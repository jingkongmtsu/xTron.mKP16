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
// BRA1 as redundant position, total RHS integrals evaluated as: 89133
// BRA2 as redundant position, total RHS integrals evaluated as: 90454
// KET1 as redundant position, total RHS integrals evaluated as: 95733
// KET2 as redundant position, total RHS integrals evaluated as: 86396
// the redundant position is: KET2
//

//
// @@@@ derivative position-direction information
// BRA1  BRA1
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####
// BRA1  BRA2
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
// BRA1  KET1
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
// KET1  KET1
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_eri_d_p_p_s_d2(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_S_Px_Px_S = 0.0E0;
  Double I_ERI_S_Py_Px_S = 0.0E0;
  Double I_ERI_S_Pz_Px_S = 0.0E0;
  Double I_ERI_S_Px_Py_S = 0.0E0;
  Double I_ERI_S_Py_Py_S = 0.0E0;
  Double I_ERI_S_Pz_Py_S = 0.0E0;
  Double I_ERI_S_Px_Pz_S = 0.0E0;
  Double I_ERI_S_Py_Pz_S = 0.0E0;
  Double I_ERI_S_Pz_Pz_S = 0.0E0;
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
  Double I_ERI_Px_S_Px_S = 0.0E0;
  Double I_ERI_Py_S_Px_S = 0.0E0;
  Double I_ERI_Pz_S_Px_S = 0.0E0;
  Double I_ERI_Px_S_Py_S = 0.0E0;
  Double I_ERI_Py_S_Py_S = 0.0E0;
  Double I_ERI_Pz_S_Py_S = 0.0E0;
  Double I_ERI_Px_S_Pz_S = 0.0E0;
  Double I_ERI_Py_S_Pz_S = 0.0E0;
  Double I_ERI_Pz_S_Pz_S = 0.0E0;
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
  Double I_ERI_D2x_S_S_S = 0.0E0;
  Double I_ERI_Dxy_S_S_S = 0.0E0;
  Double I_ERI_Dxz_S_S_S = 0.0E0;
  Double I_ERI_D2y_S_S_S = 0.0E0;
  Double I_ERI_Dyz_S_S_S = 0.0E0;
  Double I_ERI_D2z_S_S_S = 0.0E0;
  Double I_ERI_H5x_S_Px_S_aa = 0.0E0;
  Double I_ERI_H4xy_S_Px_S_aa = 0.0E0;
  Double I_ERI_H4xz_S_Px_S_aa = 0.0E0;
  Double I_ERI_H3x2y_S_Px_S_aa = 0.0E0;
  Double I_ERI_H3xyz_S_Px_S_aa = 0.0E0;
  Double I_ERI_H3x2z_S_Px_S_aa = 0.0E0;
  Double I_ERI_H2x3y_S_Px_S_aa = 0.0E0;
  Double I_ERI_H2x2yz_S_Px_S_aa = 0.0E0;
  Double I_ERI_H2xy2z_S_Px_S_aa = 0.0E0;
  Double I_ERI_H2x3z_S_Px_S_aa = 0.0E0;
  Double I_ERI_Hx4y_S_Px_S_aa = 0.0E0;
  Double I_ERI_Hx3yz_S_Px_S_aa = 0.0E0;
  Double I_ERI_Hx2y2z_S_Px_S_aa = 0.0E0;
  Double I_ERI_Hxy3z_S_Px_S_aa = 0.0E0;
  Double I_ERI_Hx4z_S_Px_S_aa = 0.0E0;
  Double I_ERI_H5y_S_Px_S_aa = 0.0E0;
  Double I_ERI_H4yz_S_Px_S_aa = 0.0E0;
  Double I_ERI_H3y2z_S_Px_S_aa = 0.0E0;
  Double I_ERI_H2y3z_S_Px_S_aa = 0.0E0;
  Double I_ERI_Hy4z_S_Px_S_aa = 0.0E0;
  Double I_ERI_H5z_S_Px_S_aa = 0.0E0;
  Double I_ERI_H5x_S_Py_S_aa = 0.0E0;
  Double I_ERI_H4xy_S_Py_S_aa = 0.0E0;
  Double I_ERI_H4xz_S_Py_S_aa = 0.0E0;
  Double I_ERI_H3x2y_S_Py_S_aa = 0.0E0;
  Double I_ERI_H3xyz_S_Py_S_aa = 0.0E0;
  Double I_ERI_H3x2z_S_Py_S_aa = 0.0E0;
  Double I_ERI_H2x3y_S_Py_S_aa = 0.0E0;
  Double I_ERI_H2x2yz_S_Py_S_aa = 0.0E0;
  Double I_ERI_H2xy2z_S_Py_S_aa = 0.0E0;
  Double I_ERI_H2x3z_S_Py_S_aa = 0.0E0;
  Double I_ERI_Hx4y_S_Py_S_aa = 0.0E0;
  Double I_ERI_Hx3yz_S_Py_S_aa = 0.0E0;
  Double I_ERI_Hx2y2z_S_Py_S_aa = 0.0E0;
  Double I_ERI_Hxy3z_S_Py_S_aa = 0.0E0;
  Double I_ERI_Hx4z_S_Py_S_aa = 0.0E0;
  Double I_ERI_H5y_S_Py_S_aa = 0.0E0;
  Double I_ERI_H4yz_S_Py_S_aa = 0.0E0;
  Double I_ERI_H3y2z_S_Py_S_aa = 0.0E0;
  Double I_ERI_H2y3z_S_Py_S_aa = 0.0E0;
  Double I_ERI_Hy4z_S_Py_S_aa = 0.0E0;
  Double I_ERI_H5z_S_Py_S_aa = 0.0E0;
  Double I_ERI_H5x_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H4xy_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H4xz_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H3x2y_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H3xyz_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H3x2z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H2x3y_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H2x2yz_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H2xy2z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H2x3z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Hx4y_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Hx3yz_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Hx2y2z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Hxy3z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Hx4z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H5y_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H4yz_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H3y2z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H2y3z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Hy4z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_H5z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G4x_S_Px_S_aa = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_aa = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_aa = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_aa = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_aa = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_aa = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_aa = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_aa = 0.0E0;
  Double I_ERI_G4y_S_Px_S_aa = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_aa = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_aa = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_aa = 0.0E0;
  Double I_ERI_G4z_S_Px_S_aa = 0.0E0;
  Double I_ERI_G4x_S_Py_S_aa = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_aa = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_aa = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_aa = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_aa = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_aa = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_aa = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_aa = 0.0E0;
  Double I_ERI_G4y_S_Py_S_aa = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_aa = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_aa = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_aa = 0.0E0;
  Double I_ERI_G4z_S_Py_S_aa = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_aa = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_aa = 0.0E0;
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
  Double I_ERI_G4x_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_ac = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_ac = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_ac = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_ac = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_ac = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_ac = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_ac = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_ac = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_ac = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_ac = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_ac = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_ac = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_ac = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_ac = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_ac = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_ac = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_ac = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_ac = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_ac = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_ac = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_ac = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_ac = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_ac = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_ac = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_ac = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_ac = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_ac = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_ac = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_ac = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_ac = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_ac = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_ac = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_ac = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_ac = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_ac = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_ac = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_ac = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_ac = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_ac = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_ac = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_ac = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_ac = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_ac = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_ac = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_ac = 0.0E0;
  Double I_ERI_G4x_S_S_S_a = 0.0E0;
  Double I_ERI_G3xy_S_S_S_a = 0.0E0;
  Double I_ERI_G3xz_S_S_S_a = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_a = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_a = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_a = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_a = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_a = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_a = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_a = 0.0E0;
  Double I_ERI_G4y_S_S_S_a = 0.0E0;
  Double I_ERI_G3yz_S_S_S_a = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_a = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_a = 0.0E0;
  Double I_ERI_G4z_S_S_S_a = 0.0E0;
  Double I_ERI_F3x_S_S_S_a = 0.0E0;
  Double I_ERI_F2xy_S_S_S_a = 0.0E0;
  Double I_ERI_F2xz_S_S_S_a = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_a = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_a = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_a = 0.0E0;
  Double I_ERI_F3y_S_S_S_a = 0.0E0;
  Double I_ERI_F2yz_S_S_S_a = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_a = 0.0E0;
  Double I_ERI_F3z_S_S_S_a = 0.0E0;
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
  Double I_ERI_Px_S_S_S = 0.0E0;
  Double I_ERI_Py_S_S_S = 0.0E0;
  Double I_ERI_Pz_S_S_S = 0.0E0;
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
  Double I_ERI_F3x_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_cc = 0.0E0;
  Double I_ERI_D2x_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S_cc = 0.0E0;
  Double I_ERI_D2y_S_F3x_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S_cc = 0.0E0;
  Double I_ERI_D2z_S_F3x_S_cc = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S_cc = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S_cc = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S_cc = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S_cc = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S_cc = 0.0E0;
  Double I_ERI_D2x_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S_cc = 0.0E0;
  Double I_ERI_D2y_S_F3y_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S_cc = 0.0E0;
  Double I_ERI_D2z_S_F3y_S_cc = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S_cc = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S_cc = 0.0E0;
  Double I_ERI_D2x_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S_cc = 0.0E0;
  Double I_ERI_D2y_S_F3z_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S_cc = 0.0E0;
  Double I_ERI_D2z_S_F3z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Px_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_c = 0.0E0;
  Double I_ERI_F3y_S_Px_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_c = 0.0E0;
  Double I_ERI_F3z_S_Px_S_c = 0.0E0;
  Double I_ERI_F3x_S_Py_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_c = 0.0E0;
  Double I_ERI_F3y_S_Py_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_c = 0.0E0;
  Double I_ERI_F3z_S_Py_S_c = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_c = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_c = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_c = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_c = 0.0E0;
  Double I_ERI_D2x_S_Px_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_c = 0.0E0;
  Double I_ERI_D2y_S_Px_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_c = 0.0E0;
  Double I_ERI_D2z_S_Px_S_c = 0.0E0;
  Double I_ERI_D2x_S_Py_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_c = 0.0E0;
  Double I_ERI_D2y_S_Py_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_c = 0.0E0;
  Double I_ERI_D2z_S_Py_S_c = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_c = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_c = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_c = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_c = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_c = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_H5x_S_Px_S_ab = 0.0E0;
  Double I_ERI_H4xy_S_Px_S_ab = 0.0E0;
  Double I_ERI_H4xz_S_Px_S_ab = 0.0E0;
  Double I_ERI_H3x2y_S_Px_S_ab = 0.0E0;
  Double I_ERI_H3xyz_S_Px_S_ab = 0.0E0;
  Double I_ERI_H3x2z_S_Px_S_ab = 0.0E0;
  Double I_ERI_H2x3y_S_Px_S_ab = 0.0E0;
  Double I_ERI_H2x2yz_S_Px_S_ab = 0.0E0;
  Double I_ERI_H2xy2z_S_Px_S_ab = 0.0E0;
  Double I_ERI_H2x3z_S_Px_S_ab = 0.0E0;
  Double I_ERI_Hx4y_S_Px_S_ab = 0.0E0;
  Double I_ERI_Hx3yz_S_Px_S_ab = 0.0E0;
  Double I_ERI_Hx2y2z_S_Px_S_ab = 0.0E0;
  Double I_ERI_Hxy3z_S_Px_S_ab = 0.0E0;
  Double I_ERI_Hx4z_S_Px_S_ab = 0.0E0;
  Double I_ERI_H5y_S_Px_S_ab = 0.0E0;
  Double I_ERI_H4yz_S_Px_S_ab = 0.0E0;
  Double I_ERI_H3y2z_S_Px_S_ab = 0.0E0;
  Double I_ERI_H2y3z_S_Px_S_ab = 0.0E0;
  Double I_ERI_Hy4z_S_Px_S_ab = 0.0E0;
  Double I_ERI_H5z_S_Px_S_ab = 0.0E0;
  Double I_ERI_H5x_S_Py_S_ab = 0.0E0;
  Double I_ERI_H4xy_S_Py_S_ab = 0.0E0;
  Double I_ERI_H4xz_S_Py_S_ab = 0.0E0;
  Double I_ERI_H3x2y_S_Py_S_ab = 0.0E0;
  Double I_ERI_H3xyz_S_Py_S_ab = 0.0E0;
  Double I_ERI_H3x2z_S_Py_S_ab = 0.0E0;
  Double I_ERI_H2x3y_S_Py_S_ab = 0.0E0;
  Double I_ERI_H2x2yz_S_Py_S_ab = 0.0E0;
  Double I_ERI_H2xy2z_S_Py_S_ab = 0.0E0;
  Double I_ERI_H2x3z_S_Py_S_ab = 0.0E0;
  Double I_ERI_Hx4y_S_Py_S_ab = 0.0E0;
  Double I_ERI_Hx3yz_S_Py_S_ab = 0.0E0;
  Double I_ERI_Hx2y2z_S_Py_S_ab = 0.0E0;
  Double I_ERI_Hxy3z_S_Py_S_ab = 0.0E0;
  Double I_ERI_Hx4z_S_Py_S_ab = 0.0E0;
  Double I_ERI_H5y_S_Py_S_ab = 0.0E0;
  Double I_ERI_H4yz_S_Py_S_ab = 0.0E0;
  Double I_ERI_H3y2z_S_Py_S_ab = 0.0E0;
  Double I_ERI_H2y3z_S_Py_S_ab = 0.0E0;
  Double I_ERI_Hy4z_S_Py_S_ab = 0.0E0;
  Double I_ERI_H5z_S_Py_S_ab = 0.0E0;
  Double I_ERI_H5x_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H4xy_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H4xz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H3x2y_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H3xyz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H3x2z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H2x3y_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H2x2yz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H2xy2z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H2x3z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Hx4y_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Hx3yz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Hx2y2z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Hxy3z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Hx4z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H5y_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H4yz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H3y2z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H2y3z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Hy4z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_H5z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G4x_S_Px_S_ab = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_ab = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_ab = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_ab = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_ab = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_ab = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_ab = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_ab = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_ab = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_ab = 0.0E0;
  Double I_ERI_G4y_S_Px_S_ab = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_ab = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_ab = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_ab = 0.0E0;
  Double I_ERI_G4z_S_Px_S_ab = 0.0E0;
  Double I_ERI_G4x_S_Py_S_ab = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_ab = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_ab = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_ab = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_ab = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_ab = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_ab = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_ab = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_ab = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_ab = 0.0E0;
  Double I_ERI_G4y_S_Py_S_ab = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_ab = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_ab = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_ab = 0.0E0;
  Double I_ERI_G4z_S_Py_S_ab = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_F3x_S_Px_S_ab = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_ab = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_ab = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_ab = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_ab = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_ab = 0.0E0;
  Double I_ERI_F3y_S_Px_S_ab = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_ab = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_ab = 0.0E0;
  Double I_ERI_F3z_S_Px_S_ab = 0.0E0;
  Double I_ERI_F3x_S_Py_S_ab = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_ab = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_ab = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_ab = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_ab = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_ab = 0.0E0;
  Double I_ERI_F3y_S_Py_S_ab = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_ab = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_ab = 0.0E0;
  Double I_ERI_F3z_S_Py_S_ab = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_ab = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_ab = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_ab = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_ab = 0.0E0;
  Double I_ERI_Px_S_Px_S_b = 0.0E0;
  Double I_ERI_Py_S_Px_S_b = 0.0E0;
  Double I_ERI_Pz_S_Px_S_b = 0.0E0;
  Double I_ERI_Px_S_Py_S_b = 0.0E0;
  Double I_ERI_Py_S_Py_S_b = 0.0E0;
  Double I_ERI_Pz_S_Py_S_b = 0.0E0;
  Double I_ERI_Px_S_Pz_S_b = 0.0E0;
  Double I_ERI_Py_S_Pz_S_b = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_b = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_bc = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_bc = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_bc = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_bc = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_bc = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_bc = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_bc = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_bc = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_bc = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_bc = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_bc = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_bc = 0.0E0;
  Double I_ERI_G4x_S_S_S_b = 0.0E0;
  Double I_ERI_G3xy_S_S_S_b = 0.0E0;
  Double I_ERI_G3xz_S_S_S_b = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_b = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_b = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_b = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_b = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_b = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_b = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_b = 0.0E0;
  Double I_ERI_G4y_S_S_S_b = 0.0E0;
  Double I_ERI_G3yz_S_S_S_b = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_b = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_b = 0.0E0;
  Double I_ERI_G4z_S_S_S_b = 0.0E0;
  Double I_ERI_F3x_S_S_S_b = 0.0E0;
  Double I_ERI_F2xy_S_S_S_b = 0.0E0;
  Double I_ERI_F2xz_S_S_S_b = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_b = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_b = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_b = 0.0E0;
  Double I_ERI_F3y_S_S_S_b = 0.0E0;
  Double I_ERI_F2yz_S_S_S_b = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_b = 0.0E0;
  Double I_ERI_F3z_S_S_S_b = 0.0E0;
  Double I_ERI_D2x_S_S_S_b = 0.0E0;
  Double I_ERI_Dxy_S_S_S_b = 0.0E0;
  Double I_ERI_Dxz_S_S_S_b = 0.0E0;
  Double I_ERI_D2y_S_S_S_b = 0.0E0;
  Double I_ERI_Dyz_S_S_S_b = 0.0E0;
  Double I_ERI_D2z_S_S_S_b = 0.0E0;
  Double I_ERI_H5x_S_Px_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_Px_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_Px_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_Px_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_Px_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H5y_S_Px_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H5z_S_Px_S_bb = 0.0E0;
  Double I_ERI_H5x_S_Py_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_Py_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_Py_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_Py_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_Py_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H5y_S_Py_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H5z_S_Py_S_bb = 0.0E0;
  Double I_ERI_H5x_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H5y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_H5z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G4x_S_Px_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_G4y_S_Px_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_G4z_S_Px_S_bb = 0.0E0;
  Double I_ERI_G4x_S_Py_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_G4y_S_Py_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_G4z_S_Py_S_bb = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F3x_S_Px_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_F3y_S_Px_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_F3z_S_Px_S_bb = 0.0E0;
  Double I_ERI_F3x_S_Py_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_F3y_S_Py_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_F3z_S_Py_S_bb = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_bb = 0.0E0;
  Double I_ERI_D2x_S_Px_S_bb = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_bb = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_bb = 0.0E0;
  Double I_ERI_D2y_S_Px_S_bb = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_bb = 0.0E0;
  Double I_ERI_D2z_S_Px_S_bb = 0.0E0;
  Double I_ERI_D2x_S_Py_S_bb = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_bb = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_bb = 0.0E0;
  Double I_ERI_D2y_S_Py_S_bb = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_bb = 0.0E0;
  Double I_ERI_D2z_S_Py_S_bb = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_bb = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_bb = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_bb = 0.0E0;

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
       * shell quartet name: SQ_ERI_G_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M2_vrr = PAX*I_ERI_F3x_S_S_S_M2_vrr+WPX*I_ERI_F3x_S_S_S_M3_vrr+3*oned2z*I_ERI_D2x_S_S_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_G3xy_S_S_S_M2_vrr = PAY*I_ERI_F3x_S_S_S_M2_vrr+WPY*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G3xz_S_S_S_M2_vrr = PAZ*I_ERI_F3x_S_S_S_M2_vrr+WPZ*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G2x2y_S_S_S_M2_vrr = PAY*I_ERI_F2xy_S_S_S_M2_vrr+WPY*I_ERI_F2xy_S_S_S_M3_vrr+oned2z*I_ERI_D2x_S_S_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_G2x2z_S_S_S_M2_vrr = PAZ*I_ERI_F2xz_S_S_S_M2_vrr+WPZ*I_ERI_F2xz_S_S_S_M3_vrr+oned2z*I_ERI_D2x_S_S_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Gx3y_S_S_S_M2_vrr = PAX*I_ERI_F3y_S_S_S_M2_vrr+WPX*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_Gx3z_S_S_S_M2_vrr = PAX*I_ERI_F3z_S_S_S_M2_vrr+WPX*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4y_S_S_S_M2_vrr = PAY*I_ERI_F3y_S_S_S_M2_vrr+WPY*I_ERI_F3y_S_S_S_M3_vrr+3*oned2z*I_ERI_D2y_S_S_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_G3yz_S_S_S_M2_vrr = PAZ*I_ERI_F3y_S_S_S_M2_vrr+WPZ*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_G2y2z_S_S_S_M2_vrr = PAZ*I_ERI_F2yz_S_S_S_M2_vrr+WPZ*I_ERI_F2yz_S_S_S_M3_vrr+oned2z*I_ERI_D2y_S_S_S_M2_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Gy3z_S_S_S_M2_vrr = PAY*I_ERI_F3z_S_S_S_M2_vrr+WPY*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4z_S_S_S_M2_vrr = PAZ*I_ERI_F3z_S_S_S_M2_vrr+WPZ*I_ERI_F3z_S_S_S_M3_vrr+3*oned2z*I_ERI_D2z_S_S_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_H_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       ************************************************************/
      Double I_ERI_H5x_S_S_S_M1_vrr = PAX*I_ERI_G4x_S_S_S_M1_vrr+WPX*I_ERI_G4x_S_S_S_M2_vrr+4*oned2z*I_ERI_F3x_S_S_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_H4xy_S_S_S_M1_vrr = PAY*I_ERI_G4x_S_S_S_M1_vrr+WPY*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H4xz_S_S_S_M1_vrr = PAZ*I_ERI_G4x_S_S_S_M1_vrr+WPZ*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H3x2y_S_S_S_M1_vrr = PAY*I_ERI_G3xy_S_S_S_M1_vrr+WPY*I_ERI_G3xy_S_S_S_M2_vrr+oned2z*I_ERI_F3x_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_H3xyz_S_S_S_M1_vrr = PAZ*I_ERI_G3xy_S_S_S_M1_vrr+WPZ*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_H3x2z_S_S_S_M1_vrr = PAZ*I_ERI_G3xz_S_S_S_M1_vrr+WPZ*I_ERI_G3xz_S_S_S_M2_vrr+oned2z*I_ERI_F3x_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_H2x3y_S_S_S_M1_vrr = PAX*I_ERI_Gx3y_S_S_S_M1_vrr+WPX*I_ERI_Gx3y_S_S_S_M2_vrr+oned2z*I_ERI_F3y_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_H2x2yz_S_S_S_M1_vrr = PAZ*I_ERI_G2x2y_S_S_S_M1_vrr+WPZ*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_H2xy2z_S_S_S_M1_vrr = PAY*I_ERI_G2x2z_S_S_S_M1_vrr+WPY*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_H2x3z_S_S_S_M1_vrr = PAX*I_ERI_Gx3z_S_S_S_M1_vrr+WPX*I_ERI_Gx3z_S_S_S_M2_vrr+oned2z*I_ERI_F3z_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4y_S_S_S_M1_vrr = PAX*I_ERI_G4y_S_S_S_M1_vrr+WPX*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_Hx3yz_S_S_S_M1_vrr = PAZ*I_ERI_Gx3y_S_S_S_M1_vrr+WPZ*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_Hx2y2z_S_S_S_M1_vrr = PAX*I_ERI_G2y2z_S_S_S_M1_vrr+WPX*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_Hxy3z_S_S_S_M1_vrr = PAY*I_ERI_Gx3z_S_S_S_M1_vrr+WPY*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4z_S_S_S_M1_vrr = PAX*I_ERI_G4z_S_S_S_M1_vrr+WPX*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_H5y_S_S_S_M1_vrr = PAY*I_ERI_G4y_S_S_S_M1_vrr+WPY*I_ERI_G4y_S_S_S_M2_vrr+4*oned2z*I_ERI_F3y_S_S_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_H4yz_S_S_S_M1_vrr = PAZ*I_ERI_G4y_S_S_S_M1_vrr+WPZ*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_H3y2z_S_S_S_M1_vrr = PAZ*I_ERI_G3yz_S_S_S_M1_vrr+WPZ*I_ERI_G3yz_S_S_S_M2_vrr+oned2z*I_ERI_F3y_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_H2y3z_S_S_S_M1_vrr = PAY*I_ERI_Gy3z_S_S_S_M1_vrr+WPY*I_ERI_Gy3z_S_S_S_M2_vrr+oned2z*I_ERI_F3z_S_S_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_Hy4z_S_S_S_M1_vrr = PAY*I_ERI_G4z_S_S_S_M1_vrr+WPY*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_H5z_S_S_S_M1_vrr = PAZ*I_ERI_G4z_S_S_S_M1_vrr+WPZ*I_ERI_G4z_S_S_S_M2_vrr+4*oned2z*I_ERI_F3z_S_S_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_H_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_S_S_vrr = PAX*I_ERI_G4x_S_S_S_vrr+WPX*I_ERI_G4x_S_S_S_M1_vrr+4*oned2z*I_ERI_F3x_S_S_S_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_S_S_vrr = PAY*I_ERI_G4x_S_S_S_vrr+WPY*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_S_S_vrr = PAZ*I_ERI_G4x_S_S_S_vrr+WPZ*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_S_S_vrr = PAY*I_ERI_G3xy_S_S_S_vrr+WPY*I_ERI_G3xy_S_S_S_M1_vrr+oned2z*I_ERI_F3x_S_S_S_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_S_S_vrr = PAZ*I_ERI_G3xy_S_S_S_vrr+WPZ*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3x2z_S_S_S_vrr = PAZ*I_ERI_G3xz_S_S_S_vrr+WPZ*I_ERI_G3xz_S_S_S_M1_vrr+oned2z*I_ERI_F3x_S_S_S_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_S_S_vrr = PAX*I_ERI_Gx3y_S_S_S_vrr+WPX*I_ERI_Gx3y_S_S_S_M1_vrr+oned2z*I_ERI_F3y_S_S_S_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_S_S_vrr = PAZ*I_ERI_G2x2y_S_S_S_vrr+WPZ*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_H2xy2z_S_S_S_vrr = PAY*I_ERI_G2x2z_S_S_S_vrr+WPY*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3z_S_S_S_vrr = PAX*I_ERI_Gx3z_S_S_S_vrr+WPX*I_ERI_Gx3z_S_S_S_M1_vrr+oned2z*I_ERI_F3z_S_S_S_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_S_S_vrr = PAX*I_ERI_G4y_S_S_S_vrr+WPX*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_S_S_vrr = PAZ*I_ERI_Gx3y_S_S_S_vrr+WPZ*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_S_S_vrr = PAX*I_ERI_G2y2z_S_S_S_vrr+WPX*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Hxy3z_S_S_S_vrr = PAY*I_ERI_Gx3z_S_S_S_vrr+WPY*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4z_S_S_S_vrr = PAX*I_ERI_G4z_S_S_S_vrr+WPX*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_S_S_vrr = PAY*I_ERI_G4y_S_S_S_vrr+WPY*I_ERI_G4y_S_S_S_M1_vrr+4*oned2z*I_ERI_F3y_S_S_S_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_S_S_vrr = PAZ*I_ERI_G4y_S_S_S_vrr+WPZ*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_S_S_vrr = PAZ*I_ERI_G3yz_S_S_S_vrr+WPZ*I_ERI_G3yz_S_S_S_M1_vrr+oned2z*I_ERI_F3y_S_S_S_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_S_S_vrr = PAY*I_ERI_Gy3z_S_S_S_vrr+WPY*I_ERI_Gy3z_S_S_S_M1_vrr+oned2z*I_ERI_F3z_S_S_S_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_S_S_vrr = PAY*I_ERI_G4z_S_S_S_vrr+WPY*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_S_S_vrr = PAZ*I_ERI_G4z_S_S_S_vrr+WPZ*I_ERI_G4z_S_S_S_M1_vrr+4*oned2z*I_ERI_F3z_S_S_S_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_H_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_Px_S_vrr = QCX*I_ERI_H5x_S_S_S_vrr+WQX*I_ERI_H5x_S_S_S_M1_vrr+5*oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_Px_S_vrr = QCX*I_ERI_H4xy_S_S_S_vrr+WQX*I_ERI_H4xy_S_S_S_M1_vrr+4*oned2k*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_Px_S_vrr = QCX*I_ERI_H4xz_S_S_S_vrr+WQX*I_ERI_H4xz_S_S_S_M1_vrr+4*oned2k*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_Px_S_vrr = QCX*I_ERI_H3x2y_S_S_S_vrr+WQX*I_ERI_H3x2y_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_Px_S_vrr = QCX*I_ERI_H3xyz_S_S_S_vrr+WQX*I_ERI_H3xyz_S_S_S_M1_vrr+3*oned2k*I_ERI_G2xyz_S_S_S_M1_vrr;
      Double I_ERI_H3x2z_S_Px_S_vrr = QCX*I_ERI_H3x2z_S_S_S_vrr+WQX*I_ERI_H3x2z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_Px_S_vrr = QCX*I_ERI_H2x3y_S_S_S_vrr+WQX*I_ERI_H2x3y_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Px_S_vrr = QCX*I_ERI_H2x2yz_S_S_S_vrr+WQX*I_ERI_H2x2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_S_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Px_S_vrr = QCX*I_ERI_H2xy2z_S_S_S_vrr+WQX*I_ERI_H2xy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gxy2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3z_S_Px_S_vrr = QCX*I_ERI_H2x3z_S_S_S_vrr+WQX*I_ERI_H2x3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_Px_S_vrr = QCX*I_ERI_Hx4y_S_S_S_vrr+WQX*I_ERI_Hx4y_S_S_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Px_S_vrr = QCX*I_ERI_Hx3yz_S_S_S_vrr+WQX*I_ERI_Hx3yz_S_S_S_M1_vrr+oned2k*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Px_S_vrr = QCX*I_ERI_Hx2y2z_S_S_S_vrr+WQX*I_ERI_Hx2y2z_S_S_S_M1_vrr+oned2k*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Px_S_vrr = QCX*I_ERI_Hxy3z_S_S_S_vrr+WQX*I_ERI_Hxy3z_S_S_S_M1_vrr+oned2k*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4z_S_Px_S_vrr = QCX*I_ERI_Hx4z_S_S_S_vrr+WQX*I_ERI_Hx4z_S_S_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_Px_S_vrr = QCX*I_ERI_H5y_S_S_S_vrr+WQX*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_Px_S_vrr = QCX*I_ERI_H4yz_S_S_S_vrr+WQX*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_Px_S_vrr = QCX*I_ERI_H3y2z_S_S_S_vrr+WQX*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_Px_S_vrr = QCX*I_ERI_H2y3z_S_S_S_vrr+WQX*I_ERI_H2y3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_Px_S_vrr = QCX*I_ERI_Hy4z_S_S_S_vrr+WQX*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_Px_S_vrr = QCX*I_ERI_H5z_S_S_S_vrr+WQX*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_H5x_S_Py_S_vrr = QCY*I_ERI_H5x_S_S_S_vrr+WQY*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_Py_S_vrr = QCY*I_ERI_H4xy_S_S_S_vrr+WQY*I_ERI_H4xy_S_S_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_Py_S_vrr = QCY*I_ERI_H4xz_S_S_S_vrr+WQY*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_Py_S_vrr = QCY*I_ERI_H3x2y_S_S_S_vrr+WQY*I_ERI_H3x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_Py_S_vrr = QCY*I_ERI_H3xyz_S_S_S_vrr+WQY*I_ERI_H3xyz_S_S_S_M1_vrr+oned2k*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_H3x2z_S_Py_S_vrr = QCY*I_ERI_H3x2z_S_S_S_vrr+WQY*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_Py_S_vrr = QCY*I_ERI_H2x3y_S_S_S_vrr+WQY*I_ERI_H2x3y_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Py_S_vrr = QCY*I_ERI_H2x2yz_S_S_S_vrr+WQY*I_ERI_H2x2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_S_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Py_S_vrr = QCY*I_ERI_H2xy2z_S_S_S_vrr+WQY*I_ERI_H2xy2z_S_S_S_M1_vrr+oned2k*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3z_S_Py_S_vrr = QCY*I_ERI_H2x3z_S_S_S_vrr+WQY*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_Py_S_vrr = QCY*I_ERI_Hx4y_S_S_S_vrr+WQY*I_ERI_Hx4y_S_S_S_M1_vrr+4*oned2k*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Py_S_vrr = QCY*I_ERI_Hx3yz_S_S_S_vrr+WQY*I_ERI_Hx3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_Gx2yz_S_S_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Py_S_vrr = QCY*I_ERI_Hx2y2z_S_S_S_vrr+WQY*I_ERI_Hx2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gxy2z_S_S_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Py_S_vrr = QCY*I_ERI_Hxy3z_S_S_S_vrr+WQY*I_ERI_Hxy3z_S_S_S_M1_vrr+oned2k*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4z_S_Py_S_vrr = QCY*I_ERI_Hx4z_S_S_S_vrr+WQY*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_Py_S_vrr = QCY*I_ERI_H5y_S_S_S_vrr+WQY*I_ERI_H5y_S_S_S_M1_vrr+5*oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_Py_S_vrr = QCY*I_ERI_H4yz_S_S_S_vrr+WQY*I_ERI_H4yz_S_S_S_M1_vrr+4*oned2k*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_Py_S_vrr = QCY*I_ERI_H3y2z_S_S_S_vrr+WQY*I_ERI_H3y2z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_Py_S_vrr = QCY*I_ERI_H2y3z_S_S_S_vrr+WQY*I_ERI_H2y3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_Py_S_vrr = QCY*I_ERI_Hy4z_S_S_S_vrr+WQY*I_ERI_Hy4z_S_S_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_Py_S_vrr = QCY*I_ERI_H5z_S_S_S_vrr+WQY*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_H5x_S_Pz_S_vrr = QCZ*I_ERI_H5x_S_S_S_vrr+WQZ*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_Pz_S_vrr = QCZ*I_ERI_H4xy_S_S_S_vrr+WQZ*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_Pz_S_vrr = QCZ*I_ERI_H4xz_S_S_S_vrr+WQZ*I_ERI_H4xz_S_S_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_Pz_S_vrr = QCZ*I_ERI_H3x2y_S_S_S_vrr+WQZ*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_Pz_S_vrr = QCZ*I_ERI_H3xyz_S_S_S_vrr+WQZ*I_ERI_H3xyz_S_S_S_M1_vrr+oned2k*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3x2z_S_Pz_S_vrr = QCZ*I_ERI_H3x2z_S_S_S_vrr+WQZ*I_ERI_H3x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_Pz_S_vrr = QCZ*I_ERI_H2x3y_S_S_S_vrr+WQZ*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Pz_S_vrr = QCZ*I_ERI_H2x2yz_S_S_S_vrr+WQZ*I_ERI_H2x2yz_S_S_S_M1_vrr+oned2k*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Pz_S_vrr = QCZ*I_ERI_H2xy2z_S_S_S_vrr+WQZ*I_ERI_H2xy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_S_S_M1_vrr;
      Double I_ERI_H2x3z_S_Pz_S_vrr = QCZ*I_ERI_H2x3z_S_S_S_vrr+WQZ*I_ERI_H2x3z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_Pz_S_vrr = QCZ*I_ERI_Hx4y_S_S_S_vrr+WQZ*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Pz_S_vrr = QCZ*I_ERI_Hx3yz_S_S_S_vrr+WQZ*I_ERI_Hx3yz_S_S_S_M1_vrr+oned2k*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Pz_S_vrr = QCZ*I_ERI_Hx2y2z_S_S_S_vrr+WQZ*I_ERI_Hx2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_S_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Pz_S_vrr = QCZ*I_ERI_Hxy3z_S_S_S_vrr+WQZ*I_ERI_Hxy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Gxy2z_S_S_S_M1_vrr;
      Double I_ERI_Hx4z_S_Pz_S_vrr = QCZ*I_ERI_Hx4z_S_S_S_vrr+WQZ*I_ERI_Hx4z_S_S_S_M1_vrr+4*oned2k*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_Pz_S_vrr = QCZ*I_ERI_H5y_S_S_S_vrr+WQZ*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_Pz_S_vrr = QCZ*I_ERI_H4yz_S_S_S_vrr+WQZ*I_ERI_H4yz_S_S_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_Pz_S_vrr = QCZ*I_ERI_H3y2z_S_S_S_vrr+WQZ*I_ERI_H3y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_Pz_S_vrr = QCZ*I_ERI_H2y3z_S_S_S_vrr+WQZ*I_ERI_H2y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_Pz_S_vrr = QCZ*I_ERI_Hy4z_S_S_S_vrr+WQZ*I_ERI_Hy4z_S_S_S_M1_vrr+4*oned2k*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_Pz_S_vrr = QCZ*I_ERI_H5z_S_S_S_vrr+WQZ*I_ERI_H5z_S_S_S_M1_vrr+5*oned2k*I_ERI_G4z_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_D_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_D2x_S_S_S += I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S += I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S += I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S += I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S += I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S += I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_P_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_P_S_aa_coefs = alpha*alpha;
      I_ERI_H5x_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H5x_S_Px_S_vrr;
      I_ERI_H4xy_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H4xy_S_Px_S_vrr;
      I_ERI_H4xz_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H4xz_S_Px_S_vrr;
      I_ERI_H3x2y_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3x2y_S_Px_S_vrr;
      I_ERI_H3xyz_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3xyz_S_Px_S_vrr;
      I_ERI_H3x2z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3x2z_S_Px_S_vrr;
      I_ERI_H2x3y_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2x3y_S_Px_S_vrr;
      I_ERI_H2x2yz_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2x2yz_S_Px_S_vrr;
      I_ERI_H2xy2z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2xy2z_S_Px_S_vrr;
      I_ERI_H2x3z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2x3z_S_Px_S_vrr;
      I_ERI_Hx4y_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx4y_S_Px_S_vrr;
      I_ERI_Hx3yz_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx3yz_S_Px_S_vrr;
      I_ERI_Hx2y2z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx2y2z_S_Px_S_vrr;
      I_ERI_Hxy3z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hxy3z_S_Px_S_vrr;
      I_ERI_Hx4z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx4z_S_Px_S_vrr;
      I_ERI_H5y_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H5y_S_Px_S_vrr;
      I_ERI_H4yz_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H4yz_S_Px_S_vrr;
      I_ERI_H3y2z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3y2z_S_Px_S_vrr;
      I_ERI_H2y3z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2y3z_S_Px_S_vrr;
      I_ERI_Hy4z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hy4z_S_Px_S_vrr;
      I_ERI_H5z_S_Px_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H5z_S_Px_S_vrr;
      I_ERI_H5x_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H5x_S_Py_S_vrr;
      I_ERI_H4xy_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H4xy_S_Py_S_vrr;
      I_ERI_H4xz_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H4xz_S_Py_S_vrr;
      I_ERI_H3x2y_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3x2y_S_Py_S_vrr;
      I_ERI_H3xyz_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3xyz_S_Py_S_vrr;
      I_ERI_H3x2z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3x2z_S_Py_S_vrr;
      I_ERI_H2x3y_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2x3y_S_Py_S_vrr;
      I_ERI_H2x2yz_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2x2yz_S_Py_S_vrr;
      I_ERI_H2xy2z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2xy2z_S_Py_S_vrr;
      I_ERI_H2x3z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2x3z_S_Py_S_vrr;
      I_ERI_Hx4y_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx4y_S_Py_S_vrr;
      I_ERI_Hx3yz_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx3yz_S_Py_S_vrr;
      I_ERI_Hx2y2z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx2y2z_S_Py_S_vrr;
      I_ERI_Hxy3z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hxy3z_S_Py_S_vrr;
      I_ERI_Hx4z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx4z_S_Py_S_vrr;
      I_ERI_H5y_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H5y_S_Py_S_vrr;
      I_ERI_H4yz_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H4yz_S_Py_S_vrr;
      I_ERI_H3y2z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3y2z_S_Py_S_vrr;
      I_ERI_H2y3z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2y3z_S_Py_S_vrr;
      I_ERI_Hy4z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hy4z_S_Py_S_vrr;
      I_ERI_H5z_S_Py_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H5z_S_Py_S_vrr;
      I_ERI_H5x_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H5x_S_Pz_S_vrr;
      I_ERI_H4xy_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H4xy_S_Pz_S_vrr;
      I_ERI_H4xz_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H4xz_S_Pz_S_vrr;
      I_ERI_H3x2y_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3x2y_S_Pz_S_vrr;
      I_ERI_H3xyz_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3xyz_S_Pz_S_vrr;
      I_ERI_H3x2z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3x2z_S_Pz_S_vrr;
      I_ERI_H2x3y_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2x3y_S_Pz_S_vrr;
      I_ERI_H2x2yz_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2x2yz_S_Pz_S_vrr;
      I_ERI_H2xy2z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2xy2z_S_Pz_S_vrr;
      I_ERI_H2x3z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2x3z_S_Pz_S_vrr;
      I_ERI_Hx4y_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx4y_S_Pz_S_vrr;
      I_ERI_Hx3yz_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx3yz_S_Pz_S_vrr;
      I_ERI_Hx2y2z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx2y2z_S_Pz_S_vrr;
      I_ERI_Hxy3z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hxy3z_S_Pz_S_vrr;
      I_ERI_Hx4z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hx4z_S_Pz_S_vrr;
      I_ERI_H5y_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H5y_S_Pz_S_vrr;
      I_ERI_H4yz_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H4yz_S_Pz_S_vrr;
      I_ERI_H3y2z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H3y2z_S_Pz_S_vrr;
      I_ERI_H2y3z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H2y3z_S_Pz_S_vrr;
      I_ERI_Hy4z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_Hy4z_S_Pz_S_vrr;
      I_ERI_H5z_S_Pz_S_aa += SQ_ERI_H_S_P_S_aa_coefs*I_ERI_H5z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_aa_coefs = alpha*alpha;
      I_ERI_G4x_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_aa += SQ_ERI_G_S_P_S_aa_coefs*I_ERI_G4z_S_Pz_S_vrr;

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
       * shell quartet name: SQ_ERI_G_S_D_S_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_ac_coefs = alpha*gamma;
      I_ERI_G4x_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_ac += SQ_ERI_G_S_D_S_ac_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_ac_coefs = alpha*gamma;
      I_ERI_F3x_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_ac += SQ_ERI_F_S_D_S_ac_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_a_coefs = alpha;
      I_ERI_G4x_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_a += SQ_ERI_G_S_S_S_a_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_a_coefs = alpha;
      I_ERI_F3x_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_a += SQ_ERI_F_S_S_S_a_coefs*I_ERI_F3z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_P_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_Px_S_S_S += I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S += I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S += I_ERI_Pz_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_F_S_F_S_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_cc_coefs = gamma*gamma;
      I_ERI_F3x_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_cc += SQ_ERI_F_S_F_S_cc_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_F_S_cc_coefs = gamma*gamma;
      I_ERI_D2x_S_F3x_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S_cc += SQ_ERI_D_S_F_S_cc_coefs*I_ERI_D2z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_c_coefs = gamma;
      I_ERI_F3x_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_c += SQ_ERI_F_S_P_S_c_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_c_coefs = gamma;
      I_ERI_D2x_S_Px_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_c += SQ_ERI_D_S_P_S_c_coefs*I_ERI_D2z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_P_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_P_S_ab_coefs = alpha*beta;
      I_ERI_H5x_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H5x_S_Px_S_vrr;
      I_ERI_H4xy_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H4xy_S_Px_S_vrr;
      I_ERI_H4xz_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H4xz_S_Px_S_vrr;
      I_ERI_H3x2y_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3x2y_S_Px_S_vrr;
      I_ERI_H3xyz_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3xyz_S_Px_S_vrr;
      I_ERI_H3x2z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3x2z_S_Px_S_vrr;
      I_ERI_H2x3y_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2x3y_S_Px_S_vrr;
      I_ERI_H2x2yz_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2x2yz_S_Px_S_vrr;
      I_ERI_H2xy2z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2xy2z_S_Px_S_vrr;
      I_ERI_H2x3z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2x3z_S_Px_S_vrr;
      I_ERI_Hx4y_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx4y_S_Px_S_vrr;
      I_ERI_Hx3yz_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx3yz_S_Px_S_vrr;
      I_ERI_Hx2y2z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx2y2z_S_Px_S_vrr;
      I_ERI_Hxy3z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hxy3z_S_Px_S_vrr;
      I_ERI_Hx4z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx4z_S_Px_S_vrr;
      I_ERI_H5y_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H5y_S_Px_S_vrr;
      I_ERI_H4yz_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H4yz_S_Px_S_vrr;
      I_ERI_H3y2z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3y2z_S_Px_S_vrr;
      I_ERI_H2y3z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2y3z_S_Px_S_vrr;
      I_ERI_Hy4z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hy4z_S_Px_S_vrr;
      I_ERI_H5z_S_Px_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H5z_S_Px_S_vrr;
      I_ERI_H5x_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H5x_S_Py_S_vrr;
      I_ERI_H4xy_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H4xy_S_Py_S_vrr;
      I_ERI_H4xz_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H4xz_S_Py_S_vrr;
      I_ERI_H3x2y_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3x2y_S_Py_S_vrr;
      I_ERI_H3xyz_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3xyz_S_Py_S_vrr;
      I_ERI_H3x2z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3x2z_S_Py_S_vrr;
      I_ERI_H2x3y_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2x3y_S_Py_S_vrr;
      I_ERI_H2x2yz_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2x2yz_S_Py_S_vrr;
      I_ERI_H2xy2z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2xy2z_S_Py_S_vrr;
      I_ERI_H2x3z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2x3z_S_Py_S_vrr;
      I_ERI_Hx4y_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx4y_S_Py_S_vrr;
      I_ERI_Hx3yz_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx3yz_S_Py_S_vrr;
      I_ERI_Hx2y2z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx2y2z_S_Py_S_vrr;
      I_ERI_Hxy3z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hxy3z_S_Py_S_vrr;
      I_ERI_Hx4z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx4z_S_Py_S_vrr;
      I_ERI_H5y_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H5y_S_Py_S_vrr;
      I_ERI_H4yz_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H4yz_S_Py_S_vrr;
      I_ERI_H3y2z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3y2z_S_Py_S_vrr;
      I_ERI_H2y3z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2y3z_S_Py_S_vrr;
      I_ERI_Hy4z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hy4z_S_Py_S_vrr;
      I_ERI_H5z_S_Py_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H5z_S_Py_S_vrr;
      I_ERI_H5x_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H5x_S_Pz_S_vrr;
      I_ERI_H4xy_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H4xy_S_Pz_S_vrr;
      I_ERI_H4xz_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H4xz_S_Pz_S_vrr;
      I_ERI_H3x2y_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3x2y_S_Pz_S_vrr;
      I_ERI_H3xyz_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3xyz_S_Pz_S_vrr;
      I_ERI_H3x2z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3x2z_S_Pz_S_vrr;
      I_ERI_H2x3y_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2x3y_S_Pz_S_vrr;
      I_ERI_H2x2yz_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2x2yz_S_Pz_S_vrr;
      I_ERI_H2xy2z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2xy2z_S_Pz_S_vrr;
      I_ERI_H2x3z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2x3z_S_Pz_S_vrr;
      I_ERI_Hx4y_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx4y_S_Pz_S_vrr;
      I_ERI_Hx3yz_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx3yz_S_Pz_S_vrr;
      I_ERI_Hx2y2z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx2y2z_S_Pz_S_vrr;
      I_ERI_Hxy3z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hxy3z_S_Pz_S_vrr;
      I_ERI_Hx4z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hx4z_S_Pz_S_vrr;
      I_ERI_H5y_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H5y_S_Pz_S_vrr;
      I_ERI_H4yz_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H4yz_S_Pz_S_vrr;
      I_ERI_H3y2z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H3y2z_S_Pz_S_vrr;
      I_ERI_H2y3z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H2y3z_S_Pz_S_vrr;
      I_ERI_Hy4z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_Hy4z_S_Pz_S_vrr;
      I_ERI_H5z_S_Pz_S_ab += SQ_ERI_H_S_P_S_ab_coefs*I_ERI_H5z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_ab_coefs = alpha*beta;
      I_ERI_G4x_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_ab += SQ_ERI_G_S_P_S_ab_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_ab_coefs = alpha*beta;
      I_ERI_F3x_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_ab += SQ_ERI_F_S_P_S_ab_coefs*I_ERI_F3z_S_Pz_S_vrr;

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

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_bc_coefs = beta*gamma;
      I_ERI_G4x_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_bc += SQ_ERI_G_S_D_S_bc_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_bc_coefs = beta*gamma;
      I_ERI_F3x_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_bc += SQ_ERI_F_S_D_S_bc_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_bc_coefs = beta*gamma;
      I_ERI_D2x_S_D2x_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_bc += SQ_ERI_D_S_D_S_bc_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_b_coefs = beta;
      I_ERI_G4x_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_b += SQ_ERI_G_S_S_S_b_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_b_coefs = beta;
      I_ERI_F3x_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_b += SQ_ERI_F_S_S_S_b_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_b_coefs = beta;
      I_ERI_D2x_S_S_S_b += SQ_ERI_D_S_S_S_b_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_b += SQ_ERI_D_S_S_S_b_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_b += SQ_ERI_D_S_S_S_b_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_b += SQ_ERI_D_S_S_S_b_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_b += SQ_ERI_D_S_S_S_b_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_b += SQ_ERI_D_S_S_S_b_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_P_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_P_S_bb_coefs = beta*beta;
      I_ERI_H5x_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5x_S_Px_S_vrr;
      I_ERI_H4xy_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xy_S_Px_S_vrr;
      I_ERI_H4xz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xz_S_Px_S_vrr;
      I_ERI_H3x2y_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2y_S_Px_S_vrr;
      I_ERI_H3xyz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3xyz_S_Px_S_vrr;
      I_ERI_H3x2z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2z_S_Px_S_vrr;
      I_ERI_H2x3y_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3y_S_Px_S_vrr;
      I_ERI_H2x2yz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x2yz_S_Px_S_vrr;
      I_ERI_H2xy2z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2xy2z_S_Px_S_vrr;
      I_ERI_H2x3z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3z_S_Px_S_vrr;
      I_ERI_Hx4y_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4y_S_Px_S_vrr;
      I_ERI_Hx3yz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx3yz_S_Px_S_vrr;
      I_ERI_Hx2y2z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx2y2z_S_Px_S_vrr;
      I_ERI_Hxy3z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hxy3z_S_Px_S_vrr;
      I_ERI_Hx4z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4z_S_Px_S_vrr;
      I_ERI_H5y_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5y_S_Px_S_vrr;
      I_ERI_H4yz_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4yz_S_Px_S_vrr;
      I_ERI_H3y2z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3y2z_S_Px_S_vrr;
      I_ERI_H2y3z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2y3z_S_Px_S_vrr;
      I_ERI_Hy4z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hy4z_S_Px_S_vrr;
      I_ERI_H5z_S_Px_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5z_S_Px_S_vrr;
      I_ERI_H5x_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5x_S_Py_S_vrr;
      I_ERI_H4xy_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xy_S_Py_S_vrr;
      I_ERI_H4xz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xz_S_Py_S_vrr;
      I_ERI_H3x2y_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2y_S_Py_S_vrr;
      I_ERI_H3xyz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3xyz_S_Py_S_vrr;
      I_ERI_H3x2z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2z_S_Py_S_vrr;
      I_ERI_H2x3y_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3y_S_Py_S_vrr;
      I_ERI_H2x2yz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x2yz_S_Py_S_vrr;
      I_ERI_H2xy2z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2xy2z_S_Py_S_vrr;
      I_ERI_H2x3z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3z_S_Py_S_vrr;
      I_ERI_Hx4y_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4y_S_Py_S_vrr;
      I_ERI_Hx3yz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx3yz_S_Py_S_vrr;
      I_ERI_Hx2y2z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx2y2z_S_Py_S_vrr;
      I_ERI_Hxy3z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hxy3z_S_Py_S_vrr;
      I_ERI_Hx4z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4z_S_Py_S_vrr;
      I_ERI_H5y_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5y_S_Py_S_vrr;
      I_ERI_H4yz_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4yz_S_Py_S_vrr;
      I_ERI_H3y2z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3y2z_S_Py_S_vrr;
      I_ERI_H2y3z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2y3z_S_Py_S_vrr;
      I_ERI_Hy4z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hy4z_S_Py_S_vrr;
      I_ERI_H5z_S_Py_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5z_S_Py_S_vrr;
      I_ERI_H5x_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5x_S_Pz_S_vrr;
      I_ERI_H4xy_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xy_S_Pz_S_vrr;
      I_ERI_H4xz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4xz_S_Pz_S_vrr;
      I_ERI_H3x2y_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2y_S_Pz_S_vrr;
      I_ERI_H3xyz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3xyz_S_Pz_S_vrr;
      I_ERI_H3x2z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3x2z_S_Pz_S_vrr;
      I_ERI_H2x3y_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3y_S_Pz_S_vrr;
      I_ERI_H2x2yz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x2yz_S_Pz_S_vrr;
      I_ERI_H2xy2z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2xy2z_S_Pz_S_vrr;
      I_ERI_H2x3z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2x3z_S_Pz_S_vrr;
      I_ERI_Hx4y_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4y_S_Pz_S_vrr;
      I_ERI_Hx3yz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx3yz_S_Pz_S_vrr;
      I_ERI_Hx2y2z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx2y2z_S_Pz_S_vrr;
      I_ERI_Hxy3z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hxy3z_S_Pz_S_vrr;
      I_ERI_Hx4z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hx4z_S_Pz_S_vrr;
      I_ERI_H5y_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5y_S_Pz_S_vrr;
      I_ERI_H4yz_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H4yz_S_Pz_S_vrr;
      I_ERI_H3y2z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H3y2z_S_Pz_S_vrr;
      I_ERI_H2y3z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H2y3z_S_Pz_S_vrr;
      I_ERI_Hy4z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_Hy4z_S_Pz_S_vrr;
      I_ERI_H5z_S_Pz_S_bb += SQ_ERI_H_S_P_S_bb_coefs*I_ERI_H5z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_bb_coefs = beta*beta;
      I_ERI_G4x_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_bb += SQ_ERI_G_S_P_S_bb_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_bb_coefs = beta*beta;
      I_ERI_F3x_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_bb += SQ_ERI_F_S_P_S_bb_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_bb_coefs = beta*beta;
      I_ERI_D2x_S_Px_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_bb += SQ_ERI_D_S_P_S_bb_coefs*I_ERI_D2z_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_P_P_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   * RHS shell quartet name: SQ_ERI_P_S_S_S
   ************************************************************/
  Double I_ERI_Px_Px_S_S = I_ERI_D2x_S_S_S+ABX*I_ERI_Px_S_S_S;
  Double I_ERI_Py_Px_S_S = I_ERI_Dxy_S_S_S+ABX*I_ERI_Py_S_S_S;
  Double I_ERI_Pz_Px_S_S = I_ERI_Dxz_S_S_S+ABX*I_ERI_Pz_S_S_S;
  Double I_ERI_Px_Py_S_S = I_ERI_Dxy_S_S_S+ABY*I_ERI_Px_S_S_S;
  Double I_ERI_Py_Py_S_S = I_ERI_D2y_S_S_S+ABY*I_ERI_Py_S_S_S;
  Double I_ERI_Pz_Py_S_S = I_ERI_Dyz_S_S_S+ABY*I_ERI_Pz_S_S_S;
  Double I_ERI_Px_Pz_S_S = I_ERI_Dxz_S_S_S+ABZ*I_ERI_Px_S_S_S;
  Double I_ERI_Py_Pz_S_S = I_ERI_Dyz_S_S_S+ABZ*I_ERI_Py_S_S_S;
  Double I_ERI_Pz_Pz_S_S = I_ERI_D2z_S_S_S+ABZ*I_ERI_Pz_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_D_S_P_S_a
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_a = I_ERI_F3x_S_Px_S_a+ABX*I_ERI_D2x_S_Px_S_a;
  Double I_ERI_Dxy_Px_Px_S_a = I_ERI_F2xy_S_Px_S_a+ABX*I_ERI_Dxy_S_Px_S_a;
  Double I_ERI_Dxz_Px_Px_S_a = I_ERI_F2xz_S_Px_S_a+ABX*I_ERI_Dxz_S_Px_S_a;
  Double I_ERI_D2y_Px_Px_S_a = I_ERI_Fx2y_S_Px_S_a+ABX*I_ERI_D2y_S_Px_S_a;
  Double I_ERI_Dyz_Px_Px_S_a = I_ERI_Fxyz_S_Px_S_a+ABX*I_ERI_Dyz_S_Px_S_a;
  Double I_ERI_D2z_Px_Px_S_a = I_ERI_Fx2z_S_Px_S_a+ABX*I_ERI_D2z_S_Px_S_a;
  Double I_ERI_D2x_Py_Px_S_a = I_ERI_F2xy_S_Px_S_a+ABY*I_ERI_D2x_S_Px_S_a;
  Double I_ERI_Dxy_Py_Px_S_a = I_ERI_Fx2y_S_Px_S_a+ABY*I_ERI_Dxy_S_Px_S_a;
  Double I_ERI_Dxz_Py_Px_S_a = I_ERI_Fxyz_S_Px_S_a+ABY*I_ERI_Dxz_S_Px_S_a;
  Double I_ERI_D2y_Py_Px_S_a = I_ERI_F3y_S_Px_S_a+ABY*I_ERI_D2y_S_Px_S_a;
  Double I_ERI_Dyz_Py_Px_S_a = I_ERI_F2yz_S_Px_S_a+ABY*I_ERI_Dyz_S_Px_S_a;
  Double I_ERI_D2z_Py_Px_S_a = I_ERI_Fy2z_S_Px_S_a+ABY*I_ERI_D2z_S_Px_S_a;
  Double I_ERI_D2x_Pz_Px_S_a = I_ERI_F2xz_S_Px_S_a+ABZ*I_ERI_D2x_S_Px_S_a;
  Double I_ERI_Dxy_Pz_Px_S_a = I_ERI_Fxyz_S_Px_S_a+ABZ*I_ERI_Dxy_S_Px_S_a;
  Double I_ERI_Dxz_Pz_Px_S_a = I_ERI_Fx2z_S_Px_S_a+ABZ*I_ERI_Dxz_S_Px_S_a;
  Double I_ERI_D2y_Pz_Px_S_a = I_ERI_F2yz_S_Px_S_a+ABZ*I_ERI_D2y_S_Px_S_a;
  Double I_ERI_Dyz_Pz_Px_S_a = I_ERI_Fy2z_S_Px_S_a+ABZ*I_ERI_Dyz_S_Px_S_a;
  Double I_ERI_D2z_Pz_Px_S_a = I_ERI_F3z_S_Px_S_a+ABZ*I_ERI_D2z_S_Px_S_a;
  Double I_ERI_D2x_Px_Py_S_a = I_ERI_F3x_S_Py_S_a+ABX*I_ERI_D2x_S_Py_S_a;
  Double I_ERI_Dxy_Px_Py_S_a = I_ERI_F2xy_S_Py_S_a+ABX*I_ERI_Dxy_S_Py_S_a;
  Double I_ERI_Dxz_Px_Py_S_a = I_ERI_F2xz_S_Py_S_a+ABX*I_ERI_Dxz_S_Py_S_a;
  Double I_ERI_D2y_Px_Py_S_a = I_ERI_Fx2y_S_Py_S_a+ABX*I_ERI_D2y_S_Py_S_a;
  Double I_ERI_Dyz_Px_Py_S_a = I_ERI_Fxyz_S_Py_S_a+ABX*I_ERI_Dyz_S_Py_S_a;
  Double I_ERI_D2z_Px_Py_S_a = I_ERI_Fx2z_S_Py_S_a+ABX*I_ERI_D2z_S_Py_S_a;
  Double I_ERI_D2x_Py_Py_S_a = I_ERI_F2xy_S_Py_S_a+ABY*I_ERI_D2x_S_Py_S_a;
  Double I_ERI_Dxy_Py_Py_S_a = I_ERI_Fx2y_S_Py_S_a+ABY*I_ERI_Dxy_S_Py_S_a;
  Double I_ERI_Dxz_Py_Py_S_a = I_ERI_Fxyz_S_Py_S_a+ABY*I_ERI_Dxz_S_Py_S_a;
  Double I_ERI_D2y_Py_Py_S_a = I_ERI_F3y_S_Py_S_a+ABY*I_ERI_D2y_S_Py_S_a;
  Double I_ERI_Dyz_Py_Py_S_a = I_ERI_F2yz_S_Py_S_a+ABY*I_ERI_Dyz_S_Py_S_a;
  Double I_ERI_D2z_Py_Py_S_a = I_ERI_Fy2z_S_Py_S_a+ABY*I_ERI_D2z_S_Py_S_a;
  Double I_ERI_D2x_Pz_Py_S_a = I_ERI_F2xz_S_Py_S_a+ABZ*I_ERI_D2x_S_Py_S_a;
  Double I_ERI_Dxy_Pz_Py_S_a = I_ERI_Fxyz_S_Py_S_a+ABZ*I_ERI_Dxy_S_Py_S_a;
  Double I_ERI_Dxz_Pz_Py_S_a = I_ERI_Fx2z_S_Py_S_a+ABZ*I_ERI_Dxz_S_Py_S_a;
  Double I_ERI_D2y_Pz_Py_S_a = I_ERI_F2yz_S_Py_S_a+ABZ*I_ERI_D2y_S_Py_S_a;
  Double I_ERI_Dyz_Pz_Py_S_a = I_ERI_Fy2z_S_Py_S_a+ABZ*I_ERI_Dyz_S_Py_S_a;
  Double I_ERI_D2z_Pz_Py_S_a = I_ERI_F3z_S_Py_S_a+ABZ*I_ERI_D2z_S_Py_S_a;
  Double I_ERI_D2x_Px_Pz_S_a = I_ERI_F3x_S_Pz_S_a+ABX*I_ERI_D2x_S_Pz_S_a;
  Double I_ERI_Dxy_Px_Pz_S_a = I_ERI_F2xy_S_Pz_S_a+ABX*I_ERI_Dxy_S_Pz_S_a;
  Double I_ERI_Dxz_Px_Pz_S_a = I_ERI_F2xz_S_Pz_S_a+ABX*I_ERI_Dxz_S_Pz_S_a;
  Double I_ERI_D2y_Px_Pz_S_a = I_ERI_Fx2y_S_Pz_S_a+ABX*I_ERI_D2y_S_Pz_S_a;
  Double I_ERI_Dyz_Px_Pz_S_a = I_ERI_Fxyz_S_Pz_S_a+ABX*I_ERI_Dyz_S_Pz_S_a;
  Double I_ERI_D2z_Px_Pz_S_a = I_ERI_Fx2z_S_Pz_S_a+ABX*I_ERI_D2z_S_Pz_S_a;
  Double I_ERI_D2x_Py_Pz_S_a = I_ERI_F2xy_S_Pz_S_a+ABY*I_ERI_D2x_S_Pz_S_a;
  Double I_ERI_Dxy_Py_Pz_S_a = I_ERI_Fx2y_S_Pz_S_a+ABY*I_ERI_Dxy_S_Pz_S_a;
  Double I_ERI_Dxz_Py_Pz_S_a = I_ERI_Fxyz_S_Pz_S_a+ABY*I_ERI_Dxz_S_Pz_S_a;
  Double I_ERI_D2y_Py_Pz_S_a = I_ERI_F3y_S_Pz_S_a+ABY*I_ERI_D2y_S_Pz_S_a;
  Double I_ERI_Dyz_Py_Pz_S_a = I_ERI_F2yz_S_Pz_S_a+ABY*I_ERI_Dyz_S_Pz_S_a;
  Double I_ERI_D2z_Py_Pz_S_a = I_ERI_Fy2z_S_Pz_S_a+ABY*I_ERI_D2z_S_Pz_S_a;
  Double I_ERI_D2x_Pz_Pz_S_a = I_ERI_F2xz_S_Pz_S_a+ABZ*I_ERI_D2x_S_Pz_S_a;
  Double I_ERI_Dxy_Pz_Pz_S_a = I_ERI_Fxyz_S_Pz_S_a+ABZ*I_ERI_Dxy_S_Pz_S_a;
  Double I_ERI_Dxz_Pz_Pz_S_a = I_ERI_Fx2z_S_Pz_S_a+ABZ*I_ERI_Dxz_S_Pz_S_a;
  Double I_ERI_D2y_Pz_Pz_S_a = I_ERI_F2yz_S_Pz_S_a+ABZ*I_ERI_D2y_S_Pz_S_a;
  Double I_ERI_Dyz_Pz_Pz_S_a = I_ERI_Fy2z_S_Pz_S_a+ABZ*I_ERI_Dyz_S_Pz_S_a;
  Double I_ERI_D2z_Pz_Pz_S_a = I_ERI_F3z_S_Pz_S_a+ABZ*I_ERI_D2z_S_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_a
   * RHS shell quartet name: SQ_ERI_F_S_S_S_a
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_a = I_ERI_G4x_S_S_S_a+ABX*I_ERI_F3x_S_S_S_a;
  Double I_ERI_F2xy_Px_S_S_a = I_ERI_G3xy_S_S_S_a+ABX*I_ERI_F2xy_S_S_S_a;
  Double I_ERI_F2xz_Px_S_S_a = I_ERI_G3xz_S_S_S_a+ABX*I_ERI_F2xz_S_S_S_a;
  Double I_ERI_Fx2y_Px_S_S_a = I_ERI_G2x2y_S_S_S_a+ABX*I_ERI_Fx2y_S_S_S_a;
  Double I_ERI_Fxyz_Px_S_S_a = I_ERI_G2xyz_S_S_S_a+ABX*I_ERI_Fxyz_S_S_S_a;
  Double I_ERI_Fx2z_Px_S_S_a = I_ERI_G2x2z_S_S_S_a+ABX*I_ERI_Fx2z_S_S_S_a;
  Double I_ERI_F3y_Px_S_S_a = I_ERI_Gx3y_S_S_S_a+ABX*I_ERI_F3y_S_S_S_a;
  Double I_ERI_F2yz_Px_S_S_a = I_ERI_Gx2yz_S_S_S_a+ABX*I_ERI_F2yz_S_S_S_a;
  Double I_ERI_Fy2z_Px_S_S_a = I_ERI_Gxy2z_S_S_S_a+ABX*I_ERI_Fy2z_S_S_S_a;
  Double I_ERI_F3z_Px_S_S_a = I_ERI_Gx3z_S_S_S_a+ABX*I_ERI_F3z_S_S_S_a;
  Double I_ERI_F3x_Py_S_S_a = I_ERI_G3xy_S_S_S_a+ABY*I_ERI_F3x_S_S_S_a;
  Double I_ERI_F2xy_Py_S_S_a = I_ERI_G2x2y_S_S_S_a+ABY*I_ERI_F2xy_S_S_S_a;
  Double I_ERI_F2xz_Py_S_S_a = I_ERI_G2xyz_S_S_S_a+ABY*I_ERI_F2xz_S_S_S_a;
  Double I_ERI_Fx2y_Py_S_S_a = I_ERI_Gx3y_S_S_S_a+ABY*I_ERI_Fx2y_S_S_S_a;
  Double I_ERI_Fxyz_Py_S_S_a = I_ERI_Gx2yz_S_S_S_a+ABY*I_ERI_Fxyz_S_S_S_a;
  Double I_ERI_Fx2z_Py_S_S_a = I_ERI_Gxy2z_S_S_S_a+ABY*I_ERI_Fx2z_S_S_S_a;
  Double I_ERI_F3y_Py_S_S_a = I_ERI_G4y_S_S_S_a+ABY*I_ERI_F3y_S_S_S_a;
  Double I_ERI_F2yz_Py_S_S_a = I_ERI_G3yz_S_S_S_a+ABY*I_ERI_F2yz_S_S_S_a;
  Double I_ERI_Fy2z_Py_S_S_a = I_ERI_G2y2z_S_S_S_a+ABY*I_ERI_Fy2z_S_S_S_a;
  Double I_ERI_F3z_Py_S_S_a = I_ERI_Gy3z_S_S_S_a+ABY*I_ERI_F3z_S_S_S_a;
  Double I_ERI_F3x_Pz_S_S_a = I_ERI_G3xz_S_S_S_a+ABZ*I_ERI_F3x_S_S_S_a;
  Double I_ERI_F2xy_Pz_S_S_a = I_ERI_G2xyz_S_S_S_a+ABZ*I_ERI_F2xy_S_S_S_a;
  Double I_ERI_F2xz_Pz_S_S_a = I_ERI_G2x2z_S_S_S_a+ABZ*I_ERI_F2xz_S_S_S_a;
  Double I_ERI_Fx2y_Pz_S_S_a = I_ERI_Gx2yz_S_S_S_a+ABZ*I_ERI_Fx2y_S_S_S_a;
  Double I_ERI_Fxyz_Pz_S_S_a = I_ERI_Gxy2z_S_S_S_a+ABZ*I_ERI_Fxyz_S_S_S_a;
  Double I_ERI_Fx2z_Pz_S_S_a = I_ERI_Gx3z_S_S_S_a+ABZ*I_ERI_Fx2z_S_S_S_a;
  Double I_ERI_F3y_Pz_S_S_a = I_ERI_G3yz_S_S_S_a+ABZ*I_ERI_F3y_S_S_S_a;
  Double I_ERI_F2yz_Pz_S_S_a = I_ERI_G2y2z_S_S_S_a+ABZ*I_ERI_F2yz_S_S_S_a;
  Double I_ERI_Fy2z_Pz_S_S_a = I_ERI_Gy3z_S_S_S_a+ABZ*I_ERI_Fy2z_S_S_S_a;
  Double I_ERI_F3z_Pz_S_S_a = I_ERI_G4z_S_S_S_a+ABZ*I_ERI_F3z_S_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_P_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S_b
   ************************************************************/
  Double I_ERI_Px_Px_Px_S_b = I_ERI_D2x_S_Px_S_b+ABX*I_ERI_Px_S_Px_S_b;
  Double I_ERI_Py_Px_Px_S_b = I_ERI_Dxy_S_Px_S_b+ABX*I_ERI_Py_S_Px_S_b;
  Double I_ERI_Pz_Px_Px_S_b = I_ERI_Dxz_S_Px_S_b+ABX*I_ERI_Pz_S_Px_S_b;
  Double I_ERI_Px_Py_Px_S_b = I_ERI_Dxy_S_Px_S_b+ABY*I_ERI_Px_S_Px_S_b;
  Double I_ERI_Py_Py_Px_S_b = I_ERI_D2y_S_Px_S_b+ABY*I_ERI_Py_S_Px_S_b;
  Double I_ERI_Pz_Py_Px_S_b = I_ERI_Dyz_S_Px_S_b+ABY*I_ERI_Pz_S_Px_S_b;
  Double I_ERI_Px_Pz_Px_S_b = I_ERI_Dxz_S_Px_S_b+ABZ*I_ERI_Px_S_Px_S_b;
  Double I_ERI_Py_Pz_Px_S_b = I_ERI_Dyz_S_Px_S_b+ABZ*I_ERI_Py_S_Px_S_b;
  Double I_ERI_Pz_Pz_Px_S_b = I_ERI_D2z_S_Px_S_b+ABZ*I_ERI_Pz_S_Px_S_b;
  Double I_ERI_Px_Px_Py_S_b = I_ERI_D2x_S_Py_S_b+ABX*I_ERI_Px_S_Py_S_b;
  Double I_ERI_Py_Px_Py_S_b = I_ERI_Dxy_S_Py_S_b+ABX*I_ERI_Py_S_Py_S_b;
  Double I_ERI_Pz_Px_Py_S_b = I_ERI_Dxz_S_Py_S_b+ABX*I_ERI_Pz_S_Py_S_b;
  Double I_ERI_Px_Py_Py_S_b = I_ERI_Dxy_S_Py_S_b+ABY*I_ERI_Px_S_Py_S_b;
  Double I_ERI_Py_Py_Py_S_b = I_ERI_D2y_S_Py_S_b+ABY*I_ERI_Py_S_Py_S_b;
  Double I_ERI_Pz_Py_Py_S_b = I_ERI_Dyz_S_Py_S_b+ABY*I_ERI_Pz_S_Py_S_b;
  Double I_ERI_Px_Pz_Py_S_b = I_ERI_Dxz_S_Py_S_b+ABZ*I_ERI_Px_S_Py_S_b;
  Double I_ERI_Py_Pz_Py_S_b = I_ERI_Dyz_S_Py_S_b+ABZ*I_ERI_Py_S_Py_S_b;
  Double I_ERI_Pz_Pz_Py_S_b = I_ERI_D2z_S_Py_S_b+ABZ*I_ERI_Pz_S_Py_S_b;
  Double I_ERI_Px_Px_Pz_S_b = I_ERI_D2x_S_Pz_S_b+ABX*I_ERI_Px_S_Pz_S_b;
  Double I_ERI_Py_Px_Pz_S_b = I_ERI_Dxy_S_Pz_S_b+ABX*I_ERI_Py_S_Pz_S_b;
  Double I_ERI_Pz_Px_Pz_S_b = I_ERI_Dxz_S_Pz_S_b+ABX*I_ERI_Pz_S_Pz_S_b;
  Double I_ERI_Px_Py_Pz_S_b = I_ERI_Dxy_S_Pz_S_b+ABY*I_ERI_Px_S_Pz_S_b;
  Double I_ERI_Py_Py_Pz_S_b = I_ERI_D2y_S_Pz_S_b+ABY*I_ERI_Py_S_Pz_S_b;
  Double I_ERI_Pz_Py_Pz_S_b = I_ERI_Dyz_S_Pz_S_b+ABY*I_ERI_Pz_S_Pz_S_b;
  Double I_ERI_Px_Pz_Pz_S_b = I_ERI_Dxz_S_Pz_S_b+ABZ*I_ERI_Px_S_Pz_S_b;
  Double I_ERI_Py_Pz_Pz_S_b = I_ERI_Dyz_S_Pz_S_b+ABZ*I_ERI_Py_S_Pz_S_b;
  Double I_ERI_Pz_Pz_Pz_S_b = I_ERI_D2z_S_Pz_S_b+ABZ*I_ERI_Pz_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_S_S_b
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_b = I_ERI_F3x_S_S_S_b+ABX*I_ERI_D2x_S_S_S_b;
  Double I_ERI_Dxy_Px_S_S_b = I_ERI_F2xy_S_S_S_b+ABX*I_ERI_Dxy_S_S_S_b;
  Double I_ERI_Dxz_Px_S_S_b = I_ERI_F2xz_S_S_S_b+ABX*I_ERI_Dxz_S_S_S_b;
  Double I_ERI_D2y_Px_S_S_b = I_ERI_Fx2y_S_S_S_b+ABX*I_ERI_D2y_S_S_S_b;
  Double I_ERI_Dyz_Px_S_S_b = I_ERI_Fxyz_S_S_S_b+ABX*I_ERI_Dyz_S_S_S_b;
  Double I_ERI_D2z_Px_S_S_b = I_ERI_Fx2z_S_S_S_b+ABX*I_ERI_D2z_S_S_S_b;
  Double I_ERI_D2x_Py_S_S_b = I_ERI_F2xy_S_S_S_b+ABY*I_ERI_D2x_S_S_S_b;
  Double I_ERI_Dxy_Py_S_S_b = I_ERI_Fx2y_S_S_S_b+ABY*I_ERI_Dxy_S_S_S_b;
  Double I_ERI_Dxz_Py_S_S_b = I_ERI_Fxyz_S_S_S_b+ABY*I_ERI_Dxz_S_S_S_b;
  Double I_ERI_D2y_Py_S_S_b = I_ERI_F3y_S_S_S_b+ABY*I_ERI_D2y_S_S_S_b;
  Double I_ERI_Dyz_Py_S_S_b = I_ERI_F2yz_S_S_S_b+ABY*I_ERI_Dyz_S_S_S_b;
  Double I_ERI_D2z_Py_S_S_b = I_ERI_Fy2z_S_S_S_b+ABY*I_ERI_D2z_S_S_S_b;
  Double I_ERI_D2x_Pz_S_S_b = I_ERI_F2xz_S_S_S_b+ABZ*I_ERI_D2x_S_S_S_b;
  Double I_ERI_Dxy_Pz_S_S_b = I_ERI_Fxyz_S_S_S_b+ABZ*I_ERI_Dxy_S_S_S_b;
  Double I_ERI_Dxz_Pz_S_S_b = I_ERI_Fx2z_S_S_S_b+ABZ*I_ERI_Dxz_S_S_S_b;
  Double I_ERI_D2y_Pz_S_S_b = I_ERI_F2yz_S_S_S_b+ABZ*I_ERI_D2y_S_S_S_b;
  Double I_ERI_Dyz_Pz_S_S_b = I_ERI_Fy2z_S_S_S_b+ABZ*I_ERI_Dyz_S_S_S_b;
  Double I_ERI_D2z_Pz_S_S_b = I_ERI_F3z_S_S_S_b+ABZ*I_ERI_D2z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_b
   * RHS shell quartet name: SQ_ERI_D_S_P_S_b
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_b = I_ERI_F3x_S_Px_S_b+ABX*I_ERI_D2x_S_Px_S_b;
  Double I_ERI_Dxy_Px_Px_S_b = I_ERI_F2xy_S_Px_S_b+ABX*I_ERI_Dxy_S_Px_S_b;
  Double I_ERI_Dxz_Px_Px_S_b = I_ERI_F2xz_S_Px_S_b+ABX*I_ERI_Dxz_S_Px_S_b;
  Double I_ERI_D2y_Px_Px_S_b = I_ERI_Fx2y_S_Px_S_b+ABX*I_ERI_D2y_S_Px_S_b;
  Double I_ERI_Dyz_Px_Px_S_b = I_ERI_Fxyz_S_Px_S_b+ABX*I_ERI_Dyz_S_Px_S_b;
  Double I_ERI_D2z_Px_Px_S_b = I_ERI_Fx2z_S_Px_S_b+ABX*I_ERI_D2z_S_Px_S_b;
  Double I_ERI_D2x_Py_Px_S_b = I_ERI_F2xy_S_Px_S_b+ABY*I_ERI_D2x_S_Px_S_b;
  Double I_ERI_Dxy_Py_Px_S_b = I_ERI_Fx2y_S_Px_S_b+ABY*I_ERI_Dxy_S_Px_S_b;
  Double I_ERI_Dxz_Py_Px_S_b = I_ERI_Fxyz_S_Px_S_b+ABY*I_ERI_Dxz_S_Px_S_b;
  Double I_ERI_D2y_Py_Px_S_b = I_ERI_F3y_S_Px_S_b+ABY*I_ERI_D2y_S_Px_S_b;
  Double I_ERI_Dyz_Py_Px_S_b = I_ERI_F2yz_S_Px_S_b+ABY*I_ERI_Dyz_S_Px_S_b;
  Double I_ERI_D2z_Py_Px_S_b = I_ERI_Fy2z_S_Px_S_b+ABY*I_ERI_D2z_S_Px_S_b;
  Double I_ERI_D2x_Pz_Px_S_b = I_ERI_F2xz_S_Px_S_b+ABZ*I_ERI_D2x_S_Px_S_b;
  Double I_ERI_Dxy_Pz_Px_S_b = I_ERI_Fxyz_S_Px_S_b+ABZ*I_ERI_Dxy_S_Px_S_b;
  Double I_ERI_Dxz_Pz_Px_S_b = I_ERI_Fx2z_S_Px_S_b+ABZ*I_ERI_Dxz_S_Px_S_b;
  Double I_ERI_D2y_Pz_Px_S_b = I_ERI_F2yz_S_Px_S_b+ABZ*I_ERI_D2y_S_Px_S_b;
  Double I_ERI_Dyz_Pz_Px_S_b = I_ERI_Fy2z_S_Px_S_b+ABZ*I_ERI_Dyz_S_Px_S_b;
  Double I_ERI_D2z_Pz_Px_S_b = I_ERI_F3z_S_Px_S_b+ABZ*I_ERI_D2z_S_Px_S_b;
  Double I_ERI_D2x_Px_Py_S_b = I_ERI_F3x_S_Py_S_b+ABX*I_ERI_D2x_S_Py_S_b;
  Double I_ERI_Dxy_Px_Py_S_b = I_ERI_F2xy_S_Py_S_b+ABX*I_ERI_Dxy_S_Py_S_b;
  Double I_ERI_Dxz_Px_Py_S_b = I_ERI_F2xz_S_Py_S_b+ABX*I_ERI_Dxz_S_Py_S_b;
  Double I_ERI_D2y_Px_Py_S_b = I_ERI_Fx2y_S_Py_S_b+ABX*I_ERI_D2y_S_Py_S_b;
  Double I_ERI_Dyz_Px_Py_S_b = I_ERI_Fxyz_S_Py_S_b+ABX*I_ERI_Dyz_S_Py_S_b;
  Double I_ERI_D2z_Px_Py_S_b = I_ERI_Fx2z_S_Py_S_b+ABX*I_ERI_D2z_S_Py_S_b;
  Double I_ERI_D2x_Py_Py_S_b = I_ERI_F2xy_S_Py_S_b+ABY*I_ERI_D2x_S_Py_S_b;
  Double I_ERI_Dxy_Py_Py_S_b = I_ERI_Fx2y_S_Py_S_b+ABY*I_ERI_Dxy_S_Py_S_b;
  Double I_ERI_Dxz_Py_Py_S_b = I_ERI_Fxyz_S_Py_S_b+ABY*I_ERI_Dxz_S_Py_S_b;
  Double I_ERI_D2y_Py_Py_S_b = I_ERI_F3y_S_Py_S_b+ABY*I_ERI_D2y_S_Py_S_b;
  Double I_ERI_Dyz_Py_Py_S_b = I_ERI_F2yz_S_Py_S_b+ABY*I_ERI_Dyz_S_Py_S_b;
  Double I_ERI_D2z_Py_Py_S_b = I_ERI_Fy2z_S_Py_S_b+ABY*I_ERI_D2z_S_Py_S_b;
  Double I_ERI_D2x_Pz_Py_S_b = I_ERI_F2xz_S_Py_S_b+ABZ*I_ERI_D2x_S_Py_S_b;
  Double I_ERI_Dxy_Pz_Py_S_b = I_ERI_Fxyz_S_Py_S_b+ABZ*I_ERI_Dxy_S_Py_S_b;
  Double I_ERI_Dxz_Pz_Py_S_b = I_ERI_Fx2z_S_Py_S_b+ABZ*I_ERI_Dxz_S_Py_S_b;
  Double I_ERI_D2y_Pz_Py_S_b = I_ERI_F2yz_S_Py_S_b+ABZ*I_ERI_D2y_S_Py_S_b;
  Double I_ERI_Dyz_Pz_Py_S_b = I_ERI_Fy2z_S_Py_S_b+ABZ*I_ERI_Dyz_S_Py_S_b;
  Double I_ERI_D2z_Pz_Py_S_b = I_ERI_F3z_S_Py_S_b+ABZ*I_ERI_D2z_S_Py_S_b;
  Double I_ERI_D2x_Px_Pz_S_b = I_ERI_F3x_S_Pz_S_b+ABX*I_ERI_D2x_S_Pz_S_b;
  Double I_ERI_Dxy_Px_Pz_S_b = I_ERI_F2xy_S_Pz_S_b+ABX*I_ERI_Dxy_S_Pz_S_b;
  Double I_ERI_Dxz_Px_Pz_S_b = I_ERI_F2xz_S_Pz_S_b+ABX*I_ERI_Dxz_S_Pz_S_b;
  Double I_ERI_D2y_Px_Pz_S_b = I_ERI_Fx2y_S_Pz_S_b+ABX*I_ERI_D2y_S_Pz_S_b;
  Double I_ERI_Dyz_Px_Pz_S_b = I_ERI_Fxyz_S_Pz_S_b+ABX*I_ERI_Dyz_S_Pz_S_b;
  Double I_ERI_D2z_Px_Pz_S_b = I_ERI_Fx2z_S_Pz_S_b+ABX*I_ERI_D2z_S_Pz_S_b;
  Double I_ERI_D2x_Py_Pz_S_b = I_ERI_F2xy_S_Pz_S_b+ABY*I_ERI_D2x_S_Pz_S_b;
  Double I_ERI_Dxy_Py_Pz_S_b = I_ERI_Fx2y_S_Pz_S_b+ABY*I_ERI_Dxy_S_Pz_S_b;
  Double I_ERI_Dxz_Py_Pz_S_b = I_ERI_Fxyz_S_Pz_S_b+ABY*I_ERI_Dxz_S_Pz_S_b;
  Double I_ERI_D2y_Py_Pz_S_b = I_ERI_F3y_S_Pz_S_b+ABY*I_ERI_D2y_S_Pz_S_b;
  Double I_ERI_Dyz_Py_Pz_S_b = I_ERI_F2yz_S_Pz_S_b+ABY*I_ERI_Dyz_S_Pz_S_b;
  Double I_ERI_D2z_Py_Pz_S_b = I_ERI_Fy2z_S_Pz_S_b+ABY*I_ERI_D2z_S_Pz_S_b;
  Double I_ERI_D2x_Pz_Pz_S_b = I_ERI_F2xz_S_Pz_S_b+ABZ*I_ERI_D2x_S_Pz_S_b;
  Double I_ERI_Dxy_Pz_Pz_S_b = I_ERI_Fxyz_S_Pz_S_b+ABZ*I_ERI_Dxy_S_Pz_S_b;
  Double I_ERI_Dxz_Pz_Pz_S_b = I_ERI_Fx2z_S_Pz_S_b+ABZ*I_ERI_Dxz_S_Pz_S_b;
  Double I_ERI_D2y_Pz_Pz_S_b = I_ERI_F2yz_S_Pz_S_b+ABZ*I_ERI_D2y_S_Pz_S_b;
  Double I_ERI_Dyz_Pz_Pz_S_b = I_ERI_Fy2z_S_Pz_S_b+ABZ*I_ERI_Dyz_S_Pz_S_b;
  Double I_ERI_D2z_Pz_Pz_S_b = I_ERI_F3z_S_Pz_S_b+ABZ*I_ERI_D2z_S_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_P_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   * RHS shell quartet name: SQ_ERI_P_P_P_S_b
   ************************************************************/
  Double I_ERI_Px_D2x_Px_S_b = I_ERI_D2x_Px_Px_S_b+ABX*I_ERI_Px_Px_Px_S_b;
  Double I_ERI_Py_D2x_Px_S_b = I_ERI_Dxy_Px_Px_S_b+ABX*I_ERI_Py_Px_Px_S_b;
  Double I_ERI_Pz_D2x_Px_S_b = I_ERI_Dxz_Px_Px_S_b+ABX*I_ERI_Pz_Px_Px_S_b;
  Double I_ERI_Px_Dxy_Px_S_b = I_ERI_Dxy_Px_Px_S_b+ABY*I_ERI_Px_Px_Px_S_b;
  Double I_ERI_Py_Dxy_Px_S_b = I_ERI_D2y_Px_Px_S_b+ABY*I_ERI_Py_Px_Px_S_b;
  Double I_ERI_Pz_Dxy_Px_S_b = I_ERI_Dyz_Px_Px_S_b+ABY*I_ERI_Pz_Px_Px_S_b;
  Double I_ERI_Px_Dxz_Px_S_b = I_ERI_Dxz_Px_Px_S_b+ABZ*I_ERI_Px_Px_Px_S_b;
  Double I_ERI_Py_Dxz_Px_S_b = I_ERI_Dyz_Px_Px_S_b+ABZ*I_ERI_Py_Px_Px_S_b;
  Double I_ERI_Pz_Dxz_Px_S_b = I_ERI_D2z_Px_Px_S_b+ABZ*I_ERI_Pz_Px_Px_S_b;
  Double I_ERI_Px_D2y_Px_S_b = I_ERI_Dxy_Py_Px_S_b+ABY*I_ERI_Px_Py_Px_S_b;
  Double I_ERI_Py_D2y_Px_S_b = I_ERI_D2y_Py_Px_S_b+ABY*I_ERI_Py_Py_Px_S_b;
  Double I_ERI_Pz_D2y_Px_S_b = I_ERI_Dyz_Py_Px_S_b+ABY*I_ERI_Pz_Py_Px_S_b;
  Double I_ERI_Px_Dyz_Px_S_b = I_ERI_Dxz_Py_Px_S_b+ABZ*I_ERI_Px_Py_Px_S_b;
  Double I_ERI_Py_Dyz_Px_S_b = I_ERI_Dyz_Py_Px_S_b+ABZ*I_ERI_Py_Py_Px_S_b;
  Double I_ERI_Pz_Dyz_Px_S_b = I_ERI_D2z_Py_Px_S_b+ABZ*I_ERI_Pz_Py_Px_S_b;
  Double I_ERI_Px_D2z_Px_S_b = I_ERI_Dxz_Pz_Px_S_b+ABZ*I_ERI_Px_Pz_Px_S_b;
  Double I_ERI_Py_D2z_Px_S_b = I_ERI_Dyz_Pz_Px_S_b+ABZ*I_ERI_Py_Pz_Px_S_b;
  Double I_ERI_Pz_D2z_Px_S_b = I_ERI_D2z_Pz_Px_S_b+ABZ*I_ERI_Pz_Pz_Px_S_b;
  Double I_ERI_Px_D2x_Py_S_b = I_ERI_D2x_Px_Py_S_b+ABX*I_ERI_Px_Px_Py_S_b;
  Double I_ERI_Py_D2x_Py_S_b = I_ERI_Dxy_Px_Py_S_b+ABX*I_ERI_Py_Px_Py_S_b;
  Double I_ERI_Pz_D2x_Py_S_b = I_ERI_Dxz_Px_Py_S_b+ABX*I_ERI_Pz_Px_Py_S_b;
  Double I_ERI_Px_Dxy_Py_S_b = I_ERI_Dxy_Px_Py_S_b+ABY*I_ERI_Px_Px_Py_S_b;
  Double I_ERI_Py_Dxy_Py_S_b = I_ERI_D2y_Px_Py_S_b+ABY*I_ERI_Py_Px_Py_S_b;
  Double I_ERI_Pz_Dxy_Py_S_b = I_ERI_Dyz_Px_Py_S_b+ABY*I_ERI_Pz_Px_Py_S_b;
  Double I_ERI_Px_Dxz_Py_S_b = I_ERI_Dxz_Px_Py_S_b+ABZ*I_ERI_Px_Px_Py_S_b;
  Double I_ERI_Py_Dxz_Py_S_b = I_ERI_Dyz_Px_Py_S_b+ABZ*I_ERI_Py_Px_Py_S_b;
  Double I_ERI_Pz_Dxz_Py_S_b = I_ERI_D2z_Px_Py_S_b+ABZ*I_ERI_Pz_Px_Py_S_b;
  Double I_ERI_Px_D2y_Py_S_b = I_ERI_Dxy_Py_Py_S_b+ABY*I_ERI_Px_Py_Py_S_b;
  Double I_ERI_Py_D2y_Py_S_b = I_ERI_D2y_Py_Py_S_b+ABY*I_ERI_Py_Py_Py_S_b;
  Double I_ERI_Pz_D2y_Py_S_b = I_ERI_Dyz_Py_Py_S_b+ABY*I_ERI_Pz_Py_Py_S_b;
  Double I_ERI_Px_Dyz_Py_S_b = I_ERI_Dxz_Py_Py_S_b+ABZ*I_ERI_Px_Py_Py_S_b;
  Double I_ERI_Py_Dyz_Py_S_b = I_ERI_Dyz_Py_Py_S_b+ABZ*I_ERI_Py_Py_Py_S_b;
  Double I_ERI_Pz_Dyz_Py_S_b = I_ERI_D2z_Py_Py_S_b+ABZ*I_ERI_Pz_Py_Py_S_b;
  Double I_ERI_Px_D2z_Py_S_b = I_ERI_Dxz_Pz_Py_S_b+ABZ*I_ERI_Px_Pz_Py_S_b;
  Double I_ERI_Py_D2z_Py_S_b = I_ERI_Dyz_Pz_Py_S_b+ABZ*I_ERI_Py_Pz_Py_S_b;
  Double I_ERI_Pz_D2z_Py_S_b = I_ERI_D2z_Pz_Py_S_b+ABZ*I_ERI_Pz_Pz_Py_S_b;
  Double I_ERI_Px_D2x_Pz_S_b = I_ERI_D2x_Px_Pz_S_b+ABX*I_ERI_Px_Px_Pz_S_b;
  Double I_ERI_Py_D2x_Pz_S_b = I_ERI_Dxy_Px_Pz_S_b+ABX*I_ERI_Py_Px_Pz_S_b;
  Double I_ERI_Pz_D2x_Pz_S_b = I_ERI_Dxz_Px_Pz_S_b+ABX*I_ERI_Pz_Px_Pz_S_b;
  Double I_ERI_Px_Dxy_Pz_S_b = I_ERI_Dxy_Px_Pz_S_b+ABY*I_ERI_Px_Px_Pz_S_b;
  Double I_ERI_Py_Dxy_Pz_S_b = I_ERI_D2y_Px_Pz_S_b+ABY*I_ERI_Py_Px_Pz_S_b;
  Double I_ERI_Pz_Dxy_Pz_S_b = I_ERI_Dyz_Px_Pz_S_b+ABY*I_ERI_Pz_Px_Pz_S_b;
  Double I_ERI_Px_Dxz_Pz_S_b = I_ERI_Dxz_Px_Pz_S_b+ABZ*I_ERI_Px_Px_Pz_S_b;
  Double I_ERI_Py_Dxz_Pz_S_b = I_ERI_Dyz_Px_Pz_S_b+ABZ*I_ERI_Py_Px_Pz_S_b;
  Double I_ERI_Pz_Dxz_Pz_S_b = I_ERI_D2z_Px_Pz_S_b+ABZ*I_ERI_Pz_Px_Pz_S_b;
  Double I_ERI_Px_D2y_Pz_S_b = I_ERI_Dxy_Py_Pz_S_b+ABY*I_ERI_Px_Py_Pz_S_b;
  Double I_ERI_Py_D2y_Pz_S_b = I_ERI_D2y_Py_Pz_S_b+ABY*I_ERI_Py_Py_Pz_S_b;
  Double I_ERI_Pz_D2y_Pz_S_b = I_ERI_Dyz_Py_Pz_S_b+ABY*I_ERI_Pz_Py_Pz_S_b;
  Double I_ERI_Px_Dyz_Pz_S_b = I_ERI_Dxz_Py_Pz_S_b+ABZ*I_ERI_Px_Py_Pz_S_b;
  Double I_ERI_Py_Dyz_Pz_S_b = I_ERI_Dyz_Py_Pz_S_b+ABZ*I_ERI_Py_Py_Pz_S_b;
  Double I_ERI_Pz_Dyz_Pz_S_b = I_ERI_D2z_Py_Pz_S_b+ABZ*I_ERI_Pz_Py_Pz_S_b;
  Double I_ERI_Px_D2z_Pz_S_b = I_ERI_Dxz_Pz_Pz_S_b+ABZ*I_ERI_Px_Pz_Pz_S_b;
  Double I_ERI_Py_D2z_Pz_S_b = I_ERI_Dyz_Pz_Pz_S_b+ABZ*I_ERI_Py_Pz_Pz_S_b;
  Double I_ERI_Pz_D2z_Pz_S_b = I_ERI_D2z_Pz_Pz_S_b+ABZ*I_ERI_Pz_Pz_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 5 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S_b
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_b = I_ERI_G4x_S_S_S_b+ABX*I_ERI_F3x_S_S_S_b;
  Double I_ERI_F2xy_Px_S_S_b = I_ERI_G3xy_S_S_S_b+ABX*I_ERI_F2xy_S_S_S_b;
  Double I_ERI_F2xz_Px_S_S_b = I_ERI_G3xz_S_S_S_b+ABX*I_ERI_F2xz_S_S_S_b;
  Double I_ERI_Fx2y_Px_S_S_b = I_ERI_G2x2y_S_S_S_b+ABX*I_ERI_Fx2y_S_S_S_b;
  Double I_ERI_Fxyz_Px_S_S_b = I_ERI_G2xyz_S_S_S_b+ABX*I_ERI_Fxyz_S_S_S_b;
  Double I_ERI_Fx2z_Px_S_S_b = I_ERI_G2x2z_S_S_S_b+ABX*I_ERI_Fx2z_S_S_S_b;
  Double I_ERI_F3y_Px_S_S_b = I_ERI_Gx3y_S_S_S_b+ABX*I_ERI_F3y_S_S_S_b;
  Double I_ERI_F2yz_Px_S_S_b = I_ERI_Gx2yz_S_S_S_b+ABX*I_ERI_F2yz_S_S_S_b;
  Double I_ERI_Fy2z_Px_S_S_b = I_ERI_Gxy2z_S_S_S_b+ABX*I_ERI_Fy2z_S_S_S_b;
  Double I_ERI_F3z_Px_S_S_b = I_ERI_Gx3z_S_S_S_b+ABX*I_ERI_F3z_S_S_S_b;
  Double I_ERI_F2xy_Py_S_S_b = I_ERI_G2x2y_S_S_S_b+ABY*I_ERI_F2xy_S_S_S_b;
  Double I_ERI_F2xz_Py_S_S_b = I_ERI_G2xyz_S_S_S_b+ABY*I_ERI_F2xz_S_S_S_b;
  Double I_ERI_Fx2y_Py_S_S_b = I_ERI_Gx3y_S_S_S_b+ABY*I_ERI_Fx2y_S_S_S_b;
  Double I_ERI_Fxyz_Py_S_S_b = I_ERI_Gx2yz_S_S_S_b+ABY*I_ERI_Fxyz_S_S_S_b;
  Double I_ERI_Fx2z_Py_S_S_b = I_ERI_Gxy2z_S_S_S_b+ABY*I_ERI_Fx2z_S_S_S_b;
  Double I_ERI_F3y_Py_S_S_b = I_ERI_G4y_S_S_S_b+ABY*I_ERI_F3y_S_S_S_b;
  Double I_ERI_F2yz_Py_S_S_b = I_ERI_G3yz_S_S_S_b+ABY*I_ERI_F2yz_S_S_S_b;
  Double I_ERI_Fy2z_Py_S_S_b = I_ERI_G2y2z_S_S_S_b+ABY*I_ERI_Fy2z_S_S_S_b;
  Double I_ERI_F3z_Py_S_S_b = I_ERI_Gy3z_S_S_S_b+ABY*I_ERI_F3z_S_S_S_b;
  Double I_ERI_F2xz_Pz_S_S_b = I_ERI_G2x2z_S_S_S_b+ABZ*I_ERI_F2xz_S_S_S_b;
  Double I_ERI_Fxyz_Pz_S_S_b = I_ERI_Gxy2z_S_S_S_b+ABZ*I_ERI_Fxyz_S_S_S_b;
  Double I_ERI_Fx2z_Pz_S_S_b = I_ERI_Gx3z_S_S_S_b+ABZ*I_ERI_Fx2z_S_S_S_b;
  Double I_ERI_F2yz_Pz_S_S_b = I_ERI_G2y2z_S_S_S_b+ABZ*I_ERI_F2yz_S_S_S_b;
  Double I_ERI_Fy2z_Pz_S_S_b = I_ERI_Gy3z_S_S_S_b+ABZ*I_ERI_Fy2z_S_S_S_b;
  Double I_ERI_F3z_Pz_S_S_b = I_ERI_G4z_S_S_S_b+ABZ*I_ERI_F3z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S_b
   ************************************************************/
  Double I_ERI_D2x_D2x_S_S_b = I_ERI_F3x_Px_S_S_b+ABX*I_ERI_D2x_Px_S_S_b;
  Double I_ERI_Dxy_D2x_S_S_b = I_ERI_F2xy_Px_S_S_b+ABX*I_ERI_Dxy_Px_S_S_b;
  Double I_ERI_Dxz_D2x_S_S_b = I_ERI_F2xz_Px_S_S_b+ABX*I_ERI_Dxz_Px_S_S_b;
  Double I_ERI_D2y_D2x_S_S_b = I_ERI_Fx2y_Px_S_S_b+ABX*I_ERI_D2y_Px_S_S_b;
  Double I_ERI_Dyz_D2x_S_S_b = I_ERI_Fxyz_Px_S_S_b+ABX*I_ERI_Dyz_Px_S_S_b;
  Double I_ERI_D2z_D2x_S_S_b = I_ERI_Fx2z_Px_S_S_b+ABX*I_ERI_D2z_Px_S_S_b;
  Double I_ERI_D2x_Dxy_S_S_b = I_ERI_F2xy_Px_S_S_b+ABY*I_ERI_D2x_Px_S_S_b;
  Double I_ERI_Dxy_Dxy_S_S_b = I_ERI_Fx2y_Px_S_S_b+ABY*I_ERI_Dxy_Px_S_S_b;
  Double I_ERI_Dxz_Dxy_S_S_b = I_ERI_Fxyz_Px_S_S_b+ABY*I_ERI_Dxz_Px_S_S_b;
  Double I_ERI_D2y_Dxy_S_S_b = I_ERI_F3y_Px_S_S_b+ABY*I_ERI_D2y_Px_S_S_b;
  Double I_ERI_Dyz_Dxy_S_S_b = I_ERI_F2yz_Px_S_S_b+ABY*I_ERI_Dyz_Px_S_S_b;
  Double I_ERI_D2z_Dxy_S_S_b = I_ERI_Fy2z_Px_S_S_b+ABY*I_ERI_D2z_Px_S_S_b;
  Double I_ERI_D2x_Dxz_S_S_b = I_ERI_F2xz_Px_S_S_b+ABZ*I_ERI_D2x_Px_S_S_b;
  Double I_ERI_Dxy_Dxz_S_S_b = I_ERI_Fxyz_Px_S_S_b+ABZ*I_ERI_Dxy_Px_S_S_b;
  Double I_ERI_Dxz_Dxz_S_S_b = I_ERI_Fx2z_Px_S_S_b+ABZ*I_ERI_Dxz_Px_S_S_b;
  Double I_ERI_D2y_Dxz_S_S_b = I_ERI_F2yz_Px_S_S_b+ABZ*I_ERI_D2y_Px_S_S_b;
  Double I_ERI_Dyz_Dxz_S_S_b = I_ERI_Fy2z_Px_S_S_b+ABZ*I_ERI_Dyz_Px_S_S_b;
  Double I_ERI_D2z_Dxz_S_S_b = I_ERI_F3z_Px_S_S_b+ABZ*I_ERI_D2z_Px_S_S_b;
  Double I_ERI_D2x_D2y_S_S_b = I_ERI_F2xy_Py_S_S_b+ABY*I_ERI_D2x_Py_S_S_b;
  Double I_ERI_Dxy_D2y_S_S_b = I_ERI_Fx2y_Py_S_S_b+ABY*I_ERI_Dxy_Py_S_S_b;
  Double I_ERI_Dxz_D2y_S_S_b = I_ERI_Fxyz_Py_S_S_b+ABY*I_ERI_Dxz_Py_S_S_b;
  Double I_ERI_D2y_D2y_S_S_b = I_ERI_F3y_Py_S_S_b+ABY*I_ERI_D2y_Py_S_S_b;
  Double I_ERI_Dyz_D2y_S_S_b = I_ERI_F2yz_Py_S_S_b+ABY*I_ERI_Dyz_Py_S_S_b;
  Double I_ERI_D2z_D2y_S_S_b = I_ERI_Fy2z_Py_S_S_b+ABY*I_ERI_D2z_Py_S_S_b;
  Double I_ERI_D2x_Dyz_S_S_b = I_ERI_F2xz_Py_S_S_b+ABZ*I_ERI_D2x_Py_S_S_b;
  Double I_ERI_Dxy_Dyz_S_S_b = I_ERI_Fxyz_Py_S_S_b+ABZ*I_ERI_Dxy_Py_S_S_b;
  Double I_ERI_Dxz_Dyz_S_S_b = I_ERI_Fx2z_Py_S_S_b+ABZ*I_ERI_Dxz_Py_S_S_b;
  Double I_ERI_D2y_Dyz_S_S_b = I_ERI_F2yz_Py_S_S_b+ABZ*I_ERI_D2y_Py_S_S_b;
  Double I_ERI_Dyz_Dyz_S_S_b = I_ERI_Fy2z_Py_S_S_b+ABZ*I_ERI_Dyz_Py_S_S_b;
  Double I_ERI_D2z_Dyz_S_S_b = I_ERI_F3z_Py_S_S_b+ABZ*I_ERI_D2z_Py_S_S_b;
  Double I_ERI_D2x_D2z_S_S_b = I_ERI_F2xz_Pz_S_S_b+ABZ*I_ERI_D2x_Pz_S_S_b;
  Double I_ERI_Dxy_D2z_S_S_b = I_ERI_Fxyz_Pz_S_S_b+ABZ*I_ERI_Dxy_Pz_S_S_b;
  Double I_ERI_Dxz_D2z_S_S_b = I_ERI_Fx2z_Pz_S_S_b+ABZ*I_ERI_Dxz_Pz_S_S_b;
  Double I_ERI_D2y_D2z_S_S_b = I_ERI_F2yz_Pz_S_S_b+ABZ*I_ERI_D2y_Pz_S_S_b;
  Double I_ERI_Dyz_D2z_S_S_b = I_ERI_Fy2z_Pz_S_S_b+ABZ*I_ERI_Dyz_Pz_S_S_b;
  Double I_ERI_D2z_D2z_S_S_b = I_ERI_F3z_Pz_S_S_b+ABZ*I_ERI_D2z_Pz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_P_P_D_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_P_S_D_S_c
   ************************************************************/
  Double I_ERI_Px_Px_D2x_S_c = I_ERI_D2x_S_D2x_S_c+ABX*I_ERI_Px_S_D2x_S_c;
  Double I_ERI_Py_Px_D2x_S_c = I_ERI_Dxy_S_D2x_S_c+ABX*I_ERI_Py_S_D2x_S_c;
  Double I_ERI_Pz_Px_D2x_S_c = I_ERI_Dxz_S_D2x_S_c+ABX*I_ERI_Pz_S_D2x_S_c;
  Double I_ERI_Px_Py_D2x_S_c = I_ERI_Dxy_S_D2x_S_c+ABY*I_ERI_Px_S_D2x_S_c;
  Double I_ERI_Py_Py_D2x_S_c = I_ERI_D2y_S_D2x_S_c+ABY*I_ERI_Py_S_D2x_S_c;
  Double I_ERI_Pz_Py_D2x_S_c = I_ERI_Dyz_S_D2x_S_c+ABY*I_ERI_Pz_S_D2x_S_c;
  Double I_ERI_Px_Pz_D2x_S_c = I_ERI_Dxz_S_D2x_S_c+ABZ*I_ERI_Px_S_D2x_S_c;
  Double I_ERI_Py_Pz_D2x_S_c = I_ERI_Dyz_S_D2x_S_c+ABZ*I_ERI_Py_S_D2x_S_c;
  Double I_ERI_Pz_Pz_D2x_S_c = I_ERI_D2z_S_D2x_S_c+ABZ*I_ERI_Pz_S_D2x_S_c;
  Double I_ERI_Px_Px_Dxy_S_c = I_ERI_D2x_S_Dxy_S_c+ABX*I_ERI_Px_S_Dxy_S_c;
  Double I_ERI_Py_Px_Dxy_S_c = I_ERI_Dxy_S_Dxy_S_c+ABX*I_ERI_Py_S_Dxy_S_c;
  Double I_ERI_Pz_Px_Dxy_S_c = I_ERI_Dxz_S_Dxy_S_c+ABX*I_ERI_Pz_S_Dxy_S_c;
  Double I_ERI_Px_Py_Dxy_S_c = I_ERI_Dxy_S_Dxy_S_c+ABY*I_ERI_Px_S_Dxy_S_c;
  Double I_ERI_Py_Py_Dxy_S_c = I_ERI_D2y_S_Dxy_S_c+ABY*I_ERI_Py_S_Dxy_S_c;
  Double I_ERI_Pz_Py_Dxy_S_c = I_ERI_Dyz_S_Dxy_S_c+ABY*I_ERI_Pz_S_Dxy_S_c;
  Double I_ERI_Px_Pz_Dxy_S_c = I_ERI_Dxz_S_Dxy_S_c+ABZ*I_ERI_Px_S_Dxy_S_c;
  Double I_ERI_Py_Pz_Dxy_S_c = I_ERI_Dyz_S_Dxy_S_c+ABZ*I_ERI_Py_S_Dxy_S_c;
  Double I_ERI_Pz_Pz_Dxy_S_c = I_ERI_D2z_S_Dxy_S_c+ABZ*I_ERI_Pz_S_Dxy_S_c;
  Double I_ERI_Px_Px_Dxz_S_c = I_ERI_D2x_S_Dxz_S_c+ABX*I_ERI_Px_S_Dxz_S_c;
  Double I_ERI_Py_Px_Dxz_S_c = I_ERI_Dxy_S_Dxz_S_c+ABX*I_ERI_Py_S_Dxz_S_c;
  Double I_ERI_Pz_Px_Dxz_S_c = I_ERI_Dxz_S_Dxz_S_c+ABX*I_ERI_Pz_S_Dxz_S_c;
  Double I_ERI_Px_Py_Dxz_S_c = I_ERI_Dxy_S_Dxz_S_c+ABY*I_ERI_Px_S_Dxz_S_c;
  Double I_ERI_Py_Py_Dxz_S_c = I_ERI_D2y_S_Dxz_S_c+ABY*I_ERI_Py_S_Dxz_S_c;
  Double I_ERI_Pz_Py_Dxz_S_c = I_ERI_Dyz_S_Dxz_S_c+ABY*I_ERI_Pz_S_Dxz_S_c;
  Double I_ERI_Px_Pz_Dxz_S_c = I_ERI_Dxz_S_Dxz_S_c+ABZ*I_ERI_Px_S_Dxz_S_c;
  Double I_ERI_Py_Pz_Dxz_S_c = I_ERI_Dyz_S_Dxz_S_c+ABZ*I_ERI_Py_S_Dxz_S_c;
  Double I_ERI_Pz_Pz_Dxz_S_c = I_ERI_D2z_S_Dxz_S_c+ABZ*I_ERI_Pz_S_Dxz_S_c;
  Double I_ERI_Px_Px_D2y_S_c = I_ERI_D2x_S_D2y_S_c+ABX*I_ERI_Px_S_D2y_S_c;
  Double I_ERI_Py_Px_D2y_S_c = I_ERI_Dxy_S_D2y_S_c+ABX*I_ERI_Py_S_D2y_S_c;
  Double I_ERI_Pz_Px_D2y_S_c = I_ERI_Dxz_S_D2y_S_c+ABX*I_ERI_Pz_S_D2y_S_c;
  Double I_ERI_Px_Py_D2y_S_c = I_ERI_Dxy_S_D2y_S_c+ABY*I_ERI_Px_S_D2y_S_c;
  Double I_ERI_Py_Py_D2y_S_c = I_ERI_D2y_S_D2y_S_c+ABY*I_ERI_Py_S_D2y_S_c;
  Double I_ERI_Pz_Py_D2y_S_c = I_ERI_Dyz_S_D2y_S_c+ABY*I_ERI_Pz_S_D2y_S_c;
  Double I_ERI_Px_Pz_D2y_S_c = I_ERI_Dxz_S_D2y_S_c+ABZ*I_ERI_Px_S_D2y_S_c;
  Double I_ERI_Py_Pz_D2y_S_c = I_ERI_Dyz_S_D2y_S_c+ABZ*I_ERI_Py_S_D2y_S_c;
  Double I_ERI_Pz_Pz_D2y_S_c = I_ERI_D2z_S_D2y_S_c+ABZ*I_ERI_Pz_S_D2y_S_c;
  Double I_ERI_Px_Px_Dyz_S_c = I_ERI_D2x_S_Dyz_S_c+ABX*I_ERI_Px_S_Dyz_S_c;
  Double I_ERI_Py_Px_Dyz_S_c = I_ERI_Dxy_S_Dyz_S_c+ABX*I_ERI_Py_S_Dyz_S_c;
  Double I_ERI_Pz_Px_Dyz_S_c = I_ERI_Dxz_S_Dyz_S_c+ABX*I_ERI_Pz_S_Dyz_S_c;
  Double I_ERI_Px_Py_Dyz_S_c = I_ERI_Dxy_S_Dyz_S_c+ABY*I_ERI_Px_S_Dyz_S_c;
  Double I_ERI_Py_Py_Dyz_S_c = I_ERI_D2y_S_Dyz_S_c+ABY*I_ERI_Py_S_Dyz_S_c;
  Double I_ERI_Pz_Py_Dyz_S_c = I_ERI_Dyz_S_Dyz_S_c+ABY*I_ERI_Pz_S_Dyz_S_c;
  Double I_ERI_Px_Pz_Dyz_S_c = I_ERI_Dxz_S_Dyz_S_c+ABZ*I_ERI_Px_S_Dyz_S_c;
  Double I_ERI_Py_Pz_Dyz_S_c = I_ERI_Dyz_S_Dyz_S_c+ABZ*I_ERI_Py_S_Dyz_S_c;
  Double I_ERI_Pz_Pz_Dyz_S_c = I_ERI_D2z_S_Dyz_S_c+ABZ*I_ERI_Pz_S_Dyz_S_c;
  Double I_ERI_Px_Px_D2z_S_c = I_ERI_D2x_S_D2z_S_c+ABX*I_ERI_Px_S_D2z_S_c;
  Double I_ERI_Py_Px_D2z_S_c = I_ERI_Dxy_S_D2z_S_c+ABX*I_ERI_Py_S_D2z_S_c;
  Double I_ERI_Pz_Px_D2z_S_c = I_ERI_Dxz_S_D2z_S_c+ABX*I_ERI_Pz_S_D2z_S_c;
  Double I_ERI_Px_Py_D2z_S_c = I_ERI_Dxy_S_D2z_S_c+ABY*I_ERI_Px_S_D2z_S_c;
  Double I_ERI_Py_Py_D2z_S_c = I_ERI_D2y_S_D2z_S_c+ABY*I_ERI_Py_S_D2z_S_c;
  Double I_ERI_Pz_Py_D2z_S_c = I_ERI_Dyz_S_D2z_S_c+ABY*I_ERI_Pz_S_D2z_S_c;
  Double I_ERI_Px_Pz_D2z_S_c = I_ERI_Dxz_S_D2z_S_c+ABZ*I_ERI_Px_S_D2z_S_c;
  Double I_ERI_Py_Pz_D2z_S_c = I_ERI_Dyz_S_D2z_S_c+ABZ*I_ERI_Py_S_D2z_S_c;
  Double I_ERI_Pz_Pz_D2z_S_c = I_ERI_D2z_S_D2z_S_c+ABZ*I_ERI_Pz_S_D2z_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   * RHS shell quartet name: SQ_ERI_D_S_P_S_c
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_c = I_ERI_F3x_S_Px_S_c+ABX*I_ERI_D2x_S_Px_S_c;
  Double I_ERI_Dxy_Px_Px_S_c = I_ERI_F2xy_S_Px_S_c+ABX*I_ERI_Dxy_S_Px_S_c;
  Double I_ERI_Dxz_Px_Px_S_c = I_ERI_F2xz_S_Px_S_c+ABX*I_ERI_Dxz_S_Px_S_c;
  Double I_ERI_D2y_Px_Px_S_c = I_ERI_Fx2y_S_Px_S_c+ABX*I_ERI_D2y_S_Px_S_c;
  Double I_ERI_Dyz_Px_Px_S_c = I_ERI_Fxyz_S_Px_S_c+ABX*I_ERI_Dyz_S_Px_S_c;
  Double I_ERI_D2z_Px_Px_S_c = I_ERI_Fx2z_S_Px_S_c+ABX*I_ERI_D2z_S_Px_S_c;
  Double I_ERI_D2x_Py_Px_S_c = I_ERI_F2xy_S_Px_S_c+ABY*I_ERI_D2x_S_Px_S_c;
  Double I_ERI_Dxy_Py_Px_S_c = I_ERI_Fx2y_S_Px_S_c+ABY*I_ERI_Dxy_S_Px_S_c;
  Double I_ERI_Dxz_Py_Px_S_c = I_ERI_Fxyz_S_Px_S_c+ABY*I_ERI_Dxz_S_Px_S_c;
  Double I_ERI_D2y_Py_Px_S_c = I_ERI_F3y_S_Px_S_c+ABY*I_ERI_D2y_S_Px_S_c;
  Double I_ERI_Dyz_Py_Px_S_c = I_ERI_F2yz_S_Px_S_c+ABY*I_ERI_Dyz_S_Px_S_c;
  Double I_ERI_D2z_Py_Px_S_c = I_ERI_Fy2z_S_Px_S_c+ABY*I_ERI_D2z_S_Px_S_c;
  Double I_ERI_D2x_Pz_Px_S_c = I_ERI_F2xz_S_Px_S_c+ABZ*I_ERI_D2x_S_Px_S_c;
  Double I_ERI_Dxy_Pz_Px_S_c = I_ERI_Fxyz_S_Px_S_c+ABZ*I_ERI_Dxy_S_Px_S_c;
  Double I_ERI_Dxz_Pz_Px_S_c = I_ERI_Fx2z_S_Px_S_c+ABZ*I_ERI_Dxz_S_Px_S_c;
  Double I_ERI_D2y_Pz_Px_S_c = I_ERI_F2yz_S_Px_S_c+ABZ*I_ERI_D2y_S_Px_S_c;
  Double I_ERI_Dyz_Pz_Px_S_c = I_ERI_Fy2z_S_Px_S_c+ABZ*I_ERI_Dyz_S_Px_S_c;
  Double I_ERI_D2z_Pz_Px_S_c = I_ERI_F3z_S_Px_S_c+ABZ*I_ERI_D2z_S_Px_S_c;
  Double I_ERI_D2x_Px_Py_S_c = I_ERI_F3x_S_Py_S_c+ABX*I_ERI_D2x_S_Py_S_c;
  Double I_ERI_Dxy_Px_Py_S_c = I_ERI_F2xy_S_Py_S_c+ABX*I_ERI_Dxy_S_Py_S_c;
  Double I_ERI_Dxz_Px_Py_S_c = I_ERI_F2xz_S_Py_S_c+ABX*I_ERI_Dxz_S_Py_S_c;
  Double I_ERI_D2y_Px_Py_S_c = I_ERI_Fx2y_S_Py_S_c+ABX*I_ERI_D2y_S_Py_S_c;
  Double I_ERI_Dyz_Px_Py_S_c = I_ERI_Fxyz_S_Py_S_c+ABX*I_ERI_Dyz_S_Py_S_c;
  Double I_ERI_D2z_Px_Py_S_c = I_ERI_Fx2z_S_Py_S_c+ABX*I_ERI_D2z_S_Py_S_c;
  Double I_ERI_D2x_Py_Py_S_c = I_ERI_F2xy_S_Py_S_c+ABY*I_ERI_D2x_S_Py_S_c;
  Double I_ERI_Dxy_Py_Py_S_c = I_ERI_Fx2y_S_Py_S_c+ABY*I_ERI_Dxy_S_Py_S_c;
  Double I_ERI_Dxz_Py_Py_S_c = I_ERI_Fxyz_S_Py_S_c+ABY*I_ERI_Dxz_S_Py_S_c;
  Double I_ERI_D2y_Py_Py_S_c = I_ERI_F3y_S_Py_S_c+ABY*I_ERI_D2y_S_Py_S_c;
  Double I_ERI_Dyz_Py_Py_S_c = I_ERI_F2yz_S_Py_S_c+ABY*I_ERI_Dyz_S_Py_S_c;
  Double I_ERI_D2z_Py_Py_S_c = I_ERI_Fy2z_S_Py_S_c+ABY*I_ERI_D2z_S_Py_S_c;
  Double I_ERI_D2x_Pz_Py_S_c = I_ERI_F2xz_S_Py_S_c+ABZ*I_ERI_D2x_S_Py_S_c;
  Double I_ERI_Dxy_Pz_Py_S_c = I_ERI_Fxyz_S_Py_S_c+ABZ*I_ERI_Dxy_S_Py_S_c;
  Double I_ERI_Dxz_Pz_Py_S_c = I_ERI_Fx2z_S_Py_S_c+ABZ*I_ERI_Dxz_S_Py_S_c;
  Double I_ERI_D2y_Pz_Py_S_c = I_ERI_F2yz_S_Py_S_c+ABZ*I_ERI_D2y_S_Py_S_c;
  Double I_ERI_Dyz_Pz_Py_S_c = I_ERI_Fy2z_S_Py_S_c+ABZ*I_ERI_Dyz_S_Py_S_c;
  Double I_ERI_D2z_Pz_Py_S_c = I_ERI_F3z_S_Py_S_c+ABZ*I_ERI_D2z_S_Py_S_c;
  Double I_ERI_D2x_Px_Pz_S_c = I_ERI_F3x_S_Pz_S_c+ABX*I_ERI_D2x_S_Pz_S_c;
  Double I_ERI_Dxy_Px_Pz_S_c = I_ERI_F2xy_S_Pz_S_c+ABX*I_ERI_Dxy_S_Pz_S_c;
  Double I_ERI_Dxz_Px_Pz_S_c = I_ERI_F2xz_S_Pz_S_c+ABX*I_ERI_Dxz_S_Pz_S_c;
  Double I_ERI_D2y_Px_Pz_S_c = I_ERI_Fx2y_S_Pz_S_c+ABX*I_ERI_D2y_S_Pz_S_c;
  Double I_ERI_Dyz_Px_Pz_S_c = I_ERI_Fxyz_S_Pz_S_c+ABX*I_ERI_Dyz_S_Pz_S_c;
  Double I_ERI_D2z_Px_Pz_S_c = I_ERI_Fx2z_S_Pz_S_c+ABX*I_ERI_D2z_S_Pz_S_c;
  Double I_ERI_D2x_Py_Pz_S_c = I_ERI_F2xy_S_Pz_S_c+ABY*I_ERI_D2x_S_Pz_S_c;
  Double I_ERI_Dxy_Py_Pz_S_c = I_ERI_Fx2y_S_Pz_S_c+ABY*I_ERI_Dxy_S_Pz_S_c;
  Double I_ERI_Dxz_Py_Pz_S_c = I_ERI_Fxyz_S_Pz_S_c+ABY*I_ERI_Dxz_S_Pz_S_c;
  Double I_ERI_D2y_Py_Pz_S_c = I_ERI_F3y_S_Pz_S_c+ABY*I_ERI_D2y_S_Pz_S_c;
  Double I_ERI_Dyz_Py_Pz_S_c = I_ERI_F2yz_S_Pz_S_c+ABY*I_ERI_Dyz_S_Pz_S_c;
  Double I_ERI_D2z_Py_Pz_S_c = I_ERI_Fy2z_S_Pz_S_c+ABY*I_ERI_D2z_S_Pz_S_c;
  Double I_ERI_D2x_Pz_Pz_S_c = I_ERI_F2xz_S_Pz_S_c+ABZ*I_ERI_D2x_S_Pz_S_c;
  Double I_ERI_Dxy_Pz_Pz_S_c = I_ERI_Fxyz_S_Pz_S_c+ABZ*I_ERI_Dxy_S_Pz_S_c;
  Double I_ERI_Dxz_Pz_Pz_S_c = I_ERI_Fx2z_S_Pz_S_c+ABZ*I_ERI_Dxz_S_Pz_S_c;
  Double I_ERI_D2y_Pz_Pz_S_c = I_ERI_F2yz_S_Pz_S_c+ABZ*I_ERI_D2y_S_Pz_S_c;
  Double I_ERI_Dyz_Pz_Pz_S_c = I_ERI_Fy2z_S_Pz_S_c+ABZ*I_ERI_Dyz_S_Pz_S_c;
  Double I_ERI_D2z_Pz_Pz_S_c = I_ERI_F3z_S_Pz_S_c+ABZ*I_ERI_D2z_S_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_P_S_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_P_S_aa
   * RHS shell quartet name: SQ_ERI_G_S_P_S_aa
   ************************************************************/
  Double I_ERI_G4x_Px_Px_S_aa = I_ERI_H5x_S_Px_S_aa+ABX*I_ERI_G4x_S_Px_S_aa;
  Double I_ERI_G3xy_Px_Px_S_aa = I_ERI_H4xy_S_Px_S_aa+ABX*I_ERI_G3xy_S_Px_S_aa;
  Double I_ERI_G3xz_Px_Px_S_aa = I_ERI_H4xz_S_Px_S_aa+ABX*I_ERI_G3xz_S_Px_S_aa;
  Double I_ERI_G2x2y_Px_Px_S_aa = I_ERI_H3x2y_S_Px_S_aa+ABX*I_ERI_G2x2y_S_Px_S_aa;
  Double I_ERI_G2xyz_Px_Px_S_aa = I_ERI_H3xyz_S_Px_S_aa+ABX*I_ERI_G2xyz_S_Px_S_aa;
  Double I_ERI_G2x2z_Px_Px_S_aa = I_ERI_H3x2z_S_Px_S_aa+ABX*I_ERI_G2x2z_S_Px_S_aa;
  Double I_ERI_Gx3y_Px_Px_S_aa = I_ERI_H2x3y_S_Px_S_aa+ABX*I_ERI_Gx3y_S_Px_S_aa;
  Double I_ERI_Gx2yz_Px_Px_S_aa = I_ERI_H2x2yz_S_Px_S_aa+ABX*I_ERI_Gx2yz_S_Px_S_aa;
  Double I_ERI_Gxy2z_Px_Px_S_aa = I_ERI_H2xy2z_S_Px_S_aa+ABX*I_ERI_Gxy2z_S_Px_S_aa;
  Double I_ERI_Gx3z_Px_Px_S_aa = I_ERI_H2x3z_S_Px_S_aa+ABX*I_ERI_Gx3z_S_Px_S_aa;
  Double I_ERI_G4y_Px_Px_S_aa = I_ERI_Hx4y_S_Px_S_aa+ABX*I_ERI_G4y_S_Px_S_aa;
  Double I_ERI_G3yz_Px_Px_S_aa = I_ERI_Hx3yz_S_Px_S_aa+ABX*I_ERI_G3yz_S_Px_S_aa;
  Double I_ERI_G2y2z_Px_Px_S_aa = I_ERI_Hx2y2z_S_Px_S_aa+ABX*I_ERI_G2y2z_S_Px_S_aa;
  Double I_ERI_Gy3z_Px_Px_S_aa = I_ERI_Hxy3z_S_Px_S_aa+ABX*I_ERI_Gy3z_S_Px_S_aa;
  Double I_ERI_G4z_Px_Px_S_aa = I_ERI_Hx4z_S_Px_S_aa+ABX*I_ERI_G4z_S_Px_S_aa;
  Double I_ERI_G4x_Py_Px_S_aa = I_ERI_H4xy_S_Px_S_aa+ABY*I_ERI_G4x_S_Px_S_aa;
  Double I_ERI_G3xy_Py_Px_S_aa = I_ERI_H3x2y_S_Px_S_aa+ABY*I_ERI_G3xy_S_Px_S_aa;
  Double I_ERI_G3xz_Py_Px_S_aa = I_ERI_H3xyz_S_Px_S_aa+ABY*I_ERI_G3xz_S_Px_S_aa;
  Double I_ERI_G2x2y_Py_Px_S_aa = I_ERI_H2x3y_S_Px_S_aa+ABY*I_ERI_G2x2y_S_Px_S_aa;
  Double I_ERI_G2xyz_Py_Px_S_aa = I_ERI_H2x2yz_S_Px_S_aa+ABY*I_ERI_G2xyz_S_Px_S_aa;
  Double I_ERI_G2x2z_Py_Px_S_aa = I_ERI_H2xy2z_S_Px_S_aa+ABY*I_ERI_G2x2z_S_Px_S_aa;
  Double I_ERI_Gx3y_Py_Px_S_aa = I_ERI_Hx4y_S_Px_S_aa+ABY*I_ERI_Gx3y_S_Px_S_aa;
  Double I_ERI_Gx2yz_Py_Px_S_aa = I_ERI_Hx3yz_S_Px_S_aa+ABY*I_ERI_Gx2yz_S_Px_S_aa;
  Double I_ERI_Gxy2z_Py_Px_S_aa = I_ERI_Hx2y2z_S_Px_S_aa+ABY*I_ERI_Gxy2z_S_Px_S_aa;
  Double I_ERI_Gx3z_Py_Px_S_aa = I_ERI_Hxy3z_S_Px_S_aa+ABY*I_ERI_Gx3z_S_Px_S_aa;
  Double I_ERI_G4y_Py_Px_S_aa = I_ERI_H5y_S_Px_S_aa+ABY*I_ERI_G4y_S_Px_S_aa;
  Double I_ERI_G3yz_Py_Px_S_aa = I_ERI_H4yz_S_Px_S_aa+ABY*I_ERI_G3yz_S_Px_S_aa;
  Double I_ERI_G2y2z_Py_Px_S_aa = I_ERI_H3y2z_S_Px_S_aa+ABY*I_ERI_G2y2z_S_Px_S_aa;
  Double I_ERI_Gy3z_Py_Px_S_aa = I_ERI_H2y3z_S_Px_S_aa+ABY*I_ERI_Gy3z_S_Px_S_aa;
  Double I_ERI_G4z_Py_Px_S_aa = I_ERI_Hy4z_S_Px_S_aa+ABY*I_ERI_G4z_S_Px_S_aa;
  Double I_ERI_G4x_Pz_Px_S_aa = I_ERI_H4xz_S_Px_S_aa+ABZ*I_ERI_G4x_S_Px_S_aa;
  Double I_ERI_G3xy_Pz_Px_S_aa = I_ERI_H3xyz_S_Px_S_aa+ABZ*I_ERI_G3xy_S_Px_S_aa;
  Double I_ERI_G3xz_Pz_Px_S_aa = I_ERI_H3x2z_S_Px_S_aa+ABZ*I_ERI_G3xz_S_Px_S_aa;
  Double I_ERI_G2x2y_Pz_Px_S_aa = I_ERI_H2x2yz_S_Px_S_aa+ABZ*I_ERI_G2x2y_S_Px_S_aa;
  Double I_ERI_G2xyz_Pz_Px_S_aa = I_ERI_H2xy2z_S_Px_S_aa+ABZ*I_ERI_G2xyz_S_Px_S_aa;
  Double I_ERI_G2x2z_Pz_Px_S_aa = I_ERI_H2x3z_S_Px_S_aa+ABZ*I_ERI_G2x2z_S_Px_S_aa;
  Double I_ERI_Gx3y_Pz_Px_S_aa = I_ERI_Hx3yz_S_Px_S_aa+ABZ*I_ERI_Gx3y_S_Px_S_aa;
  Double I_ERI_Gx2yz_Pz_Px_S_aa = I_ERI_Hx2y2z_S_Px_S_aa+ABZ*I_ERI_Gx2yz_S_Px_S_aa;
  Double I_ERI_Gxy2z_Pz_Px_S_aa = I_ERI_Hxy3z_S_Px_S_aa+ABZ*I_ERI_Gxy2z_S_Px_S_aa;
  Double I_ERI_Gx3z_Pz_Px_S_aa = I_ERI_Hx4z_S_Px_S_aa+ABZ*I_ERI_Gx3z_S_Px_S_aa;
  Double I_ERI_G4y_Pz_Px_S_aa = I_ERI_H4yz_S_Px_S_aa+ABZ*I_ERI_G4y_S_Px_S_aa;
  Double I_ERI_G3yz_Pz_Px_S_aa = I_ERI_H3y2z_S_Px_S_aa+ABZ*I_ERI_G3yz_S_Px_S_aa;
  Double I_ERI_G2y2z_Pz_Px_S_aa = I_ERI_H2y3z_S_Px_S_aa+ABZ*I_ERI_G2y2z_S_Px_S_aa;
  Double I_ERI_Gy3z_Pz_Px_S_aa = I_ERI_Hy4z_S_Px_S_aa+ABZ*I_ERI_Gy3z_S_Px_S_aa;
  Double I_ERI_G4z_Pz_Px_S_aa = I_ERI_H5z_S_Px_S_aa+ABZ*I_ERI_G4z_S_Px_S_aa;
  Double I_ERI_G4x_Px_Py_S_aa = I_ERI_H5x_S_Py_S_aa+ABX*I_ERI_G4x_S_Py_S_aa;
  Double I_ERI_G3xy_Px_Py_S_aa = I_ERI_H4xy_S_Py_S_aa+ABX*I_ERI_G3xy_S_Py_S_aa;
  Double I_ERI_G3xz_Px_Py_S_aa = I_ERI_H4xz_S_Py_S_aa+ABX*I_ERI_G3xz_S_Py_S_aa;
  Double I_ERI_G2x2y_Px_Py_S_aa = I_ERI_H3x2y_S_Py_S_aa+ABX*I_ERI_G2x2y_S_Py_S_aa;
  Double I_ERI_G2xyz_Px_Py_S_aa = I_ERI_H3xyz_S_Py_S_aa+ABX*I_ERI_G2xyz_S_Py_S_aa;
  Double I_ERI_G2x2z_Px_Py_S_aa = I_ERI_H3x2z_S_Py_S_aa+ABX*I_ERI_G2x2z_S_Py_S_aa;
  Double I_ERI_Gx3y_Px_Py_S_aa = I_ERI_H2x3y_S_Py_S_aa+ABX*I_ERI_Gx3y_S_Py_S_aa;
  Double I_ERI_Gx2yz_Px_Py_S_aa = I_ERI_H2x2yz_S_Py_S_aa+ABX*I_ERI_Gx2yz_S_Py_S_aa;
  Double I_ERI_Gxy2z_Px_Py_S_aa = I_ERI_H2xy2z_S_Py_S_aa+ABX*I_ERI_Gxy2z_S_Py_S_aa;
  Double I_ERI_Gx3z_Px_Py_S_aa = I_ERI_H2x3z_S_Py_S_aa+ABX*I_ERI_Gx3z_S_Py_S_aa;
  Double I_ERI_G4y_Px_Py_S_aa = I_ERI_Hx4y_S_Py_S_aa+ABX*I_ERI_G4y_S_Py_S_aa;
  Double I_ERI_G3yz_Px_Py_S_aa = I_ERI_Hx3yz_S_Py_S_aa+ABX*I_ERI_G3yz_S_Py_S_aa;
  Double I_ERI_G2y2z_Px_Py_S_aa = I_ERI_Hx2y2z_S_Py_S_aa+ABX*I_ERI_G2y2z_S_Py_S_aa;
  Double I_ERI_Gy3z_Px_Py_S_aa = I_ERI_Hxy3z_S_Py_S_aa+ABX*I_ERI_Gy3z_S_Py_S_aa;
  Double I_ERI_G4z_Px_Py_S_aa = I_ERI_Hx4z_S_Py_S_aa+ABX*I_ERI_G4z_S_Py_S_aa;
  Double I_ERI_G4x_Py_Py_S_aa = I_ERI_H4xy_S_Py_S_aa+ABY*I_ERI_G4x_S_Py_S_aa;
  Double I_ERI_G3xy_Py_Py_S_aa = I_ERI_H3x2y_S_Py_S_aa+ABY*I_ERI_G3xy_S_Py_S_aa;
  Double I_ERI_G3xz_Py_Py_S_aa = I_ERI_H3xyz_S_Py_S_aa+ABY*I_ERI_G3xz_S_Py_S_aa;
  Double I_ERI_G2x2y_Py_Py_S_aa = I_ERI_H2x3y_S_Py_S_aa+ABY*I_ERI_G2x2y_S_Py_S_aa;
  Double I_ERI_G2xyz_Py_Py_S_aa = I_ERI_H2x2yz_S_Py_S_aa+ABY*I_ERI_G2xyz_S_Py_S_aa;
  Double I_ERI_G2x2z_Py_Py_S_aa = I_ERI_H2xy2z_S_Py_S_aa+ABY*I_ERI_G2x2z_S_Py_S_aa;
  Double I_ERI_Gx3y_Py_Py_S_aa = I_ERI_Hx4y_S_Py_S_aa+ABY*I_ERI_Gx3y_S_Py_S_aa;
  Double I_ERI_Gx2yz_Py_Py_S_aa = I_ERI_Hx3yz_S_Py_S_aa+ABY*I_ERI_Gx2yz_S_Py_S_aa;
  Double I_ERI_Gxy2z_Py_Py_S_aa = I_ERI_Hx2y2z_S_Py_S_aa+ABY*I_ERI_Gxy2z_S_Py_S_aa;
  Double I_ERI_Gx3z_Py_Py_S_aa = I_ERI_Hxy3z_S_Py_S_aa+ABY*I_ERI_Gx3z_S_Py_S_aa;
  Double I_ERI_G4y_Py_Py_S_aa = I_ERI_H5y_S_Py_S_aa+ABY*I_ERI_G4y_S_Py_S_aa;
  Double I_ERI_G3yz_Py_Py_S_aa = I_ERI_H4yz_S_Py_S_aa+ABY*I_ERI_G3yz_S_Py_S_aa;
  Double I_ERI_G2y2z_Py_Py_S_aa = I_ERI_H3y2z_S_Py_S_aa+ABY*I_ERI_G2y2z_S_Py_S_aa;
  Double I_ERI_Gy3z_Py_Py_S_aa = I_ERI_H2y3z_S_Py_S_aa+ABY*I_ERI_Gy3z_S_Py_S_aa;
  Double I_ERI_G4z_Py_Py_S_aa = I_ERI_Hy4z_S_Py_S_aa+ABY*I_ERI_G4z_S_Py_S_aa;
  Double I_ERI_G4x_Pz_Py_S_aa = I_ERI_H4xz_S_Py_S_aa+ABZ*I_ERI_G4x_S_Py_S_aa;
  Double I_ERI_G3xy_Pz_Py_S_aa = I_ERI_H3xyz_S_Py_S_aa+ABZ*I_ERI_G3xy_S_Py_S_aa;
  Double I_ERI_G3xz_Pz_Py_S_aa = I_ERI_H3x2z_S_Py_S_aa+ABZ*I_ERI_G3xz_S_Py_S_aa;
  Double I_ERI_G2x2y_Pz_Py_S_aa = I_ERI_H2x2yz_S_Py_S_aa+ABZ*I_ERI_G2x2y_S_Py_S_aa;
  Double I_ERI_G2xyz_Pz_Py_S_aa = I_ERI_H2xy2z_S_Py_S_aa+ABZ*I_ERI_G2xyz_S_Py_S_aa;
  Double I_ERI_G2x2z_Pz_Py_S_aa = I_ERI_H2x3z_S_Py_S_aa+ABZ*I_ERI_G2x2z_S_Py_S_aa;
  Double I_ERI_Gx3y_Pz_Py_S_aa = I_ERI_Hx3yz_S_Py_S_aa+ABZ*I_ERI_Gx3y_S_Py_S_aa;
  Double I_ERI_Gx2yz_Pz_Py_S_aa = I_ERI_Hx2y2z_S_Py_S_aa+ABZ*I_ERI_Gx2yz_S_Py_S_aa;
  Double I_ERI_Gxy2z_Pz_Py_S_aa = I_ERI_Hxy3z_S_Py_S_aa+ABZ*I_ERI_Gxy2z_S_Py_S_aa;
  Double I_ERI_Gx3z_Pz_Py_S_aa = I_ERI_Hx4z_S_Py_S_aa+ABZ*I_ERI_Gx3z_S_Py_S_aa;
  Double I_ERI_G4y_Pz_Py_S_aa = I_ERI_H4yz_S_Py_S_aa+ABZ*I_ERI_G4y_S_Py_S_aa;
  Double I_ERI_G3yz_Pz_Py_S_aa = I_ERI_H3y2z_S_Py_S_aa+ABZ*I_ERI_G3yz_S_Py_S_aa;
  Double I_ERI_G2y2z_Pz_Py_S_aa = I_ERI_H2y3z_S_Py_S_aa+ABZ*I_ERI_G2y2z_S_Py_S_aa;
  Double I_ERI_Gy3z_Pz_Py_S_aa = I_ERI_Hy4z_S_Py_S_aa+ABZ*I_ERI_Gy3z_S_Py_S_aa;
  Double I_ERI_G4z_Pz_Py_S_aa = I_ERI_H5z_S_Py_S_aa+ABZ*I_ERI_G4z_S_Py_S_aa;
  Double I_ERI_G4x_Px_Pz_S_aa = I_ERI_H5x_S_Pz_S_aa+ABX*I_ERI_G4x_S_Pz_S_aa;
  Double I_ERI_G3xy_Px_Pz_S_aa = I_ERI_H4xy_S_Pz_S_aa+ABX*I_ERI_G3xy_S_Pz_S_aa;
  Double I_ERI_G3xz_Px_Pz_S_aa = I_ERI_H4xz_S_Pz_S_aa+ABX*I_ERI_G3xz_S_Pz_S_aa;
  Double I_ERI_G2x2y_Px_Pz_S_aa = I_ERI_H3x2y_S_Pz_S_aa+ABX*I_ERI_G2x2y_S_Pz_S_aa;
  Double I_ERI_G2xyz_Px_Pz_S_aa = I_ERI_H3xyz_S_Pz_S_aa+ABX*I_ERI_G2xyz_S_Pz_S_aa;
  Double I_ERI_G2x2z_Px_Pz_S_aa = I_ERI_H3x2z_S_Pz_S_aa+ABX*I_ERI_G2x2z_S_Pz_S_aa;
  Double I_ERI_Gx3y_Px_Pz_S_aa = I_ERI_H2x3y_S_Pz_S_aa+ABX*I_ERI_Gx3y_S_Pz_S_aa;
  Double I_ERI_Gx2yz_Px_Pz_S_aa = I_ERI_H2x2yz_S_Pz_S_aa+ABX*I_ERI_Gx2yz_S_Pz_S_aa;
  Double I_ERI_Gxy2z_Px_Pz_S_aa = I_ERI_H2xy2z_S_Pz_S_aa+ABX*I_ERI_Gxy2z_S_Pz_S_aa;
  Double I_ERI_Gx3z_Px_Pz_S_aa = I_ERI_H2x3z_S_Pz_S_aa+ABX*I_ERI_Gx3z_S_Pz_S_aa;
  Double I_ERI_G4y_Px_Pz_S_aa = I_ERI_Hx4y_S_Pz_S_aa+ABX*I_ERI_G4y_S_Pz_S_aa;
  Double I_ERI_G3yz_Px_Pz_S_aa = I_ERI_Hx3yz_S_Pz_S_aa+ABX*I_ERI_G3yz_S_Pz_S_aa;
  Double I_ERI_G2y2z_Px_Pz_S_aa = I_ERI_Hx2y2z_S_Pz_S_aa+ABX*I_ERI_G2y2z_S_Pz_S_aa;
  Double I_ERI_Gy3z_Px_Pz_S_aa = I_ERI_Hxy3z_S_Pz_S_aa+ABX*I_ERI_Gy3z_S_Pz_S_aa;
  Double I_ERI_G4z_Px_Pz_S_aa = I_ERI_Hx4z_S_Pz_S_aa+ABX*I_ERI_G4z_S_Pz_S_aa;
  Double I_ERI_G4x_Py_Pz_S_aa = I_ERI_H4xy_S_Pz_S_aa+ABY*I_ERI_G4x_S_Pz_S_aa;
  Double I_ERI_G3xy_Py_Pz_S_aa = I_ERI_H3x2y_S_Pz_S_aa+ABY*I_ERI_G3xy_S_Pz_S_aa;
  Double I_ERI_G3xz_Py_Pz_S_aa = I_ERI_H3xyz_S_Pz_S_aa+ABY*I_ERI_G3xz_S_Pz_S_aa;
  Double I_ERI_G2x2y_Py_Pz_S_aa = I_ERI_H2x3y_S_Pz_S_aa+ABY*I_ERI_G2x2y_S_Pz_S_aa;
  Double I_ERI_G2xyz_Py_Pz_S_aa = I_ERI_H2x2yz_S_Pz_S_aa+ABY*I_ERI_G2xyz_S_Pz_S_aa;
  Double I_ERI_G2x2z_Py_Pz_S_aa = I_ERI_H2xy2z_S_Pz_S_aa+ABY*I_ERI_G2x2z_S_Pz_S_aa;
  Double I_ERI_Gx3y_Py_Pz_S_aa = I_ERI_Hx4y_S_Pz_S_aa+ABY*I_ERI_Gx3y_S_Pz_S_aa;
  Double I_ERI_Gx2yz_Py_Pz_S_aa = I_ERI_Hx3yz_S_Pz_S_aa+ABY*I_ERI_Gx2yz_S_Pz_S_aa;
  Double I_ERI_Gxy2z_Py_Pz_S_aa = I_ERI_Hx2y2z_S_Pz_S_aa+ABY*I_ERI_Gxy2z_S_Pz_S_aa;
  Double I_ERI_Gx3z_Py_Pz_S_aa = I_ERI_Hxy3z_S_Pz_S_aa+ABY*I_ERI_Gx3z_S_Pz_S_aa;
  Double I_ERI_G4y_Py_Pz_S_aa = I_ERI_H5y_S_Pz_S_aa+ABY*I_ERI_G4y_S_Pz_S_aa;
  Double I_ERI_G3yz_Py_Pz_S_aa = I_ERI_H4yz_S_Pz_S_aa+ABY*I_ERI_G3yz_S_Pz_S_aa;
  Double I_ERI_G2y2z_Py_Pz_S_aa = I_ERI_H3y2z_S_Pz_S_aa+ABY*I_ERI_G2y2z_S_Pz_S_aa;
  Double I_ERI_Gy3z_Py_Pz_S_aa = I_ERI_H2y3z_S_Pz_S_aa+ABY*I_ERI_Gy3z_S_Pz_S_aa;
  Double I_ERI_G4z_Py_Pz_S_aa = I_ERI_Hy4z_S_Pz_S_aa+ABY*I_ERI_G4z_S_Pz_S_aa;
  Double I_ERI_G4x_Pz_Pz_S_aa = I_ERI_H4xz_S_Pz_S_aa+ABZ*I_ERI_G4x_S_Pz_S_aa;
  Double I_ERI_G3xy_Pz_Pz_S_aa = I_ERI_H3xyz_S_Pz_S_aa+ABZ*I_ERI_G3xy_S_Pz_S_aa;
  Double I_ERI_G3xz_Pz_Pz_S_aa = I_ERI_H3x2z_S_Pz_S_aa+ABZ*I_ERI_G3xz_S_Pz_S_aa;
  Double I_ERI_G2x2y_Pz_Pz_S_aa = I_ERI_H2x2yz_S_Pz_S_aa+ABZ*I_ERI_G2x2y_S_Pz_S_aa;
  Double I_ERI_G2xyz_Pz_Pz_S_aa = I_ERI_H2xy2z_S_Pz_S_aa+ABZ*I_ERI_G2xyz_S_Pz_S_aa;
  Double I_ERI_G2x2z_Pz_Pz_S_aa = I_ERI_H2x3z_S_Pz_S_aa+ABZ*I_ERI_G2x2z_S_Pz_S_aa;
  Double I_ERI_Gx3y_Pz_Pz_S_aa = I_ERI_Hx3yz_S_Pz_S_aa+ABZ*I_ERI_Gx3y_S_Pz_S_aa;
  Double I_ERI_Gx2yz_Pz_Pz_S_aa = I_ERI_Hx2y2z_S_Pz_S_aa+ABZ*I_ERI_Gx2yz_S_Pz_S_aa;
  Double I_ERI_Gxy2z_Pz_Pz_S_aa = I_ERI_Hxy3z_S_Pz_S_aa+ABZ*I_ERI_Gxy2z_S_Pz_S_aa;
  Double I_ERI_Gx3z_Pz_Pz_S_aa = I_ERI_Hx4z_S_Pz_S_aa+ABZ*I_ERI_Gx3z_S_Pz_S_aa;
  Double I_ERI_G4y_Pz_Pz_S_aa = I_ERI_H4yz_S_Pz_S_aa+ABZ*I_ERI_G4y_S_Pz_S_aa;
  Double I_ERI_G3yz_Pz_Pz_S_aa = I_ERI_H3y2z_S_Pz_S_aa+ABZ*I_ERI_G3yz_S_Pz_S_aa;
  Double I_ERI_G2y2z_Pz_Pz_S_aa = I_ERI_H2y3z_S_Pz_S_aa+ABZ*I_ERI_G2y2z_S_Pz_S_aa;
  Double I_ERI_Gy3z_Pz_Pz_S_aa = I_ERI_Hy4z_S_Pz_S_aa+ABZ*I_ERI_Gy3z_S_Pz_S_aa;
  Double I_ERI_G4z_Pz_Pz_S_aa = I_ERI_H5z_S_Pz_S_aa+ABZ*I_ERI_G4z_S_Pz_S_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ab
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_ab = I_ERI_G4x_S_Px_S_ab+ABX*I_ERI_F3x_S_Px_S_ab;
  Double I_ERI_F2xy_Px_Px_S_ab = I_ERI_G3xy_S_Px_S_ab+ABX*I_ERI_F2xy_S_Px_S_ab;
  Double I_ERI_F2xz_Px_Px_S_ab = I_ERI_G3xz_S_Px_S_ab+ABX*I_ERI_F2xz_S_Px_S_ab;
  Double I_ERI_Fx2y_Px_Px_S_ab = I_ERI_G2x2y_S_Px_S_ab+ABX*I_ERI_Fx2y_S_Px_S_ab;
  Double I_ERI_Fxyz_Px_Px_S_ab = I_ERI_G2xyz_S_Px_S_ab+ABX*I_ERI_Fxyz_S_Px_S_ab;
  Double I_ERI_Fx2z_Px_Px_S_ab = I_ERI_G2x2z_S_Px_S_ab+ABX*I_ERI_Fx2z_S_Px_S_ab;
  Double I_ERI_F3y_Px_Px_S_ab = I_ERI_Gx3y_S_Px_S_ab+ABX*I_ERI_F3y_S_Px_S_ab;
  Double I_ERI_F2yz_Px_Px_S_ab = I_ERI_Gx2yz_S_Px_S_ab+ABX*I_ERI_F2yz_S_Px_S_ab;
  Double I_ERI_Fy2z_Px_Px_S_ab = I_ERI_Gxy2z_S_Px_S_ab+ABX*I_ERI_Fy2z_S_Px_S_ab;
  Double I_ERI_F3z_Px_Px_S_ab = I_ERI_Gx3z_S_Px_S_ab+ABX*I_ERI_F3z_S_Px_S_ab;
  Double I_ERI_F3x_Py_Px_S_ab = I_ERI_G3xy_S_Px_S_ab+ABY*I_ERI_F3x_S_Px_S_ab;
  Double I_ERI_F2xy_Py_Px_S_ab = I_ERI_G2x2y_S_Px_S_ab+ABY*I_ERI_F2xy_S_Px_S_ab;
  Double I_ERI_F2xz_Py_Px_S_ab = I_ERI_G2xyz_S_Px_S_ab+ABY*I_ERI_F2xz_S_Px_S_ab;
  Double I_ERI_Fx2y_Py_Px_S_ab = I_ERI_Gx3y_S_Px_S_ab+ABY*I_ERI_Fx2y_S_Px_S_ab;
  Double I_ERI_Fxyz_Py_Px_S_ab = I_ERI_Gx2yz_S_Px_S_ab+ABY*I_ERI_Fxyz_S_Px_S_ab;
  Double I_ERI_Fx2z_Py_Px_S_ab = I_ERI_Gxy2z_S_Px_S_ab+ABY*I_ERI_Fx2z_S_Px_S_ab;
  Double I_ERI_F3y_Py_Px_S_ab = I_ERI_G4y_S_Px_S_ab+ABY*I_ERI_F3y_S_Px_S_ab;
  Double I_ERI_F2yz_Py_Px_S_ab = I_ERI_G3yz_S_Px_S_ab+ABY*I_ERI_F2yz_S_Px_S_ab;
  Double I_ERI_Fy2z_Py_Px_S_ab = I_ERI_G2y2z_S_Px_S_ab+ABY*I_ERI_Fy2z_S_Px_S_ab;
  Double I_ERI_F3z_Py_Px_S_ab = I_ERI_Gy3z_S_Px_S_ab+ABY*I_ERI_F3z_S_Px_S_ab;
  Double I_ERI_F3x_Pz_Px_S_ab = I_ERI_G3xz_S_Px_S_ab+ABZ*I_ERI_F3x_S_Px_S_ab;
  Double I_ERI_F2xy_Pz_Px_S_ab = I_ERI_G2xyz_S_Px_S_ab+ABZ*I_ERI_F2xy_S_Px_S_ab;
  Double I_ERI_F2xz_Pz_Px_S_ab = I_ERI_G2x2z_S_Px_S_ab+ABZ*I_ERI_F2xz_S_Px_S_ab;
  Double I_ERI_Fx2y_Pz_Px_S_ab = I_ERI_Gx2yz_S_Px_S_ab+ABZ*I_ERI_Fx2y_S_Px_S_ab;
  Double I_ERI_Fxyz_Pz_Px_S_ab = I_ERI_Gxy2z_S_Px_S_ab+ABZ*I_ERI_Fxyz_S_Px_S_ab;
  Double I_ERI_Fx2z_Pz_Px_S_ab = I_ERI_Gx3z_S_Px_S_ab+ABZ*I_ERI_Fx2z_S_Px_S_ab;
  Double I_ERI_F3y_Pz_Px_S_ab = I_ERI_G3yz_S_Px_S_ab+ABZ*I_ERI_F3y_S_Px_S_ab;
  Double I_ERI_F2yz_Pz_Px_S_ab = I_ERI_G2y2z_S_Px_S_ab+ABZ*I_ERI_F2yz_S_Px_S_ab;
  Double I_ERI_Fy2z_Pz_Px_S_ab = I_ERI_Gy3z_S_Px_S_ab+ABZ*I_ERI_Fy2z_S_Px_S_ab;
  Double I_ERI_F3z_Pz_Px_S_ab = I_ERI_G4z_S_Px_S_ab+ABZ*I_ERI_F3z_S_Px_S_ab;
  Double I_ERI_F3x_Px_Py_S_ab = I_ERI_G4x_S_Py_S_ab+ABX*I_ERI_F3x_S_Py_S_ab;
  Double I_ERI_F2xy_Px_Py_S_ab = I_ERI_G3xy_S_Py_S_ab+ABX*I_ERI_F2xy_S_Py_S_ab;
  Double I_ERI_F2xz_Px_Py_S_ab = I_ERI_G3xz_S_Py_S_ab+ABX*I_ERI_F2xz_S_Py_S_ab;
  Double I_ERI_Fx2y_Px_Py_S_ab = I_ERI_G2x2y_S_Py_S_ab+ABX*I_ERI_Fx2y_S_Py_S_ab;
  Double I_ERI_Fxyz_Px_Py_S_ab = I_ERI_G2xyz_S_Py_S_ab+ABX*I_ERI_Fxyz_S_Py_S_ab;
  Double I_ERI_Fx2z_Px_Py_S_ab = I_ERI_G2x2z_S_Py_S_ab+ABX*I_ERI_Fx2z_S_Py_S_ab;
  Double I_ERI_F3y_Px_Py_S_ab = I_ERI_Gx3y_S_Py_S_ab+ABX*I_ERI_F3y_S_Py_S_ab;
  Double I_ERI_F2yz_Px_Py_S_ab = I_ERI_Gx2yz_S_Py_S_ab+ABX*I_ERI_F2yz_S_Py_S_ab;
  Double I_ERI_Fy2z_Px_Py_S_ab = I_ERI_Gxy2z_S_Py_S_ab+ABX*I_ERI_Fy2z_S_Py_S_ab;
  Double I_ERI_F3z_Px_Py_S_ab = I_ERI_Gx3z_S_Py_S_ab+ABX*I_ERI_F3z_S_Py_S_ab;
  Double I_ERI_F3x_Py_Py_S_ab = I_ERI_G3xy_S_Py_S_ab+ABY*I_ERI_F3x_S_Py_S_ab;
  Double I_ERI_F2xy_Py_Py_S_ab = I_ERI_G2x2y_S_Py_S_ab+ABY*I_ERI_F2xy_S_Py_S_ab;
  Double I_ERI_F2xz_Py_Py_S_ab = I_ERI_G2xyz_S_Py_S_ab+ABY*I_ERI_F2xz_S_Py_S_ab;
  Double I_ERI_Fx2y_Py_Py_S_ab = I_ERI_Gx3y_S_Py_S_ab+ABY*I_ERI_Fx2y_S_Py_S_ab;
  Double I_ERI_Fxyz_Py_Py_S_ab = I_ERI_Gx2yz_S_Py_S_ab+ABY*I_ERI_Fxyz_S_Py_S_ab;
  Double I_ERI_Fx2z_Py_Py_S_ab = I_ERI_Gxy2z_S_Py_S_ab+ABY*I_ERI_Fx2z_S_Py_S_ab;
  Double I_ERI_F3y_Py_Py_S_ab = I_ERI_G4y_S_Py_S_ab+ABY*I_ERI_F3y_S_Py_S_ab;
  Double I_ERI_F2yz_Py_Py_S_ab = I_ERI_G3yz_S_Py_S_ab+ABY*I_ERI_F2yz_S_Py_S_ab;
  Double I_ERI_Fy2z_Py_Py_S_ab = I_ERI_G2y2z_S_Py_S_ab+ABY*I_ERI_Fy2z_S_Py_S_ab;
  Double I_ERI_F3z_Py_Py_S_ab = I_ERI_Gy3z_S_Py_S_ab+ABY*I_ERI_F3z_S_Py_S_ab;
  Double I_ERI_F3x_Pz_Py_S_ab = I_ERI_G3xz_S_Py_S_ab+ABZ*I_ERI_F3x_S_Py_S_ab;
  Double I_ERI_F2xy_Pz_Py_S_ab = I_ERI_G2xyz_S_Py_S_ab+ABZ*I_ERI_F2xy_S_Py_S_ab;
  Double I_ERI_F2xz_Pz_Py_S_ab = I_ERI_G2x2z_S_Py_S_ab+ABZ*I_ERI_F2xz_S_Py_S_ab;
  Double I_ERI_Fx2y_Pz_Py_S_ab = I_ERI_Gx2yz_S_Py_S_ab+ABZ*I_ERI_Fx2y_S_Py_S_ab;
  Double I_ERI_Fxyz_Pz_Py_S_ab = I_ERI_Gxy2z_S_Py_S_ab+ABZ*I_ERI_Fxyz_S_Py_S_ab;
  Double I_ERI_Fx2z_Pz_Py_S_ab = I_ERI_Gx3z_S_Py_S_ab+ABZ*I_ERI_Fx2z_S_Py_S_ab;
  Double I_ERI_F3y_Pz_Py_S_ab = I_ERI_G3yz_S_Py_S_ab+ABZ*I_ERI_F3y_S_Py_S_ab;
  Double I_ERI_F2yz_Pz_Py_S_ab = I_ERI_G2y2z_S_Py_S_ab+ABZ*I_ERI_F2yz_S_Py_S_ab;
  Double I_ERI_Fy2z_Pz_Py_S_ab = I_ERI_Gy3z_S_Py_S_ab+ABZ*I_ERI_Fy2z_S_Py_S_ab;
  Double I_ERI_F3z_Pz_Py_S_ab = I_ERI_G4z_S_Py_S_ab+ABZ*I_ERI_F3z_S_Py_S_ab;
  Double I_ERI_F3x_Px_Pz_S_ab = I_ERI_G4x_S_Pz_S_ab+ABX*I_ERI_F3x_S_Pz_S_ab;
  Double I_ERI_F2xy_Px_Pz_S_ab = I_ERI_G3xy_S_Pz_S_ab+ABX*I_ERI_F2xy_S_Pz_S_ab;
  Double I_ERI_F2xz_Px_Pz_S_ab = I_ERI_G3xz_S_Pz_S_ab+ABX*I_ERI_F2xz_S_Pz_S_ab;
  Double I_ERI_Fx2y_Px_Pz_S_ab = I_ERI_G2x2y_S_Pz_S_ab+ABX*I_ERI_Fx2y_S_Pz_S_ab;
  Double I_ERI_Fxyz_Px_Pz_S_ab = I_ERI_G2xyz_S_Pz_S_ab+ABX*I_ERI_Fxyz_S_Pz_S_ab;
  Double I_ERI_Fx2z_Px_Pz_S_ab = I_ERI_G2x2z_S_Pz_S_ab+ABX*I_ERI_Fx2z_S_Pz_S_ab;
  Double I_ERI_F3y_Px_Pz_S_ab = I_ERI_Gx3y_S_Pz_S_ab+ABX*I_ERI_F3y_S_Pz_S_ab;
  Double I_ERI_F2yz_Px_Pz_S_ab = I_ERI_Gx2yz_S_Pz_S_ab+ABX*I_ERI_F2yz_S_Pz_S_ab;
  Double I_ERI_Fy2z_Px_Pz_S_ab = I_ERI_Gxy2z_S_Pz_S_ab+ABX*I_ERI_Fy2z_S_Pz_S_ab;
  Double I_ERI_F3z_Px_Pz_S_ab = I_ERI_Gx3z_S_Pz_S_ab+ABX*I_ERI_F3z_S_Pz_S_ab;
  Double I_ERI_F3x_Py_Pz_S_ab = I_ERI_G3xy_S_Pz_S_ab+ABY*I_ERI_F3x_S_Pz_S_ab;
  Double I_ERI_F2xy_Py_Pz_S_ab = I_ERI_G2x2y_S_Pz_S_ab+ABY*I_ERI_F2xy_S_Pz_S_ab;
  Double I_ERI_F2xz_Py_Pz_S_ab = I_ERI_G2xyz_S_Pz_S_ab+ABY*I_ERI_F2xz_S_Pz_S_ab;
  Double I_ERI_Fx2y_Py_Pz_S_ab = I_ERI_Gx3y_S_Pz_S_ab+ABY*I_ERI_Fx2y_S_Pz_S_ab;
  Double I_ERI_Fxyz_Py_Pz_S_ab = I_ERI_Gx2yz_S_Pz_S_ab+ABY*I_ERI_Fxyz_S_Pz_S_ab;
  Double I_ERI_Fx2z_Py_Pz_S_ab = I_ERI_Gxy2z_S_Pz_S_ab+ABY*I_ERI_Fx2z_S_Pz_S_ab;
  Double I_ERI_F3y_Py_Pz_S_ab = I_ERI_G4y_S_Pz_S_ab+ABY*I_ERI_F3y_S_Pz_S_ab;
  Double I_ERI_F2yz_Py_Pz_S_ab = I_ERI_G3yz_S_Pz_S_ab+ABY*I_ERI_F2yz_S_Pz_S_ab;
  Double I_ERI_Fy2z_Py_Pz_S_ab = I_ERI_G2y2z_S_Pz_S_ab+ABY*I_ERI_Fy2z_S_Pz_S_ab;
  Double I_ERI_F3z_Py_Pz_S_ab = I_ERI_Gy3z_S_Pz_S_ab+ABY*I_ERI_F3z_S_Pz_S_ab;
  Double I_ERI_F3x_Pz_Pz_S_ab = I_ERI_G3xz_S_Pz_S_ab+ABZ*I_ERI_F3x_S_Pz_S_ab;
  Double I_ERI_F2xy_Pz_Pz_S_ab = I_ERI_G2xyz_S_Pz_S_ab+ABZ*I_ERI_F2xy_S_Pz_S_ab;
  Double I_ERI_F2xz_Pz_Pz_S_ab = I_ERI_G2x2z_S_Pz_S_ab+ABZ*I_ERI_F2xz_S_Pz_S_ab;
  Double I_ERI_Fx2y_Pz_Pz_S_ab = I_ERI_Gx2yz_S_Pz_S_ab+ABZ*I_ERI_Fx2y_S_Pz_S_ab;
  Double I_ERI_Fxyz_Pz_Pz_S_ab = I_ERI_Gxy2z_S_Pz_S_ab+ABZ*I_ERI_Fxyz_S_Pz_S_ab;
  Double I_ERI_Fx2z_Pz_Pz_S_ab = I_ERI_Gx3z_S_Pz_S_ab+ABZ*I_ERI_Fx2z_S_Pz_S_ab;
  Double I_ERI_F3y_Pz_Pz_S_ab = I_ERI_G3yz_S_Pz_S_ab+ABZ*I_ERI_F3y_S_Pz_S_ab;
  Double I_ERI_F2yz_Pz_Pz_S_ab = I_ERI_G2y2z_S_Pz_S_ab+ABZ*I_ERI_F2yz_S_Pz_S_ab;
  Double I_ERI_Fy2z_Pz_Pz_S_ab = I_ERI_Gy3z_S_Pz_S_ab+ABZ*I_ERI_Fy2z_S_Pz_S_ab;
  Double I_ERI_F3z_Pz_Pz_S_ab = I_ERI_G4z_S_Pz_S_ab+ABZ*I_ERI_F3z_S_Pz_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_P_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 18 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_P_S_ab
   * RHS shell quartet name: SQ_ERI_G_S_P_S_ab
   ************************************************************/
  Double I_ERI_G4x_Px_Px_S_ab = I_ERI_H5x_S_Px_S_ab+ABX*I_ERI_G4x_S_Px_S_ab;
  Double I_ERI_G3xy_Px_Px_S_ab = I_ERI_H4xy_S_Px_S_ab+ABX*I_ERI_G3xy_S_Px_S_ab;
  Double I_ERI_G3xz_Px_Px_S_ab = I_ERI_H4xz_S_Px_S_ab+ABX*I_ERI_G3xz_S_Px_S_ab;
  Double I_ERI_G2x2y_Px_Px_S_ab = I_ERI_H3x2y_S_Px_S_ab+ABX*I_ERI_G2x2y_S_Px_S_ab;
  Double I_ERI_G2xyz_Px_Px_S_ab = I_ERI_H3xyz_S_Px_S_ab+ABX*I_ERI_G2xyz_S_Px_S_ab;
  Double I_ERI_G2x2z_Px_Px_S_ab = I_ERI_H3x2z_S_Px_S_ab+ABX*I_ERI_G2x2z_S_Px_S_ab;
  Double I_ERI_Gx3y_Px_Px_S_ab = I_ERI_H2x3y_S_Px_S_ab+ABX*I_ERI_Gx3y_S_Px_S_ab;
  Double I_ERI_Gx2yz_Px_Px_S_ab = I_ERI_H2x2yz_S_Px_S_ab+ABX*I_ERI_Gx2yz_S_Px_S_ab;
  Double I_ERI_Gxy2z_Px_Px_S_ab = I_ERI_H2xy2z_S_Px_S_ab+ABX*I_ERI_Gxy2z_S_Px_S_ab;
  Double I_ERI_Gx3z_Px_Px_S_ab = I_ERI_H2x3z_S_Px_S_ab+ABX*I_ERI_Gx3z_S_Px_S_ab;
  Double I_ERI_G4y_Px_Px_S_ab = I_ERI_Hx4y_S_Px_S_ab+ABX*I_ERI_G4y_S_Px_S_ab;
  Double I_ERI_G3yz_Px_Px_S_ab = I_ERI_Hx3yz_S_Px_S_ab+ABX*I_ERI_G3yz_S_Px_S_ab;
  Double I_ERI_G2y2z_Px_Px_S_ab = I_ERI_Hx2y2z_S_Px_S_ab+ABX*I_ERI_G2y2z_S_Px_S_ab;
  Double I_ERI_Gy3z_Px_Px_S_ab = I_ERI_Hxy3z_S_Px_S_ab+ABX*I_ERI_Gy3z_S_Px_S_ab;
  Double I_ERI_G4z_Px_Px_S_ab = I_ERI_Hx4z_S_Px_S_ab+ABX*I_ERI_G4z_S_Px_S_ab;
  Double I_ERI_G3xy_Py_Px_S_ab = I_ERI_H3x2y_S_Px_S_ab+ABY*I_ERI_G3xy_S_Px_S_ab;
  Double I_ERI_G3xz_Py_Px_S_ab = I_ERI_H3xyz_S_Px_S_ab+ABY*I_ERI_G3xz_S_Px_S_ab;
  Double I_ERI_G2x2y_Py_Px_S_ab = I_ERI_H2x3y_S_Px_S_ab+ABY*I_ERI_G2x2y_S_Px_S_ab;
  Double I_ERI_G2xyz_Py_Px_S_ab = I_ERI_H2x2yz_S_Px_S_ab+ABY*I_ERI_G2xyz_S_Px_S_ab;
  Double I_ERI_G2x2z_Py_Px_S_ab = I_ERI_H2xy2z_S_Px_S_ab+ABY*I_ERI_G2x2z_S_Px_S_ab;
  Double I_ERI_Gx3y_Py_Px_S_ab = I_ERI_Hx4y_S_Px_S_ab+ABY*I_ERI_Gx3y_S_Px_S_ab;
  Double I_ERI_Gx2yz_Py_Px_S_ab = I_ERI_Hx3yz_S_Px_S_ab+ABY*I_ERI_Gx2yz_S_Px_S_ab;
  Double I_ERI_Gxy2z_Py_Px_S_ab = I_ERI_Hx2y2z_S_Px_S_ab+ABY*I_ERI_Gxy2z_S_Px_S_ab;
  Double I_ERI_Gx3z_Py_Px_S_ab = I_ERI_Hxy3z_S_Px_S_ab+ABY*I_ERI_Gx3z_S_Px_S_ab;
  Double I_ERI_G4y_Py_Px_S_ab = I_ERI_H5y_S_Px_S_ab+ABY*I_ERI_G4y_S_Px_S_ab;
  Double I_ERI_G3yz_Py_Px_S_ab = I_ERI_H4yz_S_Px_S_ab+ABY*I_ERI_G3yz_S_Px_S_ab;
  Double I_ERI_G2y2z_Py_Px_S_ab = I_ERI_H3y2z_S_Px_S_ab+ABY*I_ERI_G2y2z_S_Px_S_ab;
  Double I_ERI_Gy3z_Py_Px_S_ab = I_ERI_H2y3z_S_Px_S_ab+ABY*I_ERI_Gy3z_S_Px_S_ab;
  Double I_ERI_G4z_Py_Px_S_ab = I_ERI_Hy4z_S_Px_S_ab+ABY*I_ERI_G4z_S_Px_S_ab;
  Double I_ERI_G3xz_Pz_Px_S_ab = I_ERI_H3x2z_S_Px_S_ab+ABZ*I_ERI_G3xz_S_Px_S_ab;
  Double I_ERI_G2xyz_Pz_Px_S_ab = I_ERI_H2xy2z_S_Px_S_ab+ABZ*I_ERI_G2xyz_S_Px_S_ab;
  Double I_ERI_G2x2z_Pz_Px_S_ab = I_ERI_H2x3z_S_Px_S_ab+ABZ*I_ERI_G2x2z_S_Px_S_ab;
  Double I_ERI_Gx2yz_Pz_Px_S_ab = I_ERI_Hx2y2z_S_Px_S_ab+ABZ*I_ERI_Gx2yz_S_Px_S_ab;
  Double I_ERI_Gxy2z_Pz_Px_S_ab = I_ERI_Hxy3z_S_Px_S_ab+ABZ*I_ERI_Gxy2z_S_Px_S_ab;
  Double I_ERI_Gx3z_Pz_Px_S_ab = I_ERI_Hx4z_S_Px_S_ab+ABZ*I_ERI_Gx3z_S_Px_S_ab;
  Double I_ERI_G3yz_Pz_Px_S_ab = I_ERI_H3y2z_S_Px_S_ab+ABZ*I_ERI_G3yz_S_Px_S_ab;
  Double I_ERI_G2y2z_Pz_Px_S_ab = I_ERI_H2y3z_S_Px_S_ab+ABZ*I_ERI_G2y2z_S_Px_S_ab;
  Double I_ERI_Gy3z_Pz_Px_S_ab = I_ERI_Hy4z_S_Px_S_ab+ABZ*I_ERI_Gy3z_S_Px_S_ab;
  Double I_ERI_G4z_Pz_Px_S_ab = I_ERI_H5z_S_Px_S_ab+ABZ*I_ERI_G4z_S_Px_S_ab;
  Double I_ERI_G4x_Px_Py_S_ab = I_ERI_H5x_S_Py_S_ab+ABX*I_ERI_G4x_S_Py_S_ab;
  Double I_ERI_G3xy_Px_Py_S_ab = I_ERI_H4xy_S_Py_S_ab+ABX*I_ERI_G3xy_S_Py_S_ab;
  Double I_ERI_G3xz_Px_Py_S_ab = I_ERI_H4xz_S_Py_S_ab+ABX*I_ERI_G3xz_S_Py_S_ab;
  Double I_ERI_G2x2y_Px_Py_S_ab = I_ERI_H3x2y_S_Py_S_ab+ABX*I_ERI_G2x2y_S_Py_S_ab;
  Double I_ERI_G2xyz_Px_Py_S_ab = I_ERI_H3xyz_S_Py_S_ab+ABX*I_ERI_G2xyz_S_Py_S_ab;
  Double I_ERI_G2x2z_Px_Py_S_ab = I_ERI_H3x2z_S_Py_S_ab+ABX*I_ERI_G2x2z_S_Py_S_ab;
  Double I_ERI_Gx3y_Px_Py_S_ab = I_ERI_H2x3y_S_Py_S_ab+ABX*I_ERI_Gx3y_S_Py_S_ab;
  Double I_ERI_Gx2yz_Px_Py_S_ab = I_ERI_H2x2yz_S_Py_S_ab+ABX*I_ERI_Gx2yz_S_Py_S_ab;
  Double I_ERI_Gxy2z_Px_Py_S_ab = I_ERI_H2xy2z_S_Py_S_ab+ABX*I_ERI_Gxy2z_S_Py_S_ab;
  Double I_ERI_Gx3z_Px_Py_S_ab = I_ERI_H2x3z_S_Py_S_ab+ABX*I_ERI_Gx3z_S_Py_S_ab;
  Double I_ERI_G4y_Px_Py_S_ab = I_ERI_Hx4y_S_Py_S_ab+ABX*I_ERI_G4y_S_Py_S_ab;
  Double I_ERI_G3yz_Px_Py_S_ab = I_ERI_Hx3yz_S_Py_S_ab+ABX*I_ERI_G3yz_S_Py_S_ab;
  Double I_ERI_G2y2z_Px_Py_S_ab = I_ERI_Hx2y2z_S_Py_S_ab+ABX*I_ERI_G2y2z_S_Py_S_ab;
  Double I_ERI_Gy3z_Px_Py_S_ab = I_ERI_Hxy3z_S_Py_S_ab+ABX*I_ERI_Gy3z_S_Py_S_ab;
  Double I_ERI_G4z_Px_Py_S_ab = I_ERI_Hx4z_S_Py_S_ab+ABX*I_ERI_G4z_S_Py_S_ab;
  Double I_ERI_G3xy_Py_Py_S_ab = I_ERI_H3x2y_S_Py_S_ab+ABY*I_ERI_G3xy_S_Py_S_ab;
  Double I_ERI_G3xz_Py_Py_S_ab = I_ERI_H3xyz_S_Py_S_ab+ABY*I_ERI_G3xz_S_Py_S_ab;
  Double I_ERI_G2x2y_Py_Py_S_ab = I_ERI_H2x3y_S_Py_S_ab+ABY*I_ERI_G2x2y_S_Py_S_ab;
  Double I_ERI_G2xyz_Py_Py_S_ab = I_ERI_H2x2yz_S_Py_S_ab+ABY*I_ERI_G2xyz_S_Py_S_ab;
  Double I_ERI_G2x2z_Py_Py_S_ab = I_ERI_H2xy2z_S_Py_S_ab+ABY*I_ERI_G2x2z_S_Py_S_ab;
  Double I_ERI_Gx3y_Py_Py_S_ab = I_ERI_Hx4y_S_Py_S_ab+ABY*I_ERI_Gx3y_S_Py_S_ab;
  Double I_ERI_Gx2yz_Py_Py_S_ab = I_ERI_Hx3yz_S_Py_S_ab+ABY*I_ERI_Gx2yz_S_Py_S_ab;
  Double I_ERI_Gxy2z_Py_Py_S_ab = I_ERI_Hx2y2z_S_Py_S_ab+ABY*I_ERI_Gxy2z_S_Py_S_ab;
  Double I_ERI_Gx3z_Py_Py_S_ab = I_ERI_Hxy3z_S_Py_S_ab+ABY*I_ERI_Gx3z_S_Py_S_ab;
  Double I_ERI_G4y_Py_Py_S_ab = I_ERI_H5y_S_Py_S_ab+ABY*I_ERI_G4y_S_Py_S_ab;
  Double I_ERI_G3yz_Py_Py_S_ab = I_ERI_H4yz_S_Py_S_ab+ABY*I_ERI_G3yz_S_Py_S_ab;
  Double I_ERI_G2y2z_Py_Py_S_ab = I_ERI_H3y2z_S_Py_S_ab+ABY*I_ERI_G2y2z_S_Py_S_ab;
  Double I_ERI_Gy3z_Py_Py_S_ab = I_ERI_H2y3z_S_Py_S_ab+ABY*I_ERI_Gy3z_S_Py_S_ab;
  Double I_ERI_G4z_Py_Py_S_ab = I_ERI_Hy4z_S_Py_S_ab+ABY*I_ERI_G4z_S_Py_S_ab;
  Double I_ERI_G3xz_Pz_Py_S_ab = I_ERI_H3x2z_S_Py_S_ab+ABZ*I_ERI_G3xz_S_Py_S_ab;
  Double I_ERI_G2xyz_Pz_Py_S_ab = I_ERI_H2xy2z_S_Py_S_ab+ABZ*I_ERI_G2xyz_S_Py_S_ab;
  Double I_ERI_G2x2z_Pz_Py_S_ab = I_ERI_H2x3z_S_Py_S_ab+ABZ*I_ERI_G2x2z_S_Py_S_ab;
  Double I_ERI_Gx2yz_Pz_Py_S_ab = I_ERI_Hx2y2z_S_Py_S_ab+ABZ*I_ERI_Gx2yz_S_Py_S_ab;
  Double I_ERI_Gxy2z_Pz_Py_S_ab = I_ERI_Hxy3z_S_Py_S_ab+ABZ*I_ERI_Gxy2z_S_Py_S_ab;
  Double I_ERI_Gx3z_Pz_Py_S_ab = I_ERI_Hx4z_S_Py_S_ab+ABZ*I_ERI_Gx3z_S_Py_S_ab;
  Double I_ERI_G3yz_Pz_Py_S_ab = I_ERI_H3y2z_S_Py_S_ab+ABZ*I_ERI_G3yz_S_Py_S_ab;
  Double I_ERI_G2y2z_Pz_Py_S_ab = I_ERI_H2y3z_S_Py_S_ab+ABZ*I_ERI_G2y2z_S_Py_S_ab;
  Double I_ERI_Gy3z_Pz_Py_S_ab = I_ERI_Hy4z_S_Py_S_ab+ABZ*I_ERI_Gy3z_S_Py_S_ab;
  Double I_ERI_G4z_Pz_Py_S_ab = I_ERI_H5z_S_Py_S_ab+ABZ*I_ERI_G4z_S_Py_S_ab;
  Double I_ERI_G4x_Px_Pz_S_ab = I_ERI_H5x_S_Pz_S_ab+ABX*I_ERI_G4x_S_Pz_S_ab;
  Double I_ERI_G3xy_Px_Pz_S_ab = I_ERI_H4xy_S_Pz_S_ab+ABX*I_ERI_G3xy_S_Pz_S_ab;
  Double I_ERI_G3xz_Px_Pz_S_ab = I_ERI_H4xz_S_Pz_S_ab+ABX*I_ERI_G3xz_S_Pz_S_ab;
  Double I_ERI_G2x2y_Px_Pz_S_ab = I_ERI_H3x2y_S_Pz_S_ab+ABX*I_ERI_G2x2y_S_Pz_S_ab;
  Double I_ERI_G2xyz_Px_Pz_S_ab = I_ERI_H3xyz_S_Pz_S_ab+ABX*I_ERI_G2xyz_S_Pz_S_ab;
  Double I_ERI_G2x2z_Px_Pz_S_ab = I_ERI_H3x2z_S_Pz_S_ab+ABX*I_ERI_G2x2z_S_Pz_S_ab;
  Double I_ERI_Gx3y_Px_Pz_S_ab = I_ERI_H2x3y_S_Pz_S_ab+ABX*I_ERI_Gx3y_S_Pz_S_ab;
  Double I_ERI_Gx2yz_Px_Pz_S_ab = I_ERI_H2x2yz_S_Pz_S_ab+ABX*I_ERI_Gx2yz_S_Pz_S_ab;
  Double I_ERI_Gxy2z_Px_Pz_S_ab = I_ERI_H2xy2z_S_Pz_S_ab+ABX*I_ERI_Gxy2z_S_Pz_S_ab;
  Double I_ERI_Gx3z_Px_Pz_S_ab = I_ERI_H2x3z_S_Pz_S_ab+ABX*I_ERI_Gx3z_S_Pz_S_ab;
  Double I_ERI_G4y_Px_Pz_S_ab = I_ERI_Hx4y_S_Pz_S_ab+ABX*I_ERI_G4y_S_Pz_S_ab;
  Double I_ERI_G3yz_Px_Pz_S_ab = I_ERI_Hx3yz_S_Pz_S_ab+ABX*I_ERI_G3yz_S_Pz_S_ab;
  Double I_ERI_G2y2z_Px_Pz_S_ab = I_ERI_Hx2y2z_S_Pz_S_ab+ABX*I_ERI_G2y2z_S_Pz_S_ab;
  Double I_ERI_Gy3z_Px_Pz_S_ab = I_ERI_Hxy3z_S_Pz_S_ab+ABX*I_ERI_Gy3z_S_Pz_S_ab;
  Double I_ERI_G4z_Px_Pz_S_ab = I_ERI_Hx4z_S_Pz_S_ab+ABX*I_ERI_G4z_S_Pz_S_ab;
  Double I_ERI_G3xy_Py_Pz_S_ab = I_ERI_H3x2y_S_Pz_S_ab+ABY*I_ERI_G3xy_S_Pz_S_ab;
  Double I_ERI_G3xz_Py_Pz_S_ab = I_ERI_H3xyz_S_Pz_S_ab+ABY*I_ERI_G3xz_S_Pz_S_ab;
  Double I_ERI_G2x2y_Py_Pz_S_ab = I_ERI_H2x3y_S_Pz_S_ab+ABY*I_ERI_G2x2y_S_Pz_S_ab;
  Double I_ERI_G2xyz_Py_Pz_S_ab = I_ERI_H2x2yz_S_Pz_S_ab+ABY*I_ERI_G2xyz_S_Pz_S_ab;
  Double I_ERI_G2x2z_Py_Pz_S_ab = I_ERI_H2xy2z_S_Pz_S_ab+ABY*I_ERI_G2x2z_S_Pz_S_ab;
  Double I_ERI_Gx3y_Py_Pz_S_ab = I_ERI_Hx4y_S_Pz_S_ab+ABY*I_ERI_Gx3y_S_Pz_S_ab;
  Double I_ERI_Gx2yz_Py_Pz_S_ab = I_ERI_Hx3yz_S_Pz_S_ab+ABY*I_ERI_Gx2yz_S_Pz_S_ab;
  Double I_ERI_Gxy2z_Py_Pz_S_ab = I_ERI_Hx2y2z_S_Pz_S_ab+ABY*I_ERI_Gxy2z_S_Pz_S_ab;
  Double I_ERI_Gx3z_Py_Pz_S_ab = I_ERI_Hxy3z_S_Pz_S_ab+ABY*I_ERI_Gx3z_S_Pz_S_ab;
  Double I_ERI_G4y_Py_Pz_S_ab = I_ERI_H5y_S_Pz_S_ab+ABY*I_ERI_G4y_S_Pz_S_ab;
  Double I_ERI_G3yz_Py_Pz_S_ab = I_ERI_H4yz_S_Pz_S_ab+ABY*I_ERI_G3yz_S_Pz_S_ab;
  Double I_ERI_G2y2z_Py_Pz_S_ab = I_ERI_H3y2z_S_Pz_S_ab+ABY*I_ERI_G2y2z_S_Pz_S_ab;
  Double I_ERI_Gy3z_Py_Pz_S_ab = I_ERI_H2y3z_S_Pz_S_ab+ABY*I_ERI_Gy3z_S_Pz_S_ab;
  Double I_ERI_G4z_Py_Pz_S_ab = I_ERI_Hy4z_S_Pz_S_ab+ABY*I_ERI_G4z_S_Pz_S_ab;
  Double I_ERI_G3xz_Pz_Pz_S_ab = I_ERI_H3x2z_S_Pz_S_ab+ABZ*I_ERI_G3xz_S_Pz_S_ab;
  Double I_ERI_G2xyz_Pz_Pz_S_ab = I_ERI_H2xy2z_S_Pz_S_ab+ABZ*I_ERI_G2xyz_S_Pz_S_ab;
  Double I_ERI_G2x2z_Pz_Pz_S_ab = I_ERI_H2x3z_S_Pz_S_ab+ABZ*I_ERI_G2x2z_S_Pz_S_ab;
  Double I_ERI_Gx2yz_Pz_Pz_S_ab = I_ERI_Hx2y2z_S_Pz_S_ab+ABZ*I_ERI_Gx2yz_S_Pz_S_ab;
  Double I_ERI_Gxy2z_Pz_Pz_S_ab = I_ERI_Hxy3z_S_Pz_S_ab+ABZ*I_ERI_Gxy2z_S_Pz_S_ab;
  Double I_ERI_Gx3z_Pz_Pz_S_ab = I_ERI_Hx4z_S_Pz_S_ab+ABZ*I_ERI_Gx3z_S_Pz_S_ab;
  Double I_ERI_G3yz_Pz_Pz_S_ab = I_ERI_H3y2z_S_Pz_S_ab+ABZ*I_ERI_G3yz_S_Pz_S_ab;
  Double I_ERI_G2y2z_Pz_Pz_S_ab = I_ERI_H2y3z_S_Pz_S_ab+ABZ*I_ERI_G2y2z_S_Pz_S_ab;
  Double I_ERI_Gy3z_Pz_Pz_S_ab = I_ERI_Hy4z_S_Pz_S_ab+ABZ*I_ERI_Gy3z_S_Pz_S_ab;
  Double I_ERI_G4z_Pz_Pz_S_ab = I_ERI_H5z_S_Pz_S_ab+ABZ*I_ERI_G4z_S_Pz_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_P_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_P_P_S_ab
   ************************************************************/
  Double I_ERI_F3x_D2x_Px_S_ab = I_ERI_G4x_Px_Px_S_ab+ABX*I_ERI_F3x_Px_Px_S_ab;
  Double I_ERI_F2xy_D2x_Px_S_ab = I_ERI_G3xy_Px_Px_S_ab+ABX*I_ERI_F2xy_Px_Px_S_ab;
  Double I_ERI_F2xz_D2x_Px_S_ab = I_ERI_G3xz_Px_Px_S_ab+ABX*I_ERI_F2xz_Px_Px_S_ab;
  Double I_ERI_Fx2y_D2x_Px_S_ab = I_ERI_G2x2y_Px_Px_S_ab+ABX*I_ERI_Fx2y_Px_Px_S_ab;
  Double I_ERI_Fxyz_D2x_Px_S_ab = I_ERI_G2xyz_Px_Px_S_ab+ABX*I_ERI_Fxyz_Px_Px_S_ab;
  Double I_ERI_Fx2z_D2x_Px_S_ab = I_ERI_G2x2z_Px_Px_S_ab+ABX*I_ERI_Fx2z_Px_Px_S_ab;
  Double I_ERI_F3y_D2x_Px_S_ab = I_ERI_Gx3y_Px_Px_S_ab+ABX*I_ERI_F3y_Px_Px_S_ab;
  Double I_ERI_F2yz_D2x_Px_S_ab = I_ERI_Gx2yz_Px_Px_S_ab+ABX*I_ERI_F2yz_Px_Px_S_ab;
  Double I_ERI_Fy2z_D2x_Px_S_ab = I_ERI_Gxy2z_Px_Px_S_ab+ABX*I_ERI_Fy2z_Px_Px_S_ab;
  Double I_ERI_F3z_D2x_Px_S_ab = I_ERI_Gx3z_Px_Px_S_ab+ABX*I_ERI_F3z_Px_Px_S_ab;
  Double I_ERI_F3x_Dxy_Px_S_ab = I_ERI_G3xy_Px_Px_S_ab+ABY*I_ERI_F3x_Px_Px_S_ab;
  Double I_ERI_F2xy_Dxy_Px_S_ab = I_ERI_G2x2y_Px_Px_S_ab+ABY*I_ERI_F2xy_Px_Px_S_ab;
  Double I_ERI_F2xz_Dxy_Px_S_ab = I_ERI_G2xyz_Px_Px_S_ab+ABY*I_ERI_F2xz_Px_Px_S_ab;
  Double I_ERI_Fx2y_Dxy_Px_S_ab = I_ERI_Gx3y_Px_Px_S_ab+ABY*I_ERI_Fx2y_Px_Px_S_ab;
  Double I_ERI_Fxyz_Dxy_Px_S_ab = I_ERI_Gx2yz_Px_Px_S_ab+ABY*I_ERI_Fxyz_Px_Px_S_ab;
  Double I_ERI_Fx2z_Dxy_Px_S_ab = I_ERI_Gxy2z_Px_Px_S_ab+ABY*I_ERI_Fx2z_Px_Px_S_ab;
  Double I_ERI_F3y_Dxy_Px_S_ab = I_ERI_G4y_Px_Px_S_ab+ABY*I_ERI_F3y_Px_Px_S_ab;
  Double I_ERI_F2yz_Dxy_Px_S_ab = I_ERI_G3yz_Px_Px_S_ab+ABY*I_ERI_F2yz_Px_Px_S_ab;
  Double I_ERI_Fy2z_Dxy_Px_S_ab = I_ERI_G2y2z_Px_Px_S_ab+ABY*I_ERI_Fy2z_Px_Px_S_ab;
  Double I_ERI_F3z_Dxy_Px_S_ab = I_ERI_Gy3z_Px_Px_S_ab+ABY*I_ERI_F3z_Px_Px_S_ab;
  Double I_ERI_F3x_Dxz_Px_S_ab = I_ERI_G3xz_Px_Px_S_ab+ABZ*I_ERI_F3x_Px_Px_S_ab;
  Double I_ERI_F2xy_Dxz_Px_S_ab = I_ERI_G2xyz_Px_Px_S_ab+ABZ*I_ERI_F2xy_Px_Px_S_ab;
  Double I_ERI_F2xz_Dxz_Px_S_ab = I_ERI_G2x2z_Px_Px_S_ab+ABZ*I_ERI_F2xz_Px_Px_S_ab;
  Double I_ERI_Fx2y_Dxz_Px_S_ab = I_ERI_Gx2yz_Px_Px_S_ab+ABZ*I_ERI_Fx2y_Px_Px_S_ab;
  Double I_ERI_Fxyz_Dxz_Px_S_ab = I_ERI_Gxy2z_Px_Px_S_ab+ABZ*I_ERI_Fxyz_Px_Px_S_ab;
  Double I_ERI_Fx2z_Dxz_Px_S_ab = I_ERI_Gx3z_Px_Px_S_ab+ABZ*I_ERI_Fx2z_Px_Px_S_ab;
  Double I_ERI_F3y_Dxz_Px_S_ab = I_ERI_G3yz_Px_Px_S_ab+ABZ*I_ERI_F3y_Px_Px_S_ab;
  Double I_ERI_F2yz_Dxz_Px_S_ab = I_ERI_G2y2z_Px_Px_S_ab+ABZ*I_ERI_F2yz_Px_Px_S_ab;
  Double I_ERI_Fy2z_Dxz_Px_S_ab = I_ERI_Gy3z_Px_Px_S_ab+ABZ*I_ERI_Fy2z_Px_Px_S_ab;
  Double I_ERI_F3z_Dxz_Px_S_ab = I_ERI_G4z_Px_Px_S_ab+ABZ*I_ERI_F3z_Px_Px_S_ab;
  Double I_ERI_F3x_D2y_Px_S_ab = I_ERI_G3xy_Py_Px_S_ab+ABY*I_ERI_F3x_Py_Px_S_ab;
  Double I_ERI_F2xy_D2y_Px_S_ab = I_ERI_G2x2y_Py_Px_S_ab+ABY*I_ERI_F2xy_Py_Px_S_ab;
  Double I_ERI_F2xz_D2y_Px_S_ab = I_ERI_G2xyz_Py_Px_S_ab+ABY*I_ERI_F2xz_Py_Px_S_ab;
  Double I_ERI_Fx2y_D2y_Px_S_ab = I_ERI_Gx3y_Py_Px_S_ab+ABY*I_ERI_Fx2y_Py_Px_S_ab;
  Double I_ERI_Fxyz_D2y_Px_S_ab = I_ERI_Gx2yz_Py_Px_S_ab+ABY*I_ERI_Fxyz_Py_Px_S_ab;
  Double I_ERI_Fx2z_D2y_Px_S_ab = I_ERI_Gxy2z_Py_Px_S_ab+ABY*I_ERI_Fx2z_Py_Px_S_ab;
  Double I_ERI_F3y_D2y_Px_S_ab = I_ERI_G4y_Py_Px_S_ab+ABY*I_ERI_F3y_Py_Px_S_ab;
  Double I_ERI_F2yz_D2y_Px_S_ab = I_ERI_G3yz_Py_Px_S_ab+ABY*I_ERI_F2yz_Py_Px_S_ab;
  Double I_ERI_Fy2z_D2y_Px_S_ab = I_ERI_G2y2z_Py_Px_S_ab+ABY*I_ERI_Fy2z_Py_Px_S_ab;
  Double I_ERI_F3z_D2y_Px_S_ab = I_ERI_Gy3z_Py_Px_S_ab+ABY*I_ERI_F3z_Py_Px_S_ab;
  Double I_ERI_F3x_Dyz_Px_S_ab = I_ERI_G3xz_Py_Px_S_ab+ABZ*I_ERI_F3x_Py_Px_S_ab;
  Double I_ERI_F2xy_Dyz_Px_S_ab = I_ERI_G2xyz_Py_Px_S_ab+ABZ*I_ERI_F2xy_Py_Px_S_ab;
  Double I_ERI_F2xz_Dyz_Px_S_ab = I_ERI_G2x2z_Py_Px_S_ab+ABZ*I_ERI_F2xz_Py_Px_S_ab;
  Double I_ERI_Fx2y_Dyz_Px_S_ab = I_ERI_Gx2yz_Py_Px_S_ab+ABZ*I_ERI_Fx2y_Py_Px_S_ab;
  Double I_ERI_Fxyz_Dyz_Px_S_ab = I_ERI_Gxy2z_Py_Px_S_ab+ABZ*I_ERI_Fxyz_Py_Px_S_ab;
  Double I_ERI_Fx2z_Dyz_Px_S_ab = I_ERI_Gx3z_Py_Px_S_ab+ABZ*I_ERI_Fx2z_Py_Px_S_ab;
  Double I_ERI_F3y_Dyz_Px_S_ab = I_ERI_G3yz_Py_Px_S_ab+ABZ*I_ERI_F3y_Py_Px_S_ab;
  Double I_ERI_F2yz_Dyz_Px_S_ab = I_ERI_G2y2z_Py_Px_S_ab+ABZ*I_ERI_F2yz_Py_Px_S_ab;
  Double I_ERI_Fy2z_Dyz_Px_S_ab = I_ERI_Gy3z_Py_Px_S_ab+ABZ*I_ERI_Fy2z_Py_Px_S_ab;
  Double I_ERI_F3z_Dyz_Px_S_ab = I_ERI_G4z_Py_Px_S_ab+ABZ*I_ERI_F3z_Py_Px_S_ab;
  Double I_ERI_F3x_D2z_Px_S_ab = I_ERI_G3xz_Pz_Px_S_ab+ABZ*I_ERI_F3x_Pz_Px_S_ab;
  Double I_ERI_F2xy_D2z_Px_S_ab = I_ERI_G2xyz_Pz_Px_S_ab+ABZ*I_ERI_F2xy_Pz_Px_S_ab;
  Double I_ERI_F2xz_D2z_Px_S_ab = I_ERI_G2x2z_Pz_Px_S_ab+ABZ*I_ERI_F2xz_Pz_Px_S_ab;
  Double I_ERI_Fx2y_D2z_Px_S_ab = I_ERI_Gx2yz_Pz_Px_S_ab+ABZ*I_ERI_Fx2y_Pz_Px_S_ab;
  Double I_ERI_Fxyz_D2z_Px_S_ab = I_ERI_Gxy2z_Pz_Px_S_ab+ABZ*I_ERI_Fxyz_Pz_Px_S_ab;
  Double I_ERI_Fx2z_D2z_Px_S_ab = I_ERI_Gx3z_Pz_Px_S_ab+ABZ*I_ERI_Fx2z_Pz_Px_S_ab;
  Double I_ERI_F3y_D2z_Px_S_ab = I_ERI_G3yz_Pz_Px_S_ab+ABZ*I_ERI_F3y_Pz_Px_S_ab;
  Double I_ERI_F2yz_D2z_Px_S_ab = I_ERI_G2y2z_Pz_Px_S_ab+ABZ*I_ERI_F2yz_Pz_Px_S_ab;
  Double I_ERI_Fy2z_D2z_Px_S_ab = I_ERI_Gy3z_Pz_Px_S_ab+ABZ*I_ERI_Fy2z_Pz_Px_S_ab;
  Double I_ERI_F3z_D2z_Px_S_ab = I_ERI_G4z_Pz_Px_S_ab+ABZ*I_ERI_F3z_Pz_Px_S_ab;
  Double I_ERI_F3x_D2x_Py_S_ab = I_ERI_G4x_Px_Py_S_ab+ABX*I_ERI_F3x_Px_Py_S_ab;
  Double I_ERI_F2xy_D2x_Py_S_ab = I_ERI_G3xy_Px_Py_S_ab+ABX*I_ERI_F2xy_Px_Py_S_ab;
  Double I_ERI_F2xz_D2x_Py_S_ab = I_ERI_G3xz_Px_Py_S_ab+ABX*I_ERI_F2xz_Px_Py_S_ab;
  Double I_ERI_Fx2y_D2x_Py_S_ab = I_ERI_G2x2y_Px_Py_S_ab+ABX*I_ERI_Fx2y_Px_Py_S_ab;
  Double I_ERI_Fxyz_D2x_Py_S_ab = I_ERI_G2xyz_Px_Py_S_ab+ABX*I_ERI_Fxyz_Px_Py_S_ab;
  Double I_ERI_Fx2z_D2x_Py_S_ab = I_ERI_G2x2z_Px_Py_S_ab+ABX*I_ERI_Fx2z_Px_Py_S_ab;
  Double I_ERI_F3y_D2x_Py_S_ab = I_ERI_Gx3y_Px_Py_S_ab+ABX*I_ERI_F3y_Px_Py_S_ab;
  Double I_ERI_F2yz_D2x_Py_S_ab = I_ERI_Gx2yz_Px_Py_S_ab+ABX*I_ERI_F2yz_Px_Py_S_ab;
  Double I_ERI_Fy2z_D2x_Py_S_ab = I_ERI_Gxy2z_Px_Py_S_ab+ABX*I_ERI_Fy2z_Px_Py_S_ab;
  Double I_ERI_F3z_D2x_Py_S_ab = I_ERI_Gx3z_Px_Py_S_ab+ABX*I_ERI_F3z_Px_Py_S_ab;
  Double I_ERI_F3x_Dxy_Py_S_ab = I_ERI_G3xy_Px_Py_S_ab+ABY*I_ERI_F3x_Px_Py_S_ab;
  Double I_ERI_F2xy_Dxy_Py_S_ab = I_ERI_G2x2y_Px_Py_S_ab+ABY*I_ERI_F2xy_Px_Py_S_ab;
  Double I_ERI_F2xz_Dxy_Py_S_ab = I_ERI_G2xyz_Px_Py_S_ab+ABY*I_ERI_F2xz_Px_Py_S_ab;
  Double I_ERI_Fx2y_Dxy_Py_S_ab = I_ERI_Gx3y_Px_Py_S_ab+ABY*I_ERI_Fx2y_Px_Py_S_ab;
  Double I_ERI_Fxyz_Dxy_Py_S_ab = I_ERI_Gx2yz_Px_Py_S_ab+ABY*I_ERI_Fxyz_Px_Py_S_ab;
  Double I_ERI_Fx2z_Dxy_Py_S_ab = I_ERI_Gxy2z_Px_Py_S_ab+ABY*I_ERI_Fx2z_Px_Py_S_ab;
  Double I_ERI_F3y_Dxy_Py_S_ab = I_ERI_G4y_Px_Py_S_ab+ABY*I_ERI_F3y_Px_Py_S_ab;
  Double I_ERI_F2yz_Dxy_Py_S_ab = I_ERI_G3yz_Px_Py_S_ab+ABY*I_ERI_F2yz_Px_Py_S_ab;
  Double I_ERI_Fy2z_Dxy_Py_S_ab = I_ERI_G2y2z_Px_Py_S_ab+ABY*I_ERI_Fy2z_Px_Py_S_ab;
  Double I_ERI_F3z_Dxy_Py_S_ab = I_ERI_Gy3z_Px_Py_S_ab+ABY*I_ERI_F3z_Px_Py_S_ab;
  Double I_ERI_F3x_Dxz_Py_S_ab = I_ERI_G3xz_Px_Py_S_ab+ABZ*I_ERI_F3x_Px_Py_S_ab;
  Double I_ERI_F2xy_Dxz_Py_S_ab = I_ERI_G2xyz_Px_Py_S_ab+ABZ*I_ERI_F2xy_Px_Py_S_ab;
  Double I_ERI_F2xz_Dxz_Py_S_ab = I_ERI_G2x2z_Px_Py_S_ab+ABZ*I_ERI_F2xz_Px_Py_S_ab;
  Double I_ERI_Fx2y_Dxz_Py_S_ab = I_ERI_Gx2yz_Px_Py_S_ab+ABZ*I_ERI_Fx2y_Px_Py_S_ab;
  Double I_ERI_Fxyz_Dxz_Py_S_ab = I_ERI_Gxy2z_Px_Py_S_ab+ABZ*I_ERI_Fxyz_Px_Py_S_ab;
  Double I_ERI_Fx2z_Dxz_Py_S_ab = I_ERI_Gx3z_Px_Py_S_ab+ABZ*I_ERI_Fx2z_Px_Py_S_ab;
  Double I_ERI_F3y_Dxz_Py_S_ab = I_ERI_G3yz_Px_Py_S_ab+ABZ*I_ERI_F3y_Px_Py_S_ab;
  Double I_ERI_F2yz_Dxz_Py_S_ab = I_ERI_G2y2z_Px_Py_S_ab+ABZ*I_ERI_F2yz_Px_Py_S_ab;
  Double I_ERI_Fy2z_Dxz_Py_S_ab = I_ERI_Gy3z_Px_Py_S_ab+ABZ*I_ERI_Fy2z_Px_Py_S_ab;
  Double I_ERI_F3z_Dxz_Py_S_ab = I_ERI_G4z_Px_Py_S_ab+ABZ*I_ERI_F3z_Px_Py_S_ab;
  Double I_ERI_F3x_D2y_Py_S_ab = I_ERI_G3xy_Py_Py_S_ab+ABY*I_ERI_F3x_Py_Py_S_ab;
  Double I_ERI_F2xy_D2y_Py_S_ab = I_ERI_G2x2y_Py_Py_S_ab+ABY*I_ERI_F2xy_Py_Py_S_ab;
  Double I_ERI_F2xz_D2y_Py_S_ab = I_ERI_G2xyz_Py_Py_S_ab+ABY*I_ERI_F2xz_Py_Py_S_ab;
  Double I_ERI_Fx2y_D2y_Py_S_ab = I_ERI_Gx3y_Py_Py_S_ab+ABY*I_ERI_Fx2y_Py_Py_S_ab;
  Double I_ERI_Fxyz_D2y_Py_S_ab = I_ERI_Gx2yz_Py_Py_S_ab+ABY*I_ERI_Fxyz_Py_Py_S_ab;
  Double I_ERI_Fx2z_D2y_Py_S_ab = I_ERI_Gxy2z_Py_Py_S_ab+ABY*I_ERI_Fx2z_Py_Py_S_ab;
  Double I_ERI_F3y_D2y_Py_S_ab = I_ERI_G4y_Py_Py_S_ab+ABY*I_ERI_F3y_Py_Py_S_ab;
  Double I_ERI_F2yz_D2y_Py_S_ab = I_ERI_G3yz_Py_Py_S_ab+ABY*I_ERI_F2yz_Py_Py_S_ab;
  Double I_ERI_Fy2z_D2y_Py_S_ab = I_ERI_G2y2z_Py_Py_S_ab+ABY*I_ERI_Fy2z_Py_Py_S_ab;
  Double I_ERI_F3z_D2y_Py_S_ab = I_ERI_Gy3z_Py_Py_S_ab+ABY*I_ERI_F3z_Py_Py_S_ab;
  Double I_ERI_F3x_Dyz_Py_S_ab = I_ERI_G3xz_Py_Py_S_ab+ABZ*I_ERI_F3x_Py_Py_S_ab;
  Double I_ERI_F2xy_Dyz_Py_S_ab = I_ERI_G2xyz_Py_Py_S_ab+ABZ*I_ERI_F2xy_Py_Py_S_ab;
  Double I_ERI_F2xz_Dyz_Py_S_ab = I_ERI_G2x2z_Py_Py_S_ab+ABZ*I_ERI_F2xz_Py_Py_S_ab;
  Double I_ERI_Fx2y_Dyz_Py_S_ab = I_ERI_Gx2yz_Py_Py_S_ab+ABZ*I_ERI_Fx2y_Py_Py_S_ab;
  Double I_ERI_Fxyz_Dyz_Py_S_ab = I_ERI_Gxy2z_Py_Py_S_ab+ABZ*I_ERI_Fxyz_Py_Py_S_ab;
  Double I_ERI_Fx2z_Dyz_Py_S_ab = I_ERI_Gx3z_Py_Py_S_ab+ABZ*I_ERI_Fx2z_Py_Py_S_ab;
  Double I_ERI_F3y_Dyz_Py_S_ab = I_ERI_G3yz_Py_Py_S_ab+ABZ*I_ERI_F3y_Py_Py_S_ab;
  Double I_ERI_F2yz_Dyz_Py_S_ab = I_ERI_G2y2z_Py_Py_S_ab+ABZ*I_ERI_F2yz_Py_Py_S_ab;
  Double I_ERI_Fy2z_Dyz_Py_S_ab = I_ERI_Gy3z_Py_Py_S_ab+ABZ*I_ERI_Fy2z_Py_Py_S_ab;
  Double I_ERI_F3z_Dyz_Py_S_ab = I_ERI_G4z_Py_Py_S_ab+ABZ*I_ERI_F3z_Py_Py_S_ab;
  Double I_ERI_F3x_D2z_Py_S_ab = I_ERI_G3xz_Pz_Py_S_ab+ABZ*I_ERI_F3x_Pz_Py_S_ab;
  Double I_ERI_F2xy_D2z_Py_S_ab = I_ERI_G2xyz_Pz_Py_S_ab+ABZ*I_ERI_F2xy_Pz_Py_S_ab;
  Double I_ERI_F2xz_D2z_Py_S_ab = I_ERI_G2x2z_Pz_Py_S_ab+ABZ*I_ERI_F2xz_Pz_Py_S_ab;
  Double I_ERI_Fx2y_D2z_Py_S_ab = I_ERI_Gx2yz_Pz_Py_S_ab+ABZ*I_ERI_Fx2y_Pz_Py_S_ab;
  Double I_ERI_Fxyz_D2z_Py_S_ab = I_ERI_Gxy2z_Pz_Py_S_ab+ABZ*I_ERI_Fxyz_Pz_Py_S_ab;
  Double I_ERI_Fx2z_D2z_Py_S_ab = I_ERI_Gx3z_Pz_Py_S_ab+ABZ*I_ERI_Fx2z_Pz_Py_S_ab;
  Double I_ERI_F3y_D2z_Py_S_ab = I_ERI_G3yz_Pz_Py_S_ab+ABZ*I_ERI_F3y_Pz_Py_S_ab;
  Double I_ERI_F2yz_D2z_Py_S_ab = I_ERI_G2y2z_Pz_Py_S_ab+ABZ*I_ERI_F2yz_Pz_Py_S_ab;
  Double I_ERI_Fy2z_D2z_Py_S_ab = I_ERI_Gy3z_Pz_Py_S_ab+ABZ*I_ERI_Fy2z_Pz_Py_S_ab;
  Double I_ERI_F3z_D2z_Py_S_ab = I_ERI_G4z_Pz_Py_S_ab+ABZ*I_ERI_F3z_Pz_Py_S_ab;
  Double I_ERI_F3x_D2x_Pz_S_ab = I_ERI_G4x_Px_Pz_S_ab+ABX*I_ERI_F3x_Px_Pz_S_ab;
  Double I_ERI_F2xy_D2x_Pz_S_ab = I_ERI_G3xy_Px_Pz_S_ab+ABX*I_ERI_F2xy_Px_Pz_S_ab;
  Double I_ERI_F2xz_D2x_Pz_S_ab = I_ERI_G3xz_Px_Pz_S_ab+ABX*I_ERI_F2xz_Px_Pz_S_ab;
  Double I_ERI_Fx2y_D2x_Pz_S_ab = I_ERI_G2x2y_Px_Pz_S_ab+ABX*I_ERI_Fx2y_Px_Pz_S_ab;
  Double I_ERI_Fxyz_D2x_Pz_S_ab = I_ERI_G2xyz_Px_Pz_S_ab+ABX*I_ERI_Fxyz_Px_Pz_S_ab;
  Double I_ERI_Fx2z_D2x_Pz_S_ab = I_ERI_G2x2z_Px_Pz_S_ab+ABX*I_ERI_Fx2z_Px_Pz_S_ab;
  Double I_ERI_F3y_D2x_Pz_S_ab = I_ERI_Gx3y_Px_Pz_S_ab+ABX*I_ERI_F3y_Px_Pz_S_ab;
  Double I_ERI_F2yz_D2x_Pz_S_ab = I_ERI_Gx2yz_Px_Pz_S_ab+ABX*I_ERI_F2yz_Px_Pz_S_ab;
  Double I_ERI_Fy2z_D2x_Pz_S_ab = I_ERI_Gxy2z_Px_Pz_S_ab+ABX*I_ERI_Fy2z_Px_Pz_S_ab;
  Double I_ERI_F3z_D2x_Pz_S_ab = I_ERI_Gx3z_Px_Pz_S_ab+ABX*I_ERI_F3z_Px_Pz_S_ab;
  Double I_ERI_F3x_Dxy_Pz_S_ab = I_ERI_G3xy_Px_Pz_S_ab+ABY*I_ERI_F3x_Px_Pz_S_ab;
  Double I_ERI_F2xy_Dxy_Pz_S_ab = I_ERI_G2x2y_Px_Pz_S_ab+ABY*I_ERI_F2xy_Px_Pz_S_ab;
  Double I_ERI_F2xz_Dxy_Pz_S_ab = I_ERI_G2xyz_Px_Pz_S_ab+ABY*I_ERI_F2xz_Px_Pz_S_ab;
  Double I_ERI_Fx2y_Dxy_Pz_S_ab = I_ERI_Gx3y_Px_Pz_S_ab+ABY*I_ERI_Fx2y_Px_Pz_S_ab;
  Double I_ERI_Fxyz_Dxy_Pz_S_ab = I_ERI_Gx2yz_Px_Pz_S_ab+ABY*I_ERI_Fxyz_Px_Pz_S_ab;
  Double I_ERI_Fx2z_Dxy_Pz_S_ab = I_ERI_Gxy2z_Px_Pz_S_ab+ABY*I_ERI_Fx2z_Px_Pz_S_ab;
  Double I_ERI_F3y_Dxy_Pz_S_ab = I_ERI_G4y_Px_Pz_S_ab+ABY*I_ERI_F3y_Px_Pz_S_ab;
  Double I_ERI_F2yz_Dxy_Pz_S_ab = I_ERI_G3yz_Px_Pz_S_ab+ABY*I_ERI_F2yz_Px_Pz_S_ab;
  Double I_ERI_Fy2z_Dxy_Pz_S_ab = I_ERI_G2y2z_Px_Pz_S_ab+ABY*I_ERI_Fy2z_Px_Pz_S_ab;
  Double I_ERI_F3z_Dxy_Pz_S_ab = I_ERI_Gy3z_Px_Pz_S_ab+ABY*I_ERI_F3z_Px_Pz_S_ab;
  Double I_ERI_F3x_Dxz_Pz_S_ab = I_ERI_G3xz_Px_Pz_S_ab+ABZ*I_ERI_F3x_Px_Pz_S_ab;
  Double I_ERI_F2xy_Dxz_Pz_S_ab = I_ERI_G2xyz_Px_Pz_S_ab+ABZ*I_ERI_F2xy_Px_Pz_S_ab;
  Double I_ERI_F2xz_Dxz_Pz_S_ab = I_ERI_G2x2z_Px_Pz_S_ab+ABZ*I_ERI_F2xz_Px_Pz_S_ab;
  Double I_ERI_Fx2y_Dxz_Pz_S_ab = I_ERI_Gx2yz_Px_Pz_S_ab+ABZ*I_ERI_Fx2y_Px_Pz_S_ab;
  Double I_ERI_Fxyz_Dxz_Pz_S_ab = I_ERI_Gxy2z_Px_Pz_S_ab+ABZ*I_ERI_Fxyz_Px_Pz_S_ab;
  Double I_ERI_Fx2z_Dxz_Pz_S_ab = I_ERI_Gx3z_Px_Pz_S_ab+ABZ*I_ERI_Fx2z_Px_Pz_S_ab;
  Double I_ERI_F3y_Dxz_Pz_S_ab = I_ERI_G3yz_Px_Pz_S_ab+ABZ*I_ERI_F3y_Px_Pz_S_ab;
  Double I_ERI_F2yz_Dxz_Pz_S_ab = I_ERI_G2y2z_Px_Pz_S_ab+ABZ*I_ERI_F2yz_Px_Pz_S_ab;
  Double I_ERI_Fy2z_Dxz_Pz_S_ab = I_ERI_Gy3z_Px_Pz_S_ab+ABZ*I_ERI_Fy2z_Px_Pz_S_ab;
  Double I_ERI_F3z_Dxz_Pz_S_ab = I_ERI_G4z_Px_Pz_S_ab+ABZ*I_ERI_F3z_Px_Pz_S_ab;
  Double I_ERI_F3x_D2y_Pz_S_ab = I_ERI_G3xy_Py_Pz_S_ab+ABY*I_ERI_F3x_Py_Pz_S_ab;
  Double I_ERI_F2xy_D2y_Pz_S_ab = I_ERI_G2x2y_Py_Pz_S_ab+ABY*I_ERI_F2xy_Py_Pz_S_ab;
  Double I_ERI_F2xz_D2y_Pz_S_ab = I_ERI_G2xyz_Py_Pz_S_ab+ABY*I_ERI_F2xz_Py_Pz_S_ab;
  Double I_ERI_Fx2y_D2y_Pz_S_ab = I_ERI_Gx3y_Py_Pz_S_ab+ABY*I_ERI_Fx2y_Py_Pz_S_ab;
  Double I_ERI_Fxyz_D2y_Pz_S_ab = I_ERI_Gx2yz_Py_Pz_S_ab+ABY*I_ERI_Fxyz_Py_Pz_S_ab;
  Double I_ERI_Fx2z_D2y_Pz_S_ab = I_ERI_Gxy2z_Py_Pz_S_ab+ABY*I_ERI_Fx2z_Py_Pz_S_ab;
  Double I_ERI_F3y_D2y_Pz_S_ab = I_ERI_G4y_Py_Pz_S_ab+ABY*I_ERI_F3y_Py_Pz_S_ab;
  Double I_ERI_F2yz_D2y_Pz_S_ab = I_ERI_G3yz_Py_Pz_S_ab+ABY*I_ERI_F2yz_Py_Pz_S_ab;
  Double I_ERI_Fy2z_D2y_Pz_S_ab = I_ERI_G2y2z_Py_Pz_S_ab+ABY*I_ERI_Fy2z_Py_Pz_S_ab;
  Double I_ERI_F3z_D2y_Pz_S_ab = I_ERI_Gy3z_Py_Pz_S_ab+ABY*I_ERI_F3z_Py_Pz_S_ab;
  Double I_ERI_F3x_Dyz_Pz_S_ab = I_ERI_G3xz_Py_Pz_S_ab+ABZ*I_ERI_F3x_Py_Pz_S_ab;
  Double I_ERI_F2xy_Dyz_Pz_S_ab = I_ERI_G2xyz_Py_Pz_S_ab+ABZ*I_ERI_F2xy_Py_Pz_S_ab;
  Double I_ERI_F2xz_Dyz_Pz_S_ab = I_ERI_G2x2z_Py_Pz_S_ab+ABZ*I_ERI_F2xz_Py_Pz_S_ab;
  Double I_ERI_Fx2y_Dyz_Pz_S_ab = I_ERI_Gx2yz_Py_Pz_S_ab+ABZ*I_ERI_Fx2y_Py_Pz_S_ab;
  Double I_ERI_Fxyz_Dyz_Pz_S_ab = I_ERI_Gxy2z_Py_Pz_S_ab+ABZ*I_ERI_Fxyz_Py_Pz_S_ab;
  Double I_ERI_Fx2z_Dyz_Pz_S_ab = I_ERI_Gx3z_Py_Pz_S_ab+ABZ*I_ERI_Fx2z_Py_Pz_S_ab;
  Double I_ERI_F3y_Dyz_Pz_S_ab = I_ERI_G3yz_Py_Pz_S_ab+ABZ*I_ERI_F3y_Py_Pz_S_ab;
  Double I_ERI_F2yz_Dyz_Pz_S_ab = I_ERI_G2y2z_Py_Pz_S_ab+ABZ*I_ERI_F2yz_Py_Pz_S_ab;
  Double I_ERI_Fy2z_Dyz_Pz_S_ab = I_ERI_Gy3z_Py_Pz_S_ab+ABZ*I_ERI_Fy2z_Py_Pz_S_ab;
  Double I_ERI_F3z_Dyz_Pz_S_ab = I_ERI_G4z_Py_Pz_S_ab+ABZ*I_ERI_F3z_Py_Pz_S_ab;
  Double I_ERI_F3x_D2z_Pz_S_ab = I_ERI_G3xz_Pz_Pz_S_ab+ABZ*I_ERI_F3x_Pz_Pz_S_ab;
  Double I_ERI_F2xy_D2z_Pz_S_ab = I_ERI_G2xyz_Pz_Pz_S_ab+ABZ*I_ERI_F2xy_Pz_Pz_S_ab;
  Double I_ERI_F2xz_D2z_Pz_S_ab = I_ERI_G2x2z_Pz_Pz_S_ab+ABZ*I_ERI_F2xz_Pz_Pz_S_ab;
  Double I_ERI_Fx2y_D2z_Pz_S_ab = I_ERI_Gx2yz_Pz_Pz_S_ab+ABZ*I_ERI_Fx2y_Pz_Pz_S_ab;
  Double I_ERI_Fxyz_D2z_Pz_S_ab = I_ERI_Gxy2z_Pz_Pz_S_ab+ABZ*I_ERI_Fxyz_Pz_Pz_S_ab;
  Double I_ERI_Fx2z_D2z_Pz_S_ab = I_ERI_Gx3z_Pz_Pz_S_ab+ABZ*I_ERI_Fx2z_Pz_Pz_S_ab;
  Double I_ERI_F3y_D2z_Pz_S_ab = I_ERI_G3yz_Pz_Pz_S_ab+ABZ*I_ERI_F3y_Pz_Pz_S_ab;
  Double I_ERI_F2yz_D2z_Pz_S_ab = I_ERI_G2y2z_Pz_Pz_S_ab+ABZ*I_ERI_F2yz_Pz_Pz_S_ab;
  Double I_ERI_Fy2z_D2z_Pz_S_ab = I_ERI_Gy3z_Pz_Pz_S_ab+ABZ*I_ERI_Fy2z_Pz_Pz_S_ab;
  Double I_ERI_F3z_D2z_Pz_S_ab = I_ERI_G4z_Pz_Pz_S_ab+ABZ*I_ERI_F3z_Pz_Pz_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_S_ac
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_S_D_S_ac
   ************************************************************/
  Double I_ERI_F3x_Px_D2x_S_ac = I_ERI_G4x_S_D2x_S_ac+ABX*I_ERI_F3x_S_D2x_S_ac;
  Double I_ERI_F2xy_Px_D2x_S_ac = I_ERI_G3xy_S_D2x_S_ac+ABX*I_ERI_F2xy_S_D2x_S_ac;
  Double I_ERI_F2xz_Px_D2x_S_ac = I_ERI_G3xz_S_D2x_S_ac+ABX*I_ERI_F2xz_S_D2x_S_ac;
  Double I_ERI_Fx2y_Px_D2x_S_ac = I_ERI_G2x2y_S_D2x_S_ac+ABX*I_ERI_Fx2y_S_D2x_S_ac;
  Double I_ERI_Fxyz_Px_D2x_S_ac = I_ERI_G2xyz_S_D2x_S_ac+ABX*I_ERI_Fxyz_S_D2x_S_ac;
  Double I_ERI_Fx2z_Px_D2x_S_ac = I_ERI_G2x2z_S_D2x_S_ac+ABX*I_ERI_Fx2z_S_D2x_S_ac;
  Double I_ERI_F3y_Px_D2x_S_ac = I_ERI_Gx3y_S_D2x_S_ac+ABX*I_ERI_F3y_S_D2x_S_ac;
  Double I_ERI_F2yz_Px_D2x_S_ac = I_ERI_Gx2yz_S_D2x_S_ac+ABX*I_ERI_F2yz_S_D2x_S_ac;
  Double I_ERI_Fy2z_Px_D2x_S_ac = I_ERI_Gxy2z_S_D2x_S_ac+ABX*I_ERI_Fy2z_S_D2x_S_ac;
  Double I_ERI_F3z_Px_D2x_S_ac = I_ERI_Gx3z_S_D2x_S_ac+ABX*I_ERI_F3z_S_D2x_S_ac;
  Double I_ERI_F3x_Py_D2x_S_ac = I_ERI_G3xy_S_D2x_S_ac+ABY*I_ERI_F3x_S_D2x_S_ac;
  Double I_ERI_F2xy_Py_D2x_S_ac = I_ERI_G2x2y_S_D2x_S_ac+ABY*I_ERI_F2xy_S_D2x_S_ac;
  Double I_ERI_F2xz_Py_D2x_S_ac = I_ERI_G2xyz_S_D2x_S_ac+ABY*I_ERI_F2xz_S_D2x_S_ac;
  Double I_ERI_Fx2y_Py_D2x_S_ac = I_ERI_Gx3y_S_D2x_S_ac+ABY*I_ERI_Fx2y_S_D2x_S_ac;
  Double I_ERI_Fxyz_Py_D2x_S_ac = I_ERI_Gx2yz_S_D2x_S_ac+ABY*I_ERI_Fxyz_S_D2x_S_ac;
  Double I_ERI_Fx2z_Py_D2x_S_ac = I_ERI_Gxy2z_S_D2x_S_ac+ABY*I_ERI_Fx2z_S_D2x_S_ac;
  Double I_ERI_F3y_Py_D2x_S_ac = I_ERI_G4y_S_D2x_S_ac+ABY*I_ERI_F3y_S_D2x_S_ac;
  Double I_ERI_F2yz_Py_D2x_S_ac = I_ERI_G3yz_S_D2x_S_ac+ABY*I_ERI_F2yz_S_D2x_S_ac;
  Double I_ERI_Fy2z_Py_D2x_S_ac = I_ERI_G2y2z_S_D2x_S_ac+ABY*I_ERI_Fy2z_S_D2x_S_ac;
  Double I_ERI_F3z_Py_D2x_S_ac = I_ERI_Gy3z_S_D2x_S_ac+ABY*I_ERI_F3z_S_D2x_S_ac;
  Double I_ERI_F3x_Pz_D2x_S_ac = I_ERI_G3xz_S_D2x_S_ac+ABZ*I_ERI_F3x_S_D2x_S_ac;
  Double I_ERI_F2xy_Pz_D2x_S_ac = I_ERI_G2xyz_S_D2x_S_ac+ABZ*I_ERI_F2xy_S_D2x_S_ac;
  Double I_ERI_F2xz_Pz_D2x_S_ac = I_ERI_G2x2z_S_D2x_S_ac+ABZ*I_ERI_F2xz_S_D2x_S_ac;
  Double I_ERI_Fx2y_Pz_D2x_S_ac = I_ERI_Gx2yz_S_D2x_S_ac+ABZ*I_ERI_Fx2y_S_D2x_S_ac;
  Double I_ERI_Fxyz_Pz_D2x_S_ac = I_ERI_Gxy2z_S_D2x_S_ac+ABZ*I_ERI_Fxyz_S_D2x_S_ac;
  Double I_ERI_Fx2z_Pz_D2x_S_ac = I_ERI_Gx3z_S_D2x_S_ac+ABZ*I_ERI_Fx2z_S_D2x_S_ac;
  Double I_ERI_F3y_Pz_D2x_S_ac = I_ERI_G3yz_S_D2x_S_ac+ABZ*I_ERI_F3y_S_D2x_S_ac;
  Double I_ERI_F2yz_Pz_D2x_S_ac = I_ERI_G2y2z_S_D2x_S_ac+ABZ*I_ERI_F2yz_S_D2x_S_ac;
  Double I_ERI_Fy2z_Pz_D2x_S_ac = I_ERI_Gy3z_S_D2x_S_ac+ABZ*I_ERI_Fy2z_S_D2x_S_ac;
  Double I_ERI_F3z_Pz_D2x_S_ac = I_ERI_G4z_S_D2x_S_ac+ABZ*I_ERI_F3z_S_D2x_S_ac;
  Double I_ERI_F3x_Px_Dxy_S_ac = I_ERI_G4x_S_Dxy_S_ac+ABX*I_ERI_F3x_S_Dxy_S_ac;
  Double I_ERI_F2xy_Px_Dxy_S_ac = I_ERI_G3xy_S_Dxy_S_ac+ABX*I_ERI_F2xy_S_Dxy_S_ac;
  Double I_ERI_F2xz_Px_Dxy_S_ac = I_ERI_G3xz_S_Dxy_S_ac+ABX*I_ERI_F2xz_S_Dxy_S_ac;
  Double I_ERI_Fx2y_Px_Dxy_S_ac = I_ERI_G2x2y_S_Dxy_S_ac+ABX*I_ERI_Fx2y_S_Dxy_S_ac;
  Double I_ERI_Fxyz_Px_Dxy_S_ac = I_ERI_G2xyz_S_Dxy_S_ac+ABX*I_ERI_Fxyz_S_Dxy_S_ac;
  Double I_ERI_Fx2z_Px_Dxy_S_ac = I_ERI_G2x2z_S_Dxy_S_ac+ABX*I_ERI_Fx2z_S_Dxy_S_ac;
  Double I_ERI_F3y_Px_Dxy_S_ac = I_ERI_Gx3y_S_Dxy_S_ac+ABX*I_ERI_F3y_S_Dxy_S_ac;
  Double I_ERI_F2yz_Px_Dxy_S_ac = I_ERI_Gx2yz_S_Dxy_S_ac+ABX*I_ERI_F2yz_S_Dxy_S_ac;
  Double I_ERI_Fy2z_Px_Dxy_S_ac = I_ERI_Gxy2z_S_Dxy_S_ac+ABX*I_ERI_Fy2z_S_Dxy_S_ac;
  Double I_ERI_F3z_Px_Dxy_S_ac = I_ERI_Gx3z_S_Dxy_S_ac+ABX*I_ERI_F3z_S_Dxy_S_ac;
  Double I_ERI_F3x_Py_Dxy_S_ac = I_ERI_G3xy_S_Dxy_S_ac+ABY*I_ERI_F3x_S_Dxy_S_ac;
  Double I_ERI_F2xy_Py_Dxy_S_ac = I_ERI_G2x2y_S_Dxy_S_ac+ABY*I_ERI_F2xy_S_Dxy_S_ac;
  Double I_ERI_F2xz_Py_Dxy_S_ac = I_ERI_G2xyz_S_Dxy_S_ac+ABY*I_ERI_F2xz_S_Dxy_S_ac;
  Double I_ERI_Fx2y_Py_Dxy_S_ac = I_ERI_Gx3y_S_Dxy_S_ac+ABY*I_ERI_Fx2y_S_Dxy_S_ac;
  Double I_ERI_Fxyz_Py_Dxy_S_ac = I_ERI_Gx2yz_S_Dxy_S_ac+ABY*I_ERI_Fxyz_S_Dxy_S_ac;
  Double I_ERI_Fx2z_Py_Dxy_S_ac = I_ERI_Gxy2z_S_Dxy_S_ac+ABY*I_ERI_Fx2z_S_Dxy_S_ac;
  Double I_ERI_F3y_Py_Dxy_S_ac = I_ERI_G4y_S_Dxy_S_ac+ABY*I_ERI_F3y_S_Dxy_S_ac;
  Double I_ERI_F2yz_Py_Dxy_S_ac = I_ERI_G3yz_S_Dxy_S_ac+ABY*I_ERI_F2yz_S_Dxy_S_ac;
  Double I_ERI_Fy2z_Py_Dxy_S_ac = I_ERI_G2y2z_S_Dxy_S_ac+ABY*I_ERI_Fy2z_S_Dxy_S_ac;
  Double I_ERI_F3z_Py_Dxy_S_ac = I_ERI_Gy3z_S_Dxy_S_ac+ABY*I_ERI_F3z_S_Dxy_S_ac;
  Double I_ERI_F3x_Pz_Dxy_S_ac = I_ERI_G3xz_S_Dxy_S_ac+ABZ*I_ERI_F3x_S_Dxy_S_ac;
  Double I_ERI_F2xy_Pz_Dxy_S_ac = I_ERI_G2xyz_S_Dxy_S_ac+ABZ*I_ERI_F2xy_S_Dxy_S_ac;
  Double I_ERI_F2xz_Pz_Dxy_S_ac = I_ERI_G2x2z_S_Dxy_S_ac+ABZ*I_ERI_F2xz_S_Dxy_S_ac;
  Double I_ERI_Fx2y_Pz_Dxy_S_ac = I_ERI_Gx2yz_S_Dxy_S_ac+ABZ*I_ERI_Fx2y_S_Dxy_S_ac;
  Double I_ERI_Fxyz_Pz_Dxy_S_ac = I_ERI_Gxy2z_S_Dxy_S_ac+ABZ*I_ERI_Fxyz_S_Dxy_S_ac;
  Double I_ERI_Fx2z_Pz_Dxy_S_ac = I_ERI_Gx3z_S_Dxy_S_ac+ABZ*I_ERI_Fx2z_S_Dxy_S_ac;
  Double I_ERI_F3y_Pz_Dxy_S_ac = I_ERI_G3yz_S_Dxy_S_ac+ABZ*I_ERI_F3y_S_Dxy_S_ac;
  Double I_ERI_F2yz_Pz_Dxy_S_ac = I_ERI_G2y2z_S_Dxy_S_ac+ABZ*I_ERI_F2yz_S_Dxy_S_ac;
  Double I_ERI_Fy2z_Pz_Dxy_S_ac = I_ERI_Gy3z_S_Dxy_S_ac+ABZ*I_ERI_Fy2z_S_Dxy_S_ac;
  Double I_ERI_F3z_Pz_Dxy_S_ac = I_ERI_G4z_S_Dxy_S_ac+ABZ*I_ERI_F3z_S_Dxy_S_ac;
  Double I_ERI_F3x_Px_Dxz_S_ac = I_ERI_G4x_S_Dxz_S_ac+ABX*I_ERI_F3x_S_Dxz_S_ac;
  Double I_ERI_F2xy_Px_Dxz_S_ac = I_ERI_G3xy_S_Dxz_S_ac+ABX*I_ERI_F2xy_S_Dxz_S_ac;
  Double I_ERI_F2xz_Px_Dxz_S_ac = I_ERI_G3xz_S_Dxz_S_ac+ABX*I_ERI_F2xz_S_Dxz_S_ac;
  Double I_ERI_Fx2y_Px_Dxz_S_ac = I_ERI_G2x2y_S_Dxz_S_ac+ABX*I_ERI_Fx2y_S_Dxz_S_ac;
  Double I_ERI_Fxyz_Px_Dxz_S_ac = I_ERI_G2xyz_S_Dxz_S_ac+ABX*I_ERI_Fxyz_S_Dxz_S_ac;
  Double I_ERI_Fx2z_Px_Dxz_S_ac = I_ERI_G2x2z_S_Dxz_S_ac+ABX*I_ERI_Fx2z_S_Dxz_S_ac;
  Double I_ERI_F3y_Px_Dxz_S_ac = I_ERI_Gx3y_S_Dxz_S_ac+ABX*I_ERI_F3y_S_Dxz_S_ac;
  Double I_ERI_F2yz_Px_Dxz_S_ac = I_ERI_Gx2yz_S_Dxz_S_ac+ABX*I_ERI_F2yz_S_Dxz_S_ac;
  Double I_ERI_Fy2z_Px_Dxz_S_ac = I_ERI_Gxy2z_S_Dxz_S_ac+ABX*I_ERI_Fy2z_S_Dxz_S_ac;
  Double I_ERI_F3z_Px_Dxz_S_ac = I_ERI_Gx3z_S_Dxz_S_ac+ABX*I_ERI_F3z_S_Dxz_S_ac;
  Double I_ERI_F3x_Py_Dxz_S_ac = I_ERI_G3xy_S_Dxz_S_ac+ABY*I_ERI_F3x_S_Dxz_S_ac;
  Double I_ERI_F2xy_Py_Dxz_S_ac = I_ERI_G2x2y_S_Dxz_S_ac+ABY*I_ERI_F2xy_S_Dxz_S_ac;
  Double I_ERI_F2xz_Py_Dxz_S_ac = I_ERI_G2xyz_S_Dxz_S_ac+ABY*I_ERI_F2xz_S_Dxz_S_ac;
  Double I_ERI_Fx2y_Py_Dxz_S_ac = I_ERI_Gx3y_S_Dxz_S_ac+ABY*I_ERI_Fx2y_S_Dxz_S_ac;
  Double I_ERI_Fxyz_Py_Dxz_S_ac = I_ERI_Gx2yz_S_Dxz_S_ac+ABY*I_ERI_Fxyz_S_Dxz_S_ac;
  Double I_ERI_Fx2z_Py_Dxz_S_ac = I_ERI_Gxy2z_S_Dxz_S_ac+ABY*I_ERI_Fx2z_S_Dxz_S_ac;
  Double I_ERI_F3y_Py_Dxz_S_ac = I_ERI_G4y_S_Dxz_S_ac+ABY*I_ERI_F3y_S_Dxz_S_ac;
  Double I_ERI_F2yz_Py_Dxz_S_ac = I_ERI_G3yz_S_Dxz_S_ac+ABY*I_ERI_F2yz_S_Dxz_S_ac;
  Double I_ERI_Fy2z_Py_Dxz_S_ac = I_ERI_G2y2z_S_Dxz_S_ac+ABY*I_ERI_Fy2z_S_Dxz_S_ac;
  Double I_ERI_F3z_Py_Dxz_S_ac = I_ERI_Gy3z_S_Dxz_S_ac+ABY*I_ERI_F3z_S_Dxz_S_ac;
  Double I_ERI_F3x_Pz_Dxz_S_ac = I_ERI_G3xz_S_Dxz_S_ac+ABZ*I_ERI_F3x_S_Dxz_S_ac;
  Double I_ERI_F2xy_Pz_Dxz_S_ac = I_ERI_G2xyz_S_Dxz_S_ac+ABZ*I_ERI_F2xy_S_Dxz_S_ac;
  Double I_ERI_F2xz_Pz_Dxz_S_ac = I_ERI_G2x2z_S_Dxz_S_ac+ABZ*I_ERI_F2xz_S_Dxz_S_ac;
  Double I_ERI_Fx2y_Pz_Dxz_S_ac = I_ERI_Gx2yz_S_Dxz_S_ac+ABZ*I_ERI_Fx2y_S_Dxz_S_ac;
  Double I_ERI_Fxyz_Pz_Dxz_S_ac = I_ERI_Gxy2z_S_Dxz_S_ac+ABZ*I_ERI_Fxyz_S_Dxz_S_ac;
  Double I_ERI_Fx2z_Pz_Dxz_S_ac = I_ERI_Gx3z_S_Dxz_S_ac+ABZ*I_ERI_Fx2z_S_Dxz_S_ac;
  Double I_ERI_F3y_Pz_Dxz_S_ac = I_ERI_G3yz_S_Dxz_S_ac+ABZ*I_ERI_F3y_S_Dxz_S_ac;
  Double I_ERI_F2yz_Pz_Dxz_S_ac = I_ERI_G2y2z_S_Dxz_S_ac+ABZ*I_ERI_F2yz_S_Dxz_S_ac;
  Double I_ERI_Fy2z_Pz_Dxz_S_ac = I_ERI_Gy3z_S_Dxz_S_ac+ABZ*I_ERI_Fy2z_S_Dxz_S_ac;
  Double I_ERI_F3z_Pz_Dxz_S_ac = I_ERI_G4z_S_Dxz_S_ac+ABZ*I_ERI_F3z_S_Dxz_S_ac;
  Double I_ERI_F3x_Px_D2y_S_ac = I_ERI_G4x_S_D2y_S_ac+ABX*I_ERI_F3x_S_D2y_S_ac;
  Double I_ERI_F2xy_Px_D2y_S_ac = I_ERI_G3xy_S_D2y_S_ac+ABX*I_ERI_F2xy_S_D2y_S_ac;
  Double I_ERI_F2xz_Px_D2y_S_ac = I_ERI_G3xz_S_D2y_S_ac+ABX*I_ERI_F2xz_S_D2y_S_ac;
  Double I_ERI_Fx2y_Px_D2y_S_ac = I_ERI_G2x2y_S_D2y_S_ac+ABX*I_ERI_Fx2y_S_D2y_S_ac;
  Double I_ERI_Fxyz_Px_D2y_S_ac = I_ERI_G2xyz_S_D2y_S_ac+ABX*I_ERI_Fxyz_S_D2y_S_ac;
  Double I_ERI_Fx2z_Px_D2y_S_ac = I_ERI_G2x2z_S_D2y_S_ac+ABX*I_ERI_Fx2z_S_D2y_S_ac;
  Double I_ERI_F3y_Px_D2y_S_ac = I_ERI_Gx3y_S_D2y_S_ac+ABX*I_ERI_F3y_S_D2y_S_ac;
  Double I_ERI_F2yz_Px_D2y_S_ac = I_ERI_Gx2yz_S_D2y_S_ac+ABX*I_ERI_F2yz_S_D2y_S_ac;
  Double I_ERI_Fy2z_Px_D2y_S_ac = I_ERI_Gxy2z_S_D2y_S_ac+ABX*I_ERI_Fy2z_S_D2y_S_ac;
  Double I_ERI_F3z_Px_D2y_S_ac = I_ERI_Gx3z_S_D2y_S_ac+ABX*I_ERI_F3z_S_D2y_S_ac;
  Double I_ERI_F3x_Py_D2y_S_ac = I_ERI_G3xy_S_D2y_S_ac+ABY*I_ERI_F3x_S_D2y_S_ac;
  Double I_ERI_F2xy_Py_D2y_S_ac = I_ERI_G2x2y_S_D2y_S_ac+ABY*I_ERI_F2xy_S_D2y_S_ac;
  Double I_ERI_F2xz_Py_D2y_S_ac = I_ERI_G2xyz_S_D2y_S_ac+ABY*I_ERI_F2xz_S_D2y_S_ac;
  Double I_ERI_Fx2y_Py_D2y_S_ac = I_ERI_Gx3y_S_D2y_S_ac+ABY*I_ERI_Fx2y_S_D2y_S_ac;
  Double I_ERI_Fxyz_Py_D2y_S_ac = I_ERI_Gx2yz_S_D2y_S_ac+ABY*I_ERI_Fxyz_S_D2y_S_ac;
  Double I_ERI_Fx2z_Py_D2y_S_ac = I_ERI_Gxy2z_S_D2y_S_ac+ABY*I_ERI_Fx2z_S_D2y_S_ac;
  Double I_ERI_F3y_Py_D2y_S_ac = I_ERI_G4y_S_D2y_S_ac+ABY*I_ERI_F3y_S_D2y_S_ac;
  Double I_ERI_F2yz_Py_D2y_S_ac = I_ERI_G3yz_S_D2y_S_ac+ABY*I_ERI_F2yz_S_D2y_S_ac;
  Double I_ERI_Fy2z_Py_D2y_S_ac = I_ERI_G2y2z_S_D2y_S_ac+ABY*I_ERI_Fy2z_S_D2y_S_ac;
  Double I_ERI_F3z_Py_D2y_S_ac = I_ERI_Gy3z_S_D2y_S_ac+ABY*I_ERI_F3z_S_D2y_S_ac;
  Double I_ERI_F3x_Pz_D2y_S_ac = I_ERI_G3xz_S_D2y_S_ac+ABZ*I_ERI_F3x_S_D2y_S_ac;
  Double I_ERI_F2xy_Pz_D2y_S_ac = I_ERI_G2xyz_S_D2y_S_ac+ABZ*I_ERI_F2xy_S_D2y_S_ac;
  Double I_ERI_F2xz_Pz_D2y_S_ac = I_ERI_G2x2z_S_D2y_S_ac+ABZ*I_ERI_F2xz_S_D2y_S_ac;
  Double I_ERI_Fx2y_Pz_D2y_S_ac = I_ERI_Gx2yz_S_D2y_S_ac+ABZ*I_ERI_Fx2y_S_D2y_S_ac;
  Double I_ERI_Fxyz_Pz_D2y_S_ac = I_ERI_Gxy2z_S_D2y_S_ac+ABZ*I_ERI_Fxyz_S_D2y_S_ac;
  Double I_ERI_Fx2z_Pz_D2y_S_ac = I_ERI_Gx3z_S_D2y_S_ac+ABZ*I_ERI_Fx2z_S_D2y_S_ac;
  Double I_ERI_F3y_Pz_D2y_S_ac = I_ERI_G3yz_S_D2y_S_ac+ABZ*I_ERI_F3y_S_D2y_S_ac;
  Double I_ERI_F2yz_Pz_D2y_S_ac = I_ERI_G2y2z_S_D2y_S_ac+ABZ*I_ERI_F2yz_S_D2y_S_ac;
  Double I_ERI_Fy2z_Pz_D2y_S_ac = I_ERI_Gy3z_S_D2y_S_ac+ABZ*I_ERI_Fy2z_S_D2y_S_ac;
  Double I_ERI_F3z_Pz_D2y_S_ac = I_ERI_G4z_S_D2y_S_ac+ABZ*I_ERI_F3z_S_D2y_S_ac;
  Double I_ERI_F3x_Px_Dyz_S_ac = I_ERI_G4x_S_Dyz_S_ac+ABX*I_ERI_F3x_S_Dyz_S_ac;
  Double I_ERI_F2xy_Px_Dyz_S_ac = I_ERI_G3xy_S_Dyz_S_ac+ABX*I_ERI_F2xy_S_Dyz_S_ac;
  Double I_ERI_F2xz_Px_Dyz_S_ac = I_ERI_G3xz_S_Dyz_S_ac+ABX*I_ERI_F2xz_S_Dyz_S_ac;
  Double I_ERI_Fx2y_Px_Dyz_S_ac = I_ERI_G2x2y_S_Dyz_S_ac+ABX*I_ERI_Fx2y_S_Dyz_S_ac;
  Double I_ERI_Fxyz_Px_Dyz_S_ac = I_ERI_G2xyz_S_Dyz_S_ac+ABX*I_ERI_Fxyz_S_Dyz_S_ac;
  Double I_ERI_Fx2z_Px_Dyz_S_ac = I_ERI_G2x2z_S_Dyz_S_ac+ABX*I_ERI_Fx2z_S_Dyz_S_ac;
  Double I_ERI_F3y_Px_Dyz_S_ac = I_ERI_Gx3y_S_Dyz_S_ac+ABX*I_ERI_F3y_S_Dyz_S_ac;
  Double I_ERI_F2yz_Px_Dyz_S_ac = I_ERI_Gx2yz_S_Dyz_S_ac+ABX*I_ERI_F2yz_S_Dyz_S_ac;
  Double I_ERI_Fy2z_Px_Dyz_S_ac = I_ERI_Gxy2z_S_Dyz_S_ac+ABX*I_ERI_Fy2z_S_Dyz_S_ac;
  Double I_ERI_F3z_Px_Dyz_S_ac = I_ERI_Gx3z_S_Dyz_S_ac+ABX*I_ERI_F3z_S_Dyz_S_ac;
  Double I_ERI_F3x_Py_Dyz_S_ac = I_ERI_G3xy_S_Dyz_S_ac+ABY*I_ERI_F3x_S_Dyz_S_ac;
  Double I_ERI_F2xy_Py_Dyz_S_ac = I_ERI_G2x2y_S_Dyz_S_ac+ABY*I_ERI_F2xy_S_Dyz_S_ac;
  Double I_ERI_F2xz_Py_Dyz_S_ac = I_ERI_G2xyz_S_Dyz_S_ac+ABY*I_ERI_F2xz_S_Dyz_S_ac;
  Double I_ERI_Fx2y_Py_Dyz_S_ac = I_ERI_Gx3y_S_Dyz_S_ac+ABY*I_ERI_Fx2y_S_Dyz_S_ac;
  Double I_ERI_Fxyz_Py_Dyz_S_ac = I_ERI_Gx2yz_S_Dyz_S_ac+ABY*I_ERI_Fxyz_S_Dyz_S_ac;
  Double I_ERI_Fx2z_Py_Dyz_S_ac = I_ERI_Gxy2z_S_Dyz_S_ac+ABY*I_ERI_Fx2z_S_Dyz_S_ac;
  Double I_ERI_F3y_Py_Dyz_S_ac = I_ERI_G4y_S_Dyz_S_ac+ABY*I_ERI_F3y_S_Dyz_S_ac;
  Double I_ERI_F2yz_Py_Dyz_S_ac = I_ERI_G3yz_S_Dyz_S_ac+ABY*I_ERI_F2yz_S_Dyz_S_ac;
  Double I_ERI_Fy2z_Py_Dyz_S_ac = I_ERI_G2y2z_S_Dyz_S_ac+ABY*I_ERI_Fy2z_S_Dyz_S_ac;
  Double I_ERI_F3z_Py_Dyz_S_ac = I_ERI_Gy3z_S_Dyz_S_ac+ABY*I_ERI_F3z_S_Dyz_S_ac;
  Double I_ERI_F3x_Pz_Dyz_S_ac = I_ERI_G3xz_S_Dyz_S_ac+ABZ*I_ERI_F3x_S_Dyz_S_ac;
  Double I_ERI_F2xy_Pz_Dyz_S_ac = I_ERI_G2xyz_S_Dyz_S_ac+ABZ*I_ERI_F2xy_S_Dyz_S_ac;
  Double I_ERI_F2xz_Pz_Dyz_S_ac = I_ERI_G2x2z_S_Dyz_S_ac+ABZ*I_ERI_F2xz_S_Dyz_S_ac;
  Double I_ERI_Fx2y_Pz_Dyz_S_ac = I_ERI_Gx2yz_S_Dyz_S_ac+ABZ*I_ERI_Fx2y_S_Dyz_S_ac;
  Double I_ERI_Fxyz_Pz_Dyz_S_ac = I_ERI_Gxy2z_S_Dyz_S_ac+ABZ*I_ERI_Fxyz_S_Dyz_S_ac;
  Double I_ERI_Fx2z_Pz_Dyz_S_ac = I_ERI_Gx3z_S_Dyz_S_ac+ABZ*I_ERI_Fx2z_S_Dyz_S_ac;
  Double I_ERI_F3y_Pz_Dyz_S_ac = I_ERI_G3yz_S_Dyz_S_ac+ABZ*I_ERI_F3y_S_Dyz_S_ac;
  Double I_ERI_F2yz_Pz_Dyz_S_ac = I_ERI_G2y2z_S_Dyz_S_ac+ABZ*I_ERI_F2yz_S_Dyz_S_ac;
  Double I_ERI_Fy2z_Pz_Dyz_S_ac = I_ERI_Gy3z_S_Dyz_S_ac+ABZ*I_ERI_Fy2z_S_Dyz_S_ac;
  Double I_ERI_F3z_Pz_Dyz_S_ac = I_ERI_G4z_S_Dyz_S_ac+ABZ*I_ERI_F3z_S_Dyz_S_ac;
  Double I_ERI_F3x_Px_D2z_S_ac = I_ERI_G4x_S_D2z_S_ac+ABX*I_ERI_F3x_S_D2z_S_ac;
  Double I_ERI_F2xy_Px_D2z_S_ac = I_ERI_G3xy_S_D2z_S_ac+ABX*I_ERI_F2xy_S_D2z_S_ac;
  Double I_ERI_F2xz_Px_D2z_S_ac = I_ERI_G3xz_S_D2z_S_ac+ABX*I_ERI_F2xz_S_D2z_S_ac;
  Double I_ERI_Fx2y_Px_D2z_S_ac = I_ERI_G2x2y_S_D2z_S_ac+ABX*I_ERI_Fx2y_S_D2z_S_ac;
  Double I_ERI_Fxyz_Px_D2z_S_ac = I_ERI_G2xyz_S_D2z_S_ac+ABX*I_ERI_Fxyz_S_D2z_S_ac;
  Double I_ERI_Fx2z_Px_D2z_S_ac = I_ERI_G2x2z_S_D2z_S_ac+ABX*I_ERI_Fx2z_S_D2z_S_ac;
  Double I_ERI_F3y_Px_D2z_S_ac = I_ERI_Gx3y_S_D2z_S_ac+ABX*I_ERI_F3y_S_D2z_S_ac;
  Double I_ERI_F2yz_Px_D2z_S_ac = I_ERI_Gx2yz_S_D2z_S_ac+ABX*I_ERI_F2yz_S_D2z_S_ac;
  Double I_ERI_Fy2z_Px_D2z_S_ac = I_ERI_Gxy2z_S_D2z_S_ac+ABX*I_ERI_Fy2z_S_D2z_S_ac;
  Double I_ERI_F3z_Px_D2z_S_ac = I_ERI_Gx3z_S_D2z_S_ac+ABX*I_ERI_F3z_S_D2z_S_ac;
  Double I_ERI_F3x_Py_D2z_S_ac = I_ERI_G3xy_S_D2z_S_ac+ABY*I_ERI_F3x_S_D2z_S_ac;
  Double I_ERI_F2xy_Py_D2z_S_ac = I_ERI_G2x2y_S_D2z_S_ac+ABY*I_ERI_F2xy_S_D2z_S_ac;
  Double I_ERI_F2xz_Py_D2z_S_ac = I_ERI_G2xyz_S_D2z_S_ac+ABY*I_ERI_F2xz_S_D2z_S_ac;
  Double I_ERI_Fx2y_Py_D2z_S_ac = I_ERI_Gx3y_S_D2z_S_ac+ABY*I_ERI_Fx2y_S_D2z_S_ac;
  Double I_ERI_Fxyz_Py_D2z_S_ac = I_ERI_Gx2yz_S_D2z_S_ac+ABY*I_ERI_Fxyz_S_D2z_S_ac;
  Double I_ERI_Fx2z_Py_D2z_S_ac = I_ERI_Gxy2z_S_D2z_S_ac+ABY*I_ERI_Fx2z_S_D2z_S_ac;
  Double I_ERI_F3y_Py_D2z_S_ac = I_ERI_G4y_S_D2z_S_ac+ABY*I_ERI_F3y_S_D2z_S_ac;
  Double I_ERI_F2yz_Py_D2z_S_ac = I_ERI_G3yz_S_D2z_S_ac+ABY*I_ERI_F2yz_S_D2z_S_ac;
  Double I_ERI_Fy2z_Py_D2z_S_ac = I_ERI_G2y2z_S_D2z_S_ac+ABY*I_ERI_Fy2z_S_D2z_S_ac;
  Double I_ERI_F3z_Py_D2z_S_ac = I_ERI_Gy3z_S_D2z_S_ac+ABY*I_ERI_F3z_S_D2z_S_ac;
  Double I_ERI_F3x_Pz_D2z_S_ac = I_ERI_G3xz_S_D2z_S_ac+ABZ*I_ERI_F3x_S_D2z_S_ac;
  Double I_ERI_F2xy_Pz_D2z_S_ac = I_ERI_G2xyz_S_D2z_S_ac+ABZ*I_ERI_F2xy_S_D2z_S_ac;
  Double I_ERI_F2xz_Pz_D2z_S_ac = I_ERI_G2x2z_S_D2z_S_ac+ABZ*I_ERI_F2xz_S_D2z_S_ac;
  Double I_ERI_Fx2y_Pz_D2z_S_ac = I_ERI_Gx2yz_S_D2z_S_ac+ABZ*I_ERI_Fx2y_S_D2z_S_ac;
  Double I_ERI_Fxyz_Pz_D2z_S_ac = I_ERI_Gxy2z_S_D2z_S_ac+ABZ*I_ERI_Fxyz_S_D2z_S_ac;
  Double I_ERI_Fx2z_Pz_D2z_S_ac = I_ERI_Gx3z_S_D2z_S_ac+ABZ*I_ERI_Fx2z_S_D2z_S_ac;
  Double I_ERI_F3y_Pz_D2z_S_ac = I_ERI_G3yz_S_D2z_S_ac+ABZ*I_ERI_F3y_S_D2z_S_ac;
  Double I_ERI_F2yz_Pz_D2z_S_ac = I_ERI_G2y2z_S_D2z_S_ac+ABZ*I_ERI_F2yz_S_D2z_S_ac;
  Double I_ERI_Fy2z_Pz_D2z_S_ac = I_ERI_Gy3z_S_D2z_S_ac+ABZ*I_ERI_Fy2z_S_D2z_S_ac;
  Double I_ERI_F3z_Pz_D2z_S_ac = I_ERI_G4z_S_D2z_S_ac+ABZ*I_ERI_F3z_S_D2z_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_bb
   * RHS shell quartet name: SQ_ERI_D_S_P_S_bb
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_bb = I_ERI_F3x_S_Px_S_bb+ABX*I_ERI_D2x_S_Px_S_bb;
  Double I_ERI_Dxy_Px_Px_S_bb = I_ERI_F2xy_S_Px_S_bb+ABX*I_ERI_Dxy_S_Px_S_bb;
  Double I_ERI_Dxz_Px_Px_S_bb = I_ERI_F2xz_S_Px_S_bb+ABX*I_ERI_Dxz_S_Px_S_bb;
  Double I_ERI_D2y_Px_Px_S_bb = I_ERI_Fx2y_S_Px_S_bb+ABX*I_ERI_D2y_S_Px_S_bb;
  Double I_ERI_Dyz_Px_Px_S_bb = I_ERI_Fxyz_S_Px_S_bb+ABX*I_ERI_Dyz_S_Px_S_bb;
  Double I_ERI_D2z_Px_Px_S_bb = I_ERI_Fx2z_S_Px_S_bb+ABX*I_ERI_D2z_S_Px_S_bb;
  Double I_ERI_D2x_Py_Px_S_bb = I_ERI_F2xy_S_Px_S_bb+ABY*I_ERI_D2x_S_Px_S_bb;
  Double I_ERI_Dxy_Py_Px_S_bb = I_ERI_Fx2y_S_Px_S_bb+ABY*I_ERI_Dxy_S_Px_S_bb;
  Double I_ERI_Dxz_Py_Px_S_bb = I_ERI_Fxyz_S_Px_S_bb+ABY*I_ERI_Dxz_S_Px_S_bb;
  Double I_ERI_D2y_Py_Px_S_bb = I_ERI_F3y_S_Px_S_bb+ABY*I_ERI_D2y_S_Px_S_bb;
  Double I_ERI_Dyz_Py_Px_S_bb = I_ERI_F2yz_S_Px_S_bb+ABY*I_ERI_Dyz_S_Px_S_bb;
  Double I_ERI_D2z_Py_Px_S_bb = I_ERI_Fy2z_S_Px_S_bb+ABY*I_ERI_D2z_S_Px_S_bb;
  Double I_ERI_D2x_Pz_Px_S_bb = I_ERI_F2xz_S_Px_S_bb+ABZ*I_ERI_D2x_S_Px_S_bb;
  Double I_ERI_Dxy_Pz_Px_S_bb = I_ERI_Fxyz_S_Px_S_bb+ABZ*I_ERI_Dxy_S_Px_S_bb;
  Double I_ERI_Dxz_Pz_Px_S_bb = I_ERI_Fx2z_S_Px_S_bb+ABZ*I_ERI_Dxz_S_Px_S_bb;
  Double I_ERI_D2y_Pz_Px_S_bb = I_ERI_F2yz_S_Px_S_bb+ABZ*I_ERI_D2y_S_Px_S_bb;
  Double I_ERI_Dyz_Pz_Px_S_bb = I_ERI_Fy2z_S_Px_S_bb+ABZ*I_ERI_Dyz_S_Px_S_bb;
  Double I_ERI_D2z_Pz_Px_S_bb = I_ERI_F3z_S_Px_S_bb+ABZ*I_ERI_D2z_S_Px_S_bb;
  Double I_ERI_D2x_Px_Py_S_bb = I_ERI_F3x_S_Py_S_bb+ABX*I_ERI_D2x_S_Py_S_bb;
  Double I_ERI_Dxy_Px_Py_S_bb = I_ERI_F2xy_S_Py_S_bb+ABX*I_ERI_Dxy_S_Py_S_bb;
  Double I_ERI_Dxz_Px_Py_S_bb = I_ERI_F2xz_S_Py_S_bb+ABX*I_ERI_Dxz_S_Py_S_bb;
  Double I_ERI_D2y_Px_Py_S_bb = I_ERI_Fx2y_S_Py_S_bb+ABX*I_ERI_D2y_S_Py_S_bb;
  Double I_ERI_Dyz_Px_Py_S_bb = I_ERI_Fxyz_S_Py_S_bb+ABX*I_ERI_Dyz_S_Py_S_bb;
  Double I_ERI_D2z_Px_Py_S_bb = I_ERI_Fx2z_S_Py_S_bb+ABX*I_ERI_D2z_S_Py_S_bb;
  Double I_ERI_D2x_Py_Py_S_bb = I_ERI_F2xy_S_Py_S_bb+ABY*I_ERI_D2x_S_Py_S_bb;
  Double I_ERI_Dxy_Py_Py_S_bb = I_ERI_Fx2y_S_Py_S_bb+ABY*I_ERI_Dxy_S_Py_S_bb;
  Double I_ERI_Dxz_Py_Py_S_bb = I_ERI_Fxyz_S_Py_S_bb+ABY*I_ERI_Dxz_S_Py_S_bb;
  Double I_ERI_D2y_Py_Py_S_bb = I_ERI_F3y_S_Py_S_bb+ABY*I_ERI_D2y_S_Py_S_bb;
  Double I_ERI_Dyz_Py_Py_S_bb = I_ERI_F2yz_S_Py_S_bb+ABY*I_ERI_Dyz_S_Py_S_bb;
  Double I_ERI_D2z_Py_Py_S_bb = I_ERI_Fy2z_S_Py_S_bb+ABY*I_ERI_D2z_S_Py_S_bb;
  Double I_ERI_D2x_Pz_Py_S_bb = I_ERI_F2xz_S_Py_S_bb+ABZ*I_ERI_D2x_S_Py_S_bb;
  Double I_ERI_Dxy_Pz_Py_S_bb = I_ERI_Fxyz_S_Py_S_bb+ABZ*I_ERI_Dxy_S_Py_S_bb;
  Double I_ERI_Dxz_Pz_Py_S_bb = I_ERI_Fx2z_S_Py_S_bb+ABZ*I_ERI_Dxz_S_Py_S_bb;
  Double I_ERI_D2y_Pz_Py_S_bb = I_ERI_F2yz_S_Py_S_bb+ABZ*I_ERI_D2y_S_Py_S_bb;
  Double I_ERI_Dyz_Pz_Py_S_bb = I_ERI_Fy2z_S_Py_S_bb+ABZ*I_ERI_Dyz_S_Py_S_bb;
  Double I_ERI_D2z_Pz_Py_S_bb = I_ERI_F3z_S_Py_S_bb+ABZ*I_ERI_D2z_S_Py_S_bb;
  Double I_ERI_D2x_Px_Pz_S_bb = I_ERI_F3x_S_Pz_S_bb+ABX*I_ERI_D2x_S_Pz_S_bb;
  Double I_ERI_Dxy_Px_Pz_S_bb = I_ERI_F2xy_S_Pz_S_bb+ABX*I_ERI_Dxy_S_Pz_S_bb;
  Double I_ERI_Dxz_Px_Pz_S_bb = I_ERI_F2xz_S_Pz_S_bb+ABX*I_ERI_Dxz_S_Pz_S_bb;
  Double I_ERI_D2y_Px_Pz_S_bb = I_ERI_Fx2y_S_Pz_S_bb+ABX*I_ERI_D2y_S_Pz_S_bb;
  Double I_ERI_Dyz_Px_Pz_S_bb = I_ERI_Fxyz_S_Pz_S_bb+ABX*I_ERI_Dyz_S_Pz_S_bb;
  Double I_ERI_D2z_Px_Pz_S_bb = I_ERI_Fx2z_S_Pz_S_bb+ABX*I_ERI_D2z_S_Pz_S_bb;
  Double I_ERI_D2x_Py_Pz_S_bb = I_ERI_F2xy_S_Pz_S_bb+ABY*I_ERI_D2x_S_Pz_S_bb;
  Double I_ERI_Dxy_Py_Pz_S_bb = I_ERI_Fx2y_S_Pz_S_bb+ABY*I_ERI_Dxy_S_Pz_S_bb;
  Double I_ERI_Dxz_Py_Pz_S_bb = I_ERI_Fxyz_S_Pz_S_bb+ABY*I_ERI_Dxz_S_Pz_S_bb;
  Double I_ERI_D2y_Py_Pz_S_bb = I_ERI_F3y_S_Pz_S_bb+ABY*I_ERI_D2y_S_Pz_S_bb;
  Double I_ERI_Dyz_Py_Pz_S_bb = I_ERI_F2yz_S_Pz_S_bb+ABY*I_ERI_Dyz_S_Pz_S_bb;
  Double I_ERI_D2z_Py_Pz_S_bb = I_ERI_Fy2z_S_Pz_S_bb+ABY*I_ERI_D2z_S_Pz_S_bb;
  Double I_ERI_D2x_Pz_Pz_S_bb = I_ERI_F2xz_S_Pz_S_bb+ABZ*I_ERI_D2x_S_Pz_S_bb;
  Double I_ERI_Dxy_Pz_Pz_S_bb = I_ERI_Fxyz_S_Pz_S_bb+ABZ*I_ERI_Dxy_S_Pz_S_bb;
  Double I_ERI_Dxz_Pz_Pz_S_bb = I_ERI_Fx2z_S_Pz_S_bb+ABZ*I_ERI_Dxz_S_Pz_S_bb;
  Double I_ERI_D2y_Pz_Pz_S_bb = I_ERI_F2yz_S_Pz_S_bb+ABZ*I_ERI_D2y_S_Pz_S_bb;
  Double I_ERI_Dyz_Pz_Pz_S_bb = I_ERI_Fy2z_S_Pz_S_bb+ABZ*I_ERI_Dyz_S_Pz_S_bb;
  Double I_ERI_D2z_Pz_Pz_S_bb = I_ERI_F3z_S_Pz_S_bb+ABZ*I_ERI_D2z_S_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_P_S_bb
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_bb = I_ERI_G4x_S_Px_S_bb+ABX*I_ERI_F3x_S_Px_S_bb;
  Double I_ERI_F2xy_Px_Px_S_bb = I_ERI_G3xy_S_Px_S_bb+ABX*I_ERI_F2xy_S_Px_S_bb;
  Double I_ERI_F2xz_Px_Px_S_bb = I_ERI_G3xz_S_Px_S_bb+ABX*I_ERI_F2xz_S_Px_S_bb;
  Double I_ERI_Fx2y_Px_Px_S_bb = I_ERI_G2x2y_S_Px_S_bb+ABX*I_ERI_Fx2y_S_Px_S_bb;
  Double I_ERI_Fxyz_Px_Px_S_bb = I_ERI_G2xyz_S_Px_S_bb+ABX*I_ERI_Fxyz_S_Px_S_bb;
  Double I_ERI_Fx2z_Px_Px_S_bb = I_ERI_G2x2z_S_Px_S_bb+ABX*I_ERI_Fx2z_S_Px_S_bb;
  Double I_ERI_F3y_Px_Px_S_bb = I_ERI_Gx3y_S_Px_S_bb+ABX*I_ERI_F3y_S_Px_S_bb;
  Double I_ERI_F2yz_Px_Px_S_bb = I_ERI_Gx2yz_S_Px_S_bb+ABX*I_ERI_F2yz_S_Px_S_bb;
  Double I_ERI_Fy2z_Px_Px_S_bb = I_ERI_Gxy2z_S_Px_S_bb+ABX*I_ERI_Fy2z_S_Px_S_bb;
  Double I_ERI_F3z_Px_Px_S_bb = I_ERI_Gx3z_S_Px_S_bb+ABX*I_ERI_F3z_S_Px_S_bb;
  Double I_ERI_F3x_Py_Px_S_bb = I_ERI_G3xy_S_Px_S_bb+ABY*I_ERI_F3x_S_Px_S_bb;
  Double I_ERI_F2xy_Py_Px_S_bb = I_ERI_G2x2y_S_Px_S_bb+ABY*I_ERI_F2xy_S_Px_S_bb;
  Double I_ERI_F2xz_Py_Px_S_bb = I_ERI_G2xyz_S_Px_S_bb+ABY*I_ERI_F2xz_S_Px_S_bb;
  Double I_ERI_Fx2y_Py_Px_S_bb = I_ERI_Gx3y_S_Px_S_bb+ABY*I_ERI_Fx2y_S_Px_S_bb;
  Double I_ERI_Fxyz_Py_Px_S_bb = I_ERI_Gx2yz_S_Px_S_bb+ABY*I_ERI_Fxyz_S_Px_S_bb;
  Double I_ERI_Fx2z_Py_Px_S_bb = I_ERI_Gxy2z_S_Px_S_bb+ABY*I_ERI_Fx2z_S_Px_S_bb;
  Double I_ERI_F3y_Py_Px_S_bb = I_ERI_G4y_S_Px_S_bb+ABY*I_ERI_F3y_S_Px_S_bb;
  Double I_ERI_F2yz_Py_Px_S_bb = I_ERI_G3yz_S_Px_S_bb+ABY*I_ERI_F2yz_S_Px_S_bb;
  Double I_ERI_Fy2z_Py_Px_S_bb = I_ERI_G2y2z_S_Px_S_bb+ABY*I_ERI_Fy2z_S_Px_S_bb;
  Double I_ERI_F3z_Py_Px_S_bb = I_ERI_Gy3z_S_Px_S_bb+ABY*I_ERI_F3z_S_Px_S_bb;
  Double I_ERI_F3x_Pz_Px_S_bb = I_ERI_G3xz_S_Px_S_bb+ABZ*I_ERI_F3x_S_Px_S_bb;
  Double I_ERI_F2xy_Pz_Px_S_bb = I_ERI_G2xyz_S_Px_S_bb+ABZ*I_ERI_F2xy_S_Px_S_bb;
  Double I_ERI_F2xz_Pz_Px_S_bb = I_ERI_G2x2z_S_Px_S_bb+ABZ*I_ERI_F2xz_S_Px_S_bb;
  Double I_ERI_Fx2y_Pz_Px_S_bb = I_ERI_Gx2yz_S_Px_S_bb+ABZ*I_ERI_Fx2y_S_Px_S_bb;
  Double I_ERI_Fxyz_Pz_Px_S_bb = I_ERI_Gxy2z_S_Px_S_bb+ABZ*I_ERI_Fxyz_S_Px_S_bb;
  Double I_ERI_Fx2z_Pz_Px_S_bb = I_ERI_Gx3z_S_Px_S_bb+ABZ*I_ERI_Fx2z_S_Px_S_bb;
  Double I_ERI_F3y_Pz_Px_S_bb = I_ERI_G3yz_S_Px_S_bb+ABZ*I_ERI_F3y_S_Px_S_bb;
  Double I_ERI_F2yz_Pz_Px_S_bb = I_ERI_G2y2z_S_Px_S_bb+ABZ*I_ERI_F2yz_S_Px_S_bb;
  Double I_ERI_Fy2z_Pz_Px_S_bb = I_ERI_Gy3z_S_Px_S_bb+ABZ*I_ERI_Fy2z_S_Px_S_bb;
  Double I_ERI_F3z_Pz_Px_S_bb = I_ERI_G4z_S_Px_S_bb+ABZ*I_ERI_F3z_S_Px_S_bb;
  Double I_ERI_F3x_Px_Py_S_bb = I_ERI_G4x_S_Py_S_bb+ABX*I_ERI_F3x_S_Py_S_bb;
  Double I_ERI_F2xy_Px_Py_S_bb = I_ERI_G3xy_S_Py_S_bb+ABX*I_ERI_F2xy_S_Py_S_bb;
  Double I_ERI_F2xz_Px_Py_S_bb = I_ERI_G3xz_S_Py_S_bb+ABX*I_ERI_F2xz_S_Py_S_bb;
  Double I_ERI_Fx2y_Px_Py_S_bb = I_ERI_G2x2y_S_Py_S_bb+ABX*I_ERI_Fx2y_S_Py_S_bb;
  Double I_ERI_Fxyz_Px_Py_S_bb = I_ERI_G2xyz_S_Py_S_bb+ABX*I_ERI_Fxyz_S_Py_S_bb;
  Double I_ERI_Fx2z_Px_Py_S_bb = I_ERI_G2x2z_S_Py_S_bb+ABX*I_ERI_Fx2z_S_Py_S_bb;
  Double I_ERI_F3y_Px_Py_S_bb = I_ERI_Gx3y_S_Py_S_bb+ABX*I_ERI_F3y_S_Py_S_bb;
  Double I_ERI_F2yz_Px_Py_S_bb = I_ERI_Gx2yz_S_Py_S_bb+ABX*I_ERI_F2yz_S_Py_S_bb;
  Double I_ERI_Fy2z_Px_Py_S_bb = I_ERI_Gxy2z_S_Py_S_bb+ABX*I_ERI_Fy2z_S_Py_S_bb;
  Double I_ERI_F3z_Px_Py_S_bb = I_ERI_Gx3z_S_Py_S_bb+ABX*I_ERI_F3z_S_Py_S_bb;
  Double I_ERI_F3x_Py_Py_S_bb = I_ERI_G3xy_S_Py_S_bb+ABY*I_ERI_F3x_S_Py_S_bb;
  Double I_ERI_F2xy_Py_Py_S_bb = I_ERI_G2x2y_S_Py_S_bb+ABY*I_ERI_F2xy_S_Py_S_bb;
  Double I_ERI_F2xz_Py_Py_S_bb = I_ERI_G2xyz_S_Py_S_bb+ABY*I_ERI_F2xz_S_Py_S_bb;
  Double I_ERI_Fx2y_Py_Py_S_bb = I_ERI_Gx3y_S_Py_S_bb+ABY*I_ERI_Fx2y_S_Py_S_bb;
  Double I_ERI_Fxyz_Py_Py_S_bb = I_ERI_Gx2yz_S_Py_S_bb+ABY*I_ERI_Fxyz_S_Py_S_bb;
  Double I_ERI_Fx2z_Py_Py_S_bb = I_ERI_Gxy2z_S_Py_S_bb+ABY*I_ERI_Fx2z_S_Py_S_bb;
  Double I_ERI_F3y_Py_Py_S_bb = I_ERI_G4y_S_Py_S_bb+ABY*I_ERI_F3y_S_Py_S_bb;
  Double I_ERI_F2yz_Py_Py_S_bb = I_ERI_G3yz_S_Py_S_bb+ABY*I_ERI_F2yz_S_Py_S_bb;
  Double I_ERI_Fy2z_Py_Py_S_bb = I_ERI_G2y2z_S_Py_S_bb+ABY*I_ERI_Fy2z_S_Py_S_bb;
  Double I_ERI_F3z_Py_Py_S_bb = I_ERI_Gy3z_S_Py_S_bb+ABY*I_ERI_F3z_S_Py_S_bb;
  Double I_ERI_F3x_Pz_Py_S_bb = I_ERI_G3xz_S_Py_S_bb+ABZ*I_ERI_F3x_S_Py_S_bb;
  Double I_ERI_F2xy_Pz_Py_S_bb = I_ERI_G2xyz_S_Py_S_bb+ABZ*I_ERI_F2xy_S_Py_S_bb;
  Double I_ERI_F2xz_Pz_Py_S_bb = I_ERI_G2x2z_S_Py_S_bb+ABZ*I_ERI_F2xz_S_Py_S_bb;
  Double I_ERI_Fx2y_Pz_Py_S_bb = I_ERI_Gx2yz_S_Py_S_bb+ABZ*I_ERI_Fx2y_S_Py_S_bb;
  Double I_ERI_Fxyz_Pz_Py_S_bb = I_ERI_Gxy2z_S_Py_S_bb+ABZ*I_ERI_Fxyz_S_Py_S_bb;
  Double I_ERI_Fx2z_Pz_Py_S_bb = I_ERI_Gx3z_S_Py_S_bb+ABZ*I_ERI_Fx2z_S_Py_S_bb;
  Double I_ERI_F3y_Pz_Py_S_bb = I_ERI_G3yz_S_Py_S_bb+ABZ*I_ERI_F3y_S_Py_S_bb;
  Double I_ERI_F2yz_Pz_Py_S_bb = I_ERI_G2y2z_S_Py_S_bb+ABZ*I_ERI_F2yz_S_Py_S_bb;
  Double I_ERI_Fy2z_Pz_Py_S_bb = I_ERI_Gy3z_S_Py_S_bb+ABZ*I_ERI_Fy2z_S_Py_S_bb;
  Double I_ERI_F3z_Pz_Py_S_bb = I_ERI_G4z_S_Py_S_bb+ABZ*I_ERI_F3z_S_Py_S_bb;
  Double I_ERI_F3x_Px_Pz_S_bb = I_ERI_G4x_S_Pz_S_bb+ABX*I_ERI_F3x_S_Pz_S_bb;
  Double I_ERI_F2xy_Px_Pz_S_bb = I_ERI_G3xy_S_Pz_S_bb+ABX*I_ERI_F2xy_S_Pz_S_bb;
  Double I_ERI_F2xz_Px_Pz_S_bb = I_ERI_G3xz_S_Pz_S_bb+ABX*I_ERI_F2xz_S_Pz_S_bb;
  Double I_ERI_Fx2y_Px_Pz_S_bb = I_ERI_G2x2y_S_Pz_S_bb+ABX*I_ERI_Fx2y_S_Pz_S_bb;
  Double I_ERI_Fxyz_Px_Pz_S_bb = I_ERI_G2xyz_S_Pz_S_bb+ABX*I_ERI_Fxyz_S_Pz_S_bb;
  Double I_ERI_Fx2z_Px_Pz_S_bb = I_ERI_G2x2z_S_Pz_S_bb+ABX*I_ERI_Fx2z_S_Pz_S_bb;
  Double I_ERI_F3y_Px_Pz_S_bb = I_ERI_Gx3y_S_Pz_S_bb+ABX*I_ERI_F3y_S_Pz_S_bb;
  Double I_ERI_F2yz_Px_Pz_S_bb = I_ERI_Gx2yz_S_Pz_S_bb+ABX*I_ERI_F2yz_S_Pz_S_bb;
  Double I_ERI_Fy2z_Px_Pz_S_bb = I_ERI_Gxy2z_S_Pz_S_bb+ABX*I_ERI_Fy2z_S_Pz_S_bb;
  Double I_ERI_F3z_Px_Pz_S_bb = I_ERI_Gx3z_S_Pz_S_bb+ABX*I_ERI_F3z_S_Pz_S_bb;
  Double I_ERI_F3x_Py_Pz_S_bb = I_ERI_G3xy_S_Pz_S_bb+ABY*I_ERI_F3x_S_Pz_S_bb;
  Double I_ERI_F2xy_Py_Pz_S_bb = I_ERI_G2x2y_S_Pz_S_bb+ABY*I_ERI_F2xy_S_Pz_S_bb;
  Double I_ERI_F2xz_Py_Pz_S_bb = I_ERI_G2xyz_S_Pz_S_bb+ABY*I_ERI_F2xz_S_Pz_S_bb;
  Double I_ERI_Fx2y_Py_Pz_S_bb = I_ERI_Gx3y_S_Pz_S_bb+ABY*I_ERI_Fx2y_S_Pz_S_bb;
  Double I_ERI_Fxyz_Py_Pz_S_bb = I_ERI_Gx2yz_S_Pz_S_bb+ABY*I_ERI_Fxyz_S_Pz_S_bb;
  Double I_ERI_Fx2z_Py_Pz_S_bb = I_ERI_Gxy2z_S_Pz_S_bb+ABY*I_ERI_Fx2z_S_Pz_S_bb;
  Double I_ERI_F3y_Py_Pz_S_bb = I_ERI_G4y_S_Pz_S_bb+ABY*I_ERI_F3y_S_Pz_S_bb;
  Double I_ERI_F2yz_Py_Pz_S_bb = I_ERI_G3yz_S_Pz_S_bb+ABY*I_ERI_F2yz_S_Pz_S_bb;
  Double I_ERI_Fy2z_Py_Pz_S_bb = I_ERI_G2y2z_S_Pz_S_bb+ABY*I_ERI_Fy2z_S_Pz_S_bb;
  Double I_ERI_F3z_Py_Pz_S_bb = I_ERI_Gy3z_S_Pz_S_bb+ABY*I_ERI_F3z_S_Pz_S_bb;
  Double I_ERI_F3x_Pz_Pz_S_bb = I_ERI_G3xz_S_Pz_S_bb+ABZ*I_ERI_F3x_S_Pz_S_bb;
  Double I_ERI_F2xy_Pz_Pz_S_bb = I_ERI_G2xyz_S_Pz_S_bb+ABZ*I_ERI_F2xy_S_Pz_S_bb;
  Double I_ERI_F2xz_Pz_Pz_S_bb = I_ERI_G2x2z_S_Pz_S_bb+ABZ*I_ERI_F2xz_S_Pz_S_bb;
  Double I_ERI_Fx2y_Pz_Pz_S_bb = I_ERI_Gx2yz_S_Pz_S_bb+ABZ*I_ERI_Fx2y_S_Pz_S_bb;
  Double I_ERI_Fxyz_Pz_Pz_S_bb = I_ERI_Gxy2z_S_Pz_S_bb+ABZ*I_ERI_Fxyz_S_Pz_S_bb;
  Double I_ERI_Fx2z_Pz_Pz_S_bb = I_ERI_Gx3z_S_Pz_S_bb+ABZ*I_ERI_Fx2z_S_Pz_S_bb;
  Double I_ERI_F3y_Pz_Pz_S_bb = I_ERI_G3yz_S_Pz_S_bb+ABZ*I_ERI_F3y_S_Pz_S_bb;
  Double I_ERI_F2yz_Pz_Pz_S_bb = I_ERI_G2y2z_S_Pz_S_bb+ABZ*I_ERI_F2yz_S_Pz_S_bb;
  Double I_ERI_Fy2z_Pz_Pz_S_bb = I_ERI_Gy3z_S_Pz_S_bb+ABZ*I_ERI_Fy2z_S_Pz_S_bb;
  Double I_ERI_F3z_Pz_Pz_S_bb = I_ERI_G4z_S_Pz_S_bb+ABZ*I_ERI_F3z_S_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_P_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 36 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_bb
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bb
   ************************************************************/
  Double I_ERI_D2x_D2x_Px_S_bb = I_ERI_F3x_Px_Px_S_bb+ABX*I_ERI_D2x_Px_Px_S_bb;
  Double I_ERI_Dxy_D2x_Px_S_bb = I_ERI_F2xy_Px_Px_S_bb+ABX*I_ERI_Dxy_Px_Px_S_bb;
  Double I_ERI_Dxz_D2x_Px_S_bb = I_ERI_F2xz_Px_Px_S_bb+ABX*I_ERI_Dxz_Px_Px_S_bb;
  Double I_ERI_D2y_D2x_Px_S_bb = I_ERI_Fx2y_Px_Px_S_bb+ABX*I_ERI_D2y_Px_Px_S_bb;
  Double I_ERI_Dyz_D2x_Px_S_bb = I_ERI_Fxyz_Px_Px_S_bb+ABX*I_ERI_Dyz_Px_Px_S_bb;
  Double I_ERI_D2z_D2x_Px_S_bb = I_ERI_Fx2z_Px_Px_S_bb+ABX*I_ERI_D2z_Px_Px_S_bb;
  Double I_ERI_D2x_Dxy_Px_S_bb = I_ERI_F2xy_Px_Px_S_bb+ABY*I_ERI_D2x_Px_Px_S_bb;
  Double I_ERI_Dxy_Dxy_Px_S_bb = I_ERI_Fx2y_Px_Px_S_bb+ABY*I_ERI_Dxy_Px_Px_S_bb;
  Double I_ERI_Dxz_Dxy_Px_S_bb = I_ERI_Fxyz_Px_Px_S_bb+ABY*I_ERI_Dxz_Px_Px_S_bb;
  Double I_ERI_D2y_Dxy_Px_S_bb = I_ERI_F3y_Px_Px_S_bb+ABY*I_ERI_D2y_Px_Px_S_bb;
  Double I_ERI_Dyz_Dxy_Px_S_bb = I_ERI_F2yz_Px_Px_S_bb+ABY*I_ERI_Dyz_Px_Px_S_bb;
  Double I_ERI_D2z_Dxy_Px_S_bb = I_ERI_Fy2z_Px_Px_S_bb+ABY*I_ERI_D2z_Px_Px_S_bb;
  Double I_ERI_D2x_D2y_Px_S_bb = I_ERI_F2xy_Py_Px_S_bb+ABY*I_ERI_D2x_Py_Px_S_bb;
  Double I_ERI_Dxy_D2y_Px_S_bb = I_ERI_Fx2y_Py_Px_S_bb+ABY*I_ERI_Dxy_Py_Px_S_bb;
  Double I_ERI_Dxz_D2y_Px_S_bb = I_ERI_Fxyz_Py_Px_S_bb+ABY*I_ERI_Dxz_Py_Px_S_bb;
  Double I_ERI_D2y_D2y_Px_S_bb = I_ERI_F3y_Py_Px_S_bb+ABY*I_ERI_D2y_Py_Px_S_bb;
  Double I_ERI_Dyz_D2y_Px_S_bb = I_ERI_F2yz_Py_Px_S_bb+ABY*I_ERI_Dyz_Py_Px_S_bb;
  Double I_ERI_D2z_D2y_Px_S_bb = I_ERI_Fy2z_Py_Px_S_bb+ABY*I_ERI_D2z_Py_Px_S_bb;
  Double I_ERI_D2x_D2z_Px_S_bb = I_ERI_F2xz_Pz_Px_S_bb+ABZ*I_ERI_D2x_Pz_Px_S_bb;
  Double I_ERI_Dxy_D2z_Px_S_bb = I_ERI_Fxyz_Pz_Px_S_bb+ABZ*I_ERI_Dxy_Pz_Px_S_bb;
  Double I_ERI_Dxz_D2z_Px_S_bb = I_ERI_Fx2z_Pz_Px_S_bb+ABZ*I_ERI_Dxz_Pz_Px_S_bb;
  Double I_ERI_D2y_D2z_Px_S_bb = I_ERI_F2yz_Pz_Px_S_bb+ABZ*I_ERI_D2y_Pz_Px_S_bb;
  Double I_ERI_Dyz_D2z_Px_S_bb = I_ERI_Fy2z_Pz_Px_S_bb+ABZ*I_ERI_Dyz_Pz_Px_S_bb;
  Double I_ERI_D2z_D2z_Px_S_bb = I_ERI_F3z_Pz_Px_S_bb+ABZ*I_ERI_D2z_Pz_Px_S_bb;
  Double I_ERI_D2x_D2x_Py_S_bb = I_ERI_F3x_Px_Py_S_bb+ABX*I_ERI_D2x_Px_Py_S_bb;
  Double I_ERI_Dxy_D2x_Py_S_bb = I_ERI_F2xy_Px_Py_S_bb+ABX*I_ERI_Dxy_Px_Py_S_bb;
  Double I_ERI_Dxz_D2x_Py_S_bb = I_ERI_F2xz_Px_Py_S_bb+ABX*I_ERI_Dxz_Px_Py_S_bb;
  Double I_ERI_D2y_D2x_Py_S_bb = I_ERI_Fx2y_Px_Py_S_bb+ABX*I_ERI_D2y_Px_Py_S_bb;
  Double I_ERI_Dyz_D2x_Py_S_bb = I_ERI_Fxyz_Px_Py_S_bb+ABX*I_ERI_Dyz_Px_Py_S_bb;
  Double I_ERI_D2z_D2x_Py_S_bb = I_ERI_Fx2z_Px_Py_S_bb+ABX*I_ERI_D2z_Px_Py_S_bb;
  Double I_ERI_D2x_Dxy_Py_S_bb = I_ERI_F2xy_Px_Py_S_bb+ABY*I_ERI_D2x_Px_Py_S_bb;
  Double I_ERI_Dxy_Dxy_Py_S_bb = I_ERI_Fx2y_Px_Py_S_bb+ABY*I_ERI_Dxy_Px_Py_S_bb;
  Double I_ERI_Dxz_Dxy_Py_S_bb = I_ERI_Fxyz_Px_Py_S_bb+ABY*I_ERI_Dxz_Px_Py_S_bb;
  Double I_ERI_D2y_Dxy_Py_S_bb = I_ERI_F3y_Px_Py_S_bb+ABY*I_ERI_D2y_Px_Py_S_bb;
  Double I_ERI_Dyz_Dxy_Py_S_bb = I_ERI_F2yz_Px_Py_S_bb+ABY*I_ERI_Dyz_Px_Py_S_bb;
  Double I_ERI_D2z_Dxy_Py_S_bb = I_ERI_Fy2z_Px_Py_S_bb+ABY*I_ERI_D2z_Px_Py_S_bb;
  Double I_ERI_D2x_D2y_Py_S_bb = I_ERI_F2xy_Py_Py_S_bb+ABY*I_ERI_D2x_Py_Py_S_bb;
  Double I_ERI_Dxy_D2y_Py_S_bb = I_ERI_Fx2y_Py_Py_S_bb+ABY*I_ERI_Dxy_Py_Py_S_bb;
  Double I_ERI_Dxz_D2y_Py_S_bb = I_ERI_Fxyz_Py_Py_S_bb+ABY*I_ERI_Dxz_Py_Py_S_bb;
  Double I_ERI_D2y_D2y_Py_S_bb = I_ERI_F3y_Py_Py_S_bb+ABY*I_ERI_D2y_Py_Py_S_bb;
  Double I_ERI_Dyz_D2y_Py_S_bb = I_ERI_F2yz_Py_Py_S_bb+ABY*I_ERI_Dyz_Py_Py_S_bb;
  Double I_ERI_D2z_D2y_Py_S_bb = I_ERI_Fy2z_Py_Py_S_bb+ABY*I_ERI_D2z_Py_Py_S_bb;
  Double I_ERI_D2x_D2z_Py_S_bb = I_ERI_F2xz_Pz_Py_S_bb+ABZ*I_ERI_D2x_Pz_Py_S_bb;
  Double I_ERI_Dxy_D2z_Py_S_bb = I_ERI_Fxyz_Pz_Py_S_bb+ABZ*I_ERI_Dxy_Pz_Py_S_bb;
  Double I_ERI_Dxz_D2z_Py_S_bb = I_ERI_Fx2z_Pz_Py_S_bb+ABZ*I_ERI_Dxz_Pz_Py_S_bb;
  Double I_ERI_D2y_D2z_Py_S_bb = I_ERI_F2yz_Pz_Py_S_bb+ABZ*I_ERI_D2y_Pz_Py_S_bb;
  Double I_ERI_Dyz_D2z_Py_S_bb = I_ERI_Fy2z_Pz_Py_S_bb+ABZ*I_ERI_Dyz_Pz_Py_S_bb;
  Double I_ERI_D2z_D2z_Py_S_bb = I_ERI_F3z_Pz_Py_S_bb+ABZ*I_ERI_D2z_Pz_Py_S_bb;
  Double I_ERI_D2x_D2x_Pz_S_bb = I_ERI_F3x_Px_Pz_S_bb+ABX*I_ERI_D2x_Px_Pz_S_bb;
  Double I_ERI_Dxy_D2x_Pz_S_bb = I_ERI_F2xy_Px_Pz_S_bb+ABX*I_ERI_Dxy_Px_Pz_S_bb;
  Double I_ERI_Dxz_D2x_Pz_S_bb = I_ERI_F2xz_Px_Pz_S_bb+ABX*I_ERI_Dxz_Px_Pz_S_bb;
  Double I_ERI_D2y_D2x_Pz_S_bb = I_ERI_Fx2y_Px_Pz_S_bb+ABX*I_ERI_D2y_Px_Pz_S_bb;
  Double I_ERI_Dyz_D2x_Pz_S_bb = I_ERI_Fxyz_Px_Pz_S_bb+ABX*I_ERI_Dyz_Px_Pz_S_bb;
  Double I_ERI_D2z_D2x_Pz_S_bb = I_ERI_Fx2z_Px_Pz_S_bb+ABX*I_ERI_D2z_Px_Pz_S_bb;
  Double I_ERI_D2x_Dxy_Pz_S_bb = I_ERI_F2xy_Px_Pz_S_bb+ABY*I_ERI_D2x_Px_Pz_S_bb;
  Double I_ERI_Dxy_Dxy_Pz_S_bb = I_ERI_Fx2y_Px_Pz_S_bb+ABY*I_ERI_Dxy_Px_Pz_S_bb;
  Double I_ERI_Dxz_Dxy_Pz_S_bb = I_ERI_Fxyz_Px_Pz_S_bb+ABY*I_ERI_Dxz_Px_Pz_S_bb;
  Double I_ERI_D2y_Dxy_Pz_S_bb = I_ERI_F3y_Px_Pz_S_bb+ABY*I_ERI_D2y_Px_Pz_S_bb;
  Double I_ERI_Dyz_Dxy_Pz_S_bb = I_ERI_F2yz_Px_Pz_S_bb+ABY*I_ERI_Dyz_Px_Pz_S_bb;
  Double I_ERI_D2z_Dxy_Pz_S_bb = I_ERI_Fy2z_Px_Pz_S_bb+ABY*I_ERI_D2z_Px_Pz_S_bb;
  Double I_ERI_D2x_D2y_Pz_S_bb = I_ERI_F2xy_Py_Pz_S_bb+ABY*I_ERI_D2x_Py_Pz_S_bb;
  Double I_ERI_Dxy_D2y_Pz_S_bb = I_ERI_Fx2y_Py_Pz_S_bb+ABY*I_ERI_Dxy_Py_Pz_S_bb;
  Double I_ERI_Dxz_D2y_Pz_S_bb = I_ERI_Fxyz_Py_Pz_S_bb+ABY*I_ERI_Dxz_Py_Pz_S_bb;
  Double I_ERI_D2y_D2y_Pz_S_bb = I_ERI_F3y_Py_Pz_S_bb+ABY*I_ERI_D2y_Py_Pz_S_bb;
  Double I_ERI_Dyz_D2y_Pz_S_bb = I_ERI_F2yz_Py_Pz_S_bb+ABY*I_ERI_Dyz_Py_Pz_S_bb;
  Double I_ERI_D2z_D2y_Pz_S_bb = I_ERI_Fy2z_Py_Pz_S_bb+ABY*I_ERI_D2z_Py_Pz_S_bb;
  Double I_ERI_D2x_D2z_Pz_S_bb = I_ERI_F2xz_Pz_Pz_S_bb+ABZ*I_ERI_D2x_Pz_Pz_S_bb;
  Double I_ERI_Dxy_D2z_Pz_S_bb = I_ERI_Fxyz_Pz_Pz_S_bb+ABZ*I_ERI_Dxy_Pz_Pz_S_bb;
  Double I_ERI_Dxz_D2z_Pz_S_bb = I_ERI_Fx2z_Pz_Pz_S_bb+ABZ*I_ERI_Dxz_Pz_Pz_S_bb;
  Double I_ERI_D2y_D2z_Pz_S_bb = I_ERI_F2yz_Pz_Pz_S_bb+ABZ*I_ERI_D2y_Pz_Pz_S_bb;
  Double I_ERI_Dyz_D2z_Pz_S_bb = I_ERI_Fy2z_Pz_Pz_S_bb+ABZ*I_ERI_Dyz_Pz_Pz_S_bb;
  Double I_ERI_D2z_D2z_Pz_S_bb = I_ERI_F3z_Pz_Pz_S_bb+ABZ*I_ERI_D2z_Pz_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_P_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 36 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_P_S_bb
   * RHS shell quartet name: SQ_ERI_G_S_P_S_bb
   ************************************************************/
  Double I_ERI_G4x_Px_Px_S_bb = I_ERI_H5x_S_Px_S_bb+ABX*I_ERI_G4x_S_Px_S_bb;
  Double I_ERI_G3xy_Px_Px_S_bb = I_ERI_H4xy_S_Px_S_bb+ABX*I_ERI_G3xy_S_Px_S_bb;
  Double I_ERI_G3xz_Px_Px_S_bb = I_ERI_H4xz_S_Px_S_bb+ABX*I_ERI_G3xz_S_Px_S_bb;
  Double I_ERI_G2x2y_Px_Px_S_bb = I_ERI_H3x2y_S_Px_S_bb+ABX*I_ERI_G2x2y_S_Px_S_bb;
  Double I_ERI_G2xyz_Px_Px_S_bb = I_ERI_H3xyz_S_Px_S_bb+ABX*I_ERI_G2xyz_S_Px_S_bb;
  Double I_ERI_G2x2z_Px_Px_S_bb = I_ERI_H3x2z_S_Px_S_bb+ABX*I_ERI_G2x2z_S_Px_S_bb;
  Double I_ERI_Gx3y_Px_Px_S_bb = I_ERI_H2x3y_S_Px_S_bb+ABX*I_ERI_Gx3y_S_Px_S_bb;
  Double I_ERI_Gx2yz_Px_Px_S_bb = I_ERI_H2x2yz_S_Px_S_bb+ABX*I_ERI_Gx2yz_S_Px_S_bb;
  Double I_ERI_Gxy2z_Px_Px_S_bb = I_ERI_H2xy2z_S_Px_S_bb+ABX*I_ERI_Gxy2z_S_Px_S_bb;
  Double I_ERI_Gx3z_Px_Px_S_bb = I_ERI_H2x3z_S_Px_S_bb+ABX*I_ERI_Gx3z_S_Px_S_bb;
  Double I_ERI_G3yz_Px_Px_S_bb = I_ERI_Hx3yz_S_Px_S_bb+ABX*I_ERI_G3yz_S_Px_S_bb;
  Double I_ERI_G2y2z_Px_Px_S_bb = I_ERI_Hx2y2z_S_Px_S_bb+ABX*I_ERI_G2y2z_S_Px_S_bb;
  Double I_ERI_Gy3z_Px_Px_S_bb = I_ERI_Hxy3z_S_Px_S_bb+ABX*I_ERI_Gy3z_S_Px_S_bb;
  Double I_ERI_G3xy_Py_Px_S_bb = I_ERI_H3x2y_S_Px_S_bb+ABY*I_ERI_G3xy_S_Px_S_bb;
  Double I_ERI_G2x2y_Py_Px_S_bb = I_ERI_H2x3y_S_Px_S_bb+ABY*I_ERI_G2x2y_S_Px_S_bb;
  Double I_ERI_G2xyz_Py_Px_S_bb = I_ERI_H2x2yz_S_Px_S_bb+ABY*I_ERI_G2xyz_S_Px_S_bb;
  Double I_ERI_Gx3y_Py_Px_S_bb = I_ERI_Hx4y_S_Px_S_bb+ABY*I_ERI_Gx3y_S_Px_S_bb;
  Double I_ERI_Gx2yz_Py_Px_S_bb = I_ERI_Hx3yz_S_Px_S_bb+ABY*I_ERI_Gx2yz_S_Px_S_bb;
  Double I_ERI_Gxy2z_Py_Px_S_bb = I_ERI_Hx2y2z_S_Px_S_bb+ABY*I_ERI_Gxy2z_S_Px_S_bb;
  Double I_ERI_G4y_Py_Px_S_bb = I_ERI_H5y_S_Px_S_bb+ABY*I_ERI_G4y_S_Px_S_bb;
  Double I_ERI_G3yz_Py_Px_S_bb = I_ERI_H4yz_S_Px_S_bb+ABY*I_ERI_G3yz_S_Px_S_bb;
  Double I_ERI_G2y2z_Py_Px_S_bb = I_ERI_H3y2z_S_Px_S_bb+ABY*I_ERI_G2y2z_S_Px_S_bb;
  Double I_ERI_Gy3z_Py_Px_S_bb = I_ERI_H2y3z_S_Px_S_bb+ABY*I_ERI_Gy3z_S_Px_S_bb;
  Double I_ERI_G3xz_Pz_Px_S_bb = I_ERI_H3x2z_S_Px_S_bb+ABZ*I_ERI_G3xz_S_Px_S_bb;
  Double I_ERI_G2xyz_Pz_Px_S_bb = I_ERI_H2xy2z_S_Px_S_bb+ABZ*I_ERI_G2xyz_S_Px_S_bb;
  Double I_ERI_G2x2z_Pz_Px_S_bb = I_ERI_H2x3z_S_Px_S_bb+ABZ*I_ERI_G2x2z_S_Px_S_bb;
  Double I_ERI_Gx2yz_Pz_Px_S_bb = I_ERI_Hx2y2z_S_Px_S_bb+ABZ*I_ERI_Gx2yz_S_Px_S_bb;
  Double I_ERI_Gxy2z_Pz_Px_S_bb = I_ERI_Hxy3z_S_Px_S_bb+ABZ*I_ERI_Gxy2z_S_Px_S_bb;
  Double I_ERI_Gx3z_Pz_Px_S_bb = I_ERI_Hx4z_S_Px_S_bb+ABZ*I_ERI_Gx3z_S_Px_S_bb;
  Double I_ERI_G3yz_Pz_Px_S_bb = I_ERI_H3y2z_S_Px_S_bb+ABZ*I_ERI_G3yz_S_Px_S_bb;
  Double I_ERI_G2y2z_Pz_Px_S_bb = I_ERI_H2y3z_S_Px_S_bb+ABZ*I_ERI_G2y2z_S_Px_S_bb;
  Double I_ERI_Gy3z_Pz_Px_S_bb = I_ERI_Hy4z_S_Px_S_bb+ABZ*I_ERI_Gy3z_S_Px_S_bb;
  Double I_ERI_G4z_Pz_Px_S_bb = I_ERI_H5z_S_Px_S_bb+ABZ*I_ERI_G4z_S_Px_S_bb;
  Double I_ERI_G4x_Px_Py_S_bb = I_ERI_H5x_S_Py_S_bb+ABX*I_ERI_G4x_S_Py_S_bb;
  Double I_ERI_G3xy_Px_Py_S_bb = I_ERI_H4xy_S_Py_S_bb+ABX*I_ERI_G3xy_S_Py_S_bb;
  Double I_ERI_G3xz_Px_Py_S_bb = I_ERI_H4xz_S_Py_S_bb+ABX*I_ERI_G3xz_S_Py_S_bb;
  Double I_ERI_G2x2y_Px_Py_S_bb = I_ERI_H3x2y_S_Py_S_bb+ABX*I_ERI_G2x2y_S_Py_S_bb;
  Double I_ERI_G2xyz_Px_Py_S_bb = I_ERI_H3xyz_S_Py_S_bb+ABX*I_ERI_G2xyz_S_Py_S_bb;
  Double I_ERI_G2x2z_Px_Py_S_bb = I_ERI_H3x2z_S_Py_S_bb+ABX*I_ERI_G2x2z_S_Py_S_bb;
  Double I_ERI_Gx3y_Px_Py_S_bb = I_ERI_H2x3y_S_Py_S_bb+ABX*I_ERI_Gx3y_S_Py_S_bb;
  Double I_ERI_Gx2yz_Px_Py_S_bb = I_ERI_H2x2yz_S_Py_S_bb+ABX*I_ERI_Gx2yz_S_Py_S_bb;
  Double I_ERI_Gxy2z_Px_Py_S_bb = I_ERI_H2xy2z_S_Py_S_bb+ABX*I_ERI_Gxy2z_S_Py_S_bb;
  Double I_ERI_Gx3z_Px_Py_S_bb = I_ERI_H2x3z_S_Py_S_bb+ABX*I_ERI_Gx3z_S_Py_S_bb;
  Double I_ERI_G3yz_Px_Py_S_bb = I_ERI_Hx3yz_S_Py_S_bb+ABX*I_ERI_G3yz_S_Py_S_bb;
  Double I_ERI_G2y2z_Px_Py_S_bb = I_ERI_Hx2y2z_S_Py_S_bb+ABX*I_ERI_G2y2z_S_Py_S_bb;
  Double I_ERI_Gy3z_Px_Py_S_bb = I_ERI_Hxy3z_S_Py_S_bb+ABX*I_ERI_Gy3z_S_Py_S_bb;
  Double I_ERI_G3xy_Py_Py_S_bb = I_ERI_H3x2y_S_Py_S_bb+ABY*I_ERI_G3xy_S_Py_S_bb;
  Double I_ERI_G2x2y_Py_Py_S_bb = I_ERI_H2x3y_S_Py_S_bb+ABY*I_ERI_G2x2y_S_Py_S_bb;
  Double I_ERI_G2xyz_Py_Py_S_bb = I_ERI_H2x2yz_S_Py_S_bb+ABY*I_ERI_G2xyz_S_Py_S_bb;
  Double I_ERI_Gx3y_Py_Py_S_bb = I_ERI_Hx4y_S_Py_S_bb+ABY*I_ERI_Gx3y_S_Py_S_bb;
  Double I_ERI_Gx2yz_Py_Py_S_bb = I_ERI_Hx3yz_S_Py_S_bb+ABY*I_ERI_Gx2yz_S_Py_S_bb;
  Double I_ERI_Gxy2z_Py_Py_S_bb = I_ERI_Hx2y2z_S_Py_S_bb+ABY*I_ERI_Gxy2z_S_Py_S_bb;
  Double I_ERI_G4y_Py_Py_S_bb = I_ERI_H5y_S_Py_S_bb+ABY*I_ERI_G4y_S_Py_S_bb;
  Double I_ERI_G3yz_Py_Py_S_bb = I_ERI_H4yz_S_Py_S_bb+ABY*I_ERI_G3yz_S_Py_S_bb;
  Double I_ERI_G2y2z_Py_Py_S_bb = I_ERI_H3y2z_S_Py_S_bb+ABY*I_ERI_G2y2z_S_Py_S_bb;
  Double I_ERI_Gy3z_Py_Py_S_bb = I_ERI_H2y3z_S_Py_S_bb+ABY*I_ERI_Gy3z_S_Py_S_bb;
  Double I_ERI_G3xz_Pz_Py_S_bb = I_ERI_H3x2z_S_Py_S_bb+ABZ*I_ERI_G3xz_S_Py_S_bb;
  Double I_ERI_G2xyz_Pz_Py_S_bb = I_ERI_H2xy2z_S_Py_S_bb+ABZ*I_ERI_G2xyz_S_Py_S_bb;
  Double I_ERI_G2x2z_Pz_Py_S_bb = I_ERI_H2x3z_S_Py_S_bb+ABZ*I_ERI_G2x2z_S_Py_S_bb;
  Double I_ERI_Gx2yz_Pz_Py_S_bb = I_ERI_Hx2y2z_S_Py_S_bb+ABZ*I_ERI_Gx2yz_S_Py_S_bb;
  Double I_ERI_Gxy2z_Pz_Py_S_bb = I_ERI_Hxy3z_S_Py_S_bb+ABZ*I_ERI_Gxy2z_S_Py_S_bb;
  Double I_ERI_Gx3z_Pz_Py_S_bb = I_ERI_Hx4z_S_Py_S_bb+ABZ*I_ERI_Gx3z_S_Py_S_bb;
  Double I_ERI_G3yz_Pz_Py_S_bb = I_ERI_H3y2z_S_Py_S_bb+ABZ*I_ERI_G3yz_S_Py_S_bb;
  Double I_ERI_G2y2z_Pz_Py_S_bb = I_ERI_H2y3z_S_Py_S_bb+ABZ*I_ERI_G2y2z_S_Py_S_bb;
  Double I_ERI_Gy3z_Pz_Py_S_bb = I_ERI_Hy4z_S_Py_S_bb+ABZ*I_ERI_Gy3z_S_Py_S_bb;
  Double I_ERI_G4z_Pz_Py_S_bb = I_ERI_H5z_S_Py_S_bb+ABZ*I_ERI_G4z_S_Py_S_bb;
  Double I_ERI_G4x_Px_Pz_S_bb = I_ERI_H5x_S_Pz_S_bb+ABX*I_ERI_G4x_S_Pz_S_bb;
  Double I_ERI_G3xy_Px_Pz_S_bb = I_ERI_H4xy_S_Pz_S_bb+ABX*I_ERI_G3xy_S_Pz_S_bb;
  Double I_ERI_G3xz_Px_Pz_S_bb = I_ERI_H4xz_S_Pz_S_bb+ABX*I_ERI_G3xz_S_Pz_S_bb;
  Double I_ERI_G2x2y_Px_Pz_S_bb = I_ERI_H3x2y_S_Pz_S_bb+ABX*I_ERI_G2x2y_S_Pz_S_bb;
  Double I_ERI_G2xyz_Px_Pz_S_bb = I_ERI_H3xyz_S_Pz_S_bb+ABX*I_ERI_G2xyz_S_Pz_S_bb;
  Double I_ERI_G2x2z_Px_Pz_S_bb = I_ERI_H3x2z_S_Pz_S_bb+ABX*I_ERI_G2x2z_S_Pz_S_bb;
  Double I_ERI_Gx3y_Px_Pz_S_bb = I_ERI_H2x3y_S_Pz_S_bb+ABX*I_ERI_Gx3y_S_Pz_S_bb;
  Double I_ERI_Gx2yz_Px_Pz_S_bb = I_ERI_H2x2yz_S_Pz_S_bb+ABX*I_ERI_Gx2yz_S_Pz_S_bb;
  Double I_ERI_Gxy2z_Px_Pz_S_bb = I_ERI_H2xy2z_S_Pz_S_bb+ABX*I_ERI_Gxy2z_S_Pz_S_bb;
  Double I_ERI_Gx3z_Px_Pz_S_bb = I_ERI_H2x3z_S_Pz_S_bb+ABX*I_ERI_Gx3z_S_Pz_S_bb;
  Double I_ERI_G3yz_Px_Pz_S_bb = I_ERI_Hx3yz_S_Pz_S_bb+ABX*I_ERI_G3yz_S_Pz_S_bb;
  Double I_ERI_G2y2z_Px_Pz_S_bb = I_ERI_Hx2y2z_S_Pz_S_bb+ABX*I_ERI_G2y2z_S_Pz_S_bb;
  Double I_ERI_Gy3z_Px_Pz_S_bb = I_ERI_Hxy3z_S_Pz_S_bb+ABX*I_ERI_Gy3z_S_Pz_S_bb;
  Double I_ERI_G3xy_Py_Pz_S_bb = I_ERI_H3x2y_S_Pz_S_bb+ABY*I_ERI_G3xy_S_Pz_S_bb;
  Double I_ERI_G2x2y_Py_Pz_S_bb = I_ERI_H2x3y_S_Pz_S_bb+ABY*I_ERI_G2x2y_S_Pz_S_bb;
  Double I_ERI_G2xyz_Py_Pz_S_bb = I_ERI_H2x2yz_S_Pz_S_bb+ABY*I_ERI_G2xyz_S_Pz_S_bb;
  Double I_ERI_Gx3y_Py_Pz_S_bb = I_ERI_Hx4y_S_Pz_S_bb+ABY*I_ERI_Gx3y_S_Pz_S_bb;
  Double I_ERI_Gx2yz_Py_Pz_S_bb = I_ERI_Hx3yz_S_Pz_S_bb+ABY*I_ERI_Gx2yz_S_Pz_S_bb;
  Double I_ERI_Gxy2z_Py_Pz_S_bb = I_ERI_Hx2y2z_S_Pz_S_bb+ABY*I_ERI_Gxy2z_S_Pz_S_bb;
  Double I_ERI_G4y_Py_Pz_S_bb = I_ERI_H5y_S_Pz_S_bb+ABY*I_ERI_G4y_S_Pz_S_bb;
  Double I_ERI_G3yz_Py_Pz_S_bb = I_ERI_H4yz_S_Pz_S_bb+ABY*I_ERI_G3yz_S_Pz_S_bb;
  Double I_ERI_G2y2z_Py_Pz_S_bb = I_ERI_H3y2z_S_Pz_S_bb+ABY*I_ERI_G2y2z_S_Pz_S_bb;
  Double I_ERI_Gy3z_Py_Pz_S_bb = I_ERI_H2y3z_S_Pz_S_bb+ABY*I_ERI_Gy3z_S_Pz_S_bb;
  Double I_ERI_G3xz_Pz_Pz_S_bb = I_ERI_H3x2z_S_Pz_S_bb+ABZ*I_ERI_G3xz_S_Pz_S_bb;
  Double I_ERI_G2xyz_Pz_Pz_S_bb = I_ERI_H2xy2z_S_Pz_S_bb+ABZ*I_ERI_G2xyz_S_Pz_S_bb;
  Double I_ERI_G2x2z_Pz_Pz_S_bb = I_ERI_H2x3z_S_Pz_S_bb+ABZ*I_ERI_G2x2z_S_Pz_S_bb;
  Double I_ERI_Gx2yz_Pz_Pz_S_bb = I_ERI_Hx2y2z_S_Pz_S_bb+ABZ*I_ERI_Gx2yz_S_Pz_S_bb;
  Double I_ERI_Gxy2z_Pz_Pz_S_bb = I_ERI_Hxy3z_S_Pz_S_bb+ABZ*I_ERI_Gxy2z_S_Pz_S_bb;
  Double I_ERI_Gx3z_Pz_Pz_S_bb = I_ERI_Hx4z_S_Pz_S_bb+ABZ*I_ERI_Gx3z_S_Pz_S_bb;
  Double I_ERI_G3yz_Pz_Pz_S_bb = I_ERI_H3y2z_S_Pz_S_bb+ABZ*I_ERI_G3yz_S_Pz_S_bb;
  Double I_ERI_G2y2z_Pz_Pz_S_bb = I_ERI_H2y3z_S_Pz_S_bb+ABZ*I_ERI_G2y2z_S_Pz_S_bb;
  Double I_ERI_Gy3z_Pz_Pz_S_bb = I_ERI_Hy4z_S_Pz_S_bb+ABZ*I_ERI_Gy3z_S_Pz_S_bb;
  Double I_ERI_G4z_Pz_Pz_S_bb = I_ERI_H5z_S_Pz_S_bb+ABZ*I_ERI_G4z_S_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_P_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 72 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_bb
   * RHS shell quartet name: SQ_ERI_F_P_P_S_bb
   ************************************************************/
  Double I_ERI_F3x_D2x_Px_S_bb = I_ERI_G4x_Px_Px_S_bb+ABX*I_ERI_F3x_Px_Px_S_bb;
  Double I_ERI_F2xy_D2x_Px_S_bb = I_ERI_G3xy_Px_Px_S_bb+ABX*I_ERI_F2xy_Px_Px_S_bb;
  Double I_ERI_F2xz_D2x_Px_S_bb = I_ERI_G3xz_Px_Px_S_bb+ABX*I_ERI_F2xz_Px_Px_S_bb;
  Double I_ERI_Fx2y_D2x_Px_S_bb = I_ERI_G2x2y_Px_Px_S_bb+ABX*I_ERI_Fx2y_Px_Px_S_bb;
  Double I_ERI_Fxyz_D2x_Px_S_bb = I_ERI_G2xyz_Px_Px_S_bb+ABX*I_ERI_Fxyz_Px_Px_S_bb;
  Double I_ERI_Fx2z_D2x_Px_S_bb = I_ERI_G2x2z_Px_Px_S_bb+ABX*I_ERI_Fx2z_Px_Px_S_bb;
  Double I_ERI_F3y_D2x_Px_S_bb = I_ERI_Gx3y_Px_Px_S_bb+ABX*I_ERI_F3y_Px_Px_S_bb;
  Double I_ERI_F2yz_D2x_Px_S_bb = I_ERI_Gx2yz_Px_Px_S_bb+ABX*I_ERI_F2yz_Px_Px_S_bb;
  Double I_ERI_Fy2z_D2x_Px_S_bb = I_ERI_Gxy2z_Px_Px_S_bb+ABX*I_ERI_Fy2z_Px_Px_S_bb;
  Double I_ERI_F3z_D2x_Px_S_bb = I_ERI_Gx3z_Px_Px_S_bb+ABX*I_ERI_F3z_Px_Px_S_bb;
  Double I_ERI_F2xz_Dxy_Px_S_bb = I_ERI_G2xyz_Px_Px_S_bb+ABY*I_ERI_F2xz_Px_Px_S_bb;
  Double I_ERI_Fxyz_Dxy_Px_S_bb = I_ERI_Gx2yz_Px_Px_S_bb+ABY*I_ERI_Fxyz_Px_Px_S_bb;
  Double I_ERI_Fx2z_Dxy_Px_S_bb = I_ERI_Gxy2z_Px_Px_S_bb+ABY*I_ERI_Fx2z_Px_Px_S_bb;
  Double I_ERI_F2yz_Dxy_Px_S_bb = I_ERI_G3yz_Px_Px_S_bb+ABY*I_ERI_F2yz_Px_Px_S_bb;
  Double I_ERI_Fy2z_Dxy_Px_S_bb = I_ERI_G2y2z_Px_Px_S_bb+ABY*I_ERI_Fy2z_Px_Px_S_bb;
  Double I_ERI_F3z_Dxy_Px_S_bb = I_ERI_Gy3z_Px_Px_S_bb+ABY*I_ERI_F3z_Px_Px_S_bb;
  Double I_ERI_F3x_D2y_Px_S_bb = I_ERI_G3xy_Py_Px_S_bb+ABY*I_ERI_F3x_Py_Px_S_bb;
  Double I_ERI_F2xy_D2y_Px_S_bb = I_ERI_G2x2y_Py_Px_S_bb+ABY*I_ERI_F2xy_Py_Px_S_bb;
  Double I_ERI_F2xz_D2y_Px_S_bb = I_ERI_G2xyz_Py_Px_S_bb+ABY*I_ERI_F2xz_Py_Px_S_bb;
  Double I_ERI_Fx2y_D2y_Px_S_bb = I_ERI_Gx3y_Py_Px_S_bb+ABY*I_ERI_Fx2y_Py_Px_S_bb;
  Double I_ERI_Fxyz_D2y_Px_S_bb = I_ERI_Gx2yz_Py_Px_S_bb+ABY*I_ERI_Fxyz_Py_Px_S_bb;
  Double I_ERI_Fx2z_D2y_Px_S_bb = I_ERI_Gxy2z_Py_Px_S_bb+ABY*I_ERI_Fx2z_Py_Px_S_bb;
  Double I_ERI_F3y_D2y_Px_S_bb = I_ERI_G4y_Py_Px_S_bb+ABY*I_ERI_F3y_Py_Px_S_bb;
  Double I_ERI_F2yz_D2y_Px_S_bb = I_ERI_G3yz_Py_Px_S_bb+ABY*I_ERI_F2yz_Py_Px_S_bb;
  Double I_ERI_Fy2z_D2y_Px_S_bb = I_ERI_G2y2z_Py_Px_S_bb+ABY*I_ERI_Fy2z_Py_Px_S_bb;
  Double I_ERI_F3z_D2y_Px_S_bb = I_ERI_Gy3z_Py_Px_S_bb+ABY*I_ERI_F3z_Py_Px_S_bb;
  Double I_ERI_F3x_D2z_Px_S_bb = I_ERI_G3xz_Pz_Px_S_bb+ABZ*I_ERI_F3x_Pz_Px_S_bb;
  Double I_ERI_F2xy_D2z_Px_S_bb = I_ERI_G2xyz_Pz_Px_S_bb+ABZ*I_ERI_F2xy_Pz_Px_S_bb;
  Double I_ERI_F2xz_D2z_Px_S_bb = I_ERI_G2x2z_Pz_Px_S_bb+ABZ*I_ERI_F2xz_Pz_Px_S_bb;
  Double I_ERI_Fx2y_D2z_Px_S_bb = I_ERI_Gx2yz_Pz_Px_S_bb+ABZ*I_ERI_Fx2y_Pz_Px_S_bb;
  Double I_ERI_Fxyz_D2z_Px_S_bb = I_ERI_Gxy2z_Pz_Px_S_bb+ABZ*I_ERI_Fxyz_Pz_Px_S_bb;
  Double I_ERI_Fx2z_D2z_Px_S_bb = I_ERI_Gx3z_Pz_Px_S_bb+ABZ*I_ERI_Fx2z_Pz_Px_S_bb;
  Double I_ERI_F3y_D2z_Px_S_bb = I_ERI_G3yz_Pz_Px_S_bb+ABZ*I_ERI_F3y_Pz_Px_S_bb;
  Double I_ERI_F2yz_D2z_Px_S_bb = I_ERI_G2y2z_Pz_Px_S_bb+ABZ*I_ERI_F2yz_Pz_Px_S_bb;
  Double I_ERI_Fy2z_D2z_Px_S_bb = I_ERI_Gy3z_Pz_Px_S_bb+ABZ*I_ERI_Fy2z_Pz_Px_S_bb;
  Double I_ERI_F3z_D2z_Px_S_bb = I_ERI_G4z_Pz_Px_S_bb+ABZ*I_ERI_F3z_Pz_Px_S_bb;
  Double I_ERI_F3x_D2x_Py_S_bb = I_ERI_G4x_Px_Py_S_bb+ABX*I_ERI_F3x_Px_Py_S_bb;
  Double I_ERI_F2xy_D2x_Py_S_bb = I_ERI_G3xy_Px_Py_S_bb+ABX*I_ERI_F2xy_Px_Py_S_bb;
  Double I_ERI_F2xz_D2x_Py_S_bb = I_ERI_G3xz_Px_Py_S_bb+ABX*I_ERI_F2xz_Px_Py_S_bb;
  Double I_ERI_Fx2y_D2x_Py_S_bb = I_ERI_G2x2y_Px_Py_S_bb+ABX*I_ERI_Fx2y_Px_Py_S_bb;
  Double I_ERI_Fxyz_D2x_Py_S_bb = I_ERI_G2xyz_Px_Py_S_bb+ABX*I_ERI_Fxyz_Px_Py_S_bb;
  Double I_ERI_Fx2z_D2x_Py_S_bb = I_ERI_G2x2z_Px_Py_S_bb+ABX*I_ERI_Fx2z_Px_Py_S_bb;
  Double I_ERI_F3y_D2x_Py_S_bb = I_ERI_Gx3y_Px_Py_S_bb+ABX*I_ERI_F3y_Px_Py_S_bb;
  Double I_ERI_F2yz_D2x_Py_S_bb = I_ERI_Gx2yz_Px_Py_S_bb+ABX*I_ERI_F2yz_Px_Py_S_bb;
  Double I_ERI_Fy2z_D2x_Py_S_bb = I_ERI_Gxy2z_Px_Py_S_bb+ABX*I_ERI_Fy2z_Px_Py_S_bb;
  Double I_ERI_F3z_D2x_Py_S_bb = I_ERI_Gx3z_Px_Py_S_bb+ABX*I_ERI_F3z_Px_Py_S_bb;
  Double I_ERI_F2xz_Dxy_Py_S_bb = I_ERI_G2xyz_Px_Py_S_bb+ABY*I_ERI_F2xz_Px_Py_S_bb;
  Double I_ERI_Fxyz_Dxy_Py_S_bb = I_ERI_Gx2yz_Px_Py_S_bb+ABY*I_ERI_Fxyz_Px_Py_S_bb;
  Double I_ERI_Fx2z_Dxy_Py_S_bb = I_ERI_Gxy2z_Px_Py_S_bb+ABY*I_ERI_Fx2z_Px_Py_S_bb;
  Double I_ERI_F2yz_Dxy_Py_S_bb = I_ERI_G3yz_Px_Py_S_bb+ABY*I_ERI_F2yz_Px_Py_S_bb;
  Double I_ERI_Fy2z_Dxy_Py_S_bb = I_ERI_G2y2z_Px_Py_S_bb+ABY*I_ERI_Fy2z_Px_Py_S_bb;
  Double I_ERI_F3z_Dxy_Py_S_bb = I_ERI_Gy3z_Px_Py_S_bb+ABY*I_ERI_F3z_Px_Py_S_bb;
  Double I_ERI_F3x_D2y_Py_S_bb = I_ERI_G3xy_Py_Py_S_bb+ABY*I_ERI_F3x_Py_Py_S_bb;
  Double I_ERI_F2xy_D2y_Py_S_bb = I_ERI_G2x2y_Py_Py_S_bb+ABY*I_ERI_F2xy_Py_Py_S_bb;
  Double I_ERI_F2xz_D2y_Py_S_bb = I_ERI_G2xyz_Py_Py_S_bb+ABY*I_ERI_F2xz_Py_Py_S_bb;
  Double I_ERI_Fx2y_D2y_Py_S_bb = I_ERI_Gx3y_Py_Py_S_bb+ABY*I_ERI_Fx2y_Py_Py_S_bb;
  Double I_ERI_Fxyz_D2y_Py_S_bb = I_ERI_Gx2yz_Py_Py_S_bb+ABY*I_ERI_Fxyz_Py_Py_S_bb;
  Double I_ERI_Fx2z_D2y_Py_S_bb = I_ERI_Gxy2z_Py_Py_S_bb+ABY*I_ERI_Fx2z_Py_Py_S_bb;
  Double I_ERI_F3y_D2y_Py_S_bb = I_ERI_G4y_Py_Py_S_bb+ABY*I_ERI_F3y_Py_Py_S_bb;
  Double I_ERI_F2yz_D2y_Py_S_bb = I_ERI_G3yz_Py_Py_S_bb+ABY*I_ERI_F2yz_Py_Py_S_bb;
  Double I_ERI_Fy2z_D2y_Py_S_bb = I_ERI_G2y2z_Py_Py_S_bb+ABY*I_ERI_Fy2z_Py_Py_S_bb;
  Double I_ERI_F3z_D2y_Py_S_bb = I_ERI_Gy3z_Py_Py_S_bb+ABY*I_ERI_F3z_Py_Py_S_bb;
  Double I_ERI_F3x_D2z_Py_S_bb = I_ERI_G3xz_Pz_Py_S_bb+ABZ*I_ERI_F3x_Pz_Py_S_bb;
  Double I_ERI_F2xy_D2z_Py_S_bb = I_ERI_G2xyz_Pz_Py_S_bb+ABZ*I_ERI_F2xy_Pz_Py_S_bb;
  Double I_ERI_F2xz_D2z_Py_S_bb = I_ERI_G2x2z_Pz_Py_S_bb+ABZ*I_ERI_F2xz_Pz_Py_S_bb;
  Double I_ERI_Fx2y_D2z_Py_S_bb = I_ERI_Gx2yz_Pz_Py_S_bb+ABZ*I_ERI_Fx2y_Pz_Py_S_bb;
  Double I_ERI_Fxyz_D2z_Py_S_bb = I_ERI_Gxy2z_Pz_Py_S_bb+ABZ*I_ERI_Fxyz_Pz_Py_S_bb;
  Double I_ERI_Fx2z_D2z_Py_S_bb = I_ERI_Gx3z_Pz_Py_S_bb+ABZ*I_ERI_Fx2z_Pz_Py_S_bb;
  Double I_ERI_F3y_D2z_Py_S_bb = I_ERI_G3yz_Pz_Py_S_bb+ABZ*I_ERI_F3y_Pz_Py_S_bb;
  Double I_ERI_F2yz_D2z_Py_S_bb = I_ERI_G2y2z_Pz_Py_S_bb+ABZ*I_ERI_F2yz_Pz_Py_S_bb;
  Double I_ERI_Fy2z_D2z_Py_S_bb = I_ERI_Gy3z_Pz_Py_S_bb+ABZ*I_ERI_Fy2z_Pz_Py_S_bb;
  Double I_ERI_F3z_D2z_Py_S_bb = I_ERI_G4z_Pz_Py_S_bb+ABZ*I_ERI_F3z_Pz_Py_S_bb;
  Double I_ERI_F3x_D2x_Pz_S_bb = I_ERI_G4x_Px_Pz_S_bb+ABX*I_ERI_F3x_Px_Pz_S_bb;
  Double I_ERI_F2xy_D2x_Pz_S_bb = I_ERI_G3xy_Px_Pz_S_bb+ABX*I_ERI_F2xy_Px_Pz_S_bb;
  Double I_ERI_F2xz_D2x_Pz_S_bb = I_ERI_G3xz_Px_Pz_S_bb+ABX*I_ERI_F2xz_Px_Pz_S_bb;
  Double I_ERI_Fx2y_D2x_Pz_S_bb = I_ERI_G2x2y_Px_Pz_S_bb+ABX*I_ERI_Fx2y_Px_Pz_S_bb;
  Double I_ERI_Fxyz_D2x_Pz_S_bb = I_ERI_G2xyz_Px_Pz_S_bb+ABX*I_ERI_Fxyz_Px_Pz_S_bb;
  Double I_ERI_Fx2z_D2x_Pz_S_bb = I_ERI_G2x2z_Px_Pz_S_bb+ABX*I_ERI_Fx2z_Px_Pz_S_bb;
  Double I_ERI_F3y_D2x_Pz_S_bb = I_ERI_Gx3y_Px_Pz_S_bb+ABX*I_ERI_F3y_Px_Pz_S_bb;
  Double I_ERI_F2yz_D2x_Pz_S_bb = I_ERI_Gx2yz_Px_Pz_S_bb+ABX*I_ERI_F2yz_Px_Pz_S_bb;
  Double I_ERI_Fy2z_D2x_Pz_S_bb = I_ERI_Gxy2z_Px_Pz_S_bb+ABX*I_ERI_Fy2z_Px_Pz_S_bb;
  Double I_ERI_F3z_D2x_Pz_S_bb = I_ERI_Gx3z_Px_Pz_S_bb+ABX*I_ERI_F3z_Px_Pz_S_bb;
  Double I_ERI_F2xz_Dxy_Pz_S_bb = I_ERI_G2xyz_Px_Pz_S_bb+ABY*I_ERI_F2xz_Px_Pz_S_bb;
  Double I_ERI_Fxyz_Dxy_Pz_S_bb = I_ERI_Gx2yz_Px_Pz_S_bb+ABY*I_ERI_Fxyz_Px_Pz_S_bb;
  Double I_ERI_Fx2z_Dxy_Pz_S_bb = I_ERI_Gxy2z_Px_Pz_S_bb+ABY*I_ERI_Fx2z_Px_Pz_S_bb;
  Double I_ERI_F2yz_Dxy_Pz_S_bb = I_ERI_G3yz_Px_Pz_S_bb+ABY*I_ERI_F2yz_Px_Pz_S_bb;
  Double I_ERI_Fy2z_Dxy_Pz_S_bb = I_ERI_G2y2z_Px_Pz_S_bb+ABY*I_ERI_Fy2z_Px_Pz_S_bb;
  Double I_ERI_F3z_Dxy_Pz_S_bb = I_ERI_Gy3z_Px_Pz_S_bb+ABY*I_ERI_F3z_Px_Pz_S_bb;
  Double I_ERI_F3x_D2y_Pz_S_bb = I_ERI_G3xy_Py_Pz_S_bb+ABY*I_ERI_F3x_Py_Pz_S_bb;
  Double I_ERI_F2xy_D2y_Pz_S_bb = I_ERI_G2x2y_Py_Pz_S_bb+ABY*I_ERI_F2xy_Py_Pz_S_bb;
  Double I_ERI_F2xz_D2y_Pz_S_bb = I_ERI_G2xyz_Py_Pz_S_bb+ABY*I_ERI_F2xz_Py_Pz_S_bb;
  Double I_ERI_Fx2y_D2y_Pz_S_bb = I_ERI_Gx3y_Py_Pz_S_bb+ABY*I_ERI_Fx2y_Py_Pz_S_bb;
  Double I_ERI_Fxyz_D2y_Pz_S_bb = I_ERI_Gx2yz_Py_Pz_S_bb+ABY*I_ERI_Fxyz_Py_Pz_S_bb;
  Double I_ERI_Fx2z_D2y_Pz_S_bb = I_ERI_Gxy2z_Py_Pz_S_bb+ABY*I_ERI_Fx2z_Py_Pz_S_bb;
  Double I_ERI_F3y_D2y_Pz_S_bb = I_ERI_G4y_Py_Pz_S_bb+ABY*I_ERI_F3y_Py_Pz_S_bb;
  Double I_ERI_F2yz_D2y_Pz_S_bb = I_ERI_G3yz_Py_Pz_S_bb+ABY*I_ERI_F2yz_Py_Pz_S_bb;
  Double I_ERI_Fy2z_D2y_Pz_S_bb = I_ERI_G2y2z_Py_Pz_S_bb+ABY*I_ERI_Fy2z_Py_Pz_S_bb;
  Double I_ERI_F3z_D2y_Pz_S_bb = I_ERI_Gy3z_Py_Pz_S_bb+ABY*I_ERI_F3z_Py_Pz_S_bb;
  Double I_ERI_F3x_D2z_Pz_S_bb = I_ERI_G3xz_Pz_Pz_S_bb+ABZ*I_ERI_F3x_Pz_Pz_S_bb;
  Double I_ERI_F2xy_D2z_Pz_S_bb = I_ERI_G2xyz_Pz_Pz_S_bb+ABZ*I_ERI_F2xy_Pz_Pz_S_bb;
  Double I_ERI_F2xz_D2z_Pz_S_bb = I_ERI_G2x2z_Pz_Pz_S_bb+ABZ*I_ERI_F2xz_Pz_Pz_S_bb;
  Double I_ERI_Fx2y_D2z_Pz_S_bb = I_ERI_Gx2yz_Pz_Pz_S_bb+ABZ*I_ERI_Fx2y_Pz_Pz_S_bb;
  Double I_ERI_Fxyz_D2z_Pz_S_bb = I_ERI_Gxy2z_Pz_Pz_S_bb+ABZ*I_ERI_Fxyz_Pz_Pz_S_bb;
  Double I_ERI_Fx2z_D2z_Pz_S_bb = I_ERI_Gx3z_Pz_Pz_S_bb+ABZ*I_ERI_Fx2z_Pz_Pz_S_bb;
  Double I_ERI_F3y_D2z_Pz_S_bb = I_ERI_G3yz_Pz_Pz_S_bb+ABZ*I_ERI_F3y_Pz_Pz_S_bb;
  Double I_ERI_F2yz_D2z_Pz_S_bb = I_ERI_G2y2z_Pz_Pz_S_bb+ABZ*I_ERI_F2yz_Pz_Pz_S_bb;
  Double I_ERI_Fy2z_D2z_Pz_S_bb = I_ERI_Gy3z_Pz_Pz_S_bb+ABZ*I_ERI_Fy2z_Pz_Pz_S_bb;
  Double I_ERI_F3z_D2z_Pz_S_bb = I_ERI_G4z_Pz_Pz_S_bb+ABZ*I_ERI_F3z_Pz_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_F_P_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_bb
   * RHS shell quartet name: SQ_ERI_D_D_P_S_bb
   ************************************************************/
  Double I_ERI_D2x_F3x_Px_S_bb = I_ERI_F3x_D2x_Px_S_bb+ABX*I_ERI_D2x_D2x_Px_S_bb;
  Double I_ERI_Dxy_F3x_Px_S_bb = I_ERI_F2xy_D2x_Px_S_bb+ABX*I_ERI_Dxy_D2x_Px_S_bb;
  Double I_ERI_Dxz_F3x_Px_S_bb = I_ERI_F2xz_D2x_Px_S_bb+ABX*I_ERI_Dxz_D2x_Px_S_bb;
  Double I_ERI_D2y_F3x_Px_S_bb = I_ERI_Fx2y_D2x_Px_S_bb+ABX*I_ERI_D2y_D2x_Px_S_bb;
  Double I_ERI_Dyz_F3x_Px_S_bb = I_ERI_Fxyz_D2x_Px_S_bb+ABX*I_ERI_Dyz_D2x_Px_S_bb;
  Double I_ERI_D2z_F3x_Px_S_bb = I_ERI_Fx2z_D2x_Px_S_bb+ABX*I_ERI_D2z_D2x_Px_S_bb;
  Double I_ERI_D2x_F2xy_Px_S_bb = I_ERI_F2xy_D2x_Px_S_bb+ABY*I_ERI_D2x_D2x_Px_S_bb;
  Double I_ERI_Dxy_F2xy_Px_S_bb = I_ERI_Fx2y_D2x_Px_S_bb+ABY*I_ERI_Dxy_D2x_Px_S_bb;
  Double I_ERI_Dxz_F2xy_Px_S_bb = I_ERI_Fxyz_D2x_Px_S_bb+ABY*I_ERI_Dxz_D2x_Px_S_bb;
  Double I_ERI_D2y_F2xy_Px_S_bb = I_ERI_F3y_D2x_Px_S_bb+ABY*I_ERI_D2y_D2x_Px_S_bb;
  Double I_ERI_Dyz_F2xy_Px_S_bb = I_ERI_F2yz_D2x_Px_S_bb+ABY*I_ERI_Dyz_D2x_Px_S_bb;
  Double I_ERI_D2z_F2xy_Px_S_bb = I_ERI_Fy2z_D2x_Px_S_bb+ABY*I_ERI_D2z_D2x_Px_S_bb;
  Double I_ERI_D2x_F2xz_Px_S_bb = I_ERI_F2xz_D2x_Px_S_bb+ABZ*I_ERI_D2x_D2x_Px_S_bb;
  Double I_ERI_Dxy_F2xz_Px_S_bb = I_ERI_Fxyz_D2x_Px_S_bb+ABZ*I_ERI_Dxy_D2x_Px_S_bb;
  Double I_ERI_Dxz_F2xz_Px_S_bb = I_ERI_Fx2z_D2x_Px_S_bb+ABZ*I_ERI_Dxz_D2x_Px_S_bb;
  Double I_ERI_D2y_F2xz_Px_S_bb = I_ERI_F2yz_D2x_Px_S_bb+ABZ*I_ERI_D2y_D2x_Px_S_bb;
  Double I_ERI_Dyz_F2xz_Px_S_bb = I_ERI_Fy2z_D2x_Px_S_bb+ABZ*I_ERI_Dyz_D2x_Px_S_bb;
  Double I_ERI_D2z_F2xz_Px_S_bb = I_ERI_F3z_D2x_Px_S_bb+ABZ*I_ERI_D2z_D2x_Px_S_bb;
  Double I_ERI_D2x_Fx2y_Px_S_bb = I_ERI_F3x_D2y_Px_S_bb+ABX*I_ERI_D2x_D2y_Px_S_bb;
  Double I_ERI_Dxy_Fx2y_Px_S_bb = I_ERI_F2xy_D2y_Px_S_bb+ABX*I_ERI_Dxy_D2y_Px_S_bb;
  Double I_ERI_Dxz_Fx2y_Px_S_bb = I_ERI_F2xz_D2y_Px_S_bb+ABX*I_ERI_Dxz_D2y_Px_S_bb;
  Double I_ERI_D2y_Fx2y_Px_S_bb = I_ERI_Fx2y_D2y_Px_S_bb+ABX*I_ERI_D2y_D2y_Px_S_bb;
  Double I_ERI_Dyz_Fx2y_Px_S_bb = I_ERI_Fxyz_D2y_Px_S_bb+ABX*I_ERI_Dyz_D2y_Px_S_bb;
  Double I_ERI_D2z_Fx2y_Px_S_bb = I_ERI_Fx2z_D2y_Px_S_bb+ABX*I_ERI_D2z_D2y_Px_S_bb;
  Double I_ERI_D2x_Fxyz_Px_S_bb = I_ERI_F2xz_Dxy_Px_S_bb+ABZ*I_ERI_D2x_Dxy_Px_S_bb;
  Double I_ERI_Dxy_Fxyz_Px_S_bb = I_ERI_Fxyz_Dxy_Px_S_bb+ABZ*I_ERI_Dxy_Dxy_Px_S_bb;
  Double I_ERI_Dxz_Fxyz_Px_S_bb = I_ERI_Fx2z_Dxy_Px_S_bb+ABZ*I_ERI_Dxz_Dxy_Px_S_bb;
  Double I_ERI_D2y_Fxyz_Px_S_bb = I_ERI_F2yz_Dxy_Px_S_bb+ABZ*I_ERI_D2y_Dxy_Px_S_bb;
  Double I_ERI_Dyz_Fxyz_Px_S_bb = I_ERI_Fy2z_Dxy_Px_S_bb+ABZ*I_ERI_Dyz_Dxy_Px_S_bb;
  Double I_ERI_D2z_Fxyz_Px_S_bb = I_ERI_F3z_Dxy_Px_S_bb+ABZ*I_ERI_D2z_Dxy_Px_S_bb;
  Double I_ERI_D2x_Fx2z_Px_S_bb = I_ERI_F3x_D2z_Px_S_bb+ABX*I_ERI_D2x_D2z_Px_S_bb;
  Double I_ERI_Dxy_Fx2z_Px_S_bb = I_ERI_F2xy_D2z_Px_S_bb+ABX*I_ERI_Dxy_D2z_Px_S_bb;
  Double I_ERI_Dxz_Fx2z_Px_S_bb = I_ERI_F2xz_D2z_Px_S_bb+ABX*I_ERI_Dxz_D2z_Px_S_bb;
  Double I_ERI_D2y_Fx2z_Px_S_bb = I_ERI_Fx2y_D2z_Px_S_bb+ABX*I_ERI_D2y_D2z_Px_S_bb;
  Double I_ERI_Dyz_Fx2z_Px_S_bb = I_ERI_Fxyz_D2z_Px_S_bb+ABX*I_ERI_Dyz_D2z_Px_S_bb;
  Double I_ERI_D2z_Fx2z_Px_S_bb = I_ERI_Fx2z_D2z_Px_S_bb+ABX*I_ERI_D2z_D2z_Px_S_bb;
  Double I_ERI_D2x_F3y_Px_S_bb = I_ERI_F2xy_D2y_Px_S_bb+ABY*I_ERI_D2x_D2y_Px_S_bb;
  Double I_ERI_Dxy_F3y_Px_S_bb = I_ERI_Fx2y_D2y_Px_S_bb+ABY*I_ERI_Dxy_D2y_Px_S_bb;
  Double I_ERI_Dxz_F3y_Px_S_bb = I_ERI_Fxyz_D2y_Px_S_bb+ABY*I_ERI_Dxz_D2y_Px_S_bb;
  Double I_ERI_D2y_F3y_Px_S_bb = I_ERI_F3y_D2y_Px_S_bb+ABY*I_ERI_D2y_D2y_Px_S_bb;
  Double I_ERI_Dyz_F3y_Px_S_bb = I_ERI_F2yz_D2y_Px_S_bb+ABY*I_ERI_Dyz_D2y_Px_S_bb;
  Double I_ERI_D2z_F3y_Px_S_bb = I_ERI_Fy2z_D2y_Px_S_bb+ABY*I_ERI_D2z_D2y_Px_S_bb;
  Double I_ERI_D2x_F2yz_Px_S_bb = I_ERI_F2xz_D2y_Px_S_bb+ABZ*I_ERI_D2x_D2y_Px_S_bb;
  Double I_ERI_Dxy_F2yz_Px_S_bb = I_ERI_Fxyz_D2y_Px_S_bb+ABZ*I_ERI_Dxy_D2y_Px_S_bb;
  Double I_ERI_Dxz_F2yz_Px_S_bb = I_ERI_Fx2z_D2y_Px_S_bb+ABZ*I_ERI_Dxz_D2y_Px_S_bb;
  Double I_ERI_D2y_F2yz_Px_S_bb = I_ERI_F2yz_D2y_Px_S_bb+ABZ*I_ERI_D2y_D2y_Px_S_bb;
  Double I_ERI_Dyz_F2yz_Px_S_bb = I_ERI_Fy2z_D2y_Px_S_bb+ABZ*I_ERI_Dyz_D2y_Px_S_bb;
  Double I_ERI_D2z_F2yz_Px_S_bb = I_ERI_F3z_D2y_Px_S_bb+ABZ*I_ERI_D2z_D2y_Px_S_bb;
  Double I_ERI_D2x_Fy2z_Px_S_bb = I_ERI_F2xy_D2z_Px_S_bb+ABY*I_ERI_D2x_D2z_Px_S_bb;
  Double I_ERI_Dxy_Fy2z_Px_S_bb = I_ERI_Fx2y_D2z_Px_S_bb+ABY*I_ERI_Dxy_D2z_Px_S_bb;
  Double I_ERI_Dxz_Fy2z_Px_S_bb = I_ERI_Fxyz_D2z_Px_S_bb+ABY*I_ERI_Dxz_D2z_Px_S_bb;
  Double I_ERI_D2y_Fy2z_Px_S_bb = I_ERI_F3y_D2z_Px_S_bb+ABY*I_ERI_D2y_D2z_Px_S_bb;
  Double I_ERI_Dyz_Fy2z_Px_S_bb = I_ERI_F2yz_D2z_Px_S_bb+ABY*I_ERI_Dyz_D2z_Px_S_bb;
  Double I_ERI_D2z_Fy2z_Px_S_bb = I_ERI_Fy2z_D2z_Px_S_bb+ABY*I_ERI_D2z_D2z_Px_S_bb;
  Double I_ERI_D2x_F3z_Px_S_bb = I_ERI_F2xz_D2z_Px_S_bb+ABZ*I_ERI_D2x_D2z_Px_S_bb;
  Double I_ERI_Dxy_F3z_Px_S_bb = I_ERI_Fxyz_D2z_Px_S_bb+ABZ*I_ERI_Dxy_D2z_Px_S_bb;
  Double I_ERI_Dxz_F3z_Px_S_bb = I_ERI_Fx2z_D2z_Px_S_bb+ABZ*I_ERI_Dxz_D2z_Px_S_bb;
  Double I_ERI_D2y_F3z_Px_S_bb = I_ERI_F2yz_D2z_Px_S_bb+ABZ*I_ERI_D2y_D2z_Px_S_bb;
  Double I_ERI_Dyz_F3z_Px_S_bb = I_ERI_Fy2z_D2z_Px_S_bb+ABZ*I_ERI_Dyz_D2z_Px_S_bb;
  Double I_ERI_D2z_F3z_Px_S_bb = I_ERI_F3z_D2z_Px_S_bb+ABZ*I_ERI_D2z_D2z_Px_S_bb;
  Double I_ERI_D2x_F3x_Py_S_bb = I_ERI_F3x_D2x_Py_S_bb+ABX*I_ERI_D2x_D2x_Py_S_bb;
  Double I_ERI_Dxy_F3x_Py_S_bb = I_ERI_F2xy_D2x_Py_S_bb+ABX*I_ERI_Dxy_D2x_Py_S_bb;
  Double I_ERI_Dxz_F3x_Py_S_bb = I_ERI_F2xz_D2x_Py_S_bb+ABX*I_ERI_Dxz_D2x_Py_S_bb;
  Double I_ERI_D2y_F3x_Py_S_bb = I_ERI_Fx2y_D2x_Py_S_bb+ABX*I_ERI_D2y_D2x_Py_S_bb;
  Double I_ERI_Dyz_F3x_Py_S_bb = I_ERI_Fxyz_D2x_Py_S_bb+ABX*I_ERI_Dyz_D2x_Py_S_bb;
  Double I_ERI_D2z_F3x_Py_S_bb = I_ERI_Fx2z_D2x_Py_S_bb+ABX*I_ERI_D2z_D2x_Py_S_bb;
  Double I_ERI_D2x_F2xy_Py_S_bb = I_ERI_F2xy_D2x_Py_S_bb+ABY*I_ERI_D2x_D2x_Py_S_bb;
  Double I_ERI_Dxy_F2xy_Py_S_bb = I_ERI_Fx2y_D2x_Py_S_bb+ABY*I_ERI_Dxy_D2x_Py_S_bb;
  Double I_ERI_Dxz_F2xy_Py_S_bb = I_ERI_Fxyz_D2x_Py_S_bb+ABY*I_ERI_Dxz_D2x_Py_S_bb;
  Double I_ERI_D2y_F2xy_Py_S_bb = I_ERI_F3y_D2x_Py_S_bb+ABY*I_ERI_D2y_D2x_Py_S_bb;
  Double I_ERI_Dyz_F2xy_Py_S_bb = I_ERI_F2yz_D2x_Py_S_bb+ABY*I_ERI_Dyz_D2x_Py_S_bb;
  Double I_ERI_D2z_F2xy_Py_S_bb = I_ERI_Fy2z_D2x_Py_S_bb+ABY*I_ERI_D2z_D2x_Py_S_bb;
  Double I_ERI_D2x_F2xz_Py_S_bb = I_ERI_F2xz_D2x_Py_S_bb+ABZ*I_ERI_D2x_D2x_Py_S_bb;
  Double I_ERI_Dxy_F2xz_Py_S_bb = I_ERI_Fxyz_D2x_Py_S_bb+ABZ*I_ERI_Dxy_D2x_Py_S_bb;
  Double I_ERI_Dxz_F2xz_Py_S_bb = I_ERI_Fx2z_D2x_Py_S_bb+ABZ*I_ERI_Dxz_D2x_Py_S_bb;
  Double I_ERI_D2y_F2xz_Py_S_bb = I_ERI_F2yz_D2x_Py_S_bb+ABZ*I_ERI_D2y_D2x_Py_S_bb;
  Double I_ERI_Dyz_F2xz_Py_S_bb = I_ERI_Fy2z_D2x_Py_S_bb+ABZ*I_ERI_Dyz_D2x_Py_S_bb;
  Double I_ERI_D2z_F2xz_Py_S_bb = I_ERI_F3z_D2x_Py_S_bb+ABZ*I_ERI_D2z_D2x_Py_S_bb;
  Double I_ERI_D2x_Fx2y_Py_S_bb = I_ERI_F3x_D2y_Py_S_bb+ABX*I_ERI_D2x_D2y_Py_S_bb;
  Double I_ERI_Dxy_Fx2y_Py_S_bb = I_ERI_F2xy_D2y_Py_S_bb+ABX*I_ERI_Dxy_D2y_Py_S_bb;
  Double I_ERI_Dxz_Fx2y_Py_S_bb = I_ERI_F2xz_D2y_Py_S_bb+ABX*I_ERI_Dxz_D2y_Py_S_bb;
  Double I_ERI_D2y_Fx2y_Py_S_bb = I_ERI_Fx2y_D2y_Py_S_bb+ABX*I_ERI_D2y_D2y_Py_S_bb;
  Double I_ERI_Dyz_Fx2y_Py_S_bb = I_ERI_Fxyz_D2y_Py_S_bb+ABX*I_ERI_Dyz_D2y_Py_S_bb;
  Double I_ERI_D2z_Fx2y_Py_S_bb = I_ERI_Fx2z_D2y_Py_S_bb+ABX*I_ERI_D2z_D2y_Py_S_bb;
  Double I_ERI_D2x_Fxyz_Py_S_bb = I_ERI_F2xz_Dxy_Py_S_bb+ABZ*I_ERI_D2x_Dxy_Py_S_bb;
  Double I_ERI_Dxy_Fxyz_Py_S_bb = I_ERI_Fxyz_Dxy_Py_S_bb+ABZ*I_ERI_Dxy_Dxy_Py_S_bb;
  Double I_ERI_Dxz_Fxyz_Py_S_bb = I_ERI_Fx2z_Dxy_Py_S_bb+ABZ*I_ERI_Dxz_Dxy_Py_S_bb;
  Double I_ERI_D2y_Fxyz_Py_S_bb = I_ERI_F2yz_Dxy_Py_S_bb+ABZ*I_ERI_D2y_Dxy_Py_S_bb;
  Double I_ERI_Dyz_Fxyz_Py_S_bb = I_ERI_Fy2z_Dxy_Py_S_bb+ABZ*I_ERI_Dyz_Dxy_Py_S_bb;
  Double I_ERI_D2z_Fxyz_Py_S_bb = I_ERI_F3z_Dxy_Py_S_bb+ABZ*I_ERI_D2z_Dxy_Py_S_bb;
  Double I_ERI_D2x_Fx2z_Py_S_bb = I_ERI_F3x_D2z_Py_S_bb+ABX*I_ERI_D2x_D2z_Py_S_bb;
  Double I_ERI_Dxy_Fx2z_Py_S_bb = I_ERI_F2xy_D2z_Py_S_bb+ABX*I_ERI_Dxy_D2z_Py_S_bb;
  Double I_ERI_Dxz_Fx2z_Py_S_bb = I_ERI_F2xz_D2z_Py_S_bb+ABX*I_ERI_Dxz_D2z_Py_S_bb;
  Double I_ERI_D2y_Fx2z_Py_S_bb = I_ERI_Fx2y_D2z_Py_S_bb+ABX*I_ERI_D2y_D2z_Py_S_bb;
  Double I_ERI_Dyz_Fx2z_Py_S_bb = I_ERI_Fxyz_D2z_Py_S_bb+ABX*I_ERI_Dyz_D2z_Py_S_bb;
  Double I_ERI_D2z_Fx2z_Py_S_bb = I_ERI_Fx2z_D2z_Py_S_bb+ABX*I_ERI_D2z_D2z_Py_S_bb;
  Double I_ERI_D2x_F3y_Py_S_bb = I_ERI_F2xy_D2y_Py_S_bb+ABY*I_ERI_D2x_D2y_Py_S_bb;
  Double I_ERI_Dxy_F3y_Py_S_bb = I_ERI_Fx2y_D2y_Py_S_bb+ABY*I_ERI_Dxy_D2y_Py_S_bb;
  Double I_ERI_Dxz_F3y_Py_S_bb = I_ERI_Fxyz_D2y_Py_S_bb+ABY*I_ERI_Dxz_D2y_Py_S_bb;
  Double I_ERI_D2y_F3y_Py_S_bb = I_ERI_F3y_D2y_Py_S_bb+ABY*I_ERI_D2y_D2y_Py_S_bb;
  Double I_ERI_Dyz_F3y_Py_S_bb = I_ERI_F2yz_D2y_Py_S_bb+ABY*I_ERI_Dyz_D2y_Py_S_bb;
  Double I_ERI_D2z_F3y_Py_S_bb = I_ERI_Fy2z_D2y_Py_S_bb+ABY*I_ERI_D2z_D2y_Py_S_bb;
  Double I_ERI_D2x_F2yz_Py_S_bb = I_ERI_F2xz_D2y_Py_S_bb+ABZ*I_ERI_D2x_D2y_Py_S_bb;
  Double I_ERI_Dxy_F2yz_Py_S_bb = I_ERI_Fxyz_D2y_Py_S_bb+ABZ*I_ERI_Dxy_D2y_Py_S_bb;
  Double I_ERI_Dxz_F2yz_Py_S_bb = I_ERI_Fx2z_D2y_Py_S_bb+ABZ*I_ERI_Dxz_D2y_Py_S_bb;
  Double I_ERI_D2y_F2yz_Py_S_bb = I_ERI_F2yz_D2y_Py_S_bb+ABZ*I_ERI_D2y_D2y_Py_S_bb;
  Double I_ERI_Dyz_F2yz_Py_S_bb = I_ERI_Fy2z_D2y_Py_S_bb+ABZ*I_ERI_Dyz_D2y_Py_S_bb;
  Double I_ERI_D2z_F2yz_Py_S_bb = I_ERI_F3z_D2y_Py_S_bb+ABZ*I_ERI_D2z_D2y_Py_S_bb;
  Double I_ERI_D2x_Fy2z_Py_S_bb = I_ERI_F2xy_D2z_Py_S_bb+ABY*I_ERI_D2x_D2z_Py_S_bb;
  Double I_ERI_Dxy_Fy2z_Py_S_bb = I_ERI_Fx2y_D2z_Py_S_bb+ABY*I_ERI_Dxy_D2z_Py_S_bb;
  Double I_ERI_Dxz_Fy2z_Py_S_bb = I_ERI_Fxyz_D2z_Py_S_bb+ABY*I_ERI_Dxz_D2z_Py_S_bb;
  Double I_ERI_D2y_Fy2z_Py_S_bb = I_ERI_F3y_D2z_Py_S_bb+ABY*I_ERI_D2y_D2z_Py_S_bb;
  Double I_ERI_Dyz_Fy2z_Py_S_bb = I_ERI_F2yz_D2z_Py_S_bb+ABY*I_ERI_Dyz_D2z_Py_S_bb;
  Double I_ERI_D2z_Fy2z_Py_S_bb = I_ERI_Fy2z_D2z_Py_S_bb+ABY*I_ERI_D2z_D2z_Py_S_bb;
  Double I_ERI_D2x_F3z_Py_S_bb = I_ERI_F2xz_D2z_Py_S_bb+ABZ*I_ERI_D2x_D2z_Py_S_bb;
  Double I_ERI_Dxy_F3z_Py_S_bb = I_ERI_Fxyz_D2z_Py_S_bb+ABZ*I_ERI_Dxy_D2z_Py_S_bb;
  Double I_ERI_Dxz_F3z_Py_S_bb = I_ERI_Fx2z_D2z_Py_S_bb+ABZ*I_ERI_Dxz_D2z_Py_S_bb;
  Double I_ERI_D2y_F3z_Py_S_bb = I_ERI_F2yz_D2z_Py_S_bb+ABZ*I_ERI_D2y_D2z_Py_S_bb;
  Double I_ERI_Dyz_F3z_Py_S_bb = I_ERI_Fy2z_D2z_Py_S_bb+ABZ*I_ERI_Dyz_D2z_Py_S_bb;
  Double I_ERI_D2z_F3z_Py_S_bb = I_ERI_F3z_D2z_Py_S_bb+ABZ*I_ERI_D2z_D2z_Py_S_bb;
  Double I_ERI_D2x_F3x_Pz_S_bb = I_ERI_F3x_D2x_Pz_S_bb+ABX*I_ERI_D2x_D2x_Pz_S_bb;
  Double I_ERI_Dxy_F3x_Pz_S_bb = I_ERI_F2xy_D2x_Pz_S_bb+ABX*I_ERI_Dxy_D2x_Pz_S_bb;
  Double I_ERI_Dxz_F3x_Pz_S_bb = I_ERI_F2xz_D2x_Pz_S_bb+ABX*I_ERI_Dxz_D2x_Pz_S_bb;
  Double I_ERI_D2y_F3x_Pz_S_bb = I_ERI_Fx2y_D2x_Pz_S_bb+ABX*I_ERI_D2y_D2x_Pz_S_bb;
  Double I_ERI_Dyz_F3x_Pz_S_bb = I_ERI_Fxyz_D2x_Pz_S_bb+ABX*I_ERI_Dyz_D2x_Pz_S_bb;
  Double I_ERI_D2z_F3x_Pz_S_bb = I_ERI_Fx2z_D2x_Pz_S_bb+ABX*I_ERI_D2z_D2x_Pz_S_bb;
  Double I_ERI_D2x_F2xy_Pz_S_bb = I_ERI_F2xy_D2x_Pz_S_bb+ABY*I_ERI_D2x_D2x_Pz_S_bb;
  Double I_ERI_Dxy_F2xy_Pz_S_bb = I_ERI_Fx2y_D2x_Pz_S_bb+ABY*I_ERI_Dxy_D2x_Pz_S_bb;
  Double I_ERI_Dxz_F2xy_Pz_S_bb = I_ERI_Fxyz_D2x_Pz_S_bb+ABY*I_ERI_Dxz_D2x_Pz_S_bb;
  Double I_ERI_D2y_F2xy_Pz_S_bb = I_ERI_F3y_D2x_Pz_S_bb+ABY*I_ERI_D2y_D2x_Pz_S_bb;
  Double I_ERI_Dyz_F2xy_Pz_S_bb = I_ERI_F2yz_D2x_Pz_S_bb+ABY*I_ERI_Dyz_D2x_Pz_S_bb;
  Double I_ERI_D2z_F2xy_Pz_S_bb = I_ERI_Fy2z_D2x_Pz_S_bb+ABY*I_ERI_D2z_D2x_Pz_S_bb;
  Double I_ERI_D2x_F2xz_Pz_S_bb = I_ERI_F2xz_D2x_Pz_S_bb+ABZ*I_ERI_D2x_D2x_Pz_S_bb;
  Double I_ERI_Dxy_F2xz_Pz_S_bb = I_ERI_Fxyz_D2x_Pz_S_bb+ABZ*I_ERI_Dxy_D2x_Pz_S_bb;
  Double I_ERI_Dxz_F2xz_Pz_S_bb = I_ERI_Fx2z_D2x_Pz_S_bb+ABZ*I_ERI_Dxz_D2x_Pz_S_bb;
  Double I_ERI_D2y_F2xz_Pz_S_bb = I_ERI_F2yz_D2x_Pz_S_bb+ABZ*I_ERI_D2y_D2x_Pz_S_bb;
  Double I_ERI_Dyz_F2xz_Pz_S_bb = I_ERI_Fy2z_D2x_Pz_S_bb+ABZ*I_ERI_Dyz_D2x_Pz_S_bb;
  Double I_ERI_D2z_F2xz_Pz_S_bb = I_ERI_F3z_D2x_Pz_S_bb+ABZ*I_ERI_D2z_D2x_Pz_S_bb;
  Double I_ERI_D2x_Fx2y_Pz_S_bb = I_ERI_F3x_D2y_Pz_S_bb+ABX*I_ERI_D2x_D2y_Pz_S_bb;
  Double I_ERI_Dxy_Fx2y_Pz_S_bb = I_ERI_F2xy_D2y_Pz_S_bb+ABX*I_ERI_Dxy_D2y_Pz_S_bb;
  Double I_ERI_Dxz_Fx2y_Pz_S_bb = I_ERI_F2xz_D2y_Pz_S_bb+ABX*I_ERI_Dxz_D2y_Pz_S_bb;
  Double I_ERI_D2y_Fx2y_Pz_S_bb = I_ERI_Fx2y_D2y_Pz_S_bb+ABX*I_ERI_D2y_D2y_Pz_S_bb;
  Double I_ERI_Dyz_Fx2y_Pz_S_bb = I_ERI_Fxyz_D2y_Pz_S_bb+ABX*I_ERI_Dyz_D2y_Pz_S_bb;
  Double I_ERI_D2z_Fx2y_Pz_S_bb = I_ERI_Fx2z_D2y_Pz_S_bb+ABX*I_ERI_D2z_D2y_Pz_S_bb;
  Double I_ERI_D2x_Fxyz_Pz_S_bb = I_ERI_F2xz_Dxy_Pz_S_bb+ABZ*I_ERI_D2x_Dxy_Pz_S_bb;
  Double I_ERI_Dxy_Fxyz_Pz_S_bb = I_ERI_Fxyz_Dxy_Pz_S_bb+ABZ*I_ERI_Dxy_Dxy_Pz_S_bb;
  Double I_ERI_Dxz_Fxyz_Pz_S_bb = I_ERI_Fx2z_Dxy_Pz_S_bb+ABZ*I_ERI_Dxz_Dxy_Pz_S_bb;
  Double I_ERI_D2y_Fxyz_Pz_S_bb = I_ERI_F2yz_Dxy_Pz_S_bb+ABZ*I_ERI_D2y_Dxy_Pz_S_bb;
  Double I_ERI_Dyz_Fxyz_Pz_S_bb = I_ERI_Fy2z_Dxy_Pz_S_bb+ABZ*I_ERI_Dyz_Dxy_Pz_S_bb;
  Double I_ERI_D2z_Fxyz_Pz_S_bb = I_ERI_F3z_Dxy_Pz_S_bb+ABZ*I_ERI_D2z_Dxy_Pz_S_bb;
  Double I_ERI_D2x_Fx2z_Pz_S_bb = I_ERI_F3x_D2z_Pz_S_bb+ABX*I_ERI_D2x_D2z_Pz_S_bb;
  Double I_ERI_Dxy_Fx2z_Pz_S_bb = I_ERI_F2xy_D2z_Pz_S_bb+ABX*I_ERI_Dxy_D2z_Pz_S_bb;
  Double I_ERI_Dxz_Fx2z_Pz_S_bb = I_ERI_F2xz_D2z_Pz_S_bb+ABX*I_ERI_Dxz_D2z_Pz_S_bb;
  Double I_ERI_D2y_Fx2z_Pz_S_bb = I_ERI_Fx2y_D2z_Pz_S_bb+ABX*I_ERI_D2y_D2z_Pz_S_bb;
  Double I_ERI_Dyz_Fx2z_Pz_S_bb = I_ERI_Fxyz_D2z_Pz_S_bb+ABX*I_ERI_Dyz_D2z_Pz_S_bb;
  Double I_ERI_D2z_Fx2z_Pz_S_bb = I_ERI_Fx2z_D2z_Pz_S_bb+ABX*I_ERI_D2z_D2z_Pz_S_bb;
  Double I_ERI_D2x_F3y_Pz_S_bb = I_ERI_F2xy_D2y_Pz_S_bb+ABY*I_ERI_D2x_D2y_Pz_S_bb;
  Double I_ERI_Dxy_F3y_Pz_S_bb = I_ERI_Fx2y_D2y_Pz_S_bb+ABY*I_ERI_Dxy_D2y_Pz_S_bb;
  Double I_ERI_Dxz_F3y_Pz_S_bb = I_ERI_Fxyz_D2y_Pz_S_bb+ABY*I_ERI_Dxz_D2y_Pz_S_bb;
  Double I_ERI_D2y_F3y_Pz_S_bb = I_ERI_F3y_D2y_Pz_S_bb+ABY*I_ERI_D2y_D2y_Pz_S_bb;
  Double I_ERI_Dyz_F3y_Pz_S_bb = I_ERI_F2yz_D2y_Pz_S_bb+ABY*I_ERI_Dyz_D2y_Pz_S_bb;
  Double I_ERI_D2z_F3y_Pz_S_bb = I_ERI_Fy2z_D2y_Pz_S_bb+ABY*I_ERI_D2z_D2y_Pz_S_bb;
  Double I_ERI_D2x_F2yz_Pz_S_bb = I_ERI_F2xz_D2y_Pz_S_bb+ABZ*I_ERI_D2x_D2y_Pz_S_bb;
  Double I_ERI_Dxy_F2yz_Pz_S_bb = I_ERI_Fxyz_D2y_Pz_S_bb+ABZ*I_ERI_Dxy_D2y_Pz_S_bb;
  Double I_ERI_Dxz_F2yz_Pz_S_bb = I_ERI_Fx2z_D2y_Pz_S_bb+ABZ*I_ERI_Dxz_D2y_Pz_S_bb;
  Double I_ERI_D2y_F2yz_Pz_S_bb = I_ERI_F2yz_D2y_Pz_S_bb+ABZ*I_ERI_D2y_D2y_Pz_S_bb;
  Double I_ERI_Dyz_F2yz_Pz_S_bb = I_ERI_Fy2z_D2y_Pz_S_bb+ABZ*I_ERI_Dyz_D2y_Pz_S_bb;
  Double I_ERI_D2z_F2yz_Pz_S_bb = I_ERI_F3z_D2y_Pz_S_bb+ABZ*I_ERI_D2z_D2y_Pz_S_bb;
  Double I_ERI_D2x_Fy2z_Pz_S_bb = I_ERI_F2xy_D2z_Pz_S_bb+ABY*I_ERI_D2x_D2z_Pz_S_bb;
  Double I_ERI_Dxy_Fy2z_Pz_S_bb = I_ERI_Fx2y_D2z_Pz_S_bb+ABY*I_ERI_Dxy_D2z_Pz_S_bb;
  Double I_ERI_Dxz_Fy2z_Pz_S_bb = I_ERI_Fxyz_D2z_Pz_S_bb+ABY*I_ERI_Dxz_D2z_Pz_S_bb;
  Double I_ERI_D2y_Fy2z_Pz_S_bb = I_ERI_F3y_D2z_Pz_S_bb+ABY*I_ERI_D2y_D2z_Pz_S_bb;
  Double I_ERI_Dyz_Fy2z_Pz_S_bb = I_ERI_F2yz_D2z_Pz_S_bb+ABY*I_ERI_Dyz_D2z_Pz_S_bb;
  Double I_ERI_D2z_Fy2z_Pz_S_bb = I_ERI_Fy2z_D2z_Pz_S_bb+ABY*I_ERI_D2z_D2z_Pz_S_bb;
  Double I_ERI_D2x_F3z_Pz_S_bb = I_ERI_F2xz_D2z_Pz_S_bb+ABZ*I_ERI_D2x_D2z_Pz_S_bb;
  Double I_ERI_Dxy_F3z_Pz_S_bb = I_ERI_Fxyz_D2z_Pz_S_bb+ABZ*I_ERI_Dxy_D2z_Pz_S_bb;
  Double I_ERI_Dxz_F3z_Pz_S_bb = I_ERI_Fx2z_D2z_Pz_S_bb+ABZ*I_ERI_Dxz_D2z_Pz_S_bb;
  Double I_ERI_D2y_F3z_Pz_S_bb = I_ERI_F2yz_D2z_Pz_S_bb+ABZ*I_ERI_D2y_D2z_Pz_S_bb;
  Double I_ERI_Dyz_F3z_Pz_S_bb = I_ERI_Fy2z_D2z_Pz_S_bb+ABZ*I_ERI_Dyz_D2z_Pz_S_bb;
  Double I_ERI_D2z_F3z_Pz_S_bb = I_ERI_F3z_D2z_Pz_S_bb+ABZ*I_ERI_D2z_D2z_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_D_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_S_D_S_bc
   ************************************************************/
  Double I_ERI_D2x_Px_D2x_S_bc = I_ERI_F3x_S_D2x_S_bc+ABX*I_ERI_D2x_S_D2x_S_bc;
  Double I_ERI_Dxy_Px_D2x_S_bc = I_ERI_F2xy_S_D2x_S_bc+ABX*I_ERI_Dxy_S_D2x_S_bc;
  Double I_ERI_Dxz_Px_D2x_S_bc = I_ERI_F2xz_S_D2x_S_bc+ABX*I_ERI_Dxz_S_D2x_S_bc;
  Double I_ERI_D2y_Px_D2x_S_bc = I_ERI_Fx2y_S_D2x_S_bc+ABX*I_ERI_D2y_S_D2x_S_bc;
  Double I_ERI_Dyz_Px_D2x_S_bc = I_ERI_Fxyz_S_D2x_S_bc+ABX*I_ERI_Dyz_S_D2x_S_bc;
  Double I_ERI_D2z_Px_D2x_S_bc = I_ERI_Fx2z_S_D2x_S_bc+ABX*I_ERI_D2z_S_D2x_S_bc;
  Double I_ERI_D2x_Py_D2x_S_bc = I_ERI_F2xy_S_D2x_S_bc+ABY*I_ERI_D2x_S_D2x_S_bc;
  Double I_ERI_Dxy_Py_D2x_S_bc = I_ERI_Fx2y_S_D2x_S_bc+ABY*I_ERI_Dxy_S_D2x_S_bc;
  Double I_ERI_Dxz_Py_D2x_S_bc = I_ERI_Fxyz_S_D2x_S_bc+ABY*I_ERI_Dxz_S_D2x_S_bc;
  Double I_ERI_D2y_Py_D2x_S_bc = I_ERI_F3y_S_D2x_S_bc+ABY*I_ERI_D2y_S_D2x_S_bc;
  Double I_ERI_Dyz_Py_D2x_S_bc = I_ERI_F2yz_S_D2x_S_bc+ABY*I_ERI_Dyz_S_D2x_S_bc;
  Double I_ERI_D2z_Py_D2x_S_bc = I_ERI_Fy2z_S_D2x_S_bc+ABY*I_ERI_D2z_S_D2x_S_bc;
  Double I_ERI_D2x_Pz_D2x_S_bc = I_ERI_F2xz_S_D2x_S_bc+ABZ*I_ERI_D2x_S_D2x_S_bc;
  Double I_ERI_Dxy_Pz_D2x_S_bc = I_ERI_Fxyz_S_D2x_S_bc+ABZ*I_ERI_Dxy_S_D2x_S_bc;
  Double I_ERI_Dxz_Pz_D2x_S_bc = I_ERI_Fx2z_S_D2x_S_bc+ABZ*I_ERI_Dxz_S_D2x_S_bc;
  Double I_ERI_D2y_Pz_D2x_S_bc = I_ERI_F2yz_S_D2x_S_bc+ABZ*I_ERI_D2y_S_D2x_S_bc;
  Double I_ERI_Dyz_Pz_D2x_S_bc = I_ERI_Fy2z_S_D2x_S_bc+ABZ*I_ERI_Dyz_S_D2x_S_bc;
  Double I_ERI_D2z_Pz_D2x_S_bc = I_ERI_F3z_S_D2x_S_bc+ABZ*I_ERI_D2z_S_D2x_S_bc;
  Double I_ERI_D2x_Px_Dxy_S_bc = I_ERI_F3x_S_Dxy_S_bc+ABX*I_ERI_D2x_S_Dxy_S_bc;
  Double I_ERI_Dxy_Px_Dxy_S_bc = I_ERI_F2xy_S_Dxy_S_bc+ABX*I_ERI_Dxy_S_Dxy_S_bc;
  Double I_ERI_Dxz_Px_Dxy_S_bc = I_ERI_F2xz_S_Dxy_S_bc+ABX*I_ERI_Dxz_S_Dxy_S_bc;
  Double I_ERI_D2y_Px_Dxy_S_bc = I_ERI_Fx2y_S_Dxy_S_bc+ABX*I_ERI_D2y_S_Dxy_S_bc;
  Double I_ERI_Dyz_Px_Dxy_S_bc = I_ERI_Fxyz_S_Dxy_S_bc+ABX*I_ERI_Dyz_S_Dxy_S_bc;
  Double I_ERI_D2z_Px_Dxy_S_bc = I_ERI_Fx2z_S_Dxy_S_bc+ABX*I_ERI_D2z_S_Dxy_S_bc;
  Double I_ERI_D2x_Py_Dxy_S_bc = I_ERI_F2xy_S_Dxy_S_bc+ABY*I_ERI_D2x_S_Dxy_S_bc;
  Double I_ERI_Dxy_Py_Dxy_S_bc = I_ERI_Fx2y_S_Dxy_S_bc+ABY*I_ERI_Dxy_S_Dxy_S_bc;
  Double I_ERI_Dxz_Py_Dxy_S_bc = I_ERI_Fxyz_S_Dxy_S_bc+ABY*I_ERI_Dxz_S_Dxy_S_bc;
  Double I_ERI_D2y_Py_Dxy_S_bc = I_ERI_F3y_S_Dxy_S_bc+ABY*I_ERI_D2y_S_Dxy_S_bc;
  Double I_ERI_Dyz_Py_Dxy_S_bc = I_ERI_F2yz_S_Dxy_S_bc+ABY*I_ERI_Dyz_S_Dxy_S_bc;
  Double I_ERI_D2z_Py_Dxy_S_bc = I_ERI_Fy2z_S_Dxy_S_bc+ABY*I_ERI_D2z_S_Dxy_S_bc;
  Double I_ERI_D2x_Pz_Dxy_S_bc = I_ERI_F2xz_S_Dxy_S_bc+ABZ*I_ERI_D2x_S_Dxy_S_bc;
  Double I_ERI_Dxy_Pz_Dxy_S_bc = I_ERI_Fxyz_S_Dxy_S_bc+ABZ*I_ERI_Dxy_S_Dxy_S_bc;
  Double I_ERI_Dxz_Pz_Dxy_S_bc = I_ERI_Fx2z_S_Dxy_S_bc+ABZ*I_ERI_Dxz_S_Dxy_S_bc;
  Double I_ERI_D2y_Pz_Dxy_S_bc = I_ERI_F2yz_S_Dxy_S_bc+ABZ*I_ERI_D2y_S_Dxy_S_bc;
  Double I_ERI_Dyz_Pz_Dxy_S_bc = I_ERI_Fy2z_S_Dxy_S_bc+ABZ*I_ERI_Dyz_S_Dxy_S_bc;
  Double I_ERI_D2z_Pz_Dxy_S_bc = I_ERI_F3z_S_Dxy_S_bc+ABZ*I_ERI_D2z_S_Dxy_S_bc;
  Double I_ERI_D2x_Px_Dxz_S_bc = I_ERI_F3x_S_Dxz_S_bc+ABX*I_ERI_D2x_S_Dxz_S_bc;
  Double I_ERI_Dxy_Px_Dxz_S_bc = I_ERI_F2xy_S_Dxz_S_bc+ABX*I_ERI_Dxy_S_Dxz_S_bc;
  Double I_ERI_Dxz_Px_Dxz_S_bc = I_ERI_F2xz_S_Dxz_S_bc+ABX*I_ERI_Dxz_S_Dxz_S_bc;
  Double I_ERI_D2y_Px_Dxz_S_bc = I_ERI_Fx2y_S_Dxz_S_bc+ABX*I_ERI_D2y_S_Dxz_S_bc;
  Double I_ERI_Dyz_Px_Dxz_S_bc = I_ERI_Fxyz_S_Dxz_S_bc+ABX*I_ERI_Dyz_S_Dxz_S_bc;
  Double I_ERI_D2z_Px_Dxz_S_bc = I_ERI_Fx2z_S_Dxz_S_bc+ABX*I_ERI_D2z_S_Dxz_S_bc;
  Double I_ERI_D2x_Py_Dxz_S_bc = I_ERI_F2xy_S_Dxz_S_bc+ABY*I_ERI_D2x_S_Dxz_S_bc;
  Double I_ERI_Dxy_Py_Dxz_S_bc = I_ERI_Fx2y_S_Dxz_S_bc+ABY*I_ERI_Dxy_S_Dxz_S_bc;
  Double I_ERI_Dxz_Py_Dxz_S_bc = I_ERI_Fxyz_S_Dxz_S_bc+ABY*I_ERI_Dxz_S_Dxz_S_bc;
  Double I_ERI_D2y_Py_Dxz_S_bc = I_ERI_F3y_S_Dxz_S_bc+ABY*I_ERI_D2y_S_Dxz_S_bc;
  Double I_ERI_Dyz_Py_Dxz_S_bc = I_ERI_F2yz_S_Dxz_S_bc+ABY*I_ERI_Dyz_S_Dxz_S_bc;
  Double I_ERI_D2z_Py_Dxz_S_bc = I_ERI_Fy2z_S_Dxz_S_bc+ABY*I_ERI_D2z_S_Dxz_S_bc;
  Double I_ERI_D2x_Pz_Dxz_S_bc = I_ERI_F2xz_S_Dxz_S_bc+ABZ*I_ERI_D2x_S_Dxz_S_bc;
  Double I_ERI_Dxy_Pz_Dxz_S_bc = I_ERI_Fxyz_S_Dxz_S_bc+ABZ*I_ERI_Dxy_S_Dxz_S_bc;
  Double I_ERI_Dxz_Pz_Dxz_S_bc = I_ERI_Fx2z_S_Dxz_S_bc+ABZ*I_ERI_Dxz_S_Dxz_S_bc;
  Double I_ERI_D2y_Pz_Dxz_S_bc = I_ERI_F2yz_S_Dxz_S_bc+ABZ*I_ERI_D2y_S_Dxz_S_bc;
  Double I_ERI_Dyz_Pz_Dxz_S_bc = I_ERI_Fy2z_S_Dxz_S_bc+ABZ*I_ERI_Dyz_S_Dxz_S_bc;
  Double I_ERI_D2z_Pz_Dxz_S_bc = I_ERI_F3z_S_Dxz_S_bc+ABZ*I_ERI_D2z_S_Dxz_S_bc;
  Double I_ERI_D2x_Px_D2y_S_bc = I_ERI_F3x_S_D2y_S_bc+ABX*I_ERI_D2x_S_D2y_S_bc;
  Double I_ERI_Dxy_Px_D2y_S_bc = I_ERI_F2xy_S_D2y_S_bc+ABX*I_ERI_Dxy_S_D2y_S_bc;
  Double I_ERI_Dxz_Px_D2y_S_bc = I_ERI_F2xz_S_D2y_S_bc+ABX*I_ERI_Dxz_S_D2y_S_bc;
  Double I_ERI_D2y_Px_D2y_S_bc = I_ERI_Fx2y_S_D2y_S_bc+ABX*I_ERI_D2y_S_D2y_S_bc;
  Double I_ERI_Dyz_Px_D2y_S_bc = I_ERI_Fxyz_S_D2y_S_bc+ABX*I_ERI_Dyz_S_D2y_S_bc;
  Double I_ERI_D2z_Px_D2y_S_bc = I_ERI_Fx2z_S_D2y_S_bc+ABX*I_ERI_D2z_S_D2y_S_bc;
  Double I_ERI_D2x_Py_D2y_S_bc = I_ERI_F2xy_S_D2y_S_bc+ABY*I_ERI_D2x_S_D2y_S_bc;
  Double I_ERI_Dxy_Py_D2y_S_bc = I_ERI_Fx2y_S_D2y_S_bc+ABY*I_ERI_Dxy_S_D2y_S_bc;
  Double I_ERI_Dxz_Py_D2y_S_bc = I_ERI_Fxyz_S_D2y_S_bc+ABY*I_ERI_Dxz_S_D2y_S_bc;
  Double I_ERI_D2y_Py_D2y_S_bc = I_ERI_F3y_S_D2y_S_bc+ABY*I_ERI_D2y_S_D2y_S_bc;
  Double I_ERI_Dyz_Py_D2y_S_bc = I_ERI_F2yz_S_D2y_S_bc+ABY*I_ERI_Dyz_S_D2y_S_bc;
  Double I_ERI_D2z_Py_D2y_S_bc = I_ERI_Fy2z_S_D2y_S_bc+ABY*I_ERI_D2z_S_D2y_S_bc;
  Double I_ERI_D2x_Pz_D2y_S_bc = I_ERI_F2xz_S_D2y_S_bc+ABZ*I_ERI_D2x_S_D2y_S_bc;
  Double I_ERI_Dxy_Pz_D2y_S_bc = I_ERI_Fxyz_S_D2y_S_bc+ABZ*I_ERI_Dxy_S_D2y_S_bc;
  Double I_ERI_Dxz_Pz_D2y_S_bc = I_ERI_Fx2z_S_D2y_S_bc+ABZ*I_ERI_Dxz_S_D2y_S_bc;
  Double I_ERI_D2y_Pz_D2y_S_bc = I_ERI_F2yz_S_D2y_S_bc+ABZ*I_ERI_D2y_S_D2y_S_bc;
  Double I_ERI_Dyz_Pz_D2y_S_bc = I_ERI_Fy2z_S_D2y_S_bc+ABZ*I_ERI_Dyz_S_D2y_S_bc;
  Double I_ERI_D2z_Pz_D2y_S_bc = I_ERI_F3z_S_D2y_S_bc+ABZ*I_ERI_D2z_S_D2y_S_bc;
  Double I_ERI_D2x_Px_Dyz_S_bc = I_ERI_F3x_S_Dyz_S_bc+ABX*I_ERI_D2x_S_Dyz_S_bc;
  Double I_ERI_Dxy_Px_Dyz_S_bc = I_ERI_F2xy_S_Dyz_S_bc+ABX*I_ERI_Dxy_S_Dyz_S_bc;
  Double I_ERI_Dxz_Px_Dyz_S_bc = I_ERI_F2xz_S_Dyz_S_bc+ABX*I_ERI_Dxz_S_Dyz_S_bc;
  Double I_ERI_D2y_Px_Dyz_S_bc = I_ERI_Fx2y_S_Dyz_S_bc+ABX*I_ERI_D2y_S_Dyz_S_bc;
  Double I_ERI_Dyz_Px_Dyz_S_bc = I_ERI_Fxyz_S_Dyz_S_bc+ABX*I_ERI_Dyz_S_Dyz_S_bc;
  Double I_ERI_D2z_Px_Dyz_S_bc = I_ERI_Fx2z_S_Dyz_S_bc+ABX*I_ERI_D2z_S_Dyz_S_bc;
  Double I_ERI_D2x_Py_Dyz_S_bc = I_ERI_F2xy_S_Dyz_S_bc+ABY*I_ERI_D2x_S_Dyz_S_bc;
  Double I_ERI_Dxy_Py_Dyz_S_bc = I_ERI_Fx2y_S_Dyz_S_bc+ABY*I_ERI_Dxy_S_Dyz_S_bc;
  Double I_ERI_Dxz_Py_Dyz_S_bc = I_ERI_Fxyz_S_Dyz_S_bc+ABY*I_ERI_Dxz_S_Dyz_S_bc;
  Double I_ERI_D2y_Py_Dyz_S_bc = I_ERI_F3y_S_Dyz_S_bc+ABY*I_ERI_D2y_S_Dyz_S_bc;
  Double I_ERI_Dyz_Py_Dyz_S_bc = I_ERI_F2yz_S_Dyz_S_bc+ABY*I_ERI_Dyz_S_Dyz_S_bc;
  Double I_ERI_D2z_Py_Dyz_S_bc = I_ERI_Fy2z_S_Dyz_S_bc+ABY*I_ERI_D2z_S_Dyz_S_bc;
  Double I_ERI_D2x_Pz_Dyz_S_bc = I_ERI_F2xz_S_Dyz_S_bc+ABZ*I_ERI_D2x_S_Dyz_S_bc;
  Double I_ERI_Dxy_Pz_Dyz_S_bc = I_ERI_Fxyz_S_Dyz_S_bc+ABZ*I_ERI_Dxy_S_Dyz_S_bc;
  Double I_ERI_Dxz_Pz_Dyz_S_bc = I_ERI_Fx2z_S_Dyz_S_bc+ABZ*I_ERI_Dxz_S_Dyz_S_bc;
  Double I_ERI_D2y_Pz_Dyz_S_bc = I_ERI_F2yz_S_Dyz_S_bc+ABZ*I_ERI_D2y_S_Dyz_S_bc;
  Double I_ERI_Dyz_Pz_Dyz_S_bc = I_ERI_Fy2z_S_Dyz_S_bc+ABZ*I_ERI_Dyz_S_Dyz_S_bc;
  Double I_ERI_D2z_Pz_Dyz_S_bc = I_ERI_F3z_S_Dyz_S_bc+ABZ*I_ERI_D2z_S_Dyz_S_bc;
  Double I_ERI_D2x_Px_D2z_S_bc = I_ERI_F3x_S_D2z_S_bc+ABX*I_ERI_D2x_S_D2z_S_bc;
  Double I_ERI_Dxy_Px_D2z_S_bc = I_ERI_F2xy_S_D2z_S_bc+ABX*I_ERI_Dxy_S_D2z_S_bc;
  Double I_ERI_Dxz_Px_D2z_S_bc = I_ERI_F2xz_S_D2z_S_bc+ABX*I_ERI_Dxz_S_D2z_S_bc;
  Double I_ERI_D2y_Px_D2z_S_bc = I_ERI_Fx2y_S_D2z_S_bc+ABX*I_ERI_D2y_S_D2z_S_bc;
  Double I_ERI_Dyz_Px_D2z_S_bc = I_ERI_Fxyz_S_D2z_S_bc+ABX*I_ERI_Dyz_S_D2z_S_bc;
  Double I_ERI_D2z_Px_D2z_S_bc = I_ERI_Fx2z_S_D2z_S_bc+ABX*I_ERI_D2z_S_D2z_S_bc;
  Double I_ERI_D2x_Py_D2z_S_bc = I_ERI_F2xy_S_D2z_S_bc+ABY*I_ERI_D2x_S_D2z_S_bc;
  Double I_ERI_Dxy_Py_D2z_S_bc = I_ERI_Fx2y_S_D2z_S_bc+ABY*I_ERI_Dxy_S_D2z_S_bc;
  Double I_ERI_Dxz_Py_D2z_S_bc = I_ERI_Fxyz_S_D2z_S_bc+ABY*I_ERI_Dxz_S_D2z_S_bc;
  Double I_ERI_D2y_Py_D2z_S_bc = I_ERI_F3y_S_D2z_S_bc+ABY*I_ERI_D2y_S_D2z_S_bc;
  Double I_ERI_Dyz_Py_D2z_S_bc = I_ERI_F2yz_S_D2z_S_bc+ABY*I_ERI_Dyz_S_D2z_S_bc;
  Double I_ERI_D2z_Py_D2z_S_bc = I_ERI_Fy2z_S_D2z_S_bc+ABY*I_ERI_D2z_S_D2z_S_bc;
  Double I_ERI_D2x_Pz_D2z_S_bc = I_ERI_F2xz_S_D2z_S_bc+ABZ*I_ERI_D2x_S_D2z_S_bc;
  Double I_ERI_Dxy_Pz_D2z_S_bc = I_ERI_Fxyz_S_D2z_S_bc+ABZ*I_ERI_Dxy_S_D2z_S_bc;
  Double I_ERI_Dxz_Pz_D2z_S_bc = I_ERI_Fx2z_S_D2z_S_bc+ABZ*I_ERI_Dxz_S_D2z_S_bc;
  Double I_ERI_D2y_Pz_D2z_S_bc = I_ERI_F2yz_S_D2z_S_bc+ABZ*I_ERI_D2y_S_D2z_S_bc;
  Double I_ERI_Dyz_Pz_D2z_S_bc = I_ERI_Fy2z_S_D2z_S_bc+ABZ*I_ERI_Dyz_S_D2z_S_bc;
  Double I_ERI_D2z_Pz_D2z_S_bc = I_ERI_F3z_S_D2z_S_bc+ABZ*I_ERI_D2z_S_D2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_bc
   * RHS shell quartet name: SQ_ERI_F_S_D_S_bc
   ************************************************************/
  Double I_ERI_F3x_Px_D2x_S_bc = I_ERI_G4x_S_D2x_S_bc+ABX*I_ERI_F3x_S_D2x_S_bc;
  Double I_ERI_F2xy_Px_D2x_S_bc = I_ERI_G3xy_S_D2x_S_bc+ABX*I_ERI_F2xy_S_D2x_S_bc;
  Double I_ERI_F2xz_Px_D2x_S_bc = I_ERI_G3xz_S_D2x_S_bc+ABX*I_ERI_F2xz_S_D2x_S_bc;
  Double I_ERI_Fx2y_Px_D2x_S_bc = I_ERI_G2x2y_S_D2x_S_bc+ABX*I_ERI_Fx2y_S_D2x_S_bc;
  Double I_ERI_Fxyz_Px_D2x_S_bc = I_ERI_G2xyz_S_D2x_S_bc+ABX*I_ERI_Fxyz_S_D2x_S_bc;
  Double I_ERI_Fx2z_Px_D2x_S_bc = I_ERI_G2x2z_S_D2x_S_bc+ABX*I_ERI_Fx2z_S_D2x_S_bc;
  Double I_ERI_F3y_Px_D2x_S_bc = I_ERI_Gx3y_S_D2x_S_bc+ABX*I_ERI_F3y_S_D2x_S_bc;
  Double I_ERI_F2yz_Px_D2x_S_bc = I_ERI_Gx2yz_S_D2x_S_bc+ABX*I_ERI_F2yz_S_D2x_S_bc;
  Double I_ERI_Fy2z_Px_D2x_S_bc = I_ERI_Gxy2z_S_D2x_S_bc+ABX*I_ERI_Fy2z_S_D2x_S_bc;
  Double I_ERI_F3z_Px_D2x_S_bc = I_ERI_Gx3z_S_D2x_S_bc+ABX*I_ERI_F3z_S_D2x_S_bc;
  Double I_ERI_F2xy_Py_D2x_S_bc = I_ERI_G2x2y_S_D2x_S_bc+ABY*I_ERI_F2xy_S_D2x_S_bc;
  Double I_ERI_F2xz_Py_D2x_S_bc = I_ERI_G2xyz_S_D2x_S_bc+ABY*I_ERI_F2xz_S_D2x_S_bc;
  Double I_ERI_Fx2y_Py_D2x_S_bc = I_ERI_Gx3y_S_D2x_S_bc+ABY*I_ERI_Fx2y_S_D2x_S_bc;
  Double I_ERI_Fxyz_Py_D2x_S_bc = I_ERI_Gx2yz_S_D2x_S_bc+ABY*I_ERI_Fxyz_S_D2x_S_bc;
  Double I_ERI_Fx2z_Py_D2x_S_bc = I_ERI_Gxy2z_S_D2x_S_bc+ABY*I_ERI_Fx2z_S_D2x_S_bc;
  Double I_ERI_F3y_Py_D2x_S_bc = I_ERI_G4y_S_D2x_S_bc+ABY*I_ERI_F3y_S_D2x_S_bc;
  Double I_ERI_F2yz_Py_D2x_S_bc = I_ERI_G3yz_S_D2x_S_bc+ABY*I_ERI_F2yz_S_D2x_S_bc;
  Double I_ERI_Fy2z_Py_D2x_S_bc = I_ERI_G2y2z_S_D2x_S_bc+ABY*I_ERI_Fy2z_S_D2x_S_bc;
  Double I_ERI_F3z_Py_D2x_S_bc = I_ERI_Gy3z_S_D2x_S_bc+ABY*I_ERI_F3z_S_D2x_S_bc;
  Double I_ERI_F2xz_Pz_D2x_S_bc = I_ERI_G2x2z_S_D2x_S_bc+ABZ*I_ERI_F2xz_S_D2x_S_bc;
  Double I_ERI_Fxyz_Pz_D2x_S_bc = I_ERI_Gxy2z_S_D2x_S_bc+ABZ*I_ERI_Fxyz_S_D2x_S_bc;
  Double I_ERI_Fx2z_Pz_D2x_S_bc = I_ERI_Gx3z_S_D2x_S_bc+ABZ*I_ERI_Fx2z_S_D2x_S_bc;
  Double I_ERI_F2yz_Pz_D2x_S_bc = I_ERI_G2y2z_S_D2x_S_bc+ABZ*I_ERI_F2yz_S_D2x_S_bc;
  Double I_ERI_Fy2z_Pz_D2x_S_bc = I_ERI_Gy3z_S_D2x_S_bc+ABZ*I_ERI_Fy2z_S_D2x_S_bc;
  Double I_ERI_F3z_Pz_D2x_S_bc = I_ERI_G4z_S_D2x_S_bc+ABZ*I_ERI_F3z_S_D2x_S_bc;
  Double I_ERI_F3x_Px_Dxy_S_bc = I_ERI_G4x_S_Dxy_S_bc+ABX*I_ERI_F3x_S_Dxy_S_bc;
  Double I_ERI_F2xy_Px_Dxy_S_bc = I_ERI_G3xy_S_Dxy_S_bc+ABX*I_ERI_F2xy_S_Dxy_S_bc;
  Double I_ERI_F2xz_Px_Dxy_S_bc = I_ERI_G3xz_S_Dxy_S_bc+ABX*I_ERI_F2xz_S_Dxy_S_bc;
  Double I_ERI_Fx2y_Px_Dxy_S_bc = I_ERI_G2x2y_S_Dxy_S_bc+ABX*I_ERI_Fx2y_S_Dxy_S_bc;
  Double I_ERI_Fxyz_Px_Dxy_S_bc = I_ERI_G2xyz_S_Dxy_S_bc+ABX*I_ERI_Fxyz_S_Dxy_S_bc;
  Double I_ERI_Fx2z_Px_Dxy_S_bc = I_ERI_G2x2z_S_Dxy_S_bc+ABX*I_ERI_Fx2z_S_Dxy_S_bc;
  Double I_ERI_F3y_Px_Dxy_S_bc = I_ERI_Gx3y_S_Dxy_S_bc+ABX*I_ERI_F3y_S_Dxy_S_bc;
  Double I_ERI_F2yz_Px_Dxy_S_bc = I_ERI_Gx2yz_S_Dxy_S_bc+ABX*I_ERI_F2yz_S_Dxy_S_bc;
  Double I_ERI_Fy2z_Px_Dxy_S_bc = I_ERI_Gxy2z_S_Dxy_S_bc+ABX*I_ERI_Fy2z_S_Dxy_S_bc;
  Double I_ERI_F3z_Px_Dxy_S_bc = I_ERI_Gx3z_S_Dxy_S_bc+ABX*I_ERI_F3z_S_Dxy_S_bc;
  Double I_ERI_F2xy_Py_Dxy_S_bc = I_ERI_G2x2y_S_Dxy_S_bc+ABY*I_ERI_F2xy_S_Dxy_S_bc;
  Double I_ERI_F2xz_Py_Dxy_S_bc = I_ERI_G2xyz_S_Dxy_S_bc+ABY*I_ERI_F2xz_S_Dxy_S_bc;
  Double I_ERI_Fx2y_Py_Dxy_S_bc = I_ERI_Gx3y_S_Dxy_S_bc+ABY*I_ERI_Fx2y_S_Dxy_S_bc;
  Double I_ERI_Fxyz_Py_Dxy_S_bc = I_ERI_Gx2yz_S_Dxy_S_bc+ABY*I_ERI_Fxyz_S_Dxy_S_bc;
  Double I_ERI_Fx2z_Py_Dxy_S_bc = I_ERI_Gxy2z_S_Dxy_S_bc+ABY*I_ERI_Fx2z_S_Dxy_S_bc;
  Double I_ERI_F3y_Py_Dxy_S_bc = I_ERI_G4y_S_Dxy_S_bc+ABY*I_ERI_F3y_S_Dxy_S_bc;
  Double I_ERI_F2yz_Py_Dxy_S_bc = I_ERI_G3yz_S_Dxy_S_bc+ABY*I_ERI_F2yz_S_Dxy_S_bc;
  Double I_ERI_Fy2z_Py_Dxy_S_bc = I_ERI_G2y2z_S_Dxy_S_bc+ABY*I_ERI_Fy2z_S_Dxy_S_bc;
  Double I_ERI_F3z_Py_Dxy_S_bc = I_ERI_Gy3z_S_Dxy_S_bc+ABY*I_ERI_F3z_S_Dxy_S_bc;
  Double I_ERI_F2xz_Pz_Dxy_S_bc = I_ERI_G2x2z_S_Dxy_S_bc+ABZ*I_ERI_F2xz_S_Dxy_S_bc;
  Double I_ERI_Fxyz_Pz_Dxy_S_bc = I_ERI_Gxy2z_S_Dxy_S_bc+ABZ*I_ERI_Fxyz_S_Dxy_S_bc;
  Double I_ERI_Fx2z_Pz_Dxy_S_bc = I_ERI_Gx3z_S_Dxy_S_bc+ABZ*I_ERI_Fx2z_S_Dxy_S_bc;
  Double I_ERI_F2yz_Pz_Dxy_S_bc = I_ERI_G2y2z_S_Dxy_S_bc+ABZ*I_ERI_F2yz_S_Dxy_S_bc;
  Double I_ERI_Fy2z_Pz_Dxy_S_bc = I_ERI_Gy3z_S_Dxy_S_bc+ABZ*I_ERI_Fy2z_S_Dxy_S_bc;
  Double I_ERI_F3z_Pz_Dxy_S_bc = I_ERI_G4z_S_Dxy_S_bc+ABZ*I_ERI_F3z_S_Dxy_S_bc;
  Double I_ERI_F3x_Px_Dxz_S_bc = I_ERI_G4x_S_Dxz_S_bc+ABX*I_ERI_F3x_S_Dxz_S_bc;
  Double I_ERI_F2xy_Px_Dxz_S_bc = I_ERI_G3xy_S_Dxz_S_bc+ABX*I_ERI_F2xy_S_Dxz_S_bc;
  Double I_ERI_F2xz_Px_Dxz_S_bc = I_ERI_G3xz_S_Dxz_S_bc+ABX*I_ERI_F2xz_S_Dxz_S_bc;
  Double I_ERI_Fx2y_Px_Dxz_S_bc = I_ERI_G2x2y_S_Dxz_S_bc+ABX*I_ERI_Fx2y_S_Dxz_S_bc;
  Double I_ERI_Fxyz_Px_Dxz_S_bc = I_ERI_G2xyz_S_Dxz_S_bc+ABX*I_ERI_Fxyz_S_Dxz_S_bc;
  Double I_ERI_Fx2z_Px_Dxz_S_bc = I_ERI_G2x2z_S_Dxz_S_bc+ABX*I_ERI_Fx2z_S_Dxz_S_bc;
  Double I_ERI_F3y_Px_Dxz_S_bc = I_ERI_Gx3y_S_Dxz_S_bc+ABX*I_ERI_F3y_S_Dxz_S_bc;
  Double I_ERI_F2yz_Px_Dxz_S_bc = I_ERI_Gx2yz_S_Dxz_S_bc+ABX*I_ERI_F2yz_S_Dxz_S_bc;
  Double I_ERI_Fy2z_Px_Dxz_S_bc = I_ERI_Gxy2z_S_Dxz_S_bc+ABX*I_ERI_Fy2z_S_Dxz_S_bc;
  Double I_ERI_F3z_Px_Dxz_S_bc = I_ERI_Gx3z_S_Dxz_S_bc+ABX*I_ERI_F3z_S_Dxz_S_bc;
  Double I_ERI_F2xy_Py_Dxz_S_bc = I_ERI_G2x2y_S_Dxz_S_bc+ABY*I_ERI_F2xy_S_Dxz_S_bc;
  Double I_ERI_F2xz_Py_Dxz_S_bc = I_ERI_G2xyz_S_Dxz_S_bc+ABY*I_ERI_F2xz_S_Dxz_S_bc;
  Double I_ERI_Fx2y_Py_Dxz_S_bc = I_ERI_Gx3y_S_Dxz_S_bc+ABY*I_ERI_Fx2y_S_Dxz_S_bc;
  Double I_ERI_Fxyz_Py_Dxz_S_bc = I_ERI_Gx2yz_S_Dxz_S_bc+ABY*I_ERI_Fxyz_S_Dxz_S_bc;
  Double I_ERI_Fx2z_Py_Dxz_S_bc = I_ERI_Gxy2z_S_Dxz_S_bc+ABY*I_ERI_Fx2z_S_Dxz_S_bc;
  Double I_ERI_F3y_Py_Dxz_S_bc = I_ERI_G4y_S_Dxz_S_bc+ABY*I_ERI_F3y_S_Dxz_S_bc;
  Double I_ERI_F2yz_Py_Dxz_S_bc = I_ERI_G3yz_S_Dxz_S_bc+ABY*I_ERI_F2yz_S_Dxz_S_bc;
  Double I_ERI_Fy2z_Py_Dxz_S_bc = I_ERI_G2y2z_S_Dxz_S_bc+ABY*I_ERI_Fy2z_S_Dxz_S_bc;
  Double I_ERI_F3z_Py_Dxz_S_bc = I_ERI_Gy3z_S_Dxz_S_bc+ABY*I_ERI_F3z_S_Dxz_S_bc;
  Double I_ERI_F2xz_Pz_Dxz_S_bc = I_ERI_G2x2z_S_Dxz_S_bc+ABZ*I_ERI_F2xz_S_Dxz_S_bc;
  Double I_ERI_Fxyz_Pz_Dxz_S_bc = I_ERI_Gxy2z_S_Dxz_S_bc+ABZ*I_ERI_Fxyz_S_Dxz_S_bc;
  Double I_ERI_Fx2z_Pz_Dxz_S_bc = I_ERI_Gx3z_S_Dxz_S_bc+ABZ*I_ERI_Fx2z_S_Dxz_S_bc;
  Double I_ERI_F2yz_Pz_Dxz_S_bc = I_ERI_G2y2z_S_Dxz_S_bc+ABZ*I_ERI_F2yz_S_Dxz_S_bc;
  Double I_ERI_Fy2z_Pz_Dxz_S_bc = I_ERI_Gy3z_S_Dxz_S_bc+ABZ*I_ERI_Fy2z_S_Dxz_S_bc;
  Double I_ERI_F3z_Pz_Dxz_S_bc = I_ERI_G4z_S_Dxz_S_bc+ABZ*I_ERI_F3z_S_Dxz_S_bc;
  Double I_ERI_F3x_Px_D2y_S_bc = I_ERI_G4x_S_D2y_S_bc+ABX*I_ERI_F3x_S_D2y_S_bc;
  Double I_ERI_F2xy_Px_D2y_S_bc = I_ERI_G3xy_S_D2y_S_bc+ABX*I_ERI_F2xy_S_D2y_S_bc;
  Double I_ERI_F2xz_Px_D2y_S_bc = I_ERI_G3xz_S_D2y_S_bc+ABX*I_ERI_F2xz_S_D2y_S_bc;
  Double I_ERI_Fx2y_Px_D2y_S_bc = I_ERI_G2x2y_S_D2y_S_bc+ABX*I_ERI_Fx2y_S_D2y_S_bc;
  Double I_ERI_Fxyz_Px_D2y_S_bc = I_ERI_G2xyz_S_D2y_S_bc+ABX*I_ERI_Fxyz_S_D2y_S_bc;
  Double I_ERI_Fx2z_Px_D2y_S_bc = I_ERI_G2x2z_S_D2y_S_bc+ABX*I_ERI_Fx2z_S_D2y_S_bc;
  Double I_ERI_F3y_Px_D2y_S_bc = I_ERI_Gx3y_S_D2y_S_bc+ABX*I_ERI_F3y_S_D2y_S_bc;
  Double I_ERI_F2yz_Px_D2y_S_bc = I_ERI_Gx2yz_S_D2y_S_bc+ABX*I_ERI_F2yz_S_D2y_S_bc;
  Double I_ERI_Fy2z_Px_D2y_S_bc = I_ERI_Gxy2z_S_D2y_S_bc+ABX*I_ERI_Fy2z_S_D2y_S_bc;
  Double I_ERI_F3z_Px_D2y_S_bc = I_ERI_Gx3z_S_D2y_S_bc+ABX*I_ERI_F3z_S_D2y_S_bc;
  Double I_ERI_F2xy_Py_D2y_S_bc = I_ERI_G2x2y_S_D2y_S_bc+ABY*I_ERI_F2xy_S_D2y_S_bc;
  Double I_ERI_F2xz_Py_D2y_S_bc = I_ERI_G2xyz_S_D2y_S_bc+ABY*I_ERI_F2xz_S_D2y_S_bc;
  Double I_ERI_Fx2y_Py_D2y_S_bc = I_ERI_Gx3y_S_D2y_S_bc+ABY*I_ERI_Fx2y_S_D2y_S_bc;
  Double I_ERI_Fxyz_Py_D2y_S_bc = I_ERI_Gx2yz_S_D2y_S_bc+ABY*I_ERI_Fxyz_S_D2y_S_bc;
  Double I_ERI_Fx2z_Py_D2y_S_bc = I_ERI_Gxy2z_S_D2y_S_bc+ABY*I_ERI_Fx2z_S_D2y_S_bc;
  Double I_ERI_F3y_Py_D2y_S_bc = I_ERI_G4y_S_D2y_S_bc+ABY*I_ERI_F3y_S_D2y_S_bc;
  Double I_ERI_F2yz_Py_D2y_S_bc = I_ERI_G3yz_S_D2y_S_bc+ABY*I_ERI_F2yz_S_D2y_S_bc;
  Double I_ERI_Fy2z_Py_D2y_S_bc = I_ERI_G2y2z_S_D2y_S_bc+ABY*I_ERI_Fy2z_S_D2y_S_bc;
  Double I_ERI_F3z_Py_D2y_S_bc = I_ERI_Gy3z_S_D2y_S_bc+ABY*I_ERI_F3z_S_D2y_S_bc;
  Double I_ERI_F2xz_Pz_D2y_S_bc = I_ERI_G2x2z_S_D2y_S_bc+ABZ*I_ERI_F2xz_S_D2y_S_bc;
  Double I_ERI_Fxyz_Pz_D2y_S_bc = I_ERI_Gxy2z_S_D2y_S_bc+ABZ*I_ERI_Fxyz_S_D2y_S_bc;
  Double I_ERI_Fx2z_Pz_D2y_S_bc = I_ERI_Gx3z_S_D2y_S_bc+ABZ*I_ERI_Fx2z_S_D2y_S_bc;
  Double I_ERI_F2yz_Pz_D2y_S_bc = I_ERI_G2y2z_S_D2y_S_bc+ABZ*I_ERI_F2yz_S_D2y_S_bc;
  Double I_ERI_Fy2z_Pz_D2y_S_bc = I_ERI_Gy3z_S_D2y_S_bc+ABZ*I_ERI_Fy2z_S_D2y_S_bc;
  Double I_ERI_F3z_Pz_D2y_S_bc = I_ERI_G4z_S_D2y_S_bc+ABZ*I_ERI_F3z_S_D2y_S_bc;
  Double I_ERI_F3x_Px_Dyz_S_bc = I_ERI_G4x_S_Dyz_S_bc+ABX*I_ERI_F3x_S_Dyz_S_bc;
  Double I_ERI_F2xy_Px_Dyz_S_bc = I_ERI_G3xy_S_Dyz_S_bc+ABX*I_ERI_F2xy_S_Dyz_S_bc;
  Double I_ERI_F2xz_Px_Dyz_S_bc = I_ERI_G3xz_S_Dyz_S_bc+ABX*I_ERI_F2xz_S_Dyz_S_bc;
  Double I_ERI_Fx2y_Px_Dyz_S_bc = I_ERI_G2x2y_S_Dyz_S_bc+ABX*I_ERI_Fx2y_S_Dyz_S_bc;
  Double I_ERI_Fxyz_Px_Dyz_S_bc = I_ERI_G2xyz_S_Dyz_S_bc+ABX*I_ERI_Fxyz_S_Dyz_S_bc;
  Double I_ERI_Fx2z_Px_Dyz_S_bc = I_ERI_G2x2z_S_Dyz_S_bc+ABX*I_ERI_Fx2z_S_Dyz_S_bc;
  Double I_ERI_F3y_Px_Dyz_S_bc = I_ERI_Gx3y_S_Dyz_S_bc+ABX*I_ERI_F3y_S_Dyz_S_bc;
  Double I_ERI_F2yz_Px_Dyz_S_bc = I_ERI_Gx2yz_S_Dyz_S_bc+ABX*I_ERI_F2yz_S_Dyz_S_bc;
  Double I_ERI_Fy2z_Px_Dyz_S_bc = I_ERI_Gxy2z_S_Dyz_S_bc+ABX*I_ERI_Fy2z_S_Dyz_S_bc;
  Double I_ERI_F3z_Px_Dyz_S_bc = I_ERI_Gx3z_S_Dyz_S_bc+ABX*I_ERI_F3z_S_Dyz_S_bc;
  Double I_ERI_F2xy_Py_Dyz_S_bc = I_ERI_G2x2y_S_Dyz_S_bc+ABY*I_ERI_F2xy_S_Dyz_S_bc;
  Double I_ERI_F2xz_Py_Dyz_S_bc = I_ERI_G2xyz_S_Dyz_S_bc+ABY*I_ERI_F2xz_S_Dyz_S_bc;
  Double I_ERI_Fx2y_Py_Dyz_S_bc = I_ERI_Gx3y_S_Dyz_S_bc+ABY*I_ERI_Fx2y_S_Dyz_S_bc;
  Double I_ERI_Fxyz_Py_Dyz_S_bc = I_ERI_Gx2yz_S_Dyz_S_bc+ABY*I_ERI_Fxyz_S_Dyz_S_bc;
  Double I_ERI_Fx2z_Py_Dyz_S_bc = I_ERI_Gxy2z_S_Dyz_S_bc+ABY*I_ERI_Fx2z_S_Dyz_S_bc;
  Double I_ERI_F3y_Py_Dyz_S_bc = I_ERI_G4y_S_Dyz_S_bc+ABY*I_ERI_F3y_S_Dyz_S_bc;
  Double I_ERI_F2yz_Py_Dyz_S_bc = I_ERI_G3yz_S_Dyz_S_bc+ABY*I_ERI_F2yz_S_Dyz_S_bc;
  Double I_ERI_Fy2z_Py_Dyz_S_bc = I_ERI_G2y2z_S_Dyz_S_bc+ABY*I_ERI_Fy2z_S_Dyz_S_bc;
  Double I_ERI_F3z_Py_Dyz_S_bc = I_ERI_Gy3z_S_Dyz_S_bc+ABY*I_ERI_F3z_S_Dyz_S_bc;
  Double I_ERI_F2xz_Pz_Dyz_S_bc = I_ERI_G2x2z_S_Dyz_S_bc+ABZ*I_ERI_F2xz_S_Dyz_S_bc;
  Double I_ERI_Fxyz_Pz_Dyz_S_bc = I_ERI_Gxy2z_S_Dyz_S_bc+ABZ*I_ERI_Fxyz_S_Dyz_S_bc;
  Double I_ERI_Fx2z_Pz_Dyz_S_bc = I_ERI_Gx3z_S_Dyz_S_bc+ABZ*I_ERI_Fx2z_S_Dyz_S_bc;
  Double I_ERI_F2yz_Pz_Dyz_S_bc = I_ERI_G2y2z_S_Dyz_S_bc+ABZ*I_ERI_F2yz_S_Dyz_S_bc;
  Double I_ERI_Fy2z_Pz_Dyz_S_bc = I_ERI_Gy3z_S_Dyz_S_bc+ABZ*I_ERI_Fy2z_S_Dyz_S_bc;
  Double I_ERI_F3z_Pz_Dyz_S_bc = I_ERI_G4z_S_Dyz_S_bc+ABZ*I_ERI_F3z_S_Dyz_S_bc;
  Double I_ERI_F3x_Px_D2z_S_bc = I_ERI_G4x_S_D2z_S_bc+ABX*I_ERI_F3x_S_D2z_S_bc;
  Double I_ERI_F2xy_Px_D2z_S_bc = I_ERI_G3xy_S_D2z_S_bc+ABX*I_ERI_F2xy_S_D2z_S_bc;
  Double I_ERI_F2xz_Px_D2z_S_bc = I_ERI_G3xz_S_D2z_S_bc+ABX*I_ERI_F2xz_S_D2z_S_bc;
  Double I_ERI_Fx2y_Px_D2z_S_bc = I_ERI_G2x2y_S_D2z_S_bc+ABX*I_ERI_Fx2y_S_D2z_S_bc;
  Double I_ERI_Fxyz_Px_D2z_S_bc = I_ERI_G2xyz_S_D2z_S_bc+ABX*I_ERI_Fxyz_S_D2z_S_bc;
  Double I_ERI_Fx2z_Px_D2z_S_bc = I_ERI_G2x2z_S_D2z_S_bc+ABX*I_ERI_Fx2z_S_D2z_S_bc;
  Double I_ERI_F3y_Px_D2z_S_bc = I_ERI_Gx3y_S_D2z_S_bc+ABX*I_ERI_F3y_S_D2z_S_bc;
  Double I_ERI_F2yz_Px_D2z_S_bc = I_ERI_Gx2yz_S_D2z_S_bc+ABX*I_ERI_F2yz_S_D2z_S_bc;
  Double I_ERI_Fy2z_Px_D2z_S_bc = I_ERI_Gxy2z_S_D2z_S_bc+ABX*I_ERI_Fy2z_S_D2z_S_bc;
  Double I_ERI_F3z_Px_D2z_S_bc = I_ERI_Gx3z_S_D2z_S_bc+ABX*I_ERI_F3z_S_D2z_S_bc;
  Double I_ERI_F2xy_Py_D2z_S_bc = I_ERI_G2x2y_S_D2z_S_bc+ABY*I_ERI_F2xy_S_D2z_S_bc;
  Double I_ERI_F2xz_Py_D2z_S_bc = I_ERI_G2xyz_S_D2z_S_bc+ABY*I_ERI_F2xz_S_D2z_S_bc;
  Double I_ERI_Fx2y_Py_D2z_S_bc = I_ERI_Gx3y_S_D2z_S_bc+ABY*I_ERI_Fx2y_S_D2z_S_bc;
  Double I_ERI_Fxyz_Py_D2z_S_bc = I_ERI_Gx2yz_S_D2z_S_bc+ABY*I_ERI_Fxyz_S_D2z_S_bc;
  Double I_ERI_Fx2z_Py_D2z_S_bc = I_ERI_Gxy2z_S_D2z_S_bc+ABY*I_ERI_Fx2z_S_D2z_S_bc;
  Double I_ERI_F3y_Py_D2z_S_bc = I_ERI_G4y_S_D2z_S_bc+ABY*I_ERI_F3y_S_D2z_S_bc;
  Double I_ERI_F2yz_Py_D2z_S_bc = I_ERI_G3yz_S_D2z_S_bc+ABY*I_ERI_F2yz_S_D2z_S_bc;
  Double I_ERI_Fy2z_Py_D2z_S_bc = I_ERI_G2y2z_S_D2z_S_bc+ABY*I_ERI_Fy2z_S_D2z_S_bc;
  Double I_ERI_F3z_Py_D2z_S_bc = I_ERI_Gy3z_S_D2z_S_bc+ABY*I_ERI_F3z_S_D2z_S_bc;
  Double I_ERI_F2xz_Pz_D2z_S_bc = I_ERI_G2x2z_S_D2z_S_bc+ABZ*I_ERI_F2xz_S_D2z_S_bc;
  Double I_ERI_Fxyz_Pz_D2z_S_bc = I_ERI_Gxy2z_S_D2z_S_bc+ABZ*I_ERI_Fxyz_S_D2z_S_bc;
  Double I_ERI_Fx2z_Pz_D2z_S_bc = I_ERI_Gx3z_S_D2z_S_bc+ABZ*I_ERI_Fx2z_S_D2z_S_bc;
  Double I_ERI_F2yz_Pz_D2z_S_bc = I_ERI_G2y2z_S_D2z_S_bc+ABZ*I_ERI_F2yz_S_D2z_S_bc;
  Double I_ERI_Fy2z_Pz_D2z_S_bc = I_ERI_Gy3z_S_D2z_S_bc+ABZ*I_ERI_Fy2z_S_D2z_S_bc;
  Double I_ERI_F3z_Pz_D2z_S_bc = I_ERI_G4z_S_D2z_S_bc+ABZ*I_ERI_F3z_S_D2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_D_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_P_D_S_bc
   ************************************************************/
  Double I_ERI_D2x_D2x_D2x_S_bc = I_ERI_F3x_Px_D2x_S_bc+ABX*I_ERI_D2x_Px_D2x_S_bc;
  Double I_ERI_Dxy_D2x_D2x_S_bc = I_ERI_F2xy_Px_D2x_S_bc+ABX*I_ERI_Dxy_Px_D2x_S_bc;
  Double I_ERI_Dxz_D2x_D2x_S_bc = I_ERI_F2xz_Px_D2x_S_bc+ABX*I_ERI_Dxz_Px_D2x_S_bc;
  Double I_ERI_D2y_D2x_D2x_S_bc = I_ERI_Fx2y_Px_D2x_S_bc+ABX*I_ERI_D2y_Px_D2x_S_bc;
  Double I_ERI_Dyz_D2x_D2x_S_bc = I_ERI_Fxyz_Px_D2x_S_bc+ABX*I_ERI_Dyz_Px_D2x_S_bc;
  Double I_ERI_D2z_D2x_D2x_S_bc = I_ERI_Fx2z_Px_D2x_S_bc+ABX*I_ERI_D2z_Px_D2x_S_bc;
  Double I_ERI_D2x_Dxy_D2x_S_bc = I_ERI_F2xy_Px_D2x_S_bc+ABY*I_ERI_D2x_Px_D2x_S_bc;
  Double I_ERI_Dxy_Dxy_D2x_S_bc = I_ERI_Fx2y_Px_D2x_S_bc+ABY*I_ERI_Dxy_Px_D2x_S_bc;
  Double I_ERI_Dxz_Dxy_D2x_S_bc = I_ERI_Fxyz_Px_D2x_S_bc+ABY*I_ERI_Dxz_Px_D2x_S_bc;
  Double I_ERI_D2y_Dxy_D2x_S_bc = I_ERI_F3y_Px_D2x_S_bc+ABY*I_ERI_D2y_Px_D2x_S_bc;
  Double I_ERI_Dyz_Dxy_D2x_S_bc = I_ERI_F2yz_Px_D2x_S_bc+ABY*I_ERI_Dyz_Px_D2x_S_bc;
  Double I_ERI_D2z_Dxy_D2x_S_bc = I_ERI_Fy2z_Px_D2x_S_bc+ABY*I_ERI_D2z_Px_D2x_S_bc;
  Double I_ERI_D2x_Dxz_D2x_S_bc = I_ERI_F2xz_Px_D2x_S_bc+ABZ*I_ERI_D2x_Px_D2x_S_bc;
  Double I_ERI_Dxy_Dxz_D2x_S_bc = I_ERI_Fxyz_Px_D2x_S_bc+ABZ*I_ERI_Dxy_Px_D2x_S_bc;
  Double I_ERI_Dxz_Dxz_D2x_S_bc = I_ERI_Fx2z_Px_D2x_S_bc+ABZ*I_ERI_Dxz_Px_D2x_S_bc;
  Double I_ERI_D2y_Dxz_D2x_S_bc = I_ERI_F2yz_Px_D2x_S_bc+ABZ*I_ERI_D2y_Px_D2x_S_bc;
  Double I_ERI_Dyz_Dxz_D2x_S_bc = I_ERI_Fy2z_Px_D2x_S_bc+ABZ*I_ERI_Dyz_Px_D2x_S_bc;
  Double I_ERI_D2z_Dxz_D2x_S_bc = I_ERI_F3z_Px_D2x_S_bc+ABZ*I_ERI_D2z_Px_D2x_S_bc;
  Double I_ERI_D2x_D2y_D2x_S_bc = I_ERI_F2xy_Py_D2x_S_bc+ABY*I_ERI_D2x_Py_D2x_S_bc;
  Double I_ERI_Dxy_D2y_D2x_S_bc = I_ERI_Fx2y_Py_D2x_S_bc+ABY*I_ERI_Dxy_Py_D2x_S_bc;
  Double I_ERI_Dxz_D2y_D2x_S_bc = I_ERI_Fxyz_Py_D2x_S_bc+ABY*I_ERI_Dxz_Py_D2x_S_bc;
  Double I_ERI_D2y_D2y_D2x_S_bc = I_ERI_F3y_Py_D2x_S_bc+ABY*I_ERI_D2y_Py_D2x_S_bc;
  Double I_ERI_Dyz_D2y_D2x_S_bc = I_ERI_F2yz_Py_D2x_S_bc+ABY*I_ERI_Dyz_Py_D2x_S_bc;
  Double I_ERI_D2z_D2y_D2x_S_bc = I_ERI_Fy2z_Py_D2x_S_bc+ABY*I_ERI_D2z_Py_D2x_S_bc;
  Double I_ERI_D2x_Dyz_D2x_S_bc = I_ERI_F2xz_Py_D2x_S_bc+ABZ*I_ERI_D2x_Py_D2x_S_bc;
  Double I_ERI_Dxy_Dyz_D2x_S_bc = I_ERI_Fxyz_Py_D2x_S_bc+ABZ*I_ERI_Dxy_Py_D2x_S_bc;
  Double I_ERI_Dxz_Dyz_D2x_S_bc = I_ERI_Fx2z_Py_D2x_S_bc+ABZ*I_ERI_Dxz_Py_D2x_S_bc;
  Double I_ERI_D2y_Dyz_D2x_S_bc = I_ERI_F2yz_Py_D2x_S_bc+ABZ*I_ERI_D2y_Py_D2x_S_bc;
  Double I_ERI_Dyz_Dyz_D2x_S_bc = I_ERI_Fy2z_Py_D2x_S_bc+ABZ*I_ERI_Dyz_Py_D2x_S_bc;
  Double I_ERI_D2z_Dyz_D2x_S_bc = I_ERI_F3z_Py_D2x_S_bc+ABZ*I_ERI_D2z_Py_D2x_S_bc;
  Double I_ERI_D2x_D2z_D2x_S_bc = I_ERI_F2xz_Pz_D2x_S_bc+ABZ*I_ERI_D2x_Pz_D2x_S_bc;
  Double I_ERI_Dxy_D2z_D2x_S_bc = I_ERI_Fxyz_Pz_D2x_S_bc+ABZ*I_ERI_Dxy_Pz_D2x_S_bc;
  Double I_ERI_Dxz_D2z_D2x_S_bc = I_ERI_Fx2z_Pz_D2x_S_bc+ABZ*I_ERI_Dxz_Pz_D2x_S_bc;
  Double I_ERI_D2y_D2z_D2x_S_bc = I_ERI_F2yz_Pz_D2x_S_bc+ABZ*I_ERI_D2y_Pz_D2x_S_bc;
  Double I_ERI_Dyz_D2z_D2x_S_bc = I_ERI_Fy2z_Pz_D2x_S_bc+ABZ*I_ERI_Dyz_Pz_D2x_S_bc;
  Double I_ERI_D2z_D2z_D2x_S_bc = I_ERI_F3z_Pz_D2x_S_bc+ABZ*I_ERI_D2z_Pz_D2x_S_bc;
  Double I_ERI_D2x_D2x_Dxy_S_bc = I_ERI_F3x_Px_Dxy_S_bc+ABX*I_ERI_D2x_Px_Dxy_S_bc;
  Double I_ERI_Dxy_D2x_Dxy_S_bc = I_ERI_F2xy_Px_Dxy_S_bc+ABX*I_ERI_Dxy_Px_Dxy_S_bc;
  Double I_ERI_Dxz_D2x_Dxy_S_bc = I_ERI_F2xz_Px_Dxy_S_bc+ABX*I_ERI_Dxz_Px_Dxy_S_bc;
  Double I_ERI_D2y_D2x_Dxy_S_bc = I_ERI_Fx2y_Px_Dxy_S_bc+ABX*I_ERI_D2y_Px_Dxy_S_bc;
  Double I_ERI_Dyz_D2x_Dxy_S_bc = I_ERI_Fxyz_Px_Dxy_S_bc+ABX*I_ERI_Dyz_Px_Dxy_S_bc;
  Double I_ERI_D2z_D2x_Dxy_S_bc = I_ERI_Fx2z_Px_Dxy_S_bc+ABX*I_ERI_D2z_Px_Dxy_S_bc;
  Double I_ERI_D2x_Dxy_Dxy_S_bc = I_ERI_F2xy_Px_Dxy_S_bc+ABY*I_ERI_D2x_Px_Dxy_S_bc;
  Double I_ERI_Dxy_Dxy_Dxy_S_bc = I_ERI_Fx2y_Px_Dxy_S_bc+ABY*I_ERI_Dxy_Px_Dxy_S_bc;
  Double I_ERI_Dxz_Dxy_Dxy_S_bc = I_ERI_Fxyz_Px_Dxy_S_bc+ABY*I_ERI_Dxz_Px_Dxy_S_bc;
  Double I_ERI_D2y_Dxy_Dxy_S_bc = I_ERI_F3y_Px_Dxy_S_bc+ABY*I_ERI_D2y_Px_Dxy_S_bc;
  Double I_ERI_Dyz_Dxy_Dxy_S_bc = I_ERI_F2yz_Px_Dxy_S_bc+ABY*I_ERI_Dyz_Px_Dxy_S_bc;
  Double I_ERI_D2z_Dxy_Dxy_S_bc = I_ERI_Fy2z_Px_Dxy_S_bc+ABY*I_ERI_D2z_Px_Dxy_S_bc;
  Double I_ERI_D2x_Dxz_Dxy_S_bc = I_ERI_F2xz_Px_Dxy_S_bc+ABZ*I_ERI_D2x_Px_Dxy_S_bc;
  Double I_ERI_Dxy_Dxz_Dxy_S_bc = I_ERI_Fxyz_Px_Dxy_S_bc+ABZ*I_ERI_Dxy_Px_Dxy_S_bc;
  Double I_ERI_Dxz_Dxz_Dxy_S_bc = I_ERI_Fx2z_Px_Dxy_S_bc+ABZ*I_ERI_Dxz_Px_Dxy_S_bc;
  Double I_ERI_D2y_Dxz_Dxy_S_bc = I_ERI_F2yz_Px_Dxy_S_bc+ABZ*I_ERI_D2y_Px_Dxy_S_bc;
  Double I_ERI_Dyz_Dxz_Dxy_S_bc = I_ERI_Fy2z_Px_Dxy_S_bc+ABZ*I_ERI_Dyz_Px_Dxy_S_bc;
  Double I_ERI_D2z_Dxz_Dxy_S_bc = I_ERI_F3z_Px_Dxy_S_bc+ABZ*I_ERI_D2z_Px_Dxy_S_bc;
  Double I_ERI_D2x_D2y_Dxy_S_bc = I_ERI_F2xy_Py_Dxy_S_bc+ABY*I_ERI_D2x_Py_Dxy_S_bc;
  Double I_ERI_Dxy_D2y_Dxy_S_bc = I_ERI_Fx2y_Py_Dxy_S_bc+ABY*I_ERI_Dxy_Py_Dxy_S_bc;
  Double I_ERI_Dxz_D2y_Dxy_S_bc = I_ERI_Fxyz_Py_Dxy_S_bc+ABY*I_ERI_Dxz_Py_Dxy_S_bc;
  Double I_ERI_D2y_D2y_Dxy_S_bc = I_ERI_F3y_Py_Dxy_S_bc+ABY*I_ERI_D2y_Py_Dxy_S_bc;
  Double I_ERI_Dyz_D2y_Dxy_S_bc = I_ERI_F2yz_Py_Dxy_S_bc+ABY*I_ERI_Dyz_Py_Dxy_S_bc;
  Double I_ERI_D2z_D2y_Dxy_S_bc = I_ERI_Fy2z_Py_Dxy_S_bc+ABY*I_ERI_D2z_Py_Dxy_S_bc;
  Double I_ERI_D2x_Dyz_Dxy_S_bc = I_ERI_F2xz_Py_Dxy_S_bc+ABZ*I_ERI_D2x_Py_Dxy_S_bc;
  Double I_ERI_Dxy_Dyz_Dxy_S_bc = I_ERI_Fxyz_Py_Dxy_S_bc+ABZ*I_ERI_Dxy_Py_Dxy_S_bc;
  Double I_ERI_Dxz_Dyz_Dxy_S_bc = I_ERI_Fx2z_Py_Dxy_S_bc+ABZ*I_ERI_Dxz_Py_Dxy_S_bc;
  Double I_ERI_D2y_Dyz_Dxy_S_bc = I_ERI_F2yz_Py_Dxy_S_bc+ABZ*I_ERI_D2y_Py_Dxy_S_bc;
  Double I_ERI_Dyz_Dyz_Dxy_S_bc = I_ERI_Fy2z_Py_Dxy_S_bc+ABZ*I_ERI_Dyz_Py_Dxy_S_bc;
  Double I_ERI_D2z_Dyz_Dxy_S_bc = I_ERI_F3z_Py_Dxy_S_bc+ABZ*I_ERI_D2z_Py_Dxy_S_bc;
  Double I_ERI_D2x_D2z_Dxy_S_bc = I_ERI_F2xz_Pz_Dxy_S_bc+ABZ*I_ERI_D2x_Pz_Dxy_S_bc;
  Double I_ERI_Dxy_D2z_Dxy_S_bc = I_ERI_Fxyz_Pz_Dxy_S_bc+ABZ*I_ERI_Dxy_Pz_Dxy_S_bc;
  Double I_ERI_Dxz_D2z_Dxy_S_bc = I_ERI_Fx2z_Pz_Dxy_S_bc+ABZ*I_ERI_Dxz_Pz_Dxy_S_bc;
  Double I_ERI_D2y_D2z_Dxy_S_bc = I_ERI_F2yz_Pz_Dxy_S_bc+ABZ*I_ERI_D2y_Pz_Dxy_S_bc;
  Double I_ERI_Dyz_D2z_Dxy_S_bc = I_ERI_Fy2z_Pz_Dxy_S_bc+ABZ*I_ERI_Dyz_Pz_Dxy_S_bc;
  Double I_ERI_D2z_D2z_Dxy_S_bc = I_ERI_F3z_Pz_Dxy_S_bc+ABZ*I_ERI_D2z_Pz_Dxy_S_bc;
  Double I_ERI_D2x_D2x_Dxz_S_bc = I_ERI_F3x_Px_Dxz_S_bc+ABX*I_ERI_D2x_Px_Dxz_S_bc;
  Double I_ERI_Dxy_D2x_Dxz_S_bc = I_ERI_F2xy_Px_Dxz_S_bc+ABX*I_ERI_Dxy_Px_Dxz_S_bc;
  Double I_ERI_Dxz_D2x_Dxz_S_bc = I_ERI_F2xz_Px_Dxz_S_bc+ABX*I_ERI_Dxz_Px_Dxz_S_bc;
  Double I_ERI_D2y_D2x_Dxz_S_bc = I_ERI_Fx2y_Px_Dxz_S_bc+ABX*I_ERI_D2y_Px_Dxz_S_bc;
  Double I_ERI_Dyz_D2x_Dxz_S_bc = I_ERI_Fxyz_Px_Dxz_S_bc+ABX*I_ERI_Dyz_Px_Dxz_S_bc;
  Double I_ERI_D2z_D2x_Dxz_S_bc = I_ERI_Fx2z_Px_Dxz_S_bc+ABX*I_ERI_D2z_Px_Dxz_S_bc;
  Double I_ERI_D2x_Dxy_Dxz_S_bc = I_ERI_F2xy_Px_Dxz_S_bc+ABY*I_ERI_D2x_Px_Dxz_S_bc;
  Double I_ERI_Dxy_Dxy_Dxz_S_bc = I_ERI_Fx2y_Px_Dxz_S_bc+ABY*I_ERI_Dxy_Px_Dxz_S_bc;
  Double I_ERI_Dxz_Dxy_Dxz_S_bc = I_ERI_Fxyz_Px_Dxz_S_bc+ABY*I_ERI_Dxz_Px_Dxz_S_bc;
  Double I_ERI_D2y_Dxy_Dxz_S_bc = I_ERI_F3y_Px_Dxz_S_bc+ABY*I_ERI_D2y_Px_Dxz_S_bc;
  Double I_ERI_Dyz_Dxy_Dxz_S_bc = I_ERI_F2yz_Px_Dxz_S_bc+ABY*I_ERI_Dyz_Px_Dxz_S_bc;
  Double I_ERI_D2z_Dxy_Dxz_S_bc = I_ERI_Fy2z_Px_Dxz_S_bc+ABY*I_ERI_D2z_Px_Dxz_S_bc;
  Double I_ERI_D2x_Dxz_Dxz_S_bc = I_ERI_F2xz_Px_Dxz_S_bc+ABZ*I_ERI_D2x_Px_Dxz_S_bc;
  Double I_ERI_Dxy_Dxz_Dxz_S_bc = I_ERI_Fxyz_Px_Dxz_S_bc+ABZ*I_ERI_Dxy_Px_Dxz_S_bc;
  Double I_ERI_Dxz_Dxz_Dxz_S_bc = I_ERI_Fx2z_Px_Dxz_S_bc+ABZ*I_ERI_Dxz_Px_Dxz_S_bc;
  Double I_ERI_D2y_Dxz_Dxz_S_bc = I_ERI_F2yz_Px_Dxz_S_bc+ABZ*I_ERI_D2y_Px_Dxz_S_bc;
  Double I_ERI_Dyz_Dxz_Dxz_S_bc = I_ERI_Fy2z_Px_Dxz_S_bc+ABZ*I_ERI_Dyz_Px_Dxz_S_bc;
  Double I_ERI_D2z_Dxz_Dxz_S_bc = I_ERI_F3z_Px_Dxz_S_bc+ABZ*I_ERI_D2z_Px_Dxz_S_bc;
  Double I_ERI_D2x_D2y_Dxz_S_bc = I_ERI_F2xy_Py_Dxz_S_bc+ABY*I_ERI_D2x_Py_Dxz_S_bc;
  Double I_ERI_Dxy_D2y_Dxz_S_bc = I_ERI_Fx2y_Py_Dxz_S_bc+ABY*I_ERI_Dxy_Py_Dxz_S_bc;
  Double I_ERI_Dxz_D2y_Dxz_S_bc = I_ERI_Fxyz_Py_Dxz_S_bc+ABY*I_ERI_Dxz_Py_Dxz_S_bc;
  Double I_ERI_D2y_D2y_Dxz_S_bc = I_ERI_F3y_Py_Dxz_S_bc+ABY*I_ERI_D2y_Py_Dxz_S_bc;
  Double I_ERI_Dyz_D2y_Dxz_S_bc = I_ERI_F2yz_Py_Dxz_S_bc+ABY*I_ERI_Dyz_Py_Dxz_S_bc;
  Double I_ERI_D2z_D2y_Dxz_S_bc = I_ERI_Fy2z_Py_Dxz_S_bc+ABY*I_ERI_D2z_Py_Dxz_S_bc;
  Double I_ERI_D2x_Dyz_Dxz_S_bc = I_ERI_F2xz_Py_Dxz_S_bc+ABZ*I_ERI_D2x_Py_Dxz_S_bc;
  Double I_ERI_Dxy_Dyz_Dxz_S_bc = I_ERI_Fxyz_Py_Dxz_S_bc+ABZ*I_ERI_Dxy_Py_Dxz_S_bc;
  Double I_ERI_Dxz_Dyz_Dxz_S_bc = I_ERI_Fx2z_Py_Dxz_S_bc+ABZ*I_ERI_Dxz_Py_Dxz_S_bc;
  Double I_ERI_D2y_Dyz_Dxz_S_bc = I_ERI_F2yz_Py_Dxz_S_bc+ABZ*I_ERI_D2y_Py_Dxz_S_bc;
  Double I_ERI_Dyz_Dyz_Dxz_S_bc = I_ERI_Fy2z_Py_Dxz_S_bc+ABZ*I_ERI_Dyz_Py_Dxz_S_bc;
  Double I_ERI_D2z_Dyz_Dxz_S_bc = I_ERI_F3z_Py_Dxz_S_bc+ABZ*I_ERI_D2z_Py_Dxz_S_bc;
  Double I_ERI_D2x_D2z_Dxz_S_bc = I_ERI_F2xz_Pz_Dxz_S_bc+ABZ*I_ERI_D2x_Pz_Dxz_S_bc;
  Double I_ERI_Dxy_D2z_Dxz_S_bc = I_ERI_Fxyz_Pz_Dxz_S_bc+ABZ*I_ERI_Dxy_Pz_Dxz_S_bc;
  Double I_ERI_Dxz_D2z_Dxz_S_bc = I_ERI_Fx2z_Pz_Dxz_S_bc+ABZ*I_ERI_Dxz_Pz_Dxz_S_bc;
  Double I_ERI_D2y_D2z_Dxz_S_bc = I_ERI_F2yz_Pz_Dxz_S_bc+ABZ*I_ERI_D2y_Pz_Dxz_S_bc;
  Double I_ERI_Dyz_D2z_Dxz_S_bc = I_ERI_Fy2z_Pz_Dxz_S_bc+ABZ*I_ERI_Dyz_Pz_Dxz_S_bc;
  Double I_ERI_D2z_D2z_Dxz_S_bc = I_ERI_F3z_Pz_Dxz_S_bc+ABZ*I_ERI_D2z_Pz_Dxz_S_bc;
  Double I_ERI_D2x_D2x_D2y_S_bc = I_ERI_F3x_Px_D2y_S_bc+ABX*I_ERI_D2x_Px_D2y_S_bc;
  Double I_ERI_Dxy_D2x_D2y_S_bc = I_ERI_F2xy_Px_D2y_S_bc+ABX*I_ERI_Dxy_Px_D2y_S_bc;
  Double I_ERI_Dxz_D2x_D2y_S_bc = I_ERI_F2xz_Px_D2y_S_bc+ABX*I_ERI_Dxz_Px_D2y_S_bc;
  Double I_ERI_D2y_D2x_D2y_S_bc = I_ERI_Fx2y_Px_D2y_S_bc+ABX*I_ERI_D2y_Px_D2y_S_bc;
  Double I_ERI_Dyz_D2x_D2y_S_bc = I_ERI_Fxyz_Px_D2y_S_bc+ABX*I_ERI_Dyz_Px_D2y_S_bc;
  Double I_ERI_D2z_D2x_D2y_S_bc = I_ERI_Fx2z_Px_D2y_S_bc+ABX*I_ERI_D2z_Px_D2y_S_bc;
  Double I_ERI_D2x_Dxy_D2y_S_bc = I_ERI_F2xy_Px_D2y_S_bc+ABY*I_ERI_D2x_Px_D2y_S_bc;
  Double I_ERI_Dxy_Dxy_D2y_S_bc = I_ERI_Fx2y_Px_D2y_S_bc+ABY*I_ERI_Dxy_Px_D2y_S_bc;
  Double I_ERI_Dxz_Dxy_D2y_S_bc = I_ERI_Fxyz_Px_D2y_S_bc+ABY*I_ERI_Dxz_Px_D2y_S_bc;
  Double I_ERI_D2y_Dxy_D2y_S_bc = I_ERI_F3y_Px_D2y_S_bc+ABY*I_ERI_D2y_Px_D2y_S_bc;
  Double I_ERI_Dyz_Dxy_D2y_S_bc = I_ERI_F2yz_Px_D2y_S_bc+ABY*I_ERI_Dyz_Px_D2y_S_bc;
  Double I_ERI_D2z_Dxy_D2y_S_bc = I_ERI_Fy2z_Px_D2y_S_bc+ABY*I_ERI_D2z_Px_D2y_S_bc;
  Double I_ERI_D2x_Dxz_D2y_S_bc = I_ERI_F2xz_Px_D2y_S_bc+ABZ*I_ERI_D2x_Px_D2y_S_bc;
  Double I_ERI_Dxy_Dxz_D2y_S_bc = I_ERI_Fxyz_Px_D2y_S_bc+ABZ*I_ERI_Dxy_Px_D2y_S_bc;
  Double I_ERI_Dxz_Dxz_D2y_S_bc = I_ERI_Fx2z_Px_D2y_S_bc+ABZ*I_ERI_Dxz_Px_D2y_S_bc;
  Double I_ERI_D2y_Dxz_D2y_S_bc = I_ERI_F2yz_Px_D2y_S_bc+ABZ*I_ERI_D2y_Px_D2y_S_bc;
  Double I_ERI_Dyz_Dxz_D2y_S_bc = I_ERI_Fy2z_Px_D2y_S_bc+ABZ*I_ERI_Dyz_Px_D2y_S_bc;
  Double I_ERI_D2z_Dxz_D2y_S_bc = I_ERI_F3z_Px_D2y_S_bc+ABZ*I_ERI_D2z_Px_D2y_S_bc;
  Double I_ERI_D2x_D2y_D2y_S_bc = I_ERI_F2xy_Py_D2y_S_bc+ABY*I_ERI_D2x_Py_D2y_S_bc;
  Double I_ERI_Dxy_D2y_D2y_S_bc = I_ERI_Fx2y_Py_D2y_S_bc+ABY*I_ERI_Dxy_Py_D2y_S_bc;
  Double I_ERI_Dxz_D2y_D2y_S_bc = I_ERI_Fxyz_Py_D2y_S_bc+ABY*I_ERI_Dxz_Py_D2y_S_bc;
  Double I_ERI_D2y_D2y_D2y_S_bc = I_ERI_F3y_Py_D2y_S_bc+ABY*I_ERI_D2y_Py_D2y_S_bc;
  Double I_ERI_Dyz_D2y_D2y_S_bc = I_ERI_F2yz_Py_D2y_S_bc+ABY*I_ERI_Dyz_Py_D2y_S_bc;
  Double I_ERI_D2z_D2y_D2y_S_bc = I_ERI_Fy2z_Py_D2y_S_bc+ABY*I_ERI_D2z_Py_D2y_S_bc;
  Double I_ERI_D2x_Dyz_D2y_S_bc = I_ERI_F2xz_Py_D2y_S_bc+ABZ*I_ERI_D2x_Py_D2y_S_bc;
  Double I_ERI_Dxy_Dyz_D2y_S_bc = I_ERI_Fxyz_Py_D2y_S_bc+ABZ*I_ERI_Dxy_Py_D2y_S_bc;
  Double I_ERI_Dxz_Dyz_D2y_S_bc = I_ERI_Fx2z_Py_D2y_S_bc+ABZ*I_ERI_Dxz_Py_D2y_S_bc;
  Double I_ERI_D2y_Dyz_D2y_S_bc = I_ERI_F2yz_Py_D2y_S_bc+ABZ*I_ERI_D2y_Py_D2y_S_bc;
  Double I_ERI_Dyz_Dyz_D2y_S_bc = I_ERI_Fy2z_Py_D2y_S_bc+ABZ*I_ERI_Dyz_Py_D2y_S_bc;
  Double I_ERI_D2z_Dyz_D2y_S_bc = I_ERI_F3z_Py_D2y_S_bc+ABZ*I_ERI_D2z_Py_D2y_S_bc;
  Double I_ERI_D2x_D2z_D2y_S_bc = I_ERI_F2xz_Pz_D2y_S_bc+ABZ*I_ERI_D2x_Pz_D2y_S_bc;
  Double I_ERI_Dxy_D2z_D2y_S_bc = I_ERI_Fxyz_Pz_D2y_S_bc+ABZ*I_ERI_Dxy_Pz_D2y_S_bc;
  Double I_ERI_Dxz_D2z_D2y_S_bc = I_ERI_Fx2z_Pz_D2y_S_bc+ABZ*I_ERI_Dxz_Pz_D2y_S_bc;
  Double I_ERI_D2y_D2z_D2y_S_bc = I_ERI_F2yz_Pz_D2y_S_bc+ABZ*I_ERI_D2y_Pz_D2y_S_bc;
  Double I_ERI_Dyz_D2z_D2y_S_bc = I_ERI_Fy2z_Pz_D2y_S_bc+ABZ*I_ERI_Dyz_Pz_D2y_S_bc;
  Double I_ERI_D2z_D2z_D2y_S_bc = I_ERI_F3z_Pz_D2y_S_bc+ABZ*I_ERI_D2z_Pz_D2y_S_bc;
  Double I_ERI_D2x_D2x_Dyz_S_bc = I_ERI_F3x_Px_Dyz_S_bc+ABX*I_ERI_D2x_Px_Dyz_S_bc;
  Double I_ERI_Dxy_D2x_Dyz_S_bc = I_ERI_F2xy_Px_Dyz_S_bc+ABX*I_ERI_Dxy_Px_Dyz_S_bc;
  Double I_ERI_Dxz_D2x_Dyz_S_bc = I_ERI_F2xz_Px_Dyz_S_bc+ABX*I_ERI_Dxz_Px_Dyz_S_bc;
  Double I_ERI_D2y_D2x_Dyz_S_bc = I_ERI_Fx2y_Px_Dyz_S_bc+ABX*I_ERI_D2y_Px_Dyz_S_bc;
  Double I_ERI_Dyz_D2x_Dyz_S_bc = I_ERI_Fxyz_Px_Dyz_S_bc+ABX*I_ERI_Dyz_Px_Dyz_S_bc;
  Double I_ERI_D2z_D2x_Dyz_S_bc = I_ERI_Fx2z_Px_Dyz_S_bc+ABX*I_ERI_D2z_Px_Dyz_S_bc;
  Double I_ERI_D2x_Dxy_Dyz_S_bc = I_ERI_F2xy_Px_Dyz_S_bc+ABY*I_ERI_D2x_Px_Dyz_S_bc;
  Double I_ERI_Dxy_Dxy_Dyz_S_bc = I_ERI_Fx2y_Px_Dyz_S_bc+ABY*I_ERI_Dxy_Px_Dyz_S_bc;
  Double I_ERI_Dxz_Dxy_Dyz_S_bc = I_ERI_Fxyz_Px_Dyz_S_bc+ABY*I_ERI_Dxz_Px_Dyz_S_bc;
  Double I_ERI_D2y_Dxy_Dyz_S_bc = I_ERI_F3y_Px_Dyz_S_bc+ABY*I_ERI_D2y_Px_Dyz_S_bc;
  Double I_ERI_Dyz_Dxy_Dyz_S_bc = I_ERI_F2yz_Px_Dyz_S_bc+ABY*I_ERI_Dyz_Px_Dyz_S_bc;
  Double I_ERI_D2z_Dxy_Dyz_S_bc = I_ERI_Fy2z_Px_Dyz_S_bc+ABY*I_ERI_D2z_Px_Dyz_S_bc;
  Double I_ERI_D2x_Dxz_Dyz_S_bc = I_ERI_F2xz_Px_Dyz_S_bc+ABZ*I_ERI_D2x_Px_Dyz_S_bc;
  Double I_ERI_Dxy_Dxz_Dyz_S_bc = I_ERI_Fxyz_Px_Dyz_S_bc+ABZ*I_ERI_Dxy_Px_Dyz_S_bc;
  Double I_ERI_Dxz_Dxz_Dyz_S_bc = I_ERI_Fx2z_Px_Dyz_S_bc+ABZ*I_ERI_Dxz_Px_Dyz_S_bc;
  Double I_ERI_D2y_Dxz_Dyz_S_bc = I_ERI_F2yz_Px_Dyz_S_bc+ABZ*I_ERI_D2y_Px_Dyz_S_bc;
  Double I_ERI_Dyz_Dxz_Dyz_S_bc = I_ERI_Fy2z_Px_Dyz_S_bc+ABZ*I_ERI_Dyz_Px_Dyz_S_bc;
  Double I_ERI_D2z_Dxz_Dyz_S_bc = I_ERI_F3z_Px_Dyz_S_bc+ABZ*I_ERI_D2z_Px_Dyz_S_bc;
  Double I_ERI_D2x_D2y_Dyz_S_bc = I_ERI_F2xy_Py_Dyz_S_bc+ABY*I_ERI_D2x_Py_Dyz_S_bc;
  Double I_ERI_Dxy_D2y_Dyz_S_bc = I_ERI_Fx2y_Py_Dyz_S_bc+ABY*I_ERI_Dxy_Py_Dyz_S_bc;
  Double I_ERI_Dxz_D2y_Dyz_S_bc = I_ERI_Fxyz_Py_Dyz_S_bc+ABY*I_ERI_Dxz_Py_Dyz_S_bc;
  Double I_ERI_D2y_D2y_Dyz_S_bc = I_ERI_F3y_Py_Dyz_S_bc+ABY*I_ERI_D2y_Py_Dyz_S_bc;
  Double I_ERI_Dyz_D2y_Dyz_S_bc = I_ERI_F2yz_Py_Dyz_S_bc+ABY*I_ERI_Dyz_Py_Dyz_S_bc;
  Double I_ERI_D2z_D2y_Dyz_S_bc = I_ERI_Fy2z_Py_Dyz_S_bc+ABY*I_ERI_D2z_Py_Dyz_S_bc;
  Double I_ERI_D2x_Dyz_Dyz_S_bc = I_ERI_F2xz_Py_Dyz_S_bc+ABZ*I_ERI_D2x_Py_Dyz_S_bc;
  Double I_ERI_Dxy_Dyz_Dyz_S_bc = I_ERI_Fxyz_Py_Dyz_S_bc+ABZ*I_ERI_Dxy_Py_Dyz_S_bc;
  Double I_ERI_Dxz_Dyz_Dyz_S_bc = I_ERI_Fx2z_Py_Dyz_S_bc+ABZ*I_ERI_Dxz_Py_Dyz_S_bc;
  Double I_ERI_D2y_Dyz_Dyz_S_bc = I_ERI_F2yz_Py_Dyz_S_bc+ABZ*I_ERI_D2y_Py_Dyz_S_bc;
  Double I_ERI_Dyz_Dyz_Dyz_S_bc = I_ERI_Fy2z_Py_Dyz_S_bc+ABZ*I_ERI_Dyz_Py_Dyz_S_bc;
  Double I_ERI_D2z_Dyz_Dyz_S_bc = I_ERI_F3z_Py_Dyz_S_bc+ABZ*I_ERI_D2z_Py_Dyz_S_bc;
  Double I_ERI_D2x_D2z_Dyz_S_bc = I_ERI_F2xz_Pz_Dyz_S_bc+ABZ*I_ERI_D2x_Pz_Dyz_S_bc;
  Double I_ERI_Dxy_D2z_Dyz_S_bc = I_ERI_Fxyz_Pz_Dyz_S_bc+ABZ*I_ERI_Dxy_Pz_Dyz_S_bc;
  Double I_ERI_Dxz_D2z_Dyz_S_bc = I_ERI_Fx2z_Pz_Dyz_S_bc+ABZ*I_ERI_Dxz_Pz_Dyz_S_bc;
  Double I_ERI_D2y_D2z_Dyz_S_bc = I_ERI_F2yz_Pz_Dyz_S_bc+ABZ*I_ERI_D2y_Pz_Dyz_S_bc;
  Double I_ERI_Dyz_D2z_Dyz_S_bc = I_ERI_Fy2z_Pz_Dyz_S_bc+ABZ*I_ERI_Dyz_Pz_Dyz_S_bc;
  Double I_ERI_D2z_D2z_Dyz_S_bc = I_ERI_F3z_Pz_Dyz_S_bc+ABZ*I_ERI_D2z_Pz_Dyz_S_bc;
  Double I_ERI_D2x_D2x_D2z_S_bc = I_ERI_F3x_Px_D2z_S_bc+ABX*I_ERI_D2x_Px_D2z_S_bc;
  Double I_ERI_Dxy_D2x_D2z_S_bc = I_ERI_F2xy_Px_D2z_S_bc+ABX*I_ERI_Dxy_Px_D2z_S_bc;
  Double I_ERI_Dxz_D2x_D2z_S_bc = I_ERI_F2xz_Px_D2z_S_bc+ABX*I_ERI_Dxz_Px_D2z_S_bc;
  Double I_ERI_D2y_D2x_D2z_S_bc = I_ERI_Fx2y_Px_D2z_S_bc+ABX*I_ERI_D2y_Px_D2z_S_bc;
  Double I_ERI_Dyz_D2x_D2z_S_bc = I_ERI_Fxyz_Px_D2z_S_bc+ABX*I_ERI_Dyz_Px_D2z_S_bc;
  Double I_ERI_D2z_D2x_D2z_S_bc = I_ERI_Fx2z_Px_D2z_S_bc+ABX*I_ERI_D2z_Px_D2z_S_bc;
  Double I_ERI_D2x_Dxy_D2z_S_bc = I_ERI_F2xy_Px_D2z_S_bc+ABY*I_ERI_D2x_Px_D2z_S_bc;
  Double I_ERI_Dxy_Dxy_D2z_S_bc = I_ERI_Fx2y_Px_D2z_S_bc+ABY*I_ERI_Dxy_Px_D2z_S_bc;
  Double I_ERI_Dxz_Dxy_D2z_S_bc = I_ERI_Fxyz_Px_D2z_S_bc+ABY*I_ERI_Dxz_Px_D2z_S_bc;
  Double I_ERI_D2y_Dxy_D2z_S_bc = I_ERI_F3y_Px_D2z_S_bc+ABY*I_ERI_D2y_Px_D2z_S_bc;
  Double I_ERI_Dyz_Dxy_D2z_S_bc = I_ERI_F2yz_Px_D2z_S_bc+ABY*I_ERI_Dyz_Px_D2z_S_bc;
  Double I_ERI_D2z_Dxy_D2z_S_bc = I_ERI_Fy2z_Px_D2z_S_bc+ABY*I_ERI_D2z_Px_D2z_S_bc;
  Double I_ERI_D2x_Dxz_D2z_S_bc = I_ERI_F2xz_Px_D2z_S_bc+ABZ*I_ERI_D2x_Px_D2z_S_bc;
  Double I_ERI_Dxy_Dxz_D2z_S_bc = I_ERI_Fxyz_Px_D2z_S_bc+ABZ*I_ERI_Dxy_Px_D2z_S_bc;
  Double I_ERI_Dxz_Dxz_D2z_S_bc = I_ERI_Fx2z_Px_D2z_S_bc+ABZ*I_ERI_Dxz_Px_D2z_S_bc;
  Double I_ERI_D2y_Dxz_D2z_S_bc = I_ERI_F2yz_Px_D2z_S_bc+ABZ*I_ERI_D2y_Px_D2z_S_bc;
  Double I_ERI_Dyz_Dxz_D2z_S_bc = I_ERI_Fy2z_Px_D2z_S_bc+ABZ*I_ERI_Dyz_Px_D2z_S_bc;
  Double I_ERI_D2z_Dxz_D2z_S_bc = I_ERI_F3z_Px_D2z_S_bc+ABZ*I_ERI_D2z_Px_D2z_S_bc;
  Double I_ERI_D2x_D2y_D2z_S_bc = I_ERI_F2xy_Py_D2z_S_bc+ABY*I_ERI_D2x_Py_D2z_S_bc;
  Double I_ERI_Dxy_D2y_D2z_S_bc = I_ERI_Fx2y_Py_D2z_S_bc+ABY*I_ERI_Dxy_Py_D2z_S_bc;
  Double I_ERI_Dxz_D2y_D2z_S_bc = I_ERI_Fxyz_Py_D2z_S_bc+ABY*I_ERI_Dxz_Py_D2z_S_bc;
  Double I_ERI_D2y_D2y_D2z_S_bc = I_ERI_F3y_Py_D2z_S_bc+ABY*I_ERI_D2y_Py_D2z_S_bc;
  Double I_ERI_Dyz_D2y_D2z_S_bc = I_ERI_F2yz_Py_D2z_S_bc+ABY*I_ERI_Dyz_Py_D2z_S_bc;
  Double I_ERI_D2z_D2y_D2z_S_bc = I_ERI_Fy2z_Py_D2z_S_bc+ABY*I_ERI_D2z_Py_D2z_S_bc;
  Double I_ERI_D2x_Dyz_D2z_S_bc = I_ERI_F2xz_Py_D2z_S_bc+ABZ*I_ERI_D2x_Py_D2z_S_bc;
  Double I_ERI_Dxy_Dyz_D2z_S_bc = I_ERI_Fxyz_Py_D2z_S_bc+ABZ*I_ERI_Dxy_Py_D2z_S_bc;
  Double I_ERI_Dxz_Dyz_D2z_S_bc = I_ERI_Fx2z_Py_D2z_S_bc+ABZ*I_ERI_Dxz_Py_D2z_S_bc;
  Double I_ERI_D2y_Dyz_D2z_S_bc = I_ERI_F2yz_Py_D2z_S_bc+ABZ*I_ERI_D2y_Py_D2z_S_bc;
  Double I_ERI_Dyz_Dyz_D2z_S_bc = I_ERI_Fy2z_Py_D2z_S_bc+ABZ*I_ERI_Dyz_Py_D2z_S_bc;
  Double I_ERI_D2z_Dyz_D2z_S_bc = I_ERI_F3z_Py_D2z_S_bc+ABZ*I_ERI_D2z_Py_D2z_S_bc;
  Double I_ERI_D2x_D2z_D2z_S_bc = I_ERI_F2xz_Pz_D2z_S_bc+ABZ*I_ERI_D2x_Pz_D2z_S_bc;
  Double I_ERI_Dxy_D2z_D2z_S_bc = I_ERI_Fxyz_Pz_D2z_S_bc+ABZ*I_ERI_Dxy_Pz_D2z_S_bc;
  Double I_ERI_Dxz_D2z_D2z_S_bc = I_ERI_Fx2z_Pz_D2z_S_bc+ABZ*I_ERI_Dxz_Pz_D2z_S_bc;
  Double I_ERI_D2y_D2z_D2z_S_bc = I_ERI_F2yz_Pz_D2z_S_bc+ABZ*I_ERI_D2y_Pz_D2z_S_bc;
  Double I_ERI_Dyz_D2z_D2z_S_bc = I_ERI_Fy2z_Pz_D2z_S_bc+ABZ*I_ERI_Dyz_Pz_D2z_S_bc;
  Double I_ERI_D2z_D2z_D2z_S_bc = I_ERI_F3z_Pz_D2z_S_bc+ABZ*I_ERI_D2z_Pz_D2z_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_F_S_cc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_cc
   * RHS shell quartet name: SQ_ERI_D_S_F_S_cc
   ************************************************************/
  Double I_ERI_D2x_Px_F3x_S_cc = I_ERI_F3x_S_F3x_S_cc+ABX*I_ERI_D2x_S_F3x_S_cc;
  Double I_ERI_Dxy_Px_F3x_S_cc = I_ERI_F2xy_S_F3x_S_cc+ABX*I_ERI_Dxy_S_F3x_S_cc;
  Double I_ERI_Dxz_Px_F3x_S_cc = I_ERI_F2xz_S_F3x_S_cc+ABX*I_ERI_Dxz_S_F3x_S_cc;
  Double I_ERI_D2y_Px_F3x_S_cc = I_ERI_Fx2y_S_F3x_S_cc+ABX*I_ERI_D2y_S_F3x_S_cc;
  Double I_ERI_Dyz_Px_F3x_S_cc = I_ERI_Fxyz_S_F3x_S_cc+ABX*I_ERI_Dyz_S_F3x_S_cc;
  Double I_ERI_D2z_Px_F3x_S_cc = I_ERI_Fx2z_S_F3x_S_cc+ABX*I_ERI_D2z_S_F3x_S_cc;
  Double I_ERI_D2x_Py_F3x_S_cc = I_ERI_F2xy_S_F3x_S_cc+ABY*I_ERI_D2x_S_F3x_S_cc;
  Double I_ERI_Dxy_Py_F3x_S_cc = I_ERI_Fx2y_S_F3x_S_cc+ABY*I_ERI_Dxy_S_F3x_S_cc;
  Double I_ERI_Dxz_Py_F3x_S_cc = I_ERI_Fxyz_S_F3x_S_cc+ABY*I_ERI_Dxz_S_F3x_S_cc;
  Double I_ERI_D2y_Py_F3x_S_cc = I_ERI_F3y_S_F3x_S_cc+ABY*I_ERI_D2y_S_F3x_S_cc;
  Double I_ERI_Dyz_Py_F3x_S_cc = I_ERI_F2yz_S_F3x_S_cc+ABY*I_ERI_Dyz_S_F3x_S_cc;
  Double I_ERI_D2z_Py_F3x_S_cc = I_ERI_Fy2z_S_F3x_S_cc+ABY*I_ERI_D2z_S_F3x_S_cc;
  Double I_ERI_D2x_Pz_F3x_S_cc = I_ERI_F2xz_S_F3x_S_cc+ABZ*I_ERI_D2x_S_F3x_S_cc;
  Double I_ERI_Dxy_Pz_F3x_S_cc = I_ERI_Fxyz_S_F3x_S_cc+ABZ*I_ERI_Dxy_S_F3x_S_cc;
  Double I_ERI_Dxz_Pz_F3x_S_cc = I_ERI_Fx2z_S_F3x_S_cc+ABZ*I_ERI_Dxz_S_F3x_S_cc;
  Double I_ERI_D2y_Pz_F3x_S_cc = I_ERI_F2yz_S_F3x_S_cc+ABZ*I_ERI_D2y_S_F3x_S_cc;
  Double I_ERI_Dyz_Pz_F3x_S_cc = I_ERI_Fy2z_S_F3x_S_cc+ABZ*I_ERI_Dyz_S_F3x_S_cc;
  Double I_ERI_D2z_Pz_F3x_S_cc = I_ERI_F3z_S_F3x_S_cc+ABZ*I_ERI_D2z_S_F3x_S_cc;
  Double I_ERI_D2x_Px_F2xy_S_cc = I_ERI_F3x_S_F2xy_S_cc+ABX*I_ERI_D2x_S_F2xy_S_cc;
  Double I_ERI_Dxy_Px_F2xy_S_cc = I_ERI_F2xy_S_F2xy_S_cc+ABX*I_ERI_Dxy_S_F2xy_S_cc;
  Double I_ERI_Dxz_Px_F2xy_S_cc = I_ERI_F2xz_S_F2xy_S_cc+ABX*I_ERI_Dxz_S_F2xy_S_cc;
  Double I_ERI_D2y_Px_F2xy_S_cc = I_ERI_Fx2y_S_F2xy_S_cc+ABX*I_ERI_D2y_S_F2xy_S_cc;
  Double I_ERI_Dyz_Px_F2xy_S_cc = I_ERI_Fxyz_S_F2xy_S_cc+ABX*I_ERI_Dyz_S_F2xy_S_cc;
  Double I_ERI_D2z_Px_F2xy_S_cc = I_ERI_Fx2z_S_F2xy_S_cc+ABX*I_ERI_D2z_S_F2xy_S_cc;
  Double I_ERI_D2x_Py_F2xy_S_cc = I_ERI_F2xy_S_F2xy_S_cc+ABY*I_ERI_D2x_S_F2xy_S_cc;
  Double I_ERI_Dxy_Py_F2xy_S_cc = I_ERI_Fx2y_S_F2xy_S_cc+ABY*I_ERI_Dxy_S_F2xy_S_cc;
  Double I_ERI_Dxz_Py_F2xy_S_cc = I_ERI_Fxyz_S_F2xy_S_cc+ABY*I_ERI_Dxz_S_F2xy_S_cc;
  Double I_ERI_D2y_Py_F2xy_S_cc = I_ERI_F3y_S_F2xy_S_cc+ABY*I_ERI_D2y_S_F2xy_S_cc;
  Double I_ERI_Dyz_Py_F2xy_S_cc = I_ERI_F2yz_S_F2xy_S_cc+ABY*I_ERI_Dyz_S_F2xy_S_cc;
  Double I_ERI_D2z_Py_F2xy_S_cc = I_ERI_Fy2z_S_F2xy_S_cc+ABY*I_ERI_D2z_S_F2xy_S_cc;
  Double I_ERI_D2x_Pz_F2xy_S_cc = I_ERI_F2xz_S_F2xy_S_cc+ABZ*I_ERI_D2x_S_F2xy_S_cc;
  Double I_ERI_Dxy_Pz_F2xy_S_cc = I_ERI_Fxyz_S_F2xy_S_cc+ABZ*I_ERI_Dxy_S_F2xy_S_cc;
  Double I_ERI_Dxz_Pz_F2xy_S_cc = I_ERI_Fx2z_S_F2xy_S_cc+ABZ*I_ERI_Dxz_S_F2xy_S_cc;
  Double I_ERI_D2y_Pz_F2xy_S_cc = I_ERI_F2yz_S_F2xy_S_cc+ABZ*I_ERI_D2y_S_F2xy_S_cc;
  Double I_ERI_Dyz_Pz_F2xy_S_cc = I_ERI_Fy2z_S_F2xy_S_cc+ABZ*I_ERI_Dyz_S_F2xy_S_cc;
  Double I_ERI_D2z_Pz_F2xy_S_cc = I_ERI_F3z_S_F2xy_S_cc+ABZ*I_ERI_D2z_S_F2xy_S_cc;
  Double I_ERI_D2x_Px_F2xz_S_cc = I_ERI_F3x_S_F2xz_S_cc+ABX*I_ERI_D2x_S_F2xz_S_cc;
  Double I_ERI_Dxy_Px_F2xz_S_cc = I_ERI_F2xy_S_F2xz_S_cc+ABX*I_ERI_Dxy_S_F2xz_S_cc;
  Double I_ERI_Dxz_Px_F2xz_S_cc = I_ERI_F2xz_S_F2xz_S_cc+ABX*I_ERI_Dxz_S_F2xz_S_cc;
  Double I_ERI_D2y_Px_F2xz_S_cc = I_ERI_Fx2y_S_F2xz_S_cc+ABX*I_ERI_D2y_S_F2xz_S_cc;
  Double I_ERI_Dyz_Px_F2xz_S_cc = I_ERI_Fxyz_S_F2xz_S_cc+ABX*I_ERI_Dyz_S_F2xz_S_cc;
  Double I_ERI_D2z_Px_F2xz_S_cc = I_ERI_Fx2z_S_F2xz_S_cc+ABX*I_ERI_D2z_S_F2xz_S_cc;
  Double I_ERI_D2x_Py_F2xz_S_cc = I_ERI_F2xy_S_F2xz_S_cc+ABY*I_ERI_D2x_S_F2xz_S_cc;
  Double I_ERI_Dxy_Py_F2xz_S_cc = I_ERI_Fx2y_S_F2xz_S_cc+ABY*I_ERI_Dxy_S_F2xz_S_cc;
  Double I_ERI_Dxz_Py_F2xz_S_cc = I_ERI_Fxyz_S_F2xz_S_cc+ABY*I_ERI_Dxz_S_F2xz_S_cc;
  Double I_ERI_D2y_Py_F2xz_S_cc = I_ERI_F3y_S_F2xz_S_cc+ABY*I_ERI_D2y_S_F2xz_S_cc;
  Double I_ERI_Dyz_Py_F2xz_S_cc = I_ERI_F2yz_S_F2xz_S_cc+ABY*I_ERI_Dyz_S_F2xz_S_cc;
  Double I_ERI_D2z_Py_F2xz_S_cc = I_ERI_Fy2z_S_F2xz_S_cc+ABY*I_ERI_D2z_S_F2xz_S_cc;
  Double I_ERI_D2x_Pz_F2xz_S_cc = I_ERI_F2xz_S_F2xz_S_cc+ABZ*I_ERI_D2x_S_F2xz_S_cc;
  Double I_ERI_Dxy_Pz_F2xz_S_cc = I_ERI_Fxyz_S_F2xz_S_cc+ABZ*I_ERI_Dxy_S_F2xz_S_cc;
  Double I_ERI_Dxz_Pz_F2xz_S_cc = I_ERI_Fx2z_S_F2xz_S_cc+ABZ*I_ERI_Dxz_S_F2xz_S_cc;
  Double I_ERI_D2y_Pz_F2xz_S_cc = I_ERI_F2yz_S_F2xz_S_cc+ABZ*I_ERI_D2y_S_F2xz_S_cc;
  Double I_ERI_Dyz_Pz_F2xz_S_cc = I_ERI_Fy2z_S_F2xz_S_cc+ABZ*I_ERI_Dyz_S_F2xz_S_cc;
  Double I_ERI_D2z_Pz_F2xz_S_cc = I_ERI_F3z_S_F2xz_S_cc+ABZ*I_ERI_D2z_S_F2xz_S_cc;
  Double I_ERI_D2x_Px_Fx2y_S_cc = I_ERI_F3x_S_Fx2y_S_cc+ABX*I_ERI_D2x_S_Fx2y_S_cc;
  Double I_ERI_Dxy_Px_Fx2y_S_cc = I_ERI_F2xy_S_Fx2y_S_cc+ABX*I_ERI_Dxy_S_Fx2y_S_cc;
  Double I_ERI_Dxz_Px_Fx2y_S_cc = I_ERI_F2xz_S_Fx2y_S_cc+ABX*I_ERI_Dxz_S_Fx2y_S_cc;
  Double I_ERI_D2y_Px_Fx2y_S_cc = I_ERI_Fx2y_S_Fx2y_S_cc+ABX*I_ERI_D2y_S_Fx2y_S_cc;
  Double I_ERI_Dyz_Px_Fx2y_S_cc = I_ERI_Fxyz_S_Fx2y_S_cc+ABX*I_ERI_Dyz_S_Fx2y_S_cc;
  Double I_ERI_D2z_Px_Fx2y_S_cc = I_ERI_Fx2z_S_Fx2y_S_cc+ABX*I_ERI_D2z_S_Fx2y_S_cc;
  Double I_ERI_D2x_Py_Fx2y_S_cc = I_ERI_F2xy_S_Fx2y_S_cc+ABY*I_ERI_D2x_S_Fx2y_S_cc;
  Double I_ERI_Dxy_Py_Fx2y_S_cc = I_ERI_Fx2y_S_Fx2y_S_cc+ABY*I_ERI_Dxy_S_Fx2y_S_cc;
  Double I_ERI_Dxz_Py_Fx2y_S_cc = I_ERI_Fxyz_S_Fx2y_S_cc+ABY*I_ERI_Dxz_S_Fx2y_S_cc;
  Double I_ERI_D2y_Py_Fx2y_S_cc = I_ERI_F3y_S_Fx2y_S_cc+ABY*I_ERI_D2y_S_Fx2y_S_cc;
  Double I_ERI_Dyz_Py_Fx2y_S_cc = I_ERI_F2yz_S_Fx2y_S_cc+ABY*I_ERI_Dyz_S_Fx2y_S_cc;
  Double I_ERI_D2z_Py_Fx2y_S_cc = I_ERI_Fy2z_S_Fx2y_S_cc+ABY*I_ERI_D2z_S_Fx2y_S_cc;
  Double I_ERI_D2x_Pz_Fx2y_S_cc = I_ERI_F2xz_S_Fx2y_S_cc+ABZ*I_ERI_D2x_S_Fx2y_S_cc;
  Double I_ERI_Dxy_Pz_Fx2y_S_cc = I_ERI_Fxyz_S_Fx2y_S_cc+ABZ*I_ERI_Dxy_S_Fx2y_S_cc;
  Double I_ERI_Dxz_Pz_Fx2y_S_cc = I_ERI_Fx2z_S_Fx2y_S_cc+ABZ*I_ERI_Dxz_S_Fx2y_S_cc;
  Double I_ERI_D2y_Pz_Fx2y_S_cc = I_ERI_F2yz_S_Fx2y_S_cc+ABZ*I_ERI_D2y_S_Fx2y_S_cc;
  Double I_ERI_Dyz_Pz_Fx2y_S_cc = I_ERI_Fy2z_S_Fx2y_S_cc+ABZ*I_ERI_Dyz_S_Fx2y_S_cc;
  Double I_ERI_D2z_Pz_Fx2y_S_cc = I_ERI_F3z_S_Fx2y_S_cc+ABZ*I_ERI_D2z_S_Fx2y_S_cc;
  Double I_ERI_D2x_Px_Fxyz_S_cc = I_ERI_F3x_S_Fxyz_S_cc+ABX*I_ERI_D2x_S_Fxyz_S_cc;
  Double I_ERI_Dxy_Px_Fxyz_S_cc = I_ERI_F2xy_S_Fxyz_S_cc+ABX*I_ERI_Dxy_S_Fxyz_S_cc;
  Double I_ERI_Dxz_Px_Fxyz_S_cc = I_ERI_F2xz_S_Fxyz_S_cc+ABX*I_ERI_Dxz_S_Fxyz_S_cc;
  Double I_ERI_D2y_Px_Fxyz_S_cc = I_ERI_Fx2y_S_Fxyz_S_cc+ABX*I_ERI_D2y_S_Fxyz_S_cc;
  Double I_ERI_Dyz_Px_Fxyz_S_cc = I_ERI_Fxyz_S_Fxyz_S_cc+ABX*I_ERI_Dyz_S_Fxyz_S_cc;
  Double I_ERI_D2z_Px_Fxyz_S_cc = I_ERI_Fx2z_S_Fxyz_S_cc+ABX*I_ERI_D2z_S_Fxyz_S_cc;
  Double I_ERI_D2x_Py_Fxyz_S_cc = I_ERI_F2xy_S_Fxyz_S_cc+ABY*I_ERI_D2x_S_Fxyz_S_cc;
  Double I_ERI_Dxy_Py_Fxyz_S_cc = I_ERI_Fx2y_S_Fxyz_S_cc+ABY*I_ERI_Dxy_S_Fxyz_S_cc;
  Double I_ERI_Dxz_Py_Fxyz_S_cc = I_ERI_Fxyz_S_Fxyz_S_cc+ABY*I_ERI_Dxz_S_Fxyz_S_cc;
  Double I_ERI_D2y_Py_Fxyz_S_cc = I_ERI_F3y_S_Fxyz_S_cc+ABY*I_ERI_D2y_S_Fxyz_S_cc;
  Double I_ERI_Dyz_Py_Fxyz_S_cc = I_ERI_F2yz_S_Fxyz_S_cc+ABY*I_ERI_Dyz_S_Fxyz_S_cc;
  Double I_ERI_D2z_Py_Fxyz_S_cc = I_ERI_Fy2z_S_Fxyz_S_cc+ABY*I_ERI_D2z_S_Fxyz_S_cc;
  Double I_ERI_D2x_Pz_Fxyz_S_cc = I_ERI_F2xz_S_Fxyz_S_cc+ABZ*I_ERI_D2x_S_Fxyz_S_cc;
  Double I_ERI_Dxy_Pz_Fxyz_S_cc = I_ERI_Fxyz_S_Fxyz_S_cc+ABZ*I_ERI_Dxy_S_Fxyz_S_cc;
  Double I_ERI_Dxz_Pz_Fxyz_S_cc = I_ERI_Fx2z_S_Fxyz_S_cc+ABZ*I_ERI_Dxz_S_Fxyz_S_cc;
  Double I_ERI_D2y_Pz_Fxyz_S_cc = I_ERI_F2yz_S_Fxyz_S_cc+ABZ*I_ERI_D2y_S_Fxyz_S_cc;
  Double I_ERI_Dyz_Pz_Fxyz_S_cc = I_ERI_Fy2z_S_Fxyz_S_cc+ABZ*I_ERI_Dyz_S_Fxyz_S_cc;
  Double I_ERI_D2z_Pz_Fxyz_S_cc = I_ERI_F3z_S_Fxyz_S_cc+ABZ*I_ERI_D2z_S_Fxyz_S_cc;
  Double I_ERI_D2x_Px_Fx2z_S_cc = I_ERI_F3x_S_Fx2z_S_cc+ABX*I_ERI_D2x_S_Fx2z_S_cc;
  Double I_ERI_Dxy_Px_Fx2z_S_cc = I_ERI_F2xy_S_Fx2z_S_cc+ABX*I_ERI_Dxy_S_Fx2z_S_cc;
  Double I_ERI_Dxz_Px_Fx2z_S_cc = I_ERI_F2xz_S_Fx2z_S_cc+ABX*I_ERI_Dxz_S_Fx2z_S_cc;
  Double I_ERI_D2y_Px_Fx2z_S_cc = I_ERI_Fx2y_S_Fx2z_S_cc+ABX*I_ERI_D2y_S_Fx2z_S_cc;
  Double I_ERI_Dyz_Px_Fx2z_S_cc = I_ERI_Fxyz_S_Fx2z_S_cc+ABX*I_ERI_Dyz_S_Fx2z_S_cc;
  Double I_ERI_D2z_Px_Fx2z_S_cc = I_ERI_Fx2z_S_Fx2z_S_cc+ABX*I_ERI_D2z_S_Fx2z_S_cc;
  Double I_ERI_D2x_Py_Fx2z_S_cc = I_ERI_F2xy_S_Fx2z_S_cc+ABY*I_ERI_D2x_S_Fx2z_S_cc;
  Double I_ERI_Dxy_Py_Fx2z_S_cc = I_ERI_Fx2y_S_Fx2z_S_cc+ABY*I_ERI_Dxy_S_Fx2z_S_cc;
  Double I_ERI_Dxz_Py_Fx2z_S_cc = I_ERI_Fxyz_S_Fx2z_S_cc+ABY*I_ERI_Dxz_S_Fx2z_S_cc;
  Double I_ERI_D2y_Py_Fx2z_S_cc = I_ERI_F3y_S_Fx2z_S_cc+ABY*I_ERI_D2y_S_Fx2z_S_cc;
  Double I_ERI_Dyz_Py_Fx2z_S_cc = I_ERI_F2yz_S_Fx2z_S_cc+ABY*I_ERI_Dyz_S_Fx2z_S_cc;
  Double I_ERI_D2z_Py_Fx2z_S_cc = I_ERI_Fy2z_S_Fx2z_S_cc+ABY*I_ERI_D2z_S_Fx2z_S_cc;
  Double I_ERI_D2x_Pz_Fx2z_S_cc = I_ERI_F2xz_S_Fx2z_S_cc+ABZ*I_ERI_D2x_S_Fx2z_S_cc;
  Double I_ERI_Dxy_Pz_Fx2z_S_cc = I_ERI_Fxyz_S_Fx2z_S_cc+ABZ*I_ERI_Dxy_S_Fx2z_S_cc;
  Double I_ERI_Dxz_Pz_Fx2z_S_cc = I_ERI_Fx2z_S_Fx2z_S_cc+ABZ*I_ERI_Dxz_S_Fx2z_S_cc;
  Double I_ERI_D2y_Pz_Fx2z_S_cc = I_ERI_F2yz_S_Fx2z_S_cc+ABZ*I_ERI_D2y_S_Fx2z_S_cc;
  Double I_ERI_Dyz_Pz_Fx2z_S_cc = I_ERI_Fy2z_S_Fx2z_S_cc+ABZ*I_ERI_Dyz_S_Fx2z_S_cc;
  Double I_ERI_D2z_Pz_Fx2z_S_cc = I_ERI_F3z_S_Fx2z_S_cc+ABZ*I_ERI_D2z_S_Fx2z_S_cc;
  Double I_ERI_D2x_Px_F3y_S_cc = I_ERI_F3x_S_F3y_S_cc+ABX*I_ERI_D2x_S_F3y_S_cc;
  Double I_ERI_Dxy_Px_F3y_S_cc = I_ERI_F2xy_S_F3y_S_cc+ABX*I_ERI_Dxy_S_F3y_S_cc;
  Double I_ERI_Dxz_Px_F3y_S_cc = I_ERI_F2xz_S_F3y_S_cc+ABX*I_ERI_Dxz_S_F3y_S_cc;
  Double I_ERI_D2y_Px_F3y_S_cc = I_ERI_Fx2y_S_F3y_S_cc+ABX*I_ERI_D2y_S_F3y_S_cc;
  Double I_ERI_Dyz_Px_F3y_S_cc = I_ERI_Fxyz_S_F3y_S_cc+ABX*I_ERI_Dyz_S_F3y_S_cc;
  Double I_ERI_D2z_Px_F3y_S_cc = I_ERI_Fx2z_S_F3y_S_cc+ABX*I_ERI_D2z_S_F3y_S_cc;
  Double I_ERI_D2x_Py_F3y_S_cc = I_ERI_F2xy_S_F3y_S_cc+ABY*I_ERI_D2x_S_F3y_S_cc;
  Double I_ERI_Dxy_Py_F3y_S_cc = I_ERI_Fx2y_S_F3y_S_cc+ABY*I_ERI_Dxy_S_F3y_S_cc;
  Double I_ERI_Dxz_Py_F3y_S_cc = I_ERI_Fxyz_S_F3y_S_cc+ABY*I_ERI_Dxz_S_F3y_S_cc;
  Double I_ERI_D2y_Py_F3y_S_cc = I_ERI_F3y_S_F3y_S_cc+ABY*I_ERI_D2y_S_F3y_S_cc;
  Double I_ERI_Dyz_Py_F3y_S_cc = I_ERI_F2yz_S_F3y_S_cc+ABY*I_ERI_Dyz_S_F3y_S_cc;
  Double I_ERI_D2z_Py_F3y_S_cc = I_ERI_Fy2z_S_F3y_S_cc+ABY*I_ERI_D2z_S_F3y_S_cc;
  Double I_ERI_D2x_Pz_F3y_S_cc = I_ERI_F2xz_S_F3y_S_cc+ABZ*I_ERI_D2x_S_F3y_S_cc;
  Double I_ERI_Dxy_Pz_F3y_S_cc = I_ERI_Fxyz_S_F3y_S_cc+ABZ*I_ERI_Dxy_S_F3y_S_cc;
  Double I_ERI_Dxz_Pz_F3y_S_cc = I_ERI_Fx2z_S_F3y_S_cc+ABZ*I_ERI_Dxz_S_F3y_S_cc;
  Double I_ERI_D2y_Pz_F3y_S_cc = I_ERI_F2yz_S_F3y_S_cc+ABZ*I_ERI_D2y_S_F3y_S_cc;
  Double I_ERI_Dyz_Pz_F3y_S_cc = I_ERI_Fy2z_S_F3y_S_cc+ABZ*I_ERI_Dyz_S_F3y_S_cc;
  Double I_ERI_D2z_Pz_F3y_S_cc = I_ERI_F3z_S_F3y_S_cc+ABZ*I_ERI_D2z_S_F3y_S_cc;
  Double I_ERI_D2x_Px_F2yz_S_cc = I_ERI_F3x_S_F2yz_S_cc+ABX*I_ERI_D2x_S_F2yz_S_cc;
  Double I_ERI_Dxy_Px_F2yz_S_cc = I_ERI_F2xy_S_F2yz_S_cc+ABX*I_ERI_Dxy_S_F2yz_S_cc;
  Double I_ERI_Dxz_Px_F2yz_S_cc = I_ERI_F2xz_S_F2yz_S_cc+ABX*I_ERI_Dxz_S_F2yz_S_cc;
  Double I_ERI_D2y_Px_F2yz_S_cc = I_ERI_Fx2y_S_F2yz_S_cc+ABX*I_ERI_D2y_S_F2yz_S_cc;
  Double I_ERI_Dyz_Px_F2yz_S_cc = I_ERI_Fxyz_S_F2yz_S_cc+ABX*I_ERI_Dyz_S_F2yz_S_cc;
  Double I_ERI_D2z_Px_F2yz_S_cc = I_ERI_Fx2z_S_F2yz_S_cc+ABX*I_ERI_D2z_S_F2yz_S_cc;
  Double I_ERI_D2x_Py_F2yz_S_cc = I_ERI_F2xy_S_F2yz_S_cc+ABY*I_ERI_D2x_S_F2yz_S_cc;
  Double I_ERI_Dxy_Py_F2yz_S_cc = I_ERI_Fx2y_S_F2yz_S_cc+ABY*I_ERI_Dxy_S_F2yz_S_cc;
  Double I_ERI_Dxz_Py_F2yz_S_cc = I_ERI_Fxyz_S_F2yz_S_cc+ABY*I_ERI_Dxz_S_F2yz_S_cc;
  Double I_ERI_D2y_Py_F2yz_S_cc = I_ERI_F3y_S_F2yz_S_cc+ABY*I_ERI_D2y_S_F2yz_S_cc;
  Double I_ERI_Dyz_Py_F2yz_S_cc = I_ERI_F2yz_S_F2yz_S_cc+ABY*I_ERI_Dyz_S_F2yz_S_cc;
  Double I_ERI_D2z_Py_F2yz_S_cc = I_ERI_Fy2z_S_F2yz_S_cc+ABY*I_ERI_D2z_S_F2yz_S_cc;
  Double I_ERI_D2x_Pz_F2yz_S_cc = I_ERI_F2xz_S_F2yz_S_cc+ABZ*I_ERI_D2x_S_F2yz_S_cc;
  Double I_ERI_Dxy_Pz_F2yz_S_cc = I_ERI_Fxyz_S_F2yz_S_cc+ABZ*I_ERI_Dxy_S_F2yz_S_cc;
  Double I_ERI_Dxz_Pz_F2yz_S_cc = I_ERI_Fx2z_S_F2yz_S_cc+ABZ*I_ERI_Dxz_S_F2yz_S_cc;
  Double I_ERI_D2y_Pz_F2yz_S_cc = I_ERI_F2yz_S_F2yz_S_cc+ABZ*I_ERI_D2y_S_F2yz_S_cc;
  Double I_ERI_Dyz_Pz_F2yz_S_cc = I_ERI_Fy2z_S_F2yz_S_cc+ABZ*I_ERI_Dyz_S_F2yz_S_cc;
  Double I_ERI_D2z_Pz_F2yz_S_cc = I_ERI_F3z_S_F2yz_S_cc+ABZ*I_ERI_D2z_S_F2yz_S_cc;
  Double I_ERI_D2x_Px_Fy2z_S_cc = I_ERI_F3x_S_Fy2z_S_cc+ABX*I_ERI_D2x_S_Fy2z_S_cc;
  Double I_ERI_Dxy_Px_Fy2z_S_cc = I_ERI_F2xy_S_Fy2z_S_cc+ABX*I_ERI_Dxy_S_Fy2z_S_cc;
  Double I_ERI_Dxz_Px_Fy2z_S_cc = I_ERI_F2xz_S_Fy2z_S_cc+ABX*I_ERI_Dxz_S_Fy2z_S_cc;
  Double I_ERI_D2y_Px_Fy2z_S_cc = I_ERI_Fx2y_S_Fy2z_S_cc+ABX*I_ERI_D2y_S_Fy2z_S_cc;
  Double I_ERI_Dyz_Px_Fy2z_S_cc = I_ERI_Fxyz_S_Fy2z_S_cc+ABX*I_ERI_Dyz_S_Fy2z_S_cc;
  Double I_ERI_D2z_Px_Fy2z_S_cc = I_ERI_Fx2z_S_Fy2z_S_cc+ABX*I_ERI_D2z_S_Fy2z_S_cc;
  Double I_ERI_D2x_Py_Fy2z_S_cc = I_ERI_F2xy_S_Fy2z_S_cc+ABY*I_ERI_D2x_S_Fy2z_S_cc;
  Double I_ERI_Dxy_Py_Fy2z_S_cc = I_ERI_Fx2y_S_Fy2z_S_cc+ABY*I_ERI_Dxy_S_Fy2z_S_cc;
  Double I_ERI_Dxz_Py_Fy2z_S_cc = I_ERI_Fxyz_S_Fy2z_S_cc+ABY*I_ERI_Dxz_S_Fy2z_S_cc;
  Double I_ERI_D2y_Py_Fy2z_S_cc = I_ERI_F3y_S_Fy2z_S_cc+ABY*I_ERI_D2y_S_Fy2z_S_cc;
  Double I_ERI_Dyz_Py_Fy2z_S_cc = I_ERI_F2yz_S_Fy2z_S_cc+ABY*I_ERI_Dyz_S_Fy2z_S_cc;
  Double I_ERI_D2z_Py_Fy2z_S_cc = I_ERI_Fy2z_S_Fy2z_S_cc+ABY*I_ERI_D2z_S_Fy2z_S_cc;
  Double I_ERI_D2x_Pz_Fy2z_S_cc = I_ERI_F2xz_S_Fy2z_S_cc+ABZ*I_ERI_D2x_S_Fy2z_S_cc;
  Double I_ERI_Dxy_Pz_Fy2z_S_cc = I_ERI_Fxyz_S_Fy2z_S_cc+ABZ*I_ERI_Dxy_S_Fy2z_S_cc;
  Double I_ERI_Dxz_Pz_Fy2z_S_cc = I_ERI_Fx2z_S_Fy2z_S_cc+ABZ*I_ERI_Dxz_S_Fy2z_S_cc;
  Double I_ERI_D2y_Pz_Fy2z_S_cc = I_ERI_F2yz_S_Fy2z_S_cc+ABZ*I_ERI_D2y_S_Fy2z_S_cc;
  Double I_ERI_Dyz_Pz_Fy2z_S_cc = I_ERI_Fy2z_S_Fy2z_S_cc+ABZ*I_ERI_Dyz_S_Fy2z_S_cc;
  Double I_ERI_D2z_Pz_Fy2z_S_cc = I_ERI_F3z_S_Fy2z_S_cc+ABZ*I_ERI_D2z_S_Fy2z_S_cc;
  Double I_ERI_D2x_Px_F3z_S_cc = I_ERI_F3x_S_F3z_S_cc+ABX*I_ERI_D2x_S_F3z_S_cc;
  Double I_ERI_Dxy_Px_F3z_S_cc = I_ERI_F2xy_S_F3z_S_cc+ABX*I_ERI_Dxy_S_F3z_S_cc;
  Double I_ERI_Dxz_Px_F3z_S_cc = I_ERI_F2xz_S_F3z_S_cc+ABX*I_ERI_Dxz_S_F3z_S_cc;
  Double I_ERI_D2y_Px_F3z_S_cc = I_ERI_Fx2y_S_F3z_S_cc+ABX*I_ERI_D2y_S_F3z_S_cc;
  Double I_ERI_Dyz_Px_F3z_S_cc = I_ERI_Fxyz_S_F3z_S_cc+ABX*I_ERI_Dyz_S_F3z_S_cc;
  Double I_ERI_D2z_Px_F3z_S_cc = I_ERI_Fx2z_S_F3z_S_cc+ABX*I_ERI_D2z_S_F3z_S_cc;
  Double I_ERI_D2x_Py_F3z_S_cc = I_ERI_F2xy_S_F3z_S_cc+ABY*I_ERI_D2x_S_F3z_S_cc;
  Double I_ERI_Dxy_Py_F3z_S_cc = I_ERI_Fx2y_S_F3z_S_cc+ABY*I_ERI_Dxy_S_F3z_S_cc;
  Double I_ERI_Dxz_Py_F3z_S_cc = I_ERI_Fxyz_S_F3z_S_cc+ABY*I_ERI_Dxz_S_F3z_S_cc;
  Double I_ERI_D2y_Py_F3z_S_cc = I_ERI_F3y_S_F3z_S_cc+ABY*I_ERI_D2y_S_F3z_S_cc;
  Double I_ERI_Dyz_Py_F3z_S_cc = I_ERI_F2yz_S_F3z_S_cc+ABY*I_ERI_Dyz_S_F3z_S_cc;
  Double I_ERI_D2z_Py_F3z_S_cc = I_ERI_Fy2z_S_F3z_S_cc+ABY*I_ERI_D2z_S_F3z_S_cc;
  Double I_ERI_D2x_Pz_F3z_S_cc = I_ERI_F2xz_S_F3z_S_cc+ABZ*I_ERI_D2x_S_F3z_S_cc;
  Double I_ERI_Dxy_Pz_F3z_S_cc = I_ERI_Fxyz_S_F3z_S_cc+ABZ*I_ERI_Dxy_S_F3z_S_cc;
  Double I_ERI_Dxz_Pz_F3z_S_cc = I_ERI_Fx2z_S_F3z_S_cc+ABZ*I_ERI_Dxz_S_F3z_S_cc;
  Double I_ERI_D2y_Pz_F3z_S_cc = I_ERI_F2yz_S_F3z_S_cc+ABZ*I_ERI_D2y_S_F3z_S_cc;
  Double I_ERI_Dyz_Pz_F3z_S_cc = I_ERI_Fy2z_S_F3z_S_cc+ABZ*I_ERI_Dyz_S_F3z_S_cc;
  Double I_ERI_D2z_Pz_F3z_S_cc = I_ERI_F3z_S_F3z_S_cc+ABZ*I_ERI_D2z_S_F3z_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_aa
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S
   ************************************************************/
  abcd[0] = 4.0E0*I_ERI_G4x_Px_Px_S_aa-2.0E0*2*I_ERI_D2x_Px_Px_S_a-2.0E0*3*I_ERI_D2x_Px_Px_S_a+2*1*I_ERI_S_Px_Px_S;
  abcd[1] = 4.0E0*I_ERI_G3xy_Px_Px_S_aa-2.0E0*1*I_ERI_Dxy_Px_Px_S_a-2.0E0*2*I_ERI_Dxy_Px_Px_S_a;
  abcd[2] = 4.0E0*I_ERI_G3xz_Px_Px_S_aa-2.0E0*1*I_ERI_Dxz_Px_Px_S_a-2.0E0*2*I_ERI_Dxz_Px_Px_S_a;
  abcd[3] = 4.0E0*I_ERI_G2x2y_Px_Px_S_aa-2.0E0*1*I_ERI_D2y_Px_Px_S_a;
  abcd[4] = 4.0E0*I_ERI_G2xyz_Px_Px_S_aa-2.0E0*1*I_ERI_Dyz_Px_Px_S_a;
  abcd[5] = 4.0E0*I_ERI_G2x2z_Px_Px_S_aa-2.0E0*1*I_ERI_D2z_Px_Px_S_a;
  abcd[6] = 4.0E0*I_ERI_G4x_Py_Px_S_aa-2.0E0*2*I_ERI_D2x_Py_Px_S_a-2.0E0*3*I_ERI_D2x_Py_Px_S_a+2*1*I_ERI_S_Py_Px_S;
  abcd[7] = 4.0E0*I_ERI_G3xy_Py_Px_S_aa-2.0E0*1*I_ERI_Dxy_Py_Px_S_a-2.0E0*2*I_ERI_Dxy_Py_Px_S_a;
  abcd[8] = 4.0E0*I_ERI_G3xz_Py_Px_S_aa-2.0E0*1*I_ERI_Dxz_Py_Px_S_a-2.0E0*2*I_ERI_Dxz_Py_Px_S_a;
  abcd[9] = 4.0E0*I_ERI_G2x2y_Py_Px_S_aa-2.0E0*1*I_ERI_D2y_Py_Px_S_a;
  abcd[10] = 4.0E0*I_ERI_G2xyz_Py_Px_S_aa-2.0E0*1*I_ERI_Dyz_Py_Px_S_a;
  abcd[11] = 4.0E0*I_ERI_G2x2z_Py_Px_S_aa-2.0E0*1*I_ERI_D2z_Py_Px_S_a;
  abcd[12] = 4.0E0*I_ERI_G4x_Pz_Px_S_aa-2.0E0*2*I_ERI_D2x_Pz_Px_S_a-2.0E0*3*I_ERI_D2x_Pz_Px_S_a+2*1*I_ERI_S_Pz_Px_S;
  abcd[13] = 4.0E0*I_ERI_G3xy_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Px_S_a-2.0E0*2*I_ERI_Dxy_Pz_Px_S_a;
  abcd[14] = 4.0E0*I_ERI_G3xz_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Px_S_a-2.0E0*2*I_ERI_Dxz_Pz_Px_S_a;
  abcd[15] = 4.0E0*I_ERI_G2x2y_Pz_Px_S_aa-2.0E0*1*I_ERI_D2y_Pz_Px_S_a;
  abcd[16] = 4.0E0*I_ERI_G2xyz_Pz_Px_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Px_S_a;
  abcd[17] = 4.0E0*I_ERI_G2x2z_Pz_Px_S_aa-2.0E0*1*I_ERI_D2z_Pz_Px_S_a;
  abcd[18] = 4.0E0*I_ERI_G4x_Px_Py_S_aa-2.0E0*2*I_ERI_D2x_Px_Py_S_a-2.0E0*3*I_ERI_D2x_Px_Py_S_a+2*1*I_ERI_S_Px_Py_S;
  abcd[19] = 4.0E0*I_ERI_G3xy_Px_Py_S_aa-2.0E0*1*I_ERI_Dxy_Px_Py_S_a-2.0E0*2*I_ERI_Dxy_Px_Py_S_a;
  abcd[20] = 4.0E0*I_ERI_G3xz_Px_Py_S_aa-2.0E0*1*I_ERI_Dxz_Px_Py_S_a-2.0E0*2*I_ERI_Dxz_Px_Py_S_a;
  abcd[21] = 4.0E0*I_ERI_G2x2y_Px_Py_S_aa-2.0E0*1*I_ERI_D2y_Px_Py_S_a;
  abcd[22] = 4.0E0*I_ERI_G2xyz_Px_Py_S_aa-2.0E0*1*I_ERI_Dyz_Px_Py_S_a;
  abcd[23] = 4.0E0*I_ERI_G2x2z_Px_Py_S_aa-2.0E0*1*I_ERI_D2z_Px_Py_S_a;
  abcd[24] = 4.0E0*I_ERI_G4x_Py_Py_S_aa-2.0E0*2*I_ERI_D2x_Py_Py_S_a-2.0E0*3*I_ERI_D2x_Py_Py_S_a+2*1*I_ERI_S_Py_Py_S;
  abcd[25] = 4.0E0*I_ERI_G3xy_Py_Py_S_aa-2.0E0*1*I_ERI_Dxy_Py_Py_S_a-2.0E0*2*I_ERI_Dxy_Py_Py_S_a;
  abcd[26] = 4.0E0*I_ERI_G3xz_Py_Py_S_aa-2.0E0*1*I_ERI_Dxz_Py_Py_S_a-2.0E0*2*I_ERI_Dxz_Py_Py_S_a;
  abcd[27] = 4.0E0*I_ERI_G2x2y_Py_Py_S_aa-2.0E0*1*I_ERI_D2y_Py_Py_S_a;
  abcd[28] = 4.0E0*I_ERI_G2xyz_Py_Py_S_aa-2.0E0*1*I_ERI_Dyz_Py_Py_S_a;
  abcd[29] = 4.0E0*I_ERI_G2x2z_Py_Py_S_aa-2.0E0*1*I_ERI_D2z_Py_Py_S_a;
  abcd[30] = 4.0E0*I_ERI_G4x_Pz_Py_S_aa-2.0E0*2*I_ERI_D2x_Pz_Py_S_a-2.0E0*3*I_ERI_D2x_Pz_Py_S_a+2*1*I_ERI_S_Pz_Py_S;
  abcd[31] = 4.0E0*I_ERI_G3xy_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Py_S_a-2.0E0*2*I_ERI_Dxy_Pz_Py_S_a;
  abcd[32] = 4.0E0*I_ERI_G3xz_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Py_S_a-2.0E0*2*I_ERI_Dxz_Pz_Py_S_a;
  abcd[33] = 4.0E0*I_ERI_G2x2y_Pz_Py_S_aa-2.0E0*1*I_ERI_D2y_Pz_Py_S_a;
  abcd[34] = 4.0E0*I_ERI_G2xyz_Pz_Py_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Py_S_a;
  abcd[35] = 4.0E0*I_ERI_G2x2z_Pz_Py_S_aa-2.0E0*1*I_ERI_D2z_Pz_Py_S_a;
  abcd[36] = 4.0E0*I_ERI_G4x_Px_Pz_S_aa-2.0E0*2*I_ERI_D2x_Px_Pz_S_a-2.0E0*3*I_ERI_D2x_Px_Pz_S_a+2*1*I_ERI_S_Px_Pz_S;
  abcd[37] = 4.0E0*I_ERI_G3xy_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Px_Pz_S_a-2.0E0*2*I_ERI_Dxy_Px_Pz_S_a;
  abcd[38] = 4.0E0*I_ERI_G3xz_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Px_Pz_S_a-2.0E0*2*I_ERI_Dxz_Px_Pz_S_a;
  abcd[39] = 4.0E0*I_ERI_G2x2y_Px_Pz_S_aa-2.0E0*1*I_ERI_D2y_Px_Pz_S_a;
  abcd[40] = 4.0E0*I_ERI_G2xyz_Px_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Px_Pz_S_a;
  abcd[41] = 4.0E0*I_ERI_G2x2z_Px_Pz_S_aa-2.0E0*1*I_ERI_D2z_Px_Pz_S_a;
  abcd[42] = 4.0E0*I_ERI_G4x_Py_Pz_S_aa-2.0E0*2*I_ERI_D2x_Py_Pz_S_a-2.0E0*3*I_ERI_D2x_Py_Pz_S_a+2*1*I_ERI_S_Py_Pz_S;
  abcd[43] = 4.0E0*I_ERI_G3xy_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Py_Pz_S_a-2.0E0*2*I_ERI_Dxy_Py_Pz_S_a;
  abcd[44] = 4.0E0*I_ERI_G3xz_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Py_Pz_S_a-2.0E0*2*I_ERI_Dxz_Py_Pz_S_a;
  abcd[45] = 4.0E0*I_ERI_G2x2y_Py_Pz_S_aa-2.0E0*1*I_ERI_D2y_Py_Pz_S_a;
  abcd[46] = 4.0E0*I_ERI_G2xyz_Py_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Py_Pz_S_a;
  abcd[47] = 4.0E0*I_ERI_G2x2z_Py_Pz_S_aa-2.0E0*1*I_ERI_D2z_Py_Pz_S_a;
  abcd[48] = 4.0E0*I_ERI_G4x_Pz_Pz_S_aa-2.0E0*2*I_ERI_D2x_Pz_Pz_S_a-2.0E0*3*I_ERI_D2x_Pz_Pz_S_a+2*1*I_ERI_S_Pz_Pz_S;
  abcd[49] = 4.0E0*I_ERI_G3xy_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_a-2.0E0*2*I_ERI_Dxy_Pz_Pz_S_a;
  abcd[50] = 4.0E0*I_ERI_G3xz_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_a-2.0E0*2*I_ERI_Dxz_Pz_Pz_S_a;
  abcd[51] = 4.0E0*I_ERI_G2x2y_Pz_Pz_S_aa-2.0E0*1*I_ERI_D2y_Pz_Pz_S_a;
  abcd[52] = 4.0E0*I_ERI_G2xyz_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_a;
  abcd[53] = 4.0E0*I_ERI_G2x2z_Pz_Pz_S_aa-2.0E0*1*I_ERI_D2z_Pz_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_aa
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S
   ************************************************************/
  abcd[54] = 4.0E0*I_ERI_G3xy_Px_Px_S_aa-2.0E0*2*I_ERI_Dxy_Px_Px_S_a;
  abcd[55] = 4.0E0*I_ERI_G2x2y_Px_Px_S_aa-2.0E0*1*I_ERI_D2x_Px_Px_S_a-2.0E0*1*I_ERI_D2y_Px_Px_S_a+1*I_ERI_S_Px_Px_S;
  abcd[56] = 4.0E0*I_ERI_G2xyz_Px_Px_S_aa-2.0E0*1*I_ERI_Dyz_Px_Px_S_a;
  abcd[57] = 4.0E0*I_ERI_Gx3y_Px_Px_S_aa-2.0E0*2*I_ERI_Dxy_Px_Px_S_a;
  abcd[58] = 4.0E0*I_ERI_Gx2yz_Px_Px_S_aa-2.0E0*1*I_ERI_Dxz_Px_Px_S_a;
  abcd[59] = 4.0E0*I_ERI_Gxy2z_Px_Px_S_aa;
  abcd[60] = 4.0E0*I_ERI_G3xy_Py_Px_S_aa-2.0E0*2*I_ERI_Dxy_Py_Px_S_a;
  abcd[61] = 4.0E0*I_ERI_G2x2y_Py_Px_S_aa-2.0E0*1*I_ERI_D2x_Py_Px_S_a-2.0E0*1*I_ERI_D2y_Py_Px_S_a+1*I_ERI_S_Py_Px_S;
  abcd[62] = 4.0E0*I_ERI_G2xyz_Py_Px_S_aa-2.0E0*1*I_ERI_Dyz_Py_Px_S_a;
  abcd[63] = 4.0E0*I_ERI_Gx3y_Py_Px_S_aa-2.0E0*2*I_ERI_Dxy_Py_Px_S_a;
  abcd[64] = 4.0E0*I_ERI_Gx2yz_Py_Px_S_aa-2.0E0*1*I_ERI_Dxz_Py_Px_S_a;
  abcd[65] = 4.0E0*I_ERI_Gxy2z_Py_Px_S_aa;
  abcd[66] = 4.0E0*I_ERI_G3xy_Pz_Px_S_aa-2.0E0*2*I_ERI_Dxy_Pz_Px_S_a;
  abcd[67] = 4.0E0*I_ERI_G2x2y_Pz_Px_S_aa-2.0E0*1*I_ERI_D2x_Pz_Px_S_a-2.0E0*1*I_ERI_D2y_Pz_Px_S_a+1*I_ERI_S_Pz_Px_S;
  abcd[68] = 4.0E0*I_ERI_G2xyz_Pz_Px_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Px_S_a;
  abcd[69] = 4.0E0*I_ERI_Gx3y_Pz_Px_S_aa-2.0E0*2*I_ERI_Dxy_Pz_Px_S_a;
  abcd[70] = 4.0E0*I_ERI_Gx2yz_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Px_S_a;
  abcd[71] = 4.0E0*I_ERI_Gxy2z_Pz_Px_S_aa;
  abcd[72] = 4.0E0*I_ERI_G3xy_Px_Py_S_aa-2.0E0*2*I_ERI_Dxy_Px_Py_S_a;
  abcd[73] = 4.0E0*I_ERI_G2x2y_Px_Py_S_aa-2.0E0*1*I_ERI_D2x_Px_Py_S_a-2.0E0*1*I_ERI_D2y_Px_Py_S_a+1*I_ERI_S_Px_Py_S;
  abcd[74] = 4.0E0*I_ERI_G2xyz_Px_Py_S_aa-2.0E0*1*I_ERI_Dyz_Px_Py_S_a;
  abcd[75] = 4.0E0*I_ERI_Gx3y_Px_Py_S_aa-2.0E0*2*I_ERI_Dxy_Px_Py_S_a;
  abcd[76] = 4.0E0*I_ERI_Gx2yz_Px_Py_S_aa-2.0E0*1*I_ERI_Dxz_Px_Py_S_a;
  abcd[77] = 4.0E0*I_ERI_Gxy2z_Px_Py_S_aa;
  abcd[78] = 4.0E0*I_ERI_G3xy_Py_Py_S_aa-2.0E0*2*I_ERI_Dxy_Py_Py_S_a;
  abcd[79] = 4.0E0*I_ERI_G2x2y_Py_Py_S_aa-2.0E0*1*I_ERI_D2x_Py_Py_S_a-2.0E0*1*I_ERI_D2y_Py_Py_S_a+1*I_ERI_S_Py_Py_S;
  abcd[80] = 4.0E0*I_ERI_G2xyz_Py_Py_S_aa-2.0E0*1*I_ERI_Dyz_Py_Py_S_a;
  abcd[81] = 4.0E0*I_ERI_Gx3y_Py_Py_S_aa-2.0E0*2*I_ERI_Dxy_Py_Py_S_a;
  abcd[82] = 4.0E0*I_ERI_Gx2yz_Py_Py_S_aa-2.0E0*1*I_ERI_Dxz_Py_Py_S_a;
  abcd[83] = 4.0E0*I_ERI_Gxy2z_Py_Py_S_aa;
  abcd[84] = 4.0E0*I_ERI_G3xy_Pz_Py_S_aa-2.0E0*2*I_ERI_Dxy_Pz_Py_S_a;
  abcd[85] = 4.0E0*I_ERI_G2x2y_Pz_Py_S_aa-2.0E0*1*I_ERI_D2x_Pz_Py_S_a-2.0E0*1*I_ERI_D2y_Pz_Py_S_a+1*I_ERI_S_Pz_Py_S;
  abcd[86] = 4.0E0*I_ERI_G2xyz_Pz_Py_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Py_S_a;
  abcd[87] = 4.0E0*I_ERI_Gx3y_Pz_Py_S_aa-2.0E0*2*I_ERI_Dxy_Pz_Py_S_a;
  abcd[88] = 4.0E0*I_ERI_Gx2yz_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Py_S_a;
  abcd[89] = 4.0E0*I_ERI_Gxy2z_Pz_Py_S_aa;
  abcd[90] = 4.0E0*I_ERI_G3xy_Px_Pz_S_aa-2.0E0*2*I_ERI_Dxy_Px_Pz_S_a;
  abcd[91] = 4.0E0*I_ERI_G2x2y_Px_Pz_S_aa-2.0E0*1*I_ERI_D2x_Px_Pz_S_a-2.0E0*1*I_ERI_D2y_Px_Pz_S_a+1*I_ERI_S_Px_Pz_S;
  abcd[92] = 4.0E0*I_ERI_G2xyz_Px_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Px_Pz_S_a;
  abcd[93] = 4.0E0*I_ERI_Gx3y_Px_Pz_S_aa-2.0E0*2*I_ERI_Dxy_Px_Pz_S_a;
  abcd[94] = 4.0E0*I_ERI_Gx2yz_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Px_Pz_S_a;
  abcd[95] = 4.0E0*I_ERI_Gxy2z_Px_Pz_S_aa;
  abcd[96] = 4.0E0*I_ERI_G3xy_Py_Pz_S_aa-2.0E0*2*I_ERI_Dxy_Py_Pz_S_a;
  abcd[97] = 4.0E0*I_ERI_G2x2y_Py_Pz_S_aa-2.0E0*1*I_ERI_D2x_Py_Pz_S_a-2.0E0*1*I_ERI_D2y_Py_Pz_S_a+1*I_ERI_S_Py_Pz_S;
  abcd[98] = 4.0E0*I_ERI_G2xyz_Py_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Py_Pz_S_a;
  abcd[99] = 4.0E0*I_ERI_Gx3y_Py_Pz_S_aa-2.0E0*2*I_ERI_Dxy_Py_Pz_S_a;
  abcd[100] = 4.0E0*I_ERI_Gx2yz_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Py_Pz_S_a;
  abcd[101] = 4.0E0*I_ERI_Gxy2z_Py_Pz_S_aa;
  abcd[102] = 4.0E0*I_ERI_G3xy_Pz_Pz_S_aa-2.0E0*2*I_ERI_Dxy_Pz_Pz_S_a;
  abcd[103] = 4.0E0*I_ERI_G2x2y_Pz_Pz_S_aa-2.0E0*1*I_ERI_D2x_Pz_Pz_S_a-2.0E0*1*I_ERI_D2y_Pz_Pz_S_a+1*I_ERI_S_Pz_Pz_S;
  abcd[104] = 4.0E0*I_ERI_G2xyz_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_a;
  abcd[105] = 4.0E0*I_ERI_Gx3y_Pz_Pz_S_aa-2.0E0*2*I_ERI_Dxy_Pz_Pz_S_a;
  abcd[106] = 4.0E0*I_ERI_Gx2yz_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_a;
  abcd[107] = 4.0E0*I_ERI_Gxy2z_Pz_Pz_S_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_aa
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S
   ************************************************************/
  abcd[108] = 4.0E0*I_ERI_G3xz_Px_Px_S_aa-2.0E0*2*I_ERI_Dxz_Px_Px_S_a;
  abcd[109] = 4.0E0*I_ERI_G2xyz_Px_Px_S_aa-2.0E0*1*I_ERI_Dyz_Px_Px_S_a;
  abcd[110] = 4.0E0*I_ERI_G2x2z_Px_Px_S_aa-2.0E0*1*I_ERI_D2x_Px_Px_S_a-2.0E0*1*I_ERI_D2z_Px_Px_S_a+1*I_ERI_S_Px_Px_S;
  abcd[111] = 4.0E0*I_ERI_Gx2yz_Px_Px_S_aa;
  abcd[112] = 4.0E0*I_ERI_Gxy2z_Px_Px_S_aa-2.0E0*1*I_ERI_Dxy_Px_Px_S_a;
  abcd[113] = 4.0E0*I_ERI_Gx3z_Px_Px_S_aa-2.0E0*2*I_ERI_Dxz_Px_Px_S_a;
  abcd[114] = 4.0E0*I_ERI_G3xz_Py_Px_S_aa-2.0E0*2*I_ERI_Dxz_Py_Px_S_a;
  abcd[115] = 4.0E0*I_ERI_G2xyz_Py_Px_S_aa-2.0E0*1*I_ERI_Dyz_Py_Px_S_a;
  abcd[116] = 4.0E0*I_ERI_G2x2z_Py_Px_S_aa-2.0E0*1*I_ERI_D2x_Py_Px_S_a-2.0E0*1*I_ERI_D2z_Py_Px_S_a+1*I_ERI_S_Py_Px_S;
  abcd[117] = 4.0E0*I_ERI_Gx2yz_Py_Px_S_aa;
  abcd[118] = 4.0E0*I_ERI_Gxy2z_Py_Px_S_aa-2.0E0*1*I_ERI_Dxy_Py_Px_S_a;
  abcd[119] = 4.0E0*I_ERI_Gx3z_Py_Px_S_aa-2.0E0*2*I_ERI_Dxz_Py_Px_S_a;
  abcd[120] = 4.0E0*I_ERI_G3xz_Pz_Px_S_aa-2.0E0*2*I_ERI_Dxz_Pz_Px_S_a;
  abcd[121] = 4.0E0*I_ERI_G2xyz_Pz_Px_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Px_S_a;
  abcd[122] = 4.0E0*I_ERI_G2x2z_Pz_Px_S_aa-2.0E0*1*I_ERI_D2x_Pz_Px_S_a-2.0E0*1*I_ERI_D2z_Pz_Px_S_a+1*I_ERI_S_Pz_Px_S;
  abcd[123] = 4.0E0*I_ERI_Gx2yz_Pz_Px_S_aa;
  abcd[124] = 4.0E0*I_ERI_Gxy2z_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Px_S_a;
  abcd[125] = 4.0E0*I_ERI_Gx3z_Pz_Px_S_aa-2.0E0*2*I_ERI_Dxz_Pz_Px_S_a;
  abcd[126] = 4.0E0*I_ERI_G3xz_Px_Py_S_aa-2.0E0*2*I_ERI_Dxz_Px_Py_S_a;
  abcd[127] = 4.0E0*I_ERI_G2xyz_Px_Py_S_aa-2.0E0*1*I_ERI_Dyz_Px_Py_S_a;
  abcd[128] = 4.0E0*I_ERI_G2x2z_Px_Py_S_aa-2.0E0*1*I_ERI_D2x_Px_Py_S_a-2.0E0*1*I_ERI_D2z_Px_Py_S_a+1*I_ERI_S_Px_Py_S;
  abcd[129] = 4.0E0*I_ERI_Gx2yz_Px_Py_S_aa;
  abcd[130] = 4.0E0*I_ERI_Gxy2z_Px_Py_S_aa-2.0E0*1*I_ERI_Dxy_Px_Py_S_a;
  abcd[131] = 4.0E0*I_ERI_Gx3z_Px_Py_S_aa-2.0E0*2*I_ERI_Dxz_Px_Py_S_a;
  abcd[132] = 4.0E0*I_ERI_G3xz_Py_Py_S_aa-2.0E0*2*I_ERI_Dxz_Py_Py_S_a;
  abcd[133] = 4.0E0*I_ERI_G2xyz_Py_Py_S_aa-2.0E0*1*I_ERI_Dyz_Py_Py_S_a;
  abcd[134] = 4.0E0*I_ERI_G2x2z_Py_Py_S_aa-2.0E0*1*I_ERI_D2x_Py_Py_S_a-2.0E0*1*I_ERI_D2z_Py_Py_S_a+1*I_ERI_S_Py_Py_S;
  abcd[135] = 4.0E0*I_ERI_Gx2yz_Py_Py_S_aa;
  abcd[136] = 4.0E0*I_ERI_Gxy2z_Py_Py_S_aa-2.0E0*1*I_ERI_Dxy_Py_Py_S_a;
  abcd[137] = 4.0E0*I_ERI_Gx3z_Py_Py_S_aa-2.0E0*2*I_ERI_Dxz_Py_Py_S_a;
  abcd[138] = 4.0E0*I_ERI_G3xz_Pz_Py_S_aa-2.0E0*2*I_ERI_Dxz_Pz_Py_S_a;
  abcd[139] = 4.0E0*I_ERI_G2xyz_Pz_Py_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Py_S_a;
  abcd[140] = 4.0E0*I_ERI_G2x2z_Pz_Py_S_aa-2.0E0*1*I_ERI_D2x_Pz_Py_S_a-2.0E0*1*I_ERI_D2z_Pz_Py_S_a+1*I_ERI_S_Pz_Py_S;
  abcd[141] = 4.0E0*I_ERI_Gx2yz_Pz_Py_S_aa;
  abcd[142] = 4.0E0*I_ERI_Gxy2z_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Py_S_a;
  abcd[143] = 4.0E0*I_ERI_Gx3z_Pz_Py_S_aa-2.0E0*2*I_ERI_Dxz_Pz_Py_S_a;
  abcd[144] = 4.0E0*I_ERI_G3xz_Px_Pz_S_aa-2.0E0*2*I_ERI_Dxz_Px_Pz_S_a;
  abcd[145] = 4.0E0*I_ERI_G2xyz_Px_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Px_Pz_S_a;
  abcd[146] = 4.0E0*I_ERI_G2x2z_Px_Pz_S_aa-2.0E0*1*I_ERI_D2x_Px_Pz_S_a-2.0E0*1*I_ERI_D2z_Px_Pz_S_a+1*I_ERI_S_Px_Pz_S;
  abcd[147] = 4.0E0*I_ERI_Gx2yz_Px_Pz_S_aa;
  abcd[148] = 4.0E0*I_ERI_Gxy2z_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Px_Pz_S_a;
  abcd[149] = 4.0E0*I_ERI_Gx3z_Px_Pz_S_aa-2.0E0*2*I_ERI_Dxz_Px_Pz_S_a;
  abcd[150] = 4.0E0*I_ERI_G3xz_Py_Pz_S_aa-2.0E0*2*I_ERI_Dxz_Py_Pz_S_a;
  abcd[151] = 4.0E0*I_ERI_G2xyz_Py_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Py_Pz_S_a;
  abcd[152] = 4.0E0*I_ERI_G2x2z_Py_Pz_S_aa-2.0E0*1*I_ERI_D2x_Py_Pz_S_a-2.0E0*1*I_ERI_D2z_Py_Pz_S_a+1*I_ERI_S_Py_Pz_S;
  abcd[153] = 4.0E0*I_ERI_Gx2yz_Py_Pz_S_aa;
  abcd[154] = 4.0E0*I_ERI_Gxy2z_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Py_Pz_S_a;
  abcd[155] = 4.0E0*I_ERI_Gx3z_Py_Pz_S_aa-2.0E0*2*I_ERI_Dxz_Py_Pz_S_a;
  abcd[156] = 4.0E0*I_ERI_G3xz_Pz_Pz_S_aa-2.0E0*2*I_ERI_Dxz_Pz_Pz_S_a;
  abcd[157] = 4.0E0*I_ERI_G2xyz_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_a;
  abcd[158] = 4.0E0*I_ERI_G2x2z_Pz_Pz_S_aa-2.0E0*1*I_ERI_D2x_Pz_Pz_S_a-2.0E0*1*I_ERI_D2z_Pz_Pz_S_a+1*I_ERI_S_Pz_Pz_S;
  abcd[159] = 4.0E0*I_ERI_Gx2yz_Pz_Pz_S_aa;
  abcd[160] = 4.0E0*I_ERI_Gxy2z_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_a;
  abcd[161] = 4.0E0*I_ERI_Gx3z_Pz_Pz_S_aa-2.0E0*2*I_ERI_Dxz_Pz_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_aa
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S
   ************************************************************/
  abcd[162] = 4.0E0*I_ERI_G2x2y_Px_Px_S_aa-2.0E0*1*I_ERI_D2x_Px_Px_S_a;
  abcd[163] = 4.0E0*I_ERI_Gx3y_Px_Px_S_aa-2.0E0*1*I_ERI_Dxy_Px_Px_S_a-2.0E0*2*I_ERI_Dxy_Px_Px_S_a;
  abcd[164] = 4.0E0*I_ERI_Gx2yz_Px_Px_S_aa-2.0E0*1*I_ERI_Dxz_Px_Px_S_a;
  abcd[165] = 4.0E0*I_ERI_G4y_Px_Px_S_aa-2.0E0*2*I_ERI_D2y_Px_Px_S_a-2.0E0*3*I_ERI_D2y_Px_Px_S_a+2*1*I_ERI_S_Px_Px_S;
  abcd[166] = 4.0E0*I_ERI_G3yz_Px_Px_S_aa-2.0E0*1*I_ERI_Dyz_Px_Px_S_a-2.0E0*2*I_ERI_Dyz_Px_Px_S_a;
  abcd[167] = 4.0E0*I_ERI_G2y2z_Px_Px_S_aa-2.0E0*1*I_ERI_D2z_Px_Px_S_a;
  abcd[168] = 4.0E0*I_ERI_G2x2y_Py_Px_S_aa-2.0E0*1*I_ERI_D2x_Py_Px_S_a;
  abcd[169] = 4.0E0*I_ERI_Gx3y_Py_Px_S_aa-2.0E0*1*I_ERI_Dxy_Py_Px_S_a-2.0E0*2*I_ERI_Dxy_Py_Px_S_a;
  abcd[170] = 4.0E0*I_ERI_Gx2yz_Py_Px_S_aa-2.0E0*1*I_ERI_Dxz_Py_Px_S_a;
  abcd[171] = 4.0E0*I_ERI_G4y_Py_Px_S_aa-2.0E0*2*I_ERI_D2y_Py_Px_S_a-2.0E0*3*I_ERI_D2y_Py_Px_S_a+2*1*I_ERI_S_Py_Px_S;
  abcd[172] = 4.0E0*I_ERI_G3yz_Py_Px_S_aa-2.0E0*1*I_ERI_Dyz_Py_Px_S_a-2.0E0*2*I_ERI_Dyz_Py_Px_S_a;
  abcd[173] = 4.0E0*I_ERI_G2y2z_Py_Px_S_aa-2.0E0*1*I_ERI_D2z_Py_Px_S_a;
  abcd[174] = 4.0E0*I_ERI_G2x2y_Pz_Px_S_aa-2.0E0*1*I_ERI_D2x_Pz_Px_S_a;
  abcd[175] = 4.0E0*I_ERI_Gx3y_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Px_S_a-2.0E0*2*I_ERI_Dxy_Pz_Px_S_a;
  abcd[176] = 4.0E0*I_ERI_Gx2yz_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Px_S_a;
  abcd[177] = 4.0E0*I_ERI_G4y_Pz_Px_S_aa-2.0E0*2*I_ERI_D2y_Pz_Px_S_a-2.0E0*3*I_ERI_D2y_Pz_Px_S_a+2*1*I_ERI_S_Pz_Px_S;
  abcd[178] = 4.0E0*I_ERI_G3yz_Pz_Px_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Px_S_a-2.0E0*2*I_ERI_Dyz_Pz_Px_S_a;
  abcd[179] = 4.0E0*I_ERI_G2y2z_Pz_Px_S_aa-2.0E0*1*I_ERI_D2z_Pz_Px_S_a;
  abcd[180] = 4.0E0*I_ERI_G2x2y_Px_Py_S_aa-2.0E0*1*I_ERI_D2x_Px_Py_S_a;
  abcd[181] = 4.0E0*I_ERI_Gx3y_Px_Py_S_aa-2.0E0*1*I_ERI_Dxy_Px_Py_S_a-2.0E0*2*I_ERI_Dxy_Px_Py_S_a;
  abcd[182] = 4.0E0*I_ERI_Gx2yz_Px_Py_S_aa-2.0E0*1*I_ERI_Dxz_Px_Py_S_a;
  abcd[183] = 4.0E0*I_ERI_G4y_Px_Py_S_aa-2.0E0*2*I_ERI_D2y_Px_Py_S_a-2.0E0*3*I_ERI_D2y_Px_Py_S_a+2*1*I_ERI_S_Px_Py_S;
  abcd[184] = 4.0E0*I_ERI_G3yz_Px_Py_S_aa-2.0E0*1*I_ERI_Dyz_Px_Py_S_a-2.0E0*2*I_ERI_Dyz_Px_Py_S_a;
  abcd[185] = 4.0E0*I_ERI_G2y2z_Px_Py_S_aa-2.0E0*1*I_ERI_D2z_Px_Py_S_a;
  abcd[186] = 4.0E0*I_ERI_G2x2y_Py_Py_S_aa-2.0E0*1*I_ERI_D2x_Py_Py_S_a;
  abcd[187] = 4.0E0*I_ERI_Gx3y_Py_Py_S_aa-2.0E0*1*I_ERI_Dxy_Py_Py_S_a-2.0E0*2*I_ERI_Dxy_Py_Py_S_a;
  abcd[188] = 4.0E0*I_ERI_Gx2yz_Py_Py_S_aa-2.0E0*1*I_ERI_Dxz_Py_Py_S_a;
  abcd[189] = 4.0E0*I_ERI_G4y_Py_Py_S_aa-2.0E0*2*I_ERI_D2y_Py_Py_S_a-2.0E0*3*I_ERI_D2y_Py_Py_S_a+2*1*I_ERI_S_Py_Py_S;
  abcd[190] = 4.0E0*I_ERI_G3yz_Py_Py_S_aa-2.0E0*1*I_ERI_Dyz_Py_Py_S_a-2.0E0*2*I_ERI_Dyz_Py_Py_S_a;
  abcd[191] = 4.0E0*I_ERI_G2y2z_Py_Py_S_aa-2.0E0*1*I_ERI_D2z_Py_Py_S_a;
  abcd[192] = 4.0E0*I_ERI_G2x2y_Pz_Py_S_aa-2.0E0*1*I_ERI_D2x_Pz_Py_S_a;
  abcd[193] = 4.0E0*I_ERI_Gx3y_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Py_S_a-2.0E0*2*I_ERI_Dxy_Pz_Py_S_a;
  abcd[194] = 4.0E0*I_ERI_Gx2yz_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Py_S_a;
  abcd[195] = 4.0E0*I_ERI_G4y_Pz_Py_S_aa-2.0E0*2*I_ERI_D2y_Pz_Py_S_a-2.0E0*3*I_ERI_D2y_Pz_Py_S_a+2*1*I_ERI_S_Pz_Py_S;
  abcd[196] = 4.0E0*I_ERI_G3yz_Pz_Py_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Py_S_a-2.0E0*2*I_ERI_Dyz_Pz_Py_S_a;
  abcd[197] = 4.0E0*I_ERI_G2y2z_Pz_Py_S_aa-2.0E0*1*I_ERI_D2z_Pz_Py_S_a;
  abcd[198] = 4.0E0*I_ERI_G2x2y_Px_Pz_S_aa-2.0E0*1*I_ERI_D2x_Px_Pz_S_a;
  abcd[199] = 4.0E0*I_ERI_Gx3y_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Px_Pz_S_a-2.0E0*2*I_ERI_Dxy_Px_Pz_S_a;
  abcd[200] = 4.0E0*I_ERI_Gx2yz_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Px_Pz_S_a;
  abcd[201] = 4.0E0*I_ERI_G4y_Px_Pz_S_aa-2.0E0*2*I_ERI_D2y_Px_Pz_S_a-2.0E0*3*I_ERI_D2y_Px_Pz_S_a+2*1*I_ERI_S_Px_Pz_S;
  abcd[202] = 4.0E0*I_ERI_G3yz_Px_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Px_Pz_S_a-2.0E0*2*I_ERI_Dyz_Px_Pz_S_a;
  abcd[203] = 4.0E0*I_ERI_G2y2z_Px_Pz_S_aa-2.0E0*1*I_ERI_D2z_Px_Pz_S_a;
  abcd[204] = 4.0E0*I_ERI_G2x2y_Py_Pz_S_aa-2.0E0*1*I_ERI_D2x_Py_Pz_S_a;
  abcd[205] = 4.0E0*I_ERI_Gx3y_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Py_Pz_S_a-2.0E0*2*I_ERI_Dxy_Py_Pz_S_a;
  abcd[206] = 4.0E0*I_ERI_Gx2yz_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Py_Pz_S_a;
  abcd[207] = 4.0E0*I_ERI_G4y_Py_Pz_S_aa-2.0E0*2*I_ERI_D2y_Py_Pz_S_a-2.0E0*3*I_ERI_D2y_Py_Pz_S_a+2*1*I_ERI_S_Py_Pz_S;
  abcd[208] = 4.0E0*I_ERI_G3yz_Py_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Py_Pz_S_a-2.0E0*2*I_ERI_Dyz_Py_Pz_S_a;
  abcd[209] = 4.0E0*I_ERI_G2y2z_Py_Pz_S_aa-2.0E0*1*I_ERI_D2z_Py_Pz_S_a;
  abcd[210] = 4.0E0*I_ERI_G2x2y_Pz_Pz_S_aa-2.0E0*1*I_ERI_D2x_Pz_Pz_S_a;
  abcd[211] = 4.0E0*I_ERI_Gx3y_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_a-2.0E0*2*I_ERI_Dxy_Pz_Pz_S_a;
  abcd[212] = 4.0E0*I_ERI_Gx2yz_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_a;
  abcd[213] = 4.0E0*I_ERI_G4y_Pz_Pz_S_aa-2.0E0*2*I_ERI_D2y_Pz_Pz_S_a-2.0E0*3*I_ERI_D2y_Pz_Pz_S_a+2*1*I_ERI_S_Pz_Pz_S;
  abcd[214] = 4.0E0*I_ERI_G3yz_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_a-2.0E0*2*I_ERI_Dyz_Pz_Pz_S_a;
  abcd[215] = 4.0E0*I_ERI_G2y2z_Pz_Pz_S_aa-2.0E0*1*I_ERI_D2z_Pz_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_aa
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S
   ************************************************************/
  abcd[216] = 4.0E0*I_ERI_G2xyz_Px_Px_S_aa;
  abcd[217] = 4.0E0*I_ERI_Gx2yz_Px_Px_S_aa-2.0E0*1*I_ERI_Dxz_Px_Px_S_a;
  abcd[218] = 4.0E0*I_ERI_Gxy2z_Px_Px_S_aa-2.0E0*1*I_ERI_Dxy_Px_Px_S_a;
  abcd[219] = 4.0E0*I_ERI_G3yz_Px_Px_S_aa-2.0E0*2*I_ERI_Dyz_Px_Px_S_a;
  abcd[220] = 4.0E0*I_ERI_G2y2z_Px_Px_S_aa-2.0E0*1*I_ERI_D2y_Px_Px_S_a-2.0E0*1*I_ERI_D2z_Px_Px_S_a+1*I_ERI_S_Px_Px_S;
  abcd[221] = 4.0E0*I_ERI_Gy3z_Px_Px_S_aa-2.0E0*2*I_ERI_Dyz_Px_Px_S_a;
  abcd[222] = 4.0E0*I_ERI_G2xyz_Py_Px_S_aa;
  abcd[223] = 4.0E0*I_ERI_Gx2yz_Py_Px_S_aa-2.0E0*1*I_ERI_Dxz_Py_Px_S_a;
  abcd[224] = 4.0E0*I_ERI_Gxy2z_Py_Px_S_aa-2.0E0*1*I_ERI_Dxy_Py_Px_S_a;
  abcd[225] = 4.0E0*I_ERI_G3yz_Py_Px_S_aa-2.0E0*2*I_ERI_Dyz_Py_Px_S_a;
  abcd[226] = 4.0E0*I_ERI_G2y2z_Py_Px_S_aa-2.0E0*1*I_ERI_D2y_Py_Px_S_a-2.0E0*1*I_ERI_D2z_Py_Px_S_a+1*I_ERI_S_Py_Px_S;
  abcd[227] = 4.0E0*I_ERI_Gy3z_Py_Px_S_aa-2.0E0*2*I_ERI_Dyz_Py_Px_S_a;
  abcd[228] = 4.0E0*I_ERI_G2xyz_Pz_Px_S_aa;
  abcd[229] = 4.0E0*I_ERI_Gx2yz_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Px_S_a;
  abcd[230] = 4.0E0*I_ERI_Gxy2z_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Px_S_a;
  abcd[231] = 4.0E0*I_ERI_G3yz_Pz_Px_S_aa-2.0E0*2*I_ERI_Dyz_Pz_Px_S_a;
  abcd[232] = 4.0E0*I_ERI_G2y2z_Pz_Px_S_aa-2.0E0*1*I_ERI_D2y_Pz_Px_S_a-2.0E0*1*I_ERI_D2z_Pz_Px_S_a+1*I_ERI_S_Pz_Px_S;
  abcd[233] = 4.0E0*I_ERI_Gy3z_Pz_Px_S_aa-2.0E0*2*I_ERI_Dyz_Pz_Px_S_a;
  abcd[234] = 4.0E0*I_ERI_G2xyz_Px_Py_S_aa;
  abcd[235] = 4.0E0*I_ERI_Gx2yz_Px_Py_S_aa-2.0E0*1*I_ERI_Dxz_Px_Py_S_a;
  abcd[236] = 4.0E0*I_ERI_Gxy2z_Px_Py_S_aa-2.0E0*1*I_ERI_Dxy_Px_Py_S_a;
  abcd[237] = 4.0E0*I_ERI_G3yz_Px_Py_S_aa-2.0E0*2*I_ERI_Dyz_Px_Py_S_a;
  abcd[238] = 4.0E0*I_ERI_G2y2z_Px_Py_S_aa-2.0E0*1*I_ERI_D2y_Px_Py_S_a-2.0E0*1*I_ERI_D2z_Px_Py_S_a+1*I_ERI_S_Px_Py_S;
  abcd[239] = 4.0E0*I_ERI_Gy3z_Px_Py_S_aa-2.0E0*2*I_ERI_Dyz_Px_Py_S_a;
  abcd[240] = 4.0E0*I_ERI_G2xyz_Py_Py_S_aa;
  abcd[241] = 4.0E0*I_ERI_Gx2yz_Py_Py_S_aa-2.0E0*1*I_ERI_Dxz_Py_Py_S_a;
  abcd[242] = 4.0E0*I_ERI_Gxy2z_Py_Py_S_aa-2.0E0*1*I_ERI_Dxy_Py_Py_S_a;
  abcd[243] = 4.0E0*I_ERI_G3yz_Py_Py_S_aa-2.0E0*2*I_ERI_Dyz_Py_Py_S_a;
  abcd[244] = 4.0E0*I_ERI_G2y2z_Py_Py_S_aa-2.0E0*1*I_ERI_D2y_Py_Py_S_a-2.0E0*1*I_ERI_D2z_Py_Py_S_a+1*I_ERI_S_Py_Py_S;
  abcd[245] = 4.0E0*I_ERI_Gy3z_Py_Py_S_aa-2.0E0*2*I_ERI_Dyz_Py_Py_S_a;
  abcd[246] = 4.0E0*I_ERI_G2xyz_Pz_Py_S_aa;
  abcd[247] = 4.0E0*I_ERI_Gx2yz_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Py_S_a;
  abcd[248] = 4.0E0*I_ERI_Gxy2z_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Py_S_a;
  abcd[249] = 4.0E0*I_ERI_G3yz_Pz_Py_S_aa-2.0E0*2*I_ERI_Dyz_Pz_Py_S_a;
  abcd[250] = 4.0E0*I_ERI_G2y2z_Pz_Py_S_aa-2.0E0*1*I_ERI_D2y_Pz_Py_S_a-2.0E0*1*I_ERI_D2z_Pz_Py_S_a+1*I_ERI_S_Pz_Py_S;
  abcd[251] = 4.0E0*I_ERI_Gy3z_Pz_Py_S_aa-2.0E0*2*I_ERI_Dyz_Pz_Py_S_a;
  abcd[252] = 4.0E0*I_ERI_G2xyz_Px_Pz_S_aa;
  abcd[253] = 4.0E0*I_ERI_Gx2yz_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Px_Pz_S_a;
  abcd[254] = 4.0E0*I_ERI_Gxy2z_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Px_Pz_S_a;
  abcd[255] = 4.0E0*I_ERI_G3yz_Px_Pz_S_aa-2.0E0*2*I_ERI_Dyz_Px_Pz_S_a;
  abcd[256] = 4.0E0*I_ERI_G2y2z_Px_Pz_S_aa-2.0E0*1*I_ERI_D2y_Px_Pz_S_a-2.0E0*1*I_ERI_D2z_Px_Pz_S_a+1*I_ERI_S_Px_Pz_S;
  abcd[257] = 4.0E0*I_ERI_Gy3z_Px_Pz_S_aa-2.0E0*2*I_ERI_Dyz_Px_Pz_S_a;
  abcd[258] = 4.0E0*I_ERI_G2xyz_Py_Pz_S_aa;
  abcd[259] = 4.0E0*I_ERI_Gx2yz_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Py_Pz_S_a;
  abcd[260] = 4.0E0*I_ERI_Gxy2z_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Py_Pz_S_a;
  abcd[261] = 4.0E0*I_ERI_G3yz_Py_Pz_S_aa-2.0E0*2*I_ERI_Dyz_Py_Pz_S_a;
  abcd[262] = 4.0E0*I_ERI_G2y2z_Py_Pz_S_aa-2.0E0*1*I_ERI_D2y_Py_Pz_S_a-2.0E0*1*I_ERI_D2z_Py_Pz_S_a+1*I_ERI_S_Py_Pz_S;
  abcd[263] = 4.0E0*I_ERI_Gy3z_Py_Pz_S_aa-2.0E0*2*I_ERI_Dyz_Py_Pz_S_a;
  abcd[264] = 4.0E0*I_ERI_G2xyz_Pz_Pz_S_aa;
  abcd[265] = 4.0E0*I_ERI_Gx2yz_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_a;
  abcd[266] = 4.0E0*I_ERI_Gxy2z_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_a;
  abcd[267] = 4.0E0*I_ERI_G3yz_Pz_Pz_S_aa-2.0E0*2*I_ERI_Dyz_Pz_Pz_S_a;
  abcd[268] = 4.0E0*I_ERI_G2y2z_Pz_Pz_S_aa-2.0E0*1*I_ERI_D2y_Pz_Pz_S_a-2.0E0*1*I_ERI_D2z_Pz_Pz_S_a+1*I_ERI_S_Pz_Pz_S;
  abcd[269] = 4.0E0*I_ERI_Gy3z_Pz_Pz_S_aa-2.0E0*2*I_ERI_Dyz_Pz_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_aa
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_D_P_P_S_a
   * RHS shell quartet name: SQ_ERI_S_P_P_S
   ************************************************************/
  abcd[270] = 4.0E0*I_ERI_G2x2z_Px_Px_S_aa-2.0E0*1*I_ERI_D2x_Px_Px_S_a;
  abcd[271] = 4.0E0*I_ERI_Gxy2z_Px_Px_S_aa-2.0E0*1*I_ERI_Dxy_Px_Px_S_a;
  abcd[272] = 4.0E0*I_ERI_Gx3z_Px_Px_S_aa-2.0E0*1*I_ERI_Dxz_Px_Px_S_a-2.0E0*2*I_ERI_Dxz_Px_Px_S_a;
  abcd[273] = 4.0E0*I_ERI_G2y2z_Px_Px_S_aa-2.0E0*1*I_ERI_D2y_Px_Px_S_a;
  abcd[274] = 4.0E0*I_ERI_Gy3z_Px_Px_S_aa-2.0E0*1*I_ERI_Dyz_Px_Px_S_a-2.0E0*2*I_ERI_Dyz_Px_Px_S_a;
  abcd[275] = 4.0E0*I_ERI_G4z_Px_Px_S_aa-2.0E0*2*I_ERI_D2z_Px_Px_S_a-2.0E0*3*I_ERI_D2z_Px_Px_S_a+2*1*I_ERI_S_Px_Px_S;
  abcd[276] = 4.0E0*I_ERI_G2x2z_Py_Px_S_aa-2.0E0*1*I_ERI_D2x_Py_Px_S_a;
  abcd[277] = 4.0E0*I_ERI_Gxy2z_Py_Px_S_aa-2.0E0*1*I_ERI_Dxy_Py_Px_S_a;
  abcd[278] = 4.0E0*I_ERI_Gx3z_Py_Px_S_aa-2.0E0*1*I_ERI_Dxz_Py_Px_S_a-2.0E0*2*I_ERI_Dxz_Py_Px_S_a;
  abcd[279] = 4.0E0*I_ERI_G2y2z_Py_Px_S_aa-2.0E0*1*I_ERI_D2y_Py_Px_S_a;
  abcd[280] = 4.0E0*I_ERI_Gy3z_Py_Px_S_aa-2.0E0*1*I_ERI_Dyz_Py_Px_S_a-2.0E0*2*I_ERI_Dyz_Py_Px_S_a;
  abcd[281] = 4.0E0*I_ERI_G4z_Py_Px_S_aa-2.0E0*2*I_ERI_D2z_Py_Px_S_a-2.0E0*3*I_ERI_D2z_Py_Px_S_a+2*1*I_ERI_S_Py_Px_S;
  abcd[282] = 4.0E0*I_ERI_G2x2z_Pz_Px_S_aa-2.0E0*1*I_ERI_D2x_Pz_Px_S_a;
  abcd[283] = 4.0E0*I_ERI_Gxy2z_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Px_S_a;
  abcd[284] = 4.0E0*I_ERI_Gx3z_Pz_Px_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Px_S_a-2.0E0*2*I_ERI_Dxz_Pz_Px_S_a;
  abcd[285] = 4.0E0*I_ERI_G2y2z_Pz_Px_S_aa-2.0E0*1*I_ERI_D2y_Pz_Px_S_a;
  abcd[286] = 4.0E0*I_ERI_Gy3z_Pz_Px_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Px_S_a-2.0E0*2*I_ERI_Dyz_Pz_Px_S_a;
  abcd[287] = 4.0E0*I_ERI_G4z_Pz_Px_S_aa-2.0E0*2*I_ERI_D2z_Pz_Px_S_a-2.0E0*3*I_ERI_D2z_Pz_Px_S_a+2*1*I_ERI_S_Pz_Px_S;
  abcd[288] = 4.0E0*I_ERI_G2x2z_Px_Py_S_aa-2.0E0*1*I_ERI_D2x_Px_Py_S_a;
  abcd[289] = 4.0E0*I_ERI_Gxy2z_Px_Py_S_aa-2.0E0*1*I_ERI_Dxy_Px_Py_S_a;
  abcd[290] = 4.0E0*I_ERI_Gx3z_Px_Py_S_aa-2.0E0*1*I_ERI_Dxz_Px_Py_S_a-2.0E0*2*I_ERI_Dxz_Px_Py_S_a;
  abcd[291] = 4.0E0*I_ERI_G2y2z_Px_Py_S_aa-2.0E0*1*I_ERI_D2y_Px_Py_S_a;
  abcd[292] = 4.0E0*I_ERI_Gy3z_Px_Py_S_aa-2.0E0*1*I_ERI_Dyz_Px_Py_S_a-2.0E0*2*I_ERI_Dyz_Px_Py_S_a;
  abcd[293] = 4.0E0*I_ERI_G4z_Px_Py_S_aa-2.0E0*2*I_ERI_D2z_Px_Py_S_a-2.0E0*3*I_ERI_D2z_Px_Py_S_a+2*1*I_ERI_S_Px_Py_S;
  abcd[294] = 4.0E0*I_ERI_G2x2z_Py_Py_S_aa-2.0E0*1*I_ERI_D2x_Py_Py_S_a;
  abcd[295] = 4.0E0*I_ERI_Gxy2z_Py_Py_S_aa-2.0E0*1*I_ERI_Dxy_Py_Py_S_a;
  abcd[296] = 4.0E0*I_ERI_Gx3z_Py_Py_S_aa-2.0E0*1*I_ERI_Dxz_Py_Py_S_a-2.0E0*2*I_ERI_Dxz_Py_Py_S_a;
  abcd[297] = 4.0E0*I_ERI_G2y2z_Py_Py_S_aa-2.0E0*1*I_ERI_D2y_Py_Py_S_a;
  abcd[298] = 4.0E0*I_ERI_Gy3z_Py_Py_S_aa-2.0E0*1*I_ERI_Dyz_Py_Py_S_a-2.0E0*2*I_ERI_Dyz_Py_Py_S_a;
  abcd[299] = 4.0E0*I_ERI_G4z_Py_Py_S_aa-2.0E0*2*I_ERI_D2z_Py_Py_S_a-2.0E0*3*I_ERI_D2z_Py_Py_S_a+2*1*I_ERI_S_Py_Py_S;
  abcd[300] = 4.0E0*I_ERI_G2x2z_Pz_Py_S_aa-2.0E0*1*I_ERI_D2x_Pz_Py_S_a;
  abcd[301] = 4.0E0*I_ERI_Gxy2z_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Py_S_a;
  abcd[302] = 4.0E0*I_ERI_Gx3z_Pz_Py_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Py_S_a-2.0E0*2*I_ERI_Dxz_Pz_Py_S_a;
  abcd[303] = 4.0E0*I_ERI_G2y2z_Pz_Py_S_aa-2.0E0*1*I_ERI_D2y_Pz_Py_S_a;
  abcd[304] = 4.0E0*I_ERI_Gy3z_Pz_Py_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Py_S_a-2.0E0*2*I_ERI_Dyz_Pz_Py_S_a;
  abcd[305] = 4.0E0*I_ERI_G4z_Pz_Py_S_aa-2.0E0*2*I_ERI_D2z_Pz_Py_S_a-2.0E0*3*I_ERI_D2z_Pz_Py_S_a+2*1*I_ERI_S_Pz_Py_S;
  abcd[306] = 4.0E0*I_ERI_G2x2z_Px_Pz_S_aa-2.0E0*1*I_ERI_D2x_Px_Pz_S_a;
  abcd[307] = 4.0E0*I_ERI_Gxy2z_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Px_Pz_S_a;
  abcd[308] = 4.0E0*I_ERI_Gx3z_Px_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Px_Pz_S_a-2.0E0*2*I_ERI_Dxz_Px_Pz_S_a;
  abcd[309] = 4.0E0*I_ERI_G2y2z_Px_Pz_S_aa-2.0E0*1*I_ERI_D2y_Px_Pz_S_a;
  abcd[310] = 4.0E0*I_ERI_Gy3z_Px_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Px_Pz_S_a-2.0E0*2*I_ERI_Dyz_Px_Pz_S_a;
  abcd[311] = 4.0E0*I_ERI_G4z_Px_Pz_S_aa-2.0E0*2*I_ERI_D2z_Px_Pz_S_a-2.0E0*3*I_ERI_D2z_Px_Pz_S_a+2*1*I_ERI_S_Px_Pz_S;
  abcd[312] = 4.0E0*I_ERI_G2x2z_Py_Pz_S_aa-2.0E0*1*I_ERI_D2x_Py_Pz_S_a;
  abcd[313] = 4.0E0*I_ERI_Gxy2z_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Py_Pz_S_a;
  abcd[314] = 4.0E0*I_ERI_Gx3z_Py_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Py_Pz_S_a-2.0E0*2*I_ERI_Dxz_Py_Pz_S_a;
  abcd[315] = 4.0E0*I_ERI_G2y2z_Py_Pz_S_aa-2.0E0*1*I_ERI_D2y_Py_Pz_S_a;
  abcd[316] = 4.0E0*I_ERI_Gy3z_Py_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Py_Pz_S_a-2.0E0*2*I_ERI_Dyz_Py_Pz_S_a;
  abcd[317] = 4.0E0*I_ERI_G4z_Py_Pz_S_aa-2.0E0*2*I_ERI_D2z_Py_Pz_S_a-2.0E0*3*I_ERI_D2z_Py_Pz_S_a+2*1*I_ERI_S_Py_Pz_S;
  abcd[318] = 4.0E0*I_ERI_G2x2z_Pz_Pz_S_aa-2.0E0*1*I_ERI_D2x_Pz_Pz_S_a;
  abcd[319] = 4.0E0*I_ERI_Gxy2z_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_a;
  abcd[320] = 4.0E0*I_ERI_Gx3z_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_a-2.0E0*2*I_ERI_Dxz_Pz_Pz_S_a;
  abcd[321] = 4.0E0*I_ERI_G2y2z_Pz_Pz_S_aa-2.0E0*1*I_ERI_D2y_Pz_Pz_S_a;
  abcd[322] = 4.0E0*I_ERI_Gy3z_Pz_Pz_S_aa-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_a-2.0E0*2*I_ERI_Dyz_Pz_Pz_S_a;
  abcd[323] = 4.0E0*I_ERI_G4z_Pz_Pz_S_aa-2.0E0*2*I_ERI_D2z_Pz_Pz_S_a-2.0E0*3*I_ERI_D2z_Pz_Pz_S_a+2*1*I_ERI_S_Pz_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_D_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[324] = 4.0E0*I_ERI_F3x_D2x_Px_S_ab-2.0E0*1*I_ERI_F3x_S_Px_S_a-2.0E0*2*I_ERI_Px_D2x_Px_S_b+2*1*I_ERI_Px_S_Px_S;
  abcd[325] = 4.0E0*I_ERI_F2xy_D2x_Px_S_ab-2.0E0*1*I_ERI_F2xy_S_Px_S_a-2.0E0*1*I_ERI_Py_D2x_Px_S_b+1*I_ERI_Py_S_Px_S;
  abcd[326] = 4.0E0*I_ERI_F2xz_D2x_Px_S_ab-2.0E0*1*I_ERI_F2xz_S_Px_S_a-2.0E0*1*I_ERI_Pz_D2x_Px_S_b+1*I_ERI_Pz_S_Px_S;
  abcd[327] = 4.0E0*I_ERI_Fx2y_D2x_Px_S_ab-2.0E0*1*I_ERI_Fx2y_S_Px_S_a;
  abcd[328] = 4.0E0*I_ERI_Fxyz_D2x_Px_S_ab-2.0E0*1*I_ERI_Fxyz_S_Px_S_a;
  abcd[329] = 4.0E0*I_ERI_Fx2z_D2x_Px_S_ab-2.0E0*1*I_ERI_Fx2z_S_Px_S_a;
  abcd[330] = 4.0E0*I_ERI_F3x_Dxy_Px_S_ab-2.0E0*2*I_ERI_Px_Dxy_Px_S_b;
  abcd[331] = 4.0E0*I_ERI_F2xy_Dxy_Px_S_ab-2.0E0*1*I_ERI_Py_Dxy_Px_S_b;
  abcd[332] = 4.0E0*I_ERI_F2xz_Dxy_Px_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Px_S_b;
  abcd[333] = 4.0E0*I_ERI_Fx2y_Dxy_Px_S_ab;
  abcd[334] = 4.0E0*I_ERI_Fxyz_Dxy_Px_S_ab;
  abcd[335] = 4.0E0*I_ERI_Fx2z_Dxy_Px_S_ab;
  abcd[336] = 4.0E0*I_ERI_F3x_Dxz_Px_S_ab-2.0E0*2*I_ERI_Px_Dxz_Px_S_b;
  abcd[337] = 4.0E0*I_ERI_F2xy_Dxz_Px_S_ab-2.0E0*1*I_ERI_Py_Dxz_Px_S_b;
  abcd[338] = 4.0E0*I_ERI_F2xz_Dxz_Px_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Px_S_b;
  abcd[339] = 4.0E0*I_ERI_Fx2y_Dxz_Px_S_ab;
  abcd[340] = 4.0E0*I_ERI_Fxyz_Dxz_Px_S_ab;
  abcd[341] = 4.0E0*I_ERI_Fx2z_Dxz_Px_S_ab;
  abcd[342] = 4.0E0*I_ERI_F3x_D2x_Py_S_ab-2.0E0*1*I_ERI_F3x_S_Py_S_a-2.0E0*2*I_ERI_Px_D2x_Py_S_b+2*1*I_ERI_Px_S_Py_S;
  abcd[343] = 4.0E0*I_ERI_F2xy_D2x_Py_S_ab-2.0E0*1*I_ERI_F2xy_S_Py_S_a-2.0E0*1*I_ERI_Py_D2x_Py_S_b+1*I_ERI_Py_S_Py_S;
  abcd[344] = 4.0E0*I_ERI_F2xz_D2x_Py_S_ab-2.0E0*1*I_ERI_F2xz_S_Py_S_a-2.0E0*1*I_ERI_Pz_D2x_Py_S_b+1*I_ERI_Pz_S_Py_S;
  abcd[345] = 4.0E0*I_ERI_Fx2y_D2x_Py_S_ab-2.0E0*1*I_ERI_Fx2y_S_Py_S_a;
  abcd[346] = 4.0E0*I_ERI_Fxyz_D2x_Py_S_ab-2.0E0*1*I_ERI_Fxyz_S_Py_S_a;
  abcd[347] = 4.0E0*I_ERI_Fx2z_D2x_Py_S_ab-2.0E0*1*I_ERI_Fx2z_S_Py_S_a;
  abcd[348] = 4.0E0*I_ERI_F3x_Dxy_Py_S_ab-2.0E0*2*I_ERI_Px_Dxy_Py_S_b;
  abcd[349] = 4.0E0*I_ERI_F2xy_Dxy_Py_S_ab-2.0E0*1*I_ERI_Py_Dxy_Py_S_b;
  abcd[350] = 4.0E0*I_ERI_F2xz_Dxy_Py_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Py_S_b;
  abcd[351] = 4.0E0*I_ERI_Fx2y_Dxy_Py_S_ab;
  abcd[352] = 4.0E0*I_ERI_Fxyz_Dxy_Py_S_ab;
  abcd[353] = 4.0E0*I_ERI_Fx2z_Dxy_Py_S_ab;
  abcd[354] = 4.0E0*I_ERI_F3x_Dxz_Py_S_ab-2.0E0*2*I_ERI_Px_Dxz_Py_S_b;
  abcd[355] = 4.0E0*I_ERI_F2xy_Dxz_Py_S_ab-2.0E0*1*I_ERI_Py_Dxz_Py_S_b;
  abcd[356] = 4.0E0*I_ERI_F2xz_Dxz_Py_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Py_S_b;
  abcd[357] = 4.0E0*I_ERI_Fx2y_Dxz_Py_S_ab;
  abcd[358] = 4.0E0*I_ERI_Fxyz_Dxz_Py_S_ab;
  abcd[359] = 4.0E0*I_ERI_Fx2z_Dxz_Py_S_ab;
  abcd[360] = 4.0E0*I_ERI_F3x_D2x_Pz_S_ab-2.0E0*1*I_ERI_F3x_S_Pz_S_a-2.0E0*2*I_ERI_Px_D2x_Pz_S_b+2*1*I_ERI_Px_S_Pz_S;
  abcd[361] = 4.0E0*I_ERI_F2xy_D2x_Pz_S_ab-2.0E0*1*I_ERI_F2xy_S_Pz_S_a-2.0E0*1*I_ERI_Py_D2x_Pz_S_b+1*I_ERI_Py_S_Pz_S;
  abcd[362] = 4.0E0*I_ERI_F2xz_D2x_Pz_S_ab-2.0E0*1*I_ERI_F2xz_S_Pz_S_a-2.0E0*1*I_ERI_Pz_D2x_Pz_S_b+1*I_ERI_Pz_S_Pz_S;
  abcd[363] = 4.0E0*I_ERI_Fx2y_D2x_Pz_S_ab-2.0E0*1*I_ERI_Fx2y_S_Pz_S_a;
  abcd[364] = 4.0E0*I_ERI_Fxyz_D2x_Pz_S_ab-2.0E0*1*I_ERI_Fxyz_S_Pz_S_a;
  abcd[365] = 4.0E0*I_ERI_Fx2z_D2x_Pz_S_ab-2.0E0*1*I_ERI_Fx2z_S_Pz_S_a;
  abcd[366] = 4.0E0*I_ERI_F3x_Dxy_Pz_S_ab-2.0E0*2*I_ERI_Px_Dxy_Pz_S_b;
  abcd[367] = 4.0E0*I_ERI_F2xy_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Py_Dxy_Pz_S_b;
  abcd[368] = 4.0E0*I_ERI_F2xz_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Pz_S_b;
  abcd[369] = 4.0E0*I_ERI_Fx2y_Dxy_Pz_S_ab;
  abcd[370] = 4.0E0*I_ERI_Fxyz_Dxy_Pz_S_ab;
  abcd[371] = 4.0E0*I_ERI_Fx2z_Dxy_Pz_S_ab;
  abcd[372] = 4.0E0*I_ERI_F3x_Dxz_Pz_S_ab-2.0E0*2*I_ERI_Px_Dxz_Pz_S_b;
  abcd[373] = 4.0E0*I_ERI_F2xy_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Py_Dxz_Pz_S_b;
  abcd[374] = 4.0E0*I_ERI_F2xz_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Pz_S_b;
  abcd[375] = 4.0E0*I_ERI_Fx2y_Dxz_Pz_S_ab;
  abcd[376] = 4.0E0*I_ERI_Fxyz_Dxz_Pz_S_ab;
  abcd[377] = 4.0E0*I_ERI_Fx2z_Dxz_Pz_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_D_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[378] = 4.0E0*I_ERI_F3x_Dxy_Px_S_ab-2.0E0*2*I_ERI_Px_Dxy_Px_S_b;
  abcd[379] = 4.0E0*I_ERI_F2xy_Dxy_Px_S_ab-2.0E0*1*I_ERI_Py_Dxy_Px_S_b;
  abcd[380] = 4.0E0*I_ERI_F2xz_Dxy_Px_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Px_S_b;
  abcd[381] = 4.0E0*I_ERI_Fx2y_Dxy_Px_S_ab;
  abcd[382] = 4.0E0*I_ERI_Fxyz_Dxy_Px_S_ab;
  abcd[383] = 4.0E0*I_ERI_Fx2z_Dxy_Px_S_ab;
  abcd[384] = 4.0E0*I_ERI_F3x_D2y_Px_S_ab-2.0E0*1*I_ERI_F3x_S_Px_S_a-2.0E0*2*I_ERI_Px_D2y_Px_S_b+2*1*I_ERI_Px_S_Px_S;
  abcd[385] = 4.0E0*I_ERI_F2xy_D2y_Px_S_ab-2.0E0*1*I_ERI_F2xy_S_Px_S_a-2.0E0*1*I_ERI_Py_D2y_Px_S_b+1*I_ERI_Py_S_Px_S;
  abcd[386] = 4.0E0*I_ERI_F2xz_D2y_Px_S_ab-2.0E0*1*I_ERI_F2xz_S_Px_S_a-2.0E0*1*I_ERI_Pz_D2y_Px_S_b+1*I_ERI_Pz_S_Px_S;
  abcd[387] = 4.0E0*I_ERI_Fx2y_D2y_Px_S_ab-2.0E0*1*I_ERI_Fx2y_S_Px_S_a;
  abcd[388] = 4.0E0*I_ERI_Fxyz_D2y_Px_S_ab-2.0E0*1*I_ERI_Fxyz_S_Px_S_a;
  abcd[389] = 4.0E0*I_ERI_Fx2z_D2y_Px_S_ab-2.0E0*1*I_ERI_Fx2z_S_Px_S_a;
  abcd[390] = 4.0E0*I_ERI_F3x_Dyz_Px_S_ab-2.0E0*2*I_ERI_Px_Dyz_Px_S_b;
  abcd[391] = 4.0E0*I_ERI_F2xy_Dyz_Px_S_ab-2.0E0*1*I_ERI_Py_Dyz_Px_S_b;
  abcd[392] = 4.0E0*I_ERI_F2xz_Dyz_Px_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Px_S_b;
  abcd[393] = 4.0E0*I_ERI_Fx2y_Dyz_Px_S_ab;
  abcd[394] = 4.0E0*I_ERI_Fxyz_Dyz_Px_S_ab;
  abcd[395] = 4.0E0*I_ERI_Fx2z_Dyz_Px_S_ab;
  abcd[396] = 4.0E0*I_ERI_F3x_Dxy_Py_S_ab-2.0E0*2*I_ERI_Px_Dxy_Py_S_b;
  abcd[397] = 4.0E0*I_ERI_F2xy_Dxy_Py_S_ab-2.0E0*1*I_ERI_Py_Dxy_Py_S_b;
  abcd[398] = 4.0E0*I_ERI_F2xz_Dxy_Py_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Py_S_b;
  abcd[399] = 4.0E0*I_ERI_Fx2y_Dxy_Py_S_ab;
  abcd[400] = 4.0E0*I_ERI_Fxyz_Dxy_Py_S_ab;
  abcd[401] = 4.0E0*I_ERI_Fx2z_Dxy_Py_S_ab;
  abcd[402] = 4.0E0*I_ERI_F3x_D2y_Py_S_ab-2.0E0*1*I_ERI_F3x_S_Py_S_a-2.0E0*2*I_ERI_Px_D2y_Py_S_b+2*1*I_ERI_Px_S_Py_S;
  abcd[403] = 4.0E0*I_ERI_F2xy_D2y_Py_S_ab-2.0E0*1*I_ERI_F2xy_S_Py_S_a-2.0E0*1*I_ERI_Py_D2y_Py_S_b+1*I_ERI_Py_S_Py_S;
  abcd[404] = 4.0E0*I_ERI_F2xz_D2y_Py_S_ab-2.0E0*1*I_ERI_F2xz_S_Py_S_a-2.0E0*1*I_ERI_Pz_D2y_Py_S_b+1*I_ERI_Pz_S_Py_S;
  abcd[405] = 4.0E0*I_ERI_Fx2y_D2y_Py_S_ab-2.0E0*1*I_ERI_Fx2y_S_Py_S_a;
  abcd[406] = 4.0E0*I_ERI_Fxyz_D2y_Py_S_ab-2.0E0*1*I_ERI_Fxyz_S_Py_S_a;
  abcd[407] = 4.0E0*I_ERI_Fx2z_D2y_Py_S_ab-2.0E0*1*I_ERI_Fx2z_S_Py_S_a;
  abcd[408] = 4.0E0*I_ERI_F3x_Dyz_Py_S_ab-2.0E0*2*I_ERI_Px_Dyz_Py_S_b;
  abcd[409] = 4.0E0*I_ERI_F2xy_Dyz_Py_S_ab-2.0E0*1*I_ERI_Py_Dyz_Py_S_b;
  abcd[410] = 4.0E0*I_ERI_F2xz_Dyz_Py_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Py_S_b;
  abcd[411] = 4.0E0*I_ERI_Fx2y_Dyz_Py_S_ab;
  abcd[412] = 4.0E0*I_ERI_Fxyz_Dyz_Py_S_ab;
  abcd[413] = 4.0E0*I_ERI_Fx2z_Dyz_Py_S_ab;
  abcd[414] = 4.0E0*I_ERI_F3x_Dxy_Pz_S_ab-2.0E0*2*I_ERI_Px_Dxy_Pz_S_b;
  abcd[415] = 4.0E0*I_ERI_F2xy_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Py_Dxy_Pz_S_b;
  abcd[416] = 4.0E0*I_ERI_F2xz_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Pz_S_b;
  abcd[417] = 4.0E0*I_ERI_Fx2y_Dxy_Pz_S_ab;
  abcd[418] = 4.0E0*I_ERI_Fxyz_Dxy_Pz_S_ab;
  abcd[419] = 4.0E0*I_ERI_Fx2z_Dxy_Pz_S_ab;
  abcd[420] = 4.0E0*I_ERI_F3x_D2y_Pz_S_ab-2.0E0*1*I_ERI_F3x_S_Pz_S_a-2.0E0*2*I_ERI_Px_D2y_Pz_S_b+2*1*I_ERI_Px_S_Pz_S;
  abcd[421] = 4.0E0*I_ERI_F2xy_D2y_Pz_S_ab-2.0E0*1*I_ERI_F2xy_S_Pz_S_a-2.0E0*1*I_ERI_Py_D2y_Pz_S_b+1*I_ERI_Py_S_Pz_S;
  abcd[422] = 4.0E0*I_ERI_F2xz_D2y_Pz_S_ab-2.0E0*1*I_ERI_F2xz_S_Pz_S_a-2.0E0*1*I_ERI_Pz_D2y_Pz_S_b+1*I_ERI_Pz_S_Pz_S;
  abcd[423] = 4.0E0*I_ERI_Fx2y_D2y_Pz_S_ab-2.0E0*1*I_ERI_Fx2y_S_Pz_S_a;
  abcd[424] = 4.0E0*I_ERI_Fxyz_D2y_Pz_S_ab-2.0E0*1*I_ERI_Fxyz_S_Pz_S_a;
  abcd[425] = 4.0E0*I_ERI_Fx2z_D2y_Pz_S_ab-2.0E0*1*I_ERI_Fx2z_S_Pz_S_a;
  abcd[426] = 4.0E0*I_ERI_F3x_Dyz_Pz_S_ab-2.0E0*2*I_ERI_Px_Dyz_Pz_S_b;
  abcd[427] = 4.0E0*I_ERI_F2xy_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Py_Dyz_Pz_S_b;
  abcd[428] = 4.0E0*I_ERI_F2xz_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Pz_S_b;
  abcd[429] = 4.0E0*I_ERI_Fx2y_Dyz_Pz_S_ab;
  abcd[430] = 4.0E0*I_ERI_Fxyz_Dyz_Pz_S_ab;
  abcd[431] = 4.0E0*I_ERI_Fx2z_Dyz_Pz_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_D_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[432] = 4.0E0*I_ERI_F3x_Dxz_Px_S_ab-2.0E0*2*I_ERI_Px_Dxz_Px_S_b;
  abcd[433] = 4.0E0*I_ERI_F2xy_Dxz_Px_S_ab-2.0E0*1*I_ERI_Py_Dxz_Px_S_b;
  abcd[434] = 4.0E0*I_ERI_F2xz_Dxz_Px_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Px_S_b;
  abcd[435] = 4.0E0*I_ERI_Fx2y_Dxz_Px_S_ab;
  abcd[436] = 4.0E0*I_ERI_Fxyz_Dxz_Px_S_ab;
  abcd[437] = 4.0E0*I_ERI_Fx2z_Dxz_Px_S_ab;
  abcd[438] = 4.0E0*I_ERI_F3x_Dyz_Px_S_ab-2.0E0*2*I_ERI_Px_Dyz_Px_S_b;
  abcd[439] = 4.0E0*I_ERI_F2xy_Dyz_Px_S_ab-2.0E0*1*I_ERI_Py_Dyz_Px_S_b;
  abcd[440] = 4.0E0*I_ERI_F2xz_Dyz_Px_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Px_S_b;
  abcd[441] = 4.0E0*I_ERI_Fx2y_Dyz_Px_S_ab;
  abcd[442] = 4.0E0*I_ERI_Fxyz_Dyz_Px_S_ab;
  abcd[443] = 4.0E0*I_ERI_Fx2z_Dyz_Px_S_ab;
  abcd[444] = 4.0E0*I_ERI_F3x_D2z_Px_S_ab-2.0E0*1*I_ERI_F3x_S_Px_S_a-2.0E0*2*I_ERI_Px_D2z_Px_S_b+2*1*I_ERI_Px_S_Px_S;
  abcd[445] = 4.0E0*I_ERI_F2xy_D2z_Px_S_ab-2.0E0*1*I_ERI_F2xy_S_Px_S_a-2.0E0*1*I_ERI_Py_D2z_Px_S_b+1*I_ERI_Py_S_Px_S;
  abcd[446] = 4.0E0*I_ERI_F2xz_D2z_Px_S_ab-2.0E0*1*I_ERI_F2xz_S_Px_S_a-2.0E0*1*I_ERI_Pz_D2z_Px_S_b+1*I_ERI_Pz_S_Px_S;
  abcd[447] = 4.0E0*I_ERI_Fx2y_D2z_Px_S_ab-2.0E0*1*I_ERI_Fx2y_S_Px_S_a;
  abcd[448] = 4.0E0*I_ERI_Fxyz_D2z_Px_S_ab-2.0E0*1*I_ERI_Fxyz_S_Px_S_a;
  abcd[449] = 4.0E0*I_ERI_Fx2z_D2z_Px_S_ab-2.0E0*1*I_ERI_Fx2z_S_Px_S_a;
  abcd[450] = 4.0E0*I_ERI_F3x_Dxz_Py_S_ab-2.0E0*2*I_ERI_Px_Dxz_Py_S_b;
  abcd[451] = 4.0E0*I_ERI_F2xy_Dxz_Py_S_ab-2.0E0*1*I_ERI_Py_Dxz_Py_S_b;
  abcd[452] = 4.0E0*I_ERI_F2xz_Dxz_Py_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Py_S_b;
  abcd[453] = 4.0E0*I_ERI_Fx2y_Dxz_Py_S_ab;
  abcd[454] = 4.0E0*I_ERI_Fxyz_Dxz_Py_S_ab;
  abcd[455] = 4.0E0*I_ERI_Fx2z_Dxz_Py_S_ab;
  abcd[456] = 4.0E0*I_ERI_F3x_Dyz_Py_S_ab-2.0E0*2*I_ERI_Px_Dyz_Py_S_b;
  abcd[457] = 4.0E0*I_ERI_F2xy_Dyz_Py_S_ab-2.0E0*1*I_ERI_Py_Dyz_Py_S_b;
  abcd[458] = 4.0E0*I_ERI_F2xz_Dyz_Py_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Py_S_b;
  abcd[459] = 4.0E0*I_ERI_Fx2y_Dyz_Py_S_ab;
  abcd[460] = 4.0E0*I_ERI_Fxyz_Dyz_Py_S_ab;
  abcd[461] = 4.0E0*I_ERI_Fx2z_Dyz_Py_S_ab;
  abcd[462] = 4.0E0*I_ERI_F3x_D2z_Py_S_ab-2.0E0*1*I_ERI_F3x_S_Py_S_a-2.0E0*2*I_ERI_Px_D2z_Py_S_b+2*1*I_ERI_Px_S_Py_S;
  abcd[463] = 4.0E0*I_ERI_F2xy_D2z_Py_S_ab-2.0E0*1*I_ERI_F2xy_S_Py_S_a-2.0E0*1*I_ERI_Py_D2z_Py_S_b+1*I_ERI_Py_S_Py_S;
  abcd[464] = 4.0E0*I_ERI_F2xz_D2z_Py_S_ab-2.0E0*1*I_ERI_F2xz_S_Py_S_a-2.0E0*1*I_ERI_Pz_D2z_Py_S_b+1*I_ERI_Pz_S_Py_S;
  abcd[465] = 4.0E0*I_ERI_Fx2y_D2z_Py_S_ab-2.0E0*1*I_ERI_Fx2y_S_Py_S_a;
  abcd[466] = 4.0E0*I_ERI_Fxyz_D2z_Py_S_ab-2.0E0*1*I_ERI_Fxyz_S_Py_S_a;
  abcd[467] = 4.0E0*I_ERI_Fx2z_D2z_Py_S_ab-2.0E0*1*I_ERI_Fx2z_S_Py_S_a;
  abcd[468] = 4.0E0*I_ERI_F3x_Dxz_Pz_S_ab-2.0E0*2*I_ERI_Px_Dxz_Pz_S_b;
  abcd[469] = 4.0E0*I_ERI_F2xy_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Py_Dxz_Pz_S_b;
  abcd[470] = 4.0E0*I_ERI_F2xz_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Pz_S_b;
  abcd[471] = 4.0E0*I_ERI_Fx2y_Dxz_Pz_S_ab;
  abcd[472] = 4.0E0*I_ERI_Fxyz_Dxz_Pz_S_ab;
  abcd[473] = 4.0E0*I_ERI_Fx2z_Dxz_Pz_S_ab;
  abcd[474] = 4.0E0*I_ERI_F3x_Dyz_Pz_S_ab-2.0E0*2*I_ERI_Px_Dyz_Pz_S_b;
  abcd[475] = 4.0E0*I_ERI_F2xy_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Py_Dyz_Pz_S_b;
  abcd[476] = 4.0E0*I_ERI_F2xz_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Pz_S_b;
  abcd[477] = 4.0E0*I_ERI_Fx2y_Dyz_Pz_S_ab;
  abcd[478] = 4.0E0*I_ERI_Fxyz_Dyz_Pz_S_ab;
  abcd[479] = 4.0E0*I_ERI_Fx2z_Dyz_Pz_S_ab;
  abcd[480] = 4.0E0*I_ERI_F3x_D2z_Pz_S_ab-2.0E0*1*I_ERI_F3x_S_Pz_S_a-2.0E0*2*I_ERI_Px_D2z_Pz_S_b+2*1*I_ERI_Px_S_Pz_S;
  abcd[481] = 4.0E0*I_ERI_F2xy_D2z_Pz_S_ab-2.0E0*1*I_ERI_F2xy_S_Pz_S_a-2.0E0*1*I_ERI_Py_D2z_Pz_S_b+1*I_ERI_Py_S_Pz_S;
  abcd[482] = 4.0E0*I_ERI_F2xz_D2z_Pz_S_ab-2.0E0*1*I_ERI_F2xz_S_Pz_S_a-2.0E0*1*I_ERI_Pz_D2z_Pz_S_b+1*I_ERI_Pz_S_Pz_S;
  abcd[483] = 4.0E0*I_ERI_Fx2y_D2z_Pz_S_ab-2.0E0*1*I_ERI_Fx2y_S_Pz_S_a;
  abcd[484] = 4.0E0*I_ERI_Fxyz_D2z_Pz_S_ab-2.0E0*1*I_ERI_Fxyz_S_Pz_S_a;
  abcd[485] = 4.0E0*I_ERI_Fx2z_D2z_Pz_S_ab-2.0E0*1*I_ERI_Fx2z_S_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_D_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[486] = 4.0E0*I_ERI_F2xy_D2x_Px_S_ab-2.0E0*1*I_ERI_F2xy_S_Px_S_a;
  abcd[487] = 4.0E0*I_ERI_Fx2y_D2x_Px_S_ab-2.0E0*1*I_ERI_Fx2y_S_Px_S_a-2.0E0*1*I_ERI_Px_D2x_Px_S_b+1*I_ERI_Px_S_Px_S;
  abcd[488] = 4.0E0*I_ERI_Fxyz_D2x_Px_S_ab-2.0E0*1*I_ERI_Fxyz_S_Px_S_a;
  abcd[489] = 4.0E0*I_ERI_F3y_D2x_Px_S_ab-2.0E0*1*I_ERI_F3y_S_Px_S_a-2.0E0*2*I_ERI_Py_D2x_Px_S_b+2*1*I_ERI_Py_S_Px_S;
  abcd[490] = 4.0E0*I_ERI_F2yz_D2x_Px_S_ab-2.0E0*1*I_ERI_F2yz_S_Px_S_a-2.0E0*1*I_ERI_Pz_D2x_Px_S_b+1*I_ERI_Pz_S_Px_S;
  abcd[491] = 4.0E0*I_ERI_Fy2z_D2x_Px_S_ab-2.0E0*1*I_ERI_Fy2z_S_Px_S_a;
  abcd[492] = 4.0E0*I_ERI_F2xy_Dxy_Px_S_ab;
  abcd[493] = 4.0E0*I_ERI_Fx2y_Dxy_Px_S_ab-2.0E0*1*I_ERI_Px_Dxy_Px_S_b;
  abcd[494] = 4.0E0*I_ERI_Fxyz_Dxy_Px_S_ab;
  abcd[495] = 4.0E0*I_ERI_F3y_Dxy_Px_S_ab-2.0E0*2*I_ERI_Py_Dxy_Px_S_b;
  abcd[496] = 4.0E0*I_ERI_F2yz_Dxy_Px_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Px_S_b;
  abcd[497] = 4.0E0*I_ERI_Fy2z_Dxy_Px_S_ab;
  abcd[498] = 4.0E0*I_ERI_F2xy_Dxz_Px_S_ab;
  abcd[499] = 4.0E0*I_ERI_Fx2y_Dxz_Px_S_ab-2.0E0*1*I_ERI_Px_Dxz_Px_S_b;
  abcd[500] = 4.0E0*I_ERI_Fxyz_Dxz_Px_S_ab;
  abcd[501] = 4.0E0*I_ERI_F3y_Dxz_Px_S_ab-2.0E0*2*I_ERI_Py_Dxz_Px_S_b;
  abcd[502] = 4.0E0*I_ERI_F2yz_Dxz_Px_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Px_S_b;
  abcd[503] = 4.0E0*I_ERI_Fy2z_Dxz_Px_S_ab;
  abcd[504] = 4.0E0*I_ERI_F2xy_D2x_Py_S_ab-2.0E0*1*I_ERI_F2xy_S_Py_S_a;
  abcd[505] = 4.0E0*I_ERI_Fx2y_D2x_Py_S_ab-2.0E0*1*I_ERI_Fx2y_S_Py_S_a-2.0E0*1*I_ERI_Px_D2x_Py_S_b+1*I_ERI_Px_S_Py_S;
  abcd[506] = 4.0E0*I_ERI_Fxyz_D2x_Py_S_ab-2.0E0*1*I_ERI_Fxyz_S_Py_S_a;
  abcd[507] = 4.0E0*I_ERI_F3y_D2x_Py_S_ab-2.0E0*1*I_ERI_F3y_S_Py_S_a-2.0E0*2*I_ERI_Py_D2x_Py_S_b+2*1*I_ERI_Py_S_Py_S;
  abcd[508] = 4.0E0*I_ERI_F2yz_D2x_Py_S_ab-2.0E0*1*I_ERI_F2yz_S_Py_S_a-2.0E0*1*I_ERI_Pz_D2x_Py_S_b+1*I_ERI_Pz_S_Py_S;
  abcd[509] = 4.0E0*I_ERI_Fy2z_D2x_Py_S_ab-2.0E0*1*I_ERI_Fy2z_S_Py_S_a;
  abcd[510] = 4.0E0*I_ERI_F2xy_Dxy_Py_S_ab;
  abcd[511] = 4.0E0*I_ERI_Fx2y_Dxy_Py_S_ab-2.0E0*1*I_ERI_Px_Dxy_Py_S_b;
  abcd[512] = 4.0E0*I_ERI_Fxyz_Dxy_Py_S_ab;
  abcd[513] = 4.0E0*I_ERI_F3y_Dxy_Py_S_ab-2.0E0*2*I_ERI_Py_Dxy_Py_S_b;
  abcd[514] = 4.0E0*I_ERI_F2yz_Dxy_Py_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Py_S_b;
  abcd[515] = 4.0E0*I_ERI_Fy2z_Dxy_Py_S_ab;
  abcd[516] = 4.0E0*I_ERI_F2xy_Dxz_Py_S_ab;
  abcd[517] = 4.0E0*I_ERI_Fx2y_Dxz_Py_S_ab-2.0E0*1*I_ERI_Px_Dxz_Py_S_b;
  abcd[518] = 4.0E0*I_ERI_Fxyz_Dxz_Py_S_ab;
  abcd[519] = 4.0E0*I_ERI_F3y_Dxz_Py_S_ab-2.0E0*2*I_ERI_Py_Dxz_Py_S_b;
  abcd[520] = 4.0E0*I_ERI_F2yz_Dxz_Py_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Py_S_b;
  abcd[521] = 4.0E0*I_ERI_Fy2z_Dxz_Py_S_ab;
  abcd[522] = 4.0E0*I_ERI_F2xy_D2x_Pz_S_ab-2.0E0*1*I_ERI_F2xy_S_Pz_S_a;
  abcd[523] = 4.0E0*I_ERI_Fx2y_D2x_Pz_S_ab-2.0E0*1*I_ERI_Fx2y_S_Pz_S_a-2.0E0*1*I_ERI_Px_D2x_Pz_S_b+1*I_ERI_Px_S_Pz_S;
  abcd[524] = 4.0E0*I_ERI_Fxyz_D2x_Pz_S_ab-2.0E0*1*I_ERI_Fxyz_S_Pz_S_a;
  abcd[525] = 4.0E0*I_ERI_F3y_D2x_Pz_S_ab-2.0E0*1*I_ERI_F3y_S_Pz_S_a-2.0E0*2*I_ERI_Py_D2x_Pz_S_b+2*1*I_ERI_Py_S_Pz_S;
  abcd[526] = 4.0E0*I_ERI_F2yz_D2x_Pz_S_ab-2.0E0*1*I_ERI_F2yz_S_Pz_S_a-2.0E0*1*I_ERI_Pz_D2x_Pz_S_b+1*I_ERI_Pz_S_Pz_S;
  abcd[527] = 4.0E0*I_ERI_Fy2z_D2x_Pz_S_ab-2.0E0*1*I_ERI_Fy2z_S_Pz_S_a;
  abcd[528] = 4.0E0*I_ERI_F2xy_Dxy_Pz_S_ab;
  abcd[529] = 4.0E0*I_ERI_Fx2y_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Px_Dxy_Pz_S_b;
  abcd[530] = 4.0E0*I_ERI_Fxyz_Dxy_Pz_S_ab;
  abcd[531] = 4.0E0*I_ERI_F3y_Dxy_Pz_S_ab-2.0E0*2*I_ERI_Py_Dxy_Pz_S_b;
  abcd[532] = 4.0E0*I_ERI_F2yz_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Pz_S_b;
  abcd[533] = 4.0E0*I_ERI_Fy2z_Dxy_Pz_S_ab;
  abcd[534] = 4.0E0*I_ERI_F2xy_Dxz_Pz_S_ab;
  abcd[535] = 4.0E0*I_ERI_Fx2y_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Px_Dxz_Pz_S_b;
  abcd[536] = 4.0E0*I_ERI_Fxyz_Dxz_Pz_S_ab;
  abcd[537] = 4.0E0*I_ERI_F3y_Dxz_Pz_S_ab-2.0E0*2*I_ERI_Py_Dxz_Pz_S_b;
  abcd[538] = 4.0E0*I_ERI_F2yz_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Pz_S_b;
  abcd[539] = 4.0E0*I_ERI_Fy2z_Dxz_Pz_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_D_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[540] = 4.0E0*I_ERI_F2xy_Dxy_Px_S_ab;
  abcd[541] = 4.0E0*I_ERI_Fx2y_Dxy_Px_S_ab-2.0E0*1*I_ERI_Px_Dxy_Px_S_b;
  abcd[542] = 4.0E0*I_ERI_Fxyz_Dxy_Px_S_ab;
  abcd[543] = 4.0E0*I_ERI_F3y_Dxy_Px_S_ab-2.0E0*2*I_ERI_Py_Dxy_Px_S_b;
  abcd[544] = 4.0E0*I_ERI_F2yz_Dxy_Px_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Px_S_b;
  abcd[545] = 4.0E0*I_ERI_Fy2z_Dxy_Px_S_ab;
  abcd[546] = 4.0E0*I_ERI_F2xy_D2y_Px_S_ab-2.0E0*1*I_ERI_F2xy_S_Px_S_a;
  abcd[547] = 4.0E0*I_ERI_Fx2y_D2y_Px_S_ab-2.0E0*1*I_ERI_Fx2y_S_Px_S_a-2.0E0*1*I_ERI_Px_D2y_Px_S_b+1*I_ERI_Px_S_Px_S;
  abcd[548] = 4.0E0*I_ERI_Fxyz_D2y_Px_S_ab-2.0E0*1*I_ERI_Fxyz_S_Px_S_a;
  abcd[549] = 4.0E0*I_ERI_F3y_D2y_Px_S_ab-2.0E0*1*I_ERI_F3y_S_Px_S_a-2.0E0*2*I_ERI_Py_D2y_Px_S_b+2*1*I_ERI_Py_S_Px_S;
  abcd[550] = 4.0E0*I_ERI_F2yz_D2y_Px_S_ab-2.0E0*1*I_ERI_F2yz_S_Px_S_a-2.0E0*1*I_ERI_Pz_D2y_Px_S_b+1*I_ERI_Pz_S_Px_S;
  abcd[551] = 4.0E0*I_ERI_Fy2z_D2y_Px_S_ab-2.0E0*1*I_ERI_Fy2z_S_Px_S_a;
  abcd[552] = 4.0E0*I_ERI_F2xy_Dyz_Px_S_ab;
  abcd[553] = 4.0E0*I_ERI_Fx2y_Dyz_Px_S_ab-2.0E0*1*I_ERI_Px_Dyz_Px_S_b;
  abcd[554] = 4.0E0*I_ERI_Fxyz_Dyz_Px_S_ab;
  abcd[555] = 4.0E0*I_ERI_F3y_Dyz_Px_S_ab-2.0E0*2*I_ERI_Py_Dyz_Px_S_b;
  abcd[556] = 4.0E0*I_ERI_F2yz_Dyz_Px_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Px_S_b;
  abcd[557] = 4.0E0*I_ERI_Fy2z_Dyz_Px_S_ab;
  abcd[558] = 4.0E0*I_ERI_F2xy_Dxy_Py_S_ab;
  abcd[559] = 4.0E0*I_ERI_Fx2y_Dxy_Py_S_ab-2.0E0*1*I_ERI_Px_Dxy_Py_S_b;
  abcd[560] = 4.0E0*I_ERI_Fxyz_Dxy_Py_S_ab;
  abcd[561] = 4.0E0*I_ERI_F3y_Dxy_Py_S_ab-2.0E0*2*I_ERI_Py_Dxy_Py_S_b;
  abcd[562] = 4.0E0*I_ERI_F2yz_Dxy_Py_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Py_S_b;
  abcd[563] = 4.0E0*I_ERI_Fy2z_Dxy_Py_S_ab;
  abcd[564] = 4.0E0*I_ERI_F2xy_D2y_Py_S_ab-2.0E0*1*I_ERI_F2xy_S_Py_S_a;
  abcd[565] = 4.0E0*I_ERI_Fx2y_D2y_Py_S_ab-2.0E0*1*I_ERI_Fx2y_S_Py_S_a-2.0E0*1*I_ERI_Px_D2y_Py_S_b+1*I_ERI_Px_S_Py_S;
  abcd[566] = 4.0E0*I_ERI_Fxyz_D2y_Py_S_ab-2.0E0*1*I_ERI_Fxyz_S_Py_S_a;
  abcd[567] = 4.0E0*I_ERI_F3y_D2y_Py_S_ab-2.0E0*1*I_ERI_F3y_S_Py_S_a-2.0E0*2*I_ERI_Py_D2y_Py_S_b+2*1*I_ERI_Py_S_Py_S;
  abcd[568] = 4.0E0*I_ERI_F2yz_D2y_Py_S_ab-2.0E0*1*I_ERI_F2yz_S_Py_S_a-2.0E0*1*I_ERI_Pz_D2y_Py_S_b+1*I_ERI_Pz_S_Py_S;
  abcd[569] = 4.0E0*I_ERI_Fy2z_D2y_Py_S_ab-2.0E0*1*I_ERI_Fy2z_S_Py_S_a;
  abcd[570] = 4.0E0*I_ERI_F2xy_Dyz_Py_S_ab;
  abcd[571] = 4.0E0*I_ERI_Fx2y_Dyz_Py_S_ab-2.0E0*1*I_ERI_Px_Dyz_Py_S_b;
  abcd[572] = 4.0E0*I_ERI_Fxyz_Dyz_Py_S_ab;
  abcd[573] = 4.0E0*I_ERI_F3y_Dyz_Py_S_ab-2.0E0*2*I_ERI_Py_Dyz_Py_S_b;
  abcd[574] = 4.0E0*I_ERI_F2yz_Dyz_Py_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Py_S_b;
  abcd[575] = 4.0E0*I_ERI_Fy2z_Dyz_Py_S_ab;
  abcd[576] = 4.0E0*I_ERI_F2xy_Dxy_Pz_S_ab;
  abcd[577] = 4.0E0*I_ERI_Fx2y_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Px_Dxy_Pz_S_b;
  abcd[578] = 4.0E0*I_ERI_Fxyz_Dxy_Pz_S_ab;
  abcd[579] = 4.0E0*I_ERI_F3y_Dxy_Pz_S_ab-2.0E0*2*I_ERI_Py_Dxy_Pz_S_b;
  abcd[580] = 4.0E0*I_ERI_F2yz_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dxy_Pz_S_b;
  abcd[581] = 4.0E0*I_ERI_Fy2z_Dxy_Pz_S_ab;
  abcd[582] = 4.0E0*I_ERI_F2xy_D2y_Pz_S_ab-2.0E0*1*I_ERI_F2xy_S_Pz_S_a;
  abcd[583] = 4.0E0*I_ERI_Fx2y_D2y_Pz_S_ab-2.0E0*1*I_ERI_Fx2y_S_Pz_S_a-2.0E0*1*I_ERI_Px_D2y_Pz_S_b+1*I_ERI_Px_S_Pz_S;
  abcd[584] = 4.0E0*I_ERI_Fxyz_D2y_Pz_S_ab-2.0E0*1*I_ERI_Fxyz_S_Pz_S_a;
  abcd[585] = 4.0E0*I_ERI_F3y_D2y_Pz_S_ab-2.0E0*1*I_ERI_F3y_S_Pz_S_a-2.0E0*2*I_ERI_Py_D2y_Pz_S_b+2*1*I_ERI_Py_S_Pz_S;
  abcd[586] = 4.0E0*I_ERI_F2yz_D2y_Pz_S_ab-2.0E0*1*I_ERI_F2yz_S_Pz_S_a-2.0E0*1*I_ERI_Pz_D2y_Pz_S_b+1*I_ERI_Pz_S_Pz_S;
  abcd[587] = 4.0E0*I_ERI_Fy2z_D2y_Pz_S_ab-2.0E0*1*I_ERI_Fy2z_S_Pz_S_a;
  abcd[588] = 4.0E0*I_ERI_F2xy_Dyz_Pz_S_ab;
  abcd[589] = 4.0E0*I_ERI_Fx2y_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Px_Dyz_Pz_S_b;
  abcd[590] = 4.0E0*I_ERI_Fxyz_Dyz_Pz_S_ab;
  abcd[591] = 4.0E0*I_ERI_F3y_Dyz_Pz_S_ab-2.0E0*2*I_ERI_Py_Dyz_Pz_S_b;
  abcd[592] = 4.0E0*I_ERI_F2yz_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Pz_S_b;
  abcd[593] = 4.0E0*I_ERI_Fy2z_Dyz_Pz_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_D_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[594] = 4.0E0*I_ERI_F2xy_Dxz_Px_S_ab;
  abcd[595] = 4.0E0*I_ERI_Fx2y_Dxz_Px_S_ab-2.0E0*1*I_ERI_Px_Dxz_Px_S_b;
  abcd[596] = 4.0E0*I_ERI_Fxyz_Dxz_Px_S_ab;
  abcd[597] = 4.0E0*I_ERI_F3y_Dxz_Px_S_ab-2.0E0*2*I_ERI_Py_Dxz_Px_S_b;
  abcd[598] = 4.0E0*I_ERI_F2yz_Dxz_Px_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Px_S_b;
  abcd[599] = 4.0E0*I_ERI_Fy2z_Dxz_Px_S_ab;
  abcd[600] = 4.0E0*I_ERI_F2xy_Dyz_Px_S_ab;
  abcd[601] = 4.0E0*I_ERI_Fx2y_Dyz_Px_S_ab-2.0E0*1*I_ERI_Px_Dyz_Px_S_b;
  abcd[602] = 4.0E0*I_ERI_Fxyz_Dyz_Px_S_ab;
  abcd[603] = 4.0E0*I_ERI_F3y_Dyz_Px_S_ab-2.0E0*2*I_ERI_Py_Dyz_Px_S_b;
  abcd[604] = 4.0E0*I_ERI_F2yz_Dyz_Px_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Px_S_b;
  abcd[605] = 4.0E0*I_ERI_Fy2z_Dyz_Px_S_ab;
  abcd[606] = 4.0E0*I_ERI_F2xy_D2z_Px_S_ab-2.0E0*1*I_ERI_F2xy_S_Px_S_a;
  abcd[607] = 4.0E0*I_ERI_Fx2y_D2z_Px_S_ab-2.0E0*1*I_ERI_Fx2y_S_Px_S_a-2.0E0*1*I_ERI_Px_D2z_Px_S_b+1*I_ERI_Px_S_Px_S;
  abcd[608] = 4.0E0*I_ERI_Fxyz_D2z_Px_S_ab-2.0E0*1*I_ERI_Fxyz_S_Px_S_a;
  abcd[609] = 4.0E0*I_ERI_F3y_D2z_Px_S_ab-2.0E0*1*I_ERI_F3y_S_Px_S_a-2.0E0*2*I_ERI_Py_D2z_Px_S_b+2*1*I_ERI_Py_S_Px_S;
  abcd[610] = 4.0E0*I_ERI_F2yz_D2z_Px_S_ab-2.0E0*1*I_ERI_F2yz_S_Px_S_a-2.0E0*1*I_ERI_Pz_D2z_Px_S_b+1*I_ERI_Pz_S_Px_S;
  abcd[611] = 4.0E0*I_ERI_Fy2z_D2z_Px_S_ab-2.0E0*1*I_ERI_Fy2z_S_Px_S_a;
  abcd[612] = 4.0E0*I_ERI_F2xy_Dxz_Py_S_ab;
  abcd[613] = 4.0E0*I_ERI_Fx2y_Dxz_Py_S_ab-2.0E0*1*I_ERI_Px_Dxz_Py_S_b;
  abcd[614] = 4.0E0*I_ERI_Fxyz_Dxz_Py_S_ab;
  abcd[615] = 4.0E0*I_ERI_F3y_Dxz_Py_S_ab-2.0E0*2*I_ERI_Py_Dxz_Py_S_b;
  abcd[616] = 4.0E0*I_ERI_F2yz_Dxz_Py_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Py_S_b;
  abcd[617] = 4.0E0*I_ERI_Fy2z_Dxz_Py_S_ab;
  abcd[618] = 4.0E0*I_ERI_F2xy_Dyz_Py_S_ab;
  abcd[619] = 4.0E0*I_ERI_Fx2y_Dyz_Py_S_ab-2.0E0*1*I_ERI_Px_Dyz_Py_S_b;
  abcd[620] = 4.0E0*I_ERI_Fxyz_Dyz_Py_S_ab;
  abcd[621] = 4.0E0*I_ERI_F3y_Dyz_Py_S_ab-2.0E0*2*I_ERI_Py_Dyz_Py_S_b;
  abcd[622] = 4.0E0*I_ERI_F2yz_Dyz_Py_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Py_S_b;
  abcd[623] = 4.0E0*I_ERI_Fy2z_Dyz_Py_S_ab;
  abcd[624] = 4.0E0*I_ERI_F2xy_D2z_Py_S_ab-2.0E0*1*I_ERI_F2xy_S_Py_S_a;
  abcd[625] = 4.0E0*I_ERI_Fx2y_D2z_Py_S_ab-2.0E0*1*I_ERI_Fx2y_S_Py_S_a-2.0E0*1*I_ERI_Px_D2z_Py_S_b+1*I_ERI_Px_S_Py_S;
  abcd[626] = 4.0E0*I_ERI_Fxyz_D2z_Py_S_ab-2.0E0*1*I_ERI_Fxyz_S_Py_S_a;
  abcd[627] = 4.0E0*I_ERI_F3y_D2z_Py_S_ab-2.0E0*1*I_ERI_F3y_S_Py_S_a-2.0E0*2*I_ERI_Py_D2z_Py_S_b+2*1*I_ERI_Py_S_Py_S;
  abcd[628] = 4.0E0*I_ERI_F2yz_D2z_Py_S_ab-2.0E0*1*I_ERI_F2yz_S_Py_S_a-2.0E0*1*I_ERI_Pz_D2z_Py_S_b+1*I_ERI_Pz_S_Py_S;
  abcd[629] = 4.0E0*I_ERI_Fy2z_D2z_Py_S_ab-2.0E0*1*I_ERI_Fy2z_S_Py_S_a;
  abcd[630] = 4.0E0*I_ERI_F2xy_Dxz_Pz_S_ab;
  abcd[631] = 4.0E0*I_ERI_Fx2y_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Px_Dxz_Pz_S_b;
  abcd[632] = 4.0E0*I_ERI_Fxyz_Dxz_Pz_S_ab;
  abcd[633] = 4.0E0*I_ERI_F3y_Dxz_Pz_S_ab-2.0E0*2*I_ERI_Py_Dxz_Pz_S_b;
  abcd[634] = 4.0E0*I_ERI_F2yz_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dxz_Pz_S_b;
  abcd[635] = 4.0E0*I_ERI_Fy2z_Dxz_Pz_S_ab;
  abcd[636] = 4.0E0*I_ERI_F2xy_Dyz_Pz_S_ab;
  abcd[637] = 4.0E0*I_ERI_Fx2y_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Px_Dyz_Pz_S_b;
  abcd[638] = 4.0E0*I_ERI_Fxyz_Dyz_Pz_S_ab;
  abcd[639] = 4.0E0*I_ERI_F3y_Dyz_Pz_S_ab-2.0E0*2*I_ERI_Py_Dyz_Pz_S_b;
  abcd[640] = 4.0E0*I_ERI_F2yz_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Pz_Dyz_Pz_S_b;
  abcd[641] = 4.0E0*I_ERI_Fy2z_Dyz_Pz_S_ab;
  abcd[642] = 4.0E0*I_ERI_F2xy_D2z_Pz_S_ab-2.0E0*1*I_ERI_F2xy_S_Pz_S_a;
  abcd[643] = 4.0E0*I_ERI_Fx2y_D2z_Pz_S_ab-2.0E0*1*I_ERI_Fx2y_S_Pz_S_a-2.0E0*1*I_ERI_Px_D2z_Pz_S_b+1*I_ERI_Px_S_Pz_S;
  abcd[644] = 4.0E0*I_ERI_Fxyz_D2z_Pz_S_ab-2.0E0*1*I_ERI_Fxyz_S_Pz_S_a;
  abcd[645] = 4.0E0*I_ERI_F3y_D2z_Pz_S_ab-2.0E0*1*I_ERI_F3y_S_Pz_S_a-2.0E0*2*I_ERI_Py_D2z_Pz_S_b+2*1*I_ERI_Py_S_Pz_S;
  abcd[646] = 4.0E0*I_ERI_F2yz_D2z_Pz_S_ab-2.0E0*1*I_ERI_F2yz_S_Pz_S_a-2.0E0*1*I_ERI_Pz_D2z_Pz_S_b+1*I_ERI_Pz_S_Pz_S;
  abcd[647] = 4.0E0*I_ERI_Fy2z_D2z_Pz_S_ab-2.0E0*1*I_ERI_Fy2z_S_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_D_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[648] = 4.0E0*I_ERI_F2xz_D2x_Px_S_ab-2.0E0*1*I_ERI_F2xz_S_Px_S_a;
  abcd[649] = 4.0E0*I_ERI_Fxyz_D2x_Px_S_ab-2.0E0*1*I_ERI_Fxyz_S_Px_S_a;
  abcd[650] = 4.0E0*I_ERI_Fx2z_D2x_Px_S_ab-2.0E0*1*I_ERI_Fx2z_S_Px_S_a-2.0E0*1*I_ERI_Px_D2x_Px_S_b+1*I_ERI_Px_S_Px_S;
  abcd[651] = 4.0E0*I_ERI_F2yz_D2x_Px_S_ab-2.0E0*1*I_ERI_F2yz_S_Px_S_a;
  abcd[652] = 4.0E0*I_ERI_Fy2z_D2x_Px_S_ab-2.0E0*1*I_ERI_Fy2z_S_Px_S_a-2.0E0*1*I_ERI_Py_D2x_Px_S_b+1*I_ERI_Py_S_Px_S;
  abcd[653] = 4.0E0*I_ERI_F3z_D2x_Px_S_ab-2.0E0*1*I_ERI_F3z_S_Px_S_a-2.0E0*2*I_ERI_Pz_D2x_Px_S_b+2*1*I_ERI_Pz_S_Px_S;
  abcd[654] = 4.0E0*I_ERI_F2xz_Dxy_Px_S_ab;
  abcd[655] = 4.0E0*I_ERI_Fxyz_Dxy_Px_S_ab;
  abcd[656] = 4.0E0*I_ERI_Fx2z_Dxy_Px_S_ab-2.0E0*1*I_ERI_Px_Dxy_Px_S_b;
  abcd[657] = 4.0E0*I_ERI_F2yz_Dxy_Px_S_ab;
  abcd[658] = 4.0E0*I_ERI_Fy2z_Dxy_Px_S_ab-2.0E0*1*I_ERI_Py_Dxy_Px_S_b;
  abcd[659] = 4.0E0*I_ERI_F3z_Dxy_Px_S_ab-2.0E0*2*I_ERI_Pz_Dxy_Px_S_b;
  abcd[660] = 4.0E0*I_ERI_F2xz_Dxz_Px_S_ab;
  abcd[661] = 4.0E0*I_ERI_Fxyz_Dxz_Px_S_ab;
  abcd[662] = 4.0E0*I_ERI_Fx2z_Dxz_Px_S_ab-2.0E0*1*I_ERI_Px_Dxz_Px_S_b;
  abcd[663] = 4.0E0*I_ERI_F2yz_Dxz_Px_S_ab;
  abcd[664] = 4.0E0*I_ERI_Fy2z_Dxz_Px_S_ab-2.0E0*1*I_ERI_Py_Dxz_Px_S_b;
  abcd[665] = 4.0E0*I_ERI_F3z_Dxz_Px_S_ab-2.0E0*2*I_ERI_Pz_Dxz_Px_S_b;
  abcd[666] = 4.0E0*I_ERI_F2xz_D2x_Py_S_ab-2.0E0*1*I_ERI_F2xz_S_Py_S_a;
  abcd[667] = 4.0E0*I_ERI_Fxyz_D2x_Py_S_ab-2.0E0*1*I_ERI_Fxyz_S_Py_S_a;
  abcd[668] = 4.0E0*I_ERI_Fx2z_D2x_Py_S_ab-2.0E0*1*I_ERI_Fx2z_S_Py_S_a-2.0E0*1*I_ERI_Px_D2x_Py_S_b+1*I_ERI_Px_S_Py_S;
  abcd[669] = 4.0E0*I_ERI_F2yz_D2x_Py_S_ab-2.0E0*1*I_ERI_F2yz_S_Py_S_a;
  abcd[670] = 4.0E0*I_ERI_Fy2z_D2x_Py_S_ab-2.0E0*1*I_ERI_Fy2z_S_Py_S_a-2.0E0*1*I_ERI_Py_D2x_Py_S_b+1*I_ERI_Py_S_Py_S;
  abcd[671] = 4.0E0*I_ERI_F3z_D2x_Py_S_ab-2.0E0*1*I_ERI_F3z_S_Py_S_a-2.0E0*2*I_ERI_Pz_D2x_Py_S_b+2*1*I_ERI_Pz_S_Py_S;
  abcd[672] = 4.0E0*I_ERI_F2xz_Dxy_Py_S_ab;
  abcd[673] = 4.0E0*I_ERI_Fxyz_Dxy_Py_S_ab;
  abcd[674] = 4.0E0*I_ERI_Fx2z_Dxy_Py_S_ab-2.0E0*1*I_ERI_Px_Dxy_Py_S_b;
  abcd[675] = 4.0E0*I_ERI_F2yz_Dxy_Py_S_ab;
  abcd[676] = 4.0E0*I_ERI_Fy2z_Dxy_Py_S_ab-2.0E0*1*I_ERI_Py_Dxy_Py_S_b;
  abcd[677] = 4.0E0*I_ERI_F3z_Dxy_Py_S_ab-2.0E0*2*I_ERI_Pz_Dxy_Py_S_b;
  abcd[678] = 4.0E0*I_ERI_F2xz_Dxz_Py_S_ab;
  abcd[679] = 4.0E0*I_ERI_Fxyz_Dxz_Py_S_ab;
  abcd[680] = 4.0E0*I_ERI_Fx2z_Dxz_Py_S_ab-2.0E0*1*I_ERI_Px_Dxz_Py_S_b;
  abcd[681] = 4.0E0*I_ERI_F2yz_Dxz_Py_S_ab;
  abcd[682] = 4.0E0*I_ERI_Fy2z_Dxz_Py_S_ab-2.0E0*1*I_ERI_Py_Dxz_Py_S_b;
  abcd[683] = 4.0E0*I_ERI_F3z_Dxz_Py_S_ab-2.0E0*2*I_ERI_Pz_Dxz_Py_S_b;
  abcd[684] = 4.0E0*I_ERI_F2xz_D2x_Pz_S_ab-2.0E0*1*I_ERI_F2xz_S_Pz_S_a;
  abcd[685] = 4.0E0*I_ERI_Fxyz_D2x_Pz_S_ab-2.0E0*1*I_ERI_Fxyz_S_Pz_S_a;
  abcd[686] = 4.0E0*I_ERI_Fx2z_D2x_Pz_S_ab-2.0E0*1*I_ERI_Fx2z_S_Pz_S_a-2.0E0*1*I_ERI_Px_D2x_Pz_S_b+1*I_ERI_Px_S_Pz_S;
  abcd[687] = 4.0E0*I_ERI_F2yz_D2x_Pz_S_ab-2.0E0*1*I_ERI_F2yz_S_Pz_S_a;
  abcd[688] = 4.0E0*I_ERI_Fy2z_D2x_Pz_S_ab-2.0E0*1*I_ERI_Fy2z_S_Pz_S_a-2.0E0*1*I_ERI_Py_D2x_Pz_S_b+1*I_ERI_Py_S_Pz_S;
  abcd[689] = 4.0E0*I_ERI_F3z_D2x_Pz_S_ab-2.0E0*1*I_ERI_F3z_S_Pz_S_a-2.0E0*2*I_ERI_Pz_D2x_Pz_S_b+2*1*I_ERI_Pz_S_Pz_S;
  abcd[690] = 4.0E0*I_ERI_F2xz_Dxy_Pz_S_ab;
  abcd[691] = 4.0E0*I_ERI_Fxyz_Dxy_Pz_S_ab;
  abcd[692] = 4.0E0*I_ERI_Fx2z_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Px_Dxy_Pz_S_b;
  abcd[693] = 4.0E0*I_ERI_F2yz_Dxy_Pz_S_ab;
  abcd[694] = 4.0E0*I_ERI_Fy2z_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Py_Dxy_Pz_S_b;
  abcd[695] = 4.0E0*I_ERI_F3z_Dxy_Pz_S_ab-2.0E0*2*I_ERI_Pz_Dxy_Pz_S_b;
  abcd[696] = 4.0E0*I_ERI_F2xz_Dxz_Pz_S_ab;
  abcd[697] = 4.0E0*I_ERI_Fxyz_Dxz_Pz_S_ab;
  abcd[698] = 4.0E0*I_ERI_Fx2z_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Px_Dxz_Pz_S_b;
  abcd[699] = 4.0E0*I_ERI_F2yz_Dxz_Pz_S_ab;
  abcd[700] = 4.0E0*I_ERI_Fy2z_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Py_Dxz_Pz_S_b;
  abcd[701] = 4.0E0*I_ERI_F3z_Dxz_Pz_S_ab-2.0E0*2*I_ERI_Pz_Dxz_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_D_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[702] = 4.0E0*I_ERI_F2xz_Dxy_Px_S_ab;
  abcd[703] = 4.0E0*I_ERI_Fxyz_Dxy_Px_S_ab;
  abcd[704] = 4.0E0*I_ERI_Fx2z_Dxy_Px_S_ab-2.0E0*1*I_ERI_Px_Dxy_Px_S_b;
  abcd[705] = 4.0E0*I_ERI_F2yz_Dxy_Px_S_ab;
  abcd[706] = 4.0E0*I_ERI_Fy2z_Dxy_Px_S_ab-2.0E0*1*I_ERI_Py_Dxy_Px_S_b;
  abcd[707] = 4.0E0*I_ERI_F3z_Dxy_Px_S_ab-2.0E0*2*I_ERI_Pz_Dxy_Px_S_b;
  abcd[708] = 4.0E0*I_ERI_F2xz_D2y_Px_S_ab-2.0E0*1*I_ERI_F2xz_S_Px_S_a;
  abcd[709] = 4.0E0*I_ERI_Fxyz_D2y_Px_S_ab-2.0E0*1*I_ERI_Fxyz_S_Px_S_a;
  abcd[710] = 4.0E0*I_ERI_Fx2z_D2y_Px_S_ab-2.0E0*1*I_ERI_Fx2z_S_Px_S_a-2.0E0*1*I_ERI_Px_D2y_Px_S_b+1*I_ERI_Px_S_Px_S;
  abcd[711] = 4.0E0*I_ERI_F2yz_D2y_Px_S_ab-2.0E0*1*I_ERI_F2yz_S_Px_S_a;
  abcd[712] = 4.0E0*I_ERI_Fy2z_D2y_Px_S_ab-2.0E0*1*I_ERI_Fy2z_S_Px_S_a-2.0E0*1*I_ERI_Py_D2y_Px_S_b+1*I_ERI_Py_S_Px_S;
  abcd[713] = 4.0E0*I_ERI_F3z_D2y_Px_S_ab-2.0E0*1*I_ERI_F3z_S_Px_S_a-2.0E0*2*I_ERI_Pz_D2y_Px_S_b+2*1*I_ERI_Pz_S_Px_S;
  abcd[714] = 4.0E0*I_ERI_F2xz_Dyz_Px_S_ab;
  abcd[715] = 4.0E0*I_ERI_Fxyz_Dyz_Px_S_ab;
  abcd[716] = 4.0E0*I_ERI_Fx2z_Dyz_Px_S_ab-2.0E0*1*I_ERI_Px_Dyz_Px_S_b;
  abcd[717] = 4.0E0*I_ERI_F2yz_Dyz_Px_S_ab;
  abcd[718] = 4.0E0*I_ERI_Fy2z_Dyz_Px_S_ab-2.0E0*1*I_ERI_Py_Dyz_Px_S_b;
  abcd[719] = 4.0E0*I_ERI_F3z_Dyz_Px_S_ab-2.0E0*2*I_ERI_Pz_Dyz_Px_S_b;
  abcd[720] = 4.0E0*I_ERI_F2xz_Dxy_Py_S_ab;
  abcd[721] = 4.0E0*I_ERI_Fxyz_Dxy_Py_S_ab;
  abcd[722] = 4.0E0*I_ERI_Fx2z_Dxy_Py_S_ab-2.0E0*1*I_ERI_Px_Dxy_Py_S_b;
  abcd[723] = 4.0E0*I_ERI_F2yz_Dxy_Py_S_ab;
  abcd[724] = 4.0E0*I_ERI_Fy2z_Dxy_Py_S_ab-2.0E0*1*I_ERI_Py_Dxy_Py_S_b;
  abcd[725] = 4.0E0*I_ERI_F3z_Dxy_Py_S_ab-2.0E0*2*I_ERI_Pz_Dxy_Py_S_b;
  abcd[726] = 4.0E0*I_ERI_F2xz_D2y_Py_S_ab-2.0E0*1*I_ERI_F2xz_S_Py_S_a;
  abcd[727] = 4.0E0*I_ERI_Fxyz_D2y_Py_S_ab-2.0E0*1*I_ERI_Fxyz_S_Py_S_a;
  abcd[728] = 4.0E0*I_ERI_Fx2z_D2y_Py_S_ab-2.0E0*1*I_ERI_Fx2z_S_Py_S_a-2.0E0*1*I_ERI_Px_D2y_Py_S_b+1*I_ERI_Px_S_Py_S;
  abcd[729] = 4.0E0*I_ERI_F2yz_D2y_Py_S_ab-2.0E0*1*I_ERI_F2yz_S_Py_S_a;
  abcd[730] = 4.0E0*I_ERI_Fy2z_D2y_Py_S_ab-2.0E0*1*I_ERI_Fy2z_S_Py_S_a-2.0E0*1*I_ERI_Py_D2y_Py_S_b+1*I_ERI_Py_S_Py_S;
  abcd[731] = 4.0E0*I_ERI_F3z_D2y_Py_S_ab-2.0E0*1*I_ERI_F3z_S_Py_S_a-2.0E0*2*I_ERI_Pz_D2y_Py_S_b+2*1*I_ERI_Pz_S_Py_S;
  abcd[732] = 4.0E0*I_ERI_F2xz_Dyz_Py_S_ab;
  abcd[733] = 4.0E0*I_ERI_Fxyz_Dyz_Py_S_ab;
  abcd[734] = 4.0E0*I_ERI_Fx2z_Dyz_Py_S_ab-2.0E0*1*I_ERI_Px_Dyz_Py_S_b;
  abcd[735] = 4.0E0*I_ERI_F2yz_Dyz_Py_S_ab;
  abcd[736] = 4.0E0*I_ERI_Fy2z_Dyz_Py_S_ab-2.0E0*1*I_ERI_Py_Dyz_Py_S_b;
  abcd[737] = 4.0E0*I_ERI_F3z_Dyz_Py_S_ab-2.0E0*2*I_ERI_Pz_Dyz_Py_S_b;
  abcd[738] = 4.0E0*I_ERI_F2xz_Dxy_Pz_S_ab;
  abcd[739] = 4.0E0*I_ERI_Fxyz_Dxy_Pz_S_ab;
  abcd[740] = 4.0E0*I_ERI_Fx2z_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Px_Dxy_Pz_S_b;
  abcd[741] = 4.0E0*I_ERI_F2yz_Dxy_Pz_S_ab;
  abcd[742] = 4.0E0*I_ERI_Fy2z_Dxy_Pz_S_ab-2.0E0*1*I_ERI_Py_Dxy_Pz_S_b;
  abcd[743] = 4.0E0*I_ERI_F3z_Dxy_Pz_S_ab-2.0E0*2*I_ERI_Pz_Dxy_Pz_S_b;
  abcd[744] = 4.0E0*I_ERI_F2xz_D2y_Pz_S_ab-2.0E0*1*I_ERI_F2xz_S_Pz_S_a;
  abcd[745] = 4.0E0*I_ERI_Fxyz_D2y_Pz_S_ab-2.0E0*1*I_ERI_Fxyz_S_Pz_S_a;
  abcd[746] = 4.0E0*I_ERI_Fx2z_D2y_Pz_S_ab-2.0E0*1*I_ERI_Fx2z_S_Pz_S_a-2.0E0*1*I_ERI_Px_D2y_Pz_S_b+1*I_ERI_Px_S_Pz_S;
  abcd[747] = 4.0E0*I_ERI_F2yz_D2y_Pz_S_ab-2.0E0*1*I_ERI_F2yz_S_Pz_S_a;
  abcd[748] = 4.0E0*I_ERI_Fy2z_D2y_Pz_S_ab-2.0E0*1*I_ERI_Fy2z_S_Pz_S_a-2.0E0*1*I_ERI_Py_D2y_Pz_S_b+1*I_ERI_Py_S_Pz_S;
  abcd[749] = 4.0E0*I_ERI_F3z_D2y_Pz_S_ab-2.0E0*1*I_ERI_F3z_S_Pz_S_a-2.0E0*2*I_ERI_Pz_D2y_Pz_S_b+2*1*I_ERI_Pz_S_Pz_S;
  abcd[750] = 4.0E0*I_ERI_F2xz_Dyz_Pz_S_ab;
  abcd[751] = 4.0E0*I_ERI_Fxyz_Dyz_Pz_S_ab;
  abcd[752] = 4.0E0*I_ERI_Fx2z_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Px_Dyz_Pz_S_b;
  abcd[753] = 4.0E0*I_ERI_F2yz_Dyz_Pz_S_ab;
  abcd[754] = 4.0E0*I_ERI_Fy2z_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Py_Dyz_Pz_S_b;
  abcd[755] = 4.0E0*I_ERI_F3z_Dyz_Pz_S_ab-2.0E0*2*I_ERI_Pz_Dyz_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_P_S_a
   * RHS shell quartet name: SQ_ERI_P_D_P_S_b
   * RHS shell quartet name: SQ_ERI_P_S_P_S
   ************************************************************/
  abcd[756] = 4.0E0*I_ERI_F2xz_Dxz_Px_S_ab;
  abcd[757] = 4.0E0*I_ERI_Fxyz_Dxz_Px_S_ab;
  abcd[758] = 4.0E0*I_ERI_Fx2z_Dxz_Px_S_ab-2.0E0*1*I_ERI_Px_Dxz_Px_S_b;
  abcd[759] = 4.0E0*I_ERI_F2yz_Dxz_Px_S_ab;
  abcd[760] = 4.0E0*I_ERI_Fy2z_Dxz_Px_S_ab-2.0E0*1*I_ERI_Py_Dxz_Px_S_b;
  abcd[761] = 4.0E0*I_ERI_F3z_Dxz_Px_S_ab-2.0E0*2*I_ERI_Pz_Dxz_Px_S_b;
  abcd[762] = 4.0E0*I_ERI_F2xz_Dyz_Px_S_ab;
  abcd[763] = 4.0E0*I_ERI_Fxyz_Dyz_Px_S_ab;
  abcd[764] = 4.0E0*I_ERI_Fx2z_Dyz_Px_S_ab-2.0E0*1*I_ERI_Px_Dyz_Px_S_b;
  abcd[765] = 4.0E0*I_ERI_F2yz_Dyz_Px_S_ab;
  abcd[766] = 4.0E0*I_ERI_Fy2z_Dyz_Px_S_ab-2.0E0*1*I_ERI_Py_Dyz_Px_S_b;
  abcd[767] = 4.0E0*I_ERI_F3z_Dyz_Px_S_ab-2.0E0*2*I_ERI_Pz_Dyz_Px_S_b;
  abcd[768] = 4.0E0*I_ERI_F2xz_D2z_Px_S_ab-2.0E0*1*I_ERI_F2xz_S_Px_S_a;
  abcd[769] = 4.0E0*I_ERI_Fxyz_D2z_Px_S_ab-2.0E0*1*I_ERI_Fxyz_S_Px_S_a;
  abcd[770] = 4.0E0*I_ERI_Fx2z_D2z_Px_S_ab-2.0E0*1*I_ERI_Fx2z_S_Px_S_a-2.0E0*1*I_ERI_Px_D2z_Px_S_b+1*I_ERI_Px_S_Px_S;
  abcd[771] = 4.0E0*I_ERI_F2yz_D2z_Px_S_ab-2.0E0*1*I_ERI_F2yz_S_Px_S_a;
  abcd[772] = 4.0E0*I_ERI_Fy2z_D2z_Px_S_ab-2.0E0*1*I_ERI_Fy2z_S_Px_S_a-2.0E0*1*I_ERI_Py_D2z_Px_S_b+1*I_ERI_Py_S_Px_S;
  abcd[773] = 4.0E0*I_ERI_F3z_D2z_Px_S_ab-2.0E0*1*I_ERI_F3z_S_Px_S_a-2.0E0*2*I_ERI_Pz_D2z_Px_S_b+2*1*I_ERI_Pz_S_Px_S;
  abcd[774] = 4.0E0*I_ERI_F2xz_Dxz_Py_S_ab;
  abcd[775] = 4.0E0*I_ERI_Fxyz_Dxz_Py_S_ab;
  abcd[776] = 4.0E0*I_ERI_Fx2z_Dxz_Py_S_ab-2.0E0*1*I_ERI_Px_Dxz_Py_S_b;
  abcd[777] = 4.0E0*I_ERI_F2yz_Dxz_Py_S_ab;
  abcd[778] = 4.0E0*I_ERI_Fy2z_Dxz_Py_S_ab-2.0E0*1*I_ERI_Py_Dxz_Py_S_b;
  abcd[779] = 4.0E0*I_ERI_F3z_Dxz_Py_S_ab-2.0E0*2*I_ERI_Pz_Dxz_Py_S_b;
  abcd[780] = 4.0E0*I_ERI_F2xz_Dyz_Py_S_ab;
  abcd[781] = 4.0E0*I_ERI_Fxyz_Dyz_Py_S_ab;
  abcd[782] = 4.0E0*I_ERI_Fx2z_Dyz_Py_S_ab-2.0E0*1*I_ERI_Px_Dyz_Py_S_b;
  abcd[783] = 4.0E0*I_ERI_F2yz_Dyz_Py_S_ab;
  abcd[784] = 4.0E0*I_ERI_Fy2z_Dyz_Py_S_ab-2.0E0*1*I_ERI_Py_Dyz_Py_S_b;
  abcd[785] = 4.0E0*I_ERI_F3z_Dyz_Py_S_ab-2.0E0*2*I_ERI_Pz_Dyz_Py_S_b;
  abcd[786] = 4.0E0*I_ERI_F2xz_D2z_Py_S_ab-2.0E0*1*I_ERI_F2xz_S_Py_S_a;
  abcd[787] = 4.0E0*I_ERI_Fxyz_D2z_Py_S_ab-2.0E0*1*I_ERI_Fxyz_S_Py_S_a;
  abcd[788] = 4.0E0*I_ERI_Fx2z_D2z_Py_S_ab-2.0E0*1*I_ERI_Fx2z_S_Py_S_a-2.0E0*1*I_ERI_Px_D2z_Py_S_b+1*I_ERI_Px_S_Py_S;
  abcd[789] = 4.0E0*I_ERI_F2yz_D2z_Py_S_ab-2.0E0*1*I_ERI_F2yz_S_Py_S_a;
  abcd[790] = 4.0E0*I_ERI_Fy2z_D2z_Py_S_ab-2.0E0*1*I_ERI_Fy2z_S_Py_S_a-2.0E0*1*I_ERI_Py_D2z_Py_S_b+1*I_ERI_Py_S_Py_S;
  abcd[791] = 4.0E0*I_ERI_F3z_D2z_Py_S_ab-2.0E0*1*I_ERI_F3z_S_Py_S_a-2.0E0*2*I_ERI_Pz_D2z_Py_S_b+2*1*I_ERI_Pz_S_Py_S;
  abcd[792] = 4.0E0*I_ERI_F2xz_Dxz_Pz_S_ab;
  abcd[793] = 4.0E0*I_ERI_Fxyz_Dxz_Pz_S_ab;
  abcd[794] = 4.0E0*I_ERI_Fx2z_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Px_Dxz_Pz_S_b;
  abcd[795] = 4.0E0*I_ERI_F2yz_Dxz_Pz_S_ab;
  abcd[796] = 4.0E0*I_ERI_Fy2z_Dxz_Pz_S_ab-2.0E0*1*I_ERI_Py_Dxz_Pz_S_b;
  abcd[797] = 4.0E0*I_ERI_F3z_Dxz_Pz_S_ab-2.0E0*2*I_ERI_Pz_Dxz_Pz_S_b;
  abcd[798] = 4.0E0*I_ERI_F2xz_Dyz_Pz_S_ab;
  abcd[799] = 4.0E0*I_ERI_Fxyz_Dyz_Pz_S_ab;
  abcd[800] = 4.0E0*I_ERI_Fx2z_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Px_Dyz_Pz_S_b;
  abcd[801] = 4.0E0*I_ERI_F2yz_Dyz_Pz_S_ab;
  abcd[802] = 4.0E0*I_ERI_Fy2z_Dyz_Pz_S_ab-2.0E0*1*I_ERI_Py_Dyz_Pz_S_b;
  abcd[803] = 4.0E0*I_ERI_F3z_Dyz_Pz_S_ab-2.0E0*2*I_ERI_Pz_Dyz_Pz_S_b;
  abcd[804] = 4.0E0*I_ERI_F2xz_D2z_Pz_S_ab-2.0E0*1*I_ERI_F2xz_S_Pz_S_a;
  abcd[805] = 4.0E0*I_ERI_Fxyz_D2z_Pz_S_ab-2.0E0*1*I_ERI_Fxyz_S_Pz_S_a;
  abcd[806] = 4.0E0*I_ERI_Fx2z_D2z_Pz_S_ab-2.0E0*1*I_ERI_Fx2z_S_Pz_S_a-2.0E0*1*I_ERI_Px_D2z_Pz_S_b+1*I_ERI_Px_S_Pz_S;
  abcd[807] = 4.0E0*I_ERI_F2yz_D2z_Pz_S_ab-2.0E0*1*I_ERI_F2yz_S_Pz_S_a;
  abcd[808] = 4.0E0*I_ERI_Fy2z_D2z_Pz_S_ab-2.0E0*1*I_ERI_Fy2z_S_Pz_S_a-2.0E0*1*I_ERI_Py_D2z_Pz_S_b+1*I_ERI_Py_S_Pz_S;
  abcd[809] = 4.0E0*I_ERI_F3z_D2z_Pz_S_ab-2.0E0*1*I_ERI_F3z_S_Pz_S_a-2.0E0*2*I_ERI_Pz_D2z_Pz_S_b+2*1*I_ERI_Pz_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   * RHS shell quartet name: SQ_ERI_P_P_D_S_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  abcd[810] = 4.0E0*I_ERI_F3x_Px_D2x_S_ac-2.0E0*1*I_ERI_F3x_Px_S_S_a-2.0E0*2*I_ERI_Px_Px_D2x_S_c+2*1*I_ERI_Px_Px_S_S;
  abcd[811] = 4.0E0*I_ERI_F2xy_Px_D2x_S_ac-2.0E0*1*I_ERI_F2xy_Px_S_S_a-2.0E0*1*I_ERI_Py_Px_D2x_S_c+1*I_ERI_Py_Px_S_S;
  abcd[812] = 4.0E0*I_ERI_F2xz_Px_D2x_S_ac-2.0E0*1*I_ERI_F2xz_Px_S_S_a-2.0E0*1*I_ERI_Pz_Px_D2x_S_c+1*I_ERI_Pz_Px_S_S;
  abcd[813] = 4.0E0*I_ERI_Fx2y_Px_D2x_S_ac-2.0E0*1*I_ERI_Fx2y_Px_S_S_a;
  abcd[814] = 4.0E0*I_ERI_Fxyz_Px_D2x_S_ac-2.0E0*1*I_ERI_Fxyz_Px_S_S_a;
  abcd[815] = 4.0E0*I_ERI_Fx2z_Px_D2x_S_ac-2.0E0*1*I_ERI_Fx2z_Px_S_S_a;
  abcd[816] = 4.0E0*I_ERI_F3x_Py_D2x_S_ac-2.0E0*1*I_ERI_F3x_Py_S_S_a-2.0E0*2*I_ERI_Px_Py_D2x_S_c+2*1*I_ERI_Px_Py_S_S;
  abcd[817] = 4.0E0*I_ERI_F2xy_Py_D2x_S_ac-2.0E0*1*I_ERI_F2xy_Py_S_S_a-2.0E0*1*I_ERI_Py_Py_D2x_S_c+1*I_ERI_Py_Py_S_S;
  abcd[818] = 4.0E0*I_ERI_F2xz_Py_D2x_S_ac-2.0E0*1*I_ERI_F2xz_Py_S_S_a-2.0E0*1*I_ERI_Pz_Py_D2x_S_c+1*I_ERI_Pz_Py_S_S;
  abcd[819] = 4.0E0*I_ERI_Fx2y_Py_D2x_S_ac-2.0E0*1*I_ERI_Fx2y_Py_S_S_a;
  abcd[820] = 4.0E0*I_ERI_Fxyz_Py_D2x_S_ac-2.0E0*1*I_ERI_Fxyz_Py_S_S_a;
  abcd[821] = 4.0E0*I_ERI_Fx2z_Py_D2x_S_ac-2.0E0*1*I_ERI_Fx2z_Py_S_S_a;
  abcd[822] = 4.0E0*I_ERI_F3x_Pz_D2x_S_ac-2.0E0*1*I_ERI_F3x_Pz_S_S_a-2.0E0*2*I_ERI_Px_Pz_D2x_S_c+2*1*I_ERI_Px_Pz_S_S;
  abcd[823] = 4.0E0*I_ERI_F2xy_Pz_D2x_S_ac-2.0E0*1*I_ERI_F2xy_Pz_S_S_a-2.0E0*1*I_ERI_Py_Pz_D2x_S_c+1*I_ERI_Py_Pz_S_S;
  abcd[824] = 4.0E0*I_ERI_F2xz_Pz_D2x_S_ac-2.0E0*1*I_ERI_F2xz_Pz_S_S_a-2.0E0*1*I_ERI_Pz_Pz_D2x_S_c+1*I_ERI_Pz_Pz_S_S;
  abcd[825] = 4.0E0*I_ERI_Fx2y_Pz_D2x_S_ac-2.0E0*1*I_ERI_Fx2y_Pz_S_S_a;
  abcd[826] = 4.0E0*I_ERI_Fxyz_Pz_D2x_S_ac-2.0E0*1*I_ERI_Fxyz_Pz_S_S_a;
  abcd[827] = 4.0E0*I_ERI_Fx2z_Pz_D2x_S_ac-2.0E0*1*I_ERI_Fx2z_Pz_S_S_a;
  abcd[828] = 4.0E0*I_ERI_F3x_Px_Dxy_S_ac-2.0E0*2*I_ERI_Px_Px_Dxy_S_c;
  abcd[829] = 4.0E0*I_ERI_F2xy_Px_Dxy_S_ac-2.0E0*1*I_ERI_Py_Px_Dxy_S_c;
  abcd[830] = 4.0E0*I_ERI_F2xz_Px_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Px_Dxy_S_c;
  abcd[831] = 4.0E0*I_ERI_Fx2y_Px_Dxy_S_ac;
  abcd[832] = 4.0E0*I_ERI_Fxyz_Px_Dxy_S_ac;
  abcd[833] = 4.0E0*I_ERI_Fx2z_Px_Dxy_S_ac;
  abcd[834] = 4.0E0*I_ERI_F3x_Py_Dxy_S_ac-2.0E0*2*I_ERI_Px_Py_Dxy_S_c;
  abcd[835] = 4.0E0*I_ERI_F2xy_Py_Dxy_S_ac-2.0E0*1*I_ERI_Py_Py_Dxy_S_c;
  abcd[836] = 4.0E0*I_ERI_F2xz_Py_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Py_Dxy_S_c;
  abcd[837] = 4.0E0*I_ERI_Fx2y_Py_Dxy_S_ac;
  abcd[838] = 4.0E0*I_ERI_Fxyz_Py_Dxy_S_ac;
  abcd[839] = 4.0E0*I_ERI_Fx2z_Py_Dxy_S_ac;
  abcd[840] = 4.0E0*I_ERI_F3x_Pz_Dxy_S_ac-2.0E0*2*I_ERI_Px_Pz_Dxy_S_c;
  abcd[841] = 4.0E0*I_ERI_F2xy_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Py_Pz_Dxy_S_c;
  abcd[842] = 4.0E0*I_ERI_F2xz_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dxy_S_c;
  abcd[843] = 4.0E0*I_ERI_Fx2y_Pz_Dxy_S_ac;
  abcd[844] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_S_ac;
  abcd[845] = 4.0E0*I_ERI_Fx2z_Pz_Dxy_S_ac;
  abcd[846] = 4.0E0*I_ERI_F3x_Px_Dxz_S_ac-2.0E0*2*I_ERI_Px_Px_Dxz_S_c;
  abcd[847] = 4.0E0*I_ERI_F2xy_Px_Dxz_S_ac-2.0E0*1*I_ERI_Py_Px_Dxz_S_c;
  abcd[848] = 4.0E0*I_ERI_F2xz_Px_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Px_Dxz_S_c;
  abcd[849] = 4.0E0*I_ERI_Fx2y_Px_Dxz_S_ac;
  abcd[850] = 4.0E0*I_ERI_Fxyz_Px_Dxz_S_ac;
  abcd[851] = 4.0E0*I_ERI_Fx2z_Px_Dxz_S_ac;
  abcd[852] = 4.0E0*I_ERI_F3x_Py_Dxz_S_ac-2.0E0*2*I_ERI_Px_Py_Dxz_S_c;
  abcd[853] = 4.0E0*I_ERI_F2xy_Py_Dxz_S_ac-2.0E0*1*I_ERI_Py_Py_Dxz_S_c;
  abcd[854] = 4.0E0*I_ERI_F2xz_Py_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Py_Dxz_S_c;
  abcd[855] = 4.0E0*I_ERI_Fx2y_Py_Dxz_S_ac;
  abcd[856] = 4.0E0*I_ERI_Fxyz_Py_Dxz_S_ac;
  abcd[857] = 4.0E0*I_ERI_Fx2z_Py_Dxz_S_ac;
  abcd[858] = 4.0E0*I_ERI_F3x_Pz_Dxz_S_ac-2.0E0*2*I_ERI_Px_Pz_Dxz_S_c;
  abcd[859] = 4.0E0*I_ERI_F2xy_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Py_Pz_Dxz_S_c;
  abcd[860] = 4.0E0*I_ERI_F2xz_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dxz_S_c;
  abcd[861] = 4.0E0*I_ERI_Fx2y_Pz_Dxz_S_ac;
  abcd[862] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_S_ac;
  abcd[863] = 4.0E0*I_ERI_Fx2z_Pz_Dxz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   * RHS shell quartet name: SQ_ERI_P_P_D_S_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  abcd[864] = 4.0E0*I_ERI_F3x_Px_Dxy_S_ac-2.0E0*2*I_ERI_Px_Px_Dxy_S_c;
  abcd[865] = 4.0E0*I_ERI_F2xy_Px_Dxy_S_ac-2.0E0*1*I_ERI_Py_Px_Dxy_S_c;
  abcd[866] = 4.0E0*I_ERI_F2xz_Px_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Px_Dxy_S_c;
  abcd[867] = 4.0E0*I_ERI_Fx2y_Px_Dxy_S_ac;
  abcd[868] = 4.0E0*I_ERI_Fxyz_Px_Dxy_S_ac;
  abcd[869] = 4.0E0*I_ERI_Fx2z_Px_Dxy_S_ac;
  abcd[870] = 4.0E0*I_ERI_F3x_Py_Dxy_S_ac-2.0E0*2*I_ERI_Px_Py_Dxy_S_c;
  abcd[871] = 4.0E0*I_ERI_F2xy_Py_Dxy_S_ac-2.0E0*1*I_ERI_Py_Py_Dxy_S_c;
  abcd[872] = 4.0E0*I_ERI_F2xz_Py_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Py_Dxy_S_c;
  abcd[873] = 4.0E0*I_ERI_Fx2y_Py_Dxy_S_ac;
  abcd[874] = 4.0E0*I_ERI_Fxyz_Py_Dxy_S_ac;
  abcd[875] = 4.0E0*I_ERI_Fx2z_Py_Dxy_S_ac;
  abcd[876] = 4.0E0*I_ERI_F3x_Pz_Dxy_S_ac-2.0E0*2*I_ERI_Px_Pz_Dxy_S_c;
  abcd[877] = 4.0E0*I_ERI_F2xy_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Py_Pz_Dxy_S_c;
  abcd[878] = 4.0E0*I_ERI_F2xz_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dxy_S_c;
  abcd[879] = 4.0E0*I_ERI_Fx2y_Pz_Dxy_S_ac;
  abcd[880] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_S_ac;
  abcd[881] = 4.0E0*I_ERI_Fx2z_Pz_Dxy_S_ac;
  abcd[882] = 4.0E0*I_ERI_F3x_Px_D2y_S_ac-2.0E0*1*I_ERI_F3x_Px_S_S_a-2.0E0*2*I_ERI_Px_Px_D2y_S_c+2*1*I_ERI_Px_Px_S_S;
  abcd[883] = 4.0E0*I_ERI_F2xy_Px_D2y_S_ac-2.0E0*1*I_ERI_F2xy_Px_S_S_a-2.0E0*1*I_ERI_Py_Px_D2y_S_c+1*I_ERI_Py_Px_S_S;
  abcd[884] = 4.0E0*I_ERI_F2xz_Px_D2y_S_ac-2.0E0*1*I_ERI_F2xz_Px_S_S_a-2.0E0*1*I_ERI_Pz_Px_D2y_S_c+1*I_ERI_Pz_Px_S_S;
  abcd[885] = 4.0E0*I_ERI_Fx2y_Px_D2y_S_ac-2.0E0*1*I_ERI_Fx2y_Px_S_S_a;
  abcd[886] = 4.0E0*I_ERI_Fxyz_Px_D2y_S_ac-2.0E0*1*I_ERI_Fxyz_Px_S_S_a;
  abcd[887] = 4.0E0*I_ERI_Fx2z_Px_D2y_S_ac-2.0E0*1*I_ERI_Fx2z_Px_S_S_a;
  abcd[888] = 4.0E0*I_ERI_F3x_Py_D2y_S_ac-2.0E0*1*I_ERI_F3x_Py_S_S_a-2.0E0*2*I_ERI_Px_Py_D2y_S_c+2*1*I_ERI_Px_Py_S_S;
  abcd[889] = 4.0E0*I_ERI_F2xy_Py_D2y_S_ac-2.0E0*1*I_ERI_F2xy_Py_S_S_a-2.0E0*1*I_ERI_Py_Py_D2y_S_c+1*I_ERI_Py_Py_S_S;
  abcd[890] = 4.0E0*I_ERI_F2xz_Py_D2y_S_ac-2.0E0*1*I_ERI_F2xz_Py_S_S_a-2.0E0*1*I_ERI_Pz_Py_D2y_S_c+1*I_ERI_Pz_Py_S_S;
  abcd[891] = 4.0E0*I_ERI_Fx2y_Py_D2y_S_ac-2.0E0*1*I_ERI_Fx2y_Py_S_S_a;
  abcd[892] = 4.0E0*I_ERI_Fxyz_Py_D2y_S_ac-2.0E0*1*I_ERI_Fxyz_Py_S_S_a;
  abcd[893] = 4.0E0*I_ERI_Fx2z_Py_D2y_S_ac-2.0E0*1*I_ERI_Fx2z_Py_S_S_a;
  abcd[894] = 4.0E0*I_ERI_F3x_Pz_D2y_S_ac-2.0E0*1*I_ERI_F3x_Pz_S_S_a-2.0E0*2*I_ERI_Px_Pz_D2y_S_c+2*1*I_ERI_Px_Pz_S_S;
  abcd[895] = 4.0E0*I_ERI_F2xy_Pz_D2y_S_ac-2.0E0*1*I_ERI_F2xy_Pz_S_S_a-2.0E0*1*I_ERI_Py_Pz_D2y_S_c+1*I_ERI_Py_Pz_S_S;
  abcd[896] = 4.0E0*I_ERI_F2xz_Pz_D2y_S_ac-2.0E0*1*I_ERI_F2xz_Pz_S_S_a-2.0E0*1*I_ERI_Pz_Pz_D2y_S_c+1*I_ERI_Pz_Pz_S_S;
  abcd[897] = 4.0E0*I_ERI_Fx2y_Pz_D2y_S_ac-2.0E0*1*I_ERI_Fx2y_Pz_S_S_a;
  abcd[898] = 4.0E0*I_ERI_Fxyz_Pz_D2y_S_ac-2.0E0*1*I_ERI_Fxyz_Pz_S_S_a;
  abcd[899] = 4.0E0*I_ERI_Fx2z_Pz_D2y_S_ac-2.0E0*1*I_ERI_Fx2z_Pz_S_S_a;
  abcd[900] = 4.0E0*I_ERI_F3x_Px_Dyz_S_ac-2.0E0*2*I_ERI_Px_Px_Dyz_S_c;
  abcd[901] = 4.0E0*I_ERI_F2xy_Px_Dyz_S_ac-2.0E0*1*I_ERI_Py_Px_Dyz_S_c;
  abcd[902] = 4.0E0*I_ERI_F2xz_Px_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Px_Dyz_S_c;
  abcd[903] = 4.0E0*I_ERI_Fx2y_Px_Dyz_S_ac;
  abcd[904] = 4.0E0*I_ERI_Fxyz_Px_Dyz_S_ac;
  abcd[905] = 4.0E0*I_ERI_Fx2z_Px_Dyz_S_ac;
  abcd[906] = 4.0E0*I_ERI_F3x_Py_Dyz_S_ac-2.0E0*2*I_ERI_Px_Py_Dyz_S_c;
  abcd[907] = 4.0E0*I_ERI_F2xy_Py_Dyz_S_ac-2.0E0*1*I_ERI_Py_Py_Dyz_S_c;
  abcd[908] = 4.0E0*I_ERI_F2xz_Py_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Py_Dyz_S_c;
  abcd[909] = 4.0E0*I_ERI_Fx2y_Py_Dyz_S_ac;
  abcd[910] = 4.0E0*I_ERI_Fxyz_Py_Dyz_S_ac;
  abcd[911] = 4.0E0*I_ERI_Fx2z_Py_Dyz_S_ac;
  abcd[912] = 4.0E0*I_ERI_F3x_Pz_Dyz_S_ac-2.0E0*2*I_ERI_Px_Pz_Dyz_S_c;
  abcd[913] = 4.0E0*I_ERI_F2xy_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Py_Pz_Dyz_S_c;
  abcd[914] = 4.0E0*I_ERI_F2xz_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dyz_S_c;
  abcd[915] = 4.0E0*I_ERI_Fx2y_Pz_Dyz_S_ac;
  abcd[916] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_S_ac;
  abcd[917] = 4.0E0*I_ERI_Fx2z_Pz_Dyz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   * RHS shell quartet name: SQ_ERI_P_P_D_S_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  abcd[918] = 4.0E0*I_ERI_F3x_Px_Dxz_S_ac-2.0E0*2*I_ERI_Px_Px_Dxz_S_c;
  abcd[919] = 4.0E0*I_ERI_F2xy_Px_Dxz_S_ac-2.0E0*1*I_ERI_Py_Px_Dxz_S_c;
  abcd[920] = 4.0E0*I_ERI_F2xz_Px_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Px_Dxz_S_c;
  abcd[921] = 4.0E0*I_ERI_Fx2y_Px_Dxz_S_ac;
  abcd[922] = 4.0E0*I_ERI_Fxyz_Px_Dxz_S_ac;
  abcd[923] = 4.0E0*I_ERI_Fx2z_Px_Dxz_S_ac;
  abcd[924] = 4.0E0*I_ERI_F3x_Py_Dxz_S_ac-2.0E0*2*I_ERI_Px_Py_Dxz_S_c;
  abcd[925] = 4.0E0*I_ERI_F2xy_Py_Dxz_S_ac-2.0E0*1*I_ERI_Py_Py_Dxz_S_c;
  abcd[926] = 4.0E0*I_ERI_F2xz_Py_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Py_Dxz_S_c;
  abcd[927] = 4.0E0*I_ERI_Fx2y_Py_Dxz_S_ac;
  abcd[928] = 4.0E0*I_ERI_Fxyz_Py_Dxz_S_ac;
  abcd[929] = 4.0E0*I_ERI_Fx2z_Py_Dxz_S_ac;
  abcd[930] = 4.0E0*I_ERI_F3x_Pz_Dxz_S_ac-2.0E0*2*I_ERI_Px_Pz_Dxz_S_c;
  abcd[931] = 4.0E0*I_ERI_F2xy_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Py_Pz_Dxz_S_c;
  abcd[932] = 4.0E0*I_ERI_F2xz_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dxz_S_c;
  abcd[933] = 4.0E0*I_ERI_Fx2y_Pz_Dxz_S_ac;
  abcd[934] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_S_ac;
  abcd[935] = 4.0E0*I_ERI_Fx2z_Pz_Dxz_S_ac;
  abcd[936] = 4.0E0*I_ERI_F3x_Px_Dyz_S_ac-2.0E0*2*I_ERI_Px_Px_Dyz_S_c;
  abcd[937] = 4.0E0*I_ERI_F2xy_Px_Dyz_S_ac-2.0E0*1*I_ERI_Py_Px_Dyz_S_c;
  abcd[938] = 4.0E0*I_ERI_F2xz_Px_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Px_Dyz_S_c;
  abcd[939] = 4.0E0*I_ERI_Fx2y_Px_Dyz_S_ac;
  abcd[940] = 4.0E0*I_ERI_Fxyz_Px_Dyz_S_ac;
  abcd[941] = 4.0E0*I_ERI_Fx2z_Px_Dyz_S_ac;
  abcd[942] = 4.0E0*I_ERI_F3x_Py_Dyz_S_ac-2.0E0*2*I_ERI_Px_Py_Dyz_S_c;
  abcd[943] = 4.0E0*I_ERI_F2xy_Py_Dyz_S_ac-2.0E0*1*I_ERI_Py_Py_Dyz_S_c;
  abcd[944] = 4.0E0*I_ERI_F2xz_Py_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Py_Dyz_S_c;
  abcd[945] = 4.0E0*I_ERI_Fx2y_Py_Dyz_S_ac;
  abcd[946] = 4.0E0*I_ERI_Fxyz_Py_Dyz_S_ac;
  abcd[947] = 4.0E0*I_ERI_Fx2z_Py_Dyz_S_ac;
  abcd[948] = 4.0E0*I_ERI_F3x_Pz_Dyz_S_ac-2.0E0*2*I_ERI_Px_Pz_Dyz_S_c;
  abcd[949] = 4.0E0*I_ERI_F2xy_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Py_Pz_Dyz_S_c;
  abcd[950] = 4.0E0*I_ERI_F2xz_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dyz_S_c;
  abcd[951] = 4.0E0*I_ERI_Fx2y_Pz_Dyz_S_ac;
  abcd[952] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_S_ac;
  abcd[953] = 4.0E0*I_ERI_Fx2z_Pz_Dyz_S_ac;
  abcd[954] = 4.0E0*I_ERI_F3x_Px_D2z_S_ac-2.0E0*1*I_ERI_F3x_Px_S_S_a-2.0E0*2*I_ERI_Px_Px_D2z_S_c+2*1*I_ERI_Px_Px_S_S;
  abcd[955] = 4.0E0*I_ERI_F2xy_Px_D2z_S_ac-2.0E0*1*I_ERI_F2xy_Px_S_S_a-2.0E0*1*I_ERI_Py_Px_D2z_S_c+1*I_ERI_Py_Px_S_S;
  abcd[956] = 4.0E0*I_ERI_F2xz_Px_D2z_S_ac-2.0E0*1*I_ERI_F2xz_Px_S_S_a-2.0E0*1*I_ERI_Pz_Px_D2z_S_c+1*I_ERI_Pz_Px_S_S;
  abcd[957] = 4.0E0*I_ERI_Fx2y_Px_D2z_S_ac-2.0E0*1*I_ERI_Fx2y_Px_S_S_a;
  abcd[958] = 4.0E0*I_ERI_Fxyz_Px_D2z_S_ac-2.0E0*1*I_ERI_Fxyz_Px_S_S_a;
  abcd[959] = 4.0E0*I_ERI_Fx2z_Px_D2z_S_ac-2.0E0*1*I_ERI_Fx2z_Px_S_S_a;
  abcd[960] = 4.0E0*I_ERI_F3x_Py_D2z_S_ac-2.0E0*1*I_ERI_F3x_Py_S_S_a-2.0E0*2*I_ERI_Px_Py_D2z_S_c+2*1*I_ERI_Px_Py_S_S;
  abcd[961] = 4.0E0*I_ERI_F2xy_Py_D2z_S_ac-2.0E0*1*I_ERI_F2xy_Py_S_S_a-2.0E0*1*I_ERI_Py_Py_D2z_S_c+1*I_ERI_Py_Py_S_S;
  abcd[962] = 4.0E0*I_ERI_F2xz_Py_D2z_S_ac-2.0E0*1*I_ERI_F2xz_Py_S_S_a-2.0E0*1*I_ERI_Pz_Py_D2z_S_c+1*I_ERI_Pz_Py_S_S;
  abcd[963] = 4.0E0*I_ERI_Fx2y_Py_D2z_S_ac-2.0E0*1*I_ERI_Fx2y_Py_S_S_a;
  abcd[964] = 4.0E0*I_ERI_Fxyz_Py_D2z_S_ac-2.0E0*1*I_ERI_Fxyz_Py_S_S_a;
  abcd[965] = 4.0E0*I_ERI_Fx2z_Py_D2z_S_ac-2.0E0*1*I_ERI_Fx2z_Py_S_S_a;
  abcd[966] = 4.0E0*I_ERI_F3x_Pz_D2z_S_ac-2.0E0*1*I_ERI_F3x_Pz_S_S_a-2.0E0*2*I_ERI_Px_Pz_D2z_S_c+2*1*I_ERI_Px_Pz_S_S;
  abcd[967] = 4.0E0*I_ERI_F2xy_Pz_D2z_S_ac-2.0E0*1*I_ERI_F2xy_Pz_S_S_a-2.0E0*1*I_ERI_Py_Pz_D2z_S_c+1*I_ERI_Py_Pz_S_S;
  abcd[968] = 4.0E0*I_ERI_F2xz_Pz_D2z_S_ac-2.0E0*1*I_ERI_F2xz_Pz_S_S_a-2.0E0*1*I_ERI_Pz_Pz_D2z_S_c+1*I_ERI_Pz_Pz_S_S;
  abcd[969] = 4.0E0*I_ERI_Fx2y_Pz_D2z_S_ac-2.0E0*1*I_ERI_Fx2y_Pz_S_S_a;
  abcd[970] = 4.0E0*I_ERI_Fxyz_Pz_D2z_S_ac-2.0E0*1*I_ERI_Fxyz_Pz_S_S_a;
  abcd[971] = 4.0E0*I_ERI_Fx2z_Pz_D2z_S_ac-2.0E0*1*I_ERI_Fx2z_Pz_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   * RHS shell quartet name: SQ_ERI_P_P_D_S_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  abcd[972] = 4.0E0*I_ERI_F2xy_Px_D2x_S_ac-2.0E0*1*I_ERI_F2xy_Px_S_S_a;
  abcd[973] = 4.0E0*I_ERI_Fx2y_Px_D2x_S_ac-2.0E0*1*I_ERI_Fx2y_Px_S_S_a-2.0E0*1*I_ERI_Px_Px_D2x_S_c+1*I_ERI_Px_Px_S_S;
  abcd[974] = 4.0E0*I_ERI_Fxyz_Px_D2x_S_ac-2.0E0*1*I_ERI_Fxyz_Px_S_S_a;
  abcd[975] = 4.0E0*I_ERI_F3y_Px_D2x_S_ac-2.0E0*1*I_ERI_F3y_Px_S_S_a-2.0E0*2*I_ERI_Py_Px_D2x_S_c+2*1*I_ERI_Py_Px_S_S;
  abcd[976] = 4.0E0*I_ERI_F2yz_Px_D2x_S_ac-2.0E0*1*I_ERI_F2yz_Px_S_S_a-2.0E0*1*I_ERI_Pz_Px_D2x_S_c+1*I_ERI_Pz_Px_S_S;
  abcd[977] = 4.0E0*I_ERI_Fy2z_Px_D2x_S_ac-2.0E0*1*I_ERI_Fy2z_Px_S_S_a;
  abcd[978] = 4.0E0*I_ERI_F2xy_Py_D2x_S_ac-2.0E0*1*I_ERI_F2xy_Py_S_S_a;
  abcd[979] = 4.0E0*I_ERI_Fx2y_Py_D2x_S_ac-2.0E0*1*I_ERI_Fx2y_Py_S_S_a-2.0E0*1*I_ERI_Px_Py_D2x_S_c+1*I_ERI_Px_Py_S_S;
  abcd[980] = 4.0E0*I_ERI_Fxyz_Py_D2x_S_ac-2.0E0*1*I_ERI_Fxyz_Py_S_S_a;
  abcd[981] = 4.0E0*I_ERI_F3y_Py_D2x_S_ac-2.0E0*1*I_ERI_F3y_Py_S_S_a-2.0E0*2*I_ERI_Py_Py_D2x_S_c+2*1*I_ERI_Py_Py_S_S;
  abcd[982] = 4.0E0*I_ERI_F2yz_Py_D2x_S_ac-2.0E0*1*I_ERI_F2yz_Py_S_S_a-2.0E0*1*I_ERI_Pz_Py_D2x_S_c+1*I_ERI_Pz_Py_S_S;
  abcd[983] = 4.0E0*I_ERI_Fy2z_Py_D2x_S_ac-2.0E0*1*I_ERI_Fy2z_Py_S_S_a;
  abcd[984] = 4.0E0*I_ERI_F2xy_Pz_D2x_S_ac-2.0E0*1*I_ERI_F2xy_Pz_S_S_a;
  abcd[985] = 4.0E0*I_ERI_Fx2y_Pz_D2x_S_ac-2.0E0*1*I_ERI_Fx2y_Pz_S_S_a-2.0E0*1*I_ERI_Px_Pz_D2x_S_c+1*I_ERI_Px_Pz_S_S;
  abcd[986] = 4.0E0*I_ERI_Fxyz_Pz_D2x_S_ac-2.0E0*1*I_ERI_Fxyz_Pz_S_S_a;
  abcd[987] = 4.0E0*I_ERI_F3y_Pz_D2x_S_ac-2.0E0*1*I_ERI_F3y_Pz_S_S_a-2.0E0*2*I_ERI_Py_Pz_D2x_S_c+2*1*I_ERI_Py_Pz_S_S;
  abcd[988] = 4.0E0*I_ERI_F2yz_Pz_D2x_S_ac-2.0E0*1*I_ERI_F2yz_Pz_S_S_a-2.0E0*1*I_ERI_Pz_Pz_D2x_S_c+1*I_ERI_Pz_Pz_S_S;
  abcd[989] = 4.0E0*I_ERI_Fy2z_Pz_D2x_S_ac-2.0E0*1*I_ERI_Fy2z_Pz_S_S_a;
  abcd[990] = 4.0E0*I_ERI_F2xy_Px_Dxy_S_ac;
  abcd[991] = 4.0E0*I_ERI_Fx2y_Px_Dxy_S_ac-2.0E0*1*I_ERI_Px_Px_Dxy_S_c;
  abcd[992] = 4.0E0*I_ERI_Fxyz_Px_Dxy_S_ac;
  abcd[993] = 4.0E0*I_ERI_F3y_Px_Dxy_S_ac-2.0E0*2*I_ERI_Py_Px_Dxy_S_c;
  abcd[994] = 4.0E0*I_ERI_F2yz_Px_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Px_Dxy_S_c;
  abcd[995] = 4.0E0*I_ERI_Fy2z_Px_Dxy_S_ac;
  abcd[996] = 4.0E0*I_ERI_F2xy_Py_Dxy_S_ac;
  abcd[997] = 4.0E0*I_ERI_Fx2y_Py_Dxy_S_ac-2.0E0*1*I_ERI_Px_Py_Dxy_S_c;
  abcd[998] = 4.0E0*I_ERI_Fxyz_Py_Dxy_S_ac;
  abcd[999] = 4.0E0*I_ERI_F3y_Py_Dxy_S_ac-2.0E0*2*I_ERI_Py_Py_Dxy_S_c;
  abcd[1000] = 4.0E0*I_ERI_F2yz_Py_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Py_Dxy_S_c;
  abcd[1001] = 4.0E0*I_ERI_Fy2z_Py_Dxy_S_ac;
  abcd[1002] = 4.0E0*I_ERI_F2xy_Pz_Dxy_S_ac;
  abcd[1003] = 4.0E0*I_ERI_Fx2y_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Px_Pz_Dxy_S_c;
  abcd[1004] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_S_ac;
  abcd[1005] = 4.0E0*I_ERI_F3y_Pz_Dxy_S_ac-2.0E0*2*I_ERI_Py_Pz_Dxy_S_c;
  abcd[1006] = 4.0E0*I_ERI_F2yz_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dxy_S_c;
  abcd[1007] = 4.0E0*I_ERI_Fy2z_Pz_Dxy_S_ac;
  abcd[1008] = 4.0E0*I_ERI_F2xy_Px_Dxz_S_ac;
  abcd[1009] = 4.0E0*I_ERI_Fx2y_Px_Dxz_S_ac-2.0E0*1*I_ERI_Px_Px_Dxz_S_c;
  abcd[1010] = 4.0E0*I_ERI_Fxyz_Px_Dxz_S_ac;
  abcd[1011] = 4.0E0*I_ERI_F3y_Px_Dxz_S_ac-2.0E0*2*I_ERI_Py_Px_Dxz_S_c;
  abcd[1012] = 4.0E0*I_ERI_F2yz_Px_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Px_Dxz_S_c;
  abcd[1013] = 4.0E0*I_ERI_Fy2z_Px_Dxz_S_ac;
  abcd[1014] = 4.0E0*I_ERI_F2xy_Py_Dxz_S_ac;
  abcd[1015] = 4.0E0*I_ERI_Fx2y_Py_Dxz_S_ac-2.0E0*1*I_ERI_Px_Py_Dxz_S_c;
  abcd[1016] = 4.0E0*I_ERI_Fxyz_Py_Dxz_S_ac;
  abcd[1017] = 4.0E0*I_ERI_F3y_Py_Dxz_S_ac-2.0E0*2*I_ERI_Py_Py_Dxz_S_c;
  abcd[1018] = 4.0E0*I_ERI_F2yz_Py_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Py_Dxz_S_c;
  abcd[1019] = 4.0E0*I_ERI_Fy2z_Py_Dxz_S_ac;
  abcd[1020] = 4.0E0*I_ERI_F2xy_Pz_Dxz_S_ac;
  abcd[1021] = 4.0E0*I_ERI_Fx2y_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Px_Pz_Dxz_S_c;
  abcd[1022] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_S_ac;
  abcd[1023] = 4.0E0*I_ERI_F3y_Pz_Dxz_S_ac-2.0E0*2*I_ERI_Py_Pz_Dxz_S_c;
  abcd[1024] = 4.0E0*I_ERI_F2yz_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dxz_S_c;
  abcd[1025] = 4.0E0*I_ERI_Fy2z_Pz_Dxz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   * RHS shell quartet name: SQ_ERI_P_P_D_S_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  abcd[1026] = 4.0E0*I_ERI_F2xy_Px_Dxy_S_ac;
  abcd[1027] = 4.0E0*I_ERI_Fx2y_Px_Dxy_S_ac-2.0E0*1*I_ERI_Px_Px_Dxy_S_c;
  abcd[1028] = 4.0E0*I_ERI_Fxyz_Px_Dxy_S_ac;
  abcd[1029] = 4.0E0*I_ERI_F3y_Px_Dxy_S_ac-2.0E0*2*I_ERI_Py_Px_Dxy_S_c;
  abcd[1030] = 4.0E0*I_ERI_F2yz_Px_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Px_Dxy_S_c;
  abcd[1031] = 4.0E0*I_ERI_Fy2z_Px_Dxy_S_ac;
  abcd[1032] = 4.0E0*I_ERI_F2xy_Py_Dxy_S_ac;
  abcd[1033] = 4.0E0*I_ERI_Fx2y_Py_Dxy_S_ac-2.0E0*1*I_ERI_Px_Py_Dxy_S_c;
  abcd[1034] = 4.0E0*I_ERI_Fxyz_Py_Dxy_S_ac;
  abcd[1035] = 4.0E0*I_ERI_F3y_Py_Dxy_S_ac-2.0E0*2*I_ERI_Py_Py_Dxy_S_c;
  abcd[1036] = 4.0E0*I_ERI_F2yz_Py_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Py_Dxy_S_c;
  abcd[1037] = 4.0E0*I_ERI_Fy2z_Py_Dxy_S_ac;
  abcd[1038] = 4.0E0*I_ERI_F2xy_Pz_Dxy_S_ac;
  abcd[1039] = 4.0E0*I_ERI_Fx2y_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Px_Pz_Dxy_S_c;
  abcd[1040] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_S_ac;
  abcd[1041] = 4.0E0*I_ERI_F3y_Pz_Dxy_S_ac-2.0E0*2*I_ERI_Py_Pz_Dxy_S_c;
  abcd[1042] = 4.0E0*I_ERI_F2yz_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dxy_S_c;
  abcd[1043] = 4.0E0*I_ERI_Fy2z_Pz_Dxy_S_ac;
  abcd[1044] = 4.0E0*I_ERI_F2xy_Px_D2y_S_ac-2.0E0*1*I_ERI_F2xy_Px_S_S_a;
  abcd[1045] = 4.0E0*I_ERI_Fx2y_Px_D2y_S_ac-2.0E0*1*I_ERI_Fx2y_Px_S_S_a-2.0E0*1*I_ERI_Px_Px_D2y_S_c+1*I_ERI_Px_Px_S_S;
  abcd[1046] = 4.0E0*I_ERI_Fxyz_Px_D2y_S_ac-2.0E0*1*I_ERI_Fxyz_Px_S_S_a;
  abcd[1047] = 4.0E0*I_ERI_F3y_Px_D2y_S_ac-2.0E0*1*I_ERI_F3y_Px_S_S_a-2.0E0*2*I_ERI_Py_Px_D2y_S_c+2*1*I_ERI_Py_Px_S_S;
  abcd[1048] = 4.0E0*I_ERI_F2yz_Px_D2y_S_ac-2.0E0*1*I_ERI_F2yz_Px_S_S_a-2.0E0*1*I_ERI_Pz_Px_D2y_S_c+1*I_ERI_Pz_Px_S_S;
  abcd[1049] = 4.0E0*I_ERI_Fy2z_Px_D2y_S_ac-2.0E0*1*I_ERI_Fy2z_Px_S_S_a;
  abcd[1050] = 4.0E0*I_ERI_F2xy_Py_D2y_S_ac-2.0E0*1*I_ERI_F2xy_Py_S_S_a;
  abcd[1051] = 4.0E0*I_ERI_Fx2y_Py_D2y_S_ac-2.0E0*1*I_ERI_Fx2y_Py_S_S_a-2.0E0*1*I_ERI_Px_Py_D2y_S_c+1*I_ERI_Px_Py_S_S;
  abcd[1052] = 4.0E0*I_ERI_Fxyz_Py_D2y_S_ac-2.0E0*1*I_ERI_Fxyz_Py_S_S_a;
  abcd[1053] = 4.0E0*I_ERI_F3y_Py_D2y_S_ac-2.0E0*1*I_ERI_F3y_Py_S_S_a-2.0E0*2*I_ERI_Py_Py_D2y_S_c+2*1*I_ERI_Py_Py_S_S;
  abcd[1054] = 4.0E0*I_ERI_F2yz_Py_D2y_S_ac-2.0E0*1*I_ERI_F2yz_Py_S_S_a-2.0E0*1*I_ERI_Pz_Py_D2y_S_c+1*I_ERI_Pz_Py_S_S;
  abcd[1055] = 4.0E0*I_ERI_Fy2z_Py_D2y_S_ac-2.0E0*1*I_ERI_Fy2z_Py_S_S_a;
  abcd[1056] = 4.0E0*I_ERI_F2xy_Pz_D2y_S_ac-2.0E0*1*I_ERI_F2xy_Pz_S_S_a;
  abcd[1057] = 4.0E0*I_ERI_Fx2y_Pz_D2y_S_ac-2.0E0*1*I_ERI_Fx2y_Pz_S_S_a-2.0E0*1*I_ERI_Px_Pz_D2y_S_c+1*I_ERI_Px_Pz_S_S;
  abcd[1058] = 4.0E0*I_ERI_Fxyz_Pz_D2y_S_ac-2.0E0*1*I_ERI_Fxyz_Pz_S_S_a;
  abcd[1059] = 4.0E0*I_ERI_F3y_Pz_D2y_S_ac-2.0E0*1*I_ERI_F3y_Pz_S_S_a-2.0E0*2*I_ERI_Py_Pz_D2y_S_c+2*1*I_ERI_Py_Pz_S_S;
  abcd[1060] = 4.0E0*I_ERI_F2yz_Pz_D2y_S_ac-2.0E0*1*I_ERI_F2yz_Pz_S_S_a-2.0E0*1*I_ERI_Pz_Pz_D2y_S_c+1*I_ERI_Pz_Pz_S_S;
  abcd[1061] = 4.0E0*I_ERI_Fy2z_Pz_D2y_S_ac-2.0E0*1*I_ERI_Fy2z_Pz_S_S_a;
  abcd[1062] = 4.0E0*I_ERI_F2xy_Px_Dyz_S_ac;
  abcd[1063] = 4.0E0*I_ERI_Fx2y_Px_Dyz_S_ac-2.0E0*1*I_ERI_Px_Px_Dyz_S_c;
  abcd[1064] = 4.0E0*I_ERI_Fxyz_Px_Dyz_S_ac;
  abcd[1065] = 4.0E0*I_ERI_F3y_Px_Dyz_S_ac-2.0E0*2*I_ERI_Py_Px_Dyz_S_c;
  abcd[1066] = 4.0E0*I_ERI_F2yz_Px_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Px_Dyz_S_c;
  abcd[1067] = 4.0E0*I_ERI_Fy2z_Px_Dyz_S_ac;
  abcd[1068] = 4.0E0*I_ERI_F2xy_Py_Dyz_S_ac;
  abcd[1069] = 4.0E0*I_ERI_Fx2y_Py_Dyz_S_ac-2.0E0*1*I_ERI_Px_Py_Dyz_S_c;
  abcd[1070] = 4.0E0*I_ERI_Fxyz_Py_Dyz_S_ac;
  abcd[1071] = 4.0E0*I_ERI_F3y_Py_Dyz_S_ac-2.0E0*2*I_ERI_Py_Py_Dyz_S_c;
  abcd[1072] = 4.0E0*I_ERI_F2yz_Py_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Py_Dyz_S_c;
  abcd[1073] = 4.0E0*I_ERI_Fy2z_Py_Dyz_S_ac;
  abcd[1074] = 4.0E0*I_ERI_F2xy_Pz_Dyz_S_ac;
  abcd[1075] = 4.0E0*I_ERI_Fx2y_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Px_Pz_Dyz_S_c;
  abcd[1076] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_S_ac;
  abcd[1077] = 4.0E0*I_ERI_F3y_Pz_Dyz_S_ac-2.0E0*2*I_ERI_Py_Pz_Dyz_S_c;
  abcd[1078] = 4.0E0*I_ERI_F2yz_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dyz_S_c;
  abcd[1079] = 4.0E0*I_ERI_Fy2z_Pz_Dyz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   * RHS shell quartet name: SQ_ERI_P_P_D_S_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  abcd[1080] = 4.0E0*I_ERI_F2xy_Px_Dxz_S_ac;
  abcd[1081] = 4.0E0*I_ERI_Fx2y_Px_Dxz_S_ac-2.0E0*1*I_ERI_Px_Px_Dxz_S_c;
  abcd[1082] = 4.0E0*I_ERI_Fxyz_Px_Dxz_S_ac;
  abcd[1083] = 4.0E0*I_ERI_F3y_Px_Dxz_S_ac-2.0E0*2*I_ERI_Py_Px_Dxz_S_c;
  abcd[1084] = 4.0E0*I_ERI_F2yz_Px_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Px_Dxz_S_c;
  abcd[1085] = 4.0E0*I_ERI_Fy2z_Px_Dxz_S_ac;
  abcd[1086] = 4.0E0*I_ERI_F2xy_Py_Dxz_S_ac;
  abcd[1087] = 4.0E0*I_ERI_Fx2y_Py_Dxz_S_ac-2.0E0*1*I_ERI_Px_Py_Dxz_S_c;
  abcd[1088] = 4.0E0*I_ERI_Fxyz_Py_Dxz_S_ac;
  abcd[1089] = 4.0E0*I_ERI_F3y_Py_Dxz_S_ac-2.0E0*2*I_ERI_Py_Py_Dxz_S_c;
  abcd[1090] = 4.0E0*I_ERI_F2yz_Py_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Py_Dxz_S_c;
  abcd[1091] = 4.0E0*I_ERI_Fy2z_Py_Dxz_S_ac;
  abcd[1092] = 4.0E0*I_ERI_F2xy_Pz_Dxz_S_ac;
  abcd[1093] = 4.0E0*I_ERI_Fx2y_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Px_Pz_Dxz_S_c;
  abcd[1094] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_S_ac;
  abcd[1095] = 4.0E0*I_ERI_F3y_Pz_Dxz_S_ac-2.0E0*2*I_ERI_Py_Pz_Dxz_S_c;
  abcd[1096] = 4.0E0*I_ERI_F2yz_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dxz_S_c;
  abcd[1097] = 4.0E0*I_ERI_Fy2z_Pz_Dxz_S_ac;
  abcd[1098] = 4.0E0*I_ERI_F2xy_Px_Dyz_S_ac;
  abcd[1099] = 4.0E0*I_ERI_Fx2y_Px_Dyz_S_ac-2.0E0*1*I_ERI_Px_Px_Dyz_S_c;
  abcd[1100] = 4.0E0*I_ERI_Fxyz_Px_Dyz_S_ac;
  abcd[1101] = 4.0E0*I_ERI_F3y_Px_Dyz_S_ac-2.0E0*2*I_ERI_Py_Px_Dyz_S_c;
  abcd[1102] = 4.0E0*I_ERI_F2yz_Px_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Px_Dyz_S_c;
  abcd[1103] = 4.0E0*I_ERI_Fy2z_Px_Dyz_S_ac;
  abcd[1104] = 4.0E0*I_ERI_F2xy_Py_Dyz_S_ac;
  abcd[1105] = 4.0E0*I_ERI_Fx2y_Py_Dyz_S_ac-2.0E0*1*I_ERI_Px_Py_Dyz_S_c;
  abcd[1106] = 4.0E0*I_ERI_Fxyz_Py_Dyz_S_ac;
  abcd[1107] = 4.0E0*I_ERI_F3y_Py_Dyz_S_ac-2.0E0*2*I_ERI_Py_Py_Dyz_S_c;
  abcd[1108] = 4.0E0*I_ERI_F2yz_Py_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Py_Dyz_S_c;
  abcd[1109] = 4.0E0*I_ERI_Fy2z_Py_Dyz_S_ac;
  abcd[1110] = 4.0E0*I_ERI_F2xy_Pz_Dyz_S_ac;
  abcd[1111] = 4.0E0*I_ERI_Fx2y_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Px_Pz_Dyz_S_c;
  abcd[1112] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_S_ac;
  abcd[1113] = 4.0E0*I_ERI_F3y_Pz_Dyz_S_ac-2.0E0*2*I_ERI_Py_Pz_Dyz_S_c;
  abcd[1114] = 4.0E0*I_ERI_F2yz_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Pz_Pz_Dyz_S_c;
  abcd[1115] = 4.0E0*I_ERI_Fy2z_Pz_Dyz_S_ac;
  abcd[1116] = 4.0E0*I_ERI_F2xy_Px_D2z_S_ac-2.0E0*1*I_ERI_F2xy_Px_S_S_a;
  abcd[1117] = 4.0E0*I_ERI_Fx2y_Px_D2z_S_ac-2.0E0*1*I_ERI_Fx2y_Px_S_S_a-2.0E0*1*I_ERI_Px_Px_D2z_S_c+1*I_ERI_Px_Px_S_S;
  abcd[1118] = 4.0E0*I_ERI_Fxyz_Px_D2z_S_ac-2.0E0*1*I_ERI_Fxyz_Px_S_S_a;
  abcd[1119] = 4.0E0*I_ERI_F3y_Px_D2z_S_ac-2.0E0*1*I_ERI_F3y_Px_S_S_a-2.0E0*2*I_ERI_Py_Px_D2z_S_c+2*1*I_ERI_Py_Px_S_S;
  abcd[1120] = 4.0E0*I_ERI_F2yz_Px_D2z_S_ac-2.0E0*1*I_ERI_F2yz_Px_S_S_a-2.0E0*1*I_ERI_Pz_Px_D2z_S_c+1*I_ERI_Pz_Px_S_S;
  abcd[1121] = 4.0E0*I_ERI_Fy2z_Px_D2z_S_ac-2.0E0*1*I_ERI_Fy2z_Px_S_S_a;
  abcd[1122] = 4.0E0*I_ERI_F2xy_Py_D2z_S_ac-2.0E0*1*I_ERI_F2xy_Py_S_S_a;
  abcd[1123] = 4.0E0*I_ERI_Fx2y_Py_D2z_S_ac-2.0E0*1*I_ERI_Fx2y_Py_S_S_a-2.0E0*1*I_ERI_Px_Py_D2z_S_c+1*I_ERI_Px_Py_S_S;
  abcd[1124] = 4.0E0*I_ERI_Fxyz_Py_D2z_S_ac-2.0E0*1*I_ERI_Fxyz_Py_S_S_a;
  abcd[1125] = 4.0E0*I_ERI_F3y_Py_D2z_S_ac-2.0E0*1*I_ERI_F3y_Py_S_S_a-2.0E0*2*I_ERI_Py_Py_D2z_S_c+2*1*I_ERI_Py_Py_S_S;
  abcd[1126] = 4.0E0*I_ERI_F2yz_Py_D2z_S_ac-2.0E0*1*I_ERI_F2yz_Py_S_S_a-2.0E0*1*I_ERI_Pz_Py_D2z_S_c+1*I_ERI_Pz_Py_S_S;
  abcd[1127] = 4.0E0*I_ERI_Fy2z_Py_D2z_S_ac-2.0E0*1*I_ERI_Fy2z_Py_S_S_a;
  abcd[1128] = 4.0E0*I_ERI_F2xy_Pz_D2z_S_ac-2.0E0*1*I_ERI_F2xy_Pz_S_S_a;
  abcd[1129] = 4.0E0*I_ERI_Fx2y_Pz_D2z_S_ac-2.0E0*1*I_ERI_Fx2y_Pz_S_S_a-2.0E0*1*I_ERI_Px_Pz_D2z_S_c+1*I_ERI_Px_Pz_S_S;
  abcd[1130] = 4.0E0*I_ERI_Fxyz_Pz_D2z_S_ac-2.0E0*1*I_ERI_Fxyz_Pz_S_S_a;
  abcd[1131] = 4.0E0*I_ERI_F3y_Pz_D2z_S_ac-2.0E0*1*I_ERI_F3y_Pz_S_S_a-2.0E0*2*I_ERI_Py_Pz_D2z_S_c+2*1*I_ERI_Py_Pz_S_S;
  abcd[1132] = 4.0E0*I_ERI_F2yz_Pz_D2z_S_ac-2.0E0*1*I_ERI_F2yz_Pz_S_S_a-2.0E0*1*I_ERI_Pz_Pz_D2z_S_c+1*I_ERI_Pz_Pz_S_S;
  abcd[1133] = 4.0E0*I_ERI_Fy2z_Pz_D2z_S_ac-2.0E0*1*I_ERI_Fy2z_Pz_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   * RHS shell quartet name: SQ_ERI_P_P_D_S_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  abcd[1134] = 4.0E0*I_ERI_F2xz_Px_D2x_S_ac-2.0E0*1*I_ERI_F2xz_Px_S_S_a;
  abcd[1135] = 4.0E0*I_ERI_Fxyz_Px_D2x_S_ac-2.0E0*1*I_ERI_Fxyz_Px_S_S_a;
  abcd[1136] = 4.0E0*I_ERI_Fx2z_Px_D2x_S_ac-2.0E0*1*I_ERI_Fx2z_Px_S_S_a-2.0E0*1*I_ERI_Px_Px_D2x_S_c+1*I_ERI_Px_Px_S_S;
  abcd[1137] = 4.0E0*I_ERI_F2yz_Px_D2x_S_ac-2.0E0*1*I_ERI_F2yz_Px_S_S_a;
  abcd[1138] = 4.0E0*I_ERI_Fy2z_Px_D2x_S_ac-2.0E0*1*I_ERI_Fy2z_Px_S_S_a-2.0E0*1*I_ERI_Py_Px_D2x_S_c+1*I_ERI_Py_Px_S_S;
  abcd[1139] = 4.0E0*I_ERI_F3z_Px_D2x_S_ac-2.0E0*1*I_ERI_F3z_Px_S_S_a-2.0E0*2*I_ERI_Pz_Px_D2x_S_c+2*1*I_ERI_Pz_Px_S_S;
  abcd[1140] = 4.0E0*I_ERI_F2xz_Py_D2x_S_ac-2.0E0*1*I_ERI_F2xz_Py_S_S_a;
  abcd[1141] = 4.0E0*I_ERI_Fxyz_Py_D2x_S_ac-2.0E0*1*I_ERI_Fxyz_Py_S_S_a;
  abcd[1142] = 4.0E0*I_ERI_Fx2z_Py_D2x_S_ac-2.0E0*1*I_ERI_Fx2z_Py_S_S_a-2.0E0*1*I_ERI_Px_Py_D2x_S_c+1*I_ERI_Px_Py_S_S;
  abcd[1143] = 4.0E0*I_ERI_F2yz_Py_D2x_S_ac-2.0E0*1*I_ERI_F2yz_Py_S_S_a;
  abcd[1144] = 4.0E0*I_ERI_Fy2z_Py_D2x_S_ac-2.0E0*1*I_ERI_Fy2z_Py_S_S_a-2.0E0*1*I_ERI_Py_Py_D2x_S_c+1*I_ERI_Py_Py_S_S;
  abcd[1145] = 4.0E0*I_ERI_F3z_Py_D2x_S_ac-2.0E0*1*I_ERI_F3z_Py_S_S_a-2.0E0*2*I_ERI_Pz_Py_D2x_S_c+2*1*I_ERI_Pz_Py_S_S;
  abcd[1146] = 4.0E0*I_ERI_F2xz_Pz_D2x_S_ac-2.0E0*1*I_ERI_F2xz_Pz_S_S_a;
  abcd[1147] = 4.0E0*I_ERI_Fxyz_Pz_D2x_S_ac-2.0E0*1*I_ERI_Fxyz_Pz_S_S_a;
  abcd[1148] = 4.0E0*I_ERI_Fx2z_Pz_D2x_S_ac-2.0E0*1*I_ERI_Fx2z_Pz_S_S_a-2.0E0*1*I_ERI_Px_Pz_D2x_S_c+1*I_ERI_Px_Pz_S_S;
  abcd[1149] = 4.0E0*I_ERI_F2yz_Pz_D2x_S_ac-2.0E0*1*I_ERI_F2yz_Pz_S_S_a;
  abcd[1150] = 4.0E0*I_ERI_Fy2z_Pz_D2x_S_ac-2.0E0*1*I_ERI_Fy2z_Pz_S_S_a-2.0E0*1*I_ERI_Py_Pz_D2x_S_c+1*I_ERI_Py_Pz_S_S;
  abcd[1151] = 4.0E0*I_ERI_F3z_Pz_D2x_S_ac-2.0E0*1*I_ERI_F3z_Pz_S_S_a-2.0E0*2*I_ERI_Pz_Pz_D2x_S_c+2*1*I_ERI_Pz_Pz_S_S;
  abcd[1152] = 4.0E0*I_ERI_F2xz_Px_Dxy_S_ac;
  abcd[1153] = 4.0E0*I_ERI_Fxyz_Px_Dxy_S_ac;
  abcd[1154] = 4.0E0*I_ERI_Fx2z_Px_Dxy_S_ac-2.0E0*1*I_ERI_Px_Px_Dxy_S_c;
  abcd[1155] = 4.0E0*I_ERI_F2yz_Px_Dxy_S_ac;
  abcd[1156] = 4.0E0*I_ERI_Fy2z_Px_Dxy_S_ac-2.0E0*1*I_ERI_Py_Px_Dxy_S_c;
  abcd[1157] = 4.0E0*I_ERI_F3z_Px_Dxy_S_ac-2.0E0*2*I_ERI_Pz_Px_Dxy_S_c;
  abcd[1158] = 4.0E0*I_ERI_F2xz_Py_Dxy_S_ac;
  abcd[1159] = 4.0E0*I_ERI_Fxyz_Py_Dxy_S_ac;
  abcd[1160] = 4.0E0*I_ERI_Fx2z_Py_Dxy_S_ac-2.0E0*1*I_ERI_Px_Py_Dxy_S_c;
  abcd[1161] = 4.0E0*I_ERI_F2yz_Py_Dxy_S_ac;
  abcd[1162] = 4.0E0*I_ERI_Fy2z_Py_Dxy_S_ac-2.0E0*1*I_ERI_Py_Py_Dxy_S_c;
  abcd[1163] = 4.0E0*I_ERI_F3z_Py_Dxy_S_ac-2.0E0*2*I_ERI_Pz_Py_Dxy_S_c;
  abcd[1164] = 4.0E0*I_ERI_F2xz_Pz_Dxy_S_ac;
  abcd[1165] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_S_ac;
  abcd[1166] = 4.0E0*I_ERI_Fx2z_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Px_Pz_Dxy_S_c;
  abcd[1167] = 4.0E0*I_ERI_F2yz_Pz_Dxy_S_ac;
  abcd[1168] = 4.0E0*I_ERI_Fy2z_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Py_Pz_Dxy_S_c;
  abcd[1169] = 4.0E0*I_ERI_F3z_Pz_Dxy_S_ac-2.0E0*2*I_ERI_Pz_Pz_Dxy_S_c;
  abcd[1170] = 4.0E0*I_ERI_F2xz_Px_Dxz_S_ac;
  abcd[1171] = 4.0E0*I_ERI_Fxyz_Px_Dxz_S_ac;
  abcd[1172] = 4.0E0*I_ERI_Fx2z_Px_Dxz_S_ac-2.0E0*1*I_ERI_Px_Px_Dxz_S_c;
  abcd[1173] = 4.0E0*I_ERI_F2yz_Px_Dxz_S_ac;
  abcd[1174] = 4.0E0*I_ERI_Fy2z_Px_Dxz_S_ac-2.0E0*1*I_ERI_Py_Px_Dxz_S_c;
  abcd[1175] = 4.0E0*I_ERI_F3z_Px_Dxz_S_ac-2.0E0*2*I_ERI_Pz_Px_Dxz_S_c;
  abcd[1176] = 4.0E0*I_ERI_F2xz_Py_Dxz_S_ac;
  abcd[1177] = 4.0E0*I_ERI_Fxyz_Py_Dxz_S_ac;
  abcd[1178] = 4.0E0*I_ERI_Fx2z_Py_Dxz_S_ac-2.0E0*1*I_ERI_Px_Py_Dxz_S_c;
  abcd[1179] = 4.0E0*I_ERI_F2yz_Py_Dxz_S_ac;
  abcd[1180] = 4.0E0*I_ERI_Fy2z_Py_Dxz_S_ac-2.0E0*1*I_ERI_Py_Py_Dxz_S_c;
  abcd[1181] = 4.0E0*I_ERI_F3z_Py_Dxz_S_ac-2.0E0*2*I_ERI_Pz_Py_Dxz_S_c;
  abcd[1182] = 4.0E0*I_ERI_F2xz_Pz_Dxz_S_ac;
  abcd[1183] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_S_ac;
  abcd[1184] = 4.0E0*I_ERI_Fx2z_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Px_Pz_Dxz_S_c;
  abcd[1185] = 4.0E0*I_ERI_F2yz_Pz_Dxz_S_ac;
  abcd[1186] = 4.0E0*I_ERI_Fy2z_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Py_Pz_Dxz_S_c;
  abcd[1187] = 4.0E0*I_ERI_F3z_Pz_Dxz_S_ac-2.0E0*2*I_ERI_Pz_Pz_Dxz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   * RHS shell quartet name: SQ_ERI_P_P_D_S_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  abcd[1188] = 4.0E0*I_ERI_F2xz_Px_Dxy_S_ac;
  abcd[1189] = 4.0E0*I_ERI_Fxyz_Px_Dxy_S_ac;
  abcd[1190] = 4.0E0*I_ERI_Fx2z_Px_Dxy_S_ac-2.0E0*1*I_ERI_Px_Px_Dxy_S_c;
  abcd[1191] = 4.0E0*I_ERI_F2yz_Px_Dxy_S_ac;
  abcd[1192] = 4.0E0*I_ERI_Fy2z_Px_Dxy_S_ac-2.0E0*1*I_ERI_Py_Px_Dxy_S_c;
  abcd[1193] = 4.0E0*I_ERI_F3z_Px_Dxy_S_ac-2.0E0*2*I_ERI_Pz_Px_Dxy_S_c;
  abcd[1194] = 4.0E0*I_ERI_F2xz_Py_Dxy_S_ac;
  abcd[1195] = 4.0E0*I_ERI_Fxyz_Py_Dxy_S_ac;
  abcd[1196] = 4.0E0*I_ERI_Fx2z_Py_Dxy_S_ac-2.0E0*1*I_ERI_Px_Py_Dxy_S_c;
  abcd[1197] = 4.0E0*I_ERI_F2yz_Py_Dxy_S_ac;
  abcd[1198] = 4.0E0*I_ERI_Fy2z_Py_Dxy_S_ac-2.0E0*1*I_ERI_Py_Py_Dxy_S_c;
  abcd[1199] = 4.0E0*I_ERI_F3z_Py_Dxy_S_ac-2.0E0*2*I_ERI_Pz_Py_Dxy_S_c;
  abcd[1200] = 4.0E0*I_ERI_F2xz_Pz_Dxy_S_ac;
  abcd[1201] = 4.0E0*I_ERI_Fxyz_Pz_Dxy_S_ac;
  abcd[1202] = 4.0E0*I_ERI_Fx2z_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Px_Pz_Dxy_S_c;
  abcd[1203] = 4.0E0*I_ERI_F2yz_Pz_Dxy_S_ac;
  abcd[1204] = 4.0E0*I_ERI_Fy2z_Pz_Dxy_S_ac-2.0E0*1*I_ERI_Py_Pz_Dxy_S_c;
  abcd[1205] = 4.0E0*I_ERI_F3z_Pz_Dxy_S_ac-2.0E0*2*I_ERI_Pz_Pz_Dxy_S_c;
  abcd[1206] = 4.0E0*I_ERI_F2xz_Px_D2y_S_ac-2.0E0*1*I_ERI_F2xz_Px_S_S_a;
  abcd[1207] = 4.0E0*I_ERI_Fxyz_Px_D2y_S_ac-2.0E0*1*I_ERI_Fxyz_Px_S_S_a;
  abcd[1208] = 4.0E0*I_ERI_Fx2z_Px_D2y_S_ac-2.0E0*1*I_ERI_Fx2z_Px_S_S_a-2.0E0*1*I_ERI_Px_Px_D2y_S_c+1*I_ERI_Px_Px_S_S;
  abcd[1209] = 4.0E0*I_ERI_F2yz_Px_D2y_S_ac-2.0E0*1*I_ERI_F2yz_Px_S_S_a;
  abcd[1210] = 4.0E0*I_ERI_Fy2z_Px_D2y_S_ac-2.0E0*1*I_ERI_Fy2z_Px_S_S_a-2.0E0*1*I_ERI_Py_Px_D2y_S_c+1*I_ERI_Py_Px_S_S;
  abcd[1211] = 4.0E0*I_ERI_F3z_Px_D2y_S_ac-2.0E0*1*I_ERI_F3z_Px_S_S_a-2.0E0*2*I_ERI_Pz_Px_D2y_S_c+2*1*I_ERI_Pz_Px_S_S;
  abcd[1212] = 4.0E0*I_ERI_F2xz_Py_D2y_S_ac-2.0E0*1*I_ERI_F2xz_Py_S_S_a;
  abcd[1213] = 4.0E0*I_ERI_Fxyz_Py_D2y_S_ac-2.0E0*1*I_ERI_Fxyz_Py_S_S_a;
  abcd[1214] = 4.0E0*I_ERI_Fx2z_Py_D2y_S_ac-2.0E0*1*I_ERI_Fx2z_Py_S_S_a-2.0E0*1*I_ERI_Px_Py_D2y_S_c+1*I_ERI_Px_Py_S_S;
  abcd[1215] = 4.0E0*I_ERI_F2yz_Py_D2y_S_ac-2.0E0*1*I_ERI_F2yz_Py_S_S_a;
  abcd[1216] = 4.0E0*I_ERI_Fy2z_Py_D2y_S_ac-2.0E0*1*I_ERI_Fy2z_Py_S_S_a-2.0E0*1*I_ERI_Py_Py_D2y_S_c+1*I_ERI_Py_Py_S_S;
  abcd[1217] = 4.0E0*I_ERI_F3z_Py_D2y_S_ac-2.0E0*1*I_ERI_F3z_Py_S_S_a-2.0E0*2*I_ERI_Pz_Py_D2y_S_c+2*1*I_ERI_Pz_Py_S_S;
  abcd[1218] = 4.0E0*I_ERI_F2xz_Pz_D2y_S_ac-2.0E0*1*I_ERI_F2xz_Pz_S_S_a;
  abcd[1219] = 4.0E0*I_ERI_Fxyz_Pz_D2y_S_ac-2.0E0*1*I_ERI_Fxyz_Pz_S_S_a;
  abcd[1220] = 4.0E0*I_ERI_Fx2z_Pz_D2y_S_ac-2.0E0*1*I_ERI_Fx2z_Pz_S_S_a-2.0E0*1*I_ERI_Px_Pz_D2y_S_c+1*I_ERI_Px_Pz_S_S;
  abcd[1221] = 4.0E0*I_ERI_F2yz_Pz_D2y_S_ac-2.0E0*1*I_ERI_F2yz_Pz_S_S_a;
  abcd[1222] = 4.0E0*I_ERI_Fy2z_Pz_D2y_S_ac-2.0E0*1*I_ERI_Fy2z_Pz_S_S_a-2.0E0*1*I_ERI_Py_Pz_D2y_S_c+1*I_ERI_Py_Pz_S_S;
  abcd[1223] = 4.0E0*I_ERI_F3z_Pz_D2y_S_ac-2.0E0*1*I_ERI_F3z_Pz_S_S_a-2.0E0*2*I_ERI_Pz_Pz_D2y_S_c+2*1*I_ERI_Pz_Pz_S_S;
  abcd[1224] = 4.0E0*I_ERI_F2xz_Px_Dyz_S_ac;
  abcd[1225] = 4.0E0*I_ERI_Fxyz_Px_Dyz_S_ac;
  abcd[1226] = 4.0E0*I_ERI_Fx2z_Px_Dyz_S_ac-2.0E0*1*I_ERI_Px_Px_Dyz_S_c;
  abcd[1227] = 4.0E0*I_ERI_F2yz_Px_Dyz_S_ac;
  abcd[1228] = 4.0E0*I_ERI_Fy2z_Px_Dyz_S_ac-2.0E0*1*I_ERI_Py_Px_Dyz_S_c;
  abcd[1229] = 4.0E0*I_ERI_F3z_Px_Dyz_S_ac-2.0E0*2*I_ERI_Pz_Px_Dyz_S_c;
  abcd[1230] = 4.0E0*I_ERI_F2xz_Py_Dyz_S_ac;
  abcd[1231] = 4.0E0*I_ERI_Fxyz_Py_Dyz_S_ac;
  abcd[1232] = 4.0E0*I_ERI_Fx2z_Py_Dyz_S_ac-2.0E0*1*I_ERI_Px_Py_Dyz_S_c;
  abcd[1233] = 4.0E0*I_ERI_F2yz_Py_Dyz_S_ac;
  abcd[1234] = 4.0E0*I_ERI_Fy2z_Py_Dyz_S_ac-2.0E0*1*I_ERI_Py_Py_Dyz_S_c;
  abcd[1235] = 4.0E0*I_ERI_F3z_Py_Dyz_S_ac-2.0E0*2*I_ERI_Pz_Py_Dyz_S_c;
  abcd[1236] = 4.0E0*I_ERI_F2xz_Pz_Dyz_S_ac;
  abcd[1237] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_S_ac;
  abcd[1238] = 4.0E0*I_ERI_Fx2z_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Px_Pz_Dyz_S_c;
  abcd[1239] = 4.0E0*I_ERI_F2yz_Pz_Dyz_S_ac;
  abcd[1240] = 4.0E0*I_ERI_Fy2z_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Py_Pz_Dyz_S_c;
  abcd[1241] = 4.0E0*I_ERI_F3z_Pz_Dyz_S_ac-2.0E0*2*I_ERI_Pz_Pz_Dyz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_D_S_ac
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   * RHS shell quartet name: SQ_ERI_P_P_D_S_c
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  abcd[1242] = 4.0E0*I_ERI_F2xz_Px_Dxz_S_ac;
  abcd[1243] = 4.0E0*I_ERI_Fxyz_Px_Dxz_S_ac;
  abcd[1244] = 4.0E0*I_ERI_Fx2z_Px_Dxz_S_ac-2.0E0*1*I_ERI_Px_Px_Dxz_S_c;
  abcd[1245] = 4.0E0*I_ERI_F2yz_Px_Dxz_S_ac;
  abcd[1246] = 4.0E0*I_ERI_Fy2z_Px_Dxz_S_ac-2.0E0*1*I_ERI_Py_Px_Dxz_S_c;
  abcd[1247] = 4.0E0*I_ERI_F3z_Px_Dxz_S_ac-2.0E0*2*I_ERI_Pz_Px_Dxz_S_c;
  abcd[1248] = 4.0E0*I_ERI_F2xz_Py_Dxz_S_ac;
  abcd[1249] = 4.0E0*I_ERI_Fxyz_Py_Dxz_S_ac;
  abcd[1250] = 4.0E0*I_ERI_Fx2z_Py_Dxz_S_ac-2.0E0*1*I_ERI_Px_Py_Dxz_S_c;
  abcd[1251] = 4.0E0*I_ERI_F2yz_Py_Dxz_S_ac;
  abcd[1252] = 4.0E0*I_ERI_Fy2z_Py_Dxz_S_ac-2.0E0*1*I_ERI_Py_Py_Dxz_S_c;
  abcd[1253] = 4.0E0*I_ERI_F3z_Py_Dxz_S_ac-2.0E0*2*I_ERI_Pz_Py_Dxz_S_c;
  abcd[1254] = 4.0E0*I_ERI_F2xz_Pz_Dxz_S_ac;
  abcd[1255] = 4.0E0*I_ERI_Fxyz_Pz_Dxz_S_ac;
  abcd[1256] = 4.0E0*I_ERI_Fx2z_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Px_Pz_Dxz_S_c;
  abcd[1257] = 4.0E0*I_ERI_F2yz_Pz_Dxz_S_ac;
  abcd[1258] = 4.0E0*I_ERI_Fy2z_Pz_Dxz_S_ac-2.0E0*1*I_ERI_Py_Pz_Dxz_S_c;
  abcd[1259] = 4.0E0*I_ERI_F3z_Pz_Dxz_S_ac-2.0E0*2*I_ERI_Pz_Pz_Dxz_S_c;
  abcd[1260] = 4.0E0*I_ERI_F2xz_Px_Dyz_S_ac;
  abcd[1261] = 4.0E0*I_ERI_Fxyz_Px_Dyz_S_ac;
  abcd[1262] = 4.0E0*I_ERI_Fx2z_Px_Dyz_S_ac-2.0E0*1*I_ERI_Px_Px_Dyz_S_c;
  abcd[1263] = 4.0E0*I_ERI_F2yz_Px_Dyz_S_ac;
  abcd[1264] = 4.0E0*I_ERI_Fy2z_Px_Dyz_S_ac-2.0E0*1*I_ERI_Py_Px_Dyz_S_c;
  abcd[1265] = 4.0E0*I_ERI_F3z_Px_Dyz_S_ac-2.0E0*2*I_ERI_Pz_Px_Dyz_S_c;
  abcd[1266] = 4.0E0*I_ERI_F2xz_Py_Dyz_S_ac;
  abcd[1267] = 4.0E0*I_ERI_Fxyz_Py_Dyz_S_ac;
  abcd[1268] = 4.0E0*I_ERI_Fx2z_Py_Dyz_S_ac-2.0E0*1*I_ERI_Px_Py_Dyz_S_c;
  abcd[1269] = 4.0E0*I_ERI_F2yz_Py_Dyz_S_ac;
  abcd[1270] = 4.0E0*I_ERI_Fy2z_Py_Dyz_S_ac-2.0E0*1*I_ERI_Py_Py_Dyz_S_c;
  abcd[1271] = 4.0E0*I_ERI_F3z_Py_Dyz_S_ac-2.0E0*2*I_ERI_Pz_Py_Dyz_S_c;
  abcd[1272] = 4.0E0*I_ERI_F2xz_Pz_Dyz_S_ac;
  abcd[1273] = 4.0E0*I_ERI_Fxyz_Pz_Dyz_S_ac;
  abcd[1274] = 4.0E0*I_ERI_Fx2z_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Px_Pz_Dyz_S_c;
  abcd[1275] = 4.0E0*I_ERI_F2yz_Pz_Dyz_S_ac;
  abcd[1276] = 4.0E0*I_ERI_Fy2z_Pz_Dyz_S_ac-2.0E0*1*I_ERI_Py_Pz_Dyz_S_c;
  abcd[1277] = 4.0E0*I_ERI_F3z_Pz_Dyz_S_ac-2.0E0*2*I_ERI_Pz_Pz_Dyz_S_c;
  abcd[1278] = 4.0E0*I_ERI_F2xz_Px_D2z_S_ac-2.0E0*1*I_ERI_F2xz_Px_S_S_a;
  abcd[1279] = 4.0E0*I_ERI_Fxyz_Px_D2z_S_ac-2.0E0*1*I_ERI_Fxyz_Px_S_S_a;
  abcd[1280] = 4.0E0*I_ERI_Fx2z_Px_D2z_S_ac-2.0E0*1*I_ERI_Fx2z_Px_S_S_a-2.0E0*1*I_ERI_Px_Px_D2z_S_c+1*I_ERI_Px_Px_S_S;
  abcd[1281] = 4.0E0*I_ERI_F2yz_Px_D2z_S_ac-2.0E0*1*I_ERI_F2yz_Px_S_S_a;
  abcd[1282] = 4.0E0*I_ERI_Fy2z_Px_D2z_S_ac-2.0E0*1*I_ERI_Fy2z_Px_S_S_a-2.0E0*1*I_ERI_Py_Px_D2z_S_c+1*I_ERI_Py_Px_S_S;
  abcd[1283] = 4.0E0*I_ERI_F3z_Px_D2z_S_ac-2.0E0*1*I_ERI_F3z_Px_S_S_a-2.0E0*2*I_ERI_Pz_Px_D2z_S_c+2*1*I_ERI_Pz_Px_S_S;
  abcd[1284] = 4.0E0*I_ERI_F2xz_Py_D2z_S_ac-2.0E0*1*I_ERI_F2xz_Py_S_S_a;
  abcd[1285] = 4.0E0*I_ERI_Fxyz_Py_D2z_S_ac-2.0E0*1*I_ERI_Fxyz_Py_S_S_a;
  abcd[1286] = 4.0E0*I_ERI_Fx2z_Py_D2z_S_ac-2.0E0*1*I_ERI_Fx2z_Py_S_S_a-2.0E0*1*I_ERI_Px_Py_D2z_S_c+1*I_ERI_Px_Py_S_S;
  abcd[1287] = 4.0E0*I_ERI_F2yz_Py_D2z_S_ac-2.0E0*1*I_ERI_F2yz_Py_S_S_a;
  abcd[1288] = 4.0E0*I_ERI_Fy2z_Py_D2z_S_ac-2.0E0*1*I_ERI_Fy2z_Py_S_S_a-2.0E0*1*I_ERI_Py_Py_D2z_S_c+1*I_ERI_Py_Py_S_S;
  abcd[1289] = 4.0E0*I_ERI_F3z_Py_D2z_S_ac-2.0E0*1*I_ERI_F3z_Py_S_S_a-2.0E0*2*I_ERI_Pz_Py_D2z_S_c+2*1*I_ERI_Pz_Py_S_S;
  abcd[1290] = 4.0E0*I_ERI_F2xz_Pz_D2z_S_ac-2.0E0*1*I_ERI_F2xz_Pz_S_S_a;
  abcd[1291] = 4.0E0*I_ERI_Fxyz_Pz_D2z_S_ac-2.0E0*1*I_ERI_Fxyz_Pz_S_S_a;
  abcd[1292] = 4.0E0*I_ERI_Fx2z_Pz_D2z_S_ac-2.0E0*1*I_ERI_Fx2z_Pz_S_S_a-2.0E0*1*I_ERI_Px_Pz_D2z_S_c+1*I_ERI_Px_Pz_S_S;
  abcd[1293] = 4.0E0*I_ERI_F2yz_Pz_D2z_S_ac-2.0E0*1*I_ERI_F2yz_Pz_S_S_a;
  abcd[1294] = 4.0E0*I_ERI_Fy2z_Pz_D2z_S_ac-2.0E0*1*I_ERI_Fy2z_Pz_S_S_a-2.0E0*1*I_ERI_Py_Pz_D2z_S_c+1*I_ERI_Py_Pz_S_S;
  abcd[1295] = 4.0E0*I_ERI_F3z_Pz_D2z_S_ac-2.0E0*1*I_ERI_F3z_Pz_S_S_a-2.0E0*2*I_ERI_Pz_Pz_D2z_S_c+2*1*I_ERI_Pz_Pz_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_P_S_bb
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   ************************************************************/
  abcd[1296] = 4.0E0*I_ERI_D2x_F3x_Px_S_bb-2.0E0*1*I_ERI_D2x_Px_Px_S_b-2.0E0*2*I_ERI_D2x_Px_Px_S_b;
  abcd[1297] = 4.0E0*I_ERI_Dxy_F3x_Px_S_bb-2.0E0*1*I_ERI_Dxy_Px_Px_S_b-2.0E0*2*I_ERI_Dxy_Px_Px_S_b;
  abcd[1298] = 4.0E0*I_ERI_Dxz_F3x_Px_S_bb-2.0E0*1*I_ERI_Dxz_Px_Px_S_b-2.0E0*2*I_ERI_Dxz_Px_Px_S_b;
  abcd[1299] = 4.0E0*I_ERI_D2y_F3x_Px_S_bb-2.0E0*1*I_ERI_D2y_Px_Px_S_b-2.0E0*2*I_ERI_D2y_Px_Px_S_b;
  abcd[1300] = 4.0E0*I_ERI_Dyz_F3x_Px_S_bb-2.0E0*1*I_ERI_Dyz_Px_Px_S_b-2.0E0*2*I_ERI_Dyz_Px_Px_S_b;
  abcd[1301] = 4.0E0*I_ERI_D2z_F3x_Px_S_bb-2.0E0*1*I_ERI_D2z_Px_Px_S_b-2.0E0*2*I_ERI_D2z_Px_Px_S_b;
  abcd[1302] = 4.0E0*I_ERI_D2x_F2xy_Px_S_bb-2.0E0*1*I_ERI_D2x_Py_Px_S_b;
  abcd[1303] = 4.0E0*I_ERI_Dxy_F2xy_Px_S_bb-2.0E0*1*I_ERI_Dxy_Py_Px_S_b;
  abcd[1304] = 4.0E0*I_ERI_Dxz_F2xy_Px_S_bb-2.0E0*1*I_ERI_Dxz_Py_Px_S_b;
  abcd[1305] = 4.0E0*I_ERI_D2y_F2xy_Px_S_bb-2.0E0*1*I_ERI_D2y_Py_Px_S_b;
  abcd[1306] = 4.0E0*I_ERI_Dyz_F2xy_Px_S_bb-2.0E0*1*I_ERI_Dyz_Py_Px_S_b;
  abcd[1307] = 4.0E0*I_ERI_D2z_F2xy_Px_S_bb-2.0E0*1*I_ERI_D2z_Py_Px_S_b;
  abcd[1308] = 4.0E0*I_ERI_D2x_F2xz_Px_S_bb-2.0E0*1*I_ERI_D2x_Pz_Px_S_b;
  abcd[1309] = 4.0E0*I_ERI_Dxy_F2xz_Px_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Px_S_b;
  abcd[1310] = 4.0E0*I_ERI_Dxz_F2xz_Px_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Px_S_b;
  abcd[1311] = 4.0E0*I_ERI_D2y_F2xz_Px_S_bb-2.0E0*1*I_ERI_D2y_Pz_Px_S_b;
  abcd[1312] = 4.0E0*I_ERI_Dyz_F2xz_Px_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Px_S_b;
  abcd[1313] = 4.0E0*I_ERI_D2z_F2xz_Px_S_bb-2.0E0*1*I_ERI_D2z_Pz_Px_S_b;
  abcd[1314] = 4.0E0*I_ERI_D2x_F3x_Py_S_bb-2.0E0*1*I_ERI_D2x_Px_Py_S_b-2.0E0*2*I_ERI_D2x_Px_Py_S_b;
  abcd[1315] = 4.0E0*I_ERI_Dxy_F3x_Py_S_bb-2.0E0*1*I_ERI_Dxy_Px_Py_S_b-2.0E0*2*I_ERI_Dxy_Px_Py_S_b;
  abcd[1316] = 4.0E0*I_ERI_Dxz_F3x_Py_S_bb-2.0E0*1*I_ERI_Dxz_Px_Py_S_b-2.0E0*2*I_ERI_Dxz_Px_Py_S_b;
  abcd[1317] = 4.0E0*I_ERI_D2y_F3x_Py_S_bb-2.0E0*1*I_ERI_D2y_Px_Py_S_b-2.0E0*2*I_ERI_D2y_Px_Py_S_b;
  abcd[1318] = 4.0E0*I_ERI_Dyz_F3x_Py_S_bb-2.0E0*1*I_ERI_Dyz_Px_Py_S_b-2.0E0*2*I_ERI_Dyz_Px_Py_S_b;
  abcd[1319] = 4.0E0*I_ERI_D2z_F3x_Py_S_bb-2.0E0*1*I_ERI_D2z_Px_Py_S_b-2.0E0*2*I_ERI_D2z_Px_Py_S_b;
  abcd[1320] = 4.0E0*I_ERI_D2x_F2xy_Py_S_bb-2.0E0*1*I_ERI_D2x_Py_Py_S_b;
  abcd[1321] = 4.0E0*I_ERI_Dxy_F2xy_Py_S_bb-2.0E0*1*I_ERI_Dxy_Py_Py_S_b;
  abcd[1322] = 4.0E0*I_ERI_Dxz_F2xy_Py_S_bb-2.0E0*1*I_ERI_Dxz_Py_Py_S_b;
  abcd[1323] = 4.0E0*I_ERI_D2y_F2xy_Py_S_bb-2.0E0*1*I_ERI_D2y_Py_Py_S_b;
  abcd[1324] = 4.0E0*I_ERI_Dyz_F2xy_Py_S_bb-2.0E0*1*I_ERI_Dyz_Py_Py_S_b;
  abcd[1325] = 4.0E0*I_ERI_D2z_F2xy_Py_S_bb-2.0E0*1*I_ERI_D2z_Py_Py_S_b;
  abcd[1326] = 4.0E0*I_ERI_D2x_F2xz_Py_S_bb-2.0E0*1*I_ERI_D2x_Pz_Py_S_b;
  abcd[1327] = 4.0E0*I_ERI_Dxy_F2xz_Py_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Py_S_b;
  abcd[1328] = 4.0E0*I_ERI_Dxz_F2xz_Py_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Py_S_b;
  abcd[1329] = 4.0E0*I_ERI_D2y_F2xz_Py_S_bb-2.0E0*1*I_ERI_D2y_Pz_Py_S_b;
  abcd[1330] = 4.0E0*I_ERI_Dyz_F2xz_Py_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Py_S_b;
  abcd[1331] = 4.0E0*I_ERI_D2z_F2xz_Py_S_bb-2.0E0*1*I_ERI_D2z_Pz_Py_S_b;
  abcd[1332] = 4.0E0*I_ERI_D2x_F3x_Pz_S_bb-2.0E0*1*I_ERI_D2x_Px_Pz_S_b-2.0E0*2*I_ERI_D2x_Px_Pz_S_b;
  abcd[1333] = 4.0E0*I_ERI_Dxy_F3x_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Px_Pz_S_b-2.0E0*2*I_ERI_Dxy_Px_Pz_S_b;
  abcd[1334] = 4.0E0*I_ERI_Dxz_F3x_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Px_Pz_S_b-2.0E0*2*I_ERI_Dxz_Px_Pz_S_b;
  abcd[1335] = 4.0E0*I_ERI_D2y_F3x_Pz_S_bb-2.0E0*1*I_ERI_D2y_Px_Pz_S_b-2.0E0*2*I_ERI_D2y_Px_Pz_S_b;
  abcd[1336] = 4.0E0*I_ERI_Dyz_F3x_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Px_Pz_S_b-2.0E0*2*I_ERI_Dyz_Px_Pz_S_b;
  abcd[1337] = 4.0E0*I_ERI_D2z_F3x_Pz_S_bb-2.0E0*1*I_ERI_D2z_Px_Pz_S_b-2.0E0*2*I_ERI_D2z_Px_Pz_S_b;
  abcd[1338] = 4.0E0*I_ERI_D2x_F2xy_Pz_S_bb-2.0E0*1*I_ERI_D2x_Py_Pz_S_b;
  abcd[1339] = 4.0E0*I_ERI_Dxy_F2xy_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Py_Pz_S_b;
  abcd[1340] = 4.0E0*I_ERI_Dxz_F2xy_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Py_Pz_S_b;
  abcd[1341] = 4.0E0*I_ERI_D2y_F2xy_Pz_S_bb-2.0E0*1*I_ERI_D2y_Py_Pz_S_b;
  abcd[1342] = 4.0E0*I_ERI_Dyz_F2xy_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Py_Pz_S_b;
  abcd[1343] = 4.0E0*I_ERI_D2z_F2xy_Pz_S_bb-2.0E0*1*I_ERI_D2z_Py_Pz_S_b;
  abcd[1344] = 4.0E0*I_ERI_D2x_F2xz_Pz_S_bb-2.0E0*1*I_ERI_D2x_Pz_Pz_S_b;
  abcd[1345] = 4.0E0*I_ERI_Dxy_F2xz_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_b;
  abcd[1346] = 4.0E0*I_ERI_Dxz_F2xz_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_b;
  abcd[1347] = 4.0E0*I_ERI_D2y_F2xz_Pz_S_bb-2.0E0*1*I_ERI_D2y_Pz_Pz_S_b;
  abcd[1348] = 4.0E0*I_ERI_Dyz_F2xz_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_b;
  abcd[1349] = 4.0E0*I_ERI_D2z_F2xz_Pz_S_bb-2.0E0*1*I_ERI_D2z_Pz_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_P_S_bb
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   ************************************************************/
  abcd[1350] = 4.0E0*I_ERI_D2x_F2xy_Px_S_bb-2.0E0*1*I_ERI_D2x_Py_Px_S_b;
  abcd[1351] = 4.0E0*I_ERI_Dxy_F2xy_Px_S_bb-2.0E0*1*I_ERI_Dxy_Py_Px_S_b;
  abcd[1352] = 4.0E0*I_ERI_Dxz_F2xy_Px_S_bb-2.0E0*1*I_ERI_Dxz_Py_Px_S_b;
  abcd[1353] = 4.0E0*I_ERI_D2y_F2xy_Px_S_bb-2.0E0*1*I_ERI_D2y_Py_Px_S_b;
  abcd[1354] = 4.0E0*I_ERI_Dyz_F2xy_Px_S_bb-2.0E0*1*I_ERI_Dyz_Py_Px_S_b;
  abcd[1355] = 4.0E0*I_ERI_D2z_F2xy_Px_S_bb-2.0E0*1*I_ERI_D2z_Py_Px_S_b;
  abcd[1356] = 4.0E0*I_ERI_D2x_Fx2y_Px_S_bb-2.0E0*1*I_ERI_D2x_Px_Px_S_b;
  abcd[1357] = 4.0E0*I_ERI_Dxy_Fx2y_Px_S_bb-2.0E0*1*I_ERI_Dxy_Px_Px_S_b;
  abcd[1358] = 4.0E0*I_ERI_Dxz_Fx2y_Px_S_bb-2.0E0*1*I_ERI_Dxz_Px_Px_S_b;
  abcd[1359] = 4.0E0*I_ERI_D2y_Fx2y_Px_S_bb-2.0E0*1*I_ERI_D2y_Px_Px_S_b;
  abcd[1360] = 4.0E0*I_ERI_Dyz_Fx2y_Px_S_bb-2.0E0*1*I_ERI_Dyz_Px_Px_S_b;
  abcd[1361] = 4.0E0*I_ERI_D2z_Fx2y_Px_S_bb-2.0E0*1*I_ERI_D2z_Px_Px_S_b;
  abcd[1362] = 4.0E0*I_ERI_D2x_Fxyz_Px_S_bb;
  abcd[1363] = 4.0E0*I_ERI_Dxy_Fxyz_Px_S_bb;
  abcd[1364] = 4.0E0*I_ERI_Dxz_Fxyz_Px_S_bb;
  abcd[1365] = 4.0E0*I_ERI_D2y_Fxyz_Px_S_bb;
  abcd[1366] = 4.0E0*I_ERI_Dyz_Fxyz_Px_S_bb;
  abcd[1367] = 4.0E0*I_ERI_D2z_Fxyz_Px_S_bb;
  abcd[1368] = 4.0E0*I_ERI_D2x_F2xy_Py_S_bb-2.0E0*1*I_ERI_D2x_Py_Py_S_b;
  abcd[1369] = 4.0E0*I_ERI_Dxy_F2xy_Py_S_bb-2.0E0*1*I_ERI_Dxy_Py_Py_S_b;
  abcd[1370] = 4.0E0*I_ERI_Dxz_F2xy_Py_S_bb-2.0E0*1*I_ERI_Dxz_Py_Py_S_b;
  abcd[1371] = 4.0E0*I_ERI_D2y_F2xy_Py_S_bb-2.0E0*1*I_ERI_D2y_Py_Py_S_b;
  abcd[1372] = 4.0E0*I_ERI_Dyz_F2xy_Py_S_bb-2.0E0*1*I_ERI_Dyz_Py_Py_S_b;
  abcd[1373] = 4.0E0*I_ERI_D2z_F2xy_Py_S_bb-2.0E0*1*I_ERI_D2z_Py_Py_S_b;
  abcd[1374] = 4.0E0*I_ERI_D2x_Fx2y_Py_S_bb-2.0E0*1*I_ERI_D2x_Px_Py_S_b;
  abcd[1375] = 4.0E0*I_ERI_Dxy_Fx2y_Py_S_bb-2.0E0*1*I_ERI_Dxy_Px_Py_S_b;
  abcd[1376] = 4.0E0*I_ERI_Dxz_Fx2y_Py_S_bb-2.0E0*1*I_ERI_Dxz_Px_Py_S_b;
  abcd[1377] = 4.0E0*I_ERI_D2y_Fx2y_Py_S_bb-2.0E0*1*I_ERI_D2y_Px_Py_S_b;
  abcd[1378] = 4.0E0*I_ERI_Dyz_Fx2y_Py_S_bb-2.0E0*1*I_ERI_Dyz_Px_Py_S_b;
  abcd[1379] = 4.0E0*I_ERI_D2z_Fx2y_Py_S_bb-2.0E0*1*I_ERI_D2z_Px_Py_S_b;
  abcd[1380] = 4.0E0*I_ERI_D2x_Fxyz_Py_S_bb;
  abcd[1381] = 4.0E0*I_ERI_Dxy_Fxyz_Py_S_bb;
  abcd[1382] = 4.0E0*I_ERI_Dxz_Fxyz_Py_S_bb;
  abcd[1383] = 4.0E0*I_ERI_D2y_Fxyz_Py_S_bb;
  abcd[1384] = 4.0E0*I_ERI_Dyz_Fxyz_Py_S_bb;
  abcd[1385] = 4.0E0*I_ERI_D2z_Fxyz_Py_S_bb;
  abcd[1386] = 4.0E0*I_ERI_D2x_F2xy_Pz_S_bb-2.0E0*1*I_ERI_D2x_Py_Pz_S_b;
  abcd[1387] = 4.0E0*I_ERI_Dxy_F2xy_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Py_Pz_S_b;
  abcd[1388] = 4.0E0*I_ERI_Dxz_F2xy_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Py_Pz_S_b;
  abcd[1389] = 4.0E0*I_ERI_D2y_F2xy_Pz_S_bb-2.0E0*1*I_ERI_D2y_Py_Pz_S_b;
  abcd[1390] = 4.0E0*I_ERI_Dyz_F2xy_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Py_Pz_S_b;
  abcd[1391] = 4.0E0*I_ERI_D2z_F2xy_Pz_S_bb-2.0E0*1*I_ERI_D2z_Py_Pz_S_b;
  abcd[1392] = 4.0E0*I_ERI_D2x_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_D2x_Px_Pz_S_b;
  abcd[1393] = 4.0E0*I_ERI_Dxy_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Px_Pz_S_b;
  abcd[1394] = 4.0E0*I_ERI_Dxz_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Px_Pz_S_b;
  abcd[1395] = 4.0E0*I_ERI_D2y_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_D2y_Px_Pz_S_b;
  abcd[1396] = 4.0E0*I_ERI_Dyz_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Px_Pz_S_b;
  abcd[1397] = 4.0E0*I_ERI_D2z_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_D2z_Px_Pz_S_b;
  abcd[1398] = 4.0E0*I_ERI_D2x_Fxyz_Pz_S_bb;
  abcd[1399] = 4.0E0*I_ERI_Dxy_Fxyz_Pz_S_bb;
  abcd[1400] = 4.0E0*I_ERI_Dxz_Fxyz_Pz_S_bb;
  abcd[1401] = 4.0E0*I_ERI_D2y_Fxyz_Pz_S_bb;
  abcd[1402] = 4.0E0*I_ERI_Dyz_Fxyz_Pz_S_bb;
  abcd[1403] = 4.0E0*I_ERI_D2z_Fxyz_Pz_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_P_S_bb
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   ************************************************************/
  abcd[1404] = 4.0E0*I_ERI_D2x_F2xz_Px_S_bb-2.0E0*1*I_ERI_D2x_Pz_Px_S_b;
  abcd[1405] = 4.0E0*I_ERI_Dxy_F2xz_Px_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Px_S_b;
  abcd[1406] = 4.0E0*I_ERI_Dxz_F2xz_Px_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Px_S_b;
  abcd[1407] = 4.0E0*I_ERI_D2y_F2xz_Px_S_bb-2.0E0*1*I_ERI_D2y_Pz_Px_S_b;
  abcd[1408] = 4.0E0*I_ERI_Dyz_F2xz_Px_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Px_S_b;
  abcd[1409] = 4.0E0*I_ERI_D2z_F2xz_Px_S_bb-2.0E0*1*I_ERI_D2z_Pz_Px_S_b;
  abcd[1410] = 4.0E0*I_ERI_D2x_Fxyz_Px_S_bb;
  abcd[1411] = 4.0E0*I_ERI_Dxy_Fxyz_Px_S_bb;
  abcd[1412] = 4.0E0*I_ERI_Dxz_Fxyz_Px_S_bb;
  abcd[1413] = 4.0E0*I_ERI_D2y_Fxyz_Px_S_bb;
  abcd[1414] = 4.0E0*I_ERI_Dyz_Fxyz_Px_S_bb;
  abcd[1415] = 4.0E0*I_ERI_D2z_Fxyz_Px_S_bb;
  abcd[1416] = 4.0E0*I_ERI_D2x_Fx2z_Px_S_bb-2.0E0*1*I_ERI_D2x_Px_Px_S_b;
  abcd[1417] = 4.0E0*I_ERI_Dxy_Fx2z_Px_S_bb-2.0E0*1*I_ERI_Dxy_Px_Px_S_b;
  abcd[1418] = 4.0E0*I_ERI_Dxz_Fx2z_Px_S_bb-2.0E0*1*I_ERI_Dxz_Px_Px_S_b;
  abcd[1419] = 4.0E0*I_ERI_D2y_Fx2z_Px_S_bb-2.0E0*1*I_ERI_D2y_Px_Px_S_b;
  abcd[1420] = 4.0E0*I_ERI_Dyz_Fx2z_Px_S_bb-2.0E0*1*I_ERI_Dyz_Px_Px_S_b;
  abcd[1421] = 4.0E0*I_ERI_D2z_Fx2z_Px_S_bb-2.0E0*1*I_ERI_D2z_Px_Px_S_b;
  abcd[1422] = 4.0E0*I_ERI_D2x_F2xz_Py_S_bb-2.0E0*1*I_ERI_D2x_Pz_Py_S_b;
  abcd[1423] = 4.0E0*I_ERI_Dxy_F2xz_Py_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Py_S_b;
  abcd[1424] = 4.0E0*I_ERI_Dxz_F2xz_Py_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Py_S_b;
  abcd[1425] = 4.0E0*I_ERI_D2y_F2xz_Py_S_bb-2.0E0*1*I_ERI_D2y_Pz_Py_S_b;
  abcd[1426] = 4.0E0*I_ERI_Dyz_F2xz_Py_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Py_S_b;
  abcd[1427] = 4.0E0*I_ERI_D2z_F2xz_Py_S_bb-2.0E0*1*I_ERI_D2z_Pz_Py_S_b;
  abcd[1428] = 4.0E0*I_ERI_D2x_Fxyz_Py_S_bb;
  abcd[1429] = 4.0E0*I_ERI_Dxy_Fxyz_Py_S_bb;
  abcd[1430] = 4.0E0*I_ERI_Dxz_Fxyz_Py_S_bb;
  abcd[1431] = 4.0E0*I_ERI_D2y_Fxyz_Py_S_bb;
  abcd[1432] = 4.0E0*I_ERI_Dyz_Fxyz_Py_S_bb;
  abcd[1433] = 4.0E0*I_ERI_D2z_Fxyz_Py_S_bb;
  abcd[1434] = 4.0E0*I_ERI_D2x_Fx2z_Py_S_bb-2.0E0*1*I_ERI_D2x_Px_Py_S_b;
  abcd[1435] = 4.0E0*I_ERI_Dxy_Fx2z_Py_S_bb-2.0E0*1*I_ERI_Dxy_Px_Py_S_b;
  abcd[1436] = 4.0E0*I_ERI_Dxz_Fx2z_Py_S_bb-2.0E0*1*I_ERI_Dxz_Px_Py_S_b;
  abcd[1437] = 4.0E0*I_ERI_D2y_Fx2z_Py_S_bb-2.0E0*1*I_ERI_D2y_Px_Py_S_b;
  abcd[1438] = 4.0E0*I_ERI_Dyz_Fx2z_Py_S_bb-2.0E0*1*I_ERI_Dyz_Px_Py_S_b;
  abcd[1439] = 4.0E0*I_ERI_D2z_Fx2z_Py_S_bb-2.0E0*1*I_ERI_D2z_Px_Py_S_b;
  abcd[1440] = 4.0E0*I_ERI_D2x_F2xz_Pz_S_bb-2.0E0*1*I_ERI_D2x_Pz_Pz_S_b;
  abcd[1441] = 4.0E0*I_ERI_Dxy_F2xz_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_b;
  abcd[1442] = 4.0E0*I_ERI_Dxz_F2xz_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_b;
  abcd[1443] = 4.0E0*I_ERI_D2y_F2xz_Pz_S_bb-2.0E0*1*I_ERI_D2y_Pz_Pz_S_b;
  abcd[1444] = 4.0E0*I_ERI_Dyz_F2xz_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_b;
  abcd[1445] = 4.0E0*I_ERI_D2z_F2xz_Pz_S_bb-2.0E0*1*I_ERI_D2z_Pz_Pz_S_b;
  abcd[1446] = 4.0E0*I_ERI_D2x_Fxyz_Pz_S_bb;
  abcd[1447] = 4.0E0*I_ERI_Dxy_Fxyz_Pz_S_bb;
  abcd[1448] = 4.0E0*I_ERI_Dxz_Fxyz_Pz_S_bb;
  abcd[1449] = 4.0E0*I_ERI_D2y_Fxyz_Pz_S_bb;
  abcd[1450] = 4.0E0*I_ERI_Dyz_Fxyz_Pz_S_bb;
  abcd[1451] = 4.0E0*I_ERI_D2z_Fxyz_Pz_S_bb;
  abcd[1452] = 4.0E0*I_ERI_D2x_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_D2x_Px_Pz_S_b;
  abcd[1453] = 4.0E0*I_ERI_Dxy_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Px_Pz_S_b;
  abcd[1454] = 4.0E0*I_ERI_Dxz_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Px_Pz_S_b;
  abcd[1455] = 4.0E0*I_ERI_D2y_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_D2y_Px_Pz_S_b;
  abcd[1456] = 4.0E0*I_ERI_Dyz_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Px_Pz_S_b;
  abcd[1457] = 4.0E0*I_ERI_D2z_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_D2z_Px_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_P_S_bb
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   ************************************************************/
  abcd[1458] = 4.0E0*I_ERI_D2x_Fx2y_Px_S_bb-2.0E0*1*I_ERI_D2x_Px_Px_S_b;
  abcd[1459] = 4.0E0*I_ERI_Dxy_Fx2y_Px_S_bb-2.0E0*1*I_ERI_Dxy_Px_Px_S_b;
  abcd[1460] = 4.0E0*I_ERI_Dxz_Fx2y_Px_S_bb-2.0E0*1*I_ERI_Dxz_Px_Px_S_b;
  abcd[1461] = 4.0E0*I_ERI_D2y_Fx2y_Px_S_bb-2.0E0*1*I_ERI_D2y_Px_Px_S_b;
  abcd[1462] = 4.0E0*I_ERI_Dyz_Fx2y_Px_S_bb-2.0E0*1*I_ERI_Dyz_Px_Px_S_b;
  abcd[1463] = 4.0E0*I_ERI_D2z_Fx2y_Px_S_bb-2.0E0*1*I_ERI_D2z_Px_Px_S_b;
  abcd[1464] = 4.0E0*I_ERI_D2x_F3y_Px_S_bb-2.0E0*1*I_ERI_D2x_Py_Px_S_b-2.0E0*2*I_ERI_D2x_Py_Px_S_b;
  abcd[1465] = 4.0E0*I_ERI_Dxy_F3y_Px_S_bb-2.0E0*1*I_ERI_Dxy_Py_Px_S_b-2.0E0*2*I_ERI_Dxy_Py_Px_S_b;
  abcd[1466] = 4.0E0*I_ERI_Dxz_F3y_Px_S_bb-2.0E0*1*I_ERI_Dxz_Py_Px_S_b-2.0E0*2*I_ERI_Dxz_Py_Px_S_b;
  abcd[1467] = 4.0E0*I_ERI_D2y_F3y_Px_S_bb-2.0E0*1*I_ERI_D2y_Py_Px_S_b-2.0E0*2*I_ERI_D2y_Py_Px_S_b;
  abcd[1468] = 4.0E0*I_ERI_Dyz_F3y_Px_S_bb-2.0E0*1*I_ERI_Dyz_Py_Px_S_b-2.0E0*2*I_ERI_Dyz_Py_Px_S_b;
  abcd[1469] = 4.0E0*I_ERI_D2z_F3y_Px_S_bb-2.0E0*1*I_ERI_D2z_Py_Px_S_b-2.0E0*2*I_ERI_D2z_Py_Px_S_b;
  abcd[1470] = 4.0E0*I_ERI_D2x_F2yz_Px_S_bb-2.0E0*1*I_ERI_D2x_Pz_Px_S_b;
  abcd[1471] = 4.0E0*I_ERI_Dxy_F2yz_Px_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Px_S_b;
  abcd[1472] = 4.0E0*I_ERI_Dxz_F2yz_Px_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Px_S_b;
  abcd[1473] = 4.0E0*I_ERI_D2y_F2yz_Px_S_bb-2.0E0*1*I_ERI_D2y_Pz_Px_S_b;
  abcd[1474] = 4.0E0*I_ERI_Dyz_F2yz_Px_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Px_S_b;
  abcd[1475] = 4.0E0*I_ERI_D2z_F2yz_Px_S_bb-2.0E0*1*I_ERI_D2z_Pz_Px_S_b;
  abcd[1476] = 4.0E0*I_ERI_D2x_Fx2y_Py_S_bb-2.0E0*1*I_ERI_D2x_Px_Py_S_b;
  abcd[1477] = 4.0E0*I_ERI_Dxy_Fx2y_Py_S_bb-2.0E0*1*I_ERI_Dxy_Px_Py_S_b;
  abcd[1478] = 4.0E0*I_ERI_Dxz_Fx2y_Py_S_bb-2.0E0*1*I_ERI_Dxz_Px_Py_S_b;
  abcd[1479] = 4.0E0*I_ERI_D2y_Fx2y_Py_S_bb-2.0E0*1*I_ERI_D2y_Px_Py_S_b;
  abcd[1480] = 4.0E0*I_ERI_Dyz_Fx2y_Py_S_bb-2.0E0*1*I_ERI_Dyz_Px_Py_S_b;
  abcd[1481] = 4.0E0*I_ERI_D2z_Fx2y_Py_S_bb-2.0E0*1*I_ERI_D2z_Px_Py_S_b;
  abcd[1482] = 4.0E0*I_ERI_D2x_F3y_Py_S_bb-2.0E0*1*I_ERI_D2x_Py_Py_S_b-2.0E0*2*I_ERI_D2x_Py_Py_S_b;
  abcd[1483] = 4.0E0*I_ERI_Dxy_F3y_Py_S_bb-2.0E0*1*I_ERI_Dxy_Py_Py_S_b-2.0E0*2*I_ERI_Dxy_Py_Py_S_b;
  abcd[1484] = 4.0E0*I_ERI_Dxz_F3y_Py_S_bb-2.0E0*1*I_ERI_Dxz_Py_Py_S_b-2.0E0*2*I_ERI_Dxz_Py_Py_S_b;
  abcd[1485] = 4.0E0*I_ERI_D2y_F3y_Py_S_bb-2.0E0*1*I_ERI_D2y_Py_Py_S_b-2.0E0*2*I_ERI_D2y_Py_Py_S_b;
  abcd[1486] = 4.0E0*I_ERI_Dyz_F3y_Py_S_bb-2.0E0*1*I_ERI_Dyz_Py_Py_S_b-2.0E0*2*I_ERI_Dyz_Py_Py_S_b;
  abcd[1487] = 4.0E0*I_ERI_D2z_F3y_Py_S_bb-2.0E0*1*I_ERI_D2z_Py_Py_S_b-2.0E0*2*I_ERI_D2z_Py_Py_S_b;
  abcd[1488] = 4.0E0*I_ERI_D2x_F2yz_Py_S_bb-2.0E0*1*I_ERI_D2x_Pz_Py_S_b;
  abcd[1489] = 4.0E0*I_ERI_Dxy_F2yz_Py_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Py_S_b;
  abcd[1490] = 4.0E0*I_ERI_Dxz_F2yz_Py_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Py_S_b;
  abcd[1491] = 4.0E0*I_ERI_D2y_F2yz_Py_S_bb-2.0E0*1*I_ERI_D2y_Pz_Py_S_b;
  abcd[1492] = 4.0E0*I_ERI_Dyz_F2yz_Py_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Py_S_b;
  abcd[1493] = 4.0E0*I_ERI_D2z_F2yz_Py_S_bb-2.0E0*1*I_ERI_D2z_Pz_Py_S_b;
  abcd[1494] = 4.0E0*I_ERI_D2x_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_D2x_Px_Pz_S_b;
  abcd[1495] = 4.0E0*I_ERI_Dxy_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Px_Pz_S_b;
  abcd[1496] = 4.0E0*I_ERI_Dxz_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Px_Pz_S_b;
  abcd[1497] = 4.0E0*I_ERI_D2y_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_D2y_Px_Pz_S_b;
  abcd[1498] = 4.0E0*I_ERI_Dyz_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Px_Pz_S_b;
  abcd[1499] = 4.0E0*I_ERI_D2z_Fx2y_Pz_S_bb-2.0E0*1*I_ERI_D2z_Px_Pz_S_b;
  abcd[1500] = 4.0E0*I_ERI_D2x_F3y_Pz_S_bb-2.0E0*1*I_ERI_D2x_Py_Pz_S_b-2.0E0*2*I_ERI_D2x_Py_Pz_S_b;
  abcd[1501] = 4.0E0*I_ERI_Dxy_F3y_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Py_Pz_S_b-2.0E0*2*I_ERI_Dxy_Py_Pz_S_b;
  abcd[1502] = 4.0E0*I_ERI_Dxz_F3y_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Py_Pz_S_b-2.0E0*2*I_ERI_Dxz_Py_Pz_S_b;
  abcd[1503] = 4.0E0*I_ERI_D2y_F3y_Pz_S_bb-2.0E0*1*I_ERI_D2y_Py_Pz_S_b-2.0E0*2*I_ERI_D2y_Py_Pz_S_b;
  abcd[1504] = 4.0E0*I_ERI_Dyz_F3y_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Py_Pz_S_b-2.0E0*2*I_ERI_Dyz_Py_Pz_S_b;
  abcd[1505] = 4.0E0*I_ERI_D2z_F3y_Pz_S_bb-2.0E0*1*I_ERI_D2z_Py_Pz_S_b-2.0E0*2*I_ERI_D2z_Py_Pz_S_b;
  abcd[1506] = 4.0E0*I_ERI_D2x_F2yz_Pz_S_bb-2.0E0*1*I_ERI_D2x_Pz_Pz_S_b;
  abcd[1507] = 4.0E0*I_ERI_Dxy_F2yz_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_b;
  abcd[1508] = 4.0E0*I_ERI_Dxz_F2yz_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_b;
  abcd[1509] = 4.0E0*I_ERI_D2y_F2yz_Pz_S_bb-2.0E0*1*I_ERI_D2y_Pz_Pz_S_b;
  abcd[1510] = 4.0E0*I_ERI_Dyz_F2yz_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_b;
  abcd[1511] = 4.0E0*I_ERI_D2z_F2yz_Pz_S_bb-2.0E0*1*I_ERI_D2z_Pz_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_P_S_bb
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   ************************************************************/
  abcd[1512] = 4.0E0*I_ERI_D2x_Fxyz_Px_S_bb;
  abcd[1513] = 4.0E0*I_ERI_Dxy_Fxyz_Px_S_bb;
  abcd[1514] = 4.0E0*I_ERI_Dxz_Fxyz_Px_S_bb;
  abcd[1515] = 4.0E0*I_ERI_D2y_Fxyz_Px_S_bb;
  abcd[1516] = 4.0E0*I_ERI_Dyz_Fxyz_Px_S_bb;
  abcd[1517] = 4.0E0*I_ERI_D2z_Fxyz_Px_S_bb;
  abcd[1518] = 4.0E0*I_ERI_D2x_F2yz_Px_S_bb-2.0E0*1*I_ERI_D2x_Pz_Px_S_b;
  abcd[1519] = 4.0E0*I_ERI_Dxy_F2yz_Px_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Px_S_b;
  abcd[1520] = 4.0E0*I_ERI_Dxz_F2yz_Px_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Px_S_b;
  abcd[1521] = 4.0E0*I_ERI_D2y_F2yz_Px_S_bb-2.0E0*1*I_ERI_D2y_Pz_Px_S_b;
  abcd[1522] = 4.0E0*I_ERI_Dyz_F2yz_Px_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Px_S_b;
  abcd[1523] = 4.0E0*I_ERI_D2z_F2yz_Px_S_bb-2.0E0*1*I_ERI_D2z_Pz_Px_S_b;
  abcd[1524] = 4.0E0*I_ERI_D2x_Fy2z_Px_S_bb-2.0E0*1*I_ERI_D2x_Py_Px_S_b;
  abcd[1525] = 4.0E0*I_ERI_Dxy_Fy2z_Px_S_bb-2.0E0*1*I_ERI_Dxy_Py_Px_S_b;
  abcd[1526] = 4.0E0*I_ERI_Dxz_Fy2z_Px_S_bb-2.0E0*1*I_ERI_Dxz_Py_Px_S_b;
  abcd[1527] = 4.0E0*I_ERI_D2y_Fy2z_Px_S_bb-2.0E0*1*I_ERI_D2y_Py_Px_S_b;
  abcd[1528] = 4.0E0*I_ERI_Dyz_Fy2z_Px_S_bb-2.0E0*1*I_ERI_Dyz_Py_Px_S_b;
  abcd[1529] = 4.0E0*I_ERI_D2z_Fy2z_Px_S_bb-2.0E0*1*I_ERI_D2z_Py_Px_S_b;
  abcd[1530] = 4.0E0*I_ERI_D2x_Fxyz_Py_S_bb;
  abcd[1531] = 4.0E0*I_ERI_Dxy_Fxyz_Py_S_bb;
  abcd[1532] = 4.0E0*I_ERI_Dxz_Fxyz_Py_S_bb;
  abcd[1533] = 4.0E0*I_ERI_D2y_Fxyz_Py_S_bb;
  abcd[1534] = 4.0E0*I_ERI_Dyz_Fxyz_Py_S_bb;
  abcd[1535] = 4.0E0*I_ERI_D2z_Fxyz_Py_S_bb;
  abcd[1536] = 4.0E0*I_ERI_D2x_F2yz_Py_S_bb-2.0E0*1*I_ERI_D2x_Pz_Py_S_b;
  abcd[1537] = 4.0E0*I_ERI_Dxy_F2yz_Py_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Py_S_b;
  abcd[1538] = 4.0E0*I_ERI_Dxz_F2yz_Py_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Py_S_b;
  abcd[1539] = 4.0E0*I_ERI_D2y_F2yz_Py_S_bb-2.0E0*1*I_ERI_D2y_Pz_Py_S_b;
  abcd[1540] = 4.0E0*I_ERI_Dyz_F2yz_Py_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Py_S_b;
  abcd[1541] = 4.0E0*I_ERI_D2z_F2yz_Py_S_bb-2.0E0*1*I_ERI_D2z_Pz_Py_S_b;
  abcd[1542] = 4.0E0*I_ERI_D2x_Fy2z_Py_S_bb-2.0E0*1*I_ERI_D2x_Py_Py_S_b;
  abcd[1543] = 4.0E0*I_ERI_Dxy_Fy2z_Py_S_bb-2.0E0*1*I_ERI_Dxy_Py_Py_S_b;
  abcd[1544] = 4.0E0*I_ERI_Dxz_Fy2z_Py_S_bb-2.0E0*1*I_ERI_Dxz_Py_Py_S_b;
  abcd[1545] = 4.0E0*I_ERI_D2y_Fy2z_Py_S_bb-2.0E0*1*I_ERI_D2y_Py_Py_S_b;
  abcd[1546] = 4.0E0*I_ERI_Dyz_Fy2z_Py_S_bb-2.0E0*1*I_ERI_Dyz_Py_Py_S_b;
  abcd[1547] = 4.0E0*I_ERI_D2z_Fy2z_Py_S_bb-2.0E0*1*I_ERI_D2z_Py_Py_S_b;
  abcd[1548] = 4.0E0*I_ERI_D2x_Fxyz_Pz_S_bb;
  abcd[1549] = 4.0E0*I_ERI_Dxy_Fxyz_Pz_S_bb;
  abcd[1550] = 4.0E0*I_ERI_Dxz_Fxyz_Pz_S_bb;
  abcd[1551] = 4.0E0*I_ERI_D2y_Fxyz_Pz_S_bb;
  abcd[1552] = 4.0E0*I_ERI_Dyz_Fxyz_Pz_S_bb;
  abcd[1553] = 4.0E0*I_ERI_D2z_Fxyz_Pz_S_bb;
  abcd[1554] = 4.0E0*I_ERI_D2x_F2yz_Pz_S_bb-2.0E0*1*I_ERI_D2x_Pz_Pz_S_b;
  abcd[1555] = 4.0E0*I_ERI_Dxy_F2yz_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_b;
  abcd[1556] = 4.0E0*I_ERI_Dxz_F2yz_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_b;
  abcd[1557] = 4.0E0*I_ERI_D2y_F2yz_Pz_S_bb-2.0E0*1*I_ERI_D2y_Pz_Pz_S_b;
  abcd[1558] = 4.0E0*I_ERI_Dyz_F2yz_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_b;
  abcd[1559] = 4.0E0*I_ERI_D2z_F2yz_Pz_S_bb-2.0E0*1*I_ERI_D2z_Pz_Pz_S_b;
  abcd[1560] = 4.0E0*I_ERI_D2x_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_D2x_Py_Pz_S_b;
  abcd[1561] = 4.0E0*I_ERI_Dxy_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Py_Pz_S_b;
  abcd[1562] = 4.0E0*I_ERI_Dxz_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Py_Pz_S_b;
  abcd[1563] = 4.0E0*I_ERI_D2y_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_D2y_Py_Pz_S_b;
  abcd[1564] = 4.0E0*I_ERI_Dyz_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Py_Pz_S_b;
  abcd[1565] = 4.0E0*I_ERI_D2z_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_D2z_Py_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_F_P_S_bb
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   * RHS shell quartet name: SQ_ERI_D_P_P_S_b
   ************************************************************/
  abcd[1566] = 4.0E0*I_ERI_D2x_Fx2z_Px_S_bb-2.0E0*1*I_ERI_D2x_Px_Px_S_b;
  abcd[1567] = 4.0E0*I_ERI_Dxy_Fx2z_Px_S_bb-2.0E0*1*I_ERI_Dxy_Px_Px_S_b;
  abcd[1568] = 4.0E0*I_ERI_Dxz_Fx2z_Px_S_bb-2.0E0*1*I_ERI_Dxz_Px_Px_S_b;
  abcd[1569] = 4.0E0*I_ERI_D2y_Fx2z_Px_S_bb-2.0E0*1*I_ERI_D2y_Px_Px_S_b;
  abcd[1570] = 4.0E0*I_ERI_Dyz_Fx2z_Px_S_bb-2.0E0*1*I_ERI_Dyz_Px_Px_S_b;
  abcd[1571] = 4.0E0*I_ERI_D2z_Fx2z_Px_S_bb-2.0E0*1*I_ERI_D2z_Px_Px_S_b;
  abcd[1572] = 4.0E0*I_ERI_D2x_Fy2z_Px_S_bb-2.0E0*1*I_ERI_D2x_Py_Px_S_b;
  abcd[1573] = 4.0E0*I_ERI_Dxy_Fy2z_Px_S_bb-2.0E0*1*I_ERI_Dxy_Py_Px_S_b;
  abcd[1574] = 4.0E0*I_ERI_Dxz_Fy2z_Px_S_bb-2.0E0*1*I_ERI_Dxz_Py_Px_S_b;
  abcd[1575] = 4.0E0*I_ERI_D2y_Fy2z_Px_S_bb-2.0E0*1*I_ERI_D2y_Py_Px_S_b;
  abcd[1576] = 4.0E0*I_ERI_Dyz_Fy2z_Px_S_bb-2.0E0*1*I_ERI_Dyz_Py_Px_S_b;
  abcd[1577] = 4.0E0*I_ERI_D2z_Fy2z_Px_S_bb-2.0E0*1*I_ERI_D2z_Py_Px_S_b;
  abcd[1578] = 4.0E0*I_ERI_D2x_F3z_Px_S_bb-2.0E0*1*I_ERI_D2x_Pz_Px_S_b-2.0E0*2*I_ERI_D2x_Pz_Px_S_b;
  abcd[1579] = 4.0E0*I_ERI_Dxy_F3z_Px_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Px_S_b-2.0E0*2*I_ERI_Dxy_Pz_Px_S_b;
  abcd[1580] = 4.0E0*I_ERI_Dxz_F3z_Px_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Px_S_b-2.0E0*2*I_ERI_Dxz_Pz_Px_S_b;
  abcd[1581] = 4.0E0*I_ERI_D2y_F3z_Px_S_bb-2.0E0*1*I_ERI_D2y_Pz_Px_S_b-2.0E0*2*I_ERI_D2y_Pz_Px_S_b;
  abcd[1582] = 4.0E0*I_ERI_Dyz_F3z_Px_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Px_S_b-2.0E0*2*I_ERI_Dyz_Pz_Px_S_b;
  abcd[1583] = 4.0E0*I_ERI_D2z_F3z_Px_S_bb-2.0E0*1*I_ERI_D2z_Pz_Px_S_b-2.0E0*2*I_ERI_D2z_Pz_Px_S_b;
  abcd[1584] = 4.0E0*I_ERI_D2x_Fx2z_Py_S_bb-2.0E0*1*I_ERI_D2x_Px_Py_S_b;
  abcd[1585] = 4.0E0*I_ERI_Dxy_Fx2z_Py_S_bb-2.0E0*1*I_ERI_Dxy_Px_Py_S_b;
  abcd[1586] = 4.0E0*I_ERI_Dxz_Fx2z_Py_S_bb-2.0E0*1*I_ERI_Dxz_Px_Py_S_b;
  abcd[1587] = 4.0E0*I_ERI_D2y_Fx2z_Py_S_bb-2.0E0*1*I_ERI_D2y_Px_Py_S_b;
  abcd[1588] = 4.0E0*I_ERI_Dyz_Fx2z_Py_S_bb-2.0E0*1*I_ERI_Dyz_Px_Py_S_b;
  abcd[1589] = 4.0E0*I_ERI_D2z_Fx2z_Py_S_bb-2.0E0*1*I_ERI_D2z_Px_Py_S_b;
  abcd[1590] = 4.0E0*I_ERI_D2x_Fy2z_Py_S_bb-2.0E0*1*I_ERI_D2x_Py_Py_S_b;
  abcd[1591] = 4.0E0*I_ERI_Dxy_Fy2z_Py_S_bb-2.0E0*1*I_ERI_Dxy_Py_Py_S_b;
  abcd[1592] = 4.0E0*I_ERI_Dxz_Fy2z_Py_S_bb-2.0E0*1*I_ERI_Dxz_Py_Py_S_b;
  abcd[1593] = 4.0E0*I_ERI_D2y_Fy2z_Py_S_bb-2.0E0*1*I_ERI_D2y_Py_Py_S_b;
  abcd[1594] = 4.0E0*I_ERI_Dyz_Fy2z_Py_S_bb-2.0E0*1*I_ERI_Dyz_Py_Py_S_b;
  abcd[1595] = 4.0E0*I_ERI_D2z_Fy2z_Py_S_bb-2.0E0*1*I_ERI_D2z_Py_Py_S_b;
  abcd[1596] = 4.0E0*I_ERI_D2x_F3z_Py_S_bb-2.0E0*1*I_ERI_D2x_Pz_Py_S_b-2.0E0*2*I_ERI_D2x_Pz_Py_S_b;
  abcd[1597] = 4.0E0*I_ERI_Dxy_F3z_Py_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Py_S_b-2.0E0*2*I_ERI_Dxy_Pz_Py_S_b;
  abcd[1598] = 4.0E0*I_ERI_Dxz_F3z_Py_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Py_S_b-2.0E0*2*I_ERI_Dxz_Pz_Py_S_b;
  abcd[1599] = 4.0E0*I_ERI_D2y_F3z_Py_S_bb-2.0E0*1*I_ERI_D2y_Pz_Py_S_b-2.0E0*2*I_ERI_D2y_Pz_Py_S_b;
  abcd[1600] = 4.0E0*I_ERI_Dyz_F3z_Py_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Py_S_b-2.0E0*2*I_ERI_Dyz_Pz_Py_S_b;
  abcd[1601] = 4.0E0*I_ERI_D2z_F3z_Py_S_bb-2.0E0*1*I_ERI_D2z_Pz_Py_S_b-2.0E0*2*I_ERI_D2z_Pz_Py_S_b;
  abcd[1602] = 4.0E0*I_ERI_D2x_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_D2x_Px_Pz_S_b;
  abcd[1603] = 4.0E0*I_ERI_Dxy_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Px_Pz_S_b;
  abcd[1604] = 4.0E0*I_ERI_Dxz_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Px_Pz_S_b;
  abcd[1605] = 4.0E0*I_ERI_D2y_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_D2y_Px_Pz_S_b;
  abcd[1606] = 4.0E0*I_ERI_Dyz_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Px_Pz_S_b;
  abcd[1607] = 4.0E0*I_ERI_D2z_Fx2z_Pz_S_bb-2.0E0*1*I_ERI_D2z_Px_Pz_S_b;
  abcd[1608] = 4.0E0*I_ERI_D2x_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_D2x_Py_Pz_S_b;
  abcd[1609] = 4.0E0*I_ERI_Dxy_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Py_Pz_S_b;
  abcd[1610] = 4.0E0*I_ERI_Dxz_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Py_Pz_S_b;
  abcd[1611] = 4.0E0*I_ERI_D2y_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_D2y_Py_Pz_S_b;
  abcd[1612] = 4.0E0*I_ERI_Dyz_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Py_Pz_S_b;
  abcd[1613] = 4.0E0*I_ERI_D2z_Fy2z_Pz_S_bb-2.0E0*1*I_ERI_D2z_Py_Pz_S_b;
  abcd[1614] = 4.0E0*I_ERI_D2x_F3z_Pz_S_bb-2.0E0*1*I_ERI_D2x_Pz_Pz_S_b-2.0E0*2*I_ERI_D2x_Pz_Pz_S_b;
  abcd[1615] = 4.0E0*I_ERI_Dxy_F3z_Pz_S_bb-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_b-2.0E0*2*I_ERI_Dxy_Pz_Pz_S_b;
  abcd[1616] = 4.0E0*I_ERI_Dxz_F3z_Pz_S_bb-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_b-2.0E0*2*I_ERI_Dxz_Pz_Pz_S_b;
  abcd[1617] = 4.0E0*I_ERI_D2y_F3z_Pz_S_bb-2.0E0*1*I_ERI_D2y_Pz_Pz_S_b-2.0E0*2*I_ERI_D2y_Pz_Pz_S_b;
  abcd[1618] = 4.0E0*I_ERI_Dyz_F3z_Pz_S_bb-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_b-2.0E0*2*I_ERI_Dyz_Pz_Pz_S_b;
  abcd[1619] = 4.0E0*I_ERI_D2z_F3z_Pz_S_bb-2.0E0*1*I_ERI_D2z_Pz_Pz_S_b-2.0E0*2*I_ERI_D2z_Pz_Pz_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[1620] = 4.0E0*I_ERI_D2x_D2x_D2x_S_bc-2.0E0*1*I_ERI_D2x_D2x_S_S_b-2.0E0*1*I_ERI_D2x_S_D2x_S_c+1*I_ERI_D2x_S_S_S;
  abcd[1621] = 4.0E0*I_ERI_Dxy_D2x_D2x_S_bc-2.0E0*1*I_ERI_Dxy_D2x_S_S_b-2.0E0*1*I_ERI_Dxy_S_D2x_S_c+1*I_ERI_Dxy_S_S_S;
  abcd[1622] = 4.0E0*I_ERI_Dxz_D2x_D2x_S_bc-2.0E0*1*I_ERI_Dxz_D2x_S_S_b-2.0E0*1*I_ERI_Dxz_S_D2x_S_c+1*I_ERI_Dxz_S_S_S;
  abcd[1623] = 4.0E0*I_ERI_D2y_D2x_D2x_S_bc-2.0E0*1*I_ERI_D2y_D2x_S_S_b-2.0E0*1*I_ERI_D2y_S_D2x_S_c+1*I_ERI_D2y_S_S_S;
  abcd[1624] = 4.0E0*I_ERI_Dyz_D2x_D2x_S_bc-2.0E0*1*I_ERI_Dyz_D2x_S_S_b-2.0E0*1*I_ERI_Dyz_S_D2x_S_c+1*I_ERI_Dyz_S_S_S;
  abcd[1625] = 4.0E0*I_ERI_D2z_D2x_D2x_S_bc-2.0E0*1*I_ERI_D2z_D2x_S_S_b-2.0E0*1*I_ERI_D2z_S_D2x_S_c+1*I_ERI_D2z_S_S_S;
  abcd[1626] = 4.0E0*I_ERI_D2x_Dxy_D2x_S_bc-2.0E0*1*I_ERI_D2x_Dxy_S_S_b;
  abcd[1627] = 4.0E0*I_ERI_Dxy_Dxy_D2x_S_bc-2.0E0*1*I_ERI_Dxy_Dxy_S_S_b;
  abcd[1628] = 4.0E0*I_ERI_Dxz_Dxy_D2x_S_bc-2.0E0*1*I_ERI_Dxz_Dxy_S_S_b;
  abcd[1629] = 4.0E0*I_ERI_D2y_Dxy_D2x_S_bc-2.0E0*1*I_ERI_D2y_Dxy_S_S_b;
  abcd[1630] = 4.0E0*I_ERI_Dyz_Dxy_D2x_S_bc-2.0E0*1*I_ERI_Dyz_Dxy_S_S_b;
  abcd[1631] = 4.0E0*I_ERI_D2z_Dxy_D2x_S_bc-2.0E0*1*I_ERI_D2z_Dxy_S_S_b;
  abcd[1632] = 4.0E0*I_ERI_D2x_Dxz_D2x_S_bc-2.0E0*1*I_ERI_D2x_Dxz_S_S_b;
  abcd[1633] = 4.0E0*I_ERI_Dxy_Dxz_D2x_S_bc-2.0E0*1*I_ERI_Dxy_Dxz_S_S_b;
  abcd[1634] = 4.0E0*I_ERI_Dxz_Dxz_D2x_S_bc-2.0E0*1*I_ERI_Dxz_Dxz_S_S_b;
  abcd[1635] = 4.0E0*I_ERI_D2y_Dxz_D2x_S_bc-2.0E0*1*I_ERI_D2y_Dxz_S_S_b;
  abcd[1636] = 4.0E0*I_ERI_Dyz_Dxz_D2x_S_bc-2.0E0*1*I_ERI_Dyz_Dxz_S_S_b;
  abcd[1637] = 4.0E0*I_ERI_D2z_Dxz_D2x_S_bc-2.0E0*1*I_ERI_D2z_Dxz_S_S_b;
  abcd[1638] = 4.0E0*I_ERI_D2x_D2x_Dxy_S_bc-2.0E0*1*I_ERI_D2x_S_Dxy_S_c;
  abcd[1639] = 4.0E0*I_ERI_Dxy_D2x_Dxy_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxy_S_c;
  abcd[1640] = 4.0E0*I_ERI_Dxz_D2x_Dxy_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxy_S_c;
  abcd[1641] = 4.0E0*I_ERI_D2y_D2x_Dxy_S_bc-2.0E0*1*I_ERI_D2y_S_Dxy_S_c;
  abcd[1642] = 4.0E0*I_ERI_Dyz_D2x_Dxy_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxy_S_c;
  abcd[1643] = 4.0E0*I_ERI_D2z_D2x_Dxy_S_bc-2.0E0*1*I_ERI_D2z_S_Dxy_S_c;
  abcd[1644] = 4.0E0*I_ERI_D2x_Dxy_Dxy_S_bc;
  abcd[1645] = 4.0E0*I_ERI_Dxy_Dxy_Dxy_S_bc;
  abcd[1646] = 4.0E0*I_ERI_Dxz_Dxy_Dxy_S_bc;
  abcd[1647] = 4.0E0*I_ERI_D2y_Dxy_Dxy_S_bc;
  abcd[1648] = 4.0E0*I_ERI_Dyz_Dxy_Dxy_S_bc;
  abcd[1649] = 4.0E0*I_ERI_D2z_Dxy_Dxy_S_bc;
  abcd[1650] = 4.0E0*I_ERI_D2x_Dxz_Dxy_S_bc;
  abcd[1651] = 4.0E0*I_ERI_Dxy_Dxz_Dxy_S_bc;
  abcd[1652] = 4.0E0*I_ERI_Dxz_Dxz_Dxy_S_bc;
  abcd[1653] = 4.0E0*I_ERI_D2y_Dxz_Dxy_S_bc;
  abcd[1654] = 4.0E0*I_ERI_Dyz_Dxz_Dxy_S_bc;
  abcd[1655] = 4.0E0*I_ERI_D2z_Dxz_Dxy_S_bc;
  abcd[1656] = 4.0E0*I_ERI_D2x_D2x_Dxz_S_bc-2.0E0*1*I_ERI_D2x_S_Dxz_S_c;
  abcd[1657] = 4.0E0*I_ERI_Dxy_D2x_Dxz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxz_S_c;
  abcd[1658] = 4.0E0*I_ERI_Dxz_D2x_Dxz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxz_S_c;
  abcd[1659] = 4.0E0*I_ERI_D2y_D2x_Dxz_S_bc-2.0E0*1*I_ERI_D2y_S_Dxz_S_c;
  abcd[1660] = 4.0E0*I_ERI_Dyz_D2x_Dxz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxz_S_c;
  abcd[1661] = 4.0E0*I_ERI_D2z_D2x_Dxz_S_bc-2.0E0*1*I_ERI_D2z_S_Dxz_S_c;
  abcd[1662] = 4.0E0*I_ERI_D2x_Dxy_Dxz_S_bc;
  abcd[1663] = 4.0E0*I_ERI_Dxy_Dxy_Dxz_S_bc;
  abcd[1664] = 4.0E0*I_ERI_Dxz_Dxy_Dxz_S_bc;
  abcd[1665] = 4.0E0*I_ERI_D2y_Dxy_Dxz_S_bc;
  abcd[1666] = 4.0E0*I_ERI_Dyz_Dxy_Dxz_S_bc;
  abcd[1667] = 4.0E0*I_ERI_D2z_Dxy_Dxz_S_bc;
  abcd[1668] = 4.0E0*I_ERI_D2x_Dxz_Dxz_S_bc;
  abcd[1669] = 4.0E0*I_ERI_Dxy_Dxz_Dxz_S_bc;
  abcd[1670] = 4.0E0*I_ERI_Dxz_Dxz_Dxz_S_bc;
  abcd[1671] = 4.0E0*I_ERI_D2y_Dxz_Dxz_S_bc;
  abcd[1672] = 4.0E0*I_ERI_Dyz_Dxz_Dxz_S_bc;
  abcd[1673] = 4.0E0*I_ERI_D2z_Dxz_Dxz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[1674] = 4.0E0*I_ERI_D2x_D2x_Dxy_S_bc-2.0E0*1*I_ERI_D2x_S_Dxy_S_c;
  abcd[1675] = 4.0E0*I_ERI_Dxy_D2x_Dxy_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxy_S_c;
  abcd[1676] = 4.0E0*I_ERI_Dxz_D2x_Dxy_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxy_S_c;
  abcd[1677] = 4.0E0*I_ERI_D2y_D2x_Dxy_S_bc-2.0E0*1*I_ERI_D2y_S_Dxy_S_c;
  abcd[1678] = 4.0E0*I_ERI_Dyz_D2x_Dxy_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxy_S_c;
  abcd[1679] = 4.0E0*I_ERI_D2z_D2x_Dxy_S_bc-2.0E0*1*I_ERI_D2z_S_Dxy_S_c;
  abcd[1680] = 4.0E0*I_ERI_D2x_Dxy_Dxy_S_bc;
  abcd[1681] = 4.0E0*I_ERI_Dxy_Dxy_Dxy_S_bc;
  abcd[1682] = 4.0E0*I_ERI_Dxz_Dxy_Dxy_S_bc;
  abcd[1683] = 4.0E0*I_ERI_D2y_Dxy_Dxy_S_bc;
  abcd[1684] = 4.0E0*I_ERI_Dyz_Dxy_Dxy_S_bc;
  abcd[1685] = 4.0E0*I_ERI_D2z_Dxy_Dxy_S_bc;
  abcd[1686] = 4.0E0*I_ERI_D2x_Dxz_Dxy_S_bc;
  abcd[1687] = 4.0E0*I_ERI_Dxy_Dxz_Dxy_S_bc;
  abcd[1688] = 4.0E0*I_ERI_Dxz_Dxz_Dxy_S_bc;
  abcd[1689] = 4.0E0*I_ERI_D2y_Dxz_Dxy_S_bc;
  abcd[1690] = 4.0E0*I_ERI_Dyz_Dxz_Dxy_S_bc;
  abcd[1691] = 4.0E0*I_ERI_D2z_Dxz_Dxy_S_bc;
  abcd[1692] = 4.0E0*I_ERI_D2x_D2x_D2y_S_bc-2.0E0*1*I_ERI_D2x_D2x_S_S_b-2.0E0*1*I_ERI_D2x_S_D2y_S_c+1*I_ERI_D2x_S_S_S;
  abcd[1693] = 4.0E0*I_ERI_Dxy_D2x_D2y_S_bc-2.0E0*1*I_ERI_Dxy_D2x_S_S_b-2.0E0*1*I_ERI_Dxy_S_D2y_S_c+1*I_ERI_Dxy_S_S_S;
  abcd[1694] = 4.0E0*I_ERI_Dxz_D2x_D2y_S_bc-2.0E0*1*I_ERI_Dxz_D2x_S_S_b-2.0E0*1*I_ERI_Dxz_S_D2y_S_c+1*I_ERI_Dxz_S_S_S;
  abcd[1695] = 4.0E0*I_ERI_D2y_D2x_D2y_S_bc-2.0E0*1*I_ERI_D2y_D2x_S_S_b-2.0E0*1*I_ERI_D2y_S_D2y_S_c+1*I_ERI_D2y_S_S_S;
  abcd[1696] = 4.0E0*I_ERI_Dyz_D2x_D2y_S_bc-2.0E0*1*I_ERI_Dyz_D2x_S_S_b-2.0E0*1*I_ERI_Dyz_S_D2y_S_c+1*I_ERI_Dyz_S_S_S;
  abcd[1697] = 4.0E0*I_ERI_D2z_D2x_D2y_S_bc-2.0E0*1*I_ERI_D2z_D2x_S_S_b-2.0E0*1*I_ERI_D2z_S_D2y_S_c+1*I_ERI_D2z_S_S_S;
  abcd[1698] = 4.0E0*I_ERI_D2x_Dxy_D2y_S_bc-2.0E0*1*I_ERI_D2x_Dxy_S_S_b;
  abcd[1699] = 4.0E0*I_ERI_Dxy_Dxy_D2y_S_bc-2.0E0*1*I_ERI_Dxy_Dxy_S_S_b;
  abcd[1700] = 4.0E0*I_ERI_Dxz_Dxy_D2y_S_bc-2.0E0*1*I_ERI_Dxz_Dxy_S_S_b;
  abcd[1701] = 4.0E0*I_ERI_D2y_Dxy_D2y_S_bc-2.0E0*1*I_ERI_D2y_Dxy_S_S_b;
  abcd[1702] = 4.0E0*I_ERI_Dyz_Dxy_D2y_S_bc-2.0E0*1*I_ERI_Dyz_Dxy_S_S_b;
  abcd[1703] = 4.0E0*I_ERI_D2z_Dxy_D2y_S_bc-2.0E0*1*I_ERI_D2z_Dxy_S_S_b;
  abcd[1704] = 4.0E0*I_ERI_D2x_Dxz_D2y_S_bc-2.0E0*1*I_ERI_D2x_Dxz_S_S_b;
  abcd[1705] = 4.0E0*I_ERI_Dxy_Dxz_D2y_S_bc-2.0E0*1*I_ERI_Dxy_Dxz_S_S_b;
  abcd[1706] = 4.0E0*I_ERI_Dxz_Dxz_D2y_S_bc-2.0E0*1*I_ERI_Dxz_Dxz_S_S_b;
  abcd[1707] = 4.0E0*I_ERI_D2y_Dxz_D2y_S_bc-2.0E0*1*I_ERI_D2y_Dxz_S_S_b;
  abcd[1708] = 4.0E0*I_ERI_Dyz_Dxz_D2y_S_bc-2.0E0*1*I_ERI_Dyz_Dxz_S_S_b;
  abcd[1709] = 4.0E0*I_ERI_D2z_Dxz_D2y_S_bc-2.0E0*1*I_ERI_D2z_Dxz_S_S_b;
  abcd[1710] = 4.0E0*I_ERI_D2x_D2x_Dyz_S_bc-2.0E0*1*I_ERI_D2x_S_Dyz_S_c;
  abcd[1711] = 4.0E0*I_ERI_Dxy_D2x_Dyz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dyz_S_c;
  abcd[1712] = 4.0E0*I_ERI_Dxz_D2x_Dyz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dyz_S_c;
  abcd[1713] = 4.0E0*I_ERI_D2y_D2x_Dyz_S_bc-2.0E0*1*I_ERI_D2y_S_Dyz_S_c;
  abcd[1714] = 4.0E0*I_ERI_Dyz_D2x_Dyz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dyz_S_c;
  abcd[1715] = 4.0E0*I_ERI_D2z_D2x_Dyz_S_bc-2.0E0*1*I_ERI_D2z_S_Dyz_S_c;
  abcd[1716] = 4.0E0*I_ERI_D2x_Dxy_Dyz_S_bc;
  abcd[1717] = 4.0E0*I_ERI_Dxy_Dxy_Dyz_S_bc;
  abcd[1718] = 4.0E0*I_ERI_Dxz_Dxy_Dyz_S_bc;
  abcd[1719] = 4.0E0*I_ERI_D2y_Dxy_Dyz_S_bc;
  abcd[1720] = 4.0E0*I_ERI_Dyz_Dxy_Dyz_S_bc;
  abcd[1721] = 4.0E0*I_ERI_D2z_Dxy_Dyz_S_bc;
  abcd[1722] = 4.0E0*I_ERI_D2x_Dxz_Dyz_S_bc;
  abcd[1723] = 4.0E0*I_ERI_Dxy_Dxz_Dyz_S_bc;
  abcd[1724] = 4.0E0*I_ERI_Dxz_Dxz_Dyz_S_bc;
  abcd[1725] = 4.0E0*I_ERI_D2y_Dxz_Dyz_S_bc;
  abcd[1726] = 4.0E0*I_ERI_Dyz_Dxz_Dyz_S_bc;
  abcd[1727] = 4.0E0*I_ERI_D2z_Dxz_Dyz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[1728] = 4.0E0*I_ERI_D2x_D2x_Dxz_S_bc-2.0E0*1*I_ERI_D2x_S_Dxz_S_c;
  abcd[1729] = 4.0E0*I_ERI_Dxy_D2x_Dxz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxz_S_c;
  abcd[1730] = 4.0E0*I_ERI_Dxz_D2x_Dxz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxz_S_c;
  abcd[1731] = 4.0E0*I_ERI_D2y_D2x_Dxz_S_bc-2.0E0*1*I_ERI_D2y_S_Dxz_S_c;
  abcd[1732] = 4.0E0*I_ERI_Dyz_D2x_Dxz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxz_S_c;
  abcd[1733] = 4.0E0*I_ERI_D2z_D2x_Dxz_S_bc-2.0E0*1*I_ERI_D2z_S_Dxz_S_c;
  abcd[1734] = 4.0E0*I_ERI_D2x_Dxy_Dxz_S_bc;
  abcd[1735] = 4.0E0*I_ERI_Dxy_Dxy_Dxz_S_bc;
  abcd[1736] = 4.0E0*I_ERI_Dxz_Dxy_Dxz_S_bc;
  abcd[1737] = 4.0E0*I_ERI_D2y_Dxy_Dxz_S_bc;
  abcd[1738] = 4.0E0*I_ERI_Dyz_Dxy_Dxz_S_bc;
  abcd[1739] = 4.0E0*I_ERI_D2z_Dxy_Dxz_S_bc;
  abcd[1740] = 4.0E0*I_ERI_D2x_Dxz_Dxz_S_bc;
  abcd[1741] = 4.0E0*I_ERI_Dxy_Dxz_Dxz_S_bc;
  abcd[1742] = 4.0E0*I_ERI_Dxz_Dxz_Dxz_S_bc;
  abcd[1743] = 4.0E0*I_ERI_D2y_Dxz_Dxz_S_bc;
  abcd[1744] = 4.0E0*I_ERI_Dyz_Dxz_Dxz_S_bc;
  abcd[1745] = 4.0E0*I_ERI_D2z_Dxz_Dxz_S_bc;
  abcd[1746] = 4.0E0*I_ERI_D2x_D2x_Dyz_S_bc-2.0E0*1*I_ERI_D2x_S_Dyz_S_c;
  abcd[1747] = 4.0E0*I_ERI_Dxy_D2x_Dyz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dyz_S_c;
  abcd[1748] = 4.0E0*I_ERI_Dxz_D2x_Dyz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dyz_S_c;
  abcd[1749] = 4.0E0*I_ERI_D2y_D2x_Dyz_S_bc-2.0E0*1*I_ERI_D2y_S_Dyz_S_c;
  abcd[1750] = 4.0E0*I_ERI_Dyz_D2x_Dyz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dyz_S_c;
  abcd[1751] = 4.0E0*I_ERI_D2z_D2x_Dyz_S_bc-2.0E0*1*I_ERI_D2z_S_Dyz_S_c;
  abcd[1752] = 4.0E0*I_ERI_D2x_Dxy_Dyz_S_bc;
  abcd[1753] = 4.0E0*I_ERI_Dxy_Dxy_Dyz_S_bc;
  abcd[1754] = 4.0E0*I_ERI_Dxz_Dxy_Dyz_S_bc;
  abcd[1755] = 4.0E0*I_ERI_D2y_Dxy_Dyz_S_bc;
  abcd[1756] = 4.0E0*I_ERI_Dyz_Dxy_Dyz_S_bc;
  abcd[1757] = 4.0E0*I_ERI_D2z_Dxy_Dyz_S_bc;
  abcd[1758] = 4.0E0*I_ERI_D2x_Dxz_Dyz_S_bc;
  abcd[1759] = 4.0E0*I_ERI_Dxy_Dxz_Dyz_S_bc;
  abcd[1760] = 4.0E0*I_ERI_Dxz_Dxz_Dyz_S_bc;
  abcd[1761] = 4.0E0*I_ERI_D2y_Dxz_Dyz_S_bc;
  abcd[1762] = 4.0E0*I_ERI_Dyz_Dxz_Dyz_S_bc;
  abcd[1763] = 4.0E0*I_ERI_D2z_Dxz_Dyz_S_bc;
  abcd[1764] = 4.0E0*I_ERI_D2x_D2x_D2z_S_bc-2.0E0*1*I_ERI_D2x_D2x_S_S_b-2.0E0*1*I_ERI_D2x_S_D2z_S_c+1*I_ERI_D2x_S_S_S;
  abcd[1765] = 4.0E0*I_ERI_Dxy_D2x_D2z_S_bc-2.0E0*1*I_ERI_Dxy_D2x_S_S_b-2.0E0*1*I_ERI_Dxy_S_D2z_S_c+1*I_ERI_Dxy_S_S_S;
  abcd[1766] = 4.0E0*I_ERI_Dxz_D2x_D2z_S_bc-2.0E0*1*I_ERI_Dxz_D2x_S_S_b-2.0E0*1*I_ERI_Dxz_S_D2z_S_c+1*I_ERI_Dxz_S_S_S;
  abcd[1767] = 4.0E0*I_ERI_D2y_D2x_D2z_S_bc-2.0E0*1*I_ERI_D2y_D2x_S_S_b-2.0E0*1*I_ERI_D2y_S_D2z_S_c+1*I_ERI_D2y_S_S_S;
  abcd[1768] = 4.0E0*I_ERI_Dyz_D2x_D2z_S_bc-2.0E0*1*I_ERI_Dyz_D2x_S_S_b-2.0E0*1*I_ERI_Dyz_S_D2z_S_c+1*I_ERI_Dyz_S_S_S;
  abcd[1769] = 4.0E0*I_ERI_D2z_D2x_D2z_S_bc-2.0E0*1*I_ERI_D2z_D2x_S_S_b-2.0E0*1*I_ERI_D2z_S_D2z_S_c+1*I_ERI_D2z_S_S_S;
  abcd[1770] = 4.0E0*I_ERI_D2x_Dxy_D2z_S_bc-2.0E0*1*I_ERI_D2x_Dxy_S_S_b;
  abcd[1771] = 4.0E0*I_ERI_Dxy_Dxy_D2z_S_bc-2.0E0*1*I_ERI_Dxy_Dxy_S_S_b;
  abcd[1772] = 4.0E0*I_ERI_Dxz_Dxy_D2z_S_bc-2.0E0*1*I_ERI_Dxz_Dxy_S_S_b;
  abcd[1773] = 4.0E0*I_ERI_D2y_Dxy_D2z_S_bc-2.0E0*1*I_ERI_D2y_Dxy_S_S_b;
  abcd[1774] = 4.0E0*I_ERI_Dyz_Dxy_D2z_S_bc-2.0E0*1*I_ERI_Dyz_Dxy_S_S_b;
  abcd[1775] = 4.0E0*I_ERI_D2z_Dxy_D2z_S_bc-2.0E0*1*I_ERI_D2z_Dxy_S_S_b;
  abcd[1776] = 4.0E0*I_ERI_D2x_Dxz_D2z_S_bc-2.0E0*1*I_ERI_D2x_Dxz_S_S_b;
  abcd[1777] = 4.0E0*I_ERI_Dxy_Dxz_D2z_S_bc-2.0E0*1*I_ERI_Dxy_Dxz_S_S_b;
  abcd[1778] = 4.0E0*I_ERI_Dxz_Dxz_D2z_S_bc-2.0E0*1*I_ERI_Dxz_Dxz_S_S_b;
  abcd[1779] = 4.0E0*I_ERI_D2y_Dxz_D2z_S_bc-2.0E0*1*I_ERI_D2y_Dxz_S_S_b;
  abcd[1780] = 4.0E0*I_ERI_Dyz_Dxz_D2z_S_bc-2.0E0*1*I_ERI_Dyz_Dxz_S_S_b;
  abcd[1781] = 4.0E0*I_ERI_D2z_Dxz_D2z_S_bc-2.0E0*1*I_ERI_D2z_Dxz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[1782] = 4.0E0*I_ERI_D2x_Dxy_D2x_S_bc-2.0E0*1*I_ERI_D2x_Dxy_S_S_b;
  abcd[1783] = 4.0E0*I_ERI_Dxy_Dxy_D2x_S_bc-2.0E0*1*I_ERI_Dxy_Dxy_S_S_b;
  abcd[1784] = 4.0E0*I_ERI_Dxz_Dxy_D2x_S_bc-2.0E0*1*I_ERI_Dxz_Dxy_S_S_b;
  abcd[1785] = 4.0E0*I_ERI_D2y_Dxy_D2x_S_bc-2.0E0*1*I_ERI_D2y_Dxy_S_S_b;
  abcd[1786] = 4.0E0*I_ERI_Dyz_Dxy_D2x_S_bc-2.0E0*1*I_ERI_Dyz_Dxy_S_S_b;
  abcd[1787] = 4.0E0*I_ERI_D2z_Dxy_D2x_S_bc-2.0E0*1*I_ERI_D2z_Dxy_S_S_b;
  abcd[1788] = 4.0E0*I_ERI_D2x_D2y_D2x_S_bc-2.0E0*1*I_ERI_D2x_D2y_S_S_b-2.0E0*1*I_ERI_D2x_S_D2x_S_c+1*I_ERI_D2x_S_S_S;
  abcd[1789] = 4.0E0*I_ERI_Dxy_D2y_D2x_S_bc-2.0E0*1*I_ERI_Dxy_D2y_S_S_b-2.0E0*1*I_ERI_Dxy_S_D2x_S_c+1*I_ERI_Dxy_S_S_S;
  abcd[1790] = 4.0E0*I_ERI_Dxz_D2y_D2x_S_bc-2.0E0*1*I_ERI_Dxz_D2y_S_S_b-2.0E0*1*I_ERI_Dxz_S_D2x_S_c+1*I_ERI_Dxz_S_S_S;
  abcd[1791] = 4.0E0*I_ERI_D2y_D2y_D2x_S_bc-2.0E0*1*I_ERI_D2y_D2y_S_S_b-2.0E0*1*I_ERI_D2y_S_D2x_S_c+1*I_ERI_D2y_S_S_S;
  abcd[1792] = 4.0E0*I_ERI_Dyz_D2y_D2x_S_bc-2.0E0*1*I_ERI_Dyz_D2y_S_S_b-2.0E0*1*I_ERI_Dyz_S_D2x_S_c+1*I_ERI_Dyz_S_S_S;
  abcd[1793] = 4.0E0*I_ERI_D2z_D2y_D2x_S_bc-2.0E0*1*I_ERI_D2z_D2y_S_S_b-2.0E0*1*I_ERI_D2z_S_D2x_S_c+1*I_ERI_D2z_S_S_S;
  abcd[1794] = 4.0E0*I_ERI_D2x_Dyz_D2x_S_bc-2.0E0*1*I_ERI_D2x_Dyz_S_S_b;
  abcd[1795] = 4.0E0*I_ERI_Dxy_Dyz_D2x_S_bc-2.0E0*1*I_ERI_Dxy_Dyz_S_S_b;
  abcd[1796] = 4.0E0*I_ERI_Dxz_Dyz_D2x_S_bc-2.0E0*1*I_ERI_Dxz_Dyz_S_S_b;
  abcd[1797] = 4.0E0*I_ERI_D2y_Dyz_D2x_S_bc-2.0E0*1*I_ERI_D2y_Dyz_S_S_b;
  abcd[1798] = 4.0E0*I_ERI_Dyz_Dyz_D2x_S_bc-2.0E0*1*I_ERI_Dyz_Dyz_S_S_b;
  abcd[1799] = 4.0E0*I_ERI_D2z_Dyz_D2x_S_bc-2.0E0*1*I_ERI_D2z_Dyz_S_S_b;
  abcd[1800] = 4.0E0*I_ERI_D2x_Dxy_Dxy_S_bc;
  abcd[1801] = 4.0E0*I_ERI_Dxy_Dxy_Dxy_S_bc;
  abcd[1802] = 4.0E0*I_ERI_Dxz_Dxy_Dxy_S_bc;
  abcd[1803] = 4.0E0*I_ERI_D2y_Dxy_Dxy_S_bc;
  abcd[1804] = 4.0E0*I_ERI_Dyz_Dxy_Dxy_S_bc;
  abcd[1805] = 4.0E0*I_ERI_D2z_Dxy_Dxy_S_bc;
  abcd[1806] = 4.0E0*I_ERI_D2x_D2y_Dxy_S_bc-2.0E0*1*I_ERI_D2x_S_Dxy_S_c;
  abcd[1807] = 4.0E0*I_ERI_Dxy_D2y_Dxy_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxy_S_c;
  abcd[1808] = 4.0E0*I_ERI_Dxz_D2y_Dxy_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxy_S_c;
  abcd[1809] = 4.0E0*I_ERI_D2y_D2y_Dxy_S_bc-2.0E0*1*I_ERI_D2y_S_Dxy_S_c;
  abcd[1810] = 4.0E0*I_ERI_Dyz_D2y_Dxy_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxy_S_c;
  abcd[1811] = 4.0E0*I_ERI_D2z_D2y_Dxy_S_bc-2.0E0*1*I_ERI_D2z_S_Dxy_S_c;
  abcd[1812] = 4.0E0*I_ERI_D2x_Dyz_Dxy_S_bc;
  abcd[1813] = 4.0E0*I_ERI_Dxy_Dyz_Dxy_S_bc;
  abcd[1814] = 4.0E0*I_ERI_Dxz_Dyz_Dxy_S_bc;
  abcd[1815] = 4.0E0*I_ERI_D2y_Dyz_Dxy_S_bc;
  abcd[1816] = 4.0E0*I_ERI_Dyz_Dyz_Dxy_S_bc;
  abcd[1817] = 4.0E0*I_ERI_D2z_Dyz_Dxy_S_bc;
  abcd[1818] = 4.0E0*I_ERI_D2x_Dxy_Dxz_S_bc;
  abcd[1819] = 4.0E0*I_ERI_Dxy_Dxy_Dxz_S_bc;
  abcd[1820] = 4.0E0*I_ERI_Dxz_Dxy_Dxz_S_bc;
  abcd[1821] = 4.0E0*I_ERI_D2y_Dxy_Dxz_S_bc;
  abcd[1822] = 4.0E0*I_ERI_Dyz_Dxy_Dxz_S_bc;
  abcd[1823] = 4.0E0*I_ERI_D2z_Dxy_Dxz_S_bc;
  abcd[1824] = 4.0E0*I_ERI_D2x_D2y_Dxz_S_bc-2.0E0*1*I_ERI_D2x_S_Dxz_S_c;
  abcd[1825] = 4.0E0*I_ERI_Dxy_D2y_Dxz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxz_S_c;
  abcd[1826] = 4.0E0*I_ERI_Dxz_D2y_Dxz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxz_S_c;
  abcd[1827] = 4.0E0*I_ERI_D2y_D2y_Dxz_S_bc-2.0E0*1*I_ERI_D2y_S_Dxz_S_c;
  abcd[1828] = 4.0E0*I_ERI_Dyz_D2y_Dxz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxz_S_c;
  abcd[1829] = 4.0E0*I_ERI_D2z_D2y_Dxz_S_bc-2.0E0*1*I_ERI_D2z_S_Dxz_S_c;
  abcd[1830] = 4.0E0*I_ERI_D2x_Dyz_Dxz_S_bc;
  abcd[1831] = 4.0E0*I_ERI_Dxy_Dyz_Dxz_S_bc;
  abcd[1832] = 4.0E0*I_ERI_Dxz_Dyz_Dxz_S_bc;
  abcd[1833] = 4.0E0*I_ERI_D2y_Dyz_Dxz_S_bc;
  abcd[1834] = 4.0E0*I_ERI_Dyz_Dyz_Dxz_S_bc;
  abcd[1835] = 4.0E0*I_ERI_D2z_Dyz_Dxz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[1836] = 4.0E0*I_ERI_D2x_Dxy_Dxy_S_bc;
  abcd[1837] = 4.0E0*I_ERI_Dxy_Dxy_Dxy_S_bc;
  abcd[1838] = 4.0E0*I_ERI_Dxz_Dxy_Dxy_S_bc;
  abcd[1839] = 4.0E0*I_ERI_D2y_Dxy_Dxy_S_bc;
  abcd[1840] = 4.0E0*I_ERI_Dyz_Dxy_Dxy_S_bc;
  abcd[1841] = 4.0E0*I_ERI_D2z_Dxy_Dxy_S_bc;
  abcd[1842] = 4.0E0*I_ERI_D2x_D2y_Dxy_S_bc-2.0E0*1*I_ERI_D2x_S_Dxy_S_c;
  abcd[1843] = 4.0E0*I_ERI_Dxy_D2y_Dxy_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxy_S_c;
  abcd[1844] = 4.0E0*I_ERI_Dxz_D2y_Dxy_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxy_S_c;
  abcd[1845] = 4.0E0*I_ERI_D2y_D2y_Dxy_S_bc-2.0E0*1*I_ERI_D2y_S_Dxy_S_c;
  abcd[1846] = 4.0E0*I_ERI_Dyz_D2y_Dxy_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxy_S_c;
  abcd[1847] = 4.0E0*I_ERI_D2z_D2y_Dxy_S_bc-2.0E0*1*I_ERI_D2z_S_Dxy_S_c;
  abcd[1848] = 4.0E0*I_ERI_D2x_Dyz_Dxy_S_bc;
  abcd[1849] = 4.0E0*I_ERI_Dxy_Dyz_Dxy_S_bc;
  abcd[1850] = 4.0E0*I_ERI_Dxz_Dyz_Dxy_S_bc;
  abcd[1851] = 4.0E0*I_ERI_D2y_Dyz_Dxy_S_bc;
  abcd[1852] = 4.0E0*I_ERI_Dyz_Dyz_Dxy_S_bc;
  abcd[1853] = 4.0E0*I_ERI_D2z_Dyz_Dxy_S_bc;
  abcd[1854] = 4.0E0*I_ERI_D2x_Dxy_D2y_S_bc-2.0E0*1*I_ERI_D2x_Dxy_S_S_b;
  abcd[1855] = 4.0E0*I_ERI_Dxy_Dxy_D2y_S_bc-2.0E0*1*I_ERI_Dxy_Dxy_S_S_b;
  abcd[1856] = 4.0E0*I_ERI_Dxz_Dxy_D2y_S_bc-2.0E0*1*I_ERI_Dxz_Dxy_S_S_b;
  abcd[1857] = 4.0E0*I_ERI_D2y_Dxy_D2y_S_bc-2.0E0*1*I_ERI_D2y_Dxy_S_S_b;
  abcd[1858] = 4.0E0*I_ERI_Dyz_Dxy_D2y_S_bc-2.0E0*1*I_ERI_Dyz_Dxy_S_S_b;
  abcd[1859] = 4.0E0*I_ERI_D2z_Dxy_D2y_S_bc-2.0E0*1*I_ERI_D2z_Dxy_S_S_b;
  abcd[1860] = 4.0E0*I_ERI_D2x_D2y_D2y_S_bc-2.0E0*1*I_ERI_D2x_D2y_S_S_b-2.0E0*1*I_ERI_D2x_S_D2y_S_c+1*I_ERI_D2x_S_S_S;
  abcd[1861] = 4.0E0*I_ERI_Dxy_D2y_D2y_S_bc-2.0E0*1*I_ERI_Dxy_D2y_S_S_b-2.0E0*1*I_ERI_Dxy_S_D2y_S_c+1*I_ERI_Dxy_S_S_S;
  abcd[1862] = 4.0E0*I_ERI_Dxz_D2y_D2y_S_bc-2.0E0*1*I_ERI_Dxz_D2y_S_S_b-2.0E0*1*I_ERI_Dxz_S_D2y_S_c+1*I_ERI_Dxz_S_S_S;
  abcd[1863] = 4.0E0*I_ERI_D2y_D2y_D2y_S_bc-2.0E0*1*I_ERI_D2y_D2y_S_S_b-2.0E0*1*I_ERI_D2y_S_D2y_S_c+1*I_ERI_D2y_S_S_S;
  abcd[1864] = 4.0E0*I_ERI_Dyz_D2y_D2y_S_bc-2.0E0*1*I_ERI_Dyz_D2y_S_S_b-2.0E0*1*I_ERI_Dyz_S_D2y_S_c+1*I_ERI_Dyz_S_S_S;
  abcd[1865] = 4.0E0*I_ERI_D2z_D2y_D2y_S_bc-2.0E0*1*I_ERI_D2z_D2y_S_S_b-2.0E0*1*I_ERI_D2z_S_D2y_S_c+1*I_ERI_D2z_S_S_S;
  abcd[1866] = 4.0E0*I_ERI_D2x_Dyz_D2y_S_bc-2.0E0*1*I_ERI_D2x_Dyz_S_S_b;
  abcd[1867] = 4.0E0*I_ERI_Dxy_Dyz_D2y_S_bc-2.0E0*1*I_ERI_Dxy_Dyz_S_S_b;
  abcd[1868] = 4.0E0*I_ERI_Dxz_Dyz_D2y_S_bc-2.0E0*1*I_ERI_Dxz_Dyz_S_S_b;
  abcd[1869] = 4.0E0*I_ERI_D2y_Dyz_D2y_S_bc-2.0E0*1*I_ERI_D2y_Dyz_S_S_b;
  abcd[1870] = 4.0E0*I_ERI_Dyz_Dyz_D2y_S_bc-2.0E0*1*I_ERI_Dyz_Dyz_S_S_b;
  abcd[1871] = 4.0E0*I_ERI_D2z_Dyz_D2y_S_bc-2.0E0*1*I_ERI_D2z_Dyz_S_S_b;
  abcd[1872] = 4.0E0*I_ERI_D2x_Dxy_Dyz_S_bc;
  abcd[1873] = 4.0E0*I_ERI_Dxy_Dxy_Dyz_S_bc;
  abcd[1874] = 4.0E0*I_ERI_Dxz_Dxy_Dyz_S_bc;
  abcd[1875] = 4.0E0*I_ERI_D2y_Dxy_Dyz_S_bc;
  abcd[1876] = 4.0E0*I_ERI_Dyz_Dxy_Dyz_S_bc;
  abcd[1877] = 4.0E0*I_ERI_D2z_Dxy_Dyz_S_bc;
  abcd[1878] = 4.0E0*I_ERI_D2x_D2y_Dyz_S_bc-2.0E0*1*I_ERI_D2x_S_Dyz_S_c;
  abcd[1879] = 4.0E0*I_ERI_Dxy_D2y_Dyz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dyz_S_c;
  abcd[1880] = 4.0E0*I_ERI_Dxz_D2y_Dyz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dyz_S_c;
  abcd[1881] = 4.0E0*I_ERI_D2y_D2y_Dyz_S_bc-2.0E0*1*I_ERI_D2y_S_Dyz_S_c;
  abcd[1882] = 4.0E0*I_ERI_Dyz_D2y_Dyz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dyz_S_c;
  abcd[1883] = 4.0E0*I_ERI_D2z_D2y_Dyz_S_bc-2.0E0*1*I_ERI_D2z_S_Dyz_S_c;
  abcd[1884] = 4.0E0*I_ERI_D2x_Dyz_Dyz_S_bc;
  abcd[1885] = 4.0E0*I_ERI_Dxy_Dyz_Dyz_S_bc;
  abcd[1886] = 4.0E0*I_ERI_Dxz_Dyz_Dyz_S_bc;
  abcd[1887] = 4.0E0*I_ERI_D2y_Dyz_Dyz_S_bc;
  abcd[1888] = 4.0E0*I_ERI_Dyz_Dyz_Dyz_S_bc;
  abcd[1889] = 4.0E0*I_ERI_D2z_Dyz_Dyz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[1890] = 4.0E0*I_ERI_D2x_Dxy_Dxz_S_bc;
  abcd[1891] = 4.0E0*I_ERI_Dxy_Dxy_Dxz_S_bc;
  abcd[1892] = 4.0E0*I_ERI_Dxz_Dxy_Dxz_S_bc;
  abcd[1893] = 4.0E0*I_ERI_D2y_Dxy_Dxz_S_bc;
  abcd[1894] = 4.0E0*I_ERI_Dyz_Dxy_Dxz_S_bc;
  abcd[1895] = 4.0E0*I_ERI_D2z_Dxy_Dxz_S_bc;
  abcd[1896] = 4.0E0*I_ERI_D2x_D2y_Dxz_S_bc-2.0E0*1*I_ERI_D2x_S_Dxz_S_c;
  abcd[1897] = 4.0E0*I_ERI_Dxy_D2y_Dxz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxz_S_c;
  abcd[1898] = 4.0E0*I_ERI_Dxz_D2y_Dxz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxz_S_c;
  abcd[1899] = 4.0E0*I_ERI_D2y_D2y_Dxz_S_bc-2.0E0*1*I_ERI_D2y_S_Dxz_S_c;
  abcd[1900] = 4.0E0*I_ERI_Dyz_D2y_Dxz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxz_S_c;
  abcd[1901] = 4.0E0*I_ERI_D2z_D2y_Dxz_S_bc-2.0E0*1*I_ERI_D2z_S_Dxz_S_c;
  abcd[1902] = 4.0E0*I_ERI_D2x_Dyz_Dxz_S_bc;
  abcd[1903] = 4.0E0*I_ERI_Dxy_Dyz_Dxz_S_bc;
  abcd[1904] = 4.0E0*I_ERI_Dxz_Dyz_Dxz_S_bc;
  abcd[1905] = 4.0E0*I_ERI_D2y_Dyz_Dxz_S_bc;
  abcd[1906] = 4.0E0*I_ERI_Dyz_Dyz_Dxz_S_bc;
  abcd[1907] = 4.0E0*I_ERI_D2z_Dyz_Dxz_S_bc;
  abcd[1908] = 4.0E0*I_ERI_D2x_Dxy_Dyz_S_bc;
  abcd[1909] = 4.0E0*I_ERI_Dxy_Dxy_Dyz_S_bc;
  abcd[1910] = 4.0E0*I_ERI_Dxz_Dxy_Dyz_S_bc;
  abcd[1911] = 4.0E0*I_ERI_D2y_Dxy_Dyz_S_bc;
  abcd[1912] = 4.0E0*I_ERI_Dyz_Dxy_Dyz_S_bc;
  abcd[1913] = 4.0E0*I_ERI_D2z_Dxy_Dyz_S_bc;
  abcd[1914] = 4.0E0*I_ERI_D2x_D2y_Dyz_S_bc-2.0E0*1*I_ERI_D2x_S_Dyz_S_c;
  abcd[1915] = 4.0E0*I_ERI_Dxy_D2y_Dyz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dyz_S_c;
  abcd[1916] = 4.0E0*I_ERI_Dxz_D2y_Dyz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dyz_S_c;
  abcd[1917] = 4.0E0*I_ERI_D2y_D2y_Dyz_S_bc-2.0E0*1*I_ERI_D2y_S_Dyz_S_c;
  abcd[1918] = 4.0E0*I_ERI_Dyz_D2y_Dyz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dyz_S_c;
  abcd[1919] = 4.0E0*I_ERI_D2z_D2y_Dyz_S_bc-2.0E0*1*I_ERI_D2z_S_Dyz_S_c;
  abcd[1920] = 4.0E0*I_ERI_D2x_Dyz_Dyz_S_bc;
  abcd[1921] = 4.0E0*I_ERI_Dxy_Dyz_Dyz_S_bc;
  abcd[1922] = 4.0E0*I_ERI_Dxz_Dyz_Dyz_S_bc;
  abcd[1923] = 4.0E0*I_ERI_D2y_Dyz_Dyz_S_bc;
  abcd[1924] = 4.0E0*I_ERI_Dyz_Dyz_Dyz_S_bc;
  abcd[1925] = 4.0E0*I_ERI_D2z_Dyz_Dyz_S_bc;
  abcd[1926] = 4.0E0*I_ERI_D2x_Dxy_D2z_S_bc-2.0E0*1*I_ERI_D2x_Dxy_S_S_b;
  abcd[1927] = 4.0E0*I_ERI_Dxy_Dxy_D2z_S_bc-2.0E0*1*I_ERI_Dxy_Dxy_S_S_b;
  abcd[1928] = 4.0E0*I_ERI_Dxz_Dxy_D2z_S_bc-2.0E0*1*I_ERI_Dxz_Dxy_S_S_b;
  abcd[1929] = 4.0E0*I_ERI_D2y_Dxy_D2z_S_bc-2.0E0*1*I_ERI_D2y_Dxy_S_S_b;
  abcd[1930] = 4.0E0*I_ERI_Dyz_Dxy_D2z_S_bc-2.0E0*1*I_ERI_Dyz_Dxy_S_S_b;
  abcd[1931] = 4.0E0*I_ERI_D2z_Dxy_D2z_S_bc-2.0E0*1*I_ERI_D2z_Dxy_S_S_b;
  abcd[1932] = 4.0E0*I_ERI_D2x_D2y_D2z_S_bc-2.0E0*1*I_ERI_D2x_D2y_S_S_b-2.0E0*1*I_ERI_D2x_S_D2z_S_c+1*I_ERI_D2x_S_S_S;
  abcd[1933] = 4.0E0*I_ERI_Dxy_D2y_D2z_S_bc-2.0E0*1*I_ERI_Dxy_D2y_S_S_b-2.0E0*1*I_ERI_Dxy_S_D2z_S_c+1*I_ERI_Dxy_S_S_S;
  abcd[1934] = 4.0E0*I_ERI_Dxz_D2y_D2z_S_bc-2.0E0*1*I_ERI_Dxz_D2y_S_S_b-2.0E0*1*I_ERI_Dxz_S_D2z_S_c+1*I_ERI_Dxz_S_S_S;
  abcd[1935] = 4.0E0*I_ERI_D2y_D2y_D2z_S_bc-2.0E0*1*I_ERI_D2y_D2y_S_S_b-2.0E0*1*I_ERI_D2y_S_D2z_S_c+1*I_ERI_D2y_S_S_S;
  abcd[1936] = 4.0E0*I_ERI_Dyz_D2y_D2z_S_bc-2.0E0*1*I_ERI_Dyz_D2y_S_S_b-2.0E0*1*I_ERI_Dyz_S_D2z_S_c+1*I_ERI_Dyz_S_S_S;
  abcd[1937] = 4.0E0*I_ERI_D2z_D2y_D2z_S_bc-2.0E0*1*I_ERI_D2z_D2y_S_S_b-2.0E0*1*I_ERI_D2z_S_D2z_S_c+1*I_ERI_D2z_S_S_S;
  abcd[1938] = 4.0E0*I_ERI_D2x_Dyz_D2z_S_bc-2.0E0*1*I_ERI_D2x_Dyz_S_S_b;
  abcd[1939] = 4.0E0*I_ERI_Dxy_Dyz_D2z_S_bc-2.0E0*1*I_ERI_Dxy_Dyz_S_S_b;
  abcd[1940] = 4.0E0*I_ERI_Dxz_Dyz_D2z_S_bc-2.0E0*1*I_ERI_Dxz_Dyz_S_S_b;
  abcd[1941] = 4.0E0*I_ERI_D2y_Dyz_D2z_S_bc-2.0E0*1*I_ERI_D2y_Dyz_S_S_b;
  abcd[1942] = 4.0E0*I_ERI_Dyz_Dyz_D2z_S_bc-2.0E0*1*I_ERI_Dyz_Dyz_S_S_b;
  abcd[1943] = 4.0E0*I_ERI_D2z_Dyz_D2z_S_bc-2.0E0*1*I_ERI_D2z_Dyz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[1944] = 4.0E0*I_ERI_D2x_Dxz_D2x_S_bc-2.0E0*1*I_ERI_D2x_Dxz_S_S_b;
  abcd[1945] = 4.0E0*I_ERI_Dxy_Dxz_D2x_S_bc-2.0E0*1*I_ERI_Dxy_Dxz_S_S_b;
  abcd[1946] = 4.0E0*I_ERI_Dxz_Dxz_D2x_S_bc-2.0E0*1*I_ERI_Dxz_Dxz_S_S_b;
  abcd[1947] = 4.0E0*I_ERI_D2y_Dxz_D2x_S_bc-2.0E0*1*I_ERI_D2y_Dxz_S_S_b;
  abcd[1948] = 4.0E0*I_ERI_Dyz_Dxz_D2x_S_bc-2.0E0*1*I_ERI_Dyz_Dxz_S_S_b;
  abcd[1949] = 4.0E0*I_ERI_D2z_Dxz_D2x_S_bc-2.0E0*1*I_ERI_D2z_Dxz_S_S_b;
  abcd[1950] = 4.0E0*I_ERI_D2x_Dyz_D2x_S_bc-2.0E0*1*I_ERI_D2x_Dyz_S_S_b;
  abcd[1951] = 4.0E0*I_ERI_Dxy_Dyz_D2x_S_bc-2.0E0*1*I_ERI_Dxy_Dyz_S_S_b;
  abcd[1952] = 4.0E0*I_ERI_Dxz_Dyz_D2x_S_bc-2.0E0*1*I_ERI_Dxz_Dyz_S_S_b;
  abcd[1953] = 4.0E0*I_ERI_D2y_Dyz_D2x_S_bc-2.0E0*1*I_ERI_D2y_Dyz_S_S_b;
  abcd[1954] = 4.0E0*I_ERI_Dyz_Dyz_D2x_S_bc-2.0E0*1*I_ERI_Dyz_Dyz_S_S_b;
  abcd[1955] = 4.0E0*I_ERI_D2z_Dyz_D2x_S_bc-2.0E0*1*I_ERI_D2z_Dyz_S_S_b;
  abcd[1956] = 4.0E0*I_ERI_D2x_D2z_D2x_S_bc-2.0E0*1*I_ERI_D2x_D2z_S_S_b-2.0E0*1*I_ERI_D2x_S_D2x_S_c+1*I_ERI_D2x_S_S_S;
  abcd[1957] = 4.0E0*I_ERI_Dxy_D2z_D2x_S_bc-2.0E0*1*I_ERI_Dxy_D2z_S_S_b-2.0E0*1*I_ERI_Dxy_S_D2x_S_c+1*I_ERI_Dxy_S_S_S;
  abcd[1958] = 4.0E0*I_ERI_Dxz_D2z_D2x_S_bc-2.0E0*1*I_ERI_Dxz_D2z_S_S_b-2.0E0*1*I_ERI_Dxz_S_D2x_S_c+1*I_ERI_Dxz_S_S_S;
  abcd[1959] = 4.0E0*I_ERI_D2y_D2z_D2x_S_bc-2.0E0*1*I_ERI_D2y_D2z_S_S_b-2.0E0*1*I_ERI_D2y_S_D2x_S_c+1*I_ERI_D2y_S_S_S;
  abcd[1960] = 4.0E0*I_ERI_Dyz_D2z_D2x_S_bc-2.0E0*1*I_ERI_Dyz_D2z_S_S_b-2.0E0*1*I_ERI_Dyz_S_D2x_S_c+1*I_ERI_Dyz_S_S_S;
  abcd[1961] = 4.0E0*I_ERI_D2z_D2z_D2x_S_bc-2.0E0*1*I_ERI_D2z_D2z_S_S_b-2.0E0*1*I_ERI_D2z_S_D2x_S_c+1*I_ERI_D2z_S_S_S;
  abcd[1962] = 4.0E0*I_ERI_D2x_Dxz_Dxy_S_bc;
  abcd[1963] = 4.0E0*I_ERI_Dxy_Dxz_Dxy_S_bc;
  abcd[1964] = 4.0E0*I_ERI_Dxz_Dxz_Dxy_S_bc;
  abcd[1965] = 4.0E0*I_ERI_D2y_Dxz_Dxy_S_bc;
  abcd[1966] = 4.0E0*I_ERI_Dyz_Dxz_Dxy_S_bc;
  abcd[1967] = 4.0E0*I_ERI_D2z_Dxz_Dxy_S_bc;
  abcd[1968] = 4.0E0*I_ERI_D2x_Dyz_Dxy_S_bc;
  abcd[1969] = 4.0E0*I_ERI_Dxy_Dyz_Dxy_S_bc;
  abcd[1970] = 4.0E0*I_ERI_Dxz_Dyz_Dxy_S_bc;
  abcd[1971] = 4.0E0*I_ERI_D2y_Dyz_Dxy_S_bc;
  abcd[1972] = 4.0E0*I_ERI_Dyz_Dyz_Dxy_S_bc;
  abcd[1973] = 4.0E0*I_ERI_D2z_Dyz_Dxy_S_bc;
  abcd[1974] = 4.0E0*I_ERI_D2x_D2z_Dxy_S_bc-2.0E0*1*I_ERI_D2x_S_Dxy_S_c;
  abcd[1975] = 4.0E0*I_ERI_Dxy_D2z_Dxy_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxy_S_c;
  abcd[1976] = 4.0E0*I_ERI_Dxz_D2z_Dxy_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxy_S_c;
  abcd[1977] = 4.0E0*I_ERI_D2y_D2z_Dxy_S_bc-2.0E0*1*I_ERI_D2y_S_Dxy_S_c;
  abcd[1978] = 4.0E0*I_ERI_Dyz_D2z_Dxy_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxy_S_c;
  abcd[1979] = 4.0E0*I_ERI_D2z_D2z_Dxy_S_bc-2.0E0*1*I_ERI_D2z_S_Dxy_S_c;
  abcd[1980] = 4.0E0*I_ERI_D2x_Dxz_Dxz_S_bc;
  abcd[1981] = 4.0E0*I_ERI_Dxy_Dxz_Dxz_S_bc;
  abcd[1982] = 4.0E0*I_ERI_Dxz_Dxz_Dxz_S_bc;
  abcd[1983] = 4.0E0*I_ERI_D2y_Dxz_Dxz_S_bc;
  abcd[1984] = 4.0E0*I_ERI_Dyz_Dxz_Dxz_S_bc;
  abcd[1985] = 4.0E0*I_ERI_D2z_Dxz_Dxz_S_bc;
  abcd[1986] = 4.0E0*I_ERI_D2x_Dyz_Dxz_S_bc;
  abcd[1987] = 4.0E0*I_ERI_Dxy_Dyz_Dxz_S_bc;
  abcd[1988] = 4.0E0*I_ERI_Dxz_Dyz_Dxz_S_bc;
  abcd[1989] = 4.0E0*I_ERI_D2y_Dyz_Dxz_S_bc;
  abcd[1990] = 4.0E0*I_ERI_Dyz_Dyz_Dxz_S_bc;
  abcd[1991] = 4.0E0*I_ERI_D2z_Dyz_Dxz_S_bc;
  abcd[1992] = 4.0E0*I_ERI_D2x_D2z_Dxz_S_bc-2.0E0*1*I_ERI_D2x_S_Dxz_S_c;
  abcd[1993] = 4.0E0*I_ERI_Dxy_D2z_Dxz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxz_S_c;
  abcd[1994] = 4.0E0*I_ERI_Dxz_D2z_Dxz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxz_S_c;
  abcd[1995] = 4.0E0*I_ERI_D2y_D2z_Dxz_S_bc-2.0E0*1*I_ERI_D2y_S_Dxz_S_c;
  abcd[1996] = 4.0E0*I_ERI_Dyz_D2z_Dxz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxz_S_c;
  abcd[1997] = 4.0E0*I_ERI_D2z_D2z_Dxz_S_bc-2.0E0*1*I_ERI_D2z_S_Dxz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[1998] = 4.0E0*I_ERI_D2x_Dxz_Dxy_S_bc;
  abcd[1999] = 4.0E0*I_ERI_Dxy_Dxz_Dxy_S_bc;
  abcd[2000] = 4.0E0*I_ERI_Dxz_Dxz_Dxy_S_bc;
  abcd[2001] = 4.0E0*I_ERI_D2y_Dxz_Dxy_S_bc;
  abcd[2002] = 4.0E0*I_ERI_Dyz_Dxz_Dxy_S_bc;
  abcd[2003] = 4.0E0*I_ERI_D2z_Dxz_Dxy_S_bc;
  abcd[2004] = 4.0E0*I_ERI_D2x_Dyz_Dxy_S_bc;
  abcd[2005] = 4.0E0*I_ERI_Dxy_Dyz_Dxy_S_bc;
  abcd[2006] = 4.0E0*I_ERI_Dxz_Dyz_Dxy_S_bc;
  abcd[2007] = 4.0E0*I_ERI_D2y_Dyz_Dxy_S_bc;
  abcd[2008] = 4.0E0*I_ERI_Dyz_Dyz_Dxy_S_bc;
  abcd[2009] = 4.0E0*I_ERI_D2z_Dyz_Dxy_S_bc;
  abcd[2010] = 4.0E0*I_ERI_D2x_D2z_Dxy_S_bc-2.0E0*1*I_ERI_D2x_S_Dxy_S_c;
  abcd[2011] = 4.0E0*I_ERI_Dxy_D2z_Dxy_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxy_S_c;
  abcd[2012] = 4.0E0*I_ERI_Dxz_D2z_Dxy_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxy_S_c;
  abcd[2013] = 4.0E0*I_ERI_D2y_D2z_Dxy_S_bc-2.0E0*1*I_ERI_D2y_S_Dxy_S_c;
  abcd[2014] = 4.0E0*I_ERI_Dyz_D2z_Dxy_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxy_S_c;
  abcd[2015] = 4.0E0*I_ERI_D2z_D2z_Dxy_S_bc-2.0E0*1*I_ERI_D2z_S_Dxy_S_c;
  abcd[2016] = 4.0E0*I_ERI_D2x_Dxz_D2y_S_bc-2.0E0*1*I_ERI_D2x_Dxz_S_S_b;
  abcd[2017] = 4.0E0*I_ERI_Dxy_Dxz_D2y_S_bc-2.0E0*1*I_ERI_Dxy_Dxz_S_S_b;
  abcd[2018] = 4.0E0*I_ERI_Dxz_Dxz_D2y_S_bc-2.0E0*1*I_ERI_Dxz_Dxz_S_S_b;
  abcd[2019] = 4.0E0*I_ERI_D2y_Dxz_D2y_S_bc-2.0E0*1*I_ERI_D2y_Dxz_S_S_b;
  abcd[2020] = 4.0E0*I_ERI_Dyz_Dxz_D2y_S_bc-2.0E0*1*I_ERI_Dyz_Dxz_S_S_b;
  abcd[2021] = 4.0E0*I_ERI_D2z_Dxz_D2y_S_bc-2.0E0*1*I_ERI_D2z_Dxz_S_S_b;
  abcd[2022] = 4.0E0*I_ERI_D2x_Dyz_D2y_S_bc-2.0E0*1*I_ERI_D2x_Dyz_S_S_b;
  abcd[2023] = 4.0E0*I_ERI_Dxy_Dyz_D2y_S_bc-2.0E0*1*I_ERI_Dxy_Dyz_S_S_b;
  abcd[2024] = 4.0E0*I_ERI_Dxz_Dyz_D2y_S_bc-2.0E0*1*I_ERI_Dxz_Dyz_S_S_b;
  abcd[2025] = 4.0E0*I_ERI_D2y_Dyz_D2y_S_bc-2.0E0*1*I_ERI_D2y_Dyz_S_S_b;
  abcd[2026] = 4.0E0*I_ERI_Dyz_Dyz_D2y_S_bc-2.0E0*1*I_ERI_Dyz_Dyz_S_S_b;
  abcd[2027] = 4.0E0*I_ERI_D2z_Dyz_D2y_S_bc-2.0E0*1*I_ERI_D2z_Dyz_S_S_b;
  abcd[2028] = 4.0E0*I_ERI_D2x_D2z_D2y_S_bc-2.0E0*1*I_ERI_D2x_D2z_S_S_b-2.0E0*1*I_ERI_D2x_S_D2y_S_c+1*I_ERI_D2x_S_S_S;
  abcd[2029] = 4.0E0*I_ERI_Dxy_D2z_D2y_S_bc-2.0E0*1*I_ERI_Dxy_D2z_S_S_b-2.0E0*1*I_ERI_Dxy_S_D2y_S_c+1*I_ERI_Dxy_S_S_S;
  abcd[2030] = 4.0E0*I_ERI_Dxz_D2z_D2y_S_bc-2.0E0*1*I_ERI_Dxz_D2z_S_S_b-2.0E0*1*I_ERI_Dxz_S_D2y_S_c+1*I_ERI_Dxz_S_S_S;
  abcd[2031] = 4.0E0*I_ERI_D2y_D2z_D2y_S_bc-2.0E0*1*I_ERI_D2y_D2z_S_S_b-2.0E0*1*I_ERI_D2y_S_D2y_S_c+1*I_ERI_D2y_S_S_S;
  abcd[2032] = 4.0E0*I_ERI_Dyz_D2z_D2y_S_bc-2.0E0*1*I_ERI_Dyz_D2z_S_S_b-2.0E0*1*I_ERI_Dyz_S_D2y_S_c+1*I_ERI_Dyz_S_S_S;
  abcd[2033] = 4.0E0*I_ERI_D2z_D2z_D2y_S_bc-2.0E0*1*I_ERI_D2z_D2z_S_S_b-2.0E0*1*I_ERI_D2z_S_D2y_S_c+1*I_ERI_D2z_S_S_S;
  abcd[2034] = 4.0E0*I_ERI_D2x_Dxz_Dyz_S_bc;
  abcd[2035] = 4.0E0*I_ERI_Dxy_Dxz_Dyz_S_bc;
  abcd[2036] = 4.0E0*I_ERI_Dxz_Dxz_Dyz_S_bc;
  abcd[2037] = 4.0E0*I_ERI_D2y_Dxz_Dyz_S_bc;
  abcd[2038] = 4.0E0*I_ERI_Dyz_Dxz_Dyz_S_bc;
  abcd[2039] = 4.0E0*I_ERI_D2z_Dxz_Dyz_S_bc;
  abcd[2040] = 4.0E0*I_ERI_D2x_Dyz_Dyz_S_bc;
  abcd[2041] = 4.0E0*I_ERI_Dxy_Dyz_Dyz_S_bc;
  abcd[2042] = 4.0E0*I_ERI_Dxz_Dyz_Dyz_S_bc;
  abcd[2043] = 4.0E0*I_ERI_D2y_Dyz_Dyz_S_bc;
  abcd[2044] = 4.0E0*I_ERI_Dyz_Dyz_Dyz_S_bc;
  abcd[2045] = 4.0E0*I_ERI_D2z_Dyz_Dyz_S_bc;
  abcd[2046] = 4.0E0*I_ERI_D2x_D2z_Dyz_S_bc-2.0E0*1*I_ERI_D2x_S_Dyz_S_c;
  abcd[2047] = 4.0E0*I_ERI_Dxy_D2z_Dyz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dyz_S_c;
  abcd[2048] = 4.0E0*I_ERI_Dxz_D2z_Dyz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dyz_S_c;
  abcd[2049] = 4.0E0*I_ERI_D2y_D2z_Dyz_S_bc-2.0E0*1*I_ERI_D2y_S_Dyz_S_c;
  abcd[2050] = 4.0E0*I_ERI_Dyz_D2z_Dyz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dyz_S_c;
  abcd[2051] = 4.0E0*I_ERI_D2z_D2z_Dyz_S_bc-2.0E0*1*I_ERI_D2z_S_Dyz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_D_S_bc
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_S_D_S_c
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  abcd[2052] = 4.0E0*I_ERI_D2x_Dxz_Dxz_S_bc;
  abcd[2053] = 4.0E0*I_ERI_Dxy_Dxz_Dxz_S_bc;
  abcd[2054] = 4.0E0*I_ERI_Dxz_Dxz_Dxz_S_bc;
  abcd[2055] = 4.0E0*I_ERI_D2y_Dxz_Dxz_S_bc;
  abcd[2056] = 4.0E0*I_ERI_Dyz_Dxz_Dxz_S_bc;
  abcd[2057] = 4.0E0*I_ERI_D2z_Dxz_Dxz_S_bc;
  abcd[2058] = 4.0E0*I_ERI_D2x_Dyz_Dxz_S_bc;
  abcd[2059] = 4.0E0*I_ERI_Dxy_Dyz_Dxz_S_bc;
  abcd[2060] = 4.0E0*I_ERI_Dxz_Dyz_Dxz_S_bc;
  abcd[2061] = 4.0E0*I_ERI_D2y_Dyz_Dxz_S_bc;
  abcd[2062] = 4.0E0*I_ERI_Dyz_Dyz_Dxz_S_bc;
  abcd[2063] = 4.0E0*I_ERI_D2z_Dyz_Dxz_S_bc;
  abcd[2064] = 4.0E0*I_ERI_D2x_D2z_Dxz_S_bc-2.0E0*1*I_ERI_D2x_S_Dxz_S_c;
  abcd[2065] = 4.0E0*I_ERI_Dxy_D2z_Dxz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dxz_S_c;
  abcd[2066] = 4.0E0*I_ERI_Dxz_D2z_Dxz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dxz_S_c;
  abcd[2067] = 4.0E0*I_ERI_D2y_D2z_Dxz_S_bc-2.0E0*1*I_ERI_D2y_S_Dxz_S_c;
  abcd[2068] = 4.0E0*I_ERI_Dyz_D2z_Dxz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dxz_S_c;
  abcd[2069] = 4.0E0*I_ERI_D2z_D2z_Dxz_S_bc-2.0E0*1*I_ERI_D2z_S_Dxz_S_c;
  abcd[2070] = 4.0E0*I_ERI_D2x_Dxz_Dyz_S_bc;
  abcd[2071] = 4.0E0*I_ERI_Dxy_Dxz_Dyz_S_bc;
  abcd[2072] = 4.0E0*I_ERI_Dxz_Dxz_Dyz_S_bc;
  abcd[2073] = 4.0E0*I_ERI_D2y_Dxz_Dyz_S_bc;
  abcd[2074] = 4.0E0*I_ERI_Dyz_Dxz_Dyz_S_bc;
  abcd[2075] = 4.0E0*I_ERI_D2z_Dxz_Dyz_S_bc;
  abcd[2076] = 4.0E0*I_ERI_D2x_Dyz_Dyz_S_bc;
  abcd[2077] = 4.0E0*I_ERI_Dxy_Dyz_Dyz_S_bc;
  abcd[2078] = 4.0E0*I_ERI_Dxz_Dyz_Dyz_S_bc;
  abcd[2079] = 4.0E0*I_ERI_D2y_Dyz_Dyz_S_bc;
  abcd[2080] = 4.0E0*I_ERI_Dyz_Dyz_Dyz_S_bc;
  abcd[2081] = 4.0E0*I_ERI_D2z_Dyz_Dyz_S_bc;
  abcd[2082] = 4.0E0*I_ERI_D2x_D2z_Dyz_S_bc-2.0E0*1*I_ERI_D2x_S_Dyz_S_c;
  abcd[2083] = 4.0E0*I_ERI_Dxy_D2z_Dyz_S_bc-2.0E0*1*I_ERI_Dxy_S_Dyz_S_c;
  abcd[2084] = 4.0E0*I_ERI_Dxz_D2z_Dyz_S_bc-2.0E0*1*I_ERI_Dxz_S_Dyz_S_c;
  abcd[2085] = 4.0E0*I_ERI_D2y_D2z_Dyz_S_bc-2.0E0*1*I_ERI_D2y_S_Dyz_S_c;
  abcd[2086] = 4.0E0*I_ERI_Dyz_D2z_Dyz_S_bc-2.0E0*1*I_ERI_Dyz_S_Dyz_S_c;
  abcd[2087] = 4.0E0*I_ERI_D2z_D2z_Dyz_S_bc-2.0E0*1*I_ERI_D2z_S_Dyz_S_c;
  abcd[2088] = 4.0E0*I_ERI_D2x_Dxz_D2z_S_bc-2.0E0*1*I_ERI_D2x_Dxz_S_S_b;
  abcd[2089] = 4.0E0*I_ERI_Dxy_Dxz_D2z_S_bc-2.0E0*1*I_ERI_Dxy_Dxz_S_S_b;
  abcd[2090] = 4.0E0*I_ERI_Dxz_Dxz_D2z_S_bc-2.0E0*1*I_ERI_Dxz_Dxz_S_S_b;
  abcd[2091] = 4.0E0*I_ERI_D2y_Dxz_D2z_S_bc-2.0E0*1*I_ERI_D2y_Dxz_S_S_b;
  abcd[2092] = 4.0E0*I_ERI_Dyz_Dxz_D2z_S_bc-2.0E0*1*I_ERI_Dyz_Dxz_S_S_b;
  abcd[2093] = 4.0E0*I_ERI_D2z_Dxz_D2z_S_bc-2.0E0*1*I_ERI_D2z_Dxz_S_S_b;
  abcd[2094] = 4.0E0*I_ERI_D2x_Dyz_D2z_S_bc-2.0E0*1*I_ERI_D2x_Dyz_S_S_b;
  abcd[2095] = 4.0E0*I_ERI_Dxy_Dyz_D2z_S_bc-2.0E0*1*I_ERI_Dxy_Dyz_S_S_b;
  abcd[2096] = 4.0E0*I_ERI_Dxz_Dyz_D2z_S_bc-2.0E0*1*I_ERI_Dxz_Dyz_S_S_b;
  abcd[2097] = 4.0E0*I_ERI_D2y_Dyz_D2z_S_bc-2.0E0*1*I_ERI_D2y_Dyz_S_S_b;
  abcd[2098] = 4.0E0*I_ERI_Dyz_Dyz_D2z_S_bc-2.0E0*1*I_ERI_Dyz_Dyz_S_S_b;
  abcd[2099] = 4.0E0*I_ERI_D2z_Dyz_D2z_S_bc-2.0E0*1*I_ERI_D2z_Dyz_S_S_b;
  abcd[2100] = 4.0E0*I_ERI_D2x_D2z_D2z_S_bc-2.0E0*1*I_ERI_D2x_D2z_S_S_b-2.0E0*1*I_ERI_D2x_S_D2z_S_c+1*I_ERI_D2x_S_S_S;
  abcd[2101] = 4.0E0*I_ERI_Dxy_D2z_D2z_S_bc-2.0E0*1*I_ERI_Dxy_D2z_S_S_b-2.0E0*1*I_ERI_Dxy_S_D2z_S_c+1*I_ERI_Dxy_S_S_S;
  abcd[2102] = 4.0E0*I_ERI_Dxz_D2z_D2z_S_bc-2.0E0*1*I_ERI_Dxz_D2z_S_S_b-2.0E0*1*I_ERI_Dxz_S_D2z_S_c+1*I_ERI_Dxz_S_S_S;
  abcd[2103] = 4.0E0*I_ERI_D2y_D2z_D2z_S_bc-2.0E0*1*I_ERI_D2y_D2z_S_S_b-2.0E0*1*I_ERI_D2y_S_D2z_S_c+1*I_ERI_D2y_S_S_S;
  abcd[2104] = 4.0E0*I_ERI_Dyz_D2z_D2z_S_bc-2.0E0*1*I_ERI_Dyz_D2z_S_S_b-2.0E0*1*I_ERI_Dyz_S_D2z_S_c+1*I_ERI_Dyz_S_S_S;
  abcd[2105] = 4.0E0*I_ERI_D2z_D2z_D2z_S_bc-2.0E0*1*I_ERI_D2z_D2z_S_S_b-2.0E0*1*I_ERI_D2z_S_D2z_S_c+1*I_ERI_D2z_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_F_S_cc
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   ************************************************************/
  abcd[2106] = 4.0E0*I_ERI_D2x_Px_F3x_S_cc-2.0E0*1*I_ERI_D2x_Px_Px_S_c-2.0E0*2*I_ERI_D2x_Px_Px_S_c;
  abcd[2107] = 4.0E0*I_ERI_Dxy_Px_F3x_S_cc-2.0E0*1*I_ERI_Dxy_Px_Px_S_c-2.0E0*2*I_ERI_Dxy_Px_Px_S_c;
  abcd[2108] = 4.0E0*I_ERI_Dxz_Px_F3x_S_cc-2.0E0*1*I_ERI_Dxz_Px_Px_S_c-2.0E0*2*I_ERI_Dxz_Px_Px_S_c;
  abcd[2109] = 4.0E0*I_ERI_D2y_Px_F3x_S_cc-2.0E0*1*I_ERI_D2y_Px_Px_S_c-2.0E0*2*I_ERI_D2y_Px_Px_S_c;
  abcd[2110] = 4.0E0*I_ERI_Dyz_Px_F3x_S_cc-2.0E0*1*I_ERI_Dyz_Px_Px_S_c-2.0E0*2*I_ERI_Dyz_Px_Px_S_c;
  abcd[2111] = 4.0E0*I_ERI_D2z_Px_F3x_S_cc-2.0E0*1*I_ERI_D2z_Px_Px_S_c-2.0E0*2*I_ERI_D2z_Px_Px_S_c;
  abcd[2112] = 4.0E0*I_ERI_D2x_Py_F3x_S_cc-2.0E0*1*I_ERI_D2x_Py_Px_S_c-2.0E0*2*I_ERI_D2x_Py_Px_S_c;
  abcd[2113] = 4.0E0*I_ERI_Dxy_Py_F3x_S_cc-2.0E0*1*I_ERI_Dxy_Py_Px_S_c-2.0E0*2*I_ERI_Dxy_Py_Px_S_c;
  abcd[2114] = 4.0E0*I_ERI_Dxz_Py_F3x_S_cc-2.0E0*1*I_ERI_Dxz_Py_Px_S_c-2.0E0*2*I_ERI_Dxz_Py_Px_S_c;
  abcd[2115] = 4.0E0*I_ERI_D2y_Py_F3x_S_cc-2.0E0*1*I_ERI_D2y_Py_Px_S_c-2.0E0*2*I_ERI_D2y_Py_Px_S_c;
  abcd[2116] = 4.0E0*I_ERI_Dyz_Py_F3x_S_cc-2.0E0*1*I_ERI_Dyz_Py_Px_S_c-2.0E0*2*I_ERI_Dyz_Py_Px_S_c;
  abcd[2117] = 4.0E0*I_ERI_D2z_Py_F3x_S_cc-2.0E0*1*I_ERI_D2z_Py_Px_S_c-2.0E0*2*I_ERI_D2z_Py_Px_S_c;
  abcd[2118] = 4.0E0*I_ERI_D2x_Pz_F3x_S_cc-2.0E0*1*I_ERI_D2x_Pz_Px_S_c-2.0E0*2*I_ERI_D2x_Pz_Px_S_c;
  abcd[2119] = 4.0E0*I_ERI_Dxy_Pz_F3x_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Px_S_c-2.0E0*2*I_ERI_Dxy_Pz_Px_S_c;
  abcd[2120] = 4.0E0*I_ERI_Dxz_Pz_F3x_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Px_S_c-2.0E0*2*I_ERI_Dxz_Pz_Px_S_c;
  abcd[2121] = 4.0E0*I_ERI_D2y_Pz_F3x_S_cc-2.0E0*1*I_ERI_D2y_Pz_Px_S_c-2.0E0*2*I_ERI_D2y_Pz_Px_S_c;
  abcd[2122] = 4.0E0*I_ERI_Dyz_Pz_F3x_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Px_S_c-2.0E0*2*I_ERI_Dyz_Pz_Px_S_c;
  abcd[2123] = 4.0E0*I_ERI_D2z_Pz_F3x_S_cc-2.0E0*1*I_ERI_D2z_Pz_Px_S_c-2.0E0*2*I_ERI_D2z_Pz_Px_S_c;
  abcd[2124] = 4.0E0*I_ERI_D2x_Px_F2xy_S_cc-2.0E0*1*I_ERI_D2x_Px_Py_S_c;
  abcd[2125] = 4.0E0*I_ERI_Dxy_Px_F2xy_S_cc-2.0E0*1*I_ERI_Dxy_Px_Py_S_c;
  abcd[2126] = 4.0E0*I_ERI_Dxz_Px_F2xy_S_cc-2.0E0*1*I_ERI_Dxz_Px_Py_S_c;
  abcd[2127] = 4.0E0*I_ERI_D2y_Px_F2xy_S_cc-2.0E0*1*I_ERI_D2y_Px_Py_S_c;
  abcd[2128] = 4.0E0*I_ERI_Dyz_Px_F2xy_S_cc-2.0E0*1*I_ERI_Dyz_Px_Py_S_c;
  abcd[2129] = 4.0E0*I_ERI_D2z_Px_F2xy_S_cc-2.0E0*1*I_ERI_D2z_Px_Py_S_c;
  abcd[2130] = 4.0E0*I_ERI_D2x_Py_F2xy_S_cc-2.0E0*1*I_ERI_D2x_Py_Py_S_c;
  abcd[2131] = 4.0E0*I_ERI_Dxy_Py_F2xy_S_cc-2.0E0*1*I_ERI_Dxy_Py_Py_S_c;
  abcd[2132] = 4.0E0*I_ERI_Dxz_Py_F2xy_S_cc-2.0E0*1*I_ERI_Dxz_Py_Py_S_c;
  abcd[2133] = 4.0E0*I_ERI_D2y_Py_F2xy_S_cc-2.0E0*1*I_ERI_D2y_Py_Py_S_c;
  abcd[2134] = 4.0E0*I_ERI_Dyz_Py_F2xy_S_cc-2.0E0*1*I_ERI_Dyz_Py_Py_S_c;
  abcd[2135] = 4.0E0*I_ERI_D2z_Py_F2xy_S_cc-2.0E0*1*I_ERI_D2z_Py_Py_S_c;
  abcd[2136] = 4.0E0*I_ERI_D2x_Pz_F2xy_S_cc-2.0E0*1*I_ERI_D2x_Pz_Py_S_c;
  abcd[2137] = 4.0E0*I_ERI_Dxy_Pz_F2xy_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Py_S_c;
  abcd[2138] = 4.0E0*I_ERI_Dxz_Pz_F2xy_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Py_S_c;
  abcd[2139] = 4.0E0*I_ERI_D2y_Pz_F2xy_S_cc-2.0E0*1*I_ERI_D2y_Pz_Py_S_c;
  abcd[2140] = 4.0E0*I_ERI_Dyz_Pz_F2xy_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Py_S_c;
  abcd[2141] = 4.0E0*I_ERI_D2z_Pz_F2xy_S_cc-2.0E0*1*I_ERI_D2z_Pz_Py_S_c;
  abcd[2142] = 4.0E0*I_ERI_D2x_Px_F2xz_S_cc-2.0E0*1*I_ERI_D2x_Px_Pz_S_c;
  abcd[2143] = 4.0E0*I_ERI_Dxy_Px_F2xz_S_cc-2.0E0*1*I_ERI_Dxy_Px_Pz_S_c;
  abcd[2144] = 4.0E0*I_ERI_Dxz_Px_F2xz_S_cc-2.0E0*1*I_ERI_Dxz_Px_Pz_S_c;
  abcd[2145] = 4.0E0*I_ERI_D2y_Px_F2xz_S_cc-2.0E0*1*I_ERI_D2y_Px_Pz_S_c;
  abcd[2146] = 4.0E0*I_ERI_Dyz_Px_F2xz_S_cc-2.0E0*1*I_ERI_Dyz_Px_Pz_S_c;
  abcd[2147] = 4.0E0*I_ERI_D2z_Px_F2xz_S_cc-2.0E0*1*I_ERI_D2z_Px_Pz_S_c;
  abcd[2148] = 4.0E0*I_ERI_D2x_Py_F2xz_S_cc-2.0E0*1*I_ERI_D2x_Py_Pz_S_c;
  abcd[2149] = 4.0E0*I_ERI_Dxy_Py_F2xz_S_cc-2.0E0*1*I_ERI_Dxy_Py_Pz_S_c;
  abcd[2150] = 4.0E0*I_ERI_Dxz_Py_F2xz_S_cc-2.0E0*1*I_ERI_Dxz_Py_Pz_S_c;
  abcd[2151] = 4.0E0*I_ERI_D2y_Py_F2xz_S_cc-2.0E0*1*I_ERI_D2y_Py_Pz_S_c;
  abcd[2152] = 4.0E0*I_ERI_Dyz_Py_F2xz_S_cc-2.0E0*1*I_ERI_Dyz_Py_Pz_S_c;
  abcd[2153] = 4.0E0*I_ERI_D2z_Py_F2xz_S_cc-2.0E0*1*I_ERI_D2z_Py_Pz_S_c;
  abcd[2154] = 4.0E0*I_ERI_D2x_Pz_F2xz_S_cc-2.0E0*1*I_ERI_D2x_Pz_Pz_S_c;
  abcd[2155] = 4.0E0*I_ERI_Dxy_Pz_F2xz_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_c;
  abcd[2156] = 4.0E0*I_ERI_Dxz_Pz_F2xz_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_c;
  abcd[2157] = 4.0E0*I_ERI_D2y_Pz_F2xz_S_cc-2.0E0*1*I_ERI_D2y_Pz_Pz_S_c;
  abcd[2158] = 4.0E0*I_ERI_Dyz_Pz_F2xz_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_c;
  abcd[2159] = 4.0E0*I_ERI_D2z_Pz_F2xz_S_cc-2.0E0*1*I_ERI_D2z_Pz_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_F_S_cc
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   ************************************************************/
  abcd[2160] = 4.0E0*I_ERI_D2x_Px_F2xy_S_cc-2.0E0*1*I_ERI_D2x_Px_Py_S_c;
  abcd[2161] = 4.0E0*I_ERI_Dxy_Px_F2xy_S_cc-2.0E0*1*I_ERI_Dxy_Px_Py_S_c;
  abcd[2162] = 4.0E0*I_ERI_Dxz_Px_F2xy_S_cc-2.0E0*1*I_ERI_Dxz_Px_Py_S_c;
  abcd[2163] = 4.0E0*I_ERI_D2y_Px_F2xy_S_cc-2.0E0*1*I_ERI_D2y_Px_Py_S_c;
  abcd[2164] = 4.0E0*I_ERI_Dyz_Px_F2xy_S_cc-2.0E0*1*I_ERI_Dyz_Px_Py_S_c;
  abcd[2165] = 4.0E0*I_ERI_D2z_Px_F2xy_S_cc-2.0E0*1*I_ERI_D2z_Px_Py_S_c;
  abcd[2166] = 4.0E0*I_ERI_D2x_Py_F2xy_S_cc-2.0E0*1*I_ERI_D2x_Py_Py_S_c;
  abcd[2167] = 4.0E0*I_ERI_Dxy_Py_F2xy_S_cc-2.0E0*1*I_ERI_Dxy_Py_Py_S_c;
  abcd[2168] = 4.0E0*I_ERI_Dxz_Py_F2xy_S_cc-2.0E0*1*I_ERI_Dxz_Py_Py_S_c;
  abcd[2169] = 4.0E0*I_ERI_D2y_Py_F2xy_S_cc-2.0E0*1*I_ERI_D2y_Py_Py_S_c;
  abcd[2170] = 4.0E0*I_ERI_Dyz_Py_F2xy_S_cc-2.0E0*1*I_ERI_Dyz_Py_Py_S_c;
  abcd[2171] = 4.0E0*I_ERI_D2z_Py_F2xy_S_cc-2.0E0*1*I_ERI_D2z_Py_Py_S_c;
  abcd[2172] = 4.0E0*I_ERI_D2x_Pz_F2xy_S_cc-2.0E0*1*I_ERI_D2x_Pz_Py_S_c;
  abcd[2173] = 4.0E0*I_ERI_Dxy_Pz_F2xy_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Py_S_c;
  abcd[2174] = 4.0E0*I_ERI_Dxz_Pz_F2xy_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Py_S_c;
  abcd[2175] = 4.0E0*I_ERI_D2y_Pz_F2xy_S_cc-2.0E0*1*I_ERI_D2y_Pz_Py_S_c;
  abcd[2176] = 4.0E0*I_ERI_Dyz_Pz_F2xy_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Py_S_c;
  abcd[2177] = 4.0E0*I_ERI_D2z_Pz_F2xy_S_cc-2.0E0*1*I_ERI_D2z_Pz_Py_S_c;
  abcd[2178] = 4.0E0*I_ERI_D2x_Px_Fx2y_S_cc-2.0E0*1*I_ERI_D2x_Px_Px_S_c;
  abcd[2179] = 4.0E0*I_ERI_Dxy_Px_Fx2y_S_cc-2.0E0*1*I_ERI_Dxy_Px_Px_S_c;
  abcd[2180] = 4.0E0*I_ERI_Dxz_Px_Fx2y_S_cc-2.0E0*1*I_ERI_Dxz_Px_Px_S_c;
  abcd[2181] = 4.0E0*I_ERI_D2y_Px_Fx2y_S_cc-2.0E0*1*I_ERI_D2y_Px_Px_S_c;
  abcd[2182] = 4.0E0*I_ERI_Dyz_Px_Fx2y_S_cc-2.0E0*1*I_ERI_Dyz_Px_Px_S_c;
  abcd[2183] = 4.0E0*I_ERI_D2z_Px_Fx2y_S_cc-2.0E0*1*I_ERI_D2z_Px_Px_S_c;
  abcd[2184] = 4.0E0*I_ERI_D2x_Py_Fx2y_S_cc-2.0E0*1*I_ERI_D2x_Py_Px_S_c;
  abcd[2185] = 4.0E0*I_ERI_Dxy_Py_Fx2y_S_cc-2.0E0*1*I_ERI_Dxy_Py_Px_S_c;
  abcd[2186] = 4.0E0*I_ERI_Dxz_Py_Fx2y_S_cc-2.0E0*1*I_ERI_Dxz_Py_Px_S_c;
  abcd[2187] = 4.0E0*I_ERI_D2y_Py_Fx2y_S_cc-2.0E0*1*I_ERI_D2y_Py_Px_S_c;
  abcd[2188] = 4.0E0*I_ERI_Dyz_Py_Fx2y_S_cc-2.0E0*1*I_ERI_Dyz_Py_Px_S_c;
  abcd[2189] = 4.0E0*I_ERI_D2z_Py_Fx2y_S_cc-2.0E0*1*I_ERI_D2z_Py_Px_S_c;
  abcd[2190] = 4.0E0*I_ERI_D2x_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_D2x_Pz_Px_S_c;
  abcd[2191] = 4.0E0*I_ERI_Dxy_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Px_S_c;
  abcd[2192] = 4.0E0*I_ERI_Dxz_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Px_S_c;
  abcd[2193] = 4.0E0*I_ERI_D2y_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_D2y_Pz_Px_S_c;
  abcd[2194] = 4.0E0*I_ERI_Dyz_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Px_S_c;
  abcd[2195] = 4.0E0*I_ERI_D2z_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_D2z_Pz_Px_S_c;
  abcd[2196] = 4.0E0*I_ERI_D2x_Px_Fxyz_S_cc;
  abcd[2197] = 4.0E0*I_ERI_Dxy_Px_Fxyz_S_cc;
  abcd[2198] = 4.0E0*I_ERI_Dxz_Px_Fxyz_S_cc;
  abcd[2199] = 4.0E0*I_ERI_D2y_Px_Fxyz_S_cc;
  abcd[2200] = 4.0E0*I_ERI_Dyz_Px_Fxyz_S_cc;
  abcd[2201] = 4.0E0*I_ERI_D2z_Px_Fxyz_S_cc;
  abcd[2202] = 4.0E0*I_ERI_D2x_Py_Fxyz_S_cc;
  abcd[2203] = 4.0E0*I_ERI_Dxy_Py_Fxyz_S_cc;
  abcd[2204] = 4.0E0*I_ERI_Dxz_Py_Fxyz_S_cc;
  abcd[2205] = 4.0E0*I_ERI_D2y_Py_Fxyz_S_cc;
  abcd[2206] = 4.0E0*I_ERI_Dyz_Py_Fxyz_S_cc;
  abcd[2207] = 4.0E0*I_ERI_D2z_Py_Fxyz_S_cc;
  abcd[2208] = 4.0E0*I_ERI_D2x_Pz_Fxyz_S_cc;
  abcd[2209] = 4.0E0*I_ERI_Dxy_Pz_Fxyz_S_cc;
  abcd[2210] = 4.0E0*I_ERI_Dxz_Pz_Fxyz_S_cc;
  abcd[2211] = 4.0E0*I_ERI_D2y_Pz_Fxyz_S_cc;
  abcd[2212] = 4.0E0*I_ERI_Dyz_Pz_Fxyz_S_cc;
  abcd[2213] = 4.0E0*I_ERI_D2z_Pz_Fxyz_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_F_S_cc
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   ************************************************************/
  abcd[2214] = 4.0E0*I_ERI_D2x_Px_F2xz_S_cc-2.0E0*1*I_ERI_D2x_Px_Pz_S_c;
  abcd[2215] = 4.0E0*I_ERI_Dxy_Px_F2xz_S_cc-2.0E0*1*I_ERI_Dxy_Px_Pz_S_c;
  abcd[2216] = 4.0E0*I_ERI_Dxz_Px_F2xz_S_cc-2.0E0*1*I_ERI_Dxz_Px_Pz_S_c;
  abcd[2217] = 4.0E0*I_ERI_D2y_Px_F2xz_S_cc-2.0E0*1*I_ERI_D2y_Px_Pz_S_c;
  abcd[2218] = 4.0E0*I_ERI_Dyz_Px_F2xz_S_cc-2.0E0*1*I_ERI_Dyz_Px_Pz_S_c;
  abcd[2219] = 4.0E0*I_ERI_D2z_Px_F2xz_S_cc-2.0E0*1*I_ERI_D2z_Px_Pz_S_c;
  abcd[2220] = 4.0E0*I_ERI_D2x_Py_F2xz_S_cc-2.0E0*1*I_ERI_D2x_Py_Pz_S_c;
  abcd[2221] = 4.0E0*I_ERI_Dxy_Py_F2xz_S_cc-2.0E0*1*I_ERI_Dxy_Py_Pz_S_c;
  abcd[2222] = 4.0E0*I_ERI_Dxz_Py_F2xz_S_cc-2.0E0*1*I_ERI_Dxz_Py_Pz_S_c;
  abcd[2223] = 4.0E0*I_ERI_D2y_Py_F2xz_S_cc-2.0E0*1*I_ERI_D2y_Py_Pz_S_c;
  abcd[2224] = 4.0E0*I_ERI_Dyz_Py_F2xz_S_cc-2.0E0*1*I_ERI_Dyz_Py_Pz_S_c;
  abcd[2225] = 4.0E0*I_ERI_D2z_Py_F2xz_S_cc-2.0E0*1*I_ERI_D2z_Py_Pz_S_c;
  abcd[2226] = 4.0E0*I_ERI_D2x_Pz_F2xz_S_cc-2.0E0*1*I_ERI_D2x_Pz_Pz_S_c;
  abcd[2227] = 4.0E0*I_ERI_Dxy_Pz_F2xz_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_c;
  abcd[2228] = 4.0E0*I_ERI_Dxz_Pz_F2xz_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_c;
  abcd[2229] = 4.0E0*I_ERI_D2y_Pz_F2xz_S_cc-2.0E0*1*I_ERI_D2y_Pz_Pz_S_c;
  abcd[2230] = 4.0E0*I_ERI_Dyz_Pz_F2xz_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_c;
  abcd[2231] = 4.0E0*I_ERI_D2z_Pz_F2xz_S_cc-2.0E0*1*I_ERI_D2z_Pz_Pz_S_c;
  abcd[2232] = 4.0E0*I_ERI_D2x_Px_Fxyz_S_cc;
  abcd[2233] = 4.0E0*I_ERI_Dxy_Px_Fxyz_S_cc;
  abcd[2234] = 4.0E0*I_ERI_Dxz_Px_Fxyz_S_cc;
  abcd[2235] = 4.0E0*I_ERI_D2y_Px_Fxyz_S_cc;
  abcd[2236] = 4.0E0*I_ERI_Dyz_Px_Fxyz_S_cc;
  abcd[2237] = 4.0E0*I_ERI_D2z_Px_Fxyz_S_cc;
  abcd[2238] = 4.0E0*I_ERI_D2x_Py_Fxyz_S_cc;
  abcd[2239] = 4.0E0*I_ERI_Dxy_Py_Fxyz_S_cc;
  abcd[2240] = 4.0E0*I_ERI_Dxz_Py_Fxyz_S_cc;
  abcd[2241] = 4.0E0*I_ERI_D2y_Py_Fxyz_S_cc;
  abcd[2242] = 4.0E0*I_ERI_Dyz_Py_Fxyz_S_cc;
  abcd[2243] = 4.0E0*I_ERI_D2z_Py_Fxyz_S_cc;
  abcd[2244] = 4.0E0*I_ERI_D2x_Pz_Fxyz_S_cc;
  abcd[2245] = 4.0E0*I_ERI_Dxy_Pz_Fxyz_S_cc;
  abcd[2246] = 4.0E0*I_ERI_Dxz_Pz_Fxyz_S_cc;
  abcd[2247] = 4.0E0*I_ERI_D2y_Pz_Fxyz_S_cc;
  abcd[2248] = 4.0E0*I_ERI_Dyz_Pz_Fxyz_S_cc;
  abcd[2249] = 4.0E0*I_ERI_D2z_Pz_Fxyz_S_cc;
  abcd[2250] = 4.0E0*I_ERI_D2x_Px_Fx2z_S_cc-2.0E0*1*I_ERI_D2x_Px_Px_S_c;
  abcd[2251] = 4.0E0*I_ERI_Dxy_Px_Fx2z_S_cc-2.0E0*1*I_ERI_Dxy_Px_Px_S_c;
  abcd[2252] = 4.0E0*I_ERI_Dxz_Px_Fx2z_S_cc-2.0E0*1*I_ERI_Dxz_Px_Px_S_c;
  abcd[2253] = 4.0E0*I_ERI_D2y_Px_Fx2z_S_cc-2.0E0*1*I_ERI_D2y_Px_Px_S_c;
  abcd[2254] = 4.0E0*I_ERI_Dyz_Px_Fx2z_S_cc-2.0E0*1*I_ERI_Dyz_Px_Px_S_c;
  abcd[2255] = 4.0E0*I_ERI_D2z_Px_Fx2z_S_cc-2.0E0*1*I_ERI_D2z_Px_Px_S_c;
  abcd[2256] = 4.0E0*I_ERI_D2x_Py_Fx2z_S_cc-2.0E0*1*I_ERI_D2x_Py_Px_S_c;
  abcd[2257] = 4.0E0*I_ERI_Dxy_Py_Fx2z_S_cc-2.0E0*1*I_ERI_Dxy_Py_Px_S_c;
  abcd[2258] = 4.0E0*I_ERI_Dxz_Py_Fx2z_S_cc-2.0E0*1*I_ERI_Dxz_Py_Px_S_c;
  abcd[2259] = 4.0E0*I_ERI_D2y_Py_Fx2z_S_cc-2.0E0*1*I_ERI_D2y_Py_Px_S_c;
  abcd[2260] = 4.0E0*I_ERI_Dyz_Py_Fx2z_S_cc-2.0E0*1*I_ERI_Dyz_Py_Px_S_c;
  abcd[2261] = 4.0E0*I_ERI_D2z_Py_Fx2z_S_cc-2.0E0*1*I_ERI_D2z_Py_Px_S_c;
  abcd[2262] = 4.0E0*I_ERI_D2x_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_D2x_Pz_Px_S_c;
  abcd[2263] = 4.0E0*I_ERI_Dxy_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Px_S_c;
  abcd[2264] = 4.0E0*I_ERI_Dxz_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Px_S_c;
  abcd[2265] = 4.0E0*I_ERI_D2y_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_D2y_Pz_Px_S_c;
  abcd[2266] = 4.0E0*I_ERI_Dyz_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Px_S_c;
  abcd[2267] = 4.0E0*I_ERI_D2z_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_D2z_Pz_Px_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_F_S_cc
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   ************************************************************/
  abcd[2268] = 4.0E0*I_ERI_D2x_Px_Fx2y_S_cc-2.0E0*1*I_ERI_D2x_Px_Px_S_c;
  abcd[2269] = 4.0E0*I_ERI_Dxy_Px_Fx2y_S_cc-2.0E0*1*I_ERI_Dxy_Px_Px_S_c;
  abcd[2270] = 4.0E0*I_ERI_Dxz_Px_Fx2y_S_cc-2.0E0*1*I_ERI_Dxz_Px_Px_S_c;
  abcd[2271] = 4.0E0*I_ERI_D2y_Px_Fx2y_S_cc-2.0E0*1*I_ERI_D2y_Px_Px_S_c;
  abcd[2272] = 4.0E0*I_ERI_Dyz_Px_Fx2y_S_cc-2.0E0*1*I_ERI_Dyz_Px_Px_S_c;
  abcd[2273] = 4.0E0*I_ERI_D2z_Px_Fx2y_S_cc-2.0E0*1*I_ERI_D2z_Px_Px_S_c;
  abcd[2274] = 4.0E0*I_ERI_D2x_Py_Fx2y_S_cc-2.0E0*1*I_ERI_D2x_Py_Px_S_c;
  abcd[2275] = 4.0E0*I_ERI_Dxy_Py_Fx2y_S_cc-2.0E0*1*I_ERI_Dxy_Py_Px_S_c;
  abcd[2276] = 4.0E0*I_ERI_Dxz_Py_Fx2y_S_cc-2.0E0*1*I_ERI_Dxz_Py_Px_S_c;
  abcd[2277] = 4.0E0*I_ERI_D2y_Py_Fx2y_S_cc-2.0E0*1*I_ERI_D2y_Py_Px_S_c;
  abcd[2278] = 4.0E0*I_ERI_Dyz_Py_Fx2y_S_cc-2.0E0*1*I_ERI_Dyz_Py_Px_S_c;
  abcd[2279] = 4.0E0*I_ERI_D2z_Py_Fx2y_S_cc-2.0E0*1*I_ERI_D2z_Py_Px_S_c;
  abcd[2280] = 4.0E0*I_ERI_D2x_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_D2x_Pz_Px_S_c;
  abcd[2281] = 4.0E0*I_ERI_Dxy_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Px_S_c;
  abcd[2282] = 4.0E0*I_ERI_Dxz_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Px_S_c;
  abcd[2283] = 4.0E0*I_ERI_D2y_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_D2y_Pz_Px_S_c;
  abcd[2284] = 4.0E0*I_ERI_Dyz_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Px_S_c;
  abcd[2285] = 4.0E0*I_ERI_D2z_Pz_Fx2y_S_cc-2.0E0*1*I_ERI_D2z_Pz_Px_S_c;
  abcd[2286] = 4.0E0*I_ERI_D2x_Px_F3y_S_cc-2.0E0*1*I_ERI_D2x_Px_Py_S_c-2.0E0*2*I_ERI_D2x_Px_Py_S_c;
  abcd[2287] = 4.0E0*I_ERI_Dxy_Px_F3y_S_cc-2.0E0*1*I_ERI_Dxy_Px_Py_S_c-2.0E0*2*I_ERI_Dxy_Px_Py_S_c;
  abcd[2288] = 4.0E0*I_ERI_Dxz_Px_F3y_S_cc-2.0E0*1*I_ERI_Dxz_Px_Py_S_c-2.0E0*2*I_ERI_Dxz_Px_Py_S_c;
  abcd[2289] = 4.0E0*I_ERI_D2y_Px_F3y_S_cc-2.0E0*1*I_ERI_D2y_Px_Py_S_c-2.0E0*2*I_ERI_D2y_Px_Py_S_c;
  abcd[2290] = 4.0E0*I_ERI_Dyz_Px_F3y_S_cc-2.0E0*1*I_ERI_Dyz_Px_Py_S_c-2.0E0*2*I_ERI_Dyz_Px_Py_S_c;
  abcd[2291] = 4.0E0*I_ERI_D2z_Px_F3y_S_cc-2.0E0*1*I_ERI_D2z_Px_Py_S_c-2.0E0*2*I_ERI_D2z_Px_Py_S_c;
  abcd[2292] = 4.0E0*I_ERI_D2x_Py_F3y_S_cc-2.0E0*1*I_ERI_D2x_Py_Py_S_c-2.0E0*2*I_ERI_D2x_Py_Py_S_c;
  abcd[2293] = 4.0E0*I_ERI_Dxy_Py_F3y_S_cc-2.0E0*1*I_ERI_Dxy_Py_Py_S_c-2.0E0*2*I_ERI_Dxy_Py_Py_S_c;
  abcd[2294] = 4.0E0*I_ERI_Dxz_Py_F3y_S_cc-2.0E0*1*I_ERI_Dxz_Py_Py_S_c-2.0E0*2*I_ERI_Dxz_Py_Py_S_c;
  abcd[2295] = 4.0E0*I_ERI_D2y_Py_F3y_S_cc-2.0E0*1*I_ERI_D2y_Py_Py_S_c-2.0E0*2*I_ERI_D2y_Py_Py_S_c;
  abcd[2296] = 4.0E0*I_ERI_Dyz_Py_F3y_S_cc-2.0E0*1*I_ERI_Dyz_Py_Py_S_c-2.0E0*2*I_ERI_Dyz_Py_Py_S_c;
  abcd[2297] = 4.0E0*I_ERI_D2z_Py_F3y_S_cc-2.0E0*1*I_ERI_D2z_Py_Py_S_c-2.0E0*2*I_ERI_D2z_Py_Py_S_c;
  abcd[2298] = 4.0E0*I_ERI_D2x_Pz_F3y_S_cc-2.0E0*1*I_ERI_D2x_Pz_Py_S_c-2.0E0*2*I_ERI_D2x_Pz_Py_S_c;
  abcd[2299] = 4.0E0*I_ERI_Dxy_Pz_F3y_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Py_S_c-2.0E0*2*I_ERI_Dxy_Pz_Py_S_c;
  abcd[2300] = 4.0E0*I_ERI_Dxz_Pz_F3y_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Py_S_c-2.0E0*2*I_ERI_Dxz_Pz_Py_S_c;
  abcd[2301] = 4.0E0*I_ERI_D2y_Pz_F3y_S_cc-2.0E0*1*I_ERI_D2y_Pz_Py_S_c-2.0E0*2*I_ERI_D2y_Pz_Py_S_c;
  abcd[2302] = 4.0E0*I_ERI_Dyz_Pz_F3y_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Py_S_c-2.0E0*2*I_ERI_Dyz_Pz_Py_S_c;
  abcd[2303] = 4.0E0*I_ERI_D2z_Pz_F3y_S_cc-2.0E0*1*I_ERI_D2z_Pz_Py_S_c-2.0E0*2*I_ERI_D2z_Pz_Py_S_c;
  abcd[2304] = 4.0E0*I_ERI_D2x_Px_F2yz_S_cc-2.0E0*1*I_ERI_D2x_Px_Pz_S_c;
  abcd[2305] = 4.0E0*I_ERI_Dxy_Px_F2yz_S_cc-2.0E0*1*I_ERI_Dxy_Px_Pz_S_c;
  abcd[2306] = 4.0E0*I_ERI_Dxz_Px_F2yz_S_cc-2.0E0*1*I_ERI_Dxz_Px_Pz_S_c;
  abcd[2307] = 4.0E0*I_ERI_D2y_Px_F2yz_S_cc-2.0E0*1*I_ERI_D2y_Px_Pz_S_c;
  abcd[2308] = 4.0E0*I_ERI_Dyz_Px_F2yz_S_cc-2.0E0*1*I_ERI_Dyz_Px_Pz_S_c;
  abcd[2309] = 4.0E0*I_ERI_D2z_Px_F2yz_S_cc-2.0E0*1*I_ERI_D2z_Px_Pz_S_c;
  abcd[2310] = 4.0E0*I_ERI_D2x_Py_F2yz_S_cc-2.0E0*1*I_ERI_D2x_Py_Pz_S_c;
  abcd[2311] = 4.0E0*I_ERI_Dxy_Py_F2yz_S_cc-2.0E0*1*I_ERI_Dxy_Py_Pz_S_c;
  abcd[2312] = 4.0E0*I_ERI_Dxz_Py_F2yz_S_cc-2.0E0*1*I_ERI_Dxz_Py_Pz_S_c;
  abcd[2313] = 4.0E0*I_ERI_D2y_Py_F2yz_S_cc-2.0E0*1*I_ERI_D2y_Py_Pz_S_c;
  abcd[2314] = 4.0E0*I_ERI_Dyz_Py_F2yz_S_cc-2.0E0*1*I_ERI_Dyz_Py_Pz_S_c;
  abcd[2315] = 4.0E0*I_ERI_D2z_Py_F2yz_S_cc-2.0E0*1*I_ERI_D2z_Py_Pz_S_c;
  abcd[2316] = 4.0E0*I_ERI_D2x_Pz_F2yz_S_cc-2.0E0*1*I_ERI_D2x_Pz_Pz_S_c;
  abcd[2317] = 4.0E0*I_ERI_Dxy_Pz_F2yz_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_c;
  abcd[2318] = 4.0E0*I_ERI_Dxz_Pz_F2yz_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_c;
  abcd[2319] = 4.0E0*I_ERI_D2y_Pz_F2yz_S_cc-2.0E0*1*I_ERI_D2y_Pz_Pz_S_c;
  abcd[2320] = 4.0E0*I_ERI_Dyz_Pz_F2yz_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_c;
  abcd[2321] = 4.0E0*I_ERI_D2z_Pz_F2yz_S_cc-2.0E0*1*I_ERI_D2z_Pz_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_F_S_cc
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   ************************************************************/
  abcd[2322] = 4.0E0*I_ERI_D2x_Px_Fxyz_S_cc;
  abcd[2323] = 4.0E0*I_ERI_Dxy_Px_Fxyz_S_cc;
  abcd[2324] = 4.0E0*I_ERI_Dxz_Px_Fxyz_S_cc;
  abcd[2325] = 4.0E0*I_ERI_D2y_Px_Fxyz_S_cc;
  abcd[2326] = 4.0E0*I_ERI_Dyz_Px_Fxyz_S_cc;
  abcd[2327] = 4.0E0*I_ERI_D2z_Px_Fxyz_S_cc;
  abcd[2328] = 4.0E0*I_ERI_D2x_Py_Fxyz_S_cc;
  abcd[2329] = 4.0E0*I_ERI_Dxy_Py_Fxyz_S_cc;
  abcd[2330] = 4.0E0*I_ERI_Dxz_Py_Fxyz_S_cc;
  abcd[2331] = 4.0E0*I_ERI_D2y_Py_Fxyz_S_cc;
  abcd[2332] = 4.0E0*I_ERI_Dyz_Py_Fxyz_S_cc;
  abcd[2333] = 4.0E0*I_ERI_D2z_Py_Fxyz_S_cc;
  abcd[2334] = 4.0E0*I_ERI_D2x_Pz_Fxyz_S_cc;
  abcd[2335] = 4.0E0*I_ERI_Dxy_Pz_Fxyz_S_cc;
  abcd[2336] = 4.0E0*I_ERI_Dxz_Pz_Fxyz_S_cc;
  abcd[2337] = 4.0E0*I_ERI_D2y_Pz_Fxyz_S_cc;
  abcd[2338] = 4.0E0*I_ERI_Dyz_Pz_Fxyz_S_cc;
  abcd[2339] = 4.0E0*I_ERI_D2z_Pz_Fxyz_S_cc;
  abcd[2340] = 4.0E0*I_ERI_D2x_Px_F2yz_S_cc-2.0E0*1*I_ERI_D2x_Px_Pz_S_c;
  abcd[2341] = 4.0E0*I_ERI_Dxy_Px_F2yz_S_cc-2.0E0*1*I_ERI_Dxy_Px_Pz_S_c;
  abcd[2342] = 4.0E0*I_ERI_Dxz_Px_F2yz_S_cc-2.0E0*1*I_ERI_Dxz_Px_Pz_S_c;
  abcd[2343] = 4.0E0*I_ERI_D2y_Px_F2yz_S_cc-2.0E0*1*I_ERI_D2y_Px_Pz_S_c;
  abcd[2344] = 4.0E0*I_ERI_Dyz_Px_F2yz_S_cc-2.0E0*1*I_ERI_Dyz_Px_Pz_S_c;
  abcd[2345] = 4.0E0*I_ERI_D2z_Px_F2yz_S_cc-2.0E0*1*I_ERI_D2z_Px_Pz_S_c;
  abcd[2346] = 4.0E0*I_ERI_D2x_Py_F2yz_S_cc-2.0E0*1*I_ERI_D2x_Py_Pz_S_c;
  abcd[2347] = 4.0E0*I_ERI_Dxy_Py_F2yz_S_cc-2.0E0*1*I_ERI_Dxy_Py_Pz_S_c;
  abcd[2348] = 4.0E0*I_ERI_Dxz_Py_F2yz_S_cc-2.0E0*1*I_ERI_Dxz_Py_Pz_S_c;
  abcd[2349] = 4.0E0*I_ERI_D2y_Py_F2yz_S_cc-2.0E0*1*I_ERI_D2y_Py_Pz_S_c;
  abcd[2350] = 4.0E0*I_ERI_Dyz_Py_F2yz_S_cc-2.0E0*1*I_ERI_Dyz_Py_Pz_S_c;
  abcd[2351] = 4.0E0*I_ERI_D2z_Py_F2yz_S_cc-2.0E0*1*I_ERI_D2z_Py_Pz_S_c;
  abcd[2352] = 4.0E0*I_ERI_D2x_Pz_F2yz_S_cc-2.0E0*1*I_ERI_D2x_Pz_Pz_S_c;
  abcd[2353] = 4.0E0*I_ERI_Dxy_Pz_F2yz_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_c;
  abcd[2354] = 4.0E0*I_ERI_Dxz_Pz_F2yz_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_c;
  abcd[2355] = 4.0E0*I_ERI_D2y_Pz_F2yz_S_cc-2.0E0*1*I_ERI_D2y_Pz_Pz_S_c;
  abcd[2356] = 4.0E0*I_ERI_Dyz_Pz_F2yz_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_c;
  abcd[2357] = 4.0E0*I_ERI_D2z_Pz_F2yz_S_cc-2.0E0*1*I_ERI_D2z_Pz_Pz_S_c;
  abcd[2358] = 4.0E0*I_ERI_D2x_Px_Fy2z_S_cc-2.0E0*1*I_ERI_D2x_Px_Py_S_c;
  abcd[2359] = 4.0E0*I_ERI_Dxy_Px_Fy2z_S_cc-2.0E0*1*I_ERI_Dxy_Px_Py_S_c;
  abcd[2360] = 4.0E0*I_ERI_Dxz_Px_Fy2z_S_cc-2.0E0*1*I_ERI_Dxz_Px_Py_S_c;
  abcd[2361] = 4.0E0*I_ERI_D2y_Px_Fy2z_S_cc-2.0E0*1*I_ERI_D2y_Px_Py_S_c;
  abcd[2362] = 4.0E0*I_ERI_Dyz_Px_Fy2z_S_cc-2.0E0*1*I_ERI_Dyz_Px_Py_S_c;
  abcd[2363] = 4.0E0*I_ERI_D2z_Px_Fy2z_S_cc-2.0E0*1*I_ERI_D2z_Px_Py_S_c;
  abcd[2364] = 4.0E0*I_ERI_D2x_Py_Fy2z_S_cc-2.0E0*1*I_ERI_D2x_Py_Py_S_c;
  abcd[2365] = 4.0E0*I_ERI_Dxy_Py_Fy2z_S_cc-2.0E0*1*I_ERI_Dxy_Py_Py_S_c;
  abcd[2366] = 4.0E0*I_ERI_Dxz_Py_Fy2z_S_cc-2.0E0*1*I_ERI_Dxz_Py_Py_S_c;
  abcd[2367] = 4.0E0*I_ERI_D2y_Py_Fy2z_S_cc-2.0E0*1*I_ERI_D2y_Py_Py_S_c;
  abcd[2368] = 4.0E0*I_ERI_Dyz_Py_Fy2z_S_cc-2.0E0*1*I_ERI_Dyz_Py_Py_S_c;
  abcd[2369] = 4.0E0*I_ERI_D2z_Py_Fy2z_S_cc-2.0E0*1*I_ERI_D2z_Py_Py_S_c;
  abcd[2370] = 4.0E0*I_ERI_D2x_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_D2x_Pz_Py_S_c;
  abcd[2371] = 4.0E0*I_ERI_Dxy_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Py_S_c;
  abcd[2372] = 4.0E0*I_ERI_Dxz_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Py_S_c;
  abcd[2373] = 4.0E0*I_ERI_D2y_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_D2y_Pz_Py_S_c;
  abcd[2374] = 4.0E0*I_ERI_Dyz_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Py_S_c;
  abcd[2375] = 4.0E0*I_ERI_D2z_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_D2z_Pz_Py_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_F_S_cc
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   ************************************************************/
  abcd[2376] = 4.0E0*I_ERI_D2x_Px_Fx2z_S_cc-2.0E0*1*I_ERI_D2x_Px_Px_S_c;
  abcd[2377] = 4.0E0*I_ERI_Dxy_Px_Fx2z_S_cc-2.0E0*1*I_ERI_Dxy_Px_Px_S_c;
  abcd[2378] = 4.0E0*I_ERI_Dxz_Px_Fx2z_S_cc-2.0E0*1*I_ERI_Dxz_Px_Px_S_c;
  abcd[2379] = 4.0E0*I_ERI_D2y_Px_Fx2z_S_cc-2.0E0*1*I_ERI_D2y_Px_Px_S_c;
  abcd[2380] = 4.0E0*I_ERI_Dyz_Px_Fx2z_S_cc-2.0E0*1*I_ERI_Dyz_Px_Px_S_c;
  abcd[2381] = 4.0E0*I_ERI_D2z_Px_Fx2z_S_cc-2.0E0*1*I_ERI_D2z_Px_Px_S_c;
  abcd[2382] = 4.0E0*I_ERI_D2x_Py_Fx2z_S_cc-2.0E0*1*I_ERI_D2x_Py_Px_S_c;
  abcd[2383] = 4.0E0*I_ERI_Dxy_Py_Fx2z_S_cc-2.0E0*1*I_ERI_Dxy_Py_Px_S_c;
  abcd[2384] = 4.0E0*I_ERI_Dxz_Py_Fx2z_S_cc-2.0E0*1*I_ERI_Dxz_Py_Px_S_c;
  abcd[2385] = 4.0E0*I_ERI_D2y_Py_Fx2z_S_cc-2.0E0*1*I_ERI_D2y_Py_Px_S_c;
  abcd[2386] = 4.0E0*I_ERI_Dyz_Py_Fx2z_S_cc-2.0E0*1*I_ERI_Dyz_Py_Px_S_c;
  abcd[2387] = 4.0E0*I_ERI_D2z_Py_Fx2z_S_cc-2.0E0*1*I_ERI_D2z_Py_Px_S_c;
  abcd[2388] = 4.0E0*I_ERI_D2x_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_D2x_Pz_Px_S_c;
  abcd[2389] = 4.0E0*I_ERI_Dxy_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Px_S_c;
  abcd[2390] = 4.0E0*I_ERI_Dxz_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Px_S_c;
  abcd[2391] = 4.0E0*I_ERI_D2y_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_D2y_Pz_Px_S_c;
  abcd[2392] = 4.0E0*I_ERI_Dyz_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Px_S_c;
  abcd[2393] = 4.0E0*I_ERI_D2z_Pz_Fx2z_S_cc-2.0E0*1*I_ERI_D2z_Pz_Px_S_c;
  abcd[2394] = 4.0E0*I_ERI_D2x_Px_Fy2z_S_cc-2.0E0*1*I_ERI_D2x_Px_Py_S_c;
  abcd[2395] = 4.0E0*I_ERI_Dxy_Px_Fy2z_S_cc-2.0E0*1*I_ERI_Dxy_Px_Py_S_c;
  abcd[2396] = 4.0E0*I_ERI_Dxz_Px_Fy2z_S_cc-2.0E0*1*I_ERI_Dxz_Px_Py_S_c;
  abcd[2397] = 4.0E0*I_ERI_D2y_Px_Fy2z_S_cc-2.0E0*1*I_ERI_D2y_Px_Py_S_c;
  abcd[2398] = 4.0E0*I_ERI_Dyz_Px_Fy2z_S_cc-2.0E0*1*I_ERI_Dyz_Px_Py_S_c;
  abcd[2399] = 4.0E0*I_ERI_D2z_Px_Fy2z_S_cc-2.0E0*1*I_ERI_D2z_Px_Py_S_c;
  abcd[2400] = 4.0E0*I_ERI_D2x_Py_Fy2z_S_cc-2.0E0*1*I_ERI_D2x_Py_Py_S_c;
  abcd[2401] = 4.0E0*I_ERI_Dxy_Py_Fy2z_S_cc-2.0E0*1*I_ERI_Dxy_Py_Py_S_c;
  abcd[2402] = 4.0E0*I_ERI_Dxz_Py_Fy2z_S_cc-2.0E0*1*I_ERI_Dxz_Py_Py_S_c;
  abcd[2403] = 4.0E0*I_ERI_D2y_Py_Fy2z_S_cc-2.0E0*1*I_ERI_D2y_Py_Py_S_c;
  abcd[2404] = 4.0E0*I_ERI_Dyz_Py_Fy2z_S_cc-2.0E0*1*I_ERI_Dyz_Py_Py_S_c;
  abcd[2405] = 4.0E0*I_ERI_D2z_Py_Fy2z_S_cc-2.0E0*1*I_ERI_D2z_Py_Py_S_c;
  abcd[2406] = 4.0E0*I_ERI_D2x_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_D2x_Pz_Py_S_c;
  abcd[2407] = 4.0E0*I_ERI_Dxy_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Py_S_c;
  abcd[2408] = 4.0E0*I_ERI_Dxz_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Py_S_c;
  abcd[2409] = 4.0E0*I_ERI_D2y_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_D2y_Pz_Py_S_c;
  abcd[2410] = 4.0E0*I_ERI_Dyz_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Py_S_c;
  abcd[2411] = 4.0E0*I_ERI_D2z_Pz_Fy2z_S_cc-2.0E0*1*I_ERI_D2z_Pz_Py_S_c;
  abcd[2412] = 4.0E0*I_ERI_D2x_Px_F3z_S_cc-2.0E0*1*I_ERI_D2x_Px_Pz_S_c-2.0E0*2*I_ERI_D2x_Px_Pz_S_c;
  abcd[2413] = 4.0E0*I_ERI_Dxy_Px_F3z_S_cc-2.0E0*1*I_ERI_Dxy_Px_Pz_S_c-2.0E0*2*I_ERI_Dxy_Px_Pz_S_c;
  abcd[2414] = 4.0E0*I_ERI_Dxz_Px_F3z_S_cc-2.0E0*1*I_ERI_Dxz_Px_Pz_S_c-2.0E0*2*I_ERI_Dxz_Px_Pz_S_c;
  abcd[2415] = 4.0E0*I_ERI_D2y_Px_F3z_S_cc-2.0E0*1*I_ERI_D2y_Px_Pz_S_c-2.0E0*2*I_ERI_D2y_Px_Pz_S_c;
  abcd[2416] = 4.0E0*I_ERI_Dyz_Px_F3z_S_cc-2.0E0*1*I_ERI_Dyz_Px_Pz_S_c-2.0E0*2*I_ERI_Dyz_Px_Pz_S_c;
  abcd[2417] = 4.0E0*I_ERI_D2z_Px_F3z_S_cc-2.0E0*1*I_ERI_D2z_Px_Pz_S_c-2.0E0*2*I_ERI_D2z_Px_Pz_S_c;
  abcd[2418] = 4.0E0*I_ERI_D2x_Py_F3z_S_cc-2.0E0*1*I_ERI_D2x_Py_Pz_S_c-2.0E0*2*I_ERI_D2x_Py_Pz_S_c;
  abcd[2419] = 4.0E0*I_ERI_Dxy_Py_F3z_S_cc-2.0E0*1*I_ERI_Dxy_Py_Pz_S_c-2.0E0*2*I_ERI_Dxy_Py_Pz_S_c;
  abcd[2420] = 4.0E0*I_ERI_Dxz_Py_F3z_S_cc-2.0E0*1*I_ERI_Dxz_Py_Pz_S_c-2.0E0*2*I_ERI_Dxz_Py_Pz_S_c;
  abcd[2421] = 4.0E0*I_ERI_D2y_Py_F3z_S_cc-2.0E0*1*I_ERI_D2y_Py_Pz_S_c-2.0E0*2*I_ERI_D2y_Py_Pz_S_c;
  abcd[2422] = 4.0E0*I_ERI_Dyz_Py_F3z_S_cc-2.0E0*1*I_ERI_Dyz_Py_Pz_S_c-2.0E0*2*I_ERI_Dyz_Py_Pz_S_c;
  abcd[2423] = 4.0E0*I_ERI_D2z_Py_F3z_S_cc-2.0E0*1*I_ERI_D2z_Py_Pz_S_c-2.0E0*2*I_ERI_D2z_Py_Pz_S_c;
  abcd[2424] = 4.0E0*I_ERI_D2x_Pz_F3z_S_cc-2.0E0*1*I_ERI_D2x_Pz_Pz_S_c-2.0E0*2*I_ERI_D2x_Pz_Pz_S_c;
  abcd[2425] = 4.0E0*I_ERI_Dxy_Pz_F3z_S_cc-2.0E0*1*I_ERI_Dxy_Pz_Pz_S_c-2.0E0*2*I_ERI_Dxy_Pz_Pz_S_c;
  abcd[2426] = 4.0E0*I_ERI_Dxz_Pz_F3z_S_cc-2.0E0*1*I_ERI_Dxz_Pz_Pz_S_c-2.0E0*2*I_ERI_Dxz_Pz_Pz_S_c;
  abcd[2427] = 4.0E0*I_ERI_D2y_Pz_F3z_S_cc-2.0E0*1*I_ERI_D2y_Pz_Pz_S_c-2.0E0*2*I_ERI_D2y_Pz_Pz_S_c;
  abcd[2428] = 4.0E0*I_ERI_Dyz_Pz_F3z_S_cc-2.0E0*1*I_ERI_Dyz_Pz_Pz_S_c-2.0E0*2*I_ERI_Dyz_Pz_Pz_S_c;
  abcd[2429] = 4.0E0*I_ERI_D2z_Pz_F3z_S_cc-2.0E0*1*I_ERI_D2z_Pz_Pz_S_c-2.0E0*2*I_ERI_D2z_Pz_Pz_S_c;
}
