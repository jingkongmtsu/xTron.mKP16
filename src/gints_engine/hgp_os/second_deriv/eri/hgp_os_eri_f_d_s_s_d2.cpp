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
// BRA1 as redundant position, total RHS integrals evaluated as: 62068
// BRA2 as redundant position, total RHS integrals evaluated as: 61848
// KET1 as redundant position, total RHS integrals evaluated as: 41473
// KET2 as redundant position, total RHS integrals evaluated as: 41473
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

void hgp_os_eri_f_d_s_s_d2(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
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
  Double I_ERI_H5x_S_S_S_a = 0.0E0;
  Double I_ERI_H4xy_S_S_S_a = 0.0E0;
  Double I_ERI_H4xz_S_S_S_a = 0.0E0;
  Double I_ERI_H3x2y_S_S_S_a = 0.0E0;
  Double I_ERI_H3xyz_S_S_S_a = 0.0E0;
  Double I_ERI_H3x2z_S_S_S_a = 0.0E0;
  Double I_ERI_H2x3y_S_S_S_a = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S_a = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S_a = 0.0E0;
  Double I_ERI_H2x3z_S_S_S_a = 0.0E0;
  Double I_ERI_Hx4y_S_S_S_a = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S_a = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S_a = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S_a = 0.0E0;
  Double I_ERI_Hx4z_S_S_S_a = 0.0E0;
  Double I_ERI_H5y_S_S_S_a = 0.0E0;
  Double I_ERI_H4yz_S_S_S_a = 0.0E0;
  Double I_ERI_H3y2z_S_S_S_a = 0.0E0;
  Double I_ERI_H2y3z_S_S_S_a = 0.0E0;
  Double I_ERI_Hy4z_S_S_S_a = 0.0E0;
  Double I_ERI_H5z_S_S_S_a = 0.0E0;
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
  Double I_ERI_D2x_S_S_S = 0.0E0;
  Double I_ERI_Dxy_S_S_S = 0.0E0;
  Double I_ERI_Dxz_S_S_S = 0.0E0;
  Double I_ERI_D2y_S_S_S = 0.0E0;
  Double I_ERI_Dyz_S_S_S = 0.0E0;
  Double I_ERI_D2z_S_S_S = 0.0E0;
  Double I_ERI_G4x_S_Px_S_c = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_c = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_c = 0.0E0;
  Double I_ERI_G4y_S_Px_S_c = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_c = 0.0E0;
  Double I_ERI_G4z_S_Px_S_c = 0.0E0;
  Double I_ERI_G4x_S_Py_S_c = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_c = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_c = 0.0E0;
  Double I_ERI_G4y_S_Py_S_c = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_c = 0.0E0;
  Double I_ERI_G4z_S_Py_S_c = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_c = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_c = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_c = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_c = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_c = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_c = 0.0E0;
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
  Double I_ERI_K7x_S_S_S_aa = 0.0E0;
  Double I_ERI_K6xy_S_S_S_aa = 0.0E0;
  Double I_ERI_K6xz_S_S_S_aa = 0.0E0;
  Double I_ERI_K5x2y_S_S_S_aa = 0.0E0;
  Double I_ERI_K5xyz_S_S_S_aa = 0.0E0;
  Double I_ERI_K5x2z_S_S_S_aa = 0.0E0;
  Double I_ERI_K4x3y_S_S_S_aa = 0.0E0;
  Double I_ERI_K4x2yz_S_S_S_aa = 0.0E0;
  Double I_ERI_K4xy2z_S_S_S_aa = 0.0E0;
  Double I_ERI_K4x3z_S_S_S_aa = 0.0E0;
  Double I_ERI_K3x4y_S_S_S_aa = 0.0E0;
  Double I_ERI_K3x3yz_S_S_S_aa = 0.0E0;
  Double I_ERI_K3x2y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_K3xy3z_S_S_S_aa = 0.0E0;
  Double I_ERI_K3x4z_S_S_S_aa = 0.0E0;
  Double I_ERI_K2x5y_S_S_S_aa = 0.0E0;
  Double I_ERI_K2x4yz_S_S_S_aa = 0.0E0;
  Double I_ERI_K2x3y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_K2x2y3z_S_S_S_aa = 0.0E0;
  Double I_ERI_K2xy4z_S_S_S_aa = 0.0E0;
  Double I_ERI_K2x5z_S_S_S_aa = 0.0E0;
  Double I_ERI_Kx6y_S_S_S_aa = 0.0E0;
  Double I_ERI_Kx5yz_S_S_S_aa = 0.0E0;
  Double I_ERI_Kx4y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_Kx3y3z_S_S_S_aa = 0.0E0;
  Double I_ERI_Kx2y4z_S_S_S_aa = 0.0E0;
  Double I_ERI_Kxy5z_S_S_S_aa = 0.0E0;
  Double I_ERI_Kx6z_S_S_S_aa = 0.0E0;
  Double I_ERI_K7y_S_S_S_aa = 0.0E0;
  Double I_ERI_K6yz_S_S_S_aa = 0.0E0;
  Double I_ERI_K5y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_K4y3z_S_S_S_aa = 0.0E0;
  Double I_ERI_K3y4z_S_S_S_aa = 0.0E0;
  Double I_ERI_K2y5z_S_S_S_aa = 0.0E0;
  Double I_ERI_Ky6z_S_S_S_aa = 0.0E0;
  Double I_ERI_K7z_S_S_S_aa = 0.0E0;
  Double I_ERI_I6x_S_S_S_aa = 0.0E0;
  Double I_ERI_I5xy_S_S_S_aa = 0.0E0;
  Double I_ERI_I5xz_S_S_S_aa = 0.0E0;
  Double I_ERI_I4x2y_S_S_S_aa = 0.0E0;
  Double I_ERI_I4xyz_S_S_S_aa = 0.0E0;
  Double I_ERI_I4x2z_S_S_S_aa = 0.0E0;
  Double I_ERI_I3x3y_S_S_S_aa = 0.0E0;
  Double I_ERI_I3x2yz_S_S_S_aa = 0.0E0;
  Double I_ERI_I3xy2z_S_S_S_aa = 0.0E0;
  Double I_ERI_I3x3z_S_S_S_aa = 0.0E0;
  Double I_ERI_I2x4y_S_S_S_aa = 0.0E0;
  Double I_ERI_I2x3yz_S_S_S_aa = 0.0E0;
  Double I_ERI_I2x2y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_I2xy3z_S_S_S_aa = 0.0E0;
  Double I_ERI_I2x4z_S_S_S_aa = 0.0E0;
  Double I_ERI_Ix5y_S_S_S_aa = 0.0E0;
  Double I_ERI_Ix4yz_S_S_S_aa = 0.0E0;
  Double I_ERI_Ix3y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_Ix2y3z_S_S_S_aa = 0.0E0;
  Double I_ERI_Ixy4z_S_S_S_aa = 0.0E0;
  Double I_ERI_Ix5z_S_S_S_aa = 0.0E0;
  Double I_ERI_I6y_S_S_S_aa = 0.0E0;
  Double I_ERI_I5yz_S_S_S_aa = 0.0E0;
  Double I_ERI_I4y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_I3y3z_S_S_S_aa = 0.0E0;
  Double I_ERI_I2y4z_S_S_S_aa = 0.0E0;
  Double I_ERI_Iy5z_S_S_S_aa = 0.0E0;
  Double I_ERI_I6z_S_S_S_aa = 0.0E0;
  Double I_ERI_H5x_S_S_S_aa = 0.0E0;
  Double I_ERI_H4xy_S_S_S_aa = 0.0E0;
  Double I_ERI_H4xz_S_S_S_aa = 0.0E0;
  Double I_ERI_H3x2y_S_S_S_aa = 0.0E0;
  Double I_ERI_H3xyz_S_S_S_aa = 0.0E0;
  Double I_ERI_H3x2z_S_S_S_aa = 0.0E0;
  Double I_ERI_H2x3y_S_S_S_aa = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S_aa = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S_aa = 0.0E0;
  Double I_ERI_H2x3z_S_S_S_aa = 0.0E0;
  Double I_ERI_Hx4y_S_S_S_aa = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S_aa = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S_aa = 0.0E0;
  Double I_ERI_Hx4z_S_S_S_aa = 0.0E0;
  Double I_ERI_H5y_S_S_S_aa = 0.0E0;
  Double I_ERI_H4yz_S_S_S_aa = 0.0E0;
  Double I_ERI_H3y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_H2y3z_S_S_S_aa = 0.0E0;
  Double I_ERI_Hy4z_S_S_S_aa = 0.0E0;
  Double I_ERI_H5z_S_S_S_aa = 0.0E0;
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
  Double I_ERI_Px_S_S_S = 0.0E0;
  Double I_ERI_Py_S_S_S = 0.0E0;
  Double I_ERI_Pz_S_S_S = 0.0E0;
  Double I_ERI_I6x_S_Px_S_ac = 0.0E0;
  Double I_ERI_I5xy_S_Px_S_ac = 0.0E0;
  Double I_ERI_I5xz_S_Px_S_ac = 0.0E0;
  Double I_ERI_I4x2y_S_Px_S_ac = 0.0E0;
  Double I_ERI_I4xyz_S_Px_S_ac = 0.0E0;
  Double I_ERI_I4x2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I3x3y_S_Px_S_ac = 0.0E0;
  Double I_ERI_I3x2yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_I3xy2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I3x3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I2x4y_S_Px_S_ac = 0.0E0;
  Double I_ERI_I2x3yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_I2x2y2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I2xy3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I2x4z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Ix5y_S_Px_S_ac = 0.0E0;
  Double I_ERI_Ix4yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_Ix3y2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Ix2y3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Ixy4z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Ix5z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I6y_S_Px_S_ac = 0.0E0;
  Double I_ERI_I5yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_I4y2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I3y3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I2y4z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Iy5z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I6z_S_Px_S_ac = 0.0E0;
  Double I_ERI_I6x_S_Py_S_ac = 0.0E0;
  Double I_ERI_I5xy_S_Py_S_ac = 0.0E0;
  Double I_ERI_I5xz_S_Py_S_ac = 0.0E0;
  Double I_ERI_I4x2y_S_Py_S_ac = 0.0E0;
  Double I_ERI_I4xyz_S_Py_S_ac = 0.0E0;
  Double I_ERI_I4x2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I3x3y_S_Py_S_ac = 0.0E0;
  Double I_ERI_I3x2yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_I3xy2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I3x3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I2x4y_S_Py_S_ac = 0.0E0;
  Double I_ERI_I2x3yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_I2x2y2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I2xy3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I2x4z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Ix5y_S_Py_S_ac = 0.0E0;
  Double I_ERI_Ix4yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_Ix3y2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Ix2y3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Ixy4z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Ix5z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I6y_S_Py_S_ac = 0.0E0;
  Double I_ERI_I5yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_I4y2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I3y3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I2y4z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Iy5z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I6z_S_Py_S_ac = 0.0E0;
  Double I_ERI_I6x_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I5xy_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I5xz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I4x2y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I4xyz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I4x2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I3x3y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I3x2yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I3xy2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I3x3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I2x4y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I2x3yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I2x2y2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I2xy3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I2x4z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Ix5y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Ix4yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Ix3y2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Ix2y3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Ixy4z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Ix5z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I6y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I5yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I4y2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I3y3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I2y4z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Iy5z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_I6z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H5x_S_Px_S_ac = 0.0E0;
  Double I_ERI_H4xy_S_Px_S_ac = 0.0E0;
  Double I_ERI_H4xz_S_Px_S_ac = 0.0E0;
  Double I_ERI_H3x2y_S_Px_S_ac = 0.0E0;
  Double I_ERI_H3xyz_S_Px_S_ac = 0.0E0;
  Double I_ERI_H3x2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_H2x3y_S_Px_S_ac = 0.0E0;
  Double I_ERI_H2x2yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_H2xy2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_H2x3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Hx4y_S_Px_S_ac = 0.0E0;
  Double I_ERI_Hx3yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_Hx2y2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Hxy3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Hx4z_S_Px_S_ac = 0.0E0;
  Double I_ERI_H5y_S_Px_S_ac = 0.0E0;
  Double I_ERI_H4yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_H3y2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_H2y3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Hy4z_S_Px_S_ac = 0.0E0;
  Double I_ERI_H5z_S_Px_S_ac = 0.0E0;
  Double I_ERI_H5x_S_Py_S_ac = 0.0E0;
  Double I_ERI_H4xy_S_Py_S_ac = 0.0E0;
  Double I_ERI_H4xz_S_Py_S_ac = 0.0E0;
  Double I_ERI_H3x2y_S_Py_S_ac = 0.0E0;
  Double I_ERI_H3xyz_S_Py_S_ac = 0.0E0;
  Double I_ERI_H3x2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_H2x3y_S_Py_S_ac = 0.0E0;
  Double I_ERI_H2x2yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_H2xy2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_H2x3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Hx4y_S_Py_S_ac = 0.0E0;
  Double I_ERI_Hx3yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_Hx2y2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Hxy3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Hx4z_S_Py_S_ac = 0.0E0;
  Double I_ERI_H5y_S_Py_S_ac = 0.0E0;
  Double I_ERI_H4yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_H3y2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_H2y3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Hy4z_S_Py_S_ac = 0.0E0;
  Double I_ERI_H5z_S_Py_S_ac = 0.0E0;
  Double I_ERI_H5x_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H4xy_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H4xz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H3x2y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H3xyz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H3x2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H2x3y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H2x2yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H2xy2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H2x3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Hx4y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Hx3yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Hx2y2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Hxy3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Hx4z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H5y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H4yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H3y2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H2y3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Hy4z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_H5z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G4x_S_Px_S_ac = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_ac = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_ac = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_ac = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_ac = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_G4y_S_Px_S_ac = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_G4z_S_Px_S_ac = 0.0E0;
  Double I_ERI_G4x_S_Py_S_ac = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_ac = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_ac = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_ac = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_ac = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_G4y_S_Py_S_ac = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_G4z_S_Py_S_ac = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_ac = 0.0E0;
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
  Double I_ERI_H5x_S_S_S_b = 0.0E0;
  Double I_ERI_H4xy_S_S_S_b = 0.0E0;
  Double I_ERI_H4xz_S_S_S_b = 0.0E0;
  Double I_ERI_H3x2y_S_S_S_b = 0.0E0;
  Double I_ERI_H3xyz_S_S_S_b = 0.0E0;
  Double I_ERI_H3x2z_S_S_S_b = 0.0E0;
  Double I_ERI_H2x3y_S_S_S_b = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S_b = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S_b = 0.0E0;
  Double I_ERI_H2x3z_S_S_S_b = 0.0E0;
  Double I_ERI_Hx4y_S_S_S_b = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S_b = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S_b = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S_b = 0.0E0;
  Double I_ERI_Hx4z_S_S_S_b = 0.0E0;
  Double I_ERI_H5y_S_S_S_b = 0.0E0;
  Double I_ERI_H4yz_S_S_S_b = 0.0E0;
  Double I_ERI_H3y2z_S_S_S_b = 0.0E0;
  Double I_ERI_H2y3z_S_S_S_b = 0.0E0;
  Double I_ERI_Hy4z_S_S_S_b = 0.0E0;
  Double I_ERI_H5z_S_S_S_b = 0.0E0;
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
  Double I_ERI_H5x_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H4xy_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H4xz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H3x2y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H3xyz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H3x2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H2x3y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H2x2yz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H2xy2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H2x3z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Hx4y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Hx3yz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Hxy3z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Hx4z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H5y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H4yz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H3y2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H2y3z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Hy4z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H5z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_H5x_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H4xy_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H4xz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H3x2y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H3xyz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H3x2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H2x3y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H2x3z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Hx4y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Hx4z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H5y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H4yz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H3y2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H2y3z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Hy4z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H5z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_H5x_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H4xy_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H4xz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H3x2y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H3xyz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H3x2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H2x3y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H2x3z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Hx4y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Hx4z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H5y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H4yz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H3y2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H2y3z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Hy4z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H5z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_H5x_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H4xy_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H4xz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H3x2y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H3xyz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H3x2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H2x3y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H2x2yz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H2xy2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H2x3z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Hx4y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Hx3yz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Hxy3z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Hx4z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H5y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H4yz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H3y2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H2y3z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Hy4z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H5z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_H5x_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H4xy_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H4xz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H3x2y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H3xyz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H3x2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H2x3y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H2x2yz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H2xy2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H2x3z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Hx4y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Hx3yz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Hxy3z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Hx4z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H5y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H4yz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H3y2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H2y3z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Hy4z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H5z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_H5x_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H4xy_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H4xz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H3x2y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H3xyz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H3x2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H2x3y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H2x2yz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H2xy2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H2x3z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Hx4y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Hx3yz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Hxy3z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Hx4z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H5y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H4yz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H3y2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H2y3z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Hy4z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H5z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_cc = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_cc = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_H5x_S_S_S_c = 0.0E0;
  Double I_ERI_H4xy_S_S_S_c = 0.0E0;
  Double I_ERI_H4xz_S_S_S_c = 0.0E0;
  Double I_ERI_H3x2y_S_S_S_c = 0.0E0;
  Double I_ERI_H3xyz_S_S_S_c = 0.0E0;
  Double I_ERI_H3x2z_S_S_S_c = 0.0E0;
  Double I_ERI_H2x3y_S_S_S_c = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S_c = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S_c = 0.0E0;
  Double I_ERI_H2x3z_S_S_S_c = 0.0E0;
  Double I_ERI_Hx4y_S_S_S_c = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S_c = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S_c = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S_c = 0.0E0;
  Double I_ERI_Hx4z_S_S_S_c = 0.0E0;
  Double I_ERI_H5y_S_S_S_c = 0.0E0;
  Double I_ERI_H4yz_S_S_S_c = 0.0E0;
  Double I_ERI_H3y2z_S_S_S_c = 0.0E0;
  Double I_ERI_H2y3z_S_S_S_c = 0.0E0;
  Double I_ERI_Hy4z_S_S_S_c = 0.0E0;
  Double I_ERI_H5z_S_S_S_c = 0.0E0;
  Double I_ERI_G4x_S_S_S_c = 0.0E0;
  Double I_ERI_G3xy_S_S_S_c = 0.0E0;
  Double I_ERI_G3xz_S_S_S_c = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_c = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_c = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_c = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_c = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_c = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_c = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_c = 0.0E0;
  Double I_ERI_G4y_S_S_S_c = 0.0E0;
  Double I_ERI_G3yz_S_S_S_c = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_c = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_c = 0.0E0;
  Double I_ERI_G4z_S_S_S_c = 0.0E0;
  Double I_ERI_F3x_S_S_S_c = 0.0E0;
  Double I_ERI_F2xy_S_S_S_c = 0.0E0;
  Double I_ERI_F2xz_S_S_S_c = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_c = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_c = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_c = 0.0E0;
  Double I_ERI_F3y_S_S_S_c = 0.0E0;
  Double I_ERI_F2yz_S_S_S_c = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_c = 0.0E0;
  Double I_ERI_F3z_S_S_S_c = 0.0E0;
  Double I_ERI_K7x_S_S_S_ab = 0.0E0;
  Double I_ERI_K6xy_S_S_S_ab = 0.0E0;
  Double I_ERI_K6xz_S_S_S_ab = 0.0E0;
  Double I_ERI_K5x2y_S_S_S_ab = 0.0E0;
  Double I_ERI_K5xyz_S_S_S_ab = 0.0E0;
  Double I_ERI_K5x2z_S_S_S_ab = 0.0E0;
  Double I_ERI_K4x3y_S_S_S_ab = 0.0E0;
  Double I_ERI_K4x2yz_S_S_S_ab = 0.0E0;
  Double I_ERI_K4xy2z_S_S_S_ab = 0.0E0;
  Double I_ERI_K4x3z_S_S_S_ab = 0.0E0;
  Double I_ERI_K3x4y_S_S_S_ab = 0.0E0;
  Double I_ERI_K3x3yz_S_S_S_ab = 0.0E0;
  Double I_ERI_K3x2y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_K3xy3z_S_S_S_ab = 0.0E0;
  Double I_ERI_K3x4z_S_S_S_ab = 0.0E0;
  Double I_ERI_K2x5y_S_S_S_ab = 0.0E0;
  Double I_ERI_K2x4yz_S_S_S_ab = 0.0E0;
  Double I_ERI_K2x3y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_K2x2y3z_S_S_S_ab = 0.0E0;
  Double I_ERI_K2xy4z_S_S_S_ab = 0.0E0;
  Double I_ERI_K2x5z_S_S_S_ab = 0.0E0;
  Double I_ERI_Kx6y_S_S_S_ab = 0.0E0;
  Double I_ERI_Kx5yz_S_S_S_ab = 0.0E0;
  Double I_ERI_Kx4y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_Kx3y3z_S_S_S_ab = 0.0E0;
  Double I_ERI_Kx2y4z_S_S_S_ab = 0.0E0;
  Double I_ERI_Kxy5z_S_S_S_ab = 0.0E0;
  Double I_ERI_Kx6z_S_S_S_ab = 0.0E0;
  Double I_ERI_K7y_S_S_S_ab = 0.0E0;
  Double I_ERI_K6yz_S_S_S_ab = 0.0E0;
  Double I_ERI_K5y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_K4y3z_S_S_S_ab = 0.0E0;
  Double I_ERI_K3y4z_S_S_S_ab = 0.0E0;
  Double I_ERI_K2y5z_S_S_S_ab = 0.0E0;
  Double I_ERI_Ky6z_S_S_S_ab = 0.0E0;
  Double I_ERI_K7z_S_S_S_ab = 0.0E0;
  Double I_ERI_I6x_S_S_S_ab = 0.0E0;
  Double I_ERI_I5xy_S_S_S_ab = 0.0E0;
  Double I_ERI_I5xz_S_S_S_ab = 0.0E0;
  Double I_ERI_I4x2y_S_S_S_ab = 0.0E0;
  Double I_ERI_I4xyz_S_S_S_ab = 0.0E0;
  Double I_ERI_I4x2z_S_S_S_ab = 0.0E0;
  Double I_ERI_I3x3y_S_S_S_ab = 0.0E0;
  Double I_ERI_I3x2yz_S_S_S_ab = 0.0E0;
  Double I_ERI_I3xy2z_S_S_S_ab = 0.0E0;
  Double I_ERI_I3x3z_S_S_S_ab = 0.0E0;
  Double I_ERI_I2x4y_S_S_S_ab = 0.0E0;
  Double I_ERI_I2x3yz_S_S_S_ab = 0.0E0;
  Double I_ERI_I2x2y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_I2xy3z_S_S_S_ab = 0.0E0;
  Double I_ERI_I2x4z_S_S_S_ab = 0.0E0;
  Double I_ERI_Ix5y_S_S_S_ab = 0.0E0;
  Double I_ERI_Ix4yz_S_S_S_ab = 0.0E0;
  Double I_ERI_Ix3y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_Ix2y3z_S_S_S_ab = 0.0E0;
  Double I_ERI_Ixy4z_S_S_S_ab = 0.0E0;
  Double I_ERI_Ix5z_S_S_S_ab = 0.0E0;
  Double I_ERI_I6y_S_S_S_ab = 0.0E0;
  Double I_ERI_I5yz_S_S_S_ab = 0.0E0;
  Double I_ERI_I4y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_I3y3z_S_S_S_ab = 0.0E0;
  Double I_ERI_I2y4z_S_S_S_ab = 0.0E0;
  Double I_ERI_Iy5z_S_S_S_ab = 0.0E0;
  Double I_ERI_I6z_S_S_S_ab = 0.0E0;
  Double I_ERI_H5x_S_S_S_ab = 0.0E0;
  Double I_ERI_H4xy_S_S_S_ab = 0.0E0;
  Double I_ERI_H4xz_S_S_S_ab = 0.0E0;
  Double I_ERI_H3x2y_S_S_S_ab = 0.0E0;
  Double I_ERI_H3xyz_S_S_S_ab = 0.0E0;
  Double I_ERI_H3x2z_S_S_S_ab = 0.0E0;
  Double I_ERI_H2x3y_S_S_S_ab = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S_ab = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S_ab = 0.0E0;
  Double I_ERI_H2x3z_S_S_S_ab = 0.0E0;
  Double I_ERI_Hx4y_S_S_S_ab = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S_ab = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S_ab = 0.0E0;
  Double I_ERI_Hx4z_S_S_S_ab = 0.0E0;
  Double I_ERI_H5y_S_S_S_ab = 0.0E0;
  Double I_ERI_H4yz_S_S_S_ab = 0.0E0;
  Double I_ERI_H3y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_H2y3z_S_S_S_ab = 0.0E0;
  Double I_ERI_Hy4z_S_S_S_ab = 0.0E0;
  Double I_ERI_H5z_S_S_S_ab = 0.0E0;
  Double I_ERI_G4x_S_S_S_ab = 0.0E0;
  Double I_ERI_G3xy_S_S_S_ab = 0.0E0;
  Double I_ERI_G3xz_S_S_S_ab = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_ab = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_ab = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_ab = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_ab = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_ab = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_ab = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_ab = 0.0E0;
  Double I_ERI_G4y_S_S_S_ab = 0.0E0;
  Double I_ERI_G3yz_S_S_S_ab = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_ab = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_ab = 0.0E0;
  Double I_ERI_G4z_S_S_S_ab = 0.0E0;
  Double I_ERI_D2x_S_S_S_b = 0.0E0;
  Double I_ERI_Dxy_S_S_S_b = 0.0E0;
  Double I_ERI_Dxz_S_S_S_b = 0.0E0;
  Double I_ERI_D2y_S_S_S_b = 0.0E0;
  Double I_ERI_Dyz_S_S_S_b = 0.0E0;
  Double I_ERI_D2z_S_S_S_b = 0.0E0;
  Double I_ERI_I6x_S_Px_S_bc = 0.0E0;
  Double I_ERI_I5xy_S_Px_S_bc = 0.0E0;
  Double I_ERI_I5xz_S_Px_S_bc = 0.0E0;
  Double I_ERI_I4x2y_S_Px_S_bc = 0.0E0;
  Double I_ERI_I4xyz_S_Px_S_bc = 0.0E0;
  Double I_ERI_I4x2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I3x3y_S_Px_S_bc = 0.0E0;
  Double I_ERI_I3x2yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_I3xy2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I3x3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I2x4y_S_Px_S_bc = 0.0E0;
  Double I_ERI_I2x3yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_I2x2y2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I2xy3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I2x4z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Ix5y_S_Px_S_bc = 0.0E0;
  Double I_ERI_Ix4yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_Ix3y2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Ix2y3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Ixy4z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Ix5z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I6y_S_Px_S_bc = 0.0E0;
  Double I_ERI_I5yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_I4y2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I3y3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I2y4z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Iy5z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I6z_S_Px_S_bc = 0.0E0;
  Double I_ERI_I6x_S_Py_S_bc = 0.0E0;
  Double I_ERI_I5xy_S_Py_S_bc = 0.0E0;
  Double I_ERI_I5xz_S_Py_S_bc = 0.0E0;
  Double I_ERI_I4x2y_S_Py_S_bc = 0.0E0;
  Double I_ERI_I4xyz_S_Py_S_bc = 0.0E0;
  Double I_ERI_I4x2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I3x3y_S_Py_S_bc = 0.0E0;
  Double I_ERI_I3x2yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_I3xy2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I3x3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I2x4y_S_Py_S_bc = 0.0E0;
  Double I_ERI_I2x3yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_I2x2y2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I2xy3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I2x4z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Ix5y_S_Py_S_bc = 0.0E0;
  Double I_ERI_Ix4yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_Ix3y2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Ix2y3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Ixy4z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Ix5z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I6y_S_Py_S_bc = 0.0E0;
  Double I_ERI_I5yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_I4y2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I3y3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I2y4z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Iy5z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I6z_S_Py_S_bc = 0.0E0;
  Double I_ERI_I6x_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I5xy_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I5xz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I4x2y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I4xyz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I4x2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I3x3y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I3x2yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I3xy2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I3x3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I2x4y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I2x3yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I2x2y2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I2xy3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I2x4z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Ix5y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Ix4yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Ix3y2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Ix2y3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Ixy4z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Ix5z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I6y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I5yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I4y2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I3y3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I2y4z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Iy5z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_I6z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H5x_S_Px_S_bc = 0.0E0;
  Double I_ERI_H4xy_S_Px_S_bc = 0.0E0;
  Double I_ERI_H4xz_S_Px_S_bc = 0.0E0;
  Double I_ERI_H3x2y_S_Px_S_bc = 0.0E0;
  Double I_ERI_H3xyz_S_Px_S_bc = 0.0E0;
  Double I_ERI_H3x2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_H2x3y_S_Px_S_bc = 0.0E0;
  Double I_ERI_H2x2yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_H2xy2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_H2x3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Hx4y_S_Px_S_bc = 0.0E0;
  Double I_ERI_Hx3yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_Hx2y2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Hxy3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Hx4z_S_Px_S_bc = 0.0E0;
  Double I_ERI_H5y_S_Px_S_bc = 0.0E0;
  Double I_ERI_H4yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_H3y2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_H2y3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Hy4z_S_Px_S_bc = 0.0E0;
  Double I_ERI_H5z_S_Px_S_bc = 0.0E0;
  Double I_ERI_H5x_S_Py_S_bc = 0.0E0;
  Double I_ERI_H4xy_S_Py_S_bc = 0.0E0;
  Double I_ERI_H4xz_S_Py_S_bc = 0.0E0;
  Double I_ERI_H3x2y_S_Py_S_bc = 0.0E0;
  Double I_ERI_H3xyz_S_Py_S_bc = 0.0E0;
  Double I_ERI_H3x2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_H2x3y_S_Py_S_bc = 0.0E0;
  Double I_ERI_H2x2yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_H2xy2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_H2x3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Hx4y_S_Py_S_bc = 0.0E0;
  Double I_ERI_Hx3yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_Hx2y2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Hxy3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Hx4z_S_Py_S_bc = 0.0E0;
  Double I_ERI_H5y_S_Py_S_bc = 0.0E0;
  Double I_ERI_H4yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_H3y2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_H2y3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Hy4z_S_Py_S_bc = 0.0E0;
  Double I_ERI_H5z_S_Py_S_bc = 0.0E0;
  Double I_ERI_H5x_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H4xy_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H4xz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H3x2y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H3xyz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H3x2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H2x3y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H2x2yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H2xy2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H2x3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Hx4y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Hx3yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Hx2y2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Hxy3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Hx4z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H5y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H4yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H3y2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H2y3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Hy4z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_H5z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Px_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Px_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Px_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Px_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Px_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Py_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Py_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Py_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Py_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Py_S_bc = 0.0E0;
  Double I_ERI_G4x_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G4y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_G4z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Px_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Px_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Px_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Py_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Py_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Py_S_bc = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_bc = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_bc = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_bc = 0.0E0;
  Double I_ERI_K7x_S_S_S_bb = 0.0E0;
  Double I_ERI_K6xy_S_S_S_bb = 0.0E0;
  Double I_ERI_K6xz_S_S_S_bb = 0.0E0;
  Double I_ERI_K5x2y_S_S_S_bb = 0.0E0;
  Double I_ERI_K5xyz_S_S_S_bb = 0.0E0;
  Double I_ERI_K5x2z_S_S_S_bb = 0.0E0;
  Double I_ERI_K4x3y_S_S_S_bb = 0.0E0;
  Double I_ERI_K4x2yz_S_S_S_bb = 0.0E0;
  Double I_ERI_K4xy2z_S_S_S_bb = 0.0E0;
  Double I_ERI_K4x3z_S_S_S_bb = 0.0E0;
  Double I_ERI_K3x4y_S_S_S_bb = 0.0E0;
  Double I_ERI_K3x3yz_S_S_S_bb = 0.0E0;
  Double I_ERI_K3x2y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_K3xy3z_S_S_S_bb = 0.0E0;
  Double I_ERI_K3x4z_S_S_S_bb = 0.0E0;
  Double I_ERI_K2x5y_S_S_S_bb = 0.0E0;
  Double I_ERI_K2x4yz_S_S_S_bb = 0.0E0;
  Double I_ERI_K2x3y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_K2x2y3z_S_S_S_bb = 0.0E0;
  Double I_ERI_K2xy4z_S_S_S_bb = 0.0E0;
  Double I_ERI_K2x5z_S_S_S_bb = 0.0E0;
  Double I_ERI_Kx6y_S_S_S_bb = 0.0E0;
  Double I_ERI_Kx5yz_S_S_S_bb = 0.0E0;
  Double I_ERI_Kx4y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_Kx3y3z_S_S_S_bb = 0.0E0;
  Double I_ERI_Kx2y4z_S_S_S_bb = 0.0E0;
  Double I_ERI_Kxy5z_S_S_S_bb = 0.0E0;
  Double I_ERI_Kx6z_S_S_S_bb = 0.0E0;
  Double I_ERI_K7y_S_S_S_bb = 0.0E0;
  Double I_ERI_K6yz_S_S_S_bb = 0.0E0;
  Double I_ERI_K5y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_K4y3z_S_S_S_bb = 0.0E0;
  Double I_ERI_K3y4z_S_S_S_bb = 0.0E0;
  Double I_ERI_K2y5z_S_S_S_bb = 0.0E0;
  Double I_ERI_Ky6z_S_S_S_bb = 0.0E0;
  Double I_ERI_K7z_S_S_S_bb = 0.0E0;
  Double I_ERI_I6x_S_S_S_bb = 0.0E0;
  Double I_ERI_I5xy_S_S_S_bb = 0.0E0;
  Double I_ERI_I5xz_S_S_S_bb = 0.0E0;
  Double I_ERI_I4x2y_S_S_S_bb = 0.0E0;
  Double I_ERI_I4xyz_S_S_S_bb = 0.0E0;
  Double I_ERI_I4x2z_S_S_S_bb = 0.0E0;
  Double I_ERI_I3x3y_S_S_S_bb = 0.0E0;
  Double I_ERI_I3x2yz_S_S_S_bb = 0.0E0;
  Double I_ERI_I3xy2z_S_S_S_bb = 0.0E0;
  Double I_ERI_I3x3z_S_S_S_bb = 0.0E0;
  Double I_ERI_I2x4y_S_S_S_bb = 0.0E0;
  Double I_ERI_I2x3yz_S_S_S_bb = 0.0E0;
  Double I_ERI_I2x2y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_I2xy3z_S_S_S_bb = 0.0E0;
  Double I_ERI_I2x4z_S_S_S_bb = 0.0E0;
  Double I_ERI_Ix5y_S_S_S_bb = 0.0E0;
  Double I_ERI_Ix4yz_S_S_S_bb = 0.0E0;
  Double I_ERI_Ix3y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_Ix2y3z_S_S_S_bb = 0.0E0;
  Double I_ERI_Ixy4z_S_S_S_bb = 0.0E0;
  Double I_ERI_Ix5z_S_S_S_bb = 0.0E0;
  Double I_ERI_I6y_S_S_S_bb = 0.0E0;
  Double I_ERI_I5yz_S_S_S_bb = 0.0E0;
  Double I_ERI_I4y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_I3y3z_S_S_S_bb = 0.0E0;
  Double I_ERI_I2y4z_S_S_S_bb = 0.0E0;
  Double I_ERI_Iy5z_S_S_S_bb = 0.0E0;
  Double I_ERI_I6z_S_S_S_bb = 0.0E0;
  Double I_ERI_H5x_S_S_S_bb = 0.0E0;
  Double I_ERI_H4xy_S_S_S_bb = 0.0E0;
  Double I_ERI_H4xz_S_S_S_bb = 0.0E0;
  Double I_ERI_H3x2y_S_S_S_bb = 0.0E0;
  Double I_ERI_H3xyz_S_S_S_bb = 0.0E0;
  Double I_ERI_H3x2z_S_S_S_bb = 0.0E0;
  Double I_ERI_H2x3y_S_S_S_bb = 0.0E0;
  Double I_ERI_H2x2yz_S_S_S_bb = 0.0E0;
  Double I_ERI_H2xy2z_S_S_S_bb = 0.0E0;
  Double I_ERI_H2x3z_S_S_S_bb = 0.0E0;
  Double I_ERI_Hx4y_S_S_S_bb = 0.0E0;
  Double I_ERI_Hx3yz_S_S_S_bb = 0.0E0;
  Double I_ERI_Hx2y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_Hxy3z_S_S_S_bb = 0.0E0;
  Double I_ERI_Hx4z_S_S_S_bb = 0.0E0;
  Double I_ERI_H5y_S_S_S_bb = 0.0E0;
  Double I_ERI_H4yz_S_S_S_bb = 0.0E0;
  Double I_ERI_H3y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_H2y3z_S_S_S_bb = 0.0E0;
  Double I_ERI_Hy4z_S_S_S_bb = 0.0E0;
  Double I_ERI_H5z_S_S_S_bb = 0.0E0;
  Double I_ERI_G4x_S_S_S_bb = 0.0E0;
  Double I_ERI_G3xy_S_S_S_bb = 0.0E0;
  Double I_ERI_G3xz_S_S_S_bb = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_bb = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_bb = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_bb = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_bb = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_bb = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_bb = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_bb = 0.0E0;
  Double I_ERI_G4y_S_S_S_bb = 0.0E0;
  Double I_ERI_G3yz_S_S_S_bb = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_bb = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_bb = 0.0E0;
  Double I_ERI_G4z_S_S_S_bb = 0.0E0;
  Double I_ERI_F3x_S_S_S_bb = 0.0E0;
  Double I_ERI_F2xy_S_S_S_bb = 0.0E0;
  Double I_ERI_F2xz_S_S_S_bb = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_bb = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_bb = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_bb = 0.0E0;
  Double I_ERI_F3y_S_S_S_bb = 0.0E0;
  Double I_ERI_F2yz_S_S_S_bb = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_bb = 0.0E0;
  Double I_ERI_F3z_S_S_S_bb = 0.0E0;

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
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M4_vrr = PAX*I_ERI_Px_S_S_S_M4_vrr+WPX*I_ERI_Px_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_D2y_S_S_S_M4_vrr = PAY*I_ERI_Py_S_S_S_M4_vrr+WPY*I_ERI_Py_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_D2z_S_S_S_M4_vrr = PAZ*I_ERI_Pz_S_S_S_M4_vrr+WPZ*I_ERI_Pz_S_S_S_M5_vrr+oned2z*I_ERI_S_S_S_S_M4_vrr-rhod2zsq*I_ERI_S_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M4_vrr = PAX*I_ERI_D2x_S_S_S_M4_vrr+WPX*I_ERI_D2x_S_S_S_M5_vrr+2*oned2z*I_ERI_Px_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_F2xy_S_S_S_M4_vrr = PAY*I_ERI_D2x_S_S_S_M4_vrr+WPY*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_F2xz_S_S_S_M4_vrr = PAZ*I_ERI_D2x_S_S_S_M4_vrr+WPZ*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_F3y_S_S_S_M4_vrr = PAY*I_ERI_D2y_S_S_S_M4_vrr+WPY*I_ERI_D2y_S_S_S_M5_vrr+2*oned2z*I_ERI_Py_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_F2yz_S_S_S_M4_vrr = PAZ*I_ERI_D2y_S_S_S_M4_vrr+WPZ*I_ERI_D2y_S_S_S_M5_vrr;
      Double I_ERI_F3z_S_S_S_M4_vrr = PAZ*I_ERI_D2z_S_S_S_M4_vrr+WPZ*I_ERI_D2z_S_S_S_M5_vrr+2*oned2z*I_ERI_Pz_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_F_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M3_vrr = PAX*I_ERI_D2x_S_S_S_M3_vrr+WPX*I_ERI_D2x_S_S_S_M4_vrr+2*oned2z*I_ERI_Px_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_F2xy_S_S_S_M3_vrr = PAY*I_ERI_D2x_S_S_S_M3_vrr+WPY*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F2xz_S_S_S_M3_vrr = PAZ*I_ERI_D2x_S_S_S_M3_vrr+WPZ*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Fx2y_S_S_S_M3_vrr = PAX*I_ERI_D2y_S_S_S_M3_vrr+WPX*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Fx2z_S_S_S_M3_vrr = PAX*I_ERI_D2z_S_S_S_M3_vrr+WPX*I_ERI_D2z_S_S_S_M4_vrr;
      Double I_ERI_F3y_S_S_S_M3_vrr = PAY*I_ERI_D2y_S_S_S_M3_vrr+WPY*I_ERI_D2y_S_S_S_M4_vrr+2*oned2z*I_ERI_Py_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_F2yz_S_S_S_M3_vrr = PAZ*I_ERI_D2y_S_S_S_M3_vrr+WPZ*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_F3z_S_S_S_M3_vrr = PAZ*I_ERI_D2z_S_S_S_M3_vrr+WPZ*I_ERI_D2z_S_S_S_M4_vrr+2*oned2z*I_ERI_Pz_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M3_vrr = PAX*I_ERI_F3x_S_S_S_M3_vrr+WPX*I_ERI_F3x_S_S_S_M4_vrr+3*oned2z*I_ERI_D2x_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G3xy_S_S_S_M3_vrr = PAY*I_ERI_F3x_S_S_S_M3_vrr+WPY*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G3xz_S_S_S_M3_vrr = PAZ*I_ERI_F3x_S_S_S_M3_vrr+WPZ*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G2x2y_S_S_S_M3_vrr = PAY*I_ERI_F2xy_S_S_S_M3_vrr+WPY*I_ERI_F2xy_S_S_S_M4_vrr+oned2z*I_ERI_D2x_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G2x2z_S_S_S_M3_vrr = PAZ*I_ERI_F2xz_S_S_S_M3_vrr+WPZ*I_ERI_F2xz_S_S_S_M4_vrr+oned2z*I_ERI_D2x_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Gx3y_S_S_S_M3_vrr = PAX*I_ERI_F3y_S_S_S_M3_vrr+WPX*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_Gx3z_S_S_S_M3_vrr = PAX*I_ERI_F3z_S_S_S_M3_vrr+WPX*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4y_S_S_S_M3_vrr = PAY*I_ERI_F3y_S_S_S_M3_vrr+WPY*I_ERI_F3y_S_S_S_M4_vrr+3*oned2z*I_ERI_D2y_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_G3yz_S_S_S_M3_vrr = PAZ*I_ERI_F3y_S_S_S_M3_vrr+WPZ*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_G2y2z_S_S_S_M3_vrr = PAZ*I_ERI_F2yz_S_S_S_M3_vrr+WPZ*I_ERI_F2yz_S_S_S_M4_vrr+oned2z*I_ERI_D2y_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Gy3z_S_S_S_M3_vrr = PAY*I_ERI_F3z_S_S_S_M3_vrr+WPY*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4z_S_S_S_M3_vrr = PAZ*I_ERI_F3z_S_S_S_M3_vrr+WPZ*I_ERI_F3z_S_S_S_M4_vrr+3*oned2z*I_ERI_D2z_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_H_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       ************************************************************/
      Double I_ERI_H5x_S_S_S_M2_vrr = PAX*I_ERI_G4x_S_S_S_M2_vrr+WPX*I_ERI_G4x_S_S_S_M3_vrr+4*oned2z*I_ERI_F3x_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H4xy_S_S_S_M2_vrr = PAY*I_ERI_G4x_S_S_S_M2_vrr+WPY*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_H4xz_S_S_S_M2_vrr = PAZ*I_ERI_G4x_S_S_S_M2_vrr+WPZ*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_H3x2y_S_S_S_M2_vrr = PAY*I_ERI_G3xy_S_S_S_M2_vrr+WPY*I_ERI_G3xy_S_S_S_M3_vrr+oned2z*I_ERI_F3x_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H3xyz_S_S_S_M2_vrr = PAZ*I_ERI_G3xy_S_S_S_M2_vrr+WPZ*I_ERI_G3xy_S_S_S_M3_vrr;
      Double I_ERI_H3x2z_S_S_S_M2_vrr = PAZ*I_ERI_G3xz_S_S_S_M2_vrr+WPZ*I_ERI_G3xz_S_S_S_M3_vrr+oned2z*I_ERI_F3x_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H2x3y_S_S_S_M2_vrr = PAX*I_ERI_Gx3y_S_S_S_M2_vrr+WPX*I_ERI_Gx3y_S_S_S_M3_vrr+oned2z*I_ERI_F3y_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H2x2yz_S_S_S_M2_vrr = PAZ*I_ERI_G2x2y_S_S_S_M2_vrr+WPZ*I_ERI_G2x2y_S_S_S_M3_vrr;
      Double I_ERI_H2xy2z_S_S_S_M2_vrr = PAY*I_ERI_G2x2z_S_S_S_M2_vrr+WPY*I_ERI_G2x2z_S_S_S_M3_vrr;
      Double I_ERI_H2x3z_S_S_S_M2_vrr = PAX*I_ERI_Gx3z_S_S_S_M2_vrr+WPX*I_ERI_Gx3z_S_S_S_M3_vrr+oned2z*I_ERI_F3z_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_Hx4y_S_S_S_M2_vrr = PAX*I_ERI_G4y_S_S_S_M2_vrr+WPX*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_Hx3yz_S_S_S_M2_vrr = PAZ*I_ERI_Gx3y_S_S_S_M2_vrr+WPZ*I_ERI_Gx3y_S_S_S_M3_vrr;
      Double I_ERI_Hx2y2z_S_S_S_M2_vrr = PAX*I_ERI_G2y2z_S_S_S_M2_vrr+WPX*I_ERI_G2y2z_S_S_S_M3_vrr;
      Double I_ERI_Hxy3z_S_S_S_M2_vrr = PAY*I_ERI_Gx3z_S_S_S_M2_vrr+WPY*I_ERI_Gx3z_S_S_S_M3_vrr;
      Double I_ERI_Hx4z_S_S_S_M2_vrr = PAX*I_ERI_G4z_S_S_S_M2_vrr+WPX*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_H5y_S_S_S_M2_vrr = PAY*I_ERI_G4y_S_S_S_M2_vrr+WPY*I_ERI_G4y_S_S_S_M3_vrr+4*oned2z*I_ERI_F3y_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H4yz_S_S_S_M2_vrr = PAZ*I_ERI_G4y_S_S_S_M2_vrr+WPZ*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_H3y2z_S_S_S_M2_vrr = PAZ*I_ERI_G3yz_S_S_S_M2_vrr+WPZ*I_ERI_G3yz_S_S_S_M3_vrr+oned2z*I_ERI_F3y_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H2y3z_S_S_S_M2_vrr = PAY*I_ERI_Gy3z_S_S_S_M2_vrr+WPY*I_ERI_Gy3z_S_S_S_M3_vrr+oned2z*I_ERI_F3z_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_Hy4z_S_S_S_M2_vrr = PAY*I_ERI_G4z_S_S_S_M2_vrr+WPY*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_H5z_S_S_S_M2_vrr = PAZ*I_ERI_G4z_S_S_S_M2_vrr+WPZ*I_ERI_G4z_S_S_S_M3_vrr+4*oned2z*I_ERI_F3z_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_G_S_P_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_Px_S_M1_vrr = QCX*I_ERI_G4x_S_S_S_M1_vrr+WQX*I_ERI_G4x_S_S_S_M2_vrr+4*oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_Px_S_M1_vrr = QCX*I_ERI_G3xy_S_S_S_M1_vrr+WQX*I_ERI_G3xy_S_S_S_M2_vrr+3*oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_Px_S_M1_vrr = QCX*I_ERI_G3xz_S_S_S_M1_vrr+WQX*I_ERI_G3xz_S_S_S_M2_vrr+3*oned2k*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_Px_S_M1_vrr = QCX*I_ERI_G2x2y_S_S_S_M1_vrr+WQX*I_ERI_G2x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_Px_S_M1_vrr = QCX*I_ERI_G2xyz_S_S_S_M1_vrr+WQX*I_ERI_G2xyz_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_Px_S_M1_vrr = QCX*I_ERI_G2x2z_S_S_S_M1_vrr+WQX*I_ERI_G2x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_Px_S_M1_vrr = QCX*I_ERI_Gx3y_S_S_S_M1_vrr+WQX*I_ERI_Gx3y_S_S_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Px_S_M1_vrr = QCX*I_ERI_Gx2yz_S_S_S_M1_vrr+WQX*I_ERI_Gx2yz_S_S_S_M2_vrr+oned2k*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Px_S_M1_vrr = QCX*I_ERI_Gxy2z_S_S_S_M1_vrr+WQX*I_ERI_Gxy2z_S_S_S_M2_vrr+oned2k*I_ERI_Fy2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_Px_S_M1_vrr = QCX*I_ERI_Gx3z_S_S_S_M1_vrr+WQX*I_ERI_Gx3z_S_S_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_Px_S_M1_vrr = QCX*I_ERI_G4y_S_S_S_M1_vrr+WQX*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_Px_S_M1_vrr = QCX*I_ERI_G3yz_S_S_S_M1_vrr+WQX*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_Px_S_M1_vrr = QCX*I_ERI_G2y2z_S_S_S_M1_vrr+WQX*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_Px_S_M1_vrr = QCX*I_ERI_Gy3z_S_S_S_M1_vrr+WQX*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_Px_S_M1_vrr = QCX*I_ERI_G4z_S_S_S_M1_vrr+WQX*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_G4x_S_Py_S_M1_vrr = QCY*I_ERI_G4x_S_S_S_M1_vrr+WQY*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_Py_S_M1_vrr = QCY*I_ERI_G3xy_S_S_S_M1_vrr+WQY*I_ERI_G3xy_S_S_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_Py_S_M1_vrr = QCY*I_ERI_G3xz_S_S_S_M1_vrr+WQY*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_Py_S_M1_vrr = QCY*I_ERI_G2x2y_S_S_S_M1_vrr+WQY*I_ERI_G2x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_Py_S_M1_vrr = QCY*I_ERI_G2xyz_S_S_S_M1_vrr+WQY*I_ERI_G2xyz_S_S_S_M2_vrr+oned2k*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_Py_S_M1_vrr = QCY*I_ERI_G2x2z_S_S_S_M1_vrr+WQY*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_Py_S_M1_vrr = QCY*I_ERI_Gx3y_S_S_S_M1_vrr+WQY*I_ERI_Gx3y_S_S_S_M2_vrr+3*oned2k*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Py_S_M1_vrr = QCY*I_ERI_Gx2yz_S_S_S_M1_vrr+WQY*I_ERI_Gx2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Py_S_M1_vrr = QCY*I_ERI_Gxy2z_S_S_S_M1_vrr+WQY*I_ERI_Gxy2z_S_S_S_M2_vrr+oned2k*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_Py_S_M1_vrr = QCY*I_ERI_Gx3z_S_S_S_M1_vrr+WQY*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_Py_S_M1_vrr = QCY*I_ERI_G4y_S_S_S_M1_vrr+WQY*I_ERI_G4y_S_S_S_M2_vrr+4*oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_Py_S_M1_vrr = QCY*I_ERI_G3yz_S_S_S_M1_vrr+WQY*I_ERI_G3yz_S_S_S_M2_vrr+3*oned2k*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_Py_S_M1_vrr = QCY*I_ERI_G2y2z_S_S_S_M1_vrr+WQY*I_ERI_G2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fy2z_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_Py_S_M1_vrr = QCY*I_ERI_Gy3z_S_S_S_M1_vrr+WQY*I_ERI_Gy3z_S_S_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_Py_S_M1_vrr = QCY*I_ERI_G4z_S_S_S_M1_vrr+WQY*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_G4x_S_Pz_S_M1_vrr = QCZ*I_ERI_G4x_S_S_S_M1_vrr+WQZ*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_Pz_S_M1_vrr = QCZ*I_ERI_G3xy_S_S_S_M1_vrr+WQZ*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_Pz_S_M1_vrr = QCZ*I_ERI_G3xz_S_S_S_M1_vrr+WQZ*I_ERI_G3xz_S_S_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_Pz_S_M1_vrr = QCZ*I_ERI_G2x2y_S_S_S_M1_vrr+WQZ*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_Pz_S_M1_vrr = QCZ*I_ERI_G2xyz_S_S_S_M1_vrr+WQZ*I_ERI_G2xyz_S_S_S_M2_vrr+oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_Pz_S_M1_vrr = QCZ*I_ERI_G2x2z_S_S_S_M1_vrr+WQZ*I_ERI_G2x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_Pz_S_M1_vrr = QCZ*I_ERI_Gx3y_S_S_S_M1_vrr+WQZ*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Pz_S_M1_vrr = QCZ*I_ERI_Gx2yz_S_S_S_M1_vrr+WQZ*I_ERI_Gx2yz_S_S_S_M2_vrr+oned2k*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Pz_S_M1_vrr = QCZ*I_ERI_Gxy2z_S_S_S_M1_vrr+WQZ*I_ERI_Gxy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_Pz_S_M1_vrr = QCZ*I_ERI_Gx3z_S_S_S_M1_vrr+WQZ*I_ERI_Gx3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_Pz_S_M1_vrr = QCZ*I_ERI_G4y_S_S_S_M1_vrr+WQZ*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_Pz_S_M1_vrr = QCZ*I_ERI_G3yz_S_S_S_M1_vrr+WQZ*I_ERI_G3yz_S_S_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_Pz_S_M1_vrr = QCZ*I_ERI_G2y2z_S_S_S_M1_vrr+WQZ*I_ERI_G2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_Pz_S_M1_vrr = QCZ*I_ERI_Gy3z_S_S_S_M1_vrr+WQZ*I_ERI_Gy3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Fy2z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_Pz_S_M1_vrr = QCZ*I_ERI_G4z_S_S_S_M1_vrr+WQZ*I_ERI_G4z_S_S_S_M2_vrr+4*oned2k*I_ERI_F3z_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_H_S_P_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       ************************************************************/
      Double I_ERI_H5x_S_Px_S_M1_vrr = QCX*I_ERI_H5x_S_S_S_M1_vrr+WQX*I_ERI_H5x_S_S_S_M2_vrr+5*oned2k*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H4xy_S_Px_S_M1_vrr = QCX*I_ERI_H4xy_S_S_S_M1_vrr+WQX*I_ERI_H4xy_S_S_S_M2_vrr+4*oned2k*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_H4xz_S_Px_S_M1_vrr = QCX*I_ERI_H4xz_S_S_S_M1_vrr+WQX*I_ERI_H4xz_S_S_S_M2_vrr+4*oned2k*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_H3x2y_S_Px_S_M1_vrr = QCX*I_ERI_H3x2y_S_S_S_M1_vrr+WQX*I_ERI_H3x2y_S_S_S_M2_vrr+3*oned2k*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_H3xyz_S_Px_S_M1_vrr = QCX*I_ERI_H3xyz_S_S_S_M1_vrr+WQX*I_ERI_H3xyz_S_S_S_M2_vrr+3*oned2k*I_ERI_G2xyz_S_S_S_M2_vrr;
      Double I_ERI_H3x2z_S_Px_S_M1_vrr = QCX*I_ERI_H3x2z_S_S_S_M1_vrr+WQX*I_ERI_H3x2z_S_S_S_M2_vrr+3*oned2k*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_H2x3y_S_Px_S_M1_vrr = QCX*I_ERI_H2x3y_S_S_S_M1_vrr+WQX*I_ERI_H2x3y_S_S_S_M2_vrr+2*oned2k*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Px_S_M1_vrr = QCX*I_ERI_H2x2yz_S_S_S_M1_vrr+WQX*I_ERI_H2x2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_Gx2yz_S_S_S_M2_vrr;
      Double I_ERI_H2xy2z_S_Px_S_M1_vrr = QCX*I_ERI_H2xy2z_S_S_S_M1_vrr+WQX*I_ERI_H2xy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gxy2z_S_S_S_M2_vrr;
      Double I_ERI_H2x3z_S_Px_S_M1_vrr = QCX*I_ERI_H2x3z_S_S_S_M1_vrr+WQX*I_ERI_H2x3z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4y_S_Px_S_M1_vrr = QCX*I_ERI_Hx4y_S_S_S_M1_vrr+WQX*I_ERI_Hx4y_S_S_S_M2_vrr+oned2k*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_Hx3yz_S_Px_S_M1_vrr = QCX*I_ERI_Hx3yz_S_S_S_M1_vrr+WQX*I_ERI_Hx3yz_S_S_S_M2_vrr+oned2k*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_Hx2y2z_S_Px_S_M1_vrr = QCX*I_ERI_Hx2y2z_S_S_S_M1_vrr+WQX*I_ERI_Hx2y2z_S_S_S_M2_vrr+oned2k*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_Hxy3z_S_Px_S_M1_vrr = QCX*I_ERI_Hxy3z_S_S_S_M1_vrr+WQX*I_ERI_Hxy3z_S_S_S_M2_vrr+oned2k*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4z_S_Px_S_M1_vrr = QCX*I_ERI_Hx4z_S_S_S_M1_vrr+WQX*I_ERI_Hx4z_S_S_S_M2_vrr+oned2k*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_H5y_S_Px_S_M1_vrr = QCX*I_ERI_H5y_S_S_S_M1_vrr+WQX*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_H4yz_S_Px_S_M1_vrr = QCX*I_ERI_H4yz_S_S_S_M1_vrr+WQX*I_ERI_H4yz_S_S_S_M2_vrr;
      Double I_ERI_H3y2z_S_Px_S_M1_vrr = QCX*I_ERI_H3y2z_S_S_S_M1_vrr+WQX*I_ERI_H3y2z_S_S_S_M2_vrr;
      Double I_ERI_H2y3z_S_Px_S_M1_vrr = QCX*I_ERI_H2y3z_S_S_S_M1_vrr+WQX*I_ERI_H2y3z_S_S_S_M2_vrr;
      Double I_ERI_Hy4z_S_Px_S_M1_vrr = QCX*I_ERI_Hy4z_S_S_S_M1_vrr+WQX*I_ERI_Hy4z_S_S_S_M2_vrr;
      Double I_ERI_H5z_S_Px_S_M1_vrr = QCX*I_ERI_H5z_S_S_S_M1_vrr+WQX*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_H5x_S_Py_S_M1_vrr = QCY*I_ERI_H5x_S_S_S_M1_vrr+WQY*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_H4xy_S_Py_S_M1_vrr = QCY*I_ERI_H4xy_S_S_S_M1_vrr+WQY*I_ERI_H4xy_S_S_S_M2_vrr+oned2k*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H4xz_S_Py_S_M1_vrr = QCY*I_ERI_H4xz_S_S_S_M1_vrr+WQY*I_ERI_H4xz_S_S_S_M2_vrr;
      Double I_ERI_H3x2y_S_Py_S_M1_vrr = QCY*I_ERI_H3x2y_S_S_S_M1_vrr+WQY*I_ERI_H3x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_H3xyz_S_Py_S_M1_vrr = QCY*I_ERI_H3xyz_S_S_S_M1_vrr+WQY*I_ERI_H3xyz_S_S_S_M2_vrr+oned2k*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_H3x2z_S_Py_S_M1_vrr = QCY*I_ERI_H3x2z_S_S_S_M1_vrr+WQY*I_ERI_H3x2z_S_S_S_M2_vrr;
      Double I_ERI_H2x3y_S_Py_S_M1_vrr = QCY*I_ERI_H2x3y_S_S_S_M1_vrr+WQY*I_ERI_H2x3y_S_S_S_M2_vrr+3*oned2k*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Py_S_M1_vrr = QCY*I_ERI_H2x2yz_S_S_S_M1_vrr+WQY*I_ERI_H2x2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_G2xyz_S_S_S_M2_vrr;
      Double I_ERI_H2xy2z_S_Py_S_M1_vrr = QCY*I_ERI_H2xy2z_S_S_S_M1_vrr+WQY*I_ERI_H2xy2z_S_S_S_M2_vrr+oned2k*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_H2x3z_S_Py_S_M1_vrr = QCY*I_ERI_H2x3z_S_S_S_M1_vrr+WQY*I_ERI_H2x3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4y_S_Py_S_M1_vrr = QCY*I_ERI_Hx4y_S_S_S_M1_vrr+WQY*I_ERI_Hx4y_S_S_S_M2_vrr+4*oned2k*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_Hx3yz_S_Py_S_M1_vrr = QCY*I_ERI_Hx3yz_S_S_S_M1_vrr+WQY*I_ERI_Hx3yz_S_S_S_M2_vrr+3*oned2k*I_ERI_Gx2yz_S_S_S_M2_vrr;
      Double I_ERI_Hx2y2z_S_Py_S_M1_vrr = QCY*I_ERI_Hx2y2z_S_S_S_M1_vrr+WQY*I_ERI_Hx2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gxy2z_S_S_S_M2_vrr;
      Double I_ERI_Hxy3z_S_Py_S_M1_vrr = QCY*I_ERI_Hxy3z_S_S_S_M1_vrr+WQY*I_ERI_Hxy3z_S_S_S_M2_vrr+oned2k*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4z_S_Py_S_M1_vrr = QCY*I_ERI_Hx4z_S_S_S_M1_vrr+WQY*I_ERI_Hx4z_S_S_S_M2_vrr;
      Double I_ERI_H5y_S_Py_S_M1_vrr = QCY*I_ERI_H5y_S_S_S_M1_vrr+WQY*I_ERI_H5y_S_S_S_M2_vrr+5*oned2k*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_H4yz_S_Py_S_M1_vrr = QCY*I_ERI_H4yz_S_S_S_M1_vrr+WQY*I_ERI_H4yz_S_S_S_M2_vrr+4*oned2k*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_H3y2z_S_Py_S_M1_vrr = QCY*I_ERI_H3y2z_S_S_S_M1_vrr+WQY*I_ERI_H3y2z_S_S_S_M2_vrr+3*oned2k*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_H2y3z_S_Py_S_M1_vrr = QCY*I_ERI_H2y3z_S_S_S_M1_vrr+WQY*I_ERI_H2y3z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_Hy4z_S_Py_S_M1_vrr = QCY*I_ERI_Hy4z_S_S_S_M1_vrr+WQY*I_ERI_Hy4z_S_S_S_M2_vrr+oned2k*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_H5z_S_Py_S_M1_vrr = QCY*I_ERI_H5z_S_S_S_M1_vrr+WQY*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_H5x_S_Pz_S_M1_vrr = QCZ*I_ERI_H5x_S_S_S_M1_vrr+WQZ*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_H4xy_S_Pz_S_M1_vrr = QCZ*I_ERI_H4xy_S_S_S_M1_vrr+WQZ*I_ERI_H4xy_S_S_S_M2_vrr;
      Double I_ERI_H4xz_S_Pz_S_M1_vrr = QCZ*I_ERI_H4xz_S_S_S_M1_vrr+WQZ*I_ERI_H4xz_S_S_S_M2_vrr+oned2k*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H3x2y_S_Pz_S_M1_vrr = QCZ*I_ERI_H3x2y_S_S_S_M1_vrr+WQZ*I_ERI_H3x2y_S_S_S_M2_vrr;
      Double I_ERI_H3xyz_S_Pz_S_M1_vrr = QCZ*I_ERI_H3xyz_S_S_S_M1_vrr+WQZ*I_ERI_H3xyz_S_S_S_M2_vrr+oned2k*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_H3x2z_S_Pz_S_M1_vrr = QCZ*I_ERI_H3x2z_S_S_S_M1_vrr+WQZ*I_ERI_H3x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_H2x3y_S_Pz_S_M1_vrr = QCZ*I_ERI_H2x3y_S_S_S_M1_vrr+WQZ*I_ERI_H2x3y_S_S_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Pz_S_M1_vrr = QCZ*I_ERI_H2x2yz_S_S_S_M1_vrr+WQZ*I_ERI_H2x2yz_S_S_S_M2_vrr+oned2k*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_H2xy2z_S_Pz_S_M1_vrr = QCZ*I_ERI_H2xy2z_S_S_S_M1_vrr+WQZ*I_ERI_H2xy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_G2xyz_S_S_S_M2_vrr;
      Double I_ERI_H2x3z_S_Pz_S_M1_vrr = QCZ*I_ERI_H2x3z_S_S_S_M1_vrr+WQZ*I_ERI_H2x3z_S_S_S_M2_vrr+3*oned2k*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_Hx4y_S_Pz_S_M1_vrr = QCZ*I_ERI_Hx4y_S_S_S_M1_vrr+WQZ*I_ERI_Hx4y_S_S_S_M2_vrr;
      Double I_ERI_Hx3yz_S_Pz_S_M1_vrr = QCZ*I_ERI_Hx3yz_S_S_S_M1_vrr+WQZ*I_ERI_Hx3yz_S_S_S_M2_vrr+oned2k*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_Hx2y2z_S_Pz_S_M1_vrr = QCZ*I_ERI_Hx2y2z_S_S_S_M1_vrr+WQZ*I_ERI_Hx2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gx2yz_S_S_S_M2_vrr;
      Double I_ERI_Hxy3z_S_Pz_S_M1_vrr = QCZ*I_ERI_Hxy3z_S_S_S_M1_vrr+WQZ*I_ERI_Hxy3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Gxy2z_S_S_S_M2_vrr;
      Double I_ERI_Hx4z_S_Pz_S_M1_vrr = QCZ*I_ERI_Hx4z_S_S_S_M1_vrr+WQZ*I_ERI_Hx4z_S_S_S_M2_vrr+4*oned2k*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_H5y_S_Pz_S_M1_vrr = QCZ*I_ERI_H5y_S_S_S_M1_vrr+WQZ*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_H4yz_S_Pz_S_M1_vrr = QCZ*I_ERI_H4yz_S_S_S_M1_vrr+WQZ*I_ERI_H4yz_S_S_S_M2_vrr+oned2k*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_H3y2z_S_Pz_S_M1_vrr = QCZ*I_ERI_H3y2z_S_S_S_M1_vrr+WQZ*I_ERI_H3y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_H2y3z_S_Pz_S_M1_vrr = QCZ*I_ERI_H2y3z_S_S_S_M1_vrr+WQZ*I_ERI_H2y3z_S_S_S_M2_vrr+3*oned2k*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_Hy4z_S_Pz_S_M1_vrr = QCZ*I_ERI_Hy4z_S_S_S_M1_vrr+WQZ*I_ERI_Hy4z_S_S_S_M2_vrr+4*oned2k*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_H5z_S_Pz_S_M1_vrr = QCZ*I_ERI_H5z_S_S_S_M1_vrr+WQZ*I_ERI_H5z_S_S_S_M2_vrr+5*oned2k*I_ERI_G4z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       ************************************************************/
      Double I_ERI_I6x_S_S_S_M1_vrr = PAX*I_ERI_H5x_S_S_S_M1_vrr+WPX*I_ERI_H5x_S_S_S_M2_vrr+5*oned2z*I_ERI_G4x_S_S_S_M1_vrr-5*rhod2zsq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_I5xy_S_S_S_M1_vrr = PAY*I_ERI_H5x_S_S_S_M1_vrr+WPY*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_I5xz_S_S_S_M1_vrr = PAZ*I_ERI_H5x_S_S_S_M1_vrr+WPZ*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_I4x2y_S_S_S_M1_vrr = PAY*I_ERI_H4xy_S_S_S_M1_vrr+WPY*I_ERI_H4xy_S_S_S_M2_vrr+oned2z*I_ERI_G4x_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_I4x2z_S_S_S_M1_vrr = PAZ*I_ERI_H4xz_S_S_S_M1_vrr+WPZ*I_ERI_H4xz_S_S_S_M2_vrr+oned2z*I_ERI_G4x_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_I3x3y_S_S_S_M1_vrr = PAY*I_ERI_H3x2y_S_S_S_M1_vrr+WPY*I_ERI_H3x2y_S_S_S_M2_vrr+2*oned2z*I_ERI_G3xy_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_I3x2yz_S_S_S_M1_vrr = PAZ*I_ERI_H3x2y_S_S_S_M1_vrr+WPZ*I_ERI_H3x2y_S_S_S_M2_vrr;
      Double I_ERI_I3x3z_S_S_S_M1_vrr = PAZ*I_ERI_H3x2z_S_S_S_M1_vrr+WPZ*I_ERI_H3x2z_S_S_S_M2_vrr+2*oned2z*I_ERI_G3xz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_I2x4y_S_S_S_M1_vrr = PAX*I_ERI_Hx4y_S_S_S_M1_vrr+WPX*I_ERI_Hx4y_S_S_S_M2_vrr+oned2z*I_ERI_G4y_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_I2x3yz_S_S_S_M1_vrr = PAZ*I_ERI_H2x3y_S_S_S_M1_vrr+WPZ*I_ERI_H2x3y_S_S_S_M2_vrr;
      Double I_ERI_I2xy3z_S_S_S_M1_vrr = PAY*I_ERI_H2x3z_S_S_S_M1_vrr+WPY*I_ERI_H2x3z_S_S_S_M2_vrr;
      Double I_ERI_I2x4z_S_S_S_M1_vrr = PAX*I_ERI_Hx4z_S_S_S_M1_vrr+WPX*I_ERI_Hx4z_S_S_S_M2_vrr+oned2z*I_ERI_G4z_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_Ix5y_S_S_S_M1_vrr = PAX*I_ERI_H5y_S_S_S_M1_vrr+WPX*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_Ix5z_S_S_S_M1_vrr = PAX*I_ERI_H5z_S_S_S_M1_vrr+WPX*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_I6y_S_S_S_M1_vrr = PAY*I_ERI_H5y_S_S_S_M1_vrr+WPY*I_ERI_H5y_S_S_S_M2_vrr+5*oned2z*I_ERI_G4y_S_S_S_M1_vrr-5*rhod2zsq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_I5yz_S_S_S_M1_vrr = PAZ*I_ERI_H5y_S_S_S_M1_vrr+WPZ*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_I4y2z_S_S_S_M1_vrr = PAZ*I_ERI_H4yz_S_S_S_M1_vrr+WPZ*I_ERI_H4yz_S_S_S_M2_vrr+oned2z*I_ERI_G4y_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_I3y3z_S_S_S_M1_vrr = PAZ*I_ERI_H3y2z_S_S_S_M1_vrr+WPZ*I_ERI_H3y2z_S_S_S_M2_vrr+2*oned2z*I_ERI_G3yz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_I2y4z_S_S_S_M1_vrr = PAY*I_ERI_Hy4z_S_S_S_M1_vrr+WPY*I_ERI_Hy4z_S_S_S_M2_vrr+oned2z*I_ERI_G4z_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_Iy5z_S_S_S_M1_vrr = PAY*I_ERI_H5z_S_S_S_M1_vrr+WPY*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_I6z_S_S_S_M1_vrr = PAZ*I_ERI_H5z_S_S_S_M1_vrr+WPZ*I_ERI_H5z_S_S_S_M2_vrr+5*oned2z*I_ERI_G4z_S_S_S_M1_vrr-5*rhod2zsq*I_ERI_G4z_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_H_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_P_S
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_Px_S_vrr = PAX*I_ERI_G4x_S_Px_S_vrr+WPX*I_ERI_G4x_S_Px_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Px_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Px_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_Px_S_vrr = PAY*I_ERI_G4x_S_Px_S_vrr+WPY*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_Px_S_vrr = PAZ*I_ERI_G4x_S_Px_S_vrr+WPZ*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_Px_S_vrr = PAY*I_ERI_G3xy_S_Px_S_vrr+WPY*I_ERI_G3xy_S_Px_S_M1_vrr+oned2z*I_ERI_F3x_S_Px_S_vrr-rhod2zsq*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_Px_S_vrr = PAZ*I_ERI_G3xy_S_Px_S_vrr+WPZ*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_Px_S_vrr = PAZ*I_ERI_G3xz_S_Px_S_vrr+WPZ*I_ERI_G3xz_S_Px_S_M1_vrr+oned2z*I_ERI_F3x_S_Px_S_vrr-rhod2zsq*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_Px_S_vrr = PAX*I_ERI_Gx3y_S_Px_S_vrr+WPX*I_ERI_Gx3y_S_Px_S_M1_vrr+oned2z*I_ERI_F3y_S_Px_S_vrr-rhod2zsq*I_ERI_F3y_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Px_S_vrr = PAZ*I_ERI_G2x2y_S_Px_S_vrr+WPZ*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Px_S_vrr = PAY*I_ERI_G2x2z_S_Px_S_vrr+WPY*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_Px_S_vrr = PAX*I_ERI_Gx3z_S_Px_S_vrr+WPX*I_ERI_Gx3z_S_Px_S_M1_vrr+oned2z*I_ERI_F3z_S_Px_S_vrr-rhod2zsq*I_ERI_F3z_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_Px_S_vrr = PAX*I_ERI_G4y_S_Px_S_vrr+WPX*I_ERI_G4y_S_Px_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Px_S_vrr = PAZ*I_ERI_Gx3y_S_Px_S_vrr+WPZ*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Px_S_vrr = PAX*I_ERI_G2y2z_S_Px_S_vrr+WPX*I_ERI_G2y2z_S_Px_S_M1_vrr+oned2k*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Px_S_vrr = PAY*I_ERI_Gx3z_S_Px_S_vrr+WPY*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_Px_S_vrr = PAX*I_ERI_G4z_S_Px_S_vrr+WPX*I_ERI_G4z_S_Px_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_Px_S_vrr = PAY*I_ERI_G4y_S_Px_S_vrr+WPY*I_ERI_G4y_S_Px_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Px_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_H4yz_S_Px_S_vrr = PAZ*I_ERI_G4y_S_Px_S_vrr+WPZ*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_H3y2z_S_Px_S_vrr = PAZ*I_ERI_G3yz_S_Px_S_vrr+WPZ*I_ERI_G3yz_S_Px_S_M1_vrr+oned2z*I_ERI_F3y_S_Px_S_vrr-rhod2zsq*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_H2y3z_S_Px_S_vrr = PAY*I_ERI_Gy3z_S_Px_S_vrr+WPY*I_ERI_Gy3z_S_Px_S_M1_vrr+oned2z*I_ERI_F3z_S_Px_S_vrr-rhod2zsq*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_Hy4z_S_Px_S_vrr = PAY*I_ERI_G4z_S_Px_S_vrr+WPY*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5z_S_Px_S_vrr = PAZ*I_ERI_G4z_S_Px_S_vrr+WPZ*I_ERI_G4z_S_Px_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Px_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_H5x_S_Py_S_vrr = PAX*I_ERI_G4x_S_Py_S_vrr+WPX*I_ERI_G4x_S_Py_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Py_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_H4xy_S_Py_S_vrr = PAY*I_ERI_G4x_S_Py_S_vrr+WPY*I_ERI_G4x_S_Py_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_Py_S_vrr = PAZ*I_ERI_G4x_S_Py_S_vrr+WPZ*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H3x2y_S_Py_S_vrr = PAY*I_ERI_G3xy_S_Py_S_vrr+WPY*I_ERI_G3xy_S_Py_S_M1_vrr+oned2z*I_ERI_F3x_S_Py_S_vrr-rhod2zsq*I_ERI_F3x_S_Py_S_M1_vrr+oned2k*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_Py_S_vrr = PAZ*I_ERI_G3xy_S_Py_S_vrr+WPZ*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_H3x2z_S_Py_S_vrr = PAZ*I_ERI_G3xz_S_Py_S_vrr+WPZ*I_ERI_G3xz_S_Py_S_M1_vrr+oned2z*I_ERI_F3x_S_Py_S_vrr-rhod2zsq*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_H2x3y_S_Py_S_vrr = PAX*I_ERI_Gx3y_S_Py_S_vrr+WPX*I_ERI_Gx3y_S_Py_S_M1_vrr+oned2z*I_ERI_F3y_S_Py_S_vrr-rhod2zsq*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Py_S_vrr = PAZ*I_ERI_G2x2y_S_Py_S_vrr+WPZ*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Py_S_vrr = PAY*I_ERI_G2x2z_S_Py_S_vrr+WPY*I_ERI_G2x2z_S_Py_S_M1_vrr+oned2k*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3z_S_Py_S_vrr = PAX*I_ERI_Gx3z_S_Py_S_vrr+WPX*I_ERI_Gx3z_S_Py_S_M1_vrr+oned2z*I_ERI_F3z_S_Py_S_vrr-rhod2zsq*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4y_S_Py_S_vrr = PAX*I_ERI_G4y_S_Py_S_vrr+WPX*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Py_S_vrr = PAZ*I_ERI_Gx3y_S_Py_S_vrr+WPZ*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Py_S_vrr = PAX*I_ERI_G2y2z_S_Py_S_vrr+WPX*I_ERI_G2y2z_S_Py_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Py_S_vrr = PAY*I_ERI_Gx3z_S_Py_S_vrr+WPY*I_ERI_Gx3z_S_Py_S_M1_vrr+oned2k*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4z_S_Py_S_vrr = PAX*I_ERI_G4z_S_Py_S_vrr+WPX*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5y_S_Py_S_vrr = PAY*I_ERI_G4y_S_Py_S_vrr+WPY*I_ERI_G4y_S_Py_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Py_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Py_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_Py_S_vrr = PAZ*I_ERI_G4y_S_Py_S_vrr+WPZ*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_H3y2z_S_Py_S_vrr = PAZ*I_ERI_G3yz_S_Py_S_vrr+WPZ*I_ERI_G3yz_S_Py_S_M1_vrr+oned2z*I_ERI_F3y_S_Py_S_vrr-rhod2zsq*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_H2y3z_S_Py_S_vrr = PAY*I_ERI_Gy3z_S_Py_S_vrr+WPY*I_ERI_Gy3z_S_Py_S_M1_vrr+oned2z*I_ERI_F3z_S_Py_S_vrr-rhod2zsq*I_ERI_F3z_S_Py_S_M1_vrr+oned2k*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_Py_S_vrr = PAY*I_ERI_G4z_S_Py_S_vrr+WPY*I_ERI_G4z_S_Py_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_Py_S_vrr = PAZ*I_ERI_G4z_S_Py_S_vrr+WPZ*I_ERI_G4z_S_Py_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Py_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_H5x_S_Pz_S_vrr = PAX*I_ERI_G4x_S_Pz_S_vrr+WPX*I_ERI_G4x_S_Pz_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Pz_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_H4xy_S_Pz_S_vrr = PAY*I_ERI_G4x_S_Pz_S_vrr+WPY*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_H4xz_S_Pz_S_vrr = PAZ*I_ERI_G4x_S_Pz_S_vrr+WPZ*I_ERI_G4x_S_Pz_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_Pz_S_vrr = PAY*I_ERI_G3xy_S_Pz_S_vrr+WPY*I_ERI_G3xy_S_Pz_S_M1_vrr+oned2z*I_ERI_F3x_S_Pz_S_vrr-rhod2zsq*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_H3xyz_S_Pz_S_vrr = PAZ*I_ERI_G3xy_S_Pz_S_vrr+WPZ*I_ERI_G3xy_S_Pz_S_M1_vrr+oned2k*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3x2z_S_Pz_S_vrr = PAZ*I_ERI_G3xz_S_Pz_S_vrr+WPZ*I_ERI_G3xz_S_Pz_S_M1_vrr+oned2z*I_ERI_F3x_S_Pz_S_vrr-rhod2zsq*I_ERI_F3x_S_Pz_S_M1_vrr+oned2k*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_Pz_S_vrr = PAX*I_ERI_Gx3y_S_Pz_S_vrr+WPX*I_ERI_Gx3y_S_Pz_S_M1_vrr+oned2z*I_ERI_F3y_S_Pz_S_vrr-rhod2zsq*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Pz_S_vrr = PAZ*I_ERI_G2x2y_S_Pz_S_vrr+WPZ*I_ERI_G2x2y_S_Pz_S_M1_vrr+oned2k*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Pz_S_vrr = PAY*I_ERI_G2x2z_S_Pz_S_vrr+WPY*I_ERI_G2x2z_S_Pz_S_M1_vrr;
      Double I_ERI_H2x3z_S_Pz_S_vrr = PAX*I_ERI_Gx3z_S_Pz_S_vrr+WPX*I_ERI_Gx3z_S_Pz_S_M1_vrr+oned2z*I_ERI_F3z_S_Pz_S_vrr-rhod2zsq*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4y_S_Pz_S_vrr = PAX*I_ERI_G4y_S_Pz_S_vrr+WPX*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Pz_S_vrr = PAZ*I_ERI_Gx3y_S_Pz_S_vrr+WPZ*I_ERI_Gx3y_S_Pz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Pz_S_vrr = PAX*I_ERI_G2y2z_S_Pz_S_vrr+WPX*I_ERI_G2y2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Pz_S_vrr = PAY*I_ERI_Gx3z_S_Pz_S_vrr+WPY*I_ERI_Gx3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4z_S_Pz_S_vrr = PAX*I_ERI_G4z_S_Pz_S_vrr+WPX*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_H5y_S_Pz_S_vrr = PAY*I_ERI_G4y_S_Pz_S_vrr+WPY*I_ERI_G4y_S_Pz_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Pz_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_H4yz_S_Pz_S_vrr = PAZ*I_ERI_G4y_S_Pz_S_vrr+WPZ*I_ERI_G4y_S_Pz_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_Pz_S_vrr = PAZ*I_ERI_G3yz_S_Pz_S_vrr+WPZ*I_ERI_G3yz_S_Pz_S_M1_vrr+oned2z*I_ERI_F3y_S_Pz_S_vrr-rhod2zsq*I_ERI_F3y_S_Pz_S_M1_vrr+oned2k*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_Pz_S_vrr = PAY*I_ERI_Gy3z_S_Pz_S_vrr+WPY*I_ERI_Gy3z_S_Pz_S_M1_vrr+oned2z*I_ERI_F3z_S_Pz_S_vrr-rhod2zsq*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hy4z_S_Pz_S_vrr = PAY*I_ERI_G4z_S_Pz_S_vrr+WPY*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_H5z_S_Pz_S_vrr = PAZ*I_ERI_G4z_S_Pz_S_vrr+WPZ*I_ERI_G4z_S_Pz_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Pz_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Pz_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       ************************************************************/
      Double I_ERI_I6x_S_S_S_vrr = PAX*I_ERI_H5x_S_S_S_vrr+WPX*I_ERI_H5x_S_S_S_M1_vrr+5*oned2z*I_ERI_G4x_S_S_S_vrr-5*rhod2zsq*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_I5xy_S_S_S_vrr = PAY*I_ERI_H5x_S_S_S_vrr+WPY*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I5xz_S_S_S_vrr = PAZ*I_ERI_H5x_S_S_S_vrr+WPZ*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I4x2y_S_S_S_vrr = PAY*I_ERI_H4xy_S_S_S_vrr+WPY*I_ERI_H4xy_S_S_S_M1_vrr+oned2z*I_ERI_G4x_S_S_S_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_I4xyz_S_S_S_vrr = PAZ*I_ERI_H4xy_S_S_S_vrr+WPZ*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_I4x2z_S_S_S_vrr = PAZ*I_ERI_H4xz_S_S_S_vrr+WPZ*I_ERI_H4xz_S_S_S_M1_vrr+oned2z*I_ERI_G4x_S_S_S_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_I3x3y_S_S_S_vrr = PAY*I_ERI_H3x2y_S_S_S_vrr+WPY*I_ERI_H3x2y_S_S_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_S_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_I3x2yz_S_S_S_vrr = PAZ*I_ERI_H3x2y_S_S_S_vrr+WPZ*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_I3xy2z_S_S_S_vrr = PAY*I_ERI_H3x2z_S_S_S_vrr+WPY*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_I3x3z_S_S_S_vrr = PAZ*I_ERI_H3x2z_S_S_S_vrr+WPZ*I_ERI_H3x2z_S_S_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_S_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_I2x4y_S_S_S_vrr = PAX*I_ERI_Hx4y_S_S_S_vrr+WPX*I_ERI_Hx4y_S_S_S_M1_vrr+oned2z*I_ERI_G4y_S_S_S_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_I2x3yz_S_S_S_vrr = PAZ*I_ERI_H2x3y_S_S_S_vrr+WPZ*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_S_S_vrr = PAZ*I_ERI_H2x2yz_S_S_S_vrr+WPZ*I_ERI_H2x2yz_S_S_S_M1_vrr+oned2z*I_ERI_G2x2y_S_S_S_vrr-rhod2zsq*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_I2xy3z_S_S_S_vrr = PAY*I_ERI_H2x3z_S_S_S_vrr+WPY*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_I2x4z_S_S_S_vrr = PAX*I_ERI_Hx4z_S_S_S_vrr+WPX*I_ERI_Hx4z_S_S_S_M1_vrr+oned2z*I_ERI_G4z_S_S_S_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5y_S_S_S_vrr = PAX*I_ERI_H5y_S_S_S_vrr+WPX*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_Ix4yz_S_S_S_vrr = PAZ*I_ERI_Hx4y_S_S_S_vrr+WPZ*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_S_S_vrr = PAX*I_ERI_H3y2z_S_S_S_vrr+WPX*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_S_S_vrr = PAX*I_ERI_H2y3z_S_S_S_vrr+WPX*I_ERI_H2y3z_S_S_S_M1_vrr;
      Double I_ERI_Ixy4z_S_S_S_vrr = PAY*I_ERI_Hx4z_S_S_S_vrr+WPY*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5z_S_S_S_vrr = PAX*I_ERI_H5z_S_S_S_vrr+WPX*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_I6y_S_S_S_vrr = PAY*I_ERI_H5y_S_S_S_vrr+WPY*I_ERI_H5y_S_S_S_M1_vrr+5*oned2z*I_ERI_G4y_S_S_S_vrr-5*rhod2zsq*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_I5yz_S_S_S_vrr = PAZ*I_ERI_H5y_S_S_S_vrr+WPZ*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_I4y2z_S_S_S_vrr = PAZ*I_ERI_H4yz_S_S_S_vrr+WPZ*I_ERI_H4yz_S_S_S_M1_vrr+oned2z*I_ERI_G4y_S_S_S_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_I3y3z_S_S_S_vrr = PAZ*I_ERI_H3y2z_S_S_S_vrr+WPZ*I_ERI_H3y2z_S_S_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_S_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_I2y4z_S_S_S_vrr = PAY*I_ERI_Hy4z_S_S_S_vrr+WPY*I_ERI_Hy4z_S_S_S_M1_vrr+oned2z*I_ERI_G4z_S_S_S_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_Iy5z_S_S_S_vrr = PAY*I_ERI_H5z_S_S_S_vrr+WPY*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_I6z_S_S_S_vrr = PAZ*I_ERI_H5z_S_S_S_vrr+WPZ*I_ERI_H5z_S_S_S_M1_vrr+5*oned2z*I_ERI_G4z_S_S_S_vrr-5*rhod2zsq*I_ERI_G4z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_P_S
       * RHS shell quartet name: SQ_ERI_H_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_D2x_S_vrr = QCX*I_ERI_H5x_S_Px_S_vrr+WQX*I_ERI_H5x_S_Px_S_M1_vrr+oned2e*I_ERI_H5x_S_S_S_vrr-rhod2esq*I_ERI_H5x_S_S_S_M1_vrr+5*oned2k*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H4xy_S_D2x_S_vrr = QCX*I_ERI_H4xy_S_Px_S_vrr+WQX*I_ERI_H4xy_S_Px_S_M1_vrr+oned2e*I_ERI_H4xy_S_S_S_vrr-rhod2esq*I_ERI_H4xy_S_S_S_M1_vrr+4*oned2k*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_D2x_S_vrr = QCX*I_ERI_H4xz_S_Px_S_vrr+WQX*I_ERI_H4xz_S_Px_S_M1_vrr+oned2e*I_ERI_H4xz_S_S_S_vrr-rhod2esq*I_ERI_H4xz_S_S_S_M1_vrr+4*oned2k*I_ERI_G3xz_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_D2x_S_vrr = QCX*I_ERI_H3x2y_S_Px_S_vrr+WQX*I_ERI_H3x2y_S_Px_S_M1_vrr+oned2e*I_ERI_H3x2y_S_S_S_vrr-rhod2esq*I_ERI_H3x2y_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_D2x_S_vrr = QCX*I_ERI_H3xyz_S_Px_S_vrr+WQX*I_ERI_H3xyz_S_Px_S_M1_vrr+oned2e*I_ERI_H3xyz_S_S_S_vrr-rhod2esq*I_ERI_H3xyz_S_S_S_M1_vrr+3*oned2k*I_ERI_G2xyz_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_D2x_S_vrr = QCX*I_ERI_H3x2z_S_Px_S_vrr+WQX*I_ERI_H3x2z_S_Px_S_M1_vrr+oned2e*I_ERI_H3x2z_S_S_S_vrr-rhod2esq*I_ERI_H3x2z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_D2x_S_vrr = QCX*I_ERI_H2x3y_S_Px_S_vrr+WQX*I_ERI_H2x3y_S_Px_S_M1_vrr+oned2e*I_ERI_H2x3y_S_S_S_vrr-rhod2esq*I_ERI_H2x3y_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_H2x2yz_S_D2x_S_vrr = QCX*I_ERI_H2x2yz_S_Px_S_vrr+WQX*I_ERI_H2x2yz_S_Px_S_M1_vrr+oned2e*I_ERI_H2x2yz_S_S_S_vrr-rhod2esq*I_ERI_H2x2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_D2x_S_vrr = QCX*I_ERI_H2xy2z_S_Px_S_vrr+WQX*I_ERI_H2xy2z_S_Px_S_M1_vrr+oned2e*I_ERI_H2xy2z_S_S_S_vrr-rhod2esq*I_ERI_H2xy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gxy2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_D2x_S_vrr = QCX*I_ERI_H2x3z_S_Px_S_vrr+WQX*I_ERI_H2x3z_S_Px_S_M1_vrr+oned2e*I_ERI_H2x3z_S_S_S_vrr-rhod2esq*I_ERI_H2x3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4y_S_D2x_S_vrr = QCX*I_ERI_Hx4y_S_Px_S_vrr+WQX*I_ERI_Hx4y_S_Px_S_M1_vrr+oned2e*I_ERI_Hx4y_S_S_S_vrr-rhod2esq*I_ERI_Hx4y_S_S_S_M1_vrr+oned2k*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_Hx3yz_S_D2x_S_vrr = QCX*I_ERI_Hx3yz_S_Px_S_vrr+WQX*I_ERI_Hx3yz_S_Px_S_M1_vrr+oned2e*I_ERI_Hx3yz_S_S_S_vrr-rhod2esq*I_ERI_Hx3yz_S_S_S_M1_vrr+oned2k*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_D2x_S_vrr = QCX*I_ERI_Hx2y2z_S_Px_S_vrr+WQX*I_ERI_Hx2y2z_S_Px_S_M1_vrr+oned2e*I_ERI_Hx2y2z_S_S_S_vrr-rhod2esq*I_ERI_Hx2y2z_S_S_S_M1_vrr+oned2k*I_ERI_G2y2z_S_Px_S_M1_vrr;
      Double I_ERI_Hxy3z_S_D2x_S_vrr = QCX*I_ERI_Hxy3z_S_Px_S_vrr+WQX*I_ERI_Hxy3z_S_Px_S_M1_vrr+oned2e*I_ERI_Hxy3z_S_S_S_vrr-rhod2esq*I_ERI_Hxy3z_S_S_S_M1_vrr+oned2k*I_ERI_Gy3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_D2x_S_vrr = QCX*I_ERI_Hx4z_S_Px_S_vrr+WQX*I_ERI_Hx4z_S_Px_S_M1_vrr+oned2e*I_ERI_Hx4z_S_S_S_vrr-rhod2esq*I_ERI_Hx4z_S_S_S_M1_vrr+oned2k*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5y_S_D2x_S_vrr = QCX*I_ERI_H5y_S_Px_S_vrr+WQX*I_ERI_H5y_S_Px_S_M1_vrr+oned2e*I_ERI_H5y_S_S_S_vrr-rhod2esq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_D2x_S_vrr = QCX*I_ERI_H4yz_S_Px_S_vrr+WQX*I_ERI_H4yz_S_Px_S_M1_vrr+oned2e*I_ERI_H4yz_S_S_S_vrr-rhod2esq*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_D2x_S_vrr = QCX*I_ERI_H3y2z_S_Px_S_vrr+WQX*I_ERI_H3y2z_S_Px_S_M1_vrr+oned2e*I_ERI_H3y2z_S_S_S_vrr-rhod2esq*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_D2x_S_vrr = QCX*I_ERI_H2y3z_S_Px_S_vrr+WQX*I_ERI_H2y3z_S_Px_S_M1_vrr+oned2e*I_ERI_H2y3z_S_S_S_vrr-rhod2esq*I_ERI_H2y3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_D2x_S_vrr = QCX*I_ERI_Hy4z_S_Px_S_vrr+WQX*I_ERI_Hy4z_S_Px_S_M1_vrr+oned2e*I_ERI_Hy4z_S_S_S_vrr-rhod2esq*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_D2x_S_vrr = QCX*I_ERI_H5z_S_Px_S_vrr+WQX*I_ERI_H5z_S_Px_S_M1_vrr+oned2e*I_ERI_H5z_S_S_S_vrr-rhod2esq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_H5x_S_Dxy_S_vrr = QCY*I_ERI_H5x_S_Px_S_vrr+WQY*I_ERI_H5x_S_Px_S_M1_vrr;
      Double I_ERI_H4xy_S_Dxy_S_vrr = QCY*I_ERI_H4xy_S_Px_S_vrr+WQY*I_ERI_H4xy_S_Px_S_M1_vrr+oned2k*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_Dxy_S_vrr = QCY*I_ERI_H4xz_S_Px_S_vrr+WQY*I_ERI_H4xz_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_Dxy_S_vrr = QCY*I_ERI_H3x2y_S_Px_S_vrr+WQY*I_ERI_H3x2y_S_Px_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_Dxy_S_vrr = QCY*I_ERI_H3xyz_S_Px_S_vrr+WQY*I_ERI_H3xyz_S_Px_S_M1_vrr+oned2k*I_ERI_G3xz_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_Dxy_S_vrr = QCY*I_ERI_H3x2z_S_Px_S_vrr+WQY*I_ERI_H3x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_Dxy_S_vrr = QCY*I_ERI_H2x3y_S_Px_S_vrr+WQY*I_ERI_H2x3y_S_Px_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Dxy_S_vrr = QCY*I_ERI_H2x2yz_S_Px_S_vrr+WQY*I_ERI_H2x2yz_S_Px_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Dxy_S_vrr = QCY*I_ERI_H2xy2z_S_Px_S_vrr+WQY*I_ERI_H2xy2z_S_Px_S_M1_vrr+oned2k*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_Dxy_S_vrr = QCY*I_ERI_H2x3z_S_Px_S_vrr+WQY*I_ERI_H2x3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4y_S_Dxy_S_vrr = QCY*I_ERI_Hx4y_S_Px_S_vrr+WQY*I_ERI_Hx4y_S_Px_S_M1_vrr+4*oned2k*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Dxy_S_vrr = QCY*I_ERI_Hx3yz_S_Px_S_vrr+WQY*I_ERI_Hx3yz_S_Px_S_M1_vrr+3*oned2k*I_ERI_Gx2yz_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Dxy_S_vrr = QCY*I_ERI_Hx2y2z_S_Px_S_vrr+WQY*I_ERI_Hx2y2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Gxy2z_S_Px_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Dxy_S_vrr = QCY*I_ERI_Hxy3z_S_Px_S_vrr+WQY*I_ERI_Hxy3z_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_Dxy_S_vrr = QCY*I_ERI_Hx4z_S_Px_S_vrr+WQY*I_ERI_Hx4z_S_Px_S_M1_vrr;
      Double I_ERI_H5y_S_Dxy_S_vrr = QCY*I_ERI_H5y_S_Px_S_vrr+WQY*I_ERI_H5y_S_Px_S_M1_vrr+5*oned2k*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_H4yz_S_Dxy_S_vrr = QCY*I_ERI_H4yz_S_Px_S_vrr+WQY*I_ERI_H4yz_S_Px_S_M1_vrr+4*oned2k*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_H3y2z_S_Dxy_S_vrr = QCY*I_ERI_H3y2z_S_Px_S_vrr+WQY*I_ERI_H3y2z_S_Px_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Px_S_M1_vrr;
      Double I_ERI_H2y3z_S_Dxy_S_vrr = QCY*I_ERI_H2y3z_S_Px_S_vrr+WQY*I_ERI_H2y3z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Gy3z_S_Px_S_M1_vrr;
      Double I_ERI_Hy4z_S_Dxy_S_vrr = QCY*I_ERI_Hy4z_S_Px_S_vrr+WQY*I_ERI_Hy4z_S_Px_S_M1_vrr+oned2k*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5z_S_Dxy_S_vrr = QCY*I_ERI_H5z_S_Px_S_vrr+WQY*I_ERI_H5z_S_Px_S_M1_vrr;
      Double I_ERI_H5x_S_Dxz_S_vrr = QCZ*I_ERI_H5x_S_Px_S_vrr+WQZ*I_ERI_H5x_S_Px_S_M1_vrr;
      Double I_ERI_H4xy_S_Dxz_S_vrr = QCZ*I_ERI_H4xy_S_Px_S_vrr+WQZ*I_ERI_H4xy_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_Dxz_S_vrr = QCZ*I_ERI_H4xz_S_Px_S_vrr+WQZ*I_ERI_H4xz_S_Px_S_M1_vrr+oned2k*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_Dxz_S_vrr = QCZ*I_ERI_H3x2y_S_Px_S_vrr+WQZ*I_ERI_H3x2y_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_Dxz_S_vrr = QCZ*I_ERI_H3xyz_S_Px_S_vrr+WQZ*I_ERI_H3xyz_S_Px_S_M1_vrr+oned2k*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_Dxz_S_vrr = QCZ*I_ERI_H3x2z_S_Px_S_vrr+WQZ*I_ERI_H3x2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_Dxz_S_vrr = QCZ*I_ERI_H2x3y_S_Px_S_vrr+WQZ*I_ERI_H2x3y_S_Px_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Dxz_S_vrr = QCZ*I_ERI_H2x2yz_S_Px_S_vrr+WQZ*I_ERI_H2x2yz_S_Px_S_M1_vrr+oned2k*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Dxz_S_vrr = QCZ*I_ERI_H2xy2z_S_Px_S_vrr+WQZ*I_ERI_H2xy2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_Dxz_S_vrr = QCZ*I_ERI_H2x3z_S_Px_S_vrr+WQZ*I_ERI_H2x3z_S_Px_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4y_S_Dxz_S_vrr = QCZ*I_ERI_Hx4y_S_Px_S_vrr+WQZ*I_ERI_Hx4y_S_Px_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Dxz_S_vrr = QCZ*I_ERI_Hx3yz_S_Px_S_vrr+WQZ*I_ERI_Hx3yz_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Dxz_S_vrr = QCZ*I_ERI_Hx2y2z_S_Px_S_vrr+WQZ*I_ERI_Hx2y2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_Px_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Dxz_S_vrr = QCZ*I_ERI_Hxy3z_S_Px_S_vrr+WQZ*I_ERI_Hxy3z_S_Px_S_M1_vrr+3*oned2k*I_ERI_Gxy2z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_Dxz_S_vrr = QCZ*I_ERI_Hx4z_S_Px_S_vrr+WQZ*I_ERI_Hx4z_S_Px_S_M1_vrr+4*oned2k*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_H5y_S_Dxz_S_vrr = QCZ*I_ERI_H5y_S_Px_S_vrr+WQZ*I_ERI_H5y_S_Px_S_M1_vrr;
      Double I_ERI_H4yz_S_Dxz_S_vrr = QCZ*I_ERI_H4yz_S_Px_S_vrr+WQZ*I_ERI_H4yz_S_Px_S_M1_vrr+oned2k*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_H3y2z_S_Dxz_S_vrr = QCZ*I_ERI_H3y2z_S_Px_S_vrr+WQZ*I_ERI_H3y2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_H2y3z_S_Dxz_S_vrr = QCZ*I_ERI_H2y3z_S_Px_S_vrr+WQZ*I_ERI_H2y3z_S_Px_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Px_S_M1_vrr;
      Double I_ERI_Hy4z_S_Dxz_S_vrr = QCZ*I_ERI_Hy4z_S_Px_S_vrr+WQZ*I_ERI_Hy4z_S_Px_S_M1_vrr+4*oned2k*I_ERI_Gy3z_S_Px_S_M1_vrr;
      Double I_ERI_H5z_S_Dxz_S_vrr = QCZ*I_ERI_H5z_S_Px_S_vrr+WQZ*I_ERI_H5z_S_Px_S_M1_vrr+5*oned2k*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5x_S_D2y_S_vrr = QCY*I_ERI_H5x_S_Py_S_vrr+WQY*I_ERI_H5x_S_Py_S_M1_vrr+oned2e*I_ERI_H5x_S_S_S_vrr-rhod2esq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_D2y_S_vrr = QCY*I_ERI_H4xy_S_Py_S_vrr+WQY*I_ERI_H4xy_S_Py_S_M1_vrr+oned2e*I_ERI_H4xy_S_S_S_vrr-rhod2esq*I_ERI_H4xy_S_S_S_M1_vrr+oned2k*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H4xz_S_D2y_S_vrr = QCY*I_ERI_H4xz_S_Py_S_vrr+WQY*I_ERI_H4xz_S_Py_S_M1_vrr+oned2e*I_ERI_H4xz_S_S_S_vrr-rhod2esq*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_D2y_S_vrr = QCY*I_ERI_H3x2y_S_Py_S_vrr+WQY*I_ERI_H3x2y_S_Py_S_M1_vrr+oned2e*I_ERI_H3x2y_S_S_S_vrr-rhod2esq*I_ERI_H3x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_H3xyz_S_D2y_S_vrr = QCY*I_ERI_H3xyz_S_Py_S_vrr+WQY*I_ERI_H3xyz_S_Py_S_M1_vrr+oned2e*I_ERI_H3xyz_S_S_S_vrr-rhod2esq*I_ERI_H3xyz_S_S_S_M1_vrr+oned2k*I_ERI_G3xz_S_Py_S_M1_vrr;
      Double I_ERI_H3x2z_S_D2y_S_vrr = QCY*I_ERI_H3x2z_S_Py_S_vrr+WQY*I_ERI_H3x2z_S_Py_S_M1_vrr+oned2e*I_ERI_H3x2z_S_S_S_vrr-rhod2esq*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_D2y_S_vrr = QCY*I_ERI_H2x3y_S_Py_S_vrr+WQY*I_ERI_H2x3y_S_Py_S_M1_vrr+oned2e*I_ERI_H2x3y_S_S_S_vrr-rhod2esq*I_ERI_H2x3y_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_H2x2yz_S_D2y_S_vrr = QCY*I_ERI_H2x2yz_S_Py_S_vrr+WQY*I_ERI_H2x2yz_S_Py_S_M1_vrr+oned2e*I_ERI_H2x2yz_S_S_S_vrr-rhod2esq*I_ERI_H2x2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Py_S_M1_vrr;
      Double I_ERI_H2xy2z_S_D2y_S_vrr = QCY*I_ERI_H2xy2z_S_Py_S_vrr+WQY*I_ERI_H2xy2z_S_Py_S_M1_vrr+oned2e*I_ERI_H2xy2z_S_S_S_vrr-rhod2esq*I_ERI_H2xy2z_S_S_S_M1_vrr+oned2k*I_ERI_G2x2z_S_Py_S_M1_vrr;
      Double I_ERI_H2x3z_S_D2y_S_vrr = QCY*I_ERI_H2x3z_S_Py_S_vrr+WQY*I_ERI_H2x3z_S_Py_S_M1_vrr+oned2e*I_ERI_H2x3z_S_S_S_vrr-rhod2esq*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_D2y_S_vrr = QCY*I_ERI_Hx4y_S_Py_S_vrr+WQY*I_ERI_Hx4y_S_Py_S_M1_vrr+oned2e*I_ERI_Hx4y_S_S_S_vrr-rhod2esq*I_ERI_Hx4y_S_S_S_M1_vrr+4*oned2k*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_Hx3yz_S_D2y_S_vrr = QCY*I_ERI_Hx3yz_S_Py_S_vrr+WQY*I_ERI_Hx3yz_S_Py_S_M1_vrr+oned2e*I_ERI_Hx3yz_S_S_S_vrr-rhod2esq*I_ERI_Hx3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_Gx2yz_S_Py_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_D2y_S_vrr = QCY*I_ERI_Hx2y2z_S_Py_S_vrr+WQY*I_ERI_Hx2y2z_S_Py_S_M1_vrr+oned2e*I_ERI_Hx2y2z_S_S_S_vrr-rhod2esq*I_ERI_Hx2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gxy2z_S_Py_S_M1_vrr;
      Double I_ERI_Hxy3z_S_D2y_S_vrr = QCY*I_ERI_Hxy3z_S_Py_S_vrr+WQY*I_ERI_Hxy3z_S_Py_S_M1_vrr+oned2e*I_ERI_Hxy3z_S_S_S_vrr-rhod2esq*I_ERI_Hxy3z_S_S_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4z_S_D2y_S_vrr = QCY*I_ERI_Hx4z_S_Py_S_vrr+WQY*I_ERI_Hx4z_S_Py_S_M1_vrr+oned2e*I_ERI_Hx4z_S_S_S_vrr-rhod2esq*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_D2y_S_vrr = QCY*I_ERI_H5y_S_Py_S_vrr+WQY*I_ERI_H5y_S_Py_S_M1_vrr+oned2e*I_ERI_H5y_S_S_S_vrr-rhod2esq*I_ERI_H5y_S_S_S_M1_vrr+5*oned2k*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_H4yz_S_D2y_S_vrr = QCY*I_ERI_H4yz_S_Py_S_vrr+WQY*I_ERI_H4yz_S_Py_S_M1_vrr+oned2e*I_ERI_H4yz_S_S_S_vrr-rhod2esq*I_ERI_H4yz_S_S_S_M1_vrr+4*oned2k*I_ERI_G3yz_S_Py_S_M1_vrr;
      Double I_ERI_H3y2z_S_D2y_S_vrr = QCY*I_ERI_H3y2z_S_Py_S_vrr+WQY*I_ERI_H3y2z_S_Py_S_M1_vrr+oned2e*I_ERI_H3y2z_S_S_S_vrr-rhod2esq*I_ERI_H3y2z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Py_S_M1_vrr;
      Double I_ERI_H2y3z_S_D2y_S_vrr = QCY*I_ERI_H2y3z_S_Py_S_vrr+WQY*I_ERI_H2y3z_S_Py_S_M1_vrr+oned2e*I_ERI_H2y3z_S_S_S_vrr-rhod2esq*I_ERI_H2y3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gy3z_S_Py_S_M1_vrr;
      Double I_ERI_Hy4z_S_D2y_S_vrr = QCY*I_ERI_Hy4z_S_Py_S_vrr+WQY*I_ERI_Hy4z_S_Py_S_M1_vrr+oned2e*I_ERI_Hy4z_S_S_S_vrr-rhod2esq*I_ERI_Hy4z_S_S_S_M1_vrr+oned2k*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5z_S_D2y_S_vrr = QCY*I_ERI_H5z_S_Py_S_vrr+WQY*I_ERI_H5z_S_Py_S_M1_vrr+oned2e*I_ERI_H5z_S_S_S_vrr-rhod2esq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_H5x_S_Dyz_S_vrr = QCZ*I_ERI_H5x_S_Py_S_vrr+WQZ*I_ERI_H5x_S_Py_S_M1_vrr;
      Double I_ERI_H4xy_S_Dyz_S_vrr = QCZ*I_ERI_H4xy_S_Py_S_vrr+WQZ*I_ERI_H4xy_S_Py_S_M1_vrr;
      Double I_ERI_H4xz_S_Dyz_S_vrr = QCZ*I_ERI_H4xz_S_Py_S_vrr+WQZ*I_ERI_H4xz_S_Py_S_M1_vrr+oned2k*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H3x2y_S_Dyz_S_vrr = QCZ*I_ERI_H3x2y_S_Py_S_vrr+WQZ*I_ERI_H3x2y_S_Py_S_M1_vrr;
      Double I_ERI_H3xyz_S_Dyz_S_vrr = QCZ*I_ERI_H3xyz_S_Py_S_vrr+WQZ*I_ERI_H3xyz_S_Py_S_M1_vrr+oned2k*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_H3x2z_S_Dyz_S_vrr = QCZ*I_ERI_H3x2z_S_Py_S_vrr+WQZ*I_ERI_H3x2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_Py_S_M1_vrr;
      Double I_ERI_H2x3y_S_Dyz_S_vrr = QCZ*I_ERI_H2x3y_S_Py_S_vrr+WQZ*I_ERI_H2x3y_S_Py_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Dyz_S_vrr = QCZ*I_ERI_H2x2yz_S_Py_S_vrr+WQZ*I_ERI_H2x2yz_S_Py_S_M1_vrr+oned2k*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Dyz_S_vrr = QCZ*I_ERI_H2xy2z_S_Py_S_vrr+WQZ*I_ERI_H2xy2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Py_S_M1_vrr;
      Double I_ERI_H2x3z_S_Dyz_S_vrr = QCZ*I_ERI_H2x3z_S_Py_S_vrr+WQZ*I_ERI_H2x3z_S_Py_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4y_S_Dyz_S_vrr = QCZ*I_ERI_Hx4y_S_Py_S_vrr+WQZ*I_ERI_Hx4y_S_Py_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Dyz_S_vrr = QCZ*I_ERI_Hx3yz_S_Py_S_vrr+WQZ*I_ERI_Hx3yz_S_Py_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Dyz_S_vrr = QCZ*I_ERI_Hx2y2z_S_Py_S_vrr+WQZ*I_ERI_Hx2y2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_Py_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Dyz_S_vrr = QCZ*I_ERI_Hxy3z_S_Py_S_vrr+WQZ*I_ERI_Hxy3z_S_Py_S_M1_vrr+3*oned2k*I_ERI_Gxy2z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4z_S_Dyz_S_vrr = QCZ*I_ERI_Hx4z_S_Py_S_vrr+WQZ*I_ERI_Hx4z_S_Py_S_M1_vrr+4*oned2k*I_ERI_Gx3z_S_Py_S_M1_vrr;
      Double I_ERI_H5y_S_Dyz_S_vrr = QCZ*I_ERI_H5y_S_Py_S_vrr+WQZ*I_ERI_H5y_S_Py_S_M1_vrr;
      Double I_ERI_H4yz_S_Dyz_S_vrr = QCZ*I_ERI_H4yz_S_Py_S_vrr+WQZ*I_ERI_H4yz_S_Py_S_M1_vrr+oned2k*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_H3y2z_S_Dyz_S_vrr = QCZ*I_ERI_H3y2z_S_Py_S_vrr+WQZ*I_ERI_H3y2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_Py_S_M1_vrr;
      Double I_ERI_H2y3z_S_Dyz_S_vrr = QCZ*I_ERI_H2y3z_S_Py_S_vrr+WQZ*I_ERI_H2y3z_S_Py_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Py_S_M1_vrr;
      Double I_ERI_Hy4z_S_Dyz_S_vrr = QCZ*I_ERI_Hy4z_S_Py_S_vrr+WQZ*I_ERI_Hy4z_S_Py_S_M1_vrr+4*oned2k*I_ERI_Gy3z_S_Py_S_M1_vrr;
      Double I_ERI_H5z_S_Dyz_S_vrr = QCZ*I_ERI_H5z_S_Py_S_vrr+WQZ*I_ERI_H5z_S_Py_S_M1_vrr+5*oned2k*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5x_S_D2z_S_vrr = QCZ*I_ERI_H5x_S_Pz_S_vrr+WQZ*I_ERI_H5x_S_Pz_S_M1_vrr+oned2e*I_ERI_H5x_S_S_S_vrr-rhod2esq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_D2z_S_vrr = QCZ*I_ERI_H4xy_S_Pz_S_vrr+WQZ*I_ERI_H4xy_S_Pz_S_M1_vrr+oned2e*I_ERI_H4xy_S_S_S_vrr-rhod2esq*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_D2z_S_vrr = QCZ*I_ERI_H4xz_S_Pz_S_vrr+WQZ*I_ERI_H4xz_S_Pz_S_M1_vrr+oned2e*I_ERI_H4xz_S_S_S_vrr-rhod2esq*I_ERI_H4xz_S_S_S_M1_vrr+oned2k*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_H3x2y_S_D2z_S_vrr = QCZ*I_ERI_H3x2y_S_Pz_S_vrr+WQZ*I_ERI_H3x2y_S_Pz_S_M1_vrr+oned2e*I_ERI_H3x2y_S_S_S_vrr-rhod2esq*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_D2z_S_vrr = QCZ*I_ERI_H3xyz_S_Pz_S_vrr+WQZ*I_ERI_H3xyz_S_Pz_S_M1_vrr+oned2e*I_ERI_H3xyz_S_S_S_vrr-rhod2esq*I_ERI_H3xyz_S_S_S_M1_vrr+oned2k*I_ERI_G3xy_S_Pz_S_M1_vrr;
      Double I_ERI_H3x2z_S_D2z_S_vrr = QCZ*I_ERI_H3x2z_S_Pz_S_vrr+WQZ*I_ERI_H3x2z_S_Pz_S_M1_vrr+oned2e*I_ERI_H3x2z_S_S_S_vrr-rhod2esq*I_ERI_H3x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_Pz_S_M1_vrr;
      Double I_ERI_H2x3y_S_D2z_S_vrr = QCZ*I_ERI_H2x3y_S_Pz_S_vrr+WQZ*I_ERI_H2x3y_S_Pz_S_M1_vrr+oned2e*I_ERI_H2x3y_S_S_S_vrr-rhod2esq*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_D2z_S_vrr = QCZ*I_ERI_H2x2yz_S_Pz_S_vrr+WQZ*I_ERI_H2x2yz_S_Pz_S_M1_vrr+oned2e*I_ERI_H2x2yz_S_S_S_vrr-rhod2esq*I_ERI_H2x2yz_S_S_S_M1_vrr+oned2k*I_ERI_G2x2y_S_Pz_S_M1_vrr;
      Double I_ERI_H2xy2z_S_D2z_S_vrr = QCZ*I_ERI_H2xy2z_S_Pz_S_vrr+WQZ*I_ERI_H2xy2z_S_Pz_S_M1_vrr+oned2e*I_ERI_H2xy2z_S_S_S_vrr-rhod2esq*I_ERI_H2xy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Pz_S_M1_vrr;
      Double I_ERI_H2x3z_S_D2z_S_vrr = QCZ*I_ERI_H2x3z_S_Pz_S_vrr+WQZ*I_ERI_H2x3z_S_Pz_S_M1_vrr+oned2e*I_ERI_H2x3z_S_S_S_vrr-rhod2esq*I_ERI_H2x3z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4y_S_D2z_S_vrr = QCZ*I_ERI_Hx4y_S_Pz_S_vrr+WQZ*I_ERI_Hx4y_S_Pz_S_M1_vrr+oned2e*I_ERI_Hx4y_S_S_S_vrr-rhod2esq*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_D2z_S_vrr = QCZ*I_ERI_Hx3yz_S_Pz_S_vrr+WQZ*I_ERI_Hx3yz_S_Pz_S_M1_vrr+oned2e*I_ERI_Hx3yz_S_S_S_vrr-rhod2esq*I_ERI_Hx3yz_S_S_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Pz_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_D2z_S_vrr = QCZ*I_ERI_Hx2y2z_S_Pz_S_vrr+WQZ*I_ERI_Hx2y2z_S_Pz_S_M1_vrr+oned2e*I_ERI_Hx2y2z_S_S_S_vrr-rhod2esq*I_ERI_Hx2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_Pz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_D2z_S_vrr = QCZ*I_ERI_Hxy3z_S_Pz_S_vrr+WQZ*I_ERI_Hxy3z_S_Pz_S_M1_vrr+oned2e*I_ERI_Hxy3z_S_S_S_vrr-rhod2esq*I_ERI_Hxy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Gxy2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4z_S_D2z_S_vrr = QCZ*I_ERI_Hx4z_S_Pz_S_vrr+WQZ*I_ERI_Hx4z_S_Pz_S_M1_vrr+oned2e*I_ERI_Hx4z_S_S_S_vrr-rhod2esq*I_ERI_Hx4z_S_S_S_M1_vrr+4*oned2k*I_ERI_Gx3z_S_Pz_S_M1_vrr;
      Double I_ERI_H5y_S_D2z_S_vrr = QCZ*I_ERI_H5y_S_Pz_S_vrr+WQZ*I_ERI_H5y_S_Pz_S_M1_vrr+oned2e*I_ERI_H5y_S_S_S_vrr-rhod2esq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_D2z_S_vrr = QCZ*I_ERI_H4yz_S_Pz_S_vrr+WQZ*I_ERI_H4yz_S_Pz_S_M1_vrr+oned2e*I_ERI_H4yz_S_S_S_vrr-rhod2esq*I_ERI_H4yz_S_S_S_M1_vrr+oned2k*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_H3y2z_S_D2z_S_vrr = QCZ*I_ERI_H3y2z_S_Pz_S_vrr+WQZ*I_ERI_H3y2z_S_Pz_S_M1_vrr+oned2e*I_ERI_H3y2z_S_S_S_vrr-rhod2esq*I_ERI_H3y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_Pz_S_M1_vrr;
      Double I_ERI_H2y3z_S_D2z_S_vrr = QCZ*I_ERI_H2y3z_S_Pz_S_vrr+WQZ*I_ERI_H2y3z_S_Pz_S_M1_vrr+oned2e*I_ERI_H2y3z_S_S_S_vrr-rhod2esq*I_ERI_H2y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hy4z_S_D2z_S_vrr = QCZ*I_ERI_Hy4z_S_Pz_S_vrr+WQZ*I_ERI_Hy4z_S_Pz_S_M1_vrr+oned2e*I_ERI_Hy4z_S_S_S_vrr-rhod2esq*I_ERI_Hy4z_S_S_S_M1_vrr+4*oned2k*I_ERI_Gy3z_S_Pz_S_M1_vrr;
      Double I_ERI_H5z_S_D2z_S_vrr = QCZ*I_ERI_H5z_S_Pz_S_vrr+WQZ*I_ERI_H5z_S_Pz_S_M1_vrr+oned2e*I_ERI_H5z_S_S_S_vrr-rhod2esq*I_ERI_H5z_S_S_S_M1_vrr+5*oned2k*I_ERI_G4z_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_P_S
       * RHS shell quartet name: SQ_ERI_H_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_P_S
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       ************************************************************/
      Double I_ERI_I6x_S_Px_S_vrr = PAX*I_ERI_H5x_S_Px_S_vrr+WPX*I_ERI_H5x_S_Px_S_M1_vrr+5*oned2z*I_ERI_G4x_S_Px_S_vrr-5*rhod2zsq*I_ERI_G4x_S_Px_S_M1_vrr+oned2k*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I5xy_S_Px_S_vrr = PAY*I_ERI_H5x_S_Px_S_vrr+WPY*I_ERI_H5x_S_Px_S_M1_vrr;
      Double I_ERI_I5xz_S_Px_S_vrr = PAZ*I_ERI_H5x_S_Px_S_vrr+WPZ*I_ERI_H5x_S_Px_S_M1_vrr;
      Double I_ERI_I4x2y_S_Px_S_vrr = PAY*I_ERI_H4xy_S_Px_S_vrr+WPY*I_ERI_H4xy_S_Px_S_M1_vrr+oned2z*I_ERI_G4x_S_Px_S_vrr-rhod2zsq*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_I4xyz_S_Px_S_vrr = PAZ*I_ERI_H4xy_S_Px_S_vrr+WPZ*I_ERI_H4xy_S_Px_S_M1_vrr;
      Double I_ERI_I4x2z_S_Px_S_vrr = PAZ*I_ERI_H4xz_S_Px_S_vrr+WPZ*I_ERI_H4xz_S_Px_S_M1_vrr+oned2z*I_ERI_G4x_S_Px_S_vrr-rhod2zsq*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_I3x3y_S_Px_S_vrr = PAY*I_ERI_H3x2y_S_Px_S_vrr+WPY*I_ERI_H3x2y_S_Px_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_Px_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Px_S_vrr = PAZ*I_ERI_H3x2y_S_Px_S_vrr+WPZ*I_ERI_H3x2y_S_Px_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Px_S_vrr = PAY*I_ERI_H3x2z_S_Px_S_vrr+WPY*I_ERI_H3x2z_S_Px_S_M1_vrr;
      Double I_ERI_I3x3z_S_Px_S_vrr = PAZ*I_ERI_H3x2z_S_Px_S_vrr+WPZ*I_ERI_H3x2z_S_Px_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_Px_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_Px_S_M1_vrr;
      Double I_ERI_I2x4y_S_Px_S_vrr = PAX*I_ERI_Hx4y_S_Px_S_vrr+WPX*I_ERI_Hx4y_S_Px_S_M1_vrr+oned2z*I_ERI_G4y_S_Px_S_vrr-rhod2zsq*I_ERI_G4y_S_Px_S_M1_vrr+oned2k*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Px_S_vrr = PAZ*I_ERI_H2x3y_S_Px_S_vrr+WPZ*I_ERI_H2x3y_S_Px_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Px_S_vrr = PAZ*I_ERI_H2x2yz_S_Px_S_vrr+WPZ*I_ERI_H2x2yz_S_Px_S_M1_vrr+oned2z*I_ERI_G2x2y_S_Px_S_vrr-rhod2zsq*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Px_S_vrr = PAY*I_ERI_H2x3z_S_Px_S_vrr+WPY*I_ERI_H2x3z_S_Px_S_M1_vrr;
      Double I_ERI_I2x4z_S_Px_S_vrr = PAX*I_ERI_Hx4z_S_Px_S_vrr+WPX*I_ERI_Hx4z_S_Px_S_M1_vrr+oned2z*I_ERI_G4z_S_Px_S_vrr-rhod2zsq*I_ERI_G4z_S_Px_S_M1_vrr+oned2k*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5y_S_Px_S_vrr = PAX*I_ERI_H5y_S_Px_S_vrr+WPX*I_ERI_H5y_S_Px_S_M1_vrr+oned2k*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Px_S_vrr = PAZ*I_ERI_Hx4y_S_Px_S_vrr+WPZ*I_ERI_Hx4y_S_Px_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Px_S_vrr = PAX*I_ERI_H3y2z_S_Px_S_vrr+WPX*I_ERI_H3y2z_S_Px_S_M1_vrr+oned2k*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Px_S_vrr = PAX*I_ERI_H2y3z_S_Px_S_vrr+WPX*I_ERI_H2y3z_S_Px_S_M1_vrr+oned2k*I_ERI_H2y3z_S_S_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Px_S_vrr = PAY*I_ERI_Hx4z_S_Px_S_vrr+WPY*I_ERI_Hx4z_S_Px_S_M1_vrr;
      Double I_ERI_Ix5z_S_Px_S_vrr = PAX*I_ERI_H5z_S_Px_S_vrr+WPX*I_ERI_H5z_S_Px_S_M1_vrr+oned2k*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_I6y_S_Px_S_vrr = PAY*I_ERI_H5y_S_Px_S_vrr+WPY*I_ERI_H5y_S_Px_S_M1_vrr+5*oned2z*I_ERI_G4y_S_Px_S_vrr-5*rhod2zsq*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_I5yz_S_Px_S_vrr = PAZ*I_ERI_H5y_S_Px_S_vrr+WPZ*I_ERI_H5y_S_Px_S_M1_vrr;
      Double I_ERI_I4y2z_S_Px_S_vrr = PAZ*I_ERI_H4yz_S_Px_S_vrr+WPZ*I_ERI_H4yz_S_Px_S_M1_vrr+oned2z*I_ERI_G4y_S_Px_S_vrr-rhod2zsq*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_I3y3z_S_Px_S_vrr = PAZ*I_ERI_H3y2z_S_Px_S_vrr+WPZ*I_ERI_H3y2z_S_Px_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_Px_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_I2y4z_S_Px_S_vrr = PAY*I_ERI_Hy4z_S_Px_S_vrr+WPY*I_ERI_Hy4z_S_Px_S_M1_vrr+oned2z*I_ERI_G4z_S_Px_S_vrr-rhod2zsq*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_Iy5z_S_Px_S_vrr = PAY*I_ERI_H5z_S_Px_S_vrr+WPY*I_ERI_H5z_S_Px_S_M1_vrr;
      Double I_ERI_I6z_S_Px_S_vrr = PAZ*I_ERI_H5z_S_Px_S_vrr+WPZ*I_ERI_H5z_S_Px_S_M1_vrr+5*oned2z*I_ERI_G4z_S_Px_S_vrr-5*rhod2zsq*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_I6x_S_Py_S_vrr = PAX*I_ERI_H5x_S_Py_S_vrr+WPX*I_ERI_H5x_S_Py_S_M1_vrr+5*oned2z*I_ERI_G4x_S_Py_S_vrr-5*rhod2zsq*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_I5xy_S_Py_S_vrr = PAY*I_ERI_H5x_S_Py_S_vrr+WPY*I_ERI_H5x_S_Py_S_M1_vrr+oned2k*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I5xz_S_Py_S_vrr = PAZ*I_ERI_H5x_S_Py_S_vrr+WPZ*I_ERI_H5x_S_Py_S_M1_vrr;
      Double I_ERI_I4x2y_S_Py_S_vrr = PAY*I_ERI_H4xy_S_Py_S_vrr+WPY*I_ERI_H4xy_S_Py_S_M1_vrr+oned2z*I_ERI_G4x_S_Py_S_vrr-rhod2zsq*I_ERI_G4x_S_Py_S_M1_vrr+oned2k*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_I4xyz_S_Py_S_vrr = PAZ*I_ERI_H4xy_S_Py_S_vrr+WPZ*I_ERI_H4xy_S_Py_S_M1_vrr;
      Double I_ERI_I4x2z_S_Py_S_vrr = PAZ*I_ERI_H4xz_S_Py_S_vrr+WPZ*I_ERI_H4xz_S_Py_S_M1_vrr+oned2z*I_ERI_G4x_S_Py_S_vrr-rhod2zsq*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_I3x3y_S_Py_S_vrr = PAY*I_ERI_H3x2y_S_Py_S_vrr+WPY*I_ERI_H3x2y_S_Py_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_Py_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_Py_S_M1_vrr+oned2k*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Py_S_vrr = PAZ*I_ERI_H3x2y_S_Py_S_vrr+WPZ*I_ERI_H3x2y_S_Py_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Py_S_vrr = PAY*I_ERI_H3x2z_S_Py_S_vrr+WPY*I_ERI_H3x2z_S_Py_S_M1_vrr+oned2k*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_I3x3z_S_Py_S_vrr = PAZ*I_ERI_H3x2z_S_Py_S_vrr+WPZ*I_ERI_H3x2z_S_Py_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_Py_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_Py_S_M1_vrr;
      Double I_ERI_I2x4y_S_Py_S_vrr = PAX*I_ERI_Hx4y_S_Py_S_vrr+WPX*I_ERI_Hx4y_S_Py_S_M1_vrr+oned2z*I_ERI_G4y_S_Py_S_vrr-rhod2zsq*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Py_S_vrr = PAZ*I_ERI_H2x3y_S_Py_S_vrr+WPZ*I_ERI_H2x3y_S_Py_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Py_S_vrr = PAZ*I_ERI_H2x2yz_S_Py_S_vrr+WPZ*I_ERI_H2x2yz_S_Py_S_M1_vrr+oned2z*I_ERI_G2x2y_S_Py_S_vrr-rhod2zsq*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Py_S_vrr = PAY*I_ERI_H2x3z_S_Py_S_vrr+WPY*I_ERI_H2x3z_S_Py_S_M1_vrr+oned2k*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_I2x4z_S_Py_S_vrr = PAX*I_ERI_Hx4z_S_Py_S_vrr+WPX*I_ERI_Hx4z_S_Py_S_M1_vrr+oned2z*I_ERI_G4z_S_Py_S_vrr-rhod2zsq*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_Ix5y_S_Py_S_vrr = PAX*I_ERI_H5y_S_Py_S_vrr+WPX*I_ERI_H5y_S_Py_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Py_S_vrr = PAZ*I_ERI_Hx4y_S_Py_S_vrr+WPZ*I_ERI_Hx4y_S_Py_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Py_S_vrr = PAX*I_ERI_H3y2z_S_Py_S_vrr+WPX*I_ERI_H3y2z_S_Py_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Py_S_vrr = PAX*I_ERI_H2y3z_S_Py_S_vrr+WPX*I_ERI_H2y3z_S_Py_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Py_S_vrr = PAY*I_ERI_Hx4z_S_Py_S_vrr+WPY*I_ERI_Hx4z_S_Py_S_M1_vrr+oned2k*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5z_S_Py_S_vrr = PAX*I_ERI_H5z_S_Py_S_vrr+WPX*I_ERI_H5z_S_Py_S_M1_vrr;
      Double I_ERI_I6y_S_Py_S_vrr = PAY*I_ERI_H5y_S_Py_S_vrr+WPY*I_ERI_H5y_S_Py_S_M1_vrr+5*oned2z*I_ERI_G4y_S_Py_S_vrr-5*rhod2zsq*I_ERI_G4y_S_Py_S_M1_vrr+oned2k*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_I5yz_S_Py_S_vrr = PAZ*I_ERI_H5y_S_Py_S_vrr+WPZ*I_ERI_H5y_S_Py_S_M1_vrr;
      Double I_ERI_I4y2z_S_Py_S_vrr = PAZ*I_ERI_H4yz_S_Py_S_vrr+WPZ*I_ERI_H4yz_S_Py_S_M1_vrr+oned2z*I_ERI_G4y_S_Py_S_vrr-rhod2zsq*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_I3y3z_S_Py_S_vrr = PAZ*I_ERI_H3y2z_S_Py_S_vrr+WPZ*I_ERI_H3y2z_S_Py_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_Py_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_Py_S_M1_vrr;
      Double I_ERI_I2y4z_S_Py_S_vrr = PAY*I_ERI_Hy4z_S_Py_S_vrr+WPY*I_ERI_Hy4z_S_Py_S_M1_vrr+oned2z*I_ERI_G4z_S_Py_S_vrr-rhod2zsq*I_ERI_G4z_S_Py_S_M1_vrr+oned2k*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_Iy5z_S_Py_S_vrr = PAY*I_ERI_H5z_S_Py_S_vrr+WPY*I_ERI_H5z_S_Py_S_M1_vrr+oned2k*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_I6z_S_Py_S_vrr = PAZ*I_ERI_H5z_S_Py_S_vrr+WPZ*I_ERI_H5z_S_Py_S_M1_vrr+5*oned2z*I_ERI_G4z_S_Py_S_vrr-5*rhod2zsq*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_I6x_S_Pz_S_vrr = PAX*I_ERI_H5x_S_Pz_S_vrr+WPX*I_ERI_H5x_S_Pz_S_M1_vrr+5*oned2z*I_ERI_G4x_S_Pz_S_vrr-5*rhod2zsq*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_I5xy_S_Pz_S_vrr = PAY*I_ERI_H5x_S_Pz_S_vrr+WPY*I_ERI_H5x_S_Pz_S_M1_vrr;
      Double I_ERI_I5xz_S_Pz_S_vrr = PAZ*I_ERI_H5x_S_Pz_S_vrr+WPZ*I_ERI_H5x_S_Pz_S_M1_vrr+oned2k*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I4x2y_S_Pz_S_vrr = PAY*I_ERI_H4xy_S_Pz_S_vrr+WPY*I_ERI_H4xy_S_Pz_S_M1_vrr+oned2z*I_ERI_G4x_S_Pz_S_vrr-rhod2zsq*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_I4xyz_S_Pz_S_vrr = PAZ*I_ERI_H4xy_S_Pz_S_vrr+WPZ*I_ERI_H4xy_S_Pz_S_M1_vrr+oned2k*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_I4x2z_S_Pz_S_vrr = PAZ*I_ERI_H4xz_S_Pz_S_vrr+WPZ*I_ERI_H4xz_S_Pz_S_M1_vrr+oned2z*I_ERI_G4x_S_Pz_S_vrr-rhod2zsq*I_ERI_G4x_S_Pz_S_M1_vrr+oned2k*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_I3x3y_S_Pz_S_vrr = PAY*I_ERI_H3x2y_S_Pz_S_vrr+WPY*I_ERI_H3x2y_S_Pz_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_Pz_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_Pz_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Pz_S_vrr = PAZ*I_ERI_H3x2y_S_Pz_S_vrr+WPZ*I_ERI_H3x2y_S_Pz_S_M1_vrr+oned2k*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Pz_S_vrr = PAY*I_ERI_H3x2z_S_Pz_S_vrr+WPY*I_ERI_H3x2z_S_Pz_S_M1_vrr;
      Double I_ERI_I3x3z_S_Pz_S_vrr = PAZ*I_ERI_H3x2z_S_Pz_S_vrr+WPZ*I_ERI_H3x2z_S_Pz_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_Pz_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_Pz_S_M1_vrr+oned2k*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_I2x4y_S_Pz_S_vrr = PAX*I_ERI_Hx4y_S_Pz_S_vrr+WPX*I_ERI_Hx4y_S_Pz_S_M1_vrr+oned2z*I_ERI_G4y_S_Pz_S_vrr-rhod2zsq*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Pz_S_vrr = PAZ*I_ERI_H2x3y_S_Pz_S_vrr+WPZ*I_ERI_H2x3y_S_Pz_S_M1_vrr+oned2k*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Pz_S_vrr = PAZ*I_ERI_H2x2yz_S_Pz_S_vrr+WPZ*I_ERI_H2x2yz_S_Pz_S_M1_vrr+oned2z*I_ERI_G2x2y_S_Pz_S_vrr-rhod2zsq*I_ERI_G2x2y_S_Pz_S_M1_vrr+oned2k*I_ERI_H2x2yz_S_S_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Pz_S_vrr = PAY*I_ERI_H2x3z_S_Pz_S_vrr+WPY*I_ERI_H2x3z_S_Pz_S_M1_vrr;
      Double I_ERI_I2x4z_S_Pz_S_vrr = PAX*I_ERI_Hx4z_S_Pz_S_vrr+WPX*I_ERI_Hx4z_S_Pz_S_M1_vrr+oned2z*I_ERI_G4z_S_Pz_S_vrr-rhod2zsq*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_Ix5y_S_Pz_S_vrr = PAX*I_ERI_H5y_S_Pz_S_vrr+WPX*I_ERI_H5y_S_Pz_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Pz_S_vrr = PAZ*I_ERI_Hx4y_S_Pz_S_vrr+WPZ*I_ERI_Hx4y_S_Pz_S_M1_vrr+oned2k*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Pz_S_vrr = PAX*I_ERI_H3y2z_S_Pz_S_vrr+WPX*I_ERI_H3y2z_S_Pz_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Pz_S_vrr = PAX*I_ERI_H2y3z_S_Pz_S_vrr+WPX*I_ERI_H2y3z_S_Pz_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Pz_S_vrr = PAY*I_ERI_Hx4z_S_Pz_S_vrr+WPY*I_ERI_Hx4z_S_Pz_S_M1_vrr;
      Double I_ERI_Ix5z_S_Pz_S_vrr = PAX*I_ERI_H5z_S_Pz_S_vrr+WPX*I_ERI_H5z_S_Pz_S_M1_vrr;
      Double I_ERI_I6y_S_Pz_S_vrr = PAY*I_ERI_H5y_S_Pz_S_vrr+WPY*I_ERI_H5y_S_Pz_S_M1_vrr+5*oned2z*I_ERI_G4y_S_Pz_S_vrr-5*rhod2zsq*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_I5yz_S_Pz_S_vrr = PAZ*I_ERI_H5y_S_Pz_S_vrr+WPZ*I_ERI_H5y_S_Pz_S_M1_vrr+oned2k*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_I4y2z_S_Pz_S_vrr = PAZ*I_ERI_H4yz_S_Pz_S_vrr+WPZ*I_ERI_H4yz_S_Pz_S_M1_vrr+oned2z*I_ERI_G4y_S_Pz_S_vrr-rhod2zsq*I_ERI_G4y_S_Pz_S_M1_vrr+oned2k*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_I3y3z_S_Pz_S_vrr = PAZ*I_ERI_H3y2z_S_Pz_S_vrr+WPZ*I_ERI_H3y2z_S_Pz_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_Pz_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_Pz_S_M1_vrr+oned2k*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_I2y4z_S_Pz_S_vrr = PAY*I_ERI_Hy4z_S_Pz_S_vrr+WPY*I_ERI_Hy4z_S_Pz_S_M1_vrr+oned2z*I_ERI_G4z_S_Pz_S_vrr-rhod2zsq*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_Iy5z_S_Pz_S_vrr = PAY*I_ERI_H5z_S_Pz_S_vrr+WPY*I_ERI_H5z_S_Pz_S_M1_vrr;
      Double I_ERI_I6z_S_Pz_S_vrr = PAZ*I_ERI_H5z_S_Pz_S_vrr+WPZ*I_ERI_H5z_S_Pz_S_M1_vrr+5*oned2z*I_ERI_G4z_S_Pz_S_vrr-5*rhod2zsq*I_ERI_G4z_S_Pz_S_M1_vrr+oned2k*I_ERI_H5z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_I_S_S_S
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       ************************************************************/
      Double I_ERI_K7x_S_S_S_vrr = PAX*I_ERI_I6x_S_S_S_vrr+WPX*I_ERI_I6x_S_S_S_M1_vrr+6*oned2z*I_ERI_H5x_S_S_S_vrr-6*rhod2zsq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_K6xy_S_S_S_vrr = PAY*I_ERI_I6x_S_S_S_vrr+WPY*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_K6xz_S_S_S_vrr = PAZ*I_ERI_I6x_S_S_S_vrr+WPZ*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_K5x2y_S_S_S_vrr = PAY*I_ERI_I5xy_S_S_S_vrr+WPY*I_ERI_I5xy_S_S_S_M1_vrr+oned2z*I_ERI_H5x_S_S_S_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_K5xyz_S_S_S_vrr = PAZ*I_ERI_I5xy_S_S_S_vrr+WPZ*I_ERI_I5xy_S_S_S_M1_vrr;
      Double I_ERI_K5x2z_S_S_S_vrr = PAZ*I_ERI_I5xz_S_S_S_vrr+WPZ*I_ERI_I5xz_S_S_S_M1_vrr+oned2z*I_ERI_H5x_S_S_S_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_K4x3y_S_S_S_vrr = PAY*I_ERI_I4x2y_S_S_S_vrr+WPY*I_ERI_I4x2y_S_S_S_M1_vrr+2*oned2z*I_ERI_H4xy_S_S_S_vrr-2*rhod2zsq*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_K4x2yz_S_S_S_vrr = PAZ*I_ERI_I4x2y_S_S_S_vrr+WPZ*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_K4xy2z_S_S_S_vrr = PAY*I_ERI_I4x2z_S_S_S_vrr+WPY*I_ERI_I4x2z_S_S_S_M1_vrr;
      Double I_ERI_K4x3z_S_S_S_vrr = PAZ*I_ERI_I4x2z_S_S_S_vrr+WPZ*I_ERI_I4x2z_S_S_S_M1_vrr+2*oned2z*I_ERI_H4xz_S_S_S_vrr-2*rhod2zsq*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_K3x4y_S_S_S_vrr = PAX*I_ERI_I2x4y_S_S_S_vrr+WPX*I_ERI_I2x4y_S_S_S_M1_vrr+2*oned2z*I_ERI_Hx4y_S_S_S_vrr-2*rhod2zsq*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_K3x3yz_S_S_S_vrr = PAZ*I_ERI_I3x3y_S_S_S_vrr+WPZ*I_ERI_I3x3y_S_S_S_M1_vrr;
      Double I_ERI_K3x2y2z_S_S_S_vrr = PAZ*I_ERI_I3x2yz_S_S_S_vrr+WPZ*I_ERI_I3x2yz_S_S_S_M1_vrr+oned2z*I_ERI_H3x2y_S_S_S_vrr-rhod2zsq*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_K3xy3z_S_S_S_vrr = PAY*I_ERI_I3x3z_S_S_S_vrr+WPY*I_ERI_I3x3z_S_S_S_M1_vrr;
      Double I_ERI_K3x4z_S_S_S_vrr = PAX*I_ERI_I2x4z_S_S_S_vrr+WPX*I_ERI_I2x4z_S_S_S_M1_vrr+2*oned2z*I_ERI_Hx4z_S_S_S_vrr-2*rhod2zsq*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5y_S_S_S_vrr = PAX*I_ERI_Ix5y_S_S_S_vrr+WPX*I_ERI_Ix5y_S_S_S_M1_vrr+oned2z*I_ERI_H5y_S_S_S_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_K2x4yz_S_S_S_vrr = PAZ*I_ERI_I2x4y_S_S_S_vrr+WPZ*I_ERI_I2x4y_S_S_S_M1_vrr;
      Double I_ERI_K2x3y2z_S_S_S_vrr = PAZ*I_ERI_I2x3yz_S_S_S_vrr+WPZ*I_ERI_I2x3yz_S_S_S_M1_vrr+oned2z*I_ERI_H2x3y_S_S_S_vrr-rhod2zsq*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_K2x2y3z_S_S_S_vrr = PAY*I_ERI_I2xy3z_S_S_S_vrr+WPY*I_ERI_I2xy3z_S_S_S_M1_vrr+oned2z*I_ERI_H2x3z_S_S_S_vrr-rhod2zsq*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_K2xy4z_S_S_S_vrr = PAY*I_ERI_I2x4z_S_S_S_vrr+WPY*I_ERI_I2x4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5z_S_S_S_vrr = PAX*I_ERI_Ix5z_S_S_S_vrr+WPX*I_ERI_Ix5z_S_S_S_M1_vrr+oned2z*I_ERI_H5z_S_S_S_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6y_S_S_S_vrr = PAX*I_ERI_I6y_S_S_S_vrr+WPX*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_Kx5yz_S_S_S_vrr = PAZ*I_ERI_Ix5y_S_S_S_vrr+WPZ*I_ERI_Ix5y_S_S_S_M1_vrr;
      Double I_ERI_Kx4y2z_S_S_S_vrr = PAX*I_ERI_I4y2z_S_S_S_vrr+WPX*I_ERI_I4y2z_S_S_S_M1_vrr;
      Double I_ERI_Kx3y3z_S_S_S_vrr = PAX*I_ERI_I3y3z_S_S_S_vrr+WPX*I_ERI_I3y3z_S_S_S_M1_vrr;
      Double I_ERI_Kx2y4z_S_S_S_vrr = PAX*I_ERI_I2y4z_S_S_S_vrr+WPX*I_ERI_I2y4z_S_S_S_M1_vrr;
      Double I_ERI_Kxy5z_S_S_S_vrr = PAY*I_ERI_Ix5z_S_S_S_vrr+WPY*I_ERI_Ix5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6z_S_S_S_vrr = PAX*I_ERI_I6z_S_S_S_vrr+WPX*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_K7y_S_S_S_vrr = PAY*I_ERI_I6y_S_S_S_vrr+WPY*I_ERI_I6y_S_S_S_M1_vrr+6*oned2z*I_ERI_H5y_S_S_S_vrr-6*rhod2zsq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_K6yz_S_S_S_vrr = PAZ*I_ERI_I6y_S_S_S_vrr+WPZ*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_K5y2z_S_S_S_vrr = PAZ*I_ERI_I5yz_S_S_S_vrr+WPZ*I_ERI_I5yz_S_S_S_M1_vrr+oned2z*I_ERI_H5y_S_S_S_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_K4y3z_S_S_S_vrr = PAZ*I_ERI_I4y2z_S_S_S_vrr+WPZ*I_ERI_I4y2z_S_S_S_M1_vrr+2*oned2z*I_ERI_H4yz_S_S_S_vrr-2*rhod2zsq*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_K3y4z_S_S_S_vrr = PAY*I_ERI_I2y4z_S_S_S_vrr+WPY*I_ERI_I2y4z_S_S_S_M1_vrr+2*oned2z*I_ERI_Hy4z_S_S_S_vrr-2*rhod2zsq*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_K2y5z_S_S_S_vrr = PAY*I_ERI_Iy5z_S_S_S_vrr+WPY*I_ERI_Iy5z_S_S_S_M1_vrr+oned2z*I_ERI_H5z_S_S_S_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_Ky6z_S_S_S_vrr = PAY*I_ERI_I6z_S_S_S_vrr+WPY*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_K7z_S_S_S_vrr = PAZ*I_ERI_I6z_S_S_S_vrr+WPZ*I_ERI_I6z_S_S_S_M1_vrr+6*oned2z*I_ERI_H5z_S_S_S_vrr-6*rhod2zsq*I_ERI_H5z_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_H_S_S_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_S_S_a_coefs = alpha;
      I_ERI_H5x_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S_a += SQ_ERI_H_S_S_S_a_coefs*I_ERI_H5z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_G_S_P_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_c_coefs = gamma;
      I_ERI_G4x_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_c += SQ_ERI_G_S_P_S_c_coefs*I_ERI_G4z_S_Pz_S_vrr;

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
       * shell quartet name: SQ_ERI_K_S_S_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_K_S_S_S_aa_coefs = alpha*alpha;
      I_ERI_K7x_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K7x_S_S_S_vrr;
      I_ERI_K6xy_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K6xy_S_S_S_vrr;
      I_ERI_K6xz_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K6xz_S_S_S_vrr;
      I_ERI_K5x2y_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K5x2y_S_S_S_vrr;
      I_ERI_K5xyz_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K5xyz_S_S_S_vrr;
      I_ERI_K5x2z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K5x2z_S_S_S_vrr;
      I_ERI_K4x3y_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K4x3y_S_S_S_vrr;
      I_ERI_K4x2yz_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K4x2yz_S_S_S_vrr;
      I_ERI_K4xy2z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K4xy2z_S_S_S_vrr;
      I_ERI_K4x3z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K4x3z_S_S_S_vrr;
      I_ERI_K3x4y_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K3x4y_S_S_S_vrr;
      I_ERI_K3x3yz_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K3x3yz_S_S_S_vrr;
      I_ERI_K3x2y2z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K3x2y2z_S_S_S_vrr;
      I_ERI_K3xy3z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K3xy3z_S_S_S_vrr;
      I_ERI_K3x4z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K3x4z_S_S_S_vrr;
      I_ERI_K2x5y_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K2x5y_S_S_S_vrr;
      I_ERI_K2x4yz_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K2x4yz_S_S_S_vrr;
      I_ERI_K2x3y2z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K2x3y2z_S_S_S_vrr;
      I_ERI_K2x2y3z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K2x2y3z_S_S_S_vrr;
      I_ERI_K2xy4z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K2xy4z_S_S_S_vrr;
      I_ERI_K2x5z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K2x5z_S_S_S_vrr;
      I_ERI_Kx6y_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_Kx6y_S_S_S_vrr;
      I_ERI_Kx5yz_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_Kx5yz_S_S_S_vrr;
      I_ERI_Kx4y2z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_Kx4y2z_S_S_S_vrr;
      I_ERI_Kx3y3z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_Kx3y3z_S_S_S_vrr;
      I_ERI_Kx2y4z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_Kx2y4z_S_S_S_vrr;
      I_ERI_Kxy5z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_Kxy5z_S_S_S_vrr;
      I_ERI_Kx6z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_Kx6z_S_S_S_vrr;
      I_ERI_K7y_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K7y_S_S_S_vrr;
      I_ERI_K6yz_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K6yz_S_S_S_vrr;
      I_ERI_K5y2z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K5y2z_S_S_S_vrr;
      I_ERI_K4y3z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K4y3z_S_S_S_vrr;
      I_ERI_K3y4z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K3y4z_S_S_S_vrr;
      I_ERI_K2y5z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K2y5z_S_S_S_vrr;
      I_ERI_Ky6z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_Ky6z_S_S_S_vrr;
      I_ERI_K7z_S_S_S_aa += SQ_ERI_K_S_S_S_aa_coefs*I_ERI_K7z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_I_S_S_S_aa_coefs = alpha*alpha;
      I_ERI_I6x_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I6x_S_S_S_vrr;
      I_ERI_I5xy_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I5xy_S_S_S_vrr;
      I_ERI_I5xz_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I5xz_S_S_S_vrr;
      I_ERI_I4x2y_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I4x2y_S_S_S_vrr;
      I_ERI_I4xyz_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I4xyz_S_S_S_vrr;
      I_ERI_I4x2z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I4x2z_S_S_S_vrr;
      I_ERI_I3x3y_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I3x3y_S_S_S_vrr;
      I_ERI_I3x2yz_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I3x2yz_S_S_S_vrr;
      I_ERI_I3xy2z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I3xy2z_S_S_S_vrr;
      I_ERI_I3x3z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I3x3z_S_S_S_vrr;
      I_ERI_I2x4y_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I2x4y_S_S_S_vrr;
      I_ERI_I2x3yz_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I2x3yz_S_S_S_vrr;
      I_ERI_I2x2y2z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I2x2y2z_S_S_S_vrr;
      I_ERI_I2xy3z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I2xy3z_S_S_S_vrr;
      I_ERI_I2x4z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I2x4z_S_S_S_vrr;
      I_ERI_Ix5y_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_Ix5y_S_S_S_vrr;
      I_ERI_Ix4yz_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_Ix4yz_S_S_S_vrr;
      I_ERI_Ix3y2z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_Ix3y2z_S_S_S_vrr;
      I_ERI_Ix2y3z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_Ix2y3z_S_S_S_vrr;
      I_ERI_Ixy4z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_Ixy4z_S_S_S_vrr;
      I_ERI_Ix5z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_Ix5z_S_S_S_vrr;
      I_ERI_I6y_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I6y_S_S_S_vrr;
      I_ERI_I5yz_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I5yz_S_S_S_vrr;
      I_ERI_I4y2z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I4y2z_S_S_S_vrr;
      I_ERI_I3y3z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I3y3z_S_S_S_vrr;
      I_ERI_I2y4z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I2y4z_S_S_S_vrr;
      I_ERI_Iy5z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_Iy5z_S_S_S_vrr;
      I_ERI_I6z_S_S_S_aa += SQ_ERI_I_S_S_S_aa_coefs*I_ERI_I6z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_S_S_aa_coefs = alpha*alpha;
      I_ERI_H5x_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S_aa += SQ_ERI_H_S_S_S_aa_coefs*I_ERI_H5z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_P_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_Px_S_S_S += I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S += I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S += I_ERI_Pz_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_P_S_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_I_S_P_S_ac_coefs = alpha*gamma;
      I_ERI_I6x_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I6x_S_Px_S_vrr;
      I_ERI_I5xy_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I5xy_S_Px_S_vrr;
      I_ERI_I5xz_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I5xz_S_Px_S_vrr;
      I_ERI_I4x2y_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4x2y_S_Px_S_vrr;
      I_ERI_I4xyz_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4xyz_S_Px_S_vrr;
      I_ERI_I4x2z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4x2z_S_Px_S_vrr;
      I_ERI_I3x3y_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3x3y_S_Px_S_vrr;
      I_ERI_I3x2yz_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3x2yz_S_Px_S_vrr;
      I_ERI_I3xy2z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3xy2z_S_Px_S_vrr;
      I_ERI_I3x3z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3x3z_S_Px_S_vrr;
      I_ERI_I2x4y_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x4y_S_Px_S_vrr;
      I_ERI_I2x3yz_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x3yz_S_Px_S_vrr;
      I_ERI_I2x2y2z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x2y2z_S_Px_S_vrr;
      I_ERI_I2xy3z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2xy3z_S_Px_S_vrr;
      I_ERI_I2x4z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x4z_S_Px_S_vrr;
      I_ERI_Ix5y_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix5y_S_Px_S_vrr;
      I_ERI_Ix4yz_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix4yz_S_Px_S_vrr;
      I_ERI_Ix3y2z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix3y2z_S_Px_S_vrr;
      I_ERI_Ix2y3z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix2y3z_S_Px_S_vrr;
      I_ERI_Ixy4z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ixy4z_S_Px_S_vrr;
      I_ERI_Ix5z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix5z_S_Px_S_vrr;
      I_ERI_I6y_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I6y_S_Px_S_vrr;
      I_ERI_I5yz_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I5yz_S_Px_S_vrr;
      I_ERI_I4y2z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4y2z_S_Px_S_vrr;
      I_ERI_I3y3z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3y3z_S_Px_S_vrr;
      I_ERI_I2y4z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2y4z_S_Px_S_vrr;
      I_ERI_Iy5z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Iy5z_S_Px_S_vrr;
      I_ERI_I6z_S_Px_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I6z_S_Px_S_vrr;
      I_ERI_I6x_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I6x_S_Py_S_vrr;
      I_ERI_I5xy_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I5xy_S_Py_S_vrr;
      I_ERI_I5xz_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I5xz_S_Py_S_vrr;
      I_ERI_I4x2y_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4x2y_S_Py_S_vrr;
      I_ERI_I4xyz_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4xyz_S_Py_S_vrr;
      I_ERI_I4x2z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4x2z_S_Py_S_vrr;
      I_ERI_I3x3y_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3x3y_S_Py_S_vrr;
      I_ERI_I3x2yz_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3x2yz_S_Py_S_vrr;
      I_ERI_I3xy2z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3xy2z_S_Py_S_vrr;
      I_ERI_I3x3z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3x3z_S_Py_S_vrr;
      I_ERI_I2x4y_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x4y_S_Py_S_vrr;
      I_ERI_I2x3yz_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x3yz_S_Py_S_vrr;
      I_ERI_I2x2y2z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x2y2z_S_Py_S_vrr;
      I_ERI_I2xy3z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2xy3z_S_Py_S_vrr;
      I_ERI_I2x4z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x4z_S_Py_S_vrr;
      I_ERI_Ix5y_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix5y_S_Py_S_vrr;
      I_ERI_Ix4yz_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix4yz_S_Py_S_vrr;
      I_ERI_Ix3y2z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix3y2z_S_Py_S_vrr;
      I_ERI_Ix2y3z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix2y3z_S_Py_S_vrr;
      I_ERI_Ixy4z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ixy4z_S_Py_S_vrr;
      I_ERI_Ix5z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix5z_S_Py_S_vrr;
      I_ERI_I6y_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I6y_S_Py_S_vrr;
      I_ERI_I5yz_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I5yz_S_Py_S_vrr;
      I_ERI_I4y2z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4y2z_S_Py_S_vrr;
      I_ERI_I3y3z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3y3z_S_Py_S_vrr;
      I_ERI_I2y4z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2y4z_S_Py_S_vrr;
      I_ERI_Iy5z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Iy5z_S_Py_S_vrr;
      I_ERI_I6z_S_Py_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I6z_S_Py_S_vrr;
      I_ERI_I6x_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I6x_S_Pz_S_vrr;
      I_ERI_I5xy_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I5xy_S_Pz_S_vrr;
      I_ERI_I5xz_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I5xz_S_Pz_S_vrr;
      I_ERI_I4x2y_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4x2y_S_Pz_S_vrr;
      I_ERI_I4xyz_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4xyz_S_Pz_S_vrr;
      I_ERI_I4x2z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4x2z_S_Pz_S_vrr;
      I_ERI_I3x3y_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3x3y_S_Pz_S_vrr;
      I_ERI_I3x2yz_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3x2yz_S_Pz_S_vrr;
      I_ERI_I3xy2z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3xy2z_S_Pz_S_vrr;
      I_ERI_I3x3z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3x3z_S_Pz_S_vrr;
      I_ERI_I2x4y_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x4y_S_Pz_S_vrr;
      I_ERI_I2x3yz_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x3yz_S_Pz_S_vrr;
      I_ERI_I2x2y2z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x2y2z_S_Pz_S_vrr;
      I_ERI_I2xy3z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2xy3z_S_Pz_S_vrr;
      I_ERI_I2x4z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2x4z_S_Pz_S_vrr;
      I_ERI_Ix5y_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix5y_S_Pz_S_vrr;
      I_ERI_Ix4yz_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix4yz_S_Pz_S_vrr;
      I_ERI_Ix3y2z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix3y2z_S_Pz_S_vrr;
      I_ERI_Ix2y3z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix2y3z_S_Pz_S_vrr;
      I_ERI_Ixy4z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ixy4z_S_Pz_S_vrr;
      I_ERI_Ix5z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Ix5z_S_Pz_S_vrr;
      I_ERI_I6y_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I6y_S_Pz_S_vrr;
      I_ERI_I5yz_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I5yz_S_Pz_S_vrr;
      I_ERI_I4y2z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I4y2z_S_Pz_S_vrr;
      I_ERI_I3y3z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I3y3z_S_Pz_S_vrr;
      I_ERI_I2y4z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I2y4z_S_Pz_S_vrr;
      I_ERI_Iy5z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_Iy5z_S_Pz_S_vrr;
      I_ERI_I6z_S_Pz_S_ac += SQ_ERI_I_S_P_S_ac_coefs*I_ERI_I6z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_P_S_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_P_S_ac_coefs = alpha*gamma;
      I_ERI_H5x_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H5x_S_Px_S_vrr;
      I_ERI_H4xy_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H4xy_S_Px_S_vrr;
      I_ERI_H4xz_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H4xz_S_Px_S_vrr;
      I_ERI_H3x2y_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3x2y_S_Px_S_vrr;
      I_ERI_H3xyz_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3xyz_S_Px_S_vrr;
      I_ERI_H3x2z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3x2z_S_Px_S_vrr;
      I_ERI_H2x3y_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2x3y_S_Px_S_vrr;
      I_ERI_H2x2yz_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2x2yz_S_Px_S_vrr;
      I_ERI_H2xy2z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2xy2z_S_Px_S_vrr;
      I_ERI_H2x3z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2x3z_S_Px_S_vrr;
      I_ERI_Hx4y_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx4y_S_Px_S_vrr;
      I_ERI_Hx3yz_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx3yz_S_Px_S_vrr;
      I_ERI_Hx2y2z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx2y2z_S_Px_S_vrr;
      I_ERI_Hxy3z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hxy3z_S_Px_S_vrr;
      I_ERI_Hx4z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx4z_S_Px_S_vrr;
      I_ERI_H5y_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H5y_S_Px_S_vrr;
      I_ERI_H4yz_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H4yz_S_Px_S_vrr;
      I_ERI_H3y2z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3y2z_S_Px_S_vrr;
      I_ERI_H2y3z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2y3z_S_Px_S_vrr;
      I_ERI_Hy4z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hy4z_S_Px_S_vrr;
      I_ERI_H5z_S_Px_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H5z_S_Px_S_vrr;
      I_ERI_H5x_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H5x_S_Py_S_vrr;
      I_ERI_H4xy_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H4xy_S_Py_S_vrr;
      I_ERI_H4xz_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H4xz_S_Py_S_vrr;
      I_ERI_H3x2y_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3x2y_S_Py_S_vrr;
      I_ERI_H3xyz_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3xyz_S_Py_S_vrr;
      I_ERI_H3x2z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3x2z_S_Py_S_vrr;
      I_ERI_H2x3y_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2x3y_S_Py_S_vrr;
      I_ERI_H2x2yz_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2x2yz_S_Py_S_vrr;
      I_ERI_H2xy2z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2xy2z_S_Py_S_vrr;
      I_ERI_H2x3z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2x3z_S_Py_S_vrr;
      I_ERI_Hx4y_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx4y_S_Py_S_vrr;
      I_ERI_Hx3yz_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx3yz_S_Py_S_vrr;
      I_ERI_Hx2y2z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx2y2z_S_Py_S_vrr;
      I_ERI_Hxy3z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hxy3z_S_Py_S_vrr;
      I_ERI_Hx4z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx4z_S_Py_S_vrr;
      I_ERI_H5y_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H5y_S_Py_S_vrr;
      I_ERI_H4yz_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H4yz_S_Py_S_vrr;
      I_ERI_H3y2z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3y2z_S_Py_S_vrr;
      I_ERI_H2y3z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2y3z_S_Py_S_vrr;
      I_ERI_Hy4z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hy4z_S_Py_S_vrr;
      I_ERI_H5z_S_Py_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H5z_S_Py_S_vrr;
      I_ERI_H5x_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H5x_S_Pz_S_vrr;
      I_ERI_H4xy_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H4xy_S_Pz_S_vrr;
      I_ERI_H4xz_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H4xz_S_Pz_S_vrr;
      I_ERI_H3x2y_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3x2y_S_Pz_S_vrr;
      I_ERI_H3xyz_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3xyz_S_Pz_S_vrr;
      I_ERI_H3x2z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3x2z_S_Pz_S_vrr;
      I_ERI_H2x3y_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2x3y_S_Pz_S_vrr;
      I_ERI_H2x2yz_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2x2yz_S_Pz_S_vrr;
      I_ERI_H2xy2z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2xy2z_S_Pz_S_vrr;
      I_ERI_H2x3z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2x3z_S_Pz_S_vrr;
      I_ERI_Hx4y_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx4y_S_Pz_S_vrr;
      I_ERI_Hx3yz_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx3yz_S_Pz_S_vrr;
      I_ERI_Hx2y2z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx2y2z_S_Pz_S_vrr;
      I_ERI_Hxy3z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hxy3z_S_Pz_S_vrr;
      I_ERI_Hx4z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hx4z_S_Pz_S_vrr;
      I_ERI_H5y_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H5y_S_Pz_S_vrr;
      I_ERI_H4yz_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H4yz_S_Pz_S_vrr;
      I_ERI_H3y2z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H3y2z_S_Pz_S_vrr;
      I_ERI_H2y3z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H2y3z_S_Pz_S_vrr;
      I_ERI_Hy4z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_Hy4z_S_Pz_S_vrr;
      I_ERI_H5z_S_Pz_S_ac += SQ_ERI_H_S_P_S_ac_coefs*I_ERI_H5z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_ac_coefs = alpha*gamma;
      I_ERI_G4x_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_ac += SQ_ERI_G_S_P_S_ac_coefs*I_ERI_G4z_S_Pz_S_vrr;

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
       * shell quartet name: SQ_ERI_H_S_S_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_S_S_b_coefs = beta;
      I_ERI_H5x_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S_b += SQ_ERI_H_S_S_S_b_coefs*I_ERI_H5z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_H_S_D_S_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_D_S_cc_coefs = gamma*gamma;
      I_ERI_H5x_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5x_S_D2x_S_vrr;
      I_ERI_H4xy_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xy_S_D2x_S_vrr;
      I_ERI_H4xz_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xz_S_D2x_S_vrr;
      I_ERI_H3x2y_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2y_S_D2x_S_vrr;
      I_ERI_H3xyz_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3xyz_S_D2x_S_vrr;
      I_ERI_H3x2z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2z_S_D2x_S_vrr;
      I_ERI_H2x3y_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3y_S_D2x_S_vrr;
      I_ERI_H2x2yz_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x2yz_S_D2x_S_vrr;
      I_ERI_H2xy2z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2xy2z_S_D2x_S_vrr;
      I_ERI_H2x3z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3z_S_D2x_S_vrr;
      I_ERI_Hx4y_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4y_S_D2x_S_vrr;
      I_ERI_Hx3yz_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx3yz_S_D2x_S_vrr;
      I_ERI_Hx2y2z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx2y2z_S_D2x_S_vrr;
      I_ERI_Hxy3z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hxy3z_S_D2x_S_vrr;
      I_ERI_Hx4z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4z_S_D2x_S_vrr;
      I_ERI_H5y_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5y_S_D2x_S_vrr;
      I_ERI_H4yz_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4yz_S_D2x_S_vrr;
      I_ERI_H3y2z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3y2z_S_D2x_S_vrr;
      I_ERI_H2y3z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2y3z_S_D2x_S_vrr;
      I_ERI_Hy4z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hy4z_S_D2x_S_vrr;
      I_ERI_H5z_S_D2x_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5z_S_D2x_S_vrr;
      I_ERI_H5x_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5x_S_Dxy_S_vrr;
      I_ERI_H4xy_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xy_S_Dxy_S_vrr;
      I_ERI_H4xz_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xz_S_Dxy_S_vrr;
      I_ERI_H3x2y_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2y_S_Dxy_S_vrr;
      I_ERI_H3xyz_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3xyz_S_Dxy_S_vrr;
      I_ERI_H3x2z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2z_S_Dxy_S_vrr;
      I_ERI_H2x3y_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3y_S_Dxy_S_vrr;
      I_ERI_H2x2yz_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x2yz_S_Dxy_S_vrr;
      I_ERI_H2xy2z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2xy2z_S_Dxy_S_vrr;
      I_ERI_H2x3z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3z_S_Dxy_S_vrr;
      I_ERI_Hx4y_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4y_S_Dxy_S_vrr;
      I_ERI_Hx3yz_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx3yz_S_Dxy_S_vrr;
      I_ERI_Hx2y2z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx2y2z_S_Dxy_S_vrr;
      I_ERI_Hxy3z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hxy3z_S_Dxy_S_vrr;
      I_ERI_Hx4z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4z_S_Dxy_S_vrr;
      I_ERI_H5y_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5y_S_Dxy_S_vrr;
      I_ERI_H4yz_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4yz_S_Dxy_S_vrr;
      I_ERI_H3y2z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3y2z_S_Dxy_S_vrr;
      I_ERI_H2y3z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2y3z_S_Dxy_S_vrr;
      I_ERI_Hy4z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hy4z_S_Dxy_S_vrr;
      I_ERI_H5z_S_Dxy_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5z_S_Dxy_S_vrr;
      I_ERI_H5x_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5x_S_Dxz_S_vrr;
      I_ERI_H4xy_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xy_S_Dxz_S_vrr;
      I_ERI_H4xz_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xz_S_Dxz_S_vrr;
      I_ERI_H3x2y_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2y_S_Dxz_S_vrr;
      I_ERI_H3xyz_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3xyz_S_Dxz_S_vrr;
      I_ERI_H3x2z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2z_S_Dxz_S_vrr;
      I_ERI_H2x3y_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3y_S_Dxz_S_vrr;
      I_ERI_H2x2yz_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x2yz_S_Dxz_S_vrr;
      I_ERI_H2xy2z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2xy2z_S_Dxz_S_vrr;
      I_ERI_H2x3z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3z_S_Dxz_S_vrr;
      I_ERI_Hx4y_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4y_S_Dxz_S_vrr;
      I_ERI_Hx3yz_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx3yz_S_Dxz_S_vrr;
      I_ERI_Hx2y2z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx2y2z_S_Dxz_S_vrr;
      I_ERI_Hxy3z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hxy3z_S_Dxz_S_vrr;
      I_ERI_Hx4z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4z_S_Dxz_S_vrr;
      I_ERI_H5y_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5y_S_Dxz_S_vrr;
      I_ERI_H4yz_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4yz_S_Dxz_S_vrr;
      I_ERI_H3y2z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3y2z_S_Dxz_S_vrr;
      I_ERI_H2y3z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2y3z_S_Dxz_S_vrr;
      I_ERI_Hy4z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hy4z_S_Dxz_S_vrr;
      I_ERI_H5z_S_Dxz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5z_S_Dxz_S_vrr;
      I_ERI_H5x_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5x_S_D2y_S_vrr;
      I_ERI_H4xy_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xy_S_D2y_S_vrr;
      I_ERI_H4xz_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xz_S_D2y_S_vrr;
      I_ERI_H3x2y_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2y_S_D2y_S_vrr;
      I_ERI_H3xyz_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3xyz_S_D2y_S_vrr;
      I_ERI_H3x2z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2z_S_D2y_S_vrr;
      I_ERI_H2x3y_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3y_S_D2y_S_vrr;
      I_ERI_H2x2yz_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x2yz_S_D2y_S_vrr;
      I_ERI_H2xy2z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2xy2z_S_D2y_S_vrr;
      I_ERI_H2x3z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3z_S_D2y_S_vrr;
      I_ERI_Hx4y_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4y_S_D2y_S_vrr;
      I_ERI_Hx3yz_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx3yz_S_D2y_S_vrr;
      I_ERI_Hx2y2z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx2y2z_S_D2y_S_vrr;
      I_ERI_Hxy3z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hxy3z_S_D2y_S_vrr;
      I_ERI_Hx4z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4z_S_D2y_S_vrr;
      I_ERI_H5y_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5y_S_D2y_S_vrr;
      I_ERI_H4yz_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4yz_S_D2y_S_vrr;
      I_ERI_H3y2z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3y2z_S_D2y_S_vrr;
      I_ERI_H2y3z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2y3z_S_D2y_S_vrr;
      I_ERI_Hy4z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hy4z_S_D2y_S_vrr;
      I_ERI_H5z_S_D2y_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5z_S_D2y_S_vrr;
      I_ERI_H5x_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5x_S_Dyz_S_vrr;
      I_ERI_H4xy_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xy_S_Dyz_S_vrr;
      I_ERI_H4xz_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xz_S_Dyz_S_vrr;
      I_ERI_H3x2y_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2y_S_Dyz_S_vrr;
      I_ERI_H3xyz_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3xyz_S_Dyz_S_vrr;
      I_ERI_H3x2z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2z_S_Dyz_S_vrr;
      I_ERI_H2x3y_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3y_S_Dyz_S_vrr;
      I_ERI_H2x2yz_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x2yz_S_Dyz_S_vrr;
      I_ERI_H2xy2z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2xy2z_S_Dyz_S_vrr;
      I_ERI_H2x3z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3z_S_Dyz_S_vrr;
      I_ERI_Hx4y_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4y_S_Dyz_S_vrr;
      I_ERI_Hx3yz_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx3yz_S_Dyz_S_vrr;
      I_ERI_Hx2y2z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx2y2z_S_Dyz_S_vrr;
      I_ERI_Hxy3z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hxy3z_S_Dyz_S_vrr;
      I_ERI_Hx4z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4z_S_Dyz_S_vrr;
      I_ERI_H5y_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5y_S_Dyz_S_vrr;
      I_ERI_H4yz_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4yz_S_Dyz_S_vrr;
      I_ERI_H3y2z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3y2z_S_Dyz_S_vrr;
      I_ERI_H2y3z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2y3z_S_Dyz_S_vrr;
      I_ERI_Hy4z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hy4z_S_Dyz_S_vrr;
      I_ERI_H5z_S_Dyz_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5z_S_Dyz_S_vrr;
      I_ERI_H5x_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5x_S_D2z_S_vrr;
      I_ERI_H4xy_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xy_S_D2z_S_vrr;
      I_ERI_H4xz_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4xz_S_D2z_S_vrr;
      I_ERI_H3x2y_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2y_S_D2z_S_vrr;
      I_ERI_H3xyz_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3xyz_S_D2z_S_vrr;
      I_ERI_H3x2z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3x2z_S_D2z_S_vrr;
      I_ERI_H2x3y_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3y_S_D2z_S_vrr;
      I_ERI_H2x2yz_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x2yz_S_D2z_S_vrr;
      I_ERI_H2xy2z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2xy2z_S_D2z_S_vrr;
      I_ERI_H2x3z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2x3z_S_D2z_S_vrr;
      I_ERI_Hx4y_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4y_S_D2z_S_vrr;
      I_ERI_Hx3yz_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx3yz_S_D2z_S_vrr;
      I_ERI_Hx2y2z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx2y2z_S_D2z_S_vrr;
      I_ERI_Hxy3z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hxy3z_S_D2z_S_vrr;
      I_ERI_Hx4z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hx4z_S_D2z_S_vrr;
      I_ERI_H5y_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5y_S_D2z_S_vrr;
      I_ERI_H4yz_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H4yz_S_D2z_S_vrr;
      I_ERI_H3y2z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H3y2z_S_D2z_S_vrr;
      I_ERI_H2y3z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H2y3z_S_D2z_S_vrr;
      I_ERI_Hy4z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_Hy4z_S_D2z_S_vrr;
      I_ERI_H5z_S_D2z_S_cc += SQ_ERI_H_S_D_S_cc_coefs*I_ERI_H5z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_cc_coefs = gamma*gamma;
      I_ERI_G4x_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_cc += SQ_ERI_G_S_D_S_cc_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_cc_coefs = gamma*gamma;
      I_ERI_F3x_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_cc += SQ_ERI_F_S_D_S_cc_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_S_S_c_coefs = gamma;
      I_ERI_H5x_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S_c += SQ_ERI_H_S_S_S_c_coefs*I_ERI_H5z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_c_coefs = gamma;
      I_ERI_G4x_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_c += SQ_ERI_G_S_S_S_c_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_c_coefs = gamma;
      I_ERI_F3x_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_c += SQ_ERI_F_S_S_S_c_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_S_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_K_S_S_S_ab_coefs = alpha*beta;
      I_ERI_K7x_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K7x_S_S_S_vrr;
      I_ERI_K6xy_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K6xy_S_S_S_vrr;
      I_ERI_K6xz_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K6xz_S_S_S_vrr;
      I_ERI_K5x2y_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K5x2y_S_S_S_vrr;
      I_ERI_K5xyz_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K5xyz_S_S_S_vrr;
      I_ERI_K5x2z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K5x2z_S_S_S_vrr;
      I_ERI_K4x3y_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K4x3y_S_S_S_vrr;
      I_ERI_K4x2yz_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K4x2yz_S_S_S_vrr;
      I_ERI_K4xy2z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K4xy2z_S_S_S_vrr;
      I_ERI_K4x3z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K4x3z_S_S_S_vrr;
      I_ERI_K3x4y_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K3x4y_S_S_S_vrr;
      I_ERI_K3x3yz_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K3x3yz_S_S_S_vrr;
      I_ERI_K3x2y2z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K3x2y2z_S_S_S_vrr;
      I_ERI_K3xy3z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K3xy3z_S_S_S_vrr;
      I_ERI_K3x4z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K3x4z_S_S_S_vrr;
      I_ERI_K2x5y_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K2x5y_S_S_S_vrr;
      I_ERI_K2x4yz_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K2x4yz_S_S_S_vrr;
      I_ERI_K2x3y2z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K2x3y2z_S_S_S_vrr;
      I_ERI_K2x2y3z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K2x2y3z_S_S_S_vrr;
      I_ERI_K2xy4z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K2xy4z_S_S_S_vrr;
      I_ERI_K2x5z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K2x5z_S_S_S_vrr;
      I_ERI_Kx6y_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_Kx6y_S_S_S_vrr;
      I_ERI_Kx5yz_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_Kx5yz_S_S_S_vrr;
      I_ERI_Kx4y2z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_Kx4y2z_S_S_S_vrr;
      I_ERI_Kx3y3z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_Kx3y3z_S_S_S_vrr;
      I_ERI_Kx2y4z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_Kx2y4z_S_S_S_vrr;
      I_ERI_Kxy5z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_Kxy5z_S_S_S_vrr;
      I_ERI_Kx6z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_Kx6z_S_S_S_vrr;
      I_ERI_K7y_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K7y_S_S_S_vrr;
      I_ERI_K6yz_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K6yz_S_S_S_vrr;
      I_ERI_K5y2z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K5y2z_S_S_S_vrr;
      I_ERI_K4y3z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K4y3z_S_S_S_vrr;
      I_ERI_K3y4z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K3y4z_S_S_S_vrr;
      I_ERI_K2y5z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K2y5z_S_S_S_vrr;
      I_ERI_Ky6z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_Ky6z_S_S_S_vrr;
      I_ERI_K7z_S_S_S_ab += SQ_ERI_K_S_S_S_ab_coefs*I_ERI_K7z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_I_S_S_S_ab_coefs = alpha*beta;
      I_ERI_I6x_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I6x_S_S_S_vrr;
      I_ERI_I5xy_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I5xy_S_S_S_vrr;
      I_ERI_I5xz_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I5xz_S_S_S_vrr;
      I_ERI_I4x2y_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I4x2y_S_S_S_vrr;
      I_ERI_I4xyz_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I4xyz_S_S_S_vrr;
      I_ERI_I4x2z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I4x2z_S_S_S_vrr;
      I_ERI_I3x3y_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I3x3y_S_S_S_vrr;
      I_ERI_I3x2yz_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I3x2yz_S_S_S_vrr;
      I_ERI_I3xy2z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I3xy2z_S_S_S_vrr;
      I_ERI_I3x3z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I3x3z_S_S_S_vrr;
      I_ERI_I2x4y_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I2x4y_S_S_S_vrr;
      I_ERI_I2x3yz_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I2x3yz_S_S_S_vrr;
      I_ERI_I2x2y2z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I2x2y2z_S_S_S_vrr;
      I_ERI_I2xy3z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I2xy3z_S_S_S_vrr;
      I_ERI_I2x4z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I2x4z_S_S_S_vrr;
      I_ERI_Ix5y_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_Ix5y_S_S_S_vrr;
      I_ERI_Ix4yz_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_Ix4yz_S_S_S_vrr;
      I_ERI_Ix3y2z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_Ix3y2z_S_S_S_vrr;
      I_ERI_Ix2y3z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_Ix2y3z_S_S_S_vrr;
      I_ERI_Ixy4z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_Ixy4z_S_S_S_vrr;
      I_ERI_Ix5z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_Ix5z_S_S_S_vrr;
      I_ERI_I6y_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I6y_S_S_S_vrr;
      I_ERI_I5yz_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I5yz_S_S_S_vrr;
      I_ERI_I4y2z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I4y2z_S_S_S_vrr;
      I_ERI_I3y3z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I3y3z_S_S_S_vrr;
      I_ERI_I2y4z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I2y4z_S_S_S_vrr;
      I_ERI_Iy5z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_Iy5z_S_S_S_vrr;
      I_ERI_I6z_S_S_S_ab += SQ_ERI_I_S_S_S_ab_coefs*I_ERI_I6z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_S_S_ab_coefs = alpha*beta;
      I_ERI_H5x_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S_ab += SQ_ERI_H_S_S_S_ab_coefs*I_ERI_H5z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_ab_coefs = alpha*beta;
      I_ERI_G4x_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_ab += SQ_ERI_G_S_S_S_ab_coefs*I_ERI_G4z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_I_S_P_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_I_S_P_S_bc_coefs = beta*gamma;
      I_ERI_I6x_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I6x_S_Px_S_vrr;
      I_ERI_I5xy_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I5xy_S_Px_S_vrr;
      I_ERI_I5xz_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I5xz_S_Px_S_vrr;
      I_ERI_I4x2y_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4x2y_S_Px_S_vrr;
      I_ERI_I4xyz_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4xyz_S_Px_S_vrr;
      I_ERI_I4x2z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4x2z_S_Px_S_vrr;
      I_ERI_I3x3y_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3x3y_S_Px_S_vrr;
      I_ERI_I3x2yz_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3x2yz_S_Px_S_vrr;
      I_ERI_I3xy2z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3xy2z_S_Px_S_vrr;
      I_ERI_I3x3z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3x3z_S_Px_S_vrr;
      I_ERI_I2x4y_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x4y_S_Px_S_vrr;
      I_ERI_I2x3yz_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x3yz_S_Px_S_vrr;
      I_ERI_I2x2y2z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x2y2z_S_Px_S_vrr;
      I_ERI_I2xy3z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2xy3z_S_Px_S_vrr;
      I_ERI_I2x4z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x4z_S_Px_S_vrr;
      I_ERI_Ix5y_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix5y_S_Px_S_vrr;
      I_ERI_Ix4yz_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix4yz_S_Px_S_vrr;
      I_ERI_Ix3y2z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix3y2z_S_Px_S_vrr;
      I_ERI_Ix2y3z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix2y3z_S_Px_S_vrr;
      I_ERI_Ixy4z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ixy4z_S_Px_S_vrr;
      I_ERI_Ix5z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix5z_S_Px_S_vrr;
      I_ERI_I6y_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I6y_S_Px_S_vrr;
      I_ERI_I5yz_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I5yz_S_Px_S_vrr;
      I_ERI_I4y2z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4y2z_S_Px_S_vrr;
      I_ERI_I3y3z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3y3z_S_Px_S_vrr;
      I_ERI_I2y4z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2y4z_S_Px_S_vrr;
      I_ERI_Iy5z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Iy5z_S_Px_S_vrr;
      I_ERI_I6z_S_Px_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I6z_S_Px_S_vrr;
      I_ERI_I6x_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I6x_S_Py_S_vrr;
      I_ERI_I5xy_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I5xy_S_Py_S_vrr;
      I_ERI_I5xz_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I5xz_S_Py_S_vrr;
      I_ERI_I4x2y_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4x2y_S_Py_S_vrr;
      I_ERI_I4xyz_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4xyz_S_Py_S_vrr;
      I_ERI_I4x2z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4x2z_S_Py_S_vrr;
      I_ERI_I3x3y_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3x3y_S_Py_S_vrr;
      I_ERI_I3x2yz_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3x2yz_S_Py_S_vrr;
      I_ERI_I3xy2z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3xy2z_S_Py_S_vrr;
      I_ERI_I3x3z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3x3z_S_Py_S_vrr;
      I_ERI_I2x4y_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x4y_S_Py_S_vrr;
      I_ERI_I2x3yz_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x3yz_S_Py_S_vrr;
      I_ERI_I2x2y2z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x2y2z_S_Py_S_vrr;
      I_ERI_I2xy3z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2xy3z_S_Py_S_vrr;
      I_ERI_I2x4z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x4z_S_Py_S_vrr;
      I_ERI_Ix5y_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix5y_S_Py_S_vrr;
      I_ERI_Ix4yz_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix4yz_S_Py_S_vrr;
      I_ERI_Ix3y2z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix3y2z_S_Py_S_vrr;
      I_ERI_Ix2y3z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix2y3z_S_Py_S_vrr;
      I_ERI_Ixy4z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ixy4z_S_Py_S_vrr;
      I_ERI_Ix5z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix5z_S_Py_S_vrr;
      I_ERI_I6y_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I6y_S_Py_S_vrr;
      I_ERI_I5yz_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I5yz_S_Py_S_vrr;
      I_ERI_I4y2z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4y2z_S_Py_S_vrr;
      I_ERI_I3y3z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3y3z_S_Py_S_vrr;
      I_ERI_I2y4z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2y4z_S_Py_S_vrr;
      I_ERI_Iy5z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Iy5z_S_Py_S_vrr;
      I_ERI_I6z_S_Py_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I6z_S_Py_S_vrr;
      I_ERI_I6x_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I6x_S_Pz_S_vrr;
      I_ERI_I5xy_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I5xy_S_Pz_S_vrr;
      I_ERI_I5xz_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I5xz_S_Pz_S_vrr;
      I_ERI_I4x2y_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4x2y_S_Pz_S_vrr;
      I_ERI_I4xyz_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4xyz_S_Pz_S_vrr;
      I_ERI_I4x2z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4x2z_S_Pz_S_vrr;
      I_ERI_I3x3y_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3x3y_S_Pz_S_vrr;
      I_ERI_I3x2yz_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3x2yz_S_Pz_S_vrr;
      I_ERI_I3xy2z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3xy2z_S_Pz_S_vrr;
      I_ERI_I3x3z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3x3z_S_Pz_S_vrr;
      I_ERI_I2x4y_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x4y_S_Pz_S_vrr;
      I_ERI_I2x3yz_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x3yz_S_Pz_S_vrr;
      I_ERI_I2x2y2z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x2y2z_S_Pz_S_vrr;
      I_ERI_I2xy3z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2xy3z_S_Pz_S_vrr;
      I_ERI_I2x4z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2x4z_S_Pz_S_vrr;
      I_ERI_Ix5y_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix5y_S_Pz_S_vrr;
      I_ERI_Ix4yz_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix4yz_S_Pz_S_vrr;
      I_ERI_Ix3y2z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix3y2z_S_Pz_S_vrr;
      I_ERI_Ix2y3z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix2y3z_S_Pz_S_vrr;
      I_ERI_Ixy4z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ixy4z_S_Pz_S_vrr;
      I_ERI_Ix5z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Ix5z_S_Pz_S_vrr;
      I_ERI_I6y_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I6y_S_Pz_S_vrr;
      I_ERI_I5yz_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I5yz_S_Pz_S_vrr;
      I_ERI_I4y2z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I4y2z_S_Pz_S_vrr;
      I_ERI_I3y3z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I3y3z_S_Pz_S_vrr;
      I_ERI_I2y4z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I2y4z_S_Pz_S_vrr;
      I_ERI_Iy5z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_Iy5z_S_Pz_S_vrr;
      I_ERI_I6z_S_Pz_S_bc += SQ_ERI_I_S_P_S_bc_coefs*I_ERI_I6z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_P_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_P_S_bc_coefs = beta*gamma;
      I_ERI_H5x_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H5x_S_Px_S_vrr;
      I_ERI_H4xy_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H4xy_S_Px_S_vrr;
      I_ERI_H4xz_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H4xz_S_Px_S_vrr;
      I_ERI_H3x2y_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3x2y_S_Px_S_vrr;
      I_ERI_H3xyz_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3xyz_S_Px_S_vrr;
      I_ERI_H3x2z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3x2z_S_Px_S_vrr;
      I_ERI_H2x3y_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2x3y_S_Px_S_vrr;
      I_ERI_H2x2yz_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2x2yz_S_Px_S_vrr;
      I_ERI_H2xy2z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2xy2z_S_Px_S_vrr;
      I_ERI_H2x3z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2x3z_S_Px_S_vrr;
      I_ERI_Hx4y_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx4y_S_Px_S_vrr;
      I_ERI_Hx3yz_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx3yz_S_Px_S_vrr;
      I_ERI_Hx2y2z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx2y2z_S_Px_S_vrr;
      I_ERI_Hxy3z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hxy3z_S_Px_S_vrr;
      I_ERI_Hx4z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx4z_S_Px_S_vrr;
      I_ERI_H5y_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H5y_S_Px_S_vrr;
      I_ERI_H4yz_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H4yz_S_Px_S_vrr;
      I_ERI_H3y2z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3y2z_S_Px_S_vrr;
      I_ERI_H2y3z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2y3z_S_Px_S_vrr;
      I_ERI_Hy4z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hy4z_S_Px_S_vrr;
      I_ERI_H5z_S_Px_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H5z_S_Px_S_vrr;
      I_ERI_H5x_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H5x_S_Py_S_vrr;
      I_ERI_H4xy_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H4xy_S_Py_S_vrr;
      I_ERI_H4xz_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H4xz_S_Py_S_vrr;
      I_ERI_H3x2y_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3x2y_S_Py_S_vrr;
      I_ERI_H3xyz_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3xyz_S_Py_S_vrr;
      I_ERI_H3x2z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3x2z_S_Py_S_vrr;
      I_ERI_H2x3y_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2x3y_S_Py_S_vrr;
      I_ERI_H2x2yz_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2x2yz_S_Py_S_vrr;
      I_ERI_H2xy2z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2xy2z_S_Py_S_vrr;
      I_ERI_H2x3z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2x3z_S_Py_S_vrr;
      I_ERI_Hx4y_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx4y_S_Py_S_vrr;
      I_ERI_Hx3yz_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx3yz_S_Py_S_vrr;
      I_ERI_Hx2y2z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx2y2z_S_Py_S_vrr;
      I_ERI_Hxy3z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hxy3z_S_Py_S_vrr;
      I_ERI_Hx4z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx4z_S_Py_S_vrr;
      I_ERI_H5y_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H5y_S_Py_S_vrr;
      I_ERI_H4yz_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H4yz_S_Py_S_vrr;
      I_ERI_H3y2z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3y2z_S_Py_S_vrr;
      I_ERI_H2y3z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2y3z_S_Py_S_vrr;
      I_ERI_Hy4z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hy4z_S_Py_S_vrr;
      I_ERI_H5z_S_Py_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H5z_S_Py_S_vrr;
      I_ERI_H5x_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H5x_S_Pz_S_vrr;
      I_ERI_H4xy_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H4xy_S_Pz_S_vrr;
      I_ERI_H4xz_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H4xz_S_Pz_S_vrr;
      I_ERI_H3x2y_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3x2y_S_Pz_S_vrr;
      I_ERI_H3xyz_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3xyz_S_Pz_S_vrr;
      I_ERI_H3x2z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3x2z_S_Pz_S_vrr;
      I_ERI_H2x3y_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2x3y_S_Pz_S_vrr;
      I_ERI_H2x2yz_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2x2yz_S_Pz_S_vrr;
      I_ERI_H2xy2z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2xy2z_S_Pz_S_vrr;
      I_ERI_H2x3z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2x3z_S_Pz_S_vrr;
      I_ERI_Hx4y_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx4y_S_Pz_S_vrr;
      I_ERI_Hx3yz_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx3yz_S_Pz_S_vrr;
      I_ERI_Hx2y2z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx2y2z_S_Pz_S_vrr;
      I_ERI_Hxy3z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hxy3z_S_Pz_S_vrr;
      I_ERI_Hx4z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hx4z_S_Pz_S_vrr;
      I_ERI_H5y_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H5y_S_Pz_S_vrr;
      I_ERI_H4yz_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H4yz_S_Pz_S_vrr;
      I_ERI_H3y2z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H3y2z_S_Pz_S_vrr;
      I_ERI_H2y3z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H2y3z_S_Pz_S_vrr;
      I_ERI_Hy4z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_Hy4z_S_Pz_S_vrr;
      I_ERI_H5z_S_Pz_S_bc += SQ_ERI_H_S_P_S_bc_coefs*I_ERI_H5z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_P_S_bc_coefs = beta*gamma;
      I_ERI_G4x_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S_bc += SQ_ERI_G_S_P_S_bc_coefs*I_ERI_G4z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_bc_coefs = beta*gamma;
      I_ERI_F3x_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_bc += SQ_ERI_F_S_P_S_bc_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_S_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_K_S_S_S_bb_coefs = beta*beta;
      I_ERI_K7x_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K7x_S_S_S_vrr;
      I_ERI_K6xy_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K6xy_S_S_S_vrr;
      I_ERI_K6xz_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K6xz_S_S_S_vrr;
      I_ERI_K5x2y_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K5x2y_S_S_S_vrr;
      I_ERI_K5xyz_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K5xyz_S_S_S_vrr;
      I_ERI_K5x2z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K5x2z_S_S_S_vrr;
      I_ERI_K4x3y_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K4x3y_S_S_S_vrr;
      I_ERI_K4x2yz_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K4x2yz_S_S_S_vrr;
      I_ERI_K4xy2z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K4xy2z_S_S_S_vrr;
      I_ERI_K4x3z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K4x3z_S_S_S_vrr;
      I_ERI_K3x4y_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K3x4y_S_S_S_vrr;
      I_ERI_K3x3yz_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K3x3yz_S_S_S_vrr;
      I_ERI_K3x2y2z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K3x2y2z_S_S_S_vrr;
      I_ERI_K3xy3z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K3xy3z_S_S_S_vrr;
      I_ERI_K3x4z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K3x4z_S_S_S_vrr;
      I_ERI_K2x5y_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K2x5y_S_S_S_vrr;
      I_ERI_K2x4yz_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K2x4yz_S_S_S_vrr;
      I_ERI_K2x3y2z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K2x3y2z_S_S_S_vrr;
      I_ERI_K2x2y3z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K2x2y3z_S_S_S_vrr;
      I_ERI_K2xy4z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K2xy4z_S_S_S_vrr;
      I_ERI_K2x5z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K2x5z_S_S_S_vrr;
      I_ERI_Kx6y_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_Kx6y_S_S_S_vrr;
      I_ERI_Kx5yz_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_Kx5yz_S_S_S_vrr;
      I_ERI_Kx4y2z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_Kx4y2z_S_S_S_vrr;
      I_ERI_Kx3y3z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_Kx3y3z_S_S_S_vrr;
      I_ERI_Kx2y4z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_Kx2y4z_S_S_S_vrr;
      I_ERI_Kxy5z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_Kxy5z_S_S_S_vrr;
      I_ERI_Kx6z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_Kx6z_S_S_S_vrr;
      I_ERI_K7y_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K7y_S_S_S_vrr;
      I_ERI_K6yz_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K6yz_S_S_S_vrr;
      I_ERI_K5y2z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K5y2z_S_S_S_vrr;
      I_ERI_K4y3z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K4y3z_S_S_S_vrr;
      I_ERI_K3y4z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K3y4z_S_S_S_vrr;
      I_ERI_K2y5z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K2y5z_S_S_S_vrr;
      I_ERI_Ky6z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_Ky6z_S_S_S_vrr;
      I_ERI_K7z_S_S_S_bb += SQ_ERI_K_S_S_S_bb_coefs*I_ERI_K7z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_I_S_S_S_bb_coefs = beta*beta;
      I_ERI_I6x_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I6x_S_S_S_vrr;
      I_ERI_I5xy_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I5xy_S_S_S_vrr;
      I_ERI_I5xz_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I5xz_S_S_S_vrr;
      I_ERI_I4x2y_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I4x2y_S_S_S_vrr;
      I_ERI_I4xyz_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I4xyz_S_S_S_vrr;
      I_ERI_I4x2z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I4x2z_S_S_S_vrr;
      I_ERI_I3x3y_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I3x3y_S_S_S_vrr;
      I_ERI_I3x2yz_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I3x2yz_S_S_S_vrr;
      I_ERI_I3xy2z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I3xy2z_S_S_S_vrr;
      I_ERI_I3x3z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I3x3z_S_S_S_vrr;
      I_ERI_I2x4y_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I2x4y_S_S_S_vrr;
      I_ERI_I2x3yz_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I2x3yz_S_S_S_vrr;
      I_ERI_I2x2y2z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I2x2y2z_S_S_S_vrr;
      I_ERI_I2xy3z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I2xy3z_S_S_S_vrr;
      I_ERI_I2x4z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I2x4z_S_S_S_vrr;
      I_ERI_Ix5y_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_Ix5y_S_S_S_vrr;
      I_ERI_Ix4yz_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_Ix4yz_S_S_S_vrr;
      I_ERI_Ix3y2z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_Ix3y2z_S_S_S_vrr;
      I_ERI_Ix2y3z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_Ix2y3z_S_S_S_vrr;
      I_ERI_Ixy4z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_Ixy4z_S_S_S_vrr;
      I_ERI_Ix5z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_Ix5z_S_S_S_vrr;
      I_ERI_I6y_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I6y_S_S_S_vrr;
      I_ERI_I5yz_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I5yz_S_S_S_vrr;
      I_ERI_I4y2z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I4y2z_S_S_S_vrr;
      I_ERI_I3y3z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I3y3z_S_S_S_vrr;
      I_ERI_I2y4z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I2y4z_S_S_S_vrr;
      I_ERI_Iy5z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_Iy5z_S_S_S_vrr;
      I_ERI_I6z_S_S_S_bb += SQ_ERI_I_S_S_S_bb_coefs*I_ERI_I6z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_S_S_bb_coefs = beta*beta;
      I_ERI_H5x_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H5x_S_S_S_vrr;
      I_ERI_H4xy_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H4xy_S_S_S_vrr;
      I_ERI_H4xz_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H4xz_S_S_S_vrr;
      I_ERI_H3x2y_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H3x2y_S_S_S_vrr;
      I_ERI_H3xyz_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H3xyz_S_S_S_vrr;
      I_ERI_H3x2z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H3x2z_S_S_S_vrr;
      I_ERI_H2x3y_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H2x3y_S_S_S_vrr;
      I_ERI_H2x2yz_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H2x2yz_S_S_S_vrr;
      I_ERI_H2xy2z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H2xy2z_S_S_S_vrr;
      I_ERI_H2x3z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H2x3z_S_S_S_vrr;
      I_ERI_Hx4y_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_Hx4y_S_S_S_vrr;
      I_ERI_Hx3yz_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_Hx3yz_S_S_S_vrr;
      I_ERI_Hx2y2z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_Hx2y2z_S_S_S_vrr;
      I_ERI_Hxy3z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_Hxy3z_S_S_S_vrr;
      I_ERI_Hx4z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_Hx4z_S_S_S_vrr;
      I_ERI_H5y_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H5y_S_S_S_vrr;
      I_ERI_H4yz_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H4yz_S_S_S_vrr;
      I_ERI_H3y2z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H3y2z_S_S_S_vrr;
      I_ERI_H2y3z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H2y3z_S_S_S_vrr;
      I_ERI_Hy4z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_Hy4z_S_S_S_vrr;
      I_ERI_H5z_S_S_S_bb += SQ_ERI_H_S_S_S_bb_coefs*I_ERI_H5z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_bb_coefs = beta*beta;
      I_ERI_G4x_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_bb += SQ_ERI_G_S_S_S_bb_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_bb_coefs = beta*beta;
      I_ERI_F3x_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_bb += SQ_ERI_F_S_S_S_bb_coefs*I_ERI_F3z_S_S_S_vrr;
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
   * shell quartet name: SQ_ERI_D_P_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   * RHS shell quartet name: SQ_ERI_D_S_S_S
   ************************************************************/
  Double I_ERI_D2x_Px_S_S = I_ERI_F3x_S_S_S+ABX*I_ERI_D2x_S_S_S;
  Double I_ERI_Dxy_Px_S_S = I_ERI_F2xy_S_S_S+ABX*I_ERI_Dxy_S_S_S;
  Double I_ERI_Dxz_Px_S_S = I_ERI_F2xz_S_S_S+ABX*I_ERI_Dxz_S_S_S;
  Double I_ERI_D2y_Px_S_S = I_ERI_Fx2y_S_S_S+ABX*I_ERI_D2y_S_S_S;
  Double I_ERI_Dyz_Px_S_S = I_ERI_Fxyz_S_S_S+ABX*I_ERI_Dyz_S_S_S;
  Double I_ERI_D2z_Px_S_S = I_ERI_Fx2z_S_S_S+ABX*I_ERI_D2z_S_S_S;
  Double I_ERI_D2x_Py_S_S = I_ERI_F2xy_S_S_S+ABY*I_ERI_D2x_S_S_S;
  Double I_ERI_Dxy_Py_S_S = I_ERI_Fx2y_S_S_S+ABY*I_ERI_Dxy_S_S_S;
  Double I_ERI_Dxz_Py_S_S = I_ERI_Fxyz_S_S_S+ABY*I_ERI_Dxz_S_S_S;
  Double I_ERI_D2y_Py_S_S = I_ERI_F3y_S_S_S+ABY*I_ERI_D2y_S_S_S;
  Double I_ERI_Dyz_Py_S_S = I_ERI_F2yz_S_S_S+ABY*I_ERI_Dyz_S_S_S;
  Double I_ERI_D2z_Py_S_S = I_ERI_Fy2z_S_S_S+ABY*I_ERI_D2z_S_S_S;
  Double I_ERI_D2x_Pz_S_S = I_ERI_F2xz_S_S_S+ABZ*I_ERI_D2x_S_S_S;
  Double I_ERI_Dxy_Pz_S_S = I_ERI_Fxyz_S_S_S+ABZ*I_ERI_Dxy_S_S_S;
  Double I_ERI_Dxz_Pz_S_S = I_ERI_Fx2z_S_S_S+ABZ*I_ERI_Dxz_S_S_S;
  Double I_ERI_D2y_Pz_S_S = I_ERI_F2yz_S_S_S+ABZ*I_ERI_D2y_S_S_S;
  Double I_ERI_Dyz_Pz_S_S = I_ERI_Fy2z_S_S_S+ABZ*I_ERI_Dyz_S_S_S;
  Double I_ERI_D2z_Pz_S_S = I_ERI_F3z_S_S_S+ABZ*I_ERI_D2z_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_P_D_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   * RHS shell quartet name: SQ_ERI_P_P_S_S
   ************************************************************/
  Double I_ERI_Px_D2x_S_S = I_ERI_D2x_Px_S_S+ABX*I_ERI_Px_Px_S_S;
  Double I_ERI_Py_D2x_S_S = I_ERI_Dxy_Px_S_S+ABX*I_ERI_Py_Px_S_S;
  Double I_ERI_Pz_D2x_S_S = I_ERI_Dxz_Px_S_S+ABX*I_ERI_Pz_Px_S_S;
  Double I_ERI_Px_Dxy_S_S = I_ERI_Dxy_Px_S_S+ABY*I_ERI_Px_Px_S_S;
  Double I_ERI_Py_Dxy_S_S = I_ERI_D2y_Px_S_S+ABY*I_ERI_Py_Px_S_S;
  Double I_ERI_Pz_Dxy_S_S = I_ERI_Dyz_Px_S_S+ABY*I_ERI_Pz_Px_S_S;
  Double I_ERI_Px_Dxz_S_S = I_ERI_Dxz_Px_S_S+ABZ*I_ERI_Px_Px_S_S;
  Double I_ERI_Py_Dxz_S_S = I_ERI_Dyz_Px_S_S+ABZ*I_ERI_Py_Px_S_S;
  Double I_ERI_Pz_Dxz_S_S = I_ERI_D2z_Px_S_S+ABZ*I_ERI_Pz_Px_S_S;
  Double I_ERI_Px_D2y_S_S = I_ERI_Dxy_Py_S_S+ABY*I_ERI_Px_Py_S_S;
  Double I_ERI_Py_D2y_S_S = I_ERI_D2y_Py_S_S+ABY*I_ERI_Py_Py_S_S;
  Double I_ERI_Pz_D2y_S_S = I_ERI_Dyz_Py_S_S+ABY*I_ERI_Pz_Py_S_S;
  Double I_ERI_Px_Dyz_S_S = I_ERI_Dxz_Py_S_S+ABZ*I_ERI_Px_Py_S_S;
  Double I_ERI_Py_Dyz_S_S = I_ERI_Dyz_Py_S_S+ABZ*I_ERI_Py_Py_S_S;
  Double I_ERI_Pz_Dyz_S_S = I_ERI_D2z_Py_S_S+ABZ*I_ERI_Pz_Py_S_S;
  Double I_ERI_Px_D2z_S_S = I_ERI_Dxz_Pz_S_S+ABZ*I_ERI_Px_Pz_S_S;
  Double I_ERI_Py_D2z_S_S = I_ERI_Dyz_Pz_S_S+ABZ*I_ERI_Py_Pz_S_S;
  Double I_ERI_Pz_D2z_S_S = I_ERI_D2z_Pz_S_S+ABZ*I_ERI_Pz_Pz_S_S;

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
   * shell quartet name: SQ_ERI_G_P_S_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_S_S_a
   * RHS shell quartet name: SQ_ERI_G_S_S_S_a
   ************************************************************/
  Double I_ERI_G4x_Px_S_S_a = I_ERI_H5x_S_S_S_a+ABX*I_ERI_G4x_S_S_S_a;
  Double I_ERI_G3xy_Px_S_S_a = I_ERI_H4xy_S_S_S_a+ABX*I_ERI_G3xy_S_S_S_a;
  Double I_ERI_G3xz_Px_S_S_a = I_ERI_H4xz_S_S_S_a+ABX*I_ERI_G3xz_S_S_S_a;
  Double I_ERI_G2x2y_Px_S_S_a = I_ERI_H3x2y_S_S_S_a+ABX*I_ERI_G2x2y_S_S_S_a;
  Double I_ERI_G2xyz_Px_S_S_a = I_ERI_H3xyz_S_S_S_a+ABX*I_ERI_G2xyz_S_S_S_a;
  Double I_ERI_G2x2z_Px_S_S_a = I_ERI_H3x2z_S_S_S_a+ABX*I_ERI_G2x2z_S_S_S_a;
  Double I_ERI_Gx3y_Px_S_S_a = I_ERI_H2x3y_S_S_S_a+ABX*I_ERI_Gx3y_S_S_S_a;
  Double I_ERI_Gx2yz_Px_S_S_a = I_ERI_H2x2yz_S_S_S_a+ABX*I_ERI_Gx2yz_S_S_S_a;
  Double I_ERI_Gxy2z_Px_S_S_a = I_ERI_H2xy2z_S_S_S_a+ABX*I_ERI_Gxy2z_S_S_S_a;
  Double I_ERI_Gx3z_Px_S_S_a = I_ERI_H2x3z_S_S_S_a+ABX*I_ERI_Gx3z_S_S_S_a;
  Double I_ERI_G4y_Px_S_S_a = I_ERI_Hx4y_S_S_S_a+ABX*I_ERI_G4y_S_S_S_a;
  Double I_ERI_G3yz_Px_S_S_a = I_ERI_Hx3yz_S_S_S_a+ABX*I_ERI_G3yz_S_S_S_a;
  Double I_ERI_G2y2z_Px_S_S_a = I_ERI_Hx2y2z_S_S_S_a+ABX*I_ERI_G2y2z_S_S_S_a;
  Double I_ERI_Gy3z_Px_S_S_a = I_ERI_Hxy3z_S_S_S_a+ABX*I_ERI_Gy3z_S_S_S_a;
  Double I_ERI_G4z_Px_S_S_a = I_ERI_Hx4z_S_S_S_a+ABX*I_ERI_G4z_S_S_S_a;
  Double I_ERI_G4x_Py_S_S_a = I_ERI_H4xy_S_S_S_a+ABY*I_ERI_G4x_S_S_S_a;
  Double I_ERI_G3xy_Py_S_S_a = I_ERI_H3x2y_S_S_S_a+ABY*I_ERI_G3xy_S_S_S_a;
  Double I_ERI_G3xz_Py_S_S_a = I_ERI_H3xyz_S_S_S_a+ABY*I_ERI_G3xz_S_S_S_a;
  Double I_ERI_G2x2y_Py_S_S_a = I_ERI_H2x3y_S_S_S_a+ABY*I_ERI_G2x2y_S_S_S_a;
  Double I_ERI_G2xyz_Py_S_S_a = I_ERI_H2x2yz_S_S_S_a+ABY*I_ERI_G2xyz_S_S_S_a;
  Double I_ERI_G2x2z_Py_S_S_a = I_ERI_H2xy2z_S_S_S_a+ABY*I_ERI_G2x2z_S_S_S_a;
  Double I_ERI_Gx3y_Py_S_S_a = I_ERI_Hx4y_S_S_S_a+ABY*I_ERI_Gx3y_S_S_S_a;
  Double I_ERI_Gx2yz_Py_S_S_a = I_ERI_Hx3yz_S_S_S_a+ABY*I_ERI_Gx2yz_S_S_S_a;
  Double I_ERI_Gxy2z_Py_S_S_a = I_ERI_Hx2y2z_S_S_S_a+ABY*I_ERI_Gxy2z_S_S_S_a;
  Double I_ERI_Gx3z_Py_S_S_a = I_ERI_Hxy3z_S_S_S_a+ABY*I_ERI_Gx3z_S_S_S_a;
  Double I_ERI_G4y_Py_S_S_a = I_ERI_H5y_S_S_S_a+ABY*I_ERI_G4y_S_S_S_a;
  Double I_ERI_G3yz_Py_S_S_a = I_ERI_H4yz_S_S_S_a+ABY*I_ERI_G3yz_S_S_S_a;
  Double I_ERI_G2y2z_Py_S_S_a = I_ERI_H3y2z_S_S_S_a+ABY*I_ERI_G2y2z_S_S_S_a;
  Double I_ERI_Gy3z_Py_S_S_a = I_ERI_H2y3z_S_S_S_a+ABY*I_ERI_Gy3z_S_S_S_a;
  Double I_ERI_G4z_Py_S_S_a = I_ERI_Hy4z_S_S_S_a+ABY*I_ERI_G4z_S_S_S_a;
  Double I_ERI_G4x_Pz_S_S_a = I_ERI_H4xz_S_S_S_a+ABZ*I_ERI_G4x_S_S_S_a;
  Double I_ERI_G3xy_Pz_S_S_a = I_ERI_H3xyz_S_S_S_a+ABZ*I_ERI_G3xy_S_S_S_a;
  Double I_ERI_G3xz_Pz_S_S_a = I_ERI_H3x2z_S_S_S_a+ABZ*I_ERI_G3xz_S_S_S_a;
  Double I_ERI_G2x2y_Pz_S_S_a = I_ERI_H2x2yz_S_S_S_a+ABZ*I_ERI_G2x2y_S_S_S_a;
  Double I_ERI_G2xyz_Pz_S_S_a = I_ERI_H2xy2z_S_S_S_a+ABZ*I_ERI_G2xyz_S_S_S_a;
  Double I_ERI_G2x2z_Pz_S_S_a = I_ERI_H2x3z_S_S_S_a+ABZ*I_ERI_G2x2z_S_S_S_a;
  Double I_ERI_Gx3y_Pz_S_S_a = I_ERI_Hx3yz_S_S_S_a+ABZ*I_ERI_Gx3y_S_S_S_a;
  Double I_ERI_Gx2yz_Pz_S_S_a = I_ERI_Hx2y2z_S_S_S_a+ABZ*I_ERI_Gx2yz_S_S_S_a;
  Double I_ERI_Gxy2z_Pz_S_S_a = I_ERI_Hxy3z_S_S_S_a+ABZ*I_ERI_Gxy2z_S_S_S_a;
  Double I_ERI_Gx3z_Pz_S_S_a = I_ERI_Hx4z_S_S_S_a+ABZ*I_ERI_Gx3z_S_S_S_a;
  Double I_ERI_G4y_Pz_S_S_a = I_ERI_H4yz_S_S_S_a+ABZ*I_ERI_G4y_S_S_S_a;
  Double I_ERI_G3yz_Pz_S_S_a = I_ERI_H3y2z_S_S_S_a+ABZ*I_ERI_G3yz_S_S_S_a;
  Double I_ERI_G2y2z_Pz_S_S_a = I_ERI_H2y3z_S_S_S_a+ABZ*I_ERI_G2y2z_S_S_S_a;
  Double I_ERI_Gy3z_Pz_S_S_a = I_ERI_Hy4z_S_S_S_a+ABZ*I_ERI_Gy3z_S_S_S_a;
  Double I_ERI_G4z_Pz_S_S_a = I_ERI_H5z_S_S_S_a+ABZ*I_ERI_G4z_S_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_F_P_S_S_a
   ************************************************************/
  Double I_ERI_F3x_D2x_S_S_a = I_ERI_G4x_Px_S_S_a+ABX*I_ERI_F3x_Px_S_S_a;
  Double I_ERI_F2xy_D2x_S_S_a = I_ERI_G3xy_Px_S_S_a+ABX*I_ERI_F2xy_Px_S_S_a;
  Double I_ERI_F2xz_D2x_S_S_a = I_ERI_G3xz_Px_S_S_a+ABX*I_ERI_F2xz_Px_S_S_a;
  Double I_ERI_Fx2y_D2x_S_S_a = I_ERI_G2x2y_Px_S_S_a+ABX*I_ERI_Fx2y_Px_S_S_a;
  Double I_ERI_Fxyz_D2x_S_S_a = I_ERI_G2xyz_Px_S_S_a+ABX*I_ERI_Fxyz_Px_S_S_a;
  Double I_ERI_Fx2z_D2x_S_S_a = I_ERI_G2x2z_Px_S_S_a+ABX*I_ERI_Fx2z_Px_S_S_a;
  Double I_ERI_F3y_D2x_S_S_a = I_ERI_Gx3y_Px_S_S_a+ABX*I_ERI_F3y_Px_S_S_a;
  Double I_ERI_F2yz_D2x_S_S_a = I_ERI_Gx2yz_Px_S_S_a+ABX*I_ERI_F2yz_Px_S_S_a;
  Double I_ERI_Fy2z_D2x_S_S_a = I_ERI_Gxy2z_Px_S_S_a+ABX*I_ERI_Fy2z_Px_S_S_a;
  Double I_ERI_F3z_D2x_S_S_a = I_ERI_Gx3z_Px_S_S_a+ABX*I_ERI_F3z_Px_S_S_a;
  Double I_ERI_F3x_Dxy_S_S_a = I_ERI_G3xy_Px_S_S_a+ABY*I_ERI_F3x_Px_S_S_a;
  Double I_ERI_F2xy_Dxy_S_S_a = I_ERI_G2x2y_Px_S_S_a+ABY*I_ERI_F2xy_Px_S_S_a;
  Double I_ERI_F2xz_Dxy_S_S_a = I_ERI_G2xyz_Px_S_S_a+ABY*I_ERI_F2xz_Px_S_S_a;
  Double I_ERI_Fx2y_Dxy_S_S_a = I_ERI_Gx3y_Px_S_S_a+ABY*I_ERI_Fx2y_Px_S_S_a;
  Double I_ERI_Fxyz_Dxy_S_S_a = I_ERI_Gx2yz_Px_S_S_a+ABY*I_ERI_Fxyz_Px_S_S_a;
  Double I_ERI_Fx2z_Dxy_S_S_a = I_ERI_Gxy2z_Px_S_S_a+ABY*I_ERI_Fx2z_Px_S_S_a;
  Double I_ERI_F3y_Dxy_S_S_a = I_ERI_G4y_Px_S_S_a+ABY*I_ERI_F3y_Px_S_S_a;
  Double I_ERI_F2yz_Dxy_S_S_a = I_ERI_G3yz_Px_S_S_a+ABY*I_ERI_F2yz_Px_S_S_a;
  Double I_ERI_Fy2z_Dxy_S_S_a = I_ERI_G2y2z_Px_S_S_a+ABY*I_ERI_Fy2z_Px_S_S_a;
  Double I_ERI_F3z_Dxy_S_S_a = I_ERI_Gy3z_Px_S_S_a+ABY*I_ERI_F3z_Px_S_S_a;
  Double I_ERI_F3x_Dxz_S_S_a = I_ERI_G3xz_Px_S_S_a+ABZ*I_ERI_F3x_Px_S_S_a;
  Double I_ERI_F2xy_Dxz_S_S_a = I_ERI_G2xyz_Px_S_S_a+ABZ*I_ERI_F2xy_Px_S_S_a;
  Double I_ERI_F2xz_Dxz_S_S_a = I_ERI_G2x2z_Px_S_S_a+ABZ*I_ERI_F2xz_Px_S_S_a;
  Double I_ERI_Fx2y_Dxz_S_S_a = I_ERI_Gx2yz_Px_S_S_a+ABZ*I_ERI_Fx2y_Px_S_S_a;
  Double I_ERI_Fxyz_Dxz_S_S_a = I_ERI_Gxy2z_Px_S_S_a+ABZ*I_ERI_Fxyz_Px_S_S_a;
  Double I_ERI_Fx2z_Dxz_S_S_a = I_ERI_Gx3z_Px_S_S_a+ABZ*I_ERI_Fx2z_Px_S_S_a;
  Double I_ERI_F3y_Dxz_S_S_a = I_ERI_G3yz_Px_S_S_a+ABZ*I_ERI_F3y_Px_S_S_a;
  Double I_ERI_F2yz_Dxz_S_S_a = I_ERI_G2y2z_Px_S_S_a+ABZ*I_ERI_F2yz_Px_S_S_a;
  Double I_ERI_Fy2z_Dxz_S_S_a = I_ERI_Gy3z_Px_S_S_a+ABZ*I_ERI_Fy2z_Px_S_S_a;
  Double I_ERI_F3z_Dxz_S_S_a = I_ERI_G4z_Px_S_S_a+ABZ*I_ERI_F3z_Px_S_S_a;
  Double I_ERI_F3x_D2y_S_S_a = I_ERI_G3xy_Py_S_S_a+ABY*I_ERI_F3x_Py_S_S_a;
  Double I_ERI_F2xy_D2y_S_S_a = I_ERI_G2x2y_Py_S_S_a+ABY*I_ERI_F2xy_Py_S_S_a;
  Double I_ERI_F2xz_D2y_S_S_a = I_ERI_G2xyz_Py_S_S_a+ABY*I_ERI_F2xz_Py_S_S_a;
  Double I_ERI_Fx2y_D2y_S_S_a = I_ERI_Gx3y_Py_S_S_a+ABY*I_ERI_Fx2y_Py_S_S_a;
  Double I_ERI_Fxyz_D2y_S_S_a = I_ERI_Gx2yz_Py_S_S_a+ABY*I_ERI_Fxyz_Py_S_S_a;
  Double I_ERI_Fx2z_D2y_S_S_a = I_ERI_Gxy2z_Py_S_S_a+ABY*I_ERI_Fx2z_Py_S_S_a;
  Double I_ERI_F3y_D2y_S_S_a = I_ERI_G4y_Py_S_S_a+ABY*I_ERI_F3y_Py_S_S_a;
  Double I_ERI_F2yz_D2y_S_S_a = I_ERI_G3yz_Py_S_S_a+ABY*I_ERI_F2yz_Py_S_S_a;
  Double I_ERI_Fy2z_D2y_S_S_a = I_ERI_G2y2z_Py_S_S_a+ABY*I_ERI_Fy2z_Py_S_S_a;
  Double I_ERI_F3z_D2y_S_S_a = I_ERI_Gy3z_Py_S_S_a+ABY*I_ERI_F3z_Py_S_S_a;
  Double I_ERI_F3x_Dyz_S_S_a = I_ERI_G3xz_Py_S_S_a+ABZ*I_ERI_F3x_Py_S_S_a;
  Double I_ERI_F2xy_Dyz_S_S_a = I_ERI_G2xyz_Py_S_S_a+ABZ*I_ERI_F2xy_Py_S_S_a;
  Double I_ERI_F2xz_Dyz_S_S_a = I_ERI_G2x2z_Py_S_S_a+ABZ*I_ERI_F2xz_Py_S_S_a;
  Double I_ERI_Fx2y_Dyz_S_S_a = I_ERI_Gx2yz_Py_S_S_a+ABZ*I_ERI_Fx2y_Py_S_S_a;
  Double I_ERI_Fxyz_Dyz_S_S_a = I_ERI_Gxy2z_Py_S_S_a+ABZ*I_ERI_Fxyz_Py_S_S_a;
  Double I_ERI_Fx2z_Dyz_S_S_a = I_ERI_Gx3z_Py_S_S_a+ABZ*I_ERI_Fx2z_Py_S_S_a;
  Double I_ERI_F3y_Dyz_S_S_a = I_ERI_G3yz_Py_S_S_a+ABZ*I_ERI_F3y_Py_S_S_a;
  Double I_ERI_F2yz_Dyz_S_S_a = I_ERI_G2y2z_Py_S_S_a+ABZ*I_ERI_F2yz_Py_S_S_a;
  Double I_ERI_Fy2z_Dyz_S_S_a = I_ERI_Gy3z_Py_S_S_a+ABZ*I_ERI_Fy2z_Py_S_S_a;
  Double I_ERI_F3z_Dyz_S_S_a = I_ERI_G4z_Py_S_S_a+ABZ*I_ERI_F3z_Py_S_S_a;
  Double I_ERI_F3x_D2z_S_S_a = I_ERI_G3xz_Pz_S_S_a+ABZ*I_ERI_F3x_Pz_S_S_a;
  Double I_ERI_F2xy_D2z_S_S_a = I_ERI_G2xyz_Pz_S_S_a+ABZ*I_ERI_F2xy_Pz_S_S_a;
  Double I_ERI_F2xz_D2z_S_S_a = I_ERI_G2x2z_Pz_S_S_a+ABZ*I_ERI_F2xz_Pz_S_S_a;
  Double I_ERI_Fx2y_D2z_S_S_a = I_ERI_Gx2yz_Pz_S_S_a+ABZ*I_ERI_Fx2y_Pz_S_S_a;
  Double I_ERI_Fxyz_D2z_S_S_a = I_ERI_Gxy2z_Pz_S_S_a+ABZ*I_ERI_Fxyz_Pz_S_S_a;
  Double I_ERI_Fx2z_D2z_S_S_a = I_ERI_Gx3z_Pz_S_S_a+ABZ*I_ERI_Fx2z_Pz_S_S_a;
  Double I_ERI_F3y_D2z_S_S_a = I_ERI_G3yz_Pz_S_S_a+ABZ*I_ERI_F3y_Pz_S_S_a;
  Double I_ERI_F2yz_D2z_S_S_a = I_ERI_G2y2z_Pz_S_S_a+ABZ*I_ERI_F2yz_Pz_S_S_a;
  Double I_ERI_Fy2z_D2z_S_S_a = I_ERI_Gy3z_Pz_S_S_a+ABZ*I_ERI_Fy2z_Pz_S_S_a;
  Double I_ERI_F3z_D2z_S_S_a = I_ERI_G4z_Pz_S_S_a+ABZ*I_ERI_F3z_Pz_S_S_a;

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
   * shell quartet name: SQ_ERI_F_P_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
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
  Double I_ERI_F3x_Py_S_S_b = I_ERI_G3xy_S_S_S_b+ABY*I_ERI_F3x_S_S_S_b;
  Double I_ERI_F2xy_Py_S_S_b = I_ERI_G2x2y_S_S_S_b+ABY*I_ERI_F2xy_S_S_S_b;
  Double I_ERI_F2xz_Py_S_S_b = I_ERI_G2xyz_S_S_S_b+ABY*I_ERI_F2xz_S_S_S_b;
  Double I_ERI_Fx2y_Py_S_S_b = I_ERI_Gx3y_S_S_S_b+ABY*I_ERI_Fx2y_S_S_S_b;
  Double I_ERI_Fxyz_Py_S_S_b = I_ERI_Gx2yz_S_S_S_b+ABY*I_ERI_Fxyz_S_S_S_b;
  Double I_ERI_Fx2z_Py_S_S_b = I_ERI_Gxy2z_S_S_S_b+ABY*I_ERI_Fx2z_S_S_S_b;
  Double I_ERI_F3y_Py_S_S_b = I_ERI_G4y_S_S_S_b+ABY*I_ERI_F3y_S_S_S_b;
  Double I_ERI_F2yz_Py_S_S_b = I_ERI_G3yz_S_S_S_b+ABY*I_ERI_F2yz_S_S_S_b;
  Double I_ERI_Fy2z_Py_S_S_b = I_ERI_G2y2z_S_S_S_b+ABY*I_ERI_Fy2z_S_S_S_b;
  Double I_ERI_F3z_Py_S_S_b = I_ERI_Gy3z_S_S_S_b+ABY*I_ERI_F3z_S_S_S_b;
  Double I_ERI_F3x_Pz_S_S_b = I_ERI_G3xz_S_S_S_b+ABZ*I_ERI_F3x_S_S_S_b;
  Double I_ERI_F2xy_Pz_S_S_b = I_ERI_G2xyz_S_S_S_b+ABZ*I_ERI_F2xy_S_S_S_b;
  Double I_ERI_F2xz_Pz_S_S_b = I_ERI_G2x2z_S_S_S_b+ABZ*I_ERI_F2xz_S_S_S_b;
  Double I_ERI_Fx2y_Pz_S_S_b = I_ERI_Gx2yz_S_S_S_b+ABZ*I_ERI_Fx2y_S_S_S_b;
  Double I_ERI_Fxyz_Pz_S_S_b = I_ERI_Gxy2z_S_S_S_b+ABZ*I_ERI_Fxyz_S_S_S_b;
  Double I_ERI_Fx2z_Pz_S_S_b = I_ERI_Gx3z_S_S_S_b+ABZ*I_ERI_Fx2z_S_S_S_b;
  Double I_ERI_F3y_Pz_S_S_b = I_ERI_G3yz_S_S_S_b+ABZ*I_ERI_F3y_S_S_S_b;
  Double I_ERI_F2yz_Pz_S_S_b = I_ERI_G2y2z_S_S_S_b+ABZ*I_ERI_F2yz_S_S_S_b;
  Double I_ERI_Fy2z_Pz_S_S_b = I_ERI_Gy3z_S_S_S_b+ABZ*I_ERI_Fy2z_S_S_S_b;
  Double I_ERI_F3z_Pz_S_S_b = I_ERI_G4z_S_S_S_b+ABZ*I_ERI_F3z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
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
  Double I_ERI_D2x_D2y_S_S_b = I_ERI_F2xy_Py_S_S_b+ABY*I_ERI_D2x_Py_S_S_b;
  Double I_ERI_Dxy_D2y_S_S_b = I_ERI_Fx2y_Py_S_S_b+ABY*I_ERI_Dxy_Py_S_S_b;
  Double I_ERI_Dxz_D2y_S_S_b = I_ERI_Fxyz_Py_S_S_b+ABY*I_ERI_Dxz_Py_S_S_b;
  Double I_ERI_D2y_D2y_S_S_b = I_ERI_F3y_Py_S_S_b+ABY*I_ERI_D2y_Py_S_S_b;
  Double I_ERI_Dyz_D2y_S_S_b = I_ERI_F2yz_Py_S_S_b+ABY*I_ERI_Dyz_Py_S_S_b;
  Double I_ERI_D2z_D2y_S_S_b = I_ERI_Fy2z_Py_S_S_b+ABY*I_ERI_D2z_Py_S_S_b;
  Double I_ERI_D2x_D2z_S_S_b = I_ERI_F2xz_Pz_S_S_b+ABZ*I_ERI_D2x_Pz_S_S_b;
  Double I_ERI_Dxy_D2z_S_S_b = I_ERI_Fxyz_Pz_S_S_b+ABZ*I_ERI_Dxy_Pz_S_S_b;
  Double I_ERI_Dxz_D2z_S_S_b = I_ERI_Fx2z_Pz_S_S_b+ABZ*I_ERI_Dxz_Pz_S_S_b;
  Double I_ERI_D2y_D2z_S_S_b = I_ERI_F2yz_Pz_S_S_b+ABZ*I_ERI_D2y_Pz_S_S_b;
  Double I_ERI_Dyz_D2z_S_S_b = I_ERI_Fy2z_Pz_S_S_b+ABZ*I_ERI_Dyz_Pz_S_S_b;
  Double I_ERI_D2z_D2z_S_S_b = I_ERI_F3z_Pz_S_S_b+ABZ*I_ERI_D2z_Pz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_S_S_b
   * RHS shell quartet name: SQ_ERI_G_S_S_S_b
   ************************************************************/
  Double I_ERI_G4x_Px_S_S_b = I_ERI_H5x_S_S_S_b+ABX*I_ERI_G4x_S_S_S_b;
  Double I_ERI_G3xy_Px_S_S_b = I_ERI_H4xy_S_S_S_b+ABX*I_ERI_G3xy_S_S_S_b;
  Double I_ERI_G3xz_Px_S_S_b = I_ERI_H4xz_S_S_S_b+ABX*I_ERI_G3xz_S_S_S_b;
  Double I_ERI_G2x2y_Px_S_S_b = I_ERI_H3x2y_S_S_S_b+ABX*I_ERI_G2x2y_S_S_S_b;
  Double I_ERI_G2xyz_Px_S_S_b = I_ERI_H3xyz_S_S_S_b+ABX*I_ERI_G2xyz_S_S_S_b;
  Double I_ERI_G2x2z_Px_S_S_b = I_ERI_H3x2z_S_S_S_b+ABX*I_ERI_G2x2z_S_S_S_b;
  Double I_ERI_Gx3y_Px_S_S_b = I_ERI_H2x3y_S_S_S_b+ABX*I_ERI_Gx3y_S_S_S_b;
  Double I_ERI_Gx2yz_Px_S_S_b = I_ERI_H2x2yz_S_S_S_b+ABX*I_ERI_Gx2yz_S_S_S_b;
  Double I_ERI_Gxy2z_Px_S_S_b = I_ERI_H2xy2z_S_S_S_b+ABX*I_ERI_Gxy2z_S_S_S_b;
  Double I_ERI_Gx3z_Px_S_S_b = I_ERI_H2x3z_S_S_S_b+ABX*I_ERI_Gx3z_S_S_S_b;
  Double I_ERI_G4y_Px_S_S_b = I_ERI_Hx4y_S_S_S_b+ABX*I_ERI_G4y_S_S_S_b;
  Double I_ERI_G3yz_Px_S_S_b = I_ERI_Hx3yz_S_S_S_b+ABX*I_ERI_G3yz_S_S_S_b;
  Double I_ERI_G2y2z_Px_S_S_b = I_ERI_Hx2y2z_S_S_S_b+ABX*I_ERI_G2y2z_S_S_S_b;
  Double I_ERI_Gy3z_Px_S_S_b = I_ERI_Hxy3z_S_S_S_b+ABX*I_ERI_Gy3z_S_S_S_b;
  Double I_ERI_G4z_Px_S_S_b = I_ERI_Hx4z_S_S_S_b+ABX*I_ERI_G4z_S_S_S_b;
  Double I_ERI_G3xy_Py_S_S_b = I_ERI_H3x2y_S_S_S_b+ABY*I_ERI_G3xy_S_S_S_b;
  Double I_ERI_G3xz_Py_S_S_b = I_ERI_H3xyz_S_S_S_b+ABY*I_ERI_G3xz_S_S_S_b;
  Double I_ERI_G2x2y_Py_S_S_b = I_ERI_H2x3y_S_S_S_b+ABY*I_ERI_G2x2y_S_S_S_b;
  Double I_ERI_G2xyz_Py_S_S_b = I_ERI_H2x2yz_S_S_S_b+ABY*I_ERI_G2xyz_S_S_S_b;
  Double I_ERI_G2x2z_Py_S_S_b = I_ERI_H2xy2z_S_S_S_b+ABY*I_ERI_G2x2z_S_S_S_b;
  Double I_ERI_Gx3y_Py_S_S_b = I_ERI_Hx4y_S_S_S_b+ABY*I_ERI_Gx3y_S_S_S_b;
  Double I_ERI_Gx2yz_Py_S_S_b = I_ERI_Hx3yz_S_S_S_b+ABY*I_ERI_Gx2yz_S_S_S_b;
  Double I_ERI_Gxy2z_Py_S_S_b = I_ERI_Hx2y2z_S_S_S_b+ABY*I_ERI_Gxy2z_S_S_S_b;
  Double I_ERI_Gx3z_Py_S_S_b = I_ERI_Hxy3z_S_S_S_b+ABY*I_ERI_Gx3z_S_S_S_b;
  Double I_ERI_G4y_Py_S_S_b = I_ERI_H5y_S_S_S_b+ABY*I_ERI_G4y_S_S_S_b;
  Double I_ERI_G3yz_Py_S_S_b = I_ERI_H4yz_S_S_S_b+ABY*I_ERI_G3yz_S_S_S_b;
  Double I_ERI_G2y2z_Py_S_S_b = I_ERI_H3y2z_S_S_S_b+ABY*I_ERI_G2y2z_S_S_S_b;
  Double I_ERI_Gy3z_Py_S_S_b = I_ERI_H2y3z_S_S_S_b+ABY*I_ERI_Gy3z_S_S_S_b;
  Double I_ERI_G4z_Py_S_S_b = I_ERI_Hy4z_S_S_S_b+ABY*I_ERI_G4z_S_S_S_b;
  Double I_ERI_G3xz_Pz_S_S_b = I_ERI_H3x2z_S_S_S_b+ABZ*I_ERI_G3xz_S_S_S_b;
  Double I_ERI_G2xyz_Pz_S_S_b = I_ERI_H2xy2z_S_S_S_b+ABZ*I_ERI_G2xyz_S_S_S_b;
  Double I_ERI_G2x2z_Pz_S_S_b = I_ERI_H2x3z_S_S_S_b+ABZ*I_ERI_G2x2z_S_S_S_b;
  Double I_ERI_Gx2yz_Pz_S_S_b = I_ERI_Hx2y2z_S_S_S_b+ABZ*I_ERI_Gx2yz_S_S_S_b;
  Double I_ERI_Gxy2z_Pz_S_S_b = I_ERI_Hxy3z_S_S_S_b+ABZ*I_ERI_Gxy2z_S_S_S_b;
  Double I_ERI_Gx3z_Pz_S_S_b = I_ERI_Hx4z_S_S_S_b+ABZ*I_ERI_Gx3z_S_S_S_b;
  Double I_ERI_G3yz_Pz_S_S_b = I_ERI_H3y2z_S_S_S_b+ABZ*I_ERI_G3yz_S_S_S_b;
  Double I_ERI_G2y2z_Pz_S_S_b = I_ERI_H2y3z_S_S_S_b+ABZ*I_ERI_G2y2z_S_S_S_b;
  Double I_ERI_Gy3z_Pz_S_S_b = I_ERI_Hy4z_S_S_S_b+ABZ*I_ERI_Gy3z_S_S_S_b;
  Double I_ERI_G4z_Pz_S_S_b = I_ERI_H5z_S_S_S_b+ABZ*I_ERI_G4z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_b
   * RHS shell quartet name: SQ_ERI_F_P_S_S_b
   ************************************************************/
  Double I_ERI_F3x_D2x_S_S_b = I_ERI_G4x_Px_S_S_b+ABX*I_ERI_F3x_Px_S_S_b;
  Double I_ERI_F2xy_D2x_S_S_b = I_ERI_G3xy_Px_S_S_b+ABX*I_ERI_F2xy_Px_S_S_b;
  Double I_ERI_F2xz_D2x_S_S_b = I_ERI_G3xz_Px_S_S_b+ABX*I_ERI_F2xz_Px_S_S_b;
  Double I_ERI_Fx2y_D2x_S_S_b = I_ERI_G2x2y_Px_S_S_b+ABX*I_ERI_Fx2y_Px_S_S_b;
  Double I_ERI_Fxyz_D2x_S_S_b = I_ERI_G2xyz_Px_S_S_b+ABX*I_ERI_Fxyz_Px_S_S_b;
  Double I_ERI_Fx2z_D2x_S_S_b = I_ERI_G2x2z_Px_S_S_b+ABX*I_ERI_Fx2z_Px_S_S_b;
  Double I_ERI_F3y_D2x_S_S_b = I_ERI_Gx3y_Px_S_S_b+ABX*I_ERI_F3y_Px_S_S_b;
  Double I_ERI_F2yz_D2x_S_S_b = I_ERI_Gx2yz_Px_S_S_b+ABX*I_ERI_F2yz_Px_S_S_b;
  Double I_ERI_Fy2z_D2x_S_S_b = I_ERI_Gxy2z_Px_S_S_b+ABX*I_ERI_Fy2z_Px_S_S_b;
  Double I_ERI_F3z_D2x_S_S_b = I_ERI_Gx3z_Px_S_S_b+ABX*I_ERI_F3z_Px_S_S_b;
  Double I_ERI_F3x_Dxy_S_S_b = I_ERI_G3xy_Px_S_S_b+ABY*I_ERI_F3x_Px_S_S_b;
  Double I_ERI_F2xy_Dxy_S_S_b = I_ERI_G2x2y_Px_S_S_b+ABY*I_ERI_F2xy_Px_S_S_b;
  Double I_ERI_F2xz_Dxy_S_S_b = I_ERI_G2xyz_Px_S_S_b+ABY*I_ERI_F2xz_Px_S_S_b;
  Double I_ERI_Fx2y_Dxy_S_S_b = I_ERI_Gx3y_Px_S_S_b+ABY*I_ERI_Fx2y_Px_S_S_b;
  Double I_ERI_Fxyz_Dxy_S_S_b = I_ERI_Gx2yz_Px_S_S_b+ABY*I_ERI_Fxyz_Px_S_S_b;
  Double I_ERI_Fx2z_Dxy_S_S_b = I_ERI_Gxy2z_Px_S_S_b+ABY*I_ERI_Fx2z_Px_S_S_b;
  Double I_ERI_F3y_Dxy_S_S_b = I_ERI_G4y_Px_S_S_b+ABY*I_ERI_F3y_Px_S_S_b;
  Double I_ERI_F2yz_Dxy_S_S_b = I_ERI_G3yz_Px_S_S_b+ABY*I_ERI_F2yz_Px_S_S_b;
  Double I_ERI_Fy2z_Dxy_S_S_b = I_ERI_G2y2z_Px_S_S_b+ABY*I_ERI_Fy2z_Px_S_S_b;
  Double I_ERI_F3z_Dxy_S_S_b = I_ERI_Gy3z_Px_S_S_b+ABY*I_ERI_F3z_Px_S_S_b;
  Double I_ERI_F3x_Dxz_S_S_b = I_ERI_G3xz_Px_S_S_b+ABZ*I_ERI_F3x_Px_S_S_b;
  Double I_ERI_F2xy_Dxz_S_S_b = I_ERI_G2xyz_Px_S_S_b+ABZ*I_ERI_F2xy_Px_S_S_b;
  Double I_ERI_F2xz_Dxz_S_S_b = I_ERI_G2x2z_Px_S_S_b+ABZ*I_ERI_F2xz_Px_S_S_b;
  Double I_ERI_Fx2y_Dxz_S_S_b = I_ERI_Gx2yz_Px_S_S_b+ABZ*I_ERI_Fx2y_Px_S_S_b;
  Double I_ERI_Fxyz_Dxz_S_S_b = I_ERI_Gxy2z_Px_S_S_b+ABZ*I_ERI_Fxyz_Px_S_S_b;
  Double I_ERI_Fx2z_Dxz_S_S_b = I_ERI_Gx3z_Px_S_S_b+ABZ*I_ERI_Fx2z_Px_S_S_b;
  Double I_ERI_F3y_Dxz_S_S_b = I_ERI_G3yz_Px_S_S_b+ABZ*I_ERI_F3y_Px_S_S_b;
  Double I_ERI_F2yz_Dxz_S_S_b = I_ERI_G2y2z_Px_S_S_b+ABZ*I_ERI_F2yz_Px_S_S_b;
  Double I_ERI_Fy2z_Dxz_S_S_b = I_ERI_Gy3z_Px_S_S_b+ABZ*I_ERI_Fy2z_Px_S_S_b;
  Double I_ERI_F3z_Dxz_S_S_b = I_ERI_G4z_Px_S_S_b+ABZ*I_ERI_F3z_Px_S_S_b;
  Double I_ERI_F3x_D2y_S_S_b = I_ERI_G3xy_Py_S_S_b+ABY*I_ERI_F3x_Py_S_S_b;
  Double I_ERI_F2xy_D2y_S_S_b = I_ERI_G2x2y_Py_S_S_b+ABY*I_ERI_F2xy_Py_S_S_b;
  Double I_ERI_F2xz_D2y_S_S_b = I_ERI_G2xyz_Py_S_S_b+ABY*I_ERI_F2xz_Py_S_S_b;
  Double I_ERI_Fx2y_D2y_S_S_b = I_ERI_Gx3y_Py_S_S_b+ABY*I_ERI_Fx2y_Py_S_S_b;
  Double I_ERI_Fxyz_D2y_S_S_b = I_ERI_Gx2yz_Py_S_S_b+ABY*I_ERI_Fxyz_Py_S_S_b;
  Double I_ERI_Fx2z_D2y_S_S_b = I_ERI_Gxy2z_Py_S_S_b+ABY*I_ERI_Fx2z_Py_S_S_b;
  Double I_ERI_F3y_D2y_S_S_b = I_ERI_G4y_Py_S_S_b+ABY*I_ERI_F3y_Py_S_S_b;
  Double I_ERI_F2yz_D2y_S_S_b = I_ERI_G3yz_Py_S_S_b+ABY*I_ERI_F2yz_Py_S_S_b;
  Double I_ERI_Fy2z_D2y_S_S_b = I_ERI_G2y2z_Py_S_S_b+ABY*I_ERI_Fy2z_Py_S_S_b;
  Double I_ERI_F3z_D2y_S_S_b = I_ERI_Gy3z_Py_S_S_b+ABY*I_ERI_F3z_Py_S_S_b;
  Double I_ERI_F3x_Dyz_S_S_b = I_ERI_G3xz_Py_S_S_b+ABZ*I_ERI_F3x_Py_S_S_b;
  Double I_ERI_F2xy_Dyz_S_S_b = I_ERI_G2xyz_Py_S_S_b+ABZ*I_ERI_F2xy_Py_S_S_b;
  Double I_ERI_F2xz_Dyz_S_S_b = I_ERI_G2x2z_Py_S_S_b+ABZ*I_ERI_F2xz_Py_S_S_b;
  Double I_ERI_Fx2y_Dyz_S_S_b = I_ERI_Gx2yz_Py_S_S_b+ABZ*I_ERI_Fx2y_Py_S_S_b;
  Double I_ERI_Fxyz_Dyz_S_S_b = I_ERI_Gxy2z_Py_S_S_b+ABZ*I_ERI_Fxyz_Py_S_S_b;
  Double I_ERI_Fx2z_Dyz_S_S_b = I_ERI_Gx3z_Py_S_S_b+ABZ*I_ERI_Fx2z_Py_S_S_b;
  Double I_ERI_F3y_Dyz_S_S_b = I_ERI_G3yz_Py_S_S_b+ABZ*I_ERI_F3y_Py_S_S_b;
  Double I_ERI_F2yz_Dyz_S_S_b = I_ERI_G2y2z_Py_S_S_b+ABZ*I_ERI_F2yz_Py_S_S_b;
  Double I_ERI_Fy2z_Dyz_S_S_b = I_ERI_Gy3z_Py_S_S_b+ABZ*I_ERI_Fy2z_Py_S_S_b;
  Double I_ERI_F3z_Dyz_S_S_b = I_ERI_G4z_Py_S_S_b+ABZ*I_ERI_F3z_Py_S_S_b;
  Double I_ERI_F3x_D2z_S_S_b = I_ERI_G3xz_Pz_S_S_b+ABZ*I_ERI_F3x_Pz_S_S_b;
  Double I_ERI_F2xy_D2z_S_S_b = I_ERI_G2xyz_Pz_S_S_b+ABZ*I_ERI_F2xy_Pz_S_S_b;
  Double I_ERI_F2xz_D2z_S_S_b = I_ERI_G2x2z_Pz_S_S_b+ABZ*I_ERI_F2xz_Pz_S_S_b;
  Double I_ERI_Fx2y_D2z_S_S_b = I_ERI_Gx2yz_Pz_S_S_b+ABZ*I_ERI_Fx2y_Pz_S_S_b;
  Double I_ERI_Fxyz_D2z_S_S_b = I_ERI_Gxy2z_Pz_S_S_b+ABZ*I_ERI_Fxyz_Pz_S_S_b;
  Double I_ERI_Fx2z_D2z_S_S_b = I_ERI_Gx3z_Pz_S_S_b+ABZ*I_ERI_Fx2z_Pz_S_S_b;
  Double I_ERI_F3y_D2z_S_S_b = I_ERI_G3yz_Pz_S_S_b+ABZ*I_ERI_F3y_Pz_S_S_b;
  Double I_ERI_F2yz_D2z_S_S_b = I_ERI_G2y2z_Pz_S_S_b+ABZ*I_ERI_F2yz_Pz_S_S_b;
  Double I_ERI_Fy2z_D2z_S_S_b = I_ERI_Gy3z_Pz_S_S_b+ABZ*I_ERI_Fy2z_Pz_S_S_b;
  Double I_ERI_F3z_D2z_S_S_b = I_ERI_G4z_Pz_S_S_b+ABZ*I_ERI_F3z_Pz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_F_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_D_D_S_S_b
   ************************************************************/
  Double I_ERI_D2x_F3x_S_S_b = I_ERI_F3x_D2x_S_S_b+ABX*I_ERI_D2x_D2x_S_S_b;
  Double I_ERI_Dxy_F3x_S_S_b = I_ERI_F2xy_D2x_S_S_b+ABX*I_ERI_Dxy_D2x_S_S_b;
  Double I_ERI_Dxz_F3x_S_S_b = I_ERI_F2xz_D2x_S_S_b+ABX*I_ERI_Dxz_D2x_S_S_b;
  Double I_ERI_D2y_F3x_S_S_b = I_ERI_Fx2y_D2x_S_S_b+ABX*I_ERI_D2y_D2x_S_S_b;
  Double I_ERI_Dyz_F3x_S_S_b = I_ERI_Fxyz_D2x_S_S_b+ABX*I_ERI_Dyz_D2x_S_S_b;
  Double I_ERI_D2z_F3x_S_S_b = I_ERI_Fx2z_D2x_S_S_b+ABX*I_ERI_D2z_D2x_S_S_b;
  Double I_ERI_D2x_F2xy_S_S_b = I_ERI_F2xy_D2x_S_S_b+ABY*I_ERI_D2x_D2x_S_S_b;
  Double I_ERI_Dxy_F2xy_S_S_b = I_ERI_Fx2y_D2x_S_S_b+ABY*I_ERI_Dxy_D2x_S_S_b;
  Double I_ERI_Dxz_F2xy_S_S_b = I_ERI_Fxyz_D2x_S_S_b+ABY*I_ERI_Dxz_D2x_S_S_b;
  Double I_ERI_D2y_F2xy_S_S_b = I_ERI_F3y_D2x_S_S_b+ABY*I_ERI_D2y_D2x_S_S_b;
  Double I_ERI_Dyz_F2xy_S_S_b = I_ERI_F2yz_D2x_S_S_b+ABY*I_ERI_Dyz_D2x_S_S_b;
  Double I_ERI_D2z_F2xy_S_S_b = I_ERI_Fy2z_D2x_S_S_b+ABY*I_ERI_D2z_D2x_S_S_b;
  Double I_ERI_D2x_F2xz_S_S_b = I_ERI_F2xz_D2x_S_S_b+ABZ*I_ERI_D2x_D2x_S_S_b;
  Double I_ERI_Dxy_F2xz_S_S_b = I_ERI_Fxyz_D2x_S_S_b+ABZ*I_ERI_Dxy_D2x_S_S_b;
  Double I_ERI_Dxz_F2xz_S_S_b = I_ERI_Fx2z_D2x_S_S_b+ABZ*I_ERI_Dxz_D2x_S_S_b;
  Double I_ERI_D2y_F2xz_S_S_b = I_ERI_F2yz_D2x_S_S_b+ABZ*I_ERI_D2y_D2x_S_S_b;
  Double I_ERI_Dyz_F2xz_S_S_b = I_ERI_Fy2z_D2x_S_S_b+ABZ*I_ERI_Dyz_D2x_S_S_b;
  Double I_ERI_D2z_F2xz_S_S_b = I_ERI_F3z_D2x_S_S_b+ABZ*I_ERI_D2z_D2x_S_S_b;
  Double I_ERI_D2x_Fx2y_S_S_b = I_ERI_F3x_D2y_S_S_b+ABX*I_ERI_D2x_D2y_S_S_b;
  Double I_ERI_Dxy_Fx2y_S_S_b = I_ERI_F2xy_D2y_S_S_b+ABX*I_ERI_Dxy_D2y_S_S_b;
  Double I_ERI_Dxz_Fx2y_S_S_b = I_ERI_F2xz_D2y_S_S_b+ABX*I_ERI_Dxz_D2y_S_S_b;
  Double I_ERI_D2y_Fx2y_S_S_b = I_ERI_Fx2y_D2y_S_S_b+ABX*I_ERI_D2y_D2y_S_S_b;
  Double I_ERI_Dyz_Fx2y_S_S_b = I_ERI_Fxyz_D2y_S_S_b+ABX*I_ERI_Dyz_D2y_S_S_b;
  Double I_ERI_D2z_Fx2y_S_S_b = I_ERI_Fx2z_D2y_S_S_b+ABX*I_ERI_D2z_D2y_S_S_b;
  Double I_ERI_D2x_Fxyz_S_S_b = I_ERI_F2xz_Dxy_S_S_b+ABZ*I_ERI_D2x_Dxy_S_S_b;
  Double I_ERI_Dxy_Fxyz_S_S_b = I_ERI_Fxyz_Dxy_S_S_b+ABZ*I_ERI_Dxy_Dxy_S_S_b;
  Double I_ERI_Dxz_Fxyz_S_S_b = I_ERI_Fx2z_Dxy_S_S_b+ABZ*I_ERI_Dxz_Dxy_S_S_b;
  Double I_ERI_D2y_Fxyz_S_S_b = I_ERI_F2yz_Dxy_S_S_b+ABZ*I_ERI_D2y_Dxy_S_S_b;
  Double I_ERI_Dyz_Fxyz_S_S_b = I_ERI_Fy2z_Dxy_S_S_b+ABZ*I_ERI_Dyz_Dxy_S_S_b;
  Double I_ERI_D2z_Fxyz_S_S_b = I_ERI_F3z_Dxy_S_S_b+ABZ*I_ERI_D2z_Dxy_S_S_b;
  Double I_ERI_D2x_Fx2z_S_S_b = I_ERI_F3x_D2z_S_S_b+ABX*I_ERI_D2x_D2z_S_S_b;
  Double I_ERI_Dxy_Fx2z_S_S_b = I_ERI_F2xy_D2z_S_S_b+ABX*I_ERI_Dxy_D2z_S_S_b;
  Double I_ERI_Dxz_Fx2z_S_S_b = I_ERI_F2xz_D2z_S_S_b+ABX*I_ERI_Dxz_D2z_S_S_b;
  Double I_ERI_D2y_Fx2z_S_S_b = I_ERI_Fx2y_D2z_S_S_b+ABX*I_ERI_D2y_D2z_S_S_b;
  Double I_ERI_Dyz_Fx2z_S_S_b = I_ERI_Fxyz_D2z_S_S_b+ABX*I_ERI_Dyz_D2z_S_S_b;
  Double I_ERI_D2z_Fx2z_S_S_b = I_ERI_Fx2z_D2z_S_S_b+ABX*I_ERI_D2z_D2z_S_S_b;
  Double I_ERI_D2x_F3y_S_S_b = I_ERI_F2xy_D2y_S_S_b+ABY*I_ERI_D2x_D2y_S_S_b;
  Double I_ERI_Dxy_F3y_S_S_b = I_ERI_Fx2y_D2y_S_S_b+ABY*I_ERI_Dxy_D2y_S_S_b;
  Double I_ERI_Dxz_F3y_S_S_b = I_ERI_Fxyz_D2y_S_S_b+ABY*I_ERI_Dxz_D2y_S_S_b;
  Double I_ERI_D2y_F3y_S_S_b = I_ERI_F3y_D2y_S_S_b+ABY*I_ERI_D2y_D2y_S_S_b;
  Double I_ERI_Dyz_F3y_S_S_b = I_ERI_F2yz_D2y_S_S_b+ABY*I_ERI_Dyz_D2y_S_S_b;
  Double I_ERI_D2z_F3y_S_S_b = I_ERI_Fy2z_D2y_S_S_b+ABY*I_ERI_D2z_D2y_S_S_b;
  Double I_ERI_D2x_F2yz_S_S_b = I_ERI_F2xz_D2y_S_S_b+ABZ*I_ERI_D2x_D2y_S_S_b;
  Double I_ERI_Dxy_F2yz_S_S_b = I_ERI_Fxyz_D2y_S_S_b+ABZ*I_ERI_Dxy_D2y_S_S_b;
  Double I_ERI_Dxz_F2yz_S_S_b = I_ERI_Fx2z_D2y_S_S_b+ABZ*I_ERI_Dxz_D2y_S_S_b;
  Double I_ERI_D2y_F2yz_S_S_b = I_ERI_F2yz_D2y_S_S_b+ABZ*I_ERI_D2y_D2y_S_S_b;
  Double I_ERI_Dyz_F2yz_S_S_b = I_ERI_Fy2z_D2y_S_S_b+ABZ*I_ERI_Dyz_D2y_S_S_b;
  Double I_ERI_D2z_F2yz_S_S_b = I_ERI_F3z_D2y_S_S_b+ABZ*I_ERI_D2z_D2y_S_S_b;
  Double I_ERI_D2x_Fy2z_S_S_b = I_ERI_F2xy_D2z_S_S_b+ABY*I_ERI_D2x_D2z_S_S_b;
  Double I_ERI_Dxy_Fy2z_S_S_b = I_ERI_Fx2y_D2z_S_S_b+ABY*I_ERI_Dxy_D2z_S_S_b;
  Double I_ERI_Dxz_Fy2z_S_S_b = I_ERI_Fxyz_D2z_S_S_b+ABY*I_ERI_Dxz_D2z_S_S_b;
  Double I_ERI_D2y_Fy2z_S_S_b = I_ERI_F3y_D2z_S_S_b+ABY*I_ERI_D2y_D2z_S_S_b;
  Double I_ERI_Dyz_Fy2z_S_S_b = I_ERI_F2yz_D2z_S_S_b+ABY*I_ERI_Dyz_D2z_S_S_b;
  Double I_ERI_D2z_Fy2z_S_S_b = I_ERI_Fy2z_D2z_S_S_b+ABY*I_ERI_D2z_D2z_S_S_b;
  Double I_ERI_D2x_F3z_S_S_b = I_ERI_F2xz_D2z_S_S_b+ABZ*I_ERI_D2x_D2z_S_S_b;
  Double I_ERI_Dxy_F3z_S_S_b = I_ERI_Fxyz_D2z_S_S_b+ABZ*I_ERI_Dxy_D2z_S_S_b;
  Double I_ERI_Dxz_F3z_S_S_b = I_ERI_Fx2z_D2z_S_S_b+ABZ*I_ERI_Dxz_D2z_S_S_b;
  Double I_ERI_D2y_F3z_S_S_b = I_ERI_F2yz_D2z_S_S_b+ABZ*I_ERI_D2y_D2z_S_S_b;
  Double I_ERI_Dyz_F3z_S_S_b = I_ERI_Fy2z_D2z_S_S_b+ABZ*I_ERI_Dyz_D2z_S_S_b;
  Double I_ERI_D2z_F3z_S_S_b = I_ERI_F3z_D2z_S_S_b+ABZ*I_ERI_D2z_D2z_S_S_b;

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
   * shell quartet name: SQ_ERI_F_P_S_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_c
   * RHS shell quartet name: SQ_ERI_F_S_S_S_c
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_c = I_ERI_G4x_S_S_S_c+ABX*I_ERI_F3x_S_S_S_c;
  Double I_ERI_F2xy_Px_S_S_c = I_ERI_G3xy_S_S_S_c+ABX*I_ERI_F2xy_S_S_S_c;
  Double I_ERI_F2xz_Px_S_S_c = I_ERI_G3xz_S_S_S_c+ABX*I_ERI_F2xz_S_S_S_c;
  Double I_ERI_Fx2y_Px_S_S_c = I_ERI_G2x2y_S_S_S_c+ABX*I_ERI_Fx2y_S_S_S_c;
  Double I_ERI_Fxyz_Px_S_S_c = I_ERI_G2xyz_S_S_S_c+ABX*I_ERI_Fxyz_S_S_S_c;
  Double I_ERI_Fx2z_Px_S_S_c = I_ERI_G2x2z_S_S_S_c+ABX*I_ERI_Fx2z_S_S_S_c;
  Double I_ERI_F3y_Px_S_S_c = I_ERI_Gx3y_S_S_S_c+ABX*I_ERI_F3y_S_S_S_c;
  Double I_ERI_F2yz_Px_S_S_c = I_ERI_Gx2yz_S_S_S_c+ABX*I_ERI_F2yz_S_S_S_c;
  Double I_ERI_Fy2z_Px_S_S_c = I_ERI_Gxy2z_S_S_S_c+ABX*I_ERI_Fy2z_S_S_S_c;
  Double I_ERI_F3z_Px_S_S_c = I_ERI_Gx3z_S_S_S_c+ABX*I_ERI_F3z_S_S_S_c;
  Double I_ERI_F3x_Py_S_S_c = I_ERI_G3xy_S_S_S_c+ABY*I_ERI_F3x_S_S_S_c;
  Double I_ERI_F2xy_Py_S_S_c = I_ERI_G2x2y_S_S_S_c+ABY*I_ERI_F2xy_S_S_S_c;
  Double I_ERI_F2xz_Py_S_S_c = I_ERI_G2xyz_S_S_S_c+ABY*I_ERI_F2xz_S_S_S_c;
  Double I_ERI_Fx2y_Py_S_S_c = I_ERI_Gx3y_S_S_S_c+ABY*I_ERI_Fx2y_S_S_S_c;
  Double I_ERI_Fxyz_Py_S_S_c = I_ERI_Gx2yz_S_S_S_c+ABY*I_ERI_Fxyz_S_S_S_c;
  Double I_ERI_Fx2z_Py_S_S_c = I_ERI_Gxy2z_S_S_S_c+ABY*I_ERI_Fx2z_S_S_S_c;
  Double I_ERI_F3y_Py_S_S_c = I_ERI_G4y_S_S_S_c+ABY*I_ERI_F3y_S_S_S_c;
  Double I_ERI_F2yz_Py_S_S_c = I_ERI_G3yz_S_S_S_c+ABY*I_ERI_F2yz_S_S_S_c;
  Double I_ERI_Fy2z_Py_S_S_c = I_ERI_G2y2z_S_S_S_c+ABY*I_ERI_Fy2z_S_S_S_c;
  Double I_ERI_F3z_Py_S_S_c = I_ERI_Gy3z_S_S_S_c+ABY*I_ERI_F3z_S_S_S_c;
  Double I_ERI_F3x_Pz_S_S_c = I_ERI_G3xz_S_S_S_c+ABZ*I_ERI_F3x_S_S_S_c;
  Double I_ERI_F2xy_Pz_S_S_c = I_ERI_G2xyz_S_S_S_c+ABZ*I_ERI_F2xy_S_S_S_c;
  Double I_ERI_F2xz_Pz_S_S_c = I_ERI_G2x2z_S_S_S_c+ABZ*I_ERI_F2xz_S_S_S_c;
  Double I_ERI_Fx2y_Pz_S_S_c = I_ERI_Gx2yz_S_S_S_c+ABZ*I_ERI_Fx2y_S_S_S_c;
  Double I_ERI_Fxyz_Pz_S_S_c = I_ERI_Gxy2z_S_S_S_c+ABZ*I_ERI_Fxyz_S_S_S_c;
  Double I_ERI_Fx2z_Pz_S_S_c = I_ERI_Gx3z_S_S_S_c+ABZ*I_ERI_Fx2z_S_S_S_c;
  Double I_ERI_F3y_Pz_S_S_c = I_ERI_G3yz_S_S_S_c+ABZ*I_ERI_F3y_S_S_S_c;
  Double I_ERI_F2yz_Pz_S_S_c = I_ERI_G2y2z_S_S_S_c+ABZ*I_ERI_F2yz_S_S_S_c;
  Double I_ERI_Fy2z_Pz_S_S_c = I_ERI_Gy3z_S_S_S_c+ABZ*I_ERI_Fy2z_S_S_S_c;
  Double I_ERI_F3z_Pz_S_S_c = I_ERI_G4z_S_S_S_c+ABZ*I_ERI_F3z_S_S_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_c
   * RHS shell quartet name: SQ_ERI_F_S_P_S_c
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_c = I_ERI_G4x_S_Px_S_c+ABX*I_ERI_F3x_S_Px_S_c;
  Double I_ERI_F2xy_Px_Px_S_c = I_ERI_G3xy_S_Px_S_c+ABX*I_ERI_F2xy_S_Px_S_c;
  Double I_ERI_F2xz_Px_Px_S_c = I_ERI_G3xz_S_Px_S_c+ABX*I_ERI_F2xz_S_Px_S_c;
  Double I_ERI_Fx2y_Px_Px_S_c = I_ERI_G2x2y_S_Px_S_c+ABX*I_ERI_Fx2y_S_Px_S_c;
  Double I_ERI_Fxyz_Px_Px_S_c = I_ERI_G2xyz_S_Px_S_c+ABX*I_ERI_Fxyz_S_Px_S_c;
  Double I_ERI_Fx2z_Px_Px_S_c = I_ERI_G2x2z_S_Px_S_c+ABX*I_ERI_Fx2z_S_Px_S_c;
  Double I_ERI_F3y_Px_Px_S_c = I_ERI_Gx3y_S_Px_S_c+ABX*I_ERI_F3y_S_Px_S_c;
  Double I_ERI_F2yz_Px_Px_S_c = I_ERI_Gx2yz_S_Px_S_c+ABX*I_ERI_F2yz_S_Px_S_c;
  Double I_ERI_Fy2z_Px_Px_S_c = I_ERI_Gxy2z_S_Px_S_c+ABX*I_ERI_Fy2z_S_Px_S_c;
  Double I_ERI_F3z_Px_Px_S_c = I_ERI_Gx3z_S_Px_S_c+ABX*I_ERI_F3z_S_Px_S_c;
  Double I_ERI_F3x_Py_Px_S_c = I_ERI_G3xy_S_Px_S_c+ABY*I_ERI_F3x_S_Px_S_c;
  Double I_ERI_F2xy_Py_Px_S_c = I_ERI_G2x2y_S_Px_S_c+ABY*I_ERI_F2xy_S_Px_S_c;
  Double I_ERI_F2xz_Py_Px_S_c = I_ERI_G2xyz_S_Px_S_c+ABY*I_ERI_F2xz_S_Px_S_c;
  Double I_ERI_Fx2y_Py_Px_S_c = I_ERI_Gx3y_S_Px_S_c+ABY*I_ERI_Fx2y_S_Px_S_c;
  Double I_ERI_Fxyz_Py_Px_S_c = I_ERI_Gx2yz_S_Px_S_c+ABY*I_ERI_Fxyz_S_Px_S_c;
  Double I_ERI_Fx2z_Py_Px_S_c = I_ERI_Gxy2z_S_Px_S_c+ABY*I_ERI_Fx2z_S_Px_S_c;
  Double I_ERI_F3y_Py_Px_S_c = I_ERI_G4y_S_Px_S_c+ABY*I_ERI_F3y_S_Px_S_c;
  Double I_ERI_F2yz_Py_Px_S_c = I_ERI_G3yz_S_Px_S_c+ABY*I_ERI_F2yz_S_Px_S_c;
  Double I_ERI_Fy2z_Py_Px_S_c = I_ERI_G2y2z_S_Px_S_c+ABY*I_ERI_Fy2z_S_Px_S_c;
  Double I_ERI_F3z_Py_Px_S_c = I_ERI_Gy3z_S_Px_S_c+ABY*I_ERI_F3z_S_Px_S_c;
  Double I_ERI_F3x_Pz_Px_S_c = I_ERI_G3xz_S_Px_S_c+ABZ*I_ERI_F3x_S_Px_S_c;
  Double I_ERI_F2xy_Pz_Px_S_c = I_ERI_G2xyz_S_Px_S_c+ABZ*I_ERI_F2xy_S_Px_S_c;
  Double I_ERI_F2xz_Pz_Px_S_c = I_ERI_G2x2z_S_Px_S_c+ABZ*I_ERI_F2xz_S_Px_S_c;
  Double I_ERI_Fx2y_Pz_Px_S_c = I_ERI_Gx2yz_S_Px_S_c+ABZ*I_ERI_Fx2y_S_Px_S_c;
  Double I_ERI_Fxyz_Pz_Px_S_c = I_ERI_Gxy2z_S_Px_S_c+ABZ*I_ERI_Fxyz_S_Px_S_c;
  Double I_ERI_Fx2z_Pz_Px_S_c = I_ERI_Gx3z_S_Px_S_c+ABZ*I_ERI_Fx2z_S_Px_S_c;
  Double I_ERI_F3y_Pz_Px_S_c = I_ERI_G3yz_S_Px_S_c+ABZ*I_ERI_F3y_S_Px_S_c;
  Double I_ERI_F2yz_Pz_Px_S_c = I_ERI_G2y2z_S_Px_S_c+ABZ*I_ERI_F2yz_S_Px_S_c;
  Double I_ERI_Fy2z_Pz_Px_S_c = I_ERI_Gy3z_S_Px_S_c+ABZ*I_ERI_Fy2z_S_Px_S_c;
  Double I_ERI_F3z_Pz_Px_S_c = I_ERI_G4z_S_Px_S_c+ABZ*I_ERI_F3z_S_Px_S_c;
  Double I_ERI_F3x_Px_Py_S_c = I_ERI_G4x_S_Py_S_c+ABX*I_ERI_F3x_S_Py_S_c;
  Double I_ERI_F2xy_Px_Py_S_c = I_ERI_G3xy_S_Py_S_c+ABX*I_ERI_F2xy_S_Py_S_c;
  Double I_ERI_F2xz_Px_Py_S_c = I_ERI_G3xz_S_Py_S_c+ABX*I_ERI_F2xz_S_Py_S_c;
  Double I_ERI_Fx2y_Px_Py_S_c = I_ERI_G2x2y_S_Py_S_c+ABX*I_ERI_Fx2y_S_Py_S_c;
  Double I_ERI_Fxyz_Px_Py_S_c = I_ERI_G2xyz_S_Py_S_c+ABX*I_ERI_Fxyz_S_Py_S_c;
  Double I_ERI_Fx2z_Px_Py_S_c = I_ERI_G2x2z_S_Py_S_c+ABX*I_ERI_Fx2z_S_Py_S_c;
  Double I_ERI_F3y_Px_Py_S_c = I_ERI_Gx3y_S_Py_S_c+ABX*I_ERI_F3y_S_Py_S_c;
  Double I_ERI_F2yz_Px_Py_S_c = I_ERI_Gx2yz_S_Py_S_c+ABX*I_ERI_F2yz_S_Py_S_c;
  Double I_ERI_Fy2z_Px_Py_S_c = I_ERI_Gxy2z_S_Py_S_c+ABX*I_ERI_Fy2z_S_Py_S_c;
  Double I_ERI_F3z_Px_Py_S_c = I_ERI_Gx3z_S_Py_S_c+ABX*I_ERI_F3z_S_Py_S_c;
  Double I_ERI_F3x_Py_Py_S_c = I_ERI_G3xy_S_Py_S_c+ABY*I_ERI_F3x_S_Py_S_c;
  Double I_ERI_F2xy_Py_Py_S_c = I_ERI_G2x2y_S_Py_S_c+ABY*I_ERI_F2xy_S_Py_S_c;
  Double I_ERI_F2xz_Py_Py_S_c = I_ERI_G2xyz_S_Py_S_c+ABY*I_ERI_F2xz_S_Py_S_c;
  Double I_ERI_Fx2y_Py_Py_S_c = I_ERI_Gx3y_S_Py_S_c+ABY*I_ERI_Fx2y_S_Py_S_c;
  Double I_ERI_Fxyz_Py_Py_S_c = I_ERI_Gx2yz_S_Py_S_c+ABY*I_ERI_Fxyz_S_Py_S_c;
  Double I_ERI_Fx2z_Py_Py_S_c = I_ERI_Gxy2z_S_Py_S_c+ABY*I_ERI_Fx2z_S_Py_S_c;
  Double I_ERI_F3y_Py_Py_S_c = I_ERI_G4y_S_Py_S_c+ABY*I_ERI_F3y_S_Py_S_c;
  Double I_ERI_F2yz_Py_Py_S_c = I_ERI_G3yz_S_Py_S_c+ABY*I_ERI_F2yz_S_Py_S_c;
  Double I_ERI_Fy2z_Py_Py_S_c = I_ERI_G2y2z_S_Py_S_c+ABY*I_ERI_Fy2z_S_Py_S_c;
  Double I_ERI_F3z_Py_Py_S_c = I_ERI_Gy3z_S_Py_S_c+ABY*I_ERI_F3z_S_Py_S_c;
  Double I_ERI_F3x_Pz_Py_S_c = I_ERI_G3xz_S_Py_S_c+ABZ*I_ERI_F3x_S_Py_S_c;
  Double I_ERI_F2xy_Pz_Py_S_c = I_ERI_G2xyz_S_Py_S_c+ABZ*I_ERI_F2xy_S_Py_S_c;
  Double I_ERI_F2xz_Pz_Py_S_c = I_ERI_G2x2z_S_Py_S_c+ABZ*I_ERI_F2xz_S_Py_S_c;
  Double I_ERI_Fx2y_Pz_Py_S_c = I_ERI_Gx2yz_S_Py_S_c+ABZ*I_ERI_Fx2y_S_Py_S_c;
  Double I_ERI_Fxyz_Pz_Py_S_c = I_ERI_Gxy2z_S_Py_S_c+ABZ*I_ERI_Fxyz_S_Py_S_c;
  Double I_ERI_Fx2z_Pz_Py_S_c = I_ERI_Gx3z_S_Py_S_c+ABZ*I_ERI_Fx2z_S_Py_S_c;
  Double I_ERI_F3y_Pz_Py_S_c = I_ERI_G3yz_S_Py_S_c+ABZ*I_ERI_F3y_S_Py_S_c;
  Double I_ERI_F2yz_Pz_Py_S_c = I_ERI_G2y2z_S_Py_S_c+ABZ*I_ERI_F2yz_S_Py_S_c;
  Double I_ERI_Fy2z_Pz_Py_S_c = I_ERI_Gy3z_S_Py_S_c+ABZ*I_ERI_Fy2z_S_Py_S_c;
  Double I_ERI_F3z_Pz_Py_S_c = I_ERI_G4z_S_Py_S_c+ABZ*I_ERI_F3z_S_Py_S_c;
  Double I_ERI_F3x_Px_Pz_S_c = I_ERI_G4x_S_Pz_S_c+ABX*I_ERI_F3x_S_Pz_S_c;
  Double I_ERI_F2xy_Px_Pz_S_c = I_ERI_G3xy_S_Pz_S_c+ABX*I_ERI_F2xy_S_Pz_S_c;
  Double I_ERI_F2xz_Px_Pz_S_c = I_ERI_G3xz_S_Pz_S_c+ABX*I_ERI_F2xz_S_Pz_S_c;
  Double I_ERI_Fx2y_Px_Pz_S_c = I_ERI_G2x2y_S_Pz_S_c+ABX*I_ERI_Fx2y_S_Pz_S_c;
  Double I_ERI_Fxyz_Px_Pz_S_c = I_ERI_G2xyz_S_Pz_S_c+ABX*I_ERI_Fxyz_S_Pz_S_c;
  Double I_ERI_Fx2z_Px_Pz_S_c = I_ERI_G2x2z_S_Pz_S_c+ABX*I_ERI_Fx2z_S_Pz_S_c;
  Double I_ERI_F3y_Px_Pz_S_c = I_ERI_Gx3y_S_Pz_S_c+ABX*I_ERI_F3y_S_Pz_S_c;
  Double I_ERI_F2yz_Px_Pz_S_c = I_ERI_Gx2yz_S_Pz_S_c+ABX*I_ERI_F2yz_S_Pz_S_c;
  Double I_ERI_Fy2z_Px_Pz_S_c = I_ERI_Gxy2z_S_Pz_S_c+ABX*I_ERI_Fy2z_S_Pz_S_c;
  Double I_ERI_F3z_Px_Pz_S_c = I_ERI_Gx3z_S_Pz_S_c+ABX*I_ERI_F3z_S_Pz_S_c;
  Double I_ERI_F3x_Py_Pz_S_c = I_ERI_G3xy_S_Pz_S_c+ABY*I_ERI_F3x_S_Pz_S_c;
  Double I_ERI_F2xy_Py_Pz_S_c = I_ERI_G2x2y_S_Pz_S_c+ABY*I_ERI_F2xy_S_Pz_S_c;
  Double I_ERI_F2xz_Py_Pz_S_c = I_ERI_G2xyz_S_Pz_S_c+ABY*I_ERI_F2xz_S_Pz_S_c;
  Double I_ERI_Fx2y_Py_Pz_S_c = I_ERI_Gx3y_S_Pz_S_c+ABY*I_ERI_Fx2y_S_Pz_S_c;
  Double I_ERI_Fxyz_Py_Pz_S_c = I_ERI_Gx2yz_S_Pz_S_c+ABY*I_ERI_Fxyz_S_Pz_S_c;
  Double I_ERI_Fx2z_Py_Pz_S_c = I_ERI_Gxy2z_S_Pz_S_c+ABY*I_ERI_Fx2z_S_Pz_S_c;
  Double I_ERI_F3y_Py_Pz_S_c = I_ERI_G4y_S_Pz_S_c+ABY*I_ERI_F3y_S_Pz_S_c;
  Double I_ERI_F2yz_Py_Pz_S_c = I_ERI_G3yz_S_Pz_S_c+ABY*I_ERI_F2yz_S_Pz_S_c;
  Double I_ERI_Fy2z_Py_Pz_S_c = I_ERI_G2y2z_S_Pz_S_c+ABY*I_ERI_Fy2z_S_Pz_S_c;
  Double I_ERI_F3z_Py_Pz_S_c = I_ERI_Gy3z_S_Pz_S_c+ABY*I_ERI_F3z_S_Pz_S_c;
  Double I_ERI_F3x_Pz_Pz_S_c = I_ERI_G3xz_S_Pz_S_c+ABZ*I_ERI_F3x_S_Pz_S_c;
  Double I_ERI_F2xy_Pz_Pz_S_c = I_ERI_G2xyz_S_Pz_S_c+ABZ*I_ERI_F2xy_S_Pz_S_c;
  Double I_ERI_F2xz_Pz_Pz_S_c = I_ERI_G2x2z_S_Pz_S_c+ABZ*I_ERI_F2xz_S_Pz_S_c;
  Double I_ERI_Fx2y_Pz_Pz_S_c = I_ERI_Gx2yz_S_Pz_S_c+ABZ*I_ERI_Fx2y_S_Pz_S_c;
  Double I_ERI_Fxyz_Pz_Pz_S_c = I_ERI_Gxy2z_S_Pz_S_c+ABZ*I_ERI_Fxyz_S_Pz_S_c;
  Double I_ERI_Fx2z_Pz_Pz_S_c = I_ERI_Gx3z_S_Pz_S_c+ABZ*I_ERI_Fx2z_S_Pz_S_c;
  Double I_ERI_F3y_Pz_Pz_S_c = I_ERI_G3yz_S_Pz_S_c+ABZ*I_ERI_F3y_S_Pz_S_c;
  Double I_ERI_F2yz_Pz_Pz_S_c = I_ERI_G2y2z_S_Pz_S_c+ABZ*I_ERI_F2yz_S_Pz_S_c;
  Double I_ERI_Fy2z_Pz_Pz_S_c = I_ERI_Gy3z_S_Pz_S_c+ABZ*I_ERI_Fy2z_S_Pz_S_c;
  Double I_ERI_F3z_Pz_Pz_S_c = I_ERI_G4z_S_Pz_S_c+ABZ*I_ERI_F3z_S_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_P_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   * RHS shell quartet name: SQ_ERI_D_P_P_S_c
   ************************************************************/
  Double I_ERI_D2x_D2x_Px_S_c = I_ERI_F3x_Px_Px_S_c+ABX*I_ERI_D2x_Px_Px_S_c;
  Double I_ERI_Dxy_D2x_Px_S_c = I_ERI_F2xy_Px_Px_S_c+ABX*I_ERI_Dxy_Px_Px_S_c;
  Double I_ERI_Dxz_D2x_Px_S_c = I_ERI_F2xz_Px_Px_S_c+ABX*I_ERI_Dxz_Px_Px_S_c;
  Double I_ERI_D2y_D2x_Px_S_c = I_ERI_Fx2y_Px_Px_S_c+ABX*I_ERI_D2y_Px_Px_S_c;
  Double I_ERI_Dyz_D2x_Px_S_c = I_ERI_Fxyz_Px_Px_S_c+ABX*I_ERI_Dyz_Px_Px_S_c;
  Double I_ERI_D2z_D2x_Px_S_c = I_ERI_Fx2z_Px_Px_S_c+ABX*I_ERI_D2z_Px_Px_S_c;
  Double I_ERI_D2x_Dxy_Px_S_c = I_ERI_F2xy_Px_Px_S_c+ABY*I_ERI_D2x_Px_Px_S_c;
  Double I_ERI_Dxy_Dxy_Px_S_c = I_ERI_Fx2y_Px_Px_S_c+ABY*I_ERI_Dxy_Px_Px_S_c;
  Double I_ERI_Dxz_Dxy_Px_S_c = I_ERI_Fxyz_Px_Px_S_c+ABY*I_ERI_Dxz_Px_Px_S_c;
  Double I_ERI_D2y_Dxy_Px_S_c = I_ERI_F3y_Px_Px_S_c+ABY*I_ERI_D2y_Px_Px_S_c;
  Double I_ERI_Dyz_Dxy_Px_S_c = I_ERI_F2yz_Px_Px_S_c+ABY*I_ERI_Dyz_Px_Px_S_c;
  Double I_ERI_D2z_Dxy_Px_S_c = I_ERI_Fy2z_Px_Px_S_c+ABY*I_ERI_D2z_Px_Px_S_c;
  Double I_ERI_D2x_Dxz_Px_S_c = I_ERI_F2xz_Px_Px_S_c+ABZ*I_ERI_D2x_Px_Px_S_c;
  Double I_ERI_Dxy_Dxz_Px_S_c = I_ERI_Fxyz_Px_Px_S_c+ABZ*I_ERI_Dxy_Px_Px_S_c;
  Double I_ERI_Dxz_Dxz_Px_S_c = I_ERI_Fx2z_Px_Px_S_c+ABZ*I_ERI_Dxz_Px_Px_S_c;
  Double I_ERI_D2y_Dxz_Px_S_c = I_ERI_F2yz_Px_Px_S_c+ABZ*I_ERI_D2y_Px_Px_S_c;
  Double I_ERI_Dyz_Dxz_Px_S_c = I_ERI_Fy2z_Px_Px_S_c+ABZ*I_ERI_Dyz_Px_Px_S_c;
  Double I_ERI_D2z_Dxz_Px_S_c = I_ERI_F3z_Px_Px_S_c+ABZ*I_ERI_D2z_Px_Px_S_c;
  Double I_ERI_D2x_D2y_Px_S_c = I_ERI_F2xy_Py_Px_S_c+ABY*I_ERI_D2x_Py_Px_S_c;
  Double I_ERI_Dxy_D2y_Px_S_c = I_ERI_Fx2y_Py_Px_S_c+ABY*I_ERI_Dxy_Py_Px_S_c;
  Double I_ERI_Dxz_D2y_Px_S_c = I_ERI_Fxyz_Py_Px_S_c+ABY*I_ERI_Dxz_Py_Px_S_c;
  Double I_ERI_D2y_D2y_Px_S_c = I_ERI_F3y_Py_Px_S_c+ABY*I_ERI_D2y_Py_Px_S_c;
  Double I_ERI_Dyz_D2y_Px_S_c = I_ERI_F2yz_Py_Px_S_c+ABY*I_ERI_Dyz_Py_Px_S_c;
  Double I_ERI_D2z_D2y_Px_S_c = I_ERI_Fy2z_Py_Px_S_c+ABY*I_ERI_D2z_Py_Px_S_c;
  Double I_ERI_D2x_Dyz_Px_S_c = I_ERI_F2xz_Py_Px_S_c+ABZ*I_ERI_D2x_Py_Px_S_c;
  Double I_ERI_Dxy_Dyz_Px_S_c = I_ERI_Fxyz_Py_Px_S_c+ABZ*I_ERI_Dxy_Py_Px_S_c;
  Double I_ERI_Dxz_Dyz_Px_S_c = I_ERI_Fx2z_Py_Px_S_c+ABZ*I_ERI_Dxz_Py_Px_S_c;
  Double I_ERI_D2y_Dyz_Px_S_c = I_ERI_F2yz_Py_Px_S_c+ABZ*I_ERI_D2y_Py_Px_S_c;
  Double I_ERI_Dyz_Dyz_Px_S_c = I_ERI_Fy2z_Py_Px_S_c+ABZ*I_ERI_Dyz_Py_Px_S_c;
  Double I_ERI_D2z_Dyz_Px_S_c = I_ERI_F3z_Py_Px_S_c+ABZ*I_ERI_D2z_Py_Px_S_c;
  Double I_ERI_D2x_D2z_Px_S_c = I_ERI_F2xz_Pz_Px_S_c+ABZ*I_ERI_D2x_Pz_Px_S_c;
  Double I_ERI_Dxy_D2z_Px_S_c = I_ERI_Fxyz_Pz_Px_S_c+ABZ*I_ERI_Dxy_Pz_Px_S_c;
  Double I_ERI_Dxz_D2z_Px_S_c = I_ERI_Fx2z_Pz_Px_S_c+ABZ*I_ERI_Dxz_Pz_Px_S_c;
  Double I_ERI_D2y_D2z_Px_S_c = I_ERI_F2yz_Pz_Px_S_c+ABZ*I_ERI_D2y_Pz_Px_S_c;
  Double I_ERI_Dyz_D2z_Px_S_c = I_ERI_Fy2z_Pz_Px_S_c+ABZ*I_ERI_Dyz_Pz_Px_S_c;
  Double I_ERI_D2z_D2z_Px_S_c = I_ERI_F3z_Pz_Px_S_c+ABZ*I_ERI_D2z_Pz_Px_S_c;
  Double I_ERI_D2x_D2x_Py_S_c = I_ERI_F3x_Px_Py_S_c+ABX*I_ERI_D2x_Px_Py_S_c;
  Double I_ERI_Dxy_D2x_Py_S_c = I_ERI_F2xy_Px_Py_S_c+ABX*I_ERI_Dxy_Px_Py_S_c;
  Double I_ERI_Dxz_D2x_Py_S_c = I_ERI_F2xz_Px_Py_S_c+ABX*I_ERI_Dxz_Px_Py_S_c;
  Double I_ERI_D2y_D2x_Py_S_c = I_ERI_Fx2y_Px_Py_S_c+ABX*I_ERI_D2y_Px_Py_S_c;
  Double I_ERI_Dyz_D2x_Py_S_c = I_ERI_Fxyz_Px_Py_S_c+ABX*I_ERI_Dyz_Px_Py_S_c;
  Double I_ERI_D2z_D2x_Py_S_c = I_ERI_Fx2z_Px_Py_S_c+ABX*I_ERI_D2z_Px_Py_S_c;
  Double I_ERI_D2x_Dxy_Py_S_c = I_ERI_F2xy_Px_Py_S_c+ABY*I_ERI_D2x_Px_Py_S_c;
  Double I_ERI_Dxy_Dxy_Py_S_c = I_ERI_Fx2y_Px_Py_S_c+ABY*I_ERI_Dxy_Px_Py_S_c;
  Double I_ERI_Dxz_Dxy_Py_S_c = I_ERI_Fxyz_Px_Py_S_c+ABY*I_ERI_Dxz_Px_Py_S_c;
  Double I_ERI_D2y_Dxy_Py_S_c = I_ERI_F3y_Px_Py_S_c+ABY*I_ERI_D2y_Px_Py_S_c;
  Double I_ERI_Dyz_Dxy_Py_S_c = I_ERI_F2yz_Px_Py_S_c+ABY*I_ERI_Dyz_Px_Py_S_c;
  Double I_ERI_D2z_Dxy_Py_S_c = I_ERI_Fy2z_Px_Py_S_c+ABY*I_ERI_D2z_Px_Py_S_c;
  Double I_ERI_D2x_Dxz_Py_S_c = I_ERI_F2xz_Px_Py_S_c+ABZ*I_ERI_D2x_Px_Py_S_c;
  Double I_ERI_Dxy_Dxz_Py_S_c = I_ERI_Fxyz_Px_Py_S_c+ABZ*I_ERI_Dxy_Px_Py_S_c;
  Double I_ERI_Dxz_Dxz_Py_S_c = I_ERI_Fx2z_Px_Py_S_c+ABZ*I_ERI_Dxz_Px_Py_S_c;
  Double I_ERI_D2y_Dxz_Py_S_c = I_ERI_F2yz_Px_Py_S_c+ABZ*I_ERI_D2y_Px_Py_S_c;
  Double I_ERI_Dyz_Dxz_Py_S_c = I_ERI_Fy2z_Px_Py_S_c+ABZ*I_ERI_Dyz_Px_Py_S_c;
  Double I_ERI_D2z_Dxz_Py_S_c = I_ERI_F3z_Px_Py_S_c+ABZ*I_ERI_D2z_Px_Py_S_c;
  Double I_ERI_D2x_D2y_Py_S_c = I_ERI_F2xy_Py_Py_S_c+ABY*I_ERI_D2x_Py_Py_S_c;
  Double I_ERI_Dxy_D2y_Py_S_c = I_ERI_Fx2y_Py_Py_S_c+ABY*I_ERI_Dxy_Py_Py_S_c;
  Double I_ERI_Dxz_D2y_Py_S_c = I_ERI_Fxyz_Py_Py_S_c+ABY*I_ERI_Dxz_Py_Py_S_c;
  Double I_ERI_D2y_D2y_Py_S_c = I_ERI_F3y_Py_Py_S_c+ABY*I_ERI_D2y_Py_Py_S_c;
  Double I_ERI_Dyz_D2y_Py_S_c = I_ERI_F2yz_Py_Py_S_c+ABY*I_ERI_Dyz_Py_Py_S_c;
  Double I_ERI_D2z_D2y_Py_S_c = I_ERI_Fy2z_Py_Py_S_c+ABY*I_ERI_D2z_Py_Py_S_c;
  Double I_ERI_D2x_Dyz_Py_S_c = I_ERI_F2xz_Py_Py_S_c+ABZ*I_ERI_D2x_Py_Py_S_c;
  Double I_ERI_Dxy_Dyz_Py_S_c = I_ERI_Fxyz_Py_Py_S_c+ABZ*I_ERI_Dxy_Py_Py_S_c;
  Double I_ERI_Dxz_Dyz_Py_S_c = I_ERI_Fx2z_Py_Py_S_c+ABZ*I_ERI_Dxz_Py_Py_S_c;
  Double I_ERI_D2y_Dyz_Py_S_c = I_ERI_F2yz_Py_Py_S_c+ABZ*I_ERI_D2y_Py_Py_S_c;
  Double I_ERI_Dyz_Dyz_Py_S_c = I_ERI_Fy2z_Py_Py_S_c+ABZ*I_ERI_Dyz_Py_Py_S_c;
  Double I_ERI_D2z_Dyz_Py_S_c = I_ERI_F3z_Py_Py_S_c+ABZ*I_ERI_D2z_Py_Py_S_c;
  Double I_ERI_D2x_D2z_Py_S_c = I_ERI_F2xz_Pz_Py_S_c+ABZ*I_ERI_D2x_Pz_Py_S_c;
  Double I_ERI_Dxy_D2z_Py_S_c = I_ERI_Fxyz_Pz_Py_S_c+ABZ*I_ERI_Dxy_Pz_Py_S_c;
  Double I_ERI_Dxz_D2z_Py_S_c = I_ERI_Fx2z_Pz_Py_S_c+ABZ*I_ERI_Dxz_Pz_Py_S_c;
  Double I_ERI_D2y_D2z_Py_S_c = I_ERI_F2yz_Pz_Py_S_c+ABZ*I_ERI_D2y_Pz_Py_S_c;
  Double I_ERI_Dyz_D2z_Py_S_c = I_ERI_Fy2z_Pz_Py_S_c+ABZ*I_ERI_Dyz_Pz_Py_S_c;
  Double I_ERI_D2z_D2z_Py_S_c = I_ERI_F3z_Pz_Py_S_c+ABZ*I_ERI_D2z_Pz_Py_S_c;
  Double I_ERI_D2x_D2x_Pz_S_c = I_ERI_F3x_Px_Pz_S_c+ABX*I_ERI_D2x_Px_Pz_S_c;
  Double I_ERI_Dxy_D2x_Pz_S_c = I_ERI_F2xy_Px_Pz_S_c+ABX*I_ERI_Dxy_Px_Pz_S_c;
  Double I_ERI_Dxz_D2x_Pz_S_c = I_ERI_F2xz_Px_Pz_S_c+ABX*I_ERI_Dxz_Px_Pz_S_c;
  Double I_ERI_D2y_D2x_Pz_S_c = I_ERI_Fx2y_Px_Pz_S_c+ABX*I_ERI_D2y_Px_Pz_S_c;
  Double I_ERI_Dyz_D2x_Pz_S_c = I_ERI_Fxyz_Px_Pz_S_c+ABX*I_ERI_Dyz_Px_Pz_S_c;
  Double I_ERI_D2z_D2x_Pz_S_c = I_ERI_Fx2z_Px_Pz_S_c+ABX*I_ERI_D2z_Px_Pz_S_c;
  Double I_ERI_D2x_Dxy_Pz_S_c = I_ERI_F2xy_Px_Pz_S_c+ABY*I_ERI_D2x_Px_Pz_S_c;
  Double I_ERI_Dxy_Dxy_Pz_S_c = I_ERI_Fx2y_Px_Pz_S_c+ABY*I_ERI_Dxy_Px_Pz_S_c;
  Double I_ERI_Dxz_Dxy_Pz_S_c = I_ERI_Fxyz_Px_Pz_S_c+ABY*I_ERI_Dxz_Px_Pz_S_c;
  Double I_ERI_D2y_Dxy_Pz_S_c = I_ERI_F3y_Px_Pz_S_c+ABY*I_ERI_D2y_Px_Pz_S_c;
  Double I_ERI_Dyz_Dxy_Pz_S_c = I_ERI_F2yz_Px_Pz_S_c+ABY*I_ERI_Dyz_Px_Pz_S_c;
  Double I_ERI_D2z_Dxy_Pz_S_c = I_ERI_Fy2z_Px_Pz_S_c+ABY*I_ERI_D2z_Px_Pz_S_c;
  Double I_ERI_D2x_Dxz_Pz_S_c = I_ERI_F2xz_Px_Pz_S_c+ABZ*I_ERI_D2x_Px_Pz_S_c;
  Double I_ERI_Dxy_Dxz_Pz_S_c = I_ERI_Fxyz_Px_Pz_S_c+ABZ*I_ERI_Dxy_Px_Pz_S_c;
  Double I_ERI_Dxz_Dxz_Pz_S_c = I_ERI_Fx2z_Px_Pz_S_c+ABZ*I_ERI_Dxz_Px_Pz_S_c;
  Double I_ERI_D2y_Dxz_Pz_S_c = I_ERI_F2yz_Px_Pz_S_c+ABZ*I_ERI_D2y_Px_Pz_S_c;
  Double I_ERI_Dyz_Dxz_Pz_S_c = I_ERI_Fy2z_Px_Pz_S_c+ABZ*I_ERI_Dyz_Px_Pz_S_c;
  Double I_ERI_D2z_Dxz_Pz_S_c = I_ERI_F3z_Px_Pz_S_c+ABZ*I_ERI_D2z_Px_Pz_S_c;
  Double I_ERI_D2x_D2y_Pz_S_c = I_ERI_F2xy_Py_Pz_S_c+ABY*I_ERI_D2x_Py_Pz_S_c;
  Double I_ERI_Dxy_D2y_Pz_S_c = I_ERI_Fx2y_Py_Pz_S_c+ABY*I_ERI_Dxy_Py_Pz_S_c;
  Double I_ERI_Dxz_D2y_Pz_S_c = I_ERI_Fxyz_Py_Pz_S_c+ABY*I_ERI_Dxz_Py_Pz_S_c;
  Double I_ERI_D2y_D2y_Pz_S_c = I_ERI_F3y_Py_Pz_S_c+ABY*I_ERI_D2y_Py_Pz_S_c;
  Double I_ERI_Dyz_D2y_Pz_S_c = I_ERI_F2yz_Py_Pz_S_c+ABY*I_ERI_Dyz_Py_Pz_S_c;
  Double I_ERI_D2z_D2y_Pz_S_c = I_ERI_Fy2z_Py_Pz_S_c+ABY*I_ERI_D2z_Py_Pz_S_c;
  Double I_ERI_D2x_Dyz_Pz_S_c = I_ERI_F2xz_Py_Pz_S_c+ABZ*I_ERI_D2x_Py_Pz_S_c;
  Double I_ERI_Dxy_Dyz_Pz_S_c = I_ERI_Fxyz_Py_Pz_S_c+ABZ*I_ERI_Dxy_Py_Pz_S_c;
  Double I_ERI_Dxz_Dyz_Pz_S_c = I_ERI_Fx2z_Py_Pz_S_c+ABZ*I_ERI_Dxz_Py_Pz_S_c;
  Double I_ERI_D2y_Dyz_Pz_S_c = I_ERI_F2yz_Py_Pz_S_c+ABZ*I_ERI_D2y_Py_Pz_S_c;
  Double I_ERI_Dyz_Dyz_Pz_S_c = I_ERI_Fy2z_Py_Pz_S_c+ABZ*I_ERI_Dyz_Py_Pz_S_c;
  Double I_ERI_D2z_Dyz_Pz_S_c = I_ERI_F3z_Py_Pz_S_c+ABZ*I_ERI_D2z_Py_Pz_S_c;
  Double I_ERI_D2x_D2z_Pz_S_c = I_ERI_F2xz_Pz_Pz_S_c+ABZ*I_ERI_D2x_Pz_Pz_S_c;
  Double I_ERI_Dxy_D2z_Pz_S_c = I_ERI_Fxyz_Pz_Pz_S_c+ABZ*I_ERI_Dxy_Pz_Pz_S_c;
  Double I_ERI_Dxz_D2z_Pz_S_c = I_ERI_Fx2z_Pz_Pz_S_c+ABZ*I_ERI_Dxz_Pz_Pz_S_c;
  Double I_ERI_D2y_D2z_Pz_S_c = I_ERI_F2yz_Pz_Pz_S_c+ABZ*I_ERI_D2y_Pz_Pz_S_c;
  Double I_ERI_Dyz_D2z_Pz_S_c = I_ERI_Fy2z_Pz_Pz_S_c+ABZ*I_ERI_Dyz_Pz_Pz_S_c;
  Double I_ERI_D2z_D2z_Pz_S_c = I_ERI_F3z_Pz_Pz_S_c+ABZ*I_ERI_D2z_Pz_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_S_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_S_S_c
   * RHS shell quartet name: SQ_ERI_G_S_S_S_c
   ************************************************************/
  Double I_ERI_G4x_Px_S_S_c = I_ERI_H5x_S_S_S_c+ABX*I_ERI_G4x_S_S_S_c;
  Double I_ERI_G3xy_Px_S_S_c = I_ERI_H4xy_S_S_S_c+ABX*I_ERI_G3xy_S_S_S_c;
  Double I_ERI_G3xz_Px_S_S_c = I_ERI_H4xz_S_S_S_c+ABX*I_ERI_G3xz_S_S_S_c;
  Double I_ERI_G2x2y_Px_S_S_c = I_ERI_H3x2y_S_S_S_c+ABX*I_ERI_G2x2y_S_S_S_c;
  Double I_ERI_G2xyz_Px_S_S_c = I_ERI_H3xyz_S_S_S_c+ABX*I_ERI_G2xyz_S_S_S_c;
  Double I_ERI_G2x2z_Px_S_S_c = I_ERI_H3x2z_S_S_S_c+ABX*I_ERI_G2x2z_S_S_S_c;
  Double I_ERI_Gx3y_Px_S_S_c = I_ERI_H2x3y_S_S_S_c+ABX*I_ERI_Gx3y_S_S_S_c;
  Double I_ERI_Gx2yz_Px_S_S_c = I_ERI_H2x2yz_S_S_S_c+ABX*I_ERI_Gx2yz_S_S_S_c;
  Double I_ERI_Gxy2z_Px_S_S_c = I_ERI_H2xy2z_S_S_S_c+ABX*I_ERI_Gxy2z_S_S_S_c;
  Double I_ERI_Gx3z_Px_S_S_c = I_ERI_H2x3z_S_S_S_c+ABX*I_ERI_Gx3z_S_S_S_c;
  Double I_ERI_G4y_Px_S_S_c = I_ERI_Hx4y_S_S_S_c+ABX*I_ERI_G4y_S_S_S_c;
  Double I_ERI_G3yz_Px_S_S_c = I_ERI_Hx3yz_S_S_S_c+ABX*I_ERI_G3yz_S_S_S_c;
  Double I_ERI_G2y2z_Px_S_S_c = I_ERI_Hx2y2z_S_S_S_c+ABX*I_ERI_G2y2z_S_S_S_c;
  Double I_ERI_Gy3z_Px_S_S_c = I_ERI_Hxy3z_S_S_S_c+ABX*I_ERI_Gy3z_S_S_S_c;
  Double I_ERI_G4z_Px_S_S_c = I_ERI_Hx4z_S_S_S_c+ABX*I_ERI_G4z_S_S_S_c;
  Double I_ERI_G3xy_Py_S_S_c = I_ERI_H3x2y_S_S_S_c+ABY*I_ERI_G3xy_S_S_S_c;
  Double I_ERI_G3xz_Py_S_S_c = I_ERI_H3xyz_S_S_S_c+ABY*I_ERI_G3xz_S_S_S_c;
  Double I_ERI_G2x2y_Py_S_S_c = I_ERI_H2x3y_S_S_S_c+ABY*I_ERI_G2x2y_S_S_S_c;
  Double I_ERI_G2xyz_Py_S_S_c = I_ERI_H2x2yz_S_S_S_c+ABY*I_ERI_G2xyz_S_S_S_c;
  Double I_ERI_G2x2z_Py_S_S_c = I_ERI_H2xy2z_S_S_S_c+ABY*I_ERI_G2x2z_S_S_S_c;
  Double I_ERI_Gx3y_Py_S_S_c = I_ERI_Hx4y_S_S_S_c+ABY*I_ERI_Gx3y_S_S_S_c;
  Double I_ERI_Gx2yz_Py_S_S_c = I_ERI_Hx3yz_S_S_S_c+ABY*I_ERI_Gx2yz_S_S_S_c;
  Double I_ERI_Gxy2z_Py_S_S_c = I_ERI_Hx2y2z_S_S_S_c+ABY*I_ERI_Gxy2z_S_S_S_c;
  Double I_ERI_Gx3z_Py_S_S_c = I_ERI_Hxy3z_S_S_S_c+ABY*I_ERI_Gx3z_S_S_S_c;
  Double I_ERI_G4y_Py_S_S_c = I_ERI_H5y_S_S_S_c+ABY*I_ERI_G4y_S_S_S_c;
  Double I_ERI_G3yz_Py_S_S_c = I_ERI_H4yz_S_S_S_c+ABY*I_ERI_G3yz_S_S_S_c;
  Double I_ERI_G2y2z_Py_S_S_c = I_ERI_H3y2z_S_S_S_c+ABY*I_ERI_G2y2z_S_S_S_c;
  Double I_ERI_Gy3z_Py_S_S_c = I_ERI_H2y3z_S_S_S_c+ABY*I_ERI_Gy3z_S_S_S_c;
  Double I_ERI_G4z_Py_S_S_c = I_ERI_Hy4z_S_S_S_c+ABY*I_ERI_G4z_S_S_S_c;
  Double I_ERI_G3xz_Pz_S_S_c = I_ERI_H3x2z_S_S_S_c+ABZ*I_ERI_G3xz_S_S_S_c;
  Double I_ERI_G2xyz_Pz_S_S_c = I_ERI_H2xy2z_S_S_S_c+ABZ*I_ERI_G2xyz_S_S_S_c;
  Double I_ERI_G2x2z_Pz_S_S_c = I_ERI_H2x3z_S_S_S_c+ABZ*I_ERI_G2x2z_S_S_S_c;
  Double I_ERI_Gx2yz_Pz_S_S_c = I_ERI_Hx2y2z_S_S_S_c+ABZ*I_ERI_Gx2yz_S_S_S_c;
  Double I_ERI_Gxy2z_Pz_S_S_c = I_ERI_Hxy3z_S_S_S_c+ABZ*I_ERI_Gxy2z_S_S_S_c;
  Double I_ERI_Gx3z_Pz_S_S_c = I_ERI_Hx4z_S_S_S_c+ABZ*I_ERI_Gx3z_S_S_S_c;
  Double I_ERI_G3yz_Pz_S_S_c = I_ERI_H3y2z_S_S_S_c+ABZ*I_ERI_G3yz_S_S_S_c;
  Double I_ERI_G2y2z_Pz_S_S_c = I_ERI_H2y3z_S_S_S_c+ABZ*I_ERI_G2y2z_S_S_S_c;
  Double I_ERI_Gy3z_Pz_S_S_c = I_ERI_Hy4z_S_S_S_c+ABZ*I_ERI_Gy3z_S_S_S_c;
  Double I_ERI_G4z_Pz_S_S_c = I_ERI_H5z_S_S_S_c+ABZ*I_ERI_G4z_S_S_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_c
   * RHS shell quartet name: SQ_ERI_F_P_S_S_c
   ************************************************************/
  Double I_ERI_F3x_D2x_S_S_c = I_ERI_G4x_Px_S_S_c+ABX*I_ERI_F3x_Px_S_S_c;
  Double I_ERI_F2xy_D2x_S_S_c = I_ERI_G3xy_Px_S_S_c+ABX*I_ERI_F2xy_Px_S_S_c;
  Double I_ERI_F2xz_D2x_S_S_c = I_ERI_G3xz_Px_S_S_c+ABX*I_ERI_F2xz_Px_S_S_c;
  Double I_ERI_Fx2y_D2x_S_S_c = I_ERI_G2x2y_Px_S_S_c+ABX*I_ERI_Fx2y_Px_S_S_c;
  Double I_ERI_Fxyz_D2x_S_S_c = I_ERI_G2xyz_Px_S_S_c+ABX*I_ERI_Fxyz_Px_S_S_c;
  Double I_ERI_Fx2z_D2x_S_S_c = I_ERI_G2x2z_Px_S_S_c+ABX*I_ERI_Fx2z_Px_S_S_c;
  Double I_ERI_F3y_D2x_S_S_c = I_ERI_Gx3y_Px_S_S_c+ABX*I_ERI_F3y_Px_S_S_c;
  Double I_ERI_F2yz_D2x_S_S_c = I_ERI_Gx2yz_Px_S_S_c+ABX*I_ERI_F2yz_Px_S_S_c;
  Double I_ERI_Fy2z_D2x_S_S_c = I_ERI_Gxy2z_Px_S_S_c+ABX*I_ERI_Fy2z_Px_S_S_c;
  Double I_ERI_F3z_D2x_S_S_c = I_ERI_Gx3z_Px_S_S_c+ABX*I_ERI_F3z_Px_S_S_c;
  Double I_ERI_F3x_Dxy_S_S_c = I_ERI_G3xy_Px_S_S_c+ABY*I_ERI_F3x_Px_S_S_c;
  Double I_ERI_F2xy_Dxy_S_S_c = I_ERI_G2x2y_Px_S_S_c+ABY*I_ERI_F2xy_Px_S_S_c;
  Double I_ERI_F2xz_Dxy_S_S_c = I_ERI_G2xyz_Px_S_S_c+ABY*I_ERI_F2xz_Px_S_S_c;
  Double I_ERI_Fx2y_Dxy_S_S_c = I_ERI_Gx3y_Px_S_S_c+ABY*I_ERI_Fx2y_Px_S_S_c;
  Double I_ERI_Fxyz_Dxy_S_S_c = I_ERI_Gx2yz_Px_S_S_c+ABY*I_ERI_Fxyz_Px_S_S_c;
  Double I_ERI_Fx2z_Dxy_S_S_c = I_ERI_Gxy2z_Px_S_S_c+ABY*I_ERI_Fx2z_Px_S_S_c;
  Double I_ERI_F3y_Dxy_S_S_c = I_ERI_G4y_Px_S_S_c+ABY*I_ERI_F3y_Px_S_S_c;
  Double I_ERI_F2yz_Dxy_S_S_c = I_ERI_G3yz_Px_S_S_c+ABY*I_ERI_F2yz_Px_S_S_c;
  Double I_ERI_Fy2z_Dxy_S_S_c = I_ERI_G2y2z_Px_S_S_c+ABY*I_ERI_Fy2z_Px_S_S_c;
  Double I_ERI_F3z_Dxy_S_S_c = I_ERI_Gy3z_Px_S_S_c+ABY*I_ERI_F3z_Px_S_S_c;
  Double I_ERI_F3x_Dxz_S_S_c = I_ERI_G3xz_Px_S_S_c+ABZ*I_ERI_F3x_Px_S_S_c;
  Double I_ERI_F2xy_Dxz_S_S_c = I_ERI_G2xyz_Px_S_S_c+ABZ*I_ERI_F2xy_Px_S_S_c;
  Double I_ERI_F2xz_Dxz_S_S_c = I_ERI_G2x2z_Px_S_S_c+ABZ*I_ERI_F2xz_Px_S_S_c;
  Double I_ERI_Fx2y_Dxz_S_S_c = I_ERI_Gx2yz_Px_S_S_c+ABZ*I_ERI_Fx2y_Px_S_S_c;
  Double I_ERI_Fxyz_Dxz_S_S_c = I_ERI_Gxy2z_Px_S_S_c+ABZ*I_ERI_Fxyz_Px_S_S_c;
  Double I_ERI_Fx2z_Dxz_S_S_c = I_ERI_Gx3z_Px_S_S_c+ABZ*I_ERI_Fx2z_Px_S_S_c;
  Double I_ERI_F3y_Dxz_S_S_c = I_ERI_G3yz_Px_S_S_c+ABZ*I_ERI_F3y_Px_S_S_c;
  Double I_ERI_F2yz_Dxz_S_S_c = I_ERI_G2y2z_Px_S_S_c+ABZ*I_ERI_F2yz_Px_S_S_c;
  Double I_ERI_Fy2z_Dxz_S_S_c = I_ERI_Gy3z_Px_S_S_c+ABZ*I_ERI_Fy2z_Px_S_S_c;
  Double I_ERI_F3z_Dxz_S_S_c = I_ERI_G4z_Px_S_S_c+ABZ*I_ERI_F3z_Px_S_S_c;
  Double I_ERI_F3x_D2y_S_S_c = I_ERI_G3xy_Py_S_S_c+ABY*I_ERI_F3x_Py_S_S_c;
  Double I_ERI_F2xy_D2y_S_S_c = I_ERI_G2x2y_Py_S_S_c+ABY*I_ERI_F2xy_Py_S_S_c;
  Double I_ERI_F2xz_D2y_S_S_c = I_ERI_G2xyz_Py_S_S_c+ABY*I_ERI_F2xz_Py_S_S_c;
  Double I_ERI_Fx2y_D2y_S_S_c = I_ERI_Gx3y_Py_S_S_c+ABY*I_ERI_Fx2y_Py_S_S_c;
  Double I_ERI_Fxyz_D2y_S_S_c = I_ERI_Gx2yz_Py_S_S_c+ABY*I_ERI_Fxyz_Py_S_S_c;
  Double I_ERI_Fx2z_D2y_S_S_c = I_ERI_Gxy2z_Py_S_S_c+ABY*I_ERI_Fx2z_Py_S_S_c;
  Double I_ERI_F3y_D2y_S_S_c = I_ERI_G4y_Py_S_S_c+ABY*I_ERI_F3y_Py_S_S_c;
  Double I_ERI_F2yz_D2y_S_S_c = I_ERI_G3yz_Py_S_S_c+ABY*I_ERI_F2yz_Py_S_S_c;
  Double I_ERI_Fy2z_D2y_S_S_c = I_ERI_G2y2z_Py_S_S_c+ABY*I_ERI_Fy2z_Py_S_S_c;
  Double I_ERI_F3z_D2y_S_S_c = I_ERI_Gy3z_Py_S_S_c+ABY*I_ERI_F3z_Py_S_S_c;
  Double I_ERI_F3x_Dyz_S_S_c = I_ERI_G3xz_Py_S_S_c+ABZ*I_ERI_F3x_Py_S_S_c;
  Double I_ERI_F2xy_Dyz_S_S_c = I_ERI_G2xyz_Py_S_S_c+ABZ*I_ERI_F2xy_Py_S_S_c;
  Double I_ERI_F2xz_Dyz_S_S_c = I_ERI_G2x2z_Py_S_S_c+ABZ*I_ERI_F2xz_Py_S_S_c;
  Double I_ERI_Fx2y_Dyz_S_S_c = I_ERI_Gx2yz_Py_S_S_c+ABZ*I_ERI_Fx2y_Py_S_S_c;
  Double I_ERI_Fxyz_Dyz_S_S_c = I_ERI_Gxy2z_Py_S_S_c+ABZ*I_ERI_Fxyz_Py_S_S_c;
  Double I_ERI_Fx2z_Dyz_S_S_c = I_ERI_Gx3z_Py_S_S_c+ABZ*I_ERI_Fx2z_Py_S_S_c;
  Double I_ERI_F3y_Dyz_S_S_c = I_ERI_G3yz_Py_S_S_c+ABZ*I_ERI_F3y_Py_S_S_c;
  Double I_ERI_F2yz_Dyz_S_S_c = I_ERI_G2y2z_Py_S_S_c+ABZ*I_ERI_F2yz_Py_S_S_c;
  Double I_ERI_Fy2z_Dyz_S_S_c = I_ERI_Gy3z_Py_S_S_c+ABZ*I_ERI_Fy2z_Py_S_S_c;
  Double I_ERI_F3z_Dyz_S_S_c = I_ERI_G4z_Py_S_S_c+ABZ*I_ERI_F3z_Py_S_S_c;
  Double I_ERI_F3x_D2z_S_S_c = I_ERI_G3xz_Pz_S_S_c+ABZ*I_ERI_F3x_Pz_S_S_c;
  Double I_ERI_F2xy_D2z_S_S_c = I_ERI_G2xyz_Pz_S_S_c+ABZ*I_ERI_F2xy_Pz_S_S_c;
  Double I_ERI_F2xz_D2z_S_S_c = I_ERI_G2x2z_Pz_S_S_c+ABZ*I_ERI_F2xz_Pz_S_S_c;
  Double I_ERI_Fx2y_D2z_S_S_c = I_ERI_Gx2yz_Pz_S_S_c+ABZ*I_ERI_Fx2y_Pz_S_S_c;
  Double I_ERI_Fxyz_D2z_S_S_c = I_ERI_Gxy2z_Pz_S_S_c+ABZ*I_ERI_Fxyz_Pz_S_S_c;
  Double I_ERI_Fx2z_D2z_S_S_c = I_ERI_Gx3z_Pz_S_S_c+ABZ*I_ERI_Fx2z_Pz_S_S_c;
  Double I_ERI_F3y_D2z_S_S_c = I_ERI_G3yz_Pz_S_S_c+ABZ*I_ERI_F3y_Pz_S_S_c;
  Double I_ERI_F2yz_D2z_S_S_c = I_ERI_G2y2z_Pz_S_S_c+ABZ*I_ERI_F2yz_Pz_S_S_c;
  Double I_ERI_Fy2z_D2z_S_S_c = I_ERI_Gy3z_Pz_S_S_c+ABZ*I_ERI_Fy2z_Pz_S_S_c;
  Double I_ERI_F3z_D2z_S_S_c = I_ERI_G4z_Pz_S_S_c+ABZ*I_ERI_F3z_Pz_S_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_H_P_S_S_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_S_S_aa
   * RHS shell quartet name: SQ_ERI_H_S_S_S_aa
   ************************************************************/
  Double I_ERI_H5x_Px_S_S_aa = I_ERI_I6x_S_S_S_aa+ABX*I_ERI_H5x_S_S_S_aa;
  Double I_ERI_H4xy_Px_S_S_aa = I_ERI_I5xy_S_S_S_aa+ABX*I_ERI_H4xy_S_S_S_aa;
  Double I_ERI_H4xz_Px_S_S_aa = I_ERI_I5xz_S_S_S_aa+ABX*I_ERI_H4xz_S_S_S_aa;
  Double I_ERI_H3x2y_Px_S_S_aa = I_ERI_I4x2y_S_S_S_aa+ABX*I_ERI_H3x2y_S_S_S_aa;
  Double I_ERI_H3xyz_Px_S_S_aa = I_ERI_I4xyz_S_S_S_aa+ABX*I_ERI_H3xyz_S_S_S_aa;
  Double I_ERI_H3x2z_Px_S_S_aa = I_ERI_I4x2z_S_S_S_aa+ABX*I_ERI_H3x2z_S_S_S_aa;
  Double I_ERI_H2x3y_Px_S_S_aa = I_ERI_I3x3y_S_S_S_aa+ABX*I_ERI_H2x3y_S_S_S_aa;
  Double I_ERI_H2x2yz_Px_S_S_aa = I_ERI_I3x2yz_S_S_S_aa+ABX*I_ERI_H2x2yz_S_S_S_aa;
  Double I_ERI_H2xy2z_Px_S_S_aa = I_ERI_I3xy2z_S_S_S_aa+ABX*I_ERI_H2xy2z_S_S_S_aa;
  Double I_ERI_H2x3z_Px_S_S_aa = I_ERI_I3x3z_S_S_S_aa+ABX*I_ERI_H2x3z_S_S_S_aa;
  Double I_ERI_Hx4y_Px_S_S_aa = I_ERI_I2x4y_S_S_S_aa+ABX*I_ERI_Hx4y_S_S_S_aa;
  Double I_ERI_Hx3yz_Px_S_S_aa = I_ERI_I2x3yz_S_S_S_aa+ABX*I_ERI_Hx3yz_S_S_S_aa;
  Double I_ERI_Hx2y2z_Px_S_S_aa = I_ERI_I2x2y2z_S_S_S_aa+ABX*I_ERI_Hx2y2z_S_S_S_aa;
  Double I_ERI_Hxy3z_Px_S_S_aa = I_ERI_I2xy3z_S_S_S_aa+ABX*I_ERI_Hxy3z_S_S_S_aa;
  Double I_ERI_Hx4z_Px_S_S_aa = I_ERI_I2x4z_S_S_S_aa+ABX*I_ERI_Hx4z_S_S_S_aa;
  Double I_ERI_H5y_Px_S_S_aa = I_ERI_Ix5y_S_S_S_aa+ABX*I_ERI_H5y_S_S_S_aa;
  Double I_ERI_H4yz_Px_S_S_aa = I_ERI_Ix4yz_S_S_S_aa+ABX*I_ERI_H4yz_S_S_S_aa;
  Double I_ERI_H3y2z_Px_S_S_aa = I_ERI_Ix3y2z_S_S_S_aa+ABX*I_ERI_H3y2z_S_S_S_aa;
  Double I_ERI_H2y3z_Px_S_S_aa = I_ERI_Ix2y3z_S_S_S_aa+ABX*I_ERI_H2y3z_S_S_S_aa;
  Double I_ERI_Hy4z_Px_S_S_aa = I_ERI_Ixy4z_S_S_S_aa+ABX*I_ERI_Hy4z_S_S_S_aa;
  Double I_ERI_H5z_Px_S_S_aa = I_ERI_Ix5z_S_S_S_aa+ABX*I_ERI_H5z_S_S_S_aa;
  Double I_ERI_H5x_Py_S_S_aa = I_ERI_I5xy_S_S_S_aa+ABY*I_ERI_H5x_S_S_S_aa;
  Double I_ERI_H4xy_Py_S_S_aa = I_ERI_I4x2y_S_S_S_aa+ABY*I_ERI_H4xy_S_S_S_aa;
  Double I_ERI_H4xz_Py_S_S_aa = I_ERI_I4xyz_S_S_S_aa+ABY*I_ERI_H4xz_S_S_S_aa;
  Double I_ERI_H3x2y_Py_S_S_aa = I_ERI_I3x3y_S_S_S_aa+ABY*I_ERI_H3x2y_S_S_S_aa;
  Double I_ERI_H3xyz_Py_S_S_aa = I_ERI_I3x2yz_S_S_S_aa+ABY*I_ERI_H3xyz_S_S_S_aa;
  Double I_ERI_H3x2z_Py_S_S_aa = I_ERI_I3xy2z_S_S_S_aa+ABY*I_ERI_H3x2z_S_S_S_aa;
  Double I_ERI_H2x3y_Py_S_S_aa = I_ERI_I2x4y_S_S_S_aa+ABY*I_ERI_H2x3y_S_S_S_aa;
  Double I_ERI_H2x2yz_Py_S_S_aa = I_ERI_I2x3yz_S_S_S_aa+ABY*I_ERI_H2x2yz_S_S_S_aa;
  Double I_ERI_H2xy2z_Py_S_S_aa = I_ERI_I2x2y2z_S_S_S_aa+ABY*I_ERI_H2xy2z_S_S_S_aa;
  Double I_ERI_H2x3z_Py_S_S_aa = I_ERI_I2xy3z_S_S_S_aa+ABY*I_ERI_H2x3z_S_S_S_aa;
  Double I_ERI_Hx4y_Py_S_S_aa = I_ERI_Ix5y_S_S_S_aa+ABY*I_ERI_Hx4y_S_S_S_aa;
  Double I_ERI_Hx3yz_Py_S_S_aa = I_ERI_Ix4yz_S_S_S_aa+ABY*I_ERI_Hx3yz_S_S_S_aa;
  Double I_ERI_Hx2y2z_Py_S_S_aa = I_ERI_Ix3y2z_S_S_S_aa+ABY*I_ERI_Hx2y2z_S_S_S_aa;
  Double I_ERI_Hxy3z_Py_S_S_aa = I_ERI_Ix2y3z_S_S_S_aa+ABY*I_ERI_Hxy3z_S_S_S_aa;
  Double I_ERI_Hx4z_Py_S_S_aa = I_ERI_Ixy4z_S_S_S_aa+ABY*I_ERI_Hx4z_S_S_S_aa;
  Double I_ERI_H5y_Py_S_S_aa = I_ERI_I6y_S_S_S_aa+ABY*I_ERI_H5y_S_S_S_aa;
  Double I_ERI_H4yz_Py_S_S_aa = I_ERI_I5yz_S_S_S_aa+ABY*I_ERI_H4yz_S_S_S_aa;
  Double I_ERI_H3y2z_Py_S_S_aa = I_ERI_I4y2z_S_S_S_aa+ABY*I_ERI_H3y2z_S_S_S_aa;
  Double I_ERI_H2y3z_Py_S_S_aa = I_ERI_I3y3z_S_S_S_aa+ABY*I_ERI_H2y3z_S_S_S_aa;
  Double I_ERI_Hy4z_Py_S_S_aa = I_ERI_I2y4z_S_S_S_aa+ABY*I_ERI_Hy4z_S_S_S_aa;
  Double I_ERI_H5z_Py_S_S_aa = I_ERI_Iy5z_S_S_S_aa+ABY*I_ERI_H5z_S_S_S_aa;
  Double I_ERI_H5x_Pz_S_S_aa = I_ERI_I5xz_S_S_S_aa+ABZ*I_ERI_H5x_S_S_S_aa;
  Double I_ERI_H4xy_Pz_S_S_aa = I_ERI_I4xyz_S_S_S_aa+ABZ*I_ERI_H4xy_S_S_S_aa;
  Double I_ERI_H4xz_Pz_S_S_aa = I_ERI_I4x2z_S_S_S_aa+ABZ*I_ERI_H4xz_S_S_S_aa;
  Double I_ERI_H3x2y_Pz_S_S_aa = I_ERI_I3x2yz_S_S_S_aa+ABZ*I_ERI_H3x2y_S_S_S_aa;
  Double I_ERI_H3xyz_Pz_S_S_aa = I_ERI_I3xy2z_S_S_S_aa+ABZ*I_ERI_H3xyz_S_S_S_aa;
  Double I_ERI_H3x2z_Pz_S_S_aa = I_ERI_I3x3z_S_S_S_aa+ABZ*I_ERI_H3x2z_S_S_S_aa;
  Double I_ERI_H2x3y_Pz_S_S_aa = I_ERI_I2x3yz_S_S_S_aa+ABZ*I_ERI_H2x3y_S_S_S_aa;
  Double I_ERI_H2x2yz_Pz_S_S_aa = I_ERI_I2x2y2z_S_S_S_aa+ABZ*I_ERI_H2x2yz_S_S_S_aa;
  Double I_ERI_H2xy2z_Pz_S_S_aa = I_ERI_I2xy3z_S_S_S_aa+ABZ*I_ERI_H2xy2z_S_S_S_aa;
  Double I_ERI_H2x3z_Pz_S_S_aa = I_ERI_I2x4z_S_S_S_aa+ABZ*I_ERI_H2x3z_S_S_S_aa;
  Double I_ERI_Hx4y_Pz_S_S_aa = I_ERI_Ix4yz_S_S_S_aa+ABZ*I_ERI_Hx4y_S_S_S_aa;
  Double I_ERI_Hx3yz_Pz_S_S_aa = I_ERI_Ix3y2z_S_S_S_aa+ABZ*I_ERI_Hx3yz_S_S_S_aa;
  Double I_ERI_Hx2y2z_Pz_S_S_aa = I_ERI_Ix2y3z_S_S_S_aa+ABZ*I_ERI_Hx2y2z_S_S_S_aa;
  Double I_ERI_Hxy3z_Pz_S_S_aa = I_ERI_Ixy4z_S_S_S_aa+ABZ*I_ERI_Hxy3z_S_S_S_aa;
  Double I_ERI_Hx4z_Pz_S_S_aa = I_ERI_Ix5z_S_S_S_aa+ABZ*I_ERI_Hx4z_S_S_S_aa;
  Double I_ERI_H5y_Pz_S_S_aa = I_ERI_I5yz_S_S_S_aa+ABZ*I_ERI_H5y_S_S_S_aa;
  Double I_ERI_H4yz_Pz_S_S_aa = I_ERI_I4y2z_S_S_S_aa+ABZ*I_ERI_H4yz_S_S_S_aa;
  Double I_ERI_H3y2z_Pz_S_S_aa = I_ERI_I3y3z_S_S_S_aa+ABZ*I_ERI_H3y2z_S_S_S_aa;
  Double I_ERI_H2y3z_Pz_S_S_aa = I_ERI_I2y4z_S_S_S_aa+ABZ*I_ERI_H2y3z_S_S_S_aa;
  Double I_ERI_Hy4z_Pz_S_S_aa = I_ERI_Iy5z_S_S_S_aa+ABZ*I_ERI_Hy4z_S_S_S_aa;
  Double I_ERI_H5z_Pz_S_S_aa = I_ERI_I6z_S_S_S_aa+ABZ*I_ERI_H5z_S_S_S_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_I_P_S_S_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 8 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_K_S_S_S_aa
   * RHS shell quartet name: SQ_ERI_I_S_S_S_aa
   ************************************************************/
  Double I_ERI_I6x_Px_S_S_aa = I_ERI_K7x_S_S_S_aa+ABX*I_ERI_I6x_S_S_S_aa;
  Double I_ERI_I5xy_Px_S_S_aa = I_ERI_K6xy_S_S_S_aa+ABX*I_ERI_I5xy_S_S_S_aa;
  Double I_ERI_I5xz_Px_S_S_aa = I_ERI_K6xz_S_S_S_aa+ABX*I_ERI_I5xz_S_S_S_aa;
  Double I_ERI_I4x2y_Px_S_S_aa = I_ERI_K5x2y_S_S_S_aa+ABX*I_ERI_I4x2y_S_S_S_aa;
  Double I_ERI_I4xyz_Px_S_S_aa = I_ERI_K5xyz_S_S_S_aa+ABX*I_ERI_I4xyz_S_S_S_aa;
  Double I_ERI_I4x2z_Px_S_S_aa = I_ERI_K5x2z_S_S_S_aa+ABX*I_ERI_I4x2z_S_S_S_aa;
  Double I_ERI_I3x3y_Px_S_S_aa = I_ERI_K4x3y_S_S_S_aa+ABX*I_ERI_I3x3y_S_S_S_aa;
  Double I_ERI_I3x2yz_Px_S_S_aa = I_ERI_K4x2yz_S_S_S_aa+ABX*I_ERI_I3x2yz_S_S_S_aa;
  Double I_ERI_I3xy2z_Px_S_S_aa = I_ERI_K4xy2z_S_S_S_aa+ABX*I_ERI_I3xy2z_S_S_S_aa;
  Double I_ERI_I3x3z_Px_S_S_aa = I_ERI_K4x3z_S_S_S_aa+ABX*I_ERI_I3x3z_S_S_S_aa;
  Double I_ERI_I2x4y_Px_S_S_aa = I_ERI_K3x4y_S_S_S_aa+ABX*I_ERI_I2x4y_S_S_S_aa;
  Double I_ERI_I2x3yz_Px_S_S_aa = I_ERI_K3x3yz_S_S_S_aa+ABX*I_ERI_I2x3yz_S_S_S_aa;
  Double I_ERI_I2x2y2z_Px_S_S_aa = I_ERI_K3x2y2z_S_S_S_aa+ABX*I_ERI_I2x2y2z_S_S_S_aa;
  Double I_ERI_I2xy3z_Px_S_S_aa = I_ERI_K3xy3z_S_S_S_aa+ABX*I_ERI_I2xy3z_S_S_S_aa;
  Double I_ERI_I2x4z_Px_S_S_aa = I_ERI_K3x4z_S_S_S_aa+ABX*I_ERI_I2x4z_S_S_S_aa;
  Double I_ERI_Ix5y_Px_S_S_aa = I_ERI_K2x5y_S_S_S_aa+ABX*I_ERI_Ix5y_S_S_S_aa;
  Double I_ERI_Ix4yz_Px_S_S_aa = I_ERI_K2x4yz_S_S_S_aa+ABX*I_ERI_Ix4yz_S_S_S_aa;
  Double I_ERI_Ix3y2z_Px_S_S_aa = I_ERI_K2x3y2z_S_S_S_aa+ABX*I_ERI_Ix3y2z_S_S_S_aa;
  Double I_ERI_Ix2y3z_Px_S_S_aa = I_ERI_K2x2y3z_S_S_S_aa+ABX*I_ERI_Ix2y3z_S_S_S_aa;
  Double I_ERI_Ixy4z_Px_S_S_aa = I_ERI_K2xy4z_S_S_S_aa+ABX*I_ERI_Ixy4z_S_S_S_aa;
  Double I_ERI_Ix5z_Px_S_S_aa = I_ERI_K2x5z_S_S_S_aa+ABX*I_ERI_Ix5z_S_S_S_aa;
  Double I_ERI_I6y_Px_S_S_aa = I_ERI_Kx6y_S_S_S_aa+ABX*I_ERI_I6y_S_S_S_aa;
  Double I_ERI_I5yz_Px_S_S_aa = I_ERI_Kx5yz_S_S_S_aa+ABX*I_ERI_I5yz_S_S_S_aa;
  Double I_ERI_I4y2z_Px_S_S_aa = I_ERI_Kx4y2z_S_S_S_aa+ABX*I_ERI_I4y2z_S_S_S_aa;
  Double I_ERI_I3y3z_Px_S_S_aa = I_ERI_Kx3y3z_S_S_S_aa+ABX*I_ERI_I3y3z_S_S_S_aa;
  Double I_ERI_I2y4z_Px_S_S_aa = I_ERI_Kx2y4z_S_S_S_aa+ABX*I_ERI_I2y4z_S_S_S_aa;
  Double I_ERI_Iy5z_Px_S_S_aa = I_ERI_Kxy5z_S_S_S_aa+ABX*I_ERI_Iy5z_S_S_S_aa;
  Double I_ERI_I6z_Px_S_S_aa = I_ERI_Kx6z_S_S_S_aa+ABX*I_ERI_I6z_S_S_S_aa;
  Double I_ERI_I5xy_Py_S_S_aa = I_ERI_K5x2y_S_S_S_aa+ABY*I_ERI_I5xy_S_S_S_aa;
  Double I_ERI_I5xz_Py_S_S_aa = I_ERI_K5xyz_S_S_S_aa+ABY*I_ERI_I5xz_S_S_S_aa;
  Double I_ERI_I4x2y_Py_S_S_aa = I_ERI_K4x3y_S_S_S_aa+ABY*I_ERI_I4x2y_S_S_S_aa;
  Double I_ERI_I4xyz_Py_S_S_aa = I_ERI_K4x2yz_S_S_S_aa+ABY*I_ERI_I4xyz_S_S_S_aa;
  Double I_ERI_I4x2z_Py_S_S_aa = I_ERI_K4xy2z_S_S_S_aa+ABY*I_ERI_I4x2z_S_S_S_aa;
  Double I_ERI_I3x3y_Py_S_S_aa = I_ERI_K3x4y_S_S_S_aa+ABY*I_ERI_I3x3y_S_S_S_aa;
  Double I_ERI_I3x2yz_Py_S_S_aa = I_ERI_K3x3yz_S_S_S_aa+ABY*I_ERI_I3x2yz_S_S_S_aa;
  Double I_ERI_I3xy2z_Py_S_S_aa = I_ERI_K3x2y2z_S_S_S_aa+ABY*I_ERI_I3xy2z_S_S_S_aa;
  Double I_ERI_I3x3z_Py_S_S_aa = I_ERI_K3xy3z_S_S_S_aa+ABY*I_ERI_I3x3z_S_S_S_aa;
  Double I_ERI_I2x4y_Py_S_S_aa = I_ERI_K2x5y_S_S_S_aa+ABY*I_ERI_I2x4y_S_S_S_aa;
  Double I_ERI_I2x3yz_Py_S_S_aa = I_ERI_K2x4yz_S_S_S_aa+ABY*I_ERI_I2x3yz_S_S_S_aa;
  Double I_ERI_I2x2y2z_Py_S_S_aa = I_ERI_K2x3y2z_S_S_S_aa+ABY*I_ERI_I2x2y2z_S_S_S_aa;
  Double I_ERI_I2xy3z_Py_S_S_aa = I_ERI_K2x2y3z_S_S_S_aa+ABY*I_ERI_I2xy3z_S_S_S_aa;
  Double I_ERI_I2x4z_Py_S_S_aa = I_ERI_K2xy4z_S_S_S_aa+ABY*I_ERI_I2x4z_S_S_S_aa;
  Double I_ERI_Ix5y_Py_S_S_aa = I_ERI_Kx6y_S_S_S_aa+ABY*I_ERI_Ix5y_S_S_S_aa;
  Double I_ERI_Ix4yz_Py_S_S_aa = I_ERI_Kx5yz_S_S_S_aa+ABY*I_ERI_Ix4yz_S_S_S_aa;
  Double I_ERI_Ix3y2z_Py_S_S_aa = I_ERI_Kx4y2z_S_S_S_aa+ABY*I_ERI_Ix3y2z_S_S_S_aa;
  Double I_ERI_Ix2y3z_Py_S_S_aa = I_ERI_Kx3y3z_S_S_S_aa+ABY*I_ERI_Ix2y3z_S_S_S_aa;
  Double I_ERI_Ixy4z_Py_S_S_aa = I_ERI_Kx2y4z_S_S_S_aa+ABY*I_ERI_Ixy4z_S_S_S_aa;
  Double I_ERI_Ix5z_Py_S_S_aa = I_ERI_Kxy5z_S_S_S_aa+ABY*I_ERI_Ix5z_S_S_S_aa;
  Double I_ERI_I6y_Py_S_S_aa = I_ERI_K7y_S_S_S_aa+ABY*I_ERI_I6y_S_S_S_aa;
  Double I_ERI_I5yz_Py_S_S_aa = I_ERI_K6yz_S_S_S_aa+ABY*I_ERI_I5yz_S_S_S_aa;
  Double I_ERI_I4y2z_Py_S_S_aa = I_ERI_K5y2z_S_S_S_aa+ABY*I_ERI_I4y2z_S_S_S_aa;
  Double I_ERI_I3y3z_Py_S_S_aa = I_ERI_K4y3z_S_S_S_aa+ABY*I_ERI_I3y3z_S_S_S_aa;
  Double I_ERI_I2y4z_Py_S_S_aa = I_ERI_K3y4z_S_S_S_aa+ABY*I_ERI_I2y4z_S_S_S_aa;
  Double I_ERI_Iy5z_Py_S_S_aa = I_ERI_K2y5z_S_S_S_aa+ABY*I_ERI_Iy5z_S_S_S_aa;
  Double I_ERI_I6z_Py_S_S_aa = I_ERI_Ky6z_S_S_S_aa+ABY*I_ERI_I6z_S_S_S_aa;
  Double I_ERI_I5xz_Pz_S_S_aa = I_ERI_K5x2z_S_S_S_aa+ABZ*I_ERI_I5xz_S_S_S_aa;
  Double I_ERI_I4xyz_Pz_S_S_aa = I_ERI_K4xy2z_S_S_S_aa+ABZ*I_ERI_I4xyz_S_S_S_aa;
  Double I_ERI_I4x2z_Pz_S_S_aa = I_ERI_K4x3z_S_S_S_aa+ABZ*I_ERI_I4x2z_S_S_S_aa;
  Double I_ERI_I3x2yz_Pz_S_S_aa = I_ERI_K3x2y2z_S_S_S_aa+ABZ*I_ERI_I3x2yz_S_S_S_aa;
  Double I_ERI_I3xy2z_Pz_S_S_aa = I_ERI_K3xy3z_S_S_S_aa+ABZ*I_ERI_I3xy2z_S_S_S_aa;
  Double I_ERI_I3x3z_Pz_S_S_aa = I_ERI_K3x4z_S_S_S_aa+ABZ*I_ERI_I3x3z_S_S_S_aa;
  Double I_ERI_I2x3yz_Pz_S_S_aa = I_ERI_K2x3y2z_S_S_S_aa+ABZ*I_ERI_I2x3yz_S_S_S_aa;
  Double I_ERI_I2x2y2z_Pz_S_S_aa = I_ERI_K2x2y3z_S_S_S_aa+ABZ*I_ERI_I2x2y2z_S_S_S_aa;
  Double I_ERI_I2xy3z_Pz_S_S_aa = I_ERI_K2xy4z_S_S_S_aa+ABZ*I_ERI_I2xy3z_S_S_S_aa;
  Double I_ERI_I2x4z_Pz_S_S_aa = I_ERI_K2x5z_S_S_S_aa+ABZ*I_ERI_I2x4z_S_S_S_aa;
  Double I_ERI_Ix4yz_Pz_S_S_aa = I_ERI_Kx4y2z_S_S_S_aa+ABZ*I_ERI_Ix4yz_S_S_S_aa;
  Double I_ERI_Ix3y2z_Pz_S_S_aa = I_ERI_Kx3y3z_S_S_S_aa+ABZ*I_ERI_Ix3y2z_S_S_S_aa;
  Double I_ERI_Ix2y3z_Pz_S_S_aa = I_ERI_Kx2y4z_S_S_S_aa+ABZ*I_ERI_Ix2y3z_S_S_S_aa;
  Double I_ERI_Ixy4z_Pz_S_S_aa = I_ERI_Kxy5z_S_S_S_aa+ABZ*I_ERI_Ixy4z_S_S_S_aa;
  Double I_ERI_Ix5z_Pz_S_S_aa = I_ERI_Kx6z_S_S_S_aa+ABZ*I_ERI_Ix5z_S_S_S_aa;
  Double I_ERI_I5yz_Pz_S_S_aa = I_ERI_K5y2z_S_S_S_aa+ABZ*I_ERI_I5yz_S_S_S_aa;
  Double I_ERI_I4y2z_Pz_S_S_aa = I_ERI_K4y3z_S_S_S_aa+ABZ*I_ERI_I4y2z_S_S_S_aa;
  Double I_ERI_I3y3z_Pz_S_S_aa = I_ERI_K3y4z_S_S_S_aa+ABZ*I_ERI_I3y3z_S_S_S_aa;
  Double I_ERI_I2y4z_Pz_S_S_aa = I_ERI_K2y5z_S_S_S_aa+ABZ*I_ERI_I2y4z_S_S_S_aa;
  Double I_ERI_Iy5z_Pz_S_S_aa = I_ERI_Ky6z_S_S_S_aa+ABZ*I_ERI_Iy5z_S_S_S_aa;
  Double I_ERI_I6z_Pz_S_S_aa = I_ERI_K7z_S_S_S_aa+ABZ*I_ERI_I6z_S_S_S_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_H_D_S_S_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_P_S_S_aa
   * RHS shell quartet name: SQ_ERI_H_P_S_S_aa
   ************************************************************/
  Double I_ERI_H5x_D2x_S_S_aa = I_ERI_I6x_Px_S_S_aa+ABX*I_ERI_H5x_Px_S_S_aa;
  Double I_ERI_H4xy_D2x_S_S_aa = I_ERI_I5xy_Px_S_S_aa+ABX*I_ERI_H4xy_Px_S_S_aa;
  Double I_ERI_H4xz_D2x_S_S_aa = I_ERI_I5xz_Px_S_S_aa+ABX*I_ERI_H4xz_Px_S_S_aa;
  Double I_ERI_H3x2y_D2x_S_S_aa = I_ERI_I4x2y_Px_S_S_aa+ABX*I_ERI_H3x2y_Px_S_S_aa;
  Double I_ERI_H3xyz_D2x_S_S_aa = I_ERI_I4xyz_Px_S_S_aa+ABX*I_ERI_H3xyz_Px_S_S_aa;
  Double I_ERI_H3x2z_D2x_S_S_aa = I_ERI_I4x2z_Px_S_S_aa+ABX*I_ERI_H3x2z_Px_S_S_aa;
  Double I_ERI_H2x3y_D2x_S_S_aa = I_ERI_I3x3y_Px_S_S_aa+ABX*I_ERI_H2x3y_Px_S_S_aa;
  Double I_ERI_H2x2yz_D2x_S_S_aa = I_ERI_I3x2yz_Px_S_S_aa+ABX*I_ERI_H2x2yz_Px_S_S_aa;
  Double I_ERI_H2xy2z_D2x_S_S_aa = I_ERI_I3xy2z_Px_S_S_aa+ABX*I_ERI_H2xy2z_Px_S_S_aa;
  Double I_ERI_H2x3z_D2x_S_S_aa = I_ERI_I3x3z_Px_S_S_aa+ABX*I_ERI_H2x3z_Px_S_S_aa;
  Double I_ERI_Hx4y_D2x_S_S_aa = I_ERI_I2x4y_Px_S_S_aa+ABX*I_ERI_Hx4y_Px_S_S_aa;
  Double I_ERI_Hx3yz_D2x_S_S_aa = I_ERI_I2x3yz_Px_S_S_aa+ABX*I_ERI_Hx3yz_Px_S_S_aa;
  Double I_ERI_Hx2y2z_D2x_S_S_aa = I_ERI_I2x2y2z_Px_S_S_aa+ABX*I_ERI_Hx2y2z_Px_S_S_aa;
  Double I_ERI_Hxy3z_D2x_S_S_aa = I_ERI_I2xy3z_Px_S_S_aa+ABX*I_ERI_Hxy3z_Px_S_S_aa;
  Double I_ERI_Hx4z_D2x_S_S_aa = I_ERI_I2x4z_Px_S_S_aa+ABX*I_ERI_Hx4z_Px_S_S_aa;
  Double I_ERI_H5y_D2x_S_S_aa = I_ERI_Ix5y_Px_S_S_aa+ABX*I_ERI_H5y_Px_S_S_aa;
  Double I_ERI_H4yz_D2x_S_S_aa = I_ERI_Ix4yz_Px_S_S_aa+ABX*I_ERI_H4yz_Px_S_S_aa;
  Double I_ERI_H3y2z_D2x_S_S_aa = I_ERI_Ix3y2z_Px_S_S_aa+ABX*I_ERI_H3y2z_Px_S_S_aa;
  Double I_ERI_H2y3z_D2x_S_S_aa = I_ERI_Ix2y3z_Px_S_S_aa+ABX*I_ERI_H2y3z_Px_S_S_aa;
  Double I_ERI_Hy4z_D2x_S_S_aa = I_ERI_Ixy4z_Px_S_S_aa+ABX*I_ERI_Hy4z_Px_S_S_aa;
  Double I_ERI_H5z_D2x_S_S_aa = I_ERI_Ix5z_Px_S_S_aa+ABX*I_ERI_H5z_Px_S_S_aa;
  Double I_ERI_H5x_Dxy_S_S_aa = I_ERI_I5xy_Px_S_S_aa+ABY*I_ERI_H5x_Px_S_S_aa;
  Double I_ERI_H4xy_Dxy_S_S_aa = I_ERI_I4x2y_Px_S_S_aa+ABY*I_ERI_H4xy_Px_S_S_aa;
  Double I_ERI_H4xz_Dxy_S_S_aa = I_ERI_I4xyz_Px_S_S_aa+ABY*I_ERI_H4xz_Px_S_S_aa;
  Double I_ERI_H3x2y_Dxy_S_S_aa = I_ERI_I3x3y_Px_S_S_aa+ABY*I_ERI_H3x2y_Px_S_S_aa;
  Double I_ERI_H3xyz_Dxy_S_S_aa = I_ERI_I3x2yz_Px_S_S_aa+ABY*I_ERI_H3xyz_Px_S_S_aa;
  Double I_ERI_H3x2z_Dxy_S_S_aa = I_ERI_I3xy2z_Px_S_S_aa+ABY*I_ERI_H3x2z_Px_S_S_aa;
  Double I_ERI_H2x3y_Dxy_S_S_aa = I_ERI_I2x4y_Px_S_S_aa+ABY*I_ERI_H2x3y_Px_S_S_aa;
  Double I_ERI_H2x2yz_Dxy_S_S_aa = I_ERI_I2x3yz_Px_S_S_aa+ABY*I_ERI_H2x2yz_Px_S_S_aa;
  Double I_ERI_H2xy2z_Dxy_S_S_aa = I_ERI_I2x2y2z_Px_S_S_aa+ABY*I_ERI_H2xy2z_Px_S_S_aa;
  Double I_ERI_H2x3z_Dxy_S_S_aa = I_ERI_I2xy3z_Px_S_S_aa+ABY*I_ERI_H2x3z_Px_S_S_aa;
  Double I_ERI_Hx4y_Dxy_S_S_aa = I_ERI_Ix5y_Px_S_S_aa+ABY*I_ERI_Hx4y_Px_S_S_aa;
  Double I_ERI_Hx3yz_Dxy_S_S_aa = I_ERI_Ix4yz_Px_S_S_aa+ABY*I_ERI_Hx3yz_Px_S_S_aa;
  Double I_ERI_Hx2y2z_Dxy_S_S_aa = I_ERI_Ix3y2z_Px_S_S_aa+ABY*I_ERI_Hx2y2z_Px_S_S_aa;
  Double I_ERI_Hxy3z_Dxy_S_S_aa = I_ERI_Ix2y3z_Px_S_S_aa+ABY*I_ERI_Hxy3z_Px_S_S_aa;
  Double I_ERI_Hx4z_Dxy_S_S_aa = I_ERI_Ixy4z_Px_S_S_aa+ABY*I_ERI_Hx4z_Px_S_S_aa;
  Double I_ERI_H5y_Dxy_S_S_aa = I_ERI_I6y_Px_S_S_aa+ABY*I_ERI_H5y_Px_S_S_aa;
  Double I_ERI_H4yz_Dxy_S_S_aa = I_ERI_I5yz_Px_S_S_aa+ABY*I_ERI_H4yz_Px_S_S_aa;
  Double I_ERI_H3y2z_Dxy_S_S_aa = I_ERI_I4y2z_Px_S_S_aa+ABY*I_ERI_H3y2z_Px_S_S_aa;
  Double I_ERI_H2y3z_Dxy_S_S_aa = I_ERI_I3y3z_Px_S_S_aa+ABY*I_ERI_H2y3z_Px_S_S_aa;
  Double I_ERI_Hy4z_Dxy_S_S_aa = I_ERI_I2y4z_Px_S_S_aa+ABY*I_ERI_Hy4z_Px_S_S_aa;
  Double I_ERI_H5z_Dxy_S_S_aa = I_ERI_Iy5z_Px_S_S_aa+ABY*I_ERI_H5z_Px_S_S_aa;
  Double I_ERI_H5x_Dxz_S_S_aa = I_ERI_I5xz_Px_S_S_aa+ABZ*I_ERI_H5x_Px_S_S_aa;
  Double I_ERI_H4xy_Dxz_S_S_aa = I_ERI_I4xyz_Px_S_S_aa+ABZ*I_ERI_H4xy_Px_S_S_aa;
  Double I_ERI_H4xz_Dxz_S_S_aa = I_ERI_I4x2z_Px_S_S_aa+ABZ*I_ERI_H4xz_Px_S_S_aa;
  Double I_ERI_H3x2y_Dxz_S_S_aa = I_ERI_I3x2yz_Px_S_S_aa+ABZ*I_ERI_H3x2y_Px_S_S_aa;
  Double I_ERI_H3xyz_Dxz_S_S_aa = I_ERI_I3xy2z_Px_S_S_aa+ABZ*I_ERI_H3xyz_Px_S_S_aa;
  Double I_ERI_H3x2z_Dxz_S_S_aa = I_ERI_I3x3z_Px_S_S_aa+ABZ*I_ERI_H3x2z_Px_S_S_aa;
  Double I_ERI_H2x3y_Dxz_S_S_aa = I_ERI_I2x3yz_Px_S_S_aa+ABZ*I_ERI_H2x3y_Px_S_S_aa;
  Double I_ERI_H2x2yz_Dxz_S_S_aa = I_ERI_I2x2y2z_Px_S_S_aa+ABZ*I_ERI_H2x2yz_Px_S_S_aa;
  Double I_ERI_H2xy2z_Dxz_S_S_aa = I_ERI_I2xy3z_Px_S_S_aa+ABZ*I_ERI_H2xy2z_Px_S_S_aa;
  Double I_ERI_H2x3z_Dxz_S_S_aa = I_ERI_I2x4z_Px_S_S_aa+ABZ*I_ERI_H2x3z_Px_S_S_aa;
  Double I_ERI_Hx4y_Dxz_S_S_aa = I_ERI_Ix4yz_Px_S_S_aa+ABZ*I_ERI_Hx4y_Px_S_S_aa;
  Double I_ERI_Hx3yz_Dxz_S_S_aa = I_ERI_Ix3y2z_Px_S_S_aa+ABZ*I_ERI_Hx3yz_Px_S_S_aa;
  Double I_ERI_Hx2y2z_Dxz_S_S_aa = I_ERI_Ix2y3z_Px_S_S_aa+ABZ*I_ERI_Hx2y2z_Px_S_S_aa;
  Double I_ERI_Hxy3z_Dxz_S_S_aa = I_ERI_Ixy4z_Px_S_S_aa+ABZ*I_ERI_Hxy3z_Px_S_S_aa;
  Double I_ERI_Hx4z_Dxz_S_S_aa = I_ERI_Ix5z_Px_S_S_aa+ABZ*I_ERI_Hx4z_Px_S_S_aa;
  Double I_ERI_H5y_Dxz_S_S_aa = I_ERI_I5yz_Px_S_S_aa+ABZ*I_ERI_H5y_Px_S_S_aa;
  Double I_ERI_H4yz_Dxz_S_S_aa = I_ERI_I4y2z_Px_S_S_aa+ABZ*I_ERI_H4yz_Px_S_S_aa;
  Double I_ERI_H3y2z_Dxz_S_S_aa = I_ERI_I3y3z_Px_S_S_aa+ABZ*I_ERI_H3y2z_Px_S_S_aa;
  Double I_ERI_H2y3z_Dxz_S_S_aa = I_ERI_I2y4z_Px_S_S_aa+ABZ*I_ERI_H2y3z_Px_S_S_aa;
  Double I_ERI_Hy4z_Dxz_S_S_aa = I_ERI_Iy5z_Px_S_S_aa+ABZ*I_ERI_Hy4z_Px_S_S_aa;
  Double I_ERI_H5z_Dxz_S_S_aa = I_ERI_I6z_Px_S_S_aa+ABZ*I_ERI_H5z_Px_S_S_aa;
  Double I_ERI_H5x_D2y_S_S_aa = I_ERI_I5xy_Py_S_S_aa+ABY*I_ERI_H5x_Py_S_S_aa;
  Double I_ERI_H4xy_D2y_S_S_aa = I_ERI_I4x2y_Py_S_S_aa+ABY*I_ERI_H4xy_Py_S_S_aa;
  Double I_ERI_H4xz_D2y_S_S_aa = I_ERI_I4xyz_Py_S_S_aa+ABY*I_ERI_H4xz_Py_S_S_aa;
  Double I_ERI_H3x2y_D2y_S_S_aa = I_ERI_I3x3y_Py_S_S_aa+ABY*I_ERI_H3x2y_Py_S_S_aa;
  Double I_ERI_H3xyz_D2y_S_S_aa = I_ERI_I3x2yz_Py_S_S_aa+ABY*I_ERI_H3xyz_Py_S_S_aa;
  Double I_ERI_H3x2z_D2y_S_S_aa = I_ERI_I3xy2z_Py_S_S_aa+ABY*I_ERI_H3x2z_Py_S_S_aa;
  Double I_ERI_H2x3y_D2y_S_S_aa = I_ERI_I2x4y_Py_S_S_aa+ABY*I_ERI_H2x3y_Py_S_S_aa;
  Double I_ERI_H2x2yz_D2y_S_S_aa = I_ERI_I2x3yz_Py_S_S_aa+ABY*I_ERI_H2x2yz_Py_S_S_aa;
  Double I_ERI_H2xy2z_D2y_S_S_aa = I_ERI_I2x2y2z_Py_S_S_aa+ABY*I_ERI_H2xy2z_Py_S_S_aa;
  Double I_ERI_H2x3z_D2y_S_S_aa = I_ERI_I2xy3z_Py_S_S_aa+ABY*I_ERI_H2x3z_Py_S_S_aa;
  Double I_ERI_Hx4y_D2y_S_S_aa = I_ERI_Ix5y_Py_S_S_aa+ABY*I_ERI_Hx4y_Py_S_S_aa;
  Double I_ERI_Hx3yz_D2y_S_S_aa = I_ERI_Ix4yz_Py_S_S_aa+ABY*I_ERI_Hx3yz_Py_S_S_aa;
  Double I_ERI_Hx2y2z_D2y_S_S_aa = I_ERI_Ix3y2z_Py_S_S_aa+ABY*I_ERI_Hx2y2z_Py_S_S_aa;
  Double I_ERI_Hxy3z_D2y_S_S_aa = I_ERI_Ix2y3z_Py_S_S_aa+ABY*I_ERI_Hxy3z_Py_S_S_aa;
  Double I_ERI_Hx4z_D2y_S_S_aa = I_ERI_Ixy4z_Py_S_S_aa+ABY*I_ERI_Hx4z_Py_S_S_aa;
  Double I_ERI_H5y_D2y_S_S_aa = I_ERI_I6y_Py_S_S_aa+ABY*I_ERI_H5y_Py_S_S_aa;
  Double I_ERI_H4yz_D2y_S_S_aa = I_ERI_I5yz_Py_S_S_aa+ABY*I_ERI_H4yz_Py_S_S_aa;
  Double I_ERI_H3y2z_D2y_S_S_aa = I_ERI_I4y2z_Py_S_S_aa+ABY*I_ERI_H3y2z_Py_S_S_aa;
  Double I_ERI_H2y3z_D2y_S_S_aa = I_ERI_I3y3z_Py_S_S_aa+ABY*I_ERI_H2y3z_Py_S_S_aa;
  Double I_ERI_Hy4z_D2y_S_S_aa = I_ERI_I2y4z_Py_S_S_aa+ABY*I_ERI_Hy4z_Py_S_S_aa;
  Double I_ERI_H5z_D2y_S_S_aa = I_ERI_Iy5z_Py_S_S_aa+ABY*I_ERI_H5z_Py_S_S_aa;
  Double I_ERI_H5x_Dyz_S_S_aa = I_ERI_I5xz_Py_S_S_aa+ABZ*I_ERI_H5x_Py_S_S_aa;
  Double I_ERI_H4xy_Dyz_S_S_aa = I_ERI_I4xyz_Py_S_S_aa+ABZ*I_ERI_H4xy_Py_S_S_aa;
  Double I_ERI_H4xz_Dyz_S_S_aa = I_ERI_I4x2z_Py_S_S_aa+ABZ*I_ERI_H4xz_Py_S_S_aa;
  Double I_ERI_H3x2y_Dyz_S_S_aa = I_ERI_I3x2yz_Py_S_S_aa+ABZ*I_ERI_H3x2y_Py_S_S_aa;
  Double I_ERI_H3xyz_Dyz_S_S_aa = I_ERI_I3xy2z_Py_S_S_aa+ABZ*I_ERI_H3xyz_Py_S_S_aa;
  Double I_ERI_H3x2z_Dyz_S_S_aa = I_ERI_I3x3z_Py_S_S_aa+ABZ*I_ERI_H3x2z_Py_S_S_aa;
  Double I_ERI_H2x3y_Dyz_S_S_aa = I_ERI_I2x3yz_Py_S_S_aa+ABZ*I_ERI_H2x3y_Py_S_S_aa;
  Double I_ERI_H2x2yz_Dyz_S_S_aa = I_ERI_I2x2y2z_Py_S_S_aa+ABZ*I_ERI_H2x2yz_Py_S_S_aa;
  Double I_ERI_H2xy2z_Dyz_S_S_aa = I_ERI_I2xy3z_Py_S_S_aa+ABZ*I_ERI_H2xy2z_Py_S_S_aa;
  Double I_ERI_H2x3z_Dyz_S_S_aa = I_ERI_I2x4z_Py_S_S_aa+ABZ*I_ERI_H2x3z_Py_S_S_aa;
  Double I_ERI_Hx4y_Dyz_S_S_aa = I_ERI_Ix4yz_Py_S_S_aa+ABZ*I_ERI_Hx4y_Py_S_S_aa;
  Double I_ERI_Hx3yz_Dyz_S_S_aa = I_ERI_Ix3y2z_Py_S_S_aa+ABZ*I_ERI_Hx3yz_Py_S_S_aa;
  Double I_ERI_Hx2y2z_Dyz_S_S_aa = I_ERI_Ix2y3z_Py_S_S_aa+ABZ*I_ERI_Hx2y2z_Py_S_S_aa;
  Double I_ERI_Hxy3z_Dyz_S_S_aa = I_ERI_Ixy4z_Py_S_S_aa+ABZ*I_ERI_Hxy3z_Py_S_S_aa;
  Double I_ERI_Hx4z_Dyz_S_S_aa = I_ERI_Ix5z_Py_S_S_aa+ABZ*I_ERI_Hx4z_Py_S_S_aa;
  Double I_ERI_H5y_Dyz_S_S_aa = I_ERI_I5yz_Py_S_S_aa+ABZ*I_ERI_H5y_Py_S_S_aa;
  Double I_ERI_H4yz_Dyz_S_S_aa = I_ERI_I4y2z_Py_S_S_aa+ABZ*I_ERI_H4yz_Py_S_S_aa;
  Double I_ERI_H3y2z_Dyz_S_S_aa = I_ERI_I3y3z_Py_S_S_aa+ABZ*I_ERI_H3y2z_Py_S_S_aa;
  Double I_ERI_H2y3z_Dyz_S_S_aa = I_ERI_I2y4z_Py_S_S_aa+ABZ*I_ERI_H2y3z_Py_S_S_aa;
  Double I_ERI_Hy4z_Dyz_S_S_aa = I_ERI_Iy5z_Py_S_S_aa+ABZ*I_ERI_Hy4z_Py_S_S_aa;
  Double I_ERI_H5z_Dyz_S_S_aa = I_ERI_I6z_Py_S_S_aa+ABZ*I_ERI_H5z_Py_S_S_aa;
  Double I_ERI_H5x_D2z_S_S_aa = I_ERI_I5xz_Pz_S_S_aa+ABZ*I_ERI_H5x_Pz_S_S_aa;
  Double I_ERI_H4xy_D2z_S_S_aa = I_ERI_I4xyz_Pz_S_S_aa+ABZ*I_ERI_H4xy_Pz_S_S_aa;
  Double I_ERI_H4xz_D2z_S_S_aa = I_ERI_I4x2z_Pz_S_S_aa+ABZ*I_ERI_H4xz_Pz_S_S_aa;
  Double I_ERI_H3x2y_D2z_S_S_aa = I_ERI_I3x2yz_Pz_S_S_aa+ABZ*I_ERI_H3x2y_Pz_S_S_aa;
  Double I_ERI_H3xyz_D2z_S_S_aa = I_ERI_I3xy2z_Pz_S_S_aa+ABZ*I_ERI_H3xyz_Pz_S_S_aa;
  Double I_ERI_H3x2z_D2z_S_S_aa = I_ERI_I3x3z_Pz_S_S_aa+ABZ*I_ERI_H3x2z_Pz_S_S_aa;
  Double I_ERI_H2x3y_D2z_S_S_aa = I_ERI_I2x3yz_Pz_S_S_aa+ABZ*I_ERI_H2x3y_Pz_S_S_aa;
  Double I_ERI_H2x2yz_D2z_S_S_aa = I_ERI_I2x2y2z_Pz_S_S_aa+ABZ*I_ERI_H2x2yz_Pz_S_S_aa;
  Double I_ERI_H2xy2z_D2z_S_S_aa = I_ERI_I2xy3z_Pz_S_S_aa+ABZ*I_ERI_H2xy2z_Pz_S_S_aa;
  Double I_ERI_H2x3z_D2z_S_S_aa = I_ERI_I2x4z_Pz_S_S_aa+ABZ*I_ERI_H2x3z_Pz_S_S_aa;
  Double I_ERI_Hx4y_D2z_S_S_aa = I_ERI_Ix4yz_Pz_S_S_aa+ABZ*I_ERI_Hx4y_Pz_S_S_aa;
  Double I_ERI_Hx3yz_D2z_S_S_aa = I_ERI_Ix3y2z_Pz_S_S_aa+ABZ*I_ERI_Hx3yz_Pz_S_S_aa;
  Double I_ERI_Hx2y2z_D2z_S_S_aa = I_ERI_Ix2y3z_Pz_S_S_aa+ABZ*I_ERI_Hx2y2z_Pz_S_S_aa;
  Double I_ERI_Hxy3z_D2z_S_S_aa = I_ERI_Ixy4z_Pz_S_S_aa+ABZ*I_ERI_Hxy3z_Pz_S_S_aa;
  Double I_ERI_Hx4z_D2z_S_S_aa = I_ERI_Ix5z_Pz_S_S_aa+ABZ*I_ERI_Hx4z_Pz_S_S_aa;
  Double I_ERI_H5y_D2z_S_S_aa = I_ERI_I5yz_Pz_S_S_aa+ABZ*I_ERI_H5y_Pz_S_S_aa;
  Double I_ERI_H4yz_D2z_S_S_aa = I_ERI_I4y2z_Pz_S_S_aa+ABZ*I_ERI_H4yz_Pz_S_S_aa;
  Double I_ERI_H3y2z_D2z_S_S_aa = I_ERI_I3y3z_Pz_S_S_aa+ABZ*I_ERI_H3y2z_Pz_S_S_aa;
  Double I_ERI_H2y3z_D2z_S_S_aa = I_ERI_I2y4z_Pz_S_S_aa+ABZ*I_ERI_H2y3z_Pz_S_S_aa;
  Double I_ERI_Hy4z_D2z_S_S_aa = I_ERI_Iy5z_Pz_S_S_aa+ABZ*I_ERI_Hy4z_Pz_S_S_aa;
  Double I_ERI_H5z_D2z_S_S_aa = I_ERI_I6z_Pz_S_S_aa+ABZ*I_ERI_H5z_Pz_S_S_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_S_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_S_S_S_ab
   ************************************************************/
  Double I_ERI_G4x_Px_S_S_ab = I_ERI_H5x_S_S_S_ab+ABX*I_ERI_G4x_S_S_S_ab;
  Double I_ERI_G3xy_Px_S_S_ab = I_ERI_H4xy_S_S_S_ab+ABX*I_ERI_G3xy_S_S_S_ab;
  Double I_ERI_G3xz_Px_S_S_ab = I_ERI_H4xz_S_S_S_ab+ABX*I_ERI_G3xz_S_S_S_ab;
  Double I_ERI_G2x2y_Px_S_S_ab = I_ERI_H3x2y_S_S_S_ab+ABX*I_ERI_G2x2y_S_S_S_ab;
  Double I_ERI_G2xyz_Px_S_S_ab = I_ERI_H3xyz_S_S_S_ab+ABX*I_ERI_G2xyz_S_S_S_ab;
  Double I_ERI_G2x2z_Px_S_S_ab = I_ERI_H3x2z_S_S_S_ab+ABX*I_ERI_G2x2z_S_S_S_ab;
  Double I_ERI_Gx3y_Px_S_S_ab = I_ERI_H2x3y_S_S_S_ab+ABX*I_ERI_Gx3y_S_S_S_ab;
  Double I_ERI_Gx2yz_Px_S_S_ab = I_ERI_H2x2yz_S_S_S_ab+ABX*I_ERI_Gx2yz_S_S_S_ab;
  Double I_ERI_Gxy2z_Px_S_S_ab = I_ERI_H2xy2z_S_S_S_ab+ABX*I_ERI_Gxy2z_S_S_S_ab;
  Double I_ERI_Gx3z_Px_S_S_ab = I_ERI_H2x3z_S_S_S_ab+ABX*I_ERI_Gx3z_S_S_S_ab;
  Double I_ERI_G4y_Px_S_S_ab = I_ERI_Hx4y_S_S_S_ab+ABX*I_ERI_G4y_S_S_S_ab;
  Double I_ERI_G3yz_Px_S_S_ab = I_ERI_Hx3yz_S_S_S_ab+ABX*I_ERI_G3yz_S_S_S_ab;
  Double I_ERI_G2y2z_Px_S_S_ab = I_ERI_Hx2y2z_S_S_S_ab+ABX*I_ERI_G2y2z_S_S_S_ab;
  Double I_ERI_Gy3z_Px_S_S_ab = I_ERI_Hxy3z_S_S_S_ab+ABX*I_ERI_Gy3z_S_S_S_ab;
  Double I_ERI_G4z_Px_S_S_ab = I_ERI_Hx4z_S_S_S_ab+ABX*I_ERI_G4z_S_S_S_ab;
  Double I_ERI_G4x_Py_S_S_ab = I_ERI_H4xy_S_S_S_ab+ABY*I_ERI_G4x_S_S_S_ab;
  Double I_ERI_G3xy_Py_S_S_ab = I_ERI_H3x2y_S_S_S_ab+ABY*I_ERI_G3xy_S_S_S_ab;
  Double I_ERI_G3xz_Py_S_S_ab = I_ERI_H3xyz_S_S_S_ab+ABY*I_ERI_G3xz_S_S_S_ab;
  Double I_ERI_G2x2y_Py_S_S_ab = I_ERI_H2x3y_S_S_S_ab+ABY*I_ERI_G2x2y_S_S_S_ab;
  Double I_ERI_G2xyz_Py_S_S_ab = I_ERI_H2x2yz_S_S_S_ab+ABY*I_ERI_G2xyz_S_S_S_ab;
  Double I_ERI_G2x2z_Py_S_S_ab = I_ERI_H2xy2z_S_S_S_ab+ABY*I_ERI_G2x2z_S_S_S_ab;
  Double I_ERI_Gx3y_Py_S_S_ab = I_ERI_Hx4y_S_S_S_ab+ABY*I_ERI_Gx3y_S_S_S_ab;
  Double I_ERI_Gx2yz_Py_S_S_ab = I_ERI_Hx3yz_S_S_S_ab+ABY*I_ERI_Gx2yz_S_S_S_ab;
  Double I_ERI_Gxy2z_Py_S_S_ab = I_ERI_Hx2y2z_S_S_S_ab+ABY*I_ERI_Gxy2z_S_S_S_ab;
  Double I_ERI_Gx3z_Py_S_S_ab = I_ERI_Hxy3z_S_S_S_ab+ABY*I_ERI_Gx3z_S_S_S_ab;
  Double I_ERI_G4y_Py_S_S_ab = I_ERI_H5y_S_S_S_ab+ABY*I_ERI_G4y_S_S_S_ab;
  Double I_ERI_G3yz_Py_S_S_ab = I_ERI_H4yz_S_S_S_ab+ABY*I_ERI_G3yz_S_S_S_ab;
  Double I_ERI_G2y2z_Py_S_S_ab = I_ERI_H3y2z_S_S_S_ab+ABY*I_ERI_G2y2z_S_S_S_ab;
  Double I_ERI_Gy3z_Py_S_S_ab = I_ERI_H2y3z_S_S_S_ab+ABY*I_ERI_Gy3z_S_S_S_ab;
  Double I_ERI_G4z_Py_S_S_ab = I_ERI_Hy4z_S_S_S_ab+ABY*I_ERI_G4z_S_S_S_ab;
  Double I_ERI_G4x_Pz_S_S_ab = I_ERI_H4xz_S_S_S_ab+ABZ*I_ERI_G4x_S_S_S_ab;
  Double I_ERI_G3xy_Pz_S_S_ab = I_ERI_H3xyz_S_S_S_ab+ABZ*I_ERI_G3xy_S_S_S_ab;
  Double I_ERI_G3xz_Pz_S_S_ab = I_ERI_H3x2z_S_S_S_ab+ABZ*I_ERI_G3xz_S_S_S_ab;
  Double I_ERI_G2x2y_Pz_S_S_ab = I_ERI_H2x2yz_S_S_S_ab+ABZ*I_ERI_G2x2y_S_S_S_ab;
  Double I_ERI_G2xyz_Pz_S_S_ab = I_ERI_H2xy2z_S_S_S_ab+ABZ*I_ERI_G2xyz_S_S_S_ab;
  Double I_ERI_G2x2z_Pz_S_S_ab = I_ERI_H2x3z_S_S_S_ab+ABZ*I_ERI_G2x2z_S_S_S_ab;
  Double I_ERI_Gx3y_Pz_S_S_ab = I_ERI_Hx3yz_S_S_S_ab+ABZ*I_ERI_Gx3y_S_S_S_ab;
  Double I_ERI_Gx2yz_Pz_S_S_ab = I_ERI_Hx2y2z_S_S_S_ab+ABZ*I_ERI_Gx2yz_S_S_S_ab;
  Double I_ERI_Gxy2z_Pz_S_S_ab = I_ERI_Hxy3z_S_S_S_ab+ABZ*I_ERI_Gxy2z_S_S_S_ab;
  Double I_ERI_Gx3z_Pz_S_S_ab = I_ERI_Hx4z_S_S_S_ab+ABZ*I_ERI_Gx3z_S_S_S_ab;
  Double I_ERI_G4y_Pz_S_S_ab = I_ERI_H4yz_S_S_S_ab+ABZ*I_ERI_G4y_S_S_S_ab;
  Double I_ERI_G3yz_Pz_S_S_ab = I_ERI_H3y2z_S_S_S_ab+ABZ*I_ERI_G3yz_S_S_S_ab;
  Double I_ERI_G2y2z_Pz_S_S_ab = I_ERI_H2y3z_S_S_S_ab+ABZ*I_ERI_G2y2z_S_S_S_ab;
  Double I_ERI_Gy3z_Pz_S_S_ab = I_ERI_Hy4z_S_S_S_ab+ABZ*I_ERI_Gy3z_S_S_S_ab;
  Double I_ERI_G4z_Pz_S_S_ab = I_ERI_H5z_S_S_S_ab+ABZ*I_ERI_G4z_S_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_H_P_S_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_S_S_ab
   * RHS shell quartet name: SQ_ERI_H_S_S_S_ab
   ************************************************************/
  Double I_ERI_H5x_Px_S_S_ab = I_ERI_I6x_S_S_S_ab+ABX*I_ERI_H5x_S_S_S_ab;
  Double I_ERI_H4xy_Px_S_S_ab = I_ERI_I5xy_S_S_S_ab+ABX*I_ERI_H4xy_S_S_S_ab;
  Double I_ERI_H4xz_Px_S_S_ab = I_ERI_I5xz_S_S_S_ab+ABX*I_ERI_H4xz_S_S_S_ab;
  Double I_ERI_H3x2y_Px_S_S_ab = I_ERI_I4x2y_S_S_S_ab+ABX*I_ERI_H3x2y_S_S_S_ab;
  Double I_ERI_H3xyz_Px_S_S_ab = I_ERI_I4xyz_S_S_S_ab+ABX*I_ERI_H3xyz_S_S_S_ab;
  Double I_ERI_H3x2z_Px_S_S_ab = I_ERI_I4x2z_S_S_S_ab+ABX*I_ERI_H3x2z_S_S_S_ab;
  Double I_ERI_H2x3y_Px_S_S_ab = I_ERI_I3x3y_S_S_S_ab+ABX*I_ERI_H2x3y_S_S_S_ab;
  Double I_ERI_H2x2yz_Px_S_S_ab = I_ERI_I3x2yz_S_S_S_ab+ABX*I_ERI_H2x2yz_S_S_S_ab;
  Double I_ERI_H2xy2z_Px_S_S_ab = I_ERI_I3xy2z_S_S_S_ab+ABX*I_ERI_H2xy2z_S_S_S_ab;
  Double I_ERI_H2x3z_Px_S_S_ab = I_ERI_I3x3z_S_S_S_ab+ABX*I_ERI_H2x3z_S_S_S_ab;
  Double I_ERI_Hx4y_Px_S_S_ab = I_ERI_I2x4y_S_S_S_ab+ABX*I_ERI_Hx4y_S_S_S_ab;
  Double I_ERI_Hx3yz_Px_S_S_ab = I_ERI_I2x3yz_S_S_S_ab+ABX*I_ERI_Hx3yz_S_S_S_ab;
  Double I_ERI_Hx2y2z_Px_S_S_ab = I_ERI_I2x2y2z_S_S_S_ab+ABX*I_ERI_Hx2y2z_S_S_S_ab;
  Double I_ERI_Hxy3z_Px_S_S_ab = I_ERI_I2xy3z_S_S_S_ab+ABX*I_ERI_Hxy3z_S_S_S_ab;
  Double I_ERI_Hx4z_Px_S_S_ab = I_ERI_I2x4z_S_S_S_ab+ABX*I_ERI_Hx4z_S_S_S_ab;
  Double I_ERI_H5y_Px_S_S_ab = I_ERI_Ix5y_S_S_S_ab+ABX*I_ERI_H5y_S_S_S_ab;
  Double I_ERI_H4yz_Px_S_S_ab = I_ERI_Ix4yz_S_S_S_ab+ABX*I_ERI_H4yz_S_S_S_ab;
  Double I_ERI_H3y2z_Px_S_S_ab = I_ERI_Ix3y2z_S_S_S_ab+ABX*I_ERI_H3y2z_S_S_S_ab;
  Double I_ERI_H2y3z_Px_S_S_ab = I_ERI_Ix2y3z_S_S_S_ab+ABX*I_ERI_H2y3z_S_S_S_ab;
  Double I_ERI_Hy4z_Px_S_S_ab = I_ERI_Ixy4z_S_S_S_ab+ABX*I_ERI_Hy4z_S_S_S_ab;
  Double I_ERI_H5z_Px_S_S_ab = I_ERI_Ix5z_S_S_S_ab+ABX*I_ERI_H5z_S_S_S_ab;
  Double I_ERI_H5x_Py_S_S_ab = I_ERI_I5xy_S_S_S_ab+ABY*I_ERI_H5x_S_S_S_ab;
  Double I_ERI_H4xy_Py_S_S_ab = I_ERI_I4x2y_S_S_S_ab+ABY*I_ERI_H4xy_S_S_S_ab;
  Double I_ERI_H4xz_Py_S_S_ab = I_ERI_I4xyz_S_S_S_ab+ABY*I_ERI_H4xz_S_S_S_ab;
  Double I_ERI_H3x2y_Py_S_S_ab = I_ERI_I3x3y_S_S_S_ab+ABY*I_ERI_H3x2y_S_S_S_ab;
  Double I_ERI_H3xyz_Py_S_S_ab = I_ERI_I3x2yz_S_S_S_ab+ABY*I_ERI_H3xyz_S_S_S_ab;
  Double I_ERI_H3x2z_Py_S_S_ab = I_ERI_I3xy2z_S_S_S_ab+ABY*I_ERI_H3x2z_S_S_S_ab;
  Double I_ERI_H2x3y_Py_S_S_ab = I_ERI_I2x4y_S_S_S_ab+ABY*I_ERI_H2x3y_S_S_S_ab;
  Double I_ERI_H2x2yz_Py_S_S_ab = I_ERI_I2x3yz_S_S_S_ab+ABY*I_ERI_H2x2yz_S_S_S_ab;
  Double I_ERI_H2xy2z_Py_S_S_ab = I_ERI_I2x2y2z_S_S_S_ab+ABY*I_ERI_H2xy2z_S_S_S_ab;
  Double I_ERI_H2x3z_Py_S_S_ab = I_ERI_I2xy3z_S_S_S_ab+ABY*I_ERI_H2x3z_S_S_S_ab;
  Double I_ERI_Hx4y_Py_S_S_ab = I_ERI_Ix5y_S_S_S_ab+ABY*I_ERI_Hx4y_S_S_S_ab;
  Double I_ERI_Hx3yz_Py_S_S_ab = I_ERI_Ix4yz_S_S_S_ab+ABY*I_ERI_Hx3yz_S_S_S_ab;
  Double I_ERI_Hx2y2z_Py_S_S_ab = I_ERI_Ix3y2z_S_S_S_ab+ABY*I_ERI_Hx2y2z_S_S_S_ab;
  Double I_ERI_Hxy3z_Py_S_S_ab = I_ERI_Ix2y3z_S_S_S_ab+ABY*I_ERI_Hxy3z_S_S_S_ab;
  Double I_ERI_Hx4z_Py_S_S_ab = I_ERI_Ixy4z_S_S_S_ab+ABY*I_ERI_Hx4z_S_S_S_ab;
  Double I_ERI_H5y_Py_S_S_ab = I_ERI_I6y_S_S_S_ab+ABY*I_ERI_H5y_S_S_S_ab;
  Double I_ERI_H4yz_Py_S_S_ab = I_ERI_I5yz_S_S_S_ab+ABY*I_ERI_H4yz_S_S_S_ab;
  Double I_ERI_H3y2z_Py_S_S_ab = I_ERI_I4y2z_S_S_S_ab+ABY*I_ERI_H3y2z_S_S_S_ab;
  Double I_ERI_H2y3z_Py_S_S_ab = I_ERI_I3y3z_S_S_S_ab+ABY*I_ERI_H2y3z_S_S_S_ab;
  Double I_ERI_Hy4z_Py_S_S_ab = I_ERI_I2y4z_S_S_S_ab+ABY*I_ERI_Hy4z_S_S_S_ab;
  Double I_ERI_H5z_Py_S_S_ab = I_ERI_Iy5z_S_S_S_ab+ABY*I_ERI_H5z_S_S_S_ab;
  Double I_ERI_H5x_Pz_S_S_ab = I_ERI_I5xz_S_S_S_ab+ABZ*I_ERI_H5x_S_S_S_ab;
  Double I_ERI_H4xy_Pz_S_S_ab = I_ERI_I4xyz_S_S_S_ab+ABZ*I_ERI_H4xy_S_S_S_ab;
  Double I_ERI_H4xz_Pz_S_S_ab = I_ERI_I4x2z_S_S_S_ab+ABZ*I_ERI_H4xz_S_S_S_ab;
  Double I_ERI_H3x2y_Pz_S_S_ab = I_ERI_I3x2yz_S_S_S_ab+ABZ*I_ERI_H3x2y_S_S_S_ab;
  Double I_ERI_H3xyz_Pz_S_S_ab = I_ERI_I3xy2z_S_S_S_ab+ABZ*I_ERI_H3xyz_S_S_S_ab;
  Double I_ERI_H3x2z_Pz_S_S_ab = I_ERI_I3x3z_S_S_S_ab+ABZ*I_ERI_H3x2z_S_S_S_ab;
  Double I_ERI_H2x3y_Pz_S_S_ab = I_ERI_I2x3yz_S_S_S_ab+ABZ*I_ERI_H2x3y_S_S_S_ab;
  Double I_ERI_H2x2yz_Pz_S_S_ab = I_ERI_I2x2y2z_S_S_S_ab+ABZ*I_ERI_H2x2yz_S_S_S_ab;
  Double I_ERI_H2xy2z_Pz_S_S_ab = I_ERI_I2xy3z_S_S_S_ab+ABZ*I_ERI_H2xy2z_S_S_S_ab;
  Double I_ERI_H2x3z_Pz_S_S_ab = I_ERI_I2x4z_S_S_S_ab+ABZ*I_ERI_H2x3z_S_S_S_ab;
  Double I_ERI_Hx4y_Pz_S_S_ab = I_ERI_Ix4yz_S_S_S_ab+ABZ*I_ERI_Hx4y_S_S_S_ab;
  Double I_ERI_Hx3yz_Pz_S_S_ab = I_ERI_Ix3y2z_S_S_S_ab+ABZ*I_ERI_Hx3yz_S_S_S_ab;
  Double I_ERI_Hx2y2z_Pz_S_S_ab = I_ERI_Ix2y3z_S_S_S_ab+ABZ*I_ERI_Hx2y2z_S_S_S_ab;
  Double I_ERI_Hxy3z_Pz_S_S_ab = I_ERI_Ixy4z_S_S_S_ab+ABZ*I_ERI_Hxy3z_S_S_S_ab;
  Double I_ERI_Hx4z_Pz_S_S_ab = I_ERI_Ix5z_S_S_S_ab+ABZ*I_ERI_Hx4z_S_S_S_ab;
  Double I_ERI_H5y_Pz_S_S_ab = I_ERI_I5yz_S_S_S_ab+ABZ*I_ERI_H5y_S_S_S_ab;
  Double I_ERI_H4yz_Pz_S_S_ab = I_ERI_I4y2z_S_S_S_ab+ABZ*I_ERI_H4yz_S_S_S_ab;
  Double I_ERI_H3y2z_Pz_S_S_ab = I_ERI_I3y3z_S_S_S_ab+ABZ*I_ERI_H3y2z_S_S_S_ab;
  Double I_ERI_H2y3z_Pz_S_S_ab = I_ERI_I2y4z_S_S_S_ab+ABZ*I_ERI_H2y3z_S_S_S_ab;
  Double I_ERI_Hy4z_Pz_S_S_ab = I_ERI_Iy5z_S_S_S_ab+ABZ*I_ERI_Hy4z_S_S_S_ab;
  Double I_ERI_H5z_Pz_S_S_ab = I_ERI_I6z_S_S_S_ab+ABZ*I_ERI_H5z_S_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_G_D_S_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_ab
   ************************************************************/
  Double I_ERI_G4x_D2x_S_S_ab = I_ERI_H5x_Px_S_S_ab+ABX*I_ERI_G4x_Px_S_S_ab;
  Double I_ERI_G3xy_D2x_S_S_ab = I_ERI_H4xy_Px_S_S_ab+ABX*I_ERI_G3xy_Px_S_S_ab;
  Double I_ERI_G3xz_D2x_S_S_ab = I_ERI_H4xz_Px_S_S_ab+ABX*I_ERI_G3xz_Px_S_S_ab;
  Double I_ERI_G2x2y_D2x_S_S_ab = I_ERI_H3x2y_Px_S_S_ab+ABX*I_ERI_G2x2y_Px_S_S_ab;
  Double I_ERI_G2xyz_D2x_S_S_ab = I_ERI_H3xyz_Px_S_S_ab+ABX*I_ERI_G2xyz_Px_S_S_ab;
  Double I_ERI_G2x2z_D2x_S_S_ab = I_ERI_H3x2z_Px_S_S_ab+ABX*I_ERI_G2x2z_Px_S_S_ab;
  Double I_ERI_Gx3y_D2x_S_S_ab = I_ERI_H2x3y_Px_S_S_ab+ABX*I_ERI_Gx3y_Px_S_S_ab;
  Double I_ERI_Gx2yz_D2x_S_S_ab = I_ERI_H2x2yz_Px_S_S_ab+ABX*I_ERI_Gx2yz_Px_S_S_ab;
  Double I_ERI_Gxy2z_D2x_S_S_ab = I_ERI_H2xy2z_Px_S_S_ab+ABX*I_ERI_Gxy2z_Px_S_S_ab;
  Double I_ERI_Gx3z_D2x_S_S_ab = I_ERI_H2x3z_Px_S_S_ab+ABX*I_ERI_Gx3z_Px_S_S_ab;
  Double I_ERI_G4y_D2x_S_S_ab = I_ERI_Hx4y_Px_S_S_ab+ABX*I_ERI_G4y_Px_S_S_ab;
  Double I_ERI_G3yz_D2x_S_S_ab = I_ERI_Hx3yz_Px_S_S_ab+ABX*I_ERI_G3yz_Px_S_S_ab;
  Double I_ERI_G2y2z_D2x_S_S_ab = I_ERI_Hx2y2z_Px_S_S_ab+ABX*I_ERI_G2y2z_Px_S_S_ab;
  Double I_ERI_Gy3z_D2x_S_S_ab = I_ERI_Hxy3z_Px_S_S_ab+ABX*I_ERI_Gy3z_Px_S_S_ab;
  Double I_ERI_G4z_D2x_S_S_ab = I_ERI_Hx4z_Px_S_S_ab+ABX*I_ERI_G4z_Px_S_S_ab;
  Double I_ERI_G4x_Dxy_S_S_ab = I_ERI_H4xy_Px_S_S_ab+ABY*I_ERI_G4x_Px_S_S_ab;
  Double I_ERI_G3xy_Dxy_S_S_ab = I_ERI_H3x2y_Px_S_S_ab+ABY*I_ERI_G3xy_Px_S_S_ab;
  Double I_ERI_G3xz_Dxy_S_S_ab = I_ERI_H3xyz_Px_S_S_ab+ABY*I_ERI_G3xz_Px_S_S_ab;
  Double I_ERI_G2x2y_Dxy_S_S_ab = I_ERI_H2x3y_Px_S_S_ab+ABY*I_ERI_G2x2y_Px_S_S_ab;
  Double I_ERI_G2xyz_Dxy_S_S_ab = I_ERI_H2x2yz_Px_S_S_ab+ABY*I_ERI_G2xyz_Px_S_S_ab;
  Double I_ERI_G2x2z_Dxy_S_S_ab = I_ERI_H2xy2z_Px_S_S_ab+ABY*I_ERI_G2x2z_Px_S_S_ab;
  Double I_ERI_Gx3y_Dxy_S_S_ab = I_ERI_Hx4y_Px_S_S_ab+ABY*I_ERI_Gx3y_Px_S_S_ab;
  Double I_ERI_Gx2yz_Dxy_S_S_ab = I_ERI_Hx3yz_Px_S_S_ab+ABY*I_ERI_Gx2yz_Px_S_S_ab;
  Double I_ERI_Gxy2z_Dxy_S_S_ab = I_ERI_Hx2y2z_Px_S_S_ab+ABY*I_ERI_Gxy2z_Px_S_S_ab;
  Double I_ERI_Gx3z_Dxy_S_S_ab = I_ERI_Hxy3z_Px_S_S_ab+ABY*I_ERI_Gx3z_Px_S_S_ab;
  Double I_ERI_G4y_Dxy_S_S_ab = I_ERI_H5y_Px_S_S_ab+ABY*I_ERI_G4y_Px_S_S_ab;
  Double I_ERI_G3yz_Dxy_S_S_ab = I_ERI_H4yz_Px_S_S_ab+ABY*I_ERI_G3yz_Px_S_S_ab;
  Double I_ERI_G2y2z_Dxy_S_S_ab = I_ERI_H3y2z_Px_S_S_ab+ABY*I_ERI_G2y2z_Px_S_S_ab;
  Double I_ERI_Gy3z_Dxy_S_S_ab = I_ERI_H2y3z_Px_S_S_ab+ABY*I_ERI_Gy3z_Px_S_S_ab;
  Double I_ERI_G4z_Dxy_S_S_ab = I_ERI_Hy4z_Px_S_S_ab+ABY*I_ERI_G4z_Px_S_S_ab;
  Double I_ERI_G4x_D2y_S_S_ab = I_ERI_H4xy_Py_S_S_ab+ABY*I_ERI_G4x_Py_S_S_ab;
  Double I_ERI_G3xy_D2y_S_S_ab = I_ERI_H3x2y_Py_S_S_ab+ABY*I_ERI_G3xy_Py_S_S_ab;
  Double I_ERI_G3xz_D2y_S_S_ab = I_ERI_H3xyz_Py_S_S_ab+ABY*I_ERI_G3xz_Py_S_S_ab;
  Double I_ERI_G2x2y_D2y_S_S_ab = I_ERI_H2x3y_Py_S_S_ab+ABY*I_ERI_G2x2y_Py_S_S_ab;
  Double I_ERI_G2xyz_D2y_S_S_ab = I_ERI_H2x2yz_Py_S_S_ab+ABY*I_ERI_G2xyz_Py_S_S_ab;
  Double I_ERI_G2x2z_D2y_S_S_ab = I_ERI_H2xy2z_Py_S_S_ab+ABY*I_ERI_G2x2z_Py_S_S_ab;
  Double I_ERI_Gx3y_D2y_S_S_ab = I_ERI_Hx4y_Py_S_S_ab+ABY*I_ERI_Gx3y_Py_S_S_ab;
  Double I_ERI_Gx2yz_D2y_S_S_ab = I_ERI_Hx3yz_Py_S_S_ab+ABY*I_ERI_Gx2yz_Py_S_S_ab;
  Double I_ERI_Gxy2z_D2y_S_S_ab = I_ERI_Hx2y2z_Py_S_S_ab+ABY*I_ERI_Gxy2z_Py_S_S_ab;
  Double I_ERI_Gx3z_D2y_S_S_ab = I_ERI_Hxy3z_Py_S_S_ab+ABY*I_ERI_Gx3z_Py_S_S_ab;
  Double I_ERI_G4y_D2y_S_S_ab = I_ERI_H5y_Py_S_S_ab+ABY*I_ERI_G4y_Py_S_S_ab;
  Double I_ERI_G3yz_D2y_S_S_ab = I_ERI_H4yz_Py_S_S_ab+ABY*I_ERI_G3yz_Py_S_S_ab;
  Double I_ERI_G2y2z_D2y_S_S_ab = I_ERI_H3y2z_Py_S_S_ab+ABY*I_ERI_G2y2z_Py_S_S_ab;
  Double I_ERI_Gy3z_D2y_S_S_ab = I_ERI_H2y3z_Py_S_S_ab+ABY*I_ERI_Gy3z_Py_S_S_ab;
  Double I_ERI_G4z_D2y_S_S_ab = I_ERI_Hy4z_Py_S_S_ab+ABY*I_ERI_G4z_Py_S_S_ab;
  Double I_ERI_G4x_D2z_S_S_ab = I_ERI_H4xz_Pz_S_S_ab+ABZ*I_ERI_G4x_Pz_S_S_ab;
  Double I_ERI_G3xy_D2z_S_S_ab = I_ERI_H3xyz_Pz_S_S_ab+ABZ*I_ERI_G3xy_Pz_S_S_ab;
  Double I_ERI_G3xz_D2z_S_S_ab = I_ERI_H3x2z_Pz_S_S_ab+ABZ*I_ERI_G3xz_Pz_S_S_ab;
  Double I_ERI_G2x2y_D2z_S_S_ab = I_ERI_H2x2yz_Pz_S_S_ab+ABZ*I_ERI_G2x2y_Pz_S_S_ab;
  Double I_ERI_G2xyz_D2z_S_S_ab = I_ERI_H2xy2z_Pz_S_S_ab+ABZ*I_ERI_G2xyz_Pz_S_S_ab;
  Double I_ERI_G2x2z_D2z_S_S_ab = I_ERI_H2x3z_Pz_S_S_ab+ABZ*I_ERI_G2x2z_Pz_S_S_ab;
  Double I_ERI_Gx3y_D2z_S_S_ab = I_ERI_Hx3yz_Pz_S_S_ab+ABZ*I_ERI_Gx3y_Pz_S_S_ab;
  Double I_ERI_Gx2yz_D2z_S_S_ab = I_ERI_Hx2y2z_Pz_S_S_ab+ABZ*I_ERI_Gx2yz_Pz_S_S_ab;
  Double I_ERI_Gxy2z_D2z_S_S_ab = I_ERI_Hxy3z_Pz_S_S_ab+ABZ*I_ERI_Gxy2z_Pz_S_S_ab;
  Double I_ERI_Gx3z_D2z_S_S_ab = I_ERI_Hx4z_Pz_S_S_ab+ABZ*I_ERI_Gx3z_Pz_S_S_ab;
  Double I_ERI_G4y_D2z_S_S_ab = I_ERI_H4yz_Pz_S_S_ab+ABZ*I_ERI_G4y_Pz_S_S_ab;
  Double I_ERI_G3yz_D2z_S_S_ab = I_ERI_H3y2z_Pz_S_S_ab+ABZ*I_ERI_G3yz_Pz_S_S_ab;
  Double I_ERI_G2y2z_D2z_S_S_ab = I_ERI_H2y3z_Pz_S_S_ab+ABZ*I_ERI_G2y2z_Pz_S_S_ab;
  Double I_ERI_Gy3z_D2z_S_S_ab = I_ERI_Hy4z_Pz_S_S_ab+ABZ*I_ERI_Gy3z_Pz_S_S_ab;
  Double I_ERI_G4z_D2z_S_S_ab = I_ERI_H5z_Pz_S_S_ab+ABZ*I_ERI_G4z_Pz_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_I_P_S_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 16 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_K_S_S_S_ab
   * RHS shell quartet name: SQ_ERI_I_S_S_S_ab
   ************************************************************/
  Double I_ERI_I6x_Px_S_S_ab = I_ERI_K7x_S_S_S_ab+ABX*I_ERI_I6x_S_S_S_ab;
  Double I_ERI_I5xy_Px_S_S_ab = I_ERI_K6xy_S_S_S_ab+ABX*I_ERI_I5xy_S_S_S_ab;
  Double I_ERI_I5xz_Px_S_S_ab = I_ERI_K6xz_S_S_S_ab+ABX*I_ERI_I5xz_S_S_S_ab;
  Double I_ERI_I4x2y_Px_S_S_ab = I_ERI_K5x2y_S_S_S_ab+ABX*I_ERI_I4x2y_S_S_S_ab;
  Double I_ERI_I4xyz_Px_S_S_ab = I_ERI_K5xyz_S_S_S_ab+ABX*I_ERI_I4xyz_S_S_S_ab;
  Double I_ERI_I4x2z_Px_S_S_ab = I_ERI_K5x2z_S_S_S_ab+ABX*I_ERI_I4x2z_S_S_S_ab;
  Double I_ERI_I3x3y_Px_S_S_ab = I_ERI_K4x3y_S_S_S_ab+ABX*I_ERI_I3x3y_S_S_S_ab;
  Double I_ERI_I3x2yz_Px_S_S_ab = I_ERI_K4x2yz_S_S_S_ab+ABX*I_ERI_I3x2yz_S_S_S_ab;
  Double I_ERI_I3xy2z_Px_S_S_ab = I_ERI_K4xy2z_S_S_S_ab+ABX*I_ERI_I3xy2z_S_S_S_ab;
  Double I_ERI_I3x3z_Px_S_S_ab = I_ERI_K4x3z_S_S_S_ab+ABX*I_ERI_I3x3z_S_S_S_ab;
  Double I_ERI_I2x4y_Px_S_S_ab = I_ERI_K3x4y_S_S_S_ab+ABX*I_ERI_I2x4y_S_S_S_ab;
  Double I_ERI_I2x3yz_Px_S_S_ab = I_ERI_K3x3yz_S_S_S_ab+ABX*I_ERI_I2x3yz_S_S_S_ab;
  Double I_ERI_I2x2y2z_Px_S_S_ab = I_ERI_K3x2y2z_S_S_S_ab+ABX*I_ERI_I2x2y2z_S_S_S_ab;
  Double I_ERI_I2xy3z_Px_S_S_ab = I_ERI_K3xy3z_S_S_S_ab+ABX*I_ERI_I2xy3z_S_S_S_ab;
  Double I_ERI_I2x4z_Px_S_S_ab = I_ERI_K3x4z_S_S_S_ab+ABX*I_ERI_I2x4z_S_S_S_ab;
  Double I_ERI_Ix5y_Px_S_S_ab = I_ERI_K2x5y_S_S_S_ab+ABX*I_ERI_Ix5y_S_S_S_ab;
  Double I_ERI_Ix4yz_Px_S_S_ab = I_ERI_K2x4yz_S_S_S_ab+ABX*I_ERI_Ix4yz_S_S_S_ab;
  Double I_ERI_Ix3y2z_Px_S_S_ab = I_ERI_K2x3y2z_S_S_S_ab+ABX*I_ERI_Ix3y2z_S_S_S_ab;
  Double I_ERI_Ix2y3z_Px_S_S_ab = I_ERI_K2x2y3z_S_S_S_ab+ABX*I_ERI_Ix2y3z_S_S_S_ab;
  Double I_ERI_Ixy4z_Px_S_S_ab = I_ERI_K2xy4z_S_S_S_ab+ABX*I_ERI_Ixy4z_S_S_S_ab;
  Double I_ERI_Ix5z_Px_S_S_ab = I_ERI_K2x5z_S_S_S_ab+ABX*I_ERI_Ix5z_S_S_S_ab;
  Double I_ERI_I5yz_Px_S_S_ab = I_ERI_Kx5yz_S_S_S_ab+ABX*I_ERI_I5yz_S_S_S_ab;
  Double I_ERI_I4y2z_Px_S_S_ab = I_ERI_Kx4y2z_S_S_S_ab+ABX*I_ERI_I4y2z_S_S_S_ab;
  Double I_ERI_I3y3z_Px_S_S_ab = I_ERI_Kx3y3z_S_S_S_ab+ABX*I_ERI_I3y3z_S_S_S_ab;
  Double I_ERI_I2y4z_Px_S_S_ab = I_ERI_Kx2y4z_S_S_S_ab+ABX*I_ERI_I2y4z_S_S_S_ab;
  Double I_ERI_Iy5z_Px_S_S_ab = I_ERI_Kxy5z_S_S_S_ab+ABX*I_ERI_Iy5z_S_S_S_ab;
  Double I_ERI_I5xy_Py_S_S_ab = I_ERI_K5x2y_S_S_S_ab+ABY*I_ERI_I5xy_S_S_S_ab;
  Double I_ERI_I4x2y_Py_S_S_ab = I_ERI_K4x3y_S_S_S_ab+ABY*I_ERI_I4x2y_S_S_S_ab;
  Double I_ERI_I4xyz_Py_S_S_ab = I_ERI_K4x2yz_S_S_S_ab+ABY*I_ERI_I4xyz_S_S_S_ab;
  Double I_ERI_I3x3y_Py_S_S_ab = I_ERI_K3x4y_S_S_S_ab+ABY*I_ERI_I3x3y_S_S_S_ab;
  Double I_ERI_I3x2yz_Py_S_S_ab = I_ERI_K3x3yz_S_S_S_ab+ABY*I_ERI_I3x2yz_S_S_S_ab;
  Double I_ERI_I3xy2z_Py_S_S_ab = I_ERI_K3x2y2z_S_S_S_ab+ABY*I_ERI_I3xy2z_S_S_S_ab;
  Double I_ERI_I2x4y_Py_S_S_ab = I_ERI_K2x5y_S_S_S_ab+ABY*I_ERI_I2x4y_S_S_S_ab;
  Double I_ERI_I2x3yz_Py_S_S_ab = I_ERI_K2x4yz_S_S_S_ab+ABY*I_ERI_I2x3yz_S_S_S_ab;
  Double I_ERI_I2x2y2z_Py_S_S_ab = I_ERI_K2x3y2z_S_S_S_ab+ABY*I_ERI_I2x2y2z_S_S_S_ab;
  Double I_ERI_I2xy3z_Py_S_S_ab = I_ERI_K2x2y3z_S_S_S_ab+ABY*I_ERI_I2xy3z_S_S_S_ab;
  Double I_ERI_Ix5y_Py_S_S_ab = I_ERI_Kx6y_S_S_S_ab+ABY*I_ERI_Ix5y_S_S_S_ab;
  Double I_ERI_Ix4yz_Py_S_S_ab = I_ERI_Kx5yz_S_S_S_ab+ABY*I_ERI_Ix4yz_S_S_S_ab;
  Double I_ERI_Ix3y2z_Py_S_S_ab = I_ERI_Kx4y2z_S_S_S_ab+ABY*I_ERI_Ix3y2z_S_S_S_ab;
  Double I_ERI_Ix2y3z_Py_S_S_ab = I_ERI_Kx3y3z_S_S_S_ab+ABY*I_ERI_Ix2y3z_S_S_S_ab;
  Double I_ERI_Ixy4z_Py_S_S_ab = I_ERI_Kx2y4z_S_S_S_ab+ABY*I_ERI_Ixy4z_S_S_S_ab;
  Double I_ERI_I6y_Py_S_S_ab = I_ERI_K7y_S_S_S_ab+ABY*I_ERI_I6y_S_S_S_ab;
  Double I_ERI_I5yz_Py_S_S_ab = I_ERI_K6yz_S_S_S_ab+ABY*I_ERI_I5yz_S_S_S_ab;
  Double I_ERI_I4y2z_Py_S_S_ab = I_ERI_K5y2z_S_S_S_ab+ABY*I_ERI_I4y2z_S_S_S_ab;
  Double I_ERI_I3y3z_Py_S_S_ab = I_ERI_K4y3z_S_S_S_ab+ABY*I_ERI_I3y3z_S_S_S_ab;
  Double I_ERI_I2y4z_Py_S_S_ab = I_ERI_K3y4z_S_S_S_ab+ABY*I_ERI_I2y4z_S_S_S_ab;
  Double I_ERI_Iy5z_Py_S_S_ab = I_ERI_K2y5z_S_S_S_ab+ABY*I_ERI_Iy5z_S_S_S_ab;
  Double I_ERI_I5xz_Pz_S_S_ab = I_ERI_K5x2z_S_S_S_ab+ABZ*I_ERI_I5xz_S_S_S_ab;
  Double I_ERI_I4xyz_Pz_S_S_ab = I_ERI_K4xy2z_S_S_S_ab+ABZ*I_ERI_I4xyz_S_S_S_ab;
  Double I_ERI_I4x2z_Pz_S_S_ab = I_ERI_K4x3z_S_S_S_ab+ABZ*I_ERI_I4x2z_S_S_S_ab;
  Double I_ERI_I3x2yz_Pz_S_S_ab = I_ERI_K3x2y2z_S_S_S_ab+ABZ*I_ERI_I3x2yz_S_S_S_ab;
  Double I_ERI_I3xy2z_Pz_S_S_ab = I_ERI_K3xy3z_S_S_S_ab+ABZ*I_ERI_I3xy2z_S_S_S_ab;
  Double I_ERI_I3x3z_Pz_S_S_ab = I_ERI_K3x4z_S_S_S_ab+ABZ*I_ERI_I3x3z_S_S_S_ab;
  Double I_ERI_I2x3yz_Pz_S_S_ab = I_ERI_K2x3y2z_S_S_S_ab+ABZ*I_ERI_I2x3yz_S_S_S_ab;
  Double I_ERI_I2x2y2z_Pz_S_S_ab = I_ERI_K2x2y3z_S_S_S_ab+ABZ*I_ERI_I2x2y2z_S_S_S_ab;
  Double I_ERI_I2xy3z_Pz_S_S_ab = I_ERI_K2xy4z_S_S_S_ab+ABZ*I_ERI_I2xy3z_S_S_S_ab;
  Double I_ERI_I2x4z_Pz_S_S_ab = I_ERI_K2x5z_S_S_S_ab+ABZ*I_ERI_I2x4z_S_S_S_ab;
  Double I_ERI_Ix4yz_Pz_S_S_ab = I_ERI_Kx4y2z_S_S_S_ab+ABZ*I_ERI_Ix4yz_S_S_S_ab;
  Double I_ERI_Ix3y2z_Pz_S_S_ab = I_ERI_Kx3y3z_S_S_S_ab+ABZ*I_ERI_Ix3y2z_S_S_S_ab;
  Double I_ERI_Ix2y3z_Pz_S_S_ab = I_ERI_Kx2y4z_S_S_S_ab+ABZ*I_ERI_Ix2y3z_S_S_S_ab;
  Double I_ERI_Ixy4z_Pz_S_S_ab = I_ERI_Kxy5z_S_S_S_ab+ABZ*I_ERI_Ixy4z_S_S_S_ab;
  Double I_ERI_Ix5z_Pz_S_S_ab = I_ERI_Kx6z_S_S_S_ab+ABZ*I_ERI_Ix5z_S_S_S_ab;
  Double I_ERI_I5yz_Pz_S_S_ab = I_ERI_K5y2z_S_S_S_ab+ABZ*I_ERI_I5yz_S_S_S_ab;
  Double I_ERI_I4y2z_Pz_S_S_ab = I_ERI_K4y3z_S_S_S_ab+ABZ*I_ERI_I4y2z_S_S_S_ab;
  Double I_ERI_I3y3z_Pz_S_S_ab = I_ERI_K3y4z_S_S_S_ab+ABZ*I_ERI_I3y3z_S_S_S_ab;
  Double I_ERI_I2y4z_Pz_S_S_ab = I_ERI_K2y5z_S_S_S_ab+ABZ*I_ERI_I2y4z_S_S_S_ab;
  Double I_ERI_Iy5z_Pz_S_S_ab = I_ERI_Ky6z_S_S_S_ab+ABZ*I_ERI_Iy5z_S_S_S_ab;
  Double I_ERI_I6z_Pz_S_S_ab = I_ERI_K7z_S_S_S_ab+ABZ*I_ERI_I6z_S_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_H_D_S_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 48 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_H_P_S_S_ab
   ************************************************************/
  Double I_ERI_H5x_D2x_S_S_ab = I_ERI_I6x_Px_S_S_ab+ABX*I_ERI_H5x_Px_S_S_ab;
  Double I_ERI_H4xy_D2x_S_S_ab = I_ERI_I5xy_Px_S_S_ab+ABX*I_ERI_H4xy_Px_S_S_ab;
  Double I_ERI_H4xz_D2x_S_S_ab = I_ERI_I5xz_Px_S_S_ab+ABX*I_ERI_H4xz_Px_S_S_ab;
  Double I_ERI_H3x2y_D2x_S_S_ab = I_ERI_I4x2y_Px_S_S_ab+ABX*I_ERI_H3x2y_Px_S_S_ab;
  Double I_ERI_H3xyz_D2x_S_S_ab = I_ERI_I4xyz_Px_S_S_ab+ABX*I_ERI_H3xyz_Px_S_S_ab;
  Double I_ERI_H3x2z_D2x_S_S_ab = I_ERI_I4x2z_Px_S_S_ab+ABX*I_ERI_H3x2z_Px_S_S_ab;
  Double I_ERI_H2x3y_D2x_S_S_ab = I_ERI_I3x3y_Px_S_S_ab+ABX*I_ERI_H2x3y_Px_S_S_ab;
  Double I_ERI_H2x2yz_D2x_S_S_ab = I_ERI_I3x2yz_Px_S_S_ab+ABX*I_ERI_H2x2yz_Px_S_S_ab;
  Double I_ERI_H2xy2z_D2x_S_S_ab = I_ERI_I3xy2z_Px_S_S_ab+ABX*I_ERI_H2xy2z_Px_S_S_ab;
  Double I_ERI_H2x3z_D2x_S_S_ab = I_ERI_I3x3z_Px_S_S_ab+ABX*I_ERI_H2x3z_Px_S_S_ab;
  Double I_ERI_Hx4y_D2x_S_S_ab = I_ERI_I2x4y_Px_S_S_ab+ABX*I_ERI_Hx4y_Px_S_S_ab;
  Double I_ERI_Hx3yz_D2x_S_S_ab = I_ERI_I2x3yz_Px_S_S_ab+ABX*I_ERI_Hx3yz_Px_S_S_ab;
  Double I_ERI_Hx2y2z_D2x_S_S_ab = I_ERI_I2x2y2z_Px_S_S_ab+ABX*I_ERI_Hx2y2z_Px_S_S_ab;
  Double I_ERI_Hxy3z_D2x_S_S_ab = I_ERI_I2xy3z_Px_S_S_ab+ABX*I_ERI_Hxy3z_Px_S_S_ab;
  Double I_ERI_Hx4z_D2x_S_S_ab = I_ERI_I2x4z_Px_S_S_ab+ABX*I_ERI_Hx4z_Px_S_S_ab;
  Double I_ERI_H5y_D2x_S_S_ab = I_ERI_Ix5y_Px_S_S_ab+ABX*I_ERI_H5y_Px_S_S_ab;
  Double I_ERI_H4yz_D2x_S_S_ab = I_ERI_Ix4yz_Px_S_S_ab+ABX*I_ERI_H4yz_Px_S_S_ab;
  Double I_ERI_H3y2z_D2x_S_S_ab = I_ERI_Ix3y2z_Px_S_S_ab+ABX*I_ERI_H3y2z_Px_S_S_ab;
  Double I_ERI_H2y3z_D2x_S_S_ab = I_ERI_Ix2y3z_Px_S_S_ab+ABX*I_ERI_H2y3z_Px_S_S_ab;
  Double I_ERI_Hy4z_D2x_S_S_ab = I_ERI_Ixy4z_Px_S_S_ab+ABX*I_ERI_Hy4z_Px_S_S_ab;
  Double I_ERI_H5z_D2x_S_S_ab = I_ERI_Ix5z_Px_S_S_ab+ABX*I_ERI_H5z_Px_S_S_ab;
  Double I_ERI_H4xz_Dxy_S_S_ab = I_ERI_I4xyz_Px_S_S_ab+ABY*I_ERI_H4xz_Px_S_S_ab;
  Double I_ERI_H3xyz_Dxy_S_S_ab = I_ERI_I3x2yz_Px_S_S_ab+ABY*I_ERI_H3xyz_Px_S_S_ab;
  Double I_ERI_H3x2z_Dxy_S_S_ab = I_ERI_I3xy2z_Px_S_S_ab+ABY*I_ERI_H3x2z_Px_S_S_ab;
  Double I_ERI_H2x2yz_Dxy_S_S_ab = I_ERI_I2x3yz_Px_S_S_ab+ABY*I_ERI_H2x2yz_Px_S_S_ab;
  Double I_ERI_H2xy2z_Dxy_S_S_ab = I_ERI_I2x2y2z_Px_S_S_ab+ABY*I_ERI_H2xy2z_Px_S_S_ab;
  Double I_ERI_H2x3z_Dxy_S_S_ab = I_ERI_I2xy3z_Px_S_S_ab+ABY*I_ERI_H2x3z_Px_S_S_ab;
  Double I_ERI_Hx3yz_Dxy_S_S_ab = I_ERI_Ix4yz_Px_S_S_ab+ABY*I_ERI_Hx3yz_Px_S_S_ab;
  Double I_ERI_Hx2y2z_Dxy_S_S_ab = I_ERI_Ix3y2z_Px_S_S_ab+ABY*I_ERI_Hx2y2z_Px_S_S_ab;
  Double I_ERI_Hxy3z_Dxy_S_S_ab = I_ERI_Ix2y3z_Px_S_S_ab+ABY*I_ERI_Hxy3z_Px_S_S_ab;
  Double I_ERI_Hx4z_Dxy_S_S_ab = I_ERI_Ixy4z_Px_S_S_ab+ABY*I_ERI_Hx4z_Px_S_S_ab;
  Double I_ERI_H4yz_Dxy_S_S_ab = I_ERI_I5yz_Px_S_S_ab+ABY*I_ERI_H4yz_Px_S_S_ab;
  Double I_ERI_H3y2z_Dxy_S_S_ab = I_ERI_I4y2z_Px_S_S_ab+ABY*I_ERI_H3y2z_Px_S_S_ab;
  Double I_ERI_H2y3z_Dxy_S_S_ab = I_ERI_I3y3z_Px_S_S_ab+ABY*I_ERI_H2y3z_Px_S_S_ab;
  Double I_ERI_Hy4z_Dxy_S_S_ab = I_ERI_I2y4z_Px_S_S_ab+ABY*I_ERI_Hy4z_Px_S_S_ab;
  Double I_ERI_H5z_Dxy_S_S_ab = I_ERI_Iy5z_Px_S_S_ab+ABY*I_ERI_H5z_Px_S_S_ab;
  Double I_ERI_H5x_D2y_S_S_ab = I_ERI_I5xy_Py_S_S_ab+ABY*I_ERI_H5x_Py_S_S_ab;
  Double I_ERI_H4xy_D2y_S_S_ab = I_ERI_I4x2y_Py_S_S_ab+ABY*I_ERI_H4xy_Py_S_S_ab;
  Double I_ERI_H4xz_D2y_S_S_ab = I_ERI_I4xyz_Py_S_S_ab+ABY*I_ERI_H4xz_Py_S_S_ab;
  Double I_ERI_H3x2y_D2y_S_S_ab = I_ERI_I3x3y_Py_S_S_ab+ABY*I_ERI_H3x2y_Py_S_S_ab;
  Double I_ERI_H3xyz_D2y_S_S_ab = I_ERI_I3x2yz_Py_S_S_ab+ABY*I_ERI_H3xyz_Py_S_S_ab;
  Double I_ERI_H3x2z_D2y_S_S_ab = I_ERI_I3xy2z_Py_S_S_ab+ABY*I_ERI_H3x2z_Py_S_S_ab;
  Double I_ERI_H2x3y_D2y_S_S_ab = I_ERI_I2x4y_Py_S_S_ab+ABY*I_ERI_H2x3y_Py_S_S_ab;
  Double I_ERI_H2x2yz_D2y_S_S_ab = I_ERI_I2x3yz_Py_S_S_ab+ABY*I_ERI_H2x2yz_Py_S_S_ab;
  Double I_ERI_H2xy2z_D2y_S_S_ab = I_ERI_I2x2y2z_Py_S_S_ab+ABY*I_ERI_H2xy2z_Py_S_S_ab;
  Double I_ERI_H2x3z_D2y_S_S_ab = I_ERI_I2xy3z_Py_S_S_ab+ABY*I_ERI_H2x3z_Py_S_S_ab;
  Double I_ERI_Hx4y_D2y_S_S_ab = I_ERI_Ix5y_Py_S_S_ab+ABY*I_ERI_Hx4y_Py_S_S_ab;
  Double I_ERI_Hx3yz_D2y_S_S_ab = I_ERI_Ix4yz_Py_S_S_ab+ABY*I_ERI_Hx3yz_Py_S_S_ab;
  Double I_ERI_Hx2y2z_D2y_S_S_ab = I_ERI_Ix3y2z_Py_S_S_ab+ABY*I_ERI_Hx2y2z_Py_S_S_ab;
  Double I_ERI_Hxy3z_D2y_S_S_ab = I_ERI_Ix2y3z_Py_S_S_ab+ABY*I_ERI_Hxy3z_Py_S_S_ab;
  Double I_ERI_Hx4z_D2y_S_S_ab = I_ERI_Ixy4z_Py_S_S_ab+ABY*I_ERI_Hx4z_Py_S_S_ab;
  Double I_ERI_H5y_D2y_S_S_ab = I_ERI_I6y_Py_S_S_ab+ABY*I_ERI_H5y_Py_S_S_ab;
  Double I_ERI_H4yz_D2y_S_S_ab = I_ERI_I5yz_Py_S_S_ab+ABY*I_ERI_H4yz_Py_S_S_ab;
  Double I_ERI_H3y2z_D2y_S_S_ab = I_ERI_I4y2z_Py_S_S_ab+ABY*I_ERI_H3y2z_Py_S_S_ab;
  Double I_ERI_H2y3z_D2y_S_S_ab = I_ERI_I3y3z_Py_S_S_ab+ABY*I_ERI_H2y3z_Py_S_S_ab;
  Double I_ERI_Hy4z_D2y_S_S_ab = I_ERI_I2y4z_Py_S_S_ab+ABY*I_ERI_Hy4z_Py_S_S_ab;
  Double I_ERI_H5z_D2y_S_S_ab = I_ERI_Iy5z_Py_S_S_ab+ABY*I_ERI_H5z_Py_S_S_ab;
  Double I_ERI_H5x_D2z_S_S_ab = I_ERI_I5xz_Pz_S_S_ab+ABZ*I_ERI_H5x_Pz_S_S_ab;
  Double I_ERI_H4xy_D2z_S_S_ab = I_ERI_I4xyz_Pz_S_S_ab+ABZ*I_ERI_H4xy_Pz_S_S_ab;
  Double I_ERI_H4xz_D2z_S_S_ab = I_ERI_I4x2z_Pz_S_S_ab+ABZ*I_ERI_H4xz_Pz_S_S_ab;
  Double I_ERI_H3x2y_D2z_S_S_ab = I_ERI_I3x2yz_Pz_S_S_ab+ABZ*I_ERI_H3x2y_Pz_S_S_ab;
  Double I_ERI_H3xyz_D2z_S_S_ab = I_ERI_I3xy2z_Pz_S_S_ab+ABZ*I_ERI_H3xyz_Pz_S_S_ab;
  Double I_ERI_H3x2z_D2z_S_S_ab = I_ERI_I3x3z_Pz_S_S_ab+ABZ*I_ERI_H3x2z_Pz_S_S_ab;
  Double I_ERI_H2x3y_D2z_S_S_ab = I_ERI_I2x3yz_Pz_S_S_ab+ABZ*I_ERI_H2x3y_Pz_S_S_ab;
  Double I_ERI_H2x2yz_D2z_S_S_ab = I_ERI_I2x2y2z_Pz_S_S_ab+ABZ*I_ERI_H2x2yz_Pz_S_S_ab;
  Double I_ERI_H2xy2z_D2z_S_S_ab = I_ERI_I2xy3z_Pz_S_S_ab+ABZ*I_ERI_H2xy2z_Pz_S_S_ab;
  Double I_ERI_H2x3z_D2z_S_S_ab = I_ERI_I2x4z_Pz_S_S_ab+ABZ*I_ERI_H2x3z_Pz_S_S_ab;
  Double I_ERI_Hx4y_D2z_S_S_ab = I_ERI_Ix4yz_Pz_S_S_ab+ABZ*I_ERI_Hx4y_Pz_S_S_ab;
  Double I_ERI_Hx3yz_D2z_S_S_ab = I_ERI_Ix3y2z_Pz_S_S_ab+ABZ*I_ERI_Hx3yz_Pz_S_S_ab;
  Double I_ERI_Hx2y2z_D2z_S_S_ab = I_ERI_Ix2y3z_Pz_S_S_ab+ABZ*I_ERI_Hx2y2z_Pz_S_S_ab;
  Double I_ERI_Hxy3z_D2z_S_S_ab = I_ERI_Ixy4z_Pz_S_S_ab+ABZ*I_ERI_Hxy3z_Pz_S_S_ab;
  Double I_ERI_Hx4z_D2z_S_S_ab = I_ERI_Ix5z_Pz_S_S_ab+ABZ*I_ERI_Hx4z_Pz_S_S_ab;
  Double I_ERI_H5y_D2z_S_S_ab = I_ERI_I5yz_Pz_S_S_ab+ABZ*I_ERI_H5y_Pz_S_S_ab;
  Double I_ERI_H4yz_D2z_S_S_ab = I_ERI_I4y2z_Pz_S_S_ab+ABZ*I_ERI_H4yz_Pz_S_S_ab;
  Double I_ERI_H3y2z_D2z_S_S_ab = I_ERI_I3y3z_Pz_S_S_ab+ABZ*I_ERI_H3y2z_Pz_S_S_ab;
  Double I_ERI_H2y3z_D2z_S_S_ab = I_ERI_I2y4z_Pz_S_S_ab+ABZ*I_ERI_H2y3z_Pz_S_S_ab;
  Double I_ERI_Hy4z_D2z_S_S_ab = I_ERI_Iy5z_Pz_S_S_ab+ABZ*I_ERI_Hy4z_Pz_S_S_ab;
  Double I_ERI_H5z_D2z_S_S_ab = I_ERI_I6z_Pz_S_S_ab+ABZ*I_ERI_H5z_Pz_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_G_F_S_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_D_S_S_ab
   ************************************************************/
  Double I_ERI_G4x_F3x_S_S_ab = I_ERI_H5x_D2x_S_S_ab+ABX*I_ERI_G4x_D2x_S_S_ab;
  Double I_ERI_G3xy_F3x_S_S_ab = I_ERI_H4xy_D2x_S_S_ab+ABX*I_ERI_G3xy_D2x_S_S_ab;
  Double I_ERI_G3xz_F3x_S_S_ab = I_ERI_H4xz_D2x_S_S_ab+ABX*I_ERI_G3xz_D2x_S_S_ab;
  Double I_ERI_G2x2y_F3x_S_S_ab = I_ERI_H3x2y_D2x_S_S_ab+ABX*I_ERI_G2x2y_D2x_S_S_ab;
  Double I_ERI_G2xyz_F3x_S_S_ab = I_ERI_H3xyz_D2x_S_S_ab+ABX*I_ERI_G2xyz_D2x_S_S_ab;
  Double I_ERI_G2x2z_F3x_S_S_ab = I_ERI_H3x2z_D2x_S_S_ab+ABX*I_ERI_G2x2z_D2x_S_S_ab;
  Double I_ERI_Gx3y_F3x_S_S_ab = I_ERI_H2x3y_D2x_S_S_ab+ABX*I_ERI_Gx3y_D2x_S_S_ab;
  Double I_ERI_Gx2yz_F3x_S_S_ab = I_ERI_H2x2yz_D2x_S_S_ab+ABX*I_ERI_Gx2yz_D2x_S_S_ab;
  Double I_ERI_Gxy2z_F3x_S_S_ab = I_ERI_H2xy2z_D2x_S_S_ab+ABX*I_ERI_Gxy2z_D2x_S_S_ab;
  Double I_ERI_Gx3z_F3x_S_S_ab = I_ERI_H2x3z_D2x_S_S_ab+ABX*I_ERI_Gx3z_D2x_S_S_ab;
  Double I_ERI_G4y_F3x_S_S_ab = I_ERI_Hx4y_D2x_S_S_ab+ABX*I_ERI_G4y_D2x_S_S_ab;
  Double I_ERI_G3yz_F3x_S_S_ab = I_ERI_Hx3yz_D2x_S_S_ab+ABX*I_ERI_G3yz_D2x_S_S_ab;
  Double I_ERI_G2y2z_F3x_S_S_ab = I_ERI_Hx2y2z_D2x_S_S_ab+ABX*I_ERI_G2y2z_D2x_S_S_ab;
  Double I_ERI_Gy3z_F3x_S_S_ab = I_ERI_Hxy3z_D2x_S_S_ab+ABX*I_ERI_Gy3z_D2x_S_S_ab;
  Double I_ERI_G4z_F3x_S_S_ab = I_ERI_Hx4z_D2x_S_S_ab+ABX*I_ERI_G4z_D2x_S_S_ab;
  Double I_ERI_G4x_F2xy_S_S_ab = I_ERI_H4xy_D2x_S_S_ab+ABY*I_ERI_G4x_D2x_S_S_ab;
  Double I_ERI_G3xy_F2xy_S_S_ab = I_ERI_H3x2y_D2x_S_S_ab+ABY*I_ERI_G3xy_D2x_S_S_ab;
  Double I_ERI_G3xz_F2xy_S_S_ab = I_ERI_H3xyz_D2x_S_S_ab+ABY*I_ERI_G3xz_D2x_S_S_ab;
  Double I_ERI_G2x2y_F2xy_S_S_ab = I_ERI_H2x3y_D2x_S_S_ab+ABY*I_ERI_G2x2y_D2x_S_S_ab;
  Double I_ERI_G2xyz_F2xy_S_S_ab = I_ERI_H2x2yz_D2x_S_S_ab+ABY*I_ERI_G2xyz_D2x_S_S_ab;
  Double I_ERI_G2x2z_F2xy_S_S_ab = I_ERI_H2xy2z_D2x_S_S_ab+ABY*I_ERI_G2x2z_D2x_S_S_ab;
  Double I_ERI_Gx3y_F2xy_S_S_ab = I_ERI_Hx4y_D2x_S_S_ab+ABY*I_ERI_Gx3y_D2x_S_S_ab;
  Double I_ERI_Gx2yz_F2xy_S_S_ab = I_ERI_Hx3yz_D2x_S_S_ab+ABY*I_ERI_Gx2yz_D2x_S_S_ab;
  Double I_ERI_Gxy2z_F2xy_S_S_ab = I_ERI_Hx2y2z_D2x_S_S_ab+ABY*I_ERI_Gxy2z_D2x_S_S_ab;
  Double I_ERI_Gx3z_F2xy_S_S_ab = I_ERI_Hxy3z_D2x_S_S_ab+ABY*I_ERI_Gx3z_D2x_S_S_ab;
  Double I_ERI_G4y_F2xy_S_S_ab = I_ERI_H5y_D2x_S_S_ab+ABY*I_ERI_G4y_D2x_S_S_ab;
  Double I_ERI_G3yz_F2xy_S_S_ab = I_ERI_H4yz_D2x_S_S_ab+ABY*I_ERI_G3yz_D2x_S_S_ab;
  Double I_ERI_G2y2z_F2xy_S_S_ab = I_ERI_H3y2z_D2x_S_S_ab+ABY*I_ERI_G2y2z_D2x_S_S_ab;
  Double I_ERI_Gy3z_F2xy_S_S_ab = I_ERI_H2y3z_D2x_S_S_ab+ABY*I_ERI_Gy3z_D2x_S_S_ab;
  Double I_ERI_G4z_F2xy_S_S_ab = I_ERI_Hy4z_D2x_S_S_ab+ABY*I_ERI_G4z_D2x_S_S_ab;
  Double I_ERI_G4x_F2xz_S_S_ab = I_ERI_H4xz_D2x_S_S_ab+ABZ*I_ERI_G4x_D2x_S_S_ab;
  Double I_ERI_G3xy_F2xz_S_S_ab = I_ERI_H3xyz_D2x_S_S_ab+ABZ*I_ERI_G3xy_D2x_S_S_ab;
  Double I_ERI_G3xz_F2xz_S_S_ab = I_ERI_H3x2z_D2x_S_S_ab+ABZ*I_ERI_G3xz_D2x_S_S_ab;
  Double I_ERI_G2x2y_F2xz_S_S_ab = I_ERI_H2x2yz_D2x_S_S_ab+ABZ*I_ERI_G2x2y_D2x_S_S_ab;
  Double I_ERI_G2xyz_F2xz_S_S_ab = I_ERI_H2xy2z_D2x_S_S_ab+ABZ*I_ERI_G2xyz_D2x_S_S_ab;
  Double I_ERI_G2x2z_F2xz_S_S_ab = I_ERI_H2x3z_D2x_S_S_ab+ABZ*I_ERI_G2x2z_D2x_S_S_ab;
  Double I_ERI_Gx3y_F2xz_S_S_ab = I_ERI_Hx3yz_D2x_S_S_ab+ABZ*I_ERI_Gx3y_D2x_S_S_ab;
  Double I_ERI_Gx2yz_F2xz_S_S_ab = I_ERI_Hx2y2z_D2x_S_S_ab+ABZ*I_ERI_Gx2yz_D2x_S_S_ab;
  Double I_ERI_Gxy2z_F2xz_S_S_ab = I_ERI_Hxy3z_D2x_S_S_ab+ABZ*I_ERI_Gxy2z_D2x_S_S_ab;
  Double I_ERI_Gx3z_F2xz_S_S_ab = I_ERI_Hx4z_D2x_S_S_ab+ABZ*I_ERI_Gx3z_D2x_S_S_ab;
  Double I_ERI_G4y_F2xz_S_S_ab = I_ERI_H4yz_D2x_S_S_ab+ABZ*I_ERI_G4y_D2x_S_S_ab;
  Double I_ERI_G3yz_F2xz_S_S_ab = I_ERI_H3y2z_D2x_S_S_ab+ABZ*I_ERI_G3yz_D2x_S_S_ab;
  Double I_ERI_G2y2z_F2xz_S_S_ab = I_ERI_H2y3z_D2x_S_S_ab+ABZ*I_ERI_G2y2z_D2x_S_S_ab;
  Double I_ERI_Gy3z_F2xz_S_S_ab = I_ERI_Hy4z_D2x_S_S_ab+ABZ*I_ERI_Gy3z_D2x_S_S_ab;
  Double I_ERI_G4z_F2xz_S_S_ab = I_ERI_H5z_D2x_S_S_ab+ABZ*I_ERI_G4z_D2x_S_S_ab;
  Double I_ERI_G4x_Fx2y_S_S_ab = I_ERI_H5x_D2y_S_S_ab+ABX*I_ERI_G4x_D2y_S_S_ab;
  Double I_ERI_G3xy_Fx2y_S_S_ab = I_ERI_H4xy_D2y_S_S_ab+ABX*I_ERI_G3xy_D2y_S_S_ab;
  Double I_ERI_G3xz_Fx2y_S_S_ab = I_ERI_H4xz_D2y_S_S_ab+ABX*I_ERI_G3xz_D2y_S_S_ab;
  Double I_ERI_G2x2y_Fx2y_S_S_ab = I_ERI_H3x2y_D2y_S_S_ab+ABX*I_ERI_G2x2y_D2y_S_S_ab;
  Double I_ERI_G2xyz_Fx2y_S_S_ab = I_ERI_H3xyz_D2y_S_S_ab+ABX*I_ERI_G2xyz_D2y_S_S_ab;
  Double I_ERI_G2x2z_Fx2y_S_S_ab = I_ERI_H3x2z_D2y_S_S_ab+ABX*I_ERI_G2x2z_D2y_S_S_ab;
  Double I_ERI_Gx3y_Fx2y_S_S_ab = I_ERI_H2x3y_D2y_S_S_ab+ABX*I_ERI_Gx3y_D2y_S_S_ab;
  Double I_ERI_Gx2yz_Fx2y_S_S_ab = I_ERI_H2x2yz_D2y_S_S_ab+ABX*I_ERI_Gx2yz_D2y_S_S_ab;
  Double I_ERI_Gxy2z_Fx2y_S_S_ab = I_ERI_H2xy2z_D2y_S_S_ab+ABX*I_ERI_Gxy2z_D2y_S_S_ab;
  Double I_ERI_Gx3z_Fx2y_S_S_ab = I_ERI_H2x3z_D2y_S_S_ab+ABX*I_ERI_Gx3z_D2y_S_S_ab;
  Double I_ERI_G4y_Fx2y_S_S_ab = I_ERI_Hx4y_D2y_S_S_ab+ABX*I_ERI_G4y_D2y_S_S_ab;
  Double I_ERI_G3yz_Fx2y_S_S_ab = I_ERI_Hx3yz_D2y_S_S_ab+ABX*I_ERI_G3yz_D2y_S_S_ab;
  Double I_ERI_G2y2z_Fx2y_S_S_ab = I_ERI_Hx2y2z_D2y_S_S_ab+ABX*I_ERI_G2y2z_D2y_S_S_ab;
  Double I_ERI_Gy3z_Fx2y_S_S_ab = I_ERI_Hxy3z_D2y_S_S_ab+ABX*I_ERI_Gy3z_D2y_S_S_ab;
  Double I_ERI_G4z_Fx2y_S_S_ab = I_ERI_Hx4z_D2y_S_S_ab+ABX*I_ERI_G4z_D2y_S_S_ab;
  Double I_ERI_G4x_Fxyz_S_S_ab = I_ERI_H4xz_Dxy_S_S_ab+ABZ*I_ERI_G4x_Dxy_S_S_ab;
  Double I_ERI_G3xy_Fxyz_S_S_ab = I_ERI_H3xyz_Dxy_S_S_ab+ABZ*I_ERI_G3xy_Dxy_S_S_ab;
  Double I_ERI_G3xz_Fxyz_S_S_ab = I_ERI_H3x2z_Dxy_S_S_ab+ABZ*I_ERI_G3xz_Dxy_S_S_ab;
  Double I_ERI_G2x2y_Fxyz_S_S_ab = I_ERI_H2x2yz_Dxy_S_S_ab+ABZ*I_ERI_G2x2y_Dxy_S_S_ab;
  Double I_ERI_G2xyz_Fxyz_S_S_ab = I_ERI_H2xy2z_Dxy_S_S_ab+ABZ*I_ERI_G2xyz_Dxy_S_S_ab;
  Double I_ERI_G2x2z_Fxyz_S_S_ab = I_ERI_H2x3z_Dxy_S_S_ab+ABZ*I_ERI_G2x2z_Dxy_S_S_ab;
  Double I_ERI_Gx3y_Fxyz_S_S_ab = I_ERI_Hx3yz_Dxy_S_S_ab+ABZ*I_ERI_Gx3y_Dxy_S_S_ab;
  Double I_ERI_Gx2yz_Fxyz_S_S_ab = I_ERI_Hx2y2z_Dxy_S_S_ab+ABZ*I_ERI_Gx2yz_Dxy_S_S_ab;
  Double I_ERI_Gxy2z_Fxyz_S_S_ab = I_ERI_Hxy3z_Dxy_S_S_ab+ABZ*I_ERI_Gxy2z_Dxy_S_S_ab;
  Double I_ERI_Gx3z_Fxyz_S_S_ab = I_ERI_Hx4z_Dxy_S_S_ab+ABZ*I_ERI_Gx3z_Dxy_S_S_ab;
  Double I_ERI_G4y_Fxyz_S_S_ab = I_ERI_H4yz_Dxy_S_S_ab+ABZ*I_ERI_G4y_Dxy_S_S_ab;
  Double I_ERI_G3yz_Fxyz_S_S_ab = I_ERI_H3y2z_Dxy_S_S_ab+ABZ*I_ERI_G3yz_Dxy_S_S_ab;
  Double I_ERI_G2y2z_Fxyz_S_S_ab = I_ERI_H2y3z_Dxy_S_S_ab+ABZ*I_ERI_G2y2z_Dxy_S_S_ab;
  Double I_ERI_Gy3z_Fxyz_S_S_ab = I_ERI_Hy4z_Dxy_S_S_ab+ABZ*I_ERI_Gy3z_Dxy_S_S_ab;
  Double I_ERI_G4z_Fxyz_S_S_ab = I_ERI_H5z_Dxy_S_S_ab+ABZ*I_ERI_G4z_Dxy_S_S_ab;
  Double I_ERI_G4x_Fx2z_S_S_ab = I_ERI_H5x_D2z_S_S_ab+ABX*I_ERI_G4x_D2z_S_S_ab;
  Double I_ERI_G3xy_Fx2z_S_S_ab = I_ERI_H4xy_D2z_S_S_ab+ABX*I_ERI_G3xy_D2z_S_S_ab;
  Double I_ERI_G3xz_Fx2z_S_S_ab = I_ERI_H4xz_D2z_S_S_ab+ABX*I_ERI_G3xz_D2z_S_S_ab;
  Double I_ERI_G2x2y_Fx2z_S_S_ab = I_ERI_H3x2y_D2z_S_S_ab+ABX*I_ERI_G2x2y_D2z_S_S_ab;
  Double I_ERI_G2xyz_Fx2z_S_S_ab = I_ERI_H3xyz_D2z_S_S_ab+ABX*I_ERI_G2xyz_D2z_S_S_ab;
  Double I_ERI_G2x2z_Fx2z_S_S_ab = I_ERI_H3x2z_D2z_S_S_ab+ABX*I_ERI_G2x2z_D2z_S_S_ab;
  Double I_ERI_Gx3y_Fx2z_S_S_ab = I_ERI_H2x3y_D2z_S_S_ab+ABX*I_ERI_Gx3y_D2z_S_S_ab;
  Double I_ERI_Gx2yz_Fx2z_S_S_ab = I_ERI_H2x2yz_D2z_S_S_ab+ABX*I_ERI_Gx2yz_D2z_S_S_ab;
  Double I_ERI_Gxy2z_Fx2z_S_S_ab = I_ERI_H2xy2z_D2z_S_S_ab+ABX*I_ERI_Gxy2z_D2z_S_S_ab;
  Double I_ERI_Gx3z_Fx2z_S_S_ab = I_ERI_H2x3z_D2z_S_S_ab+ABX*I_ERI_Gx3z_D2z_S_S_ab;
  Double I_ERI_G4y_Fx2z_S_S_ab = I_ERI_Hx4y_D2z_S_S_ab+ABX*I_ERI_G4y_D2z_S_S_ab;
  Double I_ERI_G3yz_Fx2z_S_S_ab = I_ERI_Hx3yz_D2z_S_S_ab+ABX*I_ERI_G3yz_D2z_S_S_ab;
  Double I_ERI_G2y2z_Fx2z_S_S_ab = I_ERI_Hx2y2z_D2z_S_S_ab+ABX*I_ERI_G2y2z_D2z_S_S_ab;
  Double I_ERI_Gy3z_Fx2z_S_S_ab = I_ERI_Hxy3z_D2z_S_S_ab+ABX*I_ERI_Gy3z_D2z_S_S_ab;
  Double I_ERI_G4z_Fx2z_S_S_ab = I_ERI_Hx4z_D2z_S_S_ab+ABX*I_ERI_G4z_D2z_S_S_ab;
  Double I_ERI_G4x_F3y_S_S_ab = I_ERI_H4xy_D2y_S_S_ab+ABY*I_ERI_G4x_D2y_S_S_ab;
  Double I_ERI_G3xy_F3y_S_S_ab = I_ERI_H3x2y_D2y_S_S_ab+ABY*I_ERI_G3xy_D2y_S_S_ab;
  Double I_ERI_G3xz_F3y_S_S_ab = I_ERI_H3xyz_D2y_S_S_ab+ABY*I_ERI_G3xz_D2y_S_S_ab;
  Double I_ERI_G2x2y_F3y_S_S_ab = I_ERI_H2x3y_D2y_S_S_ab+ABY*I_ERI_G2x2y_D2y_S_S_ab;
  Double I_ERI_G2xyz_F3y_S_S_ab = I_ERI_H2x2yz_D2y_S_S_ab+ABY*I_ERI_G2xyz_D2y_S_S_ab;
  Double I_ERI_G2x2z_F3y_S_S_ab = I_ERI_H2xy2z_D2y_S_S_ab+ABY*I_ERI_G2x2z_D2y_S_S_ab;
  Double I_ERI_Gx3y_F3y_S_S_ab = I_ERI_Hx4y_D2y_S_S_ab+ABY*I_ERI_Gx3y_D2y_S_S_ab;
  Double I_ERI_Gx2yz_F3y_S_S_ab = I_ERI_Hx3yz_D2y_S_S_ab+ABY*I_ERI_Gx2yz_D2y_S_S_ab;
  Double I_ERI_Gxy2z_F3y_S_S_ab = I_ERI_Hx2y2z_D2y_S_S_ab+ABY*I_ERI_Gxy2z_D2y_S_S_ab;
  Double I_ERI_Gx3z_F3y_S_S_ab = I_ERI_Hxy3z_D2y_S_S_ab+ABY*I_ERI_Gx3z_D2y_S_S_ab;
  Double I_ERI_G4y_F3y_S_S_ab = I_ERI_H5y_D2y_S_S_ab+ABY*I_ERI_G4y_D2y_S_S_ab;
  Double I_ERI_G3yz_F3y_S_S_ab = I_ERI_H4yz_D2y_S_S_ab+ABY*I_ERI_G3yz_D2y_S_S_ab;
  Double I_ERI_G2y2z_F3y_S_S_ab = I_ERI_H3y2z_D2y_S_S_ab+ABY*I_ERI_G2y2z_D2y_S_S_ab;
  Double I_ERI_Gy3z_F3y_S_S_ab = I_ERI_H2y3z_D2y_S_S_ab+ABY*I_ERI_Gy3z_D2y_S_S_ab;
  Double I_ERI_G4z_F3y_S_S_ab = I_ERI_Hy4z_D2y_S_S_ab+ABY*I_ERI_G4z_D2y_S_S_ab;
  Double I_ERI_G4x_F2yz_S_S_ab = I_ERI_H4xz_D2y_S_S_ab+ABZ*I_ERI_G4x_D2y_S_S_ab;
  Double I_ERI_G3xy_F2yz_S_S_ab = I_ERI_H3xyz_D2y_S_S_ab+ABZ*I_ERI_G3xy_D2y_S_S_ab;
  Double I_ERI_G3xz_F2yz_S_S_ab = I_ERI_H3x2z_D2y_S_S_ab+ABZ*I_ERI_G3xz_D2y_S_S_ab;
  Double I_ERI_G2x2y_F2yz_S_S_ab = I_ERI_H2x2yz_D2y_S_S_ab+ABZ*I_ERI_G2x2y_D2y_S_S_ab;
  Double I_ERI_G2xyz_F2yz_S_S_ab = I_ERI_H2xy2z_D2y_S_S_ab+ABZ*I_ERI_G2xyz_D2y_S_S_ab;
  Double I_ERI_G2x2z_F2yz_S_S_ab = I_ERI_H2x3z_D2y_S_S_ab+ABZ*I_ERI_G2x2z_D2y_S_S_ab;
  Double I_ERI_Gx3y_F2yz_S_S_ab = I_ERI_Hx3yz_D2y_S_S_ab+ABZ*I_ERI_Gx3y_D2y_S_S_ab;
  Double I_ERI_Gx2yz_F2yz_S_S_ab = I_ERI_Hx2y2z_D2y_S_S_ab+ABZ*I_ERI_Gx2yz_D2y_S_S_ab;
  Double I_ERI_Gxy2z_F2yz_S_S_ab = I_ERI_Hxy3z_D2y_S_S_ab+ABZ*I_ERI_Gxy2z_D2y_S_S_ab;
  Double I_ERI_Gx3z_F2yz_S_S_ab = I_ERI_Hx4z_D2y_S_S_ab+ABZ*I_ERI_Gx3z_D2y_S_S_ab;
  Double I_ERI_G4y_F2yz_S_S_ab = I_ERI_H4yz_D2y_S_S_ab+ABZ*I_ERI_G4y_D2y_S_S_ab;
  Double I_ERI_G3yz_F2yz_S_S_ab = I_ERI_H3y2z_D2y_S_S_ab+ABZ*I_ERI_G3yz_D2y_S_S_ab;
  Double I_ERI_G2y2z_F2yz_S_S_ab = I_ERI_H2y3z_D2y_S_S_ab+ABZ*I_ERI_G2y2z_D2y_S_S_ab;
  Double I_ERI_Gy3z_F2yz_S_S_ab = I_ERI_Hy4z_D2y_S_S_ab+ABZ*I_ERI_Gy3z_D2y_S_S_ab;
  Double I_ERI_G4z_F2yz_S_S_ab = I_ERI_H5z_D2y_S_S_ab+ABZ*I_ERI_G4z_D2y_S_S_ab;
  Double I_ERI_G4x_Fy2z_S_S_ab = I_ERI_H4xy_D2z_S_S_ab+ABY*I_ERI_G4x_D2z_S_S_ab;
  Double I_ERI_G3xy_Fy2z_S_S_ab = I_ERI_H3x2y_D2z_S_S_ab+ABY*I_ERI_G3xy_D2z_S_S_ab;
  Double I_ERI_G3xz_Fy2z_S_S_ab = I_ERI_H3xyz_D2z_S_S_ab+ABY*I_ERI_G3xz_D2z_S_S_ab;
  Double I_ERI_G2x2y_Fy2z_S_S_ab = I_ERI_H2x3y_D2z_S_S_ab+ABY*I_ERI_G2x2y_D2z_S_S_ab;
  Double I_ERI_G2xyz_Fy2z_S_S_ab = I_ERI_H2x2yz_D2z_S_S_ab+ABY*I_ERI_G2xyz_D2z_S_S_ab;
  Double I_ERI_G2x2z_Fy2z_S_S_ab = I_ERI_H2xy2z_D2z_S_S_ab+ABY*I_ERI_G2x2z_D2z_S_S_ab;
  Double I_ERI_Gx3y_Fy2z_S_S_ab = I_ERI_Hx4y_D2z_S_S_ab+ABY*I_ERI_Gx3y_D2z_S_S_ab;
  Double I_ERI_Gx2yz_Fy2z_S_S_ab = I_ERI_Hx3yz_D2z_S_S_ab+ABY*I_ERI_Gx2yz_D2z_S_S_ab;
  Double I_ERI_Gxy2z_Fy2z_S_S_ab = I_ERI_Hx2y2z_D2z_S_S_ab+ABY*I_ERI_Gxy2z_D2z_S_S_ab;
  Double I_ERI_Gx3z_Fy2z_S_S_ab = I_ERI_Hxy3z_D2z_S_S_ab+ABY*I_ERI_Gx3z_D2z_S_S_ab;
  Double I_ERI_G4y_Fy2z_S_S_ab = I_ERI_H5y_D2z_S_S_ab+ABY*I_ERI_G4y_D2z_S_S_ab;
  Double I_ERI_G3yz_Fy2z_S_S_ab = I_ERI_H4yz_D2z_S_S_ab+ABY*I_ERI_G3yz_D2z_S_S_ab;
  Double I_ERI_G2y2z_Fy2z_S_S_ab = I_ERI_H3y2z_D2z_S_S_ab+ABY*I_ERI_G2y2z_D2z_S_S_ab;
  Double I_ERI_Gy3z_Fy2z_S_S_ab = I_ERI_H2y3z_D2z_S_S_ab+ABY*I_ERI_Gy3z_D2z_S_S_ab;
  Double I_ERI_G4z_Fy2z_S_S_ab = I_ERI_Hy4z_D2z_S_S_ab+ABY*I_ERI_G4z_D2z_S_S_ab;
  Double I_ERI_G4x_F3z_S_S_ab = I_ERI_H4xz_D2z_S_S_ab+ABZ*I_ERI_G4x_D2z_S_S_ab;
  Double I_ERI_G3xy_F3z_S_S_ab = I_ERI_H3xyz_D2z_S_S_ab+ABZ*I_ERI_G3xy_D2z_S_S_ab;
  Double I_ERI_G3xz_F3z_S_S_ab = I_ERI_H3x2z_D2z_S_S_ab+ABZ*I_ERI_G3xz_D2z_S_S_ab;
  Double I_ERI_G2x2y_F3z_S_S_ab = I_ERI_H2x2yz_D2z_S_S_ab+ABZ*I_ERI_G2x2y_D2z_S_S_ab;
  Double I_ERI_G2xyz_F3z_S_S_ab = I_ERI_H2xy2z_D2z_S_S_ab+ABZ*I_ERI_G2xyz_D2z_S_S_ab;
  Double I_ERI_G2x2z_F3z_S_S_ab = I_ERI_H2x3z_D2z_S_S_ab+ABZ*I_ERI_G2x2z_D2z_S_S_ab;
  Double I_ERI_Gx3y_F3z_S_S_ab = I_ERI_Hx3yz_D2z_S_S_ab+ABZ*I_ERI_Gx3y_D2z_S_S_ab;
  Double I_ERI_Gx2yz_F3z_S_S_ab = I_ERI_Hx2y2z_D2z_S_S_ab+ABZ*I_ERI_Gx2yz_D2z_S_S_ab;
  Double I_ERI_Gxy2z_F3z_S_S_ab = I_ERI_Hxy3z_D2z_S_S_ab+ABZ*I_ERI_Gxy2z_D2z_S_S_ab;
  Double I_ERI_Gx3z_F3z_S_S_ab = I_ERI_Hx4z_D2z_S_S_ab+ABZ*I_ERI_Gx3z_D2z_S_S_ab;
  Double I_ERI_G4y_F3z_S_S_ab = I_ERI_H4yz_D2z_S_S_ab+ABZ*I_ERI_G4y_D2z_S_S_ab;
  Double I_ERI_G3yz_F3z_S_S_ab = I_ERI_H3y2z_D2z_S_S_ab+ABZ*I_ERI_G3yz_D2z_S_S_ab;
  Double I_ERI_G2y2z_F3z_S_S_ab = I_ERI_H2y3z_D2z_S_S_ab+ABZ*I_ERI_G2y2z_D2z_S_S_ab;
  Double I_ERI_Gy3z_F3z_S_S_ab = I_ERI_Hy4z_D2z_S_S_ab+ABZ*I_ERI_Gy3z_D2z_S_S_ab;
  Double I_ERI_G4z_F3z_S_S_ab = I_ERI_H5z_D2z_S_S_ab+ABZ*I_ERI_G4z_D2z_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_P_S_ac
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_G_S_P_S_ac
   ************************************************************/
  Double I_ERI_G4x_Px_Px_S_ac = I_ERI_H5x_S_Px_S_ac+ABX*I_ERI_G4x_S_Px_S_ac;
  Double I_ERI_G3xy_Px_Px_S_ac = I_ERI_H4xy_S_Px_S_ac+ABX*I_ERI_G3xy_S_Px_S_ac;
  Double I_ERI_G3xz_Px_Px_S_ac = I_ERI_H4xz_S_Px_S_ac+ABX*I_ERI_G3xz_S_Px_S_ac;
  Double I_ERI_G2x2y_Px_Px_S_ac = I_ERI_H3x2y_S_Px_S_ac+ABX*I_ERI_G2x2y_S_Px_S_ac;
  Double I_ERI_G2xyz_Px_Px_S_ac = I_ERI_H3xyz_S_Px_S_ac+ABX*I_ERI_G2xyz_S_Px_S_ac;
  Double I_ERI_G2x2z_Px_Px_S_ac = I_ERI_H3x2z_S_Px_S_ac+ABX*I_ERI_G2x2z_S_Px_S_ac;
  Double I_ERI_Gx3y_Px_Px_S_ac = I_ERI_H2x3y_S_Px_S_ac+ABX*I_ERI_Gx3y_S_Px_S_ac;
  Double I_ERI_Gx2yz_Px_Px_S_ac = I_ERI_H2x2yz_S_Px_S_ac+ABX*I_ERI_Gx2yz_S_Px_S_ac;
  Double I_ERI_Gxy2z_Px_Px_S_ac = I_ERI_H2xy2z_S_Px_S_ac+ABX*I_ERI_Gxy2z_S_Px_S_ac;
  Double I_ERI_Gx3z_Px_Px_S_ac = I_ERI_H2x3z_S_Px_S_ac+ABX*I_ERI_Gx3z_S_Px_S_ac;
  Double I_ERI_G4y_Px_Px_S_ac = I_ERI_Hx4y_S_Px_S_ac+ABX*I_ERI_G4y_S_Px_S_ac;
  Double I_ERI_G3yz_Px_Px_S_ac = I_ERI_Hx3yz_S_Px_S_ac+ABX*I_ERI_G3yz_S_Px_S_ac;
  Double I_ERI_G2y2z_Px_Px_S_ac = I_ERI_Hx2y2z_S_Px_S_ac+ABX*I_ERI_G2y2z_S_Px_S_ac;
  Double I_ERI_Gy3z_Px_Px_S_ac = I_ERI_Hxy3z_S_Px_S_ac+ABX*I_ERI_Gy3z_S_Px_S_ac;
  Double I_ERI_G4z_Px_Px_S_ac = I_ERI_Hx4z_S_Px_S_ac+ABX*I_ERI_G4z_S_Px_S_ac;
  Double I_ERI_G4x_Py_Px_S_ac = I_ERI_H4xy_S_Px_S_ac+ABY*I_ERI_G4x_S_Px_S_ac;
  Double I_ERI_G3xy_Py_Px_S_ac = I_ERI_H3x2y_S_Px_S_ac+ABY*I_ERI_G3xy_S_Px_S_ac;
  Double I_ERI_G3xz_Py_Px_S_ac = I_ERI_H3xyz_S_Px_S_ac+ABY*I_ERI_G3xz_S_Px_S_ac;
  Double I_ERI_G2x2y_Py_Px_S_ac = I_ERI_H2x3y_S_Px_S_ac+ABY*I_ERI_G2x2y_S_Px_S_ac;
  Double I_ERI_G2xyz_Py_Px_S_ac = I_ERI_H2x2yz_S_Px_S_ac+ABY*I_ERI_G2xyz_S_Px_S_ac;
  Double I_ERI_G2x2z_Py_Px_S_ac = I_ERI_H2xy2z_S_Px_S_ac+ABY*I_ERI_G2x2z_S_Px_S_ac;
  Double I_ERI_Gx3y_Py_Px_S_ac = I_ERI_Hx4y_S_Px_S_ac+ABY*I_ERI_Gx3y_S_Px_S_ac;
  Double I_ERI_Gx2yz_Py_Px_S_ac = I_ERI_Hx3yz_S_Px_S_ac+ABY*I_ERI_Gx2yz_S_Px_S_ac;
  Double I_ERI_Gxy2z_Py_Px_S_ac = I_ERI_Hx2y2z_S_Px_S_ac+ABY*I_ERI_Gxy2z_S_Px_S_ac;
  Double I_ERI_Gx3z_Py_Px_S_ac = I_ERI_Hxy3z_S_Px_S_ac+ABY*I_ERI_Gx3z_S_Px_S_ac;
  Double I_ERI_G4y_Py_Px_S_ac = I_ERI_H5y_S_Px_S_ac+ABY*I_ERI_G4y_S_Px_S_ac;
  Double I_ERI_G3yz_Py_Px_S_ac = I_ERI_H4yz_S_Px_S_ac+ABY*I_ERI_G3yz_S_Px_S_ac;
  Double I_ERI_G2y2z_Py_Px_S_ac = I_ERI_H3y2z_S_Px_S_ac+ABY*I_ERI_G2y2z_S_Px_S_ac;
  Double I_ERI_Gy3z_Py_Px_S_ac = I_ERI_H2y3z_S_Px_S_ac+ABY*I_ERI_Gy3z_S_Px_S_ac;
  Double I_ERI_G4z_Py_Px_S_ac = I_ERI_Hy4z_S_Px_S_ac+ABY*I_ERI_G4z_S_Px_S_ac;
  Double I_ERI_G4x_Pz_Px_S_ac = I_ERI_H4xz_S_Px_S_ac+ABZ*I_ERI_G4x_S_Px_S_ac;
  Double I_ERI_G3xy_Pz_Px_S_ac = I_ERI_H3xyz_S_Px_S_ac+ABZ*I_ERI_G3xy_S_Px_S_ac;
  Double I_ERI_G3xz_Pz_Px_S_ac = I_ERI_H3x2z_S_Px_S_ac+ABZ*I_ERI_G3xz_S_Px_S_ac;
  Double I_ERI_G2x2y_Pz_Px_S_ac = I_ERI_H2x2yz_S_Px_S_ac+ABZ*I_ERI_G2x2y_S_Px_S_ac;
  Double I_ERI_G2xyz_Pz_Px_S_ac = I_ERI_H2xy2z_S_Px_S_ac+ABZ*I_ERI_G2xyz_S_Px_S_ac;
  Double I_ERI_G2x2z_Pz_Px_S_ac = I_ERI_H2x3z_S_Px_S_ac+ABZ*I_ERI_G2x2z_S_Px_S_ac;
  Double I_ERI_Gx3y_Pz_Px_S_ac = I_ERI_Hx3yz_S_Px_S_ac+ABZ*I_ERI_Gx3y_S_Px_S_ac;
  Double I_ERI_Gx2yz_Pz_Px_S_ac = I_ERI_Hx2y2z_S_Px_S_ac+ABZ*I_ERI_Gx2yz_S_Px_S_ac;
  Double I_ERI_Gxy2z_Pz_Px_S_ac = I_ERI_Hxy3z_S_Px_S_ac+ABZ*I_ERI_Gxy2z_S_Px_S_ac;
  Double I_ERI_Gx3z_Pz_Px_S_ac = I_ERI_Hx4z_S_Px_S_ac+ABZ*I_ERI_Gx3z_S_Px_S_ac;
  Double I_ERI_G4y_Pz_Px_S_ac = I_ERI_H4yz_S_Px_S_ac+ABZ*I_ERI_G4y_S_Px_S_ac;
  Double I_ERI_G3yz_Pz_Px_S_ac = I_ERI_H3y2z_S_Px_S_ac+ABZ*I_ERI_G3yz_S_Px_S_ac;
  Double I_ERI_G2y2z_Pz_Px_S_ac = I_ERI_H2y3z_S_Px_S_ac+ABZ*I_ERI_G2y2z_S_Px_S_ac;
  Double I_ERI_Gy3z_Pz_Px_S_ac = I_ERI_Hy4z_S_Px_S_ac+ABZ*I_ERI_Gy3z_S_Px_S_ac;
  Double I_ERI_G4z_Pz_Px_S_ac = I_ERI_H5z_S_Px_S_ac+ABZ*I_ERI_G4z_S_Px_S_ac;
  Double I_ERI_G4x_Px_Py_S_ac = I_ERI_H5x_S_Py_S_ac+ABX*I_ERI_G4x_S_Py_S_ac;
  Double I_ERI_G3xy_Px_Py_S_ac = I_ERI_H4xy_S_Py_S_ac+ABX*I_ERI_G3xy_S_Py_S_ac;
  Double I_ERI_G3xz_Px_Py_S_ac = I_ERI_H4xz_S_Py_S_ac+ABX*I_ERI_G3xz_S_Py_S_ac;
  Double I_ERI_G2x2y_Px_Py_S_ac = I_ERI_H3x2y_S_Py_S_ac+ABX*I_ERI_G2x2y_S_Py_S_ac;
  Double I_ERI_G2xyz_Px_Py_S_ac = I_ERI_H3xyz_S_Py_S_ac+ABX*I_ERI_G2xyz_S_Py_S_ac;
  Double I_ERI_G2x2z_Px_Py_S_ac = I_ERI_H3x2z_S_Py_S_ac+ABX*I_ERI_G2x2z_S_Py_S_ac;
  Double I_ERI_Gx3y_Px_Py_S_ac = I_ERI_H2x3y_S_Py_S_ac+ABX*I_ERI_Gx3y_S_Py_S_ac;
  Double I_ERI_Gx2yz_Px_Py_S_ac = I_ERI_H2x2yz_S_Py_S_ac+ABX*I_ERI_Gx2yz_S_Py_S_ac;
  Double I_ERI_Gxy2z_Px_Py_S_ac = I_ERI_H2xy2z_S_Py_S_ac+ABX*I_ERI_Gxy2z_S_Py_S_ac;
  Double I_ERI_Gx3z_Px_Py_S_ac = I_ERI_H2x3z_S_Py_S_ac+ABX*I_ERI_Gx3z_S_Py_S_ac;
  Double I_ERI_G4y_Px_Py_S_ac = I_ERI_Hx4y_S_Py_S_ac+ABX*I_ERI_G4y_S_Py_S_ac;
  Double I_ERI_G3yz_Px_Py_S_ac = I_ERI_Hx3yz_S_Py_S_ac+ABX*I_ERI_G3yz_S_Py_S_ac;
  Double I_ERI_G2y2z_Px_Py_S_ac = I_ERI_Hx2y2z_S_Py_S_ac+ABX*I_ERI_G2y2z_S_Py_S_ac;
  Double I_ERI_Gy3z_Px_Py_S_ac = I_ERI_Hxy3z_S_Py_S_ac+ABX*I_ERI_Gy3z_S_Py_S_ac;
  Double I_ERI_G4z_Px_Py_S_ac = I_ERI_Hx4z_S_Py_S_ac+ABX*I_ERI_G4z_S_Py_S_ac;
  Double I_ERI_G4x_Py_Py_S_ac = I_ERI_H4xy_S_Py_S_ac+ABY*I_ERI_G4x_S_Py_S_ac;
  Double I_ERI_G3xy_Py_Py_S_ac = I_ERI_H3x2y_S_Py_S_ac+ABY*I_ERI_G3xy_S_Py_S_ac;
  Double I_ERI_G3xz_Py_Py_S_ac = I_ERI_H3xyz_S_Py_S_ac+ABY*I_ERI_G3xz_S_Py_S_ac;
  Double I_ERI_G2x2y_Py_Py_S_ac = I_ERI_H2x3y_S_Py_S_ac+ABY*I_ERI_G2x2y_S_Py_S_ac;
  Double I_ERI_G2xyz_Py_Py_S_ac = I_ERI_H2x2yz_S_Py_S_ac+ABY*I_ERI_G2xyz_S_Py_S_ac;
  Double I_ERI_G2x2z_Py_Py_S_ac = I_ERI_H2xy2z_S_Py_S_ac+ABY*I_ERI_G2x2z_S_Py_S_ac;
  Double I_ERI_Gx3y_Py_Py_S_ac = I_ERI_Hx4y_S_Py_S_ac+ABY*I_ERI_Gx3y_S_Py_S_ac;
  Double I_ERI_Gx2yz_Py_Py_S_ac = I_ERI_Hx3yz_S_Py_S_ac+ABY*I_ERI_Gx2yz_S_Py_S_ac;
  Double I_ERI_Gxy2z_Py_Py_S_ac = I_ERI_Hx2y2z_S_Py_S_ac+ABY*I_ERI_Gxy2z_S_Py_S_ac;
  Double I_ERI_Gx3z_Py_Py_S_ac = I_ERI_Hxy3z_S_Py_S_ac+ABY*I_ERI_Gx3z_S_Py_S_ac;
  Double I_ERI_G4y_Py_Py_S_ac = I_ERI_H5y_S_Py_S_ac+ABY*I_ERI_G4y_S_Py_S_ac;
  Double I_ERI_G3yz_Py_Py_S_ac = I_ERI_H4yz_S_Py_S_ac+ABY*I_ERI_G3yz_S_Py_S_ac;
  Double I_ERI_G2y2z_Py_Py_S_ac = I_ERI_H3y2z_S_Py_S_ac+ABY*I_ERI_G2y2z_S_Py_S_ac;
  Double I_ERI_Gy3z_Py_Py_S_ac = I_ERI_H2y3z_S_Py_S_ac+ABY*I_ERI_Gy3z_S_Py_S_ac;
  Double I_ERI_G4z_Py_Py_S_ac = I_ERI_Hy4z_S_Py_S_ac+ABY*I_ERI_G4z_S_Py_S_ac;
  Double I_ERI_G4x_Pz_Py_S_ac = I_ERI_H4xz_S_Py_S_ac+ABZ*I_ERI_G4x_S_Py_S_ac;
  Double I_ERI_G3xy_Pz_Py_S_ac = I_ERI_H3xyz_S_Py_S_ac+ABZ*I_ERI_G3xy_S_Py_S_ac;
  Double I_ERI_G3xz_Pz_Py_S_ac = I_ERI_H3x2z_S_Py_S_ac+ABZ*I_ERI_G3xz_S_Py_S_ac;
  Double I_ERI_G2x2y_Pz_Py_S_ac = I_ERI_H2x2yz_S_Py_S_ac+ABZ*I_ERI_G2x2y_S_Py_S_ac;
  Double I_ERI_G2xyz_Pz_Py_S_ac = I_ERI_H2xy2z_S_Py_S_ac+ABZ*I_ERI_G2xyz_S_Py_S_ac;
  Double I_ERI_G2x2z_Pz_Py_S_ac = I_ERI_H2x3z_S_Py_S_ac+ABZ*I_ERI_G2x2z_S_Py_S_ac;
  Double I_ERI_Gx3y_Pz_Py_S_ac = I_ERI_Hx3yz_S_Py_S_ac+ABZ*I_ERI_Gx3y_S_Py_S_ac;
  Double I_ERI_Gx2yz_Pz_Py_S_ac = I_ERI_Hx2y2z_S_Py_S_ac+ABZ*I_ERI_Gx2yz_S_Py_S_ac;
  Double I_ERI_Gxy2z_Pz_Py_S_ac = I_ERI_Hxy3z_S_Py_S_ac+ABZ*I_ERI_Gxy2z_S_Py_S_ac;
  Double I_ERI_Gx3z_Pz_Py_S_ac = I_ERI_Hx4z_S_Py_S_ac+ABZ*I_ERI_Gx3z_S_Py_S_ac;
  Double I_ERI_G4y_Pz_Py_S_ac = I_ERI_H4yz_S_Py_S_ac+ABZ*I_ERI_G4y_S_Py_S_ac;
  Double I_ERI_G3yz_Pz_Py_S_ac = I_ERI_H3y2z_S_Py_S_ac+ABZ*I_ERI_G3yz_S_Py_S_ac;
  Double I_ERI_G2y2z_Pz_Py_S_ac = I_ERI_H2y3z_S_Py_S_ac+ABZ*I_ERI_G2y2z_S_Py_S_ac;
  Double I_ERI_Gy3z_Pz_Py_S_ac = I_ERI_Hy4z_S_Py_S_ac+ABZ*I_ERI_Gy3z_S_Py_S_ac;
  Double I_ERI_G4z_Pz_Py_S_ac = I_ERI_H5z_S_Py_S_ac+ABZ*I_ERI_G4z_S_Py_S_ac;
  Double I_ERI_G4x_Px_Pz_S_ac = I_ERI_H5x_S_Pz_S_ac+ABX*I_ERI_G4x_S_Pz_S_ac;
  Double I_ERI_G3xy_Px_Pz_S_ac = I_ERI_H4xy_S_Pz_S_ac+ABX*I_ERI_G3xy_S_Pz_S_ac;
  Double I_ERI_G3xz_Px_Pz_S_ac = I_ERI_H4xz_S_Pz_S_ac+ABX*I_ERI_G3xz_S_Pz_S_ac;
  Double I_ERI_G2x2y_Px_Pz_S_ac = I_ERI_H3x2y_S_Pz_S_ac+ABX*I_ERI_G2x2y_S_Pz_S_ac;
  Double I_ERI_G2xyz_Px_Pz_S_ac = I_ERI_H3xyz_S_Pz_S_ac+ABX*I_ERI_G2xyz_S_Pz_S_ac;
  Double I_ERI_G2x2z_Px_Pz_S_ac = I_ERI_H3x2z_S_Pz_S_ac+ABX*I_ERI_G2x2z_S_Pz_S_ac;
  Double I_ERI_Gx3y_Px_Pz_S_ac = I_ERI_H2x3y_S_Pz_S_ac+ABX*I_ERI_Gx3y_S_Pz_S_ac;
  Double I_ERI_Gx2yz_Px_Pz_S_ac = I_ERI_H2x2yz_S_Pz_S_ac+ABX*I_ERI_Gx2yz_S_Pz_S_ac;
  Double I_ERI_Gxy2z_Px_Pz_S_ac = I_ERI_H2xy2z_S_Pz_S_ac+ABX*I_ERI_Gxy2z_S_Pz_S_ac;
  Double I_ERI_Gx3z_Px_Pz_S_ac = I_ERI_H2x3z_S_Pz_S_ac+ABX*I_ERI_Gx3z_S_Pz_S_ac;
  Double I_ERI_G4y_Px_Pz_S_ac = I_ERI_Hx4y_S_Pz_S_ac+ABX*I_ERI_G4y_S_Pz_S_ac;
  Double I_ERI_G3yz_Px_Pz_S_ac = I_ERI_Hx3yz_S_Pz_S_ac+ABX*I_ERI_G3yz_S_Pz_S_ac;
  Double I_ERI_G2y2z_Px_Pz_S_ac = I_ERI_Hx2y2z_S_Pz_S_ac+ABX*I_ERI_G2y2z_S_Pz_S_ac;
  Double I_ERI_Gy3z_Px_Pz_S_ac = I_ERI_Hxy3z_S_Pz_S_ac+ABX*I_ERI_Gy3z_S_Pz_S_ac;
  Double I_ERI_G4z_Px_Pz_S_ac = I_ERI_Hx4z_S_Pz_S_ac+ABX*I_ERI_G4z_S_Pz_S_ac;
  Double I_ERI_G4x_Py_Pz_S_ac = I_ERI_H4xy_S_Pz_S_ac+ABY*I_ERI_G4x_S_Pz_S_ac;
  Double I_ERI_G3xy_Py_Pz_S_ac = I_ERI_H3x2y_S_Pz_S_ac+ABY*I_ERI_G3xy_S_Pz_S_ac;
  Double I_ERI_G3xz_Py_Pz_S_ac = I_ERI_H3xyz_S_Pz_S_ac+ABY*I_ERI_G3xz_S_Pz_S_ac;
  Double I_ERI_G2x2y_Py_Pz_S_ac = I_ERI_H2x3y_S_Pz_S_ac+ABY*I_ERI_G2x2y_S_Pz_S_ac;
  Double I_ERI_G2xyz_Py_Pz_S_ac = I_ERI_H2x2yz_S_Pz_S_ac+ABY*I_ERI_G2xyz_S_Pz_S_ac;
  Double I_ERI_G2x2z_Py_Pz_S_ac = I_ERI_H2xy2z_S_Pz_S_ac+ABY*I_ERI_G2x2z_S_Pz_S_ac;
  Double I_ERI_Gx3y_Py_Pz_S_ac = I_ERI_Hx4y_S_Pz_S_ac+ABY*I_ERI_Gx3y_S_Pz_S_ac;
  Double I_ERI_Gx2yz_Py_Pz_S_ac = I_ERI_Hx3yz_S_Pz_S_ac+ABY*I_ERI_Gx2yz_S_Pz_S_ac;
  Double I_ERI_Gxy2z_Py_Pz_S_ac = I_ERI_Hx2y2z_S_Pz_S_ac+ABY*I_ERI_Gxy2z_S_Pz_S_ac;
  Double I_ERI_Gx3z_Py_Pz_S_ac = I_ERI_Hxy3z_S_Pz_S_ac+ABY*I_ERI_Gx3z_S_Pz_S_ac;
  Double I_ERI_G4y_Py_Pz_S_ac = I_ERI_H5y_S_Pz_S_ac+ABY*I_ERI_G4y_S_Pz_S_ac;
  Double I_ERI_G3yz_Py_Pz_S_ac = I_ERI_H4yz_S_Pz_S_ac+ABY*I_ERI_G3yz_S_Pz_S_ac;
  Double I_ERI_G2y2z_Py_Pz_S_ac = I_ERI_H3y2z_S_Pz_S_ac+ABY*I_ERI_G2y2z_S_Pz_S_ac;
  Double I_ERI_Gy3z_Py_Pz_S_ac = I_ERI_H2y3z_S_Pz_S_ac+ABY*I_ERI_Gy3z_S_Pz_S_ac;
  Double I_ERI_G4z_Py_Pz_S_ac = I_ERI_Hy4z_S_Pz_S_ac+ABY*I_ERI_G4z_S_Pz_S_ac;
  Double I_ERI_G4x_Pz_Pz_S_ac = I_ERI_H4xz_S_Pz_S_ac+ABZ*I_ERI_G4x_S_Pz_S_ac;
  Double I_ERI_G3xy_Pz_Pz_S_ac = I_ERI_H3xyz_S_Pz_S_ac+ABZ*I_ERI_G3xy_S_Pz_S_ac;
  Double I_ERI_G3xz_Pz_Pz_S_ac = I_ERI_H3x2z_S_Pz_S_ac+ABZ*I_ERI_G3xz_S_Pz_S_ac;
  Double I_ERI_G2x2y_Pz_Pz_S_ac = I_ERI_H2x2yz_S_Pz_S_ac+ABZ*I_ERI_G2x2y_S_Pz_S_ac;
  Double I_ERI_G2xyz_Pz_Pz_S_ac = I_ERI_H2xy2z_S_Pz_S_ac+ABZ*I_ERI_G2xyz_S_Pz_S_ac;
  Double I_ERI_G2x2z_Pz_Pz_S_ac = I_ERI_H2x3z_S_Pz_S_ac+ABZ*I_ERI_G2x2z_S_Pz_S_ac;
  Double I_ERI_Gx3y_Pz_Pz_S_ac = I_ERI_Hx3yz_S_Pz_S_ac+ABZ*I_ERI_Gx3y_S_Pz_S_ac;
  Double I_ERI_Gx2yz_Pz_Pz_S_ac = I_ERI_Hx2y2z_S_Pz_S_ac+ABZ*I_ERI_Gx2yz_S_Pz_S_ac;
  Double I_ERI_Gxy2z_Pz_Pz_S_ac = I_ERI_Hxy3z_S_Pz_S_ac+ABZ*I_ERI_Gxy2z_S_Pz_S_ac;
  Double I_ERI_Gx3z_Pz_Pz_S_ac = I_ERI_Hx4z_S_Pz_S_ac+ABZ*I_ERI_Gx3z_S_Pz_S_ac;
  Double I_ERI_G4y_Pz_Pz_S_ac = I_ERI_H4yz_S_Pz_S_ac+ABZ*I_ERI_G4y_S_Pz_S_ac;
  Double I_ERI_G3yz_Pz_Pz_S_ac = I_ERI_H3y2z_S_Pz_S_ac+ABZ*I_ERI_G3yz_S_Pz_S_ac;
  Double I_ERI_G2y2z_Pz_Pz_S_ac = I_ERI_H2y3z_S_Pz_S_ac+ABZ*I_ERI_G2y2z_S_Pz_S_ac;
  Double I_ERI_Gy3z_Pz_Pz_S_ac = I_ERI_Hy4z_S_Pz_S_ac+ABZ*I_ERI_Gy3z_S_Pz_S_ac;
  Double I_ERI_G4z_Pz_Pz_S_ac = I_ERI_H5z_S_Pz_S_ac+ABZ*I_ERI_G4z_S_Pz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_H_P_P_S_ac
   * expanding position: BRA2
   * code section is: HRR
   * totally 21 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_H_S_P_S_ac
   ************************************************************/
  Double I_ERI_H5x_Px_Px_S_ac = I_ERI_I6x_S_Px_S_ac+ABX*I_ERI_H5x_S_Px_S_ac;
  Double I_ERI_H4xy_Px_Px_S_ac = I_ERI_I5xy_S_Px_S_ac+ABX*I_ERI_H4xy_S_Px_S_ac;
  Double I_ERI_H4xz_Px_Px_S_ac = I_ERI_I5xz_S_Px_S_ac+ABX*I_ERI_H4xz_S_Px_S_ac;
  Double I_ERI_H3x2y_Px_Px_S_ac = I_ERI_I4x2y_S_Px_S_ac+ABX*I_ERI_H3x2y_S_Px_S_ac;
  Double I_ERI_H3xyz_Px_Px_S_ac = I_ERI_I4xyz_S_Px_S_ac+ABX*I_ERI_H3xyz_S_Px_S_ac;
  Double I_ERI_H3x2z_Px_Px_S_ac = I_ERI_I4x2z_S_Px_S_ac+ABX*I_ERI_H3x2z_S_Px_S_ac;
  Double I_ERI_H2x3y_Px_Px_S_ac = I_ERI_I3x3y_S_Px_S_ac+ABX*I_ERI_H2x3y_S_Px_S_ac;
  Double I_ERI_H2x2yz_Px_Px_S_ac = I_ERI_I3x2yz_S_Px_S_ac+ABX*I_ERI_H2x2yz_S_Px_S_ac;
  Double I_ERI_H2xy2z_Px_Px_S_ac = I_ERI_I3xy2z_S_Px_S_ac+ABX*I_ERI_H2xy2z_S_Px_S_ac;
  Double I_ERI_H2x3z_Px_Px_S_ac = I_ERI_I3x3z_S_Px_S_ac+ABX*I_ERI_H2x3z_S_Px_S_ac;
  Double I_ERI_Hx4y_Px_Px_S_ac = I_ERI_I2x4y_S_Px_S_ac+ABX*I_ERI_Hx4y_S_Px_S_ac;
  Double I_ERI_Hx3yz_Px_Px_S_ac = I_ERI_I2x3yz_S_Px_S_ac+ABX*I_ERI_Hx3yz_S_Px_S_ac;
  Double I_ERI_Hx2y2z_Px_Px_S_ac = I_ERI_I2x2y2z_S_Px_S_ac+ABX*I_ERI_Hx2y2z_S_Px_S_ac;
  Double I_ERI_Hxy3z_Px_Px_S_ac = I_ERI_I2xy3z_S_Px_S_ac+ABX*I_ERI_Hxy3z_S_Px_S_ac;
  Double I_ERI_Hx4z_Px_Px_S_ac = I_ERI_I2x4z_S_Px_S_ac+ABX*I_ERI_Hx4z_S_Px_S_ac;
  Double I_ERI_H5y_Px_Px_S_ac = I_ERI_Ix5y_S_Px_S_ac+ABX*I_ERI_H5y_S_Px_S_ac;
  Double I_ERI_H4yz_Px_Px_S_ac = I_ERI_Ix4yz_S_Px_S_ac+ABX*I_ERI_H4yz_S_Px_S_ac;
  Double I_ERI_H3y2z_Px_Px_S_ac = I_ERI_Ix3y2z_S_Px_S_ac+ABX*I_ERI_H3y2z_S_Px_S_ac;
  Double I_ERI_H2y3z_Px_Px_S_ac = I_ERI_Ix2y3z_S_Px_S_ac+ABX*I_ERI_H2y3z_S_Px_S_ac;
  Double I_ERI_Hy4z_Px_Px_S_ac = I_ERI_Ixy4z_S_Px_S_ac+ABX*I_ERI_Hy4z_S_Px_S_ac;
  Double I_ERI_H5z_Px_Px_S_ac = I_ERI_Ix5z_S_Px_S_ac+ABX*I_ERI_H5z_S_Px_S_ac;
  Double I_ERI_H4xy_Py_Px_S_ac = I_ERI_I4x2y_S_Px_S_ac+ABY*I_ERI_H4xy_S_Px_S_ac;
  Double I_ERI_H4xz_Py_Px_S_ac = I_ERI_I4xyz_S_Px_S_ac+ABY*I_ERI_H4xz_S_Px_S_ac;
  Double I_ERI_H3x2y_Py_Px_S_ac = I_ERI_I3x3y_S_Px_S_ac+ABY*I_ERI_H3x2y_S_Px_S_ac;
  Double I_ERI_H3xyz_Py_Px_S_ac = I_ERI_I3x2yz_S_Px_S_ac+ABY*I_ERI_H3xyz_S_Px_S_ac;
  Double I_ERI_H3x2z_Py_Px_S_ac = I_ERI_I3xy2z_S_Px_S_ac+ABY*I_ERI_H3x2z_S_Px_S_ac;
  Double I_ERI_H2x3y_Py_Px_S_ac = I_ERI_I2x4y_S_Px_S_ac+ABY*I_ERI_H2x3y_S_Px_S_ac;
  Double I_ERI_H2x2yz_Py_Px_S_ac = I_ERI_I2x3yz_S_Px_S_ac+ABY*I_ERI_H2x2yz_S_Px_S_ac;
  Double I_ERI_H2xy2z_Py_Px_S_ac = I_ERI_I2x2y2z_S_Px_S_ac+ABY*I_ERI_H2xy2z_S_Px_S_ac;
  Double I_ERI_H2x3z_Py_Px_S_ac = I_ERI_I2xy3z_S_Px_S_ac+ABY*I_ERI_H2x3z_S_Px_S_ac;
  Double I_ERI_Hx4y_Py_Px_S_ac = I_ERI_Ix5y_S_Px_S_ac+ABY*I_ERI_Hx4y_S_Px_S_ac;
  Double I_ERI_Hx3yz_Py_Px_S_ac = I_ERI_Ix4yz_S_Px_S_ac+ABY*I_ERI_Hx3yz_S_Px_S_ac;
  Double I_ERI_Hx2y2z_Py_Px_S_ac = I_ERI_Ix3y2z_S_Px_S_ac+ABY*I_ERI_Hx2y2z_S_Px_S_ac;
  Double I_ERI_Hxy3z_Py_Px_S_ac = I_ERI_Ix2y3z_S_Px_S_ac+ABY*I_ERI_Hxy3z_S_Px_S_ac;
  Double I_ERI_Hx4z_Py_Px_S_ac = I_ERI_Ixy4z_S_Px_S_ac+ABY*I_ERI_Hx4z_S_Px_S_ac;
  Double I_ERI_H5y_Py_Px_S_ac = I_ERI_I6y_S_Px_S_ac+ABY*I_ERI_H5y_S_Px_S_ac;
  Double I_ERI_H4yz_Py_Px_S_ac = I_ERI_I5yz_S_Px_S_ac+ABY*I_ERI_H4yz_S_Px_S_ac;
  Double I_ERI_H3y2z_Py_Px_S_ac = I_ERI_I4y2z_S_Px_S_ac+ABY*I_ERI_H3y2z_S_Px_S_ac;
  Double I_ERI_H2y3z_Py_Px_S_ac = I_ERI_I3y3z_S_Px_S_ac+ABY*I_ERI_H2y3z_S_Px_S_ac;
  Double I_ERI_Hy4z_Py_Px_S_ac = I_ERI_I2y4z_S_Px_S_ac+ABY*I_ERI_Hy4z_S_Px_S_ac;
  Double I_ERI_H5z_Py_Px_S_ac = I_ERI_Iy5z_S_Px_S_ac+ABY*I_ERI_H5z_S_Px_S_ac;
  Double I_ERI_H4xz_Pz_Px_S_ac = I_ERI_I4x2z_S_Px_S_ac+ABZ*I_ERI_H4xz_S_Px_S_ac;
  Double I_ERI_H3xyz_Pz_Px_S_ac = I_ERI_I3xy2z_S_Px_S_ac+ABZ*I_ERI_H3xyz_S_Px_S_ac;
  Double I_ERI_H3x2z_Pz_Px_S_ac = I_ERI_I3x3z_S_Px_S_ac+ABZ*I_ERI_H3x2z_S_Px_S_ac;
  Double I_ERI_H2x2yz_Pz_Px_S_ac = I_ERI_I2x2y2z_S_Px_S_ac+ABZ*I_ERI_H2x2yz_S_Px_S_ac;
  Double I_ERI_H2xy2z_Pz_Px_S_ac = I_ERI_I2xy3z_S_Px_S_ac+ABZ*I_ERI_H2xy2z_S_Px_S_ac;
  Double I_ERI_H2x3z_Pz_Px_S_ac = I_ERI_I2x4z_S_Px_S_ac+ABZ*I_ERI_H2x3z_S_Px_S_ac;
  Double I_ERI_Hx3yz_Pz_Px_S_ac = I_ERI_Ix3y2z_S_Px_S_ac+ABZ*I_ERI_Hx3yz_S_Px_S_ac;
  Double I_ERI_Hx2y2z_Pz_Px_S_ac = I_ERI_Ix2y3z_S_Px_S_ac+ABZ*I_ERI_Hx2y2z_S_Px_S_ac;
  Double I_ERI_Hxy3z_Pz_Px_S_ac = I_ERI_Ixy4z_S_Px_S_ac+ABZ*I_ERI_Hxy3z_S_Px_S_ac;
  Double I_ERI_Hx4z_Pz_Px_S_ac = I_ERI_Ix5z_S_Px_S_ac+ABZ*I_ERI_Hx4z_S_Px_S_ac;
  Double I_ERI_H4yz_Pz_Px_S_ac = I_ERI_I4y2z_S_Px_S_ac+ABZ*I_ERI_H4yz_S_Px_S_ac;
  Double I_ERI_H3y2z_Pz_Px_S_ac = I_ERI_I3y3z_S_Px_S_ac+ABZ*I_ERI_H3y2z_S_Px_S_ac;
  Double I_ERI_H2y3z_Pz_Px_S_ac = I_ERI_I2y4z_S_Px_S_ac+ABZ*I_ERI_H2y3z_S_Px_S_ac;
  Double I_ERI_Hy4z_Pz_Px_S_ac = I_ERI_Iy5z_S_Px_S_ac+ABZ*I_ERI_Hy4z_S_Px_S_ac;
  Double I_ERI_H5z_Pz_Px_S_ac = I_ERI_I6z_S_Px_S_ac+ABZ*I_ERI_H5z_S_Px_S_ac;
  Double I_ERI_H5x_Px_Py_S_ac = I_ERI_I6x_S_Py_S_ac+ABX*I_ERI_H5x_S_Py_S_ac;
  Double I_ERI_H4xy_Px_Py_S_ac = I_ERI_I5xy_S_Py_S_ac+ABX*I_ERI_H4xy_S_Py_S_ac;
  Double I_ERI_H4xz_Px_Py_S_ac = I_ERI_I5xz_S_Py_S_ac+ABX*I_ERI_H4xz_S_Py_S_ac;
  Double I_ERI_H3x2y_Px_Py_S_ac = I_ERI_I4x2y_S_Py_S_ac+ABX*I_ERI_H3x2y_S_Py_S_ac;
  Double I_ERI_H3xyz_Px_Py_S_ac = I_ERI_I4xyz_S_Py_S_ac+ABX*I_ERI_H3xyz_S_Py_S_ac;
  Double I_ERI_H3x2z_Px_Py_S_ac = I_ERI_I4x2z_S_Py_S_ac+ABX*I_ERI_H3x2z_S_Py_S_ac;
  Double I_ERI_H2x3y_Px_Py_S_ac = I_ERI_I3x3y_S_Py_S_ac+ABX*I_ERI_H2x3y_S_Py_S_ac;
  Double I_ERI_H2x2yz_Px_Py_S_ac = I_ERI_I3x2yz_S_Py_S_ac+ABX*I_ERI_H2x2yz_S_Py_S_ac;
  Double I_ERI_H2xy2z_Px_Py_S_ac = I_ERI_I3xy2z_S_Py_S_ac+ABX*I_ERI_H2xy2z_S_Py_S_ac;
  Double I_ERI_H2x3z_Px_Py_S_ac = I_ERI_I3x3z_S_Py_S_ac+ABX*I_ERI_H2x3z_S_Py_S_ac;
  Double I_ERI_Hx4y_Px_Py_S_ac = I_ERI_I2x4y_S_Py_S_ac+ABX*I_ERI_Hx4y_S_Py_S_ac;
  Double I_ERI_Hx3yz_Px_Py_S_ac = I_ERI_I2x3yz_S_Py_S_ac+ABX*I_ERI_Hx3yz_S_Py_S_ac;
  Double I_ERI_Hx2y2z_Px_Py_S_ac = I_ERI_I2x2y2z_S_Py_S_ac+ABX*I_ERI_Hx2y2z_S_Py_S_ac;
  Double I_ERI_Hxy3z_Px_Py_S_ac = I_ERI_I2xy3z_S_Py_S_ac+ABX*I_ERI_Hxy3z_S_Py_S_ac;
  Double I_ERI_Hx4z_Px_Py_S_ac = I_ERI_I2x4z_S_Py_S_ac+ABX*I_ERI_Hx4z_S_Py_S_ac;
  Double I_ERI_H5y_Px_Py_S_ac = I_ERI_Ix5y_S_Py_S_ac+ABX*I_ERI_H5y_S_Py_S_ac;
  Double I_ERI_H4yz_Px_Py_S_ac = I_ERI_Ix4yz_S_Py_S_ac+ABX*I_ERI_H4yz_S_Py_S_ac;
  Double I_ERI_H3y2z_Px_Py_S_ac = I_ERI_Ix3y2z_S_Py_S_ac+ABX*I_ERI_H3y2z_S_Py_S_ac;
  Double I_ERI_H2y3z_Px_Py_S_ac = I_ERI_Ix2y3z_S_Py_S_ac+ABX*I_ERI_H2y3z_S_Py_S_ac;
  Double I_ERI_Hy4z_Px_Py_S_ac = I_ERI_Ixy4z_S_Py_S_ac+ABX*I_ERI_Hy4z_S_Py_S_ac;
  Double I_ERI_H5z_Px_Py_S_ac = I_ERI_Ix5z_S_Py_S_ac+ABX*I_ERI_H5z_S_Py_S_ac;
  Double I_ERI_H4xy_Py_Py_S_ac = I_ERI_I4x2y_S_Py_S_ac+ABY*I_ERI_H4xy_S_Py_S_ac;
  Double I_ERI_H4xz_Py_Py_S_ac = I_ERI_I4xyz_S_Py_S_ac+ABY*I_ERI_H4xz_S_Py_S_ac;
  Double I_ERI_H3x2y_Py_Py_S_ac = I_ERI_I3x3y_S_Py_S_ac+ABY*I_ERI_H3x2y_S_Py_S_ac;
  Double I_ERI_H3xyz_Py_Py_S_ac = I_ERI_I3x2yz_S_Py_S_ac+ABY*I_ERI_H3xyz_S_Py_S_ac;
  Double I_ERI_H3x2z_Py_Py_S_ac = I_ERI_I3xy2z_S_Py_S_ac+ABY*I_ERI_H3x2z_S_Py_S_ac;
  Double I_ERI_H2x3y_Py_Py_S_ac = I_ERI_I2x4y_S_Py_S_ac+ABY*I_ERI_H2x3y_S_Py_S_ac;
  Double I_ERI_H2x2yz_Py_Py_S_ac = I_ERI_I2x3yz_S_Py_S_ac+ABY*I_ERI_H2x2yz_S_Py_S_ac;
  Double I_ERI_H2xy2z_Py_Py_S_ac = I_ERI_I2x2y2z_S_Py_S_ac+ABY*I_ERI_H2xy2z_S_Py_S_ac;
  Double I_ERI_H2x3z_Py_Py_S_ac = I_ERI_I2xy3z_S_Py_S_ac+ABY*I_ERI_H2x3z_S_Py_S_ac;
  Double I_ERI_Hx4y_Py_Py_S_ac = I_ERI_Ix5y_S_Py_S_ac+ABY*I_ERI_Hx4y_S_Py_S_ac;
  Double I_ERI_Hx3yz_Py_Py_S_ac = I_ERI_Ix4yz_S_Py_S_ac+ABY*I_ERI_Hx3yz_S_Py_S_ac;
  Double I_ERI_Hx2y2z_Py_Py_S_ac = I_ERI_Ix3y2z_S_Py_S_ac+ABY*I_ERI_Hx2y2z_S_Py_S_ac;
  Double I_ERI_Hxy3z_Py_Py_S_ac = I_ERI_Ix2y3z_S_Py_S_ac+ABY*I_ERI_Hxy3z_S_Py_S_ac;
  Double I_ERI_Hx4z_Py_Py_S_ac = I_ERI_Ixy4z_S_Py_S_ac+ABY*I_ERI_Hx4z_S_Py_S_ac;
  Double I_ERI_H5y_Py_Py_S_ac = I_ERI_I6y_S_Py_S_ac+ABY*I_ERI_H5y_S_Py_S_ac;
  Double I_ERI_H4yz_Py_Py_S_ac = I_ERI_I5yz_S_Py_S_ac+ABY*I_ERI_H4yz_S_Py_S_ac;
  Double I_ERI_H3y2z_Py_Py_S_ac = I_ERI_I4y2z_S_Py_S_ac+ABY*I_ERI_H3y2z_S_Py_S_ac;
  Double I_ERI_H2y3z_Py_Py_S_ac = I_ERI_I3y3z_S_Py_S_ac+ABY*I_ERI_H2y3z_S_Py_S_ac;
  Double I_ERI_Hy4z_Py_Py_S_ac = I_ERI_I2y4z_S_Py_S_ac+ABY*I_ERI_Hy4z_S_Py_S_ac;
  Double I_ERI_H5z_Py_Py_S_ac = I_ERI_Iy5z_S_Py_S_ac+ABY*I_ERI_H5z_S_Py_S_ac;
  Double I_ERI_H4xz_Pz_Py_S_ac = I_ERI_I4x2z_S_Py_S_ac+ABZ*I_ERI_H4xz_S_Py_S_ac;
  Double I_ERI_H3xyz_Pz_Py_S_ac = I_ERI_I3xy2z_S_Py_S_ac+ABZ*I_ERI_H3xyz_S_Py_S_ac;
  Double I_ERI_H3x2z_Pz_Py_S_ac = I_ERI_I3x3z_S_Py_S_ac+ABZ*I_ERI_H3x2z_S_Py_S_ac;
  Double I_ERI_H2x2yz_Pz_Py_S_ac = I_ERI_I2x2y2z_S_Py_S_ac+ABZ*I_ERI_H2x2yz_S_Py_S_ac;
  Double I_ERI_H2xy2z_Pz_Py_S_ac = I_ERI_I2xy3z_S_Py_S_ac+ABZ*I_ERI_H2xy2z_S_Py_S_ac;
  Double I_ERI_H2x3z_Pz_Py_S_ac = I_ERI_I2x4z_S_Py_S_ac+ABZ*I_ERI_H2x3z_S_Py_S_ac;
  Double I_ERI_Hx3yz_Pz_Py_S_ac = I_ERI_Ix3y2z_S_Py_S_ac+ABZ*I_ERI_Hx3yz_S_Py_S_ac;
  Double I_ERI_Hx2y2z_Pz_Py_S_ac = I_ERI_Ix2y3z_S_Py_S_ac+ABZ*I_ERI_Hx2y2z_S_Py_S_ac;
  Double I_ERI_Hxy3z_Pz_Py_S_ac = I_ERI_Ixy4z_S_Py_S_ac+ABZ*I_ERI_Hxy3z_S_Py_S_ac;
  Double I_ERI_Hx4z_Pz_Py_S_ac = I_ERI_Ix5z_S_Py_S_ac+ABZ*I_ERI_Hx4z_S_Py_S_ac;
  Double I_ERI_H4yz_Pz_Py_S_ac = I_ERI_I4y2z_S_Py_S_ac+ABZ*I_ERI_H4yz_S_Py_S_ac;
  Double I_ERI_H3y2z_Pz_Py_S_ac = I_ERI_I3y3z_S_Py_S_ac+ABZ*I_ERI_H3y2z_S_Py_S_ac;
  Double I_ERI_H2y3z_Pz_Py_S_ac = I_ERI_I2y4z_S_Py_S_ac+ABZ*I_ERI_H2y3z_S_Py_S_ac;
  Double I_ERI_Hy4z_Pz_Py_S_ac = I_ERI_Iy5z_S_Py_S_ac+ABZ*I_ERI_Hy4z_S_Py_S_ac;
  Double I_ERI_H5z_Pz_Py_S_ac = I_ERI_I6z_S_Py_S_ac+ABZ*I_ERI_H5z_S_Py_S_ac;
  Double I_ERI_H5x_Px_Pz_S_ac = I_ERI_I6x_S_Pz_S_ac+ABX*I_ERI_H5x_S_Pz_S_ac;
  Double I_ERI_H4xy_Px_Pz_S_ac = I_ERI_I5xy_S_Pz_S_ac+ABX*I_ERI_H4xy_S_Pz_S_ac;
  Double I_ERI_H4xz_Px_Pz_S_ac = I_ERI_I5xz_S_Pz_S_ac+ABX*I_ERI_H4xz_S_Pz_S_ac;
  Double I_ERI_H3x2y_Px_Pz_S_ac = I_ERI_I4x2y_S_Pz_S_ac+ABX*I_ERI_H3x2y_S_Pz_S_ac;
  Double I_ERI_H3xyz_Px_Pz_S_ac = I_ERI_I4xyz_S_Pz_S_ac+ABX*I_ERI_H3xyz_S_Pz_S_ac;
  Double I_ERI_H3x2z_Px_Pz_S_ac = I_ERI_I4x2z_S_Pz_S_ac+ABX*I_ERI_H3x2z_S_Pz_S_ac;
  Double I_ERI_H2x3y_Px_Pz_S_ac = I_ERI_I3x3y_S_Pz_S_ac+ABX*I_ERI_H2x3y_S_Pz_S_ac;
  Double I_ERI_H2x2yz_Px_Pz_S_ac = I_ERI_I3x2yz_S_Pz_S_ac+ABX*I_ERI_H2x2yz_S_Pz_S_ac;
  Double I_ERI_H2xy2z_Px_Pz_S_ac = I_ERI_I3xy2z_S_Pz_S_ac+ABX*I_ERI_H2xy2z_S_Pz_S_ac;
  Double I_ERI_H2x3z_Px_Pz_S_ac = I_ERI_I3x3z_S_Pz_S_ac+ABX*I_ERI_H2x3z_S_Pz_S_ac;
  Double I_ERI_Hx4y_Px_Pz_S_ac = I_ERI_I2x4y_S_Pz_S_ac+ABX*I_ERI_Hx4y_S_Pz_S_ac;
  Double I_ERI_Hx3yz_Px_Pz_S_ac = I_ERI_I2x3yz_S_Pz_S_ac+ABX*I_ERI_Hx3yz_S_Pz_S_ac;
  Double I_ERI_Hx2y2z_Px_Pz_S_ac = I_ERI_I2x2y2z_S_Pz_S_ac+ABX*I_ERI_Hx2y2z_S_Pz_S_ac;
  Double I_ERI_Hxy3z_Px_Pz_S_ac = I_ERI_I2xy3z_S_Pz_S_ac+ABX*I_ERI_Hxy3z_S_Pz_S_ac;
  Double I_ERI_Hx4z_Px_Pz_S_ac = I_ERI_I2x4z_S_Pz_S_ac+ABX*I_ERI_Hx4z_S_Pz_S_ac;
  Double I_ERI_H5y_Px_Pz_S_ac = I_ERI_Ix5y_S_Pz_S_ac+ABX*I_ERI_H5y_S_Pz_S_ac;
  Double I_ERI_H4yz_Px_Pz_S_ac = I_ERI_Ix4yz_S_Pz_S_ac+ABX*I_ERI_H4yz_S_Pz_S_ac;
  Double I_ERI_H3y2z_Px_Pz_S_ac = I_ERI_Ix3y2z_S_Pz_S_ac+ABX*I_ERI_H3y2z_S_Pz_S_ac;
  Double I_ERI_H2y3z_Px_Pz_S_ac = I_ERI_Ix2y3z_S_Pz_S_ac+ABX*I_ERI_H2y3z_S_Pz_S_ac;
  Double I_ERI_Hy4z_Px_Pz_S_ac = I_ERI_Ixy4z_S_Pz_S_ac+ABX*I_ERI_Hy4z_S_Pz_S_ac;
  Double I_ERI_H5z_Px_Pz_S_ac = I_ERI_Ix5z_S_Pz_S_ac+ABX*I_ERI_H5z_S_Pz_S_ac;
  Double I_ERI_H4xy_Py_Pz_S_ac = I_ERI_I4x2y_S_Pz_S_ac+ABY*I_ERI_H4xy_S_Pz_S_ac;
  Double I_ERI_H4xz_Py_Pz_S_ac = I_ERI_I4xyz_S_Pz_S_ac+ABY*I_ERI_H4xz_S_Pz_S_ac;
  Double I_ERI_H3x2y_Py_Pz_S_ac = I_ERI_I3x3y_S_Pz_S_ac+ABY*I_ERI_H3x2y_S_Pz_S_ac;
  Double I_ERI_H3xyz_Py_Pz_S_ac = I_ERI_I3x2yz_S_Pz_S_ac+ABY*I_ERI_H3xyz_S_Pz_S_ac;
  Double I_ERI_H3x2z_Py_Pz_S_ac = I_ERI_I3xy2z_S_Pz_S_ac+ABY*I_ERI_H3x2z_S_Pz_S_ac;
  Double I_ERI_H2x3y_Py_Pz_S_ac = I_ERI_I2x4y_S_Pz_S_ac+ABY*I_ERI_H2x3y_S_Pz_S_ac;
  Double I_ERI_H2x2yz_Py_Pz_S_ac = I_ERI_I2x3yz_S_Pz_S_ac+ABY*I_ERI_H2x2yz_S_Pz_S_ac;
  Double I_ERI_H2xy2z_Py_Pz_S_ac = I_ERI_I2x2y2z_S_Pz_S_ac+ABY*I_ERI_H2xy2z_S_Pz_S_ac;
  Double I_ERI_H2x3z_Py_Pz_S_ac = I_ERI_I2xy3z_S_Pz_S_ac+ABY*I_ERI_H2x3z_S_Pz_S_ac;
  Double I_ERI_Hx4y_Py_Pz_S_ac = I_ERI_Ix5y_S_Pz_S_ac+ABY*I_ERI_Hx4y_S_Pz_S_ac;
  Double I_ERI_Hx3yz_Py_Pz_S_ac = I_ERI_Ix4yz_S_Pz_S_ac+ABY*I_ERI_Hx3yz_S_Pz_S_ac;
  Double I_ERI_Hx2y2z_Py_Pz_S_ac = I_ERI_Ix3y2z_S_Pz_S_ac+ABY*I_ERI_Hx2y2z_S_Pz_S_ac;
  Double I_ERI_Hxy3z_Py_Pz_S_ac = I_ERI_Ix2y3z_S_Pz_S_ac+ABY*I_ERI_Hxy3z_S_Pz_S_ac;
  Double I_ERI_Hx4z_Py_Pz_S_ac = I_ERI_Ixy4z_S_Pz_S_ac+ABY*I_ERI_Hx4z_S_Pz_S_ac;
  Double I_ERI_H5y_Py_Pz_S_ac = I_ERI_I6y_S_Pz_S_ac+ABY*I_ERI_H5y_S_Pz_S_ac;
  Double I_ERI_H4yz_Py_Pz_S_ac = I_ERI_I5yz_S_Pz_S_ac+ABY*I_ERI_H4yz_S_Pz_S_ac;
  Double I_ERI_H3y2z_Py_Pz_S_ac = I_ERI_I4y2z_S_Pz_S_ac+ABY*I_ERI_H3y2z_S_Pz_S_ac;
  Double I_ERI_H2y3z_Py_Pz_S_ac = I_ERI_I3y3z_S_Pz_S_ac+ABY*I_ERI_H2y3z_S_Pz_S_ac;
  Double I_ERI_Hy4z_Py_Pz_S_ac = I_ERI_I2y4z_S_Pz_S_ac+ABY*I_ERI_Hy4z_S_Pz_S_ac;
  Double I_ERI_H5z_Py_Pz_S_ac = I_ERI_Iy5z_S_Pz_S_ac+ABY*I_ERI_H5z_S_Pz_S_ac;
  Double I_ERI_H4xz_Pz_Pz_S_ac = I_ERI_I4x2z_S_Pz_S_ac+ABZ*I_ERI_H4xz_S_Pz_S_ac;
  Double I_ERI_H3xyz_Pz_Pz_S_ac = I_ERI_I3xy2z_S_Pz_S_ac+ABZ*I_ERI_H3xyz_S_Pz_S_ac;
  Double I_ERI_H3x2z_Pz_Pz_S_ac = I_ERI_I3x3z_S_Pz_S_ac+ABZ*I_ERI_H3x2z_S_Pz_S_ac;
  Double I_ERI_H2x2yz_Pz_Pz_S_ac = I_ERI_I2x2y2z_S_Pz_S_ac+ABZ*I_ERI_H2x2yz_S_Pz_S_ac;
  Double I_ERI_H2xy2z_Pz_Pz_S_ac = I_ERI_I2xy3z_S_Pz_S_ac+ABZ*I_ERI_H2xy2z_S_Pz_S_ac;
  Double I_ERI_H2x3z_Pz_Pz_S_ac = I_ERI_I2x4z_S_Pz_S_ac+ABZ*I_ERI_H2x3z_S_Pz_S_ac;
  Double I_ERI_Hx3yz_Pz_Pz_S_ac = I_ERI_Ix3y2z_S_Pz_S_ac+ABZ*I_ERI_Hx3yz_S_Pz_S_ac;
  Double I_ERI_Hx2y2z_Pz_Pz_S_ac = I_ERI_Ix2y3z_S_Pz_S_ac+ABZ*I_ERI_Hx2y2z_S_Pz_S_ac;
  Double I_ERI_Hxy3z_Pz_Pz_S_ac = I_ERI_Ixy4z_S_Pz_S_ac+ABZ*I_ERI_Hxy3z_S_Pz_S_ac;
  Double I_ERI_Hx4z_Pz_Pz_S_ac = I_ERI_Ix5z_S_Pz_S_ac+ABZ*I_ERI_Hx4z_S_Pz_S_ac;
  Double I_ERI_H4yz_Pz_Pz_S_ac = I_ERI_I4y2z_S_Pz_S_ac+ABZ*I_ERI_H4yz_S_Pz_S_ac;
  Double I_ERI_H3y2z_Pz_Pz_S_ac = I_ERI_I3y3z_S_Pz_S_ac+ABZ*I_ERI_H3y2z_S_Pz_S_ac;
  Double I_ERI_H2y3z_Pz_Pz_S_ac = I_ERI_I2y4z_S_Pz_S_ac+ABZ*I_ERI_H2y3z_S_Pz_S_ac;
  Double I_ERI_Hy4z_Pz_Pz_S_ac = I_ERI_Iy5z_S_Pz_S_ac+ABZ*I_ERI_Hy4z_S_Pz_S_ac;
  Double I_ERI_H5z_Pz_Pz_S_ac = I_ERI_I6z_S_Pz_S_ac+ABZ*I_ERI_H5z_S_Pz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_G_D_P_S_ac
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_P_P_S_ac
   * RHS shell quartet name: SQ_ERI_G_P_P_S_ac
   ************************************************************/
  Double I_ERI_G4x_D2x_Px_S_ac = I_ERI_H5x_Px_Px_S_ac+ABX*I_ERI_G4x_Px_Px_S_ac;
  Double I_ERI_G3xy_D2x_Px_S_ac = I_ERI_H4xy_Px_Px_S_ac+ABX*I_ERI_G3xy_Px_Px_S_ac;
  Double I_ERI_G3xz_D2x_Px_S_ac = I_ERI_H4xz_Px_Px_S_ac+ABX*I_ERI_G3xz_Px_Px_S_ac;
  Double I_ERI_G2x2y_D2x_Px_S_ac = I_ERI_H3x2y_Px_Px_S_ac+ABX*I_ERI_G2x2y_Px_Px_S_ac;
  Double I_ERI_G2xyz_D2x_Px_S_ac = I_ERI_H3xyz_Px_Px_S_ac+ABX*I_ERI_G2xyz_Px_Px_S_ac;
  Double I_ERI_G2x2z_D2x_Px_S_ac = I_ERI_H3x2z_Px_Px_S_ac+ABX*I_ERI_G2x2z_Px_Px_S_ac;
  Double I_ERI_Gx3y_D2x_Px_S_ac = I_ERI_H2x3y_Px_Px_S_ac+ABX*I_ERI_Gx3y_Px_Px_S_ac;
  Double I_ERI_Gx2yz_D2x_Px_S_ac = I_ERI_H2x2yz_Px_Px_S_ac+ABX*I_ERI_Gx2yz_Px_Px_S_ac;
  Double I_ERI_Gxy2z_D2x_Px_S_ac = I_ERI_H2xy2z_Px_Px_S_ac+ABX*I_ERI_Gxy2z_Px_Px_S_ac;
  Double I_ERI_Gx3z_D2x_Px_S_ac = I_ERI_H2x3z_Px_Px_S_ac+ABX*I_ERI_Gx3z_Px_Px_S_ac;
  Double I_ERI_G4y_D2x_Px_S_ac = I_ERI_Hx4y_Px_Px_S_ac+ABX*I_ERI_G4y_Px_Px_S_ac;
  Double I_ERI_G3yz_D2x_Px_S_ac = I_ERI_Hx3yz_Px_Px_S_ac+ABX*I_ERI_G3yz_Px_Px_S_ac;
  Double I_ERI_G2y2z_D2x_Px_S_ac = I_ERI_Hx2y2z_Px_Px_S_ac+ABX*I_ERI_G2y2z_Px_Px_S_ac;
  Double I_ERI_Gy3z_D2x_Px_S_ac = I_ERI_Hxy3z_Px_Px_S_ac+ABX*I_ERI_Gy3z_Px_Px_S_ac;
  Double I_ERI_G4z_D2x_Px_S_ac = I_ERI_Hx4z_Px_Px_S_ac+ABX*I_ERI_G4z_Px_Px_S_ac;
  Double I_ERI_G4x_Dxy_Px_S_ac = I_ERI_H4xy_Px_Px_S_ac+ABY*I_ERI_G4x_Px_Px_S_ac;
  Double I_ERI_G3xy_Dxy_Px_S_ac = I_ERI_H3x2y_Px_Px_S_ac+ABY*I_ERI_G3xy_Px_Px_S_ac;
  Double I_ERI_G3xz_Dxy_Px_S_ac = I_ERI_H3xyz_Px_Px_S_ac+ABY*I_ERI_G3xz_Px_Px_S_ac;
  Double I_ERI_G2x2y_Dxy_Px_S_ac = I_ERI_H2x3y_Px_Px_S_ac+ABY*I_ERI_G2x2y_Px_Px_S_ac;
  Double I_ERI_G2xyz_Dxy_Px_S_ac = I_ERI_H2x2yz_Px_Px_S_ac+ABY*I_ERI_G2xyz_Px_Px_S_ac;
  Double I_ERI_G2x2z_Dxy_Px_S_ac = I_ERI_H2xy2z_Px_Px_S_ac+ABY*I_ERI_G2x2z_Px_Px_S_ac;
  Double I_ERI_Gx3y_Dxy_Px_S_ac = I_ERI_Hx4y_Px_Px_S_ac+ABY*I_ERI_Gx3y_Px_Px_S_ac;
  Double I_ERI_Gx2yz_Dxy_Px_S_ac = I_ERI_Hx3yz_Px_Px_S_ac+ABY*I_ERI_Gx2yz_Px_Px_S_ac;
  Double I_ERI_Gxy2z_Dxy_Px_S_ac = I_ERI_Hx2y2z_Px_Px_S_ac+ABY*I_ERI_Gxy2z_Px_Px_S_ac;
  Double I_ERI_Gx3z_Dxy_Px_S_ac = I_ERI_Hxy3z_Px_Px_S_ac+ABY*I_ERI_Gx3z_Px_Px_S_ac;
  Double I_ERI_G4y_Dxy_Px_S_ac = I_ERI_H5y_Px_Px_S_ac+ABY*I_ERI_G4y_Px_Px_S_ac;
  Double I_ERI_G3yz_Dxy_Px_S_ac = I_ERI_H4yz_Px_Px_S_ac+ABY*I_ERI_G3yz_Px_Px_S_ac;
  Double I_ERI_G2y2z_Dxy_Px_S_ac = I_ERI_H3y2z_Px_Px_S_ac+ABY*I_ERI_G2y2z_Px_Px_S_ac;
  Double I_ERI_Gy3z_Dxy_Px_S_ac = I_ERI_H2y3z_Px_Px_S_ac+ABY*I_ERI_Gy3z_Px_Px_S_ac;
  Double I_ERI_G4z_Dxy_Px_S_ac = I_ERI_Hy4z_Px_Px_S_ac+ABY*I_ERI_G4z_Px_Px_S_ac;
  Double I_ERI_G4x_Dxz_Px_S_ac = I_ERI_H4xz_Px_Px_S_ac+ABZ*I_ERI_G4x_Px_Px_S_ac;
  Double I_ERI_G3xy_Dxz_Px_S_ac = I_ERI_H3xyz_Px_Px_S_ac+ABZ*I_ERI_G3xy_Px_Px_S_ac;
  Double I_ERI_G3xz_Dxz_Px_S_ac = I_ERI_H3x2z_Px_Px_S_ac+ABZ*I_ERI_G3xz_Px_Px_S_ac;
  Double I_ERI_G2x2y_Dxz_Px_S_ac = I_ERI_H2x2yz_Px_Px_S_ac+ABZ*I_ERI_G2x2y_Px_Px_S_ac;
  Double I_ERI_G2xyz_Dxz_Px_S_ac = I_ERI_H2xy2z_Px_Px_S_ac+ABZ*I_ERI_G2xyz_Px_Px_S_ac;
  Double I_ERI_G2x2z_Dxz_Px_S_ac = I_ERI_H2x3z_Px_Px_S_ac+ABZ*I_ERI_G2x2z_Px_Px_S_ac;
  Double I_ERI_Gx3y_Dxz_Px_S_ac = I_ERI_Hx3yz_Px_Px_S_ac+ABZ*I_ERI_Gx3y_Px_Px_S_ac;
  Double I_ERI_Gx2yz_Dxz_Px_S_ac = I_ERI_Hx2y2z_Px_Px_S_ac+ABZ*I_ERI_Gx2yz_Px_Px_S_ac;
  Double I_ERI_Gxy2z_Dxz_Px_S_ac = I_ERI_Hxy3z_Px_Px_S_ac+ABZ*I_ERI_Gxy2z_Px_Px_S_ac;
  Double I_ERI_Gx3z_Dxz_Px_S_ac = I_ERI_Hx4z_Px_Px_S_ac+ABZ*I_ERI_Gx3z_Px_Px_S_ac;
  Double I_ERI_G4y_Dxz_Px_S_ac = I_ERI_H4yz_Px_Px_S_ac+ABZ*I_ERI_G4y_Px_Px_S_ac;
  Double I_ERI_G3yz_Dxz_Px_S_ac = I_ERI_H3y2z_Px_Px_S_ac+ABZ*I_ERI_G3yz_Px_Px_S_ac;
  Double I_ERI_G2y2z_Dxz_Px_S_ac = I_ERI_H2y3z_Px_Px_S_ac+ABZ*I_ERI_G2y2z_Px_Px_S_ac;
  Double I_ERI_Gy3z_Dxz_Px_S_ac = I_ERI_Hy4z_Px_Px_S_ac+ABZ*I_ERI_Gy3z_Px_Px_S_ac;
  Double I_ERI_G4z_Dxz_Px_S_ac = I_ERI_H5z_Px_Px_S_ac+ABZ*I_ERI_G4z_Px_Px_S_ac;
  Double I_ERI_G4x_D2y_Px_S_ac = I_ERI_H4xy_Py_Px_S_ac+ABY*I_ERI_G4x_Py_Px_S_ac;
  Double I_ERI_G3xy_D2y_Px_S_ac = I_ERI_H3x2y_Py_Px_S_ac+ABY*I_ERI_G3xy_Py_Px_S_ac;
  Double I_ERI_G3xz_D2y_Px_S_ac = I_ERI_H3xyz_Py_Px_S_ac+ABY*I_ERI_G3xz_Py_Px_S_ac;
  Double I_ERI_G2x2y_D2y_Px_S_ac = I_ERI_H2x3y_Py_Px_S_ac+ABY*I_ERI_G2x2y_Py_Px_S_ac;
  Double I_ERI_G2xyz_D2y_Px_S_ac = I_ERI_H2x2yz_Py_Px_S_ac+ABY*I_ERI_G2xyz_Py_Px_S_ac;
  Double I_ERI_G2x2z_D2y_Px_S_ac = I_ERI_H2xy2z_Py_Px_S_ac+ABY*I_ERI_G2x2z_Py_Px_S_ac;
  Double I_ERI_Gx3y_D2y_Px_S_ac = I_ERI_Hx4y_Py_Px_S_ac+ABY*I_ERI_Gx3y_Py_Px_S_ac;
  Double I_ERI_Gx2yz_D2y_Px_S_ac = I_ERI_Hx3yz_Py_Px_S_ac+ABY*I_ERI_Gx2yz_Py_Px_S_ac;
  Double I_ERI_Gxy2z_D2y_Px_S_ac = I_ERI_Hx2y2z_Py_Px_S_ac+ABY*I_ERI_Gxy2z_Py_Px_S_ac;
  Double I_ERI_Gx3z_D2y_Px_S_ac = I_ERI_Hxy3z_Py_Px_S_ac+ABY*I_ERI_Gx3z_Py_Px_S_ac;
  Double I_ERI_G4y_D2y_Px_S_ac = I_ERI_H5y_Py_Px_S_ac+ABY*I_ERI_G4y_Py_Px_S_ac;
  Double I_ERI_G3yz_D2y_Px_S_ac = I_ERI_H4yz_Py_Px_S_ac+ABY*I_ERI_G3yz_Py_Px_S_ac;
  Double I_ERI_G2y2z_D2y_Px_S_ac = I_ERI_H3y2z_Py_Px_S_ac+ABY*I_ERI_G2y2z_Py_Px_S_ac;
  Double I_ERI_Gy3z_D2y_Px_S_ac = I_ERI_H2y3z_Py_Px_S_ac+ABY*I_ERI_Gy3z_Py_Px_S_ac;
  Double I_ERI_G4z_D2y_Px_S_ac = I_ERI_Hy4z_Py_Px_S_ac+ABY*I_ERI_G4z_Py_Px_S_ac;
  Double I_ERI_G4x_Dyz_Px_S_ac = I_ERI_H4xz_Py_Px_S_ac+ABZ*I_ERI_G4x_Py_Px_S_ac;
  Double I_ERI_G3xy_Dyz_Px_S_ac = I_ERI_H3xyz_Py_Px_S_ac+ABZ*I_ERI_G3xy_Py_Px_S_ac;
  Double I_ERI_G3xz_Dyz_Px_S_ac = I_ERI_H3x2z_Py_Px_S_ac+ABZ*I_ERI_G3xz_Py_Px_S_ac;
  Double I_ERI_G2x2y_Dyz_Px_S_ac = I_ERI_H2x2yz_Py_Px_S_ac+ABZ*I_ERI_G2x2y_Py_Px_S_ac;
  Double I_ERI_G2xyz_Dyz_Px_S_ac = I_ERI_H2xy2z_Py_Px_S_ac+ABZ*I_ERI_G2xyz_Py_Px_S_ac;
  Double I_ERI_G2x2z_Dyz_Px_S_ac = I_ERI_H2x3z_Py_Px_S_ac+ABZ*I_ERI_G2x2z_Py_Px_S_ac;
  Double I_ERI_Gx3y_Dyz_Px_S_ac = I_ERI_Hx3yz_Py_Px_S_ac+ABZ*I_ERI_Gx3y_Py_Px_S_ac;
  Double I_ERI_Gx2yz_Dyz_Px_S_ac = I_ERI_Hx2y2z_Py_Px_S_ac+ABZ*I_ERI_Gx2yz_Py_Px_S_ac;
  Double I_ERI_Gxy2z_Dyz_Px_S_ac = I_ERI_Hxy3z_Py_Px_S_ac+ABZ*I_ERI_Gxy2z_Py_Px_S_ac;
  Double I_ERI_Gx3z_Dyz_Px_S_ac = I_ERI_Hx4z_Py_Px_S_ac+ABZ*I_ERI_Gx3z_Py_Px_S_ac;
  Double I_ERI_G4y_Dyz_Px_S_ac = I_ERI_H4yz_Py_Px_S_ac+ABZ*I_ERI_G4y_Py_Px_S_ac;
  Double I_ERI_G3yz_Dyz_Px_S_ac = I_ERI_H3y2z_Py_Px_S_ac+ABZ*I_ERI_G3yz_Py_Px_S_ac;
  Double I_ERI_G2y2z_Dyz_Px_S_ac = I_ERI_H2y3z_Py_Px_S_ac+ABZ*I_ERI_G2y2z_Py_Px_S_ac;
  Double I_ERI_Gy3z_Dyz_Px_S_ac = I_ERI_Hy4z_Py_Px_S_ac+ABZ*I_ERI_Gy3z_Py_Px_S_ac;
  Double I_ERI_G4z_Dyz_Px_S_ac = I_ERI_H5z_Py_Px_S_ac+ABZ*I_ERI_G4z_Py_Px_S_ac;
  Double I_ERI_G4x_D2z_Px_S_ac = I_ERI_H4xz_Pz_Px_S_ac+ABZ*I_ERI_G4x_Pz_Px_S_ac;
  Double I_ERI_G3xy_D2z_Px_S_ac = I_ERI_H3xyz_Pz_Px_S_ac+ABZ*I_ERI_G3xy_Pz_Px_S_ac;
  Double I_ERI_G3xz_D2z_Px_S_ac = I_ERI_H3x2z_Pz_Px_S_ac+ABZ*I_ERI_G3xz_Pz_Px_S_ac;
  Double I_ERI_G2x2y_D2z_Px_S_ac = I_ERI_H2x2yz_Pz_Px_S_ac+ABZ*I_ERI_G2x2y_Pz_Px_S_ac;
  Double I_ERI_G2xyz_D2z_Px_S_ac = I_ERI_H2xy2z_Pz_Px_S_ac+ABZ*I_ERI_G2xyz_Pz_Px_S_ac;
  Double I_ERI_G2x2z_D2z_Px_S_ac = I_ERI_H2x3z_Pz_Px_S_ac+ABZ*I_ERI_G2x2z_Pz_Px_S_ac;
  Double I_ERI_Gx3y_D2z_Px_S_ac = I_ERI_Hx3yz_Pz_Px_S_ac+ABZ*I_ERI_Gx3y_Pz_Px_S_ac;
  Double I_ERI_Gx2yz_D2z_Px_S_ac = I_ERI_Hx2y2z_Pz_Px_S_ac+ABZ*I_ERI_Gx2yz_Pz_Px_S_ac;
  Double I_ERI_Gxy2z_D2z_Px_S_ac = I_ERI_Hxy3z_Pz_Px_S_ac+ABZ*I_ERI_Gxy2z_Pz_Px_S_ac;
  Double I_ERI_Gx3z_D2z_Px_S_ac = I_ERI_Hx4z_Pz_Px_S_ac+ABZ*I_ERI_Gx3z_Pz_Px_S_ac;
  Double I_ERI_G4y_D2z_Px_S_ac = I_ERI_H4yz_Pz_Px_S_ac+ABZ*I_ERI_G4y_Pz_Px_S_ac;
  Double I_ERI_G3yz_D2z_Px_S_ac = I_ERI_H3y2z_Pz_Px_S_ac+ABZ*I_ERI_G3yz_Pz_Px_S_ac;
  Double I_ERI_G2y2z_D2z_Px_S_ac = I_ERI_H2y3z_Pz_Px_S_ac+ABZ*I_ERI_G2y2z_Pz_Px_S_ac;
  Double I_ERI_Gy3z_D2z_Px_S_ac = I_ERI_Hy4z_Pz_Px_S_ac+ABZ*I_ERI_Gy3z_Pz_Px_S_ac;
  Double I_ERI_G4z_D2z_Px_S_ac = I_ERI_H5z_Pz_Px_S_ac+ABZ*I_ERI_G4z_Pz_Px_S_ac;
  Double I_ERI_G4x_D2x_Py_S_ac = I_ERI_H5x_Px_Py_S_ac+ABX*I_ERI_G4x_Px_Py_S_ac;
  Double I_ERI_G3xy_D2x_Py_S_ac = I_ERI_H4xy_Px_Py_S_ac+ABX*I_ERI_G3xy_Px_Py_S_ac;
  Double I_ERI_G3xz_D2x_Py_S_ac = I_ERI_H4xz_Px_Py_S_ac+ABX*I_ERI_G3xz_Px_Py_S_ac;
  Double I_ERI_G2x2y_D2x_Py_S_ac = I_ERI_H3x2y_Px_Py_S_ac+ABX*I_ERI_G2x2y_Px_Py_S_ac;
  Double I_ERI_G2xyz_D2x_Py_S_ac = I_ERI_H3xyz_Px_Py_S_ac+ABX*I_ERI_G2xyz_Px_Py_S_ac;
  Double I_ERI_G2x2z_D2x_Py_S_ac = I_ERI_H3x2z_Px_Py_S_ac+ABX*I_ERI_G2x2z_Px_Py_S_ac;
  Double I_ERI_Gx3y_D2x_Py_S_ac = I_ERI_H2x3y_Px_Py_S_ac+ABX*I_ERI_Gx3y_Px_Py_S_ac;
  Double I_ERI_Gx2yz_D2x_Py_S_ac = I_ERI_H2x2yz_Px_Py_S_ac+ABX*I_ERI_Gx2yz_Px_Py_S_ac;
  Double I_ERI_Gxy2z_D2x_Py_S_ac = I_ERI_H2xy2z_Px_Py_S_ac+ABX*I_ERI_Gxy2z_Px_Py_S_ac;
  Double I_ERI_Gx3z_D2x_Py_S_ac = I_ERI_H2x3z_Px_Py_S_ac+ABX*I_ERI_Gx3z_Px_Py_S_ac;
  Double I_ERI_G4y_D2x_Py_S_ac = I_ERI_Hx4y_Px_Py_S_ac+ABX*I_ERI_G4y_Px_Py_S_ac;
  Double I_ERI_G3yz_D2x_Py_S_ac = I_ERI_Hx3yz_Px_Py_S_ac+ABX*I_ERI_G3yz_Px_Py_S_ac;
  Double I_ERI_G2y2z_D2x_Py_S_ac = I_ERI_Hx2y2z_Px_Py_S_ac+ABX*I_ERI_G2y2z_Px_Py_S_ac;
  Double I_ERI_Gy3z_D2x_Py_S_ac = I_ERI_Hxy3z_Px_Py_S_ac+ABX*I_ERI_Gy3z_Px_Py_S_ac;
  Double I_ERI_G4z_D2x_Py_S_ac = I_ERI_Hx4z_Px_Py_S_ac+ABX*I_ERI_G4z_Px_Py_S_ac;
  Double I_ERI_G4x_Dxy_Py_S_ac = I_ERI_H4xy_Px_Py_S_ac+ABY*I_ERI_G4x_Px_Py_S_ac;
  Double I_ERI_G3xy_Dxy_Py_S_ac = I_ERI_H3x2y_Px_Py_S_ac+ABY*I_ERI_G3xy_Px_Py_S_ac;
  Double I_ERI_G3xz_Dxy_Py_S_ac = I_ERI_H3xyz_Px_Py_S_ac+ABY*I_ERI_G3xz_Px_Py_S_ac;
  Double I_ERI_G2x2y_Dxy_Py_S_ac = I_ERI_H2x3y_Px_Py_S_ac+ABY*I_ERI_G2x2y_Px_Py_S_ac;
  Double I_ERI_G2xyz_Dxy_Py_S_ac = I_ERI_H2x2yz_Px_Py_S_ac+ABY*I_ERI_G2xyz_Px_Py_S_ac;
  Double I_ERI_G2x2z_Dxy_Py_S_ac = I_ERI_H2xy2z_Px_Py_S_ac+ABY*I_ERI_G2x2z_Px_Py_S_ac;
  Double I_ERI_Gx3y_Dxy_Py_S_ac = I_ERI_Hx4y_Px_Py_S_ac+ABY*I_ERI_Gx3y_Px_Py_S_ac;
  Double I_ERI_Gx2yz_Dxy_Py_S_ac = I_ERI_Hx3yz_Px_Py_S_ac+ABY*I_ERI_Gx2yz_Px_Py_S_ac;
  Double I_ERI_Gxy2z_Dxy_Py_S_ac = I_ERI_Hx2y2z_Px_Py_S_ac+ABY*I_ERI_Gxy2z_Px_Py_S_ac;
  Double I_ERI_Gx3z_Dxy_Py_S_ac = I_ERI_Hxy3z_Px_Py_S_ac+ABY*I_ERI_Gx3z_Px_Py_S_ac;
  Double I_ERI_G4y_Dxy_Py_S_ac = I_ERI_H5y_Px_Py_S_ac+ABY*I_ERI_G4y_Px_Py_S_ac;
  Double I_ERI_G3yz_Dxy_Py_S_ac = I_ERI_H4yz_Px_Py_S_ac+ABY*I_ERI_G3yz_Px_Py_S_ac;
  Double I_ERI_G2y2z_Dxy_Py_S_ac = I_ERI_H3y2z_Px_Py_S_ac+ABY*I_ERI_G2y2z_Px_Py_S_ac;
  Double I_ERI_Gy3z_Dxy_Py_S_ac = I_ERI_H2y3z_Px_Py_S_ac+ABY*I_ERI_Gy3z_Px_Py_S_ac;
  Double I_ERI_G4z_Dxy_Py_S_ac = I_ERI_Hy4z_Px_Py_S_ac+ABY*I_ERI_G4z_Px_Py_S_ac;
  Double I_ERI_G4x_Dxz_Py_S_ac = I_ERI_H4xz_Px_Py_S_ac+ABZ*I_ERI_G4x_Px_Py_S_ac;
  Double I_ERI_G3xy_Dxz_Py_S_ac = I_ERI_H3xyz_Px_Py_S_ac+ABZ*I_ERI_G3xy_Px_Py_S_ac;
  Double I_ERI_G3xz_Dxz_Py_S_ac = I_ERI_H3x2z_Px_Py_S_ac+ABZ*I_ERI_G3xz_Px_Py_S_ac;
  Double I_ERI_G2x2y_Dxz_Py_S_ac = I_ERI_H2x2yz_Px_Py_S_ac+ABZ*I_ERI_G2x2y_Px_Py_S_ac;
  Double I_ERI_G2xyz_Dxz_Py_S_ac = I_ERI_H2xy2z_Px_Py_S_ac+ABZ*I_ERI_G2xyz_Px_Py_S_ac;
  Double I_ERI_G2x2z_Dxz_Py_S_ac = I_ERI_H2x3z_Px_Py_S_ac+ABZ*I_ERI_G2x2z_Px_Py_S_ac;
  Double I_ERI_Gx3y_Dxz_Py_S_ac = I_ERI_Hx3yz_Px_Py_S_ac+ABZ*I_ERI_Gx3y_Px_Py_S_ac;
  Double I_ERI_Gx2yz_Dxz_Py_S_ac = I_ERI_Hx2y2z_Px_Py_S_ac+ABZ*I_ERI_Gx2yz_Px_Py_S_ac;
  Double I_ERI_Gxy2z_Dxz_Py_S_ac = I_ERI_Hxy3z_Px_Py_S_ac+ABZ*I_ERI_Gxy2z_Px_Py_S_ac;
  Double I_ERI_Gx3z_Dxz_Py_S_ac = I_ERI_Hx4z_Px_Py_S_ac+ABZ*I_ERI_Gx3z_Px_Py_S_ac;
  Double I_ERI_G4y_Dxz_Py_S_ac = I_ERI_H4yz_Px_Py_S_ac+ABZ*I_ERI_G4y_Px_Py_S_ac;
  Double I_ERI_G3yz_Dxz_Py_S_ac = I_ERI_H3y2z_Px_Py_S_ac+ABZ*I_ERI_G3yz_Px_Py_S_ac;
  Double I_ERI_G2y2z_Dxz_Py_S_ac = I_ERI_H2y3z_Px_Py_S_ac+ABZ*I_ERI_G2y2z_Px_Py_S_ac;
  Double I_ERI_Gy3z_Dxz_Py_S_ac = I_ERI_Hy4z_Px_Py_S_ac+ABZ*I_ERI_Gy3z_Px_Py_S_ac;
  Double I_ERI_G4z_Dxz_Py_S_ac = I_ERI_H5z_Px_Py_S_ac+ABZ*I_ERI_G4z_Px_Py_S_ac;
  Double I_ERI_G4x_D2y_Py_S_ac = I_ERI_H4xy_Py_Py_S_ac+ABY*I_ERI_G4x_Py_Py_S_ac;
  Double I_ERI_G3xy_D2y_Py_S_ac = I_ERI_H3x2y_Py_Py_S_ac+ABY*I_ERI_G3xy_Py_Py_S_ac;
  Double I_ERI_G3xz_D2y_Py_S_ac = I_ERI_H3xyz_Py_Py_S_ac+ABY*I_ERI_G3xz_Py_Py_S_ac;
  Double I_ERI_G2x2y_D2y_Py_S_ac = I_ERI_H2x3y_Py_Py_S_ac+ABY*I_ERI_G2x2y_Py_Py_S_ac;
  Double I_ERI_G2xyz_D2y_Py_S_ac = I_ERI_H2x2yz_Py_Py_S_ac+ABY*I_ERI_G2xyz_Py_Py_S_ac;
  Double I_ERI_G2x2z_D2y_Py_S_ac = I_ERI_H2xy2z_Py_Py_S_ac+ABY*I_ERI_G2x2z_Py_Py_S_ac;
  Double I_ERI_Gx3y_D2y_Py_S_ac = I_ERI_Hx4y_Py_Py_S_ac+ABY*I_ERI_Gx3y_Py_Py_S_ac;
  Double I_ERI_Gx2yz_D2y_Py_S_ac = I_ERI_Hx3yz_Py_Py_S_ac+ABY*I_ERI_Gx2yz_Py_Py_S_ac;
  Double I_ERI_Gxy2z_D2y_Py_S_ac = I_ERI_Hx2y2z_Py_Py_S_ac+ABY*I_ERI_Gxy2z_Py_Py_S_ac;
  Double I_ERI_Gx3z_D2y_Py_S_ac = I_ERI_Hxy3z_Py_Py_S_ac+ABY*I_ERI_Gx3z_Py_Py_S_ac;
  Double I_ERI_G4y_D2y_Py_S_ac = I_ERI_H5y_Py_Py_S_ac+ABY*I_ERI_G4y_Py_Py_S_ac;
  Double I_ERI_G3yz_D2y_Py_S_ac = I_ERI_H4yz_Py_Py_S_ac+ABY*I_ERI_G3yz_Py_Py_S_ac;
  Double I_ERI_G2y2z_D2y_Py_S_ac = I_ERI_H3y2z_Py_Py_S_ac+ABY*I_ERI_G2y2z_Py_Py_S_ac;
  Double I_ERI_Gy3z_D2y_Py_S_ac = I_ERI_H2y3z_Py_Py_S_ac+ABY*I_ERI_Gy3z_Py_Py_S_ac;
  Double I_ERI_G4z_D2y_Py_S_ac = I_ERI_Hy4z_Py_Py_S_ac+ABY*I_ERI_G4z_Py_Py_S_ac;
  Double I_ERI_G4x_Dyz_Py_S_ac = I_ERI_H4xz_Py_Py_S_ac+ABZ*I_ERI_G4x_Py_Py_S_ac;
  Double I_ERI_G3xy_Dyz_Py_S_ac = I_ERI_H3xyz_Py_Py_S_ac+ABZ*I_ERI_G3xy_Py_Py_S_ac;
  Double I_ERI_G3xz_Dyz_Py_S_ac = I_ERI_H3x2z_Py_Py_S_ac+ABZ*I_ERI_G3xz_Py_Py_S_ac;
  Double I_ERI_G2x2y_Dyz_Py_S_ac = I_ERI_H2x2yz_Py_Py_S_ac+ABZ*I_ERI_G2x2y_Py_Py_S_ac;
  Double I_ERI_G2xyz_Dyz_Py_S_ac = I_ERI_H2xy2z_Py_Py_S_ac+ABZ*I_ERI_G2xyz_Py_Py_S_ac;
  Double I_ERI_G2x2z_Dyz_Py_S_ac = I_ERI_H2x3z_Py_Py_S_ac+ABZ*I_ERI_G2x2z_Py_Py_S_ac;
  Double I_ERI_Gx3y_Dyz_Py_S_ac = I_ERI_Hx3yz_Py_Py_S_ac+ABZ*I_ERI_Gx3y_Py_Py_S_ac;
  Double I_ERI_Gx2yz_Dyz_Py_S_ac = I_ERI_Hx2y2z_Py_Py_S_ac+ABZ*I_ERI_Gx2yz_Py_Py_S_ac;
  Double I_ERI_Gxy2z_Dyz_Py_S_ac = I_ERI_Hxy3z_Py_Py_S_ac+ABZ*I_ERI_Gxy2z_Py_Py_S_ac;
  Double I_ERI_Gx3z_Dyz_Py_S_ac = I_ERI_Hx4z_Py_Py_S_ac+ABZ*I_ERI_Gx3z_Py_Py_S_ac;
  Double I_ERI_G4y_Dyz_Py_S_ac = I_ERI_H4yz_Py_Py_S_ac+ABZ*I_ERI_G4y_Py_Py_S_ac;
  Double I_ERI_G3yz_Dyz_Py_S_ac = I_ERI_H3y2z_Py_Py_S_ac+ABZ*I_ERI_G3yz_Py_Py_S_ac;
  Double I_ERI_G2y2z_Dyz_Py_S_ac = I_ERI_H2y3z_Py_Py_S_ac+ABZ*I_ERI_G2y2z_Py_Py_S_ac;
  Double I_ERI_Gy3z_Dyz_Py_S_ac = I_ERI_Hy4z_Py_Py_S_ac+ABZ*I_ERI_Gy3z_Py_Py_S_ac;
  Double I_ERI_G4z_Dyz_Py_S_ac = I_ERI_H5z_Py_Py_S_ac+ABZ*I_ERI_G4z_Py_Py_S_ac;
  Double I_ERI_G4x_D2z_Py_S_ac = I_ERI_H4xz_Pz_Py_S_ac+ABZ*I_ERI_G4x_Pz_Py_S_ac;
  Double I_ERI_G3xy_D2z_Py_S_ac = I_ERI_H3xyz_Pz_Py_S_ac+ABZ*I_ERI_G3xy_Pz_Py_S_ac;
  Double I_ERI_G3xz_D2z_Py_S_ac = I_ERI_H3x2z_Pz_Py_S_ac+ABZ*I_ERI_G3xz_Pz_Py_S_ac;
  Double I_ERI_G2x2y_D2z_Py_S_ac = I_ERI_H2x2yz_Pz_Py_S_ac+ABZ*I_ERI_G2x2y_Pz_Py_S_ac;
  Double I_ERI_G2xyz_D2z_Py_S_ac = I_ERI_H2xy2z_Pz_Py_S_ac+ABZ*I_ERI_G2xyz_Pz_Py_S_ac;
  Double I_ERI_G2x2z_D2z_Py_S_ac = I_ERI_H2x3z_Pz_Py_S_ac+ABZ*I_ERI_G2x2z_Pz_Py_S_ac;
  Double I_ERI_Gx3y_D2z_Py_S_ac = I_ERI_Hx3yz_Pz_Py_S_ac+ABZ*I_ERI_Gx3y_Pz_Py_S_ac;
  Double I_ERI_Gx2yz_D2z_Py_S_ac = I_ERI_Hx2y2z_Pz_Py_S_ac+ABZ*I_ERI_Gx2yz_Pz_Py_S_ac;
  Double I_ERI_Gxy2z_D2z_Py_S_ac = I_ERI_Hxy3z_Pz_Py_S_ac+ABZ*I_ERI_Gxy2z_Pz_Py_S_ac;
  Double I_ERI_Gx3z_D2z_Py_S_ac = I_ERI_Hx4z_Pz_Py_S_ac+ABZ*I_ERI_Gx3z_Pz_Py_S_ac;
  Double I_ERI_G4y_D2z_Py_S_ac = I_ERI_H4yz_Pz_Py_S_ac+ABZ*I_ERI_G4y_Pz_Py_S_ac;
  Double I_ERI_G3yz_D2z_Py_S_ac = I_ERI_H3y2z_Pz_Py_S_ac+ABZ*I_ERI_G3yz_Pz_Py_S_ac;
  Double I_ERI_G2y2z_D2z_Py_S_ac = I_ERI_H2y3z_Pz_Py_S_ac+ABZ*I_ERI_G2y2z_Pz_Py_S_ac;
  Double I_ERI_Gy3z_D2z_Py_S_ac = I_ERI_Hy4z_Pz_Py_S_ac+ABZ*I_ERI_Gy3z_Pz_Py_S_ac;
  Double I_ERI_G4z_D2z_Py_S_ac = I_ERI_H5z_Pz_Py_S_ac+ABZ*I_ERI_G4z_Pz_Py_S_ac;
  Double I_ERI_G4x_D2x_Pz_S_ac = I_ERI_H5x_Px_Pz_S_ac+ABX*I_ERI_G4x_Px_Pz_S_ac;
  Double I_ERI_G3xy_D2x_Pz_S_ac = I_ERI_H4xy_Px_Pz_S_ac+ABX*I_ERI_G3xy_Px_Pz_S_ac;
  Double I_ERI_G3xz_D2x_Pz_S_ac = I_ERI_H4xz_Px_Pz_S_ac+ABX*I_ERI_G3xz_Px_Pz_S_ac;
  Double I_ERI_G2x2y_D2x_Pz_S_ac = I_ERI_H3x2y_Px_Pz_S_ac+ABX*I_ERI_G2x2y_Px_Pz_S_ac;
  Double I_ERI_G2xyz_D2x_Pz_S_ac = I_ERI_H3xyz_Px_Pz_S_ac+ABX*I_ERI_G2xyz_Px_Pz_S_ac;
  Double I_ERI_G2x2z_D2x_Pz_S_ac = I_ERI_H3x2z_Px_Pz_S_ac+ABX*I_ERI_G2x2z_Px_Pz_S_ac;
  Double I_ERI_Gx3y_D2x_Pz_S_ac = I_ERI_H2x3y_Px_Pz_S_ac+ABX*I_ERI_Gx3y_Px_Pz_S_ac;
  Double I_ERI_Gx2yz_D2x_Pz_S_ac = I_ERI_H2x2yz_Px_Pz_S_ac+ABX*I_ERI_Gx2yz_Px_Pz_S_ac;
  Double I_ERI_Gxy2z_D2x_Pz_S_ac = I_ERI_H2xy2z_Px_Pz_S_ac+ABX*I_ERI_Gxy2z_Px_Pz_S_ac;
  Double I_ERI_Gx3z_D2x_Pz_S_ac = I_ERI_H2x3z_Px_Pz_S_ac+ABX*I_ERI_Gx3z_Px_Pz_S_ac;
  Double I_ERI_G4y_D2x_Pz_S_ac = I_ERI_Hx4y_Px_Pz_S_ac+ABX*I_ERI_G4y_Px_Pz_S_ac;
  Double I_ERI_G3yz_D2x_Pz_S_ac = I_ERI_Hx3yz_Px_Pz_S_ac+ABX*I_ERI_G3yz_Px_Pz_S_ac;
  Double I_ERI_G2y2z_D2x_Pz_S_ac = I_ERI_Hx2y2z_Px_Pz_S_ac+ABX*I_ERI_G2y2z_Px_Pz_S_ac;
  Double I_ERI_Gy3z_D2x_Pz_S_ac = I_ERI_Hxy3z_Px_Pz_S_ac+ABX*I_ERI_Gy3z_Px_Pz_S_ac;
  Double I_ERI_G4z_D2x_Pz_S_ac = I_ERI_Hx4z_Px_Pz_S_ac+ABX*I_ERI_G4z_Px_Pz_S_ac;
  Double I_ERI_G4x_Dxy_Pz_S_ac = I_ERI_H4xy_Px_Pz_S_ac+ABY*I_ERI_G4x_Px_Pz_S_ac;
  Double I_ERI_G3xy_Dxy_Pz_S_ac = I_ERI_H3x2y_Px_Pz_S_ac+ABY*I_ERI_G3xy_Px_Pz_S_ac;
  Double I_ERI_G3xz_Dxy_Pz_S_ac = I_ERI_H3xyz_Px_Pz_S_ac+ABY*I_ERI_G3xz_Px_Pz_S_ac;
  Double I_ERI_G2x2y_Dxy_Pz_S_ac = I_ERI_H2x3y_Px_Pz_S_ac+ABY*I_ERI_G2x2y_Px_Pz_S_ac;
  Double I_ERI_G2xyz_Dxy_Pz_S_ac = I_ERI_H2x2yz_Px_Pz_S_ac+ABY*I_ERI_G2xyz_Px_Pz_S_ac;
  Double I_ERI_G2x2z_Dxy_Pz_S_ac = I_ERI_H2xy2z_Px_Pz_S_ac+ABY*I_ERI_G2x2z_Px_Pz_S_ac;
  Double I_ERI_Gx3y_Dxy_Pz_S_ac = I_ERI_Hx4y_Px_Pz_S_ac+ABY*I_ERI_Gx3y_Px_Pz_S_ac;
  Double I_ERI_Gx2yz_Dxy_Pz_S_ac = I_ERI_Hx3yz_Px_Pz_S_ac+ABY*I_ERI_Gx2yz_Px_Pz_S_ac;
  Double I_ERI_Gxy2z_Dxy_Pz_S_ac = I_ERI_Hx2y2z_Px_Pz_S_ac+ABY*I_ERI_Gxy2z_Px_Pz_S_ac;
  Double I_ERI_Gx3z_Dxy_Pz_S_ac = I_ERI_Hxy3z_Px_Pz_S_ac+ABY*I_ERI_Gx3z_Px_Pz_S_ac;
  Double I_ERI_G4y_Dxy_Pz_S_ac = I_ERI_H5y_Px_Pz_S_ac+ABY*I_ERI_G4y_Px_Pz_S_ac;
  Double I_ERI_G3yz_Dxy_Pz_S_ac = I_ERI_H4yz_Px_Pz_S_ac+ABY*I_ERI_G3yz_Px_Pz_S_ac;
  Double I_ERI_G2y2z_Dxy_Pz_S_ac = I_ERI_H3y2z_Px_Pz_S_ac+ABY*I_ERI_G2y2z_Px_Pz_S_ac;
  Double I_ERI_Gy3z_Dxy_Pz_S_ac = I_ERI_H2y3z_Px_Pz_S_ac+ABY*I_ERI_Gy3z_Px_Pz_S_ac;
  Double I_ERI_G4z_Dxy_Pz_S_ac = I_ERI_Hy4z_Px_Pz_S_ac+ABY*I_ERI_G4z_Px_Pz_S_ac;
  Double I_ERI_G4x_Dxz_Pz_S_ac = I_ERI_H4xz_Px_Pz_S_ac+ABZ*I_ERI_G4x_Px_Pz_S_ac;
  Double I_ERI_G3xy_Dxz_Pz_S_ac = I_ERI_H3xyz_Px_Pz_S_ac+ABZ*I_ERI_G3xy_Px_Pz_S_ac;
  Double I_ERI_G3xz_Dxz_Pz_S_ac = I_ERI_H3x2z_Px_Pz_S_ac+ABZ*I_ERI_G3xz_Px_Pz_S_ac;
  Double I_ERI_G2x2y_Dxz_Pz_S_ac = I_ERI_H2x2yz_Px_Pz_S_ac+ABZ*I_ERI_G2x2y_Px_Pz_S_ac;
  Double I_ERI_G2xyz_Dxz_Pz_S_ac = I_ERI_H2xy2z_Px_Pz_S_ac+ABZ*I_ERI_G2xyz_Px_Pz_S_ac;
  Double I_ERI_G2x2z_Dxz_Pz_S_ac = I_ERI_H2x3z_Px_Pz_S_ac+ABZ*I_ERI_G2x2z_Px_Pz_S_ac;
  Double I_ERI_Gx3y_Dxz_Pz_S_ac = I_ERI_Hx3yz_Px_Pz_S_ac+ABZ*I_ERI_Gx3y_Px_Pz_S_ac;
  Double I_ERI_Gx2yz_Dxz_Pz_S_ac = I_ERI_Hx2y2z_Px_Pz_S_ac+ABZ*I_ERI_Gx2yz_Px_Pz_S_ac;
  Double I_ERI_Gxy2z_Dxz_Pz_S_ac = I_ERI_Hxy3z_Px_Pz_S_ac+ABZ*I_ERI_Gxy2z_Px_Pz_S_ac;
  Double I_ERI_Gx3z_Dxz_Pz_S_ac = I_ERI_Hx4z_Px_Pz_S_ac+ABZ*I_ERI_Gx3z_Px_Pz_S_ac;
  Double I_ERI_G4y_Dxz_Pz_S_ac = I_ERI_H4yz_Px_Pz_S_ac+ABZ*I_ERI_G4y_Px_Pz_S_ac;
  Double I_ERI_G3yz_Dxz_Pz_S_ac = I_ERI_H3y2z_Px_Pz_S_ac+ABZ*I_ERI_G3yz_Px_Pz_S_ac;
  Double I_ERI_G2y2z_Dxz_Pz_S_ac = I_ERI_H2y3z_Px_Pz_S_ac+ABZ*I_ERI_G2y2z_Px_Pz_S_ac;
  Double I_ERI_Gy3z_Dxz_Pz_S_ac = I_ERI_Hy4z_Px_Pz_S_ac+ABZ*I_ERI_Gy3z_Px_Pz_S_ac;
  Double I_ERI_G4z_Dxz_Pz_S_ac = I_ERI_H5z_Px_Pz_S_ac+ABZ*I_ERI_G4z_Px_Pz_S_ac;
  Double I_ERI_G4x_D2y_Pz_S_ac = I_ERI_H4xy_Py_Pz_S_ac+ABY*I_ERI_G4x_Py_Pz_S_ac;
  Double I_ERI_G3xy_D2y_Pz_S_ac = I_ERI_H3x2y_Py_Pz_S_ac+ABY*I_ERI_G3xy_Py_Pz_S_ac;
  Double I_ERI_G3xz_D2y_Pz_S_ac = I_ERI_H3xyz_Py_Pz_S_ac+ABY*I_ERI_G3xz_Py_Pz_S_ac;
  Double I_ERI_G2x2y_D2y_Pz_S_ac = I_ERI_H2x3y_Py_Pz_S_ac+ABY*I_ERI_G2x2y_Py_Pz_S_ac;
  Double I_ERI_G2xyz_D2y_Pz_S_ac = I_ERI_H2x2yz_Py_Pz_S_ac+ABY*I_ERI_G2xyz_Py_Pz_S_ac;
  Double I_ERI_G2x2z_D2y_Pz_S_ac = I_ERI_H2xy2z_Py_Pz_S_ac+ABY*I_ERI_G2x2z_Py_Pz_S_ac;
  Double I_ERI_Gx3y_D2y_Pz_S_ac = I_ERI_Hx4y_Py_Pz_S_ac+ABY*I_ERI_Gx3y_Py_Pz_S_ac;
  Double I_ERI_Gx2yz_D2y_Pz_S_ac = I_ERI_Hx3yz_Py_Pz_S_ac+ABY*I_ERI_Gx2yz_Py_Pz_S_ac;
  Double I_ERI_Gxy2z_D2y_Pz_S_ac = I_ERI_Hx2y2z_Py_Pz_S_ac+ABY*I_ERI_Gxy2z_Py_Pz_S_ac;
  Double I_ERI_Gx3z_D2y_Pz_S_ac = I_ERI_Hxy3z_Py_Pz_S_ac+ABY*I_ERI_Gx3z_Py_Pz_S_ac;
  Double I_ERI_G4y_D2y_Pz_S_ac = I_ERI_H5y_Py_Pz_S_ac+ABY*I_ERI_G4y_Py_Pz_S_ac;
  Double I_ERI_G3yz_D2y_Pz_S_ac = I_ERI_H4yz_Py_Pz_S_ac+ABY*I_ERI_G3yz_Py_Pz_S_ac;
  Double I_ERI_G2y2z_D2y_Pz_S_ac = I_ERI_H3y2z_Py_Pz_S_ac+ABY*I_ERI_G2y2z_Py_Pz_S_ac;
  Double I_ERI_Gy3z_D2y_Pz_S_ac = I_ERI_H2y3z_Py_Pz_S_ac+ABY*I_ERI_Gy3z_Py_Pz_S_ac;
  Double I_ERI_G4z_D2y_Pz_S_ac = I_ERI_Hy4z_Py_Pz_S_ac+ABY*I_ERI_G4z_Py_Pz_S_ac;
  Double I_ERI_G4x_Dyz_Pz_S_ac = I_ERI_H4xz_Py_Pz_S_ac+ABZ*I_ERI_G4x_Py_Pz_S_ac;
  Double I_ERI_G3xy_Dyz_Pz_S_ac = I_ERI_H3xyz_Py_Pz_S_ac+ABZ*I_ERI_G3xy_Py_Pz_S_ac;
  Double I_ERI_G3xz_Dyz_Pz_S_ac = I_ERI_H3x2z_Py_Pz_S_ac+ABZ*I_ERI_G3xz_Py_Pz_S_ac;
  Double I_ERI_G2x2y_Dyz_Pz_S_ac = I_ERI_H2x2yz_Py_Pz_S_ac+ABZ*I_ERI_G2x2y_Py_Pz_S_ac;
  Double I_ERI_G2xyz_Dyz_Pz_S_ac = I_ERI_H2xy2z_Py_Pz_S_ac+ABZ*I_ERI_G2xyz_Py_Pz_S_ac;
  Double I_ERI_G2x2z_Dyz_Pz_S_ac = I_ERI_H2x3z_Py_Pz_S_ac+ABZ*I_ERI_G2x2z_Py_Pz_S_ac;
  Double I_ERI_Gx3y_Dyz_Pz_S_ac = I_ERI_Hx3yz_Py_Pz_S_ac+ABZ*I_ERI_Gx3y_Py_Pz_S_ac;
  Double I_ERI_Gx2yz_Dyz_Pz_S_ac = I_ERI_Hx2y2z_Py_Pz_S_ac+ABZ*I_ERI_Gx2yz_Py_Pz_S_ac;
  Double I_ERI_Gxy2z_Dyz_Pz_S_ac = I_ERI_Hxy3z_Py_Pz_S_ac+ABZ*I_ERI_Gxy2z_Py_Pz_S_ac;
  Double I_ERI_Gx3z_Dyz_Pz_S_ac = I_ERI_Hx4z_Py_Pz_S_ac+ABZ*I_ERI_Gx3z_Py_Pz_S_ac;
  Double I_ERI_G4y_Dyz_Pz_S_ac = I_ERI_H4yz_Py_Pz_S_ac+ABZ*I_ERI_G4y_Py_Pz_S_ac;
  Double I_ERI_G3yz_Dyz_Pz_S_ac = I_ERI_H3y2z_Py_Pz_S_ac+ABZ*I_ERI_G3yz_Py_Pz_S_ac;
  Double I_ERI_G2y2z_Dyz_Pz_S_ac = I_ERI_H2y3z_Py_Pz_S_ac+ABZ*I_ERI_G2y2z_Py_Pz_S_ac;
  Double I_ERI_Gy3z_Dyz_Pz_S_ac = I_ERI_Hy4z_Py_Pz_S_ac+ABZ*I_ERI_Gy3z_Py_Pz_S_ac;
  Double I_ERI_G4z_Dyz_Pz_S_ac = I_ERI_H5z_Py_Pz_S_ac+ABZ*I_ERI_G4z_Py_Pz_S_ac;
  Double I_ERI_G4x_D2z_Pz_S_ac = I_ERI_H4xz_Pz_Pz_S_ac+ABZ*I_ERI_G4x_Pz_Pz_S_ac;
  Double I_ERI_G3xy_D2z_Pz_S_ac = I_ERI_H3xyz_Pz_Pz_S_ac+ABZ*I_ERI_G3xy_Pz_Pz_S_ac;
  Double I_ERI_G3xz_D2z_Pz_S_ac = I_ERI_H3x2z_Pz_Pz_S_ac+ABZ*I_ERI_G3xz_Pz_Pz_S_ac;
  Double I_ERI_G2x2y_D2z_Pz_S_ac = I_ERI_H2x2yz_Pz_Pz_S_ac+ABZ*I_ERI_G2x2y_Pz_Pz_S_ac;
  Double I_ERI_G2xyz_D2z_Pz_S_ac = I_ERI_H2xy2z_Pz_Pz_S_ac+ABZ*I_ERI_G2xyz_Pz_Pz_S_ac;
  Double I_ERI_G2x2z_D2z_Pz_S_ac = I_ERI_H2x3z_Pz_Pz_S_ac+ABZ*I_ERI_G2x2z_Pz_Pz_S_ac;
  Double I_ERI_Gx3y_D2z_Pz_S_ac = I_ERI_Hx3yz_Pz_Pz_S_ac+ABZ*I_ERI_Gx3y_Pz_Pz_S_ac;
  Double I_ERI_Gx2yz_D2z_Pz_S_ac = I_ERI_Hx2y2z_Pz_Pz_S_ac+ABZ*I_ERI_Gx2yz_Pz_Pz_S_ac;
  Double I_ERI_Gxy2z_D2z_Pz_S_ac = I_ERI_Hxy3z_Pz_Pz_S_ac+ABZ*I_ERI_Gxy2z_Pz_Pz_S_ac;
  Double I_ERI_Gx3z_D2z_Pz_S_ac = I_ERI_Hx4z_Pz_Pz_S_ac+ABZ*I_ERI_Gx3z_Pz_Pz_S_ac;
  Double I_ERI_G4y_D2z_Pz_S_ac = I_ERI_H4yz_Pz_Pz_S_ac+ABZ*I_ERI_G4y_Pz_Pz_S_ac;
  Double I_ERI_G3yz_D2z_Pz_S_ac = I_ERI_H3y2z_Pz_Pz_S_ac+ABZ*I_ERI_G3yz_Pz_Pz_S_ac;
  Double I_ERI_G2y2z_D2z_Pz_S_ac = I_ERI_H2y3z_Pz_Pz_S_ac+ABZ*I_ERI_G2y2z_Pz_Pz_S_ac;
  Double I_ERI_Gy3z_D2z_Pz_S_ac = I_ERI_Hy4z_Pz_Pz_S_ac+ABZ*I_ERI_Gy3z_Pz_Pz_S_ac;
  Double I_ERI_G4z_D2z_Pz_S_ac = I_ERI_H5z_Pz_Pz_S_ac+ABZ*I_ERI_G4z_Pz_Pz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_S_S_S_bb
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_bb = I_ERI_G4x_S_S_S_bb+ABX*I_ERI_F3x_S_S_S_bb;
  Double I_ERI_F2xy_Px_S_S_bb = I_ERI_G3xy_S_S_S_bb+ABX*I_ERI_F2xy_S_S_S_bb;
  Double I_ERI_F2xz_Px_S_S_bb = I_ERI_G3xz_S_S_S_bb+ABX*I_ERI_F2xz_S_S_S_bb;
  Double I_ERI_Fx2y_Px_S_S_bb = I_ERI_G2x2y_S_S_S_bb+ABX*I_ERI_Fx2y_S_S_S_bb;
  Double I_ERI_Fxyz_Px_S_S_bb = I_ERI_G2xyz_S_S_S_bb+ABX*I_ERI_Fxyz_S_S_S_bb;
  Double I_ERI_Fx2z_Px_S_S_bb = I_ERI_G2x2z_S_S_S_bb+ABX*I_ERI_Fx2z_S_S_S_bb;
  Double I_ERI_F3y_Px_S_S_bb = I_ERI_Gx3y_S_S_S_bb+ABX*I_ERI_F3y_S_S_S_bb;
  Double I_ERI_F2yz_Px_S_S_bb = I_ERI_Gx2yz_S_S_S_bb+ABX*I_ERI_F2yz_S_S_S_bb;
  Double I_ERI_Fy2z_Px_S_S_bb = I_ERI_Gxy2z_S_S_S_bb+ABX*I_ERI_Fy2z_S_S_S_bb;
  Double I_ERI_F3z_Px_S_S_bb = I_ERI_Gx3z_S_S_S_bb+ABX*I_ERI_F3z_S_S_S_bb;
  Double I_ERI_F3x_Py_S_S_bb = I_ERI_G3xy_S_S_S_bb+ABY*I_ERI_F3x_S_S_S_bb;
  Double I_ERI_F2xy_Py_S_S_bb = I_ERI_G2x2y_S_S_S_bb+ABY*I_ERI_F2xy_S_S_S_bb;
  Double I_ERI_F2xz_Py_S_S_bb = I_ERI_G2xyz_S_S_S_bb+ABY*I_ERI_F2xz_S_S_S_bb;
  Double I_ERI_Fx2y_Py_S_S_bb = I_ERI_Gx3y_S_S_S_bb+ABY*I_ERI_Fx2y_S_S_S_bb;
  Double I_ERI_Fxyz_Py_S_S_bb = I_ERI_Gx2yz_S_S_S_bb+ABY*I_ERI_Fxyz_S_S_S_bb;
  Double I_ERI_Fx2z_Py_S_S_bb = I_ERI_Gxy2z_S_S_S_bb+ABY*I_ERI_Fx2z_S_S_S_bb;
  Double I_ERI_F3y_Py_S_S_bb = I_ERI_G4y_S_S_S_bb+ABY*I_ERI_F3y_S_S_S_bb;
  Double I_ERI_F2yz_Py_S_S_bb = I_ERI_G3yz_S_S_S_bb+ABY*I_ERI_F2yz_S_S_S_bb;
  Double I_ERI_Fy2z_Py_S_S_bb = I_ERI_G2y2z_S_S_S_bb+ABY*I_ERI_Fy2z_S_S_S_bb;
  Double I_ERI_F3z_Py_S_S_bb = I_ERI_Gy3z_S_S_S_bb+ABY*I_ERI_F3z_S_S_S_bb;
  Double I_ERI_F3x_Pz_S_S_bb = I_ERI_G3xz_S_S_S_bb+ABZ*I_ERI_F3x_S_S_S_bb;
  Double I_ERI_F2xy_Pz_S_S_bb = I_ERI_G2xyz_S_S_S_bb+ABZ*I_ERI_F2xy_S_S_S_bb;
  Double I_ERI_F2xz_Pz_S_S_bb = I_ERI_G2x2z_S_S_S_bb+ABZ*I_ERI_F2xz_S_S_S_bb;
  Double I_ERI_Fx2y_Pz_S_S_bb = I_ERI_Gx2yz_S_S_S_bb+ABZ*I_ERI_Fx2y_S_S_S_bb;
  Double I_ERI_Fxyz_Pz_S_S_bb = I_ERI_Gxy2z_S_S_S_bb+ABZ*I_ERI_Fxyz_S_S_S_bb;
  Double I_ERI_Fx2z_Pz_S_S_bb = I_ERI_Gx3z_S_S_S_bb+ABZ*I_ERI_Fx2z_S_S_S_bb;
  Double I_ERI_F3y_Pz_S_S_bb = I_ERI_G3yz_S_S_S_bb+ABZ*I_ERI_F3y_S_S_S_bb;
  Double I_ERI_F2yz_Pz_S_S_bb = I_ERI_G2y2z_S_S_S_bb+ABZ*I_ERI_F2yz_S_S_S_bb;
  Double I_ERI_Fy2z_Pz_S_S_bb = I_ERI_Gy3z_S_S_S_bb+ABZ*I_ERI_Fy2z_S_S_S_bb;
  Double I_ERI_F3z_Pz_S_S_bb = I_ERI_G4z_S_S_S_bb+ABZ*I_ERI_F3z_S_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_S_S_bb
   * RHS shell quartet name: SQ_ERI_G_S_S_S_bb
   ************************************************************/
  Double I_ERI_G4x_Px_S_S_bb = I_ERI_H5x_S_S_S_bb+ABX*I_ERI_G4x_S_S_S_bb;
  Double I_ERI_G3xy_Px_S_S_bb = I_ERI_H4xy_S_S_S_bb+ABX*I_ERI_G3xy_S_S_S_bb;
  Double I_ERI_G3xz_Px_S_S_bb = I_ERI_H4xz_S_S_S_bb+ABX*I_ERI_G3xz_S_S_S_bb;
  Double I_ERI_G2x2y_Px_S_S_bb = I_ERI_H3x2y_S_S_S_bb+ABX*I_ERI_G2x2y_S_S_S_bb;
  Double I_ERI_G2xyz_Px_S_S_bb = I_ERI_H3xyz_S_S_S_bb+ABX*I_ERI_G2xyz_S_S_S_bb;
  Double I_ERI_G2x2z_Px_S_S_bb = I_ERI_H3x2z_S_S_S_bb+ABX*I_ERI_G2x2z_S_S_S_bb;
  Double I_ERI_Gx3y_Px_S_S_bb = I_ERI_H2x3y_S_S_S_bb+ABX*I_ERI_Gx3y_S_S_S_bb;
  Double I_ERI_Gx2yz_Px_S_S_bb = I_ERI_H2x2yz_S_S_S_bb+ABX*I_ERI_Gx2yz_S_S_S_bb;
  Double I_ERI_Gxy2z_Px_S_S_bb = I_ERI_H2xy2z_S_S_S_bb+ABX*I_ERI_Gxy2z_S_S_S_bb;
  Double I_ERI_Gx3z_Px_S_S_bb = I_ERI_H2x3z_S_S_S_bb+ABX*I_ERI_Gx3z_S_S_S_bb;
  Double I_ERI_G4y_Px_S_S_bb = I_ERI_Hx4y_S_S_S_bb+ABX*I_ERI_G4y_S_S_S_bb;
  Double I_ERI_G3yz_Px_S_S_bb = I_ERI_Hx3yz_S_S_S_bb+ABX*I_ERI_G3yz_S_S_S_bb;
  Double I_ERI_G2y2z_Px_S_S_bb = I_ERI_Hx2y2z_S_S_S_bb+ABX*I_ERI_G2y2z_S_S_S_bb;
  Double I_ERI_Gy3z_Px_S_S_bb = I_ERI_Hxy3z_S_S_S_bb+ABX*I_ERI_Gy3z_S_S_S_bb;
  Double I_ERI_G4z_Px_S_S_bb = I_ERI_Hx4z_S_S_S_bb+ABX*I_ERI_G4z_S_S_S_bb;
  Double I_ERI_G4x_Py_S_S_bb = I_ERI_H4xy_S_S_S_bb+ABY*I_ERI_G4x_S_S_S_bb;
  Double I_ERI_G3xy_Py_S_S_bb = I_ERI_H3x2y_S_S_S_bb+ABY*I_ERI_G3xy_S_S_S_bb;
  Double I_ERI_G3xz_Py_S_S_bb = I_ERI_H3xyz_S_S_S_bb+ABY*I_ERI_G3xz_S_S_S_bb;
  Double I_ERI_G2x2y_Py_S_S_bb = I_ERI_H2x3y_S_S_S_bb+ABY*I_ERI_G2x2y_S_S_S_bb;
  Double I_ERI_G2xyz_Py_S_S_bb = I_ERI_H2x2yz_S_S_S_bb+ABY*I_ERI_G2xyz_S_S_S_bb;
  Double I_ERI_G2x2z_Py_S_S_bb = I_ERI_H2xy2z_S_S_S_bb+ABY*I_ERI_G2x2z_S_S_S_bb;
  Double I_ERI_Gx3y_Py_S_S_bb = I_ERI_Hx4y_S_S_S_bb+ABY*I_ERI_Gx3y_S_S_S_bb;
  Double I_ERI_Gx2yz_Py_S_S_bb = I_ERI_Hx3yz_S_S_S_bb+ABY*I_ERI_Gx2yz_S_S_S_bb;
  Double I_ERI_Gxy2z_Py_S_S_bb = I_ERI_Hx2y2z_S_S_S_bb+ABY*I_ERI_Gxy2z_S_S_S_bb;
  Double I_ERI_Gx3z_Py_S_S_bb = I_ERI_Hxy3z_S_S_S_bb+ABY*I_ERI_Gx3z_S_S_S_bb;
  Double I_ERI_G4y_Py_S_S_bb = I_ERI_H5y_S_S_S_bb+ABY*I_ERI_G4y_S_S_S_bb;
  Double I_ERI_G3yz_Py_S_S_bb = I_ERI_H4yz_S_S_S_bb+ABY*I_ERI_G3yz_S_S_S_bb;
  Double I_ERI_G2y2z_Py_S_S_bb = I_ERI_H3y2z_S_S_S_bb+ABY*I_ERI_G2y2z_S_S_S_bb;
  Double I_ERI_Gy3z_Py_S_S_bb = I_ERI_H2y3z_S_S_S_bb+ABY*I_ERI_Gy3z_S_S_S_bb;
  Double I_ERI_G4z_Py_S_S_bb = I_ERI_Hy4z_S_S_S_bb+ABY*I_ERI_G4z_S_S_S_bb;
  Double I_ERI_G4x_Pz_S_S_bb = I_ERI_H4xz_S_S_S_bb+ABZ*I_ERI_G4x_S_S_S_bb;
  Double I_ERI_G3xy_Pz_S_S_bb = I_ERI_H3xyz_S_S_S_bb+ABZ*I_ERI_G3xy_S_S_S_bb;
  Double I_ERI_G3xz_Pz_S_S_bb = I_ERI_H3x2z_S_S_S_bb+ABZ*I_ERI_G3xz_S_S_S_bb;
  Double I_ERI_G2x2y_Pz_S_S_bb = I_ERI_H2x2yz_S_S_S_bb+ABZ*I_ERI_G2x2y_S_S_S_bb;
  Double I_ERI_G2xyz_Pz_S_S_bb = I_ERI_H2xy2z_S_S_S_bb+ABZ*I_ERI_G2xyz_S_S_S_bb;
  Double I_ERI_G2x2z_Pz_S_S_bb = I_ERI_H2x3z_S_S_S_bb+ABZ*I_ERI_G2x2z_S_S_S_bb;
  Double I_ERI_Gx3y_Pz_S_S_bb = I_ERI_Hx3yz_S_S_S_bb+ABZ*I_ERI_Gx3y_S_S_S_bb;
  Double I_ERI_Gx2yz_Pz_S_S_bb = I_ERI_Hx2y2z_S_S_S_bb+ABZ*I_ERI_Gx2yz_S_S_S_bb;
  Double I_ERI_Gxy2z_Pz_S_S_bb = I_ERI_Hxy3z_S_S_S_bb+ABZ*I_ERI_Gxy2z_S_S_S_bb;
  Double I_ERI_Gx3z_Pz_S_S_bb = I_ERI_Hx4z_S_S_S_bb+ABZ*I_ERI_Gx3z_S_S_S_bb;
  Double I_ERI_G4y_Pz_S_S_bb = I_ERI_H4yz_S_S_S_bb+ABZ*I_ERI_G4y_S_S_S_bb;
  Double I_ERI_G3yz_Pz_S_S_bb = I_ERI_H3y2z_S_S_S_bb+ABZ*I_ERI_G3yz_S_S_S_bb;
  Double I_ERI_G2y2z_Pz_S_S_bb = I_ERI_H2y3z_S_S_S_bb+ABZ*I_ERI_G2y2z_S_S_S_bb;
  Double I_ERI_Gy3z_Pz_S_S_bb = I_ERI_Hy4z_S_S_S_bb+ABZ*I_ERI_Gy3z_S_S_S_bb;
  Double I_ERI_G4z_Pz_S_S_bb = I_ERI_H5z_S_S_S_bb+ABZ*I_ERI_G4z_S_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_P_S_S_bb
   ************************************************************/
  Double I_ERI_F3x_D2x_S_S_bb = I_ERI_G4x_Px_S_S_bb+ABX*I_ERI_F3x_Px_S_S_bb;
  Double I_ERI_F2xy_D2x_S_S_bb = I_ERI_G3xy_Px_S_S_bb+ABX*I_ERI_F2xy_Px_S_S_bb;
  Double I_ERI_F2xz_D2x_S_S_bb = I_ERI_G3xz_Px_S_S_bb+ABX*I_ERI_F2xz_Px_S_S_bb;
  Double I_ERI_Fx2y_D2x_S_S_bb = I_ERI_G2x2y_Px_S_S_bb+ABX*I_ERI_Fx2y_Px_S_S_bb;
  Double I_ERI_Fxyz_D2x_S_S_bb = I_ERI_G2xyz_Px_S_S_bb+ABX*I_ERI_Fxyz_Px_S_S_bb;
  Double I_ERI_Fx2z_D2x_S_S_bb = I_ERI_G2x2z_Px_S_S_bb+ABX*I_ERI_Fx2z_Px_S_S_bb;
  Double I_ERI_F3y_D2x_S_S_bb = I_ERI_Gx3y_Px_S_S_bb+ABX*I_ERI_F3y_Px_S_S_bb;
  Double I_ERI_F2yz_D2x_S_S_bb = I_ERI_Gx2yz_Px_S_S_bb+ABX*I_ERI_F2yz_Px_S_S_bb;
  Double I_ERI_Fy2z_D2x_S_S_bb = I_ERI_Gxy2z_Px_S_S_bb+ABX*I_ERI_Fy2z_Px_S_S_bb;
  Double I_ERI_F3z_D2x_S_S_bb = I_ERI_Gx3z_Px_S_S_bb+ABX*I_ERI_F3z_Px_S_S_bb;
  Double I_ERI_F3x_D2y_S_S_bb = I_ERI_G3xy_Py_S_S_bb+ABY*I_ERI_F3x_Py_S_S_bb;
  Double I_ERI_F2xy_D2y_S_S_bb = I_ERI_G2x2y_Py_S_S_bb+ABY*I_ERI_F2xy_Py_S_S_bb;
  Double I_ERI_F2xz_D2y_S_S_bb = I_ERI_G2xyz_Py_S_S_bb+ABY*I_ERI_F2xz_Py_S_S_bb;
  Double I_ERI_Fx2y_D2y_S_S_bb = I_ERI_Gx3y_Py_S_S_bb+ABY*I_ERI_Fx2y_Py_S_S_bb;
  Double I_ERI_Fxyz_D2y_S_S_bb = I_ERI_Gx2yz_Py_S_S_bb+ABY*I_ERI_Fxyz_Py_S_S_bb;
  Double I_ERI_Fx2z_D2y_S_S_bb = I_ERI_Gxy2z_Py_S_S_bb+ABY*I_ERI_Fx2z_Py_S_S_bb;
  Double I_ERI_F3y_D2y_S_S_bb = I_ERI_G4y_Py_S_S_bb+ABY*I_ERI_F3y_Py_S_S_bb;
  Double I_ERI_F2yz_D2y_S_S_bb = I_ERI_G3yz_Py_S_S_bb+ABY*I_ERI_F2yz_Py_S_S_bb;
  Double I_ERI_Fy2z_D2y_S_S_bb = I_ERI_G2y2z_Py_S_S_bb+ABY*I_ERI_Fy2z_Py_S_S_bb;
  Double I_ERI_F3z_D2y_S_S_bb = I_ERI_Gy3z_Py_S_S_bb+ABY*I_ERI_F3z_Py_S_S_bb;
  Double I_ERI_F3x_D2z_S_S_bb = I_ERI_G3xz_Pz_S_S_bb+ABZ*I_ERI_F3x_Pz_S_S_bb;
  Double I_ERI_F2xy_D2z_S_S_bb = I_ERI_G2xyz_Pz_S_S_bb+ABZ*I_ERI_F2xy_Pz_S_S_bb;
  Double I_ERI_F2xz_D2z_S_S_bb = I_ERI_G2x2z_Pz_S_S_bb+ABZ*I_ERI_F2xz_Pz_S_S_bb;
  Double I_ERI_Fx2y_D2z_S_S_bb = I_ERI_Gx2yz_Pz_S_S_bb+ABZ*I_ERI_Fx2y_Pz_S_S_bb;
  Double I_ERI_Fxyz_D2z_S_S_bb = I_ERI_Gxy2z_Pz_S_S_bb+ABZ*I_ERI_Fxyz_Pz_S_S_bb;
  Double I_ERI_Fx2z_D2z_S_S_bb = I_ERI_Gx3z_Pz_S_S_bb+ABZ*I_ERI_Fx2z_Pz_S_S_bb;
  Double I_ERI_F3y_D2z_S_S_bb = I_ERI_G3yz_Pz_S_S_bb+ABZ*I_ERI_F3y_Pz_S_S_bb;
  Double I_ERI_F2yz_D2z_S_S_bb = I_ERI_G2y2z_Pz_S_S_bb+ABZ*I_ERI_F2yz_Pz_S_S_bb;
  Double I_ERI_Fy2z_D2z_S_S_bb = I_ERI_Gy3z_Pz_S_S_bb+ABZ*I_ERI_Fy2z_Pz_S_S_bb;
  Double I_ERI_F3z_D2z_S_S_bb = I_ERI_G4z_Pz_S_S_bb+ABZ*I_ERI_F3z_Pz_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_H_P_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_S_S_bb
   * RHS shell quartet name: SQ_ERI_H_S_S_S_bb
   ************************************************************/
  Double I_ERI_H5x_Px_S_S_bb = I_ERI_I6x_S_S_S_bb+ABX*I_ERI_H5x_S_S_S_bb;
  Double I_ERI_H4xy_Px_S_S_bb = I_ERI_I5xy_S_S_S_bb+ABX*I_ERI_H4xy_S_S_S_bb;
  Double I_ERI_H4xz_Px_S_S_bb = I_ERI_I5xz_S_S_S_bb+ABX*I_ERI_H4xz_S_S_S_bb;
  Double I_ERI_H3x2y_Px_S_S_bb = I_ERI_I4x2y_S_S_S_bb+ABX*I_ERI_H3x2y_S_S_S_bb;
  Double I_ERI_H3xyz_Px_S_S_bb = I_ERI_I4xyz_S_S_S_bb+ABX*I_ERI_H3xyz_S_S_S_bb;
  Double I_ERI_H3x2z_Px_S_S_bb = I_ERI_I4x2z_S_S_S_bb+ABX*I_ERI_H3x2z_S_S_S_bb;
  Double I_ERI_H2x3y_Px_S_S_bb = I_ERI_I3x3y_S_S_S_bb+ABX*I_ERI_H2x3y_S_S_S_bb;
  Double I_ERI_H2x2yz_Px_S_S_bb = I_ERI_I3x2yz_S_S_S_bb+ABX*I_ERI_H2x2yz_S_S_S_bb;
  Double I_ERI_H2xy2z_Px_S_S_bb = I_ERI_I3xy2z_S_S_S_bb+ABX*I_ERI_H2xy2z_S_S_S_bb;
  Double I_ERI_H2x3z_Px_S_S_bb = I_ERI_I3x3z_S_S_S_bb+ABX*I_ERI_H2x3z_S_S_S_bb;
  Double I_ERI_Hx4y_Px_S_S_bb = I_ERI_I2x4y_S_S_S_bb+ABX*I_ERI_Hx4y_S_S_S_bb;
  Double I_ERI_Hx3yz_Px_S_S_bb = I_ERI_I2x3yz_S_S_S_bb+ABX*I_ERI_Hx3yz_S_S_S_bb;
  Double I_ERI_Hx2y2z_Px_S_S_bb = I_ERI_I2x2y2z_S_S_S_bb+ABX*I_ERI_Hx2y2z_S_S_S_bb;
  Double I_ERI_Hxy3z_Px_S_S_bb = I_ERI_I2xy3z_S_S_S_bb+ABX*I_ERI_Hxy3z_S_S_S_bb;
  Double I_ERI_Hx4z_Px_S_S_bb = I_ERI_I2x4z_S_S_S_bb+ABX*I_ERI_Hx4z_S_S_S_bb;
  Double I_ERI_H5y_Px_S_S_bb = I_ERI_Ix5y_S_S_S_bb+ABX*I_ERI_H5y_S_S_S_bb;
  Double I_ERI_H4yz_Px_S_S_bb = I_ERI_Ix4yz_S_S_S_bb+ABX*I_ERI_H4yz_S_S_S_bb;
  Double I_ERI_H3y2z_Px_S_S_bb = I_ERI_Ix3y2z_S_S_S_bb+ABX*I_ERI_H3y2z_S_S_S_bb;
  Double I_ERI_H2y3z_Px_S_S_bb = I_ERI_Ix2y3z_S_S_S_bb+ABX*I_ERI_H2y3z_S_S_S_bb;
  Double I_ERI_Hy4z_Px_S_S_bb = I_ERI_Ixy4z_S_S_S_bb+ABX*I_ERI_Hy4z_S_S_S_bb;
  Double I_ERI_H5z_Px_S_S_bb = I_ERI_Ix5z_S_S_S_bb+ABX*I_ERI_H5z_S_S_S_bb;
  Double I_ERI_H4xy_Py_S_S_bb = I_ERI_I4x2y_S_S_S_bb+ABY*I_ERI_H4xy_S_S_S_bb;
  Double I_ERI_H4xz_Py_S_S_bb = I_ERI_I4xyz_S_S_S_bb+ABY*I_ERI_H4xz_S_S_S_bb;
  Double I_ERI_H3x2y_Py_S_S_bb = I_ERI_I3x3y_S_S_S_bb+ABY*I_ERI_H3x2y_S_S_S_bb;
  Double I_ERI_H3xyz_Py_S_S_bb = I_ERI_I3x2yz_S_S_S_bb+ABY*I_ERI_H3xyz_S_S_S_bb;
  Double I_ERI_H3x2z_Py_S_S_bb = I_ERI_I3xy2z_S_S_S_bb+ABY*I_ERI_H3x2z_S_S_S_bb;
  Double I_ERI_H2x3y_Py_S_S_bb = I_ERI_I2x4y_S_S_S_bb+ABY*I_ERI_H2x3y_S_S_S_bb;
  Double I_ERI_H2x2yz_Py_S_S_bb = I_ERI_I2x3yz_S_S_S_bb+ABY*I_ERI_H2x2yz_S_S_S_bb;
  Double I_ERI_H2xy2z_Py_S_S_bb = I_ERI_I2x2y2z_S_S_S_bb+ABY*I_ERI_H2xy2z_S_S_S_bb;
  Double I_ERI_H2x3z_Py_S_S_bb = I_ERI_I2xy3z_S_S_S_bb+ABY*I_ERI_H2x3z_S_S_S_bb;
  Double I_ERI_Hx4y_Py_S_S_bb = I_ERI_Ix5y_S_S_S_bb+ABY*I_ERI_Hx4y_S_S_S_bb;
  Double I_ERI_Hx3yz_Py_S_S_bb = I_ERI_Ix4yz_S_S_S_bb+ABY*I_ERI_Hx3yz_S_S_S_bb;
  Double I_ERI_Hx2y2z_Py_S_S_bb = I_ERI_Ix3y2z_S_S_S_bb+ABY*I_ERI_Hx2y2z_S_S_S_bb;
  Double I_ERI_Hxy3z_Py_S_S_bb = I_ERI_Ix2y3z_S_S_S_bb+ABY*I_ERI_Hxy3z_S_S_S_bb;
  Double I_ERI_Hx4z_Py_S_S_bb = I_ERI_Ixy4z_S_S_S_bb+ABY*I_ERI_Hx4z_S_S_S_bb;
  Double I_ERI_H5y_Py_S_S_bb = I_ERI_I6y_S_S_S_bb+ABY*I_ERI_H5y_S_S_S_bb;
  Double I_ERI_H4yz_Py_S_S_bb = I_ERI_I5yz_S_S_S_bb+ABY*I_ERI_H4yz_S_S_S_bb;
  Double I_ERI_H3y2z_Py_S_S_bb = I_ERI_I4y2z_S_S_S_bb+ABY*I_ERI_H3y2z_S_S_S_bb;
  Double I_ERI_H2y3z_Py_S_S_bb = I_ERI_I3y3z_S_S_S_bb+ABY*I_ERI_H2y3z_S_S_S_bb;
  Double I_ERI_Hy4z_Py_S_S_bb = I_ERI_I2y4z_S_S_S_bb+ABY*I_ERI_Hy4z_S_S_S_bb;
  Double I_ERI_H5z_Py_S_S_bb = I_ERI_Iy5z_S_S_S_bb+ABY*I_ERI_H5z_S_S_S_bb;
  Double I_ERI_H4xy_Pz_S_S_bb = I_ERI_I4xyz_S_S_S_bb+ABZ*I_ERI_H4xy_S_S_S_bb;
  Double I_ERI_H4xz_Pz_S_S_bb = I_ERI_I4x2z_S_S_S_bb+ABZ*I_ERI_H4xz_S_S_S_bb;
  Double I_ERI_H3x2y_Pz_S_S_bb = I_ERI_I3x2yz_S_S_S_bb+ABZ*I_ERI_H3x2y_S_S_S_bb;
  Double I_ERI_H3xyz_Pz_S_S_bb = I_ERI_I3xy2z_S_S_S_bb+ABZ*I_ERI_H3xyz_S_S_S_bb;
  Double I_ERI_H3x2z_Pz_S_S_bb = I_ERI_I3x3z_S_S_S_bb+ABZ*I_ERI_H3x2z_S_S_S_bb;
  Double I_ERI_H2x3y_Pz_S_S_bb = I_ERI_I2x3yz_S_S_S_bb+ABZ*I_ERI_H2x3y_S_S_S_bb;
  Double I_ERI_H2x2yz_Pz_S_S_bb = I_ERI_I2x2y2z_S_S_S_bb+ABZ*I_ERI_H2x2yz_S_S_S_bb;
  Double I_ERI_H2xy2z_Pz_S_S_bb = I_ERI_I2xy3z_S_S_S_bb+ABZ*I_ERI_H2xy2z_S_S_S_bb;
  Double I_ERI_H2x3z_Pz_S_S_bb = I_ERI_I2x4z_S_S_S_bb+ABZ*I_ERI_H2x3z_S_S_S_bb;
  Double I_ERI_Hx4y_Pz_S_S_bb = I_ERI_Ix4yz_S_S_S_bb+ABZ*I_ERI_Hx4y_S_S_S_bb;
  Double I_ERI_Hx3yz_Pz_S_S_bb = I_ERI_Ix3y2z_S_S_S_bb+ABZ*I_ERI_Hx3yz_S_S_S_bb;
  Double I_ERI_Hx2y2z_Pz_S_S_bb = I_ERI_Ix2y3z_S_S_S_bb+ABZ*I_ERI_Hx2y2z_S_S_S_bb;
  Double I_ERI_Hxy3z_Pz_S_S_bb = I_ERI_Ixy4z_S_S_S_bb+ABZ*I_ERI_Hxy3z_S_S_S_bb;
  Double I_ERI_Hx4z_Pz_S_S_bb = I_ERI_Ix5z_S_S_S_bb+ABZ*I_ERI_Hx4z_S_S_S_bb;
  Double I_ERI_H4yz_Pz_S_S_bb = I_ERI_I4y2z_S_S_S_bb+ABZ*I_ERI_H4yz_S_S_S_bb;
  Double I_ERI_H3y2z_Pz_S_S_bb = I_ERI_I3y3z_S_S_S_bb+ABZ*I_ERI_H3y2z_S_S_S_bb;
  Double I_ERI_H2y3z_Pz_S_S_bb = I_ERI_I2y4z_S_S_S_bb+ABZ*I_ERI_H2y3z_S_S_S_bb;
  Double I_ERI_Hy4z_Pz_S_S_bb = I_ERI_Iy5z_S_S_S_bb+ABZ*I_ERI_Hy4z_S_S_S_bb;
  Double I_ERI_H5z_Pz_S_S_bb = I_ERI_I6z_S_S_S_bb+ABZ*I_ERI_H5z_S_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_G_D_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 45 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_P_S_S_bb
   * RHS shell quartet name: SQ_ERI_G_P_S_S_bb
   ************************************************************/
  Double I_ERI_G4x_D2x_S_S_bb = I_ERI_H5x_Px_S_S_bb+ABX*I_ERI_G4x_Px_S_S_bb;
  Double I_ERI_G3xy_D2x_S_S_bb = I_ERI_H4xy_Px_S_S_bb+ABX*I_ERI_G3xy_Px_S_S_bb;
  Double I_ERI_G3xz_D2x_S_S_bb = I_ERI_H4xz_Px_S_S_bb+ABX*I_ERI_G3xz_Px_S_S_bb;
  Double I_ERI_G2x2y_D2x_S_S_bb = I_ERI_H3x2y_Px_S_S_bb+ABX*I_ERI_G2x2y_Px_S_S_bb;
  Double I_ERI_G2xyz_D2x_S_S_bb = I_ERI_H3xyz_Px_S_S_bb+ABX*I_ERI_G2xyz_Px_S_S_bb;
  Double I_ERI_G2x2z_D2x_S_S_bb = I_ERI_H3x2z_Px_S_S_bb+ABX*I_ERI_G2x2z_Px_S_S_bb;
  Double I_ERI_Gx3y_D2x_S_S_bb = I_ERI_H2x3y_Px_S_S_bb+ABX*I_ERI_Gx3y_Px_S_S_bb;
  Double I_ERI_Gx2yz_D2x_S_S_bb = I_ERI_H2x2yz_Px_S_S_bb+ABX*I_ERI_Gx2yz_Px_S_S_bb;
  Double I_ERI_Gxy2z_D2x_S_S_bb = I_ERI_H2xy2z_Px_S_S_bb+ABX*I_ERI_Gxy2z_Px_S_S_bb;
  Double I_ERI_Gx3z_D2x_S_S_bb = I_ERI_H2x3z_Px_S_S_bb+ABX*I_ERI_Gx3z_Px_S_S_bb;
  Double I_ERI_G4y_D2x_S_S_bb = I_ERI_Hx4y_Px_S_S_bb+ABX*I_ERI_G4y_Px_S_S_bb;
  Double I_ERI_G3yz_D2x_S_S_bb = I_ERI_Hx3yz_Px_S_S_bb+ABX*I_ERI_G3yz_Px_S_S_bb;
  Double I_ERI_G2y2z_D2x_S_S_bb = I_ERI_Hx2y2z_Px_S_S_bb+ABX*I_ERI_G2y2z_Px_S_S_bb;
  Double I_ERI_Gy3z_D2x_S_S_bb = I_ERI_Hxy3z_Px_S_S_bb+ABX*I_ERI_Gy3z_Px_S_S_bb;
  Double I_ERI_G4z_D2x_S_S_bb = I_ERI_Hx4z_Px_S_S_bb+ABX*I_ERI_G4z_Px_S_S_bb;
  Double I_ERI_G4x_D2y_S_S_bb = I_ERI_H4xy_Py_S_S_bb+ABY*I_ERI_G4x_Py_S_S_bb;
  Double I_ERI_G3xy_D2y_S_S_bb = I_ERI_H3x2y_Py_S_S_bb+ABY*I_ERI_G3xy_Py_S_S_bb;
  Double I_ERI_G3xz_D2y_S_S_bb = I_ERI_H3xyz_Py_S_S_bb+ABY*I_ERI_G3xz_Py_S_S_bb;
  Double I_ERI_G2x2y_D2y_S_S_bb = I_ERI_H2x3y_Py_S_S_bb+ABY*I_ERI_G2x2y_Py_S_S_bb;
  Double I_ERI_G2xyz_D2y_S_S_bb = I_ERI_H2x2yz_Py_S_S_bb+ABY*I_ERI_G2xyz_Py_S_S_bb;
  Double I_ERI_G2x2z_D2y_S_S_bb = I_ERI_H2xy2z_Py_S_S_bb+ABY*I_ERI_G2x2z_Py_S_S_bb;
  Double I_ERI_Gx3y_D2y_S_S_bb = I_ERI_Hx4y_Py_S_S_bb+ABY*I_ERI_Gx3y_Py_S_S_bb;
  Double I_ERI_Gx2yz_D2y_S_S_bb = I_ERI_Hx3yz_Py_S_S_bb+ABY*I_ERI_Gx2yz_Py_S_S_bb;
  Double I_ERI_Gxy2z_D2y_S_S_bb = I_ERI_Hx2y2z_Py_S_S_bb+ABY*I_ERI_Gxy2z_Py_S_S_bb;
  Double I_ERI_Gx3z_D2y_S_S_bb = I_ERI_Hxy3z_Py_S_S_bb+ABY*I_ERI_Gx3z_Py_S_S_bb;
  Double I_ERI_G4y_D2y_S_S_bb = I_ERI_H5y_Py_S_S_bb+ABY*I_ERI_G4y_Py_S_S_bb;
  Double I_ERI_G3yz_D2y_S_S_bb = I_ERI_H4yz_Py_S_S_bb+ABY*I_ERI_G3yz_Py_S_S_bb;
  Double I_ERI_G2y2z_D2y_S_S_bb = I_ERI_H3y2z_Py_S_S_bb+ABY*I_ERI_G2y2z_Py_S_S_bb;
  Double I_ERI_Gy3z_D2y_S_S_bb = I_ERI_H2y3z_Py_S_S_bb+ABY*I_ERI_Gy3z_Py_S_S_bb;
  Double I_ERI_G4z_D2y_S_S_bb = I_ERI_Hy4z_Py_S_S_bb+ABY*I_ERI_G4z_Py_S_S_bb;
  Double I_ERI_G4x_D2z_S_S_bb = I_ERI_H4xz_Pz_S_S_bb+ABZ*I_ERI_G4x_Pz_S_S_bb;
  Double I_ERI_G3xy_D2z_S_S_bb = I_ERI_H3xyz_Pz_S_S_bb+ABZ*I_ERI_G3xy_Pz_S_S_bb;
  Double I_ERI_G3xz_D2z_S_S_bb = I_ERI_H3x2z_Pz_S_S_bb+ABZ*I_ERI_G3xz_Pz_S_S_bb;
  Double I_ERI_G2x2y_D2z_S_S_bb = I_ERI_H2x2yz_Pz_S_S_bb+ABZ*I_ERI_G2x2y_Pz_S_S_bb;
  Double I_ERI_G2xyz_D2z_S_S_bb = I_ERI_H2xy2z_Pz_S_S_bb+ABZ*I_ERI_G2xyz_Pz_S_S_bb;
  Double I_ERI_G2x2z_D2z_S_S_bb = I_ERI_H2x3z_Pz_S_S_bb+ABZ*I_ERI_G2x2z_Pz_S_S_bb;
  Double I_ERI_Gx3y_D2z_S_S_bb = I_ERI_Hx3yz_Pz_S_S_bb+ABZ*I_ERI_Gx3y_Pz_S_S_bb;
  Double I_ERI_Gx2yz_D2z_S_S_bb = I_ERI_Hx2y2z_Pz_S_S_bb+ABZ*I_ERI_Gx2yz_Pz_S_S_bb;
  Double I_ERI_Gxy2z_D2z_S_S_bb = I_ERI_Hxy3z_Pz_S_S_bb+ABZ*I_ERI_Gxy2z_Pz_S_S_bb;
  Double I_ERI_Gx3z_D2z_S_S_bb = I_ERI_Hx4z_Pz_S_S_bb+ABZ*I_ERI_Gx3z_Pz_S_S_bb;
  Double I_ERI_G4y_D2z_S_S_bb = I_ERI_H4yz_Pz_S_S_bb+ABZ*I_ERI_G4y_Pz_S_S_bb;
  Double I_ERI_G3yz_D2z_S_S_bb = I_ERI_H3y2z_Pz_S_S_bb+ABZ*I_ERI_G3yz_Pz_S_S_bb;
  Double I_ERI_G2y2z_D2z_S_S_bb = I_ERI_H2y3z_Pz_S_S_bb+ABZ*I_ERI_G2y2z_Pz_S_S_bb;
  Double I_ERI_Gy3z_D2z_S_S_bb = I_ERI_Hy4z_Pz_S_S_bb+ABZ*I_ERI_Gy3z_Pz_S_S_bb;
  Double I_ERI_G4z_D2z_S_S_bb = I_ERI_H5z_Pz_S_S_bb+ABZ*I_ERI_G4z_Pz_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_F_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_D_S_S_bb
   ************************************************************/
  Double I_ERI_F3x_F3x_S_S_bb = I_ERI_G4x_D2x_S_S_bb+ABX*I_ERI_F3x_D2x_S_S_bb;
  Double I_ERI_F2xy_F3x_S_S_bb = I_ERI_G3xy_D2x_S_S_bb+ABX*I_ERI_F2xy_D2x_S_S_bb;
  Double I_ERI_F2xz_F3x_S_S_bb = I_ERI_G3xz_D2x_S_S_bb+ABX*I_ERI_F2xz_D2x_S_S_bb;
  Double I_ERI_Fx2y_F3x_S_S_bb = I_ERI_G2x2y_D2x_S_S_bb+ABX*I_ERI_Fx2y_D2x_S_S_bb;
  Double I_ERI_Fxyz_F3x_S_S_bb = I_ERI_G2xyz_D2x_S_S_bb+ABX*I_ERI_Fxyz_D2x_S_S_bb;
  Double I_ERI_Fx2z_F3x_S_S_bb = I_ERI_G2x2z_D2x_S_S_bb+ABX*I_ERI_Fx2z_D2x_S_S_bb;
  Double I_ERI_F3y_F3x_S_S_bb = I_ERI_Gx3y_D2x_S_S_bb+ABX*I_ERI_F3y_D2x_S_S_bb;
  Double I_ERI_F2yz_F3x_S_S_bb = I_ERI_Gx2yz_D2x_S_S_bb+ABX*I_ERI_F2yz_D2x_S_S_bb;
  Double I_ERI_Fy2z_F3x_S_S_bb = I_ERI_Gxy2z_D2x_S_S_bb+ABX*I_ERI_Fy2z_D2x_S_S_bb;
  Double I_ERI_F3z_F3x_S_S_bb = I_ERI_Gx3z_D2x_S_S_bb+ABX*I_ERI_F3z_D2x_S_S_bb;
  Double I_ERI_F3x_F2xy_S_S_bb = I_ERI_G3xy_D2x_S_S_bb+ABY*I_ERI_F3x_D2x_S_S_bb;
  Double I_ERI_F2xy_F2xy_S_S_bb = I_ERI_G2x2y_D2x_S_S_bb+ABY*I_ERI_F2xy_D2x_S_S_bb;
  Double I_ERI_F2xz_F2xy_S_S_bb = I_ERI_G2xyz_D2x_S_S_bb+ABY*I_ERI_F2xz_D2x_S_S_bb;
  Double I_ERI_Fx2y_F2xy_S_S_bb = I_ERI_Gx3y_D2x_S_S_bb+ABY*I_ERI_Fx2y_D2x_S_S_bb;
  Double I_ERI_Fxyz_F2xy_S_S_bb = I_ERI_Gx2yz_D2x_S_S_bb+ABY*I_ERI_Fxyz_D2x_S_S_bb;
  Double I_ERI_Fx2z_F2xy_S_S_bb = I_ERI_Gxy2z_D2x_S_S_bb+ABY*I_ERI_Fx2z_D2x_S_S_bb;
  Double I_ERI_F3y_F2xy_S_S_bb = I_ERI_G4y_D2x_S_S_bb+ABY*I_ERI_F3y_D2x_S_S_bb;
  Double I_ERI_F2yz_F2xy_S_S_bb = I_ERI_G3yz_D2x_S_S_bb+ABY*I_ERI_F2yz_D2x_S_S_bb;
  Double I_ERI_Fy2z_F2xy_S_S_bb = I_ERI_G2y2z_D2x_S_S_bb+ABY*I_ERI_Fy2z_D2x_S_S_bb;
  Double I_ERI_F3z_F2xy_S_S_bb = I_ERI_Gy3z_D2x_S_S_bb+ABY*I_ERI_F3z_D2x_S_S_bb;
  Double I_ERI_F3x_F2xz_S_S_bb = I_ERI_G3xz_D2x_S_S_bb+ABZ*I_ERI_F3x_D2x_S_S_bb;
  Double I_ERI_F2xy_F2xz_S_S_bb = I_ERI_G2xyz_D2x_S_S_bb+ABZ*I_ERI_F2xy_D2x_S_S_bb;
  Double I_ERI_F2xz_F2xz_S_S_bb = I_ERI_G2x2z_D2x_S_S_bb+ABZ*I_ERI_F2xz_D2x_S_S_bb;
  Double I_ERI_Fx2y_F2xz_S_S_bb = I_ERI_Gx2yz_D2x_S_S_bb+ABZ*I_ERI_Fx2y_D2x_S_S_bb;
  Double I_ERI_Fxyz_F2xz_S_S_bb = I_ERI_Gxy2z_D2x_S_S_bb+ABZ*I_ERI_Fxyz_D2x_S_S_bb;
  Double I_ERI_Fx2z_F2xz_S_S_bb = I_ERI_Gx3z_D2x_S_S_bb+ABZ*I_ERI_Fx2z_D2x_S_S_bb;
  Double I_ERI_F3y_F2xz_S_S_bb = I_ERI_G3yz_D2x_S_S_bb+ABZ*I_ERI_F3y_D2x_S_S_bb;
  Double I_ERI_F2yz_F2xz_S_S_bb = I_ERI_G2y2z_D2x_S_S_bb+ABZ*I_ERI_F2yz_D2x_S_S_bb;
  Double I_ERI_Fy2z_F2xz_S_S_bb = I_ERI_Gy3z_D2x_S_S_bb+ABZ*I_ERI_Fy2z_D2x_S_S_bb;
  Double I_ERI_F3z_F2xz_S_S_bb = I_ERI_G4z_D2x_S_S_bb+ABZ*I_ERI_F3z_D2x_S_S_bb;
  Double I_ERI_F3x_Fx2y_S_S_bb = I_ERI_G4x_D2y_S_S_bb+ABX*I_ERI_F3x_D2y_S_S_bb;
  Double I_ERI_F2xy_Fx2y_S_S_bb = I_ERI_G3xy_D2y_S_S_bb+ABX*I_ERI_F2xy_D2y_S_S_bb;
  Double I_ERI_F2xz_Fx2y_S_S_bb = I_ERI_G3xz_D2y_S_S_bb+ABX*I_ERI_F2xz_D2y_S_S_bb;
  Double I_ERI_Fx2y_Fx2y_S_S_bb = I_ERI_G2x2y_D2y_S_S_bb+ABX*I_ERI_Fx2y_D2y_S_S_bb;
  Double I_ERI_Fxyz_Fx2y_S_S_bb = I_ERI_G2xyz_D2y_S_S_bb+ABX*I_ERI_Fxyz_D2y_S_S_bb;
  Double I_ERI_Fx2z_Fx2y_S_S_bb = I_ERI_G2x2z_D2y_S_S_bb+ABX*I_ERI_Fx2z_D2y_S_S_bb;
  Double I_ERI_F3y_Fx2y_S_S_bb = I_ERI_Gx3y_D2y_S_S_bb+ABX*I_ERI_F3y_D2y_S_S_bb;
  Double I_ERI_F2yz_Fx2y_S_S_bb = I_ERI_Gx2yz_D2y_S_S_bb+ABX*I_ERI_F2yz_D2y_S_S_bb;
  Double I_ERI_Fy2z_Fx2y_S_S_bb = I_ERI_Gxy2z_D2y_S_S_bb+ABX*I_ERI_Fy2z_D2y_S_S_bb;
  Double I_ERI_F3z_Fx2y_S_S_bb = I_ERI_Gx3z_D2y_S_S_bb+ABX*I_ERI_F3z_D2y_S_S_bb;
  Double I_ERI_F3x_Fx2z_S_S_bb = I_ERI_G4x_D2z_S_S_bb+ABX*I_ERI_F3x_D2z_S_S_bb;
  Double I_ERI_F2xy_Fx2z_S_S_bb = I_ERI_G3xy_D2z_S_S_bb+ABX*I_ERI_F2xy_D2z_S_S_bb;
  Double I_ERI_F2xz_Fx2z_S_S_bb = I_ERI_G3xz_D2z_S_S_bb+ABX*I_ERI_F2xz_D2z_S_S_bb;
  Double I_ERI_Fx2y_Fx2z_S_S_bb = I_ERI_G2x2y_D2z_S_S_bb+ABX*I_ERI_Fx2y_D2z_S_S_bb;
  Double I_ERI_Fxyz_Fx2z_S_S_bb = I_ERI_G2xyz_D2z_S_S_bb+ABX*I_ERI_Fxyz_D2z_S_S_bb;
  Double I_ERI_Fx2z_Fx2z_S_S_bb = I_ERI_G2x2z_D2z_S_S_bb+ABX*I_ERI_Fx2z_D2z_S_S_bb;
  Double I_ERI_F3y_Fx2z_S_S_bb = I_ERI_Gx3y_D2z_S_S_bb+ABX*I_ERI_F3y_D2z_S_S_bb;
  Double I_ERI_F2yz_Fx2z_S_S_bb = I_ERI_Gx2yz_D2z_S_S_bb+ABX*I_ERI_F2yz_D2z_S_S_bb;
  Double I_ERI_Fy2z_Fx2z_S_S_bb = I_ERI_Gxy2z_D2z_S_S_bb+ABX*I_ERI_Fy2z_D2z_S_S_bb;
  Double I_ERI_F3z_Fx2z_S_S_bb = I_ERI_Gx3z_D2z_S_S_bb+ABX*I_ERI_F3z_D2z_S_S_bb;
  Double I_ERI_F3x_F3y_S_S_bb = I_ERI_G3xy_D2y_S_S_bb+ABY*I_ERI_F3x_D2y_S_S_bb;
  Double I_ERI_F2xy_F3y_S_S_bb = I_ERI_G2x2y_D2y_S_S_bb+ABY*I_ERI_F2xy_D2y_S_S_bb;
  Double I_ERI_F2xz_F3y_S_S_bb = I_ERI_G2xyz_D2y_S_S_bb+ABY*I_ERI_F2xz_D2y_S_S_bb;
  Double I_ERI_Fx2y_F3y_S_S_bb = I_ERI_Gx3y_D2y_S_S_bb+ABY*I_ERI_Fx2y_D2y_S_S_bb;
  Double I_ERI_Fxyz_F3y_S_S_bb = I_ERI_Gx2yz_D2y_S_S_bb+ABY*I_ERI_Fxyz_D2y_S_S_bb;
  Double I_ERI_Fx2z_F3y_S_S_bb = I_ERI_Gxy2z_D2y_S_S_bb+ABY*I_ERI_Fx2z_D2y_S_S_bb;
  Double I_ERI_F3y_F3y_S_S_bb = I_ERI_G4y_D2y_S_S_bb+ABY*I_ERI_F3y_D2y_S_S_bb;
  Double I_ERI_F2yz_F3y_S_S_bb = I_ERI_G3yz_D2y_S_S_bb+ABY*I_ERI_F2yz_D2y_S_S_bb;
  Double I_ERI_Fy2z_F3y_S_S_bb = I_ERI_G2y2z_D2y_S_S_bb+ABY*I_ERI_Fy2z_D2y_S_S_bb;
  Double I_ERI_F3z_F3y_S_S_bb = I_ERI_Gy3z_D2y_S_S_bb+ABY*I_ERI_F3z_D2y_S_S_bb;
  Double I_ERI_F3x_F2yz_S_S_bb = I_ERI_G3xz_D2y_S_S_bb+ABZ*I_ERI_F3x_D2y_S_S_bb;
  Double I_ERI_F2xy_F2yz_S_S_bb = I_ERI_G2xyz_D2y_S_S_bb+ABZ*I_ERI_F2xy_D2y_S_S_bb;
  Double I_ERI_F2xz_F2yz_S_S_bb = I_ERI_G2x2z_D2y_S_S_bb+ABZ*I_ERI_F2xz_D2y_S_S_bb;
  Double I_ERI_Fx2y_F2yz_S_S_bb = I_ERI_Gx2yz_D2y_S_S_bb+ABZ*I_ERI_Fx2y_D2y_S_S_bb;
  Double I_ERI_Fxyz_F2yz_S_S_bb = I_ERI_Gxy2z_D2y_S_S_bb+ABZ*I_ERI_Fxyz_D2y_S_S_bb;
  Double I_ERI_Fx2z_F2yz_S_S_bb = I_ERI_Gx3z_D2y_S_S_bb+ABZ*I_ERI_Fx2z_D2y_S_S_bb;
  Double I_ERI_F3y_F2yz_S_S_bb = I_ERI_G3yz_D2y_S_S_bb+ABZ*I_ERI_F3y_D2y_S_S_bb;
  Double I_ERI_F2yz_F2yz_S_S_bb = I_ERI_G2y2z_D2y_S_S_bb+ABZ*I_ERI_F2yz_D2y_S_S_bb;
  Double I_ERI_Fy2z_F2yz_S_S_bb = I_ERI_Gy3z_D2y_S_S_bb+ABZ*I_ERI_Fy2z_D2y_S_S_bb;
  Double I_ERI_F3z_F2yz_S_S_bb = I_ERI_G4z_D2y_S_S_bb+ABZ*I_ERI_F3z_D2y_S_S_bb;
  Double I_ERI_F3x_F3z_S_S_bb = I_ERI_G3xz_D2z_S_S_bb+ABZ*I_ERI_F3x_D2z_S_S_bb;
  Double I_ERI_F2xy_F3z_S_S_bb = I_ERI_G2xyz_D2z_S_S_bb+ABZ*I_ERI_F2xy_D2z_S_S_bb;
  Double I_ERI_F2xz_F3z_S_S_bb = I_ERI_G2x2z_D2z_S_S_bb+ABZ*I_ERI_F2xz_D2z_S_S_bb;
  Double I_ERI_Fx2y_F3z_S_S_bb = I_ERI_Gx2yz_D2z_S_S_bb+ABZ*I_ERI_Fx2y_D2z_S_S_bb;
  Double I_ERI_Fxyz_F3z_S_S_bb = I_ERI_Gxy2z_D2z_S_S_bb+ABZ*I_ERI_Fxyz_D2z_S_S_bb;
  Double I_ERI_Fx2z_F3z_S_S_bb = I_ERI_Gx3z_D2z_S_S_bb+ABZ*I_ERI_Fx2z_D2z_S_S_bb;
  Double I_ERI_F3y_F3z_S_S_bb = I_ERI_G3yz_D2z_S_S_bb+ABZ*I_ERI_F3y_D2z_S_S_bb;
  Double I_ERI_F2yz_F3z_S_S_bb = I_ERI_G2y2z_D2z_S_S_bb+ABZ*I_ERI_F2yz_D2z_S_S_bb;
  Double I_ERI_Fy2z_F3z_S_S_bb = I_ERI_Gy3z_D2z_S_S_bb+ABZ*I_ERI_Fy2z_D2z_S_S_bb;
  Double I_ERI_F3z_F3z_S_S_bb = I_ERI_G4z_D2z_S_S_bb+ABZ*I_ERI_F3z_D2z_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_I_P_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 24 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_K_S_S_S_bb
   * RHS shell quartet name: SQ_ERI_I_S_S_S_bb
   ************************************************************/
  Double I_ERI_I6x_Px_S_S_bb = I_ERI_K7x_S_S_S_bb+ABX*I_ERI_I6x_S_S_S_bb;
  Double I_ERI_I5xy_Px_S_S_bb = I_ERI_K6xy_S_S_S_bb+ABX*I_ERI_I5xy_S_S_S_bb;
  Double I_ERI_I5xz_Px_S_S_bb = I_ERI_K6xz_S_S_S_bb+ABX*I_ERI_I5xz_S_S_S_bb;
  Double I_ERI_I4x2y_Px_S_S_bb = I_ERI_K5x2y_S_S_S_bb+ABX*I_ERI_I4x2y_S_S_S_bb;
  Double I_ERI_I4xyz_Px_S_S_bb = I_ERI_K5xyz_S_S_S_bb+ABX*I_ERI_I4xyz_S_S_S_bb;
  Double I_ERI_I4x2z_Px_S_S_bb = I_ERI_K5x2z_S_S_S_bb+ABX*I_ERI_I4x2z_S_S_S_bb;
  Double I_ERI_I3x3y_Px_S_S_bb = I_ERI_K4x3y_S_S_S_bb+ABX*I_ERI_I3x3y_S_S_S_bb;
  Double I_ERI_I3x2yz_Px_S_S_bb = I_ERI_K4x2yz_S_S_S_bb+ABX*I_ERI_I3x2yz_S_S_S_bb;
  Double I_ERI_I3xy2z_Px_S_S_bb = I_ERI_K4xy2z_S_S_S_bb+ABX*I_ERI_I3xy2z_S_S_S_bb;
  Double I_ERI_I3x3z_Px_S_S_bb = I_ERI_K4x3z_S_S_S_bb+ABX*I_ERI_I3x3z_S_S_S_bb;
  Double I_ERI_I2x4y_Px_S_S_bb = I_ERI_K3x4y_S_S_S_bb+ABX*I_ERI_I2x4y_S_S_S_bb;
  Double I_ERI_I2x3yz_Px_S_S_bb = I_ERI_K3x3yz_S_S_S_bb+ABX*I_ERI_I2x3yz_S_S_S_bb;
  Double I_ERI_I2x2y2z_Px_S_S_bb = I_ERI_K3x2y2z_S_S_S_bb+ABX*I_ERI_I2x2y2z_S_S_S_bb;
  Double I_ERI_I2xy3z_Px_S_S_bb = I_ERI_K3xy3z_S_S_S_bb+ABX*I_ERI_I2xy3z_S_S_S_bb;
  Double I_ERI_I2x4z_Px_S_S_bb = I_ERI_K3x4z_S_S_S_bb+ABX*I_ERI_I2x4z_S_S_S_bb;
  Double I_ERI_Ix5y_Px_S_S_bb = I_ERI_K2x5y_S_S_S_bb+ABX*I_ERI_Ix5y_S_S_S_bb;
  Double I_ERI_Ix4yz_Px_S_S_bb = I_ERI_K2x4yz_S_S_S_bb+ABX*I_ERI_Ix4yz_S_S_S_bb;
  Double I_ERI_Ix3y2z_Px_S_S_bb = I_ERI_K2x3y2z_S_S_S_bb+ABX*I_ERI_Ix3y2z_S_S_S_bb;
  Double I_ERI_Ix2y3z_Px_S_S_bb = I_ERI_K2x2y3z_S_S_S_bb+ABX*I_ERI_Ix2y3z_S_S_S_bb;
  Double I_ERI_Ixy4z_Px_S_S_bb = I_ERI_K2xy4z_S_S_S_bb+ABX*I_ERI_Ixy4z_S_S_S_bb;
  Double I_ERI_Ix5z_Px_S_S_bb = I_ERI_K2x5z_S_S_S_bb+ABX*I_ERI_Ix5z_S_S_S_bb;
  Double I_ERI_I4x2y_Py_S_S_bb = I_ERI_K4x3y_S_S_S_bb+ABY*I_ERI_I4x2y_S_S_S_bb;
  Double I_ERI_I4xyz_Py_S_S_bb = I_ERI_K4x2yz_S_S_S_bb+ABY*I_ERI_I4xyz_S_S_S_bb;
  Double I_ERI_I3x3y_Py_S_S_bb = I_ERI_K3x4y_S_S_S_bb+ABY*I_ERI_I3x3y_S_S_S_bb;
  Double I_ERI_I3x2yz_Py_S_S_bb = I_ERI_K3x3yz_S_S_S_bb+ABY*I_ERI_I3x2yz_S_S_S_bb;
  Double I_ERI_I3xy2z_Py_S_S_bb = I_ERI_K3x2y2z_S_S_S_bb+ABY*I_ERI_I3xy2z_S_S_S_bb;
  Double I_ERI_I2x4y_Py_S_S_bb = I_ERI_K2x5y_S_S_S_bb+ABY*I_ERI_I2x4y_S_S_S_bb;
  Double I_ERI_I2x3yz_Py_S_S_bb = I_ERI_K2x4yz_S_S_S_bb+ABY*I_ERI_I2x3yz_S_S_S_bb;
  Double I_ERI_I2x2y2z_Py_S_S_bb = I_ERI_K2x3y2z_S_S_S_bb+ABY*I_ERI_I2x2y2z_S_S_S_bb;
  Double I_ERI_I2xy3z_Py_S_S_bb = I_ERI_K2x2y3z_S_S_S_bb+ABY*I_ERI_I2xy3z_S_S_S_bb;
  Double I_ERI_Ix5y_Py_S_S_bb = I_ERI_Kx6y_S_S_S_bb+ABY*I_ERI_Ix5y_S_S_S_bb;
  Double I_ERI_Ix4yz_Py_S_S_bb = I_ERI_Kx5yz_S_S_S_bb+ABY*I_ERI_Ix4yz_S_S_S_bb;
  Double I_ERI_Ix3y2z_Py_S_S_bb = I_ERI_Kx4y2z_S_S_S_bb+ABY*I_ERI_Ix3y2z_S_S_S_bb;
  Double I_ERI_Ix2y3z_Py_S_S_bb = I_ERI_Kx3y3z_S_S_S_bb+ABY*I_ERI_Ix2y3z_S_S_S_bb;
  Double I_ERI_Ixy4z_Py_S_S_bb = I_ERI_Kx2y4z_S_S_S_bb+ABY*I_ERI_Ixy4z_S_S_S_bb;
  Double I_ERI_I6y_Py_S_S_bb = I_ERI_K7y_S_S_S_bb+ABY*I_ERI_I6y_S_S_S_bb;
  Double I_ERI_I5yz_Py_S_S_bb = I_ERI_K6yz_S_S_S_bb+ABY*I_ERI_I5yz_S_S_S_bb;
  Double I_ERI_I4y2z_Py_S_S_bb = I_ERI_K5y2z_S_S_S_bb+ABY*I_ERI_I4y2z_S_S_S_bb;
  Double I_ERI_I3y3z_Py_S_S_bb = I_ERI_K4y3z_S_S_S_bb+ABY*I_ERI_I3y3z_S_S_S_bb;
  Double I_ERI_I2y4z_Py_S_S_bb = I_ERI_K3y4z_S_S_S_bb+ABY*I_ERI_I2y4z_S_S_S_bb;
  Double I_ERI_Iy5z_Py_S_S_bb = I_ERI_K2y5z_S_S_S_bb+ABY*I_ERI_Iy5z_S_S_S_bb;
  Double I_ERI_I4xyz_Pz_S_S_bb = I_ERI_K4xy2z_S_S_S_bb+ABZ*I_ERI_I4xyz_S_S_S_bb;
  Double I_ERI_I4x2z_Pz_S_S_bb = I_ERI_K4x3z_S_S_S_bb+ABZ*I_ERI_I4x2z_S_S_S_bb;
  Double I_ERI_I3x2yz_Pz_S_S_bb = I_ERI_K3x2y2z_S_S_S_bb+ABZ*I_ERI_I3x2yz_S_S_S_bb;
  Double I_ERI_I3xy2z_Pz_S_S_bb = I_ERI_K3xy3z_S_S_S_bb+ABZ*I_ERI_I3xy2z_S_S_S_bb;
  Double I_ERI_I3x3z_Pz_S_S_bb = I_ERI_K3x4z_S_S_S_bb+ABZ*I_ERI_I3x3z_S_S_S_bb;
  Double I_ERI_I2x3yz_Pz_S_S_bb = I_ERI_K2x3y2z_S_S_S_bb+ABZ*I_ERI_I2x3yz_S_S_S_bb;
  Double I_ERI_I2x2y2z_Pz_S_S_bb = I_ERI_K2x2y3z_S_S_S_bb+ABZ*I_ERI_I2x2y2z_S_S_S_bb;
  Double I_ERI_I2xy3z_Pz_S_S_bb = I_ERI_K2xy4z_S_S_S_bb+ABZ*I_ERI_I2xy3z_S_S_S_bb;
  Double I_ERI_I2x4z_Pz_S_S_bb = I_ERI_K2x5z_S_S_S_bb+ABZ*I_ERI_I2x4z_S_S_S_bb;
  Double I_ERI_Ix4yz_Pz_S_S_bb = I_ERI_Kx4y2z_S_S_S_bb+ABZ*I_ERI_Ix4yz_S_S_S_bb;
  Double I_ERI_Ix3y2z_Pz_S_S_bb = I_ERI_Kx3y3z_S_S_S_bb+ABZ*I_ERI_Ix3y2z_S_S_S_bb;
  Double I_ERI_Ix2y3z_Pz_S_S_bb = I_ERI_Kx2y4z_S_S_S_bb+ABZ*I_ERI_Ix2y3z_S_S_S_bb;
  Double I_ERI_Ixy4z_Pz_S_S_bb = I_ERI_Kxy5z_S_S_S_bb+ABZ*I_ERI_Ixy4z_S_S_S_bb;
  Double I_ERI_Ix5z_Pz_S_S_bb = I_ERI_Kx6z_S_S_S_bb+ABZ*I_ERI_Ix5z_S_S_S_bb;
  Double I_ERI_I4y2z_Pz_S_S_bb = I_ERI_K4y3z_S_S_S_bb+ABZ*I_ERI_I4y2z_S_S_S_bb;
  Double I_ERI_I3y3z_Pz_S_S_bb = I_ERI_K3y4z_S_S_S_bb+ABZ*I_ERI_I3y3z_S_S_S_bb;
  Double I_ERI_I2y4z_Pz_S_S_bb = I_ERI_K2y5z_S_S_S_bb+ABZ*I_ERI_I2y4z_S_S_S_bb;
  Double I_ERI_Iy5z_Pz_S_S_bb = I_ERI_Ky6z_S_S_S_bb+ABZ*I_ERI_Iy5z_S_S_S_bb;
  Double I_ERI_I6z_Pz_S_S_bb = I_ERI_K7z_S_S_S_bb+ABZ*I_ERI_I6z_S_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_H_D_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 66 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_P_S_S_bb
   * RHS shell quartet name: SQ_ERI_H_P_S_S_bb
   ************************************************************/
  Double I_ERI_H5x_D2x_S_S_bb = I_ERI_I6x_Px_S_S_bb+ABX*I_ERI_H5x_Px_S_S_bb;
  Double I_ERI_H4xy_D2x_S_S_bb = I_ERI_I5xy_Px_S_S_bb+ABX*I_ERI_H4xy_Px_S_S_bb;
  Double I_ERI_H4xz_D2x_S_S_bb = I_ERI_I5xz_Px_S_S_bb+ABX*I_ERI_H4xz_Px_S_S_bb;
  Double I_ERI_H3x2y_D2x_S_S_bb = I_ERI_I4x2y_Px_S_S_bb+ABX*I_ERI_H3x2y_Px_S_S_bb;
  Double I_ERI_H3xyz_D2x_S_S_bb = I_ERI_I4xyz_Px_S_S_bb+ABX*I_ERI_H3xyz_Px_S_S_bb;
  Double I_ERI_H3x2z_D2x_S_S_bb = I_ERI_I4x2z_Px_S_S_bb+ABX*I_ERI_H3x2z_Px_S_S_bb;
  Double I_ERI_H2x3y_D2x_S_S_bb = I_ERI_I3x3y_Px_S_S_bb+ABX*I_ERI_H2x3y_Px_S_S_bb;
  Double I_ERI_H2x2yz_D2x_S_S_bb = I_ERI_I3x2yz_Px_S_S_bb+ABX*I_ERI_H2x2yz_Px_S_S_bb;
  Double I_ERI_H2xy2z_D2x_S_S_bb = I_ERI_I3xy2z_Px_S_S_bb+ABX*I_ERI_H2xy2z_Px_S_S_bb;
  Double I_ERI_H2x3z_D2x_S_S_bb = I_ERI_I3x3z_Px_S_S_bb+ABX*I_ERI_H2x3z_Px_S_S_bb;
  Double I_ERI_Hx4y_D2x_S_S_bb = I_ERI_I2x4y_Px_S_S_bb+ABX*I_ERI_Hx4y_Px_S_S_bb;
  Double I_ERI_Hx3yz_D2x_S_S_bb = I_ERI_I2x3yz_Px_S_S_bb+ABX*I_ERI_Hx3yz_Px_S_S_bb;
  Double I_ERI_Hx2y2z_D2x_S_S_bb = I_ERI_I2x2y2z_Px_S_S_bb+ABX*I_ERI_Hx2y2z_Px_S_S_bb;
  Double I_ERI_Hxy3z_D2x_S_S_bb = I_ERI_I2xy3z_Px_S_S_bb+ABX*I_ERI_Hxy3z_Px_S_S_bb;
  Double I_ERI_Hx4z_D2x_S_S_bb = I_ERI_I2x4z_Px_S_S_bb+ABX*I_ERI_Hx4z_Px_S_S_bb;
  Double I_ERI_H5y_D2x_S_S_bb = I_ERI_Ix5y_Px_S_S_bb+ABX*I_ERI_H5y_Px_S_S_bb;
  Double I_ERI_H4yz_D2x_S_S_bb = I_ERI_Ix4yz_Px_S_S_bb+ABX*I_ERI_H4yz_Px_S_S_bb;
  Double I_ERI_H3y2z_D2x_S_S_bb = I_ERI_Ix3y2z_Px_S_S_bb+ABX*I_ERI_H3y2z_Px_S_S_bb;
  Double I_ERI_H2y3z_D2x_S_S_bb = I_ERI_Ix2y3z_Px_S_S_bb+ABX*I_ERI_H2y3z_Px_S_S_bb;
  Double I_ERI_Hy4z_D2x_S_S_bb = I_ERI_Ixy4z_Px_S_S_bb+ABX*I_ERI_Hy4z_Px_S_S_bb;
  Double I_ERI_H5z_D2x_S_S_bb = I_ERI_Ix5z_Px_S_S_bb+ABX*I_ERI_H5z_Px_S_S_bb;
  Double I_ERI_H4xy_D2y_S_S_bb = I_ERI_I4x2y_Py_S_S_bb+ABY*I_ERI_H4xy_Py_S_S_bb;
  Double I_ERI_H4xz_D2y_S_S_bb = I_ERI_I4xyz_Py_S_S_bb+ABY*I_ERI_H4xz_Py_S_S_bb;
  Double I_ERI_H3x2y_D2y_S_S_bb = I_ERI_I3x3y_Py_S_S_bb+ABY*I_ERI_H3x2y_Py_S_S_bb;
  Double I_ERI_H3xyz_D2y_S_S_bb = I_ERI_I3x2yz_Py_S_S_bb+ABY*I_ERI_H3xyz_Py_S_S_bb;
  Double I_ERI_H3x2z_D2y_S_S_bb = I_ERI_I3xy2z_Py_S_S_bb+ABY*I_ERI_H3x2z_Py_S_S_bb;
  Double I_ERI_H2x3y_D2y_S_S_bb = I_ERI_I2x4y_Py_S_S_bb+ABY*I_ERI_H2x3y_Py_S_S_bb;
  Double I_ERI_H2x2yz_D2y_S_S_bb = I_ERI_I2x3yz_Py_S_S_bb+ABY*I_ERI_H2x2yz_Py_S_S_bb;
  Double I_ERI_H2xy2z_D2y_S_S_bb = I_ERI_I2x2y2z_Py_S_S_bb+ABY*I_ERI_H2xy2z_Py_S_S_bb;
  Double I_ERI_H2x3z_D2y_S_S_bb = I_ERI_I2xy3z_Py_S_S_bb+ABY*I_ERI_H2x3z_Py_S_S_bb;
  Double I_ERI_Hx4y_D2y_S_S_bb = I_ERI_Ix5y_Py_S_S_bb+ABY*I_ERI_Hx4y_Py_S_S_bb;
  Double I_ERI_Hx3yz_D2y_S_S_bb = I_ERI_Ix4yz_Py_S_S_bb+ABY*I_ERI_Hx3yz_Py_S_S_bb;
  Double I_ERI_Hx2y2z_D2y_S_S_bb = I_ERI_Ix3y2z_Py_S_S_bb+ABY*I_ERI_Hx2y2z_Py_S_S_bb;
  Double I_ERI_Hxy3z_D2y_S_S_bb = I_ERI_Ix2y3z_Py_S_S_bb+ABY*I_ERI_Hxy3z_Py_S_S_bb;
  Double I_ERI_Hx4z_D2y_S_S_bb = I_ERI_Ixy4z_Py_S_S_bb+ABY*I_ERI_Hx4z_Py_S_S_bb;
  Double I_ERI_H5y_D2y_S_S_bb = I_ERI_I6y_Py_S_S_bb+ABY*I_ERI_H5y_Py_S_S_bb;
  Double I_ERI_H4yz_D2y_S_S_bb = I_ERI_I5yz_Py_S_S_bb+ABY*I_ERI_H4yz_Py_S_S_bb;
  Double I_ERI_H3y2z_D2y_S_S_bb = I_ERI_I4y2z_Py_S_S_bb+ABY*I_ERI_H3y2z_Py_S_S_bb;
  Double I_ERI_H2y3z_D2y_S_S_bb = I_ERI_I3y3z_Py_S_S_bb+ABY*I_ERI_H2y3z_Py_S_S_bb;
  Double I_ERI_Hy4z_D2y_S_S_bb = I_ERI_I2y4z_Py_S_S_bb+ABY*I_ERI_Hy4z_Py_S_S_bb;
  Double I_ERI_H5z_D2y_S_S_bb = I_ERI_Iy5z_Py_S_S_bb+ABY*I_ERI_H5z_Py_S_S_bb;
  Double I_ERI_H4xy_D2z_S_S_bb = I_ERI_I4xyz_Pz_S_S_bb+ABZ*I_ERI_H4xy_Pz_S_S_bb;
  Double I_ERI_H4xz_D2z_S_S_bb = I_ERI_I4x2z_Pz_S_S_bb+ABZ*I_ERI_H4xz_Pz_S_S_bb;
  Double I_ERI_H3x2y_D2z_S_S_bb = I_ERI_I3x2yz_Pz_S_S_bb+ABZ*I_ERI_H3x2y_Pz_S_S_bb;
  Double I_ERI_H3xyz_D2z_S_S_bb = I_ERI_I3xy2z_Pz_S_S_bb+ABZ*I_ERI_H3xyz_Pz_S_S_bb;
  Double I_ERI_H3x2z_D2z_S_S_bb = I_ERI_I3x3z_Pz_S_S_bb+ABZ*I_ERI_H3x2z_Pz_S_S_bb;
  Double I_ERI_H2x3y_D2z_S_S_bb = I_ERI_I2x3yz_Pz_S_S_bb+ABZ*I_ERI_H2x3y_Pz_S_S_bb;
  Double I_ERI_H2x2yz_D2z_S_S_bb = I_ERI_I2x2y2z_Pz_S_S_bb+ABZ*I_ERI_H2x2yz_Pz_S_S_bb;
  Double I_ERI_H2xy2z_D2z_S_S_bb = I_ERI_I2xy3z_Pz_S_S_bb+ABZ*I_ERI_H2xy2z_Pz_S_S_bb;
  Double I_ERI_H2x3z_D2z_S_S_bb = I_ERI_I2x4z_Pz_S_S_bb+ABZ*I_ERI_H2x3z_Pz_S_S_bb;
  Double I_ERI_Hx4y_D2z_S_S_bb = I_ERI_Ix4yz_Pz_S_S_bb+ABZ*I_ERI_Hx4y_Pz_S_S_bb;
  Double I_ERI_Hx3yz_D2z_S_S_bb = I_ERI_Ix3y2z_Pz_S_S_bb+ABZ*I_ERI_Hx3yz_Pz_S_S_bb;
  Double I_ERI_Hx2y2z_D2z_S_S_bb = I_ERI_Ix2y3z_Pz_S_S_bb+ABZ*I_ERI_Hx2y2z_Pz_S_S_bb;
  Double I_ERI_Hxy3z_D2z_S_S_bb = I_ERI_Ixy4z_Pz_S_S_bb+ABZ*I_ERI_Hxy3z_Pz_S_S_bb;
  Double I_ERI_Hx4z_D2z_S_S_bb = I_ERI_Ix5z_Pz_S_S_bb+ABZ*I_ERI_Hx4z_Pz_S_S_bb;
  Double I_ERI_H4yz_D2z_S_S_bb = I_ERI_I4y2z_Pz_S_S_bb+ABZ*I_ERI_H4yz_Pz_S_S_bb;
  Double I_ERI_H3y2z_D2z_S_S_bb = I_ERI_I3y3z_Pz_S_S_bb+ABZ*I_ERI_H3y2z_Pz_S_S_bb;
  Double I_ERI_H2y3z_D2z_S_S_bb = I_ERI_I2y4z_Pz_S_S_bb+ABZ*I_ERI_H2y3z_Pz_S_S_bb;
  Double I_ERI_Hy4z_D2z_S_S_bb = I_ERI_Iy5z_Pz_S_S_bb+ABZ*I_ERI_Hy4z_Pz_S_S_bb;
  Double I_ERI_H5z_D2z_S_S_bb = I_ERI_I6z_Pz_S_S_bb+ABZ*I_ERI_H5z_Pz_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_G_F_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 51 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_S_S_bb
   * RHS shell quartet name: SQ_ERI_G_D_S_S_bb
   ************************************************************/
  Double I_ERI_G4x_F3x_S_S_bb = I_ERI_H5x_D2x_S_S_bb+ABX*I_ERI_G4x_D2x_S_S_bb;
  Double I_ERI_G3xy_F3x_S_S_bb = I_ERI_H4xy_D2x_S_S_bb+ABX*I_ERI_G3xy_D2x_S_S_bb;
  Double I_ERI_G3xz_F3x_S_S_bb = I_ERI_H4xz_D2x_S_S_bb+ABX*I_ERI_G3xz_D2x_S_S_bb;
  Double I_ERI_G2x2y_F3x_S_S_bb = I_ERI_H3x2y_D2x_S_S_bb+ABX*I_ERI_G2x2y_D2x_S_S_bb;
  Double I_ERI_G2xyz_F3x_S_S_bb = I_ERI_H3xyz_D2x_S_S_bb+ABX*I_ERI_G2xyz_D2x_S_S_bb;
  Double I_ERI_G2x2z_F3x_S_S_bb = I_ERI_H3x2z_D2x_S_S_bb+ABX*I_ERI_G2x2z_D2x_S_S_bb;
  Double I_ERI_Gx3y_F3x_S_S_bb = I_ERI_H2x3y_D2x_S_S_bb+ABX*I_ERI_Gx3y_D2x_S_S_bb;
  Double I_ERI_Gx2yz_F3x_S_S_bb = I_ERI_H2x2yz_D2x_S_S_bb+ABX*I_ERI_Gx2yz_D2x_S_S_bb;
  Double I_ERI_Gxy2z_F3x_S_S_bb = I_ERI_H2xy2z_D2x_S_S_bb+ABX*I_ERI_Gxy2z_D2x_S_S_bb;
  Double I_ERI_Gx3z_F3x_S_S_bb = I_ERI_H2x3z_D2x_S_S_bb+ABX*I_ERI_Gx3z_D2x_S_S_bb;
  Double I_ERI_G4y_F3x_S_S_bb = I_ERI_Hx4y_D2x_S_S_bb+ABX*I_ERI_G4y_D2x_S_S_bb;
  Double I_ERI_G3yz_F3x_S_S_bb = I_ERI_Hx3yz_D2x_S_S_bb+ABX*I_ERI_G3yz_D2x_S_S_bb;
  Double I_ERI_G2y2z_F3x_S_S_bb = I_ERI_Hx2y2z_D2x_S_S_bb+ABX*I_ERI_G2y2z_D2x_S_S_bb;
  Double I_ERI_Gy3z_F3x_S_S_bb = I_ERI_Hxy3z_D2x_S_S_bb+ABX*I_ERI_Gy3z_D2x_S_S_bb;
  Double I_ERI_G4z_F3x_S_S_bb = I_ERI_Hx4z_D2x_S_S_bb+ABX*I_ERI_G4z_D2x_S_S_bb;
  Double I_ERI_G3xy_F2xy_S_S_bb = I_ERI_H3x2y_D2x_S_S_bb+ABY*I_ERI_G3xy_D2x_S_S_bb;
  Double I_ERI_G3xz_F2xy_S_S_bb = I_ERI_H3xyz_D2x_S_S_bb+ABY*I_ERI_G3xz_D2x_S_S_bb;
  Double I_ERI_G2x2y_F2xy_S_S_bb = I_ERI_H2x3y_D2x_S_S_bb+ABY*I_ERI_G2x2y_D2x_S_S_bb;
  Double I_ERI_G2xyz_F2xy_S_S_bb = I_ERI_H2x2yz_D2x_S_S_bb+ABY*I_ERI_G2xyz_D2x_S_S_bb;
  Double I_ERI_G2x2z_F2xy_S_S_bb = I_ERI_H2xy2z_D2x_S_S_bb+ABY*I_ERI_G2x2z_D2x_S_S_bb;
  Double I_ERI_Gx3y_F2xy_S_S_bb = I_ERI_Hx4y_D2x_S_S_bb+ABY*I_ERI_Gx3y_D2x_S_S_bb;
  Double I_ERI_Gx2yz_F2xy_S_S_bb = I_ERI_Hx3yz_D2x_S_S_bb+ABY*I_ERI_Gx2yz_D2x_S_S_bb;
  Double I_ERI_Gxy2z_F2xy_S_S_bb = I_ERI_Hx2y2z_D2x_S_S_bb+ABY*I_ERI_Gxy2z_D2x_S_S_bb;
  Double I_ERI_Gx3z_F2xy_S_S_bb = I_ERI_Hxy3z_D2x_S_S_bb+ABY*I_ERI_Gx3z_D2x_S_S_bb;
  Double I_ERI_G4y_F2xy_S_S_bb = I_ERI_H5y_D2x_S_S_bb+ABY*I_ERI_G4y_D2x_S_S_bb;
  Double I_ERI_G3yz_F2xy_S_S_bb = I_ERI_H4yz_D2x_S_S_bb+ABY*I_ERI_G3yz_D2x_S_S_bb;
  Double I_ERI_G2y2z_F2xy_S_S_bb = I_ERI_H3y2z_D2x_S_S_bb+ABY*I_ERI_G2y2z_D2x_S_S_bb;
  Double I_ERI_Gy3z_F2xy_S_S_bb = I_ERI_H2y3z_D2x_S_S_bb+ABY*I_ERI_Gy3z_D2x_S_S_bb;
  Double I_ERI_G4z_F2xy_S_S_bb = I_ERI_Hy4z_D2x_S_S_bb+ABY*I_ERI_G4z_D2x_S_S_bb;
  Double I_ERI_G3xz_F2xz_S_S_bb = I_ERI_H3x2z_D2x_S_S_bb+ABZ*I_ERI_G3xz_D2x_S_S_bb;
  Double I_ERI_G2xyz_F2xz_S_S_bb = I_ERI_H2xy2z_D2x_S_S_bb+ABZ*I_ERI_G2xyz_D2x_S_S_bb;
  Double I_ERI_G2x2z_F2xz_S_S_bb = I_ERI_H2x3z_D2x_S_S_bb+ABZ*I_ERI_G2x2z_D2x_S_S_bb;
  Double I_ERI_Gx2yz_F2xz_S_S_bb = I_ERI_Hx2y2z_D2x_S_S_bb+ABZ*I_ERI_Gx2yz_D2x_S_S_bb;
  Double I_ERI_Gxy2z_F2xz_S_S_bb = I_ERI_Hxy3z_D2x_S_S_bb+ABZ*I_ERI_Gxy2z_D2x_S_S_bb;
  Double I_ERI_Gx3z_F2xz_S_S_bb = I_ERI_Hx4z_D2x_S_S_bb+ABZ*I_ERI_Gx3z_D2x_S_S_bb;
  Double I_ERI_G3yz_F2xz_S_S_bb = I_ERI_H3y2z_D2x_S_S_bb+ABZ*I_ERI_G3yz_D2x_S_S_bb;
  Double I_ERI_G2y2z_F2xz_S_S_bb = I_ERI_H2y3z_D2x_S_S_bb+ABZ*I_ERI_G2y2z_D2x_S_S_bb;
  Double I_ERI_Gy3z_F2xz_S_S_bb = I_ERI_Hy4z_D2x_S_S_bb+ABZ*I_ERI_Gy3z_D2x_S_S_bb;
  Double I_ERI_G4z_F2xz_S_S_bb = I_ERI_H5z_D2x_S_S_bb+ABZ*I_ERI_G4z_D2x_S_S_bb;
  Double I_ERI_G3xz_Fx2y_S_S_bb = I_ERI_H4xz_D2y_S_S_bb+ABX*I_ERI_G3xz_D2y_S_S_bb;
  Double I_ERI_G2xyz_Fx2y_S_S_bb = I_ERI_H3xyz_D2y_S_S_bb+ABX*I_ERI_G2xyz_D2y_S_S_bb;
  Double I_ERI_G2x2z_Fx2y_S_S_bb = I_ERI_H3x2z_D2y_S_S_bb+ABX*I_ERI_G2x2z_D2y_S_S_bb;
  Double I_ERI_Gx2yz_Fx2y_S_S_bb = I_ERI_H2x2yz_D2y_S_S_bb+ABX*I_ERI_Gx2yz_D2y_S_S_bb;
  Double I_ERI_Gxy2z_Fx2y_S_S_bb = I_ERI_H2xy2z_D2y_S_S_bb+ABX*I_ERI_Gxy2z_D2y_S_S_bb;
  Double I_ERI_Gx3z_Fx2y_S_S_bb = I_ERI_H2x3z_D2y_S_S_bb+ABX*I_ERI_Gx3z_D2y_S_S_bb;
  Double I_ERI_G3yz_Fx2y_S_S_bb = I_ERI_Hx3yz_D2y_S_S_bb+ABX*I_ERI_G3yz_D2y_S_S_bb;
  Double I_ERI_G2y2z_Fx2y_S_S_bb = I_ERI_Hx2y2z_D2y_S_S_bb+ABX*I_ERI_G2y2z_D2y_S_S_bb;
  Double I_ERI_Gy3z_Fx2y_S_S_bb = I_ERI_Hxy3z_D2y_S_S_bb+ABX*I_ERI_Gy3z_D2y_S_S_bb;
  Double I_ERI_G4z_Fx2y_S_S_bb = I_ERI_Hx4z_D2y_S_S_bb+ABX*I_ERI_G4z_D2y_S_S_bb;
  Double I_ERI_G3xy_Fx2z_S_S_bb = I_ERI_H4xy_D2z_S_S_bb+ABX*I_ERI_G3xy_D2z_S_S_bb;
  Double I_ERI_G2x2y_Fx2z_S_S_bb = I_ERI_H3x2y_D2z_S_S_bb+ABX*I_ERI_G2x2y_D2z_S_S_bb;
  Double I_ERI_G2xyz_Fx2z_S_S_bb = I_ERI_H3xyz_D2z_S_S_bb+ABX*I_ERI_G2xyz_D2z_S_S_bb;
  Double I_ERI_Gx3y_Fx2z_S_S_bb = I_ERI_H2x3y_D2z_S_S_bb+ABX*I_ERI_Gx3y_D2z_S_S_bb;
  Double I_ERI_Gx2yz_Fx2z_S_S_bb = I_ERI_H2x2yz_D2z_S_S_bb+ABX*I_ERI_Gx2yz_D2z_S_S_bb;
  Double I_ERI_Gxy2z_Fx2z_S_S_bb = I_ERI_H2xy2z_D2z_S_S_bb+ABX*I_ERI_Gxy2z_D2z_S_S_bb;
  Double I_ERI_G4y_Fx2z_S_S_bb = I_ERI_Hx4y_D2z_S_S_bb+ABX*I_ERI_G4y_D2z_S_S_bb;
  Double I_ERI_G3yz_Fx2z_S_S_bb = I_ERI_Hx3yz_D2z_S_S_bb+ABX*I_ERI_G3yz_D2z_S_S_bb;
  Double I_ERI_G2y2z_Fx2z_S_S_bb = I_ERI_Hx2y2z_D2z_S_S_bb+ABX*I_ERI_G2y2z_D2z_S_S_bb;
  Double I_ERI_Gy3z_Fx2z_S_S_bb = I_ERI_Hxy3z_D2z_S_S_bb+ABX*I_ERI_Gy3z_D2z_S_S_bb;
  Double I_ERI_G4x_F3y_S_S_bb = I_ERI_H4xy_D2y_S_S_bb+ABY*I_ERI_G4x_D2y_S_S_bb;
  Double I_ERI_G3xy_F3y_S_S_bb = I_ERI_H3x2y_D2y_S_S_bb+ABY*I_ERI_G3xy_D2y_S_S_bb;
  Double I_ERI_G3xz_F3y_S_S_bb = I_ERI_H3xyz_D2y_S_S_bb+ABY*I_ERI_G3xz_D2y_S_S_bb;
  Double I_ERI_G2x2y_F3y_S_S_bb = I_ERI_H2x3y_D2y_S_S_bb+ABY*I_ERI_G2x2y_D2y_S_S_bb;
  Double I_ERI_G2xyz_F3y_S_S_bb = I_ERI_H2x2yz_D2y_S_S_bb+ABY*I_ERI_G2xyz_D2y_S_S_bb;
  Double I_ERI_G2x2z_F3y_S_S_bb = I_ERI_H2xy2z_D2y_S_S_bb+ABY*I_ERI_G2x2z_D2y_S_S_bb;
  Double I_ERI_Gx3y_F3y_S_S_bb = I_ERI_Hx4y_D2y_S_S_bb+ABY*I_ERI_Gx3y_D2y_S_S_bb;
  Double I_ERI_Gx2yz_F3y_S_S_bb = I_ERI_Hx3yz_D2y_S_S_bb+ABY*I_ERI_Gx2yz_D2y_S_S_bb;
  Double I_ERI_Gxy2z_F3y_S_S_bb = I_ERI_Hx2y2z_D2y_S_S_bb+ABY*I_ERI_Gxy2z_D2y_S_S_bb;
  Double I_ERI_Gx3z_F3y_S_S_bb = I_ERI_Hxy3z_D2y_S_S_bb+ABY*I_ERI_Gx3z_D2y_S_S_bb;
  Double I_ERI_G4y_F3y_S_S_bb = I_ERI_H5y_D2y_S_S_bb+ABY*I_ERI_G4y_D2y_S_S_bb;
  Double I_ERI_G3yz_F3y_S_S_bb = I_ERI_H4yz_D2y_S_S_bb+ABY*I_ERI_G3yz_D2y_S_S_bb;
  Double I_ERI_G2y2z_F3y_S_S_bb = I_ERI_H3y2z_D2y_S_S_bb+ABY*I_ERI_G2y2z_D2y_S_S_bb;
  Double I_ERI_Gy3z_F3y_S_S_bb = I_ERI_H2y3z_D2y_S_S_bb+ABY*I_ERI_Gy3z_D2y_S_S_bb;
  Double I_ERI_G4z_F3y_S_S_bb = I_ERI_Hy4z_D2y_S_S_bb+ABY*I_ERI_G4z_D2y_S_S_bb;
  Double I_ERI_G3xz_F2yz_S_S_bb = I_ERI_H3x2z_D2y_S_S_bb+ABZ*I_ERI_G3xz_D2y_S_S_bb;
  Double I_ERI_G2xyz_F2yz_S_S_bb = I_ERI_H2xy2z_D2y_S_S_bb+ABZ*I_ERI_G2xyz_D2y_S_S_bb;
  Double I_ERI_G2x2z_F2yz_S_S_bb = I_ERI_H2x3z_D2y_S_S_bb+ABZ*I_ERI_G2x2z_D2y_S_S_bb;
  Double I_ERI_Gx2yz_F2yz_S_S_bb = I_ERI_Hx2y2z_D2y_S_S_bb+ABZ*I_ERI_Gx2yz_D2y_S_S_bb;
  Double I_ERI_Gxy2z_F2yz_S_S_bb = I_ERI_Hxy3z_D2y_S_S_bb+ABZ*I_ERI_Gxy2z_D2y_S_S_bb;
  Double I_ERI_Gx3z_F2yz_S_S_bb = I_ERI_Hx4z_D2y_S_S_bb+ABZ*I_ERI_Gx3z_D2y_S_S_bb;
  Double I_ERI_G3yz_F2yz_S_S_bb = I_ERI_H3y2z_D2y_S_S_bb+ABZ*I_ERI_G3yz_D2y_S_S_bb;
  Double I_ERI_G2y2z_F2yz_S_S_bb = I_ERI_H2y3z_D2y_S_S_bb+ABZ*I_ERI_G2y2z_D2y_S_S_bb;
  Double I_ERI_Gy3z_F2yz_S_S_bb = I_ERI_Hy4z_D2y_S_S_bb+ABZ*I_ERI_Gy3z_D2y_S_S_bb;
  Double I_ERI_G4z_F2yz_S_S_bb = I_ERI_H5z_D2y_S_S_bb+ABZ*I_ERI_G4z_D2y_S_S_bb;
  Double I_ERI_G4x_F3z_S_S_bb = I_ERI_H4xz_D2z_S_S_bb+ABZ*I_ERI_G4x_D2z_S_S_bb;
  Double I_ERI_G3xy_F3z_S_S_bb = I_ERI_H3xyz_D2z_S_S_bb+ABZ*I_ERI_G3xy_D2z_S_S_bb;
  Double I_ERI_G3xz_F3z_S_S_bb = I_ERI_H3x2z_D2z_S_S_bb+ABZ*I_ERI_G3xz_D2z_S_S_bb;
  Double I_ERI_G2x2y_F3z_S_S_bb = I_ERI_H2x2yz_D2z_S_S_bb+ABZ*I_ERI_G2x2y_D2z_S_S_bb;
  Double I_ERI_G2xyz_F3z_S_S_bb = I_ERI_H2xy2z_D2z_S_S_bb+ABZ*I_ERI_G2xyz_D2z_S_S_bb;
  Double I_ERI_G2x2z_F3z_S_S_bb = I_ERI_H2x3z_D2z_S_S_bb+ABZ*I_ERI_G2x2z_D2z_S_S_bb;
  Double I_ERI_Gx3y_F3z_S_S_bb = I_ERI_Hx3yz_D2z_S_S_bb+ABZ*I_ERI_Gx3y_D2z_S_S_bb;
  Double I_ERI_Gx2yz_F3z_S_S_bb = I_ERI_Hx2y2z_D2z_S_S_bb+ABZ*I_ERI_Gx2yz_D2z_S_S_bb;
  Double I_ERI_Gxy2z_F3z_S_S_bb = I_ERI_Hxy3z_D2z_S_S_bb+ABZ*I_ERI_Gxy2z_D2z_S_S_bb;
  Double I_ERI_Gx3z_F3z_S_S_bb = I_ERI_Hx4z_D2z_S_S_bb+ABZ*I_ERI_Gx3z_D2z_S_S_bb;
  Double I_ERI_G4y_F3z_S_S_bb = I_ERI_H4yz_D2z_S_S_bb+ABZ*I_ERI_G4y_D2z_S_S_bb;
  Double I_ERI_G3yz_F3z_S_S_bb = I_ERI_H3y2z_D2z_S_S_bb+ABZ*I_ERI_G3yz_D2z_S_S_bb;
  Double I_ERI_G2y2z_F3z_S_S_bb = I_ERI_H2y3z_D2z_S_S_bb+ABZ*I_ERI_G2y2z_D2z_S_S_bb;
  Double I_ERI_Gy3z_F3z_S_S_bb = I_ERI_Hy4z_D2z_S_S_bb+ABZ*I_ERI_Gy3z_D2z_S_S_bb;
  Double I_ERI_G4z_F3z_S_S_bb = I_ERI_H5z_D2z_S_S_bb+ABZ*I_ERI_G4z_D2z_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_G_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_F_S_S_bb
   ************************************************************/
  Double I_ERI_F3x_G4x_S_S_bb = I_ERI_G4x_F3x_S_S_bb+ABX*I_ERI_F3x_F3x_S_S_bb;
  Double I_ERI_F2xy_G4x_S_S_bb = I_ERI_G3xy_F3x_S_S_bb+ABX*I_ERI_F2xy_F3x_S_S_bb;
  Double I_ERI_F2xz_G4x_S_S_bb = I_ERI_G3xz_F3x_S_S_bb+ABX*I_ERI_F2xz_F3x_S_S_bb;
  Double I_ERI_Fx2y_G4x_S_S_bb = I_ERI_G2x2y_F3x_S_S_bb+ABX*I_ERI_Fx2y_F3x_S_S_bb;
  Double I_ERI_Fxyz_G4x_S_S_bb = I_ERI_G2xyz_F3x_S_S_bb+ABX*I_ERI_Fxyz_F3x_S_S_bb;
  Double I_ERI_Fx2z_G4x_S_S_bb = I_ERI_G2x2z_F3x_S_S_bb+ABX*I_ERI_Fx2z_F3x_S_S_bb;
  Double I_ERI_F3y_G4x_S_S_bb = I_ERI_Gx3y_F3x_S_S_bb+ABX*I_ERI_F3y_F3x_S_S_bb;
  Double I_ERI_F2yz_G4x_S_S_bb = I_ERI_Gx2yz_F3x_S_S_bb+ABX*I_ERI_F2yz_F3x_S_S_bb;
  Double I_ERI_Fy2z_G4x_S_S_bb = I_ERI_Gxy2z_F3x_S_S_bb+ABX*I_ERI_Fy2z_F3x_S_S_bb;
  Double I_ERI_F3z_G4x_S_S_bb = I_ERI_Gx3z_F3x_S_S_bb+ABX*I_ERI_F3z_F3x_S_S_bb;
  Double I_ERI_F3x_G3xy_S_S_bb = I_ERI_G3xy_F3x_S_S_bb+ABY*I_ERI_F3x_F3x_S_S_bb;
  Double I_ERI_F2xy_G3xy_S_S_bb = I_ERI_G2x2y_F3x_S_S_bb+ABY*I_ERI_F2xy_F3x_S_S_bb;
  Double I_ERI_F2xz_G3xy_S_S_bb = I_ERI_G2xyz_F3x_S_S_bb+ABY*I_ERI_F2xz_F3x_S_S_bb;
  Double I_ERI_Fx2y_G3xy_S_S_bb = I_ERI_Gx3y_F3x_S_S_bb+ABY*I_ERI_Fx2y_F3x_S_S_bb;
  Double I_ERI_Fxyz_G3xy_S_S_bb = I_ERI_Gx2yz_F3x_S_S_bb+ABY*I_ERI_Fxyz_F3x_S_S_bb;
  Double I_ERI_Fx2z_G3xy_S_S_bb = I_ERI_Gxy2z_F3x_S_S_bb+ABY*I_ERI_Fx2z_F3x_S_S_bb;
  Double I_ERI_F3y_G3xy_S_S_bb = I_ERI_G4y_F3x_S_S_bb+ABY*I_ERI_F3y_F3x_S_S_bb;
  Double I_ERI_F2yz_G3xy_S_S_bb = I_ERI_G3yz_F3x_S_S_bb+ABY*I_ERI_F2yz_F3x_S_S_bb;
  Double I_ERI_Fy2z_G3xy_S_S_bb = I_ERI_G2y2z_F3x_S_S_bb+ABY*I_ERI_Fy2z_F3x_S_S_bb;
  Double I_ERI_F3z_G3xy_S_S_bb = I_ERI_Gy3z_F3x_S_S_bb+ABY*I_ERI_F3z_F3x_S_S_bb;
  Double I_ERI_F3x_G3xz_S_S_bb = I_ERI_G3xz_F3x_S_S_bb+ABZ*I_ERI_F3x_F3x_S_S_bb;
  Double I_ERI_F2xy_G3xz_S_S_bb = I_ERI_G2xyz_F3x_S_S_bb+ABZ*I_ERI_F2xy_F3x_S_S_bb;
  Double I_ERI_F2xz_G3xz_S_S_bb = I_ERI_G2x2z_F3x_S_S_bb+ABZ*I_ERI_F2xz_F3x_S_S_bb;
  Double I_ERI_Fx2y_G3xz_S_S_bb = I_ERI_Gx2yz_F3x_S_S_bb+ABZ*I_ERI_Fx2y_F3x_S_S_bb;
  Double I_ERI_Fxyz_G3xz_S_S_bb = I_ERI_Gxy2z_F3x_S_S_bb+ABZ*I_ERI_Fxyz_F3x_S_S_bb;
  Double I_ERI_Fx2z_G3xz_S_S_bb = I_ERI_Gx3z_F3x_S_S_bb+ABZ*I_ERI_Fx2z_F3x_S_S_bb;
  Double I_ERI_F3y_G3xz_S_S_bb = I_ERI_G3yz_F3x_S_S_bb+ABZ*I_ERI_F3y_F3x_S_S_bb;
  Double I_ERI_F2yz_G3xz_S_S_bb = I_ERI_G2y2z_F3x_S_S_bb+ABZ*I_ERI_F2yz_F3x_S_S_bb;
  Double I_ERI_Fy2z_G3xz_S_S_bb = I_ERI_Gy3z_F3x_S_S_bb+ABZ*I_ERI_Fy2z_F3x_S_S_bb;
  Double I_ERI_F3z_G3xz_S_S_bb = I_ERI_G4z_F3x_S_S_bb+ABZ*I_ERI_F3z_F3x_S_S_bb;
  Double I_ERI_F3x_G2x2y_S_S_bb = I_ERI_G3xy_F2xy_S_S_bb+ABY*I_ERI_F3x_F2xy_S_S_bb;
  Double I_ERI_F2xy_G2x2y_S_S_bb = I_ERI_G2x2y_F2xy_S_S_bb+ABY*I_ERI_F2xy_F2xy_S_S_bb;
  Double I_ERI_F2xz_G2x2y_S_S_bb = I_ERI_G2xyz_F2xy_S_S_bb+ABY*I_ERI_F2xz_F2xy_S_S_bb;
  Double I_ERI_Fx2y_G2x2y_S_S_bb = I_ERI_Gx3y_F2xy_S_S_bb+ABY*I_ERI_Fx2y_F2xy_S_S_bb;
  Double I_ERI_Fxyz_G2x2y_S_S_bb = I_ERI_Gx2yz_F2xy_S_S_bb+ABY*I_ERI_Fxyz_F2xy_S_S_bb;
  Double I_ERI_Fx2z_G2x2y_S_S_bb = I_ERI_Gxy2z_F2xy_S_S_bb+ABY*I_ERI_Fx2z_F2xy_S_S_bb;
  Double I_ERI_F3y_G2x2y_S_S_bb = I_ERI_G4y_F2xy_S_S_bb+ABY*I_ERI_F3y_F2xy_S_S_bb;
  Double I_ERI_F2yz_G2x2y_S_S_bb = I_ERI_G3yz_F2xy_S_S_bb+ABY*I_ERI_F2yz_F2xy_S_S_bb;
  Double I_ERI_Fy2z_G2x2y_S_S_bb = I_ERI_G2y2z_F2xy_S_S_bb+ABY*I_ERI_Fy2z_F2xy_S_S_bb;
  Double I_ERI_F3z_G2x2y_S_S_bb = I_ERI_Gy3z_F2xy_S_S_bb+ABY*I_ERI_F3z_F2xy_S_S_bb;
  Double I_ERI_F3x_G2xyz_S_S_bb = I_ERI_G3xz_F2xy_S_S_bb+ABZ*I_ERI_F3x_F2xy_S_S_bb;
  Double I_ERI_F2xy_G2xyz_S_S_bb = I_ERI_G2xyz_F2xy_S_S_bb+ABZ*I_ERI_F2xy_F2xy_S_S_bb;
  Double I_ERI_F2xz_G2xyz_S_S_bb = I_ERI_G2x2z_F2xy_S_S_bb+ABZ*I_ERI_F2xz_F2xy_S_S_bb;
  Double I_ERI_Fx2y_G2xyz_S_S_bb = I_ERI_Gx2yz_F2xy_S_S_bb+ABZ*I_ERI_Fx2y_F2xy_S_S_bb;
  Double I_ERI_Fxyz_G2xyz_S_S_bb = I_ERI_Gxy2z_F2xy_S_S_bb+ABZ*I_ERI_Fxyz_F2xy_S_S_bb;
  Double I_ERI_Fx2z_G2xyz_S_S_bb = I_ERI_Gx3z_F2xy_S_S_bb+ABZ*I_ERI_Fx2z_F2xy_S_S_bb;
  Double I_ERI_F3y_G2xyz_S_S_bb = I_ERI_G3yz_F2xy_S_S_bb+ABZ*I_ERI_F3y_F2xy_S_S_bb;
  Double I_ERI_F2yz_G2xyz_S_S_bb = I_ERI_G2y2z_F2xy_S_S_bb+ABZ*I_ERI_F2yz_F2xy_S_S_bb;
  Double I_ERI_Fy2z_G2xyz_S_S_bb = I_ERI_Gy3z_F2xy_S_S_bb+ABZ*I_ERI_Fy2z_F2xy_S_S_bb;
  Double I_ERI_F3z_G2xyz_S_S_bb = I_ERI_G4z_F2xy_S_S_bb+ABZ*I_ERI_F3z_F2xy_S_S_bb;
  Double I_ERI_F3x_G2x2z_S_S_bb = I_ERI_G3xz_F2xz_S_S_bb+ABZ*I_ERI_F3x_F2xz_S_S_bb;
  Double I_ERI_F2xy_G2x2z_S_S_bb = I_ERI_G2xyz_F2xz_S_S_bb+ABZ*I_ERI_F2xy_F2xz_S_S_bb;
  Double I_ERI_F2xz_G2x2z_S_S_bb = I_ERI_G2x2z_F2xz_S_S_bb+ABZ*I_ERI_F2xz_F2xz_S_S_bb;
  Double I_ERI_Fx2y_G2x2z_S_S_bb = I_ERI_Gx2yz_F2xz_S_S_bb+ABZ*I_ERI_Fx2y_F2xz_S_S_bb;
  Double I_ERI_Fxyz_G2x2z_S_S_bb = I_ERI_Gxy2z_F2xz_S_S_bb+ABZ*I_ERI_Fxyz_F2xz_S_S_bb;
  Double I_ERI_Fx2z_G2x2z_S_S_bb = I_ERI_Gx3z_F2xz_S_S_bb+ABZ*I_ERI_Fx2z_F2xz_S_S_bb;
  Double I_ERI_F3y_G2x2z_S_S_bb = I_ERI_G3yz_F2xz_S_S_bb+ABZ*I_ERI_F3y_F2xz_S_S_bb;
  Double I_ERI_F2yz_G2x2z_S_S_bb = I_ERI_G2y2z_F2xz_S_S_bb+ABZ*I_ERI_F2yz_F2xz_S_S_bb;
  Double I_ERI_Fy2z_G2x2z_S_S_bb = I_ERI_Gy3z_F2xz_S_S_bb+ABZ*I_ERI_Fy2z_F2xz_S_S_bb;
  Double I_ERI_F3z_G2x2z_S_S_bb = I_ERI_G4z_F2xz_S_S_bb+ABZ*I_ERI_F3z_F2xz_S_S_bb;
  Double I_ERI_F3x_Gx3y_S_S_bb = I_ERI_G4x_F3y_S_S_bb+ABX*I_ERI_F3x_F3y_S_S_bb;
  Double I_ERI_F2xy_Gx3y_S_S_bb = I_ERI_G3xy_F3y_S_S_bb+ABX*I_ERI_F2xy_F3y_S_S_bb;
  Double I_ERI_F2xz_Gx3y_S_S_bb = I_ERI_G3xz_F3y_S_S_bb+ABX*I_ERI_F2xz_F3y_S_S_bb;
  Double I_ERI_Fx2y_Gx3y_S_S_bb = I_ERI_G2x2y_F3y_S_S_bb+ABX*I_ERI_Fx2y_F3y_S_S_bb;
  Double I_ERI_Fxyz_Gx3y_S_S_bb = I_ERI_G2xyz_F3y_S_S_bb+ABX*I_ERI_Fxyz_F3y_S_S_bb;
  Double I_ERI_Fx2z_Gx3y_S_S_bb = I_ERI_G2x2z_F3y_S_S_bb+ABX*I_ERI_Fx2z_F3y_S_S_bb;
  Double I_ERI_F3y_Gx3y_S_S_bb = I_ERI_Gx3y_F3y_S_S_bb+ABX*I_ERI_F3y_F3y_S_S_bb;
  Double I_ERI_F2yz_Gx3y_S_S_bb = I_ERI_Gx2yz_F3y_S_S_bb+ABX*I_ERI_F2yz_F3y_S_S_bb;
  Double I_ERI_Fy2z_Gx3y_S_S_bb = I_ERI_Gxy2z_F3y_S_S_bb+ABX*I_ERI_Fy2z_F3y_S_S_bb;
  Double I_ERI_F3z_Gx3y_S_S_bb = I_ERI_Gx3z_F3y_S_S_bb+ABX*I_ERI_F3z_F3y_S_S_bb;
  Double I_ERI_F3x_Gx2yz_S_S_bb = I_ERI_G3xz_Fx2y_S_S_bb+ABZ*I_ERI_F3x_Fx2y_S_S_bb;
  Double I_ERI_F2xy_Gx2yz_S_S_bb = I_ERI_G2xyz_Fx2y_S_S_bb+ABZ*I_ERI_F2xy_Fx2y_S_S_bb;
  Double I_ERI_F2xz_Gx2yz_S_S_bb = I_ERI_G2x2z_Fx2y_S_S_bb+ABZ*I_ERI_F2xz_Fx2y_S_S_bb;
  Double I_ERI_Fx2y_Gx2yz_S_S_bb = I_ERI_Gx2yz_Fx2y_S_S_bb+ABZ*I_ERI_Fx2y_Fx2y_S_S_bb;
  Double I_ERI_Fxyz_Gx2yz_S_S_bb = I_ERI_Gxy2z_Fx2y_S_S_bb+ABZ*I_ERI_Fxyz_Fx2y_S_S_bb;
  Double I_ERI_Fx2z_Gx2yz_S_S_bb = I_ERI_Gx3z_Fx2y_S_S_bb+ABZ*I_ERI_Fx2z_Fx2y_S_S_bb;
  Double I_ERI_F3y_Gx2yz_S_S_bb = I_ERI_G3yz_Fx2y_S_S_bb+ABZ*I_ERI_F3y_Fx2y_S_S_bb;
  Double I_ERI_F2yz_Gx2yz_S_S_bb = I_ERI_G2y2z_Fx2y_S_S_bb+ABZ*I_ERI_F2yz_Fx2y_S_S_bb;
  Double I_ERI_Fy2z_Gx2yz_S_S_bb = I_ERI_Gy3z_Fx2y_S_S_bb+ABZ*I_ERI_Fy2z_Fx2y_S_S_bb;
  Double I_ERI_F3z_Gx2yz_S_S_bb = I_ERI_G4z_Fx2y_S_S_bb+ABZ*I_ERI_F3z_Fx2y_S_S_bb;
  Double I_ERI_F3x_Gxy2z_S_S_bb = I_ERI_G3xy_Fx2z_S_S_bb+ABY*I_ERI_F3x_Fx2z_S_S_bb;
  Double I_ERI_F2xy_Gxy2z_S_S_bb = I_ERI_G2x2y_Fx2z_S_S_bb+ABY*I_ERI_F2xy_Fx2z_S_S_bb;
  Double I_ERI_F2xz_Gxy2z_S_S_bb = I_ERI_G2xyz_Fx2z_S_S_bb+ABY*I_ERI_F2xz_Fx2z_S_S_bb;
  Double I_ERI_Fx2y_Gxy2z_S_S_bb = I_ERI_Gx3y_Fx2z_S_S_bb+ABY*I_ERI_Fx2y_Fx2z_S_S_bb;
  Double I_ERI_Fxyz_Gxy2z_S_S_bb = I_ERI_Gx2yz_Fx2z_S_S_bb+ABY*I_ERI_Fxyz_Fx2z_S_S_bb;
  Double I_ERI_Fx2z_Gxy2z_S_S_bb = I_ERI_Gxy2z_Fx2z_S_S_bb+ABY*I_ERI_Fx2z_Fx2z_S_S_bb;
  Double I_ERI_F3y_Gxy2z_S_S_bb = I_ERI_G4y_Fx2z_S_S_bb+ABY*I_ERI_F3y_Fx2z_S_S_bb;
  Double I_ERI_F2yz_Gxy2z_S_S_bb = I_ERI_G3yz_Fx2z_S_S_bb+ABY*I_ERI_F2yz_Fx2z_S_S_bb;
  Double I_ERI_Fy2z_Gxy2z_S_S_bb = I_ERI_G2y2z_Fx2z_S_S_bb+ABY*I_ERI_Fy2z_Fx2z_S_S_bb;
  Double I_ERI_F3z_Gxy2z_S_S_bb = I_ERI_Gy3z_Fx2z_S_S_bb+ABY*I_ERI_F3z_Fx2z_S_S_bb;
  Double I_ERI_F3x_Gx3z_S_S_bb = I_ERI_G4x_F3z_S_S_bb+ABX*I_ERI_F3x_F3z_S_S_bb;
  Double I_ERI_F2xy_Gx3z_S_S_bb = I_ERI_G3xy_F3z_S_S_bb+ABX*I_ERI_F2xy_F3z_S_S_bb;
  Double I_ERI_F2xz_Gx3z_S_S_bb = I_ERI_G3xz_F3z_S_S_bb+ABX*I_ERI_F2xz_F3z_S_S_bb;
  Double I_ERI_Fx2y_Gx3z_S_S_bb = I_ERI_G2x2y_F3z_S_S_bb+ABX*I_ERI_Fx2y_F3z_S_S_bb;
  Double I_ERI_Fxyz_Gx3z_S_S_bb = I_ERI_G2xyz_F3z_S_S_bb+ABX*I_ERI_Fxyz_F3z_S_S_bb;
  Double I_ERI_Fx2z_Gx3z_S_S_bb = I_ERI_G2x2z_F3z_S_S_bb+ABX*I_ERI_Fx2z_F3z_S_S_bb;
  Double I_ERI_F3y_Gx3z_S_S_bb = I_ERI_Gx3y_F3z_S_S_bb+ABX*I_ERI_F3y_F3z_S_S_bb;
  Double I_ERI_F2yz_Gx3z_S_S_bb = I_ERI_Gx2yz_F3z_S_S_bb+ABX*I_ERI_F2yz_F3z_S_S_bb;
  Double I_ERI_Fy2z_Gx3z_S_S_bb = I_ERI_Gxy2z_F3z_S_S_bb+ABX*I_ERI_Fy2z_F3z_S_S_bb;
  Double I_ERI_F3z_Gx3z_S_S_bb = I_ERI_Gx3z_F3z_S_S_bb+ABX*I_ERI_F3z_F3z_S_S_bb;
  Double I_ERI_F3x_G4y_S_S_bb = I_ERI_G3xy_F3y_S_S_bb+ABY*I_ERI_F3x_F3y_S_S_bb;
  Double I_ERI_F2xy_G4y_S_S_bb = I_ERI_G2x2y_F3y_S_S_bb+ABY*I_ERI_F2xy_F3y_S_S_bb;
  Double I_ERI_F2xz_G4y_S_S_bb = I_ERI_G2xyz_F3y_S_S_bb+ABY*I_ERI_F2xz_F3y_S_S_bb;
  Double I_ERI_Fx2y_G4y_S_S_bb = I_ERI_Gx3y_F3y_S_S_bb+ABY*I_ERI_Fx2y_F3y_S_S_bb;
  Double I_ERI_Fxyz_G4y_S_S_bb = I_ERI_Gx2yz_F3y_S_S_bb+ABY*I_ERI_Fxyz_F3y_S_S_bb;
  Double I_ERI_Fx2z_G4y_S_S_bb = I_ERI_Gxy2z_F3y_S_S_bb+ABY*I_ERI_Fx2z_F3y_S_S_bb;
  Double I_ERI_F3y_G4y_S_S_bb = I_ERI_G4y_F3y_S_S_bb+ABY*I_ERI_F3y_F3y_S_S_bb;
  Double I_ERI_F2yz_G4y_S_S_bb = I_ERI_G3yz_F3y_S_S_bb+ABY*I_ERI_F2yz_F3y_S_S_bb;
  Double I_ERI_Fy2z_G4y_S_S_bb = I_ERI_G2y2z_F3y_S_S_bb+ABY*I_ERI_Fy2z_F3y_S_S_bb;
  Double I_ERI_F3z_G4y_S_S_bb = I_ERI_Gy3z_F3y_S_S_bb+ABY*I_ERI_F3z_F3y_S_S_bb;
  Double I_ERI_F3x_G3yz_S_S_bb = I_ERI_G3xz_F3y_S_S_bb+ABZ*I_ERI_F3x_F3y_S_S_bb;
  Double I_ERI_F2xy_G3yz_S_S_bb = I_ERI_G2xyz_F3y_S_S_bb+ABZ*I_ERI_F2xy_F3y_S_S_bb;
  Double I_ERI_F2xz_G3yz_S_S_bb = I_ERI_G2x2z_F3y_S_S_bb+ABZ*I_ERI_F2xz_F3y_S_S_bb;
  Double I_ERI_Fx2y_G3yz_S_S_bb = I_ERI_Gx2yz_F3y_S_S_bb+ABZ*I_ERI_Fx2y_F3y_S_S_bb;
  Double I_ERI_Fxyz_G3yz_S_S_bb = I_ERI_Gxy2z_F3y_S_S_bb+ABZ*I_ERI_Fxyz_F3y_S_S_bb;
  Double I_ERI_Fx2z_G3yz_S_S_bb = I_ERI_Gx3z_F3y_S_S_bb+ABZ*I_ERI_Fx2z_F3y_S_S_bb;
  Double I_ERI_F3y_G3yz_S_S_bb = I_ERI_G3yz_F3y_S_S_bb+ABZ*I_ERI_F3y_F3y_S_S_bb;
  Double I_ERI_F2yz_G3yz_S_S_bb = I_ERI_G2y2z_F3y_S_S_bb+ABZ*I_ERI_F2yz_F3y_S_S_bb;
  Double I_ERI_Fy2z_G3yz_S_S_bb = I_ERI_Gy3z_F3y_S_S_bb+ABZ*I_ERI_Fy2z_F3y_S_S_bb;
  Double I_ERI_F3z_G3yz_S_S_bb = I_ERI_G4z_F3y_S_S_bb+ABZ*I_ERI_F3z_F3y_S_S_bb;
  Double I_ERI_F3x_G2y2z_S_S_bb = I_ERI_G3xz_F2yz_S_S_bb+ABZ*I_ERI_F3x_F2yz_S_S_bb;
  Double I_ERI_F2xy_G2y2z_S_S_bb = I_ERI_G2xyz_F2yz_S_S_bb+ABZ*I_ERI_F2xy_F2yz_S_S_bb;
  Double I_ERI_F2xz_G2y2z_S_S_bb = I_ERI_G2x2z_F2yz_S_S_bb+ABZ*I_ERI_F2xz_F2yz_S_S_bb;
  Double I_ERI_Fx2y_G2y2z_S_S_bb = I_ERI_Gx2yz_F2yz_S_S_bb+ABZ*I_ERI_Fx2y_F2yz_S_S_bb;
  Double I_ERI_Fxyz_G2y2z_S_S_bb = I_ERI_Gxy2z_F2yz_S_S_bb+ABZ*I_ERI_Fxyz_F2yz_S_S_bb;
  Double I_ERI_Fx2z_G2y2z_S_S_bb = I_ERI_Gx3z_F2yz_S_S_bb+ABZ*I_ERI_Fx2z_F2yz_S_S_bb;
  Double I_ERI_F3y_G2y2z_S_S_bb = I_ERI_G3yz_F2yz_S_S_bb+ABZ*I_ERI_F3y_F2yz_S_S_bb;
  Double I_ERI_F2yz_G2y2z_S_S_bb = I_ERI_G2y2z_F2yz_S_S_bb+ABZ*I_ERI_F2yz_F2yz_S_S_bb;
  Double I_ERI_Fy2z_G2y2z_S_S_bb = I_ERI_Gy3z_F2yz_S_S_bb+ABZ*I_ERI_Fy2z_F2yz_S_S_bb;
  Double I_ERI_F3z_G2y2z_S_S_bb = I_ERI_G4z_F2yz_S_S_bb+ABZ*I_ERI_F3z_F2yz_S_S_bb;
  Double I_ERI_F3x_Gy3z_S_S_bb = I_ERI_G3xy_F3z_S_S_bb+ABY*I_ERI_F3x_F3z_S_S_bb;
  Double I_ERI_F2xy_Gy3z_S_S_bb = I_ERI_G2x2y_F3z_S_S_bb+ABY*I_ERI_F2xy_F3z_S_S_bb;
  Double I_ERI_F2xz_Gy3z_S_S_bb = I_ERI_G2xyz_F3z_S_S_bb+ABY*I_ERI_F2xz_F3z_S_S_bb;
  Double I_ERI_Fx2y_Gy3z_S_S_bb = I_ERI_Gx3y_F3z_S_S_bb+ABY*I_ERI_Fx2y_F3z_S_S_bb;
  Double I_ERI_Fxyz_Gy3z_S_S_bb = I_ERI_Gx2yz_F3z_S_S_bb+ABY*I_ERI_Fxyz_F3z_S_S_bb;
  Double I_ERI_Fx2z_Gy3z_S_S_bb = I_ERI_Gxy2z_F3z_S_S_bb+ABY*I_ERI_Fx2z_F3z_S_S_bb;
  Double I_ERI_F3y_Gy3z_S_S_bb = I_ERI_G4y_F3z_S_S_bb+ABY*I_ERI_F3y_F3z_S_S_bb;
  Double I_ERI_F2yz_Gy3z_S_S_bb = I_ERI_G3yz_F3z_S_S_bb+ABY*I_ERI_F2yz_F3z_S_S_bb;
  Double I_ERI_Fy2z_Gy3z_S_S_bb = I_ERI_G2y2z_F3z_S_S_bb+ABY*I_ERI_Fy2z_F3z_S_S_bb;
  Double I_ERI_F3z_Gy3z_S_S_bb = I_ERI_Gy3z_F3z_S_S_bb+ABY*I_ERI_F3z_F3z_S_S_bb;
  Double I_ERI_F3x_G4z_S_S_bb = I_ERI_G3xz_F3z_S_S_bb+ABZ*I_ERI_F3x_F3z_S_S_bb;
  Double I_ERI_F2xy_G4z_S_S_bb = I_ERI_G2xyz_F3z_S_S_bb+ABZ*I_ERI_F2xy_F3z_S_S_bb;
  Double I_ERI_F2xz_G4z_S_S_bb = I_ERI_G2x2z_F3z_S_S_bb+ABZ*I_ERI_F2xz_F3z_S_S_bb;
  Double I_ERI_Fx2y_G4z_S_S_bb = I_ERI_Gx2yz_F3z_S_S_bb+ABZ*I_ERI_Fx2y_F3z_S_S_bb;
  Double I_ERI_Fxyz_G4z_S_S_bb = I_ERI_Gxy2z_F3z_S_S_bb+ABZ*I_ERI_Fxyz_F3z_S_S_bb;
  Double I_ERI_Fx2z_G4z_S_S_bb = I_ERI_Gx3z_F3z_S_S_bb+ABZ*I_ERI_Fx2z_F3z_S_S_bb;
  Double I_ERI_F3y_G4z_S_S_bb = I_ERI_G3yz_F3z_S_S_bb+ABZ*I_ERI_F3y_F3z_S_S_bb;
  Double I_ERI_F2yz_G4z_S_S_bb = I_ERI_G2y2z_F3z_S_S_bb+ABZ*I_ERI_F2yz_F3z_S_S_bb;
  Double I_ERI_Fy2z_G4z_S_S_bb = I_ERI_Gy3z_F3z_S_S_bb+ABZ*I_ERI_Fy2z_F3z_S_S_bb;
  Double I_ERI_F3z_G4z_S_S_bb = I_ERI_G4z_F3z_S_S_bb+ABZ*I_ERI_F3z_F3z_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_P_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_S_P_S_bc
   ************************************************************/
  Double I_ERI_F3x_Px_Px_S_bc = I_ERI_G4x_S_Px_S_bc+ABX*I_ERI_F3x_S_Px_S_bc;
  Double I_ERI_F2xy_Px_Px_S_bc = I_ERI_G3xy_S_Px_S_bc+ABX*I_ERI_F2xy_S_Px_S_bc;
  Double I_ERI_F2xz_Px_Px_S_bc = I_ERI_G3xz_S_Px_S_bc+ABX*I_ERI_F2xz_S_Px_S_bc;
  Double I_ERI_Fx2y_Px_Px_S_bc = I_ERI_G2x2y_S_Px_S_bc+ABX*I_ERI_Fx2y_S_Px_S_bc;
  Double I_ERI_Fxyz_Px_Px_S_bc = I_ERI_G2xyz_S_Px_S_bc+ABX*I_ERI_Fxyz_S_Px_S_bc;
  Double I_ERI_Fx2z_Px_Px_S_bc = I_ERI_G2x2z_S_Px_S_bc+ABX*I_ERI_Fx2z_S_Px_S_bc;
  Double I_ERI_F3y_Px_Px_S_bc = I_ERI_Gx3y_S_Px_S_bc+ABX*I_ERI_F3y_S_Px_S_bc;
  Double I_ERI_F2yz_Px_Px_S_bc = I_ERI_Gx2yz_S_Px_S_bc+ABX*I_ERI_F2yz_S_Px_S_bc;
  Double I_ERI_Fy2z_Px_Px_S_bc = I_ERI_Gxy2z_S_Px_S_bc+ABX*I_ERI_Fy2z_S_Px_S_bc;
  Double I_ERI_F3z_Px_Px_S_bc = I_ERI_Gx3z_S_Px_S_bc+ABX*I_ERI_F3z_S_Px_S_bc;
  Double I_ERI_F3x_Py_Px_S_bc = I_ERI_G3xy_S_Px_S_bc+ABY*I_ERI_F3x_S_Px_S_bc;
  Double I_ERI_F2xy_Py_Px_S_bc = I_ERI_G2x2y_S_Px_S_bc+ABY*I_ERI_F2xy_S_Px_S_bc;
  Double I_ERI_F2xz_Py_Px_S_bc = I_ERI_G2xyz_S_Px_S_bc+ABY*I_ERI_F2xz_S_Px_S_bc;
  Double I_ERI_Fx2y_Py_Px_S_bc = I_ERI_Gx3y_S_Px_S_bc+ABY*I_ERI_Fx2y_S_Px_S_bc;
  Double I_ERI_Fxyz_Py_Px_S_bc = I_ERI_Gx2yz_S_Px_S_bc+ABY*I_ERI_Fxyz_S_Px_S_bc;
  Double I_ERI_Fx2z_Py_Px_S_bc = I_ERI_Gxy2z_S_Px_S_bc+ABY*I_ERI_Fx2z_S_Px_S_bc;
  Double I_ERI_F3y_Py_Px_S_bc = I_ERI_G4y_S_Px_S_bc+ABY*I_ERI_F3y_S_Px_S_bc;
  Double I_ERI_F2yz_Py_Px_S_bc = I_ERI_G3yz_S_Px_S_bc+ABY*I_ERI_F2yz_S_Px_S_bc;
  Double I_ERI_Fy2z_Py_Px_S_bc = I_ERI_G2y2z_S_Px_S_bc+ABY*I_ERI_Fy2z_S_Px_S_bc;
  Double I_ERI_F3z_Py_Px_S_bc = I_ERI_Gy3z_S_Px_S_bc+ABY*I_ERI_F3z_S_Px_S_bc;
  Double I_ERI_F3x_Pz_Px_S_bc = I_ERI_G3xz_S_Px_S_bc+ABZ*I_ERI_F3x_S_Px_S_bc;
  Double I_ERI_F2xy_Pz_Px_S_bc = I_ERI_G2xyz_S_Px_S_bc+ABZ*I_ERI_F2xy_S_Px_S_bc;
  Double I_ERI_F2xz_Pz_Px_S_bc = I_ERI_G2x2z_S_Px_S_bc+ABZ*I_ERI_F2xz_S_Px_S_bc;
  Double I_ERI_Fx2y_Pz_Px_S_bc = I_ERI_Gx2yz_S_Px_S_bc+ABZ*I_ERI_Fx2y_S_Px_S_bc;
  Double I_ERI_Fxyz_Pz_Px_S_bc = I_ERI_Gxy2z_S_Px_S_bc+ABZ*I_ERI_Fxyz_S_Px_S_bc;
  Double I_ERI_Fx2z_Pz_Px_S_bc = I_ERI_Gx3z_S_Px_S_bc+ABZ*I_ERI_Fx2z_S_Px_S_bc;
  Double I_ERI_F3y_Pz_Px_S_bc = I_ERI_G3yz_S_Px_S_bc+ABZ*I_ERI_F3y_S_Px_S_bc;
  Double I_ERI_F2yz_Pz_Px_S_bc = I_ERI_G2y2z_S_Px_S_bc+ABZ*I_ERI_F2yz_S_Px_S_bc;
  Double I_ERI_Fy2z_Pz_Px_S_bc = I_ERI_Gy3z_S_Px_S_bc+ABZ*I_ERI_Fy2z_S_Px_S_bc;
  Double I_ERI_F3z_Pz_Px_S_bc = I_ERI_G4z_S_Px_S_bc+ABZ*I_ERI_F3z_S_Px_S_bc;
  Double I_ERI_F3x_Px_Py_S_bc = I_ERI_G4x_S_Py_S_bc+ABX*I_ERI_F3x_S_Py_S_bc;
  Double I_ERI_F2xy_Px_Py_S_bc = I_ERI_G3xy_S_Py_S_bc+ABX*I_ERI_F2xy_S_Py_S_bc;
  Double I_ERI_F2xz_Px_Py_S_bc = I_ERI_G3xz_S_Py_S_bc+ABX*I_ERI_F2xz_S_Py_S_bc;
  Double I_ERI_Fx2y_Px_Py_S_bc = I_ERI_G2x2y_S_Py_S_bc+ABX*I_ERI_Fx2y_S_Py_S_bc;
  Double I_ERI_Fxyz_Px_Py_S_bc = I_ERI_G2xyz_S_Py_S_bc+ABX*I_ERI_Fxyz_S_Py_S_bc;
  Double I_ERI_Fx2z_Px_Py_S_bc = I_ERI_G2x2z_S_Py_S_bc+ABX*I_ERI_Fx2z_S_Py_S_bc;
  Double I_ERI_F3y_Px_Py_S_bc = I_ERI_Gx3y_S_Py_S_bc+ABX*I_ERI_F3y_S_Py_S_bc;
  Double I_ERI_F2yz_Px_Py_S_bc = I_ERI_Gx2yz_S_Py_S_bc+ABX*I_ERI_F2yz_S_Py_S_bc;
  Double I_ERI_Fy2z_Px_Py_S_bc = I_ERI_Gxy2z_S_Py_S_bc+ABX*I_ERI_Fy2z_S_Py_S_bc;
  Double I_ERI_F3z_Px_Py_S_bc = I_ERI_Gx3z_S_Py_S_bc+ABX*I_ERI_F3z_S_Py_S_bc;
  Double I_ERI_F3x_Py_Py_S_bc = I_ERI_G3xy_S_Py_S_bc+ABY*I_ERI_F3x_S_Py_S_bc;
  Double I_ERI_F2xy_Py_Py_S_bc = I_ERI_G2x2y_S_Py_S_bc+ABY*I_ERI_F2xy_S_Py_S_bc;
  Double I_ERI_F2xz_Py_Py_S_bc = I_ERI_G2xyz_S_Py_S_bc+ABY*I_ERI_F2xz_S_Py_S_bc;
  Double I_ERI_Fx2y_Py_Py_S_bc = I_ERI_Gx3y_S_Py_S_bc+ABY*I_ERI_Fx2y_S_Py_S_bc;
  Double I_ERI_Fxyz_Py_Py_S_bc = I_ERI_Gx2yz_S_Py_S_bc+ABY*I_ERI_Fxyz_S_Py_S_bc;
  Double I_ERI_Fx2z_Py_Py_S_bc = I_ERI_Gxy2z_S_Py_S_bc+ABY*I_ERI_Fx2z_S_Py_S_bc;
  Double I_ERI_F3y_Py_Py_S_bc = I_ERI_G4y_S_Py_S_bc+ABY*I_ERI_F3y_S_Py_S_bc;
  Double I_ERI_F2yz_Py_Py_S_bc = I_ERI_G3yz_S_Py_S_bc+ABY*I_ERI_F2yz_S_Py_S_bc;
  Double I_ERI_Fy2z_Py_Py_S_bc = I_ERI_G2y2z_S_Py_S_bc+ABY*I_ERI_Fy2z_S_Py_S_bc;
  Double I_ERI_F3z_Py_Py_S_bc = I_ERI_Gy3z_S_Py_S_bc+ABY*I_ERI_F3z_S_Py_S_bc;
  Double I_ERI_F3x_Pz_Py_S_bc = I_ERI_G3xz_S_Py_S_bc+ABZ*I_ERI_F3x_S_Py_S_bc;
  Double I_ERI_F2xy_Pz_Py_S_bc = I_ERI_G2xyz_S_Py_S_bc+ABZ*I_ERI_F2xy_S_Py_S_bc;
  Double I_ERI_F2xz_Pz_Py_S_bc = I_ERI_G2x2z_S_Py_S_bc+ABZ*I_ERI_F2xz_S_Py_S_bc;
  Double I_ERI_Fx2y_Pz_Py_S_bc = I_ERI_Gx2yz_S_Py_S_bc+ABZ*I_ERI_Fx2y_S_Py_S_bc;
  Double I_ERI_Fxyz_Pz_Py_S_bc = I_ERI_Gxy2z_S_Py_S_bc+ABZ*I_ERI_Fxyz_S_Py_S_bc;
  Double I_ERI_Fx2z_Pz_Py_S_bc = I_ERI_Gx3z_S_Py_S_bc+ABZ*I_ERI_Fx2z_S_Py_S_bc;
  Double I_ERI_F3y_Pz_Py_S_bc = I_ERI_G3yz_S_Py_S_bc+ABZ*I_ERI_F3y_S_Py_S_bc;
  Double I_ERI_F2yz_Pz_Py_S_bc = I_ERI_G2y2z_S_Py_S_bc+ABZ*I_ERI_F2yz_S_Py_S_bc;
  Double I_ERI_Fy2z_Pz_Py_S_bc = I_ERI_Gy3z_S_Py_S_bc+ABZ*I_ERI_Fy2z_S_Py_S_bc;
  Double I_ERI_F3z_Pz_Py_S_bc = I_ERI_G4z_S_Py_S_bc+ABZ*I_ERI_F3z_S_Py_S_bc;
  Double I_ERI_F3x_Px_Pz_S_bc = I_ERI_G4x_S_Pz_S_bc+ABX*I_ERI_F3x_S_Pz_S_bc;
  Double I_ERI_F2xy_Px_Pz_S_bc = I_ERI_G3xy_S_Pz_S_bc+ABX*I_ERI_F2xy_S_Pz_S_bc;
  Double I_ERI_F2xz_Px_Pz_S_bc = I_ERI_G3xz_S_Pz_S_bc+ABX*I_ERI_F2xz_S_Pz_S_bc;
  Double I_ERI_Fx2y_Px_Pz_S_bc = I_ERI_G2x2y_S_Pz_S_bc+ABX*I_ERI_Fx2y_S_Pz_S_bc;
  Double I_ERI_Fxyz_Px_Pz_S_bc = I_ERI_G2xyz_S_Pz_S_bc+ABX*I_ERI_Fxyz_S_Pz_S_bc;
  Double I_ERI_Fx2z_Px_Pz_S_bc = I_ERI_G2x2z_S_Pz_S_bc+ABX*I_ERI_Fx2z_S_Pz_S_bc;
  Double I_ERI_F3y_Px_Pz_S_bc = I_ERI_Gx3y_S_Pz_S_bc+ABX*I_ERI_F3y_S_Pz_S_bc;
  Double I_ERI_F2yz_Px_Pz_S_bc = I_ERI_Gx2yz_S_Pz_S_bc+ABX*I_ERI_F2yz_S_Pz_S_bc;
  Double I_ERI_Fy2z_Px_Pz_S_bc = I_ERI_Gxy2z_S_Pz_S_bc+ABX*I_ERI_Fy2z_S_Pz_S_bc;
  Double I_ERI_F3z_Px_Pz_S_bc = I_ERI_Gx3z_S_Pz_S_bc+ABX*I_ERI_F3z_S_Pz_S_bc;
  Double I_ERI_F3x_Py_Pz_S_bc = I_ERI_G3xy_S_Pz_S_bc+ABY*I_ERI_F3x_S_Pz_S_bc;
  Double I_ERI_F2xy_Py_Pz_S_bc = I_ERI_G2x2y_S_Pz_S_bc+ABY*I_ERI_F2xy_S_Pz_S_bc;
  Double I_ERI_F2xz_Py_Pz_S_bc = I_ERI_G2xyz_S_Pz_S_bc+ABY*I_ERI_F2xz_S_Pz_S_bc;
  Double I_ERI_Fx2y_Py_Pz_S_bc = I_ERI_Gx3y_S_Pz_S_bc+ABY*I_ERI_Fx2y_S_Pz_S_bc;
  Double I_ERI_Fxyz_Py_Pz_S_bc = I_ERI_Gx2yz_S_Pz_S_bc+ABY*I_ERI_Fxyz_S_Pz_S_bc;
  Double I_ERI_Fx2z_Py_Pz_S_bc = I_ERI_Gxy2z_S_Pz_S_bc+ABY*I_ERI_Fx2z_S_Pz_S_bc;
  Double I_ERI_F3y_Py_Pz_S_bc = I_ERI_G4y_S_Pz_S_bc+ABY*I_ERI_F3y_S_Pz_S_bc;
  Double I_ERI_F2yz_Py_Pz_S_bc = I_ERI_G3yz_S_Pz_S_bc+ABY*I_ERI_F2yz_S_Pz_S_bc;
  Double I_ERI_Fy2z_Py_Pz_S_bc = I_ERI_G2y2z_S_Pz_S_bc+ABY*I_ERI_Fy2z_S_Pz_S_bc;
  Double I_ERI_F3z_Py_Pz_S_bc = I_ERI_Gy3z_S_Pz_S_bc+ABY*I_ERI_F3z_S_Pz_S_bc;
  Double I_ERI_F3x_Pz_Pz_S_bc = I_ERI_G3xz_S_Pz_S_bc+ABZ*I_ERI_F3x_S_Pz_S_bc;
  Double I_ERI_F2xy_Pz_Pz_S_bc = I_ERI_G2xyz_S_Pz_S_bc+ABZ*I_ERI_F2xy_S_Pz_S_bc;
  Double I_ERI_F2xz_Pz_Pz_S_bc = I_ERI_G2x2z_S_Pz_S_bc+ABZ*I_ERI_F2xz_S_Pz_S_bc;
  Double I_ERI_Fx2y_Pz_Pz_S_bc = I_ERI_Gx2yz_S_Pz_S_bc+ABZ*I_ERI_Fx2y_S_Pz_S_bc;
  Double I_ERI_Fxyz_Pz_Pz_S_bc = I_ERI_Gxy2z_S_Pz_S_bc+ABZ*I_ERI_Fxyz_S_Pz_S_bc;
  Double I_ERI_Fx2z_Pz_Pz_S_bc = I_ERI_Gx3z_S_Pz_S_bc+ABZ*I_ERI_Fx2z_S_Pz_S_bc;
  Double I_ERI_F3y_Pz_Pz_S_bc = I_ERI_G3yz_S_Pz_S_bc+ABZ*I_ERI_F3y_S_Pz_S_bc;
  Double I_ERI_F2yz_Pz_Pz_S_bc = I_ERI_G2y2z_S_Pz_S_bc+ABZ*I_ERI_F2yz_S_Pz_S_bc;
  Double I_ERI_Fy2z_Pz_Pz_S_bc = I_ERI_Gy3z_S_Pz_S_bc+ABZ*I_ERI_Fy2z_S_Pz_S_bc;
  Double I_ERI_F3z_Pz_Pz_S_bc = I_ERI_G4z_S_Pz_S_bc+ABZ*I_ERI_F3z_S_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_P_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_P_S_bc
   * RHS shell quartet name: SQ_ERI_G_S_P_S_bc
   ************************************************************/
  Double I_ERI_G4x_Px_Px_S_bc = I_ERI_H5x_S_Px_S_bc+ABX*I_ERI_G4x_S_Px_S_bc;
  Double I_ERI_G3xy_Px_Px_S_bc = I_ERI_H4xy_S_Px_S_bc+ABX*I_ERI_G3xy_S_Px_S_bc;
  Double I_ERI_G3xz_Px_Px_S_bc = I_ERI_H4xz_S_Px_S_bc+ABX*I_ERI_G3xz_S_Px_S_bc;
  Double I_ERI_G2x2y_Px_Px_S_bc = I_ERI_H3x2y_S_Px_S_bc+ABX*I_ERI_G2x2y_S_Px_S_bc;
  Double I_ERI_G2xyz_Px_Px_S_bc = I_ERI_H3xyz_S_Px_S_bc+ABX*I_ERI_G2xyz_S_Px_S_bc;
  Double I_ERI_G2x2z_Px_Px_S_bc = I_ERI_H3x2z_S_Px_S_bc+ABX*I_ERI_G2x2z_S_Px_S_bc;
  Double I_ERI_Gx3y_Px_Px_S_bc = I_ERI_H2x3y_S_Px_S_bc+ABX*I_ERI_Gx3y_S_Px_S_bc;
  Double I_ERI_Gx2yz_Px_Px_S_bc = I_ERI_H2x2yz_S_Px_S_bc+ABX*I_ERI_Gx2yz_S_Px_S_bc;
  Double I_ERI_Gxy2z_Px_Px_S_bc = I_ERI_H2xy2z_S_Px_S_bc+ABX*I_ERI_Gxy2z_S_Px_S_bc;
  Double I_ERI_Gx3z_Px_Px_S_bc = I_ERI_H2x3z_S_Px_S_bc+ABX*I_ERI_Gx3z_S_Px_S_bc;
  Double I_ERI_G4y_Px_Px_S_bc = I_ERI_Hx4y_S_Px_S_bc+ABX*I_ERI_G4y_S_Px_S_bc;
  Double I_ERI_G3yz_Px_Px_S_bc = I_ERI_Hx3yz_S_Px_S_bc+ABX*I_ERI_G3yz_S_Px_S_bc;
  Double I_ERI_G2y2z_Px_Px_S_bc = I_ERI_Hx2y2z_S_Px_S_bc+ABX*I_ERI_G2y2z_S_Px_S_bc;
  Double I_ERI_Gy3z_Px_Px_S_bc = I_ERI_Hxy3z_S_Px_S_bc+ABX*I_ERI_Gy3z_S_Px_S_bc;
  Double I_ERI_G4z_Px_Px_S_bc = I_ERI_Hx4z_S_Px_S_bc+ABX*I_ERI_G4z_S_Px_S_bc;
  Double I_ERI_G4x_Py_Px_S_bc = I_ERI_H4xy_S_Px_S_bc+ABY*I_ERI_G4x_S_Px_S_bc;
  Double I_ERI_G3xy_Py_Px_S_bc = I_ERI_H3x2y_S_Px_S_bc+ABY*I_ERI_G3xy_S_Px_S_bc;
  Double I_ERI_G3xz_Py_Px_S_bc = I_ERI_H3xyz_S_Px_S_bc+ABY*I_ERI_G3xz_S_Px_S_bc;
  Double I_ERI_G2x2y_Py_Px_S_bc = I_ERI_H2x3y_S_Px_S_bc+ABY*I_ERI_G2x2y_S_Px_S_bc;
  Double I_ERI_G2xyz_Py_Px_S_bc = I_ERI_H2x2yz_S_Px_S_bc+ABY*I_ERI_G2xyz_S_Px_S_bc;
  Double I_ERI_G2x2z_Py_Px_S_bc = I_ERI_H2xy2z_S_Px_S_bc+ABY*I_ERI_G2x2z_S_Px_S_bc;
  Double I_ERI_Gx3y_Py_Px_S_bc = I_ERI_Hx4y_S_Px_S_bc+ABY*I_ERI_Gx3y_S_Px_S_bc;
  Double I_ERI_Gx2yz_Py_Px_S_bc = I_ERI_Hx3yz_S_Px_S_bc+ABY*I_ERI_Gx2yz_S_Px_S_bc;
  Double I_ERI_Gxy2z_Py_Px_S_bc = I_ERI_Hx2y2z_S_Px_S_bc+ABY*I_ERI_Gxy2z_S_Px_S_bc;
  Double I_ERI_Gx3z_Py_Px_S_bc = I_ERI_Hxy3z_S_Px_S_bc+ABY*I_ERI_Gx3z_S_Px_S_bc;
  Double I_ERI_G4y_Py_Px_S_bc = I_ERI_H5y_S_Px_S_bc+ABY*I_ERI_G4y_S_Px_S_bc;
  Double I_ERI_G3yz_Py_Px_S_bc = I_ERI_H4yz_S_Px_S_bc+ABY*I_ERI_G3yz_S_Px_S_bc;
  Double I_ERI_G2y2z_Py_Px_S_bc = I_ERI_H3y2z_S_Px_S_bc+ABY*I_ERI_G2y2z_S_Px_S_bc;
  Double I_ERI_Gy3z_Py_Px_S_bc = I_ERI_H2y3z_S_Px_S_bc+ABY*I_ERI_Gy3z_S_Px_S_bc;
  Double I_ERI_G4z_Py_Px_S_bc = I_ERI_Hy4z_S_Px_S_bc+ABY*I_ERI_G4z_S_Px_S_bc;
  Double I_ERI_G4x_Pz_Px_S_bc = I_ERI_H4xz_S_Px_S_bc+ABZ*I_ERI_G4x_S_Px_S_bc;
  Double I_ERI_G3xy_Pz_Px_S_bc = I_ERI_H3xyz_S_Px_S_bc+ABZ*I_ERI_G3xy_S_Px_S_bc;
  Double I_ERI_G3xz_Pz_Px_S_bc = I_ERI_H3x2z_S_Px_S_bc+ABZ*I_ERI_G3xz_S_Px_S_bc;
  Double I_ERI_G2x2y_Pz_Px_S_bc = I_ERI_H2x2yz_S_Px_S_bc+ABZ*I_ERI_G2x2y_S_Px_S_bc;
  Double I_ERI_G2xyz_Pz_Px_S_bc = I_ERI_H2xy2z_S_Px_S_bc+ABZ*I_ERI_G2xyz_S_Px_S_bc;
  Double I_ERI_G2x2z_Pz_Px_S_bc = I_ERI_H2x3z_S_Px_S_bc+ABZ*I_ERI_G2x2z_S_Px_S_bc;
  Double I_ERI_Gx3y_Pz_Px_S_bc = I_ERI_Hx3yz_S_Px_S_bc+ABZ*I_ERI_Gx3y_S_Px_S_bc;
  Double I_ERI_Gx2yz_Pz_Px_S_bc = I_ERI_Hx2y2z_S_Px_S_bc+ABZ*I_ERI_Gx2yz_S_Px_S_bc;
  Double I_ERI_Gxy2z_Pz_Px_S_bc = I_ERI_Hxy3z_S_Px_S_bc+ABZ*I_ERI_Gxy2z_S_Px_S_bc;
  Double I_ERI_Gx3z_Pz_Px_S_bc = I_ERI_Hx4z_S_Px_S_bc+ABZ*I_ERI_Gx3z_S_Px_S_bc;
  Double I_ERI_G4y_Pz_Px_S_bc = I_ERI_H4yz_S_Px_S_bc+ABZ*I_ERI_G4y_S_Px_S_bc;
  Double I_ERI_G3yz_Pz_Px_S_bc = I_ERI_H3y2z_S_Px_S_bc+ABZ*I_ERI_G3yz_S_Px_S_bc;
  Double I_ERI_G2y2z_Pz_Px_S_bc = I_ERI_H2y3z_S_Px_S_bc+ABZ*I_ERI_G2y2z_S_Px_S_bc;
  Double I_ERI_Gy3z_Pz_Px_S_bc = I_ERI_Hy4z_S_Px_S_bc+ABZ*I_ERI_Gy3z_S_Px_S_bc;
  Double I_ERI_G4z_Pz_Px_S_bc = I_ERI_H5z_S_Px_S_bc+ABZ*I_ERI_G4z_S_Px_S_bc;
  Double I_ERI_G4x_Px_Py_S_bc = I_ERI_H5x_S_Py_S_bc+ABX*I_ERI_G4x_S_Py_S_bc;
  Double I_ERI_G3xy_Px_Py_S_bc = I_ERI_H4xy_S_Py_S_bc+ABX*I_ERI_G3xy_S_Py_S_bc;
  Double I_ERI_G3xz_Px_Py_S_bc = I_ERI_H4xz_S_Py_S_bc+ABX*I_ERI_G3xz_S_Py_S_bc;
  Double I_ERI_G2x2y_Px_Py_S_bc = I_ERI_H3x2y_S_Py_S_bc+ABX*I_ERI_G2x2y_S_Py_S_bc;
  Double I_ERI_G2xyz_Px_Py_S_bc = I_ERI_H3xyz_S_Py_S_bc+ABX*I_ERI_G2xyz_S_Py_S_bc;
  Double I_ERI_G2x2z_Px_Py_S_bc = I_ERI_H3x2z_S_Py_S_bc+ABX*I_ERI_G2x2z_S_Py_S_bc;
  Double I_ERI_Gx3y_Px_Py_S_bc = I_ERI_H2x3y_S_Py_S_bc+ABX*I_ERI_Gx3y_S_Py_S_bc;
  Double I_ERI_Gx2yz_Px_Py_S_bc = I_ERI_H2x2yz_S_Py_S_bc+ABX*I_ERI_Gx2yz_S_Py_S_bc;
  Double I_ERI_Gxy2z_Px_Py_S_bc = I_ERI_H2xy2z_S_Py_S_bc+ABX*I_ERI_Gxy2z_S_Py_S_bc;
  Double I_ERI_Gx3z_Px_Py_S_bc = I_ERI_H2x3z_S_Py_S_bc+ABX*I_ERI_Gx3z_S_Py_S_bc;
  Double I_ERI_G4y_Px_Py_S_bc = I_ERI_Hx4y_S_Py_S_bc+ABX*I_ERI_G4y_S_Py_S_bc;
  Double I_ERI_G3yz_Px_Py_S_bc = I_ERI_Hx3yz_S_Py_S_bc+ABX*I_ERI_G3yz_S_Py_S_bc;
  Double I_ERI_G2y2z_Px_Py_S_bc = I_ERI_Hx2y2z_S_Py_S_bc+ABX*I_ERI_G2y2z_S_Py_S_bc;
  Double I_ERI_Gy3z_Px_Py_S_bc = I_ERI_Hxy3z_S_Py_S_bc+ABX*I_ERI_Gy3z_S_Py_S_bc;
  Double I_ERI_G4z_Px_Py_S_bc = I_ERI_Hx4z_S_Py_S_bc+ABX*I_ERI_G4z_S_Py_S_bc;
  Double I_ERI_G4x_Py_Py_S_bc = I_ERI_H4xy_S_Py_S_bc+ABY*I_ERI_G4x_S_Py_S_bc;
  Double I_ERI_G3xy_Py_Py_S_bc = I_ERI_H3x2y_S_Py_S_bc+ABY*I_ERI_G3xy_S_Py_S_bc;
  Double I_ERI_G3xz_Py_Py_S_bc = I_ERI_H3xyz_S_Py_S_bc+ABY*I_ERI_G3xz_S_Py_S_bc;
  Double I_ERI_G2x2y_Py_Py_S_bc = I_ERI_H2x3y_S_Py_S_bc+ABY*I_ERI_G2x2y_S_Py_S_bc;
  Double I_ERI_G2xyz_Py_Py_S_bc = I_ERI_H2x2yz_S_Py_S_bc+ABY*I_ERI_G2xyz_S_Py_S_bc;
  Double I_ERI_G2x2z_Py_Py_S_bc = I_ERI_H2xy2z_S_Py_S_bc+ABY*I_ERI_G2x2z_S_Py_S_bc;
  Double I_ERI_Gx3y_Py_Py_S_bc = I_ERI_Hx4y_S_Py_S_bc+ABY*I_ERI_Gx3y_S_Py_S_bc;
  Double I_ERI_Gx2yz_Py_Py_S_bc = I_ERI_Hx3yz_S_Py_S_bc+ABY*I_ERI_Gx2yz_S_Py_S_bc;
  Double I_ERI_Gxy2z_Py_Py_S_bc = I_ERI_Hx2y2z_S_Py_S_bc+ABY*I_ERI_Gxy2z_S_Py_S_bc;
  Double I_ERI_Gx3z_Py_Py_S_bc = I_ERI_Hxy3z_S_Py_S_bc+ABY*I_ERI_Gx3z_S_Py_S_bc;
  Double I_ERI_G4y_Py_Py_S_bc = I_ERI_H5y_S_Py_S_bc+ABY*I_ERI_G4y_S_Py_S_bc;
  Double I_ERI_G3yz_Py_Py_S_bc = I_ERI_H4yz_S_Py_S_bc+ABY*I_ERI_G3yz_S_Py_S_bc;
  Double I_ERI_G2y2z_Py_Py_S_bc = I_ERI_H3y2z_S_Py_S_bc+ABY*I_ERI_G2y2z_S_Py_S_bc;
  Double I_ERI_Gy3z_Py_Py_S_bc = I_ERI_H2y3z_S_Py_S_bc+ABY*I_ERI_Gy3z_S_Py_S_bc;
  Double I_ERI_G4z_Py_Py_S_bc = I_ERI_Hy4z_S_Py_S_bc+ABY*I_ERI_G4z_S_Py_S_bc;
  Double I_ERI_G4x_Pz_Py_S_bc = I_ERI_H4xz_S_Py_S_bc+ABZ*I_ERI_G4x_S_Py_S_bc;
  Double I_ERI_G3xy_Pz_Py_S_bc = I_ERI_H3xyz_S_Py_S_bc+ABZ*I_ERI_G3xy_S_Py_S_bc;
  Double I_ERI_G3xz_Pz_Py_S_bc = I_ERI_H3x2z_S_Py_S_bc+ABZ*I_ERI_G3xz_S_Py_S_bc;
  Double I_ERI_G2x2y_Pz_Py_S_bc = I_ERI_H2x2yz_S_Py_S_bc+ABZ*I_ERI_G2x2y_S_Py_S_bc;
  Double I_ERI_G2xyz_Pz_Py_S_bc = I_ERI_H2xy2z_S_Py_S_bc+ABZ*I_ERI_G2xyz_S_Py_S_bc;
  Double I_ERI_G2x2z_Pz_Py_S_bc = I_ERI_H2x3z_S_Py_S_bc+ABZ*I_ERI_G2x2z_S_Py_S_bc;
  Double I_ERI_Gx3y_Pz_Py_S_bc = I_ERI_Hx3yz_S_Py_S_bc+ABZ*I_ERI_Gx3y_S_Py_S_bc;
  Double I_ERI_Gx2yz_Pz_Py_S_bc = I_ERI_Hx2y2z_S_Py_S_bc+ABZ*I_ERI_Gx2yz_S_Py_S_bc;
  Double I_ERI_Gxy2z_Pz_Py_S_bc = I_ERI_Hxy3z_S_Py_S_bc+ABZ*I_ERI_Gxy2z_S_Py_S_bc;
  Double I_ERI_Gx3z_Pz_Py_S_bc = I_ERI_Hx4z_S_Py_S_bc+ABZ*I_ERI_Gx3z_S_Py_S_bc;
  Double I_ERI_G4y_Pz_Py_S_bc = I_ERI_H4yz_S_Py_S_bc+ABZ*I_ERI_G4y_S_Py_S_bc;
  Double I_ERI_G3yz_Pz_Py_S_bc = I_ERI_H3y2z_S_Py_S_bc+ABZ*I_ERI_G3yz_S_Py_S_bc;
  Double I_ERI_G2y2z_Pz_Py_S_bc = I_ERI_H2y3z_S_Py_S_bc+ABZ*I_ERI_G2y2z_S_Py_S_bc;
  Double I_ERI_Gy3z_Pz_Py_S_bc = I_ERI_Hy4z_S_Py_S_bc+ABZ*I_ERI_Gy3z_S_Py_S_bc;
  Double I_ERI_G4z_Pz_Py_S_bc = I_ERI_H5z_S_Py_S_bc+ABZ*I_ERI_G4z_S_Py_S_bc;
  Double I_ERI_G4x_Px_Pz_S_bc = I_ERI_H5x_S_Pz_S_bc+ABX*I_ERI_G4x_S_Pz_S_bc;
  Double I_ERI_G3xy_Px_Pz_S_bc = I_ERI_H4xy_S_Pz_S_bc+ABX*I_ERI_G3xy_S_Pz_S_bc;
  Double I_ERI_G3xz_Px_Pz_S_bc = I_ERI_H4xz_S_Pz_S_bc+ABX*I_ERI_G3xz_S_Pz_S_bc;
  Double I_ERI_G2x2y_Px_Pz_S_bc = I_ERI_H3x2y_S_Pz_S_bc+ABX*I_ERI_G2x2y_S_Pz_S_bc;
  Double I_ERI_G2xyz_Px_Pz_S_bc = I_ERI_H3xyz_S_Pz_S_bc+ABX*I_ERI_G2xyz_S_Pz_S_bc;
  Double I_ERI_G2x2z_Px_Pz_S_bc = I_ERI_H3x2z_S_Pz_S_bc+ABX*I_ERI_G2x2z_S_Pz_S_bc;
  Double I_ERI_Gx3y_Px_Pz_S_bc = I_ERI_H2x3y_S_Pz_S_bc+ABX*I_ERI_Gx3y_S_Pz_S_bc;
  Double I_ERI_Gx2yz_Px_Pz_S_bc = I_ERI_H2x2yz_S_Pz_S_bc+ABX*I_ERI_Gx2yz_S_Pz_S_bc;
  Double I_ERI_Gxy2z_Px_Pz_S_bc = I_ERI_H2xy2z_S_Pz_S_bc+ABX*I_ERI_Gxy2z_S_Pz_S_bc;
  Double I_ERI_Gx3z_Px_Pz_S_bc = I_ERI_H2x3z_S_Pz_S_bc+ABX*I_ERI_Gx3z_S_Pz_S_bc;
  Double I_ERI_G4y_Px_Pz_S_bc = I_ERI_Hx4y_S_Pz_S_bc+ABX*I_ERI_G4y_S_Pz_S_bc;
  Double I_ERI_G3yz_Px_Pz_S_bc = I_ERI_Hx3yz_S_Pz_S_bc+ABX*I_ERI_G3yz_S_Pz_S_bc;
  Double I_ERI_G2y2z_Px_Pz_S_bc = I_ERI_Hx2y2z_S_Pz_S_bc+ABX*I_ERI_G2y2z_S_Pz_S_bc;
  Double I_ERI_Gy3z_Px_Pz_S_bc = I_ERI_Hxy3z_S_Pz_S_bc+ABX*I_ERI_Gy3z_S_Pz_S_bc;
  Double I_ERI_G4z_Px_Pz_S_bc = I_ERI_Hx4z_S_Pz_S_bc+ABX*I_ERI_G4z_S_Pz_S_bc;
  Double I_ERI_G4x_Py_Pz_S_bc = I_ERI_H4xy_S_Pz_S_bc+ABY*I_ERI_G4x_S_Pz_S_bc;
  Double I_ERI_G3xy_Py_Pz_S_bc = I_ERI_H3x2y_S_Pz_S_bc+ABY*I_ERI_G3xy_S_Pz_S_bc;
  Double I_ERI_G3xz_Py_Pz_S_bc = I_ERI_H3xyz_S_Pz_S_bc+ABY*I_ERI_G3xz_S_Pz_S_bc;
  Double I_ERI_G2x2y_Py_Pz_S_bc = I_ERI_H2x3y_S_Pz_S_bc+ABY*I_ERI_G2x2y_S_Pz_S_bc;
  Double I_ERI_G2xyz_Py_Pz_S_bc = I_ERI_H2x2yz_S_Pz_S_bc+ABY*I_ERI_G2xyz_S_Pz_S_bc;
  Double I_ERI_G2x2z_Py_Pz_S_bc = I_ERI_H2xy2z_S_Pz_S_bc+ABY*I_ERI_G2x2z_S_Pz_S_bc;
  Double I_ERI_Gx3y_Py_Pz_S_bc = I_ERI_Hx4y_S_Pz_S_bc+ABY*I_ERI_Gx3y_S_Pz_S_bc;
  Double I_ERI_Gx2yz_Py_Pz_S_bc = I_ERI_Hx3yz_S_Pz_S_bc+ABY*I_ERI_Gx2yz_S_Pz_S_bc;
  Double I_ERI_Gxy2z_Py_Pz_S_bc = I_ERI_Hx2y2z_S_Pz_S_bc+ABY*I_ERI_Gxy2z_S_Pz_S_bc;
  Double I_ERI_Gx3z_Py_Pz_S_bc = I_ERI_Hxy3z_S_Pz_S_bc+ABY*I_ERI_Gx3z_S_Pz_S_bc;
  Double I_ERI_G4y_Py_Pz_S_bc = I_ERI_H5y_S_Pz_S_bc+ABY*I_ERI_G4y_S_Pz_S_bc;
  Double I_ERI_G3yz_Py_Pz_S_bc = I_ERI_H4yz_S_Pz_S_bc+ABY*I_ERI_G3yz_S_Pz_S_bc;
  Double I_ERI_G2y2z_Py_Pz_S_bc = I_ERI_H3y2z_S_Pz_S_bc+ABY*I_ERI_G2y2z_S_Pz_S_bc;
  Double I_ERI_Gy3z_Py_Pz_S_bc = I_ERI_H2y3z_S_Pz_S_bc+ABY*I_ERI_Gy3z_S_Pz_S_bc;
  Double I_ERI_G4z_Py_Pz_S_bc = I_ERI_Hy4z_S_Pz_S_bc+ABY*I_ERI_G4z_S_Pz_S_bc;
  Double I_ERI_G4x_Pz_Pz_S_bc = I_ERI_H4xz_S_Pz_S_bc+ABZ*I_ERI_G4x_S_Pz_S_bc;
  Double I_ERI_G3xy_Pz_Pz_S_bc = I_ERI_H3xyz_S_Pz_S_bc+ABZ*I_ERI_G3xy_S_Pz_S_bc;
  Double I_ERI_G3xz_Pz_Pz_S_bc = I_ERI_H3x2z_S_Pz_S_bc+ABZ*I_ERI_G3xz_S_Pz_S_bc;
  Double I_ERI_G2x2y_Pz_Pz_S_bc = I_ERI_H2x2yz_S_Pz_S_bc+ABZ*I_ERI_G2x2y_S_Pz_S_bc;
  Double I_ERI_G2xyz_Pz_Pz_S_bc = I_ERI_H2xy2z_S_Pz_S_bc+ABZ*I_ERI_G2xyz_S_Pz_S_bc;
  Double I_ERI_G2x2z_Pz_Pz_S_bc = I_ERI_H2x3z_S_Pz_S_bc+ABZ*I_ERI_G2x2z_S_Pz_S_bc;
  Double I_ERI_Gx3y_Pz_Pz_S_bc = I_ERI_Hx3yz_S_Pz_S_bc+ABZ*I_ERI_Gx3y_S_Pz_S_bc;
  Double I_ERI_Gx2yz_Pz_Pz_S_bc = I_ERI_Hx2y2z_S_Pz_S_bc+ABZ*I_ERI_Gx2yz_S_Pz_S_bc;
  Double I_ERI_Gxy2z_Pz_Pz_S_bc = I_ERI_Hxy3z_S_Pz_S_bc+ABZ*I_ERI_Gxy2z_S_Pz_S_bc;
  Double I_ERI_Gx3z_Pz_Pz_S_bc = I_ERI_Hx4z_S_Pz_S_bc+ABZ*I_ERI_Gx3z_S_Pz_S_bc;
  Double I_ERI_G4y_Pz_Pz_S_bc = I_ERI_H4yz_S_Pz_S_bc+ABZ*I_ERI_G4y_S_Pz_S_bc;
  Double I_ERI_G3yz_Pz_Pz_S_bc = I_ERI_H3y2z_S_Pz_S_bc+ABZ*I_ERI_G3yz_S_Pz_S_bc;
  Double I_ERI_G2y2z_Pz_Pz_S_bc = I_ERI_H2y3z_S_Pz_S_bc+ABZ*I_ERI_G2y2z_S_Pz_S_bc;
  Double I_ERI_Gy3z_Pz_Pz_S_bc = I_ERI_Hy4z_S_Pz_S_bc+ABZ*I_ERI_Gy3z_S_Pz_S_bc;
  Double I_ERI_G4z_Pz_Pz_S_bc = I_ERI_H5z_S_Pz_S_bc+ABZ*I_ERI_G4z_S_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_P_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 60 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_bc
   ************************************************************/
  Double I_ERI_F3x_D2x_Px_S_bc = I_ERI_G4x_Px_Px_S_bc+ABX*I_ERI_F3x_Px_Px_S_bc;
  Double I_ERI_F2xy_D2x_Px_S_bc = I_ERI_G3xy_Px_Px_S_bc+ABX*I_ERI_F2xy_Px_Px_S_bc;
  Double I_ERI_F2xz_D2x_Px_S_bc = I_ERI_G3xz_Px_Px_S_bc+ABX*I_ERI_F2xz_Px_Px_S_bc;
  Double I_ERI_Fx2y_D2x_Px_S_bc = I_ERI_G2x2y_Px_Px_S_bc+ABX*I_ERI_Fx2y_Px_Px_S_bc;
  Double I_ERI_Fxyz_D2x_Px_S_bc = I_ERI_G2xyz_Px_Px_S_bc+ABX*I_ERI_Fxyz_Px_Px_S_bc;
  Double I_ERI_Fx2z_D2x_Px_S_bc = I_ERI_G2x2z_Px_Px_S_bc+ABX*I_ERI_Fx2z_Px_Px_S_bc;
  Double I_ERI_F3y_D2x_Px_S_bc = I_ERI_Gx3y_Px_Px_S_bc+ABX*I_ERI_F3y_Px_Px_S_bc;
  Double I_ERI_F2yz_D2x_Px_S_bc = I_ERI_Gx2yz_Px_Px_S_bc+ABX*I_ERI_F2yz_Px_Px_S_bc;
  Double I_ERI_Fy2z_D2x_Px_S_bc = I_ERI_Gxy2z_Px_Px_S_bc+ABX*I_ERI_Fy2z_Px_Px_S_bc;
  Double I_ERI_F3z_D2x_Px_S_bc = I_ERI_Gx3z_Px_Px_S_bc+ABX*I_ERI_F3z_Px_Px_S_bc;
  Double I_ERI_F3x_Dxy_Px_S_bc = I_ERI_G3xy_Px_Px_S_bc+ABY*I_ERI_F3x_Px_Px_S_bc;
  Double I_ERI_F2xy_Dxy_Px_S_bc = I_ERI_G2x2y_Px_Px_S_bc+ABY*I_ERI_F2xy_Px_Px_S_bc;
  Double I_ERI_F2xz_Dxy_Px_S_bc = I_ERI_G2xyz_Px_Px_S_bc+ABY*I_ERI_F2xz_Px_Px_S_bc;
  Double I_ERI_Fx2y_Dxy_Px_S_bc = I_ERI_Gx3y_Px_Px_S_bc+ABY*I_ERI_Fx2y_Px_Px_S_bc;
  Double I_ERI_Fxyz_Dxy_Px_S_bc = I_ERI_Gx2yz_Px_Px_S_bc+ABY*I_ERI_Fxyz_Px_Px_S_bc;
  Double I_ERI_Fx2z_Dxy_Px_S_bc = I_ERI_Gxy2z_Px_Px_S_bc+ABY*I_ERI_Fx2z_Px_Px_S_bc;
  Double I_ERI_F3y_Dxy_Px_S_bc = I_ERI_G4y_Px_Px_S_bc+ABY*I_ERI_F3y_Px_Px_S_bc;
  Double I_ERI_F2yz_Dxy_Px_S_bc = I_ERI_G3yz_Px_Px_S_bc+ABY*I_ERI_F2yz_Px_Px_S_bc;
  Double I_ERI_Fy2z_Dxy_Px_S_bc = I_ERI_G2y2z_Px_Px_S_bc+ABY*I_ERI_Fy2z_Px_Px_S_bc;
  Double I_ERI_F3z_Dxy_Px_S_bc = I_ERI_Gy3z_Px_Px_S_bc+ABY*I_ERI_F3z_Px_Px_S_bc;
  Double I_ERI_F3x_D2y_Px_S_bc = I_ERI_G3xy_Py_Px_S_bc+ABY*I_ERI_F3x_Py_Px_S_bc;
  Double I_ERI_F2xy_D2y_Px_S_bc = I_ERI_G2x2y_Py_Px_S_bc+ABY*I_ERI_F2xy_Py_Px_S_bc;
  Double I_ERI_F2xz_D2y_Px_S_bc = I_ERI_G2xyz_Py_Px_S_bc+ABY*I_ERI_F2xz_Py_Px_S_bc;
  Double I_ERI_Fx2y_D2y_Px_S_bc = I_ERI_Gx3y_Py_Px_S_bc+ABY*I_ERI_Fx2y_Py_Px_S_bc;
  Double I_ERI_Fxyz_D2y_Px_S_bc = I_ERI_Gx2yz_Py_Px_S_bc+ABY*I_ERI_Fxyz_Py_Px_S_bc;
  Double I_ERI_Fx2z_D2y_Px_S_bc = I_ERI_Gxy2z_Py_Px_S_bc+ABY*I_ERI_Fx2z_Py_Px_S_bc;
  Double I_ERI_F3y_D2y_Px_S_bc = I_ERI_G4y_Py_Px_S_bc+ABY*I_ERI_F3y_Py_Px_S_bc;
  Double I_ERI_F2yz_D2y_Px_S_bc = I_ERI_G3yz_Py_Px_S_bc+ABY*I_ERI_F2yz_Py_Px_S_bc;
  Double I_ERI_Fy2z_D2y_Px_S_bc = I_ERI_G2y2z_Py_Px_S_bc+ABY*I_ERI_Fy2z_Py_Px_S_bc;
  Double I_ERI_F3z_D2y_Px_S_bc = I_ERI_Gy3z_Py_Px_S_bc+ABY*I_ERI_F3z_Py_Px_S_bc;
  Double I_ERI_F3x_D2z_Px_S_bc = I_ERI_G3xz_Pz_Px_S_bc+ABZ*I_ERI_F3x_Pz_Px_S_bc;
  Double I_ERI_F2xy_D2z_Px_S_bc = I_ERI_G2xyz_Pz_Px_S_bc+ABZ*I_ERI_F2xy_Pz_Px_S_bc;
  Double I_ERI_F2xz_D2z_Px_S_bc = I_ERI_G2x2z_Pz_Px_S_bc+ABZ*I_ERI_F2xz_Pz_Px_S_bc;
  Double I_ERI_Fx2y_D2z_Px_S_bc = I_ERI_Gx2yz_Pz_Px_S_bc+ABZ*I_ERI_Fx2y_Pz_Px_S_bc;
  Double I_ERI_Fxyz_D2z_Px_S_bc = I_ERI_Gxy2z_Pz_Px_S_bc+ABZ*I_ERI_Fxyz_Pz_Px_S_bc;
  Double I_ERI_Fx2z_D2z_Px_S_bc = I_ERI_Gx3z_Pz_Px_S_bc+ABZ*I_ERI_Fx2z_Pz_Px_S_bc;
  Double I_ERI_F3y_D2z_Px_S_bc = I_ERI_G3yz_Pz_Px_S_bc+ABZ*I_ERI_F3y_Pz_Px_S_bc;
  Double I_ERI_F2yz_D2z_Px_S_bc = I_ERI_G2y2z_Pz_Px_S_bc+ABZ*I_ERI_F2yz_Pz_Px_S_bc;
  Double I_ERI_Fy2z_D2z_Px_S_bc = I_ERI_Gy3z_Pz_Px_S_bc+ABZ*I_ERI_Fy2z_Pz_Px_S_bc;
  Double I_ERI_F3z_D2z_Px_S_bc = I_ERI_G4z_Pz_Px_S_bc+ABZ*I_ERI_F3z_Pz_Px_S_bc;
  Double I_ERI_F3x_D2x_Py_S_bc = I_ERI_G4x_Px_Py_S_bc+ABX*I_ERI_F3x_Px_Py_S_bc;
  Double I_ERI_F2xy_D2x_Py_S_bc = I_ERI_G3xy_Px_Py_S_bc+ABX*I_ERI_F2xy_Px_Py_S_bc;
  Double I_ERI_F2xz_D2x_Py_S_bc = I_ERI_G3xz_Px_Py_S_bc+ABX*I_ERI_F2xz_Px_Py_S_bc;
  Double I_ERI_Fx2y_D2x_Py_S_bc = I_ERI_G2x2y_Px_Py_S_bc+ABX*I_ERI_Fx2y_Px_Py_S_bc;
  Double I_ERI_Fxyz_D2x_Py_S_bc = I_ERI_G2xyz_Px_Py_S_bc+ABX*I_ERI_Fxyz_Px_Py_S_bc;
  Double I_ERI_Fx2z_D2x_Py_S_bc = I_ERI_G2x2z_Px_Py_S_bc+ABX*I_ERI_Fx2z_Px_Py_S_bc;
  Double I_ERI_F3y_D2x_Py_S_bc = I_ERI_Gx3y_Px_Py_S_bc+ABX*I_ERI_F3y_Px_Py_S_bc;
  Double I_ERI_F2yz_D2x_Py_S_bc = I_ERI_Gx2yz_Px_Py_S_bc+ABX*I_ERI_F2yz_Px_Py_S_bc;
  Double I_ERI_Fy2z_D2x_Py_S_bc = I_ERI_Gxy2z_Px_Py_S_bc+ABX*I_ERI_Fy2z_Px_Py_S_bc;
  Double I_ERI_F3z_D2x_Py_S_bc = I_ERI_Gx3z_Px_Py_S_bc+ABX*I_ERI_F3z_Px_Py_S_bc;
  Double I_ERI_F3x_Dxy_Py_S_bc = I_ERI_G3xy_Px_Py_S_bc+ABY*I_ERI_F3x_Px_Py_S_bc;
  Double I_ERI_F2xy_Dxy_Py_S_bc = I_ERI_G2x2y_Px_Py_S_bc+ABY*I_ERI_F2xy_Px_Py_S_bc;
  Double I_ERI_F2xz_Dxy_Py_S_bc = I_ERI_G2xyz_Px_Py_S_bc+ABY*I_ERI_F2xz_Px_Py_S_bc;
  Double I_ERI_Fx2y_Dxy_Py_S_bc = I_ERI_Gx3y_Px_Py_S_bc+ABY*I_ERI_Fx2y_Px_Py_S_bc;
  Double I_ERI_Fxyz_Dxy_Py_S_bc = I_ERI_Gx2yz_Px_Py_S_bc+ABY*I_ERI_Fxyz_Px_Py_S_bc;
  Double I_ERI_Fx2z_Dxy_Py_S_bc = I_ERI_Gxy2z_Px_Py_S_bc+ABY*I_ERI_Fx2z_Px_Py_S_bc;
  Double I_ERI_F3y_Dxy_Py_S_bc = I_ERI_G4y_Px_Py_S_bc+ABY*I_ERI_F3y_Px_Py_S_bc;
  Double I_ERI_F2yz_Dxy_Py_S_bc = I_ERI_G3yz_Px_Py_S_bc+ABY*I_ERI_F2yz_Px_Py_S_bc;
  Double I_ERI_Fy2z_Dxy_Py_S_bc = I_ERI_G2y2z_Px_Py_S_bc+ABY*I_ERI_Fy2z_Px_Py_S_bc;
  Double I_ERI_F3z_Dxy_Py_S_bc = I_ERI_Gy3z_Px_Py_S_bc+ABY*I_ERI_F3z_Px_Py_S_bc;
  Double I_ERI_F3x_D2y_Py_S_bc = I_ERI_G3xy_Py_Py_S_bc+ABY*I_ERI_F3x_Py_Py_S_bc;
  Double I_ERI_F2xy_D2y_Py_S_bc = I_ERI_G2x2y_Py_Py_S_bc+ABY*I_ERI_F2xy_Py_Py_S_bc;
  Double I_ERI_F2xz_D2y_Py_S_bc = I_ERI_G2xyz_Py_Py_S_bc+ABY*I_ERI_F2xz_Py_Py_S_bc;
  Double I_ERI_Fx2y_D2y_Py_S_bc = I_ERI_Gx3y_Py_Py_S_bc+ABY*I_ERI_Fx2y_Py_Py_S_bc;
  Double I_ERI_Fxyz_D2y_Py_S_bc = I_ERI_Gx2yz_Py_Py_S_bc+ABY*I_ERI_Fxyz_Py_Py_S_bc;
  Double I_ERI_Fx2z_D2y_Py_S_bc = I_ERI_Gxy2z_Py_Py_S_bc+ABY*I_ERI_Fx2z_Py_Py_S_bc;
  Double I_ERI_F3y_D2y_Py_S_bc = I_ERI_G4y_Py_Py_S_bc+ABY*I_ERI_F3y_Py_Py_S_bc;
  Double I_ERI_F2yz_D2y_Py_S_bc = I_ERI_G3yz_Py_Py_S_bc+ABY*I_ERI_F2yz_Py_Py_S_bc;
  Double I_ERI_Fy2z_D2y_Py_S_bc = I_ERI_G2y2z_Py_Py_S_bc+ABY*I_ERI_Fy2z_Py_Py_S_bc;
  Double I_ERI_F3z_D2y_Py_S_bc = I_ERI_Gy3z_Py_Py_S_bc+ABY*I_ERI_F3z_Py_Py_S_bc;
  Double I_ERI_F3x_D2z_Py_S_bc = I_ERI_G3xz_Pz_Py_S_bc+ABZ*I_ERI_F3x_Pz_Py_S_bc;
  Double I_ERI_F2xy_D2z_Py_S_bc = I_ERI_G2xyz_Pz_Py_S_bc+ABZ*I_ERI_F2xy_Pz_Py_S_bc;
  Double I_ERI_F2xz_D2z_Py_S_bc = I_ERI_G2x2z_Pz_Py_S_bc+ABZ*I_ERI_F2xz_Pz_Py_S_bc;
  Double I_ERI_Fx2y_D2z_Py_S_bc = I_ERI_Gx2yz_Pz_Py_S_bc+ABZ*I_ERI_Fx2y_Pz_Py_S_bc;
  Double I_ERI_Fxyz_D2z_Py_S_bc = I_ERI_Gxy2z_Pz_Py_S_bc+ABZ*I_ERI_Fxyz_Pz_Py_S_bc;
  Double I_ERI_Fx2z_D2z_Py_S_bc = I_ERI_Gx3z_Pz_Py_S_bc+ABZ*I_ERI_Fx2z_Pz_Py_S_bc;
  Double I_ERI_F3y_D2z_Py_S_bc = I_ERI_G3yz_Pz_Py_S_bc+ABZ*I_ERI_F3y_Pz_Py_S_bc;
  Double I_ERI_F2yz_D2z_Py_S_bc = I_ERI_G2y2z_Pz_Py_S_bc+ABZ*I_ERI_F2yz_Pz_Py_S_bc;
  Double I_ERI_Fy2z_D2z_Py_S_bc = I_ERI_Gy3z_Pz_Py_S_bc+ABZ*I_ERI_Fy2z_Pz_Py_S_bc;
  Double I_ERI_F3z_D2z_Py_S_bc = I_ERI_G4z_Pz_Py_S_bc+ABZ*I_ERI_F3z_Pz_Py_S_bc;
  Double I_ERI_F3x_D2x_Pz_S_bc = I_ERI_G4x_Px_Pz_S_bc+ABX*I_ERI_F3x_Px_Pz_S_bc;
  Double I_ERI_F2xy_D2x_Pz_S_bc = I_ERI_G3xy_Px_Pz_S_bc+ABX*I_ERI_F2xy_Px_Pz_S_bc;
  Double I_ERI_F2xz_D2x_Pz_S_bc = I_ERI_G3xz_Px_Pz_S_bc+ABX*I_ERI_F2xz_Px_Pz_S_bc;
  Double I_ERI_Fx2y_D2x_Pz_S_bc = I_ERI_G2x2y_Px_Pz_S_bc+ABX*I_ERI_Fx2y_Px_Pz_S_bc;
  Double I_ERI_Fxyz_D2x_Pz_S_bc = I_ERI_G2xyz_Px_Pz_S_bc+ABX*I_ERI_Fxyz_Px_Pz_S_bc;
  Double I_ERI_Fx2z_D2x_Pz_S_bc = I_ERI_G2x2z_Px_Pz_S_bc+ABX*I_ERI_Fx2z_Px_Pz_S_bc;
  Double I_ERI_F3y_D2x_Pz_S_bc = I_ERI_Gx3y_Px_Pz_S_bc+ABX*I_ERI_F3y_Px_Pz_S_bc;
  Double I_ERI_F2yz_D2x_Pz_S_bc = I_ERI_Gx2yz_Px_Pz_S_bc+ABX*I_ERI_F2yz_Px_Pz_S_bc;
  Double I_ERI_Fy2z_D2x_Pz_S_bc = I_ERI_Gxy2z_Px_Pz_S_bc+ABX*I_ERI_Fy2z_Px_Pz_S_bc;
  Double I_ERI_F3z_D2x_Pz_S_bc = I_ERI_Gx3z_Px_Pz_S_bc+ABX*I_ERI_F3z_Px_Pz_S_bc;
  Double I_ERI_F3x_Dxy_Pz_S_bc = I_ERI_G3xy_Px_Pz_S_bc+ABY*I_ERI_F3x_Px_Pz_S_bc;
  Double I_ERI_F2xy_Dxy_Pz_S_bc = I_ERI_G2x2y_Px_Pz_S_bc+ABY*I_ERI_F2xy_Px_Pz_S_bc;
  Double I_ERI_F2xz_Dxy_Pz_S_bc = I_ERI_G2xyz_Px_Pz_S_bc+ABY*I_ERI_F2xz_Px_Pz_S_bc;
  Double I_ERI_Fx2y_Dxy_Pz_S_bc = I_ERI_Gx3y_Px_Pz_S_bc+ABY*I_ERI_Fx2y_Px_Pz_S_bc;
  Double I_ERI_Fxyz_Dxy_Pz_S_bc = I_ERI_Gx2yz_Px_Pz_S_bc+ABY*I_ERI_Fxyz_Px_Pz_S_bc;
  Double I_ERI_Fx2z_Dxy_Pz_S_bc = I_ERI_Gxy2z_Px_Pz_S_bc+ABY*I_ERI_Fx2z_Px_Pz_S_bc;
  Double I_ERI_F3y_Dxy_Pz_S_bc = I_ERI_G4y_Px_Pz_S_bc+ABY*I_ERI_F3y_Px_Pz_S_bc;
  Double I_ERI_F2yz_Dxy_Pz_S_bc = I_ERI_G3yz_Px_Pz_S_bc+ABY*I_ERI_F2yz_Px_Pz_S_bc;
  Double I_ERI_Fy2z_Dxy_Pz_S_bc = I_ERI_G2y2z_Px_Pz_S_bc+ABY*I_ERI_Fy2z_Px_Pz_S_bc;
  Double I_ERI_F3z_Dxy_Pz_S_bc = I_ERI_Gy3z_Px_Pz_S_bc+ABY*I_ERI_F3z_Px_Pz_S_bc;
  Double I_ERI_F3x_D2y_Pz_S_bc = I_ERI_G3xy_Py_Pz_S_bc+ABY*I_ERI_F3x_Py_Pz_S_bc;
  Double I_ERI_F2xy_D2y_Pz_S_bc = I_ERI_G2x2y_Py_Pz_S_bc+ABY*I_ERI_F2xy_Py_Pz_S_bc;
  Double I_ERI_F2xz_D2y_Pz_S_bc = I_ERI_G2xyz_Py_Pz_S_bc+ABY*I_ERI_F2xz_Py_Pz_S_bc;
  Double I_ERI_Fx2y_D2y_Pz_S_bc = I_ERI_Gx3y_Py_Pz_S_bc+ABY*I_ERI_Fx2y_Py_Pz_S_bc;
  Double I_ERI_Fxyz_D2y_Pz_S_bc = I_ERI_Gx2yz_Py_Pz_S_bc+ABY*I_ERI_Fxyz_Py_Pz_S_bc;
  Double I_ERI_Fx2z_D2y_Pz_S_bc = I_ERI_Gxy2z_Py_Pz_S_bc+ABY*I_ERI_Fx2z_Py_Pz_S_bc;
  Double I_ERI_F3y_D2y_Pz_S_bc = I_ERI_G4y_Py_Pz_S_bc+ABY*I_ERI_F3y_Py_Pz_S_bc;
  Double I_ERI_F2yz_D2y_Pz_S_bc = I_ERI_G3yz_Py_Pz_S_bc+ABY*I_ERI_F2yz_Py_Pz_S_bc;
  Double I_ERI_Fy2z_D2y_Pz_S_bc = I_ERI_G2y2z_Py_Pz_S_bc+ABY*I_ERI_Fy2z_Py_Pz_S_bc;
  Double I_ERI_F3z_D2y_Pz_S_bc = I_ERI_Gy3z_Py_Pz_S_bc+ABY*I_ERI_F3z_Py_Pz_S_bc;
  Double I_ERI_F3x_D2z_Pz_S_bc = I_ERI_G3xz_Pz_Pz_S_bc+ABZ*I_ERI_F3x_Pz_Pz_S_bc;
  Double I_ERI_F2xy_D2z_Pz_S_bc = I_ERI_G2xyz_Pz_Pz_S_bc+ABZ*I_ERI_F2xy_Pz_Pz_S_bc;
  Double I_ERI_F2xz_D2z_Pz_S_bc = I_ERI_G2x2z_Pz_Pz_S_bc+ABZ*I_ERI_F2xz_Pz_Pz_S_bc;
  Double I_ERI_Fx2y_D2z_Pz_S_bc = I_ERI_Gx2yz_Pz_Pz_S_bc+ABZ*I_ERI_Fx2y_Pz_Pz_S_bc;
  Double I_ERI_Fxyz_D2z_Pz_S_bc = I_ERI_Gxy2z_Pz_Pz_S_bc+ABZ*I_ERI_Fxyz_Pz_Pz_S_bc;
  Double I_ERI_Fx2z_D2z_Pz_S_bc = I_ERI_Gx3z_Pz_Pz_S_bc+ABZ*I_ERI_Fx2z_Pz_Pz_S_bc;
  Double I_ERI_F3y_D2z_Pz_S_bc = I_ERI_G3yz_Pz_Pz_S_bc+ABZ*I_ERI_F3y_Pz_Pz_S_bc;
  Double I_ERI_F2yz_D2z_Pz_S_bc = I_ERI_G2y2z_Pz_Pz_S_bc+ABZ*I_ERI_F2yz_Pz_Pz_S_bc;
  Double I_ERI_Fy2z_D2z_Pz_S_bc = I_ERI_Gy3z_Pz_Pz_S_bc+ABZ*I_ERI_Fy2z_Pz_Pz_S_bc;
  Double I_ERI_F3z_D2z_Pz_S_bc = I_ERI_G4z_Pz_Pz_S_bc+ABZ*I_ERI_F3z_Pz_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_H_P_P_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 42 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_P_S_bc
   * RHS shell quartet name: SQ_ERI_H_S_P_S_bc
   ************************************************************/
  Double I_ERI_H5x_Px_Px_S_bc = I_ERI_I6x_S_Px_S_bc+ABX*I_ERI_H5x_S_Px_S_bc;
  Double I_ERI_H4xy_Px_Px_S_bc = I_ERI_I5xy_S_Px_S_bc+ABX*I_ERI_H4xy_S_Px_S_bc;
  Double I_ERI_H4xz_Px_Px_S_bc = I_ERI_I5xz_S_Px_S_bc+ABX*I_ERI_H4xz_S_Px_S_bc;
  Double I_ERI_H3x2y_Px_Px_S_bc = I_ERI_I4x2y_S_Px_S_bc+ABX*I_ERI_H3x2y_S_Px_S_bc;
  Double I_ERI_H3xyz_Px_Px_S_bc = I_ERI_I4xyz_S_Px_S_bc+ABX*I_ERI_H3xyz_S_Px_S_bc;
  Double I_ERI_H3x2z_Px_Px_S_bc = I_ERI_I4x2z_S_Px_S_bc+ABX*I_ERI_H3x2z_S_Px_S_bc;
  Double I_ERI_H2x3y_Px_Px_S_bc = I_ERI_I3x3y_S_Px_S_bc+ABX*I_ERI_H2x3y_S_Px_S_bc;
  Double I_ERI_H2x2yz_Px_Px_S_bc = I_ERI_I3x2yz_S_Px_S_bc+ABX*I_ERI_H2x2yz_S_Px_S_bc;
  Double I_ERI_H2xy2z_Px_Px_S_bc = I_ERI_I3xy2z_S_Px_S_bc+ABX*I_ERI_H2xy2z_S_Px_S_bc;
  Double I_ERI_H2x3z_Px_Px_S_bc = I_ERI_I3x3z_S_Px_S_bc+ABX*I_ERI_H2x3z_S_Px_S_bc;
  Double I_ERI_Hx4y_Px_Px_S_bc = I_ERI_I2x4y_S_Px_S_bc+ABX*I_ERI_Hx4y_S_Px_S_bc;
  Double I_ERI_Hx3yz_Px_Px_S_bc = I_ERI_I2x3yz_S_Px_S_bc+ABX*I_ERI_Hx3yz_S_Px_S_bc;
  Double I_ERI_Hx2y2z_Px_Px_S_bc = I_ERI_I2x2y2z_S_Px_S_bc+ABX*I_ERI_Hx2y2z_S_Px_S_bc;
  Double I_ERI_Hxy3z_Px_Px_S_bc = I_ERI_I2xy3z_S_Px_S_bc+ABX*I_ERI_Hxy3z_S_Px_S_bc;
  Double I_ERI_Hx4z_Px_Px_S_bc = I_ERI_I2x4z_S_Px_S_bc+ABX*I_ERI_Hx4z_S_Px_S_bc;
  Double I_ERI_H4yz_Px_Px_S_bc = I_ERI_Ix4yz_S_Px_S_bc+ABX*I_ERI_H4yz_S_Px_S_bc;
  Double I_ERI_H3y2z_Px_Px_S_bc = I_ERI_Ix3y2z_S_Px_S_bc+ABX*I_ERI_H3y2z_S_Px_S_bc;
  Double I_ERI_H2y3z_Px_Px_S_bc = I_ERI_Ix2y3z_S_Px_S_bc+ABX*I_ERI_H2y3z_S_Px_S_bc;
  Double I_ERI_Hy4z_Px_Px_S_bc = I_ERI_Ixy4z_S_Px_S_bc+ABX*I_ERI_Hy4z_S_Px_S_bc;
  Double I_ERI_H4xy_Py_Px_S_bc = I_ERI_I4x2y_S_Px_S_bc+ABY*I_ERI_H4xy_S_Px_S_bc;
  Double I_ERI_H3x2y_Py_Px_S_bc = I_ERI_I3x3y_S_Px_S_bc+ABY*I_ERI_H3x2y_S_Px_S_bc;
  Double I_ERI_H3xyz_Py_Px_S_bc = I_ERI_I3x2yz_S_Px_S_bc+ABY*I_ERI_H3xyz_S_Px_S_bc;
  Double I_ERI_H2x3y_Py_Px_S_bc = I_ERI_I2x4y_S_Px_S_bc+ABY*I_ERI_H2x3y_S_Px_S_bc;
  Double I_ERI_H2x2yz_Py_Px_S_bc = I_ERI_I2x3yz_S_Px_S_bc+ABY*I_ERI_H2x2yz_S_Px_S_bc;
  Double I_ERI_H2xy2z_Py_Px_S_bc = I_ERI_I2x2y2z_S_Px_S_bc+ABY*I_ERI_H2xy2z_S_Px_S_bc;
  Double I_ERI_Hx4y_Py_Px_S_bc = I_ERI_Ix5y_S_Px_S_bc+ABY*I_ERI_Hx4y_S_Px_S_bc;
  Double I_ERI_Hx3yz_Py_Px_S_bc = I_ERI_Ix4yz_S_Px_S_bc+ABY*I_ERI_Hx3yz_S_Px_S_bc;
  Double I_ERI_Hx2y2z_Py_Px_S_bc = I_ERI_Ix3y2z_S_Px_S_bc+ABY*I_ERI_Hx2y2z_S_Px_S_bc;
  Double I_ERI_Hxy3z_Py_Px_S_bc = I_ERI_Ix2y3z_S_Px_S_bc+ABY*I_ERI_Hxy3z_S_Px_S_bc;
  Double I_ERI_H5y_Py_Px_S_bc = I_ERI_I6y_S_Px_S_bc+ABY*I_ERI_H5y_S_Px_S_bc;
  Double I_ERI_H4yz_Py_Px_S_bc = I_ERI_I5yz_S_Px_S_bc+ABY*I_ERI_H4yz_S_Px_S_bc;
  Double I_ERI_H3y2z_Py_Px_S_bc = I_ERI_I4y2z_S_Px_S_bc+ABY*I_ERI_H3y2z_S_Px_S_bc;
  Double I_ERI_H2y3z_Py_Px_S_bc = I_ERI_I3y3z_S_Px_S_bc+ABY*I_ERI_H2y3z_S_Px_S_bc;
  Double I_ERI_Hy4z_Py_Px_S_bc = I_ERI_I2y4z_S_Px_S_bc+ABY*I_ERI_Hy4z_S_Px_S_bc;
  Double I_ERI_H4xz_Pz_Px_S_bc = I_ERI_I4x2z_S_Px_S_bc+ABZ*I_ERI_H4xz_S_Px_S_bc;
  Double I_ERI_H3xyz_Pz_Px_S_bc = I_ERI_I3xy2z_S_Px_S_bc+ABZ*I_ERI_H3xyz_S_Px_S_bc;
  Double I_ERI_H3x2z_Pz_Px_S_bc = I_ERI_I3x3z_S_Px_S_bc+ABZ*I_ERI_H3x2z_S_Px_S_bc;
  Double I_ERI_H2x2yz_Pz_Px_S_bc = I_ERI_I2x2y2z_S_Px_S_bc+ABZ*I_ERI_H2x2yz_S_Px_S_bc;
  Double I_ERI_H2xy2z_Pz_Px_S_bc = I_ERI_I2xy3z_S_Px_S_bc+ABZ*I_ERI_H2xy2z_S_Px_S_bc;
  Double I_ERI_H2x3z_Pz_Px_S_bc = I_ERI_I2x4z_S_Px_S_bc+ABZ*I_ERI_H2x3z_S_Px_S_bc;
  Double I_ERI_Hx3yz_Pz_Px_S_bc = I_ERI_Ix3y2z_S_Px_S_bc+ABZ*I_ERI_Hx3yz_S_Px_S_bc;
  Double I_ERI_Hx2y2z_Pz_Px_S_bc = I_ERI_Ix2y3z_S_Px_S_bc+ABZ*I_ERI_Hx2y2z_S_Px_S_bc;
  Double I_ERI_Hxy3z_Pz_Px_S_bc = I_ERI_Ixy4z_S_Px_S_bc+ABZ*I_ERI_Hxy3z_S_Px_S_bc;
  Double I_ERI_Hx4z_Pz_Px_S_bc = I_ERI_Ix5z_S_Px_S_bc+ABZ*I_ERI_Hx4z_S_Px_S_bc;
  Double I_ERI_H4yz_Pz_Px_S_bc = I_ERI_I4y2z_S_Px_S_bc+ABZ*I_ERI_H4yz_S_Px_S_bc;
  Double I_ERI_H3y2z_Pz_Px_S_bc = I_ERI_I3y3z_S_Px_S_bc+ABZ*I_ERI_H3y2z_S_Px_S_bc;
  Double I_ERI_H2y3z_Pz_Px_S_bc = I_ERI_I2y4z_S_Px_S_bc+ABZ*I_ERI_H2y3z_S_Px_S_bc;
  Double I_ERI_Hy4z_Pz_Px_S_bc = I_ERI_Iy5z_S_Px_S_bc+ABZ*I_ERI_Hy4z_S_Px_S_bc;
  Double I_ERI_H5z_Pz_Px_S_bc = I_ERI_I6z_S_Px_S_bc+ABZ*I_ERI_H5z_S_Px_S_bc;
  Double I_ERI_H5x_Px_Py_S_bc = I_ERI_I6x_S_Py_S_bc+ABX*I_ERI_H5x_S_Py_S_bc;
  Double I_ERI_H4xy_Px_Py_S_bc = I_ERI_I5xy_S_Py_S_bc+ABX*I_ERI_H4xy_S_Py_S_bc;
  Double I_ERI_H4xz_Px_Py_S_bc = I_ERI_I5xz_S_Py_S_bc+ABX*I_ERI_H4xz_S_Py_S_bc;
  Double I_ERI_H3x2y_Px_Py_S_bc = I_ERI_I4x2y_S_Py_S_bc+ABX*I_ERI_H3x2y_S_Py_S_bc;
  Double I_ERI_H3xyz_Px_Py_S_bc = I_ERI_I4xyz_S_Py_S_bc+ABX*I_ERI_H3xyz_S_Py_S_bc;
  Double I_ERI_H3x2z_Px_Py_S_bc = I_ERI_I4x2z_S_Py_S_bc+ABX*I_ERI_H3x2z_S_Py_S_bc;
  Double I_ERI_H2x3y_Px_Py_S_bc = I_ERI_I3x3y_S_Py_S_bc+ABX*I_ERI_H2x3y_S_Py_S_bc;
  Double I_ERI_H2x2yz_Px_Py_S_bc = I_ERI_I3x2yz_S_Py_S_bc+ABX*I_ERI_H2x2yz_S_Py_S_bc;
  Double I_ERI_H2xy2z_Px_Py_S_bc = I_ERI_I3xy2z_S_Py_S_bc+ABX*I_ERI_H2xy2z_S_Py_S_bc;
  Double I_ERI_H2x3z_Px_Py_S_bc = I_ERI_I3x3z_S_Py_S_bc+ABX*I_ERI_H2x3z_S_Py_S_bc;
  Double I_ERI_Hx4y_Px_Py_S_bc = I_ERI_I2x4y_S_Py_S_bc+ABX*I_ERI_Hx4y_S_Py_S_bc;
  Double I_ERI_Hx3yz_Px_Py_S_bc = I_ERI_I2x3yz_S_Py_S_bc+ABX*I_ERI_Hx3yz_S_Py_S_bc;
  Double I_ERI_Hx2y2z_Px_Py_S_bc = I_ERI_I2x2y2z_S_Py_S_bc+ABX*I_ERI_Hx2y2z_S_Py_S_bc;
  Double I_ERI_Hxy3z_Px_Py_S_bc = I_ERI_I2xy3z_S_Py_S_bc+ABX*I_ERI_Hxy3z_S_Py_S_bc;
  Double I_ERI_Hx4z_Px_Py_S_bc = I_ERI_I2x4z_S_Py_S_bc+ABX*I_ERI_Hx4z_S_Py_S_bc;
  Double I_ERI_H4yz_Px_Py_S_bc = I_ERI_Ix4yz_S_Py_S_bc+ABX*I_ERI_H4yz_S_Py_S_bc;
  Double I_ERI_H3y2z_Px_Py_S_bc = I_ERI_Ix3y2z_S_Py_S_bc+ABX*I_ERI_H3y2z_S_Py_S_bc;
  Double I_ERI_H2y3z_Px_Py_S_bc = I_ERI_Ix2y3z_S_Py_S_bc+ABX*I_ERI_H2y3z_S_Py_S_bc;
  Double I_ERI_Hy4z_Px_Py_S_bc = I_ERI_Ixy4z_S_Py_S_bc+ABX*I_ERI_Hy4z_S_Py_S_bc;
  Double I_ERI_H4xy_Py_Py_S_bc = I_ERI_I4x2y_S_Py_S_bc+ABY*I_ERI_H4xy_S_Py_S_bc;
  Double I_ERI_H3x2y_Py_Py_S_bc = I_ERI_I3x3y_S_Py_S_bc+ABY*I_ERI_H3x2y_S_Py_S_bc;
  Double I_ERI_H3xyz_Py_Py_S_bc = I_ERI_I3x2yz_S_Py_S_bc+ABY*I_ERI_H3xyz_S_Py_S_bc;
  Double I_ERI_H2x3y_Py_Py_S_bc = I_ERI_I2x4y_S_Py_S_bc+ABY*I_ERI_H2x3y_S_Py_S_bc;
  Double I_ERI_H2x2yz_Py_Py_S_bc = I_ERI_I2x3yz_S_Py_S_bc+ABY*I_ERI_H2x2yz_S_Py_S_bc;
  Double I_ERI_H2xy2z_Py_Py_S_bc = I_ERI_I2x2y2z_S_Py_S_bc+ABY*I_ERI_H2xy2z_S_Py_S_bc;
  Double I_ERI_Hx4y_Py_Py_S_bc = I_ERI_Ix5y_S_Py_S_bc+ABY*I_ERI_Hx4y_S_Py_S_bc;
  Double I_ERI_Hx3yz_Py_Py_S_bc = I_ERI_Ix4yz_S_Py_S_bc+ABY*I_ERI_Hx3yz_S_Py_S_bc;
  Double I_ERI_Hx2y2z_Py_Py_S_bc = I_ERI_Ix3y2z_S_Py_S_bc+ABY*I_ERI_Hx2y2z_S_Py_S_bc;
  Double I_ERI_Hxy3z_Py_Py_S_bc = I_ERI_Ix2y3z_S_Py_S_bc+ABY*I_ERI_Hxy3z_S_Py_S_bc;
  Double I_ERI_H5y_Py_Py_S_bc = I_ERI_I6y_S_Py_S_bc+ABY*I_ERI_H5y_S_Py_S_bc;
  Double I_ERI_H4yz_Py_Py_S_bc = I_ERI_I5yz_S_Py_S_bc+ABY*I_ERI_H4yz_S_Py_S_bc;
  Double I_ERI_H3y2z_Py_Py_S_bc = I_ERI_I4y2z_S_Py_S_bc+ABY*I_ERI_H3y2z_S_Py_S_bc;
  Double I_ERI_H2y3z_Py_Py_S_bc = I_ERI_I3y3z_S_Py_S_bc+ABY*I_ERI_H2y3z_S_Py_S_bc;
  Double I_ERI_Hy4z_Py_Py_S_bc = I_ERI_I2y4z_S_Py_S_bc+ABY*I_ERI_Hy4z_S_Py_S_bc;
  Double I_ERI_H4xz_Pz_Py_S_bc = I_ERI_I4x2z_S_Py_S_bc+ABZ*I_ERI_H4xz_S_Py_S_bc;
  Double I_ERI_H3xyz_Pz_Py_S_bc = I_ERI_I3xy2z_S_Py_S_bc+ABZ*I_ERI_H3xyz_S_Py_S_bc;
  Double I_ERI_H3x2z_Pz_Py_S_bc = I_ERI_I3x3z_S_Py_S_bc+ABZ*I_ERI_H3x2z_S_Py_S_bc;
  Double I_ERI_H2x2yz_Pz_Py_S_bc = I_ERI_I2x2y2z_S_Py_S_bc+ABZ*I_ERI_H2x2yz_S_Py_S_bc;
  Double I_ERI_H2xy2z_Pz_Py_S_bc = I_ERI_I2xy3z_S_Py_S_bc+ABZ*I_ERI_H2xy2z_S_Py_S_bc;
  Double I_ERI_H2x3z_Pz_Py_S_bc = I_ERI_I2x4z_S_Py_S_bc+ABZ*I_ERI_H2x3z_S_Py_S_bc;
  Double I_ERI_Hx3yz_Pz_Py_S_bc = I_ERI_Ix3y2z_S_Py_S_bc+ABZ*I_ERI_Hx3yz_S_Py_S_bc;
  Double I_ERI_Hx2y2z_Pz_Py_S_bc = I_ERI_Ix2y3z_S_Py_S_bc+ABZ*I_ERI_Hx2y2z_S_Py_S_bc;
  Double I_ERI_Hxy3z_Pz_Py_S_bc = I_ERI_Ixy4z_S_Py_S_bc+ABZ*I_ERI_Hxy3z_S_Py_S_bc;
  Double I_ERI_Hx4z_Pz_Py_S_bc = I_ERI_Ix5z_S_Py_S_bc+ABZ*I_ERI_Hx4z_S_Py_S_bc;
  Double I_ERI_H4yz_Pz_Py_S_bc = I_ERI_I4y2z_S_Py_S_bc+ABZ*I_ERI_H4yz_S_Py_S_bc;
  Double I_ERI_H3y2z_Pz_Py_S_bc = I_ERI_I3y3z_S_Py_S_bc+ABZ*I_ERI_H3y2z_S_Py_S_bc;
  Double I_ERI_H2y3z_Pz_Py_S_bc = I_ERI_I2y4z_S_Py_S_bc+ABZ*I_ERI_H2y3z_S_Py_S_bc;
  Double I_ERI_Hy4z_Pz_Py_S_bc = I_ERI_Iy5z_S_Py_S_bc+ABZ*I_ERI_Hy4z_S_Py_S_bc;
  Double I_ERI_H5z_Pz_Py_S_bc = I_ERI_I6z_S_Py_S_bc+ABZ*I_ERI_H5z_S_Py_S_bc;
  Double I_ERI_H5x_Px_Pz_S_bc = I_ERI_I6x_S_Pz_S_bc+ABX*I_ERI_H5x_S_Pz_S_bc;
  Double I_ERI_H4xy_Px_Pz_S_bc = I_ERI_I5xy_S_Pz_S_bc+ABX*I_ERI_H4xy_S_Pz_S_bc;
  Double I_ERI_H4xz_Px_Pz_S_bc = I_ERI_I5xz_S_Pz_S_bc+ABX*I_ERI_H4xz_S_Pz_S_bc;
  Double I_ERI_H3x2y_Px_Pz_S_bc = I_ERI_I4x2y_S_Pz_S_bc+ABX*I_ERI_H3x2y_S_Pz_S_bc;
  Double I_ERI_H3xyz_Px_Pz_S_bc = I_ERI_I4xyz_S_Pz_S_bc+ABX*I_ERI_H3xyz_S_Pz_S_bc;
  Double I_ERI_H3x2z_Px_Pz_S_bc = I_ERI_I4x2z_S_Pz_S_bc+ABX*I_ERI_H3x2z_S_Pz_S_bc;
  Double I_ERI_H2x3y_Px_Pz_S_bc = I_ERI_I3x3y_S_Pz_S_bc+ABX*I_ERI_H2x3y_S_Pz_S_bc;
  Double I_ERI_H2x2yz_Px_Pz_S_bc = I_ERI_I3x2yz_S_Pz_S_bc+ABX*I_ERI_H2x2yz_S_Pz_S_bc;
  Double I_ERI_H2xy2z_Px_Pz_S_bc = I_ERI_I3xy2z_S_Pz_S_bc+ABX*I_ERI_H2xy2z_S_Pz_S_bc;
  Double I_ERI_H2x3z_Px_Pz_S_bc = I_ERI_I3x3z_S_Pz_S_bc+ABX*I_ERI_H2x3z_S_Pz_S_bc;
  Double I_ERI_Hx4y_Px_Pz_S_bc = I_ERI_I2x4y_S_Pz_S_bc+ABX*I_ERI_Hx4y_S_Pz_S_bc;
  Double I_ERI_Hx3yz_Px_Pz_S_bc = I_ERI_I2x3yz_S_Pz_S_bc+ABX*I_ERI_Hx3yz_S_Pz_S_bc;
  Double I_ERI_Hx2y2z_Px_Pz_S_bc = I_ERI_I2x2y2z_S_Pz_S_bc+ABX*I_ERI_Hx2y2z_S_Pz_S_bc;
  Double I_ERI_Hxy3z_Px_Pz_S_bc = I_ERI_I2xy3z_S_Pz_S_bc+ABX*I_ERI_Hxy3z_S_Pz_S_bc;
  Double I_ERI_Hx4z_Px_Pz_S_bc = I_ERI_I2x4z_S_Pz_S_bc+ABX*I_ERI_Hx4z_S_Pz_S_bc;
  Double I_ERI_H4yz_Px_Pz_S_bc = I_ERI_Ix4yz_S_Pz_S_bc+ABX*I_ERI_H4yz_S_Pz_S_bc;
  Double I_ERI_H3y2z_Px_Pz_S_bc = I_ERI_Ix3y2z_S_Pz_S_bc+ABX*I_ERI_H3y2z_S_Pz_S_bc;
  Double I_ERI_H2y3z_Px_Pz_S_bc = I_ERI_Ix2y3z_S_Pz_S_bc+ABX*I_ERI_H2y3z_S_Pz_S_bc;
  Double I_ERI_Hy4z_Px_Pz_S_bc = I_ERI_Ixy4z_S_Pz_S_bc+ABX*I_ERI_Hy4z_S_Pz_S_bc;
  Double I_ERI_H4xy_Py_Pz_S_bc = I_ERI_I4x2y_S_Pz_S_bc+ABY*I_ERI_H4xy_S_Pz_S_bc;
  Double I_ERI_H3x2y_Py_Pz_S_bc = I_ERI_I3x3y_S_Pz_S_bc+ABY*I_ERI_H3x2y_S_Pz_S_bc;
  Double I_ERI_H3xyz_Py_Pz_S_bc = I_ERI_I3x2yz_S_Pz_S_bc+ABY*I_ERI_H3xyz_S_Pz_S_bc;
  Double I_ERI_H2x3y_Py_Pz_S_bc = I_ERI_I2x4y_S_Pz_S_bc+ABY*I_ERI_H2x3y_S_Pz_S_bc;
  Double I_ERI_H2x2yz_Py_Pz_S_bc = I_ERI_I2x3yz_S_Pz_S_bc+ABY*I_ERI_H2x2yz_S_Pz_S_bc;
  Double I_ERI_H2xy2z_Py_Pz_S_bc = I_ERI_I2x2y2z_S_Pz_S_bc+ABY*I_ERI_H2xy2z_S_Pz_S_bc;
  Double I_ERI_Hx4y_Py_Pz_S_bc = I_ERI_Ix5y_S_Pz_S_bc+ABY*I_ERI_Hx4y_S_Pz_S_bc;
  Double I_ERI_Hx3yz_Py_Pz_S_bc = I_ERI_Ix4yz_S_Pz_S_bc+ABY*I_ERI_Hx3yz_S_Pz_S_bc;
  Double I_ERI_Hx2y2z_Py_Pz_S_bc = I_ERI_Ix3y2z_S_Pz_S_bc+ABY*I_ERI_Hx2y2z_S_Pz_S_bc;
  Double I_ERI_Hxy3z_Py_Pz_S_bc = I_ERI_Ix2y3z_S_Pz_S_bc+ABY*I_ERI_Hxy3z_S_Pz_S_bc;
  Double I_ERI_H5y_Py_Pz_S_bc = I_ERI_I6y_S_Pz_S_bc+ABY*I_ERI_H5y_S_Pz_S_bc;
  Double I_ERI_H4yz_Py_Pz_S_bc = I_ERI_I5yz_S_Pz_S_bc+ABY*I_ERI_H4yz_S_Pz_S_bc;
  Double I_ERI_H3y2z_Py_Pz_S_bc = I_ERI_I4y2z_S_Pz_S_bc+ABY*I_ERI_H3y2z_S_Pz_S_bc;
  Double I_ERI_H2y3z_Py_Pz_S_bc = I_ERI_I3y3z_S_Pz_S_bc+ABY*I_ERI_H2y3z_S_Pz_S_bc;
  Double I_ERI_Hy4z_Py_Pz_S_bc = I_ERI_I2y4z_S_Pz_S_bc+ABY*I_ERI_Hy4z_S_Pz_S_bc;
  Double I_ERI_H4xz_Pz_Pz_S_bc = I_ERI_I4x2z_S_Pz_S_bc+ABZ*I_ERI_H4xz_S_Pz_S_bc;
  Double I_ERI_H3xyz_Pz_Pz_S_bc = I_ERI_I3xy2z_S_Pz_S_bc+ABZ*I_ERI_H3xyz_S_Pz_S_bc;
  Double I_ERI_H3x2z_Pz_Pz_S_bc = I_ERI_I3x3z_S_Pz_S_bc+ABZ*I_ERI_H3x2z_S_Pz_S_bc;
  Double I_ERI_H2x2yz_Pz_Pz_S_bc = I_ERI_I2x2y2z_S_Pz_S_bc+ABZ*I_ERI_H2x2yz_S_Pz_S_bc;
  Double I_ERI_H2xy2z_Pz_Pz_S_bc = I_ERI_I2xy3z_S_Pz_S_bc+ABZ*I_ERI_H2xy2z_S_Pz_S_bc;
  Double I_ERI_H2x3z_Pz_Pz_S_bc = I_ERI_I2x4z_S_Pz_S_bc+ABZ*I_ERI_H2x3z_S_Pz_S_bc;
  Double I_ERI_Hx3yz_Pz_Pz_S_bc = I_ERI_Ix3y2z_S_Pz_S_bc+ABZ*I_ERI_Hx3yz_S_Pz_S_bc;
  Double I_ERI_Hx2y2z_Pz_Pz_S_bc = I_ERI_Ix2y3z_S_Pz_S_bc+ABZ*I_ERI_Hx2y2z_S_Pz_S_bc;
  Double I_ERI_Hxy3z_Pz_Pz_S_bc = I_ERI_Ixy4z_S_Pz_S_bc+ABZ*I_ERI_Hxy3z_S_Pz_S_bc;
  Double I_ERI_Hx4z_Pz_Pz_S_bc = I_ERI_Ix5z_S_Pz_S_bc+ABZ*I_ERI_Hx4z_S_Pz_S_bc;
  Double I_ERI_H4yz_Pz_Pz_S_bc = I_ERI_I4y2z_S_Pz_S_bc+ABZ*I_ERI_H4yz_S_Pz_S_bc;
  Double I_ERI_H3y2z_Pz_Pz_S_bc = I_ERI_I3y3z_S_Pz_S_bc+ABZ*I_ERI_H3y2z_S_Pz_S_bc;
  Double I_ERI_H2y3z_Pz_Pz_S_bc = I_ERI_I2y4z_S_Pz_S_bc+ABZ*I_ERI_H2y3z_S_Pz_S_bc;
  Double I_ERI_Hy4z_Pz_Pz_S_bc = I_ERI_Iy5z_S_Pz_S_bc+ABZ*I_ERI_Hy4z_S_Pz_S_bc;
  Double I_ERI_H5z_Pz_Pz_S_bc = I_ERI_I6z_S_Pz_S_bc+ABZ*I_ERI_H5z_S_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_G_D_P_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 105 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_P_P_S_bc
   * RHS shell quartet name: SQ_ERI_G_P_P_S_bc
   ************************************************************/
  Double I_ERI_G4x_D2x_Px_S_bc = I_ERI_H5x_Px_Px_S_bc+ABX*I_ERI_G4x_Px_Px_S_bc;
  Double I_ERI_G3xy_D2x_Px_S_bc = I_ERI_H4xy_Px_Px_S_bc+ABX*I_ERI_G3xy_Px_Px_S_bc;
  Double I_ERI_G3xz_D2x_Px_S_bc = I_ERI_H4xz_Px_Px_S_bc+ABX*I_ERI_G3xz_Px_Px_S_bc;
  Double I_ERI_G2x2y_D2x_Px_S_bc = I_ERI_H3x2y_Px_Px_S_bc+ABX*I_ERI_G2x2y_Px_Px_S_bc;
  Double I_ERI_G2xyz_D2x_Px_S_bc = I_ERI_H3xyz_Px_Px_S_bc+ABX*I_ERI_G2xyz_Px_Px_S_bc;
  Double I_ERI_G2x2z_D2x_Px_S_bc = I_ERI_H3x2z_Px_Px_S_bc+ABX*I_ERI_G2x2z_Px_Px_S_bc;
  Double I_ERI_Gx3y_D2x_Px_S_bc = I_ERI_H2x3y_Px_Px_S_bc+ABX*I_ERI_Gx3y_Px_Px_S_bc;
  Double I_ERI_Gx2yz_D2x_Px_S_bc = I_ERI_H2x2yz_Px_Px_S_bc+ABX*I_ERI_Gx2yz_Px_Px_S_bc;
  Double I_ERI_Gxy2z_D2x_Px_S_bc = I_ERI_H2xy2z_Px_Px_S_bc+ABX*I_ERI_Gxy2z_Px_Px_S_bc;
  Double I_ERI_Gx3z_D2x_Px_S_bc = I_ERI_H2x3z_Px_Px_S_bc+ABX*I_ERI_Gx3z_Px_Px_S_bc;
  Double I_ERI_G4y_D2x_Px_S_bc = I_ERI_Hx4y_Px_Px_S_bc+ABX*I_ERI_G4y_Px_Px_S_bc;
  Double I_ERI_G3yz_D2x_Px_S_bc = I_ERI_Hx3yz_Px_Px_S_bc+ABX*I_ERI_G3yz_Px_Px_S_bc;
  Double I_ERI_G2y2z_D2x_Px_S_bc = I_ERI_Hx2y2z_Px_Px_S_bc+ABX*I_ERI_G2y2z_Px_Px_S_bc;
  Double I_ERI_Gy3z_D2x_Px_S_bc = I_ERI_Hxy3z_Px_Px_S_bc+ABX*I_ERI_Gy3z_Px_Px_S_bc;
  Double I_ERI_G4z_D2x_Px_S_bc = I_ERI_Hx4z_Px_Px_S_bc+ABX*I_ERI_G4z_Px_Px_S_bc;
  Double I_ERI_G3xz_Dxy_Px_S_bc = I_ERI_H3xyz_Px_Px_S_bc+ABY*I_ERI_G3xz_Px_Px_S_bc;
  Double I_ERI_G2xyz_Dxy_Px_S_bc = I_ERI_H2x2yz_Px_Px_S_bc+ABY*I_ERI_G2xyz_Px_Px_S_bc;
  Double I_ERI_G2x2z_Dxy_Px_S_bc = I_ERI_H2xy2z_Px_Px_S_bc+ABY*I_ERI_G2x2z_Px_Px_S_bc;
  Double I_ERI_Gx2yz_Dxy_Px_S_bc = I_ERI_Hx3yz_Px_Px_S_bc+ABY*I_ERI_Gx2yz_Px_Px_S_bc;
  Double I_ERI_Gxy2z_Dxy_Px_S_bc = I_ERI_Hx2y2z_Px_Px_S_bc+ABY*I_ERI_Gxy2z_Px_Px_S_bc;
  Double I_ERI_Gx3z_Dxy_Px_S_bc = I_ERI_Hxy3z_Px_Px_S_bc+ABY*I_ERI_Gx3z_Px_Px_S_bc;
  Double I_ERI_G3yz_Dxy_Px_S_bc = I_ERI_H4yz_Px_Px_S_bc+ABY*I_ERI_G3yz_Px_Px_S_bc;
  Double I_ERI_G2y2z_Dxy_Px_S_bc = I_ERI_H3y2z_Px_Px_S_bc+ABY*I_ERI_G2y2z_Px_Px_S_bc;
  Double I_ERI_Gy3z_Dxy_Px_S_bc = I_ERI_H2y3z_Px_Px_S_bc+ABY*I_ERI_Gy3z_Px_Px_S_bc;
  Double I_ERI_G4z_Dxy_Px_S_bc = I_ERI_Hy4z_Px_Px_S_bc+ABY*I_ERI_G4z_Px_Px_S_bc;
  Double I_ERI_G4x_D2y_Px_S_bc = I_ERI_H4xy_Py_Px_S_bc+ABY*I_ERI_G4x_Py_Px_S_bc;
  Double I_ERI_G3xy_D2y_Px_S_bc = I_ERI_H3x2y_Py_Px_S_bc+ABY*I_ERI_G3xy_Py_Px_S_bc;
  Double I_ERI_G3xz_D2y_Px_S_bc = I_ERI_H3xyz_Py_Px_S_bc+ABY*I_ERI_G3xz_Py_Px_S_bc;
  Double I_ERI_G2x2y_D2y_Px_S_bc = I_ERI_H2x3y_Py_Px_S_bc+ABY*I_ERI_G2x2y_Py_Px_S_bc;
  Double I_ERI_G2xyz_D2y_Px_S_bc = I_ERI_H2x2yz_Py_Px_S_bc+ABY*I_ERI_G2xyz_Py_Px_S_bc;
  Double I_ERI_G2x2z_D2y_Px_S_bc = I_ERI_H2xy2z_Py_Px_S_bc+ABY*I_ERI_G2x2z_Py_Px_S_bc;
  Double I_ERI_Gx3y_D2y_Px_S_bc = I_ERI_Hx4y_Py_Px_S_bc+ABY*I_ERI_Gx3y_Py_Px_S_bc;
  Double I_ERI_Gx2yz_D2y_Px_S_bc = I_ERI_Hx3yz_Py_Px_S_bc+ABY*I_ERI_Gx2yz_Py_Px_S_bc;
  Double I_ERI_Gxy2z_D2y_Px_S_bc = I_ERI_Hx2y2z_Py_Px_S_bc+ABY*I_ERI_Gxy2z_Py_Px_S_bc;
  Double I_ERI_Gx3z_D2y_Px_S_bc = I_ERI_Hxy3z_Py_Px_S_bc+ABY*I_ERI_Gx3z_Py_Px_S_bc;
  Double I_ERI_G4y_D2y_Px_S_bc = I_ERI_H5y_Py_Px_S_bc+ABY*I_ERI_G4y_Py_Px_S_bc;
  Double I_ERI_G3yz_D2y_Px_S_bc = I_ERI_H4yz_Py_Px_S_bc+ABY*I_ERI_G3yz_Py_Px_S_bc;
  Double I_ERI_G2y2z_D2y_Px_S_bc = I_ERI_H3y2z_Py_Px_S_bc+ABY*I_ERI_G2y2z_Py_Px_S_bc;
  Double I_ERI_Gy3z_D2y_Px_S_bc = I_ERI_H2y3z_Py_Px_S_bc+ABY*I_ERI_Gy3z_Py_Px_S_bc;
  Double I_ERI_G4z_D2y_Px_S_bc = I_ERI_Hy4z_Py_Px_S_bc+ABY*I_ERI_G4z_Py_Px_S_bc;
  Double I_ERI_G4x_D2z_Px_S_bc = I_ERI_H4xz_Pz_Px_S_bc+ABZ*I_ERI_G4x_Pz_Px_S_bc;
  Double I_ERI_G3xy_D2z_Px_S_bc = I_ERI_H3xyz_Pz_Px_S_bc+ABZ*I_ERI_G3xy_Pz_Px_S_bc;
  Double I_ERI_G3xz_D2z_Px_S_bc = I_ERI_H3x2z_Pz_Px_S_bc+ABZ*I_ERI_G3xz_Pz_Px_S_bc;
  Double I_ERI_G2x2y_D2z_Px_S_bc = I_ERI_H2x2yz_Pz_Px_S_bc+ABZ*I_ERI_G2x2y_Pz_Px_S_bc;
  Double I_ERI_G2xyz_D2z_Px_S_bc = I_ERI_H2xy2z_Pz_Px_S_bc+ABZ*I_ERI_G2xyz_Pz_Px_S_bc;
  Double I_ERI_G2x2z_D2z_Px_S_bc = I_ERI_H2x3z_Pz_Px_S_bc+ABZ*I_ERI_G2x2z_Pz_Px_S_bc;
  Double I_ERI_Gx3y_D2z_Px_S_bc = I_ERI_Hx3yz_Pz_Px_S_bc+ABZ*I_ERI_Gx3y_Pz_Px_S_bc;
  Double I_ERI_Gx2yz_D2z_Px_S_bc = I_ERI_Hx2y2z_Pz_Px_S_bc+ABZ*I_ERI_Gx2yz_Pz_Px_S_bc;
  Double I_ERI_Gxy2z_D2z_Px_S_bc = I_ERI_Hxy3z_Pz_Px_S_bc+ABZ*I_ERI_Gxy2z_Pz_Px_S_bc;
  Double I_ERI_Gx3z_D2z_Px_S_bc = I_ERI_Hx4z_Pz_Px_S_bc+ABZ*I_ERI_Gx3z_Pz_Px_S_bc;
  Double I_ERI_G4y_D2z_Px_S_bc = I_ERI_H4yz_Pz_Px_S_bc+ABZ*I_ERI_G4y_Pz_Px_S_bc;
  Double I_ERI_G3yz_D2z_Px_S_bc = I_ERI_H3y2z_Pz_Px_S_bc+ABZ*I_ERI_G3yz_Pz_Px_S_bc;
  Double I_ERI_G2y2z_D2z_Px_S_bc = I_ERI_H2y3z_Pz_Px_S_bc+ABZ*I_ERI_G2y2z_Pz_Px_S_bc;
  Double I_ERI_Gy3z_D2z_Px_S_bc = I_ERI_Hy4z_Pz_Px_S_bc+ABZ*I_ERI_Gy3z_Pz_Px_S_bc;
  Double I_ERI_G4z_D2z_Px_S_bc = I_ERI_H5z_Pz_Px_S_bc+ABZ*I_ERI_G4z_Pz_Px_S_bc;
  Double I_ERI_G4x_D2x_Py_S_bc = I_ERI_H5x_Px_Py_S_bc+ABX*I_ERI_G4x_Px_Py_S_bc;
  Double I_ERI_G3xy_D2x_Py_S_bc = I_ERI_H4xy_Px_Py_S_bc+ABX*I_ERI_G3xy_Px_Py_S_bc;
  Double I_ERI_G3xz_D2x_Py_S_bc = I_ERI_H4xz_Px_Py_S_bc+ABX*I_ERI_G3xz_Px_Py_S_bc;
  Double I_ERI_G2x2y_D2x_Py_S_bc = I_ERI_H3x2y_Px_Py_S_bc+ABX*I_ERI_G2x2y_Px_Py_S_bc;
  Double I_ERI_G2xyz_D2x_Py_S_bc = I_ERI_H3xyz_Px_Py_S_bc+ABX*I_ERI_G2xyz_Px_Py_S_bc;
  Double I_ERI_G2x2z_D2x_Py_S_bc = I_ERI_H3x2z_Px_Py_S_bc+ABX*I_ERI_G2x2z_Px_Py_S_bc;
  Double I_ERI_Gx3y_D2x_Py_S_bc = I_ERI_H2x3y_Px_Py_S_bc+ABX*I_ERI_Gx3y_Px_Py_S_bc;
  Double I_ERI_Gx2yz_D2x_Py_S_bc = I_ERI_H2x2yz_Px_Py_S_bc+ABX*I_ERI_Gx2yz_Px_Py_S_bc;
  Double I_ERI_Gxy2z_D2x_Py_S_bc = I_ERI_H2xy2z_Px_Py_S_bc+ABX*I_ERI_Gxy2z_Px_Py_S_bc;
  Double I_ERI_Gx3z_D2x_Py_S_bc = I_ERI_H2x3z_Px_Py_S_bc+ABX*I_ERI_Gx3z_Px_Py_S_bc;
  Double I_ERI_G4y_D2x_Py_S_bc = I_ERI_Hx4y_Px_Py_S_bc+ABX*I_ERI_G4y_Px_Py_S_bc;
  Double I_ERI_G3yz_D2x_Py_S_bc = I_ERI_Hx3yz_Px_Py_S_bc+ABX*I_ERI_G3yz_Px_Py_S_bc;
  Double I_ERI_G2y2z_D2x_Py_S_bc = I_ERI_Hx2y2z_Px_Py_S_bc+ABX*I_ERI_G2y2z_Px_Py_S_bc;
  Double I_ERI_Gy3z_D2x_Py_S_bc = I_ERI_Hxy3z_Px_Py_S_bc+ABX*I_ERI_Gy3z_Px_Py_S_bc;
  Double I_ERI_G4z_D2x_Py_S_bc = I_ERI_Hx4z_Px_Py_S_bc+ABX*I_ERI_G4z_Px_Py_S_bc;
  Double I_ERI_G3xz_Dxy_Py_S_bc = I_ERI_H3xyz_Px_Py_S_bc+ABY*I_ERI_G3xz_Px_Py_S_bc;
  Double I_ERI_G2xyz_Dxy_Py_S_bc = I_ERI_H2x2yz_Px_Py_S_bc+ABY*I_ERI_G2xyz_Px_Py_S_bc;
  Double I_ERI_G2x2z_Dxy_Py_S_bc = I_ERI_H2xy2z_Px_Py_S_bc+ABY*I_ERI_G2x2z_Px_Py_S_bc;
  Double I_ERI_Gx2yz_Dxy_Py_S_bc = I_ERI_Hx3yz_Px_Py_S_bc+ABY*I_ERI_Gx2yz_Px_Py_S_bc;
  Double I_ERI_Gxy2z_Dxy_Py_S_bc = I_ERI_Hx2y2z_Px_Py_S_bc+ABY*I_ERI_Gxy2z_Px_Py_S_bc;
  Double I_ERI_Gx3z_Dxy_Py_S_bc = I_ERI_Hxy3z_Px_Py_S_bc+ABY*I_ERI_Gx3z_Px_Py_S_bc;
  Double I_ERI_G3yz_Dxy_Py_S_bc = I_ERI_H4yz_Px_Py_S_bc+ABY*I_ERI_G3yz_Px_Py_S_bc;
  Double I_ERI_G2y2z_Dxy_Py_S_bc = I_ERI_H3y2z_Px_Py_S_bc+ABY*I_ERI_G2y2z_Px_Py_S_bc;
  Double I_ERI_Gy3z_Dxy_Py_S_bc = I_ERI_H2y3z_Px_Py_S_bc+ABY*I_ERI_Gy3z_Px_Py_S_bc;
  Double I_ERI_G4z_Dxy_Py_S_bc = I_ERI_Hy4z_Px_Py_S_bc+ABY*I_ERI_G4z_Px_Py_S_bc;
  Double I_ERI_G4x_D2y_Py_S_bc = I_ERI_H4xy_Py_Py_S_bc+ABY*I_ERI_G4x_Py_Py_S_bc;
  Double I_ERI_G3xy_D2y_Py_S_bc = I_ERI_H3x2y_Py_Py_S_bc+ABY*I_ERI_G3xy_Py_Py_S_bc;
  Double I_ERI_G3xz_D2y_Py_S_bc = I_ERI_H3xyz_Py_Py_S_bc+ABY*I_ERI_G3xz_Py_Py_S_bc;
  Double I_ERI_G2x2y_D2y_Py_S_bc = I_ERI_H2x3y_Py_Py_S_bc+ABY*I_ERI_G2x2y_Py_Py_S_bc;
  Double I_ERI_G2xyz_D2y_Py_S_bc = I_ERI_H2x2yz_Py_Py_S_bc+ABY*I_ERI_G2xyz_Py_Py_S_bc;
  Double I_ERI_G2x2z_D2y_Py_S_bc = I_ERI_H2xy2z_Py_Py_S_bc+ABY*I_ERI_G2x2z_Py_Py_S_bc;
  Double I_ERI_Gx3y_D2y_Py_S_bc = I_ERI_Hx4y_Py_Py_S_bc+ABY*I_ERI_Gx3y_Py_Py_S_bc;
  Double I_ERI_Gx2yz_D2y_Py_S_bc = I_ERI_Hx3yz_Py_Py_S_bc+ABY*I_ERI_Gx2yz_Py_Py_S_bc;
  Double I_ERI_Gxy2z_D2y_Py_S_bc = I_ERI_Hx2y2z_Py_Py_S_bc+ABY*I_ERI_Gxy2z_Py_Py_S_bc;
  Double I_ERI_Gx3z_D2y_Py_S_bc = I_ERI_Hxy3z_Py_Py_S_bc+ABY*I_ERI_Gx3z_Py_Py_S_bc;
  Double I_ERI_G4y_D2y_Py_S_bc = I_ERI_H5y_Py_Py_S_bc+ABY*I_ERI_G4y_Py_Py_S_bc;
  Double I_ERI_G3yz_D2y_Py_S_bc = I_ERI_H4yz_Py_Py_S_bc+ABY*I_ERI_G3yz_Py_Py_S_bc;
  Double I_ERI_G2y2z_D2y_Py_S_bc = I_ERI_H3y2z_Py_Py_S_bc+ABY*I_ERI_G2y2z_Py_Py_S_bc;
  Double I_ERI_Gy3z_D2y_Py_S_bc = I_ERI_H2y3z_Py_Py_S_bc+ABY*I_ERI_Gy3z_Py_Py_S_bc;
  Double I_ERI_G4z_D2y_Py_S_bc = I_ERI_Hy4z_Py_Py_S_bc+ABY*I_ERI_G4z_Py_Py_S_bc;
  Double I_ERI_G4x_D2z_Py_S_bc = I_ERI_H4xz_Pz_Py_S_bc+ABZ*I_ERI_G4x_Pz_Py_S_bc;
  Double I_ERI_G3xy_D2z_Py_S_bc = I_ERI_H3xyz_Pz_Py_S_bc+ABZ*I_ERI_G3xy_Pz_Py_S_bc;
  Double I_ERI_G3xz_D2z_Py_S_bc = I_ERI_H3x2z_Pz_Py_S_bc+ABZ*I_ERI_G3xz_Pz_Py_S_bc;
  Double I_ERI_G2x2y_D2z_Py_S_bc = I_ERI_H2x2yz_Pz_Py_S_bc+ABZ*I_ERI_G2x2y_Pz_Py_S_bc;
  Double I_ERI_G2xyz_D2z_Py_S_bc = I_ERI_H2xy2z_Pz_Py_S_bc+ABZ*I_ERI_G2xyz_Pz_Py_S_bc;
  Double I_ERI_G2x2z_D2z_Py_S_bc = I_ERI_H2x3z_Pz_Py_S_bc+ABZ*I_ERI_G2x2z_Pz_Py_S_bc;
  Double I_ERI_Gx3y_D2z_Py_S_bc = I_ERI_Hx3yz_Pz_Py_S_bc+ABZ*I_ERI_Gx3y_Pz_Py_S_bc;
  Double I_ERI_Gx2yz_D2z_Py_S_bc = I_ERI_Hx2y2z_Pz_Py_S_bc+ABZ*I_ERI_Gx2yz_Pz_Py_S_bc;
  Double I_ERI_Gxy2z_D2z_Py_S_bc = I_ERI_Hxy3z_Pz_Py_S_bc+ABZ*I_ERI_Gxy2z_Pz_Py_S_bc;
  Double I_ERI_Gx3z_D2z_Py_S_bc = I_ERI_Hx4z_Pz_Py_S_bc+ABZ*I_ERI_Gx3z_Pz_Py_S_bc;
  Double I_ERI_G4y_D2z_Py_S_bc = I_ERI_H4yz_Pz_Py_S_bc+ABZ*I_ERI_G4y_Pz_Py_S_bc;
  Double I_ERI_G3yz_D2z_Py_S_bc = I_ERI_H3y2z_Pz_Py_S_bc+ABZ*I_ERI_G3yz_Pz_Py_S_bc;
  Double I_ERI_G2y2z_D2z_Py_S_bc = I_ERI_H2y3z_Pz_Py_S_bc+ABZ*I_ERI_G2y2z_Pz_Py_S_bc;
  Double I_ERI_Gy3z_D2z_Py_S_bc = I_ERI_Hy4z_Pz_Py_S_bc+ABZ*I_ERI_Gy3z_Pz_Py_S_bc;
  Double I_ERI_G4z_D2z_Py_S_bc = I_ERI_H5z_Pz_Py_S_bc+ABZ*I_ERI_G4z_Pz_Py_S_bc;
  Double I_ERI_G4x_D2x_Pz_S_bc = I_ERI_H5x_Px_Pz_S_bc+ABX*I_ERI_G4x_Px_Pz_S_bc;
  Double I_ERI_G3xy_D2x_Pz_S_bc = I_ERI_H4xy_Px_Pz_S_bc+ABX*I_ERI_G3xy_Px_Pz_S_bc;
  Double I_ERI_G3xz_D2x_Pz_S_bc = I_ERI_H4xz_Px_Pz_S_bc+ABX*I_ERI_G3xz_Px_Pz_S_bc;
  Double I_ERI_G2x2y_D2x_Pz_S_bc = I_ERI_H3x2y_Px_Pz_S_bc+ABX*I_ERI_G2x2y_Px_Pz_S_bc;
  Double I_ERI_G2xyz_D2x_Pz_S_bc = I_ERI_H3xyz_Px_Pz_S_bc+ABX*I_ERI_G2xyz_Px_Pz_S_bc;
  Double I_ERI_G2x2z_D2x_Pz_S_bc = I_ERI_H3x2z_Px_Pz_S_bc+ABX*I_ERI_G2x2z_Px_Pz_S_bc;
  Double I_ERI_Gx3y_D2x_Pz_S_bc = I_ERI_H2x3y_Px_Pz_S_bc+ABX*I_ERI_Gx3y_Px_Pz_S_bc;
  Double I_ERI_Gx2yz_D2x_Pz_S_bc = I_ERI_H2x2yz_Px_Pz_S_bc+ABX*I_ERI_Gx2yz_Px_Pz_S_bc;
  Double I_ERI_Gxy2z_D2x_Pz_S_bc = I_ERI_H2xy2z_Px_Pz_S_bc+ABX*I_ERI_Gxy2z_Px_Pz_S_bc;
  Double I_ERI_Gx3z_D2x_Pz_S_bc = I_ERI_H2x3z_Px_Pz_S_bc+ABX*I_ERI_Gx3z_Px_Pz_S_bc;
  Double I_ERI_G4y_D2x_Pz_S_bc = I_ERI_Hx4y_Px_Pz_S_bc+ABX*I_ERI_G4y_Px_Pz_S_bc;
  Double I_ERI_G3yz_D2x_Pz_S_bc = I_ERI_Hx3yz_Px_Pz_S_bc+ABX*I_ERI_G3yz_Px_Pz_S_bc;
  Double I_ERI_G2y2z_D2x_Pz_S_bc = I_ERI_Hx2y2z_Px_Pz_S_bc+ABX*I_ERI_G2y2z_Px_Pz_S_bc;
  Double I_ERI_Gy3z_D2x_Pz_S_bc = I_ERI_Hxy3z_Px_Pz_S_bc+ABX*I_ERI_Gy3z_Px_Pz_S_bc;
  Double I_ERI_G4z_D2x_Pz_S_bc = I_ERI_Hx4z_Px_Pz_S_bc+ABX*I_ERI_G4z_Px_Pz_S_bc;
  Double I_ERI_G3xz_Dxy_Pz_S_bc = I_ERI_H3xyz_Px_Pz_S_bc+ABY*I_ERI_G3xz_Px_Pz_S_bc;
  Double I_ERI_G2xyz_Dxy_Pz_S_bc = I_ERI_H2x2yz_Px_Pz_S_bc+ABY*I_ERI_G2xyz_Px_Pz_S_bc;
  Double I_ERI_G2x2z_Dxy_Pz_S_bc = I_ERI_H2xy2z_Px_Pz_S_bc+ABY*I_ERI_G2x2z_Px_Pz_S_bc;
  Double I_ERI_Gx2yz_Dxy_Pz_S_bc = I_ERI_Hx3yz_Px_Pz_S_bc+ABY*I_ERI_Gx2yz_Px_Pz_S_bc;
  Double I_ERI_Gxy2z_Dxy_Pz_S_bc = I_ERI_Hx2y2z_Px_Pz_S_bc+ABY*I_ERI_Gxy2z_Px_Pz_S_bc;
  Double I_ERI_Gx3z_Dxy_Pz_S_bc = I_ERI_Hxy3z_Px_Pz_S_bc+ABY*I_ERI_Gx3z_Px_Pz_S_bc;
  Double I_ERI_G3yz_Dxy_Pz_S_bc = I_ERI_H4yz_Px_Pz_S_bc+ABY*I_ERI_G3yz_Px_Pz_S_bc;
  Double I_ERI_G2y2z_Dxy_Pz_S_bc = I_ERI_H3y2z_Px_Pz_S_bc+ABY*I_ERI_G2y2z_Px_Pz_S_bc;
  Double I_ERI_Gy3z_Dxy_Pz_S_bc = I_ERI_H2y3z_Px_Pz_S_bc+ABY*I_ERI_Gy3z_Px_Pz_S_bc;
  Double I_ERI_G4z_Dxy_Pz_S_bc = I_ERI_Hy4z_Px_Pz_S_bc+ABY*I_ERI_G4z_Px_Pz_S_bc;
  Double I_ERI_G4x_D2y_Pz_S_bc = I_ERI_H4xy_Py_Pz_S_bc+ABY*I_ERI_G4x_Py_Pz_S_bc;
  Double I_ERI_G3xy_D2y_Pz_S_bc = I_ERI_H3x2y_Py_Pz_S_bc+ABY*I_ERI_G3xy_Py_Pz_S_bc;
  Double I_ERI_G3xz_D2y_Pz_S_bc = I_ERI_H3xyz_Py_Pz_S_bc+ABY*I_ERI_G3xz_Py_Pz_S_bc;
  Double I_ERI_G2x2y_D2y_Pz_S_bc = I_ERI_H2x3y_Py_Pz_S_bc+ABY*I_ERI_G2x2y_Py_Pz_S_bc;
  Double I_ERI_G2xyz_D2y_Pz_S_bc = I_ERI_H2x2yz_Py_Pz_S_bc+ABY*I_ERI_G2xyz_Py_Pz_S_bc;
  Double I_ERI_G2x2z_D2y_Pz_S_bc = I_ERI_H2xy2z_Py_Pz_S_bc+ABY*I_ERI_G2x2z_Py_Pz_S_bc;
  Double I_ERI_Gx3y_D2y_Pz_S_bc = I_ERI_Hx4y_Py_Pz_S_bc+ABY*I_ERI_Gx3y_Py_Pz_S_bc;
  Double I_ERI_Gx2yz_D2y_Pz_S_bc = I_ERI_Hx3yz_Py_Pz_S_bc+ABY*I_ERI_Gx2yz_Py_Pz_S_bc;
  Double I_ERI_Gxy2z_D2y_Pz_S_bc = I_ERI_Hx2y2z_Py_Pz_S_bc+ABY*I_ERI_Gxy2z_Py_Pz_S_bc;
  Double I_ERI_Gx3z_D2y_Pz_S_bc = I_ERI_Hxy3z_Py_Pz_S_bc+ABY*I_ERI_Gx3z_Py_Pz_S_bc;
  Double I_ERI_G4y_D2y_Pz_S_bc = I_ERI_H5y_Py_Pz_S_bc+ABY*I_ERI_G4y_Py_Pz_S_bc;
  Double I_ERI_G3yz_D2y_Pz_S_bc = I_ERI_H4yz_Py_Pz_S_bc+ABY*I_ERI_G3yz_Py_Pz_S_bc;
  Double I_ERI_G2y2z_D2y_Pz_S_bc = I_ERI_H3y2z_Py_Pz_S_bc+ABY*I_ERI_G2y2z_Py_Pz_S_bc;
  Double I_ERI_Gy3z_D2y_Pz_S_bc = I_ERI_H2y3z_Py_Pz_S_bc+ABY*I_ERI_Gy3z_Py_Pz_S_bc;
  Double I_ERI_G4z_D2y_Pz_S_bc = I_ERI_Hy4z_Py_Pz_S_bc+ABY*I_ERI_G4z_Py_Pz_S_bc;
  Double I_ERI_G4x_D2z_Pz_S_bc = I_ERI_H4xz_Pz_Pz_S_bc+ABZ*I_ERI_G4x_Pz_Pz_S_bc;
  Double I_ERI_G3xy_D2z_Pz_S_bc = I_ERI_H3xyz_Pz_Pz_S_bc+ABZ*I_ERI_G3xy_Pz_Pz_S_bc;
  Double I_ERI_G3xz_D2z_Pz_S_bc = I_ERI_H3x2z_Pz_Pz_S_bc+ABZ*I_ERI_G3xz_Pz_Pz_S_bc;
  Double I_ERI_G2x2y_D2z_Pz_S_bc = I_ERI_H2x2yz_Pz_Pz_S_bc+ABZ*I_ERI_G2x2y_Pz_Pz_S_bc;
  Double I_ERI_G2xyz_D2z_Pz_S_bc = I_ERI_H2xy2z_Pz_Pz_S_bc+ABZ*I_ERI_G2xyz_Pz_Pz_S_bc;
  Double I_ERI_G2x2z_D2z_Pz_S_bc = I_ERI_H2x3z_Pz_Pz_S_bc+ABZ*I_ERI_G2x2z_Pz_Pz_S_bc;
  Double I_ERI_Gx3y_D2z_Pz_S_bc = I_ERI_Hx3yz_Pz_Pz_S_bc+ABZ*I_ERI_Gx3y_Pz_Pz_S_bc;
  Double I_ERI_Gx2yz_D2z_Pz_S_bc = I_ERI_Hx2y2z_Pz_Pz_S_bc+ABZ*I_ERI_Gx2yz_Pz_Pz_S_bc;
  Double I_ERI_Gxy2z_D2z_Pz_S_bc = I_ERI_Hxy3z_Pz_Pz_S_bc+ABZ*I_ERI_Gxy2z_Pz_Pz_S_bc;
  Double I_ERI_Gx3z_D2z_Pz_S_bc = I_ERI_Hx4z_Pz_Pz_S_bc+ABZ*I_ERI_Gx3z_Pz_Pz_S_bc;
  Double I_ERI_G4y_D2z_Pz_S_bc = I_ERI_H4yz_Pz_Pz_S_bc+ABZ*I_ERI_G4y_Pz_Pz_S_bc;
  Double I_ERI_G3yz_D2z_Pz_S_bc = I_ERI_H3y2z_Pz_Pz_S_bc+ABZ*I_ERI_G3yz_Pz_Pz_S_bc;
  Double I_ERI_G2y2z_D2z_Pz_S_bc = I_ERI_H2y3z_Pz_Pz_S_bc+ABZ*I_ERI_G2y2z_Pz_Pz_S_bc;
  Double I_ERI_Gy3z_D2z_Pz_S_bc = I_ERI_Hy4z_Pz_Pz_S_bc+ABZ*I_ERI_Gy3z_Pz_Pz_S_bc;
  Double I_ERI_G4z_D2z_Pz_S_bc = I_ERI_H5z_Pz_Pz_S_bc+ABZ*I_ERI_G4z_Pz_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_F_P_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_D_P_S_bc
   ************************************************************/
  Double I_ERI_F3x_F3x_Px_S_bc = I_ERI_G4x_D2x_Px_S_bc+ABX*I_ERI_F3x_D2x_Px_S_bc;
  Double I_ERI_F2xy_F3x_Px_S_bc = I_ERI_G3xy_D2x_Px_S_bc+ABX*I_ERI_F2xy_D2x_Px_S_bc;
  Double I_ERI_F2xz_F3x_Px_S_bc = I_ERI_G3xz_D2x_Px_S_bc+ABX*I_ERI_F2xz_D2x_Px_S_bc;
  Double I_ERI_Fx2y_F3x_Px_S_bc = I_ERI_G2x2y_D2x_Px_S_bc+ABX*I_ERI_Fx2y_D2x_Px_S_bc;
  Double I_ERI_Fxyz_F3x_Px_S_bc = I_ERI_G2xyz_D2x_Px_S_bc+ABX*I_ERI_Fxyz_D2x_Px_S_bc;
  Double I_ERI_Fx2z_F3x_Px_S_bc = I_ERI_G2x2z_D2x_Px_S_bc+ABX*I_ERI_Fx2z_D2x_Px_S_bc;
  Double I_ERI_F3y_F3x_Px_S_bc = I_ERI_Gx3y_D2x_Px_S_bc+ABX*I_ERI_F3y_D2x_Px_S_bc;
  Double I_ERI_F2yz_F3x_Px_S_bc = I_ERI_Gx2yz_D2x_Px_S_bc+ABX*I_ERI_F2yz_D2x_Px_S_bc;
  Double I_ERI_Fy2z_F3x_Px_S_bc = I_ERI_Gxy2z_D2x_Px_S_bc+ABX*I_ERI_Fy2z_D2x_Px_S_bc;
  Double I_ERI_F3z_F3x_Px_S_bc = I_ERI_Gx3z_D2x_Px_S_bc+ABX*I_ERI_F3z_D2x_Px_S_bc;
  Double I_ERI_F3x_F2xy_Px_S_bc = I_ERI_G3xy_D2x_Px_S_bc+ABY*I_ERI_F3x_D2x_Px_S_bc;
  Double I_ERI_F2xy_F2xy_Px_S_bc = I_ERI_G2x2y_D2x_Px_S_bc+ABY*I_ERI_F2xy_D2x_Px_S_bc;
  Double I_ERI_F2xz_F2xy_Px_S_bc = I_ERI_G2xyz_D2x_Px_S_bc+ABY*I_ERI_F2xz_D2x_Px_S_bc;
  Double I_ERI_Fx2y_F2xy_Px_S_bc = I_ERI_Gx3y_D2x_Px_S_bc+ABY*I_ERI_Fx2y_D2x_Px_S_bc;
  Double I_ERI_Fxyz_F2xy_Px_S_bc = I_ERI_Gx2yz_D2x_Px_S_bc+ABY*I_ERI_Fxyz_D2x_Px_S_bc;
  Double I_ERI_Fx2z_F2xy_Px_S_bc = I_ERI_Gxy2z_D2x_Px_S_bc+ABY*I_ERI_Fx2z_D2x_Px_S_bc;
  Double I_ERI_F3y_F2xy_Px_S_bc = I_ERI_G4y_D2x_Px_S_bc+ABY*I_ERI_F3y_D2x_Px_S_bc;
  Double I_ERI_F2yz_F2xy_Px_S_bc = I_ERI_G3yz_D2x_Px_S_bc+ABY*I_ERI_F2yz_D2x_Px_S_bc;
  Double I_ERI_Fy2z_F2xy_Px_S_bc = I_ERI_G2y2z_D2x_Px_S_bc+ABY*I_ERI_Fy2z_D2x_Px_S_bc;
  Double I_ERI_F3z_F2xy_Px_S_bc = I_ERI_Gy3z_D2x_Px_S_bc+ABY*I_ERI_F3z_D2x_Px_S_bc;
  Double I_ERI_F3x_F2xz_Px_S_bc = I_ERI_G3xz_D2x_Px_S_bc+ABZ*I_ERI_F3x_D2x_Px_S_bc;
  Double I_ERI_F2xy_F2xz_Px_S_bc = I_ERI_G2xyz_D2x_Px_S_bc+ABZ*I_ERI_F2xy_D2x_Px_S_bc;
  Double I_ERI_F2xz_F2xz_Px_S_bc = I_ERI_G2x2z_D2x_Px_S_bc+ABZ*I_ERI_F2xz_D2x_Px_S_bc;
  Double I_ERI_Fx2y_F2xz_Px_S_bc = I_ERI_Gx2yz_D2x_Px_S_bc+ABZ*I_ERI_Fx2y_D2x_Px_S_bc;
  Double I_ERI_Fxyz_F2xz_Px_S_bc = I_ERI_Gxy2z_D2x_Px_S_bc+ABZ*I_ERI_Fxyz_D2x_Px_S_bc;
  Double I_ERI_Fx2z_F2xz_Px_S_bc = I_ERI_Gx3z_D2x_Px_S_bc+ABZ*I_ERI_Fx2z_D2x_Px_S_bc;
  Double I_ERI_F3y_F2xz_Px_S_bc = I_ERI_G3yz_D2x_Px_S_bc+ABZ*I_ERI_F3y_D2x_Px_S_bc;
  Double I_ERI_F2yz_F2xz_Px_S_bc = I_ERI_G2y2z_D2x_Px_S_bc+ABZ*I_ERI_F2yz_D2x_Px_S_bc;
  Double I_ERI_Fy2z_F2xz_Px_S_bc = I_ERI_Gy3z_D2x_Px_S_bc+ABZ*I_ERI_Fy2z_D2x_Px_S_bc;
  Double I_ERI_F3z_F2xz_Px_S_bc = I_ERI_G4z_D2x_Px_S_bc+ABZ*I_ERI_F3z_D2x_Px_S_bc;
  Double I_ERI_F3x_Fx2y_Px_S_bc = I_ERI_G4x_D2y_Px_S_bc+ABX*I_ERI_F3x_D2y_Px_S_bc;
  Double I_ERI_F2xy_Fx2y_Px_S_bc = I_ERI_G3xy_D2y_Px_S_bc+ABX*I_ERI_F2xy_D2y_Px_S_bc;
  Double I_ERI_F2xz_Fx2y_Px_S_bc = I_ERI_G3xz_D2y_Px_S_bc+ABX*I_ERI_F2xz_D2y_Px_S_bc;
  Double I_ERI_Fx2y_Fx2y_Px_S_bc = I_ERI_G2x2y_D2y_Px_S_bc+ABX*I_ERI_Fx2y_D2y_Px_S_bc;
  Double I_ERI_Fxyz_Fx2y_Px_S_bc = I_ERI_G2xyz_D2y_Px_S_bc+ABX*I_ERI_Fxyz_D2y_Px_S_bc;
  Double I_ERI_Fx2z_Fx2y_Px_S_bc = I_ERI_G2x2z_D2y_Px_S_bc+ABX*I_ERI_Fx2z_D2y_Px_S_bc;
  Double I_ERI_F3y_Fx2y_Px_S_bc = I_ERI_Gx3y_D2y_Px_S_bc+ABX*I_ERI_F3y_D2y_Px_S_bc;
  Double I_ERI_F2yz_Fx2y_Px_S_bc = I_ERI_Gx2yz_D2y_Px_S_bc+ABX*I_ERI_F2yz_D2y_Px_S_bc;
  Double I_ERI_Fy2z_Fx2y_Px_S_bc = I_ERI_Gxy2z_D2y_Px_S_bc+ABX*I_ERI_Fy2z_D2y_Px_S_bc;
  Double I_ERI_F3z_Fx2y_Px_S_bc = I_ERI_Gx3z_D2y_Px_S_bc+ABX*I_ERI_F3z_D2y_Px_S_bc;
  Double I_ERI_F3x_Fxyz_Px_S_bc = I_ERI_G3xz_Dxy_Px_S_bc+ABZ*I_ERI_F3x_Dxy_Px_S_bc;
  Double I_ERI_F2xy_Fxyz_Px_S_bc = I_ERI_G2xyz_Dxy_Px_S_bc+ABZ*I_ERI_F2xy_Dxy_Px_S_bc;
  Double I_ERI_F2xz_Fxyz_Px_S_bc = I_ERI_G2x2z_Dxy_Px_S_bc+ABZ*I_ERI_F2xz_Dxy_Px_S_bc;
  Double I_ERI_Fx2y_Fxyz_Px_S_bc = I_ERI_Gx2yz_Dxy_Px_S_bc+ABZ*I_ERI_Fx2y_Dxy_Px_S_bc;
  Double I_ERI_Fxyz_Fxyz_Px_S_bc = I_ERI_Gxy2z_Dxy_Px_S_bc+ABZ*I_ERI_Fxyz_Dxy_Px_S_bc;
  Double I_ERI_Fx2z_Fxyz_Px_S_bc = I_ERI_Gx3z_Dxy_Px_S_bc+ABZ*I_ERI_Fx2z_Dxy_Px_S_bc;
  Double I_ERI_F3y_Fxyz_Px_S_bc = I_ERI_G3yz_Dxy_Px_S_bc+ABZ*I_ERI_F3y_Dxy_Px_S_bc;
  Double I_ERI_F2yz_Fxyz_Px_S_bc = I_ERI_G2y2z_Dxy_Px_S_bc+ABZ*I_ERI_F2yz_Dxy_Px_S_bc;
  Double I_ERI_Fy2z_Fxyz_Px_S_bc = I_ERI_Gy3z_Dxy_Px_S_bc+ABZ*I_ERI_Fy2z_Dxy_Px_S_bc;
  Double I_ERI_F3z_Fxyz_Px_S_bc = I_ERI_G4z_Dxy_Px_S_bc+ABZ*I_ERI_F3z_Dxy_Px_S_bc;
  Double I_ERI_F3x_Fx2z_Px_S_bc = I_ERI_G4x_D2z_Px_S_bc+ABX*I_ERI_F3x_D2z_Px_S_bc;
  Double I_ERI_F2xy_Fx2z_Px_S_bc = I_ERI_G3xy_D2z_Px_S_bc+ABX*I_ERI_F2xy_D2z_Px_S_bc;
  Double I_ERI_F2xz_Fx2z_Px_S_bc = I_ERI_G3xz_D2z_Px_S_bc+ABX*I_ERI_F2xz_D2z_Px_S_bc;
  Double I_ERI_Fx2y_Fx2z_Px_S_bc = I_ERI_G2x2y_D2z_Px_S_bc+ABX*I_ERI_Fx2y_D2z_Px_S_bc;
  Double I_ERI_Fxyz_Fx2z_Px_S_bc = I_ERI_G2xyz_D2z_Px_S_bc+ABX*I_ERI_Fxyz_D2z_Px_S_bc;
  Double I_ERI_Fx2z_Fx2z_Px_S_bc = I_ERI_G2x2z_D2z_Px_S_bc+ABX*I_ERI_Fx2z_D2z_Px_S_bc;
  Double I_ERI_F3y_Fx2z_Px_S_bc = I_ERI_Gx3y_D2z_Px_S_bc+ABX*I_ERI_F3y_D2z_Px_S_bc;
  Double I_ERI_F2yz_Fx2z_Px_S_bc = I_ERI_Gx2yz_D2z_Px_S_bc+ABX*I_ERI_F2yz_D2z_Px_S_bc;
  Double I_ERI_Fy2z_Fx2z_Px_S_bc = I_ERI_Gxy2z_D2z_Px_S_bc+ABX*I_ERI_Fy2z_D2z_Px_S_bc;
  Double I_ERI_F3z_Fx2z_Px_S_bc = I_ERI_Gx3z_D2z_Px_S_bc+ABX*I_ERI_F3z_D2z_Px_S_bc;
  Double I_ERI_F3x_F3y_Px_S_bc = I_ERI_G3xy_D2y_Px_S_bc+ABY*I_ERI_F3x_D2y_Px_S_bc;
  Double I_ERI_F2xy_F3y_Px_S_bc = I_ERI_G2x2y_D2y_Px_S_bc+ABY*I_ERI_F2xy_D2y_Px_S_bc;
  Double I_ERI_F2xz_F3y_Px_S_bc = I_ERI_G2xyz_D2y_Px_S_bc+ABY*I_ERI_F2xz_D2y_Px_S_bc;
  Double I_ERI_Fx2y_F3y_Px_S_bc = I_ERI_Gx3y_D2y_Px_S_bc+ABY*I_ERI_Fx2y_D2y_Px_S_bc;
  Double I_ERI_Fxyz_F3y_Px_S_bc = I_ERI_Gx2yz_D2y_Px_S_bc+ABY*I_ERI_Fxyz_D2y_Px_S_bc;
  Double I_ERI_Fx2z_F3y_Px_S_bc = I_ERI_Gxy2z_D2y_Px_S_bc+ABY*I_ERI_Fx2z_D2y_Px_S_bc;
  Double I_ERI_F3y_F3y_Px_S_bc = I_ERI_G4y_D2y_Px_S_bc+ABY*I_ERI_F3y_D2y_Px_S_bc;
  Double I_ERI_F2yz_F3y_Px_S_bc = I_ERI_G3yz_D2y_Px_S_bc+ABY*I_ERI_F2yz_D2y_Px_S_bc;
  Double I_ERI_Fy2z_F3y_Px_S_bc = I_ERI_G2y2z_D2y_Px_S_bc+ABY*I_ERI_Fy2z_D2y_Px_S_bc;
  Double I_ERI_F3z_F3y_Px_S_bc = I_ERI_Gy3z_D2y_Px_S_bc+ABY*I_ERI_F3z_D2y_Px_S_bc;
  Double I_ERI_F3x_F2yz_Px_S_bc = I_ERI_G3xz_D2y_Px_S_bc+ABZ*I_ERI_F3x_D2y_Px_S_bc;
  Double I_ERI_F2xy_F2yz_Px_S_bc = I_ERI_G2xyz_D2y_Px_S_bc+ABZ*I_ERI_F2xy_D2y_Px_S_bc;
  Double I_ERI_F2xz_F2yz_Px_S_bc = I_ERI_G2x2z_D2y_Px_S_bc+ABZ*I_ERI_F2xz_D2y_Px_S_bc;
  Double I_ERI_Fx2y_F2yz_Px_S_bc = I_ERI_Gx2yz_D2y_Px_S_bc+ABZ*I_ERI_Fx2y_D2y_Px_S_bc;
  Double I_ERI_Fxyz_F2yz_Px_S_bc = I_ERI_Gxy2z_D2y_Px_S_bc+ABZ*I_ERI_Fxyz_D2y_Px_S_bc;
  Double I_ERI_Fx2z_F2yz_Px_S_bc = I_ERI_Gx3z_D2y_Px_S_bc+ABZ*I_ERI_Fx2z_D2y_Px_S_bc;
  Double I_ERI_F3y_F2yz_Px_S_bc = I_ERI_G3yz_D2y_Px_S_bc+ABZ*I_ERI_F3y_D2y_Px_S_bc;
  Double I_ERI_F2yz_F2yz_Px_S_bc = I_ERI_G2y2z_D2y_Px_S_bc+ABZ*I_ERI_F2yz_D2y_Px_S_bc;
  Double I_ERI_Fy2z_F2yz_Px_S_bc = I_ERI_Gy3z_D2y_Px_S_bc+ABZ*I_ERI_Fy2z_D2y_Px_S_bc;
  Double I_ERI_F3z_F2yz_Px_S_bc = I_ERI_G4z_D2y_Px_S_bc+ABZ*I_ERI_F3z_D2y_Px_S_bc;
  Double I_ERI_F3x_Fy2z_Px_S_bc = I_ERI_G3xy_D2z_Px_S_bc+ABY*I_ERI_F3x_D2z_Px_S_bc;
  Double I_ERI_F2xy_Fy2z_Px_S_bc = I_ERI_G2x2y_D2z_Px_S_bc+ABY*I_ERI_F2xy_D2z_Px_S_bc;
  Double I_ERI_F2xz_Fy2z_Px_S_bc = I_ERI_G2xyz_D2z_Px_S_bc+ABY*I_ERI_F2xz_D2z_Px_S_bc;
  Double I_ERI_Fx2y_Fy2z_Px_S_bc = I_ERI_Gx3y_D2z_Px_S_bc+ABY*I_ERI_Fx2y_D2z_Px_S_bc;
  Double I_ERI_Fxyz_Fy2z_Px_S_bc = I_ERI_Gx2yz_D2z_Px_S_bc+ABY*I_ERI_Fxyz_D2z_Px_S_bc;
  Double I_ERI_Fx2z_Fy2z_Px_S_bc = I_ERI_Gxy2z_D2z_Px_S_bc+ABY*I_ERI_Fx2z_D2z_Px_S_bc;
  Double I_ERI_F3y_Fy2z_Px_S_bc = I_ERI_G4y_D2z_Px_S_bc+ABY*I_ERI_F3y_D2z_Px_S_bc;
  Double I_ERI_F2yz_Fy2z_Px_S_bc = I_ERI_G3yz_D2z_Px_S_bc+ABY*I_ERI_F2yz_D2z_Px_S_bc;
  Double I_ERI_Fy2z_Fy2z_Px_S_bc = I_ERI_G2y2z_D2z_Px_S_bc+ABY*I_ERI_Fy2z_D2z_Px_S_bc;
  Double I_ERI_F3z_Fy2z_Px_S_bc = I_ERI_Gy3z_D2z_Px_S_bc+ABY*I_ERI_F3z_D2z_Px_S_bc;
  Double I_ERI_F3x_F3z_Px_S_bc = I_ERI_G3xz_D2z_Px_S_bc+ABZ*I_ERI_F3x_D2z_Px_S_bc;
  Double I_ERI_F2xy_F3z_Px_S_bc = I_ERI_G2xyz_D2z_Px_S_bc+ABZ*I_ERI_F2xy_D2z_Px_S_bc;
  Double I_ERI_F2xz_F3z_Px_S_bc = I_ERI_G2x2z_D2z_Px_S_bc+ABZ*I_ERI_F2xz_D2z_Px_S_bc;
  Double I_ERI_Fx2y_F3z_Px_S_bc = I_ERI_Gx2yz_D2z_Px_S_bc+ABZ*I_ERI_Fx2y_D2z_Px_S_bc;
  Double I_ERI_Fxyz_F3z_Px_S_bc = I_ERI_Gxy2z_D2z_Px_S_bc+ABZ*I_ERI_Fxyz_D2z_Px_S_bc;
  Double I_ERI_Fx2z_F3z_Px_S_bc = I_ERI_Gx3z_D2z_Px_S_bc+ABZ*I_ERI_Fx2z_D2z_Px_S_bc;
  Double I_ERI_F3y_F3z_Px_S_bc = I_ERI_G3yz_D2z_Px_S_bc+ABZ*I_ERI_F3y_D2z_Px_S_bc;
  Double I_ERI_F2yz_F3z_Px_S_bc = I_ERI_G2y2z_D2z_Px_S_bc+ABZ*I_ERI_F2yz_D2z_Px_S_bc;
  Double I_ERI_Fy2z_F3z_Px_S_bc = I_ERI_Gy3z_D2z_Px_S_bc+ABZ*I_ERI_Fy2z_D2z_Px_S_bc;
  Double I_ERI_F3z_F3z_Px_S_bc = I_ERI_G4z_D2z_Px_S_bc+ABZ*I_ERI_F3z_D2z_Px_S_bc;
  Double I_ERI_F3x_F3x_Py_S_bc = I_ERI_G4x_D2x_Py_S_bc+ABX*I_ERI_F3x_D2x_Py_S_bc;
  Double I_ERI_F2xy_F3x_Py_S_bc = I_ERI_G3xy_D2x_Py_S_bc+ABX*I_ERI_F2xy_D2x_Py_S_bc;
  Double I_ERI_F2xz_F3x_Py_S_bc = I_ERI_G3xz_D2x_Py_S_bc+ABX*I_ERI_F2xz_D2x_Py_S_bc;
  Double I_ERI_Fx2y_F3x_Py_S_bc = I_ERI_G2x2y_D2x_Py_S_bc+ABX*I_ERI_Fx2y_D2x_Py_S_bc;
  Double I_ERI_Fxyz_F3x_Py_S_bc = I_ERI_G2xyz_D2x_Py_S_bc+ABX*I_ERI_Fxyz_D2x_Py_S_bc;
  Double I_ERI_Fx2z_F3x_Py_S_bc = I_ERI_G2x2z_D2x_Py_S_bc+ABX*I_ERI_Fx2z_D2x_Py_S_bc;
  Double I_ERI_F3y_F3x_Py_S_bc = I_ERI_Gx3y_D2x_Py_S_bc+ABX*I_ERI_F3y_D2x_Py_S_bc;
  Double I_ERI_F2yz_F3x_Py_S_bc = I_ERI_Gx2yz_D2x_Py_S_bc+ABX*I_ERI_F2yz_D2x_Py_S_bc;
  Double I_ERI_Fy2z_F3x_Py_S_bc = I_ERI_Gxy2z_D2x_Py_S_bc+ABX*I_ERI_Fy2z_D2x_Py_S_bc;
  Double I_ERI_F3z_F3x_Py_S_bc = I_ERI_Gx3z_D2x_Py_S_bc+ABX*I_ERI_F3z_D2x_Py_S_bc;
  Double I_ERI_F3x_F2xy_Py_S_bc = I_ERI_G3xy_D2x_Py_S_bc+ABY*I_ERI_F3x_D2x_Py_S_bc;
  Double I_ERI_F2xy_F2xy_Py_S_bc = I_ERI_G2x2y_D2x_Py_S_bc+ABY*I_ERI_F2xy_D2x_Py_S_bc;
  Double I_ERI_F2xz_F2xy_Py_S_bc = I_ERI_G2xyz_D2x_Py_S_bc+ABY*I_ERI_F2xz_D2x_Py_S_bc;
  Double I_ERI_Fx2y_F2xy_Py_S_bc = I_ERI_Gx3y_D2x_Py_S_bc+ABY*I_ERI_Fx2y_D2x_Py_S_bc;
  Double I_ERI_Fxyz_F2xy_Py_S_bc = I_ERI_Gx2yz_D2x_Py_S_bc+ABY*I_ERI_Fxyz_D2x_Py_S_bc;
  Double I_ERI_Fx2z_F2xy_Py_S_bc = I_ERI_Gxy2z_D2x_Py_S_bc+ABY*I_ERI_Fx2z_D2x_Py_S_bc;
  Double I_ERI_F3y_F2xy_Py_S_bc = I_ERI_G4y_D2x_Py_S_bc+ABY*I_ERI_F3y_D2x_Py_S_bc;
  Double I_ERI_F2yz_F2xy_Py_S_bc = I_ERI_G3yz_D2x_Py_S_bc+ABY*I_ERI_F2yz_D2x_Py_S_bc;
  Double I_ERI_Fy2z_F2xy_Py_S_bc = I_ERI_G2y2z_D2x_Py_S_bc+ABY*I_ERI_Fy2z_D2x_Py_S_bc;
  Double I_ERI_F3z_F2xy_Py_S_bc = I_ERI_Gy3z_D2x_Py_S_bc+ABY*I_ERI_F3z_D2x_Py_S_bc;
  Double I_ERI_F3x_F2xz_Py_S_bc = I_ERI_G3xz_D2x_Py_S_bc+ABZ*I_ERI_F3x_D2x_Py_S_bc;
  Double I_ERI_F2xy_F2xz_Py_S_bc = I_ERI_G2xyz_D2x_Py_S_bc+ABZ*I_ERI_F2xy_D2x_Py_S_bc;
  Double I_ERI_F2xz_F2xz_Py_S_bc = I_ERI_G2x2z_D2x_Py_S_bc+ABZ*I_ERI_F2xz_D2x_Py_S_bc;
  Double I_ERI_Fx2y_F2xz_Py_S_bc = I_ERI_Gx2yz_D2x_Py_S_bc+ABZ*I_ERI_Fx2y_D2x_Py_S_bc;
  Double I_ERI_Fxyz_F2xz_Py_S_bc = I_ERI_Gxy2z_D2x_Py_S_bc+ABZ*I_ERI_Fxyz_D2x_Py_S_bc;
  Double I_ERI_Fx2z_F2xz_Py_S_bc = I_ERI_Gx3z_D2x_Py_S_bc+ABZ*I_ERI_Fx2z_D2x_Py_S_bc;
  Double I_ERI_F3y_F2xz_Py_S_bc = I_ERI_G3yz_D2x_Py_S_bc+ABZ*I_ERI_F3y_D2x_Py_S_bc;
  Double I_ERI_F2yz_F2xz_Py_S_bc = I_ERI_G2y2z_D2x_Py_S_bc+ABZ*I_ERI_F2yz_D2x_Py_S_bc;
  Double I_ERI_Fy2z_F2xz_Py_S_bc = I_ERI_Gy3z_D2x_Py_S_bc+ABZ*I_ERI_Fy2z_D2x_Py_S_bc;
  Double I_ERI_F3z_F2xz_Py_S_bc = I_ERI_G4z_D2x_Py_S_bc+ABZ*I_ERI_F3z_D2x_Py_S_bc;
  Double I_ERI_F3x_Fx2y_Py_S_bc = I_ERI_G4x_D2y_Py_S_bc+ABX*I_ERI_F3x_D2y_Py_S_bc;
  Double I_ERI_F2xy_Fx2y_Py_S_bc = I_ERI_G3xy_D2y_Py_S_bc+ABX*I_ERI_F2xy_D2y_Py_S_bc;
  Double I_ERI_F2xz_Fx2y_Py_S_bc = I_ERI_G3xz_D2y_Py_S_bc+ABX*I_ERI_F2xz_D2y_Py_S_bc;
  Double I_ERI_Fx2y_Fx2y_Py_S_bc = I_ERI_G2x2y_D2y_Py_S_bc+ABX*I_ERI_Fx2y_D2y_Py_S_bc;
  Double I_ERI_Fxyz_Fx2y_Py_S_bc = I_ERI_G2xyz_D2y_Py_S_bc+ABX*I_ERI_Fxyz_D2y_Py_S_bc;
  Double I_ERI_Fx2z_Fx2y_Py_S_bc = I_ERI_G2x2z_D2y_Py_S_bc+ABX*I_ERI_Fx2z_D2y_Py_S_bc;
  Double I_ERI_F3y_Fx2y_Py_S_bc = I_ERI_Gx3y_D2y_Py_S_bc+ABX*I_ERI_F3y_D2y_Py_S_bc;
  Double I_ERI_F2yz_Fx2y_Py_S_bc = I_ERI_Gx2yz_D2y_Py_S_bc+ABX*I_ERI_F2yz_D2y_Py_S_bc;
  Double I_ERI_Fy2z_Fx2y_Py_S_bc = I_ERI_Gxy2z_D2y_Py_S_bc+ABX*I_ERI_Fy2z_D2y_Py_S_bc;
  Double I_ERI_F3z_Fx2y_Py_S_bc = I_ERI_Gx3z_D2y_Py_S_bc+ABX*I_ERI_F3z_D2y_Py_S_bc;
  Double I_ERI_F3x_Fxyz_Py_S_bc = I_ERI_G3xz_Dxy_Py_S_bc+ABZ*I_ERI_F3x_Dxy_Py_S_bc;
  Double I_ERI_F2xy_Fxyz_Py_S_bc = I_ERI_G2xyz_Dxy_Py_S_bc+ABZ*I_ERI_F2xy_Dxy_Py_S_bc;
  Double I_ERI_F2xz_Fxyz_Py_S_bc = I_ERI_G2x2z_Dxy_Py_S_bc+ABZ*I_ERI_F2xz_Dxy_Py_S_bc;
  Double I_ERI_Fx2y_Fxyz_Py_S_bc = I_ERI_Gx2yz_Dxy_Py_S_bc+ABZ*I_ERI_Fx2y_Dxy_Py_S_bc;
  Double I_ERI_Fxyz_Fxyz_Py_S_bc = I_ERI_Gxy2z_Dxy_Py_S_bc+ABZ*I_ERI_Fxyz_Dxy_Py_S_bc;
  Double I_ERI_Fx2z_Fxyz_Py_S_bc = I_ERI_Gx3z_Dxy_Py_S_bc+ABZ*I_ERI_Fx2z_Dxy_Py_S_bc;
  Double I_ERI_F3y_Fxyz_Py_S_bc = I_ERI_G3yz_Dxy_Py_S_bc+ABZ*I_ERI_F3y_Dxy_Py_S_bc;
  Double I_ERI_F2yz_Fxyz_Py_S_bc = I_ERI_G2y2z_Dxy_Py_S_bc+ABZ*I_ERI_F2yz_Dxy_Py_S_bc;
  Double I_ERI_Fy2z_Fxyz_Py_S_bc = I_ERI_Gy3z_Dxy_Py_S_bc+ABZ*I_ERI_Fy2z_Dxy_Py_S_bc;
  Double I_ERI_F3z_Fxyz_Py_S_bc = I_ERI_G4z_Dxy_Py_S_bc+ABZ*I_ERI_F3z_Dxy_Py_S_bc;
  Double I_ERI_F3x_Fx2z_Py_S_bc = I_ERI_G4x_D2z_Py_S_bc+ABX*I_ERI_F3x_D2z_Py_S_bc;
  Double I_ERI_F2xy_Fx2z_Py_S_bc = I_ERI_G3xy_D2z_Py_S_bc+ABX*I_ERI_F2xy_D2z_Py_S_bc;
  Double I_ERI_F2xz_Fx2z_Py_S_bc = I_ERI_G3xz_D2z_Py_S_bc+ABX*I_ERI_F2xz_D2z_Py_S_bc;
  Double I_ERI_Fx2y_Fx2z_Py_S_bc = I_ERI_G2x2y_D2z_Py_S_bc+ABX*I_ERI_Fx2y_D2z_Py_S_bc;
  Double I_ERI_Fxyz_Fx2z_Py_S_bc = I_ERI_G2xyz_D2z_Py_S_bc+ABX*I_ERI_Fxyz_D2z_Py_S_bc;
  Double I_ERI_Fx2z_Fx2z_Py_S_bc = I_ERI_G2x2z_D2z_Py_S_bc+ABX*I_ERI_Fx2z_D2z_Py_S_bc;
  Double I_ERI_F3y_Fx2z_Py_S_bc = I_ERI_Gx3y_D2z_Py_S_bc+ABX*I_ERI_F3y_D2z_Py_S_bc;
  Double I_ERI_F2yz_Fx2z_Py_S_bc = I_ERI_Gx2yz_D2z_Py_S_bc+ABX*I_ERI_F2yz_D2z_Py_S_bc;
  Double I_ERI_Fy2z_Fx2z_Py_S_bc = I_ERI_Gxy2z_D2z_Py_S_bc+ABX*I_ERI_Fy2z_D2z_Py_S_bc;
  Double I_ERI_F3z_Fx2z_Py_S_bc = I_ERI_Gx3z_D2z_Py_S_bc+ABX*I_ERI_F3z_D2z_Py_S_bc;
  Double I_ERI_F3x_F3y_Py_S_bc = I_ERI_G3xy_D2y_Py_S_bc+ABY*I_ERI_F3x_D2y_Py_S_bc;
  Double I_ERI_F2xy_F3y_Py_S_bc = I_ERI_G2x2y_D2y_Py_S_bc+ABY*I_ERI_F2xy_D2y_Py_S_bc;
  Double I_ERI_F2xz_F3y_Py_S_bc = I_ERI_G2xyz_D2y_Py_S_bc+ABY*I_ERI_F2xz_D2y_Py_S_bc;
  Double I_ERI_Fx2y_F3y_Py_S_bc = I_ERI_Gx3y_D2y_Py_S_bc+ABY*I_ERI_Fx2y_D2y_Py_S_bc;
  Double I_ERI_Fxyz_F3y_Py_S_bc = I_ERI_Gx2yz_D2y_Py_S_bc+ABY*I_ERI_Fxyz_D2y_Py_S_bc;
  Double I_ERI_Fx2z_F3y_Py_S_bc = I_ERI_Gxy2z_D2y_Py_S_bc+ABY*I_ERI_Fx2z_D2y_Py_S_bc;
  Double I_ERI_F3y_F3y_Py_S_bc = I_ERI_G4y_D2y_Py_S_bc+ABY*I_ERI_F3y_D2y_Py_S_bc;
  Double I_ERI_F2yz_F3y_Py_S_bc = I_ERI_G3yz_D2y_Py_S_bc+ABY*I_ERI_F2yz_D2y_Py_S_bc;
  Double I_ERI_Fy2z_F3y_Py_S_bc = I_ERI_G2y2z_D2y_Py_S_bc+ABY*I_ERI_Fy2z_D2y_Py_S_bc;
  Double I_ERI_F3z_F3y_Py_S_bc = I_ERI_Gy3z_D2y_Py_S_bc+ABY*I_ERI_F3z_D2y_Py_S_bc;
  Double I_ERI_F3x_F2yz_Py_S_bc = I_ERI_G3xz_D2y_Py_S_bc+ABZ*I_ERI_F3x_D2y_Py_S_bc;
  Double I_ERI_F2xy_F2yz_Py_S_bc = I_ERI_G2xyz_D2y_Py_S_bc+ABZ*I_ERI_F2xy_D2y_Py_S_bc;
  Double I_ERI_F2xz_F2yz_Py_S_bc = I_ERI_G2x2z_D2y_Py_S_bc+ABZ*I_ERI_F2xz_D2y_Py_S_bc;
  Double I_ERI_Fx2y_F2yz_Py_S_bc = I_ERI_Gx2yz_D2y_Py_S_bc+ABZ*I_ERI_Fx2y_D2y_Py_S_bc;
  Double I_ERI_Fxyz_F2yz_Py_S_bc = I_ERI_Gxy2z_D2y_Py_S_bc+ABZ*I_ERI_Fxyz_D2y_Py_S_bc;
  Double I_ERI_Fx2z_F2yz_Py_S_bc = I_ERI_Gx3z_D2y_Py_S_bc+ABZ*I_ERI_Fx2z_D2y_Py_S_bc;
  Double I_ERI_F3y_F2yz_Py_S_bc = I_ERI_G3yz_D2y_Py_S_bc+ABZ*I_ERI_F3y_D2y_Py_S_bc;
  Double I_ERI_F2yz_F2yz_Py_S_bc = I_ERI_G2y2z_D2y_Py_S_bc+ABZ*I_ERI_F2yz_D2y_Py_S_bc;
  Double I_ERI_Fy2z_F2yz_Py_S_bc = I_ERI_Gy3z_D2y_Py_S_bc+ABZ*I_ERI_Fy2z_D2y_Py_S_bc;
  Double I_ERI_F3z_F2yz_Py_S_bc = I_ERI_G4z_D2y_Py_S_bc+ABZ*I_ERI_F3z_D2y_Py_S_bc;
  Double I_ERI_F3x_Fy2z_Py_S_bc = I_ERI_G3xy_D2z_Py_S_bc+ABY*I_ERI_F3x_D2z_Py_S_bc;
  Double I_ERI_F2xy_Fy2z_Py_S_bc = I_ERI_G2x2y_D2z_Py_S_bc+ABY*I_ERI_F2xy_D2z_Py_S_bc;
  Double I_ERI_F2xz_Fy2z_Py_S_bc = I_ERI_G2xyz_D2z_Py_S_bc+ABY*I_ERI_F2xz_D2z_Py_S_bc;
  Double I_ERI_Fx2y_Fy2z_Py_S_bc = I_ERI_Gx3y_D2z_Py_S_bc+ABY*I_ERI_Fx2y_D2z_Py_S_bc;
  Double I_ERI_Fxyz_Fy2z_Py_S_bc = I_ERI_Gx2yz_D2z_Py_S_bc+ABY*I_ERI_Fxyz_D2z_Py_S_bc;
  Double I_ERI_Fx2z_Fy2z_Py_S_bc = I_ERI_Gxy2z_D2z_Py_S_bc+ABY*I_ERI_Fx2z_D2z_Py_S_bc;
  Double I_ERI_F3y_Fy2z_Py_S_bc = I_ERI_G4y_D2z_Py_S_bc+ABY*I_ERI_F3y_D2z_Py_S_bc;
  Double I_ERI_F2yz_Fy2z_Py_S_bc = I_ERI_G3yz_D2z_Py_S_bc+ABY*I_ERI_F2yz_D2z_Py_S_bc;
  Double I_ERI_Fy2z_Fy2z_Py_S_bc = I_ERI_G2y2z_D2z_Py_S_bc+ABY*I_ERI_Fy2z_D2z_Py_S_bc;
  Double I_ERI_F3z_Fy2z_Py_S_bc = I_ERI_Gy3z_D2z_Py_S_bc+ABY*I_ERI_F3z_D2z_Py_S_bc;
  Double I_ERI_F3x_F3z_Py_S_bc = I_ERI_G3xz_D2z_Py_S_bc+ABZ*I_ERI_F3x_D2z_Py_S_bc;
  Double I_ERI_F2xy_F3z_Py_S_bc = I_ERI_G2xyz_D2z_Py_S_bc+ABZ*I_ERI_F2xy_D2z_Py_S_bc;
  Double I_ERI_F2xz_F3z_Py_S_bc = I_ERI_G2x2z_D2z_Py_S_bc+ABZ*I_ERI_F2xz_D2z_Py_S_bc;
  Double I_ERI_Fx2y_F3z_Py_S_bc = I_ERI_Gx2yz_D2z_Py_S_bc+ABZ*I_ERI_Fx2y_D2z_Py_S_bc;
  Double I_ERI_Fxyz_F3z_Py_S_bc = I_ERI_Gxy2z_D2z_Py_S_bc+ABZ*I_ERI_Fxyz_D2z_Py_S_bc;
  Double I_ERI_Fx2z_F3z_Py_S_bc = I_ERI_Gx3z_D2z_Py_S_bc+ABZ*I_ERI_Fx2z_D2z_Py_S_bc;
  Double I_ERI_F3y_F3z_Py_S_bc = I_ERI_G3yz_D2z_Py_S_bc+ABZ*I_ERI_F3y_D2z_Py_S_bc;
  Double I_ERI_F2yz_F3z_Py_S_bc = I_ERI_G2y2z_D2z_Py_S_bc+ABZ*I_ERI_F2yz_D2z_Py_S_bc;
  Double I_ERI_Fy2z_F3z_Py_S_bc = I_ERI_Gy3z_D2z_Py_S_bc+ABZ*I_ERI_Fy2z_D2z_Py_S_bc;
  Double I_ERI_F3z_F3z_Py_S_bc = I_ERI_G4z_D2z_Py_S_bc+ABZ*I_ERI_F3z_D2z_Py_S_bc;
  Double I_ERI_F3x_F3x_Pz_S_bc = I_ERI_G4x_D2x_Pz_S_bc+ABX*I_ERI_F3x_D2x_Pz_S_bc;
  Double I_ERI_F2xy_F3x_Pz_S_bc = I_ERI_G3xy_D2x_Pz_S_bc+ABX*I_ERI_F2xy_D2x_Pz_S_bc;
  Double I_ERI_F2xz_F3x_Pz_S_bc = I_ERI_G3xz_D2x_Pz_S_bc+ABX*I_ERI_F2xz_D2x_Pz_S_bc;
  Double I_ERI_Fx2y_F3x_Pz_S_bc = I_ERI_G2x2y_D2x_Pz_S_bc+ABX*I_ERI_Fx2y_D2x_Pz_S_bc;
  Double I_ERI_Fxyz_F3x_Pz_S_bc = I_ERI_G2xyz_D2x_Pz_S_bc+ABX*I_ERI_Fxyz_D2x_Pz_S_bc;
  Double I_ERI_Fx2z_F3x_Pz_S_bc = I_ERI_G2x2z_D2x_Pz_S_bc+ABX*I_ERI_Fx2z_D2x_Pz_S_bc;
  Double I_ERI_F3y_F3x_Pz_S_bc = I_ERI_Gx3y_D2x_Pz_S_bc+ABX*I_ERI_F3y_D2x_Pz_S_bc;
  Double I_ERI_F2yz_F3x_Pz_S_bc = I_ERI_Gx2yz_D2x_Pz_S_bc+ABX*I_ERI_F2yz_D2x_Pz_S_bc;
  Double I_ERI_Fy2z_F3x_Pz_S_bc = I_ERI_Gxy2z_D2x_Pz_S_bc+ABX*I_ERI_Fy2z_D2x_Pz_S_bc;
  Double I_ERI_F3z_F3x_Pz_S_bc = I_ERI_Gx3z_D2x_Pz_S_bc+ABX*I_ERI_F3z_D2x_Pz_S_bc;
  Double I_ERI_F3x_F2xy_Pz_S_bc = I_ERI_G3xy_D2x_Pz_S_bc+ABY*I_ERI_F3x_D2x_Pz_S_bc;
  Double I_ERI_F2xy_F2xy_Pz_S_bc = I_ERI_G2x2y_D2x_Pz_S_bc+ABY*I_ERI_F2xy_D2x_Pz_S_bc;
  Double I_ERI_F2xz_F2xy_Pz_S_bc = I_ERI_G2xyz_D2x_Pz_S_bc+ABY*I_ERI_F2xz_D2x_Pz_S_bc;
  Double I_ERI_Fx2y_F2xy_Pz_S_bc = I_ERI_Gx3y_D2x_Pz_S_bc+ABY*I_ERI_Fx2y_D2x_Pz_S_bc;
  Double I_ERI_Fxyz_F2xy_Pz_S_bc = I_ERI_Gx2yz_D2x_Pz_S_bc+ABY*I_ERI_Fxyz_D2x_Pz_S_bc;
  Double I_ERI_Fx2z_F2xy_Pz_S_bc = I_ERI_Gxy2z_D2x_Pz_S_bc+ABY*I_ERI_Fx2z_D2x_Pz_S_bc;
  Double I_ERI_F3y_F2xy_Pz_S_bc = I_ERI_G4y_D2x_Pz_S_bc+ABY*I_ERI_F3y_D2x_Pz_S_bc;
  Double I_ERI_F2yz_F2xy_Pz_S_bc = I_ERI_G3yz_D2x_Pz_S_bc+ABY*I_ERI_F2yz_D2x_Pz_S_bc;
  Double I_ERI_Fy2z_F2xy_Pz_S_bc = I_ERI_G2y2z_D2x_Pz_S_bc+ABY*I_ERI_Fy2z_D2x_Pz_S_bc;
  Double I_ERI_F3z_F2xy_Pz_S_bc = I_ERI_Gy3z_D2x_Pz_S_bc+ABY*I_ERI_F3z_D2x_Pz_S_bc;
  Double I_ERI_F3x_F2xz_Pz_S_bc = I_ERI_G3xz_D2x_Pz_S_bc+ABZ*I_ERI_F3x_D2x_Pz_S_bc;
  Double I_ERI_F2xy_F2xz_Pz_S_bc = I_ERI_G2xyz_D2x_Pz_S_bc+ABZ*I_ERI_F2xy_D2x_Pz_S_bc;
  Double I_ERI_F2xz_F2xz_Pz_S_bc = I_ERI_G2x2z_D2x_Pz_S_bc+ABZ*I_ERI_F2xz_D2x_Pz_S_bc;
  Double I_ERI_Fx2y_F2xz_Pz_S_bc = I_ERI_Gx2yz_D2x_Pz_S_bc+ABZ*I_ERI_Fx2y_D2x_Pz_S_bc;
  Double I_ERI_Fxyz_F2xz_Pz_S_bc = I_ERI_Gxy2z_D2x_Pz_S_bc+ABZ*I_ERI_Fxyz_D2x_Pz_S_bc;
  Double I_ERI_Fx2z_F2xz_Pz_S_bc = I_ERI_Gx3z_D2x_Pz_S_bc+ABZ*I_ERI_Fx2z_D2x_Pz_S_bc;
  Double I_ERI_F3y_F2xz_Pz_S_bc = I_ERI_G3yz_D2x_Pz_S_bc+ABZ*I_ERI_F3y_D2x_Pz_S_bc;
  Double I_ERI_F2yz_F2xz_Pz_S_bc = I_ERI_G2y2z_D2x_Pz_S_bc+ABZ*I_ERI_F2yz_D2x_Pz_S_bc;
  Double I_ERI_Fy2z_F2xz_Pz_S_bc = I_ERI_Gy3z_D2x_Pz_S_bc+ABZ*I_ERI_Fy2z_D2x_Pz_S_bc;
  Double I_ERI_F3z_F2xz_Pz_S_bc = I_ERI_G4z_D2x_Pz_S_bc+ABZ*I_ERI_F3z_D2x_Pz_S_bc;
  Double I_ERI_F3x_Fx2y_Pz_S_bc = I_ERI_G4x_D2y_Pz_S_bc+ABX*I_ERI_F3x_D2y_Pz_S_bc;
  Double I_ERI_F2xy_Fx2y_Pz_S_bc = I_ERI_G3xy_D2y_Pz_S_bc+ABX*I_ERI_F2xy_D2y_Pz_S_bc;
  Double I_ERI_F2xz_Fx2y_Pz_S_bc = I_ERI_G3xz_D2y_Pz_S_bc+ABX*I_ERI_F2xz_D2y_Pz_S_bc;
  Double I_ERI_Fx2y_Fx2y_Pz_S_bc = I_ERI_G2x2y_D2y_Pz_S_bc+ABX*I_ERI_Fx2y_D2y_Pz_S_bc;
  Double I_ERI_Fxyz_Fx2y_Pz_S_bc = I_ERI_G2xyz_D2y_Pz_S_bc+ABX*I_ERI_Fxyz_D2y_Pz_S_bc;
  Double I_ERI_Fx2z_Fx2y_Pz_S_bc = I_ERI_G2x2z_D2y_Pz_S_bc+ABX*I_ERI_Fx2z_D2y_Pz_S_bc;
  Double I_ERI_F3y_Fx2y_Pz_S_bc = I_ERI_Gx3y_D2y_Pz_S_bc+ABX*I_ERI_F3y_D2y_Pz_S_bc;
  Double I_ERI_F2yz_Fx2y_Pz_S_bc = I_ERI_Gx2yz_D2y_Pz_S_bc+ABX*I_ERI_F2yz_D2y_Pz_S_bc;
  Double I_ERI_Fy2z_Fx2y_Pz_S_bc = I_ERI_Gxy2z_D2y_Pz_S_bc+ABX*I_ERI_Fy2z_D2y_Pz_S_bc;
  Double I_ERI_F3z_Fx2y_Pz_S_bc = I_ERI_Gx3z_D2y_Pz_S_bc+ABX*I_ERI_F3z_D2y_Pz_S_bc;
  Double I_ERI_F3x_Fxyz_Pz_S_bc = I_ERI_G3xz_Dxy_Pz_S_bc+ABZ*I_ERI_F3x_Dxy_Pz_S_bc;
  Double I_ERI_F2xy_Fxyz_Pz_S_bc = I_ERI_G2xyz_Dxy_Pz_S_bc+ABZ*I_ERI_F2xy_Dxy_Pz_S_bc;
  Double I_ERI_F2xz_Fxyz_Pz_S_bc = I_ERI_G2x2z_Dxy_Pz_S_bc+ABZ*I_ERI_F2xz_Dxy_Pz_S_bc;
  Double I_ERI_Fx2y_Fxyz_Pz_S_bc = I_ERI_Gx2yz_Dxy_Pz_S_bc+ABZ*I_ERI_Fx2y_Dxy_Pz_S_bc;
  Double I_ERI_Fxyz_Fxyz_Pz_S_bc = I_ERI_Gxy2z_Dxy_Pz_S_bc+ABZ*I_ERI_Fxyz_Dxy_Pz_S_bc;
  Double I_ERI_Fx2z_Fxyz_Pz_S_bc = I_ERI_Gx3z_Dxy_Pz_S_bc+ABZ*I_ERI_Fx2z_Dxy_Pz_S_bc;
  Double I_ERI_F3y_Fxyz_Pz_S_bc = I_ERI_G3yz_Dxy_Pz_S_bc+ABZ*I_ERI_F3y_Dxy_Pz_S_bc;
  Double I_ERI_F2yz_Fxyz_Pz_S_bc = I_ERI_G2y2z_Dxy_Pz_S_bc+ABZ*I_ERI_F2yz_Dxy_Pz_S_bc;
  Double I_ERI_Fy2z_Fxyz_Pz_S_bc = I_ERI_Gy3z_Dxy_Pz_S_bc+ABZ*I_ERI_Fy2z_Dxy_Pz_S_bc;
  Double I_ERI_F3z_Fxyz_Pz_S_bc = I_ERI_G4z_Dxy_Pz_S_bc+ABZ*I_ERI_F3z_Dxy_Pz_S_bc;
  Double I_ERI_F3x_Fx2z_Pz_S_bc = I_ERI_G4x_D2z_Pz_S_bc+ABX*I_ERI_F3x_D2z_Pz_S_bc;
  Double I_ERI_F2xy_Fx2z_Pz_S_bc = I_ERI_G3xy_D2z_Pz_S_bc+ABX*I_ERI_F2xy_D2z_Pz_S_bc;
  Double I_ERI_F2xz_Fx2z_Pz_S_bc = I_ERI_G3xz_D2z_Pz_S_bc+ABX*I_ERI_F2xz_D2z_Pz_S_bc;
  Double I_ERI_Fx2y_Fx2z_Pz_S_bc = I_ERI_G2x2y_D2z_Pz_S_bc+ABX*I_ERI_Fx2y_D2z_Pz_S_bc;
  Double I_ERI_Fxyz_Fx2z_Pz_S_bc = I_ERI_G2xyz_D2z_Pz_S_bc+ABX*I_ERI_Fxyz_D2z_Pz_S_bc;
  Double I_ERI_Fx2z_Fx2z_Pz_S_bc = I_ERI_G2x2z_D2z_Pz_S_bc+ABX*I_ERI_Fx2z_D2z_Pz_S_bc;
  Double I_ERI_F3y_Fx2z_Pz_S_bc = I_ERI_Gx3y_D2z_Pz_S_bc+ABX*I_ERI_F3y_D2z_Pz_S_bc;
  Double I_ERI_F2yz_Fx2z_Pz_S_bc = I_ERI_Gx2yz_D2z_Pz_S_bc+ABX*I_ERI_F2yz_D2z_Pz_S_bc;
  Double I_ERI_Fy2z_Fx2z_Pz_S_bc = I_ERI_Gxy2z_D2z_Pz_S_bc+ABX*I_ERI_Fy2z_D2z_Pz_S_bc;
  Double I_ERI_F3z_Fx2z_Pz_S_bc = I_ERI_Gx3z_D2z_Pz_S_bc+ABX*I_ERI_F3z_D2z_Pz_S_bc;
  Double I_ERI_F3x_F3y_Pz_S_bc = I_ERI_G3xy_D2y_Pz_S_bc+ABY*I_ERI_F3x_D2y_Pz_S_bc;
  Double I_ERI_F2xy_F3y_Pz_S_bc = I_ERI_G2x2y_D2y_Pz_S_bc+ABY*I_ERI_F2xy_D2y_Pz_S_bc;
  Double I_ERI_F2xz_F3y_Pz_S_bc = I_ERI_G2xyz_D2y_Pz_S_bc+ABY*I_ERI_F2xz_D2y_Pz_S_bc;
  Double I_ERI_Fx2y_F3y_Pz_S_bc = I_ERI_Gx3y_D2y_Pz_S_bc+ABY*I_ERI_Fx2y_D2y_Pz_S_bc;
  Double I_ERI_Fxyz_F3y_Pz_S_bc = I_ERI_Gx2yz_D2y_Pz_S_bc+ABY*I_ERI_Fxyz_D2y_Pz_S_bc;
  Double I_ERI_Fx2z_F3y_Pz_S_bc = I_ERI_Gxy2z_D2y_Pz_S_bc+ABY*I_ERI_Fx2z_D2y_Pz_S_bc;
  Double I_ERI_F3y_F3y_Pz_S_bc = I_ERI_G4y_D2y_Pz_S_bc+ABY*I_ERI_F3y_D2y_Pz_S_bc;
  Double I_ERI_F2yz_F3y_Pz_S_bc = I_ERI_G3yz_D2y_Pz_S_bc+ABY*I_ERI_F2yz_D2y_Pz_S_bc;
  Double I_ERI_Fy2z_F3y_Pz_S_bc = I_ERI_G2y2z_D2y_Pz_S_bc+ABY*I_ERI_Fy2z_D2y_Pz_S_bc;
  Double I_ERI_F3z_F3y_Pz_S_bc = I_ERI_Gy3z_D2y_Pz_S_bc+ABY*I_ERI_F3z_D2y_Pz_S_bc;
  Double I_ERI_F3x_F2yz_Pz_S_bc = I_ERI_G3xz_D2y_Pz_S_bc+ABZ*I_ERI_F3x_D2y_Pz_S_bc;
  Double I_ERI_F2xy_F2yz_Pz_S_bc = I_ERI_G2xyz_D2y_Pz_S_bc+ABZ*I_ERI_F2xy_D2y_Pz_S_bc;
  Double I_ERI_F2xz_F2yz_Pz_S_bc = I_ERI_G2x2z_D2y_Pz_S_bc+ABZ*I_ERI_F2xz_D2y_Pz_S_bc;
  Double I_ERI_Fx2y_F2yz_Pz_S_bc = I_ERI_Gx2yz_D2y_Pz_S_bc+ABZ*I_ERI_Fx2y_D2y_Pz_S_bc;
  Double I_ERI_Fxyz_F2yz_Pz_S_bc = I_ERI_Gxy2z_D2y_Pz_S_bc+ABZ*I_ERI_Fxyz_D2y_Pz_S_bc;
  Double I_ERI_Fx2z_F2yz_Pz_S_bc = I_ERI_Gx3z_D2y_Pz_S_bc+ABZ*I_ERI_Fx2z_D2y_Pz_S_bc;
  Double I_ERI_F3y_F2yz_Pz_S_bc = I_ERI_G3yz_D2y_Pz_S_bc+ABZ*I_ERI_F3y_D2y_Pz_S_bc;
  Double I_ERI_F2yz_F2yz_Pz_S_bc = I_ERI_G2y2z_D2y_Pz_S_bc+ABZ*I_ERI_F2yz_D2y_Pz_S_bc;
  Double I_ERI_Fy2z_F2yz_Pz_S_bc = I_ERI_Gy3z_D2y_Pz_S_bc+ABZ*I_ERI_Fy2z_D2y_Pz_S_bc;
  Double I_ERI_F3z_F2yz_Pz_S_bc = I_ERI_G4z_D2y_Pz_S_bc+ABZ*I_ERI_F3z_D2y_Pz_S_bc;
  Double I_ERI_F3x_Fy2z_Pz_S_bc = I_ERI_G3xy_D2z_Pz_S_bc+ABY*I_ERI_F3x_D2z_Pz_S_bc;
  Double I_ERI_F2xy_Fy2z_Pz_S_bc = I_ERI_G2x2y_D2z_Pz_S_bc+ABY*I_ERI_F2xy_D2z_Pz_S_bc;
  Double I_ERI_F2xz_Fy2z_Pz_S_bc = I_ERI_G2xyz_D2z_Pz_S_bc+ABY*I_ERI_F2xz_D2z_Pz_S_bc;
  Double I_ERI_Fx2y_Fy2z_Pz_S_bc = I_ERI_Gx3y_D2z_Pz_S_bc+ABY*I_ERI_Fx2y_D2z_Pz_S_bc;
  Double I_ERI_Fxyz_Fy2z_Pz_S_bc = I_ERI_Gx2yz_D2z_Pz_S_bc+ABY*I_ERI_Fxyz_D2z_Pz_S_bc;
  Double I_ERI_Fx2z_Fy2z_Pz_S_bc = I_ERI_Gxy2z_D2z_Pz_S_bc+ABY*I_ERI_Fx2z_D2z_Pz_S_bc;
  Double I_ERI_F3y_Fy2z_Pz_S_bc = I_ERI_G4y_D2z_Pz_S_bc+ABY*I_ERI_F3y_D2z_Pz_S_bc;
  Double I_ERI_F2yz_Fy2z_Pz_S_bc = I_ERI_G3yz_D2z_Pz_S_bc+ABY*I_ERI_F2yz_D2z_Pz_S_bc;
  Double I_ERI_Fy2z_Fy2z_Pz_S_bc = I_ERI_G2y2z_D2z_Pz_S_bc+ABY*I_ERI_Fy2z_D2z_Pz_S_bc;
  Double I_ERI_F3z_Fy2z_Pz_S_bc = I_ERI_Gy3z_D2z_Pz_S_bc+ABY*I_ERI_F3z_D2z_Pz_S_bc;
  Double I_ERI_F3x_F3z_Pz_S_bc = I_ERI_G3xz_D2z_Pz_S_bc+ABZ*I_ERI_F3x_D2z_Pz_S_bc;
  Double I_ERI_F2xy_F3z_Pz_S_bc = I_ERI_G2xyz_D2z_Pz_S_bc+ABZ*I_ERI_F2xy_D2z_Pz_S_bc;
  Double I_ERI_F2xz_F3z_Pz_S_bc = I_ERI_G2x2z_D2z_Pz_S_bc+ABZ*I_ERI_F2xz_D2z_Pz_S_bc;
  Double I_ERI_Fx2y_F3z_Pz_S_bc = I_ERI_Gx2yz_D2z_Pz_S_bc+ABZ*I_ERI_Fx2y_D2z_Pz_S_bc;
  Double I_ERI_Fxyz_F3z_Pz_S_bc = I_ERI_Gxy2z_D2z_Pz_S_bc+ABZ*I_ERI_Fxyz_D2z_Pz_S_bc;
  Double I_ERI_Fx2z_F3z_Pz_S_bc = I_ERI_Gx3z_D2z_Pz_S_bc+ABZ*I_ERI_Fx2z_D2z_Pz_S_bc;
  Double I_ERI_F3y_F3z_Pz_S_bc = I_ERI_G3yz_D2z_Pz_S_bc+ABZ*I_ERI_F3y_D2z_Pz_S_bc;
  Double I_ERI_F2yz_F3z_Pz_S_bc = I_ERI_G2y2z_D2z_Pz_S_bc+ABZ*I_ERI_F2yz_D2z_Pz_S_bc;
  Double I_ERI_Fy2z_F3z_Pz_S_bc = I_ERI_Gy3z_D2z_Pz_S_bc+ABZ*I_ERI_Fy2z_D2z_Pz_S_bc;
  Double I_ERI_F3z_F3z_Pz_S_bc = I_ERI_G4z_D2z_Pz_S_bc+ABZ*I_ERI_F3z_D2z_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_S_cc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_cc
   * RHS shell quartet name: SQ_ERI_F_S_D_S_cc
   ************************************************************/
  Double I_ERI_F3x_Px_D2x_S_cc = I_ERI_G4x_S_D2x_S_cc+ABX*I_ERI_F3x_S_D2x_S_cc;
  Double I_ERI_F2xy_Px_D2x_S_cc = I_ERI_G3xy_S_D2x_S_cc+ABX*I_ERI_F2xy_S_D2x_S_cc;
  Double I_ERI_F2xz_Px_D2x_S_cc = I_ERI_G3xz_S_D2x_S_cc+ABX*I_ERI_F2xz_S_D2x_S_cc;
  Double I_ERI_Fx2y_Px_D2x_S_cc = I_ERI_G2x2y_S_D2x_S_cc+ABX*I_ERI_Fx2y_S_D2x_S_cc;
  Double I_ERI_Fxyz_Px_D2x_S_cc = I_ERI_G2xyz_S_D2x_S_cc+ABX*I_ERI_Fxyz_S_D2x_S_cc;
  Double I_ERI_Fx2z_Px_D2x_S_cc = I_ERI_G2x2z_S_D2x_S_cc+ABX*I_ERI_Fx2z_S_D2x_S_cc;
  Double I_ERI_F3y_Px_D2x_S_cc = I_ERI_Gx3y_S_D2x_S_cc+ABX*I_ERI_F3y_S_D2x_S_cc;
  Double I_ERI_F2yz_Px_D2x_S_cc = I_ERI_Gx2yz_S_D2x_S_cc+ABX*I_ERI_F2yz_S_D2x_S_cc;
  Double I_ERI_Fy2z_Px_D2x_S_cc = I_ERI_Gxy2z_S_D2x_S_cc+ABX*I_ERI_Fy2z_S_D2x_S_cc;
  Double I_ERI_F3z_Px_D2x_S_cc = I_ERI_Gx3z_S_D2x_S_cc+ABX*I_ERI_F3z_S_D2x_S_cc;
  Double I_ERI_F3x_Py_D2x_S_cc = I_ERI_G3xy_S_D2x_S_cc+ABY*I_ERI_F3x_S_D2x_S_cc;
  Double I_ERI_F2xy_Py_D2x_S_cc = I_ERI_G2x2y_S_D2x_S_cc+ABY*I_ERI_F2xy_S_D2x_S_cc;
  Double I_ERI_F2xz_Py_D2x_S_cc = I_ERI_G2xyz_S_D2x_S_cc+ABY*I_ERI_F2xz_S_D2x_S_cc;
  Double I_ERI_Fx2y_Py_D2x_S_cc = I_ERI_Gx3y_S_D2x_S_cc+ABY*I_ERI_Fx2y_S_D2x_S_cc;
  Double I_ERI_Fxyz_Py_D2x_S_cc = I_ERI_Gx2yz_S_D2x_S_cc+ABY*I_ERI_Fxyz_S_D2x_S_cc;
  Double I_ERI_Fx2z_Py_D2x_S_cc = I_ERI_Gxy2z_S_D2x_S_cc+ABY*I_ERI_Fx2z_S_D2x_S_cc;
  Double I_ERI_F3y_Py_D2x_S_cc = I_ERI_G4y_S_D2x_S_cc+ABY*I_ERI_F3y_S_D2x_S_cc;
  Double I_ERI_F2yz_Py_D2x_S_cc = I_ERI_G3yz_S_D2x_S_cc+ABY*I_ERI_F2yz_S_D2x_S_cc;
  Double I_ERI_Fy2z_Py_D2x_S_cc = I_ERI_G2y2z_S_D2x_S_cc+ABY*I_ERI_Fy2z_S_D2x_S_cc;
  Double I_ERI_F3z_Py_D2x_S_cc = I_ERI_Gy3z_S_D2x_S_cc+ABY*I_ERI_F3z_S_D2x_S_cc;
  Double I_ERI_F3x_Pz_D2x_S_cc = I_ERI_G3xz_S_D2x_S_cc+ABZ*I_ERI_F3x_S_D2x_S_cc;
  Double I_ERI_F2xy_Pz_D2x_S_cc = I_ERI_G2xyz_S_D2x_S_cc+ABZ*I_ERI_F2xy_S_D2x_S_cc;
  Double I_ERI_F2xz_Pz_D2x_S_cc = I_ERI_G2x2z_S_D2x_S_cc+ABZ*I_ERI_F2xz_S_D2x_S_cc;
  Double I_ERI_Fx2y_Pz_D2x_S_cc = I_ERI_Gx2yz_S_D2x_S_cc+ABZ*I_ERI_Fx2y_S_D2x_S_cc;
  Double I_ERI_Fxyz_Pz_D2x_S_cc = I_ERI_Gxy2z_S_D2x_S_cc+ABZ*I_ERI_Fxyz_S_D2x_S_cc;
  Double I_ERI_Fx2z_Pz_D2x_S_cc = I_ERI_Gx3z_S_D2x_S_cc+ABZ*I_ERI_Fx2z_S_D2x_S_cc;
  Double I_ERI_F3y_Pz_D2x_S_cc = I_ERI_G3yz_S_D2x_S_cc+ABZ*I_ERI_F3y_S_D2x_S_cc;
  Double I_ERI_F2yz_Pz_D2x_S_cc = I_ERI_G2y2z_S_D2x_S_cc+ABZ*I_ERI_F2yz_S_D2x_S_cc;
  Double I_ERI_Fy2z_Pz_D2x_S_cc = I_ERI_Gy3z_S_D2x_S_cc+ABZ*I_ERI_Fy2z_S_D2x_S_cc;
  Double I_ERI_F3z_Pz_D2x_S_cc = I_ERI_G4z_S_D2x_S_cc+ABZ*I_ERI_F3z_S_D2x_S_cc;
  Double I_ERI_F3x_Px_Dxy_S_cc = I_ERI_G4x_S_Dxy_S_cc+ABX*I_ERI_F3x_S_Dxy_S_cc;
  Double I_ERI_F2xy_Px_Dxy_S_cc = I_ERI_G3xy_S_Dxy_S_cc+ABX*I_ERI_F2xy_S_Dxy_S_cc;
  Double I_ERI_F2xz_Px_Dxy_S_cc = I_ERI_G3xz_S_Dxy_S_cc+ABX*I_ERI_F2xz_S_Dxy_S_cc;
  Double I_ERI_Fx2y_Px_Dxy_S_cc = I_ERI_G2x2y_S_Dxy_S_cc+ABX*I_ERI_Fx2y_S_Dxy_S_cc;
  Double I_ERI_Fxyz_Px_Dxy_S_cc = I_ERI_G2xyz_S_Dxy_S_cc+ABX*I_ERI_Fxyz_S_Dxy_S_cc;
  Double I_ERI_Fx2z_Px_Dxy_S_cc = I_ERI_G2x2z_S_Dxy_S_cc+ABX*I_ERI_Fx2z_S_Dxy_S_cc;
  Double I_ERI_F3y_Px_Dxy_S_cc = I_ERI_Gx3y_S_Dxy_S_cc+ABX*I_ERI_F3y_S_Dxy_S_cc;
  Double I_ERI_F2yz_Px_Dxy_S_cc = I_ERI_Gx2yz_S_Dxy_S_cc+ABX*I_ERI_F2yz_S_Dxy_S_cc;
  Double I_ERI_Fy2z_Px_Dxy_S_cc = I_ERI_Gxy2z_S_Dxy_S_cc+ABX*I_ERI_Fy2z_S_Dxy_S_cc;
  Double I_ERI_F3z_Px_Dxy_S_cc = I_ERI_Gx3z_S_Dxy_S_cc+ABX*I_ERI_F3z_S_Dxy_S_cc;
  Double I_ERI_F3x_Py_Dxy_S_cc = I_ERI_G3xy_S_Dxy_S_cc+ABY*I_ERI_F3x_S_Dxy_S_cc;
  Double I_ERI_F2xy_Py_Dxy_S_cc = I_ERI_G2x2y_S_Dxy_S_cc+ABY*I_ERI_F2xy_S_Dxy_S_cc;
  Double I_ERI_F2xz_Py_Dxy_S_cc = I_ERI_G2xyz_S_Dxy_S_cc+ABY*I_ERI_F2xz_S_Dxy_S_cc;
  Double I_ERI_Fx2y_Py_Dxy_S_cc = I_ERI_Gx3y_S_Dxy_S_cc+ABY*I_ERI_Fx2y_S_Dxy_S_cc;
  Double I_ERI_Fxyz_Py_Dxy_S_cc = I_ERI_Gx2yz_S_Dxy_S_cc+ABY*I_ERI_Fxyz_S_Dxy_S_cc;
  Double I_ERI_Fx2z_Py_Dxy_S_cc = I_ERI_Gxy2z_S_Dxy_S_cc+ABY*I_ERI_Fx2z_S_Dxy_S_cc;
  Double I_ERI_F3y_Py_Dxy_S_cc = I_ERI_G4y_S_Dxy_S_cc+ABY*I_ERI_F3y_S_Dxy_S_cc;
  Double I_ERI_F2yz_Py_Dxy_S_cc = I_ERI_G3yz_S_Dxy_S_cc+ABY*I_ERI_F2yz_S_Dxy_S_cc;
  Double I_ERI_Fy2z_Py_Dxy_S_cc = I_ERI_G2y2z_S_Dxy_S_cc+ABY*I_ERI_Fy2z_S_Dxy_S_cc;
  Double I_ERI_F3z_Py_Dxy_S_cc = I_ERI_Gy3z_S_Dxy_S_cc+ABY*I_ERI_F3z_S_Dxy_S_cc;
  Double I_ERI_F3x_Pz_Dxy_S_cc = I_ERI_G3xz_S_Dxy_S_cc+ABZ*I_ERI_F3x_S_Dxy_S_cc;
  Double I_ERI_F2xy_Pz_Dxy_S_cc = I_ERI_G2xyz_S_Dxy_S_cc+ABZ*I_ERI_F2xy_S_Dxy_S_cc;
  Double I_ERI_F2xz_Pz_Dxy_S_cc = I_ERI_G2x2z_S_Dxy_S_cc+ABZ*I_ERI_F2xz_S_Dxy_S_cc;
  Double I_ERI_Fx2y_Pz_Dxy_S_cc = I_ERI_Gx2yz_S_Dxy_S_cc+ABZ*I_ERI_Fx2y_S_Dxy_S_cc;
  Double I_ERI_Fxyz_Pz_Dxy_S_cc = I_ERI_Gxy2z_S_Dxy_S_cc+ABZ*I_ERI_Fxyz_S_Dxy_S_cc;
  Double I_ERI_Fx2z_Pz_Dxy_S_cc = I_ERI_Gx3z_S_Dxy_S_cc+ABZ*I_ERI_Fx2z_S_Dxy_S_cc;
  Double I_ERI_F3y_Pz_Dxy_S_cc = I_ERI_G3yz_S_Dxy_S_cc+ABZ*I_ERI_F3y_S_Dxy_S_cc;
  Double I_ERI_F2yz_Pz_Dxy_S_cc = I_ERI_G2y2z_S_Dxy_S_cc+ABZ*I_ERI_F2yz_S_Dxy_S_cc;
  Double I_ERI_Fy2z_Pz_Dxy_S_cc = I_ERI_Gy3z_S_Dxy_S_cc+ABZ*I_ERI_Fy2z_S_Dxy_S_cc;
  Double I_ERI_F3z_Pz_Dxy_S_cc = I_ERI_G4z_S_Dxy_S_cc+ABZ*I_ERI_F3z_S_Dxy_S_cc;
  Double I_ERI_F3x_Px_Dxz_S_cc = I_ERI_G4x_S_Dxz_S_cc+ABX*I_ERI_F3x_S_Dxz_S_cc;
  Double I_ERI_F2xy_Px_Dxz_S_cc = I_ERI_G3xy_S_Dxz_S_cc+ABX*I_ERI_F2xy_S_Dxz_S_cc;
  Double I_ERI_F2xz_Px_Dxz_S_cc = I_ERI_G3xz_S_Dxz_S_cc+ABX*I_ERI_F2xz_S_Dxz_S_cc;
  Double I_ERI_Fx2y_Px_Dxz_S_cc = I_ERI_G2x2y_S_Dxz_S_cc+ABX*I_ERI_Fx2y_S_Dxz_S_cc;
  Double I_ERI_Fxyz_Px_Dxz_S_cc = I_ERI_G2xyz_S_Dxz_S_cc+ABX*I_ERI_Fxyz_S_Dxz_S_cc;
  Double I_ERI_Fx2z_Px_Dxz_S_cc = I_ERI_G2x2z_S_Dxz_S_cc+ABX*I_ERI_Fx2z_S_Dxz_S_cc;
  Double I_ERI_F3y_Px_Dxz_S_cc = I_ERI_Gx3y_S_Dxz_S_cc+ABX*I_ERI_F3y_S_Dxz_S_cc;
  Double I_ERI_F2yz_Px_Dxz_S_cc = I_ERI_Gx2yz_S_Dxz_S_cc+ABX*I_ERI_F2yz_S_Dxz_S_cc;
  Double I_ERI_Fy2z_Px_Dxz_S_cc = I_ERI_Gxy2z_S_Dxz_S_cc+ABX*I_ERI_Fy2z_S_Dxz_S_cc;
  Double I_ERI_F3z_Px_Dxz_S_cc = I_ERI_Gx3z_S_Dxz_S_cc+ABX*I_ERI_F3z_S_Dxz_S_cc;
  Double I_ERI_F3x_Py_Dxz_S_cc = I_ERI_G3xy_S_Dxz_S_cc+ABY*I_ERI_F3x_S_Dxz_S_cc;
  Double I_ERI_F2xy_Py_Dxz_S_cc = I_ERI_G2x2y_S_Dxz_S_cc+ABY*I_ERI_F2xy_S_Dxz_S_cc;
  Double I_ERI_F2xz_Py_Dxz_S_cc = I_ERI_G2xyz_S_Dxz_S_cc+ABY*I_ERI_F2xz_S_Dxz_S_cc;
  Double I_ERI_Fx2y_Py_Dxz_S_cc = I_ERI_Gx3y_S_Dxz_S_cc+ABY*I_ERI_Fx2y_S_Dxz_S_cc;
  Double I_ERI_Fxyz_Py_Dxz_S_cc = I_ERI_Gx2yz_S_Dxz_S_cc+ABY*I_ERI_Fxyz_S_Dxz_S_cc;
  Double I_ERI_Fx2z_Py_Dxz_S_cc = I_ERI_Gxy2z_S_Dxz_S_cc+ABY*I_ERI_Fx2z_S_Dxz_S_cc;
  Double I_ERI_F3y_Py_Dxz_S_cc = I_ERI_G4y_S_Dxz_S_cc+ABY*I_ERI_F3y_S_Dxz_S_cc;
  Double I_ERI_F2yz_Py_Dxz_S_cc = I_ERI_G3yz_S_Dxz_S_cc+ABY*I_ERI_F2yz_S_Dxz_S_cc;
  Double I_ERI_Fy2z_Py_Dxz_S_cc = I_ERI_G2y2z_S_Dxz_S_cc+ABY*I_ERI_Fy2z_S_Dxz_S_cc;
  Double I_ERI_F3z_Py_Dxz_S_cc = I_ERI_Gy3z_S_Dxz_S_cc+ABY*I_ERI_F3z_S_Dxz_S_cc;
  Double I_ERI_F3x_Pz_Dxz_S_cc = I_ERI_G3xz_S_Dxz_S_cc+ABZ*I_ERI_F3x_S_Dxz_S_cc;
  Double I_ERI_F2xy_Pz_Dxz_S_cc = I_ERI_G2xyz_S_Dxz_S_cc+ABZ*I_ERI_F2xy_S_Dxz_S_cc;
  Double I_ERI_F2xz_Pz_Dxz_S_cc = I_ERI_G2x2z_S_Dxz_S_cc+ABZ*I_ERI_F2xz_S_Dxz_S_cc;
  Double I_ERI_Fx2y_Pz_Dxz_S_cc = I_ERI_Gx2yz_S_Dxz_S_cc+ABZ*I_ERI_Fx2y_S_Dxz_S_cc;
  Double I_ERI_Fxyz_Pz_Dxz_S_cc = I_ERI_Gxy2z_S_Dxz_S_cc+ABZ*I_ERI_Fxyz_S_Dxz_S_cc;
  Double I_ERI_Fx2z_Pz_Dxz_S_cc = I_ERI_Gx3z_S_Dxz_S_cc+ABZ*I_ERI_Fx2z_S_Dxz_S_cc;
  Double I_ERI_F3y_Pz_Dxz_S_cc = I_ERI_G3yz_S_Dxz_S_cc+ABZ*I_ERI_F3y_S_Dxz_S_cc;
  Double I_ERI_F2yz_Pz_Dxz_S_cc = I_ERI_G2y2z_S_Dxz_S_cc+ABZ*I_ERI_F2yz_S_Dxz_S_cc;
  Double I_ERI_Fy2z_Pz_Dxz_S_cc = I_ERI_Gy3z_S_Dxz_S_cc+ABZ*I_ERI_Fy2z_S_Dxz_S_cc;
  Double I_ERI_F3z_Pz_Dxz_S_cc = I_ERI_G4z_S_Dxz_S_cc+ABZ*I_ERI_F3z_S_Dxz_S_cc;
  Double I_ERI_F3x_Px_D2y_S_cc = I_ERI_G4x_S_D2y_S_cc+ABX*I_ERI_F3x_S_D2y_S_cc;
  Double I_ERI_F2xy_Px_D2y_S_cc = I_ERI_G3xy_S_D2y_S_cc+ABX*I_ERI_F2xy_S_D2y_S_cc;
  Double I_ERI_F2xz_Px_D2y_S_cc = I_ERI_G3xz_S_D2y_S_cc+ABX*I_ERI_F2xz_S_D2y_S_cc;
  Double I_ERI_Fx2y_Px_D2y_S_cc = I_ERI_G2x2y_S_D2y_S_cc+ABX*I_ERI_Fx2y_S_D2y_S_cc;
  Double I_ERI_Fxyz_Px_D2y_S_cc = I_ERI_G2xyz_S_D2y_S_cc+ABX*I_ERI_Fxyz_S_D2y_S_cc;
  Double I_ERI_Fx2z_Px_D2y_S_cc = I_ERI_G2x2z_S_D2y_S_cc+ABX*I_ERI_Fx2z_S_D2y_S_cc;
  Double I_ERI_F3y_Px_D2y_S_cc = I_ERI_Gx3y_S_D2y_S_cc+ABX*I_ERI_F3y_S_D2y_S_cc;
  Double I_ERI_F2yz_Px_D2y_S_cc = I_ERI_Gx2yz_S_D2y_S_cc+ABX*I_ERI_F2yz_S_D2y_S_cc;
  Double I_ERI_Fy2z_Px_D2y_S_cc = I_ERI_Gxy2z_S_D2y_S_cc+ABX*I_ERI_Fy2z_S_D2y_S_cc;
  Double I_ERI_F3z_Px_D2y_S_cc = I_ERI_Gx3z_S_D2y_S_cc+ABX*I_ERI_F3z_S_D2y_S_cc;
  Double I_ERI_F3x_Py_D2y_S_cc = I_ERI_G3xy_S_D2y_S_cc+ABY*I_ERI_F3x_S_D2y_S_cc;
  Double I_ERI_F2xy_Py_D2y_S_cc = I_ERI_G2x2y_S_D2y_S_cc+ABY*I_ERI_F2xy_S_D2y_S_cc;
  Double I_ERI_F2xz_Py_D2y_S_cc = I_ERI_G2xyz_S_D2y_S_cc+ABY*I_ERI_F2xz_S_D2y_S_cc;
  Double I_ERI_Fx2y_Py_D2y_S_cc = I_ERI_Gx3y_S_D2y_S_cc+ABY*I_ERI_Fx2y_S_D2y_S_cc;
  Double I_ERI_Fxyz_Py_D2y_S_cc = I_ERI_Gx2yz_S_D2y_S_cc+ABY*I_ERI_Fxyz_S_D2y_S_cc;
  Double I_ERI_Fx2z_Py_D2y_S_cc = I_ERI_Gxy2z_S_D2y_S_cc+ABY*I_ERI_Fx2z_S_D2y_S_cc;
  Double I_ERI_F3y_Py_D2y_S_cc = I_ERI_G4y_S_D2y_S_cc+ABY*I_ERI_F3y_S_D2y_S_cc;
  Double I_ERI_F2yz_Py_D2y_S_cc = I_ERI_G3yz_S_D2y_S_cc+ABY*I_ERI_F2yz_S_D2y_S_cc;
  Double I_ERI_Fy2z_Py_D2y_S_cc = I_ERI_G2y2z_S_D2y_S_cc+ABY*I_ERI_Fy2z_S_D2y_S_cc;
  Double I_ERI_F3z_Py_D2y_S_cc = I_ERI_Gy3z_S_D2y_S_cc+ABY*I_ERI_F3z_S_D2y_S_cc;
  Double I_ERI_F3x_Pz_D2y_S_cc = I_ERI_G3xz_S_D2y_S_cc+ABZ*I_ERI_F3x_S_D2y_S_cc;
  Double I_ERI_F2xy_Pz_D2y_S_cc = I_ERI_G2xyz_S_D2y_S_cc+ABZ*I_ERI_F2xy_S_D2y_S_cc;
  Double I_ERI_F2xz_Pz_D2y_S_cc = I_ERI_G2x2z_S_D2y_S_cc+ABZ*I_ERI_F2xz_S_D2y_S_cc;
  Double I_ERI_Fx2y_Pz_D2y_S_cc = I_ERI_Gx2yz_S_D2y_S_cc+ABZ*I_ERI_Fx2y_S_D2y_S_cc;
  Double I_ERI_Fxyz_Pz_D2y_S_cc = I_ERI_Gxy2z_S_D2y_S_cc+ABZ*I_ERI_Fxyz_S_D2y_S_cc;
  Double I_ERI_Fx2z_Pz_D2y_S_cc = I_ERI_Gx3z_S_D2y_S_cc+ABZ*I_ERI_Fx2z_S_D2y_S_cc;
  Double I_ERI_F3y_Pz_D2y_S_cc = I_ERI_G3yz_S_D2y_S_cc+ABZ*I_ERI_F3y_S_D2y_S_cc;
  Double I_ERI_F2yz_Pz_D2y_S_cc = I_ERI_G2y2z_S_D2y_S_cc+ABZ*I_ERI_F2yz_S_D2y_S_cc;
  Double I_ERI_Fy2z_Pz_D2y_S_cc = I_ERI_Gy3z_S_D2y_S_cc+ABZ*I_ERI_Fy2z_S_D2y_S_cc;
  Double I_ERI_F3z_Pz_D2y_S_cc = I_ERI_G4z_S_D2y_S_cc+ABZ*I_ERI_F3z_S_D2y_S_cc;
  Double I_ERI_F3x_Px_Dyz_S_cc = I_ERI_G4x_S_Dyz_S_cc+ABX*I_ERI_F3x_S_Dyz_S_cc;
  Double I_ERI_F2xy_Px_Dyz_S_cc = I_ERI_G3xy_S_Dyz_S_cc+ABX*I_ERI_F2xy_S_Dyz_S_cc;
  Double I_ERI_F2xz_Px_Dyz_S_cc = I_ERI_G3xz_S_Dyz_S_cc+ABX*I_ERI_F2xz_S_Dyz_S_cc;
  Double I_ERI_Fx2y_Px_Dyz_S_cc = I_ERI_G2x2y_S_Dyz_S_cc+ABX*I_ERI_Fx2y_S_Dyz_S_cc;
  Double I_ERI_Fxyz_Px_Dyz_S_cc = I_ERI_G2xyz_S_Dyz_S_cc+ABX*I_ERI_Fxyz_S_Dyz_S_cc;
  Double I_ERI_Fx2z_Px_Dyz_S_cc = I_ERI_G2x2z_S_Dyz_S_cc+ABX*I_ERI_Fx2z_S_Dyz_S_cc;
  Double I_ERI_F3y_Px_Dyz_S_cc = I_ERI_Gx3y_S_Dyz_S_cc+ABX*I_ERI_F3y_S_Dyz_S_cc;
  Double I_ERI_F2yz_Px_Dyz_S_cc = I_ERI_Gx2yz_S_Dyz_S_cc+ABX*I_ERI_F2yz_S_Dyz_S_cc;
  Double I_ERI_Fy2z_Px_Dyz_S_cc = I_ERI_Gxy2z_S_Dyz_S_cc+ABX*I_ERI_Fy2z_S_Dyz_S_cc;
  Double I_ERI_F3z_Px_Dyz_S_cc = I_ERI_Gx3z_S_Dyz_S_cc+ABX*I_ERI_F3z_S_Dyz_S_cc;
  Double I_ERI_F3x_Py_Dyz_S_cc = I_ERI_G3xy_S_Dyz_S_cc+ABY*I_ERI_F3x_S_Dyz_S_cc;
  Double I_ERI_F2xy_Py_Dyz_S_cc = I_ERI_G2x2y_S_Dyz_S_cc+ABY*I_ERI_F2xy_S_Dyz_S_cc;
  Double I_ERI_F2xz_Py_Dyz_S_cc = I_ERI_G2xyz_S_Dyz_S_cc+ABY*I_ERI_F2xz_S_Dyz_S_cc;
  Double I_ERI_Fx2y_Py_Dyz_S_cc = I_ERI_Gx3y_S_Dyz_S_cc+ABY*I_ERI_Fx2y_S_Dyz_S_cc;
  Double I_ERI_Fxyz_Py_Dyz_S_cc = I_ERI_Gx2yz_S_Dyz_S_cc+ABY*I_ERI_Fxyz_S_Dyz_S_cc;
  Double I_ERI_Fx2z_Py_Dyz_S_cc = I_ERI_Gxy2z_S_Dyz_S_cc+ABY*I_ERI_Fx2z_S_Dyz_S_cc;
  Double I_ERI_F3y_Py_Dyz_S_cc = I_ERI_G4y_S_Dyz_S_cc+ABY*I_ERI_F3y_S_Dyz_S_cc;
  Double I_ERI_F2yz_Py_Dyz_S_cc = I_ERI_G3yz_S_Dyz_S_cc+ABY*I_ERI_F2yz_S_Dyz_S_cc;
  Double I_ERI_Fy2z_Py_Dyz_S_cc = I_ERI_G2y2z_S_Dyz_S_cc+ABY*I_ERI_Fy2z_S_Dyz_S_cc;
  Double I_ERI_F3z_Py_Dyz_S_cc = I_ERI_Gy3z_S_Dyz_S_cc+ABY*I_ERI_F3z_S_Dyz_S_cc;
  Double I_ERI_F3x_Pz_Dyz_S_cc = I_ERI_G3xz_S_Dyz_S_cc+ABZ*I_ERI_F3x_S_Dyz_S_cc;
  Double I_ERI_F2xy_Pz_Dyz_S_cc = I_ERI_G2xyz_S_Dyz_S_cc+ABZ*I_ERI_F2xy_S_Dyz_S_cc;
  Double I_ERI_F2xz_Pz_Dyz_S_cc = I_ERI_G2x2z_S_Dyz_S_cc+ABZ*I_ERI_F2xz_S_Dyz_S_cc;
  Double I_ERI_Fx2y_Pz_Dyz_S_cc = I_ERI_Gx2yz_S_Dyz_S_cc+ABZ*I_ERI_Fx2y_S_Dyz_S_cc;
  Double I_ERI_Fxyz_Pz_Dyz_S_cc = I_ERI_Gxy2z_S_Dyz_S_cc+ABZ*I_ERI_Fxyz_S_Dyz_S_cc;
  Double I_ERI_Fx2z_Pz_Dyz_S_cc = I_ERI_Gx3z_S_Dyz_S_cc+ABZ*I_ERI_Fx2z_S_Dyz_S_cc;
  Double I_ERI_F3y_Pz_Dyz_S_cc = I_ERI_G3yz_S_Dyz_S_cc+ABZ*I_ERI_F3y_S_Dyz_S_cc;
  Double I_ERI_F2yz_Pz_Dyz_S_cc = I_ERI_G2y2z_S_Dyz_S_cc+ABZ*I_ERI_F2yz_S_Dyz_S_cc;
  Double I_ERI_Fy2z_Pz_Dyz_S_cc = I_ERI_Gy3z_S_Dyz_S_cc+ABZ*I_ERI_Fy2z_S_Dyz_S_cc;
  Double I_ERI_F3z_Pz_Dyz_S_cc = I_ERI_G4z_S_Dyz_S_cc+ABZ*I_ERI_F3z_S_Dyz_S_cc;
  Double I_ERI_F3x_Px_D2z_S_cc = I_ERI_G4x_S_D2z_S_cc+ABX*I_ERI_F3x_S_D2z_S_cc;
  Double I_ERI_F2xy_Px_D2z_S_cc = I_ERI_G3xy_S_D2z_S_cc+ABX*I_ERI_F2xy_S_D2z_S_cc;
  Double I_ERI_F2xz_Px_D2z_S_cc = I_ERI_G3xz_S_D2z_S_cc+ABX*I_ERI_F2xz_S_D2z_S_cc;
  Double I_ERI_Fx2y_Px_D2z_S_cc = I_ERI_G2x2y_S_D2z_S_cc+ABX*I_ERI_Fx2y_S_D2z_S_cc;
  Double I_ERI_Fxyz_Px_D2z_S_cc = I_ERI_G2xyz_S_D2z_S_cc+ABX*I_ERI_Fxyz_S_D2z_S_cc;
  Double I_ERI_Fx2z_Px_D2z_S_cc = I_ERI_G2x2z_S_D2z_S_cc+ABX*I_ERI_Fx2z_S_D2z_S_cc;
  Double I_ERI_F3y_Px_D2z_S_cc = I_ERI_Gx3y_S_D2z_S_cc+ABX*I_ERI_F3y_S_D2z_S_cc;
  Double I_ERI_F2yz_Px_D2z_S_cc = I_ERI_Gx2yz_S_D2z_S_cc+ABX*I_ERI_F2yz_S_D2z_S_cc;
  Double I_ERI_Fy2z_Px_D2z_S_cc = I_ERI_Gxy2z_S_D2z_S_cc+ABX*I_ERI_Fy2z_S_D2z_S_cc;
  Double I_ERI_F3z_Px_D2z_S_cc = I_ERI_Gx3z_S_D2z_S_cc+ABX*I_ERI_F3z_S_D2z_S_cc;
  Double I_ERI_F3x_Py_D2z_S_cc = I_ERI_G3xy_S_D2z_S_cc+ABY*I_ERI_F3x_S_D2z_S_cc;
  Double I_ERI_F2xy_Py_D2z_S_cc = I_ERI_G2x2y_S_D2z_S_cc+ABY*I_ERI_F2xy_S_D2z_S_cc;
  Double I_ERI_F2xz_Py_D2z_S_cc = I_ERI_G2xyz_S_D2z_S_cc+ABY*I_ERI_F2xz_S_D2z_S_cc;
  Double I_ERI_Fx2y_Py_D2z_S_cc = I_ERI_Gx3y_S_D2z_S_cc+ABY*I_ERI_Fx2y_S_D2z_S_cc;
  Double I_ERI_Fxyz_Py_D2z_S_cc = I_ERI_Gx2yz_S_D2z_S_cc+ABY*I_ERI_Fxyz_S_D2z_S_cc;
  Double I_ERI_Fx2z_Py_D2z_S_cc = I_ERI_Gxy2z_S_D2z_S_cc+ABY*I_ERI_Fx2z_S_D2z_S_cc;
  Double I_ERI_F3y_Py_D2z_S_cc = I_ERI_G4y_S_D2z_S_cc+ABY*I_ERI_F3y_S_D2z_S_cc;
  Double I_ERI_F2yz_Py_D2z_S_cc = I_ERI_G3yz_S_D2z_S_cc+ABY*I_ERI_F2yz_S_D2z_S_cc;
  Double I_ERI_Fy2z_Py_D2z_S_cc = I_ERI_G2y2z_S_D2z_S_cc+ABY*I_ERI_Fy2z_S_D2z_S_cc;
  Double I_ERI_F3z_Py_D2z_S_cc = I_ERI_Gy3z_S_D2z_S_cc+ABY*I_ERI_F3z_S_D2z_S_cc;
  Double I_ERI_F3x_Pz_D2z_S_cc = I_ERI_G3xz_S_D2z_S_cc+ABZ*I_ERI_F3x_S_D2z_S_cc;
  Double I_ERI_F2xy_Pz_D2z_S_cc = I_ERI_G2xyz_S_D2z_S_cc+ABZ*I_ERI_F2xy_S_D2z_S_cc;
  Double I_ERI_F2xz_Pz_D2z_S_cc = I_ERI_G2x2z_S_D2z_S_cc+ABZ*I_ERI_F2xz_S_D2z_S_cc;
  Double I_ERI_Fx2y_Pz_D2z_S_cc = I_ERI_Gx2yz_S_D2z_S_cc+ABZ*I_ERI_Fx2y_S_D2z_S_cc;
  Double I_ERI_Fxyz_Pz_D2z_S_cc = I_ERI_Gxy2z_S_D2z_S_cc+ABZ*I_ERI_Fxyz_S_D2z_S_cc;
  Double I_ERI_Fx2z_Pz_D2z_S_cc = I_ERI_Gx3z_S_D2z_S_cc+ABZ*I_ERI_Fx2z_S_D2z_S_cc;
  Double I_ERI_F3y_Pz_D2z_S_cc = I_ERI_G3yz_S_D2z_S_cc+ABZ*I_ERI_F3y_S_D2z_S_cc;
  Double I_ERI_F2yz_Pz_D2z_S_cc = I_ERI_G2y2z_S_D2z_S_cc+ABZ*I_ERI_F2yz_S_D2z_S_cc;
  Double I_ERI_Fy2z_Pz_D2z_S_cc = I_ERI_Gy3z_S_D2z_S_cc+ABZ*I_ERI_Fy2z_S_D2z_S_cc;
  Double I_ERI_F3z_Pz_D2z_S_cc = I_ERI_G4z_S_D2z_S_cc+ABZ*I_ERI_F3z_S_D2z_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_D_S_cc
   * expanding position: BRA2
   * code section is: HRR
   * totally 36 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_D_S_cc
   * RHS shell quartet name: SQ_ERI_G_S_D_S_cc
   ************************************************************/
  Double I_ERI_G4x_Px_D2x_S_cc = I_ERI_H5x_S_D2x_S_cc+ABX*I_ERI_G4x_S_D2x_S_cc;
  Double I_ERI_G3xy_Px_D2x_S_cc = I_ERI_H4xy_S_D2x_S_cc+ABX*I_ERI_G3xy_S_D2x_S_cc;
  Double I_ERI_G3xz_Px_D2x_S_cc = I_ERI_H4xz_S_D2x_S_cc+ABX*I_ERI_G3xz_S_D2x_S_cc;
  Double I_ERI_G2x2y_Px_D2x_S_cc = I_ERI_H3x2y_S_D2x_S_cc+ABX*I_ERI_G2x2y_S_D2x_S_cc;
  Double I_ERI_G2xyz_Px_D2x_S_cc = I_ERI_H3xyz_S_D2x_S_cc+ABX*I_ERI_G2xyz_S_D2x_S_cc;
  Double I_ERI_G2x2z_Px_D2x_S_cc = I_ERI_H3x2z_S_D2x_S_cc+ABX*I_ERI_G2x2z_S_D2x_S_cc;
  Double I_ERI_Gx3y_Px_D2x_S_cc = I_ERI_H2x3y_S_D2x_S_cc+ABX*I_ERI_Gx3y_S_D2x_S_cc;
  Double I_ERI_Gx2yz_Px_D2x_S_cc = I_ERI_H2x2yz_S_D2x_S_cc+ABX*I_ERI_Gx2yz_S_D2x_S_cc;
  Double I_ERI_Gxy2z_Px_D2x_S_cc = I_ERI_H2xy2z_S_D2x_S_cc+ABX*I_ERI_Gxy2z_S_D2x_S_cc;
  Double I_ERI_Gx3z_Px_D2x_S_cc = I_ERI_H2x3z_S_D2x_S_cc+ABX*I_ERI_Gx3z_S_D2x_S_cc;
  Double I_ERI_G4y_Px_D2x_S_cc = I_ERI_Hx4y_S_D2x_S_cc+ABX*I_ERI_G4y_S_D2x_S_cc;
  Double I_ERI_G3yz_Px_D2x_S_cc = I_ERI_Hx3yz_S_D2x_S_cc+ABX*I_ERI_G3yz_S_D2x_S_cc;
  Double I_ERI_G2y2z_Px_D2x_S_cc = I_ERI_Hx2y2z_S_D2x_S_cc+ABX*I_ERI_G2y2z_S_D2x_S_cc;
  Double I_ERI_Gy3z_Px_D2x_S_cc = I_ERI_Hxy3z_S_D2x_S_cc+ABX*I_ERI_Gy3z_S_D2x_S_cc;
  Double I_ERI_G4z_Px_D2x_S_cc = I_ERI_Hx4z_S_D2x_S_cc+ABX*I_ERI_G4z_S_D2x_S_cc;
  Double I_ERI_G3xy_Py_D2x_S_cc = I_ERI_H3x2y_S_D2x_S_cc+ABY*I_ERI_G3xy_S_D2x_S_cc;
  Double I_ERI_G3xz_Py_D2x_S_cc = I_ERI_H3xyz_S_D2x_S_cc+ABY*I_ERI_G3xz_S_D2x_S_cc;
  Double I_ERI_G2x2y_Py_D2x_S_cc = I_ERI_H2x3y_S_D2x_S_cc+ABY*I_ERI_G2x2y_S_D2x_S_cc;
  Double I_ERI_G2xyz_Py_D2x_S_cc = I_ERI_H2x2yz_S_D2x_S_cc+ABY*I_ERI_G2xyz_S_D2x_S_cc;
  Double I_ERI_G2x2z_Py_D2x_S_cc = I_ERI_H2xy2z_S_D2x_S_cc+ABY*I_ERI_G2x2z_S_D2x_S_cc;
  Double I_ERI_Gx3y_Py_D2x_S_cc = I_ERI_Hx4y_S_D2x_S_cc+ABY*I_ERI_Gx3y_S_D2x_S_cc;
  Double I_ERI_Gx2yz_Py_D2x_S_cc = I_ERI_Hx3yz_S_D2x_S_cc+ABY*I_ERI_Gx2yz_S_D2x_S_cc;
  Double I_ERI_Gxy2z_Py_D2x_S_cc = I_ERI_Hx2y2z_S_D2x_S_cc+ABY*I_ERI_Gxy2z_S_D2x_S_cc;
  Double I_ERI_Gx3z_Py_D2x_S_cc = I_ERI_Hxy3z_S_D2x_S_cc+ABY*I_ERI_Gx3z_S_D2x_S_cc;
  Double I_ERI_G4y_Py_D2x_S_cc = I_ERI_H5y_S_D2x_S_cc+ABY*I_ERI_G4y_S_D2x_S_cc;
  Double I_ERI_G3yz_Py_D2x_S_cc = I_ERI_H4yz_S_D2x_S_cc+ABY*I_ERI_G3yz_S_D2x_S_cc;
  Double I_ERI_G2y2z_Py_D2x_S_cc = I_ERI_H3y2z_S_D2x_S_cc+ABY*I_ERI_G2y2z_S_D2x_S_cc;
  Double I_ERI_Gy3z_Py_D2x_S_cc = I_ERI_H2y3z_S_D2x_S_cc+ABY*I_ERI_Gy3z_S_D2x_S_cc;
  Double I_ERI_G4z_Py_D2x_S_cc = I_ERI_Hy4z_S_D2x_S_cc+ABY*I_ERI_G4z_S_D2x_S_cc;
  Double I_ERI_G3xz_Pz_D2x_S_cc = I_ERI_H3x2z_S_D2x_S_cc+ABZ*I_ERI_G3xz_S_D2x_S_cc;
  Double I_ERI_G2xyz_Pz_D2x_S_cc = I_ERI_H2xy2z_S_D2x_S_cc+ABZ*I_ERI_G2xyz_S_D2x_S_cc;
  Double I_ERI_G2x2z_Pz_D2x_S_cc = I_ERI_H2x3z_S_D2x_S_cc+ABZ*I_ERI_G2x2z_S_D2x_S_cc;
  Double I_ERI_Gx2yz_Pz_D2x_S_cc = I_ERI_Hx2y2z_S_D2x_S_cc+ABZ*I_ERI_Gx2yz_S_D2x_S_cc;
  Double I_ERI_Gxy2z_Pz_D2x_S_cc = I_ERI_Hxy3z_S_D2x_S_cc+ABZ*I_ERI_Gxy2z_S_D2x_S_cc;
  Double I_ERI_Gx3z_Pz_D2x_S_cc = I_ERI_Hx4z_S_D2x_S_cc+ABZ*I_ERI_Gx3z_S_D2x_S_cc;
  Double I_ERI_G3yz_Pz_D2x_S_cc = I_ERI_H3y2z_S_D2x_S_cc+ABZ*I_ERI_G3yz_S_D2x_S_cc;
  Double I_ERI_G2y2z_Pz_D2x_S_cc = I_ERI_H2y3z_S_D2x_S_cc+ABZ*I_ERI_G2y2z_S_D2x_S_cc;
  Double I_ERI_Gy3z_Pz_D2x_S_cc = I_ERI_Hy4z_S_D2x_S_cc+ABZ*I_ERI_Gy3z_S_D2x_S_cc;
  Double I_ERI_G4z_Pz_D2x_S_cc = I_ERI_H5z_S_D2x_S_cc+ABZ*I_ERI_G4z_S_D2x_S_cc;
  Double I_ERI_G4x_Px_Dxy_S_cc = I_ERI_H5x_S_Dxy_S_cc+ABX*I_ERI_G4x_S_Dxy_S_cc;
  Double I_ERI_G3xy_Px_Dxy_S_cc = I_ERI_H4xy_S_Dxy_S_cc+ABX*I_ERI_G3xy_S_Dxy_S_cc;
  Double I_ERI_G3xz_Px_Dxy_S_cc = I_ERI_H4xz_S_Dxy_S_cc+ABX*I_ERI_G3xz_S_Dxy_S_cc;
  Double I_ERI_G2x2y_Px_Dxy_S_cc = I_ERI_H3x2y_S_Dxy_S_cc+ABX*I_ERI_G2x2y_S_Dxy_S_cc;
  Double I_ERI_G2xyz_Px_Dxy_S_cc = I_ERI_H3xyz_S_Dxy_S_cc+ABX*I_ERI_G2xyz_S_Dxy_S_cc;
  Double I_ERI_G2x2z_Px_Dxy_S_cc = I_ERI_H3x2z_S_Dxy_S_cc+ABX*I_ERI_G2x2z_S_Dxy_S_cc;
  Double I_ERI_Gx3y_Px_Dxy_S_cc = I_ERI_H2x3y_S_Dxy_S_cc+ABX*I_ERI_Gx3y_S_Dxy_S_cc;
  Double I_ERI_Gx2yz_Px_Dxy_S_cc = I_ERI_H2x2yz_S_Dxy_S_cc+ABX*I_ERI_Gx2yz_S_Dxy_S_cc;
  Double I_ERI_Gxy2z_Px_Dxy_S_cc = I_ERI_H2xy2z_S_Dxy_S_cc+ABX*I_ERI_Gxy2z_S_Dxy_S_cc;
  Double I_ERI_Gx3z_Px_Dxy_S_cc = I_ERI_H2x3z_S_Dxy_S_cc+ABX*I_ERI_Gx3z_S_Dxy_S_cc;
  Double I_ERI_G4y_Px_Dxy_S_cc = I_ERI_Hx4y_S_Dxy_S_cc+ABX*I_ERI_G4y_S_Dxy_S_cc;
  Double I_ERI_G3yz_Px_Dxy_S_cc = I_ERI_Hx3yz_S_Dxy_S_cc+ABX*I_ERI_G3yz_S_Dxy_S_cc;
  Double I_ERI_G2y2z_Px_Dxy_S_cc = I_ERI_Hx2y2z_S_Dxy_S_cc+ABX*I_ERI_G2y2z_S_Dxy_S_cc;
  Double I_ERI_Gy3z_Px_Dxy_S_cc = I_ERI_Hxy3z_S_Dxy_S_cc+ABX*I_ERI_Gy3z_S_Dxy_S_cc;
  Double I_ERI_G4z_Px_Dxy_S_cc = I_ERI_Hx4z_S_Dxy_S_cc+ABX*I_ERI_G4z_S_Dxy_S_cc;
  Double I_ERI_G3xy_Py_Dxy_S_cc = I_ERI_H3x2y_S_Dxy_S_cc+ABY*I_ERI_G3xy_S_Dxy_S_cc;
  Double I_ERI_G3xz_Py_Dxy_S_cc = I_ERI_H3xyz_S_Dxy_S_cc+ABY*I_ERI_G3xz_S_Dxy_S_cc;
  Double I_ERI_G2x2y_Py_Dxy_S_cc = I_ERI_H2x3y_S_Dxy_S_cc+ABY*I_ERI_G2x2y_S_Dxy_S_cc;
  Double I_ERI_G2xyz_Py_Dxy_S_cc = I_ERI_H2x2yz_S_Dxy_S_cc+ABY*I_ERI_G2xyz_S_Dxy_S_cc;
  Double I_ERI_G2x2z_Py_Dxy_S_cc = I_ERI_H2xy2z_S_Dxy_S_cc+ABY*I_ERI_G2x2z_S_Dxy_S_cc;
  Double I_ERI_Gx3y_Py_Dxy_S_cc = I_ERI_Hx4y_S_Dxy_S_cc+ABY*I_ERI_Gx3y_S_Dxy_S_cc;
  Double I_ERI_Gx2yz_Py_Dxy_S_cc = I_ERI_Hx3yz_S_Dxy_S_cc+ABY*I_ERI_Gx2yz_S_Dxy_S_cc;
  Double I_ERI_Gxy2z_Py_Dxy_S_cc = I_ERI_Hx2y2z_S_Dxy_S_cc+ABY*I_ERI_Gxy2z_S_Dxy_S_cc;
  Double I_ERI_Gx3z_Py_Dxy_S_cc = I_ERI_Hxy3z_S_Dxy_S_cc+ABY*I_ERI_Gx3z_S_Dxy_S_cc;
  Double I_ERI_G4y_Py_Dxy_S_cc = I_ERI_H5y_S_Dxy_S_cc+ABY*I_ERI_G4y_S_Dxy_S_cc;
  Double I_ERI_G3yz_Py_Dxy_S_cc = I_ERI_H4yz_S_Dxy_S_cc+ABY*I_ERI_G3yz_S_Dxy_S_cc;
  Double I_ERI_G2y2z_Py_Dxy_S_cc = I_ERI_H3y2z_S_Dxy_S_cc+ABY*I_ERI_G2y2z_S_Dxy_S_cc;
  Double I_ERI_Gy3z_Py_Dxy_S_cc = I_ERI_H2y3z_S_Dxy_S_cc+ABY*I_ERI_Gy3z_S_Dxy_S_cc;
  Double I_ERI_G4z_Py_Dxy_S_cc = I_ERI_Hy4z_S_Dxy_S_cc+ABY*I_ERI_G4z_S_Dxy_S_cc;
  Double I_ERI_G3xz_Pz_Dxy_S_cc = I_ERI_H3x2z_S_Dxy_S_cc+ABZ*I_ERI_G3xz_S_Dxy_S_cc;
  Double I_ERI_G2xyz_Pz_Dxy_S_cc = I_ERI_H2xy2z_S_Dxy_S_cc+ABZ*I_ERI_G2xyz_S_Dxy_S_cc;
  Double I_ERI_G2x2z_Pz_Dxy_S_cc = I_ERI_H2x3z_S_Dxy_S_cc+ABZ*I_ERI_G2x2z_S_Dxy_S_cc;
  Double I_ERI_Gx2yz_Pz_Dxy_S_cc = I_ERI_Hx2y2z_S_Dxy_S_cc+ABZ*I_ERI_Gx2yz_S_Dxy_S_cc;
  Double I_ERI_Gxy2z_Pz_Dxy_S_cc = I_ERI_Hxy3z_S_Dxy_S_cc+ABZ*I_ERI_Gxy2z_S_Dxy_S_cc;
  Double I_ERI_Gx3z_Pz_Dxy_S_cc = I_ERI_Hx4z_S_Dxy_S_cc+ABZ*I_ERI_Gx3z_S_Dxy_S_cc;
  Double I_ERI_G3yz_Pz_Dxy_S_cc = I_ERI_H3y2z_S_Dxy_S_cc+ABZ*I_ERI_G3yz_S_Dxy_S_cc;
  Double I_ERI_G2y2z_Pz_Dxy_S_cc = I_ERI_H2y3z_S_Dxy_S_cc+ABZ*I_ERI_G2y2z_S_Dxy_S_cc;
  Double I_ERI_Gy3z_Pz_Dxy_S_cc = I_ERI_Hy4z_S_Dxy_S_cc+ABZ*I_ERI_Gy3z_S_Dxy_S_cc;
  Double I_ERI_G4z_Pz_Dxy_S_cc = I_ERI_H5z_S_Dxy_S_cc+ABZ*I_ERI_G4z_S_Dxy_S_cc;
  Double I_ERI_G4x_Px_Dxz_S_cc = I_ERI_H5x_S_Dxz_S_cc+ABX*I_ERI_G4x_S_Dxz_S_cc;
  Double I_ERI_G3xy_Px_Dxz_S_cc = I_ERI_H4xy_S_Dxz_S_cc+ABX*I_ERI_G3xy_S_Dxz_S_cc;
  Double I_ERI_G3xz_Px_Dxz_S_cc = I_ERI_H4xz_S_Dxz_S_cc+ABX*I_ERI_G3xz_S_Dxz_S_cc;
  Double I_ERI_G2x2y_Px_Dxz_S_cc = I_ERI_H3x2y_S_Dxz_S_cc+ABX*I_ERI_G2x2y_S_Dxz_S_cc;
  Double I_ERI_G2xyz_Px_Dxz_S_cc = I_ERI_H3xyz_S_Dxz_S_cc+ABX*I_ERI_G2xyz_S_Dxz_S_cc;
  Double I_ERI_G2x2z_Px_Dxz_S_cc = I_ERI_H3x2z_S_Dxz_S_cc+ABX*I_ERI_G2x2z_S_Dxz_S_cc;
  Double I_ERI_Gx3y_Px_Dxz_S_cc = I_ERI_H2x3y_S_Dxz_S_cc+ABX*I_ERI_Gx3y_S_Dxz_S_cc;
  Double I_ERI_Gx2yz_Px_Dxz_S_cc = I_ERI_H2x2yz_S_Dxz_S_cc+ABX*I_ERI_Gx2yz_S_Dxz_S_cc;
  Double I_ERI_Gxy2z_Px_Dxz_S_cc = I_ERI_H2xy2z_S_Dxz_S_cc+ABX*I_ERI_Gxy2z_S_Dxz_S_cc;
  Double I_ERI_Gx3z_Px_Dxz_S_cc = I_ERI_H2x3z_S_Dxz_S_cc+ABX*I_ERI_Gx3z_S_Dxz_S_cc;
  Double I_ERI_G4y_Px_Dxz_S_cc = I_ERI_Hx4y_S_Dxz_S_cc+ABX*I_ERI_G4y_S_Dxz_S_cc;
  Double I_ERI_G3yz_Px_Dxz_S_cc = I_ERI_Hx3yz_S_Dxz_S_cc+ABX*I_ERI_G3yz_S_Dxz_S_cc;
  Double I_ERI_G2y2z_Px_Dxz_S_cc = I_ERI_Hx2y2z_S_Dxz_S_cc+ABX*I_ERI_G2y2z_S_Dxz_S_cc;
  Double I_ERI_Gy3z_Px_Dxz_S_cc = I_ERI_Hxy3z_S_Dxz_S_cc+ABX*I_ERI_Gy3z_S_Dxz_S_cc;
  Double I_ERI_G4z_Px_Dxz_S_cc = I_ERI_Hx4z_S_Dxz_S_cc+ABX*I_ERI_G4z_S_Dxz_S_cc;
  Double I_ERI_G3xy_Py_Dxz_S_cc = I_ERI_H3x2y_S_Dxz_S_cc+ABY*I_ERI_G3xy_S_Dxz_S_cc;
  Double I_ERI_G3xz_Py_Dxz_S_cc = I_ERI_H3xyz_S_Dxz_S_cc+ABY*I_ERI_G3xz_S_Dxz_S_cc;
  Double I_ERI_G2x2y_Py_Dxz_S_cc = I_ERI_H2x3y_S_Dxz_S_cc+ABY*I_ERI_G2x2y_S_Dxz_S_cc;
  Double I_ERI_G2xyz_Py_Dxz_S_cc = I_ERI_H2x2yz_S_Dxz_S_cc+ABY*I_ERI_G2xyz_S_Dxz_S_cc;
  Double I_ERI_G2x2z_Py_Dxz_S_cc = I_ERI_H2xy2z_S_Dxz_S_cc+ABY*I_ERI_G2x2z_S_Dxz_S_cc;
  Double I_ERI_Gx3y_Py_Dxz_S_cc = I_ERI_Hx4y_S_Dxz_S_cc+ABY*I_ERI_Gx3y_S_Dxz_S_cc;
  Double I_ERI_Gx2yz_Py_Dxz_S_cc = I_ERI_Hx3yz_S_Dxz_S_cc+ABY*I_ERI_Gx2yz_S_Dxz_S_cc;
  Double I_ERI_Gxy2z_Py_Dxz_S_cc = I_ERI_Hx2y2z_S_Dxz_S_cc+ABY*I_ERI_Gxy2z_S_Dxz_S_cc;
  Double I_ERI_Gx3z_Py_Dxz_S_cc = I_ERI_Hxy3z_S_Dxz_S_cc+ABY*I_ERI_Gx3z_S_Dxz_S_cc;
  Double I_ERI_G4y_Py_Dxz_S_cc = I_ERI_H5y_S_Dxz_S_cc+ABY*I_ERI_G4y_S_Dxz_S_cc;
  Double I_ERI_G3yz_Py_Dxz_S_cc = I_ERI_H4yz_S_Dxz_S_cc+ABY*I_ERI_G3yz_S_Dxz_S_cc;
  Double I_ERI_G2y2z_Py_Dxz_S_cc = I_ERI_H3y2z_S_Dxz_S_cc+ABY*I_ERI_G2y2z_S_Dxz_S_cc;
  Double I_ERI_Gy3z_Py_Dxz_S_cc = I_ERI_H2y3z_S_Dxz_S_cc+ABY*I_ERI_Gy3z_S_Dxz_S_cc;
  Double I_ERI_G4z_Py_Dxz_S_cc = I_ERI_Hy4z_S_Dxz_S_cc+ABY*I_ERI_G4z_S_Dxz_S_cc;
  Double I_ERI_G3xz_Pz_Dxz_S_cc = I_ERI_H3x2z_S_Dxz_S_cc+ABZ*I_ERI_G3xz_S_Dxz_S_cc;
  Double I_ERI_G2xyz_Pz_Dxz_S_cc = I_ERI_H2xy2z_S_Dxz_S_cc+ABZ*I_ERI_G2xyz_S_Dxz_S_cc;
  Double I_ERI_G2x2z_Pz_Dxz_S_cc = I_ERI_H2x3z_S_Dxz_S_cc+ABZ*I_ERI_G2x2z_S_Dxz_S_cc;
  Double I_ERI_Gx2yz_Pz_Dxz_S_cc = I_ERI_Hx2y2z_S_Dxz_S_cc+ABZ*I_ERI_Gx2yz_S_Dxz_S_cc;
  Double I_ERI_Gxy2z_Pz_Dxz_S_cc = I_ERI_Hxy3z_S_Dxz_S_cc+ABZ*I_ERI_Gxy2z_S_Dxz_S_cc;
  Double I_ERI_Gx3z_Pz_Dxz_S_cc = I_ERI_Hx4z_S_Dxz_S_cc+ABZ*I_ERI_Gx3z_S_Dxz_S_cc;
  Double I_ERI_G3yz_Pz_Dxz_S_cc = I_ERI_H3y2z_S_Dxz_S_cc+ABZ*I_ERI_G3yz_S_Dxz_S_cc;
  Double I_ERI_G2y2z_Pz_Dxz_S_cc = I_ERI_H2y3z_S_Dxz_S_cc+ABZ*I_ERI_G2y2z_S_Dxz_S_cc;
  Double I_ERI_Gy3z_Pz_Dxz_S_cc = I_ERI_Hy4z_S_Dxz_S_cc+ABZ*I_ERI_Gy3z_S_Dxz_S_cc;
  Double I_ERI_G4z_Pz_Dxz_S_cc = I_ERI_H5z_S_Dxz_S_cc+ABZ*I_ERI_G4z_S_Dxz_S_cc;
  Double I_ERI_G4x_Px_D2y_S_cc = I_ERI_H5x_S_D2y_S_cc+ABX*I_ERI_G4x_S_D2y_S_cc;
  Double I_ERI_G3xy_Px_D2y_S_cc = I_ERI_H4xy_S_D2y_S_cc+ABX*I_ERI_G3xy_S_D2y_S_cc;
  Double I_ERI_G3xz_Px_D2y_S_cc = I_ERI_H4xz_S_D2y_S_cc+ABX*I_ERI_G3xz_S_D2y_S_cc;
  Double I_ERI_G2x2y_Px_D2y_S_cc = I_ERI_H3x2y_S_D2y_S_cc+ABX*I_ERI_G2x2y_S_D2y_S_cc;
  Double I_ERI_G2xyz_Px_D2y_S_cc = I_ERI_H3xyz_S_D2y_S_cc+ABX*I_ERI_G2xyz_S_D2y_S_cc;
  Double I_ERI_G2x2z_Px_D2y_S_cc = I_ERI_H3x2z_S_D2y_S_cc+ABX*I_ERI_G2x2z_S_D2y_S_cc;
  Double I_ERI_Gx3y_Px_D2y_S_cc = I_ERI_H2x3y_S_D2y_S_cc+ABX*I_ERI_Gx3y_S_D2y_S_cc;
  Double I_ERI_Gx2yz_Px_D2y_S_cc = I_ERI_H2x2yz_S_D2y_S_cc+ABX*I_ERI_Gx2yz_S_D2y_S_cc;
  Double I_ERI_Gxy2z_Px_D2y_S_cc = I_ERI_H2xy2z_S_D2y_S_cc+ABX*I_ERI_Gxy2z_S_D2y_S_cc;
  Double I_ERI_Gx3z_Px_D2y_S_cc = I_ERI_H2x3z_S_D2y_S_cc+ABX*I_ERI_Gx3z_S_D2y_S_cc;
  Double I_ERI_G4y_Px_D2y_S_cc = I_ERI_Hx4y_S_D2y_S_cc+ABX*I_ERI_G4y_S_D2y_S_cc;
  Double I_ERI_G3yz_Px_D2y_S_cc = I_ERI_Hx3yz_S_D2y_S_cc+ABX*I_ERI_G3yz_S_D2y_S_cc;
  Double I_ERI_G2y2z_Px_D2y_S_cc = I_ERI_Hx2y2z_S_D2y_S_cc+ABX*I_ERI_G2y2z_S_D2y_S_cc;
  Double I_ERI_Gy3z_Px_D2y_S_cc = I_ERI_Hxy3z_S_D2y_S_cc+ABX*I_ERI_Gy3z_S_D2y_S_cc;
  Double I_ERI_G4z_Px_D2y_S_cc = I_ERI_Hx4z_S_D2y_S_cc+ABX*I_ERI_G4z_S_D2y_S_cc;
  Double I_ERI_G3xy_Py_D2y_S_cc = I_ERI_H3x2y_S_D2y_S_cc+ABY*I_ERI_G3xy_S_D2y_S_cc;
  Double I_ERI_G3xz_Py_D2y_S_cc = I_ERI_H3xyz_S_D2y_S_cc+ABY*I_ERI_G3xz_S_D2y_S_cc;
  Double I_ERI_G2x2y_Py_D2y_S_cc = I_ERI_H2x3y_S_D2y_S_cc+ABY*I_ERI_G2x2y_S_D2y_S_cc;
  Double I_ERI_G2xyz_Py_D2y_S_cc = I_ERI_H2x2yz_S_D2y_S_cc+ABY*I_ERI_G2xyz_S_D2y_S_cc;
  Double I_ERI_G2x2z_Py_D2y_S_cc = I_ERI_H2xy2z_S_D2y_S_cc+ABY*I_ERI_G2x2z_S_D2y_S_cc;
  Double I_ERI_Gx3y_Py_D2y_S_cc = I_ERI_Hx4y_S_D2y_S_cc+ABY*I_ERI_Gx3y_S_D2y_S_cc;
  Double I_ERI_Gx2yz_Py_D2y_S_cc = I_ERI_Hx3yz_S_D2y_S_cc+ABY*I_ERI_Gx2yz_S_D2y_S_cc;
  Double I_ERI_Gxy2z_Py_D2y_S_cc = I_ERI_Hx2y2z_S_D2y_S_cc+ABY*I_ERI_Gxy2z_S_D2y_S_cc;
  Double I_ERI_Gx3z_Py_D2y_S_cc = I_ERI_Hxy3z_S_D2y_S_cc+ABY*I_ERI_Gx3z_S_D2y_S_cc;
  Double I_ERI_G4y_Py_D2y_S_cc = I_ERI_H5y_S_D2y_S_cc+ABY*I_ERI_G4y_S_D2y_S_cc;
  Double I_ERI_G3yz_Py_D2y_S_cc = I_ERI_H4yz_S_D2y_S_cc+ABY*I_ERI_G3yz_S_D2y_S_cc;
  Double I_ERI_G2y2z_Py_D2y_S_cc = I_ERI_H3y2z_S_D2y_S_cc+ABY*I_ERI_G2y2z_S_D2y_S_cc;
  Double I_ERI_Gy3z_Py_D2y_S_cc = I_ERI_H2y3z_S_D2y_S_cc+ABY*I_ERI_Gy3z_S_D2y_S_cc;
  Double I_ERI_G4z_Py_D2y_S_cc = I_ERI_Hy4z_S_D2y_S_cc+ABY*I_ERI_G4z_S_D2y_S_cc;
  Double I_ERI_G3xz_Pz_D2y_S_cc = I_ERI_H3x2z_S_D2y_S_cc+ABZ*I_ERI_G3xz_S_D2y_S_cc;
  Double I_ERI_G2xyz_Pz_D2y_S_cc = I_ERI_H2xy2z_S_D2y_S_cc+ABZ*I_ERI_G2xyz_S_D2y_S_cc;
  Double I_ERI_G2x2z_Pz_D2y_S_cc = I_ERI_H2x3z_S_D2y_S_cc+ABZ*I_ERI_G2x2z_S_D2y_S_cc;
  Double I_ERI_Gx2yz_Pz_D2y_S_cc = I_ERI_Hx2y2z_S_D2y_S_cc+ABZ*I_ERI_Gx2yz_S_D2y_S_cc;
  Double I_ERI_Gxy2z_Pz_D2y_S_cc = I_ERI_Hxy3z_S_D2y_S_cc+ABZ*I_ERI_Gxy2z_S_D2y_S_cc;
  Double I_ERI_Gx3z_Pz_D2y_S_cc = I_ERI_Hx4z_S_D2y_S_cc+ABZ*I_ERI_Gx3z_S_D2y_S_cc;
  Double I_ERI_G3yz_Pz_D2y_S_cc = I_ERI_H3y2z_S_D2y_S_cc+ABZ*I_ERI_G3yz_S_D2y_S_cc;
  Double I_ERI_G2y2z_Pz_D2y_S_cc = I_ERI_H2y3z_S_D2y_S_cc+ABZ*I_ERI_G2y2z_S_D2y_S_cc;
  Double I_ERI_Gy3z_Pz_D2y_S_cc = I_ERI_Hy4z_S_D2y_S_cc+ABZ*I_ERI_Gy3z_S_D2y_S_cc;
  Double I_ERI_G4z_Pz_D2y_S_cc = I_ERI_H5z_S_D2y_S_cc+ABZ*I_ERI_G4z_S_D2y_S_cc;
  Double I_ERI_G4x_Px_Dyz_S_cc = I_ERI_H5x_S_Dyz_S_cc+ABX*I_ERI_G4x_S_Dyz_S_cc;
  Double I_ERI_G3xy_Px_Dyz_S_cc = I_ERI_H4xy_S_Dyz_S_cc+ABX*I_ERI_G3xy_S_Dyz_S_cc;
  Double I_ERI_G3xz_Px_Dyz_S_cc = I_ERI_H4xz_S_Dyz_S_cc+ABX*I_ERI_G3xz_S_Dyz_S_cc;
  Double I_ERI_G2x2y_Px_Dyz_S_cc = I_ERI_H3x2y_S_Dyz_S_cc+ABX*I_ERI_G2x2y_S_Dyz_S_cc;
  Double I_ERI_G2xyz_Px_Dyz_S_cc = I_ERI_H3xyz_S_Dyz_S_cc+ABX*I_ERI_G2xyz_S_Dyz_S_cc;
  Double I_ERI_G2x2z_Px_Dyz_S_cc = I_ERI_H3x2z_S_Dyz_S_cc+ABX*I_ERI_G2x2z_S_Dyz_S_cc;
  Double I_ERI_Gx3y_Px_Dyz_S_cc = I_ERI_H2x3y_S_Dyz_S_cc+ABX*I_ERI_Gx3y_S_Dyz_S_cc;
  Double I_ERI_Gx2yz_Px_Dyz_S_cc = I_ERI_H2x2yz_S_Dyz_S_cc+ABX*I_ERI_Gx2yz_S_Dyz_S_cc;
  Double I_ERI_Gxy2z_Px_Dyz_S_cc = I_ERI_H2xy2z_S_Dyz_S_cc+ABX*I_ERI_Gxy2z_S_Dyz_S_cc;
  Double I_ERI_Gx3z_Px_Dyz_S_cc = I_ERI_H2x3z_S_Dyz_S_cc+ABX*I_ERI_Gx3z_S_Dyz_S_cc;
  Double I_ERI_G4y_Px_Dyz_S_cc = I_ERI_Hx4y_S_Dyz_S_cc+ABX*I_ERI_G4y_S_Dyz_S_cc;
  Double I_ERI_G3yz_Px_Dyz_S_cc = I_ERI_Hx3yz_S_Dyz_S_cc+ABX*I_ERI_G3yz_S_Dyz_S_cc;
  Double I_ERI_G2y2z_Px_Dyz_S_cc = I_ERI_Hx2y2z_S_Dyz_S_cc+ABX*I_ERI_G2y2z_S_Dyz_S_cc;
  Double I_ERI_Gy3z_Px_Dyz_S_cc = I_ERI_Hxy3z_S_Dyz_S_cc+ABX*I_ERI_Gy3z_S_Dyz_S_cc;
  Double I_ERI_G4z_Px_Dyz_S_cc = I_ERI_Hx4z_S_Dyz_S_cc+ABX*I_ERI_G4z_S_Dyz_S_cc;
  Double I_ERI_G3xy_Py_Dyz_S_cc = I_ERI_H3x2y_S_Dyz_S_cc+ABY*I_ERI_G3xy_S_Dyz_S_cc;
  Double I_ERI_G3xz_Py_Dyz_S_cc = I_ERI_H3xyz_S_Dyz_S_cc+ABY*I_ERI_G3xz_S_Dyz_S_cc;
  Double I_ERI_G2x2y_Py_Dyz_S_cc = I_ERI_H2x3y_S_Dyz_S_cc+ABY*I_ERI_G2x2y_S_Dyz_S_cc;
  Double I_ERI_G2xyz_Py_Dyz_S_cc = I_ERI_H2x2yz_S_Dyz_S_cc+ABY*I_ERI_G2xyz_S_Dyz_S_cc;
  Double I_ERI_G2x2z_Py_Dyz_S_cc = I_ERI_H2xy2z_S_Dyz_S_cc+ABY*I_ERI_G2x2z_S_Dyz_S_cc;
  Double I_ERI_Gx3y_Py_Dyz_S_cc = I_ERI_Hx4y_S_Dyz_S_cc+ABY*I_ERI_Gx3y_S_Dyz_S_cc;
  Double I_ERI_Gx2yz_Py_Dyz_S_cc = I_ERI_Hx3yz_S_Dyz_S_cc+ABY*I_ERI_Gx2yz_S_Dyz_S_cc;
  Double I_ERI_Gxy2z_Py_Dyz_S_cc = I_ERI_Hx2y2z_S_Dyz_S_cc+ABY*I_ERI_Gxy2z_S_Dyz_S_cc;
  Double I_ERI_Gx3z_Py_Dyz_S_cc = I_ERI_Hxy3z_S_Dyz_S_cc+ABY*I_ERI_Gx3z_S_Dyz_S_cc;
  Double I_ERI_G4y_Py_Dyz_S_cc = I_ERI_H5y_S_Dyz_S_cc+ABY*I_ERI_G4y_S_Dyz_S_cc;
  Double I_ERI_G3yz_Py_Dyz_S_cc = I_ERI_H4yz_S_Dyz_S_cc+ABY*I_ERI_G3yz_S_Dyz_S_cc;
  Double I_ERI_G2y2z_Py_Dyz_S_cc = I_ERI_H3y2z_S_Dyz_S_cc+ABY*I_ERI_G2y2z_S_Dyz_S_cc;
  Double I_ERI_Gy3z_Py_Dyz_S_cc = I_ERI_H2y3z_S_Dyz_S_cc+ABY*I_ERI_Gy3z_S_Dyz_S_cc;
  Double I_ERI_G4z_Py_Dyz_S_cc = I_ERI_Hy4z_S_Dyz_S_cc+ABY*I_ERI_G4z_S_Dyz_S_cc;
  Double I_ERI_G3xz_Pz_Dyz_S_cc = I_ERI_H3x2z_S_Dyz_S_cc+ABZ*I_ERI_G3xz_S_Dyz_S_cc;
  Double I_ERI_G2xyz_Pz_Dyz_S_cc = I_ERI_H2xy2z_S_Dyz_S_cc+ABZ*I_ERI_G2xyz_S_Dyz_S_cc;
  Double I_ERI_G2x2z_Pz_Dyz_S_cc = I_ERI_H2x3z_S_Dyz_S_cc+ABZ*I_ERI_G2x2z_S_Dyz_S_cc;
  Double I_ERI_Gx2yz_Pz_Dyz_S_cc = I_ERI_Hx2y2z_S_Dyz_S_cc+ABZ*I_ERI_Gx2yz_S_Dyz_S_cc;
  Double I_ERI_Gxy2z_Pz_Dyz_S_cc = I_ERI_Hxy3z_S_Dyz_S_cc+ABZ*I_ERI_Gxy2z_S_Dyz_S_cc;
  Double I_ERI_Gx3z_Pz_Dyz_S_cc = I_ERI_Hx4z_S_Dyz_S_cc+ABZ*I_ERI_Gx3z_S_Dyz_S_cc;
  Double I_ERI_G3yz_Pz_Dyz_S_cc = I_ERI_H3y2z_S_Dyz_S_cc+ABZ*I_ERI_G3yz_S_Dyz_S_cc;
  Double I_ERI_G2y2z_Pz_Dyz_S_cc = I_ERI_H2y3z_S_Dyz_S_cc+ABZ*I_ERI_G2y2z_S_Dyz_S_cc;
  Double I_ERI_Gy3z_Pz_Dyz_S_cc = I_ERI_Hy4z_S_Dyz_S_cc+ABZ*I_ERI_Gy3z_S_Dyz_S_cc;
  Double I_ERI_G4z_Pz_Dyz_S_cc = I_ERI_H5z_S_Dyz_S_cc+ABZ*I_ERI_G4z_S_Dyz_S_cc;
  Double I_ERI_G4x_Px_D2z_S_cc = I_ERI_H5x_S_D2z_S_cc+ABX*I_ERI_G4x_S_D2z_S_cc;
  Double I_ERI_G3xy_Px_D2z_S_cc = I_ERI_H4xy_S_D2z_S_cc+ABX*I_ERI_G3xy_S_D2z_S_cc;
  Double I_ERI_G3xz_Px_D2z_S_cc = I_ERI_H4xz_S_D2z_S_cc+ABX*I_ERI_G3xz_S_D2z_S_cc;
  Double I_ERI_G2x2y_Px_D2z_S_cc = I_ERI_H3x2y_S_D2z_S_cc+ABX*I_ERI_G2x2y_S_D2z_S_cc;
  Double I_ERI_G2xyz_Px_D2z_S_cc = I_ERI_H3xyz_S_D2z_S_cc+ABX*I_ERI_G2xyz_S_D2z_S_cc;
  Double I_ERI_G2x2z_Px_D2z_S_cc = I_ERI_H3x2z_S_D2z_S_cc+ABX*I_ERI_G2x2z_S_D2z_S_cc;
  Double I_ERI_Gx3y_Px_D2z_S_cc = I_ERI_H2x3y_S_D2z_S_cc+ABX*I_ERI_Gx3y_S_D2z_S_cc;
  Double I_ERI_Gx2yz_Px_D2z_S_cc = I_ERI_H2x2yz_S_D2z_S_cc+ABX*I_ERI_Gx2yz_S_D2z_S_cc;
  Double I_ERI_Gxy2z_Px_D2z_S_cc = I_ERI_H2xy2z_S_D2z_S_cc+ABX*I_ERI_Gxy2z_S_D2z_S_cc;
  Double I_ERI_Gx3z_Px_D2z_S_cc = I_ERI_H2x3z_S_D2z_S_cc+ABX*I_ERI_Gx3z_S_D2z_S_cc;
  Double I_ERI_G4y_Px_D2z_S_cc = I_ERI_Hx4y_S_D2z_S_cc+ABX*I_ERI_G4y_S_D2z_S_cc;
  Double I_ERI_G3yz_Px_D2z_S_cc = I_ERI_Hx3yz_S_D2z_S_cc+ABX*I_ERI_G3yz_S_D2z_S_cc;
  Double I_ERI_G2y2z_Px_D2z_S_cc = I_ERI_Hx2y2z_S_D2z_S_cc+ABX*I_ERI_G2y2z_S_D2z_S_cc;
  Double I_ERI_Gy3z_Px_D2z_S_cc = I_ERI_Hxy3z_S_D2z_S_cc+ABX*I_ERI_Gy3z_S_D2z_S_cc;
  Double I_ERI_G4z_Px_D2z_S_cc = I_ERI_Hx4z_S_D2z_S_cc+ABX*I_ERI_G4z_S_D2z_S_cc;
  Double I_ERI_G3xy_Py_D2z_S_cc = I_ERI_H3x2y_S_D2z_S_cc+ABY*I_ERI_G3xy_S_D2z_S_cc;
  Double I_ERI_G3xz_Py_D2z_S_cc = I_ERI_H3xyz_S_D2z_S_cc+ABY*I_ERI_G3xz_S_D2z_S_cc;
  Double I_ERI_G2x2y_Py_D2z_S_cc = I_ERI_H2x3y_S_D2z_S_cc+ABY*I_ERI_G2x2y_S_D2z_S_cc;
  Double I_ERI_G2xyz_Py_D2z_S_cc = I_ERI_H2x2yz_S_D2z_S_cc+ABY*I_ERI_G2xyz_S_D2z_S_cc;
  Double I_ERI_G2x2z_Py_D2z_S_cc = I_ERI_H2xy2z_S_D2z_S_cc+ABY*I_ERI_G2x2z_S_D2z_S_cc;
  Double I_ERI_Gx3y_Py_D2z_S_cc = I_ERI_Hx4y_S_D2z_S_cc+ABY*I_ERI_Gx3y_S_D2z_S_cc;
  Double I_ERI_Gx2yz_Py_D2z_S_cc = I_ERI_Hx3yz_S_D2z_S_cc+ABY*I_ERI_Gx2yz_S_D2z_S_cc;
  Double I_ERI_Gxy2z_Py_D2z_S_cc = I_ERI_Hx2y2z_S_D2z_S_cc+ABY*I_ERI_Gxy2z_S_D2z_S_cc;
  Double I_ERI_Gx3z_Py_D2z_S_cc = I_ERI_Hxy3z_S_D2z_S_cc+ABY*I_ERI_Gx3z_S_D2z_S_cc;
  Double I_ERI_G4y_Py_D2z_S_cc = I_ERI_H5y_S_D2z_S_cc+ABY*I_ERI_G4y_S_D2z_S_cc;
  Double I_ERI_G3yz_Py_D2z_S_cc = I_ERI_H4yz_S_D2z_S_cc+ABY*I_ERI_G3yz_S_D2z_S_cc;
  Double I_ERI_G2y2z_Py_D2z_S_cc = I_ERI_H3y2z_S_D2z_S_cc+ABY*I_ERI_G2y2z_S_D2z_S_cc;
  Double I_ERI_Gy3z_Py_D2z_S_cc = I_ERI_H2y3z_S_D2z_S_cc+ABY*I_ERI_Gy3z_S_D2z_S_cc;
  Double I_ERI_G4z_Py_D2z_S_cc = I_ERI_Hy4z_S_D2z_S_cc+ABY*I_ERI_G4z_S_D2z_S_cc;
  Double I_ERI_G3xz_Pz_D2z_S_cc = I_ERI_H3x2z_S_D2z_S_cc+ABZ*I_ERI_G3xz_S_D2z_S_cc;
  Double I_ERI_G2xyz_Pz_D2z_S_cc = I_ERI_H2xy2z_S_D2z_S_cc+ABZ*I_ERI_G2xyz_S_D2z_S_cc;
  Double I_ERI_G2x2z_Pz_D2z_S_cc = I_ERI_H2x3z_S_D2z_S_cc+ABZ*I_ERI_G2x2z_S_D2z_S_cc;
  Double I_ERI_Gx2yz_Pz_D2z_S_cc = I_ERI_Hx2y2z_S_D2z_S_cc+ABZ*I_ERI_Gx2yz_S_D2z_S_cc;
  Double I_ERI_Gxy2z_Pz_D2z_S_cc = I_ERI_Hxy3z_S_D2z_S_cc+ABZ*I_ERI_Gxy2z_S_D2z_S_cc;
  Double I_ERI_Gx3z_Pz_D2z_S_cc = I_ERI_Hx4z_S_D2z_S_cc+ABZ*I_ERI_Gx3z_S_D2z_S_cc;
  Double I_ERI_G3yz_Pz_D2z_S_cc = I_ERI_H3y2z_S_D2z_S_cc+ABZ*I_ERI_G3yz_S_D2z_S_cc;
  Double I_ERI_G2y2z_Pz_D2z_S_cc = I_ERI_H2y3z_S_D2z_S_cc+ABZ*I_ERI_G2y2z_S_D2z_S_cc;
  Double I_ERI_Gy3z_Pz_D2z_S_cc = I_ERI_Hy4z_S_D2z_S_cc+ABZ*I_ERI_Gy3z_S_D2z_S_cc;
  Double I_ERI_G4z_Pz_D2z_S_cc = I_ERI_H5z_S_D2z_S_cc+ABZ*I_ERI_G4z_S_D2z_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_D_S_cc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_D_S_cc
   * RHS shell quartet name: SQ_ERI_F_P_D_S_cc
   ************************************************************/
  Double I_ERI_F3x_D2x_D2x_S_cc = I_ERI_G4x_Px_D2x_S_cc+ABX*I_ERI_F3x_Px_D2x_S_cc;
  Double I_ERI_F2xy_D2x_D2x_S_cc = I_ERI_G3xy_Px_D2x_S_cc+ABX*I_ERI_F2xy_Px_D2x_S_cc;
  Double I_ERI_F2xz_D2x_D2x_S_cc = I_ERI_G3xz_Px_D2x_S_cc+ABX*I_ERI_F2xz_Px_D2x_S_cc;
  Double I_ERI_Fx2y_D2x_D2x_S_cc = I_ERI_G2x2y_Px_D2x_S_cc+ABX*I_ERI_Fx2y_Px_D2x_S_cc;
  Double I_ERI_Fxyz_D2x_D2x_S_cc = I_ERI_G2xyz_Px_D2x_S_cc+ABX*I_ERI_Fxyz_Px_D2x_S_cc;
  Double I_ERI_Fx2z_D2x_D2x_S_cc = I_ERI_G2x2z_Px_D2x_S_cc+ABX*I_ERI_Fx2z_Px_D2x_S_cc;
  Double I_ERI_F3y_D2x_D2x_S_cc = I_ERI_Gx3y_Px_D2x_S_cc+ABX*I_ERI_F3y_Px_D2x_S_cc;
  Double I_ERI_F2yz_D2x_D2x_S_cc = I_ERI_Gx2yz_Px_D2x_S_cc+ABX*I_ERI_F2yz_Px_D2x_S_cc;
  Double I_ERI_Fy2z_D2x_D2x_S_cc = I_ERI_Gxy2z_Px_D2x_S_cc+ABX*I_ERI_Fy2z_Px_D2x_S_cc;
  Double I_ERI_F3z_D2x_D2x_S_cc = I_ERI_Gx3z_Px_D2x_S_cc+ABX*I_ERI_F3z_Px_D2x_S_cc;
  Double I_ERI_F3x_Dxy_D2x_S_cc = I_ERI_G3xy_Px_D2x_S_cc+ABY*I_ERI_F3x_Px_D2x_S_cc;
  Double I_ERI_F2xy_Dxy_D2x_S_cc = I_ERI_G2x2y_Px_D2x_S_cc+ABY*I_ERI_F2xy_Px_D2x_S_cc;
  Double I_ERI_F2xz_Dxy_D2x_S_cc = I_ERI_G2xyz_Px_D2x_S_cc+ABY*I_ERI_F2xz_Px_D2x_S_cc;
  Double I_ERI_Fx2y_Dxy_D2x_S_cc = I_ERI_Gx3y_Px_D2x_S_cc+ABY*I_ERI_Fx2y_Px_D2x_S_cc;
  Double I_ERI_Fxyz_Dxy_D2x_S_cc = I_ERI_Gx2yz_Px_D2x_S_cc+ABY*I_ERI_Fxyz_Px_D2x_S_cc;
  Double I_ERI_Fx2z_Dxy_D2x_S_cc = I_ERI_Gxy2z_Px_D2x_S_cc+ABY*I_ERI_Fx2z_Px_D2x_S_cc;
  Double I_ERI_F3y_Dxy_D2x_S_cc = I_ERI_G4y_Px_D2x_S_cc+ABY*I_ERI_F3y_Px_D2x_S_cc;
  Double I_ERI_F2yz_Dxy_D2x_S_cc = I_ERI_G3yz_Px_D2x_S_cc+ABY*I_ERI_F2yz_Px_D2x_S_cc;
  Double I_ERI_Fy2z_Dxy_D2x_S_cc = I_ERI_G2y2z_Px_D2x_S_cc+ABY*I_ERI_Fy2z_Px_D2x_S_cc;
  Double I_ERI_F3z_Dxy_D2x_S_cc = I_ERI_Gy3z_Px_D2x_S_cc+ABY*I_ERI_F3z_Px_D2x_S_cc;
  Double I_ERI_F3x_Dxz_D2x_S_cc = I_ERI_G3xz_Px_D2x_S_cc+ABZ*I_ERI_F3x_Px_D2x_S_cc;
  Double I_ERI_F2xy_Dxz_D2x_S_cc = I_ERI_G2xyz_Px_D2x_S_cc+ABZ*I_ERI_F2xy_Px_D2x_S_cc;
  Double I_ERI_F2xz_Dxz_D2x_S_cc = I_ERI_G2x2z_Px_D2x_S_cc+ABZ*I_ERI_F2xz_Px_D2x_S_cc;
  Double I_ERI_Fx2y_Dxz_D2x_S_cc = I_ERI_Gx2yz_Px_D2x_S_cc+ABZ*I_ERI_Fx2y_Px_D2x_S_cc;
  Double I_ERI_Fxyz_Dxz_D2x_S_cc = I_ERI_Gxy2z_Px_D2x_S_cc+ABZ*I_ERI_Fxyz_Px_D2x_S_cc;
  Double I_ERI_Fx2z_Dxz_D2x_S_cc = I_ERI_Gx3z_Px_D2x_S_cc+ABZ*I_ERI_Fx2z_Px_D2x_S_cc;
  Double I_ERI_F3y_Dxz_D2x_S_cc = I_ERI_G3yz_Px_D2x_S_cc+ABZ*I_ERI_F3y_Px_D2x_S_cc;
  Double I_ERI_F2yz_Dxz_D2x_S_cc = I_ERI_G2y2z_Px_D2x_S_cc+ABZ*I_ERI_F2yz_Px_D2x_S_cc;
  Double I_ERI_Fy2z_Dxz_D2x_S_cc = I_ERI_Gy3z_Px_D2x_S_cc+ABZ*I_ERI_Fy2z_Px_D2x_S_cc;
  Double I_ERI_F3z_Dxz_D2x_S_cc = I_ERI_G4z_Px_D2x_S_cc+ABZ*I_ERI_F3z_Px_D2x_S_cc;
  Double I_ERI_F3x_D2y_D2x_S_cc = I_ERI_G3xy_Py_D2x_S_cc+ABY*I_ERI_F3x_Py_D2x_S_cc;
  Double I_ERI_F2xy_D2y_D2x_S_cc = I_ERI_G2x2y_Py_D2x_S_cc+ABY*I_ERI_F2xy_Py_D2x_S_cc;
  Double I_ERI_F2xz_D2y_D2x_S_cc = I_ERI_G2xyz_Py_D2x_S_cc+ABY*I_ERI_F2xz_Py_D2x_S_cc;
  Double I_ERI_Fx2y_D2y_D2x_S_cc = I_ERI_Gx3y_Py_D2x_S_cc+ABY*I_ERI_Fx2y_Py_D2x_S_cc;
  Double I_ERI_Fxyz_D2y_D2x_S_cc = I_ERI_Gx2yz_Py_D2x_S_cc+ABY*I_ERI_Fxyz_Py_D2x_S_cc;
  Double I_ERI_Fx2z_D2y_D2x_S_cc = I_ERI_Gxy2z_Py_D2x_S_cc+ABY*I_ERI_Fx2z_Py_D2x_S_cc;
  Double I_ERI_F3y_D2y_D2x_S_cc = I_ERI_G4y_Py_D2x_S_cc+ABY*I_ERI_F3y_Py_D2x_S_cc;
  Double I_ERI_F2yz_D2y_D2x_S_cc = I_ERI_G3yz_Py_D2x_S_cc+ABY*I_ERI_F2yz_Py_D2x_S_cc;
  Double I_ERI_Fy2z_D2y_D2x_S_cc = I_ERI_G2y2z_Py_D2x_S_cc+ABY*I_ERI_Fy2z_Py_D2x_S_cc;
  Double I_ERI_F3z_D2y_D2x_S_cc = I_ERI_Gy3z_Py_D2x_S_cc+ABY*I_ERI_F3z_Py_D2x_S_cc;
  Double I_ERI_F3x_Dyz_D2x_S_cc = I_ERI_G3xz_Py_D2x_S_cc+ABZ*I_ERI_F3x_Py_D2x_S_cc;
  Double I_ERI_F2xy_Dyz_D2x_S_cc = I_ERI_G2xyz_Py_D2x_S_cc+ABZ*I_ERI_F2xy_Py_D2x_S_cc;
  Double I_ERI_F2xz_Dyz_D2x_S_cc = I_ERI_G2x2z_Py_D2x_S_cc+ABZ*I_ERI_F2xz_Py_D2x_S_cc;
  Double I_ERI_Fx2y_Dyz_D2x_S_cc = I_ERI_Gx2yz_Py_D2x_S_cc+ABZ*I_ERI_Fx2y_Py_D2x_S_cc;
  Double I_ERI_Fxyz_Dyz_D2x_S_cc = I_ERI_Gxy2z_Py_D2x_S_cc+ABZ*I_ERI_Fxyz_Py_D2x_S_cc;
  Double I_ERI_Fx2z_Dyz_D2x_S_cc = I_ERI_Gx3z_Py_D2x_S_cc+ABZ*I_ERI_Fx2z_Py_D2x_S_cc;
  Double I_ERI_F3y_Dyz_D2x_S_cc = I_ERI_G3yz_Py_D2x_S_cc+ABZ*I_ERI_F3y_Py_D2x_S_cc;
  Double I_ERI_F2yz_Dyz_D2x_S_cc = I_ERI_G2y2z_Py_D2x_S_cc+ABZ*I_ERI_F2yz_Py_D2x_S_cc;
  Double I_ERI_Fy2z_Dyz_D2x_S_cc = I_ERI_Gy3z_Py_D2x_S_cc+ABZ*I_ERI_Fy2z_Py_D2x_S_cc;
  Double I_ERI_F3z_Dyz_D2x_S_cc = I_ERI_G4z_Py_D2x_S_cc+ABZ*I_ERI_F3z_Py_D2x_S_cc;
  Double I_ERI_F3x_D2z_D2x_S_cc = I_ERI_G3xz_Pz_D2x_S_cc+ABZ*I_ERI_F3x_Pz_D2x_S_cc;
  Double I_ERI_F2xy_D2z_D2x_S_cc = I_ERI_G2xyz_Pz_D2x_S_cc+ABZ*I_ERI_F2xy_Pz_D2x_S_cc;
  Double I_ERI_F2xz_D2z_D2x_S_cc = I_ERI_G2x2z_Pz_D2x_S_cc+ABZ*I_ERI_F2xz_Pz_D2x_S_cc;
  Double I_ERI_Fx2y_D2z_D2x_S_cc = I_ERI_Gx2yz_Pz_D2x_S_cc+ABZ*I_ERI_Fx2y_Pz_D2x_S_cc;
  Double I_ERI_Fxyz_D2z_D2x_S_cc = I_ERI_Gxy2z_Pz_D2x_S_cc+ABZ*I_ERI_Fxyz_Pz_D2x_S_cc;
  Double I_ERI_Fx2z_D2z_D2x_S_cc = I_ERI_Gx3z_Pz_D2x_S_cc+ABZ*I_ERI_Fx2z_Pz_D2x_S_cc;
  Double I_ERI_F3y_D2z_D2x_S_cc = I_ERI_G3yz_Pz_D2x_S_cc+ABZ*I_ERI_F3y_Pz_D2x_S_cc;
  Double I_ERI_F2yz_D2z_D2x_S_cc = I_ERI_G2y2z_Pz_D2x_S_cc+ABZ*I_ERI_F2yz_Pz_D2x_S_cc;
  Double I_ERI_Fy2z_D2z_D2x_S_cc = I_ERI_Gy3z_Pz_D2x_S_cc+ABZ*I_ERI_Fy2z_Pz_D2x_S_cc;
  Double I_ERI_F3z_D2z_D2x_S_cc = I_ERI_G4z_Pz_D2x_S_cc+ABZ*I_ERI_F3z_Pz_D2x_S_cc;
  Double I_ERI_F3x_D2x_Dxy_S_cc = I_ERI_G4x_Px_Dxy_S_cc+ABX*I_ERI_F3x_Px_Dxy_S_cc;
  Double I_ERI_F2xy_D2x_Dxy_S_cc = I_ERI_G3xy_Px_Dxy_S_cc+ABX*I_ERI_F2xy_Px_Dxy_S_cc;
  Double I_ERI_F2xz_D2x_Dxy_S_cc = I_ERI_G3xz_Px_Dxy_S_cc+ABX*I_ERI_F2xz_Px_Dxy_S_cc;
  Double I_ERI_Fx2y_D2x_Dxy_S_cc = I_ERI_G2x2y_Px_Dxy_S_cc+ABX*I_ERI_Fx2y_Px_Dxy_S_cc;
  Double I_ERI_Fxyz_D2x_Dxy_S_cc = I_ERI_G2xyz_Px_Dxy_S_cc+ABX*I_ERI_Fxyz_Px_Dxy_S_cc;
  Double I_ERI_Fx2z_D2x_Dxy_S_cc = I_ERI_G2x2z_Px_Dxy_S_cc+ABX*I_ERI_Fx2z_Px_Dxy_S_cc;
  Double I_ERI_F3y_D2x_Dxy_S_cc = I_ERI_Gx3y_Px_Dxy_S_cc+ABX*I_ERI_F3y_Px_Dxy_S_cc;
  Double I_ERI_F2yz_D2x_Dxy_S_cc = I_ERI_Gx2yz_Px_Dxy_S_cc+ABX*I_ERI_F2yz_Px_Dxy_S_cc;
  Double I_ERI_Fy2z_D2x_Dxy_S_cc = I_ERI_Gxy2z_Px_Dxy_S_cc+ABX*I_ERI_Fy2z_Px_Dxy_S_cc;
  Double I_ERI_F3z_D2x_Dxy_S_cc = I_ERI_Gx3z_Px_Dxy_S_cc+ABX*I_ERI_F3z_Px_Dxy_S_cc;
  Double I_ERI_F3x_Dxy_Dxy_S_cc = I_ERI_G3xy_Px_Dxy_S_cc+ABY*I_ERI_F3x_Px_Dxy_S_cc;
  Double I_ERI_F2xy_Dxy_Dxy_S_cc = I_ERI_G2x2y_Px_Dxy_S_cc+ABY*I_ERI_F2xy_Px_Dxy_S_cc;
  Double I_ERI_F2xz_Dxy_Dxy_S_cc = I_ERI_G2xyz_Px_Dxy_S_cc+ABY*I_ERI_F2xz_Px_Dxy_S_cc;
  Double I_ERI_Fx2y_Dxy_Dxy_S_cc = I_ERI_Gx3y_Px_Dxy_S_cc+ABY*I_ERI_Fx2y_Px_Dxy_S_cc;
  Double I_ERI_Fxyz_Dxy_Dxy_S_cc = I_ERI_Gx2yz_Px_Dxy_S_cc+ABY*I_ERI_Fxyz_Px_Dxy_S_cc;
  Double I_ERI_Fx2z_Dxy_Dxy_S_cc = I_ERI_Gxy2z_Px_Dxy_S_cc+ABY*I_ERI_Fx2z_Px_Dxy_S_cc;
  Double I_ERI_F3y_Dxy_Dxy_S_cc = I_ERI_G4y_Px_Dxy_S_cc+ABY*I_ERI_F3y_Px_Dxy_S_cc;
  Double I_ERI_F2yz_Dxy_Dxy_S_cc = I_ERI_G3yz_Px_Dxy_S_cc+ABY*I_ERI_F2yz_Px_Dxy_S_cc;
  Double I_ERI_Fy2z_Dxy_Dxy_S_cc = I_ERI_G2y2z_Px_Dxy_S_cc+ABY*I_ERI_Fy2z_Px_Dxy_S_cc;
  Double I_ERI_F3z_Dxy_Dxy_S_cc = I_ERI_Gy3z_Px_Dxy_S_cc+ABY*I_ERI_F3z_Px_Dxy_S_cc;
  Double I_ERI_F3x_Dxz_Dxy_S_cc = I_ERI_G3xz_Px_Dxy_S_cc+ABZ*I_ERI_F3x_Px_Dxy_S_cc;
  Double I_ERI_F2xy_Dxz_Dxy_S_cc = I_ERI_G2xyz_Px_Dxy_S_cc+ABZ*I_ERI_F2xy_Px_Dxy_S_cc;
  Double I_ERI_F2xz_Dxz_Dxy_S_cc = I_ERI_G2x2z_Px_Dxy_S_cc+ABZ*I_ERI_F2xz_Px_Dxy_S_cc;
  Double I_ERI_Fx2y_Dxz_Dxy_S_cc = I_ERI_Gx2yz_Px_Dxy_S_cc+ABZ*I_ERI_Fx2y_Px_Dxy_S_cc;
  Double I_ERI_Fxyz_Dxz_Dxy_S_cc = I_ERI_Gxy2z_Px_Dxy_S_cc+ABZ*I_ERI_Fxyz_Px_Dxy_S_cc;
  Double I_ERI_Fx2z_Dxz_Dxy_S_cc = I_ERI_Gx3z_Px_Dxy_S_cc+ABZ*I_ERI_Fx2z_Px_Dxy_S_cc;
  Double I_ERI_F3y_Dxz_Dxy_S_cc = I_ERI_G3yz_Px_Dxy_S_cc+ABZ*I_ERI_F3y_Px_Dxy_S_cc;
  Double I_ERI_F2yz_Dxz_Dxy_S_cc = I_ERI_G2y2z_Px_Dxy_S_cc+ABZ*I_ERI_F2yz_Px_Dxy_S_cc;
  Double I_ERI_Fy2z_Dxz_Dxy_S_cc = I_ERI_Gy3z_Px_Dxy_S_cc+ABZ*I_ERI_Fy2z_Px_Dxy_S_cc;
  Double I_ERI_F3z_Dxz_Dxy_S_cc = I_ERI_G4z_Px_Dxy_S_cc+ABZ*I_ERI_F3z_Px_Dxy_S_cc;
  Double I_ERI_F3x_D2y_Dxy_S_cc = I_ERI_G3xy_Py_Dxy_S_cc+ABY*I_ERI_F3x_Py_Dxy_S_cc;
  Double I_ERI_F2xy_D2y_Dxy_S_cc = I_ERI_G2x2y_Py_Dxy_S_cc+ABY*I_ERI_F2xy_Py_Dxy_S_cc;
  Double I_ERI_F2xz_D2y_Dxy_S_cc = I_ERI_G2xyz_Py_Dxy_S_cc+ABY*I_ERI_F2xz_Py_Dxy_S_cc;
  Double I_ERI_Fx2y_D2y_Dxy_S_cc = I_ERI_Gx3y_Py_Dxy_S_cc+ABY*I_ERI_Fx2y_Py_Dxy_S_cc;
  Double I_ERI_Fxyz_D2y_Dxy_S_cc = I_ERI_Gx2yz_Py_Dxy_S_cc+ABY*I_ERI_Fxyz_Py_Dxy_S_cc;
  Double I_ERI_Fx2z_D2y_Dxy_S_cc = I_ERI_Gxy2z_Py_Dxy_S_cc+ABY*I_ERI_Fx2z_Py_Dxy_S_cc;
  Double I_ERI_F3y_D2y_Dxy_S_cc = I_ERI_G4y_Py_Dxy_S_cc+ABY*I_ERI_F3y_Py_Dxy_S_cc;
  Double I_ERI_F2yz_D2y_Dxy_S_cc = I_ERI_G3yz_Py_Dxy_S_cc+ABY*I_ERI_F2yz_Py_Dxy_S_cc;
  Double I_ERI_Fy2z_D2y_Dxy_S_cc = I_ERI_G2y2z_Py_Dxy_S_cc+ABY*I_ERI_Fy2z_Py_Dxy_S_cc;
  Double I_ERI_F3z_D2y_Dxy_S_cc = I_ERI_Gy3z_Py_Dxy_S_cc+ABY*I_ERI_F3z_Py_Dxy_S_cc;
  Double I_ERI_F3x_Dyz_Dxy_S_cc = I_ERI_G3xz_Py_Dxy_S_cc+ABZ*I_ERI_F3x_Py_Dxy_S_cc;
  Double I_ERI_F2xy_Dyz_Dxy_S_cc = I_ERI_G2xyz_Py_Dxy_S_cc+ABZ*I_ERI_F2xy_Py_Dxy_S_cc;
  Double I_ERI_F2xz_Dyz_Dxy_S_cc = I_ERI_G2x2z_Py_Dxy_S_cc+ABZ*I_ERI_F2xz_Py_Dxy_S_cc;
  Double I_ERI_Fx2y_Dyz_Dxy_S_cc = I_ERI_Gx2yz_Py_Dxy_S_cc+ABZ*I_ERI_Fx2y_Py_Dxy_S_cc;
  Double I_ERI_Fxyz_Dyz_Dxy_S_cc = I_ERI_Gxy2z_Py_Dxy_S_cc+ABZ*I_ERI_Fxyz_Py_Dxy_S_cc;
  Double I_ERI_Fx2z_Dyz_Dxy_S_cc = I_ERI_Gx3z_Py_Dxy_S_cc+ABZ*I_ERI_Fx2z_Py_Dxy_S_cc;
  Double I_ERI_F3y_Dyz_Dxy_S_cc = I_ERI_G3yz_Py_Dxy_S_cc+ABZ*I_ERI_F3y_Py_Dxy_S_cc;
  Double I_ERI_F2yz_Dyz_Dxy_S_cc = I_ERI_G2y2z_Py_Dxy_S_cc+ABZ*I_ERI_F2yz_Py_Dxy_S_cc;
  Double I_ERI_Fy2z_Dyz_Dxy_S_cc = I_ERI_Gy3z_Py_Dxy_S_cc+ABZ*I_ERI_Fy2z_Py_Dxy_S_cc;
  Double I_ERI_F3z_Dyz_Dxy_S_cc = I_ERI_G4z_Py_Dxy_S_cc+ABZ*I_ERI_F3z_Py_Dxy_S_cc;
  Double I_ERI_F3x_D2z_Dxy_S_cc = I_ERI_G3xz_Pz_Dxy_S_cc+ABZ*I_ERI_F3x_Pz_Dxy_S_cc;
  Double I_ERI_F2xy_D2z_Dxy_S_cc = I_ERI_G2xyz_Pz_Dxy_S_cc+ABZ*I_ERI_F2xy_Pz_Dxy_S_cc;
  Double I_ERI_F2xz_D2z_Dxy_S_cc = I_ERI_G2x2z_Pz_Dxy_S_cc+ABZ*I_ERI_F2xz_Pz_Dxy_S_cc;
  Double I_ERI_Fx2y_D2z_Dxy_S_cc = I_ERI_Gx2yz_Pz_Dxy_S_cc+ABZ*I_ERI_Fx2y_Pz_Dxy_S_cc;
  Double I_ERI_Fxyz_D2z_Dxy_S_cc = I_ERI_Gxy2z_Pz_Dxy_S_cc+ABZ*I_ERI_Fxyz_Pz_Dxy_S_cc;
  Double I_ERI_Fx2z_D2z_Dxy_S_cc = I_ERI_Gx3z_Pz_Dxy_S_cc+ABZ*I_ERI_Fx2z_Pz_Dxy_S_cc;
  Double I_ERI_F3y_D2z_Dxy_S_cc = I_ERI_G3yz_Pz_Dxy_S_cc+ABZ*I_ERI_F3y_Pz_Dxy_S_cc;
  Double I_ERI_F2yz_D2z_Dxy_S_cc = I_ERI_G2y2z_Pz_Dxy_S_cc+ABZ*I_ERI_F2yz_Pz_Dxy_S_cc;
  Double I_ERI_Fy2z_D2z_Dxy_S_cc = I_ERI_Gy3z_Pz_Dxy_S_cc+ABZ*I_ERI_Fy2z_Pz_Dxy_S_cc;
  Double I_ERI_F3z_D2z_Dxy_S_cc = I_ERI_G4z_Pz_Dxy_S_cc+ABZ*I_ERI_F3z_Pz_Dxy_S_cc;
  Double I_ERI_F3x_D2x_Dxz_S_cc = I_ERI_G4x_Px_Dxz_S_cc+ABX*I_ERI_F3x_Px_Dxz_S_cc;
  Double I_ERI_F2xy_D2x_Dxz_S_cc = I_ERI_G3xy_Px_Dxz_S_cc+ABX*I_ERI_F2xy_Px_Dxz_S_cc;
  Double I_ERI_F2xz_D2x_Dxz_S_cc = I_ERI_G3xz_Px_Dxz_S_cc+ABX*I_ERI_F2xz_Px_Dxz_S_cc;
  Double I_ERI_Fx2y_D2x_Dxz_S_cc = I_ERI_G2x2y_Px_Dxz_S_cc+ABX*I_ERI_Fx2y_Px_Dxz_S_cc;
  Double I_ERI_Fxyz_D2x_Dxz_S_cc = I_ERI_G2xyz_Px_Dxz_S_cc+ABX*I_ERI_Fxyz_Px_Dxz_S_cc;
  Double I_ERI_Fx2z_D2x_Dxz_S_cc = I_ERI_G2x2z_Px_Dxz_S_cc+ABX*I_ERI_Fx2z_Px_Dxz_S_cc;
  Double I_ERI_F3y_D2x_Dxz_S_cc = I_ERI_Gx3y_Px_Dxz_S_cc+ABX*I_ERI_F3y_Px_Dxz_S_cc;
  Double I_ERI_F2yz_D2x_Dxz_S_cc = I_ERI_Gx2yz_Px_Dxz_S_cc+ABX*I_ERI_F2yz_Px_Dxz_S_cc;
  Double I_ERI_Fy2z_D2x_Dxz_S_cc = I_ERI_Gxy2z_Px_Dxz_S_cc+ABX*I_ERI_Fy2z_Px_Dxz_S_cc;
  Double I_ERI_F3z_D2x_Dxz_S_cc = I_ERI_Gx3z_Px_Dxz_S_cc+ABX*I_ERI_F3z_Px_Dxz_S_cc;
  Double I_ERI_F3x_Dxy_Dxz_S_cc = I_ERI_G3xy_Px_Dxz_S_cc+ABY*I_ERI_F3x_Px_Dxz_S_cc;
  Double I_ERI_F2xy_Dxy_Dxz_S_cc = I_ERI_G2x2y_Px_Dxz_S_cc+ABY*I_ERI_F2xy_Px_Dxz_S_cc;
  Double I_ERI_F2xz_Dxy_Dxz_S_cc = I_ERI_G2xyz_Px_Dxz_S_cc+ABY*I_ERI_F2xz_Px_Dxz_S_cc;
  Double I_ERI_Fx2y_Dxy_Dxz_S_cc = I_ERI_Gx3y_Px_Dxz_S_cc+ABY*I_ERI_Fx2y_Px_Dxz_S_cc;
  Double I_ERI_Fxyz_Dxy_Dxz_S_cc = I_ERI_Gx2yz_Px_Dxz_S_cc+ABY*I_ERI_Fxyz_Px_Dxz_S_cc;
  Double I_ERI_Fx2z_Dxy_Dxz_S_cc = I_ERI_Gxy2z_Px_Dxz_S_cc+ABY*I_ERI_Fx2z_Px_Dxz_S_cc;
  Double I_ERI_F3y_Dxy_Dxz_S_cc = I_ERI_G4y_Px_Dxz_S_cc+ABY*I_ERI_F3y_Px_Dxz_S_cc;
  Double I_ERI_F2yz_Dxy_Dxz_S_cc = I_ERI_G3yz_Px_Dxz_S_cc+ABY*I_ERI_F2yz_Px_Dxz_S_cc;
  Double I_ERI_Fy2z_Dxy_Dxz_S_cc = I_ERI_G2y2z_Px_Dxz_S_cc+ABY*I_ERI_Fy2z_Px_Dxz_S_cc;
  Double I_ERI_F3z_Dxy_Dxz_S_cc = I_ERI_Gy3z_Px_Dxz_S_cc+ABY*I_ERI_F3z_Px_Dxz_S_cc;
  Double I_ERI_F3x_Dxz_Dxz_S_cc = I_ERI_G3xz_Px_Dxz_S_cc+ABZ*I_ERI_F3x_Px_Dxz_S_cc;
  Double I_ERI_F2xy_Dxz_Dxz_S_cc = I_ERI_G2xyz_Px_Dxz_S_cc+ABZ*I_ERI_F2xy_Px_Dxz_S_cc;
  Double I_ERI_F2xz_Dxz_Dxz_S_cc = I_ERI_G2x2z_Px_Dxz_S_cc+ABZ*I_ERI_F2xz_Px_Dxz_S_cc;
  Double I_ERI_Fx2y_Dxz_Dxz_S_cc = I_ERI_Gx2yz_Px_Dxz_S_cc+ABZ*I_ERI_Fx2y_Px_Dxz_S_cc;
  Double I_ERI_Fxyz_Dxz_Dxz_S_cc = I_ERI_Gxy2z_Px_Dxz_S_cc+ABZ*I_ERI_Fxyz_Px_Dxz_S_cc;
  Double I_ERI_Fx2z_Dxz_Dxz_S_cc = I_ERI_Gx3z_Px_Dxz_S_cc+ABZ*I_ERI_Fx2z_Px_Dxz_S_cc;
  Double I_ERI_F3y_Dxz_Dxz_S_cc = I_ERI_G3yz_Px_Dxz_S_cc+ABZ*I_ERI_F3y_Px_Dxz_S_cc;
  Double I_ERI_F2yz_Dxz_Dxz_S_cc = I_ERI_G2y2z_Px_Dxz_S_cc+ABZ*I_ERI_F2yz_Px_Dxz_S_cc;
  Double I_ERI_Fy2z_Dxz_Dxz_S_cc = I_ERI_Gy3z_Px_Dxz_S_cc+ABZ*I_ERI_Fy2z_Px_Dxz_S_cc;
  Double I_ERI_F3z_Dxz_Dxz_S_cc = I_ERI_G4z_Px_Dxz_S_cc+ABZ*I_ERI_F3z_Px_Dxz_S_cc;
  Double I_ERI_F3x_D2y_Dxz_S_cc = I_ERI_G3xy_Py_Dxz_S_cc+ABY*I_ERI_F3x_Py_Dxz_S_cc;
  Double I_ERI_F2xy_D2y_Dxz_S_cc = I_ERI_G2x2y_Py_Dxz_S_cc+ABY*I_ERI_F2xy_Py_Dxz_S_cc;
  Double I_ERI_F2xz_D2y_Dxz_S_cc = I_ERI_G2xyz_Py_Dxz_S_cc+ABY*I_ERI_F2xz_Py_Dxz_S_cc;
  Double I_ERI_Fx2y_D2y_Dxz_S_cc = I_ERI_Gx3y_Py_Dxz_S_cc+ABY*I_ERI_Fx2y_Py_Dxz_S_cc;
  Double I_ERI_Fxyz_D2y_Dxz_S_cc = I_ERI_Gx2yz_Py_Dxz_S_cc+ABY*I_ERI_Fxyz_Py_Dxz_S_cc;
  Double I_ERI_Fx2z_D2y_Dxz_S_cc = I_ERI_Gxy2z_Py_Dxz_S_cc+ABY*I_ERI_Fx2z_Py_Dxz_S_cc;
  Double I_ERI_F3y_D2y_Dxz_S_cc = I_ERI_G4y_Py_Dxz_S_cc+ABY*I_ERI_F3y_Py_Dxz_S_cc;
  Double I_ERI_F2yz_D2y_Dxz_S_cc = I_ERI_G3yz_Py_Dxz_S_cc+ABY*I_ERI_F2yz_Py_Dxz_S_cc;
  Double I_ERI_Fy2z_D2y_Dxz_S_cc = I_ERI_G2y2z_Py_Dxz_S_cc+ABY*I_ERI_Fy2z_Py_Dxz_S_cc;
  Double I_ERI_F3z_D2y_Dxz_S_cc = I_ERI_Gy3z_Py_Dxz_S_cc+ABY*I_ERI_F3z_Py_Dxz_S_cc;
  Double I_ERI_F3x_Dyz_Dxz_S_cc = I_ERI_G3xz_Py_Dxz_S_cc+ABZ*I_ERI_F3x_Py_Dxz_S_cc;
  Double I_ERI_F2xy_Dyz_Dxz_S_cc = I_ERI_G2xyz_Py_Dxz_S_cc+ABZ*I_ERI_F2xy_Py_Dxz_S_cc;
  Double I_ERI_F2xz_Dyz_Dxz_S_cc = I_ERI_G2x2z_Py_Dxz_S_cc+ABZ*I_ERI_F2xz_Py_Dxz_S_cc;
  Double I_ERI_Fx2y_Dyz_Dxz_S_cc = I_ERI_Gx2yz_Py_Dxz_S_cc+ABZ*I_ERI_Fx2y_Py_Dxz_S_cc;
  Double I_ERI_Fxyz_Dyz_Dxz_S_cc = I_ERI_Gxy2z_Py_Dxz_S_cc+ABZ*I_ERI_Fxyz_Py_Dxz_S_cc;
  Double I_ERI_Fx2z_Dyz_Dxz_S_cc = I_ERI_Gx3z_Py_Dxz_S_cc+ABZ*I_ERI_Fx2z_Py_Dxz_S_cc;
  Double I_ERI_F3y_Dyz_Dxz_S_cc = I_ERI_G3yz_Py_Dxz_S_cc+ABZ*I_ERI_F3y_Py_Dxz_S_cc;
  Double I_ERI_F2yz_Dyz_Dxz_S_cc = I_ERI_G2y2z_Py_Dxz_S_cc+ABZ*I_ERI_F2yz_Py_Dxz_S_cc;
  Double I_ERI_Fy2z_Dyz_Dxz_S_cc = I_ERI_Gy3z_Py_Dxz_S_cc+ABZ*I_ERI_Fy2z_Py_Dxz_S_cc;
  Double I_ERI_F3z_Dyz_Dxz_S_cc = I_ERI_G4z_Py_Dxz_S_cc+ABZ*I_ERI_F3z_Py_Dxz_S_cc;
  Double I_ERI_F3x_D2z_Dxz_S_cc = I_ERI_G3xz_Pz_Dxz_S_cc+ABZ*I_ERI_F3x_Pz_Dxz_S_cc;
  Double I_ERI_F2xy_D2z_Dxz_S_cc = I_ERI_G2xyz_Pz_Dxz_S_cc+ABZ*I_ERI_F2xy_Pz_Dxz_S_cc;
  Double I_ERI_F2xz_D2z_Dxz_S_cc = I_ERI_G2x2z_Pz_Dxz_S_cc+ABZ*I_ERI_F2xz_Pz_Dxz_S_cc;
  Double I_ERI_Fx2y_D2z_Dxz_S_cc = I_ERI_Gx2yz_Pz_Dxz_S_cc+ABZ*I_ERI_Fx2y_Pz_Dxz_S_cc;
  Double I_ERI_Fxyz_D2z_Dxz_S_cc = I_ERI_Gxy2z_Pz_Dxz_S_cc+ABZ*I_ERI_Fxyz_Pz_Dxz_S_cc;
  Double I_ERI_Fx2z_D2z_Dxz_S_cc = I_ERI_Gx3z_Pz_Dxz_S_cc+ABZ*I_ERI_Fx2z_Pz_Dxz_S_cc;
  Double I_ERI_F3y_D2z_Dxz_S_cc = I_ERI_G3yz_Pz_Dxz_S_cc+ABZ*I_ERI_F3y_Pz_Dxz_S_cc;
  Double I_ERI_F2yz_D2z_Dxz_S_cc = I_ERI_G2y2z_Pz_Dxz_S_cc+ABZ*I_ERI_F2yz_Pz_Dxz_S_cc;
  Double I_ERI_Fy2z_D2z_Dxz_S_cc = I_ERI_Gy3z_Pz_Dxz_S_cc+ABZ*I_ERI_Fy2z_Pz_Dxz_S_cc;
  Double I_ERI_F3z_D2z_Dxz_S_cc = I_ERI_G4z_Pz_Dxz_S_cc+ABZ*I_ERI_F3z_Pz_Dxz_S_cc;
  Double I_ERI_F3x_D2x_D2y_S_cc = I_ERI_G4x_Px_D2y_S_cc+ABX*I_ERI_F3x_Px_D2y_S_cc;
  Double I_ERI_F2xy_D2x_D2y_S_cc = I_ERI_G3xy_Px_D2y_S_cc+ABX*I_ERI_F2xy_Px_D2y_S_cc;
  Double I_ERI_F2xz_D2x_D2y_S_cc = I_ERI_G3xz_Px_D2y_S_cc+ABX*I_ERI_F2xz_Px_D2y_S_cc;
  Double I_ERI_Fx2y_D2x_D2y_S_cc = I_ERI_G2x2y_Px_D2y_S_cc+ABX*I_ERI_Fx2y_Px_D2y_S_cc;
  Double I_ERI_Fxyz_D2x_D2y_S_cc = I_ERI_G2xyz_Px_D2y_S_cc+ABX*I_ERI_Fxyz_Px_D2y_S_cc;
  Double I_ERI_Fx2z_D2x_D2y_S_cc = I_ERI_G2x2z_Px_D2y_S_cc+ABX*I_ERI_Fx2z_Px_D2y_S_cc;
  Double I_ERI_F3y_D2x_D2y_S_cc = I_ERI_Gx3y_Px_D2y_S_cc+ABX*I_ERI_F3y_Px_D2y_S_cc;
  Double I_ERI_F2yz_D2x_D2y_S_cc = I_ERI_Gx2yz_Px_D2y_S_cc+ABX*I_ERI_F2yz_Px_D2y_S_cc;
  Double I_ERI_Fy2z_D2x_D2y_S_cc = I_ERI_Gxy2z_Px_D2y_S_cc+ABX*I_ERI_Fy2z_Px_D2y_S_cc;
  Double I_ERI_F3z_D2x_D2y_S_cc = I_ERI_Gx3z_Px_D2y_S_cc+ABX*I_ERI_F3z_Px_D2y_S_cc;
  Double I_ERI_F3x_Dxy_D2y_S_cc = I_ERI_G3xy_Px_D2y_S_cc+ABY*I_ERI_F3x_Px_D2y_S_cc;
  Double I_ERI_F2xy_Dxy_D2y_S_cc = I_ERI_G2x2y_Px_D2y_S_cc+ABY*I_ERI_F2xy_Px_D2y_S_cc;
  Double I_ERI_F2xz_Dxy_D2y_S_cc = I_ERI_G2xyz_Px_D2y_S_cc+ABY*I_ERI_F2xz_Px_D2y_S_cc;
  Double I_ERI_Fx2y_Dxy_D2y_S_cc = I_ERI_Gx3y_Px_D2y_S_cc+ABY*I_ERI_Fx2y_Px_D2y_S_cc;
  Double I_ERI_Fxyz_Dxy_D2y_S_cc = I_ERI_Gx2yz_Px_D2y_S_cc+ABY*I_ERI_Fxyz_Px_D2y_S_cc;
  Double I_ERI_Fx2z_Dxy_D2y_S_cc = I_ERI_Gxy2z_Px_D2y_S_cc+ABY*I_ERI_Fx2z_Px_D2y_S_cc;
  Double I_ERI_F3y_Dxy_D2y_S_cc = I_ERI_G4y_Px_D2y_S_cc+ABY*I_ERI_F3y_Px_D2y_S_cc;
  Double I_ERI_F2yz_Dxy_D2y_S_cc = I_ERI_G3yz_Px_D2y_S_cc+ABY*I_ERI_F2yz_Px_D2y_S_cc;
  Double I_ERI_Fy2z_Dxy_D2y_S_cc = I_ERI_G2y2z_Px_D2y_S_cc+ABY*I_ERI_Fy2z_Px_D2y_S_cc;
  Double I_ERI_F3z_Dxy_D2y_S_cc = I_ERI_Gy3z_Px_D2y_S_cc+ABY*I_ERI_F3z_Px_D2y_S_cc;
  Double I_ERI_F3x_Dxz_D2y_S_cc = I_ERI_G3xz_Px_D2y_S_cc+ABZ*I_ERI_F3x_Px_D2y_S_cc;
  Double I_ERI_F2xy_Dxz_D2y_S_cc = I_ERI_G2xyz_Px_D2y_S_cc+ABZ*I_ERI_F2xy_Px_D2y_S_cc;
  Double I_ERI_F2xz_Dxz_D2y_S_cc = I_ERI_G2x2z_Px_D2y_S_cc+ABZ*I_ERI_F2xz_Px_D2y_S_cc;
  Double I_ERI_Fx2y_Dxz_D2y_S_cc = I_ERI_Gx2yz_Px_D2y_S_cc+ABZ*I_ERI_Fx2y_Px_D2y_S_cc;
  Double I_ERI_Fxyz_Dxz_D2y_S_cc = I_ERI_Gxy2z_Px_D2y_S_cc+ABZ*I_ERI_Fxyz_Px_D2y_S_cc;
  Double I_ERI_Fx2z_Dxz_D2y_S_cc = I_ERI_Gx3z_Px_D2y_S_cc+ABZ*I_ERI_Fx2z_Px_D2y_S_cc;
  Double I_ERI_F3y_Dxz_D2y_S_cc = I_ERI_G3yz_Px_D2y_S_cc+ABZ*I_ERI_F3y_Px_D2y_S_cc;
  Double I_ERI_F2yz_Dxz_D2y_S_cc = I_ERI_G2y2z_Px_D2y_S_cc+ABZ*I_ERI_F2yz_Px_D2y_S_cc;
  Double I_ERI_Fy2z_Dxz_D2y_S_cc = I_ERI_Gy3z_Px_D2y_S_cc+ABZ*I_ERI_Fy2z_Px_D2y_S_cc;
  Double I_ERI_F3z_Dxz_D2y_S_cc = I_ERI_G4z_Px_D2y_S_cc+ABZ*I_ERI_F3z_Px_D2y_S_cc;
  Double I_ERI_F3x_D2y_D2y_S_cc = I_ERI_G3xy_Py_D2y_S_cc+ABY*I_ERI_F3x_Py_D2y_S_cc;
  Double I_ERI_F2xy_D2y_D2y_S_cc = I_ERI_G2x2y_Py_D2y_S_cc+ABY*I_ERI_F2xy_Py_D2y_S_cc;
  Double I_ERI_F2xz_D2y_D2y_S_cc = I_ERI_G2xyz_Py_D2y_S_cc+ABY*I_ERI_F2xz_Py_D2y_S_cc;
  Double I_ERI_Fx2y_D2y_D2y_S_cc = I_ERI_Gx3y_Py_D2y_S_cc+ABY*I_ERI_Fx2y_Py_D2y_S_cc;
  Double I_ERI_Fxyz_D2y_D2y_S_cc = I_ERI_Gx2yz_Py_D2y_S_cc+ABY*I_ERI_Fxyz_Py_D2y_S_cc;
  Double I_ERI_Fx2z_D2y_D2y_S_cc = I_ERI_Gxy2z_Py_D2y_S_cc+ABY*I_ERI_Fx2z_Py_D2y_S_cc;
  Double I_ERI_F3y_D2y_D2y_S_cc = I_ERI_G4y_Py_D2y_S_cc+ABY*I_ERI_F3y_Py_D2y_S_cc;
  Double I_ERI_F2yz_D2y_D2y_S_cc = I_ERI_G3yz_Py_D2y_S_cc+ABY*I_ERI_F2yz_Py_D2y_S_cc;
  Double I_ERI_Fy2z_D2y_D2y_S_cc = I_ERI_G2y2z_Py_D2y_S_cc+ABY*I_ERI_Fy2z_Py_D2y_S_cc;
  Double I_ERI_F3z_D2y_D2y_S_cc = I_ERI_Gy3z_Py_D2y_S_cc+ABY*I_ERI_F3z_Py_D2y_S_cc;
  Double I_ERI_F3x_Dyz_D2y_S_cc = I_ERI_G3xz_Py_D2y_S_cc+ABZ*I_ERI_F3x_Py_D2y_S_cc;
  Double I_ERI_F2xy_Dyz_D2y_S_cc = I_ERI_G2xyz_Py_D2y_S_cc+ABZ*I_ERI_F2xy_Py_D2y_S_cc;
  Double I_ERI_F2xz_Dyz_D2y_S_cc = I_ERI_G2x2z_Py_D2y_S_cc+ABZ*I_ERI_F2xz_Py_D2y_S_cc;
  Double I_ERI_Fx2y_Dyz_D2y_S_cc = I_ERI_Gx2yz_Py_D2y_S_cc+ABZ*I_ERI_Fx2y_Py_D2y_S_cc;
  Double I_ERI_Fxyz_Dyz_D2y_S_cc = I_ERI_Gxy2z_Py_D2y_S_cc+ABZ*I_ERI_Fxyz_Py_D2y_S_cc;
  Double I_ERI_Fx2z_Dyz_D2y_S_cc = I_ERI_Gx3z_Py_D2y_S_cc+ABZ*I_ERI_Fx2z_Py_D2y_S_cc;
  Double I_ERI_F3y_Dyz_D2y_S_cc = I_ERI_G3yz_Py_D2y_S_cc+ABZ*I_ERI_F3y_Py_D2y_S_cc;
  Double I_ERI_F2yz_Dyz_D2y_S_cc = I_ERI_G2y2z_Py_D2y_S_cc+ABZ*I_ERI_F2yz_Py_D2y_S_cc;
  Double I_ERI_Fy2z_Dyz_D2y_S_cc = I_ERI_Gy3z_Py_D2y_S_cc+ABZ*I_ERI_Fy2z_Py_D2y_S_cc;
  Double I_ERI_F3z_Dyz_D2y_S_cc = I_ERI_G4z_Py_D2y_S_cc+ABZ*I_ERI_F3z_Py_D2y_S_cc;
  Double I_ERI_F3x_D2z_D2y_S_cc = I_ERI_G3xz_Pz_D2y_S_cc+ABZ*I_ERI_F3x_Pz_D2y_S_cc;
  Double I_ERI_F2xy_D2z_D2y_S_cc = I_ERI_G2xyz_Pz_D2y_S_cc+ABZ*I_ERI_F2xy_Pz_D2y_S_cc;
  Double I_ERI_F2xz_D2z_D2y_S_cc = I_ERI_G2x2z_Pz_D2y_S_cc+ABZ*I_ERI_F2xz_Pz_D2y_S_cc;
  Double I_ERI_Fx2y_D2z_D2y_S_cc = I_ERI_Gx2yz_Pz_D2y_S_cc+ABZ*I_ERI_Fx2y_Pz_D2y_S_cc;
  Double I_ERI_Fxyz_D2z_D2y_S_cc = I_ERI_Gxy2z_Pz_D2y_S_cc+ABZ*I_ERI_Fxyz_Pz_D2y_S_cc;
  Double I_ERI_Fx2z_D2z_D2y_S_cc = I_ERI_Gx3z_Pz_D2y_S_cc+ABZ*I_ERI_Fx2z_Pz_D2y_S_cc;
  Double I_ERI_F3y_D2z_D2y_S_cc = I_ERI_G3yz_Pz_D2y_S_cc+ABZ*I_ERI_F3y_Pz_D2y_S_cc;
  Double I_ERI_F2yz_D2z_D2y_S_cc = I_ERI_G2y2z_Pz_D2y_S_cc+ABZ*I_ERI_F2yz_Pz_D2y_S_cc;
  Double I_ERI_Fy2z_D2z_D2y_S_cc = I_ERI_Gy3z_Pz_D2y_S_cc+ABZ*I_ERI_Fy2z_Pz_D2y_S_cc;
  Double I_ERI_F3z_D2z_D2y_S_cc = I_ERI_G4z_Pz_D2y_S_cc+ABZ*I_ERI_F3z_Pz_D2y_S_cc;
  Double I_ERI_F3x_D2x_Dyz_S_cc = I_ERI_G4x_Px_Dyz_S_cc+ABX*I_ERI_F3x_Px_Dyz_S_cc;
  Double I_ERI_F2xy_D2x_Dyz_S_cc = I_ERI_G3xy_Px_Dyz_S_cc+ABX*I_ERI_F2xy_Px_Dyz_S_cc;
  Double I_ERI_F2xz_D2x_Dyz_S_cc = I_ERI_G3xz_Px_Dyz_S_cc+ABX*I_ERI_F2xz_Px_Dyz_S_cc;
  Double I_ERI_Fx2y_D2x_Dyz_S_cc = I_ERI_G2x2y_Px_Dyz_S_cc+ABX*I_ERI_Fx2y_Px_Dyz_S_cc;
  Double I_ERI_Fxyz_D2x_Dyz_S_cc = I_ERI_G2xyz_Px_Dyz_S_cc+ABX*I_ERI_Fxyz_Px_Dyz_S_cc;
  Double I_ERI_Fx2z_D2x_Dyz_S_cc = I_ERI_G2x2z_Px_Dyz_S_cc+ABX*I_ERI_Fx2z_Px_Dyz_S_cc;
  Double I_ERI_F3y_D2x_Dyz_S_cc = I_ERI_Gx3y_Px_Dyz_S_cc+ABX*I_ERI_F3y_Px_Dyz_S_cc;
  Double I_ERI_F2yz_D2x_Dyz_S_cc = I_ERI_Gx2yz_Px_Dyz_S_cc+ABX*I_ERI_F2yz_Px_Dyz_S_cc;
  Double I_ERI_Fy2z_D2x_Dyz_S_cc = I_ERI_Gxy2z_Px_Dyz_S_cc+ABX*I_ERI_Fy2z_Px_Dyz_S_cc;
  Double I_ERI_F3z_D2x_Dyz_S_cc = I_ERI_Gx3z_Px_Dyz_S_cc+ABX*I_ERI_F3z_Px_Dyz_S_cc;
  Double I_ERI_F3x_Dxy_Dyz_S_cc = I_ERI_G3xy_Px_Dyz_S_cc+ABY*I_ERI_F3x_Px_Dyz_S_cc;
  Double I_ERI_F2xy_Dxy_Dyz_S_cc = I_ERI_G2x2y_Px_Dyz_S_cc+ABY*I_ERI_F2xy_Px_Dyz_S_cc;
  Double I_ERI_F2xz_Dxy_Dyz_S_cc = I_ERI_G2xyz_Px_Dyz_S_cc+ABY*I_ERI_F2xz_Px_Dyz_S_cc;
  Double I_ERI_Fx2y_Dxy_Dyz_S_cc = I_ERI_Gx3y_Px_Dyz_S_cc+ABY*I_ERI_Fx2y_Px_Dyz_S_cc;
  Double I_ERI_Fxyz_Dxy_Dyz_S_cc = I_ERI_Gx2yz_Px_Dyz_S_cc+ABY*I_ERI_Fxyz_Px_Dyz_S_cc;
  Double I_ERI_Fx2z_Dxy_Dyz_S_cc = I_ERI_Gxy2z_Px_Dyz_S_cc+ABY*I_ERI_Fx2z_Px_Dyz_S_cc;
  Double I_ERI_F3y_Dxy_Dyz_S_cc = I_ERI_G4y_Px_Dyz_S_cc+ABY*I_ERI_F3y_Px_Dyz_S_cc;
  Double I_ERI_F2yz_Dxy_Dyz_S_cc = I_ERI_G3yz_Px_Dyz_S_cc+ABY*I_ERI_F2yz_Px_Dyz_S_cc;
  Double I_ERI_Fy2z_Dxy_Dyz_S_cc = I_ERI_G2y2z_Px_Dyz_S_cc+ABY*I_ERI_Fy2z_Px_Dyz_S_cc;
  Double I_ERI_F3z_Dxy_Dyz_S_cc = I_ERI_Gy3z_Px_Dyz_S_cc+ABY*I_ERI_F3z_Px_Dyz_S_cc;
  Double I_ERI_F3x_Dxz_Dyz_S_cc = I_ERI_G3xz_Px_Dyz_S_cc+ABZ*I_ERI_F3x_Px_Dyz_S_cc;
  Double I_ERI_F2xy_Dxz_Dyz_S_cc = I_ERI_G2xyz_Px_Dyz_S_cc+ABZ*I_ERI_F2xy_Px_Dyz_S_cc;
  Double I_ERI_F2xz_Dxz_Dyz_S_cc = I_ERI_G2x2z_Px_Dyz_S_cc+ABZ*I_ERI_F2xz_Px_Dyz_S_cc;
  Double I_ERI_Fx2y_Dxz_Dyz_S_cc = I_ERI_Gx2yz_Px_Dyz_S_cc+ABZ*I_ERI_Fx2y_Px_Dyz_S_cc;
  Double I_ERI_Fxyz_Dxz_Dyz_S_cc = I_ERI_Gxy2z_Px_Dyz_S_cc+ABZ*I_ERI_Fxyz_Px_Dyz_S_cc;
  Double I_ERI_Fx2z_Dxz_Dyz_S_cc = I_ERI_Gx3z_Px_Dyz_S_cc+ABZ*I_ERI_Fx2z_Px_Dyz_S_cc;
  Double I_ERI_F3y_Dxz_Dyz_S_cc = I_ERI_G3yz_Px_Dyz_S_cc+ABZ*I_ERI_F3y_Px_Dyz_S_cc;
  Double I_ERI_F2yz_Dxz_Dyz_S_cc = I_ERI_G2y2z_Px_Dyz_S_cc+ABZ*I_ERI_F2yz_Px_Dyz_S_cc;
  Double I_ERI_Fy2z_Dxz_Dyz_S_cc = I_ERI_Gy3z_Px_Dyz_S_cc+ABZ*I_ERI_Fy2z_Px_Dyz_S_cc;
  Double I_ERI_F3z_Dxz_Dyz_S_cc = I_ERI_G4z_Px_Dyz_S_cc+ABZ*I_ERI_F3z_Px_Dyz_S_cc;
  Double I_ERI_F3x_D2y_Dyz_S_cc = I_ERI_G3xy_Py_Dyz_S_cc+ABY*I_ERI_F3x_Py_Dyz_S_cc;
  Double I_ERI_F2xy_D2y_Dyz_S_cc = I_ERI_G2x2y_Py_Dyz_S_cc+ABY*I_ERI_F2xy_Py_Dyz_S_cc;
  Double I_ERI_F2xz_D2y_Dyz_S_cc = I_ERI_G2xyz_Py_Dyz_S_cc+ABY*I_ERI_F2xz_Py_Dyz_S_cc;
  Double I_ERI_Fx2y_D2y_Dyz_S_cc = I_ERI_Gx3y_Py_Dyz_S_cc+ABY*I_ERI_Fx2y_Py_Dyz_S_cc;
  Double I_ERI_Fxyz_D2y_Dyz_S_cc = I_ERI_Gx2yz_Py_Dyz_S_cc+ABY*I_ERI_Fxyz_Py_Dyz_S_cc;
  Double I_ERI_Fx2z_D2y_Dyz_S_cc = I_ERI_Gxy2z_Py_Dyz_S_cc+ABY*I_ERI_Fx2z_Py_Dyz_S_cc;
  Double I_ERI_F3y_D2y_Dyz_S_cc = I_ERI_G4y_Py_Dyz_S_cc+ABY*I_ERI_F3y_Py_Dyz_S_cc;
  Double I_ERI_F2yz_D2y_Dyz_S_cc = I_ERI_G3yz_Py_Dyz_S_cc+ABY*I_ERI_F2yz_Py_Dyz_S_cc;
  Double I_ERI_Fy2z_D2y_Dyz_S_cc = I_ERI_G2y2z_Py_Dyz_S_cc+ABY*I_ERI_Fy2z_Py_Dyz_S_cc;
  Double I_ERI_F3z_D2y_Dyz_S_cc = I_ERI_Gy3z_Py_Dyz_S_cc+ABY*I_ERI_F3z_Py_Dyz_S_cc;
  Double I_ERI_F3x_Dyz_Dyz_S_cc = I_ERI_G3xz_Py_Dyz_S_cc+ABZ*I_ERI_F3x_Py_Dyz_S_cc;
  Double I_ERI_F2xy_Dyz_Dyz_S_cc = I_ERI_G2xyz_Py_Dyz_S_cc+ABZ*I_ERI_F2xy_Py_Dyz_S_cc;
  Double I_ERI_F2xz_Dyz_Dyz_S_cc = I_ERI_G2x2z_Py_Dyz_S_cc+ABZ*I_ERI_F2xz_Py_Dyz_S_cc;
  Double I_ERI_Fx2y_Dyz_Dyz_S_cc = I_ERI_Gx2yz_Py_Dyz_S_cc+ABZ*I_ERI_Fx2y_Py_Dyz_S_cc;
  Double I_ERI_Fxyz_Dyz_Dyz_S_cc = I_ERI_Gxy2z_Py_Dyz_S_cc+ABZ*I_ERI_Fxyz_Py_Dyz_S_cc;
  Double I_ERI_Fx2z_Dyz_Dyz_S_cc = I_ERI_Gx3z_Py_Dyz_S_cc+ABZ*I_ERI_Fx2z_Py_Dyz_S_cc;
  Double I_ERI_F3y_Dyz_Dyz_S_cc = I_ERI_G3yz_Py_Dyz_S_cc+ABZ*I_ERI_F3y_Py_Dyz_S_cc;
  Double I_ERI_F2yz_Dyz_Dyz_S_cc = I_ERI_G2y2z_Py_Dyz_S_cc+ABZ*I_ERI_F2yz_Py_Dyz_S_cc;
  Double I_ERI_Fy2z_Dyz_Dyz_S_cc = I_ERI_Gy3z_Py_Dyz_S_cc+ABZ*I_ERI_Fy2z_Py_Dyz_S_cc;
  Double I_ERI_F3z_Dyz_Dyz_S_cc = I_ERI_G4z_Py_Dyz_S_cc+ABZ*I_ERI_F3z_Py_Dyz_S_cc;
  Double I_ERI_F3x_D2z_Dyz_S_cc = I_ERI_G3xz_Pz_Dyz_S_cc+ABZ*I_ERI_F3x_Pz_Dyz_S_cc;
  Double I_ERI_F2xy_D2z_Dyz_S_cc = I_ERI_G2xyz_Pz_Dyz_S_cc+ABZ*I_ERI_F2xy_Pz_Dyz_S_cc;
  Double I_ERI_F2xz_D2z_Dyz_S_cc = I_ERI_G2x2z_Pz_Dyz_S_cc+ABZ*I_ERI_F2xz_Pz_Dyz_S_cc;
  Double I_ERI_Fx2y_D2z_Dyz_S_cc = I_ERI_Gx2yz_Pz_Dyz_S_cc+ABZ*I_ERI_Fx2y_Pz_Dyz_S_cc;
  Double I_ERI_Fxyz_D2z_Dyz_S_cc = I_ERI_Gxy2z_Pz_Dyz_S_cc+ABZ*I_ERI_Fxyz_Pz_Dyz_S_cc;
  Double I_ERI_Fx2z_D2z_Dyz_S_cc = I_ERI_Gx3z_Pz_Dyz_S_cc+ABZ*I_ERI_Fx2z_Pz_Dyz_S_cc;
  Double I_ERI_F3y_D2z_Dyz_S_cc = I_ERI_G3yz_Pz_Dyz_S_cc+ABZ*I_ERI_F3y_Pz_Dyz_S_cc;
  Double I_ERI_F2yz_D2z_Dyz_S_cc = I_ERI_G2y2z_Pz_Dyz_S_cc+ABZ*I_ERI_F2yz_Pz_Dyz_S_cc;
  Double I_ERI_Fy2z_D2z_Dyz_S_cc = I_ERI_Gy3z_Pz_Dyz_S_cc+ABZ*I_ERI_Fy2z_Pz_Dyz_S_cc;
  Double I_ERI_F3z_D2z_Dyz_S_cc = I_ERI_G4z_Pz_Dyz_S_cc+ABZ*I_ERI_F3z_Pz_Dyz_S_cc;
  Double I_ERI_F3x_D2x_D2z_S_cc = I_ERI_G4x_Px_D2z_S_cc+ABX*I_ERI_F3x_Px_D2z_S_cc;
  Double I_ERI_F2xy_D2x_D2z_S_cc = I_ERI_G3xy_Px_D2z_S_cc+ABX*I_ERI_F2xy_Px_D2z_S_cc;
  Double I_ERI_F2xz_D2x_D2z_S_cc = I_ERI_G3xz_Px_D2z_S_cc+ABX*I_ERI_F2xz_Px_D2z_S_cc;
  Double I_ERI_Fx2y_D2x_D2z_S_cc = I_ERI_G2x2y_Px_D2z_S_cc+ABX*I_ERI_Fx2y_Px_D2z_S_cc;
  Double I_ERI_Fxyz_D2x_D2z_S_cc = I_ERI_G2xyz_Px_D2z_S_cc+ABX*I_ERI_Fxyz_Px_D2z_S_cc;
  Double I_ERI_Fx2z_D2x_D2z_S_cc = I_ERI_G2x2z_Px_D2z_S_cc+ABX*I_ERI_Fx2z_Px_D2z_S_cc;
  Double I_ERI_F3y_D2x_D2z_S_cc = I_ERI_Gx3y_Px_D2z_S_cc+ABX*I_ERI_F3y_Px_D2z_S_cc;
  Double I_ERI_F2yz_D2x_D2z_S_cc = I_ERI_Gx2yz_Px_D2z_S_cc+ABX*I_ERI_F2yz_Px_D2z_S_cc;
  Double I_ERI_Fy2z_D2x_D2z_S_cc = I_ERI_Gxy2z_Px_D2z_S_cc+ABX*I_ERI_Fy2z_Px_D2z_S_cc;
  Double I_ERI_F3z_D2x_D2z_S_cc = I_ERI_Gx3z_Px_D2z_S_cc+ABX*I_ERI_F3z_Px_D2z_S_cc;
  Double I_ERI_F3x_Dxy_D2z_S_cc = I_ERI_G3xy_Px_D2z_S_cc+ABY*I_ERI_F3x_Px_D2z_S_cc;
  Double I_ERI_F2xy_Dxy_D2z_S_cc = I_ERI_G2x2y_Px_D2z_S_cc+ABY*I_ERI_F2xy_Px_D2z_S_cc;
  Double I_ERI_F2xz_Dxy_D2z_S_cc = I_ERI_G2xyz_Px_D2z_S_cc+ABY*I_ERI_F2xz_Px_D2z_S_cc;
  Double I_ERI_Fx2y_Dxy_D2z_S_cc = I_ERI_Gx3y_Px_D2z_S_cc+ABY*I_ERI_Fx2y_Px_D2z_S_cc;
  Double I_ERI_Fxyz_Dxy_D2z_S_cc = I_ERI_Gx2yz_Px_D2z_S_cc+ABY*I_ERI_Fxyz_Px_D2z_S_cc;
  Double I_ERI_Fx2z_Dxy_D2z_S_cc = I_ERI_Gxy2z_Px_D2z_S_cc+ABY*I_ERI_Fx2z_Px_D2z_S_cc;
  Double I_ERI_F3y_Dxy_D2z_S_cc = I_ERI_G4y_Px_D2z_S_cc+ABY*I_ERI_F3y_Px_D2z_S_cc;
  Double I_ERI_F2yz_Dxy_D2z_S_cc = I_ERI_G3yz_Px_D2z_S_cc+ABY*I_ERI_F2yz_Px_D2z_S_cc;
  Double I_ERI_Fy2z_Dxy_D2z_S_cc = I_ERI_G2y2z_Px_D2z_S_cc+ABY*I_ERI_Fy2z_Px_D2z_S_cc;
  Double I_ERI_F3z_Dxy_D2z_S_cc = I_ERI_Gy3z_Px_D2z_S_cc+ABY*I_ERI_F3z_Px_D2z_S_cc;
  Double I_ERI_F3x_Dxz_D2z_S_cc = I_ERI_G3xz_Px_D2z_S_cc+ABZ*I_ERI_F3x_Px_D2z_S_cc;
  Double I_ERI_F2xy_Dxz_D2z_S_cc = I_ERI_G2xyz_Px_D2z_S_cc+ABZ*I_ERI_F2xy_Px_D2z_S_cc;
  Double I_ERI_F2xz_Dxz_D2z_S_cc = I_ERI_G2x2z_Px_D2z_S_cc+ABZ*I_ERI_F2xz_Px_D2z_S_cc;
  Double I_ERI_Fx2y_Dxz_D2z_S_cc = I_ERI_Gx2yz_Px_D2z_S_cc+ABZ*I_ERI_Fx2y_Px_D2z_S_cc;
  Double I_ERI_Fxyz_Dxz_D2z_S_cc = I_ERI_Gxy2z_Px_D2z_S_cc+ABZ*I_ERI_Fxyz_Px_D2z_S_cc;
  Double I_ERI_Fx2z_Dxz_D2z_S_cc = I_ERI_Gx3z_Px_D2z_S_cc+ABZ*I_ERI_Fx2z_Px_D2z_S_cc;
  Double I_ERI_F3y_Dxz_D2z_S_cc = I_ERI_G3yz_Px_D2z_S_cc+ABZ*I_ERI_F3y_Px_D2z_S_cc;
  Double I_ERI_F2yz_Dxz_D2z_S_cc = I_ERI_G2y2z_Px_D2z_S_cc+ABZ*I_ERI_F2yz_Px_D2z_S_cc;
  Double I_ERI_Fy2z_Dxz_D2z_S_cc = I_ERI_Gy3z_Px_D2z_S_cc+ABZ*I_ERI_Fy2z_Px_D2z_S_cc;
  Double I_ERI_F3z_Dxz_D2z_S_cc = I_ERI_G4z_Px_D2z_S_cc+ABZ*I_ERI_F3z_Px_D2z_S_cc;
  Double I_ERI_F3x_D2y_D2z_S_cc = I_ERI_G3xy_Py_D2z_S_cc+ABY*I_ERI_F3x_Py_D2z_S_cc;
  Double I_ERI_F2xy_D2y_D2z_S_cc = I_ERI_G2x2y_Py_D2z_S_cc+ABY*I_ERI_F2xy_Py_D2z_S_cc;
  Double I_ERI_F2xz_D2y_D2z_S_cc = I_ERI_G2xyz_Py_D2z_S_cc+ABY*I_ERI_F2xz_Py_D2z_S_cc;
  Double I_ERI_Fx2y_D2y_D2z_S_cc = I_ERI_Gx3y_Py_D2z_S_cc+ABY*I_ERI_Fx2y_Py_D2z_S_cc;
  Double I_ERI_Fxyz_D2y_D2z_S_cc = I_ERI_Gx2yz_Py_D2z_S_cc+ABY*I_ERI_Fxyz_Py_D2z_S_cc;
  Double I_ERI_Fx2z_D2y_D2z_S_cc = I_ERI_Gxy2z_Py_D2z_S_cc+ABY*I_ERI_Fx2z_Py_D2z_S_cc;
  Double I_ERI_F3y_D2y_D2z_S_cc = I_ERI_G4y_Py_D2z_S_cc+ABY*I_ERI_F3y_Py_D2z_S_cc;
  Double I_ERI_F2yz_D2y_D2z_S_cc = I_ERI_G3yz_Py_D2z_S_cc+ABY*I_ERI_F2yz_Py_D2z_S_cc;
  Double I_ERI_Fy2z_D2y_D2z_S_cc = I_ERI_G2y2z_Py_D2z_S_cc+ABY*I_ERI_Fy2z_Py_D2z_S_cc;
  Double I_ERI_F3z_D2y_D2z_S_cc = I_ERI_Gy3z_Py_D2z_S_cc+ABY*I_ERI_F3z_Py_D2z_S_cc;
  Double I_ERI_F3x_Dyz_D2z_S_cc = I_ERI_G3xz_Py_D2z_S_cc+ABZ*I_ERI_F3x_Py_D2z_S_cc;
  Double I_ERI_F2xy_Dyz_D2z_S_cc = I_ERI_G2xyz_Py_D2z_S_cc+ABZ*I_ERI_F2xy_Py_D2z_S_cc;
  Double I_ERI_F2xz_Dyz_D2z_S_cc = I_ERI_G2x2z_Py_D2z_S_cc+ABZ*I_ERI_F2xz_Py_D2z_S_cc;
  Double I_ERI_Fx2y_Dyz_D2z_S_cc = I_ERI_Gx2yz_Py_D2z_S_cc+ABZ*I_ERI_Fx2y_Py_D2z_S_cc;
  Double I_ERI_Fxyz_Dyz_D2z_S_cc = I_ERI_Gxy2z_Py_D2z_S_cc+ABZ*I_ERI_Fxyz_Py_D2z_S_cc;
  Double I_ERI_Fx2z_Dyz_D2z_S_cc = I_ERI_Gx3z_Py_D2z_S_cc+ABZ*I_ERI_Fx2z_Py_D2z_S_cc;
  Double I_ERI_F3y_Dyz_D2z_S_cc = I_ERI_G3yz_Py_D2z_S_cc+ABZ*I_ERI_F3y_Py_D2z_S_cc;
  Double I_ERI_F2yz_Dyz_D2z_S_cc = I_ERI_G2y2z_Py_D2z_S_cc+ABZ*I_ERI_F2yz_Py_D2z_S_cc;
  Double I_ERI_Fy2z_Dyz_D2z_S_cc = I_ERI_Gy3z_Py_D2z_S_cc+ABZ*I_ERI_Fy2z_Py_D2z_S_cc;
  Double I_ERI_F3z_Dyz_D2z_S_cc = I_ERI_G4z_Py_D2z_S_cc+ABZ*I_ERI_F3z_Py_D2z_S_cc;
  Double I_ERI_F3x_D2z_D2z_S_cc = I_ERI_G3xz_Pz_D2z_S_cc+ABZ*I_ERI_F3x_Pz_D2z_S_cc;
  Double I_ERI_F2xy_D2z_D2z_S_cc = I_ERI_G2xyz_Pz_D2z_S_cc+ABZ*I_ERI_F2xy_Pz_D2z_S_cc;
  Double I_ERI_F2xz_D2z_D2z_S_cc = I_ERI_G2x2z_Pz_D2z_S_cc+ABZ*I_ERI_F2xz_Pz_D2z_S_cc;
  Double I_ERI_Fx2y_D2z_D2z_S_cc = I_ERI_Gx2yz_Pz_D2z_S_cc+ABZ*I_ERI_Fx2y_Pz_D2z_S_cc;
  Double I_ERI_Fxyz_D2z_D2z_S_cc = I_ERI_Gxy2z_Pz_D2z_S_cc+ABZ*I_ERI_Fxyz_Pz_D2z_S_cc;
  Double I_ERI_Fx2z_D2z_D2z_S_cc = I_ERI_Gx3z_Pz_D2z_S_cc+ABZ*I_ERI_Fx2z_Pz_D2z_S_cc;
  Double I_ERI_F3y_D2z_D2z_S_cc = I_ERI_G3yz_Pz_D2z_S_cc+ABZ*I_ERI_F3y_Pz_D2z_S_cc;
  Double I_ERI_F2yz_D2z_D2z_S_cc = I_ERI_G2y2z_Pz_D2z_S_cc+ABZ*I_ERI_F2yz_Pz_D2z_S_cc;
  Double I_ERI_Fy2z_D2z_D2z_S_cc = I_ERI_Gy3z_Pz_D2z_S_cc+ABZ*I_ERI_Fy2z_Pz_D2z_S_cc;
  Double I_ERI_F3z_D2z_D2z_S_cc = I_ERI_G4z_Pz_D2z_S_cc+ABZ*I_ERI_F3z_Pz_D2z_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_S_S_aa
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S
   ************************************************************/
  abcd[0] = 4.0E0*I_ERI_H5x_D2x_S_S_aa-2.0E0*3*I_ERI_F3x_D2x_S_S_a-2.0E0*4*I_ERI_F3x_D2x_S_S_a+3*2*I_ERI_Px_D2x_S_S;
  abcd[1] = 4.0E0*I_ERI_H4xy_D2x_S_S_aa-2.0E0*2*I_ERI_F2xy_D2x_S_S_a-2.0E0*3*I_ERI_F2xy_D2x_S_S_a+2*1*I_ERI_Py_D2x_S_S;
  abcd[2] = 4.0E0*I_ERI_H4xz_D2x_S_S_aa-2.0E0*2*I_ERI_F2xz_D2x_S_S_a-2.0E0*3*I_ERI_F2xz_D2x_S_S_a+2*1*I_ERI_Pz_D2x_S_S;
  abcd[3] = 4.0E0*I_ERI_H3x2y_D2x_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2x_S_S_a-2.0E0*2*I_ERI_Fx2y_D2x_S_S_a;
  abcd[4] = 4.0E0*I_ERI_H3xyz_D2x_S_S_aa-2.0E0*1*I_ERI_Fxyz_D2x_S_S_a-2.0E0*2*I_ERI_Fxyz_D2x_S_S_a;
  abcd[5] = 4.0E0*I_ERI_H3x2z_D2x_S_S_aa-2.0E0*1*I_ERI_Fx2z_D2x_S_S_a-2.0E0*2*I_ERI_Fx2z_D2x_S_S_a;
  abcd[6] = 4.0E0*I_ERI_H2x3y_D2x_S_S_aa-2.0E0*1*I_ERI_F3y_D2x_S_S_a;
  abcd[7] = 4.0E0*I_ERI_H2x2yz_D2x_S_S_aa-2.0E0*1*I_ERI_F2yz_D2x_S_S_a;
  abcd[8] = 4.0E0*I_ERI_H2xy2z_D2x_S_S_aa-2.0E0*1*I_ERI_Fy2z_D2x_S_S_a;
  abcd[9] = 4.0E0*I_ERI_H2x3z_D2x_S_S_aa-2.0E0*1*I_ERI_F3z_D2x_S_S_a;
  abcd[10] = 4.0E0*I_ERI_H5x_Dxy_S_S_aa-2.0E0*3*I_ERI_F3x_Dxy_S_S_a-2.0E0*4*I_ERI_F3x_Dxy_S_S_a+3*2*I_ERI_Px_Dxy_S_S;
  abcd[11] = 4.0E0*I_ERI_H4xy_Dxy_S_S_aa-2.0E0*2*I_ERI_F2xy_Dxy_S_S_a-2.0E0*3*I_ERI_F2xy_Dxy_S_S_a+2*1*I_ERI_Py_Dxy_S_S;
  abcd[12] = 4.0E0*I_ERI_H4xz_Dxy_S_S_aa-2.0E0*2*I_ERI_F2xz_Dxy_S_S_a-2.0E0*3*I_ERI_F2xz_Dxy_S_S_a+2*1*I_ERI_Pz_Dxy_S_S;
  abcd[13] = 4.0E0*I_ERI_H3x2y_Dxy_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_a-2.0E0*2*I_ERI_Fx2y_Dxy_S_S_a;
  abcd[14] = 4.0E0*I_ERI_H3xyz_Dxy_S_S_aa-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_a-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_a;
  abcd[15] = 4.0E0*I_ERI_H3x2z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_a-2.0E0*2*I_ERI_Fx2z_Dxy_S_S_a;
  abcd[16] = 4.0E0*I_ERI_H2x3y_Dxy_S_S_aa-2.0E0*1*I_ERI_F3y_Dxy_S_S_a;
  abcd[17] = 4.0E0*I_ERI_H2x2yz_Dxy_S_S_aa-2.0E0*1*I_ERI_F2yz_Dxy_S_S_a;
  abcd[18] = 4.0E0*I_ERI_H2xy2z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_a;
  abcd[19] = 4.0E0*I_ERI_H2x3z_Dxy_S_S_aa-2.0E0*1*I_ERI_F3z_Dxy_S_S_a;
  abcd[20] = 4.0E0*I_ERI_H5x_Dxz_S_S_aa-2.0E0*3*I_ERI_F3x_Dxz_S_S_a-2.0E0*4*I_ERI_F3x_Dxz_S_S_a+3*2*I_ERI_Px_Dxz_S_S;
  abcd[21] = 4.0E0*I_ERI_H4xy_Dxz_S_S_aa-2.0E0*2*I_ERI_F2xy_Dxz_S_S_a-2.0E0*3*I_ERI_F2xy_Dxz_S_S_a+2*1*I_ERI_Py_Dxz_S_S;
  abcd[22] = 4.0E0*I_ERI_H4xz_Dxz_S_S_aa-2.0E0*2*I_ERI_F2xz_Dxz_S_S_a-2.0E0*3*I_ERI_F2xz_Dxz_S_S_a+2*1*I_ERI_Pz_Dxz_S_S;
  abcd[23] = 4.0E0*I_ERI_H3x2y_Dxz_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_a-2.0E0*2*I_ERI_Fx2y_Dxz_S_S_a;
  abcd[24] = 4.0E0*I_ERI_H3xyz_Dxz_S_S_aa-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_a-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_a;
  abcd[25] = 4.0E0*I_ERI_H3x2z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_a-2.0E0*2*I_ERI_Fx2z_Dxz_S_S_a;
  abcd[26] = 4.0E0*I_ERI_H2x3y_Dxz_S_S_aa-2.0E0*1*I_ERI_F3y_Dxz_S_S_a;
  abcd[27] = 4.0E0*I_ERI_H2x2yz_Dxz_S_S_aa-2.0E0*1*I_ERI_F2yz_Dxz_S_S_a;
  abcd[28] = 4.0E0*I_ERI_H2xy2z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_a;
  abcd[29] = 4.0E0*I_ERI_H2x3z_Dxz_S_S_aa-2.0E0*1*I_ERI_F3z_Dxz_S_S_a;
  abcd[30] = 4.0E0*I_ERI_H5x_D2y_S_S_aa-2.0E0*3*I_ERI_F3x_D2y_S_S_a-2.0E0*4*I_ERI_F3x_D2y_S_S_a+3*2*I_ERI_Px_D2y_S_S;
  abcd[31] = 4.0E0*I_ERI_H4xy_D2y_S_S_aa-2.0E0*2*I_ERI_F2xy_D2y_S_S_a-2.0E0*3*I_ERI_F2xy_D2y_S_S_a+2*1*I_ERI_Py_D2y_S_S;
  abcd[32] = 4.0E0*I_ERI_H4xz_D2y_S_S_aa-2.0E0*2*I_ERI_F2xz_D2y_S_S_a-2.0E0*3*I_ERI_F2xz_D2y_S_S_a+2*1*I_ERI_Pz_D2y_S_S;
  abcd[33] = 4.0E0*I_ERI_H3x2y_D2y_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2y_S_S_a-2.0E0*2*I_ERI_Fx2y_D2y_S_S_a;
  abcd[34] = 4.0E0*I_ERI_H3xyz_D2y_S_S_aa-2.0E0*1*I_ERI_Fxyz_D2y_S_S_a-2.0E0*2*I_ERI_Fxyz_D2y_S_S_a;
  abcd[35] = 4.0E0*I_ERI_H3x2z_D2y_S_S_aa-2.0E0*1*I_ERI_Fx2z_D2y_S_S_a-2.0E0*2*I_ERI_Fx2z_D2y_S_S_a;
  abcd[36] = 4.0E0*I_ERI_H2x3y_D2y_S_S_aa-2.0E0*1*I_ERI_F3y_D2y_S_S_a;
  abcd[37] = 4.0E0*I_ERI_H2x2yz_D2y_S_S_aa-2.0E0*1*I_ERI_F2yz_D2y_S_S_a;
  abcd[38] = 4.0E0*I_ERI_H2xy2z_D2y_S_S_aa-2.0E0*1*I_ERI_Fy2z_D2y_S_S_a;
  abcd[39] = 4.0E0*I_ERI_H2x3z_D2y_S_S_aa-2.0E0*1*I_ERI_F3z_D2y_S_S_a;
  abcd[40] = 4.0E0*I_ERI_H5x_Dyz_S_S_aa-2.0E0*3*I_ERI_F3x_Dyz_S_S_a-2.0E0*4*I_ERI_F3x_Dyz_S_S_a+3*2*I_ERI_Px_Dyz_S_S;
  abcd[41] = 4.0E0*I_ERI_H4xy_Dyz_S_S_aa-2.0E0*2*I_ERI_F2xy_Dyz_S_S_a-2.0E0*3*I_ERI_F2xy_Dyz_S_S_a+2*1*I_ERI_Py_Dyz_S_S;
  abcd[42] = 4.0E0*I_ERI_H4xz_Dyz_S_S_aa-2.0E0*2*I_ERI_F2xz_Dyz_S_S_a-2.0E0*3*I_ERI_F2xz_Dyz_S_S_a+2*1*I_ERI_Pz_Dyz_S_S;
  abcd[43] = 4.0E0*I_ERI_H3x2y_Dyz_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_a-2.0E0*2*I_ERI_Fx2y_Dyz_S_S_a;
  abcd[44] = 4.0E0*I_ERI_H3xyz_Dyz_S_S_aa-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_a-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_a;
  abcd[45] = 4.0E0*I_ERI_H3x2z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_a-2.0E0*2*I_ERI_Fx2z_Dyz_S_S_a;
  abcd[46] = 4.0E0*I_ERI_H2x3y_Dyz_S_S_aa-2.0E0*1*I_ERI_F3y_Dyz_S_S_a;
  abcd[47] = 4.0E0*I_ERI_H2x2yz_Dyz_S_S_aa-2.0E0*1*I_ERI_F2yz_Dyz_S_S_a;
  abcd[48] = 4.0E0*I_ERI_H2xy2z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_a;
  abcd[49] = 4.0E0*I_ERI_H2x3z_Dyz_S_S_aa-2.0E0*1*I_ERI_F3z_Dyz_S_S_a;
  abcd[50] = 4.0E0*I_ERI_H5x_D2z_S_S_aa-2.0E0*3*I_ERI_F3x_D2z_S_S_a-2.0E0*4*I_ERI_F3x_D2z_S_S_a+3*2*I_ERI_Px_D2z_S_S;
  abcd[51] = 4.0E0*I_ERI_H4xy_D2z_S_S_aa-2.0E0*2*I_ERI_F2xy_D2z_S_S_a-2.0E0*3*I_ERI_F2xy_D2z_S_S_a+2*1*I_ERI_Py_D2z_S_S;
  abcd[52] = 4.0E0*I_ERI_H4xz_D2z_S_S_aa-2.0E0*2*I_ERI_F2xz_D2z_S_S_a-2.0E0*3*I_ERI_F2xz_D2z_S_S_a+2*1*I_ERI_Pz_D2z_S_S;
  abcd[53] = 4.0E0*I_ERI_H3x2y_D2z_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2z_S_S_a-2.0E0*2*I_ERI_Fx2y_D2z_S_S_a;
  abcd[54] = 4.0E0*I_ERI_H3xyz_D2z_S_S_aa-2.0E0*1*I_ERI_Fxyz_D2z_S_S_a-2.0E0*2*I_ERI_Fxyz_D2z_S_S_a;
  abcd[55] = 4.0E0*I_ERI_H3x2z_D2z_S_S_aa-2.0E0*1*I_ERI_Fx2z_D2z_S_S_a-2.0E0*2*I_ERI_Fx2z_D2z_S_S_a;
  abcd[56] = 4.0E0*I_ERI_H2x3y_D2z_S_S_aa-2.0E0*1*I_ERI_F3y_D2z_S_S_a;
  abcd[57] = 4.0E0*I_ERI_H2x2yz_D2z_S_S_aa-2.0E0*1*I_ERI_F2yz_D2z_S_S_a;
  abcd[58] = 4.0E0*I_ERI_H2xy2z_D2z_S_S_aa-2.0E0*1*I_ERI_Fy2z_D2z_S_S_a;
  abcd[59] = 4.0E0*I_ERI_H2x3z_D2z_S_S_aa-2.0E0*1*I_ERI_F3z_D2z_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_S_S_aa
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S
   ************************************************************/
  abcd[60] = 4.0E0*I_ERI_H4xy_D2x_S_S_aa-2.0E0*3*I_ERI_F2xy_D2x_S_S_a;
  abcd[61] = 4.0E0*I_ERI_H3x2y_D2x_S_S_aa-2.0E0*1*I_ERI_F3x_D2x_S_S_a-2.0E0*2*I_ERI_Fx2y_D2x_S_S_a+2*1*I_ERI_Px_D2x_S_S;
  abcd[62] = 4.0E0*I_ERI_H3xyz_D2x_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2x_S_S_a;
  abcd[63] = 4.0E0*I_ERI_H2x3y_D2x_S_S_aa-2.0E0*2*I_ERI_F2xy_D2x_S_S_a-2.0E0*1*I_ERI_F3y_D2x_S_S_a+2*I_ERI_Py_D2x_S_S;
  abcd[64] = 4.0E0*I_ERI_H2x2yz_D2x_S_S_aa-2.0E0*1*I_ERI_F2xz_D2x_S_S_a-2.0E0*1*I_ERI_F2yz_D2x_S_S_a+1*I_ERI_Pz_D2x_S_S;
  abcd[65] = 4.0E0*I_ERI_H2xy2z_D2x_S_S_aa-2.0E0*1*I_ERI_Fy2z_D2x_S_S_a;
  abcd[66] = 4.0E0*I_ERI_Hx4y_D2x_S_S_aa-2.0E0*3*I_ERI_Fx2y_D2x_S_S_a;
  abcd[67] = 4.0E0*I_ERI_Hx3yz_D2x_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2x_S_S_a;
  abcd[68] = 4.0E0*I_ERI_Hx2y2z_D2x_S_S_aa-2.0E0*1*I_ERI_Fx2z_D2x_S_S_a;
  abcd[69] = 4.0E0*I_ERI_Hxy3z_D2x_S_S_aa;
  abcd[70] = 4.0E0*I_ERI_H4xy_Dxy_S_S_aa-2.0E0*3*I_ERI_F2xy_Dxy_S_S_a;
  abcd[71] = 4.0E0*I_ERI_H3x2y_Dxy_S_S_aa-2.0E0*1*I_ERI_F3x_Dxy_S_S_a-2.0E0*2*I_ERI_Fx2y_Dxy_S_S_a+2*1*I_ERI_Px_Dxy_S_S;
  abcd[72] = 4.0E0*I_ERI_H3xyz_Dxy_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_a;
  abcd[73] = 4.0E0*I_ERI_H2x3y_Dxy_S_S_aa-2.0E0*2*I_ERI_F2xy_Dxy_S_S_a-2.0E0*1*I_ERI_F3y_Dxy_S_S_a+2*I_ERI_Py_Dxy_S_S;
  abcd[74] = 4.0E0*I_ERI_H2x2yz_Dxy_S_S_aa-2.0E0*1*I_ERI_F2xz_Dxy_S_S_a-2.0E0*1*I_ERI_F2yz_Dxy_S_S_a+1*I_ERI_Pz_Dxy_S_S;
  abcd[75] = 4.0E0*I_ERI_H2xy2z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_a;
  abcd[76] = 4.0E0*I_ERI_Hx4y_Dxy_S_S_aa-2.0E0*3*I_ERI_Fx2y_Dxy_S_S_a;
  abcd[77] = 4.0E0*I_ERI_Hx3yz_Dxy_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_a;
  abcd[78] = 4.0E0*I_ERI_Hx2y2z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_a;
  abcd[79] = 4.0E0*I_ERI_Hxy3z_Dxy_S_S_aa;
  abcd[80] = 4.0E0*I_ERI_H4xy_Dxz_S_S_aa-2.0E0*3*I_ERI_F2xy_Dxz_S_S_a;
  abcd[81] = 4.0E0*I_ERI_H3x2y_Dxz_S_S_aa-2.0E0*1*I_ERI_F3x_Dxz_S_S_a-2.0E0*2*I_ERI_Fx2y_Dxz_S_S_a+2*1*I_ERI_Px_Dxz_S_S;
  abcd[82] = 4.0E0*I_ERI_H3xyz_Dxz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_a;
  abcd[83] = 4.0E0*I_ERI_H2x3y_Dxz_S_S_aa-2.0E0*2*I_ERI_F2xy_Dxz_S_S_a-2.0E0*1*I_ERI_F3y_Dxz_S_S_a+2*I_ERI_Py_Dxz_S_S;
  abcd[84] = 4.0E0*I_ERI_H2x2yz_Dxz_S_S_aa-2.0E0*1*I_ERI_F2xz_Dxz_S_S_a-2.0E0*1*I_ERI_F2yz_Dxz_S_S_a+1*I_ERI_Pz_Dxz_S_S;
  abcd[85] = 4.0E0*I_ERI_H2xy2z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_a;
  abcd[86] = 4.0E0*I_ERI_Hx4y_Dxz_S_S_aa-2.0E0*3*I_ERI_Fx2y_Dxz_S_S_a;
  abcd[87] = 4.0E0*I_ERI_Hx3yz_Dxz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_a;
  abcd[88] = 4.0E0*I_ERI_Hx2y2z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_a;
  abcd[89] = 4.0E0*I_ERI_Hxy3z_Dxz_S_S_aa;
  abcd[90] = 4.0E0*I_ERI_H4xy_D2y_S_S_aa-2.0E0*3*I_ERI_F2xy_D2y_S_S_a;
  abcd[91] = 4.0E0*I_ERI_H3x2y_D2y_S_S_aa-2.0E0*1*I_ERI_F3x_D2y_S_S_a-2.0E0*2*I_ERI_Fx2y_D2y_S_S_a+2*1*I_ERI_Px_D2y_S_S;
  abcd[92] = 4.0E0*I_ERI_H3xyz_D2y_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2y_S_S_a;
  abcd[93] = 4.0E0*I_ERI_H2x3y_D2y_S_S_aa-2.0E0*2*I_ERI_F2xy_D2y_S_S_a-2.0E0*1*I_ERI_F3y_D2y_S_S_a+2*I_ERI_Py_D2y_S_S;
  abcd[94] = 4.0E0*I_ERI_H2x2yz_D2y_S_S_aa-2.0E0*1*I_ERI_F2xz_D2y_S_S_a-2.0E0*1*I_ERI_F2yz_D2y_S_S_a+1*I_ERI_Pz_D2y_S_S;
  abcd[95] = 4.0E0*I_ERI_H2xy2z_D2y_S_S_aa-2.0E0*1*I_ERI_Fy2z_D2y_S_S_a;
  abcd[96] = 4.0E0*I_ERI_Hx4y_D2y_S_S_aa-2.0E0*3*I_ERI_Fx2y_D2y_S_S_a;
  abcd[97] = 4.0E0*I_ERI_Hx3yz_D2y_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2y_S_S_a;
  abcd[98] = 4.0E0*I_ERI_Hx2y2z_D2y_S_S_aa-2.0E0*1*I_ERI_Fx2z_D2y_S_S_a;
  abcd[99] = 4.0E0*I_ERI_Hxy3z_D2y_S_S_aa;
  abcd[100] = 4.0E0*I_ERI_H4xy_Dyz_S_S_aa-2.0E0*3*I_ERI_F2xy_Dyz_S_S_a;
  abcd[101] = 4.0E0*I_ERI_H3x2y_Dyz_S_S_aa-2.0E0*1*I_ERI_F3x_Dyz_S_S_a-2.0E0*2*I_ERI_Fx2y_Dyz_S_S_a+2*1*I_ERI_Px_Dyz_S_S;
  abcd[102] = 4.0E0*I_ERI_H3xyz_Dyz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_a;
  abcd[103] = 4.0E0*I_ERI_H2x3y_Dyz_S_S_aa-2.0E0*2*I_ERI_F2xy_Dyz_S_S_a-2.0E0*1*I_ERI_F3y_Dyz_S_S_a+2*I_ERI_Py_Dyz_S_S;
  abcd[104] = 4.0E0*I_ERI_H2x2yz_Dyz_S_S_aa-2.0E0*1*I_ERI_F2xz_Dyz_S_S_a-2.0E0*1*I_ERI_F2yz_Dyz_S_S_a+1*I_ERI_Pz_Dyz_S_S;
  abcd[105] = 4.0E0*I_ERI_H2xy2z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_a;
  abcd[106] = 4.0E0*I_ERI_Hx4y_Dyz_S_S_aa-2.0E0*3*I_ERI_Fx2y_Dyz_S_S_a;
  abcd[107] = 4.0E0*I_ERI_Hx3yz_Dyz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_a;
  abcd[108] = 4.0E0*I_ERI_Hx2y2z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_a;
  abcd[109] = 4.0E0*I_ERI_Hxy3z_Dyz_S_S_aa;
  abcd[110] = 4.0E0*I_ERI_H4xy_D2z_S_S_aa-2.0E0*3*I_ERI_F2xy_D2z_S_S_a;
  abcd[111] = 4.0E0*I_ERI_H3x2y_D2z_S_S_aa-2.0E0*1*I_ERI_F3x_D2z_S_S_a-2.0E0*2*I_ERI_Fx2y_D2z_S_S_a+2*1*I_ERI_Px_D2z_S_S;
  abcd[112] = 4.0E0*I_ERI_H3xyz_D2z_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2z_S_S_a;
  abcd[113] = 4.0E0*I_ERI_H2x3y_D2z_S_S_aa-2.0E0*2*I_ERI_F2xy_D2z_S_S_a-2.0E0*1*I_ERI_F3y_D2z_S_S_a+2*I_ERI_Py_D2z_S_S;
  abcd[114] = 4.0E0*I_ERI_H2x2yz_D2z_S_S_aa-2.0E0*1*I_ERI_F2xz_D2z_S_S_a-2.0E0*1*I_ERI_F2yz_D2z_S_S_a+1*I_ERI_Pz_D2z_S_S;
  abcd[115] = 4.0E0*I_ERI_H2xy2z_D2z_S_S_aa-2.0E0*1*I_ERI_Fy2z_D2z_S_S_a;
  abcd[116] = 4.0E0*I_ERI_Hx4y_D2z_S_S_aa-2.0E0*3*I_ERI_Fx2y_D2z_S_S_a;
  abcd[117] = 4.0E0*I_ERI_Hx3yz_D2z_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2z_S_S_a;
  abcd[118] = 4.0E0*I_ERI_Hx2y2z_D2z_S_S_aa-2.0E0*1*I_ERI_Fx2z_D2z_S_S_a;
  abcd[119] = 4.0E0*I_ERI_Hxy3z_D2z_S_S_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_S_S_aa
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S
   ************************************************************/
  abcd[120] = 4.0E0*I_ERI_H4xz_D2x_S_S_aa-2.0E0*3*I_ERI_F2xz_D2x_S_S_a;
  abcd[121] = 4.0E0*I_ERI_H3xyz_D2x_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2x_S_S_a;
  abcd[122] = 4.0E0*I_ERI_H3x2z_D2x_S_S_aa-2.0E0*1*I_ERI_F3x_D2x_S_S_a-2.0E0*2*I_ERI_Fx2z_D2x_S_S_a+2*1*I_ERI_Px_D2x_S_S;
  abcd[123] = 4.0E0*I_ERI_H2x2yz_D2x_S_S_aa-2.0E0*1*I_ERI_F2yz_D2x_S_S_a;
  abcd[124] = 4.0E0*I_ERI_H2xy2z_D2x_S_S_aa-2.0E0*1*I_ERI_F2xy_D2x_S_S_a-2.0E0*1*I_ERI_Fy2z_D2x_S_S_a+1*I_ERI_Py_D2x_S_S;
  abcd[125] = 4.0E0*I_ERI_H2x3z_D2x_S_S_aa-2.0E0*2*I_ERI_F2xz_D2x_S_S_a-2.0E0*1*I_ERI_F3z_D2x_S_S_a+2*I_ERI_Pz_D2x_S_S;
  abcd[126] = 4.0E0*I_ERI_Hx3yz_D2x_S_S_aa;
  abcd[127] = 4.0E0*I_ERI_Hx2y2z_D2x_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2x_S_S_a;
  abcd[128] = 4.0E0*I_ERI_Hxy3z_D2x_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2x_S_S_a;
  abcd[129] = 4.0E0*I_ERI_Hx4z_D2x_S_S_aa-2.0E0*3*I_ERI_Fx2z_D2x_S_S_a;
  abcd[130] = 4.0E0*I_ERI_H4xz_Dxy_S_S_aa-2.0E0*3*I_ERI_F2xz_Dxy_S_S_a;
  abcd[131] = 4.0E0*I_ERI_H3xyz_Dxy_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_a;
  abcd[132] = 4.0E0*I_ERI_H3x2z_Dxy_S_S_aa-2.0E0*1*I_ERI_F3x_Dxy_S_S_a-2.0E0*2*I_ERI_Fx2z_Dxy_S_S_a+2*1*I_ERI_Px_Dxy_S_S;
  abcd[133] = 4.0E0*I_ERI_H2x2yz_Dxy_S_S_aa-2.0E0*1*I_ERI_F2yz_Dxy_S_S_a;
  abcd[134] = 4.0E0*I_ERI_H2xy2z_Dxy_S_S_aa-2.0E0*1*I_ERI_F2xy_Dxy_S_S_a-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_a+1*I_ERI_Py_Dxy_S_S;
  abcd[135] = 4.0E0*I_ERI_H2x3z_Dxy_S_S_aa-2.0E0*2*I_ERI_F2xz_Dxy_S_S_a-2.0E0*1*I_ERI_F3z_Dxy_S_S_a+2*I_ERI_Pz_Dxy_S_S;
  abcd[136] = 4.0E0*I_ERI_Hx3yz_Dxy_S_S_aa;
  abcd[137] = 4.0E0*I_ERI_Hx2y2z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_a;
  abcd[138] = 4.0E0*I_ERI_Hxy3z_Dxy_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_a;
  abcd[139] = 4.0E0*I_ERI_Hx4z_Dxy_S_S_aa-2.0E0*3*I_ERI_Fx2z_Dxy_S_S_a;
  abcd[140] = 4.0E0*I_ERI_H4xz_Dxz_S_S_aa-2.0E0*3*I_ERI_F2xz_Dxz_S_S_a;
  abcd[141] = 4.0E0*I_ERI_H3xyz_Dxz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_a;
  abcd[142] = 4.0E0*I_ERI_H3x2z_Dxz_S_S_aa-2.0E0*1*I_ERI_F3x_Dxz_S_S_a-2.0E0*2*I_ERI_Fx2z_Dxz_S_S_a+2*1*I_ERI_Px_Dxz_S_S;
  abcd[143] = 4.0E0*I_ERI_H2x2yz_Dxz_S_S_aa-2.0E0*1*I_ERI_F2yz_Dxz_S_S_a;
  abcd[144] = 4.0E0*I_ERI_H2xy2z_Dxz_S_S_aa-2.0E0*1*I_ERI_F2xy_Dxz_S_S_a-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_a+1*I_ERI_Py_Dxz_S_S;
  abcd[145] = 4.0E0*I_ERI_H2x3z_Dxz_S_S_aa-2.0E0*2*I_ERI_F2xz_Dxz_S_S_a-2.0E0*1*I_ERI_F3z_Dxz_S_S_a+2*I_ERI_Pz_Dxz_S_S;
  abcd[146] = 4.0E0*I_ERI_Hx3yz_Dxz_S_S_aa;
  abcd[147] = 4.0E0*I_ERI_Hx2y2z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_a;
  abcd[148] = 4.0E0*I_ERI_Hxy3z_Dxz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_a;
  abcd[149] = 4.0E0*I_ERI_Hx4z_Dxz_S_S_aa-2.0E0*3*I_ERI_Fx2z_Dxz_S_S_a;
  abcd[150] = 4.0E0*I_ERI_H4xz_D2y_S_S_aa-2.0E0*3*I_ERI_F2xz_D2y_S_S_a;
  abcd[151] = 4.0E0*I_ERI_H3xyz_D2y_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2y_S_S_a;
  abcd[152] = 4.0E0*I_ERI_H3x2z_D2y_S_S_aa-2.0E0*1*I_ERI_F3x_D2y_S_S_a-2.0E0*2*I_ERI_Fx2z_D2y_S_S_a+2*1*I_ERI_Px_D2y_S_S;
  abcd[153] = 4.0E0*I_ERI_H2x2yz_D2y_S_S_aa-2.0E0*1*I_ERI_F2yz_D2y_S_S_a;
  abcd[154] = 4.0E0*I_ERI_H2xy2z_D2y_S_S_aa-2.0E0*1*I_ERI_F2xy_D2y_S_S_a-2.0E0*1*I_ERI_Fy2z_D2y_S_S_a+1*I_ERI_Py_D2y_S_S;
  abcd[155] = 4.0E0*I_ERI_H2x3z_D2y_S_S_aa-2.0E0*2*I_ERI_F2xz_D2y_S_S_a-2.0E0*1*I_ERI_F3z_D2y_S_S_a+2*I_ERI_Pz_D2y_S_S;
  abcd[156] = 4.0E0*I_ERI_Hx3yz_D2y_S_S_aa;
  abcd[157] = 4.0E0*I_ERI_Hx2y2z_D2y_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2y_S_S_a;
  abcd[158] = 4.0E0*I_ERI_Hxy3z_D2y_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2y_S_S_a;
  abcd[159] = 4.0E0*I_ERI_Hx4z_D2y_S_S_aa-2.0E0*3*I_ERI_Fx2z_D2y_S_S_a;
  abcd[160] = 4.0E0*I_ERI_H4xz_Dyz_S_S_aa-2.0E0*3*I_ERI_F2xz_Dyz_S_S_a;
  abcd[161] = 4.0E0*I_ERI_H3xyz_Dyz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_a;
  abcd[162] = 4.0E0*I_ERI_H3x2z_Dyz_S_S_aa-2.0E0*1*I_ERI_F3x_Dyz_S_S_a-2.0E0*2*I_ERI_Fx2z_Dyz_S_S_a+2*1*I_ERI_Px_Dyz_S_S;
  abcd[163] = 4.0E0*I_ERI_H2x2yz_Dyz_S_S_aa-2.0E0*1*I_ERI_F2yz_Dyz_S_S_a;
  abcd[164] = 4.0E0*I_ERI_H2xy2z_Dyz_S_S_aa-2.0E0*1*I_ERI_F2xy_Dyz_S_S_a-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_a+1*I_ERI_Py_Dyz_S_S;
  abcd[165] = 4.0E0*I_ERI_H2x3z_Dyz_S_S_aa-2.0E0*2*I_ERI_F2xz_Dyz_S_S_a-2.0E0*1*I_ERI_F3z_Dyz_S_S_a+2*I_ERI_Pz_Dyz_S_S;
  abcd[166] = 4.0E0*I_ERI_Hx3yz_Dyz_S_S_aa;
  abcd[167] = 4.0E0*I_ERI_Hx2y2z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_a;
  abcd[168] = 4.0E0*I_ERI_Hxy3z_Dyz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_a;
  abcd[169] = 4.0E0*I_ERI_Hx4z_Dyz_S_S_aa-2.0E0*3*I_ERI_Fx2z_Dyz_S_S_a;
  abcd[170] = 4.0E0*I_ERI_H4xz_D2z_S_S_aa-2.0E0*3*I_ERI_F2xz_D2z_S_S_a;
  abcd[171] = 4.0E0*I_ERI_H3xyz_D2z_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2z_S_S_a;
  abcd[172] = 4.0E0*I_ERI_H3x2z_D2z_S_S_aa-2.0E0*1*I_ERI_F3x_D2z_S_S_a-2.0E0*2*I_ERI_Fx2z_D2z_S_S_a+2*1*I_ERI_Px_D2z_S_S;
  abcd[173] = 4.0E0*I_ERI_H2x2yz_D2z_S_S_aa-2.0E0*1*I_ERI_F2yz_D2z_S_S_a;
  abcd[174] = 4.0E0*I_ERI_H2xy2z_D2z_S_S_aa-2.0E0*1*I_ERI_F2xy_D2z_S_S_a-2.0E0*1*I_ERI_Fy2z_D2z_S_S_a+1*I_ERI_Py_D2z_S_S;
  abcd[175] = 4.0E0*I_ERI_H2x3z_D2z_S_S_aa-2.0E0*2*I_ERI_F2xz_D2z_S_S_a-2.0E0*1*I_ERI_F3z_D2z_S_S_a+2*I_ERI_Pz_D2z_S_S;
  abcd[176] = 4.0E0*I_ERI_Hx3yz_D2z_S_S_aa;
  abcd[177] = 4.0E0*I_ERI_Hx2y2z_D2z_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2z_S_S_a;
  abcd[178] = 4.0E0*I_ERI_Hxy3z_D2z_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2z_S_S_a;
  abcd[179] = 4.0E0*I_ERI_Hx4z_D2z_S_S_aa-2.0E0*3*I_ERI_Fx2z_D2z_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_S_S_aa
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S
   ************************************************************/
  abcd[180] = 4.0E0*I_ERI_H3x2y_D2x_S_S_aa-2.0E0*1*I_ERI_F3x_D2x_S_S_a;
  abcd[181] = 4.0E0*I_ERI_H2x3y_D2x_S_S_aa-2.0E0*1*I_ERI_F2xy_D2x_S_S_a-2.0E0*2*I_ERI_F2xy_D2x_S_S_a;
  abcd[182] = 4.0E0*I_ERI_H2x2yz_D2x_S_S_aa-2.0E0*1*I_ERI_F2xz_D2x_S_S_a;
  abcd[183] = 4.0E0*I_ERI_Hx4y_D2x_S_S_aa-2.0E0*2*I_ERI_Fx2y_D2x_S_S_a-2.0E0*3*I_ERI_Fx2y_D2x_S_S_a+2*1*I_ERI_Px_D2x_S_S;
  abcd[184] = 4.0E0*I_ERI_Hx3yz_D2x_S_S_aa-2.0E0*1*I_ERI_Fxyz_D2x_S_S_a-2.0E0*2*I_ERI_Fxyz_D2x_S_S_a;
  abcd[185] = 4.0E0*I_ERI_Hx2y2z_D2x_S_S_aa-2.0E0*1*I_ERI_Fx2z_D2x_S_S_a;
  abcd[186] = 4.0E0*I_ERI_H5y_D2x_S_S_aa-2.0E0*3*I_ERI_F3y_D2x_S_S_a-2.0E0*4*I_ERI_F3y_D2x_S_S_a+3*2*I_ERI_Py_D2x_S_S;
  abcd[187] = 4.0E0*I_ERI_H4yz_D2x_S_S_aa-2.0E0*2*I_ERI_F2yz_D2x_S_S_a-2.0E0*3*I_ERI_F2yz_D2x_S_S_a+2*1*I_ERI_Pz_D2x_S_S;
  abcd[188] = 4.0E0*I_ERI_H3y2z_D2x_S_S_aa-2.0E0*1*I_ERI_Fy2z_D2x_S_S_a-2.0E0*2*I_ERI_Fy2z_D2x_S_S_a;
  abcd[189] = 4.0E0*I_ERI_H2y3z_D2x_S_S_aa-2.0E0*1*I_ERI_F3z_D2x_S_S_a;
  abcd[190] = 4.0E0*I_ERI_H3x2y_Dxy_S_S_aa-2.0E0*1*I_ERI_F3x_Dxy_S_S_a;
  abcd[191] = 4.0E0*I_ERI_H2x3y_Dxy_S_S_aa-2.0E0*1*I_ERI_F2xy_Dxy_S_S_a-2.0E0*2*I_ERI_F2xy_Dxy_S_S_a;
  abcd[192] = 4.0E0*I_ERI_H2x2yz_Dxy_S_S_aa-2.0E0*1*I_ERI_F2xz_Dxy_S_S_a;
  abcd[193] = 4.0E0*I_ERI_Hx4y_Dxy_S_S_aa-2.0E0*2*I_ERI_Fx2y_Dxy_S_S_a-2.0E0*3*I_ERI_Fx2y_Dxy_S_S_a+2*1*I_ERI_Px_Dxy_S_S;
  abcd[194] = 4.0E0*I_ERI_Hx3yz_Dxy_S_S_aa-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_a-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_a;
  abcd[195] = 4.0E0*I_ERI_Hx2y2z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_a;
  abcd[196] = 4.0E0*I_ERI_H5y_Dxy_S_S_aa-2.0E0*3*I_ERI_F3y_Dxy_S_S_a-2.0E0*4*I_ERI_F3y_Dxy_S_S_a+3*2*I_ERI_Py_Dxy_S_S;
  abcd[197] = 4.0E0*I_ERI_H4yz_Dxy_S_S_aa-2.0E0*2*I_ERI_F2yz_Dxy_S_S_a-2.0E0*3*I_ERI_F2yz_Dxy_S_S_a+2*1*I_ERI_Pz_Dxy_S_S;
  abcd[198] = 4.0E0*I_ERI_H3y2z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_a-2.0E0*2*I_ERI_Fy2z_Dxy_S_S_a;
  abcd[199] = 4.0E0*I_ERI_H2y3z_Dxy_S_S_aa-2.0E0*1*I_ERI_F3z_Dxy_S_S_a;
  abcd[200] = 4.0E0*I_ERI_H3x2y_Dxz_S_S_aa-2.0E0*1*I_ERI_F3x_Dxz_S_S_a;
  abcd[201] = 4.0E0*I_ERI_H2x3y_Dxz_S_S_aa-2.0E0*1*I_ERI_F2xy_Dxz_S_S_a-2.0E0*2*I_ERI_F2xy_Dxz_S_S_a;
  abcd[202] = 4.0E0*I_ERI_H2x2yz_Dxz_S_S_aa-2.0E0*1*I_ERI_F2xz_Dxz_S_S_a;
  abcd[203] = 4.0E0*I_ERI_Hx4y_Dxz_S_S_aa-2.0E0*2*I_ERI_Fx2y_Dxz_S_S_a-2.0E0*3*I_ERI_Fx2y_Dxz_S_S_a+2*1*I_ERI_Px_Dxz_S_S;
  abcd[204] = 4.0E0*I_ERI_Hx3yz_Dxz_S_S_aa-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_a-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_a;
  abcd[205] = 4.0E0*I_ERI_Hx2y2z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_a;
  abcd[206] = 4.0E0*I_ERI_H5y_Dxz_S_S_aa-2.0E0*3*I_ERI_F3y_Dxz_S_S_a-2.0E0*4*I_ERI_F3y_Dxz_S_S_a+3*2*I_ERI_Py_Dxz_S_S;
  abcd[207] = 4.0E0*I_ERI_H4yz_Dxz_S_S_aa-2.0E0*2*I_ERI_F2yz_Dxz_S_S_a-2.0E0*3*I_ERI_F2yz_Dxz_S_S_a+2*1*I_ERI_Pz_Dxz_S_S;
  abcd[208] = 4.0E0*I_ERI_H3y2z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_a-2.0E0*2*I_ERI_Fy2z_Dxz_S_S_a;
  abcd[209] = 4.0E0*I_ERI_H2y3z_Dxz_S_S_aa-2.0E0*1*I_ERI_F3z_Dxz_S_S_a;
  abcd[210] = 4.0E0*I_ERI_H3x2y_D2y_S_S_aa-2.0E0*1*I_ERI_F3x_D2y_S_S_a;
  abcd[211] = 4.0E0*I_ERI_H2x3y_D2y_S_S_aa-2.0E0*1*I_ERI_F2xy_D2y_S_S_a-2.0E0*2*I_ERI_F2xy_D2y_S_S_a;
  abcd[212] = 4.0E0*I_ERI_H2x2yz_D2y_S_S_aa-2.0E0*1*I_ERI_F2xz_D2y_S_S_a;
  abcd[213] = 4.0E0*I_ERI_Hx4y_D2y_S_S_aa-2.0E0*2*I_ERI_Fx2y_D2y_S_S_a-2.0E0*3*I_ERI_Fx2y_D2y_S_S_a+2*1*I_ERI_Px_D2y_S_S;
  abcd[214] = 4.0E0*I_ERI_Hx3yz_D2y_S_S_aa-2.0E0*1*I_ERI_Fxyz_D2y_S_S_a-2.0E0*2*I_ERI_Fxyz_D2y_S_S_a;
  abcd[215] = 4.0E0*I_ERI_Hx2y2z_D2y_S_S_aa-2.0E0*1*I_ERI_Fx2z_D2y_S_S_a;
  abcd[216] = 4.0E0*I_ERI_H5y_D2y_S_S_aa-2.0E0*3*I_ERI_F3y_D2y_S_S_a-2.0E0*4*I_ERI_F3y_D2y_S_S_a+3*2*I_ERI_Py_D2y_S_S;
  abcd[217] = 4.0E0*I_ERI_H4yz_D2y_S_S_aa-2.0E0*2*I_ERI_F2yz_D2y_S_S_a-2.0E0*3*I_ERI_F2yz_D2y_S_S_a+2*1*I_ERI_Pz_D2y_S_S;
  abcd[218] = 4.0E0*I_ERI_H3y2z_D2y_S_S_aa-2.0E0*1*I_ERI_Fy2z_D2y_S_S_a-2.0E0*2*I_ERI_Fy2z_D2y_S_S_a;
  abcd[219] = 4.0E0*I_ERI_H2y3z_D2y_S_S_aa-2.0E0*1*I_ERI_F3z_D2y_S_S_a;
  abcd[220] = 4.0E0*I_ERI_H3x2y_Dyz_S_S_aa-2.0E0*1*I_ERI_F3x_Dyz_S_S_a;
  abcd[221] = 4.0E0*I_ERI_H2x3y_Dyz_S_S_aa-2.0E0*1*I_ERI_F2xy_Dyz_S_S_a-2.0E0*2*I_ERI_F2xy_Dyz_S_S_a;
  abcd[222] = 4.0E0*I_ERI_H2x2yz_Dyz_S_S_aa-2.0E0*1*I_ERI_F2xz_Dyz_S_S_a;
  abcd[223] = 4.0E0*I_ERI_Hx4y_Dyz_S_S_aa-2.0E0*2*I_ERI_Fx2y_Dyz_S_S_a-2.0E0*3*I_ERI_Fx2y_Dyz_S_S_a+2*1*I_ERI_Px_Dyz_S_S;
  abcd[224] = 4.0E0*I_ERI_Hx3yz_Dyz_S_S_aa-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_a-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_a;
  abcd[225] = 4.0E0*I_ERI_Hx2y2z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_a;
  abcd[226] = 4.0E0*I_ERI_H5y_Dyz_S_S_aa-2.0E0*3*I_ERI_F3y_Dyz_S_S_a-2.0E0*4*I_ERI_F3y_Dyz_S_S_a+3*2*I_ERI_Py_Dyz_S_S;
  abcd[227] = 4.0E0*I_ERI_H4yz_Dyz_S_S_aa-2.0E0*2*I_ERI_F2yz_Dyz_S_S_a-2.0E0*3*I_ERI_F2yz_Dyz_S_S_a+2*1*I_ERI_Pz_Dyz_S_S;
  abcd[228] = 4.0E0*I_ERI_H3y2z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_a-2.0E0*2*I_ERI_Fy2z_Dyz_S_S_a;
  abcd[229] = 4.0E0*I_ERI_H2y3z_Dyz_S_S_aa-2.0E0*1*I_ERI_F3z_Dyz_S_S_a;
  abcd[230] = 4.0E0*I_ERI_H3x2y_D2z_S_S_aa-2.0E0*1*I_ERI_F3x_D2z_S_S_a;
  abcd[231] = 4.0E0*I_ERI_H2x3y_D2z_S_S_aa-2.0E0*1*I_ERI_F2xy_D2z_S_S_a-2.0E0*2*I_ERI_F2xy_D2z_S_S_a;
  abcd[232] = 4.0E0*I_ERI_H2x2yz_D2z_S_S_aa-2.0E0*1*I_ERI_F2xz_D2z_S_S_a;
  abcd[233] = 4.0E0*I_ERI_Hx4y_D2z_S_S_aa-2.0E0*2*I_ERI_Fx2y_D2z_S_S_a-2.0E0*3*I_ERI_Fx2y_D2z_S_S_a+2*1*I_ERI_Px_D2z_S_S;
  abcd[234] = 4.0E0*I_ERI_Hx3yz_D2z_S_S_aa-2.0E0*1*I_ERI_Fxyz_D2z_S_S_a-2.0E0*2*I_ERI_Fxyz_D2z_S_S_a;
  abcd[235] = 4.0E0*I_ERI_Hx2y2z_D2z_S_S_aa-2.0E0*1*I_ERI_Fx2z_D2z_S_S_a;
  abcd[236] = 4.0E0*I_ERI_H5y_D2z_S_S_aa-2.0E0*3*I_ERI_F3y_D2z_S_S_a-2.0E0*4*I_ERI_F3y_D2z_S_S_a+3*2*I_ERI_Py_D2z_S_S;
  abcd[237] = 4.0E0*I_ERI_H4yz_D2z_S_S_aa-2.0E0*2*I_ERI_F2yz_D2z_S_S_a-2.0E0*3*I_ERI_F2yz_D2z_S_S_a+2*1*I_ERI_Pz_D2z_S_S;
  abcd[238] = 4.0E0*I_ERI_H3y2z_D2z_S_S_aa-2.0E0*1*I_ERI_Fy2z_D2z_S_S_a-2.0E0*2*I_ERI_Fy2z_D2z_S_S_a;
  abcd[239] = 4.0E0*I_ERI_H2y3z_D2z_S_S_aa-2.0E0*1*I_ERI_F3z_D2z_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_S_S_aa
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S
   ************************************************************/
  abcd[240] = 4.0E0*I_ERI_H3xyz_D2x_S_S_aa;
  abcd[241] = 4.0E0*I_ERI_H2x2yz_D2x_S_S_aa-2.0E0*1*I_ERI_F2xz_D2x_S_S_a;
  abcd[242] = 4.0E0*I_ERI_H2xy2z_D2x_S_S_aa-2.0E0*1*I_ERI_F2xy_D2x_S_S_a;
  abcd[243] = 4.0E0*I_ERI_Hx3yz_D2x_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2x_S_S_a;
  abcd[244] = 4.0E0*I_ERI_Hx2y2z_D2x_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2x_S_S_a-2.0E0*1*I_ERI_Fx2z_D2x_S_S_a+1*I_ERI_Px_D2x_S_S;
  abcd[245] = 4.0E0*I_ERI_Hxy3z_D2x_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2x_S_S_a;
  abcd[246] = 4.0E0*I_ERI_H4yz_D2x_S_S_aa-2.0E0*3*I_ERI_F2yz_D2x_S_S_a;
  abcd[247] = 4.0E0*I_ERI_H3y2z_D2x_S_S_aa-2.0E0*1*I_ERI_F3y_D2x_S_S_a-2.0E0*2*I_ERI_Fy2z_D2x_S_S_a+2*1*I_ERI_Py_D2x_S_S;
  abcd[248] = 4.0E0*I_ERI_H2y3z_D2x_S_S_aa-2.0E0*2*I_ERI_F2yz_D2x_S_S_a-2.0E0*1*I_ERI_F3z_D2x_S_S_a+2*I_ERI_Pz_D2x_S_S;
  abcd[249] = 4.0E0*I_ERI_Hy4z_D2x_S_S_aa-2.0E0*3*I_ERI_Fy2z_D2x_S_S_a;
  abcd[250] = 4.0E0*I_ERI_H3xyz_Dxy_S_S_aa;
  abcd[251] = 4.0E0*I_ERI_H2x2yz_Dxy_S_S_aa-2.0E0*1*I_ERI_F2xz_Dxy_S_S_a;
  abcd[252] = 4.0E0*I_ERI_H2xy2z_Dxy_S_S_aa-2.0E0*1*I_ERI_F2xy_Dxy_S_S_a;
  abcd[253] = 4.0E0*I_ERI_Hx3yz_Dxy_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_a;
  abcd[254] = 4.0E0*I_ERI_Hx2y2z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_a-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_a+1*I_ERI_Px_Dxy_S_S;
  abcd[255] = 4.0E0*I_ERI_Hxy3z_Dxy_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_a;
  abcd[256] = 4.0E0*I_ERI_H4yz_Dxy_S_S_aa-2.0E0*3*I_ERI_F2yz_Dxy_S_S_a;
  abcd[257] = 4.0E0*I_ERI_H3y2z_Dxy_S_S_aa-2.0E0*1*I_ERI_F3y_Dxy_S_S_a-2.0E0*2*I_ERI_Fy2z_Dxy_S_S_a+2*1*I_ERI_Py_Dxy_S_S;
  abcd[258] = 4.0E0*I_ERI_H2y3z_Dxy_S_S_aa-2.0E0*2*I_ERI_F2yz_Dxy_S_S_a-2.0E0*1*I_ERI_F3z_Dxy_S_S_a+2*I_ERI_Pz_Dxy_S_S;
  abcd[259] = 4.0E0*I_ERI_Hy4z_Dxy_S_S_aa-2.0E0*3*I_ERI_Fy2z_Dxy_S_S_a;
  abcd[260] = 4.0E0*I_ERI_H3xyz_Dxz_S_S_aa;
  abcd[261] = 4.0E0*I_ERI_H2x2yz_Dxz_S_S_aa-2.0E0*1*I_ERI_F2xz_Dxz_S_S_a;
  abcd[262] = 4.0E0*I_ERI_H2xy2z_Dxz_S_S_aa-2.0E0*1*I_ERI_F2xy_Dxz_S_S_a;
  abcd[263] = 4.0E0*I_ERI_Hx3yz_Dxz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_a;
  abcd[264] = 4.0E0*I_ERI_Hx2y2z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_a-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_a+1*I_ERI_Px_Dxz_S_S;
  abcd[265] = 4.0E0*I_ERI_Hxy3z_Dxz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_a;
  abcd[266] = 4.0E0*I_ERI_H4yz_Dxz_S_S_aa-2.0E0*3*I_ERI_F2yz_Dxz_S_S_a;
  abcd[267] = 4.0E0*I_ERI_H3y2z_Dxz_S_S_aa-2.0E0*1*I_ERI_F3y_Dxz_S_S_a-2.0E0*2*I_ERI_Fy2z_Dxz_S_S_a+2*1*I_ERI_Py_Dxz_S_S;
  abcd[268] = 4.0E0*I_ERI_H2y3z_Dxz_S_S_aa-2.0E0*2*I_ERI_F2yz_Dxz_S_S_a-2.0E0*1*I_ERI_F3z_Dxz_S_S_a+2*I_ERI_Pz_Dxz_S_S;
  abcd[269] = 4.0E0*I_ERI_Hy4z_Dxz_S_S_aa-2.0E0*3*I_ERI_Fy2z_Dxz_S_S_a;
  abcd[270] = 4.0E0*I_ERI_H3xyz_D2y_S_S_aa;
  abcd[271] = 4.0E0*I_ERI_H2x2yz_D2y_S_S_aa-2.0E0*1*I_ERI_F2xz_D2y_S_S_a;
  abcd[272] = 4.0E0*I_ERI_H2xy2z_D2y_S_S_aa-2.0E0*1*I_ERI_F2xy_D2y_S_S_a;
  abcd[273] = 4.0E0*I_ERI_Hx3yz_D2y_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2y_S_S_a;
  abcd[274] = 4.0E0*I_ERI_Hx2y2z_D2y_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2y_S_S_a-2.0E0*1*I_ERI_Fx2z_D2y_S_S_a+1*I_ERI_Px_D2y_S_S;
  abcd[275] = 4.0E0*I_ERI_Hxy3z_D2y_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2y_S_S_a;
  abcd[276] = 4.0E0*I_ERI_H4yz_D2y_S_S_aa-2.0E0*3*I_ERI_F2yz_D2y_S_S_a;
  abcd[277] = 4.0E0*I_ERI_H3y2z_D2y_S_S_aa-2.0E0*1*I_ERI_F3y_D2y_S_S_a-2.0E0*2*I_ERI_Fy2z_D2y_S_S_a+2*1*I_ERI_Py_D2y_S_S;
  abcd[278] = 4.0E0*I_ERI_H2y3z_D2y_S_S_aa-2.0E0*2*I_ERI_F2yz_D2y_S_S_a-2.0E0*1*I_ERI_F3z_D2y_S_S_a+2*I_ERI_Pz_D2y_S_S;
  abcd[279] = 4.0E0*I_ERI_Hy4z_D2y_S_S_aa-2.0E0*3*I_ERI_Fy2z_D2y_S_S_a;
  abcd[280] = 4.0E0*I_ERI_H3xyz_Dyz_S_S_aa;
  abcd[281] = 4.0E0*I_ERI_H2x2yz_Dyz_S_S_aa-2.0E0*1*I_ERI_F2xz_Dyz_S_S_a;
  abcd[282] = 4.0E0*I_ERI_H2xy2z_Dyz_S_S_aa-2.0E0*1*I_ERI_F2xy_Dyz_S_S_a;
  abcd[283] = 4.0E0*I_ERI_Hx3yz_Dyz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_a;
  abcd[284] = 4.0E0*I_ERI_Hx2y2z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_a-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_a+1*I_ERI_Px_Dyz_S_S;
  abcd[285] = 4.0E0*I_ERI_Hxy3z_Dyz_S_S_aa-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_a;
  abcd[286] = 4.0E0*I_ERI_H4yz_Dyz_S_S_aa-2.0E0*3*I_ERI_F2yz_Dyz_S_S_a;
  abcd[287] = 4.0E0*I_ERI_H3y2z_Dyz_S_S_aa-2.0E0*1*I_ERI_F3y_Dyz_S_S_a-2.0E0*2*I_ERI_Fy2z_Dyz_S_S_a+2*1*I_ERI_Py_Dyz_S_S;
  abcd[288] = 4.0E0*I_ERI_H2y3z_Dyz_S_S_aa-2.0E0*2*I_ERI_F2yz_Dyz_S_S_a-2.0E0*1*I_ERI_F3z_Dyz_S_S_a+2*I_ERI_Pz_Dyz_S_S;
  abcd[289] = 4.0E0*I_ERI_Hy4z_Dyz_S_S_aa-2.0E0*3*I_ERI_Fy2z_Dyz_S_S_a;
  abcd[290] = 4.0E0*I_ERI_H3xyz_D2z_S_S_aa;
  abcd[291] = 4.0E0*I_ERI_H2x2yz_D2z_S_S_aa-2.0E0*1*I_ERI_F2xz_D2z_S_S_a;
  abcd[292] = 4.0E0*I_ERI_H2xy2z_D2z_S_S_aa-2.0E0*1*I_ERI_F2xy_D2z_S_S_a;
  abcd[293] = 4.0E0*I_ERI_Hx3yz_D2z_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2z_S_S_a;
  abcd[294] = 4.0E0*I_ERI_Hx2y2z_D2z_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2z_S_S_a-2.0E0*1*I_ERI_Fx2z_D2z_S_S_a+1*I_ERI_Px_D2z_S_S;
  abcd[295] = 4.0E0*I_ERI_Hxy3z_D2z_S_S_aa-2.0E0*2*I_ERI_Fxyz_D2z_S_S_a;
  abcd[296] = 4.0E0*I_ERI_H4yz_D2z_S_S_aa-2.0E0*3*I_ERI_F2yz_D2z_S_S_a;
  abcd[297] = 4.0E0*I_ERI_H3y2z_D2z_S_S_aa-2.0E0*1*I_ERI_F3y_D2z_S_S_a-2.0E0*2*I_ERI_Fy2z_D2z_S_S_a+2*1*I_ERI_Py_D2z_S_S;
  abcd[298] = 4.0E0*I_ERI_H2y3z_D2z_S_S_aa-2.0E0*2*I_ERI_F2yz_D2z_S_S_a-2.0E0*1*I_ERI_F3z_D2z_S_S_a+2*I_ERI_Pz_D2z_S_S;
  abcd[299] = 4.0E0*I_ERI_Hy4z_D2z_S_S_aa-2.0E0*3*I_ERI_Fy2z_D2z_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_S_S_aa
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_F_D_S_S_a
   * RHS shell quartet name: SQ_ERI_P_D_S_S
   ************************************************************/
  abcd[300] = 4.0E0*I_ERI_H3x2z_D2x_S_S_aa-2.0E0*1*I_ERI_F3x_D2x_S_S_a;
  abcd[301] = 4.0E0*I_ERI_H2xy2z_D2x_S_S_aa-2.0E0*1*I_ERI_F2xy_D2x_S_S_a;
  abcd[302] = 4.0E0*I_ERI_H2x3z_D2x_S_S_aa-2.0E0*1*I_ERI_F2xz_D2x_S_S_a-2.0E0*2*I_ERI_F2xz_D2x_S_S_a;
  abcd[303] = 4.0E0*I_ERI_Hx2y2z_D2x_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2x_S_S_a;
  abcd[304] = 4.0E0*I_ERI_Hxy3z_D2x_S_S_aa-2.0E0*1*I_ERI_Fxyz_D2x_S_S_a-2.0E0*2*I_ERI_Fxyz_D2x_S_S_a;
  abcd[305] = 4.0E0*I_ERI_Hx4z_D2x_S_S_aa-2.0E0*2*I_ERI_Fx2z_D2x_S_S_a-2.0E0*3*I_ERI_Fx2z_D2x_S_S_a+2*1*I_ERI_Px_D2x_S_S;
  abcd[306] = 4.0E0*I_ERI_H3y2z_D2x_S_S_aa-2.0E0*1*I_ERI_F3y_D2x_S_S_a;
  abcd[307] = 4.0E0*I_ERI_H2y3z_D2x_S_S_aa-2.0E0*1*I_ERI_F2yz_D2x_S_S_a-2.0E0*2*I_ERI_F2yz_D2x_S_S_a;
  abcd[308] = 4.0E0*I_ERI_Hy4z_D2x_S_S_aa-2.0E0*2*I_ERI_Fy2z_D2x_S_S_a-2.0E0*3*I_ERI_Fy2z_D2x_S_S_a+2*1*I_ERI_Py_D2x_S_S;
  abcd[309] = 4.0E0*I_ERI_H5z_D2x_S_S_aa-2.0E0*3*I_ERI_F3z_D2x_S_S_a-2.0E0*4*I_ERI_F3z_D2x_S_S_a+3*2*I_ERI_Pz_D2x_S_S;
  abcd[310] = 4.0E0*I_ERI_H3x2z_Dxy_S_S_aa-2.0E0*1*I_ERI_F3x_Dxy_S_S_a;
  abcd[311] = 4.0E0*I_ERI_H2xy2z_Dxy_S_S_aa-2.0E0*1*I_ERI_F2xy_Dxy_S_S_a;
  abcd[312] = 4.0E0*I_ERI_H2x3z_Dxy_S_S_aa-2.0E0*1*I_ERI_F2xz_Dxy_S_S_a-2.0E0*2*I_ERI_F2xz_Dxy_S_S_a;
  abcd[313] = 4.0E0*I_ERI_Hx2y2z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_a;
  abcd[314] = 4.0E0*I_ERI_Hxy3z_Dxy_S_S_aa-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_a-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_a;
  abcd[315] = 4.0E0*I_ERI_Hx4z_Dxy_S_S_aa-2.0E0*2*I_ERI_Fx2z_Dxy_S_S_a-2.0E0*3*I_ERI_Fx2z_Dxy_S_S_a+2*1*I_ERI_Px_Dxy_S_S;
  abcd[316] = 4.0E0*I_ERI_H3y2z_Dxy_S_S_aa-2.0E0*1*I_ERI_F3y_Dxy_S_S_a;
  abcd[317] = 4.0E0*I_ERI_H2y3z_Dxy_S_S_aa-2.0E0*1*I_ERI_F2yz_Dxy_S_S_a-2.0E0*2*I_ERI_F2yz_Dxy_S_S_a;
  abcd[318] = 4.0E0*I_ERI_Hy4z_Dxy_S_S_aa-2.0E0*2*I_ERI_Fy2z_Dxy_S_S_a-2.0E0*3*I_ERI_Fy2z_Dxy_S_S_a+2*1*I_ERI_Py_Dxy_S_S;
  abcd[319] = 4.0E0*I_ERI_H5z_Dxy_S_S_aa-2.0E0*3*I_ERI_F3z_Dxy_S_S_a-2.0E0*4*I_ERI_F3z_Dxy_S_S_a+3*2*I_ERI_Pz_Dxy_S_S;
  abcd[320] = 4.0E0*I_ERI_H3x2z_Dxz_S_S_aa-2.0E0*1*I_ERI_F3x_Dxz_S_S_a;
  abcd[321] = 4.0E0*I_ERI_H2xy2z_Dxz_S_S_aa-2.0E0*1*I_ERI_F2xy_Dxz_S_S_a;
  abcd[322] = 4.0E0*I_ERI_H2x3z_Dxz_S_S_aa-2.0E0*1*I_ERI_F2xz_Dxz_S_S_a-2.0E0*2*I_ERI_F2xz_Dxz_S_S_a;
  abcd[323] = 4.0E0*I_ERI_Hx2y2z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_a;
  abcd[324] = 4.0E0*I_ERI_Hxy3z_Dxz_S_S_aa-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_a-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_a;
  abcd[325] = 4.0E0*I_ERI_Hx4z_Dxz_S_S_aa-2.0E0*2*I_ERI_Fx2z_Dxz_S_S_a-2.0E0*3*I_ERI_Fx2z_Dxz_S_S_a+2*1*I_ERI_Px_Dxz_S_S;
  abcd[326] = 4.0E0*I_ERI_H3y2z_Dxz_S_S_aa-2.0E0*1*I_ERI_F3y_Dxz_S_S_a;
  abcd[327] = 4.0E0*I_ERI_H2y3z_Dxz_S_S_aa-2.0E0*1*I_ERI_F2yz_Dxz_S_S_a-2.0E0*2*I_ERI_F2yz_Dxz_S_S_a;
  abcd[328] = 4.0E0*I_ERI_Hy4z_Dxz_S_S_aa-2.0E0*2*I_ERI_Fy2z_Dxz_S_S_a-2.0E0*3*I_ERI_Fy2z_Dxz_S_S_a+2*1*I_ERI_Py_Dxz_S_S;
  abcd[329] = 4.0E0*I_ERI_H5z_Dxz_S_S_aa-2.0E0*3*I_ERI_F3z_Dxz_S_S_a-2.0E0*4*I_ERI_F3z_Dxz_S_S_a+3*2*I_ERI_Pz_Dxz_S_S;
  abcd[330] = 4.0E0*I_ERI_H3x2z_D2y_S_S_aa-2.0E0*1*I_ERI_F3x_D2y_S_S_a;
  abcd[331] = 4.0E0*I_ERI_H2xy2z_D2y_S_S_aa-2.0E0*1*I_ERI_F2xy_D2y_S_S_a;
  abcd[332] = 4.0E0*I_ERI_H2x3z_D2y_S_S_aa-2.0E0*1*I_ERI_F2xz_D2y_S_S_a-2.0E0*2*I_ERI_F2xz_D2y_S_S_a;
  abcd[333] = 4.0E0*I_ERI_Hx2y2z_D2y_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2y_S_S_a;
  abcd[334] = 4.0E0*I_ERI_Hxy3z_D2y_S_S_aa-2.0E0*1*I_ERI_Fxyz_D2y_S_S_a-2.0E0*2*I_ERI_Fxyz_D2y_S_S_a;
  abcd[335] = 4.0E0*I_ERI_Hx4z_D2y_S_S_aa-2.0E0*2*I_ERI_Fx2z_D2y_S_S_a-2.0E0*3*I_ERI_Fx2z_D2y_S_S_a+2*1*I_ERI_Px_D2y_S_S;
  abcd[336] = 4.0E0*I_ERI_H3y2z_D2y_S_S_aa-2.0E0*1*I_ERI_F3y_D2y_S_S_a;
  abcd[337] = 4.0E0*I_ERI_H2y3z_D2y_S_S_aa-2.0E0*1*I_ERI_F2yz_D2y_S_S_a-2.0E0*2*I_ERI_F2yz_D2y_S_S_a;
  abcd[338] = 4.0E0*I_ERI_Hy4z_D2y_S_S_aa-2.0E0*2*I_ERI_Fy2z_D2y_S_S_a-2.0E0*3*I_ERI_Fy2z_D2y_S_S_a+2*1*I_ERI_Py_D2y_S_S;
  abcd[339] = 4.0E0*I_ERI_H5z_D2y_S_S_aa-2.0E0*3*I_ERI_F3z_D2y_S_S_a-2.0E0*4*I_ERI_F3z_D2y_S_S_a+3*2*I_ERI_Pz_D2y_S_S;
  abcd[340] = 4.0E0*I_ERI_H3x2z_Dyz_S_S_aa-2.0E0*1*I_ERI_F3x_Dyz_S_S_a;
  abcd[341] = 4.0E0*I_ERI_H2xy2z_Dyz_S_S_aa-2.0E0*1*I_ERI_F2xy_Dyz_S_S_a;
  abcd[342] = 4.0E0*I_ERI_H2x3z_Dyz_S_S_aa-2.0E0*1*I_ERI_F2xz_Dyz_S_S_a-2.0E0*2*I_ERI_F2xz_Dyz_S_S_a;
  abcd[343] = 4.0E0*I_ERI_Hx2y2z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_a;
  abcd[344] = 4.0E0*I_ERI_Hxy3z_Dyz_S_S_aa-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_a-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_a;
  abcd[345] = 4.0E0*I_ERI_Hx4z_Dyz_S_S_aa-2.0E0*2*I_ERI_Fx2z_Dyz_S_S_a-2.0E0*3*I_ERI_Fx2z_Dyz_S_S_a+2*1*I_ERI_Px_Dyz_S_S;
  abcd[346] = 4.0E0*I_ERI_H3y2z_Dyz_S_S_aa-2.0E0*1*I_ERI_F3y_Dyz_S_S_a;
  abcd[347] = 4.0E0*I_ERI_H2y3z_Dyz_S_S_aa-2.0E0*1*I_ERI_F2yz_Dyz_S_S_a-2.0E0*2*I_ERI_F2yz_Dyz_S_S_a;
  abcd[348] = 4.0E0*I_ERI_Hy4z_Dyz_S_S_aa-2.0E0*2*I_ERI_Fy2z_Dyz_S_S_a-2.0E0*3*I_ERI_Fy2z_Dyz_S_S_a+2*1*I_ERI_Py_Dyz_S_S;
  abcd[349] = 4.0E0*I_ERI_H5z_Dyz_S_S_aa-2.0E0*3*I_ERI_F3z_Dyz_S_S_a-2.0E0*4*I_ERI_F3z_Dyz_S_S_a+3*2*I_ERI_Pz_Dyz_S_S;
  abcd[350] = 4.0E0*I_ERI_H3x2z_D2z_S_S_aa-2.0E0*1*I_ERI_F3x_D2z_S_S_a;
  abcd[351] = 4.0E0*I_ERI_H2xy2z_D2z_S_S_aa-2.0E0*1*I_ERI_F2xy_D2z_S_S_a;
  abcd[352] = 4.0E0*I_ERI_H2x3z_D2z_S_S_aa-2.0E0*1*I_ERI_F2xz_D2z_S_S_a-2.0E0*2*I_ERI_F2xz_D2z_S_S_a;
  abcd[353] = 4.0E0*I_ERI_Hx2y2z_D2z_S_S_aa-2.0E0*1*I_ERI_Fx2y_D2z_S_S_a;
  abcd[354] = 4.0E0*I_ERI_Hxy3z_D2z_S_S_aa-2.0E0*1*I_ERI_Fxyz_D2z_S_S_a-2.0E0*2*I_ERI_Fxyz_D2z_S_S_a;
  abcd[355] = 4.0E0*I_ERI_Hx4z_D2z_S_S_aa-2.0E0*2*I_ERI_Fx2z_D2z_S_S_a-2.0E0*3*I_ERI_Fx2z_D2z_S_S_a+2*1*I_ERI_Px_D2z_S_S;
  abcd[356] = 4.0E0*I_ERI_H3y2z_D2z_S_S_aa-2.0E0*1*I_ERI_F3y_D2z_S_S_a;
  abcd[357] = 4.0E0*I_ERI_H2y3z_D2z_S_S_aa-2.0E0*1*I_ERI_F2yz_D2z_S_S_a-2.0E0*2*I_ERI_F2yz_D2z_S_S_a;
  abcd[358] = 4.0E0*I_ERI_Hy4z_D2z_S_S_aa-2.0E0*2*I_ERI_Fy2z_D2z_S_S_a-2.0E0*3*I_ERI_Fy2z_D2z_S_S_a+2*1*I_ERI_Py_D2z_S_S;
  abcd[359] = 4.0E0*I_ERI_H5z_D2z_S_S_aa-2.0E0*3*I_ERI_F3z_D2z_S_S_a-2.0E0*4*I_ERI_F3z_D2z_S_S_a+3*2*I_ERI_Pz_D2z_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_D_F_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  abcd[360] = 4.0E0*I_ERI_G4x_F3x_S_S_ab-2.0E0*2*I_ERI_G4x_Px_S_S_a-2.0E0*3*I_ERI_D2x_F3x_S_S_b+3*2*I_ERI_D2x_Px_S_S;
  abcd[361] = 4.0E0*I_ERI_G3xy_F3x_S_S_ab-2.0E0*2*I_ERI_G3xy_Px_S_S_a-2.0E0*2*I_ERI_Dxy_F3x_S_S_b+2*2*I_ERI_Dxy_Px_S_S;
  abcd[362] = 4.0E0*I_ERI_G3xz_F3x_S_S_ab-2.0E0*2*I_ERI_G3xz_Px_S_S_a-2.0E0*2*I_ERI_Dxz_F3x_S_S_b+2*2*I_ERI_Dxz_Px_S_S;
  abcd[363] = 4.0E0*I_ERI_G2x2y_F3x_S_S_ab-2.0E0*2*I_ERI_G2x2y_Px_S_S_a-2.0E0*1*I_ERI_D2y_F3x_S_S_b+2*I_ERI_D2y_Px_S_S;
  abcd[364] = 4.0E0*I_ERI_G2xyz_F3x_S_S_ab-2.0E0*2*I_ERI_G2xyz_Px_S_S_a-2.0E0*1*I_ERI_Dyz_F3x_S_S_b+2*I_ERI_Dyz_Px_S_S;
  abcd[365] = 4.0E0*I_ERI_G2x2z_F3x_S_S_ab-2.0E0*2*I_ERI_G2x2z_Px_S_S_a-2.0E0*1*I_ERI_D2z_F3x_S_S_b+2*I_ERI_D2z_Px_S_S;
  abcd[366] = 4.0E0*I_ERI_Gx3y_F3x_S_S_ab-2.0E0*2*I_ERI_Gx3y_Px_S_S_a;
  abcd[367] = 4.0E0*I_ERI_Gx2yz_F3x_S_S_ab-2.0E0*2*I_ERI_Gx2yz_Px_S_S_a;
  abcd[368] = 4.0E0*I_ERI_Gxy2z_F3x_S_S_ab-2.0E0*2*I_ERI_Gxy2z_Px_S_S_a;
  abcd[369] = 4.0E0*I_ERI_Gx3z_F3x_S_S_ab-2.0E0*2*I_ERI_Gx3z_Px_S_S_a;
  abcd[370] = 4.0E0*I_ERI_G4x_F2xy_S_S_ab-2.0E0*1*I_ERI_G4x_Py_S_S_a-2.0E0*3*I_ERI_D2x_F2xy_S_S_b+3*1*I_ERI_D2x_Py_S_S;
  abcd[371] = 4.0E0*I_ERI_G3xy_F2xy_S_S_ab-2.0E0*1*I_ERI_G3xy_Py_S_S_a-2.0E0*2*I_ERI_Dxy_F2xy_S_S_b+2*1*I_ERI_Dxy_Py_S_S;
  abcd[372] = 4.0E0*I_ERI_G3xz_F2xy_S_S_ab-2.0E0*1*I_ERI_G3xz_Py_S_S_a-2.0E0*2*I_ERI_Dxz_F2xy_S_S_b+2*1*I_ERI_Dxz_Py_S_S;
  abcd[373] = 4.0E0*I_ERI_G2x2y_F2xy_S_S_ab-2.0E0*1*I_ERI_G2x2y_Py_S_S_a-2.0E0*1*I_ERI_D2y_F2xy_S_S_b+1*I_ERI_D2y_Py_S_S;
  abcd[374] = 4.0E0*I_ERI_G2xyz_F2xy_S_S_ab-2.0E0*1*I_ERI_G2xyz_Py_S_S_a-2.0E0*1*I_ERI_Dyz_F2xy_S_S_b+1*I_ERI_Dyz_Py_S_S;
  abcd[375] = 4.0E0*I_ERI_G2x2z_F2xy_S_S_ab-2.0E0*1*I_ERI_G2x2z_Py_S_S_a-2.0E0*1*I_ERI_D2z_F2xy_S_S_b+1*I_ERI_D2z_Py_S_S;
  abcd[376] = 4.0E0*I_ERI_Gx3y_F2xy_S_S_ab-2.0E0*1*I_ERI_Gx3y_Py_S_S_a;
  abcd[377] = 4.0E0*I_ERI_Gx2yz_F2xy_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Py_S_S_a;
  abcd[378] = 4.0E0*I_ERI_Gxy2z_F2xy_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Py_S_S_a;
  abcd[379] = 4.0E0*I_ERI_Gx3z_F2xy_S_S_ab-2.0E0*1*I_ERI_Gx3z_Py_S_S_a;
  abcd[380] = 4.0E0*I_ERI_G4x_F2xz_S_S_ab-2.0E0*1*I_ERI_G4x_Pz_S_S_a-2.0E0*3*I_ERI_D2x_F2xz_S_S_b+3*1*I_ERI_D2x_Pz_S_S;
  abcd[381] = 4.0E0*I_ERI_G3xy_F2xz_S_S_ab-2.0E0*1*I_ERI_G3xy_Pz_S_S_a-2.0E0*2*I_ERI_Dxy_F2xz_S_S_b+2*1*I_ERI_Dxy_Pz_S_S;
  abcd[382] = 4.0E0*I_ERI_G3xz_F2xz_S_S_ab-2.0E0*1*I_ERI_G3xz_Pz_S_S_a-2.0E0*2*I_ERI_Dxz_F2xz_S_S_b+2*1*I_ERI_Dxz_Pz_S_S;
  abcd[383] = 4.0E0*I_ERI_G2x2y_F2xz_S_S_ab-2.0E0*1*I_ERI_G2x2y_Pz_S_S_a-2.0E0*1*I_ERI_D2y_F2xz_S_S_b+1*I_ERI_D2y_Pz_S_S;
  abcd[384] = 4.0E0*I_ERI_G2xyz_F2xz_S_S_ab-2.0E0*1*I_ERI_G2xyz_Pz_S_S_a-2.0E0*1*I_ERI_Dyz_F2xz_S_S_b+1*I_ERI_Dyz_Pz_S_S;
  abcd[385] = 4.0E0*I_ERI_G2x2z_F2xz_S_S_ab-2.0E0*1*I_ERI_G2x2z_Pz_S_S_a-2.0E0*1*I_ERI_D2z_F2xz_S_S_b+1*I_ERI_D2z_Pz_S_S;
  abcd[386] = 4.0E0*I_ERI_Gx3y_F2xz_S_S_ab-2.0E0*1*I_ERI_Gx3y_Pz_S_S_a;
  abcd[387] = 4.0E0*I_ERI_Gx2yz_F2xz_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Pz_S_S_a;
  abcd[388] = 4.0E0*I_ERI_Gxy2z_F2xz_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Pz_S_S_a;
  abcd[389] = 4.0E0*I_ERI_Gx3z_F2xz_S_S_ab-2.0E0*1*I_ERI_Gx3z_Pz_S_S_a;
  abcd[390] = 4.0E0*I_ERI_G4x_Fx2y_S_S_ab-2.0E0*3*I_ERI_D2x_Fx2y_S_S_b;
  abcd[391] = 4.0E0*I_ERI_G3xy_Fx2y_S_S_ab-2.0E0*2*I_ERI_Dxy_Fx2y_S_S_b;
  abcd[392] = 4.0E0*I_ERI_G3xz_Fx2y_S_S_ab-2.0E0*2*I_ERI_Dxz_Fx2y_S_S_b;
  abcd[393] = 4.0E0*I_ERI_G2x2y_Fx2y_S_S_ab-2.0E0*1*I_ERI_D2y_Fx2y_S_S_b;
  abcd[394] = 4.0E0*I_ERI_G2xyz_Fx2y_S_S_ab-2.0E0*1*I_ERI_Dyz_Fx2y_S_S_b;
  abcd[395] = 4.0E0*I_ERI_G2x2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_D2z_Fx2y_S_S_b;
  abcd[396] = 4.0E0*I_ERI_Gx3y_Fx2y_S_S_ab;
  abcd[397] = 4.0E0*I_ERI_Gx2yz_Fx2y_S_S_ab;
  abcd[398] = 4.0E0*I_ERI_Gxy2z_Fx2y_S_S_ab;
  abcd[399] = 4.0E0*I_ERI_Gx3z_Fx2y_S_S_ab;
  abcd[400] = 4.0E0*I_ERI_G4x_Fxyz_S_S_ab-2.0E0*3*I_ERI_D2x_Fxyz_S_S_b;
  abcd[401] = 4.0E0*I_ERI_G3xy_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxy_Fxyz_S_S_b;
  abcd[402] = 4.0E0*I_ERI_G3xz_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxz_Fxyz_S_S_b;
  abcd[403] = 4.0E0*I_ERI_G2x2y_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2y_Fxyz_S_S_b;
  abcd[404] = 4.0E0*I_ERI_G2xyz_Fxyz_S_S_ab-2.0E0*1*I_ERI_Dyz_Fxyz_S_S_b;
  abcd[405] = 4.0E0*I_ERI_G2x2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2z_Fxyz_S_S_b;
  abcd[406] = 4.0E0*I_ERI_Gx3y_Fxyz_S_S_ab;
  abcd[407] = 4.0E0*I_ERI_Gx2yz_Fxyz_S_S_ab;
  abcd[408] = 4.0E0*I_ERI_Gxy2z_Fxyz_S_S_ab;
  abcd[409] = 4.0E0*I_ERI_Gx3z_Fxyz_S_S_ab;
  abcd[410] = 4.0E0*I_ERI_G4x_Fx2z_S_S_ab-2.0E0*3*I_ERI_D2x_Fx2z_S_S_b;
  abcd[411] = 4.0E0*I_ERI_G3xy_Fx2z_S_S_ab-2.0E0*2*I_ERI_Dxy_Fx2z_S_S_b;
  abcd[412] = 4.0E0*I_ERI_G3xz_Fx2z_S_S_ab-2.0E0*2*I_ERI_Dxz_Fx2z_S_S_b;
  abcd[413] = 4.0E0*I_ERI_G2x2y_Fx2z_S_S_ab-2.0E0*1*I_ERI_D2y_Fx2z_S_S_b;
  abcd[414] = 4.0E0*I_ERI_G2xyz_Fx2z_S_S_ab-2.0E0*1*I_ERI_Dyz_Fx2z_S_S_b;
  abcd[415] = 4.0E0*I_ERI_G2x2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_D2z_Fx2z_S_S_b;
  abcd[416] = 4.0E0*I_ERI_Gx3y_Fx2z_S_S_ab;
  abcd[417] = 4.0E0*I_ERI_Gx2yz_Fx2z_S_S_ab;
  abcd[418] = 4.0E0*I_ERI_Gxy2z_Fx2z_S_S_ab;
  abcd[419] = 4.0E0*I_ERI_Gx3z_Fx2z_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_D_F_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  abcd[420] = 4.0E0*I_ERI_G4x_F2xy_S_S_ab-2.0E0*3*I_ERI_D2x_F2xy_S_S_b;
  abcd[421] = 4.0E0*I_ERI_G3xy_F2xy_S_S_ab-2.0E0*2*I_ERI_Dxy_F2xy_S_S_b;
  abcd[422] = 4.0E0*I_ERI_G3xz_F2xy_S_S_ab-2.0E0*2*I_ERI_Dxz_F2xy_S_S_b;
  abcd[423] = 4.0E0*I_ERI_G2x2y_F2xy_S_S_ab-2.0E0*1*I_ERI_D2y_F2xy_S_S_b;
  abcd[424] = 4.0E0*I_ERI_G2xyz_F2xy_S_S_ab-2.0E0*1*I_ERI_Dyz_F2xy_S_S_b;
  abcd[425] = 4.0E0*I_ERI_G2x2z_F2xy_S_S_ab-2.0E0*1*I_ERI_D2z_F2xy_S_S_b;
  abcd[426] = 4.0E0*I_ERI_Gx3y_F2xy_S_S_ab;
  abcd[427] = 4.0E0*I_ERI_Gx2yz_F2xy_S_S_ab;
  abcd[428] = 4.0E0*I_ERI_Gxy2z_F2xy_S_S_ab;
  abcd[429] = 4.0E0*I_ERI_Gx3z_F2xy_S_S_ab;
  abcd[430] = 4.0E0*I_ERI_G4x_Fx2y_S_S_ab-2.0E0*1*I_ERI_G4x_Px_S_S_a-2.0E0*3*I_ERI_D2x_Fx2y_S_S_b+3*1*I_ERI_D2x_Px_S_S;
  abcd[431] = 4.0E0*I_ERI_G3xy_Fx2y_S_S_ab-2.0E0*1*I_ERI_G3xy_Px_S_S_a-2.0E0*2*I_ERI_Dxy_Fx2y_S_S_b+2*1*I_ERI_Dxy_Px_S_S;
  abcd[432] = 4.0E0*I_ERI_G3xz_Fx2y_S_S_ab-2.0E0*1*I_ERI_G3xz_Px_S_S_a-2.0E0*2*I_ERI_Dxz_Fx2y_S_S_b+2*1*I_ERI_Dxz_Px_S_S;
  abcd[433] = 4.0E0*I_ERI_G2x2y_Fx2y_S_S_ab-2.0E0*1*I_ERI_G2x2y_Px_S_S_a-2.0E0*1*I_ERI_D2y_Fx2y_S_S_b+1*I_ERI_D2y_Px_S_S;
  abcd[434] = 4.0E0*I_ERI_G2xyz_Fx2y_S_S_ab-2.0E0*1*I_ERI_G2xyz_Px_S_S_a-2.0E0*1*I_ERI_Dyz_Fx2y_S_S_b+1*I_ERI_Dyz_Px_S_S;
  abcd[435] = 4.0E0*I_ERI_G2x2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_G2x2z_Px_S_S_a-2.0E0*1*I_ERI_D2z_Fx2y_S_S_b+1*I_ERI_D2z_Px_S_S;
  abcd[436] = 4.0E0*I_ERI_Gx3y_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gx3y_Px_S_S_a;
  abcd[437] = 4.0E0*I_ERI_Gx2yz_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Px_S_S_a;
  abcd[438] = 4.0E0*I_ERI_Gxy2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Px_S_S_a;
  abcd[439] = 4.0E0*I_ERI_Gx3z_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gx3z_Px_S_S_a;
  abcd[440] = 4.0E0*I_ERI_G4x_Fxyz_S_S_ab-2.0E0*3*I_ERI_D2x_Fxyz_S_S_b;
  abcd[441] = 4.0E0*I_ERI_G3xy_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxy_Fxyz_S_S_b;
  abcd[442] = 4.0E0*I_ERI_G3xz_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxz_Fxyz_S_S_b;
  abcd[443] = 4.0E0*I_ERI_G2x2y_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2y_Fxyz_S_S_b;
  abcd[444] = 4.0E0*I_ERI_G2xyz_Fxyz_S_S_ab-2.0E0*1*I_ERI_Dyz_Fxyz_S_S_b;
  abcd[445] = 4.0E0*I_ERI_G2x2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2z_Fxyz_S_S_b;
  abcd[446] = 4.0E0*I_ERI_Gx3y_Fxyz_S_S_ab;
  abcd[447] = 4.0E0*I_ERI_Gx2yz_Fxyz_S_S_ab;
  abcd[448] = 4.0E0*I_ERI_Gxy2z_Fxyz_S_S_ab;
  abcd[449] = 4.0E0*I_ERI_Gx3z_Fxyz_S_S_ab;
  abcd[450] = 4.0E0*I_ERI_G4x_F3y_S_S_ab-2.0E0*2*I_ERI_G4x_Py_S_S_a-2.0E0*3*I_ERI_D2x_F3y_S_S_b+3*2*I_ERI_D2x_Py_S_S;
  abcd[451] = 4.0E0*I_ERI_G3xy_F3y_S_S_ab-2.0E0*2*I_ERI_G3xy_Py_S_S_a-2.0E0*2*I_ERI_Dxy_F3y_S_S_b+2*2*I_ERI_Dxy_Py_S_S;
  abcd[452] = 4.0E0*I_ERI_G3xz_F3y_S_S_ab-2.0E0*2*I_ERI_G3xz_Py_S_S_a-2.0E0*2*I_ERI_Dxz_F3y_S_S_b+2*2*I_ERI_Dxz_Py_S_S;
  abcd[453] = 4.0E0*I_ERI_G2x2y_F3y_S_S_ab-2.0E0*2*I_ERI_G2x2y_Py_S_S_a-2.0E0*1*I_ERI_D2y_F3y_S_S_b+2*I_ERI_D2y_Py_S_S;
  abcd[454] = 4.0E0*I_ERI_G2xyz_F3y_S_S_ab-2.0E0*2*I_ERI_G2xyz_Py_S_S_a-2.0E0*1*I_ERI_Dyz_F3y_S_S_b+2*I_ERI_Dyz_Py_S_S;
  abcd[455] = 4.0E0*I_ERI_G2x2z_F3y_S_S_ab-2.0E0*2*I_ERI_G2x2z_Py_S_S_a-2.0E0*1*I_ERI_D2z_F3y_S_S_b+2*I_ERI_D2z_Py_S_S;
  abcd[456] = 4.0E0*I_ERI_Gx3y_F3y_S_S_ab-2.0E0*2*I_ERI_Gx3y_Py_S_S_a;
  abcd[457] = 4.0E0*I_ERI_Gx2yz_F3y_S_S_ab-2.0E0*2*I_ERI_Gx2yz_Py_S_S_a;
  abcd[458] = 4.0E0*I_ERI_Gxy2z_F3y_S_S_ab-2.0E0*2*I_ERI_Gxy2z_Py_S_S_a;
  abcd[459] = 4.0E0*I_ERI_Gx3z_F3y_S_S_ab-2.0E0*2*I_ERI_Gx3z_Py_S_S_a;
  abcd[460] = 4.0E0*I_ERI_G4x_F2yz_S_S_ab-2.0E0*1*I_ERI_G4x_Pz_S_S_a-2.0E0*3*I_ERI_D2x_F2yz_S_S_b+3*1*I_ERI_D2x_Pz_S_S;
  abcd[461] = 4.0E0*I_ERI_G3xy_F2yz_S_S_ab-2.0E0*1*I_ERI_G3xy_Pz_S_S_a-2.0E0*2*I_ERI_Dxy_F2yz_S_S_b+2*1*I_ERI_Dxy_Pz_S_S;
  abcd[462] = 4.0E0*I_ERI_G3xz_F2yz_S_S_ab-2.0E0*1*I_ERI_G3xz_Pz_S_S_a-2.0E0*2*I_ERI_Dxz_F2yz_S_S_b+2*1*I_ERI_Dxz_Pz_S_S;
  abcd[463] = 4.0E0*I_ERI_G2x2y_F2yz_S_S_ab-2.0E0*1*I_ERI_G2x2y_Pz_S_S_a-2.0E0*1*I_ERI_D2y_F2yz_S_S_b+1*I_ERI_D2y_Pz_S_S;
  abcd[464] = 4.0E0*I_ERI_G2xyz_F2yz_S_S_ab-2.0E0*1*I_ERI_G2xyz_Pz_S_S_a-2.0E0*1*I_ERI_Dyz_F2yz_S_S_b+1*I_ERI_Dyz_Pz_S_S;
  abcd[465] = 4.0E0*I_ERI_G2x2z_F2yz_S_S_ab-2.0E0*1*I_ERI_G2x2z_Pz_S_S_a-2.0E0*1*I_ERI_D2z_F2yz_S_S_b+1*I_ERI_D2z_Pz_S_S;
  abcd[466] = 4.0E0*I_ERI_Gx3y_F2yz_S_S_ab-2.0E0*1*I_ERI_Gx3y_Pz_S_S_a;
  abcd[467] = 4.0E0*I_ERI_Gx2yz_F2yz_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Pz_S_S_a;
  abcd[468] = 4.0E0*I_ERI_Gxy2z_F2yz_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Pz_S_S_a;
  abcd[469] = 4.0E0*I_ERI_Gx3z_F2yz_S_S_ab-2.0E0*1*I_ERI_Gx3z_Pz_S_S_a;
  abcd[470] = 4.0E0*I_ERI_G4x_Fy2z_S_S_ab-2.0E0*3*I_ERI_D2x_Fy2z_S_S_b;
  abcd[471] = 4.0E0*I_ERI_G3xy_Fy2z_S_S_ab-2.0E0*2*I_ERI_Dxy_Fy2z_S_S_b;
  abcd[472] = 4.0E0*I_ERI_G3xz_Fy2z_S_S_ab-2.0E0*2*I_ERI_Dxz_Fy2z_S_S_b;
  abcd[473] = 4.0E0*I_ERI_G2x2y_Fy2z_S_S_ab-2.0E0*1*I_ERI_D2y_Fy2z_S_S_b;
  abcd[474] = 4.0E0*I_ERI_G2xyz_Fy2z_S_S_ab-2.0E0*1*I_ERI_Dyz_Fy2z_S_S_b;
  abcd[475] = 4.0E0*I_ERI_G2x2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_D2z_Fy2z_S_S_b;
  abcd[476] = 4.0E0*I_ERI_Gx3y_Fy2z_S_S_ab;
  abcd[477] = 4.0E0*I_ERI_Gx2yz_Fy2z_S_S_ab;
  abcd[478] = 4.0E0*I_ERI_Gxy2z_Fy2z_S_S_ab;
  abcd[479] = 4.0E0*I_ERI_Gx3z_Fy2z_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_D_F_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  abcd[480] = 4.0E0*I_ERI_G4x_F2xz_S_S_ab-2.0E0*3*I_ERI_D2x_F2xz_S_S_b;
  abcd[481] = 4.0E0*I_ERI_G3xy_F2xz_S_S_ab-2.0E0*2*I_ERI_Dxy_F2xz_S_S_b;
  abcd[482] = 4.0E0*I_ERI_G3xz_F2xz_S_S_ab-2.0E0*2*I_ERI_Dxz_F2xz_S_S_b;
  abcd[483] = 4.0E0*I_ERI_G2x2y_F2xz_S_S_ab-2.0E0*1*I_ERI_D2y_F2xz_S_S_b;
  abcd[484] = 4.0E0*I_ERI_G2xyz_F2xz_S_S_ab-2.0E0*1*I_ERI_Dyz_F2xz_S_S_b;
  abcd[485] = 4.0E0*I_ERI_G2x2z_F2xz_S_S_ab-2.0E0*1*I_ERI_D2z_F2xz_S_S_b;
  abcd[486] = 4.0E0*I_ERI_Gx3y_F2xz_S_S_ab;
  abcd[487] = 4.0E0*I_ERI_Gx2yz_F2xz_S_S_ab;
  abcd[488] = 4.0E0*I_ERI_Gxy2z_F2xz_S_S_ab;
  abcd[489] = 4.0E0*I_ERI_Gx3z_F2xz_S_S_ab;
  abcd[490] = 4.0E0*I_ERI_G4x_Fxyz_S_S_ab-2.0E0*3*I_ERI_D2x_Fxyz_S_S_b;
  abcd[491] = 4.0E0*I_ERI_G3xy_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxy_Fxyz_S_S_b;
  abcd[492] = 4.0E0*I_ERI_G3xz_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxz_Fxyz_S_S_b;
  abcd[493] = 4.0E0*I_ERI_G2x2y_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2y_Fxyz_S_S_b;
  abcd[494] = 4.0E0*I_ERI_G2xyz_Fxyz_S_S_ab-2.0E0*1*I_ERI_Dyz_Fxyz_S_S_b;
  abcd[495] = 4.0E0*I_ERI_G2x2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2z_Fxyz_S_S_b;
  abcd[496] = 4.0E0*I_ERI_Gx3y_Fxyz_S_S_ab;
  abcd[497] = 4.0E0*I_ERI_Gx2yz_Fxyz_S_S_ab;
  abcd[498] = 4.0E0*I_ERI_Gxy2z_Fxyz_S_S_ab;
  abcd[499] = 4.0E0*I_ERI_Gx3z_Fxyz_S_S_ab;
  abcd[500] = 4.0E0*I_ERI_G4x_Fx2z_S_S_ab-2.0E0*1*I_ERI_G4x_Px_S_S_a-2.0E0*3*I_ERI_D2x_Fx2z_S_S_b+3*1*I_ERI_D2x_Px_S_S;
  abcd[501] = 4.0E0*I_ERI_G3xy_Fx2z_S_S_ab-2.0E0*1*I_ERI_G3xy_Px_S_S_a-2.0E0*2*I_ERI_Dxy_Fx2z_S_S_b+2*1*I_ERI_Dxy_Px_S_S;
  abcd[502] = 4.0E0*I_ERI_G3xz_Fx2z_S_S_ab-2.0E0*1*I_ERI_G3xz_Px_S_S_a-2.0E0*2*I_ERI_Dxz_Fx2z_S_S_b+2*1*I_ERI_Dxz_Px_S_S;
  abcd[503] = 4.0E0*I_ERI_G2x2y_Fx2z_S_S_ab-2.0E0*1*I_ERI_G2x2y_Px_S_S_a-2.0E0*1*I_ERI_D2y_Fx2z_S_S_b+1*I_ERI_D2y_Px_S_S;
  abcd[504] = 4.0E0*I_ERI_G2xyz_Fx2z_S_S_ab-2.0E0*1*I_ERI_G2xyz_Px_S_S_a-2.0E0*1*I_ERI_Dyz_Fx2z_S_S_b+1*I_ERI_Dyz_Px_S_S;
  abcd[505] = 4.0E0*I_ERI_G2x2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_G2x2z_Px_S_S_a-2.0E0*1*I_ERI_D2z_Fx2z_S_S_b+1*I_ERI_D2z_Px_S_S;
  abcd[506] = 4.0E0*I_ERI_Gx3y_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gx3y_Px_S_S_a;
  abcd[507] = 4.0E0*I_ERI_Gx2yz_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Px_S_S_a;
  abcd[508] = 4.0E0*I_ERI_Gxy2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Px_S_S_a;
  abcd[509] = 4.0E0*I_ERI_Gx3z_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gx3z_Px_S_S_a;
  abcd[510] = 4.0E0*I_ERI_G4x_F2yz_S_S_ab-2.0E0*3*I_ERI_D2x_F2yz_S_S_b;
  abcd[511] = 4.0E0*I_ERI_G3xy_F2yz_S_S_ab-2.0E0*2*I_ERI_Dxy_F2yz_S_S_b;
  abcd[512] = 4.0E0*I_ERI_G3xz_F2yz_S_S_ab-2.0E0*2*I_ERI_Dxz_F2yz_S_S_b;
  abcd[513] = 4.0E0*I_ERI_G2x2y_F2yz_S_S_ab-2.0E0*1*I_ERI_D2y_F2yz_S_S_b;
  abcd[514] = 4.0E0*I_ERI_G2xyz_F2yz_S_S_ab-2.0E0*1*I_ERI_Dyz_F2yz_S_S_b;
  abcd[515] = 4.0E0*I_ERI_G2x2z_F2yz_S_S_ab-2.0E0*1*I_ERI_D2z_F2yz_S_S_b;
  abcd[516] = 4.0E0*I_ERI_Gx3y_F2yz_S_S_ab;
  abcd[517] = 4.0E0*I_ERI_Gx2yz_F2yz_S_S_ab;
  abcd[518] = 4.0E0*I_ERI_Gxy2z_F2yz_S_S_ab;
  abcd[519] = 4.0E0*I_ERI_Gx3z_F2yz_S_S_ab;
  abcd[520] = 4.0E0*I_ERI_G4x_Fy2z_S_S_ab-2.0E0*1*I_ERI_G4x_Py_S_S_a-2.0E0*3*I_ERI_D2x_Fy2z_S_S_b+3*1*I_ERI_D2x_Py_S_S;
  abcd[521] = 4.0E0*I_ERI_G3xy_Fy2z_S_S_ab-2.0E0*1*I_ERI_G3xy_Py_S_S_a-2.0E0*2*I_ERI_Dxy_Fy2z_S_S_b+2*1*I_ERI_Dxy_Py_S_S;
  abcd[522] = 4.0E0*I_ERI_G3xz_Fy2z_S_S_ab-2.0E0*1*I_ERI_G3xz_Py_S_S_a-2.0E0*2*I_ERI_Dxz_Fy2z_S_S_b+2*1*I_ERI_Dxz_Py_S_S;
  abcd[523] = 4.0E0*I_ERI_G2x2y_Fy2z_S_S_ab-2.0E0*1*I_ERI_G2x2y_Py_S_S_a-2.0E0*1*I_ERI_D2y_Fy2z_S_S_b+1*I_ERI_D2y_Py_S_S;
  abcd[524] = 4.0E0*I_ERI_G2xyz_Fy2z_S_S_ab-2.0E0*1*I_ERI_G2xyz_Py_S_S_a-2.0E0*1*I_ERI_Dyz_Fy2z_S_S_b+1*I_ERI_Dyz_Py_S_S;
  abcd[525] = 4.0E0*I_ERI_G2x2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_G2x2z_Py_S_S_a-2.0E0*1*I_ERI_D2z_Fy2z_S_S_b+1*I_ERI_D2z_Py_S_S;
  abcd[526] = 4.0E0*I_ERI_Gx3y_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gx3y_Py_S_S_a;
  abcd[527] = 4.0E0*I_ERI_Gx2yz_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Py_S_S_a;
  abcd[528] = 4.0E0*I_ERI_Gxy2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Py_S_S_a;
  abcd[529] = 4.0E0*I_ERI_Gx3z_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gx3z_Py_S_S_a;
  abcd[530] = 4.0E0*I_ERI_G4x_F3z_S_S_ab-2.0E0*2*I_ERI_G4x_Pz_S_S_a-2.0E0*3*I_ERI_D2x_F3z_S_S_b+3*2*I_ERI_D2x_Pz_S_S;
  abcd[531] = 4.0E0*I_ERI_G3xy_F3z_S_S_ab-2.0E0*2*I_ERI_G3xy_Pz_S_S_a-2.0E0*2*I_ERI_Dxy_F3z_S_S_b+2*2*I_ERI_Dxy_Pz_S_S;
  abcd[532] = 4.0E0*I_ERI_G3xz_F3z_S_S_ab-2.0E0*2*I_ERI_G3xz_Pz_S_S_a-2.0E0*2*I_ERI_Dxz_F3z_S_S_b+2*2*I_ERI_Dxz_Pz_S_S;
  abcd[533] = 4.0E0*I_ERI_G2x2y_F3z_S_S_ab-2.0E0*2*I_ERI_G2x2y_Pz_S_S_a-2.0E0*1*I_ERI_D2y_F3z_S_S_b+2*I_ERI_D2y_Pz_S_S;
  abcd[534] = 4.0E0*I_ERI_G2xyz_F3z_S_S_ab-2.0E0*2*I_ERI_G2xyz_Pz_S_S_a-2.0E0*1*I_ERI_Dyz_F3z_S_S_b+2*I_ERI_Dyz_Pz_S_S;
  abcd[535] = 4.0E0*I_ERI_G2x2z_F3z_S_S_ab-2.0E0*2*I_ERI_G2x2z_Pz_S_S_a-2.0E0*1*I_ERI_D2z_F3z_S_S_b+2*I_ERI_D2z_Pz_S_S;
  abcd[536] = 4.0E0*I_ERI_Gx3y_F3z_S_S_ab-2.0E0*2*I_ERI_Gx3y_Pz_S_S_a;
  abcd[537] = 4.0E0*I_ERI_Gx2yz_F3z_S_S_ab-2.0E0*2*I_ERI_Gx2yz_Pz_S_S_a;
  abcd[538] = 4.0E0*I_ERI_Gxy2z_F3z_S_S_ab-2.0E0*2*I_ERI_Gxy2z_Pz_S_S_a;
  abcd[539] = 4.0E0*I_ERI_Gx3z_F3z_S_S_ab-2.0E0*2*I_ERI_Gx3z_Pz_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_D_F_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  abcd[540] = 4.0E0*I_ERI_G3xy_F3x_S_S_ab-2.0E0*2*I_ERI_G3xy_Px_S_S_a;
  abcd[541] = 4.0E0*I_ERI_G2x2y_F3x_S_S_ab-2.0E0*2*I_ERI_G2x2y_Px_S_S_a-2.0E0*1*I_ERI_D2x_F3x_S_S_b+2*I_ERI_D2x_Px_S_S;
  abcd[542] = 4.0E0*I_ERI_G2xyz_F3x_S_S_ab-2.0E0*2*I_ERI_G2xyz_Px_S_S_a;
  abcd[543] = 4.0E0*I_ERI_Gx3y_F3x_S_S_ab-2.0E0*2*I_ERI_Gx3y_Px_S_S_a-2.0E0*2*I_ERI_Dxy_F3x_S_S_b+2*2*I_ERI_Dxy_Px_S_S;
  abcd[544] = 4.0E0*I_ERI_Gx2yz_F3x_S_S_ab-2.0E0*2*I_ERI_Gx2yz_Px_S_S_a-2.0E0*1*I_ERI_Dxz_F3x_S_S_b+2*I_ERI_Dxz_Px_S_S;
  abcd[545] = 4.0E0*I_ERI_Gxy2z_F3x_S_S_ab-2.0E0*2*I_ERI_Gxy2z_Px_S_S_a;
  abcd[546] = 4.0E0*I_ERI_G4y_F3x_S_S_ab-2.0E0*2*I_ERI_G4y_Px_S_S_a-2.0E0*3*I_ERI_D2y_F3x_S_S_b+3*2*I_ERI_D2y_Px_S_S;
  abcd[547] = 4.0E0*I_ERI_G3yz_F3x_S_S_ab-2.0E0*2*I_ERI_G3yz_Px_S_S_a-2.0E0*2*I_ERI_Dyz_F3x_S_S_b+2*2*I_ERI_Dyz_Px_S_S;
  abcd[548] = 4.0E0*I_ERI_G2y2z_F3x_S_S_ab-2.0E0*2*I_ERI_G2y2z_Px_S_S_a-2.0E0*1*I_ERI_D2z_F3x_S_S_b+2*I_ERI_D2z_Px_S_S;
  abcd[549] = 4.0E0*I_ERI_Gy3z_F3x_S_S_ab-2.0E0*2*I_ERI_Gy3z_Px_S_S_a;
  abcd[550] = 4.0E0*I_ERI_G3xy_F2xy_S_S_ab-2.0E0*1*I_ERI_G3xy_Py_S_S_a;
  abcd[551] = 4.0E0*I_ERI_G2x2y_F2xy_S_S_ab-2.0E0*1*I_ERI_G2x2y_Py_S_S_a-2.0E0*1*I_ERI_D2x_F2xy_S_S_b+1*I_ERI_D2x_Py_S_S;
  abcd[552] = 4.0E0*I_ERI_G2xyz_F2xy_S_S_ab-2.0E0*1*I_ERI_G2xyz_Py_S_S_a;
  abcd[553] = 4.0E0*I_ERI_Gx3y_F2xy_S_S_ab-2.0E0*1*I_ERI_Gx3y_Py_S_S_a-2.0E0*2*I_ERI_Dxy_F2xy_S_S_b+2*1*I_ERI_Dxy_Py_S_S;
  abcd[554] = 4.0E0*I_ERI_Gx2yz_F2xy_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Py_S_S_a-2.0E0*1*I_ERI_Dxz_F2xy_S_S_b+1*I_ERI_Dxz_Py_S_S;
  abcd[555] = 4.0E0*I_ERI_Gxy2z_F2xy_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Py_S_S_a;
  abcd[556] = 4.0E0*I_ERI_G4y_F2xy_S_S_ab-2.0E0*1*I_ERI_G4y_Py_S_S_a-2.0E0*3*I_ERI_D2y_F2xy_S_S_b+3*1*I_ERI_D2y_Py_S_S;
  abcd[557] = 4.0E0*I_ERI_G3yz_F2xy_S_S_ab-2.0E0*1*I_ERI_G3yz_Py_S_S_a-2.0E0*2*I_ERI_Dyz_F2xy_S_S_b+2*1*I_ERI_Dyz_Py_S_S;
  abcd[558] = 4.0E0*I_ERI_G2y2z_F2xy_S_S_ab-2.0E0*1*I_ERI_G2y2z_Py_S_S_a-2.0E0*1*I_ERI_D2z_F2xy_S_S_b+1*I_ERI_D2z_Py_S_S;
  abcd[559] = 4.0E0*I_ERI_Gy3z_F2xy_S_S_ab-2.0E0*1*I_ERI_Gy3z_Py_S_S_a;
  abcd[560] = 4.0E0*I_ERI_G3xy_F2xz_S_S_ab-2.0E0*1*I_ERI_G3xy_Pz_S_S_a;
  abcd[561] = 4.0E0*I_ERI_G2x2y_F2xz_S_S_ab-2.0E0*1*I_ERI_G2x2y_Pz_S_S_a-2.0E0*1*I_ERI_D2x_F2xz_S_S_b+1*I_ERI_D2x_Pz_S_S;
  abcd[562] = 4.0E0*I_ERI_G2xyz_F2xz_S_S_ab-2.0E0*1*I_ERI_G2xyz_Pz_S_S_a;
  abcd[563] = 4.0E0*I_ERI_Gx3y_F2xz_S_S_ab-2.0E0*1*I_ERI_Gx3y_Pz_S_S_a-2.0E0*2*I_ERI_Dxy_F2xz_S_S_b+2*1*I_ERI_Dxy_Pz_S_S;
  abcd[564] = 4.0E0*I_ERI_Gx2yz_F2xz_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Pz_S_S_a-2.0E0*1*I_ERI_Dxz_F2xz_S_S_b+1*I_ERI_Dxz_Pz_S_S;
  abcd[565] = 4.0E0*I_ERI_Gxy2z_F2xz_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Pz_S_S_a;
  abcd[566] = 4.0E0*I_ERI_G4y_F2xz_S_S_ab-2.0E0*1*I_ERI_G4y_Pz_S_S_a-2.0E0*3*I_ERI_D2y_F2xz_S_S_b+3*1*I_ERI_D2y_Pz_S_S;
  abcd[567] = 4.0E0*I_ERI_G3yz_F2xz_S_S_ab-2.0E0*1*I_ERI_G3yz_Pz_S_S_a-2.0E0*2*I_ERI_Dyz_F2xz_S_S_b+2*1*I_ERI_Dyz_Pz_S_S;
  abcd[568] = 4.0E0*I_ERI_G2y2z_F2xz_S_S_ab-2.0E0*1*I_ERI_G2y2z_Pz_S_S_a-2.0E0*1*I_ERI_D2z_F2xz_S_S_b+1*I_ERI_D2z_Pz_S_S;
  abcd[569] = 4.0E0*I_ERI_Gy3z_F2xz_S_S_ab-2.0E0*1*I_ERI_Gy3z_Pz_S_S_a;
  abcd[570] = 4.0E0*I_ERI_G3xy_Fx2y_S_S_ab;
  abcd[571] = 4.0E0*I_ERI_G2x2y_Fx2y_S_S_ab-2.0E0*1*I_ERI_D2x_Fx2y_S_S_b;
  abcd[572] = 4.0E0*I_ERI_G2xyz_Fx2y_S_S_ab;
  abcd[573] = 4.0E0*I_ERI_Gx3y_Fx2y_S_S_ab-2.0E0*2*I_ERI_Dxy_Fx2y_S_S_b;
  abcd[574] = 4.0E0*I_ERI_Gx2yz_Fx2y_S_S_ab-2.0E0*1*I_ERI_Dxz_Fx2y_S_S_b;
  abcd[575] = 4.0E0*I_ERI_Gxy2z_Fx2y_S_S_ab;
  abcd[576] = 4.0E0*I_ERI_G4y_Fx2y_S_S_ab-2.0E0*3*I_ERI_D2y_Fx2y_S_S_b;
  abcd[577] = 4.0E0*I_ERI_G3yz_Fx2y_S_S_ab-2.0E0*2*I_ERI_Dyz_Fx2y_S_S_b;
  abcd[578] = 4.0E0*I_ERI_G2y2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_D2z_Fx2y_S_S_b;
  abcd[579] = 4.0E0*I_ERI_Gy3z_Fx2y_S_S_ab;
  abcd[580] = 4.0E0*I_ERI_G3xy_Fxyz_S_S_ab;
  abcd[581] = 4.0E0*I_ERI_G2x2y_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2x_Fxyz_S_S_b;
  abcd[582] = 4.0E0*I_ERI_G2xyz_Fxyz_S_S_ab;
  abcd[583] = 4.0E0*I_ERI_Gx3y_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxy_Fxyz_S_S_b;
  abcd[584] = 4.0E0*I_ERI_Gx2yz_Fxyz_S_S_ab-2.0E0*1*I_ERI_Dxz_Fxyz_S_S_b;
  abcd[585] = 4.0E0*I_ERI_Gxy2z_Fxyz_S_S_ab;
  abcd[586] = 4.0E0*I_ERI_G4y_Fxyz_S_S_ab-2.0E0*3*I_ERI_D2y_Fxyz_S_S_b;
  abcd[587] = 4.0E0*I_ERI_G3yz_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dyz_Fxyz_S_S_b;
  abcd[588] = 4.0E0*I_ERI_G2y2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2z_Fxyz_S_S_b;
  abcd[589] = 4.0E0*I_ERI_Gy3z_Fxyz_S_S_ab;
  abcd[590] = 4.0E0*I_ERI_G3xy_Fx2z_S_S_ab;
  abcd[591] = 4.0E0*I_ERI_G2x2y_Fx2z_S_S_ab-2.0E0*1*I_ERI_D2x_Fx2z_S_S_b;
  abcd[592] = 4.0E0*I_ERI_G2xyz_Fx2z_S_S_ab;
  abcd[593] = 4.0E0*I_ERI_Gx3y_Fx2z_S_S_ab-2.0E0*2*I_ERI_Dxy_Fx2z_S_S_b;
  abcd[594] = 4.0E0*I_ERI_Gx2yz_Fx2z_S_S_ab-2.0E0*1*I_ERI_Dxz_Fx2z_S_S_b;
  abcd[595] = 4.0E0*I_ERI_Gxy2z_Fx2z_S_S_ab;
  abcd[596] = 4.0E0*I_ERI_G4y_Fx2z_S_S_ab-2.0E0*3*I_ERI_D2y_Fx2z_S_S_b;
  abcd[597] = 4.0E0*I_ERI_G3yz_Fx2z_S_S_ab-2.0E0*2*I_ERI_Dyz_Fx2z_S_S_b;
  abcd[598] = 4.0E0*I_ERI_G2y2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_D2z_Fx2z_S_S_b;
  abcd[599] = 4.0E0*I_ERI_Gy3z_Fx2z_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_D_F_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  abcd[600] = 4.0E0*I_ERI_G3xy_F2xy_S_S_ab;
  abcd[601] = 4.0E0*I_ERI_G2x2y_F2xy_S_S_ab-2.0E0*1*I_ERI_D2x_F2xy_S_S_b;
  abcd[602] = 4.0E0*I_ERI_G2xyz_F2xy_S_S_ab;
  abcd[603] = 4.0E0*I_ERI_Gx3y_F2xy_S_S_ab-2.0E0*2*I_ERI_Dxy_F2xy_S_S_b;
  abcd[604] = 4.0E0*I_ERI_Gx2yz_F2xy_S_S_ab-2.0E0*1*I_ERI_Dxz_F2xy_S_S_b;
  abcd[605] = 4.0E0*I_ERI_Gxy2z_F2xy_S_S_ab;
  abcd[606] = 4.0E0*I_ERI_G4y_F2xy_S_S_ab-2.0E0*3*I_ERI_D2y_F2xy_S_S_b;
  abcd[607] = 4.0E0*I_ERI_G3yz_F2xy_S_S_ab-2.0E0*2*I_ERI_Dyz_F2xy_S_S_b;
  abcd[608] = 4.0E0*I_ERI_G2y2z_F2xy_S_S_ab-2.0E0*1*I_ERI_D2z_F2xy_S_S_b;
  abcd[609] = 4.0E0*I_ERI_Gy3z_F2xy_S_S_ab;
  abcd[610] = 4.0E0*I_ERI_G3xy_Fx2y_S_S_ab-2.0E0*1*I_ERI_G3xy_Px_S_S_a;
  abcd[611] = 4.0E0*I_ERI_G2x2y_Fx2y_S_S_ab-2.0E0*1*I_ERI_G2x2y_Px_S_S_a-2.0E0*1*I_ERI_D2x_Fx2y_S_S_b+1*I_ERI_D2x_Px_S_S;
  abcd[612] = 4.0E0*I_ERI_G2xyz_Fx2y_S_S_ab-2.0E0*1*I_ERI_G2xyz_Px_S_S_a;
  abcd[613] = 4.0E0*I_ERI_Gx3y_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gx3y_Px_S_S_a-2.0E0*2*I_ERI_Dxy_Fx2y_S_S_b+2*1*I_ERI_Dxy_Px_S_S;
  abcd[614] = 4.0E0*I_ERI_Gx2yz_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Px_S_S_a-2.0E0*1*I_ERI_Dxz_Fx2y_S_S_b+1*I_ERI_Dxz_Px_S_S;
  abcd[615] = 4.0E0*I_ERI_Gxy2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Px_S_S_a;
  abcd[616] = 4.0E0*I_ERI_G4y_Fx2y_S_S_ab-2.0E0*1*I_ERI_G4y_Px_S_S_a-2.0E0*3*I_ERI_D2y_Fx2y_S_S_b+3*1*I_ERI_D2y_Px_S_S;
  abcd[617] = 4.0E0*I_ERI_G3yz_Fx2y_S_S_ab-2.0E0*1*I_ERI_G3yz_Px_S_S_a-2.0E0*2*I_ERI_Dyz_Fx2y_S_S_b+2*1*I_ERI_Dyz_Px_S_S;
  abcd[618] = 4.0E0*I_ERI_G2y2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_G2y2z_Px_S_S_a-2.0E0*1*I_ERI_D2z_Fx2y_S_S_b+1*I_ERI_D2z_Px_S_S;
  abcd[619] = 4.0E0*I_ERI_Gy3z_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gy3z_Px_S_S_a;
  abcd[620] = 4.0E0*I_ERI_G3xy_Fxyz_S_S_ab;
  abcd[621] = 4.0E0*I_ERI_G2x2y_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2x_Fxyz_S_S_b;
  abcd[622] = 4.0E0*I_ERI_G2xyz_Fxyz_S_S_ab;
  abcd[623] = 4.0E0*I_ERI_Gx3y_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxy_Fxyz_S_S_b;
  abcd[624] = 4.0E0*I_ERI_Gx2yz_Fxyz_S_S_ab-2.0E0*1*I_ERI_Dxz_Fxyz_S_S_b;
  abcd[625] = 4.0E0*I_ERI_Gxy2z_Fxyz_S_S_ab;
  abcd[626] = 4.0E0*I_ERI_G4y_Fxyz_S_S_ab-2.0E0*3*I_ERI_D2y_Fxyz_S_S_b;
  abcd[627] = 4.0E0*I_ERI_G3yz_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dyz_Fxyz_S_S_b;
  abcd[628] = 4.0E0*I_ERI_G2y2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2z_Fxyz_S_S_b;
  abcd[629] = 4.0E0*I_ERI_Gy3z_Fxyz_S_S_ab;
  abcd[630] = 4.0E0*I_ERI_G3xy_F3y_S_S_ab-2.0E0*2*I_ERI_G3xy_Py_S_S_a;
  abcd[631] = 4.0E0*I_ERI_G2x2y_F3y_S_S_ab-2.0E0*2*I_ERI_G2x2y_Py_S_S_a-2.0E0*1*I_ERI_D2x_F3y_S_S_b+2*I_ERI_D2x_Py_S_S;
  abcd[632] = 4.0E0*I_ERI_G2xyz_F3y_S_S_ab-2.0E0*2*I_ERI_G2xyz_Py_S_S_a;
  abcd[633] = 4.0E0*I_ERI_Gx3y_F3y_S_S_ab-2.0E0*2*I_ERI_Gx3y_Py_S_S_a-2.0E0*2*I_ERI_Dxy_F3y_S_S_b+2*2*I_ERI_Dxy_Py_S_S;
  abcd[634] = 4.0E0*I_ERI_Gx2yz_F3y_S_S_ab-2.0E0*2*I_ERI_Gx2yz_Py_S_S_a-2.0E0*1*I_ERI_Dxz_F3y_S_S_b+2*I_ERI_Dxz_Py_S_S;
  abcd[635] = 4.0E0*I_ERI_Gxy2z_F3y_S_S_ab-2.0E0*2*I_ERI_Gxy2z_Py_S_S_a;
  abcd[636] = 4.0E0*I_ERI_G4y_F3y_S_S_ab-2.0E0*2*I_ERI_G4y_Py_S_S_a-2.0E0*3*I_ERI_D2y_F3y_S_S_b+3*2*I_ERI_D2y_Py_S_S;
  abcd[637] = 4.0E0*I_ERI_G3yz_F3y_S_S_ab-2.0E0*2*I_ERI_G3yz_Py_S_S_a-2.0E0*2*I_ERI_Dyz_F3y_S_S_b+2*2*I_ERI_Dyz_Py_S_S;
  abcd[638] = 4.0E0*I_ERI_G2y2z_F3y_S_S_ab-2.0E0*2*I_ERI_G2y2z_Py_S_S_a-2.0E0*1*I_ERI_D2z_F3y_S_S_b+2*I_ERI_D2z_Py_S_S;
  abcd[639] = 4.0E0*I_ERI_Gy3z_F3y_S_S_ab-2.0E0*2*I_ERI_Gy3z_Py_S_S_a;
  abcd[640] = 4.0E0*I_ERI_G3xy_F2yz_S_S_ab-2.0E0*1*I_ERI_G3xy_Pz_S_S_a;
  abcd[641] = 4.0E0*I_ERI_G2x2y_F2yz_S_S_ab-2.0E0*1*I_ERI_G2x2y_Pz_S_S_a-2.0E0*1*I_ERI_D2x_F2yz_S_S_b+1*I_ERI_D2x_Pz_S_S;
  abcd[642] = 4.0E0*I_ERI_G2xyz_F2yz_S_S_ab-2.0E0*1*I_ERI_G2xyz_Pz_S_S_a;
  abcd[643] = 4.0E0*I_ERI_Gx3y_F2yz_S_S_ab-2.0E0*1*I_ERI_Gx3y_Pz_S_S_a-2.0E0*2*I_ERI_Dxy_F2yz_S_S_b+2*1*I_ERI_Dxy_Pz_S_S;
  abcd[644] = 4.0E0*I_ERI_Gx2yz_F2yz_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Pz_S_S_a-2.0E0*1*I_ERI_Dxz_F2yz_S_S_b+1*I_ERI_Dxz_Pz_S_S;
  abcd[645] = 4.0E0*I_ERI_Gxy2z_F2yz_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Pz_S_S_a;
  abcd[646] = 4.0E0*I_ERI_G4y_F2yz_S_S_ab-2.0E0*1*I_ERI_G4y_Pz_S_S_a-2.0E0*3*I_ERI_D2y_F2yz_S_S_b+3*1*I_ERI_D2y_Pz_S_S;
  abcd[647] = 4.0E0*I_ERI_G3yz_F2yz_S_S_ab-2.0E0*1*I_ERI_G3yz_Pz_S_S_a-2.0E0*2*I_ERI_Dyz_F2yz_S_S_b+2*1*I_ERI_Dyz_Pz_S_S;
  abcd[648] = 4.0E0*I_ERI_G2y2z_F2yz_S_S_ab-2.0E0*1*I_ERI_G2y2z_Pz_S_S_a-2.0E0*1*I_ERI_D2z_F2yz_S_S_b+1*I_ERI_D2z_Pz_S_S;
  abcd[649] = 4.0E0*I_ERI_Gy3z_F2yz_S_S_ab-2.0E0*1*I_ERI_Gy3z_Pz_S_S_a;
  abcd[650] = 4.0E0*I_ERI_G3xy_Fy2z_S_S_ab;
  abcd[651] = 4.0E0*I_ERI_G2x2y_Fy2z_S_S_ab-2.0E0*1*I_ERI_D2x_Fy2z_S_S_b;
  abcd[652] = 4.0E0*I_ERI_G2xyz_Fy2z_S_S_ab;
  abcd[653] = 4.0E0*I_ERI_Gx3y_Fy2z_S_S_ab-2.0E0*2*I_ERI_Dxy_Fy2z_S_S_b;
  abcd[654] = 4.0E0*I_ERI_Gx2yz_Fy2z_S_S_ab-2.0E0*1*I_ERI_Dxz_Fy2z_S_S_b;
  abcd[655] = 4.0E0*I_ERI_Gxy2z_Fy2z_S_S_ab;
  abcd[656] = 4.0E0*I_ERI_G4y_Fy2z_S_S_ab-2.0E0*3*I_ERI_D2y_Fy2z_S_S_b;
  abcd[657] = 4.0E0*I_ERI_G3yz_Fy2z_S_S_ab-2.0E0*2*I_ERI_Dyz_Fy2z_S_S_b;
  abcd[658] = 4.0E0*I_ERI_G2y2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_D2z_Fy2z_S_S_b;
  abcd[659] = 4.0E0*I_ERI_Gy3z_Fy2z_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_D_F_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  abcd[660] = 4.0E0*I_ERI_G3xy_F2xz_S_S_ab;
  abcd[661] = 4.0E0*I_ERI_G2x2y_F2xz_S_S_ab-2.0E0*1*I_ERI_D2x_F2xz_S_S_b;
  abcd[662] = 4.0E0*I_ERI_G2xyz_F2xz_S_S_ab;
  abcd[663] = 4.0E0*I_ERI_Gx3y_F2xz_S_S_ab-2.0E0*2*I_ERI_Dxy_F2xz_S_S_b;
  abcd[664] = 4.0E0*I_ERI_Gx2yz_F2xz_S_S_ab-2.0E0*1*I_ERI_Dxz_F2xz_S_S_b;
  abcd[665] = 4.0E0*I_ERI_Gxy2z_F2xz_S_S_ab;
  abcd[666] = 4.0E0*I_ERI_G4y_F2xz_S_S_ab-2.0E0*3*I_ERI_D2y_F2xz_S_S_b;
  abcd[667] = 4.0E0*I_ERI_G3yz_F2xz_S_S_ab-2.0E0*2*I_ERI_Dyz_F2xz_S_S_b;
  abcd[668] = 4.0E0*I_ERI_G2y2z_F2xz_S_S_ab-2.0E0*1*I_ERI_D2z_F2xz_S_S_b;
  abcd[669] = 4.0E0*I_ERI_Gy3z_F2xz_S_S_ab;
  abcd[670] = 4.0E0*I_ERI_G3xy_Fxyz_S_S_ab;
  abcd[671] = 4.0E0*I_ERI_G2x2y_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2x_Fxyz_S_S_b;
  abcd[672] = 4.0E0*I_ERI_G2xyz_Fxyz_S_S_ab;
  abcd[673] = 4.0E0*I_ERI_Gx3y_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxy_Fxyz_S_S_b;
  abcd[674] = 4.0E0*I_ERI_Gx2yz_Fxyz_S_S_ab-2.0E0*1*I_ERI_Dxz_Fxyz_S_S_b;
  abcd[675] = 4.0E0*I_ERI_Gxy2z_Fxyz_S_S_ab;
  abcd[676] = 4.0E0*I_ERI_G4y_Fxyz_S_S_ab-2.0E0*3*I_ERI_D2y_Fxyz_S_S_b;
  abcd[677] = 4.0E0*I_ERI_G3yz_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dyz_Fxyz_S_S_b;
  abcd[678] = 4.0E0*I_ERI_G2y2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2z_Fxyz_S_S_b;
  abcd[679] = 4.0E0*I_ERI_Gy3z_Fxyz_S_S_ab;
  abcd[680] = 4.0E0*I_ERI_G3xy_Fx2z_S_S_ab-2.0E0*1*I_ERI_G3xy_Px_S_S_a;
  abcd[681] = 4.0E0*I_ERI_G2x2y_Fx2z_S_S_ab-2.0E0*1*I_ERI_G2x2y_Px_S_S_a-2.0E0*1*I_ERI_D2x_Fx2z_S_S_b+1*I_ERI_D2x_Px_S_S;
  abcd[682] = 4.0E0*I_ERI_G2xyz_Fx2z_S_S_ab-2.0E0*1*I_ERI_G2xyz_Px_S_S_a;
  abcd[683] = 4.0E0*I_ERI_Gx3y_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gx3y_Px_S_S_a-2.0E0*2*I_ERI_Dxy_Fx2z_S_S_b+2*1*I_ERI_Dxy_Px_S_S;
  abcd[684] = 4.0E0*I_ERI_Gx2yz_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Px_S_S_a-2.0E0*1*I_ERI_Dxz_Fx2z_S_S_b+1*I_ERI_Dxz_Px_S_S;
  abcd[685] = 4.0E0*I_ERI_Gxy2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Px_S_S_a;
  abcd[686] = 4.0E0*I_ERI_G4y_Fx2z_S_S_ab-2.0E0*1*I_ERI_G4y_Px_S_S_a-2.0E0*3*I_ERI_D2y_Fx2z_S_S_b+3*1*I_ERI_D2y_Px_S_S;
  abcd[687] = 4.0E0*I_ERI_G3yz_Fx2z_S_S_ab-2.0E0*1*I_ERI_G3yz_Px_S_S_a-2.0E0*2*I_ERI_Dyz_Fx2z_S_S_b+2*1*I_ERI_Dyz_Px_S_S;
  abcd[688] = 4.0E0*I_ERI_G2y2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_G2y2z_Px_S_S_a-2.0E0*1*I_ERI_D2z_Fx2z_S_S_b+1*I_ERI_D2z_Px_S_S;
  abcd[689] = 4.0E0*I_ERI_Gy3z_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gy3z_Px_S_S_a;
  abcd[690] = 4.0E0*I_ERI_G3xy_F2yz_S_S_ab;
  abcd[691] = 4.0E0*I_ERI_G2x2y_F2yz_S_S_ab-2.0E0*1*I_ERI_D2x_F2yz_S_S_b;
  abcd[692] = 4.0E0*I_ERI_G2xyz_F2yz_S_S_ab;
  abcd[693] = 4.0E0*I_ERI_Gx3y_F2yz_S_S_ab-2.0E0*2*I_ERI_Dxy_F2yz_S_S_b;
  abcd[694] = 4.0E0*I_ERI_Gx2yz_F2yz_S_S_ab-2.0E0*1*I_ERI_Dxz_F2yz_S_S_b;
  abcd[695] = 4.0E0*I_ERI_Gxy2z_F2yz_S_S_ab;
  abcd[696] = 4.0E0*I_ERI_G4y_F2yz_S_S_ab-2.0E0*3*I_ERI_D2y_F2yz_S_S_b;
  abcd[697] = 4.0E0*I_ERI_G3yz_F2yz_S_S_ab-2.0E0*2*I_ERI_Dyz_F2yz_S_S_b;
  abcd[698] = 4.0E0*I_ERI_G2y2z_F2yz_S_S_ab-2.0E0*1*I_ERI_D2z_F2yz_S_S_b;
  abcd[699] = 4.0E0*I_ERI_Gy3z_F2yz_S_S_ab;
  abcd[700] = 4.0E0*I_ERI_G3xy_Fy2z_S_S_ab-2.0E0*1*I_ERI_G3xy_Py_S_S_a;
  abcd[701] = 4.0E0*I_ERI_G2x2y_Fy2z_S_S_ab-2.0E0*1*I_ERI_G2x2y_Py_S_S_a-2.0E0*1*I_ERI_D2x_Fy2z_S_S_b+1*I_ERI_D2x_Py_S_S;
  abcd[702] = 4.0E0*I_ERI_G2xyz_Fy2z_S_S_ab-2.0E0*1*I_ERI_G2xyz_Py_S_S_a;
  abcd[703] = 4.0E0*I_ERI_Gx3y_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gx3y_Py_S_S_a-2.0E0*2*I_ERI_Dxy_Fy2z_S_S_b+2*1*I_ERI_Dxy_Py_S_S;
  abcd[704] = 4.0E0*I_ERI_Gx2yz_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Py_S_S_a-2.0E0*1*I_ERI_Dxz_Fy2z_S_S_b+1*I_ERI_Dxz_Py_S_S;
  abcd[705] = 4.0E0*I_ERI_Gxy2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Py_S_S_a;
  abcd[706] = 4.0E0*I_ERI_G4y_Fy2z_S_S_ab-2.0E0*1*I_ERI_G4y_Py_S_S_a-2.0E0*3*I_ERI_D2y_Fy2z_S_S_b+3*1*I_ERI_D2y_Py_S_S;
  abcd[707] = 4.0E0*I_ERI_G3yz_Fy2z_S_S_ab-2.0E0*1*I_ERI_G3yz_Py_S_S_a-2.0E0*2*I_ERI_Dyz_Fy2z_S_S_b+2*1*I_ERI_Dyz_Py_S_S;
  abcd[708] = 4.0E0*I_ERI_G2y2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_G2y2z_Py_S_S_a-2.0E0*1*I_ERI_D2z_Fy2z_S_S_b+1*I_ERI_D2z_Py_S_S;
  abcd[709] = 4.0E0*I_ERI_Gy3z_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gy3z_Py_S_S_a;
  abcd[710] = 4.0E0*I_ERI_G3xy_F3z_S_S_ab-2.0E0*2*I_ERI_G3xy_Pz_S_S_a;
  abcd[711] = 4.0E0*I_ERI_G2x2y_F3z_S_S_ab-2.0E0*2*I_ERI_G2x2y_Pz_S_S_a-2.0E0*1*I_ERI_D2x_F3z_S_S_b+2*I_ERI_D2x_Pz_S_S;
  abcd[712] = 4.0E0*I_ERI_G2xyz_F3z_S_S_ab-2.0E0*2*I_ERI_G2xyz_Pz_S_S_a;
  abcd[713] = 4.0E0*I_ERI_Gx3y_F3z_S_S_ab-2.0E0*2*I_ERI_Gx3y_Pz_S_S_a-2.0E0*2*I_ERI_Dxy_F3z_S_S_b+2*2*I_ERI_Dxy_Pz_S_S;
  abcd[714] = 4.0E0*I_ERI_Gx2yz_F3z_S_S_ab-2.0E0*2*I_ERI_Gx2yz_Pz_S_S_a-2.0E0*1*I_ERI_Dxz_F3z_S_S_b+2*I_ERI_Dxz_Pz_S_S;
  abcd[715] = 4.0E0*I_ERI_Gxy2z_F3z_S_S_ab-2.0E0*2*I_ERI_Gxy2z_Pz_S_S_a;
  abcd[716] = 4.0E0*I_ERI_G4y_F3z_S_S_ab-2.0E0*2*I_ERI_G4y_Pz_S_S_a-2.0E0*3*I_ERI_D2y_F3z_S_S_b+3*2*I_ERI_D2y_Pz_S_S;
  abcd[717] = 4.0E0*I_ERI_G3yz_F3z_S_S_ab-2.0E0*2*I_ERI_G3yz_Pz_S_S_a-2.0E0*2*I_ERI_Dyz_F3z_S_S_b+2*2*I_ERI_Dyz_Pz_S_S;
  abcd[718] = 4.0E0*I_ERI_G2y2z_F3z_S_S_ab-2.0E0*2*I_ERI_G2y2z_Pz_S_S_a-2.0E0*1*I_ERI_D2z_F3z_S_S_b+2*I_ERI_D2z_Pz_S_S;
  abcd[719] = 4.0E0*I_ERI_Gy3z_F3z_S_S_ab-2.0E0*2*I_ERI_Gy3z_Pz_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_D_F_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  abcd[720] = 4.0E0*I_ERI_G3xz_F3x_S_S_ab-2.0E0*2*I_ERI_G3xz_Px_S_S_a;
  abcd[721] = 4.0E0*I_ERI_G2xyz_F3x_S_S_ab-2.0E0*2*I_ERI_G2xyz_Px_S_S_a;
  abcd[722] = 4.0E0*I_ERI_G2x2z_F3x_S_S_ab-2.0E0*2*I_ERI_G2x2z_Px_S_S_a-2.0E0*1*I_ERI_D2x_F3x_S_S_b+2*I_ERI_D2x_Px_S_S;
  abcd[723] = 4.0E0*I_ERI_Gx2yz_F3x_S_S_ab-2.0E0*2*I_ERI_Gx2yz_Px_S_S_a;
  abcd[724] = 4.0E0*I_ERI_Gxy2z_F3x_S_S_ab-2.0E0*2*I_ERI_Gxy2z_Px_S_S_a-2.0E0*1*I_ERI_Dxy_F3x_S_S_b+2*I_ERI_Dxy_Px_S_S;
  abcd[725] = 4.0E0*I_ERI_Gx3z_F3x_S_S_ab-2.0E0*2*I_ERI_Gx3z_Px_S_S_a-2.0E0*2*I_ERI_Dxz_F3x_S_S_b+2*2*I_ERI_Dxz_Px_S_S;
  abcd[726] = 4.0E0*I_ERI_G3yz_F3x_S_S_ab-2.0E0*2*I_ERI_G3yz_Px_S_S_a;
  abcd[727] = 4.0E0*I_ERI_G2y2z_F3x_S_S_ab-2.0E0*2*I_ERI_G2y2z_Px_S_S_a-2.0E0*1*I_ERI_D2y_F3x_S_S_b+2*I_ERI_D2y_Px_S_S;
  abcd[728] = 4.0E0*I_ERI_Gy3z_F3x_S_S_ab-2.0E0*2*I_ERI_Gy3z_Px_S_S_a-2.0E0*2*I_ERI_Dyz_F3x_S_S_b+2*2*I_ERI_Dyz_Px_S_S;
  abcd[729] = 4.0E0*I_ERI_G4z_F3x_S_S_ab-2.0E0*2*I_ERI_G4z_Px_S_S_a-2.0E0*3*I_ERI_D2z_F3x_S_S_b+3*2*I_ERI_D2z_Px_S_S;
  abcd[730] = 4.0E0*I_ERI_G3xz_F2xy_S_S_ab-2.0E0*1*I_ERI_G3xz_Py_S_S_a;
  abcd[731] = 4.0E0*I_ERI_G2xyz_F2xy_S_S_ab-2.0E0*1*I_ERI_G2xyz_Py_S_S_a;
  abcd[732] = 4.0E0*I_ERI_G2x2z_F2xy_S_S_ab-2.0E0*1*I_ERI_G2x2z_Py_S_S_a-2.0E0*1*I_ERI_D2x_F2xy_S_S_b+1*I_ERI_D2x_Py_S_S;
  abcd[733] = 4.0E0*I_ERI_Gx2yz_F2xy_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Py_S_S_a;
  abcd[734] = 4.0E0*I_ERI_Gxy2z_F2xy_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Py_S_S_a-2.0E0*1*I_ERI_Dxy_F2xy_S_S_b+1*I_ERI_Dxy_Py_S_S;
  abcd[735] = 4.0E0*I_ERI_Gx3z_F2xy_S_S_ab-2.0E0*1*I_ERI_Gx3z_Py_S_S_a-2.0E0*2*I_ERI_Dxz_F2xy_S_S_b+2*1*I_ERI_Dxz_Py_S_S;
  abcd[736] = 4.0E0*I_ERI_G3yz_F2xy_S_S_ab-2.0E0*1*I_ERI_G3yz_Py_S_S_a;
  abcd[737] = 4.0E0*I_ERI_G2y2z_F2xy_S_S_ab-2.0E0*1*I_ERI_G2y2z_Py_S_S_a-2.0E0*1*I_ERI_D2y_F2xy_S_S_b+1*I_ERI_D2y_Py_S_S;
  abcd[738] = 4.0E0*I_ERI_Gy3z_F2xy_S_S_ab-2.0E0*1*I_ERI_Gy3z_Py_S_S_a-2.0E0*2*I_ERI_Dyz_F2xy_S_S_b+2*1*I_ERI_Dyz_Py_S_S;
  abcd[739] = 4.0E0*I_ERI_G4z_F2xy_S_S_ab-2.0E0*1*I_ERI_G4z_Py_S_S_a-2.0E0*3*I_ERI_D2z_F2xy_S_S_b+3*1*I_ERI_D2z_Py_S_S;
  abcd[740] = 4.0E0*I_ERI_G3xz_F2xz_S_S_ab-2.0E0*1*I_ERI_G3xz_Pz_S_S_a;
  abcd[741] = 4.0E0*I_ERI_G2xyz_F2xz_S_S_ab-2.0E0*1*I_ERI_G2xyz_Pz_S_S_a;
  abcd[742] = 4.0E0*I_ERI_G2x2z_F2xz_S_S_ab-2.0E0*1*I_ERI_G2x2z_Pz_S_S_a-2.0E0*1*I_ERI_D2x_F2xz_S_S_b+1*I_ERI_D2x_Pz_S_S;
  abcd[743] = 4.0E0*I_ERI_Gx2yz_F2xz_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Pz_S_S_a;
  abcd[744] = 4.0E0*I_ERI_Gxy2z_F2xz_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Pz_S_S_a-2.0E0*1*I_ERI_Dxy_F2xz_S_S_b+1*I_ERI_Dxy_Pz_S_S;
  abcd[745] = 4.0E0*I_ERI_Gx3z_F2xz_S_S_ab-2.0E0*1*I_ERI_Gx3z_Pz_S_S_a-2.0E0*2*I_ERI_Dxz_F2xz_S_S_b+2*1*I_ERI_Dxz_Pz_S_S;
  abcd[746] = 4.0E0*I_ERI_G3yz_F2xz_S_S_ab-2.0E0*1*I_ERI_G3yz_Pz_S_S_a;
  abcd[747] = 4.0E0*I_ERI_G2y2z_F2xz_S_S_ab-2.0E0*1*I_ERI_G2y2z_Pz_S_S_a-2.0E0*1*I_ERI_D2y_F2xz_S_S_b+1*I_ERI_D2y_Pz_S_S;
  abcd[748] = 4.0E0*I_ERI_Gy3z_F2xz_S_S_ab-2.0E0*1*I_ERI_Gy3z_Pz_S_S_a-2.0E0*2*I_ERI_Dyz_F2xz_S_S_b+2*1*I_ERI_Dyz_Pz_S_S;
  abcd[749] = 4.0E0*I_ERI_G4z_F2xz_S_S_ab-2.0E0*1*I_ERI_G4z_Pz_S_S_a-2.0E0*3*I_ERI_D2z_F2xz_S_S_b+3*1*I_ERI_D2z_Pz_S_S;
  abcd[750] = 4.0E0*I_ERI_G3xz_Fx2y_S_S_ab;
  abcd[751] = 4.0E0*I_ERI_G2xyz_Fx2y_S_S_ab;
  abcd[752] = 4.0E0*I_ERI_G2x2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_D2x_Fx2y_S_S_b;
  abcd[753] = 4.0E0*I_ERI_Gx2yz_Fx2y_S_S_ab;
  abcd[754] = 4.0E0*I_ERI_Gxy2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_Dxy_Fx2y_S_S_b;
  abcd[755] = 4.0E0*I_ERI_Gx3z_Fx2y_S_S_ab-2.0E0*2*I_ERI_Dxz_Fx2y_S_S_b;
  abcd[756] = 4.0E0*I_ERI_G3yz_Fx2y_S_S_ab;
  abcd[757] = 4.0E0*I_ERI_G2y2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_D2y_Fx2y_S_S_b;
  abcd[758] = 4.0E0*I_ERI_Gy3z_Fx2y_S_S_ab-2.0E0*2*I_ERI_Dyz_Fx2y_S_S_b;
  abcd[759] = 4.0E0*I_ERI_G4z_Fx2y_S_S_ab-2.0E0*3*I_ERI_D2z_Fx2y_S_S_b;
  abcd[760] = 4.0E0*I_ERI_G3xz_Fxyz_S_S_ab;
  abcd[761] = 4.0E0*I_ERI_G2xyz_Fxyz_S_S_ab;
  abcd[762] = 4.0E0*I_ERI_G2x2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2x_Fxyz_S_S_b;
  abcd[763] = 4.0E0*I_ERI_Gx2yz_Fxyz_S_S_ab;
  abcd[764] = 4.0E0*I_ERI_Gxy2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_Dxy_Fxyz_S_S_b;
  abcd[765] = 4.0E0*I_ERI_Gx3z_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxz_Fxyz_S_S_b;
  abcd[766] = 4.0E0*I_ERI_G3yz_Fxyz_S_S_ab;
  abcd[767] = 4.0E0*I_ERI_G2y2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2y_Fxyz_S_S_b;
  abcd[768] = 4.0E0*I_ERI_Gy3z_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dyz_Fxyz_S_S_b;
  abcd[769] = 4.0E0*I_ERI_G4z_Fxyz_S_S_ab-2.0E0*3*I_ERI_D2z_Fxyz_S_S_b;
  abcd[770] = 4.0E0*I_ERI_G3xz_Fx2z_S_S_ab;
  abcd[771] = 4.0E0*I_ERI_G2xyz_Fx2z_S_S_ab;
  abcd[772] = 4.0E0*I_ERI_G2x2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_D2x_Fx2z_S_S_b;
  abcd[773] = 4.0E0*I_ERI_Gx2yz_Fx2z_S_S_ab;
  abcd[774] = 4.0E0*I_ERI_Gxy2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_Dxy_Fx2z_S_S_b;
  abcd[775] = 4.0E0*I_ERI_Gx3z_Fx2z_S_S_ab-2.0E0*2*I_ERI_Dxz_Fx2z_S_S_b;
  abcd[776] = 4.0E0*I_ERI_G3yz_Fx2z_S_S_ab;
  abcd[777] = 4.0E0*I_ERI_G2y2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_D2y_Fx2z_S_S_b;
  abcd[778] = 4.0E0*I_ERI_Gy3z_Fx2z_S_S_ab-2.0E0*2*I_ERI_Dyz_Fx2z_S_S_b;
  abcd[779] = 4.0E0*I_ERI_G4z_Fx2z_S_S_ab-2.0E0*3*I_ERI_D2z_Fx2z_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_D_F_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  abcd[780] = 4.0E0*I_ERI_G3xz_F2xy_S_S_ab;
  abcd[781] = 4.0E0*I_ERI_G2xyz_F2xy_S_S_ab;
  abcd[782] = 4.0E0*I_ERI_G2x2z_F2xy_S_S_ab-2.0E0*1*I_ERI_D2x_F2xy_S_S_b;
  abcd[783] = 4.0E0*I_ERI_Gx2yz_F2xy_S_S_ab;
  abcd[784] = 4.0E0*I_ERI_Gxy2z_F2xy_S_S_ab-2.0E0*1*I_ERI_Dxy_F2xy_S_S_b;
  abcd[785] = 4.0E0*I_ERI_Gx3z_F2xy_S_S_ab-2.0E0*2*I_ERI_Dxz_F2xy_S_S_b;
  abcd[786] = 4.0E0*I_ERI_G3yz_F2xy_S_S_ab;
  abcd[787] = 4.0E0*I_ERI_G2y2z_F2xy_S_S_ab-2.0E0*1*I_ERI_D2y_F2xy_S_S_b;
  abcd[788] = 4.0E0*I_ERI_Gy3z_F2xy_S_S_ab-2.0E0*2*I_ERI_Dyz_F2xy_S_S_b;
  abcd[789] = 4.0E0*I_ERI_G4z_F2xy_S_S_ab-2.0E0*3*I_ERI_D2z_F2xy_S_S_b;
  abcd[790] = 4.0E0*I_ERI_G3xz_Fx2y_S_S_ab-2.0E0*1*I_ERI_G3xz_Px_S_S_a;
  abcd[791] = 4.0E0*I_ERI_G2xyz_Fx2y_S_S_ab-2.0E0*1*I_ERI_G2xyz_Px_S_S_a;
  abcd[792] = 4.0E0*I_ERI_G2x2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_G2x2z_Px_S_S_a-2.0E0*1*I_ERI_D2x_Fx2y_S_S_b+1*I_ERI_D2x_Px_S_S;
  abcd[793] = 4.0E0*I_ERI_Gx2yz_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Px_S_S_a;
  abcd[794] = 4.0E0*I_ERI_Gxy2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Px_S_S_a-2.0E0*1*I_ERI_Dxy_Fx2y_S_S_b+1*I_ERI_Dxy_Px_S_S;
  abcd[795] = 4.0E0*I_ERI_Gx3z_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gx3z_Px_S_S_a-2.0E0*2*I_ERI_Dxz_Fx2y_S_S_b+2*1*I_ERI_Dxz_Px_S_S;
  abcd[796] = 4.0E0*I_ERI_G3yz_Fx2y_S_S_ab-2.0E0*1*I_ERI_G3yz_Px_S_S_a;
  abcd[797] = 4.0E0*I_ERI_G2y2z_Fx2y_S_S_ab-2.0E0*1*I_ERI_G2y2z_Px_S_S_a-2.0E0*1*I_ERI_D2y_Fx2y_S_S_b+1*I_ERI_D2y_Px_S_S;
  abcd[798] = 4.0E0*I_ERI_Gy3z_Fx2y_S_S_ab-2.0E0*1*I_ERI_Gy3z_Px_S_S_a-2.0E0*2*I_ERI_Dyz_Fx2y_S_S_b+2*1*I_ERI_Dyz_Px_S_S;
  abcd[799] = 4.0E0*I_ERI_G4z_Fx2y_S_S_ab-2.0E0*1*I_ERI_G4z_Px_S_S_a-2.0E0*3*I_ERI_D2z_Fx2y_S_S_b+3*1*I_ERI_D2z_Px_S_S;
  abcd[800] = 4.0E0*I_ERI_G3xz_Fxyz_S_S_ab;
  abcd[801] = 4.0E0*I_ERI_G2xyz_Fxyz_S_S_ab;
  abcd[802] = 4.0E0*I_ERI_G2x2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2x_Fxyz_S_S_b;
  abcd[803] = 4.0E0*I_ERI_Gx2yz_Fxyz_S_S_ab;
  abcd[804] = 4.0E0*I_ERI_Gxy2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_Dxy_Fxyz_S_S_b;
  abcd[805] = 4.0E0*I_ERI_Gx3z_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxz_Fxyz_S_S_b;
  abcd[806] = 4.0E0*I_ERI_G3yz_Fxyz_S_S_ab;
  abcd[807] = 4.0E0*I_ERI_G2y2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2y_Fxyz_S_S_b;
  abcd[808] = 4.0E0*I_ERI_Gy3z_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dyz_Fxyz_S_S_b;
  abcd[809] = 4.0E0*I_ERI_G4z_Fxyz_S_S_ab-2.0E0*3*I_ERI_D2z_Fxyz_S_S_b;
  abcd[810] = 4.0E0*I_ERI_G3xz_F3y_S_S_ab-2.0E0*2*I_ERI_G3xz_Py_S_S_a;
  abcd[811] = 4.0E0*I_ERI_G2xyz_F3y_S_S_ab-2.0E0*2*I_ERI_G2xyz_Py_S_S_a;
  abcd[812] = 4.0E0*I_ERI_G2x2z_F3y_S_S_ab-2.0E0*2*I_ERI_G2x2z_Py_S_S_a-2.0E0*1*I_ERI_D2x_F3y_S_S_b+2*I_ERI_D2x_Py_S_S;
  abcd[813] = 4.0E0*I_ERI_Gx2yz_F3y_S_S_ab-2.0E0*2*I_ERI_Gx2yz_Py_S_S_a;
  abcd[814] = 4.0E0*I_ERI_Gxy2z_F3y_S_S_ab-2.0E0*2*I_ERI_Gxy2z_Py_S_S_a-2.0E0*1*I_ERI_Dxy_F3y_S_S_b+2*I_ERI_Dxy_Py_S_S;
  abcd[815] = 4.0E0*I_ERI_Gx3z_F3y_S_S_ab-2.0E0*2*I_ERI_Gx3z_Py_S_S_a-2.0E0*2*I_ERI_Dxz_F3y_S_S_b+2*2*I_ERI_Dxz_Py_S_S;
  abcd[816] = 4.0E0*I_ERI_G3yz_F3y_S_S_ab-2.0E0*2*I_ERI_G3yz_Py_S_S_a;
  abcd[817] = 4.0E0*I_ERI_G2y2z_F3y_S_S_ab-2.0E0*2*I_ERI_G2y2z_Py_S_S_a-2.0E0*1*I_ERI_D2y_F3y_S_S_b+2*I_ERI_D2y_Py_S_S;
  abcd[818] = 4.0E0*I_ERI_Gy3z_F3y_S_S_ab-2.0E0*2*I_ERI_Gy3z_Py_S_S_a-2.0E0*2*I_ERI_Dyz_F3y_S_S_b+2*2*I_ERI_Dyz_Py_S_S;
  abcd[819] = 4.0E0*I_ERI_G4z_F3y_S_S_ab-2.0E0*2*I_ERI_G4z_Py_S_S_a-2.0E0*3*I_ERI_D2z_F3y_S_S_b+3*2*I_ERI_D2z_Py_S_S;
  abcd[820] = 4.0E0*I_ERI_G3xz_F2yz_S_S_ab-2.0E0*1*I_ERI_G3xz_Pz_S_S_a;
  abcd[821] = 4.0E0*I_ERI_G2xyz_F2yz_S_S_ab-2.0E0*1*I_ERI_G2xyz_Pz_S_S_a;
  abcd[822] = 4.0E0*I_ERI_G2x2z_F2yz_S_S_ab-2.0E0*1*I_ERI_G2x2z_Pz_S_S_a-2.0E0*1*I_ERI_D2x_F2yz_S_S_b+1*I_ERI_D2x_Pz_S_S;
  abcd[823] = 4.0E0*I_ERI_Gx2yz_F2yz_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Pz_S_S_a;
  abcd[824] = 4.0E0*I_ERI_Gxy2z_F2yz_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Pz_S_S_a-2.0E0*1*I_ERI_Dxy_F2yz_S_S_b+1*I_ERI_Dxy_Pz_S_S;
  abcd[825] = 4.0E0*I_ERI_Gx3z_F2yz_S_S_ab-2.0E0*1*I_ERI_Gx3z_Pz_S_S_a-2.0E0*2*I_ERI_Dxz_F2yz_S_S_b+2*1*I_ERI_Dxz_Pz_S_S;
  abcd[826] = 4.0E0*I_ERI_G3yz_F2yz_S_S_ab-2.0E0*1*I_ERI_G3yz_Pz_S_S_a;
  abcd[827] = 4.0E0*I_ERI_G2y2z_F2yz_S_S_ab-2.0E0*1*I_ERI_G2y2z_Pz_S_S_a-2.0E0*1*I_ERI_D2y_F2yz_S_S_b+1*I_ERI_D2y_Pz_S_S;
  abcd[828] = 4.0E0*I_ERI_Gy3z_F2yz_S_S_ab-2.0E0*1*I_ERI_Gy3z_Pz_S_S_a-2.0E0*2*I_ERI_Dyz_F2yz_S_S_b+2*1*I_ERI_Dyz_Pz_S_S;
  abcd[829] = 4.0E0*I_ERI_G4z_F2yz_S_S_ab-2.0E0*1*I_ERI_G4z_Pz_S_S_a-2.0E0*3*I_ERI_D2z_F2yz_S_S_b+3*1*I_ERI_D2z_Pz_S_S;
  abcd[830] = 4.0E0*I_ERI_G3xz_Fy2z_S_S_ab;
  abcd[831] = 4.0E0*I_ERI_G2xyz_Fy2z_S_S_ab;
  abcd[832] = 4.0E0*I_ERI_G2x2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_D2x_Fy2z_S_S_b;
  abcd[833] = 4.0E0*I_ERI_Gx2yz_Fy2z_S_S_ab;
  abcd[834] = 4.0E0*I_ERI_Gxy2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_Dxy_Fy2z_S_S_b;
  abcd[835] = 4.0E0*I_ERI_Gx3z_Fy2z_S_S_ab-2.0E0*2*I_ERI_Dxz_Fy2z_S_S_b;
  abcd[836] = 4.0E0*I_ERI_G3yz_Fy2z_S_S_ab;
  abcd[837] = 4.0E0*I_ERI_G2y2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_D2y_Fy2z_S_S_b;
  abcd[838] = 4.0E0*I_ERI_Gy3z_Fy2z_S_S_ab-2.0E0*2*I_ERI_Dyz_Fy2z_S_S_b;
  abcd[839] = 4.0E0*I_ERI_G4z_Fy2z_S_S_ab-2.0E0*3*I_ERI_D2z_Fy2z_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_F_S_S_ab
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   * RHS shell quartet name: SQ_ERI_D_F_S_S_b
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  abcd[840] = 4.0E0*I_ERI_G3xz_F2xz_S_S_ab;
  abcd[841] = 4.0E0*I_ERI_G2xyz_F2xz_S_S_ab;
  abcd[842] = 4.0E0*I_ERI_G2x2z_F2xz_S_S_ab-2.0E0*1*I_ERI_D2x_F2xz_S_S_b;
  abcd[843] = 4.0E0*I_ERI_Gx2yz_F2xz_S_S_ab;
  abcd[844] = 4.0E0*I_ERI_Gxy2z_F2xz_S_S_ab-2.0E0*1*I_ERI_Dxy_F2xz_S_S_b;
  abcd[845] = 4.0E0*I_ERI_Gx3z_F2xz_S_S_ab-2.0E0*2*I_ERI_Dxz_F2xz_S_S_b;
  abcd[846] = 4.0E0*I_ERI_G3yz_F2xz_S_S_ab;
  abcd[847] = 4.0E0*I_ERI_G2y2z_F2xz_S_S_ab-2.0E0*1*I_ERI_D2y_F2xz_S_S_b;
  abcd[848] = 4.0E0*I_ERI_Gy3z_F2xz_S_S_ab-2.0E0*2*I_ERI_Dyz_F2xz_S_S_b;
  abcd[849] = 4.0E0*I_ERI_G4z_F2xz_S_S_ab-2.0E0*3*I_ERI_D2z_F2xz_S_S_b;
  abcd[850] = 4.0E0*I_ERI_G3xz_Fxyz_S_S_ab;
  abcd[851] = 4.0E0*I_ERI_G2xyz_Fxyz_S_S_ab;
  abcd[852] = 4.0E0*I_ERI_G2x2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2x_Fxyz_S_S_b;
  abcd[853] = 4.0E0*I_ERI_Gx2yz_Fxyz_S_S_ab;
  abcd[854] = 4.0E0*I_ERI_Gxy2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_Dxy_Fxyz_S_S_b;
  abcd[855] = 4.0E0*I_ERI_Gx3z_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dxz_Fxyz_S_S_b;
  abcd[856] = 4.0E0*I_ERI_G3yz_Fxyz_S_S_ab;
  abcd[857] = 4.0E0*I_ERI_G2y2z_Fxyz_S_S_ab-2.0E0*1*I_ERI_D2y_Fxyz_S_S_b;
  abcd[858] = 4.0E0*I_ERI_Gy3z_Fxyz_S_S_ab-2.0E0*2*I_ERI_Dyz_Fxyz_S_S_b;
  abcd[859] = 4.0E0*I_ERI_G4z_Fxyz_S_S_ab-2.0E0*3*I_ERI_D2z_Fxyz_S_S_b;
  abcd[860] = 4.0E0*I_ERI_G3xz_Fx2z_S_S_ab-2.0E0*1*I_ERI_G3xz_Px_S_S_a;
  abcd[861] = 4.0E0*I_ERI_G2xyz_Fx2z_S_S_ab-2.0E0*1*I_ERI_G2xyz_Px_S_S_a;
  abcd[862] = 4.0E0*I_ERI_G2x2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_G2x2z_Px_S_S_a-2.0E0*1*I_ERI_D2x_Fx2z_S_S_b+1*I_ERI_D2x_Px_S_S;
  abcd[863] = 4.0E0*I_ERI_Gx2yz_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Px_S_S_a;
  abcd[864] = 4.0E0*I_ERI_Gxy2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Px_S_S_a-2.0E0*1*I_ERI_Dxy_Fx2z_S_S_b+1*I_ERI_Dxy_Px_S_S;
  abcd[865] = 4.0E0*I_ERI_Gx3z_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gx3z_Px_S_S_a-2.0E0*2*I_ERI_Dxz_Fx2z_S_S_b+2*1*I_ERI_Dxz_Px_S_S;
  abcd[866] = 4.0E0*I_ERI_G3yz_Fx2z_S_S_ab-2.0E0*1*I_ERI_G3yz_Px_S_S_a;
  abcd[867] = 4.0E0*I_ERI_G2y2z_Fx2z_S_S_ab-2.0E0*1*I_ERI_G2y2z_Px_S_S_a-2.0E0*1*I_ERI_D2y_Fx2z_S_S_b+1*I_ERI_D2y_Px_S_S;
  abcd[868] = 4.0E0*I_ERI_Gy3z_Fx2z_S_S_ab-2.0E0*1*I_ERI_Gy3z_Px_S_S_a-2.0E0*2*I_ERI_Dyz_Fx2z_S_S_b+2*1*I_ERI_Dyz_Px_S_S;
  abcd[869] = 4.0E0*I_ERI_G4z_Fx2z_S_S_ab-2.0E0*1*I_ERI_G4z_Px_S_S_a-2.0E0*3*I_ERI_D2z_Fx2z_S_S_b+3*1*I_ERI_D2z_Px_S_S;
  abcd[870] = 4.0E0*I_ERI_G3xz_F2yz_S_S_ab;
  abcd[871] = 4.0E0*I_ERI_G2xyz_F2yz_S_S_ab;
  abcd[872] = 4.0E0*I_ERI_G2x2z_F2yz_S_S_ab-2.0E0*1*I_ERI_D2x_F2yz_S_S_b;
  abcd[873] = 4.0E0*I_ERI_Gx2yz_F2yz_S_S_ab;
  abcd[874] = 4.0E0*I_ERI_Gxy2z_F2yz_S_S_ab-2.0E0*1*I_ERI_Dxy_F2yz_S_S_b;
  abcd[875] = 4.0E0*I_ERI_Gx3z_F2yz_S_S_ab-2.0E0*2*I_ERI_Dxz_F2yz_S_S_b;
  abcd[876] = 4.0E0*I_ERI_G3yz_F2yz_S_S_ab;
  abcd[877] = 4.0E0*I_ERI_G2y2z_F2yz_S_S_ab-2.0E0*1*I_ERI_D2y_F2yz_S_S_b;
  abcd[878] = 4.0E0*I_ERI_Gy3z_F2yz_S_S_ab-2.0E0*2*I_ERI_Dyz_F2yz_S_S_b;
  abcd[879] = 4.0E0*I_ERI_G4z_F2yz_S_S_ab-2.0E0*3*I_ERI_D2z_F2yz_S_S_b;
  abcd[880] = 4.0E0*I_ERI_G3xz_Fy2z_S_S_ab-2.0E0*1*I_ERI_G3xz_Py_S_S_a;
  abcd[881] = 4.0E0*I_ERI_G2xyz_Fy2z_S_S_ab-2.0E0*1*I_ERI_G2xyz_Py_S_S_a;
  abcd[882] = 4.0E0*I_ERI_G2x2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_G2x2z_Py_S_S_a-2.0E0*1*I_ERI_D2x_Fy2z_S_S_b+1*I_ERI_D2x_Py_S_S;
  abcd[883] = 4.0E0*I_ERI_Gx2yz_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gx2yz_Py_S_S_a;
  abcd[884] = 4.0E0*I_ERI_Gxy2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gxy2z_Py_S_S_a-2.0E0*1*I_ERI_Dxy_Fy2z_S_S_b+1*I_ERI_Dxy_Py_S_S;
  abcd[885] = 4.0E0*I_ERI_Gx3z_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gx3z_Py_S_S_a-2.0E0*2*I_ERI_Dxz_Fy2z_S_S_b+2*1*I_ERI_Dxz_Py_S_S;
  abcd[886] = 4.0E0*I_ERI_G3yz_Fy2z_S_S_ab-2.0E0*1*I_ERI_G3yz_Py_S_S_a;
  abcd[887] = 4.0E0*I_ERI_G2y2z_Fy2z_S_S_ab-2.0E0*1*I_ERI_G2y2z_Py_S_S_a-2.0E0*1*I_ERI_D2y_Fy2z_S_S_b+1*I_ERI_D2y_Py_S_S;
  abcd[888] = 4.0E0*I_ERI_Gy3z_Fy2z_S_S_ab-2.0E0*1*I_ERI_Gy3z_Py_S_S_a-2.0E0*2*I_ERI_Dyz_Fy2z_S_S_b+2*1*I_ERI_Dyz_Py_S_S;
  abcd[889] = 4.0E0*I_ERI_G4z_Fy2z_S_S_ab-2.0E0*1*I_ERI_G4z_Py_S_S_a-2.0E0*3*I_ERI_D2z_Fy2z_S_S_b+3*1*I_ERI_D2z_Py_S_S;
  abcd[890] = 4.0E0*I_ERI_G3xz_F3z_S_S_ab-2.0E0*2*I_ERI_G3xz_Pz_S_S_a;
  abcd[891] = 4.0E0*I_ERI_G2xyz_F3z_S_S_ab-2.0E0*2*I_ERI_G2xyz_Pz_S_S_a;
  abcd[892] = 4.0E0*I_ERI_G2x2z_F3z_S_S_ab-2.0E0*2*I_ERI_G2x2z_Pz_S_S_a-2.0E0*1*I_ERI_D2x_F3z_S_S_b+2*I_ERI_D2x_Pz_S_S;
  abcd[893] = 4.0E0*I_ERI_Gx2yz_F3z_S_S_ab-2.0E0*2*I_ERI_Gx2yz_Pz_S_S_a;
  abcd[894] = 4.0E0*I_ERI_Gxy2z_F3z_S_S_ab-2.0E0*2*I_ERI_Gxy2z_Pz_S_S_a-2.0E0*1*I_ERI_Dxy_F3z_S_S_b+2*I_ERI_Dxy_Pz_S_S;
  abcd[895] = 4.0E0*I_ERI_Gx3z_F3z_S_S_ab-2.0E0*2*I_ERI_Gx3z_Pz_S_S_a-2.0E0*2*I_ERI_Dxz_F3z_S_S_b+2*2*I_ERI_Dxz_Pz_S_S;
  abcd[896] = 4.0E0*I_ERI_G3yz_F3z_S_S_ab-2.0E0*2*I_ERI_G3yz_Pz_S_S_a;
  abcd[897] = 4.0E0*I_ERI_G2y2z_F3z_S_S_ab-2.0E0*2*I_ERI_G2y2z_Pz_S_S_a-2.0E0*1*I_ERI_D2y_F3z_S_S_b+2*I_ERI_D2y_Pz_S_S;
  abcd[898] = 4.0E0*I_ERI_Gy3z_F3z_S_S_ab-2.0E0*2*I_ERI_Gy3z_Pz_S_S_a-2.0E0*2*I_ERI_Dyz_F3z_S_S_b+2*2*I_ERI_Dyz_Pz_S_S;
  abcd[899] = 4.0E0*I_ERI_G4z_F3z_S_S_ab-2.0E0*2*I_ERI_G4z_Pz_S_S_a-2.0E0*3*I_ERI_D2z_F3z_S_S_b+3*2*I_ERI_D2z_Pz_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_ac
   * RHS shell quartet name: SQ_ERI_D_D_P_S_c
   ************************************************************/
  abcd[900] = 4.0E0*I_ERI_G4x_D2x_Px_S_ac-2.0E0*3*I_ERI_D2x_D2x_Px_S_c;
  abcd[901] = 4.0E0*I_ERI_G3xy_D2x_Px_S_ac-2.0E0*2*I_ERI_Dxy_D2x_Px_S_c;
  abcd[902] = 4.0E0*I_ERI_G3xz_D2x_Px_S_ac-2.0E0*2*I_ERI_Dxz_D2x_Px_S_c;
  abcd[903] = 4.0E0*I_ERI_G2x2y_D2x_Px_S_ac-2.0E0*1*I_ERI_D2y_D2x_Px_S_c;
  abcd[904] = 4.0E0*I_ERI_G2xyz_D2x_Px_S_ac-2.0E0*1*I_ERI_Dyz_D2x_Px_S_c;
  abcd[905] = 4.0E0*I_ERI_G2x2z_D2x_Px_S_ac-2.0E0*1*I_ERI_D2z_D2x_Px_S_c;
  abcd[906] = 4.0E0*I_ERI_Gx3y_D2x_Px_S_ac;
  abcd[907] = 4.0E0*I_ERI_Gx2yz_D2x_Px_S_ac;
  abcd[908] = 4.0E0*I_ERI_Gxy2z_D2x_Px_S_ac;
  abcd[909] = 4.0E0*I_ERI_Gx3z_D2x_Px_S_ac;
  abcd[910] = 4.0E0*I_ERI_G4x_Dxy_Px_S_ac-2.0E0*3*I_ERI_D2x_Dxy_Px_S_c;
  abcd[911] = 4.0E0*I_ERI_G3xy_Dxy_Px_S_ac-2.0E0*2*I_ERI_Dxy_Dxy_Px_S_c;
  abcd[912] = 4.0E0*I_ERI_G3xz_Dxy_Px_S_ac-2.0E0*2*I_ERI_Dxz_Dxy_Px_S_c;
  abcd[913] = 4.0E0*I_ERI_G2x2y_Dxy_Px_S_ac-2.0E0*1*I_ERI_D2y_Dxy_Px_S_c;
  abcd[914] = 4.0E0*I_ERI_G2xyz_Dxy_Px_S_ac-2.0E0*1*I_ERI_Dyz_Dxy_Px_S_c;
  abcd[915] = 4.0E0*I_ERI_G2x2z_Dxy_Px_S_ac-2.0E0*1*I_ERI_D2z_Dxy_Px_S_c;
  abcd[916] = 4.0E0*I_ERI_Gx3y_Dxy_Px_S_ac;
  abcd[917] = 4.0E0*I_ERI_Gx2yz_Dxy_Px_S_ac;
  abcd[918] = 4.0E0*I_ERI_Gxy2z_Dxy_Px_S_ac;
  abcd[919] = 4.0E0*I_ERI_Gx3z_Dxy_Px_S_ac;
  abcd[920] = 4.0E0*I_ERI_G4x_Dxz_Px_S_ac-2.0E0*3*I_ERI_D2x_Dxz_Px_S_c;
  abcd[921] = 4.0E0*I_ERI_G3xy_Dxz_Px_S_ac-2.0E0*2*I_ERI_Dxy_Dxz_Px_S_c;
  abcd[922] = 4.0E0*I_ERI_G3xz_Dxz_Px_S_ac-2.0E0*2*I_ERI_Dxz_Dxz_Px_S_c;
  abcd[923] = 4.0E0*I_ERI_G2x2y_Dxz_Px_S_ac-2.0E0*1*I_ERI_D2y_Dxz_Px_S_c;
  abcd[924] = 4.0E0*I_ERI_G2xyz_Dxz_Px_S_ac-2.0E0*1*I_ERI_Dyz_Dxz_Px_S_c;
  abcd[925] = 4.0E0*I_ERI_G2x2z_Dxz_Px_S_ac-2.0E0*1*I_ERI_D2z_Dxz_Px_S_c;
  abcd[926] = 4.0E0*I_ERI_Gx3y_Dxz_Px_S_ac;
  abcd[927] = 4.0E0*I_ERI_Gx2yz_Dxz_Px_S_ac;
  abcd[928] = 4.0E0*I_ERI_Gxy2z_Dxz_Px_S_ac;
  abcd[929] = 4.0E0*I_ERI_Gx3z_Dxz_Px_S_ac;
  abcd[930] = 4.0E0*I_ERI_G4x_D2y_Px_S_ac-2.0E0*3*I_ERI_D2x_D2y_Px_S_c;
  abcd[931] = 4.0E0*I_ERI_G3xy_D2y_Px_S_ac-2.0E0*2*I_ERI_Dxy_D2y_Px_S_c;
  abcd[932] = 4.0E0*I_ERI_G3xz_D2y_Px_S_ac-2.0E0*2*I_ERI_Dxz_D2y_Px_S_c;
  abcd[933] = 4.0E0*I_ERI_G2x2y_D2y_Px_S_ac-2.0E0*1*I_ERI_D2y_D2y_Px_S_c;
  abcd[934] = 4.0E0*I_ERI_G2xyz_D2y_Px_S_ac-2.0E0*1*I_ERI_Dyz_D2y_Px_S_c;
  abcd[935] = 4.0E0*I_ERI_G2x2z_D2y_Px_S_ac-2.0E0*1*I_ERI_D2z_D2y_Px_S_c;
  abcd[936] = 4.0E0*I_ERI_Gx3y_D2y_Px_S_ac;
  abcd[937] = 4.0E0*I_ERI_Gx2yz_D2y_Px_S_ac;
  abcd[938] = 4.0E0*I_ERI_Gxy2z_D2y_Px_S_ac;
  abcd[939] = 4.0E0*I_ERI_Gx3z_D2y_Px_S_ac;
  abcd[940] = 4.0E0*I_ERI_G4x_Dyz_Px_S_ac-2.0E0*3*I_ERI_D2x_Dyz_Px_S_c;
  abcd[941] = 4.0E0*I_ERI_G3xy_Dyz_Px_S_ac-2.0E0*2*I_ERI_Dxy_Dyz_Px_S_c;
  abcd[942] = 4.0E0*I_ERI_G3xz_Dyz_Px_S_ac-2.0E0*2*I_ERI_Dxz_Dyz_Px_S_c;
  abcd[943] = 4.0E0*I_ERI_G2x2y_Dyz_Px_S_ac-2.0E0*1*I_ERI_D2y_Dyz_Px_S_c;
  abcd[944] = 4.0E0*I_ERI_G2xyz_Dyz_Px_S_ac-2.0E0*1*I_ERI_Dyz_Dyz_Px_S_c;
  abcd[945] = 4.0E0*I_ERI_G2x2z_Dyz_Px_S_ac-2.0E0*1*I_ERI_D2z_Dyz_Px_S_c;
  abcd[946] = 4.0E0*I_ERI_Gx3y_Dyz_Px_S_ac;
  abcd[947] = 4.0E0*I_ERI_Gx2yz_Dyz_Px_S_ac;
  abcd[948] = 4.0E0*I_ERI_Gxy2z_Dyz_Px_S_ac;
  abcd[949] = 4.0E0*I_ERI_Gx3z_Dyz_Px_S_ac;
  abcd[950] = 4.0E0*I_ERI_G4x_D2z_Px_S_ac-2.0E0*3*I_ERI_D2x_D2z_Px_S_c;
  abcd[951] = 4.0E0*I_ERI_G3xy_D2z_Px_S_ac-2.0E0*2*I_ERI_Dxy_D2z_Px_S_c;
  abcd[952] = 4.0E0*I_ERI_G3xz_D2z_Px_S_ac-2.0E0*2*I_ERI_Dxz_D2z_Px_S_c;
  abcd[953] = 4.0E0*I_ERI_G2x2y_D2z_Px_S_ac-2.0E0*1*I_ERI_D2y_D2z_Px_S_c;
  abcd[954] = 4.0E0*I_ERI_G2xyz_D2z_Px_S_ac-2.0E0*1*I_ERI_Dyz_D2z_Px_S_c;
  abcd[955] = 4.0E0*I_ERI_G2x2z_D2z_Px_S_ac-2.0E0*1*I_ERI_D2z_D2z_Px_S_c;
  abcd[956] = 4.0E0*I_ERI_Gx3y_D2z_Px_S_ac;
  abcd[957] = 4.0E0*I_ERI_Gx2yz_D2z_Px_S_ac;
  abcd[958] = 4.0E0*I_ERI_Gxy2z_D2z_Px_S_ac;
  abcd[959] = 4.0E0*I_ERI_Gx3z_D2z_Px_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_ac
   * RHS shell quartet name: SQ_ERI_D_D_P_S_c
   ************************************************************/
  abcd[960] = 4.0E0*I_ERI_G4x_D2x_Py_S_ac-2.0E0*3*I_ERI_D2x_D2x_Py_S_c;
  abcd[961] = 4.0E0*I_ERI_G3xy_D2x_Py_S_ac-2.0E0*2*I_ERI_Dxy_D2x_Py_S_c;
  abcd[962] = 4.0E0*I_ERI_G3xz_D2x_Py_S_ac-2.0E0*2*I_ERI_Dxz_D2x_Py_S_c;
  abcd[963] = 4.0E0*I_ERI_G2x2y_D2x_Py_S_ac-2.0E0*1*I_ERI_D2y_D2x_Py_S_c;
  abcd[964] = 4.0E0*I_ERI_G2xyz_D2x_Py_S_ac-2.0E0*1*I_ERI_Dyz_D2x_Py_S_c;
  abcd[965] = 4.0E0*I_ERI_G2x2z_D2x_Py_S_ac-2.0E0*1*I_ERI_D2z_D2x_Py_S_c;
  abcd[966] = 4.0E0*I_ERI_Gx3y_D2x_Py_S_ac;
  abcd[967] = 4.0E0*I_ERI_Gx2yz_D2x_Py_S_ac;
  abcd[968] = 4.0E0*I_ERI_Gxy2z_D2x_Py_S_ac;
  abcd[969] = 4.0E0*I_ERI_Gx3z_D2x_Py_S_ac;
  abcd[970] = 4.0E0*I_ERI_G4x_Dxy_Py_S_ac-2.0E0*3*I_ERI_D2x_Dxy_Py_S_c;
  abcd[971] = 4.0E0*I_ERI_G3xy_Dxy_Py_S_ac-2.0E0*2*I_ERI_Dxy_Dxy_Py_S_c;
  abcd[972] = 4.0E0*I_ERI_G3xz_Dxy_Py_S_ac-2.0E0*2*I_ERI_Dxz_Dxy_Py_S_c;
  abcd[973] = 4.0E0*I_ERI_G2x2y_Dxy_Py_S_ac-2.0E0*1*I_ERI_D2y_Dxy_Py_S_c;
  abcd[974] = 4.0E0*I_ERI_G2xyz_Dxy_Py_S_ac-2.0E0*1*I_ERI_Dyz_Dxy_Py_S_c;
  abcd[975] = 4.0E0*I_ERI_G2x2z_Dxy_Py_S_ac-2.0E0*1*I_ERI_D2z_Dxy_Py_S_c;
  abcd[976] = 4.0E0*I_ERI_Gx3y_Dxy_Py_S_ac;
  abcd[977] = 4.0E0*I_ERI_Gx2yz_Dxy_Py_S_ac;
  abcd[978] = 4.0E0*I_ERI_Gxy2z_Dxy_Py_S_ac;
  abcd[979] = 4.0E0*I_ERI_Gx3z_Dxy_Py_S_ac;
  abcd[980] = 4.0E0*I_ERI_G4x_Dxz_Py_S_ac-2.0E0*3*I_ERI_D2x_Dxz_Py_S_c;
  abcd[981] = 4.0E0*I_ERI_G3xy_Dxz_Py_S_ac-2.0E0*2*I_ERI_Dxy_Dxz_Py_S_c;
  abcd[982] = 4.0E0*I_ERI_G3xz_Dxz_Py_S_ac-2.0E0*2*I_ERI_Dxz_Dxz_Py_S_c;
  abcd[983] = 4.0E0*I_ERI_G2x2y_Dxz_Py_S_ac-2.0E0*1*I_ERI_D2y_Dxz_Py_S_c;
  abcd[984] = 4.0E0*I_ERI_G2xyz_Dxz_Py_S_ac-2.0E0*1*I_ERI_Dyz_Dxz_Py_S_c;
  abcd[985] = 4.0E0*I_ERI_G2x2z_Dxz_Py_S_ac-2.0E0*1*I_ERI_D2z_Dxz_Py_S_c;
  abcd[986] = 4.0E0*I_ERI_Gx3y_Dxz_Py_S_ac;
  abcd[987] = 4.0E0*I_ERI_Gx2yz_Dxz_Py_S_ac;
  abcd[988] = 4.0E0*I_ERI_Gxy2z_Dxz_Py_S_ac;
  abcd[989] = 4.0E0*I_ERI_Gx3z_Dxz_Py_S_ac;
  abcd[990] = 4.0E0*I_ERI_G4x_D2y_Py_S_ac-2.0E0*3*I_ERI_D2x_D2y_Py_S_c;
  abcd[991] = 4.0E0*I_ERI_G3xy_D2y_Py_S_ac-2.0E0*2*I_ERI_Dxy_D2y_Py_S_c;
  abcd[992] = 4.0E0*I_ERI_G3xz_D2y_Py_S_ac-2.0E0*2*I_ERI_Dxz_D2y_Py_S_c;
  abcd[993] = 4.0E0*I_ERI_G2x2y_D2y_Py_S_ac-2.0E0*1*I_ERI_D2y_D2y_Py_S_c;
  abcd[994] = 4.0E0*I_ERI_G2xyz_D2y_Py_S_ac-2.0E0*1*I_ERI_Dyz_D2y_Py_S_c;
  abcd[995] = 4.0E0*I_ERI_G2x2z_D2y_Py_S_ac-2.0E0*1*I_ERI_D2z_D2y_Py_S_c;
  abcd[996] = 4.0E0*I_ERI_Gx3y_D2y_Py_S_ac;
  abcd[997] = 4.0E0*I_ERI_Gx2yz_D2y_Py_S_ac;
  abcd[998] = 4.0E0*I_ERI_Gxy2z_D2y_Py_S_ac;
  abcd[999] = 4.0E0*I_ERI_Gx3z_D2y_Py_S_ac;
  abcd[1000] = 4.0E0*I_ERI_G4x_Dyz_Py_S_ac-2.0E0*3*I_ERI_D2x_Dyz_Py_S_c;
  abcd[1001] = 4.0E0*I_ERI_G3xy_Dyz_Py_S_ac-2.0E0*2*I_ERI_Dxy_Dyz_Py_S_c;
  abcd[1002] = 4.0E0*I_ERI_G3xz_Dyz_Py_S_ac-2.0E0*2*I_ERI_Dxz_Dyz_Py_S_c;
  abcd[1003] = 4.0E0*I_ERI_G2x2y_Dyz_Py_S_ac-2.0E0*1*I_ERI_D2y_Dyz_Py_S_c;
  abcd[1004] = 4.0E0*I_ERI_G2xyz_Dyz_Py_S_ac-2.0E0*1*I_ERI_Dyz_Dyz_Py_S_c;
  abcd[1005] = 4.0E0*I_ERI_G2x2z_Dyz_Py_S_ac-2.0E0*1*I_ERI_D2z_Dyz_Py_S_c;
  abcd[1006] = 4.0E0*I_ERI_Gx3y_Dyz_Py_S_ac;
  abcd[1007] = 4.0E0*I_ERI_Gx2yz_Dyz_Py_S_ac;
  abcd[1008] = 4.0E0*I_ERI_Gxy2z_Dyz_Py_S_ac;
  abcd[1009] = 4.0E0*I_ERI_Gx3z_Dyz_Py_S_ac;
  abcd[1010] = 4.0E0*I_ERI_G4x_D2z_Py_S_ac-2.0E0*3*I_ERI_D2x_D2z_Py_S_c;
  abcd[1011] = 4.0E0*I_ERI_G3xy_D2z_Py_S_ac-2.0E0*2*I_ERI_Dxy_D2z_Py_S_c;
  abcd[1012] = 4.0E0*I_ERI_G3xz_D2z_Py_S_ac-2.0E0*2*I_ERI_Dxz_D2z_Py_S_c;
  abcd[1013] = 4.0E0*I_ERI_G2x2y_D2z_Py_S_ac-2.0E0*1*I_ERI_D2y_D2z_Py_S_c;
  abcd[1014] = 4.0E0*I_ERI_G2xyz_D2z_Py_S_ac-2.0E0*1*I_ERI_Dyz_D2z_Py_S_c;
  abcd[1015] = 4.0E0*I_ERI_G2x2z_D2z_Py_S_ac-2.0E0*1*I_ERI_D2z_D2z_Py_S_c;
  abcd[1016] = 4.0E0*I_ERI_Gx3y_D2z_Py_S_ac;
  abcd[1017] = 4.0E0*I_ERI_Gx2yz_D2z_Py_S_ac;
  abcd[1018] = 4.0E0*I_ERI_Gxy2z_D2z_Py_S_ac;
  abcd[1019] = 4.0E0*I_ERI_Gx3z_D2z_Py_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_ac
   * RHS shell quartet name: SQ_ERI_D_D_P_S_c
   ************************************************************/
  abcd[1020] = 4.0E0*I_ERI_G4x_D2x_Pz_S_ac-2.0E0*3*I_ERI_D2x_D2x_Pz_S_c;
  abcd[1021] = 4.0E0*I_ERI_G3xy_D2x_Pz_S_ac-2.0E0*2*I_ERI_Dxy_D2x_Pz_S_c;
  abcd[1022] = 4.0E0*I_ERI_G3xz_D2x_Pz_S_ac-2.0E0*2*I_ERI_Dxz_D2x_Pz_S_c;
  abcd[1023] = 4.0E0*I_ERI_G2x2y_D2x_Pz_S_ac-2.0E0*1*I_ERI_D2y_D2x_Pz_S_c;
  abcd[1024] = 4.0E0*I_ERI_G2xyz_D2x_Pz_S_ac-2.0E0*1*I_ERI_Dyz_D2x_Pz_S_c;
  abcd[1025] = 4.0E0*I_ERI_G2x2z_D2x_Pz_S_ac-2.0E0*1*I_ERI_D2z_D2x_Pz_S_c;
  abcd[1026] = 4.0E0*I_ERI_Gx3y_D2x_Pz_S_ac;
  abcd[1027] = 4.0E0*I_ERI_Gx2yz_D2x_Pz_S_ac;
  abcd[1028] = 4.0E0*I_ERI_Gxy2z_D2x_Pz_S_ac;
  abcd[1029] = 4.0E0*I_ERI_Gx3z_D2x_Pz_S_ac;
  abcd[1030] = 4.0E0*I_ERI_G4x_Dxy_Pz_S_ac-2.0E0*3*I_ERI_D2x_Dxy_Pz_S_c;
  abcd[1031] = 4.0E0*I_ERI_G3xy_Dxy_Pz_S_ac-2.0E0*2*I_ERI_Dxy_Dxy_Pz_S_c;
  abcd[1032] = 4.0E0*I_ERI_G3xz_Dxy_Pz_S_ac-2.0E0*2*I_ERI_Dxz_Dxy_Pz_S_c;
  abcd[1033] = 4.0E0*I_ERI_G2x2y_Dxy_Pz_S_ac-2.0E0*1*I_ERI_D2y_Dxy_Pz_S_c;
  abcd[1034] = 4.0E0*I_ERI_G2xyz_Dxy_Pz_S_ac-2.0E0*1*I_ERI_Dyz_Dxy_Pz_S_c;
  abcd[1035] = 4.0E0*I_ERI_G2x2z_Dxy_Pz_S_ac-2.0E0*1*I_ERI_D2z_Dxy_Pz_S_c;
  abcd[1036] = 4.0E0*I_ERI_Gx3y_Dxy_Pz_S_ac;
  abcd[1037] = 4.0E0*I_ERI_Gx2yz_Dxy_Pz_S_ac;
  abcd[1038] = 4.0E0*I_ERI_Gxy2z_Dxy_Pz_S_ac;
  abcd[1039] = 4.0E0*I_ERI_Gx3z_Dxy_Pz_S_ac;
  abcd[1040] = 4.0E0*I_ERI_G4x_Dxz_Pz_S_ac-2.0E0*3*I_ERI_D2x_Dxz_Pz_S_c;
  abcd[1041] = 4.0E0*I_ERI_G3xy_Dxz_Pz_S_ac-2.0E0*2*I_ERI_Dxy_Dxz_Pz_S_c;
  abcd[1042] = 4.0E0*I_ERI_G3xz_Dxz_Pz_S_ac-2.0E0*2*I_ERI_Dxz_Dxz_Pz_S_c;
  abcd[1043] = 4.0E0*I_ERI_G2x2y_Dxz_Pz_S_ac-2.0E0*1*I_ERI_D2y_Dxz_Pz_S_c;
  abcd[1044] = 4.0E0*I_ERI_G2xyz_Dxz_Pz_S_ac-2.0E0*1*I_ERI_Dyz_Dxz_Pz_S_c;
  abcd[1045] = 4.0E0*I_ERI_G2x2z_Dxz_Pz_S_ac-2.0E0*1*I_ERI_D2z_Dxz_Pz_S_c;
  abcd[1046] = 4.0E0*I_ERI_Gx3y_Dxz_Pz_S_ac;
  abcd[1047] = 4.0E0*I_ERI_Gx2yz_Dxz_Pz_S_ac;
  abcd[1048] = 4.0E0*I_ERI_Gxy2z_Dxz_Pz_S_ac;
  abcd[1049] = 4.0E0*I_ERI_Gx3z_Dxz_Pz_S_ac;
  abcd[1050] = 4.0E0*I_ERI_G4x_D2y_Pz_S_ac-2.0E0*3*I_ERI_D2x_D2y_Pz_S_c;
  abcd[1051] = 4.0E0*I_ERI_G3xy_D2y_Pz_S_ac-2.0E0*2*I_ERI_Dxy_D2y_Pz_S_c;
  abcd[1052] = 4.0E0*I_ERI_G3xz_D2y_Pz_S_ac-2.0E0*2*I_ERI_Dxz_D2y_Pz_S_c;
  abcd[1053] = 4.0E0*I_ERI_G2x2y_D2y_Pz_S_ac-2.0E0*1*I_ERI_D2y_D2y_Pz_S_c;
  abcd[1054] = 4.0E0*I_ERI_G2xyz_D2y_Pz_S_ac-2.0E0*1*I_ERI_Dyz_D2y_Pz_S_c;
  abcd[1055] = 4.0E0*I_ERI_G2x2z_D2y_Pz_S_ac-2.0E0*1*I_ERI_D2z_D2y_Pz_S_c;
  abcd[1056] = 4.0E0*I_ERI_Gx3y_D2y_Pz_S_ac;
  abcd[1057] = 4.0E0*I_ERI_Gx2yz_D2y_Pz_S_ac;
  abcd[1058] = 4.0E0*I_ERI_Gxy2z_D2y_Pz_S_ac;
  abcd[1059] = 4.0E0*I_ERI_Gx3z_D2y_Pz_S_ac;
  abcd[1060] = 4.0E0*I_ERI_G4x_Dyz_Pz_S_ac-2.0E0*3*I_ERI_D2x_Dyz_Pz_S_c;
  abcd[1061] = 4.0E0*I_ERI_G3xy_Dyz_Pz_S_ac-2.0E0*2*I_ERI_Dxy_Dyz_Pz_S_c;
  abcd[1062] = 4.0E0*I_ERI_G3xz_Dyz_Pz_S_ac-2.0E0*2*I_ERI_Dxz_Dyz_Pz_S_c;
  abcd[1063] = 4.0E0*I_ERI_G2x2y_Dyz_Pz_S_ac-2.0E0*1*I_ERI_D2y_Dyz_Pz_S_c;
  abcd[1064] = 4.0E0*I_ERI_G2xyz_Dyz_Pz_S_ac-2.0E0*1*I_ERI_Dyz_Dyz_Pz_S_c;
  abcd[1065] = 4.0E0*I_ERI_G2x2z_Dyz_Pz_S_ac-2.0E0*1*I_ERI_D2z_Dyz_Pz_S_c;
  abcd[1066] = 4.0E0*I_ERI_Gx3y_Dyz_Pz_S_ac;
  abcd[1067] = 4.0E0*I_ERI_Gx2yz_Dyz_Pz_S_ac;
  abcd[1068] = 4.0E0*I_ERI_Gxy2z_Dyz_Pz_S_ac;
  abcd[1069] = 4.0E0*I_ERI_Gx3z_Dyz_Pz_S_ac;
  abcd[1070] = 4.0E0*I_ERI_G4x_D2z_Pz_S_ac-2.0E0*3*I_ERI_D2x_D2z_Pz_S_c;
  abcd[1071] = 4.0E0*I_ERI_G3xy_D2z_Pz_S_ac-2.0E0*2*I_ERI_Dxy_D2z_Pz_S_c;
  abcd[1072] = 4.0E0*I_ERI_G3xz_D2z_Pz_S_ac-2.0E0*2*I_ERI_Dxz_D2z_Pz_S_c;
  abcd[1073] = 4.0E0*I_ERI_G2x2y_D2z_Pz_S_ac-2.0E0*1*I_ERI_D2y_D2z_Pz_S_c;
  abcd[1074] = 4.0E0*I_ERI_G2xyz_D2z_Pz_S_ac-2.0E0*1*I_ERI_Dyz_D2z_Pz_S_c;
  abcd[1075] = 4.0E0*I_ERI_G2x2z_D2z_Pz_S_ac-2.0E0*1*I_ERI_D2z_D2z_Pz_S_c;
  abcd[1076] = 4.0E0*I_ERI_Gx3y_D2z_Pz_S_ac;
  abcd[1077] = 4.0E0*I_ERI_Gx2yz_D2z_Pz_S_ac;
  abcd[1078] = 4.0E0*I_ERI_Gxy2z_D2z_Pz_S_ac;
  abcd[1079] = 4.0E0*I_ERI_Gx3z_D2z_Pz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_ac
   * RHS shell quartet name: SQ_ERI_D_D_P_S_c
   ************************************************************/
  abcd[1080] = 4.0E0*I_ERI_G3xy_D2x_Px_S_ac;
  abcd[1081] = 4.0E0*I_ERI_G2x2y_D2x_Px_S_ac-2.0E0*1*I_ERI_D2x_D2x_Px_S_c;
  abcd[1082] = 4.0E0*I_ERI_G2xyz_D2x_Px_S_ac;
  abcd[1083] = 4.0E0*I_ERI_Gx3y_D2x_Px_S_ac-2.0E0*2*I_ERI_Dxy_D2x_Px_S_c;
  abcd[1084] = 4.0E0*I_ERI_Gx2yz_D2x_Px_S_ac-2.0E0*1*I_ERI_Dxz_D2x_Px_S_c;
  abcd[1085] = 4.0E0*I_ERI_Gxy2z_D2x_Px_S_ac;
  abcd[1086] = 4.0E0*I_ERI_G4y_D2x_Px_S_ac-2.0E0*3*I_ERI_D2y_D2x_Px_S_c;
  abcd[1087] = 4.0E0*I_ERI_G3yz_D2x_Px_S_ac-2.0E0*2*I_ERI_Dyz_D2x_Px_S_c;
  abcd[1088] = 4.0E0*I_ERI_G2y2z_D2x_Px_S_ac-2.0E0*1*I_ERI_D2z_D2x_Px_S_c;
  abcd[1089] = 4.0E0*I_ERI_Gy3z_D2x_Px_S_ac;
  abcd[1090] = 4.0E0*I_ERI_G3xy_Dxy_Px_S_ac;
  abcd[1091] = 4.0E0*I_ERI_G2x2y_Dxy_Px_S_ac-2.0E0*1*I_ERI_D2x_Dxy_Px_S_c;
  abcd[1092] = 4.0E0*I_ERI_G2xyz_Dxy_Px_S_ac;
  abcd[1093] = 4.0E0*I_ERI_Gx3y_Dxy_Px_S_ac-2.0E0*2*I_ERI_Dxy_Dxy_Px_S_c;
  abcd[1094] = 4.0E0*I_ERI_Gx2yz_Dxy_Px_S_ac-2.0E0*1*I_ERI_Dxz_Dxy_Px_S_c;
  abcd[1095] = 4.0E0*I_ERI_Gxy2z_Dxy_Px_S_ac;
  abcd[1096] = 4.0E0*I_ERI_G4y_Dxy_Px_S_ac-2.0E0*3*I_ERI_D2y_Dxy_Px_S_c;
  abcd[1097] = 4.0E0*I_ERI_G3yz_Dxy_Px_S_ac-2.0E0*2*I_ERI_Dyz_Dxy_Px_S_c;
  abcd[1098] = 4.0E0*I_ERI_G2y2z_Dxy_Px_S_ac-2.0E0*1*I_ERI_D2z_Dxy_Px_S_c;
  abcd[1099] = 4.0E0*I_ERI_Gy3z_Dxy_Px_S_ac;
  abcd[1100] = 4.0E0*I_ERI_G3xy_Dxz_Px_S_ac;
  abcd[1101] = 4.0E0*I_ERI_G2x2y_Dxz_Px_S_ac-2.0E0*1*I_ERI_D2x_Dxz_Px_S_c;
  abcd[1102] = 4.0E0*I_ERI_G2xyz_Dxz_Px_S_ac;
  abcd[1103] = 4.0E0*I_ERI_Gx3y_Dxz_Px_S_ac-2.0E0*2*I_ERI_Dxy_Dxz_Px_S_c;
  abcd[1104] = 4.0E0*I_ERI_Gx2yz_Dxz_Px_S_ac-2.0E0*1*I_ERI_Dxz_Dxz_Px_S_c;
  abcd[1105] = 4.0E0*I_ERI_Gxy2z_Dxz_Px_S_ac;
  abcd[1106] = 4.0E0*I_ERI_G4y_Dxz_Px_S_ac-2.0E0*3*I_ERI_D2y_Dxz_Px_S_c;
  abcd[1107] = 4.0E0*I_ERI_G3yz_Dxz_Px_S_ac-2.0E0*2*I_ERI_Dyz_Dxz_Px_S_c;
  abcd[1108] = 4.0E0*I_ERI_G2y2z_Dxz_Px_S_ac-2.0E0*1*I_ERI_D2z_Dxz_Px_S_c;
  abcd[1109] = 4.0E0*I_ERI_Gy3z_Dxz_Px_S_ac;
  abcd[1110] = 4.0E0*I_ERI_G3xy_D2y_Px_S_ac;
  abcd[1111] = 4.0E0*I_ERI_G2x2y_D2y_Px_S_ac-2.0E0*1*I_ERI_D2x_D2y_Px_S_c;
  abcd[1112] = 4.0E0*I_ERI_G2xyz_D2y_Px_S_ac;
  abcd[1113] = 4.0E0*I_ERI_Gx3y_D2y_Px_S_ac-2.0E0*2*I_ERI_Dxy_D2y_Px_S_c;
  abcd[1114] = 4.0E0*I_ERI_Gx2yz_D2y_Px_S_ac-2.0E0*1*I_ERI_Dxz_D2y_Px_S_c;
  abcd[1115] = 4.0E0*I_ERI_Gxy2z_D2y_Px_S_ac;
  abcd[1116] = 4.0E0*I_ERI_G4y_D2y_Px_S_ac-2.0E0*3*I_ERI_D2y_D2y_Px_S_c;
  abcd[1117] = 4.0E0*I_ERI_G3yz_D2y_Px_S_ac-2.0E0*2*I_ERI_Dyz_D2y_Px_S_c;
  abcd[1118] = 4.0E0*I_ERI_G2y2z_D2y_Px_S_ac-2.0E0*1*I_ERI_D2z_D2y_Px_S_c;
  abcd[1119] = 4.0E0*I_ERI_Gy3z_D2y_Px_S_ac;
  abcd[1120] = 4.0E0*I_ERI_G3xy_Dyz_Px_S_ac;
  abcd[1121] = 4.0E0*I_ERI_G2x2y_Dyz_Px_S_ac-2.0E0*1*I_ERI_D2x_Dyz_Px_S_c;
  abcd[1122] = 4.0E0*I_ERI_G2xyz_Dyz_Px_S_ac;
  abcd[1123] = 4.0E0*I_ERI_Gx3y_Dyz_Px_S_ac-2.0E0*2*I_ERI_Dxy_Dyz_Px_S_c;
  abcd[1124] = 4.0E0*I_ERI_Gx2yz_Dyz_Px_S_ac-2.0E0*1*I_ERI_Dxz_Dyz_Px_S_c;
  abcd[1125] = 4.0E0*I_ERI_Gxy2z_Dyz_Px_S_ac;
  abcd[1126] = 4.0E0*I_ERI_G4y_Dyz_Px_S_ac-2.0E0*3*I_ERI_D2y_Dyz_Px_S_c;
  abcd[1127] = 4.0E0*I_ERI_G3yz_Dyz_Px_S_ac-2.0E0*2*I_ERI_Dyz_Dyz_Px_S_c;
  abcd[1128] = 4.0E0*I_ERI_G2y2z_Dyz_Px_S_ac-2.0E0*1*I_ERI_D2z_Dyz_Px_S_c;
  abcd[1129] = 4.0E0*I_ERI_Gy3z_Dyz_Px_S_ac;
  abcd[1130] = 4.0E0*I_ERI_G3xy_D2z_Px_S_ac;
  abcd[1131] = 4.0E0*I_ERI_G2x2y_D2z_Px_S_ac-2.0E0*1*I_ERI_D2x_D2z_Px_S_c;
  abcd[1132] = 4.0E0*I_ERI_G2xyz_D2z_Px_S_ac;
  abcd[1133] = 4.0E0*I_ERI_Gx3y_D2z_Px_S_ac-2.0E0*2*I_ERI_Dxy_D2z_Px_S_c;
  abcd[1134] = 4.0E0*I_ERI_Gx2yz_D2z_Px_S_ac-2.0E0*1*I_ERI_Dxz_D2z_Px_S_c;
  abcd[1135] = 4.0E0*I_ERI_Gxy2z_D2z_Px_S_ac;
  abcd[1136] = 4.0E0*I_ERI_G4y_D2z_Px_S_ac-2.0E0*3*I_ERI_D2y_D2z_Px_S_c;
  abcd[1137] = 4.0E0*I_ERI_G3yz_D2z_Px_S_ac-2.0E0*2*I_ERI_Dyz_D2z_Px_S_c;
  abcd[1138] = 4.0E0*I_ERI_G2y2z_D2z_Px_S_ac-2.0E0*1*I_ERI_D2z_D2z_Px_S_c;
  abcd[1139] = 4.0E0*I_ERI_Gy3z_D2z_Px_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_ac
   * RHS shell quartet name: SQ_ERI_D_D_P_S_c
   ************************************************************/
  abcd[1140] = 4.0E0*I_ERI_G3xy_D2x_Py_S_ac;
  abcd[1141] = 4.0E0*I_ERI_G2x2y_D2x_Py_S_ac-2.0E0*1*I_ERI_D2x_D2x_Py_S_c;
  abcd[1142] = 4.0E0*I_ERI_G2xyz_D2x_Py_S_ac;
  abcd[1143] = 4.0E0*I_ERI_Gx3y_D2x_Py_S_ac-2.0E0*2*I_ERI_Dxy_D2x_Py_S_c;
  abcd[1144] = 4.0E0*I_ERI_Gx2yz_D2x_Py_S_ac-2.0E0*1*I_ERI_Dxz_D2x_Py_S_c;
  abcd[1145] = 4.0E0*I_ERI_Gxy2z_D2x_Py_S_ac;
  abcd[1146] = 4.0E0*I_ERI_G4y_D2x_Py_S_ac-2.0E0*3*I_ERI_D2y_D2x_Py_S_c;
  abcd[1147] = 4.0E0*I_ERI_G3yz_D2x_Py_S_ac-2.0E0*2*I_ERI_Dyz_D2x_Py_S_c;
  abcd[1148] = 4.0E0*I_ERI_G2y2z_D2x_Py_S_ac-2.0E0*1*I_ERI_D2z_D2x_Py_S_c;
  abcd[1149] = 4.0E0*I_ERI_Gy3z_D2x_Py_S_ac;
  abcd[1150] = 4.0E0*I_ERI_G3xy_Dxy_Py_S_ac;
  abcd[1151] = 4.0E0*I_ERI_G2x2y_Dxy_Py_S_ac-2.0E0*1*I_ERI_D2x_Dxy_Py_S_c;
  abcd[1152] = 4.0E0*I_ERI_G2xyz_Dxy_Py_S_ac;
  abcd[1153] = 4.0E0*I_ERI_Gx3y_Dxy_Py_S_ac-2.0E0*2*I_ERI_Dxy_Dxy_Py_S_c;
  abcd[1154] = 4.0E0*I_ERI_Gx2yz_Dxy_Py_S_ac-2.0E0*1*I_ERI_Dxz_Dxy_Py_S_c;
  abcd[1155] = 4.0E0*I_ERI_Gxy2z_Dxy_Py_S_ac;
  abcd[1156] = 4.0E0*I_ERI_G4y_Dxy_Py_S_ac-2.0E0*3*I_ERI_D2y_Dxy_Py_S_c;
  abcd[1157] = 4.0E0*I_ERI_G3yz_Dxy_Py_S_ac-2.0E0*2*I_ERI_Dyz_Dxy_Py_S_c;
  abcd[1158] = 4.0E0*I_ERI_G2y2z_Dxy_Py_S_ac-2.0E0*1*I_ERI_D2z_Dxy_Py_S_c;
  abcd[1159] = 4.0E0*I_ERI_Gy3z_Dxy_Py_S_ac;
  abcd[1160] = 4.0E0*I_ERI_G3xy_Dxz_Py_S_ac;
  abcd[1161] = 4.0E0*I_ERI_G2x2y_Dxz_Py_S_ac-2.0E0*1*I_ERI_D2x_Dxz_Py_S_c;
  abcd[1162] = 4.0E0*I_ERI_G2xyz_Dxz_Py_S_ac;
  abcd[1163] = 4.0E0*I_ERI_Gx3y_Dxz_Py_S_ac-2.0E0*2*I_ERI_Dxy_Dxz_Py_S_c;
  abcd[1164] = 4.0E0*I_ERI_Gx2yz_Dxz_Py_S_ac-2.0E0*1*I_ERI_Dxz_Dxz_Py_S_c;
  abcd[1165] = 4.0E0*I_ERI_Gxy2z_Dxz_Py_S_ac;
  abcd[1166] = 4.0E0*I_ERI_G4y_Dxz_Py_S_ac-2.0E0*3*I_ERI_D2y_Dxz_Py_S_c;
  abcd[1167] = 4.0E0*I_ERI_G3yz_Dxz_Py_S_ac-2.0E0*2*I_ERI_Dyz_Dxz_Py_S_c;
  abcd[1168] = 4.0E0*I_ERI_G2y2z_Dxz_Py_S_ac-2.0E0*1*I_ERI_D2z_Dxz_Py_S_c;
  abcd[1169] = 4.0E0*I_ERI_Gy3z_Dxz_Py_S_ac;
  abcd[1170] = 4.0E0*I_ERI_G3xy_D2y_Py_S_ac;
  abcd[1171] = 4.0E0*I_ERI_G2x2y_D2y_Py_S_ac-2.0E0*1*I_ERI_D2x_D2y_Py_S_c;
  abcd[1172] = 4.0E0*I_ERI_G2xyz_D2y_Py_S_ac;
  abcd[1173] = 4.0E0*I_ERI_Gx3y_D2y_Py_S_ac-2.0E0*2*I_ERI_Dxy_D2y_Py_S_c;
  abcd[1174] = 4.0E0*I_ERI_Gx2yz_D2y_Py_S_ac-2.0E0*1*I_ERI_Dxz_D2y_Py_S_c;
  abcd[1175] = 4.0E0*I_ERI_Gxy2z_D2y_Py_S_ac;
  abcd[1176] = 4.0E0*I_ERI_G4y_D2y_Py_S_ac-2.0E0*3*I_ERI_D2y_D2y_Py_S_c;
  abcd[1177] = 4.0E0*I_ERI_G3yz_D2y_Py_S_ac-2.0E0*2*I_ERI_Dyz_D2y_Py_S_c;
  abcd[1178] = 4.0E0*I_ERI_G2y2z_D2y_Py_S_ac-2.0E0*1*I_ERI_D2z_D2y_Py_S_c;
  abcd[1179] = 4.0E0*I_ERI_Gy3z_D2y_Py_S_ac;
  abcd[1180] = 4.0E0*I_ERI_G3xy_Dyz_Py_S_ac;
  abcd[1181] = 4.0E0*I_ERI_G2x2y_Dyz_Py_S_ac-2.0E0*1*I_ERI_D2x_Dyz_Py_S_c;
  abcd[1182] = 4.0E0*I_ERI_G2xyz_Dyz_Py_S_ac;
  abcd[1183] = 4.0E0*I_ERI_Gx3y_Dyz_Py_S_ac-2.0E0*2*I_ERI_Dxy_Dyz_Py_S_c;
  abcd[1184] = 4.0E0*I_ERI_Gx2yz_Dyz_Py_S_ac-2.0E0*1*I_ERI_Dxz_Dyz_Py_S_c;
  abcd[1185] = 4.0E0*I_ERI_Gxy2z_Dyz_Py_S_ac;
  abcd[1186] = 4.0E0*I_ERI_G4y_Dyz_Py_S_ac-2.0E0*3*I_ERI_D2y_Dyz_Py_S_c;
  abcd[1187] = 4.0E0*I_ERI_G3yz_Dyz_Py_S_ac-2.0E0*2*I_ERI_Dyz_Dyz_Py_S_c;
  abcd[1188] = 4.0E0*I_ERI_G2y2z_Dyz_Py_S_ac-2.0E0*1*I_ERI_D2z_Dyz_Py_S_c;
  abcd[1189] = 4.0E0*I_ERI_Gy3z_Dyz_Py_S_ac;
  abcd[1190] = 4.0E0*I_ERI_G3xy_D2z_Py_S_ac;
  abcd[1191] = 4.0E0*I_ERI_G2x2y_D2z_Py_S_ac-2.0E0*1*I_ERI_D2x_D2z_Py_S_c;
  abcd[1192] = 4.0E0*I_ERI_G2xyz_D2z_Py_S_ac;
  abcd[1193] = 4.0E0*I_ERI_Gx3y_D2z_Py_S_ac-2.0E0*2*I_ERI_Dxy_D2z_Py_S_c;
  abcd[1194] = 4.0E0*I_ERI_Gx2yz_D2z_Py_S_ac-2.0E0*1*I_ERI_Dxz_D2z_Py_S_c;
  abcd[1195] = 4.0E0*I_ERI_Gxy2z_D2z_Py_S_ac;
  abcd[1196] = 4.0E0*I_ERI_G4y_D2z_Py_S_ac-2.0E0*3*I_ERI_D2y_D2z_Py_S_c;
  abcd[1197] = 4.0E0*I_ERI_G3yz_D2z_Py_S_ac-2.0E0*2*I_ERI_Dyz_D2z_Py_S_c;
  abcd[1198] = 4.0E0*I_ERI_G2y2z_D2z_Py_S_ac-2.0E0*1*I_ERI_D2z_D2z_Py_S_c;
  abcd[1199] = 4.0E0*I_ERI_Gy3z_D2z_Py_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_ac
   * RHS shell quartet name: SQ_ERI_D_D_P_S_c
   ************************************************************/
  abcd[1200] = 4.0E0*I_ERI_G3xy_D2x_Pz_S_ac;
  abcd[1201] = 4.0E0*I_ERI_G2x2y_D2x_Pz_S_ac-2.0E0*1*I_ERI_D2x_D2x_Pz_S_c;
  abcd[1202] = 4.0E0*I_ERI_G2xyz_D2x_Pz_S_ac;
  abcd[1203] = 4.0E0*I_ERI_Gx3y_D2x_Pz_S_ac-2.0E0*2*I_ERI_Dxy_D2x_Pz_S_c;
  abcd[1204] = 4.0E0*I_ERI_Gx2yz_D2x_Pz_S_ac-2.0E0*1*I_ERI_Dxz_D2x_Pz_S_c;
  abcd[1205] = 4.0E0*I_ERI_Gxy2z_D2x_Pz_S_ac;
  abcd[1206] = 4.0E0*I_ERI_G4y_D2x_Pz_S_ac-2.0E0*3*I_ERI_D2y_D2x_Pz_S_c;
  abcd[1207] = 4.0E0*I_ERI_G3yz_D2x_Pz_S_ac-2.0E0*2*I_ERI_Dyz_D2x_Pz_S_c;
  abcd[1208] = 4.0E0*I_ERI_G2y2z_D2x_Pz_S_ac-2.0E0*1*I_ERI_D2z_D2x_Pz_S_c;
  abcd[1209] = 4.0E0*I_ERI_Gy3z_D2x_Pz_S_ac;
  abcd[1210] = 4.0E0*I_ERI_G3xy_Dxy_Pz_S_ac;
  abcd[1211] = 4.0E0*I_ERI_G2x2y_Dxy_Pz_S_ac-2.0E0*1*I_ERI_D2x_Dxy_Pz_S_c;
  abcd[1212] = 4.0E0*I_ERI_G2xyz_Dxy_Pz_S_ac;
  abcd[1213] = 4.0E0*I_ERI_Gx3y_Dxy_Pz_S_ac-2.0E0*2*I_ERI_Dxy_Dxy_Pz_S_c;
  abcd[1214] = 4.0E0*I_ERI_Gx2yz_Dxy_Pz_S_ac-2.0E0*1*I_ERI_Dxz_Dxy_Pz_S_c;
  abcd[1215] = 4.0E0*I_ERI_Gxy2z_Dxy_Pz_S_ac;
  abcd[1216] = 4.0E0*I_ERI_G4y_Dxy_Pz_S_ac-2.0E0*3*I_ERI_D2y_Dxy_Pz_S_c;
  abcd[1217] = 4.0E0*I_ERI_G3yz_Dxy_Pz_S_ac-2.0E0*2*I_ERI_Dyz_Dxy_Pz_S_c;
  abcd[1218] = 4.0E0*I_ERI_G2y2z_Dxy_Pz_S_ac-2.0E0*1*I_ERI_D2z_Dxy_Pz_S_c;
  abcd[1219] = 4.0E0*I_ERI_Gy3z_Dxy_Pz_S_ac;
  abcd[1220] = 4.0E0*I_ERI_G3xy_Dxz_Pz_S_ac;
  abcd[1221] = 4.0E0*I_ERI_G2x2y_Dxz_Pz_S_ac-2.0E0*1*I_ERI_D2x_Dxz_Pz_S_c;
  abcd[1222] = 4.0E0*I_ERI_G2xyz_Dxz_Pz_S_ac;
  abcd[1223] = 4.0E0*I_ERI_Gx3y_Dxz_Pz_S_ac-2.0E0*2*I_ERI_Dxy_Dxz_Pz_S_c;
  abcd[1224] = 4.0E0*I_ERI_Gx2yz_Dxz_Pz_S_ac-2.0E0*1*I_ERI_Dxz_Dxz_Pz_S_c;
  abcd[1225] = 4.0E0*I_ERI_Gxy2z_Dxz_Pz_S_ac;
  abcd[1226] = 4.0E0*I_ERI_G4y_Dxz_Pz_S_ac-2.0E0*3*I_ERI_D2y_Dxz_Pz_S_c;
  abcd[1227] = 4.0E0*I_ERI_G3yz_Dxz_Pz_S_ac-2.0E0*2*I_ERI_Dyz_Dxz_Pz_S_c;
  abcd[1228] = 4.0E0*I_ERI_G2y2z_Dxz_Pz_S_ac-2.0E0*1*I_ERI_D2z_Dxz_Pz_S_c;
  abcd[1229] = 4.0E0*I_ERI_Gy3z_Dxz_Pz_S_ac;
  abcd[1230] = 4.0E0*I_ERI_G3xy_D2y_Pz_S_ac;
  abcd[1231] = 4.0E0*I_ERI_G2x2y_D2y_Pz_S_ac-2.0E0*1*I_ERI_D2x_D2y_Pz_S_c;
  abcd[1232] = 4.0E0*I_ERI_G2xyz_D2y_Pz_S_ac;
  abcd[1233] = 4.0E0*I_ERI_Gx3y_D2y_Pz_S_ac-2.0E0*2*I_ERI_Dxy_D2y_Pz_S_c;
  abcd[1234] = 4.0E0*I_ERI_Gx2yz_D2y_Pz_S_ac-2.0E0*1*I_ERI_Dxz_D2y_Pz_S_c;
  abcd[1235] = 4.0E0*I_ERI_Gxy2z_D2y_Pz_S_ac;
  abcd[1236] = 4.0E0*I_ERI_G4y_D2y_Pz_S_ac-2.0E0*3*I_ERI_D2y_D2y_Pz_S_c;
  abcd[1237] = 4.0E0*I_ERI_G3yz_D2y_Pz_S_ac-2.0E0*2*I_ERI_Dyz_D2y_Pz_S_c;
  abcd[1238] = 4.0E0*I_ERI_G2y2z_D2y_Pz_S_ac-2.0E0*1*I_ERI_D2z_D2y_Pz_S_c;
  abcd[1239] = 4.0E0*I_ERI_Gy3z_D2y_Pz_S_ac;
  abcd[1240] = 4.0E0*I_ERI_G3xy_Dyz_Pz_S_ac;
  abcd[1241] = 4.0E0*I_ERI_G2x2y_Dyz_Pz_S_ac-2.0E0*1*I_ERI_D2x_Dyz_Pz_S_c;
  abcd[1242] = 4.0E0*I_ERI_G2xyz_Dyz_Pz_S_ac;
  abcd[1243] = 4.0E0*I_ERI_Gx3y_Dyz_Pz_S_ac-2.0E0*2*I_ERI_Dxy_Dyz_Pz_S_c;
  abcd[1244] = 4.0E0*I_ERI_Gx2yz_Dyz_Pz_S_ac-2.0E0*1*I_ERI_Dxz_Dyz_Pz_S_c;
  abcd[1245] = 4.0E0*I_ERI_Gxy2z_Dyz_Pz_S_ac;
  abcd[1246] = 4.0E0*I_ERI_G4y_Dyz_Pz_S_ac-2.0E0*3*I_ERI_D2y_Dyz_Pz_S_c;
  abcd[1247] = 4.0E0*I_ERI_G3yz_Dyz_Pz_S_ac-2.0E0*2*I_ERI_Dyz_Dyz_Pz_S_c;
  abcd[1248] = 4.0E0*I_ERI_G2y2z_Dyz_Pz_S_ac-2.0E0*1*I_ERI_D2z_Dyz_Pz_S_c;
  abcd[1249] = 4.0E0*I_ERI_Gy3z_Dyz_Pz_S_ac;
  abcd[1250] = 4.0E0*I_ERI_G3xy_D2z_Pz_S_ac;
  abcd[1251] = 4.0E0*I_ERI_G2x2y_D2z_Pz_S_ac-2.0E0*1*I_ERI_D2x_D2z_Pz_S_c;
  abcd[1252] = 4.0E0*I_ERI_G2xyz_D2z_Pz_S_ac;
  abcd[1253] = 4.0E0*I_ERI_Gx3y_D2z_Pz_S_ac-2.0E0*2*I_ERI_Dxy_D2z_Pz_S_c;
  abcd[1254] = 4.0E0*I_ERI_Gx2yz_D2z_Pz_S_ac-2.0E0*1*I_ERI_Dxz_D2z_Pz_S_c;
  abcd[1255] = 4.0E0*I_ERI_Gxy2z_D2z_Pz_S_ac;
  abcd[1256] = 4.0E0*I_ERI_G4y_D2z_Pz_S_ac-2.0E0*3*I_ERI_D2y_D2z_Pz_S_c;
  abcd[1257] = 4.0E0*I_ERI_G3yz_D2z_Pz_S_ac-2.0E0*2*I_ERI_Dyz_D2z_Pz_S_c;
  abcd[1258] = 4.0E0*I_ERI_G2y2z_D2z_Pz_S_ac-2.0E0*1*I_ERI_D2z_D2z_Pz_S_c;
  abcd[1259] = 4.0E0*I_ERI_Gy3z_D2z_Pz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_ac
   * RHS shell quartet name: SQ_ERI_D_D_P_S_c
   ************************************************************/
  abcd[1260] = 4.0E0*I_ERI_G3xz_D2x_Px_S_ac;
  abcd[1261] = 4.0E0*I_ERI_G2xyz_D2x_Px_S_ac;
  abcd[1262] = 4.0E0*I_ERI_G2x2z_D2x_Px_S_ac-2.0E0*1*I_ERI_D2x_D2x_Px_S_c;
  abcd[1263] = 4.0E0*I_ERI_Gx2yz_D2x_Px_S_ac;
  abcd[1264] = 4.0E0*I_ERI_Gxy2z_D2x_Px_S_ac-2.0E0*1*I_ERI_Dxy_D2x_Px_S_c;
  abcd[1265] = 4.0E0*I_ERI_Gx3z_D2x_Px_S_ac-2.0E0*2*I_ERI_Dxz_D2x_Px_S_c;
  abcd[1266] = 4.0E0*I_ERI_G3yz_D2x_Px_S_ac;
  abcd[1267] = 4.0E0*I_ERI_G2y2z_D2x_Px_S_ac-2.0E0*1*I_ERI_D2y_D2x_Px_S_c;
  abcd[1268] = 4.0E0*I_ERI_Gy3z_D2x_Px_S_ac-2.0E0*2*I_ERI_Dyz_D2x_Px_S_c;
  abcd[1269] = 4.0E0*I_ERI_G4z_D2x_Px_S_ac-2.0E0*3*I_ERI_D2z_D2x_Px_S_c;
  abcd[1270] = 4.0E0*I_ERI_G3xz_Dxy_Px_S_ac;
  abcd[1271] = 4.0E0*I_ERI_G2xyz_Dxy_Px_S_ac;
  abcd[1272] = 4.0E0*I_ERI_G2x2z_Dxy_Px_S_ac-2.0E0*1*I_ERI_D2x_Dxy_Px_S_c;
  abcd[1273] = 4.0E0*I_ERI_Gx2yz_Dxy_Px_S_ac;
  abcd[1274] = 4.0E0*I_ERI_Gxy2z_Dxy_Px_S_ac-2.0E0*1*I_ERI_Dxy_Dxy_Px_S_c;
  abcd[1275] = 4.0E0*I_ERI_Gx3z_Dxy_Px_S_ac-2.0E0*2*I_ERI_Dxz_Dxy_Px_S_c;
  abcd[1276] = 4.0E0*I_ERI_G3yz_Dxy_Px_S_ac;
  abcd[1277] = 4.0E0*I_ERI_G2y2z_Dxy_Px_S_ac-2.0E0*1*I_ERI_D2y_Dxy_Px_S_c;
  abcd[1278] = 4.0E0*I_ERI_Gy3z_Dxy_Px_S_ac-2.0E0*2*I_ERI_Dyz_Dxy_Px_S_c;
  abcd[1279] = 4.0E0*I_ERI_G4z_Dxy_Px_S_ac-2.0E0*3*I_ERI_D2z_Dxy_Px_S_c;
  abcd[1280] = 4.0E0*I_ERI_G3xz_Dxz_Px_S_ac;
  abcd[1281] = 4.0E0*I_ERI_G2xyz_Dxz_Px_S_ac;
  abcd[1282] = 4.0E0*I_ERI_G2x2z_Dxz_Px_S_ac-2.0E0*1*I_ERI_D2x_Dxz_Px_S_c;
  abcd[1283] = 4.0E0*I_ERI_Gx2yz_Dxz_Px_S_ac;
  abcd[1284] = 4.0E0*I_ERI_Gxy2z_Dxz_Px_S_ac-2.0E0*1*I_ERI_Dxy_Dxz_Px_S_c;
  abcd[1285] = 4.0E0*I_ERI_Gx3z_Dxz_Px_S_ac-2.0E0*2*I_ERI_Dxz_Dxz_Px_S_c;
  abcd[1286] = 4.0E0*I_ERI_G3yz_Dxz_Px_S_ac;
  abcd[1287] = 4.0E0*I_ERI_G2y2z_Dxz_Px_S_ac-2.0E0*1*I_ERI_D2y_Dxz_Px_S_c;
  abcd[1288] = 4.0E0*I_ERI_Gy3z_Dxz_Px_S_ac-2.0E0*2*I_ERI_Dyz_Dxz_Px_S_c;
  abcd[1289] = 4.0E0*I_ERI_G4z_Dxz_Px_S_ac-2.0E0*3*I_ERI_D2z_Dxz_Px_S_c;
  abcd[1290] = 4.0E0*I_ERI_G3xz_D2y_Px_S_ac;
  abcd[1291] = 4.0E0*I_ERI_G2xyz_D2y_Px_S_ac;
  abcd[1292] = 4.0E0*I_ERI_G2x2z_D2y_Px_S_ac-2.0E0*1*I_ERI_D2x_D2y_Px_S_c;
  abcd[1293] = 4.0E0*I_ERI_Gx2yz_D2y_Px_S_ac;
  abcd[1294] = 4.0E0*I_ERI_Gxy2z_D2y_Px_S_ac-2.0E0*1*I_ERI_Dxy_D2y_Px_S_c;
  abcd[1295] = 4.0E0*I_ERI_Gx3z_D2y_Px_S_ac-2.0E0*2*I_ERI_Dxz_D2y_Px_S_c;
  abcd[1296] = 4.0E0*I_ERI_G3yz_D2y_Px_S_ac;
  abcd[1297] = 4.0E0*I_ERI_G2y2z_D2y_Px_S_ac-2.0E0*1*I_ERI_D2y_D2y_Px_S_c;
  abcd[1298] = 4.0E0*I_ERI_Gy3z_D2y_Px_S_ac-2.0E0*2*I_ERI_Dyz_D2y_Px_S_c;
  abcd[1299] = 4.0E0*I_ERI_G4z_D2y_Px_S_ac-2.0E0*3*I_ERI_D2z_D2y_Px_S_c;
  abcd[1300] = 4.0E0*I_ERI_G3xz_Dyz_Px_S_ac;
  abcd[1301] = 4.0E0*I_ERI_G2xyz_Dyz_Px_S_ac;
  abcd[1302] = 4.0E0*I_ERI_G2x2z_Dyz_Px_S_ac-2.0E0*1*I_ERI_D2x_Dyz_Px_S_c;
  abcd[1303] = 4.0E0*I_ERI_Gx2yz_Dyz_Px_S_ac;
  abcd[1304] = 4.0E0*I_ERI_Gxy2z_Dyz_Px_S_ac-2.0E0*1*I_ERI_Dxy_Dyz_Px_S_c;
  abcd[1305] = 4.0E0*I_ERI_Gx3z_Dyz_Px_S_ac-2.0E0*2*I_ERI_Dxz_Dyz_Px_S_c;
  abcd[1306] = 4.0E0*I_ERI_G3yz_Dyz_Px_S_ac;
  abcd[1307] = 4.0E0*I_ERI_G2y2z_Dyz_Px_S_ac-2.0E0*1*I_ERI_D2y_Dyz_Px_S_c;
  abcd[1308] = 4.0E0*I_ERI_Gy3z_Dyz_Px_S_ac-2.0E0*2*I_ERI_Dyz_Dyz_Px_S_c;
  abcd[1309] = 4.0E0*I_ERI_G4z_Dyz_Px_S_ac-2.0E0*3*I_ERI_D2z_Dyz_Px_S_c;
  abcd[1310] = 4.0E0*I_ERI_G3xz_D2z_Px_S_ac;
  abcd[1311] = 4.0E0*I_ERI_G2xyz_D2z_Px_S_ac;
  abcd[1312] = 4.0E0*I_ERI_G2x2z_D2z_Px_S_ac-2.0E0*1*I_ERI_D2x_D2z_Px_S_c;
  abcd[1313] = 4.0E0*I_ERI_Gx2yz_D2z_Px_S_ac;
  abcd[1314] = 4.0E0*I_ERI_Gxy2z_D2z_Px_S_ac-2.0E0*1*I_ERI_Dxy_D2z_Px_S_c;
  abcd[1315] = 4.0E0*I_ERI_Gx3z_D2z_Px_S_ac-2.0E0*2*I_ERI_Dxz_D2z_Px_S_c;
  abcd[1316] = 4.0E0*I_ERI_G3yz_D2z_Px_S_ac;
  abcd[1317] = 4.0E0*I_ERI_G2y2z_D2z_Px_S_ac-2.0E0*1*I_ERI_D2y_D2z_Px_S_c;
  abcd[1318] = 4.0E0*I_ERI_Gy3z_D2z_Px_S_ac-2.0E0*2*I_ERI_Dyz_D2z_Px_S_c;
  abcd[1319] = 4.0E0*I_ERI_G4z_D2z_Px_S_ac-2.0E0*3*I_ERI_D2z_D2z_Px_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_ac
   * RHS shell quartet name: SQ_ERI_D_D_P_S_c
   ************************************************************/
  abcd[1320] = 4.0E0*I_ERI_G3xz_D2x_Py_S_ac;
  abcd[1321] = 4.0E0*I_ERI_G2xyz_D2x_Py_S_ac;
  abcd[1322] = 4.0E0*I_ERI_G2x2z_D2x_Py_S_ac-2.0E0*1*I_ERI_D2x_D2x_Py_S_c;
  abcd[1323] = 4.0E0*I_ERI_Gx2yz_D2x_Py_S_ac;
  abcd[1324] = 4.0E0*I_ERI_Gxy2z_D2x_Py_S_ac-2.0E0*1*I_ERI_Dxy_D2x_Py_S_c;
  abcd[1325] = 4.0E0*I_ERI_Gx3z_D2x_Py_S_ac-2.0E0*2*I_ERI_Dxz_D2x_Py_S_c;
  abcd[1326] = 4.0E0*I_ERI_G3yz_D2x_Py_S_ac;
  abcd[1327] = 4.0E0*I_ERI_G2y2z_D2x_Py_S_ac-2.0E0*1*I_ERI_D2y_D2x_Py_S_c;
  abcd[1328] = 4.0E0*I_ERI_Gy3z_D2x_Py_S_ac-2.0E0*2*I_ERI_Dyz_D2x_Py_S_c;
  abcd[1329] = 4.0E0*I_ERI_G4z_D2x_Py_S_ac-2.0E0*3*I_ERI_D2z_D2x_Py_S_c;
  abcd[1330] = 4.0E0*I_ERI_G3xz_Dxy_Py_S_ac;
  abcd[1331] = 4.0E0*I_ERI_G2xyz_Dxy_Py_S_ac;
  abcd[1332] = 4.0E0*I_ERI_G2x2z_Dxy_Py_S_ac-2.0E0*1*I_ERI_D2x_Dxy_Py_S_c;
  abcd[1333] = 4.0E0*I_ERI_Gx2yz_Dxy_Py_S_ac;
  abcd[1334] = 4.0E0*I_ERI_Gxy2z_Dxy_Py_S_ac-2.0E0*1*I_ERI_Dxy_Dxy_Py_S_c;
  abcd[1335] = 4.0E0*I_ERI_Gx3z_Dxy_Py_S_ac-2.0E0*2*I_ERI_Dxz_Dxy_Py_S_c;
  abcd[1336] = 4.0E0*I_ERI_G3yz_Dxy_Py_S_ac;
  abcd[1337] = 4.0E0*I_ERI_G2y2z_Dxy_Py_S_ac-2.0E0*1*I_ERI_D2y_Dxy_Py_S_c;
  abcd[1338] = 4.0E0*I_ERI_Gy3z_Dxy_Py_S_ac-2.0E0*2*I_ERI_Dyz_Dxy_Py_S_c;
  abcd[1339] = 4.0E0*I_ERI_G4z_Dxy_Py_S_ac-2.0E0*3*I_ERI_D2z_Dxy_Py_S_c;
  abcd[1340] = 4.0E0*I_ERI_G3xz_Dxz_Py_S_ac;
  abcd[1341] = 4.0E0*I_ERI_G2xyz_Dxz_Py_S_ac;
  abcd[1342] = 4.0E0*I_ERI_G2x2z_Dxz_Py_S_ac-2.0E0*1*I_ERI_D2x_Dxz_Py_S_c;
  abcd[1343] = 4.0E0*I_ERI_Gx2yz_Dxz_Py_S_ac;
  abcd[1344] = 4.0E0*I_ERI_Gxy2z_Dxz_Py_S_ac-2.0E0*1*I_ERI_Dxy_Dxz_Py_S_c;
  abcd[1345] = 4.0E0*I_ERI_Gx3z_Dxz_Py_S_ac-2.0E0*2*I_ERI_Dxz_Dxz_Py_S_c;
  abcd[1346] = 4.0E0*I_ERI_G3yz_Dxz_Py_S_ac;
  abcd[1347] = 4.0E0*I_ERI_G2y2z_Dxz_Py_S_ac-2.0E0*1*I_ERI_D2y_Dxz_Py_S_c;
  abcd[1348] = 4.0E0*I_ERI_Gy3z_Dxz_Py_S_ac-2.0E0*2*I_ERI_Dyz_Dxz_Py_S_c;
  abcd[1349] = 4.0E0*I_ERI_G4z_Dxz_Py_S_ac-2.0E0*3*I_ERI_D2z_Dxz_Py_S_c;
  abcd[1350] = 4.0E0*I_ERI_G3xz_D2y_Py_S_ac;
  abcd[1351] = 4.0E0*I_ERI_G2xyz_D2y_Py_S_ac;
  abcd[1352] = 4.0E0*I_ERI_G2x2z_D2y_Py_S_ac-2.0E0*1*I_ERI_D2x_D2y_Py_S_c;
  abcd[1353] = 4.0E0*I_ERI_Gx2yz_D2y_Py_S_ac;
  abcd[1354] = 4.0E0*I_ERI_Gxy2z_D2y_Py_S_ac-2.0E0*1*I_ERI_Dxy_D2y_Py_S_c;
  abcd[1355] = 4.0E0*I_ERI_Gx3z_D2y_Py_S_ac-2.0E0*2*I_ERI_Dxz_D2y_Py_S_c;
  abcd[1356] = 4.0E0*I_ERI_G3yz_D2y_Py_S_ac;
  abcd[1357] = 4.0E0*I_ERI_G2y2z_D2y_Py_S_ac-2.0E0*1*I_ERI_D2y_D2y_Py_S_c;
  abcd[1358] = 4.0E0*I_ERI_Gy3z_D2y_Py_S_ac-2.0E0*2*I_ERI_Dyz_D2y_Py_S_c;
  abcd[1359] = 4.0E0*I_ERI_G4z_D2y_Py_S_ac-2.0E0*3*I_ERI_D2z_D2y_Py_S_c;
  abcd[1360] = 4.0E0*I_ERI_G3xz_Dyz_Py_S_ac;
  abcd[1361] = 4.0E0*I_ERI_G2xyz_Dyz_Py_S_ac;
  abcd[1362] = 4.0E0*I_ERI_G2x2z_Dyz_Py_S_ac-2.0E0*1*I_ERI_D2x_Dyz_Py_S_c;
  abcd[1363] = 4.0E0*I_ERI_Gx2yz_Dyz_Py_S_ac;
  abcd[1364] = 4.0E0*I_ERI_Gxy2z_Dyz_Py_S_ac-2.0E0*1*I_ERI_Dxy_Dyz_Py_S_c;
  abcd[1365] = 4.0E0*I_ERI_Gx3z_Dyz_Py_S_ac-2.0E0*2*I_ERI_Dxz_Dyz_Py_S_c;
  abcd[1366] = 4.0E0*I_ERI_G3yz_Dyz_Py_S_ac;
  abcd[1367] = 4.0E0*I_ERI_G2y2z_Dyz_Py_S_ac-2.0E0*1*I_ERI_D2y_Dyz_Py_S_c;
  abcd[1368] = 4.0E0*I_ERI_Gy3z_Dyz_Py_S_ac-2.0E0*2*I_ERI_Dyz_Dyz_Py_S_c;
  abcd[1369] = 4.0E0*I_ERI_G4z_Dyz_Py_S_ac-2.0E0*3*I_ERI_D2z_Dyz_Py_S_c;
  abcd[1370] = 4.0E0*I_ERI_G3xz_D2z_Py_S_ac;
  abcd[1371] = 4.0E0*I_ERI_G2xyz_D2z_Py_S_ac;
  abcd[1372] = 4.0E0*I_ERI_G2x2z_D2z_Py_S_ac-2.0E0*1*I_ERI_D2x_D2z_Py_S_c;
  abcd[1373] = 4.0E0*I_ERI_Gx2yz_D2z_Py_S_ac;
  abcd[1374] = 4.0E0*I_ERI_Gxy2z_D2z_Py_S_ac-2.0E0*1*I_ERI_Dxy_D2z_Py_S_c;
  abcd[1375] = 4.0E0*I_ERI_Gx3z_D2z_Py_S_ac-2.0E0*2*I_ERI_Dxz_D2z_Py_S_c;
  abcd[1376] = 4.0E0*I_ERI_G3yz_D2z_Py_S_ac;
  abcd[1377] = 4.0E0*I_ERI_G2y2z_D2z_Py_S_ac-2.0E0*1*I_ERI_D2y_D2z_Py_S_c;
  abcd[1378] = 4.0E0*I_ERI_Gy3z_D2z_Py_S_ac-2.0E0*2*I_ERI_Dyz_D2z_Py_S_c;
  abcd[1379] = 4.0E0*I_ERI_G4z_D2z_Py_S_ac-2.0E0*3*I_ERI_D2z_D2z_Py_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_P_S_ac
   * RHS shell quartet name: SQ_ERI_D_D_P_S_c
   ************************************************************/
  abcd[1380] = 4.0E0*I_ERI_G3xz_D2x_Pz_S_ac;
  abcd[1381] = 4.0E0*I_ERI_G2xyz_D2x_Pz_S_ac;
  abcd[1382] = 4.0E0*I_ERI_G2x2z_D2x_Pz_S_ac-2.0E0*1*I_ERI_D2x_D2x_Pz_S_c;
  abcd[1383] = 4.0E0*I_ERI_Gx2yz_D2x_Pz_S_ac;
  abcd[1384] = 4.0E0*I_ERI_Gxy2z_D2x_Pz_S_ac-2.0E0*1*I_ERI_Dxy_D2x_Pz_S_c;
  abcd[1385] = 4.0E0*I_ERI_Gx3z_D2x_Pz_S_ac-2.0E0*2*I_ERI_Dxz_D2x_Pz_S_c;
  abcd[1386] = 4.0E0*I_ERI_G3yz_D2x_Pz_S_ac;
  abcd[1387] = 4.0E0*I_ERI_G2y2z_D2x_Pz_S_ac-2.0E0*1*I_ERI_D2y_D2x_Pz_S_c;
  abcd[1388] = 4.0E0*I_ERI_Gy3z_D2x_Pz_S_ac-2.0E0*2*I_ERI_Dyz_D2x_Pz_S_c;
  abcd[1389] = 4.0E0*I_ERI_G4z_D2x_Pz_S_ac-2.0E0*3*I_ERI_D2z_D2x_Pz_S_c;
  abcd[1390] = 4.0E0*I_ERI_G3xz_Dxy_Pz_S_ac;
  abcd[1391] = 4.0E0*I_ERI_G2xyz_Dxy_Pz_S_ac;
  abcd[1392] = 4.0E0*I_ERI_G2x2z_Dxy_Pz_S_ac-2.0E0*1*I_ERI_D2x_Dxy_Pz_S_c;
  abcd[1393] = 4.0E0*I_ERI_Gx2yz_Dxy_Pz_S_ac;
  abcd[1394] = 4.0E0*I_ERI_Gxy2z_Dxy_Pz_S_ac-2.0E0*1*I_ERI_Dxy_Dxy_Pz_S_c;
  abcd[1395] = 4.0E0*I_ERI_Gx3z_Dxy_Pz_S_ac-2.0E0*2*I_ERI_Dxz_Dxy_Pz_S_c;
  abcd[1396] = 4.0E0*I_ERI_G3yz_Dxy_Pz_S_ac;
  abcd[1397] = 4.0E0*I_ERI_G2y2z_Dxy_Pz_S_ac-2.0E0*1*I_ERI_D2y_Dxy_Pz_S_c;
  abcd[1398] = 4.0E0*I_ERI_Gy3z_Dxy_Pz_S_ac-2.0E0*2*I_ERI_Dyz_Dxy_Pz_S_c;
  abcd[1399] = 4.0E0*I_ERI_G4z_Dxy_Pz_S_ac-2.0E0*3*I_ERI_D2z_Dxy_Pz_S_c;
  abcd[1400] = 4.0E0*I_ERI_G3xz_Dxz_Pz_S_ac;
  abcd[1401] = 4.0E0*I_ERI_G2xyz_Dxz_Pz_S_ac;
  abcd[1402] = 4.0E0*I_ERI_G2x2z_Dxz_Pz_S_ac-2.0E0*1*I_ERI_D2x_Dxz_Pz_S_c;
  abcd[1403] = 4.0E0*I_ERI_Gx2yz_Dxz_Pz_S_ac;
  abcd[1404] = 4.0E0*I_ERI_Gxy2z_Dxz_Pz_S_ac-2.0E0*1*I_ERI_Dxy_Dxz_Pz_S_c;
  abcd[1405] = 4.0E0*I_ERI_Gx3z_Dxz_Pz_S_ac-2.0E0*2*I_ERI_Dxz_Dxz_Pz_S_c;
  abcd[1406] = 4.0E0*I_ERI_G3yz_Dxz_Pz_S_ac;
  abcd[1407] = 4.0E0*I_ERI_G2y2z_Dxz_Pz_S_ac-2.0E0*1*I_ERI_D2y_Dxz_Pz_S_c;
  abcd[1408] = 4.0E0*I_ERI_Gy3z_Dxz_Pz_S_ac-2.0E0*2*I_ERI_Dyz_Dxz_Pz_S_c;
  abcd[1409] = 4.0E0*I_ERI_G4z_Dxz_Pz_S_ac-2.0E0*3*I_ERI_D2z_Dxz_Pz_S_c;
  abcd[1410] = 4.0E0*I_ERI_G3xz_D2y_Pz_S_ac;
  abcd[1411] = 4.0E0*I_ERI_G2xyz_D2y_Pz_S_ac;
  abcd[1412] = 4.0E0*I_ERI_G2x2z_D2y_Pz_S_ac-2.0E0*1*I_ERI_D2x_D2y_Pz_S_c;
  abcd[1413] = 4.0E0*I_ERI_Gx2yz_D2y_Pz_S_ac;
  abcd[1414] = 4.0E0*I_ERI_Gxy2z_D2y_Pz_S_ac-2.0E0*1*I_ERI_Dxy_D2y_Pz_S_c;
  abcd[1415] = 4.0E0*I_ERI_Gx3z_D2y_Pz_S_ac-2.0E0*2*I_ERI_Dxz_D2y_Pz_S_c;
  abcd[1416] = 4.0E0*I_ERI_G3yz_D2y_Pz_S_ac;
  abcd[1417] = 4.0E0*I_ERI_G2y2z_D2y_Pz_S_ac-2.0E0*1*I_ERI_D2y_D2y_Pz_S_c;
  abcd[1418] = 4.0E0*I_ERI_Gy3z_D2y_Pz_S_ac-2.0E0*2*I_ERI_Dyz_D2y_Pz_S_c;
  abcd[1419] = 4.0E0*I_ERI_G4z_D2y_Pz_S_ac-2.0E0*3*I_ERI_D2z_D2y_Pz_S_c;
  abcd[1420] = 4.0E0*I_ERI_G3xz_Dyz_Pz_S_ac;
  abcd[1421] = 4.0E0*I_ERI_G2xyz_Dyz_Pz_S_ac;
  abcd[1422] = 4.0E0*I_ERI_G2x2z_Dyz_Pz_S_ac-2.0E0*1*I_ERI_D2x_Dyz_Pz_S_c;
  abcd[1423] = 4.0E0*I_ERI_Gx2yz_Dyz_Pz_S_ac;
  abcd[1424] = 4.0E0*I_ERI_Gxy2z_Dyz_Pz_S_ac-2.0E0*1*I_ERI_Dxy_Dyz_Pz_S_c;
  abcd[1425] = 4.0E0*I_ERI_Gx3z_Dyz_Pz_S_ac-2.0E0*2*I_ERI_Dxz_Dyz_Pz_S_c;
  abcd[1426] = 4.0E0*I_ERI_G3yz_Dyz_Pz_S_ac;
  abcd[1427] = 4.0E0*I_ERI_G2y2z_Dyz_Pz_S_ac-2.0E0*1*I_ERI_D2y_Dyz_Pz_S_c;
  abcd[1428] = 4.0E0*I_ERI_Gy3z_Dyz_Pz_S_ac-2.0E0*2*I_ERI_Dyz_Dyz_Pz_S_c;
  abcd[1429] = 4.0E0*I_ERI_G4z_Dyz_Pz_S_ac-2.0E0*3*I_ERI_D2z_Dyz_Pz_S_c;
  abcd[1430] = 4.0E0*I_ERI_G3xz_D2z_Pz_S_ac;
  abcd[1431] = 4.0E0*I_ERI_G2xyz_D2z_Pz_S_ac;
  abcd[1432] = 4.0E0*I_ERI_G2x2z_D2z_Pz_S_ac-2.0E0*1*I_ERI_D2x_D2z_Pz_S_c;
  abcd[1433] = 4.0E0*I_ERI_Gx2yz_D2z_Pz_S_ac;
  abcd[1434] = 4.0E0*I_ERI_Gxy2z_D2z_Pz_S_ac-2.0E0*1*I_ERI_Dxy_D2z_Pz_S_c;
  abcd[1435] = 4.0E0*I_ERI_Gx3z_D2z_Pz_S_ac-2.0E0*2*I_ERI_Dxz_D2z_Pz_S_c;
  abcd[1436] = 4.0E0*I_ERI_G3yz_D2z_Pz_S_ac;
  abcd[1437] = 4.0E0*I_ERI_G2y2z_D2z_Pz_S_ac-2.0E0*1*I_ERI_D2y_D2z_Pz_S_c;
  abcd[1438] = 4.0E0*I_ERI_Gy3z_D2z_Pz_S_ac-2.0E0*2*I_ERI_Dyz_D2z_Pz_S_c;
  abcd[1439] = 4.0E0*I_ERI_G4z_D2z_Pz_S_ac-2.0E0*3*I_ERI_D2z_D2z_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_G_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1440] = 4.0E0*I_ERI_F3x_G4x_S_S_bb-2.0E0*2*I_ERI_F3x_D2x_S_S_b-2.0E0*3*I_ERI_F3x_D2x_S_S_b+2*1*I_ERI_F3x_S_S_S;
  abcd[1441] = 4.0E0*I_ERI_F2xy_G4x_S_S_bb-2.0E0*2*I_ERI_F2xy_D2x_S_S_b-2.0E0*3*I_ERI_F2xy_D2x_S_S_b+2*1*I_ERI_F2xy_S_S_S;
  abcd[1442] = 4.0E0*I_ERI_F2xz_G4x_S_S_bb-2.0E0*2*I_ERI_F2xz_D2x_S_S_b-2.0E0*3*I_ERI_F2xz_D2x_S_S_b+2*1*I_ERI_F2xz_S_S_S;
  abcd[1443] = 4.0E0*I_ERI_Fx2y_G4x_S_S_bb-2.0E0*2*I_ERI_Fx2y_D2x_S_S_b-2.0E0*3*I_ERI_Fx2y_D2x_S_S_b+2*1*I_ERI_Fx2y_S_S_S;
  abcd[1444] = 4.0E0*I_ERI_Fxyz_G4x_S_S_bb-2.0E0*2*I_ERI_Fxyz_D2x_S_S_b-2.0E0*3*I_ERI_Fxyz_D2x_S_S_b+2*1*I_ERI_Fxyz_S_S_S;
  abcd[1445] = 4.0E0*I_ERI_Fx2z_G4x_S_S_bb-2.0E0*2*I_ERI_Fx2z_D2x_S_S_b-2.0E0*3*I_ERI_Fx2z_D2x_S_S_b+2*1*I_ERI_Fx2z_S_S_S;
  abcd[1446] = 4.0E0*I_ERI_F3y_G4x_S_S_bb-2.0E0*2*I_ERI_F3y_D2x_S_S_b-2.0E0*3*I_ERI_F3y_D2x_S_S_b+2*1*I_ERI_F3y_S_S_S;
  abcd[1447] = 4.0E0*I_ERI_F2yz_G4x_S_S_bb-2.0E0*2*I_ERI_F2yz_D2x_S_S_b-2.0E0*3*I_ERI_F2yz_D2x_S_S_b+2*1*I_ERI_F2yz_S_S_S;
  abcd[1448] = 4.0E0*I_ERI_Fy2z_G4x_S_S_bb-2.0E0*2*I_ERI_Fy2z_D2x_S_S_b-2.0E0*3*I_ERI_Fy2z_D2x_S_S_b+2*1*I_ERI_Fy2z_S_S_S;
  abcd[1449] = 4.0E0*I_ERI_F3z_G4x_S_S_bb-2.0E0*2*I_ERI_F3z_D2x_S_S_b-2.0E0*3*I_ERI_F3z_D2x_S_S_b+2*1*I_ERI_F3z_S_S_S;
  abcd[1450] = 4.0E0*I_ERI_F3x_G3xy_S_S_bb-2.0E0*1*I_ERI_F3x_Dxy_S_S_b-2.0E0*2*I_ERI_F3x_Dxy_S_S_b;
  abcd[1451] = 4.0E0*I_ERI_F2xy_G3xy_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxy_S_S_b-2.0E0*2*I_ERI_F2xy_Dxy_S_S_b;
  abcd[1452] = 4.0E0*I_ERI_F2xz_G3xy_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxy_S_S_b-2.0E0*2*I_ERI_F2xz_Dxy_S_S_b;
  abcd[1453] = 4.0E0*I_ERI_Fx2y_G3xy_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_b-2.0E0*2*I_ERI_Fx2y_Dxy_S_S_b;
  abcd[1454] = 4.0E0*I_ERI_Fxyz_G3xy_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_b-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_b;
  abcd[1455] = 4.0E0*I_ERI_Fx2z_G3xy_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_b-2.0E0*2*I_ERI_Fx2z_Dxy_S_S_b;
  abcd[1456] = 4.0E0*I_ERI_F3y_G3xy_S_S_bb-2.0E0*1*I_ERI_F3y_Dxy_S_S_b-2.0E0*2*I_ERI_F3y_Dxy_S_S_b;
  abcd[1457] = 4.0E0*I_ERI_F2yz_G3xy_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxy_S_S_b-2.0E0*2*I_ERI_F2yz_Dxy_S_S_b;
  abcd[1458] = 4.0E0*I_ERI_Fy2z_G3xy_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_b-2.0E0*2*I_ERI_Fy2z_Dxy_S_S_b;
  abcd[1459] = 4.0E0*I_ERI_F3z_G3xy_S_S_bb-2.0E0*1*I_ERI_F3z_Dxy_S_S_b-2.0E0*2*I_ERI_F3z_Dxy_S_S_b;
  abcd[1460] = 4.0E0*I_ERI_F3x_G3xz_S_S_bb-2.0E0*1*I_ERI_F3x_Dxz_S_S_b-2.0E0*2*I_ERI_F3x_Dxz_S_S_b;
  abcd[1461] = 4.0E0*I_ERI_F2xy_G3xz_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxz_S_S_b-2.0E0*2*I_ERI_F2xy_Dxz_S_S_b;
  abcd[1462] = 4.0E0*I_ERI_F2xz_G3xz_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxz_S_S_b-2.0E0*2*I_ERI_F2xz_Dxz_S_S_b;
  abcd[1463] = 4.0E0*I_ERI_Fx2y_G3xz_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_b-2.0E0*2*I_ERI_Fx2y_Dxz_S_S_b;
  abcd[1464] = 4.0E0*I_ERI_Fxyz_G3xz_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_b-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_b;
  abcd[1465] = 4.0E0*I_ERI_Fx2z_G3xz_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_b-2.0E0*2*I_ERI_Fx2z_Dxz_S_S_b;
  abcd[1466] = 4.0E0*I_ERI_F3y_G3xz_S_S_bb-2.0E0*1*I_ERI_F3y_Dxz_S_S_b-2.0E0*2*I_ERI_F3y_Dxz_S_S_b;
  abcd[1467] = 4.0E0*I_ERI_F2yz_G3xz_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxz_S_S_b-2.0E0*2*I_ERI_F2yz_Dxz_S_S_b;
  abcd[1468] = 4.0E0*I_ERI_Fy2z_G3xz_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_b-2.0E0*2*I_ERI_Fy2z_Dxz_S_S_b;
  abcd[1469] = 4.0E0*I_ERI_F3z_G3xz_S_S_bb-2.0E0*1*I_ERI_F3z_Dxz_S_S_b-2.0E0*2*I_ERI_F3z_Dxz_S_S_b;
  abcd[1470] = 4.0E0*I_ERI_F3x_G2x2y_S_S_bb-2.0E0*1*I_ERI_F3x_D2y_S_S_b;
  abcd[1471] = 4.0E0*I_ERI_F2xy_G2x2y_S_S_bb-2.0E0*1*I_ERI_F2xy_D2y_S_S_b;
  abcd[1472] = 4.0E0*I_ERI_F2xz_G2x2y_S_S_bb-2.0E0*1*I_ERI_F2xz_D2y_S_S_b;
  abcd[1473] = 4.0E0*I_ERI_Fx2y_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fx2y_D2y_S_S_b;
  abcd[1474] = 4.0E0*I_ERI_Fxyz_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fxyz_D2y_S_S_b;
  abcd[1475] = 4.0E0*I_ERI_Fx2z_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fx2z_D2y_S_S_b;
  abcd[1476] = 4.0E0*I_ERI_F3y_G2x2y_S_S_bb-2.0E0*1*I_ERI_F3y_D2y_S_S_b;
  abcd[1477] = 4.0E0*I_ERI_F2yz_G2x2y_S_S_bb-2.0E0*1*I_ERI_F2yz_D2y_S_S_b;
  abcd[1478] = 4.0E0*I_ERI_Fy2z_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fy2z_D2y_S_S_b;
  abcd[1479] = 4.0E0*I_ERI_F3z_G2x2y_S_S_bb-2.0E0*1*I_ERI_F3z_D2y_S_S_b;
  abcd[1480] = 4.0E0*I_ERI_F3x_G2xyz_S_S_bb-2.0E0*1*I_ERI_F3x_Dyz_S_S_b;
  abcd[1481] = 4.0E0*I_ERI_F2xy_G2xyz_S_S_bb-2.0E0*1*I_ERI_F2xy_Dyz_S_S_b;
  abcd[1482] = 4.0E0*I_ERI_F2xz_G2xyz_S_S_bb-2.0E0*1*I_ERI_F2xz_Dyz_S_S_b;
  abcd[1483] = 4.0E0*I_ERI_Fx2y_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_b;
  abcd[1484] = 4.0E0*I_ERI_Fxyz_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_b;
  abcd[1485] = 4.0E0*I_ERI_Fx2z_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_b;
  abcd[1486] = 4.0E0*I_ERI_F3y_G2xyz_S_S_bb-2.0E0*1*I_ERI_F3y_Dyz_S_S_b;
  abcd[1487] = 4.0E0*I_ERI_F2yz_G2xyz_S_S_bb-2.0E0*1*I_ERI_F2yz_Dyz_S_S_b;
  abcd[1488] = 4.0E0*I_ERI_Fy2z_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_b;
  abcd[1489] = 4.0E0*I_ERI_F3z_G2xyz_S_S_bb-2.0E0*1*I_ERI_F3z_Dyz_S_S_b;
  abcd[1490] = 4.0E0*I_ERI_F3x_G2x2z_S_S_bb-2.0E0*1*I_ERI_F3x_D2z_S_S_b;
  abcd[1491] = 4.0E0*I_ERI_F2xy_G2x2z_S_S_bb-2.0E0*1*I_ERI_F2xy_D2z_S_S_b;
  abcd[1492] = 4.0E0*I_ERI_F2xz_G2x2z_S_S_bb-2.0E0*1*I_ERI_F2xz_D2z_S_S_b;
  abcd[1493] = 4.0E0*I_ERI_Fx2y_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fx2y_D2z_S_S_b;
  abcd[1494] = 4.0E0*I_ERI_Fxyz_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fxyz_D2z_S_S_b;
  abcd[1495] = 4.0E0*I_ERI_Fx2z_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fx2z_D2z_S_S_b;
  abcd[1496] = 4.0E0*I_ERI_F3y_G2x2z_S_S_bb-2.0E0*1*I_ERI_F3y_D2z_S_S_b;
  abcd[1497] = 4.0E0*I_ERI_F2yz_G2x2z_S_S_bb-2.0E0*1*I_ERI_F2yz_D2z_S_S_b;
  abcd[1498] = 4.0E0*I_ERI_Fy2z_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fy2z_D2z_S_S_b;
  abcd[1499] = 4.0E0*I_ERI_F3z_G2x2z_S_S_bb-2.0E0*1*I_ERI_F3z_D2z_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_G_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1500] = 4.0E0*I_ERI_F3x_G3xy_S_S_bb-2.0E0*2*I_ERI_F3x_Dxy_S_S_b;
  abcd[1501] = 4.0E0*I_ERI_F2xy_G3xy_S_S_bb-2.0E0*2*I_ERI_F2xy_Dxy_S_S_b;
  abcd[1502] = 4.0E0*I_ERI_F2xz_G3xy_S_S_bb-2.0E0*2*I_ERI_F2xz_Dxy_S_S_b;
  abcd[1503] = 4.0E0*I_ERI_Fx2y_G3xy_S_S_bb-2.0E0*2*I_ERI_Fx2y_Dxy_S_S_b;
  abcd[1504] = 4.0E0*I_ERI_Fxyz_G3xy_S_S_bb-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_b;
  abcd[1505] = 4.0E0*I_ERI_Fx2z_G3xy_S_S_bb-2.0E0*2*I_ERI_Fx2z_Dxy_S_S_b;
  abcd[1506] = 4.0E0*I_ERI_F3y_G3xy_S_S_bb-2.0E0*2*I_ERI_F3y_Dxy_S_S_b;
  abcd[1507] = 4.0E0*I_ERI_F2yz_G3xy_S_S_bb-2.0E0*2*I_ERI_F2yz_Dxy_S_S_b;
  abcd[1508] = 4.0E0*I_ERI_Fy2z_G3xy_S_S_bb-2.0E0*2*I_ERI_Fy2z_Dxy_S_S_b;
  abcd[1509] = 4.0E0*I_ERI_F3z_G3xy_S_S_bb-2.0E0*2*I_ERI_F3z_Dxy_S_S_b;
  abcd[1510] = 4.0E0*I_ERI_F3x_G2x2y_S_S_bb-2.0E0*1*I_ERI_F3x_D2x_S_S_b-2.0E0*1*I_ERI_F3x_D2y_S_S_b+1*I_ERI_F3x_S_S_S;
  abcd[1511] = 4.0E0*I_ERI_F2xy_G2x2y_S_S_bb-2.0E0*1*I_ERI_F2xy_D2x_S_S_b-2.0E0*1*I_ERI_F2xy_D2y_S_S_b+1*I_ERI_F2xy_S_S_S;
  abcd[1512] = 4.0E0*I_ERI_F2xz_G2x2y_S_S_bb-2.0E0*1*I_ERI_F2xz_D2x_S_S_b-2.0E0*1*I_ERI_F2xz_D2y_S_S_b+1*I_ERI_F2xz_S_S_S;
  abcd[1513] = 4.0E0*I_ERI_Fx2y_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fx2y_D2x_S_S_b-2.0E0*1*I_ERI_Fx2y_D2y_S_S_b+1*I_ERI_Fx2y_S_S_S;
  abcd[1514] = 4.0E0*I_ERI_Fxyz_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fxyz_D2x_S_S_b-2.0E0*1*I_ERI_Fxyz_D2y_S_S_b+1*I_ERI_Fxyz_S_S_S;
  abcd[1515] = 4.0E0*I_ERI_Fx2z_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fx2z_D2x_S_S_b-2.0E0*1*I_ERI_Fx2z_D2y_S_S_b+1*I_ERI_Fx2z_S_S_S;
  abcd[1516] = 4.0E0*I_ERI_F3y_G2x2y_S_S_bb-2.0E0*1*I_ERI_F3y_D2x_S_S_b-2.0E0*1*I_ERI_F3y_D2y_S_S_b+1*I_ERI_F3y_S_S_S;
  abcd[1517] = 4.0E0*I_ERI_F2yz_G2x2y_S_S_bb-2.0E0*1*I_ERI_F2yz_D2x_S_S_b-2.0E0*1*I_ERI_F2yz_D2y_S_S_b+1*I_ERI_F2yz_S_S_S;
  abcd[1518] = 4.0E0*I_ERI_Fy2z_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fy2z_D2x_S_S_b-2.0E0*1*I_ERI_Fy2z_D2y_S_S_b+1*I_ERI_Fy2z_S_S_S;
  abcd[1519] = 4.0E0*I_ERI_F3z_G2x2y_S_S_bb-2.0E0*1*I_ERI_F3z_D2x_S_S_b-2.0E0*1*I_ERI_F3z_D2y_S_S_b+1*I_ERI_F3z_S_S_S;
  abcd[1520] = 4.0E0*I_ERI_F3x_G2xyz_S_S_bb-2.0E0*1*I_ERI_F3x_Dyz_S_S_b;
  abcd[1521] = 4.0E0*I_ERI_F2xy_G2xyz_S_S_bb-2.0E0*1*I_ERI_F2xy_Dyz_S_S_b;
  abcd[1522] = 4.0E0*I_ERI_F2xz_G2xyz_S_S_bb-2.0E0*1*I_ERI_F2xz_Dyz_S_S_b;
  abcd[1523] = 4.0E0*I_ERI_Fx2y_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_b;
  abcd[1524] = 4.0E0*I_ERI_Fxyz_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_b;
  abcd[1525] = 4.0E0*I_ERI_Fx2z_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_b;
  abcd[1526] = 4.0E0*I_ERI_F3y_G2xyz_S_S_bb-2.0E0*1*I_ERI_F3y_Dyz_S_S_b;
  abcd[1527] = 4.0E0*I_ERI_F2yz_G2xyz_S_S_bb-2.0E0*1*I_ERI_F2yz_Dyz_S_S_b;
  abcd[1528] = 4.0E0*I_ERI_Fy2z_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_b;
  abcd[1529] = 4.0E0*I_ERI_F3z_G2xyz_S_S_bb-2.0E0*1*I_ERI_F3z_Dyz_S_S_b;
  abcd[1530] = 4.0E0*I_ERI_F3x_Gx3y_S_S_bb-2.0E0*2*I_ERI_F3x_Dxy_S_S_b;
  abcd[1531] = 4.0E0*I_ERI_F2xy_Gx3y_S_S_bb-2.0E0*2*I_ERI_F2xy_Dxy_S_S_b;
  abcd[1532] = 4.0E0*I_ERI_F2xz_Gx3y_S_S_bb-2.0E0*2*I_ERI_F2xz_Dxy_S_S_b;
  abcd[1533] = 4.0E0*I_ERI_Fx2y_Gx3y_S_S_bb-2.0E0*2*I_ERI_Fx2y_Dxy_S_S_b;
  abcd[1534] = 4.0E0*I_ERI_Fxyz_Gx3y_S_S_bb-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_b;
  abcd[1535] = 4.0E0*I_ERI_Fx2z_Gx3y_S_S_bb-2.0E0*2*I_ERI_Fx2z_Dxy_S_S_b;
  abcd[1536] = 4.0E0*I_ERI_F3y_Gx3y_S_S_bb-2.0E0*2*I_ERI_F3y_Dxy_S_S_b;
  abcd[1537] = 4.0E0*I_ERI_F2yz_Gx3y_S_S_bb-2.0E0*2*I_ERI_F2yz_Dxy_S_S_b;
  abcd[1538] = 4.0E0*I_ERI_Fy2z_Gx3y_S_S_bb-2.0E0*2*I_ERI_Fy2z_Dxy_S_S_b;
  abcd[1539] = 4.0E0*I_ERI_F3z_Gx3y_S_S_bb-2.0E0*2*I_ERI_F3z_Dxy_S_S_b;
  abcd[1540] = 4.0E0*I_ERI_F3x_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F3x_Dxz_S_S_b;
  abcd[1541] = 4.0E0*I_ERI_F2xy_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxz_S_S_b;
  abcd[1542] = 4.0E0*I_ERI_F2xz_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxz_S_S_b;
  abcd[1543] = 4.0E0*I_ERI_Fx2y_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_b;
  abcd[1544] = 4.0E0*I_ERI_Fxyz_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_b;
  abcd[1545] = 4.0E0*I_ERI_Fx2z_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_b;
  abcd[1546] = 4.0E0*I_ERI_F3y_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F3y_Dxz_S_S_b;
  abcd[1547] = 4.0E0*I_ERI_F2yz_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxz_S_S_b;
  abcd[1548] = 4.0E0*I_ERI_Fy2z_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_b;
  abcd[1549] = 4.0E0*I_ERI_F3z_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F3z_Dxz_S_S_b;
  abcd[1550] = 4.0E0*I_ERI_F3x_Gxy2z_S_S_bb;
  abcd[1551] = 4.0E0*I_ERI_F2xy_Gxy2z_S_S_bb;
  abcd[1552] = 4.0E0*I_ERI_F2xz_Gxy2z_S_S_bb;
  abcd[1553] = 4.0E0*I_ERI_Fx2y_Gxy2z_S_S_bb;
  abcd[1554] = 4.0E0*I_ERI_Fxyz_Gxy2z_S_S_bb;
  abcd[1555] = 4.0E0*I_ERI_Fx2z_Gxy2z_S_S_bb;
  abcd[1556] = 4.0E0*I_ERI_F3y_Gxy2z_S_S_bb;
  abcd[1557] = 4.0E0*I_ERI_F2yz_Gxy2z_S_S_bb;
  abcd[1558] = 4.0E0*I_ERI_Fy2z_Gxy2z_S_S_bb;
  abcd[1559] = 4.0E0*I_ERI_F3z_Gxy2z_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_G_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1560] = 4.0E0*I_ERI_F3x_G3xz_S_S_bb-2.0E0*2*I_ERI_F3x_Dxz_S_S_b;
  abcd[1561] = 4.0E0*I_ERI_F2xy_G3xz_S_S_bb-2.0E0*2*I_ERI_F2xy_Dxz_S_S_b;
  abcd[1562] = 4.0E0*I_ERI_F2xz_G3xz_S_S_bb-2.0E0*2*I_ERI_F2xz_Dxz_S_S_b;
  abcd[1563] = 4.0E0*I_ERI_Fx2y_G3xz_S_S_bb-2.0E0*2*I_ERI_Fx2y_Dxz_S_S_b;
  abcd[1564] = 4.0E0*I_ERI_Fxyz_G3xz_S_S_bb-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_b;
  abcd[1565] = 4.0E0*I_ERI_Fx2z_G3xz_S_S_bb-2.0E0*2*I_ERI_Fx2z_Dxz_S_S_b;
  abcd[1566] = 4.0E0*I_ERI_F3y_G3xz_S_S_bb-2.0E0*2*I_ERI_F3y_Dxz_S_S_b;
  abcd[1567] = 4.0E0*I_ERI_F2yz_G3xz_S_S_bb-2.0E0*2*I_ERI_F2yz_Dxz_S_S_b;
  abcd[1568] = 4.0E0*I_ERI_Fy2z_G3xz_S_S_bb-2.0E0*2*I_ERI_Fy2z_Dxz_S_S_b;
  abcd[1569] = 4.0E0*I_ERI_F3z_G3xz_S_S_bb-2.0E0*2*I_ERI_F3z_Dxz_S_S_b;
  abcd[1570] = 4.0E0*I_ERI_F3x_G2xyz_S_S_bb-2.0E0*1*I_ERI_F3x_Dyz_S_S_b;
  abcd[1571] = 4.0E0*I_ERI_F2xy_G2xyz_S_S_bb-2.0E0*1*I_ERI_F2xy_Dyz_S_S_b;
  abcd[1572] = 4.0E0*I_ERI_F2xz_G2xyz_S_S_bb-2.0E0*1*I_ERI_F2xz_Dyz_S_S_b;
  abcd[1573] = 4.0E0*I_ERI_Fx2y_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_b;
  abcd[1574] = 4.0E0*I_ERI_Fxyz_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_b;
  abcd[1575] = 4.0E0*I_ERI_Fx2z_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_b;
  abcd[1576] = 4.0E0*I_ERI_F3y_G2xyz_S_S_bb-2.0E0*1*I_ERI_F3y_Dyz_S_S_b;
  abcd[1577] = 4.0E0*I_ERI_F2yz_G2xyz_S_S_bb-2.0E0*1*I_ERI_F2yz_Dyz_S_S_b;
  abcd[1578] = 4.0E0*I_ERI_Fy2z_G2xyz_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_b;
  abcd[1579] = 4.0E0*I_ERI_F3z_G2xyz_S_S_bb-2.0E0*1*I_ERI_F3z_Dyz_S_S_b;
  abcd[1580] = 4.0E0*I_ERI_F3x_G2x2z_S_S_bb-2.0E0*1*I_ERI_F3x_D2x_S_S_b-2.0E0*1*I_ERI_F3x_D2z_S_S_b+1*I_ERI_F3x_S_S_S;
  abcd[1581] = 4.0E0*I_ERI_F2xy_G2x2z_S_S_bb-2.0E0*1*I_ERI_F2xy_D2x_S_S_b-2.0E0*1*I_ERI_F2xy_D2z_S_S_b+1*I_ERI_F2xy_S_S_S;
  abcd[1582] = 4.0E0*I_ERI_F2xz_G2x2z_S_S_bb-2.0E0*1*I_ERI_F2xz_D2x_S_S_b-2.0E0*1*I_ERI_F2xz_D2z_S_S_b+1*I_ERI_F2xz_S_S_S;
  abcd[1583] = 4.0E0*I_ERI_Fx2y_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fx2y_D2x_S_S_b-2.0E0*1*I_ERI_Fx2y_D2z_S_S_b+1*I_ERI_Fx2y_S_S_S;
  abcd[1584] = 4.0E0*I_ERI_Fxyz_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fxyz_D2x_S_S_b-2.0E0*1*I_ERI_Fxyz_D2z_S_S_b+1*I_ERI_Fxyz_S_S_S;
  abcd[1585] = 4.0E0*I_ERI_Fx2z_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fx2z_D2x_S_S_b-2.0E0*1*I_ERI_Fx2z_D2z_S_S_b+1*I_ERI_Fx2z_S_S_S;
  abcd[1586] = 4.0E0*I_ERI_F3y_G2x2z_S_S_bb-2.0E0*1*I_ERI_F3y_D2x_S_S_b-2.0E0*1*I_ERI_F3y_D2z_S_S_b+1*I_ERI_F3y_S_S_S;
  abcd[1587] = 4.0E0*I_ERI_F2yz_G2x2z_S_S_bb-2.0E0*1*I_ERI_F2yz_D2x_S_S_b-2.0E0*1*I_ERI_F2yz_D2z_S_S_b+1*I_ERI_F2yz_S_S_S;
  abcd[1588] = 4.0E0*I_ERI_Fy2z_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fy2z_D2x_S_S_b-2.0E0*1*I_ERI_Fy2z_D2z_S_S_b+1*I_ERI_Fy2z_S_S_S;
  abcd[1589] = 4.0E0*I_ERI_F3z_G2x2z_S_S_bb-2.0E0*1*I_ERI_F3z_D2x_S_S_b-2.0E0*1*I_ERI_F3z_D2z_S_S_b+1*I_ERI_F3z_S_S_S;
  abcd[1590] = 4.0E0*I_ERI_F3x_Gx2yz_S_S_bb;
  abcd[1591] = 4.0E0*I_ERI_F2xy_Gx2yz_S_S_bb;
  abcd[1592] = 4.0E0*I_ERI_F2xz_Gx2yz_S_S_bb;
  abcd[1593] = 4.0E0*I_ERI_Fx2y_Gx2yz_S_S_bb;
  abcd[1594] = 4.0E0*I_ERI_Fxyz_Gx2yz_S_S_bb;
  abcd[1595] = 4.0E0*I_ERI_Fx2z_Gx2yz_S_S_bb;
  abcd[1596] = 4.0E0*I_ERI_F3y_Gx2yz_S_S_bb;
  abcd[1597] = 4.0E0*I_ERI_F2yz_Gx2yz_S_S_bb;
  abcd[1598] = 4.0E0*I_ERI_Fy2z_Gx2yz_S_S_bb;
  abcd[1599] = 4.0E0*I_ERI_F3z_Gx2yz_S_S_bb;
  abcd[1600] = 4.0E0*I_ERI_F3x_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F3x_Dxy_S_S_b;
  abcd[1601] = 4.0E0*I_ERI_F2xy_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxy_S_S_b;
  abcd[1602] = 4.0E0*I_ERI_F2xz_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxy_S_S_b;
  abcd[1603] = 4.0E0*I_ERI_Fx2y_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_b;
  abcd[1604] = 4.0E0*I_ERI_Fxyz_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_b;
  abcd[1605] = 4.0E0*I_ERI_Fx2z_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_b;
  abcd[1606] = 4.0E0*I_ERI_F3y_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F3y_Dxy_S_S_b;
  abcd[1607] = 4.0E0*I_ERI_F2yz_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxy_S_S_b;
  abcd[1608] = 4.0E0*I_ERI_Fy2z_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_b;
  abcd[1609] = 4.0E0*I_ERI_F3z_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F3z_Dxy_S_S_b;
  abcd[1610] = 4.0E0*I_ERI_F3x_Gx3z_S_S_bb-2.0E0*2*I_ERI_F3x_Dxz_S_S_b;
  abcd[1611] = 4.0E0*I_ERI_F2xy_Gx3z_S_S_bb-2.0E0*2*I_ERI_F2xy_Dxz_S_S_b;
  abcd[1612] = 4.0E0*I_ERI_F2xz_Gx3z_S_S_bb-2.0E0*2*I_ERI_F2xz_Dxz_S_S_b;
  abcd[1613] = 4.0E0*I_ERI_Fx2y_Gx3z_S_S_bb-2.0E0*2*I_ERI_Fx2y_Dxz_S_S_b;
  abcd[1614] = 4.0E0*I_ERI_Fxyz_Gx3z_S_S_bb-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_b;
  abcd[1615] = 4.0E0*I_ERI_Fx2z_Gx3z_S_S_bb-2.0E0*2*I_ERI_Fx2z_Dxz_S_S_b;
  abcd[1616] = 4.0E0*I_ERI_F3y_Gx3z_S_S_bb-2.0E0*2*I_ERI_F3y_Dxz_S_S_b;
  abcd[1617] = 4.0E0*I_ERI_F2yz_Gx3z_S_S_bb-2.0E0*2*I_ERI_F2yz_Dxz_S_S_b;
  abcd[1618] = 4.0E0*I_ERI_Fy2z_Gx3z_S_S_bb-2.0E0*2*I_ERI_Fy2z_Dxz_S_S_b;
  abcd[1619] = 4.0E0*I_ERI_F3z_Gx3z_S_S_bb-2.0E0*2*I_ERI_F3z_Dxz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_G_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1620] = 4.0E0*I_ERI_F3x_G2x2y_S_S_bb-2.0E0*1*I_ERI_F3x_D2x_S_S_b;
  abcd[1621] = 4.0E0*I_ERI_F2xy_G2x2y_S_S_bb-2.0E0*1*I_ERI_F2xy_D2x_S_S_b;
  abcd[1622] = 4.0E0*I_ERI_F2xz_G2x2y_S_S_bb-2.0E0*1*I_ERI_F2xz_D2x_S_S_b;
  abcd[1623] = 4.0E0*I_ERI_Fx2y_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fx2y_D2x_S_S_b;
  abcd[1624] = 4.0E0*I_ERI_Fxyz_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fxyz_D2x_S_S_b;
  abcd[1625] = 4.0E0*I_ERI_Fx2z_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fx2z_D2x_S_S_b;
  abcd[1626] = 4.0E0*I_ERI_F3y_G2x2y_S_S_bb-2.0E0*1*I_ERI_F3y_D2x_S_S_b;
  abcd[1627] = 4.0E0*I_ERI_F2yz_G2x2y_S_S_bb-2.0E0*1*I_ERI_F2yz_D2x_S_S_b;
  abcd[1628] = 4.0E0*I_ERI_Fy2z_G2x2y_S_S_bb-2.0E0*1*I_ERI_Fy2z_D2x_S_S_b;
  abcd[1629] = 4.0E0*I_ERI_F3z_G2x2y_S_S_bb-2.0E0*1*I_ERI_F3z_D2x_S_S_b;
  abcd[1630] = 4.0E0*I_ERI_F3x_Gx3y_S_S_bb-2.0E0*1*I_ERI_F3x_Dxy_S_S_b-2.0E0*2*I_ERI_F3x_Dxy_S_S_b;
  abcd[1631] = 4.0E0*I_ERI_F2xy_Gx3y_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxy_S_S_b-2.0E0*2*I_ERI_F2xy_Dxy_S_S_b;
  abcd[1632] = 4.0E0*I_ERI_F2xz_Gx3y_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxy_S_S_b-2.0E0*2*I_ERI_F2xz_Dxy_S_S_b;
  abcd[1633] = 4.0E0*I_ERI_Fx2y_Gx3y_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_b-2.0E0*2*I_ERI_Fx2y_Dxy_S_S_b;
  abcd[1634] = 4.0E0*I_ERI_Fxyz_Gx3y_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_b-2.0E0*2*I_ERI_Fxyz_Dxy_S_S_b;
  abcd[1635] = 4.0E0*I_ERI_Fx2z_Gx3y_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_b-2.0E0*2*I_ERI_Fx2z_Dxy_S_S_b;
  abcd[1636] = 4.0E0*I_ERI_F3y_Gx3y_S_S_bb-2.0E0*1*I_ERI_F3y_Dxy_S_S_b-2.0E0*2*I_ERI_F3y_Dxy_S_S_b;
  abcd[1637] = 4.0E0*I_ERI_F2yz_Gx3y_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxy_S_S_b-2.0E0*2*I_ERI_F2yz_Dxy_S_S_b;
  abcd[1638] = 4.0E0*I_ERI_Fy2z_Gx3y_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_b-2.0E0*2*I_ERI_Fy2z_Dxy_S_S_b;
  abcd[1639] = 4.0E0*I_ERI_F3z_Gx3y_S_S_bb-2.0E0*1*I_ERI_F3z_Dxy_S_S_b-2.0E0*2*I_ERI_F3z_Dxy_S_S_b;
  abcd[1640] = 4.0E0*I_ERI_F3x_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F3x_Dxz_S_S_b;
  abcd[1641] = 4.0E0*I_ERI_F2xy_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxz_S_S_b;
  abcd[1642] = 4.0E0*I_ERI_F2xz_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxz_S_S_b;
  abcd[1643] = 4.0E0*I_ERI_Fx2y_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_b;
  abcd[1644] = 4.0E0*I_ERI_Fxyz_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_b;
  abcd[1645] = 4.0E0*I_ERI_Fx2z_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_b;
  abcd[1646] = 4.0E0*I_ERI_F3y_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F3y_Dxz_S_S_b;
  abcd[1647] = 4.0E0*I_ERI_F2yz_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxz_S_S_b;
  abcd[1648] = 4.0E0*I_ERI_Fy2z_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_b;
  abcd[1649] = 4.0E0*I_ERI_F3z_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F3z_Dxz_S_S_b;
  abcd[1650] = 4.0E0*I_ERI_F3x_G4y_S_S_bb-2.0E0*2*I_ERI_F3x_D2y_S_S_b-2.0E0*3*I_ERI_F3x_D2y_S_S_b+2*1*I_ERI_F3x_S_S_S;
  abcd[1651] = 4.0E0*I_ERI_F2xy_G4y_S_S_bb-2.0E0*2*I_ERI_F2xy_D2y_S_S_b-2.0E0*3*I_ERI_F2xy_D2y_S_S_b+2*1*I_ERI_F2xy_S_S_S;
  abcd[1652] = 4.0E0*I_ERI_F2xz_G4y_S_S_bb-2.0E0*2*I_ERI_F2xz_D2y_S_S_b-2.0E0*3*I_ERI_F2xz_D2y_S_S_b+2*1*I_ERI_F2xz_S_S_S;
  abcd[1653] = 4.0E0*I_ERI_Fx2y_G4y_S_S_bb-2.0E0*2*I_ERI_Fx2y_D2y_S_S_b-2.0E0*3*I_ERI_Fx2y_D2y_S_S_b+2*1*I_ERI_Fx2y_S_S_S;
  abcd[1654] = 4.0E0*I_ERI_Fxyz_G4y_S_S_bb-2.0E0*2*I_ERI_Fxyz_D2y_S_S_b-2.0E0*3*I_ERI_Fxyz_D2y_S_S_b+2*1*I_ERI_Fxyz_S_S_S;
  abcd[1655] = 4.0E0*I_ERI_Fx2z_G4y_S_S_bb-2.0E0*2*I_ERI_Fx2z_D2y_S_S_b-2.0E0*3*I_ERI_Fx2z_D2y_S_S_b+2*1*I_ERI_Fx2z_S_S_S;
  abcd[1656] = 4.0E0*I_ERI_F3y_G4y_S_S_bb-2.0E0*2*I_ERI_F3y_D2y_S_S_b-2.0E0*3*I_ERI_F3y_D2y_S_S_b+2*1*I_ERI_F3y_S_S_S;
  abcd[1657] = 4.0E0*I_ERI_F2yz_G4y_S_S_bb-2.0E0*2*I_ERI_F2yz_D2y_S_S_b-2.0E0*3*I_ERI_F2yz_D2y_S_S_b+2*1*I_ERI_F2yz_S_S_S;
  abcd[1658] = 4.0E0*I_ERI_Fy2z_G4y_S_S_bb-2.0E0*2*I_ERI_Fy2z_D2y_S_S_b-2.0E0*3*I_ERI_Fy2z_D2y_S_S_b+2*1*I_ERI_Fy2z_S_S_S;
  abcd[1659] = 4.0E0*I_ERI_F3z_G4y_S_S_bb-2.0E0*2*I_ERI_F3z_D2y_S_S_b-2.0E0*3*I_ERI_F3z_D2y_S_S_b+2*1*I_ERI_F3z_S_S_S;
  abcd[1660] = 4.0E0*I_ERI_F3x_G3yz_S_S_bb-2.0E0*1*I_ERI_F3x_Dyz_S_S_b-2.0E0*2*I_ERI_F3x_Dyz_S_S_b;
  abcd[1661] = 4.0E0*I_ERI_F2xy_G3yz_S_S_bb-2.0E0*1*I_ERI_F2xy_Dyz_S_S_b-2.0E0*2*I_ERI_F2xy_Dyz_S_S_b;
  abcd[1662] = 4.0E0*I_ERI_F2xz_G3yz_S_S_bb-2.0E0*1*I_ERI_F2xz_Dyz_S_S_b-2.0E0*2*I_ERI_F2xz_Dyz_S_S_b;
  abcd[1663] = 4.0E0*I_ERI_Fx2y_G3yz_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_b-2.0E0*2*I_ERI_Fx2y_Dyz_S_S_b;
  abcd[1664] = 4.0E0*I_ERI_Fxyz_G3yz_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_b-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_b;
  abcd[1665] = 4.0E0*I_ERI_Fx2z_G3yz_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_b-2.0E0*2*I_ERI_Fx2z_Dyz_S_S_b;
  abcd[1666] = 4.0E0*I_ERI_F3y_G3yz_S_S_bb-2.0E0*1*I_ERI_F3y_Dyz_S_S_b-2.0E0*2*I_ERI_F3y_Dyz_S_S_b;
  abcd[1667] = 4.0E0*I_ERI_F2yz_G3yz_S_S_bb-2.0E0*1*I_ERI_F2yz_Dyz_S_S_b-2.0E0*2*I_ERI_F2yz_Dyz_S_S_b;
  abcd[1668] = 4.0E0*I_ERI_Fy2z_G3yz_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_b-2.0E0*2*I_ERI_Fy2z_Dyz_S_S_b;
  abcd[1669] = 4.0E0*I_ERI_F3z_G3yz_S_S_bb-2.0E0*1*I_ERI_F3z_Dyz_S_S_b-2.0E0*2*I_ERI_F3z_Dyz_S_S_b;
  abcd[1670] = 4.0E0*I_ERI_F3x_G2y2z_S_S_bb-2.0E0*1*I_ERI_F3x_D2z_S_S_b;
  abcd[1671] = 4.0E0*I_ERI_F2xy_G2y2z_S_S_bb-2.0E0*1*I_ERI_F2xy_D2z_S_S_b;
  abcd[1672] = 4.0E0*I_ERI_F2xz_G2y2z_S_S_bb-2.0E0*1*I_ERI_F2xz_D2z_S_S_b;
  abcd[1673] = 4.0E0*I_ERI_Fx2y_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fx2y_D2z_S_S_b;
  abcd[1674] = 4.0E0*I_ERI_Fxyz_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fxyz_D2z_S_S_b;
  abcd[1675] = 4.0E0*I_ERI_Fx2z_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fx2z_D2z_S_S_b;
  abcd[1676] = 4.0E0*I_ERI_F3y_G2y2z_S_S_bb-2.0E0*1*I_ERI_F3y_D2z_S_S_b;
  abcd[1677] = 4.0E0*I_ERI_F2yz_G2y2z_S_S_bb-2.0E0*1*I_ERI_F2yz_D2z_S_S_b;
  abcd[1678] = 4.0E0*I_ERI_Fy2z_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fy2z_D2z_S_S_b;
  abcd[1679] = 4.0E0*I_ERI_F3z_G2y2z_S_S_bb-2.0E0*1*I_ERI_F3z_D2z_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_G_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1680] = 4.0E0*I_ERI_F3x_G2xyz_S_S_bb;
  abcd[1681] = 4.0E0*I_ERI_F2xy_G2xyz_S_S_bb;
  abcd[1682] = 4.0E0*I_ERI_F2xz_G2xyz_S_S_bb;
  abcd[1683] = 4.0E0*I_ERI_Fx2y_G2xyz_S_S_bb;
  abcd[1684] = 4.0E0*I_ERI_Fxyz_G2xyz_S_S_bb;
  abcd[1685] = 4.0E0*I_ERI_Fx2z_G2xyz_S_S_bb;
  abcd[1686] = 4.0E0*I_ERI_F3y_G2xyz_S_S_bb;
  abcd[1687] = 4.0E0*I_ERI_F2yz_G2xyz_S_S_bb;
  abcd[1688] = 4.0E0*I_ERI_Fy2z_G2xyz_S_S_bb;
  abcd[1689] = 4.0E0*I_ERI_F3z_G2xyz_S_S_bb;
  abcd[1690] = 4.0E0*I_ERI_F3x_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F3x_Dxz_S_S_b;
  abcd[1691] = 4.0E0*I_ERI_F2xy_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxz_S_S_b;
  abcd[1692] = 4.0E0*I_ERI_F2xz_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxz_S_S_b;
  abcd[1693] = 4.0E0*I_ERI_Fx2y_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_b;
  abcd[1694] = 4.0E0*I_ERI_Fxyz_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_b;
  abcd[1695] = 4.0E0*I_ERI_Fx2z_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_b;
  abcd[1696] = 4.0E0*I_ERI_F3y_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F3y_Dxz_S_S_b;
  abcd[1697] = 4.0E0*I_ERI_F2yz_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxz_S_S_b;
  abcd[1698] = 4.0E0*I_ERI_Fy2z_Gx2yz_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_b;
  abcd[1699] = 4.0E0*I_ERI_F3z_Gx2yz_S_S_bb-2.0E0*1*I_ERI_F3z_Dxz_S_S_b;
  abcd[1700] = 4.0E0*I_ERI_F3x_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F3x_Dxy_S_S_b;
  abcd[1701] = 4.0E0*I_ERI_F2xy_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxy_S_S_b;
  abcd[1702] = 4.0E0*I_ERI_F2xz_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxy_S_S_b;
  abcd[1703] = 4.0E0*I_ERI_Fx2y_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_b;
  abcd[1704] = 4.0E0*I_ERI_Fxyz_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_b;
  abcd[1705] = 4.0E0*I_ERI_Fx2z_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_b;
  abcd[1706] = 4.0E0*I_ERI_F3y_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F3y_Dxy_S_S_b;
  abcd[1707] = 4.0E0*I_ERI_F2yz_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxy_S_S_b;
  abcd[1708] = 4.0E0*I_ERI_Fy2z_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_b;
  abcd[1709] = 4.0E0*I_ERI_F3z_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F3z_Dxy_S_S_b;
  abcd[1710] = 4.0E0*I_ERI_F3x_G3yz_S_S_bb-2.0E0*2*I_ERI_F3x_Dyz_S_S_b;
  abcd[1711] = 4.0E0*I_ERI_F2xy_G3yz_S_S_bb-2.0E0*2*I_ERI_F2xy_Dyz_S_S_b;
  abcd[1712] = 4.0E0*I_ERI_F2xz_G3yz_S_S_bb-2.0E0*2*I_ERI_F2xz_Dyz_S_S_b;
  abcd[1713] = 4.0E0*I_ERI_Fx2y_G3yz_S_S_bb-2.0E0*2*I_ERI_Fx2y_Dyz_S_S_b;
  abcd[1714] = 4.0E0*I_ERI_Fxyz_G3yz_S_S_bb-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_b;
  abcd[1715] = 4.0E0*I_ERI_Fx2z_G3yz_S_S_bb-2.0E0*2*I_ERI_Fx2z_Dyz_S_S_b;
  abcd[1716] = 4.0E0*I_ERI_F3y_G3yz_S_S_bb-2.0E0*2*I_ERI_F3y_Dyz_S_S_b;
  abcd[1717] = 4.0E0*I_ERI_F2yz_G3yz_S_S_bb-2.0E0*2*I_ERI_F2yz_Dyz_S_S_b;
  abcd[1718] = 4.0E0*I_ERI_Fy2z_G3yz_S_S_bb-2.0E0*2*I_ERI_Fy2z_Dyz_S_S_b;
  abcd[1719] = 4.0E0*I_ERI_F3z_G3yz_S_S_bb-2.0E0*2*I_ERI_F3z_Dyz_S_S_b;
  abcd[1720] = 4.0E0*I_ERI_F3x_G2y2z_S_S_bb-2.0E0*1*I_ERI_F3x_D2y_S_S_b-2.0E0*1*I_ERI_F3x_D2z_S_S_b+1*I_ERI_F3x_S_S_S;
  abcd[1721] = 4.0E0*I_ERI_F2xy_G2y2z_S_S_bb-2.0E0*1*I_ERI_F2xy_D2y_S_S_b-2.0E0*1*I_ERI_F2xy_D2z_S_S_b+1*I_ERI_F2xy_S_S_S;
  abcd[1722] = 4.0E0*I_ERI_F2xz_G2y2z_S_S_bb-2.0E0*1*I_ERI_F2xz_D2y_S_S_b-2.0E0*1*I_ERI_F2xz_D2z_S_S_b+1*I_ERI_F2xz_S_S_S;
  abcd[1723] = 4.0E0*I_ERI_Fx2y_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fx2y_D2y_S_S_b-2.0E0*1*I_ERI_Fx2y_D2z_S_S_b+1*I_ERI_Fx2y_S_S_S;
  abcd[1724] = 4.0E0*I_ERI_Fxyz_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fxyz_D2y_S_S_b-2.0E0*1*I_ERI_Fxyz_D2z_S_S_b+1*I_ERI_Fxyz_S_S_S;
  abcd[1725] = 4.0E0*I_ERI_Fx2z_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fx2z_D2y_S_S_b-2.0E0*1*I_ERI_Fx2z_D2z_S_S_b+1*I_ERI_Fx2z_S_S_S;
  abcd[1726] = 4.0E0*I_ERI_F3y_G2y2z_S_S_bb-2.0E0*1*I_ERI_F3y_D2y_S_S_b-2.0E0*1*I_ERI_F3y_D2z_S_S_b+1*I_ERI_F3y_S_S_S;
  abcd[1727] = 4.0E0*I_ERI_F2yz_G2y2z_S_S_bb-2.0E0*1*I_ERI_F2yz_D2y_S_S_b-2.0E0*1*I_ERI_F2yz_D2z_S_S_b+1*I_ERI_F2yz_S_S_S;
  abcd[1728] = 4.0E0*I_ERI_Fy2z_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fy2z_D2y_S_S_b-2.0E0*1*I_ERI_Fy2z_D2z_S_S_b+1*I_ERI_Fy2z_S_S_S;
  abcd[1729] = 4.0E0*I_ERI_F3z_G2y2z_S_S_bb-2.0E0*1*I_ERI_F3z_D2y_S_S_b-2.0E0*1*I_ERI_F3z_D2z_S_S_b+1*I_ERI_F3z_S_S_S;
  abcd[1730] = 4.0E0*I_ERI_F3x_Gy3z_S_S_bb-2.0E0*2*I_ERI_F3x_Dyz_S_S_b;
  abcd[1731] = 4.0E0*I_ERI_F2xy_Gy3z_S_S_bb-2.0E0*2*I_ERI_F2xy_Dyz_S_S_b;
  abcd[1732] = 4.0E0*I_ERI_F2xz_Gy3z_S_S_bb-2.0E0*2*I_ERI_F2xz_Dyz_S_S_b;
  abcd[1733] = 4.0E0*I_ERI_Fx2y_Gy3z_S_S_bb-2.0E0*2*I_ERI_Fx2y_Dyz_S_S_b;
  abcd[1734] = 4.0E0*I_ERI_Fxyz_Gy3z_S_S_bb-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_b;
  abcd[1735] = 4.0E0*I_ERI_Fx2z_Gy3z_S_S_bb-2.0E0*2*I_ERI_Fx2z_Dyz_S_S_b;
  abcd[1736] = 4.0E0*I_ERI_F3y_Gy3z_S_S_bb-2.0E0*2*I_ERI_F3y_Dyz_S_S_b;
  abcd[1737] = 4.0E0*I_ERI_F2yz_Gy3z_S_S_bb-2.0E0*2*I_ERI_F2yz_Dyz_S_S_b;
  abcd[1738] = 4.0E0*I_ERI_Fy2z_Gy3z_S_S_bb-2.0E0*2*I_ERI_Fy2z_Dyz_S_S_b;
  abcd[1739] = 4.0E0*I_ERI_F3z_Gy3z_S_S_bb-2.0E0*2*I_ERI_F3z_Dyz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_G_S_S_bb
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  abcd[1740] = 4.0E0*I_ERI_F3x_G2x2z_S_S_bb-2.0E0*1*I_ERI_F3x_D2x_S_S_b;
  abcd[1741] = 4.0E0*I_ERI_F2xy_G2x2z_S_S_bb-2.0E0*1*I_ERI_F2xy_D2x_S_S_b;
  abcd[1742] = 4.0E0*I_ERI_F2xz_G2x2z_S_S_bb-2.0E0*1*I_ERI_F2xz_D2x_S_S_b;
  abcd[1743] = 4.0E0*I_ERI_Fx2y_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fx2y_D2x_S_S_b;
  abcd[1744] = 4.0E0*I_ERI_Fxyz_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fxyz_D2x_S_S_b;
  abcd[1745] = 4.0E0*I_ERI_Fx2z_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fx2z_D2x_S_S_b;
  abcd[1746] = 4.0E0*I_ERI_F3y_G2x2z_S_S_bb-2.0E0*1*I_ERI_F3y_D2x_S_S_b;
  abcd[1747] = 4.0E0*I_ERI_F2yz_G2x2z_S_S_bb-2.0E0*1*I_ERI_F2yz_D2x_S_S_b;
  abcd[1748] = 4.0E0*I_ERI_Fy2z_G2x2z_S_S_bb-2.0E0*1*I_ERI_Fy2z_D2x_S_S_b;
  abcd[1749] = 4.0E0*I_ERI_F3z_G2x2z_S_S_bb-2.0E0*1*I_ERI_F3z_D2x_S_S_b;
  abcd[1750] = 4.0E0*I_ERI_F3x_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F3x_Dxy_S_S_b;
  abcd[1751] = 4.0E0*I_ERI_F2xy_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxy_S_S_b;
  abcd[1752] = 4.0E0*I_ERI_F2xz_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxy_S_S_b;
  abcd[1753] = 4.0E0*I_ERI_Fx2y_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_b;
  abcd[1754] = 4.0E0*I_ERI_Fxyz_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_b;
  abcd[1755] = 4.0E0*I_ERI_Fx2z_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_b;
  abcd[1756] = 4.0E0*I_ERI_F3y_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F3y_Dxy_S_S_b;
  abcd[1757] = 4.0E0*I_ERI_F2yz_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxy_S_S_b;
  abcd[1758] = 4.0E0*I_ERI_Fy2z_Gxy2z_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_b;
  abcd[1759] = 4.0E0*I_ERI_F3z_Gxy2z_S_S_bb-2.0E0*1*I_ERI_F3z_Dxy_S_S_b;
  abcd[1760] = 4.0E0*I_ERI_F3x_Gx3z_S_S_bb-2.0E0*1*I_ERI_F3x_Dxz_S_S_b-2.0E0*2*I_ERI_F3x_Dxz_S_S_b;
  abcd[1761] = 4.0E0*I_ERI_F2xy_Gx3z_S_S_bb-2.0E0*1*I_ERI_F2xy_Dxz_S_S_b-2.0E0*2*I_ERI_F2xy_Dxz_S_S_b;
  abcd[1762] = 4.0E0*I_ERI_F2xz_Gx3z_S_S_bb-2.0E0*1*I_ERI_F2xz_Dxz_S_S_b-2.0E0*2*I_ERI_F2xz_Dxz_S_S_b;
  abcd[1763] = 4.0E0*I_ERI_Fx2y_Gx3z_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_b-2.0E0*2*I_ERI_Fx2y_Dxz_S_S_b;
  abcd[1764] = 4.0E0*I_ERI_Fxyz_Gx3z_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_b-2.0E0*2*I_ERI_Fxyz_Dxz_S_S_b;
  abcd[1765] = 4.0E0*I_ERI_Fx2z_Gx3z_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_b-2.0E0*2*I_ERI_Fx2z_Dxz_S_S_b;
  abcd[1766] = 4.0E0*I_ERI_F3y_Gx3z_S_S_bb-2.0E0*1*I_ERI_F3y_Dxz_S_S_b-2.0E0*2*I_ERI_F3y_Dxz_S_S_b;
  abcd[1767] = 4.0E0*I_ERI_F2yz_Gx3z_S_S_bb-2.0E0*1*I_ERI_F2yz_Dxz_S_S_b-2.0E0*2*I_ERI_F2yz_Dxz_S_S_b;
  abcd[1768] = 4.0E0*I_ERI_Fy2z_Gx3z_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_b-2.0E0*2*I_ERI_Fy2z_Dxz_S_S_b;
  abcd[1769] = 4.0E0*I_ERI_F3z_Gx3z_S_S_bb-2.0E0*1*I_ERI_F3z_Dxz_S_S_b-2.0E0*2*I_ERI_F3z_Dxz_S_S_b;
  abcd[1770] = 4.0E0*I_ERI_F3x_G2y2z_S_S_bb-2.0E0*1*I_ERI_F3x_D2y_S_S_b;
  abcd[1771] = 4.0E0*I_ERI_F2xy_G2y2z_S_S_bb-2.0E0*1*I_ERI_F2xy_D2y_S_S_b;
  abcd[1772] = 4.0E0*I_ERI_F2xz_G2y2z_S_S_bb-2.0E0*1*I_ERI_F2xz_D2y_S_S_b;
  abcd[1773] = 4.0E0*I_ERI_Fx2y_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fx2y_D2y_S_S_b;
  abcd[1774] = 4.0E0*I_ERI_Fxyz_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fxyz_D2y_S_S_b;
  abcd[1775] = 4.0E0*I_ERI_Fx2z_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fx2z_D2y_S_S_b;
  abcd[1776] = 4.0E0*I_ERI_F3y_G2y2z_S_S_bb-2.0E0*1*I_ERI_F3y_D2y_S_S_b;
  abcd[1777] = 4.0E0*I_ERI_F2yz_G2y2z_S_S_bb-2.0E0*1*I_ERI_F2yz_D2y_S_S_b;
  abcd[1778] = 4.0E0*I_ERI_Fy2z_G2y2z_S_S_bb-2.0E0*1*I_ERI_Fy2z_D2y_S_S_b;
  abcd[1779] = 4.0E0*I_ERI_F3z_G2y2z_S_S_bb-2.0E0*1*I_ERI_F3z_D2y_S_S_b;
  abcd[1780] = 4.0E0*I_ERI_F3x_Gy3z_S_S_bb-2.0E0*1*I_ERI_F3x_Dyz_S_S_b-2.0E0*2*I_ERI_F3x_Dyz_S_S_b;
  abcd[1781] = 4.0E0*I_ERI_F2xy_Gy3z_S_S_bb-2.0E0*1*I_ERI_F2xy_Dyz_S_S_b-2.0E0*2*I_ERI_F2xy_Dyz_S_S_b;
  abcd[1782] = 4.0E0*I_ERI_F2xz_Gy3z_S_S_bb-2.0E0*1*I_ERI_F2xz_Dyz_S_S_b-2.0E0*2*I_ERI_F2xz_Dyz_S_S_b;
  abcd[1783] = 4.0E0*I_ERI_Fx2y_Gy3z_S_S_bb-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_b-2.0E0*2*I_ERI_Fx2y_Dyz_S_S_b;
  abcd[1784] = 4.0E0*I_ERI_Fxyz_Gy3z_S_S_bb-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_b-2.0E0*2*I_ERI_Fxyz_Dyz_S_S_b;
  abcd[1785] = 4.0E0*I_ERI_Fx2z_Gy3z_S_S_bb-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_b-2.0E0*2*I_ERI_Fx2z_Dyz_S_S_b;
  abcd[1786] = 4.0E0*I_ERI_F3y_Gy3z_S_S_bb-2.0E0*1*I_ERI_F3y_Dyz_S_S_b-2.0E0*2*I_ERI_F3y_Dyz_S_S_b;
  abcd[1787] = 4.0E0*I_ERI_F2yz_Gy3z_S_S_bb-2.0E0*1*I_ERI_F2yz_Dyz_S_S_b-2.0E0*2*I_ERI_F2yz_Dyz_S_S_b;
  abcd[1788] = 4.0E0*I_ERI_Fy2z_Gy3z_S_S_bb-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_b-2.0E0*2*I_ERI_Fy2z_Dyz_S_S_b;
  abcd[1789] = 4.0E0*I_ERI_F3z_Gy3z_S_S_bb-2.0E0*1*I_ERI_F3z_Dyz_S_S_b-2.0E0*2*I_ERI_F3z_Dyz_S_S_b;
  abcd[1790] = 4.0E0*I_ERI_F3x_G4z_S_S_bb-2.0E0*2*I_ERI_F3x_D2z_S_S_b-2.0E0*3*I_ERI_F3x_D2z_S_S_b+2*1*I_ERI_F3x_S_S_S;
  abcd[1791] = 4.0E0*I_ERI_F2xy_G4z_S_S_bb-2.0E0*2*I_ERI_F2xy_D2z_S_S_b-2.0E0*3*I_ERI_F2xy_D2z_S_S_b+2*1*I_ERI_F2xy_S_S_S;
  abcd[1792] = 4.0E0*I_ERI_F2xz_G4z_S_S_bb-2.0E0*2*I_ERI_F2xz_D2z_S_S_b-2.0E0*3*I_ERI_F2xz_D2z_S_S_b+2*1*I_ERI_F2xz_S_S_S;
  abcd[1793] = 4.0E0*I_ERI_Fx2y_G4z_S_S_bb-2.0E0*2*I_ERI_Fx2y_D2z_S_S_b-2.0E0*3*I_ERI_Fx2y_D2z_S_S_b+2*1*I_ERI_Fx2y_S_S_S;
  abcd[1794] = 4.0E0*I_ERI_Fxyz_G4z_S_S_bb-2.0E0*2*I_ERI_Fxyz_D2z_S_S_b-2.0E0*3*I_ERI_Fxyz_D2z_S_S_b+2*1*I_ERI_Fxyz_S_S_S;
  abcd[1795] = 4.0E0*I_ERI_Fx2z_G4z_S_S_bb-2.0E0*2*I_ERI_Fx2z_D2z_S_S_b-2.0E0*3*I_ERI_Fx2z_D2z_S_S_b+2*1*I_ERI_Fx2z_S_S_S;
  abcd[1796] = 4.0E0*I_ERI_F3y_G4z_S_S_bb-2.0E0*2*I_ERI_F3y_D2z_S_S_b-2.0E0*3*I_ERI_F3y_D2z_S_S_b+2*1*I_ERI_F3y_S_S_S;
  abcd[1797] = 4.0E0*I_ERI_F2yz_G4z_S_S_bb-2.0E0*2*I_ERI_F2yz_D2z_S_S_b-2.0E0*3*I_ERI_F2yz_D2z_S_S_b+2*1*I_ERI_F2yz_S_S_S;
  abcd[1798] = 4.0E0*I_ERI_Fy2z_G4z_S_S_bb-2.0E0*2*I_ERI_Fy2z_D2z_S_S_b-2.0E0*3*I_ERI_Fy2z_D2z_S_S_b+2*1*I_ERI_Fy2z_S_S_S;
  abcd[1799] = 4.0E0*I_ERI_F3z_G4z_S_S_bb-2.0E0*2*I_ERI_F3z_D2z_S_S_b-2.0E0*3*I_ERI_F3z_D2z_S_S_b+2*1*I_ERI_F3z_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  abcd[1800] = 4.0E0*I_ERI_F3x_F3x_Px_S_bc-2.0E0*2*I_ERI_F3x_Px_Px_S_c;
  abcd[1801] = 4.0E0*I_ERI_F2xy_F3x_Px_S_bc-2.0E0*2*I_ERI_F2xy_Px_Px_S_c;
  abcd[1802] = 4.0E0*I_ERI_F2xz_F3x_Px_S_bc-2.0E0*2*I_ERI_F2xz_Px_Px_S_c;
  abcd[1803] = 4.0E0*I_ERI_Fx2y_F3x_Px_S_bc-2.0E0*2*I_ERI_Fx2y_Px_Px_S_c;
  abcd[1804] = 4.0E0*I_ERI_Fxyz_F3x_Px_S_bc-2.0E0*2*I_ERI_Fxyz_Px_Px_S_c;
  abcd[1805] = 4.0E0*I_ERI_Fx2z_F3x_Px_S_bc-2.0E0*2*I_ERI_Fx2z_Px_Px_S_c;
  abcd[1806] = 4.0E0*I_ERI_F3y_F3x_Px_S_bc-2.0E0*2*I_ERI_F3y_Px_Px_S_c;
  abcd[1807] = 4.0E0*I_ERI_F2yz_F3x_Px_S_bc-2.0E0*2*I_ERI_F2yz_Px_Px_S_c;
  abcd[1808] = 4.0E0*I_ERI_Fy2z_F3x_Px_S_bc-2.0E0*2*I_ERI_Fy2z_Px_Px_S_c;
  abcd[1809] = 4.0E0*I_ERI_F3z_F3x_Px_S_bc-2.0E0*2*I_ERI_F3z_Px_Px_S_c;
  abcd[1810] = 4.0E0*I_ERI_F3x_F2xy_Px_S_bc-2.0E0*1*I_ERI_F3x_Py_Px_S_c;
  abcd[1811] = 4.0E0*I_ERI_F2xy_F2xy_Px_S_bc-2.0E0*1*I_ERI_F2xy_Py_Px_S_c;
  abcd[1812] = 4.0E0*I_ERI_F2xz_F2xy_Px_S_bc-2.0E0*1*I_ERI_F2xz_Py_Px_S_c;
  abcd[1813] = 4.0E0*I_ERI_Fx2y_F2xy_Px_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Px_S_c;
  abcd[1814] = 4.0E0*I_ERI_Fxyz_F2xy_Px_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Px_S_c;
  abcd[1815] = 4.0E0*I_ERI_Fx2z_F2xy_Px_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Px_S_c;
  abcd[1816] = 4.0E0*I_ERI_F3y_F2xy_Px_S_bc-2.0E0*1*I_ERI_F3y_Py_Px_S_c;
  abcd[1817] = 4.0E0*I_ERI_F2yz_F2xy_Px_S_bc-2.0E0*1*I_ERI_F2yz_Py_Px_S_c;
  abcd[1818] = 4.0E0*I_ERI_Fy2z_F2xy_Px_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Px_S_c;
  abcd[1819] = 4.0E0*I_ERI_F3z_F2xy_Px_S_bc-2.0E0*1*I_ERI_F3z_Py_Px_S_c;
  abcd[1820] = 4.0E0*I_ERI_F3x_F2xz_Px_S_bc-2.0E0*1*I_ERI_F3x_Pz_Px_S_c;
  abcd[1821] = 4.0E0*I_ERI_F2xy_F2xz_Px_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Px_S_c;
  abcd[1822] = 4.0E0*I_ERI_F2xz_F2xz_Px_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Px_S_c;
  abcd[1823] = 4.0E0*I_ERI_Fx2y_F2xz_Px_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Px_S_c;
  abcd[1824] = 4.0E0*I_ERI_Fxyz_F2xz_Px_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Px_S_c;
  abcd[1825] = 4.0E0*I_ERI_Fx2z_F2xz_Px_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Px_S_c;
  abcd[1826] = 4.0E0*I_ERI_F3y_F2xz_Px_S_bc-2.0E0*1*I_ERI_F3y_Pz_Px_S_c;
  abcd[1827] = 4.0E0*I_ERI_F2yz_F2xz_Px_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Px_S_c;
  abcd[1828] = 4.0E0*I_ERI_Fy2z_F2xz_Px_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Px_S_c;
  abcd[1829] = 4.0E0*I_ERI_F3z_F2xz_Px_S_bc-2.0E0*1*I_ERI_F3z_Pz_Px_S_c;
  abcd[1830] = 4.0E0*I_ERI_F3x_Fx2y_Px_S_bc;
  abcd[1831] = 4.0E0*I_ERI_F2xy_Fx2y_Px_S_bc;
  abcd[1832] = 4.0E0*I_ERI_F2xz_Fx2y_Px_S_bc;
  abcd[1833] = 4.0E0*I_ERI_Fx2y_Fx2y_Px_S_bc;
  abcd[1834] = 4.0E0*I_ERI_Fxyz_Fx2y_Px_S_bc;
  abcd[1835] = 4.0E0*I_ERI_Fx2z_Fx2y_Px_S_bc;
  abcd[1836] = 4.0E0*I_ERI_F3y_Fx2y_Px_S_bc;
  abcd[1837] = 4.0E0*I_ERI_F2yz_Fx2y_Px_S_bc;
  abcd[1838] = 4.0E0*I_ERI_Fy2z_Fx2y_Px_S_bc;
  abcd[1839] = 4.0E0*I_ERI_F3z_Fx2y_Px_S_bc;
  abcd[1840] = 4.0E0*I_ERI_F3x_Fxyz_Px_S_bc;
  abcd[1841] = 4.0E0*I_ERI_F2xy_Fxyz_Px_S_bc;
  abcd[1842] = 4.0E0*I_ERI_F2xz_Fxyz_Px_S_bc;
  abcd[1843] = 4.0E0*I_ERI_Fx2y_Fxyz_Px_S_bc;
  abcd[1844] = 4.0E0*I_ERI_Fxyz_Fxyz_Px_S_bc;
  abcd[1845] = 4.0E0*I_ERI_Fx2z_Fxyz_Px_S_bc;
  abcd[1846] = 4.0E0*I_ERI_F3y_Fxyz_Px_S_bc;
  abcd[1847] = 4.0E0*I_ERI_F2yz_Fxyz_Px_S_bc;
  abcd[1848] = 4.0E0*I_ERI_Fy2z_Fxyz_Px_S_bc;
  abcd[1849] = 4.0E0*I_ERI_F3z_Fxyz_Px_S_bc;
  abcd[1850] = 4.0E0*I_ERI_F3x_Fx2z_Px_S_bc;
  abcd[1851] = 4.0E0*I_ERI_F2xy_Fx2z_Px_S_bc;
  abcd[1852] = 4.0E0*I_ERI_F2xz_Fx2z_Px_S_bc;
  abcd[1853] = 4.0E0*I_ERI_Fx2y_Fx2z_Px_S_bc;
  abcd[1854] = 4.0E0*I_ERI_Fxyz_Fx2z_Px_S_bc;
  abcd[1855] = 4.0E0*I_ERI_Fx2z_Fx2z_Px_S_bc;
  abcd[1856] = 4.0E0*I_ERI_F3y_Fx2z_Px_S_bc;
  abcd[1857] = 4.0E0*I_ERI_F2yz_Fx2z_Px_S_bc;
  abcd[1858] = 4.0E0*I_ERI_Fy2z_Fx2z_Px_S_bc;
  abcd[1859] = 4.0E0*I_ERI_F3z_Fx2z_Px_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  abcd[1860] = 4.0E0*I_ERI_F3x_F3x_Py_S_bc-2.0E0*2*I_ERI_F3x_Px_Py_S_c;
  abcd[1861] = 4.0E0*I_ERI_F2xy_F3x_Py_S_bc-2.0E0*2*I_ERI_F2xy_Px_Py_S_c;
  abcd[1862] = 4.0E0*I_ERI_F2xz_F3x_Py_S_bc-2.0E0*2*I_ERI_F2xz_Px_Py_S_c;
  abcd[1863] = 4.0E0*I_ERI_Fx2y_F3x_Py_S_bc-2.0E0*2*I_ERI_Fx2y_Px_Py_S_c;
  abcd[1864] = 4.0E0*I_ERI_Fxyz_F3x_Py_S_bc-2.0E0*2*I_ERI_Fxyz_Px_Py_S_c;
  abcd[1865] = 4.0E0*I_ERI_Fx2z_F3x_Py_S_bc-2.0E0*2*I_ERI_Fx2z_Px_Py_S_c;
  abcd[1866] = 4.0E0*I_ERI_F3y_F3x_Py_S_bc-2.0E0*2*I_ERI_F3y_Px_Py_S_c;
  abcd[1867] = 4.0E0*I_ERI_F2yz_F3x_Py_S_bc-2.0E0*2*I_ERI_F2yz_Px_Py_S_c;
  abcd[1868] = 4.0E0*I_ERI_Fy2z_F3x_Py_S_bc-2.0E0*2*I_ERI_Fy2z_Px_Py_S_c;
  abcd[1869] = 4.0E0*I_ERI_F3z_F3x_Py_S_bc-2.0E0*2*I_ERI_F3z_Px_Py_S_c;
  abcd[1870] = 4.0E0*I_ERI_F3x_F2xy_Py_S_bc-2.0E0*1*I_ERI_F3x_Py_Py_S_c;
  abcd[1871] = 4.0E0*I_ERI_F2xy_F2xy_Py_S_bc-2.0E0*1*I_ERI_F2xy_Py_Py_S_c;
  abcd[1872] = 4.0E0*I_ERI_F2xz_F2xy_Py_S_bc-2.0E0*1*I_ERI_F2xz_Py_Py_S_c;
  abcd[1873] = 4.0E0*I_ERI_Fx2y_F2xy_Py_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Py_S_c;
  abcd[1874] = 4.0E0*I_ERI_Fxyz_F2xy_Py_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Py_S_c;
  abcd[1875] = 4.0E0*I_ERI_Fx2z_F2xy_Py_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Py_S_c;
  abcd[1876] = 4.0E0*I_ERI_F3y_F2xy_Py_S_bc-2.0E0*1*I_ERI_F3y_Py_Py_S_c;
  abcd[1877] = 4.0E0*I_ERI_F2yz_F2xy_Py_S_bc-2.0E0*1*I_ERI_F2yz_Py_Py_S_c;
  abcd[1878] = 4.0E0*I_ERI_Fy2z_F2xy_Py_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Py_S_c;
  abcd[1879] = 4.0E0*I_ERI_F3z_F2xy_Py_S_bc-2.0E0*1*I_ERI_F3z_Py_Py_S_c;
  abcd[1880] = 4.0E0*I_ERI_F3x_F2xz_Py_S_bc-2.0E0*1*I_ERI_F3x_Pz_Py_S_c;
  abcd[1881] = 4.0E0*I_ERI_F2xy_F2xz_Py_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Py_S_c;
  abcd[1882] = 4.0E0*I_ERI_F2xz_F2xz_Py_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Py_S_c;
  abcd[1883] = 4.0E0*I_ERI_Fx2y_F2xz_Py_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Py_S_c;
  abcd[1884] = 4.0E0*I_ERI_Fxyz_F2xz_Py_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Py_S_c;
  abcd[1885] = 4.0E0*I_ERI_Fx2z_F2xz_Py_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Py_S_c;
  abcd[1886] = 4.0E0*I_ERI_F3y_F2xz_Py_S_bc-2.0E0*1*I_ERI_F3y_Pz_Py_S_c;
  abcd[1887] = 4.0E0*I_ERI_F2yz_F2xz_Py_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Py_S_c;
  abcd[1888] = 4.0E0*I_ERI_Fy2z_F2xz_Py_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Py_S_c;
  abcd[1889] = 4.0E0*I_ERI_F3z_F2xz_Py_S_bc-2.0E0*1*I_ERI_F3z_Pz_Py_S_c;
  abcd[1890] = 4.0E0*I_ERI_F3x_Fx2y_Py_S_bc;
  abcd[1891] = 4.0E0*I_ERI_F2xy_Fx2y_Py_S_bc;
  abcd[1892] = 4.0E0*I_ERI_F2xz_Fx2y_Py_S_bc;
  abcd[1893] = 4.0E0*I_ERI_Fx2y_Fx2y_Py_S_bc;
  abcd[1894] = 4.0E0*I_ERI_Fxyz_Fx2y_Py_S_bc;
  abcd[1895] = 4.0E0*I_ERI_Fx2z_Fx2y_Py_S_bc;
  abcd[1896] = 4.0E0*I_ERI_F3y_Fx2y_Py_S_bc;
  abcd[1897] = 4.0E0*I_ERI_F2yz_Fx2y_Py_S_bc;
  abcd[1898] = 4.0E0*I_ERI_Fy2z_Fx2y_Py_S_bc;
  abcd[1899] = 4.0E0*I_ERI_F3z_Fx2y_Py_S_bc;
  abcd[1900] = 4.0E0*I_ERI_F3x_Fxyz_Py_S_bc;
  abcd[1901] = 4.0E0*I_ERI_F2xy_Fxyz_Py_S_bc;
  abcd[1902] = 4.0E0*I_ERI_F2xz_Fxyz_Py_S_bc;
  abcd[1903] = 4.0E0*I_ERI_Fx2y_Fxyz_Py_S_bc;
  abcd[1904] = 4.0E0*I_ERI_Fxyz_Fxyz_Py_S_bc;
  abcd[1905] = 4.0E0*I_ERI_Fx2z_Fxyz_Py_S_bc;
  abcd[1906] = 4.0E0*I_ERI_F3y_Fxyz_Py_S_bc;
  abcd[1907] = 4.0E0*I_ERI_F2yz_Fxyz_Py_S_bc;
  abcd[1908] = 4.0E0*I_ERI_Fy2z_Fxyz_Py_S_bc;
  abcd[1909] = 4.0E0*I_ERI_F3z_Fxyz_Py_S_bc;
  abcd[1910] = 4.0E0*I_ERI_F3x_Fx2z_Py_S_bc;
  abcd[1911] = 4.0E0*I_ERI_F2xy_Fx2z_Py_S_bc;
  abcd[1912] = 4.0E0*I_ERI_F2xz_Fx2z_Py_S_bc;
  abcd[1913] = 4.0E0*I_ERI_Fx2y_Fx2z_Py_S_bc;
  abcd[1914] = 4.0E0*I_ERI_Fxyz_Fx2z_Py_S_bc;
  abcd[1915] = 4.0E0*I_ERI_Fx2z_Fx2z_Py_S_bc;
  abcd[1916] = 4.0E0*I_ERI_F3y_Fx2z_Py_S_bc;
  abcd[1917] = 4.0E0*I_ERI_F2yz_Fx2z_Py_S_bc;
  abcd[1918] = 4.0E0*I_ERI_Fy2z_Fx2z_Py_S_bc;
  abcd[1919] = 4.0E0*I_ERI_F3z_Fx2z_Py_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  abcd[1920] = 4.0E0*I_ERI_F3x_F3x_Pz_S_bc-2.0E0*2*I_ERI_F3x_Px_Pz_S_c;
  abcd[1921] = 4.0E0*I_ERI_F2xy_F3x_Pz_S_bc-2.0E0*2*I_ERI_F2xy_Px_Pz_S_c;
  abcd[1922] = 4.0E0*I_ERI_F2xz_F3x_Pz_S_bc-2.0E0*2*I_ERI_F2xz_Px_Pz_S_c;
  abcd[1923] = 4.0E0*I_ERI_Fx2y_F3x_Pz_S_bc-2.0E0*2*I_ERI_Fx2y_Px_Pz_S_c;
  abcd[1924] = 4.0E0*I_ERI_Fxyz_F3x_Pz_S_bc-2.0E0*2*I_ERI_Fxyz_Px_Pz_S_c;
  abcd[1925] = 4.0E0*I_ERI_Fx2z_F3x_Pz_S_bc-2.0E0*2*I_ERI_Fx2z_Px_Pz_S_c;
  abcd[1926] = 4.0E0*I_ERI_F3y_F3x_Pz_S_bc-2.0E0*2*I_ERI_F3y_Px_Pz_S_c;
  abcd[1927] = 4.0E0*I_ERI_F2yz_F3x_Pz_S_bc-2.0E0*2*I_ERI_F2yz_Px_Pz_S_c;
  abcd[1928] = 4.0E0*I_ERI_Fy2z_F3x_Pz_S_bc-2.0E0*2*I_ERI_Fy2z_Px_Pz_S_c;
  abcd[1929] = 4.0E0*I_ERI_F3z_F3x_Pz_S_bc-2.0E0*2*I_ERI_F3z_Px_Pz_S_c;
  abcd[1930] = 4.0E0*I_ERI_F3x_F2xy_Pz_S_bc-2.0E0*1*I_ERI_F3x_Py_Pz_S_c;
  abcd[1931] = 4.0E0*I_ERI_F2xy_F2xy_Pz_S_bc-2.0E0*1*I_ERI_F2xy_Py_Pz_S_c;
  abcd[1932] = 4.0E0*I_ERI_F2xz_F2xy_Pz_S_bc-2.0E0*1*I_ERI_F2xz_Py_Pz_S_c;
  abcd[1933] = 4.0E0*I_ERI_Fx2y_F2xy_Pz_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Pz_S_c;
  abcd[1934] = 4.0E0*I_ERI_Fxyz_F2xy_Pz_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Pz_S_c;
  abcd[1935] = 4.0E0*I_ERI_Fx2z_F2xy_Pz_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Pz_S_c;
  abcd[1936] = 4.0E0*I_ERI_F3y_F2xy_Pz_S_bc-2.0E0*1*I_ERI_F3y_Py_Pz_S_c;
  abcd[1937] = 4.0E0*I_ERI_F2yz_F2xy_Pz_S_bc-2.0E0*1*I_ERI_F2yz_Py_Pz_S_c;
  abcd[1938] = 4.0E0*I_ERI_Fy2z_F2xy_Pz_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Pz_S_c;
  abcd[1939] = 4.0E0*I_ERI_F3z_F2xy_Pz_S_bc-2.0E0*1*I_ERI_F3z_Py_Pz_S_c;
  abcd[1940] = 4.0E0*I_ERI_F3x_F2xz_Pz_S_bc-2.0E0*1*I_ERI_F3x_Pz_Pz_S_c;
  abcd[1941] = 4.0E0*I_ERI_F2xy_F2xz_Pz_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Pz_S_c;
  abcd[1942] = 4.0E0*I_ERI_F2xz_F2xz_Pz_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Pz_S_c;
  abcd[1943] = 4.0E0*I_ERI_Fx2y_F2xz_Pz_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Pz_S_c;
  abcd[1944] = 4.0E0*I_ERI_Fxyz_F2xz_Pz_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Pz_S_c;
  abcd[1945] = 4.0E0*I_ERI_Fx2z_F2xz_Pz_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Pz_S_c;
  abcd[1946] = 4.0E0*I_ERI_F3y_F2xz_Pz_S_bc-2.0E0*1*I_ERI_F3y_Pz_Pz_S_c;
  abcd[1947] = 4.0E0*I_ERI_F2yz_F2xz_Pz_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Pz_S_c;
  abcd[1948] = 4.0E0*I_ERI_Fy2z_F2xz_Pz_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Pz_S_c;
  abcd[1949] = 4.0E0*I_ERI_F3z_F2xz_Pz_S_bc-2.0E0*1*I_ERI_F3z_Pz_Pz_S_c;
  abcd[1950] = 4.0E0*I_ERI_F3x_Fx2y_Pz_S_bc;
  abcd[1951] = 4.0E0*I_ERI_F2xy_Fx2y_Pz_S_bc;
  abcd[1952] = 4.0E0*I_ERI_F2xz_Fx2y_Pz_S_bc;
  abcd[1953] = 4.0E0*I_ERI_Fx2y_Fx2y_Pz_S_bc;
  abcd[1954] = 4.0E0*I_ERI_Fxyz_Fx2y_Pz_S_bc;
  abcd[1955] = 4.0E0*I_ERI_Fx2z_Fx2y_Pz_S_bc;
  abcd[1956] = 4.0E0*I_ERI_F3y_Fx2y_Pz_S_bc;
  abcd[1957] = 4.0E0*I_ERI_F2yz_Fx2y_Pz_S_bc;
  abcd[1958] = 4.0E0*I_ERI_Fy2z_Fx2y_Pz_S_bc;
  abcd[1959] = 4.0E0*I_ERI_F3z_Fx2y_Pz_S_bc;
  abcd[1960] = 4.0E0*I_ERI_F3x_Fxyz_Pz_S_bc;
  abcd[1961] = 4.0E0*I_ERI_F2xy_Fxyz_Pz_S_bc;
  abcd[1962] = 4.0E0*I_ERI_F2xz_Fxyz_Pz_S_bc;
  abcd[1963] = 4.0E0*I_ERI_Fx2y_Fxyz_Pz_S_bc;
  abcd[1964] = 4.0E0*I_ERI_Fxyz_Fxyz_Pz_S_bc;
  abcd[1965] = 4.0E0*I_ERI_Fx2z_Fxyz_Pz_S_bc;
  abcd[1966] = 4.0E0*I_ERI_F3y_Fxyz_Pz_S_bc;
  abcd[1967] = 4.0E0*I_ERI_F2yz_Fxyz_Pz_S_bc;
  abcd[1968] = 4.0E0*I_ERI_Fy2z_Fxyz_Pz_S_bc;
  abcd[1969] = 4.0E0*I_ERI_F3z_Fxyz_Pz_S_bc;
  abcd[1970] = 4.0E0*I_ERI_F3x_Fx2z_Pz_S_bc;
  abcd[1971] = 4.0E0*I_ERI_F2xy_Fx2z_Pz_S_bc;
  abcd[1972] = 4.0E0*I_ERI_F2xz_Fx2z_Pz_S_bc;
  abcd[1973] = 4.0E0*I_ERI_Fx2y_Fx2z_Pz_S_bc;
  abcd[1974] = 4.0E0*I_ERI_Fxyz_Fx2z_Pz_S_bc;
  abcd[1975] = 4.0E0*I_ERI_Fx2z_Fx2z_Pz_S_bc;
  abcd[1976] = 4.0E0*I_ERI_F3y_Fx2z_Pz_S_bc;
  abcd[1977] = 4.0E0*I_ERI_F2yz_Fx2z_Pz_S_bc;
  abcd[1978] = 4.0E0*I_ERI_Fy2z_Fx2z_Pz_S_bc;
  abcd[1979] = 4.0E0*I_ERI_F3z_Fx2z_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  abcd[1980] = 4.0E0*I_ERI_F3x_F2xy_Px_S_bc;
  abcd[1981] = 4.0E0*I_ERI_F2xy_F2xy_Px_S_bc;
  abcd[1982] = 4.0E0*I_ERI_F2xz_F2xy_Px_S_bc;
  abcd[1983] = 4.0E0*I_ERI_Fx2y_F2xy_Px_S_bc;
  abcd[1984] = 4.0E0*I_ERI_Fxyz_F2xy_Px_S_bc;
  abcd[1985] = 4.0E0*I_ERI_Fx2z_F2xy_Px_S_bc;
  abcd[1986] = 4.0E0*I_ERI_F3y_F2xy_Px_S_bc;
  abcd[1987] = 4.0E0*I_ERI_F2yz_F2xy_Px_S_bc;
  abcd[1988] = 4.0E0*I_ERI_Fy2z_F2xy_Px_S_bc;
  abcd[1989] = 4.0E0*I_ERI_F3z_F2xy_Px_S_bc;
  abcd[1990] = 4.0E0*I_ERI_F3x_Fx2y_Px_S_bc-2.0E0*1*I_ERI_F3x_Px_Px_S_c;
  abcd[1991] = 4.0E0*I_ERI_F2xy_Fx2y_Px_S_bc-2.0E0*1*I_ERI_F2xy_Px_Px_S_c;
  abcd[1992] = 4.0E0*I_ERI_F2xz_Fx2y_Px_S_bc-2.0E0*1*I_ERI_F2xz_Px_Px_S_c;
  abcd[1993] = 4.0E0*I_ERI_Fx2y_Fx2y_Px_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Px_S_c;
  abcd[1994] = 4.0E0*I_ERI_Fxyz_Fx2y_Px_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Px_S_c;
  abcd[1995] = 4.0E0*I_ERI_Fx2z_Fx2y_Px_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Px_S_c;
  abcd[1996] = 4.0E0*I_ERI_F3y_Fx2y_Px_S_bc-2.0E0*1*I_ERI_F3y_Px_Px_S_c;
  abcd[1997] = 4.0E0*I_ERI_F2yz_Fx2y_Px_S_bc-2.0E0*1*I_ERI_F2yz_Px_Px_S_c;
  abcd[1998] = 4.0E0*I_ERI_Fy2z_Fx2y_Px_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Px_S_c;
  abcd[1999] = 4.0E0*I_ERI_F3z_Fx2y_Px_S_bc-2.0E0*1*I_ERI_F3z_Px_Px_S_c;
  abcd[2000] = 4.0E0*I_ERI_F3x_Fxyz_Px_S_bc;
  abcd[2001] = 4.0E0*I_ERI_F2xy_Fxyz_Px_S_bc;
  abcd[2002] = 4.0E0*I_ERI_F2xz_Fxyz_Px_S_bc;
  abcd[2003] = 4.0E0*I_ERI_Fx2y_Fxyz_Px_S_bc;
  abcd[2004] = 4.0E0*I_ERI_Fxyz_Fxyz_Px_S_bc;
  abcd[2005] = 4.0E0*I_ERI_Fx2z_Fxyz_Px_S_bc;
  abcd[2006] = 4.0E0*I_ERI_F3y_Fxyz_Px_S_bc;
  abcd[2007] = 4.0E0*I_ERI_F2yz_Fxyz_Px_S_bc;
  abcd[2008] = 4.0E0*I_ERI_Fy2z_Fxyz_Px_S_bc;
  abcd[2009] = 4.0E0*I_ERI_F3z_Fxyz_Px_S_bc;
  abcd[2010] = 4.0E0*I_ERI_F3x_F3y_Px_S_bc-2.0E0*2*I_ERI_F3x_Py_Px_S_c;
  abcd[2011] = 4.0E0*I_ERI_F2xy_F3y_Px_S_bc-2.0E0*2*I_ERI_F2xy_Py_Px_S_c;
  abcd[2012] = 4.0E0*I_ERI_F2xz_F3y_Px_S_bc-2.0E0*2*I_ERI_F2xz_Py_Px_S_c;
  abcd[2013] = 4.0E0*I_ERI_Fx2y_F3y_Px_S_bc-2.0E0*2*I_ERI_Fx2y_Py_Px_S_c;
  abcd[2014] = 4.0E0*I_ERI_Fxyz_F3y_Px_S_bc-2.0E0*2*I_ERI_Fxyz_Py_Px_S_c;
  abcd[2015] = 4.0E0*I_ERI_Fx2z_F3y_Px_S_bc-2.0E0*2*I_ERI_Fx2z_Py_Px_S_c;
  abcd[2016] = 4.0E0*I_ERI_F3y_F3y_Px_S_bc-2.0E0*2*I_ERI_F3y_Py_Px_S_c;
  abcd[2017] = 4.0E0*I_ERI_F2yz_F3y_Px_S_bc-2.0E0*2*I_ERI_F2yz_Py_Px_S_c;
  abcd[2018] = 4.0E0*I_ERI_Fy2z_F3y_Px_S_bc-2.0E0*2*I_ERI_Fy2z_Py_Px_S_c;
  abcd[2019] = 4.0E0*I_ERI_F3z_F3y_Px_S_bc-2.0E0*2*I_ERI_F3z_Py_Px_S_c;
  abcd[2020] = 4.0E0*I_ERI_F3x_F2yz_Px_S_bc-2.0E0*1*I_ERI_F3x_Pz_Px_S_c;
  abcd[2021] = 4.0E0*I_ERI_F2xy_F2yz_Px_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Px_S_c;
  abcd[2022] = 4.0E0*I_ERI_F2xz_F2yz_Px_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Px_S_c;
  abcd[2023] = 4.0E0*I_ERI_Fx2y_F2yz_Px_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Px_S_c;
  abcd[2024] = 4.0E0*I_ERI_Fxyz_F2yz_Px_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Px_S_c;
  abcd[2025] = 4.0E0*I_ERI_Fx2z_F2yz_Px_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Px_S_c;
  abcd[2026] = 4.0E0*I_ERI_F3y_F2yz_Px_S_bc-2.0E0*1*I_ERI_F3y_Pz_Px_S_c;
  abcd[2027] = 4.0E0*I_ERI_F2yz_F2yz_Px_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Px_S_c;
  abcd[2028] = 4.0E0*I_ERI_Fy2z_F2yz_Px_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Px_S_c;
  abcd[2029] = 4.0E0*I_ERI_F3z_F2yz_Px_S_bc-2.0E0*1*I_ERI_F3z_Pz_Px_S_c;
  abcd[2030] = 4.0E0*I_ERI_F3x_Fy2z_Px_S_bc;
  abcd[2031] = 4.0E0*I_ERI_F2xy_Fy2z_Px_S_bc;
  abcd[2032] = 4.0E0*I_ERI_F2xz_Fy2z_Px_S_bc;
  abcd[2033] = 4.0E0*I_ERI_Fx2y_Fy2z_Px_S_bc;
  abcd[2034] = 4.0E0*I_ERI_Fxyz_Fy2z_Px_S_bc;
  abcd[2035] = 4.0E0*I_ERI_Fx2z_Fy2z_Px_S_bc;
  abcd[2036] = 4.0E0*I_ERI_F3y_Fy2z_Px_S_bc;
  abcd[2037] = 4.0E0*I_ERI_F2yz_Fy2z_Px_S_bc;
  abcd[2038] = 4.0E0*I_ERI_Fy2z_Fy2z_Px_S_bc;
  abcd[2039] = 4.0E0*I_ERI_F3z_Fy2z_Px_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  abcd[2040] = 4.0E0*I_ERI_F3x_F2xy_Py_S_bc;
  abcd[2041] = 4.0E0*I_ERI_F2xy_F2xy_Py_S_bc;
  abcd[2042] = 4.0E0*I_ERI_F2xz_F2xy_Py_S_bc;
  abcd[2043] = 4.0E0*I_ERI_Fx2y_F2xy_Py_S_bc;
  abcd[2044] = 4.0E0*I_ERI_Fxyz_F2xy_Py_S_bc;
  abcd[2045] = 4.0E0*I_ERI_Fx2z_F2xy_Py_S_bc;
  abcd[2046] = 4.0E0*I_ERI_F3y_F2xy_Py_S_bc;
  abcd[2047] = 4.0E0*I_ERI_F2yz_F2xy_Py_S_bc;
  abcd[2048] = 4.0E0*I_ERI_Fy2z_F2xy_Py_S_bc;
  abcd[2049] = 4.0E0*I_ERI_F3z_F2xy_Py_S_bc;
  abcd[2050] = 4.0E0*I_ERI_F3x_Fx2y_Py_S_bc-2.0E0*1*I_ERI_F3x_Px_Py_S_c;
  abcd[2051] = 4.0E0*I_ERI_F2xy_Fx2y_Py_S_bc-2.0E0*1*I_ERI_F2xy_Px_Py_S_c;
  abcd[2052] = 4.0E0*I_ERI_F2xz_Fx2y_Py_S_bc-2.0E0*1*I_ERI_F2xz_Px_Py_S_c;
  abcd[2053] = 4.0E0*I_ERI_Fx2y_Fx2y_Py_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Py_S_c;
  abcd[2054] = 4.0E0*I_ERI_Fxyz_Fx2y_Py_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Py_S_c;
  abcd[2055] = 4.0E0*I_ERI_Fx2z_Fx2y_Py_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Py_S_c;
  abcd[2056] = 4.0E0*I_ERI_F3y_Fx2y_Py_S_bc-2.0E0*1*I_ERI_F3y_Px_Py_S_c;
  abcd[2057] = 4.0E0*I_ERI_F2yz_Fx2y_Py_S_bc-2.0E0*1*I_ERI_F2yz_Px_Py_S_c;
  abcd[2058] = 4.0E0*I_ERI_Fy2z_Fx2y_Py_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Py_S_c;
  abcd[2059] = 4.0E0*I_ERI_F3z_Fx2y_Py_S_bc-2.0E0*1*I_ERI_F3z_Px_Py_S_c;
  abcd[2060] = 4.0E0*I_ERI_F3x_Fxyz_Py_S_bc;
  abcd[2061] = 4.0E0*I_ERI_F2xy_Fxyz_Py_S_bc;
  abcd[2062] = 4.0E0*I_ERI_F2xz_Fxyz_Py_S_bc;
  abcd[2063] = 4.0E0*I_ERI_Fx2y_Fxyz_Py_S_bc;
  abcd[2064] = 4.0E0*I_ERI_Fxyz_Fxyz_Py_S_bc;
  abcd[2065] = 4.0E0*I_ERI_Fx2z_Fxyz_Py_S_bc;
  abcd[2066] = 4.0E0*I_ERI_F3y_Fxyz_Py_S_bc;
  abcd[2067] = 4.0E0*I_ERI_F2yz_Fxyz_Py_S_bc;
  abcd[2068] = 4.0E0*I_ERI_Fy2z_Fxyz_Py_S_bc;
  abcd[2069] = 4.0E0*I_ERI_F3z_Fxyz_Py_S_bc;
  abcd[2070] = 4.0E0*I_ERI_F3x_F3y_Py_S_bc-2.0E0*2*I_ERI_F3x_Py_Py_S_c;
  abcd[2071] = 4.0E0*I_ERI_F2xy_F3y_Py_S_bc-2.0E0*2*I_ERI_F2xy_Py_Py_S_c;
  abcd[2072] = 4.0E0*I_ERI_F2xz_F3y_Py_S_bc-2.0E0*2*I_ERI_F2xz_Py_Py_S_c;
  abcd[2073] = 4.0E0*I_ERI_Fx2y_F3y_Py_S_bc-2.0E0*2*I_ERI_Fx2y_Py_Py_S_c;
  abcd[2074] = 4.0E0*I_ERI_Fxyz_F3y_Py_S_bc-2.0E0*2*I_ERI_Fxyz_Py_Py_S_c;
  abcd[2075] = 4.0E0*I_ERI_Fx2z_F3y_Py_S_bc-2.0E0*2*I_ERI_Fx2z_Py_Py_S_c;
  abcd[2076] = 4.0E0*I_ERI_F3y_F3y_Py_S_bc-2.0E0*2*I_ERI_F3y_Py_Py_S_c;
  abcd[2077] = 4.0E0*I_ERI_F2yz_F3y_Py_S_bc-2.0E0*2*I_ERI_F2yz_Py_Py_S_c;
  abcd[2078] = 4.0E0*I_ERI_Fy2z_F3y_Py_S_bc-2.0E0*2*I_ERI_Fy2z_Py_Py_S_c;
  abcd[2079] = 4.0E0*I_ERI_F3z_F3y_Py_S_bc-2.0E0*2*I_ERI_F3z_Py_Py_S_c;
  abcd[2080] = 4.0E0*I_ERI_F3x_F2yz_Py_S_bc-2.0E0*1*I_ERI_F3x_Pz_Py_S_c;
  abcd[2081] = 4.0E0*I_ERI_F2xy_F2yz_Py_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Py_S_c;
  abcd[2082] = 4.0E0*I_ERI_F2xz_F2yz_Py_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Py_S_c;
  abcd[2083] = 4.0E0*I_ERI_Fx2y_F2yz_Py_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Py_S_c;
  abcd[2084] = 4.0E0*I_ERI_Fxyz_F2yz_Py_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Py_S_c;
  abcd[2085] = 4.0E0*I_ERI_Fx2z_F2yz_Py_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Py_S_c;
  abcd[2086] = 4.0E0*I_ERI_F3y_F2yz_Py_S_bc-2.0E0*1*I_ERI_F3y_Pz_Py_S_c;
  abcd[2087] = 4.0E0*I_ERI_F2yz_F2yz_Py_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Py_S_c;
  abcd[2088] = 4.0E0*I_ERI_Fy2z_F2yz_Py_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Py_S_c;
  abcd[2089] = 4.0E0*I_ERI_F3z_F2yz_Py_S_bc-2.0E0*1*I_ERI_F3z_Pz_Py_S_c;
  abcd[2090] = 4.0E0*I_ERI_F3x_Fy2z_Py_S_bc;
  abcd[2091] = 4.0E0*I_ERI_F2xy_Fy2z_Py_S_bc;
  abcd[2092] = 4.0E0*I_ERI_F2xz_Fy2z_Py_S_bc;
  abcd[2093] = 4.0E0*I_ERI_Fx2y_Fy2z_Py_S_bc;
  abcd[2094] = 4.0E0*I_ERI_Fxyz_Fy2z_Py_S_bc;
  abcd[2095] = 4.0E0*I_ERI_Fx2z_Fy2z_Py_S_bc;
  abcd[2096] = 4.0E0*I_ERI_F3y_Fy2z_Py_S_bc;
  abcd[2097] = 4.0E0*I_ERI_F2yz_Fy2z_Py_S_bc;
  abcd[2098] = 4.0E0*I_ERI_Fy2z_Fy2z_Py_S_bc;
  abcd[2099] = 4.0E0*I_ERI_F3z_Fy2z_Py_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  abcd[2100] = 4.0E0*I_ERI_F3x_F2xy_Pz_S_bc;
  abcd[2101] = 4.0E0*I_ERI_F2xy_F2xy_Pz_S_bc;
  abcd[2102] = 4.0E0*I_ERI_F2xz_F2xy_Pz_S_bc;
  abcd[2103] = 4.0E0*I_ERI_Fx2y_F2xy_Pz_S_bc;
  abcd[2104] = 4.0E0*I_ERI_Fxyz_F2xy_Pz_S_bc;
  abcd[2105] = 4.0E0*I_ERI_Fx2z_F2xy_Pz_S_bc;
  abcd[2106] = 4.0E0*I_ERI_F3y_F2xy_Pz_S_bc;
  abcd[2107] = 4.0E0*I_ERI_F2yz_F2xy_Pz_S_bc;
  abcd[2108] = 4.0E0*I_ERI_Fy2z_F2xy_Pz_S_bc;
  abcd[2109] = 4.0E0*I_ERI_F3z_F2xy_Pz_S_bc;
  abcd[2110] = 4.0E0*I_ERI_F3x_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_F3x_Px_Pz_S_c;
  abcd[2111] = 4.0E0*I_ERI_F2xy_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_F2xy_Px_Pz_S_c;
  abcd[2112] = 4.0E0*I_ERI_F2xz_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_F2xz_Px_Pz_S_c;
  abcd[2113] = 4.0E0*I_ERI_Fx2y_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Pz_S_c;
  abcd[2114] = 4.0E0*I_ERI_Fxyz_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Pz_S_c;
  abcd[2115] = 4.0E0*I_ERI_Fx2z_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Pz_S_c;
  abcd[2116] = 4.0E0*I_ERI_F3y_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_F3y_Px_Pz_S_c;
  abcd[2117] = 4.0E0*I_ERI_F2yz_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_F2yz_Px_Pz_S_c;
  abcd[2118] = 4.0E0*I_ERI_Fy2z_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Pz_S_c;
  abcd[2119] = 4.0E0*I_ERI_F3z_Fx2y_Pz_S_bc-2.0E0*1*I_ERI_F3z_Px_Pz_S_c;
  abcd[2120] = 4.0E0*I_ERI_F3x_Fxyz_Pz_S_bc;
  abcd[2121] = 4.0E0*I_ERI_F2xy_Fxyz_Pz_S_bc;
  abcd[2122] = 4.0E0*I_ERI_F2xz_Fxyz_Pz_S_bc;
  abcd[2123] = 4.0E0*I_ERI_Fx2y_Fxyz_Pz_S_bc;
  abcd[2124] = 4.0E0*I_ERI_Fxyz_Fxyz_Pz_S_bc;
  abcd[2125] = 4.0E0*I_ERI_Fx2z_Fxyz_Pz_S_bc;
  abcd[2126] = 4.0E0*I_ERI_F3y_Fxyz_Pz_S_bc;
  abcd[2127] = 4.0E0*I_ERI_F2yz_Fxyz_Pz_S_bc;
  abcd[2128] = 4.0E0*I_ERI_Fy2z_Fxyz_Pz_S_bc;
  abcd[2129] = 4.0E0*I_ERI_F3z_Fxyz_Pz_S_bc;
  abcd[2130] = 4.0E0*I_ERI_F3x_F3y_Pz_S_bc-2.0E0*2*I_ERI_F3x_Py_Pz_S_c;
  abcd[2131] = 4.0E0*I_ERI_F2xy_F3y_Pz_S_bc-2.0E0*2*I_ERI_F2xy_Py_Pz_S_c;
  abcd[2132] = 4.0E0*I_ERI_F2xz_F3y_Pz_S_bc-2.0E0*2*I_ERI_F2xz_Py_Pz_S_c;
  abcd[2133] = 4.0E0*I_ERI_Fx2y_F3y_Pz_S_bc-2.0E0*2*I_ERI_Fx2y_Py_Pz_S_c;
  abcd[2134] = 4.0E0*I_ERI_Fxyz_F3y_Pz_S_bc-2.0E0*2*I_ERI_Fxyz_Py_Pz_S_c;
  abcd[2135] = 4.0E0*I_ERI_Fx2z_F3y_Pz_S_bc-2.0E0*2*I_ERI_Fx2z_Py_Pz_S_c;
  abcd[2136] = 4.0E0*I_ERI_F3y_F3y_Pz_S_bc-2.0E0*2*I_ERI_F3y_Py_Pz_S_c;
  abcd[2137] = 4.0E0*I_ERI_F2yz_F3y_Pz_S_bc-2.0E0*2*I_ERI_F2yz_Py_Pz_S_c;
  abcd[2138] = 4.0E0*I_ERI_Fy2z_F3y_Pz_S_bc-2.0E0*2*I_ERI_Fy2z_Py_Pz_S_c;
  abcd[2139] = 4.0E0*I_ERI_F3z_F3y_Pz_S_bc-2.0E0*2*I_ERI_F3z_Py_Pz_S_c;
  abcd[2140] = 4.0E0*I_ERI_F3x_F2yz_Pz_S_bc-2.0E0*1*I_ERI_F3x_Pz_Pz_S_c;
  abcd[2141] = 4.0E0*I_ERI_F2xy_F2yz_Pz_S_bc-2.0E0*1*I_ERI_F2xy_Pz_Pz_S_c;
  abcd[2142] = 4.0E0*I_ERI_F2xz_F2yz_Pz_S_bc-2.0E0*1*I_ERI_F2xz_Pz_Pz_S_c;
  abcd[2143] = 4.0E0*I_ERI_Fx2y_F2yz_Pz_S_bc-2.0E0*1*I_ERI_Fx2y_Pz_Pz_S_c;
  abcd[2144] = 4.0E0*I_ERI_Fxyz_F2yz_Pz_S_bc-2.0E0*1*I_ERI_Fxyz_Pz_Pz_S_c;
  abcd[2145] = 4.0E0*I_ERI_Fx2z_F2yz_Pz_S_bc-2.0E0*1*I_ERI_Fx2z_Pz_Pz_S_c;
  abcd[2146] = 4.0E0*I_ERI_F3y_F2yz_Pz_S_bc-2.0E0*1*I_ERI_F3y_Pz_Pz_S_c;
  abcd[2147] = 4.0E0*I_ERI_F2yz_F2yz_Pz_S_bc-2.0E0*1*I_ERI_F2yz_Pz_Pz_S_c;
  abcd[2148] = 4.0E0*I_ERI_Fy2z_F2yz_Pz_S_bc-2.0E0*1*I_ERI_Fy2z_Pz_Pz_S_c;
  abcd[2149] = 4.0E0*I_ERI_F3z_F2yz_Pz_S_bc-2.0E0*1*I_ERI_F3z_Pz_Pz_S_c;
  abcd[2150] = 4.0E0*I_ERI_F3x_Fy2z_Pz_S_bc;
  abcd[2151] = 4.0E0*I_ERI_F2xy_Fy2z_Pz_S_bc;
  abcd[2152] = 4.0E0*I_ERI_F2xz_Fy2z_Pz_S_bc;
  abcd[2153] = 4.0E0*I_ERI_Fx2y_Fy2z_Pz_S_bc;
  abcd[2154] = 4.0E0*I_ERI_Fxyz_Fy2z_Pz_S_bc;
  abcd[2155] = 4.0E0*I_ERI_Fx2z_Fy2z_Pz_S_bc;
  abcd[2156] = 4.0E0*I_ERI_F3y_Fy2z_Pz_S_bc;
  abcd[2157] = 4.0E0*I_ERI_F2yz_Fy2z_Pz_S_bc;
  abcd[2158] = 4.0E0*I_ERI_Fy2z_Fy2z_Pz_S_bc;
  abcd[2159] = 4.0E0*I_ERI_F3z_Fy2z_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  abcd[2160] = 4.0E0*I_ERI_F3x_F2xz_Px_S_bc;
  abcd[2161] = 4.0E0*I_ERI_F2xy_F2xz_Px_S_bc;
  abcd[2162] = 4.0E0*I_ERI_F2xz_F2xz_Px_S_bc;
  abcd[2163] = 4.0E0*I_ERI_Fx2y_F2xz_Px_S_bc;
  abcd[2164] = 4.0E0*I_ERI_Fxyz_F2xz_Px_S_bc;
  abcd[2165] = 4.0E0*I_ERI_Fx2z_F2xz_Px_S_bc;
  abcd[2166] = 4.0E0*I_ERI_F3y_F2xz_Px_S_bc;
  abcd[2167] = 4.0E0*I_ERI_F2yz_F2xz_Px_S_bc;
  abcd[2168] = 4.0E0*I_ERI_Fy2z_F2xz_Px_S_bc;
  abcd[2169] = 4.0E0*I_ERI_F3z_F2xz_Px_S_bc;
  abcd[2170] = 4.0E0*I_ERI_F3x_Fxyz_Px_S_bc;
  abcd[2171] = 4.0E0*I_ERI_F2xy_Fxyz_Px_S_bc;
  abcd[2172] = 4.0E0*I_ERI_F2xz_Fxyz_Px_S_bc;
  abcd[2173] = 4.0E0*I_ERI_Fx2y_Fxyz_Px_S_bc;
  abcd[2174] = 4.0E0*I_ERI_Fxyz_Fxyz_Px_S_bc;
  abcd[2175] = 4.0E0*I_ERI_Fx2z_Fxyz_Px_S_bc;
  abcd[2176] = 4.0E0*I_ERI_F3y_Fxyz_Px_S_bc;
  abcd[2177] = 4.0E0*I_ERI_F2yz_Fxyz_Px_S_bc;
  abcd[2178] = 4.0E0*I_ERI_Fy2z_Fxyz_Px_S_bc;
  abcd[2179] = 4.0E0*I_ERI_F3z_Fxyz_Px_S_bc;
  abcd[2180] = 4.0E0*I_ERI_F3x_Fx2z_Px_S_bc-2.0E0*1*I_ERI_F3x_Px_Px_S_c;
  abcd[2181] = 4.0E0*I_ERI_F2xy_Fx2z_Px_S_bc-2.0E0*1*I_ERI_F2xy_Px_Px_S_c;
  abcd[2182] = 4.0E0*I_ERI_F2xz_Fx2z_Px_S_bc-2.0E0*1*I_ERI_F2xz_Px_Px_S_c;
  abcd[2183] = 4.0E0*I_ERI_Fx2y_Fx2z_Px_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Px_S_c;
  abcd[2184] = 4.0E0*I_ERI_Fxyz_Fx2z_Px_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Px_S_c;
  abcd[2185] = 4.0E0*I_ERI_Fx2z_Fx2z_Px_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Px_S_c;
  abcd[2186] = 4.0E0*I_ERI_F3y_Fx2z_Px_S_bc-2.0E0*1*I_ERI_F3y_Px_Px_S_c;
  abcd[2187] = 4.0E0*I_ERI_F2yz_Fx2z_Px_S_bc-2.0E0*1*I_ERI_F2yz_Px_Px_S_c;
  abcd[2188] = 4.0E0*I_ERI_Fy2z_Fx2z_Px_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Px_S_c;
  abcd[2189] = 4.0E0*I_ERI_F3z_Fx2z_Px_S_bc-2.0E0*1*I_ERI_F3z_Px_Px_S_c;
  abcd[2190] = 4.0E0*I_ERI_F3x_F2yz_Px_S_bc;
  abcd[2191] = 4.0E0*I_ERI_F2xy_F2yz_Px_S_bc;
  abcd[2192] = 4.0E0*I_ERI_F2xz_F2yz_Px_S_bc;
  abcd[2193] = 4.0E0*I_ERI_Fx2y_F2yz_Px_S_bc;
  abcd[2194] = 4.0E0*I_ERI_Fxyz_F2yz_Px_S_bc;
  abcd[2195] = 4.0E0*I_ERI_Fx2z_F2yz_Px_S_bc;
  abcd[2196] = 4.0E0*I_ERI_F3y_F2yz_Px_S_bc;
  abcd[2197] = 4.0E0*I_ERI_F2yz_F2yz_Px_S_bc;
  abcd[2198] = 4.0E0*I_ERI_Fy2z_F2yz_Px_S_bc;
  abcd[2199] = 4.0E0*I_ERI_F3z_F2yz_Px_S_bc;
  abcd[2200] = 4.0E0*I_ERI_F3x_Fy2z_Px_S_bc-2.0E0*1*I_ERI_F3x_Py_Px_S_c;
  abcd[2201] = 4.0E0*I_ERI_F2xy_Fy2z_Px_S_bc-2.0E0*1*I_ERI_F2xy_Py_Px_S_c;
  abcd[2202] = 4.0E0*I_ERI_F2xz_Fy2z_Px_S_bc-2.0E0*1*I_ERI_F2xz_Py_Px_S_c;
  abcd[2203] = 4.0E0*I_ERI_Fx2y_Fy2z_Px_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Px_S_c;
  abcd[2204] = 4.0E0*I_ERI_Fxyz_Fy2z_Px_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Px_S_c;
  abcd[2205] = 4.0E0*I_ERI_Fx2z_Fy2z_Px_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Px_S_c;
  abcd[2206] = 4.0E0*I_ERI_F3y_Fy2z_Px_S_bc-2.0E0*1*I_ERI_F3y_Py_Px_S_c;
  abcd[2207] = 4.0E0*I_ERI_F2yz_Fy2z_Px_S_bc-2.0E0*1*I_ERI_F2yz_Py_Px_S_c;
  abcd[2208] = 4.0E0*I_ERI_Fy2z_Fy2z_Px_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Px_S_c;
  abcd[2209] = 4.0E0*I_ERI_F3z_Fy2z_Px_S_bc-2.0E0*1*I_ERI_F3z_Py_Px_S_c;
  abcd[2210] = 4.0E0*I_ERI_F3x_F3z_Px_S_bc-2.0E0*2*I_ERI_F3x_Pz_Px_S_c;
  abcd[2211] = 4.0E0*I_ERI_F2xy_F3z_Px_S_bc-2.0E0*2*I_ERI_F2xy_Pz_Px_S_c;
  abcd[2212] = 4.0E0*I_ERI_F2xz_F3z_Px_S_bc-2.0E0*2*I_ERI_F2xz_Pz_Px_S_c;
  abcd[2213] = 4.0E0*I_ERI_Fx2y_F3z_Px_S_bc-2.0E0*2*I_ERI_Fx2y_Pz_Px_S_c;
  abcd[2214] = 4.0E0*I_ERI_Fxyz_F3z_Px_S_bc-2.0E0*2*I_ERI_Fxyz_Pz_Px_S_c;
  abcd[2215] = 4.0E0*I_ERI_Fx2z_F3z_Px_S_bc-2.0E0*2*I_ERI_Fx2z_Pz_Px_S_c;
  abcd[2216] = 4.0E0*I_ERI_F3y_F3z_Px_S_bc-2.0E0*2*I_ERI_F3y_Pz_Px_S_c;
  abcd[2217] = 4.0E0*I_ERI_F2yz_F3z_Px_S_bc-2.0E0*2*I_ERI_F2yz_Pz_Px_S_c;
  abcd[2218] = 4.0E0*I_ERI_Fy2z_F3z_Px_S_bc-2.0E0*2*I_ERI_Fy2z_Pz_Px_S_c;
  abcd[2219] = 4.0E0*I_ERI_F3z_F3z_Px_S_bc-2.0E0*2*I_ERI_F3z_Pz_Px_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  abcd[2220] = 4.0E0*I_ERI_F3x_F2xz_Py_S_bc;
  abcd[2221] = 4.0E0*I_ERI_F2xy_F2xz_Py_S_bc;
  abcd[2222] = 4.0E0*I_ERI_F2xz_F2xz_Py_S_bc;
  abcd[2223] = 4.0E0*I_ERI_Fx2y_F2xz_Py_S_bc;
  abcd[2224] = 4.0E0*I_ERI_Fxyz_F2xz_Py_S_bc;
  abcd[2225] = 4.0E0*I_ERI_Fx2z_F2xz_Py_S_bc;
  abcd[2226] = 4.0E0*I_ERI_F3y_F2xz_Py_S_bc;
  abcd[2227] = 4.0E0*I_ERI_F2yz_F2xz_Py_S_bc;
  abcd[2228] = 4.0E0*I_ERI_Fy2z_F2xz_Py_S_bc;
  abcd[2229] = 4.0E0*I_ERI_F3z_F2xz_Py_S_bc;
  abcd[2230] = 4.0E0*I_ERI_F3x_Fxyz_Py_S_bc;
  abcd[2231] = 4.0E0*I_ERI_F2xy_Fxyz_Py_S_bc;
  abcd[2232] = 4.0E0*I_ERI_F2xz_Fxyz_Py_S_bc;
  abcd[2233] = 4.0E0*I_ERI_Fx2y_Fxyz_Py_S_bc;
  abcd[2234] = 4.0E0*I_ERI_Fxyz_Fxyz_Py_S_bc;
  abcd[2235] = 4.0E0*I_ERI_Fx2z_Fxyz_Py_S_bc;
  abcd[2236] = 4.0E0*I_ERI_F3y_Fxyz_Py_S_bc;
  abcd[2237] = 4.0E0*I_ERI_F2yz_Fxyz_Py_S_bc;
  abcd[2238] = 4.0E0*I_ERI_Fy2z_Fxyz_Py_S_bc;
  abcd[2239] = 4.0E0*I_ERI_F3z_Fxyz_Py_S_bc;
  abcd[2240] = 4.0E0*I_ERI_F3x_Fx2z_Py_S_bc-2.0E0*1*I_ERI_F3x_Px_Py_S_c;
  abcd[2241] = 4.0E0*I_ERI_F2xy_Fx2z_Py_S_bc-2.0E0*1*I_ERI_F2xy_Px_Py_S_c;
  abcd[2242] = 4.0E0*I_ERI_F2xz_Fx2z_Py_S_bc-2.0E0*1*I_ERI_F2xz_Px_Py_S_c;
  abcd[2243] = 4.0E0*I_ERI_Fx2y_Fx2z_Py_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Py_S_c;
  abcd[2244] = 4.0E0*I_ERI_Fxyz_Fx2z_Py_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Py_S_c;
  abcd[2245] = 4.0E0*I_ERI_Fx2z_Fx2z_Py_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Py_S_c;
  abcd[2246] = 4.0E0*I_ERI_F3y_Fx2z_Py_S_bc-2.0E0*1*I_ERI_F3y_Px_Py_S_c;
  abcd[2247] = 4.0E0*I_ERI_F2yz_Fx2z_Py_S_bc-2.0E0*1*I_ERI_F2yz_Px_Py_S_c;
  abcd[2248] = 4.0E0*I_ERI_Fy2z_Fx2z_Py_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Py_S_c;
  abcd[2249] = 4.0E0*I_ERI_F3z_Fx2z_Py_S_bc-2.0E0*1*I_ERI_F3z_Px_Py_S_c;
  abcd[2250] = 4.0E0*I_ERI_F3x_F2yz_Py_S_bc;
  abcd[2251] = 4.0E0*I_ERI_F2xy_F2yz_Py_S_bc;
  abcd[2252] = 4.0E0*I_ERI_F2xz_F2yz_Py_S_bc;
  abcd[2253] = 4.0E0*I_ERI_Fx2y_F2yz_Py_S_bc;
  abcd[2254] = 4.0E0*I_ERI_Fxyz_F2yz_Py_S_bc;
  abcd[2255] = 4.0E0*I_ERI_Fx2z_F2yz_Py_S_bc;
  abcd[2256] = 4.0E0*I_ERI_F3y_F2yz_Py_S_bc;
  abcd[2257] = 4.0E0*I_ERI_F2yz_F2yz_Py_S_bc;
  abcd[2258] = 4.0E0*I_ERI_Fy2z_F2yz_Py_S_bc;
  abcd[2259] = 4.0E0*I_ERI_F3z_F2yz_Py_S_bc;
  abcd[2260] = 4.0E0*I_ERI_F3x_Fy2z_Py_S_bc-2.0E0*1*I_ERI_F3x_Py_Py_S_c;
  abcd[2261] = 4.0E0*I_ERI_F2xy_Fy2z_Py_S_bc-2.0E0*1*I_ERI_F2xy_Py_Py_S_c;
  abcd[2262] = 4.0E0*I_ERI_F2xz_Fy2z_Py_S_bc-2.0E0*1*I_ERI_F2xz_Py_Py_S_c;
  abcd[2263] = 4.0E0*I_ERI_Fx2y_Fy2z_Py_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Py_S_c;
  abcd[2264] = 4.0E0*I_ERI_Fxyz_Fy2z_Py_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Py_S_c;
  abcd[2265] = 4.0E0*I_ERI_Fx2z_Fy2z_Py_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Py_S_c;
  abcd[2266] = 4.0E0*I_ERI_F3y_Fy2z_Py_S_bc-2.0E0*1*I_ERI_F3y_Py_Py_S_c;
  abcd[2267] = 4.0E0*I_ERI_F2yz_Fy2z_Py_S_bc-2.0E0*1*I_ERI_F2yz_Py_Py_S_c;
  abcd[2268] = 4.0E0*I_ERI_Fy2z_Fy2z_Py_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Py_S_c;
  abcd[2269] = 4.0E0*I_ERI_F3z_Fy2z_Py_S_bc-2.0E0*1*I_ERI_F3z_Py_Py_S_c;
  abcd[2270] = 4.0E0*I_ERI_F3x_F3z_Py_S_bc-2.0E0*2*I_ERI_F3x_Pz_Py_S_c;
  abcd[2271] = 4.0E0*I_ERI_F2xy_F3z_Py_S_bc-2.0E0*2*I_ERI_F2xy_Pz_Py_S_c;
  abcd[2272] = 4.0E0*I_ERI_F2xz_F3z_Py_S_bc-2.0E0*2*I_ERI_F2xz_Pz_Py_S_c;
  abcd[2273] = 4.0E0*I_ERI_Fx2y_F3z_Py_S_bc-2.0E0*2*I_ERI_Fx2y_Pz_Py_S_c;
  abcd[2274] = 4.0E0*I_ERI_Fxyz_F3z_Py_S_bc-2.0E0*2*I_ERI_Fxyz_Pz_Py_S_c;
  abcd[2275] = 4.0E0*I_ERI_Fx2z_F3z_Py_S_bc-2.0E0*2*I_ERI_Fx2z_Pz_Py_S_c;
  abcd[2276] = 4.0E0*I_ERI_F3y_F3z_Py_S_bc-2.0E0*2*I_ERI_F3y_Pz_Py_S_c;
  abcd[2277] = 4.0E0*I_ERI_F2yz_F3z_Py_S_bc-2.0E0*2*I_ERI_F2yz_Pz_Py_S_c;
  abcd[2278] = 4.0E0*I_ERI_Fy2z_F3z_Py_S_bc-2.0E0*2*I_ERI_Fy2z_Pz_Py_S_c;
  abcd[2279] = 4.0E0*I_ERI_F3z_F3z_Py_S_bc-2.0E0*2*I_ERI_F3z_Pz_Py_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_P_S_bc
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  abcd[2280] = 4.0E0*I_ERI_F3x_F2xz_Pz_S_bc;
  abcd[2281] = 4.0E0*I_ERI_F2xy_F2xz_Pz_S_bc;
  abcd[2282] = 4.0E0*I_ERI_F2xz_F2xz_Pz_S_bc;
  abcd[2283] = 4.0E0*I_ERI_Fx2y_F2xz_Pz_S_bc;
  abcd[2284] = 4.0E0*I_ERI_Fxyz_F2xz_Pz_S_bc;
  abcd[2285] = 4.0E0*I_ERI_Fx2z_F2xz_Pz_S_bc;
  abcd[2286] = 4.0E0*I_ERI_F3y_F2xz_Pz_S_bc;
  abcd[2287] = 4.0E0*I_ERI_F2yz_F2xz_Pz_S_bc;
  abcd[2288] = 4.0E0*I_ERI_Fy2z_F2xz_Pz_S_bc;
  abcd[2289] = 4.0E0*I_ERI_F3z_F2xz_Pz_S_bc;
  abcd[2290] = 4.0E0*I_ERI_F3x_Fxyz_Pz_S_bc;
  abcd[2291] = 4.0E0*I_ERI_F2xy_Fxyz_Pz_S_bc;
  abcd[2292] = 4.0E0*I_ERI_F2xz_Fxyz_Pz_S_bc;
  abcd[2293] = 4.0E0*I_ERI_Fx2y_Fxyz_Pz_S_bc;
  abcd[2294] = 4.0E0*I_ERI_Fxyz_Fxyz_Pz_S_bc;
  abcd[2295] = 4.0E0*I_ERI_Fx2z_Fxyz_Pz_S_bc;
  abcd[2296] = 4.0E0*I_ERI_F3y_Fxyz_Pz_S_bc;
  abcd[2297] = 4.0E0*I_ERI_F2yz_Fxyz_Pz_S_bc;
  abcd[2298] = 4.0E0*I_ERI_Fy2z_Fxyz_Pz_S_bc;
  abcd[2299] = 4.0E0*I_ERI_F3z_Fxyz_Pz_S_bc;
  abcd[2300] = 4.0E0*I_ERI_F3x_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_F3x_Px_Pz_S_c;
  abcd[2301] = 4.0E0*I_ERI_F2xy_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_F2xy_Px_Pz_S_c;
  abcd[2302] = 4.0E0*I_ERI_F2xz_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_F2xz_Px_Pz_S_c;
  abcd[2303] = 4.0E0*I_ERI_Fx2y_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_Fx2y_Px_Pz_S_c;
  abcd[2304] = 4.0E0*I_ERI_Fxyz_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_Fxyz_Px_Pz_S_c;
  abcd[2305] = 4.0E0*I_ERI_Fx2z_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_Fx2z_Px_Pz_S_c;
  abcd[2306] = 4.0E0*I_ERI_F3y_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_F3y_Px_Pz_S_c;
  abcd[2307] = 4.0E0*I_ERI_F2yz_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_F2yz_Px_Pz_S_c;
  abcd[2308] = 4.0E0*I_ERI_Fy2z_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_Fy2z_Px_Pz_S_c;
  abcd[2309] = 4.0E0*I_ERI_F3z_Fx2z_Pz_S_bc-2.0E0*1*I_ERI_F3z_Px_Pz_S_c;
  abcd[2310] = 4.0E0*I_ERI_F3x_F2yz_Pz_S_bc;
  abcd[2311] = 4.0E0*I_ERI_F2xy_F2yz_Pz_S_bc;
  abcd[2312] = 4.0E0*I_ERI_F2xz_F2yz_Pz_S_bc;
  abcd[2313] = 4.0E0*I_ERI_Fx2y_F2yz_Pz_S_bc;
  abcd[2314] = 4.0E0*I_ERI_Fxyz_F2yz_Pz_S_bc;
  abcd[2315] = 4.0E0*I_ERI_Fx2z_F2yz_Pz_S_bc;
  abcd[2316] = 4.0E0*I_ERI_F3y_F2yz_Pz_S_bc;
  abcd[2317] = 4.0E0*I_ERI_F2yz_F2yz_Pz_S_bc;
  abcd[2318] = 4.0E0*I_ERI_Fy2z_F2yz_Pz_S_bc;
  abcd[2319] = 4.0E0*I_ERI_F3z_F2yz_Pz_S_bc;
  abcd[2320] = 4.0E0*I_ERI_F3x_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_F3x_Py_Pz_S_c;
  abcd[2321] = 4.0E0*I_ERI_F2xy_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_F2xy_Py_Pz_S_c;
  abcd[2322] = 4.0E0*I_ERI_F2xz_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_F2xz_Py_Pz_S_c;
  abcd[2323] = 4.0E0*I_ERI_Fx2y_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_Fx2y_Py_Pz_S_c;
  abcd[2324] = 4.0E0*I_ERI_Fxyz_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_Fxyz_Py_Pz_S_c;
  abcd[2325] = 4.0E0*I_ERI_Fx2z_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_Fx2z_Py_Pz_S_c;
  abcd[2326] = 4.0E0*I_ERI_F3y_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_F3y_Py_Pz_S_c;
  abcd[2327] = 4.0E0*I_ERI_F2yz_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_F2yz_Py_Pz_S_c;
  abcd[2328] = 4.0E0*I_ERI_Fy2z_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_Fy2z_Py_Pz_S_c;
  abcd[2329] = 4.0E0*I_ERI_F3z_Fy2z_Pz_S_bc-2.0E0*1*I_ERI_F3z_Py_Pz_S_c;
  abcd[2330] = 4.0E0*I_ERI_F3x_F3z_Pz_S_bc-2.0E0*2*I_ERI_F3x_Pz_Pz_S_c;
  abcd[2331] = 4.0E0*I_ERI_F2xy_F3z_Pz_S_bc-2.0E0*2*I_ERI_F2xy_Pz_Pz_S_c;
  abcd[2332] = 4.0E0*I_ERI_F2xz_F3z_Pz_S_bc-2.0E0*2*I_ERI_F2xz_Pz_Pz_S_c;
  abcd[2333] = 4.0E0*I_ERI_Fx2y_F3z_Pz_S_bc-2.0E0*2*I_ERI_Fx2y_Pz_Pz_S_c;
  abcd[2334] = 4.0E0*I_ERI_Fxyz_F3z_Pz_S_bc-2.0E0*2*I_ERI_Fxyz_Pz_Pz_S_c;
  abcd[2335] = 4.0E0*I_ERI_Fx2z_F3z_Pz_S_bc-2.0E0*2*I_ERI_Fx2z_Pz_Pz_S_c;
  abcd[2336] = 4.0E0*I_ERI_F3y_F3z_Pz_S_bc-2.0E0*2*I_ERI_F3y_Pz_Pz_S_c;
  abcd[2337] = 4.0E0*I_ERI_F2yz_F3z_Pz_S_bc-2.0E0*2*I_ERI_F2yz_Pz_Pz_S_c;
  abcd[2338] = 4.0E0*I_ERI_Fy2z_F3z_Pz_S_bc-2.0E0*2*I_ERI_Fy2z_Pz_Pz_S_c;
  abcd[2339] = 4.0E0*I_ERI_F3z_F3z_Pz_S_bc-2.0E0*2*I_ERI_F3z_Pz_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_cc
   * RHS shell quartet name: SQ_ERI_F_D_S_S_c
   ************************************************************/
  abcd[2340] = 4.0E0*I_ERI_F3x_D2x_D2x_S_cc-2.0E0*1*I_ERI_F3x_D2x_S_S_c;
  abcd[2341] = 4.0E0*I_ERI_F2xy_D2x_D2x_S_cc-2.0E0*1*I_ERI_F2xy_D2x_S_S_c;
  abcd[2342] = 4.0E0*I_ERI_F2xz_D2x_D2x_S_cc-2.0E0*1*I_ERI_F2xz_D2x_S_S_c;
  abcd[2343] = 4.0E0*I_ERI_Fx2y_D2x_D2x_S_cc-2.0E0*1*I_ERI_Fx2y_D2x_S_S_c;
  abcd[2344] = 4.0E0*I_ERI_Fxyz_D2x_D2x_S_cc-2.0E0*1*I_ERI_Fxyz_D2x_S_S_c;
  abcd[2345] = 4.0E0*I_ERI_Fx2z_D2x_D2x_S_cc-2.0E0*1*I_ERI_Fx2z_D2x_S_S_c;
  abcd[2346] = 4.0E0*I_ERI_F3y_D2x_D2x_S_cc-2.0E0*1*I_ERI_F3y_D2x_S_S_c;
  abcd[2347] = 4.0E0*I_ERI_F2yz_D2x_D2x_S_cc-2.0E0*1*I_ERI_F2yz_D2x_S_S_c;
  abcd[2348] = 4.0E0*I_ERI_Fy2z_D2x_D2x_S_cc-2.0E0*1*I_ERI_Fy2z_D2x_S_S_c;
  abcd[2349] = 4.0E0*I_ERI_F3z_D2x_D2x_S_cc-2.0E0*1*I_ERI_F3z_D2x_S_S_c;
  abcd[2350] = 4.0E0*I_ERI_F3x_Dxy_D2x_S_cc-2.0E0*1*I_ERI_F3x_Dxy_S_S_c;
  abcd[2351] = 4.0E0*I_ERI_F2xy_Dxy_D2x_S_cc-2.0E0*1*I_ERI_F2xy_Dxy_S_S_c;
  abcd[2352] = 4.0E0*I_ERI_F2xz_Dxy_D2x_S_cc-2.0E0*1*I_ERI_F2xz_Dxy_S_S_c;
  abcd[2353] = 4.0E0*I_ERI_Fx2y_Dxy_D2x_S_cc-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_c;
  abcd[2354] = 4.0E0*I_ERI_Fxyz_Dxy_D2x_S_cc-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_c;
  abcd[2355] = 4.0E0*I_ERI_Fx2z_Dxy_D2x_S_cc-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_c;
  abcd[2356] = 4.0E0*I_ERI_F3y_Dxy_D2x_S_cc-2.0E0*1*I_ERI_F3y_Dxy_S_S_c;
  abcd[2357] = 4.0E0*I_ERI_F2yz_Dxy_D2x_S_cc-2.0E0*1*I_ERI_F2yz_Dxy_S_S_c;
  abcd[2358] = 4.0E0*I_ERI_Fy2z_Dxy_D2x_S_cc-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_c;
  abcd[2359] = 4.0E0*I_ERI_F3z_Dxy_D2x_S_cc-2.0E0*1*I_ERI_F3z_Dxy_S_S_c;
  abcd[2360] = 4.0E0*I_ERI_F3x_Dxz_D2x_S_cc-2.0E0*1*I_ERI_F3x_Dxz_S_S_c;
  abcd[2361] = 4.0E0*I_ERI_F2xy_Dxz_D2x_S_cc-2.0E0*1*I_ERI_F2xy_Dxz_S_S_c;
  abcd[2362] = 4.0E0*I_ERI_F2xz_Dxz_D2x_S_cc-2.0E0*1*I_ERI_F2xz_Dxz_S_S_c;
  abcd[2363] = 4.0E0*I_ERI_Fx2y_Dxz_D2x_S_cc-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_c;
  abcd[2364] = 4.0E0*I_ERI_Fxyz_Dxz_D2x_S_cc-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_c;
  abcd[2365] = 4.0E0*I_ERI_Fx2z_Dxz_D2x_S_cc-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_c;
  abcd[2366] = 4.0E0*I_ERI_F3y_Dxz_D2x_S_cc-2.0E0*1*I_ERI_F3y_Dxz_S_S_c;
  abcd[2367] = 4.0E0*I_ERI_F2yz_Dxz_D2x_S_cc-2.0E0*1*I_ERI_F2yz_Dxz_S_S_c;
  abcd[2368] = 4.0E0*I_ERI_Fy2z_Dxz_D2x_S_cc-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_c;
  abcd[2369] = 4.0E0*I_ERI_F3z_Dxz_D2x_S_cc-2.0E0*1*I_ERI_F3z_Dxz_S_S_c;
  abcd[2370] = 4.0E0*I_ERI_F3x_D2y_D2x_S_cc-2.0E0*1*I_ERI_F3x_D2y_S_S_c;
  abcd[2371] = 4.0E0*I_ERI_F2xy_D2y_D2x_S_cc-2.0E0*1*I_ERI_F2xy_D2y_S_S_c;
  abcd[2372] = 4.0E0*I_ERI_F2xz_D2y_D2x_S_cc-2.0E0*1*I_ERI_F2xz_D2y_S_S_c;
  abcd[2373] = 4.0E0*I_ERI_Fx2y_D2y_D2x_S_cc-2.0E0*1*I_ERI_Fx2y_D2y_S_S_c;
  abcd[2374] = 4.0E0*I_ERI_Fxyz_D2y_D2x_S_cc-2.0E0*1*I_ERI_Fxyz_D2y_S_S_c;
  abcd[2375] = 4.0E0*I_ERI_Fx2z_D2y_D2x_S_cc-2.0E0*1*I_ERI_Fx2z_D2y_S_S_c;
  abcd[2376] = 4.0E0*I_ERI_F3y_D2y_D2x_S_cc-2.0E0*1*I_ERI_F3y_D2y_S_S_c;
  abcd[2377] = 4.0E0*I_ERI_F2yz_D2y_D2x_S_cc-2.0E0*1*I_ERI_F2yz_D2y_S_S_c;
  abcd[2378] = 4.0E0*I_ERI_Fy2z_D2y_D2x_S_cc-2.0E0*1*I_ERI_Fy2z_D2y_S_S_c;
  abcd[2379] = 4.0E0*I_ERI_F3z_D2y_D2x_S_cc-2.0E0*1*I_ERI_F3z_D2y_S_S_c;
  abcd[2380] = 4.0E0*I_ERI_F3x_Dyz_D2x_S_cc-2.0E0*1*I_ERI_F3x_Dyz_S_S_c;
  abcd[2381] = 4.0E0*I_ERI_F2xy_Dyz_D2x_S_cc-2.0E0*1*I_ERI_F2xy_Dyz_S_S_c;
  abcd[2382] = 4.0E0*I_ERI_F2xz_Dyz_D2x_S_cc-2.0E0*1*I_ERI_F2xz_Dyz_S_S_c;
  abcd[2383] = 4.0E0*I_ERI_Fx2y_Dyz_D2x_S_cc-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_c;
  abcd[2384] = 4.0E0*I_ERI_Fxyz_Dyz_D2x_S_cc-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_c;
  abcd[2385] = 4.0E0*I_ERI_Fx2z_Dyz_D2x_S_cc-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_c;
  abcd[2386] = 4.0E0*I_ERI_F3y_Dyz_D2x_S_cc-2.0E0*1*I_ERI_F3y_Dyz_S_S_c;
  abcd[2387] = 4.0E0*I_ERI_F2yz_Dyz_D2x_S_cc-2.0E0*1*I_ERI_F2yz_Dyz_S_S_c;
  abcd[2388] = 4.0E0*I_ERI_Fy2z_Dyz_D2x_S_cc-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_c;
  abcd[2389] = 4.0E0*I_ERI_F3z_Dyz_D2x_S_cc-2.0E0*1*I_ERI_F3z_Dyz_S_S_c;
  abcd[2390] = 4.0E0*I_ERI_F3x_D2z_D2x_S_cc-2.0E0*1*I_ERI_F3x_D2z_S_S_c;
  abcd[2391] = 4.0E0*I_ERI_F2xy_D2z_D2x_S_cc-2.0E0*1*I_ERI_F2xy_D2z_S_S_c;
  abcd[2392] = 4.0E0*I_ERI_F2xz_D2z_D2x_S_cc-2.0E0*1*I_ERI_F2xz_D2z_S_S_c;
  abcd[2393] = 4.0E0*I_ERI_Fx2y_D2z_D2x_S_cc-2.0E0*1*I_ERI_Fx2y_D2z_S_S_c;
  abcd[2394] = 4.0E0*I_ERI_Fxyz_D2z_D2x_S_cc-2.0E0*1*I_ERI_Fxyz_D2z_S_S_c;
  abcd[2395] = 4.0E0*I_ERI_Fx2z_D2z_D2x_S_cc-2.0E0*1*I_ERI_Fx2z_D2z_S_S_c;
  abcd[2396] = 4.0E0*I_ERI_F3y_D2z_D2x_S_cc-2.0E0*1*I_ERI_F3y_D2z_S_S_c;
  abcd[2397] = 4.0E0*I_ERI_F2yz_D2z_D2x_S_cc-2.0E0*1*I_ERI_F2yz_D2z_S_S_c;
  abcd[2398] = 4.0E0*I_ERI_Fy2z_D2z_D2x_S_cc-2.0E0*1*I_ERI_Fy2z_D2z_S_S_c;
  abcd[2399] = 4.0E0*I_ERI_F3z_D2z_D2x_S_cc-2.0E0*1*I_ERI_F3z_D2z_S_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_cc
   * RHS shell quartet name: SQ_ERI_F_D_S_S_c
   ************************************************************/
  abcd[2400] = 4.0E0*I_ERI_F3x_D2x_Dxy_S_cc;
  abcd[2401] = 4.0E0*I_ERI_F2xy_D2x_Dxy_S_cc;
  abcd[2402] = 4.0E0*I_ERI_F2xz_D2x_Dxy_S_cc;
  abcd[2403] = 4.0E0*I_ERI_Fx2y_D2x_Dxy_S_cc;
  abcd[2404] = 4.0E0*I_ERI_Fxyz_D2x_Dxy_S_cc;
  abcd[2405] = 4.0E0*I_ERI_Fx2z_D2x_Dxy_S_cc;
  abcd[2406] = 4.0E0*I_ERI_F3y_D2x_Dxy_S_cc;
  abcd[2407] = 4.0E0*I_ERI_F2yz_D2x_Dxy_S_cc;
  abcd[2408] = 4.0E0*I_ERI_Fy2z_D2x_Dxy_S_cc;
  abcd[2409] = 4.0E0*I_ERI_F3z_D2x_Dxy_S_cc;
  abcd[2410] = 4.0E0*I_ERI_F3x_Dxy_Dxy_S_cc;
  abcd[2411] = 4.0E0*I_ERI_F2xy_Dxy_Dxy_S_cc;
  abcd[2412] = 4.0E0*I_ERI_F2xz_Dxy_Dxy_S_cc;
  abcd[2413] = 4.0E0*I_ERI_Fx2y_Dxy_Dxy_S_cc;
  abcd[2414] = 4.0E0*I_ERI_Fxyz_Dxy_Dxy_S_cc;
  abcd[2415] = 4.0E0*I_ERI_Fx2z_Dxy_Dxy_S_cc;
  abcd[2416] = 4.0E0*I_ERI_F3y_Dxy_Dxy_S_cc;
  abcd[2417] = 4.0E0*I_ERI_F2yz_Dxy_Dxy_S_cc;
  abcd[2418] = 4.0E0*I_ERI_Fy2z_Dxy_Dxy_S_cc;
  abcd[2419] = 4.0E0*I_ERI_F3z_Dxy_Dxy_S_cc;
  abcd[2420] = 4.0E0*I_ERI_F3x_Dxz_Dxy_S_cc;
  abcd[2421] = 4.0E0*I_ERI_F2xy_Dxz_Dxy_S_cc;
  abcd[2422] = 4.0E0*I_ERI_F2xz_Dxz_Dxy_S_cc;
  abcd[2423] = 4.0E0*I_ERI_Fx2y_Dxz_Dxy_S_cc;
  abcd[2424] = 4.0E0*I_ERI_Fxyz_Dxz_Dxy_S_cc;
  abcd[2425] = 4.0E0*I_ERI_Fx2z_Dxz_Dxy_S_cc;
  abcd[2426] = 4.0E0*I_ERI_F3y_Dxz_Dxy_S_cc;
  abcd[2427] = 4.0E0*I_ERI_F2yz_Dxz_Dxy_S_cc;
  abcd[2428] = 4.0E0*I_ERI_Fy2z_Dxz_Dxy_S_cc;
  abcd[2429] = 4.0E0*I_ERI_F3z_Dxz_Dxy_S_cc;
  abcd[2430] = 4.0E0*I_ERI_F3x_D2y_Dxy_S_cc;
  abcd[2431] = 4.0E0*I_ERI_F2xy_D2y_Dxy_S_cc;
  abcd[2432] = 4.0E0*I_ERI_F2xz_D2y_Dxy_S_cc;
  abcd[2433] = 4.0E0*I_ERI_Fx2y_D2y_Dxy_S_cc;
  abcd[2434] = 4.0E0*I_ERI_Fxyz_D2y_Dxy_S_cc;
  abcd[2435] = 4.0E0*I_ERI_Fx2z_D2y_Dxy_S_cc;
  abcd[2436] = 4.0E0*I_ERI_F3y_D2y_Dxy_S_cc;
  abcd[2437] = 4.0E0*I_ERI_F2yz_D2y_Dxy_S_cc;
  abcd[2438] = 4.0E0*I_ERI_Fy2z_D2y_Dxy_S_cc;
  abcd[2439] = 4.0E0*I_ERI_F3z_D2y_Dxy_S_cc;
  abcd[2440] = 4.0E0*I_ERI_F3x_Dyz_Dxy_S_cc;
  abcd[2441] = 4.0E0*I_ERI_F2xy_Dyz_Dxy_S_cc;
  abcd[2442] = 4.0E0*I_ERI_F2xz_Dyz_Dxy_S_cc;
  abcd[2443] = 4.0E0*I_ERI_Fx2y_Dyz_Dxy_S_cc;
  abcd[2444] = 4.0E0*I_ERI_Fxyz_Dyz_Dxy_S_cc;
  abcd[2445] = 4.0E0*I_ERI_Fx2z_Dyz_Dxy_S_cc;
  abcd[2446] = 4.0E0*I_ERI_F3y_Dyz_Dxy_S_cc;
  abcd[2447] = 4.0E0*I_ERI_F2yz_Dyz_Dxy_S_cc;
  abcd[2448] = 4.0E0*I_ERI_Fy2z_Dyz_Dxy_S_cc;
  abcd[2449] = 4.0E0*I_ERI_F3z_Dyz_Dxy_S_cc;
  abcd[2450] = 4.0E0*I_ERI_F3x_D2z_Dxy_S_cc;
  abcd[2451] = 4.0E0*I_ERI_F2xy_D2z_Dxy_S_cc;
  abcd[2452] = 4.0E0*I_ERI_F2xz_D2z_Dxy_S_cc;
  abcd[2453] = 4.0E0*I_ERI_Fx2y_D2z_Dxy_S_cc;
  abcd[2454] = 4.0E0*I_ERI_Fxyz_D2z_Dxy_S_cc;
  abcd[2455] = 4.0E0*I_ERI_Fx2z_D2z_Dxy_S_cc;
  abcd[2456] = 4.0E0*I_ERI_F3y_D2z_Dxy_S_cc;
  abcd[2457] = 4.0E0*I_ERI_F2yz_D2z_Dxy_S_cc;
  abcd[2458] = 4.0E0*I_ERI_Fy2z_D2z_Dxy_S_cc;
  abcd[2459] = 4.0E0*I_ERI_F3z_D2z_Dxy_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_cc
   * RHS shell quartet name: SQ_ERI_F_D_S_S_c
   ************************************************************/
  abcd[2460] = 4.0E0*I_ERI_F3x_D2x_Dxz_S_cc;
  abcd[2461] = 4.0E0*I_ERI_F2xy_D2x_Dxz_S_cc;
  abcd[2462] = 4.0E0*I_ERI_F2xz_D2x_Dxz_S_cc;
  abcd[2463] = 4.0E0*I_ERI_Fx2y_D2x_Dxz_S_cc;
  abcd[2464] = 4.0E0*I_ERI_Fxyz_D2x_Dxz_S_cc;
  abcd[2465] = 4.0E0*I_ERI_Fx2z_D2x_Dxz_S_cc;
  abcd[2466] = 4.0E0*I_ERI_F3y_D2x_Dxz_S_cc;
  abcd[2467] = 4.0E0*I_ERI_F2yz_D2x_Dxz_S_cc;
  abcd[2468] = 4.0E0*I_ERI_Fy2z_D2x_Dxz_S_cc;
  abcd[2469] = 4.0E0*I_ERI_F3z_D2x_Dxz_S_cc;
  abcd[2470] = 4.0E0*I_ERI_F3x_Dxy_Dxz_S_cc;
  abcd[2471] = 4.0E0*I_ERI_F2xy_Dxy_Dxz_S_cc;
  abcd[2472] = 4.0E0*I_ERI_F2xz_Dxy_Dxz_S_cc;
  abcd[2473] = 4.0E0*I_ERI_Fx2y_Dxy_Dxz_S_cc;
  abcd[2474] = 4.0E0*I_ERI_Fxyz_Dxy_Dxz_S_cc;
  abcd[2475] = 4.0E0*I_ERI_Fx2z_Dxy_Dxz_S_cc;
  abcd[2476] = 4.0E0*I_ERI_F3y_Dxy_Dxz_S_cc;
  abcd[2477] = 4.0E0*I_ERI_F2yz_Dxy_Dxz_S_cc;
  abcd[2478] = 4.0E0*I_ERI_Fy2z_Dxy_Dxz_S_cc;
  abcd[2479] = 4.0E0*I_ERI_F3z_Dxy_Dxz_S_cc;
  abcd[2480] = 4.0E0*I_ERI_F3x_Dxz_Dxz_S_cc;
  abcd[2481] = 4.0E0*I_ERI_F2xy_Dxz_Dxz_S_cc;
  abcd[2482] = 4.0E0*I_ERI_F2xz_Dxz_Dxz_S_cc;
  abcd[2483] = 4.0E0*I_ERI_Fx2y_Dxz_Dxz_S_cc;
  abcd[2484] = 4.0E0*I_ERI_Fxyz_Dxz_Dxz_S_cc;
  abcd[2485] = 4.0E0*I_ERI_Fx2z_Dxz_Dxz_S_cc;
  abcd[2486] = 4.0E0*I_ERI_F3y_Dxz_Dxz_S_cc;
  abcd[2487] = 4.0E0*I_ERI_F2yz_Dxz_Dxz_S_cc;
  abcd[2488] = 4.0E0*I_ERI_Fy2z_Dxz_Dxz_S_cc;
  abcd[2489] = 4.0E0*I_ERI_F3z_Dxz_Dxz_S_cc;
  abcd[2490] = 4.0E0*I_ERI_F3x_D2y_Dxz_S_cc;
  abcd[2491] = 4.0E0*I_ERI_F2xy_D2y_Dxz_S_cc;
  abcd[2492] = 4.0E0*I_ERI_F2xz_D2y_Dxz_S_cc;
  abcd[2493] = 4.0E0*I_ERI_Fx2y_D2y_Dxz_S_cc;
  abcd[2494] = 4.0E0*I_ERI_Fxyz_D2y_Dxz_S_cc;
  abcd[2495] = 4.0E0*I_ERI_Fx2z_D2y_Dxz_S_cc;
  abcd[2496] = 4.0E0*I_ERI_F3y_D2y_Dxz_S_cc;
  abcd[2497] = 4.0E0*I_ERI_F2yz_D2y_Dxz_S_cc;
  abcd[2498] = 4.0E0*I_ERI_Fy2z_D2y_Dxz_S_cc;
  abcd[2499] = 4.0E0*I_ERI_F3z_D2y_Dxz_S_cc;
  abcd[2500] = 4.0E0*I_ERI_F3x_Dyz_Dxz_S_cc;
  abcd[2501] = 4.0E0*I_ERI_F2xy_Dyz_Dxz_S_cc;
  abcd[2502] = 4.0E0*I_ERI_F2xz_Dyz_Dxz_S_cc;
  abcd[2503] = 4.0E0*I_ERI_Fx2y_Dyz_Dxz_S_cc;
  abcd[2504] = 4.0E0*I_ERI_Fxyz_Dyz_Dxz_S_cc;
  abcd[2505] = 4.0E0*I_ERI_Fx2z_Dyz_Dxz_S_cc;
  abcd[2506] = 4.0E0*I_ERI_F3y_Dyz_Dxz_S_cc;
  abcd[2507] = 4.0E0*I_ERI_F2yz_Dyz_Dxz_S_cc;
  abcd[2508] = 4.0E0*I_ERI_Fy2z_Dyz_Dxz_S_cc;
  abcd[2509] = 4.0E0*I_ERI_F3z_Dyz_Dxz_S_cc;
  abcd[2510] = 4.0E0*I_ERI_F3x_D2z_Dxz_S_cc;
  abcd[2511] = 4.0E0*I_ERI_F2xy_D2z_Dxz_S_cc;
  abcd[2512] = 4.0E0*I_ERI_F2xz_D2z_Dxz_S_cc;
  abcd[2513] = 4.0E0*I_ERI_Fx2y_D2z_Dxz_S_cc;
  abcd[2514] = 4.0E0*I_ERI_Fxyz_D2z_Dxz_S_cc;
  abcd[2515] = 4.0E0*I_ERI_Fx2z_D2z_Dxz_S_cc;
  abcd[2516] = 4.0E0*I_ERI_F3y_D2z_Dxz_S_cc;
  abcd[2517] = 4.0E0*I_ERI_F2yz_D2z_Dxz_S_cc;
  abcd[2518] = 4.0E0*I_ERI_Fy2z_D2z_Dxz_S_cc;
  abcd[2519] = 4.0E0*I_ERI_F3z_D2z_Dxz_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_cc
   * RHS shell quartet name: SQ_ERI_F_D_S_S_c
   ************************************************************/
  abcd[2520] = 4.0E0*I_ERI_F3x_D2x_D2y_S_cc-2.0E0*1*I_ERI_F3x_D2x_S_S_c;
  abcd[2521] = 4.0E0*I_ERI_F2xy_D2x_D2y_S_cc-2.0E0*1*I_ERI_F2xy_D2x_S_S_c;
  abcd[2522] = 4.0E0*I_ERI_F2xz_D2x_D2y_S_cc-2.0E0*1*I_ERI_F2xz_D2x_S_S_c;
  abcd[2523] = 4.0E0*I_ERI_Fx2y_D2x_D2y_S_cc-2.0E0*1*I_ERI_Fx2y_D2x_S_S_c;
  abcd[2524] = 4.0E0*I_ERI_Fxyz_D2x_D2y_S_cc-2.0E0*1*I_ERI_Fxyz_D2x_S_S_c;
  abcd[2525] = 4.0E0*I_ERI_Fx2z_D2x_D2y_S_cc-2.0E0*1*I_ERI_Fx2z_D2x_S_S_c;
  abcd[2526] = 4.0E0*I_ERI_F3y_D2x_D2y_S_cc-2.0E0*1*I_ERI_F3y_D2x_S_S_c;
  abcd[2527] = 4.0E0*I_ERI_F2yz_D2x_D2y_S_cc-2.0E0*1*I_ERI_F2yz_D2x_S_S_c;
  abcd[2528] = 4.0E0*I_ERI_Fy2z_D2x_D2y_S_cc-2.0E0*1*I_ERI_Fy2z_D2x_S_S_c;
  abcd[2529] = 4.0E0*I_ERI_F3z_D2x_D2y_S_cc-2.0E0*1*I_ERI_F3z_D2x_S_S_c;
  abcd[2530] = 4.0E0*I_ERI_F3x_Dxy_D2y_S_cc-2.0E0*1*I_ERI_F3x_Dxy_S_S_c;
  abcd[2531] = 4.0E0*I_ERI_F2xy_Dxy_D2y_S_cc-2.0E0*1*I_ERI_F2xy_Dxy_S_S_c;
  abcd[2532] = 4.0E0*I_ERI_F2xz_Dxy_D2y_S_cc-2.0E0*1*I_ERI_F2xz_Dxy_S_S_c;
  abcd[2533] = 4.0E0*I_ERI_Fx2y_Dxy_D2y_S_cc-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_c;
  abcd[2534] = 4.0E0*I_ERI_Fxyz_Dxy_D2y_S_cc-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_c;
  abcd[2535] = 4.0E0*I_ERI_Fx2z_Dxy_D2y_S_cc-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_c;
  abcd[2536] = 4.0E0*I_ERI_F3y_Dxy_D2y_S_cc-2.0E0*1*I_ERI_F3y_Dxy_S_S_c;
  abcd[2537] = 4.0E0*I_ERI_F2yz_Dxy_D2y_S_cc-2.0E0*1*I_ERI_F2yz_Dxy_S_S_c;
  abcd[2538] = 4.0E0*I_ERI_Fy2z_Dxy_D2y_S_cc-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_c;
  abcd[2539] = 4.0E0*I_ERI_F3z_Dxy_D2y_S_cc-2.0E0*1*I_ERI_F3z_Dxy_S_S_c;
  abcd[2540] = 4.0E0*I_ERI_F3x_Dxz_D2y_S_cc-2.0E0*1*I_ERI_F3x_Dxz_S_S_c;
  abcd[2541] = 4.0E0*I_ERI_F2xy_Dxz_D2y_S_cc-2.0E0*1*I_ERI_F2xy_Dxz_S_S_c;
  abcd[2542] = 4.0E0*I_ERI_F2xz_Dxz_D2y_S_cc-2.0E0*1*I_ERI_F2xz_Dxz_S_S_c;
  abcd[2543] = 4.0E0*I_ERI_Fx2y_Dxz_D2y_S_cc-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_c;
  abcd[2544] = 4.0E0*I_ERI_Fxyz_Dxz_D2y_S_cc-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_c;
  abcd[2545] = 4.0E0*I_ERI_Fx2z_Dxz_D2y_S_cc-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_c;
  abcd[2546] = 4.0E0*I_ERI_F3y_Dxz_D2y_S_cc-2.0E0*1*I_ERI_F3y_Dxz_S_S_c;
  abcd[2547] = 4.0E0*I_ERI_F2yz_Dxz_D2y_S_cc-2.0E0*1*I_ERI_F2yz_Dxz_S_S_c;
  abcd[2548] = 4.0E0*I_ERI_Fy2z_Dxz_D2y_S_cc-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_c;
  abcd[2549] = 4.0E0*I_ERI_F3z_Dxz_D2y_S_cc-2.0E0*1*I_ERI_F3z_Dxz_S_S_c;
  abcd[2550] = 4.0E0*I_ERI_F3x_D2y_D2y_S_cc-2.0E0*1*I_ERI_F3x_D2y_S_S_c;
  abcd[2551] = 4.0E0*I_ERI_F2xy_D2y_D2y_S_cc-2.0E0*1*I_ERI_F2xy_D2y_S_S_c;
  abcd[2552] = 4.0E0*I_ERI_F2xz_D2y_D2y_S_cc-2.0E0*1*I_ERI_F2xz_D2y_S_S_c;
  abcd[2553] = 4.0E0*I_ERI_Fx2y_D2y_D2y_S_cc-2.0E0*1*I_ERI_Fx2y_D2y_S_S_c;
  abcd[2554] = 4.0E0*I_ERI_Fxyz_D2y_D2y_S_cc-2.0E0*1*I_ERI_Fxyz_D2y_S_S_c;
  abcd[2555] = 4.0E0*I_ERI_Fx2z_D2y_D2y_S_cc-2.0E0*1*I_ERI_Fx2z_D2y_S_S_c;
  abcd[2556] = 4.0E0*I_ERI_F3y_D2y_D2y_S_cc-2.0E0*1*I_ERI_F3y_D2y_S_S_c;
  abcd[2557] = 4.0E0*I_ERI_F2yz_D2y_D2y_S_cc-2.0E0*1*I_ERI_F2yz_D2y_S_S_c;
  abcd[2558] = 4.0E0*I_ERI_Fy2z_D2y_D2y_S_cc-2.0E0*1*I_ERI_Fy2z_D2y_S_S_c;
  abcd[2559] = 4.0E0*I_ERI_F3z_D2y_D2y_S_cc-2.0E0*1*I_ERI_F3z_D2y_S_S_c;
  abcd[2560] = 4.0E0*I_ERI_F3x_Dyz_D2y_S_cc-2.0E0*1*I_ERI_F3x_Dyz_S_S_c;
  abcd[2561] = 4.0E0*I_ERI_F2xy_Dyz_D2y_S_cc-2.0E0*1*I_ERI_F2xy_Dyz_S_S_c;
  abcd[2562] = 4.0E0*I_ERI_F2xz_Dyz_D2y_S_cc-2.0E0*1*I_ERI_F2xz_Dyz_S_S_c;
  abcd[2563] = 4.0E0*I_ERI_Fx2y_Dyz_D2y_S_cc-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_c;
  abcd[2564] = 4.0E0*I_ERI_Fxyz_Dyz_D2y_S_cc-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_c;
  abcd[2565] = 4.0E0*I_ERI_Fx2z_Dyz_D2y_S_cc-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_c;
  abcd[2566] = 4.0E0*I_ERI_F3y_Dyz_D2y_S_cc-2.0E0*1*I_ERI_F3y_Dyz_S_S_c;
  abcd[2567] = 4.0E0*I_ERI_F2yz_Dyz_D2y_S_cc-2.0E0*1*I_ERI_F2yz_Dyz_S_S_c;
  abcd[2568] = 4.0E0*I_ERI_Fy2z_Dyz_D2y_S_cc-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_c;
  abcd[2569] = 4.0E0*I_ERI_F3z_Dyz_D2y_S_cc-2.0E0*1*I_ERI_F3z_Dyz_S_S_c;
  abcd[2570] = 4.0E0*I_ERI_F3x_D2z_D2y_S_cc-2.0E0*1*I_ERI_F3x_D2z_S_S_c;
  abcd[2571] = 4.0E0*I_ERI_F2xy_D2z_D2y_S_cc-2.0E0*1*I_ERI_F2xy_D2z_S_S_c;
  abcd[2572] = 4.0E0*I_ERI_F2xz_D2z_D2y_S_cc-2.0E0*1*I_ERI_F2xz_D2z_S_S_c;
  abcd[2573] = 4.0E0*I_ERI_Fx2y_D2z_D2y_S_cc-2.0E0*1*I_ERI_Fx2y_D2z_S_S_c;
  abcd[2574] = 4.0E0*I_ERI_Fxyz_D2z_D2y_S_cc-2.0E0*1*I_ERI_Fxyz_D2z_S_S_c;
  abcd[2575] = 4.0E0*I_ERI_Fx2z_D2z_D2y_S_cc-2.0E0*1*I_ERI_Fx2z_D2z_S_S_c;
  abcd[2576] = 4.0E0*I_ERI_F3y_D2z_D2y_S_cc-2.0E0*1*I_ERI_F3y_D2z_S_S_c;
  abcd[2577] = 4.0E0*I_ERI_F2yz_D2z_D2y_S_cc-2.0E0*1*I_ERI_F2yz_D2z_S_S_c;
  abcd[2578] = 4.0E0*I_ERI_Fy2z_D2z_D2y_S_cc-2.0E0*1*I_ERI_Fy2z_D2z_S_S_c;
  abcd[2579] = 4.0E0*I_ERI_F3z_D2z_D2y_S_cc-2.0E0*1*I_ERI_F3z_D2z_S_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_cc
   * RHS shell quartet name: SQ_ERI_F_D_S_S_c
   ************************************************************/
  abcd[2580] = 4.0E0*I_ERI_F3x_D2x_Dyz_S_cc;
  abcd[2581] = 4.0E0*I_ERI_F2xy_D2x_Dyz_S_cc;
  abcd[2582] = 4.0E0*I_ERI_F2xz_D2x_Dyz_S_cc;
  abcd[2583] = 4.0E0*I_ERI_Fx2y_D2x_Dyz_S_cc;
  abcd[2584] = 4.0E0*I_ERI_Fxyz_D2x_Dyz_S_cc;
  abcd[2585] = 4.0E0*I_ERI_Fx2z_D2x_Dyz_S_cc;
  abcd[2586] = 4.0E0*I_ERI_F3y_D2x_Dyz_S_cc;
  abcd[2587] = 4.0E0*I_ERI_F2yz_D2x_Dyz_S_cc;
  abcd[2588] = 4.0E0*I_ERI_Fy2z_D2x_Dyz_S_cc;
  abcd[2589] = 4.0E0*I_ERI_F3z_D2x_Dyz_S_cc;
  abcd[2590] = 4.0E0*I_ERI_F3x_Dxy_Dyz_S_cc;
  abcd[2591] = 4.0E0*I_ERI_F2xy_Dxy_Dyz_S_cc;
  abcd[2592] = 4.0E0*I_ERI_F2xz_Dxy_Dyz_S_cc;
  abcd[2593] = 4.0E0*I_ERI_Fx2y_Dxy_Dyz_S_cc;
  abcd[2594] = 4.0E0*I_ERI_Fxyz_Dxy_Dyz_S_cc;
  abcd[2595] = 4.0E0*I_ERI_Fx2z_Dxy_Dyz_S_cc;
  abcd[2596] = 4.0E0*I_ERI_F3y_Dxy_Dyz_S_cc;
  abcd[2597] = 4.0E0*I_ERI_F2yz_Dxy_Dyz_S_cc;
  abcd[2598] = 4.0E0*I_ERI_Fy2z_Dxy_Dyz_S_cc;
  abcd[2599] = 4.0E0*I_ERI_F3z_Dxy_Dyz_S_cc;
  abcd[2600] = 4.0E0*I_ERI_F3x_Dxz_Dyz_S_cc;
  abcd[2601] = 4.0E0*I_ERI_F2xy_Dxz_Dyz_S_cc;
  abcd[2602] = 4.0E0*I_ERI_F2xz_Dxz_Dyz_S_cc;
  abcd[2603] = 4.0E0*I_ERI_Fx2y_Dxz_Dyz_S_cc;
  abcd[2604] = 4.0E0*I_ERI_Fxyz_Dxz_Dyz_S_cc;
  abcd[2605] = 4.0E0*I_ERI_Fx2z_Dxz_Dyz_S_cc;
  abcd[2606] = 4.0E0*I_ERI_F3y_Dxz_Dyz_S_cc;
  abcd[2607] = 4.0E0*I_ERI_F2yz_Dxz_Dyz_S_cc;
  abcd[2608] = 4.0E0*I_ERI_Fy2z_Dxz_Dyz_S_cc;
  abcd[2609] = 4.0E0*I_ERI_F3z_Dxz_Dyz_S_cc;
  abcd[2610] = 4.0E0*I_ERI_F3x_D2y_Dyz_S_cc;
  abcd[2611] = 4.0E0*I_ERI_F2xy_D2y_Dyz_S_cc;
  abcd[2612] = 4.0E0*I_ERI_F2xz_D2y_Dyz_S_cc;
  abcd[2613] = 4.0E0*I_ERI_Fx2y_D2y_Dyz_S_cc;
  abcd[2614] = 4.0E0*I_ERI_Fxyz_D2y_Dyz_S_cc;
  abcd[2615] = 4.0E0*I_ERI_Fx2z_D2y_Dyz_S_cc;
  abcd[2616] = 4.0E0*I_ERI_F3y_D2y_Dyz_S_cc;
  abcd[2617] = 4.0E0*I_ERI_F2yz_D2y_Dyz_S_cc;
  abcd[2618] = 4.0E0*I_ERI_Fy2z_D2y_Dyz_S_cc;
  abcd[2619] = 4.0E0*I_ERI_F3z_D2y_Dyz_S_cc;
  abcd[2620] = 4.0E0*I_ERI_F3x_Dyz_Dyz_S_cc;
  abcd[2621] = 4.0E0*I_ERI_F2xy_Dyz_Dyz_S_cc;
  abcd[2622] = 4.0E0*I_ERI_F2xz_Dyz_Dyz_S_cc;
  abcd[2623] = 4.0E0*I_ERI_Fx2y_Dyz_Dyz_S_cc;
  abcd[2624] = 4.0E0*I_ERI_Fxyz_Dyz_Dyz_S_cc;
  abcd[2625] = 4.0E0*I_ERI_Fx2z_Dyz_Dyz_S_cc;
  abcd[2626] = 4.0E0*I_ERI_F3y_Dyz_Dyz_S_cc;
  abcd[2627] = 4.0E0*I_ERI_F2yz_Dyz_Dyz_S_cc;
  abcd[2628] = 4.0E0*I_ERI_Fy2z_Dyz_Dyz_S_cc;
  abcd[2629] = 4.0E0*I_ERI_F3z_Dyz_Dyz_S_cc;
  abcd[2630] = 4.0E0*I_ERI_F3x_D2z_Dyz_S_cc;
  abcd[2631] = 4.0E0*I_ERI_F2xy_D2z_Dyz_S_cc;
  abcd[2632] = 4.0E0*I_ERI_F2xz_D2z_Dyz_S_cc;
  abcd[2633] = 4.0E0*I_ERI_Fx2y_D2z_Dyz_S_cc;
  abcd[2634] = 4.0E0*I_ERI_Fxyz_D2z_Dyz_S_cc;
  abcd[2635] = 4.0E0*I_ERI_Fx2z_D2z_Dyz_S_cc;
  abcd[2636] = 4.0E0*I_ERI_F3y_D2z_Dyz_S_cc;
  abcd[2637] = 4.0E0*I_ERI_F2yz_D2z_Dyz_S_cc;
  abcd[2638] = 4.0E0*I_ERI_Fy2z_D2z_Dyz_S_cc;
  abcd[2639] = 4.0E0*I_ERI_F3z_D2z_Dyz_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_D_S_cc
   * RHS shell quartet name: SQ_ERI_F_D_S_S_c
   ************************************************************/
  abcd[2640] = 4.0E0*I_ERI_F3x_D2x_D2z_S_cc-2.0E0*1*I_ERI_F3x_D2x_S_S_c;
  abcd[2641] = 4.0E0*I_ERI_F2xy_D2x_D2z_S_cc-2.0E0*1*I_ERI_F2xy_D2x_S_S_c;
  abcd[2642] = 4.0E0*I_ERI_F2xz_D2x_D2z_S_cc-2.0E0*1*I_ERI_F2xz_D2x_S_S_c;
  abcd[2643] = 4.0E0*I_ERI_Fx2y_D2x_D2z_S_cc-2.0E0*1*I_ERI_Fx2y_D2x_S_S_c;
  abcd[2644] = 4.0E0*I_ERI_Fxyz_D2x_D2z_S_cc-2.0E0*1*I_ERI_Fxyz_D2x_S_S_c;
  abcd[2645] = 4.0E0*I_ERI_Fx2z_D2x_D2z_S_cc-2.0E0*1*I_ERI_Fx2z_D2x_S_S_c;
  abcd[2646] = 4.0E0*I_ERI_F3y_D2x_D2z_S_cc-2.0E0*1*I_ERI_F3y_D2x_S_S_c;
  abcd[2647] = 4.0E0*I_ERI_F2yz_D2x_D2z_S_cc-2.0E0*1*I_ERI_F2yz_D2x_S_S_c;
  abcd[2648] = 4.0E0*I_ERI_Fy2z_D2x_D2z_S_cc-2.0E0*1*I_ERI_Fy2z_D2x_S_S_c;
  abcd[2649] = 4.0E0*I_ERI_F3z_D2x_D2z_S_cc-2.0E0*1*I_ERI_F3z_D2x_S_S_c;
  abcd[2650] = 4.0E0*I_ERI_F3x_Dxy_D2z_S_cc-2.0E0*1*I_ERI_F3x_Dxy_S_S_c;
  abcd[2651] = 4.0E0*I_ERI_F2xy_Dxy_D2z_S_cc-2.0E0*1*I_ERI_F2xy_Dxy_S_S_c;
  abcd[2652] = 4.0E0*I_ERI_F2xz_Dxy_D2z_S_cc-2.0E0*1*I_ERI_F2xz_Dxy_S_S_c;
  abcd[2653] = 4.0E0*I_ERI_Fx2y_Dxy_D2z_S_cc-2.0E0*1*I_ERI_Fx2y_Dxy_S_S_c;
  abcd[2654] = 4.0E0*I_ERI_Fxyz_Dxy_D2z_S_cc-2.0E0*1*I_ERI_Fxyz_Dxy_S_S_c;
  abcd[2655] = 4.0E0*I_ERI_Fx2z_Dxy_D2z_S_cc-2.0E0*1*I_ERI_Fx2z_Dxy_S_S_c;
  abcd[2656] = 4.0E0*I_ERI_F3y_Dxy_D2z_S_cc-2.0E0*1*I_ERI_F3y_Dxy_S_S_c;
  abcd[2657] = 4.0E0*I_ERI_F2yz_Dxy_D2z_S_cc-2.0E0*1*I_ERI_F2yz_Dxy_S_S_c;
  abcd[2658] = 4.0E0*I_ERI_Fy2z_Dxy_D2z_S_cc-2.0E0*1*I_ERI_Fy2z_Dxy_S_S_c;
  abcd[2659] = 4.0E0*I_ERI_F3z_Dxy_D2z_S_cc-2.0E0*1*I_ERI_F3z_Dxy_S_S_c;
  abcd[2660] = 4.0E0*I_ERI_F3x_Dxz_D2z_S_cc-2.0E0*1*I_ERI_F3x_Dxz_S_S_c;
  abcd[2661] = 4.0E0*I_ERI_F2xy_Dxz_D2z_S_cc-2.0E0*1*I_ERI_F2xy_Dxz_S_S_c;
  abcd[2662] = 4.0E0*I_ERI_F2xz_Dxz_D2z_S_cc-2.0E0*1*I_ERI_F2xz_Dxz_S_S_c;
  abcd[2663] = 4.0E0*I_ERI_Fx2y_Dxz_D2z_S_cc-2.0E0*1*I_ERI_Fx2y_Dxz_S_S_c;
  abcd[2664] = 4.0E0*I_ERI_Fxyz_Dxz_D2z_S_cc-2.0E0*1*I_ERI_Fxyz_Dxz_S_S_c;
  abcd[2665] = 4.0E0*I_ERI_Fx2z_Dxz_D2z_S_cc-2.0E0*1*I_ERI_Fx2z_Dxz_S_S_c;
  abcd[2666] = 4.0E0*I_ERI_F3y_Dxz_D2z_S_cc-2.0E0*1*I_ERI_F3y_Dxz_S_S_c;
  abcd[2667] = 4.0E0*I_ERI_F2yz_Dxz_D2z_S_cc-2.0E0*1*I_ERI_F2yz_Dxz_S_S_c;
  abcd[2668] = 4.0E0*I_ERI_Fy2z_Dxz_D2z_S_cc-2.0E0*1*I_ERI_Fy2z_Dxz_S_S_c;
  abcd[2669] = 4.0E0*I_ERI_F3z_Dxz_D2z_S_cc-2.0E0*1*I_ERI_F3z_Dxz_S_S_c;
  abcd[2670] = 4.0E0*I_ERI_F3x_D2y_D2z_S_cc-2.0E0*1*I_ERI_F3x_D2y_S_S_c;
  abcd[2671] = 4.0E0*I_ERI_F2xy_D2y_D2z_S_cc-2.0E0*1*I_ERI_F2xy_D2y_S_S_c;
  abcd[2672] = 4.0E0*I_ERI_F2xz_D2y_D2z_S_cc-2.0E0*1*I_ERI_F2xz_D2y_S_S_c;
  abcd[2673] = 4.0E0*I_ERI_Fx2y_D2y_D2z_S_cc-2.0E0*1*I_ERI_Fx2y_D2y_S_S_c;
  abcd[2674] = 4.0E0*I_ERI_Fxyz_D2y_D2z_S_cc-2.0E0*1*I_ERI_Fxyz_D2y_S_S_c;
  abcd[2675] = 4.0E0*I_ERI_Fx2z_D2y_D2z_S_cc-2.0E0*1*I_ERI_Fx2z_D2y_S_S_c;
  abcd[2676] = 4.0E0*I_ERI_F3y_D2y_D2z_S_cc-2.0E0*1*I_ERI_F3y_D2y_S_S_c;
  abcd[2677] = 4.0E0*I_ERI_F2yz_D2y_D2z_S_cc-2.0E0*1*I_ERI_F2yz_D2y_S_S_c;
  abcd[2678] = 4.0E0*I_ERI_Fy2z_D2y_D2z_S_cc-2.0E0*1*I_ERI_Fy2z_D2y_S_S_c;
  abcd[2679] = 4.0E0*I_ERI_F3z_D2y_D2z_S_cc-2.0E0*1*I_ERI_F3z_D2y_S_S_c;
  abcd[2680] = 4.0E0*I_ERI_F3x_Dyz_D2z_S_cc-2.0E0*1*I_ERI_F3x_Dyz_S_S_c;
  abcd[2681] = 4.0E0*I_ERI_F2xy_Dyz_D2z_S_cc-2.0E0*1*I_ERI_F2xy_Dyz_S_S_c;
  abcd[2682] = 4.0E0*I_ERI_F2xz_Dyz_D2z_S_cc-2.0E0*1*I_ERI_F2xz_Dyz_S_S_c;
  abcd[2683] = 4.0E0*I_ERI_Fx2y_Dyz_D2z_S_cc-2.0E0*1*I_ERI_Fx2y_Dyz_S_S_c;
  abcd[2684] = 4.0E0*I_ERI_Fxyz_Dyz_D2z_S_cc-2.0E0*1*I_ERI_Fxyz_Dyz_S_S_c;
  abcd[2685] = 4.0E0*I_ERI_Fx2z_Dyz_D2z_S_cc-2.0E0*1*I_ERI_Fx2z_Dyz_S_S_c;
  abcd[2686] = 4.0E0*I_ERI_F3y_Dyz_D2z_S_cc-2.0E0*1*I_ERI_F3y_Dyz_S_S_c;
  abcd[2687] = 4.0E0*I_ERI_F2yz_Dyz_D2z_S_cc-2.0E0*1*I_ERI_F2yz_Dyz_S_S_c;
  abcd[2688] = 4.0E0*I_ERI_Fy2z_Dyz_D2z_S_cc-2.0E0*1*I_ERI_Fy2z_Dyz_S_S_c;
  abcd[2689] = 4.0E0*I_ERI_F3z_Dyz_D2z_S_cc-2.0E0*1*I_ERI_F3z_Dyz_S_S_c;
  abcd[2690] = 4.0E0*I_ERI_F3x_D2z_D2z_S_cc-2.0E0*1*I_ERI_F3x_D2z_S_S_c;
  abcd[2691] = 4.0E0*I_ERI_F2xy_D2z_D2z_S_cc-2.0E0*1*I_ERI_F2xy_D2z_S_S_c;
  abcd[2692] = 4.0E0*I_ERI_F2xz_D2z_D2z_S_cc-2.0E0*1*I_ERI_F2xz_D2z_S_S_c;
  abcd[2693] = 4.0E0*I_ERI_Fx2y_D2z_D2z_S_cc-2.0E0*1*I_ERI_Fx2y_D2z_S_S_c;
  abcd[2694] = 4.0E0*I_ERI_Fxyz_D2z_D2z_S_cc-2.0E0*1*I_ERI_Fxyz_D2z_S_S_c;
  abcd[2695] = 4.0E0*I_ERI_Fx2z_D2z_D2z_S_cc-2.0E0*1*I_ERI_Fx2z_D2z_S_S_c;
  abcd[2696] = 4.0E0*I_ERI_F3y_D2z_D2z_S_cc-2.0E0*1*I_ERI_F3y_D2z_S_S_c;
  abcd[2697] = 4.0E0*I_ERI_F2yz_D2z_D2z_S_cc-2.0E0*1*I_ERI_F2yz_D2z_S_S_c;
  abcd[2698] = 4.0E0*I_ERI_Fy2z_D2z_D2z_S_cc-2.0E0*1*I_ERI_Fy2z_D2z_S_S_c;
  abcd[2699] = 4.0E0*I_ERI_F3z_D2z_D2z_S_cc-2.0E0*1*I_ERI_F3z_D2z_S_S_c;
}
