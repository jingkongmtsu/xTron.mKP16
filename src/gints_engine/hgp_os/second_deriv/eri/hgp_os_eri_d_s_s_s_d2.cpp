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
// BRA1 as redundant position, total RHS integrals evaluated as: 23108
// BRA2 as redundant position, total RHS integrals evaluated as: 24603
// KET1 as redundant position, total RHS integrals evaluated as: 15119
// KET2 as redundant position, total RHS integrals evaluated as: 15119
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

void hgp_os_eri_d_s_s_s_d2(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_G4x_S_S_S_aa = 0.0E0;
  Double I_ERI_G3xy_S_S_S_aa = 0.0E0;
  Double I_ERI_G3xz_S_S_S_aa = 0.0E0;
  Double I_ERI_G2x2y_S_S_S_aa = 0.0E0;
  Double I_ERI_G2xyz_S_S_S_aa = 0.0E0;
  Double I_ERI_G2x2z_S_S_S_aa = 0.0E0;
  Double I_ERI_Gx3y_S_S_S_aa = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S_aa = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S_aa = 0.0E0;
  Double I_ERI_Gx3z_S_S_S_aa = 0.0E0;
  Double I_ERI_G4y_S_S_S_aa = 0.0E0;
  Double I_ERI_G3yz_S_S_S_aa = 0.0E0;
  Double I_ERI_G2y2z_S_S_S_aa = 0.0E0;
  Double I_ERI_Gy3z_S_S_S_aa = 0.0E0;
  Double I_ERI_G4z_S_S_S_aa = 0.0E0;
  Double I_ERI_D2x_S_S_S_a = 0.0E0;
  Double I_ERI_Dxy_S_S_S_a = 0.0E0;
  Double I_ERI_Dxz_S_S_S_a = 0.0E0;
  Double I_ERI_D2y_S_S_S_a = 0.0E0;
  Double I_ERI_Dyz_S_S_S_a = 0.0E0;
  Double I_ERI_D2z_S_S_S_a = 0.0E0;
  Double I_ERI_S_S_S_S = 0.0E0;
  Double I_ERI_F3x_S_Px_S_ac = 0.0E0;
  Double I_ERI_F2xy_S_Px_S_ac = 0.0E0;
  Double I_ERI_F2xz_S_Px_S_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_F3y_S_Px_S_ac = 0.0E0;
  Double I_ERI_F2yz_S_Px_S_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S_ac = 0.0E0;
  Double I_ERI_F3z_S_Px_S_ac = 0.0E0;
  Double I_ERI_F3x_S_Py_S_ac = 0.0E0;
  Double I_ERI_F2xy_S_Py_S_ac = 0.0E0;
  Double I_ERI_F2xz_S_Py_S_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_F3y_S_Py_S_ac = 0.0E0;
  Double I_ERI_F2yz_S_Py_S_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S_ac = 0.0E0;
  Double I_ERI_F3z_S_Py_S_ac = 0.0E0;
  Double I_ERI_F3x_S_Pz_S_ac = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S_ac = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_F3y_S_Pz_S_ac = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_F3z_S_Pz_S_ac = 0.0E0;
  Double I_ERI_Px_S_Px_S_c = 0.0E0;
  Double I_ERI_Py_S_Px_S_c = 0.0E0;
  Double I_ERI_Pz_S_Px_S_c = 0.0E0;
  Double I_ERI_Px_S_Py_S_c = 0.0E0;
  Double I_ERI_Py_S_Py_S_c = 0.0E0;
  Double I_ERI_Pz_S_Py_S_c = 0.0E0;
  Double I_ERI_Px_S_Pz_S_c = 0.0E0;
  Double I_ERI_Py_S_Pz_S_c = 0.0E0;
  Double I_ERI_Pz_S_Pz_S_c = 0.0E0;
  Double I_ERI_D2x_S_S_S_b = 0.0E0;
  Double I_ERI_Dxy_S_S_S_b = 0.0E0;
  Double I_ERI_Dxz_S_S_S_b = 0.0E0;
  Double I_ERI_D2y_S_S_S_b = 0.0E0;
  Double I_ERI_Dyz_S_S_S_b = 0.0E0;
  Double I_ERI_D2z_S_S_S_b = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_cc = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_cc = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_cc = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_cc = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_cc = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_cc = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_cc = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_cc = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_cc = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_cc = 0.0E0;
  Double I_ERI_D2x_S_S_S_c = 0.0E0;
  Double I_ERI_Dxy_S_S_S_c = 0.0E0;
  Double I_ERI_Dxz_S_S_S_c = 0.0E0;
  Double I_ERI_D2y_S_S_S_c = 0.0E0;
  Double I_ERI_Dyz_S_S_S_c = 0.0E0;
  Double I_ERI_D2z_S_S_S_c = 0.0E0;
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
  Double I_ERI_F3x_S_S_S_ab = 0.0E0;
  Double I_ERI_F2xy_S_S_S_ab = 0.0E0;
  Double I_ERI_F2xz_S_S_S_ab = 0.0E0;
  Double I_ERI_Fx2y_S_S_S_ab = 0.0E0;
  Double I_ERI_Fxyz_S_S_S_ab = 0.0E0;
  Double I_ERI_Fx2z_S_S_S_ab = 0.0E0;
  Double I_ERI_F3y_S_S_S_ab = 0.0E0;
  Double I_ERI_F2yz_S_S_S_ab = 0.0E0;
  Double I_ERI_Fy2z_S_S_S_ab = 0.0E0;
  Double I_ERI_F3z_S_S_S_ab = 0.0E0;
  Double I_ERI_Px_S_S_S_b = 0.0E0;
  Double I_ERI_Py_S_S_S_b = 0.0E0;
  Double I_ERI_Pz_S_S_S_b = 0.0E0;
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
  Double I_ERI_D2x_S_Px_S_bc = 0.0E0;
  Double I_ERI_Dxy_S_Px_S_bc = 0.0E0;
  Double I_ERI_Dxz_S_Px_S_bc = 0.0E0;
  Double I_ERI_D2y_S_Px_S_bc = 0.0E0;
  Double I_ERI_Dyz_S_Px_S_bc = 0.0E0;
  Double I_ERI_D2z_S_Px_S_bc = 0.0E0;
  Double I_ERI_D2x_S_Py_S_bc = 0.0E0;
  Double I_ERI_Dxy_S_Py_S_bc = 0.0E0;
  Double I_ERI_Dxz_S_Py_S_bc = 0.0E0;
  Double I_ERI_D2y_S_Py_S_bc = 0.0E0;
  Double I_ERI_Dyz_S_Py_S_bc = 0.0E0;
  Double I_ERI_D2z_S_Py_S_bc = 0.0E0;
  Double I_ERI_D2x_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Dxy_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Dxz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_D2y_S_Pz_S_bc = 0.0E0;
  Double I_ERI_Dyz_S_Pz_S_bc = 0.0E0;
  Double I_ERI_D2z_S_Pz_S_bc = 0.0E0;
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
  Double I_ERI_D2x_S_S_S_bb = 0.0E0;
  Double I_ERI_Dxy_S_S_S_bb = 0.0E0;
  Double I_ERI_Dxz_S_S_S_bb = 0.0E0;
  Double I_ERI_D2y_S_S_S_bb = 0.0E0;
  Double I_ERI_Dyz_S_S_S_bb = 0.0E0;
  Double I_ERI_D2z_S_S_S_bb = 0.0E0;

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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M1_vrr = PAX*I_ERI_D2x_S_S_S_M1_vrr+WPX*I_ERI_D2x_S_S_S_M2_vrr+2*oned2z*I_ERI_Px_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_F2xy_S_S_S_M1_vrr = PAY*I_ERI_D2x_S_S_S_M1_vrr+WPY*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_F2xz_S_S_S_M1_vrr = PAZ*I_ERI_D2x_S_S_S_M1_vrr+WPZ*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Fx2y_S_S_S_M1_vrr = PAX*I_ERI_D2y_S_S_S_M1_vrr+WPX*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Fx2z_S_S_S_M1_vrr = PAX*I_ERI_D2z_S_S_S_M1_vrr+WPX*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_F3y_S_S_S_M1_vrr = PAY*I_ERI_D2y_S_S_S_M1_vrr+WPY*I_ERI_D2y_S_S_S_M2_vrr+2*oned2z*I_ERI_Py_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_F2yz_S_S_S_M1_vrr = PAZ*I_ERI_D2y_S_S_S_M1_vrr+WPZ*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_F3z_S_S_S_M1_vrr = PAZ*I_ERI_D2z_S_S_S_M1_vrr+WPZ*I_ERI_D2z_S_S_S_M2_vrr+2*oned2z*I_ERI_Pz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_G_S_S_S_aa
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_S_S_aa_coefs = alpha*alpha;
      I_ERI_G4x_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S_aa += SQ_ERI_G_S_S_S_aa_coefs*I_ERI_G4z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_a_coefs = alpha;
      I_ERI_D2x_S_S_S_a += SQ_ERI_D_S_S_S_a_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_a += SQ_ERI_D_S_S_S_a_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_a += SQ_ERI_D_S_S_S_a_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_a += SQ_ERI_D_S_S_S_a_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_a += SQ_ERI_D_S_S_S_a_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_a += SQ_ERI_D_S_S_S_a_coefs*I_ERI_D2z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_S_S_S_S += I_ERI_S_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_ac
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_P_S_ac_coefs = alpha*gamma;
      I_ERI_F3x_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S_ac += SQ_ERI_F_S_P_S_ac_coefs*I_ERI_F3z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_P_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_P_S_c_coefs = gamma;
      I_ERI_Px_S_Px_S_c += SQ_ERI_P_S_P_S_c_coefs*I_ERI_Px_S_Px_S_vrr;
      I_ERI_Py_S_Px_S_c += SQ_ERI_P_S_P_S_c_coefs*I_ERI_Py_S_Px_S_vrr;
      I_ERI_Pz_S_Px_S_c += SQ_ERI_P_S_P_S_c_coefs*I_ERI_Pz_S_Px_S_vrr;
      I_ERI_Px_S_Py_S_c += SQ_ERI_P_S_P_S_c_coefs*I_ERI_Px_S_Py_S_vrr;
      I_ERI_Py_S_Py_S_c += SQ_ERI_P_S_P_S_c_coefs*I_ERI_Py_S_Py_S_vrr;
      I_ERI_Pz_S_Py_S_c += SQ_ERI_P_S_P_S_c_coefs*I_ERI_Pz_S_Py_S_vrr;
      I_ERI_Px_S_Pz_S_c += SQ_ERI_P_S_P_S_c_coefs*I_ERI_Px_S_Pz_S_vrr;
      I_ERI_Py_S_Pz_S_c += SQ_ERI_P_S_P_S_c_coefs*I_ERI_Py_S_Pz_S_vrr;
      I_ERI_Pz_S_Pz_S_c += SQ_ERI_P_S_P_S_c_coefs*I_ERI_Pz_S_Pz_S_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S_cc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_cc_coefs = gamma*gamma;
      I_ERI_D2x_S_D2x_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_cc += SQ_ERI_D_S_D_S_cc_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_c_coefs = gamma;
      I_ERI_D2x_S_S_S_c += SQ_ERI_D_S_S_S_c_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_c += SQ_ERI_D_S_S_S_c_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_c += SQ_ERI_D_S_S_S_c_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_c += SQ_ERI_D_S_S_S_c_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_c += SQ_ERI_D_S_S_S_c_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_c += SQ_ERI_D_S_S_S_c_coefs*I_ERI_D2z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_F_S_S_S_ab
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_S_S_ab_coefs = alpha*beta;
      I_ERI_F3x_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_F3x_S_S_S_vrr;
      I_ERI_F2xy_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_F2xy_S_S_S_vrr;
      I_ERI_F2xz_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_F2xz_S_S_S_vrr;
      I_ERI_Fx2y_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_Fx2y_S_S_S_vrr;
      I_ERI_Fxyz_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_Fxyz_S_S_S_vrr;
      I_ERI_Fx2z_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_Fx2z_S_S_S_vrr;
      I_ERI_F3y_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_F3y_S_S_S_vrr;
      I_ERI_F2yz_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_F2yz_S_S_S_vrr;
      I_ERI_Fy2z_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_Fy2z_S_S_S_vrr;
      I_ERI_F3z_S_S_S_ab += SQ_ERI_F_S_S_S_ab_coefs*I_ERI_F3z_S_S_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_P_S_S_S_b_coefs = beta;
      I_ERI_Px_S_S_S_b += SQ_ERI_P_S_S_S_b_coefs*I_ERI_Px_S_S_S_vrr;
      I_ERI_Py_S_S_S_b += SQ_ERI_P_S_S_S_b_coefs*I_ERI_Py_S_S_S_vrr;
      I_ERI_Pz_S_S_S_b += SQ_ERI_P_S_S_S_b_coefs*I_ERI_Pz_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_D_S_P_S_bc
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_P_S_bc_coefs = beta*gamma;
      I_ERI_D2x_S_Px_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_D2x_S_Px_S_vrr;
      I_ERI_Dxy_S_Px_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_Dxy_S_Px_S_vrr;
      I_ERI_Dxz_S_Px_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_Dxz_S_Px_S_vrr;
      I_ERI_D2y_S_Px_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_D2y_S_Px_S_vrr;
      I_ERI_Dyz_S_Px_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_Dyz_S_Px_S_vrr;
      I_ERI_D2z_S_Px_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_D2z_S_Px_S_vrr;
      I_ERI_D2x_S_Py_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_D2x_S_Py_S_vrr;
      I_ERI_Dxy_S_Py_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_Dxy_S_Py_S_vrr;
      I_ERI_Dxz_S_Py_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_Dxz_S_Py_S_vrr;
      I_ERI_D2y_S_Py_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_D2y_S_Py_S_vrr;
      I_ERI_Dyz_S_Py_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_Dyz_S_Py_S_vrr;
      I_ERI_D2z_S_Py_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_D2z_S_Py_S_vrr;
      I_ERI_D2x_S_Pz_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_D2x_S_Pz_S_vrr;
      I_ERI_Dxy_S_Pz_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_Dxy_S_Pz_S_vrr;
      I_ERI_Dxz_S_Pz_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_Dxz_S_Pz_S_vrr;
      I_ERI_D2y_S_Pz_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_D2y_S_Pz_S_vrr;
      I_ERI_Dyz_S_Pz_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_Dyz_S_Pz_S_vrr;
      I_ERI_D2z_S_Pz_S_bc += SQ_ERI_D_S_P_S_bc_coefs*I_ERI_D2z_S_Pz_S_vrr;

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

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_bb
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_S_S_bb_coefs = beta*beta;
      I_ERI_D2x_S_S_S_bb += SQ_ERI_D_S_S_S_bb_coefs*I_ERI_D2x_S_S_S_vrr;
      I_ERI_Dxy_S_S_S_bb += SQ_ERI_D_S_S_S_bb_coefs*I_ERI_Dxy_S_S_S_vrr;
      I_ERI_Dxz_S_S_S_bb += SQ_ERI_D_S_S_S_bb_coefs*I_ERI_Dxz_S_S_S_vrr;
      I_ERI_D2y_S_S_S_bb += SQ_ERI_D_S_S_S_bb_coefs*I_ERI_D2y_S_S_S_vrr;
      I_ERI_Dyz_S_S_S_bb += SQ_ERI_D_S_S_S_bb_coefs*I_ERI_Dyz_S_S_S_vrr;
      I_ERI_D2z_S_S_S_bb += SQ_ERI_D_S_S_S_bb_coefs*I_ERI_D2z_S_S_S_vrr;
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
   * shell quartet name: SQ_ERI_P_P_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_S_S_b
   * RHS shell quartet name: SQ_ERI_P_S_S_S_b
   ************************************************************/
  Double I_ERI_Px_Px_S_S_b = I_ERI_D2x_S_S_S_b+ABX*I_ERI_Px_S_S_S_b;
  Double I_ERI_Py_Px_S_S_b = I_ERI_Dxy_S_S_S_b+ABX*I_ERI_Py_S_S_S_b;
  Double I_ERI_Pz_Px_S_S_b = I_ERI_Dxz_S_S_S_b+ABX*I_ERI_Pz_S_S_S_b;
  Double I_ERI_Px_Py_S_S_b = I_ERI_Dxy_S_S_S_b+ABY*I_ERI_Px_S_S_S_b;
  Double I_ERI_Py_Py_S_S_b = I_ERI_D2y_S_S_S_b+ABY*I_ERI_Py_S_S_S_b;
  Double I_ERI_Pz_Py_S_S_b = I_ERI_Dyz_S_S_S_b+ABY*I_ERI_Pz_S_S_S_b;
  Double I_ERI_Px_Pz_S_S_b = I_ERI_Dxz_S_S_S_b+ABZ*I_ERI_Px_S_S_S_b;
  Double I_ERI_Py_Pz_S_S_b = I_ERI_Dyz_S_S_S_b+ABZ*I_ERI_Py_S_S_S_b;
  Double I_ERI_Pz_Pz_S_S_b = I_ERI_D2z_S_S_S_b+ABZ*I_ERI_Pz_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_ab
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_ab
   * RHS shell quartet name: SQ_ERI_F_S_S_S_ab
   ************************************************************/
  Double I_ERI_F3x_Px_S_S_ab = I_ERI_G4x_S_S_S_ab+ABX*I_ERI_F3x_S_S_S_ab;
  Double I_ERI_F2xy_Px_S_S_ab = I_ERI_G3xy_S_S_S_ab+ABX*I_ERI_F2xy_S_S_S_ab;
  Double I_ERI_F2xz_Px_S_S_ab = I_ERI_G3xz_S_S_S_ab+ABX*I_ERI_F2xz_S_S_S_ab;
  Double I_ERI_Fx2y_Px_S_S_ab = I_ERI_G2x2y_S_S_S_ab+ABX*I_ERI_Fx2y_S_S_S_ab;
  Double I_ERI_Fxyz_Px_S_S_ab = I_ERI_G2xyz_S_S_S_ab+ABX*I_ERI_Fxyz_S_S_S_ab;
  Double I_ERI_Fx2z_Px_S_S_ab = I_ERI_G2x2z_S_S_S_ab+ABX*I_ERI_Fx2z_S_S_S_ab;
  Double I_ERI_F3y_Px_S_S_ab = I_ERI_Gx3y_S_S_S_ab+ABX*I_ERI_F3y_S_S_S_ab;
  Double I_ERI_F2yz_Px_S_S_ab = I_ERI_Gx2yz_S_S_S_ab+ABX*I_ERI_F2yz_S_S_S_ab;
  Double I_ERI_Fy2z_Px_S_S_ab = I_ERI_Gxy2z_S_S_S_ab+ABX*I_ERI_Fy2z_S_S_S_ab;
  Double I_ERI_F3z_Px_S_S_ab = I_ERI_Gx3z_S_S_S_ab+ABX*I_ERI_F3z_S_S_S_ab;
  Double I_ERI_F3x_Py_S_S_ab = I_ERI_G3xy_S_S_S_ab+ABY*I_ERI_F3x_S_S_S_ab;
  Double I_ERI_F2xy_Py_S_S_ab = I_ERI_G2x2y_S_S_S_ab+ABY*I_ERI_F2xy_S_S_S_ab;
  Double I_ERI_F2xz_Py_S_S_ab = I_ERI_G2xyz_S_S_S_ab+ABY*I_ERI_F2xz_S_S_S_ab;
  Double I_ERI_Fx2y_Py_S_S_ab = I_ERI_Gx3y_S_S_S_ab+ABY*I_ERI_Fx2y_S_S_S_ab;
  Double I_ERI_Fxyz_Py_S_S_ab = I_ERI_Gx2yz_S_S_S_ab+ABY*I_ERI_Fxyz_S_S_S_ab;
  Double I_ERI_Fx2z_Py_S_S_ab = I_ERI_Gxy2z_S_S_S_ab+ABY*I_ERI_Fx2z_S_S_S_ab;
  Double I_ERI_F3y_Py_S_S_ab = I_ERI_G4y_S_S_S_ab+ABY*I_ERI_F3y_S_S_S_ab;
  Double I_ERI_F2yz_Py_S_S_ab = I_ERI_G3yz_S_S_S_ab+ABY*I_ERI_F2yz_S_S_S_ab;
  Double I_ERI_Fy2z_Py_S_S_ab = I_ERI_G2y2z_S_S_S_ab+ABY*I_ERI_Fy2z_S_S_S_ab;
  Double I_ERI_F3z_Py_S_S_ab = I_ERI_Gy3z_S_S_S_ab+ABY*I_ERI_F3z_S_S_S_ab;
  Double I_ERI_F3x_Pz_S_S_ab = I_ERI_G3xz_S_S_S_ab+ABZ*I_ERI_F3x_S_S_S_ab;
  Double I_ERI_F2xy_Pz_S_S_ab = I_ERI_G2xyz_S_S_S_ab+ABZ*I_ERI_F2xy_S_S_S_ab;
  Double I_ERI_F2xz_Pz_S_S_ab = I_ERI_G2x2z_S_S_S_ab+ABZ*I_ERI_F2xz_S_S_S_ab;
  Double I_ERI_Fx2y_Pz_S_S_ab = I_ERI_Gx2yz_S_S_S_ab+ABZ*I_ERI_Fx2y_S_S_S_ab;
  Double I_ERI_Fxyz_Pz_S_S_ab = I_ERI_Gxy2z_S_S_S_ab+ABZ*I_ERI_Fxyz_S_S_S_ab;
  Double I_ERI_Fx2z_Pz_S_S_ab = I_ERI_Gx3z_S_S_S_ab+ABZ*I_ERI_Fx2z_S_S_S_ab;
  Double I_ERI_F3y_Pz_S_S_ab = I_ERI_G3yz_S_S_S_ab+ABZ*I_ERI_F3y_S_S_S_ab;
  Double I_ERI_F2yz_Pz_S_S_ab = I_ERI_G2y2z_S_S_S_ab+ABZ*I_ERI_F2yz_S_S_S_ab;
  Double I_ERI_Fy2z_Pz_S_S_ab = I_ERI_Gy3z_S_S_S_ab+ABZ*I_ERI_Fy2z_S_S_S_ab;
  Double I_ERI_F3z_Pz_S_S_ab = I_ERI_G4z_S_S_S_ab+ABZ*I_ERI_F3z_S_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_S_S_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_bb
   ************************************************************/
  Double I_ERI_D2x_Px_S_S_bb = I_ERI_F3x_S_S_S_bb+ABX*I_ERI_D2x_S_S_S_bb;
  Double I_ERI_Dxy_Px_S_S_bb = I_ERI_F2xy_S_S_S_bb+ABX*I_ERI_Dxy_S_S_S_bb;
  Double I_ERI_Dxz_Px_S_S_bb = I_ERI_F2xz_S_S_S_bb+ABX*I_ERI_Dxz_S_S_S_bb;
  Double I_ERI_D2y_Px_S_S_bb = I_ERI_Fx2y_S_S_S_bb+ABX*I_ERI_D2y_S_S_S_bb;
  Double I_ERI_Dyz_Px_S_S_bb = I_ERI_Fxyz_S_S_S_bb+ABX*I_ERI_Dyz_S_S_S_bb;
  Double I_ERI_D2z_Px_S_S_bb = I_ERI_Fx2z_S_S_S_bb+ABX*I_ERI_D2z_S_S_S_bb;
  Double I_ERI_D2x_Py_S_S_bb = I_ERI_F2xy_S_S_S_bb+ABY*I_ERI_D2x_S_S_S_bb;
  Double I_ERI_Dxy_Py_S_S_bb = I_ERI_Fx2y_S_S_S_bb+ABY*I_ERI_Dxy_S_S_S_bb;
  Double I_ERI_Dxz_Py_S_S_bb = I_ERI_Fxyz_S_S_S_bb+ABY*I_ERI_Dxz_S_S_S_bb;
  Double I_ERI_D2y_Py_S_S_bb = I_ERI_F3y_S_S_S_bb+ABY*I_ERI_D2y_S_S_S_bb;
  Double I_ERI_Dyz_Py_S_S_bb = I_ERI_F2yz_S_S_S_bb+ABY*I_ERI_Dyz_S_S_S_bb;
  Double I_ERI_D2z_Py_S_S_bb = I_ERI_Fy2z_S_S_S_bb+ABY*I_ERI_D2z_S_S_S_bb;
  Double I_ERI_D2x_Pz_S_S_bb = I_ERI_F2xz_S_S_S_bb+ABZ*I_ERI_D2x_S_S_S_bb;
  Double I_ERI_Dxy_Pz_S_S_bb = I_ERI_Fxyz_S_S_S_bb+ABZ*I_ERI_Dxy_S_S_S_bb;
  Double I_ERI_Dxz_Pz_S_S_bb = I_ERI_Fx2z_S_S_S_bb+ABZ*I_ERI_Dxz_S_S_S_bb;
  Double I_ERI_D2y_Pz_S_S_bb = I_ERI_F2yz_S_S_S_bb+ABZ*I_ERI_D2y_S_S_S_bb;
  Double I_ERI_Dyz_Pz_S_S_bb = I_ERI_Fy2z_S_S_S_bb+ABZ*I_ERI_Dyz_S_S_S_bb;
  Double I_ERI_D2z_Pz_S_S_bb = I_ERI_F3z_S_S_S_bb+ABZ*I_ERI_D2z_S_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 5 integrals are omitted 
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
  Double I_ERI_F2xy_Py_S_S_bb = I_ERI_G2x2y_S_S_S_bb+ABY*I_ERI_F2xy_S_S_S_bb;
  Double I_ERI_F2xz_Py_S_S_bb = I_ERI_G2xyz_S_S_S_bb+ABY*I_ERI_F2xz_S_S_S_bb;
  Double I_ERI_Fx2y_Py_S_S_bb = I_ERI_Gx3y_S_S_S_bb+ABY*I_ERI_Fx2y_S_S_S_bb;
  Double I_ERI_Fxyz_Py_S_S_bb = I_ERI_Gx2yz_S_S_S_bb+ABY*I_ERI_Fxyz_S_S_S_bb;
  Double I_ERI_Fx2z_Py_S_S_bb = I_ERI_Gxy2z_S_S_S_bb+ABY*I_ERI_Fx2z_S_S_S_bb;
  Double I_ERI_F3y_Py_S_S_bb = I_ERI_G4y_S_S_S_bb+ABY*I_ERI_F3y_S_S_S_bb;
  Double I_ERI_F2yz_Py_S_S_bb = I_ERI_G3yz_S_S_S_bb+ABY*I_ERI_F2yz_S_S_S_bb;
  Double I_ERI_Fy2z_Py_S_S_bb = I_ERI_G2y2z_S_S_S_bb+ABY*I_ERI_Fy2z_S_S_S_bb;
  Double I_ERI_F3z_Py_S_S_bb = I_ERI_Gy3z_S_S_S_bb+ABY*I_ERI_F3z_S_S_S_bb;
  Double I_ERI_F2xz_Pz_S_S_bb = I_ERI_G2x2z_S_S_S_bb+ABZ*I_ERI_F2xz_S_S_S_bb;
  Double I_ERI_Fxyz_Pz_S_S_bb = I_ERI_Gxy2z_S_S_S_bb+ABZ*I_ERI_Fxyz_S_S_S_bb;
  Double I_ERI_Fx2z_Pz_S_S_bb = I_ERI_Gx3z_S_S_S_bb+ABZ*I_ERI_Fx2z_S_S_S_bb;
  Double I_ERI_F2yz_Pz_S_S_bb = I_ERI_G2y2z_S_S_S_bb+ABZ*I_ERI_F2yz_S_S_S_bb;
  Double I_ERI_Fy2z_Pz_S_S_bb = I_ERI_Gy3z_S_S_S_bb+ABZ*I_ERI_Fy2z_S_S_S_bb;
  Double I_ERI_F3z_Pz_S_S_bb = I_ERI_G4z_S_S_S_bb+ABZ*I_ERI_F3z_S_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_S_S_bb
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_bb
   * RHS shell quartet name: SQ_ERI_D_P_S_S_bb
   ************************************************************/
  Double I_ERI_D2x_D2x_S_S_bb = I_ERI_F3x_Px_S_S_bb+ABX*I_ERI_D2x_Px_S_S_bb;
  Double I_ERI_Dxy_D2x_S_S_bb = I_ERI_F2xy_Px_S_S_bb+ABX*I_ERI_Dxy_Px_S_S_bb;
  Double I_ERI_Dxz_D2x_S_S_bb = I_ERI_F2xz_Px_S_S_bb+ABX*I_ERI_Dxz_Px_S_S_bb;
  Double I_ERI_D2y_D2x_S_S_bb = I_ERI_Fx2y_Px_S_S_bb+ABX*I_ERI_D2y_Px_S_S_bb;
  Double I_ERI_Dyz_D2x_S_S_bb = I_ERI_Fxyz_Px_S_S_bb+ABX*I_ERI_Dyz_Px_S_S_bb;
  Double I_ERI_D2z_D2x_S_S_bb = I_ERI_Fx2z_Px_S_S_bb+ABX*I_ERI_D2z_Px_S_S_bb;
  Double I_ERI_D2x_Dxy_S_S_bb = I_ERI_F2xy_Px_S_S_bb+ABY*I_ERI_D2x_Px_S_S_bb;
  Double I_ERI_Dxy_Dxy_S_S_bb = I_ERI_Fx2y_Px_S_S_bb+ABY*I_ERI_Dxy_Px_S_S_bb;
  Double I_ERI_Dxz_Dxy_S_S_bb = I_ERI_Fxyz_Px_S_S_bb+ABY*I_ERI_Dxz_Px_S_S_bb;
  Double I_ERI_D2y_Dxy_S_S_bb = I_ERI_F3y_Px_S_S_bb+ABY*I_ERI_D2y_Px_S_S_bb;
  Double I_ERI_Dyz_Dxy_S_S_bb = I_ERI_F2yz_Px_S_S_bb+ABY*I_ERI_Dyz_Px_S_S_bb;
  Double I_ERI_D2z_Dxy_S_S_bb = I_ERI_Fy2z_Px_S_S_bb+ABY*I_ERI_D2z_Px_S_S_bb;
  Double I_ERI_D2x_Dxz_S_S_bb = I_ERI_F2xz_Px_S_S_bb+ABZ*I_ERI_D2x_Px_S_S_bb;
  Double I_ERI_Dxy_Dxz_S_S_bb = I_ERI_Fxyz_Px_S_S_bb+ABZ*I_ERI_Dxy_Px_S_S_bb;
  Double I_ERI_Dxz_Dxz_S_S_bb = I_ERI_Fx2z_Px_S_S_bb+ABZ*I_ERI_Dxz_Px_S_S_bb;
  Double I_ERI_D2y_Dxz_S_S_bb = I_ERI_F2yz_Px_S_S_bb+ABZ*I_ERI_D2y_Px_S_S_bb;
  Double I_ERI_Dyz_Dxz_S_S_bb = I_ERI_Fy2z_Px_S_S_bb+ABZ*I_ERI_Dyz_Px_S_S_bb;
  Double I_ERI_D2z_Dxz_S_S_bb = I_ERI_F3z_Px_S_S_bb+ABZ*I_ERI_D2z_Px_S_S_bb;
  Double I_ERI_D2x_D2y_S_S_bb = I_ERI_F2xy_Py_S_S_bb+ABY*I_ERI_D2x_Py_S_S_bb;
  Double I_ERI_Dxy_D2y_S_S_bb = I_ERI_Fx2y_Py_S_S_bb+ABY*I_ERI_Dxy_Py_S_S_bb;
  Double I_ERI_Dxz_D2y_S_S_bb = I_ERI_Fxyz_Py_S_S_bb+ABY*I_ERI_Dxz_Py_S_S_bb;
  Double I_ERI_D2y_D2y_S_S_bb = I_ERI_F3y_Py_S_S_bb+ABY*I_ERI_D2y_Py_S_S_bb;
  Double I_ERI_Dyz_D2y_S_S_bb = I_ERI_F2yz_Py_S_S_bb+ABY*I_ERI_Dyz_Py_S_S_bb;
  Double I_ERI_D2z_D2y_S_S_bb = I_ERI_Fy2z_Py_S_S_bb+ABY*I_ERI_D2z_Py_S_S_bb;
  Double I_ERI_D2x_Dyz_S_S_bb = I_ERI_F2xz_Py_S_S_bb+ABZ*I_ERI_D2x_Py_S_S_bb;
  Double I_ERI_Dxy_Dyz_S_S_bb = I_ERI_Fxyz_Py_S_S_bb+ABZ*I_ERI_Dxy_Py_S_S_bb;
  Double I_ERI_Dxz_Dyz_S_S_bb = I_ERI_Fx2z_Py_S_S_bb+ABZ*I_ERI_Dxz_Py_S_S_bb;
  Double I_ERI_D2y_Dyz_S_S_bb = I_ERI_F2yz_Py_S_S_bb+ABZ*I_ERI_D2y_Py_S_S_bb;
  Double I_ERI_Dyz_Dyz_S_S_bb = I_ERI_Fy2z_Py_S_S_bb+ABZ*I_ERI_Dyz_Py_S_S_bb;
  Double I_ERI_D2z_Dyz_S_S_bb = I_ERI_F3z_Py_S_S_bb+ABZ*I_ERI_D2z_Py_S_S_bb;
  Double I_ERI_D2x_D2z_S_S_bb = I_ERI_F2xz_Pz_S_S_bb+ABZ*I_ERI_D2x_Pz_S_S_bb;
  Double I_ERI_Dxy_D2z_S_S_bb = I_ERI_Fxyz_Pz_S_S_bb+ABZ*I_ERI_Dxy_Pz_S_S_bb;
  Double I_ERI_Dxz_D2z_S_S_bb = I_ERI_Fx2z_Pz_S_S_bb+ABZ*I_ERI_Dxz_Pz_S_S_bb;
  Double I_ERI_D2y_D2z_S_S_bb = I_ERI_F2yz_Pz_S_S_bb+ABZ*I_ERI_D2y_Pz_S_S_bb;
  Double I_ERI_Dyz_D2z_S_S_bb = I_ERI_Fy2z_Pz_S_S_bb+ABZ*I_ERI_Dyz_Pz_S_S_bb;
  Double I_ERI_D2z_D2z_S_S_bb = I_ERI_F3z_Pz_S_S_bb+ABZ*I_ERI_D2z_Pz_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_P_S_bc
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_bc
   * RHS shell quartet name: SQ_ERI_D_S_P_S_bc
   ************************************************************/
  Double I_ERI_D2x_Px_Px_S_bc = I_ERI_F3x_S_Px_S_bc+ABX*I_ERI_D2x_S_Px_S_bc;
  Double I_ERI_Dxy_Px_Px_S_bc = I_ERI_F2xy_S_Px_S_bc+ABX*I_ERI_Dxy_S_Px_S_bc;
  Double I_ERI_Dxz_Px_Px_S_bc = I_ERI_F2xz_S_Px_S_bc+ABX*I_ERI_Dxz_S_Px_S_bc;
  Double I_ERI_D2y_Px_Px_S_bc = I_ERI_Fx2y_S_Px_S_bc+ABX*I_ERI_D2y_S_Px_S_bc;
  Double I_ERI_Dyz_Px_Px_S_bc = I_ERI_Fxyz_S_Px_S_bc+ABX*I_ERI_Dyz_S_Px_S_bc;
  Double I_ERI_D2z_Px_Px_S_bc = I_ERI_Fx2z_S_Px_S_bc+ABX*I_ERI_D2z_S_Px_S_bc;
  Double I_ERI_D2x_Py_Px_S_bc = I_ERI_F2xy_S_Px_S_bc+ABY*I_ERI_D2x_S_Px_S_bc;
  Double I_ERI_Dxy_Py_Px_S_bc = I_ERI_Fx2y_S_Px_S_bc+ABY*I_ERI_Dxy_S_Px_S_bc;
  Double I_ERI_Dxz_Py_Px_S_bc = I_ERI_Fxyz_S_Px_S_bc+ABY*I_ERI_Dxz_S_Px_S_bc;
  Double I_ERI_D2y_Py_Px_S_bc = I_ERI_F3y_S_Px_S_bc+ABY*I_ERI_D2y_S_Px_S_bc;
  Double I_ERI_Dyz_Py_Px_S_bc = I_ERI_F2yz_S_Px_S_bc+ABY*I_ERI_Dyz_S_Px_S_bc;
  Double I_ERI_D2z_Py_Px_S_bc = I_ERI_Fy2z_S_Px_S_bc+ABY*I_ERI_D2z_S_Px_S_bc;
  Double I_ERI_D2x_Pz_Px_S_bc = I_ERI_F2xz_S_Px_S_bc+ABZ*I_ERI_D2x_S_Px_S_bc;
  Double I_ERI_Dxy_Pz_Px_S_bc = I_ERI_Fxyz_S_Px_S_bc+ABZ*I_ERI_Dxy_S_Px_S_bc;
  Double I_ERI_Dxz_Pz_Px_S_bc = I_ERI_Fx2z_S_Px_S_bc+ABZ*I_ERI_Dxz_S_Px_S_bc;
  Double I_ERI_D2y_Pz_Px_S_bc = I_ERI_F2yz_S_Px_S_bc+ABZ*I_ERI_D2y_S_Px_S_bc;
  Double I_ERI_Dyz_Pz_Px_S_bc = I_ERI_Fy2z_S_Px_S_bc+ABZ*I_ERI_Dyz_S_Px_S_bc;
  Double I_ERI_D2z_Pz_Px_S_bc = I_ERI_F3z_S_Px_S_bc+ABZ*I_ERI_D2z_S_Px_S_bc;
  Double I_ERI_D2x_Px_Py_S_bc = I_ERI_F3x_S_Py_S_bc+ABX*I_ERI_D2x_S_Py_S_bc;
  Double I_ERI_Dxy_Px_Py_S_bc = I_ERI_F2xy_S_Py_S_bc+ABX*I_ERI_Dxy_S_Py_S_bc;
  Double I_ERI_Dxz_Px_Py_S_bc = I_ERI_F2xz_S_Py_S_bc+ABX*I_ERI_Dxz_S_Py_S_bc;
  Double I_ERI_D2y_Px_Py_S_bc = I_ERI_Fx2y_S_Py_S_bc+ABX*I_ERI_D2y_S_Py_S_bc;
  Double I_ERI_Dyz_Px_Py_S_bc = I_ERI_Fxyz_S_Py_S_bc+ABX*I_ERI_Dyz_S_Py_S_bc;
  Double I_ERI_D2z_Px_Py_S_bc = I_ERI_Fx2z_S_Py_S_bc+ABX*I_ERI_D2z_S_Py_S_bc;
  Double I_ERI_D2x_Py_Py_S_bc = I_ERI_F2xy_S_Py_S_bc+ABY*I_ERI_D2x_S_Py_S_bc;
  Double I_ERI_Dxy_Py_Py_S_bc = I_ERI_Fx2y_S_Py_S_bc+ABY*I_ERI_Dxy_S_Py_S_bc;
  Double I_ERI_Dxz_Py_Py_S_bc = I_ERI_Fxyz_S_Py_S_bc+ABY*I_ERI_Dxz_S_Py_S_bc;
  Double I_ERI_D2y_Py_Py_S_bc = I_ERI_F3y_S_Py_S_bc+ABY*I_ERI_D2y_S_Py_S_bc;
  Double I_ERI_Dyz_Py_Py_S_bc = I_ERI_F2yz_S_Py_S_bc+ABY*I_ERI_Dyz_S_Py_S_bc;
  Double I_ERI_D2z_Py_Py_S_bc = I_ERI_Fy2z_S_Py_S_bc+ABY*I_ERI_D2z_S_Py_S_bc;
  Double I_ERI_D2x_Pz_Py_S_bc = I_ERI_F2xz_S_Py_S_bc+ABZ*I_ERI_D2x_S_Py_S_bc;
  Double I_ERI_Dxy_Pz_Py_S_bc = I_ERI_Fxyz_S_Py_S_bc+ABZ*I_ERI_Dxy_S_Py_S_bc;
  Double I_ERI_Dxz_Pz_Py_S_bc = I_ERI_Fx2z_S_Py_S_bc+ABZ*I_ERI_Dxz_S_Py_S_bc;
  Double I_ERI_D2y_Pz_Py_S_bc = I_ERI_F2yz_S_Py_S_bc+ABZ*I_ERI_D2y_S_Py_S_bc;
  Double I_ERI_Dyz_Pz_Py_S_bc = I_ERI_Fy2z_S_Py_S_bc+ABZ*I_ERI_Dyz_S_Py_S_bc;
  Double I_ERI_D2z_Pz_Py_S_bc = I_ERI_F3z_S_Py_S_bc+ABZ*I_ERI_D2z_S_Py_S_bc;
  Double I_ERI_D2x_Px_Pz_S_bc = I_ERI_F3x_S_Pz_S_bc+ABX*I_ERI_D2x_S_Pz_S_bc;
  Double I_ERI_Dxy_Px_Pz_S_bc = I_ERI_F2xy_S_Pz_S_bc+ABX*I_ERI_Dxy_S_Pz_S_bc;
  Double I_ERI_Dxz_Px_Pz_S_bc = I_ERI_F2xz_S_Pz_S_bc+ABX*I_ERI_Dxz_S_Pz_S_bc;
  Double I_ERI_D2y_Px_Pz_S_bc = I_ERI_Fx2y_S_Pz_S_bc+ABX*I_ERI_D2y_S_Pz_S_bc;
  Double I_ERI_Dyz_Px_Pz_S_bc = I_ERI_Fxyz_S_Pz_S_bc+ABX*I_ERI_Dyz_S_Pz_S_bc;
  Double I_ERI_D2z_Px_Pz_S_bc = I_ERI_Fx2z_S_Pz_S_bc+ABX*I_ERI_D2z_S_Pz_S_bc;
  Double I_ERI_D2x_Py_Pz_S_bc = I_ERI_F2xy_S_Pz_S_bc+ABY*I_ERI_D2x_S_Pz_S_bc;
  Double I_ERI_Dxy_Py_Pz_S_bc = I_ERI_Fx2y_S_Pz_S_bc+ABY*I_ERI_Dxy_S_Pz_S_bc;
  Double I_ERI_Dxz_Py_Pz_S_bc = I_ERI_Fxyz_S_Pz_S_bc+ABY*I_ERI_Dxz_S_Pz_S_bc;
  Double I_ERI_D2y_Py_Pz_S_bc = I_ERI_F3y_S_Pz_S_bc+ABY*I_ERI_D2y_S_Pz_S_bc;
  Double I_ERI_Dyz_Py_Pz_S_bc = I_ERI_F2yz_S_Pz_S_bc+ABY*I_ERI_Dyz_S_Pz_S_bc;
  Double I_ERI_D2z_Py_Pz_S_bc = I_ERI_Fy2z_S_Pz_S_bc+ABY*I_ERI_D2z_S_Pz_S_bc;
  Double I_ERI_D2x_Pz_Pz_S_bc = I_ERI_F2xz_S_Pz_S_bc+ABZ*I_ERI_D2x_S_Pz_S_bc;
  Double I_ERI_Dxy_Pz_Pz_S_bc = I_ERI_Fxyz_S_Pz_S_bc+ABZ*I_ERI_Dxy_S_Pz_S_bc;
  Double I_ERI_Dxz_Pz_Pz_S_bc = I_ERI_Fx2z_S_Pz_S_bc+ABZ*I_ERI_Dxz_S_Pz_S_bc;
  Double I_ERI_D2y_Pz_Pz_S_bc = I_ERI_F2yz_S_Pz_S_bc+ABZ*I_ERI_D2y_S_Pz_S_bc;
  Double I_ERI_Dyz_Pz_Pz_S_bc = I_ERI_Fy2z_S_Pz_S_bc+ABZ*I_ERI_Dyz_S_Pz_S_bc;
  Double I_ERI_D2z_Pz_Pz_S_bc = I_ERI_F3z_S_Pz_S_bc+ABZ*I_ERI_D2z_S_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S
   ************************************************************/
  abcd[0] = 4.0E0*I_ERI_G4x_S_S_S_aa-2.0E0*2*I_ERI_D2x_S_S_S_a-2.0E0*3*I_ERI_D2x_S_S_S_a+2*1*I_ERI_S_S_S_S;
  abcd[1] = 4.0E0*I_ERI_G3xy_S_S_S_aa-2.0E0*1*I_ERI_Dxy_S_S_S_a-2.0E0*2*I_ERI_Dxy_S_S_S_a;
  abcd[2] = 4.0E0*I_ERI_G3xz_S_S_S_aa-2.0E0*1*I_ERI_Dxz_S_S_S_a-2.0E0*2*I_ERI_Dxz_S_S_S_a;
  abcd[3] = 4.0E0*I_ERI_G2x2y_S_S_S_aa-2.0E0*1*I_ERI_D2y_S_S_S_a;
  abcd[4] = 4.0E0*I_ERI_G2xyz_S_S_S_aa-2.0E0*1*I_ERI_Dyz_S_S_S_a;
  abcd[5] = 4.0E0*I_ERI_G2x2z_S_S_S_aa-2.0E0*1*I_ERI_D2z_S_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S
   ************************************************************/
  abcd[6] = 4.0E0*I_ERI_G3xy_S_S_S_aa-2.0E0*2*I_ERI_Dxy_S_S_S_a;
  abcd[7] = 4.0E0*I_ERI_G2x2y_S_S_S_aa-2.0E0*1*I_ERI_D2x_S_S_S_a-2.0E0*1*I_ERI_D2y_S_S_S_a+1*I_ERI_S_S_S_S;
  abcd[8] = 4.0E0*I_ERI_G2xyz_S_S_S_aa-2.0E0*1*I_ERI_Dyz_S_S_S_a;
  abcd[9] = 4.0E0*I_ERI_Gx3y_S_S_S_aa-2.0E0*2*I_ERI_Dxy_S_S_S_a;
  abcd[10] = 4.0E0*I_ERI_Gx2yz_S_S_S_aa-2.0E0*1*I_ERI_Dxz_S_S_S_a;
  abcd[11] = 4.0E0*I_ERI_Gxy2z_S_S_S_aa;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S
   ************************************************************/
  abcd[12] = 4.0E0*I_ERI_G3xz_S_S_S_aa-2.0E0*2*I_ERI_Dxz_S_S_S_a;
  abcd[13] = 4.0E0*I_ERI_G2xyz_S_S_S_aa-2.0E0*1*I_ERI_Dyz_S_S_S_a;
  abcd[14] = 4.0E0*I_ERI_G2x2z_S_S_S_aa-2.0E0*1*I_ERI_D2x_S_S_S_a-2.0E0*1*I_ERI_D2z_S_S_S_a+1*I_ERI_S_S_S_S;
  abcd[15] = 4.0E0*I_ERI_Gx2yz_S_S_S_aa;
  abcd[16] = 4.0E0*I_ERI_Gxy2z_S_S_S_aa-2.0E0*1*I_ERI_Dxy_S_S_S_a;
  abcd[17] = 4.0E0*I_ERI_Gx3z_S_S_S_aa-2.0E0*2*I_ERI_Dxz_S_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S
   ************************************************************/
  abcd[18] = 4.0E0*I_ERI_G2x2y_S_S_S_aa-2.0E0*1*I_ERI_D2x_S_S_S_a;
  abcd[19] = 4.0E0*I_ERI_Gx3y_S_S_S_aa-2.0E0*1*I_ERI_Dxy_S_S_S_a-2.0E0*2*I_ERI_Dxy_S_S_S_a;
  abcd[20] = 4.0E0*I_ERI_Gx2yz_S_S_S_aa-2.0E0*1*I_ERI_Dxz_S_S_S_a;
  abcd[21] = 4.0E0*I_ERI_G4y_S_S_S_aa-2.0E0*2*I_ERI_D2y_S_S_S_a-2.0E0*3*I_ERI_D2y_S_S_S_a+2*1*I_ERI_S_S_S_S;
  abcd[22] = 4.0E0*I_ERI_G3yz_S_S_S_aa-2.0E0*1*I_ERI_Dyz_S_S_S_a-2.0E0*2*I_ERI_Dyz_S_S_S_a;
  abcd[23] = 4.0E0*I_ERI_G2y2z_S_S_S_aa-2.0E0*1*I_ERI_D2z_S_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S
   ************************************************************/
  abcd[24] = 4.0E0*I_ERI_G2xyz_S_S_S_aa;
  abcd[25] = 4.0E0*I_ERI_Gx2yz_S_S_S_aa-2.0E0*1*I_ERI_Dxz_S_S_S_a;
  abcd[26] = 4.0E0*I_ERI_Gxy2z_S_S_S_aa-2.0E0*1*I_ERI_Dxy_S_S_S_a;
  abcd[27] = 4.0E0*I_ERI_G3yz_S_S_S_aa-2.0E0*2*I_ERI_Dyz_S_S_S_a;
  abcd[28] = 4.0E0*I_ERI_G2y2z_S_S_S_aa-2.0E0*1*I_ERI_D2y_S_S_S_a-2.0E0*1*I_ERI_D2z_S_S_S_a+1*I_ERI_S_S_S_S;
  abcd[29] = 4.0E0*I_ERI_Gy3z_S_S_S_aa-2.0E0*2*I_ERI_Dyz_S_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S_aa
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_D_S_S_S_a
   * RHS shell quartet name: SQ_ERI_S_S_S_S
   ************************************************************/
  abcd[30] = 4.0E0*I_ERI_G2x2z_S_S_S_aa-2.0E0*1*I_ERI_D2x_S_S_S_a;
  abcd[31] = 4.0E0*I_ERI_Gxy2z_S_S_S_aa-2.0E0*1*I_ERI_Dxy_S_S_S_a;
  abcd[32] = 4.0E0*I_ERI_Gx3z_S_S_S_aa-2.0E0*1*I_ERI_Dxz_S_S_S_a-2.0E0*2*I_ERI_Dxz_S_S_S_a;
  abcd[33] = 4.0E0*I_ERI_G2y2z_S_S_S_aa-2.0E0*1*I_ERI_D2y_S_S_S_a;
  abcd[34] = 4.0E0*I_ERI_Gy3z_S_S_S_aa-2.0E0*1*I_ERI_Dyz_S_S_S_a-2.0E0*2*I_ERI_Dyz_S_S_S_a;
  abcd[35] = 4.0E0*I_ERI_G4z_S_S_S_aa-2.0E0*2*I_ERI_D2z_S_S_S_a-2.0E0*3*I_ERI_D2z_S_S_S_a+2*1*I_ERI_S_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dax_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_b
   ************************************************************/
  abcd[36] = 4.0E0*I_ERI_F3x_Px_S_S_ab-2.0E0*2*I_ERI_Px_Px_S_S_b;
  abcd[37] = 4.0E0*I_ERI_F2xy_Px_S_S_ab-2.0E0*1*I_ERI_Py_Px_S_S_b;
  abcd[38] = 4.0E0*I_ERI_F2xz_Px_S_S_ab-2.0E0*1*I_ERI_Pz_Px_S_S_b;
  abcd[39] = 4.0E0*I_ERI_Fx2y_Px_S_S_ab;
  abcd[40] = 4.0E0*I_ERI_Fxyz_Px_S_S_ab;
  abcd[41] = 4.0E0*I_ERI_Fx2z_Px_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dax_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_b
   ************************************************************/
  abcd[42] = 4.0E0*I_ERI_F3x_Py_S_S_ab-2.0E0*2*I_ERI_Px_Py_S_S_b;
  abcd[43] = 4.0E0*I_ERI_F2xy_Py_S_S_ab-2.0E0*1*I_ERI_Py_Py_S_S_b;
  abcd[44] = 4.0E0*I_ERI_F2xz_Py_S_S_ab-2.0E0*1*I_ERI_Pz_Py_S_S_b;
  abcd[45] = 4.0E0*I_ERI_Fx2y_Py_S_S_ab;
  abcd[46] = 4.0E0*I_ERI_Fxyz_Py_S_S_ab;
  abcd[47] = 4.0E0*I_ERI_Fx2z_Py_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dax_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_b
   ************************************************************/
  abcd[48] = 4.0E0*I_ERI_F3x_Pz_S_S_ab-2.0E0*2*I_ERI_Px_Pz_S_S_b;
  abcd[49] = 4.0E0*I_ERI_F2xy_Pz_S_S_ab-2.0E0*1*I_ERI_Py_Pz_S_S_b;
  abcd[50] = 4.0E0*I_ERI_F2xz_Pz_S_S_ab-2.0E0*1*I_ERI_Pz_Pz_S_S_b;
  abcd[51] = 4.0E0*I_ERI_Fx2y_Pz_S_S_ab;
  abcd[52] = 4.0E0*I_ERI_Fxyz_Pz_S_S_ab;
  abcd[53] = 4.0E0*I_ERI_Fx2z_Pz_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_day_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_b
   ************************************************************/
  abcd[54] = 4.0E0*I_ERI_F2xy_Px_S_S_ab;
  abcd[55] = 4.0E0*I_ERI_Fx2y_Px_S_S_ab-2.0E0*1*I_ERI_Px_Px_S_S_b;
  abcd[56] = 4.0E0*I_ERI_Fxyz_Px_S_S_ab;
  abcd[57] = 4.0E0*I_ERI_F3y_Px_S_S_ab-2.0E0*2*I_ERI_Py_Px_S_S_b;
  abcd[58] = 4.0E0*I_ERI_F2yz_Px_S_S_ab-2.0E0*1*I_ERI_Pz_Px_S_S_b;
  abcd[59] = 4.0E0*I_ERI_Fy2z_Px_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_day_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_b
   ************************************************************/
  abcd[60] = 4.0E0*I_ERI_F2xy_Py_S_S_ab;
  abcd[61] = 4.0E0*I_ERI_Fx2y_Py_S_S_ab-2.0E0*1*I_ERI_Px_Py_S_S_b;
  abcd[62] = 4.0E0*I_ERI_Fxyz_Py_S_S_ab;
  abcd[63] = 4.0E0*I_ERI_F3y_Py_S_S_ab-2.0E0*2*I_ERI_Py_Py_S_S_b;
  abcd[64] = 4.0E0*I_ERI_F2yz_Py_S_S_ab-2.0E0*1*I_ERI_Pz_Py_S_S_b;
  abcd[65] = 4.0E0*I_ERI_Fy2z_Py_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_day_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_b
   ************************************************************/
  abcd[66] = 4.0E0*I_ERI_F2xy_Pz_S_S_ab;
  abcd[67] = 4.0E0*I_ERI_Fx2y_Pz_S_S_ab-2.0E0*1*I_ERI_Px_Pz_S_S_b;
  abcd[68] = 4.0E0*I_ERI_Fxyz_Pz_S_S_ab;
  abcd[69] = 4.0E0*I_ERI_F3y_Pz_S_S_ab-2.0E0*2*I_ERI_Py_Pz_S_S_b;
  abcd[70] = 4.0E0*I_ERI_F2yz_Pz_S_S_ab-2.0E0*1*I_ERI_Pz_Pz_S_S_b;
  abcd[71] = 4.0E0*I_ERI_Fy2z_Pz_S_S_ab;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_daz_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_b
   ************************************************************/
  abcd[72] = 4.0E0*I_ERI_F2xz_Px_S_S_ab;
  abcd[73] = 4.0E0*I_ERI_Fxyz_Px_S_S_ab;
  abcd[74] = 4.0E0*I_ERI_Fx2z_Px_S_S_ab-2.0E0*1*I_ERI_Px_Px_S_S_b;
  abcd[75] = 4.0E0*I_ERI_F2yz_Px_S_S_ab;
  abcd[76] = 4.0E0*I_ERI_Fy2z_Px_S_S_ab-2.0E0*1*I_ERI_Py_Px_S_S_b;
  abcd[77] = 4.0E0*I_ERI_F3z_Px_S_S_ab-2.0E0*2*I_ERI_Pz_Px_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_daz_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_b
   ************************************************************/
  abcd[78] = 4.0E0*I_ERI_F2xz_Py_S_S_ab;
  abcd[79] = 4.0E0*I_ERI_Fxyz_Py_S_S_ab;
  abcd[80] = 4.0E0*I_ERI_Fx2z_Py_S_S_ab-2.0E0*1*I_ERI_Px_Py_S_S_b;
  abcd[81] = 4.0E0*I_ERI_F2yz_Py_S_S_ab;
  abcd[82] = 4.0E0*I_ERI_Fy2z_Py_S_S_ab-2.0E0*1*I_ERI_Py_Py_S_S_b;
  abcd[83] = 4.0E0*I_ERI_F3z_Py_S_S_ab-2.0E0*2*I_ERI_Pz_Py_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_daz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S_ab
   * RHS shell quartet name: SQ_ERI_P_P_S_S_b
   ************************************************************/
  abcd[84] = 4.0E0*I_ERI_F2xz_Pz_S_S_ab;
  abcd[85] = 4.0E0*I_ERI_Fxyz_Pz_S_S_ab;
  abcd[86] = 4.0E0*I_ERI_Fx2z_Pz_S_S_ab-2.0E0*1*I_ERI_Px_Pz_S_S_b;
  abcd[87] = 4.0E0*I_ERI_F2yz_Pz_S_S_ab;
  abcd[88] = 4.0E0*I_ERI_Fy2z_Pz_S_S_ab-2.0E0*1*I_ERI_Py_Pz_S_S_b;
  abcd[89] = 4.0E0*I_ERI_F3z_Pz_S_S_ab-2.0E0*2*I_ERI_Pz_Pz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dax_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_c
   ************************************************************/
  abcd[90] = 4.0E0*I_ERI_F3x_S_Px_S_ac-2.0E0*2*I_ERI_Px_S_Px_S_c;
  abcd[91] = 4.0E0*I_ERI_F2xy_S_Px_S_ac-2.0E0*1*I_ERI_Py_S_Px_S_c;
  abcd[92] = 4.0E0*I_ERI_F2xz_S_Px_S_ac-2.0E0*1*I_ERI_Pz_S_Px_S_c;
  abcd[93] = 4.0E0*I_ERI_Fx2y_S_Px_S_ac;
  abcd[94] = 4.0E0*I_ERI_Fxyz_S_Px_S_ac;
  abcd[95] = 4.0E0*I_ERI_Fx2z_S_Px_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dax_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_c
   ************************************************************/
  abcd[96] = 4.0E0*I_ERI_F3x_S_Py_S_ac-2.0E0*2*I_ERI_Px_S_Py_S_c;
  abcd[97] = 4.0E0*I_ERI_F2xy_S_Py_S_ac-2.0E0*1*I_ERI_Py_S_Py_S_c;
  abcd[98] = 4.0E0*I_ERI_F2xz_S_Py_S_ac-2.0E0*1*I_ERI_Pz_S_Py_S_c;
  abcd[99] = 4.0E0*I_ERI_Fx2y_S_Py_S_ac;
  abcd[100] = 4.0E0*I_ERI_Fxyz_S_Py_S_ac;
  abcd[101] = 4.0E0*I_ERI_Fx2z_S_Py_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dax_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_c
   ************************************************************/
  abcd[102] = 4.0E0*I_ERI_F3x_S_Pz_S_ac-2.0E0*2*I_ERI_Px_S_Pz_S_c;
  abcd[103] = 4.0E0*I_ERI_F2xy_S_Pz_S_ac-2.0E0*1*I_ERI_Py_S_Pz_S_c;
  abcd[104] = 4.0E0*I_ERI_F2xz_S_Pz_S_ac-2.0E0*1*I_ERI_Pz_S_Pz_S_c;
  abcd[105] = 4.0E0*I_ERI_Fx2y_S_Pz_S_ac;
  abcd[106] = 4.0E0*I_ERI_Fxyz_S_Pz_S_ac;
  abcd[107] = 4.0E0*I_ERI_Fx2z_S_Pz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_day_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_c
   ************************************************************/
  abcd[108] = 4.0E0*I_ERI_F2xy_S_Px_S_ac;
  abcd[109] = 4.0E0*I_ERI_Fx2y_S_Px_S_ac-2.0E0*1*I_ERI_Px_S_Px_S_c;
  abcd[110] = 4.0E0*I_ERI_Fxyz_S_Px_S_ac;
  abcd[111] = 4.0E0*I_ERI_F3y_S_Px_S_ac-2.0E0*2*I_ERI_Py_S_Px_S_c;
  abcd[112] = 4.0E0*I_ERI_F2yz_S_Px_S_ac-2.0E0*1*I_ERI_Pz_S_Px_S_c;
  abcd[113] = 4.0E0*I_ERI_Fy2z_S_Px_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_day_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_c
   ************************************************************/
  abcd[114] = 4.0E0*I_ERI_F2xy_S_Py_S_ac;
  abcd[115] = 4.0E0*I_ERI_Fx2y_S_Py_S_ac-2.0E0*1*I_ERI_Px_S_Py_S_c;
  abcd[116] = 4.0E0*I_ERI_Fxyz_S_Py_S_ac;
  abcd[117] = 4.0E0*I_ERI_F3y_S_Py_S_ac-2.0E0*2*I_ERI_Py_S_Py_S_c;
  abcd[118] = 4.0E0*I_ERI_F2yz_S_Py_S_ac-2.0E0*1*I_ERI_Pz_S_Py_S_c;
  abcd[119] = 4.0E0*I_ERI_Fy2z_S_Py_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_day_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_c
   ************************************************************/
  abcd[120] = 4.0E0*I_ERI_F2xy_S_Pz_S_ac;
  abcd[121] = 4.0E0*I_ERI_Fx2y_S_Pz_S_ac-2.0E0*1*I_ERI_Px_S_Pz_S_c;
  abcd[122] = 4.0E0*I_ERI_Fxyz_S_Pz_S_ac;
  abcd[123] = 4.0E0*I_ERI_F3y_S_Pz_S_ac-2.0E0*2*I_ERI_Py_S_Pz_S_c;
  abcd[124] = 4.0E0*I_ERI_F2yz_S_Pz_S_ac-2.0E0*1*I_ERI_Pz_S_Pz_S_c;
  abcd[125] = 4.0E0*I_ERI_Fy2z_S_Pz_S_ac;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_daz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_c
   ************************************************************/
  abcd[126] = 4.0E0*I_ERI_F2xz_S_Px_S_ac;
  abcd[127] = 4.0E0*I_ERI_Fxyz_S_Px_S_ac;
  abcd[128] = 4.0E0*I_ERI_Fx2z_S_Px_S_ac-2.0E0*1*I_ERI_Px_S_Px_S_c;
  abcd[129] = 4.0E0*I_ERI_F2yz_S_Px_S_ac;
  abcd[130] = 4.0E0*I_ERI_Fy2z_S_Px_S_ac-2.0E0*1*I_ERI_Py_S_Px_S_c;
  abcd[131] = 4.0E0*I_ERI_F3z_S_Px_S_ac-2.0E0*2*I_ERI_Pz_S_Px_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_daz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_c
   ************************************************************/
  abcd[132] = 4.0E0*I_ERI_F2xz_S_Py_S_ac;
  abcd[133] = 4.0E0*I_ERI_Fxyz_S_Py_S_ac;
  abcd[134] = 4.0E0*I_ERI_Fx2z_S_Py_S_ac-2.0E0*1*I_ERI_Px_S_Py_S_c;
  abcd[135] = 4.0E0*I_ERI_F2yz_S_Py_S_ac;
  abcd[136] = 4.0E0*I_ERI_Fy2z_S_Py_S_ac-2.0E0*1*I_ERI_Py_S_Py_S_c;
  abcd[137] = 4.0E0*I_ERI_F3z_S_Py_S_ac-2.0E0*2*I_ERI_Pz_S_Py_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_daz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_P_S_ac
   * RHS shell quartet name: SQ_ERI_P_S_P_S_c
   ************************************************************/
  abcd[138] = 4.0E0*I_ERI_F2xz_S_Pz_S_ac;
  abcd[139] = 4.0E0*I_ERI_Fxyz_S_Pz_S_ac;
  abcd[140] = 4.0E0*I_ERI_Fx2z_S_Pz_S_ac-2.0E0*1*I_ERI_Px_S_Pz_S_c;
  abcd[141] = 4.0E0*I_ERI_F2yz_S_Pz_S_ac;
  abcd[142] = 4.0E0*I_ERI_Fy2z_S_Pz_S_ac-2.0E0*1*I_ERI_Py_S_Pz_S_c;
  abcd[143] = 4.0E0*I_ERI_F3z_S_Pz_S_ac-2.0E0*2*I_ERI_Pz_S_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbx_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_b
   ************************************************************/
  abcd[144] = 4.0E0*I_ERI_D2x_D2x_S_S_bb-2.0E0*1*I_ERI_D2x_S_S_S_b;
  abcd[145] = 4.0E0*I_ERI_Dxy_D2x_S_S_bb-2.0E0*1*I_ERI_Dxy_S_S_S_b;
  abcd[146] = 4.0E0*I_ERI_Dxz_D2x_S_S_bb-2.0E0*1*I_ERI_Dxz_S_S_S_b;
  abcd[147] = 4.0E0*I_ERI_D2y_D2x_S_S_bb-2.0E0*1*I_ERI_D2y_S_S_S_b;
  abcd[148] = 4.0E0*I_ERI_Dyz_D2x_S_S_bb-2.0E0*1*I_ERI_Dyz_S_S_S_b;
  abcd[149] = 4.0E0*I_ERI_D2z_D2x_S_S_bb-2.0E0*1*I_ERI_D2z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbx_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_b
   ************************************************************/
  abcd[150] = 4.0E0*I_ERI_D2x_Dxy_S_S_bb;
  abcd[151] = 4.0E0*I_ERI_Dxy_Dxy_S_S_bb;
  abcd[152] = 4.0E0*I_ERI_Dxz_Dxy_S_S_bb;
  abcd[153] = 4.0E0*I_ERI_D2y_Dxy_S_S_bb;
  abcd[154] = 4.0E0*I_ERI_Dyz_Dxy_S_S_bb;
  abcd[155] = 4.0E0*I_ERI_D2z_Dxy_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbx_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_b
   ************************************************************/
  abcd[156] = 4.0E0*I_ERI_D2x_Dxz_S_S_bb;
  abcd[157] = 4.0E0*I_ERI_Dxy_Dxz_S_S_bb;
  abcd[158] = 4.0E0*I_ERI_Dxz_Dxz_S_S_bb;
  abcd[159] = 4.0E0*I_ERI_D2y_Dxz_S_S_bb;
  abcd[160] = 4.0E0*I_ERI_Dyz_Dxz_S_S_bb;
  abcd[161] = 4.0E0*I_ERI_D2z_Dxz_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dby_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_b
   ************************************************************/
  abcd[162] = 4.0E0*I_ERI_D2x_D2y_S_S_bb-2.0E0*1*I_ERI_D2x_S_S_S_b;
  abcd[163] = 4.0E0*I_ERI_Dxy_D2y_S_S_bb-2.0E0*1*I_ERI_Dxy_S_S_S_b;
  abcd[164] = 4.0E0*I_ERI_Dxz_D2y_S_S_bb-2.0E0*1*I_ERI_Dxz_S_S_S_b;
  abcd[165] = 4.0E0*I_ERI_D2y_D2y_S_S_bb-2.0E0*1*I_ERI_D2y_S_S_S_b;
  abcd[166] = 4.0E0*I_ERI_Dyz_D2y_S_S_bb-2.0E0*1*I_ERI_Dyz_S_S_S_b;
  abcd[167] = 4.0E0*I_ERI_D2z_D2y_S_S_bb-2.0E0*1*I_ERI_D2z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dby_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_b
   ************************************************************/
  abcd[168] = 4.0E0*I_ERI_D2x_Dyz_S_S_bb;
  abcd[169] = 4.0E0*I_ERI_Dxy_Dyz_S_S_bb;
  abcd[170] = 4.0E0*I_ERI_Dxz_Dyz_S_S_bb;
  abcd[171] = 4.0E0*I_ERI_D2y_Dyz_S_S_bb;
  abcd[172] = 4.0E0*I_ERI_Dyz_Dyz_S_S_bb;
  abcd[173] = 4.0E0*I_ERI_D2z_Dyz_S_S_bb;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbz_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_D_S_S_bb
   * RHS shell quartet name: SQ_ERI_D_S_S_S_b
   ************************************************************/
  abcd[174] = 4.0E0*I_ERI_D2x_D2z_S_S_bb-2.0E0*1*I_ERI_D2x_S_S_S_b;
  abcd[175] = 4.0E0*I_ERI_Dxy_D2z_S_S_bb-2.0E0*1*I_ERI_Dxy_S_S_S_b;
  abcd[176] = 4.0E0*I_ERI_Dxz_D2z_S_S_bb-2.0E0*1*I_ERI_Dxz_S_S_S_b;
  abcd[177] = 4.0E0*I_ERI_D2y_D2z_S_S_bb-2.0E0*1*I_ERI_D2y_S_S_S_b;
  abcd[178] = 4.0E0*I_ERI_Dyz_D2z_S_S_bb-2.0E0*1*I_ERI_Dyz_S_S_S_b;
  abcd[179] = 4.0E0*I_ERI_D2z_D2z_S_S_bb-2.0E0*1*I_ERI_D2z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bc
   ************************************************************/
  abcd[180] = 4.0E0*I_ERI_D2x_Px_Px_S_bc;
  abcd[181] = 4.0E0*I_ERI_Dxy_Px_Px_S_bc;
  abcd[182] = 4.0E0*I_ERI_Dxz_Px_Px_S_bc;
  abcd[183] = 4.0E0*I_ERI_D2y_Px_Px_S_bc;
  abcd[184] = 4.0E0*I_ERI_Dyz_Px_Px_S_bc;
  abcd[185] = 4.0E0*I_ERI_D2z_Px_Px_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bc
   ************************************************************/
  abcd[186] = 4.0E0*I_ERI_D2x_Px_Py_S_bc;
  abcd[187] = 4.0E0*I_ERI_Dxy_Px_Py_S_bc;
  abcd[188] = 4.0E0*I_ERI_Dxz_Px_Py_S_bc;
  abcd[189] = 4.0E0*I_ERI_D2y_Px_Py_S_bc;
  abcd[190] = 4.0E0*I_ERI_Dyz_Px_Py_S_bc;
  abcd[191] = 4.0E0*I_ERI_D2z_Px_Py_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bc
   ************************************************************/
  abcd[192] = 4.0E0*I_ERI_D2x_Px_Pz_S_bc;
  abcd[193] = 4.0E0*I_ERI_Dxy_Px_Pz_S_bc;
  abcd[194] = 4.0E0*I_ERI_Dxz_Px_Pz_S_bc;
  abcd[195] = 4.0E0*I_ERI_D2y_Px_Pz_S_bc;
  abcd[196] = 4.0E0*I_ERI_Dyz_Px_Pz_S_bc;
  abcd[197] = 4.0E0*I_ERI_D2z_Px_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dby_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bc
   ************************************************************/
  abcd[198] = 4.0E0*I_ERI_D2x_Py_Px_S_bc;
  abcd[199] = 4.0E0*I_ERI_Dxy_Py_Px_S_bc;
  abcd[200] = 4.0E0*I_ERI_Dxz_Py_Px_S_bc;
  abcd[201] = 4.0E0*I_ERI_D2y_Py_Px_S_bc;
  abcd[202] = 4.0E0*I_ERI_Dyz_Py_Px_S_bc;
  abcd[203] = 4.0E0*I_ERI_D2z_Py_Px_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dby_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bc
   ************************************************************/
  abcd[204] = 4.0E0*I_ERI_D2x_Py_Py_S_bc;
  abcd[205] = 4.0E0*I_ERI_Dxy_Py_Py_S_bc;
  abcd[206] = 4.0E0*I_ERI_Dxz_Py_Py_S_bc;
  abcd[207] = 4.0E0*I_ERI_D2y_Py_Py_S_bc;
  abcd[208] = 4.0E0*I_ERI_Dyz_Py_Py_S_bc;
  abcd[209] = 4.0E0*I_ERI_D2z_Py_Py_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dby_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bc
   ************************************************************/
  abcd[210] = 4.0E0*I_ERI_D2x_Py_Pz_S_bc;
  abcd[211] = 4.0E0*I_ERI_Dxy_Py_Pz_S_bc;
  abcd[212] = 4.0E0*I_ERI_Dxz_Py_Pz_S_bc;
  abcd[213] = 4.0E0*I_ERI_D2y_Py_Pz_S_bc;
  abcd[214] = 4.0E0*I_ERI_Dyz_Py_Pz_S_bc;
  abcd[215] = 4.0E0*I_ERI_D2z_Py_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbz_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bc
   ************************************************************/
  abcd[216] = 4.0E0*I_ERI_D2x_Pz_Px_S_bc;
  abcd[217] = 4.0E0*I_ERI_Dxy_Pz_Px_S_bc;
  abcd[218] = 4.0E0*I_ERI_Dxz_Pz_Px_S_bc;
  abcd[219] = 4.0E0*I_ERI_D2y_Pz_Px_S_bc;
  abcd[220] = 4.0E0*I_ERI_Dyz_Pz_Px_S_bc;
  abcd[221] = 4.0E0*I_ERI_D2z_Pz_Px_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbz_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bc
   ************************************************************/
  abcd[222] = 4.0E0*I_ERI_D2x_Pz_Py_S_bc;
  abcd[223] = 4.0E0*I_ERI_Dxy_Pz_Py_S_bc;
  abcd[224] = 4.0E0*I_ERI_Dxz_Pz_Py_S_bc;
  abcd[225] = 4.0E0*I_ERI_D2y_Pz_Py_S_bc;
  abcd[226] = 4.0E0*I_ERI_Dyz_Pz_Py_S_bc;
  abcd[227] = 4.0E0*I_ERI_D2z_Pz_Py_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dbz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_P_P_S_bc
   ************************************************************/
  abcd[228] = 4.0E0*I_ERI_D2x_Pz_Pz_S_bc;
  abcd[229] = 4.0E0*I_ERI_Dxy_Pz_Pz_S_bc;
  abcd[230] = 4.0E0*I_ERI_Dxz_Pz_Pz_S_bc;
  abcd[231] = 4.0E0*I_ERI_D2y_Pz_Pz_S_bc;
  abcd[232] = 4.0E0*I_ERI_Dyz_Pz_Pz_S_bc;
  abcd[233] = 4.0E0*I_ERI_D2z_Pz_Pz_S_bc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dcx_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_c
   ************************************************************/
  abcd[234] = 4.0E0*I_ERI_D2x_S_D2x_S_cc-2.0E0*1*I_ERI_D2x_S_S_S_c;
  abcd[235] = 4.0E0*I_ERI_Dxy_S_D2x_S_cc-2.0E0*1*I_ERI_Dxy_S_S_S_c;
  abcd[236] = 4.0E0*I_ERI_Dxz_S_D2x_S_cc-2.0E0*1*I_ERI_Dxz_S_S_S_c;
  abcd[237] = 4.0E0*I_ERI_D2y_S_D2x_S_cc-2.0E0*1*I_ERI_D2y_S_S_S_c;
  abcd[238] = 4.0E0*I_ERI_Dyz_S_D2x_S_cc-2.0E0*1*I_ERI_Dyz_S_S_S_c;
  abcd[239] = 4.0E0*I_ERI_D2z_S_D2x_S_cc-2.0E0*1*I_ERI_D2z_S_S_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dcx_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_c
   ************************************************************/
  abcd[240] = 4.0E0*I_ERI_D2x_S_Dxy_S_cc;
  abcd[241] = 4.0E0*I_ERI_Dxy_S_Dxy_S_cc;
  abcd[242] = 4.0E0*I_ERI_Dxz_S_Dxy_S_cc;
  abcd[243] = 4.0E0*I_ERI_D2y_S_Dxy_S_cc;
  abcd[244] = 4.0E0*I_ERI_Dyz_S_Dxy_S_cc;
  abcd[245] = 4.0E0*I_ERI_D2z_S_Dxy_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dcx_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_c
   ************************************************************/
  abcd[246] = 4.0E0*I_ERI_D2x_S_Dxz_S_cc;
  abcd[247] = 4.0E0*I_ERI_Dxy_S_Dxz_S_cc;
  abcd[248] = 4.0E0*I_ERI_Dxz_S_Dxz_S_cc;
  abcd[249] = 4.0E0*I_ERI_D2y_S_Dxz_S_cc;
  abcd[250] = 4.0E0*I_ERI_Dyz_S_Dxz_S_cc;
  abcd[251] = 4.0E0*I_ERI_D2z_S_Dxz_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dcy_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_c
   ************************************************************/
  abcd[252] = 4.0E0*I_ERI_D2x_S_D2y_S_cc-2.0E0*1*I_ERI_D2x_S_S_S_c;
  abcd[253] = 4.0E0*I_ERI_Dxy_S_D2y_S_cc-2.0E0*1*I_ERI_Dxy_S_S_S_c;
  abcd[254] = 4.0E0*I_ERI_Dxz_S_D2y_S_cc-2.0E0*1*I_ERI_Dxz_S_S_S_c;
  abcd[255] = 4.0E0*I_ERI_D2y_S_D2y_S_cc-2.0E0*1*I_ERI_D2y_S_S_S_c;
  abcd[256] = 4.0E0*I_ERI_Dyz_S_D2y_S_cc-2.0E0*1*I_ERI_Dyz_S_S_S_c;
  abcd[257] = 4.0E0*I_ERI_D2z_S_D2y_S_cc-2.0E0*1*I_ERI_D2z_S_S_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dcy_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_c
   ************************************************************/
  abcd[258] = 4.0E0*I_ERI_D2x_S_Dyz_S_cc;
  abcd[259] = 4.0E0*I_ERI_Dxy_S_Dyz_S_cc;
  abcd[260] = 4.0E0*I_ERI_Dxz_S_Dyz_S_cc;
  abcd[261] = 4.0E0*I_ERI_D2y_S_Dyz_S_cc;
  abcd[262] = 4.0E0*I_ERI_Dyz_S_Dyz_S_cc;
  abcd[263] = 4.0E0*I_ERI_D2z_S_Dyz_S_cc;

  /************************************************************
   * shell quartet name: SQ_ERI_D_S_S_S_dcz_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_D_S_cc
   * RHS shell quartet name: SQ_ERI_D_S_S_S_c
   ************************************************************/
  abcd[264] = 4.0E0*I_ERI_D2x_S_D2z_S_cc-2.0E0*1*I_ERI_D2x_S_S_S_c;
  abcd[265] = 4.0E0*I_ERI_Dxy_S_D2z_S_cc-2.0E0*1*I_ERI_Dxy_S_S_S_c;
  abcd[266] = 4.0E0*I_ERI_Dxz_S_D2z_S_cc-2.0E0*1*I_ERI_Dxz_S_S_S_c;
  abcd[267] = 4.0E0*I_ERI_D2y_S_D2z_S_cc-2.0E0*1*I_ERI_D2y_S_S_S_c;
  abcd[268] = 4.0E0*I_ERI_Dyz_S_D2z_S_cc-2.0E0*1*I_ERI_Dyz_S_S_S_c;
  abcd[269] = 4.0E0*I_ERI_D2z_S_D2z_S_cc-2.0E0*1*I_ERI_D2z_S_S_S_c;
}
