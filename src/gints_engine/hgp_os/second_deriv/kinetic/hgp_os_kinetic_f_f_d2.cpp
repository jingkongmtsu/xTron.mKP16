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
// BRA1  BRA1
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_kinetic_f_f_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_H5x_F3x_aa = 0.0E0;
  Double I_KINETIC_H4xy_F3x_aa = 0.0E0;
  Double I_KINETIC_H4xz_F3x_aa = 0.0E0;
  Double I_KINETIC_H3x2y_F3x_aa = 0.0E0;
  Double I_KINETIC_H3xyz_F3x_aa = 0.0E0;
  Double I_KINETIC_H3x2z_F3x_aa = 0.0E0;
  Double I_KINETIC_H2x3y_F3x_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_F3x_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_F3x_aa = 0.0E0;
  Double I_KINETIC_H2x3z_F3x_aa = 0.0E0;
  Double I_KINETIC_Hx4y_F3x_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_F3x_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_F3x_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_F3x_aa = 0.0E0;
  Double I_KINETIC_Hx4z_F3x_aa = 0.0E0;
  Double I_KINETIC_H5y_F3x_aa = 0.0E0;
  Double I_KINETIC_H4yz_F3x_aa = 0.0E0;
  Double I_KINETIC_H3y2z_F3x_aa = 0.0E0;
  Double I_KINETIC_H2y3z_F3x_aa = 0.0E0;
  Double I_KINETIC_Hy4z_F3x_aa = 0.0E0;
  Double I_KINETIC_H5z_F3x_aa = 0.0E0;
  Double I_KINETIC_H5x_F2xy_aa = 0.0E0;
  Double I_KINETIC_H4xy_F2xy_aa = 0.0E0;
  Double I_KINETIC_H4xz_F2xy_aa = 0.0E0;
  Double I_KINETIC_H3x2y_F2xy_aa = 0.0E0;
  Double I_KINETIC_H3xyz_F2xy_aa = 0.0E0;
  Double I_KINETIC_H3x2z_F2xy_aa = 0.0E0;
  Double I_KINETIC_H2x3y_F2xy_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_F2xy_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_F2xy_aa = 0.0E0;
  Double I_KINETIC_H2x3z_F2xy_aa = 0.0E0;
  Double I_KINETIC_Hx4y_F2xy_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_F2xy_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_F2xy_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_F2xy_aa = 0.0E0;
  Double I_KINETIC_Hx4z_F2xy_aa = 0.0E0;
  Double I_KINETIC_H5y_F2xy_aa = 0.0E0;
  Double I_KINETIC_H4yz_F2xy_aa = 0.0E0;
  Double I_KINETIC_H3y2z_F2xy_aa = 0.0E0;
  Double I_KINETIC_H2y3z_F2xy_aa = 0.0E0;
  Double I_KINETIC_Hy4z_F2xy_aa = 0.0E0;
  Double I_KINETIC_H5z_F2xy_aa = 0.0E0;
  Double I_KINETIC_H5x_F2xz_aa = 0.0E0;
  Double I_KINETIC_H4xy_F2xz_aa = 0.0E0;
  Double I_KINETIC_H4xz_F2xz_aa = 0.0E0;
  Double I_KINETIC_H3x2y_F2xz_aa = 0.0E0;
  Double I_KINETIC_H3xyz_F2xz_aa = 0.0E0;
  Double I_KINETIC_H3x2z_F2xz_aa = 0.0E0;
  Double I_KINETIC_H2x3y_F2xz_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_F2xz_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_F2xz_aa = 0.0E0;
  Double I_KINETIC_H2x3z_F2xz_aa = 0.0E0;
  Double I_KINETIC_Hx4y_F2xz_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_F2xz_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_F2xz_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_F2xz_aa = 0.0E0;
  Double I_KINETIC_Hx4z_F2xz_aa = 0.0E0;
  Double I_KINETIC_H5y_F2xz_aa = 0.0E0;
  Double I_KINETIC_H4yz_F2xz_aa = 0.0E0;
  Double I_KINETIC_H3y2z_F2xz_aa = 0.0E0;
  Double I_KINETIC_H2y3z_F2xz_aa = 0.0E0;
  Double I_KINETIC_Hy4z_F2xz_aa = 0.0E0;
  Double I_KINETIC_H5z_F2xz_aa = 0.0E0;
  Double I_KINETIC_H5x_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H4xy_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H4xz_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H3x2y_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H3xyz_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H3x2z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H2x3y_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H2x3z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_Hx4y_Fx2y_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_Fx2y_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_Hx4z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H5y_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H4yz_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H3y2z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H2y3z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_Hy4z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H5z_Fx2y_aa = 0.0E0;
  Double I_KINETIC_H5x_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H4xy_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H4xz_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H3x2y_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H3xyz_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H3x2z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H2x3y_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H2x3z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_Hx4y_Fxyz_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_Fxyz_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_Hx4z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H5y_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H4yz_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H3y2z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H2y3z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_Hy4z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H5z_Fxyz_aa = 0.0E0;
  Double I_KINETIC_H5x_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H4xy_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H4xz_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H3x2y_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H3xyz_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H3x2z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H2x3y_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H2x3z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_Hx4y_Fx2z_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_Fx2z_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_Hx4z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H5y_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H4yz_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H3y2z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H2y3z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_Hy4z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H5z_Fx2z_aa = 0.0E0;
  Double I_KINETIC_H5x_F3y_aa = 0.0E0;
  Double I_KINETIC_H4xy_F3y_aa = 0.0E0;
  Double I_KINETIC_H4xz_F3y_aa = 0.0E0;
  Double I_KINETIC_H3x2y_F3y_aa = 0.0E0;
  Double I_KINETIC_H3xyz_F3y_aa = 0.0E0;
  Double I_KINETIC_H3x2z_F3y_aa = 0.0E0;
  Double I_KINETIC_H2x3y_F3y_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_F3y_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_F3y_aa = 0.0E0;
  Double I_KINETIC_H2x3z_F3y_aa = 0.0E0;
  Double I_KINETIC_Hx4y_F3y_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_F3y_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_F3y_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_F3y_aa = 0.0E0;
  Double I_KINETIC_Hx4z_F3y_aa = 0.0E0;
  Double I_KINETIC_H5y_F3y_aa = 0.0E0;
  Double I_KINETIC_H4yz_F3y_aa = 0.0E0;
  Double I_KINETIC_H3y2z_F3y_aa = 0.0E0;
  Double I_KINETIC_H2y3z_F3y_aa = 0.0E0;
  Double I_KINETIC_Hy4z_F3y_aa = 0.0E0;
  Double I_KINETIC_H5z_F3y_aa = 0.0E0;
  Double I_KINETIC_H5x_F2yz_aa = 0.0E0;
  Double I_KINETIC_H4xy_F2yz_aa = 0.0E0;
  Double I_KINETIC_H4xz_F2yz_aa = 0.0E0;
  Double I_KINETIC_H3x2y_F2yz_aa = 0.0E0;
  Double I_KINETIC_H3xyz_F2yz_aa = 0.0E0;
  Double I_KINETIC_H3x2z_F2yz_aa = 0.0E0;
  Double I_KINETIC_H2x3y_F2yz_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_F2yz_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_F2yz_aa = 0.0E0;
  Double I_KINETIC_H2x3z_F2yz_aa = 0.0E0;
  Double I_KINETIC_Hx4y_F2yz_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_F2yz_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_F2yz_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_F2yz_aa = 0.0E0;
  Double I_KINETIC_Hx4z_F2yz_aa = 0.0E0;
  Double I_KINETIC_H5y_F2yz_aa = 0.0E0;
  Double I_KINETIC_H4yz_F2yz_aa = 0.0E0;
  Double I_KINETIC_H3y2z_F2yz_aa = 0.0E0;
  Double I_KINETIC_H2y3z_F2yz_aa = 0.0E0;
  Double I_KINETIC_Hy4z_F2yz_aa = 0.0E0;
  Double I_KINETIC_H5z_F2yz_aa = 0.0E0;
  Double I_KINETIC_H5x_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H4xy_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H4xz_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H3x2y_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H3xyz_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H3x2z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H2x3y_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H2x3z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_Hx4y_Fy2z_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_Fy2z_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_Hx4z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H5y_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H4yz_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H3y2z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H2y3z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_Hy4z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H5z_Fy2z_aa = 0.0E0;
  Double I_KINETIC_H5x_F3z_aa = 0.0E0;
  Double I_KINETIC_H4xy_F3z_aa = 0.0E0;
  Double I_KINETIC_H4xz_F3z_aa = 0.0E0;
  Double I_KINETIC_H3x2y_F3z_aa = 0.0E0;
  Double I_KINETIC_H3xyz_F3z_aa = 0.0E0;
  Double I_KINETIC_H3x2z_F3z_aa = 0.0E0;
  Double I_KINETIC_H2x3y_F3z_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_F3z_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_F3z_aa = 0.0E0;
  Double I_KINETIC_H2x3z_F3z_aa = 0.0E0;
  Double I_KINETIC_Hx4y_F3z_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_F3z_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_F3z_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_F3z_aa = 0.0E0;
  Double I_KINETIC_Hx4z_F3z_aa = 0.0E0;
  Double I_KINETIC_H5y_F3z_aa = 0.0E0;
  Double I_KINETIC_H4yz_F3z_aa = 0.0E0;
  Double I_KINETIC_H3y2z_F3z_aa = 0.0E0;
  Double I_KINETIC_H2y3z_F3z_aa = 0.0E0;
  Double I_KINETIC_Hy4z_F3z_aa = 0.0E0;
  Double I_KINETIC_H5z_F3z_aa = 0.0E0;
  Double I_KINETIC_F3x_F3x_a = 0.0E0;
  Double I_KINETIC_F2xy_F3x_a = 0.0E0;
  Double I_KINETIC_F2xz_F3x_a = 0.0E0;
  Double I_KINETIC_Fx2y_F3x_a = 0.0E0;
  Double I_KINETIC_Fxyz_F3x_a = 0.0E0;
  Double I_KINETIC_Fx2z_F3x_a = 0.0E0;
  Double I_KINETIC_F3y_F3x_a = 0.0E0;
  Double I_KINETIC_F2yz_F3x_a = 0.0E0;
  Double I_KINETIC_Fy2z_F3x_a = 0.0E0;
  Double I_KINETIC_F3z_F3x_a = 0.0E0;
  Double I_KINETIC_F3x_F2xy_a = 0.0E0;
  Double I_KINETIC_F2xy_F2xy_a = 0.0E0;
  Double I_KINETIC_F2xz_F2xy_a = 0.0E0;
  Double I_KINETIC_Fx2y_F2xy_a = 0.0E0;
  Double I_KINETIC_Fxyz_F2xy_a = 0.0E0;
  Double I_KINETIC_Fx2z_F2xy_a = 0.0E0;
  Double I_KINETIC_F3y_F2xy_a = 0.0E0;
  Double I_KINETIC_F2yz_F2xy_a = 0.0E0;
  Double I_KINETIC_Fy2z_F2xy_a = 0.0E0;
  Double I_KINETIC_F3z_F2xy_a = 0.0E0;
  Double I_KINETIC_F3x_F2xz_a = 0.0E0;
  Double I_KINETIC_F2xy_F2xz_a = 0.0E0;
  Double I_KINETIC_F2xz_F2xz_a = 0.0E0;
  Double I_KINETIC_Fx2y_F2xz_a = 0.0E0;
  Double I_KINETIC_Fxyz_F2xz_a = 0.0E0;
  Double I_KINETIC_Fx2z_F2xz_a = 0.0E0;
  Double I_KINETIC_F3y_F2xz_a = 0.0E0;
  Double I_KINETIC_F2yz_F2xz_a = 0.0E0;
  Double I_KINETIC_Fy2z_F2xz_a = 0.0E0;
  Double I_KINETIC_F3z_F2xz_a = 0.0E0;
  Double I_KINETIC_F3x_Fx2y_a = 0.0E0;
  Double I_KINETIC_F2xy_Fx2y_a = 0.0E0;
  Double I_KINETIC_F2xz_Fx2y_a = 0.0E0;
  Double I_KINETIC_Fx2y_Fx2y_a = 0.0E0;
  Double I_KINETIC_Fxyz_Fx2y_a = 0.0E0;
  Double I_KINETIC_Fx2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_F3y_Fx2y_a = 0.0E0;
  Double I_KINETIC_F2yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_Fy2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_F3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_F3x_Fxyz_a = 0.0E0;
  Double I_KINETIC_F2xy_Fxyz_a = 0.0E0;
  Double I_KINETIC_F2xz_Fxyz_a = 0.0E0;
  Double I_KINETIC_Fx2y_Fxyz_a = 0.0E0;
  Double I_KINETIC_Fxyz_Fxyz_a = 0.0E0;
  Double I_KINETIC_Fx2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_F3y_Fxyz_a = 0.0E0;
  Double I_KINETIC_F2yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_Fy2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_F3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_F3x_Fx2z_a = 0.0E0;
  Double I_KINETIC_F2xy_Fx2z_a = 0.0E0;
  Double I_KINETIC_F2xz_Fx2z_a = 0.0E0;
  Double I_KINETIC_Fx2y_Fx2z_a = 0.0E0;
  Double I_KINETIC_Fxyz_Fx2z_a = 0.0E0;
  Double I_KINETIC_Fx2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_F3y_Fx2z_a = 0.0E0;
  Double I_KINETIC_F2yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_Fy2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_F3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_F3x_F3y_a = 0.0E0;
  Double I_KINETIC_F2xy_F3y_a = 0.0E0;
  Double I_KINETIC_F2xz_F3y_a = 0.0E0;
  Double I_KINETIC_Fx2y_F3y_a = 0.0E0;
  Double I_KINETIC_Fxyz_F3y_a = 0.0E0;
  Double I_KINETIC_Fx2z_F3y_a = 0.0E0;
  Double I_KINETIC_F3y_F3y_a = 0.0E0;
  Double I_KINETIC_F2yz_F3y_a = 0.0E0;
  Double I_KINETIC_Fy2z_F3y_a = 0.0E0;
  Double I_KINETIC_F3z_F3y_a = 0.0E0;
  Double I_KINETIC_F3x_F2yz_a = 0.0E0;
  Double I_KINETIC_F2xy_F2yz_a = 0.0E0;
  Double I_KINETIC_F2xz_F2yz_a = 0.0E0;
  Double I_KINETIC_Fx2y_F2yz_a = 0.0E0;
  Double I_KINETIC_Fxyz_F2yz_a = 0.0E0;
  Double I_KINETIC_Fx2z_F2yz_a = 0.0E0;
  Double I_KINETIC_F3y_F2yz_a = 0.0E0;
  Double I_KINETIC_F2yz_F2yz_a = 0.0E0;
  Double I_KINETIC_Fy2z_F2yz_a = 0.0E0;
  Double I_KINETIC_F3z_F2yz_a = 0.0E0;
  Double I_KINETIC_F3x_Fy2z_a = 0.0E0;
  Double I_KINETIC_F2xy_Fy2z_a = 0.0E0;
  Double I_KINETIC_F2xz_Fy2z_a = 0.0E0;
  Double I_KINETIC_Fx2y_Fy2z_a = 0.0E0;
  Double I_KINETIC_Fxyz_Fy2z_a = 0.0E0;
  Double I_KINETIC_Fx2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_F3y_Fy2z_a = 0.0E0;
  Double I_KINETIC_F2yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_Fy2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_F3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_F3x_F3z_a = 0.0E0;
  Double I_KINETIC_F2xy_F3z_a = 0.0E0;
  Double I_KINETIC_F2xz_F3z_a = 0.0E0;
  Double I_KINETIC_Fx2y_F3z_a = 0.0E0;
  Double I_KINETIC_Fxyz_F3z_a = 0.0E0;
  Double I_KINETIC_Fx2z_F3z_a = 0.0E0;
  Double I_KINETIC_F3y_F3z_a = 0.0E0;
  Double I_KINETIC_F2yz_F3z_a = 0.0E0;
  Double I_KINETIC_Fy2z_F3z_a = 0.0E0;
  Double I_KINETIC_F3z_F3z_a = 0.0E0;
  Double I_KINETIC_Px_F3x = 0.0E0;
  Double I_KINETIC_Py_F3x = 0.0E0;
  Double I_KINETIC_Pz_F3x = 0.0E0;
  Double I_KINETIC_Px_F2xy = 0.0E0;
  Double I_KINETIC_Py_F2xy = 0.0E0;
  Double I_KINETIC_Pz_F2xy = 0.0E0;
  Double I_KINETIC_Px_F2xz = 0.0E0;
  Double I_KINETIC_Py_F2xz = 0.0E0;
  Double I_KINETIC_Pz_F2xz = 0.0E0;
  Double I_KINETIC_Px_Fx2y = 0.0E0;
  Double I_KINETIC_Py_Fx2y = 0.0E0;
  Double I_KINETIC_Pz_Fx2y = 0.0E0;
  Double I_KINETIC_Px_Fxyz = 0.0E0;
  Double I_KINETIC_Py_Fxyz = 0.0E0;
  Double I_KINETIC_Pz_Fxyz = 0.0E0;
  Double I_KINETIC_Px_Fx2z = 0.0E0;
  Double I_KINETIC_Py_Fx2z = 0.0E0;
  Double I_KINETIC_Pz_Fx2z = 0.0E0;
  Double I_KINETIC_Px_F3y = 0.0E0;
  Double I_KINETIC_Py_F3y = 0.0E0;
  Double I_KINETIC_Pz_F3y = 0.0E0;
  Double I_KINETIC_Px_F2yz = 0.0E0;
  Double I_KINETIC_Py_F2yz = 0.0E0;
  Double I_KINETIC_Pz_F2yz = 0.0E0;
  Double I_KINETIC_Px_Fy2z = 0.0E0;
  Double I_KINETIC_Py_Fy2z = 0.0E0;
  Double I_KINETIC_Pz_Fy2z = 0.0E0;
  Double I_KINETIC_Px_F3z = 0.0E0;
  Double I_KINETIC_Py_F3z = 0.0E0;
  Double I_KINETIC_Pz_F3z = 0.0E0;

  Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double onedz = iexp[ip2];
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double xi    = alpha*beta*onedz;
    Double twoxi = 2.0E0*xi;
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    Double adz   = alpha*onedz;
    Double bdz   = beta*onedz;
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
    Double I_KINETIC_S_S_vrr = ic2*fbra*xi*(3.0E0-twoxi*AB2);
    Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;
    if(fabs(I_KINETIC_S_S_vrr)<THRESHOLD_MATH) continue;


    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_S_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_S_Px_vrr = PBX*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_S_Py_vrr = PBY*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_S_Pz_vrr = PBZ*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_S_vrr = PAX*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Py_S_vrr = PAY*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_S_vrr = PAZ*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_Px_vrr = PBX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Py_Px_vrr = PBX*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_Px_vrr = PBX*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_Px_Py_vrr = PBY*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Py_Py_vrr = PBY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_Py_vrr = PBY*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_Px_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Py_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_S_D2x_vrr = PBX*I_TWOBODYOVERLAP_S_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_S_Dxy_vrr = PBY*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_S_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_S_D2y_vrr = PBY*I_TWOBODYOVERLAP_S_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_S_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_S_D2z_vrr = PBZ*I_TWOBODYOVERLAP_S_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_S_vrr = PAZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_D2x_vrr = PAX*I_TWOBODYOVERLAP_S_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_Py_D2x_vrr = PAY*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Pz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Px_Dxy_vrr = PAX*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Py_Dxy_vrr = PAY*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_Pz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Px_Dxz_vrr = PAX*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_Py_Dxz_vrr = PAY*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Pz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_Px_D2y_vrr = PAX*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Py_D2y_vrr = PAY*I_TWOBODYOVERLAP_S_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Pz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Px_Dyz_vrr = PAX*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Py_Dyz_vrr = PAY*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_Pz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Px_D2z_vrr = PAX*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Py_D2z_vrr = PAY*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Pz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_S_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_Px_vrr = PAX*I_TWOBODYOVERLAP_Px_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Px_vrr = PAY*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_Px_vrr = PAY*I_TWOBODYOVERLAP_Py_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_Py_vrr = PAX*I_TWOBODYOVERLAP_Px_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Py_vrr = PAY*I_TWOBODYOVERLAP_Px_Py_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Py_vrr = PAY*I_TWOBODYOVERLAP_Py_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_Pz_vrr = PAX*I_TWOBODYOVERLAP_Px_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Pz_vrr = PAY*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_Py_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_D2x_vrr = PAX*I_TWOBODYOVERLAP_Px_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2x_vrr = PAY*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2x_vrr = PAY*I_TWOBODYOVERLAP_Py_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Py_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Px_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Py_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2y_vrr = PAX*I_TWOBODYOVERLAP_Px_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2y_vrr = PAY*I_TWOBODYOVERLAP_Px_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2y_vrr = PAY*I_TWOBODYOVERLAP_Py_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Px_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Px_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Py_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2z_vrr = PAX*I_TWOBODYOVERLAP_Px_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2z_vrr = PAY*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2z_vrr = PAY*I_TWOBODYOVERLAP_Py_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_F3x_vrr = PBX*I_TWOBODYOVERLAP_Px_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_Py_F3x_vrr = PBX*I_TWOBODYOVERLAP_Py_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_Pz_F3x_vrr = PBX*I_TWOBODYOVERLAP_Pz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_Px_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_Py_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Py_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Pz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Px_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_Py_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_Pz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Pz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Px_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Px_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Py_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_TWOBODYOVERLAP_Pz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Px_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Py_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Pz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Pz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Px_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Px_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Py_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_Pz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Pz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Px_F3y_vrr = PBY*I_TWOBODYOVERLAP_Px_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_Py_F3y_vrr = PBY*I_TWOBODYOVERLAP_Py_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_TWOBODYOVERLAP_Pz_F3y_vrr = PBY*I_TWOBODYOVERLAP_Pz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_TWOBODYOVERLAP_Px_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_Py_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_TWOBODYOVERLAP_Pz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Pz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Px_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_Py_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Py_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Pz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Pz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Px_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Px_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_Py_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Py_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_TWOBODYOVERLAP_Pz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Pz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_Px_vrr = PAX*I_TWOBODYOVERLAP_D2x_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Px_vrr = PAY*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_Px_vrr = PAY*I_TWOBODYOVERLAP_D2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Px_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_Py_vrr = PAX*I_TWOBODYOVERLAP_D2x_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Py_vrr = PAY*I_TWOBODYOVERLAP_D2x_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Py_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_Py_vrr = PAY*I_TWOBODYOVERLAP_D2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_F3z_Py_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_TWOBODYOVERLAP_F3x_Pz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Pz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Pz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 20 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_F3x_vrr = PBX*I_TWOBODYOVERLAP_D2x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F3x_vrr = PBX*I_TWOBODYOVERLAP_Dxy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3x_vrr = PBX*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2xy_vrr = PBY*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Dxy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Dxy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Dxy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2x_F3y_vrr = PBY*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F3y_vrr = PBY*I_TWOBODYOVERLAP_Dxy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3y_vrr = PBY*I_TWOBODYOVERLAP_D2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Dxy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2x_F3z_vrr = PBZ*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_D2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 24 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2y_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2y_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2y_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2z_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2z_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2z_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_F3x_vrr = PAX*I_TWOBODYOVERLAP_D2x_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F3x_vrr = PAY*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F3x_vrr = PAX*I_TWOBODYOVERLAP_D2y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F3x_vrr = PAX*I_TWOBODYOVERLAP_D2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3y_F3x_vrr = PAY*I_TWOBODYOVERLAP_D2y_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F3x_vrr = PAY*I_TWOBODYOVERLAP_D2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3x_F2xy_vrr = PAX*I_TWOBODYOVERLAP_D2x_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F2xy_vrr = PAY*I_TWOBODYOVERLAP_D2x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F2xy_vrr = PAX*I_TWOBODYOVERLAP_D2y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_D2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_D2y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_D2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_F3x_F2xz_vrr = PAX*I_TWOBODYOVERLAP_D2x_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F2xz_vrr = PAY*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F2xz_vrr = PAX*I_TWOBODYOVERLAP_D2y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_D2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F3y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_D2y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F2xz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_D2z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_F3z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_D2x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_D2x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_D2y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_D2z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_D2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_D2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_D2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_D2x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_D2y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_D2z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_D2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F3x_F3y_vrr = PAX*I_TWOBODYOVERLAP_D2x_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F3y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F3y_vrr = PAY*I_TWOBODYOVERLAP_D2x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F3y_vrr = PAX*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F3y_vrr = PAX*I_TWOBODYOVERLAP_D2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3y_F3y_vrr = PAY*I_TWOBODYOVERLAP_D2y_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F3y_vrr = PAY*I_TWOBODYOVERLAP_D2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3x_F2yz_vrr = PAX*I_TWOBODYOVERLAP_D2x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F2yz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F2yz_vrr = PAY*I_TWOBODYOVERLAP_D2x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F2yz_vrr = PAX*I_TWOBODYOVERLAP_D2y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_D2z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_F3y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_D2y_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_D2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F3z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_D2y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_D2z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_D2z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F3x_F3z_vrr = PAX*I_TWOBODYOVERLAP_D2x_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F3z_vrr = PAY*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F3z_vrr = PAX*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Dxy_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F3z_vrr = PAX*I_TWOBODYOVERLAP_D2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_F3y_F3z_vrr = PAY*I_TWOBODYOVERLAP_D2y_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F3z_vrr = PAY*I_TWOBODYOVERLAP_D2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_F3z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 18 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_D2x_vrr = PAX*I_TWOBODYOVERLAP_F3x_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2x_vrr = PAY*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2x_vrr = PAY*I_TWOBODYOVERLAP_F2xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2x_vrr = PAX*I_TWOBODYOVERLAP_F3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2x_vrr = PAX*I_TWOBODYOVERLAP_F3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2x_vrr = PAY*I_TWOBODYOVERLAP_F3y_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2x_vrr = PAY*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_F3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_F3x_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_F3x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_F3y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_F3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_F3y_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_F3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_F3x_Dxz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_F3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_F3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_F3y_Dxz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Dxz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_D2y_vrr = PAX*I_TWOBODYOVERLAP_F3x_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2y_vrr = PAY*I_TWOBODYOVERLAP_F3x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2y_vrr = PAY*I_TWOBODYOVERLAP_F2xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2y_vrr = PAX*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2y_vrr = PAX*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2y_vrr = PAY*I_TWOBODYOVERLAP_F3y_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2y_vrr = PAY*I_TWOBODYOVERLAP_F3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_F3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_F3x_Dyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_F3x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_F3y_Dyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_F3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Dyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4x_D2z_vrr = PAX*I_TWOBODYOVERLAP_F3x_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2z_vrr = PAY*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2z_vrr = PAX*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2z_vrr = PAX*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2z_vrr = PAY*I_TWOBODYOVERLAP_F3y_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2z_vrr = PAY*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_F3x_vrr = PAX*I_TWOBODYOVERLAP_F3x_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3x_vrr = PAY*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3x_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3x_vrr = PAX*I_TWOBODYOVERLAP_F3y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3x_vrr = PAX*I_TWOBODYOVERLAP_F3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3x_vrr = PAY*I_TWOBODYOVERLAP_F3y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3x_vrr = PAY*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2xy_vrr = PAX*I_TWOBODYOVERLAP_F3x_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2xy_vrr = PAY*I_TWOBODYOVERLAP_F3x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2xy_vrr = PAX*I_TWOBODYOVERLAP_F3y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_F3z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_F3y_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_F3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2xz_vrr = PAX*I_TWOBODYOVERLAP_F3x_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2xz_vrr = PAY*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2xz_vrr = PAX*I_TWOBODYOVERLAP_F3y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_F3z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_F3y_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_F3x_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_F3x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_F3y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_F3z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_F3y_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_F3z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_F3x_Fxyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_F3x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_F3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_F3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_F3y_Fxyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_F3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Fxyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_F3x_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_F3y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_F3z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_F3y_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G4x_F3y_vrr = PAX*I_TWOBODYOVERLAP_F3x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3y_vrr = PAY*I_TWOBODYOVERLAP_F3x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3y_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3y_vrr = PAX*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3y_vrr = PAX*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3y_vrr = PAY*I_TWOBODYOVERLAP_F3y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3y_vrr = PAY*I_TWOBODYOVERLAP_F3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2yz_vrr = PAX*I_TWOBODYOVERLAP_F3x_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2yz_vrr = PAY*I_TWOBODYOVERLAP_F3x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2yz_vrr = PAX*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_F3y_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_F3z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_F3x_Fy2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_F3x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_F3y_Fy2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_F3z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Fy2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G4x_F3z_vrr = PAX*I_TWOBODYOVERLAP_F3x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3z_vrr = PAY*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3z_vrr = PAX*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3z_vrr = PAX*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3z_vrr = PAY*I_TWOBODYOVERLAP_F3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3z_vrr = PAY*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_F3x_vrr = PAX*I_TWOBODYOVERLAP_G4x_F3x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F3x_vrr = PAY*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F3x_vrr = PAY*I_TWOBODYOVERLAP_G3xy_F3x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F3x_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_F3x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F3x_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F3x_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F3x_vrr = PAX*I_TWOBODYOVERLAP_G4y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F3x_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F3x_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F3x_vrr = PAX*I_TWOBODYOVERLAP_G4z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5y_F3x_vrr = PAY*I_TWOBODYOVERLAP_G4y_F3x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F3x_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F3x_vrr = PAY*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_G4z_F3x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5x_F2xy_vrr = PAX*I_TWOBODYOVERLAP_G4x_F2xy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F2xy_vrr = PAY*I_TWOBODYOVERLAP_G4x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_G3xy_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F2xy_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F2xy_vrr = PAX*I_TWOBODYOVERLAP_G4y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_G4z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H5y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_G4y_F2xy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_G4z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_G4z_F2xy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H5x_F2xz_vrr = PAX*I_TWOBODYOVERLAP_G4x_F2xz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F2xz_vrr = PAY*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F2xz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F2xz_vrr = PAX*I_TWOBODYOVERLAP_G4y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_G4z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H5y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_G4y_F2xz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H5z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_F2xz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_G4x_Fx2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_G4x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_G4y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_G4z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_G4y_Fx2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_G4z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Fx2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_G4x_Fxyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_G4x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_G4y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_G4z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_G4y_Fxyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_G4z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Fxyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_G4x_Fx2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_G4y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_G4z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_G4y_Fx2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Fx2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H5x_F3y_vrr = PAX*I_TWOBODYOVERLAP_G4x_F3y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F3y_vrr = PAY*I_TWOBODYOVERLAP_G4x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F3y_vrr = PAY*I_TWOBODYOVERLAP_G3xy_F3y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F3y_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_F3y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F3y_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F3y_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F3y_vrr = PAX*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F3y_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F3y_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F3y_vrr = PAX*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5y_F3y_vrr = PAY*I_TWOBODYOVERLAP_G4y_F3y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F3y_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Gy3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F3y_vrr = PAY*I_TWOBODYOVERLAP_G4z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_G4z_F3y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5x_F2yz_vrr = PAX*I_TWOBODYOVERLAP_G4x_F2yz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F2yz_vrr = PAY*I_TWOBODYOVERLAP_G4x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F2yz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F2yz_vrr = PAX*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H5y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_G4y_F2yz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_G4z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H5z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_F2yz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_G4x_Fy2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_G4x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_G4y_Fy2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_G4z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Fy2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H5x_F3z_vrr = PAX*I_TWOBODYOVERLAP_G4x_F3z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F3z_vrr = PAY*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F3z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_F3z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G3xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G3xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F3z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_F3z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F3z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F3z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F3z_vrr = PAX*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F3z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F3z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F3z_vrr = PAX*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5y_F3z_vrr = PAY*I_TWOBODYOVERLAP_G4y_F3z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G3yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F3z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F3z_vrr = PAY*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_F3z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_S_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_KINETIC_S_Px_vrr = PBX*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_S_Py_vrr = PBY*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_S_Pz_vrr = PBZ*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_Px_S_vrr = PAX*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_Py_S_vrr = PAY*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_Pz_S_vrr = PAZ*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_KINETIC_Px_Px_vrr = PBX*I_KINETIC_Px_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_KINETIC_Py_Px_vrr = PBX*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_KINETIC_Pz_Px_vrr = PBX*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_KINETIC_Px_Py_vrr = PBY*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_KINETIC_Py_Py_vrr = PBY*I_KINETIC_Py_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_KINETIC_Pz_Py_vrr = PBY*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_KINETIC_Px_Pz_vrr = PBZ*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_KINETIC_Py_Pz_vrr = PBZ*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_KINETIC_Pz_Pz_vrr = PBZ*I_KINETIC_Pz_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_S_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_KINETIC_S_D2x_vrr = PBX*I_KINETIC_S_Px_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_D2x_vrr-adz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_S_Dxy_vrr = PBY*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_S_Dxz_vrr = PBZ*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_S_D2y_vrr = PBY*I_KINETIC_S_Py_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_D2y_vrr-adz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_S_Dyz_vrr = PBZ*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_S_D2z_vrr = PBZ*I_KINETIC_S_Pz_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_D2z_vrr-adz*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_KINETIC_D2x_S_vrr = PAX*I_KINETIC_Px_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_D2y_S_vrr = PAY*I_KINETIC_Py_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_D2z_S_vrr = PAZ*I_KINETIC_Pz_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_D
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     ************************************************************/
    Double I_KINETIC_Px_D2x_vrr = PAX*I_KINETIC_S_D2x_vrr+2*oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_KINETIC_Py_D2x_vrr = PAY*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_KINETIC_Pz_D2x_vrr = PAZ*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_KINETIC_Px_Dxy_vrr = PAX*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_KINETIC_Py_Dxy_vrr = PAY*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_KINETIC_Pz_Dxy_vrr = PAZ*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_KINETIC_Px_Dxz_vrr = PAX*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_KINETIC_Py_Dxz_vrr = PAY*I_KINETIC_S_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_KINETIC_Pz_Dxz_vrr = PAZ*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Dxz_vrr;
    Double I_KINETIC_Px_D2y_vrr = PAX*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_KINETIC_Py_D2y_vrr = PAY*I_KINETIC_S_D2y_vrr+2*oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_KINETIC_Pz_D2y_vrr = PAZ*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_KINETIC_Px_Dyz_vrr = PAX*I_KINETIC_S_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_KINETIC_Py_Dyz_vrr = PAY*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Py_Dyz_vrr;
    Double I_KINETIC_Pz_Dyz_vrr = PAZ*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Dyz_vrr;
    Double I_KINETIC_Px_D2z_vrr = PAX*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_KINETIC_Py_D2z_vrr = PAY*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_KINETIC_Pz_D2z_vrr = PAZ*I_KINETIC_S_D2z_vrr+2*oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Pz_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_KINETIC_D2x_Px_vrr = PAX*I_KINETIC_Px_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_Dxy_Px_vrr = PAY*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_KINETIC_D2y_Px_vrr = PAY*I_KINETIC_Py_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_D2z_Px_vrr = PAZ*I_KINETIC_Pz_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_D2x_Py_vrr = PAX*I_KINETIC_Px_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_Dxy_Py_vrr = PAY*I_KINETIC_Px_Py_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_KINETIC_D2y_Py_vrr = PAY*I_KINETIC_Py_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_D2z_Py_vrr = PAZ*I_KINETIC_Pz_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_D2x_Pz_vrr = PAX*I_KINETIC_Px_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_KINETIC_Dxy_Pz_vrr = PAY*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_KINETIC_D2y_Pz_vrr = PAY*I_KINETIC_Py_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_KINETIC_D2z_Pz_vrr = PAZ*I_KINETIC_Pz_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_S_D
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     ************************************************************/
    Double I_KINETIC_D2x_D2x_vrr = PAX*I_KINETIC_Px_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+2*oned2z*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_Dxy_D2x_vrr = PAY*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_KINETIC_D2y_D2x_vrr = PAY*I_KINETIC_Py_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_D2z_D2x_vrr = PAZ*I_KINETIC_Pz_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_D2x_Dxy_vrr = PAX*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_Dxy_Dxy_vrr = PAY*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Dxy_vrr;
    Double I_KINETIC_D2y_Dxy_vrr = PAY*I_KINETIC_Py_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2z_Dxy_vrr = PAZ*I_KINETIC_Pz_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2x_Dxz_vrr = PAX*I_KINETIC_Px_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_Dxy_Dxz_vrr = PAY*I_KINETIC_Px_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Dxz_vrr;
    Double I_KINETIC_D2y_Dxz_vrr = PAY*I_KINETIC_Py_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_D2z_Dxz_vrr = PAZ*I_KINETIC_Pz_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_D2x_D2y_vrr = PAX*I_KINETIC_Px_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_Dxy_D2y_vrr = PAY*I_KINETIC_Px_D2y_vrr+2*oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_KINETIC_D2y_D2y_vrr = PAY*I_KINETIC_Py_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+2*oned2z*I_KINETIC_Py_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_D2z_D2y_vrr = PAZ*I_KINETIC_Pz_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_D2x_Dyz_vrr = PAX*I_KINETIC_Px_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_Dxy_Dyz_vrr = PAY*I_KINETIC_Px_Dyz_vrr+oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Dyz_vrr;
    Double I_KINETIC_D2y_Dyz_vrr = PAY*I_KINETIC_Py_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_Py_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2z_Dyz_vrr = PAZ*I_KINETIC_Pz_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_Pz_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2x_D2z_vrr = PAX*I_KINETIC_Px_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_KINETIC_Dxy_D2z_vrr = PAY*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2z_vrr;
    Double I_KINETIC_D2y_D2z_vrr = PAY*I_KINETIC_Py_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_KINETIC_D2z_D2z_vrr = PAZ*I_KINETIC_Pz_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+2*oned2z*I_KINETIC_Pz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_S_D
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_KINETIC_Px_F3x_vrr = PBX*I_KINETIC_Px_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+2*oned2z*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Px_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_KINETIC_Py_F3x_vrr = PBX*I_KINETIC_Py_D2x_vrr+2*oned2z*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_Py_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_KINETIC_Pz_F3x_vrr = PBX*I_KINETIC_Pz_D2x_vrr+2*oned2z*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_KINETIC_Px_F2xy_vrr = PBY*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Px_F2xy_vrr;
    Double I_KINETIC_Py_F2xy_vrr = PBY*I_KINETIC_Py_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Py_F2xy_vrr;
    Double I_KINETIC_Pz_F2xy_vrr = PBY*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F2xy_vrr;
    Double I_KINETIC_Px_F2xz_vrr = PBZ*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Px_F2xz_vrr;
    Double I_KINETIC_Py_F2xz_vrr = PBZ*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Py_F2xz_vrr;
    Double I_KINETIC_Pz_F2xz_vrr = PBZ*I_KINETIC_Pz_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F2xz_vrr;
    Double I_KINETIC_Px_Fx2y_vrr = PBX*I_KINETIC_Px_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Px_Fx2y_vrr;
    Double I_KINETIC_Py_Fx2y_vrr = PBX*I_KINETIC_Py_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Py_Fx2y_vrr;
    Double I_KINETIC_Pz_Fx2y_vrr = PBX*I_KINETIC_Pz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Fx2y_vrr;
    Double I_KINETIC_Px_Fxyz_vrr = PBZ*I_KINETIC_Px_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Px_Fxyz_vrr;
    Double I_KINETIC_Py_Fxyz_vrr = PBZ*I_KINETIC_Py_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Py_Fxyz_vrr;
    Double I_KINETIC_Pz_Fxyz_vrr = PBZ*I_KINETIC_Pz_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Fxyz_vrr;
    Double I_KINETIC_Px_Fx2z_vrr = PBX*I_KINETIC_Px_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Px_Fx2z_vrr;
    Double I_KINETIC_Py_Fx2z_vrr = PBX*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Py_Fx2z_vrr;
    Double I_KINETIC_Pz_Fx2z_vrr = PBX*I_KINETIC_Pz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Fx2z_vrr;
    Double I_KINETIC_Px_F3y_vrr = PBY*I_KINETIC_Px_D2y_vrr+2*oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_Px_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_KINETIC_Py_F3y_vrr = PBY*I_KINETIC_Py_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+2*oned2z*I_KINETIC_Py_Py_vrr+twoxi*I_TWOBODYOVERLAP_Py_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_KINETIC_Pz_F3y_vrr = PBY*I_KINETIC_Pz_D2y_vrr+2*oned2z*I_KINETIC_Pz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_KINETIC_Px_F2yz_vrr = PBZ*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Px_F2yz_vrr;
    Double I_KINETIC_Py_F2yz_vrr = PBZ*I_KINETIC_Py_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Py_F2yz_vrr;
    Double I_KINETIC_Pz_F2yz_vrr = PBZ*I_KINETIC_Pz_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F2yz_vrr;
    Double I_KINETIC_Px_Fy2z_vrr = PBY*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Px_Fy2z_vrr;
    Double I_KINETIC_Py_Fy2z_vrr = PBY*I_KINETIC_Py_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Py_Fy2z_vrr;
    Double I_KINETIC_Pz_Fy2z_vrr = PBY*I_KINETIC_Pz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Fy2z_vrr;
    Double I_KINETIC_Px_F3z_vrr = PBZ*I_KINETIC_Px_D2z_vrr+2*oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Px_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_KINETIC_Py_F3z_vrr = PBZ*I_KINETIC_Py_D2z_vrr+2*oned2z*I_KINETIC_Py_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Py_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_KINETIC_Pz_F3z_vrr = PBZ*I_KINETIC_Pz_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+2*oned2z*I_KINETIC_Pz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Pz_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_KINETIC_F3x_Px_vrr = PAX*I_KINETIC_D2x_Px_vrr+2*oned2z*I_KINETIC_Px_Px_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Px_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_KINETIC_F2xy_Px_vrr = PAY*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_KINETIC_F2xz_Px_vrr = PAZ*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_KINETIC_F3y_Px_vrr = PAY*I_KINETIC_D2y_Px_vrr+2*oned2z*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Px_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_KINETIC_F2yz_Px_vrr = PAZ*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_KINETIC_F3z_Px_vrr = PAZ*I_KINETIC_D2z_Px_vrr+2*oned2z*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Px_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_KINETIC_F3x_Py_vrr = PAX*I_KINETIC_D2x_Py_vrr+2*oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Py_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_KINETIC_F2xy_Py_vrr = PAY*I_KINETIC_D2x_Py_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_KINETIC_F2xz_Py_vrr = PAZ*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_KINETIC_F3y_Py_vrr = PAY*I_KINETIC_D2y_Py_vrr+2*oned2z*I_KINETIC_Py_Py_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Py_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_KINETIC_F2yz_Py_vrr = PAZ*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Py_vrr;
    Double I_KINETIC_F3z_Py_vrr = PAZ*I_KINETIC_D2z_Py_vrr+2*oned2z*I_KINETIC_Pz_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Py_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_KINETIC_F3x_Pz_vrr = PAX*I_KINETIC_D2x_Pz_vrr+2*oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Pz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_KINETIC_F2xy_Pz_vrr = PAY*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Pz_vrr;
    Double I_KINETIC_F2xz_Pz_vrr = PAZ*I_KINETIC_D2x_Pz_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Pz_vrr;
    Double I_KINETIC_F3y_Pz_vrr = PAY*I_KINETIC_D2y_Pz_vrr+2*oned2z*I_KINETIC_Py_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Pz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_KINETIC_F2yz_Pz_vrr = PAZ*I_KINETIC_D2y_Pz_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Pz_vrr;
    Double I_KINETIC_F3z_Pz_vrr = PAZ*I_KINETIC_D2z_Pz_vrr+2*oned2z*I_KINETIC_Pz_Pz_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Pz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 20 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_KINETIC_D2x_F3x_vrr = PBX*I_KINETIC_D2x_D2x_vrr+2*oned2z*I_KINETIC_Px_D2x_vrr+2*oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3x_vrr-2*adz*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_Dxy_F3x_vrr = PBX*I_KINETIC_Dxy_D2x_vrr+oned2z*I_KINETIC_Py_D2x_vrr+2*oned2z*I_KINETIC_Dxy_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_KINETIC_D2y_F3x_vrr = PBX*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_D2z_F3x_vrr = PBX*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_KINETIC_D2x_F2xy_vrr = PBY*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_KINETIC_Dxy_F2xy_vrr = PBY*I_KINETIC_Dxy_D2x_vrr+oned2z*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F2xy_vrr;
    Double I_KINETIC_D2y_F2xy_vrr = PBY*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_KINETIC_D2z_F2xy_vrr = PBY*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2xy_vrr;
    Double I_KINETIC_D2x_F2xz_vrr = PBZ*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_KINETIC_Dxy_F2xz_vrr = PBZ*I_KINETIC_Dxy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F2xz_vrr;
    Double I_KINETIC_D2y_F2xz_vrr = PBZ*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2xz_vrr;
    Double I_KINETIC_D2z_F2xz_vrr = PBZ*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2xz_vrr;
    Double I_KINETIC_D2x_Fx2y_vrr = PBX*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_KINETIC_Dxy_Fx2y_vrr = PBX*I_KINETIC_Dxy_D2y_vrr+oned2z*I_KINETIC_Py_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fx2y_vrr;
    Double I_KINETIC_D2y_Fx2y_vrr = PBX*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_KINETIC_D2z_Fx2y_vrr = PBX*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fx2y_vrr;
    Double I_KINETIC_D2x_Fxyz_vrr = PBZ*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fxyz_vrr;
    Double I_KINETIC_Dxy_Fxyz_vrr = PBZ*I_KINETIC_Dxy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fxyz_vrr;
    Double I_KINETIC_D2y_Fxyz_vrr = PBZ*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fxyz_vrr;
    Double I_KINETIC_D2z_Fxyz_vrr = PBZ*I_KINETIC_D2z_Dxy_vrr+2*oned2z*I_KINETIC_Pz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fxyz_vrr;
    Double I_KINETIC_D2x_Fx2z_vrr = PBX*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_KINETIC_Dxy_Fx2z_vrr = PBX*I_KINETIC_Dxy_D2z_vrr+oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fx2z_vrr;
    Double I_KINETIC_D2y_Fx2z_vrr = PBX*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_KINETIC_D2z_Fx2z_vrr = PBX*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fx2z_vrr;
    Double I_KINETIC_D2x_F3y_vrr = PBY*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3y_vrr-2*adz*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_Dxy_F3y_vrr = PBY*I_KINETIC_Dxy_D2y_vrr+oned2z*I_KINETIC_Px_D2y_vrr+2*oned2z*I_KINETIC_Dxy_Py_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_KINETIC_D2y_F3y_vrr = PBY*I_KINETIC_D2y_D2y_vrr+2*oned2z*I_KINETIC_Py_D2y_vrr+2*oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_D2z_F3y_vrr = PBY*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_KINETIC_D2x_F2yz_vrr = PBZ*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_KINETIC_Dxy_F2yz_vrr = PBZ*I_KINETIC_Dxy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F2yz_vrr;
    Double I_KINETIC_D2y_F2yz_vrr = PBZ*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2yz_vrr;
    Double I_KINETIC_D2z_F2yz_vrr = PBZ*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_Pz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2yz_vrr;
    Double I_KINETIC_D2x_Fy2z_vrr = PBY*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_KINETIC_Dxy_Fy2z_vrr = PBY*I_KINETIC_Dxy_D2z_vrr+oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fy2z_vrr;
    Double I_KINETIC_D2y_Fy2z_vrr = PBY*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fy2z_vrr;
    Double I_KINETIC_D2z_Fy2z_vrr = PBY*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fy2z_vrr;
    Double I_KINETIC_D2x_F3z_vrr = PBZ*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3z_vrr-2*adz*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_Dxy_F3z_vrr = PBZ*I_KINETIC_Dxy_D2z_vrr+2*oned2z*I_KINETIC_Dxy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_KINETIC_D2y_F3z_vrr = PBZ*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_D2z_F3z_vrr = PBZ*I_KINETIC_D2z_D2z_vrr+2*oned2z*I_KINETIC_Pz_D2z_vrr+2*oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 24 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     ************************************************************/
    Double I_KINETIC_F3x_D2x_vrr = PAX*I_KINETIC_D2x_D2x_vrr+2*oned2z*I_KINETIC_Px_D2x_vrr+2*oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_KINETIC_F2xy_D2x_vrr = PAY*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_KINETIC_F2xz_D2x_vrr = PAZ*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_KINETIC_F3y_D2x_vrr = PAY*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_KINETIC_F2yz_D2x_vrr = PAZ*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_KINETIC_F3z_D2x_vrr = PAZ*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_KINETIC_F3x_Dxy_vrr = PAX*I_KINETIC_D2x_Dxy_vrr+2*oned2z*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_KINETIC_F2xy_Dxy_vrr = PAY*I_KINETIC_D2x_Dxy_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_KINETIC_F2xz_Dxy_vrr = PAZ*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_KINETIC_F3y_Dxy_vrr = PAY*I_KINETIC_D2y_Dxy_vrr+2*oned2z*I_KINETIC_Py_Dxy_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_KINETIC_F2yz_Dxy_vrr = PAZ*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_KINETIC_F3z_Dxy_vrr = PAZ*I_KINETIC_D2z_Dxy_vrr+2*oned2z*I_KINETIC_Pz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_KINETIC_F3x_Dxz_vrr = PAX*I_KINETIC_D2x_Dxz_vrr+2*oned2z*I_KINETIC_Px_Dxz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_KINETIC_F2xy_Dxz_vrr = PAY*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxz_vrr;
    Double I_KINETIC_F2xz_Dxz_vrr = PAZ*I_KINETIC_D2x_Dxz_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxz_vrr;
    Double I_KINETIC_F3y_Dxz_vrr = PAY*I_KINETIC_D2y_Dxz_vrr+2*oned2z*I_KINETIC_Py_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_KINETIC_F2yz_Dxz_vrr = PAZ*I_KINETIC_D2y_Dxz_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxz_vrr;
    Double I_KINETIC_F3z_Dxz_vrr = PAZ*I_KINETIC_D2z_Dxz_vrr+2*oned2z*I_KINETIC_Pz_Dxz_vrr+oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Dxz_vrr;
    Double I_KINETIC_F3x_D2y_vrr = PAX*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_KINETIC_F2xy_D2y_vrr = PAY*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_KINETIC_F2xz_D2y_vrr = PAZ*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_KINETIC_F3y_D2y_vrr = PAY*I_KINETIC_D2y_D2y_vrr+2*oned2z*I_KINETIC_Py_D2y_vrr+2*oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_KINETIC_F2yz_D2y_vrr = PAZ*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_KINETIC_F3z_D2y_vrr = PAZ*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_Pz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_KINETIC_F3x_Dyz_vrr = PAX*I_KINETIC_D2x_Dyz_vrr+2*oned2z*I_KINETIC_Px_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_KINETIC_F2xy_Dyz_vrr = PAY*I_KINETIC_D2x_Dyz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dyz_vrr;
    Double I_KINETIC_F2xz_Dyz_vrr = PAZ*I_KINETIC_D2x_Dyz_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dyz_vrr;
    Double I_KINETIC_F3y_Dyz_vrr = PAY*I_KINETIC_D2y_Dyz_vrr+2*oned2z*I_KINETIC_Py_Dyz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dyz_vrr;
    Double I_KINETIC_F2yz_Dyz_vrr = PAZ*I_KINETIC_D2y_Dyz_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dyz_vrr;
    Double I_KINETIC_F3z_Dyz_vrr = PAZ*I_KINETIC_D2z_Dyz_vrr+2*oned2z*I_KINETIC_Pz_Dyz_vrr+oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Dyz_vrr;
    Double I_KINETIC_F3x_D2z_vrr = PAX*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_KINETIC_F2xy_D2z_vrr = PAY*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_KINETIC_F2xz_D2z_vrr = PAZ*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_KINETIC_F3y_D2z_vrr = PAY*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_KINETIC_F2yz_D2z_vrr = PAZ*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_KINETIC_F3z_D2z_vrr = PAZ*I_KINETIC_D2z_D2z_vrr+2*oned2z*I_KINETIC_Pz_D2z_vrr+2*oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_F
     * RHS shell quartet name: SQ_KINETIC_P_F
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_F
     ************************************************************/
    Double I_KINETIC_F3x_F3x_vrr = PAX*I_KINETIC_D2x_F3x_vrr+2*oned2z*I_KINETIC_Px_F3x_vrr+3*oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F3x_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F3x_vrr;
    Double I_KINETIC_F2xy_F3x_vrr = PAY*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F3x_vrr;
    Double I_KINETIC_F2xz_F3x_vrr = PAZ*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F3x_vrr;
    Double I_KINETIC_Fx2y_F3x_vrr = PAX*I_KINETIC_D2y_F3x_vrr+3*oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F3x_vrr;
    Double I_KINETIC_Fxyz_F3x_vrr = PAZ*I_KINETIC_Dxy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F3x_vrr;
    Double I_KINETIC_Fx2z_F3x_vrr = PAX*I_KINETIC_D2z_F3x_vrr+3*oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F3x_vrr;
    Double I_KINETIC_F3y_F3x_vrr = PAY*I_KINETIC_D2y_F3x_vrr+2*oned2z*I_KINETIC_Py_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F3x_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F3x_vrr;
    Double I_KINETIC_F2yz_F3x_vrr = PAZ*I_KINETIC_D2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F3x_vrr;
    Double I_KINETIC_Fy2z_F3x_vrr = PAY*I_KINETIC_D2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F3x_vrr;
    Double I_KINETIC_F3z_F3x_vrr = PAZ*I_KINETIC_D2z_F3x_vrr+2*oned2z*I_KINETIC_Pz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F3x_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F3x_vrr;
    Double I_KINETIC_F3x_F2xy_vrr = PAX*I_KINETIC_D2x_F2xy_vrr+2*oned2z*I_KINETIC_Px_F2xy_vrr+2*oned2z*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F2xy_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F2xy_vrr;
    Double I_KINETIC_F2xy_F2xy_vrr = PAY*I_KINETIC_D2x_F2xy_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F2xy_vrr;
    Double I_KINETIC_F2xz_F2xy_vrr = PAZ*I_KINETIC_D2x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F2xy_vrr;
    Double I_KINETIC_Fx2y_F2xy_vrr = PAX*I_KINETIC_D2y_F2xy_vrr+2*oned2z*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F2xy_vrr;
    Double I_KINETIC_Fxyz_F2xy_vrr = PAZ*I_KINETIC_Dxy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F2xy_vrr;
    Double I_KINETIC_Fx2z_F2xy_vrr = PAX*I_KINETIC_D2z_F2xy_vrr+2*oned2z*I_KINETIC_D2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F2xy_vrr;
    Double I_KINETIC_F3y_F2xy_vrr = PAY*I_KINETIC_D2y_F2xy_vrr+2*oned2z*I_KINETIC_Py_F2xy_vrr+oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F2xy_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F2xy_vrr;
    Double I_KINETIC_F2yz_F2xy_vrr = PAZ*I_KINETIC_D2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F2xy_vrr;
    Double I_KINETIC_Fy2z_F2xy_vrr = PAY*I_KINETIC_D2z_F2xy_vrr+oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F2xy_vrr;
    Double I_KINETIC_F3z_F2xy_vrr = PAZ*I_KINETIC_D2z_F2xy_vrr+2*oned2z*I_KINETIC_Pz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F2xy_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F2xy_vrr;
    Double I_KINETIC_F3x_F2xz_vrr = PAX*I_KINETIC_D2x_F2xz_vrr+2*oned2z*I_KINETIC_Px_F2xz_vrr+2*oned2z*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F2xz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F2xz_vrr;
    Double I_KINETIC_F2xy_F2xz_vrr = PAY*I_KINETIC_D2x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F2xz_vrr;
    Double I_KINETIC_F2xz_F2xz_vrr = PAZ*I_KINETIC_D2x_F2xz_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F2xz_vrr;
    Double I_KINETIC_Fx2y_F2xz_vrr = PAX*I_KINETIC_D2y_F2xz_vrr+2*oned2z*I_KINETIC_D2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F2xz_vrr;
    Double I_KINETIC_Fxyz_F2xz_vrr = PAZ*I_KINETIC_Dxy_F2xz_vrr+oned2z*I_KINETIC_Dxy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F2xz_vrr;
    Double I_KINETIC_Fx2z_F2xz_vrr = PAX*I_KINETIC_D2z_F2xz_vrr+2*oned2z*I_KINETIC_D2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F2xz_vrr;
    Double I_KINETIC_F3y_F2xz_vrr = PAY*I_KINETIC_D2y_F2xz_vrr+2*oned2z*I_KINETIC_Py_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F2xz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F2xz_vrr;
    Double I_KINETIC_F2yz_F2xz_vrr = PAZ*I_KINETIC_D2y_F2xz_vrr+oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F2xz_vrr;
    Double I_KINETIC_Fy2z_F2xz_vrr = PAY*I_KINETIC_D2z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F2xz_vrr;
    Double I_KINETIC_F3z_F2xz_vrr = PAZ*I_KINETIC_D2z_F2xz_vrr+2*oned2z*I_KINETIC_Pz_F2xz_vrr+oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F2xz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F2xz_vrr;
    Double I_KINETIC_F3x_Fx2y_vrr = PAX*I_KINETIC_D2x_Fx2y_vrr+2*oned2z*I_KINETIC_Px_Fx2y_vrr+oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fx2y_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Fx2y_vrr;
    Double I_KINETIC_F2xy_Fx2y_vrr = PAY*I_KINETIC_D2x_Fx2y_vrr+2*oned2z*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fx2y_vrr;
    Double I_KINETIC_F2xz_Fx2y_vrr = PAZ*I_KINETIC_D2x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fx2y_vrr;
    Double I_KINETIC_Fx2y_Fx2y_vrr = PAX*I_KINETIC_D2y_Fx2y_vrr+oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fx2y_vrr;
    Double I_KINETIC_Fxyz_Fx2y_vrr = PAZ*I_KINETIC_Dxy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Fx2y_vrr;
    Double I_KINETIC_Fx2z_Fx2y_vrr = PAX*I_KINETIC_D2z_Fx2y_vrr+oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fx2y_vrr;
    Double I_KINETIC_F3y_Fx2y_vrr = PAY*I_KINETIC_D2y_Fx2y_vrr+2*oned2z*I_KINETIC_Py_Fx2y_vrr+2*oned2z*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fx2y_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fx2y_vrr;
    Double I_KINETIC_F2yz_Fx2y_vrr = PAZ*I_KINETIC_D2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fx2y_vrr;
    Double I_KINETIC_Fy2z_Fx2y_vrr = PAY*I_KINETIC_D2z_Fx2y_vrr+2*oned2z*I_KINETIC_D2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Fx2y_vrr;
    Double I_KINETIC_F3z_Fx2y_vrr = PAZ*I_KINETIC_D2z_Fx2y_vrr+2*oned2z*I_KINETIC_Pz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fx2y_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Fx2y_vrr;
    Double I_KINETIC_F3x_Fxyz_vrr = PAX*I_KINETIC_D2x_Fxyz_vrr+2*oned2z*I_KINETIC_Px_Fxyz_vrr+oned2z*I_KINETIC_D2x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Fxyz_vrr;
    Double I_KINETIC_F2xy_Fxyz_vrr = PAY*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fxyz_vrr;
    Double I_KINETIC_F2xz_Fxyz_vrr = PAZ*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fxyz_vrr;
    Double I_KINETIC_Fx2y_Fxyz_vrr = PAX*I_KINETIC_D2y_Fxyz_vrr+oned2z*I_KINETIC_D2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fxyz_vrr;
    Double I_KINETIC_Fxyz_Fxyz_vrr = PAZ*I_KINETIC_Dxy_Fxyz_vrr+oned2z*I_KINETIC_Dxy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Fxyz_vrr;
    Double I_KINETIC_Fx2z_Fxyz_vrr = PAX*I_KINETIC_D2z_Fxyz_vrr+oned2z*I_KINETIC_D2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fxyz_vrr;
    Double I_KINETIC_F3y_Fxyz_vrr = PAY*I_KINETIC_D2y_Fxyz_vrr+2*oned2z*I_KINETIC_Py_Fxyz_vrr+oned2z*I_KINETIC_D2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fxyz_vrr;
    Double I_KINETIC_F2yz_Fxyz_vrr = PAZ*I_KINETIC_D2y_Fxyz_vrr+oned2z*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fxyz_vrr;
    Double I_KINETIC_Fy2z_Fxyz_vrr = PAY*I_KINETIC_D2z_Fxyz_vrr+oned2z*I_KINETIC_D2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Fxyz_vrr;
    Double I_KINETIC_F3z_Fxyz_vrr = PAZ*I_KINETIC_D2z_Fxyz_vrr+2*oned2z*I_KINETIC_Pz_Fxyz_vrr+oned2z*I_KINETIC_D2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Fxyz_vrr;
    Double I_KINETIC_F3x_Fx2z_vrr = PAX*I_KINETIC_D2x_Fx2z_vrr+2*oned2z*I_KINETIC_Px_Fx2z_vrr+oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fx2z_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Fx2z_vrr;
    Double I_KINETIC_F2xy_Fx2z_vrr = PAY*I_KINETIC_D2x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fx2z_vrr;
    Double I_KINETIC_F2xz_Fx2z_vrr = PAZ*I_KINETIC_D2x_Fx2z_vrr+2*oned2z*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fx2z_vrr;
    Double I_KINETIC_Fx2y_Fx2z_vrr = PAX*I_KINETIC_D2y_Fx2z_vrr+oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fx2z_vrr;
    Double I_KINETIC_Fxyz_Fx2z_vrr = PAZ*I_KINETIC_Dxy_Fx2z_vrr+2*oned2z*I_KINETIC_Dxy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Fx2z_vrr;
    Double I_KINETIC_Fx2z_Fx2z_vrr = PAX*I_KINETIC_D2z_Fx2z_vrr+oned2z*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fx2z_vrr;
    Double I_KINETIC_F3y_Fx2z_vrr = PAY*I_KINETIC_D2y_Fx2z_vrr+2*oned2z*I_KINETIC_Py_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fx2z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fx2z_vrr;
    Double I_KINETIC_F2yz_Fx2z_vrr = PAZ*I_KINETIC_D2y_Fx2z_vrr+2*oned2z*I_KINETIC_D2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fx2z_vrr;
    Double I_KINETIC_Fy2z_Fx2z_vrr = PAY*I_KINETIC_D2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Fx2z_vrr;
    Double I_KINETIC_F3z_Fx2z_vrr = PAZ*I_KINETIC_D2z_Fx2z_vrr+2*oned2z*I_KINETIC_Pz_Fx2z_vrr+2*oned2z*I_KINETIC_D2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fx2z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Fx2z_vrr;
    Double I_KINETIC_F3x_F3y_vrr = PAX*I_KINETIC_D2x_F3y_vrr+2*oned2z*I_KINETIC_Px_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F3y_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F3y_vrr;
    Double I_KINETIC_F2xy_F3y_vrr = PAY*I_KINETIC_D2x_F3y_vrr+3*oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F3y_vrr;
    Double I_KINETIC_F2xz_F3y_vrr = PAZ*I_KINETIC_D2x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F3y_vrr;
    Double I_KINETIC_Fx2y_F3y_vrr = PAX*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F3y_vrr;
    Double I_KINETIC_Fxyz_F3y_vrr = PAZ*I_KINETIC_Dxy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F3y_vrr;
    Double I_KINETIC_Fx2z_F3y_vrr = PAX*I_KINETIC_D2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F3y_vrr;
    Double I_KINETIC_F3y_F3y_vrr = PAY*I_KINETIC_D2y_F3y_vrr+2*oned2z*I_KINETIC_Py_F3y_vrr+3*oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F3y_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F3y_vrr;
    Double I_KINETIC_F2yz_F3y_vrr = PAZ*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F3y_vrr;
    Double I_KINETIC_Fy2z_F3y_vrr = PAY*I_KINETIC_D2z_F3y_vrr+3*oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F3y_vrr;
    Double I_KINETIC_F3z_F3y_vrr = PAZ*I_KINETIC_D2z_F3y_vrr+2*oned2z*I_KINETIC_Pz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F3y_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F3y_vrr;
    Double I_KINETIC_F3x_F2yz_vrr = PAX*I_KINETIC_D2x_F2yz_vrr+2*oned2z*I_KINETIC_Px_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F2yz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F2yz_vrr;
    Double I_KINETIC_F2xy_F2yz_vrr = PAY*I_KINETIC_D2x_F2yz_vrr+2*oned2z*I_KINETIC_D2x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F2yz_vrr;
    Double I_KINETIC_F2xz_F2yz_vrr = PAZ*I_KINETIC_D2x_F2yz_vrr+oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F2yz_vrr;
    Double I_KINETIC_Fx2y_F2yz_vrr = PAX*I_KINETIC_D2y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F2yz_vrr;
    Double I_KINETIC_Fxyz_F2yz_vrr = PAZ*I_KINETIC_Dxy_F2yz_vrr+oned2z*I_KINETIC_Dxy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F2yz_vrr;
    Double I_KINETIC_Fx2z_F2yz_vrr = PAX*I_KINETIC_D2z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F2yz_vrr;
    Double I_KINETIC_F3y_F2yz_vrr = PAY*I_KINETIC_D2y_F2yz_vrr+2*oned2z*I_KINETIC_Py_F2yz_vrr+2*oned2z*I_KINETIC_D2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F2yz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F2yz_vrr;
    Double I_KINETIC_F2yz_F2yz_vrr = PAZ*I_KINETIC_D2y_F2yz_vrr+oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F2yz_vrr;
    Double I_KINETIC_Fy2z_F2yz_vrr = PAY*I_KINETIC_D2z_F2yz_vrr+2*oned2z*I_KINETIC_D2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F2yz_vrr;
    Double I_KINETIC_F3z_F2yz_vrr = PAZ*I_KINETIC_D2z_F2yz_vrr+2*oned2z*I_KINETIC_Pz_F2yz_vrr+oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F2yz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F2yz_vrr;
    Double I_KINETIC_F3x_Fy2z_vrr = PAX*I_KINETIC_D2x_Fy2z_vrr+2*oned2z*I_KINETIC_Px_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Fy2z_vrr;
    Double I_KINETIC_F2xy_Fy2z_vrr = PAY*I_KINETIC_D2x_Fy2z_vrr+oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fy2z_vrr;
    Double I_KINETIC_F2xz_Fy2z_vrr = PAZ*I_KINETIC_D2x_Fy2z_vrr+2*oned2z*I_KINETIC_D2x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fy2z_vrr;
    Double I_KINETIC_Fx2y_Fy2z_vrr = PAX*I_KINETIC_D2y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fy2z_vrr;
    Double I_KINETIC_Fxyz_Fy2z_vrr = PAZ*I_KINETIC_Dxy_Fy2z_vrr+2*oned2z*I_KINETIC_Dxy_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Fy2z_vrr;
    Double I_KINETIC_Fx2z_Fy2z_vrr = PAX*I_KINETIC_D2z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fy2z_vrr;
    Double I_KINETIC_F3y_Fy2z_vrr = PAY*I_KINETIC_D2y_Fy2z_vrr+2*oned2z*I_KINETIC_Py_Fy2z_vrr+oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fy2z_vrr;
    Double I_KINETIC_F2yz_Fy2z_vrr = PAZ*I_KINETIC_D2y_Fy2z_vrr+2*oned2z*I_KINETIC_D2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fy2z_vrr;
    Double I_KINETIC_Fy2z_Fy2z_vrr = PAY*I_KINETIC_D2z_Fy2z_vrr+oned2z*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Fy2z_vrr;
    Double I_KINETIC_F3z_Fy2z_vrr = PAZ*I_KINETIC_D2z_Fy2z_vrr+2*oned2z*I_KINETIC_Pz_Fy2z_vrr+2*oned2z*I_KINETIC_D2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Fy2z_vrr;
    Double I_KINETIC_F3x_F3z_vrr = PAX*I_KINETIC_D2x_F3z_vrr+2*oned2z*I_KINETIC_Px_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F3z_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F3z_vrr;
    Double I_KINETIC_F2xy_F3z_vrr = PAY*I_KINETIC_D2x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F3z_vrr;
    Double I_KINETIC_F2xz_F3z_vrr = PAZ*I_KINETIC_D2x_F3z_vrr+3*oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F3z_vrr;
    Double I_KINETIC_Fx2y_F3z_vrr = PAX*I_KINETIC_D2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F3z_vrr;
    Double I_KINETIC_Fxyz_F3z_vrr = PAZ*I_KINETIC_Dxy_F3z_vrr+3*oned2z*I_KINETIC_Dxy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F3z_vrr;
    Double I_KINETIC_Fx2z_F3z_vrr = PAX*I_KINETIC_D2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F3z_vrr;
    Double I_KINETIC_F3y_F3z_vrr = PAY*I_KINETIC_D2y_F3z_vrr+2*oned2z*I_KINETIC_Py_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F3z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F3z_vrr;
    Double I_KINETIC_F2yz_F3z_vrr = PAZ*I_KINETIC_D2y_F3z_vrr+3*oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F3z_vrr;
    Double I_KINETIC_Fy2z_F3z_vrr = PAY*I_KINETIC_D2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F3z_vrr;
    Double I_KINETIC_F3z_F3z_vrr = PAZ*I_KINETIC_D2z_F3z_vrr+2*oned2z*I_KINETIC_Pz_F3z_vrr+3*oned2z*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F3z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 18 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_D
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     ************************************************************/
    Double I_KINETIC_G4x_D2x_vrr = PAX*I_KINETIC_F3x_D2x_vrr+3*oned2z*I_KINETIC_D2x_D2x_vrr+2*oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2x_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_KINETIC_G3xy_D2x_vrr = PAY*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_KINETIC_G3xz_D2x_vrr = PAZ*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2x_vrr;
    Double I_KINETIC_G2x2y_D2x_vrr = PAY*I_KINETIC_F2xy_D2x_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2x_vrr-bdz*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_KINETIC_G2x2z_D2x_vrr = PAZ*I_KINETIC_F2xz_D2x_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_KINETIC_Gx3y_D2x_vrr = PAX*I_KINETIC_F3y_D2x_vrr+2*oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2x_vrr;
    Double I_KINETIC_Gx3z_D2x_vrr = PAX*I_KINETIC_F3z_D2x_vrr+2*oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2x_vrr;
    Double I_KINETIC_G4y_D2x_vrr = PAY*I_KINETIC_F3y_D2x_vrr+3*oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2x_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_KINETIC_G3yz_D2x_vrr = PAZ*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2x_vrr;
    Double I_KINETIC_G2y2z_D2x_vrr = PAZ*I_KINETIC_F2yz_D2x_vrr+oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_KINETIC_Gy3z_D2x_vrr = PAY*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2x_vrr;
    Double I_KINETIC_G4z_D2x_vrr = PAZ*I_KINETIC_F3z_D2x_vrr+3*oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2x_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_KINETIC_G4x_Dxy_vrr = PAX*I_KINETIC_F3x_Dxy_vrr+3*oned2z*I_KINETIC_D2x_Dxy_vrr+oned2z*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Dxy_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_KINETIC_G3xy_Dxy_vrr = PAY*I_KINETIC_F3x_Dxy_vrr+oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Dxy_vrr;
    Double I_KINETIC_G3xz_Dxy_vrr = PAZ*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Dxy_vrr;
    Double I_KINETIC_G2x2y_Dxy_vrr = PAY*I_KINETIC_F2xy_Dxy_vrr+oned2z*I_KINETIC_D2x_Dxy_vrr+oned2z*I_KINETIC_F2xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_KINETIC_G2x2z_Dxy_vrr = PAZ*I_KINETIC_F2xz_Dxy_vrr+oned2z*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_KINETIC_Gx3y_Dxy_vrr = PAX*I_KINETIC_F3y_Dxy_vrr+oned2z*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr;
    Double I_KINETIC_Gx3z_Dxy_vrr = PAX*I_KINETIC_F3z_Dxy_vrr+oned2z*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr;
    Double I_KINETIC_G4y_Dxy_vrr = PAY*I_KINETIC_F3y_Dxy_vrr+3*oned2z*I_KINETIC_D2y_Dxy_vrr+oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dxy_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_KINETIC_G3yz_Dxy_vrr = PAZ*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dxy_vrr;
    Double I_KINETIC_G2y2z_Dxy_vrr = PAZ*I_KINETIC_F2yz_Dxy_vrr+oned2z*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_KINETIC_Gy3z_Dxy_vrr = PAY*I_KINETIC_F3z_Dxy_vrr+oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr;
    Double I_KINETIC_G4z_Dxy_vrr = PAZ*I_KINETIC_F3z_Dxy_vrr+3*oned2z*I_KINETIC_D2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dxy_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Dxy_vrr;
    Double I_KINETIC_G4x_Dxz_vrr = PAX*I_KINETIC_F3x_Dxz_vrr+3*oned2z*I_KINETIC_D2x_Dxz_vrr+oned2z*I_KINETIC_F3x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Dxz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_KINETIC_G3xy_Dxz_vrr = PAY*I_KINETIC_F3x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Dxz_vrr;
    Double I_KINETIC_G3xz_Dxz_vrr = PAZ*I_KINETIC_F3x_Dxz_vrr+oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Dxz_vrr;
    Double I_KINETIC_G2x2y_Dxz_vrr = PAY*I_KINETIC_F2xy_Dxz_vrr+oned2z*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_KINETIC_G2x2z_Dxz_vrr = PAZ*I_KINETIC_F2xz_Dxz_vrr+oned2z*I_KINETIC_D2x_Dxz_vrr+oned2z*I_KINETIC_F2xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_KINETIC_Gx3y_Dxz_vrr = PAX*I_KINETIC_F3y_Dxz_vrr+oned2z*I_KINETIC_F3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dxz_vrr;
    Double I_KINETIC_Gx3z_Dxz_vrr = PAX*I_KINETIC_F3z_Dxz_vrr+oned2z*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dxz_vrr;
    Double I_KINETIC_G4y_Dxz_vrr = PAY*I_KINETIC_F3y_Dxz_vrr+3*oned2z*I_KINETIC_D2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dxz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_KINETIC_G3yz_Dxz_vrr = PAZ*I_KINETIC_F3y_Dxz_vrr+oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dxz_vrr;
    Double I_KINETIC_G2y2z_Dxz_vrr = PAZ*I_KINETIC_F2yz_Dxz_vrr+oned2z*I_KINETIC_D2y_Dxz_vrr+oned2z*I_KINETIC_F2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_KINETIC_Gy3z_Dxz_vrr = PAY*I_KINETIC_F3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dxz_vrr;
    Double I_KINETIC_G4z_Dxz_vrr = PAZ*I_KINETIC_F3z_Dxz_vrr+3*oned2z*I_KINETIC_D2z_Dxz_vrr+oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dxz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Dxz_vrr;
    Double I_KINETIC_G4x_D2y_vrr = PAX*I_KINETIC_F3x_D2y_vrr+3*oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_KINETIC_G3xy_D2y_vrr = PAY*I_KINETIC_F3x_D2y_vrr+2*oned2z*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_KINETIC_G3xz_D2y_vrr = PAZ*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2y_vrr;
    Double I_KINETIC_G2x2y_D2y_vrr = PAY*I_KINETIC_F2xy_D2y_vrr+oned2z*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_F2xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2y_vrr-bdz*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_KINETIC_G2x2z_D2y_vrr = PAZ*I_KINETIC_F2xz_D2y_vrr+oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_KINETIC_Gx3y_D2y_vrr = PAX*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2y_vrr;
    Double I_KINETIC_Gx3z_D2y_vrr = PAX*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2y_vrr;
    Double I_KINETIC_G4y_D2y_vrr = PAY*I_KINETIC_F3y_D2y_vrr+3*oned2z*I_KINETIC_D2y_D2y_vrr+2*oned2z*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_KINETIC_G3yz_D2y_vrr = PAZ*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2y_vrr;
    Double I_KINETIC_G2y2z_D2y_vrr = PAZ*I_KINETIC_F2yz_D2y_vrr+oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_KINETIC_Gy3z_D2y_vrr = PAY*I_KINETIC_F3z_D2y_vrr+2*oned2z*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2y_vrr;
    Double I_KINETIC_G4z_D2y_vrr = PAZ*I_KINETIC_F3z_D2y_vrr+3*oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_KINETIC_G4x_Dyz_vrr = PAX*I_KINETIC_F3x_Dyz_vrr+3*oned2z*I_KINETIC_D2x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Dyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_KINETIC_G3xy_Dyz_vrr = PAY*I_KINETIC_F3x_Dyz_vrr+oned2z*I_KINETIC_F3x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Dyz_vrr;
    Double I_KINETIC_G3xz_Dyz_vrr = PAZ*I_KINETIC_F3x_Dyz_vrr+oned2z*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Dyz_vrr;
    Double I_KINETIC_G2x2y_Dyz_vrr = PAY*I_KINETIC_F2xy_Dyz_vrr+oned2z*I_KINETIC_D2x_Dyz_vrr+oned2z*I_KINETIC_F2xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_KINETIC_G2x2z_Dyz_vrr = PAZ*I_KINETIC_F2xz_Dyz_vrr+oned2z*I_KINETIC_D2x_Dyz_vrr+oned2z*I_KINETIC_F2xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_KINETIC_Gx3y_Dyz_vrr = PAX*I_KINETIC_F3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dyz_vrr;
    Double I_KINETIC_Gx3z_Dyz_vrr = PAX*I_KINETIC_F3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dyz_vrr;
    Double I_KINETIC_G4y_Dyz_vrr = PAY*I_KINETIC_F3y_Dyz_vrr+3*oned2z*I_KINETIC_D2y_Dyz_vrr+oned2z*I_KINETIC_F3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_KINETIC_G3yz_Dyz_vrr = PAZ*I_KINETIC_F3y_Dyz_vrr+oned2z*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dyz_vrr;
    Double I_KINETIC_G2y2z_Dyz_vrr = PAZ*I_KINETIC_F2yz_Dyz_vrr+oned2z*I_KINETIC_D2y_Dyz_vrr+oned2z*I_KINETIC_F2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_KINETIC_Gy3z_Dyz_vrr = PAY*I_KINETIC_F3z_Dyz_vrr+oned2z*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dyz_vrr;
    Double I_KINETIC_G4z_Dyz_vrr = PAZ*I_KINETIC_F3z_Dyz_vrr+3*oned2z*I_KINETIC_D2z_Dyz_vrr+oned2z*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Dyz_vrr;
    Double I_KINETIC_G4x_D2z_vrr = PAX*I_KINETIC_F3x_D2z_vrr+3*oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_KINETIC_G3xy_D2z_vrr = PAY*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2z_vrr;
    Double I_KINETIC_G3xz_D2z_vrr = PAZ*I_KINETIC_F3x_D2z_vrr+2*oned2z*I_KINETIC_F3x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2z_vrr;
    Double I_KINETIC_G2x2y_D2z_vrr = PAY*I_KINETIC_F2xy_D2z_vrr+oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_KINETIC_G2x2z_D2z_vrr = PAZ*I_KINETIC_F2xz_D2z_vrr+oned2z*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_F2xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_KINETIC_Gx3y_D2z_vrr = PAX*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2z_vrr;
    Double I_KINETIC_Gx3z_D2z_vrr = PAX*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2z_vrr;
    Double I_KINETIC_G4y_D2z_vrr = PAY*I_KINETIC_F3y_D2z_vrr+3*oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_KINETIC_G3yz_D2z_vrr = PAZ*I_KINETIC_F3y_D2z_vrr+2*oned2z*I_KINETIC_F3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2z_vrr;
    Double I_KINETIC_G2y2z_D2z_vrr = PAZ*I_KINETIC_F2yz_D2z_vrr+oned2z*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_F2yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_KINETIC_Gy3z_D2z_vrr = PAY*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2z_vrr;
    Double I_KINETIC_G4z_D2z_vrr = PAZ*I_KINETIC_F3z_D2z_vrr+3*oned2z*I_KINETIC_D2z_D2z_vrr+2*oned2z*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_F
     * RHS shell quartet name: SQ_KINETIC_D_F
     * RHS shell quartet name: SQ_KINETIC_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     ************************************************************/
    Double I_KINETIC_G4x_F3x_vrr = PAX*I_KINETIC_F3x_F3x_vrr+3*oned2z*I_KINETIC_D2x_F3x_vrr+3*oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3x_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_KINETIC_G3xy_F3x_vrr = PAY*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3x_vrr;
    Double I_KINETIC_G3xz_F3x_vrr = PAZ*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3x_vrr;
    Double I_KINETIC_G2x2y_F3x_vrr = PAY*I_KINETIC_F2xy_F3x_vrr+oned2z*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3x_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_KINETIC_G2x2z_F3x_vrr = PAZ*I_KINETIC_F2xz_F3x_vrr+oned2z*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_KINETIC_Gx3y_F3x_vrr = PAX*I_KINETIC_F3y_F3x_vrr+3*oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3x_vrr;
    Double I_KINETIC_Gx3z_F3x_vrr = PAX*I_KINETIC_F3z_F3x_vrr+3*oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3x_vrr;
    Double I_KINETIC_G4y_F3x_vrr = PAY*I_KINETIC_F3y_F3x_vrr+3*oned2z*I_KINETIC_D2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3x_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_KINETIC_G3yz_F3x_vrr = PAZ*I_KINETIC_F3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3x_vrr;
    Double I_KINETIC_G2y2z_F3x_vrr = PAZ*I_KINETIC_F2yz_F3x_vrr+oned2z*I_KINETIC_D2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_KINETIC_Gy3z_F3x_vrr = PAY*I_KINETIC_F3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3x_vrr;
    Double I_KINETIC_G4z_F3x_vrr = PAZ*I_KINETIC_F3z_F3x_vrr+3*oned2z*I_KINETIC_D2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3x_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F3x_vrr;
    Double I_KINETIC_G4x_F2xy_vrr = PAX*I_KINETIC_F3x_F2xy_vrr+3*oned2z*I_KINETIC_D2x_F2xy_vrr+2*oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2xy_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_KINETIC_G3xy_F2xy_vrr = PAY*I_KINETIC_F3x_F2xy_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2xy_vrr;
    Double I_KINETIC_G3xz_F2xy_vrr = PAZ*I_KINETIC_F3x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2xy_vrr;
    Double I_KINETIC_G2x2y_F2xy_vrr = PAY*I_KINETIC_F2xy_F2xy_vrr+oned2z*I_KINETIC_D2x_F2xy_vrr+oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_KINETIC_G2x2z_F2xy_vrr = PAZ*I_KINETIC_F2xz_F2xy_vrr+oned2z*I_KINETIC_D2x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_KINETIC_Gx3y_F2xy_vrr = PAX*I_KINETIC_F3y_F2xy_vrr+2*oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr;
    Double I_KINETIC_Gx3z_F2xy_vrr = PAX*I_KINETIC_F3z_F2xy_vrr+2*oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr;
    Double I_KINETIC_G4y_F2xy_vrr = PAY*I_KINETIC_F3y_F2xy_vrr+3*oned2z*I_KINETIC_D2y_F2xy_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2xy_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_KINETIC_G3yz_F2xy_vrr = PAZ*I_KINETIC_F3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2xy_vrr;
    Double I_KINETIC_G2y2z_F2xy_vrr = PAZ*I_KINETIC_F2yz_F2xy_vrr+oned2z*I_KINETIC_D2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_KINETIC_Gy3z_F2xy_vrr = PAY*I_KINETIC_F3z_F2xy_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2xy_vrr;
    Double I_KINETIC_G4z_F2xy_vrr = PAZ*I_KINETIC_F3z_F2xy_vrr+3*oned2z*I_KINETIC_D2z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2xy_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F2xy_vrr;
    Double I_KINETIC_G4x_F2xz_vrr = PAX*I_KINETIC_F3x_F2xz_vrr+3*oned2z*I_KINETIC_D2x_F2xz_vrr+2*oned2z*I_KINETIC_F3x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2xz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_KINETIC_G3xy_F2xz_vrr = PAY*I_KINETIC_F3x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2xz_vrr;
    Double I_KINETIC_G3xz_F2xz_vrr = PAZ*I_KINETIC_F3x_F2xz_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2xz_vrr;
    Double I_KINETIC_G2x2y_F2xz_vrr = PAY*I_KINETIC_F2xy_F2xz_vrr+oned2z*I_KINETIC_D2x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2xz_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_KINETIC_G2x2z_F2xz_vrr = PAZ*I_KINETIC_F2xz_F2xz_vrr+oned2z*I_KINETIC_D2x_F2xz_vrr+oned2z*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_KINETIC_Gx3y_F2xz_vrr = PAX*I_KINETIC_F3y_F2xz_vrr+2*oned2z*I_KINETIC_F3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2xz_vrr;
    Double I_KINETIC_Gx3z_F2xz_vrr = PAX*I_KINETIC_F3z_F2xz_vrr+2*oned2z*I_KINETIC_F3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2xz_vrr;
    Double I_KINETIC_G4y_F2xz_vrr = PAY*I_KINETIC_F3y_F2xz_vrr+3*oned2z*I_KINETIC_D2y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2xz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F2xz_vrr;
    Double I_KINETIC_G3yz_F2xz_vrr = PAZ*I_KINETIC_F3y_F2xz_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2xz_vrr;
    Double I_KINETIC_G2y2z_F2xz_vrr = PAZ*I_KINETIC_F2yz_F2xz_vrr+oned2z*I_KINETIC_D2y_F2xz_vrr+oned2z*I_KINETIC_F2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_D2y_F2xz_vrr;
    Double I_KINETIC_Gy3z_F2xz_vrr = PAY*I_KINETIC_F3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2xz_vrr;
    Double I_KINETIC_G4z_F2xz_vrr = PAZ*I_KINETIC_F3z_F2xz_vrr+3*oned2z*I_KINETIC_D2z_F2xz_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2xz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F2xz_vrr;
    Double I_KINETIC_G4x_Fx2y_vrr = PAX*I_KINETIC_F3x_Fx2y_vrr+3*oned2z*I_KINETIC_D2x_Fx2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fx2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_KINETIC_G3xy_Fx2y_vrr = PAY*I_KINETIC_F3x_Fx2y_vrr+2*oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr;
    Double I_KINETIC_G3xz_Fx2y_vrr = PAZ*I_KINETIC_F3x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fx2y_vrr;
    Double I_KINETIC_G2x2y_Fx2y_vrr = PAY*I_KINETIC_F2xy_Fx2y_vrr+oned2z*I_KINETIC_D2x_Fx2y_vrr+2*oned2z*I_KINETIC_F2xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_KINETIC_G2x2z_Fx2y_vrr = PAZ*I_KINETIC_F2xz_Fx2y_vrr+oned2z*I_KINETIC_D2x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_KINETIC_Gx3y_Fx2y_vrr = PAX*I_KINETIC_F3y_Fx2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr;
    Double I_KINETIC_Gx3z_Fx2y_vrr = PAX*I_KINETIC_F3z_Fx2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr;
    Double I_KINETIC_G4y_Fx2y_vrr = PAY*I_KINETIC_F3y_Fx2y_vrr+3*oned2z*I_KINETIC_D2y_Fx2y_vrr+2*oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fx2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_KINETIC_G3yz_Fx2y_vrr = PAZ*I_KINETIC_F3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fx2y_vrr;
    Double I_KINETIC_G2y2z_Fx2y_vrr = PAZ*I_KINETIC_F2yz_Fx2y_vrr+oned2z*I_KINETIC_D2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_KINETIC_Gy3z_Fx2y_vrr = PAY*I_KINETIC_F3z_Fx2y_vrr+2*oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr;
    Double I_KINETIC_G4z_Fx2y_vrr = PAZ*I_KINETIC_F3z_Fx2y_vrr+3*oned2z*I_KINETIC_D2z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fx2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Fx2y_vrr;
    Double I_KINETIC_G4x_Fxyz_vrr = PAX*I_KINETIC_F3x_Fxyz_vrr+3*oned2z*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_F3x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fxyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Fxyz_vrr;
    Double I_KINETIC_G3xy_Fxyz_vrr = PAY*I_KINETIC_F3x_Fxyz_vrr+oned2z*I_KINETIC_F3x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr;
    Double I_KINETIC_G3xz_Fxyz_vrr = PAZ*I_KINETIC_F3x_Fxyz_vrr+oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fxyz_vrr;
    Double I_KINETIC_G2x2y_Fxyz_vrr = PAY*I_KINETIC_F2xy_Fxyz_vrr+oned2z*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_F2xy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fxyz_vrr;
    Double I_KINETIC_G2x2z_Fxyz_vrr = PAZ*I_KINETIC_F2xz_Fxyz_vrr+oned2z*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_F2xz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fxyz_vrr;
    Double I_KINETIC_Gx3y_Fxyz_vrr = PAX*I_KINETIC_F3y_Fxyz_vrr+oned2z*I_KINETIC_F3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr;
    Double I_KINETIC_Gx3z_Fxyz_vrr = PAX*I_KINETIC_F3z_Fxyz_vrr+oned2z*I_KINETIC_F3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr;
    Double I_KINETIC_G4y_Fxyz_vrr = PAY*I_KINETIC_F3y_Fxyz_vrr+3*oned2z*I_KINETIC_D2y_Fxyz_vrr+oned2z*I_KINETIC_F3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fxyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Fxyz_vrr;
    Double I_KINETIC_G3yz_Fxyz_vrr = PAZ*I_KINETIC_F3y_Fxyz_vrr+oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fxyz_vrr;
    Double I_KINETIC_G2y2z_Fxyz_vrr = PAZ*I_KINETIC_F2yz_Fxyz_vrr+oned2z*I_KINETIC_D2y_Fxyz_vrr+oned2z*I_KINETIC_F2yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_D2y_Fxyz_vrr;
    Double I_KINETIC_Gy3z_Fxyz_vrr = PAY*I_KINETIC_F3z_Fxyz_vrr+oned2z*I_KINETIC_F3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr;
    Double I_KINETIC_G4z_Fxyz_vrr = PAZ*I_KINETIC_F3z_Fxyz_vrr+3*oned2z*I_KINETIC_D2z_Fxyz_vrr+oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fxyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Fxyz_vrr;
    Double I_KINETIC_G4x_Fx2z_vrr = PAX*I_KINETIC_F3x_Fx2z_vrr+3*oned2z*I_KINETIC_D2x_Fx2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fx2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_KINETIC_G3xy_Fx2z_vrr = PAY*I_KINETIC_F3x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr;
    Double I_KINETIC_G3xz_Fx2z_vrr = PAZ*I_KINETIC_F3x_Fx2z_vrr+2*oned2z*I_KINETIC_F3x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fx2z_vrr;
    Double I_KINETIC_G2x2y_Fx2z_vrr = PAY*I_KINETIC_F2xy_Fx2z_vrr+oned2z*I_KINETIC_D2x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_KINETIC_G2x2z_Fx2z_vrr = PAZ*I_KINETIC_F2xz_Fx2z_vrr+oned2z*I_KINETIC_D2x_Fx2z_vrr+2*oned2z*I_KINETIC_F2xz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_KINETIC_Gx3y_Fx2z_vrr = PAX*I_KINETIC_F3y_Fx2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr;
    Double I_KINETIC_Gx3z_Fx2z_vrr = PAX*I_KINETIC_F3z_Fx2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr;
    Double I_KINETIC_G4y_Fx2z_vrr = PAY*I_KINETIC_F3y_Fx2z_vrr+3*oned2z*I_KINETIC_D2y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fx2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_KINETIC_G3yz_Fx2z_vrr = PAZ*I_KINETIC_F3y_Fx2z_vrr+2*oned2z*I_KINETIC_F3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fx2z_vrr;
    Double I_KINETIC_G2y2z_Fx2z_vrr = PAZ*I_KINETIC_F2yz_Fx2z_vrr+oned2z*I_KINETIC_D2y_Fx2z_vrr+2*oned2z*I_KINETIC_F2yz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_KINETIC_Gy3z_Fx2z_vrr = PAY*I_KINETIC_F3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr;
    Double I_KINETIC_G4z_Fx2z_vrr = PAZ*I_KINETIC_F3z_Fx2z_vrr+3*oned2z*I_KINETIC_D2z_Fx2z_vrr+2*oned2z*I_KINETIC_F3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fx2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Fx2z_vrr;
    Double I_KINETIC_G4x_F3y_vrr = PAX*I_KINETIC_F3x_F3y_vrr+3*oned2z*I_KINETIC_D2x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3y_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_KINETIC_G3xy_F3y_vrr = PAY*I_KINETIC_F3x_F3y_vrr+3*oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3y_vrr;
    Double I_KINETIC_G3xz_F3y_vrr = PAZ*I_KINETIC_F3x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3y_vrr;
    Double I_KINETIC_G2x2y_F3y_vrr = PAY*I_KINETIC_F2xy_F3y_vrr+oned2z*I_KINETIC_D2x_F3y_vrr+3*oned2z*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3y_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_KINETIC_G2x2z_F3y_vrr = PAZ*I_KINETIC_F2xz_F3y_vrr+oned2z*I_KINETIC_D2x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_KINETIC_Gx3y_F3y_vrr = PAX*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3y_vrr;
    Double I_KINETIC_Gx3z_F3y_vrr = PAX*I_KINETIC_F3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3y_vrr;
    Double I_KINETIC_G4y_F3y_vrr = PAY*I_KINETIC_F3y_F3y_vrr+3*oned2z*I_KINETIC_D2y_F3y_vrr+3*oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3y_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_KINETIC_G3yz_F3y_vrr = PAZ*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3y_vrr;
    Double I_KINETIC_G2y2z_F3y_vrr = PAZ*I_KINETIC_F2yz_F3y_vrr+oned2z*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_KINETIC_Gy3z_F3y_vrr = PAY*I_KINETIC_F3z_F3y_vrr+3*oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3y_vrr;
    Double I_KINETIC_G4z_F3y_vrr = PAZ*I_KINETIC_F3z_F3y_vrr+3*oned2z*I_KINETIC_D2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3y_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F3y_vrr;
    Double I_KINETIC_G4x_F2yz_vrr = PAX*I_KINETIC_F3x_F2yz_vrr+3*oned2z*I_KINETIC_D2x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_KINETIC_G3xy_F2yz_vrr = PAY*I_KINETIC_F3x_F2yz_vrr+2*oned2z*I_KINETIC_F3x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2yz_vrr;
    Double I_KINETIC_G3xz_F2yz_vrr = PAZ*I_KINETIC_F3x_F2yz_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2yz_vrr;
    Double I_KINETIC_G2x2y_F2yz_vrr = PAY*I_KINETIC_F2xy_F2yz_vrr+oned2z*I_KINETIC_D2x_F2yz_vrr+2*oned2z*I_KINETIC_F2xy_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2yz_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_KINETIC_G2x2z_F2yz_vrr = PAZ*I_KINETIC_F2xz_F2yz_vrr+oned2z*I_KINETIC_D2x_F2yz_vrr+oned2z*I_KINETIC_F2xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_KINETIC_Gx3y_F2yz_vrr = PAX*I_KINETIC_F3y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2yz_vrr;
    Double I_KINETIC_Gx3z_F2yz_vrr = PAX*I_KINETIC_F3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2yz_vrr;
    Double I_KINETIC_G4y_F2yz_vrr = PAY*I_KINETIC_F3y_F2yz_vrr+3*oned2z*I_KINETIC_D2y_F2yz_vrr+2*oned2z*I_KINETIC_F3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F2yz_vrr;
    Double I_KINETIC_G3yz_F2yz_vrr = PAZ*I_KINETIC_F3y_F2yz_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2yz_vrr;
    Double I_KINETIC_G2y2z_F2yz_vrr = PAZ*I_KINETIC_F2yz_F2yz_vrr+oned2z*I_KINETIC_D2y_F2yz_vrr+oned2z*I_KINETIC_F2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_D2y_F2yz_vrr;
    Double I_KINETIC_Gy3z_F2yz_vrr = PAY*I_KINETIC_F3z_F2yz_vrr+2*oned2z*I_KINETIC_F3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2yz_vrr;
    Double I_KINETIC_G4z_F2yz_vrr = PAZ*I_KINETIC_F3z_F2yz_vrr+3*oned2z*I_KINETIC_D2z_F2yz_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F2yz_vrr;
    Double I_KINETIC_G4x_Fy2z_vrr = PAX*I_KINETIC_F3x_Fy2z_vrr+3*oned2z*I_KINETIC_D2x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fy2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_KINETIC_G3xy_Fy2z_vrr = PAY*I_KINETIC_F3x_Fy2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr;
    Double I_KINETIC_G3xz_Fy2z_vrr = PAZ*I_KINETIC_F3x_Fy2z_vrr+2*oned2z*I_KINETIC_F3x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fy2z_vrr;
    Double I_KINETIC_G2x2y_Fy2z_vrr = PAY*I_KINETIC_F2xy_Fy2z_vrr+oned2z*I_KINETIC_D2x_Fy2z_vrr+oned2z*I_KINETIC_F2xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_KINETIC_G2x2z_Fy2z_vrr = PAZ*I_KINETIC_F2xz_Fy2z_vrr+oned2z*I_KINETIC_D2x_Fy2z_vrr+2*oned2z*I_KINETIC_F2xz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_KINETIC_Gx3y_Fy2z_vrr = PAX*I_KINETIC_F3y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr;
    Double I_KINETIC_Gx3z_Fy2z_vrr = PAX*I_KINETIC_F3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr;
    Double I_KINETIC_G4y_Fy2z_vrr = PAY*I_KINETIC_F3y_Fy2z_vrr+3*oned2z*I_KINETIC_D2y_Fy2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fy2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Fy2z_vrr;
    Double I_KINETIC_G3yz_Fy2z_vrr = PAZ*I_KINETIC_F3y_Fy2z_vrr+2*oned2z*I_KINETIC_F3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fy2z_vrr;
    Double I_KINETIC_G2y2z_Fy2z_vrr = PAZ*I_KINETIC_F2yz_Fy2z_vrr+oned2z*I_KINETIC_D2y_Fy2z_vrr+2*oned2z*I_KINETIC_F2yz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_D2y_Fy2z_vrr;
    Double I_KINETIC_Gy3z_Fy2z_vrr = PAY*I_KINETIC_F3z_Fy2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr;
    Double I_KINETIC_G4z_Fy2z_vrr = PAZ*I_KINETIC_F3z_Fy2z_vrr+3*oned2z*I_KINETIC_D2z_Fy2z_vrr+2*oned2z*I_KINETIC_F3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fy2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Fy2z_vrr;
    Double I_KINETIC_G4x_F3z_vrr = PAX*I_KINETIC_F3x_F3z_vrr+3*oned2z*I_KINETIC_D2x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_KINETIC_G3xy_F3z_vrr = PAY*I_KINETIC_F3x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3z_vrr;
    Double I_KINETIC_G3xz_F3z_vrr = PAZ*I_KINETIC_F3x_F3z_vrr+3*oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3z_vrr;
    Double I_KINETIC_G2x2y_F3z_vrr = PAY*I_KINETIC_F2xy_F3z_vrr+oned2z*I_KINETIC_D2x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3z_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_KINETIC_G2x2z_F3z_vrr = PAZ*I_KINETIC_F2xz_F3z_vrr+oned2z*I_KINETIC_D2x_F3z_vrr+3*oned2z*I_KINETIC_F2xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_KINETIC_Gx3y_F3z_vrr = PAX*I_KINETIC_F3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3z_vrr;
    Double I_KINETIC_Gx3z_F3z_vrr = PAX*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3z_vrr;
    Double I_KINETIC_G4y_F3z_vrr = PAY*I_KINETIC_F3y_F3z_vrr+3*oned2z*I_KINETIC_D2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_KINETIC_G3yz_F3z_vrr = PAZ*I_KINETIC_F3y_F3z_vrr+3*oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3z_vrr;
    Double I_KINETIC_G2y2z_F3z_vrr = PAZ*I_KINETIC_F2yz_F3z_vrr+oned2z*I_KINETIC_D2y_F3z_vrr+3*oned2z*I_KINETIC_F2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_KINETIC_Gy3z_F3z_vrr = PAY*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3z_vrr;
    Double I_KINETIC_G4z_F3z_vrr = PAZ*I_KINETIC_F3z_F3z_vrr+3*oned2z*I_KINETIC_D2z_F3z_vrr+3*oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_F
     * RHS shell quartet name: SQ_KINETIC_F_F
     * RHS shell quartet name: SQ_KINETIC_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
     ************************************************************/
    Double I_KINETIC_H5x_F3x_vrr = PAX*I_KINETIC_G4x_F3x_vrr+4*oned2z*I_KINETIC_F3x_F3x_vrr+3*oned2z*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F3x_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_KINETIC_H4xy_F3x_vrr = PAY*I_KINETIC_G4x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F3x_vrr;
    Double I_KINETIC_H4xz_F3x_vrr = PAZ*I_KINETIC_G4x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F3x_vrr;
    Double I_KINETIC_H3x2y_F3x_vrr = PAY*I_KINETIC_G3xy_F3x_vrr+oned2z*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F3x_vrr-bdz*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_KINETIC_H3xyz_F3x_vrr = PAZ*I_KINETIC_G3xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F3x_vrr;
    Double I_KINETIC_H3x2z_F3x_vrr = PAZ*I_KINETIC_G3xz_F3x_vrr+oned2z*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_KINETIC_H2x3y_F3x_vrr = PAX*I_KINETIC_Gx3y_F3x_vrr+oned2z*I_KINETIC_F3y_F3x_vrr+3*oned2z*I_KINETIC_Gx3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F3x_vrr-bdz*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_KINETIC_H2x2yz_F3x_vrr = PAZ*I_KINETIC_G2x2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F3x_vrr;
    Double I_KINETIC_H2xy2z_F3x_vrr = PAY*I_KINETIC_G2x2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F3x_vrr;
    Double I_KINETIC_H2x3z_F3x_vrr = PAX*I_KINETIC_Gx3z_F3x_vrr+oned2z*I_KINETIC_F3z_F3x_vrr+3*oned2z*I_KINETIC_Gx3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F3x_vrr-bdz*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_KINETIC_Hx4y_F3x_vrr = PAX*I_KINETIC_G4y_F3x_vrr+3*oned2z*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F3x_vrr;
    Double I_KINETIC_Hx3yz_F3x_vrr = PAZ*I_KINETIC_Gx3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F3x_vrr;
    Double I_KINETIC_Hx2y2z_F3x_vrr = PAX*I_KINETIC_G2y2z_F3x_vrr+3*oned2z*I_KINETIC_G2y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F3x_vrr;
    Double I_KINETIC_Hxy3z_F3x_vrr = PAY*I_KINETIC_Gx3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F3x_vrr;
    Double I_KINETIC_Hx4z_F3x_vrr = PAX*I_KINETIC_G4z_F3x_vrr+3*oned2z*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F3x_vrr;
    Double I_KINETIC_H5y_F3x_vrr = PAY*I_KINETIC_G4y_F3x_vrr+4*oned2z*I_KINETIC_F3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F3x_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_KINETIC_H4yz_F3x_vrr = PAZ*I_KINETIC_G4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F3x_vrr;
    Double I_KINETIC_H3y2z_F3x_vrr = PAZ*I_KINETIC_G3yz_F3x_vrr+oned2z*I_KINETIC_F3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_KINETIC_H2y3z_F3x_vrr = PAY*I_KINETIC_Gy3z_F3x_vrr+oned2z*I_KINETIC_F3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F3x_vrr-bdz*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_KINETIC_Hy4z_F3x_vrr = PAY*I_KINETIC_G4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F3x_vrr;
    Double I_KINETIC_H5z_F3x_vrr = PAZ*I_KINETIC_G4z_F3x_vrr+4*oned2z*I_KINETIC_F3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F3x_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_KINETIC_H5x_F2xy_vrr = PAX*I_KINETIC_G4x_F2xy_vrr+4*oned2z*I_KINETIC_F3x_F2xy_vrr+2*oned2z*I_KINETIC_G4x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F2xy_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_KINETIC_H4xy_F2xy_vrr = PAY*I_KINETIC_G4x_F2xy_vrr+oned2z*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F2xy_vrr;
    Double I_KINETIC_H4xz_F2xy_vrr = PAZ*I_KINETIC_G4x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F2xy_vrr;
    Double I_KINETIC_H3x2y_F2xy_vrr = PAY*I_KINETIC_G3xy_F2xy_vrr+oned2z*I_KINETIC_F3x_F2xy_vrr+oned2z*I_KINETIC_G3xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F2xy_vrr-bdz*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_KINETIC_H3xyz_F2xy_vrr = PAZ*I_KINETIC_G3xy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F2xy_vrr;
    Double I_KINETIC_H3x2z_F2xy_vrr = PAZ*I_KINETIC_G3xz_F2xy_vrr+oned2z*I_KINETIC_F3x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_KINETIC_H2x3y_F2xy_vrr = PAX*I_KINETIC_Gx3y_F2xy_vrr+oned2z*I_KINETIC_F3y_F2xy_vrr+2*oned2z*I_KINETIC_Gx3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F2xy_vrr-bdz*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_KINETIC_H2x2yz_F2xy_vrr = PAZ*I_KINETIC_G2x2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F2xy_vrr;
    Double I_KINETIC_H2xy2z_F2xy_vrr = PAY*I_KINETIC_G2x2z_F2xy_vrr+oned2z*I_KINETIC_G2x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F2xy_vrr;
    Double I_KINETIC_H2x3z_F2xy_vrr = PAX*I_KINETIC_Gx3z_F2xy_vrr+oned2z*I_KINETIC_F3z_F2xy_vrr+2*oned2z*I_KINETIC_Gx3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_F3z_F2xy_vrr;
    Double I_KINETIC_Hx4y_F2xy_vrr = PAX*I_KINETIC_G4y_F2xy_vrr+2*oned2z*I_KINETIC_G4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F2xy_vrr;
    Double I_KINETIC_Hx3yz_F2xy_vrr = PAZ*I_KINETIC_Gx3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F2xy_vrr;
    Double I_KINETIC_Hx2y2z_F2xy_vrr = PAX*I_KINETIC_G2y2z_F2xy_vrr+2*oned2z*I_KINETIC_G2y2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F2xy_vrr;
    Double I_KINETIC_Hxy3z_F2xy_vrr = PAY*I_KINETIC_Gx3z_F2xy_vrr+oned2z*I_KINETIC_Gx3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F2xy_vrr;
    Double I_KINETIC_Hx4z_F2xy_vrr = PAX*I_KINETIC_G4z_F2xy_vrr+2*oned2z*I_KINETIC_G4z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F2xy_vrr;
    Double I_KINETIC_H5y_F2xy_vrr = PAY*I_KINETIC_G4y_F2xy_vrr+4*oned2z*I_KINETIC_F3y_F2xy_vrr+oned2z*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F2xy_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_KINETIC_H4yz_F2xy_vrr = PAZ*I_KINETIC_G4y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F2xy_vrr;
    Double I_KINETIC_H3y2z_F2xy_vrr = PAZ*I_KINETIC_G3yz_F2xy_vrr+oned2z*I_KINETIC_F3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_KINETIC_H2y3z_F2xy_vrr = PAY*I_KINETIC_Gy3z_F2xy_vrr+oned2z*I_KINETIC_F3z_F2xy_vrr+oned2z*I_KINETIC_Gy3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_F3z_F2xy_vrr;
    Double I_KINETIC_Hy4z_F2xy_vrr = PAY*I_KINETIC_G4z_F2xy_vrr+oned2z*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F2xy_vrr;
    Double I_KINETIC_H5z_F2xy_vrr = PAZ*I_KINETIC_G4z_F2xy_vrr+4*oned2z*I_KINETIC_F3z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F2xy_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_F2xy_vrr;
    Double I_KINETIC_H5x_F2xz_vrr = PAX*I_KINETIC_G4x_F2xz_vrr+4*oned2z*I_KINETIC_F3x_F2xz_vrr+2*oned2z*I_KINETIC_G4x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F2xz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_KINETIC_H4xy_F2xz_vrr = PAY*I_KINETIC_G4x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F2xz_vrr;
    Double I_KINETIC_H4xz_F2xz_vrr = PAZ*I_KINETIC_G4x_F2xz_vrr+oned2z*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F2xz_vrr;
    Double I_KINETIC_H3x2y_F2xz_vrr = PAY*I_KINETIC_G3xy_F2xz_vrr+oned2z*I_KINETIC_F3x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F2xz_vrr-bdz*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_KINETIC_H3xyz_F2xz_vrr = PAZ*I_KINETIC_G3xy_F2xz_vrr+oned2z*I_KINETIC_G3xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F2xz_vrr;
    Double I_KINETIC_H3x2z_F2xz_vrr = PAZ*I_KINETIC_G3xz_F2xz_vrr+oned2z*I_KINETIC_F3x_F2xz_vrr+oned2z*I_KINETIC_G3xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_KINETIC_H2x3y_F2xz_vrr = PAX*I_KINETIC_Gx3y_F2xz_vrr+oned2z*I_KINETIC_F3y_F2xz_vrr+2*oned2z*I_KINETIC_Gx3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F2xz_vrr-bdz*I_TWOBODYOVERLAP_F3y_F2xz_vrr;
    Double I_KINETIC_H2x2yz_F2xz_vrr = PAZ*I_KINETIC_G2x2y_F2xz_vrr+oned2z*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F2xz_vrr;
    Double I_KINETIC_H2xy2z_F2xz_vrr = PAY*I_KINETIC_G2x2z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F2xz_vrr;
    Double I_KINETIC_H2x3z_F2xz_vrr = PAX*I_KINETIC_Gx3z_F2xz_vrr+oned2z*I_KINETIC_F3z_F2xz_vrr+2*oned2z*I_KINETIC_Gx3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_KINETIC_Hx4y_F2xz_vrr = PAX*I_KINETIC_G4y_F2xz_vrr+2*oned2z*I_KINETIC_G4y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F2xz_vrr;
    Double I_KINETIC_Hx3yz_F2xz_vrr = PAZ*I_KINETIC_Gx3y_F2xz_vrr+oned2z*I_KINETIC_Gx3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F2xz_vrr;
    Double I_KINETIC_Hx2y2z_F2xz_vrr = PAX*I_KINETIC_G2y2z_F2xz_vrr+2*oned2z*I_KINETIC_G2y2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F2xz_vrr;
    Double I_KINETIC_Hxy3z_F2xz_vrr = PAY*I_KINETIC_Gx3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F2xz_vrr;
    Double I_KINETIC_Hx4z_F2xz_vrr = PAX*I_KINETIC_G4z_F2xz_vrr+2*oned2z*I_KINETIC_G4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F2xz_vrr;
    Double I_KINETIC_H5y_F2xz_vrr = PAY*I_KINETIC_G4y_F2xz_vrr+4*oned2z*I_KINETIC_F3y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F2xz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_F2xz_vrr;
    Double I_KINETIC_H4yz_F2xz_vrr = PAZ*I_KINETIC_G4y_F2xz_vrr+oned2z*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F2xz_vrr;
    Double I_KINETIC_H3y2z_F2xz_vrr = PAZ*I_KINETIC_G3yz_F2xz_vrr+oned2z*I_KINETIC_F3y_F2xz_vrr+oned2z*I_KINETIC_G3yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_F3y_F2xz_vrr;
    Double I_KINETIC_H2y3z_F2xz_vrr = PAY*I_KINETIC_Gy3z_F2xz_vrr+oned2z*I_KINETIC_F3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_KINETIC_Hy4z_F2xz_vrr = PAY*I_KINETIC_G4z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F2xz_vrr;
    Double I_KINETIC_H5z_F2xz_vrr = PAZ*I_KINETIC_G4z_F2xz_vrr+4*oned2z*I_KINETIC_F3z_F2xz_vrr+oned2z*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F2xz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_KINETIC_H5x_Fx2y_vrr = PAX*I_KINETIC_G4x_Fx2y_vrr+4*oned2z*I_KINETIC_F3x_Fx2y_vrr+oned2z*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fx2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_KINETIC_H4xy_Fx2y_vrr = PAY*I_KINETIC_G4x_Fx2y_vrr+2*oned2z*I_KINETIC_G4x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fx2y_vrr;
    Double I_KINETIC_H4xz_Fx2y_vrr = PAZ*I_KINETIC_G4x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fx2y_vrr;
    Double I_KINETIC_H3x2y_Fx2y_vrr = PAY*I_KINETIC_G3xy_Fx2y_vrr+oned2z*I_KINETIC_F3x_Fx2y_vrr+2*oned2z*I_KINETIC_G3xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_KINETIC_H3xyz_Fx2y_vrr = PAZ*I_KINETIC_G3xy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Fx2y_vrr;
    Double I_KINETIC_H3x2z_Fx2y_vrr = PAZ*I_KINETIC_G3xz_Fx2y_vrr+oned2z*I_KINETIC_F3x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_KINETIC_H2x3y_Fx2y_vrr = PAX*I_KINETIC_Gx3y_Fx2y_vrr+oned2z*I_KINETIC_F3y_Fx2y_vrr+oned2z*I_KINETIC_Gx3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_KINETIC_H2x2yz_Fx2y_vrr = PAZ*I_KINETIC_G2x2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fx2y_vrr;
    Double I_KINETIC_H2xy2z_Fx2y_vrr = PAY*I_KINETIC_G2x2z_Fx2y_vrr+2*oned2z*I_KINETIC_G2x2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Fx2y_vrr;
    Double I_KINETIC_H2x3z_Fx2y_vrr = PAX*I_KINETIC_Gx3z_Fx2y_vrr+oned2z*I_KINETIC_F3z_Fx2y_vrr+oned2z*I_KINETIC_Gx3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_F3z_Fx2y_vrr;
    Double I_KINETIC_Hx4y_Fx2y_vrr = PAX*I_KINETIC_G4y_Fx2y_vrr+oned2z*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fx2y_vrr;
    Double I_KINETIC_Hx3yz_Fx2y_vrr = PAZ*I_KINETIC_Gx3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Fx2y_vrr;
    Double I_KINETIC_Hx2y2z_Fx2y_vrr = PAX*I_KINETIC_G2y2z_Fx2y_vrr+oned2z*I_KINETIC_G2y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Fx2y_vrr;
    Double I_KINETIC_Hxy3z_Fx2y_vrr = PAY*I_KINETIC_Gx3z_Fx2y_vrr+2*oned2z*I_KINETIC_Gx3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Fx2y_vrr;
    Double I_KINETIC_Hx4z_Fx2y_vrr = PAX*I_KINETIC_G4z_Fx2y_vrr+oned2z*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fx2y_vrr;
    Double I_KINETIC_H5y_Fx2y_vrr = PAY*I_KINETIC_G4y_Fx2y_vrr+4*oned2z*I_KINETIC_F3y_Fx2y_vrr+2*oned2z*I_KINETIC_G4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fx2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_KINETIC_H4yz_Fx2y_vrr = PAZ*I_KINETIC_G4y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fx2y_vrr;
    Double I_KINETIC_H3y2z_Fx2y_vrr = PAZ*I_KINETIC_G3yz_Fx2y_vrr+oned2z*I_KINETIC_F3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_KINETIC_H2y3z_Fx2y_vrr = PAY*I_KINETIC_Gy3z_Fx2y_vrr+oned2z*I_KINETIC_F3z_Fx2y_vrr+2*oned2z*I_KINETIC_Gy3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_F3z_Fx2y_vrr;
    Double I_KINETIC_Hy4z_Fx2y_vrr = PAY*I_KINETIC_G4z_Fx2y_vrr+2*oned2z*I_KINETIC_G4z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fx2y_vrr;
    Double I_KINETIC_H5z_Fx2y_vrr = PAZ*I_KINETIC_G4z_Fx2y_vrr+4*oned2z*I_KINETIC_F3z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fx2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Fx2y_vrr;
    Double I_KINETIC_H5x_Fxyz_vrr = PAX*I_KINETIC_G4x_Fxyz_vrr+4*oned2z*I_KINETIC_F3x_Fxyz_vrr+oned2z*I_KINETIC_G4x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fxyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Fxyz_vrr;
    Double I_KINETIC_H4xy_Fxyz_vrr = PAY*I_KINETIC_G4x_Fxyz_vrr+oned2z*I_KINETIC_G4x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fxyz_vrr;
    Double I_KINETIC_H4xz_Fxyz_vrr = PAZ*I_KINETIC_G4x_Fxyz_vrr+oned2z*I_KINETIC_G4x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fxyz_vrr;
    Double I_KINETIC_H3x2y_Fxyz_vrr = PAY*I_KINETIC_G3xy_Fxyz_vrr+oned2z*I_KINETIC_F3x_Fxyz_vrr+oned2z*I_KINETIC_G3xy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Fxyz_vrr;
    Double I_KINETIC_H3xyz_Fxyz_vrr = PAZ*I_KINETIC_G3xy_Fxyz_vrr+oned2z*I_KINETIC_G3xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Fxyz_vrr;
    Double I_KINETIC_H3x2z_Fxyz_vrr = PAZ*I_KINETIC_G3xz_Fxyz_vrr+oned2z*I_KINETIC_F3x_Fxyz_vrr+oned2z*I_KINETIC_G3xz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Fxyz_vrr;
    Double I_KINETIC_H2x3y_Fxyz_vrr = PAX*I_KINETIC_Gx3y_Fxyz_vrr+oned2z*I_KINETIC_F3y_Fxyz_vrr+oned2z*I_KINETIC_Gx3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Fxyz_vrr;
    Double I_KINETIC_H2x2yz_Fxyz_vrr = PAZ*I_KINETIC_G2x2y_Fxyz_vrr+oned2z*I_KINETIC_G2x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fxyz_vrr;
    Double I_KINETIC_H2xy2z_Fxyz_vrr = PAY*I_KINETIC_G2x2z_Fxyz_vrr+oned2z*I_KINETIC_G2x2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Fxyz_vrr;
    Double I_KINETIC_H2x3z_Fxyz_vrr = PAX*I_KINETIC_Gx3z_Fxyz_vrr+oned2z*I_KINETIC_F3z_Fxyz_vrr+oned2z*I_KINETIC_Gx3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Fxyz_vrr;
    Double I_KINETIC_Hx4y_Fxyz_vrr = PAX*I_KINETIC_G4y_Fxyz_vrr+oned2z*I_KINETIC_G4y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fxyz_vrr;
    Double I_KINETIC_Hx3yz_Fxyz_vrr = PAZ*I_KINETIC_Gx3y_Fxyz_vrr+oned2z*I_KINETIC_Gx3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Fxyz_vrr;
    Double I_KINETIC_Hx2y2z_Fxyz_vrr = PAX*I_KINETIC_G2y2z_Fxyz_vrr+oned2z*I_KINETIC_G2y2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Fxyz_vrr;
    Double I_KINETIC_Hxy3z_Fxyz_vrr = PAY*I_KINETIC_Gx3z_Fxyz_vrr+oned2z*I_KINETIC_Gx3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Fxyz_vrr;
    Double I_KINETIC_Hx4z_Fxyz_vrr = PAX*I_KINETIC_G4z_Fxyz_vrr+oned2z*I_KINETIC_G4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fxyz_vrr;
    Double I_KINETIC_H5y_Fxyz_vrr = PAY*I_KINETIC_G4y_Fxyz_vrr+4*oned2z*I_KINETIC_F3y_Fxyz_vrr+oned2z*I_KINETIC_G4y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fxyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Fxyz_vrr;
    Double I_KINETIC_H4yz_Fxyz_vrr = PAZ*I_KINETIC_G4y_Fxyz_vrr+oned2z*I_KINETIC_G4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fxyz_vrr;
    Double I_KINETIC_H3y2z_Fxyz_vrr = PAZ*I_KINETIC_G3yz_Fxyz_vrr+oned2z*I_KINETIC_F3y_Fxyz_vrr+oned2z*I_KINETIC_G3yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Fxyz_vrr;
    Double I_KINETIC_H2y3z_Fxyz_vrr = PAY*I_KINETIC_Gy3z_Fxyz_vrr+oned2z*I_KINETIC_F3z_Fxyz_vrr+oned2z*I_KINETIC_Gy3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Fxyz_vrr;
    Double I_KINETIC_Hy4z_Fxyz_vrr = PAY*I_KINETIC_G4z_Fxyz_vrr+oned2z*I_KINETIC_G4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fxyz_vrr;
    Double I_KINETIC_H5z_Fxyz_vrr = PAZ*I_KINETIC_G4z_Fxyz_vrr+4*oned2z*I_KINETIC_F3z_Fxyz_vrr+oned2z*I_KINETIC_G4z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fxyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Fxyz_vrr;
    Double I_KINETIC_H5x_Fx2z_vrr = PAX*I_KINETIC_G4x_Fx2z_vrr+4*oned2z*I_KINETIC_F3x_Fx2z_vrr+oned2z*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fx2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_KINETIC_H4xy_Fx2z_vrr = PAY*I_KINETIC_G4x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fx2z_vrr;
    Double I_KINETIC_H4xz_Fx2z_vrr = PAZ*I_KINETIC_G4x_Fx2z_vrr+2*oned2z*I_KINETIC_G4x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fx2z_vrr;
    Double I_KINETIC_H3x2y_Fx2z_vrr = PAY*I_KINETIC_G3xy_Fx2z_vrr+oned2z*I_KINETIC_F3x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_KINETIC_H3xyz_Fx2z_vrr = PAZ*I_KINETIC_G3xy_Fx2z_vrr+2*oned2z*I_KINETIC_G3xy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Fx2z_vrr;
    Double I_KINETIC_H3x2z_Fx2z_vrr = PAZ*I_KINETIC_G3xz_Fx2z_vrr+oned2z*I_KINETIC_F3x_Fx2z_vrr+2*oned2z*I_KINETIC_G3xz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_KINETIC_H2x3y_Fx2z_vrr = PAX*I_KINETIC_Gx3y_Fx2z_vrr+oned2z*I_KINETIC_F3y_Fx2z_vrr+oned2z*I_KINETIC_Gx3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Fx2z_vrr;
    Double I_KINETIC_H2x2yz_Fx2z_vrr = PAZ*I_KINETIC_G2x2y_Fx2z_vrr+2*oned2z*I_KINETIC_G2x2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fx2z_vrr;
    Double I_KINETIC_H2xy2z_Fx2z_vrr = PAY*I_KINETIC_G2x2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Fx2z_vrr;
    Double I_KINETIC_H2x3z_Fx2z_vrr = PAX*I_KINETIC_Gx3z_Fx2z_vrr+oned2z*I_KINETIC_F3z_Fx2z_vrr+oned2z*I_KINETIC_Gx3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_KINETIC_Hx4y_Fx2z_vrr = PAX*I_KINETIC_G4y_Fx2z_vrr+oned2z*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fx2z_vrr;
    Double I_KINETIC_Hx3yz_Fx2z_vrr = PAZ*I_KINETIC_Gx3y_Fx2z_vrr+2*oned2z*I_KINETIC_Gx3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Fx2z_vrr;
    Double I_KINETIC_Hx2y2z_Fx2z_vrr = PAX*I_KINETIC_G2y2z_Fx2z_vrr+oned2z*I_KINETIC_G2y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Fx2z_vrr;
    Double I_KINETIC_Hxy3z_Fx2z_vrr = PAY*I_KINETIC_Gx3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Fx2z_vrr;
    Double I_KINETIC_Hx4z_Fx2z_vrr = PAX*I_KINETIC_G4z_Fx2z_vrr+oned2z*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fx2z_vrr;
    Double I_KINETIC_H5y_Fx2z_vrr = PAY*I_KINETIC_G4y_Fx2z_vrr+4*oned2z*I_KINETIC_F3y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fx2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Fx2z_vrr;
    Double I_KINETIC_H4yz_Fx2z_vrr = PAZ*I_KINETIC_G4y_Fx2z_vrr+2*oned2z*I_KINETIC_G4y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fx2z_vrr;
    Double I_KINETIC_H3y2z_Fx2z_vrr = PAZ*I_KINETIC_G3yz_Fx2z_vrr+oned2z*I_KINETIC_F3y_Fx2z_vrr+2*oned2z*I_KINETIC_G3yz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Fx2z_vrr;
    Double I_KINETIC_H2y3z_Fx2z_vrr = PAY*I_KINETIC_Gy3z_Fx2z_vrr+oned2z*I_KINETIC_F3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_KINETIC_Hy4z_Fx2z_vrr = PAY*I_KINETIC_G4z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fx2z_vrr;
    Double I_KINETIC_H5z_Fx2z_vrr = PAZ*I_KINETIC_G4z_Fx2z_vrr+4*oned2z*I_KINETIC_F3z_Fx2z_vrr+2*oned2z*I_KINETIC_G4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fx2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_KINETIC_H5x_F3y_vrr = PAX*I_KINETIC_G4x_F3y_vrr+4*oned2z*I_KINETIC_F3x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F3y_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_KINETIC_H4xy_F3y_vrr = PAY*I_KINETIC_G4x_F3y_vrr+3*oned2z*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F3y_vrr;
    Double I_KINETIC_H4xz_F3y_vrr = PAZ*I_KINETIC_G4x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F3y_vrr;
    Double I_KINETIC_H3x2y_F3y_vrr = PAY*I_KINETIC_G3xy_F3y_vrr+oned2z*I_KINETIC_F3x_F3y_vrr+3*oned2z*I_KINETIC_G3xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F3y_vrr-bdz*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_KINETIC_H3xyz_F3y_vrr = PAZ*I_KINETIC_G3xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F3y_vrr;
    Double I_KINETIC_H3x2z_F3y_vrr = PAZ*I_KINETIC_G3xz_F3y_vrr+oned2z*I_KINETIC_F3x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_KINETIC_H2x3y_F3y_vrr = PAX*I_KINETIC_Gx3y_F3y_vrr+oned2z*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F3y_vrr-bdz*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_KINETIC_H2x2yz_F3y_vrr = PAZ*I_KINETIC_G2x2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F3y_vrr;
    Double I_KINETIC_H2xy2z_F3y_vrr = PAY*I_KINETIC_G2x2z_F3y_vrr+3*oned2z*I_KINETIC_G2x2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F3y_vrr;
    Double I_KINETIC_H2x3z_F3y_vrr = PAX*I_KINETIC_Gx3z_F3y_vrr+oned2z*I_KINETIC_F3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F3y_vrr-bdz*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_KINETIC_Hx4y_F3y_vrr = PAX*I_KINETIC_G4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F3y_vrr;
    Double I_KINETIC_Hx3yz_F3y_vrr = PAZ*I_KINETIC_Gx3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F3y_vrr;
    Double I_KINETIC_Hx2y2z_F3y_vrr = PAX*I_KINETIC_G2y2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F3y_vrr;
    Double I_KINETIC_Hxy3z_F3y_vrr = PAY*I_KINETIC_Gx3z_F3y_vrr+3*oned2z*I_KINETIC_Gx3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F3y_vrr;
    Double I_KINETIC_Hx4z_F3y_vrr = PAX*I_KINETIC_G4z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F3y_vrr;
    Double I_KINETIC_H5y_F3y_vrr = PAY*I_KINETIC_G4y_F3y_vrr+4*oned2z*I_KINETIC_F3y_F3y_vrr+3*oned2z*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F3y_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_KINETIC_H4yz_F3y_vrr = PAZ*I_KINETIC_G4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F3y_vrr;
    Double I_KINETIC_H3y2z_F3y_vrr = PAZ*I_KINETIC_G3yz_F3y_vrr+oned2z*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_KINETIC_H2y3z_F3y_vrr = PAY*I_KINETIC_Gy3z_F3y_vrr+oned2z*I_KINETIC_F3z_F3y_vrr+3*oned2z*I_KINETIC_Gy3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F3y_vrr-bdz*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_KINETIC_Hy4z_F3y_vrr = PAY*I_KINETIC_G4z_F3y_vrr+3*oned2z*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F3y_vrr;
    Double I_KINETIC_H5z_F3y_vrr = PAZ*I_KINETIC_G4z_F3y_vrr+4*oned2z*I_KINETIC_F3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F3y_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_KINETIC_H5x_F2yz_vrr = PAX*I_KINETIC_G4x_F2yz_vrr+4*oned2z*I_KINETIC_F3x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F2yz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_F2yz_vrr;
    Double I_KINETIC_H4xy_F2yz_vrr = PAY*I_KINETIC_G4x_F2yz_vrr+2*oned2z*I_KINETIC_G4x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F2yz_vrr;
    Double I_KINETIC_H4xz_F2yz_vrr = PAZ*I_KINETIC_G4x_F2yz_vrr+oned2z*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F2yz_vrr;
    Double I_KINETIC_H3x2y_F2yz_vrr = PAY*I_KINETIC_G3xy_F2yz_vrr+oned2z*I_KINETIC_F3x_F2yz_vrr+2*oned2z*I_KINETIC_G3xy_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F2yz_vrr-bdz*I_TWOBODYOVERLAP_F3x_F2yz_vrr;
    Double I_KINETIC_H3xyz_F2yz_vrr = PAZ*I_KINETIC_G3xy_F2yz_vrr+oned2z*I_KINETIC_G3xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F2yz_vrr;
    Double I_KINETIC_H3x2z_F2yz_vrr = PAZ*I_KINETIC_G3xz_F2yz_vrr+oned2z*I_KINETIC_F3x_F2yz_vrr+oned2z*I_KINETIC_G3xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_F3x_F2yz_vrr;
    Double I_KINETIC_H2x3y_F2yz_vrr = PAX*I_KINETIC_Gx3y_F2yz_vrr+oned2z*I_KINETIC_F3y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F2yz_vrr-bdz*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_KINETIC_H2x2yz_F2yz_vrr = PAZ*I_KINETIC_G2x2y_F2yz_vrr+oned2z*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F2yz_vrr;
    Double I_KINETIC_H2xy2z_F2yz_vrr = PAY*I_KINETIC_G2x2z_F2yz_vrr+2*oned2z*I_KINETIC_G2x2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F2yz_vrr;
    Double I_KINETIC_H2x3z_F2yz_vrr = PAX*I_KINETIC_Gx3z_F2yz_vrr+oned2z*I_KINETIC_F3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_KINETIC_Hx4y_F2yz_vrr = PAX*I_KINETIC_G4y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F2yz_vrr;
    Double I_KINETIC_Hx3yz_F2yz_vrr = PAZ*I_KINETIC_Gx3y_F2yz_vrr+oned2z*I_KINETIC_Gx3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F2yz_vrr;
    Double I_KINETIC_Hx2y2z_F2yz_vrr = PAX*I_KINETIC_G2y2z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F2yz_vrr;
    Double I_KINETIC_Hxy3z_F2yz_vrr = PAY*I_KINETIC_Gx3z_F2yz_vrr+2*oned2z*I_KINETIC_Gx3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F2yz_vrr;
    Double I_KINETIC_Hx4z_F2yz_vrr = PAX*I_KINETIC_G4z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F2yz_vrr;
    Double I_KINETIC_H5y_F2yz_vrr = PAY*I_KINETIC_G4y_F2yz_vrr+4*oned2z*I_KINETIC_F3y_F2yz_vrr+2*oned2z*I_KINETIC_G4y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F2yz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_KINETIC_H4yz_F2yz_vrr = PAZ*I_KINETIC_G4y_F2yz_vrr+oned2z*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F2yz_vrr;
    Double I_KINETIC_H3y2z_F2yz_vrr = PAZ*I_KINETIC_G3yz_F2yz_vrr+oned2z*I_KINETIC_F3y_F2yz_vrr+oned2z*I_KINETIC_G3yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_KINETIC_H2y3z_F2yz_vrr = PAY*I_KINETIC_Gy3z_F2yz_vrr+oned2z*I_KINETIC_F3z_F2yz_vrr+2*oned2z*I_KINETIC_Gy3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_KINETIC_Hy4z_F2yz_vrr = PAY*I_KINETIC_G4z_F2yz_vrr+2*oned2z*I_KINETIC_G4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F2yz_vrr;
    Double I_KINETIC_H5z_F2yz_vrr = PAZ*I_KINETIC_G4z_F2yz_vrr+4*oned2z*I_KINETIC_F3z_F2yz_vrr+oned2z*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F2yz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_KINETIC_H5x_Fy2z_vrr = PAX*I_KINETIC_G4x_Fy2z_vrr+4*oned2z*I_KINETIC_F3x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fy2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Fy2z_vrr;
    Double I_KINETIC_H4xy_Fy2z_vrr = PAY*I_KINETIC_G4x_Fy2z_vrr+oned2z*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fy2z_vrr;
    Double I_KINETIC_H4xz_Fy2z_vrr = PAZ*I_KINETIC_G4x_Fy2z_vrr+2*oned2z*I_KINETIC_G4x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fy2z_vrr;
    Double I_KINETIC_H3x2y_Fy2z_vrr = PAY*I_KINETIC_G3xy_Fy2z_vrr+oned2z*I_KINETIC_F3x_Fy2z_vrr+oned2z*I_KINETIC_G3xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Fy2z_vrr;
    Double I_KINETIC_H3xyz_Fy2z_vrr = PAZ*I_KINETIC_G3xy_Fy2z_vrr+2*oned2z*I_KINETIC_G3xy_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Fy2z_vrr;
    Double I_KINETIC_H3x2z_Fy2z_vrr = PAZ*I_KINETIC_G3xz_Fy2z_vrr+oned2z*I_KINETIC_F3x_Fy2z_vrr+2*oned2z*I_KINETIC_G3xz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Fy2z_vrr;
    Double I_KINETIC_H2x3y_Fy2z_vrr = PAX*I_KINETIC_Gx3y_Fy2z_vrr+oned2z*I_KINETIC_F3y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_KINETIC_H2x2yz_Fy2z_vrr = PAZ*I_KINETIC_G2x2y_Fy2z_vrr+2*oned2z*I_KINETIC_G2x2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fy2z_vrr;
    Double I_KINETIC_H2xy2z_Fy2z_vrr = PAY*I_KINETIC_G2x2z_Fy2z_vrr+oned2z*I_KINETIC_G2x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Fy2z_vrr;
    Double I_KINETIC_H2x3z_Fy2z_vrr = PAX*I_KINETIC_Gx3z_Fy2z_vrr+oned2z*I_KINETIC_F3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_KINETIC_Hx4y_Fy2z_vrr = PAX*I_KINETIC_G4y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fy2z_vrr;
    Double I_KINETIC_Hx3yz_Fy2z_vrr = PAZ*I_KINETIC_Gx3y_Fy2z_vrr+2*oned2z*I_KINETIC_Gx3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Fy2z_vrr;
    Double I_KINETIC_Hx2y2z_Fy2z_vrr = PAX*I_KINETIC_G2y2z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Fy2z_vrr;
    Double I_KINETIC_Hxy3z_Fy2z_vrr = PAY*I_KINETIC_Gx3z_Fy2z_vrr+oned2z*I_KINETIC_Gx3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Fy2z_vrr;
    Double I_KINETIC_Hx4z_Fy2z_vrr = PAX*I_KINETIC_G4z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fy2z_vrr;
    Double I_KINETIC_H5y_Fy2z_vrr = PAY*I_KINETIC_G4y_Fy2z_vrr+4*oned2z*I_KINETIC_F3y_Fy2z_vrr+oned2z*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fy2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_KINETIC_H4yz_Fy2z_vrr = PAZ*I_KINETIC_G4y_Fy2z_vrr+2*oned2z*I_KINETIC_G4y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fy2z_vrr;
    Double I_KINETIC_H3y2z_Fy2z_vrr = PAZ*I_KINETIC_G3yz_Fy2z_vrr+oned2z*I_KINETIC_F3y_Fy2z_vrr+2*oned2z*I_KINETIC_G3yz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_KINETIC_H2y3z_Fy2z_vrr = PAY*I_KINETIC_Gy3z_Fy2z_vrr+oned2z*I_KINETIC_F3z_Fy2z_vrr+oned2z*I_KINETIC_Gy3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_KINETIC_Hy4z_Fy2z_vrr = PAY*I_KINETIC_G4z_Fy2z_vrr+oned2z*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fy2z_vrr;
    Double I_KINETIC_H5z_Fy2z_vrr = PAZ*I_KINETIC_G4z_Fy2z_vrr+4*oned2z*I_KINETIC_F3z_Fy2z_vrr+2*oned2z*I_KINETIC_G4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fy2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_KINETIC_H5x_F3z_vrr = PAX*I_KINETIC_G4x_F3z_vrr+4*oned2z*I_KINETIC_F3x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F3z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_KINETIC_H4xy_F3z_vrr = PAY*I_KINETIC_G4x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F3z_vrr;
    Double I_KINETIC_H4xz_F3z_vrr = PAZ*I_KINETIC_G4x_F3z_vrr+3*oned2z*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F3z_vrr;
    Double I_KINETIC_H3x2y_F3z_vrr = PAY*I_KINETIC_G3xy_F3z_vrr+oned2z*I_KINETIC_F3x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F3z_vrr-bdz*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_KINETIC_H3xyz_F3z_vrr = PAZ*I_KINETIC_G3xy_F3z_vrr+3*oned2z*I_KINETIC_G3xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F3z_vrr;
    Double I_KINETIC_H3x2z_F3z_vrr = PAZ*I_KINETIC_G3xz_F3z_vrr+oned2z*I_KINETIC_F3x_F3z_vrr+3*oned2z*I_KINETIC_G3xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_KINETIC_H2x3y_F3z_vrr = PAX*I_KINETIC_Gx3y_F3z_vrr+oned2z*I_KINETIC_F3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F3z_vrr-bdz*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_KINETIC_H2x2yz_F3z_vrr = PAZ*I_KINETIC_G2x2y_F3z_vrr+3*oned2z*I_KINETIC_G2x2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F3z_vrr;
    Double I_KINETIC_H2xy2z_F3z_vrr = PAY*I_KINETIC_G2x2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F3z_vrr;
    Double I_KINETIC_H2x3z_F3z_vrr = PAX*I_KINETIC_Gx3z_F3z_vrr+oned2z*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F3z_vrr-bdz*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_KINETIC_Hx4y_F3z_vrr = PAX*I_KINETIC_G4y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F3z_vrr;
    Double I_KINETIC_Hx3yz_F3z_vrr = PAZ*I_KINETIC_Gx3y_F3z_vrr+3*oned2z*I_KINETIC_Gx3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F3z_vrr;
    Double I_KINETIC_Hx2y2z_F3z_vrr = PAX*I_KINETIC_G2y2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F3z_vrr;
    Double I_KINETIC_Hxy3z_F3z_vrr = PAY*I_KINETIC_Gx3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F3z_vrr;
    Double I_KINETIC_Hx4z_F3z_vrr = PAX*I_KINETIC_G4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F3z_vrr;
    Double I_KINETIC_H5y_F3z_vrr = PAY*I_KINETIC_G4y_F3z_vrr+4*oned2z*I_KINETIC_F3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F3z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_KINETIC_H4yz_F3z_vrr = PAZ*I_KINETIC_G4y_F3z_vrr+3*oned2z*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F3z_vrr;
    Double I_KINETIC_H3y2z_F3z_vrr = PAZ*I_KINETIC_G3yz_F3z_vrr+oned2z*I_KINETIC_F3y_F3z_vrr+3*oned2z*I_KINETIC_G3yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_KINETIC_H2y3z_F3z_vrr = PAY*I_KINETIC_Gy3z_F3z_vrr+oned2z*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F3z_vrr-bdz*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_KINETIC_Hy4z_F3z_vrr = PAY*I_KINETIC_G4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F3z_vrr;
    Double I_KINETIC_H5z_F3z_vrr = PAZ*I_KINETIC_G4z_F3z_vrr+4*oned2z*I_KINETIC_F3z_F3z_vrr+3*oned2z*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F3z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_F_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_H_F_aa_coefs = alpha*alpha;
    I_KINETIC_H5x_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_F3x_vrr;
    I_KINETIC_H4xy_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_F3x_vrr;
    I_KINETIC_H4xz_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_F3x_vrr;
    I_KINETIC_H3x2y_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_F3x_vrr;
    I_KINETIC_H3xyz_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_F3x_vrr;
    I_KINETIC_H3x2z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_F3x_vrr;
    I_KINETIC_H2x3y_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_F3x_vrr;
    I_KINETIC_H2x2yz_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_F3x_vrr;
    I_KINETIC_H2xy2z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_F3x_vrr;
    I_KINETIC_H2x3z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_F3x_vrr;
    I_KINETIC_Hx4y_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_F3x_vrr;
    I_KINETIC_Hx3yz_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_F3x_vrr;
    I_KINETIC_Hx2y2z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_F3x_vrr;
    I_KINETIC_Hxy3z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_F3x_vrr;
    I_KINETIC_Hx4z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_F3x_vrr;
    I_KINETIC_H5y_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_F3x_vrr;
    I_KINETIC_H4yz_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_F3x_vrr;
    I_KINETIC_H3y2z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_F3x_vrr;
    I_KINETIC_H2y3z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_F3x_vrr;
    I_KINETIC_Hy4z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_F3x_vrr;
    I_KINETIC_H5z_F3x_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_F3x_vrr;
    I_KINETIC_H5x_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_F2xy_vrr;
    I_KINETIC_H4xy_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_F2xy_vrr;
    I_KINETIC_H4xz_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_F2xy_vrr;
    I_KINETIC_H3x2y_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_F2xy_vrr;
    I_KINETIC_H3xyz_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_F2xy_vrr;
    I_KINETIC_H3x2z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_F2xy_vrr;
    I_KINETIC_H2x3y_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_F2xy_vrr;
    I_KINETIC_H2x2yz_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_F2xy_vrr;
    I_KINETIC_H2xy2z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_F2xy_vrr;
    I_KINETIC_H2x3z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_F2xy_vrr;
    I_KINETIC_Hx4y_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_F2xy_vrr;
    I_KINETIC_Hx3yz_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_F2xy_vrr;
    I_KINETIC_Hx2y2z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_F2xy_vrr;
    I_KINETIC_Hxy3z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_F2xy_vrr;
    I_KINETIC_Hx4z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_F2xy_vrr;
    I_KINETIC_H5y_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_F2xy_vrr;
    I_KINETIC_H4yz_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_F2xy_vrr;
    I_KINETIC_H3y2z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_F2xy_vrr;
    I_KINETIC_H2y3z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_F2xy_vrr;
    I_KINETIC_Hy4z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_F2xy_vrr;
    I_KINETIC_H5z_F2xy_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_F2xy_vrr;
    I_KINETIC_H5x_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_F2xz_vrr;
    I_KINETIC_H4xy_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_F2xz_vrr;
    I_KINETIC_H4xz_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_F2xz_vrr;
    I_KINETIC_H3x2y_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_F2xz_vrr;
    I_KINETIC_H3xyz_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_F2xz_vrr;
    I_KINETIC_H3x2z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_F2xz_vrr;
    I_KINETIC_H2x3y_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_F2xz_vrr;
    I_KINETIC_H2x2yz_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_F2xz_vrr;
    I_KINETIC_H2xy2z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_F2xz_vrr;
    I_KINETIC_H2x3z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_F2xz_vrr;
    I_KINETIC_Hx4y_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_F2xz_vrr;
    I_KINETIC_Hx3yz_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_F2xz_vrr;
    I_KINETIC_Hx2y2z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_F2xz_vrr;
    I_KINETIC_Hxy3z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_F2xz_vrr;
    I_KINETIC_Hx4z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_F2xz_vrr;
    I_KINETIC_H5y_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_F2xz_vrr;
    I_KINETIC_H4yz_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_F2xz_vrr;
    I_KINETIC_H3y2z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_F2xz_vrr;
    I_KINETIC_H2y3z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_F2xz_vrr;
    I_KINETIC_Hy4z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_F2xz_vrr;
    I_KINETIC_H5z_F2xz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_F2xz_vrr;
    I_KINETIC_H5x_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_Fx2y_vrr;
    I_KINETIC_H4xy_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_Fx2y_vrr;
    I_KINETIC_H4xz_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_Fx2y_vrr;
    I_KINETIC_H3x2y_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_Fx2y_vrr;
    I_KINETIC_H3xyz_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_Fx2y_vrr;
    I_KINETIC_H3x2z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_Fx2y_vrr;
    I_KINETIC_H2x3y_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_Fx2y_vrr;
    I_KINETIC_H2x2yz_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_Fx2y_vrr;
    I_KINETIC_H2xy2z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_Fx2y_vrr;
    I_KINETIC_H2x3z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_Fx2y_vrr;
    I_KINETIC_Hx4y_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_Fx2y_vrr;
    I_KINETIC_Hx3yz_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_Fx2y_vrr;
    I_KINETIC_Hx2y2z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_Fx2y_vrr;
    I_KINETIC_Hxy3z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_Fx2y_vrr;
    I_KINETIC_Hx4z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_Fx2y_vrr;
    I_KINETIC_H5y_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_Fx2y_vrr;
    I_KINETIC_H4yz_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_Fx2y_vrr;
    I_KINETIC_H3y2z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_Fx2y_vrr;
    I_KINETIC_H2y3z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_Fx2y_vrr;
    I_KINETIC_Hy4z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_Fx2y_vrr;
    I_KINETIC_H5z_Fx2y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_Fx2y_vrr;
    I_KINETIC_H5x_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_Fxyz_vrr;
    I_KINETIC_H4xy_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_Fxyz_vrr;
    I_KINETIC_H4xz_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_Fxyz_vrr;
    I_KINETIC_H3x2y_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_Fxyz_vrr;
    I_KINETIC_H3xyz_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_Fxyz_vrr;
    I_KINETIC_H3x2z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_Fxyz_vrr;
    I_KINETIC_H2x3y_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_Fxyz_vrr;
    I_KINETIC_H2x2yz_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_Fxyz_vrr;
    I_KINETIC_H2xy2z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_Fxyz_vrr;
    I_KINETIC_H2x3z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_Fxyz_vrr;
    I_KINETIC_Hx4y_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_Fxyz_vrr;
    I_KINETIC_Hx3yz_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_Fxyz_vrr;
    I_KINETIC_Hx2y2z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_Fxyz_vrr;
    I_KINETIC_Hxy3z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_Fxyz_vrr;
    I_KINETIC_Hx4z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_Fxyz_vrr;
    I_KINETIC_H5y_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_Fxyz_vrr;
    I_KINETIC_H4yz_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_Fxyz_vrr;
    I_KINETIC_H3y2z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_Fxyz_vrr;
    I_KINETIC_H2y3z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_Fxyz_vrr;
    I_KINETIC_Hy4z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_Fxyz_vrr;
    I_KINETIC_H5z_Fxyz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_Fxyz_vrr;
    I_KINETIC_H5x_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_Fx2z_vrr;
    I_KINETIC_H4xy_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_Fx2z_vrr;
    I_KINETIC_H4xz_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_Fx2z_vrr;
    I_KINETIC_H3x2y_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_Fx2z_vrr;
    I_KINETIC_H3xyz_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_Fx2z_vrr;
    I_KINETIC_H3x2z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_Fx2z_vrr;
    I_KINETIC_H2x3y_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_Fx2z_vrr;
    I_KINETIC_H2x2yz_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_Fx2z_vrr;
    I_KINETIC_H2xy2z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_Fx2z_vrr;
    I_KINETIC_H2x3z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_Fx2z_vrr;
    I_KINETIC_Hx4y_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_Fx2z_vrr;
    I_KINETIC_Hx3yz_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_Fx2z_vrr;
    I_KINETIC_Hx2y2z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_Fx2z_vrr;
    I_KINETIC_Hxy3z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_Fx2z_vrr;
    I_KINETIC_Hx4z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_Fx2z_vrr;
    I_KINETIC_H5y_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_Fx2z_vrr;
    I_KINETIC_H4yz_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_Fx2z_vrr;
    I_KINETIC_H3y2z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_Fx2z_vrr;
    I_KINETIC_H2y3z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_Fx2z_vrr;
    I_KINETIC_Hy4z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_Fx2z_vrr;
    I_KINETIC_H5z_Fx2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_Fx2z_vrr;
    I_KINETIC_H5x_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_F3y_vrr;
    I_KINETIC_H4xy_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_F3y_vrr;
    I_KINETIC_H4xz_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_F3y_vrr;
    I_KINETIC_H3x2y_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_F3y_vrr;
    I_KINETIC_H3xyz_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_F3y_vrr;
    I_KINETIC_H3x2z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_F3y_vrr;
    I_KINETIC_H2x3y_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_F3y_vrr;
    I_KINETIC_H2x2yz_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_F3y_vrr;
    I_KINETIC_H2xy2z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_F3y_vrr;
    I_KINETIC_H2x3z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_F3y_vrr;
    I_KINETIC_Hx4y_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_F3y_vrr;
    I_KINETIC_Hx3yz_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_F3y_vrr;
    I_KINETIC_Hx2y2z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_F3y_vrr;
    I_KINETIC_Hxy3z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_F3y_vrr;
    I_KINETIC_Hx4z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_F3y_vrr;
    I_KINETIC_H5y_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_F3y_vrr;
    I_KINETIC_H4yz_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_F3y_vrr;
    I_KINETIC_H3y2z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_F3y_vrr;
    I_KINETIC_H2y3z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_F3y_vrr;
    I_KINETIC_Hy4z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_F3y_vrr;
    I_KINETIC_H5z_F3y_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_F3y_vrr;
    I_KINETIC_H5x_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_F2yz_vrr;
    I_KINETIC_H4xy_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_F2yz_vrr;
    I_KINETIC_H4xz_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_F2yz_vrr;
    I_KINETIC_H3x2y_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_F2yz_vrr;
    I_KINETIC_H3xyz_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_F2yz_vrr;
    I_KINETIC_H3x2z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_F2yz_vrr;
    I_KINETIC_H2x3y_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_F2yz_vrr;
    I_KINETIC_H2x2yz_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_F2yz_vrr;
    I_KINETIC_H2xy2z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_F2yz_vrr;
    I_KINETIC_H2x3z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_F2yz_vrr;
    I_KINETIC_Hx4y_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_F2yz_vrr;
    I_KINETIC_Hx3yz_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_F2yz_vrr;
    I_KINETIC_Hx2y2z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_F2yz_vrr;
    I_KINETIC_Hxy3z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_F2yz_vrr;
    I_KINETIC_Hx4z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_F2yz_vrr;
    I_KINETIC_H5y_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_F2yz_vrr;
    I_KINETIC_H4yz_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_F2yz_vrr;
    I_KINETIC_H3y2z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_F2yz_vrr;
    I_KINETIC_H2y3z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_F2yz_vrr;
    I_KINETIC_Hy4z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_F2yz_vrr;
    I_KINETIC_H5z_F2yz_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_F2yz_vrr;
    I_KINETIC_H5x_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_Fy2z_vrr;
    I_KINETIC_H4xy_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_Fy2z_vrr;
    I_KINETIC_H4xz_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_Fy2z_vrr;
    I_KINETIC_H3x2y_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_Fy2z_vrr;
    I_KINETIC_H3xyz_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_Fy2z_vrr;
    I_KINETIC_H3x2z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_Fy2z_vrr;
    I_KINETIC_H2x3y_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_Fy2z_vrr;
    I_KINETIC_H2x2yz_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_Fy2z_vrr;
    I_KINETIC_H2xy2z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_Fy2z_vrr;
    I_KINETIC_H2x3z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_Fy2z_vrr;
    I_KINETIC_Hx4y_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_Fy2z_vrr;
    I_KINETIC_Hx3yz_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_Fy2z_vrr;
    I_KINETIC_Hx2y2z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_Fy2z_vrr;
    I_KINETIC_Hxy3z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_Fy2z_vrr;
    I_KINETIC_Hx4z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_Fy2z_vrr;
    I_KINETIC_H5y_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_Fy2z_vrr;
    I_KINETIC_H4yz_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_Fy2z_vrr;
    I_KINETIC_H3y2z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_Fy2z_vrr;
    I_KINETIC_H2y3z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_Fy2z_vrr;
    I_KINETIC_Hy4z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_Fy2z_vrr;
    I_KINETIC_H5z_Fy2z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_Fy2z_vrr;
    I_KINETIC_H5x_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5x_F3z_vrr;
    I_KINETIC_H4xy_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xy_F3z_vrr;
    I_KINETIC_H4xz_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4xz_F3z_vrr;
    I_KINETIC_H3x2y_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2y_F3z_vrr;
    I_KINETIC_H3xyz_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3xyz_F3z_vrr;
    I_KINETIC_H3x2z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3x2z_F3z_vrr;
    I_KINETIC_H2x3y_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3y_F3z_vrr;
    I_KINETIC_H2x2yz_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x2yz_F3z_vrr;
    I_KINETIC_H2xy2z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2xy2z_F3z_vrr;
    I_KINETIC_H2x3z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2x3z_F3z_vrr;
    I_KINETIC_Hx4y_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4y_F3z_vrr;
    I_KINETIC_Hx3yz_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx3yz_F3z_vrr;
    I_KINETIC_Hx2y2z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx2y2z_F3z_vrr;
    I_KINETIC_Hxy3z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hxy3z_F3z_vrr;
    I_KINETIC_Hx4z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hx4z_F3z_vrr;
    I_KINETIC_H5y_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5y_F3z_vrr;
    I_KINETIC_H4yz_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H4yz_F3z_vrr;
    I_KINETIC_H3y2z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H3y2z_F3z_vrr;
    I_KINETIC_H2y3z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H2y3z_F3z_vrr;
    I_KINETIC_Hy4z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_Hy4z_F3z_vrr;
    I_KINETIC_H5z_F3z_aa += SQ_KINETIC_H_F_aa_coefs*I_KINETIC_H5z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_F_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_F_F_a_coefs = alpha;
    I_KINETIC_F3x_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_F3x_vrr;
    I_KINETIC_F2xy_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_F3x_vrr;
    I_KINETIC_F2xz_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_F3x_vrr;
    I_KINETIC_Fx2y_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_F3x_vrr;
    I_KINETIC_Fxyz_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_F3x_vrr;
    I_KINETIC_Fx2z_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_F3x_vrr;
    I_KINETIC_F3y_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_F3x_vrr;
    I_KINETIC_F2yz_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_F3x_vrr;
    I_KINETIC_Fy2z_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_F3x_vrr;
    I_KINETIC_F3z_F3x_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_F3x_vrr;
    I_KINETIC_F3x_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_F2xy_vrr;
    I_KINETIC_F2xy_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_F2xy_vrr;
    I_KINETIC_F2xz_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_F2xy_vrr;
    I_KINETIC_Fx2y_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_F2xy_vrr;
    I_KINETIC_Fxyz_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_F2xy_vrr;
    I_KINETIC_Fx2z_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_F2xy_vrr;
    I_KINETIC_F3y_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_F2xy_vrr;
    I_KINETIC_F2yz_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_F2xy_vrr;
    I_KINETIC_Fy2z_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_F2xy_vrr;
    I_KINETIC_F3z_F2xy_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_F2xy_vrr;
    I_KINETIC_F3x_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_F2xz_vrr;
    I_KINETIC_F2xy_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_F2xz_vrr;
    I_KINETIC_F2xz_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_F2xz_vrr;
    I_KINETIC_Fx2y_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_F2xz_vrr;
    I_KINETIC_Fxyz_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_F2xz_vrr;
    I_KINETIC_Fx2z_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_F2xz_vrr;
    I_KINETIC_F3y_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_F2xz_vrr;
    I_KINETIC_F2yz_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_F2xz_vrr;
    I_KINETIC_Fy2z_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_F2xz_vrr;
    I_KINETIC_F3z_F2xz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_F2xz_vrr;
    I_KINETIC_F3x_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_Fx2y_vrr;
    I_KINETIC_F2xy_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_Fx2y_vrr;
    I_KINETIC_F2xz_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_Fx2y_vrr;
    I_KINETIC_Fx2y_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_Fx2y_vrr;
    I_KINETIC_Fxyz_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_Fx2y_vrr;
    I_KINETIC_Fx2z_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_Fx2y_vrr;
    I_KINETIC_F3y_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_Fx2y_vrr;
    I_KINETIC_F2yz_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_Fx2y_vrr;
    I_KINETIC_Fy2z_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_Fx2y_vrr;
    I_KINETIC_F3z_Fx2y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_Fx2y_vrr;
    I_KINETIC_F3x_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_Fxyz_vrr;
    I_KINETIC_F2xy_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_Fxyz_vrr;
    I_KINETIC_F2xz_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_Fxyz_vrr;
    I_KINETIC_Fx2y_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_Fxyz_vrr;
    I_KINETIC_Fxyz_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_Fxyz_vrr;
    I_KINETIC_Fx2z_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_Fxyz_vrr;
    I_KINETIC_F3y_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_Fxyz_vrr;
    I_KINETIC_F2yz_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_Fxyz_vrr;
    I_KINETIC_Fy2z_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_Fxyz_vrr;
    I_KINETIC_F3z_Fxyz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_Fxyz_vrr;
    I_KINETIC_F3x_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_Fx2z_vrr;
    I_KINETIC_F2xy_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_Fx2z_vrr;
    I_KINETIC_F2xz_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_Fx2z_vrr;
    I_KINETIC_Fx2y_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_Fx2z_vrr;
    I_KINETIC_Fxyz_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_Fx2z_vrr;
    I_KINETIC_Fx2z_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_Fx2z_vrr;
    I_KINETIC_F3y_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_Fx2z_vrr;
    I_KINETIC_F2yz_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_Fx2z_vrr;
    I_KINETIC_Fy2z_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_Fx2z_vrr;
    I_KINETIC_F3z_Fx2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_Fx2z_vrr;
    I_KINETIC_F3x_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_F3y_vrr;
    I_KINETIC_F2xy_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_F3y_vrr;
    I_KINETIC_F2xz_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_F3y_vrr;
    I_KINETIC_Fx2y_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_F3y_vrr;
    I_KINETIC_Fxyz_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_F3y_vrr;
    I_KINETIC_Fx2z_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_F3y_vrr;
    I_KINETIC_F3y_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_F3y_vrr;
    I_KINETIC_F2yz_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_F3y_vrr;
    I_KINETIC_Fy2z_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_F3y_vrr;
    I_KINETIC_F3z_F3y_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_F3y_vrr;
    I_KINETIC_F3x_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_F2yz_vrr;
    I_KINETIC_F2xy_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_F2yz_vrr;
    I_KINETIC_F2xz_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_F2yz_vrr;
    I_KINETIC_Fx2y_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_F2yz_vrr;
    I_KINETIC_Fxyz_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_F2yz_vrr;
    I_KINETIC_Fx2z_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_F2yz_vrr;
    I_KINETIC_F3y_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_F2yz_vrr;
    I_KINETIC_F2yz_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_F2yz_vrr;
    I_KINETIC_Fy2z_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_F2yz_vrr;
    I_KINETIC_F3z_F2yz_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_F2yz_vrr;
    I_KINETIC_F3x_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_Fy2z_vrr;
    I_KINETIC_F2xy_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_Fy2z_vrr;
    I_KINETIC_F2xz_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_Fy2z_vrr;
    I_KINETIC_Fx2y_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_Fy2z_vrr;
    I_KINETIC_Fxyz_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_Fy2z_vrr;
    I_KINETIC_Fx2z_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_Fy2z_vrr;
    I_KINETIC_F3y_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_Fy2z_vrr;
    I_KINETIC_F2yz_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_Fy2z_vrr;
    I_KINETIC_Fy2z_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_Fy2z_vrr;
    I_KINETIC_F3z_Fy2z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_Fy2z_vrr;
    I_KINETIC_F3x_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3x_F3z_vrr;
    I_KINETIC_F2xy_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xy_F3z_vrr;
    I_KINETIC_F2xz_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2xz_F3z_vrr;
    I_KINETIC_Fx2y_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2y_F3z_vrr;
    I_KINETIC_Fxyz_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fxyz_F3z_vrr;
    I_KINETIC_Fx2z_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fx2z_F3z_vrr;
    I_KINETIC_F3y_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3y_F3z_vrr;
    I_KINETIC_F2yz_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F2yz_F3z_vrr;
    I_KINETIC_Fy2z_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_Fy2z_F3z_vrr;
    I_KINETIC_F3z_F3z_a += SQ_KINETIC_F_F_a_coefs*I_KINETIC_F3z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_F
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_Px_F3x += I_KINETIC_Px_F3x_vrr;
    I_KINETIC_Py_F3x += I_KINETIC_Py_F3x_vrr;
    I_KINETIC_Pz_F3x += I_KINETIC_Pz_F3x_vrr;
    I_KINETIC_Px_F2xy += I_KINETIC_Px_F2xy_vrr;
    I_KINETIC_Py_F2xy += I_KINETIC_Py_F2xy_vrr;
    I_KINETIC_Pz_F2xy += I_KINETIC_Pz_F2xy_vrr;
    I_KINETIC_Px_F2xz += I_KINETIC_Px_F2xz_vrr;
    I_KINETIC_Py_F2xz += I_KINETIC_Py_F2xz_vrr;
    I_KINETIC_Pz_F2xz += I_KINETIC_Pz_F2xz_vrr;
    I_KINETIC_Px_Fx2y += I_KINETIC_Px_Fx2y_vrr;
    I_KINETIC_Py_Fx2y += I_KINETIC_Py_Fx2y_vrr;
    I_KINETIC_Pz_Fx2y += I_KINETIC_Pz_Fx2y_vrr;
    I_KINETIC_Px_Fxyz += I_KINETIC_Px_Fxyz_vrr;
    I_KINETIC_Py_Fxyz += I_KINETIC_Py_Fxyz_vrr;
    I_KINETIC_Pz_Fxyz += I_KINETIC_Pz_Fxyz_vrr;
    I_KINETIC_Px_Fx2z += I_KINETIC_Px_Fx2z_vrr;
    I_KINETIC_Py_Fx2z += I_KINETIC_Py_Fx2z_vrr;
    I_KINETIC_Pz_Fx2z += I_KINETIC_Pz_Fx2z_vrr;
    I_KINETIC_Px_F3y += I_KINETIC_Px_F3y_vrr;
    I_KINETIC_Py_F3y += I_KINETIC_Py_F3y_vrr;
    I_KINETIC_Pz_F3y += I_KINETIC_Pz_F3y_vrr;
    I_KINETIC_Px_F2yz += I_KINETIC_Px_F2yz_vrr;
    I_KINETIC_Py_F2yz += I_KINETIC_Py_F2yz_vrr;
    I_KINETIC_Pz_F2yz += I_KINETIC_Pz_F2yz_vrr;
    I_KINETIC_Px_Fy2z += I_KINETIC_Px_Fy2z_vrr;
    I_KINETIC_Py_Fy2z += I_KINETIC_Py_Fy2z_vrr;
    I_KINETIC_Pz_Fy2z += I_KINETIC_Pz_Fy2z_vrr;
    I_KINETIC_Px_F3z += I_KINETIC_Px_F3z_vrr;
    I_KINETIC_Py_F3z += I_KINETIC_Py_F3z_vrr;
    I_KINETIC_Pz_F3z += I_KINETIC_Pz_F3z_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_F_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_F_aa
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_P_F
   ************************************************************/
  abcd[0] = 4.0E0*I_KINETIC_H5x_F3x_aa-2.0E0*3*I_KINETIC_F3x_F3x_a-2.0E0*4*I_KINETIC_F3x_F3x_a+3*2*I_KINETIC_Px_F3x;
  abcd[1] = 4.0E0*I_KINETIC_H4xy_F3x_aa-2.0E0*2*I_KINETIC_F2xy_F3x_a-2.0E0*3*I_KINETIC_F2xy_F3x_a+2*1*I_KINETIC_Py_F3x;
  abcd[2] = 4.0E0*I_KINETIC_H4xz_F3x_aa-2.0E0*2*I_KINETIC_F2xz_F3x_a-2.0E0*3*I_KINETIC_F2xz_F3x_a+2*1*I_KINETIC_Pz_F3x;
  abcd[3] = 4.0E0*I_KINETIC_H3x2y_F3x_aa-2.0E0*1*I_KINETIC_Fx2y_F3x_a-2.0E0*2*I_KINETIC_Fx2y_F3x_a;
  abcd[4] = 4.0E0*I_KINETIC_H3xyz_F3x_aa-2.0E0*1*I_KINETIC_Fxyz_F3x_a-2.0E0*2*I_KINETIC_Fxyz_F3x_a;
  abcd[5] = 4.0E0*I_KINETIC_H3x2z_F3x_aa-2.0E0*1*I_KINETIC_Fx2z_F3x_a-2.0E0*2*I_KINETIC_Fx2z_F3x_a;
  abcd[6] = 4.0E0*I_KINETIC_H2x3y_F3x_aa-2.0E0*1*I_KINETIC_F3y_F3x_a;
  abcd[7] = 4.0E0*I_KINETIC_H2x2yz_F3x_aa-2.0E0*1*I_KINETIC_F2yz_F3x_a;
  abcd[8] = 4.0E0*I_KINETIC_H2xy2z_F3x_aa-2.0E0*1*I_KINETIC_Fy2z_F3x_a;
  abcd[9] = 4.0E0*I_KINETIC_H2x3z_F3x_aa-2.0E0*1*I_KINETIC_F3z_F3x_a;
  abcd[10] = 4.0E0*I_KINETIC_H5x_F2xy_aa-2.0E0*3*I_KINETIC_F3x_F2xy_a-2.0E0*4*I_KINETIC_F3x_F2xy_a+3*2*I_KINETIC_Px_F2xy;
  abcd[11] = 4.0E0*I_KINETIC_H4xy_F2xy_aa-2.0E0*2*I_KINETIC_F2xy_F2xy_a-2.0E0*3*I_KINETIC_F2xy_F2xy_a+2*1*I_KINETIC_Py_F2xy;
  abcd[12] = 4.0E0*I_KINETIC_H4xz_F2xy_aa-2.0E0*2*I_KINETIC_F2xz_F2xy_a-2.0E0*3*I_KINETIC_F2xz_F2xy_a+2*1*I_KINETIC_Pz_F2xy;
  abcd[13] = 4.0E0*I_KINETIC_H3x2y_F2xy_aa-2.0E0*1*I_KINETIC_Fx2y_F2xy_a-2.0E0*2*I_KINETIC_Fx2y_F2xy_a;
  abcd[14] = 4.0E0*I_KINETIC_H3xyz_F2xy_aa-2.0E0*1*I_KINETIC_Fxyz_F2xy_a-2.0E0*2*I_KINETIC_Fxyz_F2xy_a;
  abcd[15] = 4.0E0*I_KINETIC_H3x2z_F2xy_aa-2.0E0*1*I_KINETIC_Fx2z_F2xy_a-2.0E0*2*I_KINETIC_Fx2z_F2xy_a;
  abcd[16] = 4.0E0*I_KINETIC_H2x3y_F2xy_aa-2.0E0*1*I_KINETIC_F3y_F2xy_a;
  abcd[17] = 4.0E0*I_KINETIC_H2x2yz_F2xy_aa-2.0E0*1*I_KINETIC_F2yz_F2xy_a;
  abcd[18] = 4.0E0*I_KINETIC_H2xy2z_F2xy_aa-2.0E0*1*I_KINETIC_Fy2z_F2xy_a;
  abcd[19] = 4.0E0*I_KINETIC_H2x3z_F2xy_aa-2.0E0*1*I_KINETIC_F3z_F2xy_a;
  abcd[20] = 4.0E0*I_KINETIC_H5x_F2xz_aa-2.0E0*3*I_KINETIC_F3x_F2xz_a-2.0E0*4*I_KINETIC_F3x_F2xz_a+3*2*I_KINETIC_Px_F2xz;
  abcd[21] = 4.0E0*I_KINETIC_H4xy_F2xz_aa-2.0E0*2*I_KINETIC_F2xy_F2xz_a-2.0E0*3*I_KINETIC_F2xy_F2xz_a+2*1*I_KINETIC_Py_F2xz;
  abcd[22] = 4.0E0*I_KINETIC_H4xz_F2xz_aa-2.0E0*2*I_KINETIC_F2xz_F2xz_a-2.0E0*3*I_KINETIC_F2xz_F2xz_a+2*1*I_KINETIC_Pz_F2xz;
  abcd[23] = 4.0E0*I_KINETIC_H3x2y_F2xz_aa-2.0E0*1*I_KINETIC_Fx2y_F2xz_a-2.0E0*2*I_KINETIC_Fx2y_F2xz_a;
  abcd[24] = 4.0E0*I_KINETIC_H3xyz_F2xz_aa-2.0E0*1*I_KINETIC_Fxyz_F2xz_a-2.0E0*2*I_KINETIC_Fxyz_F2xz_a;
  abcd[25] = 4.0E0*I_KINETIC_H3x2z_F2xz_aa-2.0E0*1*I_KINETIC_Fx2z_F2xz_a-2.0E0*2*I_KINETIC_Fx2z_F2xz_a;
  abcd[26] = 4.0E0*I_KINETIC_H2x3y_F2xz_aa-2.0E0*1*I_KINETIC_F3y_F2xz_a;
  abcd[27] = 4.0E0*I_KINETIC_H2x2yz_F2xz_aa-2.0E0*1*I_KINETIC_F2yz_F2xz_a;
  abcd[28] = 4.0E0*I_KINETIC_H2xy2z_F2xz_aa-2.0E0*1*I_KINETIC_Fy2z_F2xz_a;
  abcd[29] = 4.0E0*I_KINETIC_H2x3z_F2xz_aa-2.0E0*1*I_KINETIC_F3z_F2xz_a;
  abcd[30] = 4.0E0*I_KINETIC_H5x_Fx2y_aa-2.0E0*3*I_KINETIC_F3x_Fx2y_a-2.0E0*4*I_KINETIC_F3x_Fx2y_a+3*2*I_KINETIC_Px_Fx2y;
  abcd[31] = 4.0E0*I_KINETIC_H4xy_Fx2y_aa-2.0E0*2*I_KINETIC_F2xy_Fx2y_a-2.0E0*3*I_KINETIC_F2xy_Fx2y_a+2*1*I_KINETIC_Py_Fx2y;
  abcd[32] = 4.0E0*I_KINETIC_H4xz_Fx2y_aa-2.0E0*2*I_KINETIC_F2xz_Fx2y_a-2.0E0*3*I_KINETIC_F2xz_Fx2y_a+2*1*I_KINETIC_Pz_Fx2y;
  abcd[33] = 4.0E0*I_KINETIC_H3x2y_Fx2y_aa-2.0E0*1*I_KINETIC_Fx2y_Fx2y_a-2.0E0*2*I_KINETIC_Fx2y_Fx2y_a;
  abcd[34] = 4.0E0*I_KINETIC_H3xyz_Fx2y_aa-2.0E0*1*I_KINETIC_Fxyz_Fx2y_a-2.0E0*2*I_KINETIC_Fxyz_Fx2y_a;
  abcd[35] = 4.0E0*I_KINETIC_H3x2z_Fx2y_aa-2.0E0*1*I_KINETIC_Fx2z_Fx2y_a-2.0E0*2*I_KINETIC_Fx2z_Fx2y_a;
  abcd[36] = 4.0E0*I_KINETIC_H2x3y_Fx2y_aa-2.0E0*1*I_KINETIC_F3y_Fx2y_a;
  abcd[37] = 4.0E0*I_KINETIC_H2x2yz_Fx2y_aa-2.0E0*1*I_KINETIC_F2yz_Fx2y_a;
  abcd[38] = 4.0E0*I_KINETIC_H2xy2z_Fx2y_aa-2.0E0*1*I_KINETIC_Fy2z_Fx2y_a;
  abcd[39] = 4.0E0*I_KINETIC_H2x3z_Fx2y_aa-2.0E0*1*I_KINETIC_F3z_Fx2y_a;
  abcd[40] = 4.0E0*I_KINETIC_H5x_Fxyz_aa-2.0E0*3*I_KINETIC_F3x_Fxyz_a-2.0E0*4*I_KINETIC_F3x_Fxyz_a+3*2*I_KINETIC_Px_Fxyz;
  abcd[41] = 4.0E0*I_KINETIC_H4xy_Fxyz_aa-2.0E0*2*I_KINETIC_F2xy_Fxyz_a-2.0E0*3*I_KINETIC_F2xy_Fxyz_a+2*1*I_KINETIC_Py_Fxyz;
  abcd[42] = 4.0E0*I_KINETIC_H4xz_Fxyz_aa-2.0E0*2*I_KINETIC_F2xz_Fxyz_a-2.0E0*3*I_KINETIC_F2xz_Fxyz_a+2*1*I_KINETIC_Pz_Fxyz;
  abcd[43] = 4.0E0*I_KINETIC_H3x2y_Fxyz_aa-2.0E0*1*I_KINETIC_Fx2y_Fxyz_a-2.0E0*2*I_KINETIC_Fx2y_Fxyz_a;
  abcd[44] = 4.0E0*I_KINETIC_H3xyz_Fxyz_aa-2.0E0*1*I_KINETIC_Fxyz_Fxyz_a-2.0E0*2*I_KINETIC_Fxyz_Fxyz_a;
  abcd[45] = 4.0E0*I_KINETIC_H3x2z_Fxyz_aa-2.0E0*1*I_KINETIC_Fx2z_Fxyz_a-2.0E0*2*I_KINETIC_Fx2z_Fxyz_a;
  abcd[46] = 4.0E0*I_KINETIC_H2x3y_Fxyz_aa-2.0E0*1*I_KINETIC_F3y_Fxyz_a;
  abcd[47] = 4.0E0*I_KINETIC_H2x2yz_Fxyz_aa-2.0E0*1*I_KINETIC_F2yz_Fxyz_a;
  abcd[48] = 4.0E0*I_KINETIC_H2xy2z_Fxyz_aa-2.0E0*1*I_KINETIC_Fy2z_Fxyz_a;
  abcd[49] = 4.0E0*I_KINETIC_H2x3z_Fxyz_aa-2.0E0*1*I_KINETIC_F3z_Fxyz_a;
  abcd[50] = 4.0E0*I_KINETIC_H5x_Fx2z_aa-2.0E0*3*I_KINETIC_F3x_Fx2z_a-2.0E0*4*I_KINETIC_F3x_Fx2z_a+3*2*I_KINETIC_Px_Fx2z;
  abcd[51] = 4.0E0*I_KINETIC_H4xy_Fx2z_aa-2.0E0*2*I_KINETIC_F2xy_Fx2z_a-2.0E0*3*I_KINETIC_F2xy_Fx2z_a+2*1*I_KINETIC_Py_Fx2z;
  abcd[52] = 4.0E0*I_KINETIC_H4xz_Fx2z_aa-2.0E0*2*I_KINETIC_F2xz_Fx2z_a-2.0E0*3*I_KINETIC_F2xz_Fx2z_a+2*1*I_KINETIC_Pz_Fx2z;
  abcd[53] = 4.0E0*I_KINETIC_H3x2y_Fx2z_aa-2.0E0*1*I_KINETIC_Fx2y_Fx2z_a-2.0E0*2*I_KINETIC_Fx2y_Fx2z_a;
  abcd[54] = 4.0E0*I_KINETIC_H3xyz_Fx2z_aa-2.0E0*1*I_KINETIC_Fxyz_Fx2z_a-2.0E0*2*I_KINETIC_Fxyz_Fx2z_a;
  abcd[55] = 4.0E0*I_KINETIC_H3x2z_Fx2z_aa-2.0E0*1*I_KINETIC_Fx2z_Fx2z_a-2.0E0*2*I_KINETIC_Fx2z_Fx2z_a;
  abcd[56] = 4.0E0*I_KINETIC_H2x3y_Fx2z_aa-2.0E0*1*I_KINETIC_F3y_Fx2z_a;
  abcd[57] = 4.0E0*I_KINETIC_H2x2yz_Fx2z_aa-2.0E0*1*I_KINETIC_F2yz_Fx2z_a;
  abcd[58] = 4.0E0*I_KINETIC_H2xy2z_Fx2z_aa-2.0E0*1*I_KINETIC_Fy2z_Fx2z_a;
  abcd[59] = 4.0E0*I_KINETIC_H2x3z_Fx2z_aa-2.0E0*1*I_KINETIC_F3z_Fx2z_a;
  abcd[60] = 4.0E0*I_KINETIC_H5x_F3y_aa-2.0E0*3*I_KINETIC_F3x_F3y_a-2.0E0*4*I_KINETIC_F3x_F3y_a+3*2*I_KINETIC_Px_F3y;
  abcd[61] = 4.0E0*I_KINETIC_H4xy_F3y_aa-2.0E0*2*I_KINETIC_F2xy_F3y_a-2.0E0*3*I_KINETIC_F2xy_F3y_a+2*1*I_KINETIC_Py_F3y;
  abcd[62] = 4.0E0*I_KINETIC_H4xz_F3y_aa-2.0E0*2*I_KINETIC_F2xz_F3y_a-2.0E0*3*I_KINETIC_F2xz_F3y_a+2*1*I_KINETIC_Pz_F3y;
  abcd[63] = 4.0E0*I_KINETIC_H3x2y_F3y_aa-2.0E0*1*I_KINETIC_Fx2y_F3y_a-2.0E0*2*I_KINETIC_Fx2y_F3y_a;
  abcd[64] = 4.0E0*I_KINETIC_H3xyz_F3y_aa-2.0E0*1*I_KINETIC_Fxyz_F3y_a-2.0E0*2*I_KINETIC_Fxyz_F3y_a;
  abcd[65] = 4.0E0*I_KINETIC_H3x2z_F3y_aa-2.0E0*1*I_KINETIC_Fx2z_F3y_a-2.0E0*2*I_KINETIC_Fx2z_F3y_a;
  abcd[66] = 4.0E0*I_KINETIC_H2x3y_F3y_aa-2.0E0*1*I_KINETIC_F3y_F3y_a;
  abcd[67] = 4.0E0*I_KINETIC_H2x2yz_F3y_aa-2.0E0*1*I_KINETIC_F2yz_F3y_a;
  abcd[68] = 4.0E0*I_KINETIC_H2xy2z_F3y_aa-2.0E0*1*I_KINETIC_Fy2z_F3y_a;
  abcd[69] = 4.0E0*I_KINETIC_H2x3z_F3y_aa-2.0E0*1*I_KINETIC_F3z_F3y_a;
  abcd[70] = 4.0E0*I_KINETIC_H5x_F2yz_aa-2.0E0*3*I_KINETIC_F3x_F2yz_a-2.0E0*4*I_KINETIC_F3x_F2yz_a+3*2*I_KINETIC_Px_F2yz;
  abcd[71] = 4.0E0*I_KINETIC_H4xy_F2yz_aa-2.0E0*2*I_KINETIC_F2xy_F2yz_a-2.0E0*3*I_KINETIC_F2xy_F2yz_a+2*1*I_KINETIC_Py_F2yz;
  abcd[72] = 4.0E0*I_KINETIC_H4xz_F2yz_aa-2.0E0*2*I_KINETIC_F2xz_F2yz_a-2.0E0*3*I_KINETIC_F2xz_F2yz_a+2*1*I_KINETIC_Pz_F2yz;
  abcd[73] = 4.0E0*I_KINETIC_H3x2y_F2yz_aa-2.0E0*1*I_KINETIC_Fx2y_F2yz_a-2.0E0*2*I_KINETIC_Fx2y_F2yz_a;
  abcd[74] = 4.0E0*I_KINETIC_H3xyz_F2yz_aa-2.0E0*1*I_KINETIC_Fxyz_F2yz_a-2.0E0*2*I_KINETIC_Fxyz_F2yz_a;
  abcd[75] = 4.0E0*I_KINETIC_H3x2z_F2yz_aa-2.0E0*1*I_KINETIC_Fx2z_F2yz_a-2.0E0*2*I_KINETIC_Fx2z_F2yz_a;
  abcd[76] = 4.0E0*I_KINETIC_H2x3y_F2yz_aa-2.0E0*1*I_KINETIC_F3y_F2yz_a;
  abcd[77] = 4.0E0*I_KINETIC_H2x2yz_F2yz_aa-2.0E0*1*I_KINETIC_F2yz_F2yz_a;
  abcd[78] = 4.0E0*I_KINETIC_H2xy2z_F2yz_aa-2.0E0*1*I_KINETIC_Fy2z_F2yz_a;
  abcd[79] = 4.0E0*I_KINETIC_H2x3z_F2yz_aa-2.0E0*1*I_KINETIC_F3z_F2yz_a;
  abcd[80] = 4.0E0*I_KINETIC_H5x_Fy2z_aa-2.0E0*3*I_KINETIC_F3x_Fy2z_a-2.0E0*4*I_KINETIC_F3x_Fy2z_a+3*2*I_KINETIC_Px_Fy2z;
  abcd[81] = 4.0E0*I_KINETIC_H4xy_Fy2z_aa-2.0E0*2*I_KINETIC_F2xy_Fy2z_a-2.0E0*3*I_KINETIC_F2xy_Fy2z_a+2*1*I_KINETIC_Py_Fy2z;
  abcd[82] = 4.0E0*I_KINETIC_H4xz_Fy2z_aa-2.0E0*2*I_KINETIC_F2xz_Fy2z_a-2.0E0*3*I_KINETIC_F2xz_Fy2z_a+2*1*I_KINETIC_Pz_Fy2z;
  abcd[83] = 4.0E0*I_KINETIC_H3x2y_Fy2z_aa-2.0E0*1*I_KINETIC_Fx2y_Fy2z_a-2.0E0*2*I_KINETIC_Fx2y_Fy2z_a;
  abcd[84] = 4.0E0*I_KINETIC_H3xyz_Fy2z_aa-2.0E0*1*I_KINETIC_Fxyz_Fy2z_a-2.0E0*2*I_KINETIC_Fxyz_Fy2z_a;
  abcd[85] = 4.0E0*I_KINETIC_H3x2z_Fy2z_aa-2.0E0*1*I_KINETIC_Fx2z_Fy2z_a-2.0E0*2*I_KINETIC_Fx2z_Fy2z_a;
  abcd[86] = 4.0E0*I_KINETIC_H2x3y_Fy2z_aa-2.0E0*1*I_KINETIC_F3y_Fy2z_a;
  abcd[87] = 4.0E0*I_KINETIC_H2x2yz_Fy2z_aa-2.0E0*1*I_KINETIC_F2yz_Fy2z_a;
  abcd[88] = 4.0E0*I_KINETIC_H2xy2z_Fy2z_aa-2.0E0*1*I_KINETIC_Fy2z_Fy2z_a;
  abcd[89] = 4.0E0*I_KINETIC_H2x3z_Fy2z_aa-2.0E0*1*I_KINETIC_F3z_Fy2z_a;
  abcd[90] = 4.0E0*I_KINETIC_H5x_F3z_aa-2.0E0*3*I_KINETIC_F3x_F3z_a-2.0E0*4*I_KINETIC_F3x_F3z_a+3*2*I_KINETIC_Px_F3z;
  abcd[91] = 4.0E0*I_KINETIC_H4xy_F3z_aa-2.0E0*2*I_KINETIC_F2xy_F3z_a-2.0E0*3*I_KINETIC_F2xy_F3z_a+2*1*I_KINETIC_Py_F3z;
  abcd[92] = 4.0E0*I_KINETIC_H4xz_F3z_aa-2.0E0*2*I_KINETIC_F2xz_F3z_a-2.0E0*3*I_KINETIC_F2xz_F3z_a+2*1*I_KINETIC_Pz_F3z;
  abcd[93] = 4.0E0*I_KINETIC_H3x2y_F3z_aa-2.0E0*1*I_KINETIC_Fx2y_F3z_a-2.0E0*2*I_KINETIC_Fx2y_F3z_a;
  abcd[94] = 4.0E0*I_KINETIC_H3xyz_F3z_aa-2.0E0*1*I_KINETIC_Fxyz_F3z_a-2.0E0*2*I_KINETIC_Fxyz_F3z_a;
  abcd[95] = 4.0E0*I_KINETIC_H3x2z_F3z_aa-2.0E0*1*I_KINETIC_Fx2z_F3z_a-2.0E0*2*I_KINETIC_Fx2z_F3z_a;
  abcd[96] = 4.0E0*I_KINETIC_H2x3y_F3z_aa-2.0E0*1*I_KINETIC_F3y_F3z_a;
  abcd[97] = 4.0E0*I_KINETIC_H2x2yz_F3z_aa-2.0E0*1*I_KINETIC_F2yz_F3z_a;
  abcd[98] = 4.0E0*I_KINETIC_H2xy2z_F3z_aa-2.0E0*1*I_KINETIC_Fy2z_F3z_a;
  abcd[99] = 4.0E0*I_KINETIC_H2x3z_F3z_aa-2.0E0*1*I_KINETIC_F3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_F_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_F_aa
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_P_F
   ************************************************************/
  abcd[100] = 4.0E0*I_KINETIC_H4xy_F3x_aa-2.0E0*3*I_KINETIC_F2xy_F3x_a;
  abcd[101] = 4.0E0*I_KINETIC_H3x2y_F3x_aa-2.0E0*1*I_KINETIC_F3x_F3x_a-2.0E0*2*I_KINETIC_Fx2y_F3x_a+2*1*I_KINETIC_Px_F3x;
  abcd[102] = 4.0E0*I_KINETIC_H3xyz_F3x_aa-2.0E0*2*I_KINETIC_Fxyz_F3x_a;
  abcd[103] = 4.0E0*I_KINETIC_H2x3y_F3x_aa-2.0E0*2*I_KINETIC_F2xy_F3x_a-2.0E0*1*I_KINETIC_F3y_F3x_a+2*I_KINETIC_Py_F3x;
  abcd[104] = 4.0E0*I_KINETIC_H2x2yz_F3x_aa-2.0E0*1*I_KINETIC_F2xz_F3x_a-2.0E0*1*I_KINETIC_F2yz_F3x_a+1*I_KINETIC_Pz_F3x;
  abcd[105] = 4.0E0*I_KINETIC_H2xy2z_F3x_aa-2.0E0*1*I_KINETIC_Fy2z_F3x_a;
  abcd[106] = 4.0E0*I_KINETIC_Hx4y_F3x_aa-2.0E0*3*I_KINETIC_Fx2y_F3x_a;
  abcd[107] = 4.0E0*I_KINETIC_Hx3yz_F3x_aa-2.0E0*2*I_KINETIC_Fxyz_F3x_a;
  abcd[108] = 4.0E0*I_KINETIC_Hx2y2z_F3x_aa-2.0E0*1*I_KINETIC_Fx2z_F3x_a;
  abcd[109] = 4.0E0*I_KINETIC_Hxy3z_F3x_aa;
  abcd[110] = 4.0E0*I_KINETIC_H4xy_F2xy_aa-2.0E0*3*I_KINETIC_F2xy_F2xy_a;
  abcd[111] = 4.0E0*I_KINETIC_H3x2y_F2xy_aa-2.0E0*1*I_KINETIC_F3x_F2xy_a-2.0E0*2*I_KINETIC_Fx2y_F2xy_a+2*1*I_KINETIC_Px_F2xy;
  abcd[112] = 4.0E0*I_KINETIC_H3xyz_F2xy_aa-2.0E0*2*I_KINETIC_Fxyz_F2xy_a;
  abcd[113] = 4.0E0*I_KINETIC_H2x3y_F2xy_aa-2.0E0*2*I_KINETIC_F2xy_F2xy_a-2.0E0*1*I_KINETIC_F3y_F2xy_a+2*I_KINETIC_Py_F2xy;
  abcd[114] = 4.0E0*I_KINETIC_H2x2yz_F2xy_aa-2.0E0*1*I_KINETIC_F2xz_F2xy_a-2.0E0*1*I_KINETIC_F2yz_F2xy_a+1*I_KINETIC_Pz_F2xy;
  abcd[115] = 4.0E0*I_KINETIC_H2xy2z_F2xy_aa-2.0E0*1*I_KINETIC_Fy2z_F2xy_a;
  abcd[116] = 4.0E0*I_KINETIC_Hx4y_F2xy_aa-2.0E0*3*I_KINETIC_Fx2y_F2xy_a;
  abcd[117] = 4.0E0*I_KINETIC_Hx3yz_F2xy_aa-2.0E0*2*I_KINETIC_Fxyz_F2xy_a;
  abcd[118] = 4.0E0*I_KINETIC_Hx2y2z_F2xy_aa-2.0E0*1*I_KINETIC_Fx2z_F2xy_a;
  abcd[119] = 4.0E0*I_KINETIC_Hxy3z_F2xy_aa;
  abcd[120] = 4.0E0*I_KINETIC_H4xy_F2xz_aa-2.0E0*3*I_KINETIC_F2xy_F2xz_a;
  abcd[121] = 4.0E0*I_KINETIC_H3x2y_F2xz_aa-2.0E0*1*I_KINETIC_F3x_F2xz_a-2.0E0*2*I_KINETIC_Fx2y_F2xz_a+2*1*I_KINETIC_Px_F2xz;
  abcd[122] = 4.0E0*I_KINETIC_H3xyz_F2xz_aa-2.0E0*2*I_KINETIC_Fxyz_F2xz_a;
  abcd[123] = 4.0E0*I_KINETIC_H2x3y_F2xz_aa-2.0E0*2*I_KINETIC_F2xy_F2xz_a-2.0E0*1*I_KINETIC_F3y_F2xz_a+2*I_KINETIC_Py_F2xz;
  abcd[124] = 4.0E0*I_KINETIC_H2x2yz_F2xz_aa-2.0E0*1*I_KINETIC_F2xz_F2xz_a-2.0E0*1*I_KINETIC_F2yz_F2xz_a+1*I_KINETIC_Pz_F2xz;
  abcd[125] = 4.0E0*I_KINETIC_H2xy2z_F2xz_aa-2.0E0*1*I_KINETIC_Fy2z_F2xz_a;
  abcd[126] = 4.0E0*I_KINETIC_Hx4y_F2xz_aa-2.0E0*3*I_KINETIC_Fx2y_F2xz_a;
  abcd[127] = 4.0E0*I_KINETIC_Hx3yz_F2xz_aa-2.0E0*2*I_KINETIC_Fxyz_F2xz_a;
  abcd[128] = 4.0E0*I_KINETIC_Hx2y2z_F2xz_aa-2.0E0*1*I_KINETIC_Fx2z_F2xz_a;
  abcd[129] = 4.0E0*I_KINETIC_Hxy3z_F2xz_aa;
  abcd[130] = 4.0E0*I_KINETIC_H4xy_Fx2y_aa-2.0E0*3*I_KINETIC_F2xy_Fx2y_a;
  abcd[131] = 4.0E0*I_KINETIC_H3x2y_Fx2y_aa-2.0E0*1*I_KINETIC_F3x_Fx2y_a-2.0E0*2*I_KINETIC_Fx2y_Fx2y_a+2*1*I_KINETIC_Px_Fx2y;
  abcd[132] = 4.0E0*I_KINETIC_H3xyz_Fx2y_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2y_a;
  abcd[133] = 4.0E0*I_KINETIC_H2x3y_Fx2y_aa-2.0E0*2*I_KINETIC_F2xy_Fx2y_a-2.0E0*1*I_KINETIC_F3y_Fx2y_a+2*I_KINETIC_Py_Fx2y;
  abcd[134] = 4.0E0*I_KINETIC_H2x2yz_Fx2y_aa-2.0E0*1*I_KINETIC_F2xz_Fx2y_a-2.0E0*1*I_KINETIC_F2yz_Fx2y_a+1*I_KINETIC_Pz_Fx2y;
  abcd[135] = 4.0E0*I_KINETIC_H2xy2z_Fx2y_aa-2.0E0*1*I_KINETIC_Fy2z_Fx2y_a;
  abcd[136] = 4.0E0*I_KINETIC_Hx4y_Fx2y_aa-2.0E0*3*I_KINETIC_Fx2y_Fx2y_a;
  abcd[137] = 4.0E0*I_KINETIC_Hx3yz_Fx2y_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2y_a;
  abcd[138] = 4.0E0*I_KINETIC_Hx2y2z_Fx2y_aa-2.0E0*1*I_KINETIC_Fx2z_Fx2y_a;
  abcd[139] = 4.0E0*I_KINETIC_Hxy3z_Fx2y_aa;
  abcd[140] = 4.0E0*I_KINETIC_H4xy_Fxyz_aa-2.0E0*3*I_KINETIC_F2xy_Fxyz_a;
  abcd[141] = 4.0E0*I_KINETIC_H3x2y_Fxyz_aa-2.0E0*1*I_KINETIC_F3x_Fxyz_a-2.0E0*2*I_KINETIC_Fx2y_Fxyz_a+2*1*I_KINETIC_Px_Fxyz;
  abcd[142] = 4.0E0*I_KINETIC_H3xyz_Fxyz_aa-2.0E0*2*I_KINETIC_Fxyz_Fxyz_a;
  abcd[143] = 4.0E0*I_KINETIC_H2x3y_Fxyz_aa-2.0E0*2*I_KINETIC_F2xy_Fxyz_a-2.0E0*1*I_KINETIC_F3y_Fxyz_a+2*I_KINETIC_Py_Fxyz;
  abcd[144] = 4.0E0*I_KINETIC_H2x2yz_Fxyz_aa-2.0E0*1*I_KINETIC_F2xz_Fxyz_a-2.0E0*1*I_KINETIC_F2yz_Fxyz_a+1*I_KINETIC_Pz_Fxyz;
  abcd[145] = 4.0E0*I_KINETIC_H2xy2z_Fxyz_aa-2.0E0*1*I_KINETIC_Fy2z_Fxyz_a;
  abcd[146] = 4.0E0*I_KINETIC_Hx4y_Fxyz_aa-2.0E0*3*I_KINETIC_Fx2y_Fxyz_a;
  abcd[147] = 4.0E0*I_KINETIC_Hx3yz_Fxyz_aa-2.0E0*2*I_KINETIC_Fxyz_Fxyz_a;
  abcd[148] = 4.0E0*I_KINETIC_Hx2y2z_Fxyz_aa-2.0E0*1*I_KINETIC_Fx2z_Fxyz_a;
  abcd[149] = 4.0E0*I_KINETIC_Hxy3z_Fxyz_aa;
  abcd[150] = 4.0E0*I_KINETIC_H4xy_Fx2z_aa-2.0E0*3*I_KINETIC_F2xy_Fx2z_a;
  abcd[151] = 4.0E0*I_KINETIC_H3x2y_Fx2z_aa-2.0E0*1*I_KINETIC_F3x_Fx2z_a-2.0E0*2*I_KINETIC_Fx2y_Fx2z_a+2*1*I_KINETIC_Px_Fx2z;
  abcd[152] = 4.0E0*I_KINETIC_H3xyz_Fx2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2z_a;
  abcd[153] = 4.0E0*I_KINETIC_H2x3y_Fx2z_aa-2.0E0*2*I_KINETIC_F2xy_Fx2z_a-2.0E0*1*I_KINETIC_F3y_Fx2z_a+2*I_KINETIC_Py_Fx2z;
  abcd[154] = 4.0E0*I_KINETIC_H2x2yz_Fx2z_aa-2.0E0*1*I_KINETIC_F2xz_Fx2z_a-2.0E0*1*I_KINETIC_F2yz_Fx2z_a+1*I_KINETIC_Pz_Fx2z;
  abcd[155] = 4.0E0*I_KINETIC_H2xy2z_Fx2z_aa-2.0E0*1*I_KINETIC_Fy2z_Fx2z_a;
  abcd[156] = 4.0E0*I_KINETIC_Hx4y_Fx2z_aa-2.0E0*3*I_KINETIC_Fx2y_Fx2z_a;
  abcd[157] = 4.0E0*I_KINETIC_Hx3yz_Fx2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2z_a;
  abcd[158] = 4.0E0*I_KINETIC_Hx2y2z_Fx2z_aa-2.0E0*1*I_KINETIC_Fx2z_Fx2z_a;
  abcd[159] = 4.0E0*I_KINETIC_Hxy3z_Fx2z_aa;
  abcd[160] = 4.0E0*I_KINETIC_H4xy_F3y_aa-2.0E0*3*I_KINETIC_F2xy_F3y_a;
  abcd[161] = 4.0E0*I_KINETIC_H3x2y_F3y_aa-2.0E0*1*I_KINETIC_F3x_F3y_a-2.0E0*2*I_KINETIC_Fx2y_F3y_a+2*1*I_KINETIC_Px_F3y;
  abcd[162] = 4.0E0*I_KINETIC_H3xyz_F3y_aa-2.0E0*2*I_KINETIC_Fxyz_F3y_a;
  abcd[163] = 4.0E0*I_KINETIC_H2x3y_F3y_aa-2.0E0*2*I_KINETIC_F2xy_F3y_a-2.0E0*1*I_KINETIC_F3y_F3y_a+2*I_KINETIC_Py_F3y;
  abcd[164] = 4.0E0*I_KINETIC_H2x2yz_F3y_aa-2.0E0*1*I_KINETIC_F2xz_F3y_a-2.0E0*1*I_KINETIC_F2yz_F3y_a+1*I_KINETIC_Pz_F3y;
  abcd[165] = 4.0E0*I_KINETIC_H2xy2z_F3y_aa-2.0E0*1*I_KINETIC_Fy2z_F3y_a;
  abcd[166] = 4.0E0*I_KINETIC_Hx4y_F3y_aa-2.0E0*3*I_KINETIC_Fx2y_F3y_a;
  abcd[167] = 4.0E0*I_KINETIC_Hx3yz_F3y_aa-2.0E0*2*I_KINETIC_Fxyz_F3y_a;
  abcd[168] = 4.0E0*I_KINETIC_Hx2y2z_F3y_aa-2.0E0*1*I_KINETIC_Fx2z_F3y_a;
  abcd[169] = 4.0E0*I_KINETIC_Hxy3z_F3y_aa;
  abcd[170] = 4.0E0*I_KINETIC_H4xy_F2yz_aa-2.0E0*3*I_KINETIC_F2xy_F2yz_a;
  abcd[171] = 4.0E0*I_KINETIC_H3x2y_F2yz_aa-2.0E0*1*I_KINETIC_F3x_F2yz_a-2.0E0*2*I_KINETIC_Fx2y_F2yz_a+2*1*I_KINETIC_Px_F2yz;
  abcd[172] = 4.0E0*I_KINETIC_H3xyz_F2yz_aa-2.0E0*2*I_KINETIC_Fxyz_F2yz_a;
  abcd[173] = 4.0E0*I_KINETIC_H2x3y_F2yz_aa-2.0E0*2*I_KINETIC_F2xy_F2yz_a-2.0E0*1*I_KINETIC_F3y_F2yz_a+2*I_KINETIC_Py_F2yz;
  abcd[174] = 4.0E0*I_KINETIC_H2x2yz_F2yz_aa-2.0E0*1*I_KINETIC_F2xz_F2yz_a-2.0E0*1*I_KINETIC_F2yz_F2yz_a+1*I_KINETIC_Pz_F2yz;
  abcd[175] = 4.0E0*I_KINETIC_H2xy2z_F2yz_aa-2.0E0*1*I_KINETIC_Fy2z_F2yz_a;
  abcd[176] = 4.0E0*I_KINETIC_Hx4y_F2yz_aa-2.0E0*3*I_KINETIC_Fx2y_F2yz_a;
  abcd[177] = 4.0E0*I_KINETIC_Hx3yz_F2yz_aa-2.0E0*2*I_KINETIC_Fxyz_F2yz_a;
  abcd[178] = 4.0E0*I_KINETIC_Hx2y2z_F2yz_aa-2.0E0*1*I_KINETIC_Fx2z_F2yz_a;
  abcd[179] = 4.0E0*I_KINETIC_Hxy3z_F2yz_aa;
  abcd[180] = 4.0E0*I_KINETIC_H4xy_Fy2z_aa-2.0E0*3*I_KINETIC_F2xy_Fy2z_a;
  abcd[181] = 4.0E0*I_KINETIC_H3x2y_Fy2z_aa-2.0E0*1*I_KINETIC_F3x_Fy2z_a-2.0E0*2*I_KINETIC_Fx2y_Fy2z_a+2*1*I_KINETIC_Px_Fy2z;
  abcd[182] = 4.0E0*I_KINETIC_H3xyz_Fy2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fy2z_a;
  abcd[183] = 4.0E0*I_KINETIC_H2x3y_Fy2z_aa-2.0E0*2*I_KINETIC_F2xy_Fy2z_a-2.0E0*1*I_KINETIC_F3y_Fy2z_a+2*I_KINETIC_Py_Fy2z;
  abcd[184] = 4.0E0*I_KINETIC_H2x2yz_Fy2z_aa-2.0E0*1*I_KINETIC_F2xz_Fy2z_a-2.0E0*1*I_KINETIC_F2yz_Fy2z_a+1*I_KINETIC_Pz_Fy2z;
  abcd[185] = 4.0E0*I_KINETIC_H2xy2z_Fy2z_aa-2.0E0*1*I_KINETIC_Fy2z_Fy2z_a;
  abcd[186] = 4.0E0*I_KINETIC_Hx4y_Fy2z_aa-2.0E0*3*I_KINETIC_Fx2y_Fy2z_a;
  abcd[187] = 4.0E0*I_KINETIC_Hx3yz_Fy2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fy2z_a;
  abcd[188] = 4.0E0*I_KINETIC_Hx2y2z_Fy2z_aa-2.0E0*1*I_KINETIC_Fx2z_Fy2z_a;
  abcd[189] = 4.0E0*I_KINETIC_Hxy3z_Fy2z_aa;
  abcd[190] = 4.0E0*I_KINETIC_H4xy_F3z_aa-2.0E0*3*I_KINETIC_F2xy_F3z_a;
  abcd[191] = 4.0E0*I_KINETIC_H3x2y_F3z_aa-2.0E0*1*I_KINETIC_F3x_F3z_a-2.0E0*2*I_KINETIC_Fx2y_F3z_a+2*1*I_KINETIC_Px_F3z;
  abcd[192] = 4.0E0*I_KINETIC_H3xyz_F3z_aa-2.0E0*2*I_KINETIC_Fxyz_F3z_a;
  abcd[193] = 4.0E0*I_KINETIC_H2x3y_F3z_aa-2.0E0*2*I_KINETIC_F2xy_F3z_a-2.0E0*1*I_KINETIC_F3y_F3z_a+2*I_KINETIC_Py_F3z;
  abcd[194] = 4.0E0*I_KINETIC_H2x2yz_F3z_aa-2.0E0*1*I_KINETIC_F2xz_F3z_a-2.0E0*1*I_KINETIC_F2yz_F3z_a+1*I_KINETIC_Pz_F3z;
  abcd[195] = 4.0E0*I_KINETIC_H2xy2z_F3z_aa-2.0E0*1*I_KINETIC_Fy2z_F3z_a;
  abcd[196] = 4.0E0*I_KINETIC_Hx4y_F3z_aa-2.0E0*3*I_KINETIC_Fx2y_F3z_a;
  abcd[197] = 4.0E0*I_KINETIC_Hx3yz_F3z_aa-2.0E0*2*I_KINETIC_Fxyz_F3z_a;
  abcd[198] = 4.0E0*I_KINETIC_Hx2y2z_F3z_aa-2.0E0*1*I_KINETIC_Fx2z_F3z_a;
  abcd[199] = 4.0E0*I_KINETIC_Hxy3z_F3z_aa;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_F_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_F_aa
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_P_F
   ************************************************************/
  abcd[200] = 4.0E0*I_KINETIC_H4xz_F3x_aa-2.0E0*3*I_KINETIC_F2xz_F3x_a;
  abcd[201] = 4.0E0*I_KINETIC_H3xyz_F3x_aa-2.0E0*2*I_KINETIC_Fxyz_F3x_a;
  abcd[202] = 4.0E0*I_KINETIC_H3x2z_F3x_aa-2.0E0*1*I_KINETIC_F3x_F3x_a-2.0E0*2*I_KINETIC_Fx2z_F3x_a+2*1*I_KINETIC_Px_F3x;
  abcd[203] = 4.0E0*I_KINETIC_H2x2yz_F3x_aa-2.0E0*1*I_KINETIC_F2yz_F3x_a;
  abcd[204] = 4.0E0*I_KINETIC_H2xy2z_F3x_aa-2.0E0*1*I_KINETIC_F2xy_F3x_a-2.0E0*1*I_KINETIC_Fy2z_F3x_a+1*I_KINETIC_Py_F3x;
  abcd[205] = 4.0E0*I_KINETIC_H2x3z_F3x_aa-2.0E0*2*I_KINETIC_F2xz_F3x_a-2.0E0*1*I_KINETIC_F3z_F3x_a+2*I_KINETIC_Pz_F3x;
  abcd[206] = 4.0E0*I_KINETIC_Hx3yz_F3x_aa;
  abcd[207] = 4.0E0*I_KINETIC_Hx2y2z_F3x_aa-2.0E0*1*I_KINETIC_Fx2y_F3x_a;
  abcd[208] = 4.0E0*I_KINETIC_Hxy3z_F3x_aa-2.0E0*2*I_KINETIC_Fxyz_F3x_a;
  abcd[209] = 4.0E0*I_KINETIC_Hx4z_F3x_aa-2.0E0*3*I_KINETIC_Fx2z_F3x_a;
  abcd[210] = 4.0E0*I_KINETIC_H4xz_F2xy_aa-2.0E0*3*I_KINETIC_F2xz_F2xy_a;
  abcd[211] = 4.0E0*I_KINETIC_H3xyz_F2xy_aa-2.0E0*2*I_KINETIC_Fxyz_F2xy_a;
  abcd[212] = 4.0E0*I_KINETIC_H3x2z_F2xy_aa-2.0E0*1*I_KINETIC_F3x_F2xy_a-2.0E0*2*I_KINETIC_Fx2z_F2xy_a+2*1*I_KINETIC_Px_F2xy;
  abcd[213] = 4.0E0*I_KINETIC_H2x2yz_F2xy_aa-2.0E0*1*I_KINETIC_F2yz_F2xy_a;
  abcd[214] = 4.0E0*I_KINETIC_H2xy2z_F2xy_aa-2.0E0*1*I_KINETIC_F2xy_F2xy_a-2.0E0*1*I_KINETIC_Fy2z_F2xy_a+1*I_KINETIC_Py_F2xy;
  abcd[215] = 4.0E0*I_KINETIC_H2x3z_F2xy_aa-2.0E0*2*I_KINETIC_F2xz_F2xy_a-2.0E0*1*I_KINETIC_F3z_F2xy_a+2*I_KINETIC_Pz_F2xy;
  abcd[216] = 4.0E0*I_KINETIC_Hx3yz_F2xy_aa;
  abcd[217] = 4.0E0*I_KINETIC_Hx2y2z_F2xy_aa-2.0E0*1*I_KINETIC_Fx2y_F2xy_a;
  abcd[218] = 4.0E0*I_KINETIC_Hxy3z_F2xy_aa-2.0E0*2*I_KINETIC_Fxyz_F2xy_a;
  abcd[219] = 4.0E0*I_KINETIC_Hx4z_F2xy_aa-2.0E0*3*I_KINETIC_Fx2z_F2xy_a;
  abcd[220] = 4.0E0*I_KINETIC_H4xz_F2xz_aa-2.0E0*3*I_KINETIC_F2xz_F2xz_a;
  abcd[221] = 4.0E0*I_KINETIC_H3xyz_F2xz_aa-2.0E0*2*I_KINETIC_Fxyz_F2xz_a;
  abcd[222] = 4.0E0*I_KINETIC_H3x2z_F2xz_aa-2.0E0*1*I_KINETIC_F3x_F2xz_a-2.0E0*2*I_KINETIC_Fx2z_F2xz_a+2*1*I_KINETIC_Px_F2xz;
  abcd[223] = 4.0E0*I_KINETIC_H2x2yz_F2xz_aa-2.0E0*1*I_KINETIC_F2yz_F2xz_a;
  abcd[224] = 4.0E0*I_KINETIC_H2xy2z_F2xz_aa-2.0E0*1*I_KINETIC_F2xy_F2xz_a-2.0E0*1*I_KINETIC_Fy2z_F2xz_a+1*I_KINETIC_Py_F2xz;
  abcd[225] = 4.0E0*I_KINETIC_H2x3z_F2xz_aa-2.0E0*2*I_KINETIC_F2xz_F2xz_a-2.0E0*1*I_KINETIC_F3z_F2xz_a+2*I_KINETIC_Pz_F2xz;
  abcd[226] = 4.0E0*I_KINETIC_Hx3yz_F2xz_aa;
  abcd[227] = 4.0E0*I_KINETIC_Hx2y2z_F2xz_aa-2.0E0*1*I_KINETIC_Fx2y_F2xz_a;
  abcd[228] = 4.0E0*I_KINETIC_Hxy3z_F2xz_aa-2.0E0*2*I_KINETIC_Fxyz_F2xz_a;
  abcd[229] = 4.0E0*I_KINETIC_Hx4z_F2xz_aa-2.0E0*3*I_KINETIC_Fx2z_F2xz_a;
  abcd[230] = 4.0E0*I_KINETIC_H4xz_Fx2y_aa-2.0E0*3*I_KINETIC_F2xz_Fx2y_a;
  abcd[231] = 4.0E0*I_KINETIC_H3xyz_Fx2y_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2y_a;
  abcd[232] = 4.0E0*I_KINETIC_H3x2z_Fx2y_aa-2.0E0*1*I_KINETIC_F3x_Fx2y_a-2.0E0*2*I_KINETIC_Fx2z_Fx2y_a+2*1*I_KINETIC_Px_Fx2y;
  abcd[233] = 4.0E0*I_KINETIC_H2x2yz_Fx2y_aa-2.0E0*1*I_KINETIC_F2yz_Fx2y_a;
  abcd[234] = 4.0E0*I_KINETIC_H2xy2z_Fx2y_aa-2.0E0*1*I_KINETIC_F2xy_Fx2y_a-2.0E0*1*I_KINETIC_Fy2z_Fx2y_a+1*I_KINETIC_Py_Fx2y;
  abcd[235] = 4.0E0*I_KINETIC_H2x3z_Fx2y_aa-2.0E0*2*I_KINETIC_F2xz_Fx2y_a-2.0E0*1*I_KINETIC_F3z_Fx2y_a+2*I_KINETIC_Pz_Fx2y;
  abcd[236] = 4.0E0*I_KINETIC_Hx3yz_Fx2y_aa;
  abcd[237] = 4.0E0*I_KINETIC_Hx2y2z_Fx2y_aa-2.0E0*1*I_KINETIC_Fx2y_Fx2y_a;
  abcd[238] = 4.0E0*I_KINETIC_Hxy3z_Fx2y_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2y_a;
  abcd[239] = 4.0E0*I_KINETIC_Hx4z_Fx2y_aa-2.0E0*3*I_KINETIC_Fx2z_Fx2y_a;
  abcd[240] = 4.0E0*I_KINETIC_H4xz_Fxyz_aa-2.0E0*3*I_KINETIC_F2xz_Fxyz_a;
  abcd[241] = 4.0E0*I_KINETIC_H3xyz_Fxyz_aa-2.0E0*2*I_KINETIC_Fxyz_Fxyz_a;
  abcd[242] = 4.0E0*I_KINETIC_H3x2z_Fxyz_aa-2.0E0*1*I_KINETIC_F3x_Fxyz_a-2.0E0*2*I_KINETIC_Fx2z_Fxyz_a+2*1*I_KINETIC_Px_Fxyz;
  abcd[243] = 4.0E0*I_KINETIC_H2x2yz_Fxyz_aa-2.0E0*1*I_KINETIC_F2yz_Fxyz_a;
  abcd[244] = 4.0E0*I_KINETIC_H2xy2z_Fxyz_aa-2.0E0*1*I_KINETIC_F2xy_Fxyz_a-2.0E0*1*I_KINETIC_Fy2z_Fxyz_a+1*I_KINETIC_Py_Fxyz;
  abcd[245] = 4.0E0*I_KINETIC_H2x3z_Fxyz_aa-2.0E0*2*I_KINETIC_F2xz_Fxyz_a-2.0E0*1*I_KINETIC_F3z_Fxyz_a+2*I_KINETIC_Pz_Fxyz;
  abcd[246] = 4.0E0*I_KINETIC_Hx3yz_Fxyz_aa;
  abcd[247] = 4.0E0*I_KINETIC_Hx2y2z_Fxyz_aa-2.0E0*1*I_KINETIC_Fx2y_Fxyz_a;
  abcd[248] = 4.0E0*I_KINETIC_Hxy3z_Fxyz_aa-2.0E0*2*I_KINETIC_Fxyz_Fxyz_a;
  abcd[249] = 4.0E0*I_KINETIC_Hx4z_Fxyz_aa-2.0E0*3*I_KINETIC_Fx2z_Fxyz_a;
  abcd[250] = 4.0E0*I_KINETIC_H4xz_Fx2z_aa-2.0E0*3*I_KINETIC_F2xz_Fx2z_a;
  abcd[251] = 4.0E0*I_KINETIC_H3xyz_Fx2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2z_a;
  abcd[252] = 4.0E0*I_KINETIC_H3x2z_Fx2z_aa-2.0E0*1*I_KINETIC_F3x_Fx2z_a-2.0E0*2*I_KINETIC_Fx2z_Fx2z_a+2*1*I_KINETIC_Px_Fx2z;
  abcd[253] = 4.0E0*I_KINETIC_H2x2yz_Fx2z_aa-2.0E0*1*I_KINETIC_F2yz_Fx2z_a;
  abcd[254] = 4.0E0*I_KINETIC_H2xy2z_Fx2z_aa-2.0E0*1*I_KINETIC_F2xy_Fx2z_a-2.0E0*1*I_KINETIC_Fy2z_Fx2z_a+1*I_KINETIC_Py_Fx2z;
  abcd[255] = 4.0E0*I_KINETIC_H2x3z_Fx2z_aa-2.0E0*2*I_KINETIC_F2xz_Fx2z_a-2.0E0*1*I_KINETIC_F3z_Fx2z_a+2*I_KINETIC_Pz_Fx2z;
  abcd[256] = 4.0E0*I_KINETIC_Hx3yz_Fx2z_aa;
  abcd[257] = 4.0E0*I_KINETIC_Hx2y2z_Fx2z_aa-2.0E0*1*I_KINETIC_Fx2y_Fx2z_a;
  abcd[258] = 4.0E0*I_KINETIC_Hxy3z_Fx2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2z_a;
  abcd[259] = 4.0E0*I_KINETIC_Hx4z_Fx2z_aa-2.0E0*3*I_KINETIC_Fx2z_Fx2z_a;
  abcd[260] = 4.0E0*I_KINETIC_H4xz_F3y_aa-2.0E0*3*I_KINETIC_F2xz_F3y_a;
  abcd[261] = 4.0E0*I_KINETIC_H3xyz_F3y_aa-2.0E0*2*I_KINETIC_Fxyz_F3y_a;
  abcd[262] = 4.0E0*I_KINETIC_H3x2z_F3y_aa-2.0E0*1*I_KINETIC_F3x_F3y_a-2.0E0*2*I_KINETIC_Fx2z_F3y_a+2*1*I_KINETIC_Px_F3y;
  abcd[263] = 4.0E0*I_KINETIC_H2x2yz_F3y_aa-2.0E0*1*I_KINETIC_F2yz_F3y_a;
  abcd[264] = 4.0E0*I_KINETIC_H2xy2z_F3y_aa-2.0E0*1*I_KINETIC_F2xy_F3y_a-2.0E0*1*I_KINETIC_Fy2z_F3y_a+1*I_KINETIC_Py_F3y;
  abcd[265] = 4.0E0*I_KINETIC_H2x3z_F3y_aa-2.0E0*2*I_KINETIC_F2xz_F3y_a-2.0E0*1*I_KINETIC_F3z_F3y_a+2*I_KINETIC_Pz_F3y;
  abcd[266] = 4.0E0*I_KINETIC_Hx3yz_F3y_aa;
  abcd[267] = 4.0E0*I_KINETIC_Hx2y2z_F3y_aa-2.0E0*1*I_KINETIC_Fx2y_F3y_a;
  abcd[268] = 4.0E0*I_KINETIC_Hxy3z_F3y_aa-2.0E0*2*I_KINETIC_Fxyz_F3y_a;
  abcd[269] = 4.0E0*I_KINETIC_Hx4z_F3y_aa-2.0E0*3*I_KINETIC_Fx2z_F3y_a;
  abcd[270] = 4.0E0*I_KINETIC_H4xz_F2yz_aa-2.0E0*3*I_KINETIC_F2xz_F2yz_a;
  abcd[271] = 4.0E0*I_KINETIC_H3xyz_F2yz_aa-2.0E0*2*I_KINETIC_Fxyz_F2yz_a;
  abcd[272] = 4.0E0*I_KINETIC_H3x2z_F2yz_aa-2.0E0*1*I_KINETIC_F3x_F2yz_a-2.0E0*2*I_KINETIC_Fx2z_F2yz_a+2*1*I_KINETIC_Px_F2yz;
  abcd[273] = 4.0E0*I_KINETIC_H2x2yz_F2yz_aa-2.0E0*1*I_KINETIC_F2yz_F2yz_a;
  abcd[274] = 4.0E0*I_KINETIC_H2xy2z_F2yz_aa-2.0E0*1*I_KINETIC_F2xy_F2yz_a-2.0E0*1*I_KINETIC_Fy2z_F2yz_a+1*I_KINETIC_Py_F2yz;
  abcd[275] = 4.0E0*I_KINETIC_H2x3z_F2yz_aa-2.0E0*2*I_KINETIC_F2xz_F2yz_a-2.0E0*1*I_KINETIC_F3z_F2yz_a+2*I_KINETIC_Pz_F2yz;
  abcd[276] = 4.0E0*I_KINETIC_Hx3yz_F2yz_aa;
  abcd[277] = 4.0E0*I_KINETIC_Hx2y2z_F2yz_aa-2.0E0*1*I_KINETIC_Fx2y_F2yz_a;
  abcd[278] = 4.0E0*I_KINETIC_Hxy3z_F2yz_aa-2.0E0*2*I_KINETIC_Fxyz_F2yz_a;
  abcd[279] = 4.0E0*I_KINETIC_Hx4z_F2yz_aa-2.0E0*3*I_KINETIC_Fx2z_F2yz_a;
  abcd[280] = 4.0E0*I_KINETIC_H4xz_Fy2z_aa-2.0E0*3*I_KINETIC_F2xz_Fy2z_a;
  abcd[281] = 4.0E0*I_KINETIC_H3xyz_Fy2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fy2z_a;
  abcd[282] = 4.0E0*I_KINETIC_H3x2z_Fy2z_aa-2.0E0*1*I_KINETIC_F3x_Fy2z_a-2.0E0*2*I_KINETIC_Fx2z_Fy2z_a+2*1*I_KINETIC_Px_Fy2z;
  abcd[283] = 4.0E0*I_KINETIC_H2x2yz_Fy2z_aa-2.0E0*1*I_KINETIC_F2yz_Fy2z_a;
  abcd[284] = 4.0E0*I_KINETIC_H2xy2z_Fy2z_aa-2.0E0*1*I_KINETIC_F2xy_Fy2z_a-2.0E0*1*I_KINETIC_Fy2z_Fy2z_a+1*I_KINETIC_Py_Fy2z;
  abcd[285] = 4.0E0*I_KINETIC_H2x3z_Fy2z_aa-2.0E0*2*I_KINETIC_F2xz_Fy2z_a-2.0E0*1*I_KINETIC_F3z_Fy2z_a+2*I_KINETIC_Pz_Fy2z;
  abcd[286] = 4.0E0*I_KINETIC_Hx3yz_Fy2z_aa;
  abcd[287] = 4.0E0*I_KINETIC_Hx2y2z_Fy2z_aa-2.0E0*1*I_KINETIC_Fx2y_Fy2z_a;
  abcd[288] = 4.0E0*I_KINETIC_Hxy3z_Fy2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fy2z_a;
  abcd[289] = 4.0E0*I_KINETIC_Hx4z_Fy2z_aa-2.0E0*3*I_KINETIC_Fx2z_Fy2z_a;
  abcd[290] = 4.0E0*I_KINETIC_H4xz_F3z_aa-2.0E0*3*I_KINETIC_F2xz_F3z_a;
  abcd[291] = 4.0E0*I_KINETIC_H3xyz_F3z_aa-2.0E0*2*I_KINETIC_Fxyz_F3z_a;
  abcd[292] = 4.0E0*I_KINETIC_H3x2z_F3z_aa-2.0E0*1*I_KINETIC_F3x_F3z_a-2.0E0*2*I_KINETIC_Fx2z_F3z_a+2*1*I_KINETIC_Px_F3z;
  abcd[293] = 4.0E0*I_KINETIC_H2x2yz_F3z_aa-2.0E0*1*I_KINETIC_F2yz_F3z_a;
  abcd[294] = 4.0E0*I_KINETIC_H2xy2z_F3z_aa-2.0E0*1*I_KINETIC_F2xy_F3z_a-2.0E0*1*I_KINETIC_Fy2z_F3z_a+1*I_KINETIC_Py_F3z;
  abcd[295] = 4.0E0*I_KINETIC_H2x3z_F3z_aa-2.0E0*2*I_KINETIC_F2xz_F3z_a-2.0E0*1*I_KINETIC_F3z_F3z_a+2*I_KINETIC_Pz_F3z;
  abcd[296] = 4.0E0*I_KINETIC_Hx3yz_F3z_aa;
  abcd[297] = 4.0E0*I_KINETIC_Hx2y2z_F3z_aa-2.0E0*1*I_KINETIC_Fx2y_F3z_a;
  abcd[298] = 4.0E0*I_KINETIC_Hxy3z_F3z_aa-2.0E0*2*I_KINETIC_Fxyz_F3z_a;
  abcd[299] = 4.0E0*I_KINETIC_Hx4z_F3z_aa-2.0E0*3*I_KINETIC_Fx2z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_F_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_F_aa
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_P_F
   ************************************************************/
  abcd[300] = 4.0E0*I_KINETIC_H3x2y_F3x_aa-2.0E0*1*I_KINETIC_F3x_F3x_a;
  abcd[301] = 4.0E0*I_KINETIC_H2x3y_F3x_aa-2.0E0*1*I_KINETIC_F2xy_F3x_a-2.0E0*2*I_KINETIC_F2xy_F3x_a;
  abcd[302] = 4.0E0*I_KINETIC_H2x2yz_F3x_aa-2.0E0*1*I_KINETIC_F2xz_F3x_a;
  abcd[303] = 4.0E0*I_KINETIC_Hx4y_F3x_aa-2.0E0*2*I_KINETIC_Fx2y_F3x_a-2.0E0*3*I_KINETIC_Fx2y_F3x_a+2*1*I_KINETIC_Px_F3x;
  abcd[304] = 4.0E0*I_KINETIC_Hx3yz_F3x_aa-2.0E0*1*I_KINETIC_Fxyz_F3x_a-2.0E0*2*I_KINETIC_Fxyz_F3x_a;
  abcd[305] = 4.0E0*I_KINETIC_Hx2y2z_F3x_aa-2.0E0*1*I_KINETIC_Fx2z_F3x_a;
  abcd[306] = 4.0E0*I_KINETIC_H5y_F3x_aa-2.0E0*3*I_KINETIC_F3y_F3x_a-2.0E0*4*I_KINETIC_F3y_F3x_a+3*2*I_KINETIC_Py_F3x;
  abcd[307] = 4.0E0*I_KINETIC_H4yz_F3x_aa-2.0E0*2*I_KINETIC_F2yz_F3x_a-2.0E0*3*I_KINETIC_F2yz_F3x_a+2*1*I_KINETIC_Pz_F3x;
  abcd[308] = 4.0E0*I_KINETIC_H3y2z_F3x_aa-2.0E0*1*I_KINETIC_Fy2z_F3x_a-2.0E0*2*I_KINETIC_Fy2z_F3x_a;
  abcd[309] = 4.0E0*I_KINETIC_H2y3z_F3x_aa-2.0E0*1*I_KINETIC_F3z_F3x_a;
  abcd[310] = 4.0E0*I_KINETIC_H3x2y_F2xy_aa-2.0E0*1*I_KINETIC_F3x_F2xy_a;
  abcd[311] = 4.0E0*I_KINETIC_H2x3y_F2xy_aa-2.0E0*1*I_KINETIC_F2xy_F2xy_a-2.0E0*2*I_KINETIC_F2xy_F2xy_a;
  abcd[312] = 4.0E0*I_KINETIC_H2x2yz_F2xy_aa-2.0E0*1*I_KINETIC_F2xz_F2xy_a;
  abcd[313] = 4.0E0*I_KINETIC_Hx4y_F2xy_aa-2.0E0*2*I_KINETIC_Fx2y_F2xy_a-2.0E0*3*I_KINETIC_Fx2y_F2xy_a+2*1*I_KINETIC_Px_F2xy;
  abcd[314] = 4.0E0*I_KINETIC_Hx3yz_F2xy_aa-2.0E0*1*I_KINETIC_Fxyz_F2xy_a-2.0E0*2*I_KINETIC_Fxyz_F2xy_a;
  abcd[315] = 4.0E0*I_KINETIC_Hx2y2z_F2xy_aa-2.0E0*1*I_KINETIC_Fx2z_F2xy_a;
  abcd[316] = 4.0E0*I_KINETIC_H5y_F2xy_aa-2.0E0*3*I_KINETIC_F3y_F2xy_a-2.0E0*4*I_KINETIC_F3y_F2xy_a+3*2*I_KINETIC_Py_F2xy;
  abcd[317] = 4.0E0*I_KINETIC_H4yz_F2xy_aa-2.0E0*2*I_KINETIC_F2yz_F2xy_a-2.0E0*3*I_KINETIC_F2yz_F2xy_a+2*1*I_KINETIC_Pz_F2xy;
  abcd[318] = 4.0E0*I_KINETIC_H3y2z_F2xy_aa-2.0E0*1*I_KINETIC_Fy2z_F2xy_a-2.0E0*2*I_KINETIC_Fy2z_F2xy_a;
  abcd[319] = 4.0E0*I_KINETIC_H2y3z_F2xy_aa-2.0E0*1*I_KINETIC_F3z_F2xy_a;
  abcd[320] = 4.0E0*I_KINETIC_H3x2y_F2xz_aa-2.0E0*1*I_KINETIC_F3x_F2xz_a;
  abcd[321] = 4.0E0*I_KINETIC_H2x3y_F2xz_aa-2.0E0*1*I_KINETIC_F2xy_F2xz_a-2.0E0*2*I_KINETIC_F2xy_F2xz_a;
  abcd[322] = 4.0E0*I_KINETIC_H2x2yz_F2xz_aa-2.0E0*1*I_KINETIC_F2xz_F2xz_a;
  abcd[323] = 4.0E0*I_KINETIC_Hx4y_F2xz_aa-2.0E0*2*I_KINETIC_Fx2y_F2xz_a-2.0E0*3*I_KINETIC_Fx2y_F2xz_a+2*1*I_KINETIC_Px_F2xz;
  abcd[324] = 4.0E0*I_KINETIC_Hx3yz_F2xz_aa-2.0E0*1*I_KINETIC_Fxyz_F2xz_a-2.0E0*2*I_KINETIC_Fxyz_F2xz_a;
  abcd[325] = 4.0E0*I_KINETIC_Hx2y2z_F2xz_aa-2.0E0*1*I_KINETIC_Fx2z_F2xz_a;
  abcd[326] = 4.0E0*I_KINETIC_H5y_F2xz_aa-2.0E0*3*I_KINETIC_F3y_F2xz_a-2.0E0*4*I_KINETIC_F3y_F2xz_a+3*2*I_KINETIC_Py_F2xz;
  abcd[327] = 4.0E0*I_KINETIC_H4yz_F2xz_aa-2.0E0*2*I_KINETIC_F2yz_F2xz_a-2.0E0*3*I_KINETIC_F2yz_F2xz_a+2*1*I_KINETIC_Pz_F2xz;
  abcd[328] = 4.0E0*I_KINETIC_H3y2z_F2xz_aa-2.0E0*1*I_KINETIC_Fy2z_F2xz_a-2.0E0*2*I_KINETIC_Fy2z_F2xz_a;
  abcd[329] = 4.0E0*I_KINETIC_H2y3z_F2xz_aa-2.0E0*1*I_KINETIC_F3z_F2xz_a;
  abcd[330] = 4.0E0*I_KINETIC_H3x2y_Fx2y_aa-2.0E0*1*I_KINETIC_F3x_Fx2y_a;
  abcd[331] = 4.0E0*I_KINETIC_H2x3y_Fx2y_aa-2.0E0*1*I_KINETIC_F2xy_Fx2y_a-2.0E0*2*I_KINETIC_F2xy_Fx2y_a;
  abcd[332] = 4.0E0*I_KINETIC_H2x2yz_Fx2y_aa-2.0E0*1*I_KINETIC_F2xz_Fx2y_a;
  abcd[333] = 4.0E0*I_KINETIC_Hx4y_Fx2y_aa-2.0E0*2*I_KINETIC_Fx2y_Fx2y_a-2.0E0*3*I_KINETIC_Fx2y_Fx2y_a+2*1*I_KINETIC_Px_Fx2y;
  abcd[334] = 4.0E0*I_KINETIC_Hx3yz_Fx2y_aa-2.0E0*1*I_KINETIC_Fxyz_Fx2y_a-2.0E0*2*I_KINETIC_Fxyz_Fx2y_a;
  abcd[335] = 4.0E0*I_KINETIC_Hx2y2z_Fx2y_aa-2.0E0*1*I_KINETIC_Fx2z_Fx2y_a;
  abcd[336] = 4.0E0*I_KINETIC_H5y_Fx2y_aa-2.0E0*3*I_KINETIC_F3y_Fx2y_a-2.0E0*4*I_KINETIC_F3y_Fx2y_a+3*2*I_KINETIC_Py_Fx2y;
  abcd[337] = 4.0E0*I_KINETIC_H4yz_Fx2y_aa-2.0E0*2*I_KINETIC_F2yz_Fx2y_a-2.0E0*3*I_KINETIC_F2yz_Fx2y_a+2*1*I_KINETIC_Pz_Fx2y;
  abcd[338] = 4.0E0*I_KINETIC_H3y2z_Fx2y_aa-2.0E0*1*I_KINETIC_Fy2z_Fx2y_a-2.0E0*2*I_KINETIC_Fy2z_Fx2y_a;
  abcd[339] = 4.0E0*I_KINETIC_H2y3z_Fx2y_aa-2.0E0*1*I_KINETIC_F3z_Fx2y_a;
  abcd[340] = 4.0E0*I_KINETIC_H3x2y_Fxyz_aa-2.0E0*1*I_KINETIC_F3x_Fxyz_a;
  abcd[341] = 4.0E0*I_KINETIC_H2x3y_Fxyz_aa-2.0E0*1*I_KINETIC_F2xy_Fxyz_a-2.0E0*2*I_KINETIC_F2xy_Fxyz_a;
  abcd[342] = 4.0E0*I_KINETIC_H2x2yz_Fxyz_aa-2.0E0*1*I_KINETIC_F2xz_Fxyz_a;
  abcd[343] = 4.0E0*I_KINETIC_Hx4y_Fxyz_aa-2.0E0*2*I_KINETIC_Fx2y_Fxyz_a-2.0E0*3*I_KINETIC_Fx2y_Fxyz_a+2*1*I_KINETIC_Px_Fxyz;
  abcd[344] = 4.0E0*I_KINETIC_Hx3yz_Fxyz_aa-2.0E0*1*I_KINETIC_Fxyz_Fxyz_a-2.0E0*2*I_KINETIC_Fxyz_Fxyz_a;
  abcd[345] = 4.0E0*I_KINETIC_Hx2y2z_Fxyz_aa-2.0E0*1*I_KINETIC_Fx2z_Fxyz_a;
  abcd[346] = 4.0E0*I_KINETIC_H5y_Fxyz_aa-2.0E0*3*I_KINETIC_F3y_Fxyz_a-2.0E0*4*I_KINETIC_F3y_Fxyz_a+3*2*I_KINETIC_Py_Fxyz;
  abcd[347] = 4.0E0*I_KINETIC_H4yz_Fxyz_aa-2.0E0*2*I_KINETIC_F2yz_Fxyz_a-2.0E0*3*I_KINETIC_F2yz_Fxyz_a+2*1*I_KINETIC_Pz_Fxyz;
  abcd[348] = 4.0E0*I_KINETIC_H3y2z_Fxyz_aa-2.0E0*1*I_KINETIC_Fy2z_Fxyz_a-2.0E0*2*I_KINETIC_Fy2z_Fxyz_a;
  abcd[349] = 4.0E0*I_KINETIC_H2y3z_Fxyz_aa-2.0E0*1*I_KINETIC_F3z_Fxyz_a;
  abcd[350] = 4.0E0*I_KINETIC_H3x2y_Fx2z_aa-2.0E0*1*I_KINETIC_F3x_Fx2z_a;
  abcd[351] = 4.0E0*I_KINETIC_H2x3y_Fx2z_aa-2.0E0*1*I_KINETIC_F2xy_Fx2z_a-2.0E0*2*I_KINETIC_F2xy_Fx2z_a;
  abcd[352] = 4.0E0*I_KINETIC_H2x2yz_Fx2z_aa-2.0E0*1*I_KINETIC_F2xz_Fx2z_a;
  abcd[353] = 4.0E0*I_KINETIC_Hx4y_Fx2z_aa-2.0E0*2*I_KINETIC_Fx2y_Fx2z_a-2.0E0*3*I_KINETIC_Fx2y_Fx2z_a+2*1*I_KINETIC_Px_Fx2z;
  abcd[354] = 4.0E0*I_KINETIC_Hx3yz_Fx2z_aa-2.0E0*1*I_KINETIC_Fxyz_Fx2z_a-2.0E0*2*I_KINETIC_Fxyz_Fx2z_a;
  abcd[355] = 4.0E0*I_KINETIC_Hx2y2z_Fx2z_aa-2.0E0*1*I_KINETIC_Fx2z_Fx2z_a;
  abcd[356] = 4.0E0*I_KINETIC_H5y_Fx2z_aa-2.0E0*3*I_KINETIC_F3y_Fx2z_a-2.0E0*4*I_KINETIC_F3y_Fx2z_a+3*2*I_KINETIC_Py_Fx2z;
  abcd[357] = 4.0E0*I_KINETIC_H4yz_Fx2z_aa-2.0E0*2*I_KINETIC_F2yz_Fx2z_a-2.0E0*3*I_KINETIC_F2yz_Fx2z_a+2*1*I_KINETIC_Pz_Fx2z;
  abcd[358] = 4.0E0*I_KINETIC_H3y2z_Fx2z_aa-2.0E0*1*I_KINETIC_Fy2z_Fx2z_a-2.0E0*2*I_KINETIC_Fy2z_Fx2z_a;
  abcd[359] = 4.0E0*I_KINETIC_H2y3z_Fx2z_aa-2.0E0*1*I_KINETIC_F3z_Fx2z_a;
  abcd[360] = 4.0E0*I_KINETIC_H3x2y_F3y_aa-2.0E0*1*I_KINETIC_F3x_F3y_a;
  abcd[361] = 4.0E0*I_KINETIC_H2x3y_F3y_aa-2.0E0*1*I_KINETIC_F2xy_F3y_a-2.0E0*2*I_KINETIC_F2xy_F3y_a;
  abcd[362] = 4.0E0*I_KINETIC_H2x2yz_F3y_aa-2.0E0*1*I_KINETIC_F2xz_F3y_a;
  abcd[363] = 4.0E0*I_KINETIC_Hx4y_F3y_aa-2.0E0*2*I_KINETIC_Fx2y_F3y_a-2.0E0*3*I_KINETIC_Fx2y_F3y_a+2*1*I_KINETIC_Px_F3y;
  abcd[364] = 4.0E0*I_KINETIC_Hx3yz_F3y_aa-2.0E0*1*I_KINETIC_Fxyz_F3y_a-2.0E0*2*I_KINETIC_Fxyz_F3y_a;
  abcd[365] = 4.0E0*I_KINETIC_Hx2y2z_F3y_aa-2.0E0*1*I_KINETIC_Fx2z_F3y_a;
  abcd[366] = 4.0E0*I_KINETIC_H5y_F3y_aa-2.0E0*3*I_KINETIC_F3y_F3y_a-2.0E0*4*I_KINETIC_F3y_F3y_a+3*2*I_KINETIC_Py_F3y;
  abcd[367] = 4.0E0*I_KINETIC_H4yz_F3y_aa-2.0E0*2*I_KINETIC_F2yz_F3y_a-2.0E0*3*I_KINETIC_F2yz_F3y_a+2*1*I_KINETIC_Pz_F3y;
  abcd[368] = 4.0E0*I_KINETIC_H3y2z_F3y_aa-2.0E0*1*I_KINETIC_Fy2z_F3y_a-2.0E0*2*I_KINETIC_Fy2z_F3y_a;
  abcd[369] = 4.0E0*I_KINETIC_H2y3z_F3y_aa-2.0E0*1*I_KINETIC_F3z_F3y_a;
  abcd[370] = 4.0E0*I_KINETIC_H3x2y_F2yz_aa-2.0E0*1*I_KINETIC_F3x_F2yz_a;
  abcd[371] = 4.0E0*I_KINETIC_H2x3y_F2yz_aa-2.0E0*1*I_KINETIC_F2xy_F2yz_a-2.0E0*2*I_KINETIC_F2xy_F2yz_a;
  abcd[372] = 4.0E0*I_KINETIC_H2x2yz_F2yz_aa-2.0E0*1*I_KINETIC_F2xz_F2yz_a;
  abcd[373] = 4.0E0*I_KINETIC_Hx4y_F2yz_aa-2.0E0*2*I_KINETIC_Fx2y_F2yz_a-2.0E0*3*I_KINETIC_Fx2y_F2yz_a+2*1*I_KINETIC_Px_F2yz;
  abcd[374] = 4.0E0*I_KINETIC_Hx3yz_F2yz_aa-2.0E0*1*I_KINETIC_Fxyz_F2yz_a-2.0E0*2*I_KINETIC_Fxyz_F2yz_a;
  abcd[375] = 4.0E0*I_KINETIC_Hx2y2z_F2yz_aa-2.0E0*1*I_KINETIC_Fx2z_F2yz_a;
  abcd[376] = 4.0E0*I_KINETIC_H5y_F2yz_aa-2.0E0*3*I_KINETIC_F3y_F2yz_a-2.0E0*4*I_KINETIC_F3y_F2yz_a+3*2*I_KINETIC_Py_F2yz;
  abcd[377] = 4.0E0*I_KINETIC_H4yz_F2yz_aa-2.0E0*2*I_KINETIC_F2yz_F2yz_a-2.0E0*3*I_KINETIC_F2yz_F2yz_a+2*1*I_KINETIC_Pz_F2yz;
  abcd[378] = 4.0E0*I_KINETIC_H3y2z_F2yz_aa-2.0E0*1*I_KINETIC_Fy2z_F2yz_a-2.0E0*2*I_KINETIC_Fy2z_F2yz_a;
  abcd[379] = 4.0E0*I_KINETIC_H2y3z_F2yz_aa-2.0E0*1*I_KINETIC_F3z_F2yz_a;
  abcd[380] = 4.0E0*I_KINETIC_H3x2y_Fy2z_aa-2.0E0*1*I_KINETIC_F3x_Fy2z_a;
  abcd[381] = 4.0E0*I_KINETIC_H2x3y_Fy2z_aa-2.0E0*1*I_KINETIC_F2xy_Fy2z_a-2.0E0*2*I_KINETIC_F2xy_Fy2z_a;
  abcd[382] = 4.0E0*I_KINETIC_H2x2yz_Fy2z_aa-2.0E0*1*I_KINETIC_F2xz_Fy2z_a;
  abcd[383] = 4.0E0*I_KINETIC_Hx4y_Fy2z_aa-2.0E0*2*I_KINETIC_Fx2y_Fy2z_a-2.0E0*3*I_KINETIC_Fx2y_Fy2z_a+2*1*I_KINETIC_Px_Fy2z;
  abcd[384] = 4.0E0*I_KINETIC_Hx3yz_Fy2z_aa-2.0E0*1*I_KINETIC_Fxyz_Fy2z_a-2.0E0*2*I_KINETIC_Fxyz_Fy2z_a;
  abcd[385] = 4.0E0*I_KINETIC_Hx2y2z_Fy2z_aa-2.0E0*1*I_KINETIC_Fx2z_Fy2z_a;
  abcd[386] = 4.0E0*I_KINETIC_H5y_Fy2z_aa-2.0E0*3*I_KINETIC_F3y_Fy2z_a-2.0E0*4*I_KINETIC_F3y_Fy2z_a+3*2*I_KINETIC_Py_Fy2z;
  abcd[387] = 4.0E0*I_KINETIC_H4yz_Fy2z_aa-2.0E0*2*I_KINETIC_F2yz_Fy2z_a-2.0E0*3*I_KINETIC_F2yz_Fy2z_a+2*1*I_KINETIC_Pz_Fy2z;
  abcd[388] = 4.0E0*I_KINETIC_H3y2z_Fy2z_aa-2.0E0*1*I_KINETIC_Fy2z_Fy2z_a-2.0E0*2*I_KINETIC_Fy2z_Fy2z_a;
  abcd[389] = 4.0E0*I_KINETIC_H2y3z_Fy2z_aa-2.0E0*1*I_KINETIC_F3z_Fy2z_a;
  abcd[390] = 4.0E0*I_KINETIC_H3x2y_F3z_aa-2.0E0*1*I_KINETIC_F3x_F3z_a;
  abcd[391] = 4.0E0*I_KINETIC_H2x3y_F3z_aa-2.0E0*1*I_KINETIC_F2xy_F3z_a-2.0E0*2*I_KINETIC_F2xy_F3z_a;
  abcd[392] = 4.0E0*I_KINETIC_H2x2yz_F3z_aa-2.0E0*1*I_KINETIC_F2xz_F3z_a;
  abcd[393] = 4.0E0*I_KINETIC_Hx4y_F3z_aa-2.0E0*2*I_KINETIC_Fx2y_F3z_a-2.0E0*3*I_KINETIC_Fx2y_F3z_a+2*1*I_KINETIC_Px_F3z;
  abcd[394] = 4.0E0*I_KINETIC_Hx3yz_F3z_aa-2.0E0*1*I_KINETIC_Fxyz_F3z_a-2.0E0*2*I_KINETIC_Fxyz_F3z_a;
  abcd[395] = 4.0E0*I_KINETIC_Hx2y2z_F3z_aa-2.0E0*1*I_KINETIC_Fx2z_F3z_a;
  abcd[396] = 4.0E0*I_KINETIC_H5y_F3z_aa-2.0E0*3*I_KINETIC_F3y_F3z_a-2.0E0*4*I_KINETIC_F3y_F3z_a+3*2*I_KINETIC_Py_F3z;
  abcd[397] = 4.0E0*I_KINETIC_H4yz_F3z_aa-2.0E0*2*I_KINETIC_F2yz_F3z_a-2.0E0*3*I_KINETIC_F2yz_F3z_a+2*1*I_KINETIC_Pz_F3z;
  abcd[398] = 4.0E0*I_KINETIC_H3y2z_F3z_aa-2.0E0*1*I_KINETIC_Fy2z_F3z_a-2.0E0*2*I_KINETIC_Fy2z_F3z_a;
  abcd[399] = 4.0E0*I_KINETIC_H2y3z_F3z_aa-2.0E0*1*I_KINETIC_F3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_F_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_F_aa
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_P_F
   ************************************************************/
  abcd[400] = 4.0E0*I_KINETIC_H3xyz_F3x_aa;
  abcd[401] = 4.0E0*I_KINETIC_H2x2yz_F3x_aa-2.0E0*1*I_KINETIC_F2xz_F3x_a;
  abcd[402] = 4.0E0*I_KINETIC_H2xy2z_F3x_aa-2.0E0*1*I_KINETIC_F2xy_F3x_a;
  abcd[403] = 4.0E0*I_KINETIC_Hx3yz_F3x_aa-2.0E0*2*I_KINETIC_Fxyz_F3x_a;
  abcd[404] = 4.0E0*I_KINETIC_Hx2y2z_F3x_aa-2.0E0*1*I_KINETIC_Fx2y_F3x_a-2.0E0*1*I_KINETIC_Fx2z_F3x_a+1*I_KINETIC_Px_F3x;
  abcd[405] = 4.0E0*I_KINETIC_Hxy3z_F3x_aa-2.0E0*2*I_KINETIC_Fxyz_F3x_a;
  abcd[406] = 4.0E0*I_KINETIC_H4yz_F3x_aa-2.0E0*3*I_KINETIC_F2yz_F3x_a;
  abcd[407] = 4.0E0*I_KINETIC_H3y2z_F3x_aa-2.0E0*1*I_KINETIC_F3y_F3x_a-2.0E0*2*I_KINETIC_Fy2z_F3x_a+2*1*I_KINETIC_Py_F3x;
  abcd[408] = 4.0E0*I_KINETIC_H2y3z_F3x_aa-2.0E0*2*I_KINETIC_F2yz_F3x_a-2.0E0*1*I_KINETIC_F3z_F3x_a+2*I_KINETIC_Pz_F3x;
  abcd[409] = 4.0E0*I_KINETIC_Hy4z_F3x_aa-2.0E0*3*I_KINETIC_Fy2z_F3x_a;
  abcd[410] = 4.0E0*I_KINETIC_H3xyz_F2xy_aa;
  abcd[411] = 4.0E0*I_KINETIC_H2x2yz_F2xy_aa-2.0E0*1*I_KINETIC_F2xz_F2xy_a;
  abcd[412] = 4.0E0*I_KINETIC_H2xy2z_F2xy_aa-2.0E0*1*I_KINETIC_F2xy_F2xy_a;
  abcd[413] = 4.0E0*I_KINETIC_Hx3yz_F2xy_aa-2.0E0*2*I_KINETIC_Fxyz_F2xy_a;
  abcd[414] = 4.0E0*I_KINETIC_Hx2y2z_F2xy_aa-2.0E0*1*I_KINETIC_Fx2y_F2xy_a-2.0E0*1*I_KINETIC_Fx2z_F2xy_a+1*I_KINETIC_Px_F2xy;
  abcd[415] = 4.0E0*I_KINETIC_Hxy3z_F2xy_aa-2.0E0*2*I_KINETIC_Fxyz_F2xy_a;
  abcd[416] = 4.0E0*I_KINETIC_H4yz_F2xy_aa-2.0E0*3*I_KINETIC_F2yz_F2xy_a;
  abcd[417] = 4.0E0*I_KINETIC_H3y2z_F2xy_aa-2.0E0*1*I_KINETIC_F3y_F2xy_a-2.0E0*2*I_KINETIC_Fy2z_F2xy_a+2*1*I_KINETIC_Py_F2xy;
  abcd[418] = 4.0E0*I_KINETIC_H2y3z_F2xy_aa-2.0E0*2*I_KINETIC_F2yz_F2xy_a-2.0E0*1*I_KINETIC_F3z_F2xy_a+2*I_KINETIC_Pz_F2xy;
  abcd[419] = 4.0E0*I_KINETIC_Hy4z_F2xy_aa-2.0E0*3*I_KINETIC_Fy2z_F2xy_a;
  abcd[420] = 4.0E0*I_KINETIC_H3xyz_F2xz_aa;
  abcd[421] = 4.0E0*I_KINETIC_H2x2yz_F2xz_aa-2.0E0*1*I_KINETIC_F2xz_F2xz_a;
  abcd[422] = 4.0E0*I_KINETIC_H2xy2z_F2xz_aa-2.0E0*1*I_KINETIC_F2xy_F2xz_a;
  abcd[423] = 4.0E0*I_KINETIC_Hx3yz_F2xz_aa-2.0E0*2*I_KINETIC_Fxyz_F2xz_a;
  abcd[424] = 4.0E0*I_KINETIC_Hx2y2z_F2xz_aa-2.0E0*1*I_KINETIC_Fx2y_F2xz_a-2.0E0*1*I_KINETIC_Fx2z_F2xz_a+1*I_KINETIC_Px_F2xz;
  abcd[425] = 4.0E0*I_KINETIC_Hxy3z_F2xz_aa-2.0E0*2*I_KINETIC_Fxyz_F2xz_a;
  abcd[426] = 4.0E0*I_KINETIC_H4yz_F2xz_aa-2.0E0*3*I_KINETIC_F2yz_F2xz_a;
  abcd[427] = 4.0E0*I_KINETIC_H3y2z_F2xz_aa-2.0E0*1*I_KINETIC_F3y_F2xz_a-2.0E0*2*I_KINETIC_Fy2z_F2xz_a+2*1*I_KINETIC_Py_F2xz;
  abcd[428] = 4.0E0*I_KINETIC_H2y3z_F2xz_aa-2.0E0*2*I_KINETIC_F2yz_F2xz_a-2.0E0*1*I_KINETIC_F3z_F2xz_a+2*I_KINETIC_Pz_F2xz;
  abcd[429] = 4.0E0*I_KINETIC_Hy4z_F2xz_aa-2.0E0*3*I_KINETIC_Fy2z_F2xz_a;
  abcd[430] = 4.0E0*I_KINETIC_H3xyz_Fx2y_aa;
  abcd[431] = 4.0E0*I_KINETIC_H2x2yz_Fx2y_aa-2.0E0*1*I_KINETIC_F2xz_Fx2y_a;
  abcd[432] = 4.0E0*I_KINETIC_H2xy2z_Fx2y_aa-2.0E0*1*I_KINETIC_F2xy_Fx2y_a;
  abcd[433] = 4.0E0*I_KINETIC_Hx3yz_Fx2y_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2y_a;
  abcd[434] = 4.0E0*I_KINETIC_Hx2y2z_Fx2y_aa-2.0E0*1*I_KINETIC_Fx2y_Fx2y_a-2.0E0*1*I_KINETIC_Fx2z_Fx2y_a+1*I_KINETIC_Px_Fx2y;
  abcd[435] = 4.0E0*I_KINETIC_Hxy3z_Fx2y_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2y_a;
  abcd[436] = 4.0E0*I_KINETIC_H4yz_Fx2y_aa-2.0E0*3*I_KINETIC_F2yz_Fx2y_a;
  abcd[437] = 4.0E0*I_KINETIC_H3y2z_Fx2y_aa-2.0E0*1*I_KINETIC_F3y_Fx2y_a-2.0E0*2*I_KINETIC_Fy2z_Fx2y_a+2*1*I_KINETIC_Py_Fx2y;
  abcd[438] = 4.0E0*I_KINETIC_H2y3z_Fx2y_aa-2.0E0*2*I_KINETIC_F2yz_Fx2y_a-2.0E0*1*I_KINETIC_F3z_Fx2y_a+2*I_KINETIC_Pz_Fx2y;
  abcd[439] = 4.0E0*I_KINETIC_Hy4z_Fx2y_aa-2.0E0*3*I_KINETIC_Fy2z_Fx2y_a;
  abcd[440] = 4.0E0*I_KINETIC_H3xyz_Fxyz_aa;
  abcd[441] = 4.0E0*I_KINETIC_H2x2yz_Fxyz_aa-2.0E0*1*I_KINETIC_F2xz_Fxyz_a;
  abcd[442] = 4.0E0*I_KINETIC_H2xy2z_Fxyz_aa-2.0E0*1*I_KINETIC_F2xy_Fxyz_a;
  abcd[443] = 4.0E0*I_KINETIC_Hx3yz_Fxyz_aa-2.0E0*2*I_KINETIC_Fxyz_Fxyz_a;
  abcd[444] = 4.0E0*I_KINETIC_Hx2y2z_Fxyz_aa-2.0E0*1*I_KINETIC_Fx2y_Fxyz_a-2.0E0*1*I_KINETIC_Fx2z_Fxyz_a+1*I_KINETIC_Px_Fxyz;
  abcd[445] = 4.0E0*I_KINETIC_Hxy3z_Fxyz_aa-2.0E0*2*I_KINETIC_Fxyz_Fxyz_a;
  abcd[446] = 4.0E0*I_KINETIC_H4yz_Fxyz_aa-2.0E0*3*I_KINETIC_F2yz_Fxyz_a;
  abcd[447] = 4.0E0*I_KINETIC_H3y2z_Fxyz_aa-2.0E0*1*I_KINETIC_F3y_Fxyz_a-2.0E0*2*I_KINETIC_Fy2z_Fxyz_a+2*1*I_KINETIC_Py_Fxyz;
  abcd[448] = 4.0E0*I_KINETIC_H2y3z_Fxyz_aa-2.0E0*2*I_KINETIC_F2yz_Fxyz_a-2.0E0*1*I_KINETIC_F3z_Fxyz_a+2*I_KINETIC_Pz_Fxyz;
  abcd[449] = 4.0E0*I_KINETIC_Hy4z_Fxyz_aa-2.0E0*3*I_KINETIC_Fy2z_Fxyz_a;
  abcd[450] = 4.0E0*I_KINETIC_H3xyz_Fx2z_aa;
  abcd[451] = 4.0E0*I_KINETIC_H2x2yz_Fx2z_aa-2.0E0*1*I_KINETIC_F2xz_Fx2z_a;
  abcd[452] = 4.0E0*I_KINETIC_H2xy2z_Fx2z_aa-2.0E0*1*I_KINETIC_F2xy_Fx2z_a;
  abcd[453] = 4.0E0*I_KINETIC_Hx3yz_Fx2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2z_a;
  abcd[454] = 4.0E0*I_KINETIC_Hx2y2z_Fx2z_aa-2.0E0*1*I_KINETIC_Fx2y_Fx2z_a-2.0E0*1*I_KINETIC_Fx2z_Fx2z_a+1*I_KINETIC_Px_Fx2z;
  abcd[455] = 4.0E0*I_KINETIC_Hxy3z_Fx2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fx2z_a;
  abcd[456] = 4.0E0*I_KINETIC_H4yz_Fx2z_aa-2.0E0*3*I_KINETIC_F2yz_Fx2z_a;
  abcd[457] = 4.0E0*I_KINETIC_H3y2z_Fx2z_aa-2.0E0*1*I_KINETIC_F3y_Fx2z_a-2.0E0*2*I_KINETIC_Fy2z_Fx2z_a+2*1*I_KINETIC_Py_Fx2z;
  abcd[458] = 4.0E0*I_KINETIC_H2y3z_Fx2z_aa-2.0E0*2*I_KINETIC_F2yz_Fx2z_a-2.0E0*1*I_KINETIC_F3z_Fx2z_a+2*I_KINETIC_Pz_Fx2z;
  abcd[459] = 4.0E0*I_KINETIC_Hy4z_Fx2z_aa-2.0E0*3*I_KINETIC_Fy2z_Fx2z_a;
  abcd[460] = 4.0E0*I_KINETIC_H3xyz_F3y_aa;
  abcd[461] = 4.0E0*I_KINETIC_H2x2yz_F3y_aa-2.0E0*1*I_KINETIC_F2xz_F3y_a;
  abcd[462] = 4.0E0*I_KINETIC_H2xy2z_F3y_aa-2.0E0*1*I_KINETIC_F2xy_F3y_a;
  abcd[463] = 4.0E0*I_KINETIC_Hx3yz_F3y_aa-2.0E0*2*I_KINETIC_Fxyz_F3y_a;
  abcd[464] = 4.0E0*I_KINETIC_Hx2y2z_F3y_aa-2.0E0*1*I_KINETIC_Fx2y_F3y_a-2.0E0*1*I_KINETIC_Fx2z_F3y_a+1*I_KINETIC_Px_F3y;
  abcd[465] = 4.0E0*I_KINETIC_Hxy3z_F3y_aa-2.0E0*2*I_KINETIC_Fxyz_F3y_a;
  abcd[466] = 4.0E0*I_KINETIC_H4yz_F3y_aa-2.0E0*3*I_KINETIC_F2yz_F3y_a;
  abcd[467] = 4.0E0*I_KINETIC_H3y2z_F3y_aa-2.0E0*1*I_KINETIC_F3y_F3y_a-2.0E0*2*I_KINETIC_Fy2z_F3y_a+2*1*I_KINETIC_Py_F3y;
  abcd[468] = 4.0E0*I_KINETIC_H2y3z_F3y_aa-2.0E0*2*I_KINETIC_F2yz_F3y_a-2.0E0*1*I_KINETIC_F3z_F3y_a+2*I_KINETIC_Pz_F3y;
  abcd[469] = 4.0E0*I_KINETIC_Hy4z_F3y_aa-2.0E0*3*I_KINETIC_Fy2z_F3y_a;
  abcd[470] = 4.0E0*I_KINETIC_H3xyz_F2yz_aa;
  abcd[471] = 4.0E0*I_KINETIC_H2x2yz_F2yz_aa-2.0E0*1*I_KINETIC_F2xz_F2yz_a;
  abcd[472] = 4.0E0*I_KINETIC_H2xy2z_F2yz_aa-2.0E0*1*I_KINETIC_F2xy_F2yz_a;
  abcd[473] = 4.0E0*I_KINETIC_Hx3yz_F2yz_aa-2.0E0*2*I_KINETIC_Fxyz_F2yz_a;
  abcd[474] = 4.0E0*I_KINETIC_Hx2y2z_F2yz_aa-2.0E0*1*I_KINETIC_Fx2y_F2yz_a-2.0E0*1*I_KINETIC_Fx2z_F2yz_a+1*I_KINETIC_Px_F2yz;
  abcd[475] = 4.0E0*I_KINETIC_Hxy3z_F2yz_aa-2.0E0*2*I_KINETIC_Fxyz_F2yz_a;
  abcd[476] = 4.0E0*I_KINETIC_H4yz_F2yz_aa-2.0E0*3*I_KINETIC_F2yz_F2yz_a;
  abcd[477] = 4.0E0*I_KINETIC_H3y2z_F2yz_aa-2.0E0*1*I_KINETIC_F3y_F2yz_a-2.0E0*2*I_KINETIC_Fy2z_F2yz_a+2*1*I_KINETIC_Py_F2yz;
  abcd[478] = 4.0E0*I_KINETIC_H2y3z_F2yz_aa-2.0E0*2*I_KINETIC_F2yz_F2yz_a-2.0E0*1*I_KINETIC_F3z_F2yz_a+2*I_KINETIC_Pz_F2yz;
  abcd[479] = 4.0E0*I_KINETIC_Hy4z_F2yz_aa-2.0E0*3*I_KINETIC_Fy2z_F2yz_a;
  abcd[480] = 4.0E0*I_KINETIC_H3xyz_Fy2z_aa;
  abcd[481] = 4.0E0*I_KINETIC_H2x2yz_Fy2z_aa-2.0E0*1*I_KINETIC_F2xz_Fy2z_a;
  abcd[482] = 4.0E0*I_KINETIC_H2xy2z_Fy2z_aa-2.0E0*1*I_KINETIC_F2xy_Fy2z_a;
  abcd[483] = 4.0E0*I_KINETIC_Hx3yz_Fy2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fy2z_a;
  abcd[484] = 4.0E0*I_KINETIC_Hx2y2z_Fy2z_aa-2.0E0*1*I_KINETIC_Fx2y_Fy2z_a-2.0E0*1*I_KINETIC_Fx2z_Fy2z_a+1*I_KINETIC_Px_Fy2z;
  abcd[485] = 4.0E0*I_KINETIC_Hxy3z_Fy2z_aa-2.0E0*2*I_KINETIC_Fxyz_Fy2z_a;
  abcd[486] = 4.0E0*I_KINETIC_H4yz_Fy2z_aa-2.0E0*3*I_KINETIC_F2yz_Fy2z_a;
  abcd[487] = 4.0E0*I_KINETIC_H3y2z_Fy2z_aa-2.0E0*1*I_KINETIC_F3y_Fy2z_a-2.0E0*2*I_KINETIC_Fy2z_Fy2z_a+2*1*I_KINETIC_Py_Fy2z;
  abcd[488] = 4.0E0*I_KINETIC_H2y3z_Fy2z_aa-2.0E0*2*I_KINETIC_F2yz_Fy2z_a-2.0E0*1*I_KINETIC_F3z_Fy2z_a+2*I_KINETIC_Pz_Fy2z;
  abcd[489] = 4.0E0*I_KINETIC_Hy4z_Fy2z_aa-2.0E0*3*I_KINETIC_Fy2z_Fy2z_a;
  abcd[490] = 4.0E0*I_KINETIC_H3xyz_F3z_aa;
  abcd[491] = 4.0E0*I_KINETIC_H2x2yz_F3z_aa-2.0E0*1*I_KINETIC_F2xz_F3z_a;
  abcd[492] = 4.0E0*I_KINETIC_H2xy2z_F3z_aa-2.0E0*1*I_KINETIC_F2xy_F3z_a;
  abcd[493] = 4.0E0*I_KINETIC_Hx3yz_F3z_aa-2.0E0*2*I_KINETIC_Fxyz_F3z_a;
  abcd[494] = 4.0E0*I_KINETIC_Hx2y2z_F3z_aa-2.0E0*1*I_KINETIC_Fx2y_F3z_a-2.0E0*1*I_KINETIC_Fx2z_F3z_a+1*I_KINETIC_Px_F3z;
  abcd[495] = 4.0E0*I_KINETIC_Hxy3z_F3z_aa-2.0E0*2*I_KINETIC_Fxyz_F3z_a;
  abcd[496] = 4.0E0*I_KINETIC_H4yz_F3z_aa-2.0E0*3*I_KINETIC_F2yz_F3z_a;
  abcd[497] = 4.0E0*I_KINETIC_H3y2z_F3z_aa-2.0E0*1*I_KINETIC_F3y_F3z_a-2.0E0*2*I_KINETIC_Fy2z_F3z_a+2*1*I_KINETIC_Py_F3z;
  abcd[498] = 4.0E0*I_KINETIC_H2y3z_F3z_aa-2.0E0*2*I_KINETIC_F2yz_F3z_a-2.0E0*1*I_KINETIC_F3z_F3z_a+2*I_KINETIC_Pz_F3z;
  abcd[499] = 4.0E0*I_KINETIC_Hy4z_F3z_aa-2.0E0*3*I_KINETIC_Fy2z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_F_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_F_aa
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_F_F_a
   * RHS shell quartet name: SQ_KINETIC_P_F
   ************************************************************/
  abcd[500] = 4.0E0*I_KINETIC_H3x2z_F3x_aa-2.0E0*1*I_KINETIC_F3x_F3x_a;
  abcd[501] = 4.0E0*I_KINETIC_H2xy2z_F3x_aa-2.0E0*1*I_KINETIC_F2xy_F3x_a;
  abcd[502] = 4.0E0*I_KINETIC_H2x3z_F3x_aa-2.0E0*1*I_KINETIC_F2xz_F3x_a-2.0E0*2*I_KINETIC_F2xz_F3x_a;
  abcd[503] = 4.0E0*I_KINETIC_Hx2y2z_F3x_aa-2.0E0*1*I_KINETIC_Fx2y_F3x_a;
  abcd[504] = 4.0E0*I_KINETIC_Hxy3z_F3x_aa-2.0E0*1*I_KINETIC_Fxyz_F3x_a-2.0E0*2*I_KINETIC_Fxyz_F3x_a;
  abcd[505] = 4.0E0*I_KINETIC_Hx4z_F3x_aa-2.0E0*2*I_KINETIC_Fx2z_F3x_a-2.0E0*3*I_KINETIC_Fx2z_F3x_a+2*1*I_KINETIC_Px_F3x;
  abcd[506] = 4.0E0*I_KINETIC_H3y2z_F3x_aa-2.0E0*1*I_KINETIC_F3y_F3x_a;
  abcd[507] = 4.0E0*I_KINETIC_H2y3z_F3x_aa-2.0E0*1*I_KINETIC_F2yz_F3x_a-2.0E0*2*I_KINETIC_F2yz_F3x_a;
  abcd[508] = 4.0E0*I_KINETIC_Hy4z_F3x_aa-2.0E0*2*I_KINETIC_Fy2z_F3x_a-2.0E0*3*I_KINETIC_Fy2z_F3x_a+2*1*I_KINETIC_Py_F3x;
  abcd[509] = 4.0E0*I_KINETIC_H5z_F3x_aa-2.0E0*3*I_KINETIC_F3z_F3x_a-2.0E0*4*I_KINETIC_F3z_F3x_a+3*2*I_KINETIC_Pz_F3x;
  abcd[510] = 4.0E0*I_KINETIC_H3x2z_F2xy_aa-2.0E0*1*I_KINETIC_F3x_F2xy_a;
  abcd[511] = 4.0E0*I_KINETIC_H2xy2z_F2xy_aa-2.0E0*1*I_KINETIC_F2xy_F2xy_a;
  abcd[512] = 4.0E0*I_KINETIC_H2x3z_F2xy_aa-2.0E0*1*I_KINETIC_F2xz_F2xy_a-2.0E0*2*I_KINETIC_F2xz_F2xy_a;
  abcd[513] = 4.0E0*I_KINETIC_Hx2y2z_F2xy_aa-2.0E0*1*I_KINETIC_Fx2y_F2xy_a;
  abcd[514] = 4.0E0*I_KINETIC_Hxy3z_F2xy_aa-2.0E0*1*I_KINETIC_Fxyz_F2xy_a-2.0E0*2*I_KINETIC_Fxyz_F2xy_a;
  abcd[515] = 4.0E0*I_KINETIC_Hx4z_F2xy_aa-2.0E0*2*I_KINETIC_Fx2z_F2xy_a-2.0E0*3*I_KINETIC_Fx2z_F2xy_a+2*1*I_KINETIC_Px_F2xy;
  abcd[516] = 4.0E0*I_KINETIC_H3y2z_F2xy_aa-2.0E0*1*I_KINETIC_F3y_F2xy_a;
  abcd[517] = 4.0E0*I_KINETIC_H2y3z_F2xy_aa-2.0E0*1*I_KINETIC_F2yz_F2xy_a-2.0E0*2*I_KINETIC_F2yz_F2xy_a;
  abcd[518] = 4.0E0*I_KINETIC_Hy4z_F2xy_aa-2.0E0*2*I_KINETIC_Fy2z_F2xy_a-2.0E0*3*I_KINETIC_Fy2z_F2xy_a+2*1*I_KINETIC_Py_F2xy;
  abcd[519] = 4.0E0*I_KINETIC_H5z_F2xy_aa-2.0E0*3*I_KINETIC_F3z_F2xy_a-2.0E0*4*I_KINETIC_F3z_F2xy_a+3*2*I_KINETIC_Pz_F2xy;
  abcd[520] = 4.0E0*I_KINETIC_H3x2z_F2xz_aa-2.0E0*1*I_KINETIC_F3x_F2xz_a;
  abcd[521] = 4.0E0*I_KINETIC_H2xy2z_F2xz_aa-2.0E0*1*I_KINETIC_F2xy_F2xz_a;
  abcd[522] = 4.0E0*I_KINETIC_H2x3z_F2xz_aa-2.0E0*1*I_KINETIC_F2xz_F2xz_a-2.0E0*2*I_KINETIC_F2xz_F2xz_a;
  abcd[523] = 4.0E0*I_KINETIC_Hx2y2z_F2xz_aa-2.0E0*1*I_KINETIC_Fx2y_F2xz_a;
  abcd[524] = 4.0E0*I_KINETIC_Hxy3z_F2xz_aa-2.0E0*1*I_KINETIC_Fxyz_F2xz_a-2.0E0*2*I_KINETIC_Fxyz_F2xz_a;
  abcd[525] = 4.0E0*I_KINETIC_Hx4z_F2xz_aa-2.0E0*2*I_KINETIC_Fx2z_F2xz_a-2.0E0*3*I_KINETIC_Fx2z_F2xz_a+2*1*I_KINETIC_Px_F2xz;
  abcd[526] = 4.0E0*I_KINETIC_H3y2z_F2xz_aa-2.0E0*1*I_KINETIC_F3y_F2xz_a;
  abcd[527] = 4.0E0*I_KINETIC_H2y3z_F2xz_aa-2.0E0*1*I_KINETIC_F2yz_F2xz_a-2.0E0*2*I_KINETIC_F2yz_F2xz_a;
  abcd[528] = 4.0E0*I_KINETIC_Hy4z_F2xz_aa-2.0E0*2*I_KINETIC_Fy2z_F2xz_a-2.0E0*3*I_KINETIC_Fy2z_F2xz_a+2*1*I_KINETIC_Py_F2xz;
  abcd[529] = 4.0E0*I_KINETIC_H5z_F2xz_aa-2.0E0*3*I_KINETIC_F3z_F2xz_a-2.0E0*4*I_KINETIC_F3z_F2xz_a+3*2*I_KINETIC_Pz_F2xz;
  abcd[530] = 4.0E0*I_KINETIC_H3x2z_Fx2y_aa-2.0E0*1*I_KINETIC_F3x_Fx2y_a;
  abcd[531] = 4.0E0*I_KINETIC_H2xy2z_Fx2y_aa-2.0E0*1*I_KINETIC_F2xy_Fx2y_a;
  abcd[532] = 4.0E0*I_KINETIC_H2x3z_Fx2y_aa-2.0E0*1*I_KINETIC_F2xz_Fx2y_a-2.0E0*2*I_KINETIC_F2xz_Fx2y_a;
  abcd[533] = 4.0E0*I_KINETIC_Hx2y2z_Fx2y_aa-2.0E0*1*I_KINETIC_Fx2y_Fx2y_a;
  abcd[534] = 4.0E0*I_KINETIC_Hxy3z_Fx2y_aa-2.0E0*1*I_KINETIC_Fxyz_Fx2y_a-2.0E0*2*I_KINETIC_Fxyz_Fx2y_a;
  abcd[535] = 4.0E0*I_KINETIC_Hx4z_Fx2y_aa-2.0E0*2*I_KINETIC_Fx2z_Fx2y_a-2.0E0*3*I_KINETIC_Fx2z_Fx2y_a+2*1*I_KINETIC_Px_Fx2y;
  abcd[536] = 4.0E0*I_KINETIC_H3y2z_Fx2y_aa-2.0E0*1*I_KINETIC_F3y_Fx2y_a;
  abcd[537] = 4.0E0*I_KINETIC_H2y3z_Fx2y_aa-2.0E0*1*I_KINETIC_F2yz_Fx2y_a-2.0E0*2*I_KINETIC_F2yz_Fx2y_a;
  abcd[538] = 4.0E0*I_KINETIC_Hy4z_Fx2y_aa-2.0E0*2*I_KINETIC_Fy2z_Fx2y_a-2.0E0*3*I_KINETIC_Fy2z_Fx2y_a+2*1*I_KINETIC_Py_Fx2y;
  abcd[539] = 4.0E0*I_KINETIC_H5z_Fx2y_aa-2.0E0*3*I_KINETIC_F3z_Fx2y_a-2.0E0*4*I_KINETIC_F3z_Fx2y_a+3*2*I_KINETIC_Pz_Fx2y;
  abcd[540] = 4.0E0*I_KINETIC_H3x2z_Fxyz_aa-2.0E0*1*I_KINETIC_F3x_Fxyz_a;
  abcd[541] = 4.0E0*I_KINETIC_H2xy2z_Fxyz_aa-2.0E0*1*I_KINETIC_F2xy_Fxyz_a;
  abcd[542] = 4.0E0*I_KINETIC_H2x3z_Fxyz_aa-2.0E0*1*I_KINETIC_F2xz_Fxyz_a-2.0E0*2*I_KINETIC_F2xz_Fxyz_a;
  abcd[543] = 4.0E0*I_KINETIC_Hx2y2z_Fxyz_aa-2.0E0*1*I_KINETIC_Fx2y_Fxyz_a;
  abcd[544] = 4.0E0*I_KINETIC_Hxy3z_Fxyz_aa-2.0E0*1*I_KINETIC_Fxyz_Fxyz_a-2.0E0*2*I_KINETIC_Fxyz_Fxyz_a;
  abcd[545] = 4.0E0*I_KINETIC_Hx4z_Fxyz_aa-2.0E0*2*I_KINETIC_Fx2z_Fxyz_a-2.0E0*3*I_KINETIC_Fx2z_Fxyz_a+2*1*I_KINETIC_Px_Fxyz;
  abcd[546] = 4.0E0*I_KINETIC_H3y2z_Fxyz_aa-2.0E0*1*I_KINETIC_F3y_Fxyz_a;
  abcd[547] = 4.0E0*I_KINETIC_H2y3z_Fxyz_aa-2.0E0*1*I_KINETIC_F2yz_Fxyz_a-2.0E0*2*I_KINETIC_F2yz_Fxyz_a;
  abcd[548] = 4.0E0*I_KINETIC_Hy4z_Fxyz_aa-2.0E0*2*I_KINETIC_Fy2z_Fxyz_a-2.0E0*3*I_KINETIC_Fy2z_Fxyz_a+2*1*I_KINETIC_Py_Fxyz;
  abcd[549] = 4.0E0*I_KINETIC_H5z_Fxyz_aa-2.0E0*3*I_KINETIC_F3z_Fxyz_a-2.0E0*4*I_KINETIC_F3z_Fxyz_a+3*2*I_KINETIC_Pz_Fxyz;
  abcd[550] = 4.0E0*I_KINETIC_H3x2z_Fx2z_aa-2.0E0*1*I_KINETIC_F3x_Fx2z_a;
  abcd[551] = 4.0E0*I_KINETIC_H2xy2z_Fx2z_aa-2.0E0*1*I_KINETIC_F2xy_Fx2z_a;
  abcd[552] = 4.0E0*I_KINETIC_H2x3z_Fx2z_aa-2.0E0*1*I_KINETIC_F2xz_Fx2z_a-2.0E0*2*I_KINETIC_F2xz_Fx2z_a;
  abcd[553] = 4.0E0*I_KINETIC_Hx2y2z_Fx2z_aa-2.0E0*1*I_KINETIC_Fx2y_Fx2z_a;
  abcd[554] = 4.0E0*I_KINETIC_Hxy3z_Fx2z_aa-2.0E0*1*I_KINETIC_Fxyz_Fx2z_a-2.0E0*2*I_KINETIC_Fxyz_Fx2z_a;
  abcd[555] = 4.0E0*I_KINETIC_Hx4z_Fx2z_aa-2.0E0*2*I_KINETIC_Fx2z_Fx2z_a-2.0E0*3*I_KINETIC_Fx2z_Fx2z_a+2*1*I_KINETIC_Px_Fx2z;
  abcd[556] = 4.0E0*I_KINETIC_H3y2z_Fx2z_aa-2.0E0*1*I_KINETIC_F3y_Fx2z_a;
  abcd[557] = 4.0E0*I_KINETIC_H2y3z_Fx2z_aa-2.0E0*1*I_KINETIC_F2yz_Fx2z_a-2.0E0*2*I_KINETIC_F2yz_Fx2z_a;
  abcd[558] = 4.0E0*I_KINETIC_Hy4z_Fx2z_aa-2.0E0*2*I_KINETIC_Fy2z_Fx2z_a-2.0E0*3*I_KINETIC_Fy2z_Fx2z_a+2*1*I_KINETIC_Py_Fx2z;
  abcd[559] = 4.0E0*I_KINETIC_H5z_Fx2z_aa-2.0E0*3*I_KINETIC_F3z_Fx2z_a-2.0E0*4*I_KINETIC_F3z_Fx2z_a+3*2*I_KINETIC_Pz_Fx2z;
  abcd[560] = 4.0E0*I_KINETIC_H3x2z_F3y_aa-2.0E0*1*I_KINETIC_F3x_F3y_a;
  abcd[561] = 4.0E0*I_KINETIC_H2xy2z_F3y_aa-2.0E0*1*I_KINETIC_F2xy_F3y_a;
  abcd[562] = 4.0E0*I_KINETIC_H2x3z_F3y_aa-2.0E0*1*I_KINETIC_F2xz_F3y_a-2.0E0*2*I_KINETIC_F2xz_F3y_a;
  abcd[563] = 4.0E0*I_KINETIC_Hx2y2z_F3y_aa-2.0E0*1*I_KINETIC_Fx2y_F3y_a;
  abcd[564] = 4.0E0*I_KINETIC_Hxy3z_F3y_aa-2.0E0*1*I_KINETIC_Fxyz_F3y_a-2.0E0*2*I_KINETIC_Fxyz_F3y_a;
  abcd[565] = 4.0E0*I_KINETIC_Hx4z_F3y_aa-2.0E0*2*I_KINETIC_Fx2z_F3y_a-2.0E0*3*I_KINETIC_Fx2z_F3y_a+2*1*I_KINETIC_Px_F3y;
  abcd[566] = 4.0E0*I_KINETIC_H3y2z_F3y_aa-2.0E0*1*I_KINETIC_F3y_F3y_a;
  abcd[567] = 4.0E0*I_KINETIC_H2y3z_F3y_aa-2.0E0*1*I_KINETIC_F2yz_F3y_a-2.0E0*2*I_KINETIC_F2yz_F3y_a;
  abcd[568] = 4.0E0*I_KINETIC_Hy4z_F3y_aa-2.0E0*2*I_KINETIC_Fy2z_F3y_a-2.0E0*3*I_KINETIC_Fy2z_F3y_a+2*1*I_KINETIC_Py_F3y;
  abcd[569] = 4.0E0*I_KINETIC_H5z_F3y_aa-2.0E0*3*I_KINETIC_F3z_F3y_a-2.0E0*4*I_KINETIC_F3z_F3y_a+3*2*I_KINETIC_Pz_F3y;
  abcd[570] = 4.0E0*I_KINETIC_H3x2z_F2yz_aa-2.0E0*1*I_KINETIC_F3x_F2yz_a;
  abcd[571] = 4.0E0*I_KINETIC_H2xy2z_F2yz_aa-2.0E0*1*I_KINETIC_F2xy_F2yz_a;
  abcd[572] = 4.0E0*I_KINETIC_H2x3z_F2yz_aa-2.0E0*1*I_KINETIC_F2xz_F2yz_a-2.0E0*2*I_KINETIC_F2xz_F2yz_a;
  abcd[573] = 4.0E0*I_KINETIC_Hx2y2z_F2yz_aa-2.0E0*1*I_KINETIC_Fx2y_F2yz_a;
  abcd[574] = 4.0E0*I_KINETIC_Hxy3z_F2yz_aa-2.0E0*1*I_KINETIC_Fxyz_F2yz_a-2.0E0*2*I_KINETIC_Fxyz_F2yz_a;
  abcd[575] = 4.0E0*I_KINETIC_Hx4z_F2yz_aa-2.0E0*2*I_KINETIC_Fx2z_F2yz_a-2.0E0*3*I_KINETIC_Fx2z_F2yz_a+2*1*I_KINETIC_Px_F2yz;
  abcd[576] = 4.0E0*I_KINETIC_H3y2z_F2yz_aa-2.0E0*1*I_KINETIC_F3y_F2yz_a;
  abcd[577] = 4.0E0*I_KINETIC_H2y3z_F2yz_aa-2.0E0*1*I_KINETIC_F2yz_F2yz_a-2.0E0*2*I_KINETIC_F2yz_F2yz_a;
  abcd[578] = 4.0E0*I_KINETIC_Hy4z_F2yz_aa-2.0E0*2*I_KINETIC_Fy2z_F2yz_a-2.0E0*3*I_KINETIC_Fy2z_F2yz_a+2*1*I_KINETIC_Py_F2yz;
  abcd[579] = 4.0E0*I_KINETIC_H5z_F2yz_aa-2.0E0*3*I_KINETIC_F3z_F2yz_a-2.0E0*4*I_KINETIC_F3z_F2yz_a+3*2*I_KINETIC_Pz_F2yz;
  abcd[580] = 4.0E0*I_KINETIC_H3x2z_Fy2z_aa-2.0E0*1*I_KINETIC_F3x_Fy2z_a;
  abcd[581] = 4.0E0*I_KINETIC_H2xy2z_Fy2z_aa-2.0E0*1*I_KINETIC_F2xy_Fy2z_a;
  abcd[582] = 4.0E0*I_KINETIC_H2x3z_Fy2z_aa-2.0E0*1*I_KINETIC_F2xz_Fy2z_a-2.0E0*2*I_KINETIC_F2xz_Fy2z_a;
  abcd[583] = 4.0E0*I_KINETIC_Hx2y2z_Fy2z_aa-2.0E0*1*I_KINETIC_Fx2y_Fy2z_a;
  abcd[584] = 4.0E0*I_KINETIC_Hxy3z_Fy2z_aa-2.0E0*1*I_KINETIC_Fxyz_Fy2z_a-2.0E0*2*I_KINETIC_Fxyz_Fy2z_a;
  abcd[585] = 4.0E0*I_KINETIC_Hx4z_Fy2z_aa-2.0E0*2*I_KINETIC_Fx2z_Fy2z_a-2.0E0*3*I_KINETIC_Fx2z_Fy2z_a+2*1*I_KINETIC_Px_Fy2z;
  abcd[586] = 4.0E0*I_KINETIC_H3y2z_Fy2z_aa-2.0E0*1*I_KINETIC_F3y_Fy2z_a;
  abcd[587] = 4.0E0*I_KINETIC_H2y3z_Fy2z_aa-2.0E0*1*I_KINETIC_F2yz_Fy2z_a-2.0E0*2*I_KINETIC_F2yz_Fy2z_a;
  abcd[588] = 4.0E0*I_KINETIC_Hy4z_Fy2z_aa-2.0E0*2*I_KINETIC_Fy2z_Fy2z_a-2.0E0*3*I_KINETIC_Fy2z_Fy2z_a+2*1*I_KINETIC_Py_Fy2z;
  abcd[589] = 4.0E0*I_KINETIC_H5z_Fy2z_aa-2.0E0*3*I_KINETIC_F3z_Fy2z_a-2.0E0*4*I_KINETIC_F3z_Fy2z_a+3*2*I_KINETIC_Pz_Fy2z;
  abcd[590] = 4.0E0*I_KINETIC_H3x2z_F3z_aa-2.0E0*1*I_KINETIC_F3x_F3z_a;
  abcd[591] = 4.0E0*I_KINETIC_H2xy2z_F3z_aa-2.0E0*1*I_KINETIC_F2xy_F3z_a;
  abcd[592] = 4.0E0*I_KINETIC_H2x3z_F3z_aa-2.0E0*1*I_KINETIC_F2xz_F3z_a-2.0E0*2*I_KINETIC_F2xz_F3z_a;
  abcd[593] = 4.0E0*I_KINETIC_Hx2y2z_F3z_aa-2.0E0*1*I_KINETIC_Fx2y_F3z_a;
  abcd[594] = 4.0E0*I_KINETIC_Hxy3z_F3z_aa-2.0E0*1*I_KINETIC_Fxyz_F3z_a-2.0E0*2*I_KINETIC_Fxyz_F3z_a;
  abcd[595] = 4.0E0*I_KINETIC_Hx4z_F3z_aa-2.0E0*2*I_KINETIC_Fx2z_F3z_a-2.0E0*3*I_KINETIC_Fx2z_F3z_a+2*1*I_KINETIC_Px_F3z;
  abcd[596] = 4.0E0*I_KINETIC_H3y2z_F3z_aa-2.0E0*1*I_KINETIC_F3y_F3z_a;
  abcd[597] = 4.0E0*I_KINETIC_H2y3z_F3z_aa-2.0E0*1*I_KINETIC_F2yz_F3z_a-2.0E0*2*I_KINETIC_F2yz_F3z_a;
  abcd[598] = 4.0E0*I_KINETIC_Hy4z_F3z_aa-2.0E0*2*I_KINETIC_Fy2z_F3z_a-2.0E0*3*I_KINETIC_Fy2z_F3z_a+2*1*I_KINETIC_Py_F3z;
  abcd[599] = 4.0E0*I_KINETIC_H5z_F3z_aa-2.0E0*3*I_KINETIC_F3z_F3z_a-2.0E0*4*I_KINETIC_F3z_F3z_a+3*2*I_KINETIC_Pz_F3z;
}
