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
// BRA1 as redundant position, total RHS integrals evaluated as: 4504
// BRA2 as redundant position, total RHS integrals evaluated as: 4462
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: BRA2
//

//
// @@@@ derivative position-direction information
// BRA1
// X
// Y
// Z
// ####

void hgp_os_kinetic_g_f_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_H5x_F3x_a = 0.0E0;
  Double I_KINETIC_H4xy_F3x_a = 0.0E0;
  Double I_KINETIC_H4xz_F3x_a = 0.0E0;
  Double I_KINETIC_H3x2y_F3x_a = 0.0E0;
  Double I_KINETIC_H3xyz_F3x_a = 0.0E0;
  Double I_KINETIC_H3x2z_F3x_a = 0.0E0;
  Double I_KINETIC_H2x3y_F3x_a = 0.0E0;
  Double I_KINETIC_H2x2yz_F3x_a = 0.0E0;
  Double I_KINETIC_H2xy2z_F3x_a = 0.0E0;
  Double I_KINETIC_H2x3z_F3x_a = 0.0E0;
  Double I_KINETIC_Hx4y_F3x_a = 0.0E0;
  Double I_KINETIC_Hx3yz_F3x_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_F3x_a = 0.0E0;
  Double I_KINETIC_Hxy3z_F3x_a = 0.0E0;
  Double I_KINETIC_Hx4z_F3x_a = 0.0E0;
  Double I_KINETIC_H5y_F3x_a = 0.0E0;
  Double I_KINETIC_H4yz_F3x_a = 0.0E0;
  Double I_KINETIC_H3y2z_F3x_a = 0.0E0;
  Double I_KINETIC_H2y3z_F3x_a = 0.0E0;
  Double I_KINETIC_Hy4z_F3x_a = 0.0E0;
  Double I_KINETIC_H5z_F3x_a = 0.0E0;
  Double I_KINETIC_H5x_F2xy_a = 0.0E0;
  Double I_KINETIC_H4xy_F2xy_a = 0.0E0;
  Double I_KINETIC_H4xz_F2xy_a = 0.0E0;
  Double I_KINETIC_H3x2y_F2xy_a = 0.0E0;
  Double I_KINETIC_H3xyz_F2xy_a = 0.0E0;
  Double I_KINETIC_H3x2z_F2xy_a = 0.0E0;
  Double I_KINETIC_H2x3y_F2xy_a = 0.0E0;
  Double I_KINETIC_H2x2yz_F2xy_a = 0.0E0;
  Double I_KINETIC_H2xy2z_F2xy_a = 0.0E0;
  Double I_KINETIC_H2x3z_F2xy_a = 0.0E0;
  Double I_KINETIC_Hx4y_F2xy_a = 0.0E0;
  Double I_KINETIC_Hx3yz_F2xy_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_F2xy_a = 0.0E0;
  Double I_KINETIC_Hxy3z_F2xy_a = 0.0E0;
  Double I_KINETIC_Hx4z_F2xy_a = 0.0E0;
  Double I_KINETIC_H5y_F2xy_a = 0.0E0;
  Double I_KINETIC_H4yz_F2xy_a = 0.0E0;
  Double I_KINETIC_H3y2z_F2xy_a = 0.0E0;
  Double I_KINETIC_H2y3z_F2xy_a = 0.0E0;
  Double I_KINETIC_Hy4z_F2xy_a = 0.0E0;
  Double I_KINETIC_H5z_F2xy_a = 0.0E0;
  Double I_KINETIC_H5x_F2xz_a = 0.0E0;
  Double I_KINETIC_H4xy_F2xz_a = 0.0E0;
  Double I_KINETIC_H4xz_F2xz_a = 0.0E0;
  Double I_KINETIC_H3x2y_F2xz_a = 0.0E0;
  Double I_KINETIC_H3xyz_F2xz_a = 0.0E0;
  Double I_KINETIC_H3x2z_F2xz_a = 0.0E0;
  Double I_KINETIC_H2x3y_F2xz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_F2xz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_F2xz_a = 0.0E0;
  Double I_KINETIC_H2x3z_F2xz_a = 0.0E0;
  Double I_KINETIC_Hx4y_F2xz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_F2xz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_F2xz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_F2xz_a = 0.0E0;
  Double I_KINETIC_Hx4z_F2xz_a = 0.0E0;
  Double I_KINETIC_H5y_F2xz_a = 0.0E0;
  Double I_KINETIC_H4yz_F2xz_a = 0.0E0;
  Double I_KINETIC_H3y2z_F2xz_a = 0.0E0;
  Double I_KINETIC_H2y3z_F2xz_a = 0.0E0;
  Double I_KINETIC_Hy4z_F2xz_a = 0.0E0;
  Double I_KINETIC_H5z_F2xz_a = 0.0E0;
  Double I_KINETIC_H5x_Fx2y_a = 0.0E0;
  Double I_KINETIC_H4xy_Fx2y_a = 0.0E0;
  Double I_KINETIC_H4xz_Fx2y_a = 0.0E0;
  Double I_KINETIC_H3x2y_Fx2y_a = 0.0E0;
  Double I_KINETIC_H3xyz_Fx2y_a = 0.0E0;
  Double I_KINETIC_H3x2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_H2x3y_Fx2y_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_H2x3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Hx4y_Fx2y_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Hx4z_Fx2y_a = 0.0E0;
  Double I_KINETIC_H5y_Fx2y_a = 0.0E0;
  Double I_KINETIC_H4yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_H3y2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_H2y3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Hy4z_Fx2y_a = 0.0E0;
  Double I_KINETIC_H5z_Fx2y_a = 0.0E0;
  Double I_KINETIC_H5x_Fxyz_a = 0.0E0;
  Double I_KINETIC_H4xy_Fxyz_a = 0.0E0;
  Double I_KINETIC_H4xz_Fxyz_a = 0.0E0;
  Double I_KINETIC_H3x2y_Fxyz_a = 0.0E0;
  Double I_KINETIC_H3xyz_Fxyz_a = 0.0E0;
  Double I_KINETIC_H3x2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_H2x3y_Fxyz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_H2x3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Hx4y_Fxyz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Hx4z_Fxyz_a = 0.0E0;
  Double I_KINETIC_H5y_Fxyz_a = 0.0E0;
  Double I_KINETIC_H4yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_H3y2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_H2y3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Hy4z_Fxyz_a = 0.0E0;
  Double I_KINETIC_H5z_Fxyz_a = 0.0E0;
  Double I_KINETIC_H5x_Fx2z_a = 0.0E0;
  Double I_KINETIC_H4xy_Fx2z_a = 0.0E0;
  Double I_KINETIC_H4xz_Fx2z_a = 0.0E0;
  Double I_KINETIC_H3x2y_Fx2z_a = 0.0E0;
  Double I_KINETIC_H3xyz_Fx2z_a = 0.0E0;
  Double I_KINETIC_H3x2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_H2x3y_Fx2z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_H2x3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Hx4y_Fx2z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Hx4z_Fx2z_a = 0.0E0;
  Double I_KINETIC_H5y_Fx2z_a = 0.0E0;
  Double I_KINETIC_H4yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_H3y2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_H2y3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Hy4z_Fx2z_a = 0.0E0;
  Double I_KINETIC_H5z_Fx2z_a = 0.0E0;
  Double I_KINETIC_H5x_F3y_a = 0.0E0;
  Double I_KINETIC_H4xy_F3y_a = 0.0E0;
  Double I_KINETIC_H4xz_F3y_a = 0.0E0;
  Double I_KINETIC_H3x2y_F3y_a = 0.0E0;
  Double I_KINETIC_H3xyz_F3y_a = 0.0E0;
  Double I_KINETIC_H3x2z_F3y_a = 0.0E0;
  Double I_KINETIC_H2x3y_F3y_a = 0.0E0;
  Double I_KINETIC_H2x2yz_F3y_a = 0.0E0;
  Double I_KINETIC_H2xy2z_F3y_a = 0.0E0;
  Double I_KINETIC_H2x3z_F3y_a = 0.0E0;
  Double I_KINETIC_Hx4y_F3y_a = 0.0E0;
  Double I_KINETIC_Hx3yz_F3y_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_F3y_a = 0.0E0;
  Double I_KINETIC_Hxy3z_F3y_a = 0.0E0;
  Double I_KINETIC_Hx4z_F3y_a = 0.0E0;
  Double I_KINETIC_H5y_F3y_a = 0.0E0;
  Double I_KINETIC_H4yz_F3y_a = 0.0E0;
  Double I_KINETIC_H3y2z_F3y_a = 0.0E0;
  Double I_KINETIC_H2y3z_F3y_a = 0.0E0;
  Double I_KINETIC_Hy4z_F3y_a = 0.0E0;
  Double I_KINETIC_H5z_F3y_a = 0.0E0;
  Double I_KINETIC_H5x_F2yz_a = 0.0E0;
  Double I_KINETIC_H4xy_F2yz_a = 0.0E0;
  Double I_KINETIC_H4xz_F2yz_a = 0.0E0;
  Double I_KINETIC_H3x2y_F2yz_a = 0.0E0;
  Double I_KINETIC_H3xyz_F2yz_a = 0.0E0;
  Double I_KINETIC_H3x2z_F2yz_a = 0.0E0;
  Double I_KINETIC_H2x3y_F2yz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_F2yz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_F2yz_a = 0.0E0;
  Double I_KINETIC_H2x3z_F2yz_a = 0.0E0;
  Double I_KINETIC_Hx4y_F2yz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_F2yz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_F2yz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_F2yz_a = 0.0E0;
  Double I_KINETIC_Hx4z_F2yz_a = 0.0E0;
  Double I_KINETIC_H5y_F2yz_a = 0.0E0;
  Double I_KINETIC_H4yz_F2yz_a = 0.0E0;
  Double I_KINETIC_H3y2z_F2yz_a = 0.0E0;
  Double I_KINETIC_H2y3z_F2yz_a = 0.0E0;
  Double I_KINETIC_Hy4z_F2yz_a = 0.0E0;
  Double I_KINETIC_H5z_F2yz_a = 0.0E0;
  Double I_KINETIC_H5x_Fy2z_a = 0.0E0;
  Double I_KINETIC_H4xy_Fy2z_a = 0.0E0;
  Double I_KINETIC_H4xz_Fy2z_a = 0.0E0;
  Double I_KINETIC_H3x2y_Fy2z_a = 0.0E0;
  Double I_KINETIC_H3xyz_Fy2z_a = 0.0E0;
  Double I_KINETIC_H3x2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_H2x3y_Fy2z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_H2x3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Hx4y_Fy2z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Hx4z_Fy2z_a = 0.0E0;
  Double I_KINETIC_H5y_Fy2z_a = 0.0E0;
  Double I_KINETIC_H4yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_H3y2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_H2y3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Hy4z_Fy2z_a = 0.0E0;
  Double I_KINETIC_H5z_Fy2z_a = 0.0E0;
  Double I_KINETIC_H5x_F3z_a = 0.0E0;
  Double I_KINETIC_H4xy_F3z_a = 0.0E0;
  Double I_KINETIC_H4xz_F3z_a = 0.0E0;
  Double I_KINETIC_H3x2y_F3z_a = 0.0E0;
  Double I_KINETIC_H3xyz_F3z_a = 0.0E0;
  Double I_KINETIC_H3x2z_F3z_a = 0.0E0;
  Double I_KINETIC_H2x3y_F3z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_F3z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_F3z_a = 0.0E0;
  Double I_KINETIC_H2x3z_F3z_a = 0.0E0;
  Double I_KINETIC_Hx4y_F3z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_F3z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_F3z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_F3z_a = 0.0E0;
  Double I_KINETIC_Hx4z_F3z_a = 0.0E0;
  Double I_KINETIC_H5y_F3z_a = 0.0E0;
  Double I_KINETIC_H4yz_F3z_a = 0.0E0;
  Double I_KINETIC_H3y2z_F3z_a = 0.0E0;
  Double I_KINETIC_H2y3z_F3z_a = 0.0E0;
  Double I_KINETIC_Hy4z_F3z_a = 0.0E0;
  Double I_KINETIC_H5z_F3z_a = 0.0E0;
  Double I_KINETIC_F3x_F3x = 0.0E0;
  Double I_KINETIC_F2xy_F3x = 0.0E0;
  Double I_KINETIC_F2xz_F3x = 0.0E0;
  Double I_KINETIC_Fx2y_F3x = 0.0E0;
  Double I_KINETIC_Fxyz_F3x = 0.0E0;
  Double I_KINETIC_Fx2z_F3x = 0.0E0;
  Double I_KINETIC_F3y_F3x = 0.0E0;
  Double I_KINETIC_F2yz_F3x = 0.0E0;
  Double I_KINETIC_Fy2z_F3x = 0.0E0;
  Double I_KINETIC_F3z_F3x = 0.0E0;
  Double I_KINETIC_F3x_F2xy = 0.0E0;
  Double I_KINETIC_F2xy_F2xy = 0.0E0;
  Double I_KINETIC_F2xz_F2xy = 0.0E0;
  Double I_KINETIC_Fx2y_F2xy = 0.0E0;
  Double I_KINETIC_Fxyz_F2xy = 0.0E0;
  Double I_KINETIC_Fx2z_F2xy = 0.0E0;
  Double I_KINETIC_F3y_F2xy = 0.0E0;
  Double I_KINETIC_F2yz_F2xy = 0.0E0;
  Double I_KINETIC_Fy2z_F2xy = 0.0E0;
  Double I_KINETIC_F3z_F2xy = 0.0E0;
  Double I_KINETIC_F3x_F2xz = 0.0E0;
  Double I_KINETIC_F2xy_F2xz = 0.0E0;
  Double I_KINETIC_F2xz_F2xz = 0.0E0;
  Double I_KINETIC_Fx2y_F2xz = 0.0E0;
  Double I_KINETIC_Fxyz_F2xz = 0.0E0;
  Double I_KINETIC_Fx2z_F2xz = 0.0E0;
  Double I_KINETIC_F3y_F2xz = 0.0E0;
  Double I_KINETIC_F2yz_F2xz = 0.0E0;
  Double I_KINETIC_Fy2z_F2xz = 0.0E0;
  Double I_KINETIC_F3z_F2xz = 0.0E0;
  Double I_KINETIC_F3x_Fx2y = 0.0E0;
  Double I_KINETIC_F2xy_Fx2y = 0.0E0;
  Double I_KINETIC_F2xz_Fx2y = 0.0E0;
  Double I_KINETIC_Fx2y_Fx2y = 0.0E0;
  Double I_KINETIC_Fxyz_Fx2y = 0.0E0;
  Double I_KINETIC_Fx2z_Fx2y = 0.0E0;
  Double I_KINETIC_F3y_Fx2y = 0.0E0;
  Double I_KINETIC_F2yz_Fx2y = 0.0E0;
  Double I_KINETIC_Fy2z_Fx2y = 0.0E0;
  Double I_KINETIC_F3z_Fx2y = 0.0E0;
  Double I_KINETIC_F3x_Fxyz = 0.0E0;
  Double I_KINETIC_F2xy_Fxyz = 0.0E0;
  Double I_KINETIC_F2xz_Fxyz = 0.0E0;
  Double I_KINETIC_Fx2y_Fxyz = 0.0E0;
  Double I_KINETIC_Fxyz_Fxyz = 0.0E0;
  Double I_KINETIC_Fx2z_Fxyz = 0.0E0;
  Double I_KINETIC_F3y_Fxyz = 0.0E0;
  Double I_KINETIC_F2yz_Fxyz = 0.0E0;
  Double I_KINETIC_Fy2z_Fxyz = 0.0E0;
  Double I_KINETIC_F3z_Fxyz = 0.0E0;
  Double I_KINETIC_F3x_Fx2z = 0.0E0;
  Double I_KINETIC_F2xy_Fx2z = 0.0E0;
  Double I_KINETIC_F2xz_Fx2z = 0.0E0;
  Double I_KINETIC_Fx2y_Fx2z = 0.0E0;
  Double I_KINETIC_Fxyz_Fx2z = 0.0E0;
  Double I_KINETIC_Fx2z_Fx2z = 0.0E0;
  Double I_KINETIC_F3y_Fx2z = 0.0E0;
  Double I_KINETIC_F2yz_Fx2z = 0.0E0;
  Double I_KINETIC_Fy2z_Fx2z = 0.0E0;
  Double I_KINETIC_F3z_Fx2z = 0.0E0;
  Double I_KINETIC_F3x_F3y = 0.0E0;
  Double I_KINETIC_F2xy_F3y = 0.0E0;
  Double I_KINETIC_F2xz_F3y = 0.0E0;
  Double I_KINETIC_Fx2y_F3y = 0.0E0;
  Double I_KINETIC_Fxyz_F3y = 0.0E0;
  Double I_KINETIC_Fx2z_F3y = 0.0E0;
  Double I_KINETIC_F3y_F3y = 0.0E0;
  Double I_KINETIC_F2yz_F3y = 0.0E0;
  Double I_KINETIC_Fy2z_F3y = 0.0E0;
  Double I_KINETIC_F3z_F3y = 0.0E0;
  Double I_KINETIC_F3x_F2yz = 0.0E0;
  Double I_KINETIC_F2xy_F2yz = 0.0E0;
  Double I_KINETIC_F2xz_F2yz = 0.0E0;
  Double I_KINETIC_Fx2y_F2yz = 0.0E0;
  Double I_KINETIC_Fxyz_F2yz = 0.0E0;
  Double I_KINETIC_Fx2z_F2yz = 0.0E0;
  Double I_KINETIC_F3y_F2yz = 0.0E0;
  Double I_KINETIC_F2yz_F2yz = 0.0E0;
  Double I_KINETIC_Fy2z_F2yz = 0.0E0;
  Double I_KINETIC_F3z_F2yz = 0.0E0;
  Double I_KINETIC_F3x_Fy2z = 0.0E0;
  Double I_KINETIC_F2xy_Fy2z = 0.0E0;
  Double I_KINETIC_F2xz_Fy2z = 0.0E0;
  Double I_KINETIC_Fx2y_Fy2z = 0.0E0;
  Double I_KINETIC_Fxyz_Fy2z = 0.0E0;
  Double I_KINETIC_Fx2z_Fy2z = 0.0E0;
  Double I_KINETIC_F3y_Fy2z = 0.0E0;
  Double I_KINETIC_F2yz_Fy2z = 0.0E0;
  Double I_KINETIC_Fy2z_Fy2z = 0.0E0;
  Double I_KINETIC_F3z_Fy2z = 0.0E0;
  Double I_KINETIC_F3x_F3z = 0.0E0;
  Double I_KINETIC_F2xy_F3z = 0.0E0;
  Double I_KINETIC_F2xz_F3z = 0.0E0;
  Double I_KINETIC_Fx2y_F3z = 0.0E0;
  Double I_KINETIC_Fxyz_F3z = 0.0E0;
  Double I_KINETIC_Fx2z_F3z = 0.0E0;
  Double I_KINETIC_F3y_F3z = 0.0E0;
  Double I_KINETIC_F2yz_F3z = 0.0E0;
  Double I_KINETIC_Fy2z_F3z = 0.0E0;
  Double I_KINETIC_F3z_F3z = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_S_vrr = PAY*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_S_vrr = PAZ*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_S_vrr = PAZ*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_S_vrr = PAZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_D2x_vrr = PBX*I_TWOBODYOVERLAP_Px_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Py_D2x_vrr = PBX*I_TWOBODYOVERLAP_Py_Px_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_D2x_vrr = PBX*I_TWOBODYOVERLAP_Pz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_Px_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_Py_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Py_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_Pz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_Px_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_Py_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_Pz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Pz_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_Px_D2y_vrr = PBY*I_TWOBODYOVERLAP_Px_Py_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Py_D2y_vrr = PBY*I_TWOBODYOVERLAP_Py_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_D2y_vrr = PBY*I_TWOBODYOVERLAP_Pz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_Px_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_Py_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_TWOBODYOVERLAP_Pz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Pz_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Px_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Px_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Py_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Py_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Pz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_Px_vrr = PAX*I_TWOBODYOVERLAP_Px_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Px_vrr = PAY*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Px_vrr = PAZ*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_Px_vrr = PAY*I_TWOBODYOVERLAP_Py_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Px_vrr = PAZ*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_Py_vrr = PAX*I_TWOBODYOVERLAP_Px_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Py_vrr = PAY*I_TWOBODYOVERLAP_Px_Py_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Py_vrr = PAZ*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_Py_vrr = PAY*I_TWOBODYOVERLAP_Py_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Py_vrr = PAZ*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_TWOBODYOVERLAP_D2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_Pz_vrr = PAX*I_TWOBODYOVERLAP_Px_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Pz_vrr = PAY*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Px_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_Py_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Py_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_D2x_vrr = PBX*I_TWOBODYOVERLAP_D2x_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2x_vrr = PBX*I_TWOBODYOVERLAP_Dxy_Px_vrr+oned2z*I_TWOBODYOVERLAP_Py_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_D2x_vrr = PBX*I_TWOBODYOVERLAP_Dxz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_D2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_D2x_vrr = PBX*I_TWOBODYOVERLAP_Dyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_D2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxy_vrr = PBY*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Dxy_Px_vrr+oned2z*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_D2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Dyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2y_vrr = PBY*I_TWOBODYOVERLAP_D2x_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2y_vrr = PBY*I_TWOBODYOVERLAP_Dxy_Py_vrr+oned2z*I_TWOBODYOVERLAP_Px_Py_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_D2y_vrr = PBY*I_TWOBODYOVERLAP_Dxz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_D2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_D2y_vrr = PBY*I_TWOBODYOVERLAP_Dyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_D2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Dxz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_D2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Dyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Py_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_D2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_Px_vrr = PAX*I_TWOBODYOVERLAP_D2x_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Px_vrr = PAY*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Px_vrr = PAX*I_TWOBODYOVERLAP_D2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Px_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Px_vrr = PAX*I_TWOBODYOVERLAP_D2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Px_vrr = PAY*I_TWOBODYOVERLAP_D2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Px_vrr = PAY*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Px_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_Py_vrr = PAX*I_TWOBODYOVERLAP_D2x_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Py_vrr = PAY*I_TWOBODYOVERLAP_D2x_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Py_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Py_vrr = PAX*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Py_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Py_vrr = PAX*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_Py_vrr = PAY*I_TWOBODYOVERLAP_D2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Py_vrr = PAY*I_TWOBODYOVERLAP_D2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Py_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_TWOBODYOVERLAP_F3x_Pz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Pz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Pz_vrr = PAX*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Pz_vrr = PAX*I_TWOBODYOVERLAP_D2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3y_Pz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Pz_vrr = PAY*I_TWOBODYOVERLAP_D2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_F3x_vrr = PBX*I_TWOBODYOVERLAP_D2x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3x_vrr = PBX*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2xy_vrr = PBY*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2x_F3y_vrr = PBY*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3y_vrr = PBY*I_TWOBODYOVERLAP_D2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2x_F3z_vrr = PBZ*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_D2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 8 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
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
    Double I_TWOBODYOVERLAP_Fx2y_D2y_vrr = PAX*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2y_vrr = PAX*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2y_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2y_vrr = PAY*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
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
    Double I_TWOBODYOVERLAP_Fx2y_D2z_vrr = PAX*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2z_vrr = PAX*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2z_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2z_vrr = PAY*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_F3x_vrr = PBX*I_TWOBODYOVERLAP_F3x_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F3x_vrr = PBX*I_TWOBODYOVERLAP_F2xy_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F3x_vrr = PBX*I_TWOBODYOVERLAP_F2xz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F3x_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_D2x_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F3x_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_F3x_vrr = PBX*I_TWOBODYOVERLAP_F3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_F2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_F3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_F2xy_vrr = PBY*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F2xy_vrr = PBY*I_TWOBODYOVERLAP_F2xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_F3y_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_F2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_F3x_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_F2xy_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_F2xz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_D2y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_F3x_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_F2xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_F2xz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_D2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3x_F3y_vrr = PBY*I_TWOBODYOVERLAP_F3x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F3y_vrr = PBY*I_TWOBODYOVERLAP_F2xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F3y_vrr = PBY*I_TWOBODYOVERLAP_F2xz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F3y_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F3y_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_F3y_vrr = PBY*I_TWOBODYOVERLAP_F3y_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_F2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_F3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3x_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_F2xy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_F3y_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_F2yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3x_F3z_vrr = PBZ*I_TWOBODYOVERLAP_F3x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F3z_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_F3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_F3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;

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
     * shell quartet name: SQ_KINETIC_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_KINETIC_D2x_S_vrr = PAX*I_KINETIC_Px_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_Dxy_S_vrr = PAY*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_KINETIC_Dxz_S_vrr = PAZ*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_KINETIC_D2y_S_vrr = PAY*I_KINETIC_Py_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_Dyz_S_vrr = PAZ*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_KINETIC_D2z_S_vrr = PAZ*I_KINETIC_Pz_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_Px_D2x_vrr = PBX*I_KINETIC_Px_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_D2x_vrr-adz*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_Py_D2x_vrr = PBX*I_KINETIC_Py_Px_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_D2x_vrr-adz*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_Pz_D2x_vrr = PBX*I_KINETIC_Pz_Px_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_D2x_vrr-adz*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_KINETIC_Px_Dxy_vrr = PBY*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_KINETIC_Py_Dxy_vrr = PBY*I_KINETIC_Py_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_KINETIC_Pz_Dxy_vrr = PBY*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_KINETIC_Px_Dxz_vrr = PBZ*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_KINETIC_Py_Dxz_vrr = PBZ*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_KINETIC_Pz_Dxz_vrr = PBZ*I_KINETIC_Pz_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Dxz_vrr;
    Double I_KINETIC_Px_D2y_vrr = PBY*I_KINETIC_Px_Py_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_D2y_vrr-adz*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_Py_D2y_vrr = PBY*I_KINETIC_Py_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_D2y_vrr-adz*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_Pz_D2y_vrr = PBY*I_KINETIC_Pz_Py_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_D2y_vrr-adz*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_KINETIC_Px_Dyz_vrr = PBZ*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_KINETIC_Py_Dyz_vrr = PBZ*I_KINETIC_Py_Py_vrr+twoxi*I_TWOBODYOVERLAP_Py_Dyz_vrr;
    Double I_KINETIC_Pz_Dyz_vrr = PBZ*I_KINETIC_Pz_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Dyz_vrr;
    Double I_KINETIC_Px_D2z_vrr = PBZ*I_KINETIC_Px_Pz_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_D2z_vrr-adz*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_Py_D2z_vrr = PBZ*I_KINETIC_Py_Pz_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_D2z_vrr-adz*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_Pz_D2z_vrr = PBZ*I_KINETIC_Pz_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_D2z_vrr-adz*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_KINETIC_D2x_Px_vrr = PAX*I_KINETIC_Px_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_Dxy_Px_vrr = PAY*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_KINETIC_Dxz_Px_vrr = PAZ*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_KINETIC_D2y_Px_vrr = PAY*I_KINETIC_Py_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_Dyz_Px_vrr = PAZ*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_KINETIC_D2z_Px_vrr = PAZ*I_KINETIC_Pz_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_D2x_Py_vrr = PAX*I_KINETIC_Px_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_Dxy_Py_vrr = PAY*I_KINETIC_Px_Py_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_KINETIC_Dxz_Py_vrr = PAZ*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Py_vrr;
    Double I_KINETIC_D2y_Py_vrr = PAY*I_KINETIC_Py_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_Dyz_Py_vrr = PAZ*I_KINETIC_Py_Py_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Py_vrr;
    Double I_KINETIC_D2z_Py_vrr = PAZ*I_KINETIC_Pz_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_D2x_Pz_vrr = PAX*I_KINETIC_Px_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_KINETIC_Dxy_Pz_vrr = PAY*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_KINETIC_Dxz_Pz_vrr = PAZ*I_KINETIC_Px_Pz_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Pz_vrr;
    Double I_KINETIC_D2y_Pz_vrr = PAY*I_KINETIC_Py_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_KINETIC_Dyz_Pz_vrr = PAZ*I_KINETIC_Py_Pz_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Pz_vrr;
    Double I_KINETIC_D2z_Pz_vrr = PAZ*I_KINETIC_Pz_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_KINETIC_D2x_D2x_vrr = PBX*I_KINETIC_D2x_Px_vrr+2*oned2z*I_KINETIC_Px_Px_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2x_vrr-adz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_Dxy_D2x_vrr = PBX*I_KINETIC_Dxy_Px_vrr+oned2z*I_KINETIC_Py_Px_vrr+oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2x_vrr-adz*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_KINETIC_Dxz_D2x_vrr = PBX*I_KINETIC_Dxz_Px_vrr+oned2z*I_KINETIC_Pz_Px_vrr+oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_D2x_vrr-adz*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_KINETIC_D2y_D2x_vrr = PBX*I_KINETIC_D2y_Px_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2x_vrr-adz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_Dyz_D2x_vrr = PBX*I_KINETIC_Dyz_Px_vrr+oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_D2x_vrr-adz*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_KINETIC_D2z_D2x_vrr = PBX*I_KINETIC_D2z_Px_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2x_vrr-adz*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_KINETIC_D2x_Dxy_vrr = PBY*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_KINETIC_Dxy_Dxy_vrr = PBY*I_KINETIC_Dxy_Px_vrr+oned2z*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Dxy_vrr;
    Double I_KINETIC_Dxz_Dxy_vrr = PBY*I_KINETIC_Dxz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Dxy_vrr;
    Double I_KINETIC_D2y_Dxy_vrr = PBY*I_KINETIC_D2y_Px_vrr+2*oned2z*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_KINETIC_Dyz_Dxy_vrr = PBY*I_KINETIC_Dyz_Px_vrr+oned2z*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Dxy_vrr;
    Double I_KINETIC_D2z_Dxy_vrr = PBY*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxy_vrr;
    Double I_KINETIC_D2x_Dxz_vrr = PBZ*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_KINETIC_D2y_Dxz_vrr = PBZ*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_KINETIC_D2z_Dxz_vrr = PBZ*I_KINETIC_D2z_Px_vrr+2*oned2z*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxz_vrr;
    Double I_KINETIC_D2x_D2y_vrr = PBY*I_KINETIC_D2x_Py_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2y_vrr-adz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_Dxy_D2y_vrr = PBY*I_KINETIC_Dxy_Py_vrr+oned2z*I_KINETIC_Px_Py_vrr+oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2y_vrr-adz*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_KINETIC_Dxz_D2y_vrr = PBY*I_KINETIC_Dxz_Py_vrr+oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_D2y_vrr-adz*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_KINETIC_D2y_D2y_vrr = PBY*I_KINETIC_D2y_Py_vrr+2*oned2z*I_KINETIC_Py_Py_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2y_vrr-adz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_Dyz_D2y_vrr = PBY*I_KINETIC_Dyz_Py_vrr+oned2z*I_KINETIC_Pz_Py_vrr+oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_D2y_vrr-adz*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_KINETIC_D2z_D2y_vrr = PBY*I_KINETIC_D2z_Py_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2y_vrr-adz*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_KINETIC_D2x_Dyz_vrr = PBZ*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_KINETIC_D2y_Dyz_vrr = PBZ*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_KINETIC_D2z_Dyz_vrr = PBZ*I_KINETIC_D2z_Py_vrr+2*oned2z*I_KINETIC_Pz_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dyz_vrr;
    Double I_KINETIC_D2x_D2z_vrr = PBZ*I_KINETIC_D2x_Pz_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2z_vrr-adz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_Dxy_D2z_vrr = PBZ*I_KINETIC_Dxy_Pz_vrr+oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2z_vrr-adz*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_KINETIC_Dxz_D2z_vrr = PBZ*I_KINETIC_Dxz_Pz_vrr+oned2z*I_KINETIC_Px_Pz_vrr+oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_D2z_vrr-adz*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_KINETIC_D2y_D2z_vrr = PBZ*I_KINETIC_D2y_Pz_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2z_vrr-adz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_Dyz_D2z_vrr = PBZ*I_KINETIC_Dyz_Pz_vrr+oned2z*I_KINETIC_Py_Pz_vrr+oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_D2z_vrr-adz*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_KINETIC_D2z_D2z_vrr = PBZ*I_KINETIC_D2z_Pz_vrr+2*oned2z*I_KINETIC_Pz_Pz_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2z_vrr-adz*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_KINETIC_F3x_Px_vrr = PAX*I_KINETIC_D2x_Px_vrr+2*oned2z*I_KINETIC_Px_Px_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Px_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_KINETIC_F2xy_Px_vrr = PAY*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_KINETIC_F2xz_Px_vrr = PAZ*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_KINETIC_Fx2y_Px_vrr = PAX*I_KINETIC_D2y_Px_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_KINETIC_Fxyz_Px_vrr = PAZ*I_KINETIC_Dxy_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Px_vrr;
    Double I_KINETIC_Fx2z_Px_vrr = PAX*I_KINETIC_D2z_Px_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_KINETIC_F3y_Px_vrr = PAY*I_KINETIC_D2y_Px_vrr+2*oned2z*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Px_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_KINETIC_F2yz_Px_vrr = PAZ*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_KINETIC_Fy2z_Px_vrr = PAY*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_KINETIC_F3z_Px_vrr = PAZ*I_KINETIC_D2z_Px_vrr+2*oned2z*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Px_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_KINETIC_F3x_Py_vrr = PAX*I_KINETIC_D2x_Py_vrr+2*oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Py_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_KINETIC_F2xy_Py_vrr = PAY*I_KINETIC_D2x_Py_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_KINETIC_F2xz_Py_vrr = PAZ*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_KINETIC_Fx2y_Py_vrr = PAX*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_KINETIC_Fxyz_Py_vrr = PAZ*I_KINETIC_Dxy_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Py_vrr;
    Double I_KINETIC_Fx2z_Py_vrr = PAX*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Py_vrr;
    Double I_KINETIC_F3y_Py_vrr = PAY*I_KINETIC_D2y_Py_vrr+2*oned2z*I_KINETIC_Py_Py_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Py_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_KINETIC_F2yz_Py_vrr = PAZ*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Py_vrr;
    Double I_KINETIC_Fy2z_Py_vrr = PAY*I_KINETIC_D2z_Py_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Py_vrr;
    Double I_KINETIC_F3z_Py_vrr = PAZ*I_KINETIC_D2z_Py_vrr+2*oned2z*I_KINETIC_Pz_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Py_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_KINETIC_F3x_Pz_vrr = PAX*I_KINETIC_D2x_Pz_vrr+2*oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Pz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_KINETIC_F2xy_Pz_vrr = PAY*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Pz_vrr;
    Double I_KINETIC_F2xz_Pz_vrr = PAZ*I_KINETIC_D2x_Pz_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Pz_vrr;
    Double I_KINETIC_Fx2y_Pz_vrr = PAX*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Pz_vrr;
    Double I_KINETIC_Fxyz_Pz_vrr = PAZ*I_KINETIC_Dxy_Pz_vrr+oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Pz_vrr;
    Double I_KINETIC_Fx2z_Pz_vrr = PAX*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Pz_vrr;
    Double I_KINETIC_F3y_Pz_vrr = PAY*I_KINETIC_D2y_Pz_vrr+2*oned2z*I_KINETIC_Py_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Pz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_KINETIC_F2yz_Pz_vrr = PAZ*I_KINETIC_D2y_Pz_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Pz_vrr;
    Double I_KINETIC_Fy2z_Pz_vrr = PAY*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Pz_vrr;
    Double I_KINETIC_F3z_Pz_vrr = PAZ*I_KINETIC_D2z_Pz_vrr+2*oned2z*I_KINETIC_Pz_Pz_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Pz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_KINETIC_D2x_F3x_vrr = PBX*I_KINETIC_D2x_D2x_vrr+2*oned2z*I_KINETIC_Px_D2x_vrr+2*oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3x_vrr-2*adz*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_D2y_F3x_vrr = PBX*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_D2z_F3x_vrr = PBX*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_KINETIC_D2x_F2xy_vrr = PBY*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_KINETIC_D2y_F2xy_vrr = PBY*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_KINETIC_D2z_F2xy_vrr = PBY*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2xy_vrr;
    Double I_KINETIC_D2x_F2xz_vrr = PBZ*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_KINETIC_D2y_F2xz_vrr = PBZ*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2xz_vrr;
    Double I_KINETIC_D2z_F2xz_vrr = PBZ*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2xz_vrr;
    Double I_KINETIC_D2x_Fx2y_vrr = PBX*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_KINETIC_D2y_Fx2y_vrr = PBX*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_KINETIC_D2z_Fx2y_vrr = PBX*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fx2y_vrr;
    Double I_KINETIC_D2x_Fxyz_vrr = PBZ*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fxyz_vrr;
    Double I_KINETIC_D2y_Fxyz_vrr = PBZ*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fxyz_vrr;
    Double I_KINETIC_D2z_Fxyz_vrr = PBZ*I_KINETIC_D2z_Dxy_vrr+2*oned2z*I_KINETIC_Pz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fxyz_vrr;
    Double I_KINETIC_D2x_Fx2z_vrr = PBX*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_KINETIC_D2y_Fx2z_vrr = PBX*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_KINETIC_D2z_Fx2z_vrr = PBX*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fx2z_vrr;
    Double I_KINETIC_D2x_F3y_vrr = PBY*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3y_vrr-2*adz*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_D2y_F3y_vrr = PBY*I_KINETIC_D2y_D2y_vrr+2*oned2z*I_KINETIC_Py_D2y_vrr+2*oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_D2z_F3y_vrr = PBY*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_KINETIC_D2x_F2yz_vrr = PBZ*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_KINETIC_D2y_F2yz_vrr = PBZ*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2yz_vrr;
    Double I_KINETIC_D2z_F2yz_vrr = PBZ*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_Pz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2yz_vrr;
    Double I_KINETIC_D2x_Fy2z_vrr = PBY*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_KINETIC_D2y_Fy2z_vrr = PBY*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fy2z_vrr;
    Double I_KINETIC_D2z_Fy2z_vrr = PBY*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fy2z_vrr;
    Double I_KINETIC_D2x_F3z_vrr = PBZ*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3z_vrr-2*adz*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_D2y_F3z_vrr = PBZ*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_D2z_F3z_vrr = PBZ*I_KINETIC_D2z_D2z_vrr+2*oned2z*I_KINETIC_Pz_D2z_vrr+2*oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 8 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     ************************************************************/
    Double I_KINETIC_F3x_D2x_vrr = PAX*I_KINETIC_D2x_D2x_vrr+2*oned2z*I_KINETIC_Px_D2x_vrr+2*oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_KINETIC_F2xy_D2x_vrr = PAY*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_KINETIC_F2xz_D2x_vrr = PAZ*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_KINETIC_Fx2y_D2x_vrr = PAX*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_KINETIC_Fxyz_D2x_vrr = PAZ*I_KINETIC_Dxy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_KINETIC_Fx2z_D2x_vrr = PAX*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_KINETIC_F3y_D2x_vrr = PAY*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_KINETIC_F2yz_D2x_vrr = PAZ*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_KINETIC_Fy2z_D2x_vrr = PAY*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_KINETIC_F3z_D2x_vrr = PAZ*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_KINETIC_F3x_Dxy_vrr = PAX*I_KINETIC_D2x_Dxy_vrr+2*oned2z*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_KINETIC_F2xy_Dxy_vrr = PAY*I_KINETIC_D2x_Dxy_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_KINETIC_F2xz_Dxy_vrr = PAZ*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_KINETIC_Fx2y_Dxy_vrr = PAX*I_KINETIC_D2y_Dxy_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_KINETIC_Fxyz_Dxy_vrr = PAZ*I_KINETIC_Dxy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dxy_vrr;
    Double I_KINETIC_Fx2z_Dxy_vrr = PAX*I_KINETIC_D2z_Dxy_vrr+oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr;
    Double I_KINETIC_F3y_Dxy_vrr = PAY*I_KINETIC_D2y_Dxy_vrr+2*oned2z*I_KINETIC_Py_Dxy_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_KINETIC_F2yz_Dxy_vrr = PAZ*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_KINETIC_Fy2z_Dxy_vrr = PAY*I_KINETIC_D2z_Dxy_vrr+oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dxy_vrr;
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
    Double I_KINETIC_Fx2y_D2y_vrr = PAX*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_KINETIC_Fxyz_D2y_vrr = PAZ*I_KINETIC_Dxy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_KINETIC_Fx2z_D2y_vrr = PAX*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_KINETIC_F3y_D2y_vrr = PAY*I_KINETIC_D2y_D2y_vrr+2*oned2z*I_KINETIC_Py_D2y_vrr+2*oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_KINETIC_F2yz_D2y_vrr = PAZ*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_KINETIC_Fy2z_D2y_vrr = PAY*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
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
    Double I_KINETIC_Fx2y_D2z_vrr = PAX*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_KINETIC_Fxyz_D2z_vrr = PAZ*I_KINETIC_Dxy_D2z_vrr+2*oned2z*I_KINETIC_Dxy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2z_vrr;
    Double I_KINETIC_Fx2z_D2z_vrr = PAX*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_KINETIC_F3y_D2z_vrr = PAY*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_KINETIC_F2yz_D2z_vrr = PAZ*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_KINETIC_Fy2z_D2z_vrr = PAY*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2z_vrr;
    Double I_KINETIC_F3z_D2z_vrr = PAZ*I_KINETIC_D2z_D2z_vrr+2*oned2z*I_KINETIC_Pz_D2z_vrr+2*oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_D
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     ************************************************************/
    Double I_KINETIC_F3x_F3x_vrr = PBX*I_KINETIC_F3x_D2x_vrr+3*oned2z*I_KINETIC_D2x_D2x_vrr+2*oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F3x_vrr-2*adz*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_KINETIC_F2xy_F3x_vrr = PBX*I_KINETIC_F2xy_D2x_vrr+2*oned2z*I_KINETIC_Dxy_D2x_vrr+2*oned2z*I_KINETIC_F2xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F3x_vrr-2*adz*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_KINETIC_F2xz_F3x_vrr = PBX*I_KINETIC_F2xz_D2x_vrr+2*oned2z*I_KINETIC_Dxz_D2x_vrr+2*oned2z*I_KINETIC_F2xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_KINETIC_Fx2y_F3x_vrr = PBX*I_KINETIC_Fx2y_D2x_vrr+oned2z*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_Fx2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_KINETIC_Fxyz_F3x_vrr = PBX*I_KINETIC_Fxyz_D2x_vrr+oned2z*I_KINETIC_Dyz_D2x_vrr+2*oned2z*I_KINETIC_Fxyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Fxyz_Px_vrr;
    Double I_KINETIC_Fx2z_F3x_vrr = PBX*I_KINETIC_Fx2z_D2x_vrr+oned2z*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_KINETIC_F3y_F3x_vrr = PBX*I_KINETIC_F3y_D2x_vrr+2*oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_KINETIC_F2yz_F3x_vrr = PBX*I_KINETIC_F2yz_D2x_vrr+2*oned2z*I_KINETIC_F2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_KINETIC_Fy2z_F3x_vrr = PBX*I_KINETIC_Fy2z_D2x_vrr+2*oned2z*I_KINETIC_Fy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_KINETIC_F3z_F3x_vrr = PBX*I_KINETIC_F3z_D2x_vrr+2*oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_KINETIC_F3x_F2xy_vrr = PBY*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_KINETIC_F2xy_F2xy_vrr = PBY*I_KINETIC_F2xy_D2x_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F2xy_vrr;
    Double I_KINETIC_F2xz_F2xy_vrr = PBY*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F2xy_vrr;
    Double I_KINETIC_Fx2y_F2xy_vrr = PBY*I_KINETIC_Fx2y_D2x_vrr+2*oned2z*I_KINETIC_Dxy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F2xy_vrr;
    Double I_KINETIC_Fxyz_F2xy_vrr = PBY*I_KINETIC_Fxyz_D2x_vrr+oned2z*I_KINETIC_Dxz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F2xy_vrr;
    Double I_KINETIC_Fx2z_F2xy_vrr = PBY*I_KINETIC_Fx2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F2xy_vrr;
    Double I_KINETIC_F3y_F2xy_vrr = PBY*I_KINETIC_F3y_D2x_vrr+3*oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_KINETIC_F2yz_F2xy_vrr = PBY*I_KINETIC_F2yz_D2x_vrr+2*oned2z*I_KINETIC_Dyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F2xy_vrr;
    Double I_KINETIC_Fy2z_F2xy_vrr = PBY*I_KINETIC_Fy2z_D2x_vrr+oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F2xy_vrr;
    Double I_KINETIC_F3z_F2xy_vrr = PBY*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F2xy_vrr;
    Double I_KINETIC_F3x_F2xz_vrr = PBZ*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_KINETIC_F2xy_F2xz_vrr = PBZ*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F2xz_vrr;
    Double I_KINETIC_F2xz_F2xz_vrr = PBZ*I_KINETIC_F2xz_D2x_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F2xz_vrr;
    Double I_KINETIC_Fx2y_F2xz_vrr = PBZ*I_KINETIC_Fx2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F2xz_vrr;
    Double I_KINETIC_Fxyz_F2xz_vrr = PBZ*I_KINETIC_Fxyz_D2x_vrr+oned2z*I_KINETIC_Dxy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F2xz_vrr;
    Double I_KINETIC_Fx2z_F2xz_vrr = PBZ*I_KINETIC_Fx2z_D2x_vrr+2*oned2z*I_KINETIC_Dxz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F2xz_vrr;
    Double I_KINETIC_F3y_F2xz_vrr = PBZ*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F2xz_vrr;
    Double I_KINETIC_F2yz_F2xz_vrr = PBZ*I_KINETIC_F2yz_D2x_vrr+oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F2xz_vrr;
    Double I_KINETIC_Fy2z_F2xz_vrr = PBZ*I_KINETIC_Fy2z_D2x_vrr+2*oned2z*I_KINETIC_Dyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F2xz_vrr;
    Double I_KINETIC_F3z_F2xz_vrr = PBZ*I_KINETIC_F3z_D2x_vrr+3*oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_KINETIC_F3x_Fx2y_vrr = PBX*I_KINETIC_F3x_D2y_vrr+3*oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_KINETIC_F2xy_Fx2y_vrr = PBX*I_KINETIC_F2xy_D2y_vrr+2*oned2z*I_KINETIC_Dxy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fx2y_vrr;
    Double I_KINETIC_F2xz_Fx2y_vrr = PBX*I_KINETIC_F2xz_D2y_vrr+2*oned2z*I_KINETIC_Dxz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fx2y_vrr;
    Double I_KINETIC_Fx2y_Fx2y_vrr = PBX*I_KINETIC_Fx2y_D2y_vrr+oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fx2y_vrr;
    Double I_KINETIC_Fxyz_Fx2y_vrr = PBX*I_KINETIC_Fxyz_D2y_vrr+oned2z*I_KINETIC_Dyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Fx2y_vrr;
    Double I_KINETIC_Fx2z_Fx2y_vrr = PBX*I_KINETIC_Fx2z_D2y_vrr+oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fx2y_vrr;
    Double I_KINETIC_F3y_Fx2y_vrr = PBX*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_KINETIC_F2yz_Fx2y_vrr = PBX*I_KINETIC_F2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fx2y_vrr;
    Double I_KINETIC_Fy2z_Fx2y_vrr = PBX*I_KINETIC_Fy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Fx2y_vrr;
    Double I_KINETIC_F3z_Fx2y_vrr = PBX*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fx2y_vrr;
    Double I_KINETIC_F3x_Fxyz_vrr = PBZ*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fxyz_vrr;
    Double I_KINETIC_F2xy_Fxyz_vrr = PBZ*I_KINETIC_F2xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fxyz_vrr;
    Double I_KINETIC_F2xz_Fxyz_vrr = PBZ*I_KINETIC_F2xz_Dxy_vrr+oned2z*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fxyz_vrr;
    Double I_KINETIC_Fx2y_Fxyz_vrr = PBZ*I_KINETIC_Fx2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fxyz_vrr;
    Double I_KINETIC_Fxyz_Fxyz_vrr = PBZ*I_KINETIC_Fxyz_Dxy_vrr+oned2z*I_KINETIC_Dxy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Fxyz_vrr;
    Double I_KINETIC_Fx2z_Fxyz_vrr = PBZ*I_KINETIC_Fx2z_Dxy_vrr+2*oned2z*I_KINETIC_Dxz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fxyz_vrr;
    Double I_KINETIC_F3y_Fxyz_vrr = PBZ*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fxyz_vrr;
    Double I_KINETIC_F2yz_Fxyz_vrr = PBZ*I_KINETIC_F2yz_Dxy_vrr+oned2z*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fxyz_vrr;
    Double I_KINETIC_Fy2z_Fxyz_vrr = PBZ*I_KINETIC_Fy2z_Dxy_vrr+2*oned2z*I_KINETIC_Dyz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Fxyz_vrr;
    Double I_KINETIC_F3z_Fxyz_vrr = PBZ*I_KINETIC_F3z_Dxy_vrr+3*oned2z*I_KINETIC_D2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fxyz_vrr;
    Double I_KINETIC_F3x_Fx2z_vrr = PBX*I_KINETIC_F3x_D2z_vrr+3*oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_KINETIC_F2xy_Fx2z_vrr = PBX*I_KINETIC_F2xy_D2z_vrr+2*oned2z*I_KINETIC_Dxy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fx2z_vrr;
    Double I_KINETIC_F2xz_Fx2z_vrr = PBX*I_KINETIC_F2xz_D2z_vrr+2*oned2z*I_KINETIC_Dxz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fx2z_vrr;
    Double I_KINETIC_Fx2y_Fx2z_vrr = PBX*I_KINETIC_Fx2y_D2z_vrr+oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fx2z_vrr;
    Double I_KINETIC_Fxyz_Fx2z_vrr = PBX*I_KINETIC_Fxyz_D2z_vrr+oned2z*I_KINETIC_Dyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Fx2z_vrr;
    Double I_KINETIC_Fx2z_Fx2z_vrr = PBX*I_KINETIC_Fx2z_D2z_vrr+oned2z*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fx2z_vrr;
    Double I_KINETIC_F3y_Fx2z_vrr = PBX*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fx2z_vrr;
    Double I_KINETIC_F2yz_Fx2z_vrr = PBX*I_KINETIC_F2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fx2z_vrr;
    Double I_KINETIC_Fy2z_Fx2z_vrr = PBX*I_KINETIC_Fy2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Fx2z_vrr;
    Double I_KINETIC_F3z_Fx2z_vrr = PBX*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_KINETIC_F3x_F3y_vrr = PBY*I_KINETIC_F3x_D2y_vrr+2*oned2z*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F3y_vrr-2*adz*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_KINETIC_F2xy_F3y_vrr = PBY*I_KINETIC_F2xy_D2y_vrr+oned2z*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_F2xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F3y_vrr-2*adz*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_KINETIC_F2xz_F3y_vrr = PBY*I_KINETIC_F2xz_D2y_vrr+2*oned2z*I_KINETIC_F2xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_KINETIC_Fx2y_F3y_vrr = PBY*I_KINETIC_Fx2y_D2y_vrr+2*oned2z*I_KINETIC_Dxy_D2y_vrr+2*oned2z*I_KINETIC_Fx2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_KINETIC_Fxyz_F3y_vrr = PBY*I_KINETIC_Fxyz_D2y_vrr+oned2z*I_KINETIC_Dxz_D2y_vrr+2*oned2z*I_KINETIC_Fxyz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Fxyz_Py_vrr;
    Double I_KINETIC_Fx2z_F3y_vrr = PBY*I_KINETIC_Fx2z_D2y_vrr+2*oned2z*I_KINETIC_Fx2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Fx2z_Py_vrr;
    Double I_KINETIC_F3y_F3y_vrr = PBY*I_KINETIC_F3y_D2y_vrr+3*oned2z*I_KINETIC_D2y_D2y_vrr+2*oned2z*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_KINETIC_F2yz_F3y_vrr = PBY*I_KINETIC_F2yz_D2y_vrr+2*oned2z*I_KINETIC_Dyz_D2y_vrr+2*oned2z*I_KINETIC_F2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_F2yz_Py_vrr;
    Double I_KINETIC_Fy2z_F3y_vrr = PBY*I_KINETIC_Fy2z_D2y_vrr+oned2z*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_Fy2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Fy2z_Py_vrr;
    Double I_KINETIC_F3z_F3y_vrr = PBY*I_KINETIC_F3z_D2y_vrr+2*oned2z*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_KINETIC_F3x_F2yz_vrr = PBZ*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F2yz_vrr;
    Double I_KINETIC_F2xy_F2yz_vrr = PBZ*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F2yz_vrr;
    Double I_KINETIC_F2xz_F2yz_vrr = PBZ*I_KINETIC_F2xz_D2y_vrr+oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F2yz_vrr;
    Double I_KINETIC_Fx2y_F2yz_vrr = PBZ*I_KINETIC_Fx2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F2yz_vrr;
    Double I_KINETIC_Fxyz_F2yz_vrr = PBZ*I_KINETIC_Fxyz_D2y_vrr+oned2z*I_KINETIC_Dxy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F2yz_vrr;
    Double I_KINETIC_Fx2z_F2yz_vrr = PBZ*I_KINETIC_Fx2z_D2y_vrr+2*oned2z*I_KINETIC_Dxz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F2yz_vrr;
    Double I_KINETIC_F3y_F2yz_vrr = PBZ*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_KINETIC_F2yz_F2yz_vrr = PBZ*I_KINETIC_F2yz_D2y_vrr+oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F2yz_vrr;
    Double I_KINETIC_Fy2z_F2yz_vrr = PBZ*I_KINETIC_Fy2z_D2y_vrr+2*oned2z*I_KINETIC_Dyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F2yz_vrr;
    Double I_KINETIC_F3z_F2yz_vrr = PBZ*I_KINETIC_F3z_D2y_vrr+3*oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_KINETIC_F3x_Fy2z_vrr = PBY*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fy2z_vrr;
    Double I_KINETIC_F2xy_Fy2z_vrr = PBY*I_KINETIC_F2xy_D2z_vrr+oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fy2z_vrr;
    Double I_KINETIC_F2xz_Fy2z_vrr = PBY*I_KINETIC_F2xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fy2z_vrr;
    Double I_KINETIC_Fx2y_Fy2z_vrr = PBY*I_KINETIC_Fx2y_D2z_vrr+2*oned2z*I_KINETIC_Dxy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fy2z_vrr;
    Double I_KINETIC_Fxyz_Fy2z_vrr = PBY*I_KINETIC_Fxyz_D2z_vrr+oned2z*I_KINETIC_Dxz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Fy2z_vrr;
    Double I_KINETIC_Fx2z_Fy2z_vrr = PBY*I_KINETIC_Fx2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fy2z_vrr;
    Double I_KINETIC_F3y_Fy2z_vrr = PBY*I_KINETIC_F3y_D2z_vrr+3*oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_KINETIC_F2yz_Fy2z_vrr = PBY*I_KINETIC_F2yz_D2z_vrr+2*oned2z*I_KINETIC_Dyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fy2z_vrr;
    Double I_KINETIC_Fy2z_Fy2z_vrr = PBY*I_KINETIC_Fy2z_D2z_vrr+oned2z*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Fy2z_vrr;
    Double I_KINETIC_F3z_Fy2z_vrr = PBY*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_KINETIC_F3x_F3z_vrr = PBZ*I_KINETIC_F3x_D2z_vrr+2*oned2z*I_KINETIC_F3x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F3z_vrr-2*adz*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_KINETIC_F2xy_F3z_vrr = PBZ*I_KINETIC_F2xy_D2z_vrr+2*oned2z*I_KINETIC_F2xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F3z_vrr-2*adz*I_TWOBODYOVERLAP_F2xy_Pz_vrr;
    Double I_KINETIC_F2xz_F3z_vrr = PBZ*I_KINETIC_F2xz_D2z_vrr+oned2z*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_F2xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_F2xz_Pz_vrr;
    Double I_KINETIC_Fx2y_F3z_vrr = PBZ*I_KINETIC_Fx2y_D2z_vrr+2*oned2z*I_KINETIC_Fx2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Fx2y_Pz_vrr;
    Double I_KINETIC_Fxyz_F3z_vrr = PBZ*I_KINETIC_Fxyz_D2z_vrr+oned2z*I_KINETIC_Dxy_D2z_vrr+2*oned2z*I_KINETIC_Fxyz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Fxyz_Pz_vrr;
    Double I_KINETIC_Fx2z_F3z_vrr = PBZ*I_KINETIC_Fx2z_D2z_vrr+2*oned2z*I_KINETIC_Dxz_D2z_vrr+2*oned2z*I_KINETIC_Fx2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Fx2z_Pz_vrr;
    Double I_KINETIC_F3y_F3z_vrr = PBZ*I_KINETIC_F3y_D2z_vrr+2*oned2z*I_KINETIC_F3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_KINETIC_F2yz_F3z_vrr = PBZ*I_KINETIC_F2yz_D2z_vrr+oned2z*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_F2yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_F2yz_Pz_vrr;
    Double I_KINETIC_Fy2z_F3z_vrr = PBZ*I_KINETIC_Fy2z_D2z_vrr+2*oned2z*I_KINETIC_Dyz_D2z_vrr+2*oned2z*I_KINETIC_Fy2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Fy2z_Pz_vrr;
    Double I_KINETIC_F3z_F3z_vrr = PBZ*I_KINETIC_F3z_D2z_vrr+3*oned2z*I_KINETIC_D2z_D2z_vrr+2*oned2z*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_F3z_Pz_vrr;

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
     * shell quartet name: SQ_KINETIC_H_F_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_H_F_a_coefs = alpha;
    I_KINETIC_H5x_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_F3x_vrr;
    I_KINETIC_H4xy_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_F3x_vrr;
    I_KINETIC_H4xz_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_F3x_vrr;
    I_KINETIC_H3x2y_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_F3x_vrr;
    I_KINETIC_H3xyz_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_F3x_vrr;
    I_KINETIC_H3x2z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_F3x_vrr;
    I_KINETIC_H2x3y_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_F3x_vrr;
    I_KINETIC_H2x2yz_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_F3x_vrr;
    I_KINETIC_H2xy2z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_F3x_vrr;
    I_KINETIC_H2x3z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_F3x_vrr;
    I_KINETIC_Hx4y_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_F3x_vrr;
    I_KINETIC_Hx3yz_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_F3x_vrr;
    I_KINETIC_Hx2y2z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_F3x_vrr;
    I_KINETIC_Hxy3z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_F3x_vrr;
    I_KINETIC_Hx4z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_F3x_vrr;
    I_KINETIC_H5y_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_F3x_vrr;
    I_KINETIC_H4yz_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_F3x_vrr;
    I_KINETIC_H3y2z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_F3x_vrr;
    I_KINETIC_H2y3z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_F3x_vrr;
    I_KINETIC_Hy4z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_F3x_vrr;
    I_KINETIC_H5z_F3x_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_F3x_vrr;
    I_KINETIC_H5x_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_F2xy_vrr;
    I_KINETIC_H4xy_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_F2xy_vrr;
    I_KINETIC_H4xz_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_F2xy_vrr;
    I_KINETIC_H3x2y_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_F2xy_vrr;
    I_KINETIC_H3xyz_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_F2xy_vrr;
    I_KINETIC_H3x2z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_F2xy_vrr;
    I_KINETIC_H2x3y_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_F2xy_vrr;
    I_KINETIC_H2x2yz_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_F2xy_vrr;
    I_KINETIC_H2xy2z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_F2xy_vrr;
    I_KINETIC_H2x3z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_F2xy_vrr;
    I_KINETIC_Hx4y_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_F2xy_vrr;
    I_KINETIC_Hx3yz_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_F2xy_vrr;
    I_KINETIC_Hx2y2z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_F2xy_vrr;
    I_KINETIC_Hxy3z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_F2xy_vrr;
    I_KINETIC_Hx4z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_F2xy_vrr;
    I_KINETIC_H5y_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_F2xy_vrr;
    I_KINETIC_H4yz_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_F2xy_vrr;
    I_KINETIC_H3y2z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_F2xy_vrr;
    I_KINETIC_H2y3z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_F2xy_vrr;
    I_KINETIC_Hy4z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_F2xy_vrr;
    I_KINETIC_H5z_F2xy_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_F2xy_vrr;
    I_KINETIC_H5x_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_F2xz_vrr;
    I_KINETIC_H4xy_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_F2xz_vrr;
    I_KINETIC_H4xz_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_F2xz_vrr;
    I_KINETIC_H3x2y_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_F2xz_vrr;
    I_KINETIC_H3xyz_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_F2xz_vrr;
    I_KINETIC_H3x2z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_F2xz_vrr;
    I_KINETIC_H2x3y_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_F2xz_vrr;
    I_KINETIC_H2x2yz_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_F2xz_vrr;
    I_KINETIC_H2xy2z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_F2xz_vrr;
    I_KINETIC_H2x3z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_F2xz_vrr;
    I_KINETIC_Hx4y_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_F2xz_vrr;
    I_KINETIC_Hx3yz_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_F2xz_vrr;
    I_KINETIC_Hx2y2z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_F2xz_vrr;
    I_KINETIC_Hxy3z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_F2xz_vrr;
    I_KINETIC_Hx4z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_F2xz_vrr;
    I_KINETIC_H5y_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_F2xz_vrr;
    I_KINETIC_H4yz_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_F2xz_vrr;
    I_KINETIC_H3y2z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_F2xz_vrr;
    I_KINETIC_H2y3z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_F2xz_vrr;
    I_KINETIC_Hy4z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_F2xz_vrr;
    I_KINETIC_H5z_F2xz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_F2xz_vrr;
    I_KINETIC_H5x_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_Fx2y_vrr;
    I_KINETIC_H4xy_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_Fx2y_vrr;
    I_KINETIC_H4xz_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_Fx2y_vrr;
    I_KINETIC_H3x2y_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_Fx2y_vrr;
    I_KINETIC_H3xyz_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_Fx2y_vrr;
    I_KINETIC_H3x2z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_Fx2y_vrr;
    I_KINETIC_H2x3y_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_Fx2y_vrr;
    I_KINETIC_H2x2yz_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_Fx2y_vrr;
    I_KINETIC_H2xy2z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_Fx2y_vrr;
    I_KINETIC_H2x3z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_Fx2y_vrr;
    I_KINETIC_Hx4y_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_Fx2y_vrr;
    I_KINETIC_Hx3yz_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_Fx2y_vrr;
    I_KINETIC_Hx2y2z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_Fx2y_vrr;
    I_KINETIC_Hxy3z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_Fx2y_vrr;
    I_KINETIC_Hx4z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_Fx2y_vrr;
    I_KINETIC_H5y_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_Fx2y_vrr;
    I_KINETIC_H4yz_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_Fx2y_vrr;
    I_KINETIC_H3y2z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_Fx2y_vrr;
    I_KINETIC_H2y3z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_Fx2y_vrr;
    I_KINETIC_Hy4z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_Fx2y_vrr;
    I_KINETIC_H5z_Fx2y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_Fx2y_vrr;
    I_KINETIC_H5x_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_Fxyz_vrr;
    I_KINETIC_H4xy_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_Fxyz_vrr;
    I_KINETIC_H4xz_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_Fxyz_vrr;
    I_KINETIC_H3x2y_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_Fxyz_vrr;
    I_KINETIC_H3xyz_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_Fxyz_vrr;
    I_KINETIC_H3x2z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_Fxyz_vrr;
    I_KINETIC_H2x3y_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_Fxyz_vrr;
    I_KINETIC_H2x2yz_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_Fxyz_vrr;
    I_KINETIC_H2xy2z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_Fxyz_vrr;
    I_KINETIC_H2x3z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_Fxyz_vrr;
    I_KINETIC_Hx4y_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_Fxyz_vrr;
    I_KINETIC_Hx3yz_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_Fxyz_vrr;
    I_KINETIC_Hx2y2z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_Fxyz_vrr;
    I_KINETIC_Hxy3z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_Fxyz_vrr;
    I_KINETIC_Hx4z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_Fxyz_vrr;
    I_KINETIC_H5y_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_Fxyz_vrr;
    I_KINETIC_H4yz_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_Fxyz_vrr;
    I_KINETIC_H3y2z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_Fxyz_vrr;
    I_KINETIC_H2y3z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_Fxyz_vrr;
    I_KINETIC_Hy4z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_Fxyz_vrr;
    I_KINETIC_H5z_Fxyz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_Fxyz_vrr;
    I_KINETIC_H5x_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_Fx2z_vrr;
    I_KINETIC_H4xy_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_Fx2z_vrr;
    I_KINETIC_H4xz_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_Fx2z_vrr;
    I_KINETIC_H3x2y_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_Fx2z_vrr;
    I_KINETIC_H3xyz_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_Fx2z_vrr;
    I_KINETIC_H3x2z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_Fx2z_vrr;
    I_KINETIC_H2x3y_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_Fx2z_vrr;
    I_KINETIC_H2x2yz_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_Fx2z_vrr;
    I_KINETIC_H2xy2z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_Fx2z_vrr;
    I_KINETIC_H2x3z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_Fx2z_vrr;
    I_KINETIC_Hx4y_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_Fx2z_vrr;
    I_KINETIC_Hx3yz_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_Fx2z_vrr;
    I_KINETIC_Hx2y2z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_Fx2z_vrr;
    I_KINETIC_Hxy3z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_Fx2z_vrr;
    I_KINETIC_Hx4z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_Fx2z_vrr;
    I_KINETIC_H5y_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_Fx2z_vrr;
    I_KINETIC_H4yz_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_Fx2z_vrr;
    I_KINETIC_H3y2z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_Fx2z_vrr;
    I_KINETIC_H2y3z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_Fx2z_vrr;
    I_KINETIC_Hy4z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_Fx2z_vrr;
    I_KINETIC_H5z_Fx2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_Fx2z_vrr;
    I_KINETIC_H5x_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_F3y_vrr;
    I_KINETIC_H4xy_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_F3y_vrr;
    I_KINETIC_H4xz_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_F3y_vrr;
    I_KINETIC_H3x2y_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_F3y_vrr;
    I_KINETIC_H3xyz_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_F3y_vrr;
    I_KINETIC_H3x2z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_F3y_vrr;
    I_KINETIC_H2x3y_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_F3y_vrr;
    I_KINETIC_H2x2yz_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_F3y_vrr;
    I_KINETIC_H2xy2z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_F3y_vrr;
    I_KINETIC_H2x3z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_F3y_vrr;
    I_KINETIC_Hx4y_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_F3y_vrr;
    I_KINETIC_Hx3yz_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_F3y_vrr;
    I_KINETIC_Hx2y2z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_F3y_vrr;
    I_KINETIC_Hxy3z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_F3y_vrr;
    I_KINETIC_Hx4z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_F3y_vrr;
    I_KINETIC_H5y_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_F3y_vrr;
    I_KINETIC_H4yz_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_F3y_vrr;
    I_KINETIC_H3y2z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_F3y_vrr;
    I_KINETIC_H2y3z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_F3y_vrr;
    I_KINETIC_Hy4z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_F3y_vrr;
    I_KINETIC_H5z_F3y_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_F3y_vrr;
    I_KINETIC_H5x_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_F2yz_vrr;
    I_KINETIC_H4xy_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_F2yz_vrr;
    I_KINETIC_H4xz_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_F2yz_vrr;
    I_KINETIC_H3x2y_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_F2yz_vrr;
    I_KINETIC_H3xyz_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_F2yz_vrr;
    I_KINETIC_H3x2z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_F2yz_vrr;
    I_KINETIC_H2x3y_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_F2yz_vrr;
    I_KINETIC_H2x2yz_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_F2yz_vrr;
    I_KINETIC_H2xy2z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_F2yz_vrr;
    I_KINETIC_H2x3z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_F2yz_vrr;
    I_KINETIC_Hx4y_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_F2yz_vrr;
    I_KINETIC_Hx3yz_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_F2yz_vrr;
    I_KINETIC_Hx2y2z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_F2yz_vrr;
    I_KINETIC_Hxy3z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_F2yz_vrr;
    I_KINETIC_Hx4z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_F2yz_vrr;
    I_KINETIC_H5y_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_F2yz_vrr;
    I_KINETIC_H4yz_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_F2yz_vrr;
    I_KINETIC_H3y2z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_F2yz_vrr;
    I_KINETIC_H2y3z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_F2yz_vrr;
    I_KINETIC_Hy4z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_F2yz_vrr;
    I_KINETIC_H5z_F2yz_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_F2yz_vrr;
    I_KINETIC_H5x_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_Fy2z_vrr;
    I_KINETIC_H4xy_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_Fy2z_vrr;
    I_KINETIC_H4xz_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_Fy2z_vrr;
    I_KINETIC_H3x2y_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_Fy2z_vrr;
    I_KINETIC_H3xyz_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_Fy2z_vrr;
    I_KINETIC_H3x2z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_Fy2z_vrr;
    I_KINETIC_H2x3y_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_Fy2z_vrr;
    I_KINETIC_H2x2yz_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_Fy2z_vrr;
    I_KINETIC_H2xy2z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_Fy2z_vrr;
    I_KINETIC_H2x3z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_Fy2z_vrr;
    I_KINETIC_Hx4y_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_Fy2z_vrr;
    I_KINETIC_Hx3yz_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_Fy2z_vrr;
    I_KINETIC_Hx2y2z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_Fy2z_vrr;
    I_KINETIC_Hxy3z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_Fy2z_vrr;
    I_KINETIC_Hx4z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_Fy2z_vrr;
    I_KINETIC_H5y_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_Fy2z_vrr;
    I_KINETIC_H4yz_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_Fy2z_vrr;
    I_KINETIC_H3y2z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_Fy2z_vrr;
    I_KINETIC_H2y3z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_Fy2z_vrr;
    I_KINETIC_Hy4z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_Fy2z_vrr;
    I_KINETIC_H5z_Fy2z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_Fy2z_vrr;
    I_KINETIC_H5x_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5x_F3z_vrr;
    I_KINETIC_H4xy_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xy_F3z_vrr;
    I_KINETIC_H4xz_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4xz_F3z_vrr;
    I_KINETIC_H3x2y_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2y_F3z_vrr;
    I_KINETIC_H3xyz_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3xyz_F3z_vrr;
    I_KINETIC_H3x2z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3x2z_F3z_vrr;
    I_KINETIC_H2x3y_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3y_F3z_vrr;
    I_KINETIC_H2x2yz_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x2yz_F3z_vrr;
    I_KINETIC_H2xy2z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2xy2z_F3z_vrr;
    I_KINETIC_H2x3z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2x3z_F3z_vrr;
    I_KINETIC_Hx4y_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4y_F3z_vrr;
    I_KINETIC_Hx3yz_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx3yz_F3z_vrr;
    I_KINETIC_Hx2y2z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx2y2z_F3z_vrr;
    I_KINETIC_Hxy3z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hxy3z_F3z_vrr;
    I_KINETIC_Hx4z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hx4z_F3z_vrr;
    I_KINETIC_H5y_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5y_F3z_vrr;
    I_KINETIC_H4yz_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H4yz_F3z_vrr;
    I_KINETIC_H3y2z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H3y2z_F3z_vrr;
    I_KINETIC_H2y3z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H2y3z_F3z_vrr;
    I_KINETIC_Hy4z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_Hy4z_F3z_vrr;
    I_KINETIC_H5z_F3z_a += SQ_KINETIC_H_F_a_coefs*I_KINETIC_H5z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_F
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_F3x_F3x += I_KINETIC_F3x_F3x_vrr;
    I_KINETIC_F2xy_F3x += I_KINETIC_F2xy_F3x_vrr;
    I_KINETIC_F2xz_F3x += I_KINETIC_F2xz_F3x_vrr;
    I_KINETIC_Fx2y_F3x += I_KINETIC_Fx2y_F3x_vrr;
    I_KINETIC_Fxyz_F3x += I_KINETIC_Fxyz_F3x_vrr;
    I_KINETIC_Fx2z_F3x += I_KINETIC_Fx2z_F3x_vrr;
    I_KINETIC_F3y_F3x += I_KINETIC_F3y_F3x_vrr;
    I_KINETIC_F2yz_F3x += I_KINETIC_F2yz_F3x_vrr;
    I_KINETIC_Fy2z_F3x += I_KINETIC_Fy2z_F3x_vrr;
    I_KINETIC_F3z_F3x += I_KINETIC_F3z_F3x_vrr;
    I_KINETIC_F3x_F2xy += I_KINETIC_F3x_F2xy_vrr;
    I_KINETIC_F2xy_F2xy += I_KINETIC_F2xy_F2xy_vrr;
    I_KINETIC_F2xz_F2xy += I_KINETIC_F2xz_F2xy_vrr;
    I_KINETIC_Fx2y_F2xy += I_KINETIC_Fx2y_F2xy_vrr;
    I_KINETIC_Fxyz_F2xy += I_KINETIC_Fxyz_F2xy_vrr;
    I_KINETIC_Fx2z_F2xy += I_KINETIC_Fx2z_F2xy_vrr;
    I_KINETIC_F3y_F2xy += I_KINETIC_F3y_F2xy_vrr;
    I_KINETIC_F2yz_F2xy += I_KINETIC_F2yz_F2xy_vrr;
    I_KINETIC_Fy2z_F2xy += I_KINETIC_Fy2z_F2xy_vrr;
    I_KINETIC_F3z_F2xy += I_KINETIC_F3z_F2xy_vrr;
    I_KINETIC_F3x_F2xz += I_KINETIC_F3x_F2xz_vrr;
    I_KINETIC_F2xy_F2xz += I_KINETIC_F2xy_F2xz_vrr;
    I_KINETIC_F2xz_F2xz += I_KINETIC_F2xz_F2xz_vrr;
    I_KINETIC_Fx2y_F2xz += I_KINETIC_Fx2y_F2xz_vrr;
    I_KINETIC_Fxyz_F2xz += I_KINETIC_Fxyz_F2xz_vrr;
    I_KINETIC_Fx2z_F2xz += I_KINETIC_Fx2z_F2xz_vrr;
    I_KINETIC_F3y_F2xz += I_KINETIC_F3y_F2xz_vrr;
    I_KINETIC_F2yz_F2xz += I_KINETIC_F2yz_F2xz_vrr;
    I_KINETIC_Fy2z_F2xz += I_KINETIC_Fy2z_F2xz_vrr;
    I_KINETIC_F3z_F2xz += I_KINETIC_F3z_F2xz_vrr;
    I_KINETIC_F3x_Fx2y += I_KINETIC_F3x_Fx2y_vrr;
    I_KINETIC_F2xy_Fx2y += I_KINETIC_F2xy_Fx2y_vrr;
    I_KINETIC_F2xz_Fx2y += I_KINETIC_F2xz_Fx2y_vrr;
    I_KINETIC_Fx2y_Fx2y += I_KINETIC_Fx2y_Fx2y_vrr;
    I_KINETIC_Fxyz_Fx2y += I_KINETIC_Fxyz_Fx2y_vrr;
    I_KINETIC_Fx2z_Fx2y += I_KINETIC_Fx2z_Fx2y_vrr;
    I_KINETIC_F3y_Fx2y += I_KINETIC_F3y_Fx2y_vrr;
    I_KINETIC_F2yz_Fx2y += I_KINETIC_F2yz_Fx2y_vrr;
    I_KINETIC_Fy2z_Fx2y += I_KINETIC_Fy2z_Fx2y_vrr;
    I_KINETIC_F3z_Fx2y += I_KINETIC_F3z_Fx2y_vrr;
    I_KINETIC_F3x_Fxyz += I_KINETIC_F3x_Fxyz_vrr;
    I_KINETIC_F2xy_Fxyz += I_KINETIC_F2xy_Fxyz_vrr;
    I_KINETIC_F2xz_Fxyz += I_KINETIC_F2xz_Fxyz_vrr;
    I_KINETIC_Fx2y_Fxyz += I_KINETIC_Fx2y_Fxyz_vrr;
    I_KINETIC_Fxyz_Fxyz += I_KINETIC_Fxyz_Fxyz_vrr;
    I_KINETIC_Fx2z_Fxyz += I_KINETIC_Fx2z_Fxyz_vrr;
    I_KINETIC_F3y_Fxyz += I_KINETIC_F3y_Fxyz_vrr;
    I_KINETIC_F2yz_Fxyz += I_KINETIC_F2yz_Fxyz_vrr;
    I_KINETIC_Fy2z_Fxyz += I_KINETIC_Fy2z_Fxyz_vrr;
    I_KINETIC_F3z_Fxyz += I_KINETIC_F3z_Fxyz_vrr;
    I_KINETIC_F3x_Fx2z += I_KINETIC_F3x_Fx2z_vrr;
    I_KINETIC_F2xy_Fx2z += I_KINETIC_F2xy_Fx2z_vrr;
    I_KINETIC_F2xz_Fx2z += I_KINETIC_F2xz_Fx2z_vrr;
    I_KINETIC_Fx2y_Fx2z += I_KINETIC_Fx2y_Fx2z_vrr;
    I_KINETIC_Fxyz_Fx2z += I_KINETIC_Fxyz_Fx2z_vrr;
    I_KINETIC_Fx2z_Fx2z += I_KINETIC_Fx2z_Fx2z_vrr;
    I_KINETIC_F3y_Fx2z += I_KINETIC_F3y_Fx2z_vrr;
    I_KINETIC_F2yz_Fx2z += I_KINETIC_F2yz_Fx2z_vrr;
    I_KINETIC_Fy2z_Fx2z += I_KINETIC_Fy2z_Fx2z_vrr;
    I_KINETIC_F3z_Fx2z += I_KINETIC_F3z_Fx2z_vrr;
    I_KINETIC_F3x_F3y += I_KINETIC_F3x_F3y_vrr;
    I_KINETIC_F2xy_F3y += I_KINETIC_F2xy_F3y_vrr;
    I_KINETIC_F2xz_F3y += I_KINETIC_F2xz_F3y_vrr;
    I_KINETIC_Fx2y_F3y += I_KINETIC_Fx2y_F3y_vrr;
    I_KINETIC_Fxyz_F3y += I_KINETIC_Fxyz_F3y_vrr;
    I_KINETIC_Fx2z_F3y += I_KINETIC_Fx2z_F3y_vrr;
    I_KINETIC_F3y_F3y += I_KINETIC_F3y_F3y_vrr;
    I_KINETIC_F2yz_F3y += I_KINETIC_F2yz_F3y_vrr;
    I_KINETIC_Fy2z_F3y += I_KINETIC_Fy2z_F3y_vrr;
    I_KINETIC_F3z_F3y += I_KINETIC_F3z_F3y_vrr;
    I_KINETIC_F3x_F2yz += I_KINETIC_F3x_F2yz_vrr;
    I_KINETIC_F2xy_F2yz += I_KINETIC_F2xy_F2yz_vrr;
    I_KINETIC_F2xz_F2yz += I_KINETIC_F2xz_F2yz_vrr;
    I_KINETIC_Fx2y_F2yz += I_KINETIC_Fx2y_F2yz_vrr;
    I_KINETIC_Fxyz_F2yz += I_KINETIC_Fxyz_F2yz_vrr;
    I_KINETIC_Fx2z_F2yz += I_KINETIC_Fx2z_F2yz_vrr;
    I_KINETIC_F3y_F2yz += I_KINETIC_F3y_F2yz_vrr;
    I_KINETIC_F2yz_F2yz += I_KINETIC_F2yz_F2yz_vrr;
    I_KINETIC_Fy2z_F2yz += I_KINETIC_Fy2z_F2yz_vrr;
    I_KINETIC_F3z_F2yz += I_KINETIC_F3z_F2yz_vrr;
    I_KINETIC_F3x_Fy2z += I_KINETIC_F3x_Fy2z_vrr;
    I_KINETIC_F2xy_Fy2z += I_KINETIC_F2xy_Fy2z_vrr;
    I_KINETIC_F2xz_Fy2z += I_KINETIC_F2xz_Fy2z_vrr;
    I_KINETIC_Fx2y_Fy2z += I_KINETIC_Fx2y_Fy2z_vrr;
    I_KINETIC_Fxyz_Fy2z += I_KINETIC_Fxyz_Fy2z_vrr;
    I_KINETIC_Fx2z_Fy2z += I_KINETIC_Fx2z_Fy2z_vrr;
    I_KINETIC_F3y_Fy2z += I_KINETIC_F3y_Fy2z_vrr;
    I_KINETIC_F2yz_Fy2z += I_KINETIC_F2yz_Fy2z_vrr;
    I_KINETIC_Fy2z_Fy2z += I_KINETIC_Fy2z_Fy2z_vrr;
    I_KINETIC_F3z_Fy2z += I_KINETIC_F3z_Fy2z_vrr;
    I_KINETIC_F3x_F3z += I_KINETIC_F3x_F3z_vrr;
    I_KINETIC_F2xy_F3z += I_KINETIC_F2xy_F3z_vrr;
    I_KINETIC_F2xz_F3z += I_KINETIC_F2xz_F3z_vrr;
    I_KINETIC_Fx2y_F3z += I_KINETIC_Fx2y_F3z_vrr;
    I_KINETIC_Fxyz_F3z += I_KINETIC_Fxyz_F3z_vrr;
    I_KINETIC_Fx2z_F3z += I_KINETIC_Fx2z_F3z_vrr;
    I_KINETIC_F3y_F3z += I_KINETIC_F3y_F3z_vrr;
    I_KINETIC_F2yz_F3z += I_KINETIC_F2yz_F3z_vrr;
    I_KINETIC_Fy2z_F3z += I_KINETIC_Fy2z_F3z_vrr;
    I_KINETIC_F3z_F3z += I_KINETIC_F3z_F3z_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_F_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_F_a
   * RHS shell quartet name: SQ_KINETIC_F_F
   ************************************************************/
  abcd[0] = 2.0E0*I_KINETIC_H5x_F3x_a-4*I_KINETIC_F3x_F3x;
  abcd[1] = 2.0E0*I_KINETIC_H4xy_F3x_a-3*I_KINETIC_F2xy_F3x;
  abcd[2] = 2.0E0*I_KINETIC_H4xz_F3x_a-3*I_KINETIC_F2xz_F3x;
  abcd[3] = 2.0E0*I_KINETIC_H3x2y_F3x_a-2*I_KINETIC_Fx2y_F3x;
  abcd[4] = 2.0E0*I_KINETIC_H3xyz_F3x_a-2*I_KINETIC_Fxyz_F3x;
  abcd[5] = 2.0E0*I_KINETIC_H3x2z_F3x_a-2*I_KINETIC_Fx2z_F3x;
  abcd[6] = 2.0E0*I_KINETIC_H2x3y_F3x_a-1*I_KINETIC_F3y_F3x;
  abcd[7] = 2.0E0*I_KINETIC_H2x2yz_F3x_a-1*I_KINETIC_F2yz_F3x;
  abcd[8] = 2.0E0*I_KINETIC_H2xy2z_F3x_a-1*I_KINETIC_Fy2z_F3x;
  abcd[9] = 2.0E0*I_KINETIC_H2x3z_F3x_a-1*I_KINETIC_F3z_F3x;
  abcd[10] = 2.0E0*I_KINETIC_Hx4y_F3x_a;
  abcd[11] = 2.0E0*I_KINETIC_Hx3yz_F3x_a;
  abcd[12] = 2.0E0*I_KINETIC_Hx2y2z_F3x_a;
  abcd[13] = 2.0E0*I_KINETIC_Hxy3z_F3x_a;
  abcd[14] = 2.0E0*I_KINETIC_Hx4z_F3x_a;
  abcd[15] = 2.0E0*I_KINETIC_H5x_F2xy_a-4*I_KINETIC_F3x_F2xy;
  abcd[16] = 2.0E0*I_KINETIC_H4xy_F2xy_a-3*I_KINETIC_F2xy_F2xy;
  abcd[17] = 2.0E0*I_KINETIC_H4xz_F2xy_a-3*I_KINETIC_F2xz_F2xy;
  abcd[18] = 2.0E0*I_KINETIC_H3x2y_F2xy_a-2*I_KINETIC_Fx2y_F2xy;
  abcd[19] = 2.0E0*I_KINETIC_H3xyz_F2xy_a-2*I_KINETIC_Fxyz_F2xy;
  abcd[20] = 2.0E0*I_KINETIC_H3x2z_F2xy_a-2*I_KINETIC_Fx2z_F2xy;
  abcd[21] = 2.0E0*I_KINETIC_H2x3y_F2xy_a-1*I_KINETIC_F3y_F2xy;
  abcd[22] = 2.0E0*I_KINETIC_H2x2yz_F2xy_a-1*I_KINETIC_F2yz_F2xy;
  abcd[23] = 2.0E0*I_KINETIC_H2xy2z_F2xy_a-1*I_KINETIC_Fy2z_F2xy;
  abcd[24] = 2.0E0*I_KINETIC_H2x3z_F2xy_a-1*I_KINETIC_F3z_F2xy;
  abcd[25] = 2.0E0*I_KINETIC_Hx4y_F2xy_a;
  abcd[26] = 2.0E0*I_KINETIC_Hx3yz_F2xy_a;
  abcd[27] = 2.0E0*I_KINETIC_Hx2y2z_F2xy_a;
  abcd[28] = 2.0E0*I_KINETIC_Hxy3z_F2xy_a;
  abcd[29] = 2.0E0*I_KINETIC_Hx4z_F2xy_a;
  abcd[30] = 2.0E0*I_KINETIC_H5x_F2xz_a-4*I_KINETIC_F3x_F2xz;
  abcd[31] = 2.0E0*I_KINETIC_H4xy_F2xz_a-3*I_KINETIC_F2xy_F2xz;
  abcd[32] = 2.0E0*I_KINETIC_H4xz_F2xz_a-3*I_KINETIC_F2xz_F2xz;
  abcd[33] = 2.0E0*I_KINETIC_H3x2y_F2xz_a-2*I_KINETIC_Fx2y_F2xz;
  abcd[34] = 2.0E0*I_KINETIC_H3xyz_F2xz_a-2*I_KINETIC_Fxyz_F2xz;
  abcd[35] = 2.0E0*I_KINETIC_H3x2z_F2xz_a-2*I_KINETIC_Fx2z_F2xz;
  abcd[36] = 2.0E0*I_KINETIC_H2x3y_F2xz_a-1*I_KINETIC_F3y_F2xz;
  abcd[37] = 2.0E0*I_KINETIC_H2x2yz_F2xz_a-1*I_KINETIC_F2yz_F2xz;
  abcd[38] = 2.0E0*I_KINETIC_H2xy2z_F2xz_a-1*I_KINETIC_Fy2z_F2xz;
  abcd[39] = 2.0E0*I_KINETIC_H2x3z_F2xz_a-1*I_KINETIC_F3z_F2xz;
  abcd[40] = 2.0E0*I_KINETIC_Hx4y_F2xz_a;
  abcd[41] = 2.0E0*I_KINETIC_Hx3yz_F2xz_a;
  abcd[42] = 2.0E0*I_KINETIC_Hx2y2z_F2xz_a;
  abcd[43] = 2.0E0*I_KINETIC_Hxy3z_F2xz_a;
  abcd[44] = 2.0E0*I_KINETIC_Hx4z_F2xz_a;
  abcd[45] = 2.0E0*I_KINETIC_H5x_Fx2y_a-4*I_KINETIC_F3x_Fx2y;
  abcd[46] = 2.0E0*I_KINETIC_H4xy_Fx2y_a-3*I_KINETIC_F2xy_Fx2y;
  abcd[47] = 2.0E0*I_KINETIC_H4xz_Fx2y_a-3*I_KINETIC_F2xz_Fx2y;
  abcd[48] = 2.0E0*I_KINETIC_H3x2y_Fx2y_a-2*I_KINETIC_Fx2y_Fx2y;
  abcd[49] = 2.0E0*I_KINETIC_H3xyz_Fx2y_a-2*I_KINETIC_Fxyz_Fx2y;
  abcd[50] = 2.0E0*I_KINETIC_H3x2z_Fx2y_a-2*I_KINETIC_Fx2z_Fx2y;
  abcd[51] = 2.0E0*I_KINETIC_H2x3y_Fx2y_a-1*I_KINETIC_F3y_Fx2y;
  abcd[52] = 2.0E0*I_KINETIC_H2x2yz_Fx2y_a-1*I_KINETIC_F2yz_Fx2y;
  abcd[53] = 2.0E0*I_KINETIC_H2xy2z_Fx2y_a-1*I_KINETIC_Fy2z_Fx2y;
  abcd[54] = 2.0E0*I_KINETIC_H2x3z_Fx2y_a-1*I_KINETIC_F3z_Fx2y;
  abcd[55] = 2.0E0*I_KINETIC_Hx4y_Fx2y_a;
  abcd[56] = 2.0E0*I_KINETIC_Hx3yz_Fx2y_a;
  abcd[57] = 2.0E0*I_KINETIC_Hx2y2z_Fx2y_a;
  abcd[58] = 2.0E0*I_KINETIC_Hxy3z_Fx2y_a;
  abcd[59] = 2.0E0*I_KINETIC_Hx4z_Fx2y_a;
  abcd[60] = 2.0E0*I_KINETIC_H5x_Fxyz_a-4*I_KINETIC_F3x_Fxyz;
  abcd[61] = 2.0E0*I_KINETIC_H4xy_Fxyz_a-3*I_KINETIC_F2xy_Fxyz;
  abcd[62] = 2.0E0*I_KINETIC_H4xz_Fxyz_a-3*I_KINETIC_F2xz_Fxyz;
  abcd[63] = 2.0E0*I_KINETIC_H3x2y_Fxyz_a-2*I_KINETIC_Fx2y_Fxyz;
  abcd[64] = 2.0E0*I_KINETIC_H3xyz_Fxyz_a-2*I_KINETIC_Fxyz_Fxyz;
  abcd[65] = 2.0E0*I_KINETIC_H3x2z_Fxyz_a-2*I_KINETIC_Fx2z_Fxyz;
  abcd[66] = 2.0E0*I_KINETIC_H2x3y_Fxyz_a-1*I_KINETIC_F3y_Fxyz;
  abcd[67] = 2.0E0*I_KINETIC_H2x2yz_Fxyz_a-1*I_KINETIC_F2yz_Fxyz;
  abcd[68] = 2.0E0*I_KINETIC_H2xy2z_Fxyz_a-1*I_KINETIC_Fy2z_Fxyz;
  abcd[69] = 2.0E0*I_KINETIC_H2x3z_Fxyz_a-1*I_KINETIC_F3z_Fxyz;
  abcd[70] = 2.0E0*I_KINETIC_Hx4y_Fxyz_a;
  abcd[71] = 2.0E0*I_KINETIC_Hx3yz_Fxyz_a;
  abcd[72] = 2.0E0*I_KINETIC_Hx2y2z_Fxyz_a;
  abcd[73] = 2.0E0*I_KINETIC_Hxy3z_Fxyz_a;
  abcd[74] = 2.0E0*I_KINETIC_Hx4z_Fxyz_a;
  abcd[75] = 2.0E0*I_KINETIC_H5x_Fx2z_a-4*I_KINETIC_F3x_Fx2z;
  abcd[76] = 2.0E0*I_KINETIC_H4xy_Fx2z_a-3*I_KINETIC_F2xy_Fx2z;
  abcd[77] = 2.0E0*I_KINETIC_H4xz_Fx2z_a-3*I_KINETIC_F2xz_Fx2z;
  abcd[78] = 2.0E0*I_KINETIC_H3x2y_Fx2z_a-2*I_KINETIC_Fx2y_Fx2z;
  abcd[79] = 2.0E0*I_KINETIC_H3xyz_Fx2z_a-2*I_KINETIC_Fxyz_Fx2z;
  abcd[80] = 2.0E0*I_KINETIC_H3x2z_Fx2z_a-2*I_KINETIC_Fx2z_Fx2z;
  abcd[81] = 2.0E0*I_KINETIC_H2x3y_Fx2z_a-1*I_KINETIC_F3y_Fx2z;
  abcd[82] = 2.0E0*I_KINETIC_H2x2yz_Fx2z_a-1*I_KINETIC_F2yz_Fx2z;
  abcd[83] = 2.0E0*I_KINETIC_H2xy2z_Fx2z_a-1*I_KINETIC_Fy2z_Fx2z;
  abcd[84] = 2.0E0*I_KINETIC_H2x3z_Fx2z_a-1*I_KINETIC_F3z_Fx2z;
  abcd[85] = 2.0E0*I_KINETIC_Hx4y_Fx2z_a;
  abcd[86] = 2.0E0*I_KINETIC_Hx3yz_Fx2z_a;
  abcd[87] = 2.0E0*I_KINETIC_Hx2y2z_Fx2z_a;
  abcd[88] = 2.0E0*I_KINETIC_Hxy3z_Fx2z_a;
  abcd[89] = 2.0E0*I_KINETIC_Hx4z_Fx2z_a;
  abcd[90] = 2.0E0*I_KINETIC_H5x_F3y_a-4*I_KINETIC_F3x_F3y;
  abcd[91] = 2.0E0*I_KINETIC_H4xy_F3y_a-3*I_KINETIC_F2xy_F3y;
  abcd[92] = 2.0E0*I_KINETIC_H4xz_F3y_a-3*I_KINETIC_F2xz_F3y;
  abcd[93] = 2.0E0*I_KINETIC_H3x2y_F3y_a-2*I_KINETIC_Fx2y_F3y;
  abcd[94] = 2.0E0*I_KINETIC_H3xyz_F3y_a-2*I_KINETIC_Fxyz_F3y;
  abcd[95] = 2.0E0*I_KINETIC_H3x2z_F3y_a-2*I_KINETIC_Fx2z_F3y;
  abcd[96] = 2.0E0*I_KINETIC_H2x3y_F3y_a-1*I_KINETIC_F3y_F3y;
  abcd[97] = 2.0E0*I_KINETIC_H2x2yz_F3y_a-1*I_KINETIC_F2yz_F3y;
  abcd[98] = 2.0E0*I_KINETIC_H2xy2z_F3y_a-1*I_KINETIC_Fy2z_F3y;
  abcd[99] = 2.0E0*I_KINETIC_H2x3z_F3y_a-1*I_KINETIC_F3z_F3y;
  abcd[100] = 2.0E0*I_KINETIC_Hx4y_F3y_a;
  abcd[101] = 2.0E0*I_KINETIC_Hx3yz_F3y_a;
  abcd[102] = 2.0E0*I_KINETIC_Hx2y2z_F3y_a;
  abcd[103] = 2.0E0*I_KINETIC_Hxy3z_F3y_a;
  abcd[104] = 2.0E0*I_KINETIC_Hx4z_F3y_a;
  abcd[105] = 2.0E0*I_KINETIC_H5x_F2yz_a-4*I_KINETIC_F3x_F2yz;
  abcd[106] = 2.0E0*I_KINETIC_H4xy_F2yz_a-3*I_KINETIC_F2xy_F2yz;
  abcd[107] = 2.0E0*I_KINETIC_H4xz_F2yz_a-3*I_KINETIC_F2xz_F2yz;
  abcd[108] = 2.0E0*I_KINETIC_H3x2y_F2yz_a-2*I_KINETIC_Fx2y_F2yz;
  abcd[109] = 2.0E0*I_KINETIC_H3xyz_F2yz_a-2*I_KINETIC_Fxyz_F2yz;
  abcd[110] = 2.0E0*I_KINETIC_H3x2z_F2yz_a-2*I_KINETIC_Fx2z_F2yz;
  abcd[111] = 2.0E0*I_KINETIC_H2x3y_F2yz_a-1*I_KINETIC_F3y_F2yz;
  abcd[112] = 2.0E0*I_KINETIC_H2x2yz_F2yz_a-1*I_KINETIC_F2yz_F2yz;
  abcd[113] = 2.0E0*I_KINETIC_H2xy2z_F2yz_a-1*I_KINETIC_Fy2z_F2yz;
  abcd[114] = 2.0E0*I_KINETIC_H2x3z_F2yz_a-1*I_KINETIC_F3z_F2yz;
  abcd[115] = 2.0E0*I_KINETIC_Hx4y_F2yz_a;
  abcd[116] = 2.0E0*I_KINETIC_Hx3yz_F2yz_a;
  abcd[117] = 2.0E0*I_KINETIC_Hx2y2z_F2yz_a;
  abcd[118] = 2.0E0*I_KINETIC_Hxy3z_F2yz_a;
  abcd[119] = 2.0E0*I_KINETIC_Hx4z_F2yz_a;
  abcd[120] = 2.0E0*I_KINETIC_H5x_Fy2z_a-4*I_KINETIC_F3x_Fy2z;
  abcd[121] = 2.0E0*I_KINETIC_H4xy_Fy2z_a-3*I_KINETIC_F2xy_Fy2z;
  abcd[122] = 2.0E0*I_KINETIC_H4xz_Fy2z_a-3*I_KINETIC_F2xz_Fy2z;
  abcd[123] = 2.0E0*I_KINETIC_H3x2y_Fy2z_a-2*I_KINETIC_Fx2y_Fy2z;
  abcd[124] = 2.0E0*I_KINETIC_H3xyz_Fy2z_a-2*I_KINETIC_Fxyz_Fy2z;
  abcd[125] = 2.0E0*I_KINETIC_H3x2z_Fy2z_a-2*I_KINETIC_Fx2z_Fy2z;
  abcd[126] = 2.0E0*I_KINETIC_H2x3y_Fy2z_a-1*I_KINETIC_F3y_Fy2z;
  abcd[127] = 2.0E0*I_KINETIC_H2x2yz_Fy2z_a-1*I_KINETIC_F2yz_Fy2z;
  abcd[128] = 2.0E0*I_KINETIC_H2xy2z_Fy2z_a-1*I_KINETIC_Fy2z_Fy2z;
  abcd[129] = 2.0E0*I_KINETIC_H2x3z_Fy2z_a-1*I_KINETIC_F3z_Fy2z;
  abcd[130] = 2.0E0*I_KINETIC_Hx4y_Fy2z_a;
  abcd[131] = 2.0E0*I_KINETIC_Hx3yz_Fy2z_a;
  abcd[132] = 2.0E0*I_KINETIC_Hx2y2z_Fy2z_a;
  abcd[133] = 2.0E0*I_KINETIC_Hxy3z_Fy2z_a;
  abcd[134] = 2.0E0*I_KINETIC_Hx4z_Fy2z_a;
  abcd[135] = 2.0E0*I_KINETIC_H5x_F3z_a-4*I_KINETIC_F3x_F3z;
  abcd[136] = 2.0E0*I_KINETIC_H4xy_F3z_a-3*I_KINETIC_F2xy_F3z;
  abcd[137] = 2.0E0*I_KINETIC_H4xz_F3z_a-3*I_KINETIC_F2xz_F3z;
  abcd[138] = 2.0E0*I_KINETIC_H3x2y_F3z_a-2*I_KINETIC_Fx2y_F3z;
  abcd[139] = 2.0E0*I_KINETIC_H3xyz_F3z_a-2*I_KINETIC_Fxyz_F3z;
  abcd[140] = 2.0E0*I_KINETIC_H3x2z_F3z_a-2*I_KINETIC_Fx2z_F3z;
  abcd[141] = 2.0E0*I_KINETIC_H2x3y_F3z_a-1*I_KINETIC_F3y_F3z;
  abcd[142] = 2.0E0*I_KINETIC_H2x2yz_F3z_a-1*I_KINETIC_F2yz_F3z;
  abcd[143] = 2.0E0*I_KINETIC_H2xy2z_F3z_a-1*I_KINETIC_Fy2z_F3z;
  abcd[144] = 2.0E0*I_KINETIC_H2x3z_F3z_a-1*I_KINETIC_F3z_F3z;
  abcd[145] = 2.0E0*I_KINETIC_Hx4y_F3z_a;
  abcd[146] = 2.0E0*I_KINETIC_Hx3yz_F3z_a;
  abcd[147] = 2.0E0*I_KINETIC_Hx2y2z_F3z_a;
  abcd[148] = 2.0E0*I_KINETIC_Hxy3z_F3z_a;
  abcd[149] = 2.0E0*I_KINETIC_Hx4z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_F_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_F_a
   * RHS shell quartet name: SQ_KINETIC_F_F
   ************************************************************/
  abcd[150] = 2.0E0*I_KINETIC_H4xy_F3x_a;
  abcd[151] = 2.0E0*I_KINETIC_H3x2y_F3x_a-1*I_KINETIC_F3x_F3x;
  abcd[152] = 2.0E0*I_KINETIC_H3xyz_F3x_a;
  abcd[153] = 2.0E0*I_KINETIC_H2x3y_F3x_a-2*I_KINETIC_F2xy_F3x;
  abcd[154] = 2.0E0*I_KINETIC_H2x2yz_F3x_a-1*I_KINETIC_F2xz_F3x;
  abcd[155] = 2.0E0*I_KINETIC_H2xy2z_F3x_a;
  abcd[156] = 2.0E0*I_KINETIC_Hx4y_F3x_a-3*I_KINETIC_Fx2y_F3x;
  abcd[157] = 2.0E0*I_KINETIC_Hx3yz_F3x_a-2*I_KINETIC_Fxyz_F3x;
  abcd[158] = 2.0E0*I_KINETIC_Hx2y2z_F3x_a-1*I_KINETIC_Fx2z_F3x;
  abcd[159] = 2.0E0*I_KINETIC_Hxy3z_F3x_a;
  abcd[160] = 2.0E0*I_KINETIC_H5y_F3x_a-4*I_KINETIC_F3y_F3x;
  abcd[161] = 2.0E0*I_KINETIC_H4yz_F3x_a-3*I_KINETIC_F2yz_F3x;
  abcd[162] = 2.0E0*I_KINETIC_H3y2z_F3x_a-2*I_KINETIC_Fy2z_F3x;
  abcd[163] = 2.0E0*I_KINETIC_H2y3z_F3x_a-1*I_KINETIC_F3z_F3x;
  abcd[164] = 2.0E0*I_KINETIC_Hy4z_F3x_a;
  abcd[165] = 2.0E0*I_KINETIC_H4xy_F2xy_a;
  abcd[166] = 2.0E0*I_KINETIC_H3x2y_F2xy_a-1*I_KINETIC_F3x_F2xy;
  abcd[167] = 2.0E0*I_KINETIC_H3xyz_F2xy_a;
  abcd[168] = 2.0E0*I_KINETIC_H2x3y_F2xy_a-2*I_KINETIC_F2xy_F2xy;
  abcd[169] = 2.0E0*I_KINETIC_H2x2yz_F2xy_a-1*I_KINETIC_F2xz_F2xy;
  abcd[170] = 2.0E0*I_KINETIC_H2xy2z_F2xy_a;
  abcd[171] = 2.0E0*I_KINETIC_Hx4y_F2xy_a-3*I_KINETIC_Fx2y_F2xy;
  abcd[172] = 2.0E0*I_KINETIC_Hx3yz_F2xy_a-2*I_KINETIC_Fxyz_F2xy;
  abcd[173] = 2.0E0*I_KINETIC_Hx2y2z_F2xy_a-1*I_KINETIC_Fx2z_F2xy;
  abcd[174] = 2.0E0*I_KINETIC_Hxy3z_F2xy_a;
  abcd[175] = 2.0E0*I_KINETIC_H5y_F2xy_a-4*I_KINETIC_F3y_F2xy;
  abcd[176] = 2.0E0*I_KINETIC_H4yz_F2xy_a-3*I_KINETIC_F2yz_F2xy;
  abcd[177] = 2.0E0*I_KINETIC_H3y2z_F2xy_a-2*I_KINETIC_Fy2z_F2xy;
  abcd[178] = 2.0E0*I_KINETIC_H2y3z_F2xy_a-1*I_KINETIC_F3z_F2xy;
  abcd[179] = 2.0E0*I_KINETIC_Hy4z_F2xy_a;
  abcd[180] = 2.0E0*I_KINETIC_H4xy_F2xz_a;
  abcd[181] = 2.0E0*I_KINETIC_H3x2y_F2xz_a-1*I_KINETIC_F3x_F2xz;
  abcd[182] = 2.0E0*I_KINETIC_H3xyz_F2xz_a;
  abcd[183] = 2.0E0*I_KINETIC_H2x3y_F2xz_a-2*I_KINETIC_F2xy_F2xz;
  abcd[184] = 2.0E0*I_KINETIC_H2x2yz_F2xz_a-1*I_KINETIC_F2xz_F2xz;
  abcd[185] = 2.0E0*I_KINETIC_H2xy2z_F2xz_a;
  abcd[186] = 2.0E0*I_KINETIC_Hx4y_F2xz_a-3*I_KINETIC_Fx2y_F2xz;
  abcd[187] = 2.0E0*I_KINETIC_Hx3yz_F2xz_a-2*I_KINETIC_Fxyz_F2xz;
  abcd[188] = 2.0E0*I_KINETIC_Hx2y2z_F2xz_a-1*I_KINETIC_Fx2z_F2xz;
  abcd[189] = 2.0E0*I_KINETIC_Hxy3z_F2xz_a;
  abcd[190] = 2.0E0*I_KINETIC_H5y_F2xz_a-4*I_KINETIC_F3y_F2xz;
  abcd[191] = 2.0E0*I_KINETIC_H4yz_F2xz_a-3*I_KINETIC_F2yz_F2xz;
  abcd[192] = 2.0E0*I_KINETIC_H3y2z_F2xz_a-2*I_KINETIC_Fy2z_F2xz;
  abcd[193] = 2.0E0*I_KINETIC_H2y3z_F2xz_a-1*I_KINETIC_F3z_F2xz;
  abcd[194] = 2.0E0*I_KINETIC_Hy4z_F2xz_a;
  abcd[195] = 2.0E0*I_KINETIC_H4xy_Fx2y_a;
  abcd[196] = 2.0E0*I_KINETIC_H3x2y_Fx2y_a-1*I_KINETIC_F3x_Fx2y;
  abcd[197] = 2.0E0*I_KINETIC_H3xyz_Fx2y_a;
  abcd[198] = 2.0E0*I_KINETIC_H2x3y_Fx2y_a-2*I_KINETIC_F2xy_Fx2y;
  abcd[199] = 2.0E0*I_KINETIC_H2x2yz_Fx2y_a-1*I_KINETIC_F2xz_Fx2y;
  abcd[200] = 2.0E0*I_KINETIC_H2xy2z_Fx2y_a;
  abcd[201] = 2.0E0*I_KINETIC_Hx4y_Fx2y_a-3*I_KINETIC_Fx2y_Fx2y;
  abcd[202] = 2.0E0*I_KINETIC_Hx3yz_Fx2y_a-2*I_KINETIC_Fxyz_Fx2y;
  abcd[203] = 2.0E0*I_KINETIC_Hx2y2z_Fx2y_a-1*I_KINETIC_Fx2z_Fx2y;
  abcd[204] = 2.0E0*I_KINETIC_Hxy3z_Fx2y_a;
  abcd[205] = 2.0E0*I_KINETIC_H5y_Fx2y_a-4*I_KINETIC_F3y_Fx2y;
  abcd[206] = 2.0E0*I_KINETIC_H4yz_Fx2y_a-3*I_KINETIC_F2yz_Fx2y;
  abcd[207] = 2.0E0*I_KINETIC_H3y2z_Fx2y_a-2*I_KINETIC_Fy2z_Fx2y;
  abcd[208] = 2.0E0*I_KINETIC_H2y3z_Fx2y_a-1*I_KINETIC_F3z_Fx2y;
  abcd[209] = 2.0E0*I_KINETIC_Hy4z_Fx2y_a;
  abcd[210] = 2.0E0*I_KINETIC_H4xy_Fxyz_a;
  abcd[211] = 2.0E0*I_KINETIC_H3x2y_Fxyz_a-1*I_KINETIC_F3x_Fxyz;
  abcd[212] = 2.0E0*I_KINETIC_H3xyz_Fxyz_a;
  abcd[213] = 2.0E0*I_KINETIC_H2x3y_Fxyz_a-2*I_KINETIC_F2xy_Fxyz;
  abcd[214] = 2.0E0*I_KINETIC_H2x2yz_Fxyz_a-1*I_KINETIC_F2xz_Fxyz;
  abcd[215] = 2.0E0*I_KINETIC_H2xy2z_Fxyz_a;
  abcd[216] = 2.0E0*I_KINETIC_Hx4y_Fxyz_a-3*I_KINETIC_Fx2y_Fxyz;
  abcd[217] = 2.0E0*I_KINETIC_Hx3yz_Fxyz_a-2*I_KINETIC_Fxyz_Fxyz;
  abcd[218] = 2.0E0*I_KINETIC_Hx2y2z_Fxyz_a-1*I_KINETIC_Fx2z_Fxyz;
  abcd[219] = 2.0E0*I_KINETIC_Hxy3z_Fxyz_a;
  abcd[220] = 2.0E0*I_KINETIC_H5y_Fxyz_a-4*I_KINETIC_F3y_Fxyz;
  abcd[221] = 2.0E0*I_KINETIC_H4yz_Fxyz_a-3*I_KINETIC_F2yz_Fxyz;
  abcd[222] = 2.0E0*I_KINETIC_H3y2z_Fxyz_a-2*I_KINETIC_Fy2z_Fxyz;
  abcd[223] = 2.0E0*I_KINETIC_H2y3z_Fxyz_a-1*I_KINETIC_F3z_Fxyz;
  abcd[224] = 2.0E0*I_KINETIC_Hy4z_Fxyz_a;
  abcd[225] = 2.0E0*I_KINETIC_H4xy_Fx2z_a;
  abcd[226] = 2.0E0*I_KINETIC_H3x2y_Fx2z_a-1*I_KINETIC_F3x_Fx2z;
  abcd[227] = 2.0E0*I_KINETIC_H3xyz_Fx2z_a;
  abcd[228] = 2.0E0*I_KINETIC_H2x3y_Fx2z_a-2*I_KINETIC_F2xy_Fx2z;
  abcd[229] = 2.0E0*I_KINETIC_H2x2yz_Fx2z_a-1*I_KINETIC_F2xz_Fx2z;
  abcd[230] = 2.0E0*I_KINETIC_H2xy2z_Fx2z_a;
  abcd[231] = 2.0E0*I_KINETIC_Hx4y_Fx2z_a-3*I_KINETIC_Fx2y_Fx2z;
  abcd[232] = 2.0E0*I_KINETIC_Hx3yz_Fx2z_a-2*I_KINETIC_Fxyz_Fx2z;
  abcd[233] = 2.0E0*I_KINETIC_Hx2y2z_Fx2z_a-1*I_KINETIC_Fx2z_Fx2z;
  abcd[234] = 2.0E0*I_KINETIC_Hxy3z_Fx2z_a;
  abcd[235] = 2.0E0*I_KINETIC_H5y_Fx2z_a-4*I_KINETIC_F3y_Fx2z;
  abcd[236] = 2.0E0*I_KINETIC_H4yz_Fx2z_a-3*I_KINETIC_F2yz_Fx2z;
  abcd[237] = 2.0E0*I_KINETIC_H3y2z_Fx2z_a-2*I_KINETIC_Fy2z_Fx2z;
  abcd[238] = 2.0E0*I_KINETIC_H2y3z_Fx2z_a-1*I_KINETIC_F3z_Fx2z;
  abcd[239] = 2.0E0*I_KINETIC_Hy4z_Fx2z_a;
  abcd[240] = 2.0E0*I_KINETIC_H4xy_F3y_a;
  abcd[241] = 2.0E0*I_KINETIC_H3x2y_F3y_a-1*I_KINETIC_F3x_F3y;
  abcd[242] = 2.0E0*I_KINETIC_H3xyz_F3y_a;
  abcd[243] = 2.0E0*I_KINETIC_H2x3y_F3y_a-2*I_KINETIC_F2xy_F3y;
  abcd[244] = 2.0E0*I_KINETIC_H2x2yz_F3y_a-1*I_KINETIC_F2xz_F3y;
  abcd[245] = 2.0E0*I_KINETIC_H2xy2z_F3y_a;
  abcd[246] = 2.0E0*I_KINETIC_Hx4y_F3y_a-3*I_KINETIC_Fx2y_F3y;
  abcd[247] = 2.0E0*I_KINETIC_Hx3yz_F3y_a-2*I_KINETIC_Fxyz_F3y;
  abcd[248] = 2.0E0*I_KINETIC_Hx2y2z_F3y_a-1*I_KINETIC_Fx2z_F3y;
  abcd[249] = 2.0E0*I_KINETIC_Hxy3z_F3y_a;
  abcd[250] = 2.0E0*I_KINETIC_H5y_F3y_a-4*I_KINETIC_F3y_F3y;
  abcd[251] = 2.0E0*I_KINETIC_H4yz_F3y_a-3*I_KINETIC_F2yz_F3y;
  abcd[252] = 2.0E0*I_KINETIC_H3y2z_F3y_a-2*I_KINETIC_Fy2z_F3y;
  abcd[253] = 2.0E0*I_KINETIC_H2y3z_F3y_a-1*I_KINETIC_F3z_F3y;
  abcd[254] = 2.0E0*I_KINETIC_Hy4z_F3y_a;
  abcd[255] = 2.0E0*I_KINETIC_H4xy_F2yz_a;
  abcd[256] = 2.0E0*I_KINETIC_H3x2y_F2yz_a-1*I_KINETIC_F3x_F2yz;
  abcd[257] = 2.0E0*I_KINETIC_H3xyz_F2yz_a;
  abcd[258] = 2.0E0*I_KINETIC_H2x3y_F2yz_a-2*I_KINETIC_F2xy_F2yz;
  abcd[259] = 2.0E0*I_KINETIC_H2x2yz_F2yz_a-1*I_KINETIC_F2xz_F2yz;
  abcd[260] = 2.0E0*I_KINETIC_H2xy2z_F2yz_a;
  abcd[261] = 2.0E0*I_KINETIC_Hx4y_F2yz_a-3*I_KINETIC_Fx2y_F2yz;
  abcd[262] = 2.0E0*I_KINETIC_Hx3yz_F2yz_a-2*I_KINETIC_Fxyz_F2yz;
  abcd[263] = 2.0E0*I_KINETIC_Hx2y2z_F2yz_a-1*I_KINETIC_Fx2z_F2yz;
  abcd[264] = 2.0E0*I_KINETIC_Hxy3z_F2yz_a;
  abcd[265] = 2.0E0*I_KINETIC_H5y_F2yz_a-4*I_KINETIC_F3y_F2yz;
  abcd[266] = 2.0E0*I_KINETIC_H4yz_F2yz_a-3*I_KINETIC_F2yz_F2yz;
  abcd[267] = 2.0E0*I_KINETIC_H3y2z_F2yz_a-2*I_KINETIC_Fy2z_F2yz;
  abcd[268] = 2.0E0*I_KINETIC_H2y3z_F2yz_a-1*I_KINETIC_F3z_F2yz;
  abcd[269] = 2.0E0*I_KINETIC_Hy4z_F2yz_a;
  abcd[270] = 2.0E0*I_KINETIC_H4xy_Fy2z_a;
  abcd[271] = 2.0E0*I_KINETIC_H3x2y_Fy2z_a-1*I_KINETIC_F3x_Fy2z;
  abcd[272] = 2.0E0*I_KINETIC_H3xyz_Fy2z_a;
  abcd[273] = 2.0E0*I_KINETIC_H2x3y_Fy2z_a-2*I_KINETIC_F2xy_Fy2z;
  abcd[274] = 2.0E0*I_KINETIC_H2x2yz_Fy2z_a-1*I_KINETIC_F2xz_Fy2z;
  abcd[275] = 2.0E0*I_KINETIC_H2xy2z_Fy2z_a;
  abcd[276] = 2.0E0*I_KINETIC_Hx4y_Fy2z_a-3*I_KINETIC_Fx2y_Fy2z;
  abcd[277] = 2.0E0*I_KINETIC_Hx3yz_Fy2z_a-2*I_KINETIC_Fxyz_Fy2z;
  abcd[278] = 2.0E0*I_KINETIC_Hx2y2z_Fy2z_a-1*I_KINETIC_Fx2z_Fy2z;
  abcd[279] = 2.0E0*I_KINETIC_Hxy3z_Fy2z_a;
  abcd[280] = 2.0E0*I_KINETIC_H5y_Fy2z_a-4*I_KINETIC_F3y_Fy2z;
  abcd[281] = 2.0E0*I_KINETIC_H4yz_Fy2z_a-3*I_KINETIC_F2yz_Fy2z;
  abcd[282] = 2.0E0*I_KINETIC_H3y2z_Fy2z_a-2*I_KINETIC_Fy2z_Fy2z;
  abcd[283] = 2.0E0*I_KINETIC_H2y3z_Fy2z_a-1*I_KINETIC_F3z_Fy2z;
  abcd[284] = 2.0E0*I_KINETIC_Hy4z_Fy2z_a;
  abcd[285] = 2.0E0*I_KINETIC_H4xy_F3z_a;
  abcd[286] = 2.0E0*I_KINETIC_H3x2y_F3z_a-1*I_KINETIC_F3x_F3z;
  abcd[287] = 2.0E0*I_KINETIC_H3xyz_F3z_a;
  abcd[288] = 2.0E0*I_KINETIC_H2x3y_F3z_a-2*I_KINETIC_F2xy_F3z;
  abcd[289] = 2.0E0*I_KINETIC_H2x2yz_F3z_a-1*I_KINETIC_F2xz_F3z;
  abcd[290] = 2.0E0*I_KINETIC_H2xy2z_F3z_a;
  abcd[291] = 2.0E0*I_KINETIC_Hx4y_F3z_a-3*I_KINETIC_Fx2y_F3z;
  abcd[292] = 2.0E0*I_KINETIC_Hx3yz_F3z_a-2*I_KINETIC_Fxyz_F3z;
  abcd[293] = 2.0E0*I_KINETIC_Hx2y2z_F3z_a-1*I_KINETIC_Fx2z_F3z;
  abcd[294] = 2.0E0*I_KINETIC_Hxy3z_F3z_a;
  abcd[295] = 2.0E0*I_KINETIC_H5y_F3z_a-4*I_KINETIC_F3y_F3z;
  abcd[296] = 2.0E0*I_KINETIC_H4yz_F3z_a-3*I_KINETIC_F2yz_F3z;
  abcd[297] = 2.0E0*I_KINETIC_H3y2z_F3z_a-2*I_KINETIC_Fy2z_F3z;
  abcd[298] = 2.0E0*I_KINETIC_H2y3z_F3z_a-1*I_KINETIC_F3z_F3z;
  abcd[299] = 2.0E0*I_KINETIC_Hy4z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_F_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_F_a
   * RHS shell quartet name: SQ_KINETIC_F_F
   ************************************************************/
  abcd[300] = 2.0E0*I_KINETIC_H4xz_F3x_a;
  abcd[301] = 2.0E0*I_KINETIC_H3xyz_F3x_a;
  abcd[302] = 2.0E0*I_KINETIC_H3x2z_F3x_a-1*I_KINETIC_F3x_F3x;
  abcd[303] = 2.0E0*I_KINETIC_H2x2yz_F3x_a;
  abcd[304] = 2.0E0*I_KINETIC_H2xy2z_F3x_a-1*I_KINETIC_F2xy_F3x;
  abcd[305] = 2.0E0*I_KINETIC_H2x3z_F3x_a-2*I_KINETIC_F2xz_F3x;
  abcd[306] = 2.0E0*I_KINETIC_Hx3yz_F3x_a;
  abcd[307] = 2.0E0*I_KINETIC_Hx2y2z_F3x_a-1*I_KINETIC_Fx2y_F3x;
  abcd[308] = 2.0E0*I_KINETIC_Hxy3z_F3x_a-2*I_KINETIC_Fxyz_F3x;
  abcd[309] = 2.0E0*I_KINETIC_Hx4z_F3x_a-3*I_KINETIC_Fx2z_F3x;
  abcd[310] = 2.0E0*I_KINETIC_H4yz_F3x_a;
  abcd[311] = 2.0E0*I_KINETIC_H3y2z_F3x_a-1*I_KINETIC_F3y_F3x;
  abcd[312] = 2.0E0*I_KINETIC_H2y3z_F3x_a-2*I_KINETIC_F2yz_F3x;
  abcd[313] = 2.0E0*I_KINETIC_Hy4z_F3x_a-3*I_KINETIC_Fy2z_F3x;
  abcd[314] = 2.0E0*I_KINETIC_H5z_F3x_a-4*I_KINETIC_F3z_F3x;
  abcd[315] = 2.0E0*I_KINETIC_H4xz_F2xy_a;
  abcd[316] = 2.0E0*I_KINETIC_H3xyz_F2xy_a;
  abcd[317] = 2.0E0*I_KINETIC_H3x2z_F2xy_a-1*I_KINETIC_F3x_F2xy;
  abcd[318] = 2.0E0*I_KINETIC_H2x2yz_F2xy_a;
  abcd[319] = 2.0E0*I_KINETIC_H2xy2z_F2xy_a-1*I_KINETIC_F2xy_F2xy;
  abcd[320] = 2.0E0*I_KINETIC_H2x3z_F2xy_a-2*I_KINETIC_F2xz_F2xy;
  abcd[321] = 2.0E0*I_KINETIC_Hx3yz_F2xy_a;
  abcd[322] = 2.0E0*I_KINETIC_Hx2y2z_F2xy_a-1*I_KINETIC_Fx2y_F2xy;
  abcd[323] = 2.0E0*I_KINETIC_Hxy3z_F2xy_a-2*I_KINETIC_Fxyz_F2xy;
  abcd[324] = 2.0E0*I_KINETIC_Hx4z_F2xy_a-3*I_KINETIC_Fx2z_F2xy;
  abcd[325] = 2.0E0*I_KINETIC_H4yz_F2xy_a;
  abcd[326] = 2.0E0*I_KINETIC_H3y2z_F2xy_a-1*I_KINETIC_F3y_F2xy;
  abcd[327] = 2.0E0*I_KINETIC_H2y3z_F2xy_a-2*I_KINETIC_F2yz_F2xy;
  abcd[328] = 2.0E0*I_KINETIC_Hy4z_F2xy_a-3*I_KINETIC_Fy2z_F2xy;
  abcd[329] = 2.0E0*I_KINETIC_H5z_F2xy_a-4*I_KINETIC_F3z_F2xy;
  abcd[330] = 2.0E0*I_KINETIC_H4xz_F2xz_a;
  abcd[331] = 2.0E0*I_KINETIC_H3xyz_F2xz_a;
  abcd[332] = 2.0E0*I_KINETIC_H3x2z_F2xz_a-1*I_KINETIC_F3x_F2xz;
  abcd[333] = 2.0E0*I_KINETIC_H2x2yz_F2xz_a;
  abcd[334] = 2.0E0*I_KINETIC_H2xy2z_F2xz_a-1*I_KINETIC_F2xy_F2xz;
  abcd[335] = 2.0E0*I_KINETIC_H2x3z_F2xz_a-2*I_KINETIC_F2xz_F2xz;
  abcd[336] = 2.0E0*I_KINETIC_Hx3yz_F2xz_a;
  abcd[337] = 2.0E0*I_KINETIC_Hx2y2z_F2xz_a-1*I_KINETIC_Fx2y_F2xz;
  abcd[338] = 2.0E0*I_KINETIC_Hxy3z_F2xz_a-2*I_KINETIC_Fxyz_F2xz;
  abcd[339] = 2.0E0*I_KINETIC_Hx4z_F2xz_a-3*I_KINETIC_Fx2z_F2xz;
  abcd[340] = 2.0E0*I_KINETIC_H4yz_F2xz_a;
  abcd[341] = 2.0E0*I_KINETIC_H3y2z_F2xz_a-1*I_KINETIC_F3y_F2xz;
  abcd[342] = 2.0E0*I_KINETIC_H2y3z_F2xz_a-2*I_KINETIC_F2yz_F2xz;
  abcd[343] = 2.0E0*I_KINETIC_Hy4z_F2xz_a-3*I_KINETIC_Fy2z_F2xz;
  abcd[344] = 2.0E0*I_KINETIC_H5z_F2xz_a-4*I_KINETIC_F3z_F2xz;
  abcd[345] = 2.0E0*I_KINETIC_H4xz_Fx2y_a;
  abcd[346] = 2.0E0*I_KINETIC_H3xyz_Fx2y_a;
  abcd[347] = 2.0E0*I_KINETIC_H3x2z_Fx2y_a-1*I_KINETIC_F3x_Fx2y;
  abcd[348] = 2.0E0*I_KINETIC_H2x2yz_Fx2y_a;
  abcd[349] = 2.0E0*I_KINETIC_H2xy2z_Fx2y_a-1*I_KINETIC_F2xy_Fx2y;
  abcd[350] = 2.0E0*I_KINETIC_H2x3z_Fx2y_a-2*I_KINETIC_F2xz_Fx2y;
  abcd[351] = 2.0E0*I_KINETIC_Hx3yz_Fx2y_a;
  abcd[352] = 2.0E0*I_KINETIC_Hx2y2z_Fx2y_a-1*I_KINETIC_Fx2y_Fx2y;
  abcd[353] = 2.0E0*I_KINETIC_Hxy3z_Fx2y_a-2*I_KINETIC_Fxyz_Fx2y;
  abcd[354] = 2.0E0*I_KINETIC_Hx4z_Fx2y_a-3*I_KINETIC_Fx2z_Fx2y;
  abcd[355] = 2.0E0*I_KINETIC_H4yz_Fx2y_a;
  abcd[356] = 2.0E0*I_KINETIC_H3y2z_Fx2y_a-1*I_KINETIC_F3y_Fx2y;
  abcd[357] = 2.0E0*I_KINETIC_H2y3z_Fx2y_a-2*I_KINETIC_F2yz_Fx2y;
  abcd[358] = 2.0E0*I_KINETIC_Hy4z_Fx2y_a-3*I_KINETIC_Fy2z_Fx2y;
  abcd[359] = 2.0E0*I_KINETIC_H5z_Fx2y_a-4*I_KINETIC_F3z_Fx2y;
  abcd[360] = 2.0E0*I_KINETIC_H4xz_Fxyz_a;
  abcd[361] = 2.0E0*I_KINETIC_H3xyz_Fxyz_a;
  abcd[362] = 2.0E0*I_KINETIC_H3x2z_Fxyz_a-1*I_KINETIC_F3x_Fxyz;
  abcd[363] = 2.0E0*I_KINETIC_H2x2yz_Fxyz_a;
  abcd[364] = 2.0E0*I_KINETIC_H2xy2z_Fxyz_a-1*I_KINETIC_F2xy_Fxyz;
  abcd[365] = 2.0E0*I_KINETIC_H2x3z_Fxyz_a-2*I_KINETIC_F2xz_Fxyz;
  abcd[366] = 2.0E0*I_KINETIC_Hx3yz_Fxyz_a;
  abcd[367] = 2.0E0*I_KINETIC_Hx2y2z_Fxyz_a-1*I_KINETIC_Fx2y_Fxyz;
  abcd[368] = 2.0E0*I_KINETIC_Hxy3z_Fxyz_a-2*I_KINETIC_Fxyz_Fxyz;
  abcd[369] = 2.0E0*I_KINETIC_Hx4z_Fxyz_a-3*I_KINETIC_Fx2z_Fxyz;
  abcd[370] = 2.0E0*I_KINETIC_H4yz_Fxyz_a;
  abcd[371] = 2.0E0*I_KINETIC_H3y2z_Fxyz_a-1*I_KINETIC_F3y_Fxyz;
  abcd[372] = 2.0E0*I_KINETIC_H2y3z_Fxyz_a-2*I_KINETIC_F2yz_Fxyz;
  abcd[373] = 2.0E0*I_KINETIC_Hy4z_Fxyz_a-3*I_KINETIC_Fy2z_Fxyz;
  abcd[374] = 2.0E0*I_KINETIC_H5z_Fxyz_a-4*I_KINETIC_F3z_Fxyz;
  abcd[375] = 2.0E0*I_KINETIC_H4xz_Fx2z_a;
  abcd[376] = 2.0E0*I_KINETIC_H3xyz_Fx2z_a;
  abcd[377] = 2.0E0*I_KINETIC_H3x2z_Fx2z_a-1*I_KINETIC_F3x_Fx2z;
  abcd[378] = 2.0E0*I_KINETIC_H2x2yz_Fx2z_a;
  abcd[379] = 2.0E0*I_KINETIC_H2xy2z_Fx2z_a-1*I_KINETIC_F2xy_Fx2z;
  abcd[380] = 2.0E0*I_KINETIC_H2x3z_Fx2z_a-2*I_KINETIC_F2xz_Fx2z;
  abcd[381] = 2.0E0*I_KINETIC_Hx3yz_Fx2z_a;
  abcd[382] = 2.0E0*I_KINETIC_Hx2y2z_Fx2z_a-1*I_KINETIC_Fx2y_Fx2z;
  abcd[383] = 2.0E0*I_KINETIC_Hxy3z_Fx2z_a-2*I_KINETIC_Fxyz_Fx2z;
  abcd[384] = 2.0E0*I_KINETIC_Hx4z_Fx2z_a-3*I_KINETIC_Fx2z_Fx2z;
  abcd[385] = 2.0E0*I_KINETIC_H4yz_Fx2z_a;
  abcd[386] = 2.0E0*I_KINETIC_H3y2z_Fx2z_a-1*I_KINETIC_F3y_Fx2z;
  abcd[387] = 2.0E0*I_KINETIC_H2y3z_Fx2z_a-2*I_KINETIC_F2yz_Fx2z;
  abcd[388] = 2.0E0*I_KINETIC_Hy4z_Fx2z_a-3*I_KINETIC_Fy2z_Fx2z;
  abcd[389] = 2.0E0*I_KINETIC_H5z_Fx2z_a-4*I_KINETIC_F3z_Fx2z;
  abcd[390] = 2.0E0*I_KINETIC_H4xz_F3y_a;
  abcd[391] = 2.0E0*I_KINETIC_H3xyz_F3y_a;
  abcd[392] = 2.0E0*I_KINETIC_H3x2z_F3y_a-1*I_KINETIC_F3x_F3y;
  abcd[393] = 2.0E0*I_KINETIC_H2x2yz_F3y_a;
  abcd[394] = 2.0E0*I_KINETIC_H2xy2z_F3y_a-1*I_KINETIC_F2xy_F3y;
  abcd[395] = 2.0E0*I_KINETIC_H2x3z_F3y_a-2*I_KINETIC_F2xz_F3y;
  abcd[396] = 2.0E0*I_KINETIC_Hx3yz_F3y_a;
  abcd[397] = 2.0E0*I_KINETIC_Hx2y2z_F3y_a-1*I_KINETIC_Fx2y_F3y;
  abcd[398] = 2.0E0*I_KINETIC_Hxy3z_F3y_a-2*I_KINETIC_Fxyz_F3y;
  abcd[399] = 2.0E0*I_KINETIC_Hx4z_F3y_a-3*I_KINETIC_Fx2z_F3y;
  abcd[400] = 2.0E0*I_KINETIC_H4yz_F3y_a;
  abcd[401] = 2.0E0*I_KINETIC_H3y2z_F3y_a-1*I_KINETIC_F3y_F3y;
  abcd[402] = 2.0E0*I_KINETIC_H2y3z_F3y_a-2*I_KINETIC_F2yz_F3y;
  abcd[403] = 2.0E0*I_KINETIC_Hy4z_F3y_a-3*I_KINETIC_Fy2z_F3y;
  abcd[404] = 2.0E0*I_KINETIC_H5z_F3y_a-4*I_KINETIC_F3z_F3y;
  abcd[405] = 2.0E0*I_KINETIC_H4xz_F2yz_a;
  abcd[406] = 2.0E0*I_KINETIC_H3xyz_F2yz_a;
  abcd[407] = 2.0E0*I_KINETIC_H3x2z_F2yz_a-1*I_KINETIC_F3x_F2yz;
  abcd[408] = 2.0E0*I_KINETIC_H2x2yz_F2yz_a;
  abcd[409] = 2.0E0*I_KINETIC_H2xy2z_F2yz_a-1*I_KINETIC_F2xy_F2yz;
  abcd[410] = 2.0E0*I_KINETIC_H2x3z_F2yz_a-2*I_KINETIC_F2xz_F2yz;
  abcd[411] = 2.0E0*I_KINETIC_Hx3yz_F2yz_a;
  abcd[412] = 2.0E0*I_KINETIC_Hx2y2z_F2yz_a-1*I_KINETIC_Fx2y_F2yz;
  abcd[413] = 2.0E0*I_KINETIC_Hxy3z_F2yz_a-2*I_KINETIC_Fxyz_F2yz;
  abcd[414] = 2.0E0*I_KINETIC_Hx4z_F2yz_a-3*I_KINETIC_Fx2z_F2yz;
  abcd[415] = 2.0E0*I_KINETIC_H4yz_F2yz_a;
  abcd[416] = 2.0E0*I_KINETIC_H3y2z_F2yz_a-1*I_KINETIC_F3y_F2yz;
  abcd[417] = 2.0E0*I_KINETIC_H2y3z_F2yz_a-2*I_KINETIC_F2yz_F2yz;
  abcd[418] = 2.0E0*I_KINETIC_Hy4z_F2yz_a-3*I_KINETIC_Fy2z_F2yz;
  abcd[419] = 2.0E0*I_KINETIC_H5z_F2yz_a-4*I_KINETIC_F3z_F2yz;
  abcd[420] = 2.0E0*I_KINETIC_H4xz_Fy2z_a;
  abcd[421] = 2.0E0*I_KINETIC_H3xyz_Fy2z_a;
  abcd[422] = 2.0E0*I_KINETIC_H3x2z_Fy2z_a-1*I_KINETIC_F3x_Fy2z;
  abcd[423] = 2.0E0*I_KINETIC_H2x2yz_Fy2z_a;
  abcd[424] = 2.0E0*I_KINETIC_H2xy2z_Fy2z_a-1*I_KINETIC_F2xy_Fy2z;
  abcd[425] = 2.0E0*I_KINETIC_H2x3z_Fy2z_a-2*I_KINETIC_F2xz_Fy2z;
  abcd[426] = 2.0E0*I_KINETIC_Hx3yz_Fy2z_a;
  abcd[427] = 2.0E0*I_KINETIC_Hx2y2z_Fy2z_a-1*I_KINETIC_Fx2y_Fy2z;
  abcd[428] = 2.0E0*I_KINETIC_Hxy3z_Fy2z_a-2*I_KINETIC_Fxyz_Fy2z;
  abcd[429] = 2.0E0*I_KINETIC_Hx4z_Fy2z_a-3*I_KINETIC_Fx2z_Fy2z;
  abcd[430] = 2.0E0*I_KINETIC_H4yz_Fy2z_a;
  abcd[431] = 2.0E0*I_KINETIC_H3y2z_Fy2z_a-1*I_KINETIC_F3y_Fy2z;
  abcd[432] = 2.0E0*I_KINETIC_H2y3z_Fy2z_a-2*I_KINETIC_F2yz_Fy2z;
  abcd[433] = 2.0E0*I_KINETIC_Hy4z_Fy2z_a-3*I_KINETIC_Fy2z_Fy2z;
  abcd[434] = 2.0E0*I_KINETIC_H5z_Fy2z_a-4*I_KINETIC_F3z_Fy2z;
  abcd[435] = 2.0E0*I_KINETIC_H4xz_F3z_a;
  abcd[436] = 2.0E0*I_KINETIC_H3xyz_F3z_a;
  abcd[437] = 2.0E0*I_KINETIC_H3x2z_F3z_a-1*I_KINETIC_F3x_F3z;
  abcd[438] = 2.0E0*I_KINETIC_H2x2yz_F3z_a;
  abcd[439] = 2.0E0*I_KINETIC_H2xy2z_F3z_a-1*I_KINETIC_F2xy_F3z;
  abcd[440] = 2.0E0*I_KINETIC_H2x3z_F3z_a-2*I_KINETIC_F2xz_F3z;
  abcd[441] = 2.0E0*I_KINETIC_Hx3yz_F3z_a;
  abcd[442] = 2.0E0*I_KINETIC_Hx2y2z_F3z_a-1*I_KINETIC_Fx2y_F3z;
  abcd[443] = 2.0E0*I_KINETIC_Hxy3z_F3z_a-2*I_KINETIC_Fxyz_F3z;
  abcd[444] = 2.0E0*I_KINETIC_Hx4z_F3z_a-3*I_KINETIC_Fx2z_F3z;
  abcd[445] = 2.0E0*I_KINETIC_H4yz_F3z_a;
  abcd[446] = 2.0E0*I_KINETIC_H3y2z_F3z_a-1*I_KINETIC_F3y_F3z;
  abcd[447] = 2.0E0*I_KINETIC_H2y3z_F3z_a-2*I_KINETIC_F2yz_F3z;
  abcd[448] = 2.0E0*I_KINETIC_Hy4z_F3z_a-3*I_KINETIC_Fy2z_F3z;
  abcd[449] = 2.0E0*I_KINETIC_H5z_F3z_a-4*I_KINETIC_F3z_F3z;
}
