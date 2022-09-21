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
// BRA1 as redundant position, total RHS integrals evaluated as: 6601
// BRA2 as redundant position, total RHS integrals evaluated as: 6368
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

void hgp_os_kinetic_h_f_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_I6x_F3x_a = 0.0E0;
  Double I_KINETIC_I5xy_F3x_a = 0.0E0;
  Double I_KINETIC_I5xz_F3x_a = 0.0E0;
  Double I_KINETIC_I4x2y_F3x_a = 0.0E0;
  Double I_KINETIC_I4xyz_F3x_a = 0.0E0;
  Double I_KINETIC_I4x2z_F3x_a = 0.0E0;
  Double I_KINETIC_I3x3y_F3x_a = 0.0E0;
  Double I_KINETIC_I3x2yz_F3x_a = 0.0E0;
  Double I_KINETIC_I3xy2z_F3x_a = 0.0E0;
  Double I_KINETIC_I3x3z_F3x_a = 0.0E0;
  Double I_KINETIC_I2x4y_F3x_a = 0.0E0;
  Double I_KINETIC_I2x3yz_F3x_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_F3x_a = 0.0E0;
  Double I_KINETIC_I2xy3z_F3x_a = 0.0E0;
  Double I_KINETIC_I2x4z_F3x_a = 0.0E0;
  Double I_KINETIC_Ix5y_F3x_a = 0.0E0;
  Double I_KINETIC_Ix4yz_F3x_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_F3x_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_F3x_a = 0.0E0;
  Double I_KINETIC_Ixy4z_F3x_a = 0.0E0;
  Double I_KINETIC_Ix5z_F3x_a = 0.0E0;
  Double I_KINETIC_I6y_F3x_a = 0.0E0;
  Double I_KINETIC_I5yz_F3x_a = 0.0E0;
  Double I_KINETIC_I4y2z_F3x_a = 0.0E0;
  Double I_KINETIC_I3y3z_F3x_a = 0.0E0;
  Double I_KINETIC_I2y4z_F3x_a = 0.0E0;
  Double I_KINETIC_Iy5z_F3x_a = 0.0E0;
  Double I_KINETIC_I6z_F3x_a = 0.0E0;
  Double I_KINETIC_I6x_F2xy_a = 0.0E0;
  Double I_KINETIC_I5xy_F2xy_a = 0.0E0;
  Double I_KINETIC_I5xz_F2xy_a = 0.0E0;
  Double I_KINETIC_I4x2y_F2xy_a = 0.0E0;
  Double I_KINETIC_I4xyz_F2xy_a = 0.0E0;
  Double I_KINETIC_I4x2z_F2xy_a = 0.0E0;
  Double I_KINETIC_I3x3y_F2xy_a = 0.0E0;
  Double I_KINETIC_I3x2yz_F2xy_a = 0.0E0;
  Double I_KINETIC_I3xy2z_F2xy_a = 0.0E0;
  Double I_KINETIC_I3x3z_F2xy_a = 0.0E0;
  Double I_KINETIC_I2x4y_F2xy_a = 0.0E0;
  Double I_KINETIC_I2x3yz_F2xy_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_F2xy_a = 0.0E0;
  Double I_KINETIC_I2xy3z_F2xy_a = 0.0E0;
  Double I_KINETIC_I2x4z_F2xy_a = 0.0E0;
  Double I_KINETIC_Ix5y_F2xy_a = 0.0E0;
  Double I_KINETIC_Ix4yz_F2xy_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_F2xy_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_F2xy_a = 0.0E0;
  Double I_KINETIC_Ixy4z_F2xy_a = 0.0E0;
  Double I_KINETIC_Ix5z_F2xy_a = 0.0E0;
  Double I_KINETIC_I6y_F2xy_a = 0.0E0;
  Double I_KINETIC_I5yz_F2xy_a = 0.0E0;
  Double I_KINETIC_I4y2z_F2xy_a = 0.0E0;
  Double I_KINETIC_I3y3z_F2xy_a = 0.0E0;
  Double I_KINETIC_I2y4z_F2xy_a = 0.0E0;
  Double I_KINETIC_Iy5z_F2xy_a = 0.0E0;
  Double I_KINETIC_I6z_F2xy_a = 0.0E0;
  Double I_KINETIC_I6x_F2xz_a = 0.0E0;
  Double I_KINETIC_I5xy_F2xz_a = 0.0E0;
  Double I_KINETIC_I5xz_F2xz_a = 0.0E0;
  Double I_KINETIC_I4x2y_F2xz_a = 0.0E0;
  Double I_KINETIC_I4xyz_F2xz_a = 0.0E0;
  Double I_KINETIC_I4x2z_F2xz_a = 0.0E0;
  Double I_KINETIC_I3x3y_F2xz_a = 0.0E0;
  Double I_KINETIC_I3x2yz_F2xz_a = 0.0E0;
  Double I_KINETIC_I3xy2z_F2xz_a = 0.0E0;
  Double I_KINETIC_I3x3z_F2xz_a = 0.0E0;
  Double I_KINETIC_I2x4y_F2xz_a = 0.0E0;
  Double I_KINETIC_I2x3yz_F2xz_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_F2xz_a = 0.0E0;
  Double I_KINETIC_I2xy3z_F2xz_a = 0.0E0;
  Double I_KINETIC_I2x4z_F2xz_a = 0.0E0;
  Double I_KINETIC_Ix5y_F2xz_a = 0.0E0;
  Double I_KINETIC_Ix4yz_F2xz_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_F2xz_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_F2xz_a = 0.0E0;
  Double I_KINETIC_Ixy4z_F2xz_a = 0.0E0;
  Double I_KINETIC_Ix5z_F2xz_a = 0.0E0;
  Double I_KINETIC_I6y_F2xz_a = 0.0E0;
  Double I_KINETIC_I5yz_F2xz_a = 0.0E0;
  Double I_KINETIC_I4y2z_F2xz_a = 0.0E0;
  Double I_KINETIC_I3y3z_F2xz_a = 0.0E0;
  Double I_KINETIC_I2y4z_F2xz_a = 0.0E0;
  Double I_KINETIC_Iy5z_F2xz_a = 0.0E0;
  Double I_KINETIC_I6z_F2xz_a = 0.0E0;
  Double I_KINETIC_I6x_Fx2y_a = 0.0E0;
  Double I_KINETIC_I5xy_Fx2y_a = 0.0E0;
  Double I_KINETIC_I5xz_Fx2y_a = 0.0E0;
  Double I_KINETIC_I4x2y_Fx2y_a = 0.0E0;
  Double I_KINETIC_I4xyz_Fx2y_a = 0.0E0;
  Double I_KINETIC_I4x2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I3x3y_Fx2y_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I3x3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I2x4y_Fx2y_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I2x4z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Ix5y_Fx2y_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Ix5z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I6y_Fx2y_a = 0.0E0;
  Double I_KINETIC_I5yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_I4y2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I3y3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I2y4z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Iy5z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I6z_Fx2y_a = 0.0E0;
  Double I_KINETIC_I6x_Fxyz_a = 0.0E0;
  Double I_KINETIC_I5xy_Fxyz_a = 0.0E0;
  Double I_KINETIC_I5xz_Fxyz_a = 0.0E0;
  Double I_KINETIC_I4x2y_Fxyz_a = 0.0E0;
  Double I_KINETIC_I4xyz_Fxyz_a = 0.0E0;
  Double I_KINETIC_I4x2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I3x3y_Fxyz_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I3x3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I2x4y_Fxyz_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I2x4z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Ix5y_Fxyz_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Ix5z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I6y_Fxyz_a = 0.0E0;
  Double I_KINETIC_I5yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_I4y2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I3y3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I2y4z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Iy5z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I6z_Fxyz_a = 0.0E0;
  Double I_KINETIC_I6x_Fx2z_a = 0.0E0;
  Double I_KINETIC_I5xy_Fx2z_a = 0.0E0;
  Double I_KINETIC_I5xz_Fx2z_a = 0.0E0;
  Double I_KINETIC_I4x2y_Fx2z_a = 0.0E0;
  Double I_KINETIC_I4xyz_Fx2z_a = 0.0E0;
  Double I_KINETIC_I4x2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I3x3y_Fx2z_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I3x3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I2x4y_Fx2z_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I2x4z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Ix5y_Fx2z_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Ix5z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I6y_Fx2z_a = 0.0E0;
  Double I_KINETIC_I5yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_I4y2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I3y3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I2y4z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Iy5z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I6z_Fx2z_a = 0.0E0;
  Double I_KINETIC_I6x_F3y_a = 0.0E0;
  Double I_KINETIC_I5xy_F3y_a = 0.0E0;
  Double I_KINETIC_I5xz_F3y_a = 0.0E0;
  Double I_KINETIC_I4x2y_F3y_a = 0.0E0;
  Double I_KINETIC_I4xyz_F3y_a = 0.0E0;
  Double I_KINETIC_I4x2z_F3y_a = 0.0E0;
  Double I_KINETIC_I3x3y_F3y_a = 0.0E0;
  Double I_KINETIC_I3x2yz_F3y_a = 0.0E0;
  Double I_KINETIC_I3xy2z_F3y_a = 0.0E0;
  Double I_KINETIC_I3x3z_F3y_a = 0.0E0;
  Double I_KINETIC_I2x4y_F3y_a = 0.0E0;
  Double I_KINETIC_I2x3yz_F3y_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_F3y_a = 0.0E0;
  Double I_KINETIC_I2xy3z_F3y_a = 0.0E0;
  Double I_KINETIC_I2x4z_F3y_a = 0.0E0;
  Double I_KINETIC_Ix5y_F3y_a = 0.0E0;
  Double I_KINETIC_Ix4yz_F3y_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_F3y_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_F3y_a = 0.0E0;
  Double I_KINETIC_Ixy4z_F3y_a = 0.0E0;
  Double I_KINETIC_Ix5z_F3y_a = 0.0E0;
  Double I_KINETIC_I6y_F3y_a = 0.0E0;
  Double I_KINETIC_I5yz_F3y_a = 0.0E0;
  Double I_KINETIC_I4y2z_F3y_a = 0.0E0;
  Double I_KINETIC_I3y3z_F3y_a = 0.0E0;
  Double I_KINETIC_I2y4z_F3y_a = 0.0E0;
  Double I_KINETIC_Iy5z_F3y_a = 0.0E0;
  Double I_KINETIC_I6z_F3y_a = 0.0E0;
  Double I_KINETIC_I6x_F2yz_a = 0.0E0;
  Double I_KINETIC_I5xy_F2yz_a = 0.0E0;
  Double I_KINETIC_I5xz_F2yz_a = 0.0E0;
  Double I_KINETIC_I4x2y_F2yz_a = 0.0E0;
  Double I_KINETIC_I4xyz_F2yz_a = 0.0E0;
  Double I_KINETIC_I4x2z_F2yz_a = 0.0E0;
  Double I_KINETIC_I3x3y_F2yz_a = 0.0E0;
  Double I_KINETIC_I3x2yz_F2yz_a = 0.0E0;
  Double I_KINETIC_I3xy2z_F2yz_a = 0.0E0;
  Double I_KINETIC_I3x3z_F2yz_a = 0.0E0;
  Double I_KINETIC_I2x4y_F2yz_a = 0.0E0;
  Double I_KINETIC_I2x3yz_F2yz_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_F2yz_a = 0.0E0;
  Double I_KINETIC_I2xy3z_F2yz_a = 0.0E0;
  Double I_KINETIC_I2x4z_F2yz_a = 0.0E0;
  Double I_KINETIC_Ix5y_F2yz_a = 0.0E0;
  Double I_KINETIC_Ix4yz_F2yz_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_F2yz_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_F2yz_a = 0.0E0;
  Double I_KINETIC_Ixy4z_F2yz_a = 0.0E0;
  Double I_KINETIC_Ix5z_F2yz_a = 0.0E0;
  Double I_KINETIC_I6y_F2yz_a = 0.0E0;
  Double I_KINETIC_I5yz_F2yz_a = 0.0E0;
  Double I_KINETIC_I4y2z_F2yz_a = 0.0E0;
  Double I_KINETIC_I3y3z_F2yz_a = 0.0E0;
  Double I_KINETIC_I2y4z_F2yz_a = 0.0E0;
  Double I_KINETIC_Iy5z_F2yz_a = 0.0E0;
  Double I_KINETIC_I6z_F2yz_a = 0.0E0;
  Double I_KINETIC_I6x_Fy2z_a = 0.0E0;
  Double I_KINETIC_I5xy_Fy2z_a = 0.0E0;
  Double I_KINETIC_I5xz_Fy2z_a = 0.0E0;
  Double I_KINETIC_I4x2y_Fy2z_a = 0.0E0;
  Double I_KINETIC_I4xyz_Fy2z_a = 0.0E0;
  Double I_KINETIC_I4x2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I3x3y_Fy2z_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I3x3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I2x4y_Fy2z_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I2x4z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Ix5y_Fy2z_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Ix5z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I6y_Fy2z_a = 0.0E0;
  Double I_KINETIC_I5yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_I4y2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I3y3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I2y4z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Iy5z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I6z_Fy2z_a = 0.0E0;
  Double I_KINETIC_I6x_F3z_a = 0.0E0;
  Double I_KINETIC_I5xy_F3z_a = 0.0E0;
  Double I_KINETIC_I5xz_F3z_a = 0.0E0;
  Double I_KINETIC_I4x2y_F3z_a = 0.0E0;
  Double I_KINETIC_I4xyz_F3z_a = 0.0E0;
  Double I_KINETIC_I4x2z_F3z_a = 0.0E0;
  Double I_KINETIC_I3x3y_F3z_a = 0.0E0;
  Double I_KINETIC_I3x2yz_F3z_a = 0.0E0;
  Double I_KINETIC_I3xy2z_F3z_a = 0.0E0;
  Double I_KINETIC_I3x3z_F3z_a = 0.0E0;
  Double I_KINETIC_I2x4y_F3z_a = 0.0E0;
  Double I_KINETIC_I2x3yz_F3z_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_F3z_a = 0.0E0;
  Double I_KINETIC_I2xy3z_F3z_a = 0.0E0;
  Double I_KINETIC_I2x4z_F3z_a = 0.0E0;
  Double I_KINETIC_Ix5y_F3z_a = 0.0E0;
  Double I_KINETIC_Ix4yz_F3z_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_F3z_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_F3z_a = 0.0E0;
  Double I_KINETIC_Ixy4z_F3z_a = 0.0E0;
  Double I_KINETIC_Ix5z_F3z_a = 0.0E0;
  Double I_KINETIC_I6y_F3z_a = 0.0E0;
  Double I_KINETIC_I5yz_F3z_a = 0.0E0;
  Double I_KINETIC_I4y2z_F3z_a = 0.0E0;
  Double I_KINETIC_I3y3z_F3z_a = 0.0E0;
  Double I_KINETIC_I2y4z_F3z_a = 0.0E0;
  Double I_KINETIC_Iy5z_F3z_a = 0.0E0;
  Double I_KINETIC_I6z_F3z_a = 0.0E0;
  Double I_KINETIC_G4x_F3x = 0.0E0;
  Double I_KINETIC_G3xy_F3x = 0.0E0;
  Double I_KINETIC_G3xz_F3x = 0.0E0;
  Double I_KINETIC_G2x2y_F3x = 0.0E0;
  Double I_KINETIC_G2xyz_F3x = 0.0E0;
  Double I_KINETIC_G2x2z_F3x = 0.0E0;
  Double I_KINETIC_Gx3y_F3x = 0.0E0;
  Double I_KINETIC_Gx2yz_F3x = 0.0E0;
  Double I_KINETIC_Gxy2z_F3x = 0.0E0;
  Double I_KINETIC_Gx3z_F3x = 0.0E0;
  Double I_KINETIC_G4y_F3x = 0.0E0;
  Double I_KINETIC_G3yz_F3x = 0.0E0;
  Double I_KINETIC_G2y2z_F3x = 0.0E0;
  Double I_KINETIC_Gy3z_F3x = 0.0E0;
  Double I_KINETIC_G4z_F3x = 0.0E0;
  Double I_KINETIC_G4x_F2xy = 0.0E0;
  Double I_KINETIC_G3xy_F2xy = 0.0E0;
  Double I_KINETIC_G3xz_F2xy = 0.0E0;
  Double I_KINETIC_G2x2y_F2xy = 0.0E0;
  Double I_KINETIC_G2xyz_F2xy = 0.0E0;
  Double I_KINETIC_G2x2z_F2xy = 0.0E0;
  Double I_KINETIC_Gx3y_F2xy = 0.0E0;
  Double I_KINETIC_Gx2yz_F2xy = 0.0E0;
  Double I_KINETIC_Gxy2z_F2xy = 0.0E0;
  Double I_KINETIC_Gx3z_F2xy = 0.0E0;
  Double I_KINETIC_G4y_F2xy = 0.0E0;
  Double I_KINETIC_G3yz_F2xy = 0.0E0;
  Double I_KINETIC_G2y2z_F2xy = 0.0E0;
  Double I_KINETIC_Gy3z_F2xy = 0.0E0;
  Double I_KINETIC_G4z_F2xy = 0.0E0;
  Double I_KINETIC_G4x_F2xz = 0.0E0;
  Double I_KINETIC_G3xy_F2xz = 0.0E0;
  Double I_KINETIC_G3xz_F2xz = 0.0E0;
  Double I_KINETIC_G2x2y_F2xz = 0.0E0;
  Double I_KINETIC_G2xyz_F2xz = 0.0E0;
  Double I_KINETIC_G2x2z_F2xz = 0.0E0;
  Double I_KINETIC_Gx3y_F2xz = 0.0E0;
  Double I_KINETIC_Gx2yz_F2xz = 0.0E0;
  Double I_KINETIC_Gxy2z_F2xz = 0.0E0;
  Double I_KINETIC_Gx3z_F2xz = 0.0E0;
  Double I_KINETIC_G4y_F2xz = 0.0E0;
  Double I_KINETIC_G3yz_F2xz = 0.0E0;
  Double I_KINETIC_G2y2z_F2xz = 0.0E0;
  Double I_KINETIC_Gy3z_F2xz = 0.0E0;
  Double I_KINETIC_G4z_F2xz = 0.0E0;
  Double I_KINETIC_G4x_Fx2y = 0.0E0;
  Double I_KINETIC_G3xy_Fx2y = 0.0E0;
  Double I_KINETIC_G3xz_Fx2y = 0.0E0;
  Double I_KINETIC_G2x2y_Fx2y = 0.0E0;
  Double I_KINETIC_G2xyz_Fx2y = 0.0E0;
  Double I_KINETIC_G2x2z_Fx2y = 0.0E0;
  Double I_KINETIC_Gx3y_Fx2y = 0.0E0;
  Double I_KINETIC_Gx2yz_Fx2y = 0.0E0;
  Double I_KINETIC_Gxy2z_Fx2y = 0.0E0;
  Double I_KINETIC_Gx3z_Fx2y = 0.0E0;
  Double I_KINETIC_G4y_Fx2y = 0.0E0;
  Double I_KINETIC_G3yz_Fx2y = 0.0E0;
  Double I_KINETIC_G2y2z_Fx2y = 0.0E0;
  Double I_KINETIC_Gy3z_Fx2y = 0.0E0;
  Double I_KINETIC_G4z_Fx2y = 0.0E0;
  Double I_KINETIC_G4x_Fxyz = 0.0E0;
  Double I_KINETIC_G3xy_Fxyz = 0.0E0;
  Double I_KINETIC_G3xz_Fxyz = 0.0E0;
  Double I_KINETIC_G2x2y_Fxyz = 0.0E0;
  Double I_KINETIC_G2xyz_Fxyz = 0.0E0;
  Double I_KINETIC_G2x2z_Fxyz = 0.0E0;
  Double I_KINETIC_Gx3y_Fxyz = 0.0E0;
  Double I_KINETIC_Gx2yz_Fxyz = 0.0E0;
  Double I_KINETIC_Gxy2z_Fxyz = 0.0E0;
  Double I_KINETIC_Gx3z_Fxyz = 0.0E0;
  Double I_KINETIC_G4y_Fxyz = 0.0E0;
  Double I_KINETIC_G3yz_Fxyz = 0.0E0;
  Double I_KINETIC_G2y2z_Fxyz = 0.0E0;
  Double I_KINETIC_Gy3z_Fxyz = 0.0E0;
  Double I_KINETIC_G4z_Fxyz = 0.0E0;
  Double I_KINETIC_G4x_Fx2z = 0.0E0;
  Double I_KINETIC_G3xy_Fx2z = 0.0E0;
  Double I_KINETIC_G3xz_Fx2z = 0.0E0;
  Double I_KINETIC_G2x2y_Fx2z = 0.0E0;
  Double I_KINETIC_G2xyz_Fx2z = 0.0E0;
  Double I_KINETIC_G2x2z_Fx2z = 0.0E0;
  Double I_KINETIC_Gx3y_Fx2z = 0.0E0;
  Double I_KINETIC_Gx2yz_Fx2z = 0.0E0;
  Double I_KINETIC_Gxy2z_Fx2z = 0.0E0;
  Double I_KINETIC_Gx3z_Fx2z = 0.0E0;
  Double I_KINETIC_G4y_Fx2z = 0.0E0;
  Double I_KINETIC_G3yz_Fx2z = 0.0E0;
  Double I_KINETIC_G2y2z_Fx2z = 0.0E0;
  Double I_KINETIC_Gy3z_Fx2z = 0.0E0;
  Double I_KINETIC_G4z_Fx2z = 0.0E0;
  Double I_KINETIC_G4x_F3y = 0.0E0;
  Double I_KINETIC_G3xy_F3y = 0.0E0;
  Double I_KINETIC_G3xz_F3y = 0.0E0;
  Double I_KINETIC_G2x2y_F3y = 0.0E0;
  Double I_KINETIC_G2xyz_F3y = 0.0E0;
  Double I_KINETIC_G2x2z_F3y = 0.0E0;
  Double I_KINETIC_Gx3y_F3y = 0.0E0;
  Double I_KINETIC_Gx2yz_F3y = 0.0E0;
  Double I_KINETIC_Gxy2z_F3y = 0.0E0;
  Double I_KINETIC_Gx3z_F3y = 0.0E0;
  Double I_KINETIC_G4y_F3y = 0.0E0;
  Double I_KINETIC_G3yz_F3y = 0.0E0;
  Double I_KINETIC_G2y2z_F3y = 0.0E0;
  Double I_KINETIC_Gy3z_F3y = 0.0E0;
  Double I_KINETIC_G4z_F3y = 0.0E0;
  Double I_KINETIC_G4x_F2yz = 0.0E0;
  Double I_KINETIC_G3xy_F2yz = 0.0E0;
  Double I_KINETIC_G3xz_F2yz = 0.0E0;
  Double I_KINETIC_G2x2y_F2yz = 0.0E0;
  Double I_KINETIC_G2xyz_F2yz = 0.0E0;
  Double I_KINETIC_G2x2z_F2yz = 0.0E0;
  Double I_KINETIC_Gx3y_F2yz = 0.0E0;
  Double I_KINETIC_Gx2yz_F2yz = 0.0E0;
  Double I_KINETIC_Gxy2z_F2yz = 0.0E0;
  Double I_KINETIC_Gx3z_F2yz = 0.0E0;
  Double I_KINETIC_G4y_F2yz = 0.0E0;
  Double I_KINETIC_G3yz_F2yz = 0.0E0;
  Double I_KINETIC_G2y2z_F2yz = 0.0E0;
  Double I_KINETIC_Gy3z_F2yz = 0.0E0;
  Double I_KINETIC_G4z_F2yz = 0.0E0;
  Double I_KINETIC_G4x_Fy2z = 0.0E0;
  Double I_KINETIC_G3xy_Fy2z = 0.0E0;
  Double I_KINETIC_G3xz_Fy2z = 0.0E0;
  Double I_KINETIC_G2x2y_Fy2z = 0.0E0;
  Double I_KINETIC_G2xyz_Fy2z = 0.0E0;
  Double I_KINETIC_G2x2z_Fy2z = 0.0E0;
  Double I_KINETIC_Gx3y_Fy2z = 0.0E0;
  Double I_KINETIC_Gx2yz_Fy2z = 0.0E0;
  Double I_KINETIC_Gxy2z_Fy2z = 0.0E0;
  Double I_KINETIC_Gx3z_Fy2z = 0.0E0;
  Double I_KINETIC_G4y_Fy2z = 0.0E0;
  Double I_KINETIC_G3yz_Fy2z = 0.0E0;
  Double I_KINETIC_G2y2z_Fy2z = 0.0E0;
  Double I_KINETIC_Gy3z_Fy2z = 0.0E0;
  Double I_KINETIC_G4z_Fy2z = 0.0E0;
  Double I_KINETIC_G4x_F3z = 0.0E0;
  Double I_KINETIC_G3xy_F3z = 0.0E0;
  Double I_KINETIC_G3xz_F3z = 0.0E0;
  Double I_KINETIC_G2x2y_F3z = 0.0E0;
  Double I_KINETIC_G2xyz_F3z = 0.0E0;
  Double I_KINETIC_G2x2z_F3z = 0.0E0;
  Double I_KINETIC_Gx3y_F3z = 0.0E0;
  Double I_KINETIC_Gx2yz_F3z = 0.0E0;
  Double I_KINETIC_Gxy2z_F3z = 0.0E0;
  Double I_KINETIC_Gx3z_F3z = 0.0E0;
  Double I_KINETIC_G4y_F3z = 0.0E0;
  Double I_KINETIC_G3yz_F3z = 0.0E0;
  Double I_KINETIC_G2y2z_F3z = 0.0E0;
  Double I_KINETIC_Gy3z_F3z = 0.0E0;
  Double I_KINETIC_G4z_F3z = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_Px_vrr = PBX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Px_vrr = PBX*I_TWOBODYOVERLAP_Dxy_S_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Px_vrr = PBX*I_TWOBODYOVERLAP_Dxz_S_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Px_vrr = PBX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Px_vrr = PBX*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Px_vrr = PBX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_D2x_Py_vrr = PBY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Py_vrr = PBY*I_TWOBODYOVERLAP_Dxy_S_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Py_vrr = PBY*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Py_vrr = PBY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Py_vrr = PBY*I_TWOBODYOVERLAP_Dyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Py_vrr = PBY*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_D2x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Dxz_S_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Dyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_S_vrr = PAX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_S_vrr = PAY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_S_vrr = PAX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_S_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_S_vrr = PAX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_S_vrr = PAY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_S_vrr = PAY*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_S_vrr = PAZ*I_TWOBODYOVERLAP_D2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_Px_vrr = PBX*I_TWOBODYOVERLAP_F3x_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Px_vrr = PBX*I_TWOBODYOVERLAP_F2xy_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Px_vrr = PBX*I_TWOBODYOVERLAP_F2xz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Px_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Px_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_S_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Px_vrr = PBX*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Px_vrr = PBX*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Py_vrr = PBY*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Py_vrr = PBY*I_TWOBODYOVERLAP_F2xy_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Py_vrr = PBY*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Py_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Py_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Py_vrr = PBY*I_TWOBODYOVERLAP_F3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_F2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Py_vrr = PBY*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_S_vrr = PAX*I_TWOBODYOVERLAP_F3x_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_S_vrr = PAY*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_S_vrr = PAY*I_TWOBODYOVERLAP_F2xy_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_S_vrr = PAX*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_S_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_S_vrr = PAX*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_S_vrr = PAY*I_TWOBODYOVERLAP_F3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_S_vrr = PAY*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_S_vrr = PAZ*I_TWOBODYOVERLAP_F3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 14 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_D2x_vrr = PBX*I_TWOBODYOVERLAP_F3x_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2x_vrr = PBX*I_TWOBODYOVERLAP_F2xy_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2x_vrr = PBX*I_TWOBODYOVERLAP_F2xz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2x_vrr = PBX*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2y_vrr = PBY*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2y_vrr = PBY*I_TWOBODYOVERLAP_F2xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2y_vrr = PBY*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2y_vrr = PBY*I_TWOBODYOVERLAP_F3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_F2yz_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_Px_vrr = PAX*I_TWOBODYOVERLAP_F3x_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Px_vrr = PAY*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Px_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Px_vrr = PAX*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Px_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Px_vrr = PAX*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Px_vrr = PAY*I_TWOBODYOVERLAP_F3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Px_vrr = PAY*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_Py_vrr = PAX*I_TWOBODYOVERLAP_F3x_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Py_vrr = PAY*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Py_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Py_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Py_vrr = PAX*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Py_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Py_vrr = PAX*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4y_Py_vrr = PAY*I_TWOBODYOVERLAP_F3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Py_vrr = PAY*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4x_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3x_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Pz_vrr = PAY*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Pz_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4y_Pz_vrr = PAY*I_TWOBODYOVERLAP_F3y_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Pz_vrr = PAY*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 10 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_D2x_vrr = PBX*I_TWOBODYOVERLAP_G4x_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2x_vrr = PBX*I_TWOBODYOVERLAP_G3xy_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2x_vrr = PBX*I_TWOBODYOVERLAP_G3xz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2x_vrr = PBX*I_TWOBODYOVERLAP_G4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_G3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2x_vrr = PBX*I_TWOBODYOVERLAP_G4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G3xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G4y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G3yz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_D2y_vrr = PBY*I_TWOBODYOVERLAP_G4x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2y_vrr = PBY*I_TWOBODYOVERLAP_G3xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2y_vrr = PBY*I_TWOBODYOVERLAP_G3xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2y_vrr = PBY*I_TWOBODYOVERLAP_G4y_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_G3yz_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2y_vrr = PBY*I_TWOBODYOVERLAP_G4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G4x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 15 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_Px_vrr = PAX*I_TWOBODYOVERLAP_G4x_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Px_vrr = PAY*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Px_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Px_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Px_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Px_vrr = PAX*I_TWOBODYOVERLAP_G4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Px_vrr = PAX*I_TWOBODYOVERLAP_G4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_Px_vrr = PAY*I_TWOBODYOVERLAP_G4y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Px_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Px_vrr = PAY*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_Px_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_Py_vrr = PAX*I_TWOBODYOVERLAP_G4x_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Py_vrr = PAY*I_TWOBODYOVERLAP_G4x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Py_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Py_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Py_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Py_vrr = PAX*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Py_vrr = PAX*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5y_Py_vrr = PAY*I_TWOBODYOVERLAP_G4y_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Py_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Py_vrr = PAY*I_TWOBODYOVERLAP_G4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_Py_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5x_Pz_vrr = PAX*I_TWOBODYOVERLAP_G4x_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Pz_vrr = PAY*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Pz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Pz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Pz_vrr = PAX*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Pz_vrr = PAX*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5y_Pz_vrr = PAY*I_TWOBODYOVERLAP_G4y_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Pz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Pz_vrr = PAY*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_F3x_vrr = PBX*I_TWOBODYOVERLAP_G4x_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3x_vrr = PBX*I_TWOBODYOVERLAP_G3xy_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3x_vrr = PBX*I_TWOBODYOVERLAP_G3xz_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3x_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F3x_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3x_vrr = PBX*I_TWOBODYOVERLAP_G4y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_G3yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3x_vrr = PBX*I_TWOBODYOVERLAP_G4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G3xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G3xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G4y_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G3yz_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G4x_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G3xy_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G3xz_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G3yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G4x_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G3xy_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G3xz_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G3yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4x_F3y_vrr = PBY*I_TWOBODYOVERLAP_G4x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3y_vrr = PBY*I_TWOBODYOVERLAP_G3xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3y_vrr = PBY*I_TWOBODYOVERLAP_G3xz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3y_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F3y_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3y_vrr = PBY*I_TWOBODYOVERLAP_G4y_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_G3yz_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3y_vrr = PBY*I_TWOBODYOVERLAP_G4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G3xy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G3xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G4y_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G3yz_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4x_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G4x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G4z_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_D2x_vrr = PAX*I_TWOBODYOVERLAP_G4x_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2x_vrr = PAY*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2x_vrr = PAY*I_TWOBODYOVERLAP_G3xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2x_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2x_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2x_vrr = PAX*I_TWOBODYOVERLAP_G4y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2x_vrr = PAX*I_TWOBODYOVERLAP_G4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2x_vrr = PAY*I_TWOBODYOVERLAP_G4y_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2x_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2x_vrr = PAY*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G4z_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_G4x_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_G4x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_G4y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_G4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_G4y_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_G4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_G4x_Dxz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_G4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_G4z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_G4y_Dxz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Dxz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_D2y_vrr = PAX*I_TWOBODYOVERLAP_G4x_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2y_vrr = PAY*I_TWOBODYOVERLAP_G4x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2y_vrr = PAY*I_TWOBODYOVERLAP_G3xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2y_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2y_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2y_vrr = PAX*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2y_vrr = PAX*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2y_vrr = PAY*I_TWOBODYOVERLAP_G4y_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2y_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2y_vrr = PAY*I_TWOBODYOVERLAP_G4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G4z_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_G4x_Dyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_G4x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_G4y_Dyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_G4z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Dyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5x_D2z_vrr = PAX*I_TWOBODYOVERLAP_G4x_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2z_vrr = PAY*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2z_vrr = PAX*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2z_vrr = PAX*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2z_vrr = PAY*I_TWOBODYOVERLAP_G4y_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2z_vrr = PAY*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 50 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_F3x_vrr = PBX*I_TWOBODYOVERLAP_H5x_D2x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F3x_vrr = PBX*I_TWOBODYOVERLAP_H4xy_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F3x_vrr = PBX*I_TWOBODYOVERLAP_H4xz_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F3x_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F3x_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F3x_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5y_F3x_vrr = PBX*I_TWOBODYOVERLAP_H5y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_H4yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H5z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H4xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H5y_D2x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H4yz_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_D2x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H5x_D2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H4xy_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H4xz_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_Dxy_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H5x_D2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H4xy_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H4xz_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H4yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5x_F3y_vrr = PBY*I_TWOBODYOVERLAP_H5x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F3y_vrr = PBY*I_TWOBODYOVERLAP_H4xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F3y_vrr = PBY*I_TWOBODYOVERLAP_H4xz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F3y_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F3y_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F3y_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5y_F3y_vrr = PBY*I_TWOBODYOVERLAP_H5y_D2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_H4yz_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2y3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H5z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5x_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_D2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H4xy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H4xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H5y_D2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H4yz_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5x_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H5x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H5y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2y3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H5z_D2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_F3x_vrr = PAX*I_TWOBODYOVERLAP_H5x_F3x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_I5xy_F3x_vrr = PAY*I_TWOBODYOVERLAP_H5x_F3x_vrr;
    Double I_TWOBODYOVERLAP_I5xz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H5x_F3x_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_F3x_vrr = PAY*I_TWOBODYOVERLAP_H4xy_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_F3x_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_F3x_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F3x_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_F3x_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_F3x_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_F3x_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_F3x_vrr = PAX*I_TWOBODYOVERLAP_H5y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_F3x_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_F3x_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H2y3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_F3x_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_F3x_vrr = PAX*I_TWOBODYOVERLAP_H5z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I6y_F3x_vrr = PAY*I_TWOBODYOVERLAP_H5y_F3x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_I5yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H5y_F3x_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F3x_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_F3x_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_F3x_vrr = PAY*I_TWOBODYOVERLAP_H5z_F3x_vrr;
    Double I_TWOBODYOVERLAP_I6z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_H5z_F3x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_I6x_F2xy_vrr = PAX*I_TWOBODYOVERLAP_H5x_F2xy_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I5xy_F2xy_vrr = PAY*I_TWOBODYOVERLAP_H5x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_I5xz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H5x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_H4xy_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_F2xy_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_F2xy_vrr = PAX*I_TWOBODYOVERLAP_H5y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_H2y3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_H5z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I6y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_H5y_F2xy_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I5yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H5y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_H5z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I6z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_H5z_F2xy_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_I6x_F2xz_vrr = PAX*I_TWOBODYOVERLAP_H5x_F2xz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I5xy_F2xz_vrr = PAY*I_TWOBODYOVERLAP_H5x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_I5xz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H5x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_H4xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F2xz_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_F2xz_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_F2xz_vrr = PAX*I_TWOBODYOVERLAP_H5y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_H2y3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_H5z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I6y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_H5y_F2xz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_I5yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H5y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_H5z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_I6z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_H5z_F2xz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I6x_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_H5x_Fx2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_H5x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_H5y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_H5z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_I6y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_H5y_Fx2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_H5z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I6z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Fx2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_I6x_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_H5x_Fxyz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_H5x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_H5y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_H5z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I6y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_H5y_Fxyz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_H5z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I6z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Fxyz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I6x_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_H5x_Fx2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_H5x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_H5y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_H5z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I6y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_H5y_Fx2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_H5z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_I6z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Fx2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I6x_F3y_vrr = PAX*I_TWOBODYOVERLAP_H5x_F3y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_TWOBODYOVERLAP_I5xy_F3y_vrr = PAY*I_TWOBODYOVERLAP_H5x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_I5xz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H5x_F3y_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_F3y_vrr = PAY*I_TWOBODYOVERLAP_H4xy_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_F3y_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_F3y_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F3y_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_F3y_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_F3y_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_F3y_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_F3y_vrr = PAX*I_TWOBODYOVERLAP_H5y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_F3y_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_F3y_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_F3y_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_F3y_vrr = PAX*I_TWOBODYOVERLAP_H5z_F3y_vrr;
    Double I_TWOBODYOVERLAP_I6y_F3y_vrr = PAY*I_TWOBODYOVERLAP_H5y_F3y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I5yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H5y_F3y_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F3y_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_F3y_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_F3y_vrr = PAY*I_TWOBODYOVERLAP_H5z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_I6z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_H5z_F3y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_TWOBODYOVERLAP_I6x_F2yz_vrr = PAX*I_TWOBODYOVERLAP_H5x_F2yz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_I5xy_F2yz_vrr = PAY*I_TWOBODYOVERLAP_H5x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I5xz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H5x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_H4xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_F2yz_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_F2yz_vrr = PAX*I_TWOBODYOVERLAP_H5y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_H5z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_I6y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_H5y_F2yz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I5yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H5y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_H5z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I6z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_H5z_F2yz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_I6x_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_H5x_Fy2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_H5x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_H5y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_H5z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_I6y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_H5y_Fy2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_H5z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I6z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Fy2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I6x_F3z_vrr = PAX*I_TWOBODYOVERLAP_H5x_F3z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_TWOBODYOVERLAP_I5xy_F3z_vrr = PAY*I_TWOBODYOVERLAP_H5x_F3z_vrr;
    Double I_TWOBODYOVERLAP_I5xz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H5x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_F3z_vrr = PAY*I_TWOBODYOVERLAP_H4xy_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H4xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H4xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_F3z_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F3z_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_F3z_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_F3z_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_F3z_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_F3z_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_F3z_vrr = PAX*I_TWOBODYOVERLAP_H5y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_F3z_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_F3z_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_F3z_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_F3z_vrr = PAX*I_TWOBODYOVERLAP_H5z_F3z_vrr;
    Double I_TWOBODYOVERLAP_I6y_F3z_vrr = PAY*I_TWOBODYOVERLAP_H5y_F3z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_TWOBODYOVERLAP_I5yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H5y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H4yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_F3z_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_F3z_vrr = PAY*I_TWOBODYOVERLAP_H5z_F3z_vrr;
    Double I_TWOBODYOVERLAP_I6z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_H5z_F3z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H5z_D2z_vrr;

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
     * shell quartet name: SQ_KINETIC_D_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_KINETIC_D2x_Px_vrr = PBX*I_KINETIC_D2x_S_vrr+2*oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_Dxy_Px_vrr = PBX*I_KINETIC_Dxy_S_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_KINETIC_Dxz_Px_vrr = PBX*I_KINETIC_Dxz_S_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_KINETIC_D2y_Px_vrr = PBX*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_Dyz_Px_vrr = PBX*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_KINETIC_D2z_Px_vrr = PBX*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_KINETIC_D2x_Py_vrr = PBY*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_Dxy_Py_vrr = PBY*I_KINETIC_Dxy_S_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_KINETIC_Dxz_Py_vrr = PBY*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Py_vrr;
    Double I_KINETIC_D2y_Py_vrr = PBY*I_KINETIC_D2y_S_vrr+2*oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_Dyz_Py_vrr = PBY*I_KINETIC_Dyz_S_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Py_vrr;
    Double I_KINETIC_D2z_Py_vrr = PBY*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_KINETIC_D2x_Pz_vrr = PBZ*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_Dxy_Pz_vrr = PBZ*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_KINETIC_Dxz_Pz_vrr = PBZ*I_KINETIC_Dxz_S_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Pz_vrr;
    Double I_KINETIC_D2y_Pz_vrr = PBZ*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_Dyz_Pz_vrr = PBZ*I_KINETIC_Dyz_S_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Pz_vrr;
    Double I_KINETIC_D2z_Pz_vrr = PBZ*I_KINETIC_D2z_S_vrr+2*oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_F3x_S_vrr = PAX*I_KINETIC_D2x_S_vrr+2*oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_S_vrr-2*bdz*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_F2xy_S_vrr = PAY*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_S_vrr = PAZ*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_S_vrr = PAX*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_S_vrr = PAZ*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_S_vrr = PAX*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_S_vrr = PAY*I_KINETIC_D2y_S_vrr+2*oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_F2yz_S_vrr = PAZ*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_S_vrr = PAY*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_S_vrr = PAZ*I_KINETIC_D2z_S_vrr+2*oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     ************************************************************/
    Double I_KINETIC_F3x_Px_vrr = PBX*I_KINETIC_F3x_S_vrr+3*oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_KINETIC_F2xy_Px_vrr = PBX*I_KINETIC_F2xy_S_vrr+2*oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_KINETIC_F2xz_Px_vrr = PBX*I_KINETIC_F2xz_S_vrr+2*oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_KINETIC_Fx2y_Px_vrr = PBX*I_KINETIC_Fx2y_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_KINETIC_Fxyz_Px_vrr = PBX*I_KINETIC_Fxyz_S_vrr+oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Px_vrr;
    Double I_KINETIC_Fx2z_Px_vrr = PBX*I_KINETIC_Fx2z_S_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_KINETIC_F3y_Px_vrr = PBX*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_KINETIC_F2yz_Px_vrr = PBX*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_KINETIC_Fy2z_Px_vrr = PBX*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_KINETIC_F3z_Px_vrr = PBX*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_KINETIC_F3x_Py_vrr = PBY*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_KINETIC_F2xy_Py_vrr = PBY*I_KINETIC_F2xy_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_KINETIC_F2xz_Py_vrr = PBY*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_KINETIC_Fx2y_Py_vrr = PBY*I_KINETIC_Fx2y_S_vrr+2*oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_KINETIC_Fxyz_Py_vrr = PBY*I_KINETIC_Fxyz_S_vrr+oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Py_vrr;
    Double I_KINETIC_Fx2z_Py_vrr = PBY*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Py_vrr;
    Double I_KINETIC_F3y_Py_vrr = PBY*I_KINETIC_F3y_S_vrr+3*oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_KINETIC_F2yz_Py_vrr = PBY*I_KINETIC_F2yz_S_vrr+2*oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Py_vrr;
    Double I_KINETIC_Fy2z_Py_vrr = PBY*I_KINETIC_Fy2z_S_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Py_vrr;
    Double I_KINETIC_F3z_Py_vrr = PBY*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_KINETIC_F3x_Pz_vrr = PBZ*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_KINETIC_F2xy_Pz_vrr = PBZ*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Pz_vrr;
    Double I_KINETIC_F2xz_Pz_vrr = PBZ*I_KINETIC_F2xz_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Pz_vrr;
    Double I_KINETIC_Fx2y_Pz_vrr = PBZ*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Pz_vrr;
    Double I_KINETIC_Fxyz_Pz_vrr = PBZ*I_KINETIC_Fxyz_S_vrr+oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Pz_vrr;
    Double I_KINETIC_Fx2z_Pz_vrr = PBZ*I_KINETIC_Fx2z_S_vrr+2*oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Pz_vrr;
    Double I_KINETIC_F3y_Pz_vrr = PBZ*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_KINETIC_F2yz_Pz_vrr = PBZ*I_KINETIC_F2yz_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Pz_vrr;
    Double I_KINETIC_Fy2z_Pz_vrr = PBZ*I_KINETIC_Fy2z_S_vrr+2*oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Pz_vrr;
    Double I_KINETIC_F3z_Pz_vrr = PBZ*I_KINETIC_F3z_S_vrr+3*oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_KINETIC_G4x_S_vrr = PAX*I_KINETIC_F3x_S_vrr+3*oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G3xy_S_vrr = PAY*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_S_vrr = PAZ*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_S_vrr = PAY*I_KINETIC_F2xy_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G2xyz_S_vrr = PAZ*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_S_vrr = PAZ*I_KINETIC_F2xz_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_Gx3y_S_vrr = PAX*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_S_vrr = PAZ*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_S_vrr = PAY*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_KINETIC_Gx3z_S_vrr = PAX*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_S_vrr = PAY*I_KINETIC_F3y_S_vrr+3*oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_G3yz_S_vrr = PAZ*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_S_vrr = PAZ*I_KINETIC_F2yz_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_Gy3z_S_vrr = PAY*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_S_vrr = PAZ*I_KINETIC_F3z_S_vrr+3*oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 14 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_KINETIC_F3x_D2x_vrr = PBX*I_KINETIC_F3x_Px_vrr+3*oned2z*I_KINETIC_D2x_Px_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2x_vrr-adz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_F2xy_D2x_vrr = PBX*I_KINETIC_F2xy_Px_vrr+2*oned2z*I_KINETIC_Dxy_Px_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2x_vrr-adz*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_D2x_vrr = PBX*I_KINETIC_F2xz_Px_vrr+2*oned2z*I_KINETIC_Dxz_Px_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2x_vrr-adz*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_D2x_vrr = PBX*I_KINETIC_Fx2y_Px_vrr+oned2z*I_KINETIC_D2y_Px_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2x_vrr-adz*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_D2x_vrr = PBX*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_Dyz_Px_vrr+oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2x_vrr-adz*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_D2x_vrr = PBX*I_KINETIC_Fx2z_Px_vrr+oned2z*I_KINETIC_D2z_Px_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_D2x_vrr = PBX*I_KINETIC_F3y_Px_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2x_vrr-adz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_F2yz_D2x_vrr = PBX*I_KINETIC_F2yz_Px_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2x_vrr-adz*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_D2x_vrr = PBX*I_KINETIC_Fy2z_Px_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_D2x_vrr = PBX*I_KINETIC_F3z_Px_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2x_vrr-adz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_F3x_Dxy_vrr = PBY*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_KINETIC_F2xy_Dxy_vrr = PBY*I_KINETIC_F2xy_Px_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_KINETIC_F2xz_Dxy_vrr = PBY*I_KINETIC_F2xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_KINETIC_Fx2y_Dxy_vrr = PBY*I_KINETIC_Fx2y_Px_vrr+2*oned2z*I_KINETIC_Dxy_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_KINETIC_Fxyz_Dxy_vrr = PBY*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_Dxz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dxy_vrr;
    Double I_KINETIC_Fx2z_Dxy_vrr = PBY*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr;
    Double I_KINETIC_F3y_Dxy_vrr = PBY*I_KINETIC_F3y_Px_vrr+3*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_KINETIC_F2yz_Dxy_vrr = PBY*I_KINETIC_F2yz_Px_vrr+2*oned2z*I_KINETIC_Dyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_KINETIC_Fy2z_Dxy_vrr = PBY*I_KINETIC_Fy2z_Px_vrr+oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dxy_vrr;
    Double I_KINETIC_F3z_Dxy_vrr = PBY*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_KINETIC_F3x_Dxz_vrr = PBZ*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_KINETIC_F3y_Dxz_vrr = PBZ*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_KINETIC_F3z_Dxz_vrr = PBZ*I_KINETIC_F3z_Px_vrr+3*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_KINETIC_F3x_D2y_vrr = PBY*I_KINETIC_F3x_Py_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2y_vrr-adz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_F2xy_D2y_vrr = PBY*I_KINETIC_F2xy_Py_vrr+oned2z*I_KINETIC_D2x_Py_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2y_vrr-adz*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_D2y_vrr = PBY*I_KINETIC_F2xz_Py_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2y_vrr-adz*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_D2y_vrr = PBY*I_KINETIC_Fx2y_Py_vrr+2*oned2z*I_KINETIC_Dxy_Py_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2y_vrr-adz*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_D2y_vrr = PBY*I_KINETIC_Fxyz_Py_vrr+oned2z*I_KINETIC_Dxz_Py_vrr+oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2y_vrr-adz*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_D2y_vrr = PBY*I_KINETIC_Fx2z_Py_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_D2y_vrr = PBY*I_KINETIC_F3y_Py_vrr+3*oned2z*I_KINETIC_D2y_Py_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2y_vrr-adz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_F2yz_D2y_vrr = PBY*I_KINETIC_F2yz_Py_vrr+2*oned2z*I_KINETIC_Dyz_Py_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2y_vrr-adz*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_D2y_vrr = PBY*I_KINETIC_Fy2z_Py_vrr+oned2z*I_KINETIC_D2z_Py_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_D2y_vrr = PBY*I_KINETIC_F3z_Py_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2y_vrr-adz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_F3x_Dyz_vrr = PBZ*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_KINETIC_F3y_Dyz_vrr = PBZ*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_KINETIC_F3z_Dyz_vrr = PBZ*I_KINETIC_F3z_Py_vrr+3*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_KINETIC_F3x_D2z_vrr = PBZ*I_KINETIC_F3x_Pz_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2z_vrr-adz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_F2xy_D2z_vrr = PBZ*I_KINETIC_F2xy_Pz_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2z_vrr-adz*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_D2z_vrr = PBZ*I_KINETIC_F2xz_Pz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2z_vrr-adz*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_D2z_vrr = PBZ*I_KINETIC_Fx2y_Pz_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2z_vrr-adz*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_D2z_vrr = PBZ*I_KINETIC_Fxyz_Pz_vrr+oned2z*I_KINETIC_Dxy_Pz_vrr+oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2z_vrr-adz*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_D2z_vrr = PBZ*I_KINETIC_Fx2z_Pz_vrr+2*oned2z*I_KINETIC_Dxz_Pz_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_D2z_vrr = PBZ*I_KINETIC_F3y_Pz_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2z_vrr-adz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_F2yz_D2z_vrr = PBZ*I_KINETIC_F2yz_Pz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2z_vrr-adz*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_D2z_vrr = PBZ*I_KINETIC_Fy2z_Pz_vrr+2*oned2z*I_KINETIC_Dyz_Pz_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_D2z_vrr = PBZ*I_KINETIC_F3z_Pz_vrr+3*oned2z*I_KINETIC_D2z_Pz_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2z_vrr-adz*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_KINETIC_G4x_Px_vrr = PAX*I_KINETIC_F3x_Px_vrr+3*oned2z*I_KINETIC_D2x_Px_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Px_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_G3xy_Px_vrr = PAY*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_KINETIC_G3xz_Px_vrr = PAZ*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_KINETIC_G2x2y_Px_vrr = PAY*I_KINETIC_F2xy_Px_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Px_vrr-bdz*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_G2xyz_Px_vrr = PAZ*I_KINETIC_F2xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_KINETIC_G2x2z_Px_vrr = PAZ*I_KINETIC_F2xz_Px_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Px_vrr-bdz*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_Gx3y_Px_vrr = PAX*I_KINETIC_F3y_Px_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_KINETIC_Gx2yz_Px_vrr = PAZ*I_KINETIC_Fx2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_KINETIC_Gxy2z_Px_vrr = PAY*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_KINETIC_Gx3z_Px_vrr = PAX*I_KINETIC_F3z_Px_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_KINETIC_G4y_Px_vrr = PAY*I_KINETIC_F3y_Px_vrr+3*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Px_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_G3yz_Px_vrr = PAZ*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_KINETIC_G2y2z_Px_vrr = PAZ*I_KINETIC_F2yz_Px_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Px_vrr-bdz*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_Gy3z_Px_vrr = PAY*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_KINETIC_G4z_Px_vrr = PAZ*I_KINETIC_F3z_Px_vrr+3*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Px_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_KINETIC_G4x_Py_vrr = PAX*I_KINETIC_F3x_Py_vrr+3*oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Py_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_G3xy_Py_vrr = PAY*I_KINETIC_F3x_Py_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_KINETIC_G3xz_Py_vrr = PAZ*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_KINETIC_G2x2y_Py_vrr = PAY*I_KINETIC_F2xy_Py_vrr+oned2z*I_KINETIC_D2x_Py_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Py_vrr-bdz*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_G2xyz_Py_vrr = PAZ*I_KINETIC_F2xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_KINETIC_G2x2z_Py_vrr = PAZ*I_KINETIC_F2xz_Py_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Py_vrr-bdz*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_Gx3y_Py_vrr = PAX*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_KINETIC_Gx2yz_Py_vrr = PAZ*I_KINETIC_Fx2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_KINETIC_Gxy2z_Py_vrr = PAY*I_KINETIC_Fx2z_Py_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_KINETIC_Gx3z_Py_vrr = PAX*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_KINETIC_G4y_Py_vrr = PAY*I_KINETIC_F3y_Py_vrr+3*oned2z*I_KINETIC_D2y_Py_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Py_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_G3yz_Py_vrr = PAZ*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_KINETIC_G2y2z_Py_vrr = PAZ*I_KINETIC_F2yz_Py_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Py_vrr-bdz*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_Gy3z_Py_vrr = PAY*I_KINETIC_F3z_Py_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_KINETIC_G4z_Py_vrr = PAZ*I_KINETIC_F3z_Py_vrr+3*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Py_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_KINETIC_G4x_Pz_vrr = PAX*I_KINETIC_F3x_Pz_vrr+3*oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Pz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_G3xy_Pz_vrr = PAY*I_KINETIC_F3x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_KINETIC_G3xz_Pz_vrr = PAZ*I_KINETIC_F3x_Pz_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_KINETIC_G2x2y_Pz_vrr = PAY*I_KINETIC_F2xy_Pz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Pz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_G2xyz_Pz_vrr = PAZ*I_KINETIC_F2xy_Pz_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Pz_vrr;
    Double I_KINETIC_G2x2z_Pz_vrr = PAZ*I_KINETIC_F2xz_Pz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_Gx3y_Pz_vrr = PAX*I_KINETIC_F3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_KINETIC_Gx2yz_Pz_vrr = PAZ*I_KINETIC_Fx2y_Pz_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr;
    Double I_KINETIC_Gxy2z_Pz_vrr = PAY*I_KINETIC_Fx2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr;
    Double I_KINETIC_Gx3z_Pz_vrr = PAX*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_KINETIC_G4y_Pz_vrr = PAY*I_KINETIC_F3y_Pz_vrr+3*oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Pz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_G3yz_Pz_vrr = PAZ*I_KINETIC_F3y_Pz_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_KINETIC_G2y2z_Pz_vrr = PAZ*I_KINETIC_F2yz_Pz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_Gy3z_Pz_vrr = PAY*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_KINETIC_G4z_Pz_vrr = PAZ*I_KINETIC_F3z_Pz_vrr+3*oned2z*I_KINETIC_D2z_Pz_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Pz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 10 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_KINETIC_G4x_D2x_vrr = PBX*I_KINETIC_G4x_Px_vrr+4*oned2z*I_KINETIC_F3x_Px_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2x_vrr-adz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_G3xy_D2x_vrr = PBX*I_KINETIC_G3xy_Px_vrr+3*oned2z*I_KINETIC_F2xy_Px_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2x_vrr-adz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_D2x_vrr = PBX*I_KINETIC_G3xz_Px_vrr+3*oned2z*I_KINETIC_F2xz_Px_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2x_vrr-adz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_D2x_vrr = PBX*I_KINETIC_G2x2y_Px_vrr+2*oned2z*I_KINETIC_Fx2y_Px_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2x_vrr-adz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_G2xyz_D2x_vrr = PBX*I_KINETIC_G2xyz_Px_vrr+2*oned2z*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_D2x_vrr-adz*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_D2x_vrr = PBX*I_KINETIC_G2x2z_Px_vrr+2*oned2z*I_KINETIC_Fx2z_Px_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2x_vrr-adz*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_KINETIC_Gx3y_D2x_vrr = PBX*I_KINETIC_Gx3y_Px_vrr+oned2z*I_KINETIC_F3y_Px_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2x_vrr-adz*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_D2x_vrr = PBX*I_KINETIC_Gx2yz_Px_vrr+oned2z*I_KINETIC_F2yz_Px_vrr+oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr-adz*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_D2x_vrr = PBX*I_KINETIC_Gxy2z_Px_vrr+oned2z*I_KINETIC_Fy2z_Px_vrr+oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_KINETIC_Gx3z_D2x_vrr = PBX*I_KINETIC_Gx3z_Px_vrr+oned2z*I_KINETIC_F3z_Px_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2x_vrr-adz*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_D2x_vrr = PBX*I_KINETIC_G4y_Px_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2x_vrr-adz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_G3yz_D2x_vrr = PBX*I_KINETIC_G3yz_Px_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2x_vrr-adz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_D2x_vrr = PBX*I_KINETIC_G2y2z_Px_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2x_vrr-adz*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_KINETIC_Gy3z_D2x_vrr = PBX*I_KINETIC_Gy3z_Px_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2x_vrr-adz*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_D2x_vrr = PBX*I_KINETIC_G4z_Px_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2x_vrr-adz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_G4x_Dxy_vrr = PBY*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_KINETIC_G3xy_Dxy_vrr = PBY*I_KINETIC_G3xy_Px_vrr+oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Dxy_vrr;
    Double I_KINETIC_G3xz_Dxy_vrr = PBY*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Dxy_vrr;
    Double I_KINETIC_G2x2y_Dxy_vrr = PBY*I_KINETIC_G2x2y_Px_vrr+2*oned2z*I_KINETIC_F2xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_KINETIC_G2xyz_Dxy_vrr = PBY*I_KINETIC_G2xyz_Px_vrr+oned2z*I_KINETIC_F2xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Dxy_vrr;
    Double I_KINETIC_G2x2z_Dxy_vrr = PBY*I_KINETIC_G2x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dxy_vrr;
    Double I_KINETIC_Gx3y_Dxy_vrr = PBY*I_KINETIC_Gx3y_Px_vrr+3*oned2z*I_KINETIC_Fx2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr;
    Double I_KINETIC_Gx2yz_Dxy_vrr = PBY*I_KINETIC_Gx2yz_Px_vrr+2*oned2z*I_KINETIC_Fxyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Dxy_vrr;
    Double I_KINETIC_Gxy2z_Dxy_vrr = PBY*I_KINETIC_Gxy2z_Px_vrr+oned2z*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Dxy_vrr;
    Double I_KINETIC_Gx3z_Dxy_vrr = PBY*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr;
    Double I_KINETIC_G4y_Dxy_vrr = PBY*I_KINETIC_G4y_Px_vrr+4*oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_KINETIC_G3yz_Dxy_vrr = PBY*I_KINETIC_G3yz_Px_vrr+3*oned2z*I_KINETIC_F2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dxy_vrr;
    Double I_KINETIC_G2y2z_Dxy_vrr = PBY*I_KINETIC_G2y2z_Px_vrr+2*oned2z*I_KINETIC_Fy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dxy_vrr;
    Double I_KINETIC_Gy3z_Dxy_vrr = PBY*I_KINETIC_Gy3z_Px_vrr+oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr;
    Double I_KINETIC_G4z_Dxy_vrr = PBY*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_KINETIC_G4x_Dxz_vrr = PBZ*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_KINETIC_G3xy_Dxz_vrr = PBZ*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Dxz_vrr;
    Double I_KINETIC_G3xz_Dxz_vrr = PBZ*I_KINETIC_G3xz_Px_vrr+oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Dxz_vrr;
    Double I_KINETIC_G2x2y_Dxz_vrr = PBZ*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Dxz_vrr;
    Double I_KINETIC_Gx3y_Dxz_vrr = PBZ*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dxz_vrr;
    Double I_KINETIC_Gx3z_Dxz_vrr = PBZ*I_KINETIC_Gx3z_Px_vrr+3*oned2z*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dxz_vrr;
    Double I_KINETIC_G4y_Dxz_vrr = PBZ*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_KINETIC_G3yz_Dxz_vrr = PBZ*I_KINETIC_G3yz_Px_vrr+oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dxz_vrr;
    Double I_KINETIC_Gy3z_Dxz_vrr = PBZ*I_KINETIC_Gy3z_Px_vrr+3*oned2z*I_KINETIC_Fy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dxz_vrr;
    Double I_KINETIC_G4z_Dxz_vrr = PBZ*I_KINETIC_G4z_Px_vrr+4*oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_KINETIC_G4x_D2y_vrr = PBY*I_KINETIC_G4x_Py_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2y_vrr-adz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_G3xy_D2y_vrr = PBY*I_KINETIC_G3xy_Py_vrr+oned2z*I_KINETIC_F3x_Py_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2y_vrr-adz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_D2y_vrr = PBY*I_KINETIC_G3xz_Py_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2y_vrr-adz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_D2y_vrr = PBY*I_KINETIC_G2x2y_Py_vrr+2*oned2z*I_KINETIC_F2xy_Py_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2y_vrr-adz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_G2xyz_D2y_vrr = PBY*I_KINETIC_G2xyz_Py_vrr+oned2z*I_KINETIC_F2xz_Py_vrr+oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_D2y_vrr-adz*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_D2y_vrr = PBY*I_KINETIC_G2x2z_Py_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2y_vrr-adz*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_KINETIC_Gx3y_D2y_vrr = PBY*I_KINETIC_Gx3y_Py_vrr+3*oned2z*I_KINETIC_Fx2y_Py_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2y_vrr-adz*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_D2y_vrr = PBY*I_KINETIC_Gx2yz_Py_vrr+2*oned2z*I_KINETIC_Fxyz_Py_vrr+oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr-adz*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_D2y_vrr = PBY*I_KINETIC_Gxy2z_Py_vrr+oned2z*I_KINETIC_Fx2z_Py_vrr+oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_KINETIC_Gx3z_D2y_vrr = PBY*I_KINETIC_Gx3z_Py_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2y_vrr-adz*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_D2y_vrr = PBY*I_KINETIC_G4y_Py_vrr+4*oned2z*I_KINETIC_F3y_Py_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2y_vrr-adz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_G3yz_D2y_vrr = PBY*I_KINETIC_G3yz_Py_vrr+3*oned2z*I_KINETIC_F2yz_Py_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2y_vrr-adz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_D2y_vrr = PBY*I_KINETIC_G2y2z_Py_vrr+2*oned2z*I_KINETIC_Fy2z_Py_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2y_vrr-adz*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_KINETIC_Gy3z_D2y_vrr = PBY*I_KINETIC_Gy3z_Py_vrr+oned2z*I_KINETIC_F3z_Py_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2y_vrr-adz*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_D2y_vrr = PBY*I_KINETIC_G4z_Py_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2y_vrr-adz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_G4x_Dyz_vrr = PBZ*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Dyz_vrr;
    Double I_KINETIC_G3xy_Dyz_vrr = PBZ*I_KINETIC_G3xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Dyz_vrr;
    Double I_KINETIC_G3xz_Dyz_vrr = PBZ*I_KINETIC_G3xz_Py_vrr+oned2z*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Dyz_vrr;
    Double I_KINETIC_G2x2y_Dyz_vrr = PBZ*I_KINETIC_G2x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Dyz_vrr;
    Double I_KINETIC_Gx3y_Dyz_vrr = PBZ*I_KINETIC_Gx3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dyz_vrr;
    Double I_KINETIC_Gx3z_Dyz_vrr = PBZ*I_KINETIC_Gx3z_Py_vrr+3*oned2z*I_KINETIC_Fx2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dyz_vrr;
    Double I_KINETIC_G4y_Dyz_vrr = PBZ*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_KINETIC_G3yz_Dyz_vrr = PBZ*I_KINETIC_G3yz_Py_vrr+oned2z*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dyz_vrr;
    Double I_KINETIC_Gy3z_Dyz_vrr = PBZ*I_KINETIC_Gy3z_Py_vrr+3*oned2z*I_KINETIC_Fy2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dyz_vrr;
    Double I_KINETIC_G4z_Dyz_vrr = PBZ*I_KINETIC_G4z_Py_vrr+4*oned2z*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_KINETIC_G4x_D2z_vrr = PBZ*I_KINETIC_G4x_Pz_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2z_vrr-adz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_G3xy_D2z_vrr = PBZ*I_KINETIC_G3xy_Pz_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2z_vrr-adz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_D2z_vrr = PBZ*I_KINETIC_G3xz_Pz_vrr+oned2z*I_KINETIC_F3x_Pz_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2z_vrr-adz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_D2z_vrr = PBZ*I_KINETIC_G2x2y_Pz_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2z_vrr-adz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_G2xyz_D2z_vrr = PBZ*I_KINETIC_G2xyz_Pz_vrr+oned2z*I_KINETIC_F2xy_Pz_vrr+oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_D2z_vrr-adz*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_D2z_vrr = PBZ*I_KINETIC_G2x2z_Pz_vrr+2*oned2z*I_KINETIC_F2xz_Pz_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2z_vrr-adz*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_KINETIC_Gx3y_D2z_vrr = PBZ*I_KINETIC_Gx3y_Pz_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2z_vrr-adz*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_D2z_vrr = PBZ*I_KINETIC_Gx2yz_Pz_vrr+oned2z*I_KINETIC_Fx2y_Pz_vrr+oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr-adz*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_D2z_vrr = PBZ*I_KINETIC_Gxy2z_Pz_vrr+2*oned2z*I_KINETIC_Fxyz_Pz_vrr+oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_KINETIC_Gx3z_D2z_vrr = PBZ*I_KINETIC_Gx3z_Pz_vrr+3*oned2z*I_KINETIC_Fx2z_Pz_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2z_vrr-adz*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_D2z_vrr = PBZ*I_KINETIC_G4y_Pz_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2z_vrr-adz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_G3yz_D2z_vrr = PBZ*I_KINETIC_G3yz_Pz_vrr+oned2z*I_KINETIC_F3y_Pz_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2z_vrr-adz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_D2z_vrr = PBZ*I_KINETIC_G2y2z_Pz_vrr+2*oned2z*I_KINETIC_F2yz_Pz_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2z_vrr-adz*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_KINETIC_Gy3z_D2z_vrr = PBZ*I_KINETIC_Gy3z_Pz_vrr+3*oned2z*I_KINETIC_Fy2z_Pz_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2z_vrr-adz*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_D2z_vrr = PBZ*I_KINETIC_G4z_Pz_vrr+4*oned2z*I_KINETIC_F3z_Pz_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2z_vrr-adz*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 15 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     ************************************************************/
    Double I_KINETIC_H5x_Px_vrr = PAX*I_KINETIC_G4x_Px_vrr+4*oned2z*I_KINETIC_F3x_Px_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Px_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_KINETIC_H4xy_Px_vrr = PAY*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_KINETIC_H4xz_Px_vrr = PAZ*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_KINETIC_H3x2y_Px_vrr = PAY*I_KINETIC_G3xy_Px_vrr+oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Px_vrr-bdz*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_KINETIC_H3x2z_Px_vrr = PAZ*I_KINETIC_G3xz_Px_vrr+oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Px_vrr-bdz*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_KINETIC_H2x3y_Px_vrr = PAX*I_KINETIC_Gx3y_Px_vrr+oned2z*I_KINETIC_F3y_Px_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Px_vrr-bdz*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_KINETIC_H2x2yz_Px_vrr = PAZ*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_KINETIC_H2x3z_Px_vrr = PAX*I_KINETIC_Gx3z_Px_vrr+oned2z*I_KINETIC_F3z_Px_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Px_vrr-bdz*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_KINETIC_Hx4y_Px_vrr = PAX*I_KINETIC_G4y_Px_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_KINETIC_Hx4z_Px_vrr = PAX*I_KINETIC_G4z_Px_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_KINETIC_H5y_Px_vrr = PAY*I_KINETIC_G4y_Px_vrr+4*oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Px_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_KINETIC_H4yz_Px_vrr = PAZ*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_KINETIC_H3y2z_Px_vrr = PAZ*I_KINETIC_G3yz_Px_vrr+oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Px_vrr-bdz*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_KINETIC_H2y3z_Px_vrr = PAY*I_KINETIC_Gy3z_Px_vrr+oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Px_vrr-bdz*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_KINETIC_Hy4z_Px_vrr = PAY*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_KINETIC_H5z_Px_vrr = PAZ*I_KINETIC_G4z_Px_vrr+4*oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Px_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_KINETIC_H5x_Py_vrr = PAX*I_KINETIC_G4x_Py_vrr+4*oned2z*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Py_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_KINETIC_H4xy_Py_vrr = PAY*I_KINETIC_G4x_Py_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_KINETIC_H4xz_Py_vrr = PAZ*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_KINETIC_H3x2y_Py_vrr = PAY*I_KINETIC_G3xy_Py_vrr+oned2z*I_KINETIC_F3x_Py_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Py_vrr-bdz*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_KINETIC_H3x2z_Py_vrr = PAZ*I_KINETIC_G3xz_Py_vrr+oned2z*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Py_vrr-bdz*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_KINETIC_H2x3y_Py_vrr = PAX*I_KINETIC_Gx3y_Py_vrr+oned2z*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Py_vrr-bdz*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_KINETIC_H2x2yz_Py_vrr = PAZ*I_KINETIC_G2x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Py_vrr;
    Double I_KINETIC_H2x3z_Py_vrr = PAX*I_KINETIC_Gx3z_Py_vrr+oned2z*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Py_vrr-bdz*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_KINETIC_Hx4y_Py_vrr = PAX*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_KINETIC_Hx4z_Py_vrr = PAX*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_KINETIC_H5y_Py_vrr = PAY*I_KINETIC_G4y_Py_vrr+4*oned2z*I_KINETIC_F3y_Py_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Py_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_KINETIC_H4yz_Py_vrr = PAZ*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_KINETIC_H3y2z_Py_vrr = PAZ*I_KINETIC_G3yz_Py_vrr+oned2z*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Py_vrr-bdz*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_KINETIC_H2y3z_Py_vrr = PAY*I_KINETIC_Gy3z_Py_vrr+oned2z*I_KINETIC_F3z_Py_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Py_vrr-bdz*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_KINETIC_Hy4z_Py_vrr = PAY*I_KINETIC_G4z_Py_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_KINETIC_H5z_Py_vrr = PAZ*I_KINETIC_G4z_Py_vrr+4*oned2z*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Py_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_KINETIC_H5x_Pz_vrr = PAX*I_KINETIC_G4x_Pz_vrr+4*oned2z*I_KINETIC_F3x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Pz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_KINETIC_H4xy_Pz_vrr = PAY*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_KINETIC_H4xz_Pz_vrr = PAZ*I_KINETIC_G4x_Pz_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Pz_vrr;
    Double I_KINETIC_H3x2y_Pz_vrr = PAY*I_KINETIC_G3xy_Pz_vrr+oned2z*I_KINETIC_F3x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Pz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_KINETIC_H3x2z_Pz_vrr = PAZ*I_KINETIC_G3xz_Pz_vrr+oned2z*I_KINETIC_F3x_Pz_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_KINETIC_H2x3y_Pz_vrr = PAX*I_KINETIC_Gx3y_Pz_vrr+oned2z*I_KINETIC_F3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Pz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_KINETIC_H2x2yz_Pz_vrr = PAZ*I_KINETIC_G2x2y_Pz_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr;
    Double I_KINETIC_H2x3z_Pz_vrr = PAX*I_KINETIC_Gx3z_Pz_vrr+oned2z*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Pz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_KINETIC_Hx4y_Pz_vrr = PAX*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_KINETIC_Hx4z_Pz_vrr = PAX*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_KINETIC_H5y_Pz_vrr = PAY*I_KINETIC_G4y_Pz_vrr+4*oned2z*I_KINETIC_F3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Pz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_KINETIC_H4yz_Pz_vrr = PAZ*I_KINETIC_G4y_Pz_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Pz_vrr;
    Double I_KINETIC_H3y2z_Pz_vrr = PAZ*I_KINETIC_G3yz_Pz_vrr+oned2z*I_KINETIC_F3y_Pz_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_KINETIC_H2y3z_Pz_vrr = PAY*I_KINETIC_Gy3z_Pz_vrr+oned2z*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Pz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_KINETIC_Hy4z_Pz_vrr = PAY*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Pz_vrr;
    Double I_KINETIC_H5z_Pz_vrr = PAZ*I_KINETIC_G4z_Pz_vrr+4*oned2z*I_KINETIC_F3z_Pz_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Pz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_D
     * RHS shell quartet name: SQ_KINETIC_F_D
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_KINETIC_G4x_F3x_vrr = PBX*I_KINETIC_G4x_D2x_vrr+4*oned2z*I_KINETIC_F3x_D2x_vrr+2*oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_KINETIC_G3xy_F3x_vrr = PBX*I_KINETIC_G3xy_D2x_vrr+3*oned2z*I_KINETIC_F2xy_D2x_vrr+2*oned2z*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_KINETIC_G3xz_F3x_vrr = PBX*I_KINETIC_G3xz_D2x_vrr+3*oned2z*I_KINETIC_F2xz_D2x_vrr+2*oned2z*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_KINETIC_G2x2y_F3x_vrr = PBX*I_KINETIC_G2x2y_D2x_vrr+2*oned2z*I_KINETIC_Fx2y_D2x_vrr+2*oned2z*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_KINETIC_G2xyz_F3x_vrr = PBX*I_KINETIC_G2xyz_D2x_vrr+2*oned2z*I_KINETIC_Fxyz_D2x_vrr+2*oned2z*I_KINETIC_G2xyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_KINETIC_G2x2z_F3x_vrr = PBX*I_KINETIC_G2x2z_D2x_vrr+2*oned2z*I_KINETIC_Fx2z_D2x_vrr+2*oned2z*I_KINETIC_G2x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_KINETIC_Gx3y_F3x_vrr = PBX*I_KINETIC_Gx3y_D2x_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+2*oned2z*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_KINETIC_Gx2yz_F3x_vrr = PBX*I_KINETIC_Gx2yz_D2x_vrr+oned2z*I_KINETIC_F2yz_D2x_vrr+2*oned2z*I_KINETIC_Gx2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_KINETIC_Gxy2z_F3x_vrr = PBX*I_KINETIC_Gxy2z_D2x_vrr+oned2z*I_KINETIC_Fy2z_D2x_vrr+2*oned2z*I_KINETIC_Gxy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_KINETIC_Gx3z_F3x_vrr = PBX*I_KINETIC_Gx3z_D2x_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+2*oned2z*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_KINETIC_G4y_F3x_vrr = PBX*I_KINETIC_G4y_D2x_vrr+2*oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_KINETIC_G3yz_F3x_vrr = PBX*I_KINETIC_G3yz_D2x_vrr+2*oned2z*I_KINETIC_G3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_KINETIC_G2y2z_F3x_vrr = PBX*I_KINETIC_G2y2z_D2x_vrr+2*oned2z*I_KINETIC_G2y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_KINETIC_Gy3z_F3x_vrr = PBX*I_KINETIC_Gy3z_D2x_vrr+2*oned2z*I_KINETIC_Gy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_KINETIC_G4z_F3x_vrr = PBX*I_KINETIC_G4z_D2x_vrr+2*oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_KINETIC_G4x_F2xy_vrr = PBY*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_KINETIC_G3xy_F2xy_vrr = PBY*I_KINETIC_G3xy_D2x_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2xy_vrr;
    Double I_KINETIC_G3xz_F2xy_vrr = PBY*I_KINETIC_G3xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2xy_vrr;
    Double I_KINETIC_G2x2y_F2xy_vrr = PBY*I_KINETIC_G2x2y_D2x_vrr+2*oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr;
    Double I_KINETIC_G2xyz_F2xy_vrr = PBY*I_KINETIC_G2xyz_D2x_vrr+oned2z*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F2xy_vrr;
    Double I_KINETIC_G2x2z_F2xy_vrr = PBY*I_KINETIC_G2x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2xy_vrr;
    Double I_KINETIC_Gx3y_F2xy_vrr = PBY*I_KINETIC_Gx3y_D2x_vrr+3*oned2z*I_KINETIC_Fx2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr;
    Double I_KINETIC_Gx2yz_F2xy_vrr = PBY*I_KINETIC_Gx2yz_D2x_vrr+2*oned2z*I_KINETIC_Fxyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F2xy_vrr;
    Double I_KINETIC_Gxy2z_F2xy_vrr = PBY*I_KINETIC_Gxy2z_D2x_vrr+oned2z*I_KINETIC_Fx2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F2xy_vrr;
    Double I_KINETIC_Gx3z_F2xy_vrr = PBY*I_KINETIC_Gx3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr;
    Double I_KINETIC_G4y_F2xy_vrr = PBY*I_KINETIC_G4y_D2x_vrr+4*oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_KINETIC_G3yz_F2xy_vrr = PBY*I_KINETIC_G3yz_D2x_vrr+3*oned2z*I_KINETIC_F2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2xy_vrr;
    Double I_KINETIC_G2y2z_F2xy_vrr = PBY*I_KINETIC_G2y2z_D2x_vrr+2*oned2z*I_KINETIC_Fy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2xy_vrr;
    Double I_KINETIC_Gy3z_F2xy_vrr = PBY*I_KINETIC_Gy3z_D2x_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2xy_vrr;
    Double I_KINETIC_G4z_F2xy_vrr = PBY*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_KINETIC_G4x_F2xz_vrr = PBZ*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_KINETIC_G3xy_F2xz_vrr = PBZ*I_KINETIC_G3xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2xz_vrr;
    Double I_KINETIC_G3xz_F2xz_vrr = PBZ*I_KINETIC_G3xz_D2x_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2xz_vrr;
    Double I_KINETIC_G2x2y_F2xz_vrr = PBZ*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2xz_vrr;
    Double I_KINETIC_G2xyz_F2xz_vrr = PBZ*I_KINETIC_G2xyz_D2x_vrr+oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F2xz_vrr;
    Double I_KINETIC_G2x2z_F2xz_vrr = PBZ*I_KINETIC_G2x2z_D2x_vrr+2*oned2z*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2xz_vrr;
    Double I_KINETIC_Gx3y_F2xz_vrr = PBZ*I_KINETIC_Gx3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2xz_vrr;
    Double I_KINETIC_Gx2yz_F2xz_vrr = PBZ*I_KINETIC_Gx2yz_D2x_vrr+oned2z*I_KINETIC_Fx2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F2xz_vrr;
    Double I_KINETIC_Gxy2z_F2xz_vrr = PBZ*I_KINETIC_Gxy2z_D2x_vrr+2*oned2z*I_KINETIC_Fxyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F2xz_vrr;
    Double I_KINETIC_Gx3z_F2xz_vrr = PBZ*I_KINETIC_Gx3z_D2x_vrr+3*oned2z*I_KINETIC_Fx2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2xz_vrr;
    Double I_KINETIC_G4y_F2xz_vrr = PBZ*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2xz_vrr;
    Double I_KINETIC_G3yz_F2xz_vrr = PBZ*I_KINETIC_G3yz_D2x_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2xz_vrr;
    Double I_KINETIC_G2y2z_F2xz_vrr = PBZ*I_KINETIC_G2y2z_D2x_vrr+2*oned2z*I_KINETIC_F2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2xz_vrr;
    Double I_KINETIC_Gy3z_F2xz_vrr = PBZ*I_KINETIC_Gy3z_D2x_vrr+3*oned2z*I_KINETIC_Fy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2xz_vrr;
    Double I_KINETIC_G4z_F2xz_vrr = PBZ*I_KINETIC_G4z_D2x_vrr+4*oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_KINETIC_G4x_Fx2y_vrr = PBX*I_KINETIC_G4x_D2y_vrr+4*oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_KINETIC_G3xy_Fx2y_vrr = PBX*I_KINETIC_G3xy_D2y_vrr+3*oned2z*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr;
    Double I_KINETIC_G3xz_Fx2y_vrr = PBX*I_KINETIC_G3xz_D2y_vrr+3*oned2z*I_KINETIC_F2xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fx2y_vrr;
    Double I_KINETIC_G2x2y_Fx2y_vrr = PBX*I_KINETIC_G2x2y_D2y_vrr+2*oned2z*I_KINETIC_Fx2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr;
    Double I_KINETIC_G2xyz_Fx2y_vrr = PBX*I_KINETIC_G2xyz_D2y_vrr+2*oned2z*I_KINETIC_Fxyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fx2y_vrr;
    Double I_KINETIC_G2x2z_Fx2y_vrr = PBX*I_KINETIC_G2x2z_D2y_vrr+2*oned2z*I_KINETIC_Fx2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr;
    Double I_KINETIC_Gx3y_Fx2y_vrr = PBX*I_KINETIC_Gx3y_D2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr;
    Double I_KINETIC_Gx2yz_Fx2y_vrr = PBX*I_KINETIC_Gx2yz_D2y_vrr+oned2z*I_KINETIC_F2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fx2y_vrr;
    Double I_KINETIC_Gxy2z_Fx2y_vrr = PBX*I_KINETIC_Gxy2z_D2y_vrr+oned2z*I_KINETIC_Fy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fx2y_vrr;
    Double I_KINETIC_Gx3z_Fx2y_vrr = PBX*I_KINETIC_Gx3z_D2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr;
    Double I_KINETIC_G4y_Fx2y_vrr = PBX*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_KINETIC_G3yz_Fx2y_vrr = PBX*I_KINETIC_G3yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fx2y_vrr;
    Double I_KINETIC_G2y2z_Fx2y_vrr = PBX*I_KINETIC_G2y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr;
    Double I_KINETIC_Gy3z_Fx2y_vrr = PBX*I_KINETIC_Gy3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr;
    Double I_KINETIC_G4z_Fx2y_vrr = PBX*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_KINETIC_G4x_Fxyz_vrr = PBZ*I_KINETIC_G4x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fxyz_vrr;
    Double I_KINETIC_G3xy_Fxyz_vrr = PBZ*I_KINETIC_G3xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr;
    Double I_KINETIC_G3xz_Fxyz_vrr = PBZ*I_KINETIC_G3xz_Dxy_vrr+oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fxyz_vrr;
    Double I_KINETIC_G2x2y_Fxyz_vrr = PBZ*I_KINETIC_G2x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr;
    Double I_KINETIC_G2xyz_Fxyz_vrr = PBZ*I_KINETIC_G2xyz_Dxy_vrr+oned2z*I_KINETIC_F2xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fxyz_vrr;
    Double I_KINETIC_G2x2z_Fxyz_vrr = PBZ*I_KINETIC_G2x2z_Dxy_vrr+2*oned2z*I_KINETIC_F2xz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr;
    Double I_KINETIC_Gx3y_Fxyz_vrr = PBZ*I_KINETIC_Gx3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr;
    Double I_KINETIC_Gx2yz_Fxyz_vrr = PBZ*I_KINETIC_Gx2yz_Dxy_vrr+oned2z*I_KINETIC_Fx2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fxyz_vrr;
    Double I_KINETIC_Gxy2z_Fxyz_vrr = PBZ*I_KINETIC_Gxy2z_Dxy_vrr+2*oned2z*I_KINETIC_Fxyz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fxyz_vrr;
    Double I_KINETIC_Gx3z_Fxyz_vrr = PBZ*I_KINETIC_Gx3z_Dxy_vrr+3*oned2z*I_KINETIC_Fx2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr;
    Double I_KINETIC_G4y_Fxyz_vrr = PBZ*I_KINETIC_G4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fxyz_vrr;
    Double I_KINETIC_G3yz_Fxyz_vrr = PBZ*I_KINETIC_G3yz_Dxy_vrr+oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fxyz_vrr;
    Double I_KINETIC_G2y2z_Fxyz_vrr = PBZ*I_KINETIC_G2y2z_Dxy_vrr+2*oned2z*I_KINETIC_F2yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr;
    Double I_KINETIC_Gy3z_Fxyz_vrr = PBZ*I_KINETIC_Gy3z_Dxy_vrr+3*oned2z*I_KINETIC_Fy2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr;
    Double I_KINETIC_G4z_Fxyz_vrr = PBZ*I_KINETIC_G4z_Dxy_vrr+4*oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fxyz_vrr;
    Double I_KINETIC_G4x_Fx2z_vrr = PBX*I_KINETIC_G4x_D2z_vrr+4*oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_KINETIC_G3xy_Fx2z_vrr = PBX*I_KINETIC_G3xy_D2z_vrr+3*oned2z*I_KINETIC_F2xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr;
    Double I_KINETIC_G3xz_Fx2z_vrr = PBX*I_KINETIC_G3xz_D2z_vrr+3*oned2z*I_KINETIC_F2xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fx2z_vrr;
    Double I_KINETIC_G2x2y_Fx2z_vrr = PBX*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_Fx2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr;
    Double I_KINETIC_G2xyz_Fx2z_vrr = PBX*I_KINETIC_G2xyz_D2z_vrr+2*oned2z*I_KINETIC_Fxyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fx2z_vrr;
    Double I_KINETIC_G2x2z_Fx2z_vrr = PBX*I_KINETIC_G2x2z_D2z_vrr+2*oned2z*I_KINETIC_Fx2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr;
    Double I_KINETIC_Gx3y_Fx2z_vrr = PBX*I_KINETIC_Gx3y_D2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr;
    Double I_KINETIC_Gx2yz_Fx2z_vrr = PBX*I_KINETIC_Gx2yz_D2z_vrr+oned2z*I_KINETIC_F2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fx2z_vrr;
    Double I_KINETIC_Gxy2z_Fx2z_vrr = PBX*I_KINETIC_Gxy2z_D2z_vrr+oned2z*I_KINETIC_Fy2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fx2z_vrr;
    Double I_KINETIC_Gx3z_Fx2z_vrr = PBX*I_KINETIC_Gx3z_D2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr;
    Double I_KINETIC_G4y_Fx2z_vrr = PBX*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_KINETIC_G3yz_Fx2z_vrr = PBX*I_KINETIC_G3yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fx2z_vrr;
    Double I_KINETIC_G2y2z_Fx2z_vrr = PBX*I_KINETIC_G2y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr;
    Double I_KINETIC_Gy3z_Fx2z_vrr = PBX*I_KINETIC_Gy3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr;
    Double I_KINETIC_G4z_Fx2z_vrr = PBX*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_KINETIC_G4x_F3y_vrr = PBY*I_KINETIC_G4x_D2y_vrr+2*oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_KINETIC_G3xy_F3y_vrr = PBY*I_KINETIC_G3xy_D2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+2*oned2z*I_KINETIC_G3xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_KINETIC_G3xz_F3y_vrr = PBY*I_KINETIC_G3xz_D2y_vrr+2*oned2z*I_KINETIC_G3xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_KINETIC_G2x2y_F3y_vrr = PBY*I_KINETIC_G2x2y_D2y_vrr+2*oned2z*I_KINETIC_F2xy_D2y_vrr+2*oned2z*I_KINETIC_G2x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_KINETIC_G2xyz_F3y_vrr = PBY*I_KINETIC_G2xyz_D2y_vrr+oned2z*I_KINETIC_F2xz_D2y_vrr+2*oned2z*I_KINETIC_G2xyz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_KINETIC_G2x2z_F3y_vrr = PBY*I_KINETIC_G2x2z_D2y_vrr+2*oned2z*I_KINETIC_G2x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_KINETIC_Gx3y_F3y_vrr = PBY*I_KINETIC_Gx3y_D2y_vrr+3*oned2z*I_KINETIC_Fx2y_D2y_vrr+2*oned2z*I_KINETIC_Gx3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_KINETIC_Gx2yz_F3y_vrr = PBY*I_KINETIC_Gx2yz_D2y_vrr+2*oned2z*I_KINETIC_Fxyz_D2y_vrr+2*oned2z*I_KINETIC_Gx2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_KINETIC_Gxy2z_F3y_vrr = PBY*I_KINETIC_Gxy2z_D2y_vrr+oned2z*I_KINETIC_Fx2z_D2y_vrr+2*oned2z*I_KINETIC_Gxy2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_KINETIC_Gx3z_F3y_vrr = PBY*I_KINETIC_Gx3z_D2y_vrr+2*oned2z*I_KINETIC_Gx3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_KINETIC_G4y_F3y_vrr = PBY*I_KINETIC_G4y_D2y_vrr+4*oned2z*I_KINETIC_F3y_D2y_vrr+2*oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_KINETIC_G3yz_F3y_vrr = PBY*I_KINETIC_G3yz_D2y_vrr+3*oned2z*I_KINETIC_F2yz_D2y_vrr+2*oned2z*I_KINETIC_G3yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_KINETIC_G2y2z_F3y_vrr = PBY*I_KINETIC_G2y2z_D2y_vrr+2*oned2z*I_KINETIC_Fy2z_D2y_vrr+2*oned2z*I_KINETIC_G2y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_KINETIC_Gy3z_F3y_vrr = PBY*I_KINETIC_Gy3z_D2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+2*oned2z*I_KINETIC_Gy3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_KINETIC_G4z_F3y_vrr = PBY*I_KINETIC_G4z_D2y_vrr+2*oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_KINETIC_G4x_F2yz_vrr = PBZ*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2yz_vrr;
    Double I_KINETIC_G3xy_F2yz_vrr = PBZ*I_KINETIC_G3xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2yz_vrr;
    Double I_KINETIC_G3xz_F2yz_vrr = PBZ*I_KINETIC_G3xz_D2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2yz_vrr;
    Double I_KINETIC_G2x2y_F2yz_vrr = PBZ*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2yz_vrr;
    Double I_KINETIC_G2xyz_F2yz_vrr = PBZ*I_KINETIC_G2xyz_D2y_vrr+oned2z*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F2yz_vrr;
    Double I_KINETIC_G2x2z_F2yz_vrr = PBZ*I_KINETIC_G2x2z_D2y_vrr+2*oned2z*I_KINETIC_F2xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2yz_vrr;
    Double I_KINETIC_Gx3y_F2yz_vrr = PBZ*I_KINETIC_Gx3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2yz_vrr;
    Double I_KINETIC_Gx2yz_F2yz_vrr = PBZ*I_KINETIC_Gx2yz_D2y_vrr+oned2z*I_KINETIC_Fx2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F2yz_vrr;
    Double I_KINETIC_Gxy2z_F2yz_vrr = PBZ*I_KINETIC_Gxy2z_D2y_vrr+2*oned2z*I_KINETIC_Fxyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F2yz_vrr;
    Double I_KINETIC_Gx3z_F2yz_vrr = PBZ*I_KINETIC_Gx3z_D2y_vrr+3*oned2z*I_KINETIC_Fx2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2yz_vrr;
    Double I_KINETIC_G4y_F2yz_vrr = PBZ*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_KINETIC_G3yz_F2yz_vrr = PBZ*I_KINETIC_G3yz_D2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2yz_vrr;
    Double I_KINETIC_G2y2z_F2yz_vrr = PBZ*I_KINETIC_G2y2z_D2y_vrr+2*oned2z*I_KINETIC_F2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2yz_vrr;
    Double I_KINETIC_Gy3z_F2yz_vrr = PBZ*I_KINETIC_Gy3z_D2y_vrr+3*oned2z*I_KINETIC_Fy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2yz_vrr;
    Double I_KINETIC_G4z_F2yz_vrr = PBZ*I_KINETIC_G4z_D2y_vrr+4*oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_KINETIC_G4x_Fy2z_vrr = PBY*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fy2z_vrr;
    Double I_KINETIC_G3xy_Fy2z_vrr = PBY*I_KINETIC_G3xy_D2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr;
    Double I_KINETIC_G3xz_Fy2z_vrr = PBY*I_KINETIC_G3xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fy2z_vrr;
    Double I_KINETIC_G2x2y_Fy2z_vrr = PBY*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_F2xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr;
    Double I_KINETIC_G2xyz_Fy2z_vrr = PBY*I_KINETIC_G2xyz_D2z_vrr+oned2z*I_KINETIC_F2xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fy2z_vrr;
    Double I_KINETIC_G2x2z_Fy2z_vrr = PBY*I_KINETIC_G2x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr;
    Double I_KINETIC_Gx3y_Fy2z_vrr = PBY*I_KINETIC_Gx3y_D2z_vrr+3*oned2z*I_KINETIC_Fx2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr;
    Double I_KINETIC_Gx2yz_Fy2z_vrr = PBY*I_KINETIC_Gx2yz_D2z_vrr+2*oned2z*I_KINETIC_Fxyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fy2z_vrr;
    Double I_KINETIC_Gxy2z_Fy2z_vrr = PBY*I_KINETIC_Gxy2z_D2z_vrr+oned2z*I_KINETIC_Fx2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fy2z_vrr;
    Double I_KINETIC_Gx3z_Fy2z_vrr = PBY*I_KINETIC_Gx3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr;
    Double I_KINETIC_G4y_Fy2z_vrr = PBY*I_KINETIC_G4y_D2z_vrr+4*oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_KINETIC_G3yz_Fy2z_vrr = PBY*I_KINETIC_G3yz_D2z_vrr+3*oned2z*I_KINETIC_F2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fy2z_vrr;
    Double I_KINETIC_G2y2z_Fy2z_vrr = PBY*I_KINETIC_G2y2z_D2z_vrr+2*oned2z*I_KINETIC_Fy2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr;
    Double I_KINETIC_Gy3z_Fy2z_vrr = PBY*I_KINETIC_Gy3z_D2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr;
    Double I_KINETIC_G4z_Fy2z_vrr = PBY*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_KINETIC_G4x_F3z_vrr = PBZ*I_KINETIC_G4x_D2z_vrr+2*oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_KINETIC_G3xy_F3z_vrr = PBZ*I_KINETIC_G3xy_D2z_vrr+2*oned2z*I_KINETIC_G3xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_KINETIC_G3xz_F3z_vrr = PBZ*I_KINETIC_G3xz_D2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+2*oned2z*I_KINETIC_G3xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_KINETIC_G2x2y_F3z_vrr = PBZ*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_G2x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_KINETIC_G2xyz_F3z_vrr = PBZ*I_KINETIC_G2xyz_D2z_vrr+oned2z*I_KINETIC_F2xy_D2z_vrr+2*oned2z*I_KINETIC_G2xyz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G2xyz_Pz_vrr;
    Double I_KINETIC_G2x2z_F3z_vrr = PBZ*I_KINETIC_G2x2z_D2z_vrr+2*oned2z*I_KINETIC_F2xz_D2z_vrr+2*oned2z*I_KINETIC_G2x2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_KINETIC_Gx3y_F3z_vrr = PBZ*I_KINETIC_Gx3y_D2z_vrr+2*oned2z*I_KINETIC_Gx3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_KINETIC_Gx2yz_F3z_vrr = PBZ*I_KINETIC_Gx2yz_D2z_vrr+oned2z*I_KINETIC_Fx2y_D2z_vrr+2*oned2z*I_KINETIC_Gx2yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr;
    Double I_KINETIC_Gxy2z_F3z_vrr = PBZ*I_KINETIC_Gxy2z_D2z_vrr+2*oned2z*I_KINETIC_Fxyz_D2z_vrr+2*oned2z*I_KINETIC_Gxy2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr;
    Double I_KINETIC_Gx3z_F3z_vrr = PBZ*I_KINETIC_Gx3z_D2z_vrr+3*oned2z*I_KINETIC_Fx2z_D2z_vrr+2*oned2z*I_KINETIC_Gx3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_KINETIC_G4y_F3z_vrr = PBZ*I_KINETIC_G4y_D2z_vrr+2*oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_KINETIC_G3yz_F3z_vrr = PBZ*I_KINETIC_G3yz_D2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+2*oned2z*I_KINETIC_G3yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_KINETIC_G2y2z_F3z_vrr = PBZ*I_KINETIC_G2y2z_D2z_vrr+2*oned2z*I_KINETIC_F2yz_D2z_vrr+2*oned2z*I_KINETIC_G2y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_KINETIC_Gy3z_F3z_vrr = PBZ*I_KINETIC_Gy3z_D2z_vrr+3*oned2z*I_KINETIC_Fy2z_D2z_vrr+2*oned2z*I_KINETIC_Gy3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_KINETIC_G4z_F3z_vrr = PBZ*I_KINETIC_G4z_D2z_vrr+4*oned2z*I_KINETIC_F3z_D2z_vrr+2*oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_D
     * RHS shell quartet name: SQ_KINETIC_F_D
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     ************************************************************/
    Double I_KINETIC_H5x_D2x_vrr = PAX*I_KINETIC_G4x_D2x_vrr+4*oned2z*I_KINETIC_F3x_D2x_vrr+2*oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2x_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_KINETIC_H4xy_D2x_vrr = PAY*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_KINETIC_H4xz_D2x_vrr = PAZ*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_KINETIC_H3x2y_D2x_vrr = PAY*I_KINETIC_G3xy_D2x_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_KINETIC_H3x2z_D2x_vrr = PAZ*I_KINETIC_G3xz_D2x_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_KINETIC_H2x3y_D2x_vrr = PAX*I_KINETIC_Gx3y_D2x_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+2*oned2z*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_KINETIC_H2x2yz_D2x_vrr = PAZ*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr;
    Double I_KINETIC_H2x3z_D2x_vrr = PAX*I_KINETIC_Gx3z_D2x_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+2*oned2z*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_KINETIC_Hx4y_D2x_vrr = PAX*I_KINETIC_G4y_D2x_vrr+2*oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_KINETIC_Hx4z_D2x_vrr = PAX*I_KINETIC_G4z_D2x_vrr+2*oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_KINETIC_H5y_D2x_vrr = PAY*I_KINETIC_G4y_D2x_vrr+4*oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2x_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_KINETIC_H4yz_D2x_vrr = PAZ*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_KINETIC_H3y2z_D2x_vrr = PAZ*I_KINETIC_G3yz_D2x_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_KINETIC_H2y3z_D2x_vrr = PAY*I_KINETIC_Gy3z_D2x_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_KINETIC_Hy4z_D2x_vrr = PAY*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_KINETIC_H5z_D2x_vrr = PAZ*I_KINETIC_G4z_D2x_vrr+4*oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2x_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_KINETIC_H5x_Dxy_vrr = PAX*I_KINETIC_G4x_Dxy_vrr+4*oned2z*I_KINETIC_F3x_Dxy_vrr+oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dxy_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_KINETIC_H4xy_Dxy_vrr = PAY*I_KINETIC_G4x_Dxy_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dxy_vrr;
    Double I_KINETIC_H4xz_Dxy_vrr = PAZ*I_KINETIC_G4x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dxy_vrr;
    Double I_KINETIC_H3x2y_Dxy_vrr = PAY*I_KINETIC_G3xy_Dxy_vrr+oned2z*I_KINETIC_F3x_Dxy_vrr+oned2z*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_KINETIC_H3x2z_Dxy_vrr = PAZ*I_KINETIC_G3xz_Dxy_vrr+oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_KINETIC_H2x3y_Dxy_vrr = PAX*I_KINETIC_Gx3y_Dxy_vrr+oned2z*I_KINETIC_F3y_Dxy_vrr+oned2z*I_KINETIC_Gx3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_KINETIC_H2x2yz_Dxy_vrr = PAZ*I_KINETIC_G2x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr;
    Double I_KINETIC_H2x3z_Dxy_vrr = PAX*I_KINETIC_Gx3z_Dxy_vrr+oned2z*I_KINETIC_F3z_Dxy_vrr+oned2z*I_KINETIC_Gx3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_KINETIC_Hx4y_Dxy_vrr = PAX*I_KINETIC_G4y_Dxy_vrr+oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr;
    Double I_KINETIC_Hx4z_Dxy_vrr = PAX*I_KINETIC_G4z_Dxy_vrr+oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr;
    Double I_KINETIC_H5y_Dxy_vrr = PAY*I_KINETIC_G4y_Dxy_vrr+4*oned2z*I_KINETIC_F3y_Dxy_vrr+oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dxy_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_KINETIC_H4yz_Dxy_vrr = PAZ*I_KINETIC_G4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dxy_vrr;
    Double I_KINETIC_H3y2z_Dxy_vrr = PAZ*I_KINETIC_G3yz_Dxy_vrr+oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_KINETIC_H2y3z_Dxy_vrr = PAY*I_KINETIC_Gy3z_Dxy_vrr+oned2z*I_KINETIC_F3z_Dxy_vrr+oned2z*I_KINETIC_Gy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_KINETIC_Hy4z_Dxy_vrr = PAY*I_KINETIC_G4z_Dxy_vrr+oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dxy_vrr;
    Double I_KINETIC_H5z_Dxy_vrr = PAZ*I_KINETIC_G4z_Dxy_vrr+4*oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dxy_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_KINETIC_H5x_Dxz_vrr = PAX*I_KINETIC_G4x_Dxz_vrr+4*oned2z*I_KINETIC_F3x_Dxz_vrr+oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dxz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_KINETIC_H4xy_Dxz_vrr = PAY*I_KINETIC_G4x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dxz_vrr;
    Double I_KINETIC_H4xz_Dxz_vrr = PAZ*I_KINETIC_G4x_Dxz_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dxz_vrr;
    Double I_KINETIC_H3x2y_Dxz_vrr = PAY*I_KINETIC_G3xy_Dxz_vrr+oned2z*I_KINETIC_F3x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_KINETIC_H3x2z_Dxz_vrr = PAZ*I_KINETIC_G3xz_Dxz_vrr+oned2z*I_KINETIC_F3x_Dxz_vrr+oned2z*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_KINETIC_H2x3y_Dxz_vrr = PAX*I_KINETIC_Gx3y_Dxz_vrr+oned2z*I_KINETIC_F3y_Dxz_vrr+oned2z*I_KINETIC_Gx3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_KINETIC_H2x2yz_Dxz_vrr = PAZ*I_KINETIC_G2x2y_Dxz_vrr+oned2z*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr;
    Double I_KINETIC_H2x3z_Dxz_vrr = PAX*I_KINETIC_Gx3z_Dxz_vrr+oned2z*I_KINETIC_F3z_Dxz_vrr+oned2z*I_KINETIC_Gx3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_KINETIC_Hx4y_Dxz_vrr = PAX*I_KINETIC_G4y_Dxz_vrr+oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr;
    Double I_KINETIC_Hx4z_Dxz_vrr = PAX*I_KINETIC_G4z_Dxz_vrr+oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr;
    Double I_KINETIC_H5y_Dxz_vrr = PAY*I_KINETIC_G4y_Dxz_vrr+4*oned2z*I_KINETIC_F3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dxz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_KINETIC_H4yz_Dxz_vrr = PAZ*I_KINETIC_G4y_Dxz_vrr+oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dxz_vrr;
    Double I_KINETIC_H3y2z_Dxz_vrr = PAZ*I_KINETIC_G3yz_Dxz_vrr+oned2z*I_KINETIC_F3y_Dxz_vrr+oned2z*I_KINETIC_G3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_KINETIC_H2y3z_Dxz_vrr = PAY*I_KINETIC_Gy3z_Dxz_vrr+oned2z*I_KINETIC_F3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_KINETIC_Hy4z_Dxz_vrr = PAY*I_KINETIC_G4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dxz_vrr;
    Double I_KINETIC_H5z_Dxz_vrr = PAZ*I_KINETIC_G4z_Dxz_vrr+4*oned2z*I_KINETIC_F3z_Dxz_vrr+oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dxz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_KINETIC_H5x_D2y_vrr = PAX*I_KINETIC_G4x_D2y_vrr+4*oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_KINETIC_H4xy_D2y_vrr = PAY*I_KINETIC_G4x_D2y_vrr+2*oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_KINETIC_H4xz_D2y_vrr = PAZ*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2y_vrr;
    Double I_KINETIC_H3x2y_D2y_vrr = PAY*I_KINETIC_G3xy_D2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+2*oned2z*I_KINETIC_G3xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_KINETIC_H3x2z_D2y_vrr = PAZ*I_KINETIC_G3xz_D2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_KINETIC_H2x3y_D2y_vrr = PAX*I_KINETIC_Gx3y_D2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_KINETIC_H2x2yz_D2y_vrr = PAZ*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr;
    Double I_KINETIC_H2x3z_D2y_vrr = PAX*I_KINETIC_Gx3z_D2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_KINETIC_Hx4y_D2y_vrr = PAX*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_KINETIC_Hx4z_D2y_vrr = PAX*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_KINETIC_H5y_D2y_vrr = PAY*I_KINETIC_G4y_D2y_vrr+4*oned2z*I_KINETIC_F3y_D2y_vrr+2*oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_KINETIC_H4yz_D2y_vrr = PAZ*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_KINETIC_H3y2z_D2y_vrr = PAZ*I_KINETIC_G3yz_D2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_KINETIC_H2y3z_D2y_vrr = PAY*I_KINETIC_Gy3z_D2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+2*oned2z*I_KINETIC_Gy3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_KINETIC_Hy4z_D2y_vrr = PAY*I_KINETIC_G4z_D2y_vrr+2*oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_KINETIC_H5z_D2y_vrr = PAZ*I_KINETIC_G4z_D2y_vrr+4*oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_KINETIC_H5x_Dyz_vrr = PAX*I_KINETIC_G4x_Dyz_vrr+4*oned2z*I_KINETIC_F3x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_KINETIC_H4xy_Dyz_vrr = PAY*I_KINETIC_G4x_Dyz_vrr+oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dyz_vrr;
    Double I_KINETIC_H4xz_Dyz_vrr = PAZ*I_KINETIC_G4x_Dyz_vrr+oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dyz_vrr;
    Double I_KINETIC_H3x2y_Dyz_vrr = PAY*I_KINETIC_G3xy_Dyz_vrr+oned2z*I_KINETIC_F3x_Dyz_vrr+oned2z*I_KINETIC_G3xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_KINETIC_H3x2z_Dyz_vrr = PAZ*I_KINETIC_G3xz_Dyz_vrr+oned2z*I_KINETIC_F3x_Dyz_vrr+oned2z*I_KINETIC_G3xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_KINETIC_H2x3y_Dyz_vrr = PAX*I_KINETIC_Gx3y_Dyz_vrr+oned2z*I_KINETIC_F3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_KINETIC_H2x2yz_Dyz_vrr = PAZ*I_KINETIC_G2x2y_Dyz_vrr+oned2z*I_KINETIC_G2x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr;
    Double I_KINETIC_H2x3z_Dyz_vrr = PAX*I_KINETIC_Gx3z_Dyz_vrr+oned2z*I_KINETIC_F3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_KINETIC_Hx4y_Dyz_vrr = PAX*I_KINETIC_G4y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr;
    Double I_KINETIC_Hx4z_Dyz_vrr = PAX*I_KINETIC_G4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr;
    Double I_KINETIC_H5y_Dyz_vrr = PAY*I_KINETIC_G4y_Dyz_vrr+4*oned2z*I_KINETIC_F3y_Dyz_vrr+oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_KINETIC_H4yz_Dyz_vrr = PAZ*I_KINETIC_G4y_Dyz_vrr+oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dyz_vrr;
    Double I_KINETIC_H3y2z_Dyz_vrr = PAZ*I_KINETIC_G3yz_Dyz_vrr+oned2z*I_KINETIC_F3y_Dyz_vrr+oned2z*I_KINETIC_G3yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_KINETIC_H2y3z_Dyz_vrr = PAY*I_KINETIC_Gy3z_Dyz_vrr+oned2z*I_KINETIC_F3z_Dyz_vrr+oned2z*I_KINETIC_Gy3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_KINETIC_Hy4z_Dyz_vrr = PAY*I_KINETIC_G4z_Dyz_vrr+oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dyz_vrr;
    Double I_KINETIC_H5z_Dyz_vrr = PAZ*I_KINETIC_G4z_Dyz_vrr+4*oned2z*I_KINETIC_F3z_Dyz_vrr+oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_KINETIC_H5x_D2z_vrr = PAX*I_KINETIC_G4x_D2z_vrr+4*oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_KINETIC_H4xy_D2z_vrr = PAY*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2z_vrr;
    Double I_KINETIC_H4xz_D2z_vrr = PAZ*I_KINETIC_G4x_D2z_vrr+2*oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2z_vrr;
    Double I_KINETIC_H3x2y_D2z_vrr = PAY*I_KINETIC_G3xy_D2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_KINETIC_H3x2z_D2z_vrr = PAZ*I_KINETIC_G3xz_D2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+2*oned2z*I_KINETIC_G3xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_KINETIC_H2x3y_D2z_vrr = PAX*I_KINETIC_Gx3y_D2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_KINETIC_H2x2yz_D2z_vrr = PAZ*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_G2x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr;
    Double I_KINETIC_H2x3z_D2z_vrr = PAX*I_KINETIC_Gx3z_D2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_KINETIC_Hx4y_D2z_vrr = PAX*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2z_vrr;
    Double I_KINETIC_Hx4z_D2z_vrr = PAX*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_KINETIC_H5y_D2z_vrr = PAY*I_KINETIC_G4y_D2z_vrr+4*oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_KINETIC_H4yz_D2z_vrr = PAZ*I_KINETIC_G4y_D2z_vrr+2*oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2z_vrr;
    Double I_KINETIC_H3y2z_D2z_vrr = PAZ*I_KINETIC_G3yz_D2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+2*oned2z*I_KINETIC_G3yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_KINETIC_H2y3z_D2z_vrr = PAY*I_KINETIC_Gy3z_D2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_KINETIC_Hy4z_D2z_vrr = PAY*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2z_vrr;
    Double I_KINETIC_H5z_D2z_vrr = PAZ*I_KINETIC_G4z_D2z_vrr+4*oned2z*I_KINETIC_F3z_D2z_vrr+2*oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 50 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_D
     * RHS shell quartet name: SQ_KINETIC_G_D
     * RHS shell quartet name: SQ_KINETIC_H_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     ************************************************************/
    Double I_KINETIC_H5x_F3x_vrr = PBX*I_KINETIC_H5x_D2x_vrr+5*oned2z*I_KINETIC_G4x_D2x_vrr+2*oned2z*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_KINETIC_H4xy_F3x_vrr = PBX*I_KINETIC_H4xy_D2x_vrr+4*oned2z*I_KINETIC_G3xy_D2x_vrr+2*oned2z*I_KINETIC_H4xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_KINETIC_H4xz_F3x_vrr = PBX*I_KINETIC_H4xz_D2x_vrr+4*oned2z*I_KINETIC_G3xz_D2x_vrr+2*oned2z*I_KINETIC_H4xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_KINETIC_H3x2y_F3x_vrr = PBX*I_KINETIC_H3x2y_D2x_vrr+3*oned2z*I_KINETIC_G2x2y_D2x_vrr+2*oned2z*I_KINETIC_H3x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_KINETIC_H3x2z_F3x_vrr = PBX*I_KINETIC_H3x2z_D2x_vrr+3*oned2z*I_KINETIC_G2x2z_D2x_vrr+2*oned2z*I_KINETIC_H3x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_KINETIC_H2x3y_F3x_vrr = PBX*I_KINETIC_H2x3y_D2x_vrr+2*oned2z*I_KINETIC_Gx3y_D2x_vrr+2*oned2z*I_KINETIC_H2x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_KINETIC_H2x2yz_F3x_vrr = PBX*I_KINETIC_H2x2yz_D2x_vrr+2*oned2z*I_KINETIC_Gx2yz_D2x_vrr+2*oned2z*I_KINETIC_H2x2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_KINETIC_H2x3z_F3x_vrr = PBX*I_KINETIC_H2x3z_D2x_vrr+2*oned2z*I_KINETIC_Gx3z_D2x_vrr+2*oned2z*I_KINETIC_H2x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_KINETIC_Hx4y_F3x_vrr = PBX*I_KINETIC_Hx4y_D2x_vrr+oned2z*I_KINETIC_G4y_D2x_vrr+2*oned2z*I_KINETIC_Hx4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_KINETIC_Hx4z_F3x_vrr = PBX*I_KINETIC_Hx4z_D2x_vrr+oned2z*I_KINETIC_G4z_D2x_vrr+2*oned2z*I_KINETIC_Hx4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_KINETIC_H5y_F3x_vrr = PBX*I_KINETIC_H5y_D2x_vrr+2*oned2z*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_KINETIC_H4yz_F3x_vrr = PBX*I_KINETIC_H4yz_D2x_vrr+2*oned2z*I_KINETIC_H4yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_KINETIC_H3y2z_F3x_vrr = PBX*I_KINETIC_H3y2z_D2x_vrr+2*oned2z*I_KINETIC_H3y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_KINETIC_H2y3z_F3x_vrr = PBX*I_KINETIC_H2y3z_D2x_vrr+2*oned2z*I_KINETIC_H2y3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_KINETIC_Hy4z_F3x_vrr = PBX*I_KINETIC_Hy4z_D2x_vrr+2*oned2z*I_KINETIC_Hy4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_KINETIC_H5z_F3x_vrr = PBX*I_KINETIC_H5z_D2x_vrr+2*oned2z*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_KINETIC_H5x_F2xy_vrr = PBY*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F2xy_vrr;
    Double I_KINETIC_H4xy_F2xy_vrr = PBY*I_KINETIC_H4xy_D2x_vrr+oned2z*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F2xy_vrr;
    Double I_KINETIC_H4xz_F2xy_vrr = PBY*I_KINETIC_H4xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F2xy_vrr;
    Double I_KINETIC_H3x2y_F2xy_vrr = PBY*I_KINETIC_H3x2y_D2x_vrr+2*oned2z*I_KINETIC_G3xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F2xy_vrr;
    Double I_KINETIC_H3x2z_F2xy_vrr = PBY*I_KINETIC_H3x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F2xy_vrr;
    Double I_KINETIC_H2x3y_F2xy_vrr = PBY*I_KINETIC_H2x3y_D2x_vrr+3*oned2z*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F2xy_vrr;
    Double I_KINETIC_H2x2yz_F2xy_vrr = PBY*I_KINETIC_H2x2yz_D2x_vrr+2*oned2z*I_KINETIC_G2xyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F2xy_vrr;
    Double I_KINETIC_H2x3z_F2xy_vrr = PBY*I_KINETIC_H2x3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F2xy_vrr;
    Double I_KINETIC_Hx4y_F2xy_vrr = PBY*I_KINETIC_Hx4y_D2x_vrr+4*oned2z*I_KINETIC_Gx3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F2xy_vrr;
    Double I_KINETIC_Hx4z_F2xy_vrr = PBY*I_KINETIC_Hx4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F2xy_vrr;
    Double I_KINETIC_H5y_F2xy_vrr = PBY*I_KINETIC_H5y_D2x_vrr+5*oned2z*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F2xy_vrr;
    Double I_KINETIC_H4yz_F2xy_vrr = PBY*I_KINETIC_H4yz_D2x_vrr+4*oned2z*I_KINETIC_G3yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F2xy_vrr;
    Double I_KINETIC_H3y2z_F2xy_vrr = PBY*I_KINETIC_H3y2z_D2x_vrr+3*oned2z*I_KINETIC_G2y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F2xy_vrr;
    Double I_KINETIC_H2y3z_F2xy_vrr = PBY*I_KINETIC_H2y3z_D2x_vrr+2*oned2z*I_KINETIC_Gy3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F2xy_vrr;
    Double I_KINETIC_Hy4z_F2xy_vrr = PBY*I_KINETIC_Hy4z_D2x_vrr+oned2z*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F2xy_vrr;
    Double I_KINETIC_H5z_F2xy_vrr = PBY*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F2xy_vrr;
    Double I_KINETIC_H5x_F2xz_vrr = PBZ*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F2xz_vrr;
    Double I_KINETIC_H4xy_F2xz_vrr = PBZ*I_KINETIC_H4xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F2xz_vrr;
    Double I_KINETIC_H4xz_F2xz_vrr = PBZ*I_KINETIC_H4xz_D2x_vrr+oned2z*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F2xz_vrr;
    Double I_KINETIC_H3x2y_F2xz_vrr = PBZ*I_KINETIC_H3x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F2xz_vrr;
    Double I_KINETIC_H3x2z_F2xz_vrr = PBZ*I_KINETIC_H3x2z_D2x_vrr+2*oned2z*I_KINETIC_G3xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F2xz_vrr;
    Double I_KINETIC_H2x3y_F2xz_vrr = PBZ*I_KINETIC_H2x3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F2xz_vrr;
    Double I_KINETIC_H2x2yz_F2xz_vrr = PBZ*I_KINETIC_H2x2yz_D2x_vrr+oned2z*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F2xz_vrr;
    Double I_KINETIC_H2x3z_F2xz_vrr = PBZ*I_KINETIC_H2x3z_D2x_vrr+3*oned2z*I_KINETIC_G2x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F2xz_vrr;
    Double I_KINETIC_Hx4y_F2xz_vrr = PBZ*I_KINETIC_Hx4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F2xz_vrr;
    Double I_KINETIC_Hx4z_F2xz_vrr = PBZ*I_KINETIC_Hx4z_D2x_vrr+4*oned2z*I_KINETIC_Gx3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F2xz_vrr;
    Double I_KINETIC_H5y_F2xz_vrr = PBZ*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F2xz_vrr;
    Double I_KINETIC_H4yz_F2xz_vrr = PBZ*I_KINETIC_H4yz_D2x_vrr+oned2z*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F2xz_vrr;
    Double I_KINETIC_H3y2z_F2xz_vrr = PBZ*I_KINETIC_H3y2z_D2x_vrr+2*oned2z*I_KINETIC_G3yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F2xz_vrr;
    Double I_KINETIC_H2y3z_F2xz_vrr = PBZ*I_KINETIC_H2y3z_D2x_vrr+3*oned2z*I_KINETIC_G2y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F2xz_vrr;
    Double I_KINETIC_Hy4z_F2xz_vrr = PBZ*I_KINETIC_Hy4z_D2x_vrr+4*oned2z*I_KINETIC_Gy3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F2xz_vrr;
    Double I_KINETIC_H5z_F2xz_vrr = PBZ*I_KINETIC_H5z_D2x_vrr+5*oned2z*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F2xz_vrr;
    Double I_KINETIC_H5x_Fx2y_vrr = PBX*I_KINETIC_H5x_D2y_vrr+5*oned2z*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fx2y_vrr;
    Double I_KINETIC_H4xy_Fx2y_vrr = PBX*I_KINETIC_H4xy_D2y_vrr+4*oned2z*I_KINETIC_G3xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fx2y_vrr;
    Double I_KINETIC_H4xz_Fx2y_vrr = PBX*I_KINETIC_H4xz_D2y_vrr+4*oned2z*I_KINETIC_G3xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fx2y_vrr;
    Double I_KINETIC_H3x2y_Fx2y_vrr = PBX*I_KINETIC_H3x2y_D2y_vrr+3*oned2z*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fx2y_vrr;
    Double I_KINETIC_H3x2z_Fx2y_vrr = PBX*I_KINETIC_H3x2z_D2y_vrr+3*oned2z*I_KINETIC_G2x2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fx2y_vrr;
    Double I_KINETIC_H2x3y_Fx2y_vrr = PBX*I_KINETIC_H2x3y_D2y_vrr+2*oned2z*I_KINETIC_Gx3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fx2y_vrr;
    Double I_KINETIC_H2x2yz_Fx2y_vrr = PBX*I_KINETIC_H2x2yz_D2y_vrr+2*oned2z*I_KINETIC_Gx2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fx2y_vrr;
    Double I_KINETIC_H2x3z_Fx2y_vrr = PBX*I_KINETIC_H2x3z_D2y_vrr+2*oned2z*I_KINETIC_Gx3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fx2y_vrr;
    Double I_KINETIC_Hx4y_Fx2y_vrr = PBX*I_KINETIC_Hx4y_D2y_vrr+oned2z*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fx2y_vrr;
    Double I_KINETIC_Hx4z_Fx2y_vrr = PBX*I_KINETIC_Hx4z_D2y_vrr+oned2z*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fx2y_vrr;
    Double I_KINETIC_H5y_Fx2y_vrr = PBX*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fx2y_vrr;
    Double I_KINETIC_H4yz_Fx2y_vrr = PBX*I_KINETIC_H4yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fx2y_vrr;
    Double I_KINETIC_H3y2z_Fx2y_vrr = PBX*I_KINETIC_H3y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fx2y_vrr;
    Double I_KINETIC_H2y3z_Fx2y_vrr = PBX*I_KINETIC_H2y3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fx2y_vrr;
    Double I_KINETIC_Hy4z_Fx2y_vrr = PBX*I_KINETIC_Hy4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fx2y_vrr;
    Double I_KINETIC_H5z_Fx2y_vrr = PBX*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fx2y_vrr;
    Double I_KINETIC_H5x_Fxyz_vrr = PBZ*I_KINETIC_H5x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fxyz_vrr;
    Double I_KINETIC_H4xy_Fxyz_vrr = PBZ*I_KINETIC_H4xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fxyz_vrr;
    Double I_KINETIC_H4xz_Fxyz_vrr = PBZ*I_KINETIC_H4xz_Dxy_vrr+oned2z*I_KINETIC_G4x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fxyz_vrr;
    Double I_KINETIC_H3x2y_Fxyz_vrr = PBZ*I_KINETIC_H3x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fxyz_vrr;
    Double I_KINETIC_H3x2z_Fxyz_vrr = PBZ*I_KINETIC_H3x2z_Dxy_vrr+2*oned2z*I_KINETIC_G3xz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fxyz_vrr;
    Double I_KINETIC_H2x3y_Fxyz_vrr = PBZ*I_KINETIC_H2x3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fxyz_vrr;
    Double I_KINETIC_H2x2yz_Fxyz_vrr = PBZ*I_KINETIC_H2x2yz_Dxy_vrr+oned2z*I_KINETIC_G2x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fxyz_vrr;
    Double I_KINETIC_H2x3z_Fxyz_vrr = PBZ*I_KINETIC_H2x3z_Dxy_vrr+3*oned2z*I_KINETIC_G2x2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fxyz_vrr;
    Double I_KINETIC_Hx4y_Fxyz_vrr = PBZ*I_KINETIC_Hx4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fxyz_vrr;
    Double I_KINETIC_Hx4z_Fxyz_vrr = PBZ*I_KINETIC_Hx4z_Dxy_vrr+4*oned2z*I_KINETIC_Gx3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fxyz_vrr;
    Double I_KINETIC_H5y_Fxyz_vrr = PBZ*I_KINETIC_H5y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fxyz_vrr;
    Double I_KINETIC_H4yz_Fxyz_vrr = PBZ*I_KINETIC_H4yz_Dxy_vrr+oned2z*I_KINETIC_G4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fxyz_vrr;
    Double I_KINETIC_H3y2z_Fxyz_vrr = PBZ*I_KINETIC_H3y2z_Dxy_vrr+2*oned2z*I_KINETIC_G3yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fxyz_vrr;
    Double I_KINETIC_H2y3z_Fxyz_vrr = PBZ*I_KINETIC_H2y3z_Dxy_vrr+3*oned2z*I_KINETIC_G2y2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fxyz_vrr;
    Double I_KINETIC_Hy4z_Fxyz_vrr = PBZ*I_KINETIC_Hy4z_Dxy_vrr+4*oned2z*I_KINETIC_Gy3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fxyz_vrr;
    Double I_KINETIC_H5z_Fxyz_vrr = PBZ*I_KINETIC_H5z_Dxy_vrr+5*oned2z*I_KINETIC_G4z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fxyz_vrr;
    Double I_KINETIC_H5x_Fx2z_vrr = PBX*I_KINETIC_H5x_D2z_vrr+5*oned2z*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fx2z_vrr;
    Double I_KINETIC_H4xy_Fx2z_vrr = PBX*I_KINETIC_H4xy_D2z_vrr+4*oned2z*I_KINETIC_G3xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fx2z_vrr;
    Double I_KINETIC_H4xz_Fx2z_vrr = PBX*I_KINETIC_H4xz_D2z_vrr+4*oned2z*I_KINETIC_G3xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fx2z_vrr;
    Double I_KINETIC_H3x2y_Fx2z_vrr = PBX*I_KINETIC_H3x2y_D2z_vrr+3*oned2z*I_KINETIC_G2x2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fx2z_vrr;
    Double I_KINETIC_H3x2z_Fx2z_vrr = PBX*I_KINETIC_H3x2z_D2z_vrr+3*oned2z*I_KINETIC_G2x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fx2z_vrr;
    Double I_KINETIC_H2x3y_Fx2z_vrr = PBX*I_KINETIC_H2x3y_D2z_vrr+2*oned2z*I_KINETIC_Gx3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fx2z_vrr;
    Double I_KINETIC_H2x2yz_Fx2z_vrr = PBX*I_KINETIC_H2x2yz_D2z_vrr+2*oned2z*I_KINETIC_Gx2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fx2z_vrr;
    Double I_KINETIC_H2x3z_Fx2z_vrr = PBX*I_KINETIC_H2x3z_D2z_vrr+2*oned2z*I_KINETIC_Gx3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fx2z_vrr;
    Double I_KINETIC_Hx4y_Fx2z_vrr = PBX*I_KINETIC_Hx4y_D2z_vrr+oned2z*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fx2z_vrr;
    Double I_KINETIC_Hx4z_Fx2z_vrr = PBX*I_KINETIC_Hx4z_D2z_vrr+oned2z*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fx2z_vrr;
    Double I_KINETIC_H5y_Fx2z_vrr = PBX*I_KINETIC_H5y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fx2z_vrr;
    Double I_KINETIC_H4yz_Fx2z_vrr = PBX*I_KINETIC_H4yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fx2z_vrr;
    Double I_KINETIC_H3y2z_Fx2z_vrr = PBX*I_KINETIC_H3y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fx2z_vrr;
    Double I_KINETIC_H2y3z_Fx2z_vrr = PBX*I_KINETIC_H2y3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fx2z_vrr;
    Double I_KINETIC_Hy4z_Fx2z_vrr = PBX*I_KINETIC_Hy4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fx2z_vrr;
    Double I_KINETIC_H5z_Fx2z_vrr = PBX*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fx2z_vrr;
    Double I_KINETIC_H5x_F3y_vrr = PBY*I_KINETIC_H5x_D2y_vrr+2*oned2z*I_KINETIC_H5x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_KINETIC_H4xy_F3y_vrr = PBY*I_KINETIC_H4xy_D2y_vrr+oned2z*I_KINETIC_G4x_D2y_vrr+2*oned2z*I_KINETIC_H4xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_KINETIC_H4xz_F3y_vrr = PBY*I_KINETIC_H4xz_D2y_vrr+2*oned2z*I_KINETIC_H4xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_KINETIC_H3x2y_F3y_vrr = PBY*I_KINETIC_H3x2y_D2y_vrr+2*oned2z*I_KINETIC_G3xy_D2y_vrr+2*oned2z*I_KINETIC_H3x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_KINETIC_H3x2z_F3y_vrr = PBY*I_KINETIC_H3x2z_D2y_vrr+2*oned2z*I_KINETIC_H3x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_KINETIC_H2x3y_F3y_vrr = PBY*I_KINETIC_H2x3y_D2y_vrr+3*oned2z*I_KINETIC_G2x2y_D2y_vrr+2*oned2z*I_KINETIC_H2x3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_KINETIC_H2x2yz_F3y_vrr = PBY*I_KINETIC_H2x2yz_D2y_vrr+2*oned2z*I_KINETIC_G2xyz_D2y_vrr+2*oned2z*I_KINETIC_H2x2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H2x2yz_Py_vrr;
    Double I_KINETIC_H2x3z_F3y_vrr = PBY*I_KINETIC_H2x3z_D2y_vrr+2*oned2z*I_KINETIC_H2x3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H2x3z_Py_vrr;
    Double I_KINETIC_Hx4y_F3y_vrr = PBY*I_KINETIC_Hx4y_D2y_vrr+4*oned2z*I_KINETIC_Gx3y_D2y_vrr+2*oned2z*I_KINETIC_Hx4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_KINETIC_Hx4z_F3y_vrr = PBY*I_KINETIC_Hx4z_D2y_vrr+2*oned2z*I_KINETIC_Hx4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_KINETIC_H5y_F3y_vrr = PBY*I_KINETIC_H5y_D2y_vrr+5*oned2z*I_KINETIC_G4y_D2y_vrr+2*oned2z*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_KINETIC_H4yz_F3y_vrr = PBY*I_KINETIC_H4yz_D2y_vrr+4*oned2z*I_KINETIC_G3yz_D2y_vrr+2*oned2z*I_KINETIC_H4yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_KINETIC_H3y2z_F3y_vrr = PBY*I_KINETIC_H3y2z_D2y_vrr+3*oned2z*I_KINETIC_G2y2z_D2y_vrr+2*oned2z*I_KINETIC_H3y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_KINETIC_H2y3z_F3y_vrr = PBY*I_KINETIC_H2y3z_D2y_vrr+2*oned2z*I_KINETIC_Gy3z_D2y_vrr+2*oned2z*I_KINETIC_H2y3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H2y3z_Py_vrr;
    Double I_KINETIC_Hy4z_F3y_vrr = PBY*I_KINETIC_Hy4z_D2y_vrr+oned2z*I_KINETIC_G4z_D2y_vrr+2*oned2z*I_KINETIC_Hy4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_KINETIC_H5z_F3y_vrr = PBY*I_KINETIC_H5z_D2y_vrr+2*oned2z*I_KINETIC_H5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_KINETIC_H5x_F2yz_vrr = PBZ*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F2yz_vrr;
    Double I_KINETIC_H4xy_F2yz_vrr = PBZ*I_KINETIC_H4xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F2yz_vrr;
    Double I_KINETIC_H4xz_F2yz_vrr = PBZ*I_KINETIC_H4xz_D2y_vrr+oned2z*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F2yz_vrr;
    Double I_KINETIC_H3x2y_F2yz_vrr = PBZ*I_KINETIC_H3x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F2yz_vrr;
    Double I_KINETIC_H3x2z_F2yz_vrr = PBZ*I_KINETIC_H3x2z_D2y_vrr+2*oned2z*I_KINETIC_G3xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F2yz_vrr;
    Double I_KINETIC_H2x3y_F2yz_vrr = PBZ*I_KINETIC_H2x3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F2yz_vrr;
    Double I_KINETIC_H2x2yz_F2yz_vrr = PBZ*I_KINETIC_H2x2yz_D2y_vrr+oned2z*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F2yz_vrr;
    Double I_KINETIC_H2x3z_F2yz_vrr = PBZ*I_KINETIC_H2x3z_D2y_vrr+3*oned2z*I_KINETIC_G2x2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F2yz_vrr;
    Double I_KINETIC_Hx4y_F2yz_vrr = PBZ*I_KINETIC_Hx4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F2yz_vrr;
    Double I_KINETIC_Hx4z_F2yz_vrr = PBZ*I_KINETIC_Hx4z_D2y_vrr+4*oned2z*I_KINETIC_Gx3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F2yz_vrr;
    Double I_KINETIC_H5y_F2yz_vrr = PBZ*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F2yz_vrr;
    Double I_KINETIC_H4yz_F2yz_vrr = PBZ*I_KINETIC_H4yz_D2y_vrr+oned2z*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F2yz_vrr;
    Double I_KINETIC_H3y2z_F2yz_vrr = PBZ*I_KINETIC_H3y2z_D2y_vrr+2*oned2z*I_KINETIC_G3yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F2yz_vrr;
    Double I_KINETIC_H2y3z_F2yz_vrr = PBZ*I_KINETIC_H2y3z_D2y_vrr+3*oned2z*I_KINETIC_G2y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F2yz_vrr;
    Double I_KINETIC_Hy4z_F2yz_vrr = PBZ*I_KINETIC_Hy4z_D2y_vrr+4*oned2z*I_KINETIC_Gy3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F2yz_vrr;
    Double I_KINETIC_H5z_F2yz_vrr = PBZ*I_KINETIC_H5z_D2y_vrr+5*oned2z*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F2yz_vrr;
    Double I_KINETIC_H5x_Fy2z_vrr = PBY*I_KINETIC_H5x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fy2z_vrr;
    Double I_KINETIC_H4xy_Fy2z_vrr = PBY*I_KINETIC_H4xy_D2z_vrr+oned2z*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fy2z_vrr;
    Double I_KINETIC_H4xz_Fy2z_vrr = PBY*I_KINETIC_H4xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fy2z_vrr;
    Double I_KINETIC_H3x2y_Fy2z_vrr = PBY*I_KINETIC_H3x2y_D2z_vrr+2*oned2z*I_KINETIC_G3xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fy2z_vrr;
    Double I_KINETIC_H3x2z_Fy2z_vrr = PBY*I_KINETIC_H3x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fy2z_vrr;
    Double I_KINETIC_H2x3y_Fy2z_vrr = PBY*I_KINETIC_H2x3y_D2z_vrr+3*oned2z*I_KINETIC_G2x2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fy2z_vrr;
    Double I_KINETIC_H2x2yz_Fy2z_vrr = PBY*I_KINETIC_H2x2yz_D2z_vrr+2*oned2z*I_KINETIC_G2xyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fy2z_vrr;
    Double I_KINETIC_H2x3z_Fy2z_vrr = PBY*I_KINETIC_H2x3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fy2z_vrr;
    Double I_KINETIC_Hx4y_Fy2z_vrr = PBY*I_KINETIC_Hx4y_D2z_vrr+4*oned2z*I_KINETIC_Gx3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fy2z_vrr;
    Double I_KINETIC_Hx4z_Fy2z_vrr = PBY*I_KINETIC_Hx4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fy2z_vrr;
    Double I_KINETIC_H5y_Fy2z_vrr = PBY*I_KINETIC_H5y_D2z_vrr+5*oned2z*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fy2z_vrr;
    Double I_KINETIC_H4yz_Fy2z_vrr = PBY*I_KINETIC_H4yz_D2z_vrr+4*oned2z*I_KINETIC_G3yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fy2z_vrr;
    Double I_KINETIC_H3y2z_Fy2z_vrr = PBY*I_KINETIC_H3y2z_D2z_vrr+3*oned2z*I_KINETIC_G2y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fy2z_vrr;
    Double I_KINETIC_H2y3z_Fy2z_vrr = PBY*I_KINETIC_H2y3z_D2z_vrr+2*oned2z*I_KINETIC_Gy3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fy2z_vrr;
    Double I_KINETIC_Hy4z_Fy2z_vrr = PBY*I_KINETIC_Hy4z_D2z_vrr+oned2z*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fy2z_vrr;
    Double I_KINETIC_H5z_Fy2z_vrr = PBY*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fy2z_vrr;
    Double I_KINETIC_H5x_F3z_vrr = PBZ*I_KINETIC_H5x_D2z_vrr+2*oned2z*I_KINETIC_H5x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_KINETIC_H4xy_F3z_vrr = PBZ*I_KINETIC_H4xy_D2z_vrr+2*oned2z*I_KINETIC_H4xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_KINETIC_H4xz_F3z_vrr = PBZ*I_KINETIC_H4xz_D2z_vrr+oned2z*I_KINETIC_G4x_D2z_vrr+2*oned2z*I_KINETIC_H4xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H4xz_Pz_vrr;
    Double I_KINETIC_H3x2y_F3z_vrr = PBZ*I_KINETIC_H3x2y_D2z_vrr+2*oned2z*I_KINETIC_H3x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H3x2y_Pz_vrr;
    Double I_KINETIC_H3x2z_F3z_vrr = PBZ*I_KINETIC_H3x2z_D2z_vrr+2*oned2z*I_KINETIC_G3xz_D2z_vrr+2*oned2z*I_KINETIC_H3x2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H3x2z_Pz_vrr;
    Double I_KINETIC_H2x3y_F3z_vrr = PBZ*I_KINETIC_H2x3y_D2z_vrr+2*oned2z*I_KINETIC_H2x3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H2x3y_Pz_vrr;
    Double I_KINETIC_H2x2yz_F3z_vrr = PBZ*I_KINETIC_H2x2yz_D2z_vrr+oned2z*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_H2x2yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr;
    Double I_KINETIC_H2x3z_F3z_vrr = PBZ*I_KINETIC_H2x3z_D2z_vrr+3*oned2z*I_KINETIC_G2x2z_D2z_vrr+2*oned2z*I_KINETIC_H2x3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H2x3z_Pz_vrr;
    Double I_KINETIC_Hx4y_F3z_vrr = PBZ*I_KINETIC_Hx4y_D2z_vrr+2*oned2z*I_KINETIC_Hx4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_KINETIC_Hx4z_F3z_vrr = PBZ*I_KINETIC_Hx4z_D2z_vrr+4*oned2z*I_KINETIC_Gx3z_D2z_vrr+2*oned2z*I_KINETIC_Hx4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_KINETIC_H5y_F3z_vrr = PBZ*I_KINETIC_H5y_D2z_vrr+2*oned2z*I_KINETIC_H5y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_KINETIC_H4yz_F3z_vrr = PBZ*I_KINETIC_H4yz_D2z_vrr+oned2z*I_KINETIC_G4y_D2z_vrr+2*oned2z*I_KINETIC_H4yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H4yz_Pz_vrr;
    Double I_KINETIC_H3y2z_F3z_vrr = PBZ*I_KINETIC_H3y2z_D2z_vrr+2*oned2z*I_KINETIC_G3yz_D2z_vrr+2*oned2z*I_KINETIC_H3y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H3y2z_Pz_vrr;
    Double I_KINETIC_H2y3z_F3z_vrr = PBZ*I_KINETIC_H2y3z_D2z_vrr+3*oned2z*I_KINETIC_G2y2z_D2z_vrr+2*oned2z*I_KINETIC_H2y3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H2y3z_Pz_vrr;
    Double I_KINETIC_Hy4z_F3z_vrr = PBZ*I_KINETIC_Hy4z_D2z_vrr+4*oned2z*I_KINETIC_Gy3z_D2z_vrr+2*oned2z*I_KINETIC_Hy4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Hy4z_Pz_vrr;
    Double I_KINETIC_H5z_F3z_vrr = PBZ*I_KINETIC_H5z_D2z_vrr+5*oned2z*I_KINETIC_G4z_D2z_vrr+2*oned2z*I_KINETIC_H5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H5z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_F
     * RHS shell quartet name: SQ_KINETIC_G_F
     * RHS shell quartet name: SQ_KINETIC_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
     ************************************************************/
    Double I_KINETIC_I6x_F3x_vrr = PAX*I_KINETIC_H5x_F3x_vrr+5*oned2z*I_KINETIC_G4x_F3x_vrr+3*oned2z*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I6x_F3x_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_KINETIC_I5xy_F3x_vrr = PAY*I_KINETIC_H5x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_F3x_vrr;
    Double I_KINETIC_I5xz_F3x_vrr = PAZ*I_KINETIC_H5x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_F3x_vrr;
    Double I_KINETIC_I4x2y_F3x_vrr = PAY*I_KINETIC_H4xy_F3x_vrr+oned2z*I_KINETIC_G4x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_F3x_vrr-bdz*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_KINETIC_I4xyz_F3x_vrr = PAZ*I_KINETIC_H4xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_F3x_vrr;
    Double I_KINETIC_I4x2z_F3x_vrr = PAZ*I_KINETIC_H4xz_F3x_vrr+oned2z*I_KINETIC_G4x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_KINETIC_I3x3y_F3x_vrr = PAY*I_KINETIC_H3x2y_F3x_vrr+2*oned2z*I_KINETIC_G3xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_F3x_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_F3x_vrr;
    Double I_KINETIC_I3x2yz_F3x_vrr = PAZ*I_KINETIC_H3x2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_F3x_vrr;
    Double I_KINETIC_I3xy2z_F3x_vrr = PAY*I_KINETIC_H3x2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_F3x_vrr;
    Double I_KINETIC_I3x3z_F3x_vrr = PAZ*I_KINETIC_H3x2z_F3x_vrr+2*oned2z*I_KINETIC_G3xz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_F3x_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_F3x_vrr;
    Double I_KINETIC_I2x4y_F3x_vrr = PAX*I_KINETIC_Hx4y_F3x_vrr+oned2z*I_KINETIC_G4y_F3x_vrr+3*oned2z*I_KINETIC_Hx4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_F3x_vrr-bdz*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_KINETIC_I2x3yz_F3x_vrr = PAZ*I_KINETIC_H2x3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_F3x_vrr;
    Double I_KINETIC_I2x2y2z_F3x_vrr = PAZ*I_KINETIC_H2x2yz_F3x_vrr+oned2z*I_KINETIC_G2x2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_F3x_vrr;
    Double I_KINETIC_I2xy3z_F3x_vrr = PAY*I_KINETIC_H2x3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_F3x_vrr;
    Double I_KINETIC_I2x4z_F3x_vrr = PAX*I_KINETIC_Hx4z_F3x_vrr+oned2z*I_KINETIC_G4z_F3x_vrr+3*oned2z*I_KINETIC_Hx4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_F3x_vrr-bdz*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_KINETIC_Ix5y_F3x_vrr = PAX*I_KINETIC_H5y_F3x_vrr+3*oned2z*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_F3x_vrr;
    Double I_KINETIC_Ix4yz_F3x_vrr = PAZ*I_KINETIC_Hx4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_F3x_vrr;
    Double I_KINETIC_Ix3y2z_F3x_vrr = PAX*I_KINETIC_H3y2z_F3x_vrr+3*oned2z*I_KINETIC_H3y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_F3x_vrr;
    Double I_KINETIC_Ix2y3z_F3x_vrr = PAX*I_KINETIC_H2y3z_F3x_vrr+3*oned2z*I_KINETIC_H2y3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_F3x_vrr;
    Double I_KINETIC_Ixy4z_F3x_vrr = PAY*I_KINETIC_Hx4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_F3x_vrr;
    Double I_KINETIC_Ix5z_F3x_vrr = PAX*I_KINETIC_H5z_F3x_vrr+3*oned2z*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_F3x_vrr;
    Double I_KINETIC_I6y_F3x_vrr = PAY*I_KINETIC_H5y_F3x_vrr+5*oned2z*I_KINETIC_G4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I6y_F3x_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_KINETIC_I5yz_F3x_vrr = PAZ*I_KINETIC_H5y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_F3x_vrr;
    Double I_KINETIC_I4y2z_F3x_vrr = PAZ*I_KINETIC_H4yz_F3x_vrr+oned2z*I_KINETIC_G4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_KINETIC_I3y3z_F3x_vrr = PAZ*I_KINETIC_H3y2z_F3x_vrr+2*oned2z*I_KINETIC_G3yz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_F3x_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_F3x_vrr;
    Double I_KINETIC_I2y4z_F3x_vrr = PAY*I_KINETIC_Hy4z_F3x_vrr+oned2z*I_KINETIC_G4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_F3x_vrr-bdz*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_KINETIC_Iy5z_F3x_vrr = PAY*I_KINETIC_H5z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_F3x_vrr;
    Double I_KINETIC_I6z_F3x_vrr = PAZ*I_KINETIC_H5z_F3x_vrr+5*oned2z*I_KINETIC_G4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_I6z_F3x_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_KINETIC_I6x_F2xy_vrr = PAX*I_KINETIC_H5x_F2xy_vrr+5*oned2z*I_KINETIC_G4x_F2xy_vrr+2*oned2z*I_KINETIC_H5x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I6x_F2xy_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_KINETIC_I5xy_F2xy_vrr = PAY*I_KINETIC_H5x_F2xy_vrr+oned2z*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_F2xy_vrr;
    Double I_KINETIC_I5xz_F2xy_vrr = PAZ*I_KINETIC_H5x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_F2xy_vrr;
    Double I_KINETIC_I4x2y_F2xy_vrr = PAY*I_KINETIC_H4xy_F2xy_vrr+oned2z*I_KINETIC_G4x_F2xy_vrr+oned2z*I_KINETIC_H4xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_F2xy_vrr-bdz*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_KINETIC_I4xyz_F2xy_vrr = PAZ*I_KINETIC_H4xy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_F2xy_vrr;
    Double I_KINETIC_I4x2z_F2xy_vrr = PAZ*I_KINETIC_H4xz_F2xy_vrr+oned2z*I_KINETIC_G4x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_KINETIC_I3x3y_F2xy_vrr = PAY*I_KINETIC_H3x2y_F2xy_vrr+2*oned2z*I_KINETIC_G3xy_F2xy_vrr+oned2z*I_KINETIC_H3x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_F2xy_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_F2xy_vrr;
    Double I_KINETIC_I3x2yz_F2xy_vrr = PAZ*I_KINETIC_H3x2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_F2xy_vrr;
    Double I_KINETIC_I3xy2z_F2xy_vrr = PAY*I_KINETIC_H3x2z_F2xy_vrr+oned2z*I_KINETIC_H3x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_F2xy_vrr;
    Double I_KINETIC_I3x3z_F2xy_vrr = PAZ*I_KINETIC_H3x2z_F2xy_vrr+2*oned2z*I_KINETIC_G3xz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_F2xy_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_F2xy_vrr;
    Double I_KINETIC_I2x4y_F2xy_vrr = PAX*I_KINETIC_Hx4y_F2xy_vrr+oned2z*I_KINETIC_G4y_F2xy_vrr+2*oned2z*I_KINETIC_Hx4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_F2xy_vrr-bdz*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_KINETIC_I2x3yz_F2xy_vrr = PAZ*I_KINETIC_H2x3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_F2xy_vrr;
    Double I_KINETIC_I2x2y2z_F2xy_vrr = PAZ*I_KINETIC_H2x2yz_F2xy_vrr+oned2z*I_KINETIC_G2x2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr;
    Double I_KINETIC_I2xy3z_F2xy_vrr = PAY*I_KINETIC_H2x3z_F2xy_vrr+oned2z*I_KINETIC_H2x3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_F2xy_vrr;
    Double I_KINETIC_I2x4z_F2xy_vrr = PAX*I_KINETIC_Hx4z_F2xy_vrr+oned2z*I_KINETIC_G4z_F2xy_vrr+2*oned2z*I_KINETIC_Hx4z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_KINETIC_Ix5y_F2xy_vrr = PAX*I_KINETIC_H5y_F2xy_vrr+2*oned2z*I_KINETIC_H5y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_F2xy_vrr;
    Double I_KINETIC_Ix4yz_F2xy_vrr = PAZ*I_KINETIC_Hx4y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_F2xy_vrr;
    Double I_KINETIC_Ix3y2z_F2xy_vrr = PAX*I_KINETIC_H3y2z_F2xy_vrr+2*oned2z*I_KINETIC_H3y2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_F2xy_vrr;
    Double I_KINETIC_Ix2y3z_F2xy_vrr = PAX*I_KINETIC_H2y3z_F2xy_vrr+2*oned2z*I_KINETIC_H2y3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_F2xy_vrr;
    Double I_KINETIC_Ixy4z_F2xy_vrr = PAY*I_KINETIC_Hx4z_F2xy_vrr+oned2z*I_KINETIC_Hx4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_F2xy_vrr;
    Double I_KINETIC_Ix5z_F2xy_vrr = PAX*I_KINETIC_H5z_F2xy_vrr+2*oned2z*I_KINETIC_H5z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_F2xy_vrr;
    Double I_KINETIC_I6y_F2xy_vrr = PAY*I_KINETIC_H5y_F2xy_vrr+5*oned2z*I_KINETIC_G4y_F2xy_vrr+oned2z*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I6y_F2xy_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_KINETIC_I5yz_F2xy_vrr = PAZ*I_KINETIC_H5y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_F2xy_vrr;
    Double I_KINETIC_I4y2z_F2xy_vrr = PAZ*I_KINETIC_H4yz_F2xy_vrr+oned2z*I_KINETIC_G4y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_KINETIC_I3y3z_F2xy_vrr = PAZ*I_KINETIC_H3y2z_F2xy_vrr+2*oned2z*I_KINETIC_G3yz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_F2xy_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_F2xy_vrr;
    Double I_KINETIC_I2y4z_F2xy_vrr = PAY*I_KINETIC_Hy4z_F2xy_vrr+oned2z*I_KINETIC_G4z_F2xy_vrr+oned2z*I_KINETIC_Hy4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_KINETIC_Iy5z_F2xy_vrr = PAY*I_KINETIC_H5z_F2xy_vrr+oned2z*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_F2xy_vrr;
    Double I_KINETIC_I6z_F2xy_vrr = PAZ*I_KINETIC_H5z_F2xy_vrr+5*oned2z*I_KINETIC_G4z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_I6z_F2xy_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_KINETIC_I6x_F2xz_vrr = PAX*I_KINETIC_H5x_F2xz_vrr+5*oned2z*I_KINETIC_G4x_F2xz_vrr+2*oned2z*I_KINETIC_H5x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I6x_F2xz_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_KINETIC_I5xy_F2xz_vrr = PAY*I_KINETIC_H5x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_F2xz_vrr;
    Double I_KINETIC_I5xz_F2xz_vrr = PAZ*I_KINETIC_H5x_F2xz_vrr+oned2z*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_F2xz_vrr;
    Double I_KINETIC_I4x2y_F2xz_vrr = PAY*I_KINETIC_H4xy_F2xz_vrr+oned2z*I_KINETIC_G4x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_F2xz_vrr-bdz*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_KINETIC_I4xyz_F2xz_vrr = PAZ*I_KINETIC_H4xy_F2xz_vrr+oned2z*I_KINETIC_H4xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_F2xz_vrr;
    Double I_KINETIC_I4x2z_F2xz_vrr = PAZ*I_KINETIC_H4xz_F2xz_vrr+oned2z*I_KINETIC_G4x_F2xz_vrr+oned2z*I_KINETIC_H4xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_KINETIC_I3x3y_F2xz_vrr = PAY*I_KINETIC_H3x2y_F2xz_vrr+2*oned2z*I_KINETIC_G3xy_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_F2xz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_F2xz_vrr;
    Double I_KINETIC_I3x2yz_F2xz_vrr = PAZ*I_KINETIC_H3x2y_F2xz_vrr+oned2z*I_KINETIC_H3x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_F2xz_vrr;
    Double I_KINETIC_I3xy2z_F2xz_vrr = PAY*I_KINETIC_H3x2z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_F2xz_vrr;
    Double I_KINETIC_I3x3z_F2xz_vrr = PAZ*I_KINETIC_H3x2z_F2xz_vrr+2*oned2z*I_KINETIC_G3xz_F2xz_vrr+oned2z*I_KINETIC_H3x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_F2xz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_F2xz_vrr;
    Double I_KINETIC_I2x4y_F2xz_vrr = PAX*I_KINETIC_Hx4y_F2xz_vrr+oned2z*I_KINETIC_G4y_F2xz_vrr+2*oned2z*I_KINETIC_Hx4y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_F2xz_vrr-bdz*I_TWOBODYOVERLAP_G4y_F2xz_vrr;
    Double I_KINETIC_I2x3yz_F2xz_vrr = PAZ*I_KINETIC_H2x3y_F2xz_vrr+oned2z*I_KINETIC_H2x3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_F2xz_vrr;
    Double I_KINETIC_I2x2y2z_F2xz_vrr = PAZ*I_KINETIC_H2x2yz_F2xz_vrr+oned2z*I_KINETIC_G2x2y_F2xz_vrr+oned2z*I_KINETIC_H2x2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_F2xz_vrr;
    Double I_KINETIC_I2xy3z_F2xz_vrr = PAY*I_KINETIC_H2x3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_F2xz_vrr;
    Double I_KINETIC_I2x4z_F2xz_vrr = PAX*I_KINETIC_Hx4z_F2xz_vrr+oned2z*I_KINETIC_G4z_F2xz_vrr+2*oned2z*I_KINETIC_Hx4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_KINETIC_Ix5y_F2xz_vrr = PAX*I_KINETIC_H5y_F2xz_vrr+2*oned2z*I_KINETIC_H5y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_F2xz_vrr;
    Double I_KINETIC_Ix4yz_F2xz_vrr = PAZ*I_KINETIC_Hx4y_F2xz_vrr+oned2z*I_KINETIC_Hx4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_F2xz_vrr;
    Double I_KINETIC_Ix3y2z_F2xz_vrr = PAX*I_KINETIC_H3y2z_F2xz_vrr+2*oned2z*I_KINETIC_H3y2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_F2xz_vrr;
    Double I_KINETIC_Ix2y3z_F2xz_vrr = PAX*I_KINETIC_H2y3z_F2xz_vrr+2*oned2z*I_KINETIC_H2y3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_F2xz_vrr;
    Double I_KINETIC_Ixy4z_F2xz_vrr = PAY*I_KINETIC_Hx4z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_F2xz_vrr;
    Double I_KINETIC_Ix5z_F2xz_vrr = PAX*I_KINETIC_H5z_F2xz_vrr+2*oned2z*I_KINETIC_H5z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_F2xz_vrr;
    Double I_KINETIC_I6y_F2xz_vrr = PAY*I_KINETIC_H5y_F2xz_vrr+5*oned2z*I_KINETIC_G4y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_I6y_F2xz_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_F2xz_vrr;
    Double I_KINETIC_I5yz_F2xz_vrr = PAZ*I_KINETIC_H5y_F2xz_vrr+oned2z*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_F2xz_vrr;
    Double I_KINETIC_I4y2z_F2xz_vrr = PAZ*I_KINETIC_H4yz_F2xz_vrr+oned2z*I_KINETIC_G4y_F2xz_vrr+oned2z*I_KINETIC_H4yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_G4y_F2xz_vrr;
    Double I_KINETIC_I3y3z_F2xz_vrr = PAZ*I_KINETIC_H3y2z_F2xz_vrr+2*oned2z*I_KINETIC_G3yz_F2xz_vrr+oned2z*I_KINETIC_H3y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_F2xz_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_F2xz_vrr;
    Double I_KINETIC_I2y4z_F2xz_vrr = PAY*I_KINETIC_Hy4z_F2xz_vrr+oned2z*I_KINETIC_G4z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_KINETIC_Iy5z_F2xz_vrr = PAY*I_KINETIC_H5z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_F2xz_vrr;
    Double I_KINETIC_I6z_F2xz_vrr = PAZ*I_KINETIC_H5z_F2xz_vrr+5*oned2z*I_KINETIC_G4z_F2xz_vrr+oned2z*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I6z_F2xz_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_KINETIC_I6x_Fx2y_vrr = PAX*I_KINETIC_H5x_Fx2y_vrr+5*oned2z*I_KINETIC_G4x_Fx2y_vrr+oned2z*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Fx2y_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_KINETIC_I5xy_Fx2y_vrr = PAY*I_KINETIC_H5x_Fx2y_vrr+2*oned2z*I_KINETIC_H5x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Fx2y_vrr;
    Double I_KINETIC_I5xz_Fx2y_vrr = PAZ*I_KINETIC_H5x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Fx2y_vrr;
    Double I_KINETIC_I4x2y_Fx2y_vrr = PAY*I_KINETIC_H4xy_Fx2y_vrr+oned2z*I_KINETIC_G4x_Fx2y_vrr+2*oned2z*I_KINETIC_H4xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_KINETIC_I4xyz_Fx2y_vrr = PAZ*I_KINETIC_H4xy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Fx2y_vrr;
    Double I_KINETIC_I4x2z_Fx2y_vrr = PAZ*I_KINETIC_H4xz_Fx2y_vrr+oned2z*I_KINETIC_G4x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_KINETIC_I3x3y_Fx2y_vrr = PAY*I_KINETIC_H3x2y_Fx2y_vrr+2*oned2z*I_KINETIC_G3xy_Fx2y_vrr+2*oned2z*I_KINETIC_H3x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Fx2y_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr;
    Double I_KINETIC_I3x2yz_Fx2y_vrr = PAZ*I_KINETIC_H3x2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Fx2y_vrr;
    Double I_KINETIC_I3xy2z_Fx2y_vrr = PAY*I_KINETIC_H3x2z_Fx2y_vrr+2*oned2z*I_KINETIC_H3x2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Fx2y_vrr;
    Double I_KINETIC_I3x3z_Fx2y_vrr = PAZ*I_KINETIC_H3x2z_Fx2y_vrr+2*oned2z*I_KINETIC_G3xz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Fx2y_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Fx2y_vrr;
    Double I_KINETIC_I2x4y_Fx2y_vrr = PAX*I_KINETIC_Hx4y_Fx2y_vrr+oned2z*I_KINETIC_G4y_Fx2y_vrr+oned2z*I_KINETIC_Hx4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_KINETIC_I2x3yz_Fx2y_vrr = PAZ*I_KINETIC_H2x3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Fx2y_vrr;
    Double I_KINETIC_I2x2y2z_Fx2y_vrr = PAZ*I_KINETIC_H2x2yz_Fx2y_vrr+oned2z*I_KINETIC_G2x2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr;
    Double I_KINETIC_I2xy3z_Fx2y_vrr = PAY*I_KINETIC_H2x3z_Fx2y_vrr+2*oned2z*I_KINETIC_H2x3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Fx2y_vrr;
    Double I_KINETIC_I2x4z_Fx2y_vrr = PAX*I_KINETIC_Hx4z_Fx2y_vrr+oned2z*I_KINETIC_G4z_Fx2y_vrr+oned2z*I_KINETIC_Hx4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_KINETIC_Ix5y_Fx2y_vrr = PAX*I_KINETIC_H5y_Fx2y_vrr+oned2z*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Fx2y_vrr;
    Double I_KINETIC_Ix4yz_Fx2y_vrr = PAZ*I_KINETIC_Hx4y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Fx2y_vrr;
    Double I_KINETIC_Ix3y2z_Fx2y_vrr = PAX*I_KINETIC_H3y2z_Fx2y_vrr+oned2z*I_KINETIC_H3y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Fx2y_vrr;
    Double I_KINETIC_Ix2y3z_Fx2y_vrr = PAX*I_KINETIC_H2y3z_Fx2y_vrr+oned2z*I_KINETIC_H2y3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Fx2y_vrr;
    Double I_KINETIC_Ixy4z_Fx2y_vrr = PAY*I_KINETIC_Hx4z_Fx2y_vrr+2*oned2z*I_KINETIC_Hx4z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Fx2y_vrr;
    Double I_KINETIC_Ix5z_Fx2y_vrr = PAX*I_KINETIC_H5z_Fx2y_vrr+oned2z*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Fx2y_vrr;
    Double I_KINETIC_I6y_Fx2y_vrr = PAY*I_KINETIC_H5y_Fx2y_vrr+5*oned2z*I_KINETIC_G4y_Fx2y_vrr+2*oned2z*I_KINETIC_H5y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Fx2y_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_KINETIC_I5yz_Fx2y_vrr = PAZ*I_KINETIC_H5y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Fx2y_vrr;
    Double I_KINETIC_I4y2z_Fx2y_vrr = PAZ*I_KINETIC_H4yz_Fx2y_vrr+oned2z*I_KINETIC_G4y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_KINETIC_I3y3z_Fx2y_vrr = PAZ*I_KINETIC_H3y2z_Fx2y_vrr+2*oned2z*I_KINETIC_G3yz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Fx2y_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Fx2y_vrr;
    Double I_KINETIC_I2y4z_Fx2y_vrr = PAY*I_KINETIC_Hy4z_Fx2y_vrr+oned2z*I_KINETIC_G4z_Fx2y_vrr+2*oned2z*I_KINETIC_Hy4z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_KINETIC_Iy5z_Fx2y_vrr = PAY*I_KINETIC_H5z_Fx2y_vrr+2*oned2z*I_KINETIC_H5z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Fx2y_vrr;
    Double I_KINETIC_I6z_Fx2y_vrr = PAZ*I_KINETIC_H5z_Fx2y_vrr+5*oned2z*I_KINETIC_G4z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Fx2y_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_KINETIC_I6x_Fxyz_vrr = PAX*I_KINETIC_H5x_Fxyz_vrr+5*oned2z*I_KINETIC_G4x_Fxyz_vrr+oned2z*I_KINETIC_H5x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Fxyz_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Fxyz_vrr;
    Double I_KINETIC_I5xy_Fxyz_vrr = PAY*I_KINETIC_H5x_Fxyz_vrr+oned2z*I_KINETIC_H5x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Fxyz_vrr;
    Double I_KINETIC_I5xz_Fxyz_vrr = PAZ*I_KINETIC_H5x_Fxyz_vrr+oned2z*I_KINETIC_H5x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Fxyz_vrr;
    Double I_KINETIC_I4x2y_Fxyz_vrr = PAY*I_KINETIC_H4xy_Fxyz_vrr+oned2z*I_KINETIC_G4x_Fxyz_vrr+oned2z*I_KINETIC_H4xy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_G4x_Fxyz_vrr;
    Double I_KINETIC_I4xyz_Fxyz_vrr = PAZ*I_KINETIC_H4xy_Fxyz_vrr+oned2z*I_KINETIC_H4xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Fxyz_vrr;
    Double I_KINETIC_I4x2z_Fxyz_vrr = PAZ*I_KINETIC_H4xz_Fxyz_vrr+oned2z*I_KINETIC_G4x_Fxyz_vrr+oned2z*I_KINETIC_H4xz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_G4x_Fxyz_vrr;
    Double I_KINETIC_I3x3y_Fxyz_vrr = PAY*I_KINETIC_H3x2y_Fxyz_vrr+2*oned2z*I_KINETIC_G3xy_Fxyz_vrr+oned2z*I_KINETIC_H3x2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr;
    Double I_KINETIC_I3x2yz_Fxyz_vrr = PAZ*I_KINETIC_H3x2y_Fxyz_vrr+oned2z*I_KINETIC_H3x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Fxyz_vrr;
    Double I_KINETIC_I3xy2z_Fxyz_vrr = PAY*I_KINETIC_H3x2z_Fxyz_vrr+oned2z*I_KINETIC_H3x2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Fxyz_vrr;
    Double I_KINETIC_I3x3z_Fxyz_vrr = PAZ*I_KINETIC_H3x2z_Fxyz_vrr+2*oned2z*I_KINETIC_G3xz_Fxyz_vrr+oned2z*I_KINETIC_H3x2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Fxyz_vrr;
    Double I_KINETIC_I2x4y_Fxyz_vrr = PAX*I_KINETIC_Hx4y_Fxyz_vrr+oned2z*I_KINETIC_G4y_Fxyz_vrr+oned2z*I_KINETIC_Hx4y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_G4y_Fxyz_vrr;
    Double I_KINETIC_I2x3yz_Fxyz_vrr = PAZ*I_KINETIC_H2x3y_Fxyz_vrr+oned2z*I_KINETIC_H2x3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Fxyz_vrr;
    Double I_KINETIC_I2x2y2z_Fxyz_vrr = PAZ*I_KINETIC_H2x2yz_Fxyz_vrr+oned2z*I_KINETIC_G2x2y_Fxyz_vrr+oned2z*I_KINETIC_H2x2yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr;
    Double I_KINETIC_I2xy3z_Fxyz_vrr = PAY*I_KINETIC_H2x3z_Fxyz_vrr+oned2z*I_KINETIC_H2x3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Fxyz_vrr;
    Double I_KINETIC_I2x4z_Fxyz_vrr = PAX*I_KINETIC_Hx4z_Fxyz_vrr+oned2z*I_KINETIC_G4z_Fxyz_vrr+oned2z*I_KINETIC_Hx4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_G4z_Fxyz_vrr;
    Double I_KINETIC_Ix5y_Fxyz_vrr = PAX*I_KINETIC_H5y_Fxyz_vrr+oned2z*I_KINETIC_H5y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Fxyz_vrr;
    Double I_KINETIC_Ix4yz_Fxyz_vrr = PAZ*I_KINETIC_Hx4y_Fxyz_vrr+oned2z*I_KINETIC_Hx4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Fxyz_vrr;
    Double I_KINETIC_Ix3y2z_Fxyz_vrr = PAX*I_KINETIC_H3y2z_Fxyz_vrr+oned2z*I_KINETIC_H3y2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Fxyz_vrr;
    Double I_KINETIC_Ix2y3z_Fxyz_vrr = PAX*I_KINETIC_H2y3z_Fxyz_vrr+oned2z*I_KINETIC_H2y3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Fxyz_vrr;
    Double I_KINETIC_Ixy4z_Fxyz_vrr = PAY*I_KINETIC_Hx4z_Fxyz_vrr+oned2z*I_KINETIC_Hx4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Fxyz_vrr;
    Double I_KINETIC_Ix5z_Fxyz_vrr = PAX*I_KINETIC_H5z_Fxyz_vrr+oned2z*I_KINETIC_H5z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Fxyz_vrr;
    Double I_KINETIC_I6y_Fxyz_vrr = PAY*I_KINETIC_H5y_Fxyz_vrr+5*oned2z*I_KINETIC_G4y_Fxyz_vrr+oned2z*I_KINETIC_H5y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Fxyz_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Fxyz_vrr;
    Double I_KINETIC_I5yz_Fxyz_vrr = PAZ*I_KINETIC_H5y_Fxyz_vrr+oned2z*I_KINETIC_H5y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Fxyz_vrr;
    Double I_KINETIC_I4y2z_Fxyz_vrr = PAZ*I_KINETIC_H4yz_Fxyz_vrr+oned2z*I_KINETIC_G4y_Fxyz_vrr+oned2z*I_KINETIC_H4yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_G4y_Fxyz_vrr;
    Double I_KINETIC_I3y3z_Fxyz_vrr = PAZ*I_KINETIC_H3y2z_Fxyz_vrr+2*oned2z*I_KINETIC_G3yz_Fxyz_vrr+oned2z*I_KINETIC_H3y2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Fxyz_vrr;
    Double I_KINETIC_I2y4z_Fxyz_vrr = PAY*I_KINETIC_Hy4z_Fxyz_vrr+oned2z*I_KINETIC_G4z_Fxyz_vrr+oned2z*I_KINETIC_Hy4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_G4z_Fxyz_vrr;
    Double I_KINETIC_Iy5z_Fxyz_vrr = PAY*I_KINETIC_H5z_Fxyz_vrr+oned2z*I_KINETIC_H5z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Fxyz_vrr;
    Double I_KINETIC_I6z_Fxyz_vrr = PAZ*I_KINETIC_H5z_Fxyz_vrr+5*oned2z*I_KINETIC_G4z_Fxyz_vrr+oned2z*I_KINETIC_H5z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Fxyz_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Fxyz_vrr;
    Double I_KINETIC_I6x_Fx2z_vrr = PAX*I_KINETIC_H5x_Fx2z_vrr+5*oned2z*I_KINETIC_G4x_Fx2z_vrr+oned2z*I_KINETIC_H5x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Fx2z_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_KINETIC_I5xy_Fx2z_vrr = PAY*I_KINETIC_H5x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Fx2z_vrr;
    Double I_KINETIC_I5xz_Fx2z_vrr = PAZ*I_KINETIC_H5x_Fx2z_vrr+2*oned2z*I_KINETIC_H5x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Fx2z_vrr;
    Double I_KINETIC_I4x2y_Fx2z_vrr = PAY*I_KINETIC_H4xy_Fx2z_vrr+oned2z*I_KINETIC_G4x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_KINETIC_I4xyz_Fx2z_vrr = PAZ*I_KINETIC_H4xy_Fx2z_vrr+2*oned2z*I_KINETIC_H4xy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Fx2z_vrr;
    Double I_KINETIC_I4x2z_Fx2z_vrr = PAZ*I_KINETIC_H4xz_Fx2z_vrr+oned2z*I_KINETIC_G4x_Fx2z_vrr+2*oned2z*I_KINETIC_H4xz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_KINETIC_I3x3y_Fx2z_vrr = PAY*I_KINETIC_H3x2y_Fx2z_vrr+2*oned2z*I_KINETIC_G3xy_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Fx2z_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr;
    Double I_KINETIC_I3x2yz_Fx2z_vrr = PAZ*I_KINETIC_H3x2y_Fx2z_vrr+2*oned2z*I_KINETIC_H3x2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Fx2z_vrr;
    Double I_KINETIC_I3xy2z_Fx2z_vrr = PAY*I_KINETIC_H3x2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Fx2z_vrr;
    Double I_KINETIC_I3x3z_Fx2z_vrr = PAZ*I_KINETIC_H3x2z_Fx2z_vrr+2*oned2z*I_KINETIC_G3xz_Fx2z_vrr+2*oned2z*I_KINETIC_H3x2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Fx2z_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Fx2z_vrr;
    Double I_KINETIC_I2x4y_Fx2z_vrr = PAX*I_KINETIC_Hx4y_Fx2z_vrr+oned2z*I_KINETIC_G4y_Fx2z_vrr+oned2z*I_KINETIC_Hx4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_KINETIC_I2x3yz_Fx2z_vrr = PAZ*I_KINETIC_H2x3y_Fx2z_vrr+2*oned2z*I_KINETIC_H2x3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Fx2z_vrr;
    Double I_KINETIC_I2x2y2z_Fx2z_vrr = PAZ*I_KINETIC_H2x2yz_Fx2z_vrr+oned2z*I_KINETIC_G2x2y_Fx2z_vrr+2*oned2z*I_KINETIC_H2x2yz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr;
    Double I_KINETIC_I2xy3z_Fx2z_vrr = PAY*I_KINETIC_H2x3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Fx2z_vrr;
    Double I_KINETIC_I2x4z_Fx2z_vrr = PAX*I_KINETIC_Hx4z_Fx2z_vrr+oned2z*I_KINETIC_G4z_Fx2z_vrr+oned2z*I_KINETIC_Hx4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_KINETIC_Ix5y_Fx2z_vrr = PAX*I_KINETIC_H5y_Fx2z_vrr+oned2z*I_KINETIC_H5y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Fx2z_vrr;
    Double I_KINETIC_Ix4yz_Fx2z_vrr = PAZ*I_KINETIC_Hx4y_Fx2z_vrr+2*oned2z*I_KINETIC_Hx4y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Fx2z_vrr;
    Double I_KINETIC_Ix3y2z_Fx2z_vrr = PAX*I_KINETIC_H3y2z_Fx2z_vrr+oned2z*I_KINETIC_H3y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Fx2z_vrr;
    Double I_KINETIC_Ix2y3z_Fx2z_vrr = PAX*I_KINETIC_H2y3z_Fx2z_vrr+oned2z*I_KINETIC_H2y3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Fx2z_vrr;
    Double I_KINETIC_Ixy4z_Fx2z_vrr = PAY*I_KINETIC_Hx4z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Fx2z_vrr;
    Double I_KINETIC_Ix5z_Fx2z_vrr = PAX*I_KINETIC_H5z_Fx2z_vrr+oned2z*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Fx2z_vrr;
    Double I_KINETIC_I6y_Fx2z_vrr = PAY*I_KINETIC_H5y_Fx2z_vrr+5*oned2z*I_KINETIC_G4y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Fx2z_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_KINETIC_I5yz_Fx2z_vrr = PAZ*I_KINETIC_H5y_Fx2z_vrr+2*oned2z*I_KINETIC_H5y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Fx2z_vrr;
    Double I_KINETIC_I4y2z_Fx2z_vrr = PAZ*I_KINETIC_H4yz_Fx2z_vrr+oned2z*I_KINETIC_G4y_Fx2z_vrr+2*oned2z*I_KINETIC_H4yz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_KINETIC_I3y3z_Fx2z_vrr = PAZ*I_KINETIC_H3y2z_Fx2z_vrr+2*oned2z*I_KINETIC_G3yz_Fx2z_vrr+2*oned2z*I_KINETIC_H3y2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Fx2z_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Fx2z_vrr;
    Double I_KINETIC_I2y4z_Fx2z_vrr = PAY*I_KINETIC_Hy4z_Fx2z_vrr+oned2z*I_KINETIC_G4z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_KINETIC_Iy5z_Fx2z_vrr = PAY*I_KINETIC_H5z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Fx2z_vrr;
    Double I_KINETIC_I6z_Fx2z_vrr = PAZ*I_KINETIC_H5z_Fx2z_vrr+5*oned2z*I_KINETIC_G4z_Fx2z_vrr+2*oned2z*I_KINETIC_H5z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Fx2z_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_KINETIC_I6x_F3y_vrr = PAX*I_KINETIC_H5x_F3y_vrr+5*oned2z*I_KINETIC_G4x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I6x_F3y_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_KINETIC_I5xy_F3y_vrr = PAY*I_KINETIC_H5x_F3y_vrr+3*oned2z*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_F3y_vrr;
    Double I_KINETIC_I5xz_F3y_vrr = PAZ*I_KINETIC_H5x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_F3y_vrr;
    Double I_KINETIC_I4x2y_F3y_vrr = PAY*I_KINETIC_H4xy_F3y_vrr+oned2z*I_KINETIC_G4x_F3y_vrr+3*oned2z*I_KINETIC_H4xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_F3y_vrr-bdz*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_KINETIC_I4xyz_F3y_vrr = PAZ*I_KINETIC_H4xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_F3y_vrr;
    Double I_KINETIC_I4x2z_F3y_vrr = PAZ*I_KINETIC_H4xz_F3y_vrr+oned2z*I_KINETIC_G4x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_KINETIC_I3x3y_F3y_vrr = PAY*I_KINETIC_H3x2y_F3y_vrr+2*oned2z*I_KINETIC_G3xy_F3y_vrr+3*oned2z*I_KINETIC_H3x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_F3y_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_F3y_vrr;
    Double I_KINETIC_I3x2yz_F3y_vrr = PAZ*I_KINETIC_H3x2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_F3y_vrr;
    Double I_KINETIC_I3xy2z_F3y_vrr = PAY*I_KINETIC_H3x2z_F3y_vrr+3*oned2z*I_KINETIC_H3x2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_F3y_vrr;
    Double I_KINETIC_I3x3z_F3y_vrr = PAZ*I_KINETIC_H3x2z_F3y_vrr+2*oned2z*I_KINETIC_G3xz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_F3y_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_F3y_vrr;
    Double I_KINETIC_I2x4y_F3y_vrr = PAX*I_KINETIC_Hx4y_F3y_vrr+oned2z*I_KINETIC_G4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_F3y_vrr-bdz*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_KINETIC_I2x3yz_F3y_vrr = PAZ*I_KINETIC_H2x3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_F3y_vrr;
    Double I_KINETIC_I2x2y2z_F3y_vrr = PAZ*I_KINETIC_H2x2yz_F3y_vrr+oned2z*I_KINETIC_G2x2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_F3y_vrr;
    Double I_KINETIC_I2xy3z_F3y_vrr = PAY*I_KINETIC_H2x3z_F3y_vrr+3*oned2z*I_KINETIC_H2x3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_F3y_vrr;
    Double I_KINETIC_I2x4z_F3y_vrr = PAX*I_KINETIC_Hx4z_F3y_vrr+oned2z*I_KINETIC_G4z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_F3y_vrr-bdz*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_KINETIC_Ix5y_F3y_vrr = PAX*I_KINETIC_H5y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_F3y_vrr;
    Double I_KINETIC_Ix4yz_F3y_vrr = PAZ*I_KINETIC_Hx4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_F3y_vrr;
    Double I_KINETIC_Ix3y2z_F3y_vrr = PAX*I_KINETIC_H3y2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_F3y_vrr;
    Double I_KINETIC_Ix2y3z_F3y_vrr = PAX*I_KINETIC_H2y3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_F3y_vrr;
    Double I_KINETIC_Ixy4z_F3y_vrr = PAY*I_KINETIC_Hx4z_F3y_vrr+3*oned2z*I_KINETIC_Hx4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_F3y_vrr;
    Double I_KINETIC_Ix5z_F3y_vrr = PAX*I_KINETIC_H5z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_F3y_vrr;
    Double I_KINETIC_I6y_F3y_vrr = PAY*I_KINETIC_H5y_F3y_vrr+5*oned2z*I_KINETIC_G4y_F3y_vrr+3*oned2z*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I6y_F3y_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_KINETIC_I5yz_F3y_vrr = PAZ*I_KINETIC_H5y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_F3y_vrr;
    Double I_KINETIC_I4y2z_F3y_vrr = PAZ*I_KINETIC_H4yz_F3y_vrr+oned2z*I_KINETIC_G4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_KINETIC_I3y3z_F3y_vrr = PAZ*I_KINETIC_H3y2z_F3y_vrr+2*oned2z*I_KINETIC_G3yz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_F3y_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_F3y_vrr;
    Double I_KINETIC_I2y4z_F3y_vrr = PAY*I_KINETIC_Hy4z_F3y_vrr+oned2z*I_KINETIC_G4z_F3y_vrr+3*oned2z*I_KINETIC_Hy4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_F3y_vrr-bdz*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_KINETIC_Iy5z_F3y_vrr = PAY*I_KINETIC_H5z_F3y_vrr+3*oned2z*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_F3y_vrr;
    Double I_KINETIC_I6z_F3y_vrr = PAZ*I_KINETIC_H5z_F3y_vrr+5*oned2z*I_KINETIC_G4z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_I6z_F3y_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_KINETIC_I6x_F2yz_vrr = PAX*I_KINETIC_H5x_F2yz_vrr+5*oned2z*I_KINETIC_G4x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_I6x_F2yz_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_F2yz_vrr;
    Double I_KINETIC_I5xy_F2yz_vrr = PAY*I_KINETIC_H5x_F2yz_vrr+2*oned2z*I_KINETIC_H5x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_F2yz_vrr;
    Double I_KINETIC_I5xz_F2yz_vrr = PAZ*I_KINETIC_H5x_F2yz_vrr+oned2z*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_F2yz_vrr;
    Double I_KINETIC_I4x2y_F2yz_vrr = PAY*I_KINETIC_H4xy_F2yz_vrr+oned2z*I_KINETIC_G4x_F2yz_vrr+2*oned2z*I_KINETIC_H4xy_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_F2yz_vrr-bdz*I_TWOBODYOVERLAP_G4x_F2yz_vrr;
    Double I_KINETIC_I4xyz_F2yz_vrr = PAZ*I_KINETIC_H4xy_F2yz_vrr+oned2z*I_KINETIC_H4xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_F2yz_vrr;
    Double I_KINETIC_I4x2z_F2yz_vrr = PAZ*I_KINETIC_H4xz_F2yz_vrr+oned2z*I_KINETIC_G4x_F2yz_vrr+oned2z*I_KINETIC_H4xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_G4x_F2yz_vrr;
    Double I_KINETIC_I3x3y_F2yz_vrr = PAY*I_KINETIC_H3x2y_F2yz_vrr+2*oned2z*I_KINETIC_G3xy_F2yz_vrr+2*oned2z*I_KINETIC_H3x2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_F2yz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_F2yz_vrr;
    Double I_KINETIC_I3x2yz_F2yz_vrr = PAZ*I_KINETIC_H3x2y_F2yz_vrr+oned2z*I_KINETIC_H3x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_F2yz_vrr;
    Double I_KINETIC_I3xy2z_F2yz_vrr = PAY*I_KINETIC_H3x2z_F2yz_vrr+2*oned2z*I_KINETIC_H3x2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_F2yz_vrr;
    Double I_KINETIC_I3x3z_F2yz_vrr = PAZ*I_KINETIC_H3x2z_F2yz_vrr+2*oned2z*I_KINETIC_G3xz_F2yz_vrr+oned2z*I_KINETIC_H3x2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_F2yz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_F2yz_vrr;
    Double I_KINETIC_I2x4y_F2yz_vrr = PAX*I_KINETIC_Hx4y_F2yz_vrr+oned2z*I_KINETIC_G4y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_F2yz_vrr-bdz*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_KINETIC_I2x3yz_F2yz_vrr = PAZ*I_KINETIC_H2x3y_F2yz_vrr+oned2z*I_KINETIC_H2x3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_F2yz_vrr;
    Double I_KINETIC_I2x2y2z_F2yz_vrr = PAZ*I_KINETIC_H2x2yz_F2yz_vrr+oned2z*I_KINETIC_G2x2y_F2yz_vrr+oned2z*I_KINETIC_H2x2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_F2yz_vrr;
    Double I_KINETIC_I2xy3z_F2yz_vrr = PAY*I_KINETIC_H2x3z_F2yz_vrr+2*oned2z*I_KINETIC_H2x3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_F2yz_vrr;
    Double I_KINETIC_I2x4z_F2yz_vrr = PAX*I_KINETIC_Hx4z_F2yz_vrr+oned2z*I_KINETIC_G4z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_KINETIC_Ix5y_F2yz_vrr = PAX*I_KINETIC_H5y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_F2yz_vrr;
    Double I_KINETIC_Ix4yz_F2yz_vrr = PAZ*I_KINETIC_Hx4y_F2yz_vrr+oned2z*I_KINETIC_Hx4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_F2yz_vrr;
    Double I_KINETIC_Ix3y2z_F2yz_vrr = PAX*I_KINETIC_H3y2z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_F2yz_vrr;
    Double I_KINETIC_Ix2y3z_F2yz_vrr = PAX*I_KINETIC_H2y3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_F2yz_vrr;
    Double I_KINETIC_Ixy4z_F2yz_vrr = PAY*I_KINETIC_Hx4z_F2yz_vrr+2*oned2z*I_KINETIC_Hx4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_F2yz_vrr;
    Double I_KINETIC_Ix5z_F2yz_vrr = PAX*I_KINETIC_H5z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_F2yz_vrr;
    Double I_KINETIC_I6y_F2yz_vrr = PAY*I_KINETIC_H5y_F2yz_vrr+5*oned2z*I_KINETIC_G4y_F2yz_vrr+2*oned2z*I_KINETIC_H5y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I6y_F2yz_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_KINETIC_I5yz_F2yz_vrr = PAZ*I_KINETIC_H5y_F2yz_vrr+oned2z*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_F2yz_vrr;
    Double I_KINETIC_I4y2z_F2yz_vrr = PAZ*I_KINETIC_H4yz_F2yz_vrr+oned2z*I_KINETIC_G4y_F2yz_vrr+oned2z*I_KINETIC_H4yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_KINETIC_I3y3z_F2yz_vrr = PAZ*I_KINETIC_H3y2z_F2yz_vrr+2*oned2z*I_KINETIC_G3yz_F2yz_vrr+oned2z*I_KINETIC_H3y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_F2yz_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_F2yz_vrr;
    Double I_KINETIC_I2y4z_F2yz_vrr = PAY*I_KINETIC_Hy4z_F2yz_vrr+oned2z*I_KINETIC_G4z_F2yz_vrr+2*oned2z*I_KINETIC_Hy4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_KINETIC_Iy5z_F2yz_vrr = PAY*I_KINETIC_H5z_F2yz_vrr+2*oned2z*I_KINETIC_H5z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_F2yz_vrr;
    Double I_KINETIC_I6z_F2yz_vrr = PAZ*I_KINETIC_H5z_F2yz_vrr+5*oned2z*I_KINETIC_G4z_F2yz_vrr+oned2z*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I6z_F2yz_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_KINETIC_I6x_Fy2z_vrr = PAX*I_KINETIC_H5x_Fy2z_vrr+5*oned2z*I_KINETIC_G4x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Fy2z_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Fy2z_vrr;
    Double I_KINETIC_I5xy_Fy2z_vrr = PAY*I_KINETIC_H5x_Fy2z_vrr+oned2z*I_KINETIC_H5x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Fy2z_vrr;
    Double I_KINETIC_I5xz_Fy2z_vrr = PAZ*I_KINETIC_H5x_Fy2z_vrr+2*oned2z*I_KINETIC_H5x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Fy2z_vrr;
    Double I_KINETIC_I4x2y_Fy2z_vrr = PAY*I_KINETIC_H4xy_Fy2z_vrr+oned2z*I_KINETIC_G4x_Fy2z_vrr+oned2z*I_KINETIC_H4xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_G4x_Fy2z_vrr;
    Double I_KINETIC_I4xyz_Fy2z_vrr = PAZ*I_KINETIC_H4xy_Fy2z_vrr+2*oned2z*I_KINETIC_H4xy_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Fy2z_vrr;
    Double I_KINETIC_I4x2z_Fy2z_vrr = PAZ*I_KINETIC_H4xz_Fy2z_vrr+oned2z*I_KINETIC_G4x_Fy2z_vrr+2*oned2z*I_KINETIC_H4xz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_G4x_Fy2z_vrr;
    Double I_KINETIC_I3x3y_Fy2z_vrr = PAY*I_KINETIC_H3x2y_Fy2z_vrr+2*oned2z*I_KINETIC_G3xy_Fy2z_vrr+oned2z*I_KINETIC_H3x2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr;
    Double I_KINETIC_I3x2yz_Fy2z_vrr = PAZ*I_KINETIC_H3x2y_Fy2z_vrr+2*oned2z*I_KINETIC_H3x2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Fy2z_vrr;
    Double I_KINETIC_I3xy2z_Fy2z_vrr = PAY*I_KINETIC_H3x2z_Fy2z_vrr+oned2z*I_KINETIC_H3x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Fy2z_vrr;
    Double I_KINETIC_I3x3z_Fy2z_vrr = PAZ*I_KINETIC_H3x2z_Fy2z_vrr+2*oned2z*I_KINETIC_G3xz_Fy2z_vrr+2*oned2z*I_KINETIC_H3x2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Fy2z_vrr;
    Double I_KINETIC_I2x4y_Fy2z_vrr = PAX*I_KINETIC_Hx4y_Fy2z_vrr+oned2z*I_KINETIC_G4y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_KINETIC_I2x3yz_Fy2z_vrr = PAZ*I_KINETIC_H2x3y_Fy2z_vrr+2*oned2z*I_KINETIC_H2x3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Fy2z_vrr;
    Double I_KINETIC_I2x2y2z_Fy2z_vrr = PAZ*I_KINETIC_H2x2yz_Fy2z_vrr+oned2z*I_KINETIC_G2x2y_Fy2z_vrr+2*oned2z*I_KINETIC_H2x2yz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr;
    Double I_KINETIC_I2xy3z_Fy2z_vrr = PAY*I_KINETIC_H2x3z_Fy2z_vrr+oned2z*I_KINETIC_H2x3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Fy2z_vrr;
    Double I_KINETIC_I2x4z_Fy2z_vrr = PAX*I_KINETIC_Hx4z_Fy2z_vrr+oned2z*I_KINETIC_G4z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_KINETIC_Ix5y_Fy2z_vrr = PAX*I_KINETIC_H5y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Fy2z_vrr;
    Double I_KINETIC_Ix4yz_Fy2z_vrr = PAZ*I_KINETIC_Hx4y_Fy2z_vrr+2*oned2z*I_KINETIC_Hx4y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Fy2z_vrr;
    Double I_KINETIC_Ix3y2z_Fy2z_vrr = PAX*I_KINETIC_H3y2z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Fy2z_vrr;
    Double I_KINETIC_Ix2y3z_Fy2z_vrr = PAX*I_KINETIC_H2y3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Fy2z_vrr;
    Double I_KINETIC_Ixy4z_Fy2z_vrr = PAY*I_KINETIC_Hx4z_Fy2z_vrr+oned2z*I_KINETIC_Hx4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Fy2z_vrr;
    Double I_KINETIC_Ix5z_Fy2z_vrr = PAX*I_KINETIC_H5z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Fy2z_vrr;
    Double I_KINETIC_I6y_Fy2z_vrr = PAY*I_KINETIC_H5y_Fy2z_vrr+5*oned2z*I_KINETIC_G4y_Fy2z_vrr+oned2z*I_KINETIC_H5y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Fy2z_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_KINETIC_I5yz_Fy2z_vrr = PAZ*I_KINETIC_H5y_Fy2z_vrr+2*oned2z*I_KINETIC_H5y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Fy2z_vrr;
    Double I_KINETIC_I4y2z_Fy2z_vrr = PAZ*I_KINETIC_H4yz_Fy2z_vrr+oned2z*I_KINETIC_G4y_Fy2z_vrr+2*oned2z*I_KINETIC_H4yz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_KINETIC_I3y3z_Fy2z_vrr = PAZ*I_KINETIC_H3y2z_Fy2z_vrr+2*oned2z*I_KINETIC_G3yz_Fy2z_vrr+2*oned2z*I_KINETIC_H3y2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Fy2z_vrr;
    Double I_KINETIC_I2y4z_Fy2z_vrr = PAY*I_KINETIC_Hy4z_Fy2z_vrr+oned2z*I_KINETIC_G4z_Fy2z_vrr+oned2z*I_KINETIC_Hy4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_KINETIC_Iy5z_Fy2z_vrr = PAY*I_KINETIC_H5z_Fy2z_vrr+oned2z*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Fy2z_vrr;
    Double I_KINETIC_I6z_Fy2z_vrr = PAZ*I_KINETIC_H5z_Fy2z_vrr+5*oned2z*I_KINETIC_G4z_Fy2z_vrr+2*oned2z*I_KINETIC_H5z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Fy2z_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_KINETIC_I6x_F3z_vrr = PAX*I_KINETIC_H5x_F3z_vrr+5*oned2z*I_KINETIC_G4x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I6x_F3z_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_KINETIC_I5xy_F3z_vrr = PAY*I_KINETIC_H5x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_F3z_vrr;
    Double I_KINETIC_I5xz_F3z_vrr = PAZ*I_KINETIC_H5x_F3z_vrr+3*oned2z*I_KINETIC_H5x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_F3z_vrr;
    Double I_KINETIC_I4x2y_F3z_vrr = PAY*I_KINETIC_H4xy_F3z_vrr+oned2z*I_KINETIC_G4x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_F3z_vrr-bdz*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_KINETIC_I4xyz_F3z_vrr = PAZ*I_KINETIC_H4xy_F3z_vrr+3*oned2z*I_KINETIC_H4xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_F3z_vrr;
    Double I_KINETIC_I4x2z_F3z_vrr = PAZ*I_KINETIC_H4xz_F3z_vrr+oned2z*I_KINETIC_G4x_F3z_vrr+3*oned2z*I_KINETIC_H4xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_KINETIC_I3x3y_F3z_vrr = PAY*I_KINETIC_H3x2y_F3z_vrr+2*oned2z*I_KINETIC_G3xy_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_F3z_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_F3z_vrr;
    Double I_KINETIC_I3x2yz_F3z_vrr = PAZ*I_KINETIC_H3x2y_F3z_vrr+3*oned2z*I_KINETIC_H3x2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_F3z_vrr;
    Double I_KINETIC_I3xy2z_F3z_vrr = PAY*I_KINETIC_H3x2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_F3z_vrr;
    Double I_KINETIC_I3x3z_F3z_vrr = PAZ*I_KINETIC_H3x2z_F3z_vrr+2*oned2z*I_KINETIC_G3xz_F3z_vrr+3*oned2z*I_KINETIC_H3x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_F3z_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_F3z_vrr;
    Double I_KINETIC_I2x4y_F3z_vrr = PAX*I_KINETIC_Hx4y_F3z_vrr+oned2z*I_KINETIC_G4y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_F3z_vrr-bdz*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_KINETIC_I2x3yz_F3z_vrr = PAZ*I_KINETIC_H2x3y_F3z_vrr+3*oned2z*I_KINETIC_H2x3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_F3z_vrr;
    Double I_KINETIC_I2x2y2z_F3z_vrr = PAZ*I_KINETIC_H2x2yz_F3z_vrr+oned2z*I_KINETIC_G2x2y_F3z_vrr+3*oned2z*I_KINETIC_H2x2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_F3z_vrr;
    Double I_KINETIC_I2xy3z_F3z_vrr = PAY*I_KINETIC_H2x3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_F3z_vrr;
    Double I_KINETIC_I2x4z_F3z_vrr = PAX*I_KINETIC_Hx4z_F3z_vrr+oned2z*I_KINETIC_G4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_F3z_vrr-bdz*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_KINETIC_Ix5y_F3z_vrr = PAX*I_KINETIC_H5y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_F3z_vrr;
    Double I_KINETIC_Ix4yz_F3z_vrr = PAZ*I_KINETIC_Hx4y_F3z_vrr+3*oned2z*I_KINETIC_Hx4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_F3z_vrr;
    Double I_KINETIC_Ix3y2z_F3z_vrr = PAX*I_KINETIC_H3y2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_F3z_vrr;
    Double I_KINETIC_Ix2y3z_F3z_vrr = PAX*I_KINETIC_H2y3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_F3z_vrr;
    Double I_KINETIC_Ixy4z_F3z_vrr = PAY*I_KINETIC_Hx4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_F3z_vrr;
    Double I_KINETIC_Ix5z_F3z_vrr = PAX*I_KINETIC_H5z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_F3z_vrr;
    Double I_KINETIC_I6y_F3z_vrr = PAY*I_KINETIC_H5y_F3z_vrr+5*oned2z*I_KINETIC_G4y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I6y_F3z_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_KINETIC_I5yz_F3z_vrr = PAZ*I_KINETIC_H5y_F3z_vrr+3*oned2z*I_KINETIC_H5y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_F3z_vrr;
    Double I_KINETIC_I4y2z_F3z_vrr = PAZ*I_KINETIC_H4yz_F3z_vrr+oned2z*I_KINETIC_G4y_F3z_vrr+3*oned2z*I_KINETIC_H4yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_KINETIC_I3y3z_F3z_vrr = PAZ*I_KINETIC_H3y2z_F3z_vrr+2*oned2z*I_KINETIC_G3yz_F3z_vrr+3*oned2z*I_KINETIC_H3y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_F3z_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_F3z_vrr;
    Double I_KINETIC_I2y4z_F3z_vrr = PAY*I_KINETIC_Hy4z_F3z_vrr+oned2z*I_KINETIC_G4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_F3z_vrr-bdz*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_KINETIC_Iy5z_F3z_vrr = PAY*I_KINETIC_H5z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_F3z_vrr;
    Double I_KINETIC_I6z_F3z_vrr = PAZ*I_KINETIC_H5z_F3z_vrr+5*oned2z*I_KINETIC_G4z_F3z_vrr+3*oned2z*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I6z_F3z_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_F_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_I_F_a_coefs = alpha;
    I_KINETIC_I6x_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_F3x_vrr;
    I_KINETIC_I5xy_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_F3x_vrr;
    I_KINETIC_I5xz_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_F3x_vrr;
    I_KINETIC_I4x2y_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_F3x_vrr;
    I_KINETIC_I4xyz_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_F3x_vrr;
    I_KINETIC_I4x2z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_F3x_vrr;
    I_KINETIC_I3x3y_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_F3x_vrr;
    I_KINETIC_I3x2yz_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_F3x_vrr;
    I_KINETIC_I3xy2z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_F3x_vrr;
    I_KINETIC_I3x3z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_F3x_vrr;
    I_KINETIC_I2x4y_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_F3x_vrr;
    I_KINETIC_I2x3yz_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_F3x_vrr;
    I_KINETIC_I2x2y2z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_F3x_vrr;
    I_KINETIC_I2xy3z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_F3x_vrr;
    I_KINETIC_I2x4z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_F3x_vrr;
    I_KINETIC_Ix5y_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_F3x_vrr;
    I_KINETIC_Ix4yz_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_F3x_vrr;
    I_KINETIC_Ix3y2z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_F3x_vrr;
    I_KINETIC_Ix2y3z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_F3x_vrr;
    I_KINETIC_Ixy4z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_F3x_vrr;
    I_KINETIC_Ix5z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_F3x_vrr;
    I_KINETIC_I6y_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_F3x_vrr;
    I_KINETIC_I5yz_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_F3x_vrr;
    I_KINETIC_I4y2z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_F3x_vrr;
    I_KINETIC_I3y3z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_F3x_vrr;
    I_KINETIC_I2y4z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_F3x_vrr;
    I_KINETIC_Iy5z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_F3x_vrr;
    I_KINETIC_I6z_F3x_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_F3x_vrr;
    I_KINETIC_I6x_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_F2xy_vrr;
    I_KINETIC_I5xy_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_F2xy_vrr;
    I_KINETIC_I5xz_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_F2xy_vrr;
    I_KINETIC_I4x2y_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_F2xy_vrr;
    I_KINETIC_I4xyz_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_F2xy_vrr;
    I_KINETIC_I4x2z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_F2xy_vrr;
    I_KINETIC_I3x3y_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_F2xy_vrr;
    I_KINETIC_I3x2yz_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_F2xy_vrr;
    I_KINETIC_I3xy2z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_F2xy_vrr;
    I_KINETIC_I3x3z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_F2xy_vrr;
    I_KINETIC_I2x4y_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_F2xy_vrr;
    I_KINETIC_I2x3yz_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_F2xy_vrr;
    I_KINETIC_I2x2y2z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_F2xy_vrr;
    I_KINETIC_I2xy3z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_F2xy_vrr;
    I_KINETIC_I2x4z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_F2xy_vrr;
    I_KINETIC_Ix5y_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_F2xy_vrr;
    I_KINETIC_Ix4yz_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_F2xy_vrr;
    I_KINETIC_Ix3y2z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_F2xy_vrr;
    I_KINETIC_Ix2y3z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_F2xy_vrr;
    I_KINETIC_Ixy4z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_F2xy_vrr;
    I_KINETIC_Ix5z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_F2xy_vrr;
    I_KINETIC_I6y_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_F2xy_vrr;
    I_KINETIC_I5yz_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_F2xy_vrr;
    I_KINETIC_I4y2z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_F2xy_vrr;
    I_KINETIC_I3y3z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_F2xy_vrr;
    I_KINETIC_I2y4z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_F2xy_vrr;
    I_KINETIC_Iy5z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_F2xy_vrr;
    I_KINETIC_I6z_F2xy_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_F2xy_vrr;
    I_KINETIC_I6x_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_F2xz_vrr;
    I_KINETIC_I5xy_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_F2xz_vrr;
    I_KINETIC_I5xz_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_F2xz_vrr;
    I_KINETIC_I4x2y_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_F2xz_vrr;
    I_KINETIC_I4xyz_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_F2xz_vrr;
    I_KINETIC_I4x2z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_F2xz_vrr;
    I_KINETIC_I3x3y_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_F2xz_vrr;
    I_KINETIC_I3x2yz_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_F2xz_vrr;
    I_KINETIC_I3xy2z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_F2xz_vrr;
    I_KINETIC_I3x3z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_F2xz_vrr;
    I_KINETIC_I2x4y_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_F2xz_vrr;
    I_KINETIC_I2x3yz_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_F2xz_vrr;
    I_KINETIC_I2x2y2z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_F2xz_vrr;
    I_KINETIC_I2xy3z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_F2xz_vrr;
    I_KINETIC_I2x4z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_F2xz_vrr;
    I_KINETIC_Ix5y_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_F2xz_vrr;
    I_KINETIC_Ix4yz_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_F2xz_vrr;
    I_KINETIC_Ix3y2z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_F2xz_vrr;
    I_KINETIC_Ix2y3z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_F2xz_vrr;
    I_KINETIC_Ixy4z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_F2xz_vrr;
    I_KINETIC_Ix5z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_F2xz_vrr;
    I_KINETIC_I6y_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_F2xz_vrr;
    I_KINETIC_I5yz_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_F2xz_vrr;
    I_KINETIC_I4y2z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_F2xz_vrr;
    I_KINETIC_I3y3z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_F2xz_vrr;
    I_KINETIC_I2y4z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_F2xz_vrr;
    I_KINETIC_Iy5z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_F2xz_vrr;
    I_KINETIC_I6z_F2xz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_F2xz_vrr;
    I_KINETIC_I6x_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_Fx2y_vrr;
    I_KINETIC_I5xy_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_Fx2y_vrr;
    I_KINETIC_I5xz_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_Fx2y_vrr;
    I_KINETIC_I4x2y_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_Fx2y_vrr;
    I_KINETIC_I4xyz_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_Fx2y_vrr;
    I_KINETIC_I4x2z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_Fx2y_vrr;
    I_KINETIC_I3x3y_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_Fx2y_vrr;
    I_KINETIC_I3x2yz_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_Fx2y_vrr;
    I_KINETIC_I3xy2z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_Fx2y_vrr;
    I_KINETIC_I3x3z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_Fx2y_vrr;
    I_KINETIC_I2x4y_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_Fx2y_vrr;
    I_KINETIC_I2x3yz_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_Fx2y_vrr;
    I_KINETIC_I2x2y2z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_Fx2y_vrr;
    I_KINETIC_I2xy3z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_Fx2y_vrr;
    I_KINETIC_I2x4z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_Fx2y_vrr;
    I_KINETIC_Ix5y_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_Fx2y_vrr;
    I_KINETIC_Ix4yz_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_Fx2y_vrr;
    I_KINETIC_Ix3y2z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_Fx2y_vrr;
    I_KINETIC_Ix2y3z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_Fx2y_vrr;
    I_KINETIC_Ixy4z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_Fx2y_vrr;
    I_KINETIC_Ix5z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_Fx2y_vrr;
    I_KINETIC_I6y_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_Fx2y_vrr;
    I_KINETIC_I5yz_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_Fx2y_vrr;
    I_KINETIC_I4y2z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_Fx2y_vrr;
    I_KINETIC_I3y3z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_Fx2y_vrr;
    I_KINETIC_I2y4z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_Fx2y_vrr;
    I_KINETIC_Iy5z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_Fx2y_vrr;
    I_KINETIC_I6z_Fx2y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_Fx2y_vrr;
    I_KINETIC_I6x_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_Fxyz_vrr;
    I_KINETIC_I5xy_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_Fxyz_vrr;
    I_KINETIC_I5xz_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_Fxyz_vrr;
    I_KINETIC_I4x2y_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_Fxyz_vrr;
    I_KINETIC_I4xyz_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_Fxyz_vrr;
    I_KINETIC_I4x2z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_Fxyz_vrr;
    I_KINETIC_I3x3y_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_Fxyz_vrr;
    I_KINETIC_I3x2yz_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_Fxyz_vrr;
    I_KINETIC_I3xy2z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_Fxyz_vrr;
    I_KINETIC_I3x3z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_Fxyz_vrr;
    I_KINETIC_I2x4y_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_Fxyz_vrr;
    I_KINETIC_I2x3yz_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_Fxyz_vrr;
    I_KINETIC_I2x2y2z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_Fxyz_vrr;
    I_KINETIC_I2xy3z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_Fxyz_vrr;
    I_KINETIC_I2x4z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_Fxyz_vrr;
    I_KINETIC_Ix5y_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_Fxyz_vrr;
    I_KINETIC_Ix4yz_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_Fxyz_vrr;
    I_KINETIC_Ix3y2z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_Fxyz_vrr;
    I_KINETIC_Ix2y3z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_Fxyz_vrr;
    I_KINETIC_Ixy4z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_Fxyz_vrr;
    I_KINETIC_Ix5z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_Fxyz_vrr;
    I_KINETIC_I6y_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_Fxyz_vrr;
    I_KINETIC_I5yz_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_Fxyz_vrr;
    I_KINETIC_I4y2z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_Fxyz_vrr;
    I_KINETIC_I3y3z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_Fxyz_vrr;
    I_KINETIC_I2y4z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_Fxyz_vrr;
    I_KINETIC_Iy5z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_Fxyz_vrr;
    I_KINETIC_I6z_Fxyz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_Fxyz_vrr;
    I_KINETIC_I6x_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_Fx2z_vrr;
    I_KINETIC_I5xy_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_Fx2z_vrr;
    I_KINETIC_I5xz_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_Fx2z_vrr;
    I_KINETIC_I4x2y_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_Fx2z_vrr;
    I_KINETIC_I4xyz_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_Fx2z_vrr;
    I_KINETIC_I4x2z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_Fx2z_vrr;
    I_KINETIC_I3x3y_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_Fx2z_vrr;
    I_KINETIC_I3x2yz_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_Fx2z_vrr;
    I_KINETIC_I3xy2z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_Fx2z_vrr;
    I_KINETIC_I3x3z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_Fx2z_vrr;
    I_KINETIC_I2x4y_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_Fx2z_vrr;
    I_KINETIC_I2x3yz_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_Fx2z_vrr;
    I_KINETIC_I2x2y2z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_Fx2z_vrr;
    I_KINETIC_I2xy3z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_Fx2z_vrr;
    I_KINETIC_I2x4z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_Fx2z_vrr;
    I_KINETIC_Ix5y_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_Fx2z_vrr;
    I_KINETIC_Ix4yz_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_Fx2z_vrr;
    I_KINETIC_Ix3y2z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_Fx2z_vrr;
    I_KINETIC_Ix2y3z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_Fx2z_vrr;
    I_KINETIC_Ixy4z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_Fx2z_vrr;
    I_KINETIC_Ix5z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_Fx2z_vrr;
    I_KINETIC_I6y_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_Fx2z_vrr;
    I_KINETIC_I5yz_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_Fx2z_vrr;
    I_KINETIC_I4y2z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_Fx2z_vrr;
    I_KINETIC_I3y3z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_Fx2z_vrr;
    I_KINETIC_I2y4z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_Fx2z_vrr;
    I_KINETIC_Iy5z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_Fx2z_vrr;
    I_KINETIC_I6z_Fx2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_Fx2z_vrr;
    I_KINETIC_I6x_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_F3y_vrr;
    I_KINETIC_I5xy_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_F3y_vrr;
    I_KINETIC_I5xz_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_F3y_vrr;
    I_KINETIC_I4x2y_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_F3y_vrr;
    I_KINETIC_I4xyz_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_F3y_vrr;
    I_KINETIC_I4x2z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_F3y_vrr;
    I_KINETIC_I3x3y_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_F3y_vrr;
    I_KINETIC_I3x2yz_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_F3y_vrr;
    I_KINETIC_I3xy2z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_F3y_vrr;
    I_KINETIC_I3x3z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_F3y_vrr;
    I_KINETIC_I2x4y_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_F3y_vrr;
    I_KINETIC_I2x3yz_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_F3y_vrr;
    I_KINETIC_I2x2y2z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_F3y_vrr;
    I_KINETIC_I2xy3z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_F3y_vrr;
    I_KINETIC_I2x4z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_F3y_vrr;
    I_KINETIC_Ix5y_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_F3y_vrr;
    I_KINETIC_Ix4yz_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_F3y_vrr;
    I_KINETIC_Ix3y2z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_F3y_vrr;
    I_KINETIC_Ix2y3z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_F3y_vrr;
    I_KINETIC_Ixy4z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_F3y_vrr;
    I_KINETIC_Ix5z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_F3y_vrr;
    I_KINETIC_I6y_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_F3y_vrr;
    I_KINETIC_I5yz_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_F3y_vrr;
    I_KINETIC_I4y2z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_F3y_vrr;
    I_KINETIC_I3y3z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_F3y_vrr;
    I_KINETIC_I2y4z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_F3y_vrr;
    I_KINETIC_Iy5z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_F3y_vrr;
    I_KINETIC_I6z_F3y_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_F3y_vrr;
    I_KINETIC_I6x_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_F2yz_vrr;
    I_KINETIC_I5xy_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_F2yz_vrr;
    I_KINETIC_I5xz_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_F2yz_vrr;
    I_KINETIC_I4x2y_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_F2yz_vrr;
    I_KINETIC_I4xyz_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_F2yz_vrr;
    I_KINETIC_I4x2z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_F2yz_vrr;
    I_KINETIC_I3x3y_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_F2yz_vrr;
    I_KINETIC_I3x2yz_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_F2yz_vrr;
    I_KINETIC_I3xy2z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_F2yz_vrr;
    I_KINETIC_I3x3z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_F2yz_vrr;
    I_KINETIC_I2x4y_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_F2yz_vrr;
    I_KINETIC_I2x3yz_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_F2yz_vrr;
    I_KINETIC_I2x2y2z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_F2yz_vrr;
    I_KINETIC_I2xy3z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_F2yz_vrr;
    I_KINETIC_I2x4z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_F2yz_vrr;
    I_KINETIC_Ix5y_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_F2yz_vrr;
    I_KINETIC_Ix4yz_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_F2yz_vrr;
    I_KINETIC_Ix3y2z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_F2yz_vrr;
    I_KINETIC_Ix2y3z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_F2yz_vrr;
    I_KINETIC_Ixy4z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_F2yz_vrr;
    I_KINETIC_Ix5z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_F2yz_vrr;
    I_KINETIC_I6y_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_F2yz_vrr;
    I_KINETIC_I5yz_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_F2yz_vrr;
    I_KINETIC_I4y2z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_F2yz_vrr;
    I_KINETIC_I3y3z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_F2yz_vrr;
    I_KINETIC_I2y4z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_F2yz_vrr;
    I_KINETIC_Iy5z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_F2yz_vrr;
    I_KINETIC_I6z_F2yz_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_F2yz_vrr;
    I_KINETIC_I6x_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_Fy2z_vrr;
    I_KINETIC_I5xy_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_Fy2z_vrr;
    I_KINETIC_I5xz_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_Fy2z_vrr;
    I_KINETIC_I4x2y_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_Fy2z_vrr;
    I_KINETIC_I4xyz_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_Fy2z_vrr;
    I_KINETIC_I4x2z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_Fy2z_vrr;
    I_KINETIC_I3x3y_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_Fy2z_vrr;
    I_KINETIC_I3x2yz_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_Fy2z_vrr;
    I_KINETIC_I3xy2z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_Fy2z_vrr;
    I_KINETIC_I3x3z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_Fy2z_vrr;
    I_KINETIC_I2x4y_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_Fy2z_vrr;
    I_KINETIC_I2x3yz_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_Fy2z_vrr;
    I_KINETIC_I2x2y2z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_Fy2z_vrr;
    I_KINETIC_I2xy3z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_Fy2z_vrr;
    I_KINETIC_I2x4z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_Fy2z_vrr;
    I_KINETIC_Ix5y_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_Fy2z_vrr;
    I_KINETIC_Ix4yz_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_Fy2z_vrr;
    I_KINETIC_Ix3y2z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_Fy2z_vrr;
    I_KINETIC_Ix2y3z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_Fy2z_vrr;
    I_KINETIC_Ixy4z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_Fy2z_vrr;
    I_KINETIC_Ix5z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_Fy2z_vrr;
    I_KINETIC_I6y_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_Fy2z_vrr;
    I_KINETIC_I5yz_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_Fy2z_vrr;
    I_KINETIC_I4y2z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_Fy2z_vrr;
    I_KINETIC_I3y3z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_Fy2z_vrr;
    I_KINETIC_I2y4z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_Fy2z_vrr;
    I_KINETIC_Iy5z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_Fy2z_vrr;
    I_KINETIC_I6z_Fy2z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_Fy2z_vrr;
    I_KINETIC_I6x_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6x_F3z_vrr;
    I_KINETIC_I5xy_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xy_F3z_vrr;
    I_KINETIC_I5xz_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5xz_F3z_vrr;
    I_KINETIC_I4x2y_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2y_F3z_vrr;
    I_KINETIC_I4xyz_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4xyz_F3z_vrr;
    I_KINETIC_I4x2z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4x2z_F3z_vrr;
    I_KINETIC_I3x3y_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3y_F3z_vrr;
    I_KINETIC_I3x2yz_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x2yz_F3z_vrr;
    I_KINETIC_I3xy2z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3xy2z_F3z_vrr;
    I_KINETIC_I3x3z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3x3z_F3z_vrr;
    I_KINETIC_I2x4y_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4y_F3z_vrr;
    I_KINETIC_I2x3yz_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x3yz_F3z_vrr;
    I_KINETIC_I2x2y2z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x2y2z_F3z_vrr;
    I_KINETIC_I2xy3z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2xy3z_F3z_vrr;
    I_KINETIC_I2x4z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2x4z_F3z_vrr;
    I_KINETIC_Ix5y_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5y_F3z_vrr;
    I_KINETIC_Ix4yz_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix4yz_F3z_vrr;
    I_KINETIC_Ix3y2z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix3y2z_F3z_vrr;
    I_KINETIC_Ix2y3z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix2y3z_F3z_vrr;
    I_KINETIC_Ixy4z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ixy4z_F3z_vrr;
    I_KINETIC_Ix5z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Ix5z_F3z_vrr;
    I_KINETIC_I6y_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6y_F3z_vrr;
    I_KINETIC_I5yz_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I5yz_F3z_vrr;
    I_KINETIC_I4y2z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I4y2z_F3z_vrr;
    I_KINETIC_I3y3z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I3y3z_F3z_vrr;
    I_KINETIC_I2y4z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I2y4z_F3z_vrr;
    I_KINETIC_Iy5z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_Iy5z_F3z_vrr;
    I_KINETIC_I6z_F3z_a += SQ_KINETIC_I_F_a_coefs*I_KINETIC_I6z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_F
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_G4x_F3x += I_KINETIC_G4x_F3x_vrr;
    I_KINETIC_G3xy_F3x += I_KINETIC_G3xy_F3x_vrr;
    I_KINETIC_G3xz_F3x += I_KINETIC_G3xz_F3x_vrr;
    I_KINETIC_G2x2y_F3x += I_KINETIC_G2x2y_F3x_vrr;
    I_KINETIC_G2xyz_F3x += I_KINETIC_G2xyz_F3x_vrr;
    I_KINETIC_G2x2z_F3x += I_KINETIC_G2x2z_F3x_vrr;
    I_KINETIC_Gx3y_F3x += I_KINETIC_Gx3y_F3x_vrr;
    I_KINETIC_Gx2yz_F3x += I_KINETIC_Gx2yz_F3x_vrr;
    I_KINETIC_Gxy2z_F3x += I_KINETIC_Gxy2z_F3x_vrr;
    I_KINETIC_Gx3z_F3x += I_KINETIC_Gx3z_F3x_vrr;
    I_KINETIC_G4y_F3x += I_KINETIC_G4y_F3x_vrr;
    I_KINETIC_G3yz_F3x += I_KINETIC_G3yz_F3x_vrr;
    I_KINETIC_G2y2z_F3x += I_KINETIC_G2y2z_F3x_vrr;
    I_KINETIC_Gy3z_F3x += I_KINETIC_Gy3z_F3x_vrr;
    I_KINETIC_G4z_F3x += I_KINETIC_G4z_F3x_vrr;
    I_KINETIC_G4x_F2xy += I_KINETIC_G4x_F2xy_vrr;
    I_KINETIC_G3xy_F2xy += I_KINETIC_G3xy_F2xy_vrr;
    I_KINETIC_G3xz_F2xy += I_KINETIC_G3xz_F2xy_vrr;
    I_KINETIC_G2x2y_F2xy += I_KINETIC_G2x2y_F2xy_vrr;
    I_KINETIC_G2xyz_F2xy += I_KINETIC_G2xyz_F2xy_vrr;
    I_KINETIC_G2x2z_F2xy += I_KINETIC_G2x2z_F2xy_vrr;
    I_KINETIC_Gx3y_F2xy += I_KINETIC_Gx3y_F2xy_vrr;
    I_KINETIC_Gx2yz_F2xy += I_KINETIC_Gx2yz_F2xy_vrr;
    I_KINETIC_Gxy2z_F2xy += I_KINETIC_Gxy2z_F2xy_vrr;
    I_KINETIC_Gx3z_F2xy += I_KINETIC_Gx3z_F2xy_vrr;
    I_KINETIC_G4y_F2xy += I_KINETIC_G4y_F2xy_vrr;
    I_KINETIC_G3yz_F2xy += I_KINETIC_G3yz_F2xy_vrr;
    I_KINETIC_G2y2z_F2xy += I_KINETIC_G2y2z_F2xy_vrr;
    I_KINETIC_Gy3z_F2xy += I_KINETIC_Gy3z_F2xy_vrr;
    I_KINETIC_G4z_F2xy += I_KINETIC_G4z_F2xy_vrr;
    I_KINETIC_G4x_F2xz += I_KINETIC_G4x_F2xz_vrr;
    I_KINETIC_G3xy_F2xz += I_KINETIC_G3xy_F2xz_vrr;
    I_KINETIC_G3xz_F2xz += I_KINETIC_G3xz_F2xz_vrr;
    I_KINETIC_G2x2y_F2xz += I_KINETIC_G2x2y_F2xz_vrr;
    I_KINETIC_G2xyz_F2xz += I_KINETIC_G2xyz_F2xz_vrr;
    I_KINETIC_G2x2z_F2xz += I_KINETIC_G2x2z_F2xz_vrr;
    I_KINETIC_Gx3y_F2xz += I_KINETIC_Gx3y_F2xz_vrr;
    I_KINETIC_Gx2yz_F2xz += I_KINETIC_Gx2yz_F2xz_vrr;
    I_KINETIC_Gxy2z_F2xz += I_KINETIC_Gxy2z_F2xz_vrr;
    I_KINETIC_Gx3z_F2xz += I_KINETIC_Gx3z_F2xz_vrr;
    I_KINETIC_G4y_F2xz += I_KINETIC_G4y_F2xz_vrr;
    I_KINETIC_G3yz_F2xz += I_KINETIC_G3yz_F2xz_vrr;
    I_KINETIC_G2y2z_F2xz += I_KINETIC_G2y2z_F2xz_vrr;
    I_KINETIC_Gy3z_F2xz += I_KINETIC_Gy3z_F2xz_vrr;
    I_KINETIC_G4z_F2xz += I_KINETIC_G4z_F2xz_vrr;
    I_KINETIC_G4x_Fx2y += I_KINETIC_G4x_Fx2y_vrr;
    I_KINETIC_G3xy_Fx2y += I_KINETIC_G3xy_Fx2y_vrr;
    I_KINETIC_G3xz_Fx2y += I_KINETIC_G3xz_Fx2y_vrr;
    I_KINETIC_G2x2y_Fx2y += I_KINETIC_G2x2y_Fx2y_vrr;
    I_KINETIC_G2xyz_Fx2y += I_KINETIC_G2xyz_Fx2y_vrr;
    I_KINETIC_G2x2z_Fx2y += I_KINETIC_G2x2z_Fx2y_vrr;
    I_KINETIC_Gx3y_Fx2y += I_KINETIC_Gx3y_Fx2y_vrr;
    I_KINETIC_Gx2yz_Fx2y += I_KINETIC_Gx2yz_Fx2y_vrr;
    I_KINETIC_Gxy2z_Fx2y += I_KINETIC_Gxy2z_Fx2y_vrr;
    I_KINETIC_Gx3z_Fx2y += I_KINETIC_Gx3z_Fx2y_vrr;
    I_KINETIC_G4y_Fx2y += I_KINETIC_G4y_Fx2y_vrr;
    I_KINETIC_G3yz_Fx2y += I_KINETIC_G3yz_Fx2y_vrr;
    I_KINETIC_G2y2z_Fx2y += I_KINETIC_G2y2z_Fx2y_vrr;
    I_KINETIC_Gy3z_Fx2y += I_KINETIC_Gy3z_Fx2y_vrr;
    I_KINETIC_G4z_Fx2y += I_KINETIC_G4z_Fx2y_vrr;
    I_KINETIC_G4x_Fxyz += I_KINETIC_G4x_Fxyz_vrr;
    I_KINETIC_G3xy_Fxyz += I_KINETIC_G3xy_Fxyz_vrr;
    I_KINETIC_G3xz_Fxyz += I_KINETIC_G3xz_Fxyz_vrr;
    I_KINETIC_G2x2y_Fxyz += I_KINETIC_G2x2y_Fxyz_vrr;
    I_KINETIC_G2xyz_Fxyz += I_KINETIC_G2xyz_Fxyz_vrr;
    I_KINETIC_G2x2z_Fxyz += I_KINETIC_G2x2z_Fxyz_vrr;
    I_KINETIC_Gx3y_Fxyz += I_KINETIC_Gx3y_Fxyz_vrr;
    I_KINETIC_Gx2yz_Fxyz += I_KINETIC_Gx2yz_Fxyz_vrr;
    I_KINETIC_Gxy2z_Fxyz += I_KINETIC_Gxy2z_Fxyz_vrr;
    I_KINETIC_Gx3z_Fxyz += I_KINETIC_Gx3z_Fxyz_vrr;
    I_KINETIC_G4y_Fxyz += I_KINETIC_G4y_Fxyz_vrr;
    I_KINETIC_G3yz_Fxyz += I_KINETIC_G3yz_Fxyz_vrr;
    I_KINETIC_G2y2z_Fxyz += I_KINETIC_G2y2z_Fxyz_vrr;
    I_KINETIC_Gy3z_Fxyz += I_KINETIC_Gy3z_Fxyz_vrr;
    I_KINETIC_G4z_Fxyz += I_KINETIC_G4z_Fxyz_vrr;
    I_KINETIC_G4x_Fx2z += I_KINETIC_G4x_Fx2z_vrr;
    I_KINETIC_G3xy_Fx2z += I_KINETIC_G3xy_Fx2z_vrr;
    I_KINETIC_G3xz_Fx2z += I_KINETIC_G3xz_Fx2z_vrr;
    I_KINETIC_G2x2y_Fx2z += I_KINETIC_G2x2y_Fx2z_vrr;
    I_KINETIC_G2xyz_Fx2z += I_KINETIC_G2xyz_Fx2z_vrr;
    I_KINETIC_G2x2z_Fx2z += I_KINETIC_G2x2z_Fx2z_vrr;
    I_KINETIC_Gx3y_Fx2z += I_KINETIC_Gx3y_Fx2z_vrr;
    I_KINETIC_Gx2yz_Fx2z += I_KINETIC_Gx2yz_Fx2z_vrr;
    I_KINETIC_Gxy2z_Fx2z += I_KINETIC_Gxy2z_Fx2z_vrr;
    I_KINETIC_Gx3z_Fx2z += I_KINETIC_Gx3z_Fx2z_vrr;
    I_KINETIC_G4y_Fx2z += I_KINETIC_G4y_Fx2z_vrr;
    I_KINETIC_G3yz_Fx2z += I_KINETIC_G3yz_Fx2z_vrr;
    I_KINETIC_G2y2z_Fx2z += I_KINETIC_G2y2z_Fx2z_vrr;
    I_KINETIC_Gy3z_Fx2z += I_KINETIC_Gy3z_Fx2z_vrr;
    I_KINETIC_G4z_Fx2z += I_KINETIC_G4z_Fx2z_vrr;
    I_KINETIC_G4x_F3y += I_KINETIC_G4x_F3y_vrr;
    I_KINETIC_G3xy_F3y += I_KINETIC_G3xy_F3y_vrr;
    I_KINETIC_G3xz_F3y += I_KINETIC_G3xz_F3y_vrr;
    I_KINETIC_G2x2y_F3y += I_KINETIC_G2x2y_F3y_vrr;
    I_KINETIC_G2xyz_F3y += I_KINETIC_G2xyz_F3y_vrr;
    I_KINETIC_G2x2z_F3y += I_KINETIC_G2x2z_F3y_vrr;
    I_KINETIC_Gx3y_F3y += I_KINETIC_Gx3y_F3y_vrr;
    I_KINETIC_Gx2yz_F3y += I_KINETIC_Gx2yz_F3y_vrr;
    I_KINETIC_Gxy2z_F3y += I_KINETIC_Gxy2z_F3y_vrr;
    I_KINETIC_Gx3z_F3y += I_KINETIC_Gx3z_F3y_vrr;
    I_KINETIC_G4y_F3y += I_KINETIC_G4y_F3y_vrr;
    I_KINETIC_G3yz_F3y += I_KINETIC_G3yz_F3y_vrr;
    I_KINETIC_G2y2z_F3y += I_KINETIC_G2y2z_F3y_vrr;
    I_KINETIC_Gy3z_F3y += I_KINETIC_Gy3z_F3y_vrr;
    I_KINETIC_G4z_F3y += I_KINETIC_G4z_F3y_vrr;
    I_KINETIC_G4x_F2yz += I_KINETIC_G4x_F2yz_vrr;
    I_KINETIC_G3xy_F2yz += I_KINETIC_G3xy_F2yz_vrr;
    I_KINETIC_G3xz_F2yz += I_KINETIC_G3xz_F2yz_vrr;
    I_KINETIC_G2x2y_F2yz += I_KINETIC_G2x2y_F2yz_vrr;
    I_KINETIC_G2xyz_F2yz += I_KINETIC_G2xyz_F2yz_vrr;
    I_KINETIC_G2x2z_F2yz += I_KINETIC_G2x2z_F2yz_vrr;
    I_KINETIC_Gx3y_F2yz += I_KINETIC_Gx3y_F2yz_vrr;
    I_KINETIC_Gx2yz_F2yz += I_KINETIC_Gx2yz_F2yz_vrr;
    I_KINETIC_Gxy2z_F2yz += I_KINETIC_Gxy2z_F2yz_vrr;
    I_KINETIC_Gx3z_F2yz += I_KINETIC_Gx3z_F2yz_vrr;
    I_KINETIC_G4y_F2yz += I_KINETIC_G4y_F2yz_vrr;
    I_KINETIC_G3yz_F2yz += I_KINETIC_G3yz_F2yz_vrr;
    I_KINETIC_G2y2z_F2yz += I_KINETIC_G2y2z_F2yz_vrr;
    I_KINETIC_Gy3z_F2yz += I_KINETIC_Gy3z_F2yz_vrr;
    I_KINETIC_G4z_F2yz += I_KINETIC_G4z_F2yz_vrr;
    I_KINETIC_G4x_Fy2z += I_KINETIC_G4x_Fy2z_vrr;
    I_KINETIC_G3xy_Fy2z += I_KINETIC_G3xy_Fy2z_vrr;
    I_KINETIC_G3xz_Fy2z += I_KINETIC_G3xz_Fy2z_vrr;
    I_KINETIC_G2x2y_Fy2z += I_KINETIC_G2x2y_Fy2z_vrr;
    I_KINETIC_G2xyz_Fy2z += I_KINETIC_G2xyz_Fy2z_vrr;
    I_KINETIC_G2x2z_Fy2z += I_KINETIC_G2x2z_Fy2z_vrr;
    I_KINETIC_Gx3y_Fy2z += I_KINETIC_Gx3y_Fy2z_vrr;
    I_KINETIC_Gx2yz_Fy2z += I_KINETIC_Gx2yz_Fy2z_vrr;
    I_KINETIC_Gxy2z_Fy2z += I_KINETIC_Gxy2z_Fy2z_vrr;
    I_KINETIC_Gx3z_Fy2z += I_KINETIC_Gx3z_Fy2z_vrr;
    I_KINETIC_G4y_Fy2z += I_KINETIC_G4y_Fy2z_vrr;
    I_KINETIC_G3yz_Fy2z += I_KINETIC_G3yz_Fy2z_vrr;
    I_KINETIC_G2y2z_Fy2z += I_KINETIC_G2y2z_Fy2z_vrr;
    I_KINETIC_Gy3z_Fy2z += I_KINETIC_Gy3z_Fy2z_vrr;
    I_KINETIC_G4z_Fy2z += I_KINETIC_G4z_Fy2z_vrr;
    I_KINETIC_G4x_F3z += I_KINETIC_G4x_F3z_vrr;
    I_KINETIC_G3xy_F3z += I_KINETIC_G3xy_F3z_vrr;
    I_KINETIC_G3xz_F3z += I_KINETIC_G3xz_F3z_vrr;
    I_KINETIC_G2x2y_F3z += I_KINETIC_G2x2y_F3z_vrr;
    I_KINETIC_G2xyz_F3z += I_KINETIC_G2xyz_F3z_vrr;
    I_KINETIC_G2x2z_F3z += I_KINETIC_G2x2z_F3z_vrr;
    I_KINETIC_Gx3y_F3z += I_KINETIC_Gx3y_F3z_vrr;
    I_KINETIC_Gx2yz_F3z += I_KINETIC_Gx2yz_F3z_vrr;
    I_KINETIC_Gxy2z_F3z += I_KINETIC_Gxy2z_F3z_vrr;
    I_KINETIC_Gx3z_F3z += I_KINETIC_Gx3z_F3z_vrr;
    I_KINETIC_G4y_F3z += I_KINETIC_G4y_F3z_vrr;
    I_KINETIC_G3yz_F3z += I_KINETIC_G3yz_F3z_vrr;
    I_KINETIC_G2y2z_F3z += I_KINETIC_G2y2z_F3z_vrr;
    I_KINETIC_Gy3z_F3z += I_KINETIC_Gy3z_F3z_vrr;
    I_KINETIC_G4z_F3z += I_KINETIC_G4z_F3z_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_F_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_F_a
   * RHS shell quartet name: SQ_KINETIC_G_F
   ************************************************************/
  abcd[0] = 2.0E0*I_KINETIC_I6x_F3x_a-5*I_KINETIC_G4x_F3x;
  abcd[1] = 2.0E0*I_KINETIC_I5xy_F3x_a-4*I_KINETIC_G3xy_F3x;
  abcd[2] = 2.0E0*I_KINETIC_I5xz_F3x_a-4*I_KINETIC_G3xz_F3x;
  abcd[3] = 2.0E0*I_KINETIC_I4x2y_F3x_a-3*I_KINETIC_G2x2y_F3x;
  abcd[4] = 2.0E0*I_KINETIC_I4xyz_F3x_a-3*I_KINETIC_G2xyz_F3x;
  abcd[5] = 2.0E0*I_KINETIC_I4x2z_F3x_a-3*I_KINETIC_G2x2z_F3x;
  abcd[6] = 2.0E0*I_KINETIC_I3x3y_F3x_a-2*I_KINETIC_Gx3y_F3x;
  abcd[7] = 2.0E0*I_KINETIC_I3x2yz_F3x_a-2*I_KINETIC_Gx2yz_F3x;
  abcd[8] = 2.0E0*I_KINETIC_I3xy2z_F3x_a-2*I_KINETIC_Gxy2z_F3x;
  abcd[9] = 2.0E0*I_KINETIC_I3x3z_F3x_a-2*I_KINETIC_Gx3z_F3x;
  abcd[10] = 2.0E0*I_KINETIC_I2x4y_F3x_a-1*I_KINETIC_G4y_F3x;
  abcd[11] = 2.0E0*I_KINETIC_I2x3yz_F3x_a-1*I_KINETIC_G3yz_F3x;
  abcd[12] = 2.0E0*I_KINETIC_I2x2y2z_F3x_a-1*I_KINETIC_G2y2z_F3x;
  abcd[13] = 2.0E0*I_KINETIC_I2xy3z_F3x_a-1*I_KINETIC_Gy3z_F3x;
  abcd[14] = 2.0E0*I_KINETIC_I2x4z_F3x_a-1*I_KINETIC_G4z_F3x;
  abcd[15] = 2.0E0*I_KINETIC_Ix5y_F3x_a;
  abcd[16] = 2.0E0*I_KINETIC_Ix4yz_F3x_a;
  abcd[17] = 2.0E0*I_KINETIC_Ix3y2z_F3x_a;
  abcd[18] = 2.0E0*I_KINETIC_Ix2y3z_F3x_a;
  abcd[19] = 2.0E0*I_KINETIC_Ixy4z_F3x_a;
  abcd[20] = 2.0E0*I_KINETIC_Ix5z_F3x_a;
  abcd[21] = 2.0E0*I_KINETIC_I6x_F2xy_a-5*I_KINETIC_G4x_F2xy;
  abcd[22] = 2.0E0*I_KINETIC_I5xy_F2xy_a-4*I_KINETIC_G3xy_F2xy;
  abcd[23] = 2.0E0*I_KINETIC_I5xz_F2xy_a-4*I_KINETIC_G3xz_F2xy;
  abcd[24] = 2.0E0*I_KINETIC_I4x2y_F2xy_a-3*I_KINETIC_G2x2y_F2xy;
  abcd[25] = 2.0E0*I_KINETIC_I4xyz_F2xy_a-3*I_KINETIC_G2xyz_F2xy;
  abcd[26] = 2.0E0*I_KINETIC_I4x2z_F2xy_a-3*I_KINETIC_G2x2z_F2xy;
  abcd[27] = 2.0E0*I_KINETIC_I3x3y_F2xy_a-2*I_KINETIC_Gx3y_F2xy;
  abcd[28] = 2.0E0*I_KINETIC_I3x2yz_F2xy_a-2*I_KINETIC_Gx2yz_F2xy;
  abcd[29] = 2.0E0*I_KINETIC_I3xy2z_F2xy_a-2*I_KINETIC_Gxy2z_F2xy;
  abcd[30] = 2.0E0*I_KINETIC_I3x3z_F2xy_a-2*I_KINETIC_Gx3z_F2xy;
  abcd[31] = 2.0E0*I_KINETIC_I2x4y_F2xy_a-1*I_KINETIC_G4y_F2xy;
  abcd[32] = 2.0E0*I_KINETIC_I2x3yz_F2xy_a-1*I_KINETIC_G3yz_F2xy;
  abcd[33] = 2.0E0*I_KINETIC_I2x2y2z_F2xy_a-1*I_KINETIC_G2y2z_F2xy;
  abcd[34] = 2.0E0*I_KINETIC_I2xy3z_F2xy_a-1*I_KINETIC_Gy3z_F2xy;
  abcd[35] = 2.0E0*I_KINETIC_I2x4z_F2xy_a-1*I_KINETIC_G4z_F2xy;
  abcd[36] = 2.0E0*I_KINETIC_Ix5y_F2xy_a;
  abcd[37] = 2.0E0*I_KINETIC_Ix4yz_F2xy_a;
  abcd[38] = 2.0E0*I_KINETIC_Ix3y2z_F2xy_a;
  abcd[39] = 2.0E0*I_KINETIC_Ix2y3z_F2xy_a;
  abcd[40] = 2.0E0*I_KINETIC_Ixy4z_F2xy_a;
  abcd[41] = 2.0E0*I_KINETIC_Ix5z_F2xy_a;
  abcd[42] = 2.0E0*I_KINETIC_I6x_F2xz_a-5*I_KINETIC_G4x_F2xz;
  abcd[43] = 2.0E0*I_KINETIC_I5xy_F2xz_a-4*I_KINETIC_G3xy_F2xz;
  abcd[44] = 2.0E0*I_KINETIC_I5xz_F2xz_a-4*I_KINETIC_G3xz_F2xz;
  abcd[45] = 2.0E0*I_KINETIC_I4x2y_F2xz_a-3*I_KINETIC_G2x2y_F2xz;
  abcd[46] = 2.0E0*I_KINETIC_I4xyz_F2xz_a-3*I_KINETIC_G2xyz_F2xz;
  abcd[47] = 2.0E0*I_KINETIC_I4x2z_F2xz_a-3*I_KINETIC_G2x2z_F2xz;
  abcd[48] = 2.0E0*I_KINETIC_I3x3y_F2xz_a-2*I_KINETIC_Gx3y_F2xz;
  abcd[49] = 2.0E0*I_KINETIC_I3x2yz_F2xz_a-2*I_KINETIC_Gx2yz_F2xz;
  abcd[50] = 2.0E0*I_KINETIC_I3xy2z_F2xz_a-2*I_KINETIC_Gxy2z_F2xz;
  abcd[51] = 2.0E0*I_KINETIC_I3x3z_F2xz_a-2*I_KINETIC_Gx3z_F2xz;
  abcd[52] = 2.0E0*I_KINETIC_I2x4y_F2xz_a-1*I_KINETIC_G4y_F2xz;
  abcd[53] = 2.0E0*I_KINETIC_I2x3yz_F2xz_a-1*I_KINETIC_G3yz_F2xz;
  abcd[54] = 2.0E0*I_KINETIC_I2x2y2z_F2xz_a-1*I_KINETIC_G2y2z_F2xz;
  abcd[55] = 2.0E0*I_KINETIC_I2xy3z_F2xz_a-1*I_KINETIC_Gy3z_F2xz;
  abcd[56] = 2.0E0*I_KINETIC_I2x4z_F2xz_a-1*I_KINETIC_G4z_F2xz;
  abcd[57] = 2.0E0*I_KINETIC_Ix5y_F2xz_a;
  abcd[58] = 2.0E0*I_KINETIC_Ix4yz_F2xz_a;
  abcd[59] = 2.0E0*I_KINETIC_Ix3y2z_F2xz_a;
  abcd[60] = 2.0E0*I_KINETIC_Ix2y3z_F2xz_a;
  abcd[61] = 2.0E0*I_KINETIC_Ixy4z_F2xz_a;
  abcd[62] = 2.0E0*I_KINETIC_Ix5z_F2xz_a;
  abcd[63] = 2.0E0*I_KINETIC_I6x_Fx2y_a-5*I_KINETIC_G4x_Fx2y;
  abcd[64] = 2.0E0*I_KINETIC_I5xy_Fx2y_a-4*I_KINETIC_G3xy_Fx2y;
  abcd[65] = 2.0E0*I_KINETIC_I5xz_Fx2y_a-4*I_KINETIC_G3xz_Fx2y;
  abcd[66] = 2.0E0*I_KINETIC_I4x2y_Fx2y_a-3*I_KINETIC_G2x2y_Fx2y;
  abcd[67] = 2.0E0*I_KINETIC_I4xyz_Fx2y_a-3*I_KINETIC_G2xyz_Fx2y;
  abcd[68] = 2.0E0*I_KINETIC_I4x2z_Fx2y_a-3*I_KINETIC_G2x2z_Fx2y;
  abcd[69] = 2.0E0*I_KINETIC_I3x3y_Fx2y_a-2*I_KINETIC_Gx3y_Fx2y;
  abcd[70] = 2.0E0*I_KINETIC_I3x2yz_Fx2y_a-2*I_KINETIC_Gx2yz_Fx2y;
  abcd[71] = 2.0E0*I_KINETIC_I3xy2z_Fx2y_a-2*I_KINETIC_Gxy2z_Fx2y;
  abcd[72] = 2.0E0*I_KINETIC_I3x3z_Fx2y_a-2*I_KINETIC_Gx3z_Fx2y;
  abcd[73] = 2.0E0*I_KINETIC_I2x4y_Fx2y_a-1*I_KINETIC_G4y_Fx2y;
  abcd[74] = 2.0E0*I_KINETIC_I2x3yz_Fx2y_a-1*I_KINETIC_G3yz_Fx2y;
  abcd[75] = 2.0E0*I_KINETIC_I2x2y2z_Fx2y_a-1*I_KINETIC_G2y2z_Fx2y;
  abcd[76] = 2.0E0*I_KINETIC_I2xy3z_Fx2y_a-1*I_KINETIC_Gy3z_Fx2y;
  abcd[77] = 2.0E0*I_KINETIC_I2x4z_Fx2y_a-1*I_KINETIC_G4z_Fx2y;
  abcd[78] = 2.0E0*I_KINETIC_Ix5y_Fx2y_a;
  abcd[79] = 2.0E0*I_KINETIC_Ix4yz_Fx2y_a;
  abcd[80] = 2.0E0*I_KINETIC_Ix3y2z_Fx2y_a;
  abcd[81] = 2.0E0*I_KINETIC_Ix2y3z_Fx2y_a;
  abcd[82] = 2.0E0*I_KINETIC_Ixy4z_Fx2y_a;
  abcd[83] = 2.0E0*I_KINETIC_Ix5z_Fx2y_a;
  abcd[84] = 2.0E0*I_KINETIC_I6x_Fxyz_a-5*I_KINETIC_G4x_Fxyz;
  abcd[85] = 2.0E0*I_KINETIC_I5xy_Fxyz_a-4*I_KINETIC_G3xy_Fxyz;
  abcd[86] = 2.0E0*I_KINETIC_I5xz_Fxyz_a-4*I_KINETIC_G3xz_Fxyz;
  abcd[87] = 2.0E0*I_KINETIC_I4x2y_Fxyz_a-3*I_KINETIC_G2x2y_Fxyz;
  abcd[88] = 2.0E0*I_KINETIC_I4xyz_Fxyz_a-3*I_KINETIC_G2xyz_Fxyz;
  abcd[89] = 2.0E0*I_KINETIC_I4x2z_Fxyz_a-3*I_KINETIC_G2x2z_Fxyz;
  abcd[90] = 2.0E0*I_KINETIC_I3x3y_Fxyz_a-2*I_KINETIC_Gx3y_Fxyz;
  abcd[91] = 2.0E0*I_KINETIC_I3x2yz_Fxyz_a-2*I_KINETIC_Gx2yz_Fxyz;
  abcd[92] = 2.0E0*I_KINETIC_I3xy2z_Fxyz_a-2*I_KINETIC_Gxy2z_Fxyz;
  abcd[93] = 2.0E0*I_KINETIC_I3x3z_Fxyz_a-2*I_KINETIC_Gx3z_Fxyz;
  abcd[94] = 2.0E0*I_KINETIC_I2x4y_Fxyz_a-1*I_KINETIC_G4y_Fxyz;
  abcd[95] = 2.0E0*I_KINETIC_I2x3yz_Fxyz_a-1*I_KINETIC_G3yz_Fxyz;
  abcd[96] = 2.0E0*I_KINETIC_I2x2y2z_Fxyz_a-1*I_KINETIC_G2y2z_Fxyz;
  abcd[97] = 2.0E0*I_KINETIC_I2xy3z_Fxyz_a-1*I_KINETIC_Gy3z_Fxyz;
  abcd[98] = 2.0E0*I_KINETIC_I2x4z_Fxyz_a-1*I_KINETIC_G4z_Fxyz;
  abcd[99] = 2.0E0*I_KINETIC_Ix5y_Fxyz_a;
  abcd[100] = 2.0E0*I_KINETIC_Ix4yz_Fxyz_a;
  abcd[101] = 2.0E0*I_KINETIC_Ix3y2z_Fxyz_a;
  abcd[102] = 2.0E0*I_KINETIC_Ix2y3z_Fxyz_a;
  abcd[103] = 2.0E0*I_KINETIC_Ixy4z_Fxyz_a;
  abcd[104] = 2.0E0*I_KINETIC_Ix5z_Fxyz_a;
  abcd[105] = 2.0E0*I_KINETIC_I6x_Fx2z_a-5*I_KINETIC_G4x_Fx2z;
  abcd[106] = 2.0E0*I_KINETIC_I5xy_Fx2z_a-4*I_KINETIC_G3xy_Fx2z;
  abcd[107] = 2.0E0*I_KINETIC_I5xz_Fx2z_a-4*I_KINETIC_G3xz_Fx2z;
  abcd[108] = 2.0E0*I_KINETIC_I4x2y_Fx2z_a-3*I_KINETIC_G2x2y_Fx2z;
  abcd[109] = 2.0E0*I_KINETIC_I4xyz_Fx2z_a-3*I_KINETIC_G2xyz_Fx2z;
  abcd[110] = 2.0E0*I_KINETIC_I4x2z_Fx2z_a-3*I_KINETIC_G2x2z_Fx2z;
  abcd[111] = 2.0E0*I_KINETIC_I3x3y_Fx2z_a-2*I_KINETIC_Gx3y_Fx2z;
  abcd[112] = 2.0E0*I_KINETIC_I3x2yz_Fx2z_a-2*I_KINETIC_Gx2yz_Fx2z;
  abcd[113] = 2.0E0*I_KINETIC_I3xy2z_Fx2z_a-2*I_KINETIC_Gxy2z_Fx2z;
  abcd[114] = 2.0E0*I_KINETIC_I3x3z_Fx2z_a-2*I_KINETIC_Gx3z_Fx2z;
  abcd[115] = 2.0E0*I_KINETIC_I2x4y_Fx2z_a-1*I_KINETIC_G4y_Fx2z;
  abcd[116] = 2.0E0*I_KINETIC_I2x3yz_Fx2z_a-1*I_KINETIC_G3yz_Fx2z;
  abcd[117] = 2.0E0*I_KINETIC_I2x2y2z_Fx2z_a-1*I_KINETIC_G2y2z_Fx2z;
  abcd[118] = 2.0E0*I_KINETIC_I2xy3z_Fx2z_a-1*I_KINETIC_Gy3z_Fx2z;
  abcd[119] = 2.0E0*I_KINETIC_I2x4z_Fx2z_a-1*I_KINETIC_G4z_Fx2z;
  abcd[120] = 2.0E0*I_KINETIC_Ix5y_Fx2z_a;
  abcd[121] = 2.0E0*I_KINETIC_Ix4yz_Fx2z_a;
  abcd[122] = 2.0E0*I_KINETIC_Ix3y2z_Fx2z_a;
  abcd[123] = 2.0E0*I_KINETIC_Ix2y3z_Fx2z_a;
  abcd[124] = 2.0E0*I_KINETIC_Ixy4z_Fx2z_a;
  abcd[125] = 2.0E0*I_KINETIC_Ix5z_Fx2z_a;
  abcd[126] = 2.0E0*I_KINETIC_I6x_F3y_a-5*I_KINETIC_G4x_F3y;
  abcd[127] = 2.0E0*I_KINETIC_I5xy_F3y_a-4*I_KINETIC_G3xy_F3y;
  abcd[128] = 2.0E0*I_KINETIC_I5xz_F3y_a-4*I_KINETIC_G3xz_F3y;
  abcd[129] = 2.0E0*I_KINETIC_I4x2y_F3y_a-3*I_KINETIC_G2x2y_F3y;
  abcd[130] = 2.0E0*I_KINETIC_I4xyz_F3y_a-3*I_KINETIC_G2xyz_F3y;
  abcd[131] = 2.0E0*I_KINETIC_I4x2z_F3y_a-3*I_KINETIC_G2x2z_F3y;
  abcd[132] = 2.0E0*I_KINETIC_I3x3y_F3y_a-2*I_KINETIC_Gx3y_F3y;
  abcd[133] = 2.0E0*I_KINETIC_I3x2yz_F3y_a-2*I_KINETIC_Gx2yz_F3y;
  abcd[134] = 2.0E0*I_KINETIC_I3xy2z_F3y_a-2*I_KINETIC_Gxy2z_F3y;
  abcd[135] = 2.0E0*I_KINETIC_I3x3z_F3y_a-2*I_KINETIC_Gx3z_F3y;
  abcd[136] = 2.0E0*I_KINETIC_I2x4y_F3y_a-1*I_KINETIC_G4y_F3y;
  abcd[137] = 2.0E0*I_KINETIC_I2x3yz_F3y_a-1*I_KINETIC_G3yz_F3y;
  abcd[138] = 2.0E0*I_KINETIC_I2x2y2z_F3y_a-1*I_KINETIC_G2y2z_F3y;
  abcd[139] = 2.0E0*I_KINETIC_I2xy3z_F3y_a-1*I_KINETIC_Gy3z_F3y;
  abcd[140] = 2.0E0*I_KINETIC_I2x4z_F3y_a-1*I_KINETIC_G4z_F3y;
  abcd[141] = 2.0E0*I_KINETIC_Ix5y_F3y_a;
  abcd[142] = 2.0E0*I_KINETIC_Ix4yz_F3y_a;
  abcd[143] = 2.0E0*I_KINETIC_Ix3y2z_F3y_a;
  abcd[144] = 2.0E0*I_KINETIC_Ix2y3z_F3y_a;
  abcd[145] = 2.0E0*I_KINETIC_Ixy4z_F3y_a;
  abcd[146] = 2.0E0*I_KINETIC_Ix5z_F3y_a;
  abcd[147] = 2.0E0*I_KINETIC_I6x_F2yz_a-5*I_KINETIC_G4x_F2yz;
  abcd[148] = 2.0E0*I_KINETIC_I5xy_F2yz_a-4*I_KINETIC_G3xy_F2yz;
  abcd[149] = 2.0E0*I_KINETIC_I5xz_F2yz_a-4*I_KINETIC_G3xz_F2yz;
  abcd[150] = 2.0E0*I_KINETIC_I4x2y_F2yz_a-3*I_KINETIC_G2x2y_F2yz;
  abcd[151] = 2.0E0*I_KINETIC_I4xyz_F2yz_a-3*I_KINETIC_G2xyz_F2yz;
  abcd[152] = 2.0E0*I_KINETIC_I4x2z_F2yz_a-3*I_KINETIC_G2x2z_F2yz;
  abcd[153] = 2.0E0*I_KINETIC_I3x3y_F2yz_a-2*I_KINETIC_Gx3y_F2yz;
  abcd[154] = 2.0E0*I_KINETIC_I3x2yz_F2yz_a-2*I_KINETIC_Gx2yz_F2yz;
  abcd[155] = 2.0E0*I_KINETIC_I3xy2z_F2yz_a-2*I_KINETIC_Gxy2z_F2yz;
  abcd[156] = 2.0E0*I_KINETIC_I3x3z_F2yz_a-2*I_KINETIC_Gx3z_F2yz;
  abcd[157] = 2.0E0*I_KINETIC_I2x4y_F2yz_a-1*I_KINETIC_G4y_F2yz;
  abcd[158] = 2.0E0*I_KINETIC_I2x3yz_F2yz_a-1*I_KINETIC_G3yz_F2yz;
  abcd[159] = 2.0E0*I_KINETIC_I2x2y2z_F2yz_a-1*I_KINETIC_G2y2z_F2yz;
  abcd[160] = 2.0E0*I_KINETIC_I2xy3z_F2yz_a-1*I_KINETIC_Gy3z_F2yz;
  abcd[161] = 2.0E0*I_KINETIC_I2x4z_F2yz_a-1*I_KINETIC_G4z_F2yz;
  abcd[162] = 2.0E0*I_KINETIC_Ix5y_F2yz_a;
  abcd[163] = 2.0E0*I_KINETIC_Ix4yz_F2yz_a;
  abcd[164] = 2.0E0*I_KINETIC_Ix3y2z_F2yz_a;
  abcd[165] = 2.0E0*I_KINETIC_Ix2y3z_F2yz_a;
  abcd[166] = 2.0E0*I_KINETIC_Ixy4z_F2yz_a;
  abcd[167] = 2.0E0*I_KINETIC_Ix5z_F2yz_a;
  abcd[168] = 2.0E0*I_KINETIC_I6x_Fy2z_a-5*I_KINETIC_G4x_Fy2z;
  abcd[169] = 2.0E0*I_KINETIC_I5xy_Fy2z_a-4*I_KINETIC_G3xy_Fy2z;
  abcd[170] = 2.0E0*I_KINETIC_I5xz_Fy2z_a-4*I_KINETIC_G3xz_Fy2z;
  abcd[171] = 2.0E0*I_KINETIC_I4x2y_Fy2z_a-3*I_KINETIC_G2x2y_Fy2z;
  abcd[172] = 2.0E0*I_KINETIC_I4xyz_Fy2z_a-3*I_KINETIC_G2xyz_Fy2z;
  abcd[173] = 2.0E0*I_KINETIC_I4x2z_Fy2z_a-3*I_KINETIC_G2x2z_Fy2z;
  abcd[174] = 2.0E0*I_KINETIC_I3x3y_Fy2z_a-2*I_KINETIC_Gx3y_Fy2z;
  abcd[175] = 2.0E0*I_KINETIC_I3x2yz_Fy2z_a-2*I_KINETIC_Gx2yz_Fy2z;
  abcd[176] = 2.0E0*I_KINETIC_I3xy2z_Fy2z_a-2*I_KINETIC_Gxy2z_Fy2z;
  abcd[177] = 2.0E0*I_KINETIC_I3x3z_Fy2z_a-2*I_KINETIC_Gx3z_Fy2z;
  abcd[178] = 2.0E0*I_KINETIC_I2x4y_Fy2z_a-1*I_KINETIC_G4y_Fy2z;
  abcd[179] = 2.0E0*I_KINETIC_I2x3yz_Fy2z_a-1*I_KINETIC_G3yz_Fy2z;
  abcd[180] = 2.0E0*I_KINETIC_I2x2y2z_Fy2z_a-1*I_KINETIC_G2y2z_Fy2z;
  abcd[181] = 2.0E0*I_KINETIC_I2xy3z_Fy2z_a-1*I_KINETIC_Gy3z_Fy2z;
  abcd[182] = 2.0E0*I_KINETIC_I2x4z_Fy2z_a-1*I_KINETIC_G4z_Fy2z;
  abcd[183] = 2.0E0*I_KINETIC_Ix5y_Fy2z_a;
  abcd[184] = 2.0E0*I_KINETIC_Ix4yz_Fy2z_a;
  abcd[185] = 2.0E0*I_KINETIC_Ix3y2z_Fy2z_a;
  abcd[186] = 2.0E0*I_KINETIC_Ix2y3z_Fy2z_a;
  abcd[187] = 2.0E0*I_KINETIC_Ixy4z_Fy2z_a;
  abcd[188] = 2.0E0*I_KINETIC_Ix5z_Fy2z_a;
  abcd[189] = 2.0E0*I_KINETIC_I6x_F3z_a-5*I_KINETIC_G4x_F3z;
  abcd[190] = 2.0E0*I_KINETIC_I5xy_F3z_a-4*I_KINETIC_G3xy_F3z;
  abcd[191] = 2.0E0*I_KINETIC_I5xz_F3z_a-4*I_KINETIC_G3xz_F3z;
  abcd[192] = 2.0E0*I_KINETIC_I4x2y_F3z_a-3*I_KINETIC_G2x2y_F3z;
  abcd[193] = 2.0E0*I_KINETIC_I4xyz_F3z_a-3*I_KINETIC_G2xyz_F3z;
  abcd[194] = 2.0E0*I_KINETIC_I4x2z_F3z_a-3*I_KINETIC_G2x2z_F3z;
  abcd[195] = 2.0E0*I_KINETIC_I3x3y_F3z_a-2*I_KINETIC_Gx3y_F3z;
  abcd[196] = 2.0E0*I_KINETIC_I3x2yz_F3z_a-2*I_KINETIC_Gx2yz_F3z;
  abcd[197] = 2.0E0*I_KINETIC_I3xy2z_F3z_a-2*I_KINETIC_Gxy2z_F3z;
  abcd[198] = 2.0E0*I_KINETIC_I3x3z_F3z_a-2*I_KINETIC_Gx3z_F3z;
  abcd[199] = 2.0E0*I_KINETIC_I2x4y_F3z_a-1*I_KINETIC_G4y_F3z;
  abcd[200] = 2.0E0*I_KINETIC_I2x3yz_F3z_a-1*I_KINETIC_G3yz_F3z;
  abcd[201] = 2.0E0*I_KINETIC_I2x2y2z_F3z_a-1*I_KINETIC_G2y2z_F3z;
  abcd[202] = 2.0E0*I_KINETIC_I2xy3z_F3z_a-1*I_KINETIC_Gy3z_F3z;
  abcd[203] = 2.0E0*I_KINETIC_I2x4z_F3z_a-1*I_KINETIC_G4z_F3z;
  abcd[204] = 2.0E0*I_KINETIC_Ix5y_F3z_a;
  abcd[205] = 2.0E0*I_KINETIC_Ix4yz_F3z_a;
  abcd[206] = 2.0E0*I_KINETIC_Ix3y2z_F3z_a;
  abcd[207] = 2.0E0*I_KINETIC_Ix2y3z_F3z_a;
  abcd[208] = 2.0E0*I_KINETIC_Ixy4z_F3z_a;
  abcd[209] = 2.0E0*I_KINETIC_Ix5z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_F_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_F_a
   * RHS shell quartet name: SQ_KINETIC_G_F
   ************************************************************/
  abcd[210] = 2.0E0*I_KINETIC_I5xy_F3x_a;
  abcd[211] = 2.0E0*I_KINETIC_I4x2y_F3x_a-1*I_KINETIC_G4x_F3x;
  abcd[212] = 2.0E0*I_KINETIC_I4xyz_F3x_a;
  abcd[213] = 2.0E0*I_KINETIC_I3x3y_F3x_a-2*I_KINETIC_G3xy_F3x;
  abcd[214] = 2.0E0*I_KINETIC_I3x2yz_F3x_a-1*I_KINETIC_G3xz_F3x;
  abcd[215] = 2.0E0*I_KINETIC_I3xy2z_F3x_a;
  abcd[216] = 2.0E0*I_KINETIC_I2x4y_F3x_a-3*I_KINETIC_G2x2y_F3x;
  abcd[217] = 2.0E0*I_KINETIC_I2x3yz_F3x_a-2*I_KINETIC_G2xyz_F3x;
  abcd[218] = 2.0E0*I_KINETIC_I2x2y2z_F3x_a-1*I_KINETIC_G2x2z_F3x;
  abcd[219] = 2.0E0*I_KINETIC_I2xy3z_F3x_a;
  abcd[220] = 2.0E0*I_KINETIC_Ix5y_F3x_a-4*I_KINETIC_Gx3y_F3x;
  abcd[221] = 2.0E0*I_KINETIC_Ix4yz_F3x_a-3*I_KINETIC_Gx2yz_F3x;
  abcd[222] = 2.0E0*I_KINETIC_Ix3y2z_F3x_a-2*I_KINETIC_Gxy2z_F3x;
  abcd[223] = 2.0E0*I_KINETIC_Ix2y3z_F3x_a-1*I_KINETIC_Gx3z_F3x;
  abcd[224] = 2.0E0*I_KINETIC_Ixy4z_F3x_a;
  abcd[225] = 2.0E0*I_KINETIC_I6y_F3x_a-5*I_KINETIC_G4y_F3x;
  abcd[226] = 2.0E0*I_KINETIC_I5yz_F3x_a-4*I_KINETIC_G3yz_F3x;
  abcd[227] = 2.0E0*I_KINETIC_I4y2z_F3x_a-3*I_KINETIC_G2y2z_F3x;
  abcd[228] = 2.0E0*I_KINETIC_I3y3z_F3x_a-2*I_KINETIC_Gy3z_F3x;
  abcd[229] = 2.0E0*I_KINETIC_I2y4z_F3x_a-1*I_KINETIC_G4z_F3x;
  abcd[230] = 2.0E0*I_KINETIC_Iy5z_F3x_a;
  abcd[231] = 2.0E0*I_KINETIC_I5xy_F2xy_a;
  abcd[232] = 2.0E0*I_KINETIC_I4x2y_F2xy_a-1*I_KINETIC_G4x_F2xy;
  abcd[233] = 2.0E0*I_KINETIC_I4xyz_F2xy_a;
  abcd[234] = 2.0E0*I_KINETIC_I3x3y_F2xy_a-2*I_KINETIC_G3xy_F2xy;
  abcd[235] = 2.0E0*I_KINETIC_I3x2yz_F2xy_a-1*I_KINETIC_G3xz_F2xy;
  abcd[236] = 2.0E0*I_KINETIC_I3xy2z_F2xy_a;
  abcd[237] = 2.0E0*I_KINETIC_I2x4y_F2xy_a-3*I_KINETIC_G2x2y_F2xy;
  abcd[238] = 2.0E0*I_KINETIC_I2x3yz_F2xy_a-2*I_KINETIC_G2xyz_F2xy;
  abcd[239] = 2.0E0*I_KINETIC_I2x2y2z_F2xy_a-1*I_KINETIC_G2x2z_F2xy;
  abcd[240] = 2.0E0*I_KINETIC_I2xy3z_F2xy_a;
  abcd[241] = 2.0E0*I_KINETIC_Ix5y_F2xy_a-4*I_KINETIC_Gx3y_F2xy;
  abcd[242] = 2.0E0*I_KINETIC_Ix4yz_F2xy_a-3*I_KINETIC_Gx2yz_F2xy;
  abcd[243] = 2.0E0*I_KINETIC_Ix3y2z_F2xy_a-2*I_KINETIC_Gxy2z_F2xy;
  abcd[244] = 2.0E0*I_KINETIC_Ix2y3z_F2xy_a-1*I_KINETIC_Gx3z_F2xy;
  abcd[245] = 2.0E0*I_KINETIC_Ixy4z_F2xy_a;
  abcd[246] = 2.0E0*I_KINETIC_I6y_F2xy_a-5*I_KINETIC_G4y_F2xy;
  abcd[247] = 2.0E0*I_KINETIC_I5yz_F2xy_a-4*I_KINETIC_G3yz_F2xy;
  abcd[248] = 2.0E0*I_KINETIC_I4y2z_F2xy_a-3*I_KINETIC_G2y2z_F2xy;
  abcd[249] = 2.0E0*I_KINETIC_I3y3z_F2xy_a-2*I_KINETIC_Gy3z_F2xy;
  abcd[250] = 2.0E0*I_KINETIC_I2y4z_F2xy_a-1*I_KINETIC_G4z_F2xy;
  abcd[251] = 2.0E0*I_KINETIC_Iy5z_F2xy_a;
  abcd[252] = 2.0E0*I_KINETIC_I5xy_F2xz_a;
  abcd[253] = 2.0E0*I_KINETIC_I4x2y_F2xz_a-1*I_KINETIC_G4x_F2xz;
  abcd[254] = 2.0E0*I_KINETIC_I4xyz_F2xz_a;
  abcd[255] = 2.0E0*I_KINETIC_I3x3y_F2xz_a-2*I_KINETIC_G3xy_F2xz;
  abcd[256] = 2.0E0*I_KINETIC_I3x2yz_F2xz_a-1*I_KINETIC_G3xz_F2xz;
  abcd[257] = 2.0E0*I_KINETIC_I3xy2z_F2xz_a;
  abcd[258] = 2.0E0*I_KINETIC_I2x4y_F2xz_a-3*I_KINETIC_G2x2y_F2xz;
  abcd[259] = 2.0E0*I_KINETIC_I2x3yz_F2xz_a-2*I_KINETIC_G2xyz_F2xz;
  abcd[260] = 2.0E0*I_KINETIC_I2x2y2z_F2xz_a-1*I_KINETIC_G2x2z_F2xz;
  abcd[261] = 2.0E0*I_KINETIC_I2xy3z_F2xz_a;
  abcd[262] = 2.0E0*I_KINETIC_Ix5y_F2xz_a-4*I_KINETIC_Gx3y_F2xz;
  abcd[263] = 2.0E0*I_KINETIC_Ix4yz_F2xz_a-3*I_KINETIC_Gx2yz_F2xz;
  abcd[264] = 2.0E0*I_KINETIC_Ix3y2z_F2xz_a-2*I_KINETIC_Gxy2z_F2xz;
  abcd[265] = 2.0E0*I_KINETIC_Ix2y3z_F2xz_a-1*I_KINETIC_Gx3z_F2xz;
  abcd[266] = 2.0E0*I_KINETIC_Ixy4z_F2xz_a;
  abcd[267] = 2.0E0*I_KINETIC_I6y_F2xz_a-5*I_KINETIC_G4y_F2xz;
  abcd[268] = 2.0E0*I_KINETIC_I5yz_F2xz_a-4*I_KINETIC_G3yz_F2xz;
  abcd[269] = 2.0E0*I_KINETIC_I4y2z_F2xz_a-3*I_KINETIC_G2y2z_F2xz;
  abcd[270] = 2.0E0*I_KINETIC_I3y3z_F2xz_a-2*I_KINETIC_Gy3z_F2xz;
  abcd[271] = 2.0E0*I_KINETIC_I2y4z_F2xz_a-1*I_KINETIC_G4z_F2xz;
  abcd[272] = 2.0E0*I_KINETIC_Iy5z_F2xz_a;
  abcd[273] = 2.0E0*I_KINETIC_I5xy_Fx2y_a;
  abcd[274] = 2.0E0*I_KINETIC_I4x2y_Fx2y_a-1*I_KINETIC_G4x_Fx2y;
  abcd[275] = 2.0E0*I_KINETIC_I4xyz_Fx2y_a;
  abcd[276] = 2.0E0*I_KINETIC_I3x3y_Fx2y_a-2*I_KINETIC_G3xy_Fx2y;
  abcd[277] = 2.0E0*I_KINETIC_I3x2yz_Fx2y_a-1*I_KINETIC_G3xz_Fx2y;
  abcd[278] = 2.0E0*I_KINETIC_I3xy2z_Fx2y_a;
  abcd[279] = 2.0E0*I_KINETIC_I2x4y_Fx2y_a-3*I_KINETIC_G2x2y_Fx2y;
  abcd[280] = 2.0E0*I_KINETIC_I2x3yz_Fx2y_a-2*I_KINETIC_G2xyz_Fx2y;
  abcd[281] = 2.0E0*I_KINETIC_I2x2y2z_Fx2y_a-1*I_KINETIC_G2x2z_Fx2y;
  abcd[282] = 2.0E0*I_KINETIC_I2xy3z_Fx2y_a;
  abcd[283] = 2.0E0*I_KINETIC_Ix5y_Fx2y_a-4*I_KINETIC_Gx3y_Fx2y;
  abcd[284] = 2.0E0*I_KINETIC_Ix4yz_Fx2y_a-3*I_KINETIC_Gx2yz_Fx2y;
  abcd[285] = 2.0E0*I_KINETIC_Ix3y2z_Fx2y_a-2*I_KINETIC_Gxy2z_Fx2y;
  abcd[286] = 2.0E0*I_KINETIC_Ix2y3z_Fx2y_a-1*I_KINETIC_Gx3z_Fx2y;
  abcd[287] = 2.0E0*I_KINETIC_Ixy4z_Fx2y_a;
  abcd[288] = 2.0E0*I_KINETIC_I6y_Fx2y_a-5*I_KINETIC_G4y_Fx2y;
  abcd[289] = 2.0E0*I_KINETIC_I5yz_Fx2y_a-4*I_KINETIC_G3yz_Fx2y;
  abcd[290] = 2.0E0*I_KINETIC_I4y2z_Fx2y_a-3*I_KINETIC_G2y2z_Fx2y;
  abcd[291] = 2.0E0*I_KINETIC_I3y3z_Fx2y_a-2*I_KINETIC_Gy3z_Fx2y;
  abcd[292] = 2.0E0*I_KINETIC_I2y4z_Fx2y_a-1*I_KINETIC_G4z_Fx2y;
  abcd[293] = 2.0E0*I_KINETIC_Iy5z_Fx2y_a;
  abcd[294] = 2.0E0*I_KINETIC_I5xy_Fxyz_a;
  abcd[295] = 2.0E0*I_KINETIC_I4x2y_Fxyz_a-1*I_KINETIC_G4x_Fxyz;
  abcd[296] = 2.0E0*I_KINETIC_I4xyz_Fxyz_a;
  abcd[297] = 2.0E0*I_KINETIC_I3x3y_Fxyz_a-2*I_KINETIC_G3xy_Fxyz;
  abcd[298] = 2.0E0*I_KINETIC_I3x2yz_Fxyz_a-1*I_KINETIC_G3xz_Fxyz;
  abcd[299] = 2.0E0*I_KINETIC_I3xy2z_Fxyz_a;
  abcd[300] = 2.0E0*I_KINETIC_I2x4y_Fxyz_a-3*I_KINETIC_G2x2y_Fxyz;
  abcd[301] = 2.0E0*I_KINETIC_I2x3yz_Fxyz_a-2*I_KINETIC_G2xyz_Fxyz;
  abcd[302] = 2.0E0*I_KINETIC_I2x2y2z_Fxyz_a-1*I_KINETIC_G2x2z_Fxyz;
  abcd[303] = 2.0E0*I_KINETIC_I2xy3z_Fxyz_a;
  abcd[304] = 2.0E0*I_KINETIC_Ix5y_Fxyz_a-4*I_KINETIC_Gx3y_Fxyz;
  abcd[305] = 2.0E0*I_KINETIC_Ix4yz_Fxyz_a-3*I_KINETIC_Gx2yz_Fxyz;
  abcd[306] = 2.0E0*I_KINETIC_Ix3y2z_Fxyz_a-2*I_KINETIC_Gxy2z_Fxyz;
  abcd[307] = 2.0E0*I_KINETIC_Ix2y3z_Fxyz_a-1*I_KINETIC_Gx3z_Fxyz;
  abcd[308] = 2.0E0*I_KINETIC_Ixy4z_Fxyz_a;
  abcd[309] = 2.0E0*I_KINETIC_I6y_Fxyz_a-5*I_KINETIC_G4y_Fxyz;
  abcd[310] = 2.0E0*I_KINETIC_I5yz_Fxyz_a-4*I_KINETIC_G3yz_Fxyz;
  abcd[311] = 2.0E0*I_KINETIC_I4y2z_Fxyz_a-3*I_KINETIC_G2y2z_Fxyz;
  abcd[312] = 2.0E0*I_KINETIC_I3y3z_Fxyz_a-2*I_KINETIC_Gy3z_Fxyz;
  abcd[313] = 2.0E0*I_KINETIC_I2y4z_Fxyz_a-1*I_KINETIC_G4z_Fxyz;
  abcd[314] = 2.0E0*I_KINETIC_Iy5z_Fxyz_a;
  abcd[315] = 2.0E0*I_KINETIC_I5xy_Fx2z_a;
  abcd[316] = 2.0E0*I_KINETIC_I4x2y_Fx2z_a-1*I_KINETIC_G4x_Fx2z;
  abcd[317] = 2.0E0*I_KINETIC_I4xyz_Fx2z_a;
  abcd[318] = 2.0E0*I_KINETIC_I3x3y_Fx2z_a-2*I_KINETIC_G3xy_Fx2z;
  abcd[319] = 2.0E0*I_KINETIC_I3x2yz_Fx2z_a-1*I_KINETIC_G3xz_Fx2z;
  abcd[320] = 2.0E0*I_KINETIC_I3xy2z_Fx2z_a;
  abcd[321] = 2.0E0*I_KINETIC_I2x4y_Fx2z_a-3*I_KINETIC_G2x2y_Fx2z;
  abcd[322] = 2.0E0*I_KINETIC_I2x3yz_Fx2z_a-2*I_KINETIC_G2xyz_Fx2z;
  abcd[323] = 2.0E0*I_KINETIC_I2x2y2z_Fx2z_a-1*I_KINETIC_G2x2z_Fx2z;
  abcd[324] = 2.0E0*I_KINETIC_I2xy3z_Fx2z_a;
  abcd[325] = 2.0E0*I_KINETIC_Ix5y_Fx2z_a-4*I_KINETIC_Gx3y_Fx2z;
  abcd[326] = 2.0E0*I_KINETIC_Ix4yz_Fx2z_a-3*I_KINETIC_Gx2yz_Fx2z;
  abcd[327] = 2.0E0*I_KINETIC_Ix3y2z_Fx2z_a-2*I_KINETIC_Gxy2z_Fx2z;
  abcd[328] = 2.0E0*I_KINETIC_Ix2y3z_Fx2z_a-1*I_KINETIC_Gx3z_Fx2z;
  abcd[329] = 2.0E0*I_KINETIC_Ixy4z_Fx2z_a;
  abcd[330] = 2.0E0*I_KINETIC_I6y_Fx2z_a-5*I_KINETIC_G4y_Fx2z;
  abcd[331] = 2.0E0*I_KINETIC_I5yz_Fx2z_a-4*I_KINETIC_G3yz_Fx2z;
  abcd[332] = 2.0E0*I_KINETIC_I4y2z_Fx2z_a-3*I_KINETIC_G2y2z_Fx2z;
  abcd[333] = 2.0E0*I_KINETIC_I3y3z_Fx2z_a-2*I_KINETIC_Gy3z_Fx2z;
  abcd[334] = 2.0E0*I_KINETIC_I2y4z_Fx2z_a-1*I_KINETIC_G4z_Fx2z;
  abcd[335] = 2.0E0*I_KINETIC_Iy5z_Fx2z_a;
  abcd[336] = 2.0E0*I_KINETIC_I5xy_F3y_a;
  abcd[337] = 2.0E0*I_KINETIC_I4x2y_F3y_a-1*I_KINETIC_G4x_F3y;
  abcd[338] = 2.0E0*I_KINETIC_I4xyz_F3y_a;
  abcd[339] = 2.0E0*I_KINETIC_I3x3y_F3y_a-2*I_KINETIC_G3xy_F3y;
  abcd[340] = 2.0E0*I_KINETIC_I3x2yz_F3y_a-1*I_KINETIC_G3xz_F3y;
  abcd[341] = 2.0E0*I_KINETIC_I3xy2z_F3y_a;
  abcd[342] = 2.0E0*I_KINETIC_I2x4y_F3y_a-3*I_KINETIC_G2x2y_F3y;
  abcd[343] = 2.0E0*I_KINETIC_I2x3yz_F3y_a-2*I_KINETIC_G2xyz_F3y;
  abcd[344] = 2.0E0*I_KINETIC_I2x2y2z_F3y_a-1*I_KINETIC_G2x2z_F3y;
  abcd[345] = 2.0E0*I_KINETIC_I2xy3z_F3y_a;
  abcd[346] = 2.0E0*I_KINETIC_Ix5y_F3y_a-4*I_KINETIC_Gx3y_F3y;
  abcd[347] = 2.0E0*I_KINETIC_Ix4yz_F3y_a-3*I_KINETIC_Gx2yz_F3y;
  abcd[348] = 2.0E0*I_KINETIC_Ix3y2z_F3y_a-2*I_KINETIC_Gxy2z_F3y;
  abcd[349] = 2.0E0*I_KINETIC_Ix2y3z_F3y_a-1*I_KINETIC_Gx3z_F3y;
  abcd[350] = 2.0E0*I_KINETIC_Ixy4z_F3y_a;
  abcd[351] = 2.0E0*I_KINETIC_I6y_F3y_a-5*I_KINETIC_G4y_F3y;
  abcd[352] = 2.0E0*I_KINETIC_I5yz_F3y_a-4*I_KINETIC_G3yz_F3y;
  abcd[353] = 2.0E0*I_KINETIC_I4y2z_F3y_a-3*I_KINETIC_G2y2z_F3y;
  abcd[354] = 2.0E0*I_KINETIC_I3y3z_F3y_a-2*I_KINETIC_Gy3z_F3y;
  abcd[355] = 2.0E0*I_KINETIC_I2y4z_F3y_a-1*I_KINETIC_G4z_F3y;
  abcd[356] = 2.0E0*I_KINETIC_Iy5z_F3y_a;
  abcd[357] = 2.0E0*I_KINETIC_I5xy_F2yz_a;
  abcd[358] = 2.0E0*I_KINETIC_I4x2y_F2yz_a-1*I_KINETIC_G4x_F2yz;
  abcd[359] = 2.0E0*I_KINETIC_I4xyz_F2yz_a;
  abcd[360] = 2.0E0*I_KINETIC_I3x3y_F2yz_a-2*I_KINETIC_G3xy_F2yz;
  abcd[361] = 2.0E0*I_KINETIC_I3x2yz_F2yz_a-1*I_KINETIC_G3xz_F2yz;
  abcd[362] = 2.0E0*I_KINETIC_I3xy2z_F2yz_a;
  abcd[363] = 2.0E0*I_KINETIC_I2x4y_F2yz_a-3*I_KINETIC_G2x2y_F2yz;
  abcd[364] = 2.0E0*I_KINETIC_I2x3yz_F2yz_a-2*I_KINETIC_G2xyz_F2yz;
  abcd[365] = 2.0E0*I_KINETIC_I2x2y2z_F2yz_a-1*I_KINETIC_G2x2z_F2yz;
  abcd[366] = 2.0E0*I_KINETIC_I2xy3z_F2yz_a;
  abcd[367] = 2.0E0*I_KINETIC_Ix5y_F2yz_a-4*I_KINETIC_Gx3y_F2yz;
  abcd[368] = 2.0E0*I_KINETIC_Ix4yz_F2yz_a-3*I_KINETIC_Gx2yz_F2yz;
  abcd[369] = 2.0E0*I_KINETIC_Ix3y2z_F2yz_a-2*I_KINETIC_Gxy2z_F2yz;
  abcd[370] = 2.0E0*I_KINETIC_Ix2y3z_F2yz_a-1*I_KINETIC_Gx3z_F2yz;
  abcd[371] = 2.0E0*I_KINETIC_Ixy4z_F2yz_a;
  abcd[372] = 2.0E0*I_KINETIC_I6y_F2yz_a-5*I_KINETIC_G4y_F2yz;
  abcd[373] = 2.0E0*I_KINETIC_I5yz_F2yz_a-4*I_KINETIC_G3yz_F2yz;
  abcd[374] = 2.0E0*I_KINETIC_I4y2z_F2yz_a-3*I_KINETIC_G2y2z_F2yz;
  abcd[375] = 2.0E0*I_KINETIC_I3y3z_F2yz_a-2*I_KINETIC_Gy3z_F2yz;
  abcd[376] = 2.0E0*I_KINETIC_I2y4z_F2yz_a-1*I_KINETIC_G4z_F2yz;
  abcd[377] = 2.0E0*I_KINETIC_Iy5z_F2yz_a;
  abcd[378] = 2.0E0*I_KINETIC_I5xy_Fy2z_a;
  abcd[379] = 2.0E0*I_KINETIC_I4x2y_Fy2z_a-1*I_KINETIC_G4x_Fy2z;
  abcd[380] = 2.0E0*I_KINETIC_I4xyz_Fy2z_a;
  abcd[381] = 2.0E0*I_KINETIC_I3x3y_Fy2z_a-2*I_KINETIC_G3xy_Fy2z;
  abcd[382] = 2.0E0*I_KINETIC_I3x2yz_Fy2z_a-1*I_KINETIC_G3xz_Fy2z;
  abcd[383] = 2.0E0*I_KINETIC_I3xy2z_Fy2z_a;
  abcd[384] = 2.0E0*I_KINETIC_I2x4y_Fy2z_a-3*I_KINETIC_G2x2y_Fy2z;
  abcd[385] = 2.0E0*I_KINETIC_I2x3yz_Fy2z_a-2*I_KINETIC_G2xyz_Fy2z;
  abcd[386] = 2.0E0*I_KINETIC_I2x2y2z_Fy2z_a-1*I_KINETIC_G2x2z_Fy2z;
  abcd[387] = 2.0E0*I_KINETIC_I2xy3z_Fy2z_a;
  abcd[388] = 2.0E0*I_KINETIC_Ix5y_Fy2z_a-4*I_KINETIC_Gx3y_Fy2z;
  abcd[389] = 2.0E0*I_KINETIC_Ix4yz_Fy2z_a-3*I_KINETIC_Gx2yz_Fy2z;
  abcd[390] = 2.0E0*I_KINETIC_Ix3y2z_Fy2z_a-2*I_KINETIC_Gxy2z_Fy2z;
  abcd[391] = 2.0E0*I_KINETIC_Ix2y3z_Fy2z_a-1*I_KINETIC_Gx3z_Fy2z;
  abcd[392] = 2.0E0*I_KINETIC_Ixy4z_Fy2z_a;
  abcd[393] = 2.0E0*I_KINETIC_I6y_Fy2z_a-5*I_KINETIC_G4y_Fy2z;
  abcd[394] = 2.0E0*I_KINETIC_I5yz_Fy2z_a-4*I_KINETIC_G3yz_Fy2z;
  abcd[395] = 2.0E0*I_KINETIC_I4y2z_Fy2z_a-3*I_KINETIC_G2y2z_Fy2z;
  abcd[396] = 2.0E0*I_KINETIC_I3y3z_Fy2z_a-2*I_KINETIC_Gy3z_Fy2z;
  abcd[397] = 2.0E0*I_KINETIC_I2y4z_Fy2z_a-1*I_KINETIC_G4z_Fy2z;
  abcd[398] = 2.0E0*I_KINETIC_Iy5z_Fy2z_a;
  abcd[399] = 2.0E0*I_KINETIC_I5xy_F3z_a;
  abcd[400] = 2.0E0*I_KINETIC_I4x2y_F3z_a-1*I_KINETIC_G4x_F3z;
  abcd[401] = 2.0E0*I_KINETIC_I4xyz_F3z_a;
  abcd[402] = 2.0E0*I_KINETIC_I3x3y_F3z_a-2*I_KINETIC_G3xy_F3z;
  abcd[403] = 2.0E0*I_KINETIC_I3x2yz_F3z_a-1*I_KINETIC_G3xz_F3z;
  abcd[404] = 2.0E0*I_KINETIC_I3xy2z_F3z_a;
  abcd[405] = 2.0E0*I_KINETIC_I2x4y_F3z_a-3*I_KINETIC_G2x2y_F3z;
  abcd[406] = 2.0E0*I_KINETIC_I2x3yz_F3z_a-2*I_KINETIC_G2xyz_F3z;
  abcd[407] = 2.0E0*I_KINETIC_I2x2y2z_F3z_a-1*I_KINETIC_G2x2z_F3z;
  abcd[408] = 2.0E0*I_KINETIC_I2xy3z_F3z_a;
  abcd[409] = 2.0E0*I_KINETIC_Ix5y_F3z_a-4*I_KINETIC_Gx3y_F3z;
  abcd[410] = 2.0E0*I_KINETIC_Ix4yz_F3z_a-3*I_KINETIC_Gx2yz_F3z;
  abcd[411] = 2.0E0*I_KINETIC_Ix3y2z_F3z_a-2*I_KINETIC_Gxy2z_F3z;
  abcd[412] = 2.0E0*I_KINETIC_Ix2y3z_F3z_a-1*I_KINETIC_Gx3z_F3z;
  abcd[413] = 2.0E0*I_KINETIC_Ixy4z_F3z_a;
  abcd[414] = 2.0E0*I_KINETIC_I6y_F3z_a-5*I_KINETIC_G4y_F3z;
  abcd[415] = 2.0E0*I_KINETIC_I5yz_F3z_a-4*I_KINETIC_G3yz_F3z;
  abcd[416] = 2.0E0*I_KINETIC_I4y2z_F3z_a-3*I_KINETIC_G2y2z_F3z;
  abcd[417] = 2.0E0*I_KINETIC_I3y3z_F3z_a-2*I_KINETIC_Gy3z_F3z;
  abcd[418] = 2.0E0*I_KINETIC_I2y4z_F3z_a-1*I_KINETIC_G4z_F3z;
  abcd[419] = 2.0E0*I_KINETIC_Iy5z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_F_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_F_a
   * RHS shell quartet name: SQ_KINETIC_G_F
   ************************************************************/
  abcd[420] = 2.0E0*I_KINETIC_I5xz_F3x_a;
  abcd[421] = 2.0E0*I_KINETIC_I4xyz_F3x_a;
  abcd[422] = 2.0E0*I_KINETIC_I4x2z_F3x_a-1*I_KINETIC_G4x_F3x;
  abcd[423] = 2.0E0*I_KINETIC_I3x2yz_F3x_a;
  abcd[424] = 2.0E0*I_KINETIC_I3xy2z_F3x_a-1*I_KINETIC_G3xy_F3x;
  abcd[425] = 2.0E0*I_KINETIC_I3x3z_F3x_a-2*I_KINETIC_G3xz_F3x;
  abcd[426] = 2.0E0*I_KINETIC_I2x3yz_F3x_a;
  abcd[427] = 2.0E0*I_KINETIC_I2x2y2z_F3x_a-1*I_KINETIC_G2x2y_F3x;
  abcd[428] = 2.0E0*I_KINETIC_I2xy3z_F3x_a-2*I_KINETIC_G2xyz_F3x;
  abcd[429] = 2.0E0*I_KINETIC_I2x4z_F3x_a-3*I_KINETIC_G2x2z_F3x;
  abcd[430] = 2.0E0*I_KINETIC_Ix4yz_F3x_a;
  abcd[431] = 2.0E0*I_KINETIC_Ix3y2z_F3x_a-1*I_KINETIC_Gx3y_F3x;
  abcd[432] = 2.0E0*I_KINETIC_Ix2y3z_F3x_a-2*I_KINETIC_Gx2yz_F3x;
  abcd[433] = 2.0E0*I_KINETIC_Ixy4z_F3x_a-3*I_KINETIC_Gxy2z_F3x;
  abcd[434] = 2.0E0*I_KINETIC_Ix5z_F3x_a-4*I_KINETIC_Gx3z_F3x;
  abcd[435] = 2.0E0*I_KINETIC_I5yz_F3x_a;
  abcd[436] = 2.0E0*I_KINETIC_I4y2z_F3x_a-1*I_KINETIC_G4y_F3x;
  abcd[437] = 2.0E0*I_KINETIC_I3y3z_F3x_a-2*I_KINETIC_G3yz_F3x;
  abcd[438] = 2.0E0*I_KINETIC_I2y4z_F3x_a-3*I_KINETIC_G2y2z_F3x;
  abcd[439] = 2.0E0*I_KINETIC_Iy5z_F3x_a-4*I_KINETIC_Gy3z_F3x;
  abcd[440] = 2.0E0*I_KINETIC_I6z_F3x_a-5*I_KINETIC_G4z_F3x;
  abcd[441] = 2.0E0*I_KINETIC_I5xz_F2xy_a;
  abcd[442] = 2.0E0*I_KINETIC_I4xyz_F2xy_a;
  abcd[443] = 2.0E0*I_KINETIC_I4x2z_F2xy_a-1*I_KINETIC_G4x_F2xy;
  abcd[444] = 2.0E0*I_KINETIC_I3x2yz_F2xy_a;
  abcd[445] = 2.0E0*I_KINETIC_I3xy2z_F2xy_a-1*I_KINETIC_G3xy_F2xy;
  abcd[446] = 2.0E0*I_KINETIC_I3x3z_F2xy_a-2*I_KINETIC_G3xz_F2xy;
  abcd[447] = 2.0E0*I_KINETIC_I2x3yz_F2xy_a;
  abcd[448] = 2.0E0*I_KINETIC_I2x2y2z_F2xy_a-1*I_KINETIC_G2x2y_F2xy;
  abcd[449] = 2.0E0*I_KINETIC_I2xy3z_F2xy_a-2*I_KINETIC_G2xyz_F2xy;
  abcd[450] = 2.0E0*I_KINETIC_I2x4z_F2xy_a-3*I_KINETIC_G2x2z_F2xy;
  abcd[451] = 2.0E0*I_KINETIC_Ix4yz_F2xy_a;
  abcd[452] = 2.0E0*I_KINETIC_Ix3y2z_F2xy_a-1*I_KINETIC_Gx3y_F2xy;
  abcd[453] = 2.0E0*I_KINETIC_Ix2y3z_F2xy_a-2*I_KINETIC_Gx2yz_F2xy;
  abcd[454] = 2.0E0*I_KINETIC_Ixy4z_F2xy_a-3*I_KINETIC_Gxy2z_F2xy;
  abcd[455] = 2.0E0*I_KINETIC_Ix5z_F2xy_a-4*I_KINETIC_Gx3z_F2xy;
  abcd[456] = 2.0E0*I_KINETIC_I5yz_F2xy_a;
  abcd[457] = 2.0E0*I_KINETIC_I4y2z_F2xy_a-1*I_KINETIC_G4y_F2xy;
  abcd[458] = 2.0E0*I_KINETIC_I3y3z_F2xy_a-2*I_KINETIC_G3yz_F2xy;
  abcd[459] = 2.0E0*I_KINETIC_I2y4z_F2xy_a-3*I_KINETIC_G2y2z_F2xy;
  abcd[460] = 2.0E0*I_KINETIC_Iy5z_F2xy_a-4*I_KINETIC_Gy3z_F2xy;
  abcd[461] = 2.0E0*I_KINETIC_I6z_F2xy_a-5*I_KINETIC_G4z_F2xy;
  abcd[462] = 2.0E0*I_KINETIC_I5xz_F2xz_a;
  abcd[463] = 2.0E0*I_KINETIC_I4xyz_F2xz_a;
  abcd[464] = 2.0E0*I_KINETIC_I4x2z_F2xz_a-1*I_KINETIC_G4x_F2xz;
  abcd[465] = 2.0E0*I_KINETIC_I3x2yz_F2xz_a;
  abcd[466] = 2.0E0*I_KINETIC_I3xy2z_F2xz_a-1*I_KINETIC_G3xy_F2xz;
  abcd[467] = 2.0E0*I_KINETIC_I3x3z_F2xz_a-2*I_KINETIC_G3xz_F2xz;
  abcd[468] = 2.0E0*I_KINETIC_I2x3yz_F2xz_a;
  abcd[469] = 2.0E0*I_KINETIC_I2x2y2z_F2xz_a-1*I_KINETIC_G2x2y_F2xz;
  abcd[470] = 2.0E0*I_KINETIC_I2xy3z_F2xz_a-2*I_KINETIC_G2xyz_F2xz;
  abcd[471] = 2.0E0*I_KINETIC_I2x4z_F2xz_a-3*I_KINETIC_G2x2z_F2xz;
  abcd[472] = 2.0E0*I_KINETIC_Ix4yz_F2xz_a;
  abcd[473] = 2.0E0*I_KINETIC_Ix3y2z_F2xz_a-1*I_KINETIC_Gx3y_F2xz;
  abcd[474] = 2.0E0*I_KINETIC_Ix2y3z_F2xz_a-2*I_KINETIC_Gx2yz_F2xz;
  abcd[475] = 2.0E0*I_KINETIC_Ixy4z_F2xz_a-3*I_KINETIC_Gxy2z_F2xz;
  abcd[476] = 2.0E0*I_KINETIC_Ix5z_F2xz_a-4*I_KINETIC_Gx3z_F2xz;
  abcd[477] = 2.0E0*I_KINETIC_I5yz_F2xz_a;
  abcd[478] = 2.0E0*I_KINETIC_I4y2z_F2xz_a-1*I_KINETIC_G4y_F2xz;
  abcd[479] = 2.0E0*I_KINETIC_I3y3z_F2xz_a-2*I_KINETIC_G3yz_F2xz;
  abcd[480] = 2.0E0*I_KINETIC_I2y4z_F2xz_a-3*I_KINETIC_G2y2z_F2xz;
  abcd[481] = 2.0E0*I_KINETIC_Iy5z_F2xz_a-4*I_KINETIC_Gy3z_F2xz;
  abcd[482] = 2.0E0*I_KINETIC_I6z_F2xz_a-5*I_KINETIC_G4z_F2xz;
  abcd[483] = 2.0E0*I_KINETIC_I5xz_Fx2y_a;
  abcd[484] = 2.0E0*I_KINETIC_I4xyz_Fx2y_a;
  abcd[485] = 2.0E0*I_KINETIC_I4x2z_Fx2y_a-1*I_KINETIC_G4x_Fx2y;
  abcd[486] = 2.0E0*I_KINETIC_I3x2yz_Fx2y_a;
  abcd[487] = 2.0E0*I_KINETIC_I3xy2z_Fx2y_a-1*I_KINETIC_G3xy_Fx2y;
  abcd[488] = 2.0E0*I_KINETIC_I3x3z_Fx2y_a-2*I_KINETIC_G3xz_Fx2y;
  abcd[489] = 2.0E0*I_KINETIC_I2x3yz_Fx2y_a;
  abcd[490] = 2.0E0*I_KINETIC_I2x2y2z_Fx2y_a-1*I_KINETIC_G2x2y_Fx2y;
  abcd[491] = 2.0E0*I_KINETIC_I2xy3z_Fx2y_a-2*I_KINETIC_G2xyz_Fx2y;
  abcd[492] = 2.0E0*I_KINETIC_I2x4z_Fx2y_a-3*I_KINETIC_G2x2z_Fx2y;
  abcd[493] = 2.0E0*I_KINETIC_Ix4yz_Fx2y_a;
  abcd[494] = 2.0E0*I_KINETIC_Ix3y2z_Fx2y_a-1*I_KINETIC_Gx3y_Fx2y;
  abcd[495] = 2.0E0*I_KINETIC_Ix2y3z_Fx2y_a-2*I_KINETIC_Gx2yz_Fx2y;
  abcd[496] = 2.0E0*I_KINETIC_Ixy4z_Fx2y_a-3*I_KINETIC_Gxy2z_Fx2y;
  abcd[497] = 2.0E0*I_KINETIC_Ix5z_Fx2y_a-4*I_KINETIC_Gx3z_Fx2y;
  abcd[498] = 2.0E0*I_KINETIC_I5yz_Fx2y_a;
  abcd[499] = 2.0E0*I_KINETIC_I4y2z_Fx2y_a-1*I_KINETIC_G4y_Fx2y;
  abcd[500] = 2.0E0*I_KINETIC_I3y3z_Fx2y_a-2*I_KINETIC_G3yz_Fx2y;
  abcd[501] = 2.0E0*I_KINETIC_I2y4z_Fx2y_a-3*I_KINETIC_G2y2z_Fx2y;
  abcd[502] = 2.0E0*I_KINETIC_Iy5z_Fx2y_a-4*I_KINETIC_Gy3z_Fx2y;
  abcd[503] = 2.0E0*I_KINETIC_I6z_Fx2y_a-5*I_KINETIC_G4z_Fx2y;
  abcd[504] = 2.0E0*I_KINETIC_I5xz_Fxyz_a;
  abcd[505] = 2.0E0*I_KINETIC_I4xyz_Fxyz_a;
  abcd[506] = 2.0E0*I_KINETIC_I4x2z_Fxyz_a-1*I_KINETIC_G4x_Fxyz;
  abcd[507] = 2.0E0*I_KINETIC_I3x2yz_Fxyz_a;
  abcd[508] = 2.0E0*I_KINETIC_I3xy2z_Fxyz_a-1*I_KINETIC_G3xy_Fxyz;
  abcd[509] = 2.0E0*I_KINETIC_I3x3z_Fxyz_a-2*I_KINETIC_G3xz_Fxyz;
  abcd[510] = 2.0E0*I_KINETIC_I2x3yz_Fxyz_a;
  abcd[511] = 2.0E0*I_KINETIC_I2x2y2z_Fxyz_a-1*I_KINETIC_G2x2y_Fxyz;
  abcd[512] = 2.0E0*I_KINETIC_I2xy3z_Fxyz_a-2*I_KINETIC_G2xyz_Fxyz;
  abcd[513] = 2.0E0*I_KINETIC_I2x4z_Fxyz_a-3*I_KINETIC_G2x2z_Fxyz;
  abcd[514] = 2.0E0*I_KINETIC_Ix4yz_Fxyz_a;
  abcd[515] = 2.0E0*I_KINETIC_Ix3y2z_Fxyz_a-1*I_KINETIC_Gx3y_Fxyz;
  abcd[516] = 2.0E0*I_KINETIC_Ix2y3z_Fxyz_a-2*I_KINETIC_Gx2yz_Fxyz;
  abcd[517] = 2.0E0*I_KINETIC_Ixy4z_Fxyz_a-3*I_KINETIC_Gxy2z_Fxyz;
  abcd[518] = 2.0E0*I_KINETIC_Ix5z_Fxyz_a-4*I_KINETIC_Gx3z_Fxyz;
  abcd[519] = 2.0E0*I_KINETIC_I5yz_Fxyz_a;
  abcd[520] = 2.0E0*I_KINETIC_I4y2z_Fxyz_a-1*I_KINETIC_G4y_Fxyz;
  abcd[521] = 2.0E0*I_KINETIC_I3y3z_Fxyz_a-2*I_KINETIC_G3yz_Fxyz;
  abcd[522] = 2.0E0*I_KINETIC_I2y4z_Fxyz_a-3*I_KINETIC_G2y2z_Fxyz;
  abcd[523] = 2.0E0*I_KINETIC_Iy5z_Fxyz_a-4*I_KINETIC_Gy3z_Fxyz;
  abcd[524] = 2.0E0*I_KINETIC_I6z_Fxyz_a-5*I_KINETIC_G4z_Fxyz;
  abcd[525] = 2.0E0*I_KINETIC_I5xz_Fx2z_a;
  abcd[526] = 2.0E0*I_KINETIC_I4xyz_Fx2z_a;
  abcd[527] = 2.0E0*I_KINETIC_I4x2z_Fx2z_a-1*I_KINETIC_G4x_Fx2z;
  abcd[528] = 2.0E0*I_KINETIC_I3x2yz_Fx2z_a;
  abcd[529] = 2.0E0*I_KINETIC_I3xy2z_Fx2z_a-1*I_KINETIC_G3xy_Fx2z;
  abcd[530] = 2.0E0*I_KINETIC_I3x3z_Fx2z_a-2*I_KINETIC_G3xz_Fx2z;
  abcd[531] = 2.0E0*I_KINETIC_I2x3yz_Fx2z_a;
  abcd[532] = 2.0E0*I_KINETIC_I2x2y2z_Fx2z_a-1*I_KINETIC_G2x2y_Fx2z;
  abcd[533] = 2.0E0*I_KINETIC_I2xy3z_Fx2z_a-2*I_KINETIC_G2xyz_Fx2z;
  abcd[534] = 2.0E0*I_KINETIC_I2x4z_Fx2z_a-3*I_KINETIC_G2x2z_Fx2z;
  abcd[535] = 2.0E0*I_KINETIC_Ix4yz_Fx2z_a;
  abcd[536] = 2.0E0*I_KINETIC_Ix3y2z_Fx2z_a-1*I_KINETIC_Gx3y_Fx2z;
  abcd[537] = 2.0E0*I_KINETIC_Ix2y3z_Fx2z_a-2*I_KINETIC_Gx2yz_Fx2z;
  abcd[538] = 2.0E0*I_KINETIC_Ixy4z_Fx2z_a-3*I_KINETIC_Gxy2z_Fx2z;
  abcd[539] = 2.0E0*I_KINETIC_Ix5z_Fx2z_a-4*I_KINETIC_Gx3z_Fx2z;
  abcd[540] = 2.0E0*I_KINETIC_I5yz_Fx2z_a;
  abcd[541] = 2.0E0*I_KINETIC_I4y2z_Fx2z_a-1*I_KINETIC_G4y_Fx2z;
  abcd[542] = 2.0E0*I_KINETIC_I3y3z_Fx2z_a-2*I_KINETIC_G3yz_Fx2z;
  abcd[543] = 2.0E0*I_KINETIC_I2y4z_Fx2z_a-3*I_KINETIC_G2y2z_Fx2z;
  abcd[544] = 2.0E0*I_KINETIC_Iy5z_Fx2z_a-4*I_KINETIC_Gy3z_Fx2z;
  abcd[545] = 2.0E0*I_KINETIC_I6z_Fx2z_a-5*I_KINETIC_G4z_Fx2z;
  abcd[546] = 2.0E0*I_KINETIC_I5xz_F3y_a;
  abcd[547] = 2.0E0*I_KINETIC_I4xyz_F3y_a;
  abcd[548] = 2.0E0*I_KINETIC_I4x2z_F3y_a-1*I_KINETIC_G4x_F3y;
  abcd[549] = 2.0E0*I_KINETIC_I3x2yz_F3y_a;
  abcd[550] = 2.0E0*I_KINETIC_I3xy2z_F3y_a-1*I_KINETIC_G3xy_F3y;
  abcd[551] = 2.0E0*I_KINETIC_I3x3z_F3y_a-2*I_KINETIC_G3xz_F3y;
  abcd[552] = 2.0E0*I_KINETIC_I2x3yz_F3y_a;
  abcd[553] = 2.0E0*I_KINETIC_I2x2y2z_F3y_a-1*I_KINETIC_G2x2y_F3y;
  abcd[554] = 2.0E0*I_KINETIC_I2xy3z_F3y_a-2*I_KINETIC_G2xyz_F3y;
  abcd[555] = 2.0E0*I_KINETIC_I2x4z_F3y_a-3*I_KINETIC_G2x2z_F3y;
  abcd[556] = 2.0E0*I_KINETIC_Ix4yz_F3y_a;
  abcd[557] = 2.0E0*I_KINETIC_Ix3y2z_F3y_a-1*I_KINETIC_Gx3y_F3y;
  abcd[558] = 2.0E0*I_KINETIC_Ix2y3z_F3y_a-2*I_KINETIC_Gx2yz_F3y;
  abcd[559] = 2.0E0*I_KINETIC_Ixy4z_F3y_a-3*I_KINETIC_Gxy2z_F3y;
  abcd[560] = 2.0E0*I_KINETIC_Ix5z_F3y_a-4*I_KINETIC_Gx3z_F3y;
  abcd[561] = 2.0E0*I_KINETIC_I5yz_F3y_a;
  abcd[562] = 2.0E0*I_KINETIC_I4y2z_F3y_a-1*I_KINETIC_G4y_F3y;
  abcd[563] = 2.0E0*I_KINETIC_I3y3z_F3y_a-2*I_KINETIC_G3yz_F3y;
  abcd[564] = 2.0E0*I_KINETIC_I2y4z_F3y_a-3*I_KINETIC_G2y2z_F3y;
  abcd[565] = 2.0E0*I_KINETIC_Iy5z_F3y_a-4*I_KINETIC_Gy3z_F3y;
  abcd[566] = 2.0E0*I_KINETIC_I6z_F3y_a-5*I_KINETIC_G4z_F3y;
  abcd[567] = 2.0E0*I_KINETIC_I5xz_F2yz_a;
  abcd[568] = 2.0E0*I_KINETIC_I4xyz_F2yz_a;
  abcd[569] = 2.0E0*I_KINETIC_I4x2z_F2yz_a-1*I_KINETIC_G4x_F2yz;
  abcd[570] = 2.0E0*I_KINETIC_I3x2yz_F2yz_a;
  abcd[571] = 2.0E0*I_KINETIC_I3xy2z_F2yz_a-1*I_KINETIC_G3xy_F2yz;
  abcd[572] = 2.0E0*I_KINETIC_I3x3z_F2yz_a-2*I_KINETIC_G3xz_F2yz;
  abcd[573] = 2.0E0*I_KINETIC_I2x3yz_F2yz_a;
  abcd[574] = 2.0E0*I_KINETIC_I2x2y2z_F2yz_a-1*I_KINETIC_G2x2y_F2yz;
  abcd[575] = 2.0E0*I_KINETIC_I2xy3z_F2yz_a-2*I_KINETIC_G2xyz_F2yz;
  abcd[576] = 2.0E0*I_KINETIC_I2x4z_F2yz_a-3*I_KINETIC_G2x2z_F2yz;
  abcd[577] = 2.0E0*I_KINETIC_Ix4yz_F2yz_a;
  abcd[578] = 2.0E0*I_KINETIC_Ix3y2z_F2yz_a-1*I_KINETIC_Gx3y_F2yz;
  abcd[579] = 2.0E0*I_KINETIC_Ix2y3z_F2yz_a-2*I_KINETIC_Gx2yz_F2yz;
  abcd[580] = 2.0E0*I_KINETIC_Ixy4z_F2yz_a-3*I_KINETIC_Gxy2z_F2yz;
  abcd[581] = 2.0E0*I_KINETIC_Ix5z_F2yz_a-4*I_KINETIC_Gx3z_F2yz;
  abcd[582] = 2.0E0*I_KINETIC_I5yz_F2yz_a;
  abcd[583] = 2.0E0*I_KINETIC_I4y2z_F2yz_a-1*I_KINETIC_G4y_F2yz;
  abcd[584] = 2.0E0*I_KINETIC_I3y3z_F2yz_a-2*I_KINETIC_G3yz_F2yz;
  abcd[585] = 2.0E0*I_KINETIC_I2y4z_F2yz_a-3*I_KINETIC_G2y2z_F2yz;
  abcd[586] = 2.0E0*I_KINETIC_Iy5z_F2yz_a-4*I_KINETIC_Gy3z_F2yz;
  abcd[587] = 2.0E0*I_KINETIC_I6z_F2yz_a-5*I_KINETIC_G4z_F2yz;
  abcd[588] = 2.0E0*I_KINETIC_I5xz_Fy2z_a;
  abcd[589] = 2.0E0*I_KINETIC_I4xyz_Fy2z_a;
  abcd[590] = 2.0E0*I_KINETIC_I4x2z_Fy2z_a-1*I_KINETIC_G4x_Fy2z;
  abcd[591] = 2.0E0*I_KINETIC_I3x2yz_Fy2z_a;
  abcd[592] = 2.0E0*I_KINETIC_I3xy2z_Fy2z_a-1*I_KINETIC_G3xy_Fy2z;
  abcd[593] = 2.0E0*I_KINETIC_I3x3z_Fy2z_a-2*I_KINETIC_G3xz_Fy2z;
  abcd[594] = 2.0E0*I_KINETIC_I2x3yz_Fy2z_a;
  abcd[595] = 2.0E0*I_KINETIC_I2x2y2z_Fy2z_a-1*I_KINETIC_G2x2y_Fy2z;
  abcd[596] = 2.0E0*I_KINETIC_I2xy3z_Fy2z_a-2*I_KINETIC_G2xyz_Fy2z;
  abcd[597] = 2.0E0*I_KINETIC_I2x4z_Fy2z_a-3*I_KINETIC_G2x2z_Fy2z;
  abcd[598] = 2.0E0*I_KINETIC_Ix4yz_Fy2z_a;
  abcd[599] = 2.0E0*I_KINETIC_Ix3y2z_Fy2z_a-1*I_KINETIC_Gx3y_Fy2z;
  abcd[600] = 2.0E0*I_KINETIC_Ix2y3z_Fy2z_a-2*I_KINETIC_Gx2yz_Fy2z;
  abcd[601] = 2.0E0*I_KINETIC_Ixy4z_Fy2z_a-3*I_KINETIC_Gxy2z_Fy2z;
  abcd[602] = 2.0E0*I_KINETIC_Ix5z_Fy2z_a-4*I_KINETIC_Gx3z_Fy2z;
  abcd[603] = 2.0E0*I_KINETIC_I5yz_Fy2z_a;
  abcd[604] = 2.0E0*I_KINETIC_I4y2z_Fy2z_a-1*I_KINETIC_G4y_Fy2z;
  abcd[605] = 2.0E0*I_KINETIC_I3y3z_Fy2z_a-2*I_KINETIC_G3yz_Fy2z;
  abcd[606] = 2.0E0*I_KINETIC_I2y4z_Fy2z_a-3*I_KINETIC_G2y2z_Fy2z;
  abcd[607] = 2.0E0*I_KINETIC_Iy5z_Fy2z_a-4*I_KINETIC_Gy3z_Fy2z;
  abcd[608] = 2.0E0*I_KINETIC_I6z_Fy2z_a-5*I_KINETIC_G4z_Fy2z;
  abcd[609] = 2.0E0*I_KINETIC_I5xz_F3z_a;
  abcd[610] = 2.0E0*I_KINETIC_I4xyz_F3z_a;
  abcd[611] = 2.0E0*I_KINETIC_I4x2z_F3z_a-1*I_KINETIC_G4x_F3z;
  abcd[612] = 2.0E0*I_KINETIC_I3x2yz_F3z_a;
  abcd[613] = 2.0E0*I_KINETIC_I3xy2z_F3z_a-1*I_KINETIC_G3xy_F3z;
  abcd[614] = 2.0E0*I_KINETIC_I3x3z_F3z_a-2*I_KINETIC_G3xz_F3z;
  abcd[615] = 2.0E0*I_KINETIC_I2x3yz_F3z_a;
  abcd[616] = 2.0E0*I_KINETIC_I2x2y2z_F3z_a-1*I_KINETIC_G2x2y_F3z;
  abcd[617] = 2.0E0*I_KINETIC_I2xy3z_F3z_a-2*I_KINETIC_G2xyz_F3z;
  abcd[618] = 2.0E0*I_KINETIC_I2x4z_F3z_a-3*I_KINETIC_G2x2z_F3z;
  abcd[619] = 2.0E0*I_KINETIC_Ix4yz_F3z_a;
  abcd[620] = 2.0E0*I_KINETIC_Ix3y2z_F3z_a-1*I_KINETIC_Gx3y_F3z;
  abcd[621] = 2.0E0*I_KINETIC_Ix2y3z_F3z_a-2*I_KINETIC_Gx2yz_F3z;
  abcd[622] = 2.0E0*I_KINETIC_Ixy4z_F3z_a-3*I_KINETIC_Gxy2z_F3z;
  abcd[623] = 2.0E0*I_KINETIC_Ix5z_F3z_a-4*I_KINETIC_Gx3z_F3z;
  abcd[624] = 2.0E0*I_KINETIC_I5yz_F3z_a;
  abcd[625] = 2.0E0*I_KINETIC_I4y2z_F3z_a-1*I_KINETIC_G4y_F3z;
  abcd[626] = 2.0E0*I_KINETIC_I3y3z_F3z_a-2*I_KINETIC_G3yz_F3z;
  abcd[627] = 2.0E0*I_KINETIC_I2y4z_F3z_a-3*I_KINETIC_G2y2z_F3z;
  abcd[628] = 2.0E0*I_KINETIC_Iy5z_F3z_a-4*I_KINETIC_Gy3z_F3z;
  abcd[629] = 2.0E0*I_KINETIC_I6z_F3z_a-5*I_KINETIC_G4z_F3z;
}
