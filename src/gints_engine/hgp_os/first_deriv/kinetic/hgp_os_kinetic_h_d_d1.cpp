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
// BRA1 as redundant position, total RHS integrals evaluated as: 3473
// BRA2 as redundant position, total RHS integrals evaluated as: 3422
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

void hgp_os_kinetic_h_d_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_I6x_D2x_a = 0.0E0;
  Double I_KINETIC_I5xy_D2x_a = 0.0E0;
  Double I_KINETIC_I5xz_D2x_a = 0.0E0;
  Double I_KINETIC_I4x2y_D2x_a = 0.0E0;
  Double I_KINETIC_I4xyz_D2x_a = 0.0E0;
  Double I_KINETIC_I4x2z_D2x_a = 0.0E0;
  Double I_KINETIC_I3x3y_D2x_a = 0.0E0;
  Double I_KINETIC_I3x2yz_D2x_a = 0.0E0;
  Double I_KINETIC_I3xy2z_D2x_a = 0.0E0;
  Double I_KINETIC_I3x3z_D2x_a = 0.0E0;
  Double I_KINETIC_I2x4y_D2x_a = 0.0E0;
  Double I_KINETIC_I2x3yz_D2x_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_D2x_a = 0.0E0;
  Double I_KINETIC_I2xy3z_D2x_a = 0.0E0;
  Double I_KINETIC_I2x4z_D2x_a = 0.0E0;
  Double I_KINETIC_Ix5y_D2x_a = 0.0E0;
  Double I_KINETIC_Ix4yz_D2x_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_D2x_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_D2x_a = 0.0E0;
  Double I_KINETIC_Ixy4z_D2x_a = 0.0E0;
  Double I_KINETIC_Ix5z_D2x_a = 0.0E0;
  Double I_KINETIC_I6y_D2x_a = 0.0E0;
  Double I_KINETIC_I5yz_D2x_a = 0.0E0;
  Double I_KINETIC_I4y2z_D2x_a = 0.0E0;
  Double I_KINETIC_I3y3z_D2x_a = 0.0E0;
  Double I_KINETIC_I2y4z_D2x_a = 0.0E0;
  Double I_KINETIC_Iy5z_D2x_a = 0.0E0;
  Double I_KINETIC_I6z_D2x_a = 0.0E0;
  Double I_KINETIC_I6x_Dxy_a = 0.0E0;
  Double I_KINETIC_I5xy_Dxy_a = 0.0E0;
  Double I_KINETIC_I5xz_Dxy_a = 0.0E0;
  Double I_KINETIC_I4x2y_Dxy_a = 0.0E0;
  Double I_KINETIC_I4xyz_Dxy_a = 0.0E0;
  Double I_KINETIC_I4x2z_Dxy_a = 0.0E0;
  Double I_KINETIC_I3x3y_Dxy_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Dxy_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Dxy_a = 0.0E0;
  Double I_KINETIC_I3x3z_Dxy_a = 0.0E0;
  Double I_KINETIC_I2x4y_Dxy_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Dxy_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Dxy_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Dxy_a = 0.0E0;
  Double I_KINETIC_I2x4z_Dxy_a = 0.0E0;
  Double I_KINETIC_Ix5y_Dxy_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Dxy_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Dxy_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Dxy_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Dxy_a = 0.0E0;
  Double I_KINETIC_Ix5z_Dxy_a = 0.0E0;
  Double I_KINETIC_I6y_Dxy_a = 0.0E0;
  Double I_KINETIC_I5yz_Dxy_a = 0.0E0;
  Double I_KINETIC_I4y2z_Dxy_a = 0.0E0;
  Double I_KINETIC_I3y3z_Dxy_a = 0.0E0;
  Double I_KINETIC_I2y4z_Dxy_a = 0.0E0;
  Double I_KINETIC_Iy5z_Dxy_a = 0.0E0;
  Double I_KINETIC_I6z_Dxy_a = 0.0E0;
  Double I_KINETIC_I6x_Dxz_a = 0.0E0;
  Double I_KINETIC_I5xy_Dxz_a = 0.0E0;
  Double I_KINETIC_I5xz_Dxz_a = 0.0E0;
  Double I_KINETIC_I4x2y_Dxz_a = 0.0E0;
  Double I_KINETIC_I4xyz_Dxz_a = 0.0E0;
  Double I_KINETIC_I4x2z_Dxz_a = 0.0E0;
  Double I_KINETIC_I3x3y_Dxz_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Dxz_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Dxz_a = 0.0E0;
  Double I_KINETIC_I3x3z_Dxz_a = 0.0E0;
  Double I_KINETIC_I2x4y_Dxz_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Dxz_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Dxz_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Dxz_a = 0.0E0;
  Double I_KINETIC_I2x4z_Dxz_a = 0.0E0;
  Double I_KINETIC_Ix5y_Dxz_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Dxz_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Dxz_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Dxz_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Dxz_a = 0.0E0;
  Double I_KINETIC_Ix5z_Dxz_a = 0.0E0;
  Double I_KINETIC_I6y_Dxz_a = 0.0E0;
  Double I_KINETIC_I5yz_Dxz_a = 0.0E0;
  Double I_KINETIC_I4y2z_Dxz_a = 0.0E0;
  Double I_KINETIC_I3y3z_Dxz_a = 0.0E0;
  Double I_KINETIC_I2y4z_Dxz_a = 0.0E0;
  Double I_KINETIC_Iy5z_Dxz_a = 0.0E0;
  Double I_KINETIC_I6z_Dxz_a = 0.0E0;
  Double I_KINETIC_I6x_D2y_a = 0.0E0;
  Double I_KINETIC_I5xy_D2y_a = 0.0E0;
  Double I_KINETIC_I5xz_D2y_a = 0.0E0;
  Double I_KINETIC_I4x2y_D2y_a = 0.0E0;
  Double I_KINETIC_I4xyz_D2y_a = 0.0E0;
  Double I_KINETIC_I4x2z_D2y_a = 0.0E0;
  Double I_KINETIC_I3x3y_D2y_a = 0.0E0;
  Double I_KINETIC_I3x2yz_D2y_a = 0.0E0;
  Double I_KINETIC_I3xy2z_D2y_a = 0.0E0;
  Double I_KINETIC_I3x3z_D2y_a = 0.0E0;
  Double I_KINETIC_I2x4y_D2y_a = 0.0E0;
  Double I_KINETIC_I2x3yz_D2y_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_D2y_a = 0.0E0;
  Double I_KINETIC_I2xy3z_D2y_a = 0.0E0;
  Double I_KINETIC_I2x4z_D2y_a = 0.0E0;
  Double I_KINETIC_Ix5y_D2y_a = 0.0E0;
  Double I_KINETIC_Ix4yz_D2y_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_D2y_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_D2y_a = 0.0E0;
  Double I_KINETIC_Ixy4z_D2y_a = 0.0E0;
  Double I_KINETIC_Ix5z_D2y_a = 0.0E0;
  Double I_KINETIC_I6y_D2y_a = 0.0E0;
  Double I_KINETIC_I5yz_D2y_a = 0.0E0;
  Double I_KINETIC_I4y2z_D2y_a = 0.0E0;
  Double I_KINETIC_I3y3z_D2y_a = 0.0E0;
  Double I_KINETIC_I2y4z_D2y_a = 0.0E0;
  Double I_KINETIC_Iy5z_D2y_a = 0.0E0;
  Double I_KINETIC_I6z_D2y_a = 0.0E0;
  Double I_KINETIC_I6x_Dyz_a = 0.0E0;
  Double I_KINETIC_I5xy_Dyz_a = 0.0E0;
  Double I_KINETIC_I5xz_Dyz_a = 0.0E0;
  Double I_KINETIC_I4x2y_Dyz_a = 0.0E0;
  Double I_KINETIC_I4xyz_Dyz_a = 0.0E0;
  Double I_KINETIC_I4x2z_Dyz_a = 0.0E0;
  Double I_KINETIC_I3x3y_Dyz_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Dyz_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Dyz_a = 0.0E0;
  Double I_KINETIC_I3x3z_Dyz_a = 0.0E0;
  Double I_KINETIC_I2x4y_Dyz_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Dyz_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Dyz_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Dyz_a = 0.0E0;
  Double I_KINETIC_I2x4z_Dyz_a = 0.0E0;
  Double I_KINETIC_Ix5y_Dyz_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Dyz_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Dyz_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Dyz_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Dyz_a = 0.0E0;
  Double I_KINETIC_Ix5z_Dyz_a = 0.0E0;
  Double I_KINETIC_I6y_Dyz_a = 0.0E0;
  Double I_KINETIC_I5yz_Dyz_a = 0.0E0;
  Double I_KINETIC_I4y2z_Dyz_a = 0.0E0;
  Double I_KINETIC_I3y3z_Dyz_a = 0.0E0;
  Double I_KINETIC_I2y4z_Dyz_a = 0.0E0;
  Double I_KINETIC_Iy5z_Dyz_a = 0.0E0;
  Double I_KINETIC_I6z_Dyz_a = 0.0E0;
  Double I_KINETIC_I6x_D2z_a = 0.0E0;
  Double I_KINETIC_I5xy_D2z_a = 0.0E0;
  Double I_KINETIC_I5xz_D2z_a = 0.0E0;
  Double I_KINETIC_I4x2y_D2z_a = 0.0E0;
  Double I_KINETIC_I4xyz_D2z_a = 0.0E0;
  Double I_KINETIC_I4x2z_D2z_a = 0.0E0;
  Double I_KINETIC_I3x3y_D2z_a = 0.0E0;
  Double I_KINETIC_I3x2yz_D2z_a = 0.0E0;
  Double I_KINETIC_I3xy2z_D2z_a = 0.0E0;
  Double I_KINETIC_I3x3z_D2z_a = 0.0E0;
  Double I_KINETIC_I2x4y_D2z_a = 0.0E0;
  Double I_KINETIC_I2x3yz_D2z_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_D2z_a = 0.0E0;
  Double I_KINETIC_I2xy3z_D2z_a = 0.0E0;
  Double I_KINETIC_I2x4z_D2z_a = 0.0E0;
  Double I_KINETIC_Ix5y_D2z_a = 0.0E0;
  Double I_KINETIC_Ix4yz_D2z_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_D2z_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_D2z_a = 0.0E0;
  Double I_KINETIC_Ixy4z_D2z_a = 0.0E0;
  Double I_KINETIC_Ix5z_D2z_a = 0.0E0;
  Double I_KINETIC_I6y_D2z_a = 0.0E0;
  Double I_KINETIC_I5yz_D2z_a = 0.0E0;
  Double I_KINETIC_I4y2z_D2z_a = 0.0E0;
  Double I_KINETIC_I3y3z_D2z_a = 0.0E0;
  Double I_KINETIC_I2y4z_D2z_a = 0.0E0;
  Double I_KINETIC_Iy5z_D2z_a = 0.0E0;
  Double I_KINETIC_I6z_D2z_a = 0.0E0;
  Double I_KINETIC_G4x_D2x = 0.0E0;
  Double I_KINETIC_G3xy_D2x = 0.0E0;
  Double I_KINETIC_G3xz_D2x = 0.0E0;
  Double I_KINETIC_G2x2y_D2x = 0.0E0;
  Double I_KINETIC_G2xyz_D2x = 0.0E0;
  Double I_KINETIC_G2x2z_D2x = 0.0E0;
  Double I_KINETIC_Gx3y_D2x = 0.0E0;
  Double I_KINETIC_Gx2yz_D2x = 0.0E0;
  Double I_KINETIC_Gxy2z_D2x = 0.0E0;
  Double I_KINETIC_Gx3z_D2x = 0.0E0;
  Double I_KINETIC_G4y_D2x = 0.0E0;
  Double I_KINETIC_G3yz_D2x = 0.0E0;
  Double I_KINETIC_G2y2z_D2x = 0.0E0;
  Double I_KINETIC_Gy3z_D2x = 0.0E0;
  Double I_KINETIC_G4z_D2x = 0.0E0;
  Double I_KINETIC_G4x_Dxy = 0.0E0;
  Double I_KINETIC_G3xy_Dxy = 0.0E0;
  Double I_KINETIC_G3xz_Dxy = 0.0E0;
  Double I_KINETIC_G2x2y_Dxy = 0.0E0;
  Double I_KINETIC_G2xyz_Dxy = 0.0E0;
  Double I_KINETIC_G2x2z_Dxy = 0.0E0;
  Double I_KINETIC_Gx3y_Dxy = 0.0E0;
  Double I_KINETIC_Gx2yz_Dxy = 0.0E0;
  Double I_KINETIC_Gxy2z_Dxy = 0.0E0;
  Double I_KINETIC_Gx3z_Dxy = 0.0E0;
  Double I_KINETIC_G4y_Dxy = 0.0E0;
  Double I_KINETIC_G3yz_Dxy = 0.0E0;
  Double I_KINETIC_G2y2z_Dxy = 0.0E0;
  Double I_KINETIC_Gy3z_Dxy = 0.0E0;
  Double I_KINETIC_G4z_Dxy = 0.0E0;
  Double I_KINETIC_G4x_Dxz = 0.0E0;
  Double I_KINETIC_G3xy_Dxz = 0.0E0;
  Double I_KINETIC_G3xz_Dxz = 0.0E0;
  Double I_KINETIC_G2x2y_Dxz = 0.0E0;
  Double I_KINETIC_G2xyz_Dxz = 0.0E0;
  Double I_KINETIC_G2x2z_Dxz = 0.0E0;
  Double I_KINETIC_Gx3y_Dxz = 0.0E0;
  Double I_KINETIC_Gx2yz_Dxz = 0.0E0;
  Double I_KINETIC_Gxy2z_Dxz = 0.0E0;
  Double I_KINETIC_Gx3z_Dxz = 0.0E0;
  Double I_KINETIC_G4y_Dxz = 0.0E0;
  Double I_KINETIC_G3yz_Dxz = 0.0E0;
  Double I_KINETIC_G2y2z_Dxz = 0.0E0;
  Double I_KINETIC_Gy3z_Dxz = 0.0E0;
  Double I_KINETIC_G4z_Dxz = 0.0E0;
  Double I_KINETIC_G4x_D2y = 0.0E0;
  Double I_KINETIC_G3xy_D2y = 0.0E0;
  Double I_KINETIC_G3xz_D2y = 0.0E0;
  Double I_KINETIC_G2x2y_D2y = 0.0E0;
  Double I_KINETIC_G2xyz_D2y = 0.0E0;
  Double I_KINETIC_G2x2z_D2y = 0.0E0;
  Double I_KINETIC_Gx3y_D2y = 0.0E0;
  Double I_KINETIC_Gx2yz_D2y = 0.0E0;
  Double I_KINETIC_Gxy2z_D2y = 0.0E0;
  Double I_KINETIC_Gx3z_D2y = 0.0E0;
  Double I_KINETIC_G4y_D2y = 0.0E0;
  Double I_KINETIC_G3yz_D2y = 0.0E0;
  Double I_KINETIC_G2y2z_D2y = 0.0E0;
  Double I_KINETIC_Gy3z_D2y = 0.0E0;
  Double I_KINETIC_G4z_D2y = 0.0E0;
  Double I_KINETIC_G4x_Dyz = 0.0E0;
  Double I_KINETIC_G3xy_Dyz = 0.0E0;
  Double I_KINETIC_G3xz_Dyz = 0.0E0;
  Double I_KINETIC_G2x2y_Dyz = 0.0E0;
  Double I_KINETIC_G2xyz_Dyz = 0.0E0;
  Double I_KINETIC_G2x2z_Dyz = 0.0E0;
  Double I_KINETIC_Gx3y_Dyz = 0.0E0;
  Double I_KINETIC_Gx2yz_Dyz = 0.0E0;
  Double I_KINETIC_Gxy2z_Dyz = 0.0E0;
  Double I_KINETIC_Gx3z_Dyz = 0.0E0;
  Double I_KINETIC_G4y_Dyz = 0.0E0;
  Double I_KINETIC_G3yz_Dyz = 0.0E0;
  Double I_KINETIC_G2y2z_Dyz = 0.0E0;
  Double I_KINETIC_Gy3z_Dyz = 0.0E0;
  Double I_KINETIC_G4z_Dyz = 0.0E0;
  Double I_KINETIC_G4x_D2z = 0.0E0;
  Double I_KINETIC_G3xy_D2z = 0.0E0;
  Double I_KINETIC_G3xz_D2z = 0.0E0;
  Double I_KINETIC_G2x2y_D2z = 0.0E0;
  Double I_KINETIC_G2xyz_D2z = 0.0E0;
  Double I_KINETIC_G2x2z_D2z = 0.0E0;
  Double I_KINETIC_Gx3y_D2z = 0.0E0;
  Double I_KINETIC_Gx2yz_D2z = 0.0E0;
  Double I_KINETIC_Gxy2z_D2z = 0.0E0;
  Double I_KINETIC_Gx3z_D2z = 0.0E0;
  Double I_KINETIC_G4y_D2z = 0.0E0;
  Double I_KINETIC_G3yz_D2z = 0.0E0;
  Double I_KINETIC_G2y2z_D2z = 0.0E0;
  Double I_KINETIC_Gy3z_D2z = 0.0E0;
  Double I_KINETIC_G4z_D2z = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_Px_vrr = PBX*I_TWOBODYOVERLAP_G4x_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Px_vrr = PBX*I_TWOBODYOVERLAP_G3xy_S_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Px_vrr = PBX*I_TWOBODYOVERLAP_G3xz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Px_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Px_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Px_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Px_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Px_vrr = PBX*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Px_vrr = PBX*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Px_vrr = PBX*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Py_vrr = PBY*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Py_vrr = PBY*I_TWOBODYOVERLAP_G3xy_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Py_vrr = PBY*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Py_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Py_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Py_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Py_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Py_vrr = PBY*I_TWOBODYOVERLAP_G4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Py_vrr = PBY*I_TWOBODYOVERLAP_G3yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Py_vrr = PBY*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 5 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_S_vrr = PAX*I_TWOBODYOVERLAP_G4x_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_S_vrr = PAY*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_S_vrr = PAZ*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_S_vrr = PAY*I_TWOBODYOVERLAP_G3xy_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_S_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_S_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_S_vrr = PAX*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_S_vrr = PAX*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_S_vrr = PAY*I_TWOBODYOVERLAP_G4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_S_vrr = PAY*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_S_vrr = PAZ*I_TWOBODYOVERLAP_G4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
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
    Double I_TWOBODYOVERLAP_G2xyz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Px_vrr;
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
    Double I_TWOBODYOVERLAP_G2xyz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Py_vrr;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_D2x_vrr = PBX*I_TWOBODYOVERLAP_H5x_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2x_vrr = PBX*I_TWOBODYOVERLAP_H4xy_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2x_vrr = PBX*I_TWOBODYOVERLAP_H4xz_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2x_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2x_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2x_vrr = PBX*I_TWOBODYOVERLAP_H5y_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_H4yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2x_vrr = PBX*I_TWOBODYOVERLAP_H5z_Px_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H4xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H5y_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H4yz_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_Px_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_D2y_vrr = PBY*I_TWOBODYOVERLAP_H5x_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2y_vrr = PBY*I_TWOBODYOVERLAP_H4xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2y_vrr = PBY*I_TWOBODYOVERLAP_H4xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2y_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2y_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2y_vrr = PBY*I_TWOBODYOVERLAP_H5y_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_H4yz_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2y_vrr = PBY*I_TWOBODYOVERLAP_H5z_Py_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_Py_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H5x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H5y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_H5z_Pz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_D2x_vrr = PAX*I_TWOBODYOVERLAP_H5x_D2x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_I5xy_D2x_vrr = PAY*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_I5xz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_D2x_vrr = PAY*I_TWOBODYOVERLAP_H4xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_D2x_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_D2x_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_D2x_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_D2x_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_D2x_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_D2x_vrr = PAX*I_TWOBODYOVERLAP_H5y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_D2x_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_D2x_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_D2x_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_D2x_vrr = PAX*I_TWOBODYOVERLAP_H5z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6y_D2x_vrr = PAY*I_TWOBODYOVERLAP_H5y_D2x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I5yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_D2x_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_D2x_vrr = PAY*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I6z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_H5z_D2x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_I6x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_H5x_Dxy_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_H5x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_H5y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_H5z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_TWOBODYOVERLAP_I6y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_H5y_Dxy_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_H5z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Dxy_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_I6x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_H5x_Dxz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_H5y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_H5z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Pz_vrr;
    Double I_TWOBODYOVERLAP_I6y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_H5y_Dxz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_I6z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Dxz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_I6x_D2y_vrr = PAX*I_TWOBODYOVERLAP_H5x_D2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_I5xy_D2y_vrr = PAY*I_TWOBODYOVERLAP_H5x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_TWOBODYOVERLAP_I5xz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_D2y_vrr = PAY*I_TWOBODYOVERLAP_H4xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_D2y_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_D2y_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_D2y_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_D2y_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3z_Py_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_D2y_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_D2y_vrr = PAX*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_D2y_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_D2y_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_D2y_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_D2y_vrr = PAX*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_I6y_D2y_vrr = PAY*I_TWOBODYOVERLAP_H5y_D2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_I5yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_D2y_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_D2y_vrr = PAY*I_TWOBODYOVERLAP_H5z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_TWOBODYOVERLAP_I6z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_H5z_D2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_I6x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_H5x_Dyz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_H5x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H5x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_H4xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_I6y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_H5y_Dyz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H5y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_H5z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Pz_vrr;
    Double I_TWOBODYOVERLAP_I6z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_H5z_Dyz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_TWOBODYOVERLAP_I6x_D2z_vrr = PAX*I_TWOBODYOVERLAP_H5x_D2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_I5xy_D2z_vrr = PAY*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_TWOBODYOVERLAP_I5xz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H5x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_D2z_vrr = PAY*I_TWOBODYOVERLAP_H4xy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_D2z_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_D2z_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_D2z_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_D2z_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_D2z_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_D2z_vrr = PAX*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_D2z_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_D2z_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_D2z_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_D2z_vrr = PAX*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I6y_D2z_vrr = PAY*I_TWOBODYOVERLAP_H5y_D2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_I5yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H5y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_D2z_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_D2z_vrr = PAY*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_I6z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_H5z_D2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Pz_vrr;

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
     * shell quartet name: SQ_KINETIC_G_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_KINETIC_G4x_Px_vrr = PBX*I_KINETIC_G4x_S_vrr+4*oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_KINETIC_G3xy_Px_vrr = PBX*I_KINETIC_G3xy_S_vrr+3*oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_KINETIC_G3xz_Px_vrr = PBX*I_KINETIC_G3xz_S_vrr+3*oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_KINETIC_G2x2y_Px_vrr = PBX*I_KINETIC_G2x2y_S_vrr+2*oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_KINETIC_G2xyz_Px_vrr = PBX*I_KINETIC_G2xyz_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_KINETIC_G2x2z_Px_vrr = PBX*I_KINETIC_G2x2z_S_vrr+2*oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_KINETIC_Gx3y_Px_vrr = PBX*I_KINETIC_Gx3y_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_KINETIC_Gx2yz_Px_vrr = PBX*I_KINETIC_Gx2yz_S_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_KINETIC_Gxy2z_Px_vrr = PBX*I_KINETIC_Gxy2z_S_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_KINETIC_Gx3z_Px_vrr = PBX*I_KINETIC_Gx3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_KINETIC_G4y_Px_vrr = PBX*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_KINETIC_G3yz_Px_vrr = PBX*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_KINETIC_G2y2z_Px_vrr = PBX*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_KINETIC_Gy3z_Px_vrr = PBX*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_KINETIC_G4z_Px_vrr = PBX*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_KINETIC_G4x_Py_vrr = PBY*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_KINETIC_G3xy_Py_vrr = PBY*I_KINETIC_G3xy_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_KINETIC_G3xz_Py_vrr = PBY*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_KINETIC_G2x2y_Py_vrr = PBY*I_KINETIC_G2x2y_S_vrr+2*oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_KINETIC_G2xyz_Py_vrr = PBY*I_KINETIC_G2xyz_S_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_KINETIC_G2x2z_Py_vrr = PBY*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_KINETIC_Gx3y_Py_vrr = PBY*I_KINETIC_Gx3y_S_vrr+3*oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_KINETIC_Gx2yz_Py_vrr = PBY*I_KINETIC_Gx2yz_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_KINETIC_Gxy2z_Py_vrr = PBY*I_KINETIC_Gxy2z_S_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_KINETIC_Gx3z_Py_vrr = PBY*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_KINETIC_G4y_Py_vrr = PBY*I_KINETIC_G4y_S_vrr+4*oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_KINETIC_G3yz_Py_vrr = PBY*I_KINETIC_G3yz_S_vrr+3*oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_KINETIC_G2y2z_Py_vrr = PBY*I_KINETIC_G2y2z_S_vrr+2*oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_KINETIC_Gy3z_Py_vrr = PBY*I_KINETIC_Gy3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_KINETIC_G4z_Py_vrr = PBY*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_KINETIC_G4x_Pz_vrr = PBZ*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_KINETIC_G3xy_Pz_vrr = PBZ*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_KINETIC_G3xz_Pz_vrr = PBZ*I_KINETIC_G3xz_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_KINETIC_G2x2y_Pz_vrr = PBZ*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_KINETIC_G2xyz_Pz_vrr = PBZ*I_KINETIC_G2xyz_S_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Pz_vrr;
    Double I_KINETIC_G2x2z_Pz_vrr = PBZ*I_KINETIC_G2x2z_S_vrr+2*oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_KINETIC_Gx3y_Pz_vrr = PBZ*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_KINETIC_Gx2yz_Pz_vrr = PBZ*I_KINETIC_Gx2yz_S_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr;
    Double I_KINETIC_Gxy2z_Pz_vrr = PBZ*I_KINETIC_Gxy2z_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr;
    Double I_KINETIC_Gx3z_Pz_vrr = PBZ*I_KINETIC_Gx3z_S_vrr+3*oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_KINETIC_G4y_Pz_vrr = PBZ*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_KINETIC_G3yz_Pz_vrr = PBZ*I_KINETIC_G3yz_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_KINETIC_G2y2z_Pz_vrr = PBZ*I_KINETIC_G2y2z_S_vrr+2*oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_KINETIC_Gy3z_Pz_vrr = PBZ*I_KINETIC_Gy3z_S_vrr+3*oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_KINETIC_G4z_Pz_vrr = PBZ*I_KINETIC_G4z_S_vrr+4*oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 5 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_KINETIC_H5x_S_vrr = PAX*I_KINETIC_G4x_S_vrr+4*oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H4xy_S_vrr = PAY*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_S_vrr = PAZ*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_S_vrr = PAY*I_KINETIC_G3xy_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_S_vrr-bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H3x2z_S_vrr = PAZ*I_KINETIC_G3xz_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_S_vrr-bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H2x3y_S_vrr = PAX*I_KINETIC_Gx3y_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_S_vrr-bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H2x2yz_S_vrr = PAZ*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2x3z_S_vrr = PAX*I_KINETIC_Gx3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_S_vrr-bdz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_Hx4y_S_vrr = PAX*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx4z_S_vrr = PAX*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_S_vrr = PAY*I_KINETIC_G4y_S_vrr+4*oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H4yz_S_vrr = PAZ*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_S_vrr = PAZ*I_KINETIC_G3yz_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_S_vrr-bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H2y3z_S_vrr = PAY*I_KINETIC_Gy3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_S_vrr-bdz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_Hy4z_S_vrr = PAY*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_S_vrr = PAZ*I_KINETIC_G4z_S_vrr+4*oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
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
    Double I_KINETIC_G2xyz_Dxz_vrr = PBZ*I_KINETIC_G2xyz_Px_vrr+oned2z*I_KINETIC_F2xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Dxz_vrr;
    Double I_KINETIC_G2x2z_Dxz_vrr = PBZ*I_KINETIC_G2x2z_Px_vrr+2*oned2z*I_KINETIC_F2xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dxz_vrr;
    Double I_KINETIC_Gx3y_Dxz_vrr = PBZ*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dxz_vrr;
    Double I_KINETIC_Gx2yz_Dxz_vrr = PBZ*I_KINETIC_Gx2yz_Px_vrr+oned2z*I_KINETIC_Fx2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Dxz_vrr;
    Double I_KINETIC_Gxy2z_Dxz_vrr = PBZ*I_KINETIC_Gxy2z_Px_vrr+2*oned2z*I_KINETIC_Fxyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Dxz_vrr;
    Double I_KINETIC_Gx3z_Dxz_vrr = PBZ*I_KINETIC_Gx3z_Px_vrr+3*oned2z*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dxz_vrr;
    Double I_KINETIC_G4y_Dxz_vrr = PBZ*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_KINETIC_G3yz_Dxz_vrr = PBZ*I_KINETIC_G3yz_Px_vrr+oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dxz_vrr;
    Double I_KINETIC_G2y2z_Dxz_vrr = PBZ*I_KINETIC_G2y2z_Px_vrr+2*oned2z*I_KINETIC_F2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dxz_vrr;
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
    Double I_KINETIC_G2xyz_Dyz_vrr = PBZ*I_KINETIC_G2xyz_Py_vrr+oned2z*I_KINETIC_F2xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Dyz_vrr;
    Double I_KINETIC_G2x2z_Dyz_vrr = PBZ*I_KINETIC_G2x2z_Py_vrr+2*oned2z*I_KINETIC_F2xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dyz_vrr;
    Double I_KINETIC_Gx3y_Dyz_vrr = PBZ*I_KINETIC_Gx3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dyz_vrr;
    Double I_KINETIC_Gx2yz_Dyz_vrr = PBZ*I_KINETIC_Gx2yz_Py_vrr+oned2z*I_KINETIC_Fx2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Dyz_vrr;
    Double I_KINETIC_Gxy2z_Dyz_vrr = PBZ*I_KINETIC_Gxy2z_Py_vrr+2*oned2z*I_KINETIC_Fxyz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Dyz_vrr;
    Double I_KINETIC_Gx3z_Dyz_vrr = PBZ*I_KINETIC_Gx3z_Py_vrr+3*oned2z*I_KINETIC_Fx2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dyz_vrr;
    Double I_KINETIC_G4y_Dyz_vrr = PBZ*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_KINETIC_G3yz_Dyz_vrr = PBZ*I_KINETIC_G3yz_Py_vrr+oned2z*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dyz_vrr;
    Double I_KINETIC_G2y2z_Dyz_vrr = PBZ*I_KINETIC_G2y2z_Py_vrr+2*oned2z*I_KINETIC_F2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dyz_vrr;
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
     * shell quartet name: SQ_KINETIC_H_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_P
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_KINETIC_H5x_D2x_vrr = PBX*I_KINETIC_H5x_Px_vrr+5*oned2z*I_KINETIC_G4x_Px_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2x_vrr-adz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_H4xy_D2x_vrr = PBX*I_KINETIC_H4xy_Px_vrr+4*oned2z*I_KINETIC_G3xy_Px_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2x_vrr-adz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_D2x_vrr = PBX*I_KINETIC_H4xz_Px_vrr+4*oned2z*I_KINETIC_G3xz_Px_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2x_vrr-adz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_D2x_vrr = PBX*I_KINETIC_H3x2y_Px_vrr+3*oned2z*I_KINETIC_G2x2y_Px_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2x_vrr-adz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_H3x2z_D2x_vrr = PBX*I_KINETIC_H3x2z_Px_vrr+3*oned2z*I_KINETIC_G2x2z_Px_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2x_vrr-adz*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_KINETIC_H2x3y_D2x_vrr = PBX*I_KINETIC_H2x3y_Px_vrr+2*oned2z*I_KINETIC_Gx3y_Px_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2x_vrr-adz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_H2x2yz_D2x_vrr = PBX*I_KINETIC_H2x2yz_Px_vrr+2*oned2z*I_KINETIC_Gx2yz_Px_vrr+oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2x3z_D2x_vrr = PBX*I_KINETIC_H2x3z_Px_vrr+2*oned2z*I_KINETIC_Gx3z_Px_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2x_vrr-adz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_Hx4y_D2x_vrr = PBX*I_KINETIC_Hx4y_Px_vrr+oned2z*I_KINETIC_G4y_Px_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx4z_D2x_vrr = PBX*I_KINETIC_Hx4z_Px_vrr+oned2z*I_KINETIC_G4z_Px_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_D2x_vrr = PBX*I_KINETIC_H5y_Px_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2x_vrr-adz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_H4yz_D2x_vrr = PBX*I_KINETIC_H4yz_Px_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2x_vrr-adz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_D2x_vrr = PBX*I_KINETIC_H3y2z_Px_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2x_vrr-adz*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_KINETIC_H2y3z_D2x_vrr = PBX*I_KINETIC_H2y3z_Px_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2x_vrr-adz*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_KINETIC_Hy4z_D2x_vrr = PBX*I_KINETIC_Hy4z_Px_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_D2x_vrr = PBX*I_KINETIC_H5z_Px_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2x_vrr-adz*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_KINETIC_H5x_Dxy_vrr = PBY*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dxy_vrr;
    Double I_KINETIC_H4xy_Dxy_vrr = PBY*I_KINETIC_H4xy_Px_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dxy_vrr;
    Double I_KINETIC_H4xz_Dxy_vrr = PBY*I_KINETIC_H4xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dxy_vrr;
    Double I_KINETIC_H3x2y_Dxy_vrr = PBY*I_KINETIC_H3x2y_Px_vrr+2*oned2z*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr;
    Double I_KINETIC_H3x2z_Dxy_vrr = PBY*I_KINETIC_H3x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dxy_vrr;
    Double I_KINETIC_H2x3y_Dxy_vrr = PBY*I_KINETIC_H2x3y_Px_vrr+3*oned2z*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dxy_vrr;
    Double I_KINETIC_H2x2yz_Dxy_vrr = PBY*I_KINETIC_H2x2yz_Px_vrr+2*oned2z*I_KINETIC_G2xyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr;
    Double I_KINETIC_H2x3z_Dxy_vrr = PBY*I_KINETIC_H2x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dxy_vrr;
    Double I_KINETIC_Hx4y_Dxy_vrr = PBY*I_KINETIC_Hx4y_Px_vrr+4*oned2z*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr;
    Double I_KINETIC_Hx4z_Dxy_vrr = PBY*I_KINETIC_Hx4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr;
    Double I_KINETIC_H5y_Dxy_vrr = PBY*I_KINETIC_H5y_Px_vrr+5*oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dxy_vrr;
    Double I_KINETIC_H4yz_Dxy_vrr = PBY*I_KINETIC_H4yz_Px_vrr+4*oned2z*I_KINETIC_G3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dxy_vrr;
    Double I_KINETIC_H3y2z_Dxy_vrr = PBY*I_KINETIC_H3y2z_Px_vrr+3*oned2z*I_KINETIC_G2y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dxy_vrr;
    Double I_KINETIC_H2y3z_Dxy_vrr = PBY*I_KINETIC_H2y3z_Px_vrr+2*oned2z*I_KINETIC_Gy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dxy_vrr;
    Double I_KINETIC_Hy4z_Dxy_vrr = PBY*I_KINETIC_Hy4z_Px_vrr+oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dxy_vrr;
    Double I_KINETIC_H5z_Dxy_vrr = PBY*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dxy_vrr;
    Double I_KINETIC_H5x_Dxz_vrr = PBZ*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dxz_vrr;
    Double I_KINETIC_H4xy_Dxz_vrr = PBZ*I_KINETIC_H4xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dxz_vrr;
    Double I_KINETIC_H4xz_Dxz_vrr = PBZ*I_KINETIC_H4xz_Px_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dxz_vrr;
    Double I_KINETIC_H3x2y_Dxz_vrr = PBZ*I_KINETIC_H3x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr;
    Double I_KINETIC_H3x2z_Dxz_vrr = PBZ*I_KINETIC_H3x2z_Px_vrr+2*oned2z*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dxz_vrr;
    Double I_KINETIC_H2x3y_Dxz_vrr = PBZ*I_KINETIC_H2x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dxz_vrr;
    Double I_KINETIC_H2x2yz_Dxz_vrr = PBZ*I_KINETIC_H2x2yz_Px_vrr+oned2z*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr;
    Double I_KINETIC_H2x3z_Dxz_vrr = PBZ*I_KINETIC_H2x3z_Px_vrr+3*oned2z*I_KINETIC_G2x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dxz_vrr;
    Double I_KINETIC_Hx4y_Dxz_vrr = PBZ*I_KINETIC_Hx4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr;
    Double I_KINETIC_Hx4z_Dxz_vrr = PBZ*I_KINETIC_Hx4z_Px_vrr+4*oned2z*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr;
    Double I_KINETIC_H5y_Dxz_vrr = PBZ*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dxz_vrr;
    Double I_KINETIC_H4yz_Dxz_vrr = PBZ*I_KINETIC_H4yz_Px_vrr+oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dxz_vrr;
    Double I_KINETIC_H3y2z_Dxz_vrr = PBZ*I_KINETIC_H3y2z_Px_vrr+2*oned2z*I_KINETIC_G3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dxz_vrr;
    Double I_KINETIC_H2y3z_Dxz_vrr = PBZ*I_KINETIC_H2y3z_Px_vrr+3*oned2z*I_KINETIC_G2y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dxz_vrr;
    Double I_KINETIC_Hy4z_Dxz_vrr = PBZ*I_KINETIC_Hy4z_Px_vrr+4*oned2z*I_KINETIC_Gy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dxz_vrr;
    Double I_KINETIC_H5z_Dxz_vrr = PBZ*I_KINETIC_H5z_Px_vrr+5*oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dxz_vrr;
    Double I_KINETIC_H5x_D2y_vrr = PBY*I_KINETIC_H5x_Py_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2y_vrr-adz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_H4xy_D2y_vrr = PBY*I_KINETIC_H4xy_Py_vrr+oned2z*I_KINETIC_G4x_Py_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2y_vrr-adz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_D2y_vrr = PBY*I_KINETIC_H4xz_Py_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2y_vrr-adz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_D2y_vrr = PBY*I_KINETIC_H3x2y_Py_vrr+2*oned2z*I_KINETIC_G3xy_Py_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2y_vrr-adz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_H3x2z_D2y_vrr = PBY*I_KINETIC_H3x2z_Py_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2y_vrr-adz*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_KINETIC_H2x3y_D2y_vrr = PBY*I_KINETIC_H2x3y_Py_vrr+3*oned2z*I_KINETIC_G2x2y_Py_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2y_vrr-adz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_H2x2yz_D2y_vrr = PBY*I_KINETIC_H2x2yz_Py_vrr+2*oned2z*I_KINETIC_G2xyz_Py_vrr+oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2x3z_D2y_vrr = PBY*I_KINETIC_H2x3z_Py_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2y_vrr-adz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_Hx4y_D2y_vrr = PBY*I_KINETIC_Hx4y_Py_vrr+4*oned2z*I_KINETIC_Gx3y_Py_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx4z_D2y_vrr = PBY*I_KINETIC_Hx4z_Py_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_D2y_vrr = PBY*I_KINETIC_H5y_Py_vrr+5*oned2z*I_KINETIC_G4y_Py_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2y_vrr-adz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_H4yz_D2y_vrr = PBY*I_KINETIC_H4yz_Py_vrr+4*oned2z*I_KINETIC_G3yz_Py_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2y_vrr-adz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_D2y_vrr = PBY*I_KINETIC_H3y2z_Py_vrr+3*oned2z*I_KINETIC_G2y2z_Py_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2y_vrr-adz*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_KINETIC_H2y3z_D2y_vrr = PBY*I_KINETIC_H2y3z_Py_vrr+2*oned2z*I_KINETIC_Gy3z_Py_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2y_vrr-adz*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_KINETIC_Hy4z_D2y_vrr = PBY*I_KINETIC_Hy4z_Py_vrr+oned2z*I_KINETIC_G4z_Py_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_D2y_vrr = PBY*I_KINETIC_H5z_Py_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2y_vrr-adz*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_KINETIC_H5x_Dyz_vrr = PBZ*I_KINETIC_H5x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dyz_vrr;
    Double I_KINETIC_H4xy_Dyz_vrr = PBZ*I_KINETIC_H4xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dyz_vrr;
    Double I_KINETIC_H4xz_Dyz_vrr = PBZ*I_KINETIC_H4xz_Py_vrr+oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dyz_vrr;
    Double I_KINETIC_H3x2y_Dyz_vrr = PBZ*I_KINETIC_H3x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr;
    Double I_KINETIC_H3x2z_Dyz_vrr = PBZ*I_KINETIC_H3x2z_Py_vrr+2*oned2z*I_KINETIC_G3xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dyz_vrr;
    Double I_KINETIC_H2x3y_Dyz_vrr = PBZ*I_KINETIC_H2x3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dyz_vrr;
    Double I_KINETIC_H2x2yz_Dyz_vrr = PBZ*I_KINETIC_H2x2yz_Py_vrr+oned2z*I_KINETIC_G2x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr;
    Double I_KINETIC_H2x3z_Dyz_vrr = PBZ*I_KINETIC_H2x3z_Py_vrr+3*oned2z*I_KINETIC_G2x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dyz_vrr;
    Double I_KINETIC_Hx4y_Dyz_vrr = PBZ*I_KINETIC_Hx4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr;
    Double I_KINETIC_Hx4z_Dyz_vrr = PBZ*I_KINETIC_Hx4z_Py_vrr+4*oned2z*I_KINETIC_Gx3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr;
    Double I_KINETIC_H5y_Dyz_vrr = PBZ*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dyz_vrr;
    Double I_KINETIC_H4yz_Dyz_vrr = PBZ*I_KINETIC_H4yz_Py_vrr+oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dyz_vrr;
    Double I_KINETIC_H3y2z_Dyz_vrr = PBZ*I_KINETIC_H3y2z_Py_vrr+2*oned2z*I_KINETIC_G3yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dyz_vrr;
    Double I_KINETIC_H2y3z_Dyz_vrr = PBZ*I_KINETIC_H2y3z_Py_vrr+3*oned2z*I_KINETIC_G2y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dyz_vrr;
    Double I_KINETIC_Hy4z_Dyz_vrr = PBZ*I_KINETIC_Hy4z_Py_vrr+4*oned2z*I_KINETIC_Gy3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dyz_vrr;
    Double I_KINETIC_H5z_Dyz_vrr = PBZ*I_KINETIC_H5z_Py_vrr+5*oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dyz_vrr;
    Double I_KINETIC_H5x_D2z_vrr = PBZ*I_KINETIC_H5x_Pz_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2z_vrr-adz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_H4xy_D2z_vrr = PBZ*I_KINETIC_H4xy_Pz_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2z_vrr-adz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_D2z_vrr = PBZ*I_KINETIC_H4xz_Pz_vrr+oned2z*I_KINETIC_G4x_Pz_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2z_vrr-adz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_D2z_vrr = PBZ*I_KINETIC_H3x2y_Pz_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2z_vrr-adz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_H3x2z_D2z_vrr = PBZ*I_KINETIC_H3x2z_Pz_vrr+2*oned2z*I_KINETIC_G3xz_Pz_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2z_vrr-adz*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_KINETIC_H2x3y_D2z_vrr = PBZ*I_KINETIC_H2x3y_Pz_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2z_vrr-adz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_H2x2yz_D2z_vrr = PBZ*I_KINETIC_H2x2yz_Pz_vrr+oned2z*I_KINETIC_G2x2y_Pz_vrr+oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2x3z_D2z_vrr = PBZ*I_KINETIC_H2x3z_Pz_vrr+3*oned2z*I_KINETIC_G2x2z_Pz_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2z_vrr-adz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_Hx4y_D2z_vrr = PBZ*I_KINETIC_Hx4y_Pz_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx4z_D2z_vrr = PBZ*I_KINETIC_Hx4z_Pz_vrr+4*oned2z*I_KINETIC_Gx3z_Pz_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_D2z_vrr = PBZ*I_KINETIC_H5y_Pz_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2z_vrr-adz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_H4yz_D2z_vrr = PBZ*I_KINETIC_H4yz_Pz_vrr+oned2z*I_KINETIC_G4y_Pz_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2z_vrr-adz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_D2z_vrr = PBZ*I_KINETIC_H3y2z_Pz_vrr+2*oned2z*I_KINETIC_G3yz_Pz_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2z_vrr-adz*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_KINETIC_H2y3z_D2z_vrr = PBZ*I_KINETIC_H2y3z_Pz_vrr+3*oned2z*I_KINETIC_G2y2z_Pz_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2z_vrr-adz*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_KINETIC_Hy4z_D2z_vrr = PBZ*I_KINETIC_Hy4z_Pz_vrr+4*oned2z*I_KINETIC_Gy3z_Pz_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_D2z_vrr = PBZ*I_KINETIC_H5z_Pz_vrr+5*oned2z*I_KINETIC_G4z_Pz_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2z_vrr-adz*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_D
     * RHS shell quartet name: SQ_KINETIC_G_D
     * RHS shell quartet name: SQ_KINETIC_H_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     ************************************************************/
    Double I_KINETIC_I6x_D2x_vrr = PAX*I_KINETIC_H5x_D2x_vrr+5*oned2z*I_KINETIC_G4x_D2x_vrr+2*oned2z*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6x_D2x_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_KINETIC_I5xy_D2x_vrr = PAY*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_D2x_vrr;
    Double I_KINETIC_I5xz_D2x_vrr = PAZ*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_D2x_vrr;
    Double I_KINETIC_I4x2y_D2x_vrr = PAY*I_KINETIC_H4xy_D2x_vrr+oned2z*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_D2x_vrr-bdz*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_KINETIC_I4xyz_D2x_vrr = PAZ*I_KINETIC_H4xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_D2x_vrr;
    Double I_KINETIC_I4x2z_D2x_vrr = PAZ*I_KINETIC_H4xz_D2x_vrr+oned2z*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_KINETIC_I3x3y_D2x_vrr = PAY*I_KINETIC_H3x2y_D2x_vrr+2*oned2z*I_KINETIC_G3xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_KINETIC_I3x2yz_D2x_vrr = PAZ*I_KINETIC_H3x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_D2x_vrr;
    Double I_KINETIC_I3xy2z_D2x_vrr = PAY*I_KINETIC_H3x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_D2x_vrr;
    Double I_KINETIC_I3x3z_D2x_vrr = PAZ*I_KINETIC_H3x2z_D2x_vrr+2*oned2z*I_KINETIC_G3xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_D2x_vrr;
    Double I_KINETIC_I2x4y_D2x_vrr = PAX*I_KINETIC_Hx4y_D2x_vrr+oned2z*I_KINETIC_G4y_D2x_vrr+2*oned2z*I_KINETIC_Hx4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_D2x_vrr-bdz*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_KINETIC_I2x3yz_D2x_vrr = PAZ*I_KINETIC_H2x3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_D2x_vrr;
    Double I_KINETIC_I2x2y2z_D2x_vrr = PAZ*I_KINETIC_H2x2yz_D2x_vrr+oned2z*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_KINETIC_I2xy3z_D2x_vrr = PAY*I_KINETIC_H2x3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_D2x_vrr;
    Double I_KINETIC_I2x4z_D2x_vrr = PAX*I_KINETIC_Hx4z_D2x_vrr+oned2z*I_KINETIC_G4z_D2x_vrr+2*oned2z*I_KINETIC_Hx4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_D2x_vrr-bdz*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_KINETIC_Ix5y_D2x_vrr = PAX*I_KINETIC_H5y_D2x_vrr+2*oned2z*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_D2x_vrr;
    Double I_KINETIC_Ix4yz_D2x_vrr = PAZ*I_KINETIC_Hx4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_D2x_vrr;
    Double I_KINETIC_Ix3y2z_D2x_vrr = PAX*I_KINETIC_H3y2z_D2x_vrr+2*oned2z*I_KINETIC_H3y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_D2x_vrr;
    Double I_KINETIC_Ix2y3z_D2x_vrr = PAX*I_KINETIC_H2y3z_D2x_vrr+2*oned2z*I_KINETIC_H2y3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_D2x_vrr;
    Double I_KINETIC_Ixy4z_D2x_vrr = PAY*I_KINETIC_Hx4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_D2x_vrr;
    Double I_KINETIC_Ix5z_D2x_vrr = PAX*I_KINETIC_H5z_D2x_vrr+2*oned2z*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_D2x_vrr;
    Double I_KINETIC_I6y_D2x_vrr = PAY*I_KINETIC_H5y_D2x_vrr+5*oned2z*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I6y_D2x_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_KINETIC_I5yz_D2x_vrr = PAZ*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_D2x_vrr;
    Double I_KINETIC_I4y2z_D2x_vrr = PAZ*I_KINETIC_H4yz_D2x_vrr+oned2z*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_KINETIC_I3y3z_D2x_vrr = PAZ*I_KINETIC_H3y2z_D2x_vrr+2*oned2z*I_KINETIC_G3yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_D2x_vrr;
    Double I_KINETIC_I2y4z_D2x_vrr = PAY*I_KINETIC_Hy4z_D2x_vrr+oned2z*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_D2x_vrr-bdz*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_KINETIC_Iy5z_D2x_vrr = PAY*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_D2x_vrr;
    Double I_KINETIC_I6z_D2x_vrr = PAZ*I_KINETIC_H5z_D2x_vrr+5*oned2z*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_I6z_D2x_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_KINETIC_I6x_Dxy_vrr = PAX*I_KINETIC_H5x_Dxy_vrr+5*oned2z*I_KINETIC_G4x_Dxy_vrr+oned2z*I_KINETIC_H5x_Py_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Dxy_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_KINETIC_I5xy_Dxy_vrr = PAY*I_KINETIC_H5x_Dxy_vrr+oned2z*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Dxy_vrr;
    Double I_KINETIC_I5xz_Dxy_vrr = PAZ*I_KINETIC_H5x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Dxy_vrr;
    Double I_KINETIC_I4x2y_Dxy_vrr = PAY*I_KINETIC_H4xy_Dxy_vrr+oned2z*I_KINETIC_G4x_Dxy_vrr+oned2z*I_KINETIC_H4xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_KINETIC_I4xyz_Dxy_vrr = PAZ*I_KINETIC_H4xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Dxy_vrr;
    Double I_KINETIC_I4x2z_Dxy_vrr = PAZ*I_KINETIC_H4xz_Dxy_vrr+oned2z*I_KINETIC_G4x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_KINETIC_I3x3y_Dxy_vrr = PAY*I_KINETIC_H3x2y_Dxy_vrr+2*oned2z*I_KINETIC_G3xy_Dxy_vrr+oned2z*I_KINETIC_H3x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Dxy_vrr;
    Double I_KINETIC_I3x2yz_Dxy_vrr = PAZ*I_KINETIC_H3x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Dxy_vrr;
    Double I_KINETIC_I3xy2z_Dxy_vrr = PAY*I_KINETIC_H3x2z_Dxy_vrr+oned2z*I_KINETIC_H3x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Dxy_vrr;
    Double I_KINETIC_I3x3z_Dxy_vrr = PAZ*I_KINETIC_H3x2z_Dxy_vrr+2*oned2z*I_KINETIC_G3xz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Dxy_vrr;
    Double I_KINETIC_I2x4y_Dxy_vrr = PAX*I_KINETIC_Hx4y_Dxy_vrr+oned2z*I_KINETIC_G4y_Dxy_vrr+oned2z*I_KINETIC_Hx4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_KINETIC_I2x3yz_Dxy_vrr = PAZ*I_KINETIC_H2x3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Dxy_vrr;
    Double I_KINETIC_I2x2y2z_Dxy_vrr = PAZ*I_KINETIC_H2x2yz_Dxy_vrr+oned2z*I_KINETIC_G2x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_KINETIC_I2xy3z_Dxy_vrr = PAY*I_KINETIC_H2x3z_Dxy_vrr+oned2z*I_KINETIC_H2x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Dxy_vrr;
    Double I_KINETIC_I2x4z_Dxy_vrr = PAX*I_KINETIC_Hx4z_Dxy_vrr+oned2z*I_KINETIC_G4z_Dxy_vrr+oned2z*I_KINETIC_Hx4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_KINETIC_Ix5y_Dxy_vrr = PAX*I_KINETIC_H5y_Dxy_vrr+oned2z*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Dxy_vrr;
    Double I_KINETIC_Ix4yz_Dxy_vrr = PAZ*I_KINETIC_Hx4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Dxy_vrr;
    Double I_KINETIC_Ix3y2z_Dxy_vrr = PAX*I_KINETIC_H3y2z_Dxy_vrr+oned2z*I_KINETIC_H3y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Dxy_vrr;
    Double I_KINETIC_Ix2y3z_Dxy_vrr = PAX*I_KINETIC_H2y3z_Dxy_vrr+oned2z*I_KINETIC_H2y3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Dxy_vrr;
    Double I_KINETIC_Ixy4z_Dxy_vrr = PAY*I_KINETIC_Hx4z_Dxy_vrr+oned2z*I_KINETIC_Hx4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Dxy_vrr;
    Double I_KINETIC_Ix5z_Dxy_vrr = PAX*I_KINETIC_H5z_Dxy_vrr+oned2z*I_KINETIC_H5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Dxy_vrr;
    Double I_KINETIC_I6y_Dxy_vrr = PAY*I_KINETIC_H5y_Dxy_vrr+5*oned2z*I_KINETIC_G4y_Dxy_vrr+oned2z*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Dxy_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_KINETIC_I5yz_Dxy_vrr = PAZ*I_KINETIC_H5y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Dxy_vrr;
    Double I_KINETIC_I4y2z_Dxy_vrr = PAZ*I_KINETIC_H4yz_Dxy_vrr+oned2z*I_KINETIC_G4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_KINETIC_I3y3z_Dxy_vrr = PAZ*I_KINETIC_H3y2z_Dxy_vrr+2*oned2z*I_KINETIC_G3yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Dxy_vrr;
    Double I_KINETIC_I2y4z_Dxy_vrr = PAY*I_KINETIC_Hy4z_Dxy_vrr+oned2z*I_KINETIC_G4z_Dxy_vrr+oned2z*I_KINETIC_Hy4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_KINETIC_Iy5z_Dxy_vrr = PAY*I_KINETIC_H5z_Dxy_vrr+oned2z*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Dxy_vrr;
    Double I_KINETIC_I6z_Dxy_vrr = PAZ*I_KINETIC_H5z_Dxy_vrr+5*oned2z*I_KINETIC_G4z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Dxy_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_KINETIC_I6x_Dxz_vrr = PAX*I_KINETIC_H5x_Dxz_vrr+5*oned2z*I_KINETIC_G4x_Dxz_vrr+oned2z*I_KINETIC_H5x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Dxz_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_KINETIC_I5xy_Dxz_vrr = PAY*I_KINETIC_H5x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Dxz_vrr;
    Double I_KINETIC_I5xz_Dxz_vrr = PAZ*I_KINETIC_H5x_Dxz_vrr+oned2z*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Dxz_vrr;
    Double I_KINETIC_I4x2y_Dxz_vrr = PAY*I_KINETIC_H4xy_Dxz_vrr+oned2z*I_KINETIC_G4x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_KINETIC_I4xyz_Dxz_vrr = PAZ*I_KINETIC_H4xy_Dxz_vrr+oned2z*I_KINETIC_H4xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Dxz_vrr;
    Double I_KINETIC_I4x2z_Dxz_vrr = PAZ*I_KINETIC_H4xz_Dxz_vrr+oned2z*I_KINETIC_G4x_Dxz_vrr+oned2z*I_KINETIC_H4xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_KINETIC_I3x3y_Dxz_vrr = PAY*I_KINETIC_H3x2y_Dxz_vrr+2*oned2z*I_KINETIC_G3xy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Dxz_vrr;
    Double I_KINETIC_I3x2yz_Dxz_vrr = PAZ*I_KINETIC_H3x2y_Dxz_vrr+oned2z*I_KINETIC_H3x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Dxz_vrr;
    Double I_KINETIC_I3xy2z_Dxz_vrr = PAY*I_KINETIC_H3x2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Dxz_vrr;
    Double I_KINETIC_I3x3z_Dxz_vrr = PAZ*I_KINETIC_H3x2z_Dxz_vrr+2*oned2z*I_KINETIC_G3xz_Dxz_vrr+oned2z*I_KINETIC_H3x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Dxz_vrr;
    Double I_KINETIC_I2x4y_Dxz_vrr = PAX*I_KINETIC_Hx4y_Dxz_vrr+oned2z*I_KINETIC_G4y_Dxz_vrr+oned2z*I_KINETIC_Hx4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_KINETIC_I2x3yz_Dxz_vrr = PAZ*I_KINETIC_H2x3y_Dxz_vrr+oned2z*I_KINETIC_H2x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Dxz_vrr;
    Double I_KINETIC_I2x2y2z_Dxz_vrr = PAZ*I_KINETIC_H2x2yz_Dxz_vrr+oned2z*I_KINETIC_G2x2y_Dxz_vrr+oned2z*I_KINETIC_H2x2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_Dxz_vrr;
    Double I_KINETIC_I2xy3z_Dxz_vrr = PAY*I_KINETIC_H2x3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Dxz_vrr;
    Double I_KINETIC_I2x4z_Dxz_vrr = PAX*I_KINETIC_Hx4z_Dxz_vrr+oned2z*I_KINETIC_G4z_Dxz_vrr+oned2z*I_KINETIC_Hx4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_KINETIC_Ix5y_Dxz_vrr = PAX*I_KINETIC_H5y_Dxz_vrr+oned2z*I_KINETIC_H5y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Dxz_vrr;
    Double I_KINETIC_Ix4yz_Dxz_vrr = PAZ*I_KINETIC_Hx4y_Dxz_vrr+oned2z*I_KINETIC_Hx4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Dxz_vrr;
    Double I_KINETIC_Ix3y2z_Dxz_vrr = PAX*I_KINETIC_H3y2z_Dxz_vrr+oned2z*I_KINETIC_H3y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Dxz_vrr;
    Double I_KINETIC_Ix2y3z_Dxz_vrr = PAX*I_KINETIC_H2y3z_Dxz_vrr+oned2z*I_KINETIC_H2y3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Dxz_vrr;
    Double I_KINETIC_Ixy4z_Dxz_vrr = PAY*I_KINETIC_Hx4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Dxz_vrr;
    Double I_KINETIC_Ix5z_Dxz_vrr = PAX*I_KINETIC_H5z_Dxz_vrr+oned2z*I_KINETIC_H5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Dxz_vrr;
    Double I_KINETIC_I6y_Dxz_vrr = PAY*I_KINETIC_H5y_Dxz_vrr+5*oned2z*I_KINETIC_G4y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Dxz_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_KINETIC_I5yz_Dxz_vrr = PAZ*I_KINETIC_H5y_Dxz_vrr+oned2z*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Dxz_vrr;
    Double I_KINETIC_I4y2z_Dxz_vrr = PAZ*I_KINETIC_H4yz_Dxz_vrr+oned2z*I_KINETIC_G4y_Dxz_vrr+oned2z*I_KINETIC_H4yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_KINETIC_I3y3z_Dxz_vrr = PAZ*I_KINETIC_H3y2z_Dxz_vrr+2*oned2z*I_KINETIC_G3yz_Dxz_vrr+oned2z*I_KINETIC_H3y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Dxz_vrr;
    Double I_KINETIC_I2y4z_Dxz_vrr = PAY*I_KINETIC_Hy4z_Dxz_vrr+oned2z*I_KINETIC_G4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_KINETIC_Iy5z_Dxz_vrr = PAY*I_KINETIC_H5z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Dxz_vrr;
    Double I_KINETIC_I6z_Dxz_vrr = PAZ*I_KINETIC_H5z_Dxz_vrr+5*oned2z*I_KINETIC_G4z_Dxz_vrr+oned2z*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Dxz_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_KINETIC_I6x_D2y_vrr = PAX*I_KINETIC_H5x_D2y_vrr+5*oned2z*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I6x_D2y_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_KINETIC_I5xy_D2y_vrr = PAY*I_KINETIC_H5x_D2y_vrr+2*oned2z*I_KINETIC_H5x_Py_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_D2y_vrr;
    Double I_KINETIC_I5xz_D2y_vrr = PAZ*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_D2y_vrr;
    Double I_KINETIC_I4x2y_D2y_vrr = PAY*I_KINETIC_H4xy_D2y_vrr+oned2z*I_KINETIC_G4x_D2y_vrr+2*oned2z*I_KINETIC_H4xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_D2y_vrr-bdz*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_KINETIC_I4xyz_D2y_vrr = PAZ*I_KINETIC_H4xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_D2y_vrr;
    Double I_KINETIC_I4x2z_D2y_vrr = PAZ*I_KINETIC_H4xz_D2y_vrr+oned2z*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_KINETIC_I3x3y_D2y_vrr = PAY*I_KINETIC_H3x2y_D2y_vrr+2*oned2z*I_KINETIC_G3xy_D2y_vrr+2*oned2z*I_KINETIC_H3x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_KINETIC_I3x2yz_D2y_vrr = PAZ*I_KINETIC_H3x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_D2y_vrr;
    Double I_KINETIC_I3xy2z_D2y_vrr = PAY*I_KINETIC_H3x2z_D2y_vrr+2*oned2z*I_KINETIC_H3x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_D2y_vrr;
    Double I_KINETIC_I3x3z_D2y_vrr = PAZ*I_KINETIC_H3x2z_D2y_vrr+2*oned2z*I_KINETIC_G3xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_D2y_vrr;
    Double I_KINETIC_I2x4y_D2y_vrr = PAX*I_KINETIC_Hx4y_D2y_vrr+oned2z*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_D2y_vrr-bdz*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_KINETIC_I2x3yz_D2y_vrr = PAZ*I_KINETIC_H2x3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_D2y_vrr;
    Double I_KINETIC_I2x2y2z_D2y_vrr = PAZ*I_KINETIC_H2x2yz_D2y_vrr+oned2z*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_KINETIC_I2xy3z_D2y_vrr = PAY*I_KINETIC_H2x3z_D2y_vrr+2*oned2z*I_KINETIC_H2x3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_D2y_vrr;
    Double I_KINETIC_I2x4z_D2y_vrr = PAX*I_KINETIC_Hx4z_D2y_vrr+oned2z*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_D2y_vrr-bdz*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_KINETIC_Ix5y_D2y_vrr = PAX*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_D2y_vrr;
    Double I_KINETIC_Ix4yz_D2y_vrr = PAZ*I_KINETIC_Hx4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_D2y_vrr;
    Double I_KINETIC_Ix3y2z_D2y_vrr = PAX*I_KINETIC_H3y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_D2y_vrr;
    Double I_KINETIC_Ix2y3z_D2y_vrr = PAX*I_KINETIC_H2y3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_D2y_vrr;
    Double I_KINETIC_Ixy4z_D2y_vrr = PAY*I_KINETIC_Hx4z_D2y_vrr+2*oned2z*I_KINETIC_Hx4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_D2y_vrr;
    Double I_KINETIC_Ix5z_D2y_vrr = PAX*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_D2y_vrr;
    Double I_KINETIC_I6y_D2y_vrr = PAY*I_KINETIC_H5y_D2y_vrr+5*oned2z*I_KINETIC_G4y_D2y_vrr+2*oned2z*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I6y_D2y_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_KINETIC_I5yz_D2y_vrr = PAZ*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_D2y_vrr;
    Double I_KINETIC_I4y2z_D2y_vrr = PAZ*I_KINETIC_H4yz_D2y_vrr+oned2z*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_KINETIC_I3y3z_D2y_vrr = PAZ*I_KINETIC_H3y2z_D2y_vrr+2*oned2z*I_KINETIC_G3yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_D2y_vrr;
    Double I_KINETIC_I2y4z_D2y_vrr = PAY*I_KINETIC_Hy4z_D2y_vrr+oned2z*I_KINETIC_G4z_D2y_vrr+2*oned2z*I_KINETIC_Hy4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_D2y_vrr-bdz*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_KINETIC_Iy5z_D2y_vrr = PAY*I_KINETIC_H5z_D2y_vrr+2*oned2z*I_KINETIC_H5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_D2y_vrr;
    Double I_KINETIC_I6z_D2y_vrr = PAZ*I_KINETIC_H5z_D2y_vrr+5*oned2z*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_I6z_D2y_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_KINETIC_I6x_Dyz_vrr = PAX*I_KINETIC_H5x_Dyz_vrr+5*oned2z*I_KINETIC_G4x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Dyz_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_Dyz_vrr;
    Double I_KINETIC_I5xy_Dyz_vrr = PAY*I_KINETIC_H5x_Dyz_vrr+oned2z*I_KINETIC_H5x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Dyz_vrr;
    Double I_KINETIC_I5xz_Dyz_vrr = PAZ*I_KINETIC_H5x_Dyz_vrr+oned2z*I_KINETIC_H5x_Py_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Dyz_vrr;
    Double I_KINETIC_I4x2y_Dyz_vrr = PAY*I_KINETIC_H4xy_Dyz_vrr+oned2z*I_KINETIC_G4x_Dyz_vrr+oned2z*I_KINETIC_H4xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_G4x_Dyz_vrr;
    Double I_KINETIC_I4xyz_Dyz_vrr = PAZ*I_KINETIC_H4xy_Dyz_vrr+oned2z*I_KINETIC_H4xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Dyz_vrr;
    Double I_KINETIC_I4x2z_Dyz_vrr = PAZ*I_KINETIC_H4xz_Dyz_vrr+oned2z*I_KINETIC_G4x_Dyz_vrr+oned2z*I_KINETIC_H4xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_G4x_Dyz_vrr;
    Double I_KINETIC_I3x3y_Dyz_vrr = PAY*I_KINETIC_H3x2y_Dyz_vrr+2*oned2z*I_KINETIC_G3xy_Dyz_vrr+oned2z*I_KINETIC_H3x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_Dyz_vrr;
    Double I_KINETIC_I3x2yz_Dyz_vrr = PAZ*I_KINETIC_H3x2y_Dyz_vrr+oned2z*I_KINETIC_H3x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Dyz_vrr;
    Double I_KINETIC_I3xy2z_Dyz_vrr = PAY*I_KINETIC_H3x2z_Dyz_vrr+oned2z*I_KINETIC_H3x2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Dyz_vrr;
    Double I_KINETIC_I3x3z_Dyz_vrr = PAZ*I_KINETIC_H3x2z_Dyz_vrr+2*oned2z*I_KINETIC_G3xz_Dyz_vrr+oned2z*I_KINETIC_H3x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_Dyz_vrr;
    Double I_KINETIC_I2x4y_Dyz_vrr = PAX*I_KINETIC_Hx4y_Dyz_vrr+oned2z*I_KINETIC_G4y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_KINETIC_I2x3yz_Dyz_vrr = PAZ*I_KINETIC_H2x3y_Dyz_vrr+oned2z*I_KINETIC_H2x3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Dyz_vrr;
    Double I_KINETIC_I2x2y2z_Dyz_vrr = PAZ*I_KINETIC_H2x2yz_Dyz_vrr+oned2z*I_KINETIC_G2x2y_Dyz_vrr+oned2z*I_KINETIC_H2x2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_Dyz_vrr;
    Double I_KINETIC_I2xy3z_Dyz_vrr = PAY*I_KINETIC_H2x3z_Dyz_vrr+oned2z*I_KINETIC_H2x3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Dyz_vrr;
    Double I_KINETIC_I2x4z_Dyz_vrr = PAX*I_KINETIC_Hx4z_Dyz_vrr+oned2z*I_KINETIC_G4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_KINETIC_Ix5y_Dyz_vrr = PAX*I_KINETIC_H5y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Dyz_vrr;
    Double I_KINETIC_Ix4yz_Dyz_vrr = PAZ*I_KINETIC_Hx4y_Dyz_vrr+oned2z*I_KINETIC_Hx4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Dyz_vrr;
    Double I_KINETIC_Ix3y2z_Dyz_vrr = PAX*I_KINETIC_H3y2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Dyz_vrr;
    Double I_KINETIC_Ix2y3z_Dyz_vrr = PAX*I_KINETIC_H2y3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Dyz_vrr;
    Double I_KINETIC_Ixy4z_Dyz_vrr = PAY*I_KINETIC_Hx4z_Dyz_vrr+oned2z*I_KINETIC_Hx4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Dyz_vrr;
    Double I_KINETIC_Ix5z_Dyz_vrr = PAX*I_KINETIC_H5z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Dyz_vrr;
    Double I_KINETIC_I6y_Dyz_vrr = PAY*I_KINETIC_H5y_Dyz_vrr+5*oned2z*I_KINETIC_G4y_Dyz_vrr+oned2z*I_KINETIC_H5y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Dyz_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_KINETIC_I5yz_Dyz_vrr = PAZ*I_KINETIC_H5y_Dyz_vrr+oned2z*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Dyz_vrr;
    Double I_KINETIC_I4y2z_Dyz_vrr = PAZ*I_KINETIC_H4yz_Dyz_vrr+oned2z*I_KINETIC_G4y_Dyz_vrr+oned2z*I_KINETIC_H4yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_KINETIC_I3y3z_Dyz_vrr = PAZ*I_KINETIC_H3y2z_Dyz_vrr+2*oned2z*I_KINETIC_G3yz_Dyz_vrr+oned2z*I_KINETIC_H3y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_Dyz_vrr;
    Double I_KINETIC_I2y4z_Dyz_vrr = PAY*I_KINETIC_Hy4z_Dyz_vrr+oned2z*I_KINETIC_G4z_Dyz_vrr+oned2z*I_KINETIC_Hy4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_KINETIC_Iy5z_Dyz_vrr = PAY*I_KINETIC_H5z_Dyz_vrr+oned2z*I_KINETIC_H5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Dyz_vrr;
    Double I_KINETIC_I6z_Dyz_vrr = PAZ*I_KINETIC_H5z_Dyz_vrr+5*oned2z*I_KINETIC_G4z_Dyz_vrr+oned2z*I_KINETIC_H5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Dyz_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_KINETIC_I6x_D2z_vrr = PAX*I_KINETIC_H5x_D2z_vrr+5*oned2z*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I6x_D2z_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_KINETIC_I5xy_D2z_vrr = PAY*I_KINETIC_H5x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_D2z_vrr;
    Double I_KINETIC_I5xz_D2z_vrr = PAZ*I_KINETIC_H5x_D2z_vrr+2*oned2z*I_KINETIC_H5x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_D2z_vrr;
    Double I_KINETIC_I4x2y_D2z_vrr = PAY*I_KINETIC_H4xy_D2z_vrr+oned2z*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_D2z_vrr-bdz*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_KINETIC_I4xyz_D2z_vrr = PAZ*I_KINETIC_H4xy_D2z_vrr+2*oned2z*I_KINETIC_H4xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_D2z_vrr;
    Double I_KINETIC_I4x2z_D2z_vrr = PAZ*I_KINETIC_H4xz_D2z_vrr+oned2z*I_KINETIC_G4x_D2z_vrr+2*oned2z*I_KINETIC_H4xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_KINETIC_I3x3y_D2z_vrr = PAY*I_KINETIC_H3x2y_D2z_vrr+2*oned2z*I_KINETIC_G3xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_D2z_vrr;
    Double I_KINETIC_I3x2yz_D2z_vrr = PAZ*I_KINETIC_H3x2y_D2z_vrr+2*oned2z*I_KINETIC_H3x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_D2z_vrr;
    Double I_KINETIC_I3xy2z_D2z_vrr = PAY*I_KINETIC_H3x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_D2z_vrr;
    Double I_KINETIC_I3x3z_D2z_vrr = PAZ*I_KINETIC_H3x2z_D2z_vrr+2*oned2z*I_KINETIC_G3xz_D2z_vrr+2*oned2z*I_KINETIC_H3x2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_D2z_vrr;
    Double I_KINETIC_I2x4y_D2z_vrr = PAX*I_KINETIC_Hx4y_D2z_vrr+oned2z*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_D2z_vrr-bdz*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_KINETIC_I2x3yz_D2z_vrr = PAZ*I_KINETIC_H2x3y_D2z_vrr+2*oned2z*I_KINETIC_H2x3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_D2z_vrr;
    Double I_KINETIC_I2x2y2z_D2z_vrr = PAZ*I_KINETIC_H2x2yz_D2z_vrr+oned2z*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_H2x2yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_D2z_vrr;
    Double I_KINETIC_I2xy3z_D2z_vrr = PAY*I_KINETIC_H2x3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_D2z_vrr;
    Double I_KINETIC_I2x4z_D2z_vrr = PAX*I_KINETIC_Hx4z_D2z_vrr+oned2z*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_D2z_vrr-bdz*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_KINETIC_Ix5y_D2z_vrr = PAX*I_KINETIC_H5y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_D2z_vrr;
    Double I_KINETIC_Ix4yz_D2z_vrr = PAZ*I_KINETIC_Hx4y_D2z_vrr+2*oned2z*I_KINETIC_Hx4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_D2z_vrr;
    Double I_KINETIC_Ix3y2z_D2z_vrr = PAX*I_KINETIC_H3y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_D2z_vrr;
    Double I_KINETIC_Ix2y3z_D2z_vrr = PAX*I_KINETIC_H2y3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_D2z_vrr;
    Double I_KINETIC_Ixy4z_D2z_vrr = PAY*I_KINETIC_Hx4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_D2z_vrr;
    Double I_KINETIC_Ix5z_D2z_vrr = PAX*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_D2z_vrr;
    Double I_KINETIC_I6y_D2z_vrr = PAY*I_KINETIC_H5y_D2z_vrr+5*oned2z*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I6y_D2z_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_KINETIC_I5yz_D2z_vrr = PAZ*I_KINETIC_H5y_D2z_vrr+2*oned2z*I_KINETIC_H5y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_D2z_vrr;
    Double I_KINETIC_I4y2z_D2z_vrr = PAZ*I_KINETIC_H4yz_D2z_vrr+oned2z*I_KINETIC_G4y_D2z_vrr+2*oned2z*I_KINETIC_H4yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_KINETIC_I3y3z_D2z_vrr = PAZ*I_KINETIC_H3y2z_D2z_vrr+2*oned2z*I_KINETIC_G3yz_D2z_vrr+2*oned2z*I_KINETIC_H3y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_D2z_vrr;
    Double I_KINETIC_I2y4z_D2z_vrr = PAY*I_KINETIC_Hy4z_D2z_vrr+oned2z*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_D2z_vrr-bdz*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_KINETIC_Iy5z_D2z_vrr = PAY*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_D2z_vrr;
    Double I_KINETIC_I6z_D2z_vrr = PAZ*I_KINETIC_H5z_D2z_vrr+5*oned2z*I_KINETIC_G4z_D2z_vrr+2*oned2z*I_KINETIC_H5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_I6z_D2z_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_D_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_I_D_a_coefs = alpha;
    I_KINETIC_I6x_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6x_D2x_vrr;
    I_KINETIC_I5xy_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xy_D2x_vrr;
    I_KINETIC_I5xz_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xz_D2x_vrr;
    I_KINETIC_I4x2y_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2y_D2x_vrr;
    I_KINETIC_I4xyz_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4xyz_D2x_vrr;
    I_KINETIC_I4x2z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2z_D2x_vrr;
    I_KINETIC_I3x3y_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3y_D2x_vrr;
    I_KINETIC_I3x2yz_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x2yz_D2x_vrr;
    I_KINETIC_I3xy2z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3xy2z_D2x_vrr;
    I_KINETIC_I3x3z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3z_D2x_vrr;
    I_KINETIC_I2x4y_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4y_D2x_vrr;
    I_KINETIC_I2x3yz_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x3yz_D2x_vrr;
    I_KINETIC_I2x2y2z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x2y2z_D2x_vrr;
    I_KINETIC_I2xy3z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2xy3z_D2x_vrr;
    I_KINETIC_I2x4z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4z_D2x_vrr;
    I_KINETIC_Ix5y_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5y_D2x_vrr;
    I_KINETIC_Ix4yz_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix4yz_D2x_vrr;
    I_KINETIC_Ix3y2z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix3y2z_D2x_vrr;
    I_KINETIC_Ix2y3z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix2y3z_D2x_vrr;
    I_KINETIC_Ixy4z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ixy4z_D2x_vrr;
    I_KINETIC_Ix5z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5z_D2x_vrr;
    I_KINETIC_I6y_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6y_D2x_vrr;
    I_KINETIC_I5yz_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5yz_D2x_vrr;
    I_KINETIC_I4y2z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4y2z_D2x_vrr;
    I_KINETIC_I3y3z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3y3z_D2x_vrr;
    I_KINETIC_I2y4z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2y4z_D2x_vrr;
    I_KINETIC_Iy5z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Iy5z_D2x_vrr;
    I_KINETIC_I6z_D2x_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6z_D2x_vrr;
    I_KINETIC_I6x_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6x_Dxy_vrr;
    I_KINETIC_I5xy_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xy_Dxy_vrr;
    I_KINETIC_I5xz_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xz_Dxy_vrr;
    I_KINETIC_I4x2y_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2y_Dxy_vrr;
    I_KINETIC_I4xyz_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4xyz_Dxy_vrr;
    I_KINETIC_I4x2z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2z_Dxy_vrr;
    I_KINETIC_I3x3y_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3y_Dxy_vrr;
    I_KINETIC_I3x2yz_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x2yz_Dxy_vrr;
    I_KINETIC_I3xy2z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3xy2z_Dxy_vrr;
    I_KINETIC_I3x3z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3z_Dxy_vrr;
    I_KINETIC_I2x4y_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4y_Dxy_vrr;
    I_KINETIC_I2x3yz_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x3yz_Dxy_vrr;
    I_KINETIC_I2x2y2z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x2y2z_Dxy_vrr;
    I_KINETIC_I2xy3z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2xy3z_Dxy_vrr;
    I_KINETIC_I2x4z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4z_Dxy_vrr;
    I_KINETIC_Ix5y_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5y_Dxy_vrr;
    I_KINETIC_Ix4yz_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix4yz_Dxy_vrr;
    I_KINETIC_Ix3y2z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix3y2z_Dxy_vrr;
    I_KINETIC_Ix2y3z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix2y3z_Dxy_vrr;
    I_KINETIC_Ixy4z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ixy4z_Dxy_vrr;
    I_KINETIC_Ix5z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5z_Dxy_vrr;
    I_KINETIC_I6y_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6y_Dxy_vrr;
    I_KINETIC_I5yz_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5yz_Dxy_vrr;
    I_KINETIC_I4y2z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4y2z_Dxy_vrr;
    I_KINETIC_I3y3z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3y3z_Dxy_vrr;
    I_KINETIC_I2y4z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2y4z_Dxy_vrr;
    I_KINETIC_Iy5z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Iy5z_Dxy_vrr;
    I_KINETIC_I6z_Dxy_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6z_Dxy_vrr;
    I_KINETIC_I6x_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6x_Dxz_vrr;
    I_KINETIC_I5xy_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xy_Dxz_vrr;
    I_KINETIC_I5xz_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xz_Dxz_vrr;
    I_KINETIC_I4x2y_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2y_Dxz_vrr;
    I_KINETIC_I4xyz_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4xyz_Dxz_vrr;
    I_KINETIC_I4x2z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2z_Dxz_vrr;
    I_KINETIC_I3x3y_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3y_Dxz_vrr;
    I_KINETIC_I3x2yz_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x2yz_Dxz_vrr;
    I_KINETIC_I3xy2z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3xy2z_Dxz_vrr;
    I_KINETIC_I3x3z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3z_Dxz_vrr;
    I_KINETIC_I2x4y_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4y_Dxz_vrr;
    I_KINETIC_I2x3yz_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x3yz_Dxz_vrr;
    I_KINETIC_I2x2y2z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x2y2z_Dxz_vrr;
    I_KINETIC_I2xy3z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2xy3z_Dxz_vrr;
    I_KINETIC_I2x4z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4z_Dxz_vrr;
    I_KINETIC_Ix5y_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5y_Dxz_vrr;
    I_KINETIC_Ix4yz_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix4yz_Dxz_vrr;
    I_KINETIC_Ix3y2z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix3y2z_Dxz_vrr;
    I_KINETIC_Ix2y3z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix2y3z_Dxz_vrr;
    I_KINETIC_Ixy4z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ixy4z_Dxz_vrr;
    I_KINETIC_Ix5z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5z_Dxz_vrr;
    I_KINETIC_I6y_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6y_Dxz_vrr;
    I_KINETIC_I5yz_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5yz_Dxz_vrr;
    I_KINETIC_I4y2z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4y2z_Dxz_vrr;
    I_KINETIC_I3y3z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3y3z_Dxz_vrr;
    I_KINETIC_I2y4z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2y4z_Dxz_vrr;
    I_KINETIC_Iy5z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Iy5z_Dxz_vrr;
    I_KINETIC_I6z_Dxz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6z_Dxz_vrr;
    I_KINETIC_I6x_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6x_D2y_vrr;
    I_KINETIC_I5xy_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xy_D2y_vrr;
    I_KINETIC_I5xz_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xz_D2y_vrr;
    I_KINETIC_I4x2y_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2y_D2y_vrr;
    I_KINETIC_I4xyz_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4xyz_D2y_vrr;
    I_KINETIC_I4x2z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2z_D2y_vrr;
    I_KINETIC_I3x3y_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3y_D2y_vrr;
    I_KINETIC_I3x2yz_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x2yz_D2y_vrr;
    I_KINETIC_I3xy2z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3xy2z_D2y_vrr;
    I_KINETIC_I3x3z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3z_D2y_vrr;
    I_KINETIC_I2x4y_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4y_D2y_vrr;
    I_KINETIC_I2x3yz_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x3yz_D2y_vrr;
    I_KINETIC_I2x2y2z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x2y2z_D2y_vrr;
    I_KINETIC_I2xy3z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2xy3z_D2y_vrr;
    I_KINETIC_I2x4z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4z_D2y_vrr;
    I_KINETIC_Ix5y_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5y_D2y_vrr;
    I_KINETIC_Ix4yz_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix4yz_D2y_vrr;
    I_KINETIC_Ix3y2z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix3y2z_D2y_vrr;
    I_KINETIC_Ix2y3z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix2y3z_D2y_vrr;
    I_KINETIC_Ixy4z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ixy4z_D2y_vrr;
    I_KINETIC_Ix5z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5z_D2y_vrr;
    I_KINETIC_I6y_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6y_D2y_vrr;
    I_KINETIC_I5yz_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5yz_D2y_vrr;
    I_KINETIC_I4y2z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4y2z_D2y_vrr;
    I_KINETIC_I3y3z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3y3z_D2y_vrr;
    I_KINETIC_I2y4z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2y4z_D2y_vrr;
    I_KINETIC_Iy5z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Iy5z_D2y_vrr;
    I_KINETIC_I6z_D2y_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6z_D2y_vrr;
    I_KINETIC_I6x_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6x_Dyz_vrr;
    I_KINETIC_I5xy_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xy_Dyz_vrr;
    I_KINETIC_I5xz_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xz_Dyz_vrr;
    I_KINETIC_I4x2y_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2y_Dyz_vrr;
    I_KINETIC_I4xyz_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4xyz_Dyz_vrr;
    I_KINETIC_I4x2z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2z_Dyz_vrr;
    I_KINETIC_I3x3y_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3y_Dyz_vrr;
    I_KINETIC_I3x2yz_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x2yz_Dyz_vrr;
    I_KINETIC_I3xy2z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3xy2z_Dyz_vrr;
    I_KINETIC_I3x3z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3z_Dyz_vrr;
    I_KINETIC_I2x4y_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4y_Dyz_vrr;
    I_KINETIC_I2x3yz_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x3yz_Dyz_vrr;
    I_KINETIC_I2x2y2z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x2y2z_Dyz_vrr;
    I_KINETIC_I2xy3z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2xy3z_Dyz_vrr;
    I_KINETIC_I2x4z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4z_Dyz_vrr;
    I_KINETIC_Ix5y_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5y_Dyz_vrr;
    I_KINETIC_Ix4yz_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix4yz_Dyz_vrr;
    I_KINETIC_Ix3y2z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix3y2z_Dyz_vrr;
    I_KINETIC_Ix2y3z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix2y3z_Dyz_vrr;
    I_KINETIC_Ixy4z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ixy4z_Dyz_vrr;
    I_KINETIC_Ix5z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5z_Dyz_vrr;
    I_KINETIC_I6y_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6y_Dyz_vrr;
    I_KINETIC_I5yz_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5yz_Dyz_vrr;
    I_KINETIC_I4y2z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4y2z_Dyz_vrr;
    I_KINETIC_I3y3z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3y3z_Dyz_vrr;
    I_KINETIC_I2y4z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2y4z_Dyz_vrr;
    I_KINETIC_Iy5z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Iy5z_Dyz_vrr;
    I_KINETIC_I6z_Dyz_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6z_Dyz_vrr;
    I_KINETIC_I6x_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6x_D2z_vrr;
    I_KINETIC_I5xy_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xy_D2z_vrr;
    I_KINETIC_I5xz_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5xz_D2z_vrr;
    I_KINETIC_I4x2y_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2y_D2z_vrr;
    I_KINETIC_I4xyz_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4xyz_D2z_vrr;
    I_KINETIC_I4x2z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4x2z_D2z_vrr;
    I_KINETIC_I3x3y_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3y_D2z_vrr;
    I_KINETIC_I3x2yz_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x2yz_D2z_vrr;
    I_KINETIC_I3xy2z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3xy2z_D2z_vrr;
    I_KINETIC_I3x3z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3x3z_D2z_vrr;
    I_KINETIC_I2x4y_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4y_D2z_vrr;
    I_KINETIC_I2x3yz_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x3yz_D2z_vrr;
    I_KINETIC_I2x2y2z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x2y2z_D2z_vrr;
    I_KINETIC_I2xy3z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2xy3z_D2z_vrr;
    I_KINETIC_I2x4z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2x4z_D2z_vrr;
    I_KINETIC_Ix5y_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5y_D2z_vrr;
    I_KINETIC_Ix4yz_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix4yz_D2z_vrr;
    I_KINETIC_Ix3y2z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix3y2z_D2z_vrr;
    I_KINETIC_Ix2y3z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix2y3z_D2z_vrr;
    I_KINETIC_Ixy4z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ixy4z_D2z_vrr;
    I_KINETIC_Ix5z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Ix5z_D2z_vrr;
    I_KINETIC_I6y_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6y_D2z_vrr;
    I_KINETIC_I5yz_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I5yz_D2z_vrr;
    I_KINETIC_I4y2z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I4y2z_D2z_vrr;
    I_KINETIC_I3y3z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I3y3z_D2z_vrr;
    I_KINETIC_I2y4z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I2y4z_D2z_vrr;
    I_KINETIC_Iy5z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_Iy5z_D2z_vrr;
    I_KINETIC_I6z_D2z_a += SQ_KINETIC_I_D_a_coefs*I_KINETIC_I6z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_D
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_G4x_D2x += I_KINETIC_G4x_D2x_vrr;
    I_KINETIC_G3xy_D2x += I_KINETIC_G3xy_D2x_vrr;
    I_KINETIC_G3xz_D2x += I_KINETIC_G3xz_D2x_vrr;
    I_KINETIC_G2x2y_D2x += I_KINETIC_G2x2y_D2x_vrr;
    I_KINETIC_G2xyz_D2x += I_KINETIC_G2xyz_D2x_vrr;
    I_KINETIC_G2x2z_D2x += I_KINETIC_G2x2z_D2x_vrr;
    I_KINETIC_Gx3y_D2x += I_KINETIC_Gx3y_D2x_vrr;
    I_KINETIC_Gx2yz_D2x += I_KINETIC_Gx2yz_D2x_vrr;
    I_KINETIC_Gxy2z_D2x += I_KINETIC_Gxy2z_D2x_vrr;
    I_KINETIC_Gx3z_D2x += I_KINETIC_Gx3z_D2x_vrr;
    I_KINETIC_G4y_D2x += I_KINETIC_G4y_D2x_vrr;
    I_KINETIC_G3yz_D2x += I_KINETIC_G3yz_D2x_vrr;
    I_KINETIC_G2y2z_D2x += I_KINETIC_G2y2z_D2x_vrr;
    I_KINETIC_Gy3z_D2x += I_KINETIC_Gy3z_D2x_vrr;
    I_KINETIC_G4z_D2x += I_KINETIC_G4z_D2x_vrr;
    I_KINETIC_G4x_Dxy += I_KINETIC_G4x_Dxy_vrr;
    I_KINETIC_G3xy_Dxy += I_KINETIC_G3xy_Dxy_vrr;
    I_KINETIC_G3xz_Dxy += I_KINETIC_G3xz_Dxy_vrr;
    I_KINETIC_G2x2y_Dxy += I_KINETIC_G2x2y_Dxy_vrr;
    I_KINETIC_G2xyz_Dxy += I_KINETIC_G2xyz_Dxy_vrr;
    I_KINETIC_G2x2z_Dxy += I_KINETIC_G2x2z_Dxy_vrr;
    I_KINETIC_Gx3y_Dxy += I_KINETIC_Gx3y_Dxy_vrr;
    I_KINETIC_Gx2yz_Dxy += I_KINETIC_Gx2yz_Dxy_vrr;
    I_KINETIC_Gxy2z_Dxy += I_KINETIC_Gxy2z_Dxy_vrr;
    I_KINETIC_Gx3z_Dxy += I_KINETIC_Gx3z_Dxy_vrr;
    I_KINETIC_G4y_Dxy += I_KINETIC_G4y_Dxy_vrr;
    I_KINETIC_G3yz_Dxy += I_KINETIC_G3yz_Dxy_vrr;
    I_KINETIC_G2y2z_Dxy += I_KINETIC_G2y2z_Dxy_vrr;
    I_KINETIC_Gy3z_Dxy += I_KINETIC_Gy3z_Dxy_vrr;
    I_KINETIC_G4z_Dxy += I_KINETIC_G4z_Dxy_vrr;
    I_KINETIC_G4x_Dxz += I_KINETIC_G4x_Dxz_vrr;
    I_KINETIC_G3xy_Dxz += I_KINETIC_G3xy_Dxz_vrr;
    I_KINETIC_G3xz_Dxz += I_KINETIC_G3xz_Dxz_vrr;
    I_KINETIC_G2x2y_Dxz += I_KINETIC_G2x2y_Dxz_vrr;
    I_KINETIC_G2xyz_Dxz += I_KINETIC_G2xyz_Dxz_vrr;
    I_KINETIC_G2x2z_Dxz += I_KINETIC_G2x2z_Dxz_vrr;
    I_KINETIC_Gx3y_Dxz += I_KINETIC_Gx3y_Dxz_vrr;
    I_KINETIC_Gx2yz_Dxz += I_KINETIC_Gx2yz_Dxz_vrr;
    I_KINETIC_Gxy2z_Dxz += I_KINETIC_Gxy2z_Dxz_vrr;
    I_KINETIC_Gx3z_Dxz += I_KINETIC_Gx3z_Dxz_vrr;
    I_KINETIC_G4y_Dxz += I_KINETIC_G4y_Dxz_vrr;
    I_KINETIC_G3yz_Dxz += I_KINETIC_G3yz_Dxz_vrr;
    I_KINETIC_G2y2z_Dxz += I_KINETIC_G2y2z_Dxz_vrr;
    I_KINETIC_Gy3z_Dxz += I_KINETIC_Gy3z_Dxz_vrr;
    I_KINETIC_G4z_Dxz += I_KINETIC_G4z_Dxz_vrr;
    I_KINETIC_G4x_D2y += I_KINETIC_G4x_D2y_vrr;
    I_KINETIC_G3xy_D2y += I_KINETIC_G3xy_D2y_vrr;
    I_KINETIC_G3xz_D2y += I_KINETIC_G3xz_D2y_vrr;
    I_KINETIC_G2x2y_D2y += I_KINETIC_G2x2y_D2y_vrr;
    I_KINETIC_G2xyz_D2y += I_KINETIC_G2xyz_D2y_vrr;
    I_KINETIC_G2x2z_D2y += I_KINETIC_G2x2z_D2y_vrr;
    I_KINETIC_Gx3y_D2y += I_KINETIC_Gx3y_D2y_vrr;
    I_KINETIC_Gx2yz_D2y += I_KINETIC_Gx2yz_D2y_vrr;
    I_KINETIC_Gxy2z_D2y += I_KINETIC_Gxy2z_D2y_vrr;
    I_KINETIC_Gx3z_D2y += I_KINETIC_Gx3z_D2y_vrr;
    I_KINETIC_G4y_D2y += I_KINETIC_G4y_D2y_vrr;
    I_KINETIC_G3yz_D2y += I_KINETIC_G3yz_D2y_vrr;
    I_KINETIC_G2y2z_D2y += I_KINETIC_G2y2z_D2y_vrr;
    I_KINETIC_Gy3z_D2y += I_KINETIC_Gy3z_D2y_vrr;
    I_KINETIC_G4z_D2y += I_KINETIC_G4z_D2y_vrr;
    I_KINETIC_G4x_Dyz += I_KINETIC_G4x_Dyz_vrr;
    I_KINETIC_G3xy_Dyz += I_KINETIC_G3xy_Dyz_vrr;
    I_KINETIC_G3xz_Dyz += I_KINETIC_G3xz_Dyz_vrr;
    I_KINETIC_G2x2y_Dyz += I_KINETIC_G2x2y_Dyz_vrr;
    I_KINETIC_G2xyz_Dyz += I_KINETIC_G2xyz_Dyz_vrr;
    I_KINETIC_G2x2z_Dyz += I_KINETIC_G2x2z_Dyz_vrr;
    I_KINETIC_Gx3y_Dyz += I_KINETIC_Gx3y_Dyz_vrr;
    I_KINETIC_Gx2yz_Dyz += I_KINETIC_Gx2yz_Dyz_vrr;
    I_KINETIC_Gxy2z_Dyz += I_KINETIC_Gxy2z_Dyz_vrr;
    I_KINETIC_Gx3z_Dyz += I_KINETIC_Gx3z_Dyz_vrr;
    I_KINETIC_G4y_Dyz += I_KINETIC_G4y_Dyz_vrr;
    I_KINETIC_G3yz_Dyz += I_KINETIC_G3yz_Dyz_vrr;
    I_KINETIC_G2y2z_Dyz += I_KINETIC_G2y2z_Dyz_vrr;
    I_KINETIC_Gy3z_Dyz += I_KINETIC_Gy3z_Dyz_vrr;
    I_KINETIC_G4z_Dyz += I_KINETIC_G4z_Dyz_vrr;
    I_KINETIC_G4x_D2z += I_KINETIC_G4x_D2z_vrr;
    I_KINETIC_G3xy_D2z += I_KINETIC_G3xy_D2z_vrr;
    I_KINETIC_G3xz_D2z += I_KINETIC_G3xz_D2z_vrr;
    I_KINETIC_G2x2y_D2z += I_KINETIC_G2x2y_D2z_vrr;
    I_KINETIC_G2xyz_D2z += I_KINETIC_G2xyz_D2z_vrr;
    I_KINETIC_G2x2z_D2z += I_KINETIC_G2x2z_D2z_vrr;
    I_KINETIC_Gx3y_D2z += I_KINETIC_Gx3y_D2z_vrr;
    I_KINETIC_Gx2yz_D2z += I_KINETIC_Gx2yz_D2z_vrr;
    I_KINETIC_Gxy2z_D2z += I_KINETIC_Gxy2z_D2z_vrr;
    I_KINETIC_Gx3z_D2z += I_KINETIC_Gx3z_D2z_vrr;
    I_KINETIC_G4y_D2z += I_KINETIC_G4y_D2z_vrr;
    I_KINETIC_G3yz_D2z += I_KINETIC_G3yz_D2z_vrr;
    I_KINETIC_G2y2z_D2z += I_KINETIC_G2y2z_D2z_vrr;
    I_KINETIC_Gy3z_D2z += I_KINETIC_Gy3z_D2z_vrr;
    I_KINETIC_G4z_D2z += I_KINETIC_G4z_D2z_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_D_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_D_a
   * RHS shell quartet name: SQ_KINETIC_G_D
   ************************************************************/
  abcd[0] = 2.0E0*I_KINETIC_I6x_D2x_a-5*I_KINETIC_G4x_D2x;
  abcd[1] = 2.0E0*I_KINETIC_I5xy_D2x_a-4*I_KINETIC_G3xy_D2x;
  abcd[2] = 2.0E0*I_KINETIC_I5xz_D2x_a-4*I_KINETIC_G3xz_D2x;
  abcd[3] = 2.0E0*I_KINETIC_I4x2y_D2x_a-3*I_KINETIC_G2x2y_D2x;
  abcd[4] = 2.0E0*I_KINETIC_I4xyz_D2x_a-3*I_KINETIC_G2xyz_D2x;
  abcd[5] = 2.0E0*I_KINETIC_I4x2z_D2x_a-3*I_KINETIC_G2x2z_D2x;
  abcd[6] = 2.0E0*I_KINETIC_I3x3y_D2x_a-2*I_KINETIC_Gx3y_D2x;
  abcd[7] = 2.0E0*I_KINETIC_I3x2yz_D2x_a-2*I_KINETIC_Gx2yz_D2x;
  abcd[8] = 2.0E0*I_KINETIC_I3xy2z_D2x_a-2*I_KINETIC_Gxy2z_D2x;
  abcd[9] = 2.0E0*I_KINETIC_I3x3z_D2x_a-2*I_KINETIC_Gx3z_D2x;
  abcd[10] = 2.0E0*I_KINETIC_I2x4y_D2x_a-1*I_KINETIC_G4y_D2x;
  abcd[11] = 2.0E0*I_KINETIC_I2x3yz_D2x_a-1*I_KINETIC_G3yz_D2x;
  abcd[12] = 2.0E0*I_KINETIC_I2x2y2z_D2x_a-1*I_KINETIC_G2y2z_D2x;
  abcd[13] = 2.0E0*I_KINETIC_I2xy3z_D2x_a-1*I_KINETIC_Gy3z_D2x;
  abcd[14] = 2.0E0*I_KINETIC_I2x4z_D2x_a-1*I_KINETIC_G4z_D2x;
  abcd[15] = 2.0E0*I_KINETIC_Ix5y_D2x_a;
  abcd[16] = 2.0E0*I_KINETIC_Ix4yz_D2x_a;
  abcd[17] = 2.0E0*I_KINETIC_Ix3y2z_D2x_a;
  abcd[18] = 2.0E0*I_KINETIC_Ix2y3z_D2x_a;
  abcd[19] = 2.0E0*I_KINETIC_Ixy4z_D2x_a;
  abcd[20] = 2.0E0*I_KINETIC_Ix5z_D2x_a;
  abcd[21] = 2.0E0*I_KINETIC_I6x_Dxy_a-5*I_KINETIC_G4x_Dxy;
  abcd[22] = 2.0E0*I_KINETIC_I5xy_Dxy_a-4*I_KINETIC_G3xy_Dxy;
  abcd[23] = 2.0E0*I_KINETIC_I5xz_Dxy_a-4*I_KINETIC_G3xz_Dxy;
  abcd[24] = 2.0E0*I_KINETIC_I4x2y_Dxy_a-3*I_KINETIC_G2x2y_Dxy;
  abcd[25] = 2.0E0*I_KINETIC_I4xyz_Dxy_a-3*I_KINETIC_G2xyz_Dxy;
  abcd[26] = 2.0E0*I_KINETIC_I4x2z_Dxy_a-3*I_KINETIC_G2x2z_Dxy;
  abcd[27] = 2.0E0*I_KINETIC_I3x3y_Dxy_a-2*I_KINETIC_Gx3y_Dxy;
  abcd[28] = 2.0E0*I_KINETIC_I3x2yz_Dxy_a-2*I_KINETIC_Gx2yz_Dxy;
  abcd[29] = 2.0E0*I_KINETIC_I3xy2z_Dxy_a-2*I_KINETIC_Gxy2z_Dxy;
  abcd[30] = 2.0E0*I_KINETIC_I3x3z_Dxy_a-2*I_KINETIC_Gx3z_Dxy;
  abcd[31] = 2.0E0*I_KINETIC_I2x4y_Dxy_a-1*I_KINETIC_G4y_Dxy;
  abcd[32] = 2.0E0*I_KINETIC_I2x3yz_Dxy_a-1*I_KINETIC_G3yz_Dxy;
  abcd[33] = 2.0E0*I_KINETIC_I2x2y2z_Dxy_a-1*I_KINETIC_G2y2z_Dxy;
  abcd[34] = 2.0E0*I_KINETIC_I2xy3z_Dxy_a-1*I_KINETIC_Gy3z_Dxy;
  abcd[35] = 2.0E0*I_KINETIC_I2x4z_Dxy_a-1*I_KINETIC_G4z_Dxy;
  abcd[36] = 2.0E0*I_KINETIC_Ix5y_Dxy_a;
  abcd[37] = 2.0E0*I_KINETIC_Ix4yz_Dxy_a;
  abcd[38] = 2.0E0*I_KINETIC_Ix3y2z_Dxy_a;
  abcd[39] = 2.0E0*I_KINETIC_Ix2y3z_Dxy_a;
  abcd[40] = 2.0E0*I_KINETIC_Ixy4z_Dxy_a;
  abcd[41] = 2.0E0*I_KINETIC_Ix5z_Dxy_a;
  abcd[42] = 2.0E0*I_KINETIC_I6x_Dxz_a-5*I_KINETIC_G4x_Dxz;
  abcd[43] = 2.0E0*I_KINETIC_I5xy_Dxz_a-4*I_KINETIC_G3xy_Dxz;
  abcd[44] = 2.0E0*I_KINETIC_I5xz_Dxz_a-4*I_KINETIC_G3xz_Dxz;
  abcd[45] = 2.0E0*I_KINETIC_I4x2y_Dxz_a-3*I_KINETIC_G2x2y_Dxz;
  abcd[46] = 2.0E0*I_KINETIC_I4xyz_Dxz_a-3*I_KINETIC_G2xyz_Dxz;
  abcd[47] = 2.0E0*I_KINETIC_I4x2z_Dxz_a-3*I_KINETIC_G2x2z_Dxz;
  abcd[48] = 2.0E0*I_KINETIC_I3x3y_Dxz_a-2*I_KINETIC_Gx3y_Dxz;
  abcd[49] = 2.0E0*I_KINETIC_I3x2yz_Dxz_a-2*I_KINETIC_Gx2yz_Dxz;
  abcd[50] = 2.0E0*I_KINETIC_I3xy2z_Dxz_a-2*I_KINETIC_Gxy2z_Dxz;
  abcd[51] = 2.0E0*I_KINETIC_I3x3z_Dxz_a-2*I_KINETIC_Gx3z_Dxz;
  abcd[52] = 2.0E0*I_KINETIC_I2x4y_Dxz_a-1*I_KINETIC_G4y_Dxz;
  abcd[53] = 2.0E0*I_KINETIC_I2x3yz_Dxz_a-1*I_KINETIC_G3yz_Dxz;
  abcd[54] = 2.0E0*I_KINETIC_I2x2y2z_Dxz_a-1*I_KINETIC_G2y2z_Dxz;
  abcd[55] = 2.0E0*I_KINETIC_I2xy3z_Dxz_a-1*I_KINETIC_Gy3z_Dxz;
  abcd[56] = 2.0E0*I_KINETIC_I2x4z_Dxz_a-1*I_KINETIC_G4z_Dxz;
  abcd[57] = 2.0E0*I_KINETIC_Ix5y_Dxz_a;
  abcd[58] = 2.0E0*I_KINETIC_Ix4yz_Dxz_a;
  abcd[59] = 2.0E0*I_KINETIC_Ix3y2z_Dxz_a;
  abcd[60] = 2.0E0*I_KINETIC_Ix2y3z_Dxz_a;
  abcd[61] = 2.0E0*I_KINETIC_Ixy4z_Dxz_a;
  abcd[62] = 2.0E0*I_KINETIC_Ix5z_Dxz_a;
  abcd[63] = 2.0E0*I_KINETIC_I6x_D2y_a-5*I_KINETIC_G4x_D2y;
  abcd[64] = 2.0E0*I_KINETIC_I5xy_D2y_a-4*I_KINETIC_G3xy_D2y;
  abcd[65] = 2.0E0*I_KINETIC_I5xz_D2y_a-4*I_KINETIC_G3xz_D2y;
  abcd[66] = 2.0E0*I_KINETIC_I4x2y_D2y_a-3*I_KINETIC_G2x2y_D2y;
  abcd[67] = 2.0E0*I_KINETIC_I4xyz_D2y_a-3*I_KINETIC_G2xyz_D2y;
  abcd[68] = 2.0E0*I_KINETIC_I4x2z_D2y_a-3*I_KINETIC_G2x2z_D2y;
  abcd[69] = 2.0E0*I_KINETIC_I3x3y_D2y_a-2*I_KINETIC_Gx3y_D2y;
  abcd[70] = 2.0E0*I_KINETIC_I3x2yz_D2y_a-2*I_KINETIC_Gx2yz_D2y;
  abcd[71] = 2.0E0*I_KINETIC_I3xy2z_D2y_a-2*I_KINETIC_Gxy2z_D2y;
  abcd[72] = 2.0E0*I_KINETIC_I3x3z_D2y_a-2*I_KINETIC_Gx3z_D2y;
  abcd[73] = 2.0E0*I_KINETIC_I2x4y_D2y_a-1*I_KINETIC_G4y_D2y;
  abcd[74] = 2.0E0*I_KINETIC_I2x3yz_D2y_a-1*I_KINETIC_G3yz_D2y;
  abcd[75] = 2.0E0*I_KINETIC_I2x2y2z_D2y_a-1*I_KINETIC_G2y2z_D2y;
  abcd[76] = 2.0E0*I_KINETIC_I2xy3z_D2y_a-1*I_KINETIC_Gy3z_D2y;
  abcd[77] = 2.0E0*I_KINETIC_I2x4z_D2y_a-1*I_KINETIC_G4z_D2y;
  abcd[78] = 2.0E0*I_KINETIC_Ix5y_D2y_a;
  abcd[79] = 2.0E0*I_KINETIC_Ix4yz_D2y_a;
  abcd[80] = 2.0E0*I_KINETIC_Ix3y2z_D2y_a;
  abcd[81] = 2.0E0*I_KINETIC_Ix2y3z_D2y_a;
  abcd[82] = 2.0E0*I_KINETIC_Ixy4z_D2y_a;
  abcd[83] = 2.0E0*I_KINETIC_Ix5z_D2y_a;
  abcd[84] = 2.0E0*I_KINETIC_I6x_Dyz_a-5*I_KINETIC_G4x_Dyz;
  abcd[85] = 2.0E0*I_KINETIC_I5xy_Dyz_a-4*I_KINETIC_G3xy_Dyz;
  abcd[86] = 2.0E0*I_KINETIC_I5xz_Dyz_a-4*I_KINETIC_G3xz_Dyz;
  abcd[87] = 2.0E0*I_KINETIC_I4x2y_Dyz_a-3*I_KINETIC_G2x2y_Dyz;
  abcd[88] = 2.0E0*I_KINETIC_I4xyz_Dyz_a-3*I_KINETIC_G2xyz_Dyz;
  abcd[89] = 2.0E0*I_KINETIC_I4x2z_Dyz_a-3*I_KINETIC_G2x2z_Dyz;
  abcd[90] = 2.0E0*I_KINETIC_I3x3y_Dyz_a-2*I_KINETIC_Gx3y_Dyz;
  abcd[91] = 2.0E0*I_KINETIC_I3x2yz_Dyz_a-2*I_KINETIC_Gx2yz_Dyz;
  abcd[92] = 2.0E0*I_KINETIC_I3xy2z_Dyz_a-2*I_KINETIC_Gxy2z_Dyz;
  abcd[93] = 2.0E0*I_KINETIC_I3x3z_Dyz_a-2*I_KINETIC_Gx3z_Dyz;
  abcd[94] = 2.0E0*I_KINETIC_I2x4y_Dyz_a-1*I_KINETIC_G4y_Dyz;
  abcd[95] = 2.0E0*I_KINETIC_I2x3yz_Dyz_a-1*I_KINETIC_G3yz_Dyz;
  abcd[96] = 2.0E0*I_KINETIC_I2x2y2z_Dyz_a-1*I_KINETIC_G2y2z_Dyz;
  abcd[97] = 2.0E0*I_KINETIC_I2xy3z_Dyz_a-1*I_KINETIC_Gy3z_Dyz;
  abcd[98] = 2.0E0*I_KINETIC_I2x4z_Dyz_a-1*I_KINETIC_G4z_Dyz;
  abcd[99] = 2.0E0*I_KINETIC_Ix5y_Dyz_a;
  abcd[100] = 2.0E0*I_KINETIC_Ix4yz_Dyz_a;
  abcd[101] = 2.0E0*I_KINETIC_Ix3y2z_Dyz_a;
  abcd[102] = 2.0E0*I_KINETIC_Ix2y3z_Dyz_a;
  abcd[103] = 2.0E0*I_KINETIC_Ixy4z_Dyz_a;
  abcd[104] = 2.0E0*I_KINETIC_Ix5z_Dyz_a;
  abcd[105] = 2.0E0*I_KINETIC_I6x_D2z_a-5*I_KINETIC_G4x_D2z;
  abcd[106] = 2.0E0*I_KINETIC_I5xy_D2z_a-4*I_KINETIC_G3xy_D2z;
  abcd[107] = 2.0E0*I_KINETIC_I5xz_D2z_a-4*I_KINETIC_G3xz_D2z;
  abcd[108] = 2.0E0*I_KINETIC_I4x2y_D2z_a-3*I_KINETIC_G2x2y_D2z;
  abcd[109] = 2.0E0*I_KINETIC_I4xyz_D2z_a-3*I_KINETIC_G2xyz_D2z;
  abcd[110] = 2.0E0*I_KINETIC_I4x2z_D2z_a-3*I_KINETIC_G2x2z_D2z;
  abcd[111] = 2.0E0*I_KINETIC_I3x3y_D2z_a-2*I_KINETIC_Gx3y_D2z;
  abcd[112] = 2.0E0*I_KINETIC_I3x2yz_D2z_a-2*I_KINETIC_Gx2yz_D2z;
  abcd[113] = 2.0E0*I_KINETIC_I3xy2z_D2z_a-2*I_KINETIC_Gxy2z_D2z;
  abcd[114] = 2.0E0*I_KINETIC_I3x3z_D2z_a-2*I_KINETIC_Gx3z_D2z;
  abcd[115] = 2.0E0*I_KINETIC_I2x4y_D2z_a-1*I_KINETIC_G4y_D2z;
  abcd[116] = 2.0E0*I_KINETIC_I2x3yz_D2z_a-1*I_KINETIC_G3yz_D2z;
  abcd[117] = 2.0E0*I_KINETIC_I2x2y2z_D2z_a-1*I_KINETIC_G2y2z_D2z;
  abcd[118] = 2.0E0*I_KINETIC_I2xy3z_D2z_a-1*I_KINETIC_Gy3z_D2z;
  abcd[119] = 2.0E0*I_KINETIC_I2x4z_D2z_a-1*I_KINETIC_G4z_D2z;
  abcd[120] = 2.0E0*I_KINETIC_Ix5y_D2z_a;
  abcd[121] = 2.0E0*I_KINETIC_Ix4yz_D2z_a;
  abcd[122] = 2.0E0*I_KINETIC_Ix3y2z_D2z_a;
  abcd[123] = 2.0E0*I_KINETIC_Ix2y3z_D2z_a;
  abcd[124] = 2.0E0*I_KINETIC_Ixy4z_D2z_a;
  abcd[125] = 2.0E0*I_KINETIC_Ix5z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_D_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_D_a
   * RHS shell quartet name: SQ_KINETIC_G_D
   ************************************************************/
  abcd[126] = 2.0E0*I_KINETIC_I5xy_D2x_a;
  abcd[127] = 2.0E0*I_KINETIC_I4x2y_D2x_a-1*I_KINETIC_G4x_D2x;
  abcd[128] = 2.0E0*I_KINETIC_I4xyz_D2x_a;
  abcd[129] = 2.0E0*I_KINETIC_I3x3y_D2x_a-2*I_KINETIC_G3xy_D2x;
  abcd[130] = 2.0E0*I_KINETIC_I3x2yz_D2x_a-1*I_KINETIC_G3xz_D2x;
  abcd[131] = 2.0E0*I_KINETIC_I3xy2z_D2x_a;
  abcd[132] = 2.0E0*I_KINETIC_I2x4y_D2x_a-3*I_KINETIC_G2x2y_D2x;
  abcd[133] = 2.0E0*I_KINETIC_I2x3yz_D2x_a-2*I_KINETIC_G2xyz_D2x;
  abcd[134] = 2.0E0*I_KINETIC_I2x2y2z_D2x_a-1*I_KINETIC_G2x2z_D2x;
  abcd[135] = 2.0E0*I_KINETIC_I2xy3z_D2x_a;
  abcd[136] = 2.0E0*I_KINETIC_Ix5y_D2x_a-4*I_KINETIC_Gx3y_D2x;
  abcd[137] = 2.0E0*I_KINETIC_Ix4yz_D2x_a-3*I_KINETIC_Gx2yz_D2x;
  abcd[138] = 2.0E0*I_KINETIC_Ix3y2z_D2x_a-2*I_KINETIC_Gxy2z_D2x;
  abcd[139] = 2.0E0*I_KINETIC_Ix2y3z_D2x_a-1*I_KINETIC_Gx3z_D2x;
  abcd[140] = 2.0E0*I_KINETIC_Ixy4z_D2x_a;
  abcd[141] = 2.0E0*I_KINETIC_I6y_D2x_a-5*I_KINETIC_G4y_D2x;
  abcd[142] = 2.0E0*I_KINETIC_I5yz_D2x_a-4*I_KINETIC_G3yz_D2x;
  abcd[143] = 2.0E0*I_KINETIC_I4y2z_D2x_a-3*I_KINETIC_G2y2z_D2x;
  abcd[144] = 2.0E0*I_KINETIC_I3y3z_D2x_a-2*I_KINETIC_Gy3z_D2x;
  abcd[145] = 2.0E0*I_KINETIC_I2y4z_D2x_a-1*I_KINETIC_G4z_D2x;
  abcd[146] = 2.0E0*I_KINETIC_Iy5z_D2x_a;
  abcd[147] = 2.0E0*I_KINETIC_I5xy_Dxy_a;
  abcd[148] = 2.0E0*I_KINETIC_I4x2y_Dxy_a-1*I_KINETIC_G4x_Dxy;
  abcd[149] = 2.0E0*I_KINETIC_I4xyz_Dxy_a;
  abcd[150] = 2.0E0*I_KINETIC_I3x3y_Dxy_a-2*I_KINETIC_G3xy_Dxy;
  abcd[151] = 2.0E0*I_KINETIC_I3x2yz_Dxy_a-1*I_KINETIC_G3xz_Dxy;
  abcd[152] = 2.0E0*I_KINETIC_I3xy2z_Dxy_a;
  abcd[153] = 2.0E0*I_KINETIC_I2x4y_Dxy_a-3*I_KINETIC_G2x2y_Dxy;
  abcd[154] = 2.0E0*I_KINETIC_I2x3yz_Dxy_a-2*I_KINETIC_G2xyz_Dxy;
  abcd[155] = 2.0E0*I_KINETIC_I2x2y2z_Dxy_a-1*I_KINETIC_G2x2z_Dxy;
  abcd[156] = 2.0E0*I_KINETIC_I2xy3z_Dxy_a;
  abcd[157] = 2.0E0*I_KINETIC_Ix5y_Dxy_a-4*I_KINETIC_Gx3y_Dxy;
  abcd[158] = 2.0E0*I_KINETIC_Ix4yz_Dxy_a-3*I_KINETIC_Gx2yz_Dxy;
  abcd[159] = 2.0E0*I_KINETIC_Ix3y2z_Dxy_a-2*I_KINETIC_Gxy2z_Dxy;
  abcd[160] = 2.0E0*I_KINETIC_Ix2y3z_Dxy_a-1*I_KINETIC_Gx3z_Dxy;
  abcd[161] = 2.0E0*I_KINETIC_Ixy4z_Dxy_a;
  abcd[162] = 2.0E0*I_KINETIC_I6y_Dxy_a-5*I_KINETIC_G4y_Dxy;
  abcd[163] = 2.0E0*I_KINETIC_I5yz_Dxy_a-4*I_KINETIC_G3yz_Dxy;
  abcd[164] = 2.0E0*I_KINETIC_I4y2z_Dxy_a-3*I_KINETIC_G2y2z_Dxy;
  abcd[165] = 2.0E0*I_KINETIC_I3y3z_Dxy_a-2*I_KINETIC_Gy3z_Dxy;
  abcd[166] = 2.0E0*I_KINETIC_I2y4z_Dxy_a-1*I_KINETIC_G4z_Dxy;
  abcd[167] = 2.0E0*I_KINETIC_Iy5z_Dxy_a;
  abcd[168] = 2.0E0*I_KINETIC_I5xy_Dxz_a;
  abcd[169] = 2.0E0*I_KINETIC_I4x2y_Dxz_a-1*I_KINETIC_G4x_Dxz;
  abcd[170] = 2.0E0*I_KINETIC_I4xyz_Dxz_a;
  abcd[171] = 2.0E0*I_KINETIC_I3x3y_Dxz_a-2*I_KINETIC_G3xy_Dxz;
  abcd[172] = 2.0E0*I_KINETIC_I3x2yz_Dxz_a-1*I_KINETIC_G3xz_Dxz;
  abcd[173] = 2.0E0*I_KINETIC_I3xy2z_Dxz_a;
  abcd[174] = 2.0E0*I_KINETIC_I2x4y_Dxz_a-3*I_KINETIC_G2x2y_Dxz;
  abcd[175] = 2.0E0*I_KINETIC_I2x3yz_Dxz_a-2*I_KINETIC_G2xyz_Dxz;
  abcd[176] = 2.0E0*I_KINETIC_I2x2y2z_Dxz_a-1*I_KINETIC_G2x2z_Dxz;
  abcd[177] = 2.0E0*I_KINETIC_I2xy3z_Dxz_a;
  abcd[178] = 2.0E0*I_KINETIC_Ix5y_Dxz_a-4*I_KINETIC_Gx3y_Dxz;
  abcd[179] = 2.0E0*I_KINETIC_Ix4yz_Dxz_a-3*I_KINETIC_Gx2yz_Dxz;
  abcd[180] = 2.0E0*I_KINETIC_Ix3y2z_Dxz_a-2*I_KINETIC_Gxy2z_Dxz;
  abcd[181] = 2.0E0*I_KINETIC_Ix2y3z_Dxz_a-1*I_KINETIC_Gx3z_Dxz;
  abcd[182] = 2.0E0*I_KINETIC_Ixy4z_Dxz_a;
  abcd[183] = 2.0E0*I_KINETIC_I6y_Dxz_a-5*I_KINETIC_G4y_Dxz;
  abcd[184] = 2.0E0*I_KINETIC_I5yz_Dxz_a-4*I_KINETIC_G3yz_Dxz;
  abcd[185] = 2.0E0*I_KINETIC_I4y2z_Dxz_a-3*I_KINETIC_G2y2z_Dxz;
  abcd[186] = 2.0E0*I_KINETIC_I3y3z_Dxz_a-2*I_KINETIC_Gy3z_Dxz;
  abcd[187] = 2.0E0*I_KINETIC_I2y4z_Dxz_a-1*I_KINETIC_G4z_Dxz;
  abcd[188] = 2.0E0*I_KINETIC_Iy5z_Dxz_a;
  abcd[189] = 2.0E0*I_KINETIC_I5xy_D2y_a;
  abcd[190] = 2.0E0*I_KINETIC_I4x2y_D2y_a-1*I_KINETIC_G4x_D2y;
  abcd[191] = 2.0E0*I_KINETIC_I4xyz_D2y_a;
  abcd[192] = 2.0E0*I_KINETIC_I3x3y_D2y_a-2*I_KINETIC_G3xy_D2y;
  abcd[193] = 2.0E0*I_KINETIC_I3x2yz_D2y_a-1*I_KINETIC_G3xz_D2y;
  abcd[194] = 2.0E0*I_KINETIC_I3xy2z_D2y_a;
  abcd[195] = 2.0E0*I_KINETIC_I2x4y_D2y_a-3*I_KINETIC_G2x2y_D2y;
  abcd[196] = 2.0E0*I_KINETIC_I2x3yz_D2y_a-2*I_KINETIC_G2xyz_D2y;
  abcd[197] = 2.0E0*I_KINETIC_I2x2y2z_D2y_a-1*I_KINETIC_G2x2z_D2y;
  abcd[198] = 2.0E0*I_KINETIC_I2xy3z_D2y_a;
  abcd[199] = 2.0E0*I_KINETIC_Ix5y_D2y_a-4*I_KINETIC_Gx3y_D2y;
  abcd[200] = 2.0E0*I_KINETIC_Ix4yz_D2y_a-3*I_KINETIC_Gx2yz_D2y;
  abcd[201] = 2.0E0*I_KINETIC_Ix3y2z_D2y_a-2*I_KINETIC_Gxy2z_D2y;
  abcd[202] = 2.0E0*I_KINETIC_Ix2y3z_D2y_a-1*I_KINETIC_Gx3z_D2y;
  abcd[203] = 2.0E0*I_KINETIC_Ixy4z_D2y_a;
  abcd[204] = 2.0E0*I_KINETIC_I6y_D2y_a-5*I_KINETIC_G4y_D2y;
  abcd[205] = 2.0E0*I_KINETIC_I5yz_D2y_a-4*I_KINETIC_G3yz_D2y;
  abcd[206] = 2.0E0*I_KINETIC_I4y2z_D2y_a-3*I_KINETIC_G2y2z_D2y;
  abcd[207] = 2.0E0*I_KINETIC_I3y3z_D2y_a-2*I_KINETIC_Gy3z_D2y;
  abcd[208] = 2.0E0*I_KINETIC_I2y4z_D2y_a-1*I_KINETIC_G4z_D2y;
  abcd[209] = 2.0E0*I_KINETIC_Iy5z_D2y_a;
  abcd[210] = 2.0E0*I_KINETIC_I5xy_Dyz_a;
  abcd[211] = 2.0E0*I_KINETIC_I4x2y_Dyz_a-1*I_KINETIC_G4x_Dyz;
  abcd[212] = 2.0E0*I_KINETIC_I4xyz_Dyz_a;
  abcd[213] = 2.0E0*I_KINETIC_I3x3y_Dyz_a-2*I_KINETIC_G3xy_Dyz;
  abcd[214] = 2.0E0*I_KINETIC_I3x2yz_Dyz_a-1*I_KINETIC_G3xz_Dyz;
  abcd[215] = 2.0E0*I_KINETIC_I3xy2z_Dyz_a;
  abcd[216] = 2.0E0*I_KINETIC_I2x4y_Dyz_a-3*I_KINETIC_G2x2y_Dyz;
  abcd[217] = 2.0E0*I_KINETIC_I2x3yz_Dyz_a-2*I_KINETIC_G2xyz_Dyz;
  abcd[218] = 2.0E0*I_KINETIC_I2x2y2z_Dyz_a-1*I_KINETIC_G2x2z_Dyz;
  abcd[219] = 2.0E0*I_KINETIC_I2xy3z_Dyz_a;
  abcd[220] = 2.0E0*I_KINETIC_Ix5y_Dyz_a-4*I_KINETIC_Gx3y_Dyz;
  abcd[221] = 2.0E0*I_KINETIC_Ix4yz_Dyz_a-3*I_KINETIC_Gx2yz_Dyz;
  abcd[222] = 2.0E0*I_KINETIC_Ix3y2z_Dyz_a-2*I_KINETIC_Gxy2z_Dyz;
  abcd[223] = 2.0E0*I_KINETIC_Ix2y3z_Dyz_a-1*I_KINETIC_Gx3z_Dyz;
  abcd[224] = 2.0E0*I_KINETIC_Ixy4z_Dyz_a;
  abcd[225] = 2.0E0*I_KINETIC_I6y_Dyz_a-5*I_KINETIC_G4y_Dyz;
  abcd[226] = 2.0E0*I_KINETIC_I5yz_Dyz_a-4*I_KINETIC_G3yz_Dyz;
  abcd[227] = 2.0E0*I_KINETIC_I4y2z_Dyz_a-3*I_KINETIC_G2y2z_Dyz;
  abcd[228] = 2.0E0*I_KINETIC_I3y3z_Dyz_a-2*I_KINETIC_Gy3z_Dyz;
  abcd[229] = 2.0E0*I_KINETIC_I2y4z_Dyz_a-1*I_KINETIC_G4z_Dyz;
  abcd[230] = 2.0E0*I_KINETIC_Iy5z_Dyz_a;
  abcd[231] = 2.0E0*I_KINETIC_I5xy_D2z_a;
  abcd[232] = 2.0E0*I_KINETIC_I4x2y_D2z_a-1*I_KINETIC_G4x_D2z;
  abcd[233] = 2.0E0*I_KINETIC_I4xyz_D2z_a;
  abcd[234] = 2.0E0*I_KINETIC_I3x3y_D2z_a-2*I_KINETIC_G3xy_D2z;
  abcd[235] = 2.0E0*I_KINETIC_I3x2yz_D2z_a-1*I_KINETIC_G3xz_D2z;
  abcd[236] = 2.0E0*I_KINETIC_I3xy2z_D2z_a;
  abcd[237] = 2.0E0*I_KINETIC_I2x4y_D2z_a-3*I_KINETIC_G2x2y_D2z;
  abcd[238] = 2.0E0*I_KINETIC_I2x3yz_D2z_a-2*I_KINETIC_G2xyz_D2z;
  abcd[239] = 2.0E0*I_KINETIC_I2x2y2z_D2z_a-1*I_KINETIC_G2x2z_D2z;
  abcd[240] = 2.0E0*I_KINETIC_I2xy3z_D2z_a;
  abcd[241] = 2.0E0*I_KINETIC_Ix5y_D2z_a-4*I_KINETIC_Gx3y_D2z;
  abcd[242] = 2.0E0*I_KINETIC_Ix4yz_D2z_a-3*I_KINETIC_Gx2yz_D2z;
  abcd[243] = 2.0E0*I_KINETIC_Ix3y2z_D2z_a-2*I_KINETIC_Gxy2z_D2z;
  abcd[244] = 2.0E0*I_KINETIC_Ix2y3z_D2z_a-1*I_KINETIC_Gx3z_D2z;
  abcd[245] = 2.0E0*I_KINETIC_Ixy4z_D2z_a;
  abcd[246] = 2.0E0*I_KINETIC_I6y_D2z_a-5*I_KINETIC_G4y_D2z;
  abcd[247] = 2.0E0*I_KINETIC_I5yz_D2z_a-4*I_KINETIC_G3yz_D2z;
  abcd[248] = 2.0E0*I_KINETIC_I4y2z_D2z_a-3*I_KINETIC_G2y2z_D2z;
  abcd[249] = 2.0E0*I_KINETIC_I3y3z_D2z_a-2*I_KINETIC_Gy3z_D2z;
  abcd[250] = 2.0E0*I_KINETIC_I2y4z_D2z_a-1*I_KINETIC_G4z_D2z;
  abcd[251] = 2.0E0*I_KINETIC_Iy5z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_D_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_D_a
   * RHS shell quartet name: SQ_KINETIC_G_D
   ************************************************************/
  abcd[252] = 2.0E0*I_KINETIC_I5xz_D2x_a;
  abcd[253] = 2.0E0*I_KINETIC_I4xyz_D2x_a;
  abcd[254] = 2.0E0*I_KINETIC_I4x2z_D2x_a-1*I_KINETIC_G4x_D2x;
  abcd[255] = 2.0E0*I_KINETIC_I3x2yz_D2x_a;
  abcd[256] = 2.0E0*I_KINETIC_I3xy2z_D2x_a-1*I_KINETIC_G3xy_D2x;
  abcd[257] = 2.0E0*I_KINETIC_I3x3z_D2x_a-2*I_KINETIC_G3xz_D2x;
  abcd[258] = 2.0E0*I_KINETIC_I2x3yz_D2x_a;
  abcd[259] = 2.0E0*I_KINETIC_I2x2y2z_D2x_a-1*I_KINETIC_G2x2y_D2x;
  abcd[260] = 2.0E0*I_KINETIC_I2xy3z_D2x_a-2*I_KINETIC_G2xyz_D2x;
  abcd[261] = 2.0E0*I_KINETIC_I2x4z_D2x_a-3*I_KINETIC_G2x2z_D2x;
  abcd[262] = 2.0E0*I_KINETIC_Ix4yz_D2x_a;
  abcd[263] = 2.0E0*I_KINETIC_Ix3y2z_D2x_a-1*I_KINETIC_Gx3y_D2x;
  abcd[264] = 2.0E0*I_KINETIC_Ix2y3z_D2x_a-2*I_KINETIC_Gx2yz_D2x;
  abcd[265] = 2.0E0*I_KINETIC_Ixy4z_D2x_a-3*I_KINETIC_Gxy2z_D2x;
  abcd[266] = 2.0E0*I_KINETIC_Ix5z_D2x_a-4*I_KINETIC_Gx3z_D2x;
  abcd[267] = 2.0E0*I_KINETIC_I5yz_D2x_a;
  abcd[268] = 2.0E0*I_KINETIC_I4y2z_D2x_a-1*I_KINETIC_G4y_D2x;
  abcd[269] = 2.0E0*I_KINETIC_I3y3z_D2x_a-2*I_KINETIC_G3yz_D2x;
  abcd[270] = 2.0E0*I_KINETIC_I2y4z_D2x_a-3*I_KINETIC_G2y2z_D2x;
  abcd[271] = 2.0E0*I_KINETIC_Iy5z_D2x_a-4*I_KINETIC_Gy3z_D2x;
  abcd[272] = 2.0E0*I_KINETIC_I6z_D2x_a-5*I_KINETIC_G4z_D2x;
  abcd[273] = 2.0E0*I_KINETIC_I5xz_Dxy_a;
  abcd[274] = 2.0E0*I_KINETIC_I4xyz_Dxy_a;
  abcd[275] = 2.0E0*I_KINETIC_I4x2z_Dxy_a-1*I_KINETIC_G4x_Dxy;
  abcd[276] = 2.0E0*I_KINETIC_I3x2yz_Dxy_a;
  abcd[277] = 2.0E0*I_KINETIC_I3xy2z_Dxy_a-1*I_KINETIC_G3xy_Dxy;
  abcd[278] = 2.0E0*I_KINETIC_I3x3z_Dxy_a-2*I_KINETIC_G3xz_Dxy;
  abcd[279] = 2.0E0*I_KINETIC_I2x3yz_Dxy_a;
  abcd[280] = 2.0E0*I_KINETIC_I2x2y2z_Dxy_a-1*I_KINETIC_G2x2y_Dxy;
  abcd[281] = 2.0E0*I_KINETIC_I2xy3z_Dxy_a-2*I_KINETIC_G2xyz_Dxy;
  abcd[282] = 2.0E0*I_KINETIC_I2x4z_Dxy_a-3*I_KINETIC_G2x2z_Dxy;
  abcd[283] = 2.0E0*I_KINETIC_Ix4yz_Dxy_a;
  abcd[284] = 2.0E0*I_KINETIC_Ix3y2z_Dxy_a-1*I_KINETIC_Gx3y_Dxy;
  abcd[285] = 2.0E0*I_KINETIC_Ix2y3z_Dxy_a-2*I_KINETIC_Gx2yz_Dxy;
  abcd[286] = 2.0E0*I_KINETIC_Ixy4z_Dxy_a-3*I_KINETIC_Gxy2z_Dxy;
  abcd[287] = 2.0E0*I_KINETIC_Ix5z_Dxy_a-4*I_KINETIC_Gx3z_Dxy;
  abcd[288] = 2.0E0*I_KINETIC_I5yz_Dxy_a;
  abcd[289] = 2.0E0*I_KINETIC_I4y2z_Dxy_a-1*I_KINETIC_G4y_Dxy;
  abcd[290] = 2.0E0*I_KINETIC_I3y3z_Dxy_a-2*I_KINETIC_G3yz_Dxy;
  abcd[291] = 2.0E0*I_KINETIC_I2y4z_Dxy_a-3*I_KINETIC_G2y2z_Dxy;
  abcd[292] = 2.0E0*I_KINETIC_Iy5z_Dxy_a-4*I_KINETIC_Gy3z_Dxy;
  abcd[293] = 2.0E0*I_KINETIC_I6z_Dxy_a-5*I_KINETIC_G4z_Dxy;
  abcd[294] = 2.0E0*I_KINETIC_I5xz_Dxz_a;
  abcd[295] = 2.0E0*I_KINETIC_I4xyz_Dxz_a;
  abcd[296] = 2.0E0*I_KINETIC_I4x2z_Dxz_a-1*I_KINETIC_G4x_Dxz;
  abcd[297] = 2.0E0*I_KINETIC_I3x2yz_Dxz_a;
  abcd[298] = 2.0E0*I_KINETIC_I3xy2z_Dxz_a-1*I_KINETIC_G3xy_Dxz;
  abcd[299] = 2.0E0*I_KINETIC_I3x3z_Dxz_a-2*I_KINETIC_G3xz_Dxz;
  abcd[300] = 2.0E0*I_KINETIC_I2x3yz_Dxz_a;
  abcd[301] = 2.0E0*I_KINETIC_I2x2y2z_Dxz_a-1*I_KINETIC_G2x2y_Dxz;
  abcd[302] = 2.0E0*I_KINETIC_I2xy3z_Dxz_a-2*I_KINETIC_G2xyz_Dxz;
  abcd[303] = 2.0E0*I_KINETIC_I2x4z_Dxz_a-3*I_KINETIC_G2x2z_Dxz;
  abcd[304] = 2.0E0*I_KINETIC_Ix4yz_Dxz_a;
  abcd[305] = 2.0E0*I_KINETIC_Ix3y2z_Dxz_a-1*I_KINETIC_Gx3y_Dxz;
  abcd[306] = 2.0E0*I_KINETIC_Ix2y3z_Dxz_a-2*I_KINETIC_Gx2yz_Dxz;
  abcd[307] = 2.0E0*I_KINETIC_Ixy4z_Dxz_a-3*I_KINETIC_Gxy2z_Dxz;
  abcd[308] = 2.0E0*I_KINETIC_Ix5z_Dxz_a-4*I_KINETIC_Gx3z_Dxz;
  abcd[309] = 2.0E0*I_KINETIC_I5yz_Dxz_a;
  abcd[310] = 2.0E0*I_KINETIC_I4y2z_Dxz_a-1*I_KINETIC_G4y_Dxz;
  abcd[311] = 2.0E0*I_KINETIC_I3y3z_Dxz_a-2*I_KINETIC_G3yz_Dxz;
  abcd[312] = 2.0E0*I_KINETIC_I2y4z_Dxz_a-3*I_KINETIC_G2y2z_Dxz;
  abcd[313] = 2.0E0*I_KINETIC_Iy5z_Dxz_a-4*I_KINETIC_Gy3z_Dxz;
  abcd[314] = 2.0E0*I_KINETIC_I6z_Dxz_a-5*I_KINETIC_G4z_Dxz;
  abcd[315] = 2.0E0*I_KINETIC_I5xz_D2y_a;
  abcd[316] = 2.0E0*I_KINETIC_I4xyz_D2y_a;
  abcd[317] = 2.0E0*I_KINETIC_I4x2z_D2y_a-1*I_KINETIC_G4x_D2y;
  abcd[318] = 2.0E0*I_KINETIC_I3x2yz_D2y_a;
  abcd[319] = 2.0E0*I_KINETIC_I3xy2z_D2y_a-1*I_KINETIC_G3xy_D2y;
  abcd[320] = 2.0E0*I_KINETIC_I3x3z_D2y_a-2*I_KINETIC_G3xz_D2y;
  abcd[321] = 2.0E0*I_KINETIC_I2x3yz_D2y_a;
  abcd[322] = 2.0E0*I_KINETIC_I2x2y2z_D2y_a-1*I_KINETIC_G2x2y_D2y;
  abcd[323] = 2.0E0*I_KINETIC_I2xy3z_D2y_a-2*I_KINETIC_G2xyz_D2y;
  abcd[324] = 2.0E0*I_KINETIC_I2x4z_D2y_a-3*I_KINETIC_G2x2z_D2y;
  abcd[325] = 2.0E0*I_KINETIC_Ix4yz_D2y_a;
  abcd[326] = 2.0E0*I_KINETIC_Ix3y2z_D2y_a-1*I_KINETIC_Gx3y_D2y;
  abcd[327] = 2.0E0*I_KINETIC_Ix2y3z_D2y_a-2*I_KINETIC_Gx2yz_D2y;
  abcd[328] = 2.0E0*I_KINETIC_Ixy4z_D2y_a-3*I_KINETIC_Gxy2z_D2y;
  abcd[329] = 2.0E0*I_KINETIC_Ix5z_D2y_a-4*I_KINETIC_Gx3z_D2y;
  abcd[330] = 2.0E0*I_KINETIC_I5yz_D2y_a;
  abcd[331] = 2.0E0*I_KINETIC_I4y2z_D2y_a-1*I_KINETIC_G4y_D2y;
  abcd[332] = 2.0E0*I_KINETIC_I3y3z_D2y_a-2*I_KINETIC_G3yz_D2y;
  abcd[333] = 2.0E0*I_KINETIC_I2y4z_D2y_a-3*I_KINETIC_G2y2z_D2y;
  abcd[334] = 2.0E0*I_KINETIC_Iy5z_D2y_a-4*I_KINETIC_Gy3z_D2y;
  abcd[335] = 2.0E0*I_KINETIC_I6z_D2y_a-5*I_KINETIC_G4z_D2y;
  abcd[336] = 2.0E0*I_KINETIC_I5xz_Dyz_a;
  abcd[337] = 2.0E0*I_KINETIC_I4xyz_Dyz_a;
  abcd[338] = 2.0E0*I_KINETIC_I4x2z_Dyz_a-1*I_KINETIC_G4x_Dyz;
  abcd[339] = 2.0E0*I_KINETIC_I3x2yz_Dyz_a;
  abcd[340] = 2.0E0*I_KINETIC_I3xy2z_Dyz_a-1*I_KINETIC_G3xy_Dyz;
  abcd[341] = 2.0E0*I_KINETIC_I3x3z_Dyz_a-2*I_KINETIC_G3xz_Dyz;
  abcd[342] = 2.0E0*I_KINETIC_I2x3yz_Dyz_a;
  abcd[343] = 2.0E0*I_KINETIC_I2x2y2z_Dyz_a-1*I_KINETIC_G2x2y_Dyz;
  abcd[344] = 2.0E0*I_KINETIC_I2xy3z_Dyz_a-2*I_KINETIC_G2xyz_Dyz;
  abcd[345] = 2.0E0*I_KINETIC_I2x4z_Dyz_a-3*I_KINETIC_G2x2z_Dyz;
  abcd[346] = 2.0E0*I_KINETIC_Ix4yz_Dyz_a;
  abcd[347] = 2.0E0*I_KINETIC_Ix3y2z_Dyz_a-1*I_KINETIC_Gx3y_Dyz;
  abcd[348] = 2.0E0*I_KINETIC_Ix2y3z_Dyz_a-2*I_KINETIC_Gx2yz_Dyz;
  abcd[349] = 2.0E0*I_KINETIC_Ixy4z_Dyz_a-3*I_KINETIC_Gxy2z_Dyz;
  abcd[350] = 2.0E0*I_KINETIC_Ix5z_Dyz_a-4*I_KINETIC_Gx3z_Dyz;
  abcd[351] = 2.0E0*I_KINETIC_I5yz_Dyz_a;
  abcd[352] = 2.0E0*I_KINETIC_I4y2z_Dyz_a-1*I_KINETIC_G4y_Dyz;
  abcd[353] = 2.0E0*I_KINETIC_I3y3z_Dyz_a-2*I_KINETIC_G3yz_Dyz;
  abcd[354] = 2.0E0*I_KINETIC_I2y4z_Dyz_a-3*I_KINETIC_G2y2z_Dyz;
  abcd[355] = 2.0E0*I_KINETIC_Iy5z_Dyz_a-4*I_KINETIC_Gy3z_Dyz;
  abcd[356] = 2.0E0*I_KINETIC_I6z_Dyz_a-5*I_KINETIC_G4z_Dyz;
  abcd[357] = 2.0E0*I_KINETIC_I5xz_D2z_a;
  abcd[358] = 2.0E0*I_KINETIC_I4xyz_D2z_a;
  abcd[359] = 2.0E0*I_KINETIC_I4x2z_D2z_a-1*I_KINETIC_G4x_D2z;
  abcd[360] = 2.0E0*I_KINETIC_I3x2yz_D2z_a;
  abcd[361] = 2.0E0*I_KINETIC_I3xy2z_D2z_a-1*I_KINETIC_G3xy_D2z;
  abcd[362] = 2.0E0*I_KINETIC_I3x3z_D2z_a-2*I_KINETIC_G3xz_D2z;
  abcd[363] = 2.0E0*I_KINETIC_I2x3yz_D2z_a;
  abcd[364] = 2.0E0*I_KINETIC_I2x2y2z_D2z_a-1*I_KINETIC_G2x2y_D2z;
  abcd[365] = 2.0E0*I_KINETIC_I2xy3z_D2z_a-2*I_KINETIC_G2xyz_D2z;
  abcd[366] = 2.0E0*I_KINETIC_I2x4z_D2z_a-3*I_KINETIC_G2x2z_D2z;
  abcd[367] = 2.0E0*I_KINETIC_Ix4yz_D2z_a;
  abcd[368] = 2.0E0*I_KINETIC_Ix3y2z_D2z_a-1*I_KINETIC_Gx3y_D2z;
  abcd[369] = 2.0E0*I_KINETIC_Ix2y3z_D2z_a-2*I_KINETIC_Gx2yz_D2z;
  abcd[370] = 2.0E0*I_KINETIC_Ixy4z_D2z_a-3*I_KINETIC_Gxy2z_D2z;
  abcd[371] = 2.0E0*I_KINETIC_Ix5z_D2z_a-4*I_KINETIC_Gx3z_D2z;
  abcd[372] = 2.0E0*I_KINETIC_I5yz_D2z_a;
  abcd[373] = 2.0E0*I_KINETIC_I4y2z_D2z_a-1*I_KINETIC_G4y_D2z;
  abcd[374] = 2.0E0*I_KINETIC_I3y3z_D2z_a-2*I_KINETIC_G3yz_D2z;
  abcd[375] = 2.0E0*I_KINETIC_I2y4z_D2z_a-3*I_KINETIC_G2y2z_D2z;
  abcd[376] = 2.0E0*I_KINETIC_Iy5z_D2z_a-4*I_KINETIC_Gy3z_D2z;
  abcd[377] = 2.0E0*I_KINETIC_I6z_D2z_a-5*I_KINETIC_G4z_D2z;
}
