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
// BRA1
// X
// Y
// Z
// ####

void hgp_os_kinetic_g_g_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_H5x_G4x_a = 0.0E0;
  Double I_KINETIC_H4xy_G4x_a = 0.0E0;
  Double I_KINETIC_H4xz_G4x_a = 0.0E0;
  Double I_KINETIC_H3x2y_G4x_a = 0.0E0;
  Double I_KINETIC_H3xyz_G4x_a = 0.0E0;
  Double I_KINETIC_H3x2z_G4x_a = 0.0E0;
  Double I_KINETIC_H2x3y_G4x_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G4x_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G4x_a = 0.0E0;
  Double I_KINETIC_H2x3z_G4x_a = 0.0E0;
  Double I_KINETIC_Hx4y_G4x_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G4x_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G4x_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G4x_a = 0.0E0;
  Double I_KINETIC_Hx4z_G4x_a = 0.0E0;
  Double I_KINETIC_H5y_G4x_a = 0.0E0;
  Double I_KINETIC_H4yz_G4x_a = 0.0E0;
  Double I_KINETIC_H3y2z_G4x_a = 0.0E0;
  Double I_KINETIC_H2y3z_G4x_a = 0.0E0;
  Double I_KINETIC_Hy4z_G4x_a = 0.0E0;
  Double I_KINETIC_H5z_G4x_a = 0.0E0;
  Double I_KINETIC_H5x_G3xy_a = 0.0E0;
  Double I_KINETIC_H4xy_G3xy_a = 0.0E0;
  Double I_KINETIC_H4xz_G3xy_a = 0.0E0;
  Double I_KINETIC_H3x2y_G3xy_a = 0.0E0;
  Double I_KINETIC_H3xyz_G3xy_a = 0.0E0;
  Double I_KINETIC_H3x2z_G3xy_a = 0.0E0;
  Double I_KINETIC_H2x3y_G3xy_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G3xy_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G3xy_a = 0.0E0;
  Double I_KINETIC_H2x3z_G3xy_a = 0.0E0;
  Double I_KINETIC_Hx4y_G3xy_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G3xy_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G3xy_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G3xy_a = 0.0E0;
  Double I_KINETIC_Hx4z_G3xy_a = 0.0E0;
  Double I_KINETIC_H5y_G3xy_a = 0.0E0;
  Double I_KINETIC_H4yz_G3xy_a = 0.0E0;
  Double I_KINETIC_H3y2z_G3xy_a = 0.0E0;
  Double I_KINETIC_H2y3z_G3xy_a = 0.0E0;
  Double I_KINETIC_Hy4z_G3xy_a = 0.0E0;
  Double I_KINETIC_H5z_G3xy_a = 0.0E0;
  Double I_KINETIC_H5x_G3xz_a = 0.0E0;
  Double I_KINETIC_H4xy_G3xz_a = 0.0E0;
  Double I_KINETIC_H4xz_G3xz_a = 0.0E0;
  Double I_KINETIC_H3x2y_G3xz_a = 0.0E0;
  Double I_KINETIC_H3xyz_G3xz_a = 0.0E0;
  Double I_KINETIC_H3x2z_G3xz_a = 0.0E0;
  Double I_KINETIC_H2x3y_G3xz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G3xz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G3xz_a = 0.0E0;
  Double I_KINETIC_H2x3z_G3xz_a = 0.0E0;
  Double I_KINETIC_Hx4y_G3xz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G3xz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G3xz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G3xz_a = 0.0E0;
  Double I_KINETIC_Hx4z_G3xz_a = 0.0E0;
  Double I_KINETIC_H5y_G3xz_a = 0.0E0;
  Double I_KINETIC_H4yz_G3xz_a = 0.0E0;
  Double I_KINETIC_H3y2z_G3xz_a = 0.0E0;
  Double I_KINETIC_H2y3z_G3xz_a = 0.0E0;
  Double I_KINETIC_Hy4z_G3xz_a = 0.0E0;
  Double I_KINETIC_H5z_G3xz_a = 0.0E0;
  Double I_KINETIC_H5x_G2x2y_a = 0.0E0;
  Double I_KINETIC_H4xy_G2x2y_a = 0.0E0;
  Double I_KINETIC_H4xz_G2x2y_a = 0.0E0;
  Double I_KINETIC_H3x2y_G2x2y_a = 0.0E0;
  Double I_KINETIC_H3xyz_G2x2y_a = 0.0E0;
  Double I_KINETIC_H3x2z_G2x2y_a = 0.0E0;
  Double I_KINETIC_H2x3y_G2x2y_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G2x2y_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G2x2y_a = 0.0E0;
  Double I_KINETIC_H2x3z_G2x2y_a = 0.0E0;
  Double I_KINETIC_Hx4y_G2x2y_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G2x2y_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G2x2y_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G2x2y_a = 0.0E0;
  Double I_KINETIC_Hx4z_G2x2y_a = 0.0E0;
  Double I_KINETIC_H5y_G2x2y_a = 0.0E0;
  Double I_KINETIC_H4yz_G2x2y_a = 0.0E0;
  Double I_KINETIC_H3y2z_G2x2y_a = 0.0E0;
  Double I_KINETIC_H2y3z_G2x2y_a = 0.0E0;
  Double I_KINETIC_Hy4z_G2x2y_a = 0.0E0;
  Double I_KINETIC_H5z_G2x2y_a = 0.0E0;
  Double I_KINETIC_H5x_G2xyz_a = 0.0E0;
  Double I_KINETIC_H4xy_G2xyz_a = 0.0E0;
  Double I_KINETIC_H4xz_G2xyz_a = 0.0E0;
  Double I_KINETIC_H3x2y_G2xyz_a = 0.0E0;
  Double I_KINETIC_H3xyz_G2xyz_a = 0.0E0;
  Double I_KINETIC_H3x2z_G2xyz_a = 0.0E0;
  Double I_KINETIC_H2x3y_G2xyz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G2xyz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G2xyz_a = 0.0E0;
  Double I_KINETIC_H2x3z_G2xyz_a = 0.0E0;
  Double I_KINETIC_Hx4y_G2xyz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G2xyz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G2xyz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G2xyz_a = 0.0E0;
  Double I_KINETIC_Hx4z_G2xyz_a = 0.0E0;
  Double I_KINETIC_H5y_G2xyz_a = 0.0E0;
  Double I_KINETIC_H4yz_G2xyz_a = 0.0E0;
  Double I_KINETIC_H3y2z_G2xyz_a = 0.0E0;
  Double I_KINETIC_H2y3z_G2xyz_a = 0.0E0;
  Double I_KINETIC_Hy4z_G2xyz_a = 0.0E0;
  Double I_KINETIC_H5z_G2xyz_a = 0.0E0;
  Double I_KINETIC_H5x_G2x2z_a = 0.0E0;
  Double I_KINETIC_H4xy_G2x2z_a = 0.0E0;
  Double I_KINETIC_H4xz_G2x2z_a = 0.0E0;
  Double I_KINETIC_H3x2y_G2x2z_a = 0.0E0;
  Double I_KINETIC_H3xyz_G2x2z_a = 0.0E0;
  Double I_KINETIC_H3x2z_G2x2z_a = 0.0E0;
  Double I_KINETIC_H2x3y_G2x2z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G2x2z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G2x2z_a = 0.0E0;
  Double I_KINETIC_H2x3z_G2x2z_a = 0.0E0;
  Double I_KINETIC_Hx4y_G2x2z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G2x2z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G2x2z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G2x2z_a = 0.0E0;
  Double I_KINETIC_Hx4z_G2x2z_a = 0.0E0;
  Double I_KINETIC_H5y_G2x2z_a = 0.0E0;
  Double I_KINETIC_H4yz_G2x2z_a = 0.0E0;
  Double I_KINETIC_H3y2z_G2x2z_a = 0.0E0;
  Double I_KINETIC_H2y3z_G2x2z_a = 0.0E0;
  Double I_KINETIC_Hy4z_G2x2z_a = 0.0E0;
  Double I_KINETIC_H5z_G2x2z_a = 0.0E0;
  Double I_KINETIC_H5x_Gx3y_a = 0.0E0;
  Double I_KINETIC_H4xy_Gx3y_a = 0.0E0;
  Double I_KINETIC_H4xz_Gx3y_a = 0.0E0;
  Double I_KINETIC_H3x2y_Gx3y_a = 0.0E0;
  Double I_KINETIC_H3xyz_Gx3y_a = 0.0E0;
  Double I_KINETIC_H3x2z_Gx3y_a = 0.0E0;
  Double I_KINETIC_H2x3y_Gx3y_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Gx3y_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Gx3y_a = 0.0E0;
  Double I_KINETIC_H2x3z_Gx3y_a = 0.0E0;
  Double I_KINETIC_Hx4y_Gx3y_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Gx3y_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Gx3y_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Gx3y_a = 0.0E0;
  Double I_KINETIC_Hx4z_Gx3y_a = 0.0E0;
  Double I_KINETIC_H5y_Gx3y_a = 0.0E0;
  Double I_KINETIC_H4yz_Gx3y_a = 0.0E0;
  Double I_KINETIC_H3y2z_Gx3y_a = 0.0E0;
  Double I_KINETIC_H2y3z_Gx3y_a = 0.0E0;
  Double I_KINETIC_Hy4z_Gx3y_a = 0.0E0;
  Double I_KINETIC_H5z_Gx3y_a = 0.0E0;
  Double I_KINETIC_H5x_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H4xy_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H4xz_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H3x2y_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H3xyz_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H3x2z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H2x3y_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H2x3z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_Hx4y_Gx2yz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Gx2yz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_Hx4z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H5y_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H4yz_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H3y2z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H2y3z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_Hy4z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H5z_Gx2yz_a = 0.0E0;
  Double I_KINETIC_H5x_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H4xy_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H4xz_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H3x2y_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H3xyz_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H3x2z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H2x3y_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H2x3z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_Hx4y_Gxy2z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Gxy2z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_Hx4z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H5y_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H4yz_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H3y2z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H2y3z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_Hy4z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H5z_Gxy2z_a = 0.0E0;
  Double I_KINETIC_H5x_Gx3z_a = 0.0E0;
  Double I_KINETIC_H4xy_Gx3z_a = 0.0E0;
  Double I_KINETIC_H4xz_Gx3z_a = 0.0E0;
  Double I_KINETIC_H3x2y_Gx3z_a = 0.0E0;
  Double I_KINETIC_H3xyz_Gx3z_a = 0.0E0;
  Double I_KINETIC_H3x2z_Gx3z_a = 0.0E0;
  Double I_KINETIC_H2x3y_Gx3z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Gx3z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Gx3z_a = 0.0E0;
  Double I_KINETIC_H2x3z_Gx3z_a = 0.0E0;
  Double I_KINETIC_Hx4y_Gx3z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Gx3z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Gx3z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Gx3z_a = 0.0E0;
  Double I_KINETIC_Hx4z_Gx3z_a = 0.0E0;
  Double I_KINETIC_H5y_Gx3z_a = 0.0E0;
  Double I_KINETIC_H4yz_Gx3z_a = 0.0E0;
  Double I_KINETIC_H3y2z_Gx3z_a = 0.0E0;
  Double I_KINETIC_H2y3z_Gx3z_a = 0.0E0;
  Double I_KINETIC_Hy4z_Gx3z_a = 0.0E0;
  Double I_KINETIC_H5z_Gx3z_a = 0.0E0;
  Double I_KINETIC_H5x_G4y_a = 0.0E0;
  Double I_KINETIC_H4xy_G4y_a = 0.0E0;
  Double I_KINETIC_H4xz_G4y_a = 0.0E0;
  Double I_KINETIC_H3x2y_G4y_a = 0.0E0;
  Double I_KINETIC_H3xyz_G4y_a = 0.0E0;
  Double I_KINETIC_H3x2z_G4y_a = 0.0E0;
  Double I_KINETIC_H2x3y_G4y_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G4y_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G4y_a = 0.0E0;
  Double I_KINETIC_H2x3z_G4y_a = 0.0E0;
  Double I_KINETIC_Hx4y_G4y_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G4y_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G4y_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G4y_a = 0.0E0;
  Double I_KINETIC_Hx4z_G4y_a = 0.0E0;
  Double I_KINETIC_H5y_G4y_a = 0.0E0;
  Double I_KINETIC_H4yz_G4y_a = 0.0E0;
  Double I_KINETIC_H3y2z_G4y_a = 0.0E0;
  Double I_KINETIC_H2y3z_G4y_a = 0.0E0;
  Double I_KINETIC_Hy4z_G4y_a = 0.0E0;
  Double I_KINETIC_H5z_G4y_a = 0.0E0;
  Double I_KINETIC_H5x_G3yz_a = 0.0E0;
  Double I_KINETIC_H4xy_G3yz_a = 0.0E0;
  Double I_KINETIC_H4xz_G3yz_a = 0.0E0;
  Double I_KINETIC_H3x2y_G3yz_a = 0.0E0;
  Double I_KINETIC_H3xyz_G3yz_a = 0.0E0;
  Double I_KINETIC_H3x2z_G3yz_a = 0.0E0;
  Double I_KINETIC_H2x3y_G3yz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G3yz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G3yz_a = 0.0E0;
  Double I_KINETIC_H2x3z_G3yz_a = 0.0E0;
  Double I_KINETIC_Hx4y_G3yz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G3yz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G3yz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G3yz_a = 0.0E0;
  Double I_KINETIC_Hx4z_G3yz_a = 0.0E0;
  Double I_KINETIC_H5y_G3yz_a = 0.0E0;
  Double I_KINETIC_H4yz_G3yz_a = 0.0E0;
  Double I_KINETIC_H3y2z_G3yz_a = 0.0E0;
  Double I_KINETIC_H2y3z_G3yz_a = 0.0E0;
  Double I_KINETIC_Hy4z_G3yz_a = 0.0E0;
  Double I_KINETIC_H5z_G3yz_a = 0.0E0;
  Double I_KINETIC_H5x_G2y2z_a = 0.0E0;
  Double I_KINETIC_H4xy_G2y2z_a = 0.0E0;
  Double I_KINETIC_H4xz_G2y2z_a = 0.0E0;
  Double I_KINETIC_H3x2y_G2y2z_a = 0.0E0;
  Double I_KINETIC_H3xyz_G2y2z_a = 0.0E0;
  Double I_KINETIC_H3x2z_G2y2z_a = 0.0E0;
  Double I_KINETIC_H2x3y_G2y2z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G2y2z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G2y2z_a = 0.0E0;
  Double I_KINETIC_H2x3z_G2y2z_a = 0.0E0;
  Double I_KINETIC_Hx4y_G2y2z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G2y2z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G2y2z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G2y2z_a = 0.0E0;
  Double I_KINETIC_Hx4z_G2y2z_a = 0.0E0;
  Double I_KINETIC_H5y_G2y2z_a = 0.0E0;
  Double I_KINETIC_H4yz_G2y2z_a = 0.0E0;
  Double I_KINETIC_H3y2z_G2y2z_a = 0.0E0;
  Double I_KINETIC_H2y3z_G2y2z_a = 0.0E0;
  Double I_KINETIC_Hy4z_G2y2z_a = 0.0E0;
  Double I_KINETIC_H5z_G2y2z_a = 0.0E0;
  Double I_KINETIC_H5x_Gy3z_a = 0.0E0;
  Double I_KINETIC_H4xy_Gy3z_a = 0.0E0;
  Double I_KINETIC_H4xz_Gy3z_a = 0.0E0;
  Double I_KINETIC_H3x2y_Gy3z_a = 0.0E0;
  Double I_KINETIC_H3xyz_Gy3z_a = 0.0E0;
  Double I_KINETIC_H3x2z_Gy3z_a = 0.0E0;
  Double I_KINETIC_H2x3y_Gy3z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Gy3z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Gy3z_a = 0.0E0;
  Double I_KINETIC_H2x3z_Gy3z_a = 0.0E0;
  Double I_KINETIC_Hx4y_Gy3z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Gy3z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Gy3z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Gy3z_a = 0.0E0;
  Double I_KINETIC_Hx4z_Gy3z_a = 0.0E0;
  Double I_KINETIC_H5y_Gy3z_a = 0.0E0;
  Double I_KINETIC_H4yz_Gy3z_a = 0.0E0;
  Double I_KINETIC_H3y2z_Gy3z_a = 0.0E0;
  Double I_KINETIC_H2y3z_Gy3z_a = 0.0E0;
  Double I_KINETIC_Hy4z_Gy3z_a = 0.0E0;
  Double I_KINETIC_H5z_Gy3z_a = 0.0E0;
  Double I_KINETIC_H5x_G4z_a = 0.0E0;
  Double I_KINETIC_H4xy_G4z_a = 0.0E0;
  Double I_KINETIC_H4xz_G4z_a = 0.0E0;
  Double I_KINETIC_H3x2y_G4z_a = 0.0E0;
  Double I_KINETIC_H3xyz_G4z_a = 0.0E0;
  Double I_KINETIC_H3x2z_G4z_a = 0.0E0;
  Double I_KINETIC_H2x3y_G4z_a = 0.0E0;
  Double I_KINETIC_H2x2yz_G4z_a = 0.0E0;
  Double I_KINETIC_H2xy2z_G4z_a = 0.0E0;
  Double I_KINETIC_H2x3z_G4z_a = 0.0E0;
  Double I_KINETIC_Hx4y_G4z_a = 0.0E0;
  Double I_KINETIC_Hx3yz_G4z_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_G4z_a = 0.0E0;
  Double I_KINETIC_Hxy3z_G4z_a = 0.0E0;
  Double I_KINETIC_Hx4z_G4z_a = 0.0E0;
  Double I_KINETIC_H5y_G4z_a = 0.0E0;
  Double I_KINETIC_H4yz_G4z_a = 0.0E0;
  Double I_KINETIC_H3y2z_G4z_a = 0.0E0;
  Double I_KINETIC_H2y3z_G4z_a = 0.0E0;
  Double I_KINETIC_Hy4z_G4z_a = 0.0E0;
  Double I_KINETIC_H5z_G4z_a = 0.0E0;
  Double I_KINETIC_F3x_G4x = 0.0E0;
  Double I_KINETIC_F2xy_G4x = 0.0E0;
  Double I_KINETIC_F2xz_G4x = 0.0E0;
  Double I_KINETIC_Fx2y_G4x = 0.0E0;
  Double I_KINETIC_Fxyz_G4x = 0.0E0;
  Double I_KINETIC_Fx2z_G4x = 0.0E0;
  Double I_KINETIC_F3y_G4x = 0.0E0;
  Double I_KINETIC_F2yz_G4x = 0.0E0;
  Double I_KINETIC_Fy2z_G4x = 0.0E0;
  Double I_KINETIC_F3z_G4x = 0.0E0;
  Double I_KINETIC_F3x_G3xy = 0.0E0;
  Double I_KINETIC_F2xy_G3xy = 0.0E0;
  Double I_KINETIC_F2xz_G3xy = 0.0E0;
  Double I_KINETIC_Fx2y_G3xy = 0.0E0;
  Double I_KINETIC_Fxyz_G3xy = 0.0E0;
  Double I_KINETIC_Fx2z_G3xy = 0.0E0;
  Double I_KINETIC_F3y_G3xy = 0.0E0;
  Double I_KINETIC_F2yz_G3xy = 0.0E0;
  Double I_KINETIC_Fy2z_G3xy = 0.0E0;
  Double I_KINETIC_F3z_G3xy = 0.0E0;
  Double I_KINETIC_F3x_G3xz = 0.0E0;
  Double I_KINETIC_F2xy_G3xz = 0.0E0;
  Double I_KINETIC_F2xz_G3xz = 0.0E0;
  Double I_KINETIC_Fx2y_G3xz = 0.0E0;
  Double I_KINETIC_Fxyz_G3xz = 0.0E0;
  Double I_KINETIC_Fx2z_G3xz = 0.0E0;
  Double I_KINETIC_F3y_G3xz = 0.0E0;
  Double I_KINETIC_F2yz_G3xz = 0.0E0;
  Double I_KINETIC_Fy2z_G3xz = 0.0E0;
  Double I_KINETIC_F3z_G3xz = 0.0E0;
  Double I_KINETIC_F3x_G2x2y = 0.0E0;
  Double I_KINETIC_F2xy_G2x2y = 0.0E0;
  Double I_KINETIC_F2xz_G2x2y = 0.0E0;
  Double I_KINETIC_Fx2y_G2x2y = 0.0E0;
  Double I_KINETIC_Fxyz_G2x2y = 0.0E0;
  Double I_KINETIC_Fx2z_G2x2y = 0.0E0;
  Double I_KINETIC_F3y_G2x2y = 0.0E0;
  Double I_KINETIC_F2yz_G2x2y = 0.0E0;
  Double I_KINETIC_Fy2z_G2x2y = 0.0E0;
  Double I_KINETIC_F3z_G2x2y = 0.0E0;
  Double I_KINETIC_F3x_G2xyz = 0.0E0;
  Double I_KINETIC_F2xy_G2xyz = 0.0E0;
  Double I_KINETIC_F2xz_G2xyz = 0.0E0;
  Double I_KINETIC_Fx2y_G2xyz = 0.0E0;
  Double I_KINETIC_Fxyz_G2xyz = 0.0E0;
  Double I_KINETIC_Fx2z_G2xyz = 0.0E0;
  Double I_KINETIC_F3y_G2xyz = 0.0E0;
  Double I_KINETIC_F2yz_G2xyz = 0.0E0;
  Double I_KINETIC_Fy2z_G2xyz = 0.0E0;
  Double I_KINETIC_F3z_G2xyz = 0.0E0;
  Double I_KINETIC_F3x_G2x2z = 0.0E0;
  Double I_KINETIC_F2xy_G2x2z = 0.0E0;
  Double I_KINETIC_F2xz_G2x2z = 0.0E0;
  Double I_KINETIC_Fx2y_G2x2z = 0.0E0;
  Double I_KINETIC_Fxyz_G2x2z = 0.0E0;
  Double I_KINETIC_Fx2z_G2x2z = 0.0E0;
  Double I_KINETIC_F3y_G2x2z = 0.0E0;
  Double I_KINETIC_F2yz_G2x2z = 0.0E0;
  Double I_KINETIC_Fy2z_G2x2z = 0.0E0;
  Double I_KINETIC_F3z_G2x2z = 0.0E0;
  Double I_KINETIC_F3x_Gx3y = 0.0E0;
  Double I_KINETIC_F2xy_Gx3y = 0.0E0;
  Double I_KINETIC_F2xz_Gx3y = 0.0E0;
  Double I_KINETIC_Fx2y_Gx3y = 0.0E0;
  Double I_KINETIC_Fxyz_Gx3y = 0.0E0;
  Double I_KINETIC_Fx2z_Gx3y = 0.0E0;
  Double I_KINETIC_F3y_Gx3y = 0.0E0;
  Double I_KINETIC_F2yz_Gx3y = 0.0E0;
  Double I_KINETIC_Fy2z_Gx3y = 0.0E0;
  Double I_KINETIC_F3z_Gx3y = 0.0E0;
  Double I_KINETIC_F3x_Gx2yz = 0.0E0;
  Double I_KINETIC_F2xy_Gx2yz = 0.0E0;
  Double I_KINETIC_F2xz_Gx2yz = 0.0E0;
  Double I_KINETIC_Fx2y_Gx2yz = 0.0E0;
  Double I_KINETIC_Fxyz_Gx2yz = 0.0E0;
  Double I_KINETIC_Fx2z_Gx2yz = 0.0E0;
  Double I_KINETIC_F3y_Gx2yz = 0.0E0;
  Double I_KINETIC_F2yz_Gx2yz = 0.0E0;
  Double I_KINETIC_Fy2z_Gx2yz = 0.0E0;
  Double I_KINETIC_F3z_Gx2yz = 0.0E0;
  Double I_KINETIC_F3x_Gxy2z = 0.0E0;
  Double I_KINETIC_F2xy_Gxy2z = 0.0E0;
  Double I_KINETIC_F2xz_Gxy2z = 0.0E0;
  Double I_KINETIC_Fx2y_Gxy2z = 0.0E0;
  Double I_KINETIC_Fxyz_Gxy2z = 0.0E0;
  Double I_KINETIC_Fx2z_Gxy2z = 0.0E0;
  Double I_KINETIC_F3y_Gxy2z = 0.0E0;
  Double I_KINETIC_F2yz_Gxy2z = 0.0E0;
  Double I_KINETIC_Fy2z_Gxy2z = 0.0E0;
  Double I_KINETIC_F3z_Gxy2z = 0.0E0;
  Double I_KINETIC_F3x_Gx3z = 0.0E0;
  Double I_KINETIC_F2xy_Gx3z = 0.0E0;
  Double I_KINETIC_F2xz_Gx3z = 0.0E0;
  Double I_KINETIC_Fx2y_Gx3z = 0.0E0;
  Double I_KINETIC_Fxyz_Gx3z = 0.0E0;
  Double I_KINETIC_Fx2z_Gx3z = 0.0E0;
  Double I_KINETIC_F3y_Gx3z = 0.0E0;
  Double I_KINETIC_F2yz_Gx3z = 0.0E0;
  Double I_KINETIC_Fy2z_Gx3z = 0.0E0;
  Double I_KINETIC_F3z_Gx3z = 0.0E0;
  Double I_KINETIC_F3x_G4y = 0.0E0;
  Double I_KINETIC_F2xy_G4y = 0.0E0;
  Double I_KINETIC_F2xz_G4y = 0.0E0;
  Double I_KINETIC_Fx2y_G4y = 0.0E0;
  Double I_KINETIC_Fxyz_G4y = 0.0E0;
  Double I_KINETIC_Fx2z_G4y = 0.0E0;
  Double I_KINETIC_F3y_G4y = 0.0E0;
  Double I_KINETIC_F2yz_G4y = 0.0E0;
  Double I_KINETIC_Fy2z_G4y = 0.0E0;
  Double I_KINETIC_F3z_G4y = 0.0E0;
  Double I_KINETIC_F3x_G3yz = 0.0E0;
  Double I_KINETIC_F2xy_G3yz = 0.0E0;
  Double I_KINETIC_F2xz_G3yz = 0.0E0;
  Double I_KINETIC_Fx2y_G3yz = 0.0E0;
  Double I_KINETIC_Fxyz_G3yz = 0.0E0;
  Double I_KINETIC_Fx2z_G3yz = 0.0E0;
  Double I_KINETIC_F3y_G3yz = 0.0E0;
  Double I_KINETIC_F2yz_G3yz = 0.0E0;
  Double I_KINETIC_Fy2z_G3yz = 0.0E0;
  Double I_KINETIC_F3z_G3yz = 0.0E0;
  Double I_KINETIC_F3x_G2y2z = 0.0E0;
  Double I_KINETIC_F2xy_G2y2z = 0.0E0;
  Double I_KINETIC_F2xz_G2y2z = 0.0E0;
  Double I_KINETIC_Fx2y_G2y2z = 0.0E0;
  Double I_KINETIC_Fxyz_G2y2z = 0.0E0;
  Double I_KINETIC_Fx2z_G2y2z = 0.0E0;
  Double I_KINETIC_F3y_G2y2z = 0.0E0;
  Double I_KINETIC_F2yz_G2y2z = 0.0E0;
  Double I_KINETIC_Fy2z_G2y2z = 0.0E0;
  Double I_KINETIC_F3z_G2y2z = 0.0E0;
  Double I_KINETIC_F3x_Gy3z = 0.0E0;
  Double I_KINETIC_F2xy_Gy3z = 0.0E0;
  Double I_KINETIC_F2xz_Gy3z = 0.0E0;
  Double I_KINETIC_Fx2y_Gy3z = 0.0E0;
  Double I_KINETIC_Fxyz_Gy3z = 0.0E0;
  Double I_KINETIC_Fx2z_Gy3z = 0.0E0;
  Double I_KINETIC_F3y_Gy3z = 0.0E0;
  Double I_KINETIC_F2yz_Gy3z = 0.0E0;
  Double I_KINETIC_Fy2z_Gy3z = 0.0E0;
  Double I_KINETIC_F3z_Gy3z = 0.0E0;
  Double I_KINETIC_F3x_G4z = 0.0E0;
  Double I_KINETIC_F2xy_G4z = 0.0E0;
  Double I_KINETIC_F2xz_G4z = 0.0E0;
  Double I_KINETIC_Fx2y_G4z = 0.0E0;
  Double I_KINETIC_Fxyz_G4z = 0.0E0;
  Double I_KINETIC_Fx2z_G4z = 0.0E0;
  Double I_KINETIC_F3y_G4z = 0.0E0;
  Double I_KINETIC_F2yz_G4z = 0.0E0;
  Double I_KINETIC_Fy2z_G4z = 0.0E0;
  Double I_KINETIC_F3z_G4z = 0.0E0;

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
     * expanding position: BRA1
     * code section is: VRR
     * totally 8 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_D2x_vrr = PAX*I_TWOBODYOVERLAP_Px_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2x_vrr = PAY*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dxz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2x_vrr = PAY*I_TWOBODYOVERLAP_Py_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dyz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Py_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Px_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Py_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2y_vrr = PAX*I_TWOBODYOVERLAP_Px_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2y_vrr = PAY*I_TWOBODYOVERLAP_Px_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_Dxz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2y_vrr = PAY*I_TWOBODYOVERLAP_Py_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_TWOBODYOVERLAP_Dyz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Px_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Py_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2z_vrr = PAX*I_TWOBODYOVERLAP_Px_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2z_vrr = PAY*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Px_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2z_vrr = PAY*I_TWOBODYOVERLAP_Py_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dyz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Py_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Pz_vrr;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_F3x_vrr = PBX*I_TWOBODYOVERLAP_D2x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F3x_vrr = PBX*I_TWOBODYOVERLAP_Dxy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F3x_vrr = PBX*I_TWOBODYOVERLAP_Dxz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3x_vrr = PBX*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F3x_vrr = PBX*I_TWOBODYOVERLAP_Dyz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2xy_vrr = PBY*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Dxy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Dxz_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Dyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Dxz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Dyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Dxy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Dxz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Dyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Dxy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Dxz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Pz_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Dyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2x_F3y_vrr = PBY*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F3y_vrr = PBY*I_TWOBODYOVERLAP_Dxy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F3y_vrr = PBY*I_TWOBODYOVERLAP_Dxz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3y_vrr = PBY*I_TWOBODYOVERLAP_D2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F3y_vrr = PBY*I_TWOBODYOVERLAP_Dyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Pz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Py_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Dxz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Dyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2x_F3z_vrr = PBZ*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Dxz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Dyz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_D2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
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
     * expanding position: BRA1
     * code section is: VRR
     * totally 8 integrals are omitted 
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
    Double I_TWOBODYOVERLAP_F3y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
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
    Double I_TWOBODYOVERLAP_F3y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_D_G
     * expanding position: BRA2
     * code section is: VRR
     * totally 45 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_G4x_vrr = PBX*I_TWOBODYOVERLAP_D2x_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_G4x_vrr = PBX*I_TWOBODYOVERLAP_D2y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_G4x_vrr = PBX*I_TWOBODYOVERLAP_D2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_G3xy_vrr = PBY*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_D2y_G3xy_vrr = PBY*I_TWOBODYOVERLAP_D2y_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F3x_vrr;
    Double I_TWOBODYOVERLAP_D2z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_D2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_D2x_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_D2y_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_D2z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F3x_vrr;
    Double I_TWOBODYOVERLAP_D2x_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_D2x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_D2y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_D2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_D2y_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_D2z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_D2x_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_D2x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_D2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_D2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_D2x_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F3y_vrr;
    Double I_TWOBODYOVERLAP_D2y_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_D2z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_D2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_D2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_D2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_D2x_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_D2x_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F3z_vrr;
    Double I_TWOBODYOVERLAP_D2y_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_D2z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_D2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_D2x_G4y_vrr = PBY*I_TWOBODYOVERLAP_D2x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_G4y_vrr = PBY*I_TWOBODYOVERLAP_D2y_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_G4y_vrr = PBY*I_TWOBODYOVERLAP_D2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_TWOBODYOVERLAP_D2y_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_D2z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F3y_vrr;
    Double I_TWOBODYOVERLAP_D2x_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_D2x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_D2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_D2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_TWOBODYOVERLAP_D2y_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_D2y_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F3z_vrr;
    Double I_TWOBODYOVERLAP_D2z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_D2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_D2x_G4z_vrr = PBZ*I_TWOBODYOVERLAP_D2x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_G4z_vrr = PBZ*I_TWOBODYOVERLAP_D2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_D2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_G
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_G4x_vrr = PBX*I_TWOBODYOVERLAP_F3x_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G4x_vrr = PBX*I_TWOBODYOVERLAP_F2xy_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G4x_vrr = PBX*I_TWOBODYOVERLAP_F2xz_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G4x_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G4x_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G4x_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3y_G4x_vrr = PBX*I_TWOBODYOVERLAP_F3y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G4x_vrr = PBX*I_TWOBODYOVERLAP_F2yz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G4x_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_G4x_vrr = PBX*I_TWOBODYOVERLAP_F3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_G3xy_vrr = PBY*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G3xy_vrr = PBY*I_TWOBODYOVERLAP_F2xy_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G3xy_vrr = PBY*I_TWOBODYOVERLAP_F2xz_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3y_G3xy_vrr = PBY*I_TWOBODYOVERLAP_F3y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G3xy_vrr = PBY*I_TWOBODYOVERLAP_F2yz_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3x_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3y_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3x_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_F3x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_F2xy_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_F2xz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3y_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_F3y_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_F2yz_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_F3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_F3y_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_F3z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_F3x_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_F3x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3y_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_F3y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_F3z_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_F3x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_F2xy_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_F3y_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_F2xz_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3y_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_F2yz_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_F2xy_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_F2xz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_F3y_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_F2yz_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F3z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F3x_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_F3x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_F2xy_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_F2xz_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_F3y_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_F2yz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_F3z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_F3x_G4y_vrr = PBY*I_TWOBODYOVERLAP_F3x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G4y_vrr = PBY*I_TWOBODYOVERLAP_F2xy_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G4y_vrr = PBY*I_TWOBODYOVERLAP_F2xz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G4y_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G4y_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G4y_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_G4y_vrr = PBY*I_TWOBODYOVERLAP_F3y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G4y_vrr = PBY*I_TWOBODYOVERLAP_F2yz_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G4y_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_G4y_vrr = PBY*I_TWOBODYOVERLAP_F3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3y_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3x_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_F3x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_F3y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_F3z_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_F2xy_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_F2xz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_F3y_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_F3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_F2yz_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_F3z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_F3x_G4z_vrr = PBZ*I_TWOBODYOVERLAP_F3x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_G4z_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_G4z_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fxyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_G4z_vrr = PBZ*I_TWOBODYOVERLAP_F3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_G4z_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_F3z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_G_G
     * expanding position: BRA1
     * code section is: VRR
     * totally 45 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_G
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_G
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_G4x_vrr = PAX*I_TWOBODYOVERLAP_F3x_G4x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G4x_vrr = PAY*I_TWOBODYOVERLAP_F3x_G4x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G4x_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G4x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G4x_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G4x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G4x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G4x_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G4x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G4x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G4x_vrr = PAX*I_TWOBODYOVERLAP_F3y_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G4x_vrr = PAX*I_TWOBODYOVERLAP_F3z_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_G4y_G4x_vrr = PAY*I_TWOBODYOVERLAP_F3y_G4x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G4x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G4x_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G4x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G4x_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G4x_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G4x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G4x_vrr = PAY*I_TWOBODYOVERLAP_F3z_G4x_vrr;
    Double I_TWOBODYOVERLAP_G4z_G4x_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G4x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G4x_vrr;
    Double I_TWOBODYOVERLAP_G4x_G3xy_vrr = PAX*I_TWOBODYOVERLAP_F3x_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G3xy_vrr = PAY*I_TWOBODYOVERLAP_F3x_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G3xy_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G3xy_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G3xy_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G3xy_vrr = PAX*I_TWOBODYOVERLAP_F3y_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G3xy_vrr = PAX*I_TWOBODYOVERLAP_F3z_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G4y_G3xy_vrr = PAY*I_TWOBODYOVERLAP_F3y_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G3xy_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G3xy_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G3xy_vrr = PAY*I_TWOBODYOVERLAP_F3z_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_G4z_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G3xy_vrr;
    Double I_TWOBODYOVERLAP_G4x_G3xz_vrr = PAX*I_TWOBODYOVERLAP_F3x_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G3xz_vrr = PAY*I_TWOBODYOVERLAP_F3x_G3xz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G3xz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G3xz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_F3x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G3xz_vrr = PAX*I_TWOBODYOVERLAP_F3y_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G3xz_vrr = PAX*I_TWOBODYOVERLAP_F3z_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G4y_G3xz_vrr = PAY*I_TWOBODYOVERLAP_F3y_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G3xz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_F3x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G3xz_vrr = PAY*I_TWOBODYOVERLAP_F3z_G3xz_vrr;
    Double I_TWOBODYOVERLAP_G4z_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_G4x_G2x2y_vrr = PAX*I_TWOBODYOVERLAP_F3x_G2x2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_F3x_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G2x2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G2x2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G2x2y_vrr = PAX*I_TWOBODYOVERLAP_F3y_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G2x2y_vrr = PAX*I_TWOBODYOVERLAP_F3z_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G4y_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_F3y_G2x2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G2x2y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_F3z_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G4z_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G2x2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_G2xyz_vrr = PAX*I_TWOBODYOVERLAP_F3x_G2xyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G2xyz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_F3x_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G2xyz_vrr = PAX*I_TWOBODYOVERLAP_F3y_G2xyz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G2xyz_vrr = PAX*I_TWOBODYOVERLAP_F3z_G2xyz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_G4y_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_F3y_G2xyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_F3z_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G4z_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G2xyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G4x_G2x2z_vrr = PAX*I_TWOBODYOVERLAP_F3x_G2x2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_F3x_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G2x2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G2x2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G2x2z_vrr = PAX*I_TWOBODYOVERLAP_F3y_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G2x2z_vrr = PAX*I_TWOBODYOVERLAP_F3z_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_F3y_G2x2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G2x2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_F3z_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G2x2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G4x_Gx3y_vrr = PAX*I_TWOBODYOVERLAP_F3x_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_F3x_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Gx3y_vrr = PAX*I_TWOBODYOVERLAP_F3y_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Gx3y_vrr = PAX*I_TWOBODYOVERLAP_F3z_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_G4y_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_F3y_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_F3z_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G4z_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Gx2yz_vrr = PAX*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Gx2yz_vrr = PAX*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Gx2yz_vrr = PAX*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G4y_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_G4z_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Gxy2z_vrr = PAX*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Gxy2z_vrr = PAX*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Gxy2z_vrr = PAX*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_G4x_Gx3z_vrr = PAX*I_TWOBODYOVERLAP_F3x_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_F3x_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Gx3z_vrr = PAX*I_TWOBODYOVERLAP_F3y_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Gx3z_vrr = PAX*I_TWOBODYOVERLAP_F3z_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_F3y_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_F3z_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G4x_G4y_vrr = PAX*I_TWOBODYOVERLAP_F3x_G4y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G4y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G4y_vrr = PAY*I_TWOBODYOVERLAP_F3x_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G4y_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G4y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G4y_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G4y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_F2xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G4y_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G4y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G4y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G4y_vrr = PAX*I_TWOBODYOVERLAP_F3y_G4y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G4y_vrr = PAX*I_TWOBODYOVERLAP_F3z_G4y_vrr;
    Double I_TWOBODYOVERLAP_G4y_G4y_vrr = PAY*I_TWOBODYOVERLAP_F3y_G4y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G4y_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G4y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G4y_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G4y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G4y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G4y_vrr = PAY*I_TWOBODYOVERLAP_F3z_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_G4z_G4y_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G4y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G4y_vrr;
    Double I_TWOBODYOVERLAP_G4x_G3yz_vrr = PAX*I_TWOBODYOVERLAP_F3x_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G3yz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G3yz_vrr = PAY*I_TWOBODYOVERLAP_F3x_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G3yz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G3yz_vrr = PAX*I_TWOBODYOVERLAP_F3y_G3yz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G3yz_vrr = PAX*I_TWOBODYOVERLAP_F3z_G3yz_vrr;
    Double I_TWOBODYOVERLAP_G4y_G3yz_vrr = PAY*I_TWOBODYOVERLAP_F3y_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G3yz_vrr = PAY*I_TWOBODYOVERLAP_F3z_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G4z_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_G4x_G2y2z_vrr = PAX*I_TWOBODYOVERLAP_F3x_G2y2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G2y2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_F3x_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G2y2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G2y2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G2y2z_vrr = PAX*I_TWOBODYOVERLAP_F3y_G2y2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G2y2z_vrr = PAX*I_TWOBODYOVERLAP_F3z_G2y2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_F3y_G2y2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G2y2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_F3z_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G2y2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G4x_Gy3z_vrr = PAX*I_TWOBODYOVERLAP_F3x_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Gy3z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_F3x_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_F3z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Gy3z_vrr = PAX*I_TWOBODYOVERLAP_F3y_Gy3z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Gy3z_vrr = PAX*I_TWOBODYOVERLAP_F3z_Gy3z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_F3y_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_F3z_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G4x_G4z_vrr = PAX*I_TWOBODYOVERLAP_F3x_G4z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_G4z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_G4z_vrr = PAY*I_TWOBODYOVERLAP_F3x_G4z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_G4z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_G4z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_G4z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G4z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_G4z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_G4z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_F2xz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_G4z_vrr = PAX*I_TWOBODYOVERLAP_F3y_G4z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_G4z_vrr = PAX*I_TWOBODYOVERLAP_F3z_G4z_vrr;
    Double I_TWOBODYOVERLAP_G4y_G4z_vrr = PAY*I_TWOBODYOVERLAP_F3y_G4z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_G4z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_G4z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_G4z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_G4z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_F2yz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_G4z_vrr = PAY*I_TWOBODYOVERLAP_F3z_G4z_vrr;
    Double I_TWOBODYOVERLAP_G4z_G4z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_G4z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_G
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_G
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_G
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_G4x_vrr = PAX*I_TWOBODYOVERLAP_G4x_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G4x_vrr = PAY*I_TWOBODYOVERLAP_G4x_G4x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G4x_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G4x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G4x_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G4x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G4x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G4x_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G4x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G4x_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G4x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G4x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G4x_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G4x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G4x_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G4x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G4x_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G4x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G4x_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G4x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G4x_vrr = PAX*I_TWOBODYOVERLAP_G4y_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G4x_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G4x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G4x_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_G2y2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G4x_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G4x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G4x_vrr = PAX*I_TWOBODYOVERLAP_G4z_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5y_G4x_vrr = PAY*I_TWOBODYOVERLAP_G4y_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G4x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G4x_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G4x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G4x_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G4x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G4x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G4x_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G4x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G4x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G4x_vrr = PAY*I_TWOBODYOVERLAP_G4z_G4x_vrr;
    Double I_TWOBODYOVERLAP_H5z_G4x_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G4x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G4x_vrr;
    Double I_TWOBODYOVERLAP_H5x_G3xy_vrr = PAX*I_TWOBODYOVERLAP_G4x_G3xy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G3xy_vrr = PAY*I_TWOBODYOVERLAP_G4x_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G3xy_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G3xy_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G3xy_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G3xy_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G3xy_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G3xy_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G3xy_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G3xy_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G3xy_vrr = PAX*I_TWOBODYOVERLAP_G4y_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G3xy_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G3xy_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G3xy_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G3xy_vrr = PAX*I_TWOBODYOVERLAP_G4z_G3xy_vrr+3*oned2z*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H5y_G3xy_vrr = PAY*I_TWOBODYOVERLAP_G4y_G3xy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G3xy_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G3xy_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G3xy_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G3xy_vrr = PAY*I_TWOBODYOVERLAP_G4z_G3xy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5z_G3xy_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G3xy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G3xy_vrr;
    Double I_TWOBODYOVERLAP_H5x_G3xz_vrr = PAX*I_TWOBODYOVERLAP_G4x_G3xz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G3xz_vrr = PAY*I_TWOBODYOVERLAP_G4x_G3xz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G3xz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G3xz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G3xz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G3xz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G3xz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G3xz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G3xz_vrr = PAX*I_TWOBODYOVERLAP_G4y_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_G4y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G3xz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G3xz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G3xz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G3xz_vrr = PAX*I_TWOBODYOVERLAP_G4z_G3xz_vrr+3*oned2z*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H5y_G3xz_vrr = PAY*I_TWOBODYOVERLAP_G4y_G3xz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G3xz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G3xz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G3xz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G3xz_vrr = PAY*I_TWOBODYOVERLAP_G4z_G3xz_vrr;
    Double I_TWOBODYOVERLAP_H5z_G3xz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G3xz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G3xz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5x_G2x2y_vrr = PAX*I_TWOBODYOVERLAP_G4x_G2x2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_G4x_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G2x2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G2x2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G2x2y_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G2x2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G2x2y_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G2x2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G2x2y_vrr = PAX*I_TWOBODYOVERLAP_G4y_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G2x2y_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G2x2y_vrr = PAX*I_TWOBODYOVERLAP_G4z_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_G4y_G2x2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G2x2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G2x2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G2x2y_vrr = PAY*I_TWOBODYOVERLAP_G4z_G2x2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H5z_G2x2y_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G2x2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G2x2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_G2xyz_vrr = PAX*I_TWOBODYOVERLAP_G4x_G2xyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G2xyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_G4x_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G2xyz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G2xyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G2xyz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G2xyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G2xyz_vrr = PAX*I_TWOBODYOVERLAP_G4y_G2xyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G2xyz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G2xyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G2xyz_vrr = PAX*I_TWOBODYOVERLAP_G4z_G2xyz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H5y_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_G4y_G2xyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G2xyz_vrr = PAY*I_TWOBODYOVERLAP_G4z_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H5z_G2xyz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G2xyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G2xyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H5x_G2x2z_vrr = PAX*I_TWOBODYOVERLAP_G4x_G2x2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_G4x_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G2x2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G2x2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G2x2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G2x2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G2x2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G2x2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G2x2z_vrr = PAX*I_TWOBODYOVERLAP_G4y_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G2x2z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G2x2z_vrr = PAX*I_TWOBODYOVERLAP_G4z_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_G4y_G2x2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G2x2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G2x2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G2x2z_vrr = PAY*I_TWOBODYOVERLAP_G4z_G2x2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_G2x2z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G2x2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G2x2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gx3y_vrr = PAX*I_TWOBODYOVERLAP_G4x_Gx3y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_G4x_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gx3y_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gx3y_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gx3y_vrr = PAX*I_TWOBODYOVERLAP_G4y_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gx3y_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gx3y_vrr = PAX*I_TWOBODYOVERLAP_G4z_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_G4y_Gx3y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Gx3y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gx3y_vrr = PAY*I_TWOBODYOVERLAP_G4z_Gx3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gx3y_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Gx3y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Gx3y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gx2yz_vrr = PAX*I_TWOBODYOVERLAP_G4x_Gx2yz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_G4x_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gx2yz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gx2yz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gx2yz_vrr = PAX*I_TWOBODYOVERLAP_G4y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gx2yz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gx2yz_vrr = PAX*I_TWOBODYOVERLAP_G4z_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_G4y_Gx2yz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gx2yz_vrr = PAY*I_TWOBODYOVERLAP_G4z_Gx2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gx2yz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Gx2yz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gxy2z_vrr = PAX*I_TWOBODYOVERLAP_G4x_Gxy2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_G4x_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gxy2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gxy2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gxy2z_vrr = PAX*I_TWOBODYOVERLAP_G4y_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gxy2z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gxy2z_vrr = PAX*I_TWOBODYOVERLAP_G4z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_G4y_Gxy2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gxy2z_vrr = PAY*I_TWOBODYOVERLAP_G4z_Gxy2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gxy2z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Gxy2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Fxyz_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gx3z_vrr = PAX*I_TWOBODYOVERLAP_G4x_Gx3z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_G4x_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G3xz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gx3z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gx3z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gx3z_vrr = PAX*I_TWOBODYOVERLAP_G4y_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gx3z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gx3z_vrr = PAX*I_TWOBODYOVERLAP_G4z_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_G4y_Gx3z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G3yz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Gx3z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gx3z_vrr = PAY*I_TWOBODYOVERLAP_G4z_Gx3z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gx3z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Gx3z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Gx3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H5x_G4y_vrr = PAX*I_TWOBODYOVERLAP_G4x_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G4y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G4y_vrr = PAY*I_TWOBODYOVERLAP_G4x_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G4y_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G4y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G4y_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G4y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G4y_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G4y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G4y_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G4y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G4y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G4y_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G4y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G4y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G4y_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G4y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G4y_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_G2x2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G4y_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G4y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G4y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G4y_vrr = PAX*I_TWOBODYOVERLAP_G4y_G4y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G4y_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G4y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G4y_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G4y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G4y_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G4y_vrr = PAX*I_TWOBODYOVERLAP_G4z_G4y_vrr;
    Double I_TWOBODYOVERLAP_H5y_G4y_vrr = PAY*I_TWOBODYOVERLAP_G4y_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G4y_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G4y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G4y_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G4y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G4y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G4y_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G4y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G4y_vrr = PAY*I_TWOBODYOVERLAP_G4z_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5z_G4y_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G4y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G4y_vrr;
    Double I_TWOBODYOVERLAP_H5x_G3yz_vrr = PAX*I_TWOBODYOVERLAP_G4x_G3yz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G3yz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G3yz_vrr = PAY*I_TWOBODYOVERLAP_G4x_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_G4x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G3yz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_G3xy_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G3yz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G3yz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G3yz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G3yz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G3yz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G3yz_vrr = PAX*I_TWOBODYOVERLAP_G4y_G3yz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G3yz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G3yz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G3yz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G3yz_vrr = PAX*I_TWOBODYOVERLAP_G4z_G3yz_vrr;
    Double I_TWOBODYOVERLAP_H5y_G3yz_vrr = PAY*I_TWOBODYOVERLAP_G4y_G3yz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G3yz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_Gy3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G3yz_vrr = PAY*I_TWOBODYOVERLAP_G4z_G3yz_vrr+3*oned2z*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H5z_G3yz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G3yz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G3yz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5x_G2y2z_vrr = PAX*I_TWOBODYOVERLAP_G4x_G2y2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G2y2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_G4x_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G2y2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G2y2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G2y2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G2y2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G2y2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G2y2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G2y2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G2y2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G2y2z_vrr = PAX*I_TWOBODYOVERLAP_G4y_G2y2z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G2y2z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G2y2z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G2y2z_vrr = PAX*I_TWOBODYOVERLAP_G4z_G2y2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_G4y_G2y2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G2y2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G2y2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G2y2z_vrr = PAY*I_TWOBODYOVERLAP_G4z_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_G2y2z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G2y2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G2y2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gy3z_vrr = PAX*I_TWOBODYOVERLAP_G4x_Gy3z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Gy3z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_G4x_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G4x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G3xz_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gy3z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gy3z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gy3z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gy3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gy3z_vrr = PAX*I_TWOBODYOVERLAP_G4y_Gy3z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gy3z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Gy3z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gy3z_vrr = PAX*I_TWOBODYOVERLAP_G4z_Gy3z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_G4y_Gy3z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G3yz_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gy3z_vrr = PAY*I_TWOBODYOVERLAP_G4z_Gy3z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gy3z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Gy3z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Gy3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_H5x_G4z_vrr = PAX*I_TWOBODYOVERLAP_G4x_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_G4z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G4z_vrr = PAY*I_TWOBODYOVERLAP_G4x_G4z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G4z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G4z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_G4z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G4z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G4z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G4z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_G4z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G4z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_G4z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G4z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G4z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_G2x2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G4z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_G4z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G4z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_G4z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G4z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G4z_vrr = PAX*I_TWOBODYOVERLAP_G4y_G4z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G4z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G4z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_G4z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G4z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_G4z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G4z_vrr = PAX*I_TWOBODYOVERLAP_G4z_G4z_vrr;
    Double I_TWOBODYOVERLAP_H5y_G4z_vrr = PAY*I_TWOBODYOVERLAP_G4y_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_G4z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G4z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G4z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_G4z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G4z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_G4z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_G4z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G4z_vrr = PAY*I_TWOBODYOVERLAP_G4z_G4z_vrr;
    Double I_TWOBODYOVERLAP_H5z_G4z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_G4z_vrr+4*oned2z*I_TWOBODYOVERLAP_G4z_F3z_vrr;

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
     * expanding position: BRA1
     * code section is: VRR
     * totally 8 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_S_D
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     ************************************************************/
    Double I_KINETIC_D2x_D2x_vrr = PAX*I_KINETIC_Px_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+2*oned2z*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_Dxy_D2x_vrr = PAY*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_KINETIC_Dxz_D2x_vrr = PAZ*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_D2x_vrr;
    Double I_KINETIC_D2y_D2x_vrr = PAY*I_KINETIC_Py_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_Dyz_D2x_vrr = PAZ*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_D2x_vrr;
    Double I_KINETIC_D2z_D2x_vrr = PAZ*I_KINETIC_Pz_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_D2x_Dxy_vrr = PAX*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2y_Dxy_vrr = PAY*I_KINETIC_Py_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2z_Dxy_vrr = PAZ*I_KINETIC_Pz_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2x_Dxz_vrr = PAX*I_KINETIC_Px_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_Dxy_Dxz_vrr = PAY*I_KINETIC_Px_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Dxz_vrr;
    Double I_KINETIC_D2y_Dxz_vrr = PAY*I_KINETIC_Py_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_D2z_Dxz_vrr = PAZ*I_KINETIC_Pz_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_D2x_D2y_vrr = PAX*I_KINETIC_Px_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_Dxy_D2y_vrr = PAY*I_KINETIC_Px_D2y_vrr+2*oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_KINETIC_Dxz_D2y_vrr = PAZ*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_D2y_vrr;
    Double I_KINETIC_D2y_D2y_vrr = PAY*I_KINETIC_Py_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+2*oned2z*I_KINETIC_Py_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_Dyz_D2y_vrr = PAZ*I_KINETIC_Py_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_D2y_vrr;
    Double I_KINETIC_D2z_D2y_vrr = PAZ*I_KINETIC_Pz_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_D2x_Dyz_vrr = PAX*I_KINETIC_Px_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2y_Dyz_vrr = PAY*I_KINETIC_Py_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_Py_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2z_Dyz_vrr = PAZ*I_KINETIC_Pz_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_Pz_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2x_D2z_vrr = PAX*I_KINETIC_Px_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_KINETIC_Dxy_D2z_vrr = PAY*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2z_vrr;
    Double I_KINETIC_Dxz_D2z_vrr = PAZ*I_KINETIC_Px_D2z_vrr+2*oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_D2z_vrr;
    Double I_KINETIC_D2y_D2z_vrr = PAY*I_KINETIC_Py_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_KINETIC_Dyz_D2z_vrr = PAZ*I_KINETIC_Py_D2z_vrr+2*oned2z*I_KINETIC_Py_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_D2z_vrr;
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
     * shell quartet name: SQ_KINETIC_D_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_KINETIC_D2x_F3x_vrr = PBX*I_KINETIC_D2x_D2x_vrr+2*oned2z*I_KINETIC_Px_D2x_vrr+2*oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3x_vrr-2*adz*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_Dxy_F3x_vrr = PBX*I_KINETIC_Dxy_D2x_vrr+oned2z*I_KINETIC_Py_D2x_vrr+2*oned2z*I_KINETIC_Dxy_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_KINETIC_Dxz_F3x_vrr = PBX*I_KINETIC_Dxz_D2x_vrr+oned2z*I_KINETIC_Pz_D2x_vrr+2*oned2z*I_KINETIC_Dxz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_KINETIC_D2y_F3x_vrr = PBX*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_Dyz_F3x_vrr = PBX*I_KINETIC_Dyz_D2x_vrr+2*oned2z*I_KINETIC_Dyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_KINETIC_D2z_F3x_vrr = PBX*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_KINETIC_D2x_F2xy_vrr = PBY*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_KINETIC_Dxy_F2xy_vrr = PBY*I_KINETIC_Dxy_D2x_vrr+oned2z*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F2xy_vrr;
    Double I_KINETIC_Dxz_F2xy_vrr = PBY*I_KINETIC_Dxz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F2xy_vrr;
    Double I_KINETIC_D2y_F2xy_vrr = PBY*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_KINETIC_Dyz_F2xy_vrr = PBY*I_KINETIC_Dyz_D2x_vrr+oned2z*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F2xy_vrr;
    Double I_KINETIC_D2z_F2xy_vrr = PBY*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2xy_vrr;
    Double I_KINETIC_D2x_F2xz_vrr = PBZ*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_KINETIC_Dxy_F2xz_vrr = PBZ*I_KINETIC_Dxy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F2xz_vrr;
    Double I_KINETIC_Dxz_F2xz_vrr = PBZ*I_KINETIC_Dxz_D2x_vrr+oned2z*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F2xz_vrr;
    Double I_KINETIC_D2y_F2xz_vrr = PBZ*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2xz_vrr;
    Double I_KINETIC_Dyz_F2xz_vrr = PBZ*I_KINETIC_Dyz_D2x_vrr+oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F2xz_vrr;
    Double I_KINETIC_D2z_F2xz_vrr = PBZ*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2xz_vrr;
    Double I_KINETIC_D2x_Fx2y_vrr = PBX*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_KINETIC_Dxy_Fx2y_vrr = PBX*I_KINETIC_Dxy_D2y_vrr+oned2z*I_KINETIC_Py_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fx2y_vrr;
    Double I_KINETIC_Dxz_Fx2y_vrr = PBX*I_KINETIC_Dxz_D2y_vrr+oned2z*I_KINETIC_Pz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Fx2y_vrr;
    Double I_KINETIC_D2y_Fx2y_vrr = PBX*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_KINETIC_Dyz_Fx2y_vrr = PBX*I_KINETIC_Dyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Fx2y_vrr;
    Double I_KINETIC_D2z_Fx2y_vrr = PBX*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fx2y_vrr;
    Double I_KINETIC_D2x_Fxyz_vrr = PBZ*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fxyz_vrr;
    Double I_KINETIC_D2y_Fxyz_vrr = PBZ*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fxyz_vrr;
    Double I_KINETIC_D2z_Fxyz_vrr = PBZ*I_KINETIC_D2z_Dxy_vrr+2*oned2z*I_KINETIC_Pz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fxyz_vrr;
    Double I_KINETIC_D2x_Fx2z_vrr = PBX*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_KINETIC_Dxy_Fx2z_vrr = PBX*I_KINETIC_Dxy_D2z_vrr+oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fx2z_vrr;
    Double I_KINETIC_Dxz_Fx2z_vrr = PBX*I_KINETIC_Dxz_D2z_vrr+oned2z*I_KINETIC_Pz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Fx2z_vrr;
    Double I_KINETIC_D2y_Fx2z_vrr = PBX*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_KINETIC_Dyz_Fx2z_vrr = PBX*I_KINETIC_Dyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Fx2z_vrr;
    Double I_KINETIC_D2z_Fx2z_vrr = PBX*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fx2z_vrr;
    Double I_KINETIC_D2x_F3y_vrr = PBY*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3y_vrr-2*adz*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_Dxy_F3y_vrr = PBY*I_KINETIC_Dxy_D2y_vrr+oned2z*I_KINETIC_Px_D2y_vrr+2*oned2z*I_KINETIC_Dxy_Py_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_KINETIC_Dxz_F3y_vrr = PBY*I_KINETIC_Dxz_D2y_vrr+2*oned2z*I_KINETIC_Dxz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Dxz_Py_vrr;
    Double I_KINETIC_D2y_F3y_vrr = PBY*I_KINETIC_D2y_D2y_vrr+2*oned2z*I_KINETIC_Py_D2y_vrr+2*oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_Dyz_F3y_vrr = PBY*I_KINETIC_Dyz_D2y_vrr+oned2z*I_KINETIC_Pz_D2y_vrr+2*oned2z*I_KINETIC_Dyz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Dyz_Py_vrr;
    Double I_KINETIC_D2z_F3y_vrr = PBY*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_KINETIC_D2x_F2yz_vrr = PBZ*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_KINETIC_Dxy_F2yz_vrr = PBZ*I_KINETIC_Dxy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F2yz_vrr;
    Double I_KINETIC_Dxz_F2yz_vrr = PBZ*I_KINETIC_Dxz_D2y_vrr+oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F2yz_vrr;
    Double I_KINETIC_D2y_F2yz_vrr = PBZ*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2yz_vrr;
    Double I_KINETIC_Dyz_F2yz_vrr = PBZ*I_KINETIC_Dyz_D2y_vrr+oned2z*I_KINETIC_Py_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F2yz_vrr;
    Double I_KINETIC_D2z_F2yz_vrr = PBZ*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_Pz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2yz_vrr;
    Double I_KINETIC_D2x_Fy2z_vrr = PBY*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_KINETIC_D2y_Fy2z_vrr = PBY*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fy2z_vrr;
    Double I_KINETIC_D2z_Fy2z_vrr = PBY*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fy2z_vrr;
    Double I_KINETIC_D2x_F3z_vrr = PBZ*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3z_vrr-2*adz*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_Dxy_F3z_vrr = PBZ*I_KINETIC_Dxy_D2z_vrr+2*oned2z*I_KINETIC_Dxy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_KINETIC_Dxz_F3z_vrr = PBZ*I_KINETIC_Dxz_D2z_vrr+oned2z*I_KINETIC_Px_D2z_vrr+2*oned2z*I_KINETIC_Dxz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Dxz_Pz_vrr;
    Double I_KINETIC_D2y_F3z_vrr = PBZ*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_Dyz_F3z_vrr = PBZ*I_KINETIC_Dyz_D2z_vrr+oned2z*I_KINETIC_Py_D2z_vrr+2*oned2z*I_KINETIC_Dyz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Dyz_Pz_vrr;
    Double I_KINETIC_D2z_F3z_vrr = PBZ*I_KINETIC_D2z_D2z_vrr+2*oned2z*I_KINETIC_Pz_D2z_vrr+2*oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
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
     * expanding position: BRA1
     * code section is: VRR
     * totally 8 integrals are omitted 
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
    Double I_KINETIC_F3y_Fxyz_vrr = PAY*I_KINETIC_D2y_Fxyz_vrr+2*oned2z*I_KINETIC_Py_Fxyz_vrr+oned2z*I_KINETIC_D2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fxyz_vrr;
    Double I_KINETIC_F2yz_Fxyz_vrr = PAZ*I_KINETIC_D2y_Fxyz_vrr+oned2z*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fxyz_vrr;
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
    Double I_KINETIC_F3y_Fy2z_vrr = PAY*I_KINETIC_D2y_Fy2z_vrr+2*oned2z*I_KINETIC_Py_Fy2z_vrr+oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fy2z_vrr;
    Double I_KINETIC_F2yz_Fy2z_vrr = PAZ*I_KINETIC_D2y_Fy2z_vrr+2*oned2z*I_KINETIC_D2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fy2z_vrr;
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
     * shell quartet name: SQ_KINETIC_D_G
     * expanding position: BRA2
     * code section is: VRR
     * totally 45 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_F
     * RHS shell quartet name: SQ_KINETIC_P_F
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_G
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     ************************************************************/
    Double I_KINETIC_D2x_G4x_vrr = PBX*I_KINETIC_D2x_F3x_vrr+2*oned2z*I_KINETIC_Px_F3x_vrr+3*oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G4x_vrr-3*adz*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_KINETIC_D2y_G4x_vrr = PBX*I_KINETIC_D2y_F3x_vrr+3*oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G4x_vrr-3*adz*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_KINETIC_D2z_G4x_vrr = PBX*I_KINETIC_D2z_F3x_vrr+3*oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_KINETIC_D2x_G3xy_vrr = PBY*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G3xy_vrr;
    Double I_KINETIC_D2y_G3xy_vrr = PBY*I_KINETIC_D2y_F3x_vrr+2*oned2z*I_KINETIC_Py_F3x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G3xy_vrr;
    Double I_KINETIC_D2z_G3xy_vrr = PBY*I_KINETIC_D2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G3xy_vrr;
    Double I_KINETIC_D2x_G3xz_vrr = PBZ*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G3xz_vrr;
    Double I_KINETIC_D2y_G3xz_vrr = PBZ*I_KINETIC_D2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G3xz_vrr;
    Double I_KINETIC_D2z_G3xz_vrr = PBZ*I_KINETIC_D2z_F3x_vrr+2*oned2z*I_KINETIC_Pz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G3xz_vrr;
    Double I_KINETIC_D2x_G2x2y_vrr = PBY*I_KINETIC_D2x_F2xy_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G2x2y_vrr-adz*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_KINETIC_D2y_G2x2y_vrr = PBY*I_KINETIC_D2y_F2xy_vrr+2*oned2z*I_KINETIC_Py_F2xy_vrr+oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G2x2y_vrr-adz*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_KINETIC_D2z_G2x2y_vrr = PBY*I_KINETIC_D2z_F2xy_vrr+oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_KINETIC_D2x_G2xyz_vrr = PBZ*I_KINETIC_D2x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G2xyz_vrr;
    Double I_KINETIC_D2y_G2xyz_vrr = PBZ*I_KINETIC_D2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G2xyz_vrr;
    Double I_KINETIC_D2z_G2xyz_vrr = PBZ*I_KINETIC_D2z_F2xy_vrr+2*oned2z*I_KINETIC_Pz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G2xyz_vrr;
    Double I_KINETIC_D2x_G2x2z_vrr = PBZ*I_KINETIC_D2x_F2xz_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G2x2z_vrr-adz*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_KINETIC_D2y_G2x2z_vrr = PBZ*I_KINETIC_D2y_F2xz_vrr+oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G2x2z_vrr-adz*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_KINETIC_D2z_G2x2z_vrr = PBZ*I_KINETIC_D2z_F2xz_vrr+2*oned2z*I_KINETIC_Pz_F2xz_vrr+oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_KINETIC_D2x_Gx3y_vrr = PBX*I_KINETIC_D2x_F3y_vrr+2*oned2z*I_KINETIC_Px_F3y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Gx3y_vrr;
    Double I_KINETIC_D2y_Gx3y_vrr = PBX*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Gx3y_vrr;
    Double I_KINETIC_D2z_Gx3y_vrr = PBX*I_KINETIC_D2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Gx3y_vrr;
    Double I_KINETIC_D2x_Gx2yz_vrr = PBZ*I_KINETIC_D2x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Gx2yz_vrr;
    Double I_KINETIC_D2y_Gx2yz_vrr = PBZ*I_KINETIC_D2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Gx2yz_vrr;
    Double I_KINETIC_D2z_Gx2yz_vrr = PBZ*I_KINETIC_D2z_Fx2y_vrr+2*oned2z*I_KINETIC_Pz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Gx2yz_vrr;
    Double I_KINETIC_D2x_Gxy2z_vrr = PBY*I_KINETIC_D2x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Gxy2z_vrr;
    Double I_KINETIC_D2y_Gxy2z_vrr = PBY*I_KINETIC_D2y_Fx2z_vrr+2*oned2z*I_KINETIC_Py_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Gxy2z_vrr;
    Double I_KINETIC_D2z_Gxy2z_vrr = PBY*I_KINETIC_D2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Gxy2z_vrr;
    Double I_KINETIC_D2x_Gx3z_vrr = PBX*I_KINETIC_D2x_F3z_vrr+2*oned2z*I_KINETIC_Px_F3z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Gx3z_vrr;
    Double I_KINETIC_D2y_Gx3z_vrr = PBX*I_KINETIC_D2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Gx3z_vrr;
    Double I_KINETIC_D2z_Gx3z_vrr = PBX*I_KINETIC_D2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Gx3z_vrr;
    Double I_KINETIC_D2x_G4y_vrr = PBY*I_KINETIC_D2x_F3y_vrr+3*oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G4y_vrr-3*adz*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_KINETIC_D2y_G4y_vrr = PBY*I_KINETIC_D2y_F3y_vrr+2*oned2z*I_KINETIC_Py_F3y_vrr+3*oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G4y_vrr-3*adz*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_KINETIC_D2z_G4y_vrr = PBY*I_KINETIC_D2z_F3y_vrr+3*oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_KINETIC_D2x_G3yz_vrr = PBZ*I_KINETIC_D2x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G3yz_vrr;
    Double I_KINETIC_D2y_G3yz_vrr = PBZ*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G3yz_vrr;
    Double I_KINETIC_D2z_G3yz_vrr = PBZ*I_KINETIC_D2z_F3y_vrr+2*oned2z*I_KINETIC_Pz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G3yz_vrr;
    Double I_KINETIC_D2x_G2y2z_vrr = PBZ*I_KINETIC_D2x_F2yz_vrr+oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G2y2z_vrr-adz*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_KINETIC_D2y_G2y2z_vrr = PBZ*I_KINETIC_D2y_F2yz_vrr+oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G2y2z_vrr-adz*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_KINETIC_D2z_G2y2z_vrr = PBZ*I_KINETIC_D2z_F2yz_vrr+2*oned2z*I_KINETIC_Pz_F2yz_vrr+oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_KINETIC_D2x_Gy3z_vrr = PBY*I_KINETIC_D2x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Gy3z_vrr;
    Double I_KINETIC_D2y_Gy3z_vrr = PBY*I_KINETIC_D2y_F3z_vrr+2*oned2z*I_KINETIC_Py_F3z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Gy3z_vrr;
    Double I_KINETIC_D2z_Gy3z_vrr = PBY*I_KINETIC_D2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Gy3z_vrr;
    Double I_KINETIC_D2x_G4z_vrr = PBZ*I_KINETIC_D2x_F3z_vrr+3*oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_G4z_vrr-3*adz*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_KINETIC_D2y_G4z_vrr = PBZ*I_KINETIC_D2y_F3z_vrr+3*oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_G4z_vrr-3*adz*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_KINETIC_D2z_G4z_vrr = PBZ*I_KINETIC_D2z_F3z_vrr+2*oned2z*I_KINETIC_Pz_F3z_vrr+3*oned2z*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_D2z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_G
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_F
     * RHS shell quartet name: SQ_KINETIC_D_F
     * RHS shell quartet name: SQ_KINETIC_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_G
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     ************************************************************/
    Double I_KINETIC_F3x_G4x_vrr = PBX*I_KINETIC_F3x_F3x_vrr+3*oned2z*I_KINETIC_D2x_F3x_vrr+3*oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G4x_vrr-3*adz*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_KINETIC_F2xy_G4x_vrr = PBX*I_KINETIC_F2xy_F3x_vrr+2*oned2z*I_KINETIC_Dxy_F3x_vrr+3*oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G4x_vrr-3*adz*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_KINETIC_F2xz_G4x_vrr = PBX*I_KINETIC_F2xz_F3x_vrr+2*oned2z*I_KINETIC_Dxz_F3x_vrr+3*oned2z*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G4x_vrr-3*adz*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_KINETIC_Fx2y_G4x_vrr = PBX*I_KINETIC_Fx2y_F3x_vrr+oned2z*I_KINETIC_D2y_F3x_vrr+3*oned2z*I_KINETIC_Fx2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_KINETIC_Fxyz_G4x_vrr = PBX*I_KINETIC_Fxyz_F3x_vrr+oned2z*I_KINETIC_Dyz_F3x_vrr+3*oned2z*I_KINETIC_Fxyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_KINETIC_Fx2z_G4x_vrr = PBX*I_KINETIC_Fx2z_F3x_vrr+oned2z*I_KINETIC_D2z_F3x_vrr+3*oned2z*I_KINETIC_Fx2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_KINETIC_F3y_G4x_vrr = PBX*I_KINETIC_F3y_F3x_vrr+3*oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G4x_vrr-3*adz*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_KINETIC_F2yz_G4x_vrr = PBX*I_KINETIC_F2yz_F3x_vrr+3*oned2z*I_KINETIC_F2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G4x_vrr-3*adz*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_KINETIC_Fy2z_G4x_vrr = PBX*I_KINETIC_Fy2z_F3x_vrr+3*oned2z*I_KINETIC_Fy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_KINETIC_F3z_G4x_vrr = PBX*I_KINETIC_F3z_F3x_vrr+3*oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_KINETIC_F3x_G3xy_vrr = PBY*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G3xy_vrr;
    Double I_KINETIC_F2xy_G3xy_vrr = PBY*I_KINETIC_F2xy_F3x_vrr+oned2z*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G3xy_vrr;
    Double I_KINETIC_F2xz_G3xy_vrr = PBY*I_KINETIC_F2xz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G3xy_vrr;
    Double I_KINETIC_Fx2y_G3xy_vrr = PBY*I_KINETIC_Fx2y_F3x_vrr+2*oned2z*I_KINETIC_Dxy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G3xy_vrr;
    Double I_KINETIC_Fxyz_G3xy_vrr = PBY*I_KINETIC_Fxyz_F3x_vrr+oned2z*I_KINETIC_Dxz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G3xy_vrr;
    Double I_KINETIC_Fx2z_G3xy_vrr = PBY*I_KINETIC_Fx2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G3xy_vrr;
    Double I_KINETIC_F3y_G3xy_vrr = PBY*I_KINETIC_F3y_F3x_vrr+3*oned2z*I_KINETIC_D2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G3xy_vrr;
    Double I_KINETIC_F2yz_G3xy_vrr = PBY*I_KINETIC_F2yz_F3x_vrr+2*oned2z*I_KINETIC_Dyz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G3xy_vrr;
    Double I_KINETIC_Fy2z_G3xy_vrr = PBY*I_KINETIC_Fy2z_F3x_vrr+oned2z*I_KINETIC_D2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G3xy_vrr;
    Double I_KINETIC_F3z_G3xy_vrr = PBY*I_KINETIC_F3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G3xy_vrr;
    Double I_KINETIC_F3x_G3xz_vrr = PBZ*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G3xz_vrr;
    Double I_KINETIC_F2xy_G3xz_vrr = PBZ*I_KINETIC_F2xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G3xz_vrr;
    Double I_KINETIC_F2xz_G3xz_vrr = PBZ*I_KINETIC_F2xz_F3x_vrr+oned2z*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G3xz_vrr;
    Double I_KINETIC_Fx2y_G3xz_vrr = PBZ*I_KINETIC_Fx2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G3xz_vrr;
    Double I_KINETIC_Fxyz_G3xz_vrr = PBZ*I_KINETIC_Fxyz_F3x_vrr+oned2z*I_KINETIC_Dxy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G3xz_vrr;
    Double I_KINETIC_Fx2z_G3xz_vrr = PBZ*I_KINETIC_Fx2z_F3x_vrr+2*oned2z*I_KINETIC_Dxz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G3xz_vrr;
    Double I_KINETIC_F3y_G3xz_vrr = PBZ*I_KINETIC_F3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G3xz_vrr;
    Double I_KINETIC_F2yz_G3xz_vrr = PBZ*I_KINETIC_F2yz_F3x_vrr+oned2z*I_KINETIC_D2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G3xz_vrr;
    Double I_KINETIC_Fy2z_G3xz_vrr = PBZ*I_KINETIC_Fy2z_F3x_vrr+2*oned2z*I_KINETIC_Dyz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G3xz_vrr;
    Double I_KINETIC_F3z_G3xz_vrr = PBZ*I_KINETIC_F3z_F3x_vrr+3*oned2z*I_KINETIC_D2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G3xz_vrr;
    Double I_KINETIC_F3x_G2x2y_vrr = PBY*I_KINETIC_F3x_F2xy_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G2x2y_vrr-adz*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_KINETIC_F2xy_G2x2y_vrr = PBY*I_KINETIC_F2xy_F2xy_vrr+oned2z*I_KINETIC_D2x_F2xy_vrr+oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G2x2y_vrr-adz*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_KINETIC_F2xz_G2x2y_vrr = PBY*I_KINETIC_F2xz_F2xy_vrr+oned2z*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G2x2y_vrr-adz*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_KINETIC_Fx2y_G2x2y_vrr = PBY*I_KINETIC_Fx2y_F2xy_vrr+2*oned2z*I_KINETIC_Dxy_F2xy_vrr+oned2z*I_KINETIC_Fx2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_KINETIC_Fxyz_G2x2y_vrr = PBY*I_KINETIC_Fxyz_F2xy_vrr+oned2z*I_KINETIC_Dxz_F2xy_vrr+oned2z*I_KINETIC_Fxyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_KINETIC_Fx2z_G2x2y_vrr = PBY*I_KINETIC_Fx2z_F2xy_vrr+oned2z*I_KINETIC_Fx2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_KINETIC_F3y_G2x2y_vrr = PBY*I_KINETIC_F3y_F2xy_vrr+3*oned2z*I_KINETIC_D2y_F2xy_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G2x2y_vrr-adz*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_KINETIC_F2yz_G2x2y_vrr = PBY*I_KINETIC_F2yz_F2xy_vrr+2*oned2z*I_KINETIC_Dyz_F2xy_vrr+oned2z*I_KINETIC_F2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G2x2y_vrr-adz*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_KINETIC_Fy2z_G2x2y_vrr = PBY*I_KINETIC_Fy2z_F2xy_vrr+oned2z*I_KINETIC_D2z_F2xy_vrr+oned2z*I_KINETIC_Fy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_KINETIC_F3z_G2x2y_vrr = PBY*I_KINETIC_F3z_F2xy_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_KINETIC_F3x_G2xyz_vrr = PBZ*I_KINETIC_F3x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G2xyz_vrr;
    Double I_KINETIC_F2xy_G2xyz_vrr = PBZ*I_KINETIC_F2xy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G2xyz_vrr;
    Double I_KINETIC_F2xz_G2xyz_vrr = PBZ*I_KINETIC_F2xz_F2xy_vrr+oned2z*I_KINETIC_D2x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G2xyz_vrr;
    Double I_KINETIC_Fx2y_G2xyz_vrr = PBZ*I_KINETIC_Fx2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G2xyz_vrr;
    Double I_KINETIC_Fxyz_G2xyz_vrr = PBZ*I_KINETIC_Fxyz_F2xy_vrr+oned2z*I_KINETIC_Dxy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G2xyz_vrr;
    Double I_KINETIC_Fx2z_G2xyz_vrr = PBZ*I_KINETIC_Fx2z_F2xy_vrr+2*oned2z*I_KINETIC_Dxz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G2xyz_vrr;
    Double I_KINETIC_F3y_G2xyz_vrr = PBZ*I_KINETIC_F3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G2xyz_vrr;
    Double I_KINETIC_F2yz_G2xyz_vrr = PBZ*I_KINETIC_F2yz_F2xy_vrr+oned2z*I_KINETIC_D2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G2xyz_vrr;
    Double I_KINETIC_Fy2z_G2xyz_vrr = PBZ*I_KINETIC_Fy2z_F2xy_vrr+2*oned2z*I_KINETIC_Dyz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G2xyz_vrr;
    Double I_KINETIC_F3z_G2xyz_vrr = PBZ*I_KINETIC_F3z_F2xy_vrr+3*oned2z*I_KINETIC_D2z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G2xyz_vrr;
    Double I_KINETIC_F3x_G2x2z_vrr = PBZ*I_KINETIC_F3x_F2xz_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G2x2z_vrr-adz*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_KINETIC_F2xy_G2x2z_vrr = PBZ*I_KINETIC_F2xy_F2xz_vrr+oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G2x2z_vrr-adz*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_KINETIC_F2xz_G2x2z_vrr = PBZ*I_KINETIC_F2xz_F2xz_vrr+oned2z*I_KINETIC_D2x_F2xz_vrr+oned2z*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G2x2z_vrr-adz*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_KINETIC_Fx2y_G2x2z_vrr = PBZ*I_KINETIC_Fx2y_F2xz_vrr+oned2z*I_KINETIC_Fx2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_KINETIC_Fxyz_G2x2z_vrr = PBZ*I_KINETIC_Fxyz_F2xz_vrr+oned2z*I_KINETIC_Dxy_F2xz_vrr+oned2z*I_KINETIC_Fxyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_KINETIC_Fx2z_G2x2z_vrr = PBZ*I_KINETIC_Fx2z_F2xz_vrr+2*oned2z*I_KINETIC_Dxz_F2xz_vrr+oned2z*I_KINETIC_Fx2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_KINETIC_F3y_G2x2z_vrr = PBZ*I_KINETIC_F3y_F2xz_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G2x2z_vrr-adz*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_KINETIC_F2yz_G2x2z_vrr = PBZ*I_KINETIC_F2yz_F2xz_vrr+oned2z*I_KINETIC_D2y_F2xz_vrr+oned2z*I_KINETIC_F2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G2x2z_vrr-adz*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_KINETIC_Fy2z_G2x2z_vrr = PBZ*I_KINETIC_Fy2z_F2xz_vrr+2*oned2z*I_KINETIC_Dyz_F2xz_vrr+oned2z*I_KINETIC_Fy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_KINETIC_F3z_G2x2z_vrr = PBZ*I_KINETIC_F3z_F2xz_vrr+3*oned2z*I_KINETIC_D2z_F2xz_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_KINETIC_F3x_Gx3y_vrr = PBX*I_KINETIC_F3x_F3y_vrr+3*oned2z*I_KINETIC_D2x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Gx3y_vrr;
    Double I_KINETIC_F2xy_Gx3y_vrr = PBX*I_KINETIC_F2xy_F3y_vrr+2*oned2z*I_KINETIC_Dxy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Gx3y_vrr;
    Double I_KINETIC_F2xz_Gx3y_vrr = PBX*I_KINETIC_F2xz_F3y_vrr+2*oned2z*I_KINETIC_Dxz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Gx3y_vrr;
    Double I_KINETIC_Fx2y_Gx3y_vrr = PBX*I_KINETIC_Fx2y_F3y_vrr+oned2z*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Gx3y_vrr;
    Double I_KINETIC_Fxyz_Gx3y_vrr = PBX*I_KINETIC_Fxyz_F3y_vrr+oned2z*I_KINETIC_Dyz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Gx3y_vrr;
    Double I_KINETIC_Fx2z_Gx3y_vrr = PBX*I_KINETIC_Fx2z_F3y_vrr+oned2z*I_KINETIC_D2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Gx3y_vrr;
    Double I_KINETIC_F3y_Gx3y_vrr = PBX*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Gx3y_vrr;
    Double I_KINETIC_F2yz_Gx3y_vrr = PBX*I_KINETIC_F2yz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Gx3y_vrr;
    Double I_KINETIC_Fy2z_Gx3y_vrr = PBX*I_KINETIC_Fy2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Gx3y_vrr;
    Double I_KINETIC_F3z_Gx3y_vrr = PBX*I_KINETIC_F3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Gx3y_vrr;
    Double I_KINETIC_F3x_Gx2yz_vrr = PBZ*I_KINETIC_F3x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr;
    Double I_KINETIC_F2xy_Gx2yz_vrr = PBZ*I_KINETIC_F2xy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Gx2yz_vrr;
    Double I_KINETIC_F2xz_Gx2yz_vrr = PBZ*I_KINETIC_F2xz_Fx2y_vrr+oned2z*I_KINETIC_D2x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Gx2yz_vrr;
    Double I_KINETIC_Fx2y_Gx2yz_vrr = PBZ*I_KINETIC_Fx2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Gx2yz_vrr;
    Double I_KINETIC_Fxyz_Gx2yz_vrr = PBZ*I_KINETIC_Fxyz_Fx2y_vrr+oned2z*I_KINETIC_Dxy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Gx2yz_vrr;
    Double I_KINETIC_Fx2z_Gx2yz_vrr = PBZ*I_KINETIC_Fx2z_Fx2y_vrr+2*oned2z*I_KINETIC_Dxz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Gx2yz_vrr;
    Double I_KINETIC_F3y_Gx2yz_vrr = PBZ*I_KINETIC_F3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr;
    Double I_KINETIC_F2yz_Gx2yz_vrr = PBZ*I_KINETIC_F2yz_Fx2y_vrr+oned2z*I_KINETIC_D2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Gx2yz_vrr;
    Double I_KINETIC_Fy2z_Gx2yz_vrr = PBZ*I_KINETIC_Fy2z_Fx2y_vrr+2*oned2z*I_KINETIC_Dyz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Gx2yz_vrr;
    Double I_KINETIC_F3z_Gx2yz_vrr = PBZ*I_KINETIC_F3z_Fx2y_vrr+3*oned2z*I_KINETIC_D2z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr;
    Double I_KINETIC_F3x_Gxy2z_vrr = PBY*I_KINETIC_F3x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr;
    Double I_KINETIC_F2xy_Gxy2z_vrr = PBY*I_KINETIC_F2xy_Fx2z_vrr+oned2z*I_KINETIC_D2x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Gxy2z_vrr;
    Double I_KINETIC_F2xz_Gxy2z_vrr = PBY*I_KINETIC_F2xz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Gxy2z_vrr;
    Double I_KINETIC_Fx2y_Gxy2z_vrr = PBY*I_KINETIC_Fx2y_Fx2z_vrr+2*oned2z*I_KINETIC_Dxy_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Gxy2z_vrr;
    Double I_KINETIC_Fxyz_Gxy2z_vrr = PBY*I_KINETIC_Fxyz_Fx2z_vrr+oned2z*I_KINETIC_Dxz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Gxy2z_vrr;
    Double I_KINETIC_Fx2z_Gxy2z_vrr = PBY*I_KINETIC_Fx2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Gxy2z_vrr;
    Double I_KINETIC_F3y_Gxy2z_vrr = PBY*I_KINETIC_F3y_Fx2z_vrr+3*oned2z*I_KINETIC_D2y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr;
    Double I_KINETIC_F2yz_Gxy2z_vrr = PBY*I_KINETIC_F2yz_Fx2z_vrr+2*oned2z*I_KINETIC_Dyz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Gxy2z_vrr;
    Double I_KINETIC_Fy2z_Gxy2z_vrr = PBY*I_KINETIC_Fy2z_Fx2z_vrr+oned2z*I_KINETIC_D2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Gxy2z_vrr;
    Double I_KINETIC_F3z_Gxy2z_vrr = PBY*I_KINETIC_F3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr;
    Double I_KINETIC_F3x_Gx3z_vrr = PBX*I_KINETIC_F3x_F3z_vrr+3*oned2z*I_KINETIC_D2x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Gx3z_vrr;
    Double I_KINETIC_F2xy_Gx3z_vrr = PBX*I_KINETIC_F2xy_F3z_vrr+2*oned2z*I_KINETIC_Dxy_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Gx3z_vrr;
    Double I_KINETIC_F2xz_Gx3z_vrr = PBX*I_KINETIC_F2xz_F3z_vrr+2*oned2z*I_KINETIC_Dxz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Gx3z_vrr;
    Double I_KINETIC_Fx2y_Gx3z_vrr = PBX*I_KINETIC_Fx2y_F3z_vrr+oned2z*I_KINETIC_D2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Gx3z_vrr;
    Double I_KINETIC_Fxyz_Gx3z_vrr = PBX*I_KINETIC_Fxyz_F3z_vrr+oned2z*I_KINETIC_Dyz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Gx3z_vrr;
    Double I_KINETIC_Fx2z_Gx3z_vrr = PBX*I_KINETIC_Fx2z_F3z_vrr+oned2z*I_KINETIC_D2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Gx3z_vrr;
    Double I_KINETIC_F3y_Gx3z_vrr = PBX*I_KINETIC_F3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Gx3z_vrr;
    Double I_KINETIC_F2yz_Gx3z_vrr = PBX*I_KINETIC_F2yz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Gx3z_vrr;
    Double I_KINETIC_Fy2z_Gx3z_vrr = PBX*I_KINETIC_Fy2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Gx3z_vrr;
    Double I_KINETIC_F3z_Gx3z_vrr = PBX*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Gx3z_vrr;
    Double I_KINETIC_F3x_G4y_vrr = PBY*I_KINETIC_F3x_F3y_vrr+3*oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G4y_vrr-3*adz*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_KINETIC_F2xy_G4y_vrr = PBY*I_KINETIC_F2xy_F3y_vrr+oned2z*I_KINETIC_D2x_F3y_vrr+3*oned2z*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G4y_vrr-3*adz*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_KINETIC_F2xz_G4y_vrr = PBY*I_KINETIC_F2xz_F3y_vrr+3*oned2z*I_KINETIC_F2xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G4y_vrr-3*adz*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_KINETIC_Fx2y_G4y_vrr = PBY*I_KINETIC_Fx2y_F3y_vrr+2*oned2z*I_KINETIC_Dxy_F3y_vrr+3*oned2z*I_KINETIC_Fx2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_KINETIC_Fxyz_G4y_vrr = PBY*I_KINETIC_Fxyz_F3y_vrr+oned2z*I_KINETIC_Dxz_F3y_vrr+3*oned2z*I_KINETIC_Fxyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_KINETIC_Fx2z_G4y_vrr = PBY*I_KINETIC_Fx2z_F3y_vrr+3*oned2z*I_KINETIC_Fx2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_KINETIC_F3y_G4y_vrr = PBY*I_KINETIC_F3y_F3y_vrr+3*oned2z*I_KINETIC_D2y_F3y_vrr+3*oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G4y_vrr-3*adz*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_KINETIC_F2yz_G4y_vrr = PBY*I_KINETIC_F2yz_F3y_vrr+2*oned2z*I_KINETIC_Dyz_F3y_vrr+3*oned2z*I_KINETIC_F2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G4y_vrr-3*adz*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_KINETIC_Fy2z_G4y_vrr = PBY*I_KINETIC_Fy2z_F3y_vrr+oned2z*I_KINETIC_D2z_F3y_vrr+3*oned2z*I_KINETIC_Fy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_KINETIC_F3z_G4y_vrr = PBY*I_KINETIC_F3z_F3y_vrr+3*oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_KINETIC_F3x_G3yz_vrr = PBZ*I_KINETIC_F3x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G3yz_vrr;
    Double I_KINETIC_F2xy_G3yz_vrr = PBZ*I_KINETIC_F2xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G3yz_vrr;
    Double I_KINETIC_F2xz_G3yz_vrr = PBZ*I_KINETIC_F2xz_F3y_vrr+oned2z*I_KINETIC_D2x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G3yz_vrr;
    Double I_KINETIC_Fx2y_G3yz_vrr = PBZ*I_KINETIC_Fx2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G3yz_vrr;
    Double I_KINETIC_Fxyz_G3yz_vrr = PBZ*I_KINETIC_Fxyz_F3y_vrr+oned2z*I_KINETIC_Dxy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G3yz_vrr;
    Double I_KINETIC_Fx2z_G3yz_vrr = PBZ*I_KINETIC_Fx2z_F3y_vrr+2*oned2z*I_KINETIC_Dxz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G3yz_vrr;
    Double I_KINETIC_F3y_G3yz_vrr = PBZ*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G3yz_vrr;
    Double I_KINETIC_F2yz_G3yz_vrr = PBZ*I_KINETIC_F2yz_F3y_vrr+oned2z*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G3yz_vrr;
    Double I_KINETIC_Fy2z_G3yz_vrr = PBZ*I_KINETIC_Fy2z_F3y_vrr+2*oned2z*I_KINETIC_Dyz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G3yz_vrr;
    Double I_KINETIC_F3z_G3yz_vrr = PBZ*I_KINETIC_F3z_F3y_vrr+3*oned2z*I_KINETIC_D2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G3yz_vrr;
    Double I_KINETIC_F3x_G2y2z_vrr = PBZ*I_KINETIC_F3x_F2yz_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G2y2z_vrr-adz*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_KINETIC_F2xy_G2y2z_vrr = PBZ*I_KINETIC_F2xy_F2yz_vrr+oned2z*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G2y2z_vrr-adz*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_KINETIC_F2xz_G2y2z_vrr = PBZ*I_KINETIC_F2xz_F2yz_vrr+oned2z*I_KINETIC_D2x_F2yz_vrr+oned2z*I_KINETIC_F2xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G2y2z_vrr-adz*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_KINETIC_Fx2y_G2y2z_vrr = PBZ*I_KINETIC_Fx2y_F2yz_vrr+oned2z*I_KINETIC_Fx2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_KINETIC_Fxyz_G2y2z_vrr = PBZ*I_KINETIC_Fxyz_F2yz_vrr+oned2z*I_KINETIC_Dxy_F2yz_vrr+oned2z*I_KINETIC_Fxyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_KINETIC_Fx2z_G2y2z_vrr = PBZ*I_KINETIC_Fx2z_F2yz_vrr+2*oned2z*I_KINETIC_Dxz_F2yz_vrr+oned2z*I_KINETIC_Fx2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_KINETIC_F3y_G2y2z_vrr = PBZ*I_KINETIC_F3y_F2yz_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G2y2z_vrr-adz*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_KINETIC_F2yz_G2y2z_vrr = PBZ*I_KINETIC_F2yz_F2yz_vrr+oned2z*I_KINETIC_D2y_F2yz_vrr+oned2z*I_KINETIC_F2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G2y2z_vrr-adz*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_KINETIC_Fy2z_G2y2z_vrr = PBZ*I_KINETIC_Fy2z_F2yz_vrr+2*oned2z*I_KINETIC_Dyz_F2yz_vrr+oned2z*I_KINETIC_Fy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_KINETIC_F3z_G2y2z_vrr = PBZ*I_KINETIC_F3z_F2yz_vrr+3*oned2z*I_KINETIC_D2z_F2yz_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_KINETIC_F3x_Gy3z_vrr = PBY*I_KINETIC_F3x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Gy3z_vrr;
    Double I_KINETIC_F2xy_Gy3z_vrr = PBY*I_KINETIC_F2xy_F3z_vrr+oned2z*I_KINETIC_D2x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Gy3z_vrr;
    Double I_KINETIC_F2xz_Gy3z_vrr = PBY*I_KINETIC_F2xz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Gy3z_vrr;
    Double I_KINETIC_Fx2y_Gy3z_vrr = PBY*I_KINETIC_Fx2y_F3z_vrr+2*oned2z*I_KINETIC_Dxy_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Gy3z_vrr;
    Double I_KINETIC_Fxyz_Gy3z_vrr = PBY*I_KINETIC_Fxyz_F3z_vrr+oned2z*I_KINETIC_Dxz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Gy3z_vrr;
    Double I_KINETIC_Fx2z_Gy3z_vrr = PBY*I_KINETIC_Fx2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Gy3z_vrr;
    Double I_KINETIC_F3y_Gy3z_vrr = PBY*I_KINETIC_F3y_F3z_vrr+3*oned2z*I_KINETIC_D2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Gy3z_vrr;
    Double I_KINETIC_F2yz_Gy3z_vrr = PBY*I_KINETIC_F2yz_F3z_vrr+2*oned2z*I_KINETIC_Dyz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Gy3z_vrr;
    Double I_KINETIC_Fy2z_Gy3z_vrr = PBY*I_KINETIC_Fy2z_F3z_vrr+oned2z*I_KINETIC_D2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Gy3z_vrr;
    Double I_KINETIC_F3z_Gy3z_vrr = PBY*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Gy3z_vrr;
    Double I_KINETIC_F3x_G4z_vrr = PBZ*I_KINETIC_F3x_F3z_vrr+3*oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_G4z_vrr-3*adz*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_KINETIC_F2xy_G4z_vrr = PBZ*I_KINETIC_F2xy_F3z_vrr+3*oned2z*I_KINETIC_F2xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_G4z_vrr-3*adz*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_KINETIC_F2xz_G4z_vrr = PBZ*I_KINETIC_F2xz_F3z_vrr+oned2z*I_KINETIC_D2x_F3z_vrr+3*oned2z*I_KINETIC_F2xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_G4z_vrr-3*adz*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_KINETIC_Fx2y_G4z_vrr = PBZ*I_KINETIC_Fx2y_F3z_vrr+3*oned2z*I_KINETIC_Fx2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_KINETIC_Fxyz_G4z_vrr = PBZ*I_KINETIC_Fxyz_F3z_vrr+oned2z*I_KINETIC_Dxy_F3z_vrr+3*oned2z*I_KINETIC_Fxyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Fxyz_D2z_vrr;
    Double I_KINETIC_Fx2z_G4z_vrr = PBZ*I_KINETIC_Fx2z_F3z_vrr+2*oned2z*I_KINETIC_Dxz_F3z_vrr+3*oned2z*I_KINETIC_Fx2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_KINETIC_F3y_G4z_vrr = PBZ*I_KINETIC_F3y_F3z_vrr+3*oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_G4z_vrr-3*adz*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_KINETIC_F2yz_G4z_vrr = PBZ*I_KINETIC_F2yz_F3z_vrr+oned2z*I_KINETIC_D2y_F3z_vrr+3*oned2z*I_KINETIC_F2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_G4z_vrr-3*adz*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_KINETIC_Fy2z_G4z_vrr = PBZ*I_KINETIC_Fy2z_F3z_vrr+2*oned2z*I_KINETIC_Dyz_F3z_vrr+3*oned2z*I_KINETIC_Fy2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Fy2z_D2z_vrr;
    Double I_KINETIC_F3z_G4z_vrr = PBZ*I_KINETIC_F3z_F3z_vrr+3*oned2z*I_KINETIC_D2z_F3z_vrr+3*oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_F3z_D2z_vrr;

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
     * shell quartet name: SQ_KINETIC_G_G
     * expanding position: BRA1
     * code section is: VRR
     * totally 45 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_G
     * RHS shell quartet name: SQ_KINETIC_D_G
     * RHS shell quartet name: SQ_KINETIC_F_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_G
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_G
     ************************************************************/
    Double I_KINETIC_G4x_G4x_vrr = PAX*I_KINETIC_F3x_G4x_vrr+3*oned2z*I_KINETIC_D2x_G4x_vrr+4*oned2z*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G4x_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G4x_vrr;
    Double I_KINETIC_G3xy_G4x_vrr = PAY*I_KINETIC_F3x_G4x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G4x_vrr;
    Double I_KINETIC_G3xz_G4x_vrr = PAZ*I_KINETIC_F3x_G4x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G4x_vrr;
    Double I_KINETIC_G2x2y_G4x_vrr = PAY*I_KINETIC_F2xy_G4x_vrr+oned2z*I_KINETIC_D2x_G4x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G4x_vrr-bdz*I_TWOBODYOVERLAP_D2x_G4x_vrr;
    Double I_KINETIC_G2x2z_G4x_vrr = PAZ*I_KINETIC_F2xz_G4x_vrr+oned2z*I_KINETIC_D2x_G4x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G4x_vrr-bdz*I_TWOBODYOVERLAP_D2x_G4x_vrr;
    Double I_KINETIC_Gx3y_G4x_vrr = PAX*I_KINETIC_F3y_G4x_vrr+4*oned2z*I_KINETIC_F3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G4x_vrr;
    Double I_KINETIC_Gx3z_G4x_vrr = PAX*I_KINETIC_F3z_G4x_vrr+4*oned2z*I_KINETIC_F3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G4x_vrr;
    Double I_KINETIC_G4y_G4x_vrr = PAY*I_KINETIC_F3y_G4x_vrr+3*oned2z*I_KINETIC_D2y_G4x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G4x_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G4x_vrr;
    Double I_KINETIC_G3yz_G4x_vrr = PAZ*I_KINETIC_F3y_G4x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G4x_vrr;
    Double I_KINETIC_G2y2z_G4x_vrr = PAZ*I_KINETIC_F2yz_G4x_vrr+oned2z*I_KINETIC_D2y_G4x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G4x_vrr-bdz*I_TWOBODYOVERLAP_D2y_G4x_vrr;
    Double I_KINETIC_Gy3z_G4x_vrr = PAY*I_KINETIC_F3z_G4x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G4x_vrr;
    Double I_KINETIC_G4z_G4x_vrr = PAZ*I_KINETIC_F3z_G4x_vrr+3*oned2z*I_KINETIC_D2z_G4x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G4x_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G4x_vrr;
    Double I_KINETIC_G4x_G3xy_vrr = PAX*I_KINETIC_F3x_G3xy_vrr+3*oned2z*I_KINETIC_D2x_G3xy_vrr+3*oned2z*I_KINETIC_F3x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G3xy_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G3xy_vrr;
    Double I_KINETIC_G3xy_G3xy_vrr = PAY*I_KINETIC_F3x_G3xy_vrr+oned2z*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G3xy_vrr;
    Double I_KINETIC_G3xz_G3xy_vrr = PAZ*I_KINETIC_F3x_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G3xy_vrr;
    Double I_KINETIC_G2x2y_G3xy_vrr = PAY*I_KINETIC_F2xy_G3xy_vrr+oned2z*I_KINETIC_D2x_G3xy_vrr+oned2z*I_KINETIC_F2xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G3xy_vrr-bdz*I_TWOBODYOVERLAP_D2x_G3xy_vrr;
    Double I_KINETIC_G2x2z_G3xy_vrr = PAZ*I_KINETIC_F2xz_G3xy_vrr+oned2z*I_KINETIC_D2x_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G3xy_vrr-bdz*I_TWOBODYOVERLAP_D2x_G3xy_vrr;
    Double I_KINETIC_Gx3y_G3xy_vrr = PAX*I_KINETIC_F3y_G3xy_vrr+3*oned2z*I_KINETIC_F3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G3xy_vrr;
    Double I_KINETIC_Gx3z_G3xy_vrr = PAX*I_KINETIC_F3z_G3xy_vrr+3*oned2z*I_KINETIC_F3z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G3xy_vrr;
    Double I_KINETIC_G4y_G3xy_vrr = PAY*I_KINETIC_F3y_G3xy_vrr+3*oned2z*I_KINETIC_D2y_G3xy_vrr+oned2z*I_KINETIC_F3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G3xy_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G3xy_vrr;
    Double I_KINETIC_G3yz_G3xy_vrr = PAZ*I_KINETIC_F3y_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G3xy_vrr;
    Double I_KINETIC_G2y2z_G3xy_vrr = PAZ*I_KINETIC_F2yz_G3xy_vrr+oned2z*I_KINETIC_D2y_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G3xy_vrr-bdz*I_TWOBODYOVERLAP_D2y_G3xy_vrr;
    Double I_KINETIC_Gy3z_G3xy_vrr = PAY*I_KINETIC_F3z_G3xy_vrr+oned2z*I_KINETIC_F3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G3xy_vrr;
    Double I_KINETIC_G4z_G3xy_vrr = PAZ*I_KINETIC_F3z_G3xy_vrr+3*oned2z*I_KINETIC_D2z_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G3xy_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G3xy_vrr;
    Double I_KINETIC_G4x_G3xz_vrr = PAX*I_KINETIC_F3x_G3xz_vrr+3*oned2z*I_KINETIC_D2x_G3xz_vrr+3*oned2z*I_KINETIC_F3x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G3xz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G3xz_vrr;
    Double I_KINETIC_G3xy_G3xz_vrr = PAY*I_KINETIC_F3x_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G3xz_vrr;
    Double I_KINETIC_G3xz_G3xz_vrr = PAZ*I_KINETIC_F3x_G3xz_vrr+oned2z*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G3xz_vrr;
    Double I_KINETIC_G2x2y_G3xz_vrr = PAY*I_KINETIC_F2xy_G3xz_vrr+oned2z*I_KINETIC_D2x_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G3xz_vrr-bdz*I_TWOBODYOVERLAP_D2x_G3xz_vrr;
    Double I_KINETIC_G2x2z_G3xz_vrr = PAZ*I_KINETIC_F2xz_G3xz_vrr+oned2z*I_KINETIC_D2x_G3xz_vrr+oned2z*I_KINETIC_F2xz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G3xz_vrr-bdz*I_TWOBODYOVERLAP_D2x_G3xz_vrr;
    Double I_KINETIC_Gx3y_G3xz_vrr = PAX*I_KINETIC_F3y_G3xz_vrr+3*oned2z*I_KINETIC_F3y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G3xz_vrr;
    Double I_KINETIC_Gx3z_G3xz_vrr = PAX*I_KINETIC_F3z_G3xz_vrr+3*oned2z*I_KINETIC_F3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G3xz_vrr;
    Double I_KINETIC_G4y_G3xz_vrr = PAY*I_KINETIC_F3y_G3xz_vrr+3*oned2z*I_KINETIC_D2y_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G3xz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G3xz_vrr;
    Double I_KINETIC_G3yz_G3xz_vrr = PAZ*I_KINETIC_F3y_G3xz_vrr+oned2z*I_KINETIC_F3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G3xz_vrr;
    Double I_KINETIC_G2y2z_G3xz_vrr = PAZ*I_KINETIC_F2yz_G3xz_vrr+oned2z*I_KINETIC_D2y_G3xz_vrr+oned2z*I_KINETIC_F2yz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G3xz_vrr-bdz*I_TWOBODYOVERLAP_D2y_G3xz_vrr;
    Double I_KINETIC_Gy3z_G3xz_vrr = PAY*I_KINETIC_F3z_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G3xz_vrr;
    Double I_KINETIC_G4z_G3xz_vrr = PAZ*I_KINETIC_F3z_G3xz_vrr+3*oned2z*I_KINETIC_D2z_G3xz_vrr+oned2z*I_KINETIC_F3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G3xz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G3xz_vrr;
    Double I_KINETIC_G4x_G2x2y_vrr = PAX*I_KINETIC_F3x_G2x2y_vrr+3*oned2z*I_KINETIC_D2x_G2x2y_vrr+2*oned2z*I_KINETIC_F3x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G2x2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G2x2y_vrr;
    Double I_KINETIC_G3xy_G2x2y_vrr = PAY*I_KINETIC_F3x_G2x2y_vrr+2*oned2z*I_KINETIC_F3x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G2x2y_vrr;
    Double I_KINETIC_G3xz_G2x2y_vrr = PAZ*I_KINETIC_F3x_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G2x2y_vrr;
    Double I_KINETIC_G2x2y_G2x2y_vrr = PAY*I_KINETIC_F2xy_G2x2y_vrr+oned2z*I_KINETIC_D2x_G2x2y_vrr+2*oned2z*I_KINETIC_F2xy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G2x2y_vrr-bdz*I_TWOBODYOVERLAP_D2x_G2x2y_vrr;
    Double I_KINETIC_G2x2z_G2x2y_vrr = PAZ*I_KINETIC_F2xz_G2x2y_vrr+oned2z*I_KINETIC_D2x_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G2x2y_vrr-bdz*I_TWOBODYOVERLAP_D2x_G2x2y_vrr;
    Double I_KINETIC_Gx3y_G2x2y_vrr = PAX*I_KINETIC_F3y_G2x2y_vrr+2*oned2z*I_KINETIC_F3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G2x2y_vrr;
    Double I_KINETIC_Gx3z_G2x2y_vrr = PAX*I_KINETIC_F3z_G2x2y_vrr+2*oned2z*I_KINETIC_F3z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G2x2y_vrr;
    Double I_KINETIC_G4y_G2x2y_vrr = PAY*I_KINETIC_F3y_G2x2y_vrr+3*oned2z*I_KINETIC_D2y_G2x2y_vrr+2*oned2z*I_KINETIC_F3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G2x2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G2x2y_vrr;
    Double I_KINETIC_G3yz_G2x2y_vrr = PAZ*I_KINETIC_F3y_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G2x2y_vrr;
    Double I_KINETIC_G2y2z_G2x2y_vrr = PAZ*I_KINETIC_F2yz_G2x2y_vrr+oned2z*I_KINETIC_D2y_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G2x2y_vrr-bdz*I_TWOBODYOVERLAP_D2y_G2x2y_vrr;
    Double I_KINETIC_Gy3z_G2x2y_vrr = PAY*I_KINETIC_F3z_G2x2y_vrr+2*oned2z*I_KINETIC_F3z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G2x2y_vrr;
    Double I_KINETIC_G4z_G2x2y_vrr = PAZ*I_KINETIC_F3z_G2x2y_vrr+3*oned2z*I_KINETIC_D2z_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G2x2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G2x2y_vrr;
    Double I_KINETIC_G4x_G2xyz_vrr = PAX*I_KINETIC_F3x_G2xyz_vrr+3*oned2z*I_KINETIC_D2x_G2xyz_vrr+2*oned2z*I_KINETIC_F3x_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G2xyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G2xyz_vrr;
    Double I_KINETIC_G3xy_G2xyz_vrr = PAY*I_KINETIC_F3x_G2xyz_vrr+oned2z*I_KINETIC_F3x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G2xyz_vrr;
    Double I_KINETIC_G3xz_G2xyz_vrr = PAZ*I_KINETIC_F3x_G2xyz_vrr+oned2z*I_KINETIC_F3x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G2xyz_vrr;
    Double I_KINETIC_G2x2y_G2xyz_vrr = PAY*I_KINETIC_F2xy_G2xyz_vrr+oned2z*I_KINETIC_D2x_G2xyz_vrr+oned2z*I_KINETIC_F2xy_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G2xyz_vrr-bdz*I_TWOBODYOVERLAP_D2x_G2xyz_vrr;
    Double I_KINETIC_G2x2z_G2xyz_vrr = PAZ*I_KINETIC_F2xz_G2xyz_vrr+oned2z*I_KINETIC_D2x_G2xyz_vrr+oned2z*I_KINETIC_F2xz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G2xyz_vrr-bdz*I_TWOBODYOVERLAP_D2x_G2xyz_vrr;
    Double I_KINETIC_Gx3y_G2xyz_vrr = PAX*I_KINETIC_F3y_G2xyz_vrr+2*oned2z*I_KINETIC_F3y_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G2xyz_vrr;
    Double I_KINETIC_Gx3z_G2xyz_vrr = PAX*I_KINETIC_F3z_G2xyz_vrr+2*oned2z*I_KINETIC_F3z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G2xyz_vrr;
    Double I_KINETIC_G4y_G2xyz_vrr = PAY*I_KINETIC_F3y_G2xyz_vrr+3*oned2z*I_KINETIC_D2y_G2xyz_vrr+oned2z*I_KINETIC_F3y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G2xyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G2xyz_vrr;
    Double I_KINETIC_G3yz_G2xyz_vrr = PAZ*I_KINETIC_F3y_G2xyz_vrr+oned2z*I_KINETIC_F3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G2xyz_vrr;
    Double I_KINETIC_G2y2z_G2xyz_vrr = PAZ*I_KINETIC_F2yz_G2xyz_vrr+oned2z*I_KINETIC_D2y_G2xyz_vrr+oned2z*I_KINETIC_F2yz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G2xyz_vrr-bdz*I_TWOBODYOVERLAP_D2y_G2xyz_vrr;
    Double I_KINETIC_Gy3z_G2xyz_vrr = PAY*I_KINETIC_F3z_G2xyz_vrr+oned2z*I_KINETIC_F3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G2xyz_vrr;
    Double I_KINETIC_G4z_G2xyz_vrr = PAZ*I_KINETIC_F3z_G2xyz_vrr+3*oned2z*I_KINETIC_D2z_G2xyz_vrr+oned2z*I_KINETIC_F3z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G2xyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G2xyz_vrr;
    Double I_KINETIC_G4x_G2x2z_vrr = PAX*I_KINETIC_F3x_G2x2z_vrr+3*oned2z*I_KINETIC_D2x_G2x2z_vrr+2*oned2z*I_KINETIC_F3x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G2x2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G2x2z_vrr;
    Double I_KINETIC_G3xy_G2x2z_vrr = PAY*I_KINETIC_F3x_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G2x2z_vrr;
    Double I_KINETIC_G3xz_G2x2z_vrr = PAZ*I_KINETIC_F3x_G2x2z_vrr+2*oned2z*I_KINETIC_F3x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G2x2z_vrr;
    Double I_KINETIC_G2x2y_G2x2z_vrr = PAY*I_KINETIC_F2xy_G2x2z_vrr+oned2z*I_KINETIC_D2x_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G2x2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_G2x2z_vrr;
    Double I_KINETIC_G2x2z_G2x2z_vrr = PAZ*I_KINETIC_F2xz_G2x2z_vrr+oned2z*I_KINETIC_D2x_G2x2z_vrr+2*oned2z*I_KINETIC_F2xz_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G2x2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_G2x2z_vrr;
    Double I_KINETIC_Gx3y_G2x2z_vrr = PAX*I_KINETIC_F3y_G2x2z_vrr+2*oned2z*I_KINETIC_F3y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G2x2z_vrr;
    Double I_KINETIC_Gx3z_G2x2z_vrr = PAX*I_KINETIC_F3z_G2x2z_vrr+2*oned2z*I_KINETIC_F3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G2x2z_vrr;
    Double I_KINETIC_G4y_G2x2z_vrr = PAY*I_KINETIC_F3y_G2x2z_vrr+3*oned2z*I_KINETIC_D2y_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G2x2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G2x2z_vrr;
    Double I_KINETIC_G3yz_G2x2z_vrr = PAZ*I_KINETIC_F3y_G2x2z_vrr+2*oned2z*I_KINETIC_F3y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G2x2z_vrr;
    Double I_KINETIC_G2y2z_G2x2z_vrr = PAZ*I_KINETIC_F2yz_G2x2z_vrr+oned2z*I_KINETIC_D2y_G2x2z_vrr+2*oned2z*I_KINETIC_F2yz_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G2x2z_vrr-bdz*I_TWOBODYOVERLAP_D2y_G2x2z_vrr;
    Double I_KINETIC_Gy3z_G2x2z_vrr = PAY*I_KINETIC_F3z_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G2x2z_vrr;
    Double I_KINETIC_G4z_G2x2z_vrr = PAZ*I_KINETIC_F3z_G2x2z_vrr+3*oned2z*I_KINETIC_D2z_G2x2z_vrr+2*oned2z*I_KINETIC_F3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G2x2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G2x2z_vrr;
    Double I_KINETIC_G4x_Gx3y_vrr = PAX*I_KINETIC_F3x_Gx3y_vrr+3*oned2z*I_KINETIC_D2x_Gx3y_vrr+oned2z*I_KINETIC_F3x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Gx3y_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Gx3y_vrr;
    Double I_KINETIC_G3xy_Gx3y_vrr = PAY*I_KINETIC_F3x_Gx3y_vrr+3*oned2z*I_KINETIC_F3x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Gx3y_vrr;
    Double I_KINETIC_G3xz_Gx3y_vrr = PAZ*I_KINETIC_F3x_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Gx3y_vrr;
    Double I_KINETIC_G2x2y_Gx3y_vrr = PAY*I_KINETIC_F2xy_Gx3y_vrr+oned2z*I_KINETIC_D2x_Gx3y_vrr+3*oned2z*I_KINETIC_F2xy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Gx3y_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gx3y_vrr;
    Double I_KINETIC_G2x2z_Gx3y_vrr = PAZ*I_KINETIC_F2xz_Gx3y_vrr+oned2z*I_KINETIC_D2x_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Gx3y_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gx3y_vrr;
    Double I_KINETIC_Gx3y_Gx3y_vrr = PAX*I_KINETIC_F3y_Gx3y_vrr+oned2z*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Gx3y_vrr;
    Double I_KINETIC_Gx3z_Gx3y_vrr = PAX*I_KINETIC_F3z_Gx3y_vrr+oned2z*I_KINETIC_F3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Gx3y_vrr;
    Double I_KINETIC_G4y_Gx3y_vrr = PAY*I_KINETIC_F3y_Gx3y_vrr+3*oned2z*I_KINETIC_D2y_Gx3y_vrr+3*oned2z*I_KINETIC_F3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Gx3y_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Gx3y_vrr;
    Double I_KINETIC_G3yz_Gx3y_vrr = PAZ*I_KINETIC_F3y_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Gx3y_vrr;
    Double I_KINETIC_G2y2z_Gx3y_vrr = PAZ*I_KINETIC_F2yz_Gx3y_vrr+oned2z*I_KINETIC_D2y_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Gx3y_vrr-bdz*I_TWOBODYOVERLAP_D2y_Gx3y_vrr;
    Double I_KINETIC_Gy3z_Gx3y_vrr = PAY*I_KINETIC_F3z_Gx3y_vrr+3*oned2z*I_KINETIC_F3z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Gx3y_vrr;
    Double I_KINETIC_G4z_Gx3y_vrr = PAZ*I_KINETIC_F3z_Gx3y_vrr+3*oned2z*I_KINETIC_D2z_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Gx3y_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Gx3y_vrr;
    Double I_KINETIC_G4x_Gx2yz_vrr = PAX*I_KINETIC_F3x_Gx2yz_vrr+3*oned2z*I_KINETIC_D2x_Gx2yz_vrr+oned2z*I_KINETIC_F3x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Gx2yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Gx2yz_vrr;
    Double I_KINETIC_G3xy_Gx2yz_vrr = PAY*I_KINETIC_F3x_Gx2yz_vrr+2*oned2z*I_KINETIC_F3x_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Gx2yz_vrr;
    Double I_KINETIC_G3xz_Gx2yz_vrr = PAZ*I_KINETIC_F3x_Gx2yz_vrr+oned2z*I_KINETIC_F3x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Gx2yz_vrr;
    Double I_KINETIC_G2x2y_Gx2yz_vrr = PAY*I_KINETIC_F2xy_Gx2yz_vrr+oned2z*I_KINETIC_D2x_Gx2yz_vrr+2*oned2z*I_KINETIC_F2xy_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Gx2yz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gx2yz_vrr;
    Double I_KINETIC_G2x2z_Gx2yz_vrr = PAZ*I_KINETIC_F2xz_Gx2yz_vrr+oned2z*I_KINETIC_D2x_Gx2yz_vrr+oned2z*I_KINETIC_F2xz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Gx2yz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gx2yz_vrr;
    Double I_KINETIC_Gx3y_Gx2yz_vrr = PAX*I_KINETIC_F3y_Gx2yz_vrr+oned2z*I_KINETIC_F3y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Gx2yz_vrr;
    Double I_KINETIC_Gx3z_Gx2yz_vrr = PAX*I_KINETIC_F3z_Gx2yz_vrr+oned2z*I_KINETIC_F3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Gx2yz_vrr;
    Double I_KINETIC_G4y_Gx2yz_vrr = PAY*I_KINETIC_F3y_Gx2yz_vrr+3*oned2z*I_KINETIC_D2y_Gx2yz_vrr+2*oned2z*I_KINETIC_F3y_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Gx2yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Gx2yz_vrr;
    Double I_KINETIC_G3yz_Gx2yz_vrr = PAZ*I_KINETIC_F3y_Gx2yz_vrr+oned2z*I_KINETIC_F3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Gx2yz_vrr;
    Double I_KINETIC_G2y2z_Gx2yz_vrr = PAZ*I_KINETIC_F2yz_Gx2yz_vrr+oned2z*I_KINETIC_D2y_Gx2yz_vrr+oned2z*I_KINETIC_F2yz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Gx2yz_vrr-bdz*I_TWOBODYOVERLAP_D2y_Gx2yz_vrr;
    Double I_KINETIC_Gy3z_Gx2yz_vrr = PAY*I_KINETIC_F3z_Gx2yz_vrr+2*oned2z*I_KINETIC_F3z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Gx2yz_vrr;
    Double I_KINETIC_G4z_Gx2yz_vrr = PAZ*I_KINETIC_F3z_Gx2yz_vrr+3*oned2z*I_KINETIC_D2z_Gx2yz_vrr+oned2z*I_KINETIC_F3z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Gx2yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Gx2yz_vrr;
    Double I_KINETIC_G4x_Gxy2z_vrr = PAX*I_KINETIC_F3x_Gxy2z_vrr+3*oned2z*I_KINETIC_D2x_Gxy2z_vrr+oned2z*I_KINETIC_F3x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Gxy2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Gxy2z_vrr;
    Double I_KINETIC_G3xy_Gxy2z_vrr = PAY*I_KINETIC_F3x_Gxy2z_vrr+oned2z*I_KINETIC_F3x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Gxy2z_vrr;
    Double I_KINETIC_G3xz_Gxy2z_vrr = PAZ*I_KINETIC_F3x_Gxy2z_vrr+2*oned2z*I_KINETIC_F3x_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Gxy2z_vrr;
    Double I_KINETIC_G2x2y_Gxy2z_vrr = PAY*I_KINETIC_F2xy_Gxy2z_vrr+oned2z*I_KINETIC_D2x_Gxy2z_vrr+oned2z*I_KINETIC_F2xy_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Gxy2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gxy2z_vrr;
    Double I_KINETIC_G2x2z_Gxy2z_vrr = PAZ*I_KINETIC_F2xz_Gxy2z_vrr+oned2z*I_KINETIC_D2x_Gxy2z_vrr+2*oned2z*I_KINETIC_F2xz_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Gxy2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gxy2z_vrr;
    Double I_KINETIC_Gx3y_Gxy2z_vrr = PAX*I_KINETIC_F3y_Gxy2z_vrr+oned2z*I_KINETIC_F3y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Gxy2z_vrr;
    Double I_KINETIC_Gx3z_Gxy2z_vrr = PAX*I_KINETIC_F3z_Gxy2z_vrr+oned2z*I_KINETIC_F3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Gxy2z_vrr;
    Double I_KINETIC_G4y_Gxy2z_vrr = PAY*I_KINETIC_F3y_Gxy2z_vrr+3*oned2z*I_KINETIC_D2y_Gxy2z_vrr+oned2z*I_KINETIC_F3y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Gxy2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Gxy2z_vrr;
    Double I_KINETIC_G3yz_Gxy2z_vrr = PAZ*I_KINETIC_F3y_Gxy2z_vrr+2*oned2z*I_KINETIC_F3y_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Gxy2z_vrr;
    Double I_KINETIC_G2y2z_Gxy2z_vrr = PAZ*I_KINETIC_F2yz_Gxy2z_vrr+oned2z*I_KINETIC_D2y_Gxy2z_vrr+2*oned2z*I_KINETIC_F2yz_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Gxy2z_vrr-bdz*I_TWOBODYOVERLAP_D2y_Gxy2z_vrr;
    Double I_KINETIC_Gy3z_Gxy2z_vrr = PAY*I_KINETIC_F3z_Gxy2z_vrr+oned2z*I_KINETIC_F3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Gxy2z_vrr;
    Double I_KINETIC_G4z_Gxy2z_vrr = PAZ*I_KINETIC_F3z_Gxy2z_vrr+3*oned2z*I_KINETIC_D2z_Gxy2z_vrr+2*oned2z*I_KINETIC_F3z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Gxy2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Gxy2z_vrr;
    Double I_KINETIC_G4x_Gx3z_vrr = PAX*I_KINETIC_F3x_Gx3z_vrr+3*oned2z*I_KINETIC_D2x_Gx3z_vrr+oned2z*I_KINETIC_F3x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Gx3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Gx3z_vrr;
    Double I_KINETIC_G3xy_Gx3z_vrr = PAY*I_KINETIC_F3x_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Gx3z_vrr;
    Double I_KINETIC_G3xz_Gx3z_vrr = PAZ*I_KINETIC_F3x_Gx3z_vrr+3*oned2z*I_KINETIC_F3x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Gx3z_vrr;
    Double I_KINETIC_G2x2y_Gx3z_vrr = PAY*I_KINETIC_F2xy_Gx3z_vrr+oned2z*I_KINETIC_D2x_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Gx3z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gx3z_vrr;
    Double I_KINETIC_G2x2z_Gx3z_vrr = PAZ*I_KINETIC_F2xz_Gx3z_vrr+oned2z*I_KINETIC_D2x_Gx3z_vrr+3*oned2z*I_KINETIC_F2xz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Gx3z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gx3z_vrr;
    Double I_KINETIC_Gx3y_Gx3z_vrr = PAX*I_KINETIC_F3y_Gx3z_vrr+oned2z*I_KINETIC_F3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Gx3z_vrr;
    Double I_KINETIC_Gx3z_Gx3z_vrr = PAX*I_KINETIC_F3z_Gx3z_vrr+oned2z*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Gx3z_vrr;
    Double I_KINETIC_G4y_Gx3z_vrr = PAY*I_KINETIC_F3y_Gx3z_vrr+3*oned2z*I_KINETIC_D2y_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Gx3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Gx3z_vrr;
    Double I_KINETIC_G3yz_Gx3z_vrr = PAZ*I_KINETIC_F3y_Gx3z_vrr+3*oned2z*I_KINETIC_F3y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Gx3z_vrr;
    Double I_KINETIC_G2y2z_Gx3z_vrr = PAZ*I_KINETIC_F2yz_Gx3z_vrr+oned2z*I_KINETIC_D2y_Gx3z_vrr+3*oned2z*I_KINETIC_F2yz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Gx3z_vrr-bdz*I_TWOBODYOVERLAP_D2y_Gx3z_vrr;
    Double I_KINETIC_Gy3z_Gx3z_vrr = PAY*I_KINETIC_F3z_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Gx3z_vrr;
    Double I_KINETIC_G4z_Gx3z_vrr = PAZ*I_KINETIC_F3z_Gx3z_vrr+3*oned2z*I_KINETIC_D2z_Gx3z_vrr+3*oned2z*I_KINETIC_F3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Gx3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Gx3z_vrr;
    Double I_KINETIC_G4x_G4y_vrr = PAX*I_KINETIC_F3x_G4y_vrr+3*oned2z*I_KINETIC_D2x_G4y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G4y_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G4y_vrr;
    Double I_KINETIC_G3xy_G4y_vrr = PAY*I_KINETIC_F3x_G4y_vrr+4*oned2z*I_KINETIC_F3x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G4y_vrr;
    Double I_KINETIC_G3xz_G4y_vrr = PAZ*I_KINETIC_F3x_G4y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G4y_vrr;
    Double I_KINETIC_G2x2y_G4y_vrr = PAY*I_KINETIC_F2xy_G4y_vrr+oned2z*I_KINETIC_D2x_G4y_vrr+4*oned2z*I_KINETIC_F2xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G4y_vrr-bdz*I_TWOBODYOVERLAP_D2x_G4y_vrr;
    Double I_KINETIC_G2x2z_G4y_vrr = PAZ*I_KINETIC_F2xz_G4y_vrr+oned2z*I_KINETIC_D2x_G4y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G4y_vrr-bdz*I_TWOBODYOVERLAP_D2x_G4y_vrr;
    Double I_KINETIC_Gx3y_G4y_vrr = PAX*I_KINETIC_F3y_G4y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G4y_vrr;
    Double I_KINETIC_Gx3z_G4y_vrr = PAX*I_KINETIC_F3z_G4y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G4y_vrr;
    Double I_KINETIC_G4y_G4y_vrr = PAY*I_KINETIC_F3y_G4y_vrr+3*oned2z*I_KINETIC_D2y_G4y_vrr+4*oned2z*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G4y_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G4y_vrr;
    Double I_KINETIC_G3yz_G4y_vrr = PAZ*I_KINETIC_F3y_G4y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G4y_vrr;
    Double I_KINETIC_G2y2z_G4y_vrr = PAZ*I_KINETIC_F2yz_G4y_vrr+oned2z*I_KINETIC_D2y_G4y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G4y_vrr-bdz*I_TWOBODYOVERLAP_D2y_G4y_vrr;
    Double I_KINETIC_Gy3z_G4y_vrr = PAY*I_KINETIC_F3z_G4y_vrr+4*oned2z*I_KINETIC_F3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G4y_vrr;
    Double I_KINETIC_G4z_G4y_vrr = PAZ*I_KINETIC_F3z_G4y_vrr+3*oned2z*I_KINETIC_D2z_G4y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G4y_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G4y_vrr;
    Double I_KINETIC_G4x_G3yz_vrr = PAX*I_KINETIC_F3x_G3yz_vrr+3*oned2z*I_KINETIC_D2x_G3yz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G3yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G3yz_vrr;
    Double I_KINETIC_G3xy_G3yz_vrr = PAY*I_KINETIC_F3x_G3yz_vrr+3*oned2z*I_KINETIC_F3x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G3yz_vrr;
    Double I_KINETIC_G3xz_G3yz_vrr = PAZ*I_KINETIC_F3x_G3yz_vrr+oned2z*I_KINETIC_F3x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G3yz_vrr;
    Double I_KINETIC_G2x2y_G3yz_vrr = PAY*I_KINETIC_F2xy_G3yz_vrr+oned2z*I_KINETIC_D2x_G3yz_vrr+3*oned2z*I_KINETIC_F2xy_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G3yz_vrr-bdz*I_TWOBODYOVERLAP_D2x_G3yz_vrr;
    Double I_KINETIC_G2x2z_G3yz_vrr = PAZ*I_KINETIC_F2xz_G3yz_vrr+oned2z*I_KINETIC_D2x_G3yz_vrr+oned2z*I_KINETIC_F2xz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G3yz_vrr-bdz*I_TWOBODYOVERLAP_D2x_G3yz_vrr;
    Double I_KINETIC_Gx3y_G3yz_vrr = PAX*I_KINETIC_F3y_G3yz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G3yz_vrr;
    Double I_KINETIC_Gx3z_G3yz_vrr = PAX*I_KINETIC_F3z_G3yz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G3yz_vrr;
    Double I_KINETIC_G4y_G3yz_vrr = PAY*I_KINETIC_F3y_G3yz_vrr+3*oned2z*I_KINETIC_D2y_G3yz_vrr+3*oned2z*I_KINETIC_F3y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G3yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G3yz_vrr;
    Double I_KINETIC_G3yz_G3yz_vrr = PAZ*I_KINETIC_F3y_G3yz_vrr+oned2z*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G3yz_vrr;
    Double I_KINETIC_G2y2z_G3yz_vrr = PAZ*I_KINETIC_F2yz_G3yz_vrr+oned2z*I_KINETIC_D2y_G3yz_vrr+oned2z*I_KINETIC_F2yz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G3yz_vrr-bdz*I_TWOBODYOVERLAP_D2y_G3yz_vrr;
    Double I_KINETIC_Gy3z_G3yz_vrr = PAY*I_KINETIC_F3z_G3yz_vrr+3*oned2z*I_KINETIC_F3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G3yz_vrr;
    Double I_KINETIC_G4z_G3yz_vrr = PAZ*I_KINETIC_F3z_G3yz_vrr+3*oned2z*I_KINETIC_D2z_G3yz_vrr+oned2z*I_KINETIC_F3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G3yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G3yz_vrr;
    Double I_KINETIC_G4x_G2y2z_vrr = PAX*I_KINETIC_F3x_G2y2z_vrr+3*oned2z*I_KINETIC_D2x_G2y2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G2y2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G2y2z_vrr;
    Double I_KINETIC_G3xy_G2y2z_vrr = PAY*I_KINETIC_F3x_G2y2z_vrr+2*oned2z*I_KINETIC_F3x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G2y2z_vrr;
    Double I_KINETIC_G3xz_G2y2z_vrr = PAZ*I_KINETIC_F3x_G2y2z_vrr+2*oned2z*I_KINETIC_F3x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G2y2z_vrr;
    Double I_KINETIC_G2x2y_G2y2z_vrr = PAY*I_KINETIC_F2xy_G2y2z_vrr+oned2z*I_KINETIC_D2x_G2y2z_vrr+2*oned2z*I_KINETIC_F2xy_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G2y2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_G2y2z_vrr;
    Double I_KINETIC_G2x2z_G2y2z_vrr = PAZ*I_KINETIC_F2xz_G2y2z_vrr+oned2z*I_KINETIC_D2x_G2y2z_vrr+2*oned2z*I_KINETIC_F2xz_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G2y2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_G2y2z_vrr;
    Double I_KINETIC_Gx3y_G2y2z_vrr = PAX*I_KINETIC_F3y_G2y2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G2y2z_vrr;
    Double I_KINETIC_Gx3z_G2y2z_vrr = PAX*I_KINETIC_F3z_G2y2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G2y2z_vrr;
    Double I_KINETIC_G4y_G2y2z_vrr = PAY*I_KINETIC_F3y_G2y2z_vrr+3*oned2z*I_KINETIC_D2y_G2y2z_vrr+2*oned2z*I_KINETIC_F3y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G2y2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G2y2z_vrr;
    Double I_KINETIC_G3yz_G2y2z_vrr = PAZ*I_KINETIC_F3y_G2y2z_vrr+2*oned2z*I_KINETIC_F3y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G2y2z_vrr;
    Double I_KINETIC_G2y2z_G2y2z_vrr = PAZ*I_KINETIC_F2yz_G2y2z_vrr+oned2z*I_KINETIC_D2y_G2y2z_vrr+2*oned2z*I_KINETIC_F2yz_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G2y2z_vrr-bdz*I_TWOBODYOVERLAP_D2y_G2y2z_vrr;
    Double I_KINETIC_Gy3z_G2y2z_vrr = PAY*I_KINETIC_F3z_G2y2z_vrr+2*oned2z*I_KINETIC_F3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G2y2z_vrr;
    Double I_KINETIC_G4z_G2y2z_vrr = PAZ*I_KINETIC_F3z_G2y2z_vrr+3*oned2z*I_KINETIC_D2z_G2y2z_vrr+2*oned2z*I_KINETIC_F3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G2y2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G2y2z_vrr;
    Double I_KINETIC_G4x_Gy3z_vrr = PAX*I_KINETIC_F3x_Gy3z_vrr+3*oned2z*I_KINETIC_D2x_Gy3z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Gy3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Gy3z_vrr;
    Double I_KINETIC_G3xy_Gy3z_vrr = PAY*I_KINETIC_F3x_Gy3z_vrr+oned2z*I_KINETIC_F3x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Gy3z_vrr;
    Double I_KINETIC_G3xz_Gy3z_vrr = PAZ*I_KINETIC_F3x_Gy3z_vrr+3*oned2z*I_KINETIC_F3x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Gy3z_vrr;
    Double I_KINETIC_G2x2y_Gy3z_vrr = PAY*I_KINETIC_F2xy_Gy3z_vrr+oned2z*I_KINETIC_D2x_Gy3z_vrr+oned2z*I_KINETIC_F2xy_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Gy3z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gy3z_vrr;
    Double I_KINETIC_G2x2z_Gy3z_vrr = PAZ*I_KINETIC_F2xz_Gy3z_vrr+oned2z*I_KINETIC_D2x_Gy3z_vrr+3*oned2z*I_KINETIC_F2xz_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Gy3z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Gy3z_vrr;
    Double I_KINETIC_Gx3y_Gy3z_vrr = PAX*I_KINETIC_F3y_Gy3z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Gy3z_vrr;
    Double I_KINETIC_Gx3z_Gy3z_vrr = PAX*I_KINETIC_F3z_Gy3z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Gy3z_vrr;
    Double I_KINETIC_G4y_Gy3z_vrr = PAY*I_KINETIC_F3y_Gy3z_vrr+3*oned2z*I_KINETIC_D2y_Gy3z_vrr+oned2z*I_KINETIC_F3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Gy3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Gy3z_vrr;
    Double I_KINETIC_G3yz_Gy3z_vrr = PAZ*I_KINETIC_F3y_Gy3z_vrr+3*oned2z*I_KINETIC_F3y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Gy3z_vrr;
    Double I_KINETIC_G2y2z_Gy3z_vrr = PAZ*I_KINETIC_F2yz_Gy3z_vrr+oned2z*I_KINETIC_D2y_Gy3z_vrr+3*oned2z*I_KINETIC_F2yz_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Gy3z_vrr-bdz*I_TWOBODYOVERLAP_D2y_Gy3z_vrr;
    Double I_KINETIC_Gy3z_Gy3z_vrr = PAY*I_KINETIC_F3z_Gy3z_vrr+oned2z*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Gy3z_vrr;
    Double I_KINETIC_G4z_Gy3z_vrr = PAZ*I_KINETIC_F3z_Gy3z_vrr+3*oned2z*I_KINETIC_D2z_Gy3z_vrr+3*oned2z*I_KINETIC_F3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Gy3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Gy3z_vrr;
    Double I_KINETIC_G4x_G4z_vrr = PAX*I_KINETIC_F3x_G4z_vrr+3*oned2z*I_KINETIC_D2x_G4z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_G4z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_G4z_vrr;
    Double I_KINETIC_G3xy_G4z_vrr = PAY*I_KINETIC_F3x_G4z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_G4z_vrr;
    Double I_KINETIC_G3xz_G4z_vrr = PAZ*I_KINETIC_F3x_G4z_vrr+4*oned2z*I_KINETIC_F3x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_G4z_vrr;
    Double I_KINETIC_G2x2y_G4z_vrr = PAY*I_KINETIC_F2xy_G4z_vrr+oned2z*I_KINETIC_D2x_G4z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_G4z_vrr-bdz*I_TWOBODYOVERLAP_D2x_G4z_vrr;
    Double I_KINETIC_G2x2z_G4z_vrr = PAZ*I_KINETIC_F2xz_G4z_vrr+oned2z*I_KINETIC_D2x_G4z_vrr+4*oned2z*I_KINETIC_F2xz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_G4z_vrr-bdz*I_TWOBODYOVERLAP_D2x_G4z_vrr;
    Double I_KINETIC_Gx3y_G4z_vrr = PAX*I_KINETIC_F3y_G4z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_G4z_vrr;
    Double I_KINETIC_Gx3z_G4z_vrr = PAX*I_KINETIC_F3z_G4z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_G4z_vrr;
    Double I_KINETIC_G4y_G4z_vrr = PAY*I_KINETIC_F3y_G4z_vrr+3*oned2z*I_KINETIC_D2y_G4z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_G4z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_G4z_vrr;
    Double I_KINETIC_G3yz_G4z_vrr = PAZ*I_KINETIC_F3y_G4z_vrr+4*oned2z*I_KINETIC_F3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_G4z_vrr;
    Double I_KINETIC_G2y2z_G4z_vrr = PAZ*I_KINETIC_F2yz_G4z_vrr+oned2z*I_KINETIC_D2y_G4z_vrr+4*oned2z*I_KINETIC_F2yz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_G4z_vrr-bdz*I_TWOBODYOVERLAP_D2y_G4z_vrr;
    Double I_KINETIC_Gy3z_G4z_vrr = PAY*I_KINETIC_F3z_G4z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_G4z_vrr;
    Double I_KINETIC_G4z_G4z_vrr = PAZ*I_KINETIC_F3z_G4z_vrr+3*oned2z*I_KINETIC_D2z_G4z_vrr+4*oned2z*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G4z_G4z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_G4z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_G
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_G
     * RHS shell quartet name: SQ_KINETIC_F_G
     * RHS shell quartet name: SQ_KINETIC_G_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_G
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_G
     ************************************************************/
    Double I_KINETIC_H5x_G4x_vrr = PAX*I_KINETIC_G4x_G4x_vrr+4*oned2z*I_KINETIC_F3x_G4x_vrr+4*oned2z*I_KINETIC_G4x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G4x_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G4x_vrr;
    Double I_KINETIC_H4xy_G4x_vrr = PAY*I_KINETIC_G4x_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G4x_vrr;
    Double I_KINETIC_H4xz_G4x_vrr = PAZ*I_KINETIC_G4x_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G4x_vrr;
    Double I_KINETIC_H3x2y_G4x_vrr = PAY*I_KINETIC_G3xy_G4x_vrr+oned2z*I_KINETIC_F3x_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G4x_vrr-bdz*I_TWOBODYOVERLAP_F3x_G4x_vrr;
    Double I_KINETIC_H3xyz_G4x_vrr = PAZ*I_KINETIC_G3xy_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G4x_vrr;
    Double I_KINETIC_H3x2z_G4x_vrr = PAZ*I_KINETIC_G3xz_G4x_vrr+oned2z*I_KINETIC_F3x_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G4x_vrr-bdz*I_TWOBODYOVERLAP_F3x_G4x_vrr;
    Double I_KINETIC_H2x3y_G4x_vrr = PAX*I_KINETIC_Gx3y_G4x_vrr+oned2z*I_KINETIC_F3y_G4x_vrr+4*oned2z*I_KINETIC_Gx3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G4x_vrr-bdz*I_TWOBODYOVERLAP_F3y_G4x_vrr;
    Double I_KINETIC_H2x2yz_G4x_vrr = PAZ*I_KINETIC_G2x2y_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G4x_vrr;
    Double I_KINETIC_H2xy2z_G4x_vrr = PAY*I_KINETIC_G2x2z_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G4x_vrr;
    Double I_KINETIC_H2x3z_G4x_vrr = PAX*I_KINETIC_Gx3z_G4x_vrr+oned2z*I_KINETIC_F3z_G4x_vrr+4*oned2z*I_KINETIC_Gx3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G4x_vrr-bdz*I_TWOBODYOVERLAP_F3z_G4x_vrr;
    Double I_KINETIC_Hx4y_G4x_vrr = PAX*I_KINETIC_G4y_G4x_vrr+4*oned2z*I_KINETIC_G4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G4x_vrr;
    Double I_KINETIC_Hx3yz_G4x_vrr = PAZ*I_KINETIC_Gx3y_G4x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G4x_vrr;
    Double I_KINETIC_Hx2y2z_G4x_vrr = PAX*I_KINETIC_G2y2z_G4x_vrr+4*oned2z*I_KINETIC_G2y2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G4x_vrr;
    Double I_KINETIC_Hxy3z_G4x_vrr = PAY*I_KINETIC_Gx3z_G4x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G4x_vrr;
    Double I_KINETIC_Hx4z_G4x_vrr = PAX*I_KINETIC_G4z_G4x_vrr+4*oned2z*I_KINETIC_G4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G4x_vrr;
    Double I_KINETIC_H5y_G4x_vrr = PAY*I_KINETIC_G4y_G4x_vrr+4*oned2z*I_KINETIC_F3y_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G4x_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G4x_vrr;
    Double I_KINETIC_H4yz_G4x_vrr = PAZ*I_KINETIC_G4y_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G4x_vrr;
    Double I_KINETIC_H3y2z_G4x_vrr = PAZ*I_KINETIC_G3yz_G4x_vrr+oned2z*I_KINETIC_F3y_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G4x_vrr-bdz*I_TWOBODYOVERLAP_F3y_G4x_vrr;
    Double I_KINETIC_H2y3z_G4x_vrr = PAY*I_KINETIC_Gy3z_G4x_vrr+oned2z*I_KINETIC_F3z_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G4x_vrr-bdz*I_TWOBODYOVERLAP_F3z_G4x_vrr;
    Double I_KINETIC_Hy4z_G4x_vrr = PAY*I_KINETIC_G4z_G4x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G4x_vrr;
    Double I_KINETIC_H5z_G4x_vrr = PAZ*I_KINETIC_G4z_G4x_vrr+4*oned2z*I_KINETIC_F3z_G4x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G4x_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G4x_vrr;
    Double I_KINETIC_H5x_G3xy_vrr = PAX*I_KINETIC_G4x_G3xy_vrr+4*oned2z*I_KINETIC_F3x_G3xy_vrr+3*oned2z*I_KINETIC_G4x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G3xy_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G3xy_vrr;
    Double I_KINETIC_H4xy_G3xy_vrr = PAY*I_KINETIC_G4x_G3xy_vrr+oned2z*I_KINETIC_G4x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G3xy_vrr;
    Double I_KINETIC_H4xz_G3xy_vrr = PAZ*I_KINETIC_G4x_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G3xy_vrr;
    Double I_KINETIC_H3x2y_G3xy_vrr = PAY*I_KINETIC_G3xy_G3xy_vrr+oned2z*I_KINETIC_F3x_G3xy_vrr+oned2z*I_KINETIC_G3xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G3xy_vrr-bdz*I_TWOBODYOVERLAP_F3x_G3xy_vrr;
    Double I_KINETIC_H3xyz_G3xy_vrr = PAZ*I_KINETIC_G3xy_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G3xy_vrr;
    Double I_KINETIC_H3x2z_G3xy_vrr = PAZ*I_KINETIC_G3xz_G3xy_vrr+oned2z*I_KINETIC_F3x_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G3xy_vrr-bdz*I_TWOBODYOVERLAP_F3x_G3xy_vrr;
    Double I_KINETIC_H2x3y_G3xy_vrr = PAX*I_KINETIC_Gx3y_G3xy_vrr+oned2z*I_KINETIC_F3y_G3xy_vrr+3*oned2z*I_KINETIC_Gx3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G3xy_vrr-bdz*I_TWOBODYOVERLAP_F3y_G3xy_vrr;
    Double I_KINETIC_H2x2yz_G3xy_vrr = PAZ*I_KINETIC_G2x2y_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G3xy_vrr;
    Double I_KINETIC_H2xy2z_G3xy_vrr = PAY*I_KINETIC_G2x2z_G3xy_vrr+oned2z*I_KINETIC_G2x2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G3xy_vrr;
    Double I_KINETIC_H2x3z_G3xy_vrr = PAX*I_KINETIC_Gx3z_G3xy_vrr+oned2z*I_KINETIC_F3z_G3xy_vrr+3*oned2z*I_KINETIC_Gx3z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G3xy_vrr-bdz*I_TWOBODYOVERLAP_F3z_G3xy_vrr;
    Double I_KINETIC_Hx4y_G3xy_vrr = PAX*I_KINETIC_G4y_G3xy_vrr+3*oned2z*I_KINETIC_G4y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G3xy_vrr;
    Double I_KINETIC_Hx3yz_G3xy_vrr = PAZ*I_KINETIC_Gx3y_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G3xy_vrr;
    Double I_KINETIC_Hx2y2z_G3xy_vrr = PAX*I_KINETIC_G2y2z_G3xy_vrr+3*oned2z*I_KINETIC_G2y2z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G3xy_vrr;
    Double I_KINETIC_Hxy3z_G3xy_vrr = PAY*I_KINETIC_Gx3z_G3xy_vrr+oned2z*I_KINETIC_Gx3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G3xy_vrr;
    Double I_KINETIC_Hx4z_G3xy_vrr = PAX*I_KINETIC_G4z_G3xy_vrr+3*oned2z*I_KINETIC_G4z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G3xy_vrr;
    Double I_KINETIC_H5y_G3xy_vrr = PAY*I_KINETIC_G4y_G3xy_vrr+4*oned2z*I_KINETIC_F3y_G3xy_vrr+oned2z*I_KINETIC_G4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G3xy_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G3xy_vrr;
    Double I_KINETIC_H4yz_G3xy_vrr = PAZ*I_KINETIC_G4y_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G3xy_vrr;
    Double I_KINETIC_H3y2z_G3xy_vrr = PAZ*I_KINETIC_G3yz_G3xy_vrr+oned2z*I_KINETIC_F3y_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G3xy_vrr-bdz*I_TWOBODYOVERLAP_F3y_G3xy_vrr;
    Double I_KINETIC_H2y3z_G3xy_vrr = PAY*I_KINETIC_Gy3z_G3xy_vrr+oned2z*I_KINETIC_F3z_G3xy_vrr+oned2z*I_KINETIC_Gy3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G3xy_vrr-bdz*I_TWOBODYOVERLAP_F3z_G3xy_vrr;
    Double I_KINETIC_Hy4z_G3xy_vrr = PAY*I_KINETIC_G4z_G3xy_vrr+oned2z*I_KINETIC_G4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G3xy_vrr;
    Double I_KINETIC_H5z_G3xy_vrr = PAZ*I_KINETIC_G4z_G3xy_vrr+4*oned2z*I_KINETIC_F3z_G3xy_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G3xy_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G3xy_vrr;
    Double I_KINETIC_H5x_G3xz_vrr = PAX*I_KINETIC_G4x_G3xz_vrr+4*oned2z*I_KINETIC_F3x_G3xz_vrr+3*oned2z*I_KINETIC_G4x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G3xz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G3xz_vrr;
    Double I_KINETIC_H4xy_G3xz_vrr = PAY*I_KINETIC_G4x_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G3xz_vrr;
    Double I_KINETIC_H4xz_G3xz_vrr = PAZ*I_KINETIC_G4x_G3xz_vrr+oned2z*I_KINETIC_G4x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G3xz_vrr;
    Double I_KINETIC_H3x2y_G3xz_vrr = PAY*I_KINETIC_G3xy_G3xz_vrr+oned2z*I_KINETIC_F3x_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G3xz_vrr-bdz*I_TWOBODYOVERLAP_F3x_G3xz_vrr;
    Double I_KINETIC_H3xyz_G3xz_vrr = PAZ*I_KINETIC_G3xy_G3xz_vrr+oned2z*I_KINETIC_G3xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G3xz_vrr;
    Double I_KINETIC_H3x2z_G3xz_vrr = PAZ*I_KINETIC_G3xz_G3xz_vrr+oned2z*I_KINETIC_F3x_G3xz_vrr+oned2z*I_KINETIC_G3xz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G3xz_vrr-bdz*I_TWOBODYOVERLAP_F3x_G3xz_vrr;
    Double I_KINETIC_H2x3y_G3xz_vrr = PAX*I_KINETIC_Gx3y_G3xz_vrr+oned2z*I_KINETIC_F3y_G3xz_vrr+3*oned2z*I_KINETIC_Gx3y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G3xz_vrr-bdz*I_TWOBODYOVERLAP_F3y_G3xz_vrr;
    Double I_KINETIC_H2x2yz_G3xz_vrr = PAZ*I_KINETIC_G2x2y_G3xz_vrr+oned2z*I_KINETIC_G2x2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G3xz_vrr;
    Double I_KINETIC_H2xy2z_G3xz_vrr = PAY*I_KINETIC_G2x2z_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G3xz_vrr;
    Double I_KINETIC_H2x3z_G3xz_vrr = PAX*I_KINETIC_Gx3z_G3xz_vrr+oned2z*I_KINETIC_F3z_G3xz_vrr+3*oned2z*I_KINETIC_Gx3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G3xz_vrr-bdz*I_TWOBODYOVERLAP_F3z_G3xz_vrr;
    Double I_KINETIC_Hx4y_G3xz_vrr = PAX*I_KINETIC_G4y_G3xz_vrr+3*oned2z*I_KINETIC_G4y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G3xz_vrr;
    Double I_KINETIC_Hx3yz_G3xz_vrr = PAZ*I_KINETIC_Gx3y_G3xz_vrr+oned2z*I_KINETIC_Gx3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G3xz_vrr;
    Double I_KINETIC_Hx2y2z_G3xz_vrr = PAX*I_KINETIC_G2y2z_G3xz_vrr+3*oned2z*I_KINETIC_G2y2z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G3xz_vrr;
    Double I_KINETIC_Hxy3z_G3xz_vrr = PAY*I_KINETIC_Gx3z_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G3xz_vrr;
    Double I_KINETIC_Hx4z_G3xz_vrr = PAX*I_KINETIC_G4z_G3xz_vrr+3*oned2z*I_KINETIC_G4z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G3xz_vrr;
    Double I_KINETIC_H5y_G3xz_vrr = PAY*I_KINETIC_G4y_G3xz_vrr+4*oned2z*I_KINETIC_F3y_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G3xz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G3xz_vrr;
    Double I_KINETIC_H4yz_G3xz_vrr = PAZ*I_KINETIC_G4y_G3xz_vrr+oned2z*I_KINETIC_G4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G3xz_vrr;
    Double I_KINETIC_H3y2z_G3xz_vrr = PAZ*I_KINETIC_G3yz_G3xz_vrr+oned2z*I_KINETIC_F3y_G3xz_vrr+oned2z*I_KINETIC_G3yz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G3xz_vrr-bdz*I_TWOBODYOVERLAP_F3y_G3xz_vrr;
    Double I_KINETIC_H2y3z_G3xz_vrr = PAY*I_KINETIC_Gy3z_G3xz_vrr+oned2z*I_KINETIC_F3z_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G3xz_vrr-bdz*I_TWOBODYOVERLAP_F3z_G3xz_vrr;
    Double I_KINETIC_Hy4z_G3xz_vrr = PAY*I_KINETIC_G4z_G3xz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G3xz_vrr;
    Double I_KINETIC_H5z_G3xz_vrr = PAZ*I_KINETIC_G4z_G3xz_vrr+4*oned2z*I_KINETIC_F3z_G3xz_vrr+oned2z*I_KINETIC_G4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G3xz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G3xz_vrr;
    Double I_KINETIC_H5x_G2x2y_vrr = PAX*I_KINETIC_G4x_G2x2y_vrr+4*oned2z*I_KINETIC_F3x_G2x2y_vrr+2*oned2z*I_KINETIC_G4x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G2x2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G2x2y_vrr;
    Double I_KINETIC_H4xy_G2x2y_vrr = PAY*I_KINETIC_G4x_G2x2y_vrr+2*oned2z*I_KINETIC_G4x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G2x2y_vrr;
    Double I_KINETIC_H4xz_G2x2y_vrr = PAZ*I_KINETIC_G4x_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G2x2y_vrr;
    Double I_KINETIC_H3x2y_G2x2y_vrr = PAY*I_KINETIC_G3xy_G2x2y_vrr+oned2z*I_KINETIC_F3x_G2x2y_vrr+2*oned2z*I_KINETIC_G3xy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G2x2y_vrr-bdz*I_TWOBODYOVERLAP_F3x_G2x2y_vrr;
    Double I_KINETIC_H3xyz_G2x2y_vrr = PAZ*I_KINETIC_G3xy_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G2x2y_vrr;
    Double I_KINETIC_H3x2z_G2x2y_vrr = PAZ*I_KINETIC_G3xz_G2x2y_vrr+oned2z*I_KINETIC_F3x_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G2x2y_vrr-bdz*I_TWOBODYOVERLAP_F3x_G2x2y_vrr;
    Double I_KINETIC_H2x3y_G2x2y_vrr = PAX*I_KINETIC_Gx3y_G2x2y_vrr+oned2z*I_KINETIC_F3y_G2x2y_vrr+2*oned2z*I_KINETIC_Gx3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G2x2y_vrr-bdz*I_TWOBODYOVERLAP_F3y_G2x2y_vrr;
    Double I_KINETIC_H2x2yz_G2x2y_vrr = PAZ*I_KINETIC_G2x2y_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G2x2y_vrr;
    Double I_KINETIC_H2xy2z_G2x2y_vrr = PAY*I_KINETIC_G2x2z_G2x2y_vrr+2*oned2z*I_KINETIC_G2x2z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G2x2y_vrr;
    Double I_KINETIC_H2x3z_G2x2y_vrr = PAX*I_KINETIC_Gx3z_G2x2y_vrr+oned2z*I_KINETIC_F3z_G2x2y_vrr+2*oned2z*I_KINETIC_Gx3z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G2x2y_vrr-bdz*I_TWOBODYOVERLAP_F3z_G2x2y_vrr;
    Double I_KINETIC_Hx4y_G2x2y_vrr = PAX*I_KINETIC_G4y_G2x2y_vrr+2*oned2z*I_KINETIC_G4y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G2x2y_vrr;
    Double I_KINETIC_Hx3yz_G2x2y_vrr = PAZ*I_KINETIC_Gx3y_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G2x2y_vrr;
    Double I_KINETIC_Hx2y2z_G2x2y_vrr = PAX*I_KINETIC_G2y2z_G2x2y_vrr+2*oned2z*I_KINETIC_G2y2z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G2x2y_vrr;
    Double I_KINETIC_Hxy3z_G2x2y_vrr = PAY*I_KINETIC_Gx3z_G2x2y_vrr+2*oned2z*I_KINETIC_Gx3z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G2x2y_vrr;
    Double I_KINETIC_Hx4z_G2x2y_vrr = PAX*I_KINETIC_G4z_G2x2y_vrr+2*oned2z*I_KINETIC_G4z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G2x2y_vrr;
    Double I_KINETIC_H5y_G2x2y_vrr = PAY*I_KINETIC_G4y_G2x2y_vrr+4*oned2z*I_KINETIC_F3y_G2x2y_vrr+2*oned2z*I_KINETIC_G4y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G2x2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G2x2y_vrr;
    Double I_KINETIC_H4yz_G2x2y_vrr = PAZ*I_KINETIC_G4y_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G2x2y_vrr;
    Double I_KINETIC_H3y2z_G2x2y_vrr = PAZ*I_KINETIC_G3yz_G2x2y_vrr+oned2z*I_KINETIC_F3y_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G2x2y_vrr-bdz*I_TWOBODYOVERLAP_F3y_G2x2y_vrr;
    Double I_KINETIC_H2y3z_G2x2y_vrr = PAY*I_KINETIC_Gy3z_G2x2y_vrr+oned2z*I_KINETIC_F3z_G2x2y_vrr+2*oned2z*I_KINETIC_Gy3z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G2x2y_vrr-bdz*I_TWOBODYOVERLAP_F3z_G2x2y_vrr;
    Double I_KINETIC_Hy4z_G2x2y_vrr = PAY*I_KINETIC_G4z_G2x2y_vrr+2*oned2z*I_KINETIC_G4z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G2x2y_vrr;
    Double I_KINETIC_H5z_G2x2y_vrr = PAZ*I_KINETIC_G4z_G2x2y_vrr+4*oned2z*I_KINETIC_F3z_G2x2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G2x2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G2x2y_vrr;
    Double I_KINETIC_H5x_G2xyz_vrr = PAX*I_KINETIC_G4x_G2xyz_vrr+4*oned2z*I_KINETIC_F3x_G2xyz_vrr+2*oned2z*I_KINETIC_G4x_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G2xyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G2xyz_vrr;
    Double I_KINETIC_H4xy_G2xyz_vrr = PAY*I_KINETIC_G4x_G2xyz_vrr+oned2z*I_KINETIC_G4x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G2xyz_vrr;
    Double I_KINETIC_H4xz_G2xyz_vrr = PAZ*I_KINETIC_G4x_G2xyz_vrr+oned2z*I_KINETIC_G4x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G2xyz_vrr;
    Double I_KINETIC_H3x2y_G2xyz_vrr = PAY*I_KINETIC_G3xy_G2xyz_vrr+oned2z*I_KINETIC_F3x_G2xyz_vrr+oned2z*I_KINETIC_G3xy_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G2xyz_vrr-bdz*I_TWOBODYOVERLAP_F3x_G2xyz_vrr;
    Double I_KINETIC_H3xyz_G2xyz_vrr = PAZ*I_KINETIC_G3xy_G2xyz_vrr+oned2z*I_KINETIC_G3xy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G2xyz_vrr;
    Double I_KINETIC_H3x2z_G2xyz_vrr = PAZ*I_KINETIC_G3xz_G2xyz_vrr+oned2z*I_KINETIC_F3x_G2xyz_vrr+oned2z*I_KINETIC_G3xz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G2xyz_vrr-bdz*I_TWOBODYOVERLAP_F3x_G2xyz_vrr;
    Double I_KINETIC_H2x3y_G2xyz_vrr = PAX*I_KINETIC_Gx3y_G2xyz_vrr+oned2z*I_KINETIC_F3y_G2xyz_vrr+2*oned2z*I_KINETIC_Gx3y_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G2xyz_vrr-bdz*I_TWOBODYOVERLAP_F3y_G2xyz_vrr;
    Double I_KINETIC_H2x2yz_G2xyz_vrr = PAZ*I_KINETIC_G2x2y_G2xyz_vrr+oned2z*I_KINETIC_G2x2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G2xyz_vrr;
    Double I_KINETIC_H2xy2z_G2xyz_vrr = PAY*I_KINETIC_G2x2z_G2xyz_vrr+oned2z*I_KINETIC_G2x2z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G2xyz_vrr;
    Double I_KINETIC_H2x3z_G2xyz_vrr = PAX*I_KINETIC_Gx3z_G2xyz_vrr+oned2z*I_KINETIC_F3z_G2xyz_vrr+2*oned2z*I_KINETIC_Gx3z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G2xyz_vrr-bdz*I_TWOBODYOVERLAP_F3z_G2xyz_vrr;
    Double I_KINETIC_Hx4y_G2xyz_vrr = PAX*I_KINETIC_G4y_G2xyz_vrr+2*oned2z*I_KINETIC_G4y_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G2xyz_vrr;
    Double I_KINETIC_Hx3yz_G2xyz_vrr = PAZ*I_KINETIC_Gx3y_G2xyz_vrr+oned2z*I_KINETIC_Gx3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G2xyz_vrr;
    Double I_KINETIC_Hx2y2z_G2xyz_vrr = PAX*I_KINETIC_G2y2z_G2xyz_vrr+2*oned2z*I_KINETIC_G2y2z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G2xyz_vrr;
    Double I_KINETIC_Hxy3z_G2xyz_vrr = PAY*I_KINETIC_Gx3z_G2xyz_vrr+oned2z*I_KINETIC_Gx3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G2xyz_vrr;
    Double I_KINETIC_Hx4z_G2xyz_vrr = PAX*I_KINETIC_G4z_G2xyz_vrr+2*oned2z*I_KINETIC_G4z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G2xyz_vrr;
    Double I_KINETIC_H5y_G2xyz_vrr = PAY*I_KINETIC_G4y_G2xyz_vrr+4*oned2z*I_KINETIC_F3y_G2xyz_vrr+oned2z*I_KINETIC_G4y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G2xyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G2xyz_vrr;
    Double I_KINETIC_H4yz_G2xyz_vrr = PAZ*I_KINETIC_G4y_G2xyz_vrr+oned2z*I_KINETIC_G4y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G2xyz_vrr;
    Double I_KINETIC_H3y2z_G2xyz_vrr = PAZ*I_KINETIC_G3yz_G2xyz_vrr+oned2z*I_KINETIC_F3y_G2xyz_vrr+oned2z*I_KINETIC_G3yz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G2xyz_vrr-bdz*I_TWOBODYOVERLAP_F3y_G2xyz_vrr;
    Double I_KINETIC_H2y3z_G2xyz_vrr = PAY*I_KINETIC_Gy3z_G2xyz_vrr+oned2z*I_KINETIC_F3z_G2xyz_vrr+oned2z*I_KINETIC_Gy3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G2xyz_vrr-bdz*I_TWOBODYOVERLAP_F3z_G2xyz_vrr;
    Double I_KINETIC_Hy4z_G2xyz_vrr = PAY*I_KINETIC_G4z_G2xyz_vrr+oned2z*I_KINETIC_G4z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G2xyz_vrr;
    Double I_KINETIC_H5z_G2xyz_vrr = PAZ*I_KINETIC_G4z_G2xyz_vrr+4*oned2z*I_KINETIC_F3z_G2xyz_vrr+oned2z*I_KINETIC_G4z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G2xyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G2xyz_vrr;
    Double I_KINETIC_H5x_G2x2z_vrr = PAX*I_KINETIC_G4x_G2x2z_vrr+4*oned2z*I_KINETIC_F3x_G2x2z_vrr+2*oned2z*I_KINETIC_G4x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G2x2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G2x2z_vrr;
    Double I_KINETIC_H4xy_G2x2z_vrr = PAY*I_KINETIC_G4x_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G2x2z_vrr;
    Double I_KINETIC_H4xz_G2x2z_vrr = PAZ*I_KINETIC_G4x_G2x2z_vrr+2*oned2z*I_KINETIC_G4x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G2x2z_vrr;
    Double I_KINETIC_H3x2y_G2x2z_vrr = PAY*I_KINETIC_G3xy_G2x2z_vrr+oned2z*I_KINETIC_F3x_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G2x2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_G2x2z_vrr;
    Double I_KINETIC_H3xyz_G2x2z_vrr = PAZ*I_KINETIC_G3xy_G2x2z_vrr+2*oned2z*I_KINETIC_G3xy_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G2x2z_vrr;
    Double I_KINETIC_H3x2z_G2x2z_vrr = PAZ*I_KINETIC_G3xz_G2x2z_vrr+oned2z*I_KINETIC_F3x_G2x2z_vrr+2*oned2z*I_KINETIC_G3xz_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G2x2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_G2x2z_vrr;
    Double I_KINETIC_H2x3y_G2x2z_vrr = PAX*I_KINETIC_Gx3y_G2x2z_vrr+oned2z*I_KINETIC_F3y_G2x2z_vrr+2*oned2z*I_KINETIC_Gx3y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G2x2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_G2x2z_vrr;
    Double I_KINETIC_H2x2yz_G2x2z_vrr = PAZ*I_KINETIC_G2x2y_G2x2z_vrr+2*oned2z*I_KINETIC_G2x2y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G2x2z_vrr;
    Double I_KINETIC_H2xy2z_G2x2z_vrr = PAY*I_KINETIC_G2x2z_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G2x2z_vrr;
    Double I_KINETIC_H2x3z_G2x2z_vrr = PAX*I_KINETIC_Gx3z_G2x2z_vrr+oned2z*I_KINETIC_F3z_G2x2z_vrr+2*oned2z*I_KINETIC_Gx3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G2x2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_G2x2z_vrr;
    Double I_KINETIC_Hx4y_G2x2z_vrr = PAX*I_KINETIC_G4y_G2x2z_vrr+2*oned2z*I_KINETIC_G4y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G2x2z_vrr;
    Double I_KINETIC_Hx3yz_G2x2z_vrr = PAZ*I_KINETIC_Gx3y_G2x2z_vrr+2*oned2z*I_KINETIC_Gx3y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G2x2z_vrr;
    Double I_KINETIC_Hx2y2z_G2x2z_vrr = PAX*I_KINETIC_G2y2z_G2x2z_vrr+2*oned2z*I_KINETIC_G2y2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G2x2z_vrr;
    Double I_KINETIC_Hxy3z_G2x2z_vrr = PAY*I_KINETIC_Gx3z_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G2x2z_vrr;
    Double I_KINETIC_Hx4z_G2x2z_vrr = PAX*I_KINETIC_G4z_G2x2z_vrr+2*oned2z*I_KINETIC_G4z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G2x2z_vrr;
    Double I_KINETIC_H5y_G2x2z_vrr = PAY*I_KINETIC_G4y_G2x2z_vrr+4*oned2z*I_KINETIC_F3y_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G2x2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G2x2z_vrr;
    Double I_KINETIC_H4yz_G2x2z_vrr = PAZ*I_KINETIC_G4y_G2x2z_vrr+2*oned2z*I_KINETIC_G4y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G2x2z_vrr;
    Double I_KINETIC_H3y2z_G2x2z_vrr = PAZ*I_KINETIC_G3yz_G2x2z_vrr+oned2z*I_KINETIC_F3y_G2x2z_vrr+2*oned2z*I_KINETIC_G3yz_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G2x2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_G2x2z_vrr;
    Double I_KINETIC_H2y3z_G2x2z_vrr = PAY*I_KINETIC_Gy3z_G2x2z_vrr+oned2z*I_KINETIC_F3z_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G2x2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_G2x2z_vrr;
    Double I_KINETIC_Hy4z_G2x2z_vrr = PAY*I_KINETIC_G4z_G2x2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G2x2z_vrr;
    Double I_KINETIC_H5z_G2x2z_vrr = PAZ*I_KINETIC_G4z_G2x2z_vrr+4*oned2z*I_KINETIC_F3z_G2x2z_vrr+2*oned2z*I_KINETIC_G4z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G2x2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G2x2z_vrr;
    Double I_KINETIC_H5x_Gx3y_vrr = PAX*I_KINETIC_G4x_Gx3y_vrr+4*oned2z*I_KINETIC_F3x_Gx3y_vrr+oned2z*I_KINETIC_G4x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gx3y_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Gx3y_vrr;
    Double I_KINETIC_H4xy_Gx3y_vrr = PAY*I_KINETIC_G4x_Gx3y_vrr+3*oned2z*I_KINETIC_G4x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gx3y_vrr;
    Double I_KINETIC_H4xz_Gx3y_vrr = PAZ*I_KINETIC_G4x_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gx3y_vrr;
    Double I_KINETIC_H3x2y_Gx3y_vrr = PAY*I_KINETIC_G3xy_Gx3y_vrr+oned2z*I_KINETIC_F3x_Gx3y_vrr+3*oned2z*I_KINETIC_G3xy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gx3y_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gx3y_vrr;
    Double I_KINETIC_H3xyz_Gx3y_vrr = PAZ*I_KINETIC_G3xy_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gx3y_vrr;
    Double I_KINETIC_H3x2z_Gx3y_vrr = PAZ*I_KINETIC_G3xz_Gx3y_vrr+oned2z*I_KINETIC_F3x_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gx3y_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gx3y_vrr;
    Double I_KINETIC_H2x3y_Gx3y_vrr = PAX*I_KINETIC_Gx3y_Gx3y_vrr+oned2z*I_KINETIC_F3y_Gx3y_vrr+oned2z*I_KINETIC_Gx3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gx3y_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gx3y_vrr;
    Double I_KINETIC_H2x2yz_Gx3y_vrr = PAZ*I_KINETIC_G2x2y_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gx3y_vrr;
    Double I_KINETIC_H2xy2z_Gx3y_vrr = PAY*I_KINETIC_G2x2z_Gx3y_vrr+3*oned2z*I_KINETIC_G2x2z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gx3y_vrr;
    Double I_KINETIC_H2x3z_Gx3y_vrr = PAX*I_KINETIC_Gx3z_Gx3y_vrr+oned2z*I_KINETIC_F3z_Gx3y_vrr+oned2z*I_KINETIC_Gx3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gx3y_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gx3y_vrr;
    Double I_KINETIC_Hx4y_Gx3y_vrr = PAX*I_KINETIC_G4y_Gx3y_vrr+oned2z*I_KINETIC_G4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gx3y_vrr;
    Double I_KINETIC_Hx3yz_Gx3y_vrr = PAZ*I_KINETIC_Gx3y_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gx3y_vrr;
    Double I_KINETIC_Hx2y2z_Gx3y_vrr = PAX*I_KINETIC_G2y2z_Gx3y_vrr+oned2z*I_KINETIC_G2y2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gx3y_vrr;
    Double I_KINETIC_Hxy3z_Gx3y_vrr = PAY*I_KINETIC_Gx3z_Gx3y_vrr+3*oned2z*I_KINETIC_Gx3z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gx3y_vrr;
    Double I_KINETIC_Hx4z_Gx3y_vrr = PAX*I_KINETIC_G4z_Gx3y_vrr+oned2z*I_KINETIC_G4z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gx3y_vrr;
    Double I_KINETIC_H5y_Gx3y_vrr = PAY*I_KINETIC_G4y_Gx3y_vrr+4*oned2z*I_KINETIC_F3y_Gx3y_vrr+3*oned2z*I_KINETIC_G4y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gx3y_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Gx3y_vrr;
    Double I_KINETIC_H4yz_Gx3y_vrr = PAZ*I_KINETIC_G4y_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gx3y_vrr;
    Double I_KINETIC_H3y2z_Gx3y_vrr = PAZ*I_KINETIC_G3yz_Gx3y_vrr+oned2z*I_KINETIC_F3y_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gx3y_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gx3y_vrr;
    Double I_KINETIC_H2y3z_Gx3y_vrr = PAY*I_KINETIC_Gy3z_Gx3y_vrr+oned2z*I_KINETIC_F3z_Gx3y_vrr+3*oned2z*I_KINETIC_Gy3z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gx3y_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gx3y_vrr;
    Double I_KINETIC_Hy4z_Gx3y_vrr = PAY*I_KINETIC_G4z_Gx3y_vrr+3*oned2z*I_KINETIC_G4z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gx3y_vrr;
    Double I_KINETIC_H5z_Gx3y_vrr = PAZ*I_KINETIC_G4z_Gx3y_vrr+4*oned2z*I_KINETIC_F3z_Gx3y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gx3y_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Gx3y_vrr;
    Double I_KINETIC_H5x_Gx2yz_vrr = PAX*I_KINETIC_G4x_Gx2yz_vrr+4*oned2z*I_KINETIC_F3x_Gx2yz_vrr+oned2z*I_KINETIC_G4x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gx2yz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr;
    Double I_KINETIC_H4xy_Gx2yz_vrr = PAY*I_KINETIC_G4x_Gx2yz_vrr+2*oned2z*I_KINETIC_G4x_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gx2yz_vrr;
    Double I_KINETIC_H4xz_Gx2yz_vrr = PAZ*I_KINETIC_G4x_Gx2yz_vrr+oned2z*I_KINETIC_G4x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gx2yz_vrr;
    Double I_KINETIC_H3x2y_Gx2yz_vrr = PAY*I_KINETIC_G3xy_Gx2yz_vrr+oned2z*I_KINETIC_F3x_Gx2yz_vrr+2*oned2z*I_KINETIC_G3xy_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gx2yz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr;
    Double I_KINETIC_H3xyz_Gx2yz_vrr = PAZ*I_KINETIC_G3xy_Gx2yz_vrr+oned2z*I_KINETIC_G3xy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gx2yz_vrr;
    Double I_KINETIC_H3x2z_Gx2yz_vrr = PAZ*I_KINETIC_G3xz_Gx2yz_vrr+oned2z*I_KINETIC_F3x_Gx2yz_vrr+oned2z*I_KINETIC_G3xz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gx2yz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gx2yz_vrr;
    Double I_KINETIC_H2x3y_Gx2yz_vrr = PAX*I_KINETIC_Gx3y_Gx2yz_vrr+oned2z*I_KINETIC_F3y_Gx2yz_vrr+oned2z*I_KINETIC_Gx3y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gx2yz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr;
    Double I_KINETIC_H2x2yz_Gx2yz_vrr = PAZ*I_KINETIC_G2x2y_Gx2yz_vrr+oned2z*I_KINETIC_G2x2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gx2yz_vrr;
    Double I_KINETIC_H2xy2z_Gx2yz_vrr = PAY*I_KINETIC_G2x2z_Gx2yz_vrr+2*oned2z*I_KINETIC_G2x2z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gx2yz_vrr;
    Double I_KINETIC_H2x3z_Gx2yz_vrr = PAX*I_KINETIC_Gx3z_Gx2yz_vrr+oned2z*I_KINETIC_F3z_Gx2yz_vrr+oned2z*I_KINETIC_Gx3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gx2yz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr;
    Double I_KINETIC_Hx4y_Gx2yz_vrr = PAX*I_KINETIC_G4y_Gx2yz_vrr+oned2z*I_KINETIC_G4y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gx2yz_vrr;
    Double I_KINETIC_Hx3yz_Gx2yz_vrr = PAZ*I_KINETIC_Gx3y_Gx2yz_vrr+oned2z*I_KINETIC_Gx3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gx2yz_vrr;
    Double I_KINETIC_Hx2y2z_Gx2yz_vrr = PAX*I_KINETIC_G2y2z_Gx2yz_vrr+oned2z*I_KINETIC_G2y2z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gx2yz_vrr;
    Double I_KINETIC_Hxy3z_Gx2yz_vrr = PAY*I_KINETIC_Gx3z_Gx2yz_vrr+2*oned2z*I_KINETIC_Gx3z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gx2yz_vrr;
    Double I_KINETIC_Hx4z_Gx2yz_vrr = PAX*I_KINETIC_G4z_Gx2yz_vrr+oned2z*I_KINETIC_G4z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gx2yz_vrr;
    Double I_KINETIC_H5y_Gx2yz_vrr = PAY*I_KINETIC_G4y_Gx2yz_vrr+4*oned2z*I_KINETIC_F3y_Gx2yz_vrr+2*oned2z*I_KINETIC_G4y_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gx2yz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr;
    Double I_KINETIC_H4yz_Gx2yz_vrr = PAZ*I_KINETIC_G4y_Gx2yz_vrr+oned2z*I_KINETIC_G4y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gx2yz_vrr;
    Double I_KINETIC_H3y2z_Gx2yz_vrr = PAZ*I_KINETIC_G3yz_Gx2yz_vrr+oned2z*I_KINETIC_F3y_Gx2yz_vrr+oned2z*I_KINETIC_G3yz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gx2yz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gx2yz_vrr;
    Double I_KINETIC_H2y3z_Gx2yz_vrr = PAY*I_KINETIC_Gy3z_Gx2yz_vrr+oned2z*I_KINETIC_F3z_Gx2yz_vrr+2*oned2z*I_KINETIC_Gy3z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gx2yz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr;
    Double I_KINETIC_Hy4z_Gx2yz_vrr = PAY*I_KINETIC_G4z_Gx2yz_vrr+2*oned2z*I_KINETIC_G4z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gx2yz_vrr;
    Double I_KINETIC_H5z_Gx2yz_vrr = PAZ*I_KINETIC_G4z_Gx2yz_vrr+4*oned2z*I_KINETIC_F3z_Gx2yz_vrr+oned2z*I_KINETIC_G4z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gx2yz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Gx2yz_vrr;
    Double I_KINETIC_H5x_Gxy2z_vrr = PAX*I_KINETIC_G4x_Gxy2z_vrr+4*oned2z*I_KINETIC_F3x_Gxy2z_vrr+oned2z*I_KINETIC_G4x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gxy2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr;
    Double I_KINETIC_H4xy_Gxy2z_vrr = PAY*I_KINETIC_G4x_Gxy2z_vrr+oned2z*I_KINETIC_G4x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gxy2z_vrr;
    Double I_KINETIC_H4xz_Gxy2z_vrr = PAZ*I_KINETIC_G4x_Gxy2z_vrr+2*oned2z*I_KINETIC_G4x_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gxy2z_vrr;
    Double I_KINETIC_H3x2y_Gxy2z_vrr = PAY*I_KINETIC_G3xy_Gxy2z_vrr+oned2z*I_KINETIC_F3x_Gxy2z_vrr+oned2z*I_KINETIC_G3xy_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gxy2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr;
    Double I_KINETIC_H3xyz_Gxy2z_vrr = PAZ*I_KINETIC_G3xy_Gxy2z_vrr+2*oned2z*I_KINETIC_G3xy_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gxy2z_vrr;
    Double I_KINETIC_H3x2z_Gxy2z_vrr = PAZ*I_KINETIC_G3xz_Gxy2z_vrr+oned2z*I_KINETIC_F3x_Gxy2z_vrr+2*oned2z*I_KINETIC_G3xz_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gxy2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gxy2z_vrr;
    Double I_KINETIC_H2x3y_Gxy2z_vrr = PAX*I_KINETIC_Gx3y_Gxy2z_vrr+oned2z*I_KINETIC_F3y_Gxy2z_vrr+oned2z*I_KINETIC_Gx3y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gxy2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr;
    Double I_KINETIC_H2x2yz_Gxy2z_vrr = PAZ*I_KINETIC_G2x2y_Gxy2z_vrr+2*oned2z*I_KINETIC_G2x2y_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gxy2z_vrr;
    Double I_KINETIC_H2xy2z_Gxy2z_vrr = PAY*I_KINETIC_G2x2z_Gxy2z_vrr+oned2z*I_KINETIC_G2x2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gxy2z_vrr;
    Double I_KINETIC_H2x3z_Gxy2z_vrr = PAX*I_KINETIC_Gx3z_Gxy2z_vrr+oned2z*I_KINETIC_F3z_Gxy2z_vrr+oned2z*I_KINETIC_Gx3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gxy2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr;
    Double I_KINETIC_Hx4y_Gxy2z_vrr = PAX*I_KINETIC_G4y_Gxy2z_vrr+oned2z*I_KINETIC_G4y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gxy2z_vrr;
    Double I_KINETIC_Hx3yz_Gxy2z_vrr = PAZ*I_KINETIC_Gx3y_Gxy2z_vrr+2*oned2z*I_KINETIC_Gx3y_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gxy2z_vrr;
    Double I_KINETIC_Hx2y2z_Gxy2z_vrr = PAX*I_KINETIC_G2y2z_Gxy2z_vrr+oned2z*I_KINETIC_G2y2z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gxy2z_vrr;
    Double I_KINETIC_Hxy3z_Gxy2z_vrr = PAY*I_KINETIC_Gx3z_Gxy2z_vrr+oned2z*I_KINETIC_Gx3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gxy2z_vrr;
    Double I_KINETIC_Hx4z_Gxy2z_vrr = PAX*I_KINETIC_G4z_Gxy2z_vrr+oned2z*I_KINETIC_G4z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gxy2z_vrr;
    Double I_KINETIC_H5y_Gxy2z_vrr = PAY*I_KINETIC_G4y_Gxy2z_vrr+4*oned2z*I_KINETIC_F3y_Gxy2z_vrr+oned2z*I_KINETIC_G4y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gxy2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr;
    Double I_KINETIC_H4yz_Gxy2z_vrr = PAZ*I_KINETIC_G4y_Gxy2z_vrr+2*oned2z*I_KINETIC_G4y_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gxy2z_vrr;
    Double I_KINETIC_H3y2z_Gxy2z_vrr = PAZ*I_KINETIC_G3yz_Gxy2z_vrr+oned2z*I_KINETIC_F3y_Gxy2z_vrr+2*oned2z*I_KINETIC_G3yz_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gxy2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gxy2z_vrr;
    Double I_KINETIC_H2y3z_Gxy2z_vrr = PAY*I_KINETIC_Gy3z_Gxy2z_vrr+oned2z*I_KINETIC_F3z_Gxy2z_vrr+oned2z*I_KINETIC_Gy3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gxy2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr;
    Double I_KINETIC_Hy4z_Gxy2z_vrr = PAY*I_KINETIC_G4z_Gxy2z_vrr+oned2z*I_KINETIC_G4z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gxy2z_vrr;
    Double I_KINETIC_H5z_Gxy2z_vrr = PAZ*I_KINETIC_G4z_Gxy2z_vrr+4*oned2z*I_KINETIC_F3z_Gxy2z_vrr+2*oned2z*I_KINETIC_G4z_Fxyz_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gxy2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Gxy2z_vrr;
    Double I_KINETIC_H5x_Gx3z_vrr = PAX*I_KINETIC_G4x_Gx3z_vrr+4*oned2z*I_KINETIC_F3x_Gx3z_vrr+oned2z*I_KINETIC_G4x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gx3z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Gx3z_vrr;
    Double I_KINETIC_H4xy_Gx3z_vrr = PAY*I_KINETIC_G4x_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gx3z_vrr;
    Double I_KINETIC_H4xz_Gx3z_vrr = PAZ*I_KINETIC_G4x_Gx3z_vrr+3*oned2z*I_KINETIC_G4x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gx3z_vrr;
    Double I_KINETIC_H3x2y_Gx3z_vrr = PAY*I_KINETIC_G3xy_Gx3z_vrr+oned2z*I_KINETIC_F3x_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gx3z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gx3z_vrr;
    Double I_KINETIC_H3xyz_Gx3z_vrr = PAZ*I_KINETIC_G3xy_Gx3z_vrr+3*oned2z*I_KINETIC_G3xy_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gx3z_vrr;
    Double I_KINETIC_H3x2z_Gx3z_vrr = PAZ*I_KINETIC_G3xz_Gx3z_vrr+oned2z*I_KINETIC_F3x_Gx3z_vrr+3*oned2z*I_KINETIC_G3xz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gx3z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gx3z_vrr;
    Double I_KINETIC_H2x3y_Gx3z_vrr = PAX*I_KINETIC_Gx3y_Gx3z_vrr+oned2z*I_KINETIC_F3y_Gx3z_vrr+oned2z*I_KINETIC_Gx3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gx3z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gx3z_vrr;
    Double I_KINETIC_H2x2yz_Gx3z_vrr = PAZ*I_KINETIC_G2x2y_Gx3z_vrr+3*oned2z*I_KINETIC_G2x2y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gx3z_vrr;
    Double I_KINETIC_H2xy2z_Gx3z_vrr = PAY*I_KINETIC_G2x2z_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gx3z_vrr;
    Double I_KINETIC_H2x3z_Gx3z_vrr = PAX*I_KINETIC_Gx3z_Gx3z_vrr+oned2z*I_KINETIC_F3z_Gx3z_vrr+oned2z*I_KINETIC_Gx3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gx3z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gx3z_vrr;
    Double I_KINETIC_Hx4y_Gx3z_vrr = PAX*I_KINETIC_G4y_Gx3z_vrr+oned2z*I_KINETIC_G4y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gx3z_vrr;
    Double I_KINETIC_Hx3yz_Gx3z_vrr = PAZ*I_KINETIC_Gx3y_Gx3z_vrr+3*oned2z*I_KINETIC_Gx3y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gx3z_vrr;
    Double I_KINETIC_Hx2y2z_Gx3z_vrr = PAX*I_KINETIC_G2y2z_Gx3z_vrr+oned2z*I_KINETIC_G2y2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gx3z_vrr;
    Double I_KINETIC_Hxy3z_Gx3z_vrr = PAY*I_KINETIC_Gx3z_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gx3z_vrr;
    Double I_KINETIC_Hx4z_Gx3z_vrr = PAX*I_KINETIC_G4z_Gx3z_vrr+oned2z*I_KINETIC_G4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gx3z_vrr;
    Double I_KINETIC_H5y_Gx3z_vrr = PAY*I_KINETIC_G4y_Gx3z_vrr+4*oned2z*I_KINETIC_F3y_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gx3z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Gx3z_vrr;
    Double I_KINETIC_H4yz_Gx3z_vrr = PAZ*I_KINETIC_G4y_Gx3z_vrr+3*oned2z*I_KINETIC_G4y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gx3z_vrr;
    Double I_KINETIC_H3y2z_Gx3z_vrr = PAZ*I_KINETIC_G3yz_Gx3z_vrr+oned2z*I_KINETIC_F3y_Gx3z_vrr+3*oned2z*I_KINETIC_G3yz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gx3z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gx3z_vrr;
    Double I_KINETIC_H2y3z_Gx3z_vrr = PAY*I_KINETIC_Gy3z_Gx3z_vrr+oned2z*I_KINETIC_F3z_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gx3z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gx3z_vrr;
    Double I_KINETIC_Hy4z_Gx3z_vrr = PAY*I_KINETIC_G4z_Gx3z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gx3z_vrr;
    Double I_KINETIC_H5z_Gx3z_vrr = PAZ*I_KINETIC_G4z_Gx3z_vrr+4*oned2z*I_KINETIC_F3z_Gx3z_vrr+3*oned2z*I_KINETIC_G4z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gx3z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Gx3z_vrr;
    Double I_KINETIC_H5x_G4y_vrr = PAX*I_KINETIC_G4x_G4y_vrr+4*oned2z*I_KINETIC_F3x_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G4y_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G4y_vrr;
    Double I_KINETIC_H4xy_G4y_vrr = PAY*I_KINETIC_G4x_G4y_vrr+4*oned2z*I_KINETIC_G4x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G4y_vrr;
    Double I_KINETIC_H4xz_G4y_vrr = PAZ*I_KINETIC_G4x_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G4y_vrr;
    Double I_KINETIC_H3x2y_G4y_vrr = PAY*I_KINETIC_G3xy_G4y_vrr+oned2z*I_KINETIC_F3x_G4y_vrr+4*oned2z*I_KINETIC_G3xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G4y_vrr-bdz*I_TWOBODYOVERLAP_F3x_G4y_vrr;
    Double I_KINETIC_H3xyz_G4y_vrr = PAZ*I_KINETIC_G3xy_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G4y_vrr;
    Double I_KINETIC_H3x2z_G4y_vrr = PAZ*I_KINETIC_G3xz_G4y_vrr+oned2z*I_KINETIC_F3x_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G4y_vrr-bdz*I_TWOBODYOVERLAP_F3x_G4y_vrr;
    Double I_KINETIC_H2x3y_G4y_vrr = PAX*I_KINETIC_Gx3y_G4y_vrr+oned2z*I_KINETIC_F3y_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G4y_vrr-bdz*I_TWOBODYOVERLAP_F3y_G4y_vrr;
    Double I_KINETIC_H2x2yz_G4y_vrr = PAZ*I_KINETIC_G2x2y_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G4y_vrr;
    Double I_KINETIC_H2xy2z_G4y_vrr = PAY*I_KINETIC_G2x2z_G4y_vrr+4*oned2z*I_KINETIC_G2x2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G4y_vrr;
    Double I_KINETIC_H2x3z_G4y_vrr = PAX*I_KINETIC_Gx3z_G4y_vrr+oned2z*I_KINETIC_F3z_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G4y_vrr-bdz*I_TWOBODYOVERLAP_F3z_G4y_vrr;
    Double I_KINETIC_Hx4y_G4y_vrr = PAX*I_KINETIC_G4y_G4y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G4y_vrr;
    Double I_KINETIC_Hx3yz_G4y_vrr = PAZ*I_KINETIC_Gx3y_G4y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G4y_vrr;
    Double I_KINETIC_Hx2y2z_G4y_vrr = PAX*I_KINETIC_G2y2z_G4y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G4y_vrr;
    Double I_KINETIC_Hxy3z_G4y_vrr = PAY*I_KINETIC_Gx3z_G4y_vrr+4*oned2z*I_KINETIC_Gx3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G4y_vrr;
    Double I_KINETIC_Hx4z_G4y_vrr = PAX*I_KINETIC_G4z_G4y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G4y_vrr;
    Double I_KINETIC_H5y_G4y_vrr = PAY*I_KINETIC_G4y_G4y_vrr+4*oned2z*I_KINETIC_F3y_G4y_vrr+4*oned2z*I_KINETIC_G4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G4y_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G4y_vrr;
    Double I_KINETIC_H4yz_G4y_vrr = PAZ*I_KINETIC_G4y_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G4y_vrr;
    Double I_KINETIC_H3y2z_G4y_vrr = PAZ*I_KINETIC_G3yz_G4y_vrr+oned2z*I_KINETIC_F3y_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G4y_vrr-bdz*I_TWOBODYOVERLAP_F3y_G4y_vrr;
    Double I_KINETIC_H2y3z_G4y_vrr = PAY*I_KINETIC_Gy3z_G4y_vrr+oned2z*I_KINETIC_F3z_G4y_vrr+4*oned2z*I_KINETIC_Gy3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G4y_vrr-bdz*I_TWOBODYOVERLAP_F3z_G4y_vrr;
    Double I_KINETIC_Hy4z_G4y_vrr = PAY*I_KINETIC_G4z_G4y_vrr+4*oned2z*I_KINETIC_G4z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G4y_vrr;
    Double I_KINETIC_H5z_G4y_vrr = PAZ*I_KINETIC_G4z_G4y_vrr+4*oned2z*I_KINETIC_F3z_G4y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G4y_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G4y_vrr;
    Double I_KINETIC_H5x_G3yz_vrr = PAX*I_KINETIC_G4x_G3yz_vrr+4*oned2z*I_KINETIC_F3x_G3yz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G3yz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G3yz_vrr;
    Double I_KINETIC_H4xy_G3yz_vrr = PAY*I_KINETIC_G4x_G3yz_vrr+3*oned2z*I_KINETIC_G4x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G3yz_vrr;
    Double I_KINETIC_H4xz_G3yz_vrr = PAZ*I_KINETIC_G4x_G3yz_vrr+oned2z*I_KINETIC_G4x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G3yz_vrr;
    Double I_KINETIC_H3x2y_G3yz_vrr = PAY*I_KINETIC_G3xy_G3yz_vrr+oned2z*I_KINETIC_F3x_G3yz_vrr+3*oned2z*I_KINETIC_G3xy_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G3yz_vrr-bdz*I_TWOBODYOVERLAP_F3x_G3yz_vrr;
    Double I_KINETIC_H3xyz_G3yz_vrr = PAZ*I_KINETIC_G3xy_G3yz_vrr+oned2z*I_KINETIC_G3xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G3yz_vrr;
    Double I_KINETIC_H3x2z_G3yz_vrr = PAZ*I_KINETIC_G3xz_G3yz_vrr+oned2z*I_KINETIC_F3x_G3yz_vrr+oned2z*I_KINETIC_G3xz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G3yz_vrr-bdz*I_TWOBODYOVERLAP_F3x_G3yz_vrr;
    Double I_KINETIC_H2x3y_G3yz_vrr = PAX*I_KINETIC_Gx3y_G3yz_vrr+oned2z*I_KINETIC_F3y_G3yz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G3yz_vrr-bdz*I_TWOBODYOVERLAP_F3y_G3yz_vrr;
    Double I_KINETIC_H2x2yz_G3yz_vrr = PAZ*I_KINETIC_G2x2y_G3yz_vrr+oned2z*I_KINETIC_G2x2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G3yz_vrr;
    Double I_KINETIC_H2xy2z_G3yz_vrr = PAY*I_KINETIC_G2x2z_G3yz_vrr+3*oned2z*I_KINETIC_G2x2z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G3yz_vrr;
    Double I_KINETIC_H2x3z_G3yz_vrr = PAX*I_KINETIC_Gx3z_G3yz_vrr+oned2z*I_KINETIC_F3z_G3yz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G3yz_vrr-bdz*I_TWOBODYOVERLAP_F3z_G3yz_vrr;
    Double I_KINETIC_Hx4y_G3yz_vrr = PAX*I_KINETIC_G4y_G3yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G3yz_vrr;
    Double I_KINETIC_Hx3yz_G3yz_vrr = PAZ*I_KINETIC_Gx3y_G3yz_vrr+oned2z*I_KINETIC_Gx3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G3yz_vrr;
    Double I_KINETIC_Hx2y2z_G3yz_vrr = PAX*I_KINETIC_G2y2z_G3yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G3yz_vrr;
    Double I_KINETIC_Hxy3z_G3yz_vrr = PAY*I_KINETIC_Gx3z_G3yz_vrr+3*oned2z*I_KINETIC_Gx3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G3yz_vrr;
    Double I_KINETIC_Hx4z_G3yz_vrr = PAX*I_KINETIC_G4z_G3yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G3yz_vrr;
    Double I_KINETIC_H5y_G3yz_vrr = PAY*I_KINETIC_G4y_G3yz_vrr+4*oned2z*I_KINETIC_F3y_G3yz_vrr+3*oned2z*I_KINETIC_G4y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G3yz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G3yz_vrr;
    Double I_KINETIC_H4yz_G3yz_vrr = PAZ*I_KINETIC_G4y_G3yz_vrr+oned2z*I_KINETIC_G4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G3yz_vrr;
    Double I_KINETIC_H3y2z_G3yz_vrr = PAZ*I_KINETIC_G3yz_G3yz_vrr+oned2z*I_KINETIC_F3y_G3yz_vrr+oned2z*I_KINETIC_G3yz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G3yz_vrr-bdz*I_TWOBODYOVERLAP_F3y_G3yz_vrr;
    Double I_KINETIC_H2y3z_G3yz_vrr = PAY*I_KINETIC_Gy3z_G3yz_vrr+oned2z*I_KINETIC_F3z_G3yz_vrr+3*oned2z*I_KINETIC_Gy3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G3yz_vrr-bdz*I_TWOBODYOVERLAP_F3z_G3yz_vrr;
    Double I_KINETIC_Hy4z_G3yz_vrr = PAY*I_KINETIC_G4z_G3yz_vrr+3*oned2z*I_KINETIC_G4z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G3yz_vrr;
    Double I_KINETIC_H5z_G3yz_vrr = PAZ*I_KINETIC_G4z_G3yz_vrr+4*oned2z*I_KINETIC_F3z_G3yz_vrr+oned2z*I_KINETIC_G4z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G3yz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G3yz_vrr;
    Double I_KINETIC_H5x_G2y2z_vrr = PAX*I_KINETIC_G4x_G2y2z_vrr+4*oned2z*I_KINETIC_F3x_G2y2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G2y2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G2y2z_vrr;
    Double I_KINETIC_H4xy_G2y2z_vrr = PAY*I_KINETIC_G4x_G2y2z_vrr+2*oned2z*I_KINETIC_G4x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G2y2z_vrr;
    Double I_KINETIC_H4xz_G2y2z_vrr = PAZ*I_KINETIC_G4x_G2y2z_vrr+2*oned2z*I_KINETIC_G4x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G2y2z_vrr;
    Double I_KINETIC_H3x2y_G2y2z_vrr = PAY*I_KINETIC_G3xy_G2y2z_vrr+oned2z*I_KINETIC_F3x_G2y2z_vrr+2*oned2z*I_KINETIC_G3xy_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G2y2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_G2y2z_vrr;
    Double I_KINETIC_H3xyz_G2y2z_vrr = PAZ*I_KINETIC_G3xy_G2y2z_vrr+2*oned2z*I_KINETIC_G3xy_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G2y2z_vrr;
    Double I_KINETIC_H3x2z_G2y2z_vrr = PAZ*I_KINETIC_G3xz_G2y2z_vrr+oned2z*I_KINETIC_F3x_G2y2z_vrr+2*oned2z*I_KINETIC_G3xz_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G2y2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_G2y2z_vrr;
    Double I_KINETIC_H2x3y_G2y2z_vrr = PAX*I_KINETIC_Gx3y_G2y2z_vrr+oned2z*I_KINETIC_F3y_G2y2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G2y2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_G2y2z_vrr;
    Double I_KINETIC_H2x2yz_G2y2z_vrr = PAZ*I_KINETIC_G2x2y_G2y2z_vrr+2*oned2z*I_KINETIC_G2x2y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G2y2z_vrr;
    Double I_KINETIC_H2xy2z_G2y2z_vrr = PAY*I_KINETIC_G2x2z_G2y2z_vrr+2*oned2z*I_KINETIC_G2x2z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G2y2z_vrr;
    Double I_KINETIC_H2x3z_G2y2z_vrr = PAX*I_KINETIC_Gx3z_G2y2z_vrr+oned2z*I_KINETIC_F3z_G2y2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G2y2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_G2y2z_vrr;
    Double I_KINETIC_Hx4y_G2y2z_vrr = PAX*I_KINETIC_G4y_G2y2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G2y2z_vrr;
    Double I_KINETIC_Hx3yz_G2y2z_vrr = PAZ*I_KINETIC_Gx3y_G2y2z_vrr+2*oned2z*I_KINETIC_Gx3y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G2y2z_vrr;
    Double I_KINETIC_Hx2y2z_G2y2z_vrr = PAX*I_KINETIC_G2y2z_G2y2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G2y2z_vrr;
    Double I_KINETIC_Hxy3z_G2y2z_vrr = PAY*I_KINETIC_Gx3z_G2y2z_vrr+2*oned2z*I_KINETIC_Gx3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G2y2z_vrr;
    Double I_KINETIC_Hx4z_G2y2z_vrr = PAX*I_KINETIC_G4z_G2y2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G2y2z_vrr;
    Double I_KINETIC_H5y_G2y2z_vrr = PAY*I_KINETIC_G4y_G2y2z_vrr+4*oned2z*I_KINETIC_F3y_G2y2z_vrr+2*oned2z*I_KINETIC_G4y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G2y2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G2y2z_vrr;
    Double I_KINETIC_H4yz_G2y2z_vrr = PAZ*I_KINETIC_G4y_G2y2z_vrr+2*oned2z*I_KINETIC_G4y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G2y2z_vrr;
    Double I_KINETIC_H3y2z_G2y2z_vrr = PAZ*I_KINETIC_G3yz_G2y2z_vrr+oned2z*I_KINETIC_F3y_G2y2z_vrr+2*oned2z*I_KINETIC_G3yz_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G2y2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_G2y2z_vrr;
    Double I_KINETIC_H2y3z_G2y2z_vrr = PAY*I_KINETIC_Gy3z_G2y2z_vrr+oned2z*I_KINETIC_F3z_G2y2z_vrr+2*oned2z*I_KINETIC_Gy3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G2y2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_G2y2z_vrr;
    Double I_KINETIC_Hy4z_G2y2z_vrr = PAY*I_KINETIC_G4z_G2y2z_vrr+2*oned2z*I_KINETIC_G4z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G2y2z_vrr;
    Double I_KINETIC_H5z_G2y2z_vrr = PAZ*I_KINETIC_G4z_G2y2z_vrr+4*oned2z*I_KINETIC_F3z_G2y2z_vrr+2*oned2z*I_KINETIC_G4z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G2y2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G2y2z_vrr;
    Double I_KINETIC_H5x_Gy3z_vrr = PAX*I_KINETIC_G4x_Gy3z_vrr+4*oned2z*I_KINETIC_F3x_Gy3z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gy3z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Gy3z_vrr;
    Double I_KINETIC_H4xy_Gy3z_vrr = PAY*I_KINETIC_G4x_Gy3z_vrr+oned2z*I_KINETIC_G4x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gy3z_vrr;
    Double I_KINETIC_H4xz_Gy3z_vrr = PAZ*I_KINETIC_G4x_Gy3z_vrr+3*oned2z*I_KINETIC_G4x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gy3z_vrr;
    Double I_KINETIC_H3x2y_Gy3z_vrr = PAY*I_KINETIC_G3xy_Gy3z_vrr+oned2z*I_KINETIC_F3x_Gy3z_vrr+oned2z*I_KINETIC_G3xy_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gy3z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gy3z_vrr;
    Double I_KINETIC_H3xyz_Gy3z_vrr = PAZ*I_KINETIC_G3xy_Gy3z_vrr+3*oned2z*I_KINETIC_G3xy_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gy3z_vrr;
    Double I_KINETIC_H3x2z_Gy3z_vrr = PAZ*I_KINETIC_G3xz_Gy3z_vrr+oned2z*I_KINETIC_F3x_Gy3z_vrr+3*oned2z*I_KINETIC_G3xz_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gy3z_vrr-bdz*I_TWOBODYOVERLAP_F3x_Gy3z_vrr;
    Double I_KINETIC_H2x3y_Gy3z_vrr = PAX*I_KINETIC_Gx3y_Gy3z_vrr+oned2z*I_KINETIC_F3y_Gy3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gy3z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gy3z_vrr;
    Double I_KINETIC_H2x2yz_Gy3z_vrr = PAZ*I_KINETIC_G2x2y_Gy3z_vrr+3*oned2z*I_KINETIC_G2x2y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gy3z_vrr;
    Double I_KINETIC_H2xy2z_Gy3z_vrr = PAY*I_KINETIC_G2x2z_Gy3z_vrr+oned2z*I_KINETIC_G2x2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gy3z_vrr;
    Double I_KINETIC_H2x3z_Gy3z_vrr = PAX*I_KINETIC_Gx3z_Gy3z_vrr+oned2z*I_KINETIC_F3z_Gy3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gy3z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gy3z_vrr;
    Double I_KINETIC_Hx4y_Gy3z_vrr = PAX*I_KINETIC_G4y_Gy3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gy3z_vrr;
    Double I_KINETIC_Hx3yz_Gy3z_vrr = PAZ*I_KINETIC_Gx3y_Gy3z_vrr+3*oned2z*I_KINETIC_Gx3y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gy3z_vrr;
    Double I_KINETIC_Hx2y2z_Gy3z_vrr = PAX*I_KINETIC_G2y2z_Gy3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gy3z_vrr;
    Double I_KINETIC_Hxy3z_Gy3z_vrr = PAY*I_KINETIC_Gx3z_Gy3z_vrr+oned2z*I_KINETIC_Gx3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gy3z_vrr;
    Double I_KINETIC_Hx4z_Gy3z_vrr = PAX*I_KINETIC_G4z_Gy3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gy3z_vrr;
    Double I_KINETIC_H5y_Gy3z_vrr = PAY*I_KINETIC_G4y_Gy3z_vrr+4*oned2z*I_KINETIC_F3y_Gy3z_vrr+oned2z*I_KINETIC_G4y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gy3z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Gy3z_vrr;
    Double I_KINETIC_H4yz_Gy3z_vrr = PAZ*I_KINETIC_G4y_Gy3z_vrr+3*oned2z*I_KINETIC_G4y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gy3z_vrr;
    Double I_KINETIC_H3y2z_Gy3z_vrr = PAZ*I_KINETIC_G3yz_Gy3z_vrr+oned2z*I_KINETIC_F3y_Gy3z_vrr+3*oned2z*I_KINETIC_G3yz_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gy3z_vrr-bdz*I_TWOBODYOVERLAP_F3y_Gy3z_vrr;
    Double I_KINETIC_H2y3z_Gy3z_vrr = PAY*I_KINETIC_Gy3z_Gy3z_vrr+oned2z*I_KINETIC_F3z_Gy3z_vrr+oned2z*I_KINETIC_Gy3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gy3z_vrr-bdz*I_TWOBODYOVERLAP_F3z_Gy3z_vrr;
    Double I_KINETIC_Hy4z_Gy3z_vrr = PAY*I_KINETIC_G4z_Gy3z_vrr+oned2z*I_KINETIC_G4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gy3z_vrr;
    Double I_KINETIC_H5z_Gy3z_vrr = PAZ*I_KINETIC_G4z_Gy3z_vrr+4*oned2z*I_KINETIC_F3z_Gy3z_vrr+3*oned2z*I_KINETIC_G4z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gy3z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Gy3z_vrr;
    Double I_KINETIC_H5x_G4z_vrr = PAX*I_KINETIC_G4x_G4z_vrr+4*oned2z*I_KINETIC_F3x_G4z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G4z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_G4z_vrr;
    Double I_KINETIC_H4xy_G4z_vrr = PAY*I_KINETIC_G4x_G4z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G4z_vrr;
    Double I_KINETIC_H4xz_G4z_vrr = PAZ*I_KINETIC_G4x_G4z_vrr+4*oned2z*I_KINETIC_G4x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G4z_vrr;
    Double I_KINETIC_H3x2y_G4z_vrr = PAY*I_KINETIC_G3xy_G4z_vrr+oned2z*I_KINETIC_F3x_G4z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G4z_vrr-bdz*I_TWOBODYOVERLAP_F3x_G4z_vrr;
    Double I_KINETIC_H3xyz_G4z_vrr = PAZ*I_KINETIC_G3xy_G4z_vrr+4*oned2z*I_KINETIC_G3xy_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G4z_vrr;
    Double I_KINETIC_H3x2z_G4z_vrr = PAZ*I_KINETIC_G3xz_G4z_vrr+oned2z*I_KINETIC_F3x_G4z_vrr+4*oned2z*I_KINETIC_G3xz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G4z_vrr-bdz*I_TWOBODYOVERLAP_F3x_G4z_vrr;
    Double I_KINETIC_H2x3y_G4z_vrr = PAX*I_KINETIC_Gx3y_G4z_vrr+oned2z*I_KINETIC_F3y_G4z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G4z_vrr-bdz*I_TWOBODYOVERLAP_F3y_G4z_vrr;
    Double I_KINETIC_H2x2yz_G4z_vrr = PAZ*I_KINETIC_G2x2y_G4z_vrr+4*oned2z*I_KINETIC_G2x2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G4z_vrr;
    Double I_KINETIC_H2xy2z_G4z_vrr = PAY*I_KINETIC_G2x2z_G4z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G4z_vrr;
    Double I_KINETIC_H2x3z_G4z_vrr = PAX*I_KINETIC_Gx3z_G4z_vrr+oned2z*I_KINETIC_F3z_G4z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G4z_vrr-bdz*I_TWOBODYOVERLAP_F3z_G4z_vrr;
    Double I_KINETIC_Hx4y_G4z_vrr = PAX*I_KINETIC_G4y_G4z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G4z_vrr;
    Double I_KINETIC_Hx3yz_G4z_vrr = PAZ*I_KINETIC_Gx3y_G4z_vrr+4*oned2z*I_KINETIC_Gx3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G4z_vrr;
    Double I_KINETIC_Hx2y2z_G4z_vrr = PAX*I_KINETIC_G2y2z_G4z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G4z_vrr;
    Double I_KINETIC_Hxy3z_G4z_vrr = PAY*I_KINETIC_Gx3z_G4z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G4z_vrr;
    Double I_KINETIC_Hx4z_G4z_vrr = PAX*I_KINETIC_G4z_G4z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G4z_vrr;
    Double I_KINETIC_H5y_G4z_vrr = PAY*I_KINETIC_G4y_G4z_vrr+4*oned2z*I_KINETIC_F3y_G4z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G4z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_G4z_vrr;
    Double I_KINETIC_H4yz_G4z_vrr = PAZ*I_KINETIC_G4y_G4z_vrr+4*oned2z*I_KINETIC_G4y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G4z_vrr;
    Double I_KINETIC_H3y2z_G4z_vrr = PAZ*I_KINETIC_G3yz_G4z_vrr+oned2z*I_KINETIC_F3y_G4z_vrr+4*oned2z*I_KINETIC_G3yz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G4z_vrr-bdz*I_TWOBODYOVERLAP_F3y_G4z_vrr;
    Double I_KINETIC_H2y3z_G4z_vrr = PAY*I_KINETIC_Gy3z_G4z_vrr+oned2z*I_KINETIC_F3z_G4z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G4z_vrr-bdz*I_TWOBODYOVERLAP_F3z_G4z_vrr;
    Double I_KINETIC_Hy4z_G4z_vrr = PAY*I_KINETIC_G4z_G4z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G4z_vrr;
    Double I_KINETIC_H5z_G4z_vrr = PAZ*I_KINETIC_G4z_G4z_vrr+4*oned2z*I_KINETIC_F3z_G4z_vrr+4*oned2z*I_KINETIC_G4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G4z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_G4z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_G_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_H_G_a_coefs = alpha;
    I_KINETIC_H5x_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G4x_vrr;
    I_KINETIC_H4xy_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G4x_vrr;
    I_KINETIC_H4xz_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G4x_vrr;
    I_KINETIC_H3x2y_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G4x_vrr;
    I_KINETIC_H3xyz_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G4x_vrr;
    I_KINETIC_H3x2z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G4x_vrr;
    I_KINETIC_H2x3y_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G4x_vrr;
    I_KINETIC_H2x2yz_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G4x_vrr;
    I_KINETIC_H2xy2z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G4x_vrr;
    I_KINETIC_H2x3z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G4x_vrr;
    I_KINETIC_Hx4y_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G4x_vrr;
    I_KINETIC_Hx3yz_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G4x_vrr;
    I_KINETIC_Hx2y2z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G4x_vrr;
    I_KINETIC_Hxy3z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G4x_vrr;
    I_KINETIC_Hx4z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G4x_vrr;
    I_KINETIC_H5y_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G4x_vrr;
    I_KINETIC_H4yz_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G4x_vrr;
    I_KINETIC_H3y2z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G4x_vrr;
    I_KINETIC_H2y3z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G4x_vrr;
    I_KINETIC_Hy4z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G4x_vrr;
    I_KINETIC_H5z_G4x_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G4x_vrr;
    I_KINETIC_H5x_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G3xy_vrr;
    I_KINETIC_H4xy_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G3xy_vrr;
    I_KINETIC_H4xz_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G3xy_vrr;
    I_KINETIC_H3x2y_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G3xy_vrr;
    I_KINETIC_H3xyz_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G3xy_vrr;
    I_KINETIC_H3x2z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G3xy_vrr;
    I_KINETIC_H2x3y_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G3xy_vrr;
    I_KINETIC_H2x2yz_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G3xy_vrr;
    I_KINETIC_H2xy2z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G3xy_vrr;
    I_KINETIC_H2x3z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G3xy_vrr;
    I_KINETIC_Hx4y_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G3xy_vrr;
    I_KINETIC_Hx3yz_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G3xy_vrr;
    I_KINETIC_Hx2y2z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G3xy_vrr;
    I_KINETIC_Hxy3z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G3xy_vrr;
    I_KINETIC_Hx4z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G3xy_vrr;
    I_KINETIC_H5y_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G3xy_vrr;
    I_KINETIC_H4yz_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G3xy_vrr;
    I_KINETIC_H3y2z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G3xy_vrr;
    I_KINETIC_H2y3z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G3xy_vrr;
    I_KINETIC_Hy4z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G3xy_vrr;
    I_KINETIC_H5z_G3xy_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G3xy_vrr;
    I_KINETIC_H5x_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G3xz_vrr;
    I_KINETIC_H4xy_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G3xz_vrr;
    I_KINETIC_H4xz_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G3xz_vrr;
    I_KINETIC_H3x2y_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G3xz_vrr;
    I_KINETIC_H3xyz_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G3xz_vrr;
    I_KINETIC_H3x2z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G3xz_vrr;
    I_KINETIC_H2x3y_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G3xz_vrr;
    I_KINETIC_H2x2yz_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G3xz_vrr;
    I_KINETIC_H2xy2z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G3xz_vrr;
    I_KINETIC_H2x3z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G3xz_vrr;
    I_KINETIC_Hx4y_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G3xz_vrr;
    I_KINETIC_Hx3yz_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G3xz_vrr;
    I_KINETIC_Hx2y2z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G3xz_vrr;
    I_KINETIC_Hxy3z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G3xz_vrr;
    I_KINETIC_Hx4z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G3xz_vrr;
    I_KINETIC_H5y_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G3xz_vrr;
    I_KINETIC_H4yz_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G3xz_vrr;
    I_KINETIC_H3y2z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G3xz_vrr;
    I_KINETIC_H2y3z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G3xz_vrr;
    I_KINETIC_Hy4z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G3xz_vrr;
    I_KINETIC_H5z_G3xz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G3xz_vrr;
    I_KINETIC_H5x_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G2x2y_vrr;
    I_KINETIC_H4xy_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G2x2y_vrr;
    I_KINETIC_H4xz_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G2x2y_vrr;
    I_KINETIC_H3x2y_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G2x2y_vrr;
    I_KINETIC_H3xyz_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G2x2y_vrr;
    I_KINETIC_H3x2z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G2x2y_vrr;
    I_KINETIC_H2x3y_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G2x2y_vrr;
    I_KINETIC_H2x2yz_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G2x2y_vrr;
    I_KINETIC_H2xy2z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G2x2y_vrr;
    I_KINETIC_H2x3z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G2x2y_vrr;
    I_KINETIC_Hx4y_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G2x2y_vrr;
    I_KINETIC_Hx3yz_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G2x2y_vrr;
    I_KINETIC_Hx2y2z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G2x2y_vrr;
    I_KINETIC_Hxy3z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G2x2y_vrr;
    I_KINETIC_Hx4z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G2x2y_vrr;
    I_KINETIC_H5y_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G2x2y_vrr;
    I_KINETIC_H4yz_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G2x2y_vrr;
    I_KINETIC_H3y2z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G2x2y_vrr;
    I_KINETIC_H2y3z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G2x2y_vrr;
    I_KINETIC_Hy4z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G2x2y_vrr;
    I_KINETIC_H5z_G2x2y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G2x2y_vrr;
    I_KINETIC_H5x_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G2xyz_vrr;
    I_KINETIC_H4xy_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G2xyz_vrr;
    I_KINETIC_H4xz_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G2xyz_vrr;
    I_KINETIC_H3x2y_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G2xyz_vrr;
    I_KINETIC_H3xyz_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G2xyz_vrr;
    I_KINETIC_H3x2z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G2xyz_vrr;
    I_KINETIC_H2x3y_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G2xyz_vrr;
    I_KINETIC_H2x2yz_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G2xyz_vrr;
    I_KINETIC_H2xy2z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G2xyz_vrr;
    I_KINETIC_H2x3z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G2xyz_vrr;
    I_KINETIC_Hx4y_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G2xyz_vrr;
    I_KINETIC_Hx3yz_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G2xyz_vrr;
    I_KINETIC_Hx2y2z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G2xyz_vrr;
    I_KINETIC_Hxy3z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G2xyz_vrr;
    I_KINETIC_Hx4z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G2xyz_vrr;
    I_KINETIC_H5y_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G2xyz_vrr;
    I_KINETIC_H4yz_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G2xyz_vrr;
    I_KINETIC_H3y2z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G2xyz_vrr;
    I_KINETIC_H2y3z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G2xyz_vrr;
    I_KINETIC_Hy4z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G2xyz_vrr;
    I_KINETIC_H5z_G2xyz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G2xyz_vrr;
    I_KINETIC_H5x_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G2x2z_vrr;
    I_KINETIC_H4xy_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G2x2z_vrr;
    I_KINETIC_H4xz_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G2x2z_vrr;
    I_KINETIC_H3x2y_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G2x2z_vrr;
    I_KINETIC_H3xyz_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G2x2z_vrr;
    I_KINETIC_H3x2z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G2x2z_vrr;
    I_KINETIC_H2x3y_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G2x2z_vrr;
    I_KINETIC_H2x2yz_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G2x2z_vrr;
    I_KINETIC_H2xy2z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G2x2z_vrr;
    I_KINETIC_H2x3z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G2x2z_vrr;
    I_KINETIC_Hx4y_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G2x2z_vrr;
    I_KINETIC_Hx3yz_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G2x2z_vrr;
    I_KINETIC_Hx2y2z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G2x2z_vrr;
    I_KINETIC_Hxy3z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G2x2z_vrr;
    I_KINETIC_Hx4z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G2x2z_vrr;
    I_KINETIC_H5y_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G2x2z_vrr;
    I_KINETIC_H4yz_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G2x2z_vrr;
    I_KINETIC_H3y2z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G2x2z_vrr;
    I_KINETIC_H2y3z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G2x2z_vrr;
    I_KINETIC_Hy4z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G2x2z_vrr;
    I_KINETIC_H5z_G2x2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G2x2z_vrr;
    I_KINETIC_H5x_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_Gx3y_vrr;
    I_KINETIC_H4xy_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_Gx3y_vrr;
    I_KINETIC_H4xz_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_Gx3y_vrr;
    I_KINETIC_H3x2y_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_Gx3y_vrr;
    I_KINETIC_H3xyz_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_Gx3y_vrr;
    I_KINETIC_H3x2z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_Gx3y_vrr;
    I_KINETIC_H2x3y_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_Gx3y_vrr;
    I_KINETIC_H2x2yz_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_Gx3y_vrr;
    I_KINETIC_H2xy2z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_Gx3y_vrr;
    I_KINETIC_H2x3z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_Gx3y_vrr;
    I_KINETIC_Hx4y_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_Gx3y_vrr;
    I_KINETIC_Hx3yz_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_Gx3y_vrr;
    I_KINETIC_Hx2y2z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_Gx3y_vrr;
    I_KINETIC_Hxy3z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_Gx3y_vrr;
    I_KINETIC_Hx4z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_Gx3y_vrr;
    I_KINETIC_H5y_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_Gx3y_vrr;
    I_KINETIC_H4yz_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_Gx3y_vrr;
    I_KINETIC_H3y2z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_Gx3y_vrr;
    I_KINETIC_H2y3z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_Gx3y_vrr;
    I_KINETIC_Hy4z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_Gx3y_vrr;
    I_KINETIC_H5z_Gx3y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_Gx3y_vrr;
    I_KINETIC_H5x_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_Gx2yz_vrr;
    I_KINETIC_H4xy_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_Gx2yz_vrr;
    I_KINETIC_H4xz_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_Gx2yz_vrr;
    I_KINETIC_H3x2y_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_Gx2yz_vrr;
    I_KINETIC_H3xyz_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_Gx2yz_vrr;
    I_KINETIC_H3x2z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_Gx2yz_vrr;
    I_KINETIC_H2x3y_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_Gx2yz_vrr;
    I_KINETIC_H2x2yz_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_Gx2yz_vrr;
    I_KINETIC_H2xy2z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_Gx2yz_vrr;
    I_KINETIC_H2x3z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_Gx2yz_vrr;
    I_KINETIC_Hx4y_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_Gx2yz_vrr;
    I_KINETIC_Hx3yz_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_Gx2yz_vrr;
    I_KINETIC_Hx2y2z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_Gx2yz_vrr;
    I_KINETIC_Hxy3z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_Gx2yz_vrr;
    I_KINETIC_Hx4z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_Gx2yz_vrr;
    I_KINETIC_H5y_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_Gx2yz_vrr;
    I_KINETIC_H4yz_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_Gx2yz_vrr;
    I_KINETIC_H3y2z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_Gx2yz_vrr;
    I_KINETIC_H2y3z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_Gx2yz_vrr;
    I_KINETIC_Hy4z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_Gx2yz_vrr;
    I_KINETIC_H5z_Gx2yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_Gx2yz_vrr;
    I_KINETIC_H5x_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_Gxy2z_vrr;
    I_KINETIC_H4xy_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_Gxy2z_vrr;
    I_KINETIC_H4xz_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_Gxy2z_vrr;
    I_KINETIC_H3x2y_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_Gxy2z_vrr;
    I_KINETIC_H3xyz_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_Gxy2z_vrr;
    I_KINETIC_H3x2z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_Gxy2z_vrr;
    I_KINETIC_H2x3y_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_Gxy2z_vrr;
    I_KINETIC_H2x2yz_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_Gxy2z_vrr;
    I_KINETIC_H2xy2z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_Gxy2z_vrr;
    I_KINETIC_H2x3z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_Gxy2z_vrr;
    I_KINETIC_Hx4y_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_Gxy2z_vrr;
    I_KINETIC_Hx3yz_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_Gxy2z_vrr;
    I_KINETIC_Hx2y2z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_Gxy2z_vrr;
    I_KINETIC_Hxy3z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_Gxy2z_vrr;
    I_KINETIC_Hx4z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_Gxy2z_vrr;
    I_KINETIC_H5y_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_Gxy2z_vrr;
    I_KINETIC_H4yz_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_Gxy2z_vrr;
    I_KINETIC_H3y2z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_Gxy2z_vrr;
    I_KINETIC_H2y3z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_Gxy2z_vrr;
    I_KINETIC_Hy4z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_Gxy2z_vrr;
    I_KINETIC_H5z_Gxy2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_Gxy2z_vrr;
    I_KINETIC_H5x_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_Gx3z_vrr;
    I_KINETIC_H4xy_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_Gx3z_vrr;
    I_KINETIC_H4xz_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_Gx3z_vrr;
    I_KINETIC_H3x2y_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_Gx3z_vrr;
    I_KINETIC_H3xyz_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_Gx3z_vrr;
    I_KINETIC_H3x2z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_Gx3z_vrr;
    I_KINETIC_H2x3y_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_Gx3z_vrr;
    I_KINETIC_H2x2yz_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_Gx3z_vrr;
    I_KINETIC_H2xy2z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_Gx3z_vrr;
    I_KINETIC_H2x3z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_Gx3z_vrr;
    I_KINETIC_Hx4y_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_Gx3z_vrr;
    I_KINETIC_Hx3yz_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_Gx3z_vrr;
    I_KINETIC_Hx2y2z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_Gx3z_vrr;
    I_KINETIC_Hxy3z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_Gx3z_vrr;
    I_KINETIC_Hx4z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_Gx3z_vrr;
    I_KINETIC_H5y_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_Gx3z_vrr;
    I_KINETIC_H4yz_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_Gx3z_vrr;
    I_KINETIC_H3y2z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_Gx3z_vrr;
    I_KINETIC_H2y3z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_Gx3z_vrr;
    I_KINETIC_Hy4z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_Gx3z_vrr;
    I_KINETIC_H5z_Gx3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_Gx3z_vrr;
    I_KINETIC_H5x_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G4y_vrr;
    I_KINETIC_H4xy_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G4y_vrr;
    I_KINETIC_H4xz_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G4y_vrr;
    I_KINETIC_H3x2y_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G4y_vrr;
    I_KINETIC_H3xyz_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G4y_vrr;
    I_KINETIC_H3x2z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G4y_vrr;
    I_KINETIC_H2x3y_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G4y_vrr;
    I_KINETIC_H2x2yz_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G4y_vrr;
    I_KINETIC_H2xy2z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G4y_vrr;
    I_KINETIC_H2x3z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G4y_vrr;
    I_KINETIC_Hx4y_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G4y_vrr;
    I_KINETIC_Hx3yz_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G4y_vrr;
    I_KINETIC_Hx2y2z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G4y_vrr;
    I_KINETIC_Hxy3z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G4y_vrr;
    I_KINETIC_Hx4z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G4y_vrr;
    I_KINETIC_H5y_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G4y_vrr;
    I_KINETIC_H4yz_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G4y_vrr;
    I_KINETIC_H3y2z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G4y_vrr;
    I_KINETIC_H2y3z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G4y_vrr;
    I_KINETIC_Hy4z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G4y_vrr;
    I_KINETIC_H5z_G4y_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G4y_vrr;
    I_KINETIC_H5x_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G3yz_vrr;
    I_KINETIC_H4xy_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G3yz_vrr;
    I_KINETIC_H4xz_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G3yz_vrr;
    I_KINETIC_H3x2y_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G3yz_vrr;
    I_KINETIC_H3xyz_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G3yz_vrr;
    I_KINETIC_H3x2z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G3yz_vrr;
    I_KINETIC_H2x3y_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G3yz_vrr;
    I_KINETIC_H2x2yz_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G3yz_vrr;
    I_KINETIC_H2xy2z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G3yz_vrr;
    I_KINETIC_H2x3z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G3yz_vrr;
    I_KINETIC_Hx4y_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G3yz_vrr;
    I_KINETIC_Hx3yz_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G3yz_vrr;
    I_KINETIC_Hx2y2z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G3yz_vrr;
    I_KINETIC_Hxy3z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G3yz_vrr;
    I_KINETIC_Hx4z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G3yz_vrr;
    I_KINETIC_H5y_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G3yz_vrr;
    I_KINETIC_H4yz_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G3yz_vrr;
    I_KINETIC_H3y2z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G3yz_vrr;
    I_KINETIC_H2y3z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G3yz_vrr;
    I_KINETIC_Hy4z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G3yz_vrr;
    I_KINETIC_H5z_G3yz_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G3yz_vrr;
    I_KINETIC_H5x_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G2y2z_vrr;
    I_KINETIC_H4xy_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G2y2z_vrr;
    I_KINETIC_H4xz_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G2y2z_vrr;
    I_KINETIC_H3x2y_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G2y2z_vrr;
    I_KINETIC_H3xyz_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G2y2z_vrr;
    I_KINETIC_H3x2z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G2y2z_vrr;
    I_KINETIC_H2x3y_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G2y2z_vrr;
    I_KINETIC_H2x2yz_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G2y2z_vrr;
    I_KINETIC_H2xy2z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G2y2z_vrr;
    I_KINETIC_H2x3z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G2y2z_vrr;
    I_KINETIC_Hx4y_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G2y2z_vrr;
    I_KINETIC_Hx3yz_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G2y2z_vrr;
    I_KINETIC_Hx2y2z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G2y2z_vrr;
    I_KINETIC_Hxy3z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G2y2z_vrr;
    I_KINETIC_Hx4z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G2y2z_vrr;
    I_KINETIC_H5y_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G2y2z_vrr;
    I_KINETIC_H4yz_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G2y2z_vrr;
    I_KINETIC_H3y2z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G2y2z_vrr;
    I_KINETIC_H2y3z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G2y2z_vrr;
    I_KINETIC_Hy4z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G2y2z_vrr;
    I_KINETIC_H5z_G2y2z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G2y2z_vrr;
    I_KINETIC_H5x_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_Gy3z_vrr;
    I_KINETIC_H4xy_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_Gy3z_vrr;
    I_KINETIC_H4xz_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_Gy3z_vrr;
    I_KINETIC_H3x2y_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_Gy3z_vrr;
    I_KINETIC_H3xyz_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_Gy3z_vrr;
    I_KINETIC_H3x2z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_Gy3z_vrr;
    I_KINETIC_H2x3y_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_Gy3z_vrr;
    I_KINETIC_H2x2yz_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_Gy3z_vrr;
    I_KINETIC_H2xy2z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_Gy3z_vrr;
    I_KINETIC_H2x3z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_Gy3z_vrr;
    I_KINETIC_Hx4y_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_Gy3z_vrr;
    I_KINETIC_Hx3yz_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_Gy3z_vrr;
    I_KINETIC_Hx2y2z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_Gy3z_vrr;
    I_KINETIC_Hxy3z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_Gy3z_vrr;
    I_KINETIC_Hx4z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_Gy3z_vrr;
    I_KINETIC_H5y_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_Gy3z_vrr;
    I_KINETIC_H4yz_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_Gy3z_vrr;
    I_KINETIC_H3y2z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_Gy3z_vrr;
    I_KINETIC_H2y3z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_Gy3z_vrr;
    I_KINETIC_Hy4z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_Gy3z_vrr;
    I_KINETIC_H5z_Gy3z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_Gy3z_vrr;
    I_KINETIC_H5x_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5x_G4z_vrr;
    I_KINETIC_H4xy_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xy_G4z_vrr;
    I_KINETIC_H4xz_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4xz_G4z_vrr;
    I_KINETIC_H3x2y_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2y_G4z_vrr;
    I_KINETIC_H3xyz_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3xyz_G4z_vrr;
    I_KINETIC_H3x2z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3x2z_G4z_vrr;
    I_KINETIC_H2x3y_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3y_G4z_vrr;
    I_KINETIC_H2x2yz_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x2yz_G4z_vrr;
    I_KINETIC_H2xy2z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2xy2z_G4z_vrr;
    I_KINETIC_H2x3z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2x3z_G4z_vrr;
    I_KINETIC_Hx4y_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4y_G4z_vrr;
    I_KINETIC_Hx3yz_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx3yz_G4z_vrr;
    I_KINETIC_Hx2y2z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx2y2z_G4z_vrr;
    I_KINETIC_Hxy3z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hxy3z_G4z_vrr;
    I_KINETIC_Hx4z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hx4z_G4z_vrr;
    I_KINETIC_H5y_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5y_G4z_vrr;
    I_KINETIC_H4yz_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H4yz_G4z_vrr;
    I_KINETIC_H3y2z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H3y2z_G4z_vrr;
    I_KINETIC_H2y3z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H2y3z_G4z_vrr;
    I_KINETIC_Hy4z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_Hy4z_G4z_vrr;
    I_KINETIC_H5z_G4z_a += SQ_KINETIC_H_G_a_coefs*I_KINETIC_H5z_G4z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_G
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_F3x_G4x += I_KINETIC_F3x_G4x_vrr;
    I_KINETIC_F2xy_G4x += I_KINETIC_F2xy_G4x_vrr;
    I_KINETIC_F2xz_G4x += I_KINETIC_F2xz_G4x_vrr;
    I_KINETIC_Fx2y_G4x += I_KINETIC_Fx2y_G4x_vrr;
    I_KINETIC_Fxyz_G4x += I_KINETIC_Fxyz_G4x_vrr;
    I_KINETIC_Fx2z_G4x += I_KINETIC_Fx2z_G4x_vrr;
    I_KINETIC_F3y_G4x += I_KINETIC_F3y_G4x_vrr;
    I_KINETIC_F2yz_G4x += I_KINETIC_F2yz_G4x_vrr;
    I_KINETIC_Fy2z_G4x += I_KINETIC_Fy2z_G4x_vrr;
    I_KINETIC_F3z_G4x += I_KINETIC_F3z_G4x_vrr;
    I_KINETIC_F3x_G3xy += I_KINETIC_F3x_G3xy_vrr;
    I_KINETIC_F2xy_G3xy += I_KINETIC_F2xy_G3xy_vrr;
    I_KINETIC_F2xz_G3xy += I_KINETIC_F2xz_G3xy_vrr;
    I_KINETIC_Fx2y_G3xy += I_KINETIC_Fx2y_G3xy_vrr;
    I_KINETIC_Fxyz_G3xy += I_KINETIC_Fxyz_G3xy_vrr;
    I_KINETIC_Fx2z_G3xy += I_KINETIC_Fx2z_G3xy_vrr;
    I_KINETIC_F3y_G3xy += I_KINETIC_F3y_G3xy_vrr;
    I_KINETIC_F2yz_G3xy += I_KINETIC_F2yz_G3xy_vrr;
    I_KINETIC_Fy2z_G3xy += I_KINETIC_Fy2z_G3xy_vrr;
    I_KINETIC_F3z_G3xy += I_KINETIC_F3z_G3xy_vrr;
    I_KINETIC_F3x_G3xz += I_KINETIC_F3x_G3xz_vrr;
    I_KINETIC_F2xy_G3xz += I_KINETIC_F2xy_G3xz_vrr;
    I_KINETIC_F2xz_G3xz += I_KINETIC_F2xz_G3xz_vrr;
    I_KINETIC_Fx2y_G3xz += I_KINETIC_Fx2y_G3xz_vrr;
    I_KINETIC_Fxyz_G3xz += I_KINETIC_Fxyz_G3xz_vrr;
    I_KINETIC_Fx2z_G3xz += I_KINETIC_Fx2z_G3xz_vrr;
    I_KINETIC_F3y_G3xz += I_KINETIC_F3y_G3xz_vrr;
    I_KINETIC_F2yz_G3xz += I_KINETIC_F2yz_G3xz_vrr;
    I_KINETIC_Fy2z_G3xz += I_KINETIC_Fy2z_G3xz_vrr;
    I_KINETIC_F3z_G3xz += I_KINETIC_F3z_G3xz_vrr;
    I_KINETIC_F3x_G2x2y += I_KINETIC_F3x_G2x2y_vrr;
    I_KINETIC_F2xy_G2x2y += I_KINETIC_F2xy_G2x2y_vrr;
    I_KINETIC_F2xz_G2x2y += I_KINETIC_F2xz_G2x2y_vrr;
    I_KINETIC_Fx2y_G2x2y += I_KINETIC_Fx2y_G2x2y_vrr;
    I_KINETIC_Fxyz_G2x2y += I_KINETIC_Fxyz_G2x2y_vrr;
    I_KINETIC_Fx2z_G2x2y += I_KINETIC_Fx2z_G2x2y_vrr;
    I_KINETIC_F3y_G2x2y += I_KINETIC_F3y_G2x2y_vrr;
    I_KINETIC_F2yz_G2x2y += I_KINETIC_F2yz_G2x2y_vrr;
    I_KINETIC_Fy2z_G2x2y += I_KINETIC_Fy2z_G2x2y_vrr;
    I_KINETIC_F3z_G2x2y += I_KINETIC_F3z_G2x2y_vrr;
    I_KINETIC_F3x_G2xyz += I_KINETIC_F3x_G2xyz_vrr;
    I_KINETIC_F2xy_G2xyz += I_KINETIC_F2xy_G2xyz_vrr;
    I_KINETIC_F2xz_G2xyz += I_KINETIC_F2xz_G2xyz_vrr;
    I_KINETIC_Fx2y_G2xyz += I_KINETIC_Fx2y_G2xyz_vrr;
    I_KINETIC_Fxyz_G2xyz += I_KINETIC_Fxyz_G2xyz_vrr;
    I_KINETIC_Fx2z_G2xyz += I_KINETIC_Fx2z_G2xyz_vrr;
    I_KINETIC_F3y_G2xyz += I_KINETIC_F3y_G2xyz_vrr;
    I_KINETIC_F2yz_G2xyz += I_KINETIC_F2yz_G2xyz_vrr;
    I_KINETIC_Fy2z_G2xyz += I_KINETIC_Fy2z_G2xyz_vrr;
    I_KINETIC_F3z_G2xyz += I_KINETIC_F3z_G2xyz_vrr;
    I_KINETIC_F3x_G2x2z += I_KINETIC_F3x_G2x2z_vrr;
    I_KINETIC_F2xy_G2x2z += I_KINETIC_F2xy_G2x2z_vrr;
    I_KINETIC_F2xz_G2x2z += I_KINETIC_F2xz_G2x2z_vrr;
    I_KINETIC_Fx2y_G2x2z += I_KINETIC_Fx2y_G2x2z_vrr;
    I_KINETIC_Fxyz_G2x2z += I_KINETIC_Fxyz_G2x2z_vrr;
    I_KINETIC_Fx2z_G2x2z += I_KINETIC_Fx2z_G2x2z_vrr;
    I_KINETIC_F3y_G2x2z += I_KINETIC_F3y_G2x2z_vrr;
    I_KINETIC_F2yz_G2x2z += I_KINETIC_F2yz_G2x2z_vrr;
    I_KINETIC_Fy2z_G2x2z += I_KINETIC_Fy2z_G2x2z_vrr;
    I_KINETIC_F3z_G2x2z += I_KINETIC_F3z_G2x2z_vrr;
    I_KINETIC_F3x_Gx3y += I_KINETIC_F3x_Gx3y_vrr;
    I_KINETIC_F2xy_Gx3y += I_KINETIC_F2xy_Gx3y_vrr;
    I_KINETIC_F2xz_Gx3y += I_KINETIC_F2xz_Gx3y_vrr;
    I_KINETIC_Fx2y_Gx3y += I_KINETIC_Fx2y_Gx3y_vrr;
    I_KINETIC_Fxyz_Gx3y += I_KINETIC_Fxyz_Gx3y_vrr;
    I_KINETIC_Fx2z_Gx3y += I_KINETIC_Fx2z_Gx3y_vrr;
    I_KINETIC_F3y_Gx3y += I_KINETIC_F3y_Gx3y_vrr;
    I_KINETIC_F2yz_Gx3y += I_KINETIC_F2yz_Gx3y_vrr;
    I_KINETIC_Fy2z_Gx3y += I_KINETIC_Fy2z_Gx3y_vrr;
    I_KINETIC_F3z_Gx3y += I_KINETIC_F3z_Gx3y_vrr;
    I_KINETIC_F3x_Gx2yz += I_KINETIC_F3x_Gx2yz_vrr;
    I_KINETIC_F2xy_Gx2yz += I_KINETIC_F2xy_Gx2yz_vrr;
    I_KINETIC_F2xz_Gx2yz += I_KINETIC_F2xz_Gx2yz_vrr;
    I_KINETIC_Fx2y_Gx2yz += I_KINETIC_Fx2y_Gx2yz_vrr;
    I_KINETIC_Fxyz_Gx2yz += I_KINETIC_Fxyz_Gx2yz_vrr;
    I_KINETIC_Fx2z_Gx2yz += I_KINETIC_Fx2z_Gx2yz_vrr;
    I_KINETIC_F3y_Gx2yz += I_KINETIC_F3y_Gx2yz_vrr;
    I_KINETIC_F2yz_Gx2yz += I_KINETIC_F2yz_Gx2yz_vrr;
    I_KINETIC_Fy2z_Gx2yz += I_KINETIC_Fy2z_Gx2yz_vrr;
    I_KINETIC_F3z_Gx2yz += I_KINETIC_F3z_Gx2yz_vrr;
    I_KINETIC_F3x_Gxy2z += I_KINETIC_F3x_Gxy2z_vrr;
    I_KINETIC_F2xy_Gxy2z += I_KINETIC_F2xy_Gxy2z_vrr;
    I_KINETIC_F2xz_Gxy2z += I_KINETIC_F2xz_Gxy2z_vrr;
    I_KINETIC_Fx2y_Gxy2z += I_KINETIC_Fx2y_Gxy2z_vrr;
    I_KINETIC_Fxyz_Gxy2z += I_KINETIC_Fxyz_Gxy2z_vrr;
    I_KINETIC_Fx2z_Gxy2z += I_KINETIC_Fx2z_Gxy2z_vrr;
    I_KINETIC_F3y_Gxy2z += I_KINETIC_F3y_Gxy2z_vrr;
    I_KINETIC_F2yz_Gxy2z += I_KINETIC_F2yz_Gxy2z_vrr;
    I_KINETIC_Fy2z_Gxy2z += I_KINETIC_Fy2z_Gxy2z_vrr;
    I_KINETIC_F3z_Gxy2z += I_KINETIC_F3z_Gxy2z_vrr;
    I_KINETIC_F3x_Gx3z += I_KINETIC_F3x_Gx3z_vrr;
    I_KINETIC_F2xy_Gx3z += I_KINETIC_F2xy_Gx3z_vrr;
    I_KINETIC_F2xz_Gx3z += I_KINETIC_F2xz_Gx3z_vrr;
    I_KINETIC_Fx2y_Gx3z += I_KINETIC_Fx2y_Gx3z_vrr;
    I_KINETIC_Fxyz_Gx3z += I_KINETIC_Fxyz_Gx3z_vrr;
    I_KINETIC_Fx2z_Gx3z += I_KINETIC_Fx2z_Gx3z_vrr;
    I_KINETIC_F3y_Gx3z += I_KINETIC_F3y_Gx3z_vrr;
    I_KINETIC_F2yz_Gx3z += I_KINETIC_F2yz_Gx3z_vrr;
    I_KINETIC_Fy2z_Gx3z += I_KINETIC_Fy2z_Gx3z_vrr;
    I_KINETIC_F3z_Gx3z += I_KINETIC_F3z_Gx3z_vrr;
    I_KINETIC_F3x_G4y += I_KINETIC_F3x_G4y_vrr;
    I_KINETIC_F2xy_G4y += I_KINETIC_F2xy_G4y_vrr;
    I_KINETIC_F2xz_G4y += I_KINETIC_F2xz_G4y_vrr;
    I_KINETIC_Fx2y_G4y += I_KINETIC_Fx2y_G4y_vrr;
    I_KINETIC_Fxyz_G4y += I_KINETIC_Fxyz_G4y_vrr;
    I_KINETIC_Fx2z_G4y += I_KINETIC_Fx2z_G4y_vrr;
    I_KINETIC_F3y_G4y += I_KINETIC_F3y_G4y_vrr;
    I_KINETIC_F2yz_G4y += I_KINETIC_F2yz_G4y_vrr;
    I_KINETIC_Fy2z_G4y += I_KINETIC_Fy2z_G4y_vrr;
    I_KINETIC_F3z_G4y += I_KINETIC_F3z_G4y_vrr;
    I_KINETIC_F3x_G3yz += I_KINETIC_F3x_G3yz_vrr;
    I_KINETIC_F2xy_G3yz += I_KINETIC_F2xy_G3yz_vrr;
    I_KINETIC_F2xz_G3yz += I_KINETIC_F2xz_G3yz_vrr;
    I_KINETIC_Fx2y_G3yz += I_KINETIC_Fx2y_G3yz_vrr;
    I_KINETIC_Fxyz_G3yz += I_KINETIC_Fxyz_G3yz_vrr;
    I_KINETIC_Fx2z_G3yz += I_KINETIC_Fx2z_G3yz_vrr;
    I_KINETIC_F3y_G3yz += I_KINETIC_F3y_G3yz_vrr;
    I_KINETIC_F2yz_G3yz += I_KINETIC_F2yz_G3yz_vrr;
    I_KINETIC_Fy2z_G3yz += I_KINETIC_Fy2z_G3yz_vrr;
    I_KINETIC_F3z_G3yz += I_KINETIC_F3z_G3yz_vrr;
    I_KINETIC_F3x_G2y2z += I_KINETIC_F3x_G2y2z_vrr;
    I_KINETIC_F2xy_G2y2z += I_KINETIC_F2xy_G2y2z_vrr;
    I_KINETIC_F2xz_G2y2z += I_KINETIC_F2xz_G2y2z_vrr;
    I_KINETIC_Fx2y_G2y2z += I_KINETIC_Fx2y_G2y2z_vrr;
    I_KINETIC_Fxyz_G2y2z += I_KINETIC_Fxyz_G2y2z_vrr;
    I_KINETIC_Fx2z_G2y2z += I_KINETIC_Fx2z_G2y2z_vrr;
    I_KINETIC_F3y_G2y2z += I_KINETIC_F3y_G2y2z_vrr;
    I_KINETIC_F2yz_G2y2z += I_KINETIC_F2yz_G2y2z_vrr;
    I_KINETIC_Fy2z_G2y2z += I_KINETIC_Fy2z_G2y2z_vrr;
    I_KINETIC_F3z_G2y2z += I_KINETIC_F3z_G2y2z_vrr;
    I_KINETIC_F3x_Gy3z += I_KINETIC_F3x_Gy3z_vrr;
    I_KINETIC_F2xy_Gy3z += I_KINETIC_F2xy_Gy3z_vrr;
    I_KINETIC_F2xz_Gy3z += I_KINETIC_F2xz_Gy3z_vrr;
    I_KINETIC_Fx2y_Gy3z += I_KINETIC_Fx2y_Gy3z_vrr;
    I_KINETIC_Fxyz_Gy3z += I_KINETIC_Fxyz_Gy3z_vrr;
    I_KINETIC_Fx2z_Gy3z += I_KINETIC_Fx2z_Gy3z_vrr;
    I_KINETIC_F3y_Gy3z += I_KINETIC_F3y_Gy3z_vrr;
    I_KINETIC_F2yz_Gy3z += I_KINETIC_F2yz_Gy3z_vrr;
    I_KINETIC_Fy2z_Gy3z += I_KINETIC_Fy2z_Gy3z_vrr;
    I_KINETIC_F3z_Gy3z += I_KINETIC_F3z_Gy3z_vrr;
    I_KINETIC_F3x_G4z += I_KINETIC_F3x_G4z_vrr;
    I_KINETIC_F2xy_G4z += I_KINETIC_F2xy_G4z_vrr;
    I_KINETIC_F2xz_G4z += I_KINETIC_F2xz_G4z_vrr;
    I_KINETIC_Fx2y_G4z += I_KINETIC_Fx2y_G4z_vrr;
    I_KINETIC_Fxyz_G4z += I_KINETIC_Fxyz_G4z_vrr;
    I_KINETIC_Fx2z_G4z += I_KINETIC_Fx2z_G4z_vrr;
    I_KINETIC_F3y_G4z += I_KINETIC_F3y_G4z_vrr;
    I_KINETIC_F2yz_G4z += I_KINETIC_F2yz_G4z_vrr;
    I_KINETIC_Fy2z_G4z += I_KINETIC_Fy2z_G4z_vrr;
    I_KINETIC_F3z_G4z += I_KINETIC_F3z_G4z_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_G_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_G_a
   * RHS shell quartet name: SQ_KINETIC_F_G
   ************************************************************/
  abcd[0] = 2.0E0*I_KINETIC_H5x_G4x_a-4*I_KINETIC_F3x_G4x;
  abcd[1] = 2.0E0*I_KINETIC_H4xy_G4x_a-3*I_KINETIC_F2xy_G4x;
  abcd[2] = 2.0E0*I_KINETIC_H4xz_G4x_a-3*I_KINETIC_F2xz_G4x;
  abcd[3] = 2.0E0*I_KINETIC_H3x2y_G4x_a-2*I_KINETIC_Fx2y_G4x;
  abcd[4] = 2.0E0*I_KINETIC_H3xyz_G4x_a-2*I_KINETIC_Fxyz_G4x;
  abcd[5] = 2.0E0*I_KINETIC_H3x2z_G4x_a-2*I_KINETIC_Fx2z_G4x;
  abcd[6] = 2.0E0*I_KINETIC_H2x3y_G4x_a-1*I_KINETIC_F3y_G4x;
  abcd[7] = 2.0E0*I_KINETIC_H2x2yz_G4x_a-1*I_KINETIC_F2yz_G4x;
  abcd[8] = 2.0E0*I_KINETIC_H2xy2z_G4x_a-1*I_KINETIC_Fy2z_G4x;
  abcd[9] = 2.0E0*I_KINETIC_H2x3z_G4x_a-1*I_KINETIC_F3z_G4x;
  abcd[10] = 2.0E0*I_KINETIC_Hx4y_G4x_a;
  abcd[11] = 2.0E0*I_KINETIC_Hx3yz_G4x_a;
  abcd[12] = 2.0E0*I_KINETIC_Hx2y2z_G4x_a;
  abcd[13] = 2.0E0*I_KINETIC_Hxy3z_G4x_a;
  abcd[14] = 2.0E0*I_KINETIC_Hx4z_G4x_a;
  abcd[15] = 2.0E0*I_KINETIC_H5x_G3xy_a-4*I_KINETIC_F3x_G3xy;
  abcd[16] = 2.0E0*I_KINETIC_H4xy_G3xy_a-3*I_KINETIC_F2xy_G3xy;
  abcd[17] = 2.0E0*I_KINETIC_H4xz_G3xy_a-3*I_KINETIC_F2xz_G3xy;
  abcd[18] = 2.0E0*I_KINETIC_H3x2y_G3xy_a-2*I_KINETIC_Fx2y_G3xy;
  abcd[19] = 2.0E0*I_KINETIC_H3xyz_G3xy_a-2*I_KINETIC_Fxyz_G3xy;
  abcd[20] = 2.0E0*I_KINETIC_H3x2z_G3xy_a-2*I_KINETIC_Fx2z_G3xy;
  abcd[21] = 2.0E0*I_KINETIC_H2x3y_G3xy_a-1*I_KINETIC_F3y_G3xy;
  abcd[22] = 2.0E0*I_KINETIC_H2x2yz_G3xy_a-1*I_KINETIC_F2yz_G3xy;
  abcd[23] = 2.0E0*I_KINETIC_H2xy2z_G3xy_a-1*I_KINETIC_Fy2z_G3xy;
  abcd[24] = 2.0E0*I_KINETIC_H2x3z_G3xy_a-1*I_KINETIC_F3z_G3xy;
  abcd[25] = 2.0E0*I_KINETIC_Hx4y_G3xy_a;
  abcd[26] = 2.0E0*I_KINETIC_Hx3yz_G3xy_a;
  abcd[27] = 2.0E0*I_KINETIC_Hx2y2z_G3xy_a;
  abcd[28] = 2.0E0*I_KINETIC_Hxy3z_G3xy_a;
  abcd[29] = 2.0E0*I_KINETIC_Hx4z_G3xy_a;
  abcd[30] = 2.0E0*I_KINETIC_H5x_G3xz_a-4*I_KINETIC_F3x_G3xz;
  abcd[31] = 2.0E0*I_KINETIC_H4xy_G3xz_a-3*I_KINETIC_F2xy_G3xz;
  abcd[32] = 2.0E0*I_KINETIC_H4xz_G3xz_a-3*I_KINETIC_F2xz_G3xz;
  abcd[33] = 2.0E0*I_KINETIC_H3x2y_G3xz_a-2*I_KINETIC_Fx2y_G3xz;
  abcd[34] = 2.0E0*I_KINETIC_H3xyz_G3xz_a-2*I_KINETIC_Fxyz_G3xz;
  abcd[35] = 2.0E0*I_KINETIC_H3x2z_G3xz_a-2*I_KINETIC_Fx2z_G3xz;
  abcd[36] = 2.0E0*I_KINETIC_H2x3y_G3xz_a-1*I_KINETIC_F3y_G3xz;
  abcd[37] = 2.0E0*I_KINETIC_H2x2yz_G3xz_a-1*I_KINETIC_F2yz_G3xz;
  abcd[38] = 2.0E0*I_KINETIC_H2xy2z_G3xz_a-1*I_KINETIC_Fy2z_G3xz;
  abcd[39] = 2.0E0*I_KINETIC_H2x3z_G3xz_a-1*I_KINETIC_F3z_G3xz;
  abcd[40] = 2.0E0*I_KINETIC_Hx4y_G3xz_a;
  abcd[41] = 2.0E0*I_KINETIC_Hx3yz_G3xz_a;
  abcd[42] = 2.0E0*I_KINETIC_Hx2y2z_G3xz_a;
  abcd[43] = 2.0E0*I_KINETIC_Hxy3z_G3xz_a;
  abcd[44] = 2.0E0*I_KINETIC_Hx4z_G3xz_a;
  abcd[45] = 2.0E0*I_KINETIC_H5x_G2x2y_a-4*I_KINETIC_F3x_G2x2y;
  abcd[46] = 2.0E0*I_KINETIC_H4xy_G2x2y_a-3*I_KINETIC_F2xy_G2x2y;
  abcd[47] = 2.0E0*I_KINETIC_H4xz_G2x2y_a-3*I_KINETIC_F2xz_G2x2y;
  abcd[48] = 2.0E0*I_KINETIC_H3x2y_G2x2y_a-2*I_KINETIC_Fx2y_G2x2y;
  abcd[49] = 2.0E0*I_KINETIC_H3xyz_G2x2y_a-2*I_KINETIC_Fxyz_G2x2y;
  abcd[50] = 2.0E0*I_KINETIC_H3x2z_G2x2y_a-2*I_KINETIC_Fx2z_G2x2y;
  abcd[51] = 2.0E0*I_KINETIC_H2x3y_G2x2y_a-1*I_KINETIC_F3y_G2x2y;
  abcd[52] = 2.0E0*I_KINETIC_H2x2yz_G2x2y_a-1*I_KINETIC_F2yz_G2x2y;
  abcd[53] = 2.0E0*I_KINETIC_H2xy2z_G2x2y_a-1*I_KINETIC_Fy2z_G2x2y;
  abcd[54] = 2.0E0*I_KINETIC_H2x3z_G2x2y_a-1*I_KINETIC_F3z_G2x2y;
  abcd[55] = 2.0E0*I_KINETIC_Hx4y_G2x2y_a;
  abcd[56] = 2.0E0*I_KINETIC_Hx3yz_G2x2y_a;
  abcd[57] = 2.0E0*I_KINETIC_Hx2y2z_G2x2y_a;
  abcd[58] = 2.0E0*I_KINETIC_Hxy3z_G2x2y_a;
  abcd[59] = 2.0E0*I_KINETIC_Hx4z_G2x2y_a;
  abcd[60] = 2.0E0*I_KINETIC_H5x_G2xyz_a-4*I_KINETIC_F3x_G2xyz;
  abcd[61] = 2.0E0*I_KINETIC_H4xy_G2xyz_a-3*I_KINETIC_F2xy_G2xyz;
  abcd[62] = 2.0E0*I_KINETIC_H4xz_G2xyz_a-3*I_KINETIC_F2xz_G2xyz;
  abcd[63] = 2.0E0*I_KINETIC_H3x2y_G2xyz_a-2*I_KINETIC_Fx2y_G2xyz;
  abcd[64] = 2.0E0*I_KINETIC_H3xyz_G2xyz_a-2*I_KINETIC_Fxyz_G2xyz;
  abcd[65] = 2.0E0*I_KINETIC_H3x2z_G2xyz_a-2*I_KINETIC_Fx2z_G2xyz;
  abcd[66] = 2.0E0*I_KINETIC_H2x3y_G2xyz_a-1*I_KINETIC_F3y_G2xyz;
  abcd[67] = 2.0E0*I_KINETIC_H2x2yz_G2xyz_a-1*I_KINETIC_F2yz_G2xyz;
  abcd[68] = 2.0E0*I_KINETIC_H2xy2z_G2xyz_a-1*I_KINETIC_Fy2z_G2xyz;
  abcd[69] = 2.0E0*I_KINETIC_H2x3z_G2xyz_a-1*I_KINETIC_F3z_G2xyz;
  abcd[70] = 2.0E0*I_KINETIC_Hx4y_G2xyz_a;
  abcd[71] = 2.0E0*I_KINETIC_Hx3yz_G2xyz_a;
  abcd[72] = 2.0E0*I_KINETIC_Hx2y2z_G2xyz_a;
  abcd[73] = 2.0E0*I_KINETIC_Hxy3z_G2xyz_a;
  abcd[74] = 2.0E0*I_KINETIC_Hx4z_G2xyz_a;
  abcd[75] = 2.0E0*I_KINETIC_H5x_G2x2z_a-4*I_KINETIC_F3x_G2x2z;
  abcd[76] = 2.0E0*I_KINETIC_H4xy_G2x2z_a-3*I_KINETIC_F2xy_G2x2z;
  abcd[77] = 2.0E0*I_KINETIC_H4xz_G2x2z_a-3*I_KINETIC_F2xz_G2x2z;
  abcd[78] = 2.0E0*I_KINETIC_H3x2y_G2x2z_a-2*I_KINETIC_Fx2y_G2x2z;
  abcd[79] = 2.0E0*I_KINETIC_H3xyz_G2x2z_a-2*I_KINETIC_Fxyz_G2x2z;
  abcd[80] = 2.0E0*I_KINETIC_H3x2z_G2x2z_a-2*I_KINETIC_Fx2z_G2x2z;
  abcd[81] = 2.0E0*I_KINETIC_H2x3y_G2x2z_a-1*I_KINETIC_F3y_G2x2z;
  abcd[82] = 2.0E0*I_KINETIC_H2x2yz_G2x2z_a-1*I_KINETIC_F2yz_G2x2z;
  abcd[83] = 2.0E0*I_KINETIC_H2xy2z_G2x2z_a-1*I_KINETIC_Fy2z_G2x2z;
  abcd[84] = 2.0E0*I_KINETIC_H2x3z_G2x2z_a-1*I_KINETIC_F3z_G2x2z;
  abcd[85] = 2.0E0*I_KINETIC_Hx4y_G2x2z_a;
  abcd[86] = 2.0E0*I_KINETIC_Hx3yz_G2x2z_a;
  abcd[87] = 2.0E0*I_KINETIC_Hx2y2z_G2x2z_a;
  abcd[88] = 2.0E0*I_KINETIC_Hxy3z_G2x2z_a;
  abcd[89] = 2.0E0*I_KINETIC_Hx4z_G2x2z_a;
  abcd[90] = 2.0E0*I_KINETIC_H5x_Gx3y_a-4*I_KINETIC_F3x_Gx3y;
  abcd[91] = 2.0E0*I_KINETIC_H4xy_Gx3y_a-3*I_KINETIC_F2xy_Gx3y;
  abcd[92] = 2.0E0*I_KINETIC_H4xz_Gx3y_a-3*I_KINETIC_F2xz_Gx3y;
  abcd[93] = 2.0E0*I_KINETIC_H3x2y_Gx3y_a-2*I_KINETIC_Fx2y_Gx3y;
  abcd[94] = 2.0E0*I_KINETIC_H3xyz_Gx3y_a-2*I_KINETIC_Fxyz_Gx3y;
  abcd[95] = 2.0E0*I_KINETIC_H3x2z_Gx3y_a-2*I_KINETIC_Fx2z_Gx3y;
  abcd[96] = 2.0E0*I_KINETIC_H2x3y_Gx3y_a-1*I_KINETIC_F3y_Gx3y;
  abcd[97] = 2.0E0*I_KINETIC_H2x2yz_Gx3y_a-1*I_KINETIC_F2yz_Gx3y;
  abcd[98] = 2.0E0*I_KINETIC_H2xy2z_Gx3y_a-1*I_KINETIC_Fy2z_Gx3y;
  abcd[99] = 2.0E0*I_KINETIC_H2x3z_Gx3y_a-1*I_KINETIC_F3z_Gx3y;
  abcd[100] = 2.0E0*I_KINETIC_Hx4y_Gx3y_a;
  abcd[101] = 2.0E0*I_KINETIC_Hx3yz_Gx3y_a;
  abcd[102] = 2.0E0*I_KINETIC_Hx2y2z_Gx3y_a;
  abcd[103] = 2.0E0*I_KINETIC_Hxy3z_Gx3y_a;
  abcd[104] = 2.0E0*I_KINETIC_Hx4z_Gx3y_a;
  abcd[105] = 2.0E0*I_KINETIC_H5x_Gx2yz_a-4*I_KINETIC_F3x_Gx2yz;
  abcd[106] = 2.0E0*I_KINETIC_H4xy_Gx2yz_a-3*I_KINETIC_F2xy_Gx2yz;
  abcd[107] = 2.0E0*I_KINETIC_H4xz_Gx2yz_a-3*I_KINETIC_F2xz_Gx2yz;
  abcd[108] = 2.0E0*I_KINETIC_H3x2y_Gx2yz_a-2*I_KINETIC_Fx2y_Gx2yz;
  abcd[109] = 2.0E0*I_KINETIC_H3xyz_Gx2yz_a-2*I_KINETIC_Fxyz_Gx2yz;
  abcd[110] = 2.0E0*I_KINETIC_H3x2z_Gx2yz_a-2*I_KINETIC_Fx2z_Gx2yz;
  abcd[111] = 2.0E0*I_KINETIC_H2x3y_Gx2yz_a-1*I_KINETIC_F3y_Gx2yz;
  abcd[112] = 2.0E0*I_KINETIC_H2x2yz_Gx2yz_a-1*I_KINETIC_F2yz_Gx2yz;
  abcd[113] = 2.0E0*I_KINETIC_H2xy2z_Gx2yz_a-1*I_KINETIC_Fy2z_Gx2yz;
  abcd[114] = 2.0E0*I_KINETIC_H2x3z_Gx2yz_a-1*I_KINETIC_F3z_Gx2yz;
  abcd[115] = 2.0E0*I_KINETIC_Hx4y_Gx2yz_a;
  abcd[116] = 2.0E0*I_KINETIC_Hx3yz_Gx2yz_a;
  abcd[117] = 2.0E0*I_KINETIC_Hx2y2z_Gx2yz_a;
  abcd[118] = 2.0E0*I_KINETIC_Hxy3z_Gx2yz_a;
  abcd[119] = 2.0E0*I_KINETIC_Hx4z_Gx2yz_a;
  abcd[120] = 2.0E0*I_KINETIC_H5x_Gxy2z_a-4*I_KINETIC_F3x_Gxy2z;
  abcd[121] = 2.0E0*I_KINETIC_H4xy_Gxy2z_a-3*I_KINETIC_F2xy_Gxy2z;
  abcd[122] = 2.0E0*I_KINETIC_H4xz_Gxy2z_a-3*I_KINETIC_F2xz_Gxy2z;
  abcd[123] = 2.0E0*I_KINETIC_H3x2y_Gxy2z_a-2*I_KINETIC_Fx2y_Gxy2z;
  abcd[124] = 2.0E0*I_KINETIC_H3xyz_Gxy2z_a-2*I_KINETIC_Fxyz_Gxy2z;
  abcd[125] = 2.0E0*I_KINETIC_H3x2z_Gxy2z_a-2*I_KINETIC_Fx2z_Gxy2z;
  abcd[126] = 2.0E0*I_KINETIC_H2x3y_Gxy2z_a-1*I_KINETIC_F3y_Gxy2z;
  abcd[127] = 2.0E0*I_KINETIC_H2x2yz_Gxy2z_a-1*I_KINETIC_F2yz_Gxy2z;
  abcd[128] = 2.0E0*I_KINETIC_H2xy2z_Gxy2z_a-1*I_KINETIC_Fy2z_Gxy2z;
  abcd[129] = 2.0E0*I_KINETIC_H2x3z_Gxy2z_a-1*I_KINETIC_F3z_Gxy2z;
  abcd[130] = 2.0E0*I_KINETIC_Hx4y_Gxy2z_a;
  abcd[131] = 2.0E0*I_KINETIC_Hx3yz_Gxy2z_a;
  abcd[132] = 2.0E0*I_KINETIC_Hx2y2z_Gxy2z_a;
  abcd[133] = 2.0E0*I_KINETIC_Hxy3z_Gxy2z_a;
  abcd[134] = 2.0E0*I_KINETIC_Hx4z_Gxy2z_a;
  abcd[135] = 2.0E0*I_KINETIC_H5x_Gx3z_a-4*I_KINETIC_F3x_Gx3z;
  abcd[136] = 2.0E0*I_KINETIC_H4xy_Gx3z_a-3*I_KINETIC_F2xy_Gx3z;
  abcd[137] = 2.0E0*I_KINETIC_H4xz_Gx3z_a-3*I_KINETIC_F2xz_Gx3z;
  abcd[138] = 2.0E0*I_KINETIC_H3x2y_Gx3z_a-2*I_KINETIC_Fx2y_Gx3z;
  abcd[139] = 2.0E0*I_KINETIC_H3xyz_Gx3z_a-2*I_KINETIC_Fxyz_Gx3z;
  abcd[140] = 2.0E0*I_KINETIC_H3x2z_Gx3z_a-2*I_KINETIC_Fx2z_Gx3z;
  abcd[141] = 2.0E0*I_KINETIC_H2x3y_Gx3z_a-1*I_KINETIC_F3y_Gx3z;
  abcd[142] = 2.0E0*I_KINETIC_H2x2yz_Gx3z_a-1*I_KINETIC_F2yz_Gx3z;
  abcd[143] = 2.0E0*I_KINETIC_H2xy2z_Gx3z_a-1*I_KINETIC_Fy2z_Gx3z;
  abcd[144] = 2.0E0*I_KINETIC_H2x3z_Gx3z_a-1*I_KINETIC_F3z_Gx3z;
  abcd[145] = 2.0E0*I_KINETIC_Hx4y_Gx3z_a;
  abcd[146] = 2.0E0*I_KINETIC_Hx3yz_Gx3z_a;
  abcd[147] = 2.0E0*I_KINETIC_Hx2y2z_Gx3z_a;
  abcd[148] = 2.0E0*I_KINETIC_Hxy3z_Gx3z_a;
  abcd[149] = 2.0E0*I_KINETIC_Hx4z_Gx3z_a;
  abcd[150] = 2.0E0*I_KINETIC_H5x_G4y_a-4*I_KINETIC_F3x_G4y;
  abcd[151] = 2.0E0*I_KINETIC_H4xy_G4y_a-3*I_KINETIC_F2xy_G4y;
  abcd[152] = 2.0E0*I_KINETIC_H4xz_G4y_a-3*I_KINETIC_F2xz_G4y;
  abcd[153] = 2.0E0*I_KINETIC_H3x2y_G4y_a-2*I_KINETIC_Fx2y_G4y;
  abcd[154] = 2.0E0*I_KINETIC_H3xyz_G4y_a-2*I_KINETIC_Fxyz_G4y;
  abcd[155] = 2.0E0*I_KINETIC_H3x2z_G4y_a-2*I_KINETIC_Fx2z_G4y;
  abcd[156] = 2.0E0*I_KINETIC_H2x3y_G4y_a-1*I_KINETIC_F3y_G4y;
  abcd[157] = 2.0E0*I_KINETIC_H2x2yz_G4y_a-1*I_KINETIC_F2yz_G4y;
  abcd[158] = 2.0E0*I_KINETIC_H2xy2z_G4y_a-1*I_KINETIC_Fy2z_G4y;
  abcd[159] = 2.0E0*I_KINETIC_H2x3z_G4y_a-1*I_KINETIC_F3z_G4y;
  abcd[160] = 2.0E0*I_KINETIC_Hx4y_G4y_a;
  abcd[161] = 2.0E0*I_KINETIC_Hx3yz_G4y_a;
  abcd[162] = 2.0E0*I_KINETIC_Hx2y2z_G4y_a;
  abcd[163] = 2.0E0*I_KINETIC_Hxy3z_G4y_a;
  abcd[164] = 2.0E0*I_KINETIC_Hx4z_G4y_a;
  abcd[165] = 2.0E0*I_KINETIC_H5x_G3yz_a-4*I_KINETIC_F3x_G3yz;
  abcd[166] = 2.0E0*I_KINETIC_H4xy_G3yz_a-3*I_KINETIC_F2xy_G3yz;
  abcd[167] = 2.0E0*I_KINETIC_H4xz_G3yz_a-3*I_KINETIC_F2xz_G3yz;
  abcd[168] = 2.0E0*I_KINETIC_H3x2y_G3yz_a-2*I_KINETIC_Fx2y_G3yz;
  abcd[169] = 2.0E0*I_KINETIC_H3xyz_G3yz_a-2*I_KINETIC_Fxyz_G3yz;
  abcd[170] = 2.0E0*I_KINETIC_H3x2z_G3yz_a-2*I_KINETIC_Fx2z_G3yz;
  abcd[171] = 2.0E0*I_KINETIC_H2x3y_G3yz_a-1*I_KINETIC_F3y_G3yz;
  abcd[172] = 2.0E0*I_KINETIC_H2x2yz_G3yz_a-1*I_KINETIC_F2yz_G3yz;
  abcd[173] = 2.0E0*I_KINETIC_H2xy2z_G3yz_a-1*I_KINETIC_Fy2z_G3yz;
  abcd[174] = 2.0E0*I_KINETIC_H2x3z_G3yz_a-1*I_KINETIC_F3z_G3yz;
  abcd[175] = 2.0E0*I_KINETIC_Hx4y_G3yz_a;
  abcd[176] = 2.0E0*I_KINETIC_Hx3yz_G3yz_a;
  abcd[177] = 2.0E0*I_KINETIC_Hx2y2z_G3yz_a;
  abcd[178] = 2.0E0*I_KINETIC_Hxy3z_G3yz_a;
  abcd[179] = 2.0E0*I_KINETIC_Hx4z_G3yz_a;
  abcd[180] = 2.0E0*I_KINETIC_H5x_G2y2z_a-4*I_KINETIC_F3x_G2y2z;
  abcd[181] = 2.0E0*I_KINETIC_H4xy_G2y2z_a-3*I_KINETIC_F2xy_G2y2z;
  abcd[182] = 2.0E0*I_KINETIC_H4xz_G2y2z_a-3*I_KINETIC_F2xz_G2y2z;
  abcd[183] = 2.0E0*I_KINETIC_H3x2y_G2y2z_a-2*I_KINETIC_Fx2y_G2y2z;
  abcd[184] = 2.0E0*I_KINETIC_H3xyz_G2y2z_a-2*I_KINETIC_Fxyz_G2y2z;
  abcd[185] = 2.0E0*I_KINETIC_H3x2z_G2y2z_a-2*I_KINETIC_Fx2z_G2y2z;
  abcd[186] = 2.0E0*I_KINETIC_H2x3y_G2y2z_a-1*I_KINETIC_F3y_G2y2z;
  abcd[187] = 2.0E0*I_KINETIC_H2x2yz_G2y2z_a-1*I_KINETIC_F2yz_G2y2z;
  abcd[188] = 2.0E0*I_KINETIC_H2xy2z_G2y2z_a-1*I_KINETIC_Fy2z_G2y2z;
  abcd[189] = 2.0E0*I_KINETIC_H2x3z_G2y2z_a-1*I_KINETIC_F3z_G2y2z;
  abcd[190] = 2.0E0*I_KINETIC_Hx4y_G2y2z_a;
  abcd[191] = 2.0E0*I_KINETIC_Hx3yz_G2y2z_a;
  abcd[192] = 2.0E0*I_KINETIC_Hx2y2z_G2y2z_a;
  abcd[193] = 2.0E0*I_KINETIC_Hxy3z_G2y2z_a;
  abcd[194] = 2.0E0*I_KINETIC_Hx4z_G2y2z_a;
  abcd[195] = 2.0E0*I_KINETIC_H5x_Gy3z_a-4*I_KINETIC_F3x_Gy3z;
  abcd[196] = 2.0E0*I_KINETIC_H4xy_Gy3z_a-3*I_KINETIC_F2xy_Gy3z;
  abcd[197] = 2.0E0*I_KINETIC_H4xz_Gy3z_a-3*I_KINETIC_F2xz_Gy3z;
  abcd[198] = 2.0E0*I_KINETIC_H3x2y_Gy3z_a-2*I_KINETIC_Fx2y_Gy3z;
  abcd[199] = 2.0E0*I_KINETIC_H3xyz_Gy3z_a-2*I_KINETIC_Fxyz_Gy3z;
  abcd[200] = 2.0E0*I_KINETIC_H3x2z_Gy3z_a-2*I_KINETIC_Fx2z_Gy3z;
  abcd[201] = 2.0E0*I_KINETIC_H2x3y_Gy3z_a-1*I_KINETIC_F3y_Gy3z;
  abcd[202] = 2.0E0*I_KINETIC_H2x2yz_Gy3z_a-1*I_KINETIC_F2yz_Gy3z;
  abcd[203] = 2.0E0*I_KINETIC_H2xy2z_Gy3z_a-1*I_KINETIC_Fy2z_Gy3z;
  abcd[204] = 2.0E0*I_KINETIC_H2x3z_Gy3z_a-1*I_KINETIC_F3z_Gy3z;
  abcd[205] = 2.0E0*I_KINETIC_Hx4y_Gy3z_a;
  abcd[206] = 2.0E0*I_KINETIC_Hx3yz_Gy3z_a;
  abcd[207] = 2.0E0*I_KINETIC_Hx2y2z_Gy3z_a;
  abcd[208] = 2.0E0*I_KINETIC_Hxy3z_Gy3z_a;
  abcd[209] = 2.0E0*I_KINETIC_Hx4z_Gy3z_a;
  abcd[210] = 2.0E0*I_KINETIC_H5x_G4z_a-4*I_KINETIC_F3x_G4z;
  abcd[211] = 2.0E0*I_KINETIC_H4xy_G4z_a-3*I_KINETIC_F2xy_G4z;
  abcd[212] = 2.0E0*I_KINETIC_H4xz_G4z_a-3*I_KINETIC_F2xz_G4z;
  abcd[213] = 2.0E0*I_KINETIC_H3x2y_G4z_a-2*I_KINETIC_Fx2y_G4z;
  abcd[214] = 2.0E0*I_KINETIC_H3xyz_G4z_a-2*I_KINETIC_Fxyz_G4z;
  abcd[215] = 2.0E0*I_KINETIC_H3x2z_G4z_a-2*I_KINETIC_Fx2z_G4z;
  abcd[216] = 2.0E0*I_KINETIC_H2x3y_G4z_a-1*I_KINETIC_F3y_G4z;
  abcd[217] = 2.0E0*I_KINETIC_H2x2yz_G4z_a-1*I_KINETIC_F2yz_G4z;
  abcd[218] = 2.0E0*I_KINETIC_H2xy2z_G4z_a-1*I_KINETIC_Fy2z_G4z;
  abcd[219] = 2.0E0*I_KINETIC_H2x3z_G4z_a-1*I_KINETIC_F3z_G4z;
  abcd[220] = 2.0E0*I_KINETIC_Hx4y_G4z_a;
  abcd[221] = 2.0E0*I_KINETIC_Hx3yz_G4z_a;
  abcd[222] = 2.0E0*I_KINETIC_Hx2y2z_G4z_a;
  abcd[223] = 2.0E0*I_KINETIC_Hxy3z_G4z_a;
  abcd[224] = 2.0E0*I_KINETIC_Hx4z_G4z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_G_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_G_a
   * RHS shell quartet name: SQ_KINETIC_F_G
   ************************************************************/
  abcd[225] = 2.0E0*I_KINETIC_H4xy_G4x_a;
  abcd[226] = 2.0E0*I_KINETIC_H3x2y_G4x_a-1*I_KINETIC_F3x_G4x;
  abcd[227] = 2.0E0*I_KINETIC_H3xyz_G4x_a;
  abcd[228] = 2.0E0*I_KINETIC_H2x3y_G4x_a-2*I_KINETIC_F2xy_G4x;
  abcd[229] = 2.0E0*I_KINETIC_H2x2yz_G4x_a-1*I_KINETIC_F2xz_G4x;
  abcd[230] = 2.0E0*I_KINETIC_H2xy2z_G4x_a;
  abcd[231] = 2.0E0*I_KINETIC_Hx4y_G4x_a-3*I_KINETIC_Fx2y_G4x;
  abcd[232] = 2.0E0*I_KINETIC_Hx3yz_G4x_a-2*I_KINETIC_Fxyz_G4x;
  abcd[233] = 2.0E0*I_KINETIC_Hx2y2z_G4x_a-1*I_KINETIC_Fx2z_G4x;
  abcd[234] = 2.0E0*I_KINETIC_Hxy3z_G4x_a;
  abcd[235] = 2.0E0*I_KINETIC_H5y_G4x_a-4*I_KINETIC_F3y_G4x;
  abcd[236] = 2.0E0*I_KINETIC_H4yz_G4x_a-3*I_KINETIC_F2yz_G4x;
  abcd[237] = 2.0E0*I_KINETIC_H3y2z_G4x_a-2*I_KINETIC_Fy2z_G4x;
  abcd[238] = 2.0E0*I_KINETIC_H2y3z_G4x_a-1*I_KINETIC_F3z_G4x;
  abcd[239] = 2.0E0*I_KINETIC_Hy4z_G4x_a;
  abcd[240] = 2.0E0*I_KINETIC_H4xy_G3xy_a;
  abcd[241] = 2.0E0*I_KINETIC_H3x2y_G3xy_a-1*I_KINETIC_F3x_G3xy;
  abcd[242] = 2.0E0*I_KINETIC_H3xyz_G3xy_a;
  abcd[243] = 2.0E0*I_KINETIC_H2x3y_G3xy_a-2*I_KINETIC_F2xy_G3xy;
  abcd[244] = 2.0E0*I_KINETIC_H2x2yz_G3xy_a-1*I_KINETIC_F2xz_G3xy;
  abcd[245] = 2.0E0*I_KINETIC_H2xy2z_G3xy_a;
  abcd[246] = 2.0E0*I_KINETIC_Hx4y_G3xy_a-3*I_KINETIC_Fx2y_G3xy;
  abcd[247] = 2.0E0*I_KINETIC_Hx3yz_G3xy_a-2*I_KINETIC_Fxyz_G3xy;
  abcd[248] = 2.0E0*I_KINETIC_Hx2y2z_G3xy_a-1*I_KINETIC_Fx2z_G3xy;
  abcd[249] = 2.0E0*I_KINETIC_Hxy3z_G3xy_a;
  abcd[250] = 2.0E0*I_KINETIC_H5y_G3xy_a-4*I_KINETIC_F3y_G3xy;
  abcd[251] = 2.0E0*I_KINETIC_H4yz_G3xy_a-3*I_KINETIC_F2yz_G3xy;
  abcd[252] = 2.0E0*I_KINETIC_H3y2z_G3xy_a-2*I_KINETIC_Fy2z_G3xy;
  abcd[253] = 2.0E0*I_KINETIC_H2y3z_G3xy_a-1*I_KINETIC_F3z_G3xy;
  abcd[254] = 2.0E0*I_KINETIC_Hy4z_G3xy_a;
  abcd[255] = 2.0E0*I_KINETIC_H4xy_G3xz_a;
  abcd[256] = 2.0E0*I_KINETIC_H3x2y_G3xz_a-1*I_KINETIC_F3x_G3xz;
  abcd[257] = 2.0E0*I_KINETIC_H3xyz_G3xz_a;
  abcd[258] = 2.0E0*I_KINETIC_H2x3y_G3xz_a-2*I_KINETIC_F2xy_G3xz;
  abcd[259] = 2.0E0*I_KINETIC_H2x2yz_G3xz_a-1*I_KINETIC_F2xz_G3xz;
  abcd[260] = 2.0E0*I_KINETIC_H2xy2z_G3xz_a;
  abcd[261] = 2.0E0*I_KINETIC_Hx4y_G3xz_a-3*I_KINETIC_Fx2y_G3xz;
  abcd[262] = 2.0E0*I_KINETIC_Hx3yz_G3xz_a-2*I_KINETIC_Fxyz_G3xz;
  abcd[263] = 2.0E0*I_KINETIC_Hx2y2z_G3xz_a-1*I_KINETIC_Fx2z_G3xz;
  abcd[264] = 2.0E0*I_KINETIC_Hxy3z_G3xz_a;
  abcd[265] = 2.0E0*I_KINETIC_H5y_G3xz_a-4*I_KINETIC_F3y_G3xz;
  abcd[266] = 2.0E0*I_KINETIC_H4yz_G3xz_a-3*I_KINETIC_F2yz_G3xz;
  abcd[267] = 2.0E0*I_KINETIC_H3y2z_G3xz_a-2*I_KINETIC_Fy2z_G3xz;
  abcd[268] = 2.0E0*I_KINETIC_H2y3z_G3xz_a-1*I_KINETIC_F3z_G3xz;
  abcd[269] = 2.0E0*I_KINETIC_Hy4z_G3xz_a;
  abcd[270] = 2.0E0*I_KINETIC_H4xy_G2x2y_a;
  abcd[271] = 2.0E0*I_KINETIC_H3x2y_G2x2y_a-1*I_KINETIC_F3x_G2x2y;
  abcd[272] = 2.0E0*I_KINETIC_H3xyz_G2x2y_a;
  abcd[273] = 2.0E0*I_KINETIC_H2x3y_G2x2y_a-2*I_KINETIC_F2xy_G2x2y;
  abcd[274] = 2.0E0*I_KINETIC_H2x2yz_G2x2y_a-1*I_KINETIC_F2xz_G2x2y;
  abcd[275] = 2.0E0*I_KINETIC_H2xy2z_G2x2y_a;
  abcd[276] = 2.0E0*I_KINETIC_Hx4y_G2x2y_a-3*I_KINETIC_Fx2y_G2x2y;
  abcd[277] = 2.0E0*I_KINETIC_Hx3yz_G2x2y_a-2*I_KINETIC_Fxyz_G2x2y;
  abcd[278] = 2.0E0*I_KINETIC_Hx2y2z_G2x2y_a-1*I_KINETIC_Fx2z_G2x2y;
  abcd[279] = 2.0E0*I_KINETIC_Hxy3z_G2x2y_a;
  abcd[280] = 2.0E0*I_KINETIC_H5y_G2x2y_a-4*I_KINETIC_F3y_G2x2y;
  abcd[281] = 2.0E0*I_KINETIC_H4yz_G2x2y_a-3*I_KINETIC_F2yz_G2x2y;
  abcd[282] = 2.0E0*I_KINETIC_H3y2z_G2x2y_a-2*I_KINETIC_Fy2z_G2x2y;
  abcd[283] = 2.0E0*I_KINETIC_H2y3z_G2x2y_a-1*I_KINETIC_F3z_G2x2y;
  abcd[284] = 2.0E0*I_KINETIC_Hy4z_G2x2y_a;
  abcd[285] = 2.0E0*I_KINETIC_H4xy_G2xyz_a;
  abcd[286] = 2.0E0*I_KINETIC_H3x2y_G2xyz_a-1*I_KINETIC_F3x_G2xyz;
  abcd[287] = 2.0E0*I_KINETIC_H3xyz_G2xyz_a;
  abcd[288] = 2.0E0*I_KINETIC_H2x3y_G2xyz_a-2*I_KINETIC_F2xy_G2xyz;
  abcd[289] = 2.0E0*I_KINETIC_H2x2yz_G2xyz_a-1*I_KINETIC_F2xz_G2xyz;
  abcd[290] = 2.0E0*I_KINETIC_H2xy2z_G2xyz_a;
  abcd[291] = 2.0E0*I_KINETIC_Hx4y_G2xyz_a-3*I_KINETIC_Fx2y_G2xyz;
  abcd[292] = 2.0E0*I_KINETIC_Hx3yz_G2xyz_a-2*I_KINETIC_Fxyz_G2xyz;
  abcd[293] = 2.0E0*I_KINETIC_Hx2y2z_G2xyz_a-1*I_KINETIC_Fx2z_G2xyz;
  abcd[294] = 2.0E0*I_KINETIC_Hxy3z_G2xyz_a;
  abcd[295] = 2.0E0*I_KINETIC_H5y_G2xyz_a-4*I_KINETIC_F3y_G2xyz;
  abcd[296] = 2.0E0*I_KINETIC_H4yz_G2xyz_a-3*I_KINETIC_F2yz_G2xyz;
  abcd[297] = 2.0E0*I_KINETIC_H3y2z_G2xyz_a-2*I_KINETIC_Fy2z_G2xyz;
  abcd[298] = 2.0E0*I_KINETIC_H2y3z_G2xyz_a-1*I_KINETIC_F3z_G2xyz;
  abcd[299] = 2.0E0*I_KINETIC_Hy4z_G2xyz_a;
  abcd[300] = 2.0E0*I_KINETIC_H4xy_G2x2z_a;
  abcd[301] = 2.0E0*I_KINETIC_H3x2y_G2x2z_a-1*I_KINETIC_F3x_G2x2z;
  abcd[302] = 2.0E0*I_KINETIC_H3xyz_G2x2z_a;
  abcd[303] = 2.0E0*I_KINETIC_H2x3y_G2x2z_a-2*I_KINETIC_F2xy_G2x2z;
  abcd[304] = 2.0E0*I_KINETIC_H2x2yz_G2x2z_a-1*I_KINETIC_F2xz_G2x2z;
  abcd[305] = 2.0E0*I_KINETIC_H2xy2z_G2x2z_a;
  abcd[306] = 2.0E0*I_KINETIC_Hx4y_G2x2z_a-3*I_KINETIC_Fx2y_G2x2z;
  abcd[307] = 2.0E0*I_KINETIC_Hx3yz_G2x2z_a-2*I_KINETIC_Fxyz_G2x2z;
  abcd[308] = 2.0E0*I_KINETIC_Hx2y2z_G2x2z_a-1*I_KINETIC_Fx2z_G2x2z;
  abcd[309] = 2.0E0*I_KINETIC_Hxy3z_G2x2z_a;
  abcd[310] = 2.0E0*I_KINETIC_H5y_G2x2z_a-4*I_KINETIC_F3y_G2x2z;
  abcd[311] = 2.0E0*I_KINETIC_H4yz_G2x2z_a-3*I_KINETIC_F2yz_G2x2z;
  abcd[312] = 2.0E0*I_KINETIC_H3y2z_G2x2z_a-2*I_KINETIC_Fy2z_G2x2z;
  abcd[313] = 2.0E0*I_KINETIC_H2y3z_G2x2z_a-1*I_KINETIC_F3z_G2x2z;
  abcd[314] = 2.0E0*I_KINETIC_Hy4z_G2x2z_a;
  abcd[315] = 2.0E0*I_KINETIC_H4xy_Gx3y_a;
  abcd[316] = 2.0E0*I_KINETIC_H3x2y_Gx3y_a-1*I_KINETIC_F3x_Gx3y;
  abcd[317] = 2.0E0*I_KINETIC_H3xyz_Gx3y_a;
  abcd[318] = 2.0E0*I_KINETIC_H2x3y_Gx3y_a-2*I_KINETIC_F2xy_Gx3y;
  abcd[319] = 2.0E0*I_KINETIC_H2x2yz_Gx3y_a-1*I_KINETIC_F2xz_Gx3y;
  abcd[320] = 2.0E0*I_KINETIC_H2xy2z_Gx3y_a;
  abcd[321] = 2.0E0*I_KINETIC_Hx4y_Gx3y_a-3*I_KINETIC_Fx2y_Gx3y;
  abcd[322] = 2.0E0*I_KINETIC_Hx3yz_Gx3y_a-2*I_KINETIC_Fxyz_Gx3y;
  abcd[323] = 2.0E0*I_KINETIC_Hx2y2z_Gx3y_a-1*I_KINETIC_Fx2z_Gx3y;
  abcd[324] = 2.0E0*I_KINETIC_Hxy3z_Gx3y_a;
  abcd[325] = 2.0E0*I_KINETIC_H5y_Gx3y_a-4*I_KINETIC_F3y_Gx3y;
  abcd[326] = 2.0E0*I_KINETIC_H4yz_Gx3y_a-3*I_KINETIC_F2yz_Gx3y;
  abcd[327] = 2.0E0*I_KINETIC_H3y2z_Gx3y_a-2*I_KINETIC_Fy2z_Gx3y;
  abcd[328] = 2.0E0*I_KINETIC_H2y3z_Gx3y_a-1*I_KINETIC_F3z_Gx3y;
  abcd[329] = 2.0E0*I_KINETIC_Hy4z_Gx3y_a;
  abcd[330] = 2.0E0*I_KINETIC_H4xy_Gx2yz_a;
  abcd[331] = 2.0E0*I_KINETIC_H3x2y_Gx2yz_a-1*I_KINETIC_F3x_Gx2yz;
  abcd[332] = 2.0E0*I_KINETIC_H3xyz_Gx2yz_a;
  abcd[333] = 2.0E0*I_KINETIC_H2x3y_Gx2yz_a-2*I_KINETIC_F2xy_Gx2yz;
  abcd[334] = 2.0E0*I_KINETIC_H2x2yz_Gx2yz_a-1*I_KINETIC_F2xz_Gx2yz;
  abcd[335] = 2.0E0*I_KINETIC_H2xy2z_Gx2yz_a;
  abcd[336] = 2.0E0*I_KINETIC_Hx4y_Gx2yz_a-3*I_KINETIC_Fx2y_Gx2yz;
  abcd[337] = 2.0E0*I_KINETIC_Hx3yz_Gx2yz_a-2*I_KINETIC_Fxyz_Gx2yz;
  abcd[338] = 2.0E0*I_KINETIC_Hx2y2z_Gx2yz_a-1*I_KINETIC_Fx2z_Gx2yz;
  abcd[339] = 2.0E0*I_KINETIC_Hxy3z_Gx2yz_a;
  abcd[340] = 2.0E0*I_KINETIC_H5y_Gx2yz_a-4*I_KINETIC_F3y_Gx2yz;
  abcd[341] = 2.0E0*I_KINETIC_H4yz_Gx2yz_a-3*I_KINETIC_F2yz_Gx2yz;
  abcd[342] = 2.0E0*I_KINETIC_H3y2z_Gx2yz_a-2*I_KINETIC_Fy2z_Gx2yz;
  abcd[343] = 2.0E0*I_KINETIC_H2y3z_Gx2yz_a-1*I_KINETIC_F3z_Gx2yz;
  abcd[344] = 2.0E0*I_KINETIC_Hy4z_Gx2yz_a;
  abcd[345] = 2.0E0*I_KINETIC_H4xy_Gxy2z_a;
  abcd[346] = 2.0E0*I_KINETIC_H3x2y_Gxy2z_a-1*I_KINETIC_F3x_Gxy2z;
  abcd[347] = 2.0E0*I_KINETIC_H3xyz_Gxy2z_a;
  abcd[348] = 2.0E0*I_KINETIC_H2x3y_Gxy2z_a-2*I_KINETIC_F2xy_Gxy2z;
  abcd[349] = 2.0E0*I_KINETIC_H2x2yz_Gxy2z_a-1*I_KINETIC_F2xz_Gxy2z;
  abcd[350] = 2.0E0*I_KINETIC_H2xy2z_Gxy2z_a;
  abcd[351] = 2.0E0*I_KINETIC_Hx4y_Gxy2z_a-3*I_KINETIC_Fx2y_Gxy2z;
  abcd[352] = 2.0E0*I_KINETIC_Hx3yz_Gxy2z_a-2*I_KINETIC_Fxyz_Gxy2z;
  abcd[353] = 2.0E0*I_KINETIC_Hx2y2z_Gxy2z_a-1*I_KINETIC_Fx2z_Gxy2z;
  abcd[354] = 2.0E0*I_KINETIC_Hxy3z_Gxy2z_a;
  abcd[355] = 2.0E0*I_KINETIC_H5y_Gxy2z_a-4*I_KINETIC_F3y_Gxy2z;
  abcd[356] = 2.0E0*I_KINETIC_H4yz_Gxy2z_a-3*I_KINETIC_F2yz_Gxy2z;
  abcd[357] = 2.0E0*I_KINETIC_H3y2z_Gxy2z_a-2*I_KINETIC_Fy2z_Gxy2z;
  abcd[358] = 2.0E0*I_KINETIC_H2y3z_Gxy2z_a-1*I_KINETIC_F3z_Gxy2z;
  abcd[359] = 2.0E0*I_KINETIC_Hy4z_Gxy2z_a;
  abcd[360] = 2.0E0*I_KINETIC_H4xy_Gx3z_a;
  abcd[361] = 2.0E0*I_KINETIC_H3x2y_Gx3z_a-1*I_KINETIC_F3x_Gx3z;
  abcd[362] = 2.0E0*I_KINETIC_H3xyz_Gx3z_a;
  abcd[363] = 2.0E0*I_KINETIC_H2x3y_Gx3z_a-2*I_KINETIC_F2xy_Gx3z;
  abcd[364] = 2.0E0*I_KINETIC_H2x2yz_Gx3z_a-1*I_KINETIC_F2xz_Gx3z;
  abcd[365] = 2.0E0*I_KINETIC_H2xy2z_Gx3z_a;
  abcd[366] = 2.0E0*I_KINETIC_Hx4y_Gx3z_a-3*I_KINETIC_Fx2y_Gx3z;
  abcd[367] = 2.0E0*I_KINETIC_Hx3yz_Gx3z_a-2*I_KINETIC_Fxyz_Gx3z;
  abcd[368] = 2.0E0*I_KINETIC_Hx2y2z_Gx3z_a-1*I_KINETIC_Fx2z_Gx3z;
  abcd[369] = 2.0E0*I_KINETIC_Hxy3z_Gx3z_a;
  abcd[370] = 2.0E0*I_KINETIC_H5y_Gx3z_a-4*I_KINETIC_F3y_Gx3z;
  abcd[371] = 2.0E0*I_KINETIC_H4yz_Gx3z_a-3*I_KINETIC_F2yz_Gx3z;
  abcd[372] = 2.0E0*I_KINETIC_H3y2z_Gx3z_a-2*I_KINETIC_Fy2z_Gx3z;
  abcd[373] = 2.0E0*I_KINETIC_H2y3z_Gx3z_a-1*I_KINETIC_F3z_Gx3z;
  abcd[374] = 2.0E0*I_KINETIC_Hy4z_Gx3z_a;
  abcd[375] = 2.0E0*I_KINETIC_H4xy_G4y_a;
  abcd[376] = 2.0E0*I_KINETIC_H3x2y_G4y_a-1*I_KINETIC_F3x_G4y;
  abcd[377] = 2.0E0*I_KINETIC_H3xyz_G4y_a;
  abcd[378] = 2.0E0*I_KINETIC_H2x3y_G4y_a-2*I_KINETIC_F2xy_G4y;
  abcd[379] = 2.0E0*I_KINETIC_H2x2yz_G4y_a-1*I_KINETIC_F2xz_G4y;
  abcd[380] = 2.0E0*I_KINETIC_H2xy2z_G4y_a;
  abcd[381] = 2.0E0*I_KINETIC_Hx4y_G4y_a-3*I_KINETIC_Fx2y_G4y;
  abcd[382] = 2.0E0*I_KINETIC_Hx3yz_G4y_a-2*I_KINETIC_Fxyz_G4y;
  abcd[383] = 2.0E0*I_KINETIC_Hx2y2z_G4y_a-1*I_KINETIC_Fx2z_G4y;
  abcd[384] = 2.0E0*I_KINETIC_Hxy3z_G4y_a;
  abcd[385] = 2.0E0*I_KINETIC_H5y_G4y_a-4*I_KINETIC_F3y_G4y;
  abcd[386] = 2.0E0*I_KINETIC_H4yz_G4y_a-3*I_KINETIC_F2yz_G4y;
  abcd[387] = 2.0E0*I_KINETIC_H3y2z_G4y_a-2*I_KINETIC_Fy2z_G4y;
  abcd[388] = 2.0E0*I_KINETIC_H2y3z_G4y_a-1*I_KINETIC_F3z_G4y;
  abcd[389] = 2.0E0*I_KINETIC_Hy4z_G4y_a;
  abcd[390] = 2.0E0*I_KINETIC_H4xy_G3yz_a;
  abcd[391] = 2.0E0*I_KINETIC_H3x2y_G3yz_a-1*I_KINETIC_F3x_G3yz;
  abcd[392] = 2.0E0*I_KINETIC_H3xyz_G3yz_a;
  abcd[393] = 2.0E0*I_KINETIC_H2x3y_G3yz_a-2*I_KINETIC_F2xy_G3yz;
  abcd[394] = 2.0E0*I_KINETIC_H2x2yz_G3yz_a-1*I_KINETIC_F2xz_G3yz;
  abcd[395] = 2.0E0*I_KINETIC_H2xy2z_G3yz_a;
  abcd[396] = 2.0E0*I_KINETIC_Hx4y_G3yz_a-3*I_KINETIC_Fx2y_G3yz;
  abcd[397] = 2.0E0*I_KINETIC_Hx3yz_G3yz_a-2*I_KINETIC_Fxyz_G3yz;
  abcd[398] = 2.0E0*I_KINETIC_Hx2y2z_G3yz_a-1*I_KINETIC_Fx2z_G3yz;
  abcd[399] = 2.0E0*I_KINETIC_Hxy3z_G3yz_a;
  abcd[400] = 2.0E0*I_KINETIC_H5y_G3yz_a-4*I_KINETIC_F3y_G3yz;
  abcd[401] = 2.0E0*I_KINETIC_H4yz_G3yz_a-3*I_KINETIC_F2yz_G3yz;
  abcd[402] = 2.0E0*I_KINETIC_H3y2z_G3yz_a-2*I_KINETIC_Fy2z_G3yz;
  abcd[403] = 2.0E0*I_KINETIC_H2y3z_G3yz_a-1*I_KINETIC_F3z_G3yz;
  abcd[404] = 2.0E0*I_KINETIC_Hy4z_G3yz_a;
  abcd[405] = 2.0E0*I_KINETIC_H4xy_G2y2z_a;
  abcd[406] = 2.0E0*I_KINETIC_H3x2y_G2y2z_a-1*I_KINETIC_F3x_G2y2z;
  abcd[407] = 2.0E0*I_KINETIC_H3xyz_G2y2z_a;
  abcd[408] = 2.0E0*I_KINETIC_H2x3y_G2y2z_a-2*I_KINETIC_F2xy_G2y2z;
  abcd[409] = 2.0E0*I_KINETIC_H2x2yz_G2y2z_a-1*I_KINETIC_F2xz_G2y2z;
  abcd[410] = 2.0E0*I_KINETIC_H2xy2z_G2y2z_a;
  abcd[411] = 2.0E0*I_KINETIC_Hx4y_G2y2z_a-3*I_KINETIC_Fx2y_G2y2z;
  abcd[412] = 2.0E0*I_KINETIC_Hx3yz_G2y2z_a-2*I_KINETIC_Fxyz_G2y2z;
  abcd[413] = 2.0E0*I_KINETIC_Hx2y2z_G2y2z_a-1*I_KINETIC_Fx2z_G2y2z;
  abcd[414] = 2.0E0*I_KINETIC_Hxy3z_G2y2z_a;
  abcd[415] = 2.0E0*I_KINETIC_H5y_G2y2z_a-4*I_KINETIC_F3y_G2y2z;
  abcd[416] = 2.0E0*I_KINETIC_H4yz_G2y2z_a-3*I_KINETIC_F2yz_G2y2z;
  abcd[417] = 2.0E0*I_KINETIC_H3y2z_G2y2z_a-2*I_KINETIC_Fy2z_G2y2z;
  abcd[418] = 2.0E0*I_KINETIC_H2y3z_G2y2z_a-1*I_KINETIC_F3z_G2y2z;
  abcd[419] = 2.0E0*I_KINETIC_Hy4z_G2y2z_a;
  abcd[420] = 2.0E0*I_KINETIC_H4xy_Gy3z_a;
  abcd[421] = 2.0E0*I_KINETIC_H3x2y_Gy3z_a-1*I_KINETIC_F3x_Gy3z;
  abcd[422] = 2.0E0*I_KINETIC_H3xyz_Gy3z_a;
  abcd[423] = 2.0E0*I_KINETIC_H2x3y_Gy3z_a-2*I_KINETIC_F2xy_Gy3z;
  abcd[424] = 2.0E0*I_KINETIC_H2x2yz_Gy3z_a-1*I_KINETIC_F2xz_Gy3z;
  abcd[425] = 2.0E0*I_KINETIC_H2xy2z_Gy3z_a;
  abcd[426] = 2.0E0*I_KINETIC_Hx4y_Gy3z_a-3*I_KINETIC_Fx2y_Gy3z;
  abcd[427] = 2.0E0*I_KINETIC_Hx3yz_Gy3z_a-2*I_KINETIC_Fxyz_Gy3z;
  abcd[428] = 2.0E0*I_KINETIC_Hx2y2z_Gy3z_a-1*I_KINETIC_Fx2z_Gy3z;
  abcd[429] = 2.0E0*I_KINETIC_Hxy3z_Gy3z_a;
  abcd[430] = 2.0E0*I_KINETIC_H5y_Gy3z_a-4*I_KINETIC_F3y_Gy3z;
  abcd[431] = 2.0E0*I_KINETIC_H4yz_Gy3z_a-3*I_KINETIC_F2yz_Gy3z;
  abcd[432] = 2.0E0*I_KINETIC_H3y2z_Gy3z_a-2*I_KINETIC_Fy2z_Gy3z;
  abcd[433] = 2.0E0*I_KINETIC_H2y3z_Gy3z_a-1*I_KINETIC_F3z_Gy3z;
  abcd[434] = 2.0E0*I_KINETIC_Hy4z_Gy3z_a;
  abcd[435] = 2.0E0*I_KINETIC_H4xy_G4z_a;
  abcd[436] = 2.0E0*I_KINETIC_H3x2y_G4z_a-1*I_KINETIC_F3x_G4z;
  abcd[437] = 2.0E0*I_KINETIC_H3xyz_G4z_a;
  abcd[438] = 2.0E0*I_KINETIC_H2x3y_G4z_a-2*I_KINETIC_F2xy_G4z;
  abcd[439] = 2.0E0*I_KINETIC_H2x2yz_G4z_a-1*I_KINETIC_F2xz_G4z;
  abcd[440] = 2.0E0*I_KINETIC_H2xy2z_G4z_a;
  abcd[441] = 2.0E0*I_KINETIC_Hx4y_G4z_a-3*I_KINETIC_Fx2y_G4z;
  abcd[442] = 2.0E0*I_KINETIC_Hx3yz_G4z_a-2*I_KINETIC_Fxyz_G4z;
  abcd[443] = 2.0E0*I_KINETIC_Hx2y2z_G4z_a-1*I_KINETIC_Fx2z_G4z;
  abcd[444] = 2.0E0*I_KINETIC_Hxy3z_G4z_a;
  abcd[445] = 2.0E0*I_KINETIC_H5y_G4z_a-4*I_KINETIC_F3y_G4z;
  abcd[446] = 2.0E0*I_KINETIC_H4yz_G4z_a-3*I_KINETIC_F2yz_G4z;
  abcd[447] = 2.0E0*I_KINETIC_H3y2z_G4z_a-2*I_KINETIC_Fy2z_G4z;
  abcd[448] = 2.0E0*I_KINETIC_H2y3z_G4z_a-1*I_KINETIC_F3z_G4z;
  abcd[449] = 2.0E0*I_KINETIC_Hy4z_G4z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_G_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_G_a
   * RHS shell quartet name: SQ_KINETIC_F_G
   ************************************************************/
  abcd[450] = 2.0E0*I_KINETIC_H4xz_G4x_a;
  abcd[451] = 2.0E0*I_KINETIC_H3xyz_G4x_a;
  abcd[452] = 2.0E0*I_KINETIC_H3x2z_G4x_a-1*I_KINETIC_F3x_G4x;
  abcd[453] = 2.0E0*I_KINETIC_H2x2yz_G4x_a;
  abcd[454] = 2.0E0*I_KINETIC_H2xy2z_G4x_a-1*I_KINETIC_F2xy_G4x;
  abcd[455] = 2.0E0*I_KINETIC_H2x3z_G4x_a-2*I_KINETIC_F2xz_G4x;
  abcd[456] = 2.0E0*I_KINETIC_Hx3yz_G4x_a;
  abcd[457] = 2.0E0*I_KINETIC_Hx2y2z_G4x_a-1*I_KINETIC_Fx2y_G4x;
  abcd[458] = 2.0E0*I_KINETIC_Hxy3z_G4x_a-2*I_KINETIC_Fxyz_G4x;
  abcd[459] = 2.0E0*I_KINETIC_Hx4z_G4x_a-3*I_KINETIC_Fx2z_G4x;
  abcd[460] = 2.0E0*I_KINETIC_H4yz_G4x_a;
  abcd[461] = 2.0E0*I_KINETIC_H3y2z_G4x_a-1*I_KINETIC_F3y_G4x;
  abcd[462] = 2.0E0*I_KINETIC_H2y3z_G4x_a-2*I_KINETIC_F2yz_G4x;
  abcd[463] = 2.0E0*I_KINETIC_Hy4z_G4x_a-3*I_KINETIC_Fy2z_G4x;
  abcd[464] = 2.0E0*I_KINETIC_H5z_G4x_a-4*I_KINETIC_F3z_G4x;
  abcd[465] = 2.0E0*I_KINETIC_H4xz_G3xy_a;
  abcd[466] = 2.0E0*I_KINETIC_H3xyz_G3xy_a;
  abcd[467] = 2.0E0*I_KINETIC_H3x2z_G3xy_a-1*I_KINETIC_F3x_G3xy;
  abcd[468] = 2.0E0*I_KINETIC_H2x2yz_G3xy_a;
  abcd[469] = 2.0E0*I_KINETIC_H2xy2z_G3xy_a-1*I_KINETIC_F2xy_G3xy;
  abcd[470] = 2.0E0*I_KINETIC_H2x3z_G3xy_a-2*I_KINETIC_F2xz_G3xy;
  abcd[471] = 2.0E0*I_KINETIC_Hx3yz_G3xy_a;
  abcd[472] = 2.0E0*I_KINETIC_Hx2y2z_G3xy_a-1*I_KINETIC_Fx2y_G3xy;
  abcd[473] = 2.0E0*I_KINETIC_Hxy3z_G3xy_a-2*I_KINETIC_Fxyz_G3xy;
  abcd[474] = 2.0E0*I_KINETIC_Hx4z_G3xy_a-3*I_KINETIC_Fx2z_G3xy;
  abcd[475] = 2.0E0*I_KINETIC_H4yz_G3xy_a;
  abcd[476] = 2.0E0*I_KINETIC_H3y2z_G3xy_a-1*I_KINETIC_F3y_G3xy;
  abcd[477] = 2.0E0*I_KINETIC_H2y3z_G3xy_a-2*I_KINETIC_F2yz_G3xy;
  abcd[478] = 2.0E0*I_KINETIC_Hy4z_G3xy_a-3*I_KINETIC_Fy2z_G3xy;
  abcd[479] = 2.0E0*I_KINETIC_H5z_G3xy_a-4*I_KINETIC_F3z_G3xy;
  abcd[480] = 2.0E0*I_KINETIC_H4xz_G3xz_a;
  abcd[481] = 2.0E0*I_KINETIC_H3xyz_G3xz_a;
  abcd[482] = 2.0E0*I_KINETIC_H3x2z_G3xz_a-1*I_KINETIC_F3x_G3xz;
  abcd[483] = 2.0E0*I_KINETIC_H2x2yz_G3xz_a;
  abcd[484] = 2.0E0*I_KINETIC_H2xy2z_G3xz_a-1*I_KINETIC_F2xy_G3xz;
  abcd[485] = 2.0E0*I_KINETIC_H2x3z_G3xz_a-2*I_KINETIC_F2xz_G3xz;
  abcd[486] = 2.0E0*I_KINETIC_Hx3yz_G3xz_a;
  abcd[487] = 2.0E0*I_KINETIC_Hx2y2z_G3xz_a-1*I_KINETIC_Fx2y_G3xz;
  abcd[488] = 2.0E0*I_KINETIC_Hxy3z_G3xz_a-2*I_KINETIC_Fxyz_G3xz;
  abcd[489] = 2.0E0*I_KINETIC_Hx4z_G3xz_a-3*I_KINETIC_Fx2z_G3xz;
  abcd[490] = 2.0E0*I_KINETIC_H4yz_G3xz_a;
  abcd[491] = 2.0E0*I_KINETIC_H3y2z_G3xz_a-1*I_KINETIC_F3y_G3xz;
  abcd[492] = 2.0E0*I_KINETIC_H2y3z_G3xz_a-2*I_KINETIC_F2yz_G3xz;
  abcd[493] = 2.0E0*I_KINETIC_Hy4z_G3xz_a-3*I_KINETIC_Fy2z_G3xz;
  abcd[494] = 2.0E0*I_KINETIC_H5z_G3xz_a-4*I_KINETIC_F3z_G3xz;
  abcd[495] = 2.0E0*I_KINETIC_H4xz_G2x2y_a;
  abcd[496] = 2.0E0*I_KINETIC_H3xyz_G2x2y_a;
  abcd[497] = 2.0E0*I_KINETIC_H3x2z_G2x2y_a-1*I_KINETIC_F3x_G2x2y;
  abcd[498] = 2.0E0*I_KINETIC_H2x2yz_G2x2y_a;
  abcd[499] = 2.0E0*I_KINETIC_H2xy2z_G2x2y_a-1*I_KINETIC_F2xy_G2x2y;
  abcd[500] = 2.0E0*I_KINETIC_H2x3z_G2x2y_a-2*I_KINETIC_F2xz_G2x2y;
  abcd[501] = 2.0E0*I_KINETIC_Hx3yz_G2x2y_a;
  abcd[502] = 2.0E0*I_KINETIC_Hx2y2z_G2x2y_a-1*I_KINETIC_Fx2y_G2x2y;
  abcd[503] = 2.0E0*I_KINETIC_Hxy3z_G2x2y_a-2*I_KINETIC_Fxyz_G2x2y;
  abcd[504] = 2.0E0*I_KINETIC_Hx4z_G2x2y_a-3*I_KINETIC_Fx2z_G2x2y;
  abcd[505] = 2.0E0*I_KINETIC_H4yz_G2x2y_a;
  abcd[506] = 2.0E0*I_KINETIC_H3y2z_G2x2y_a-1*I_KINETIC_F3y_G2x2y;
  abcd[507] = 2.0E0*I_KINETIC_H2y3z_G2x2y_a-2*I_KINETIC_F2yz_G2x2y;
  abcd[508] = 2.0E0*I_KINETIC_Hy4z_G2x2y_a-3*I_KINETIC_Fy2z_G2x2y;
  abcd[509] = 2.0E0*I_KINETIC_H5z_G2x2y_a-4*I_KINETIC_F3z_G2x2y;
  abcd[510] = 2.0E0*I_KINETIC_H4xz_G2xyz_a;
  abcd[511] = 2.0E0*I_KINETIC_H3xyz_G2xyz_a;
  abcd[512] = 2.0E0*I_KINETIC_H3x2z_G2xyz_a-1*I_KINETIC_F3x_G2xyz;
  abcd[513] = 2.0E0*I_KINETIC_H2x2yz_G2xyz_a;
  abcd[514] = 2.0E0*I_KINETIC_H2xy2z_G2xyz_a-1*I_KINETIC_F2xy_G2xyz;
  abcd[515] = 2.0E0*I_KINETIC_H2x3z_G2xyz_a-2*I_KINETIC_F2xz_G2xyz;
  abcd[516] = 2.0E0*I_KINETIC_Hx3yz_G2xyz_a;
  abcd[517] = 2.0E0*I_KINETIC_Hx2y2z_G2xyz_a-1*I_KINETIC_Fx2y_G2xyz;
  abcd[518] = 2.0E0*I_KINETIC_Hxy3z_G2xyz_a-2*I_KINETIC_Fxyz_G2xyz;
  abcd[519] = 2.0E0*I_KINETIC_Hx4z_G2xyz_a-3*I_KINETIC_Fx2z_G2xyz;
  abcd[520] = 2.0E0*I_KINETIC_H4yz_G2xyz_a;
  abcd[521] = 2.0E0*I_KINETIC_H3y2z_G2xyz_a-1*I_KINETIC_F3y_G2xyz;
  abcd[522] = 2.0E0*I_KINETIC_H2y3z_G2xyz_a-2*I_KINETIC_F2yz_G2xyz;
  abcd[523] = 2.0E0*I_KINETIC_Hy4z_G2xyz_a-3*I_KINETIC_Fy2z_G2xyz;
  abcd[524] = 2.0E0*I_KINETIC_H5z_G2xyz_a-4*I_KINETIC_F3z_G2xyz;
  abcd[525] = 2.0E0*I_KINETIC_H4xz_G2x2z_a;
  abcd[526] = 2.0E0*I_KINETIC_H3xyz_G2x2z_a;
  abcd[527] = 2.0E0*I_KINETIC_H3x2z_G2x2z_a-1*I_KINETIC_F3x_G2x2z;
  abcd[528] = 2.0E0*I_KINETIC_H2x2yz_G2x2z_a;
  abcd[529] = 2.0E0*I_KINETIC_H2xy2z_G2x2z_a-1*I_KINETIC_F2xy_G2x2z;
  abcd[530] = 2.0E0*I_KINETIC_H2x3z_G2x2z_a-2*I_KINETIC_F2xz_G2x2z;
  abcd[531] = 2.0E0*I_KINETIC_Hx3yz_G2x2z_a;
  abcd[532] = 2.0E0*I_KINETIC_Hx2y2z_G2x2z_a-1*I_KINETIC_Fx2y_G2x2z;
  abcd[533] = 2.0E0*I_KINETIC_Hxy3z_G2x2z_a-2*I_KINETIC_Fxyz_G2x2z;
  abcd[534] = 2.0E0*I_KINETIC_Hx4z_G2x2z_a-3*I_KINETIC_Fx2z_G2x2z;
  abcd[535] = 2.0E0*I_KINETIC_H4yz_G2x2z_a;
  abcd[536] = 2.0E0*I_KINETIC_H3y2z_G2x2z_a-1*I_KINETIC_F3y_G2x2z;
  abcd[537] = 2.0E0*I_KINETIC_H2y3z_G2x2z_a-2*I_KINETIC_F2yz_G2x2z;
  abcd[538] = 2.0E0*I_KINETIC_Hy4z_G2x2z_a-3*I_KINETIC_Fy2z_G2x2z;
  abcd[539] = 2.0E0*I_KINETIC_H5z_G2x2z_a-4*I_KINETIC_F3z_G2x2z;
  abcd[540] = 2.0E0*I_KINETIC_H4xz_Gx3y_a;
  abcd[541] = 2.0E0*I_KINETIC_H3xyz_Gx3y_a;
  abcd[542] = 2.0E0*I_KINETIC_H3x2z_Gx3y_a-1*I_KINETIC_F3x_Gx3y;
  abcd[543] = 2.0E0*I_KINETIC_H2x2yz_Gx3y_a;
  abcd[544] = 2.0E0*I_KINETIC_H2xy2z_Gx3y_a-1*I_KINETIC_F2xy_Gx3y;
  abcd[545] = 2.0E0*I_KINETIC_H2x3z_Gx3y_a-2*I_KINETIC_F2xz_Gx3y;
  abcd[546] = 2.0E0*I_KINETIC_Hx3yz_Gx3y_a;
  abcd[547] = 2.0E0*I_KINETIC_Hx2y2z_Gx3y_a-1*I_KINETIC_Fx2y_Gx3y;
  abcd[548] = 2.0E0*I_KINETIC_Hxy3z_Gx3y_a-2*I_KINETIC_Fxyz_Gx3y;
  abcd[549] = 2.0E0*I_KINETIC_Hx4z_Gx3y_a-3*I_KINETIC_Fx2z_Gx3y;
  abcd[550] = 2.0E0*I_KINETIC_H4yz_Gx3y_a;
  abcd[551] = 2.0E0*I_KINETIC_H3y2z_Gx3y_a-1*I_KINETIC_F3y_Gx3y;
  abcd[552] = 2.0E0*I_KINETIC_H2y3z_Gx3y_a-2*I_KINETIC_F2yz_Gx3y;
  abcd[553] = 2.0E0*I_KINETIC_Hy4z_Gx3y_a-3*I_KINETIC_Fy2z_Gx3y;
  abcd[554] = 2.0E0*I_KINETIC_H5z_Gx3y_a-4*I_KINETIC_F3z_Gx3y;
  abcd[555] = 2.0E0*I_KINETIC_H4xz_Gx2yz_a;
  abcd[556] = 2.0E0*I_KINETIC_H3xyz_Gx2yz_a;
  abcd[557] = 2.0E0*I_KINETIC_H3x2z_Gx2yz_a-1*I_KINETIC_F3x_Gx2yz;
  abcd[558] = 2.0E0*I_KINETIC_H2x2yz_Gx2yz_a;
  abcd[559] = 2.0E0*I_KINETIC_H2xy2z_Gx2yz_a-1*I_KINETIC_F2xy_Gx2yz;
  abcd[560] = 2.0E0*I_KINETIC_H2x3z_Gx2yz_a-2*I_KINETIC_F2xz_Gx2yz;
  abcd[561] = 2.0E0*I_KINETIC_Hx3yz_Gx2yz_a;
  abcd[562] = 2.0E0*I_KINETIC_Hx2y2z_Gx2yz_a-1*I_KINETIC_Fx2y_Gx2yz;
  abcd[563] = 2.0E0*I_KINETIC_Hxy3z_Gx2yz_a-2*I_KINETIC_Fxyz_Gx2yz;
  abcd[564] = 2.0E0*I_KINETIC_Hx4z_Gx2yz_a-3*I_KINETIC_Fx2z_Gx2yz;
  abcd[565] = 2.0E0*I_KINETIC_H4yz_Gx2yz_a;
  abcd[566] = 2.0E0*I_KINETIC_H3y2z_Gx2yz_a-1*I_KINETIC_F3y_Gx2yz;
  abcd[567] = 2.0E0*I_KINETIC_H2y3z_Gx2yz_a-2*I_KINETIC_F2yz_Gx2yz;
  abcd[568] = 2.0E0*I_KINETIC_Hy4z_Gx2yz_a-3*I_KINETIC_Fy2z_Gx2yz;
  abcd[569] = 2.0E0*I_KINETIC_H5z_Gx2yz_a-4*I_KINETIC_F3z_Gx2yz;
  abcd[570] = 2.0E0*I_KINETIC_H4xz_Gxy2z_a;
  abcd[571] = 2.0E0*I_KINETIC_H3xyz_Gxy2z_a;
  abcd[572] = 2.0E0*I_KINETIC_H3x2z_Gxy2z_a-1*I_KINETIC_F3x_Gxy2z;
  abcd[573] = 2.0E0*I_KINETIC_H2x2yz_Gxy2z_a;
  abcd[574] = 2.0E0*I_KINETIC_H2xy2z_Gxy2z_a-1*I_KINETIC_F2xy_Gxy2z;
  abcd[575] = 2.0E0*I_KINETIC_H2x3z_Gxy2z_a-2*I_KINETIC_F2xz_Gxy2z;
  abcd[576] = 2.0E0*I_KINETIC_Hx3yz_Gxy2z_a;
  abcd[577] = 2.0E0*I_KINETIC_Hx2y2z_Gxy2z_a-1*I_KINETIC_Fx2y_Gxy2z;
  abcd[578] = 2.0E0*I_KINETIC_Hxy3z_Gxy2z_a-2*I_KINETIC_Fxyz_Gxy2z;
  abcd[579] = 2.0E0*I_KINETIC_Hx4z_Gxy2z_a-3*I_KINETIC_Fx2z_Gxy2z;
  abcd[580] = 2.0E0*I_KINETIC_H4yz_Gxy2z_a;
  abcd[581] = 2.0E0*I_KINETIC_H3y2z_Gxy2z_a-1*I_KINETIC_F3y_Gxy2z;
  abcd[582] = 2.0E0*I_KINETIC_H2y3z_Gxy2z_a-2*I_KINETIC_F2yz_Gxy2z;
  abcd[583] = 2.0E0*I_KINETIC_Hy4z_Gxy2z_a-3*I_KINETIC_Fy2z_Gxy2z;
  abcd[584] = 2.0E0*I_KINETIC_H5z_Gxy2z_a-4*I_KINETIC_F3z_Gxy2z;
  abcd[585] = 2.0E0*I_KINETIC_H4xz_Gx3z_a;
  abcd[586] = 2.0E0*I_KINETIC_H3xyz_Gx3z_a;
  abcd[587] = 2.0E0*I_KINETIC_H3x2z_Gx3z_a-1*I_KINETIC_F3x_Gx3z;
  abcd[588] = 2.0E0*I_KINETIC_H2x2yz_Gx3z_a;
  abcd[589] = 2.0E0*I_KINETIC_H2xy2z_Gx3z_a-1*I_KINETIC_F2xy_Gx3z;
  abcd[590] = 2.0E0*I_KINETIC_H2x3z_Gx3z_a-2*I_KINETIC_F2xz_Gx3z;
  abcd[591] = 2.0E0*I_KINETIC_Hx3yz_Gx3z_a;
  abcd[592] = 2.0E0*I_KINETIC_Hx2y2z_Gx3z_a-1*I_KINETIC_Fx2y_Gx3z;
  abcd[593] = 2.0E0*I_KINETIC_Hxy3z_Gx3z_a-2*I_KINETIC_Fxyz_Gx3z;
  abcd[594] = 2.0E0*I_KINETIC_Hx4z_Gx3z_a-3*I_KINETIC_Fx2z_Gx3z;
  abcd[595] = 2.0E0*I_KINETIC_H4yz_Gx3z_a;
  abcd[596] = 2.0E0*I_KINETIC_H3y2z_Gx3z_a-1*I_KINETIC_F3y_Gx3z;
  abcd[597] = 2.0E0*I_KINETIC_H2y3z_Gx3z_a-2*I_KINETIC_F2yz_Gx3z;
  abcd[598] = 2.0E0*I_KINETIC_Hy4z_Gx3z_a-3*I_KINETIC_Fy2z_Gx3z;
  abcd[599] = 2.0E0*I_KINETIC_H5z_Gx3z_a-4*I_KINETIC_F3z_Gx3z;
  abcd[600] = 2.0E0*I_KINETIC_H4xz_G4y_a;
  abcd[601] = 2.0E0*I_KINETIC_H3xyz_G4y_a;
  abcd[602] = 2.0E0*I_KINETIC_H3x2z_G4y_a-1*I_KINETIC_F3x_G4y;
  abcd[603] = 2.0E0*I_KINETIC_H2x2yz_G4y_a;
  abcd[604] = 2.0E0*I_KINETIC_H2xy2z_G4y_a-1*I_KINETIC_F2xy_G4y;
  abcd[605] = 2.0E0*I_KINETIC_H2x3z_G4y_a-2*I_KINETIC_F2xz_G4y;
  abcd[606] = 2.0E0*I_KINETIC_Hx3yz_G4y_a;
  abcd[607] = 2.0E0*I_KINETIC_Hx2y2z_G4y_a-1*I_KINETIC_Fx2y_G4y;
  abcd[608] = 2.0E0*I_KINETIC_Hxy3z_G4y_a-2*I_KINETIC_Fxyz_G4y;
  abcd[609] = 2.0E0*I_KINETIC_Hx4z_G4y_a-3*I_KINETIC_Fx2z_G4y;
  abcd[610] = 2.0E0*I_KINETIC_H4yz_G4y_a;
  abcd[611] = 2.0E0*I_KINETIC_H3y2z_G4y_a-1*I_KINETIC_F3y_G4y;
  abcd[612] = 2.0E0*I_KINETIC_H2y3z_G4y_a-2*I_KINETIC_F2yz_G4y;
  abcd[613] = 2.0E0*I_KINETIC_Hy4z_G4y_a-3*I_KINETIC_Fy2z_G4y;
  abcd[614] = 2.0E0*I_KINETIC_H5z_G4y_a-4*I_KINETIC_F3z_G4y;
  abcd[615] = 2.0E0*I_KINETIC_H4xz_G3yz_a;
  abcd[616] = 2.0E0*I_KINETIC_H3xyz_G3yz_a;
  abcd[617] = 2.0E0*I_KINETIC_H3x2z_G3yz_a-1*I_KINETIC_F3x_G3yz;
  abcd[618] = 2.0E0*I_KINETIC_H2x2yz_G3yz_a;
  abcd[619] = 2.0E0*I_KINETIC_H2xy2z_G3yz_a-1*I_KINETIC_F2xy_G3yz;
  abcd[620] = 2.0E0*I_KINETIC_H2x3z_G3yz_a-2*I_KINETIC_F2xz_G3yz;
  abcd[621] = 2.0E0*I_KINETIC_Hx3yz_G3yz_a;
  abcd[622] = 2.0E0*I_KINETIC_Hx2y2z_G3yz_a-1*I_KINETIC_Fx2y_G3yz;
  abcd[623] = 2.0E0*I_KINETIC_Hxy3z_G3yz_a-2*I_KINETIC_Fxyz_G3yz;
  abcd[624] = 2.0E0*I_KINETIC_Hx4z_G3yz_a-3*I_KINETIC_Fx2z_G3yz;
  abcd[625] = 2.0E0*I_KINETIC_H4yz_G3yz_a;
  abcd[626] = 2.0E0*I_KINETIC_H3y2z_G3yz_a-1*I_KINETIC_F3y_G3yz;
  abcd[627] = 2.0E0*I_KINETIC_H2y3z_G3yz_a-2*I_KINETIC_F2yz_G3yz;
  abcd[628] = 2.0E0*I_KINETIC_Hy4z_G3yz_a-3*I_KINETIC_Fy2z_G3yz;
  abcd[629] = 2.0E0*I_KINETIC_H5z_G3yz_a-4*I_KINETIC_F3z_G3yz;
  abcd[630] = 2.0E0*I_KINETIC_H4xz_G2y2z_a;
  abcd[631] = 2.0E0*I_KINETIC_H3xyz_G2y2z_a;
  abcd[632] = 2.0E0*I_KINETIC_H3x2z_G2y2z_a-1*I_KINETIC_F3x_G2y2z;
  abcd[633] = 2.0E0*I_KINETIC_H2x2yz_G2y2z_a;
  abcd[634] = 2.0E0*I_KINETIC_H2xy2z_G2y2z_a-1*I_KINETIC_F2xy_G2y2z;
  abcd[635] = 2.0E0*I_KINETIC_H2x3z_G2y2z_a-2*I_KINETIC_F2xz_G2y2z;
  abcd[636] = 2.0E0*I_KINETIC_Hx3yz_G2y2z_a;
  abcd[637] = 2.0E0*I_KINETIC_Hx2y2z_G2y2z_a-1*I_KINETIC_Fx2y_G2y2z;
  abcd[638] = 2.0E0*I_KINETIC_Hxy3z_G2y2z_a-2*I_KINETIC_Fxyz_G2y2z;
  abcd[639] = 2.0E0*I_KINETIC_Hx4z_G2y2z_a-3*I_KINETIC_Fx2z_G2y2z;
  abcd[640] = 2.0E0*I_KINETIC_H4yz_G2y2z_a;
  abcd[641] = 2.0E0*I_KINETIC_H3y2z_G2y2z_a-1*I_KINETIC_F3y_G2y2z;
  abcd[642] = 2.0E0*I_KINETIC_H2y3z_G2y2z_a-2*I_KINETIC_F2yz_G2y2z;
  abcd[643] = 2.0E0*I_KINETIC_Hy4z_G2y2z_a-3*I_KINETIC_Fy2z_G2y2z;
  abcd[644] = 2.0E0*I_KINETIC_H5z_G2y2z_a-4*I_KINETIC_F3z_G2y2z;
  abcd[645] = 2.0E0*I_KINETIC_H4xz_Gy3z_a;
  abcd[646] = 2.0E0*I_KINETIC_H3xyz_Gy3z_a;
  abcd[647] = 2.0E0*I_KINETIC_H3x2z_Gy3z_a-1*I_KINETIC_F3x_Gy3z;
  abcd[648] = 2.0E0*I_KINETIC_H2x2yz_Gy3z_a;
  abcd[649] = 2.0E0*I_KINETIC_H2xy2z_Gy3z_a-1*I_KINETIC_F2xy_Gy3z;
  abcd[650] = 2.0E0*I_KINETIC_H2x3z_Gy3z_a-2*I_KINETIC_F2xz_Gy3z;
  abcd[651] = 2.0E0*I_KINETIC_Hx3yz_Gy3z_a;
  abcd[652] = 2.0E0*I_KINETIC_Hx2y2z_Gy3z_a-1*I_KINETIC_Fx2y_Gy3z;
  abcd[653] = 2.0E0*I_KINETIC_Hxy3z_Gy3z_a-2*I_KINETIC_Fxyz_Gy3z;
  abcd[654] = 2.0E0*I_KINETIC_Hx4z_Gy3z_a-3*I_KINETIC_Fx2z_Gy3z;
  abcd[655] = 2.0E0*I_KINETIC_H4yz_Gy3z_a;
  abcd[656] = 2.0E0*I_KINETIC_H3y2z_Gy3z_a-1*I_KINETIC_F3y_Gy3z;
  abcd[657] = 2.0E0*I_KINETIC_H2y3z_Gy3z_a-2*I_KINETIC_F2yz_Gy3z;
  abcd[658] = 2.0E0*I_KINETIC_Hy4z_Gy3z_a-3*I_KINETIC_Fy2z_Gy3z;
  abcd[659] = 2.0E0*I_KINETIC_H5z_Gy3z_a-4*I_KINETIC_F3z_Gy3z;
  abcd[660] = 2.0E0*I_KINETIC_H4xz_G4z_a;
  abcd[661] = 2.0E0*I_KINETIC_H3xyz_G4z_a;
  abcd[662] = 2.0E0*I_KINETIC_H3x2z_G4z_a-1*I_KINETIC_F3x_G4z;
  abcd[663] = 2.0E0*I_KINETIC_H2x2yz_G4z_a;
  abcd[664] = 2.0E0*I_KINETIC_H2xy2z_G4z_a-1*I_KINETIC_F2xy_G4z;
  abcd[665] = 2.0E0*I_KINETIC_H2x3z_G4z_a-2*I_KINETIC_F2xz_G4z;
  abcd[666] = 2.0E0*I_KINETIC_Hx3yz_G4z_a;
  abcd[667] = 2.0E0*I_KINETIC_Hx2y2z_G4z_a-1*I_KINETIC_Fx2y_G4z;
  abcd[668] = 2.0E0*I_KINETIC_Hxy3z_G4z_a-2*I_KINETIC_Fxyz_G4z;
  abcd[669] = 2.0E0*I_KINETIC_Hx4z_G4z_a-3*I_KINETIC_Fx2z_G4z;
  abcd[670] = 2.0E0*I_KINETIC_H4yz_G4z_a;
  abcd[671] = 2.0E0*I_KINETIC_H3y2z_G4z_a-1*I_KINETIC_F3y_G4z;
  abcd[672] = 2.0E0*I_KINETIC_H2y3z_G4z_a-2*I_KINETIC_F2yz_G4z;
  abcd[673] = 2.0E0*I_KINETIC_Hy4z_G4z_a-3*I_KINETIC_Fy2z_G4z;
  abcd[674] = 2.0E0*I_KINETIC_H5z_G4z_a-4*I_KINETIC_F3z_G4z;
}
