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
// BRA1 as redundant position, total RHS integrals evaluated as: 14198
// BRA2 as redundant position, total RHS integrals evaluated as: 14094
// KET1 as redundant position, total RHS integrals evaluated as: 10737
// KET2 as redundant position, total RHS integrals evaluated as: 10737
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

void hgp_os_eri_f_d_s_s_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_G4x_S_S_S = 0.0E0;
  Double I_ERI_G3xy_S_S_S = 0.0E0;
  Double I_ERI_G3xz_S_S_S = 0.0E0;
  Double I_ERI_G2x2y_S_S_S = 0.0E0;
  Double I_ERI_G2xyz_S_S_S = 0.0E0;
  Double I_ERI_G2x2z_S_S_S = 0.0E0;
  Double I_ERI_Gx3y_S_S_S = 0.0E0;
  Double I_ERI_Gx2yz_S_S_S = 0.0E0;
  Double I_ERI_Gxy2z_S_S_S = 0.0E0;
  Double I_ERI_Gx3z_S_S_S = 0.0E0;
  Double I_ERI_G4y_S_S_S = 0.0E0;
  Double I_ERI_G3yz_S_S_S = 0.0E0;
  Double I_ERI_G2y2z_S_S_S = 0.0E0;
  Double I_ERI_Gy3z_S_S_S = 0.0E0;
  Double I_ERI_G4z_S_S_S = 0.0E0;
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
  Double I_ERI_I6x_S_S_S_a = 0.0E0;
  Double I_ERI_I5xy_S_S_S_a = 0.0E0;
  Double I_ERI_I5xz_S_S_S_a = 0.0E0;
  Double I_ERI_I4x2y_S_S_S_a = 0.0E0;
  Double I_ERI_I4xyz_S_S_S_a = 0.0E0;
  Double I_ERI_I4x2z_S_S_S_a = 0.0E0;
  Double I_ERI_I3x3y_S_S_S_a = 0.0E0;
  Double I_ERI_I3x2yz_S_S_S_a = 0.0E0;
  Double I_ERI_I3xy2z_S_S_S_a = 0.0E0;
  Double I_ERI_I3x3z_S_S_S_a = 0.0E0;
  Double I_ERI_I2x4y_S_S_S_a = 0.0E0;
  Double I_ERI_I2x3yz_S_S_S_a = 0.0E0;
  Double I_ERI_I2x2y2z_S_S_S_a = 0.0E0;
  Double I_ERI_I2xy3z_S_S_S_a = 0.0E0;
  Double I_ERI_I2x4z_S_S_S_a = 0.0E0;
  Double I_ERI_Ix5y_S_S_S_a = 0.0E0;
  Double I_ERI_Ix4yz_S_S_S_a = 0.0E0;
  Double I_ERI_Ix3y2z_S_S_S_a = 0.0E0;
  Double I_ERI_Ix2y3z_S_S_S_a = 0.0E0;
  Double I_ERI_Ixy4z_S_S_S_a = 0.0E0;
  Double I_ERI_Ix5z_S_S_S_a = 0.0E0;
  Double I_ERI_I6y_S_S_S_a = 0.0E0;
  Double I_ERI_I5yz_S_S_S_a = 0.0E0;
  Double I_ERI_I4y2z_S_S_S_a = 0.0E0;
  Double I_ERI_I3y3z_S_S_S_a = 0.0E0;
  Double I_ERI_I2y4z_S_S_S_a = 0.0E0;
  Double I_ERI_Iy5z_S_S_S_a = 0.0E0;
  Double I_ERI_I6z_S_S_S_a = 0.0E0;
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
  Double I_ERI_H5x_S_Px_S_c = 0.0E0;
  Double I_ERI_H4xy_S_Px_S_c = 0.0E0;
  Double I_ERI_H4xz_S_Px_S_c = 0.0E0;
  Double I_ERI_H3x2y_S_Px_S_c = 0.0E0;
  Double I_ERI_H3xyz_S_Px_S_c = 0.0E0;
  Double I_ERI_H3x2z_S_Px_S_c = 0.0E0;
  Double I_ERI_H2x3y_S_Px_S_c = 0.0E0;
  Double I_ERI_H2x2yz_S_Px_S_c = 0.0E0;
  Double I_ERI_H2xy2z_S_Px_S_c = 0.0E0;
  Double I_ERI_H2x3z_S_Px_S_c = 0.0E0;
  Double I_ERI_Hx4y_S_Px_S_c = 0.0E0;
  Double I_ERI_Hx3yz_S_Px_S_c = 0.0E0;
  Double I_ERI_Hx2y2z_S_Px_S_c = 0.0E0;
  Double I_ERI_Hxy3z_S_Px_S_c = 0.0E0;
  Double I_ERI_Hx4z_S_Px_S_c = 0.0E0;
  Double I_ERI_H5y_S_Px_S_c = 0.0E0;
  Double I_ERI_H4yz_S_Px_S_c = 0.0E0;
  Double I_ERI_H3y2z_S_Px_S_c = 0.0E0;
  Double I_ERI_H2y3z_S_Px_S_c = 0.0E0;
  Double I_ERI_Hy4z_S_Px_S_c = 0.0E0;
  Double I_ERI_H5z_S_Px_S_c = 0.0E0;
  Double I_ERI_H5x_S_Py_S_c = 0.0E0;
  Double I_ERI_H4xy_S_Py_S_c = 0.0E0;
  Double I_ERI_H4xz_S_Py_S_c = 0.0E0;
  Double I_ERI_H3x2y_S_Py_S_c = 0.0E0;
  Double I_ERI_H3xyz_S_Py_S_c = 0.0E0;
  Double I_ERI_H3x2z_S_Py_S_c = 0.0E0;
  Double I_ERI_H2x3y_S_Py_S_c = 0.0E0;
  Double I_ERI_H2x2yz_S_Py_S_c = 0.0E0;
  Double I_ERI_H2xy2z_S_Py_S_c = 0.0E0;
  Double I_ERI_H2x3z_S_Py_S_c = 0.0E0;
  Double I_ERI_Hx4y_S_Py_S_c = 0.0E0;
  Double I_ERI_Hx3yz_S_Py_S_c = 0.0E0;
  Double I_ERI_Hx2y2z_S_Py_S_c = 0.0E0;
  Double I_ERI_Hxy3z_S_Py_S_c = 0.0E0;
  Double I_ERI_Hx4z_S_Py_S_c = 0.0E0;
  Double I_ERI_H5y_S_Py_S_c = 0.0E0;
  Double I_ERI_H4yz_S_Py_S_c = 0.0E0;
  Double I_ERI_H3y2z_S_Py_S_c = 0.0E0;
  Double I_ERI_H2y3z_S_Py_S_c = 0.0E0;
  Double I_ERI_Hy4z_S_Py_S_c = 0.0E0;
  Double I_ERI_H5z_S_Py_S_c = 0.0E0;
  Double I_ERI_H5x_S_Pz_S_c = 0.0E0;
  Double I_ERI_H4xy_S_Pz_S_c = 0.0E0;
  Double I_ERI_H4xz_S_Pz_S_c = 0.0E0;
  Double I_ERI_H3x2y_S_Pz_S_c = 0.0E0;
  Double I_ERI_H3xyz_S_Pz_S_c = 0.0E0;
  Double I_ERI_H3x2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_H2x3y_S_Pz_S_c = 0.0E0;
  Double I_ERI_H2x2yz_S_Pz_S_c = 0.0E0;
  Double I_ERI_H2xy2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_H2x3z_S_Pz_S_c = 0.0E0;
  Double I_ERI_Hx4y_S_Pz_S_c = 0.0E0;
  Double I_ERI_Hx3yz_S_Pz_S_c = 0.0E0;
  Double I_ERI_Hx2y2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_Hxy3z_S_Pz_S_c = 0.0E0;
  Double I_ERI_Hx4z_S_Pz_S_c = 0.0E0;
  Double I_ERI_H5y_S_Pz_S_c = 0.0E0;
  Double I_ERI_H4yz_S_Pz_S_c = 0.0E0;
  Double I_ERI_H3y2z_S_Pz_S_c = 0.0E0;
  Double I_ERI_H2y3z_S_Pz_S_c = 0.0E0;
  Double I_ERI_Hy4z_S_Pz_S_c = 0.0E0;
  Double I_ERI_H5z_S_Pz_S_c = 0.0E0;
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
  Double I_ERI_I6x_S_S_S_b = 0.0E0;
  Double I_ERI_I5xy_S_S_S_b = 0.0E0;
  Double I_ERI_I5xz_S_S_S_b = 0.0E0;
  Double I_ERI_I4x2y_S_S_S_b = 0.0E0;
  Double I_ERI_I4xyz_S_S_S_b = 0.0E0;
  Double I_ERI_I4x2z_S_S_S_b = 0.0E0;
  Double I_ERI_I3x3y_S_S_S_b = 0.0E0;
  Double I_ERI_I3x2yz_S_S_S_b = 0.0E0;
  Double I_ERI_I3xy2z_S_S_S_b = 0.0E0;
  Double I_ERI_I3x3z_S_S_S_b = 0.0E0;
  Double I_ERI_I2x4y_S_S_S_b = 0.0E0;
  Double I_ERI_I2x3yz_S_S_S_b = 0.0E0;
  Double I_ERI_I2x2y2z_S_S_S_b = 0.0E0;
  Double I_ERI_I2xy3z_S_S_S_b = 0.0E0;
  Double I_ERI_I2x4z_S_S_S_b = 0.0E0;
  Double I_ERI_Ix5y_S_S_S_b = 0.0E0;
  Double I_ERI_Ix4yz_S_S_S_b = 0.0E0;
  Double I_ERI_Ix3y2z_S_S_S_b = 0.0E0;
  Double I_ERI_Ix2y3z_S_S_S_b = 0.0E0;
  Double I_ERI_Ixy4z_S_S_S_b = 0.0E0;
  Double I_ERI_Ix5z_S_S_S_b = 0.0E0;
  Double I_ERI_I6y_S_S_S_b = 0.0E0;
  Double I_ERI_I5yz_S_S_S_b = 0.0E0;
  Double I_ERI_I4y2z_S_S_S_b = 0.0E0;
  Double I_ERI_I3y3z_S_S_S_b = 0.0E0;
  Double I_ERI_I2y4z_S_S_S_b = 0.0E0;
  Double I_ERI_Iy5z_S_S_S_b = 0.0E0;
  Double I_ERI_I6z_S_S_S_b = 0.0E0;
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
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M3_vrr = PAX*I_ERI_Px_S_S_S_M3_vrr+WPX*I_ERI_Px_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_D2y_S_S_S_M3_vrr = PAY*I_ERI_Py_S_S_S_M3_vrr+WPY*I_ERI_Py_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_D2z_S_S_S_M3_vrr = PAZ*I_ERI_Pz_S_S_S_M3_vrr+WPZ*I_ERI_Pz_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M3_vrr = PAX*I_ERI_D2x_S_S_S_M3_vrr+WPX*I_ERI_D2x_S_S_S_M4_vrr+2*oned2z*I_ERI_Px_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_F2xy_S_S_S_M3_vrr = PAY*I_ERI_D2x_S_S_S_M3_vrr+WPY*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F2xz_S_S_S_M3_vrr = PAZ*I_ERI_D2x_S_S_S_M3_vrr+WPZ*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F3y_S_S_S_M3_vrr = PAY*I_ERI_D2y_S_S_S_M3_vrr+WPY*I_ERI_D2y_S_S_S_M4_vrr+2*oned2z*I_ERI_Py_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_F2yz_S_S_S_M3_vrr = PAZ*I_ERI_D2y_S_S_S_M3_vrr+WPZ*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_F3z_S_S_S_M3_vrr = PAZ*I_ERI_D2z_S_S_S_M3_vrr+WPZ*I_ERI_D2z_S_S_S_M4_vrr+2*oned2z*I_ERI_Pz_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M4_vrr;

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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M2_vrr = PAX*I_ERI_Px_S_S_S_M2_vrr+WPX*I_ERI_Px_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_Dxy_S_S_S_M2_vrr = PAY*I_ERI_Px_S_S_S_M2_vrr+WPY*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_S_S_M2_vrr = PAY*I_ERI_Py_S_S_S_M2_vrr+WPY*I_ERI_Py_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_S_S_M2_vrr = PAZ*I_ERI_Pz_S_S_S_M2_vrr+WPZ*I_ERI_Pz_S_S_S_M3_vrr+oned2z*I_ERI_S_S_S_S_M2_vrr-rhod2zsq*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M2_vrr = PAX*I_ERI_D2x_S_S_S_M2_vrr+WPX*I_ERI_D2x_S_S_S_M3_vrr+2*oned2z*I_ERI_Px_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_F2xy_S_S_S_M2_vrr = PAY*I_ERI_D2x_S_S_S_M2_vrr+WPY*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_F2xz_S_S_S_M2_vrr = PAZ*I_ERI_D2x_S_S_S_M2_vrr+WPZ*I_ERI_D2x_S_S_S_M3_vrr;
      Double I_ERI_Fx2y_S_S_S_M2_vrr = PAX*I_ERI_D2y_S_S_S_M2_vrr+WPX*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_Fx2z_S_S_S_M2_vrr = PAX*I_ERI_D2z_S_S_S_M2_vrr+WPX*I_ERI_D2z_S_S_S_M3_vrr;
      Double I_ERI_F3y_S_S_S_M2_vrr = PAY*I_ERI_D2y_S_S_S_M2_vrr+WPY*I_ERI_D2y_S_S_S_M3_vrr+2*oned2z*I_ERI_Py_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_F2yz_S_S_S_M2_vrr = PAZ*I_ERI_D2y_S_S_S_M2_vrr+WPZ*I_ERI_D2y_S_S_S_M3_vrr;
      Double I_ERI_F3z_S_S_S_M2_vrr = PAZ*I_ERI_D2z_S_S_S_M2_vrr+WPZ*I_ERI_D2z_S_S_S_M3_vrr+2*oned2z*I_ERI_Pz_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_G_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_Px_S_vrr = QCX*I_ERI_G4x_S_S_S_vrr+WQX*I_ERI_G4x_S_S_S_M1_vrr+4*oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_Px_S_vrr = QCX*I_ERI_G3xy_S_S_S_vrr+WQX*I_ERI_G3xy_S_S_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_Px_S_vrr = QCX*I_ERI_G3xz_S_S_S_vrr+WQX*I_ERI_G3xz_S_S_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_Px_S_vrr = QCX*I_ERI_G2x2y_S_S_S_vrr+WQX*I_ERI_G2x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_Px_S_vrr = QCX*I_ERI_G2xyz_S_S_S_vrr+WQX*I_ERI_G2xyz_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_Px_S_vrr = QCX*I_ERI_G2x2z_S_S_S_vrr+WQX*I_ERI_G2x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_Px_S_vrr = QCX*I_ERI_Gx3y_S_S_S_vrr+WQX*I_ERI_Gx3y_S_S_S_M1_vrr+oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Px_S_vrr = QCX*I_ERI_Gx2yz_S_S_S_vrr+WQX*I_ERI_Gx2yz_S_S_S_M1_vrr+oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Px_S_vrr = QCX*I_ERI_Gxy2z_S_S_S_vrr+WQX*I_ERI_Gxy2z_S_S_S_M1_vrr+oned2k*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_Px_S_vrr = QCX*I_ERI_Gx3z_S_S_S_vrr+WQX*I_ERI_Gx3z_S_S_S_M1_vrr+oned2k*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_Px_S_vrr = QCX*I_ERI_G4y_S_S_S_vrr+WQX*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_Px_S_vrr = QCX*I_ERI_G3yz_S_S_S_vrr+WQX*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_Px_S_vrr = QCX*I_ERI_G2y2z_S_S_S_vrr+WQX*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_Px_S_vrr = QCX*I_ERI_Gy3z_S_S_S_vrr+WQX*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_Px_S_vrr = QCX*I_ERI_G4z_S_S_S_vrr+WQX*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_G4x_S_Py_S_vrr = QCY*I_ERI_G4x_S_S_S_vrr+WQY*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_Py_S_vrr = QCY*I_ERI_G3xy_S_S_S_vrr+WQY*I_ERI_G3xy_S_S_S_M1_vrr+oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_Py_S_vrr = QCY*I_ERI_G3xz_S_S_S_vrr+WQY*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_Py_S_vrr = QCY*I_ERI_G2x2y_S_S_S_vrr+WQY*I_ERI_G2x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_Py_S_vrr = QCY*I_ERI_G2xyz_S_S_S_vrr+WQY*I_ERI_G2xyz_S_S_S_M1_vrr+oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_Py_S_vrr = QCY*I_ERI_G2x2z_S_S_S_vrr+WQY*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_Py_S_vrr = QCY*I_ERI_Gx3y_S_S_S_vrr+WQY*I_ERI_Gx3y_S_S_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Py_S_vrr = QCY*I_ERI_Gx2yz_S_S_S_vrr+WQY*I_ERI_Gx2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Py_S_vrr = QCY*I_ERI_Gxy2z_S_S_S_vrr+WQY*I_ERI_Gxy2z_S_S_S_M1_vrr+oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_Py_S_vrr = QCY*I_ERI_Gx3z_S_S_S_vrr+WQY*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_Py_S_vrr = QCY*I_ERI_G4y_S_S_S_vrr+WQY*I_ERI_G4y_S_S_S_M1_vrr+4*oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_Py_S_vrr = QCY*I_ERI_G3yz_S_S_S_vrr+WQY*I_ERI_G3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_Py_S_vrr = QCY*I_ERI_G2y2z_S_S_S_vrr+WQY*I_ERI_G2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_Py_S_vrr = QCY*I_ERI_Gy3z_S_S_S_vrr+WQY*I_ERI_Gy3z_S_S_S_M1_vrr+oned2k*I_ERI_F3z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_Py_S_vrr = QCY*I_ERI_G4z_S_S_S_vrr+WQY*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_G4x_S_Pz_S_vrr = QCZ*I_ERI_G4x_S_S_S_vrr+WQZ*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_G3xy_S_Pz_S_vrr = QCZ*I_ERI_G3xy_S_S_S_vrr+WQZ*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_G3xz_S_Pz_S_vrr = QCZ*I_ERI_G3xz_S_S_S_vrr+WQZ*I_ERI_G3xz_S_S_S_M1_vrr+oned2k*I_ERI_F3x_S_S_S_M1_vrr;
      Double I_ERI_G2x2y_S_Pz_S_vrr = QCZ*I_ERI_G2x2y_S_S_S_vrr+WQZ*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_G2xyz_S_Pz_S_vrr = QCZ*I_ERI_G2xyz_S_S_S_vrr+WQZ*I_ERI_G2xyz_S_S_S_M1_vrr+oned2k*I_ERI_F2xy_S_S_S_M1_vrr;
      Double I_ERI_G2x2z_S_Pz_S_vrr = QCZ*I_ERI_G2x2z_S_S_S_vrr+WQZ*I_ERI_G2x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_S_S_M1_vrr;
      Double I_ERI_Gx3y_S_Pz_S_vrr = QCZ*I_ERI_Gx3y_S_S_S_vrr+WQZ*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Pz_S_vrr = QCZ*I_ERI_Gx2yz_S_S_S_vrr+WQZ*I_ERI_Gx2yz_S_S_S_M1_vrr+oned2k*I_ERI_Fx2y_S_S_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Pz_S_vrr = QCZ*I_ERI_Gxy2z_S_S_S_vrr+WQZ*I_ERI_Gxy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M1_vrr;
      Double I_ERI_Gx3z_S_Pz_S_vrr = QCZ*I_ERI_Gx3z_S_S_S_vrr+WQZ*I_ERI_Gx3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_S_S_M1_vrr;
      Double I_ERI_G4y_S_Pz_S_vrr = QCZ*I_ERI_G4y_S_S_S_vrr+WQZ*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_G3yz_S_Pz_S_vrr = QCZ*I_ERI_G3yz_S_S_S_vrr+WQZ*I_ERI_G3yz_S_S_S_M1_vrr+oned2k*I_ERI_F3y_S_S_S_M1_vrr;
      Double I_ERI_G2y2z_S_Pz_S_vrr = QCZ*I_ERI_G2y2z_S_S_S_vrr+WQZ*I_ERI_G2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_S_S_M1_vrr;
      Double I_ERI_Gy3z_S_Pz_S_vrr = QCZ*I_ERI_Gy3z_S_S_S_vrr+WQZ*I_ERI_Gy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_S_S_M1_vrr;
      Double I_ERI_G4z_S_Pz_S_vrr = QCZ*I_ERI_G4z_S_S_S_vrr+WQZ*I_ERI_G4z_S_S_S_M1_vrr+4*oned2k*I_ERI_F3z_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_G_S_S_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_G4x_S_S_S += I_ERI_G4x_S_S_S_vrr;
      I_ERI_G3xy_S_S_S += I_ERI_G3xy_S_S_S_vrr;
      I_ERI_G3xz_S_S_S += I_ERI_G3xz_S_S_S_vrr;
      I_ERI_G2x2y_S_S_S += I_ERI_G2x2y_S_S_S_vrr;
      I_ERI_G2xyz_S_S_S += I_ERI_G2xyz_S_S_S_vrr;
      I_ERI_G2x2z_S_S_S += I_ERI_G2x2z_S_S_S_vrr;
      I_ERI_Gx3y_S_S_S += I_ERI_Gx3y_S_S_S_vrr;
      I_ERI_Gx2yz_S_S_S += I_ERI_Gx2yz_S_S_S_vrr;
      I_ERI_Gxy2z_S_S_S += I_ERI_Gxy2z_S_S_S_vrr;
      I_ERI_Gx3z_S_S_S += I_ERI_Gx3z_S_S_S_vrr;
      I_ERI_G4y_S_S_S += I_ERI_G4y_S_S_S_vrr;
      I_ERI_G3yz_S_S_S += I_ERI_G3yz_S_S_S_vrr;
      I_ERI_G2y2z_S_S_S += I_ERI_G2y2z_S_S_S_vrr;
      I_ERI_Gy3z_S_S_S += I_ERI_Gy3z_S_S_S_vrr;
      I_ERI_G4z_S_S_S += I_ERI_G4z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_I_S_S_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_I_S_S_S_a_coefs = alpha;
      I_ERI_I6x_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I6x_S_S_S_vrr;
      I_ERI_I5xy_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I5xy_S_S_S_vrr;
      I_ERI_I5xz_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I5xz_S_S_S_vrr;
      I_ERI_I4x2y_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I4x2y_S_S_S_vrr;
      I_ERI_I4xyz_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I4xyz_S_S_S_vrr;
      I_ERI_I4x2z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I4x2z_S_S_S_vrr;
      I_ERI_I3x3y_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I3x3y_S_S_S_vrr;
      I_ERI_I3x2yz_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I3x2yz_S_S_S_vrr;
      I_ERI_I3xy2z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I3xy2z_S_S_S_vrr;
      I_ERI_I3x3z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I3x3z_S_S_S_vrr;
      I_ERI_I2x4y_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I2x4y_S_S_S_vrr;
      I_ERI_I2x3yz_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I2x3yz_S_S_S_vrr;
      I_ERI_I2x2y2z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I2x2y2z_S_S_S_vrr;
      I_ERI_I2xy3z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I2xy3z_S_S_S_vrr;
      I_ERI_I2x4z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I2x4z_S_S_S_vrr;
      I_ERI_Ix5y_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_Ix5y_S_S_S_vrr;
      I_ERI_Ix4yz_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_Ix4yz_S_S_S_vrr;
      I_ERI_Ix3y2z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_Ix3y2z_S_S_S_vrr;
      I_ERI_Ix2y3z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_Ix2y3z_S_S_S_vrr;
      I_ERI_Ixy4z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_Ixy4z_S_S_S_vrr;
      I_ERI_Ix5z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_Ix5z_S_S_S_vrr;
      I_ERI_I6y_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I6y_S_S_S_vrr;
      I_ERI_I5yz_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I5yz_S_S_S_vrr;
      I_ERI_I4y2z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I4y2z_S_S_S_vrr;
      I_ERI_I3y3z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I3y3z_S_S_S_vrr;
      I_ERI_I2y4z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I2y4z_S_S_S_vrr;
      I_ERI_Iy5z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_Iy5z_S_S_S_vrr;
      I_ERI_I6z_S_S_S_a += SQ_ERI_I_S_S_S_a_coefs*I_ERI_I6z_S_S_S_vrr;

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
       * shell quartet name: SQ_ERI_H_S_P_S_c
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_H_S_P_S_c_coefs = gamma;
      I_ERI_H5x_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H5x_S_Px_S_vrr;
      I_ERI_H4xy_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H4xy_S_Px_S_vrr;
      I_ERI_H4xz_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H4xz_S_Px_S_vrr;
      I_ERI_H3x2y_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3x2y_S_Px_S_vrr;
      I_ERI_H3xyz_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3xyz_S_Px_S_vrr;
      I_ERI_H3x2z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3x2z_S_Px_S_vrr;
      I_ERI_H2x3y_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2x3y_S_Px_S_vrr;
      I_ERI_H2x2yz_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2x2yz_S_Px_S_vrr;
      I_ERI_H2xy2z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2xy2z_S_Px_S_vrr;
      I_ERI_H2x3z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2x3z_S_Px_S_vrr;
      I_ERI_Hx4y_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx4y_S_Px_S_vrr;
      I_ERI_Hx3yz_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx3yz_S_Px_S_vrr;
      I_ERI_Hx2y2z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx2y2z_S_Px_S_vrr;
      I_ERI_Hxy3z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hxy3z_S_Px_S_vrr;
      I_ERI_Hx4z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx4z_S_Px_S_vrr;
      I_ERI_H5y_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H5y_S_Px_S_vrr;
      I_ERI_H4yz_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H4yz_S_Px_S_vrr;
      I_ERI_H3y2z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3y2z_S_Px_S_vrr;
      I_ERI_H2y3z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2y3z_S_Px_S_vrr;
      I_ERI_Hy4z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hy4z_S_Px_S_vrr;
      I_ERI_H5z_S_Px_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H5z_S_Px_S_vrr;
      I_ERI_H5x_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H5x_S_Py_S_vrr;
      I_ERI_H4xy_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H4xy_S_Py_S_vrr;
      I_ERI_H4xz_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H4xz_S_Py_S_vrr;
      I_ERI_H3x2y_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3x2y_S_Py_S_vrr;
      I_ERI_H3xyz_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3xyz_S_Py_S_vrr;
      I_ERI_H3x2z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3x2z_S_Py_S_vrr;
      I_ERI_H2x3y_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2x3y_S_Py_S_vrr;
      I_ERI_H2x2yz_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2x2yz_S_Py_S_vrr;
      I_ERI_H2xy2z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2xy2z_S_Py_S_vrr;
      I_ERI_H2x3z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2x3z_S_Py_S_vrr;
      I_ERI_Hx4y_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx4y_S_Py_S_vrr;
      I_ERI_Hx3yz_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx3yz_S_Py_S_vrr;
      I_ERI_Hx2y2z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx2y2z_S_Py_S_vrr;
      I_ERI_Hxy3z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hxy3z_S_Py_S_vrr;
      I_ERI_Hx4z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx4z_S_Py_S_vrr;
      I_ERI_H5y_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H5y_S_Py_S_vrr;
      I_ERI_H4yz_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H4yz_S_Py_S_vrr;
      I_ERI_H3y2z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3y2z_S_Py_S_vrr;
      I_ERI_H2y3z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2y3z_S_Py_S_vrr;
      I_ERI_Hy4z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hy4z_S_Py_S_vrr;
      I_ERI_H5z_S_Py_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H5z_S_Py_S_vrr;
      I_ERI_H5x_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H5x_S_Pz_S_vrr;
      I_ERI_H4xy_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H4xy_S_Pz_S_vrr;
      I_ERI_H4xz_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H4xz_S_Pz_S_vrr;
      I_ERI_H3x2y_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3x2y_S_Pz_S_vrr;
      I_ERI_H3xyz_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3xyz_S_Pz_S_vrr;
      I_ERI_H3x2z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3x2z_S_Pz_S_vrr;
      I_ERI_H2x3y_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2x3y_S_Pz_S_vrr;
      I_ERI_H2x2yz_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2x2yz_S_Pz_S_vrr;
      I_ERI_H2xy2z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2xy2z_S_Pz_S_vrr;
      I_ERI_H2x3z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2x3z_S_Pz_S_vrr;
      I_ERI_Hx4y_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx4y_S_Pz_S_vrr;
      I_ERI_Hx3yz_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx3yz_S_Pz_S_vrr;
      I_ERI_Hx2y2z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx2y2z_S_Pz_S_vrr;
      I_ERI_Hxy3z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hxy3z_S_Pz_S_vrr;
      I_ERI_Hx4z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hx4z_S_Pz_S_vrr;
      I_ERI_H5y_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H5y_S_Pz_S_vrr;
      I_ERI_H4yz_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H4yz_S_Pz_S_vrr;
      I_ERI_H3y2z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H3y2z_S_Pz_S_vrr;
      I_ERI_H2y3z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H2y3z_S_Pz_S_vrr;
      I_ERI_Hy4z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_Hy4z_S_Pz_S_vrr;
      I_ERI_H5z_S_Pz_S_c += SQ_ERI_H_S_P_S_c_coefs*I_ERI_H5z_S_Pz_S_vrr;

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
       * shell quartet name: SQ_ERI_I_S_S_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_I_S_S_S_b_coefs = beta;
      I_ERI_I6x_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I6x_S_S_S_vrr;
      I_ERI_I5xy_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I5xy_S_S_S_vrr;
      I_ERI_I5xz_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I5xz_S_S_S_vrr;
      I_ERI_I4x2y_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I4x2y_S_S_S_vrr;
      I_ERI_I4xyz_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I4xyz_S_S_S_vrr;
      I_ERI_I4x2z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I4x2z_S_S_S_vrr;
      I_ERI_I3x3y_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I3x3y_S_S_S_vrr;
      I_ERI_I3x2yz_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I3x2yz_S_S_S_vrr;
      I_ERI_I3xy2z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I3xy2z_S_S_S_vrr;
      I_ERI_I3x3z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I3x3z_S_S_S_vrr;
      I_ERI_I2x4y_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I2x4y_S_S_S_vrr;
      I_ERI_I2x3yz_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I2x3yz_S_S_S_vrr;
      I_ERI_I2x2y2z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I2x2y2z_S_S_S_vrr;
      I_ERI_I2xy3z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I2xy3z_S_S_S_vrr;
      I_ERI_I2x4z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I2x4z_S_S_S_vrr;
      I_ERI_Ix5y_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_Ix5y_S_S_S_vrr;
      I_ERI_Ix4yz_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_Ix4yz_S_S_S_vrr;
      I_ERI_Ix3y2z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_Ix3y2z_S_S_S_vrr;
      I_ERI_Ix2y3z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_Ix2y3z_S_S_S_vrr;
      I_ERI_Ixy4z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_Ixy4z_S_S_S_vrr;
      I_ERI_Ix5z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_Ix5z_S_S_S_vrr;
      I_ERI_I6y_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I6y_S_S_S_vrr;
      I_ERI_I5yz_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I5yz_S_S_S_vrr;
      I_ERI_I4y2z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I4y2z_S_S_S_vrr;
      I_ERI_I3y3z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I3y3z_S_S_S_vrr;
      I_ERI_I2y4z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I2y4z_S_S_S_vrr;
      I_ERI_Iy5z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_Iy5z_S_S_S_vrr;
      I_ERI_I6z_S_S_S_b += SQ_ERI_I_S_S_S_b_coefs*I_ERI_I6z_S_S_S_vrr;

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
   * shell quartet name: SQ_ERI_F_P_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_S_S
   * RHS shell quartet name: SQ_ERI_F_S_S_S
   ************************************************************/
  Double I_ERI_F3x_Px_S_S = I_ERI_G4x_S_S_S+ABX*I_ERI_F3x_S_S_S;
  Double I_ERI_F2xy_Px_S_S = I_ERI_G3xy_S_S_S+ABX*I_ERI_F2xy_S_S_S;
  Double I_ERI_F2xz_Px_S_S = I_ERI_G3xz_S_S_S+ABX*I_ERI_F2xz_S_S_S;
  Double I_ERI_Fx2y_Px_S_S = I_ERI_G2x2y_S_S_S+ABX*I_ERI_Fx2y_S_S_S;
  Double I_ERI_Fxyz_Px_S_S = I_ERI_G2xyz_S_S_S+ABX*I_ERI_Fxyz_S_S_S;
  Double I_ERI_Fx2z_Px_S_S = I_ERI_G2x2z_S_S_S+ABX*I_ERI_Fx2z_S_S_S;
  Double I_ERI_F3y_Px_S_S = I_ERI_Gx3y_S_S_S+ABX*I_ERI_F3y_S_S_S;
  Double I_ERI_F2yz_Px_S_S = I_ERI_Gx2yz_S_S_S+ABX*I_ERI_F2yz_S_S_S;
  Double I_ERI_Fy2z_Px_S_S = I_ERI_Gxy2z_S_S_S+ABX*I_ERI_Fy2z_S_S_S;
  Double I_ERI_F3z_Px_S_S = I_ERI_Gx3z_S_S_S+ABX*I_ERI_F3z_S_S_S;
  Double I_ERI_F3x_Py_S_S = I_ERI_G3xy_S_S_S+ABY*I_ERI_F3x_S_S_S;
  Double I_ERI_F2xy_Py_S_S = I_ERI_G2x2y_S_S_S+ABY*I_ERI_F2xy_S_S_S;
  Double I_ERI_F2xz_Py_S_S = I_ERI_G2xyz_S_S_S+ABY*I_ERI_F2xz_S_S_S;
  Double I_ERI_Fx2y_Py_S_S = I_ERI_Gx3y_S_S_S+ABY*I_ERI_Fx2y_S_S_S;
  Double I_ERI_Fxyz_Py_S_S = I_ERI_Gx2yz_S_S_S+ABY*I_ERI_Fxyz_S_S_S;
  Double I_ERI_Fx2z_Py_S_S = I_ERI_Gxy2z_S_S_S+ABY*I_ERI_Fx2z_S_S_S;
  Double I_ERI_F3y_Py_S_S = I_ERI_G4y_S_S_S+ABY*I_ERI_F3y_S_S_S;
  Double I_ERI_F2yz_Py_S_S = I_ERI_G3yz_S_S_S+ABY*I_ERI_F2yz_S_S_S;
  Double I_ERI_Fy2z_Py_S_S = I_ERI_G2y2z_S_S_S+ABY*I_ERI_Fy2z_S_S_S;
  Double I_ERI_F3z_Py_S_S = I_ERI_Gy3z_S_S_S+ABY*I_ERI_F3z_S_S_S;
  Double I_ERI_F3x_Pz_S_S = I_ERI_G3xz_S_S_S+ABZ*I_ERI_F3x_S_S_S;
  Double I_ERI_F2xy_Pz_S_S = I_ERI_G2xyz_S_S_S+ABZ*I_ERI_F2xy_S_S_S;
  Double I_ERI_F2xz_Pz_S_S = I_ERI_G2x2z_S_S_S+ABZ*I_ERI_F2xz_S_S_S;
  Double I_ERI_Fx2y_Pz_S_S = I_ERI_Gx2yz_S_S_S+ABZ*I_ERI_Fx2y_S_S_S;
  Double I_ERI_Fxyz_Pz_S_S = I_ERI_Gxy2z_S_S_S+ABZ*I_ERI_Fxyz_S_S_S;
  Double I_ERI_Fx2z_Pz_S_S = I_ERI_Gx3z_S_S_S+ABZ*I_ERI_Fx2z_S_S_S;
  Double I_ERI_F3y_Pz_S_S = I_ERI_G3yz_S_S_S+ABZ*I_ERI_F3y_S_S_S;
  Double I_ERI_F2yz_Pz_S_S = I_ERI_G2y2z_S_S_S+ABZ*I_ERI_F2yz_S_S_S;
  Double I_ERI_Fy2z_Pz_S_S = I_ERI_Gy3z_S_S_S+ABZ*I_ERI_Fy2z_S_S_S;
  Double I_ERI_F3z_Pz_S_S = I_ERI_G4z_S_S_S+ABZ*I_ERI_F3z_S_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_S_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_S_S
   * RHS shell quartet name: SQ_ERI_D_P_S_S
   ************************************************************/
  Double I_ERI_D2x_D2x_S_S = I_ERI_F3x_Px_S_S+ABX*I_ERI_D2x_Px_S_S;
  Double I_ERI_Dxy_D2x_S_S = I_ERI_F2xy_Px_S_S+ABX*I_ERI_Dxy_Px_S_S;
  Double I_ERI_Dxz_D2x_S_S = I_ERI_F2xz_Px_S_S+ABX*I_ERI_Dxz_Px_S_S;
  Double I_ERI_D2y_D2x_S_S = I_ERI_Fx2y_Px_S_S+ABX*I_ERI_D2y_Px_S_S;
  Double I_ERI_Dyz_D2x_S_S = I_ERI_Fxyz_Px_S_S+ABX*I_ERI_Dyz_Px_S_S;
  Double I_ERI_D2z_D2x_S_S = I_ERI_Fx2z_Px_S_S+ABX*I_ERI_D2z_Px_S_S;
  Double I_ERI_D2x_Dxy_S_S = I_ERI_F2xy_Px_S_S+ABY*I_ERI_D2x_Px_S_S;
  Double I_ERI_Dxy_Dxy_S_S = I_ERI_Fx2y_Px_S_S+ABY*I_ERI_Dxy_Px_S_S;
  Double I_ERI_Dxz_Dxy_S_S = I_ERI_Fxyz_Px_S_S+ABY*I_ERI_Dxz_Px_S_S;
  Double I_ERI_D2y_Dxy_S_S = I_ERI_F3y_Px_S_S+ABY*I_ERI_D2y_Px_S_S;
  Double I_ERI_Dyz_Dxy_S_S = I_ERI_F2yz_Px_S_S+ABY*I_ERI_Dyz_Px_S_S;
  Double I_ERI_D2z_Dxy_S_S = I_ERI_Fy2z_Px_S_S+ABY*I_ERI_D2z_Px_S_S;
  Double I_ERI_D2x_Dxz_S_S = I_ERI_F2xz_Px_S_S+ABZ*I_ERI_D2x_Px_S_S;
  Double I_ERI_Dxy_Dxz_S_S = I_ERI_Fxyz_Px_S_S+ABZ*I_ERI_Dxy_Px_S_S;
  Double I_ERI_Dxz_Dxz_S_S = I_ERI_Fx2z_Px_S_S+ABZ*I_ERI_Dxz_Px_S_S;
  Double I_ERI_D2y_Dxz_S_S = I_ERI_F2yz_Px_S_S+ABZ*I_ERI_D2y_Px_S_S;
  Double I_ERI_Dyz_Dxz_S_S = I_ERI_Fy2z_Px_S_S+ABZ*I_ERI_Dyz_Px_S_S;
  Double I_ERI_D2z_Dxz_S_S = I_ERI_F3z_Px_S_S+ABZ*I_ERI_D2z_Px_S_S;
  Double I_ERI_D2x_D2y_S_S = I_ERI_F2xy_Py_S_S+ABY*I_ERI_D2x_Py_S_S;
  Double I_ERI_Dxy_D2y_S_S = I_ERI_Fx2y_Py_S_S+ABY*I_ERI_Dxy_Py_S_S;
  Double I_ERI_Dxz_D2y_S_S = I_ERI_Fxyz_Py_S_S+ABY*I_ERI_Dxz_Py_S_S;
  Double I_ERI_D2y_D2y_S_S = I_ERI_F3y_Py_S_S+ABY*I_ERI_D2y_Py_S_S;
  Double I_ERI_Dyz_D2y_S_S = I_ERI_F2yz_Py_S_S+ABY*I_ERI_Dyz_Py_S_S;
  Double I_ERI_D2z_D2y_S_S = I_ERI_Fy2z_Py_S_S+ABY*I_ERI_D2z_Py_S_S;
  Double I_ERI_D2x_Dyz_S_S = I_ERI_F2xz_Py_S_S+ABZ*I_ERI_D2x_Py_S_S;
  Double I_ERI_Dxy_Dyz_S_S = I_ERI_Fxyz_Py_S_S+ABZ*I_ERI_Dxy_Py_S_S;
  Double I_ERI_Dxz_Dyz_S_S = I_ERI_Fx2z_Py_S_S+ABZ*I_ERI_Dxz_Py_S_S;
  Double I_ERI_D2y_Dyz_S_S = I_ERI_F2yz_Py_S_S+ABZ*I_ERI_D2y_Py_S_S;
  Double I_ERI_Dyz_Dyz_S_S = I_ERI_Fy2z_Py_S_S+ABZ*I_ERI_Dyz_Py_S_S;
  Double I_ERI_D2z_Dyz_S_S = I_ERI_F3z_Py_S_S+ABZ*I_ERI_D2z_Py_S_S;
  Double I_ERI_D2x_D2z_S_S = I_ERI_F2xz_Pz_S_S+ABZ*I_ERI_D2x_Pz_S_S;
  Double I_ERI_Dxy_D2z_S_S = I_ERI_Fxyz_Pz_S_S+ABZ*I_ERI_Dxy_Pz_S_S;
  Double I_ERI_Dxz_D2z_S_S = I_ERI_Fx2z_Pz_S_S+ABZ*I_ERI_Dxz_Pz_S_S;
  Double I_ERI_D2y_D2z_S_S = I_ERI_F2yz_Pz_S_S+ABZ*I_ERI_D2y_Pz_S_S;
  Double I_ERI_Dyz_D2z_S_S = I_ERI_Fy2z_Pz_S_S+ABZ*I_ERI_Dyz_Pz_S_S;
  Double I_ERI_D2z_D2z_S_S = I_ERI_F3z_Pz_S_S+ABZ*I_ERI_D2z_Pz_S_S;

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
   * shell quartet name: SQ_ERI_H_P_S_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 7 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_S_S_a
   * RHS shell quartet name: SQ_ERI_H_S_S_S_a
   ************************************************************/
  Double I_ERI_H5x_Px_S_S_a = I_ERI_I6x_S_S_S_a+ABX*I_ERI_H5x_S_S_S_a;
  Double I_ERI_H4xy_Px_S_S_a = I_ERI_I5xy_S_S_S_a+ABX*I_ERI_H4xy_S_S_S_a;
  Double I_ERI_H4xz_Px_S_S_a = I_ERI_I5xz_S_S_S_a+ABX*I_ERI_H4xz_S_S_S_a;
  Double I_ERI_H3x2y_Px_S_S_a = I_ERI_I4x2y_S_S_S_a+ABX*I_ERI_H3x2y_S_S_S_a;
  Double I_ERI_H3xyz_Px_S_S_a = I_ERI_I4xyz_S_S_S_a+ABX*I_ERI_H3xyz_S_S_S_a;
  Double I_ERI_H3x2z_Px_S_S_a = I_ERI_I4x2z_S_S_S_a+ABX*I_ERI_H3x2z_S_S_S_a;
  Double I_ERI_H2x3y_Px_S_S_a = I_ERI_I3x3y_S_S_S_a+ABX*I_ERI_H2x3y_S_S_S_a;
  Double I_ERI_H2x2yz_Px_S_S_a = I_ERI_I3x2yz_S_S_S_a+ABX*I_ERI_H2x2yz_S_S_S_a;
  Double I_ERI_H2xy2z_Px_S_S_a = I_ERI_I3xy2z_S_S_S_a+ABX*I_ERI_H2xy2z_S_S_S_a;
  Double I_ERI_H2x3z_Px_S_S_a = I_ERI_I3x3z_S_S_S_a+ABX*I_ERI_H2x3z_S_S_S_a;
  Double I_ERI_Hx4y_Px_S_S_a = I_ERI_I2x4y_S_S_S_a+ABX*I_ERI_Hx4y_S_S_S_a;
  Double I_ERI_Hx3yz_Px_S_S_a = I_ERI_I2x3yz_S_S_S_a+ABX*I_ERI_Hx3yz_S_S_S_a;
  Double I_ERI_Hx2y2z_Px_S_S_a = I_ERI_I2x2y2z_S_S_S_a+ABX*I_ERI_Hx2y2z_S_S_S_a;
  Double I_ERI_Hxy3z_Px_S_S_a = I_ERI_I2xy3z_S_S_S_a+ABX*I_ERI_Hxy3z_S_S_S_a;
  Double I_ERI_Hx4z_Px_S_S_a = I_ERI_I2x4z_S_S_S_a+ABX*I_ERI_Hx4z_S_S_S_a;
  Double I_ERI_H5y_Px_S_S_a = I_ERI_Ix5y_S_S_S_a+ABX*I_ERI_H5y_S_S_S_a;
  Double I_ERI_H4yz_Px_S_S_a = I_ERI_Ix4yz_S_S_S_a+ABX*I_ERI_H4yz_S_S_S_a;
  Double I_ERI_H3y2z_Px_S_S_a = I_ERI_Ix3y2z_S_S_S_a+ABX*I_ERI_H3y2z_S_S_S_a;
  Double I_ERI_H2y3z_Px_S_S_a = I_ERI_Ix2y3z_S_S_S_a+ABX*I_ERI_H2y3z_S_S_S_a;
  Double I_ERI_Hy4z_Px_S_S_a = I_ERI_Ixy4z_S_S_S_a+ABX*I_ERI_Hy4z_S_S_S_a;
  Double I_ERI_H5z_Px_S_S_a = I_ERI_Ix5z_S_S_S_a+ABX*I_ERI_H5z_S_S_S_a;
  Double I_ERI_H4xy_Py_S_S_a = I_ERI_I4x2y_S_S_S_a+ABY*I_ERI_H4xy_S_S_S_a;
  Double I_ERI_H4xz_Py_S_S_a = I_ERI_I4xyz_S_S_S_a+ABY*I_ERI_H4xz_S_S_S_a;
  Double I_ERI_H3x2y_Py_S_S_a = I_ERI_I3x3y_S_S_S_a+ABY*I_ERI_H3x2y_S_S_S_a;
  Double I_ERI_H3xyz_Py_S_S_a = I_ERI_I3x2yz_S_S_S_a+ABY*I_ERI_H3xyz_S_S_S_a;
  Double I_ERI_H3x2z_Py_S_S_a = I_ERI_I3xy2z_S_S_S_a+ABY*I_ERI_H3x2z_S_S_S_a;
  Double I_ERI_H2x3y_Py_S_S_a = I_ERI_I2x4y_S_S_S_a+ABY*I_ERI_H2x3y_S_S_S_a;
  Double I_ERI_H2x2yz_Py_S_S_a = I_ERI_I2x3yz_S_S_S_a+ABY*I_ERI_H2x2yz_S_S_S_a;
  Double I_ERI_H2xy2z_Py_S_S_a = I_ERI_I2x2y2z_S_S_S_a+ABY*I_ERI_H2xy2z_S_S_S_a;
  Double I_ERI_H2x3z_Py_S_S_a = I_ERI_I2xy3z_S_S_S_a+ABY*I_ERI_H2x3z_S_S_S_a;
  Double I_ERI_Hx4y_Py_S_S_a = I_ERI_Ix5y_S_S_S_a+ABY*I_ERI_Hx4y_S_S_S_a;
  Double I_ERI_Hx3yz_Py_S_S_a = I_ERI_Ix4yz_S_S_S_a+ABY*I_ERI_Hx3yz_S_S_S_a;
  Double I_ERI_Hx2y2z_Py_S_S_a = I_ERI_Ix3y2z_S_S_S_a+ABY*I_ERI_Hx2y2z_S_S_S_a;
  Double I_ERI_Hxy3z_Py_S_S_a = I_ERI_Ix2y3z_S_S_S_a+ABY*I_ERI_Hxy3z_S_S_S_a;
  Double I_ERI_Hx4z_Py_S_S_a = I_ERI_Ixy4z_S_S_S_a+ABY*I_ERI_Hx4z_S_S_S_a;
  Double I_ERI_H5y_Py_S_S_a = I_ERI_I6y_S_S_S_a+ABY*I_ERI_H5y_S_S_S_a;
  Double I_ERI_H4yz_Py_S_S_a = I_ERI_I5yz_S_S_S_a+ABY*I_ERI_H4yz_S_S_S_a;
  Double I_ERI_H3y2z_Py_S_S_a = I_ERI_I4y2z_S_S_S_a+ABY*I_ERI_H3y2z_S_S_S_a;
  Double I_ERI_H2y3z_Py_S_S_a = I_ERI_I3y3z_S_S_S_a+ABY*I_ERI_H2y3z_S_S_S_a;
  Double I_ERI_Hy4z_Py_S_S_a = I_ERI_I2y4z_S_S_S_a+ABY*I_ERI_Hy4z_S_S_S_a;
  Double I_ERI_H5z_Py_S_S_a = I_ERI_Iy5z_S_S_S_a+ABY*I_ERI_H5z_S_S_S_a;
  Double I_ERI_H4xz_Pz_S_S_a = I_ERI_I4x2z_S_S_S_a+ABZ*I_ERI_H4xz_S_S_S_a;
  Double I_ERI_H3xyz_Pz_S_S_a = I_ERI_I3xy2z_S_S_S_a+ABZ*I_ERI_H3xyz_S_S_S_a;
  Double I_ERI_H3x2z_Pz_S_S_a = I_ERI_I3x3z_S_S_S_a+ABZ*I_ERI_H3x2z_S_S_S_a;
  Double I_ERI_H2x2yz_Pz_S_S_a = I_ERI_I2x2y2z_S_S_S_a+ABZ*I_ERI_H2x2yz_S_S_S_a;
  Double I_ERI_H2xy2z_Pz_S_S_a = I_ERI_I2xy3z_S_S_S_a+ABZ*I_ERI_H2xy2z_S_S_S_a;
  Double I_ERI_H2x3z_Pz_S_S_a = I_ERI_I2x4z_S_S_S_a+ABZ*I_ERI_H2x3z_S_S_S_a;
  Double I_ERI_Hx3yz_Pz_S_S_a = I_ERI_Ix3y2z_S_S_S_a+ABZ*I_ERI_Hx3yz_S_S_S_a;
  Double I_ERI_Hx2y2z_Pz_S_S_a = I_ERI_Ix2y3z_S_S_S_a+ABZ*I_ERI_Hx2y2z_S_S_S_a;
  Double I_ERI_Hxy3z_Pz_S_S_a = I_ERI_Ixy4z_S_S_S_a+ABZ*I_ERI_Hxy3z_S_S_S_a;
  Double I_ERI_Hx4z_Pz_S_S_a = I_ERI_Ix5z_S_S_S_a+ABZ*I_ERI_Hx4z_S_S_S_a;
  Double I_ERI_H4yz_Pz_S_S_a = I_ERI_I4y2z_S_S_S_a+ABZ*I_ERI_H4yz_S_S_S_a;
  Double I_ERI_H3y2z_Pz_S_S_a = I_ERI_I3y3z_S_S_S_a+ABZ*I_ERI_H3y2z_S_S_S_a;
  Double I_ERI_H2y3z_Pz_S_S_a = I_ERI_I2y4z_S_S_S_a+ABZ*I_ERI_H2y3z_S_S_S_a;
  Double I_ERI_Hy4z_Pz_S_S_a = I_ERI_Iy5z_S_S_S_a+ABZ*I_ERI_Hy4z_S_S_S_a;
  Double I_ERI_H5z_Pz_S_S_a = I_ERI_I6z_S_S_S_a+ABZ*I_ERI_H5z_S_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_G_D_S_S_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_P_S_S_a
   * RHS shell quartet name: SQ_ERI_G_P_S_S_a
   ************************************************************/
  Double I_ERI_G4x_D2x_S_S_a = I_ERI_H5x_Px_S_S_a+ABX*I_ERI_G4x_Px_S_S_a;
  Double I_ERI_G3xy_D2x_S_S_a = I_ERI_H4xy_Px_S_S_a+ABX*I_ERI_G3xy_Px_S_S_a;
  Double I_ERI_G3xz_D2x_S_S_a = I_ERI_H4xz_Px_S_S_a+ABX*I_ERI_G3xz_Px_S_S_a;
  Double I_ERI_G2x2y_D2x_S_S_a = I_ERI_H3x2y_Px_S_S_a+ABX*I_ERI_G2x2y_Px_S_S_a;
  Double I_ERI_G2xyz_D2x_S_S_a = I_ERI_H3xyz_Px_S_S_a+ABX*I_ERI_G2xyz_Px_S_S_a;
  Double I_ERI_G2x2z_D2x_S_S_a = I_ERI_H3x2z_Px_S_S_a+ABX*I_ERI_G2x2z_Px_S_S_a;
  Double I_ERI_Gx3y_D2x_S_S_a = I_ERI_H2x3y_Px_S_S_a+ABX*I_ERI_Gx3y_Px_S_S_a;
  Double I_ERI_Gx2yz_D2x_S_S_a = I_ERI_H2x2yz_Px_S_S_a+ABX*I_ERI_Gx2yz_Px_S_S_a;
  Double I_ERI_Gxy2z_D2x_S_S_a = I_ERI_H2xy2z_Px_S_S_a+ABX*I_ERI_Gxy2z_Px_S_S_a;
  Double I_ERI_Gx3z_D2x_S_S_a = I_ERI_H2x3z_Px_S_S_a+ABX*I_ERI_Gx3z_Px_S_S_a;
  Double I_ERI_G4y_D2x_S_S_a = I_ERI_Hx4y_Px_S_S_a+ABX*I_ERI_G4y_Px_S_S_a;
  Double I_ERI_G3yz_D2x_S_S_a = I_ERI_Hx3yz_Px_S_S_a+ABX*I_ERI_G3yz_Px_S_S_a;
  Double I_ERI_G2y2z_D2x_S_S_a = I_ERI_Hx2y2z_Px_S_S_a+ABX*I_ERI_G2y2z_Px_S_S_a;
  Double I_ERI_Gy3z_D2x_S_S_a = I_ERI_Hxy3z_Px_S_S_a+ABX*I_ERI_Gy3z_Px_S_S_a;
  Double I_ERI_G4z_D2x_S_S_a = I_ERI_Hx4z_Px_S_S_a+ABX*I_ERI_G4z_Px_S_S_a;
  Double I_ERI_G4x_Dxy_S_S_a = I_ERI_H4xy_Px_S_S_a+ABY*I_ERI_G4x_Px_S_S_a;
  Double I_ERI_G3xy_Dxy_S_S_a = I_ERI_H3x2y_Px_S_S_a+ABY*I_ERI_G3xy_Px_S_S_a;
  Double I_ERI_G3xz_Dxy_S_S_a = I_ERI_H3xyz_Px_S_S_a+ABY*I_ERI_G3xz_Px_S_S_a;
  Double I_ERI_G2x2y_Dxy_S_S_a = I_ERI_H2x3y_Px_S_S_a+ABY*I_ERI_G2x2y_Px_S_S_a;
  Double I_ERI_G2xyz_Dxy_S_S_a = I_ERI_H2x2yz_Px_S_S_a+ABY*I_ERI_G2xyz_Px_S_S_a;
  Double I_ERI_G2x2z_Dxy_S_S_a = I_ERI_H2xy2z_Px_S_S_a+ABY*I_ERI_G2x2z_Px_S_S_a;
  Double I_ERI_Gx3y_Dxy_S_S_a = I_ERI_Hx4y_Px_S_S_a+ABY*I_ERI_Gx3y_Px_S_S_a;
  Double I_ERI_Gx2yz_Dxy_S_S_a = I_ERI_Hx3yz_Px_S_S_a+ABY*I_ERI_Gx2yz_Px_S_S_a;
  Double I_ERI_Gxy2z_Dxy_S_S_a = I_ERI_Hx2y2z_Px_S_S_a+ABY*I_ERI_Gxy2z_Px_S_S_a;
  Double I_ERI_Gx3z_Dxy_S_S_a = I_ERI_Hxy3z_Px_S_S_a+ABY*I_ERI_Gx3z_Px_S_S_a;
  Double I_ERI_G4y_Dxy_S_S_a = I_ERI_H5y_Px_S_S_a+ABY*I_ERI_G4y_Px_S_S_a;
  Double I_ERI_G3yz_Dxy_S_S_a = I_ERI_H4yz_Px_S_S_a+ABY*I_ERI_G3yz_Px_S_S_a;
  Double I_ERI_G2y2z_Dxy_S_S_a = I_ERI_H3y2z_Px_S_S_a+ABY*I_ERI_G2y2z_Px_S_S_a;
  Double I_ERI_Gy3z_Dxy_S_S_a = I_ERI_H2y3z_Px_S_S_a+ABY*I_ERI_Gy3z_Px_S_S_a;
  Double I_ERI_G4z_Dxy_S_S_a = I_ERI_Hy4z_Px_S_S_a+ABY*I_ERI_G4z_Px_S_S_a;
  Double I_ERI_G4x_Dxz_S_S_a = I_ERI_H4xz_Px_S_S_a+ABZ*I_ERI_G4x_Px_S_S_a;
  Double I_ERI_G3xy_Dxz_S_S_a = I_ERI_H3xyz_Px_S_S_a+ABZ*I_ERI_G3xy_Px_S_S_a;
  Double I_ERI_G3xz_Dxz_S_S_a = I_ERI_H3x2z_Px_S_S_a+ABZ*I_ERI_G3xz_Px_S_S_a;
  Double I_ERI_G2x2y_Dxz_S_S_a = I_ERI_H2x2yz_Px_S_S_a+ABZ*I_ERI_G2x2y_Px_S_S_a;
  Double I_ERI_G2xyz_Dxz_S_S_a = I_ERI_H2xy2z_Px_S_S_a+ABZ*I_ERI_G2xyz_Px_S_S_a;
  Double I_ERI_G2x2z_Dxz_S_S_a = I_ERI_H2x3z_Px_S_S_a+ABZ*I_ERI_G2x2z_Px_S_S_a;
  Double I_ERI_Gx3y_Dxz_S_S_a = I_ERI_Hx3yz_Px_S_S_a+ABZ*I_ERI_Gx3y_Px_S_S_a;
  Double I_ERI_Gx2yz_Dxz_S_S_a = I_ERI_Hx2y2z_Px_S_S_a+ABZ*I_ERI_Gx2yz_Px_S_S_a;
  Double I_ERI_Gxy2z_Dxz_S_S_a = I_ERI_Hxy3z_Px_S_S_a+ABZ*I_ERI_Gxy2z_Px_S_S_a;
  Double I_ERI_Gx3z_Dxz_S_S_a = I_ERI_Hx4z_Px_S_S_a+ABZ*I_ERI_Gx3z_Px_S_S_a;
  Double I_ERI_G4y_Dxz_S_S_a = I_ERI_H4yz_Px_S_S_a+ABZ*I_ERI_G4y_Px_S_S_a;
  Double I_ERI_G3yz_Dxz_S_S_a = I_ERI_H3y2z_Px_S_S_a+ABZ*I_ERI_G3yz_Px_S_S_a;
  Double I_ERI_G2y2z_Dxz_S_S_a = I_ERI_H2y3z_Px_S_S_a+ABZ*I_ERI_G2y2z_Px_S_S_a;
  Double I_ERI_Gy3z_Dxz_S_S_a = I_ERI_Hy4z_Px_S_S_a+ABZ*I_ERI_Gy3z_Px_S_S_a;
  Double I_ERI_G4z_Dxz_S_S_a = I_ERI_H5z_Px_S_S_a+ABZ*I_ERI_G4z_Px_S_S_a;
  Double I_ERI_G4x_D2y_S_S_a = I_ERI_H4xy_Py_S_S_a+ABY*I_ERI_G4x_Py_S_S_a;
  Double I_ERI_G3xy_D2y_S_S_a = I_ERI_H3x2y_Py_S_S_a+ABY*I_ERI_G3xy_Py_S_S_a;
  Double I_ERI_G3xz_D2y_S_S_a = I_ERI_H3xyz_Py_S_S_a+ABY*I_ERI_G3xz_Py_S_S_a;
  Double I_ERI_G2x2y_D2y_S_S_a = I_ERI_H2x3y_Py_S_S_a+ABY*I_ERI_G2x2y_Py_S_S_a;
  Double I_ERI_G2xyz_D2y_S_S_a = I_ERI_H2x2yz_Py_S_S_a+ABY*I_ERI_G2xyz_Py_S_S_a;
  Double I_ERI_G2x2z_D2y_S_S_a = I_ERI_H2xy2z_Py_S_S_a+ABY*I_ERI_G2x2z_Py_S_S_a;
  Double I_ERI_Gx3y_D2y_S_S_a = I_ERI_Hx4y_Py_S_S_a+ABY*I_ERI_Gx3y_Py_S_S_a;
  Double I_ERI_Gx2yz_D2y_S_S_a = I_ERI_Hx3yz_Py_S_S_a+ABY*I_ERI_Gx2yz_Py_S_S_a;
  Double I_ERI_Gxy2z_D2y_S_S_a = I_ERI_Hx2y2z_Py_S_S_a+ABY*I_ERI_Gxy2z_Py_S_S_a;
  Double I_ERI_Gx3z_D2y_S_S_a = I_ERI_Hxy3z_Py_S_S_a+ABY*I_ERI_Gx3z_Py_S_S_a;
  Double I_ERI_G4y_D2y_S_S_a = I_ERI_H5y_Py_S_S_a+ABY*I_ERI_G4y_Py_S_S_a;
  Double I_ERI_G3yz_D2y_S_S_a = I_ERI_H4yz_Py_S_S_a+ABY*I_ERI_G3yz_Py_S_S_a;
  Double I_ERI_G2y2z_D2y_S_S_a = I_ERI_H3y2z_Py_S_S_a+ABY*I_ERI_G2y2z_Py_S_S_a;
  Double I_ERI_Gy3z_D2y_S_S_a = I_ERI_H2y3z_Py_S_S_a+ABY*I_ERI_Gy3z_Py_S_S_a;
  Double I_ERI_G4z_D2y_S_S_a = I_ERI_Hy4z_Py_S_S_a+ABY*I_ERI_G4z_Py_S_S_a;
  Double I_ERI_G4x_Dyz_S_S_a = I_ERI_H4xz_Py_S_S_a+ABZ*I_ERI_G4x_Py_S_S_a;
  Double I_ERI_G3xy_Dyz_S_S_a = I_ERI_H3xyz_Py_S_S_a+ABZ*I_ERI_G3xy_Py_S_S_a;
  Double I_ERI_G3xz_Dyz_S_S_a = I_ERI_H3x2z_Py_S_S_a+ABZ*I_ERI_G3xz_Py_S_S_a;
  Double I_ERI_G2x2y_Dyz_S_S_a = I_ERI_H2x2yz_Py_S_S_a+ABZ*I_ERI_G2x2y_Py_S_S_a;
  Double I_ERI_G2xyz_Dyz_S_S_a = I_ERI_H2xy2z_Py_S_S_a+ABZ*I_ERI_G2xyz_Py_S_S_a;
  Double I_ERI_G2x2z_Dyz_S_S_a = I_ERI_H2x3z_Py_S_S_a+ABZ*I_ERI_G2x2z_Py_S_S_a;
  Double I_ERI_Gx3y_Dyz_S_S_a = I_ERI_Hx3yz_Py_S_S_a+ABZ*I_ERI_Gx3y_Py_S_S_a;
  Double I_ERI_Gx2yz_Dyz_S_S_a = I_ERI_Hx2y2z_Py_S_S_a+ABZ*I_ERI_Gx2yz_Py_S_S_a;
  Double I_ERI_Gxy2z_Dyz_S_S_a = I_ERI_Hxy3z_Py_S_S_a+ABZ*I_ERI_Gxy2z_Py_S_S_a;
  Double I_ERI_Gx3z_Dyz_S_S_a = I_ERI_Hx4z_Py_S_S_a+ABZ*I_ERI_Gx3z_Py_S_S_a;
  Double I_ERI_G4y_Dyz_S_S_a = I_ERI_H4yz_Py_S_S_a+ABZ*I_ERI_G4y_Py_S_S_a;
  Double I_ERI_G3yz_Dyz_S_S_a = I_ERI_H3y2z_Py_S_S_a+ABZ*I_ERI_G3yz_Py_S_S_a;
  Double I_ERI_G2y2z_Dyz_S_S_a = I_ERI_H2y3z_Py_S_S_a+ABZ*I_ERI_G2y2z_Py_S_S_a;
  Double I_ERI_Gy3z_Dyz_S_S_a = I_ERI_Hy4z_Py_S_S_a+ABZ*I_ERI_Gy3z_Py_S_S_a;
  Double I_ERI_G4z_Dyz_S_S_a = I_ERI_H5z_Py_S_S_a+ABZ*I_ERI_G4z_Py_S_S_a;
  Double I_ERI_G4x_D2z_S_S_a = I_ERI_H4xz_Pz_S_S_a+ABZ*I_ERI_G4x_Pz_S_S_a;
  Double I_ERI_G3xy_D2z_S_S_a = I_ERI_H3xyz_Pz_S_S_a+ABZ*I_ERI_G3xy_Pz_S_S_a;
  Double I_ERI_G3xz_D2z_S_S_a = I_ERI_H3x2z_Pz_S_S_a+ABZ*I_ERI_G3xz_Pz_S_S_a;
  Double I_ERI_G2x2y_D2z_S_S_a = I_ERI_H2x2yz_Pz_S_S_a+ABZ*I_ERI_G2x2y_Pz_S_S_a;
  Double I_ERI_G2xyz_D2z_S_S_a = I_ERI_H2xy2z_Pz_S_S_a+ABZ*I_ERI_G2xyz_Pz_S_S_a;
  Double I_ERI_G2x2z_D2z_S_S_a = I_ERI_H2x3z_Pz_S_S_a+ABZ*I_ERI_G2x2z_Pz_S_S_a;
  Double I_ERI_Gx3y_D2z_S_S_a = I_ERI_Hx3yz_Pz_S_S_a+ABZ*I_ERI_Gx3y_Pz_S_S_a;
  Double I_ERI_Gx2yz_D2z_S_S_a = I_ERI_Hx2y2z_Pz_S_S_a+ABZ*I_ERI_Gx2yz_Pz_S_S_a;
  Double I_ERI_Gxy2z_D2z_S_S_a = I_ERI_Hxy3z_Pz_S_S_a+ABZ*I_ERI_Gxy2z_Pz_S_S_a;
  Double I_ERI_Gx3z_D2z_S_S_a = I_ERI_Hx4z_Pz_S_S_a+ABZ*I_ERI_Gx3z_Pz_S_S_a;
  Double I_ERI_G4y_D2z_S_S_a = I_ERI_H4yz_Pz_S_S_a+ABZ*I_ERI_G4y_Pz_S_S_a;
  Double I_ERI_G3yz_D2z_S_S_a = I_ERI_H3y2z_Pz_S_S_a+ABZ*I_ERI_G3yz_Pz_S_S_a;
  Double I_ERI_G2y2z_D2z_S_S_a = I_ERI_H2y3z_Pz_S_S_a+ABZ*I_ERI_G2y2z_Pz_S_S_a;
  Double I_ERI_Gy3z_D2z_S_S_a = I_ERI_Hy4z_Pz_S_S_a+ABZ*I_ERI_Gy3z_Pz_S_S_a;
  Double I_ERI_G4z_D2z_S_S_a = I_ERI_H5z_Pz_S_S_a+ABZ*I_ERI_G4z_Pz_S_S_a;

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
   * shell quartet name: SQ_ERI_G_P_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
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
  Double I_ERI_G4x_Py_S_S_b = I_ERI_H4xy_S_S_S_b+ABY*I_ERI_G4x_S_S_S_b;
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
  Double I_ERI_G4x_Pz_S_S_b = I_ERI_H4xz_S_S_S_b+ABZ*I_ERI_G4x_S_S_S_b;
  Double I_ERI_G3xy_Pz_S_S_b = I_ERI_H3xyz_S_S_S_b+ABZ*I_ERI_G3xy_S_S_S_b;
  Double I_ERI_G3xz_Pz_S_S_b = I_ERI_H3x2z_S_S_S_b+ABZ*I_ERI_G3xz_S_S_S_b;
  Double I_ERI_G2x2y_Pz_S_S_b = I_ERI_H2x2yz_S_S_S_b+ABZ*I_ERI_G2x2y_S_S_S_b;
  Double I_ERI_G2xyz_Pz_S_S_b = I_ERI_H2xy2z_S_S_S_b+ABZ*I_ERI_G2xyz_S_S_S_b;
  Double I_ERI_G2x2z_Pz_S_S_b = I_ERI_H2x3z_S_S_S_b+ABZ*I_ERI_G2x2z_S_S_S_b;
  Double I_ERI_Gx3y_Pz_S_S_b = I_ERI_Hx3yz_S_S_S_b+ABZ*I_ERI_Gx3y_S_S_S_b;
  Double I_ERI_Gx2yz_Pz_S_S_b = I_ERI_Hx2y2z_S_S_S_b+ABZ*I_ERI_Gx2yz_S_S_S_b;
  Double I_ERI_Gxy2z_Pz_S_S_b = I_ERI_Hxy3z_S_S_S_b+ABZ*I_ERI_Gxy2z_S_S_S_b;
  Double I_ERI_Gx3z_Pz_S_S_b = I_ERI_Hx4z_S_S_S_b+ABZ*I_ERI_Gx3z_S_S_S_b;
  Double I_ERI_G4y_Pz_S_S_b = I_ERI_H4yz_S_S_S_b+ABZ*I_ERI_G4y_S_S_S_b;
  Double I_ERI_G3yz_Pz_S_S_b = I_ERI_H3y2z_S_S_S_b+ABZ*I_ERI_G3yz_S_S_S_b;
  Double I_ERI_G2y2z_Pz_S_S_b = I_ERI_H2y3z_S_S_S_b+ABZ*I_ERI_G2y2z_S_S_S_b;
  Double I_ERI_Gy3z_Pz_S_S_b = I_ERI_Hy4z_S_S_S_b+ABZ*I_ERI_Gy3z_S_S_S_b;
  Double I_ERI_G4z_Pz_S_S_b = I_ERI_H5z_S_S_S_b+ABZ*I_ERI_G4z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
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
   * shell quartet name: SQ_ERI_H_P_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 14 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_S_S_b
   * RHS shell quartet name: SQ_ERI_H_S_S_S_b
   ************************************************************/
  Double I_ERI_H5x_Px_S_S_b = I_ERI_I6x_S_S_S_b+ABX*I_ERI_H5x_S_S_S_b;
  Double I_ERI_H4xy_Px_S_S_b = I_ERI_I5xy_S_S_S_b+ABX*I_ERI_H4xy_S_S_S_b;
  Double I_ERI_H4xz_Px_S_S_b = I_ERI_I5xz_S_S_S_b+ABX*I_ERI_H4xz_S_S_S_b;
  Double I_ERI_H3x2y_Px_S_S_b = I_ERI_I4x2y_S_S_S_b+ABX*I_ERI_H3x2y_S_S_S_b;
  Double I_ERI_H3xyz_Px_S_S_b = I_ERI_I4xyz_S_S_S_b+ABX*I_ERI_H3xyz_S_S_S_b;
  Double I_ERI_H3x2z_Px_S_S_b = I_ERI_I4x2z_S_S_S_b+ABX*I_ERI_H3x2z_S_S_S_b;
  Double I_ERI_H2x3y_Px_S_S_b = I_ERI_I3x3y_S_S_S_b+ABX*I_ERI_H2x3y_S_S_S_b;
  Double I_ERI_H2x2yz_Px_S_S_b = I_ERI_I3x2yz_S_S_S_b+ABX*I_ERI_H2x2yz_S_S_S_b;
  Double I_ERI_H2xy2z_Px_S_S_b = I_ERI_I3xy2z_S_S_S_b+ABX*I_ERI_H2xy2z_S_S_S_b;
  Double I_ERI_H2x3z_Px_S_S_b = I_ERI_I3x3z_S_S_S_b+ABX*I_ERI_H2x3z_S_S_S_b;
  Double I_ERI_Hx4y_Px_S_S_b = I_ERI_I2x4y_S_S_S_b+ABX*I_ERI_Hx4y_S_S_S_b;
  Double I_ERI_Hx3yz_Px_S_S_b = I_ERI_I2x3yz_S_S_S_b+ABX*I_ERI_Hx3yz_S_S_S_b;
  Double I_ERI_Hx2y2z_Px_S_S_b = I_ERI_I2x2y2z_S_S_S_b+ABX*I_ERI_Hx2y2z_S_S_S_b;
  Double I_ERI_Hxy3z_Px_S_S_b = I_ERI_I2xy3z_S_S_S_b+ABX*I_ERI_Hxy3z_S_S_S_b;
  Double I_ERI_Hx4z_Px_S_S_b = I_ERI_I2x4z_S_S_S_b+ABX*I_ERI_Hx4z_S_S_S_b;
  Double I_ERI_H4yz_Px_S_S_b = I_ERI_Ix4yz_S_S_S_b+ABX*I_ERI_H4yz_S_S_S_b;
  Double I_ERI_H3y2z_Px_S_S_b = I_ERI_Ix3y2z_S_S_S_b+ABX*I_ERI_H3y2z_S_S_S_b;
  Double I_ERI_H2y3z_Px_S_S_b = I_ERI_Ix2y3z_S_S_S_b+ABX*I_ERI_H2y3z_S_S_S_b;
  Double I_ERI_Hy4z_Px_S_S_b = I_ERI_Ixy4z_S_S_S_b+ABX*I_ERI_Hy4z_S_S_S_b;
  Double I_ERI_H4xy_Py_S_S_b = I_ERI_I4x2y_S_S_S_b+ABY*I_ERI_H4xy_S_S_S_b;
  Double I_ERI_H3x2y_Py_S_S_b = I_ERI_I3x3y_S_S_S_b+ABY*I_ERI_H3x2y_S_S_S_b;
  Double I_ERI_H3xyz_Py_S_S_b = I_ERI_I3x2yz_S_S_S_b+ABY*I_ERI_H3xyz_S_S_S_b;
  Double I_ERI_H2x3y_Py_S_S_b = I_ERI_I2x4y_S_S_S_b+ABY*I_ERI_H2x3y_S_S_S_b;
  Double I_ERI_H2x2yz_Py_S_S_b = I_ERI_I2x3yz_S_S_S_b+ABY*I_ERI_H2x2yz_S_S_S_b;
  Double I_ERI_H2xy2z_Py_S_S_b = I_ERI_I2x2y2z_S_S_S_b+ABY*I_ERI_H2xy2z_S_S_S_b;
  Double I_ERI_Hx4y_Py_S_S_b = I_ERI_Ix5y_S_S_S_b+ABY*I_ERI_Hx4y_S_S_S_b;
  Double I_ERI_Hx3yz_Py_S_S_b = I_ERI_Ix4yz_S_S_S_b+ABY*I_ERI_Hx3yz_S_S_S_b;
  Double I_ERI_Hx2y2z_Py_S_S_b = I_ERI_Ix3y2z_S_S_S_b+ABY*I_ERI_Hx2y2z_S_S_S_b;
  Double I_ERI_Hxy3z_Py_S_S_b = I_ERI_Ix2y3z_S_S_S_b+ABY*I_ERI_Hxy3z_S_S_S_b;
  Double I_ERI_H5y_Py_S_S_b = I_ERI_I6y_S_S_S_b+ABY*I_ERI_H5y_S_S_S_b;
  Double I_ERI_H4yz_Py_S_S_b = I_ERI_I5yz_S_S_S_b+ABY*I_ERI_H4yz_S_S_S_b;
  Double I_ERI_H3y2z_Py_S_S_b = I_ERI_I4y2z_S_S_S_b+ABY*I_ERI_H3y2z_S_S_S_b;
  Double I_ERI_H2y3z_Py_S_S_b = I_ERI_I3y3z_S_S_S_b+ABY*I_ERI_H2y3z_S_S_S_b;
  Double I_ERI_Hy4z_Py_S_S_b = I_ERI_I2y4z_S_S_S_b+ABY*I_ERI_Hy4z_S_S_S_b;
  Double I_ERI_H4xz_Pz_S_S_b = I_ERI_I4x2z_S_S_S_b+ABZ*I_ERI_H4xz_S_S_S_b;
  Double I_ERI_H3xyz_Pz_S_S_b = I_ERI_I3xy2z_S_S_S_b+ABZ*I_ERI_H3xyz_S_S_S_b;
  Double I_ERI_H3x2z_Pz_S_S_b = I_ERI_I3x3z_S_S_S_b+ABZ*I_ERI_H3x2z_S_S_S_b;
  Double I_ERI_H2x2yz_Pz_S_S_b = I_ERI_I2x2y2z_S_S_S_b+ABZ*I_ERI_H2x2yz_S_S_S_b;
  Double I_ERI_H2xy2z_Pz_S_S_b = I_ERI_I2xy3z_S_S_S_b+ABZ*I_ERI_H2xy2z_S_S_S_b;
  Double I_ERI_H2x3z_Pz_S_S_b = I_ERI_I2x4z_S_S_S_b+ABZ*I_ERI_H2x3z_S_S_S_b;
  Double I_ERI_Hx3yz_Pz_S_S_b = I_ERI_Ix3y2z_S_S_S_b+ABZ*I_ERI_Hx3yz_S_S_S_b;
  Double I_ERI_Hx2y2z_Pz_S_S_b = I_ERI_Ix2y3z_S_S_S_b+ABZ*I_ERI_Hx2y2z_S_S_S_b;
  Double I_ERI_Hxy3z_Pz_S_S_b = I_ERI_Ixy4z_S_S_S_b+ABZ*I_ERI_Hxy3z_S_S_S_b;
  Double I_ERI_Hx4z_Pz_S_S_b = I_ERI_Ix5z_S_S_S_b+ABZ*I_ERI_Hx4z_S_S_S_b;
  Double I_ERI_H4yz_Pz_S_S_b = I_ERI_I4y2z_S_S_S_b+ABZ*I_ERI_H4yz_S_S_S_b;
  Double I_ERI_H3y2z_Pz_S_S_b = I_ERI_I3y3z_S_S_S_b+ABZ*I_ERI_H3y2z_S_S_S_b;
  Double I_ERI_H2y3z_Pz_S_S_b = I_ERI_I2y4z_S_S_S_b+ABZ*I_ERI_H2y3z_S_S_S_b;
  Double I_ERI_Hy4z_Pz_S_S_b = I_ERI_Iy5z_S_S_S_b+ABZ*I_ERI_Hy4z_S_S_S_b;
  Double I_ERI_H5z_Pz_S_S_b = I_ERI_I6z_S_S_S_b+ABZ*I_ERI_H5z_S_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_G_D_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 35 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_P_S_S_b
   * RHS shell quartet name: SQ_ERI_G_P_S_S_b
   ************************************************************/
  Double I_ERI_G4x_D2x_S_S_b = I_ERI_H5x_Px_S_S_b+ABX*I_ERI_G4x_Px_S_S_b;
  Double I_ERI_G3xy_D2x_S_S_b = I_ERI_H4xy_Px_S_S_b+ABX*I_ERI_G3xy_Px_S_S_b;
  Double I_ERI_G3xz_D2x_S_S_b = I_ERI_H4xz_Px_S_S_b+ABX*I_ERI_G3xz_Px_S_S_b;
  Double I_ERI_G2x2y_D2x_S_S_b = I_ERI_H3x2y_Px_S_S_b+ABX*I_ERI_G2x2y_Px_S_S_b;
  Double I_ERI_G2xyz_D2x_S_S_b = I_ERI_H3xyz_Px_S_S_b+ABX*I_ERI_G2xyz_Px_S_S_b;
  Double I_ERI_G2x2z_D2x_S_S_b = I_ERI_H3x2z_Px_S_S_b+ABX*I_ERI_G2x2z_Px_S_S_b;
  Double I_ERI_Gx3y_D2x_S_S_b = I_ERI_H2x3y_Px_S_S_b+ABX*I_ERI_Gx3y_Px_S_S_b;
  Double I_ERI_Gx2yz_D2x_S_S_b = I_ERI_H2x2yz_Px_S_S_b+ABX*I_ERI_Gx2yz_Px_S_S_b;
  Double I_ERI_Gxy2z_D2x_S_S_b = I_ERI_H2xy2z_Px_S_S_b+ABX*I_ERI_Gxy2z_Px_S_S_b;
  Double I_ERI_Gx3z_D2x_S_S_b = I_ERI_H2x3z_Px_S_S_b+ABX*I_ERI_Gx3z_Px_S_S_b;
  Double I_ERI_G4y_D2x_S_S_b = I_ERI_Hx4y_Px_S_S_b+ABX*I_ERI_G4y_Px_S_S_b;
  Double I_ERI_G3yz_D2x_S_S_b = I_ERI_Hx3yz_Px_S_S_b+ABX*I_ERI_G3yz_Px_S_S_b;
  Double I_ERI_G2y2z_D2x_S_S_b = I_ERI_Hx2y2z_Px_S_S_b+ABX*I_ERI_G2y2z_Px_S_S_b;
  Double I_ERI_Gy3z_D2x_S_S_b = I_ERI_Hxy3z_Px_S_S_b+ABX*I_ERI_Gy3z_Px_S_S_b;
  Double I_ERI_G4z_D2x_S_S_b = I_ERI_Hx4z_Px_S_S_b+ABX*I_ERI_G4z_Px_S_S_b;
  Double I_ERI_G3xz_Dxy_S_S_b = I_ERI_H3xyz_Px_S_S_b+ABY*I_ERI_G3xz_Px_S_S_b;
  Double I_ERI_G2xyz_Dxy_S_S_b = I_ERI_H2x2yz_Px_S_S_b+ABY*I_ERI_G2xyz_Px_S_S_b;
  Double I_ERI_G2x2z_Dxy_S_S_b = I_ERI_H2xy2z_Px_S_S_b+ABY*I_ERI_G2x2z_Px_S_S_b;
  Double I_ERI_Gx2yz_Dxy_S_S_b = I_ERI_Hx3yz_Px_S_S_b+ABY*I_ERI_Gx2yz_Px_S_S_b;
  Double I_ERI_Gxy2z_Dxy_S_S_b = I_ERI_Hx2y2z_Px_S_S_b+ABY*I_ERI_Gxy2z_Px_S_S_b;
  Double I_ERI_Gx3z_Dxy_S_S_b = I_ERI_Hxy3z_Px_S_S_b+ABY*I_ERI_Gx3z_Px_S_S_b;
  Double I_ERI_G3yz_Dxy_S_S_b = I_ERI_H4yz_Px_S_S_b+ABY*I_ERI_G3yz_Px_S_S_b;
  Double I_ERI_G2y2z_Dxy_S_S_b = I_ERI_H3y2z_Px_S_S_b+ABY*I_ERI_G2y2z_Px_S_S_b;
  Double I_ERI_Gy3z_Dxy_S_S_b = I_ERI_H2y3z_Px_S_S_b+ABY*I_ERI_Gy3z_Px_S_S_b;
  Double I_ERI_G4z_Dxy_S_S_b = I_ERI_Hy4z_Px_S_S_b+ABY*I_ERI_G4z_Px_S_S_b;
  Double I_ERI_G4x_D2y_S_S_b = I_ERI_H4xy_Py_S_S_b+ABY*I_ERI_G4x_Py_S_S_b;
  Double I_ERI_G3xy_D2y_S_S_b = I_ERI_H3x2y_Py_S_S_b+ABY*I_ERI_G3xy_Py_S_S_b;
  Double I_ERI_G3xz_D2y_S_S_b = I_ERI_H3xyz_Py_S_S_b+ABY*I_ERI_G3xz_Py_S_S_b;
  Double I_ERI_G2x2y_D2y_S_S_b = I_ERI_H2x3y_Py_S_S_b+ABY*I_ERI_G2x2y_Py_S_S_b;
  Double I_ERI_G2xyz_D2y_S_S_b = I_ERI_H2x2yz_Py_S_S_b+ABY*I_ERI_G2xyz_Py_S_S_b;
  Double I_ERI_G2x2z_D2y_S_S_b = I_ERI_H2xy2z_Py_S_S_b+ABY*I_ERI_G2x2z_Py_S_S_b;
  Double I_ERI_Gx3y_D2y_S_S_b = I_ERI_Hx4y_Py_S_S_b+ABY*I_ERI_Gx3y_Py_S_S_b;
  Double I_ERI_Gx2yz_D2y_S_S_b = I_ERI_Hx3yz_Py_S_S_b+ABY*I_ERI_Gx2yz_Py_S_S_b;
  Double I_ERI_Gxy2z_D2y_S_S_b = I_ERI_Hx2y2z_Py_S_S_b+ABY*I_ERI_Gxy2z_Py_S_S_b;
  Double I_ERI_Gx3z_D2y_S_S_b = I_ERI_Hxy3z_Py_S_S_b+ABY*I_ERI_Gx3z_Py_S_S_b;
  Double I_ERI_G4y_D2y_S_S_b = I_ERI_H5y_Py_S_S_b+ABY*I_ERI_G4y_Py_S_S_b;
  Double I_ERI_G3yz_D2y_S_S_b = I_ERI_H4yz_Py_S_S_b+ABY*I_ERI_G3yz_Py_S_S_b;
  Double I_ERI_G2y2z_D2y_S_S_b = I_ERI_H3y2z_Py_S_S_b+ABY*I_ERI_G2y2z_Py_S_S_b;
  Double I_ERI_Gy3z_D2y_S_S_b = I_ERI_H2y3z_Py_S_S_b+ABY*I_ERI_Gy3z_Py_S_S_b;
  Double I_ERI_G4z_D2y_S_S_b = I_ERI_Hy4z_Py_S_S_b+ABY*I_ERI_G4z_Py_S_S_b;
  Double I_ERI_G4x_D2z_S_S_b = I_ERI_H4xz_Pz_S_S_b+ABZ*I_ERI_G4x_Pz_S_S_b;
  Double I_ERI_G3xy_D2z_S_S_b = I_ERI_H3xyz_Pz_S_S_b+ABZ*I_ERI_G3xy_Pz_S_S_b;
  Double I_ERI_G3xz_D2z_S_S_b = I_ERI_H3x2z_Pz_S_S_b+ABZ*I_ERI_G3xz_Pz_S_S_b;
  Double I_ERI_G2x2y_D2z_S_S_b = I_ERI_H2x2yz_Pz_S_S_b+ABZ*I_ERI_G2x2y_Pz_S_S_b;
  Double I_ERI_G2xyz_D2z_S_S_b = I_ERI_H2xy2z_Pz_S_S_b+ABZ*I_ERI_G2xyz_Pz_S_S_b;
  Double I_ERI_G2x2z_D2z_S_S_b = I_ERI_H2x3z_Pz_S_S_b+ABZ*I_ERI_G2x2z_Pz_S_S_b;
  Double I_ERI_Gx3y_D2z_S_S_b = I_ERI_Hx3yz_Pz_S_S_b+ABZ*I_ERI_Gx3y_Pz_S_S_b;
  Double I_ERI_Gx2yz_D2z_S_S_b = I_ERI_Hx2y2z_Pz_S_S_b+ABZ*I_ERI_Gx2yz_Pz_S_S_b;
  Double I_ERI_Gxy2z_D2z_S_S_b = I_ERI_Hxy3z_Pz_S_S_b+ABZ*I_ERI_Gxy2z_Pz_S_S_b;
  Double I_ERI_Gx3z_D2z_S_S_b = I_ERI_Hx4z_Pz_S_S_b+ABZ*I_ERI_Gx3z_Pz_S_S_b;
  Double I_ERI_G4y_D2z_S_S_b = I_ERI_H4yz_Pz_S_S_b+ABZ*I_ERI_G4y_Pz_S_S_b;
  Double I_ERI_G3yz_D2z_S_S_b = I_ERI_H3y2z_Pz_S_S_b+ABZ*I_ERI_G3yz_Pz_S_S_b;
  Double I_ERI_G2y2z_D2z_S_S_b = I_ERI_H2y3z_Pz_S_S_b+ABZ*I_ERI_G2y2z_Pz_S_S_b;
  Double I_ERI_Gy3z_D2z_S_S_b = I_ERI_Hy4z_Pz_S_S_b+ABZ*I_ERI_Gy3z_Pz_S_S_b;
  Double I_ERI_G4z_D2z_S_S_b = I_ERI_H5z_Pz_S_S_b+ABZ*I_ERI_G4z_Pz_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_F_S_S_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_S_S_b
   * RHS shell quartet name: SQ_ERI_F_D_S_S_b
   ************************************************************/
  Double I_ERI_F3x_F3x_S_S_b = I_ERI_G4x_D2x_S_S_b+ABX*I_ERI_F3x_D2x_S_S_b;
  Double I_ERI_F2xy_F3x_S_S_b = I_ERI_G3xy_D2x_S_S_b+ABX*I_ERI_F2xy_D2x_S_S_b;
  Double I_ERI_F2xz_F3x_S_S_b = I_ERI_G3xz_D2x_S_S_b+ABX*I_ERI_F2xz_D2x_S_S_b;
  Double I_ERI_Fx2y_F3x_S_S_b = I_ERI_G2x2y_D2x_S_S_b+ABX*I_ERI_Fx2y_D2x_S_S_b;
  Double I_ERI_Fxyz_F3x_S_S_b = I_ERI_G2xyz_D2x_S_S_b+ABX*I_ERI_Fxyz_D2x_S_S_b;
  Double I_ERI_Fx2z_F3x_S_S_b = I_ERI_G2x2z_D2x_S_S_b+ABX*I_ERI_Fx2z_D2x_S_S_b;
  Double I_ERI_F3y_F3x_S_S_b = I_ERI_Gx3y_D2x_S_S_b+ABX*I_ERI_F3y_D2x_S_S_b;
  Double I_ERI_F2yz_F3x_S_S_b = I_ERI_Gx2yz_D2x_S_S_b+ABX*I_ERI_F2yz_D2x_S_S_b;
  Double I_ERI_Fy2z_F3x_S_S_b = I_ERI_Gxy2z_D2x_S_S_b+ABX*I_ERI_Fy2z_D2x_S_S_b;
  Double I_ERI_F3z_F3x_S_S_b = I_ERI_Gx3z_D2x_S_S_b+ABX*I_ERI_F3z_D2x_S_S_b;
  Double I_ERI_F3x_F2xy_S_S_b = I_ERI_G3xy_D2x_S_S_b+ABY*I_ERI_F3x_D2x_S_S_b;
  Double I_ERI_F2xy_F2xy_S_S_b = I_ERI_G2x2y_D2x_S_S_b+ABY*I_ERI_F2xy_D2x_S_S_b;
  Double I_ERI_F2xz_F2xy_S_S_b = I_ERI_G2xyz_D2x_S_S_b+ABY*I_ERI_F2xz_D2x_S_S_b;
  Double I_ERI_Fx2y_F2xy_S_S_b = I_ERI_Gx3y_D2x_S_S_b+ABY*I_ERI_Fx2y_D2x_S_S_b;
  Double I_ERI_Fxyz_F2xy_S_S_b = I_ERI_Gx2yz_D2x_S_S_b+ABY*I_ERI_Fxyz_D2x_S_S_b;
  Double I_ERI_Fx2z_F2xy_S_S_b = I_ERI_Gxy2z_D2x_S_S_b+ABY*I_ERI_Fx2z_D2x_S_S_b;
  Double I_ERI_F3y_F2xy_S_S_b = I_ERI_G4y_D2x_S_S_b+ABY*I_ERI_F3y_D2x_S_S_b;
  Double I_ERI_F2yz_F2xy_S_S_b = I_ERI_G3yz_D2x_S_S_b+ABY*I_ERI_F2yz_D2x_S_S_b;
  Double I_ERI_Fy2z_F2xy_S_S_b = I_ERI_G2y2z_D2x_S_S_b+ABY*I_ERI_Fy2z_D2x_S_S_b;
  Double I_ERI_F3z_F2xy_S_S_b = I_ERI_Gy3z_D2x_S_S_b+ABY*I_ERI_F3z_D2x_S_S_b;
  Double I_ERI_F3x_F2xz_S_S_b = I_ERI_G3xz_D2x_S_S_b+ABZ*I_ERI_F3x_D2x_S_S_b;
  Double I_ERI_F2xy_F2xz_S_S_b = I_ERI_G2xyz_D2x_S_S_b+ABZ*I_ERI_F2xy_D2x_S_S_b;
  Double I_ERI_F2xz_F2xz_S_S_b = I_ERI_G2x2z_D2x_S_S_b+ABZ*I_ERI_F2xz_D2x_S_S_b;
  Double I_ERI_Fx2y_F2xz_S_S_b = I_ERI_Gx2yz_D2x_S_S_b+ABZ*I_ERI_Fx2y_D2x_S_S_b;
  Double I_ERI_Fxyz_F2xz_S_S_b = I_ERI_Gxy2z_D2x_S_S_b+ABZ*I_ERI_Fxyz_D2x_S_S_b;
  Double I_ERI_Fx2z_F2xz_S_S_b = I_ERI_Gx3z_D2x_S_S_b+ABZ*I_ERI_Fx2z_D2x_S_S_b;
  Double I_ERI_F3y_F2xz_S_S_b = I_ERI_G3yz_D2x_S_S_b+ABZ*I_ERI_F3y_D2x_S_S_b;
  Double I_ERI_F2yz_F2xz_S_S_b = I_ERI_G2y2z_D2x_S_S_b+ABZ*I_ERI_F2yz_D2x_S_S_b;
  Double I_ERI_Fy2z_F2xz_S_S_b = I_ERI_Gy3z_D2x_S_S_b+ABZ*I_ERI_Fy2z_D2x_S_S_b;
  Double I_ERI_F3z_F2xz_S_S_b = I_ERI_G4z_D2x_S_S_b+ABZ*I_ERI_F3z_D2x_S_S_b;
  Double I_ERI_F3x_Fx2y_S_S_b = I_ERI_G4x_D2y_S_S_b+ABX*I_ERI_F3x_D2y_S_S_b;
  Double I_ERI_F2xy_Fx2y_S_S_b = I_ERI_G3xy_D2y_S_S_b+ABX*I_ERI_F2xy_D2y_S_S_b;
  Double I_ERI_F2xz_Fx2y_S_S_b = I_ERI_G3xz_D2y_S_S_b+ABX*I_ERI_F2xz_D2y_S_S_b;
  Double I_ERI_Fx2y_Fx2y_S_S_b = I_ERI_G2x2y_D2y_S_S_b+ABX*I_ERI_Fx2y_D2y_S_S_b;
  Double I_ERI_Fxyz_Fx2y_S_S_b = I_ERI_G2xyz_D2y_S_S_b+ABX*I_ERI_Fxyz_D2y_S_S_b;
  Double I_ERI_Fx2z_Fx2y_S_S_b = I_ERI_G2x2z_D2y_S_S_b+ABX*I_ERI_Fx2z_D2y_S_S_b;
  Double I_ERI_F3y_Fx2y_S_S_b = I_ERI_Gx3y_D2y_S_S_b+ABX*I_ERI_F3y_D2y_S_S_b;
  Double I_ERI_F2yz_Fx2y_S_S_b = I_ERI_Gx2yz_D2y_S_S_b+ABX*I_ERI_F2yz_D2y_S_S_b;
  Double I_ERI_Fy2z_Fx2y_S_S_b = I_ERI_Gxy2z_D2y_S_S_b+ABX*I_ERI_Fy2z_D2y_S_S_b;
  Double I_ERI_F3z_Fx2y_S_S_b = I_ERI_Gx3z_D2y_S_S_b+ABX*I_ERI_F3z_D2y_S_S_b;
  Double I_ERI_F3x_Fxyz_S_S_b = I_ERI_G3xz_Dxy_S_S_b+ABZ*I_ERI_F3x_Dxy_S_S_b;
  Double I_ERI_F2xy_Fxyz_S_S_b = I_ERI_G2xyz_Dxy_S_S_b+ABZ*I_ERI_F2xy_Dxy_S_S_b;
  Double I_ERI_F2xz_Fxyz_S_S_b = I_ERI_G2x2z_Dxy_S_S_b+ABZ*I_ERI_F2xz_Dxy_S_S_b;
  Double I_ERI_Fx2y_Fxyz_S_S_b = I_ERI_Gx2yz_Dxy_S_S_b+ABZ*I_ERI_Fx2y_Dxy_S_S_b;
  Double I_ERI_Fxyz_Fxyz_S_S_b = I_ERI_Gxy2z_Dxy_S_S_b+ABZ*I_ERI_Fxyz_Dxy_S_S_b;
  Double I_ERI_Fx2z_Fxyz_S_S_b = I_ERI_Gx3z_Dxy_S_S_b+ABZ*I_ERI_Fx2z_Dxy_S_S_b;
  Double I_ERI_F3y_Fxyz_S_S_b = I_ERI_G3yz_Dxy_S_S_b+ABZ*I_ERI_F3y_Dxy_S_S_b;
  Double I_ERI_F2yz_Fxyz_S_S_b = I_ERI_G2y2z_Dxy_S_S_b+ABZ*I_ERI_F2yz_Dxy_S_S_b;
  Double I_ERI_Fy2z_Fxyz_S_S_b = I_ERI_Gy3z_Dxy_S_S_b+ABZ*I_ERI_Fy2z_Dxy_S_S_b;
  Double I_ERI_F3z_Fxyz_S_S_b = I_ERI_G4z_Dxy_S_S_b+ABZ*I_ERI_F3z_Dxy_S_S_b;
  Double I_ERI_F3x_Fx2z_S_S_b = I_ERI_G4x_D2z_S_S_b+ABX*I_ERI_F3x_D2z_S_S_b;
  Double I_ERI_F2xy_Fx2z_S_S_b = I_ERI_G3xy_D2z_S_S_b+ABX*I_ERI_F2xy_D2z_S_S_b;
  Double I_ERI_F2xz_Fx2z_S_S_b = I_ERI_G3xz_D2z_S_S_b+ABX*I_ERI_F2xz_D2z_S_S_b;
  Double I_ERI_Fx2y_Fx2z_S_S_b = I_ERI_G2x2y_D2z_S_S_b+ABX*I_ERI_Fx2y_D2z_S_S_b;
  Double I_ERI_Fxyz_Fx2z_S_S_b = I_ERI_G2xyz_D2z_S_S_b+ABX*I_ERI_Fxyz_D2z_S_S_b;
  Double I_ERI_Fx2z_Fx2z_S_S_b = I_ERI_G2x2z_D2z_S_S_b+ABX*I_ERI_Fx2z_D2z_S_S_b;
  Double I_ERI_F3y_Fx2z_S_S_b = I_ERI_Gx3y_D2z_S_S_b+ABX*I_ERI_F3y_D2z_S_S_b;
  Double I_ERI_F2yz_Fx2z_S_S_b = I_ERI_Gx2yz_D2z_S_S_b+ABX*I_ERI_F2yz_D2z_S_S_b;
  Double I_ERI_Fy2z_Fx2z_S_S_b = I_ERI_Gxy2z_D2z_S_S_b+ABX*I_ERI_Fy2z_D2z_S_S_b;
  Double I_ERI_F3z_Fx2z_S_S_b = I_ERI_Gx3z_D2z_S_S_b+ABX*I_ERI_F3z_D2z_S_S_b;
  Double I_ERI_F3x_F3y_S_S_b = I_ERI_G3xy_D2y_S_S_b+ABY*I_ERI_F3x_D2y_S_S_b;
  Double I_ERI_F2xy_F3y_S_S_b = I_ERI_G2x2y_D2y_S_S_b+ABY*I_ERI_F2xy_D2y_S_S_b;
  Double I_ERI_F2xz_F3y_S_S_b = I_ERI_G2xyz_D2y_S_S_b+ABY*I_ERI_F2xz_D2y_S_S_b;
  Double I_ERI_Fx2y_F3y_S_S_b = I_ERI_Gx3y_D2y_S_S_b+ABY*I_ERI_Fx2y_D2y_S_S_b;
  Double I_ERI_Fxyz_F3y_S_S_b = I_ERI_Gx2yz_D2y_S_S_b+ABY*I_ERI_Fxyz_D2y_S_S_b;
  Double I_ERI_Fx2z_F3y_S_S_b = I_ERI_Gxy2z_D2y_S_S_b+ABY*I_ERI_Fx2z_D2y_S_S_b;
  Double I_ERI_F3y_F3y_S_S_b = I_ERI_G4y_D2y_S_S_b+ABY*I_ERI_F3y_D2y_S_S_b;
  Double I_ERI_F2yz_F3y_S_S_b = I_ERI_G3yz_D2y_S_S_b+ABY*I_ERI_F2yz_D2y_S_S_b;
  Double I_ERI_Fy2z_F3y_S_S_b = I_ERI_G2y2z_D2y_S_S_b+ABY*I_ERI_Fy2z_D2y_S_S_b;
  Double I_ERI_F3z_F3y_S_S_b = I_ERI_Gy3z_D2y_S_S_b+ABY*I_ERI_F3z_D2y_S_S_b;
  Double I_ERI_F3x_F2yz_S_S_b = I_ERI_G3xz_D2y_S_S_b+ABZ*I_ERI_F3x_D2y_S_S_b;
  Double I_ERI_F2xy_F2yz_S_S_b = I_ERI_G2xyz_D2y_S_S_b+ABZ*I_ERI_F2xy_D2y_S_S_b;
  Double I_ERI_F2xz_F2yz_S_S_b = I_ERI_G2x2z_D2y_S_S_b+ABZ*I_ERI_F2xz_D2y_S_S_b;
  Double I_ERI_Fx2y_F2yz_S_S_b = I_ERI_Gx2yz_D2y_S_S_b+ABZ*I_ERI_Fx2y_D2y_S_S_b;
  Double I_ERI_Fxyz_F2yz_S_S_b = I_ERI_Gxy2z_D2y_S_S_b+ABZ*I_ERI_Fxyz_D2y_S_S_b;
  Double I_ERI_Fx2z_F2yz_S_S_b = I_ERI_Gx3z_D2y_S_S_b+ABZ*I_ERI_Fx2z_D2y_S_S_b;
  Double I_ERI_F3y_F2yz_S_S_b = I_ERI_G3yz_D2y_S_S_b+ABZ*I_ERI_F3y_D2y_S_S_b;
  Double I_ERI_F2yz_F2yz_S_S_b = I_ERI_G2y2z_D2y_S_S_b+ABZ*I_ERI_F2yz_D2y_S_S_b;
  Double I_ERI_Fy2z_F2yz_S_S_b = I_ERI_Gy3z_D2y_S_S_b+ABZ*I_ERI_Fy2z_D2y_S_S_b;
  Double I_ERI_F3z_F2yz_S_S_b = I_ERI_G4z_D2y_S_S_b+ABZ*I_ERI_F3z_D2y_S_S_b;
  Double I_ERI_F3x_Fy2z_S_S_b = I_ERI_G3xy_D2z_S_S_b+ABY*I_ERI_F3x_D2z_S_S_b;
  Double I_ERI_F2xy_Fy2z_S_S_b = I_ERI_G2x2y_D2z_S_S_b+ABY*I_ERI_F2xy_D2z_S_S_b;
  Double I_ERI_F2xz_Fy2z_S_S_b = I_ERI_G2xyz_D2z_S_S_b+ABY*I_ERI_F2xz_D2z_S_S_b;
  Double I_ERI_Fx2y_Fy2z_S_S_b = I_ERI_Gx3y_D2z_S_S_b+ABY*I_ERI_Fx2y_D2z_S_S_b;
  Double I_ERI_Fxyz_Fy2z_S_S_b = I_ERI_Gx2yz_D2z_S_S_b+ABY*I_ERI_Fxyz_D2z_S_S_b;
  Double I_ERI_Fx2z_Fy2z_S_S_b = I_ERI_Gxy2z_D2z_S_S_b+ABY*I_ERI_Fx2z_D2z_S_S_b;
  Double I_ERI_F3y_Fy2z_S_S_b = I_ERI_G4y_D2z_S_S_b+ABY*I_ERI_F3y_D2z_S_S_b;
  Double I_ERI_F2yz_Fy2z_S_S_b = I_ERI_G3yz_D2z_S_S_b+ABY*I_ERI_F2yz_D2z_S_S_b;
  Double I_ERI_Fy2z_Fy2z_S_S_b = I_ERI_G2y2z_D2z_S_S_b+ABY*I_ERI_Fy2z_D2z_S_S_b;
  Double I_ERI_F3z_Fy2z_S_S_b = I_ERI_Gy3z_D2z_S_S_b+ABY*I_ERI_F3z_D2z_S_S_b;
  Double I_ERI_F3x_F3z_S_S_b = I_ERI_G3xz_D2z_S_S_b+ABZ*I_ERI_F3x_D2z_S_S_b;
  Double I_ERI_F2xy_F3z_S_S_b = I_ERI_G2xyz_D2z_S_S_b+ABZ*I_ERI_F2xy_D2z_S_S_b;
  Double I_ERI_F2xz_F3z_S_S_b = I_ERI_G2x2z_D2z_S_S_b+ABZ*I_ERI_F2xz_D2z_S_S_b;
  Double I_ERI_Fx2y_F3z_S_S_b = I_ERI_Gx2yz_D2z_S_S_b+ABZ*I_ERI_Fx2y_D2z_S_S_b;
  Double I_ERI_Fxyz_F3z_S_S_b = I_ERI_Gxy2z_D2z_S_S_b+ABZ*I_ERI_Fxyz_D2z_S_S_b;
  Double I_ERI_Fx2z_F3z_S_S_b = I_ERI_Gx3z_D2z_S_S_b+ABZ*I_ERI_Fx2z_D2z_S_S_b;
  Double I_ERI_F3y_F3z_S_S_b = I_ERI_G3yz_D2z_S_S_b+ABZ*I_ERI_F3y_D2z_S_S_b;
  Double I_ERI_F2yz_F3z_S_S_b = I_ERI_G2y2z_D2z_S_S_b+ABZ*I_ERI_F2yz_D2z_S_S_b;
  Double I_ERI_Fy2z_F3z_S_S_b = I_ERI_Gy3z_D2z_S_S_b+ABZ*I_ERI_Fy2z_D2z_S_S_b;
  Double I_ERI_F3z_F3z_S_S_b = I_ERI_G4z_D2z_S_S_b+ABZ*I_ERI_F3z_D2z_S_S_b;

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
   * shell quartet name: SQ_ERI_G_P_P_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 18 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_P_S_c
   * RHS shell quartet name: SQ_ERI_G_S_P_S_c
   ************************************************************/
  Double I_ERI_G4x_Px_Px_S_c = I_ERI_H5x_S_Px_S_c+ABX*I_ERI_G4x_S_Px_S_c;
  Double I_ERI_G3xy_Px_Px_S_c = I_ERI_H4xy_S_Px_S_c+ABX*I_ERI_G3xy_S_Px_S_c;
  Double I_ERI_G3xz_Px_Px_S_c = I_ERI_H4xz_S_Px_S_c+ABX*I_ERI_G3xz_S_Px_S_c;
  Double I_ERI_G2x2y_Px_Px_S_c = I_ERI_H3x2y_S_Px_S_c+ABX*I_ERI_G2x2y_S_Px_S_c;
  Double I_ERI_G2xyz_Px_Px_S_c = I_ERI_H3xyz_S_Px_S_c+ABX*I_ERI_G2xyz_S_Px_S_c;
  Double I_ERI_G2x2z_Px_Px_S_c = I_ERI_H3x2z_S_Px_S_c+ABX*I_ERI_G2x2z_S_Px_S_c;
  Double I_ERI_Gx3y_Px_Px_S_c = I_ERI_H2x3y_S_Px_S_c+ABX*I_ERI_Gx3y_S_Px_S_c;
  Double I_ERI_Gx2yz_Px_Px_S_c = I_ERI_H2x2yz_S_Px_S_c+ABX*I_ERI_Gx2yz_S_Px_S_c;
  Double I_ERI_Gxy2z_Px_Px_S_c = I_ERI_H2xy2z_S_Px_S_c+ABX*I_ERI_Gxy2z_S_Px_S_c;
  Double I_ERI_Gx3z_Px_Px_S_c = I_ERI_H2x3z_S_Px_S_c+ABX*I_ERI_Gx3z_S_Px_S_c;
  Double I_ERI_G4y_Px_Px_S_c = I_ERI_Hx4y_S_Px_S_c+ABX*I_ERI_G4y_S_Px_S_c;
  Double I_ERI_G3yz_Px_Px_S_c = I_ERI_Hx3yz_S_Px_S_c+ABX*I_ERI_G3yz_S_Px_S_c;
  Double I_ERI_G2y2z_Px_Px_S_c = I_ERI_Hx2y2z_S_Px_S_c+ABX*I_ERI_G2y2z_S_Px_S_c;
  Double I_ERI_Gy3z_Px_Px_S_c = I_ERI_Hxy3z_S_Px_S_c+ABX*I_ERI_Gy3z_S_Px_S_c;
  Double I_ERI_G4z_Px_Px_S_c = I_ERI_Hx4z_S_Px_S_c+ABX*I_ERI_G4z_S_Px_S_c;
  Double I_ERI_G3xy_Py_Px_S_c = I_ERI_H3x2y_S_Px_S_c+ABY*I_ERI_G3xy_S_Px_S_c;
  Double I_ERI_G3xz_Py_Px_S_c = I_ERI_H3xyz_S_Px_S_c+ABY*I_ERI_G3xz_S_Px_S_c;
  Double I_ERI_G2x2y_Py_Px_S_c = I_ERI_H2x3y_S_Px_S_c+ABY*I_ERI_G2x2y_S_Px_S_c;
  Double I_ERI_G2xyz_Py_Px_S_c = I_ERI_H2x2yz_S_Px_S_c+ABY*I_ERI_G2xyz_S_Px_S_c;
  Double I_ERI_G2x2z_Py_Px_S_c = I_ERI_H2xy2z_S_Px_S_c+ABY*I_ERI_G2x2z_S_Px_S_c;
  Double I_ERI_Gx3y_Py_Px_S_c = I_ERI_Hx4y_S_Px_S_c+ABY*I_ERI_Gx3y_S_Px_S_c;
  Double I_ERI_Gx2yz_Py_Px_S_c = I_ERI_Hx3yz_S_Px_S_c+ABY*I_ERI_Gx2yz_S_Px_S_c;
  Double I_ERI_Gxy2z_Py_Px_S_c = I_ERI_Hx2y2z_S_Px_S_c+ABY*I_ERI_Gxy2z_S_Px_S_c;
  Double I_ERI_Gx3z_Py_Px_S_c = I_ERI_Hxy3z_S_Px_S_c+ABY*I_ERI_Gx3z_S_Px_S_c;
  Double I_ERI_G4y_Py_Px_S_c = I_ERI_H5y_S_Px_S_c+ABY*I_ERI_G4y_S_Px_S_c;
  Double I_ERI_G3yz_Py_Px_S_c = I_ERI_H4yz_S_Px_S_c+ABY*I_ERI_G3yz_S_Px_S_c;
  Double I_ERI_G2y2z_Py_Px_S_c = I_ERI_H3y2z_S_Px_S_c+ABY*I_ERI_G2y2z_S_Px_S_c;
  Double I_ERI_Gy3z_Py_Px_S_c = I_ERI_H2y3z_S_Px_S_c+ABY*I_ERI_Gy3z_S_Px_S_c;
  Double I_ERI_G4z_Py_Px_S_c = I_ERI_Hy4z_S_Px_S_c+ABY*I_ERI_G4z_S_Px_S_c;
  Double I_ERI_G3xz_Pz_Px_S_c = I_ERI_H3x2z_S_Px_S_c+ABZ*I_ERI_G3xz_S_Px_S_c;
  Double I_ERI_G2xyz_Pz_Px_S_c = I_ERI_H2xy2z_S_Px_S_c+ABZ*I_ERI_G2xyz_S_Px_S_c;
  Double I_ERI_G2x2z_Pz_Px_S_c = I_ERI_H2x3z_S_Px_S_c+ABZ*I_ERI_G2x2z_S_Px_S_c;
  Double I_ERI_Gx2yz_Pz_Px_S_c = I_ERI_Hx2y2z_S_Px_S_c+ABZ*I_ERI_Gx2yz_S_Px_S_c;
  Double I_ERI_Gxy2z_Pz_Px_S_c = I_ERI_Hxy3z_S_Px_S_c+ABZ*I_ERI_Gxy2z_S_Px_S_c;
  Double I_ERI_Gx3z_Pz_Px_S_c = I_ERI_Hx4z_S_Px_S_c+ABZ*I_ERI_Gx3z_S_Px_S_c;
  Double I_ERI_G3yz_Pz_Px_S_c = I_ERI_H3y2z_S_Px_S_c+ABZ*I_ERI_G3yz_S_Px_S_c;
  Double I_ERI_G2y2z_Pz_Px_S_c = I_ERI_H2y3z_S_Px_S_c+ABZ*I_ERI_G2y2z_S_Px_S_c;
  Double I_ERI_Gy3z_Pz_Px_S_c = I_ERI_Hy4z_S_Px_S_c+ABZ*I_ERI_Gy3z_S_Px_S_c;
  Double I_ERI_G4z_Pz_Px_S_c = I_ERI_H5z_S_Px_S_c+ABZ*I_ERI_G4z_S_Px_S_c;
  Double I_ERI_G4x_Px_Py_S_c = I_ERI_H5x_S_Py_S_c+ABX*I_ERI_G4x_S_Py_S_c;
  Double I_ERI_G3xy_Px_Py_S_c = I_ERI_H4xy_S_Py_S_c+ABX*I_ERI_G3xy_S_Py_S_c;
  Double I_ERI_G3xz_Px_Py_S_c = I_ERI_H4xz_S_Py_S_c+ABX*I_ERI_G3xz_S_Py_S_c;
  Double I_ERI_G2x2y_Px_Py_S_c = I_ERI_H3x2y_S_Py_S_c+ABX*I_ERI_G2x2y_S_Py_S_c;
  Double I_ERI_G2xyz_Px_Py_S_c = I_ERI_H3xyz_S_Py_S_c+ABX*I_ERI_G2xyz_S_Py_S_c;
  Double I_ERI_G2x2z_Px_Py_S_c = I_ERI_H3x2z_S_Py_S_c+ABX*I_ERI_G2x2z_S_Py_S_c;
  Double I_ERI_Gx3y_Px_Py_S_c = I_ERI_H2x3y_S_Py_S_c+ABX*I_ERI_Gx3y_S_Py_S_c;
  Double I_ERI_Gx2yz_Px_Py_S_c = I_ERI_H2x2yz_S_Py_S_c+ABX*I_ERI_Gx2yz_S_Py_S_c;
  Double I_ERI_Gxy2z_Px_Py_S_c = I_ERI_H2xy2z_S_Py_S_c+ABX*I_ERI_Gxy2z_S_Py_S_c;
  Double I_ERI_Gx3z_Px_Py_S_c = I_ERI_H2x3z_S_Py_S_c+ABX*I_ERI_Gx3z_S_Py_S_c;
  Double I_ERI_G4y_Px_Py_S_c = I_ERI_Hx4y_S_Py_S_c+ABX*I_ERI_G4y_S_Py_S_c;
  Double I_ERI_G3yz_Px_Py_S_c = I_ERI_Hx3yz_S_Py_S_c+ABX*I_ERI_G3yz_S_Py_S_c;
  Double I_ERI_G2y2z_Px_Py_S_c = I_ERI_Hx2y2z_S_Py_S_c+ABX*I_ERI_G2y2z_S_Py_S_c;
  Double I_ERI_Gy3z_Px_Py_S_c = I_ERI_Hxy3z_S_Py_S_c+ABX*I_ERI_Gy3z_S_Py_S_c;
  Double I_ERI_G4z_Px_Py_S_c = I_ERI_Hx4z_S_Py_S_c+ABX*I_ERI_G4z_S_Py_S_c;
  Double I_ERI_G3xy_Py_Py_S_c = I_ERI_H3x2y_S_Py_S_c+ABY*I_ERI_G3xy_S_Py_S_c;
  Double I_ERI_G3xz_Py_Py_S_c = I_ERI_H3xyz_S_Py_S_c+ABY*I_ERI_G3xz_S_Py_S_c;
  Double I_ERI_G2x2y_Py_Py_S_c = I_ERI_H2x3y_S_Py_S_c+ABY*I_ERI_G2x2y_S_Py_S_c;
  Double I_ERI_G2xyz_Py_Py_S_c = I_ERI_H2x2yz_S_Py_S_c+ABY*I_ERI_G2xyz_S_Py_S_c;
  Double I_ERI_G2x2z_Py_Py_S_c = I_ERI_H2xy2z_S_Py_S_c+ABY*I_ERI_G2x2z_S_Py_S_c;
  Double I_ERI_Gx3y_Py_Py_S_c = I_ERI_Hx4y_S_Py_S_c+ABY*I_ERI_Gx3y_S_Py_S_c;
  Double I_ERI_Gx2yz_Py_Py_S_c = I_ERI_Hx3yz_S_Py_S_c+ABY*I_ERI_Gx2yz_S_Py_S_c;
  Double I_ERI_Gxy2z_Py_Py_S_c = I_ERI_Hx2y2z_S_Py_S_c+ABY*I_ERI_Gxy2z_S_Py_S_c;
  Double I_ERI_Gx3z_Py_Py_S_c = I_ERI_Hxy3z_S_Py_S_c+ABY*I_ERI_Gx3z_S_Py_S_c;
  Double I_ERI_G4y_Py_Py_S_c = I_ERI_H5y_S_Py_S_c+ABY*I_ERI_G4y_S_Py_S_c;
  Double I_ERI_G3yz_Py_Py_S_c = I_ERI_H4yz_S_Py_S_c+ABY*I_ERI_G3yz_S_Py_S_c;
  Double I_ERI_G2y2z_Py_Py_S_c = I_ERI_H3y2z_S_Py_S_c+ABY*I_ERI_G2y2z_S_Py_S_c;
  Double I_ERI_Gy3z_Py_Py_S_c = I_ERI_H2y3z_S_Py_S_c+ABY*I_ERI_Gy3z_S_Py_S_c;
  Double I_ERI_G4z_Py_Py_S_c = I_ERI_Hy4z_S_Py_S_c+ABY*I_ERI_G4z_S_Py_S_c;
  Double I_ERI_G3xz_Pz_Py_S_c = I_ERI_H3x2z_S_Py_S_c+ABZ*I_ERI_G3xz_S_Py_S_c;
  Double I_ERI_G2xyz_Pz_Py_S_c = I_ERI_H2xy2z_S_Py_S_c+ABZ*I_ERI_G2xyz_S_Py_S_c;
  Double I_ERI_G2x2z_Pz_Py_S_c = I_ERI_H2x3z_S_Py_S_c+ABZ*I_ERI_G2x2z_S_Py_S_c;
  Double I_ERI_Gx2yz_Pz_Py_S_c = I_ERI_Hx2y2z_S_Py_S_c+ABZ*I_ERI_Gx2yz_S_Py_S_c;
  Double I_ERI_Gxy2z_Pz_Py_S_c = I_ERI_Hxy3z_S_Py_S_c+ABZ*I_ERI_Gxy2z_S_Py_S_c;
  Double I_ERI_Gx3z_Pz_Py_S_c = I_ERI_Hx4z_S_Py_S_c+ABZ*I_ERI_Gx3z_S_Py_S_c;
  Double I_ERI_G3yz_Pz_Py_S_c = I_ERI_H3y2z_S_Py_S_c+ABZ*I_ERI_G3yz_S_Py_S_c;
  Double I_ERI_G2y2z_Pz_Py_S_c = I_ERI_H2y3z_S_Py_S_c+ABZ*I_ERI_G2y2z_S_Py_S_c;
  Double I_ERI_Gy3z_Pz_Py_S_c = I_ERI_Hy4z_S_Py_S_c+ABZ*I_ERI_Gy3z_S_Py_S_c;
  Double I_ERI_G4z_Pz_Py_S_c = I_ERI_H5z_S_Py_S_c+ABZ*I_ERI_G4z_S_Py_S_c;
  Double I_ERI_G4x_Px_Pz_S_c = I_ERI_H5x_S_Pz_S_c+ABX*I_ERI_G4x_S_Pz_S_c;
  Double I_ERI_G3xy_Px_Pz_S_c = I_ERI_H4xy_S_Pz_S_c+ABX*I_ERI_G3xy_S_Pz_S_c;
  Double I_ERI_G3xz_Px_Pz_S_c = I_ERI_H4xz_S_Pz_S_c+ABX*I_ERI_G3xz_S_Pz_S_c;
  Double I_ERI_G2x2y_Px_Pz_S_c = I_ERI_H3x2y_S_Pz_S_c+ABX*I_ERI_G2x2y_S_Pz_S_c;
  Double I_ERI_G2xyz_Px_Pz_S_c = I_ERI_H3xyz_S_Pz_S_c+ABX*I_ERI_G2xyz_S_Pz_S_c;
  Double I_ERI_G2x2z_Px_Pz_S_c = I_ERI_H3x2z_S_Pz_S_c+ABX*I_ERI_G2x2z_S_Pz_S_c;
  Double I_ERI_Gx3y_Px_Pz_S_c = I_ERI_H2x3y_S_Pz_S_c+ABX*I_ERI_Gx3y_S_Pz_S_c;
  Double I_ERI_Gx2yz_Px_Pz_S_c = I_ERI_H2x2yz_S_Pz_S_c+ABX*I_ERI_Gx2yz_S_Pz_S_c;
  Double I_ERI_Gxy2z_Px_Pz_S_c = I_ERI_H2xy2z_S_Pz_S_c+ABX*I_ERI_Gxy2z_S_Pz_S_c;
  Double I_ERI_Gx3z_Px_Pz_S_c = I_ERI_H2x3z_S_Pz_S_c+ABX*I_ERI_Gx3z_S_Pz_S_c;
  Double I_ERI_G4y_Px_Pz_S_c = I_ERI_Hx4y_S_Pz_S_c+ABX*I_ERI_G4y_S_Pz_S_c;
  Double I_ERI_G3yz_Px_Pz_S_c = I_ERI_Hx3yz_S_Pz_S_c+ABX*I_ERI_G3yz_S_Pz_S_c;
  Double I_ERI_G2y2z_Px_Pz_S_c = I_ERI_Hx2y2z_S_Pz_S_c+ABX*I_ERI_G2y2z_S_Pz_S_c;
  Double I_ERI_Gy3z_Px_Pz_S_c = I_ERI_Hxy3z_S_Pz_S_c+ABX*I_ERI_Gy3z_S_Pz_S_c;
  Double I_ERI_G4z_Px_Pz_S_c = I_ERI_Hx4z_S_Pz_S_c+ABX*I_ERI_G4z_S_Pz_S_c;
  Double I_ERI_G3xy_Py_Pz_S_c = I_ERI_H3x2y_S_Pz_S_c+ABY*I_ERI_G3xy_S_Pz_S_c;
  Double I_ERI_G3xz_Py_Pz_S_c = I_ERI_H3xyz_S_Pz_S_c+ABY*I_ERI_G3xz_S_Pz_S_c;
  Double I_ERI_G2x2y_Py_Pz_S_c = I_ERI_H2x3y_S_Pz_S_c+ABY*I_ERI_G2x2y_S_Pz_S_c;
  Double I_ERI_G2xyz_Py_Pz_S_c = I_ERI_H2x2yz_S_Pz_S_c+ABY*I_ERI_G2xyz_S_Pz_S_c;
  Double I_ERI_G2x2z_Py_Pz_S_c = I_ERI_H2xy2z_S_Pz_S_c+ABY*I_ERI_G2x2z_S_Pz_S_c;
  Double I_ERI_Gx3y_Py_Pz_S_c = I_ERI_Hx4y_S_Pz_S_c+ABY*I_ERI_Gx3y_S_Pz_S_c;
  Double I_ERI_Gx2yz_Py_Pz_S_c = I_ERI_Hx3yz_S_Pz_S_c+ABY*I_ERI_Gx2yz_S_Pz_S_c;
  Double I_ERI_Gxy2z_Py_Pz_S_c = I_ERI_Hx2y2z_S_Pz_S_c+ABY*I_ERI_Gxy2z_S_Pz_S_c;
  Double I_ERI_Gx3z_Py_Pz_S_c = I_ERI_Hxy3z_S_Pz_S_c+ABY*I_ERI_Gx3z_S_Pz_S_c;
  Double I_ERI_G4y_Py_Pz_S_c = I_ERI_H5y_S_Pz_S_c+ABY*I_ERI_G4y_S_Pz_S_c;
  Double I_ERI_G3yz_Py_Pz_S_c = I_ERI_H4yz_S_Pz_S_c+ABY*I_ERI_G3yz_S_Pz_S_c;
  Double I_ERI_G2y2z_Py_Pz_S_c = I_ERI_H3y2z_S_Pz_S_c+ABY*I_ERI_G2y2z_S_Pz_S_c;
  Double I_ERI_Gy3z_Py_Pz_S_c = I_ERI_H2y3z_S_Pz_S_c+ABY*I_ERI_Gy3z_S_Pz_S_c;
  Double I_ERI_G4z_Py_Pz_S_c = I_ERI_Hy4z_S_Pz_S_c+ABY*I_ERI_G4z_S_Pz_S_c;
  Double I_ERI_G3xz_Pz_Pz_S_c = I_ERI_H3x2z_S_Pz_S_c+ABZ*I_ERI_G3xz_S_Pz_S_c;
  Double I_ERI_G2xyz_Pz_Pz_S_c = I_ERI_H2xy2z_S_Pz_S_c+ABZ*I_ERI_G2xyz_S_Pz_S_c;
  Double I_ERI_G2x2z_Pz_Pz_S_c = I_ERI_H2x3z_S_Pz_S_c+ABZ*I_ERI_G2x2z_S_Pz_S_c;
  Double I_ERI_Gx2yz_Pz_Pz_S_c = I_ERI_Hx2y2z_S_Pz_S_c+ABZ*I_ERI_Gx2yz_S_Pz_S_c;
  Double I_ERI_Gxy2z_Pz_Pz_S_c = I_ERI_Hxy3z_S_Pz_S_c+ABZ*I_ERI_Gxy2z_S_Pz_S_c;
  Double I_ERI_Gx3z_Pz_Pz_S_c = I_ERI_Hx4z_S_Pz_S_c+ABZ*I_ERI_Gx3z_S_Pz_S_c;
  Double I_ERI_G3yz_Pz_Pz_S_c = I_ERI_H3y2z_S_Pz_S_c+ABZ*I_ERI_G3yz_S_Pz_S_c;
  Double I_ERI_G2y2z_Pz_Pz_S_c = I_ERI_H2y3z_S_Pz_S_c+ABZ*I_ERI_G2y2z_S_Pz_S_c;
  Double I_ERI_Gy3z_Pz_Pz_S_c = I_ERI_Hy4z_S_Pz_S_c+ABZ*I_ERI_Gy3z_S_Pz_S_c;
  Double I_ERI_G4z_Pz_Pz_S_c = I_ERI_H5z_S_Pz_S_c+ABZ*I_ERI_G4z_S_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_P_S_c
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_S_c
   * RHS shell quartet name: SQ_ERI_F_P_P_S_c
   ************************************************************/
  Double I_ERI_F3x_D2x_Px_S_c = I_ERI_G4x_Px_Px_S_c+ABX*I_ERI_F3x_Px_Px_S_c;
  Double I_ERI_F2xy_D2x_Px_S_c = I_ERI_G3xy_Px_Px_S_c+ABX*I_ERI_F2xy_Px_Px_S_c;
  Double I_ERI_F2xz_D2x_Px_S_c = I_ERI_G3xz_Px_Px_S_c+ABX*I_ERI_F2xz_Px_Px_S_c;
  Double I_ERI_Fx2y_D2x_Px_S_c = I_ERI_G2x2y_Px_Px_S_c+ABX*I_ERI_Fx2y_Px_Px_S_c;
  Double I_ERI_Fxyz_D2x_Px_S_c = I_ERI_G2xyz_Px_Px_S_c+ABX*I_ERI_Fxyz_Px_Px_S_c;
  Double I_ERI_Fx2z_D2x_Px_S_c = I_ERI_G2x2z_Px_Px_S_c+ABX*I_ERI_Fx2z_Px_Px_S_c;
  Double I_ERI_F3y_D2x_Px_S_c = I_ERI_Gx3y_Px_Px_S_c+ABX*I_ERI_F3y_Px_Px_S_c;
  Double I_ERI_F2yz_D2x_Px_S_c = I_ERI_Gx2yz_Px_Px_S_c+ABX*I_ERI_F2yz_Px_Px_S_c;
  Double I_ERI_Fy2z_D2x_Px_S_c = I_ERI_Gxy2z_Px_Px_S_c+ABX*I_ERI_Fy2z_Px_Px_S_c;
  Double I_ERI_F3z_D2x_Px_S_c = I_ERI_Gx3z_Px_Px_S_c+ABX*I_ERI_F3z_Px_Px_S_c;
  Double I_ERI_F3x_Dxy_Px_S_c = I_ERI_G3xy_Px_Px_S_c+ABY*I_ERI_F3x_Px_Px_S_c;
  Double I_ERI_F2xy_Dxy_Px_S_c = I_ERI_G2x2y_Px_Px_S_c+ABY*I_ERI_F2xy_Px_Px_S_c;
  Double I_ERI_F2xz_Dxy_Px_S_c = I_ERI_G2xyz_Px_Px_S_c+ABY*I_ERI_F2xz_Px_Px_S_c;
  Double I_ERI_Fx2y_Dxy_Px_S_c = I_ERI_Gx3y_Px_Px_S_c+ABY*I_ERI_Fx2y_Px_Px_S_c;
  Double I_ERI_Fxyz_Dxy_Px_S_c = I_ERI_Gx2yz_Px_Px_S_c+ABY*I_ERI_Fxyz_Px_Px_S_c;
  Double I_ERI_Fx2z_Dxy_Px_S_c = I_ERI_Gxy2z_Px_Px_S_c+ABY*I_ERI_Fx2z_Px_Px_S_c;
  Double I_ERI_F3y_Dxy_Px_S_c = I_ERI_G4y_Px_Px_S_c+ABY*I_ERI_F3y_Px_Px_S_c;
  Double I_ERI_F2yz_Dxy_Px_S_c = I_ERI_G3yz_Px_Px_S_c+ABY*I_ERI_F2yz_Px_Px_S_c;
  Double I_ERI_Fy2z_Dxy_Px_S_c = I_ERI_G2y2z_Px_Px_S_c+ABY*I_ERI_Fy2z_Px_Px_S_c;
  Double I_ERI_F3z_Dxy_Px_S_c = I_ERI_Gy3z_Px_Px_S_c+ABY*I_ERI_F3z_Px_Px_S_c;
  Double I_ERI_F3x_Dxz_Px_S_c = I_ERI_G3xz_Px_Px_S_c+ABZ*I_ERI_F3x_Px_Px_S_c;
  Double I_ERI_F2xy_Dxz_Px_S_c = I_ERI_G2xyz_Px_Px_S_c+ABZ*I_ERI_F2xy_Px_Px_S_c;
  Double I_ERI_F2xz_Dxz_Px_S_c = I_ERI_G2x2z_Px_Px_S_c+ABZ*I_ERI_F2xz_Px_Px_S_c;
  Double I_ERI_Fx2y_Dxz_Px_S_c = I_ERI_Gx2yz_Px_Px_S_c+ABZ*I_ERI_Fx2y_Px_Px_S_c;
  Double I_ERI_Fxyz_Dxz_Px_S_c = I_ERI_Gxy2z_Px_Px_S_c+ABZ*I_ERI_Fxyz_Px_Px_S_c;
  Double I_ERI_Fx2z_Dxz_Px_S_c = I_ERI_Gx3z_Px_Px_S_c+ABZ*I_ERI_Fx2z_Px_Px_S_c;
  Double I_ERI_F3y_Dxz_Px_S_c = I_ERI_G3yz_Px_Px_S_c+ABZ*I_ERI_F3y_Px_Px_S_c;
  Double I_ERI_F2yz_Dxz_Px_S_c = I_ERI_G2y2z_Px_Px_S_c+ABZ*I_ERI_F2yz_Px_Px_S_c;
  Double I_ERI_Fy2z_Dxz_Px_S_c = I_ERI_Gy3z_Px_Px_S_c+ABZ*I_ERI_Fy2z_Px_Px_S_c;
  Double I_ERI_F3z_Dxz_Px_S_c = I_ERI_G4z_Px_Px_S_c+ABZ*I_ERI_F3z_Px_Px_S_c;
  Double I_ERI_F3x_D2y_Px_S_c = I_ERI_G3xy_Py_Px_S_c+ABY*I_ERI_F3x_Py_Px_S_c;
  Double I_ERI_F2xy_D2y_Px_S_c = I_ERI_G2x2y_Py_Px_S_c+ABY*I_ERI_F2xy_Py_Px_S_c;
  Double I_ERI_F2xz_D2y_Px_S_c = I_ERI_G2xyz_Py_Px_S_c+ABY*I_ERI_F2xz_Py_Px_S_c;
  Double I_ERI_Fx2y_D2y_Px_S_c = I_ERI_Gx3y_Py_Px_S_c+ABY*I_ERI_Fx2y_Py_Px_S_c;
  Double I_ERI_Fxyz_D2y_Px_S_c = I_ERI_Gx2yz_Py_Px_S_c+ABY*I_ERI_Fxyz_Py_Px_S_c;
  Double I_ERI_Fx2z_D2y_Px_S_c = I_ERI_Gxy2z_Py_Px_S_c+ABY*I_ERI_Fx2z_Py_Px_S_c;
  Double I_ERI_F3y_D2y_Px_S_c = I_ERI_G4y_Py_Px_S_c+ABY*I_ERI_F3y_Py_Px_S_c;
  Double I_ERI_F2yz_D2y_Px_S_c = I_ERI_G3yz_Py_Px_S_c+ABY*I_ERI_F2yz_Py_Px_S_c;
  Double I_ERI_Fy2z_D2y_Px_S_c = I_ERI_G2y2z_Py_Px_S_c+ABY*I_ERI_Fy2z_Py_Px_S_c;
  Double I_ERI_F3z_D2y_Px_S_c = I_ERI_Gy3z_Py_Px_S_c+ABY*I_ERI_F3z_Py_Px_S_c;
  Double I_ERI_F3x_Dyz_Px_S_c = I_ERI_G3xz_Py_Px_S_c+ABZ*I_ERI_F3x_Py_Px_S_c;
  Double I_ERI_F2xy_Dyz_Px_S_c = I_ERI_G2xyz_Py_Px_S_c+ABZ*I_ERI_F2xy_Py_Px_S_c;
  Double I_ERI_F2xz_Dyz_Px_S_c = I_ERI_G2x2z_Py_Px_S_c+ABZ*I_ERI_F2xz_Py_Px_S_c;
  Double I_ERI_Fx2y_Dyz_Px_S_c = I_ERI_Gx2yz_Py_Px_S_c+ABZ*I_ERI_Fx2y_Py_Px_S_c;
  Double I_ERI_Fxyz_Dyz_Px_S_c = I_ERI_Gxy2z_Py_Px_S_c+ABZ*I_ERI_Fxyz_Py_Px_S_c;
  Double I_ERI_Fx2z_Dyz_Px_S_c = I_ERI_Gx3z_Py_Px_S_c+ABZ*I_ERI_Fx2z_Py_Px_S_c;
  Double I_ERI_F3y_Dyz_Px_S_c = I_ERI_G3yz_Py_Px_S_c+ABZ*I_ERI_F3y_Py_Px_S_c;
  Double I_ERI_F2yz_Dyz_Px_S_c = I_ERI_G2y2z_Py_Px_S_c+ABZ*I_ERI_F2yz_Py_Px_S_c;
  Double I_ERI_Fy2z_Dyz_Px_S_c = I_ERI_Gy3z_Py_Px_S_c+ABZ*I_ERI_Fy2z_Py_Px_S_c;
  Double I_ERI_F3z_Dyz_Px_S_c = I_ERI_G4z_Py_Px_S_c+ABZ*I_ERI_F3z_Py_Px_S_c;
  Double I_ERI_F3x_D2z_Px_S_c = I_ERI_G3xz_Pz_Px_S_c+ABZ*I_ERI_F3x_Pz_Px_S_c;
  Double I_ERI_F2xy_D2z_Px_S_c = I_ERI_G2xyz_Pz_Px_S_c+ABZ*I_ERI_F2xy_Pz_Px_S_c;
  Double I_ERI_F2xz_D2z_Px_S_c = I_ERI_G2x2z_Pz_Px_S_c+ABZ*I_ERI_F2xz_Pz_Px_S_c;
  Double I_ERI_Fx2y_D2z_Px_S_c = I_ERI_Gx2yz_Pz_Px_S_c+ABZ*I_ERI_Fx2y_Pz_Px_S_c;
  Double I_ERI_Fxyz_D2z_Px_S_c = I_ERI_Gxy2z_Pz_Px_S_c+ABZ*I_ERI_Fxyz_Pz_Px_S_c;
  Double I_ERI_Fx2z_D2z_Px_S_c = I_ERI_Gx3z_Pz_Px_S_c+ABZ*I_ERI_Fx2z_Pz_Px_S_c;
  Double I_ERI_F3y_D2z_Px_S_c = I_ERI_G3yz_Pz_Px_S_c+ABZ*I_ERI_F3y_Pz_Px_S_c;
  Double I_ERI_F2yz_D2z_Px_S_c = I_ERI_G2y2z_Pz_Px_S_c+ABZ*I_ERI_F2yz_Pz_Px_S_c;
  Double I_ERI_Fy2z_D2z_Px_S_c = I_ERI_Gy3z_Pz_Px_S_c+ABZ*I_ERI_Fy2z_Pz_Px_S_c;
  Double I_ERI_F3z_D2z_Px_S_c = I_ERI_G4z_Pz_Px_S_c+ABZ*I_ERI_F3z_Pz_Px_S_c;
  Double I_ERI_F3x_D2x_Py_S_c = I_ERI_G4x_Px_Py_S_c+ABX*I_ERI_F3x_Px_Py_S_c;
  Double I_ERI_F2xy_D2x_Py_S_c = I_ERI_G3xy_Px_Py_S_c+ABX*I_ERI_F2xy_Px_Py_S_c;
  Double I_ERI_F2xz_D2x_Py_S_c = I_ERI_G3xz_Px_Py_S_c+ABX*I_ERI_F2xz_Px_Py_S_c;
  Double I_ERI_Fx2y_D2x_Py_S_c = I_ERI_G2x2y_Px_Py_S_c+ABX*I_ERI_Fx2y_Px_Py_S_c;
  Double I_ERI_Fxyz_D2x_Py_S_c = I_ERI_G2xyz_Px_Py_S_c+ABX*I_ERI_Fxyz_Px_Py_S_c;
  Double I_ERI_Fx2z_D2x_Py_S_c = I_ERI_G2x2z_Px_Py_S_c+ABX*I_ERI_Fx2z_Px_Py_S_c;
  Double I_ERI_F3y_D2x_Py_S_c = I_ERI_Gx3y_Px_Py_S_c+ABX*I_ERI_F3y_Px_Py_S_c;
  Double I_ERI_F2yz_D2x_Py_S_c = I_ERI_Gx2yz_Px_Py_S_c+ABX*I_ERI_F2yz_Px_Py_S_c;
  Double I_ERI_Fy2z_D2x_Py_S_c = I_ERI_Gxy2z_Px_Py_S_c+ABX*I_ERI_Fy2z_Px_Py_S_c;
  Double I_ERI_F3z_D2x_Py_S_c = I_ERI_Gx3z_Px_Py_S_c+ABX*I_ERI_F3z_Px_Py_S_c;
  Double I_ERI_F3x_Dxy_Py_S_c = I_ERI_G3xy_Px_Py_S_c+ABY*I_ERI_F3x_Px_Py_S_c;
  Double I_ERI_F2xy_Dxy_Py_S_c = I_ERI_G2x2y_Px_Py_S_c+ABY*I_ERI_F2xy_Px_Py_S_c;
  Double I_ERI_F2xz_Dxy_Py_S_c = I_ERI_G2xyz_Px_Py_S_c+ABY*I_ERI_F2xz_Px_Py_S_c;
  Double I_ERI_Fx2y_Dxy_Py_S_c = I_ERI_Gx3y_Px_Py_S_c+ABY*I_ERI_Fx2y_Px_Py_S_c;
  Double I_ERI_Fxyz_Dxy_Py_S_c = I_ERI_Gx2yz_Px_Py_S_c+ABY*I_ERI_Fxyz_Px_Py_S_c;
  Double I_ERI_Fx2z_Dxy_Py_S_c = I_ERI_Gxy2z_Px_Py_S_c+ABY*I_ERI_Fx2z_Px_Py_S_c;
  Double I_ERI_F3y_Dxy_Py_S_c = I_ERI_G4y_Px_Py_S_c+ABY*I_ERI_F3y_Px_Py_S_c;
  Double I_ERI_F2yz_Dxy_Py_S_c = I_ERI_G3yz_Px_Py_S_c+ABY*I_ERI_F2yz_Px_Py_S_c;
  Double I_ERI_Fy2z_Dxy_Py_S_c = I_ERI_G2y2z_Px_Py_S_c+ABY*I_ERI_Fy2z_Px_Py_S_c;
  Double I_ERI_F3z_Dxy_Py_S_c = I_ERI_Gy3z_Px_Py_S_c+ABY*I_ERI_F3z_Px_Py_S_c;
  Double I_ERI_F3x_Dxz_Py_S_c = I_ERI_G3xz_Px_Py_S_c+ABZ*I_ERI_F3x_Px_Py_S_c;
  Double I_ERI_F2xy_Dxz_Py_S_c = I_ERI_G2xyz_Px_Py_S_c+ABZ*I_ERI_F2xy_Px_Py_S_c;
  Double I_ERI_F2xz_Dxz_Py_S_c = I_ERI_G2x2z_Px_Py_S_c+ABZ*I_ERI_F2xz_Px_Py_S_c;
  Double I_ERI_Fx2y_Dxz_Py_S_c = I_ERI_Gx2yz_Px_Py_S_c+ABZ*I_ERI_Fx2y_Px_Py_S_c;
  Double I_ERI_Fxyz_Dxz_Py_S_c = I_ERI_Gxy2z_Px_Py_S_c+ABZ*I_ERI_Fxyz_Px_Py_S_c;
  Double I_ERI_Fx2z_Dxz_Py_S_c = I_ERI_Gx3z_Px_Py_S_c+ABZ*I_ERI_Fx2z_Px_Py_S_c;
  Double I_ERI_F3y_Dxz_Py_S_c = I_ERI_G3yz_Px_Py_S_c+ABZ*I_ERI_F3y_Px_Py_S_c;
  Double I_ERI_F2yz_Dxz_Py_S_c = I_ERI_G2y2z_Px_Py_S_c+ABZ*I_ERI_F2yz_Px_Py_S_c;
  Double I_ERI_Fy2z_Dxz_Py_S_c = I_ERI_Gy3z_Px_Py_S_c+ABZ*I_ERI_Fy2z_Px_Py_S_c;
  Double I_ERI_F3z_Dxz_Py_S_c = I_ERI_G4z_Px_Py_S_c+ABZ*I_ERI_F3z_Px_Py_S_c;
  Double I_ERI_F3x_D2y_Py_S_c = I_ERI_G3xy_Py_Py_S_c+ABY*I_ERI_F3x_Py_Py_S_c;
  Double I_ERI_F2xy_D2y_Py_S_c = I_ERI_G2x2y_Py_Py_S_c+ABY*I_ERI_F2xy_Py_Py_S_c;
  Double I_ERI_F2xz_D2y_Py_S_c = I_ERI_G2xyz_Py_Py_S_c+ABY*I_ERI_F2xz_Py_Py_S_c;
  Double I_ERI_Fx2y_D2y_Py_S_c = I_ERI_Gx3y_Py_Py_S_c+ABY*I_ERI_Fx2y_Py_Py_S_c;
  Double I_ERI_Fxyz_D2y_Py_S_c = I_ERI_Gx2yz_Py_Py_S_c+ABY*I_ERI_Fxyz_Py_Py_S_c;
  Double I_ERI_Fx2z_D2y_Py_S_c = I_ERI_Gxy2z_Py_Py_S_c+ABY*I_ERI_Fx2z_Py_Py_S_c;
  Double I_ERI_F3y_D2y_Py_S_c = I_ERI_G4y_Py_Py_S_c+ABY*I_ERI_F3y_Py_Py_S_c;
  Double I_ERI_F2yz_D2y_Py_S_c = I_ERI_G3yz_Py_Py_S_c+ABY*I_ERI_F2yz_Py_Py_S_c;
  Double I_ERI_Fy2z_D2y_Py_S_c = I_ERI_G2y2z_Py_Py_S_c+ABY*I_ERI_Fy2z_Py_Py_S_c;
  Double I_ERI_F3z_D2y_Py_S_c = I_ERI_Gy3z_Py_Py_S_c+ABY*I_ERI_F3z_Py_Py_S_c;
  Double I_ERI_F3x_Dyz_Py_S_c = I_ERI_G3xz_Py_Py_S_c+ABZ*I_ERI_F3x_Py_Py_S_c;
  Double I_ERI_F2xy_Dyz_Py_S_c = I_ERI_G2xyz_Py_Py_S_c+ABZ*I_ERI_F2xy_Py_Py_S_c;
  Double I_ERI_F2xz_Dyz_Py_S_c = I_ERI_G2x2z_Py_Py_S_c+ABZ*I_ERI_F2xz_Py_Py_S_c;
  Double I_ERI_Fx2y_Dyz_Py_S_c = I_ERI_Gx2yz_Py_Py_S_c+ABZ*I_ERI_Fx2y_Py_Py_S_c;
  Double I_ERI_Fxyz_Dyz_Py_S_c = I_ERI_Gxy2z_Py_Py_S_c+ABZ*I_ERI_Fxyz_Py_Py_S_c;
  Double I_ERI_Fx2z_Dyz_Py_S_c = I_ERI_Gx3z_Py_Py_S_c+ABZ*I_ERI_Fx2z_Py_Py_S_c;
  Double I_ERI_F3y_Dyz_Py_S_c = I_ERI_G3yz_Py_Py_S_c+ABZ*I_ERI_F3y_Py_Py_S_c;
  Double I_ERI_F2yz_Dyz_Py_S_c = I_ERI_G2y2z_Py_Py_S_c+ABZ*I_ERI_F2yz_Py_Py_S_c;
  Double I_ERI_Fy2z_Dyz_Py_S_c = I_ERI_Gy3z_Py_Py_S_c+ABZ*I_ERI_Fy2z_Py_Py_S_c;
  Double I_ERI_F3z_Dyz_Py_S_c = I_ERI_G4z_Py_Py_S_c+ABZ*I_ERI_F3z_Py_Py_S_c;
  Double I_ERI_F3x_D2z_Py_S_c = I_ERI_G3xz_Pz_Py_S_c+ABZ*I_ERI_F3x_Pz_Py_S_c;
  Double I_ERI_F2xy_D2z_Py_S_c = I_ERI_G2xyz_Pz_Py_S_c+ABZ*I_ERI_F2xy_Pz_Py_S_c;
  Double I_ERI_F2xz_D2z_Py_S_c = I_ERI_G2x2z_Pz_Py_S_c+ABZ*I_ERI_F2xz_Pz_Py_S_c;
  Double I_ERI_Fx2y_D2z_Py_S_c = I_ERI_Gx2yz_Pz_Py_S_c+ABZ*I_ERI_Fx2y_Pz_Py_S_c;
  Double I_ERI_Fxyz_D2z_Py_S_c = I_ERI_Gxy2z_Pz_Py_S_c+ABZ*I_ERI_Fxyz_Pz_Py_S_c;
  Double I_ERI_Fx2z_D2z_Py_S_c = I_ERI_Gx3z_Pz_Py_S_c+ABZ*I_ERI_Fx2z_Pz_Py_S_c;
  Double I_ERI_F3y_D2z_Py_S_c = I_ERI_G3yz_Pz_Py_S_c+ABZ*I_ERI_F3y_Pz_Py_S_c;
  Double I_ERI_F2yz_D2z_Py_S_c = I_ERI_G2y2z_Pz_Py_S_c+ABZ*I_ERI_F2yz_Pz_Py_S_c;
  Double I_ERI_Fy2z_D2z_Py_S_c = I_ERI_Gy3z_Pz_Py_S_c+ABZ*I_ERI_Fy2z_Pz_Py_S_c;
  Double I_ERI_F3z_D2z_Py_S_c = I_ERI_G4z_Pz_Py_S_c+ABZ*I_ERI_F3z_Pz_Py_S_c;
  Double I_ERI_F3x_D2x_Pz_S_c = I_ERI_G4x_Px_Pz_S_c+ABX*I_ERI_F3x_Px_Pz_S_c;
  Double I_ERI_F2xy_D2x_Pz_S_c = I_ERI_G3xy_Px_Pz_S_c+ABX*I_ERI_F2xy_Px_Pz_S_c;
  Double I_ERI_F2xz_D2x_Pz_S_c = I_ERI_G3xz_Px_Pz_S_c+ABX*I_ERI_F2xz_Px_Pz_S_c;
  Double I_ERI_Fx2y_D2x_Pz_S_c = I_ERI_G2x2y_Px_Pz_S_c+ABX*I_ERI_Fx2y_Px_Pz_S_c;
  Double I_ERI_Fxyz_D2x_Pz_S_c = I_ERI_G2xyz_Px_Pz_S_c+ABX*I_ERI_Fxyz_Px_Pz_S_c;
  Double I_ERI_Fx2z_D2x_Pz_S_c = I_ERI_G2x2z_Px_Pz_S_c+ABX*I_ERI_Fx2z_Px_Pz_S_c;
  Double I_ERI_F3y_D2x_Pz_S_c = I_ERI_Gx3y_Px_Pz_S_c+ABX*I_ERI_F3y_Px_Pz_S_c;
  Double I_ERI_F2yz_D2x_Pz_S_c = I_ERI_Gx2yz_Px_Pz_S_c+ABX*I_ERI_F2yz_Px_Pz_S_c;
  Double I_ERI_Fy2z_D2x_Pz_S_c = I_ERI_Gxy2z_Px_Pz_S_c+ABX*I_ERI_Fy2z_Px_Pz_S_c;
  Double I_ERI_F3z_D2x_Pz_S_c = I_ERI_Gx3z_Px_Pz_S_c+ABX*I_ERI_F3z_Px_Pz_S_c;
  Double I_ERI_F3x_Dxy_Pz_S_c = I_ERI_G3xy_Px_Pz_S_c+ABY*I_ERI_F3x_Px_Pz_S_c;
  Double I_ERI_F2xy_Dxy_Pz_S_c = I_ERI_G2x2y_Px_Pz_S_c+ABY*I_ERI_F2xy_Px_Pz_S_c;
  Double I_ERI_F2xz_Dxy_Pz_S_c = I_ERI_G2xyz_Px_Pz_S_c+ABY*I_ERI_F2xz_Px_Pz_S_c;
  Double I_ERI_Fx2y_Dxy_Pz_S_c = I_ERI_Gx3y_Px_Pz_S_c+ABY*I_ERI_Fx2y_Px_Pz_S_c;
  Double I_ERI_Fxyz_Dxy_Pz_S_c = I_ERI_Gx2yz_Px_Pz_S_c+ABY*I_ERI_Fxyz_Px_Pz_S_c;
  Double I_ERI_Fx2z_Dxy_Pz_S_c = I_ERI_Gxy2z_Px_Pz_S_c+ABY*I_ERI_Fx2z_Px_Pz_S_c;
  Double I_ERI_F3y_Dxy_Pz_S_c = I_ERI_G4y_Px_Pz_S_c+ABY*I_ERI_F3y_Px_Pz_S_c;
  Double I_ERI_F2yz_Dxy_Pz_S_c = I_ERI_G3yz_Px_Pz_S_c+ABY*I_ERI_F2yz_Px_Pz_S_c;
  Double I_ERI_Fy2z_Dxy_Pz_S_c = I_ERI_G2y2z_Px_Pz_S_c+ABY*I_ERI_Fy2z_Px_Pz_S_c;
  Double I_ERI_F3z_Dxy_Pz_S_c = I_ERI_Gy3z_Px_Pz_S_c+ABY*I_ERI_F3z_Px_Pz_S_c;
  Double I_ERI_F3x_Dxz_Pz_S_c = I_ERI_G3xz_Px_Pz_S_c+ABZ*I_ERI_F3x_Px_Pz_S_c;
  Double I_ERI_F2xy_Dxz_Pz_S_c = I_ERI_G2xyz_Px_Pz_S_c+ABZ*I_ERI_F2xy_Px_Pz_S_c;
  Double I_ERI_F2xz_Dxz_Pz_S_c = I_ERI_G2x2z_Px_Pz_S_c+ABZ*I_ERI_F2xz_Px_Pz_S_c;
  Double I_ERI_Fx2y_Dxz_Pz_S_c = I_ERI_Gx2yz_Px_Pz_S_c+ABZ*I_ERI_Fx2y_Px_Pz_S_c;
  Double I_ERI_Fxyz_Dxz_Pz_S_c = I_ERI_Gxy2z_Px_Pz_S_c+ABZ*I_ERI_Fxyz_Px_Pz_S_c;
  Double I_ERI_Fx2z_Dxz_Pz_S_c = I_ERI_Gx3z_Px_Pz_S_c+ABZ*I_ERI_Fx2z_Px_Pz_S_c;
  Double I_ERI_F3y_Dxz_Pz_S_c = I_ERI_G3yz_Px_Pz_S_c+ABZ*I_ERI_F3y_Px_Pz_S_c;
  Double I_ERI_F2yz_Dxz_Pz_S_c = I_ERI_G2y2z_Px_Pz_S_c+ABZ*I_ERI_F2yz_Px_Pz_S_c;
  Double I_ERI_Fy2z_Dxz_Pz_S_c = I_ERI_Gy3z_Px_Pz_S_c+ABZ*I_ERI_Fy2z_Px_Pz_S_c;
  Double I_ERI_F3z_Dxz_Pz_S_c = I_ERI_G4z_Px_Pz_S_c+ABZ*I_ERI_F3z_Px_Pz_S_c;
  Double I_ERI_F3x_D2y_Pz_S_c = I_ERI_G3xy_Py_Pz_S_c+ABY*I_ERI_F3x_Py_Pz_S_c;
  Double I_ERI_F2xy_D2y_Pz_S_c = I_ERI_G2x2y_Py_Pz_S_c+ABY*I_ERI_F2xy_Py_Pz_S_c;
  Double I_ERI_F2xz_D2y_Pz_S_c = I_ERI_G2xyz_Py_Pz_S_c+ABY*I_ERI_F2xz_Py_Pz_S_c;
  Double I_ERI_Fx2y_D2y_Pz_S_c = I_ERI_Gx3y_Py_Pz_S_c+ABY*I_ERI_Fx2y_Py_Pz_S_c;
  Double I_ERI_Fxyz_D2y_Pz_S_c = I_ERI_Gx2yz_Py_Pz_S_c+ABY*I_ERI_Fxyz_Py_Pz_S_c;
  Double I_ERI_Fx2z_D2y_Pz_S_c = I_ERI_Gxy2z_Py_Pz_S_c+ABY*I_ERI_Fx2z_Py_Pz_S_c;
  Double I_ERI_F3y_D2y_Pz_S_c = I_ERI_G4y_Py_Pz_S_c+ABY*I_ERI_F3y_Py_Pz_S_c;
  Double I_ERI_F2yz_D2y_Pz_S_c = I_ERI_G3yz_Py_Pz_S_c+ABY*I_ERI_F2yz_Py_Pz_S_c;
  Double I_ERI_Fy2z_D2y_Pz_S_c = I_ERI_G2y2z_Py_Pz_S_c+ABY*I_ERI_Fy2z_Py_Pz_S_c;
  Double I_ERI_F3z_D2y_Pz_S_c = I_ERI_Gy3z_Py_Pz_S_c+ABY*I_ERI_F3z_Py_Pz_S_c;
  Double I_ERI_F3x_Dyz_Pz_S_c = I_ERI_G3xz_Py_Pz_S_c+ABZ*I_ERI_F3x_Py_Pz_S_c;
  Double I_ERI_F2xy_Dyz_Pz_S_c = I_ERI_G2xyz_Py_Pz_S_c+ABZ*I_ERI_F2xy_Py_Pz_S_c;
  Double I_ERI_F2xz_Dyz_Pz_S_c = I_ERI_G2x2z_Py_Pz_S_c+ABZ*I_ERI_F2xz_Py_Pz_S_c;
  Double I_ERI_Fx2y_Dyz_Pz_S_c = I_ERI_Gx2yz_Py_Pz_S_c+ABZ*I_ERI_Fx2y_Py_Pz_S_c;
  Double I_ERI_Fxyz_Dyz_Pz_S_c = I_ERI_Gxy2z_Py_Pz_S_c+ABZ*I_ERI_Fxyz_Py_Pz_S_c;
  Double I_ERI_Fx2z_Dyz_Pz_S_c = I_ERI_Gx3z_Py_Pz_S_c+ABZ*I_ERI_Fx2z_Py_Pz_S_c;
  Double I_ERI_F3y_Dyz_Pz_S_c = I_ERI_G3yz_Py_Pz_S_c+ABZ*I_ERI_F3y_Py_Pz_S_c;
  Double I_ERI_F2yz_Dyz_Pz_S_c = I_ERI_G2y2z_Py_Pz_S_c+ABZ*I_ERI_F2yz_Py_Pz_S_c;
  Double I_ERI_Fy2z_Dyz_Pz_S_c = I_ERI_Gy3z_Py_Pz_S_c+ABZ*I_ERI_Fy2z_Py_Pz_S_c;
  Double I_ERI_F3z_Dyz_Pz_S_c = I_ERI_G4z_Py_Pz_S_c+ABZ*I_ERI_F3z_Py_Pz_S_c;
  Double I_ERI_F3x_D2z_Pz_S_c = I_ERI_G3xz_Pz_Pz_S_c+ABZ*I_ERI_F3x_Pz_Pz_S_c;
  Double I_ERI_F2xy_D2z_Pz_S_c = I_ERI_G2xyz_Pz_Pz_S_c+ABZ*I_ERI_F2xy_Pz_Pz_S_c;
  Double I_ERI_F2xz_D2z_Pz_S_c = I_ERI_G2x2z_Pz_Pz_S_c+ABZ*I_ERI_F2xz_Pz_Pz_S_c;
  Double I_ERI_Fx2y_D2z_Pz_S_c = I_ERI_Gx2yz_Pz_Pz_S_c+ABZ*I_ERI_Fx2y_Pz_Pz_S_c;
  Double I_ERI_Fxyz_D2z_Pz_S_c = I_ERI_Gxy2z_Pz_Pz_S_c+ABZ*I_ERI_Fxyz_Pz_Pz_S_c;
  Double I_ERI_Fx2z_D2z_Pz_S_c = I_ERI_Gx3z_Pz_Pz_S_c+ABZ*I_ERI_Fx2z_Pz_Pz_S_c;
  Double I_ERI_F3y_D2z_Pz_S_c = I_ERI_G3yz_Pz_Pz_S_c+ABZ*I_ERI_F3y_Pz_Pz_S_c;
  Double I_ERI_F2yz_D2z_Pz_S_c = I_ERI_G2y2z_Pz_Pz_S_c+ABZ*I_ERI_F2yz_Pz_Pz_S_c;
  Double I_ERI_Fy2z_D2z_Pz_S_c = I_ERI_Gy3z_Pz_Pz_S_c+ABZ*I_ERI_Fy2z_Pz_Pz_S_c;
  Double I_ERI_F3z_D2z_Pz_S_c = I_ERI_G4z_Pz_Pz_S_c+ABZ*I_ERI_F3z_Pz_Pz_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_S_S_a
   * RHS shell quartet name: SQ_ERI_D_D_S_S
   ************************************************************/
  abcd[0] = 2.0E0*I_ERI_G4x_D2x_S_S_a-3*I_ERI_D2x_D2x_S_S;
  abcd[1] = 2.0E0*I_ERI_G3xy_D2x_S_S_a-2*I_ERI_Dxy_D2x_S_S;
  abcd[2] = 2.0E0*I_ERI_G3xz_D2x_S_S_a-2*I_ERI_Dxz_D2x_S_S;
  abcd[3] = 2.0E0*I_ERI_G2x2y_D2x_S_S_a-1*I_ERI_D2y_D2x_S_S;
  abcd[4] = 2.0E0*I_ERI_G2xyz_D2x_S_S_a-1*I_ERI_Dyz_D2x_S_S;
  abcd[5] = 2.0E0*I_ERI_G2x2z_D2x_S_S_a-1*I_ERI_D2z_D2x_S_S;
  abcd[6] = 2.0E0*I_ERI_Gx3y_D2x_S_S_a;
  abcd[7] = 2.0E0*I_ERI_Gx2yz_D2x_S_S_a;
  abcd[8] = 2.0E0*I_ERI_Gxy2z_D2x_S_S_a;
  abcd[9] = 2.0E0*I_ERI_Gx3z_D2x_S_S_a;
  abcd[10] = 2.0E0*I_ERI_G4x_Dxy_S_S_a-3*I_ERI_D2x_Dxy_S_S;
  abcd[11] = 2.0E0*I_ERI_G3xy_Dxy_S_S_a-2*I_ERI_Dxy_Dxy_S_S;
  abcd[12] = 2.0E0*I_ERI_G3xz_Dxy_S_S_a-2*I_ERI_Dxz_Dxy_S_S;
  abcd[13] = 2.0E0*I_ERI_G2x2y_Dxy_S_S_a-1*I_ERI_D2y_Dxy_S_S;
  abcd[14] = 2.0E0*I_ERI_G2xyz_Dxy_S_S_a-1*I_ERI_Dyz_Dxy_S_S;
  abcd[15] = 2.0E0*I_ERI_G2x2z_Dxy_S_S_a-1*I_ERI_D2z_Dxy_S_S;
  abcd[16] = 2.0E0*I_ERI_Gx3y_Dxy_S_S_a;
  abcd[17] = 2.0E0*I_ERI_Gx2yz_Dxy_S_S_a;
  abcd[18] = 2.0E0*I_ERI_Gxy2z_Dxy_S_S_a;
  abcd[19] = 2.0E0*I_ERI_Gx3z_Dxy_S_S_a;
  abcd[20] = 2.0E0*I_ERI_G4x_Dxz_S_S_a-3*I_ERI_D2x_Dxz_S_S;
  abcd[21] = 2.0E0*I_ERI_G3xy_Dxz_S_S_a-2*I_ERI_Dxy_Dxz_S_S;
  abcd[22] = 2.0E0*I_ERI_G3xz_Dxz_S_S_a-2*I_ERI_Dxz_Dxz_S_S;
  abcd[23] = 2.0E0*I_ERI_G2x2y_Dxz_S_S_a-1*I_ERI_D2y_Dxz_S_S;
  abcd[24] = 2.0E0*I_ERI_G2xyz_Dxz_S_S_a-1*I_ERI_Dyz_Dxz_S_S;
  abcd[25] = 2.0E0*I_ERI_G2x2z_Dxz_S_S_a-1*I_ERI_D2z_Dxz_S_S;
  abcd[26] = 2.0E0*I_ERI_Gx3y_Dxz_S_S_a;
  abcd[27] = 2.0E0*I_ERI_Gx2yz_Dxz_S_S_a;
  abcd[28] = 2.0E0*I_ERI_Gxy2z_Dxz_S_S_a;
  abcd[29] = 2.0E0*I_ERI_Gx3z_Dxz_S_S_a;
  abcd[30] = 2.0E0*I_ERI_G4x_D2y_S_S_a-3*I_ERI_D2x_D2y_S_S;
  abcd[31] = 2.0E0*I_ERI_G3xy_D2y_S_S_a-2*I_ERI_Dxy_D2y_S_S;
  abcd[32] = 2.0E0*I_ERI_G3xz_D2y_S_S_a-2*I_ERI_Dxz_D2y_S_S;
  abcd[33] = 2.0E0*I_ERI_G2x2y_D2y_S_S_a-1*I_ERI_D2y_D2y_S_S;
  abcd[34] = 2.0E0*I_ERI_G2xyz_D2y_S_S_a-1*I_ERI_Dyz_D2y_S_S;
  abcd[35] = 2.0E0*I_ERI_G2x2z_D2y_S_S_a-1*I_ERI_D2z_D2y_S_S;
  abcd[36] = 2.0E0*I_ERI_Gx3y_D2y_S_S_a;
  abcd[37] = 2.0E0*I_ERI_Gx2yz_D2y_S_S_a;
  abcd[38] = 2.0E0*I_ERI_Gxy2z_D2y_S_S_a;
  abcd[39] = 2.0E0*I_ERI_Gx3z_D2y_S_S_a;
  abcd[40] = 2.0E0*I_ERI_G4x_Dyz_S_S_a-3*I_ERI_D2x_Dyz_S_S;
  abcd[41] = 2.0E0*I_ERI_G3xy_Dyz_S_S_a-2*I_ERI_Dxy_Dyz_S_S;
  abcd[42] = 2.0E0*I_ERI_G3xz_Dyz_S_S_a-2*I_ERI_Dxz_Dyz_S_S;
  abcd[43] = 2.0E0*I_ERI_G2x2y_Dyz_S_S_a-1*I_ERI_D2y_Dyz_S_S;
  abcd[44] = 2.0E0*I_ERI_G2xyz_Dyz_S_S_a-1*I_ERI_Dyz_Dyz_S_S;
  abcd[45] = 2.0E0*I_ERI_G2x2z_Dyz_S_S_a-1*I_ERI_D2z_Dyz_S_S;
  abcd[46] = 2.0E0*I_ERI_Gx3y_Dyz_S_S_a;
  abcd[47] = 2.0E0*I_ERI_Gx2yz_Dyz_S_S_a;
  abcd[48] = 2.0E0*I_ERI_Gxy2z_Dyz_S_S_a;
  abcd[49] = 2.0E0*I_ERI_Gx3z_Dyz_S_S_a;
  abcd[50] = 2.0E0*I_ERI_G4x_D2z_S_S_a-3*I_ERI_D2x_D2z_S_S;
  abcd[51] = 2.0E0*I_ERI_G3xy_D2z_S_S_a-2*I_ERI_Dxy_D2z_S_S;
  abcd[52] = 2.0E0*I_ERI_G3xz_D2z_S_S_a-2*I_ERI_Dxz_D2z_S_S;
  abcd[53] = 2.0E0*I_ERI_G2x2y_D2z_S_S_a-1*I_ERI_D2y_D2z_S_S;
  abcd[54] = 2.0E0*I_ERI_G2xyz_D2z_S_S_a-1*I_ERI_Dyz_D2z_S_S;
  abcd[55] = 2.0E0*I_ERI_G2x2z_D2z_S_S_a-1*I_ERI_D2z_D2z_S_S;
  abcd[56] = 2.0E0*I_ERI_Gx3y_D2z_S_S_a;
  abcd[57] = 2.0E0*I_ERI_Gx2yz_D2z_S_S_a;
  abcd[58] = 2.0E0*I_ERI_Gxy2z_D2z_S_S_a;
  abcd[59] = 2.0E0*I_ERI_Gx3z_D2z_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_S_S_a
   * RHS shell quartet name: SQ_ERI_D_D_S_S
   ************************************************************/
  abcd[60] = 2.0E0*I_ERI_G3xy_D2x_S_S_a;
  abcd[61] = 2.0E0*I_ERI_G2x2y_D2x_S_S_a-1*I_ERI_D2x_D2x_S_S;
  abcd[62] = 2.0E0*I_ERI_G2xyz_D2x_S_S_a;
  abcd[63] = 2.0E0*I_ERI_Gx3y_D2x_S_S_a-2*I_ERI_Dxy_D2x_S_S;
  abcd[64] = 2.0E0*I_ERI_Gx2yz_D2x_S_S_a-1*I_ERI_Dxz_D2x_S_S;
  abcd[65] = 2.0E0*I_ERI_Gxy2z_D2x_S_S_a;
  abcd[66] = 2.0E0*I_ERI_G4y_D2x_S_S_a-3*I_ERI_D2y_D2x_S_S;
  abcd[67] = 2.0E0*I_ERI_G3yz_D2x_S_S_a-2*I_ERI_Dyz_D2x_S_S;
  abcd[68] = 2.0E0*I_ERI_G2y2z_D2x_S_S_a-1*I_ERI_D2z_D2x_S_S;
  abcd[69] = 2.0E0*I_ERI_Gy3z_D2x_S_S_a;
  abcd[70] = 2.0E0*I_ERI_G3xy_Dxy_S_S_a;
  abcd[71] = 2.0E0*I_ERI_G2x2y_Dxy_S_S_a-1*I_ERI_D2x_Dxy_S_S;
  abcd[72] = 2.0E0*I_ERI_G2xyz_Dxy_S_S_a;
  abcd[73] = 2.0E0*I_ERI_Gx3y_Dxy_S_S_a-2*I_ERI_Dxy_Dxy_S_S;
  abcd[74] = 2.0E0*I_ERI_Gx2yz_Dxy_S_S_a-1*I_ERI_Dxz_Dxy_S_S;
  abcd[75] = 2.0E0*I_ERI_Gxy2z_Dxy_S_S_a;
  abcd[76] = 2.0E0*I_ERI_G4y_Dxy_S_S_a-3*I_ERI_D2y_Dxy_S_S;
  abcd[77] = 2.0E0*I_ERI_G3yz_Dxy_S_S_a-2*I_ERI_Dyz_Dxy_S_S;
  abcd[78] = 2.0E0*I_ERI_G2y2z_Dxy_S_S_a-1*I_ERI_D2z_Dxy_S_S;
  abcd[79] = 2.0E0*I_ERI_Gy3z_Dxy_S_S_a;
  abcd[80] = 2.0E0*I_ERI_G3xy_Dxz_S_S_a;
  abcd[81] = 2.0E0*I_ERI_G2x2y_Dxz_S_S_a-1*I_ERI_D2x_Dxz_S_S;
  abcd[82] = 2.0E0*I_ERI_G2xyz_Dxz_S_S_a;
  abcd[83] = 2.0E0*I_ERI_Gx3y_Dxz_S_S_a-2*I_ERI_Dxy_Dxz_S_S;
  abcd[84] = 2.0E0*I_ERI_Gx2yz_Dxz_S_S_a-1*I_ERI_Dxz_Dxz_S_S;
  abcd[85] = 2.0E0*I_ERI_Gxy2z_Dxz_S_S_a;
  abcd[86] = 2.0E0*I_ERI_G4y_Dxz_S_S_a-3*I_ERI_D2y_Dxz_S_S;
  abcd[87] = 2.0E0*I_ERI_G3yz_Dxz_S_S_a-2*I_ERI_Dyz_Dxz_S_S;
  abcd[88] = 2.0E0*I_ERI_G2y2z_Dxz_S_S_a-1*I_ERI_D2z_Dxz_S_S;
  abcd[89] = 2.0E0*I_ERI_Gy3z_Dxz_S_S_a;
  abcd[90] = 2.0E0*I_ERI_G3xy_D2y_S_S_a;
  abcd[91] = 2.0E0*I_ERI_G2x2y_D2y_S_S_a-1*I_ERI_D2x_D2y_S_S;
  abcd[92] = 2.0E0*I_ERI_G2xyz_D2y_S_S_a;
  abcd[93] = 2.0E0*I_ERI_Gx3y_D2y_S_S_a-2*I_ERI_Dxy_D2y_S_S;
  abcd[94] = 2.0E0*I_ERI_Gx2yz_D2y_S_S_a-1*I_ERI_Dxz_D2y_S_S;
  abcd[95] = 2.0E0*I_ERI_Gxy2z_D2y_S_S_a;
  abcd[96] = 2.0E0*I_ERI_G4y_D2y_S_S_a-3*I_ERI_D2y_D2y_S_S;
  abcd[97] = 2.0E0*I_ERI_G3yz_D2y_S_S_a-2*I_ERI_Dyz_D2y_S_S;
  abcd[98] = 2.0E0*I_ERI_G2y2z_D2y_S_S_a-1*I_ERI_D2z_D2y_S_S;
  abcd[99] = 2.0E0*I_ERI_Gy3z_D2y_S_S_a;
  abcd[100] = 2.0E0*I_ERI_G3xy_Dyz_S_S_a;
  abcd[101] = 2.0E0*I_ERI_G2x2y_Dyz_S_S_a-1*I_ERI_D2x_Dyz_S_S;
  abcd[102] = 2.0E0*I_ERI_G2xyz_Dyz_S_S_a;
  abcd[103] = 2.0E0*I_ERI_Gx3y_Dyz_S_S_a-2*I_ERI_Dxy_Dyz_S_S;
  abcd[104] = 2.0E0*I_ERI_Gx2yz_Dyz_S_S_a-1*I_ERI_Dxz_Dyz_S_S;
  abcd[105] = 2.0E0*I_ERI_Gxy2z_Dyz_S_S_a;
  abcd[106] = 2.0E0*I_ERI_G4y_Dyz_S_S_a-3*I_ERI_D2y_Dyz_S_S;
  abcd[107] = 2.0E0*I_ERI_G3yz_Dyz_S_S_a-2*I_ERI_Dyz_Dyz_S_S;
  abcd[108] = 2.0E0*I_ERI_G2y2z_Dyz_S_S_a-1*I_ERI_D2z_Dyz_S_S;
  abcd[109] = 2.0E0*I_ERI_Gy3z_Dyz_S_S_a;
  abcd[110] = 2.0E0*I_ERI_G3xy_D2z_S_S_a;
  abcd[111] = 2.0E0*I_ERI_G2x2y_D2z_S_S_a-1*I_ERI_D2x_D2z_S_S;
  abcd[112] = 2.0E0*I_ERI_G2xyz_D2z_S_S_a;
  abcd[113] = 2.0E0*I_ERI_Gx3y_D2z_S_S_a-2*I_ERI_Dxy_D2z_S_S;
  abcd[114] = 2.0E0*I_ERI_Gx2yz_D2z_S_S_a-1*I_ERI_Dxz_D2z_S_S;
  abcd[115] = 2.0E0*I_ERI_Gxy2z_D2z_S_S_a;
  abcd[116] = 2.0E0*I_ERI_G4y_D2z_S_S_a-3*I_ERI_D2y_D2z_S_S;
  abcd[117] = 2.0E0*I_ERI_G3yz_D2z_S_S_a-2*I_ERI_Dyz_D2z_S_S;
  abcd[118] = 2.0E0*I_ERI_G2y2z_D2z_S_S_a-1*I_ERI_D2z_D2z_S_S;
  abcd[119] = 2.0E0*I_ERI_Gy3z_D2z_S_S_a;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_D_S_S_a
   * RHS shell quartet name: SQ_ERI_D_D_S_S
   ************************************************************/
  abcd[120] = 2.0E0*I_ERI_G3xz_D2x_S_S_a;
  abcd[121] = 2.0E0*I_ERI_G2xyz_D2x_S_S_a;
  abcd[122] = 2.0E0*I_ERI_G2x2z_D2x_S_S_a-1*I_ERI_D2x_D2x_S_S;
  abcd[123] = 2.0E0*I_ERI_Gx2yz_D2x_S_S_a;
  abcd[124] = 2.0E0*I_ERI_Gxy2z_D2x_S_S_a-1*I_ERI_Dxy_D2x_S_S;
  abcd[125] = 2.0E0*I_ERI_Gx3z_D2x_S_S_a-2*I_ERI_Dxz_D2x_S_S;
  abcd[126] = 2.0E0*I_ERI_G3yz_D2x_S_S_a;
  abcd[127] = 2.0E0*I_ERI_G2y2z_D2x_S_S_a-1*I_ERI_D2y_D2x_S_S;
  abcd[128] = 2.0E0*I_ERI_Gy3z_D2x_S_S_a-2*I_ERI_Dyz_D2x_S_S;
  abcd[129] = 2.0E0*I_ERI_G4z_D2x_S_S_a-3*I_ERI_D2z_D2x_S_S;
  abcd[130] = 2.0E0*I_ERI_G3xz_Dxy_S_S_a;
  abcd[131] = 2.0E0*I_ERI_G2xyz_Dxy_S_S_a;
  abcd[132] = 2.0E0*I_ERI_G2x2z_Dxy_S_S_a-1*I_ERI_D2x_Dxy_S_S;
  abcd[133] = 2.0E0*I_ERI_Gx2yz_Dxy_S_S_a;
  abcd[134] = 2.0E0*I_ERI_Gxy2z_Dxy_S_S_a-1*I_ERI_Dxy_Dxy_S_S;
  abcd[135] = 2.0E0*I_ERI_Gx3z_Dxy_S_S_a-2*I_ERI_Dxz_Dxy_S_S;
  abcd[136] = 2.0E0*I_ERI_G3yz_Dxy_S_S_a;
  abcd[137] = 2.0E0*I_ERI_G2y2z_Dxy_S_S_a-1*I_ERI_D2y_Dxy_S_S;
  abcd[138] = 2.0E0*I_ERI_Gy3z_Dxy_S_S_a-2*I_ERI_Dyz_Dxy_S_S;
  abcd[139] = 2.0E0*I_ERI_G4z_Dxy_S_S_a-3*I_ERI_D2z_Dxy_S_S;
  abcd[140] = 2.0E0*I_ERI_G3xz_Dxz_S_S_a;
  abcd[141] = 2.0E0*I_ERI_G2xyz_Dxz_S_S_a;
  abcd[142] = 2.0E0*I_ERI_G2x2z_Dxz_S_S_a-1*I_ERI_D2x_Dxz_S_S;
  abcd[143] = 2.0E0*I_ERI_Gx2yz_Dxz_S_S_a;
  abcd[144] = 2.0E0*I_ERI_Gxy2z_Dxz_S_S_a-1*I_ERI_Dxy_Dxz_S_S;
  abcd[145] = 2.0E0*I_ERI_Gx3z_Dxz_S_S_a-2*I_ERI_Dxz_Dxz_S_S;
  abcd[146] = 2.0E0*I_ERI_G3yz_Dxz_S_S_a;
  abcd[147] = 2.0E0*I_ERI_G2y2z_Dxz_S_S_a-1*I_ERI_D2y_Dxz_S_S;
  abcd[148] = 2.0E0*I_ERI_Gy3z_Dxz_S_S_a-2*I_ERI_Dyz_Dxz_S_S;
  abcd[149] = 2.0E0*I_ERI_G4z_Dxz_S_S_a-3*I_ERI_D2z_Dxz_S_S;
  abcd[150] = 2.0E0*I_ERI_G3xz_D2y_S_S_a;
  abcd[151] = 2.0E0*I_ERI_G2xyz_D2y_S_S_a;
  abcd[152] = 2.0E0*I_ERI_G2x2z_D2y_S_S_a-1*I_ERI_D2x_D2y_S_S;
  abcd[153] = 2.0E0*I_ERI_Gx2yz_D2y_S_S_a;
  abcd[154] = 2.0E0*I_ERI_Gxy2z_D2y_S_S_a-1*I_ERI_Dxy_D2y_S_S;
  abcd[155] = 2.0E0*I_ERI_Gx3z_D2y_S_S_a-2*I_ERI_Dxz_D2y_S_S;
  abcd[156] = 2.0E0*I_ERI_G3yz_D2y_S_S_a;
  abcd[157] = 2.0E0*I_ERI_G2y2z_D2y_S_S_a-1*I_ERI_D2y_D2y_S_S;
  abcd[158] = 2.0E0*I_ERI_Gy3z_D2y_S_S_a-2*I_ERI_Dyz_D2y_S_S;
  abcd[159] = 2.0E0*I_ERI_G4z_D2y_S_S_a-3*I_ERI_D2z_D2y_S_S;
  abcd[160] = 2.0E0*I_ERI_G3xz_Dyz_S_S_a;
  abcd[161] = 2.0E0*I_ERI_G2xyz_Dyz_S_S_a;
  abcd[162] = 2.0E0*I_ERI_G2x2z_Dyz_S_S_a-1*I_ERI_D2x_Dyz_S_S;
  abcd[163] = 2.0E0*I_ERI_Gx2yz_Dyz_S_S_a;
  abcd[164] = 2.0E0*I_ERI_Gxy2z_Dyz_S_S_a-1*I_ERI_Dxy_Dyz_S_S;
  abcd[165] = 2.0E0*I_ERI_Gx3z_Dyz_S_S_a-2*I_ERI_Dxz_Dyz_S_S;
  abcd[166] = 2.0E0*I_ERI_G3yz_Dyz_S_S_a;
  abcd[167] = 2.0E0*I_ERI_G2y2z_Dyz_S_S_a-1*I_ERI_D2y_Dyz_S_S;
  abcd[168] = 2.0E0*I_ERI_Gy3z_Dyz_S_S_a-2*I_ERI_Dyz_Dyz_S_S;
  abcd[169] = 2.0E0*I_ERI_G4z_Dyz_S_S_a-3*I_ERI_D2z_Dyz_S_S;
  abcd[170] = 2.0E0*I_ERI_G3xz_D2z_S_S_a;
  abcd[171] = 2.0E0*I_ERI_G2xyz_D2z_S_S_a;
  abcd[172] = 2.0E0*I_ERI_G2x2z_D2z_S_S_a-1*I_ERI_D2x_D2z_S_S;
  abcd[173] = 2.0E0*I_ERI_Gx2yz_D2z_S_S_a;
  abcd[174] = 2.0E0*I_ERI_Gxy2z_D2z_S_S_a-1*I_ERI_Dxy_D2z_S_S;
  abcd[175] = 2.0E0*I_ERI_Gx3z_D2z_S_S_a-2*I_ERI_Dxz_D2z_S_S;
  abcd[176] = 2.0E0*I_ERI_G3yz_D2z_S_S_a;
  abcd[177] = 2.0E0*I_ERI_G2y2z_D2z_S_S_a-1*I_ERI_D2y_D2z_S_S;
  abcd[178] = 2.0E0*I_ERI_Gy3z_D2z_S_S_a-2*I_ERI_Dyz_D2z_S_S;
  abcd[179] = 2.0E0*I_ERI_G4z_D2z_S_S_a-3*I_ERI_D2z_D2z_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_S_S_b
   * RHS shell quartet name: SQ_ERI_F_P_S_S
   ************************************************************/
  abcd[180] = 2.0E0*I_ERI_F3x_F3x_S_S_b-2*I_ERI_F3x_Px_S_S;
  abcd[181] = 2.0E0*I_ERI_F2xy_F3x_S_S_b-2*I_ERI_F2xy_Px_S_S;
  abcd[182] = 2.0E0*I_ERI_F2xz_F3x_S_S_b-2*I_ERI_F2xz_Px_S_S;
  abcd[183] = 2.0E0*I_ERI_Fx2y_F3x_S_S_b-2*I_ERI_Fx2y_Px_S_S;
  abcd[184] = 2.0E0*I_ERI_Fxyz_F3x_S_S_b-2*I_ERI_Fxyz_Px_S_S;
  abcd[185] = 2.0E0*I_ERI_Fx2z_F3x_S_S_b-2*I_ERI_Fx2z_Px_S_S;
  abcd[186] = 2.0E0*I_ERI_F3y_F3x_S_S_b-2*I_ERI_F3y_Px_S_S;
  abcd[187] = 2.0E0*I_ERI_F2yz_F3x_S_S_b-2*I_ERI_F2yz_Px_S_S;
  abcd[188] = 2.0E0*I_ERI_Fy2z_F3x_S_S_b-2*I_ERI_Fy2z_Px_S_S;
  abcd[189] = 2.0E0*I_ERI_F3z_F3x_S_S_b-2*I_ERI_F3z_Px_S_S;
  abcd[190] = 2.0E0*I_ERI_F3x_F2xy_S_S_b-1*I_ERI_F3x_Py_S_S;
  abcd[191] = 2.0E0*I_ERI_F2xy_F2xy_S_S_b-1*I_ERI_F2xy_Py_S_S;
  abcd[192] = 2.0E0*I_ERI_F2xz_F2xy_S_S_b-1*I_ERI_F2xz_Py_S_S;
  abcd[193] = 2.0E0*I_ERI_Fx2y_F2xy_S_S_b-1*I_ERI_Fx2y_Py_S_S;
  abcd[194] = 2.0E0*I_ERI_Fxyz_F2xy_S_S_b-1*I_ERI_Fxyz_Py_S_S;
  abcd[195] = 2.0E0*I_ERI_Fx2z_F2xy_S_S_b-1*I_ERI_Fx2z_Py_S_S;
  abcd[196] = 2.0E0*I_ERI_F3y_F2xy_S_S_b-1*I_ERI_F3y_Py_S_S;
  abcd[197] = 2.0E0*I_ERI_F2yz_F2xy_S_S_b-1*I_ERI_F2yz_Py_S_S;
  abcd[198] = 2.0E0*I_ERI_Fy2z_F2xy_S_S_b-1*I_ERI_Fy2z_Py_S_S;
  abcd[199] = 2.0E0*I_ERI_F3z_F2xy_S_S_b-1*I_ERI_F3z_Py_S_S;
  abcd[200] = 2.0E0*I_ERI_F3x_F2xz_S_S_b-1*I_ERI_F3x_Pz_S_S;
  abcd[201] = 2.0E0*I_ERI_F2xy_F2xz_S_S_b-1*I_ERI_F2xy_Pz_S_S;
  abcd[202] = 2.0E0*I_ERI_F2xz_F2xz_S_S_b-1*I_ERI_F2xz_Pz_S_S;
  abcd[203] = 2.0E0*I_ERI_Fx2y_F2xz_S_S_b-1*I_ERI_Fx2y_Pz_S_S;
  abcd[204] = 2.0E0*I_ERI_Fxyz_F2xz_S_S_b-1*I_ERI_Fxyz_Pz_S_S;
  abcd[205] = 2.0E0*I_ERI_Fx2z_F2xz_S_S_b-1*I_ERI_Fx2z_Pz_S_S;
  abcd[206] = 2.0E0*I_ERI_F3y_F2xz_S_S_b-1*I_ERI_F3y_Pz_S_S;
  abcd[207] = 2.0E0*I_ERI_F2yz_F2xz_S_S_b-1*I_ERI_F2yz_Pz_S_S;
  abcd[208] = 2.0E0*I_ERI_Fy2z_F2xz_S_S_b-1*I_ERI_Fy2z_Pz_S_S;
  abcd[209] = 2.0E0*I_ERI_F3z_F2xz_S_S_b-1*I_ERI_F3z_Pz_S_S;
  abcd[210] = 2.0E0*I_ERI_F3x_Fx2y_S_S_b;
  abcd[211] = 2.0E0*I_ERI_F2xy_Fx2y_S_S_b;
  abcd[212] = 2.0E0*I_ERI_F2xz_Fx2y_S_S_b;
  abcd[213] = 2.0E0*I_ERI_Fx2y_Fx2y_S_S_b;
  abcd[214] = 2.0E0*I_ERI_Fxyz_Fx2y_S_S_b;
  abcd[215] = 2.0E0*I_ERI_Fx2z_Fx2y_S_S_b;
  abcd[216] = 2.0E0*I_ERI_F3y_Fx2y_S_S_b;
  abcd[217] = 2.0E0*I_ERI_F2yz_Fx2y_S_S_b;
  abcd[218] = 2.0E0*I_ERI_Fy2z_Fx2y_S_S_b;
  abcd[219] = 2.0E0*I_ERI_F3z_Fx2y_S_S_b;
  abcd[220] = 2.0E0*I_ERI_F3x_Fxyz_S_S_b;
  abcd[221] = 2.0E0*I_ERI_F2xy_Fxyz_S_S_b;
  abcd[222] = 2.0E0*I_ERI_F2xz_Fxyz_S_S_b;
  abcd[223] = 2.0E0*I_ERI_Fx2y_Fxyz_S_S_b;
  abcd[224] = 2.0E0*I_ERI_Fxyz_Fxyz_S_S_b;
  abcd[225] = 2.0E0*I_ERI_Fx2z_Fxyz_S_S_b;
  abcd[226] = 2.0E0*I_ERI_F3y_Fxyz_S_S_b;
  abcd[227] = 2.0E0*I_ERI_F2yz_Fxyz_S_S_b;
  abcd[228] = 2.0E0*I_ERI_Fy2z_Fxyz_S_S_b;
  abcd[229] = 2.0E0*I_ERI_F3z_Fxyz_S_S_b;
  abcd[230] = 2.0E0*I_ERI_F3x_Fx2z_S_S_b;
  abcd[231] = 2.0E0*I_ERI_F2xy_Fx2z_S_S_b;
  abcd[232] = 2.0E0*I_ERI_F2xz_Fx2z_S_S_b;
  abcd[233] = 2.0E0*I_ERI_Fx2y_Fx2z_S_S_b;
  abcd[234] = 2.0E0*I_ERI_Fxyz_Fx2z_S_S_b;
  abcd[235] = 2.0E0*I_ERI_Fx2z_Fx2z_S_S_b;
  abcd[236] = 2.0E0*I_ERI_F3y_Fx2z_S_S_b;
  abcd[237] = 2.0E0*I_ERI_F2yz_Fx2z_S_S_b;
  abcd[238] = 2.0E0*I_ERI_Fy2z_Fx2z_S_S_b;
  abcd[239] = 2.0E0*I_ERI_F3z_Fx2z_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_S_S_b
   * RHS shell quartet name: SQ_ERI_F_P_S_S
   ************************************************************/
  abcd[240] = 2.0E0*I_ERI_F3x_F2xy_S_S_b;
  abcd[241] = 2.0E0*I_ERI_F2xy_F2xy_S_S_b;
  abcd[242] = 2.0E0*I_ERI_F2xz_F2xy_S_S_b;
  abcd[243] = 2.0E0*I_ERI_Fx2y_F2xy_S_S_b;
  abcd[244] = 2.0E0*I_ERI_Fxyz_F2xy_S_S_b;
  abcd[245] = 2.0E0*I_ERI_Fx2z_F2xy_S_S_b;
  abcd[246] = 2.0E0*I_ERI_F3y_F2xy_S_S_b;
  abcd[247] = 2.0E0*I_ERI_F2yz_F2xy_S_S_b;
  abcd[248] = 2.0E0*I_ERI_Fy2z_F2xy_S_S_b;
  abcd[249] = 2.0E0*I_ERI_F3z_F2xy_S_S_b;
  abcd[250] = 2.0E0*I_ERI_F3x_Fx2y_S_S_b-1*I_ERI_F3x_Px_S_S;
  abcd[251] = 2.0E0*I_ERI_F2xy_Fx2y_S_S_b-1*I_ERI_F2xy_Px_S_S;
  abcd[252] = 2.0E0*I_ERI_F2xz_Fx2y_S_S_b-1*I_ERI_F2xz_Px_S_S;
  abcd[253] = 2.0E0*I_ERI_Fx2y_Fx2y_S_S_b-1*I_ERI_Fx2y_Px_S_S;
  abcd[254] = 2.0E0*I_ERI_Fxyz_Fx2y_S_S_b-1*I_ERI_Fxyz_Px_S_S;
  abcd[255] = 2.0E0*I_ERI_Fx2z_Fx2y_S_S_b-1*I_ERI_Fx2z_Px_S_S;
  abcd[256] = 2.0E0*I_ERI_F3y_Fx2y_S_S_b-1*I_ERI_F3y_Px_S_S;
  abcd[257] = 2.0E0*I_ERI_F2yz_Fx2y_S_S_b-1*I_ERI_F2yz_Px_S_S;
  abcd[258] = 2.0E0*I_ERI_Fy2z_Fx2y_S_S_b-1*I_ERI_Fy2z_Px_S_S;
  abcd[259] = 2.0E0*I_ERI_F3z_Fx2y_S_S_b-1*I_ERI_F3z_Px_S_S;
  abcd[260] = 2.0E0*I_ERI_F3x_Fxyz_S_S_b;
  abcd[261] = 2.0E0*I_ERI_F2xy_Fxyz_S_S_b;
  abcd[262] = 2.0E0*I_ERI_F2xz_Fxyz_S_S_b;
  abcd[263] = 2.0E0*I_ERI_Fx2y_Fxyz_S_S_b;
  abcd[264] = 2.0E0*I_ERI_Fxyz_Fxyz_S_S_b;
  abcd[265] = 2.0E0*I_ERI_Fx2z_Fxyz_S_S_b;
  abcd[266] = 2.0E0*I_ERI_F3y_Fxyz_S_S_b;
  abcd[267] = 2.0E0*I_ERI_F2yz_Fxyz_S_S_b;
  abcd[268] = 2.0E0*I_ERI_Fy2z_Fxyz_S_S_b;
  abcd[269] = 2.0E0*I_ERI_F3z_Fxyz_S_S_b;
  abcd[270] = 2.0E0*I_ERI_F3x_F3y_S_S_b-2*I_ERI_F3x_Py_S_S;
  abcd[271] = 2.0E0*I_ERI_F2xy_F3y_S_S_b-2*I_ERI_F2xy_Py_S_S;
  abcd[272] = 2.0E0*I_ERI_F2xz_F3y_S_S_b-2*I_ERI_F2xz_Py_S_S;
  abcd[273] = 2.0E0*I_ERI_Fx2y_F3y_S_S_b-2*I_ERI_Fx2y_Py_S_S;
  abcd[274] = 2.0E0*I_ERI_Fxyz_F3y_S_S_b-2*I_ERI_Fxyz_Py_S_S;
  abcd[275] = 2.0E0*I_ERI_Fx2z_F3y_S_S_b-2*I_ERI_Fx2z_Py_S_S;
  abcd[276] = 2.0E0*I_ERI_F3y_F3y_S_S_b-2*I_ERI_F3y_Py_S_S;
  abcd[277] = 2.0E0*I_ERI_F2yz_F3y_S_S_b-2*I_ERI_F2yz_Py_S_S;
  abcd[278] = 2.0E0*I_ERI_Fy2z_F3y_S_S_b-2*I_ERI_Fy2z_Py_S_S;
  abcd[279] = 2.0E0*I_ERI_F3z_F3y_S_S_b-2*I_ERI_F3z_Py_S_S;
  abcd[280] = 2.0E0*I_ERI_F3x_F2yz_S_S_b-1*I_ERI_F3x_Pz_S_S;
  abcd[281] = 2.0E0*I_ERI_F2xy_F2yz_S_S_b-1*I_ERI_F2xy_Pz_S_S;
  abcd[282] = 2.0E0*I_ERI_F2xz_F2yz_S_S_b-1*I_ERI_F2xz_Pz_S_S;
  abcd[283] = 2.0E0*I_ERI_Fx2y_F2yz_S_S_b-1*I_ERI_Fx2y_Pz_S_S;
  abcd[284] = 2.0E0*I_ERI_Fxyz_F2yz_S_S_b-1*I_ERI_Fxyz_Pz_S_S;
  abcd[285] = 2.0E0*I_ERI_Fx2z_F2yz_S_S_b-1*I_ERI_Fx2z_Pz_S_S;
  abcd[286] = 2.0E0*I_ERI_F3y_F2yz_S_S_b-1*I_ERI_F3y_Pz_S_S;
  abcd[287] = 2.0E0*I_ERI_F2yz_F2yz_S_S_b-1*I_ERI_F2yz_Pz_S_S;
  abcd[288] = 2.0E0*I_ERI_Fy2z_F2yz_S_S_b-1*I_ERI_Fy2z_Pz_S_S;
  abcd[289] = 2.0E0*I_ERI_F3z_F2yz_S_S_b-1*I_ERI_F3z_Pz_S_S;
  abcd[290] = 2.0E0*I_ERI_F3x_Fy2z_S_S_b;
  abcd[291] = 2.0E0*I_ERI_F2xy_Fy2z_S_S_b;
  abcd[292] = 2.0E0*I_ERI_F2xz_Fy2z_S_S_b;
  abcd[293] = 2.0E0*I_ERI_Fx2y_Fy2z_S_S_b;
  abcd[294] = 2.0E0*I_ERI_Fxyz_Fy2z_S_S_b;
  abcd[295] = 2.0E0*I_ERI_Fx2z_Fy2z_S_S_b;
  abcd[296] = 2.0E0*I_ERI_F3y_Fy2z_S_S_b;
  abcd[297] = 2.0E0*I_ERI_F2yz_Fy2z_S_S_b;
  abcd[298] = 2.0E0*I_ERI_Fy2z_Fy2z_S_S_b;
  abcd[299] = 2.0E0*I_ERI_F3z_Fy2z_S_S_b;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_F_S_S_b
   * RHS shell quartet name: SQ_ERI_F_P_S_S
   ************************************************************/
  abcd[300] = 2.0E0*I_ERI_F3x_F2xz_S_S_b;
  abcd[301] = 2.0E0*I_ERI_F2xy_F2xz_S_S_b;
  abcd[302] = 2.0E0*I_ERI_F2xz_F2xz_S_S_b;
  abcd[303] = 2.0E0*I_ERI_Fx2y_F2xz_S_S_b;
  abcd[304] = 2.0E0*I_ERI_Fxyz_F2xz_S_S_b;
  abcd[305] = 2.0E0*I_ERI_Fx2z_F2xz_S_S_b;
  abcd[306] = 2.0E0*I_ERI_F3y_F2xz_S_S_b;
  abcd[307] = 2.0E0*I_ERI_F2yz_F2xz_S_S_b;
  abcd[308] = 2.0E0*I_ERI_Fy2z_F2xz_S_S_b;
  abcd[309] = 2.0E0*I_ERI_F3z_F2xz_S_S_b;
  abcd[310] = 2.0E0*I_ERI_F3x_Fxyz_S_S_b;
  abcd[311] = 2.0E0*I_ERI_F2xy_Fxyz_S_S_b;
  abcd[312] = 2.0E0*I_ERI_F2xz_Fxyz_S_S_b;
  abcd[313] = 2.0E0*I_ERI_Fx2y_Fxyz_S_S_b;
  abcd[314] = 2.0E0*I_ERI_Fxyz_Fxyz_S_S_b;
  abcd[315] = 2.0E0*I_ERI_Fx2z_Fxyz_S_S_b;
  abcd[316] = 2.0E0*I_ERI_F3y_Fxyz_S_S_b;
  abcd[317] = 2.0E0*I_ERI_F2yz_Fxyz_S_S_b;
  abcd[318] = 2.0E0*I_ERI_Fy2z_Fxyz_S_S_b;
  abcd[319] = 2.0E0*I_ERI_F3z_Fxyz_S_S_b;
  abcd[320] = 2.0E0*I_ERI_F3x_Fx2z_S_S_b-1*I_ERI_F3x_Px_S_S;
  abcd[321] = 2.0E0*I_ERI_F2xy_Fx2z_S_S_b-1*I_ERI_F2xy_Px_S_S;
  abcd[322] = 2.0E0*I_ERI_F2xz_Fx2z_S_S_b-1*I_ERI_F2xz_Px_S_S;
  abcd[323] = 2.0E0*I_ERI_Fx2y_Fx2z_S_S_b-1*I_ERI_Fx2y_Px_S_S;
  abcd[324] = 2.0E0*I_ERI_Fxyz_Fx2z_S_S_b-1*I_ERI_Fxyz_Px_S_S;
  abcd[325] = 2.0E0*I_ERI_Fx2z_Fx2z_S_S_b-1*I_ERI_Fx2z_Px_S_S;
  abcd[326] = 2.0E0*I_ERI_F3y_Fx2z_S_S_b-1*I_ERI_F3y_Px_S_S;
  abcd[327] = 2.0E0*I_ERI_F2yz_Fx2z_S_S_b-1*I_ERI_F2yz_Px_S_S;
  abcd[328] = 2.0E0*I_ERI_Fy2z_Fx2z_S_S_b-1*I_ERI_Fy2z_Px_S_S;
  abcd[329] = 2.0E0*I_ERI_F3z_Fx2z_S_S_b-1*I_ERI_F3z_Px_S_S;
  abcd[330] = 2.0E0*I_ERI_F3x_F2yz_S_S_b;
  abcd[331] = 2.0E0*I_ERI_F2xy_F2yz_S_S_b;
  abcd[332] = 2.0E0*I_ERI_F2xz_F2yz_S_S_b;
  abcd[333] = 2.0E0*I_ERI_Fx2y_F2yz_S_S_b;
  abcd[334] = 2.0E0*I_ERI_Fxyz_F2yz_S_S_b;
  abcd[335] = 2.0E0*I_ERI_Fx2z_F2yz_S_S_b;
  abcd[336] = 2.0E0*I_ERI_F3y_F2yz_S_S_b;
  abcd[337] = 2.0E0*I_ERI_F2yz_F2yz_S_S_b;
  abcd[338] = 2.0E0*I_ERI_Fy2z_F2yz_S_S_b;
  abcd[339] = 2.0E0*I_ERI_F3z_F2yz_S_S_b;
  abcd[340] = 2.0E0*I_ERI_F3x_Fy2z_S_S_b-1*I_ERI_F3x_Py_S_S;
  abcd[341] = 2.0E0*I_ERI_F2xy_Fy2z_S_S_b-1*I_ERI_F2xy_Py_S_S;
  abcd[342] = 2.0E0*I_ERI_F2xz_Fy2z_S_S_b-1*I_ERI_F2xz_Py_S_S;
  abcd[343] = 2.0E0*I_ERI_Fx2y_Fy2z_S_S_b-1*I_ERI_Fx2y_Py_S_S;
  abcd[344] = 2.0E0*I_ERI_Fxyz_Fy2z_S_S_b-1*I_ERI_Fxyz_Py_S_S;
  abcd[345] = 2.0E0*I_ERI_Fx2z_Fy2z_S_S_b-1*I_ERI_Fx2z_Py_S_S;
  abcd[346] = 2.0E0*I_ERI_F3y_Fy2z_S_S_b-1*I_ERI_F3y_Py_S_S;
  abcd[347] = 2.0E0*I_ERI_F2yz_Fy2z_S_S_b-1*I_ERI_F2yz_Py_S_S;
  abcd[348] = 2.0E0*I_ERI_Fy2z_Fy2z_S_S_b-1*I_ERI_Fy2z_Py_S_S;
  abcd[349] = 2.0E0*I_ERI_F3z_Fy2z_S_S_b-1*I_ERI_F3z_Py_S_S;
  abcd[350] = 2.0E0*I_ERI_F3x_F3z_S_S_b-2*I_ERI_F3x_Pz_S_S;
  abcd[351] = 2.0E0*I_ERI_F2xy_F3z_S_S_b-2*I_ERI_F2xy_Pz_S_S;
  abcd[352] = 2.0E0*I_ERI_F2xz_F3z_S_S_b-2*I_ERI_F2xz_Pz_S_S;
  abcd[353] = 2.0E0*I_ERI_Fx2y_F3z_S_S_b-2*I_ERI_Fx2y_Pz_S_S;
  abcd[354] = 2.0E0*I_ERI_Fxyz_F3z_S_S_b-2*I_ERI_Fxyz_Pz_S_S;
  abcd[355] = 2.0E0*I_ERI_Fx2z_F3z_S_S_b-2*I_ERI_Fx2z_Pz_S_S;
  abcd[356] = 2.0E0*I_ERI_F3y_F3z_S_S_b-2*I_ERI_F3y_Pz_S_S;
  abcd[357] = 2.0E0*I_ERI_F2yz_F3z_S_S_b-2*I_ERI_F2yz_Pz_S_S;
  abcd[358] = 2.0E0*I_ERI_Fy2z_F3z_S_S_b-2*I_ERI_Fy2z_Pz_S_S;
  abcd[359] = 2.0E0*I_ERI_F3z_F3z_S_S_b-2*I_ERI_F3z_Pz_S_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dcx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_c
   ************************************************************/
  abcd[360] = 2.0E0*I_ERI_F3x_D2x_Px_S_c;
  abcd[361] = 2.0E0*I_ERI_F2xy_D2x_Px_S_c;
  abcd[362] = 2.0E0*I_ERI_F2xz_D2x_Px_S_c;
  abcd[363] = 2.0E0*I_ERI_Fx2y_D2x_Px_S_c;
  abcd[364] = 2.0E0*I_ERI_Fxyz_D2x_Px_S_c;
  abcd[365] = 2.0E0*I_ERI_Fx2z_D2x_Px_S_c;
  abcd[366] = 2.0E0*I_ERI_F3y_D2x_Px_S_c;
  abcd[367] = 2.0E0*I_ERI_F2yz_D2x_Px_S_c;
  abcd[368] = 2.0E0*I_ERI_Fy2z_D2x_Px_S_c;
  abcd[369] = 2.0E0*I_ERI_F3z_D2x_Px_S_c;
  abcd[370] = 2.0E0*I_ERI_F3x_Dxy_Px_S_c;
  abcd[371] = 2.0E0*I_ERI_F2xy_Dxy_Px_S_c;
  abcd[372] = 2.0E0*I_ERI_F2xz_Dxy_Px_S_c;
  abcd[373] = 2.0E0*I_ERI_Fx2y_Dxy_Px_S_c;
  abcd[374] = 2.0E0*I_ERI_Fxyz_Dxy_Px_S_c;
  abcd[375] = 2.0E0*I_ERI_Fx2z_Dxy_Px_S_c;
  abcd[376] = 2.0E0*I_ERI_F3y_Dxy_Px_S_c;
  abcd[377] = 2.0E0*I_ERI_F2yz_Dxy_Px_S_c;
  abcd[378] = 2.0E0*I_ERI_Fy2z_Dxy_Px_S_c;
  abcd[379] = 2.0E0*I_ERI_F3z_Dxy_Px_S_c;
  abcd[380] = 2.0E0*I_ERI_F3x_Dxz_Px_S_c;
  abcd[381] = 2.0E0*I_ERI_F2xy_Dxz_Px_S_c;
  abcd[382] = 2.0E0*I_ERI_F2xz_Dxz_Px_S_c;
  abcd[383] = 2.0E0*I_ERI_Fx2y_Dxz_Px_S_c;
  abcd[384] = 2.0E0*I_ERI_Fxyz_Dxz_Px_S_c;
  abcd[385] = 2.0E0*I_ERI_Fx2z_Dxz_Px_S_c;
  abcd[386] = 2.0E0*I_ERI_F3y_Dxz_Px_S_c;
  abcd[387] = 2.0E0*I_ERI_F2yz_Dxz_Px_S_c;
  abcd[388] = 2.0E0*I_ERI_Fy2z_Dxz_Px_S_c;
  abcd[389] = 2.0E0*I_ERI_F3z_Dxz_Px_S_c;
  abcd[390] = 2.0E0*I_ERI_F3x_D2y_Px_S_c;
  abcd[391] = 2.0E0*I_ERI_F2xy_D2y_Px_S_c;
  abcd[392] = 2.0E0*I_ERI_F2xz_D2y_Px_S_c;
  abcd[393] = 2.0E0*I_ERI_Fx2y_D2y_Px_S_c;
  abcd[394] = 2.0E0*I_ERI_Fxyz_D2y_Px_S_c;
  abcd[395] = 2.0E0*I_ERI_Fx2z_D2y_Px_S_c;
  abcd[396] = 2.0E0*I_ERI_F3y_D2y_Px_S_c;
  abcd[397] = 2.0E0*I_ERI_F2yz_D2y_Px_S_c;
  abcd[398] = 2.0E0*I_ERI_Fy2z_D2y_Px_S_c;
  abcd[399] = 2.0E0*I_ERI_F3z_D2y_Px_S_c;
  abcd[400] = 2.0E0*I_ERI_F3x_Dyz_Px_S_c;
  abcd[401] = 2.0E0*I_ERI_F2xy_Dyz_Px_S_c;
  abcd[402] = 2.0E0*I_ERI_F2xz_Dyz_Px_S_c;
  abcd[403] = 2.0E0*I_ERI_Fx2y_Dyz_Px_S_c;
  abcd[404] = 2.0E0*I_ERI_Fxyz_Dyz_Px_S_c;
  abcd[405] = 2.0E0*I_ERI_Fx2z_Dyz_Px_S_c;
  abcd[406] = 2.0E0*I_ERI_F3y_Dyz_Px_S_c;
  abcd[407] = 2.0E0*I_ERI_F2yz_Dyz_Px_S_c;
  abcd[408] = 2.0E0*I_ERI_Fy2z_Dyz_Px_S_c;
  abcd[409] = 2.0E0*I_ERI_F3z_Dyz_Px_S_c;
  abcd[410] = 2.0E0*I_ERI_F3x_D2z_Px_S_c;
  abcd[411] = 2.0E0*I_ERI_F2xy_D2z_Px_S_c;
  abcd[412] = 2.0E0*I_ERI_F2xz_D2z_Px_S_c;
  abcd[413] = 2.0E0*I_ERI_Fx2y_D2z_Px_S_c;
  abcd[414] = 2.0E0*I_ERI_Fxyz_D2z_Px_S_c;
  abcd[415] = 2.0E0*I_ERI_Fx2z_D2z_Px_S_c;
  abcd[416] = 2.0E0*I_ERI_F3y_D2z_Px_S_c;
  abcd[417] = 2.0E0*I_ERI_F2yz_D2z_Px_S_c;
  abcd[418] = 2.0E0*I_ERI_Fy2z_D2z_Px_S_c;
  abcd[419] = 2.0E0*I_ERI_F3z_D2z_Px_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dcy
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_c
   ************************************************************/
  abcd[420] = 2.0E0*I_ERI_F3x_D2x_Py_S_c;
  abcd[421] = 2.0E0*I_ERI_F2xy_D2x_Py_S_c;
  abcd[422] = 2.0E0*I_ERI_F2xz_D2x_Py_S_c;
  abcd[423] = 2.0E0*I_ERI_Fx2y_D2x_Py_S_c;
  abcd[424] = 2.0E0*I_ERI_Fxyz_D2x_Py_S_c;
  abcd[425] = 2.0E0*I_ERI_Fx2z_D2x_Py_S_c;
  abcd[426] = 2.0E0*I_ERI_F3y_D2x_Py_S_c;
  abcd[427] = 2.0E0*I_ERI_F2yz_D2x_Py_S_c;
  abcd[428] = 2.0E0*I_ERI_Fy2z_D2x_Py_S_c;
  abcd[429] = 2.0E0*I_ERI_F3z_D2x_Py_S_c;
  abcd[430] = 2.0E0*I_ERI_F3x_Dxy_Py_S_c;
  abcd[431] = 2.0E0*I_ERI_F2xy_Dxy_Py_S_c;
  abcd[432] = 2.0E0*I_ERI_F2xz_Dxy_Py_S_c;
  abcd[433] = 2.0E0*I_ERI_Fx2y_Dxy_Py_S_c;
  abcd[434] = 2.0E0*I_ERI_Fxyz_Dxy_Py_S_c;
  abcd[435] = 2.0E0*I_ERI_Fx2z_Dxy_Py_S_c;
  abcd[436] = 2.0E0*I_ERI_F3y_Dxy_Py_S_c;
  abcd[437] = 2.0E0*I_ERI_F2yz_Dxy_Py_S_c;
  abcd[438] = 2.0E0*I_ERI_Fy2z_Dxy_Py_S_c;
  abcd[439] = 2.0E0*I_ERI_F3z_Dxy_Py_S_c;
  abcd[440] = 2.0E0*I_ERI_F3x_Dxz_Py_S_c;
  abcd[441] = 2.0E0*I_ERI_F2xy_Dxz_Py_S_c;
  abcd[442] = 2.0E0*I_ERI_F2xz_Dxz_Py_S_c;
  abcd[443] = 2.0E0*I_ERI_Fx2y_Dxz_Py_S_c;
  abcd[444] = 2.0E0*I_ERI_Fxyz_Dxz_Py_S_c;
  abcd[445] = 2.0E0*I_ERI_Fx2z_Dxz_Py_S_c;
  abcd[446] = 2.0E0*I_ERI_F3y_Dxz_Py_S_c;
  abcd[447] = 2.0E0*I_ERI_F2yz_Dxz_Py_S_c;
  abcd[448] = 2.0E0*I_ERI_Fy2z_Dxz_Py_S_c;
  abcd[449] = 2.0E0*I_ERI_F3z_Dxz_Py_S_c;
  abcd[450] = 2.0E0*I_ERI_F3x_D2y_Py_S_c;
  abcd[451] = 2.0E0*I_ERI_F2xy_D2y_Py_S_c;
  abcd[452] = 2.0E0*I_ERI_F2xz_D2y_Py_S_c;
  abcd[453] = 2.0E0*I_ERI_Fx2y_D2y_Py_S_c;
  abcd[454] = 2.0E0*I_ERI_Fxyz_D2y_Py_S_c;
  abcd[455] = 2.0E0*I_ERI_Fx2z_D2y_Py_S_c;
  abcd[456] = 2.0E0*I_ERI_F3y_D2y_Py_S_c;
  abcd[457] = 2.0E0*I_ERI_F2yz_D2y_Py_S_c;
  abcd[458] = 2.0E0*I_ERI_Fy2z_D2y_Py_S_c;
  abcd[459] = 2.0E0*I_ERI_F3z_D2y_Py_S_c;
  abcd[460] = 2.0E0*I_ERI_F3x_Dyz_Py_S_c;
  abcd[461] = 2.0E0*I_ERI_F2xy_Dyz_Py_S_c;
  abcd[462] = 2.0E0*I_ERI_F2xz_Dyz_Py_S_c;
  abcd[463] = 2.0E0*I_ERI_Fx2y_Dyz_Py_S_c;
  abcd[464] = 2.0E0*I_ERI_Fxyz_Dyz_Py_S_c;
  abcd[465] = 2.0E0*I_ERI_Fx2z_Dyz_Py_S_c;
  abcd[466] = 2.0E0*I_ERI_F3y_Dyz_Py_S_c;
  abcd[467] = 2.0E0*I_ERI_F2yz_Dyz_Py_S_c;
  abcd[468] = 2.0E0*I_ERI_Fy2z_Dyz_Py_S_c;
  abcd[469] = 2.0E0*I_ERI_F3z_Dyz_Py_S_c;
  abcd[470] = 2.0E0*I_ERI_F3x_D2z_Py_S_c;
  abcd[471] = 2.0E0*I_ERI_F2xy_D2z_Py_S_c;
  abcd[472] = 2.0E0*I_ERI_F2xz_D2z_Py_S_c;
  abcd[473] = 2.0E0*I_ERI_Fx2y_D2z_Py_S_c;
  abcd[474] = 2.0E0*I_ERI_Fxyz_D2z_Py_S_c;
  abcd[475] = 2.0E0*I_ERI_Fx2z_D2z_Py_S_c;
  abcd[476] = 2.0E0*I_ERI_F3y_D2z_Py_S_c;
  abcd[477] = 2.0E0*I_ERI_F2yz_D2z_Py_S_c;
  abcd[478] = 2.0E0*I_ERI_Fy2z_D2z_Py_S_c;
  abcd[479] = 2.0E0*I_ERI_F3z_D2z_Py_S_c;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_S_S_dcz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_D_P_S_c
   ************************************************************/
  abcd[480] = 2.0E0*I_ERI_F3x_D2x_Pz_S_c;
  abcd[481] = 2.0E0*I_ERI_F2xy_D2x_Pz_S_c;
  abcd[482] = 2.0E0*I_ERI_F2xz_D2x_Pz_S_c;
  abcd[483] = 2.0E0*I_ERI_Fx2y_D2x_Pz_S_c;
  abcd[484] = 2.0E0*I_ERI_Fxyz_D2x_Pz_S_c;
  abcd[485] = 2.0E0*I_ERI_Fx2z_D2x_Pz_S_c;
  abcd[486] = 2.0E0*I_ERI_F3y_D2x_Pz_S_c;
  abcd[487] = 2.0E0*I_ERI_F2yz_D2x_Pz_S_c;
  abcd[488] = 2.0E0*I_ERI_Fy2z_D2x_Pz_S_c;
  abcd[489] = 2.0E0*I_ERI_F3z_D2x_Pz_S_c;
  abcd[490] = 2.0E0*I_ERI_F3x_Dxy_Pz_S_c;
  abcd[491] = 2.0E0*I_ERI_F2xy_Dxy_Pz_S_c;
  abcd[492] = 2.0E0*I_ERI_F2xz_Dxy_Pz_S_c;
  abcd[493] = 2.0E0*I_ERI_Fx2y_Dxy_Pz_S_c;
  abcd[494] = 2.0E0*I_ERI_Fxyz_Dxy_Pz_S_c;
  abcd[495] = 2.0E0*I_ERI_Fx2z_Dxy_Pz_S_c;
  abcd[496] = 2.0E0*I_ERI_F3y_Dxy_Pz_S_c;
  abcd[497] = 2.0E0*I_ERI_F2yz_Dxy_Pz_S_c;
  abcd[498] = 2.0E0*I_ERI_Fy2z_Dxy_Pz_S_c;
  abcd[499] = 2.0E0*I_ERI_F3z_Dxy_Pz_S_c;
  abcd[500] = 2.0E0*I_ERI_F3x_Dxz_Pz_S_c;
  abcd[501] = 2.0E0*I_ERI_F2xy_Dxz_Pz_S_c;
  abcd[502] = 2.0E0*I_ERI_F2xz_Dxz_Pz_S_c;
  abcd[503] = 2.0E0*I_ERI_Fx2y_Dxz_Pz_S_c;
  abcd[504] = 2.0E0*I_ERI_Fxyz_Dxz_Pz_S_c;
  abcd[505] = 2.0E0*I_ERI_Fx2z_Dxz_Pz_S_c;
  abcd[506] = 2.0E0*I_ERI_F3y_Dxz_Pz_S_c;
  abcd[507] = 2.0E0*I_ERI_F2yz_Dxz_Pz_S_c;
  abcd[508] = 2.0E0*I_ERI_Fy2z_Dxz_Pz_S_c;
  abcd[509] = 2.0E0*I_ERI_F3z_Dxz_Pz_S_c;
  abcd[510] = 2.0E0*I_ERI_F3x_D2y_Pz_S_c;
  abcd[511] = 2.0E0*I_ERI_F2xy_D2y_Pz_S_c;
  abcd[512] = 2.0E0*I_ERI_F2xz_D2y_Pz_S_c;
  abcd[513] = 2.0E0*I_ERI_Fx2y_D2y_Pz_S_c;
  abcd[514] = 2.0E0*I_ERI_Fxyz_D2y_Pz_S_c;
  abcd[515] = 2.0E0*I_ERI_Fx2z_D2y_Pz_S_c;
  abcd[516] = 2.0E0*I_ERI_F3y_D2y_Pz_S_c;
  abcd[517] = 2.0E0*I_ERI_F2yz_D2y_Pz_S_c;
  abcd[518] = 2.0E0*I_ERI_Fy2z_D2y_Pz_S_c;
  abcd[519] = 2.0E0*I_ERI_F3z_D2y_Pz_S_c;
  abcd[520] = 2.0E0*I_ERI_F3x_Dyz_Pz_S_c;
  abcd[521] = 2.0E0*I_ERI_F2xy_Dyz_Pz_S_c;
  abcd[522] = 2.0E0*I_ERI_F2xz_Dyz_Pz_S_c;
  abcd[523] = 2.0E0*I_ERI_Fx2y_Dyz_Pz_S_c;
  abcd[524] = 2.0E0*I_ERI_Fxyz_Dyz_Pz_S_c;
  abcd[525] = 2.0E0*I_ERI_Fx2z_Dyz_Pz_S_c;
  abcd[526] = 2.0E0*I_ERI_F3y_Dyz_Pz_S_c;
  abcd[527] = 2.0E0*I_ERI_F2yz_Dyz_Pz_S_c;
  abcd[528] = 2.0E0*I_ERI_Fy2z_Dyz_Pz_S_c;
  abcd[529] = 2.0E0*I_ERI_F3z_Dyz_Pz_S_c;
  abcd[530] = 2.0E0*I_ERI_F3x_D2z_Pz_S_c;
  abcd[531] = 2.0E0*I_ERI_F2xy_D2z_Pz_S_c;
  abcd[532] = 2.0E0*I_ERI_F2xz_D2z_Pz_S_c;
  abcd[533] = 2.0E0*I_ERI_Fx2y_D2z_Pz_S_c;
  abcd[534] = 2.0E0*I_ERI_Fxyz_D2z_Pz_S_c;
  abcd[535] = 2.0E0*I_ERI_Fx2z_D2z_Pz_S_c;
  abcd[536] = 2.0E0*I_ERI_F3y_D2z_Pz_S_c;
  abcd[537] = 2.0E0*I_ERI_F2yz_D2z_Pz_S_c;
  abcd[538] = 2.0E0*I_ERI_Fy2z_D2z_Pz_S_c;
  abcd[539] = 2.0E0*I_ERI_F3z_D2z_Pz_S_c;
}
