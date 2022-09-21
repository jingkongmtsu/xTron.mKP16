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

void hgp_os_eri_g_d_g_s(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_I6x_S_G4x_S = 0.0E0;
  Double I_ERI_I5xy_S_G4x_S = 0.0E0;
  Double I_ERI_I5xz_S_G4x_S = 0.0E0;
  Double I_ERI_I4x2y_S_G4x_S = 0.0E0;
  Double I_ERI_I4xyz_S_G4x_S = 0.0E0;
  Double I_ERI_I4x2z_S_G4x_S = 0.0E0;
  Double I_ERI_I3x3y_S_G4x_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G4x_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G4x_S = 0.0E0;
  Double I_ERI_I3x3z_S_G4x_S = 0.0E0;
  Double I_ERI_I2x4y_S_G4x_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G4x_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G4x_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G4x_S = 0.0E0;
  Double I_ERI_I2x4z_S_G4x_S = 0.0E0;
  Double I_ERI_Ix5y_S_G4x_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G4x_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G4x_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G4x_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G4x_S = 0.0E0;
  Double I_ERI_Ix5z_S_G4x_S = 0.0E0;
  Double I_ERI_I6y_S_G4x_S = 0.0E0;
  Double I_ERI_I5yz_S_G4x_S = 0.0E0;
  Double I_ERI_I4y2z_S_G4x_S = 0.0E0;
  Double I_ERI_I3y3z_S_G4x_S = 0.0E0;
  Double I_ERI_I2y4z_S_G4x_S = 0.0E0;
  Double I_ERI_Iy5z_S_G4x_S = 0.0E0;
  Double I_ERI_I6z_S_G4x_S = 0.0E0;
  Double I_ERI_I6x_S_G3xy_S = 0.0E0;
  Double I_ERI_I5xy_S_G3xy_S = 0.0E0;
  Double I_ERI_I5xz_S_G3xy_S = 0.0E0;
  Double I_ERI_I4x2y_S_G3xy_S = 0.0E0;
  Double I_ERI_I4xyz_S_G3xy_S = 0.0E0;
  Double I_ERI_I4x2z_S_G3xy_S = 0.0E0;
  Double I_ERI_I3x3y_S_G3xy_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G3xy_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G3xy_S = 0.0E0;
  Double I_ERI_I3x3z_S_G3xy_S = 0.0E0;
  Double I_ERI_I2x4y_S_G3xy_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G3xy_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G3xy_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G3xy_S = 0.0E0;
  Double I_ERI_I2x4z_S_G3xy_S = 0.0E0;
  Double I_ERI_Ix5y_S_G3xy_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G3xy_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G3xy_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G3xy_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G3xy_S = 0.0E0;
  Double I_ERI_Ix5z_S_G3xy_S = 0.0E0;
  Double I_ERI_I6y_S_G3xy_S = 0.0E0;
  Double I_ERI_I5yz_S_G3xy_S = 0.0E0;
  Double I_ERI_I4y2z_S_G3xy_S = 0.0E0;
  Double I_ERI_I3y3z_S_G3xy_S = 0.0E0;
  Double I_ERI_I2y4z_S_G3xy_S = 0.0E0;
  Double I_ERI_Iy5z_S_G3xy_S = 0.0E0;
  Double I_ERI_I6z_S_G3xy_S = 0.0E0;
  Double I_ERI_I6x_S_G3xz_S = 0.0E0;
  Double I_ERI_I5xy_S_G3xz_S = 0.0E0;
  Double I_ERI_I5xz_S_G3xz_S = 0.0E0;
  Double I_ERI_I4x2y_S_G3xz_S = 0.0E0;
  Double I_ERI_I4xyz_S_G3xz_S = 0.0E0;
  Double I_ERI_I4x2z_S_G3xz_S = 0.0E0;
  Double I_ERI_I3x3y_S_G3xz_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G3xz_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G3xz_S = 0.0E0;
  Double I_ERI_I3x3z_S_G3xz_S = 0.0E0;
  Double I_ERI_I2x4y_S_G3xz_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G3xz_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G3xz_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G3xz_S = 0.0E0;
  Double I_ERI_I2x4z_S_G3xz_S = 0.0E0;
  Double I_ERI_Ix5y_S_G3xz_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G3xz_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G3xz_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G3xz_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G3xz_S = 0.0E0;
  Double I_ERI_Ix5z_S_G3xz_S = 0.0E0;
  Double I_ERI_I6y_S_G3xz_S = 0.0E0;
  Double I_ERI_I5yz_S_G3xz_S = 0.0E0;
  Double I_ERI_I4y2z_S_G3xz_S = 0.0E0;
  Double I_ERI_I3y3z_S_G3xz_S = 0.0E0;
  Double I_ERI_I2y4z_S_G3xz_S = 0.0E0;
  Double I_ERI_Iy5z_S_G3xz_S = 0.0E0;
  Double I_ERI_I6z_S_G3xz_S = 0.0E0;
  Double I_ERI_I6x_S_G2x2y_S = 0.0E0;
  Double I_ERI_I5xy_S_G2x2y_S = 0.0E0;
  Double I_ERI_I5xz_S_G2x2y_S = 0.0E0;
  Double I_ERI_I4x2y_S_G2x2y_S = 0.0E0;
  Double I_ERI_I4xyz_S_G2x2y_S = 0.0E0;
  Double I_ERI_I4x2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I3x3y_S_G2x2y_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I3x3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I2x4y_S_G2x2y_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I2x4z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Ix5y_S_G2x2y_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Ix5z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I6y_S_G2x2y_S = 0.0E0;
  Double I_ERI_I5yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_I4y2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I3y3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I2y4z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Iy5z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I6z_S_G2x2y_S = 0.0E0;
  Double I_ERI_I6x_S_G2xyz_S = 0.0E0;
  Double I_ERI_I5xy_S_G2xyz_S = 0.0E0;
  Double I_ERI_I5xz_S_G2xyz_S = 0.0E0;
  Double I_ERI_I4x2y_S_G2xyz_S = 0.0E0;
  Double I_ERI_I4xyz_S_G2xyz_S = 0.0E0;
  Double I_ERI_I4x2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I3x3y_S_G2xyz_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I3x3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I2x4y_S_G2xyz_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I2x4z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Ix5y_S_G2xyz_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Ix5z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I6y_S_G2xyz_S = 0.0E0;
  Double I_ERI_I5yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_I4y2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I3y3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I2y4z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Iy5z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I6z_S_G2xyz_S = 0.0E0;
  Double I_ERI_I6x_S_G2x2z_S = 0.0E0;
  Double I_ERI_I5xy_S_G2x2z_S = 0.0E0;
  Double I_ERI_I5xz_S_G2x2z_S = 0.0E0;
  Double I_ERI_I4x2y_S_G2x2z_S = 0.0E0;
  Double I_ERI_I4xyz_S_G2x2z_S = 0.0E0;
  Double I_ERI_I4x2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I3x3y_S_G2x2z_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I3x3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I2x4y_S_G2x2z_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I2x4z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Ix5y_S_G2x2z_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Ix5z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I6y_S_G2x2z_S = 0.0E0;
  Double I_ERI_I5yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_I4y2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I3y3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I2y4z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Iy5z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I6z_S_G2x2z_S = 0.0E0;
  Double I_ERI_I6x_S_Gx3y_S = 0.0E0;
  Double I_ERI_I5xy_S_Gx3y_S = 0.0E0;
  Double I_ERI_I5xz_S_Gx3y_S = 0.0E0;
  Double I_ERI_I4x2y_S_Gx3y_S = 0.0E0;
  Double I_ERI_I4xyz_S_Gx3y_S = 0.0E0;
  Double I_ERI_I4x2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I3x3y_S_Gx3y_S = 0.0E0;
  Double I_ERI_I3x2yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_I3xy2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I3x3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I2x4y_S_Gx3y_S = 0.0E0;
  Double I_ERI_I2x3yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I2xy3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I2x4z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Ix5y_S_Gx3y_S = 0.0E0;
  Double I_ERI_Ix4yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Ixy4z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Ix5z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I6y_S_Gx3y_S = 0.0E0;
  Double I_ERI_I5yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_I4y2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I3y3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I2y4z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Iy5z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I6z_S_Gx3y_S = 0.0E0;
  Double I_ERI_I6x_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I5xy_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I5xz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I4x2y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I4xyz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I4x2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I3x3y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I3x2yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I3xy2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I3x3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I2x4y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I2x3yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I2xy3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I2x4z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Ix5y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Ix4yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Ixy4z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Ix5z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I6y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I5yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I4y2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I3y3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I2y4z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Iy5z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I6z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_I6x_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I5xy_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I5xz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I4x2y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I4xyz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I4x2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I3x3y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I3x2yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I3xy2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I3x3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I2x4y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I2x3yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I2xy3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I2x4z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Ix5y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Ix4yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Ixy4z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Ix5z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I6y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I5yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I4y2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I3y3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I2y4z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Iy5z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I6z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_I6x_S_Gx3z_S = 0.0E0;
  Double I_ERI_I5xy_S_Gx3z_S = 0.0E0;
  Double I_ERI_I5xz_S_Gx3z_S = 0.0E0;
  Double I_ERI_I4x2y_S_Gx3z_S = 0.0E0;
  Double I_ERI_I4xyz_S_Gx3z_S = 0.0E0;
  Double I_ERI_I4x2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I3x3y_S_Gx3z_S = 0.0E0;
  Double I_ERI_I3x2yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_I3xy2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I3x3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I2x4y_S_Gx3z_S = 0.0E0;
  Double I_ERI_I2x3yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I2xy3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I2x4z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Ix5y_S_Gx3z_S = 0.0E0;
  Double I_ERI_Ix4yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Ixy4z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Ix5z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I6y_S_Gx3z_S = 0.0E0;
  Double I_ERI_I5yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_I4y2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I3y3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I2y4z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Iy5z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I6z_S_Gx3z_S = 0.0E0;
  Double I_ERI_I6x_S_G4y_S = 0.0E0;
  Double I_ERI_I5xy_S_G4y_S = 0.0E0;
  Double I_ERI_I5xz_S_G4y_S = 0.0E0;
  Double I_ERI_I4x2y_S_G4y_S = 0.0E0;
  Double I_ERI_I4xyz_S_G4y_S = 0.0E0;
  Double I_ERI_I4x2z_S_G4y_S = 0.0E0;
  Double I_ERI_I3x3y_S_G4y_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G4y_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G4y_S = 0.0E0;
  Double I_ERI_I3x3z_S_G4y_S = 0.0E0;
  Double I_ERI_I2x4y_S_G4y_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G4y_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G4y_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G4y_S = 0.0E0;
  Double I_ERI_I2x4z_S_G4y_S = 0.0E0;
  Double I_ERI_Ix5y_S_G4y_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G4y_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G4y_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G4y_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G4y_S = 0.0E0;
  Double I_ERI_Ix5z_S_G4y_S = 0.0E0;
  Double I_ERI_I6y_S_G4y_S = 0.0E0;
  Double I_ERI_I5yz_S_G4y_S = 0.0E0;
  Double I_ERI_I4y2z_S_G4y_S = 0.0E0;
  Double I_ERI_I3y3z_S_G4y_S = 0.0E0;
  Double I_ERI_I2y4z_S_G4y_S = 0.0E0;
  Double I_ERI_Iy5z_S_G4y_S = 0.0E0;
  Double I_ERI_I6z_S_G4y_S = 0.0E0;
  Double I_ERI_I6x_S_G3yz_S = 0.0E0;
  Double I_ERI_I5xy_S_G3yz_S = 0.0E0;
  Double I_ERI_I5xz_S_G3yz_S = 0.0E0;
  Double I_ERI_I4x2y_S_G3yz_S = 0.0E0;
  Double I_ERI_I4xyz_S_G3yz_S = 0.0E0;
  Double I_ERI_I4x2z_S_G3yz_S = 0.0E0;
  Double I_ERI_I3x3y_S_G3yz_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G3yz_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G3yz_S = 0.0E0;
  Double I_ERI_I3x3z_S_G3yz_S = 0.0E0;
  Double I_ERI_I2x4y_S_G3yz_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G3yz_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G3yz_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G3yz_S = 0.0E0;
  Double I_ERI_I2x4z_S_G3yz_S = 0.0E0;
  Double I_ERI_Ix5y_S_G3yz_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G3yz_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G3yz_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G3yz_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G3yz_S = 0.0E0;
  Double I_ERI_Ix5z_S_G3yz_S = 0.0E0;
  Double I_ERI_I6y_S_G3yz_S = 0.0E0;
  Double I_ERI_I5yz_S_G3yz_S = 0.0E0;
  Double I_ERI_I4y2z_S_G3yz_S = 0.0E0;
  Double I_ERI_I3y3z_S_G3yz_S = 0.0E0;
  Double I_ERI_I2y4z_S_G3yz_S = 0.0E0;
  Double I_ERI_Iy5z_S_G3yz_S = 0.0E0;
  Double I_ERI_I6z_S_G3yz_S = 0.0E0;
  Double I_ERI_I6x_S_G2y2z_S = 0.0E0;
  Double I_ERI_I5xy_S_G2y2z_S = 0.0E0;
  Double I_ERI_I5xz_S_G2y2z_S = 0.0E0;
  Double I_ERI_I4x2y_S_G2y2z_S = 0.0E0;
  Double I_ERI_I4xyz_S_G2y2z_S = 0.0E0;
  Double I_ERI_I4x2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I3x3y_S_G2y2z_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I3x3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I2x4y_S_G2y2z_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I2x4z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Ix5y_S_G2y2z_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Ix5z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I6y_S_G2y2z_S = 0.0E0;
  Double I_ERI_I5yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_I4y2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I3y3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I2y4z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Iy5z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I6z_S_G2y2z_S = 0.0E0;
  Double I_ERI_I6x_S_Gy3z_S = 0.0E0;
  Double I_ERI_I5xy_S_Gy3z_S = 0.0E0;
  Double I_ERI_I5xz_S_Gy3z_S = 0.0E0;
  Double I_ERI_I4x2y_S_Gy3z_S = 0.0E0;
  Double I_ERI_I4xyz_S_Gy3z_S = 0.0E0;
  Double I_ERI_I4x2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I3x3y_S_Gy3z_S = 0.0E0;
  Double I_ERI_I3x2yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_I3xy2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I3x3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I2x4y_S_Gy3z_S = 0.0E0;
  Double I_ERI_I2x3yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I2xy3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I2x4z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Ix5y_S_Gy3z_S = 0.0E0;
  Double I_ERI_Ix4yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Ixy4z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Ix5z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I6y_S_Gy3z_S = 0.0E0;
  Double I_ERI_I5yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_I4y2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I3y3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I2y4z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Iy5z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I6z_S_Gy3z_S = 0.0E0;
  Double I_ERI_I6x_S_G4z_S = 0.0E0;
  Double I_ERI_I5xy_S_G4z_S = 0.0E0;
  Double I_ERI_I5xz_S_G4z_S = 0.0E0;
  Double I_ERI_I4x2y_S_G4z_S = 0.0E0;
  Double I_ERI_I4xyz_S_G4z_S = 0.0E0;
  Double I_ERI_I4x2z_S_G4z_S = 0.0E0;
  Double I_ERI_I3x3y_S_G4z_S = 0.0E0;
  Double I_ERI_I3x2yz_S_G4z_S = 0.0E0;
  Double I_ERI_I3xy2z_S_G4z_S = 0.0E0;
  Double I_ERI_I3x3z_S_G4z_S = 0.0E0;
  Double I_ERI_I2x4y_S_G4z_S = 0.0E0;
  Double I_ERI_I2x3yz_S_G4z_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_G4z_S = 0.0E0;
  Double I_ERI_I2xy3z_S_G4z_S = 0.0E0;
  Double I_ERI_I2x4z_S_G4z_S = 0.0E0;
  Double I_ERI_Ix5y_S_G4z_S = 0.0E0;
  Double I_ERI_Ix4yz_S_G4z_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_G4z_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_G4z_S = 0.0E0;
  Double I_ERI_Ixy4z_S_G4z_S = 0.0E0;
  Double I_ERI_Ix5z_S_G4z_S = 0.0E0;
  Double I_ERI_I6y_S_G4z_S = 0.0E0;
  Double I_ERI_I5yz_S_G4z_S = 0.0E0;
  Double I_ERI_I4y2z_S_G4z_S = 0.0E0;
  Double I_ERI_I3y3z_S_G4z_S = 0.0E0;
  Double I_ERI_I2y4z_S_G4z_S = 0.0E0;
  Double I_ERI_Iy5z_S_G4z_S = 0.0E0;
  Double I_ERI_I6z_S_G4z_S = 0.0E0;
  Double I_ERI_H5x_S_G4x_S = 0.0E0;
  Double I_ERI_H4xy_S_G4x_S = 0.0E0;
  Double I_ERI_H4xz_S_G4x_S = 0.0E0;
  Double I_ERI_H3x2y_S_G4x_S = 0.0E0;
  Double I_ERI_H3xyz_S_G4x_S = 0.0E0;
  Double I_ERI_H3x2z_S_G4x_S = 0.0E0;
  Double I_ERI_H2x3y_S_G4x_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G4x_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G4x_S = 0.0E0;
  Double I_ERI_H2x3z_S_G4x_S = 0.0E0;
  Double I_ERI_Hx4y_S_G4x_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G4x_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G4x_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G4x_S = 0.0E0;
  Double I_ERI_Hx4z_S_G4x_S = 0.0E0;
  Double I_ERI_H5y_S_G4x_S = 0.0E0;
  Double I_ERI_H4yz_S_G4x_S = 0.0E0;
  Double I_ERI_H3y2z_S_G4x_S = 0.0E0;
  Double I_ERI_H2y3z_S_G4x_S = 0.0E0;
  Double I_ERI_Hy4z_S_G4x_S = 0.0E0;
  Double I_ERI_H5z_S_G4x_S = 0.0E0;
  Double I_ERI_H5x_S_G3xy_S = 0.0E0;
  Double I_ERI_H4xy_S_G3xy_S = 0.0E0;
  Double I_ERI_H4xz_S_G3xy_S = 0.0E0;
  Double I_ERI_H3x2y_S_G3xy_S = 0.0E0;
  Double I_ERI_H3xyz_S_G3xy_S = 0.0E0;
  Double I_ERI_H3x2z_S_G3xy_S = 0.0E0;
  Double I_ERI_H2x3y_S_G3xy_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G3xy_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G3xy_S = 0.0E0;
  Double I_ERI_H2x3z_S_G3xy_S = 0.0E0;
  Double I_ERI_Hx4y_S_G3xy_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G3xy_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G3xy_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G3xy_S = 0.0E0;
  Double I_ERI_Hx4z_S_G3xy_S = 0.0E0;
  Double I_ERI_H5y_S_G3xy_S = 0.0E0;
  Double I_ERI_H4yz_S_G3xy_S = 0.0E0;
  Double I_ERI_H3y2z_S_G3xy_S = 0.0E0;
  Double I_ERI_H2y3z_S_G3xy_S = 0.0E0;
  Double I_ERI_Hy4z_S_G3xy_S = 0.0E0;
  Double I_ERI_H5z_S_G3xy_S = 0.0E0;
  Double I_ERI_H5x_S_G3xz_S = 0.0E0;
  Double I_ERI_H4xy_S_G3xz_S = 0.0E0;
  Double I_ERI_H4xz_S_G3xz_S = 0.0E0;
  Double I_ERI_H3x2y_S_G3xz_S = 0.0E0;
  Double I_ERI_H3xyz_S_G3xz_S = 0.0E0;
  Double I_ERI_H3x2z_S_G3xz_S = 0.0E0;
  Double I_ERI_H2x3y_S_G3xz_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G3xz_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G3xz_S = 0.0E0;
  Double I_ERI_H2x3z_S_G3xz_S = 0.0E0;
  Double I_ERI_Hx4y_S_G3xz_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G3xz_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G3xz_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G3xz_S = 0.0E0;
  Double I_ERI_Hx4z_S_G3xz_S = 0.0E0;
  Double I_ERI_H5y_S_G3xz_S = 0.0E0;
  Double I_ERI_H4yz_S_G3xz_S = 0.0E0;
  Double I_ERI_H3y2z_S_G3xz_S = 0.0E0;
  Double I_ERI_H2y3z_S_G3xz_S = 0.0E0;
  Double I_ERI_Hy4z_S_G3xz_S = 0.0E0;
  Double I_ERI_H5z_S_G3xz_S = 0.0E0;
  Double I_ERI_H5x_S_G2x2y_S = 0.0E0;
  Double I_ERI_H4xy_S_G2x2y_S = 0.0E0;
  Double I_ERI_H4xz_S_G2x2y_S = 0.0E0;
  Double I_ERI_H3x2y_S_G2x2y_S = 0.0E0;
  Double I_ERI_H3xyz_S_G2x2y_S = 0.0E0;
  Double I_ERI_H3x2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_H2x3y_S_G2x2y_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_H2x3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Hx4y_S_G2x2y_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Hx4z_S_G2x2y_S = 0.0E0;
  Double I_ERI_H5y_S_G2x2y_S = 0.0E0;
  Double I_ERI_H4yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_H3y2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_H2y3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Hy4z_S_G2x2y_S = 0.0E0;
  Double I_ERI_H5z_S_G2x2y_S = 0.0E0;
  Double I_ERI_H5x_S_G2xyz_S = 0.0E0;
  Double I_ERI_H4xy_S_G2xyz_S = 0.0E0;
  Double I_ERI_H4xz_S_G2xyz_S = 0.0E0;
  Double I_ERI_H3x2y_S_G2xyz_S = 0.0E0;
  Double I_ERI_H3xyz_S_G2xyz_S = 0.0E0;
  Double I_ERI_H3x2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_H2x3y_S_G2xyz_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_H2x3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Hx4y_S_G2xyz_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Hx4z_S_G2xyz_S = 0.0E0;
  Double I_ERI_H5y_S_G2xyz_S = 0.0E0;
  Double I_ERI_H4yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_H3y2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_H2y3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Hy4z_S_G2xyz_S = 0.0E0;
  Double I_ERI_H5z_S_G2xyz_S = 0.0E0;
  Double I_ERI_H5x_S_G2x2z_S = 0.0E0;
  Double I_ERI_H4xy_S_G2x2z_S = 0.0E0;
  Double I_ERI_H4xz_S_G2x2z_S = 0.0E0;
  Double I_ERI_H3x2y_S_G2x2z_S = 0.0E0;
  Double I_ERI_H3xyz_S_G2x2z_S = 0.0E0;
  Double I_ERI_H3x2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_H2x3y_S_G2x2z_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_H2x3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Hx4y_S_G2x2z_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Hx4z_S_G2x2z_S = 0.0E0;
  Double I_ERI_H5y_S_G2x2z_S = 0.0E0;
  Double I_ERI_H4yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_H3y2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_H2y3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Hy4z_S_G2x2z_S = 0.0E0;
  Double I_ERI_H5z_S_G2x2z_S = 0.0E0;
  Double I_ERI_H5x_S_Gx3y_S = 0.0E0;
  Double I_ERI_H4xy_S_Gx3y_S = 0.0E0;
  Double I_ERI_H4xz_S_Gx3y_S = 0.0E0;
  Double I_ERI_H3x2y_S_Gx3y_S = 0.0E0;
  Double I_ERI_H3xyz_S_Gx3y_S = 0.0E0;
  Double I_ERI_H3x2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_H2x3y_S_Gx3y_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_H2x3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Hx4y_S_Gx3y_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Hx4z_S_Gx3y_S = 0.0E0;
  Double I_ERI_H5y_S_Gx3y_S = 0.0E0;
  Double I_ERI_H4yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_H3y2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_H2y3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Hy4z_S_Gx3y_S = 0.0E0;
  Double I_ERI_H5z_S_Gx3y_S = 0.0E0;
  Double I_ERI_H5x_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H4xy_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H4xz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H3x2y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H3xyz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H3x2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H2x3y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H2x3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Hx4y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Hx4z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H5y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H4yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H3y2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H2y3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Hy4z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H5z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_H5x_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H4xy_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H4xz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H3x2y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H3xyz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H3x2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H2x3y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H2x3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Hx4y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Hx4z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H5y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H4yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H3y2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H2y3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Hy4z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H5z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_H5x_S_Gx3z_S = 0.0E0;
  Double I_ERI_H4xy_S_Gx3z_S = 0.0E0;
  Double I_ERI_H4xz_S_Gx3z_S = 0.0E0;
  Double I_ERI_H3x2y_S_Gx3z_S = 0.0E0;
  Double I_ERI_H3xyz_S_Gx3z_S = 0.0E0;
  Double I_ERI_H3x2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_H2x3y_S_Gx3z_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_H2x3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Hx4y_S_Gx3z_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Hx4z_S_Gx3z_S = 0.0E0;
  Double I_ERI_H5y_S_Gx3z_S = 0.0E0;
  Double I_ERI_H4yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_H3y2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_H2y3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Hy4z_S_Gx3z_S = 0.0E0;
  Double I_ERI_H5z_S_Gx3z_S = 0.0E0;
  Double I_ERI_H5x_S_G4y_S = 0.0E0;
  Double I_ERI_H4xy_S_G4y_S = 0.0E0;
  Double I_ERI_H4xz_S_G4y_S = 0.0E0;
  Double I_ERI_H3x2y_S_G4y_S = 0.0E0;
  Double I_ERI_H3xyz_S_G4y_S = 0.0E0;
  Double I_ERI_H3x2z_S_G4y_S = 0.0E0;
  Double I_ERI_H2x3y_S_G4y_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G4y_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G4y_S = 0.0E0;
  Double I_ERI_H2x3z_S_G4y_S = 0.0E0;
  Double I_ERI_Hx4y_S_G4y_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G4y_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G4y_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G4y_S = 0.0E0;
  Double I_ERI_Hx4z_S_G4y_S = 0.0E0;
  Double I_ERI_H5y_S_G4y_S = 0.0E0;
  Double I_ERI_H4yz_S_G4y_S = 0.0E0;
  Double I_ERI_H3y2z_S_G4y_S = 0.0E0;
  Double I_ERI_H2y3z_S_G4y_S = 0.0E0;
  Double I_ERI_Hy4z_S_G4y_S = 0.0E0;
  Double I_ERI_H5z_S_G4y_S = 0.0E0;
  Double I_ERI_H5x_S_G3yz_S = 0.0E0;
  Double I_ERI_H4xy_S_G3yz_S = 0.0E0;
  Double I_ERI_H4xz_S_G3yz_S = 0.0E0;
  Double I_ERI_H3x2y_S_G3yz_S = 0.0E0;
  Double I_ERI_H3xyz_S_G3yz_S = 0.0E0;
  Double I_ERI_H3x2z_S_G3yz_S = 0.0E0;
  Double I_ERI_H2x3y_S_G3yz_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G3yz_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G3yz_S = 0.0E0;
  Double I_ERI_H2x3z_S_G3yz_S = 0.0E0;
  Double I_ERI_Hx4y_S_G3yz_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G3yz_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G3yz_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G3yz_S = 0.0E0;
  Double I_ERI_Hx4z_S_G3yz_S = 0.0E0;
  Double I_ERI_H5y_S_G3yz_S = 0.0E0;
  Double I_ERI_H4yz_S_G3yz_S = 0.0E0;
  Double I_ERI_H3y2z_S_G3yz_S = 0.0E0;
  Double I_ERI_H2y3z_S_G3yz_S = 0.0E0;
  Double I_ERI_Hy4z_S_G3yz_S = 0.0E0;
  Double I_ERI_H5z_S_G3yz_S = 0.0E0;
  Double I_ERI_H5x_S_G2y2z_S = 0.0E0;
  Double I_ERI_H4xy_S_G2y2z_S = 0.0E0;
  Double I_ERI_H4xz_S_G2y2z_S = 0.0E0;
  Double I_ERI_H3x2y_S_G2y2z_S = 0.0E0;
  Double I_ERI_H3xyz_S_G2y2z_S = 0.0E0;
  Double I_ERI_H3x2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_H2x3y_S_G2y2z_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_H2x3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Hx4y_S_G2y2z_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Hx4z_S_G2y2z_S = 0.0E0;
  Double I_ERI_H5y_S_G2y2z_S = 0.0E0;
  Double I_ERI_H4yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_H3y2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_H2y3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Hy4z_S_G2y2z_S = 0.0E0;
  Double I_ERI_H5z_S_G2y2z_S = 0.0E0;
  Double I_ERI_H5x_S_Gy3z_S = 0.0E0;
  Double I_ERI_H4xy_S_Gy3z_S = 0.0E0;
  Double I_ERI_H4xz_S_Gy3z_S = 0.0E0;
  Double I_ERI_H3x2y_S_Gy3z_S = 0.0E0;
  Double I_ERI_H3xyz_S_Gy3z_S = 0.0E0;
  Double I_ERI_H3x2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_H2x3y_S_Gy3z_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_H2x3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Hx4y_S_Gy3z_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Hx4z_S_Gy3z_S = 0.0E0;
  Double I_ERI_H5y_S_Gy3z_S = 0.0E0;
  Double I_ERI_H4yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_H3y2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_H2y3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Hy4z_S_Gy3z_S = 0.0E0;
  Double I_ERI_H5z_S_Gy3z_S = 0.0E0;
  Double I_ERI_H5x_S_G4z_S = 0.0E0;
  Double I_ERI_H4xy_S_G4z_S = 0.0E0;
  Double I_ERI_H4xz_S_G4z_S = 0.0E0;
  Double I_ERI_H3x2y_S_G4z_S = 0.0E0;
  Double I_ERI_H3xyz_S_G4z_S = 0.0E0;
  Double I_ERI_H3x2z_S_G4z_S = 0.0E0;
  Double I_ERI_H2x3y_S_G4z_S = 0.0E0;
  Double I_ERI_H2x2yz_S_G4z_S = 0.0E0;
  Double I_ERI_H2xy2z_S_G4z_S = 0.0E0;
  Double I_ERI_H2x3z_S_G4z_S = 0.0E0;
  Double I_ERI_Hx4y_S_G4z_S = 0.0E0;
  Double I_ERI_Hx3yz_S_G4z_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_G4z_S = 0.0E0;
  Double I_ERI_Hxy3z_S_G4z_S = 0.0E0;
  Double I_ERI_Hx4z_S_G4z_S = 0.0E0;
  Double I_ERI_H5y_S_G4z_S = 0.0E0;
  Double I_ERI_H4yz_S_G4z_S = 0.0E0;
  Double I_ERI_H3y2z_S_G4z_S = 0.0E0;
  Double I_ERI_H2y3z_S_G4z_S = 0.0E0;
  Double I_ERI_Hy4z_S_G4z_S = 0.0E0;
  Double I_ERI_H5z_S_G4z_S = 0.0E0;
  Double I_ERI_G4x_S_G4x_S = 0.0E0;
  Double I_ERI_G3xy_S_G4x_S = 0.0E0;
  Double I_ERI_G3xz_S_G4x_S = 0.0E0;
  Double I_ERI_G2x2y_S_G4x_S = 0.0E0;
  Double I_ERI_G2xyz_S_G4x_S = 0.0E0;
  Double I_ERI_G2x2z_S_G4x_S = 0.0E0;
  Double I_ERI_Gx3y_S_G4x_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G4x_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G4x_S = 0.0E0;
  Double I_ERI_Gx3z_S_G4x_S = 0.0E0;
  Double I_ERI_G4y_S_G4x_S = 0.0E0;
  Double I_ERI_G3yz_S_G4x_S = 0.0E0;
  Double I_ERI_G2y2z_S_G4x_S = 0.0E0;
  Double I_ERI_Gy3z_S_G4x_S = 0.0E0;
  Double I_ERI_G4z_S_G4x_S = 0.0E0;
  Double I_ERI_G4x_S_G3xy_S = 0.0E0;
  Double I_ERI_G3xy_S_G3xy_S = 0.0E0;
  Double I_ERI_G3xz_S_G3xy_S = 0.0E0;
  Double I_ERI_G2x2y_S_G3xy_S = 0.0E0;
  Double I_ERI_G2xyz_S_G3xy_S = 0.0E0;
  Double I_ERI_G2x2z_S_G3xy_S = 0.0E0;
  Double I_ERI_Gx3y_S_G3xy_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G3xy_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G3xy_S = 0.0E0;
  Double I_ERI_Gx3z_S_G3xy_S = 0.0E0;
  Double I_ERI_G4y_S_G3xy_S = 0.0E0;
  Double I_ERI_G3yz_S_G3xy_S = 0.0E0;
  Double I_ERI_G2y2z_S_G3xy_S = 0.0E0;
  Double I_ERI_Gy3z_S_G3xy_S = 0.0E0;
  Double I_ERI_G4z_S_G3xy_S = 0.0E0;
  Double I_ERI_G4x_S_G3xz_S = 0.0E0;
  Double I_ERI_G3xy_S_G3xz_S = 0.0E0;
  Double I_ERI_G3xz_S_G3xz_S = 0.0E0;
  Double I_ERI_G2x2y_S_G3xz_S = 0.0E0;
  Double I_ERI_G2xyz_S_G3xz_S = 0.0E0;
  Double I_ERI_G2x2z_S_G3xz_S = 0.0E0;
  Double I_ERI_Gx3y_S_G3xz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G3xz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G3xz_S = 0.0E0;
  Double I_ERI_Gx3z_S_G3xz_S = 0.0E0;
  Double I_ERI_G4y_S_G3xz_S = 0.0E0;
  Double I_ERI_G3yz_S_G3xz_S = 0.0E0;
  Double I_ERI_G2y2z_S_G3xz_S = 0.0E0;
  Double I_ERI_Gy3z_S_G3xz_S = 0.0E0;
  Double I_ERI_G4z_S_G3xz_S = 0.0E0;
  Double I_ERI_G4x_S_G2x2y_S = 0.0E0;
  Double I_ERI_G3xy_S_G2x2y_S = 0.0E0;
  Double I_ERI_G3xz_S_G2x2y_S = 0.0E0;
  Double I_ERI_G2x2y_S_G2x2y_S = 0.0E0;
  Double I_ERI_G2xyz_S_G2x2y_S = 0.0E0;
  Double I_ERI_G2x2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Gx3y_S_G2x2y_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Gx3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_G4y_S_G2x2y_S = 0.0E0;
  Double I_ERI_G3yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_G2y2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_Gy3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_G4z_S_G2x2y_S = 0.0E0;
  Double I_ERI_G4x_S_G2xyz_S = 0.0E0;
  Double I_ERI_G3xy_S_G2xyz_S = 0.0E0;
  Double I_ERI_G3xz_S_G2xyz_S = 0.0E0;
  Double I_ERI_G2x2y_S_G2xyz_S = 0.0E0;
  Double I_ERI_G2xyz_S_G2xyz_S = 0.0E0;
  Double I_ERI_G2x2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Gx3y_S_G2xyz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Gx3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_G4y_S_G2xyz_S = 0.0E0;
  Double I_ERI_G3yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_G2y2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_Gy3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_G4z_S_G2xyz_S = 0.0E0;
  Double I_ERI_G4x_S_G2x2z_S = 0.0E0;
  Double I_ERI_G3xy_S_G2x2z_S = 0.0E0;
  Double I_ERI_G3xz_S_G2x2z_S = 0.0E0;
  Double I_ERI_G2x2y_S_G2x2z_S = 0.0E0;
  Double I_ERI_G2xyz_S_G2x2z_S = 0.0E0;
  Double I_ERI_G2x2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Gx3y_S_G2x2z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Gx3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_G4y_S_G2x2z_S = 0.0E0;
  Double I_ERI_G3yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_G2y2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_Gy3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_G4z_S_G2x2z_S = 0.0E0;
  Double I_ERI_G4x_S_Gx3y_S = 0.0E0;
  Double I_ERI_G3xy_S_Gx3y_S = 0.0E0;
  Double I_ERI_G3xz_S_Gx3y_S = 0.0E0;
  Double I_ERI_G2x2y_S_Gx3y_S = 0.0E0;
  Double I_ERI_G2xyz_S_Gx3y_S = 0.0E0;
  Double I_ERI_G2x2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Gx3y_S_Gx3y_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Gx3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_G4y_S_Gx3y_S = 0.0E0;
  Double I_ERI_G3yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_G2y2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_Gy3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_G4z_S_Gx3y_S = 0.0E0;
  Double I_ERI_G4x_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G3xy_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G3xz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G2x2y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G2xyz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G2x2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Gx3y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Gx3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G4y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G3yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G2y2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Gy3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G4z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_G4x_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G3xy_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G3xz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G2x2y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G2xyz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G2x2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Gx3y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Gx3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G4y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G3yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G2y2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Gy3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G4z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_G4x_S_Gx3z_S = 0.0E0;
  Double I_ERI_G3xy_S_Gx3z_S = 0.0E0;
  Double I_ERI_G3xz_S_Gx3z_S = 0.0E0;
  Double I_ERI_G2x2y_S_Gx3z_S = 0.0E0;
  Double I_ERI_G2xyz_S_Gx3z_S = 0.0E0;
  Double I_ERI_G2x2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Gx3y_S_Gx3z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Gx3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_G4y_S_Gx3z_S = 0.0E0;
  Double I_ERI_G3yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_G2y2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_Gy3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_G4z_S_Gx3z_S = 0.0E0;
  Double I_ERI_G4x_S_G4y_S = 0.0E0;
  Double I_ERI_G3xy_S_G4y_S = 0.0E0;
  Double I_ERI_G3xz_S_G4y_S = 0.0E0;
  Double I_ERI_G2x2y_S_G4y_S = 0.0E0;
  Double I_ERI_G2xyz_S_G4y_S = 0.0E0;
  Double I_ERI_G2x2z_S_G4y_S = 0.0E0;
  Double I_ERI_Gx3y_S_G4y_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G4y_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G4y_S = 0.0E0;
  Double I_ERI_Gx3z_S_G4y_S = 0.0E0;
  Double I_ERI_G4y_S_G4y_S = 0.0E0;
  Double I_ERI_G3yz_S_G4y_S = 0.0E0;
  Double I_ERI_G2y2z_S_G4y_S = 0.0E0;
  Double I_ERI_Gy3z_S_G4y_S = 0.0E0;
  Double I_ERI_G4z_S_G4y_S = 0.0E0;
  Double I_ERI_G4x_S_G3yz_S = 0.0E0;
  Double I_ERI_G3xy_S_G3yz_S = 0.0E0;
  Double I_ERI_G3xz_S_G3yz_S = 0.0E0;
  Double I_ERI_G2x2y_S_G3yz_S = 0.0E0;
  Double I_ERI_G2xyz_S_G3yz_S = 0.0E0;
  Double I_ERI_G2x2z_S_G3yz_S = 0.0E0;
  Double I_ERI_Gx3y_S_G3yz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G3yz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G3yz_S = 0.0E0;
  Double I_ERI_Gx3z_S_G3yz_S = 0.0E0;
  Double I_ERI_G4y_S_G3yz_S = 0.0E0;
  Double I_ERI_G3yz_S_G3yz_S = 0.0E0;
  Double I_ERI_G2y2z_S_G3yz_S = 0.0E0;
  Double I_ERI_Gy3z_S_G3yz_S = 0.0E0;
  Double I_ERI_G4z_S_G3yz_S = 0.0E0;
  Double I_ERI_G4x_S_G2y2z_S = 0.0E0;
  Double I_ERI_G3xy_S_G2y2z_S = 0.0E0;
  Double I_ERI_G3xz_S_G2y2z_S = 0.0E0;
  Double I_ERI_G2x2y_S_G2y2z_S = 0.0E0;
  Double I_ERI_G2xyz_S_G2y2z_S = 0.0E0;
  Double I_ERI_G2x2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Gx3y_S_G2y2z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Gx3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_G4y_S_G2y2z_S = 0.0E0;
  Double I_ERI_G3yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_G2y2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_Gy3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_G4z_S_G2y2z_S = 0.0E0;
  Double I_ERI_G4x_S_Gy3z_S = 0.0E0;
  Double I_ERI_G3xy_S_Gy3z_S = 0.0E0;
  Double I_ERI_G3xz_S_Gy3z_S = 0.0E0;
  Double I_ERI_G2x2y_S_Gy3z_S = 0.0E0;
  Double I_ERI_G2xyz_S_Gy3z_S = 0.0E0;
  Double I_ERI_G2x2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Gx3y_S_Gy3z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Gx3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_G4y_S_Gy3z_S = 0.0E0;
  Double I_ERI_G3yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_G2y2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_Gy3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_G4z_S_Gy3z_S = 0.0E0;
  Double I_ERI_G4x_S_G4z_S = 0.0E0;
  Double I_ERI_G3xy_S_G4z_S = 0.0E0;
  Double I_ERI_G3xz_S_G4z_S = 0.0E0;
  Double I_ERI_G2x2y_S_G4z_S = 0.0E0;
  Double I_ERI_G2xyz_S_G4z_S = 0.0E0;
  Double I_ERI_G2x2z_S_G4z_S = 0.0E0;
  Double I_ERI_Gx3y_S_G4z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_G4z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_G4z_S = 0.0E0;
  Double I_ERI_Gx3z_S_G4z_S = 0.0E0;
  Double I_ERI_G4y_S_G4z_S = 0.0E0;
  Double I_ERI_G3yz_S_G4z_S = 0.0E0;
  Double I_ERI_G2y2z_S_G4z_S = 0.0E0;
  Double I_ERI_Gy3z_S_G4z_S = 0.0E0;
  Double I_ERI_G4z_S_G4z_S = 0.0E0;

  // initialize the significance check for VRR part 
  // this will determine that whether we skip the following part 
  bool isSignificant = false;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double onedz = iexp[ip2];
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
      Double I_ERI_S_S_S_S_M8_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M9_vrr  = 0.0E0;
      Double I_ERI_S_S_S_S_M10_vrr  = 0.0E0;

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
      double I_ERI_S_S_S_S_M8_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M9_vrr_d  = 0.0E0;
      double I_ERI_S_S_S_S_M10_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER55;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER53*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER51*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER49*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER47*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER45*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER43*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER25*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = 1.0E0+u2*ONEOVER23*I_ERI_S_S_S_S_M10_vrr;
        I_ERI_S_S_S_S_M10_vrr = ONEOVER21*I_ERI_S_S_S_S_M10_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M10_vrr  = f*I_ERI_S_S_S_S_M10_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ERI_S_S_S_S_M9_vrr  = ONEOVER19*(u2*I_ERI_S_S_S_S_M10_vrr+f);
        I_ERI_S_S_S_S_M8_vrr  = ONEOVER17*(u2*I_ERI_S_S_S_S_M9_vrr+f);
        I_ERI_S_S_S_S_M7_vrr  = ONEOVER15*(u2*I_ERI_S_S_S_S_M8_vrr+f);
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
        I_ERI_S_S_S_S_M8_vrr_d = oneO2u*(15.0E0*I_ERI_S_S_S_S_M7_vrr_d-f);
        I_ERI_S_S_S_S_M9_vrr_d = oneO2u*(17.0E0*I_ERI_S_S_S_S_M8_vrr_d-f);
        I_ERI_S_S_S_S_M10_vrr_d = oneO2u*(19.0E0*I_ERI_S_S_S_S_M9_vrr_d-f);

        // write the double result back to the float var
        I_ERI_S_S_S_S_vrr = static_cast<Double>(I_ERI_S_S_S_S_vrr_d);
        I_ERI_S_S_S_S_M1_vrr = static_cast<Double>(I_ERI_S_S_S_S_M1_vrr_d);
        I_ERI_S_S_S_S_M2_vrr = static_cast<Double>(I_ERI_S_S_S_S_M2_vrr_d);
        I_ERI_S_S_S_S_M3_vrr = static_cast<Double>(I_ERI_S_S_S_S_M3_vrr_d);
        I_ERI_S_S_S_S_M4_vrr = static_cast<Double>(I_ERI_S_S_S_S_M4_vrr_d);
        I_ERI_S_S_S_S_M5_vrr = static_cast<Double>(I_ERI_S_S_S_S_M5_vrr_d);
        I_ERI_S_S_S_S_M6_vrr = static_cast<Double>(I_ERI_S_S_S_S_M6_vrr_d);
        I_ERI_S_S_S_S_M7_vrr = static_cast<Double>(I_ERI_S_S_S_S_M7_vrr_d);
        I_ERI_S_S_S_S_M8_vrr = static_cast<Double>(I_ERI_S_S_S_S_M8_vrr_d);
        I_ERI_S_S_S_S_M9_vrr = static_cast<Double>(I_ERI_S_S_S_S_M9_vrr_d);
        I_ERI_S_S_S_S_M10_vrr = static_cast<Double>(I_ERI_S_S_S_S_M10_vrr_d);

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
        I_ERI_S_S_S_S_M8_vrr = oneO2u*(15.0E0*I_ERI_S_S_S_S_M7_vrr-f);
        I_ERI_S_S_S_S_M9_vrr = oneO2u*(17.0E0*I_ERI_S_S_S_S_M8_vrr-f);
        I_ERI_S_S_S_S_M10_vrr = oneO2u*(19.0E0*I_ERI_S_S_S_S_M9_vrr-f);

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
        Double erfPref_17 = erfPref_15*erfp2;
        Double erfPref_19 = erfPref_17*erfp2;
        Double erfPref_21 = erfPref_19*erfp2;
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
        I_ERI_S_S_S_S_M5_vrr = I_ERI_S_S_S_S_M5_vrr*erfPref_11;
        I_ERI_S_S_S_S_M6_vrr = I_ERI_S_S_S_S_M6_vrr*erfPref_13;
        I_ERI_S_S_S_S_M7_vrr = I_ERI_S_S_S_S_M7_vrr*erfPref_15;
        I_ERI_S_S_S_S_M8_vrr = I_ERI_S_S_S_S_M8_vrr*erfPref_17;
        I_ERI_S_S_S_S_M9_vrr = I_ERI_S_S_S_S_M9_vrr*erfPref_19;
        I_ERI_S_S_S_S_M10_vrr = I_ERI_S_S_S_S_M10_vrr*erfPref_21;
      }

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M9
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M9
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M10
       ************************************************************/
      Double I_ERI_S_S_Px_S_M9_vrr = QCX*I_ERI_S_S_S_S_M9_vrr+WQX*I_ERI_S_S_S_S_M10_vrr;
      Double I_ERI_S_S_Py_S_M9_vrr = QCY*I_ERI_S_S_S_S_M9_vrr+WQY*I_ERI_S_S_S_S_M10_vrr;
      Double I_ERI_S_S_Pz_S_M9_vrr = QCZ*I_ERI_S_S_S_S_M9_vrr+WQZ*I_ERI_S_S_S_S_M10_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M8
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M8
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M9
       ************************************************************/
      Double I_ERI_S_S_Px_S_M8_vrr = QCX*I_ERI_S_S_S_S_M8_vrr+WQX*I_ERI_S_S_S_S_M9_vrr;
      Double I_ERI_S_S_Py_S_M8_vrr = QCY*I_ERI_S_S_S_S_M8_vrr+WQY*I_ERI_S_S_S_S_M9_vrr;
      Double I_ERI_S_S_Pz_S_M8_vrr = QCZ*I_ERI_S_S_S_S_M8_vrr+WQZ*I_ERI_S_S_S_S_M9_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M8
       * expanding position: KET1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M8
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M9
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M8
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M9
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M8_vrr = QCX*I_ERI_S_S_Px_S_M8_vrr+WQX*I_ERI_S_S_Px_S_M9_vrr+oned2e*I_ERI_S_S_S_S_M8_vrr-rhod2esq*I_ERI_S_S_S_S_M9_vrr;
      Double I_ERI_S_S_D2y_S_M8_vrr = QCY*I_ERI_S_S_Py_S_M8_vrr+WQY*I_ERI_S_S_Py_S_M9_vrr+oned2e*I_ERI_S_S_S_S_M8_vrr-rhod2esq*I_ERI_S_S_S_S_M9_vrr;
      Double I_ERI_S_S_D2z_S_M8_vrr = QCZ*I_ERI_S_S_Pz_S_M8_vrr+WQZ*I_ERI_S_S_Pz_S_M9_vrr+oned2e*I_ERI_S_S_S_S_M8_vrr-rhod2esq*I_ERI_S_S_S_S_M9_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M7
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M8
       ************************************************************/
      Double I_ERI_S_S_Px_S_M7_vrr = QCX*I_ERI_S_S_S_S_M7_vrr+WQX*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_S_S_Py_S_M7_vrr = QCY*I_ERI_S_S_S_S_M7_vrr+WQY*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_S_S_Pz_S_M7_vrr = QCZ*I_ERI_S_S_S_S_M7_vrr+WQZ*I_ERI_S_S_S_S_M8_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M7
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M8
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M8
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M7_vrr = QCX*I_ERI_S_S_Px_S_M7_vrr+WQX*I_ERI_S_S_Px_S_M8_vrr+oned2e*I_ERI_S_S_S_S_M7_vrr-rhod2esq*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_S_S_Dxy_S_M7_vrr = QCY*I_ERI_S_S_Px_S_M7_vrr+WQY*I_ERI_S_S_Px_S_M8_vrr;
      Double I_ERI_S_S_Dxz_S_M7_vrr = QCZ*I_ERI_S_S_Px_S_M7_vrr+WQZ*I_ERI_S_S_Px_S_M8_vrr;
      Double I_ERI_S_S_D2y_S_M7_vrr = QCY*I_ERI_S_S_Py_S_M7_vrr+WQY*I_ERI_S_S_Py_S_M8_vrr+oned2e*I_ERI_S_S_S_S_M7_vrr-rhod2esq*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_S_S_Dyz_S_M7_vrr = QCZ*I_ERI_S_S_Py_S_M7_vrr+WQZ*I_ERI_S_S_Py_S_M8_vrr;
      Double I_ERI_S_S_D2z_S_M7_vrr = QCZ*I_ERI_S_S_Pz_S_M7_vrr+WQZ*I_ERI_S_S_Pz_S_M8_vrr+oned2e*I_ERI_S_S_S_S_M7_vrr-rhod2esq*I_ERI_S_S_S_S_M8_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M7
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M8
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M8
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M7_vrr = QCX*I_ERI_S_S_D2x_S_M7_vrr+WQX*I_ERI_S_S_D2x_S_M8_vrr+2*oned2e*I_ERI_S_S_Px_S_M7_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M8_vrr;
      Double I_ERI_S_S_F2xy_S_M7_vrr = QCY*I_ERI_S_S_D2x_S_M7_vrr+WQY*I_ERI_S_S_D2x_S_M8_vrr;
      Double I_ERI_S_S_F2xz_S_M7_vrr = QCZ*I_ERI_S_S_D2x_S_M7_vrr+WQZ*I_ERI_S_S_D2x_S_M8_vrr;
      Double I_ERI_S_S_Fx2y_S_M7_vrr = QCX*I_ERI_S_S_D2y_S_M7_vrr+WQX*I_ERI_S_S_D2y_S_M8_vrr;
      Double I_ERI_S_S_Fx2z_S_M7_vrr = QCX*I_ERI_S_S_D2z_S_M7_vrr+WQX*I_ERI_S_S_D2z_S_M8_vrr;
      Double I_ERI_S_S_F3y_S_M7_vrr = QCY*I_ERI_S_S_D2y_S_M7_vrr+WQY*I_ERI_S_S_D2y_S_M8_vrr+2*oned2e*I_ERI_S_S_Py_S_M7_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M8_vrr;
      Double I_ERI_S_S_F2yz_S_M7_vrr = QCZ*I_ERI_S_S_D2y_S_M7_vrr+WQZ*I_ERI_S_S_D2y_S_M8_vrr;
      Double I_ERI_S_S_F3z_S_M7_vrr = QCZ*I_ERI_S_S_D2z_S_M7_vrr+WQZ*I_ERI_S_S_D2z_S_M8_vrr+2*oned2e*I_ERI_S_S_Pz_S_M7_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M8_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M6
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_S_S_Px_S_M6_vrr = QCX*I_ERI_S_S_S_S_M6_vrr+WQX*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_Py_S_M6_vrr = QCY*I_ERI_S_S_S_S_M6_vrr+WQY*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_Pz_S_M6_vrr = QCZ*I_ERI_S_S_S_S_M6_vrr+WQZ*I_ERI_S_S_S_S_M7_vrr;

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
       * shell quartet name: SQ_ERI_P_S_P_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_Px_S_Px_S_M6_vrr = PAX*I_ERI_S_S_Px_S_M6_vrr+WPX*I_ERI_S_S_Px_S_M7_vrr+oned2k*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_Py_S_Px_S_M6_vrr = PAY*I_ERI_S_S_Px_S_M6_vrr+WPY*I_ERI_S_S_Px_S_M7_vrr;
      Double I_ERI_Pz_S_Px_S_M6_vrr = PAZ*I_ERI_S_S_Px_S_M6_vrr+WPZ*I_ERI_S_S_Px_S_M7_vrr;
      Double I_ERI_Px_S_Py_S_M6_vrr = PAX*I_ERI_S_S_Py_S_M6_vrr+WPX*I_ERI_S_S_Py_S_M7_vrr;
      Double I_ERI_Py_S_Py_S_M6_vrr = PAY*I_ERI_S_S_Py_S_M6_vrr+WPY*I_ERI_S_S_Py_S_M7_vrr+oned2k*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_Pz_S_Py_S_M6_vrr = PAZ*I_ERI_S_S_Py_S_M6_vrr+WPZ*I_ERI_S_S_Py_S_M7_vrr;
      Double I_ERI_Px_S_Pz_S_M6_vrr = PAX*I_ERI_S_S_Pz_S_M6_vrr+WPX*I_ERI_S_S_Pz_S_M7_vrr;
      Double I_ERI_Py_S_Pz_S_M6_vrr = PAY*I_ERI_S_S_Pz_S_M6_vrr+WPY*I_ERI_S_S_Pz_S_M7_vrr;
      Double I_ERI_Pz_S_Pz_S_M6_vrr = PAZ*I_ERI_S_S_Pz_S_M6_vrr+WPZ*I_ERI_S_S_Pz_S_M7_vrr+oned2k*I_ERI_S_S_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M6
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M6_vrr = QCX*I_ERI_S_S_Px_S_M6_vrr+WQX*I_ERI_S_S_Px_S_M7_vrr+oned2e*I_ERI_S_S_S_S_M6_vrr-rhod2esq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_Dxy_S_M6_vrr = QCY*I_ERI_S_S_Px_S_M6_vrr+WQY*I_ERI_S_S_Px_S_M7_vrr;
      Double I_ERI_S_S_Dxz_S_M6_vrr = QCZ*I_ERI_S_S_Px_S_M6_vrr+WQZ*I_ERI_S_S_Px_S_M7_vrr;
      Double I_ERI_S_S_D2y_S_M6_vrr = QCY*I_ERI_S_S_Py_S_M6_vrr+WQY*I_ERI_S_S_Py_S_M7_vrr+oned2e*I_ERI_S_S_S_S_M6_vrr-rhod2esq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_Dyz_S_M6_vrr = QCZ*I_ERI_S_S_Py_S_M6_vrr+WQZ*I_ERI_S_S_Py_S_M7_vrr;
      Double I_ERI_S_S_D2z_S_M6_vrr = QCZ*I_ERI_S_S_Pz_S_M6_vrr+WQZ*I_ERI_S_S_Pz_S_M7_vrr+oned2e*I_ERI_S_S_S_S_M6_vrr-rhod2esq*I_ERI_S_S_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M7
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M6_vrr = PAX*I_ERI_S_S_D2x_S_M6_vrr+WPX*I_ERI_S_S_D2x_S_M7_vrr+2*oned2k*I_ERI_S_S_Px_S_M7_vrr;
      Double I_ERI_Py_S_D2x_S_M6_vrr = PAY*I_ERI_S_S_D2x_S_M6_vrr+WPY*I_ERI_S_S_D2x_S_M7_vrr;
      Double I_ERI_Pz_S_D2x_S_M6_vrr = PAZ*I_ERI_S_S_D2x_S_M6_vrr+WPZ*I_ERI_S_S_D2x_S_M7_vrr;
      Double I_ERI_Px_S_Dxy_S_M6_vrr = PAX*I_ERI_S_S_Dxy_S_M6_vrr+WPX*I_ERI_S_S_Dxy_S_M7_vrr+oned2k*I_ERI_S_S_Py_S_M7_vrr;
      Double I_ERI_Py_S_Dxy_S_M6_vrr = PAY*I_ERI_S_S_Dxy_S_M6_vrr+WPY*I_ERI_S_S_Dxy_S_M7_vrr+oned2k*I_ERI_S_S_Px_S_M7_vrr;
      Double I_ERI_Pz_S_Dxy_S_M6_vrr = PAZ*I_ERI_S_S_Dxy_S_M6_vrr+WPZ*I_ERI_S_S_Dxy_S_M7_vrr;
      Double I_ERI_Px_S_Dxz_S_M6_vrr = PAX*I_ERI_S_S_Dxz_S_M6_vrr+WPX*I_ERI_S_S_Dxz_S_M7_vrr+oned2k*I_ERI_S_S_Pz_S_M7_vrr;
      Double I_ERI_Pz_S_Dxz_S_M6_vrr = PAZ*I_ERI_S_S_Dxz_S_M6_vrr+WPZ*I_ERI_S_S_Dxz_S_M7_vrr+oned2k*I_ERI_S_S_Px_S_M7_vrr;
      Double I_ERI_Px_S_D2y_S_M6_vrr = PAX*I_ERI_S_S_D2y_S_M6_vrr+WPX*I_ERI_S_S_D2y_S_M7_vrr;
      Double I_ERI_Py_S_D2y_S_M6_vrr = PAY*I_ERI_S_S_D2y_S_M6_vrr+WPY*I_ERI_S_S_D2y_S_M7_vrr+2*oned2k*I_ERI_S_S_Py_S_M7_vrr;
      Double I_ERI_Pz_S_D2y_S_M6_vrr = PAZ*I_ERI_S_S_D2y_S_M6_vrr+WPZ*I_ERI_S_S_D2y_S_M7_vrr;
      Double I_ERI_Px_S_Dyz_S_M6_vrr = PAX*I_ERI_S_S_Dyz_S_M6_vrr+WPX*I_ERI_S_S_Dyz_S_M7_vrr;
      Double I_ERI_Py_S_Dyz_S_M6_vrr = PAY*I_ERI_S_S_Dyz_S_M6_vrr+WPY*I_ERI_S_S_Dyz_S_M7_vrr+oned2k*I_ERI_S_S_Pz_S_M7_vrr;
      Double I_ERI_Px_S_D2z_S_M6_vrr = PAX*I_ERI_S_S_D2z_S_M6_vrr+WPX*I_ERI_S_S_D2z_S_M7_vrr;
      Double I_ERI_Py_S_D2z_S_M6_vrr = PAY*I_ERI_S_S_D2z_S_M6_vrr+WPY*I_ERI_S_S_D2z_S_M7_vrr;
      Double I_ERI_Pz_S_D2z_S_M6_vrr = PAZ*I_ERI_S_S_D2z_S_M6_vrr+WPZ*I_ERI_S_S_D2z_S_M7_vrr+2*oned2k*I_ERI_S_S_Pz_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M6
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M7
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M6_vrr = QCX*I_ERI_S_S_D2x_S_M6_vrr+WQX*I_ERI_S_S_D2x_S_M7_vrr+2*oned2e*I_ERI_S_S_Px_S_M6_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M7_vrr;
      Double I_ERI_S_S_F2xy_S_M6_vrr = QCY*I_ERI_S_S_D2x_S_M6_vrr+WQY*I_ERI_S_S_D2x_S_M7_vrr;
      Double I_ERI_S_S_F2xz_S_M6_vrr = QCZ*I_ERI_S_S_D2x_S_M6_vrr+WQZ*I_ERI_S_S_D2x_S_M7_vrr;
      Double I_ERI_S_S_Fx2y_S_M6_vrr = QCX*I_ERI_S_S_D2y_S_M6_vrr+WQX*I_ERI_S_S_D2y_S_M7_vrr;
      Double I_ERI_S_S_Fx2z_S_M6_vrr = QCX*I_ERI_S_S_D2z_S_M6_vrr+WQX*I_ERI_S_S_D2z_S_M7_vrr;
      Double I_ERI_S_S_F3y_S_M6_vrr = QCY*I_ERI_S_S_D2y_S_M6_vrr+WQY*I_ERI_S_S_D2y_S_M7_vrr+2*oned2e*I_ERI_S_S_Py_S_M6_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M7_vrr;
      Double I_ERI_S_S_F2yz_S_M6_vrr = QCZ*I_ERI_S_S_D2y_S_M6_vrr+WQZ*I_ERI_S_S_D2y_S_M7_vrr;
      Double I_ERI_S_S_F3z_S_M6_vrr = QCZ*I_ERI_S_S_D2z_S_M6_vrr+WQZ*I_ERI_S_S_D2z_S_M7_vrr+2*oned2e*I_ERI_S_S_Pz_S_M6_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M7
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M6_vrr = PAX*I_ERI_S_S_F3x_S_M6_vrr+WPX*I_ERI_S_S_F3x_S_M7_vrr+3*oned2k*I_ERI_S_S_D2x_S_M7_vrr;
      Double I_ERI_Py_S_F3x_S_M6_vrr = PAY*I_ERI_S_S_F3x_S_M6_vrr+WPY*I_ERI_S_S_F3x_S_M7_vrr;
      Double I_ERI_Pz_S_F3x_S_M6_vrr = PAZ*I_ERI_S_S_F3x_S_M6_vrr+WPZ*I_ERI_S_S_F3x_S_M7_vrr;
      Double I_ERI_Px_S_F2xy_S_M6_vrr = PAX*I_ERI_S_S_F2xy_S_M6_vrr+WPX*I_ERI_S_S_F2xy_S_M7_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M7_vrr;
      Double I_ERI_Py_S_F2xy_S_M6_vrr = PAY*I_ERI_S_S_F2xy_S_M6_vrr+WPY*I_ERI_S_S_F2xy_S_M7_vrr+oned2k*I_ERI_S_S_D2x_S_M7_vrr;
      Double I_ERI_Pz_S_F2xy_S_M6_vrr = PAZ*I_ERI_S_S_F2xy_S_M6_vrr+WPZ*I_ERI_S_S_F2xy_S_M7_vrr;
      Double I_ERI_Px_S_F2xz_S_M6_vrr = PAX*I_ERI_S_S_F2xz_S_M6_vrr+WPX*I_ERI_S_S_F2xz_S_M7_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M7_vrr;
      Double I_ERI_Py_S_F2xz_S_M6_vrr = PAY*I_ERI_S_S_F2xz_S_M6_vrr+WPY*I_ERI_S_S_F2xz_S_M7_vrr;
      Double I_ERI_Pz_S_F2xz_S_M6_vrr = PAZ*I_ERI_S_S_F2xz_S_M6_vrr+WPZ*I_ERI_S_S_F2xz_S_M7_vrr+oned2k*I_ERI_S_S_D2x_S_M7_vrr;
      Double I_ERI_Px_S_Fx2y_S_M6_vrr = PAX*I_ERI_S_S_Fx2y_S_M6_vrr+WPX*I_ERI_S_S_Fx2y_S_M7_vrr+oned2k*I_ERI_S_S_D2y_S_M7_vrr;
      Double I_ERI_Py_S_Fx2y_S_M6_vrr = PAY*I_ERI_S_S_Fx2y_S_M6_vrr+WPY*I_ERI_S_S_Fx2y_S_M7_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M7_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M6_vrr = PAZ*I_ERI_S_S_Fx2y_S_M6_vrr+WPZ*I_ERI_S_S_Fx2y_S_M7_vrr;
      Double I_ERI_Px_S_Fx2z_S_M6_vrr = PAX*I_ERI_S_S_Fx2z_S_M6_vrr+WPX*I_ERI_S_S_Fx2z_S_M7_vrr+oned2k*I_ERI_S_S_D2z_S_M7_vrr;
      Double I_ERI_Py_S_Fx2z_S_M6_vrr = PAY*I_ERI_S_S_Fx2z_S_M6_vrr+WPY*I_ERI_S_S_Fx2z_S_M7_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M6_vrr = PAZ*I_ERI_S_S_Fx2z_S_M6_vrr+WPZ*I_ERI_S_S_Fx2z_S_M7_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M7_vrr;
      Double I_ERI_Px_S_F3y_S_M6_vrr = PAX*I_ERI_S_S_F3y_S_M6_vrr+WPX*I_ERI_S_S_F3y_S_M7_vrr;
      Double I_ERI_Py_S_F3y_S_M6_vrr = PAY*I_ERI_S_S_F3y_S_M6_vrr+WPY*I_ERI_S_S_F3y_S_M7_vrr+3*oned2k*I_ERI_S_S_D2y_S_M7_vrr;
      Double I_ERI_Pz_S_F3y_S_M6_vrr = PAZ*I_ERI_S_S_F3y_S_M6_vrr+WPZ*I_ERI_S_S_F3y_S_M7_vrr;
      Double I_ERI_Px_S_F2yz_S_M6_vrr = PAX*I_ERI_S_S_F2yz_S_M6_vrr+WPX*I_ERI_S_S_F2yz_S_M7_vrr;
      Double I_ERI_Py_S_F2yz_S_M6_vrr = PAY*I_ERI_S_S_F2yz_S_M6_vrr+WPY*I_ERI_S_S_F2yz_S_M7_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M7_vrr;
      Double I_ERI_Pz_S_F2yz_S_M6_vrr = PAZ*I_ERI_S_S_F2yz_S_M6_vrr+WPZ*I_ERI_S_S_F2yz_S_M7_vrr+oned2k*I_ERI_S_S_D2y_S_M7_vrr;
      Double I_ERI_Px_S_F3z_S_M6_vrr = PAX*I_ERI_S_S_F3z_S_M6_vrr+WPX*I_ERI_S_S_F3z_S_M7_vrr;
      Double I_ERI_Py_S_F3z_S_M6_vrr = PAY*I_ERI_S_S_F3z_S_M6_vrr+WPY*I_ERI_S_S_F3z_S_M7_vrr;
      Double I_ERI_Pz_S_F3z_S_M6_vrr = PAZ*I_ERI_S_S_F3z_S_M6_vrr+WPZ*I_ERI_S_S_F3z_S_M7_vrr+3*oned2k*I_ERI_S_S_D2z_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M5
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_S_S_Px_S_M5_vrr = QCX*I_ERI_S_S_S_S_M5_vrr+WQX*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Py_S_M5_vrr = QCY*I_ERI_S_S_S_S_M5_vrr+WQY*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Pz_S_M5_vrr = QCZ*I_ERI_S_S_S_S_M5_vrr+WQZ*I_ERI_S_S_S_S_M6_vrr;

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
       * shell quartet name: SQ_ERI_P_S_P_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_Px_S_Px_S_M5_vrr = PAX*I_ERI_S_S_Px_S_M5_vrr+WPX*I_ERI_S_S_Px_S_M6_vrr+oned2k*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_Py_S_Px_S_M5_vrr = PAY*I_ERI_S_S_Px_S_M5_vrr+WPY*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_Pz_S_Px_S_M5_vrr = PAZ*I_ERI_S_S_Px_S_M5_vrr+WPZ*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_Px_S_Py_S_M5_vrr = PAX*I_ERI_S_S_Py_S_M5_vrr+WPX*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_Py_S_Py_S_M5_vrr = PAY*I_ERI_S_S_Py_S_M5_vrr+WPY*I_ERI_S_S_Py_S_M6_vrr+oned2k*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_Pz_S_Py_S_M5_vrr = PAZ*I_ERI_S_S_Py_S_M5_vrr+WPZ*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_Px_S_Pz_S_M5_vrr = PAX*I_ERI_S_S_Pz_S_M5_vrr+WPX*I_ERI_S_S_Pz_S_M6_vrr;
      Double I_ERI_Py_S_Pz_S_M5_vrr = PAY*I_ERI_S_S_Pz_S_M5_vrr+WPY*I_ERI_S_S_Pz_S_M6_vrr;
      Double I_ERI_Pz_S_Pz_S_M5_vrr = PAZ*I_ERI_S_S_Pz_S_M5_vrr+WPZ*I_ERI_S_S_Pz_S_M6_vrr+oned2k*I_ERI_S_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M5
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M5_vrr = QCX*I_ERI_S_S_Px_S_M5_vrr+WQX*I_ERI_S_S_Px_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Dxy_S_M5_vrr = QCY*I_ERI_S_S_Px_S_M5_vrr+WQY*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_S_S_Dxz_S_M5_vrr = QCZ*I_ERI_S_S_Px_S_M5_vrr+WQZ*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_S_S_D2y_S_M5_vrr = QCY*I_ERI_S_S_Py_S_M5_vrr+WQY*I_ERI_S_S_Py_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Dyz_S_M5_vrr = QCZ*I_ERI_S_S_Py_S_M5_vrr+WQZ*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_S_S_D2z_S_M5_vrr = QCZ*I_ERI_S_S_Pz_S_M5_vrr+WQZ*I_ERI_S_S_Pz_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M5_vrr = PAX*I_ERI_S_S_D2x_S_M5_vrr+WPX*I_ERI_S_S_D2x_S_M6_vrr+2*oned2k*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_Py_S_D2x_S_M5_vrr = PAY*I_ERI_S_S_D2x_S_M5_vrr+WPY*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_Pz_S_D2x_S_M5_vrr = PAZ*I_ERI_S_S_D2x_S_M5_vrr+WPZ*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_Px_S_Dxy_S_M5_vrr = PAX*I_ERI_S_S_Dxy_S_M5_vrr+WPX*I_ERI_S_S_Dxy_S_M6_vrr+oned2k*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_Py_S_Dxy_S_M5_vrr = PAY*I_ERI_S_S_Dxy_S_M5_vrr+WPY*I_ERI_S_S_Dxy_S_M6_vrr+oned2k*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_Pz_S_Dxy_S_M5_vrr = PAZ*I_ERI_S_S_Dxy_S_M5_vrr+WPZ*I_ERI_S_S_Dxy_S_M6_vrr;
      Double I_ERI_Px_S_Dxz_S_M5_vrr = PAX*I_ERI_S_S_Dxz_S_M5_vrr+WPX*I_ERI_S_S_Dxz_S_M6_vrr+oned2k*I_ERI_S_S_Pz_S_M6_vrr;
      Double I_ERI_Py_S_Dxz_S_M5_vrr = PAY*I_ERI_S_S_Dxz_S_M5_vrr+WPY*I_ERI_S_S_Dxz_S_M6_vrr;
      Double I_ERI_Pz_S_Dxz_S_M5_vrr = PAZ*I_ERI_S_S_Dxz_S_M5_vrr+WPZ*I_ERI_S_S_Dxz_S_M6_vrr+oned2k*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_Px_S_D2y_S_M5_vrr = PAX*I_ERI_S_S_D2y_S_M5_vrr+WPX*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_Py_S_D2y_S_M5_vrr = PAY*I_ERI_S_S_D2y_S_M5_vrr+WPY*I_ERI_S_S_D2y_S_M6_vrr+2*oned2k*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_Pz_S_D2y_S_M5_vrr = PAZ*I_ERI_S_S_D2y_S_M5_vrr+WPZ*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_Px_S_Dyz_S_M5_vrr = PAX*I_ERI_S_S_Dyz_S_M5_vrr+WPX*I_ERI_S_S_Dyz_S_M6_vrr;
      Double I_ERI_Py_S_Dyz_S_M5_vrr = PAY*I_ERI_S_S_Dyz_S_M5_vrr+WPY*I_ERI_S_S_Dyz_S_M6_vrr+oned2k*I_ERI_S_S_Pz_S_M6_vrr;
      Double I_ERI_Pz_S_Dyz_S_M5_vrr = PAZ*I_ERI_S_S_Dyz_S_M5_vrr+WPZ*I_ERI_S_S_Dyz_S_M6_vrr+oned2k*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_Px_S_D2z_S_M5_vrr = PAX*I_ERI_S_S_D2z_S_M5_vrr+WPX*I_ERI_S_S_D2z_S_M6_vrr;
      Double I_ERI_Py_S_D2z_S_M5_vrr = PAY*I_ERI_S_S_D2z_S_M5_vrr+WPY*I_ERI_S_S_D2z_S_M6_vrr;
      Double I_ERI_Pz_S_D2z_S_M5_vrr = PAZ*I_ERI_S_S_D2z_S_M5_vrr+WPZ*I_ERI_S_S_D2z_S_M6_vrr+2*oned2k*I_ERI_S_S_Pz_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M6
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M5_vrr = PAX*I_ERI_Px_S_Px_S_M5_vrr+WPX*I_ERI_Px_S_Px_S_M6_vrr+oned2z*I_ERI_S_S_Px_S_M5_vrr-rhod2zsq*I_ERI_S_S_Px_S_M6_vrr+oned2k*I_ERI_Px_S_S_S_M6_vrr;
      Double I_ERI_Dxy_S_Px_S_M5_vrr = PAY*I_ERI_Px_S_Px_S_M5_vrr+WPY*I_ERI_Px_S_Px_S_M6_vrr;
      Double I_ERI_D2y_S_Px_S_M5_vrr = PAY*I_ERI_Py_S_Px_S_M5_vrr+WPY*I_ERI_Py_S_Px_S_M6_vrr+oned2z*I_ERI_S_S_Px_S_M5_vrr-rhod2zsq*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_D2z_S_Px_S_M5_vrr = PAZ*I_ERI_Pz_S_Px_S_M5_vrr+WPZ*I_ERI_Pz_S_Px_S_M6_vrr+oned2z*I_ERI_S_S_Px_S_M5_vrr-rhod2zsq*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_D2x_S_Py_S_M5_vrr = PAX*I_ERI_Px_S_Py_S_M5_vrr+WPX*I_ERI_Px_S_Py_S_M6_vrr+oned2z*I_ERI_S_S_Py_S_M5_vrr-rhod2zsq*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_Dxy_S_Py_S_M5_vrr = PAY*I_ERI_Px_S_Py_S_M5_vrr+WPY*I_ERI_Px_S_Py_S_M6_vrr+oned2k*I_ERI_Px_S_S_S_M6_vrr;
      Double I_ERI_D2y_S_Py_S_M5_vrr = PAY*I_ERI_Py_S_Py_S_M5_vrr+WPY*I_ERI_Py_S_Py_S_M6_vrr+oned2z*I_ERI_S_S_Py_S_M5_vrr-rhod2zsq*I_ERI_S_S_Py_S_M6_vrr+oned2k*I_ERI_Py_S_S_S_M6_vrr;
      Double I_ERI_D2z_S_Py_S_M5_vrr = PAZ*I_ERI_Pz_S_Py_S_M5_vrr+WPZ*I_ERI_Pz_S_Py_S_M6_vrr+oned2z*I_ERI_S_S_Py_S_M5_vrr-rhod2zsq*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_D2x_S_Pz_S_M5_vrr = PAX*I_ERI_Px_S_Pz_S_M5_vrr+WPX*I_ERI_Px_S_Pz_S_M6_vrr+oned2z*I_ERI_S_S_Pz_S_M5_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M6_vrr;
      Double I_ERI_Dxy_S_Pz_S_M5_vrr = PAY*I_ERI_Px_S_Pz_S_M5_vrr+WPY*I_ERI_Px_S_Pz_S_M6_vrr;
      Double I_ERI_D2y_S_Pz_S_M5_vrr = PAY*I_ERI_Py_S_Pz_S_M5_vrr+WPY*I_ERI_Py_S_Pz_S_M6_vrr+oned2z*I_ERI_S_S_Pz_S_M5_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M6_vrr;
      Double I_ERI_D2z_S_Pz_S_M5_vrr = PAZ*I_ERI_Pz_S_Pz_S_M5_vrr+WPZ*I_ERI_Pz_S_Pz_S_M6_vrr+oned2z*I_ERI_S_S_Pz_S_M5_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M6_vrr+oned2k*I_ERI_Pz_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M5
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M5_vrr = QCX*I_ERI_S_S_D2x_S_M5_vrr+WQX*I_ERI_S_S_D2x_S_M6_vrr+2*oned2e*I_ERI_S_S_Px_S_M5_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_S_S_F2xy_S_M5_vrr = QCY*I_ERI_S_S_D2x_S_M5_vrr+WQY*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_S_S_F2xz_S_M5_vrr = QCZ*I_ERI_S_S_D2x_S_M5_vrr+WQZ*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_S_S_Fx2y_S_M5_vrr = QCX*I_ERI_S_S_D2y_S_M5_vrr+WQX*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_S_S_Fx2z_S_M5_vrr = QCX*I_ERI_S_S_D2z_S_M5_vrr+WQX*I_ERI_S_S_D2z_S_M6_vrr;
      Double I_ERI_S_S_F3y_S_M5_vrr = QCY*I_ERI_S_S_D2y_S_M5_vrr+WQY*I_ERI_S_S_D2y_S_M6_vrr+2*oned2e*I_ERI_S_S_Py_S_M5_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M6_vrr;
      Double I_ERI_S_S_F2yz_S_M5_vrr = QCZ*I_ERI_S_S_D2y_S_M5_vrr+WQZ*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_S_S_F3z_S_M5_vrr = QCZ*I_ERI_S_S_D2z_S_M5_vrr+WQZ*I_ERI_S_S_D2z_S_M6_vrr+2*oned2e*I_ERI_S_S_Pz_S_M5_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 17 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M6
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M6
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M5_vrr = PAX*I_ERI_Px_S_D2x_S_M5_vrr+WPX*I_ERI_Px_S_D2x_S_M6_vrr+oned2z*I_ERI_S_S_D2x_S_M5_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M6_vrr+2*oned2k*I_ERI_Px_S_Px_S_M6_vrr;
      Double I_ERI_Dxy_S_D2x_S_M5_vrr = PAY*I_ERI_Px_S_D2x_S_M5_vrr+WPY*I_ERI_Px_S_D2x_S_M6_vrr;
      Double I_ERI_D2y_S_D2x_S_M5_vrr = PAY*I_ERI_Py_S_D2x_S_M5_vrr+WPY*I_ERI_Py_S_D2x_S_M6_vrr+oned2z*I_ERI_S_S_D2x_S_M5_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_D2z_S_D2x_S_M5_vrr = PAZ*I_ERI_Pz_S_D2x_S_M5_vrr+WPZ*I_ERI_Pz_S_D2x_S_M6_vrr+oned2z*I_ERI_S_S_D2x_S_M5_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_D2x_S_Dxy_S_M5_vrr = PAX*I_ERI_Px_S_Dxy_S_M5_vrr+WPX*I_ERI_Px_S_Dxy_S_M6_vrr+oned2z*I_ERI_S_S_Dxy_S_M5_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M6_vrr+oned2k*I_ERI_Px_S_Py_S_M6_vrr;
      Double I_ERI_D2y_S_Dxy_S_M5_vrr = PAY*I_ERI_Py_S_Dxy_S_M5_vrr+WPY*I_ERI_Py_S_Dxy_S_M6_vrr+oned2z*I_ERI_S_S_Dxy_S_M5_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M6_vrr+oned2k*I_ERI_Py_S_Px_S_M6_vrr;
      Double I_ERI_D2z_S_Dxy_S_M5_vrr = PAZ*I_ERI_Pz_S_Dxy_S_M5_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M6_vrr+oned2z*I_ERI_S_S_Dxy_S_M5_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M6_vrr;
      Double I_ERI_D2x_S_Dxz_S_M5_vrr = PAX*I_ERI_Px_S_Dxz_S_M5_vrr+WPX*I_ERI_Px_S_Dxz_S_M6_vrr+oned2z*I_ERI_S_S_Dxz_S_M5_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M6_vrr+oned2k*I_ERI_Px_S_Pz_S_M6_vrr;
      Double I_ERI_D2z_S_Dxz_S_M5_vrr = PAZ*I_ERI_Pz_S_Dxz_S_M5_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M6_vrr+oned2z*I_ERI_S_S_Dxz_S_M5_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M6_vrr+oned2k*I_ERI_Pz_S_Px_S_M6_vrr;
      Double I_ERI_D2x_S_D2y_S_M5_vrr = PAX*I_ERI_Px_S_D2y_S_M5_vrr+WPX*I_ERI_Px_S_D2y_S_M6_vrr+oned2z*I_ERI_S_S_D2y_S_M5_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_Dxy_S_D2y_S_M5_vrr = PAY*I_ERI_Px_S_D2y_S_M5_vrr+WPY*I_ERI_Px_S_D2y_S_M6_vrr+2*oned2k*I_ERI_Px_S_Py_S_M6_vrr;
      Double I_ERI_D2y_S_D2y_S_M5_vrr = PAY*I_ERI_Py_S_D2y_S_M5_vrr+WPY*I_ERI_Py_S_D2y_S_M6_vrr+oned2z*I_ERI_S_S_D2y_S_M5_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M6_vrr+2*oned2k*I_ERI_Py_S_Py_S_M6_vrr;
      Double I_ERI_D2z_S_D2y_S_M5_vrr = PAZ*I_ERI_Pz_S_D2y_S_M5_vrr+WPZ*I_ERI_Pz_S_D2y_S_M6_vrr+oned2z*I_ERI_S_S_D2y_S_M5_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_D2x_S_Dyz_S_M5_vrr = PAX*I_ERI_Px_S_Dyz_S_M5_vrr+WPX*I_ERI_Px_S_Dyz_S_M6_vrr+oned2z*I_ERI_S_S_Dyz_S_M5_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M6_vrr;
      Double I_ERI_D2y_S_Dyz_S_M5_vrr = PAY*I_ERI_Py_S_Dyz_S_M5_vrr+WPY*I_ERI_Py_S_Dyz_S_M6_vrr+oned2z*I_ERI_S_S_Dyz_S_M5_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M6_vrr+oned2k*I_ERI_Py_S_Pz_S_M6_vrr;
      Double I_ERI_D2x_S_D2z_S_M5_vrr = PAX*I_ERI_Px_S_D2z_S_M5_vrr+WPX*I_ERI_Px_S_D2z_S_M6_vrr+oned2z*I_ERI_S_S_D2z_S_M5_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M6_vrr;
      Double I_ERI_Dxy_S_D2z_S_M5_vrr = PAY*I_ERI_Px_S_D2z_S_M5_vrr+WPY*I_ERI_Px_S_D2z_S_M6_vrr;
      Double I_ERI_D2y_S_D2z_S_M5_vrr = PAY*I_ERI_Py_S_D2z_S_M5_vrr+WPY*I_ERI_Py_S_D2z_S_M6_vrr+oned2z*I_ERI_S_S_D2z_S_M5_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M6_vrr;
      Double I_ERI_D2z_S_D2z_S_M5_vrr = PAZ*I_ERI_Pz_S_D2z_S_M5_vrr+WPZ*I_ERI_Pz_S_D2z_S_M6_vrr+oned2z*I_ERI_S_S_D2z_S_M5_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M6_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M6
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M5_vrr = PAX*I_ERI_S_S_F3x_S_M5_vrr+WPX*I_ERI_S_S_F3x_S_M6_vrr+3*oned2k*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_Py_S_F3x_S_M5_vrr = PAY*I_ERI_S_S_F3x_S_M5_vrr+WPY*I_ERI_S_S_F3x_S_M6_vrr;
      Double I_ERI_Pz_S_F3x_S_M5_vrr = PAZ*I_ERI_S_S_F3x_S_M5_vrr+WPZ*I_ERI_S_S_F3x_S_M6_vrr;
      Double I_ERI_Px_S_F2xy_S_M5_vrr = PAX*I_ERI_S_S_F2xy_S_M5_vrr+WPX*I_ERI_S_S_F2xy_S_M6_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M6_vrr;
      Double I_ERI_Py_S_F2xy_S_M5_vrr = PAY*I_ERI_S_S_F2xy_S_M5_vrr+WPY*I_ERI_S_S_F2xy_S_M6_vrr+oned2k*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_Pz_S_F2xy_S_M5_vrr = PAZ*I_ERI_S_S_F2xy_S_M5_vrr+WPZ*I_ERI_S_S_F2xy_S_M6_vrr;
      Double I_ERI_Px_S_F2xz_S_M5_vrr = PAX*I_ERI_S_S_F2xz_S_M5_vrr+WPX*I_ERI_S_S_F2xz_S_M6_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M6_vrr;
      Double I_ERI_Py_S_F2xz_S_M5_vrr = PAY*I_ERI_S_S_F2xz_S_M5_vrr+WPY*I_ERI_S_S_F2xz_S_M6_vrr;
      Double I_ERI_Pz_S_F2xz_S_M5_vrr = PAZ*I_ERI_S_S_F2xz_S_M5_vrr+WPZ*I_ERI_S_S_F2xz_S_M6_vrr+oned2k*I_ERI_S_S_D2x_S_M6_vrr;
      Double I_ERI_Px_S_Fx2y_S_M5_vrr = PAX*I_ERI_S_S_Fx2y_S_M5_vrr+WPX*I_ERI_S_S_Fx2y_S_M6_vrr+oned2k*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_Py_S_Fx2y_S_M5_vrr = PAY*I_ERI_S_S_Fx2y_S_M5_vrr+WPY*I_ERI_S_S_Fx2y_S_M6_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M6_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M5_vrr = PAZ*I_ERI_S_S_Fx2y_S_M5_vrr+WPZ*I_ERI_S_S_Fx2y_S_M6_vrr;
      Double I_ERI_Px_S_Fx2z_S_M5_vrr = PAX*I_ERI_S_S_Fx2z_S_M5_vrr+WPX*I_ERI_S_S_Fx2z_S_M6_vrr+oned2k*I_ERI_S_S_D2z_S_M6_vrr;
      Double I_ERI_Py_S_Fx2z_S_M5_vrr = PAY*I_ERI_S_S_Fx2z_S_M5_vrr+WPY*I_ERI_S_S_Fx2z_S_M6_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M5_vrr = PAZ*I_ERI_S_S_Fx2z_S_M5_vrr+WPZ*I_ERI_S_S_Fx2z_S_M6_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M6_vrr;
      Double I_ERI_Px_S_F3y_S_M5_vrr = PAX*I_ERI_S_S_F3y_S_M5_vrr+WPX*I_ERI_S_S_F3y_S_M6_vrr;
      Double I_ERI_Py_S_F3y_S_M5_vrr = PAY*I_ERI_S_S_F3y_S_M5_vrr+WPY*I_ERI_S_S_F3y_S_M6_vrr+3*oned2k*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_Pz_S_F3y_S_M5_vrr = PAZ*I_ERI_S_S_F3y_S_M5_vrr+WPZ*I_ERI_S_S_F3y_S_M6_vrr;
      Double I_ERI_Px_S_F2yz_S_M5_vrr = PAX*I_ERI_S_S_F2yz_S_M5_vrr+WPX*I_ERI_S_S_F2yz_S_M6_vrr;
      Double I_ERI_Py_S_F2yz_S_M5_vrr = PAY*I_ERI_S_S_F2yz_S_M5_vrr+WPY*I_ERI_S_S_F2yz_S_M6_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M6_vrr;
      Double I_ERI_Pz_S_F2yz_S_M5_vrr = PAZ*I_ERI_S_S_F2yz_S_M5_vrr+WPZ*I_ERI_S_S_F2yz_S_M6_vrr+oned2k*I_ERI_S_S_D2y_S_M6_vrr;
      Double I_ERI_Px_S_F3z_S_M5_vrr = PAX*I_ERI_S_S_F3z_S_M5_vrr+WPX*I_ERI_S_S_F3z_S_M6_vrr;
      Double I_ERI_Py_S_F3z_S_M5_vrr = PAY*I_ERI_S_S_F3z_S_M5_vrr+WPY*I_ERI_S_S_F3z_S_M6_vrr;
      Double I_ERI_Pz_S_F3z_S_M5_vrr = PAZ*I_ERI_S_S_F3z_S_M5_vrr+WPZ*I_ERI_S_S_F3z_S_M6_vrr+3*oned2k*I_ERI_S_S_D2z_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 36 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M6
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M6
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M5_vrr = PAX*I_ERI_Px_S_F3x_S_M5_vrr+WPX*I_ERI_Px_S_F3x_S_M6_vrr+oned2z*I_ERI_S_S_F3x_S_M5_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M6_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M6_vrr;
      Double I_ERI_D2y_S_F3x_S_M5_vrr = PAY*I_ERI_Py_S_F3x_S_M5_vrr+WPY*I_ERI_Py_S_F3x_S_M6_vrr+oned2z*I_ERI_S_S_F3x_S_M5_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M6_vrr;
      Double I_ERI_D2z_S_F3x_S_M5_vrr = PAZ*I_ERI_Pz_S_F3x_S_M5_vrr+WPZ*I_ERI_Pz_S_F3x_S_M6_vrr+oned2z*I_ERI_S_S_F3x_S_M5_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M6_vrr;
      Double I_ERI_D2x_S_F2xy_S_M5_vrr = PAX*I_ERI_Px_S_F2xy_S_M5_vrr+WPX*I_ERI_Px_S_F2xy_S_M6_vrr+oned2z*I_ERI_S_S_F2xy_S_M5_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M6_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M6_vrr;
      Double I_ERI_D2y_S_F2xy_S_M5_vrr = PAY*I_ERI_Py_S_F2xy_S_M5_vrr+WPY*I_ERI_Py_S_F2xy_S_M6_vrr+oned2z*I_ERI_S_S_F2xy_S_M5_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M6_vrr+oned2k*I_ERI_Py_S_D2x_S_M6_vrr;
      Double I_ERI_D2z_S_F2xy_S_M5_vrr = PAZ*I_ERI_Pz_S_F2xy_S_M5_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M6_vrr+oned2z*I_ERI_S_S_F2xy_S_M5_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M6_vrr;
      Double I_ERI_D2x_S_F2xz_S_M5_vrr = PAX*I_ERI_Px_S_F2xz_S_M5_vrr+WPX*I_ERI_Px_S_F2xz_S_M6_vrr+oned2z*I_ERI_S_S_F2xz_S_M5_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M6_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M6_vrr;
      Double I_ERI_D2y_S_F2xz_S_M5_vrr = PAY*I_ERI_Py_S_F2xz_S_M5_vrr+WPY*I_ERI_Py_S_F2xz_S_M6_vrr+oned2z*I_ERI_S_S_F2xz_S_M5_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M6_vrr;
      Double I_ERI_D2z_S_F2xz_S_M5_vrr = PAZ*I_ERI_Pz_S_F2xz_S_M5_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M6_vrr+oned2z*I_ERI_S_S_F2xz_S_M5_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M6_vrr+oned2k*I_ERI_Pz_S_D2x_S_M6_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M5_vrr = PAX*I_ERI_Px_S_Fx2y_S_M5_vrr+WPX*I_ERI_Px_S_Fx2y_S_M6_vrr+oned2z*I_ERI_S_S_Fx2y_S_M5_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M6_vrr+oned2k*I_ERI_Px_S_D2y_S_M6_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M5_vrr = PAY*I_ERI_Py_S_Fx2y_S_M5_vrr+WPY*I_ERI_Py_S_Fx2y_S_M6_vrr+oned2z*I_ERI_S_S_Fx2y_S_M5_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M6_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M6_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M5_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_M5_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M6_vrr+oned2z*I_ERI_S_S_Fx2y_S_M5_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M6_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M5_vrr = PAX*I_ERI_Px_S_Fx2z_S_M5_vrr+WPX*I_ERI_Px_S_Fx2z_S_M6_vrr+oned2z*I_ERI_S_S_Fx2z_S_M5_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M6_vrr+oned2k*I_ERI_Px_S_D2z_S_M6_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M5_vrr = PAY*I_ERI_Py_S_Fx2z_S_M5_vrr+WPY*I_ERI_Py_S_Fx2z_S_M6_vrr+oned2z*I_ERI_S_S_Fx2z_S_M5_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M6_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M5_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_M5_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M6_vrr+oned2z*I_ERI_S_S_Fx2z_S_M5_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M6_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M6_vrr;
      Double I_ERI_D2x_S_F3y_S_M5_vrr = PAX*I_ERI_Px_S_F3y_S_M5_vrr+WPX*I_ERI_Px_S_F3y_S_M6_vrr+oned2z*I_ERI_S_S_F3y_S_M5_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M6_vrr;
      Double I_ERI_D2y_S_F3y_S_M5_vrr = PAY*I_ERI_Py_S_F3y_S_M5_vrr+WPY*I_ERI_Py_S_F3y_S_M6_vrr+oned2z*I_ERI_S_S_F3y_S_M5_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M6_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M6_vrr;
      Double I_ERI_D2z_S_F3y_S_M5_vrr = PAZ*I_ERI_Pz_S_F3y_S_M5_vrr+WPZ*I_ERI_Pz_S_F3y_S_M6_vrr+oned2z*I_ERI_S_S_F3y_S_M5_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M6_vrr;
      Double I_ERI_D2x_S_F2yz_S_M5_vrr = PAX*I_ERI_Px_S_F2yz_S_M5_vrr+WPX*I_ERI_Px_S_F2yz_S_M6_vrr+oned2z*I_ERI_S_S_F2yz_S_M5_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M6_vrr;
      Double I_ERI_D2y_S_F2yz_S_M5_vrr = PAY*I_ERI_Py_S_F2yz_S_M5_vrr+WPY*I_ERI_Py_S_F2yz_S_M6_vrr+oned2z*I_ERI_S_S_F2yz_S_M5_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M6_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M6_vrr;
      Double I_ERI_D2z_S_F2yz_S_M5_vrr = PAZ*I_ERI_Pz_S_F2yz_S_M5_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M6_vrr+oned2z*I_ERI_S_S_F2yz_S_M5_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M6_vrr+oned2k*I_ERI_Pz_S_D2y_S_M6_vrr;
      Double I_ERI_D2x_S_F3z_S_M5_vrr = PAX*I_ERI_Px_S_F3z_S_M5_vrr+WPX*I_ERI_Px_S_F3z_S_M6_vrr+oned2z*I_ERI_S_S_F3z_S_M5_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M6_vrr;
      Double I_ERI_D2y_S_F3z_S_M5_vrr = PAY*I_ERI_Py_S_F3z_S_M5_vrr+WPY*I_ERI_Py_S_F3z_S_M6_vrr+oned2z*I_ERI_S_S_F3z_S_M5_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M6_vrr;
      Double I_ERI_D2z_S_F3z_S_M5_vrr = PAZ*I_ERI_Pz_S_F3z_S_M5_vrr+WPZ*I_ERI_Pz_S_F3z_S_M6_vrr+oned2z*I_ERI_S_S_F3z_S_M5_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M6_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_P_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_S_S_Px_S_M4_vrr = QCX*I_ERI_S_S_S_S_M4_vrr+WQX*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Py_S_M4_vrr = QCY*I_ERI_S_S_S_S_M4_vrr+WQY*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Pz_S_M4_vrr = QCZ*I_ERI_S_S_S_S_M4_vrr+WQZ*I_ERI_S_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_P_S_P_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_Px_S_Px_S_M4_vrr = PAX*I_ERI_S_S_Px_S_M4_vrr+WPX*I_ERI_S_S_Px_S_M5_vrr+oned2k*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_Py_S_Px_S_M4_vrr = PAY*I_ERI_S_S_Px_S_M4_vrr+WPY*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_Pz_S_Px_S_M4_vrr = PAZ*I_ERI_S_S_Px_S_M4_vrr+WPZ*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_Px_S_Py_S_M4_vrr = PAX*I_ERI_S_S_Py_S_M4_vrr+WPX*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_Py_S_Py_S_M4_vrr = PAY*I_ERI_S_S_Py_S_M4_vrr+WPY*I_ERI_S_S_Py_S_M5_vrr+oned2k*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_Pz_S_Py_S_M4_vrr = PAZ*I_ERI_S_S_Py_S_M4_vrr+WPZ*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_Px_S_Pz_S_M4_vrr = PAX*I_ERI_S_S_Pz_S_M4_vrr+WPX*I_ERI_S_S_Pz_S_M5_vrr;
      Double I_ERI_Py_S_Pz_S_M4_vrr = PAY*I_ERI_S_S_Pz_S_M4_vrr+WPY*I_ERI_S_S_Pz_S_M5_vrr;
      Double I_ERI_Pz_S_Pz_S_M4_vrr = PAZ*I_ERI_S_S_Pz_S_M4_vrr+WPZ*I_ERI_S_S_Pz_S_M5_vrr+oned2k*I_ERI_S_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M4_vrr = QCX*I_ERI_S_S_Px_S_M4_vrr+WQX*I_ERI_S_S_Px_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Dxy_S_M4_vrr = QCY*I_ERI_S_S_Px_S_M4_vrr+WQY*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_Dxz_S_M4_vrr = QCZ*I_ERI_S_S_Px_S_M4_vrr+WQZ*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_D2y_S_M4_vrr = QCY*I_ERI_S_S_Py_S_M4_vrr+WQY*I_ERI_S_S_Py_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Dyz_S_M4_vrr = QCZ*I_ERI_S_S_Py_S_M4_vrr+WQZ*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_S_S_D2z_S_M4_vrr = QCZ*I_ERI_S_S_Pz_S_M4_vrr+WQZ*I_ERI_S_S_Pz_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_P_S_D_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M4_vrr = PAX*I_ERI_S_S_D2x_S_M4_vrr+WPX*I_ERI_S_S_D2x_S_M5_vrr+2*oned2k*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_Py_S_D2x_S_M4_vrr = PAY*I_ERI_S_S_D2x_S_M4_vrr+WPY*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_Pz_S_D2x_S_M4_vrr = PAZ*I_ERI_S_S_D2x_S_M4_vrr+WPZ*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_Px_S_Dxy_S_M4_vrr = PAX*I_ERI_S_S_Dxy_S_M4_vrr+WPX*I_ERI_S_S_Dxy_S_M5_vrr+oned2k*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_Py_S_Dxy_S_M4_vrr = PAY*I_ERI_S_S_Dxy_S_M4_vrr+WPY*I_ERI_S_S_Dxy_S_M5_vrr+oned2k*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_Pz_S_Dxy_S_M4_vrr = PAZ*I_ERI_S_S_Dxy_S_M4_vrr+WPZ*I_ERI_S_S_Dxy_S_M5_vrr;
      Double I_ERI_Px_S_Dxz_S_M4_vrr = PAX*I_ERI_S_S_Dxz_S_M4_vrr+WPX*I_ERI_S_S_Dxz_S_M5_vrr+oned2k*I_ERI_S_S_Pz_S_M5_vrr;
      Double I_ERI_Py_S_Dxz_S_M4_vrr = PAY*I_ERI_S_S_Dxz_S_M4_vrr+WPY*I_ERI_S_S_Dxz_S_M5_vrr;
      Double I_ERI_Pz_S_Dxz_S_M4_vrr = PAZ*I_ERI_S_S_Dxz_S_M4_vrr+WPZ*I_ERI_S_S_Dxz_S_M5_vrr+oned2k*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_Px_S_D2y_S_M4_vrr = PAX*I_ERI_S_S_D2y_S_M4_vrr+WPX*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_Py_S_D2y_S_M4_vrr = PAY*I_ERI_S_S_D2y_S_M4_vrr+WPY*I_ERI_S_S_D2y_S_M5_vrr+2*oned2k*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_Pz_S_D2y_S_M4_vrr = PAZ*I_ERI_S_S_D2y_S_M4_vrr+WPZ*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_Px_S_Dyz_S_M4_vrr = PAX*I_ERI_S_S_Dyz_S_M4_vrr+WPX*I_ERI_S_S_Dyz_S_M5_vrr;
      Double I_ERI_Py_S_Dyz_S_M4_vrr = PAY*I_ERI_S_S_Dyz_S_M4_vrr+WPY*I_ERI_S_S_Dyz_S_M5_vrr+oned2k*I_ERI_S_S_Pz_S_M5_vrr;
      Double I_ERI_Pz_S_Dyz_S_M4_vrr = PAZ*I_ERI_S_S_Dyz_S_M4_vrr+WPZ*I_ERI_S_S_Dyz_S_M5_vrr+oned2k*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_Px_S_D2z_S_M4_vrr = PAX*I_ERI_S_S_D2z_S_M4_vrr+WPX*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_Py_S_D2z_S_M4_vrr = PAY*I_ERI_S_S_D2z_S_M4_vrr+WPY*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_Pz_S_D2z_S_M4_vrr = PAZ*I_ERI_S_S_D2z_S_M4_vrr+WPZ*I_ERI_S_S_D2z_S_M5_vrr+2*oned2k*I_ERI_S_S_Pz_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M4_vrr = PAX*I_ERI_Px_S_Px_S_M4_vrr+WPX*I_ERI_Px_S_Px_S_M5_vrr+oned2z*I_ERI_S_S_Px_S_M4_vrr-rhod2zsq*I_ERI_S_S_Px_S_M5_vrr+oned2k*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_Dxy_S_Px_S_M4_vrr = PAY*I_ERI_Px_S_Px_S_M4_vrr+WPY*I_ERI_Px_S_Px_S_M5_vrr;
      Double I_ERI_Dxz_S_Px_S_M4_vrr = PAZ*I_ERI_Px_S_Px_S_M4_vrr+WPZ*I_ERI_Px_S_Px_S_M5_vrr;
      Double I_ERI_D2y_S_Px_S_M4_vrr = PAY*I_ERI_Py_S_Px_S_M4_vrr+WPY*I_ERI_Py_S_Px_S_M5_vrr+oned2z*I_ERI_S_S_Px_S_M4_vrr-rhod2zsq*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_Dyz_S_Px_S_M4_vrr = PAZ*I_ERI_Py_S_Px_S_M4_vrr+WPZ*I_ERI_Py_S_Px_S_M5_vrr;
      Double I_ERI_D2z_S_Px_S_M4_vrr = PAZ*I_ERI_Pz_S_Px_S_M4_vrr+WPZ*I_ERI_Pz_S_Px_S_M5_vrr+oned2z*I_ERI_S_S_Px_S_M4_vrr-rhod2zsq*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_D2x_S_Py_S_M4_vrr = PAX*I_ERI_Px_S_Py_S_M4_vrr+WPX*I_ERI_Px_S_Py_S_M5_vrr+oned2z*I_ERI_S_S_Py_S_M4_vrr-rhod2zsq*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_Dxy_S_Py_S_M4_vrr = PAY*I_ERI_Px_S_Py_S_M4_vrr+WPY*I_ERI_Px_S_Py_S_M5_vrr+oned2k*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_Dxz_S_Py_S_M4_vrr = PAZ*I_ERI_Px_S_Py_S_M4_vrr+WPZ*I_ERI_Px_S_Py_S_M5_vrr;
      Double I_ERI_D2y_S_Py_S_M4_vrr = PAY*I_ERI_Py_S_Py_S_M4_vrr+WPY*I_ERI_Py_S_Py_S_M5_vrr+oned2z*I_ERI_S_S_Py_S_M4_vrr-rhod2zsq*I_ERI_S_S_Py_S_M5_vrr+oned2k*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_Dyz_S_Py_S_M4_vrr = PAZ*I_ERI_Py_S_Py_S_M4_vrr+WPZ*I_ERI_Py_S_Py_S_M5_vrr;
      Double I_ERI_D2z_S_Py_S_M4_vrr = PAZ*I_ERI_Pz_S_Py_S_M4_vrr+WPZ*I_ERI_Pz_S_Py_S_M5_vrr+oned2z*I_ERI_S_S_Py_S_M4_vrr-rhod2zsq*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_D2x_S_Pz_S_M4_vrr = PAX*I_ERI_Px_S_Pz_S_M4_vrr+WPX*I_ERI_Px_S_Pz_S_M5_vrr+oned2z*I_ERI_S_S_Pz_S_M4_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M5_vrr;
      Double I_ERI_Dxy_S_Pz_S_M4_vrr = PAY*I_ERI_Px_S_Pz_S_M4_vrr+WPY*I_ERI_Px_S_Pz_S_M5_vrr;
      Double I_ERI_Dxz_S_Pz_S_M4_vrr = PAZ*I_ERI_Px_S_Pz_S_M4_vrr+WPZ*I_ERI_Px_S_Pz_S_M5_vrr+oned2k*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_D2y_S_Pz_S_M4_vrr = PAY*I_ERI_Py_S_Pz_S_M4_vrr+WPY*I_ERI_Py_S_Pz_S_M5_vrr+oned2z*I_ERI_S_S_Pz_S_M4_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M5_vrr;
      Double I_ERI_Dyz_S_Pz_S_M4_vrr = PAZ*I_ERI_Py_S_Pz_S_M4_vrr+WPZ*I_ERI_Py_S_Pz_S_M5_vrr+oned2k*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_D2z_S_Pz_S_M4_vrr = PAZ*I_ERI_Pz_S_Pz_S_M4_vrr+WPZ*I_ERI_Pz_S_Pz_S_M5_vrr+oned2z*I_ERI_S_S_Pz_S_M4_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M5_vrr+oned2k*I_ERI_Pz_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M5
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M4_vrr = PAX*I_ERI_Px_S_D2x_S_M4_vrr+WPX*I_ERI_Px_S_D2x_S_M5_vrr+oned2z*I_ERI_S_S_D2x_S_M4_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M5_vrr+2*oned2k*I_ERI_Px_S_Px_S_M5_vrr;
      Double I_ERI_Dxy_S_D2x_S_M4_vrr = PAY*I_ERI_Px_S_D2x_S_M4_vrr+WPY*I_ERI_Px_S_D2x_S_M5_vrr;
      Double I_ERI_Dxz_S_D2x_S_M4_vrr = PAZ*I_ERI_Px_S_D2x_S_M4_vrr+WPZ*I_ERI_Px_S_D2x_S_M5_vrr;
      Double I_ERI_D2y_S_D2x_S_M4_vrr = PAY*I_ERI_Py_S_D2x_S_M4_vrr+WPY*I_ERI_Py_S_D2x_S_M5_vrr+oned2z*I_ERI_S_S_D2x_S_M4_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_Dyz_S_D2x_S_M4_vrr = PAZ*I_ERI_Py_S_D2x_S_M4_vrr+WPZ*I_ERI_Py_S_D2x_S_M5_vrr;
      Double I_ERI_D2z_S_D2x_S_M4_vrr = PAZ*I_ERI_Pz_S_D2x_S_M4_vrr+WPZ*I_ERI_Pz_S_D2x_S_M5_vrr+oned2z*I_ERI_S_S_D2x_S_M4_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_D2x_S_Dxy_S_M4_vrr = PAX*I_ERI_Px_S_Dxy_S_M4_vrr+WPX*I_ERI_Px_S_Dxy_S_M5_vrr+oned2z*I_ERI_S_S_Dxy_S_M4_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M5_vrr+oned2k*I_ERI_Px_S_Py_S_M5_vrr;
      Double I_ERI_D2y_S_Dxy_S_M4_vrr = PAY*I_ERI_Py_S_Dxy_S_M4_vrr+WPY*I_ERI_Py_S_Dxy_S_M5_vrr+oned2z*I_ERI_S_S_Dxy_S_M4_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M5_vrr+oned2k*I_ERI_Py_S_Px_S_M5_vrr;
      Double I_ERI_D2z_S_Dxy_S_M4_vrr = PAZ*I_ERI_Pz_S_Dxy_S_M4_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M5_vrr+oned2z*I_ERI_S_S_Dxy_S_M4_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M5_vrr;
      Double I_ERI_D2x_S_Dxz_S_M4_vrr = PAX*I_ERI_Px_S_Dxz_S_M4_vrr+WPX*I_ERI_Px_S_Dxz_S_M5_vrr+oned2z*I_ERI_S_S_Dxz_S_M4_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M5_vrr+oned2k*I_ERI_Px_S_Pz_S_M5_vrr;
      Double I_ERI_D2y_S_Dxz_S_M4_vrr = PAY*I_ERI_Py_S_Dxz_S_M4_vrr+WPY*I_ERI_Py_S_Dxz_S_M5_vrr+oned2z*I_ERI_S_S_Dxz_S_M4_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M5_vrr;
      Double I_ERI_D2z_S_Dxz_S_M4_vrr = PAZ*I_ERI_Pz_S_Dxz_S_M4_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M5_vrr+oned2z*I_ERI_S_S_Dxz_S_M4_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M5_vrr+oned2k*I_ERI_Pz_S_Px_S_M5_vrr;
      Double I_ERI_D2x_S_D2y_S_M4_vrr = PAX*I_ERI_Px_S_D2y_S_M4_vrr+WPX*I_ERI_Px_S_D2y_S_M5_vrr+oned2z*I_ERI_S_S_D2y_S_M4_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_Dxy_S_D2y_S_M4_vrr = PAY*I_ERI_Px_S_D2y_S_M4_vrr+WPY*I_ERI_Px_S_D2y_S_M5_vrr+2*oned2k*I_ERI_Px_S_Py_S_M5_vrr;
      Double I_ERI_Dxz_S_D2y_S_M4_vrr = PAZ*I_ERI_Px_S_D2y_S_M4_vrr+WPZ*I_ERI_Px_S_D2y_S_M5_vrr;
      Double I_ERI_D2y_S_D2y_S_M4_vrr = PAY*I_ERI_Py_S_D2y_S_M4_vrr+WPY*I_ERI_Py_S_D2y_S_M5_vrr+oned2z*I_ERI_S_S_D2y_S_M4_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M5_vrr+2*oned2k*I_ERI_Py_S_Py_S_M5_vrr;
      Double I_ERI_Dyz_S_D2y_S_M4_vrr = PAZ*I_ERI_Py_S_D2y_S_M4_vrr+WPZ*I_ERI_Py_S_D2y_S_M5_vrr;
      Double I_ERI_D2z_S_D2y_S_M4_vrr = PAZ*I_ERI_Pz_S_D2y_S_M4_vrr+WPZ*I_ERI_Pz_S_D2y_S_M5_vrr+oned2z*I_ERI_S_S_D2y_S_M4_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_D2x_S_Dyz_S_M4_vrr = PAX*I_ERI_Px_S_Dyz_S_M4_vrr+WPX*I_ERI_Px_S_Dyz_S_M5_vrr+oned2z*I_ERI_S_S_Dyz_S_M4_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M5_vrr;
      Double I_ERI_D2y_S_Dyz_S_M4_vrr = PAY*I_ERI_Py_S_Dyz_S_M4_vrr+WPY*I_ERI_Py_S_Dyz_S_M5_vrr+oned2z*I_ERI_S_S_Dyz_S_M4_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M5_vrr+oned2k*I_ERI_Py_S_Pz_S_M5_vrr;
      Double I_ERI_D2z_S_Dyz_S_M4_vrr = PAZ*I_ERI_Pz_S_Dyz_S_M4_vrr+WPZ*I_ERI_Pz_S_Dyz_S_M5_vrr+oned2z*I_ERI_S_S_Dyz_S_M4_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M5_vrr+oned2k*I_ERI_Pz_S_Py_S_M5_vrr;
      Double I_ERI_D2x_S_D2z_S_M4_vrr = PAX*I_ERI_Px_S_D2z_S_M4_vrr+WPX*I_ERI_Px_S_D2z_S_M5_vrr+oned2z*I_ERI_S_S_D2z_S_M4_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_Dxy_S_D2z_S_M4_vrr = PAY*I_ERI_Px_S_D2z_S_M4_vrr+WPY*I_ERI_Px_S_D2z_S_M5_vrr;
      Double I_ERI_Dxz_S_D2z_S_M4_vrr = PAZ*I_ERI_Px_S_D2z_S_M4_vrr+WPZ*I_ERI_Px_S_D2z_S_M5_vrr+2*oned2k*I_ERI_Px_S_Pz_S_M5_vrr;
      Double I_ERI_D2y_S_D2z_S_M4_vrr = PAY*I_ERI_Py_S_D2z_S_M4_vrr+WPY*I_ERI_Py_S_D2z_S_M5_vrr+oned2z*I_ERI_S_S_D2z_S_M4_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_Dyz_S_D2z_S_M4_vrr = PAZ*I_ERI_Py_S_D2z_S_M4_vrr+WPZ*I_ERI_Py_S_D2z_S_M5_vrr+2*oned2k*I_ERI_Py_S_Pz_S_M5_vrr;
      Double I_ERI_D2z_S_D2z_S_M4_vrr = PAZ*I_ERI_Pz_S_D2z_S_M4_vrr+WPZ*I_ERI_Pz_S_D2z_S_M5_vrr+oned2z*I_ERI_S_S_D2z_S_M4_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M5_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M4_vrr = QCX*I_ERI_Px_S_D2x_S_M4_vrr+WQX*I_ERI_Px_S_D2x_S_M5_vrr+2*oned2e*I_ERI_Px_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_Px_S_Px_S_M5_vrr+oned2k*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_Py_S_F3x_S_M4_vrr = QCX*I_ERI_Py_S_D2x_S_M4_vrr+WQX*I_ERI_Py_S_D2x_S_M5_vrr+2*oned2e*I_ERI_Py_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_Py_S_Px_S_M5_vrr;
      Double I_ERI_Pz_S_F3x_S_M4_vrr = QCX*I_ERI_Pz_S_D2x_S_M4_vrr+WQX*I_ERI_Pz_S_D2x_S_M5_vrr+2*oned2e*I_ERI_Pz_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_Pz_S_Px_S_M5_vrr;
      Double I_ERI_Px_S_F2xy_S_M4_vrr = QCY*I_ERI_Px_S_D2x_S_M4_vrr+WQY*I_ERI_Px_S_D2x_S_M5_vrr;
      Double I_ERI_Py_S_F2xy_S_M4_vrr = QCY*I_ERI_Py_S_D2x_S_M4_vrr+WQY*I_ERI_Py_S_D2x_S_M5_vrr+oned2k*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_Pz_S_F2xy_S_M4_vrr = QCY*I_ERI_Pz_S_D2x_S_M4_vrr+WQY*I_ERI_Pz_S_D2x_S_M5_vrr;
      Double I_ERI_Px_S_F2xz_S_M4_vrr = QCZ*I_ERI_Px_S_D2x_S_M4_vrr+WQZ*I_ERI_Px_S_D2x_S_M5_vrr;
      Double I_ERI_Py_S_F2xz_S_M4_vrr = QCZ*I_ERI_Py_S_D2x_S_M4_vrr+WQZ*I_ERI_Py_S_D2x_S_M5_vrr;
      Double I_ERI_Pz_S_F2xz_S_M4_vrr = QCZ*I_ERI_Pz_S_D2x_S_M4_vrr+WQZ*I_ERI_Pz_S_D2x_S_M5_vrr+oned2k*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_Px_S_Fx2y_S_M4_vrr = QCX*I_ERI_Px_S_D2y_S_M4_vrr+WQX*I_ERI_Px_S_D2y_S_M5_vrr+oned2k*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_Py_S_Fx2y_S_M4_vrr = QCX*I_ERI_Py_S_D2y_S_M4_vrr+WQX*I_ERI_Py_S_D2y_S_M5_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M4_vrr = QCX*I_ERI_Pz_S_D2y_S_M4_vrr+WQX*I_ERI_Pz_S_D2y_S_M5_vrr;
      Double I_ERI_Px_S_Fxyz_S_M4_vrr = QCZ*I_ERI_Px_S_Dxy_S_M4_vrr+WQZ*I_ERI_Px_S_Dxy_S_M5_vrr;
      Double I_ERI_Py_S_Fxyz_S_M4_vrr = QCZ*I_ERI_Py_S_Dxy_S_M4_vrr+WQZ*I_ERI_Py_S_Dxy_S_M5_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M4_vrr = QCZ*I_ERI_Pz_S_Dxy_S_M4_vrr+WQZ*I_ERI_Pz_S_Dxy_S_M5_vrr+oned2k*I_ERI_S_S_Dxy_S_M5_vrr;
      Double I_ERI_Px_S_Fx2z_S_M4_vrr = QCX*I_ERI_Px_S_D2z_S_M4_vrr+WQX*I_ERI_Px_S_D2z_S_M5_vrr+oned2k*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_Py_S_Fx2z_S_M4_vrr = QCX*I_ERI_Py_S_D2z_S_M4_vrr+WQX*I_ERI_Py_S_D2z_S_M5_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M4_vrr = QCX*I_ERI_Pz_S_D2z_S_M4_vrr+WQX*I_ERI_Pz_S_D2z_S_M5_vrr;
      Double I_ERI_Px_S_F3y_S_M4_vrr = QCY*I_ERI_Px_S_D2y_S_M4_vrr+WQY*I_ERI_Px_S_D2y_S_M5_vrr+2*oned2e*I_ERI_Px_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_Px_S_Py_S_M5_vrr;
      Double I_ERI_Py_S_F3y_S_M4_vrr = QCY*I_ERI_Py_S_D2y_S_M4_vrr+WQY*I_ERI_Py_S_D2y_S_M5_vrr+2*oned2e*I_ERI_Py_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_Py_S_Py_S_M5_vrr+oned2k*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_Pz_S_F3y_S_M4_vrr = QCY*I_ERI_Pz_S_D2y_S_M4_vrr+WQY*I_ERI_Pz_S_D2y_S_M5_vrr+2*oned2e*I_ERI_Pz_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_Pz_S_Py_S_M5_vrr;
      Double I_ERI_Px_S_F2yz_S_M4_vrr = QCZ*I_ERI_Px_S_D2y_S_M4_vrr+WQZ*I_ERI_Px_S_D2y_S_M5_vrr;
      Double I_ERI_Py_S_F2yz_S_M4_vrr = QCZ*I_ERI_Py_S_D2y_S_M4_vrr+WQZ*I_ERI_Py_S_D2y_S_M5_vrr;
      Double I_ERI_Pz_S_F2yz_S_M4_vrr = QCZ*I_ERI_Pz_S_D2y_S_M4_vrr+WQZ*I_ERI_Pz_S_D2y_S_M5_vrr+oned2k*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_Px_S_Fy2z_S_M4_vrr = QCY*I_ERI_Px_S_D2z_S_M4_vrr+WQY*I_ERI_Px_S_D2z_S_M5_vrr;
      Double I_ERI_Py_S_Fy2z_S_M4_vrr = QCY*I_ERI_Py_S_D2z_S_M4_vrr+WQY*I_ERI_Py_S_D2z_S_M5_vrr+oned2k*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M4_vrr = QCY*I_ERI_Pz_S_D2z_S_M4_vrr+WQY*I_ERI_Pz_S_D2z_S_M5_vrr;
      Double I_ERI_Px_S_F3z_S_M4_vrr = QCZ*I_ERI_Px_S_D2z_S_M4_vrr+WQZ*I_ERI_Px_S_D2z_S_M5_vrr+2*oned2e*I_ERI_Px_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_Px_S_Pz_S_M5_vrr;
      Double I_ERI_Py_S_F3z_S_M4_vrr = QCZ*I_ERI_Py_S_D2z_S_M4_vrr+WQZ*I_ERI_Py_S_D2z_S_M5_vrr+2*oned2e*I_ERI_Py_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_Py_S_Pz_S_M5_vrr;
      Double I_ERI_Pz_S_F3z_S_M4_vrr = QCZ*I_ERI_Pz_S_D2z_S_M4_vrr+WQZ*I_ERI_Pz_S_D2z_S_M5_vrr+2*oned2e*I_ERI_Pz_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_Pz_S_Pz_S_M5_vrr+oned2k*I_ERI_S_S_D2z_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 27 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M5
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M4_vrr = QCX*I_ERI_D2x_S_D2x_S_M4_vrr+WQX*I_ERI_D2x_S_D2x_S_M5_vrr+2*oned2e*I_ERI_D2x_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_D2x_S_Px_S_M5_vrr+2*oned2k*I_ERI_Px_S_D2x_S_M5_vrr;
      Double I_ERI_Dxy_S_F3x_S_M4_vrr = QCX*I_ERI_Dxy_S_D2x_S_M4_vrr+WQX*I_ERI_Dxy_S_D2x_S_M5_vrr+2*oned2e*I_ERI_Dxy_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_Dxy_S_Px_S_M5_vrr+oned2k*I_ERI_Py_S_D2x_S_M5_vrr;
      Double I_ERI_D2y_S_F3x_S_M4_vrr = QCX*I_ERI_D2y_S_D2x_S_M4_vrr+WQX*I_ERI_D2y_S_D2x_S_M5_vrr+2*oned2e*I_ERI_D2y_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_D2y_S_Px_S_M5_vrr;
      Double I_ERI_D2z_S_F3x_S_M4_vrr = QCX*I_ERI_D2z_S_D2x_S_M4_vrr+WQX*I_ERI_D2z_S_D2x_S_M5_vrr+2*oned2e*I_ERI_D2z_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_D2z_S_Px_S_M5_vrr;
      Double I_ERI_D2x_S_F2xy_S_M4_vrr = QCY*I_ERI_D2x_S_D2x_S_M4_vrr+WQY*I_ERI_D2x_S_D2x_S_M5_vrr;
      Double I_ERI_D2y_S_F2xy_S_M4_vrr = QCY*I_ERI_D2y_S_D2x_S_M4_vrr+WQY*I_ERI_D2y_S_D2x_S_M5_vrr+2*oned2k*I_ERI_Py_S_D2x_S_M5_vrr;
      Double I_ERI_D2z_S_F2xy_S_M4_vrr = QCY*I_ERI_D2z_S_D2x_S_M4_vrr+WQY*I_ERI_D2z_S_D2x_S_M5_vrr;
      Double I_ERI_D2x_S_F2xz_S_M4_vrr = QCZ*I_ERI_D2x_S_D2x_S_M4_vrr+WQZ*I_ERI_D2x_S_D2x_S_M5_vrr;
      Double I_ERI_D2y_S_F2xz_S_M4_vrr = QCZ*I_ERI_D2y_S_D2x_S_M4_vrr+WQZ*I_ERI_D2y_S_D2x_S_M5_vrr;
      Double I_ERI_D2z_S_F2xz_S_M4_vrr = QCZ*I_ERI_D2z_S_D2x_S_M4_vrr+WQZ*I_ERI_D2z_S_D2x_S_M5_vrr+2*oned2k*I_ERI_Pz_S_D2x_S_M5_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M4_vrr = QCX*I_ERI_D2x_S_D2y_S_M4_vrr+WQX*I_ERI_D2x_S_D2y_S_M5_vrr+2*oned2k*I_ERI_Px_S_D2y_S_M5_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M4_vrr = QCX*I_ERI_D2y_S_D2y_S_M4_vrr+WQX*I_ERI_D2y_S_D2y_S_M5_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M4_vrr = QCX*I_ERI_D2z_S_D2y_S_M4_vrr+WQX*I_ERI_D2z_S_D2y_S_M5_vrr;
      Double I_ERI_D2x_S_Fxyz_S_M4_vrr = QCZ*I_ERI_D2x_S_Dxy_S_M4_vrr+WQZ*I_ERI_D2x_S_Dxy_S_M5_vrr;
      Double I_ERI_D2y_S_Fxyz_S_M4_vrr = QCZ*I_ERI_D2y_S_Dxy_S_M4_vrr+WQZ*I_ERI_D2y_S_Dxy_S_M5_vrr;
      Double I_ERI_D2z_S_Fxyz_S_M4_vrr = QCZ*I_ERI_D2z_S_Dxy_S_M4_vrr+WQZ*I_ERI_D2z_S_Dxy_S_M5_vrr+2*oned2k*I_ERI_Pz_S_Dxy_S_M5_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M4_vrr = QCX*I_ERI_D2x_S_D2z_S_M4_vrr+WQX*I_ERI_D2x_S_D2z_S_M5_vrr+2*oned2k*I_ERI_Px_S_D2z_S_M5_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M4_vrr = QCX*I_ERI_D2y_S_D2z_S_M4_vrr+WQX*I_ERI_D2y_S_D2z_S_M5_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M4_vrr = QCX*I_ERI_D2z_S_D2z_S_M4_vrr+WQX*I_ERI_D2z_S_D2z_S_M5_vrr;
      Double I_ERI_D2x_S_F3y_S_M4_vrr = QCY*I_ERI_D2x_S_D2y_S_M4_vrr+WQY*I_ERI_D2x_S_D2y_S_M5_vrr+2*oned2e*I_ERI_D2x_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_D2x_S_Py_S_M5_vrr;
      Double I_ERI_Dxy_S_F3y_S_M4_vrr = QCY*I_ERI_Dxy_S_D2y_S_M4_vrr+WQY*I_ERI_Dxy_S_D2y_S_M5_vrr+2*oned2e*I_ERI_Dxy_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_Dxy_S_Py_S_M5_vrr+oned2k*I_ERI_Px_S_D2y_S_M5_vrr;
      Double I_ERI_D2y_S_F3y_S_M4_vrr = QCY*I_ERI_D2y_S_D2y_S_M4_vrr+WQY*I_ERI_D2y_S_D2y_S_M5_vrr+2*oned2e*I_ERI_D2y_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_D2y_S_Py_S_M5_vrr+2*oned2k*I_ERI_Py_S_D2y_S_M5_vrr;
      Double I_ERI_D2z_S_F3y_S_M4_vrr = QCY*I_ERI_D2z_S_D2y_S_M4_vrr+WQY*I_ERI_D2z_S_D2y_S_M5_vrr+2*oned2e*I_ERI_D2z_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_D2z_S_Py_S_M5_vrr;
      Double I_ERI_D2x_S_F2yz_S_M4_vrr = QCZ*I_ERI_D2x_S_D2y_S_M4_vrr+WQZ*I_ERI_D2x_S_D2y_S_M5_vrr;
      Double I_ERI_D2y_S_F2yz_S_M4_vrr = QCZ*I_ERI_D2y_S_D2y_S_M4_vrr+WQZ*I_ERI_D2y_S_D2y_S_M5_vrr;
      Double I_ERI_D2z_S_F2yz_S_M4_vrr = QCZ*I_ERI_D2z_S_D2y_S_M4_vrr+WQZ*I_ERI_D2z_S_D2y_S_M5_vrr+2*oned2k*I_ERI_Pz_S_D2y_S_M5_vrr;
      Double I_ERI_D2x_S_Fy2z_S_M4_vrr = QCY*I_ERI_D2x_S_D2z_S_M4_vrr+WQY*I_ERI_D2x_S_D2z_S_M5_vrr;
      Double I_ERI_D2y_S_Fy2z_S_M4_vrr = QCY*I_ERI_D2y_S_D2z_S_M4_vrr+WQY*I_ERI_D2y_S_D2z_S_M5_vrr+2*oned2k*I_ERI_Py_S_D2z_S_M5_vrr;
      Double I_ERI_D2z_S_Fy2z_S_M4_vrr = QCY*I_ERI_D2z_S_D2z_S_M4_vrr+WQY*I_ERI_D2z_S_D2z_S_M5_vrr;
      Double I_ERI_D2x_S_F3z_S_M4_vrr = QCZ*I_ERI_D2x_S_D2z_S_M4_vrr+WQZ*I_ERI_D2x_S_D2z_S_M5_vrr+2*oned2e*I_ERI_D2x_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_D2x_S_Pz_S_M5_vrr;
      Double I_ERI_Dxy_S_F3z_S_M4_vrr = QCZ*I_ERI_Dxy_S_D2z_S_M4_vrr+WQZ*I_ERI_Dxy_S_D2z_S_M5_vrr+2*oned2e*I_ERI_Dxy_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_Dxy_S_Pz_S_M5_vrr;
      Double I_ERI_D2y_S_F3z_S_M4_vrr = QCZ*I_ERI_D2y_S_D2z_S_M4_vrr+WQZ*I_ERI_D2y_S_D2z_S_M5_vrr+2*oned2e*I_ERI_D2y_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_D2y_S_Pz_S_M5_vrr;
      Double I_ERI_D2z_S_F3z_S_M4_vrr = QCZ*I_ERI_D2z_S_D2z_S_M4_vrr+WQZ*I_ERI_D2z_S_D2z_S_M5_vrr+2*oned2e*I_ERI_D2z_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_D2z_S_Pz_S_M5_vrr+2*oned2k*I_ERI_Pz_S_D2z_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 48 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M5
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_M4_vrr = PAX*I_ERI_D2x_S_D2x_S_M4_vrr+WPX*I_ERI_D2x_S_D2x_S_M5_vrr+2*oned2z*I_ERI_Px_S_D2x_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_D2x_S_M5_vrr+2*oned2k*I_ERI_D2x_S_Px_S_M5_vrr;
      Double I_ERI_F2xy_S_D2x_S_M4_vrr = PAY*I_ERI_D2x_S_D2x_S_M4_vrr+WPY*I_ERI_D2x_S_D2x_S_M5_vrr;
      Double I_ERI_F3y_S_D2x_S_M4_vrr = PAY*I_ERI_D2y_S_D2x_S_M4_vrr+WPY*I_ERI_D2y_S_D2x_S_M5_vrr+2*oned2z*I_ERI_Py_S_D2x_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_D2x_S_M5_vrr;
      Double I_ERI_F3z_S_D2x_S_M4_vrr = PAZ*I_ERI_D2z_S_D2x_S_M4_vrr+WPZ*I_ERI_D2z_S_D2x_S_M5_vrr+2*oned2z*I_ERI_Pz_S_D2x_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_D2x_S_M5_vrr;
      Double I_ERI_F3x_S_D2y_S_M4_vrr = PAX*I_ERI_D2x_S_D2y_S_M4_vrr+WPX*I_ERI_D2x_S_D2y_S_M5_vrr+2*oned2z*I_ERI_Px_S_D2y_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_D2y_S_M5_vrr;
      Double I_ERI_F2xy_S_D2y_S_M4_vrr = PAY*I_ERI_D2x_S_D2y_S_M4_vrr+WPY*I_ERI_D2x_S_D2y_S_M5_vrr+2*oned2k*I_ERI_D2x_S_Py_S_M5_vrr;
      Double I_ERI_F3y_S_D2y_S_M4_vrr = PAY*I_ERI_D2y_S_D2y_S_M4_vrr+WPY*I_ERI_D2y_S_D2y_S_M5_vrr+2*oned2z*I_ERI_Py_S_D2y_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_D2y_S_M5_vrr+2*oned2k*I_ERI_D2y_S_Py_S_M5_vrr;
      Double I_ERI_F3z_S_D2y_S_M4_vrr = PAZ*I_ERI_D2z_S_D2y_S_M4_vrr+WPZ*I_ERI_D2z_S_D2y_S_M5_vrr+2*oned2z*I_ERI_Pz_S_D2y_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_D2y_S_M5_vrr;
      Double I_ERI_F3x_S_D2z_S_M4_vrr = PAX*I_ERI_D2x_S_D2z_S_M4_vrr+WPX*I_ERI_D2x_S_D2z_S_M5_vrr+2*oned2z*I_ERI_Px_S_D2z_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_D2z_S_M5_vrr;
      Double I_ERI_F2xy_S_D2z_S_M4_vrr = PAY*I_ERI_D2x_S_D2z_S_M4_vrr+WPY*I_ERI_D2x_S_D2z_S_M5_vrr;
      Double I_ERI_F3y_S_D2z_S_M4_vrr = PAY*I_ERI_D2y_S_D2z_S_M4_vrr+WPY*I_ERI_D2y_S_D2z_S_M5_vrr+2*oned2z*I_ERI_Py_S_D2z_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_D2z_S_M5_vrr;
      Double I_ERI_F3z_S_D2z_S_M4_vrr = PAZ*I_ERI_D2z_S_D2z_S_M4_vrr+WPZ*I_ERI_D2z_S_D2z_S_M5_vrr+2*oned2z*I_ERI_Pz_S_D2z_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_D2z_S_M5_vrr+2*oned2k*I_ERI_D2z_S_Pz_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 68 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M5
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M5
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_M4_vrr = PAX*I_ERI_D2x_S_F3x_S_M4_vrr+WPX*I_ERI_D2x_S_F3x_S_M5_vrr+2*oned2z*I_ERI_Px_S_F3x_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_F3x_S_M5_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M5_vrr;
      Double I_ERI_F2xy_S_F3x_S_M4_vrr = PAY*I_ERI_D2x_S_F3x_S_M4_vrr+WPY*I_ERI_D2x_S_F3x_S_M5_vrr;
      Double I_ERI_F3y_S_F3x_S_M4_vrr = PAY*I_ERI_D2y_S_F3x_S_M4_vrr+WPY*I_ERI_D2y_S_F3x_S_M5_vrr+2*oned2z*I_ERI_Py_S_F3x_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_F3x_S_M5_vrr;
      Double I_ERI_F3z_S_F3x_S_M4_vrr = PAZ*I_ERI_D2z_S_F3x_S_M4_vrr+WPZ*I_ERI_D2z_S_F3x_S_M5_vrr+2*oned2z*I_ERI_Pz_S_F3x_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_F3x_S_M5_vrr;
      Double I_ERI_F3x_S_F2xy_S_M4_vrr = PAX*I_ERI_D2x_S_F2xy_S_M4_vrr+WPX*I_ERI_D2x_S_F2xy_S_M5_vrr+2*oned2z*I_ERI_Px_S_F2xy_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_F2xy_S_M5_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M5_vrr;
      Double I_ERI_F2xy_S_F2xy_S_M4_vrr = PAY*I_ERI_D2x_S_F2xy_S_M4_vrr+WPY*I_ERI_D2x_S_F2xy_S_M5_vrr+oned2k*I_ERI_D2x_S_D2x_S_M5_vrr;
      Double I_ERI_F3y_S_F2xy_S_M4_vrr = PAY*I_ERI_D2y_S_F2xy_S_M4_vrr+WPY*I_ERI_D2y_S_F2xy_S_M5_vrr+2*oned2z*I_ERI_Py_S_F2xy_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_F2xy_S_M5_vrr+oned2k*I_ERI_D2y_S_D2x_S_M5_vrr;
      Double I_ERI_F3z_S_F2xy_S_M4_vrr = PAZ*I_ERI_D2z_S_F2xy_S_M4_vrr+WPZ*I_ERI_D2z_S_F2xy_S_M5_vrr+2*oned2z*I_ERI_Pz_S_F2xy_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xy_S_M5_vrr;
      Double I_ERI_F3x_S_F2xz_S_M4_vrr = PAX*I_ERI_D2x_S_F2xz_S_M4_vrr+WPX*I_ERI_D2x_S_F2xz_S_M5_vrr+2*oned2z*I_ERI_Px_S_F2xz_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_F2xz_S_M5_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M5_vrr;
      Double I_ERI_F2xy_S_F2xz_S_M4_vrr = PAY*I_ERI_D2x_S_F2xz_S_M4_vrr+WPY*I_ERI_D2x_S_F2xz_S_M5_vrr;
      Double I_ERI_F3y_S_F2xz_S_M4_vrr = PAY*I_ERI_D2y_S_F2xz_S_M4_vrr+WPY*I_ERI_D2y_S_F2xz_S_M5_vrr+2*oned2z*I_ERI_Py_S_F2xz_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_F2xz_S_M5_vrr;
      Double I_ERI_F3z_S_F2xz_S_M4_vrr = PAZ*I_ERI_D2z_S_F2xz_S_M4_vrr+WPZ*I_ERI_D2z_S_F2xz_S_M5_vrr+2*oned2z*I_ERI_Pz_S_F2xz_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xz_S_M5_vrr+oned2k*I_ERI_D2z_S_D2x_S_M5_vrr;
      Double I_ERI_F3x_S_Fx2y_S_M4_vrr = PAX*I_ERI_D2x_S_Fx2y_S_M4_vrr+WPX*I_ERI_D2x_S_Fx2y_S_M5_vrr+2*oned2z*I_ERI_Px_S_Fx2y_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2y_S_M5_vrr+oned2k*I_ERI_D2x_S_D2y_S_M5_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_M4_vrr = PAY*I_ERI_D2x_S_Fx2y_S_M4_vrr+WPY*I_ERI_D2x_S_Fx2y_S_M5_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M5_vrr;
      Double I_ERI_F3y_S_Fx2y_S_M4_vrr = PAY*I_ERI_D2y_S_Fx2y_S_M4_vrr+WPY*I_ERI_D2y_S_Fx2y_S_M5_vrr+2*oned2z*I_ERI_Py_S_Fx2y_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2y_S_M5_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M5_vrr;
      Double I_ERI_F3z_S_Fx2y_S_M4_vrr = PAZ*I_ERI_D2z_S_Fx2y_S_M4_vrr+WPZ*I_ERI_D2z_S_Fx2y_S_M5_vrr+2*oned2z*I_ERI_Pz_S_Fx2y_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2y_S_M5_vrr;
      Double I_ERI_F3x_S_Fx2z_S_M4_vrr = PAX*I_ERI_D2x_S_Fx2z_S_M4_vrr+WPX*I_ERI_D2x_S_Fx2z_S_M5_vrr+2*oned2z*I_ERI_Px_S_Fx2z_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2z_S_M5_vrr+oned2k*I_ERI_D2x_S_D2z_S_M5_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_M4_vrr = PAY*I_ERI_D2x_S_Fx2z_S_M4_vrr+WPY*I_ERI_D2x_S_Fx2z_S_M5_vrr;
      Double I_ERI_F3y_S_Fx2z_S_M4_vrr = PAY*I_ERI_D2y_S_Fx2z_S_M4_vrr+WPY*I_ERI_D2y_S_Fx2z_S_M5_vrr+2*oned2z*I_ERI_Py_S_Fx2z_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2z_S_M5_vrr;
      Double I_ERI_F3z_S_Fx2z_S_M4_vrr = PAZ*I_ERI_D2z_S_Fx2z_S_M4_vrr+WPZ*I_ERI_D2z_S_Fx2z_S_M5_vrr+2*oned2z*I_ERI_Pz_S_Fx2z_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2z_S_M5_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M5_vrr;
      Double I_ERI_F3x_S_F3y_S_M4_vrr = PAX*I_ERI_D2x_S_F3y_S_M4_vrr+WPX*I_ERI_D2x_S_F3y_S_M5_vrr+2*oned2z*I_ERI_Px_S_F3y_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_F3y_S_M5_vrr;
      Double I_ERI_F2xy_S_F3y_S_M4_vrr = PAY*I_ERI_D2x_S_F3y_S_M4_vrr+WPY*I_ERI_D2x_S_F3y_S_M5_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M5_vrr;
      Double I_ERI_F3y_S_F3y_S_M4_vrr = PAY*I_ERI_D2y_S_F3y_S_M4_vrr+WPY*I_ERI_D2y_S_F3y_S_M5_vrr+2*oned2z*I_ERI_Py_S_F3y_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_F3y_S_M5_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M5_vrr;
      Double I_ERI_F3z_S_F3y_S_M4_vrr = PAZ*I_ERI_D2z_S_F3y_S_M4_vrr+WPZ*I_ERI_D2z_S_F3y_S_M5_vrr+2*oned2z*I_ERI_Pz_S_F3y_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_F3y_S_M5_vrr;
      Double I_ERI_F3x_S_F2yz_S_M4_vrr = PAX*I_ERI_D2x_S_F2yz_S_M4_vrr+WPX*I_ERI_D2x_S_F2yz_S_M5_vrr+2*oned2z*I_ERI_Px_S_F2yz_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_F2yz_S_M5_vrr;
      Double I_ERI_F2xy_S_F2yz_S_M4_vrr = PAY*I_ERI_D2x_S_F2yz_S_M4_vrr+WPY*I_ERI_D2x_S_F2yz_S_M5_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M5_vrr;
      Double I_ERI_F3y_S_F2yz_S_M4_vrr = PAY*I_ERI_D2y_S_F2yz_S_M4_vrr+WPY*I_ERI_D2y_S_F2yz_S_M5_vrr+2*oned2z*I_ERI_Py_S_F2yz_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_F2yz_S_M5_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M5_vrr;
      Double I_ERI_F3z_S_F2yz_S_M4_vrr = PAZ*I_ERI_D2z_S_F2yz_S_M4_vrr+WPZ*I_ERI_D2z_S_F2yz_S_M5_vrr+2*oned2z*I_ERI_Pz_S_F2yz_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_F2yz_S_M5_vrr+oned2k*I_ERI_D2z_S_D2y_S_M5_vrr;
      Double I_ERI_F3x_S_F3z_S_M4_vrr = PAX*I_ERI_D2x_S_F3z_S_M4_vrr+WPX*I_ERI_D2x_S_F3z_S_M5_vrr+2*oned2z*I_ERI_Px_S_F3z_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_F3z_S_M5_vrr;
      Double I_ERI_F2xy_S_F3z_S_M4_vrr = PAY*I_ERI_D2x_S_F3z_S_M4_vrr+WPY*I_ERI_D2x_S_F3z_S_M5_vrr;
      Double I_ERI_F3y_S_F3z_S_M4_vrr = PAY*I_ERI_D2y_S_F3z_S_M4_vrr+WPY*I_ERI_D2y_S_F3z_S_M5_vrr+2*oned2z*I_ERI_Py_S_F3z_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_F3z_S_M5_vrr;
      Double I_ERI_F3z_S_F3z_S_M4_vrr = PAZ*I_ERI_D2z_S_F3z_S_M4_vrr+WPZ*I_ERI_D2z_S_F3z_S_M5_vrr+2*oned2z*I_ERI_Pz_S_F3z_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_F3z_S_M5_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_P_S_P_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_Px_S_Px_S_M3_vrr = PAX*I_ERI_S_S_Px_S_M3_vrr+WPX*I_ERI_S_S_Px_S_M4_vrr+oned2k*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Py_S_Px_S_M3_vrr = PAY*I_ERI_S_S_Px_S_M3_vrr+WPY*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Pz_S_Px_S_M3_vrr = PAZ*I_ERI_S_S_Px_S_M3_vrr+WPZ*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Px_S_Py_S_M3_vrr = PAX*I_ERI_S_S_Py_S_M3_vrr+WPX*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Py_S_Py_S_M3_vrr = PAY*I_ERI_S_S_Py_S_M3_vrr+WPY*I_ERI_S_S_Py_S_M4_vrr+oned2k*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Pz_S_Py_S_M3_vrr = PAZ*I_ERI_S_S_Py_S_M3_vrr+WPZ*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Px_S_Pz_S_M3_vrr = PAX*I_ERI_S_S_Pz_S_M3_vrr+WPX*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Py_S_Pz_S_M3_vrr = PAY*I_ERI_S_S_Pz_S_M3_vrr+WPY*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Pz_S_Pz_S_M3_vrr = PAZ*I_ERI_S_S_Pz_S_M3_vrr+WPZ*I_ERI_S_S_Pz_S_M4_vrr+oned2k*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_D_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M3_vrr = QCX*I_ERI_S_S_Px_S_M3_vrr+WQX*I_ERI_S_S_Px_S_M4_vrr+oned2e*I_ERI_S_S_S_S_M3_vrr-rhod2esq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Dxy_S_M3_vrr = QCY*I_ERI_S_S_Px_S_M3_vrr+WQY*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_S_S_Dxz_S_M3_vrr = QCZ*I_ERI_S_S_Px_S_M3_vrr+WQZ*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_S_S_D2y_S_M3_vrr = QCY*I_ERI_S_S_Py_S_M3_vrr+WQY*I_ERI_S_S_Py_S_M4_vrr+oned2e*I_ERI_S_S_S_S_M3_vrr-rhod2esq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_S_S_Dyz_S_M3_vrr = QCZ*I_ERI_S_S_Py_S_M3_vrr+WQZ*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_S_S_D2z_S_M3_vrr = QCZ*I_ERI_S_S_Pz_S_M3_vrr+WQZ*I_ERI_S_S_Pz_S_M4_vrr+oned2e*I_ERI_S_S_S_S_M3_vrr-rhod2esq*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M3_vrr = PAX*I_ERI_S_S_D2x_S_M3_vrr+WPX*I_ERI_S_S_D2x_S_M4_vrr+2*oned2k*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Py_S_D2x_S_M3_vrr = PAY*I_ERI_S_S_D2x_S_M3_vrr+WPY*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Pz_S_D2x_S_M3_vrr = PAZ*I_ERI_S_S_D2x_S_M3_vrr+WPZ*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Px_S_Dxy_S_M3_vrr = PAX*I_ERI_S_S_Dxy_S_M3_vrr+WPX*I_ERI_S_S_Dxy_S_M4_vrr+oned2k*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Py_S_Dxy_S_M3_vrr = PAY*I_ERI_S_S_Dxy_S_M3_vrr+WPY*I_ERI_S_S_Dxy_S_M4_vrr+oned2k*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Pz_S_Dxy_S_M3_vrr = PAZ*I_ERI_S_S_Dxy_S_M3_vrr+WPZ*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Px_S_Dxz_S_M3_vrr = PAX*I_ERI_S_S_Dxz_S_M3_vrr+WPX*I_ERI_S_S_Dxz_S_M4_vrr+oned2k*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Py_S_Dxz_S_M3_vrr = PAY*I_ERI_S_S_Dxz_S_M3_vrr+WPY*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Pz_S_Dxz_S_M3_vrr = PAZ*I_ERI_S_S_Dxz_S_M3_vrr+WPZ*I_ERI_S_S_Dxz_S_M4_vrr+oned2k*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Px_S_D2y_S_M3_vrr = PAX*I_ERI_S_S_D2y_S_M3_vrr+WPX*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Py_S_D2y_S_M3_vrr = PAY*I_ERI_S_S_D2y_S_M3_vrr+WPY*I_ERI_S_S_D2y_S_M4_vrr+2*oned2k*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Pz_S_D2y_S_M3_vrr = PAZ*I_ERI_S_S_D2y_S_M3_vrr+WPZ*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Px_S_Dyz_S_M3_vrr = PAX*I_ERI_S_S_Dyz_S_M3_vrr+WPX*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Py_S_Dyz_S_M3_vrr = PAY*I_ERI_S_S_Dyz_S_M3_vrr+WPY*I_ERI_S_S_Dyz_S_M4_vrr+oned2k*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Pz_S_Dyz_S_M3_vrr = PAZ*I_ERI_S_S_Dyz_S_M3_vrr+WPZ*I_ERI_S_S_Dyz_S_M4_vrr+oned2k*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Px_S_D2z_S_M3_vrr = PAX*I_ERI_S_S_D2z_S_M3_vrr+WPX*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Py_S_D2z_S_M3_vrr = PAY*I_ERI_S_S_D2z_S_M3_vrr+WPY*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Pz_S_D2z_S_M3_vrr = PAZ*I_ERI_S_S_D2z_S_M3_vrr+WPZ*I_ERI_S_S_D2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_P_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M3_vrr = PAX*I_ERI_Px_S_Px_S_M3_vrr+WPX*I_ERI_Px_S_Px_S_M4_vrr+oned2z*I_ERI_S_S_Px_S_M3_vrr-rhod2zsq*I_ERI_S_S_Px_S_M4_vrr+oned2k*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_Dxy_S_Px_S_M3_vrr = PAY*I_ERI_Px_S_Px_S_M3_vrr+WPY*I_ERI_Px_S_Px_S_M4_vrr;
      Double I_ERI_Dxz_S_Px_S_M3_vrr = PAZ*I_ERI_Px_S_Px_S_M3_vrr+WPZ*I_ERI_Px_S_Px_S_M4_vrr;
      Double I_ERI_D2y_S_Px_S_M3_vrr = PAY*I_ERI_Py_S_Px_S_M3_vrr+WPY*I_ERI_Py_S_Px_S_M4_vrr+oned2z*I_ERI_S_S_Px_S_M3_vrr-rhod2zsq*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_Dyz_S_Px_S_M3_vrr = PAZ*I_ERI_Py_S_Px_S_M3_vrr+WPZ*I_ERI_Py_S_Px_S_M4_vrr;
      Double I_ERI_D2z_S_Px_S_M3_vrr = PAZ*I_ERI_Pz_S_Px_S_M3_vrr+WPZ*I_ERI_Pz_S_Px_S_M4_vrr+oned2z*I_ERI_S_S_Px_S_M3_vrr-rhod2zsq*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_D2x_S_Py_S_M3_vrr = PAX*I_ERI_Px_S_Py_S_M3_vrr+WPX*I_ERI_Px_S_Py_S_M4_vrr+oned2z*I_ERI_S_S_Py_S_M3_vrr-rhod2zsq*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_Dxy_S_Py_S_M3_vrr = PAY*I_ERI_Px_S_Py_S_M3_vrr+WPY*I_ERI_Px_S_Py_S_M4_vrr+oned2k*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_Dxz_S_Py_S_M3_vrr = PAZ*I_ERI_Px_S_Py_S_M3_vrr+WPZ*I_ERI_Px_S_Py_S_M4_vrr;
      Double I_ERI_D2y_S_Py_S_M3_vrr = PAY*I_ERI_Py_S_Py_S_M3_vrr+WPY*I_ERI_Py_S_Py_S_M4_vrr+oned2z*I_ERI_S_S_Py_S_M3_vrr-rhod2zsq*I_ERI_S_S_Py_S_M4_vrr+oned2k*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_Dyz_S_Py_S_M3_vrr = PAZ*I_ERI_Py_S_Py_S_M3_vrr+WPZ*I_ERI_Py_S_Py_S_M4_vrr;
      Double I_ERI_D2z_S_Py_S_M3_vrr = PAZ*I_ERI_Pz_S_Py_S_M3_vrr+WPZ*I_ERI_Pz_S_Py_S_M4_vrr+oned2z*I_ERI_S_S_Py_S_M3_vrr-rhod2zsq*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_D2x_S_Pz_S_M3_vrr = PAX*I_ERI_Px_S_Pz_S_M3_vrr+WPX*I_ERI_Px_S_Pz_S_M4_vrr+oned2z*I_ERI_S_S_Pz_S_M3_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Dxy_S_Pz_S_M3_vrr = PAY*I_ERI_Px_S_Pz_S_M3_vrr+WPY*I_ERI_Px_S_Pz_S_M4_vrr;
      Double I_ERI_Dxz_S_Pz_S_M3_vrr = PAZ*I_ERI_Px_S_Pz_S_M3_vrr+WPZ*I_ERI_Px_S_Pz_S_M4_vrr+oned2k*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_D2y_S_Pz_S_M3_vrr = PAY*I_ERI_Py_S_Pz_S_M3_vrr+WPY*I_ERI_Py_S_Pz_S_M4_vrr+oned2z*I_ERI_S_S_Pz_S_M3_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M4_vrr;
      Double I_ERI_Dyz_S_Pz_S_M3_vrr = PAZ*I_ERI_Py_S_Pz_S_M3_vrr+WPZ*I_ERI_Py_S_Pz_S_M4_vrr+oned2k*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_D2z_S_Pz_S_M3_vrr = PAZ*I_ERI_Pz_S_Pz_S_M3_vrr+WPZ*I_ERI_Pz_S_Pz_S_M4_vrr+oned2z*I_ERI_S_S_Pz_S_M3_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M4_vrr+oned2k*I_ERI_Pz_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M3_vrr = QCX*I_ERI_S_S_D2x_S_M3_vrr+WQX*I_ERI_S_S_D2x_S_M4_vrr+2*oned2e*I_ERI_S_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M4_vrr;
      Double I_ERI_S_S_F2xy_S_M3_vrr = QCY*I_ERI_S_S_D2x_S_M3_vrr+WQY*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_F2xz_S_M3_vrr = QCZ*I_ERI_S_S_D2x_S_M3_vrr+WQZ*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_Fx2y_S_M3_vrr = QCX*I_ERI_S_S_D2y_S_M3_vrr+WQX*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_Fxyz_S_M3_vrr = QCZ*I_ERI_S_S_Dxy_S_M3_vrr+WQZ*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_S_S_Fx2z_S_M3_vrr = QCX*I_ERI_S_S_D2z_S_M3_vrr+WQX*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_S_S_F3y_S_M3_vrr = QCY*I_ERI_S_S_D2y_S_M3_vrr+WQY*I_ERI_S_S_D2y_S_M4_vrr+2*oned2e*I_ERI_S_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M4_vrr;
      Double I_ERI_S_S_F2yz_S_M3_vrr = QCZ*I_ERI_S_S_D2y_S_M3_vrr+WQZ*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_Fy2z_S_M3_vrr = QCY*I_ERI_S_S_D2z_S_M3_vrr+WQY*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_S_S_F3z_S_M3_vrr = QCZ*I_ERI_S_S_D2z_S_M3_vrr+WQZ*I_ERI_S_S_D2z_S_M4_vrr+2*oned2e*I_ERI_S_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M3_vrr = PAX*I_ERI_Px_S_D2x_S_M3_vrr+WPX*I_ERI_Px_S_D2x_S_M4_vrr+oned2z*I_ERI_S_S_D2x_S_M3_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M4_vrr+2*oned2k*I_ERI_Px_S_Px_S_M4_vrr;
      Double I_ERI_Dxy_S_D2x_S_M3_vrr = PAY*I_ERI_Px_S_D2x_S_M3_vrr+WPY*I_ERI_Px_S_D2x_S_M4_vrr;
      Double I_ERI_Dxz_S_D2x_S_M3_vrr = PAZ*I_ERI_Px_S_D2x_S_M3_vrr+WPZ*I_ERI_Px_S_D2x_S_M4_vrr;
      Double I_ERI_D2y_S_D2x_S_M3_vrr = PAY*I_ERI_Py_S_D2x_S_M3_vrr+WPY*I_ERI_Py_S_D2x_S_M4_vrr+oned2z*I_ERI_S_S_D2x_S_M3_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Dyz_S_D2x_S_M3_vrr = PAZ*I_ERI_Py_S_D2x_S_M3_vrr+WPZ*I_ERI_Py_S_D2x_S_M4_vrr;
      Double I_ERI_D2z_S_D2x_S_M3_vrr = PAZ*I_ERI_Pz_S_D2x_S_M3_vrr+WPZ*I_ERI_Pz_S_D2x_S_M4_vrr+oned2z*I_ERI_S_S_D2x_S_M3_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_D2x_S_Dxy_S_M3_vrr = PAX*I_ERI_Px_S_Dxy_S_M3_vrr+WPX*I_ERI_Px_S_Dxy_S_M4_vrr+oned2z*I_ERI_S_S_Dxy_S_M3_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M4_vrr+oned2k*I_ERI_Px_S_Py_S_M4_vrr;
      Double I_ERI_D2y_S_Dxy_S_M3_vrr = PAY*I_ERI_Py_S_Dxy_S_M3_vrr+WPY*I_ERI_Py_S_Dxy_S_M4_vrr+oned2z*I_ERI_S_S_Dxy_S_M3_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M4_vrr+oned2k*I_ERI_Py_S_Px_S_M4_vrr;
      Double I_ERI_D2z_S_Dxy_S_M3_vrr = PAZ*I_ERI_Pz_S_Dxy_S_M3_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M4_vrr+oned2z*I_ERI_S_S_Dxy_S_M3_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_D2x_S_Dxz_S_M3_vrr = PAX*I_ERI_Px_S_Dxz_S_M3_vrr+WPX*I_ERI_Px_S_Dxz_S_M4_vrr+oned2z*I_ERI_S_S_Dxz_S_M3_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M4_vrr+oned2k*I_ERI_Px_S_Pz_S_M4_vrr;
      Double I_ERI_D2y_S_Dxz_S_M3_vrr = PAY*I_ERI_Py_S_Dxz_S_M3_vrr+WPY*I_ERI_Py_S_Dxz_S_M4_vrr+oned2z*I_ERI_S_S_Dxz_S_M3_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_D2z_S_Dxz_S_M3_vrr = PAZ*I_ERI_Pz_S_Dxz_S_M3_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M4_vrr+oned2z*I_ERI_S_S_Dxz_S_M3_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M4_vrr+oned2k*I_ERI_Pz_S_Px_S_M4_vrr;
      Double I_ERI_D2x_S_D2y_S_M3_vrr = PAX*I_ERI_Px_S_D2y_S_M3_vrr+WPX*I_ERI_Px_S_D2y_S_M4_vrr+oned2z*I_ERI_S_S_D2y_S_M3_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Dxy_S_D2y_S_M3_vrr = PAY*I_ERI_Px_S_D2y_S_M3_vrr+WPY*I_ERI_Px_S_D2y_S_M4_vrr+2*oned2k*I_ERI_Px_S_Py_S_M4_vrr;
      Double I_ERI_Dxz_S_D2y_S_M3_vrr = PAZ*I_ERI_Px_S_D2y_S_M3_vrr+WPZ*I_ERI_Px_S_D2y_S_M4_vrr;
      Double I_ERI_D2y_S_D2y_S_M3_vrr = PAY*I_ERI_Py_S_D2y_S_M3_vrr+WPY*I_ERI_Py_S_D2y_S_M4_vrr+oned2z*I_ERI_S_S_D2y_S_M3_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M4_vrr+2*oned2k*I_ERI_Py_S_Py_S_M4_vrr;
      Double I_ERI_Dyz_S_D2y_S_M3_vrr = PAZ*I_ERI_Py_S_D2y_S_M3_vrr+WPZ*I_ERI_Py_S_D2y_S_M4_vrr;
      Double I_ERI_D2z_S_D2y_S_M3_vrr = PAZ*I_ERI_Pz_S_D2y_S_M3_vrr+WPZ*I_ERI_Pz_S_D2y_S_M4_vrr+oned2z*I_ERI_S_S_D2y_S_M3_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_D2x_S_Dyz_S_M3_vrr = PAX*I_ERI_Px_S_Dyz_S_M3_vrr+WPX*I_ERI_Px_S_Dyz_S_M4_vrr+oned2z*I_ERI_S_S_Dyz_S_M3_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_D2y_S_Dyz_S_M3_vrr = PAY*I_ERI_Py_S_Dyz_S_M3_vrr+WPY*I_ERI_Py_S_Dyz_S_M4_vrr+oned2z*I_ERI_S_S_Dyz_S_M3_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M4_vrr+oned2k*I_ERI_Py_S_Pz_S_M4_vrr;
      Double I_ERI_D2z_S_Dyz_S_M3_vrr = PAZ*I_ERI_Pz_S_Dyz_S_M3_vrr+WPZ*I_ERI_Pz_S_Dyz_S_M4_vrr+oned2z*I_ERI_S_S_Dyz_S_M3_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M4_vrr+oned2k*I_ERI_Pz_S_Py_S_M4_vrr;
      Double I_ERI_D2x_S_D2z_S_M3_vrr = PAX*I_ERI_Px_S_D2z_S_M3_vrr+WPX*I_ERI_Px_S_D2z_S_M4_vrr+oned2z*I_ERI_S_S_D2z_S_M3_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Dxy_S_D2z_S_M3_vrr = PAY*I_ERI_Px_S_D2z_S_M3_vrr+WPY*I_ERI_Px_S_D2z_S_M4_vrr;
      Double I_ERI_Dxz_S_D2z_S_M3_vrr = PAZ*I_ERI_Px_S_D2z_S_M3_vrr+WPZ*I_ERI_Px_S_D2z_S_M4_vrr+2*oned2k*I_ERI_Px_S_Pz_S_M4_vrr;
      Double I_ERI_D2y_S_D2z_S_M3_vrr = PAY*I_ERI_Py_S_D2z_S_M3_vrr+WPY*I_ERI_Py_S_D2z_S_M4_vrr+oned2z*I_ERI_S_S_D2z_S_M3_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Dyz_S_D2z_S_M3_vrr = PAZ*I_ERI_Py_S_D2z_S_M3_vrr+WPZ*I_ERI_Py_S_D2z_S_M4_vrr+2*oned2k*I_ERI_Py_S_Pz_S_M4_vrr;
      Double I_ERI_D2z_S_D2z_S_M3_vrr = PAZ*I_ERI_Pz_S_D2z_S_M3_vrr+WPZ*I_ERI_Pz_S_D2z_S_M4_vrr+oned2z*I_ERI_S_S_D2z_S_M3_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M4_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M3_vrr = QCX*I_ERI_Px_S_D2x_S_M3_vrr+WQX*I_ERI_Px_S_D2x_S_M4_vrr+2*oned2e*I_ERI_Px_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_Px_S_Px_S_M4_vrr+oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Py_S_F3x_S_M3_vrr = QCX*I_ERI_Py_S_D2x_S_M3_vrr+WQX*I_ERI_Py_S_D2x_S_M4_vrr+2*oned2e*I_ERI_Py_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_Py_S_Px_S_M4_vrr;
      Double I_ERI_Pz_S_F3x_S_M3_vrr = QCX*I_ERI_Pz_S_D2x_S_M3_vrr+WQX*I_ERI_Pz_S_D2x_S_M4_vrr+2*oned2e*I_ERI_Pz_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_Pz_S_Px_S_M4_vrr;
      Double I_ERI_Px_S_F2xy_S_M3_vrr = QCY*I_ERI_Px_S_D2x_S_M3_vrr+WQY*I_ERI_Px_S_D2x_S_M4_vrr;
      Double I_ERI_Py_S_F2xy_S_M3_vrr = QCY*I_ERI_Py_S_D2x_S_M3_vrr+WQY*I_ERI_Py_S_D2x_S_M4_vrr+oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Pz_S_F2xy_S_M3_vrr = QCY*I_ERI_Pz_S_D2x_S_M3_vrr+WQY*I_ERI_Pz_S_D2x_S_M4_vrr;
      Double I_ERI_Px_S_F2xz_S_M3_vrr = QCZ*I_ERI_Px_S_D2x_S_M3_vrr+WQZ*I_ERI_Px_S_D2x_S_M4_vrr;
      Double I_ERI_Py_S_F2xz_S_M3_vrr = QCZ*I_ERI_Py_S_D2x_S_M3_vrr+WQZ*I_ERI_Py_S_D2x_S_M4_vrr;
      Double I_ERI_Pz_S_F2xz_S_M3_vrr = QCZ*I_ERI_Pz_S_D2x_S_M3_vrr+WQZ*I_ERI_Pz_S_D2x_S_M4_vrr+oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Px_S_Fx2y_S_M3_vrr = QCX*I_ERI_Px_S_D2y_S_M3_vrr+WQX*I_ERI_Px_S_D2y_S_M4_vrr+oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Py_S_Fx2y_S_M3_vrr = QCX*I_ERI_Py_S_D2y_S_M3_vrr+WQX*I_ERI_Py_S_D2y_S_M4_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M3_vrr = QCX*I_ERI_Pz_S_D2y_S_M3_vrr+WQX*I_ERI_Pz_S_D2y_S_M4_vrr;
      Double I_ERI_Px_S_Fxyz_S_M3_vrr = QCZ*I_ERI_Px_S_Dxy_S_M3_vrr+WQZ*I_ERI_Px_S_Dxy_S_M4_vrr;
      Double I_ERI_Py_S_Fxyz_S_M3_vrr = QCZ*I_ERI_Py_S_Dxy_S_M3_vrr+WQZ*I_ERI_Py_S_Dxy_S_M4_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M3_vrr = QCZ*I_ERI_Pz_S_Dxy_S_M3_vrr+WQZ*I_ERI_Pz_S_Dxy_S_M4_vrr+oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Px_S_Fx2z_S_M3_vrr = QCX*I_ERI_Px_S_D2z_S_M3_vrr+WQX*I_ERI_Px_S_D2z_S_M4_vrr+oned2k*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Py_S_Fx2z_S_M3_vrr = QCX*I_ERI_Py_S_D2z_S_M3_vrr+WQX*I_ERI_Py_S_D2z_S_M4_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M3_vrr = QCX*I_ERI_Pz_S_D2z_S_M3_vrr+WQX*I_ERI_Pz_S_D2z_S_M4_vrr;
      Double I_ERI_Px_S_F3y_S_M3_vrr = QCY*I_ERI_Px_S_D2y_S_M3_vrr+WQY*I_ERI_Px_S_D2y_S_M4_vrr+2*oned2e*I_ERI_Px_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_Px_S_Py_S_M4_vrr;
      Double I_ERI_Py_S_F3y_S_M3_vrr = QCY*I_ERI_Py_S_D2y_S_M3_vrr+WQY*I_ERI_Py_S_D2y_S_M4_vrr+2*oned2e*I_ERI_Py_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_Py_S_Py_S_M4_vrr+oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Pz_S_F3y_S_M3_vrr = QCY*I_ERI_Pz_S_D2y_S_M3_vrr+WQY*I_ERI_Pz_S_D2y_S_M4_vrr+2*oned2e*I_ERI_Pz_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_Pz_S_Py_S_M4_vrr;
      Double I_ERI_Px_S_F2yz_S_M3_vrr = QCZ*I_ERI_Px_S_D2y_S_M3_vrr+WQZ*I_ERI_Px_S_D2y_S_M4_vrr;
      Double I_ERI_Py_S_F2yz_S_M3_vrr = QCZ*I_ERI_Py_S_D2y_S_M3_vrr+WQZ*I_ERI_Py_S_D2y_S_M4_vrr;
      Double I_ERI_Pz_S_F2yz_S_M3_vrr = QCZ*I_ERI_Pz_S_D2y_S_M3_vrr+WQZ*I_ERI_Pz_S_D2y_S_M4_vrr+oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Px_S_Fy2z_S_M3_vrr = QCY*I_ERI_Px_S_D2z_S_M3_vrr+WQY*I_ERI_Px_S_D2z_S_M4_vrr;
      Double I_ERI_Py_S_Fy2z_S_M3_vrr = QCY*I_ERI_Py_S_D2z_S_M3_vrr+WQY*I_ERI_Py_S_D2z_S_M4_vrr+oned2k*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M3_vrr = QCY*I_ERI_Pz_S_D2z_S_M3_vrr+WQY*I_ERI_Pz_S_D2z_S_M4_vrr;
      Double I_ERI_Px_S_F3z_S_M3_vrr = QCZ*I_ERI_Px_S_D2z_S_M3_vrr+WQZ*I_ERI_Px_S_D2z_S_M4_vrr+2*oned2e*I_ERI_Px_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_Px_S_Pz_S_M4_vrr;
      Double I_ERI_Py_S_F3z_S_M3_vrr = QCZ*I_ERI_Py_S_D2z_S_M3_vrr+WQZ*I_ERI_Py_S_D2z_S_M4_vrr+2*oned2e*I_ERI_Py_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_Py_S_Pz_S_M4_vrr;
      Double I_ERI_Pz_S_F3z_S_M3_vrr = QCZ*I_ERI_Pz_S_D2z_S_M3_vrr+WQZ*I_ERI_Pz_S_D2z_S_M4_vrr+2*oned2e*I_ERI_Pz_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_Pz_S_Pz_S_M4_vrr+oned2k*I_ERI_S_S_D2z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 18 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_Px_S_M3_vrr = PAX*I_ERI_D2x_S_Px_S_M3_vrr+WPX*I_ERI_D2x_S_Px_S_M4_vrr+2*oned2z*I_ERI_Px_S_Px_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Px_S_M4_vrr+oned2k*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F2xy_S_Px_S_M3_vrr = PAY*I_ERI_D2x_S_Px_S_M3_vrr+WPY*I_ERI_D2x_S_Px_S_M4_vrr;
      Double I_ERI_F3y_S_Px_S_M3_vrr = PAY*I_ERI_D2y_S_Px_S_M3_vrr+WPY*I_ERI_D2y_S_Px_S_M4_vrr+2*oned2z*I_ERI_Py_S_Px_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Px_S_M4_vrr;
      Double I_ERI_F3z_S_Px_S_M3_vrr = PAZ*I_ERI_D2z_S_Px_S_M3_vrr+WPZ*I_ERI_D2z_S_Px_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Px_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Px_S_M4_vrr;
      Double I_ERI_F3x_S_Py_S_M3_vrr = PAX*I_ERI_D2x_S_Py_S_M3_vrr+WPX*I_ERI_D2x_S_Py_S_M4_vrr+2*oned2z*I_ERI_Px_S_Py_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Py_S_M4_vrr;
      Double I_ERI_F2xy_S_Py_S_M3_vrr = PAY*I_ERI_D2x_S_Py_S_M3_vrr+WPY*I_ERI_D2x_S_Py_S_M4_vrr+oned2k*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F3y_S_Py_S_M3_vrr = PAY*I_ERI_D2y_S_Py_S_M3_vrr+WPY*I_ERI_D2y_S_Py_S_M4_vrr+2*oned2z*I_ERI_Py_S_Py_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Py_S_M4_vrr+oned2k*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_F3z_S_Py_S_M3_vrr = PAZ*I_ERI_D2z_S_Py_S_M3_vrr+WPZ*I_ERI_D2z_S_Py_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Py_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Py_S_M4_vrr;
      Double I_ERI_F3x_S_Pz_S_M3_vrr = PAX*I_ERI_D2x_S_Pz_S_M3_vrr+WPX*I_ERI_D2x_S_Pz_S_M4_vrr+2*oned2z*I_ERI_Px_S_Pz_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Pz_S_M4_vrr;
      Double I_ERI_F2xy_S_Pz_S_M3_vrr = PAY*I_ERI_D2x_S_Pz_S_M3_vrr+WPY*I_ERI_D2x_S_Pz_S_M4_vrr;
      Double I_ERI_F3y_S_Pz_S_M3_vrr = PAY*I_ERI_D2y_S_Pz_S_M3_vrr+WPY*I_ERI_D2y_S_Pz_S_M4_vrr+2*oned2z*I_ERI_Py_S_Pz_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Pz_S_M4_vrr;
      Double I_ERI_F3z_S_Pz_S_M3_vrr = PAZ*I_ERI_D2z_S_Pz_S_M3_vrr+WPZ*I_ERI_D2z_S_Pz_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Pz_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Pz_S_M4_vrr+oned2k*I_ERI_D2z_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 19 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M3_vrr = QCX*I_ERI_D2x_S_D2x_S_M3_vrr+WQX*I_ERI_D2x_S_D2x_S_M4_vrr+2*oned2e*I_ERI_D2x_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_D2x_S_Px_S_M4_vrr+2*oned2k*I_ERI_Px_S_D2x_S_M4_vrr;
      Double I_ERI_Dxy_S_F3x_S_M3_vrr = QCX*I_ERI_Dxy_S_D2x_S_M3_vrr+WQX*I_ERI_Dxy_S_D2x_S_M4_vrr+2*oned2e*I_ERI_Dxy_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_Dxy_S_Px_S_M4_vrr+oned2k*I_ERI_Py_S_D2x_S_M4_vrr;
      Double I_ERI_Dxz_S_F3x_S_M3_vrr = QCX*I_ERI_Dxz_S_D2x_S_M3_vrr+WQX*I_ERI_Dxz_S_D2x_S_M4_vrr+2*oned2e*I_ERI_Dxz_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_Dxz_S_Px_S_M4_vrr+oned2k*I_ERI_Pz_S_D2x_S_M4_vrr;
      Double I_ERI_D2y_S_F3x_S_M3_vrr = QCX*I_ERI_D2y_S_D2x_S_M3_vrr+WQX*I_ERI_D2y_S_D2x_S_M4_vrr+2*oned2e*I_ERI_D2y_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_D2y_S_Px_S_M4_vrr;
      Double I_ERI_Dyz_S_F3x_S_M3_vrr = QCX*I_ERI_Dyz_S_D2x_S_M3_vrr+WQX*I_ERI_Dyz_S_D2x_S_M4_vrr+2*oned2e*I_ERI_Dyz_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_Dyz_S_Px_S_M4_vrr;
      Double I_ERI_D2z_S_F3x_S_M3_vrr = QCX*I_ERI_D2z_S_D2x_S_M3_vrr+WQX*I_ERI_D2z_S_D2x_S_M4_vrr+2*oned2e*I_ERI_D2z_S_Px_S_M3_vrr-2*rhod2esq*I_ERI_D2z_S_Px_S_M4_vrr;
      Double I_ERI_D2x_S_F2xy_S_M3_vrr = QCY*I_ERI_D2x_S_D2x_S_M3_vrr+WQY*I_ERI_D2x_S_D2x_S_M4_vrr;
      Double I_ERI_D2y_S_F2xy_S_M3_vrr = QCY*I_ERI_D2y_S_D2x_S_M3_vrr+WQY*I_ERI_D2y_S_D2x_S_M4_vrr+2*oned2k*I_ERI_Py_S_D2x_S_M4_vrr;
      Double I_ERI_Dyz_S_F2xy_S_M3_vrr = QCY*I_ERI_Dyz_S_D2x_S_M3_vrr+WQY*I_ERI_Dyz_S_D2x_S_M4_vrr+oned2k*I_ERI_Pz_S_D2x_S_M4_vrr;
      Double I_ERI_D2z_S_F2xy_S_M3_vrr = QCY*I_ERI_D2z_S_D2x_S_M3_vrr+WQY*I_ERI_D2z_S_D2x_S_M4_vrr;
      Double I_ERI_D2x_S_F2xz_S_M3_vrr = QCZ*I_ERI_D2x_S_D2x_S_M3_vrr+WQZ*I_ERI_D2x_S_D2x_S_M4_vrr;
      Double I_ERI_D2y_S_F2xz_S_M3_vrr = QCZ*I_ERI_D2y_S_D2x_S_M3_vrr+WQZ*I_ERI_D2y_S_D2x_S_M4_vrr;
      Double I_ERI_D2z_S_F2xz_S_M3_vrr = QCZ*I_ERI_D2z_S_D2x_S_M3_vrr+WQZ*I_ERI_D2z_S_D2x_S_M4_vrr+2*oned2k*I_ERI_Pz_S_D2x_S_M4_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M3_vrr = QCX*I_ERI_D2x_S_D2y_S_M3_vrr+WQX*I_ERI_D2x_S_D2y_S_M4_vrr+2*oned2k*I_ERI_Px_S_D2y_S_M4_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M3_vrr = QCX*I_ERI_D2y_S_D2y_S_M3_vrr+WQX*I_ERI_D2y_S_D2y_S_M4_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M3_vrr = QCX*I_ERI_D2z_S_D2y_S_M3_vrr+WQX*I_ERI_D2z_S_D2y_S_M4_vrr;
      Double I_ERI_D2x_S_Fxyz_S_M3_vrr = QCZ*I_ERI_D2x_S_Dxy_S_M3_vrr+WQZ*I_ERI_D2x_S_Dxy_S_M4_vrr;
      Double I_ERI_D2y_S_Fxyz_S_M3_vrr = QCZ*I_ERI_D2y_S_Dxy_S_M3_vrr+WQZ*I_ERI_D2y_S_Dxy_S_M4_vrr;
      Double I_ERI_D2z_S_Fxyz_S_M3_vrr = QCZ*I_ERI_D2z_S_Dxy_S_M3_vrr+WQZ*I_ERI_D2z_S_Dxy_S_M4_vrr+2*oned2k*I_ERI_Pz_S_Dxy_S_M4_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M3_vrr = QCX*I_ERI_D2x_S_D2z_S_M3_vrr+WQX*I_ERI_D2x_S_D2z_S_M4_vrr+2*oned2k*I_ERI_Px_S_D2z_S_M4_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M3_vrr = QCX*I_ERI_D2y_S_D2z_S_M3_vrr+WQX*I_ERI_D2y_S_D2z_S_M4_vrr;
      Double I_ERI_Dyz_S_Fx2z_S_M3_vrr = QCX*I_ERI_Dyz_S_D2z_S_M3_vrr+WQX*I_ERI_Dyz_S_D2z_S_M4_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M3_vrr = QCX*I_ERI_D2z_S_D2z_S_M3_vrr+WQX*I_ERI_D2z_S_D2z_S_M4_vrr;
      Double I_ERI_D2x_S_F3y_S_M3_vrr = QCY*I_ERI_D2x_S_D2y_S_M3_vrr+WQY*I_ERI_D2x_S_D2y_S_M4_vrr+2*oned2e*I_ERI_D2x_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_D2x_S_Py_S_M4_vrr;
      Double I_ERI_Dxy_S_F3y_S_M3_vrr = QCY*I_ERI_Dxy_S_D2y_S_M3_vrr+WQY*I_ERI_Dxy_S_D2y_S_M4_vrr+2*oned2e*I_ERI_Dxy_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_Dxy_S_Py_S_M4_vrr+oned2k*I_ERI_Px_S_D2y_S_M4_vrr;
      Double I_ERI_Dxz_S_F3y_S_M3_vrr = QCY*I_ERI_Dxz_S_D2y_S_M3_vrr+WQY*I_ERI_Dxz_S_D2y_S_M4_vrr+2*oned2e*I_ERI_Dxz_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_Dxz_S_Py_S_M4_vrr;
      Double I_ERI_D2y_S_F3y_S_M3_vrr = QCY*I_ERI_D2y_S_D2y_S_M3_vrr+WQY*I_ERI_D2y_S_D2y_S_M4_vrr+2*oned2e*I_ERI_D2y_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_D2y_S_Py_S_M4_vrr+2*oned2k*I_ERI_Py_S_D2y_S_M4_vrr;
      Double I_ERI_Dyz_S_F3y_S_M3_vrr = QCY*I_ERI_Dyz_S_D2y_S_M3_vrr+WQY*I_ERI_Dyz_S_D2y_S_M4_vrr+2*oned2e*I_ERI_Dyz_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_Dyz_S_Py_S_M4_vrr+oned2k*I_ERI_Pz_S_D2y_S_M4_vrr;
      Double I_ERI_D2z_S_F3y_S_M3_vrr = QCY*I_ERI_D2z_S_D2y_S_M3_vrr+WQY*I_ERI_D2z_S_D2y_S_M4_vrr+2*oned2e*I_ERI_D2z_S_Py_S_M3_vrr-2*rhod2esq*I_ERI_D2z_S_Py_S_M4_vrr;
      Double I_ERI_D2x_S_F2yz_S_M3_vrr = QCZ*I_ERI_D2x_S_D2y_S_M3_vrr+WQZ*I_ERI_D2x_S_D2y_S_M4_vrr;
      Double I_ERI_D2y_S_F2yz_S_M3_vrr = QCZ*I_ERI_D2y_S_D2y_S_M3_vrr+WQZ*I_ERI_D2y_S_D2y_S_M4_vrr;
      Double I_ERI_D2z_S_F2yz_S_M3_vrr = QCZ*I_ERI_D2z_S_D2y_S_M3_vrr+WQZ*I_ERI_D2z_S_D2y_S_M4_vrr+2*oned2k*I_ERI_Pz_S_D2y_S_M4_vrr;
      Double I_ERI_D2x_S_Fy2z_S_M3_vrr = QCY*I_ERI_D2x_S_D2z_S_M3_vrr+WQY*I_ERI_D2x_S_D2z_S_M4_vrr;
      Double I_ERI_D2y_S_Fy2z_S_M3_vrr = QCY*I_ERI_D2y_S_D2z_S_M3_vrr+WQY*I_ERI_D2y_S_D2z_S_M4_vrr+2*oned2k*I_ERI_Py_S_D2z_S_M4_vrr;
      Double I_ERI_D2z_S_Fy2z_S_M3_vrr = QCY*I_ERI_D2z_S_D2z_S_M3_vrr+WQY*I_ERI_D2z_S_D2z_S_M4_vrr;
      Double I_ERI_D2x_S_F3z_S_M3_vrr = QCZ*I_ERI_D2x_S_D2z_S_M3_vrr+WQZ*I_ERI_D2x_S_D2z_S_M4_vrr+2*oned2e*I_ERI_D2x_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_D2x_S_Pz_S_M4_vrr;
      Double I_ERI_Dxy_S_F3z_S_M3_vrr = QCZ*I_ERI_Dxy_S_D2z_S_M3_vrr+WQZ*I_ERI_Dxy_S_D2z_S_M4_vrr+2*oned2e*I_ERI_Dxy_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_Dxy_S_Pz_S_M4_vrr;
      Double I_ERI_Dxz_S_F3z_S_M3_vrr = QCZ*I_ERI_Dxz_S_D2z_S_M3_vrr+WQZ*I_ERI_Dxz_S_D2z_S_M4_vrr+2*oned2e*I_ERI_Dxz_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_Dxz_S_Pz_S_M4_vrr+oned2k*I_ERI_Px_S_D2z_S_M4_vrr;
      Double I_ERI_D2y_S_F3z_S_M3_vrr = QCZ*I_ERI_D2y_S_D2z_S_M3_vrr+WQZ*I_ERI_D2y_S_D2z_S_M4_vrr+2*oned2e*I_ERI_D2y_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_D2y_S_Pz_S_M4_vrr;
      Double I_ERI_Dyz_S_F3z_S_M3_vrr = QCZ*I_ERI_Dyz_S_D2z_S_M3_vrr+WQZ*I_ERI_Dyz_S_D2z_S_M4_vrr+2*oned2e*I_ERI_Dyz_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_Dyz_S_Pz_S_M4_vrr+oned2k*I_ERI_Py_S_D2z_S_M4_vrr;
      Double I_ERI_D2z_S_F3z_S_M3_vrr = QCZ*I_ERI_D2z_S_D2z_S_M3_vrr+WQZ*I_ERI_D2z_S_D2z_S_M4_vrr+2*oned2e*I_ERI_D2z_S_Pz_S_M3_vrr-2*rhod2esq*I_ERI_D2z_S_Pz_S_M4_vrr+2*oned2k*I_ERI_Pz_S_D2z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_M3_vrr = PAX*I_ERI_D2x_S_D2x_S_M3_vrr+WPX*I_ERI_D2x_S_D2x_S_M4_vrr+2*oned2z*I_ERI_Px_S_D2x_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_D2x_S_M4_vrr+2*oned2k*I_ERI_D2x_S_Px_S_M4_vrr;
      Double I_ERI_F2xy_S_D2x_S_M3_vrr = PAY*I_ERI_D2x_S_D2x_S_M3_vrr+WPY*I_ERI_D2x_S_D2x_S_M4_vrr;
      Double I_ERI_F2xz_S_D2x_S_M3_vrr = PAZ*I_ERI_D2x_S_D2x_S_M3_vrr+WPZ*I_ERI_D2x_S_D2x_S_M4_vrr;
      Double I_ERI_F3y_S_D2x_S_M3_vrr = PAY*I_ERI_D2y_S_D2x_S_M3_vrr+WPY*I_ERI_D2y_S_D2x_S_M4_vrr+2*oned2z*I_ERI_Py_S_D2x_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_D2x_S_M4_vrr;
      Double I_ERI_F2yz_S_D2x_S_M3_vrr = PAZ*I_ERI_D2y_S_D2x_S_M3_vrr+WPZ*I_ERI_D2y_S_D2x_S_M4_vrr;
      Double I_ERI_F3z_S_D2x_S_M3_vrr = PAZ*I_ERI_D2z_S_D2x_S_M3_vrr+WPZ*I_ERI_D2z_S_D2x_S_M4_vrr+2*oned2z*I_ERI_Pz_S_D2x_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_D2x_S_M4_vrr;
      Double I_ERI_F3x_S_Dxy_S_M3_vrr = PAX*I_ERI_D2x_S_Dxy_S_M3_vrr+WPX*I_ERI_D2x_S_Dxy_S_M4_vrr+2*oned2z*I_ERI_Px_S_Dxy_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Dxy_S_M4_vrr+oned2k*I_ERI_D2x_S_Py_S_M4_vrr;
      Double I_ERI_F2xy_S_Dxy_S_M3_vrr = PAY*I_ERI_D2x_S_Dxy_S_M3_vrr+WPY*I_ERI_D2x_S_Dxy_S_M4_vrr+oned2k*I_ERI_D2x_S_Px_S_M4_vrr;
      Double I_ERI_F3y_S_Dxy_S_M3_vrr = PAY*I_ERI_D2y_S_Dxy_S_M3_vrr+WPY*I_ERI_D2y_S_Dxy_S_M4_vrr+2*oned2z*I_ERI_Py_S_Dxy_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Dxy_S_M4_vrr+oned2k*I_ERI_D2y_S_Px_S_M4_vrr;
      Double I_ERI_F3z_S_Dxy_S_M3_vrr = PAZ*I_ERI_D2z_S_Dxy_S_M3_vrr+WPZ*I_ERI_D2z_S_Dxy_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Dxy_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxy_S_M4_vrr;
      Double I_ERI_F3x_S_Dxz_S_M3_vrr = PAX*I_ERI_D2x_S_Dxz_S_M3_vrr+WPX*I_ERI_D2x_S_Dxz_S_M4_vrr+2*oned2z*I_ERI_Px_S_Dxz_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Dxz_S_M4_vrr+oned2k*I_ERI_D2x_S_Pz_S_M4_vrr;
      Double I_ERI_F2xy_S_Dxz_S_M3_vrr = PAY*I_ERI_D2x_S_Dxz_S_M3_vrr+WPY*I_ERI_D2x_S_Dxz_S_M4_vrr;
      Double I_ERI_F3y_S_Dxz_S_M3_vrr = PAY*I_ERI_D2y_S_Dxz_S_M3_vrr+WPY*I_ERI_D2y_S_Dxz_S_M4_vrr+2*oned2z*I_ERI_Py_S_Dxz_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Dxz_S_M4_vrr;
      Double I_ERI_F3z_S_Dxz_S_M3_vrr = PAZ*I_ERI_D2z_S_Dxz_S_M3_vrr+WPZ*I_ERI_D2z_S_Dxz_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Dxz_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxz_S_M4_vrr+oned2k*I_ERI_D2z_S_Px_S_M4_vrr;
      Double I_ERI_F3x_S_D2y_S_M3_vrr = PAX*I_ERI_D2x_S_D2y_S_M3_vrr+WPX*I_ERI_D2x_S_D2y_S_M4_vrr+2*oned2z*I_ERI_Px_S_D2y_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_D2y_S_M4_vrr;
      Double I_ERI_F2xy_S_D2y_S_M3_vrr = PAY*I_ERI_D2x_S_D2y_S_M3_vrr+WPY*I_ERI_D2x_S_D2y_S_M4_vrr+2*oned2k*I_ERI_D2x_S_Py_S_M4_vrr;
      Double I_ERI_F2xz_S_D2y_S_M3_vrr = PAZ*I_ERI_D2x_S_D2y_S_M3_vrr+WPZ*I_ERI_D2x_S_D2y_S_M4_vrr;
      Double I_ERI_F3y_S_D2y_S_M3_vrr = PAY*I_ERI_D2y_S_D2y_S_M3_vrr+WPY*I_ERI_D2y_S_D2y_S_M4_vrr+2*oned2z*I_ERI_Py_S_D2y_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_D2y_S_M4_vrr+2*oned2k*I_ERI_D2y_S_Py_S_M4_vrr;
      Double I_ERI_F2yz_S_D2y_S_M3_vrr = PAZ*I_ERI_D2y_S_D2y_S_M3_vrr+WPZ*I_ERI_D2y_S_D2y_S_M4_vrr;
      Double I_ERI_F3z_S_D2y_S_M3_vrr = PAZ*I_ERI_D2z_S_D2y_S_M3_vrr+WPZ*I_ERI_D2z_S_D2y_S_M4_vrr+2*oned2z*I_ERI_Pz_S_D2y_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_D2y_S_M4_vrr;
      Double I_ERI_F3x_S_Dyz_S_M3_vrr = PAX*I_ERI_D2x_S_Dyz_S_M3_vrr+WPX*I_ERI_D2x_S_Dyz_S_M4_vrr+2*oned2z*I_ERI_Px_S_Dyz_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Dyz_S_M4_vrr;
      Double I_ERI_F2xy_S_Dyz_S_M3_vrr = PAY*I_ERI_D2x_S_Dyz_S_M3_vrr+WPY*I_ERI_D2x_S_Dyz_S_M4_vrr+oned2k*I_ERI_D2x_S_Pz_S_M4_vrr;
      Double I_ERI_F3y_S_Dyz_S_M3_vrr = PAY*I_ERI_D2y_S_Dyz_S_M3_vrr+WPY*I_ERI_D2y_S_Dyz_S_M4_vrr+2*oned2z*I_ERI_Py_S_Dyz_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Dyz_S_M4_vrr+oned2k*I_ERI_D2y_S_Pz_S_M4_vrr;
      Double I_ERI_F3z_S_Dyz_S_M3_vrr = PAZ*I_ERI_D2z_S_Dyz_S_M3_vrr+WPZ*I_ERI_D2z_S_Dyz_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Dyz_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Dyz_S_M4_vrr+oned2k*I_ERI_D2z_S_Py_S_M4_vrr;
      Double I_ERI_F3x_S_D2z_S_M3_vrr = PAX*I_ERI_D2x_S_D2z_S_M3_vrr+WPX*I_ERI_D2x_S_D2z_S_M4_vrr+2*oned2z*I_ERI_Px_S_D2z_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_D2z_S_M4_vrr;
      Double I_ERI_F2xy_S_D2z_S_M3_vrr = PAY*I_ERI_D2x_S_D2z_S_M3_vrr+WPY*I_ERI_D2x_S_D2z_S_M4_vrr;
      Double I_ERI_F2xz_S_D2z_S_M3_vrr = PAZ*I_ERI_D2x_S_D2z_S_M3_vrr+WPZ*I_ERI_D2x_S_D2z_S_M4_vrr+2*oned2k*I_ERI_D2x_S_Pz_S_M4_vrr;
      Double I_ERI_F3y_S_D2z_S_M3_vrr = PAY*I_ERI_D2y_S_D2z_S_M3_vrr+WPY*I_ERI_D2y_S_D2z_S_M4_vrr+2*oned2z*I_ERI_Py_S_D2z_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_D2z_S_M4_vrr;
      Double I_ERI_F2yz_S_D2z_S_M3_vrr = PAZ*I_ERI_D2y_S_D2z_S_M3_vrr+WPZ*I_ERI_D2y_S_D2z_S_M4_vrr+2*oned2k*I_ERI_D2y_S_Pz_S_M4_vrr;
      Double I_ERI_F3z_S_D2z_S_M3_vrr = PAZ*I_ERI_D2z_S_D2z_S_M3_vrr+WPZ*I_ERI_D2z_S_D2z_S_M4_vrr+2*oned2z*I_ERI_Pz_S_D2z_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_D2z_S_M4_vrr+2*oned2k*I_ERI_D2z_S_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 44 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_M3_vrr = PAX*I_ERI_D2x_S_F3x_S_M3_vrr+WPX*I_ERI_D2x_S_F3x_S_M4_vrr+2*oned2z*I_ERI_Px_S_F3x_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_F3x_S_M4_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M4_vrr;
      Double I_ERI_F2xy_S_F3x_S_M3_vrr = PAY*I_ERI_D2x_S_F3x_S_M3_vrr+WPY*I_ERI_D2x_S_F3x_S_M4_vrr;
      Double I_ERI_F2xz_S_F3x_S_M3_vrr = PAZ*I_ERI_D2x_S_F3x_S_M3_vrr+WPZ*I_ERI_D2x_S_F3x_S_M4_vrr;
      Double I_ERI_F3y_S_F3x_S_M3_vrr = PAY*I_ERI_D2y_S_F3x_S_M3_vrr+WPY*I_ERI_D2y_S_F3x_S_M4_vrr+2*oned2z*I_ERI_Py_S_F3x_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_F3x_S_M4_vrr;
      Double I_ERI_F2yz_S_F3x_S_M3_vrr = PAZ*I_ERI_D2y_S_F3x_S_M3_vrr+WPZ*I_ERI_D2y_S_F3x_S_M4_vrr;
      Double I_ERI_F3z_S_F3x_S_M3_vrr = PAZ*I_ERI_D2z_S_F3x_S_M3_vrr+WPZ*I_ERI_D2z_S_F3x_S_M4_vrr+2*oned2z*I_ERI_Pz_S_F3x_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_F3x_S_M4_vrr;
      Double I_ERI_F3x_S_F2xy_S_M3_vrr = PAX*I_ERI_D2x_S_F2xy_S_M3_vrr+WPX*I_ERI_D2x_S_F2xy_S_M4_vrr+2*oned2z*I_ERI_Px_S_F2xy_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_F2xy_S_M4_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M4_vrr;
      Double I_ERI_F2xy_S_F2xy_S_M3_vrr = PAY*I_ERI_D2x_S_F2xy_S_M3_vrr+WPY*I_ERI_D2x_S_F2xy_S_M4_vrr+oned2k*I_ERI_D2x_S_D2x_S_M4_vrr;
      Double I_ERI_F2xz_S_F2xy_S_M3_vrr = PAZ*I_ERI_D2x_S_F2xy_S_M3_vrr+WPZ*I_ERI_D2x_S_F2xy_S_M4_vrr;
      Double I_ERI_F3y_S_F2xy_S_M3_vrr = PAY*I_ERI_D2y_S_F2xy_S_M3_vrr+WPY*I_ERI_D2y_S_F2xy_S_M4_vrr+2*oned2z*I_ERI_Py_S_F2xy_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_F2xy_S_M4_vrr+oned2k*I_ERI_D2y_S_D2x_S_M4_vrr;
      Double I_ERI_F2yz_S_F2xy_S_M3_vrr = PAZ*I_ERI_D2y_S_F2xy_S_M3_vrr+WPZ*I_ERI_D2y_S_F2xy_S_M4_vrr;
      Double I_ERI_F3z_S_F2xy_S_M3_vrr = PAZ*I_ERI_D2z_S_F2xy_S_M3_vrr+WPZ*I_ERI_D2z_S_F2xy_S_M4_vrr+2*oned2z*I_ERI_Pz_S_F2xy_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xy_S_M4_vrr;
      Double I_ERI_F3x_S_F2xz_S_M3_vrr = PAX*I_ERI_D2x_S_F2xz_S_M3_vrr+WPX*I_ERI_D2x_S_F2xz_S_M4_vrr+2*oned2z*I_ERI_Px_S_F2xz_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_F2xz_S_M4_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M4_vrr;
      Double I_ERI_F2xy_S_F2xz_S_M3_vrr = PAY*I_ERI_D2x_S_F2xz_S_M3_vrr+WPY*I_ERI_D2x_S_F2xz_S_M4_vrr;
      Double I_ERI_F2xz_S_F2xz_S_M3_vrr = PAZ*I_ERI_D2x_S_F2xz_S_M3_vrr+WPZ*I_ERI_D2x_S_F2xz_S_M4_vrr+oned2k*I_ERI_D2x_S_D2x_S_M4_vrr;
      Double I_ERI_F3y_S_F2xz_S_M3_vrr = PAY*I_ERI_D2y_S_F2xz_S_M3_vrr+WPY*I_ERI_D2y_S_F2xz_S_M4_vrr+2*oned2z*I_ERI_Py_S_F2xz_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_F2xz_S_M4_vrr;
      Double I_ERI_F2yz_S_F2xz_S_M3_vrr = PAZ*I_ERI_D2y_S_F2xz_S_M3_vrr+WPZ*I_ERI_D2y_S_F2xz_S_M4_vrr+oned2k*I_ERI_D2y_S_D2x_S_M4_vrr;
      Double I_ERI_F3z_S_F2xz_S_M3_vrr = PAZ*I_ERI_D2z_S_F2xz_S_M3_vrr+WPZ*I_ERI_D2z_S_F2xz_S_M4_vrr+2*oned2z*I_ERI_Pz_S_F2xz_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xz_S_M4_vrr+oned2k*I_ERI_D2z_S_D2x_S_M4_vrr;
      Double I_ERI_F3x_S_Fx2y_S_M3_vrr = PAX*I_ERI_D2x_S_Fx2y_S_M3_vrr+WPX*I_ERI_D2x_S_Fx2y_S_M4_vrr+2*oned2z*I_ERI_Px_S_Fx2y_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2y_S_M4_vrr+oned2k*I_ERI_D2x_S_D2y_S_M4_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_M3_vrr = PAY*I_ERI_D2x_S_Fx2y_S_M3_vrr+WPY*I_ERI_D2x_S_Fx2y_S_M4_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M4_vrr;
      Double I_ERI_F2xz_S_Fx2y_S_M3_vrr = PAZ*I_ERI_D2x_S_Fx2y_S_M3_vrr+WPZ*I_ERI_D2x_S_Fx2y_S_M4_vrr;
      Double I_ERI_F3y_S_Fx2y_S_M3_vrr = PAY*I_ERI_D2y_S_Fx2y_S_M3_vrr+WPY*I_ERI_D2y_S_Fx2y_S_M4_vrr+2*oned2z*I_ERI_Py_S_Fx2y_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2y_S_M4_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M4_vrr;
      Double I_ERI_F2yz_S_Fx2y_S_M3_vrr = PAZ*I_ERI_D2y_S_Fx2y_S_M3_vrr+WPZ*I_ERI_D2y_S_Fx2y_S_M4_vrr;
      Double I_ERI_F3z_S_Fx2y_S_M3_vrr = PAZ*I_ERI_D2z_S_Fx2y_S_M3_vrr+WPZ*I_ERI_D2z_S_Fx2y_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Fx2y_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2y_S_M4_vrr;
      Double I_ERI_F3x_S_Fxyz_S_M3_vrr = PAX*I_ERI_D2x_S_Fxyz_S_M3_vrr+WPX*I_ERI_D2x_S_Fxyz_S_M4_vrr+2*oned2z*I_ERI_Px_S_Fxyz_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Fxyz_S_M4_vrr+oned2k*I_ERI_D2x_S_Dyz_S_M4_vrr;
      Double I_ERI_F2xy_S_Fxyz_S_M3_vrr = PAY*I_ERI_D2x_S_Fxyz_S_M3_vrr+WPY*I_ERI_D2x_S_Fxyz_S_M4_vrr+oned2k*I_ERI_D2x_S_Dxz_S_M4_vrr;
      Double I_ERI_F3y_S_Fxyz_S_M3_vrr = PAY*I_ERI_D2y_S_Fxyz_S_M3_vrr+WPY*I_ERI_D2y_S_Fxyz_S_M4_vrr+2*oned2z*I_ERI_Py_S_Fxyz_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Fxyz_S_M4_vrr+oned2k*I_ERI_D2y_S_Dxz_S_M4_vrr;
      Double I_ERI_F3z_S_Fxyz_S_M3_vrr = PAZ*I_ERI_D2z_S_Fxyz_S_M3_vrr+WPZ*I_ERI_D2z_S_Fxyz_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Fxyz_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Fxyz_S_M4_vrr+oned2k*I_ERI_D2z_S_Dxy_S_M4_vrr;
      Double I_ERI_F3x_S_Fx2z_S_M3_vrr = PAX*I_ERI_D2x_S_Fx2z_S_M3_vrr+WPX*I_ERI_D2x_S_Fx2z_S_M4_vrr+2*oned2z*I_ERI_Px_S_Fx2z_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2z_S_M4_vrr+oned2k*I_ERI_D2x_S_D2z_S_M4_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_M3_vrr = PAY*I_ERI_D2x_S_Fx2z_S_M3_vrr+WPY*I_ERI_D2x_S_Fx2z_S_M4_vrr;
      Double I_ERI_F2xz_S_Fx2z_S_M3_vrr = PAZ*I_ERI_D2x_S_Fx2z_S_M3_vrr+WPZ*I_ERI_D2x_S_Fx2z_S_M4_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M4_vrr;
      Double I_ERI_F3y_S_Fx2z_S_M3_vrr = PAY*I_ERI_D2y_S_Fx2z_S_M3_vrr+WPY*I_ERI_D2y_S_Fx2z_S_M4_vrr+2*oned2z*I_ERI_Py_S_Fx2z_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2z_S_M4_vrr;
      Double I_ERI_F2yz_S_Fx2z_S_M3_vrr = PAZ*I_ERI_D2y_S_Fx2z_S_M3_vrr+WPZ*I_ERI_D2y_S_Fx2z_S_M4_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M4_vrr;
      Double I_ERI_F3z_S_Fx2z_S_M3_vrr = PAZ*I_ERI_D2z_S_Fx2z_S_M3_vrr+WPZ*I_ERI_D2z_S_Fx2z_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Fx2z_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2z_S_M4_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M4_vrr;
      Double I_ERI_F3x_S_F3y_S_M3_vrr = PAX*I_ERI_D2x_S_F3y_S_M3_vrr+WPX*I_ERI_D2x_S_F3y_S_M4_vrr+2*oned2z*I_ERI_Px_S_F3y_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_F3y_S_M4_vrr;
      Double I_ERI_F2xy_S_F3y_S_M3_vrr = PAY*I_ERI_D2x_S_F3y_S_M3_vrr+WPY*I_ERI_D2x_S_F3y_S_M4_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M4_vrr;
      Double I_ERI_F2xz_S_F3y_S_M3_vrr = PAZ*I_ERI_D2x_S_F3y_S_M3_vrr+WPZ*I_ERI_D2x_S_F3y_S_M4_vrr;
      Double I_ERI_F3y_S_F3y_S_M3_vrr = PAY*I_ERI_D2y_S_F3y_S_M3_vrr+WPY*I_ERI_D2y_S_F3y_S_M4_vrr+2*oned2z*I_ERI_Py_S_F3y_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_F3y_S_M4_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M4_vrr;
      Double I_ERI_F2yz_S_F3y_S_M3_vrr = PAZ*I_ERI_D2y_S_F3y_S_M3_vrr+WPZ*I_ERI_D2y_S_F3y_S_M4_vrr;
      Double I_ERI_F3z_S_F3y_S_M3_vrr = PAZ*I_ERI_D2z_S_F3y_S_M3_vrr+WPZ*I_ERI_D2z_S_F3y_S_M4_vrr+2*oned2z*I_ERI_Pz_S_F3y_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_F3y_S_M4_vrr;
      Double I_ERI_F3x_S_F2yz_S_M3_vrr = PAX*I_ERI_D2x_S_F2yz_S_M3_vrr+WPX*I_ERI_D2x_S_F2yz_S_M4_vrr+2*oned2z*I_ERI_Px_S_F2yz_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_F2yz_S_M4_vrr;
      Double I_ERI_F2xy_S_F2yz_S_M3_vrr = PAY*I_ERI_D2x_S_F2yz_S_M3_vrr+WPY*I_ERI_D2x_S_F2yz_S_M4_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M4_vrr;
      Double I_ERI_F2xz_S_F2yz_S_M3_vrr = PAZ*I_ERI_D2x_S_F2yz_S_M3_vrr+WPZ*I_ERI_D2x_S_F2yz_S_M4_vrr+oned2k*I_ERI_D2x_S_D2y_S_M4_vrr;
      Double I_ERI_F3y_S_F2yz_S_M3_vrr = PAY*I_ERI_D2y_S_F2yz_S_M3_vrr+WPY*I_ERI_D2y_S_F2yz_S_M4_vrr+2*oned2z*I_ERI_Py_S_F2yz_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_F2yz_S_M4_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M4_vrr;
      Double I_ERI_F2yz_S_F2yz_S_M3_vrr = PAZ*I_ERI_D2y_S_F2yz_S_M3_vrr+WPZ*I_ERI_D2y_S_F2yz_S_M4_vrr+oned2k*I_ERI_D2y_S_D2y_S_M4_vrr;
      Double I_ERI_F3z_S_F2yz_S_M3_vrr = PAZ*I_ERI_D2z_S_F2yz_S_M3_vrr+WPZ*I_ERI_D2z_S_F2yz_S_M4_vrr+2*oned2z*I_ERI_Pz_S_F2yz_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_F2yz_S_M4_vrr+oned2k*I_ERI_D2z_S_D2y_S_M4_vrr;
      Double I_ERI_F3x_S_Fy2z_S_M3_vrr = PAX*I_ERI_D2x_S_Fy2z_S_M3_vrr+WPX*I_ERI_D2x_S_Fy2z_S_M4_vrr+2*oned2z*I_ERI_Px_S_Fy2z_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_Fy2z_S_M4_vrr;
      Double I_ERI_F2xy_S_Fy2z_S_M3_vrr = PAY*I_ERI_D2x_S_Fy2z_S_M3_vrr+WPY*I_ERI_D2x_S_Fy2z_S_M4_vrr+oned2k*I_ERI_D2x_S_D2z_S_M4_vrr;
      Double I_ERI_F3y_S_Fy2z_S_M3_vrr = PAY*I_ERI_D2y_S_Fy2z_S_M3_vrr+WPY*I_ERI_D2y_S_Fy2z_S_M4_vrr+2*oned2z*I_ERI_Py_S_Fy2z_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_Fy2z_S_M4_vrr+oned2k*I_ERI_D2y_S_D2z_S_M4_vrr;
      Double I_ERI_F3z_S_Fy2z_S_M3_vrr = PAZ*I_ERI_D2z_S_Fy2z_S_M3_vrr+WPZ*I_ERI_D2z_S_Fy2z_S_M4_vrr+2*oned2z*I_ERI_Pz_S_Fy2z_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_Fy2z_S_M4_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M4_vrr;
      Double I_ERI_F3x_S_F3z_S_M3_vrr = PAX*I_ERI_D2x_S_F3z_S_M3_vrr+WPX*I_ERI_D2x_S_F3z_S_M4_vrr+2*oned2z*I_ERI_Px_S_F3z_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_F3z_S_M4_vrr;
      Double I_ERI_F2xy_S_F3z_S_M3_vrr = PAY*I_ERI_D2x_S_F3z_S_M3_vrr+WPY*I_ERI_D2x_S_F3z_S_M4_vrr;
      Double I_ERI_F2xz_S_F3z_S_M3_vrr = PAZ*I_ERI_D2x_S_F3z_S_M3_vrr+WPZ*I_ERI_D2x_S_F3z_S_M4_vrr+3*oned2k*I_ERI_D2x_S_D2z_S_M4_vrr;
      Double I_ERI_F3y_S_F3z_S_M3_vrr = PAY*I_ERI_D2y_S_F3z_S_M3_vrr+WPY*I_ERI_D2y_S_F3z_S_M4_vrr+2*oned2z*I_ERI_Py_S_F3z_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_F3z_S_M4_vrr;
      Double I_ERI_F2yz_S_F3z_S_M3_vrr = PAZ*I_ERI_D2y_S_F3z_S_M3_vrr+WPZ*I_ERI_D2y_S_F3z_S_M4_vrr+3*oned2k*I_ERI_D2y_S_D2z_S_M4_vrr;
      Double I_ERI_F3z_S_F3z_S_M3_vrr = PAZ*I_ERI_D2z_S_F3z_S_M3_vrr+WPZ*I_ERI_D2z_S_F3z_S_M4_vrr+2*oned2z*I_ERI_Pz_S_F3z_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_F3z_S_M4_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 45 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_M3_vrr = QCX*I_ERI_D2x_S_F3x_S_M3_vrr+WQX*I_ERI_D2x_S_F3x_S_M4_vrr+3*oned2e*I_ERI_D2x_S_D2x_S_M3_vrr-3*rhod2esq*I_ERI_D2x_S_D2x_S_M4_vrr+2*oned2k*I_ERI_Px_S_F3x_S_M4_vrr;
      Double I_ERI_D2y_S_G4x_S_M3_vrr = QCX*I_ERI_D2y_S_F3x_S_M3_vrr+WQX*I_ERI_D2y_S_F3x_S_M4_vrr+3*oned2e*I_ERI_D2y_S_D2x_S_M3_vrr-3*rhod2esq*I_ERI_D2y_S_D2x_S_M4_vrr;
      Double I_ERI_D2z_S_G4x_S_M3_vrr = QCX*I_ERI_D2z_S_F3x_S_M3_vrr+WQX*I_ERI_D2z_S_F3x_S_M4_vrr+3*oned2e*I_ERI_D2z_S_D2x_S_M3_vrr-3*rhod2esq*I_ERI_D2z_S_D2x_S_M4_vrr;
      Double I_ERI_D2x_S_G3xy_S_M3_vrr = QCY*I_ERI_D2x_S_F3x_S_M3_vrr+WQY*I_ERI_D2x_S_F3x_S_M4_vrr;
      Double I_ERI_D2y_S_G3xy_S_M3_vrr = QCY*I_ERI_D2y_S_F3x_S_M3_vrr+WQY*I_ERI_D2y_S_F3x_S_M4_vrr+2*oned2k*I_ERI_Py_S_F3x_S_M4_vrr;
      Double I_ERI_D2z_S_G3xy_S_M3_vrr = QCY*I_ERI_D2z_S_F3x_S_M3_vrr+WQY*I_ERI_D2z_S_F3x_S_M4_vrr;
      Double I_ERI_D2x_S_G3xz_S_M3_vrr = QCZ*I_ERI_D2x_S_F3x_S_M3_vrr+WQZ*I_ERI_D2x_S_F3x_S_M4_vrr;
      Double I_ERI_D2y_S_G3xz_S_M3_vrr = QCZ*I_ERI_D2y_S_F3x_S_M3_vrr+WQZ*I_ERI_D2y_S_F3x_S_M4_vrr;
      Double I_ERI_D2z_S_G3xz_S_M3_vrr = QCZ*I_ERI_D2z_S_F3x_S_M3_vrr+WQZ*I_ERI_D2z_S_F3x_S_M4_vrr+2*oned2k*I_ERI_Pz_S_F3x_S_M4_vrr;
      Double I_ERI_D2x_S_G2x2y_S_M3_vrr = QCY*I_ERI_D2x_S_F2xy_S_M3_vrr+WQY*I_ERI_D2x_S_F2xy_S_M4_vrr+oned2e*I_ERI_D2x_S_D2x_S_M3_vrr-rhod2esq*I_ERI_D2x_S_D2x_S_M4_vrr;
      Double I_ERI_D2y_S_G2x2y_S_M3_vrr = QCY*I_ERI_D2y_S_F2xy_S_M3_vrr+WQY*I_ERI_D2y_S_F2xy_S_M4_vrr+oned2e*I_ERI_D2y_S_D2x_S_M3_vrr-rhod2esq*I_ERI_D2y_S_D2x_S_M4_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M4_vrr;
      Double I_ERI_D2z_S_G2x2y_S_M3_vrr = QCY*I_ERI_D2z_S_F2xy_S_M3_vrr+WQY*I_ERI_D2z_S_F2xy_S_M4_vrr+oned2e*I_ERI_D2z_S_D2x_S_M3_vrr-rhod2esq*I_ERI_D2z_S_D2x_S_M4_vrr;
      Double I_ERI_D2x_S_G2xyz_S_M3_vrr = QCZ*I_ERI_D2x_S_F2xy_S_M3_vrr+WQZ*I_ERI_D2x_S_F2xy_S_M4_vrr;
      Double I_ERI_D2y_S_G2xyz_S_M3_vrr = QCZ*I_ERI_D2y_S_F2xy_S_M3_vrr+WQZ*I_ERI_D2y_S_F2xy_S_M4_vrr;
      Double I_ERI_D2z_S_G2xyz_S_M3_vrr = QCZ*I_ERI_D2z_S_F2xy_S_M3_vrr+WQZ*I_ERI_D2z_S_F2xy_S_M4_vrr+2*oned2k*I_ERI_Pz_S_F2xy_S_M4_vrr;
      Double I_ERI_D2x_S_G2x2z_S_M3_vrr = QCZ*I_ERI_D2x_S_F2xz_S_M3_vrr+WQZ*I_ERI_D2x_S_F2xz_S_M4_vrr+oned2e*I_ERI_D2x_S_D2x_S_M3_vrr-rhod2esq*I_ERI_D2x_S_D2x_S_M4_vrr;
      Double I_ERI_D2y_S_G2x2z_S_M3_vrr = QCZ*I_ERI_D2y_S_F2xz_S_M3_vrr+WQZ*I_ERI_D2y_S_F2xz_S_M4_vrr+oned2e*I_ERI_D2y_S_D2x_S_M3_vrr-rhod2esq*I_ERI_D2y_S_D2x_S_M4_vrr;
      Double I_ERI_D2z_S_G2x2z_S_M3_vrr = QCZ*I_ERI_D2z_S_F2xz_S_M3_vrr+WQZ*I_ERI_D2z_S_F2xz_S_M4_vrr+oned2e*I_ERI_D2z_S_D2x_S_M3_vrr-rhod2esq*I_ERI_D2z_S_D2x_S_M4_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M4_vrr;
      Double I_ERI_D2x_S_Gx3y_S_M3_vrr = QCX*I_ERI_D2x_S_F3y_S_M3_vrr+WQX*I_ERI_D2x_S_F3y_S_M4_vrr+2*oned2k*I_ERI_Px_S_F3y_S_M4_vrr;
      Double I_ERI_D2y_S_Gx3y_S_M3_vrr = QCX*I_ERI_D2y_S_F3y_S_M3_vrr+WQX*I_ERI_D2y_S_F3y_S_M4_vrr;
      Double I_ERI_D2z_S_Gx3y_S_M3_vrr = QCX*I_ERI_D2z_S_F3y_S_M3_vrr+WQX*I_ERI_D2z_S_F3y_S_M4_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_M3_vrr = QCZ*I_ERI_D2x_S_Fx2y_S_M3_vrr+WQZ*I_ERI_D2x_S_Fx2y_S_M4_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_M3_vrr = QCZ*I_ERI_D2y_S_Fx2y_S_M3_vrr+WQZ*I_ERI_D2y_S_Fx2y_S_M4_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_M3_vrr = QCZ*I_ERI_D2z_S_Fx2y_S_M3_vrr+WQZ*I_ERI_D2z_S_Fx2y_S_M4_vrr+2*oned2k*I_ERI_Pz_S_Fx2y_S_M4_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_M3_vrr = QCY*I_ERI_D2x_S_Fx2z_S_M3_vrr+WQY*I_ERI_D2x_S_Fx2z_S_M4_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_M3_vrr = QCY*I_ERI_D2y_S_Fx2z_S_M3_vrr+WQY*I_ERI_D2y_S_Fx2z_S_M4_vrr+2*oned2k*I_ERI_Py_S_Fx2z_S_M4_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_M3_vrr = QCY*I_ERI_D2z_S_Fx2z_S_M3_vrr+WQY*I_ERI_D2z_S_Fx2z_S_M4_vrr;
      Double I_ERI_D2x_S_Gx3z_S_M3_vrr = QCX*I_ERI_D2x_S_F3z_S_M3_vrr+WQX*I_ERI_D2x_S_F3z_S_M4_vrr+2*oned2k*I_ERI_Px_S_F3z_S_M4_vrr;
      Double I_ERI_D2y_S_Gx3z_S_M3_vrr = QCX*I_ERI_D2y_S_F3z_S_M3_vrr+WQX*I_ERI_D2y_S_F3z_S_M4_vrr;
      Double I_ERI_D2z_S_Gx3z_S_M3_vrr = QCX*I_ERI_D2z_S_F3z_S_M3_vrr+WQX*I_ERI_D2z_S_F3z_S_M4_vrr;
      Double I_ERI_D2x_S_G4y_S_M3_vrr = QCY*I_ERI_D2x_S_F3y_S_M3_vrr+WQY*I_ERI_D2x_S_F3y_S_M4_vrr+3*oned2e*I_ERI_D2x_S_D2y_S_M3_vrr-3*rhod2esq*I_ERI_D2x_S_D2y_S_M4_vrr;
      Double I_ERI_D2y_S_G4y_S_M3_vrr = QCY*I_ERI_D2y_S_F3y_S_M3_vrr+WQY*I_ERI_D2y_S_F3y_S_M4_vrr+3*oned2e*I_ERI_D2y_S_D2y_S_M3_vrr-3*rhod2esq*I_ERI_D2y_S_D2y_S_M4_vrr+2*oned2k*I_ERI_Py_S_F3y_S_M4_vrr;
      Double I_ERI_D2z_S_G4y_S_M3_vrr = QCY*I_ERI_D2z_S_F3y_S_M3_vrr+WQY*I_ERI_D2z_S_F3y_S_M4_vrr+3*oned2e*I_ERI_D2z_S_D2y_S_M3_vrr-3*rhod2esq*I_ERI_D2z_S_D2y_S_M4_vrr;
      Double I_ERI_D2x_S_G3yz_S_M3_vrr = QCZ*I_ERI_D2x_S_F3y_S_M3_vrr+WQZ*I_ERI_D2x_S_F3y_S_M4_vrr;
      Double I_ERI_D2y_S_G3yz_S_M3_vrr = QCZ*I_ERI_D2y_S_F3y_S_M3_vrr+WQZ*I_ERI_D2y_S_F3y_S_M4_vrr;
      Double I_ERI_D2z_S_G3yz_S_M3_vrr = QCZ*I_ERI_D2z_S_F3y_S_M3_vrr+WQZ*I_ERI_D2z_S_F3y_S_M4_vrr+2*oned2k*I_ERI_Pz_S_F3y_S_M4_vrr;
      Double I_ERI_D2x_S_G2y2z_S_M3_vrr = QCZ*I_ERI_D2x_S_F2yz_S_M3_vrr+WQZ*I_ERI_D2x_S_F2yz_S_M4_vrr+oned2e*I_ERI_D2x_S_D2y_S_M3_vrr-rhod2esq*I_ERI_D2x_S_D2y_S_M4_vrr;
      Double I_ERI_D2y_S_G2y2z_S_M3_vrr = QCZ*I_ERI_D2y_S_F2yz_S_M3_vrr+WQZ*I_ERI_D2y_S_F2yz_S_M4_vrr+oned2e*I_ERI_D2y_S_D2y_S_M3_vrr-rhod2esq*I_ERI_D2y_S_D2y_S_M4_vrr;
      Double I_ERI_D2z_S_G2y2z_S_M3_vrr = QCZ*I_ERI_D2z_S_F2yz_S_M3_vrr+WQZ*I_ERI_D2z_S_F2yz_S_M4_vrr+oned2e*I_ERI_D2z_S_D2y_S_M3_vrr-rhod2esq*I_ERI_D2z_S_D2y_S_M4_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M4_vrr;
      Double I_ERI_D2x_S_Gy3z_S_M3_vrr = QCY*I_ERI_D2x_S_F3z_S_M3_vrr+WQY*I_ERI_D2x_S_F3z_S_M4_vrr;
      Double I_ERI_D2y_S_Gy3z_S_M3_vrr = QCY*I_ERI_D2y_S_F3z_S_M3_vrr+WQY*I_ERI_D2y_S_F3z_S_M4_vrr+2*oned2k*I_ERI_Py_S_F3z_S_M4_vrr;
      Double I_ERI_D2z_S_Gy3z_S_M3_vrr = QCY*I_ERI_D2z_S_F3z_S_M3_vrr+WQY*I_ERI_D2z_S_F3z_S_M4_vrr;
      Double I_ERI_D2x_S_G4z_S_M3_vrr = QCZ*I_ERI_D2x_S_F3z_S_M3_vrr+WQZ*I_ERI_D2x_S_F3z_S_M4_vrr+3*oned2e*I_ERI_D2x_S_D2z_S_M3_vrr-3*rhod2esq*I_ERI_D2x_S_D2z_S_M4_vrr;
      Double I_ERI_D2y_S_G4z_S_M3_vrr = QCZ*I_ERI_D2y_S_F3z_S_M3_vrr+WQZ*I_ERI_D2y_S_F3z_S_M4_vrr+3*oned2e*I_ERI_D2y_S_D2z_S_M3_vrr-3*rhod2esq*I_ERI_D2y_S_D2z_S_M4_vrr;
      Double I_ERI_D2z_S_G4z_S_M3_vrr = QCZ*I_ERI_D2z_S_F3z_S_M3_vrr+WQZ*I_ERI_D2z_S_F3z_S_M4_vrr+3*oned2e*I_ERI_D2z_S_D2z_S_M3_vrr-3*rhod2esq*I_ERI_D2z_S_D2z_S_M4_vrr+2*oned2k*I_ERI_Pz_S_F3z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 90 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_G4x_S_M3_vrr = QCX*I_ERI_F3x_S_F3x_S_M3_vrr+WQX*I_ERI_F3x_S_F3x_S_M4_vrr+3*oned2e*I_ERI_F3x_S_D2x_S_M3_vrr-3*rhod2esq*I_ERI_F3x_S_D2x_S_M4_vrr+3*oned2k*I_ERI_D2x_S_F3x_S_M4_vrr;
      Double I_ERI_F2xy_S_G4x_S_M3_vrr = QCX*I_ERI_F2xy_S_F3x_S_M3_vrr+WQX*I_ERI_F2xy_S_F3x_S_M4_vrr+3*oned2e*I_ERI_F2xy_S_D2x_S_M3_vrr-3*rhod2esq*I_ERI_F2xy_S_D2x_S_M4_vrr+2*oned2k*I_ERI_Dxy_S_F3x_S_M4_vrr;
      Double I_ERI_F3y_S_G4x_S_M3_vrr = QCX*I_ERI_F3y_S_F3x_S_M3_vrr+WQX*I_ERI_F3y_S_F3x_S_M4_vrr+3*oned2e*I_ERI_F3y_S_D2x_S_M3_vrr-3*rhod2esq*I_ERI_F3y_S_D2x_S_M4_vrr;
      Double I_ERI_F3z_S_G4x_S_M3_vrr = QCX*I_ERI_F3z_S_F3x_S_M3_vrr+WQX*I_ERI_F3z_S_F3x_S_M4_vrr+3*oned2e*I_ERI_F3z_S_D2x_S_M3_vrr-3*rhod2esq*I_ERI_F3z_S_D2x_S_M4_vrr;
      Double I_ERI_F3x_S_G3xy_S_M3_vrr = QCY*I_ERI_F3x_S_F3x_S_M3_vrr+WQY*I_ERI_F3x_S_F3x_S_M4_vrr;
      Double I_ERI_F2xy_S_G3xy_S_M3_vrr = QCY*I_ERI_F2xy_S_F3x_S_M3_vrr+WQY*I_ERI_F2xy_S_F3x_S_M4_vrr+oned2k*I_ERI_D2x_S_F3x_S_M4_vrr;
      Double I_ERI_F3y_S_G3xy_S_M3_vrr = QCY*I_ERI_F3y_S_F3x_S_M3_vrr+WQY*I_ERI_F3y_S_F3x_S_M4_vrr+3*oned2k*I_ERI_D2y_S_F3x_S_M4_vrr;
      Double I_ERI_F3z_S_G3xy_S_M3_vrr = QCY*I_ERI_F3z_S_F3x_S_M3_vrr+WQY*I_ERI_F3z_S_F3x_S_M4_vrr;
      Double I_ERI_F3x_S_G3xz_S_M3_vrr = QCZ*I_ERI_F3x_S_F3x_S_M3_vrr+WQZ*I_ERI_F3x_S_F3x_S_M4_vrr;
      Double I_ERI_F2xy_S_G3xz_S_M3_vrr = QCZ*I_ERI_F2xy_S_F3x_S_M3_vrr+WQZ*I_ERI_F2xy_S_F3x_S_M4_vrr;
      Double I_ERI_F3y_S_G3xz_S_M3_vrr = QCZ*I_ERI_F3y_S_F3x_S_M3_vrr+WQZ*I_ERI_F3y_S_F3x_S_M4_vrr;
      Double I_ERI_F3z_S_G3xz_S_M3_vrr = QCZ*I_ERI_F3z_S_F3x_S_M3_vrr+WQZ*I_ERI_F3z_S_F3x_S_M4_vrr+3*oned2k*I_ERI_D2z_S_F3x_S_M4_vrr;
      Double I_ERI_F3x_S_G2x2y_S_M3_vrr = QCY*I_ERI_F3x_S_F2xy_S_M3_vrr+WQY*I_ERI_F3x_S_F2xy_S_M4_vrr+oned2e*I_ERI_F3x_S_D2x_S_M3_vrr-rhod2esq*I_ERI_F3x_S_D2x_S_M4_vrr;
      Double I_ERI_F2xy_S_G2x2y_S_M3_vrr = QCY*I_ERI_F2xy_S_F2xy_S_M3_vrr+WQY*I_ERI_F2xy_S_F2xy_S_M4_vrr+oned2e*I_ERI_F2xy_S_D2x_S_M3_vrr-rhod2esq*I_ERI_F2xy_S_D2x_S_M4_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M4_vrr;
      Double I_ERI_F3y_S_G2x2y_S_M3_vrr = QCY*I_ERI_F3y_S_F2xy_S_M3_vrr+WQY*I_ERI_F3y_S_F2xy_S_M4_vrr+oned2e*I_ERI_F3y_S_D2x_S_M3_vrr-rhod2esq*I_ERI_F3y_S_D2x_S_M4_vrr+3*oned2k*I_ERI_D2y_S_F2xy_S_M4_vrr;
      Double I_ERI_F3z_S_G2x2y_S_M3_vrr = QCY*I_ERI_F3z_S_F2xy_S_M3_vrr+WQY*I_ERI_F3z_S_F2xy_S_M4_vrr+oned2e*I_ERI_F3z_S_D2x_S_M3_vrr-rhod2esq*I_ERI_F3z_S_D2x_S_M4_vrr;
      Double I_ERI_F3x_S_G2xyz_S_M3_vrr = QCZ*I_ERI_F3x_S_F2xy_S_M3_vrr+WQZ*I_ERI_F3x_S_F2xy_S_M4_vrr;
      Double I_ERI_F2xy_S_G2xyz_S_M3_vrr = QCZ*I_ERI_F2xy_S_F2xy_S_M3_vrr+WQZ*I_ERI_F2xy_S_F2xy_S_M4_vrr;
      Double I_ERI_F3y_S_G2xyz_S_M3_vrr = QCZ*I_ERI_F3y_S_F2xy_S_M3_vrr+WQZ*I_ERI_F3y_S_F2xy_S_M4_vrr;
      Double I_ERI_F3z_S_G2xyz_S_M3_vrr = QCZ*I_ERI_F3z_S_F2xy_S_M3_vrr+WQZ*I_ERI_F3z_S_F2xy_S_M4_vrr+3*oned2k*I_ERI_D2z_S_F2xy_S_M4_vrr;
      Double I_ERI_F3x_S_G2x2z_S_M3_vrr = QCZ*I_ERI_F3x_S_F2xz_S_M3_vrr+WQZ*I_ERI_F3x_S_F2xz_S_M4_vrr+oned2e*I_ERI_F3x_S_D2x_S_M3_vrr-rhod2esq*I_ERI_F3x_S_D2x_S_M4_vrr;
      Double I_ERI_F2xy_S_G2x2z_S_M3_vrr = QCZ*I_ERI_F2xy_S_F2xz_S_M3_vrr+WQZ*I_ERI_F2xy_S_F2xz_S_M4_vrr+oned2e*I_ERI_F2xy_S_D2x_S_M3_vrr-rhod2esq*I_ERI_F2xy_S_D2x_S_M4_vrr;
      Double I_ERI_F3y_S_G2x2z_S_M3_vrr = QCZ*I_ERI_F3y_S_F2xz_S_M3_vrr+WQZ*I_ERI_F3y_S_F2xz_S_M4_vrr+oned2e*I_ERI_F3y_S_D2x_S_M3_vrr-rhod2esq*I_ERI_F3y_S_D2x_S_M4_vrr;
      Double I_ERI_F3z_S_G2x2z_S_M3_vrr = QCZ*I_ERI_F3z_S_F2xz_S_M3_vrr+WQZ*I_ERI_F3z_S_F2xz_S_M4_vrr+oned2e*I_ERI_F3z_S_D2x_S_M3_vrr-rhod2esq*I_ERI_F3z_S_D2x_S_M4_vrr+3*oned2k*I_ERI_D2z_S_F2xz_S_M4_vrr;
      Double I_ERI_F3x_S_Gx3y_S_M3_vrr = QCX*I_ERI_F3x_S_F3y_S_M3_vrr+WQX*I_ERI_F3x_S_F3y_S_M4_vrr+3*oned2k*I_ERI_D2x_S_F3y_S_M4_vrr;
      Double I_ERI_F2xy_S_Gx3y_S_M3_vrr = QCX*I_ERI_F2xy_S_F3y_S_M3_vrr+WQX*I_ERI_F2xy_S_F3y_S_M4_vrr+2*oned2k*I_ERI_Dxy_S_F3y_S_M4_vrr;
      Double I_ERI_F3y_S_Gx3y_S_M3_vrr = QCX*I_ERI_F3y_S_F3y_S_M3_vrr+WQX*I_ERI_F3y_S_F3y_S_M4_vrr;
      Double I_ERI_F3z_S_Gx3y_S_M3_vrr = QCX*I_ERI_F3z_S_F3y_S_M3_vrr+WQX*I_ERI_F3z_S_F3y_S_M4_vrr;
      Double I_ERI_F3x_S_Gx2yz_S_M3_vrr = QCZ*I_ERI_F3x_S_Fx2y_S_M3_vrr+WQZ*I_ERI_F3x_S_Fx2y_S_M4_vrr;
      Double I_ERI_F2xy_S_Gx2yz_S_M3_vrr = QCZ*I_ERI_F2xy_S_Fx2y_S_M3_vrr+WQZ*I_ERI_F2xy_S_Fx2y_S_M4_vrr;
      Double I_ERI_F3y_S_Gx2yz_S_M3_vrr = QCZ*I_ERI_F3y_S_Fx2y_S_M3_vrr+WQZ*I_ERI_F3y_S_Fx2y_S_M4_vrr;
      Double I_ERI_F3z_S_Gx2yz_S_M3_vrr = QCZ*I_ERI_F3z_S_Fx2y_S_M3_vrr+WQZ*I_ERI_F3z_S_Fx2y_S_M4_vrr+3*oned2k*I_ERI_D2z_S_Fx2y_S_M4_vrr;
      Double I_ERI_F3x_S_Gxy2z_S_M3_vrr = QCY*I_ERI_F3x_S_Fx2z_S_M3_vrr+WQY*I_ERI_F3x_S_Fx2z_S_M4_vrr;
      Double I_ERI_F2xy_S_Gxy2z_S_M3_vrr = QCY*I_ERI_F2xy_S_Fx2z_S_M3_vrr+WQY*I_ERI_F2xy_S_Fx2z_S_M4_vrr+oned2k*I_ERI_D2x_S_Fx2z_S_M4_vrr;
      Double I_ERI_F3y_S_Gxy2z_S_M3_vrr = QCY*I_ERI_F3y_S_Fx2z_S_M3_vrr+WQY*I_ERI_F3y_S_Fx2z_S_M4_vrr+3*oned2k*I_ERI_D2y_S_Fx2z_S_M4_vrr;
      Double I_ERI_F3z_S_Gxy2z_S_M3_vrr = QCY*I_ERI_F3z_S_Fx2z_S_M3_vrr+WQY*I_ERI_F3z_S_Fx2z_S_M4_vrr;
      Double I_ERI_F3x_S_Gx3z_S_M3_vrr = QCX*I_ERI_F3x_S_F3z_S_M3_vrr+WQX*I_ERI_F3x_S_F3z_S_M4_vrr+3*oned2k*I_ERI_D2x_S_F3z_S_M4_vrr;
      Double I_ERI_F2xy_S_Gx3z_S_M3_vrr = QCX*I_ERI_F2xy_S_F3z_S_M3_vrr+WQX*I_ERI_F2xy_S_F3z_S_M4_vrr+2*oned2k*I_ERI_Dxy_S_F3z_S_M4_vrr;
      Double I_ERI_F3y_S_Gx3z_S_M3_vrr = QCX*I_ERI_F3y_S_F3z_S_M3_vrr+WQX*I_ERI_F3y_S_F3z_S_M4_vrr;
      Double I_ERI_F3z_S_Gx3z_S_M3_vrr = QCX*I_ERI_F3z_S_F3z_S_M3_vrr+WQX*I_ERI_F3z_S_F3z_S_M4_vrr;
      Double I_ERI_F3x_S_G4y_S_M3_vrr = QCY*I_ERI_F3x_S_F3y_S_M3_vrr+WQY*I_ERI_F3x_S_F3y_S_M4_vrr+3*oned2e*I_ERI_F3x_S_D2y_S_M3_vrr-3*rhod2esq*I_ERI_F3x_S_D2y_S_M4_vrr;
      Double I_ERI_F2xy_S_G4y_S_M3_vrr = QCY*I_ERI_F2xy_S_F3y_S_M3_vrr+WQY*I_ERI_F2xy_S_F3y_S_M4_vrr+3*oned2e*I_ERI_F2xy_S_D2y_S_M3_vrr-3*rhod2esq*I_ERI_F2xy_S_D2y_S_M4_vrr+oned2k*I_ERI_D2x_S_F3y_S_M4_vrr;
      Double I_ERI_F3y_S_G4y_S_M3_vrr = QCY*I_ERI_F3y_S_F3y_S_M3_vrr+WQY*I_ERI_F3y_S_F3y_S_M4_vrr+3*oned2e*I_ERI_F3y_S_D2y_S_M3_vrr-3*rhod2esq*I_ERI_F3y_S_D2y_S_M4_vrr+3*oned2k*I_ERI_D2y_S_F3y_S_M4_vrr;
      Double I_ERI_F3z_S_G4y_S_M3_vrr = QCY*I_ERI_F3z_S_F3y_S_M3_vrr+WQY*I_ERI_F3z_S_F3y_S_M4_vrr+3*oned2e*I_ERI_F3z_S_D2y_S_M3_vrr-3*rhod2esq*I_ERI_F3z_S_D2y_S_M4_vrr;
      Double I_ERI_F3x_S_G3yz_S_M3_vrr = QCZ*I_ERI_F3x_S_F3y_S_M3_vrr+WQZ*I_ERI_F3x_S_F3y_S_M4_vrr;
      Double I_ERI_F2xy_S_G3yz_S_M3_vrr = QCZ*I_ERI_F2xy_S_F3y_S_M3_vrr+WQZ*I_ERI_F2xy_S_F3y_S_M4_vrr;
      Double I_ERI_F3y_S_G3yz_S_M3_vrr = QCZ*I_ERI_F3y_S_F3y_S_M3_vrr+WQZ*I_ERI_F3y_S_F3y_S_M4_vrr;
      Double I_ERI_F3z_S_G3yz_S_M3_vrr = QCZ*I_ERI_F3z_S_F3y_S_M3_vrr+WQZ*I_ERI_F3z_S_F3y_S_M4_vrr+3*oned2k*I_ERI_D2z_S_F3y_S_M4_vrr;
      Double I_ERI_F3x_S_G2y2z_S_M3_vrr = QCZ*I_ERI_F3x_S_F2yz_S_M3_vrr+WQZ*I_ERI_F3x_S_F2yz_S_M4_vrr+oned2e*I_ERI_F3x_S_D2y_S_M3_vrr-rhod2esq*I_ERI_F3x_S_D2y_S_M4_vrr;
      Double I_ERI_F2xy_S_G2y2z_S_M3_vrr = QCZ*I_ERI_F2xy_S_F2yz_S_M3_vrr+WQZ*I_ERI_F2xy_S_F2yz_S_M4_vrr+oned2e*I_ERI_F2xy_S_D2y_S_M3_vrr-rhod2esq*I_ERI_F2xy_S_D2y_S_M4_vrr;
      Double I_ERI_F3y_S_G2y2z_S_M3_vrr = QCZ*I_ERI_F3y_S_F2yz_S_M3_vrr+WQZ*I_ERI_F3y_S_F2yz_S_M4_vrr+oned2e*I_ERI_F3y_S_D2y_S_M3_vrr-rhod2esq*I_ERI_F3y_S_D2y_S_M4_vrr;
      Double I_ERI_F3z_S_G2y2z_S_M3_vrr = QCZ*I_ERI_F3z_S_F2yz_S_M3_vrr+WQZ*I_ERI_F3z_S_F2yz_S_M4_vrr+oned2e*I_ERI_F3z_S_D2y_S_M3_vrr-rhod2esq*I_ERI_F3z_S_D2y_S_M4_vrr+3*oned2k*I_ERI_D2z_S_F2yz_S_M4_vrr;
      Double I_ERI_F3x_S_Gy3z_S_M3_vrr = QCY*I_ERI_F3x_S_F3z_S_M3_vrr+WQY*I_ERI_F3x_S_F3z_S_M4_vrr;
      Double I_ERI_F2xy_S_Gy3z_S_M3_vrr = QCY*I_ERI_F2xy_S_F3z_S_M3_vrr+WQY*I_ERI_F2xy_S_F3z_S_M4_vrr+oned2k*I_ERI_D2x_S_F3z_S_M4_vrr;
      Double I_ERI_F3y_S_Gy3z_S_M3_vrr = QCY*I_ERI_F3y_S_F3z_S_M3_vrr+WQY*I_ERI_F3y_S_F3z_S_M4_vrr+3*oned2k*I_ERI_D2y_S_F3z_S_M4_vrr;
      Double I_ERI_F3z_S_Gy3z_S_M3_vrr = QCY*I_ERI_F3z_S_F3z_S_M3_vrr+WQY*I_ERI_F3z_S_F3z_S_M4_vrr;
      Double I_ERI_F3x_S_G4z_S_M3_vrr = QCZ*I_ERI_F3x_S_F3z_S_M3_vrr+WQZ*I_ERI_F3x_S_F3z_S_M4_vrr+3*oned2e*I_ERI_F3x_S_D2z_S_M3_vrr-3*rhod2esq*I_ERI_F3x_S_D2z_S_M4_vrr;
      Double I_ERI_F2xy_S_G4z_S_M3_vrr = QCZ*I_ERI_F2xy_S_F3z_S_M3_vrr+WQZ*I_ERI_F2xy_S_F3z_S_M4_vrr+3*oned2e*I_ERI_F2xy_S_D2z_S_M3_vrr-3*rhod2esq*I_ERI_F2xy_S_D2z_S_M4_vrr;
      Double I_ERI_F3y_S_G4z_S_M3_vrr = QCZ*I_ERI_F3y_S_F3z_S_M3_vrr+WQZ*I_ERI_F3y_S_F3z_S_M4_vrr+3*oned2e*I_ERI_F3y_S_D2z_S_M3_vrr-3*rhod2esq*I_ERI_F3y_S_D2z_S_M4_vrr;
      Double I_ERI_F3z_S_G4z_S_M3_vrr = QCZ*I_ERI_F3z_S_F3z_S_M3_vrr+WQZ*I_ERI_F3z_S_F3z_S_M4_vrr+3*oned2e*I_ERI_F3z_S_D2z_S_M3_vrr-3*rhod2esq*I_ERI_F3z_S_D2z_S_M4_vrr+3*oned2k*I_ERI_D2z_S_F3z_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M2_vrr = QCX*I_ERI_S_S_Px_S_M2_vrr+WQX*I_ERI_S_S_Px_S_M3_vrr+oned2e*I_ERI_S_S_S_S_M2_vrr-rhod2esq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Dxy_S_M2_vrr = QCY*I_ERI_S_S_Px_S_M2_vrr+WQY*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_S_S_Dxz_S_M2_vrr = QCZ*I_ERI_S_S_Px_S_M2_vrr+WQZ*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_S_S_D2y_S_M2_vrr = QCY*I_ERI_S_S_Py_S_M2_vrr+WQY*I_ERI_S_S_Py_S_M3_vrr+oned2e*I_ERI_S_S_S_S_M2_vrr-rhod2esq*I_ERI_S_S_S_S_M3_vrr;
      Double I_ERI_S_S_Dyz_S_M2_vrr = QCZ*I_ERI_S_S_Py_S_M2_vrr+WQZ*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_S_S_D2z_S_M2_vrr = QCZ*I_ERI_S_S_Pz_S_M2_vrr+WQZ*I_ERI_S_S_Pz_S_M3_vrr+oned2e*I_ERI_S_S_S_S_M2_vrr-rhod2esq*I_ERI_S_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M2_vrr = PAX*I_ERI_S_S_D2x_S_M2_vrr+WPX*I_ERI_S_S_D2x_S_M3_vrr+2*oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Py_S_D2x_S_M2_vrr = PAY*I_ERI_S_S_D2x_S_M2_vrr+WPY*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Pz_S_D2x_S_M2_vrr = PAZ*I_ERI_S_S_D2x_S_M2_vrr+WPZ*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Px_S_Dxy_S_M2_vrr = PAX*I_ERI_S_S_Dxy_S_M2_vrr+WPX*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Py_S_Dxy_S_M2_vrr = PAY*I_ERI_S_S_Dxy_S_M2_vrr+WPY*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Pz_S_Dxy_S_M2_vrr = PAZ*I_ERI_S_S_Dxy_S_M2_vrr+WPZ*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Px_S_Dxz_S_M2_vrr = PAX*I_ERI_S_S_Dxz_S_M2_vrr+WPX*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_Py_S_Dxz_S_M2_vrr = PAY*I_ERI_S_S_Dxz_S_M2_vrr+WPY*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Pz_S_Dxz_S_M2_vrr = PAZ*I_ERI_S_S_Dxz_S_M2_vrr+WPZ*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_Px_S_D2y_S_M2_vrr = PAX*I_ERI_S_S_D2y_S_M2_vrr+WPX*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Py_S_D2y_S_M2_vrr = PAY*I_ERI_S_S_D2y_S_M2_vrr+WPY*I_ERI_S_S_D2y_S_M3_vrr+2*oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Pz_S_D2y_S_M2_vrr = PAZ*I_ERI_S_S_D2y_S_M2_vrr+WPZ*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Px_S_Dyz_S_M2_vrr = PAX*I_ERI_S_S_Dyz_S_M2_vrr+WPX*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Py_S_Dyz_S_M2_vrr = PAY*I_ERI_S_S_Dyz_S_M2_vrr+WPY*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_Pz_S_Dyz_S_M2_vrr = PAZ*I_ERI_S_S_Dyz_S_M2_vrr+WPZ*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_Px_S_D2z_S_M2_vrr = PAX*I_ERI_S_S_D2z_S_M2_vrr+WPX*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Py_S_D2z_S_M2_vrr = PAY*I_ERI_S_S_D2z_S_M2_vrr+WPY*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Pz_S_D2z_S_M2_vrr = PAZ*I_ERI_S_S_D2z_S_M2_vrr+WPZ*I_ERI_S_S_D2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M2_vrr = QCX*I_ERI_S_S_D2x_S_M2_vrr+WQX*I_ERI_S_S_D2x_S_M3_vrr+2*oned2e*I_ERI_S_S_Px_S_M2_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_S_S_F2xy_S_M2_vrr = QCY*I_ERI_S_S_D2x_S_M2_vrr+WQY*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_F2xz_S_M2_vrr = QCZ*I_ERI_S_S_D2x_S_M2_vrr+WQZ*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_Fx2y_S_M2_vrr = QCX*I_ERI_S_S_D2y_S_M2_vrr+WQX*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_Fxyz_S_M2_vrr = QCZ*I_ERI_S_S_Dxy_S_M2_vrr+WQZ*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_S_S_Fx2z_S_M2_vrr = QCX*I_ERI_S_S_D2z_S_M2_vrr+WQX*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_S_S_F3y_S_M2_vrr = QCY*I_ERI_S_S_D2y_S_M2_vrr+WQY*I_ERI_S_S_D2y_S_M3_vrr+2*oned2e*I_ERI_S_S_Py_S_M2_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_S_S_F2yz_S_M2_vrr = QCZ*I_ERI_S_S_D2y_S_M2_vrr+WQZ*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_Fy2z_S_M2_vrr = QCY*I_ERI_S_S_D2z_S_M2_vrr+WQY*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_S_S_F3z_S_M2_vrr = QCZ*I_ERI_S_S_D2z_S_M2_vrr+WQZ*I_ERI_S_S_D2z_S_M3_vrr+2*oned2e*I_ERI_S_S_Pz_S_M2_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 18 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M2_vrr = PAX*I_ERI_Px_S_D2x_S_M2_vrr+WPX*I_ERI_Px_S_D2x_S_M3_vrr+oned2z*I_ERI_S_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M3_vrr+2*oned2k*I_ERI_Px_S_Px_S_M3_vrr;
      Double I_ERI_D2y_S_D2x_S_M2_vrr = PAY*I_ERI_Py_S_D2x_S_M2_vrr+WPY*I_ERI_Py_S_D2x_S_M3_vrr+oned2z*I_ERI_S_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_D2z_S_D2x_S_M2_vrr = PAZ*I_ERI_Pz_S_D2x_S_M2_vrr+WPZ*I_ERI_Pz_S_D2x_S_M3_vrr+oned2z*I_ERI_S_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_D2x_S_Dxy_S_M2_vrr = PAX*I_ERI_Px_S_Dxy_S_M2_vrr+WPX*I_ERI_Px_S_Dxy_S_M3_vrr+oned2z*I_ERI_S_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_Px_S_Py_S_M3_vrr;
      Double I_ERI_D2y_S_Dxy_S_M2_vrr = PAY*I_ERI_Py_S_Dxy_S_M2_vrr+WPY*I_ERI_Py_S_Dxy_S_M3_vrr+oned2z*I_ERI_S_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M3_vrr+oned2k*I_ERI_Py_S_Px_S_M3_vrr;
      Double I_ERI_D2z_S_Dxy_S_M2_vrr = PAZ*I_ERI_Pz_S_Dxy_S_M2_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M3_vrr+oned2z*I_ERI_S_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_D2x_S_Dxz_S_M2_vrr = PAX*I_ERI_Px_S_Dxz_S_M2_vrr+WPX*I_ERI_Px_S_Dxz_S_M3_vrr+oned2z*I_ERI_S_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_Px_S_Pz_S_M3_vrr;
      Double I_ERI_D2y_S_Dxz_S_M2_vrr = PAY*I_ERI_Py_S_Dxz_S_M2_vrr+WPY*I_ERI_Py_S_Dxz_S_M3_vrr+oned2z*I_ERI_S_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_D2z_S_Dxz_S_M2_vrr = PAZ*I_ERI_Pz_S_Dxz_S_M2_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M3_vrr+oned2z*I_ERI_S_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M3_vrr+oned2k*I_ERI_Pz_S_Px_S_M3_vrr;
      Double I_ERI_D2x_S_D2y_S_M2_vrr = PAX*I_ERI_Px_S_D2y_S_M2_vrr+WPX*I_ERI_Px_S_D2y_S_M3_vrr+oned2z*I_ERI_S_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_D2y_S_D2y_S_M2_vrr = PAY*I_ERI_Py_S_D2y_S_M2_vrr+WPY*I_ERI_Py_S_D2y_S_M3_vrr+oned2z*I_ERI_S_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M3_vrr+2*oned2k*I_ERI_Py_S_Py_S_M3_vrr;
      Double I_ERI_D2z_S_D2y_S_M2_vrr = PAZ*I_ERI_Pz_S_D2y_S_M2_vrr+WPZ*I_ERI_Pz_S_D2y_S_M3_vrr+oned2z*I_ERI_S_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_D2x_S_Dyz_S_M2_vrr = PAX*I_ERI_Px_S_Dyz_S_M2_vrr+WPX*I_ERI_Px_S_Dyz_S_M3_vrr+oned2z*I_ERI_S_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_D2y_S_Dyz_S_M2_vrr = PAY*I_ERI_Py_S_Dyz_S_M2_vrr+WPY*I_ERI_Py_S_Dyz_S_M3_vrr+oned2z*I_ERI_S_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_Py_S_Pz_S_M3_vrr;
      Double I_ERI_D2z_S_Dyz_S_M2_vrr = PAZ*I_ERI_Pz_S_Dyz_S_M2_vrr+WPZ*I_ERI_Pz_S_Dyz_S_M3_vrr+oned2z*I_ERI_S_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M3_vrr+oned2k*I_ERI_Pz_S_Py_S_M3_vrr;
      Double I_ERI_D2x_S_D2z_S_M2_vrr = PAX*I_ERI_Px_S_D2z_S_M2_vrr+WPX*I_ERI_Px_S_D2z_S_M3_vrr+oned2z*I_ERI_S_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_D2y_S_D2z_S_M2_vrr = PAY*I_ERI_Py_S_D2z_S_M2_vrr+WPY*I_ERI_Py_S_D2z_S_M3_vrr+oned2z*I_ERI_S_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_D2z_S_D2z_S_M2_vrr = PAZ*I_ERI_Pz_S_D2z_S_M2_vrr+WPZ*I_ERI_Pz_S_D2z_S_M3_vrr+oned2z*I_ERI_S_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M2_vrr = PAX*I_ERI_S_S_F3x_S_M2_vrr+WPX*I_ERI_S_S_F3x_S_M3_vrr+3*oned2k*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Py_S_F3x_S_M2_vrr = PAY*I_ERI_S_S_F3x_S_M2_vrr+WPY*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Pz_S_F3x_S_M2_vrr = PAZ*I_ERI_S_S_F3x_S_M2_vrr+WPZ*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Px_S_F2xy_S_M2_vrr = PAX*I_ERI_S_S_F2xy_S_M2_vrr+WPX*I_ERI_S_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Py_S_F2xy_S_M2_vrr = PAY*I_ERI_S_S_F2xy_S_M2_vrr+WPY*I_ERI_S_S_F2xy_S_M3_vrr+oned2k*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Pz_S_F2xy_S_M2_vrr = PAZ*I_ERI_S_S_F2xy_S_M2_vrr+WPZ*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Px_S_F2xz_S_M2_vrr = PAX*I_ERI_S_S_F2xz_S_M2_vrr+WPX*I_ERI_S_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Py_S_F2xz_S_M2_vrr = PAY*I_ERI_S_S_F2xz_S_M2_vrr+WPY*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Pz_S_F2xz_S_M2_vrr = PAZ*I_ERI_S_S_F2xz_S_M2_vrr+WPZ*I_ERI_S_S_F2xz_S_M3_vrr+oned2k*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_Px_S_Fx2y_S_M2_vrr = PAX*I_ERI_S_S_Fx2y_S_M2_vrr+WPX*I_ERI_S_S_Fx2y_S_M3_vrr+oned2k*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Py_S_Fx2y_S_M2_vrr = PAY*I_ERI_S_S_Fx2y_S_M2_vrr+WPY*I_ERI_S_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M2_vrr = PAZ*I_ERI_S_S_Fx2y_S_M2_vrr+WPZ*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Px_S_Fxyz_S_M2_vrr = PAX*I_ERI_S_S_Fxyz_S_M2_vrr+WPX*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Py_S_Fxyz_S_M2_vrr = PAY*I_ERI_S_S_Fxyz_S_M2_vrr+WPY*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M2_vrr = PAZ*I_ERI_S_S_Fxyz_S_M2_vrr+WPZ*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_S_S_Dxy_S_M3_vrr;
      Double I_ERI_Px_S_Fx2z_S_M2_vrr = PAX*I_ERI_S_S_Fx2z_S_M2_vrr+WPX*I_ERI_S_S_Fx2z_S_M3_vrr+oned2k*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Py_S_Fx2z_S_M2_vrr = PAY*I_ERI_S_S_Fx2z_S_M2_vrr+WPY*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M2_vrr = PAZ*I_ERI_S_S_Fx2z_S_M2_vrr+WPZ*I_ERI_S_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M3_vrr;
      Double I_ERI_Px_S_F3y_S_M2_vrr = PAX*I_ERI_S_S_F3y_S_M2_vrr+WPX*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Py_S_F3y_S_M2_vrr = PAY*I_ERI_S_S_F3y_S_M2_vrr+WPY*I_ERI_S_S_F3y_S_M3_vrr+3*oned2k*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Pz_S_F3y_S_M2_vrr = PAZ*I_ERI_S_S_F3y_S_M2_vrr+WPZ*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Px_S_F2yz_S_M2_vrr = PAX*I_ERI_S_S_F2yz_S_M2_vrr+WPX*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Py_S_F2yz_S_M2_vrr = PAY*I_ERI_S_S_F2yz_S_M2_vrr+WPY*I_ERI_S_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Pz_S_F2yz_S_M2_vrr = PAZ*I_ERI_S_S_F2yz_S_M2_vrr+WPZ*I_ERI_S_S_F2yz_S_M3_vrr+oned2k*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_Px_S_Fy2z_S_M2_vrr = PAX*I_ERI_S_S_Fy2z_S_M2_vrr+WPX*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_Py_S_Fy2z_S_M2_vrr = PAY*I_ERI_S_S_Fy2z_S_M2_vrr+WPY*I_ERI_S_S_Fy2z_S_M3_vrr+oned2k*I_ERI_S_S_D2z_S_M3_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M2_vrr = PAZ*I_ERI_S_S_Fy2z_S_M2_vrr+WPZ*I_ERI_S_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M3_vrr;
      Double I_ERI_Px_S_F3z_S_M2_vrr = PAX*I_ERI_S_S_F3z_S_M2_vrr+WPX*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Py_S_F3z_S_M2_vrr = PAY*I_ERI_S_S_F3z_S_M2_vrr+WPY*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Pz_S_F3z_S_M2_vrr = PAZ*I_ERI_S_S_F3z_S_M2_vrr+WPZ*I_ERI_S_S_F3z_S_M3_vrr+3*oned2k*I_ERI_S_S_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       ************************************************************/
      Double I_ERI_S_S_G4x_S_M2_vrr = QCX*I_ERI_S_S_F3x_S_M2_vrr+WQX*I_ERI_S_S_F3x_S_M3_vrr+3*oned2e*I_ERI_S_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_G3xy_S_M2_vrr = QCY*I_ERI_S_S_F3x_S_M2_vrr+WQY*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_S_S_G3xz_S_M2_vrr = QCZ*I_ERI_S_S_F3x_S_M2_vrr+WQZ*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_S_S_G2x2y_S_M2_vrr = QCY*I_ERI_S_S_F2xy_S_M2_vrr+WQY*I_ERI_S_S_F2xy_S_M3_vrr+oned2e*I_ERI_S_S_D2x_S_M2_vrr-rhod2esq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_G2xyz_S_M2_vrr = QCZ*I_ERI_S_S_F2xy_S_M2_vrr+WQZ*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_S_S_G2x2z_S_M2_vrr = QCZ*I_ERI_S_S_F2xz_S_M2_vrr+WQZ*I_ERI_S_S_F2xz_S_M3_vrr+oned2e*I_ERI_S_S_D2x_S_M2_vrr-rhod2esq*I_ERI_S_S_D2x_S_M3_vrr;
      Double I_ERI_S_S_Gx3y_S_M2_vrr = QCX*I_ERI_S_S_F3y_S_M2_vrr+WQX*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_S_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_S_S_Fx2y_S_M2_vrr+WQZ*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_S_S_Gxy2z_S_M2_vrr = QCY*I_ERI_S_S_Fx2z_S_M2_vrr+WQY*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_S_S_Gx3z_S_M2_vrr = QCX*I_ERI_S_S_F3z_S_M2_vrr+WQX*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_S_S_G4y_S_M2_vrr = QCY*I_ERI_S_S_F3y_S_M2_vrr+WQY*I_ERI_S_S_F3y_S_M3_vrr+3*oned2e*I_ERI_S_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_G3yz_S_M2_vrr = QCZ*I_ERI_S_S_F3y_S_M2_vrr+WQZ*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_S_S_G2y2z_S_M2_vrr = QCZ*I_ERI_S_S_F2yz_S_M2_vrr+WQZ*I_ERI_S_S_F2yz_S_M3_vrr+oned2e*I_ERI_S_S_D2y_S_M2_vrr-rhod2esq*I_ERI_S_S_D2y_S_M3_vrr;
      Double I_ERI_S_S_Gy3z_S_M2_vrr = QCY*I_ERI_S_S_F3z_S_M2_vrr+WQY*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_S_S_G4z_S_M2_vrr = QCZ*I_ERI_S_S_F3z_S_M2_vrr+WQZ*I_ERI_S_S_F3z_S_M3_vrr+3*oned2e*I_ERI_S_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M2_vrr = PAX*I_ERI_Px_S_F3x_S_M2_vrr+WPX*I_ERI_Px_S_F3x_S_M3_vrr+oned2z*I_ERI_S_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M3_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M3_vrr;
      Double I_ERI_D2y_S_F3x_S_M2_vrr = PAY*I_ERI_Py_S_F3x_S_M2_vrr+WPY*I_ERI_Py_S_F3x_S_M3_vrr+oned2z*I_ERI_S_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_D2z_S_F3x_S_M2_vrr = PAZ*I_ERI_Pz_S_F3x_S_M2_vrr+WPZ*I_ERI_Pz_S_F3x_S_M3_vrr+oned2z*I_ERI_S_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_D2x_S_F2xy_S_M2_vrr = PAX*I_ERI_Px_S_F2xy_S_M2_vrr+WPX*I_ERI_Px_S_F2xy_S_M3_vrr+oned2z*I_ERI_S_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M3_vrr;
      Double I_ERI_D2y_S_F2xy_S_M2_vrr = PAY*I_ERI_Py_S_F2xy_S_M2_vrr+WPY*I_ERI_Py_S_F2xy_S_M3_vrr+oned2z*I_ERI_S_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M3_vrr+oned2k*I_ERI_Py_S_D2x_S_M3_vrr;
      Double I_ERI_D2z_S_F2xy_S_M2_vrr = PAZ*I_ERI_Pz_S_F2xy_S_M2_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M3_vrr+oned2z*I_ERI_S_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_D2x_S_F2xz_S_M2_vrr = PAX*I_ERI_Px_S_F2xz_S_M2_vrr+WPX*I_ERI_Px_S_F2xz_S_M3_vrr+oned2z*I_ERI_S_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M3_vrr;
      Double I_ERI_D2y_S_F2xz_S_M2_vrr = PAY*I_ERI_Py_S_F2xz_S_M2_vrr+WPY*I_ERI_Py_S_F2xz_S_M3_vrr+oned2z*I_ERI_S_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_D2z_S_F2xz_S_M2_vrr = PAZ*I_ERI_Pz_S_F2xz_S_M2_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M3_vrr+oned2z*I_ERI_S_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M3_vrr+oned2k*I_ERI_Pz_S_D2x_S_M3_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M2_vrr = PAX*I_ERI_Px_S_Fx2y_S_M2_vrr+WPX*I_ERI_Px_S_Fx2y_S_M3_vrr+oned2z*I_ERI_S_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M3_vrr+oned2k*I_ERI_Px_S_D2y_S_M3_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M2_vrr = PAY*I_ERI_Py_S_Fx2y_S_M2_vrr+WPY*I_ERI_Py_S_Fx2y_S_M3_vrr+oned2z*I_ERI_S_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M3_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M2_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_M2_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M3_vrr+oned2z*I_ERI_S_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2x_S_Fxyz_S_M2_vrr = PAX*I_ERI_Px_S_Fxyz_S_M2_vrr+WPX*I_ERI_Px_S_Fxyz_S_M3_vrr+oned2z*I_ERI_S_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_Px_S_Dyz_S_M3_vrr;
      Double I_ERI_D2y_S_Fxyz_S_M2_vrr = PAY*I_ERI_Py_S_Fxyz_S_M2_vrr+WPY*I_ERI_Py_S_Fxyz_S_M3_vrr+oned2z*I_ERI_S_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_Py_S_Dxz_S_M3_vrr;
      Double I_ERI_D2z_S_Fxyz_S_M2_vrr = PAZ*I_ERI_Pz_S_Fxyz_S_M2_vrr+WPZ*I_ERI_Pz_S_Fxyz_S_M3_vrr+oned2z*I_ERI_S_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M3_vrr+oned2k*I_ERI_Pz_S_Dxy_S_M3_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M2_vrr = PAX*I_ERI_Px_S_Fx2z_S_M2_vrr+WPX*I_ERI_Px_S_Fx2z_S_M3_vrr+oned2z*I_ERI_S_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M3_vrr+oned2k*I_ERI_Px_S_D2z_S_M3_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M2_vrr = PAY*I_ERI_Py_S_Fx2z_S_M2_vrr+WPY*I_ERI_Py_S_Fx2z_S_M3_vrr+oned2z*I_ERI_S_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M2_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_M2_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M3_vrr+oned2z*I_ERI_S_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M3_vrr;
      Double I_ERI_D2x_S_F3y_S_M2_vrr = PAX*I_ERI_Px_S_F3y_S_M2_vrr+WPX*I_ERI_Px_S_F3y_S_M3_vrr+oned2z*I_ERI_S_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_D2y_S_F3y_S_M2_vrr = PAY*I_ERI_Py_S_F3y_S_M2_vrr+WPY*I_ERI_Py_S_F3y_S_M3_vrr+oned2z*I_ERI_S_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M3_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M3_vrr;
      Double I_ERI_D2z_S_F3y_S_M2_vrr = PAZ*I_ERI_Pz_S_F3y_S_M2_vrr+WPZ*I_ERI_Pz_S_F3y_S_M3_vrr+oned2z*I_ERI_S_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_D2x_S_F2yz_S_M2_vrr = PAX*I_ERI_Px_S_F2yz_S_M2_vrr+WPX*I_ERI_Px_S_F2yz_S_M3_vrr+oned2z*I_ERI_S_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_D2y_S_F2yz_S_M2_vrr = PAY*I_ERI_Py_S_F2yz_S_M2_vrr+WPY*I_ERI_Py_S_F2yz_S_M3_vrr+oned2z*I_ERI_S_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M3_vrr;
      Double I_ERI_D2z_S_F2yz_S_M2_vrr = PAZ*I_ERI_Pz_S_F2yz_S_M2_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M3_vrr+oned2z*I_ERI_S_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M3_vrr+oned2k*I_ERI_Pz_S_D2y_S_M3_vrr;
      Double I_ERI_D2x_S_Fy2z_S_M2_vrr = PAX*I_ERI_Px_S_Fy2z_S_M2_vrr+WPX*I_ERI_Px_S_Fy2z_S_M3_vrr+oned2z*I_ERI_S_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_D2y_S_Fy2z_S_M2_vrr = PAY*I_ERI_Py_S_Fy2z_S_M2_vrr+WPY*I_ERI_Py_S_Fy2z_S_M3_vrr+oned2z*I_ERI_S_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M3_vrr+oned2k*I_ERI_Py_S_D2z_S_M3_vrr;
      Double I_ERI_D2z_S_Fy2z_S_M2_vrr = PAZ*I_ERI_Pz_S_Fy2z_S_M2_vrr+WPZ*I_ERI_Pz_S_Fy2z_S_M3_vrr+oned2z*I_ERI_S_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Dyz_S_M3_vrr;
      Double I_ERI_D2x_S_F3z_S_M2_vrr = PAX*I_ERI_Px_S_F3z_S_M2_vrr+WPX*I_ERI_Px_S_F3z_S_M3_vrr+oned2z*I_ERI_S_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_D2y_S_F3z_S_M2_vrr = PAY*I_ERI_Py_S_F3z_S_M2_vrr+WPY*I_ERI_Py_S_F3z_S_M3_vrr+oned2z*I_ERI_S_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_D2z_S_F3z_S_M2_vrr = PAZ*I_ERI_Pz_S_F3z_S_M2_vrr+WPZ*I_ERI_Pz_S_F3z_S_M3_vrr+oned2z*I_ERI_S_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M3_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 24 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M3
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_M2_vrr = PAX*I_ERI_D2x_S_D2x_S_M2_vrr+WPX*I_ERI_D2x_S_D2x_S_M3_vrr+2*oned2z*I_ERI_Px_S_D2x_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_D2x_S_M3_vrr+2*oned2k*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_F2xy_S_D2x_S_M2_vrr = PAY*I_ERI_D2x_S_D2x_S_M2_vrr+WPY*I_ERI_D2x_S_D2x_S_M3_vrr;
      Double I_ERI_F2xz_S_D2x_S_M2_vrr = PAZ*I_ERI_D2x_S_D2x_S_M2_vrr+WPZ*I_ERI_D2x_S_D2x_S_M3_vrr;
      Double I_ERI_F3y_S_D2x_S_M2_vrr = PAY*I_ERI_D2y_S_D2x_S_M2_vrr+WPY*I_ERI_D2y_S_D2x_S_M3_vrr+2*oned2z*I_ERI_Py_S_D2x_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_D2x_S_M3_vrr;
      Double I_ERI_F2yz_S_D2x_S_M2_vrr = PAZ*I_ERI_D2y_S_D2x_S_M2_vrr+WPZ*I_ERI_D2y_S_D2x_S_M3_vrr;
      Double I_ERI_F3z_S_D2x_S_M2_vrr = PAZ*I_ERI_D2z_S_D2x_S_M2_vrr+WPZ*I_ERI_D2z_S_D2x_S_M3_vrr+2*oned2z*I_ERI_Pz_S_D2x_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_D2x_S_M3_vrr;
      Double I_ERI_F3x_S_Dxy_S_M2_vrr = PAX*I_ERI_D2x_S_Dxy_S_M2_vrr+WPX*I_ERI_D2x_S_Dxy_S_M3_vrr+2*oned2z*I_ERI_Px_S_Dxy_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Dxy_S_M3_vrr+oned2k*I_ERI_D2x_S_Py_S_M3_vrr;
      Double I_ERI_F2xy_S_Dxy_S_M2_vrr = PAY*I_ERI_D2x_S_Dxy_S_M2_vrr+WPY*I_ERI_D2x_S_Dxy_S_M3_vrr+oned2k*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_F2xz_S_Dxy_S_M2_vrr = PAZ*I_ERI_D2x_S_Dxy_S_M2_vrr+WPZ*I_ERI_D2x_S_Dxy_S_M3_vrr;
      Double I_ERI_F3y_S_Dxy_S_M2_vrr = PAY*I_ERI_D2y_S_Dxy_S_M2_vrr+WPY*I_ERI_D2y_S_Dxy_S_M3_vrr+2*oned2z*I_ERI_Py_S_Dxy_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Dxy_S_M3_vrr+oned2k*I_ERI_D2y_S_Px_S_M3_vrr;
      Double I_ERI_F2yz_S_Dxy_S_M2_vrr = PAZ*I_ERI_D2y_S_Dxy_S_M2_vrr+WPZ*I_ERI_D2y_S_Dxy_S_M3_vrr;
      Double I_ERI_F3z_S_Dxy_S_M2_vrr = PAZ*I_ERI_D2z_S_Dxy_S_M2_vrr+WPZ*I_ERI_D2z_S_Dxy_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Dxy_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxy_S_M3_vrr;
      Double I_ERI_F3x_S_Dxz_S_M2_vrr = PAX*I_ERI_D2x_S_Dxz_S_M2_vrr+WPX*I_ERI_D2x_S_Dxz_S_M3_vrr+2*oned2z*I_ERI_Px_S_Dxz_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Dxz_S_M3_vrr+oned2k*I_ERI_D2x_S_Pz_S_M3_vrr;
      Double I_ERI_F2xy_S_Dxz_S_M2_vrr = PAY*I_ERI_D2x_S_Dxz_S_M2_vrr+WPY*I_ERI_D2x_S_Dxz_S_M3_vrr;
      Double I_ERI_F2xz_S_Dxz_S_M2_vrr = PAZ*I_ERI_D2x_S_Dxz_S_M2_vrr+WPZ*I_ERI_D2x_S_Dxz_S_M3_vrr+oned2k*I_ERI_D2x_S_Px_S_M3_vrr;
      Double I_ERI_F3y_S_Dxz_S_M2_vrr = PAY*I_ERI_D2y_S_Dxz_S_M2_vrr+WPY*I_ERI_D2y_S_Dxz_S_M3_vrr+2*oned2z*I_ERI_Py_S_Dxz_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Dxz_S_M3_vrr;
      Double I_ERI_F2yz_S_Dxz_S_M2_vrr = PAZ*I_ERI_D2y_S_Dxz_S_M2_vrr+WPZ*I_ERI_D2y_S_Dxz_S_M3_vrr+oned2k*I_ERI_D2y_S_Px_S_M3_vrr;
      Double I_ERI_F3z_S_Dxz_S_M2_vrr = PAZ*I_ERI_D2z_S_Dxz_S_M2_vrr+WPZ*I_ERI_D2z_S_Dxz_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Dxz_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxz_S_M3_vrr+oned2k*I_ERI_D2z_S_Px_S_M3_vrr;
      Double I_ERI_F3x_S_D2y_S_M2_vrr = PAX*I_ERI_D2x_S_D2y_S_M2_vrr+WPX*I_ERI_D2x_S_D2y_S_M3_vrr+2*oned2z*I_ERI_Px_S_D2y_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_D2y_S_M3_vrr;
      Double I_ERI_F2xy_S_D2y_S_M2_vrr = PAY*I_ERI_D2x_S_D2y_S_M2_vrr+WPY*I_ERI_D2x_S_D2y_S_M3_vrr+2*oned2k*I_ERI_D2x_S_Py_S_M3_vrr;
      Double I_ERI_F2xz_S_D2y_S_M2_vrr = PAZ*I_ERI_D2x_S_D2y_S_M2_vrr+WPZ*I_ERI_D2x_S_D2y_S_M3_vrr;
      Double I_ERI_F3y_S_D2y_S_M2_vrr = PAY*I_ERI_D2y_S_D2y_S_M2_vrr+WPY*I_ERI_D2y_S_D2y_S_M3_vrr+2*oned2z*I_ERI_Py_S_D2y_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_D2y_S_M3_vrr+2*oned2k*I_ERI_D2y_S_Py_S_M3_vrr;
      Double I_ERI_F2yz_S_D2y_S_M2_vrr = PAZ*I_ERI_D2y_S_D2y_S_M2_vrr+WPZ*I_ERI_D2y_S_D2y_S_M3_vrr;
      Double I_ERI_F3z_S_D2y_S_M2_vrr = PAZ*I_ERI_D2z_S_D2y_S_M2_vrr+WPZ*I_ERI_D2z_S_D2y_S_M3_vrr+2*oned2z*I_ERI_Pz_S_D2y_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_D2y_S_M3_vrr;
      Double I_ERI_F3x_S_Dyz_S_M2_vrr = PAX*I_ERI_D2x_S_Dyz_S_M2_vrr+WPX*I_ERI_D2x_S_Dyz_S_M3_vrr+2*oned2z*I_ERI_Px_S_Dyz_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Dyz_S_M3_vrr;
      Double I_ERI_F2xy_S_Dyz_S_M2_vrr = PAY*I_ERI_D2x_S_Dyz_S_M2_vrr+WPY*I_ERI_D2x_S_Dyz_S_M3_vrr+oned2k*I_ERI_D2x_S_Pz_S_M3_vrr;
      Double I_ERI_F2xz_S_Dyz_S_M2_vrr = PAZ*I_ERI_D2x_S_Dyz_S_M2_vrr+WPZ*I_ERI_D2x_S_Dyz_S_M3_vrr+oned2k*I_ERI_D2x_S_Py_S_M3_vrr;
      Double I_ERI_F3y_S_Dyz_S_M2_vrr = PAY*I_ERI_D2y_S_Dyz_S_M2_vrr+WPY*I_ERI_D2y_S_Dyz_S_M3_vrr+2*oned2z*I_ERI_Py_S_Dyz_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Dyz_S_M3_vrr+oned2k*I_ERI_D2y_S_Pz_S_M3_vrr;
      Double I_ERI_F2yz_S_Dyz_S_M2_vrr = PAZ*I_ERI_D2y_S_Dyz_S_M2_vrr+WPZ*I_ERI_D2y_S_Dyz_S_M3_vrr+oned2k*I_ERI_D2y_S_Py_S_M3_vrr;
      Double I_ERI_F3z_S_Dyz_S_M2_vrr = PAZ*I_ERI_D2z_S_Dyz_S_M2_vrr+WPZ*I_ERI_D2z_S_Dyz_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Dyz_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Dyz_S_M3_vrr+oned2k*I_ERI_D2z_S_Py_S_M3_vrr;
      Double I_ERI_F3x_S_D2z_S_M2_vrr = PAX*I_ERI_D2x_S_D2z_S_M2_vrr+WPX*I_ERI_D2x_S_D2z_S_M3_vrr+2*oned2z*I_ERI_Px_S_D2z_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_D2z_S_M3_vrr;
      Double I_ERI_F2xy_S_D2z_S_M2_vrr = PAY*I_ERI_D2x_S_D2z_S_M2_vrr+WPY*I_ERI_D2x_S_D2z_S_M3_vrr;
      Double I_ERI_F2xz_S_D2z_S_M2_vrr = PAZ*I_ERI_D2x_S_D2z_S_M2_vrr+WPZ*I_ERI_D2x_S_D2z_S_M3_vrr+2*oned2k*I_ERI_D2x_S_Pz_S_M3_vrr;
      Double I_ERI_F3y_S_D2z_S_M2_vrr = PAY*I_ERI_D2y_S_D2z_S_M2_vrr+WPY*I_ERI_D2y_S_D2z_S_M3_vrr+2*oned2z*I_ERI_Py_S_D2z_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_D2z_S_M3_vrr;
      Double I_ERI_F2yz_S_D2z_S_M2_vrr = PAZ*I_ERI_D2y_S_D2z_S_M2_vrr+WPZ*I_ERI_D2y_S_D2z_S_M3_vrr+2*oned2k*I_ERI_D2y_S_Pz_S_M3_vrr;
      Double I_ERI_F3z_S_D2z_S_M2_vrr = PAZ*I_ERI_D2z_S_D2z_S_M2_vrr+WPZ*I_ERI_D2z_S_D2z_S_M3_vrr+2*oned2z*I_ERI_Pz_S_D2z_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_D2z_S_M3_vrr+2*oned2k*I_ERI_D2z_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_G_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       ************************************************************/
      Double I_ERI_Px_S_G4x_S_M2_vrr = QCX*I_ERI_Px_S_F3x_S_M2_vrr+WQX*I_ERI_Px_S_F3x_S_M3_vrr+3*oned2e*I_ERI_Px_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_Px_S_D2x_S_M3_vrr+oned2k*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Py_S_G4x_S_M2_vrr = QCX*I_ERI_Py_S_F3x_S_M2_vrr+WQX*I_ERI_Py_S_F3x_S_M3_vrr+3*oned2e*I_ERI_Py_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_Py_S_D2x_S_M3_vrr;
      Double I_ERI_Pz_S_G4x_S_M2_vrr = QCX*I_ERI_Pz_S_F3x_S_M2_vrr+WQX*I_ERI_Pz_S_F3x_S_M3_vrr+3*oned2e*I_ERI_Pz_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_Pz_S_D2x_S_M3_vrr;
      Double I_ERI_Px_S_G3xy_S_M2_vrr = QCY*I_ERI_Px_S_F3x_S_M2_vrr+WQY*I_ERI_Px_S_F3x_S_M3_vrr;
      Double I_ERI_Py_S_G3xy_S_M2_vrr = QCY*I_ERI_Py_S_F3x_S_M2_vrr+WQY*I_ERI_Py_S_F3x_S_M3_vrr+oned2k*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Pz_S_G3xy_S_M2_vrr = QCY*I_ERI_Pz_S_F3x_S_M2_vrr+WQY*I_ERI_Pz_S_F3x_S_M3_vrr;
      Double I_ERI_Px_S_G3xz_S_M2_vrr = QCZ*I_ERI_Px_S_F3x_S_M2_vrr+WQZ*I_ERI_Px_S_F3x_S_M3_vrr;
      Double I_ERI_Py_S_G3xz_S_M2_vrr = QCZ*I_ERI_Py_S_F3x_S_M2_vrr+WQZ*I_ERI_Py_S_F3x_S_M3_vrr;
      Double I_ERI_Pz_S_G3xz_S_M2_vrr = QCZ*I_ERI_Pz_S_F3x_S_M2_vrr+WQZ*I_ERI_Pz_S_F3x_S_M3_vrr+oned2k*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Px_S_G2x2y_S_M2_vrr = QCY*I_ERI_Px_S_F2xy_S_M2_vrr+WQY*I_ERI_Px_S_F2xy_S_M3_vrr+oned2e*I_ERI_Px_S_D2x_S_M2_vrr-rhod2esq*I_ERI_Px_S_D2x_S_M3_vrr;
      Double I_ERI_Py_S_G2x2y_S_M2_vrr = QCY*I_ERI_Py_S_F2xy_S_M2_vrr+WQY*I_ERI_Py_S_F2xy_S_M3_vrr+oned2e*I_ERI_Py_S_D2x_S_M2_vrr-rhod2esq*I_ERI_Py_S_D2x_S_M3_vrr+oned2k*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Pz_S_G2x2y_S_M2_vrr = QCY*I_ERI_Pz_S_F2xy_S_M2_vrr+WQY*I_ERI_Pz_S_F2xy_S_M3_vrr+oned2e*I_ERI_Pz_S_D2x_S_M2_vrr-rhod2esq*I_ERI_Pz_S_D2x_S_M3_vrr;
      Double I_ERI_Px_S_G2xyz_S_M2_vrr = QCZ*I_ERI_Px_S_F2xy_S_M2_vrr+WQZ*I_ERI_Px_S_F2xy_S_M3_vrr;
      Double I_ERI_Py_S_G2xyz_S_M2_vrr = QCZ*I_ERI_Py_S_F2xy_S_M2_vrr+WQZ*I_ERI_Py_S_F2xy_S_M3_vrr;
      Double I_ERI_Pz_S_G2xyz_S_M2_vrr = QCZ*I_ERI_Pz_S_F2xy_S_M2_vrr+WQZ*I_ERI_Pz_S_F2xy_S_M3_vrr+oned2k*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Px_S_G2x2z_S_M2_vrr = QCZ*I_ERI_Px_S_F2xz_S_M2_vrr+WQZ*I_ERI_Px_S_F2xz_S_M3_vrr+oned2e*I_ERI_Px_S_D2x_S_M2_vrr-rhod2esq*I_ERI_Px_S_D2x_S_M3_vrr;
      Double I_ERI_Py_S_G2x2z_S_M2_vrr = QCZ*I_ERI_Py_S_F2xz_S_M2_vrr+WQZ*I_ERI_Py_S_F2xz_S_M3_vrr+oned2e*I_ERI_Py_S_D2x_S_M2_vrr-rhod2esq*I_ERI_Py_S_D2x_S_M3_vrr;
      Double I_ERI_Pz_S_G2x2z_S_M2_vrr = QCZ*I_ERI_Pz_S_F2xz_S_M2_vrr+WQZ*I_ERI_Pz_S_F2xz_S_M3_vrr+oned2e*I_ERI_Pz_S_D2x_S_M2_vrr-rhod2esq*I_ERI_Pz_S_D2x_S_M3_vrr+oned2k*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Px_S_Gx3y_S_M2_vrr = QCX*I_ERI_Px_S_F3y_S_M2_vrr+WQX*I_ERI_Px_S_F3y_S_M3_vrr+oned2k*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Py_S_Gx3y_S_M2_vrr = QCX*I_ERI_Py_S_F3y_S_M2_vrr+WQX*I_ERI_Py_S_F3y_S_M3_vrr;
      Double I_ERI_Pz_S_Gx3y_S_M2_vrr = QCX*I_ERI_Pz_S_F3y_S_M2_vrr+WQX*I_ERI_Pz_S_F3y_S_M3_vrr;
      Double I_ERI_Px_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_Px_S_Fx2y_S_M2_vrr+WQZ*I_ERI_Px_S_Fx2y_S_M3_vrr;
      Double I_ERI_Py_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_Py_S_Fx2y_S_M2_vrr+WQZ*I_ERI_Py_S_Fx2y_S_M3_vrr;
      Double I_ERI_Pz_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_Pz_S_Fx2y_S_M2_vrr+WQZ*I_ERI_Pz_S_Fx2y_S_M3_vrr+oned2k*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Px_S_Gxy2z_S_M2_vrr = QCY*I_ERI_Px_S_Fx2z_S_M2_vrr+WQY*I_ERI_Px_S_Fx2z_S_M3_vrr;
      Double I_ERI_Py_S_Gxy2z_S_M2_vrr = QCY*I_ERI_Py_S_Fx2z_S_M2_vrr+WQY*I_ERI_Py_S_Fx2z_S_M3_vrr+oned2k*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Pz_S_Gxy2z_S_M2_vrr = QCY*I_ERI_Pz_S_Fx2z_S_M2_vrr+WQY*I_ERI_Pz_S_Fx2z_S_M3_vrr;
      Double I_ERI_Px_S_Gx3z_S_M2_vrr = QCX*I_ERI_Px_S_F3z_S_M2_vrr+WQX*I_ERI_Px_S_F3z_S_M3_vrr+oned2k*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Py_S_Gx3z_S_M2_vrr = QCX*I_ERI_Py_S_F3z_S_M2_vrr+WQX*I_ERI_Py_S_F3z_S_M3_vrr;
      Double I_ERI_Pz_S_Gx3z_S_M2_vrr = QCX*I_ERI_Pz_S_F3z_S_M2_vrr+WQX*I_ERI_Pz_S_F3z_S_M3_vrr;
      Double I_ERI_Px_S_G4y_S_M2_vrr = QCY*I_ERI_Px_S_F3y_S_M2_vrr+WQY*I_ERI_Px_S_F3y_S_M3_vrr+3*oned2e*I_ERI_Px_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_Px_S_D2y_S_M3_vrr;
      Double I_ERI_Py_S_G4y_S_M2_vrr = QCY*I_ERI_Py_S_F3y_S_M2_vrr+WQY*I_ERI_Py_S_F3y_S_M3_vrr+3*oned2e*I_ERI_Py_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_Py_S_D2y_S_M3_vrr+oned2k*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Pz_S_G4y_S_M2_vrr = QCY*I_ERI_Pz_S_F3y_S_M2_vrr+WQY*I_ERI_Pz_S_F3y_S_M3_vrr+3*oned2e*I_ERI_Pz_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_Pz_S_D2y_S_M3_vrr;
      Double I_ERI_Px_S_G3yz_S_M2_vrr = QCZ*I_ERI_Px_S_F3y_S_M2_vrr+WQZ*I_ERI_Px_S_F3y_S_M3_vrr;
      Double I_ERI_Py_S_G3yz_S_M2_vrr = QCZ*I_ERI_Py_S_F3y_S_M2_vrr+WQZ*I_ERI_Py_S_F3y_S_M3_vrr;
      Double I_ERI_Pz_S_G3yz_S_M2_vrr = QCZ*I_ERI_Pz_S_F3y_S_M2_vrr+WQZ*I_ERI_Pz_S_F3y_S_M3_vrr+oned2k*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Px_S_G2y2z_S_M2_vrr = QCZ*I_ERI_Px_S_F2yz_S_M2_vrr+WQZ*I_ERI_Px_S_F2yz_S_M3_vrr+oned2e*I_ERI_Px_S_D2y_S_M2_vrr-rhod2esq*I_ERI_Px_S_D2y_S_M3_vrr;
      Double I_ERI_Py_S_G2y2z_S_M2_vrr = QCZ*I_ERI_Py_S_F2yz_S_M2_vrr+WQZ*I_ERI_Py_S_F2yz_S_M3_vrr+oned2e*I_ERI_Py_S_D2y_S_M2_vrr-rhod2esq*I_ERI_Py_S_D2y_S_M3_vrr;
      Double I_ERI_Pz_S_G2y2z_S_M2_vrr = QCZ*I_ERI_Pz_S_F2yz_S_M2_vrr+WQZ*I_ERI_Pz_S_F2yz_S_M3_vrr+oned2e*I_ERI_Pz_S_D2y_S_M2_vrr-rhod2esq*I_ERI_Pz_S_D2y_S_M3_vrr+oned2k*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Px_S_Gy3z_S_M2_vrr = QCY*I_ERI_Px_S_F3z_S_M2_vrr+WQY*I_ERI_Px_S_F3z_S_M3_vrr;
      Double I_ERI_Py_S_Gy3z_S_M2_vrr = QCY*I_ERI_Py_S_F3z_S_M2_vrr+WQY*I_ERI_Py_S_F3z_S_M3_vrr+oned2k*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Pz_S_Gy3z_S_M2_vrr = QCY*I_ERI_Pz_S_F3z_S_M2_vrr+WQY*I_ERI_Pz_S_F3z_S_M3_vrr;
      Double I_ERI_Px_S_G4z_S_M2_vrr = QCZ*I_ERI_Px_S_F3z_S_M2_vrr+WQZ*I_ERI_Px_S_F3z_S_M3_vrr+3*oned2e*I_ERI_Px_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_Px_S_D2z_S_M3_vrr;
      Double I_ERI_Py_S_G4z_S_M2_vrr = QCZ*I_ERI_Py_S_F3z_S_M2_vrr+WQZ*I_ERI_Py_S_F3z_S_M3_vrr+3*oned2e*I_ERI_Py_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_Py_S_D2z_S_M3_vrr;
      Double I_ERI_Pz_S_G4z_S_M2_vrr = QCZ*I_ERI_Pz_S_F3z_S_M2_vrr+WQZ*I_ERI_Pz_S_F3z_S_M3_vrr+3*oned2e*I_ERI_Pz_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_Pz_S_D2z_S_M3_vrr+oned2k*I_ERI_S_S_F3z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 40 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M3
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_M2_vrr = PAX*I_ERI_D2x_S_F3x_S_M2_vrr+WPX*I_ERI_D2x_S_F3x_S_M3_vrr+2*oned2z*I_ERI_Px_S_F3x_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_F3x_S_M3_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M3_vrr;
      Double I_ERI_F2xy_S_F3x_S_M2_vrr = PAY*I_ERI_D2x_S_F3x_S_M2_vrr+WPY*I_ERI_D2x_S_F3x_S_M3_vrr;
      Double I_ERI_F2xz_S_F3x_S_M2_vrr = PAZ*I_ERI_D2x_S_F3x_S_M2_vrr+WPZ*I_ERI_D2x_S_F3x_S_M3_vrr;
      Double I_ERI_F3y_S_F3x_S_M2_vrr = PAY*I_ERI_D2y_S_F3x_S_M2_vrr+WPY*I_ERI_D2y_S_F3x_S_M3_vrr+2*oned2z*I_ERI_Py_S_F3x_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_F3x_S_M3_vrr;
      Double I_ERI_F2yz_S_F3x_S_M2_vrr = PAZ*I_ERI_D2y_S_F3x_S_M2_vrr+WPZ*I_ERI_D2y_S_F3x_S_M3_vrr;
      Double I_ERI_F3z_S_F3x_S_M2_vrr = PAZ*I_ERI_D2z_S_F3x_S_M2_vrr+WPZ*I_ERI_D2z_S_F3x_S_M3_vrr+2*oned2z*I_ERI_Pz_S_F3x_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_F3x_S_M3_vrr;
      Double I_ERI_F3x_S_F2xy_S_M2_vrr = PAX*I_ERI_D2x_S_F2xy_S_M2_vrr+WPX*I_ERI_D2x_S_F2xy_S_M3_vrr+2*oned2z*I_ERI_Px_S_F2xy_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M3_vrr;
      Double I_ERI_F2xy_S_F2xy_S_M2_vrr = PAY*I_ERI_D2x_S_F2xy_S_M2_vrr+WPY*I_ERI_D2x_S_F2xy_S_M3_vrr+oned2k*I_ERI_D2x_S_D2x_S_M3_vrr;
      Double I_ERI_F2xz_S_F2xy_S_M2_vrr = PAZ*I_ERI_D2x_S_F2xy_S_M2_vrr+WPZ*I_ERI_D2x_S_F2xy_S_M3_vrr;
      Double I_ERI_F3y_S_F2xy_S_M2_vrr = PAY*I_ERI_D2y_S_F2xy_S_M2_vrr+WPY*I_ERI_D2y_S_F2xy_S_M3_vrr+2*oned2z*I_ERI_Py_S_F2xy_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_F2xy_S_M3_vrr+oned2k*I_ERI_D2y_S_D2x_S_M3_vrr;
      Double I_ERI_F2yz_S_F2xy_S_M2_vrr = PAZ*I_ERI_D2y_S_F2xy_S_M2_vrr+WPZ*I_ERI_D2y_S_F2xy_S_M3_vrr;
      Double I_ERI_F3z_S_F2xy_S_M2_vrr = PAZ*I_ERI_D2z_S_F2xy_S_M2_vrr+WPZ*I_ERI_D2z_S_F2xy_S_M3_vrr+2*oned2z*I_ERI_Pz_S_F2xy_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xy_S_M3_vrr;
      Double I_ERI_F3x_S_F2xz_S_M2_vrr = PAX*I_ERI_D2x_S_F2xz_S_M2_vrr+WPX*I_ERI_D2x_S_F2xz_S_M3_vrr+2*oned2z*I_ERI_Px_S_F2xz_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M3_vrr;
      Double I_ERI_F2xy_S_F2xz_S_M2_vrr = PAY*I_ERI_D2x_S_F2xz_S_M2_vrr+WPY*I_ERI_D2x_S_F2xz_S_M3_vrr;
      Double I_ERI_F2xz_S_F2xz_S_M2_vrr = PAZ*I_ERI_D2x_S_F2xz_S_M2_vrr+WPZ*I_ERI_D2x_S_F2xz_S_M3_vrr+oned2k*I_ERI_D2x_S_D2x_S_M3_vrr;
      Double I_ERI_F3y_S_F2xz_S_M2_vrr = PAY*I_ERI_D2y_S_F2xz_S_M2_vrr+WPY*I_ERI_D2y_S_F2xz_S_M3_vrr+2*oned2z*I_ERI_Py_S_F2xz_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_F2xz_S_M3_vrr;
      Double I_ERI_F2yz_S_F2xz_S_M2_vrr = PAZ*I_ERI_D2y_S_F2xz_S_M2_vrr+WPZ*I_ERI_D2y_S_F2xz_S_M3_vrr+oned2k*I_ERI_D2y_S_D2x_S_M3_vrr;
      Double I_ERI_F3z_S_F2xz_S_M2_vrr = PAZ*I_ERI_D2z_S_F2xz_S_M2_vrr+WPZ*I_ERI_D2z_S_F2xz_S_M3_vrr+2*oned2z*I_ERI_Pz_S_F2xz_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xz_S_M3_vrr+oned2k*I_ERI_D2z_S_D2x_S_M3_vrr;
      Double I_ERI_F3x_S_Fx2y_S_M2_vrr = PAX*I_ERI_D2x_S_Fx2y_S_M2_vrr+WPX*I_ERI_D2x_S_Fx2y_S_M3_vrr+2*oned2z*I_ERI_Px_S_Fx2y_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2y_S_M3_vrr+oned2k*I_ERI_D2x_S_D2y_S_M3_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_M2_vrr = PAY*I_ERI_D2x_S_Fx2y_S_M2_vrr+WPY*I_ERI_D2x_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M3_vrr;
      Double I_ERI_F2xz_S_Fx2y_S_M2_vrr = PAZ*I_ERI_D2x_S_Fx2y_S_M2_vrr+WPZ*I_ERI_D2x_S_Fx2y_S_M3_vrr;
      Double I_ERI_F3y_S_Fx2y_S_M2_vrr = PAY*I_ERI_D2y_S_Fx2y_S_M2_vrr+WPY*I_ERI_D2y_S_Fx2y_S_M3_vrr+2*oned2z*I_ERI_Py_S_Fx2y_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M3_vrr;
      Double I_ERI_F2yz_S_Fx2y_S_M2_vrr = PAZ*I_ERI_D2y_S_Fx2y_S_M2_vrr+WPZ*I_ERI_D2y_S_Fx2y_S_M3_vrr;
      Double I_ERI_F3z_S_Fx2y_S_M2_vrr = PAZ*I_ERI_D2z_S_Fx2y_S_M2_vrr+WPZ*I_ERI_D2z_S_Fx2y_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Fx2y_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2y_S_M3_vrr;
      Double I_ERI_F3x_S_Fxyz_S_M2_vrr = PAX*I_ERI_D2x_S_Fxyz_S_M2_vrr+WPX*I_ERI_D2x_S_Fxyz_S_M3_vrr+2*oned2z*I_ERI_Px_S_Fxyz_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Fxyz_S_M3_vrr+oned2k*I_ERI_D2x_S_Dyz_S_M3_vrr;
      Double I_ERI_F2xy_S_Fxyz_S_M2_vrr = PAY*I_ERI_D2x_S_Fxyz_S_M2_vrr+WPY*I_ERI_D2x_S_Fxyz_S_M3_vrr+oned2k*I_ERI_D2x_S_Dxz_S_M3_vrr;
      Double I_ERI_F2xz_S_Fxyz_S_M2_vrr = PAZ*I_ERI_D2x_S_Fxyz_S_M2_vrr+WPZ*I_ERI_D2x_S_Fxyz_S_M3_vrr+oned2k*I_ERI_D2x_S_Dxy_S_M3_vrr;
      Double I_ERI_F3y_S_Fxyz_S_M2_vrr = PAY*I_ERI_D2y_S_Fxyz_S_M2_vrr+WPY*I_ERI_D2y_S_Fxyz_S_M3_vrr+2*oned2z*I_ERI_Py_S_Fxyz_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Fxyz_S_M3_vrr+oned2k*I_ERI_D2y_S_Dxz_S_M3_vrr;
      Double I_ERI_F2yz_S_Fxyz_S_M2_vrr = PAZ*I_ERI_D2y_S_Fxyz_S_M2_vrr+WPZ*I_ERI_D2y_S_Fxyz_S_M3_vrr+oned2k*I_ERI_D2y_S_Dxy_S_M3_vrr;
      Double I_ERI_F3z_S_Fxyz_S_M2_vrr = PAZ*I_ERI_D2z_S_Fxyz_S_M2_vrr+WPZ*I_ERI_D2z_S_Fxyz_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Fxyz_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Fxyz_S_M3_vrr+oned2k*I_ERI_D2z_S_Dxy_S_M3_vrr;
      Double I_ERI_F3x_S_Fx2z_S_M2_vrr = PAX*I_ERI_D2x_S_Fx2z_S_M2_vrr+WPX*I_ERI_D2x_S_Fx2z_S_M3_vrr+2*oned2z*I_ERI_Px_S_Fx2z_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2z_S_M3_vrr+oned2k*I_ERI_D2x_S_D2z_S_M3_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_M2_vrr = PAY*I_ERI_D2x_S_Fx2z_S_M2_vrr+WPY*I_ERI_D2x_S_Fx2z_S_M3_vrr;
      Double I_ERI_F2xz_S_Fx2z_S_M2_vrr = PAZ*I_ERI_D2x_S_Fx2z_S_M2_vrr+WPZ*I_ERI_D2x_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M3_vrr;
      Double I_ERI_F3y_S_Fx2z_S_M2_vrr = PAY*I_ERI_D2y_S_Fx2z_S_M2_vrr+WPY*I_ERI_D2y_S_Fx2z_S_M3_vrr+2*oned2z*I_ERI_Py_S_Fx2z_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2z_S_M3_vrr;
      Double I_ERI_F2yz_S_Fx2z_S_M2_vrr = PAZ*I_ERI_D2y_S_Fx2z_S_M2_vrr+WPZ*I_ERI_D2y_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M3_vrr;
      Double I_ERI_F3z_S_Fx2z_S_M2_vrr = PAZ*I_ERI_D2z_S_Fx2z_S_M2_vrr+WPZ*I_ERI_D2z_S_Fx2z_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Fx2z_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M3_vrr;
      Double I_ERI_F3x_S_F3y_S_M2_vrr = PAX*I_ERI_D2x_S_F3y_S_M2_vrr+WPX*I_ERI_D2x_S_F3y_S_M3_vrr+2*oned2z*I_ERI_Px_S_F3y_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_F3y_S_M3_vrr;
      Double I_ERI_F2xy_S_F3y_S_M2_vrr = PAY*I_ERI_D2x_S_F3y_S_M2_vrr+WPY*I_ERI_D2x_S_F3y_S_M3_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M3_vrr;
      Double I_ERI_F2xz_S_F3y_S_M2_vrr = PAZ*I_ERI_D2x_S_F3y_S_M2_vrr+WPZ*I_ERI_D2x_S_F3y_S_M3_vrr;
      Double I_ERI_F3y_S_F3y_S_M2_vrr = PAY*I_ERI_D2y_S_F3y_S_M2_vrr+WPY*I_ERI_D2y_S_F3y_S_M3_vrr+2*oned2z*I_ERI_Py_S_F3y_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_F3y_S_M3_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M3_vrr;
      Double I_ERI_F2yz_S_F3y_S_M2_vrr = PAZ*I_ERI_D2y_S_F3y_S_M2_vrr+WPZ*I_ERI_D2y_S_F3y_S_M3_vrr;
      Double I_ERI_F3z_S_F3y_S_M2_vrr = PAZ*I_ERI_D2z_S_F3y_S_M2_vrr+WPZ*I_ERI_D2z_S_F3y_S_M3_vrr+2*oned2z*I_ERI_Pz_S_F3y_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_F3y_S_M3_vrr;
      Double I_ERI_F3x_S_F2yz_S_M2_vrr = PAX*I_ERI_D2x_S_F2yz_S_M2_vrr+WPX*I_ERI_D2x_S_F2yz_S_M3_vrr+2*oned2z*I_ERI_Px_S_F2yz_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_F2yz_S_M3_vrr;
      Double I_ERI_F2xy_S_F2yz_S_M2_vrr = PAY*I_ERI_D2x_S_F2yz_S_M2_vrr+WPY*I_ERI_D2x_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M3_vrr;
      Double I_ERI_F2xz_S_F2yz_S_M2_vrr = PAZ*I_ERI_D2x_S_F2yz_S_M2_vrr+WPZ*I_ERI_D2x_S_F2yz_S_M3_vrr+oned2k*I_ERI_D2x_S_D2y_S_M3_vrr;
      Double I_ERI_F3y_S_F2yz_S_M2_vrr = PAY*I_ERI_D2y_S_F2yz_S_M2_vrr+WPY*I_ERI_D2y_S_F2yz_S_M3_vrr+2*oned2z*I_ERI_Py_S_F2yz_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M3_vrr;
      Double I_ERI_F2yz_S_F2yz_S_M2_vrr = PAZ*I_ERI_D2y_S_F2yz_S_M2_vrr+WPZ*I_ERI_D2y_S_F2yz_S_M3_vrr+oned2k*I_ERI_D2y_S_D2y_S_M3_vrr;
      Double I_ERI_F3z_S_F2yz_S_M2_vrr = PAZ*I_ERI_D2z_S_F2yz_S_M2_vrr+WPZ*I_ERI_D2z_S_F2yz_S_M3_vrr+2*oned2z*I_ERI_Pz_S_F2yz_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_F2yz_S_M3_vrr+oned2k*I_ERI_D2z_S_D2y_S_M3_vrr;
      Double I_ERI_F3x_S_Fy2z_S_M2_vrr = PAX*I_ERI_D2x_S_Fy2z_S_M2_vrr+WPX*I_ERI_D2x_S_Fy2z_S_M3_vrr+2*oned2z*I_ERI_Px_S_Fy2z_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_Fy2z_S_M3_vrr;
      Double I_ERI_F2xy_S_Fy2z_S_M2_vrr = PAY*I_ERI_D2x_S_Fy2z_S_M2_vrr+WPY*I_ERI_D2x_S_Fy2z_S_M3_vrr+oned2k*I_ERI_D2x_S_D2z_S_M3_vrr;
      Double I_ERI_F2xz_S_Fy2z_S_M2_vrr = PAZ*I_ERI_D2x_S_Fy2z_S_M2_vrr+WPZ*I_ERI_D2x_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M3_vrr;
      Double I_ERI_F3y_S_Fy2z_S_M2_vrr = PAY*I_ERI_D2y_S_Fy2z_S_M2_vrr+WPY*I_ERI_D2y_S_Fy2z_S_M3_vrr+2*oned2z*I_ERI_Py_S_Fy2z_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_Fy2z_S_M3_vrr+oned2k*I_ERI_D2y_S_D2z_S_M3_vrr;
      Double I_ERI_F2yz_S_Fy2z_S_M2_vrr = PAZ*I_ERI_D2y_S_Fy2z_S_M2_vrr+WPZ*I_ERI_D2y_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M3_vrr;
      Double I_ERI_F3z_S_Fy2z_S_M2_vrr = PAZ*I_ERI_D2z_S_Fy2z_S_M2_vrr+WPZ*I_ERI_D2z_S_Fy2z_S_M3_vrr+2*oned2z*I_ERI_Pz_S_Fy2z_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M3_vrr;
      Double I_ERI_F3x_S_F3z_S_M2_vrr = PAX*I_ERI_D2x_S_F3z_S_M2_vrr+WPX*I_ERI_D2x_S_F3z_S_M3_vrr+2*oned2z*I_ERI_Px_S_F3z_S_M2_vrr-2*rhod2zsq*I_ERI_Px_S_F3z_S_M3_vrr;
      Double I_ERI_F2xy_S_F3z_S_M2_vrr = PAY*I_ERI_D2x_S_F3z_S_M2_vrr+WPY*I_ERI_D2x_S_F3z_S_M3_vrr;
      Double I_ERI_F2xz_S_F3z_S_M2_vrr = PAZ*I_ERI_D2x_S_F3z_S_M2_vrr+WPZ*I_ERI_D2x_S_F3z_S_M3_vrr+3*oned2k*I_ERI_D2x_S_D2z_S_M3_vrr;
      Double I_ERI_F3y_S_F3z_S_M2_vrr = PAY*I_ERI_D2y_S_F3z_S_M2_vrr+WPY*I_ERI_D2y_S_F3z_S_M3_vrr+2*oned2z*I_ERI_Py_S_F3z_S_M2_vrr-2*rhod2zsq*I_ERI_Py_S_F3z_S_M3_vrr;
      Double I_ERI_F2yz_S_F3z_S_M2_vrr = PAZ*I_ERI_D2y_S_F3z_S_M2_vrr+WPZ*I_ERI_D2y_S_F3z_S_M3_vrr+3*oned2k*I_ERI_D2y_S_D2z_S_M3_vrr;
      Double I_ERI_F3z_S_F3z_S_M2_vrr = PAZ*I_ERI_D2z_S_F3z_S_M2_vrr+WPZ*I_ERI_D2z_S_F3z_S_M3_vrr+2*oned2z*I_ERI_Pz_S_F3z_S_M2_vrr-2*rhod2zsq*I_ERI_Pz_S_F3z_S_M3_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 45 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_M2_vrr = QCX*I_ERI_D2x_S_F3x_S_M2_vrr+WQX*I_ERI_D2x_S_F3x_S_M3_vrr+3*oned2e*I_ERI_D2x_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_D2x_S_D2x_S_M3_vrr+2*oned2k*I_ERI_Px_S_F3x_S_M3_vrr;
      Double I_ERI_D2y_S_G4x_S_M2_vrr = QCX*I_ERI_D2y_S_F3x_S_M2_vrr+WQX*I_ERI_D2y_S_F3x_S_M3_vrr+3*oned2e*I_ERI_D2y_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_D2y_S_D2x_S_M3_vrr;
      Double I_ERI_D2z_S_G4x_S_M2_vrr = QCX*I_ERI_D2z_S_F3x_S_M2_vrr+WQX*I_ERI_D2z_S_F3x_S_M3_vrr+3*oned2e*I_ERI_D2z_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_D2z_S_D2x_S_M3_vrr;
      Double I_ERI_D2x_S_G3xy_S_M2_vrr = QCY*I_ERI_D2x_S_F3x_S_M2_vrr+WQY*I_ERI_D2x_S_F3x_S_M3_vrr;
      Double I_ERI_D2y_S_G3xy_S_M2_vrr = QCY*I_ERI_D2y_S_F3x_S_M2_vrr+WQY*I_ERI_D2y_S_F3x_S_M3_vrr+2*oned2k*I_ERI_Py_S_F3x_S_M3_vrr;
      Double I_ERI_D2z_S_G3xy_S_M2_vrr = QCY*I_ERI_D2z_S_F3x_S_M2_vrr+WQY*I_ERI_D2z_S_F3x_S_M3_vrr;
      Double I_ERI_D2x_S_G3xz_S_M2_vrr = QCZ*I_ERI_D2x_S_F3x_S_M2_vrr+WQZ*I_ERI_D2x_S_F3x_S_M3_vrr;
      Double I_ERI_D2y_S_G3xz_S_M2_vrr = QCZ*I_ERI_D2y_S_F3x_S_M2_vrr+WQZ*I_ERI_D2y_S_F3x_S_M3_vrr;
      Double I_ERI_D2z_S_G3xz_S_M2_vrr = QCZ*I_ERI_D2z_S_F3x_S_M2_vrr+WQZ*I_ERI_D2z_S_F3x_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F3x_S_M3_vrr;
      Double I_ERI_D2x_S_G2x2y_S_M2_vrr = QCY*I_ERI_D2x_S_F2xy_S_M2_vrr+WQY*I_ERI_D2x_S_F2xy_S_M3_vrr+oned2e*I_ERI_D2x_S_D2x_S_M2_vrr-rhod2esq*I_ERI_D2x_S_D2x_S_M3_vrr;
      Double I_ERI_D2y_S_G2x2y_S_M2_vrr = QCY*I_ERI_D2y_S_F2xy_S_M2_vrr+WQY*I_ERI_D2y_S_F2xy_S_M3_vrr+oned2e*I_ERI_D2y_S_D2x_S_M2_vrr-rhod2esq*I_ERI_D2y_S_D2x_S_M3_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M3_vrr;
      Double I_ERI_D2z_S_G2x2y_S_M2_vrr = QCY*I_ERI_D2z_S_F2xy_S_M2_vrr+WQY*I_ERI_D2z_S_F2xy_S_M3_vrr+oned2e*I_ERI_D2z_S_D2x_S_M2_vrr-rhod2esq*I_ERI_D2z_S_D2x_S_M3_vrr;
      Double I_ERI_D2x_S_G2xyz_S_M2_vrr = QCZ*I_ERI_D2x_S_F2xy_S_M2_vrr+WQZ*I_ERI_D2x_S_F2xy_S_M3_vrr;
      Double I_ERI_D2y_S_G2xyz_S_M2_vrr = QCZ*I_ERI_D2y_S_F2xy_S_M2_vrr+WQZ*I_ERI_D2y_S_F2xy_S_M3_vrr;
      Double I_ERI_D2z_S_G2xyz_S_M2_vrr = QCZ*I_ERI_D2z_S_F2xy_S_M2_vrr+WQZ*I_ERI_D2z_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F2xy_S_M3_vrr;
      Double I_ERI_D2x_S_G2x2z_S_M2_vrr = QCZ*I_ERI_D2x_S_F2xz_S_M2_vrr+WQZ*I_ERI_D2x_S_F2xz_S_M3_vrr+oned2e*I_ERI_D2x_S_D2x_S_M2_vrr-rhod2esq*I_ERI_D2x_S_D2x_S_M3_vrr;
      Double I_ERI_D2y_S_G2x2z_S_M2_vrr = QCZ*I_ERI_D2y_S_F2xz_S_M2_vrr+WQZ*I_ERI_D2y_S_F2xz_S_M3_vrr+oned2e*I_ERI_D2y_S_D2x_S_M2_vrr-rhod2esq*I_ERI_D2y_S_D2x_S_M3_vrr;
      Double I_ERI_D2z_S_G2x2z_S_M2_vrr = QCZ*I_ERI_D2z_S_F2xz_S_M2_vrr+WQZ*I_ERI_D2z_S_F2xz_S_M3_vrr+oned2e*I_ERI_D2z_S_D2x_S_M2_vrr-rhod2esq*I_ERI_D2z_S_D2x_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M3_vrr;
      Double I_ERI_D2x_S_Gx3y_S_M2_vrr = QCX*I_ERI_D2x_S_F3y_S_M2_vrr+WQX*I_ERI_D2x_S_F3y_S_M3_vrr+2*oned2k*I_ERI_Px_S_F3y_S_M3_vrr;
      Double I_ERI_D2y_S_Gx3y_S_M2_vrr = QCX*I_ERI_D2y_S_F3y_S_M2_vrr+WQX*I_ERI_D2y_S_F3y_S_M3_vrr;
      Double I_ERI_D2z_S_Gx3y_S_M2_vrr = QCX*I_ERI_D2z_S_F3y_S_M2_vrr+WQX*I_ERI_D2z_S_F3y_S_M3_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_D2x_S_Fx2y_S_M2_vrr+WQZ*I_ERI_D2x_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_D2y_S_Fx2y_S_M2_vrr+WQZ*I_ERI_D2y_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_D2z_S_Fx2y_S_M2_vrr+WQZ*I_ERI_D2z_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_M2_vrr = QCY*I_ERI_D2x_S_Fx2z_S_M2_vrr+WQY*I_ERI_D2x_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_M2_vrr = QCY*I_ERI_D2y_S_Fx2z_S_M2_vrr+WQY*I_ERI_D2y_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_Py_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_M2_vrr = QCY*I_ERI_D2z_S_Fx2z_S_M2_vrr+WQY*I_ERI_D2z_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2x_S_Gx3z_S_M2_vrr = QCX*I_ERI_D2x_S_F3z_S_M2_vrr+WQX*I_ERI_D2x_S_F3z_S_M3_vrr+2*oned2k*I_ERI_Px_S_F3z_S_M3_vrr;
      Double I_ERI_D2y_S_Gx3z_S_M2_vrr = QCX*I_ERI_D2y_S_F3z_S_M2_vrr+WQX*I_ERI_D2y_S_F3z_S_M3_vrr;
      Double I_ERI_D2z_S_Gx3z_S_M2_vrr = QCX*I_ERI_D2z_S_F3z_S_M2_vrr+WQX*I_ERI_D2z_S_F3z_S_M3_vrr;
      Double I_ERI_D2x_S_G4y_S_M2_vrr = QCY*I_ERI_D2x_S_F3y_S_M2_vrr+WQY*I_ERI_D2x_S_F3y_S_M3_vrr+3*oned2e*I_ERI_D2x_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_D2x_S_D2y_S_M3_vrr;
      Double I_ERI_D2y_S_G4y_S_M2_vrr = QCY*I_ERI_D2y_S_F3y_S_M2_vrr+WQY*I_ERI_D2y_S_F3y_S_M3_vrr+3*oned2e*I_ERI_D2y_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_D2y_S_D2y_S_M3_vrr+2*oned2k*I_ERI_Py_S_F3y_S_M3_vrr;
      Double I_ERI_D2z_S_G4y_S_M2_vrr = QCY*I_ERI_D2z_S_F3y_S_M2_vrr+WQY*I_ERI_D2z_S_F3y_S_M3_vrr+3*oned2e*I_ERI_D2z_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_D2z_S_D2y_S_M3_vrr;
      Double I_ERI_D2x_S_G3yz_S_M2_vrr = QCZ*I_ERI_D2x_S_F3y_S_M2_vrr+WQZ*I_ERI_D2x_S_F3y_S_M3_vrr;
      Double I_ERI_D2y_S_G3yz_S_M2_vrr = QCZ*I_ERI_D2y_S_F3y_S_M2_vrr+WQZ*I_ERI_D2y_S_F3y_S_M3_vrr;
      Double I_ERI_D2z_S_G3yz_S_M2_vrr = QCZ*I_ERI_D2z_S_F3y_S_M2_vrr+WQZ*I_ERI_D2z_S_F3y_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F3y_S_M3_vrr;
      Double I_ERI_D2x_S_G2y2z_S_M2_vrr = QCZ*I_ERI_D2x_S_F2yz_S_M2_vrr+WQZ*I_ERI_D2x_S_F2yz_S_M3_vrr+oned2e*I_ERI_D2x_S_D2y_S_M2_vrr-rhod2esq*I_ERI_D2x_S_D2y_S_M3_vrr;
      Double I_ERI_D2y_S_G2y2z_S_M2_vrr = QCZ*I_ERI_D2y_S_F2yz_S_M2_vrr+WQZ*I_ERI_D2y_S_F2yz_S_M3_vrr+oned2e*I_ERI_D2y_S_D2y_S_M2_vrr-rhod2esq*I_ERI_D2y_S_D2y_S_M3_vrr;
      Double I_ERI_D2z_S_G2y2z_S_M2_vrr = QCZ*I_ERI_D2z_S_F2yz_S_M2_vrr+WQZ*I_ERI_D2z_S_F2yz_S_M3_vrr+oned2e*I_ERI_D2z_S_D2y_S_M2_vrr-rhod2esq*I_ERI_D2z_S_D2y_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M3_vrr;
      Double I_ERI_D2x_S_Gy3z_S_M2_vrr = QCY*I_ERI_D2x_S_F3z_S_M2_vrr+WQY*I_ERI_D2x_S_F3z_S_M3_vrr;
      Double I_ERI_D2y_S_Gy3z_S_M2_vrr = QCY*I_ERI_D2y_S_F3z_S_M2_vrr+WQY*I_ERI_D2y_S_F3z_S_M3_vrr+2*oned2k*I_ERI_Py_S_F3z_S_M3_vrr;
      Double I_ERI_D2z_S_Gy3z_S_M2_vrr = QCY*I_ERI_D2z_S_F3z_S_M2_vrr+WQY*I_ERI_D2z_S_F3z_S_M3_vrr;
      Double I_ERI_D2x_S_G4z_S_M2_vrr = QCZ*I_ERI_D2x_S_F3z_S_M2_vrr+WQZ*I_ERI_D2x_S_F3z_S_M3_vrr+3*oned2e*I_ERI_D2x_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_D2x_S_D2z_S_M3_vrr;
      Double I_ERI_D2y_S_G4z_S_M2_vrr = QCZ*I_ERI_D2y_S_F3z_S_M2_vrr+WQZ*I_ERI_D2y_S_F3z_S_M3_vrr+3*oned2e*I_ERI_D2y_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_D2y_S_D2z_S_M3_vrr;
      Double I_ERI_D2z_S_G4z_S_M2_vrr = QCZ*I_ERI_D2z_S_F3z_S_M2_vrr+WQZ*I_ERI_D2z_S_F3z_S_M3_vrr+3*oned2e*I_ERI_D2z_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_D2z_S_D2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F3z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M3
       ************************************************************/
      Double I_ERI_G4x_S_D2x_S_M2_vrr = PAX*I_ERI_F3x_S_D2x_S_M2_vrr+WPX*I_ERI_F3x_S_D2x_S_M3_vrr+3*oned2z*I_ERI_D2x_S_D2x_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_D2x_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Px_S_M3_vrr;
      Double I_ERI_G3xy_S_D2x_S_M2_vrr = PAY*I_ERI_F3x_S_D2x_S_M2_vrr+WPY*I_ERI_F3x_S_D2x_S_M3_vrr;
      Double I_ERI_G3xz_S_D2x_S_M2_vrr = PAZ*I_ERI_F3x_S_D2x_S_M2_vrr+WPZ*I_ERI_F3x_S_D2x_S_M3_vrr;
      Double I_ERI_G2x2y_S_D2x_S_M2_vrr = PAY*I_ERI_F2xy_S_D2x_S_M2_vrr+WPY*I_ERI_F2xy_S_D2x_S_M3_vrr+oned2z*I_ERI_D2x_S_D2x_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_D2x_S_M3_vrr;
      Double I_ERI_Gx3y_S_D2x_S_M2_vrr = PAX*I_ERI_F3y_S_D2x_S_M2_vrr+WPX*I_ERI_F3y_S_D2x_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Px_S_M3_vrr;
      Double I_ERI_Gx3z_S_D2x_S_M2_vrr = PAX*I_ERI_F3z_S_D2x_S_M2_vrr+WPX*I_ERI_F3z_S_D2x_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Px_S_M3_vrr;
      Double I_ERI_G4y_S_D2x_S_M2_vrr = PAY*I_ERI_F3y_S_D2x_S_M2_vrr+WPY*I_ERI_F3y_S_D2x_S_M3_vrr+3*oned2z*I_ERI_D2y_S_D2x_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_D2x_S_M3_vrr;
      Double I_ERI_G3yz_S_D2x_S_M2_vrr = PAZ*I_ERI_F3y_S_D2x_S_M2_vrr+WPZ*I_ERI_F3y_S_D2x_S_M3_vrr;
      Double I_ERI_Gy3z_S_D2x_S_M2_vrr = PAY*I_ERI_F3z_S_D2x_S_M2_vrr+WPY*I_ERI_F3z_S_D2x_S_M3_vrr;
      Double I_ERI_G4z_S_D2x_S_M2_vrr = PAZ*I_ERI_F3z_S_D2x_S_M2_vrr+WPZ*I_ERI_F3z_S_D2x_S_M3_vrr+3*oned2z*I_ERI_D2z_S_D2x_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_D2x_S_M3_vrr;
      Double I_ERI_G4x_S_Dxy_S_M2_vrr = PAX*I_ERI_F3x_S_Dxy_S_M2_vrr+WPX*I_ERI_F3x_S_Dxy_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Dxy_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Dxy_S_M3_vrr+oned2k*I_ERI_F3x_S_Py_S_M3_vrr;
      Double I_ERI_G3xy_S_Dxy_S_M2_vrr = PAY*I_ERI_F3x_S_Dxy_S_M2_vrr+WPY*I_ERI_F3x_S_Dxy_S_M3_vrr+oned2k*I_ERI_F3x_S_Px_S_M3_vrr;
      Double I_ERI_G3xz_S_Dxy_S_M2_vrr = PAZ*I_ERI_F3x_S_Dxy_S_M2_vrr+WPZ*I_ERI_F3x_S_Dxy_S_M3_vrr;
      Double I_ERI_G2x2y_S_Dxy_S_M2_vrr = PAY*I_ERI_F2xy_S_Dxy_S_M2_vrr+WPY*I_ERI_F2xy_S_Dxy_S_M3_vrr+oned2z*I_ERI_D2x_S_Dxy_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Dxy_S_M3_vrr+oned2k*I_ERI_F2xy_S_Px_S_M3_vrr;
      Double I_ERI_Gx3y_S_Dxy_S_M2_vrr = PAX*I_ERI_F3y_S_Dxy_S_M2_vrr+WPX*I_ERI_F3y_S_Dxy_S_M3_vrr+oned2k*I_ERI_F3y_S_Py_S_M3_vrr;
      Double I_ERI_Gx3z_S_Dxy_S_M2_vrr = PAX*I_ERI_F3z_S_Dxy_S_M2_vrr+WPX*I_ERI_F3z_S_Dxy_S_M3_vrr+oned2k*I_ERI_F3z_S_Py_S_M3_vrr;
      Double I_ERI_G4y_S_Dxy_S_M2_vrr = PAY*I_ERI_F3y_S_Dxy_S_M2_vrr+WPY*I_ERI_F3y_S_Dxy_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Dxy_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Dxy_S_M3_vrr+oned2k*I_ERI_F3y_S_Px_S_M3_vrr;
      Double I_ERI_G3yz_S_Dxy_S_M2_vrr = PAZ*I_ERI_F3y_S_Dxy_S_M2_vrr+WPZ*I_ERI_F3y_S_Dxy_S_M3_vrr;
      Double I_ERI_Gy3z_S_Dxy_S_M2_vrr = PAY*I_ERI_F3z_S_Dxy_S_M2_vrr+WPY*I_ERI_F3z_S_Dxy_S_M3_vrr+oned2k*I_ERI_F3z_S_Px_S_M3_vrr;
      Double I_ERI_G4z_S_Dxy_S_M2_vrr = PAZ*I_ERI_F3z_S_Dxy_S_M2_vrr+WPZ*I_ERI_F3z_S_Dxy_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Dxy_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Dxy_S_M3_vrr;
      Double I_ERI_G4x_S_Dxz_S_M2_vrr = PAX*I_ERI_F3x_S_Dxz_S_M2_vrr+WPX*I_ERI_F3x_S_Dxz_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Dxz_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Dxz_S_M3_vrr+oned2k*I_ERI_F3x_S_Pz_S_M3_vrr;
      Double I_ERI_G3xy_S_Dxz_S_M2_vrr = PAY*I_ERI_F3x_S_Dxz_S_M2_vrr+WPY*I_ERI_F3x_S_Dxz_S_M3_vrr;
      Double I_ERI_G3xz_S_Dxz_S_M2_vrr = PAZ*I_ERI_F3x_S_Dxz_S_M2_vrr+WPZ*I_ERI_F3x_S_Dxz_S_M3_vrr+oned2k*I_ERI_F3x_S_Px_S_M3_vrr;
      Double I_ERI_G2x2y_S_Dxz_S_M2_vrr = PAY*I_ERI_F2xy_S_Dxz_S_M2_vrr+WPY*I_ERI_F2xy_S_Dxz_S_M3_vrr+oned2z*I_ERI_D2x_S_Dxz_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Dxz_S_M3_vrr;
      Double I_ERI_Gx3y_S_Dxz_S_M2_vrr = PAX*I_ERI_F3y_S_Dxz_S_M2_vrr+WPX*I_ERI_F3y_S_Dxz_S_M3_vrr+oned2k*I_ERI_F3y_S_Pz_S_M3_vrr;
      Double I_ERI_Gx3z_S_Dxz_S_M2_vrr = PAX*I_ERI_F3z_S_Dxz_S_M2_vrr+WPX*I_ERI_F3z_S_Dxz_S_M3_vrr+oned2k*I_ERI_F3z_S_Pz_S_M3_vrr;
      Double I_ERI_G4y_S_Dxz_S_M2_vrr = PAY*I_ERI_F3y_S_Dxz_S_M2_vrr+WPY*I_ERI_F3y_S_Dxz_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Dxz_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Dxz_S_M3_vrr;
      Double I_ERI_G3yz_S_Dxz_S_M2_vrr = PAZ*I_ERI_F3y_S_Dxz_S_M2_vrr+WPZ*I_ERI_F3y_S_Dxz_S_M3_vrr+oned2k*I_ERI_F3y_S_Px_S_M3_vrr;
      Double I_ERI_Gy3z_S_Dxz_S_M2_vrr = PAY*I_ERI_F3z_S_Dxz_S_M2_vrr+WPY*I_ERI_F3z_S_Dxz_S_M3_vrr;
      Double I_ERI_G4z_S_Dxz_S_M2_vrr = PAZ*I_ERI_F3z_S_Dxz_S_M2_vrr+WPZ*I_ERI_F3z_S_Dxz_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Dxz_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Dxz_S_M3_vrr+oned2k*I_ERI_F3z_S_Px_S_M3_vrr;
      Double I_ERI_G4x_S_D2y_S_M2_vrr = PAX*I_ERI_F3x_S_D2y_S_M2_vrr+WPX*I_ERI_F3x_S_D2y_S_M3_vrr+3*oned2z*I_ERI_D2x_S_D2y_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_D2y_S_M3_vrr;
      Double I_ERI_G3xy_S_D2y_S_M2_vrr = PAY*I_ERI_F3x_S_D2y_S_M2_vrr+WPY*I_ERI_F3x_S_D2y_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Py_S_M3_vrr;
      Double I_ERI_G3xz_S_D2y_S_M2_vrr = PAZ*I_ERI_F3x_S_D2y_S_M2_vrr+WPZ*I_ERI_F3x_S_D2y_S_M3_vrr;
      Double I_ERI_G2x2y_S_D2y_S_M2_vrr = PAY*I_ERI_F2xy_S_D2y_S_M2_vrr+WPY*I_ERI_F2xy_S_D2y_S_M3_vrr+oned2z*I_ERI_D2x_S_D2y_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_D2y_S_M3_vrr+2*oned2k*I_ERI_F2xy_S_Py_S_M3_vrr;
      Double I_ERI_Gx3y_S_D2y_S_M2_vrr = PAX*I_ERI_F3y_S_D2y_S_M2_vrr+WPX*I_ERI_F3y_S_D2y_S_M3_vrr;
      Double I_ERI_Gx3z_S_D2y_S_M2_vrr = PAX*I_ERI_F3z_S_D2y_S_M2_vrr+WPX*I_ERI_F3z_S_D2y_S_M3_vrr;
      Double I_ERI_G4y_S_D2y_S_M2_vrr = PAY*I_ERI_F3y_S_D2y_S_M2_vrr+WPY*I_ERI_F3y_S_D2y_S_M3_vrr+3*oned2z*I_ERI_D2y_S_D2y_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_D2y_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Py_S_M3_vrr;
      Double I_ERI_G3yz_S_D2y_S_M2_vrr = PAZ*I_ERI_F3y_S_D2y_S_M2_vrr+WPZ*I_ERI_F3y_S_D2y_S_M3_vrr;
      Double I_ERI_Gy3z_S_D2y_S_M2_vrr = PAY*I_ERI_F3z_S_D2y_S_M2_vrr+WPY*I_ERI_F3z_S_D2y_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Py_S_M3_vrr;
      Double I_ERI_G4z_S_D2y_S_M2_vrr = PAZ*I_ERI_F3z_S_D2y_S_M2_vrr+WPZ*I_ERI_F3z_S_D2y_S_M3_vrr+3*oned2z*I_ERI_D2z_S_D2y_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_D2y_S_M3_vrr;
      Double I_ERI_G4x_S_Dyz_S_M2_vrr = PAX*I_ERI_F3x_S_Dyz_S_M2_vrr+WPX*I_ERI_F3x_S_Dyz_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Dyz_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Dyz_S_M3_vrr;
      Double I_ERI_G3xy_S_Dyz_S_M2_vrr = PAY*I_ERI_F3x_S_Dyz_S_M2_vrr+WPY*I_ERI_F3x_S_Dyz_S_M3_vrr+oned2k*I_ERI_F3x_S_Pz_S_M3_vrr;
      Double I_ERI_G3xz_S_Dyz_S_M2_vrr = PAZ*I_ERI_F3x_S_Dyz_S_M2_vrr+WPZ*I_ERI_F3x_S_Dyz_S_M3_vrr+oned2k*I_ERI_F3x_S_Py_S_M3_vrr;
      Double I_ERI_G2x2y_S_Dyz_S_M2_vrr = PAY*I_ERI_F2xy_S_Dyz_S_M2_vrr+WPY*I_ERI_F2xy_S_Dyz_S_M3_vrr+oned2z*I_ERI_D2x_S_Dyz_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Dyz_S_M3_vrr+oned2k*I_ERI_F2xy_S_Pz_S_M3_vrr;
      Double I_ERI_Gx3y_S_Dyz_S_M2_vrr = PAX*I_ERI_F3y_S_Dyz_S_M2_vrr+WPX*I_ERI_F3y_S_Dyz_S_M3_vrr;
      Double I_ERI_Gx3z_S_Dyz_S_M2_vrr = PAX*I_ERI_F3z_S_Dyz_S_M2_vrr+WPX*I_ERI_F3z_S_Dyz_S_M3_vrr;
      Double I_ERI_G4y_S_Dyz_S_M2_vrr = PAY*I_ERI_F3y_S_Dyz_S_M2_vrr+WPY*I_ERI_F3y_S_Dyz_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Dyz_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Dyz_S_M3_vrr+oned2k*I_ERI_F3y_S_Pz_S_M3_vrr;
      Double I_ERI_G3yz_S_Dyz_S_M2_vrr = PAZ*I_ERI_F3y_S_Dyz_S_M2_vrr+WPZ*I_ERI_F3y_S_Dyz_S_M3_vrr+oned2k*I_ERI_F3y_S_Py_S_M3_vrr;
      Double I_ERI_Gy3z_S_Dyz_S_M2_vrr = PAY*I_ERI_F3z_S_Dyz_S_M2_vrr+WPY*I_ERI_F3z_S_Dyz_S_M3_vrr+oned2k*I_ERI_F3z_S_Pz_S_M3_vrr;
      Double I_ERI_G4z_S_Dyz_S_M2_vrr = PAZ*I_ERI_F3z_S_Dyz_S_M2_vrr+WPZ*I_ERI_F3z_S_Dyz_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Dyz_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Dyz_S_M3_vrr+oned2k*I_ERI_F3z_S_Py_S_M3_vrr;
      Double I_ERI_G4x_S_D2z_S_M2_vrr = PAX*I_ERI_F3x_S_D2z_S_M2_vrr+WPX*I_ERI_F3x_S_D2z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_D2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_D2z_S_M3_vrr;
      Double I_ERI_G3xy_S_D2z_S_M2_vrr = PAY*I_ERI_F3x_S_D2z_S_M2_vrr+WPY*I_ERI_F3x_S_D2z_S_M3_vrr;
      Double I_ERI_G3xz_S_D2z_S_M2_vrr = PAZ*I_ERI_F3x_S_D2z_S_M2_vrr+WPZ*I_ERI_F3x_S_D2z_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Pz_S_M3_vrr;
      Double I_ERI_G2x2y_S_D2z_S_M2_vrr = PAY*I_ERI_F2xy_S_D2z_S_M2_vrr+WPY*I_ERI_F2xy_S_D2z_S_M3_vrr+oned2z*I_ERI_D2x_S_D2z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_D2z_S_M3_vrr;
      Double I_ERI_Gx3y_S_D2z_S_M2_vrr = PAX*I_ERI_F3y_S_D2z_S_M2_vrr+WPX*I_ERI_F3y_S_D2z_S_M3_vrr;
      Double I_ERI_Gx3z_S_D2z_S_M2_vrr = PAX*I_ERI_F3z_S_D2z_S_M2_vrr+WPX*I_ERI_F3z_S_D2z_S_M3_vrr;
      Double I_ERI_G4y_S_D2z_S_M2_vrr = PAY*I_ERI_F3y_S_D2z_S_M2_vrr+WPY*I_ERI_F3y_S_D2z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_D2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_D2z_S_M3_vrr;
      Double I_ERI_G3yz_S_D2z_S_M2_vrr = PAZ*I_ERI_F3y_S_D2z_S_M2_vrr+WPZ*I_ERI_F3y_S_D2z_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Pz_S_M3_vrr;
      Double I_ERI_Gy3z_S_D2z_S_M2_vrr = PAY*I_ERI_F3z_S_D2z_S_M2_vrr+WPY*I_ERI_F3z_S_D2z_S_M3_vrr;
      Double I_ERI_G4z_S_D2z_S_M2_vrr = PAZ*I_ERI_F3z_S_D2z_S_M2_vrr+WPZ*I_ERI_F3z_S_D2z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_D2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_D2z_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 60 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M3
       ************************************************************/
      Double I_ERI_F3x_S_G4x_S_M2_vrr = QCX*I_ERI_F3x_S_F3x_S_M2_vrr+WQX*I_ERI_F3x_S_F3x_S_M3_vrr+3*oned2e*I_ERI_F3x_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_F3x_S_D2x_S_M3_vrr+3*oned2k*I_ERI_D2x_S_F3x_S_M3_vrr;
      Double I_ERI_F2xy_S_G4x_S_M2_vrr = QCX*I_ERI_F2xy_S_F3x_S_M2_vrr+WQX*I_ERI_F2xy_S_F3x_S_M3_vrr+3*oned2e*I_ERI_F2xy_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_F2xy_S_D2x_S_M3_vrr+2*oned2k*I_ERI_Dxy_S_F3x_S_M3_vrr;
      Double I_ERI_F2xz_S_G4x_S_M2_vrr = QCX*I_ERI_F2xz_S_F3x_S_M2_vrr+WQX*I_ERI_F2xz_S_F3x_S_M3_vrr+3*oned2e*I_ERI_F2xz_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_F2xz_S_D2x_S_M3_vrr+2*oned2k*I_ERI_Dxz_S_F3x_S_M3_vrr;
      Double I_ERI_F3y_S_G4x_S_M2_vrr = QCX*I_ERI_F3y_S_F3x_S_M2_vrr+WQX*I_ERI_F3y_S_F3x_S_M3_vrr+3*oned2e*I_ERI_F3y_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_F3y_S_D2x_S_M3_vrr;
      Double I_ERI_F2yz_S_G4x_S_M2_vrr = QCX*I_ERI_F2yz_S_F3x_S_M2_vrr+WQX*I_ERI_F2yz_S_F3x_S_M3_vrr+3*oned2e*I_ERI_F2yz_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_F2yz_S_D2x_S_M3_vrr;
      Double I_ERI_F3z_S_G4x_S_M2_vrr = QCX*I_ERI_F3z_S_F3x_S_M2_vrr+WQX*I_ERI_F3z_S_F3x_S_M3_vrr+3*oned2e*I_ERI_F3z_S_D2x_S_M2_vrr-3*rhod2esq*I_ERI_F3z_S_D2x_S_M3_vrr;
      Double I_ERI_F3x_S_G3xy_S_M2_vrr = QCY*I_ERI_F3x_S_F3x_S_M2_vrr+WQY*I_ERI_F3x_S_F3x_S_M3_vrr;
      Double I_ERI_F2xy_S_G3xy_S_M2_vrr = QCY*I_ERI_F2xy_S_F3x_S_M2_vrr+WQY*I_ERI_F2xy_S_F3x_S_M3_vrr+oned2k*I_ERI_D2x_S_F3x_S_M3_vrr;
      Double I_ERI_F2xz_S_G3xy_S_M2_vrr = QCY*I_ERI_F2xz_S_F3x_S_M2_vrr+WQY*I_ERI_F2xz_S_F3x_S_M3_vrr;
      Double I_ERI_F3y_S_G3xy_S_M2_vrr = QCY*I_ERI_F3y_S_F3x_S_M2_vrr+WQY*I_ERI_F3y_S_F3x_S_M3_vrr+3*oned2k*I_ERI_D2y_S_F3x_S_M3_vrr;
      Double I_ERI_F2yz_S_G3xy_S_M2_vrr = QCY*I_ERI_F2yz_S_F3x_S_M2_vrr+WQY*I_ERI_F2yz_S_F3x_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_F3x_S_M3_vrr;
      Double I_ERI_F3z_S_G3xy_S_M2_vrr = QCY*I_ERI_F3z_S_F3x_S_M2_vrr+WQY*I_ERI_F3z_S_F3x_S_M3_vrr;
      Double I_ERI_F3x_S_G3xz_S_M2_vrr = QCZ*I_ERI_F3x_S_F3x_S_M2_vrr+WQZ*I_ERI_F3x_S_F3x_S_M3_vrr;
      Double I_ERI_F2xy_S_G3xz_S_M2_vrr = QCZ*I_ERI_F2xy_S_F3x_S_M2_vrr+WQZ*I_ERI_F2xy_S_F3x_S_M3_vrr;
      Double I_ERI_F2xz_S_G3xz_S_M2_vrr = QCZ*I_ERI_F2xz_S_F3x_S_M2_vrr+WQZ*I_ERI_F2xz_S_F3x_S_M3_vrr+oned2k*I_ERI_D2x_S_F3x_S_M3_vrr;
      Double I_ERI_F3y_S_G3xz_S_M2_vrr = QCZ*I_ERI_F3y_S_F3x_S_M2_vrr+WQZ*I_ERI_F3y_S_F3x_S_M3_vrr;
      Double I_ERI_F2yz_S_G3xz_S_M2_vrr = QCZ*I_ERI_F2yz_S_F3x_S_M2_vrr+WQZ*I_ERI_F2yz_S_F3x_S_M3_vrr+oned2k*I_ERI_D2y_S_F3x_S_M3_vrr;
      Double I_ERI_F3z_S_G3xz_S_M2_vrr = QCZ*I_ERI_F3z_S_F3x_S_M2_vrr+WQZ*I_ERI_F3z_S_F3x_S_M3_vrr+3*oned2k*I_ERI_D2z_S_F3x_S_M3_vrr;
      Double I_ERI_F3x_S_G2x2y_S_M2_vrr = QCY*I_ERI_F3x_S_F2xy_S_M2_vrr+WQY*I_ERI_F3x_S_F2xy_S_M3_vrr+oned2e*I_ERI_F3x_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F3x_S_D2x_S_M3_vrr;
      Double I_ERI_F2xy_S_G2x2y_S_M2_vrr = QCY*I_ERI_F2xy_S_F2xy_S_M2_vrr+WQY*I_ERI_F2xy_S_F2xy_S_M3_vrr+oned2e*I_ERI_F2xy_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F2xy_S_D2x_S_M3_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M3_vrr;
      Double I_ERI_F2xz_S_G2x2y_S_M2_vrr = QCY*I_ERI_F2xz_S_F2xy_S_M2_vrr+WQY*I_ERI_F2xz_S_F2xy_S_M3_vrr+oned2e*I_ERI_F2xz_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F2xz_S_D2x_S_M3_vrr;
      Double I_ERI_F3y_S_G2x2y_S_M2_vrr = QCY*I_ERI_F3y_S_F2xy_S_M2_vrr+WQY*I_ERI_F3y_S_F2xy_S_M3_vrr+oned2e*I_ERI_F3y_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F3y_S_D2x_S_M3_vrr+3*oned2k*I_ERI_D2y_S_F2xy_S_M3_vrr;
      Double I_ERI_F2yz_S_G2x2y_S_M2_vrr = QCY*I_ERI_F2yz_S_F2xy_S_M2_vrr+WQY*I_ERI_F2yz_S_F2xy_S_M3_vrr+oned2e*I_ERI_F2yz_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F2yz_S_D2x_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_F2xy_S_M3_vrr;
      Double I_ERI_F3z_S_G2x2y_S_M2_vrr = QCY*I_ERI_F3z_S_F2xy_S_M2_vrr+WQY*I_ERI_F3z_S_F2xy_S_M3_vrr+oned2e*I_ERI_F3z_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F3z_S_D2x_S_M3_vrr;
      Double I_ERI_F3x_S_G2xyz_S_M2_vrr = QCZ*I_ERI_F3x_S_F2xy_S_M2_vrr+WQZ*I_ERI_F3x_S_F2xy_S_M3_vrr;
      Double I_ERI_F2xy_S_G2xyz_S_M2_vrr = QCZ*I_ERI_F2xy_S_F2xy_S_M2_vrr+WQZ*I_ERI_F2xy_S_F2xy_S_M3_vrr;
      Double I_ERI_F2xz_S_G2xyz_S_M2_vrr = QCZ*I_ERI_F2xz_S_F2xy_S_M2_vrr+WQZ*I_ERI_F2xz_S_F2xy_S_M3_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M3_vrr;
      Double I_ERI_F3y_S_G2xyz_S_M2_vrr = QCZ*I_ERI_F3y_S_F2xy_S_M2_vrr+WQZ*I_ERI_F3y_S_F2xy_S_M3_vrr;
      Double I_ERI_F2yz_S_G2xyz_S_M2_vrr = QCZ*I_ERI_F2yz_S_F2xy_S_M2_vrr+WQZ*I_ERI_F2yz_S_F2xy_S_M3_vrr+oned2k*I_ERI_D2y_S_F2xy_S_M3_vrr;
      Double I_ERI_F3z_S_G2xyz_S_M2_vrr = QCZ*I_ERI_F3z_S_F2xy_S_M2_vrr+WQZ*I_ERI_F3z_S_F2xy_S_M3_vrr+3*oned2k*I_ERI_D2z_S_F2xy_S_M3_vrr;
      Double I_ERI_F3x_S_G2x2z_S_M2_vrr = QCZ*I_ERI_F3x_S_F2xz_S_M2_vrr+WQZ*I_ERI_F3x_S_F2xz_S_M3_vrr+oned2e*I_ERI_F3x_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F3x_S_D2x_S_M3_vrr;
      Double I_ERI_F2xy_S_G2x2z_S_M2_vrr = QCZ*I_ERI_F2xy_S_F2xz_S_M2_vrr+WQZ*I_ERI_F2xy_S_F2xz_S_M3_vrr+oned2e*I_ERI_F2xy_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F2xy_S_D2x_S_M3_vrr;
      Double I_ERI_F2xz_S_G2x2z_S_M2_vrr = QCZ*I_ERI_F2xz_S_F2xz_S_M2_vrr+WQZ*I_ERI_F2xz_S_F2xz_S_M3_vrr+oned2e*I_ERI_F2xz_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F2xz_S_D2x_S_M3_vrr+oned2k*I_ERI_D2x_S_F2xz_S_M3_vrr;
      Double I_ERI_F3y_S_G2x2z_S_M2_vrr = QCZ*I_ERI_F3y_S_F2xz_S_M2_vrr+WQZ*I_ERI_F3y_S_F2xz_S_M3_vrr+oned2e*I_ERI_F3y_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F3y_S_D2x_S_M3_vrr;
      Double I_ERI_F2yz_S_G2x2z_S_M2_vrr = QCZ*I_ERI_F2yz_S_F2xz_S_M2_vrr+WQZ*I_ERI_F2yz_S_F2xz_S_M3_vrr+oned2e*I_ERI_F2yz_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F2yz_S_D2x_S_M3_vrr+oned2k*I_ERI_D2y_S_F2xz_S_M3_vrr;
      Double I_ERI_F3z_S_G2x2z_S_M2_vrr = QCZ*I_ERI_F3z_S_F2xz_S_M2_vrr+WQZ*I_ERI_F3z_S_F2xz_S_M3_vrr+oned2e*I_ERI_F3z_S_D2x_S_M2_vrr-rhod2esq*I_ERI_F3z_S_D2x_S_M3_vrr+3*oned2k*I_ERI_D2z_S_F2xz_S_M3_vrr;
      Double I_ERI_F3x_S_Gx3y_S_M2_vrr = QCX*I_ERI_F3x_S_F3y_S_M2_vrr+WQX*I_ERI_F3x_S_F3y_S_M3_vrr+3*oned2k*I_ERI_D2x_S_F3y_S_M3_vrr;
      Double I_ERI_F2xy_S_Gx3y_S_M2_vrr = QCX*I_ERI_F2xy_S_F3y_S_M2_vrr+WQX*I_ERI_F2xy_S_F3y_S_M3_vrr+2*oned2k*I_ERI_Dxy_S_F3y_S_M3_vrr;
      Double I_ERI_F2xz_S_Gx3y_S_M2_vrr = QCX*I_ERI_F2xz_S_F3y_S_M2_vrr+WQX*I_ERI_F2xz_S_F3y_S_M3_vrr+2*oned2k*I_ERI_Dxz_S_F3y_S_M3_vrr;
      Double I_ERI_F3y_S_Gx3y_S_M2_vrr = QCX*I_ERI_F3y_S_F3y_S_M2_vrr+WQX*I_ERI_F3y_S_F3y_S_M3_vrr;
      Double I_ERI_F2yz_S_Gx3y_S_M2_vrr = QCX*I_ERI_F2yz_S_F3y_S_M2_vrr+WQX*I_ERI_F2yz_S_F3y_S_M3_vrr;
      Double I_ERI_F3z_S_Gx3y_S_M2_vrr = QCX*I_ERI_F3z_S_F3y_S_M2_vrr+WQX*I_ERI_F3z_S_F3y_S_M3_vrr;
      Double I_ERI_F3x_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_F3x_S_Fx2y_S_M2_vrr+WQZ*I_ERI_F3x_S_Fx2y_S_M3_vrr;
      Double I_ERI_F2xy_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_F2xy_S_Fx2y_S_M2_vrr+WQZ*I_ERI_F2xy_S_Fx2y_S_M3_vrr;
      Double I_ERI_F2xz_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_F2xz_S_Fx2y_S_M2_vrr+WQZ*I_ERI_F2xz_S_Fx2y_S_M3_vrr+oned2k*I_ERI_D2x_S_Fx2y_S_M3_vrr;
      Double I_ERI_F3y_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_F3y_S_Fx2y_S_M2_vrr+WQZ*I_ERI_F3y_S_Fx2y_S_M3_vrr;
      Double I_ERI_F2yz_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_F2yz_S_Fx2y_S_M2_vrr+WQZ*I_ERI_F2yz_S_Fx2y_S_M3_vrr+oned2k*I_ERI_D2y_S_Fx2y_S_M3_vrr;
      Double I_ERI_F3z_S_Gx2yz_S_M2_vrr = QCZ*I_ERI_F3z_S_Fx2y_S_M2_vrr+WQZ*I_ERI_F3z_S_Fx2y_S_M3_vrr+3*oned2k*I_ERI_D2z_S_Fx2y_S_M3_vrr;
      Double I_ERI_F3x_S_Gxy2z_S_M2_vrr = QCY*I_ERI_F3x_S_Fx2z_S_M2_vrr+WQY*I_ERI_F3x_S_Fx2z_S_M3_vrr;
      Double I_ERI_F2xy_S_Gxy2z_S_M2_vrr = QCY*I_ERI_F2xy_S_Fx2z_S_M2_vrr+WQY*I_ERI_F2xy_S_Fx2z_S_M3_vrr+oned2k*I_ERI_D2x_S_Fx2z_S_M3_vrr;
      Double I_ERI_F2xz_S_Gxy2z_S_M2_vrr = QCY*I_ERI_F2xz_S_Fx2z_S_M2_vrr+WQY*I_ERI_F2xz_S_Fx2z_S_M3_vrr;
      Double I_ERI_F3y_S_Gxy2z_S_M2_vrr = QCY*I_ERI_F3y_S_Fx2z_S_M2_vrr+WQY*I_ERI_F3y_S_Fx2z_S_M3_vrr+3*oned2k*I_ERI_D2y_S_Fx2z_S_M3_vrr;
      Double I_ERI_F2yz_S_Gxy2z_S_M2_vrr = QCY*I_ERI_F2yz_S_Fx2z_S_M2_vrr+WQY*I_ERI_F2yz_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_Fx2z_S_M3_vrr;
      Double I_ERI_F3z_S_Gxy2z_S_M2_vrr = QCY*I_ERI_F3z_S_Fx2z_S_M2_vrr+WQY*I_ERI_F3z_S_Fx2z_S_M3_vrr;
      Double I_ERI_F3x_S_Gx3z_S_M2_vrr = QCX*I_ERI_F3x_S_F3z_S_M2_vrr+WQX*I_ERI_F3x_S_F3z_S_M3_vrr+3*oned2k*I_ERI_D2x_S_F3z_S_M3_vrr;
      Double I_ERI_F2xy_S_Gx3z_S_M2_vrr = QCX*I_ERI_F2xy_S_F3z_S_M2_vrr+WQX*I_ERI_F2xy_S_F3z_S_M3_vrr+2*oned2k*I_ERI_Dxy_S_F3z_S_M3_vrr;
      Double I_ERI_F2xz_S_Gx3z_S_M2_vrr = QCX*I_ERI_F2xz_S_F3z_S_M2_vrr+WQX*I_ERI_F2xz_S_F3z_S_M3_vrr+2*oned2k*I_ERI_Dxz_S_F3z_S_M3_vrr;
      Double I_ERI_F3y_S_Gx3z_S_M2_vrr = QCX*I_ERI_F3y_S_F3z_S_M2_vrr+WQX*I_ERI_F3y_S_F3z_S_M3_vrr;
      Double I_ERI_F2yz_S_Gx3z_S_M2_vrr = QCX*I_ERI_F2yz_S_F3z_S_M2_vrr+WQX*I_ERI_F2yz_S_F3z_S_M3_vrr;
      Double I_ERI_F3z_S_Gx3z_S_M2_vrr = QCX*I_ERI_F3z_S_F3z_S_M2_vrr+WQX*I_ERI_F3z_S_F3z_S_M3_vrr;
      Double I_ERI_F3x_S_G4y_S_M2_vrr = QCY*I_ERI_F3x_S_F3y_S_M2_vrr+WQY*I_ERI_F3x_S_F3y_S_M3_vrr+3*oned2e*I_ERI_F3x_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_F3x_S_D2y_S_M3_vrr;
      Double I_ERI_F2xy_S_G4y_S_M2_vrr = QCY*I_ERI_F2xy_S_F3y_S_M2_vrr+WQY*I_ERI_F2xy_S_F3y_S_M3_vrr+3*oned2e*I_ERI_F2xy_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_F2xy_S_D2y_S_M3_vrr+oned2k*I_ERI_D2x_S_F3y_S_M3_vrr;
      Double I_ERI_F2xz_S_G4y_S_M2_vrr = QCY*I_ERI_F2xz_S_F3y_S_M2_vrr+WQY*I_ERI_F2xz_S_F3y_S_M3_vrr+3*oned2e*I_ERI_F2xz_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_F2xz_S_D2y_S_M3_vrr;
      Double I_ERI_F3y_S_G4y_S_M2_vrr = QCY*I_ERI_F3y_S_F3y_S_M2_vrr+WQY*I_ERI_F3y_S_F3y_S_M3_vrr+3*oned2e*I_ERI_F3y_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_F3y_S_D2y_S_M3_vrr+3*oned2k*I_ERI_D2y_S_F3y_S_M3_vrr;
      Double I_ERI_F2yz_S_G4y_S_M2_vrr = QCY*I_ERI_F2yz_S_F3y_S_M2_vrr+WQY*I_ERI_F2yz_S_F3y_S_M3_vrr+3*oned2e*I_ERI_F2yz_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_F2yz_S_D2y_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_F3y_S_M3_vrr;
      Double I_ERI_F3z_S_G4y_S_M2_vrr = QCY*I_ERI_F3z_S_F3y_S_M2_vrr+WQY*I_ERI_F3z_S_F3y_S_M3_vrr+3*oned2e*I_ERI_F3z_S_D2y_S_M2_vrr-3*rhod2esq*I_ERI_F3z_S_D2y_S_M3_vrr;
      Double I_ERI_F3x_S_G3yz_S_M2_vrr = QCZ*I_ERI_F3x_S_F3y_S_M2_vrr+WQZ*I_ERI_F3x_S_F3y_S_M3_vrr;
      Double I_ERI_F2xy_S_G3yz_S_M2_vrr = QCZ*I_ERI_F2xy_S_F3y_S_M2_vrr+WQZ*I_ERI_F2xy_S_F3y_S_M3_vrr;
      Double I_ERI_F2xz_S_G3yz_S_M2_vrr = QCZ*I_ERI_F2xz_S_F3y_S_M2_vrr+WQZ*I_ERI_F2xz_S_F3y_S_M3_vrr+oned2k*I_ERI_D2x_S_F3y_S_M3_vrr;
      Double I_ERI_F3y_S_G3yz_S_M2_vrr = QCZ*I_ERI_F3y_S_F3y_S_M2_vrr+WQZ*I_ERI_F3y_S_F3y_S_M3_vrr;
      Double I_ERI_F2yz_S_G3yz_S_M2_vrr = QCZ*I_ERI_F2yz_S_F3y_S_M2_vrr+WQZ*I_ERI_F2yz_S_F3y_S_M3_vrr+oned2k*I_ERI_D2y_S_F3y_S_M3_vrr;
      Double I_ERI_F3z_S_G3yz_S_M2_vrr = QCZ*I_ERI_F3z_S_F3y_S_M2_vrr+WQZ*I_ERI_F3z_S_F3y_S_M3_vrr+3*oned2k*I_ERI_D2z_S_F3y_S_M3_vrr;
      Double I_ERI_F3x_S_G2y2z_S_M2_vrr = QCZ*I_ERI_F3x_S_F2yz_S_M2_vrr+WQZ*I_ERI_F3x_S_F2yz_S_M3_vrr+oned2e*I_ERI_F3x_S_D2y_S_M2_vrr-rhod2esq*I_ERI_F3x_S_D2y_S_M3_vrr;
      Double I_ERI_F2xy_S_G2y2z_S_M2_vrr = QCZ*I_ERI_F2xy_S_F2yz_S_M2_vrr+WQZ*I_ERI_F2xy_S_F2yz_S_M3_vrr+oned2e*I_ERI_F2xy_S_D2y_S_M2_vrr-rhod2esq*I_ERI_F2xy_S_D2y_S_M3_vrr;
      Double I_ERI_F2xz_S_G2y2z_S_M2_vrr = QCZ*I_ERI_F2xz_S_F2yz_S_M2_vrr+WQZ*I_ERI_F2xz_S_F2yz_S_M3_vrr+oned2e*I_ERI_F2xz_S_D2y_S_M2_vrr-rhod2esq*I_ERI_F2xz_S_D2y_S_M3_vrr+oned2k*I_ERI_D2x_S_F2yz_S_M3_vrr;
      Double I_ERI_F3y_S_G2y2z_S_M2_vrr = QCZ*I_ERI_F3y_S_F2yz_S_M2_vrr+WQZ*I_ERI_F3y_S_F2yz_S_M3_vrr+oned2e*I_ERI_F3y_S_D2y_S_M2_vrr-rhod2esq*I_ERI_F3y_S_D2y_S_M3_vrr;
      Double I_ERI_F2yz_S_G2y2z_S_M2_vrr = QCZ*I_ERI_F2yz_S_F2yz_S_M2_vrr+WQZ*I_ERI_F2yz_S_F2yz_S_M3_vrr+oned2e*I_ERI_F2yz_S_D2y_S_M2_vrr-rhod2esq*I_ERI_F2yz_S_D2y_S_M3_vrr+oned2k*I_ERI_D2y_S_F2yz_S_M3_vrr;
      Double I_ERI_F3z_S_G2y2z_S_M2_vrr = QCZ*I_ERI_F3z_S_F2yz_S_M2_vrr+WQZ*I_ERI_F3z_S_F2yz_S_M3_vrr+oned2e*I_ERI_F3z_S_D2y_S_M2_vrr-rhod2esq*I_ERI_F3z_S_D2y_S_M3_vrr+3*oned2k*I_ERI_D2z_S_F2yz_S_M3_vrr;
      Double I_ERI_F3x_S_Gy3z_S_M2_vrr = QCY*I_ERI_F3x_S_F3z_S_M2_vrr+WQY*I_ERI_F3x_S_F3z_S_M3_vrr;
      Double I_ERI_F2xy_S_Gy3z_S_M2_vrr = QCY*I_ERI_F2xy_S_F3z_S_M2_vrr+WQY*I_ERI_F2xy_S_F3z_S_M3_vrr+oned2k*I_ERI_D2x_S_F3z_S_M3_vrr;
      Double I_ERI_F2xz_S_Gy3z_S_M2_vrr = QCY*I_ERI_F2xz_S_F3z_S_M2_vrr+WQY*I_ERI_F2xz_S_F3z_S_M3_vrr;
      Double I_ERI_F3y_S_Gy3z_S_M2_vrr = QCY*I_ERI_F3y_S_F3z_S_M2_vrr+WQY*I_ERI_F3y_S_F3z_S_M3_vrr+3*oned2k*I_ERI_D2y_S_F3z_S_M3_vrr;
      Double I_ERI_F2yz_S_Gy3z_S_M2_vrr = QCY*I_ERI_F2yz_S_F3z_S_M2_vrr+WQY*I_ERI_F2yz_S_F3z_S_M3_vrr+2*oned2k*I_ERI_Dyz_S_F3z_S_M3_vrr;
      Double I_ERI_F3z_S_Gy3z_S_M2_vrr = QCY*I_ERI_F3z_S_F3z_S_M2_vrr+WQY*I_ERI_F3z_S_F3z_S_M3_vrr;
      Double I_ERI_F3x_S_G4z_S_M2_vrr = QCZ*I_ERI_F3x_S_F3z_S_M2_vrr+WQZ*I_ERI_F3x_S_F3z_S_M3_vrr+3*oned2e*I_ERI_F3x_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_F3x_S_D2z_S_M3_vrr;
      Double I_ERI_F2xy_S_G4z_S_M2_vrr = QCZ*I_ERI_F2xy_S_F3z_S_M2_vrr+WQZ*I_ERI_F2xy_S_F3z_S_M3_vrr+3*oned2e*I_ERI_F2xy_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_F2xy_S_D2z_S_M3_vrr;
      Double I_ERI_F2xz_S_G4z_S_M2_vrr = QCZ*I_ERI_F2xz_S_F3z_S_M2_vrr+WQZ*I_ERI_F2xz_S_F3z_S_M3_vrr+3*oned2e*I_ERI_F2xz_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_F2xz_S_D2z_S_M3_vrr+oned2k*I_ERI_D2x_S_F3z_S_M3_vrr;
      Double I_ERI_F3y_S_G4z_S_M2_vrr = QCZ*I_ERI_F3y_S_F3z_S_M2_vrr+WQZ*I_ERI_F3y_S_F3z_S_M3_vrr+3*oned2e*I_ERI_F3y_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_F3y_S_D2z_S_M3_vrr;
      Double I_ERI_F2yz_S_G4z_S_M2_vrr = QCZ*I_ERI_F2yz_S_F3z_S_M2_vrr+WQZ*I_ERI_F2yz_S_F3z_S_M3_vrr+3*oned2e*I_ERI_F2yz_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_F2yz_S_D2z_S_M3_vrr+oned2k*I_ERI_D2y_S_F3z_S_M3_vrr;
      Double I_ERI_F3z_S_G4z_S_M2_vrr = QCZ*I_ERI_F3z_S_F3z_S_M2_vrr+WQZ*I_ERI_F3z_S_F3z_S_M3_vrr+3*oned2e*I_ERI_F3z_S_D2z_S_M2_vrr-3*rhod2esq*I_ERI_F3z_S_D2z_S_M3_vrr+3*oned2k*I_ERI_D2z_S_F3z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 50 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M3
       ************************************************************/
      Double I_ERI_G4x_S_F3x_S_M2_vrr = PAX*I_ERI_F3x_S_F3x_S_M2_vrr+WPX*I_ERI_F3x_S_F3x_S_M3_vrr+3*oned2z*I_ERI_D2x_S_F3x_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_F3x_S_M3_vrr+3*oned2k*I_ERI_F3x_S_D2x_S_M3_vrr;
      Double I_ERI_G3xy_S_F3x_S_M2_vrr = PAY*I_ERI_F3x_S_F3x_S_M2_vrr+WPY*I_ERI_F3x_S_F3x_S_M3_vrr;
      Double I_ERI_G3xz_S_F3x_S_M2_vrr = PAZ*I_ERI_F3x_S_F3x_S_M2_vrr+WPZ*I_ERI_F3x_S_F3x_S_M3_vrr;
      Double I_ERI_G2x2y_S_F3x_S_M2_vrr = PAY*I_ERI_F2xy_S_F3x_S_M2_vrr+WPY*I_ERI_F2xy_S_F3x_S_M3_vrr+oned2z*I_ERI_D2x_S_F3x_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_F3x_S_M3_vrr;
      Double I_ERI_Gx3y_S_F3x_S_M2_vrr = PAX*I_ERI_F3y_S_F3x_S_M2_vrr+WPX*I_ERI_F3y_S_F3x_S_M3_vrr+3*oned2k*I_ERI_F3y_S_D2x_S_M3_vrr;
      Double I_ERI_Gx3z_S_F3x_S_M2_vrr = PAX*I_ERI_F3z_S_F3x_S_M2_vrr+WPX*I_ERI_F3z_S_F3x_S_M3_vrr+3*oned2k*I_ERI_F3z_S_D2x_S_M3_vrr;
      Double I_ERI_G4y_S_F3x_S_M2_vrr = PAY*I_ERI_F3y_S_F3x_S_M2_vrr+WPY*I_ERI_F3y_S_F3x_S_M3_vrr+3*oned2z*I_ERI_D2y_S_F3x_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_F3x_S_M3_vrr;
      Double I_ERI_G3yz_S_F3x_S_M2_vrr = PAZ*I_ERI_F3y_S_F3x_S_M2_vrr+WPZ*I_ERI_F3y_S_F3x_S_M3_vrr;
      Double I_ERI_Gy3z_S_F3x_S_M2_vrr = PAY*I_ERI_F3z_S_F3x_S_M2_vrr+WPY*I_ERI_F3z_S_F3x_S_M3_vrr;
      Double I_ERI_G4z_S_F3x_S_M2_vrr = PAZ*I_ERI_F3z_S_F3x_S_M2_vrr+WPZ*I_ERI_F3z_S_F3x_S_M3_vrr+3*oned2z*I_ERI_D2z_S_F3x_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_F3x_S_M3_vrr;
      Double I_ERI_G4x_S_F2xy_S_M2_vrr = PAX*I_ERI_F3x_S_F2xy_S_M2_vrr+WPX*I_ERI_F3x_S_F2xy_S_M3_vrr+3*oned2z*I_ERI_D2x_S_F2xy_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Dxy_S_M3_vrr;
      Double I_ERI_G3xy_S_F2xy_S_M2_vrr = PAY*I_ERI_F3x_S_F2xy_S_M2_vrr+WPY*I_ERI_F3x_S_F2xy_S_M3_vrr+oned2k*I_ERI_F3x_S_D2x_S_M3_vrr;
      Double I_ERI_G3xz_S_F2xy_S_M2_vrr = PAZ*I_ERI_F3x_S_F2xy_S_M2_vrr+WPZ*I_ERI_F3x_S_F2xy_S_M3_vrr;
      Double I_ERI_G2x2y_S_F2xy_S_M2_vrr = PAY*I_ERI_F2xy_S_F2xy_S_M2_vrr+WPY*I_ERI_F2xy_S_F2xy_S_M3_vrr+oned2z*I_ERI_D2x_S_F2xy_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_F2xy_S_M3_vrr+oned2k*I_ERI_F2xy_S_D2x_S_M3_vrr;
      Double I_ERI_Gx3y_S_F2xy_S_M2_vrr = PAX*I_ERI_F3y_S_F2xy_S_M2_vrr+WPX*I_ERI_F3y_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Dxy_S_M3_vrr;
      Double I_ERI_Gx3z_S_F2xy_S_M2_vrr = PAX*I_ERI_F3z_S_F2xy_S_M2_vrr+WPX*I_ERI_F3z_S_F2xy_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Dxy_S_M3_vrr;
      Double I_ERI_G4y_S_F2xy_S_M2_vrr = PAY*I_ERI_F3y_S_F2xy_S_M2_vrr+WPY*I_ERI_F3y_S_F2xy_S_M3_vrr+3*oned2z*I_ERI_D2y_S_F2xy_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_F2xy_S_M3_vrr+oned2k*I_ERI_F3y_S_D2x_S_M3_vrr;
      Double I_ERI_G3yz_S_F2xy_S_M2_vrr = PAZ*I_ERI_F3y_S_F2xy_S_M2_vrr+WPZ*I_ERI_F3y_S_F2xy_S_M3_vrr;
      Double I_ERI_Gy3z_S_F2xy_S_M2_vrr = PAY*I_ERI_F3z_S_F2xy_S_M2_vrr+WPY*I_ERI_F3z_S_F2xy_S_M3_vrr+oned2k*I_ERI_F3z_S_D2x_S_M3_vrr;
      Double I_ERI_G4z_S_F2xy_S_M2_vrr = PAZ*I_ERI_F3z_S_F2xy_S_M2_vrr+WPZ*I_ERI_F3z_S_F2xy_S_M3_vrr+3*oned2z*I_ERI_D2z_S_F2xy_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_F2xy_S_M3_vrr;
      Double I_ERI_G4x_S_F2xz_S_M2_vrr = PAX*I_ERI_F3x_S_F2xz_S_M2_vrr+WPX*I_ERI_F3x_S_F2xz_S_M3_vrr+3*oned2z*I_ERI_D2x_S_F2xz_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Dxz_S_M3_vrr;
      Double I_ERI_G3xy_S_F2xz_S_M2_vrr = PAY*I_ERI_F3x_S_F2xz_S_M2_vrr+WPY*I_ERI_F3x_S_F2xz_S_M3_vrr;
      Double I_ERI_G3xz_S_F2xz_S_M2_vrr = PAZ*I_ERI_F3x_S_F2xz_S_M2_vrr+WPZ*I_ERI_F3x_S_F2xz_S_M3_vrr+oned2k*I_ERI_F3x_S_D2x_S_M3_vrr;
      Double I_ERI_G2x2y_S_F2xz_S_M2_vrr = PAY*I_ERI_F2xy_S_F2xz_S_M2_vrr+WPY*I_ERI_F2xy_S_F2xz_S_M3_vrr+oned2z*I_ERI_D2x_S_F2xz_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_F2xz_S_M3_vrr;
      Double I_ERI_Gx3y_S_F2xz_S_M2_vrr = PAX*I_ERI_F3y_S_F2xz_S_M2_vrr+WPX*I_ERI_F3y_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Dxz_S_M3_vrr;
      Double I_ERI_Gx3z_S_F2xz_S_M2_vrr = PAX*I_ERI_F3z_S_F2xz_S_M2_vrr+WPX*I_ERI_F3z_S_F2xz_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Dxz_S_M3_vrr;
      Double I_ERI_G4y_S_F2xz_S_M2_vrr = PAY*I_ERI_F3y_S_F2xz_S_M2_vrr+WPY*I_ERI_F3y_S_F2xz_S_M3_vrr+3*oned2z*I_ERI_D2y_S_F2xz_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_F2xz_S_M3_vrr;
      Double I_ERI_G3yz_S_F2xz_S_M2_vrr = PAZ*I_ERI_F3y_S_F2xz_S_M2_vrr+WPZ*I_ERI_F3y_S_F2xz_S_M3_vrr+oned2k*I_ERI_F3y_S_D2x_S_M3_vrr;
      Double I_ERI_Gy3z_S_F2xz_S_M2_vrr = PAY*I_ERI_F3z_S_F2xz_S_M2_vrr+WPY*I_ERI_F3z_S_F2xz_S_M3_vrr;
      Double I_ERI_G4z_S_F2xz_S_M2_vrr = PAZ*I_ERI_F3z_S_F2xz_S_M2_vrr+WPZ*I_ERI_F3z_S_F2xz_S_M3_vrr+3*oned2z*I_ERI_D2z_S_F2xz_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_F2xz_S_M3_vrr+oned2k*I_ERI_F3z_S_D2x_S_M3_vrr;
      Double I_ERI_G4x_S_Fx2y_S_M2_vrr = PAX*I_ERI_F3x_S_Fx2y_S_M2_vrr+WPX*I_ERI_F3x_S_Fx2y_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Fx2y_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Fx2y_S_M3_vrr+oned2k*I_ERI_F3x_S_D2y_S_M3_vrr;
      Double I_ERI_G3xy_S_Fx2y_S_M2_vrr = PAY*I_ERI_F3x_S_Fx2y_S_M2_vrr+WPY*I_ERI_F3x_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Dxy_S_M3_vrr;
      Double I_ERI_G3xz_S_Fx2y_S_M2_vrr = PAZ*I_ERI_F3x_S_Fx2y_S_M2_vrr+WPZ*I_ERI_F3x_S_Fx2y_S_M3_vrr;
      Double I_ERI_G2x2y_S_Fx2y_S_M2_vrr = PAY*I_ERI_F2xy_S_Fx2y_S_M2_vrr+WPY*I_ERI_F2xy_S_Fx2y_S_M3_vrr+oned2z*I_ERI_D2x_S_Fx2y_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_F2xy_S_Dxy_S_M3_vrr;
      Double I_ERI_Gx3y_S_Fx2y_S_M2_vrr = PAX*I_ERI_F3y_S_Fx2y_S_M2_vrr+WPX*I_ERI_F3y_S_Fx2y_S_M3_vrr+oned2k*I_ERI_F3y_S_D2y_S_M3_vrr;
      Double I_ERI_Gx3z_S_Fx2y_S_M2_vrr = PAX*I_ERI_F3z_S_Fx2y_S_M2_vrr+WPX*I_ERI_F3z_S_Fx2y_S_M3_vrr+oned2k*I_ERI_F3z_S_D2y_S_M3_vrr;
      Double I_ERI_G4y_S_Fx2y_S_M2_vrr = PAY*I_ERI_F3y_S_Fx2y_S_M2_vrr+WPY*I_ERI_F3y_S_Fx2y_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Fx2y_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Dxy_S_M3_vrr;
      Double I_ERI_G3yz_S_Fx2y_S_M2_vrr = PAZ*I_ERI_F3y_S_Fx2y_S_M2_vrr+WPZ*I_ERI_F3y_S_Fx2y_S_M3_vrr;
      Double I_ERI_Gy3z_S_Fx2y_S_M2_vrr = PAY*I_ERI_F3z_S_Fx2y_S_M2_vrr+WPY*I_ERI_F3z_S_Fx2y_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Dxy_S_M3_vrr;
      Double I_ERI_G4z_S_Fx2y_S_M2_vrr = PAZ*I_ERI_F3z_S_Fx2y_S_M2_vrr+WPZ*I_ERI_F3z_S_Fx2y_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Fx2y_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Fx2y_S_M3_vrr;
      Double I_ERI_G4x_S_Fxyz_S_M2_vrr = PAX*I_ERI_F3x_S_Fxyz_S_M2_vrr+WPX*I_ERI_F3x_S_Fxyz_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Fxyz_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F3x_S_Dyz_S_M3_vrr;
      Double I_ERI_G3xy_S_Fxyz_S_M2_vrr = PAY*I_ERI_F3x_S_Fxyz_S_M2_vrr+WPY*I_ERI_F3x_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F3x_S_Dxz_S_M3_vrr;
      Double I_ERI_G3xz_S_Fxyz_S_M2_vrr = PAZ*I_ERI_F3x_S_Fxyz_S_M2_vrr+WPZ*I_ERI_F3x_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F3x_S_Dxy_S_M3_vrr;
      Double I_ERI_G2x2y_S_Fxyz_S_M2_vrr = PAY*I_ERI_F2xy_S_Fxyz_S_M2_vrr+WPY*I_ERI_F2xy_S_Fxyz_S_M3_vrr+oned2z*I_ERI_D2x_S_Fxyz_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F2xy_S_Dxz_S_M3_vrr;
      Double I_ERI_Gx3y_S_Fxyz_S_M2_vrr = PAX*I_ERI_F3y_S_Fxyz_S_M2_vrr+WPX*I_ERI_F3y_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F3y_S_Dyz_S_M3_vrr;
      Double I_ERI_Gx3z_S_Fxyz_S_M2_vrr = PAX*I_ERI_F3z_S_Fxyz_S_M2_vrr+WPX*I_ERI_F3z_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F3z_S_Dyz_S_M3_vrr;
      Double I_ERI_G4y_S_Fxyz_S_M2_vrr = PAY*I_ERI_F3y_S_Fxyz_S_M2_vrr+WPY*I_ERI_F3y_S_Fxyz_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Fxyz_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F3y_S_Dxz_S_M3_vrr;
      Double I_ERI_G3yz_S_Fxyz_S_M2_vrr = PAZ*I_ERI_F3y_S_Fxyz_S_M2_vrr+WPZ*I_ERI_F3y_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F3y_S_Dxy_S_M3_vrr;
      Double I_ERI_Gy3z_S_Fxyz_S_M2_vrr = PAY*I_ERI_F3z_S_Fxyz_S_M2_vrr+WPY*I_ERI_F3z_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F3z_S_Dxz_S_M3_vrr;
      Double I_ERI_G4z_S_Fxyz_S_M2_vrr = PAZ*I_ERI_F3z_S_Fxyz_S_M2_vrr+WPZ*I_ERI_F3z_S_Fxyz_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Fxyz_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Fxyz_S_M3_vrr+oned2k*I_ERI_F3z_S_Dxy_S_M3_vrr;
      Double I_ERI_G4x_S_Fx2z_S_M2_vrr = PAX*I_ERI_F3x_S_Fx2z_S_M2_vrr+WPX*I_ERI_F3x_S_Fx2z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Fx2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Fx2z_S_M3_vrr+oned2k*I_ERI_F3x_S_D2z_S_M3_vrr;
      Double I_ERI_G3xy_S_Fx2z_S_M2_vrr = PAY*I_ERI_F3x_S_Fx2z_S_M2_vrr+WPY*I_ERI_F3x_S_Fx2z_S_M3_vrr;
      Double I_ERI_G3xz_S_Fx2z_S_M2_vrr = PAZ*I_ERI_F3x_S_Fx2z_S_M2_vrr+WPZ*I_ERI_F3x_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Dxz_S_M3_vrr;
      Double I_ERI_G2x2y_S_Fx2z_S_M2_vrr = PAY*I_ERI_F2xy_S_Fx2z_S_M2_vrr+WPY*I_ERI_F2xy_S_Fx2z_S_M3_vrr+oned2z*I_ERI_D2x_S_Fx2z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Fx2z_S_M3_vrr;
      Double I_ERI_Gx3y_S_Fx2z_S_M2_vrr = PAX*I_ERI_F3y_S_Fx2z_S_M2_vrr+WPX*I_ERI_F3y_S_Fx2z_S_M3_vrr+oned2k*I_ERI_F3y_S_D2z_S_M3_vrr;
      Double I_ERI_Gx3z_S_Fx2z_S_M2_vrr = PAX*I_ERI_F3z_S_Fx2z_S_M2_vrr+WPX*I_ERI_F3z_S_Fx2z_S_M3_vrr+oned2k*I_ERI_F3z_S_D2z_S_M3_vrr;
      Double I_ERI_G4y_S_Fx2z_S_M2_vrr = PAY*I_ERI_F3y_S_Fx2z_S_M2_vrr+WPY*I_ERI_F3y_S_Fx2z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Fx2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Fx2z_S_M3_vrr;
      Double I_ERI_G3yz_S_Fx2z_S_M2_vrr = PAZ*I_ERI_F3y_S_Fx2z_S_M2_vrr+WPZ*I_ERI_F3y_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Dxz_S_M3_vrr;
      Double I_ERI_Gy3z_S_Fx2z_S_M2_vrr = PAY*I_ERI_F3z_S_Fx2z_S_M2_vrr+WPY*I_ERI_F3z_S_Fx2z_S_M3_vrr;
      Double I_ERI_G4z_S_Fx2z_S_M2_vrr = PAZ*I_ERI_F3z_S_Fx2z_S_M2_vrr+WPZ*I_ERI_F3z_S_Fx2z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Fx2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Fx2z_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Dxz_S_M3_vrr;
      Double I_ERI_G4x_S_F3y_S_M2_vrr = PAX*I_ERI_F3x_S_F3y_S_M2_vrr+WPX*I_ERI_F3x_S_F3y_S_M3_vrr+3*oned2z*I_ERI_D2x_S_F3y_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_F3y_S_M3_vrr;
      Double I_ERI_G3xy_S_F3y_S_M2_vrr = PAY*I_ERI_F3x_S_F3y_S_M2_vrr+WPY*I_ERI_F3x_S_F3y_S_M3_vrr+3*oned2k*I_ERI_F3x_S_D2y_S_M3_vrr;
      Double I_ERI_G3xz_S_F3y_S_M2_vrr = PAZ*I_ERI_F3x_S_F3y_S_M2_vrr+WPZ*I_ERI_F3x_S_F3y_S_M3_vrr;
      Double I_ERI_G2x2y_S_F3y_S_M2_vrr = PAY*I_ERI_F2xy_S_F3y_S_M2_vrr+WPY*I_ERI_F2xy_S_F3y_S_M3_vrr+oned2z*I_ERI_D2x_S_F3y_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_F3y_S_M3_vrr+3*oned2k*I_ERI_F2xy_S_D2y_S_M3_vrr;
      Double I_ERI_Gx3y_S_F3y_S_M2_vrr = PAX*I_ERI_F3y_S_F3y_S_M2_vrr+WPX*I_ERI_F3y_S_F3y_S_M3_vrr;
      Double I_ERI_Gx3z_S_F3y_S_M2_vrr = PAX*I_ERI_F3z_S_F3y_S_M2_vrr+WPX*I_ERI_F3z_S_F3y_S_M3_vrr;
      Double I_ERI_G4y_S_F3y_S_M2_vrr = PAY*I_ERI_F3y_S_F3y_S_M2_vrr+WPY*I_ERI_F3y_S_F3y_S_M3_vrr+3*oned2z*I_ERI_D2y_S_F3y_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_F3y_S_M3_vrr+3*oned2k*I_ERI_F3y_S_D2y_S_M3_vrr;
      Double I_ERI_G3yz_S_F3y_S_M2_vrr = PAZ*I_ERI_F3y_S_F3y_S_M2_vrr+WPZ*I_ERI_F3y_S_F3y_S_M3_vrr;
      Double I_ERI_Gy3z_S_F3y_S_M2_vrr = PAY*I_ERI_F3z_S_F3y_S_M2_vrr+WPY*I_ERI_F3z_S_F3y_S_M3_vrr+3*oned2k*I_ERI_F3z_S_D2y_S_M3_vrr;
      Double I_ERI_G4z_S_F3y_S_M2_vrr = PAZ*I_ERI_F3z_S_F3y_S_M2_vrr+WPZ*I_ERI_F3z_S_F3y_S_M3_vrr+3*oned2z*I_ERI_D2z_S_F3y_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_F3y_S_M3_vrr;
      Double I_ERI_G4x_S_F2yz_S_M2_vrr = PAX*I_ERI_F3x_S_F2yz_S_M2_vrr+WPX*I_ERI_F3x_S_F2yz_S_M3_vrr+3*oned2z*I_ERI_D2x_S_F2yz_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_F2yz_S_M3_vrr;
      Double I_ERI_G3xy_S_F2yz_S_M2_vrr = PAY*I_ERI_F3x_S_F2yz_S_M2_vrr+WPY*I_ERI_F3x_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Dyz_S_M3_vrr;
      Double I_ERI_G3xz_S_F2yz_S_M2_vrr = PAZ*I_ERI_F3x_S_F2yz_S_M2_vrr+WPZ*I_ERI_F3x_S_F2yz_S_M3_vrr+oned2k*I_ERI_F3x_S_D2y_S_M3_vrr;
      Double I_ERI_G2x2y_S_F2yz_S_M2_vrr = PAY*I_ERI_F2xy_S_F2yz_S_M2_vrr+WPY*I_ERI_F2xy_S_F2yz_S_M3_vrr+oned2z*I_ERI_D2x_S_F2yz_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_F2xy_S_Dyz_S_M3_vrr;
      Double I_ERI_Gx3y_S_F2yz_S_M2_vrr = PAX*I_ERI_F3y_S_F2yz_S_M2_vrr+WPX*I_ERI_F3y_S_F2yz_S_M3_vrr;
      Double I_ERI_Gx3z_S_F2yz_S_M2_vrr = PAX*I_ERI_F3z_S_F2yz_S_M2_vrr+WPX*I_ERI_F3z_S_F2yz_S_M3_vrr;
      Double I_ERI_G4y_S_F2yz_S_M2_vrr = PAY*I_ERI_F3y_S_F2yz_S_M2_vrr+WPY*I_ERI_F3y_S_F2yz_S_M3_vrr+3*oned2z*I_ERI_D2y_S_F2yz_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Dyz_S_M3_vrr;
      Double I_ERI_G3yz_S_F2yz_S_M2_vrr = PAZ*I_ERI_F3y_S_F2yz_S_M2_vrr+WPZ*I_ERI_F3y_S_F2yz_S_M3_vrr+oned2k*I_ERI_F3y_S_D2y_S_M3_vrr;
      Double I_ERI_Gy3z_S_F2yz_S_M2_vrr = PAY*I_ERI_F3z_S_F2yz_S_M2_vrr+WPY*I_ERI_F3z_S_F2yz_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Dyz_S_M3_vrr;
      Double I_ERI_G4z_S_F2yz_S_M2_vrr = PAZ*I_ERI_F3z_S_F2yz_S_M2_vrr+WPZ*I_ERI_F3z_S_F2yz_S_M3_vrr+3*oned2z*I_ERI_D2z_S_F2yz_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_F2yz_S_M3_vrr+oned2k*I_ERI_F3z_S_D2y_S_M3_vrr;
      Double I_ERI_G4x_S_Fy2z_S_M2_vrr = PAX*I_ERI_F3x_S_Fy2z_S_M2_vrr+WPX*I_ERI_F3x_S_Fy2z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Fy2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Fy2z_S_M3_vrr;
      Double I_ERI_G3xy_S_Fy2z_S_M2_vrr = PAY*I_ERI_F3x_S_Fy2z_S_M2_vrr+WPY*I_ERI_F3x_S_Fy2z_S_M3_vrr+oned2k*I_ERI_F3x_S_D2z_S_M3_vrr;
      Double I_ERI_G3xz_S_Fy2z_S_M2_vrr = PAZ*I_ERI_F3x_S_Fy2z_S_M2_vrr+WPZ*I_ERI_F3x_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Dyz_S_M3_vrr;
      Double I_ERI_G2x2y_S_Fy2z_S_M2_vrr = PAY*I_ERI_F2xy_S_Fy2z_S_M2_vrr+WPY*I_ERI_F2xy_S_Fy2z_S_M3_vrr+oned2z*I_ERI_D2x_S_Fy2z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Fy2z_S_M3_vrr+oned2k*I_ERI_F2xy_S_D2z_S_M3_vrr;
      Double I_ERI_Gx3y_S_Fy2z_S_M2_vrr = PAX*I_ERI_F3y_S_Fy2z_S_M2_vrr+WPX*I_ERI_F3y_S_Fy2z_S_M3_vrr;
      Double I_ERI_Gx3z_S_Fy2z_S_M2_vrr = PAX*I_ERI_F3z_S_Fy2z_S_M2_vrr+WPX*I_ERI_F3z_S_Fy2z_S_M3_vrr;
      Double I_ERI_G4y_S_Fy2z_S_M2_vrr = PAY*I_ERI_F3y_S_Fy2z_S_M2_vrr+WPY*I_ERI_F3y_S_Fy2z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Fy2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Fy2z_S_M3_vrr+oned2k*I_ERI_F3y_S_D2z_S_M3_vrr;
      Double I_ERI_G3yz_S_Fy2z_S_M2_vrr = PAZ*I_ERI_F3y_S_Fy2z_S_M2_vrr+WPZ*I_ERI_F3y_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Dyz_S_M3_vrr;
      Double I_ERI_Gy3z_S_Fy2z_S_M2_vrr = PAY*I_ERI_F3z_S_Fy2z_S_M2_vrr+WPY*I_ERI_F3z_S_Fy2z_S_M3_vrr+oned2k*I_ERI_F3z_S_D2z_S_M3_vrr;
      Double I_ERI_G4z_S_Fy2z_S_M2_vrr = PAZ*I_ERI_F3z_S_Fy2z_S_M2_vrr+WPZ*I_ERI_F3z_S_Fy2z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Fy2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Fy2z_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Dyz_S_M3_vrr;
      Double I_ERI_G4x_S_F3z_S_M2_vrr = PAX*I_ERI_F3x_S_F3z_S_M2_vrr+WPX*I_ERI_F3x_S_F3z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_F3z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_F3z_S_M3_vrr;
      Double I_ERI_G3xy_S_F3z_S_M2_vrr = PAY*I_ERI_F3x_S_F3z_S_M2_vrr+WPY*I_ERI_F3x_S_F3z_S_M3_vrr;
      Double I_ERI_G3xz_S_F3z_S_M2_vrr = PAZ*I_ERI_F3x_S_F3z_S_M2_vrr+WPZ*I_ERI_F3x_S_F3z_S_M3_vrr+3*oned2k*I_ERI_F3x_S_D2z_S_M3_vrr;
      Double I_ERI_G2x2y_S_F3z_S_M2_vrr = PAY*I_ERI_F2xy_S_F3z_S_M2_vrr+WPY*I_ERI_F2xy_S_F3z_S_M3_vrr+oned2z*I_ERI_D2x_S_F3z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_F3z_S_M3_vrr;
      Double I_ERI_Gx3y_S_F3z_S_M2_vrr = PAX*I_ERI_F3y_S_F3z_S_M2_vrr+WPX*I_ERI_F3y_S_F3z_S_M3_vrr;
      Double I_ERI_Gx3z_S_F3z_S_M2_vrr = PAX*I_ERI_F3z_S_F3z_S_M2_vrr+WPX*I_ERI_F3z_S_F3z_S_M3_vrr;
      Double I_ERI_G4y_S_F3z_S_M2_vrr = PAY*I_ERI_F3y_S_F3z_S_M2_vrr+WPY*I_ERI_F3y_S_F3z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_F3z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_F3z_S_M3_vrr;
      Double I_ERI_G3yz_S_F3z_S_M2_vrr = PAZ*I_ERI_F3y_S_F3z_S_M2_vrr+WPZ*I_ERI_F3y_S_F3z_S_M3_vrr+3*oned2k*I_ERI_F3y_S_D2z_S_M3_vrr;
      Double I_ERI_Gy3z_S_F3z_S_M2_vrr = PAY*I_ERI_F3z_S_F3z_S_M2_vrr+WPY*I_ERI_F3z_S_F3z_S_M3_vrr;
      Double I_ERI_G4z_S_F3z_S_M2_vrr = PAZ*I_ERI_F3z_S_F3z_S_M2_vrr+WPZ*I_ERI_F3z_S_F3z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_F3z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_F3z_S_M3_vrr+3*oned2k*I_ERI_F3z_S_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_G_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 75 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M3
       ************************************************************/
      Double I_ERI_G4x_S_G4x_S_M2_vrr = PAX*I_ERI_F3x_S_G4x_S_M2_vrr+WPX*I_ERI_F3x_S_G4x_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G4x_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G4x_S_M3_vrr+4*oned2k*I_ERI_F3x_S_F3x_S_M3_vrr;
      Double I_ERI_G3xy_S_G4x_S_M2_vrr = PAY*I_ERI_F3x_S_G4x_S_M2_vrr+WPY*I_ERI_F3x_S_G4x_S_M3_vrr;
      Double I_ERI_G3xz_S_G4x_S_M2_vrr = PAZ*I_ERI_F3x_S_G4x_S_M2_vrr+WPZ*I_ERI_F3x_S_G4x_S_M3_vrr;
      Double I_ERI_G2x2y_S_G4x_S_M2_vrr = PAY*I_ERI_F2xy_S_G4x_S_M2_vrr+WPY*I_ERI_F2xy_S_G4x_S_M3_vrr+oned2z*I_ERI_D2x_S_G4x_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G4x_S_M3_vrr;
      Double I_ERI_Gx3y_S_G4x_S_M2_vrr = PAX*I_ERI_F3y_S_G4x_S_M2_vrr+WPX*I_ERI_F3y_S_G4x_S_M3_vrr+4*oned2k*I_ERI_F3y_S_F3x_S_M3_vrr;
      Double I_ERI_Gx3z_S_G4x_S_M2_vrr = PAX*I_ERI_F3z_S_G4x_S_M2_vrr+WPX*I_ERI_F3z_S_G4x_S_M3_vrr+4*oned2k*I_ERI_F3z_S_F3x_S_M3_vrr;
      Double I_ERI_G4y_S_G4x_S_M2_vrr = PAY*I_ERI_F3y_S_G4x_S_M2_vrr+WPY*I_ERI_F3y_S_G4x_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G4x_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G4x_S_M3_vrr;
      Double I_ERI_G3yz_S_G4x_S_M2_vrr = PAZ*I_ERI_F3y_S_G4x_S_M2_vrr+WPZ*I_ERI_F3y_S_G4x_S_M3_vrr;
      Double I_ERI_Gy3z_S_G4x_S_M2_vrr = PAY*I_ERI_F3z_S_G4x_S_M2_vrr+WPY*I_ERI_F3z_S_G4x_S_M3_vrr;
      Double I_ERI_G4z_S_G4x_S_M2_vrr = PAZ*I_ERI_F3z_S_G4x_S_M2_vrr+WPZ*I_ERI_F3z_S_G4x_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G4x_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G4x_S_M3_vrr;
      Double I_ERI_G4x_S_G3xy_S_M2_vrr = PAX*I_ERI_F3x_S_G3xy_S_M2_vrr+WPX*I_ERI_F3x_S_G3xy_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G3xy_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G3xy_S_M3_vrr+3*oned2k*I_ERI_F3x_S_F2xy_S_M3_vrr;
      Double I_ERI_G3xy_S_G3xy_S_M2_vrr = PAY*I_ERI_F3x_S_G3xy_S_M2_vrr+WPY*I_ERI_F3x_S_G3xy_S_M3_vrr+oned2k*I_ERI_F3x_S_F3x_S_M3_vrr;
      Double I_ERI_G3xz_S_G3xy_S_M2_vrr = PAZ*I_ERI_F3x_S_G3xy_S_M2_vrr+WPZ*I_ERI_F3x_S_G3xy_S_M3_vrr;
      Double I_ERI_G2x2y_S_G3xy_S_M2_vrr = PAY*I_ERI_F2xy_S_G3xy_S_M2_vrr+WPY*I_ERI_F2xy_S_G3xy_S_M3_vrr+oned2z*I_ERI_D2x_S_G3xy_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G3xy_S_M3_vrr+oned2k*I_ERI_F2xy_S_F3x_S_M3_vrr;
      Double I_ERI_Gx3y_S_G3xy_S_M2_vrr = PAX*I_ERI_F3y_S_G3xy_S_M2_vrr+WPX*I_ERI_F3y_S_G3xy_S_M3_vrr+3*oned2k*I_ERI_F3y_S_F2xy_S_M3_vrr;
      Double I_ERI_Gx3z_S_G3xy_S_M2_vrr = PAX*I_ERI_F3z_S_G3xy_S_M2_vrr+WPX*I_ERI_F3z_S_G3xy_S_M3_vrr+3*oned2k*I_ERI_F3z_S_F2xy_S_M3_vrr;
      Double I_ERI_G4y_S_G3xy_S_M2_vrr = PAY*I_ERI_F3y_S_G3xy_S_M2_vrr+WPY*I_ERI_F3y_S_G3xy_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G3xy_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G3xy_S_M3_vrr+oned2k*I_ERI_F3y_S_F3x_S_M3_vrr;
      Double I_ERI_G3yz_S_G3xy_S_M2_vrr = PAZ*I_ERI_F3y_S_G3xy_S_M2_vrr+WPZ*I_ERI_F3y_S_G3xy_S_M3_vrr;
      Double I_ERI_Gy3z_S_G3xy_S_M2_vrr = PAY*I_ERI_F3z_S_G3xy_S_M2_vrr+WPY*I_ERI_F3z_S_G3xy_S_M3_vrr+oned2k*I_ERI_F3z_S_F3x_S_M3_vrr;
      Double I_ERI_G4z_S_G3xy_S_M2_vrr = PAZ*I_ERI_F3z_S_G3xy_S_M2_vrr+WPZ*I_ERI_F3z_S_G3xy_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G3xy_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G3xy_S_M3_vrr;
      Double I_ERI_G4x_S_G3xz_S_M2_vrr = PAX*I_ERI_F3x_S_G3xz_S_M2_vrr+WPX*I_ERI_F3x_S_G3xz_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G3xz_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G3xz_S_M3_vrr+3*oned2k*I_ERI_F3x_S_F2xz_S_M3_vrr;
      Double I_ERI_G3xy_S_G3xz_S_M2_vrr = PAY*I_ERI_F3x_S_G3xz_S_M2_vrr+WPY*I_ERI_F3x_S_G3xz_S_M3_vrr;
      Double I_ERI_G3xz_S_G3xz_S_M2_vrr = PAZ*I_ERI_F3x_S_G3xz_S_M2_vrr+WPZ*I_ERI_F3x_S_G3xz_S_M3_vrr+oned2k*I_ERI_F3x_S_F3x_S_M3_vrr;
      Double I_ERI_G2x2y_S_G3xz_S_M2_vrr = PAY*I_ERI_F2xy_S_G3xz_S_M2_vrr+WPY*I_ERI_F2xy_S_G3xz_S_M3_vrr+oned2z*I_ERI_D2x_S_G3xz_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G3xz_S_M3_vrr;
      Double I_ERI_Gx3y_S_G3xz_S_M2_vrr = PAX*I_ERI_F3y_S_G3xz_S_M2_vrr+WPX*I_ERI_F3y_S_G3xz_S_M3_vrr+3*oned2k*I_ERI_F3y_S_F2xz_S_M3_vrr;
      Double I_ERI_Gx3z_S_G3xz_S_M2_vrr = PAX*I_ERI_F3z_S_G3xz_S_M2_vrr+WPX*I_ERI_F3z_S_G3xz_S_M3_vrr+3*oned2k*I_ERI_F3z_S_F2xz_S_M3_vrr;
      Double I_ERI_G4y_S_G3xz_S_M2_vrr = PAY*I_ERI_F3y_S_G3xz_S_M2_vrr+WPY*I_ERI_F3y_S_G3xz_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G3xz_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G3xz_S_M3_vrr;
      Double I_ERI_G3yz_S_G3xz_S_M2_vrr = PAZ*I_ERI_F3y_S_G3xz_S_M2_vrr+WPZ*I_ERI_F3y_S_G3xz_S_M3_vrr+oned2k*I_ERI_F3y_S_F3x_S_M3_vrr;
      Double I_ERI_Gy3z_S_G3xz_S_M2_vrr = PAY*I_ERI_F3z_S_G3xz_S_M2_vrr+WPY*I_ERI_F3z_S_G3xz_S_M3_vrr;
      Double I_ERI_G4z_S_G3xz_S_M2_vrr = PAZ*I_ERI_F3z_S_G3xz_S_M2_vrr+WPZ*I_ERI_F3z_S_G3xz_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G3xz_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G3xz_S_M3_vrr+oned2k*I_ERI_F3z_S_F3x_S_M3_vrr;
      Double I_ERI_G4x_S_G2x2y_S_M2_vrr = PAX*I_ERI_F3x_S_G2x2y_S_M2_vrr+WPX*I_ERI_F3x_S_G2x2y_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G2x2y_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Fx2y_S_M3_vrr;
      Double I_ERI_G3xy_S_G2x2y_S_M2_vrr = PAY*I_ERI_F3x_S_G2x2y_S_M2_vrr+WPY*I_ERI_F3x_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_F3x_S_F2xy_S_M3_vrr;
      Double I_ERI_G3xz_S_G2x2y_S_M2_vrr = PAZ*I_ERI_F3x_S_G2x2y_S_M2_vrr+WPZ*I_ERI_F3x_S_G2x2y_S_M3_vrr;
      Double I_ERI_G2x2y_S_G2x2y_S_M2_vrr = PAY*I_ERI_F2xy_S_G2x2y_S_M2_vrr+WPY*I_ERI_F2xy_S_G2x2y_S_M3_vrr+oned2z*I_ERI_D2x_S_G2x2y_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_F2xy_S_F2xy_S_M3_vrr;
      Double I_ERI_Gx3y_S_G2x2y_S_M2_vrr = PAX*I_ERI_F3y_S_G2x2y_S_M2_vrr+WPX*I_ERI_F3y_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Fx2y_S_M3_vrr;
      Double I_ERI_Gx3z_S_G2x2y_S_M2_vrr = PAX*I_ERI_F3z_S_G2x2y_S_M2_vrr+WPX*I_ERI_F3z_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Fx2y_S_M3_vrr;
      Double I_ERI_G4y_S_G2x2y_S_M2_vrr = PAY*I_ERI_F3y_S_G2x2y_S_M2_vrr+WPY*I_ERI_F3y_S_G2x2y_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G2x2y_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_F3y_S_F2xy_S_M3_vrr;
      Double I_ERI_G3yz_S_G2x2y_S_M2_vrr = PAZ*I_ERI_F3y_S_G2x2y_S_M2_vrr+WPZ*I_ERI_F3y_S_G2x2y_S_M3_vrr;
      Double I_ERI_Gy3z_S_G2x2y_S_M2_vrr = PAY*I_ERI_F3z_S_G2x2y_S_M2_vrr+WPY*I_ERI_F3z_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_F3z_S_F2xy_S_M3_vrr;
      Double I_ERI_G4z_S_G2x2y_S_M2_vrr = PAZ*I_ERI_F3z_S_G2x2y_S_M2_vrr+WPZ*I_ERI_F3z_S_G2x2y_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G2x2y_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G2x2y_S_M3_vrr;
      Double I_ERI_G4x_S_G2xyz_S_M2_vrr = PAX*I_ERI_F3x_S_G2xyz_S_M2_vrr+WPX*I_ERI_F3x_S_G2xyz_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G2xyz_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G2xyz_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M3_vrr;
      Double I_ERI_G3xy_S_G2xyz_S_M2_vrr = PAY*I_ERI_F3x_S_G2xyz_S_M2_vrr+WPY*I_ERI_F3x_S_G2xyz_S_M3_vrr+oned2k*I_ERI_F3x_S_F2xz_S_M3_vrr;
      Double I_ERI_G3xz_S_G2xyz_S_M2_vrr = PAZ*I_ERI_F3x_S_G2xyz_S_M2_vrr+WPZ*I_ERI_F3x_S_G2xyz_S_M3_vrr+oned2k*I_ERI_F3x_S_F2xy_S_M3_vrr;
      Double I_ERI_G2x2y_S_G2xyz_S_M2_vrr = PAY*I_ERI_F2xy_S_G2xyz_S_M2_vrr+WPY*I_ERI_F2xy_S_G2xyz_S_M3_vrr+oned2z*I_ERI_D2x_S_G2xyz_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G2xyz_S_M3_vrr+oned2k*I_ERI_F2xy_S_F2xz_S_M3_vrr;
      Double I_ERI_Gx3y_S_G2xyz_S_M2_vrr = PAX*I_ERI_F3y_S_G2xyz_S_M2_vrr+WPX*I_ERI_F3y_S_G2xyz_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M3_vrr;
      Double I_ERI_Gx3z_S_G2xyz_S_M2_vrr = PAX*I_ERI_F3z_S_G2xyz_S_M2_vrr+WPX*I_ERI_F3z_S_G2xyz_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M3_vrr;
      Double I_ERI_G4y_S_G2xyz_S_M2_vrr = PAY*I_ERI_F3y_S_G2xyz_S_M2_vrr+WPY*I_ERI_F3y_S_G2xyz_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G2xyz_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G2xyz_S_M3_vrr+oned2k*I_ERI_F3y_S_F2xz_S_M3_vrr;
      Double I_ERI_G3yz_S_G2xyz_S_M2_vrr = PAZ*I_ERI_F3y_S_G2xyz_S_M2_vrr+WPZ*I_ERI_F3y_S_G2xyz_S_M3_vrr+oned2k*I_ERI_F3y_S_F2xy_S_M3_vrr;
      Double I_ERI_Gy3z_S_G2xyz_S_M2_vrr = PAY*I_ERI_F3z_S_G2xyz_S_M2_vrr+WPY*I_ERI_F3z_S_G2xyz_S_M3_vrr+oned2k*I_ERI_F3z_S_F2xz_S_M3_vrr;
      Double I_ERI_G4z_S_G2xyz_S_M2_vrr = PAZ*I_ERI_F3z_S_G2xyz_S_M2_vrr+WPZ*I_ERI_F3z_S_G2xyz_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G2xyz_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G2xyz_S_M3_vrr+oned2k*I_ERI_F3z_S_F2xy_S_M3_vrr;
      Double I_ERI_G4x_S_G2x2z_S_M2_vrr = PAX*I_ERI_F3x_S_G2x2z_S_M2_vrr+WPX*I_ERI_F3x_S_G2x2z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G2x2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Fx2z_S_M3_vrr;
      Double I_ERI_G3xy_S_G2x2z_S_M2_vrr = PAY*I_ERI_F3x_S_G2x2z_S_M2_vrr+WPY*I_ERI_F3x_S_G2x2z_S_M3_vrr;
      Double I_ERI_G3xz_S_G2x2z_S_M2_vrr = PAZ*I_ERI_F3x_S_G2x2z_S_M2_vrr+WPZ*I_ERI_F3x_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_F3x_S_F2xz_S_M3_vrr;
      Double I_ERI_G2x2y_S_G2x2z_S_M2_vrr = PAY*I_ERI_F2xy_S_G2x2z_S_M2_vrr+WPY*I_ERI_F2xy_S_G2x2z_S_M3_vrr+oned2z*I_ERI_D2x_S_G2x2z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G2x2z_S_M3_vrr;
      Double I_ERI_Gx3y_S_G2x2z_S_M2_vrr = PAX*I_ERI_F3y_S_G2x2z_S_M2_vrr+WPX*I_ERI_F3y_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Fx2z_S_M3_vrr;
      Double I_ERI_Gx3z_S_G2x2z_S_M2_vrr = PAX*I_ERI_F3z_S_G2x2z_S_M2_vrr+WPX*I_ERI_F3z_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Fx2z_S_M3_vrr;
      Double I_ERI_G4y_S_G2x2z_S_M2_vrr = PAY*I_ERI_F3y_S_G2x2z_S_M2_vrr+WPY*I_ERI_F3y_S_G2x2z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G2x2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G2x2z_S_M3_vrr;
      Double I_ERI_G3yz_S_G2x2z_S_M2_vrr = PAZ*I_ERI_F3y_S_G2x2z_S_M2_vrr+WPZ*I_ERI_F3y_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_F3y_S_F2xz_S_M3_vrr;
      Double I_ERI_Gy3z_S_G2x2z_S_M2_vrr = PAY*I_ERI_F3z_S_G2x2z_S_M2_vrr+WPY*I_ERI_F3z_S_G2x2z_S_M3_vrr;
      Double I_ERI_G4z_S_G2x2z_S_M2_vrr = PAZ*I_ERI_F3z_S_G2x2z_S_M2_vrr+WPZ*I_ERI_F3z_S_G2x2z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G2x2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_F3z_S_F2xz_S_M3_vrr;
      Double I_ERI_G4x_S_Gx3y_S_M2_vrr = PAX*I_ERI_F3x_S_Gx3y_S_M2_vrr+WPX*I_ERI_F3x_S_Gx3y_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Gx3y_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx3y_S_M3_vrr+oned2k*I_ERI_F3x_S_F3y_S_M3_vrr;
      Double I_ERI_G3xy_S_Gx3y_S_M2_vrr = PAY*I_ERI_F3x_S_Gx3y_S_M2_vrr+WPY*I_ERI_F3x_S_Gx3y_S_M3_vrr+3*oned2k*I_ERI_F3x_S_Fx2y_S_M3_vrr;
      Double I_ERI_G3xz_S_Gx3y_S_M2_vrr = PAZ*I_ERI_F3x_S_Gx3y_S_M2_vrr+WPZ*I_ERI_F3x_S_Gx3y_S_M3_vrr;
      Double I_ERI_G2x2y_S_Gx3y_S_M2_vrr = PAY*I_ERI_F2xy_S_Gx3y_S_M2_vrr+WPY*I_ERI_F2xy_S_Gx3y_S_M3_vrr+oned2z*I_ERI_D2x_S_Gx3y_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Gx3y_S_M3_vrr+3*oned2k*I_ERI_F2xy_S_Fx2y_S_M3_vrr;
      Double I_ERI_Gx3y_S_Gx3y_S_M2_vrr = PAX*I_ERI_F3y_S_Gx3y_S_M2_vrr+WPX*I_ERI_F3y_S_Gx3y_S_M3_vrr+oned2k*I_ERI_F3y_S_F3y_S_M3_vrr;
      Double I_ERI_Gx3z_S_Gx3y_S_M2_vrr = PAX*I_ERI_F3z_S_Gx3y_S_M2_vrr+WPX*I_ERI_F3z_S_Gx3y_S_M3_vrr+oned2k*I_ERI_F3z_S_F3y_S_M3_vrr;
      Double I_ERI_G4y_S_Gx3y_S_M2_vrr = PAY*I_ERI_F3y_S_Gx3y_S_M2_vrr+WPY*I_ERI_F3y_S_Gx3y_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Gx3y_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx3y_S_M3_vrr+3*oned2k*I_ERI_F3y_S_Fx2y_S_M3_vrr;
      Double I_ERI_G3yz_S_Gx3y_S_M2_vrr = PAZ*I_ERI_F3y_S_Gx3y_S_M2_vrr+WPZ*I_ERI_F3y_S_Gx3y_S_M3_vrr;
      Double I_ERI_Gy3z_S_Gx3y_S_M2_vrr = PAY*I_ERI_F3z_S_Gx3y_S_M2_vrr+WPY*I_ERI_F3z_S_Gx3y_S_M3_vrr+3*oned2k*I_ERI_F3z_S_Fx2y_S_M3_vrr;
      Double I_ERI_G4z_S_Gx3y_S_M2_vrr = PAZ*I_ERI_F3z_S_Gx3y_S_M2_vrr+WPZ*I_ERI_F3z_S_Gx3y_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Gx3y_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx3y_S_M3_vrr;
      Double I_ERI_G4x_S_Gx2yz_S_M2_vrr = PAX*I_ERI_F3x_S_Gx2yz_S_M2_vrr+WPX*I_ERI_F3x_S_Gx2yz_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Gx2yz_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_F3x_S_F2yz_S_M3_vrr;
      Double I_ERI_G3xy_S_Gx2yz_S_M2_vrr = PAY*I_ERI_F3x_S_Gx2yz_S_M2_vrr+WPY*I_ERI_F3x_S_Gx2yz_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M3_vrr;
      Double I_ERI_G3xz_S_Gx2yz_S_M2_vrr = PAZ*I_ERI_F3x_S_Gx2yz_S_M2_vrr+WPZ*I_ERI_F3x_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_F3x_S_Fx2y_S_M3_vrr;
      Double I_ERI_G2x2y_S_Gx2yz_S_M2_vrr = PAY*I_ERI_F2xy_S_Gx2yz_S_M2_vrr+WPY*I_ERI_F2xy_S_Gx2yz_S_M3_vrr+oned2z*I_ERI_D2x_S_Gx2yz_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M3_vrr+2*oned2k*I_ERI_F2xy_S_Fxyz_S_M3_vrr;
      Double I_ERI_Gx3y_S_Gx2yz_S_M2_vrr = PAX*I_ERI_F3y_S_Gx2yz_S_M2_vrr+WPX*I_ERI_F3y_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_F3y_S_F2yz_S_M3_vrr;
      Double I_ERI_Gx3z_S_Gx2yz_S_M2_vrr = PAX*I_ERI_F3z_S_Gx2yz_S_M2_vrr+WPX*I_ERI_F3z_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_F3z_S_F2yz_S_M3_vrr;
      Double I_ERI_G4y_S_Gx2yz_S_M2_vrr = PAY*I_ERI_F3y_S_Gx2yz_S_M2_vrr+WPY*I_ERI_F3y_S_Gx2yz_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Gx2yz_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx2yz_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M3_vrr;
      Double I_ERI_G3yz_S_Gx2yz_S_M2_vrr = PAZ*I_ERI_F3y_S_Gx2yz_S_M2_vrr+WPZ*I_ERI_F3y_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_F3y_S_Fx2y_S_M3_vrr;
      Double I_ERI_Gy3z_S_Gx2yz_S_M2_vrr = PAY*I_ERI_F3z_S_Gx2yz_S_M2_vrr+WPY*I_ERI_F3z_S_Gx2yz_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M3_vrr;
      Double I_ERI_G4z_S_Gx2yz_S_M2_vrr = PAZ*I_ERI_F3z_S_Gx2yz_S_M2_vrr+WPZ*I_ERI_F3z_S_Gx2yz_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Gx2yz_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_F3z_S_Fx2y_S_M3_vrr;
      Double I_ERI_G4x_S_Gxy2z_S_M2_vrr = PAX*I_ERI_F3x_S_Gxy2z_S_M2_vrr+WPX*I_ERI_F3x_S_Gxy2z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Gxy2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_F3x_S_Fy2z_S_M3_vrr;
      Double I_ERI_G3xy_S_Gxy2z_S_M2_vrr = PAY*I_ERI_F3x_S_Gxy2z_S_M2_vrr+WPY*I_ERI_F3x_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_F3x_S_Fx2z_S_M3_vrr;
      Double I_ERI_G3xz_S_Gxy2z_S_M2_vrr = PAZ*I_ERI_F3x_S_Gxy2z_S_M2_vrr+WPZ*I_ERI_F3x_S_Gxy2z_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M3_vrr;
      Double I_ERI_G2x2y_S_Gxy2z_S_M2_vrr = PAY*I_ERI_F2xy_S_Gxy2z_S_M2_vrr+WPY*I_ERI_F2xy_S_Gxy2z_S_M3_vrr+oned2z*I_ERI_D2x_S_Gxy2z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_F2xy_S_Fx2z_S_M3_vrr;
      Double I_ERI_Gx3y_S_Gxy2z_S_M2_vrr = PAX*I_ERI_F3y_S_Gxy2z_S_M2_vrr+WPX*I_ERI_F3y_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_F3y_S_Fy2z_S_M3_vrr;
      Double I_ERI_Gx3z_S_Gxy2z_S_M2_vrr = PAX*I_ERI_F3z_S_Gxy2z_S_M2_vrr+WPX*I_ERI_F3z_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_F3z_S_Fy2z_S_M3_vrr;
      Double I_ERI_G4y_S_Gxy2z_S_M2_vrr = PAY*I_ERI_F3y_S_Gxy2z_S_M2_vrr+WPY*I_ERI_F3y_S_Gxy2z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Gxy2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_F3y_S_Fx2z_S_M3_vrr;
      Double I_ERI_G3yz_S_Gxy2z_S_M2_vrr = PAZ*I_ERI_F3y_S_Gxy2z_S_M2_vrr+WPZ*I_ERI_F3y_S_Gxy2z_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M3_vrr;
      Double I_ERI_Gy3z_S_Gxy2z_S_M2_vrr = PAY*I_ERI_F3z_S_Gxy2z_S_M2_vrr+WPY*I_ERI_F3z_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_F3z_S_Fx2z_S_M3_vrr;
      Double I_ERI_G4z_S_Gxy2z_S_M2_vrr = PAZ*I_ERI_F3z_S_Gxy2z_S_M2_vrr+WPZ*I_ERI_F3z_S_Gxy2z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Gxy2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Gxy2z_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M3_vrr;
      Double I_ERI_G4x_S_Gx3z_S_M2_vrr = PAX*I_ERI_F3x_S_Gx3z_S_M2_vrr+WPX*I_ERI_F3x_S_Gx3z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Gx3z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx3z_S_M3_vrr+oned2k*I_ERI_F3x_S_F3z_S_M3_vrr;
      Double I_ERI_G3xy_S_Gx3z_S_M2_vrr = PAY*I_ERI_F3x_S_Gx3z_S_M2_vrr+WPY*I_ERI_F3x_S_Gx3z_S_M3_vrr;
      Double I_ERI_G3xz_S_Gx3z_S_M2_vrr = PAZ*I_ERI_F3x_S_Gx3z_S_M2_vrr+WPZ*I_ERI_F3x_S_Gx3z_S_M3_vrr+3*oned2k*I_ERI_F3x_S_Fx2z_S_M3_vrr;
      Double I_ERI_G2x2y_S_Gx3z_S_M2_vrr = PAY*I_ERI_F2xy_S_Gx3z_S_M2_vrr+WPY*I_ERI_F2xy_S_Gx3z_S_M3_vrr+oned2z*I_ERI_D2x_S_Gx3z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Gx3z_S_M3_vrr;
      Double I_ERI_Gx3y_S_Gx3z_S_M2_vrr = PAX*I_ERI_F3y_S_Gx3z_S_M2_vrr+WPX*I_ERI_F3y_S_Gx3z_S_M3_vrr+oned2k*I_ERI_F3y_S_F3z_S_M3_vrr;
      Double I_ERI_Gx3z_S_Gx3z_S_M2_vrr = PAX*I_ERI_F3z_S_Gx3z_S_M2_vrr+WPX*I_ERI_F3z_S_Gx3z_S_M3_vrr+oned2k*I_ERI_F3z_S_F3z_S_M3_vrr;
      Double I_ERI_G4y_S_Gx3z_S_M2_vrr = PAY*I_ERI_F3y_S_Gx3z_S_M2_vrr+WPY*I_ERI_F3y_S_Gx3z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Gx3z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx3z_S_M3_vrr;
      Double I_ERI_G3yz_S_Gx3z_S_M2_vrr = PAZ*I_ERI_F3y_S_Gx3z_S_M2_vrr+WPZ*I_ERI_F3y_S_Gx3z_S_M3_vrr+3*oned2k*I_ERI_F3y_S_Fx2z_S_M3_vrr;
      Double I_ERI_Gy3z_S_Gx3z_S_M2_vrr = PAY*I_ERI_F3z_S_Gx3z_S_M2_vrr+WPY*I_ERI_F3z_S_Gx3z_S_M3_vrr;
      Double I_ERI_G4z_S_Gx3z_S_M2_vrr = PAZ*I_ERI_F3z_S_Gx3z_S_M2_vrr+WPZ*I_ERI_F3z_S_Gx3z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Gx3z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx3z_S_M3_vrr+3*oned2k*I_ERI_F3z_S_Fx2z_S_M3_vrr;
      Double I_ERI_G4x_S_G4y_S_M2_vrr = PAX*I_ERI_F3x_S_G4y_S_M2_vrr+WPX*I_ERI_F3x_S_G4y_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G4y_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G4y_S_M3_vrr;
      Double I_ERI_G3xy_S_G4y_S_M2_vrr = PAY*I_ERI_F3x_S_G4y_S_M2_vrr+WPY*I_ERI_F3x_S_G4y_S_M3_vrr+4*oned2k*I_ERI_F3x_S_F3y_S_M3_vrr;
      Double I_ERI_G3xz_S_G4y_S_M2_vrr = PAZ*I_ERI_F3x_S_G4y_S_M2_vrr+WPZ*I_ERI_F3x_S_G4y_S_M3_vrr;
      Double I_ERI_G2x2y_S_G4y_S_M2_vrr = PAY*I_ERI_F2xy_S_G4y_S_M2_vrr+WPY*I_ERI_F2xy_S_G4y_S_M3_vrr+oned2z*I_ERI_D2x_S_G4y_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G4y_S_M3_vrr+4*oned2k*I_ERI_F2xy_S_F3y_S_M3_vrr;
      Double I_ERI_Gx3y_S_G4y_S_M2_vrr = PAX*I_ERI_F3y_S_G4y_S_M2_vrr+WPX*I_ERI_F3y_S_G4y_S_M3_vrr;
      Double I_ERI_Gx3z_S_G4y_S_M2_vrr = PAX*I_ERI_F3z_S_G4y_S_M2_vrr+WPX*I_ERI_F3z_S_G4y_S_M3_vrr;
      Double I_ERI_G4y_S_G4y_S_M2_vrr = PAY*I_ERI_F3y_S_G4y_S_M2_vrr+WPY*I_ERI_F3y_S_G4y_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G4y_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G4y_S_M3_vrr+4*oned2k*I_ERI_F3y_S_F3y_S_M3_vrr;
      Double I_ERI_G3yz_S_G4y_S_M2_vrr = PAZ*I_ERI_F3y_S_G4y_S_M2_vrr+WPZ*I_ERI_F3y_S_G4y_S_M3_vrr;
      Double I_ERI_Gy3z_S_G4y_S_M2_vrr = PAY*I_ERI_F3z_S_G4y_S_M2_vrr+WPY*I_ERI_F3z_S_G4y_S_M3_vrr+4*oned2k*I_ERI_F3z_S_F3y_S_M3_vrr;
      Double I_ERI_G4z_S_G4y_S_M2_vrr = PAZ*I_ERI_F3z_S_G4y_S_M2_vrr+WPZ*I_ERI_F3z_S_G4y_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G4y_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G4y_S_M3_vrr;
      Double I_ERI_G4x_S_G3yz_S_M2_vrr = PAX*I_ERI_F3x_S_G3yz_S_M2_vrr+WPX*I_ERI_F3x_S_G3yz_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G3yz_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G3yz_S_M3_vrr;
      Double I_ERI_G3xy_S_G3yz_S_M2_vrr = PAY*I_ERI_F3x_S_G3yz_S_M2_vrr+WPY*I_ERI_F3x_S_G3yz_S_M3_vrr+3*oned2k*I_ERI_F3x_S_F2yz_S_M3_vrr;
      Double I_ERI_G3xz_S_G3yz_S_M2_vrr = PAZ*I_ERI_F3x_S_G3yz_S_M2_vrr+WPZ*I_ERI_F3x_S_G3yz_S_M3_vrr+oned2k*I_ERI_F3x_S_F3y_S_M3_vrr;
      Double I_ERI_G2x2y_S_G3yz_S_M2_vrr = PAY*I_ERI_F2xy_S_G3yz_S_M2_vrr+WPY*I_ERI_F2xy_S_G3yz_S_M3_vrr+oned2z*I_ERI_D2x_S_G3yz_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G3yz_S_M3_vrr+3*oned2k*I_ERI_F2xy_S_F2yz_S_M3_vrr;
      Double I_ERI_Gx3y_S_G3yz_S_M2_vrr = PAX*I_ERI_F3y_S_G3yz_S_M2_vrr+WPX*I_ERI_F3y_S_G3yz_S_M3_vrr;
      Double I_ERI_Gx3z_S_G3yz_S_M2_vrr = PAX*I_ERI_F3z_S_G3yz_S_M2_vrr+WPX*I_ERI_F3z_S_G3yz_S_M3_vrr;
      Double I_ERI_G4y_S_G3yz_S_M2_vrr = PAY*I_ERI_F3y_S_G3yz_S_M2_vrr+WPY*I_ERI_F3y_S_G3yz_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G3yz_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G3yz_S_M3_vrr+3*oned2k*I_ERI_F3y_S_F2yz_S_M3_vrr;
      Double I_ERI_G3yz_S_G3yz_S_M2_vrr = PAZ*I_ERI_F3y_S_G3yz_S_M2_vrr+WPZ*I_ERI_F3y_S_G3yz_S_M3_vrr+oned2k*I_ERI_F3y_S_F3y_S_M3_vrr;
      Double I_ERI_Gy3z_S_G3yz_S_M2_vrr = PAY*I_ERI_F3z_S_G3yz_S_M2_vrr+WPY*I_ERI_F3z_S_G3yz_S_M3_vrr+3*oned2k*I_ERI_F3z_S_F2yz_S_M3_vrr;
      Double I_ERI_G4z_S_G3yz_S_M2_vrr = PAZ*I_ERI_F3z_S_G3yz_S_M2_vrr+WPZ*I_ERI_F3z_S_G3yz_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G3yz_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G3yz_S_M3_vrr+oned2k*I_ERI_F3z_S_F3y_S_M3_vrr;
      Double I_ERI_G4x_S_G2y2z_S_M2_vrr = PAX*I_ERI_F3x_S_G2y2z_S_M2_vrr+WPX*I_ERI_F3x_S_G2y2z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G2y2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G2y2z_S_M3_vrr;
      Double I_ERI_G3xy_S_G2y2z_S_M2_vrr = PAY*I_ERI_F3x_S_G2y2z_S_M2_vrr+WPY*I_ERI_F3x_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_F3x_S_Fy2z_S_M3_vrr;
      Double I_ERI_G3xz_S_G2y2z_S_M2_vrr = PAZ*I_ERI_F3x_S_G2y2z_S_M2_vrr+WPZ*I_ERI_F3x_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_F3x_S_F2yz_S_M3_vrr;
      Double I_ERI_G2x2y_S_G2y2z_S_M2_vrr = PAY*I_ERI_F2xy_S_G2y2z_S_M2_vrr+WPY*I_ERI_F2xy_S_G2y2z_S_M3_vrr+oned2z*I_ERI_D2x_S_G2y2z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_F2xy_S_Fy2z_S_M3_vrr;
      Double I_ERI_Gx3y_S_G2y2z_S_M2_vrr = PAX*I_ERI_F3y_S_G2y2z_S_M2_vrr+WPX*I_ERI_F3y_S_G2y2z_S_M3_vrr;
      Double I_ERI_Gx3z_S_G2y2z_S_M2_vrr = PAX*I_ERI_F3z_S_G2y2z_S_M2_vrr+WPX*I_ERI_F3z_S_G2y2z_S_M3_vrr;
      Double I_ERI_G4y_S_G2y2z_S_M2_vrr = PAY*I_ERI_F3y_S_G2y2z_S_M2_vrr+WPY*I_ERI_F3y_S_G2y2z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G2y2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_F3y_S_Fy2z_S_M3_vrr;
      Double I_ERI_G3yz_S_G2y2z_S_M2_vrr = PAZ*I_ERI_F3y_S_G2y2z_S_M2_vrr+WPZ*I_ERI_F3y_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_F3y_S_F2yz_S_M3_vrr;
      Double I_ERI_Gy3z_S_G2y2z_S_M2_vrr = PAY*I_ERI_F3z_S_G2y2z_S_M2_vrr+WPY*I_ERI_F3z_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_F3z_S_Fy2z_S_M3_vrr;
      Double I_ERI_G4z_S_G2y2z_S_M2_vrr = PAZ*I_ERI_F3z_S_G2y2z_S_M2_vrr+WPZ*I_ERI_F3z_S_G2y2z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G2y2z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_F3z_S_F2yz_S_M3_vrr;
      Double I_ERI_G4x_S_Gy3z_S_M2_vrr = PAX*I_ERI_F3x_S_Gy3z_S_M2_vrr+WPX*I_ERI_F3x_S_Gy3z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_Gy3z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_Gy3z_S_M3_vrr;
      Double I_ERI_G3xy_S_Gy3z_S_M2_vrr = PAY*I_ERI_F3x_S_Gy3z_S_M2_vrr+WPY*I_ERI_F3x_S_Gy3z_S_M3_vrr+oned2k*I_ERI_F3x_S_F3z_S_M3_vrr;
      Double I_ERI_G3xz_S_Gy3z_S_M2_vrr = PAZ*I_ERI_F3x_S_Gy3z_S_M2_vrr+WPZ*I_ERI_F3x_S_Gy3z_S_M3_vrr+3*oned2k*I_ERI_F3x_S_Fy2z_S_M3_vrr;
      Double I_ERI_G2x2y_S_Gy3z_S_M2_vrr = PAY*I_ERI_F2xy_S_Gy3z_S_M2_vrr+WPY*I_ERI_F2xy_S_Gy3z_S_M3_vrr+oned2z*I_ERI_D2x_S_Gy3z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_Gy3z_S_M3_vrr+oned2k*I_ERI_F2xy_S_F3z_S_M3_vrr;
      Double I_ERI_Gx3y_S_Gy3z_S_M2_vrr = PAX*I_ERI_F3y_S_Gy3z_S_M2_vrr+WPX*I_ERI_F3y_S_Gy3z_S_M3_vrr;
      Double I_ERI_Gx3z_S_Gy3z_S_M2_vrr = PAX*I_ERI_F3z_S_Gy3z_S_M2_vrr+WPX*I_ERI_F3z_S_Gy3z_S_M3_vrr;
      Double I_ERI_G4y_S_Gy3z_S_M2_vrr = PAY*I_ERI_F3y_S_Gy3z_S_M2_vrr+WPY*I_ERI_F3y_S_Gy3z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_Gy3z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_Gy3z_S_M3_vrr+oned2k*I_ERI_F3y_S_F3z_S_M3_vrr;
      Double I_ERI_G3yz_S_Gy3z_S_M2_vrr = PAZ*I_ERI_F3y_S_Gy3z_S_M2_vrr+WPZ*I_ERI_F3y_S_Gy3z_S_M3_vrr+3*oned2k*I_ERI_F3y_S_Fy2z_S_M3_vrr;
      Double I_ERI_Gy3z_S_Gy3z_S_M2_vrr = PAY*I_ERI_F3z_S_Gy3z_S_M2_vrr+WPY*I_ERI_F3z_S_Gy3z_S_M3_vrr+oned2k*I_ERI_F3z_S_F3z_S_M3_vrr;
      Double I_ERI_G4z_S_Gy3z_S_M2_vrr = PAZ*I_ERI_F3z_S_Gy3z_S_M2_vrr+WPZ*I_ERI_F3z_S_Gy3z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_Gy3z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_Gy3z_S_M3_vrr+3*oned2k*I_ERI_F3z_S_Fy2z_S_M3_vrr;
      Double I_ERI_G4x_S_G4z_S_M2_vrr = PAX*I_ERI_F3x_S_G4z_S_M2_vrr+WPX*I_ERI_F3x_S_G4z_S_M3_vrr+3*oned2z*I_ERI_D2x_S_G4z_S_M2_vrr-3*rhod2zsq*I_ERI_D2x_S_G4z_S_M3_vrr;
      Double I_ERI_G3xy_S_G4z_S_M2_vrr = PAY*I_ERI_F3x_S_G4z_S_M2_vrr+WPY*I_ERI_F3x_S_G4z_S_M3_vrr;
      Double I_ERI_G3xz_S_G4z_S_M2_vrr = PAZ*I_ERI_F3x_S_G4z_S_M2_vrr+WPZ*I_ERI_F3x_S_G4z_S_M3_vrr+4*oned2k*I_ERI_F3x_S_F3z_S_M3_vrr;
      Double I_ERI_G2x2y_S_G4z_S_M2_vrr = PAY*I_ERI_F2xy_S_G4z_S_M2_vrr+WPY*I_ERI_F2xy_S_G4z_S_M3_vrr+oned2z*I_ERI_D2x_S_G4z_S_M2_vrr-rhod2zsq*I_ERI_D2x_S_G4z_S_M3_vrr;
      Double I_ERI_Gx3y_S_G4z_S_M2_vrr = PAX*I_ERI_F3y_S_G4z_S_M2_vrr+WPX*I_ERI_F3y_S_G4z_S_M3_vrr;
      Double I_ERI_Gx3z_S_G4z_S_M2_vrr = PAX*I_ERI_F3z_S_G4z_S_M2_vrr+WPX*I_ERI_F3z_S_G4z_S_M3_vrr;
      Double I_ERI_G4y_S_G4z_S_M2_vrr = PAY*I_ERI_F3y_S_G4z_S_M2_vrr+WPY*I_ERI_F3y_S_G4z_S_M3_vrr+3*oned2z*I_ERI_D2y_S_G4z_S_M2_vrr-3*rhod2zsq*I_ERI_D2y_S_G4z_S_M3_vrr;
      Double I_ERI_G3yz_S_G4z_S_M2_vrr = PAZ*I_ERI_F3y_S_G4z_S_M2_vrr+WPZ*I_ERI_F3y_S_G4z_S_M3_vrr+4*oned2k*I_ERI_F3y_S_F3z_S_M3_vrr;
      Double I_ERI_Gy3z_S_G4z_S_M2_vrr = PAY*I_ERI_F3z_S_G4z_S_M2_vrr+WPY*I_ERI_F3z_S_G4z_S_M3_vrr;
      Double I_ERI_G4z_S_G4z_S_M2_vrr = PAZ*I_ERI_F3z_S_G4z_S_M2_vrr+WPZ*I_ERI_F3z_S_G4z_S_M3_vrr+3*oned2z*I_ERI_D2z_S_G4z_S_M2_vrr-3*rhod2zsq*I_ERI_D2z_S_G4z_S_M3_vrr+4*oned2k*I_ERI_F3z_S_F3z_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M1_vrr = QCX*I_ERI_S_S_Px_S_M1_vrr+WQX*I_ERI_S_S_Px_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Dxy_S_M1_vrr = QCY*I_ERI_S_S_Px_S_M1_vrr+WQY*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_S_D2y_S_M1_vrr = QCY*I_ERI_S_S_Py_S_M1_vrr+WQY*I_ERI_S_S_Py_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_D2z_S_M1_vrr = QCZ*I_ERI_S_S_Pz_S_M1_vrr+WQZ*I_ERI_S_S_Pz_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M1_vrr = QCX*I_ERI_S_S_D2x_S_M1_vrr+WQX*I_ERI_S_S_D2x_S_M2_vrr+2*oned2e*I_ERI_S_S_Px_S_M1_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_S_F2xy_S_M1_vrr = QCY*I_ERI_S_S_D2x_S_M1_vrr+WQY*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_F2xz_S_M1_vrr = QCZ*I_ERI_S_S_D2x_S_M1_vrr+WQZ*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_Fx2y_S_M1_vrr = QCX*I_ERI_S_S_D2y_S_M1_vrr+WQX*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_Fxyz_S_M1_vrr = QCZ*I_ERI_S_S_Dxy_S_M1_vrr+WQZ*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_S_S_Fx2z_S_M1_vrr = QCX*I_ERI_S_S_D2z_S_M1_vrr+WQX*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_S_S_F3y_S_M1_vrr = QCY*I_ERI_S_S_D2y_S_M1_vrr+WQY*I_ERI_S_S_D2y_S_M2_vrr+2*oned2e*I_ERI_S_S_Py_S_M1_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_S_F2yz_S_M1_vrr = QCZ*I_ERI_S_S_D2y_S_M1_vrr+WQZ*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_Fy2z_S_M1_vrr = QCY*I_ERI_S_S_D2z_S_M1_vrr+WQY*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_S_S_F3z_S_M1_vrr = QCZ*I_ERI_S_S_D2z_S_M1_vrr+WQZ*I_ERI_S_S_D2z_S_M2_vrr+2*oned2e*I_ERI_S_S_Pz_S_M1_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M1_vrr = PAX*I_ERI_S_S_F3x_S_M1_vrr+WPX*I_ERI_S_S_F3x_S_M2_vrr+3*oned2k*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Py_S_F3x_S_M1_vrr = PAY*I_ERI_S_S_F3x_S_M1_vrr+WPY*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Pz_S_F3x_S_M1_vrr = PAZ*I_ERI_S_S_F3x_S_M1_vrr+WPZ*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Px_S_F2xy_S_M1_vrr = PAX*I_ERI_S_S_F2xy_S_M1_vrr+WPX*I_ERI_S_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Py_S_F2xy_S_M1_vrr = PAY*I_ERI_S_S_F2xy_S_M1_vrr+WPY*I_ERI_S_S_F2xy_S_M2_vrr+oned2k*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Pz_S_F2xy_S_M1_vrr = PAZ*I_ERI_S_S_F2xy_S_M1_vrr+WPZ*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_Px_S_F2xz_S_M1_vrr = PAX*I_ERI_S_S_F2xz_S_M1_vrr+WPX*I_ERI_S_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Py_S_F2xz_S_M1_vrr = PAY*I_ERI_S_S_F2xz_S_M1_vrr+WPY*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Pz_S_F2xz_S_M1_vrr = PAZ*I_ERI_S_S_F2xz_S_M1_vrr+WPZ*I_ERI_S_S_F2xz_S_M2_vrr+oned2k*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Px_S_Fx2y_S_M1_vrr = PAX*I_ERI_S_S_Fx2y_S_M1_vrr+WPX*I_ERI_S_S_Fx2y_S_M2_vrr+oned2k*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Py_S_Fx2y_S_M1_vrr = PAY*I_ERI_S_S_Fx2y_S_M1_vrr+WPY*I_ERI_S_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_S_S_Fx2y_S_M1_vrr+WPZ*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_Px_S_Fxyz_S_M1_vrr = PAX*I_ERI_S_S_Fxyz_S_M1_vrr+WPX*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Py_S_Fxyz_S_M1_vrr = PAY*I_ERI_S_S_Fxyz_S_M1_vrr+WPY*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_S_S_Fxyz_S_M1_vrr+WPZ*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Px_S_Fx2z_S_M1_vrr = PAX*I_ERI_S_S_Fx2z_S_M1_vrr+WPX*I_ERI_S_S_Fx2z_S_M2_vrr+oned2k*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Py_S_Fx2z_S_M1_vrr = PAY*I_ERI_S_S_Fx2z_S_M1_vrr+WPY*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_S_S_Fx2z_S_M1_vrr+WPZ*I_ERI_S_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Px_S_F3y_S_M1_vrr = PAX*I_ERI_S_S_F3y_S_M1_vrr+WPX*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Py_S_F3y_S_M1_vrr = PAY*I_ERI_S_S_F3y_S_M1_vrr+WPY*I_ERI_S_S_F3y_S_M2_vrr+3*oned2k*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Pz_S_F3y_S_M1_vrr = PAZ*I_ERI_S_S_F3y_S_M1_vrr+WPZ*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Px_S_F2yz_S_M1_vrr = PAX*I_ERI_S_S_F2yz_S_M1_vrr+WPX*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Py_S_F2yz_S_M1_vrr = PAY*I_ERI_S_S_F2yz_S_M1_vrr+WPY*I_ERI_S_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Pz_S_F2yz_S_M1_vrr = PAZ*I_ERI_S_S_F2yz_S_M1_vrr+WPZ*I_ERI_S_S_F2yz_S_M2_vrr+oned2k*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Px_S_Fy2z_S_M1_vrr = PAX*I_ERI_S_S_Fy2z_S_M1_vrr+WPX*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Py_S_Fy2z_S_M1_vrr = PAY*I_ERI_S_S_Fy2z_S_M1_vrr+WPY*I_ERI_S_S_Fy2z_S_M2_vrr+oned2k*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_S_S_Fy2z_S_M1_vrr+WPZ*I_ERI_S_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Px_S_F3z_S_M1_vrr = PAX*I_ERI_S_S_F3z_S_M1_vrr+WPX*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Py_S_F3z_S_M1_vrr = PAY*I_ERI_S_S_F3z_S_M1_vrr+WPY*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Pz_S_F3z_S_M1_vrr = PAZ*I_ERI_S_S_F3z_S_M1_vrr+WPZ*I_ERI_S_S_F3z_S_M2_vrr+3*oned2k*I_ERI_S_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       ************************************************************/
      Double I_ERI_S_S_G4x_S_M1_vrr = QCX*I_ERI_S_S_F3x_S_M1_vrr+WQX*I_ERI_S_S_F3x_S_M2_vrr+3*oned2e*I_ERI_S_S_D2x_S_M1_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_G3xy_S_M1_vrr = QCY*I_ERI_S_S_F3x_S_M1_vrr+WQY*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_S_S_G3xz_S_M1_vrr = QCZ*I_ERI_S_S_F3x_S_M1_vrr+WQZ*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_S_S_G2x2y_S_M1_vrr = QCY*I_ERI_S_S_F2xy_S_M1_vrr+WQY*I_ERI_S_S_F2xy_S_M2_vrr+oned2e*I_ERI_S_S_D2x_S_M1_vrr-rhod2esq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_G2xyz_S_M1_vrr = QCZ*I_ERI_S_S_F2xy_S_M1_vrr+WQZ*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_S_S_G2x2z_S_M1_vrr = QCZ*I_ERI_S_S_F2xz_S_M1_vrr+WQZ*I_ERI_S_S_F2xz_S_M2_vrr+oned2e*I_ERI_S_S_D2x_S_M1_vrr-rhod2esq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_S_S_Gx3y_S_M1_vrr = QCX*I_ERI_S_S_F3y_S_M1_vrr+WQX*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_S_S_Gx2yz_S_M1_vrr = QCZ*I_ERI_S_S_Fx2y_S_M1_vrr+WQZ*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_S_S_Gxy2z_S_M1_vrr = QCY*I_ERI_S_S_Fx2z_S_M1_vrr+WQY*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_S_S_Gx3z_S_M1_vrr = QCX*I_ERI_S_S_F3z_S_M1_vrr+WQX*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_S_S_G4y_S_M1_vrr = QCY*I_ERI_S_S_F3y_S_M1_vrr+WQY*I_ERI_S_S_F3y_S_M2_vrr+3*oned2e*I_ERI_S_S_D2y_S_M1_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_G3yz_S_M1_vrr = QCZ*I_ERI_S_S_F3y_S_M1_vrr+WQZ*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_S_S_G2y2z_S_M1_vrr = QCZ*I_ERI_S_S_F2yz_S_M1_vrr+WQZ*I_ERI_S_S_F2yz_S_M2_vrr+oned2e*I_ERI_S_S_D2y_S_M1_vrr-rhod2esq*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_S_S_Gy3z_S_M1_vrr = QCY*I_ERI_S_S_F3z_S_M1_vrr+WQY*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_S_S_G4z_S_M1_vrr = QCZ*I_ERI_S_S_F3z_S_M1_vrr+WQZ*I_ERI_S_S_F3z_S_M2_vrr+3*oned2e*I_ERI_S_S_D2z_S_M1_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M1_vrr = PAX*I_ERI_Px_S_F3x_S_M1_vrr+WPX*I_ERI_Px_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_D2y_S_F3x_S_M1_vrr = PAY*I_ERI_Py_S_F3x_S_M1_vrr+WPY*I_ERI_Py_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_D2z_S_F3x_S_M1_vrr = PAZ*I_ERI_Pz_S_F3x_S_M1_vrr+WPZ*I_ERI_Pz_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_D2x_S_F2xy_S_M1_vrr = PAX*I_ERI_Px_S_F2xy_S_M1_vrr+WPX*I_ERI_Px_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M2_vrr;
      Double I_ERI_D2y_S_F2xy_S_M1_vrr = PAY*I_ERI_Py_S_F2xy_S_M1_vrr+WPY*I_ERI_Py_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr+oned2k*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_D2z_S_F2xy_S_M1_vrr = PAZ*I_ERI_Pz_S_F2xy_S_M1_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_D2x_S_F2xz_S_M1_vrr = PAX*I_ERI_Px_S_F2xz_S_M1_vrr+WPX*I_ERI_Px_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M2_vrr;
      Double I_ERI_D2y_S_F2xz_S_M1_vrr = PAY*I_ERI_Py_S_F2xz_S_M1_vrr+WPY*I_ERI_Py_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_D2z_S_F2xz_S_M1_vrr = PAZ*I_ERI_Pz_S_F2xz_S_M1_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr+oned2k*I_ERI_Pz_S_D2x_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M1_vrr = PAX*I_ERI_Px_S_Fx2y_S_M1_vrr+WPX*I_ERI_Px_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr+oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M1_vrr = PAY*I_ERI_Py_S_Fx2y_S_M1_vrr+WPY*I_ERI_Py_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_M1_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2x_S_Fxyz_S_M1_vrr = PAX*I_ERI_Px_S_Fxyz_S_M1_vrr+WPX*I_ERI_Px_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Px_S_Dyz_S_M2_vrr;
      Double I_ERI_D2y_S_Fxyz_S_M1_vrr = PAY*I_ERI_Py_S_Fxyz_S_M1_vrr+WPY*I_ERI_Py_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Py_S_Dxz_S_M2_vrr;
      Double I_ERI_D2z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_Pz_S_Fxyz_S_M1_vrr+WPZ*I_ERI_Pz_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Pz_S_Dxy_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M1_vrr = PAX*I_ERI_Px_S_Fx2z_S_M1_vrr+WPX*I_ERI_Px_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr+oned2k*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M1_vrr = PAY*I_ERI_Py_S_Fx2z_S_M1_vrr+WPY*I_ERI_Py_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M2_vrr;
      Double I_ERI_D2x_S_F3y_S_M1_vrr = PAX*I_ERI_Px_S_F3y_S_M1_vrr+WPX*I_ERI_Px_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_D2y_S_F3y_S_M1_vrr = PAY*I_ERI_Py_S_F3y_S_M1_vrr+WPY*I_ERI_Py_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M2_vrr;
      Double I_ERI_D2z_S_F3y_S_M1_vrr = PAZ*I_ERI_Pz_S_F3y_S_M1_vrr+WPZ*I_ERI_Pz_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_D2x_S_F2yz_S_M1_vrr = PAX*I_ERI_Px_S_F2yz_S_M1_vrr+WPX*I_ERI_Px_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_D2y_S_F2yz_S_M1_vrr = PAY*I_ERI_Py_S_F2yz_S_M1_vrr+WPY*I_ERI_Py_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M2_vrr;
      Double I_ERI_D2z_S_F2yz_S_M1_vrr = PAZ*I_ERI_Pz_S_F2yz_S_M1_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr+oned2k*I_ERI_Pz_S_D2y_S_M2_vrr;
      Double I_ERI_D2x_S_Fy2z_S_M1_vrr = PAX*I_ERI_Px_S_Fy2z_S_M1_vrr+WPX*I_ERI_Px_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2y_S_Fy2z_S_M1_vrr = PAY*I_ERI_Py_S_Fy2z_S_M1_vrr+WPY*I_ERI_Py_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr+oned2k*I_ERI_Py_S_D2z_S_M2_vrr;
      Double I_ERI_D2z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Fy2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Dyz_S_M2_vrr;
      Double I_ERI_D2x_S_F3z_S_M1_vrr = PAX*I_ERI_Px_S_F3z_S_M1_vrr+WPX*I_ERI_Px_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_D2y_S_F3z_S_M1_vrr = PAY*I_ERI_Py_S_F3z_S_M1_vrr+WPY*I_ERI_Py_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_D2z_S_F3z_S_M1_vrr = PAZ*I_ERI_Pz_S_F3z_S_M1_vrr+WPZ*I_ERI_Pz_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       ************************************************************/
      Double I_ERI_Px_S_G4x_S_M1_vrr = PAX*I_ERI_S_S_G4x_S_M1_vrr+WPX*I_ERI_S_S_G4x_S_M2_vrr+4*oned2k*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Py_S_G4x_S_M1_vrr = PAY*I_ERI_S_S_G4x_S_M1_vrr+WPY*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_Pz_S_G4x_S_M1_vrr = PAZ*I_ERI_S_S_G4x_S_M1_vrr+WPZ*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_Px_S_G3xy_S_M1_vrr = PAX*I_ERI_S_S_G3xy_S_M1_vrr+WPX*I_ERI_S_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_Py_S_G3xy_S_M1_vrr = PAY*I_ERI_S_S_G3xy_S_M1_vrr+WPY*I_ERI_S_S_G3xy_S_M2_vrr+oned2k*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Pz_S_G3xy_S_M1_vrr = PAZ*I_ERI_S_S_G3xy_S_M1_vrr+WPZ*I_ERI_S_S_G3xy_S_M2_vrr;
      Double I_ERI_Px_S_G3xz_S_M1_vrr = PAX*I_ERI_S_S_G3xz_S_M1_vrr+WPX*I_ERI_S_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Py_S_G3xz_S_M1_vrr = PAY*I_ERI_S_S_G3xz_S_M1_vrr+WPY*I_ERI_S_S_G3xz_S_M2_vrr;
      Double I_ERI_Pz_S_G3xz_S_M1_vrr = PAZ*I_ERI_S_S_G3xz_S_M1_vrr+WPZ*I_ERI_S_S_G3xz_S_M2_vrr+oned2k*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_Px_S_G2x2y_S_M1_vrr = PAX*I_ERI_S_S_G2x2y_S_M1_vrr+WPX*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_Py_S_G2x2y_S_M1_vrr = PAY*I_ERI_S_S_G2x2y_S_M1_vrr+WPY*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_Pz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_S_S_G2x2y_S_M1_vrr+WPZ*I_ERI_S_S_G2x2y_S_M2_vrr;
      Double I_ERI_Px_S_G2xyz_S_M1_vrr = PAX*I_ERI_S_S_G2xyz_S_M1_vrr+WPX*I_ERI_S_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M2_vrr;
      Double I_ERI_Py_S_G2xyz_S_M1_vrr = PAY*I_ERI_S_S_G2xyz_S_M1_vrr+WPY*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Pz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_S_S_G2xyz_S_M1_vrr+WPZ*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_Px_S_G2x2z_S_M1_vrr = PAX*I_ERI_S_S_G2x2z_S_M1_vrr+WPX*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Py_S_G2x2z_S_M1_vrr = PAY*I_ERI_S_S_G2x2z_S_M1_vrr+WPY*I_ERI_S_S_G2x2z_S_M2_vrr;
      Double I_ERI_Pz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_S_S_G2x2z_S_M1_vrr+WPZ*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_Px_S_Gx3y_S_M1_vrr = PAX*I_ERI_S_S_Gx3y_S_M1_vrr+WPX*I_ERI_S_S_Gx3y_S_M2_vrr+oned2k*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Py_S_Gx3y_S_M1_vrr = PAY*I_ERI_S_S_Gx3y_S_M1_vrr+WPY*I_ERI_S_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_Pz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_S_S_Gx3y_S_M1_vrr+WPZ*I_ERI_S_S_Gx3y_S_M2_vrr;
      Double I_ERI_Px_S_Gx2yz_S_M1_vrr = PAX*I_ERI_S_S_Gx2yz_S_M1_vrr+WPX*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Py_S_Gx2yz_S_M1_vrr = PAY*I_ERI_S_S_Gx2yz_S_M1_vrr+WPY*I_ERI_S_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M2_vrr;
      Double I_ERI_Pz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_S_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_Px_S_Gxy2z_S_M1_vrr = PAX*I_ERI_S_S_Gxy2z_S_M1_vrr+WPX*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Py_S_Gxy2z_S_M1_vrr = PAY*I_ERI_S_S_Gxy2z_S_M1_vrr+WPY*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Pz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_S_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_S_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M2_vrr;
      Double I_ERI_Px_S_Gx3z_S_M1_vrr = PAX*I_ERI_S_S_Gx3z_S_M1_vrr+WPX*I_ERI_S_S_Gx3z_S_M2_vrr+oned2k*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Py_S_Gx3z_S_M1_vrr = PAY*I_ERI_S_S_Gx3z_S_M1_vrr+WPY*I_ERI_S_S_Gx3z_S_M2_vrr;
      Double I_ERI_Pz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_S_S_Gx3z_S_M1_vrr+WPZ*I_ERI_S_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_Px_S_G4y_S_M1_vrr = PAX*I_ERI_S_S_G4y_S_M1_vrr+WPX*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_Py_S_G4y_S_M1_vrr = PAY*I_ERI_S_S_G4y_S_M1_vrr+WPY*I_ERI_S_S_G4y_S_M2_vrr+4*oned2k*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Pz_S_G4y_S_M1_vrr = PAZ*I_ERI_S_S_G4y_S_M1_vrr+WPZ*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_Px_S_G3yz_S_M1_vrr = PAX*I_ERI_S_S_G3yz_S_M1_vrr+WPX*I_ERI_S_S_G3yz_S_M2_vrr;
      Double I_ERI_Py_S_G3yz_S_M1_vrr = PAY*I_ERI_S_S_G3yz_S_M1_vrr+WPY*I_ERI_S_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Pz_S_G3yz_S_M1_vrr = PAZ*I_ERI_S_S_G3yz_S_M1_vrr+WPZ*I_ERI_S_S_G3yz_S_M2_vrr+oned2k*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Px_S_G2y2z_S_M1_vrr = PAX*I_ERI_S_S_G2y2z_S_M1_vrr+WPX*I_ERI_S_S_G2y2z_S_M2_vrr;
      Double I_ERI_Py_S_G2y2z_S_M1_vrr = PAY*I_ERI_S_S_G2y2z_S_M1_vrr+WPY*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Pz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_S_S_G2y2z_S_M1_vrr+WPZ*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Px_S_Gy3z_S_M1_vrr = PAX*I_ERI_S_S_Gy3z_S_M1_vrr+WPX*I_ERI_S_S_Gy3z_S_M2_vrr;
      Double I_ERI_Py_S_Gy3z_S_M1_vrr = PAY*I_ERI_S_S_Gy3z_S_M1_vrr+WPY*I_ERI_S_S_Gy3z_S_M2_vrr+oned2k*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Pz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_S_S_Gy3z_S_M1_vrr+WPZ*I_ERI_S_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Px_S_G4z_S_M1_vrr = PAX*I_ERI_S_S_G4z_S_M1_vrr+WPX*I_ERI_S_S_G4z_S_M2_vrr;
      Double I_ERI_Py_S_G4z_S_M1_vrr = PAY*I_ERI_S_S_G4z_S_M1_vrr+WPY*I_ERI_S_S_G4z_S_M2_vrr;
      Double I_ERI_Pz_S_G4z_S_M1_vrr = PAZ*I_ERI_S_S_G4z_S_M1_vrr+WPZ*I_ERI_S_S_G4z_S_M2_vrr+4*oned2k*I_ERI_S_S_F3z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 20 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_M1_vrr = PAX*I_ERI_D2x_S_F3x_S_M1_vrr+WPX*I_ERI_D2x_S_F3x_S_M2_vrr+2*oned2z*I_ERI_Px_S_F3x_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F3x_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xy_S_F3x_S_M1_vrr = PAY*I_ERI_D2x_S_F3x_S_M1_vrr+WPY*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_F2xz_S_F3x_S_M1_vrr = PAZ*I_ERI_D2x_S_F3x_S_M1_vrr+WPZ*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3x_S_M1_vrr = PAX*I_ERI_D2y_S_F3x_S_M1_vrr+WPX*I_ERI_D2y_S_F3x_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3x_S_M1_vrr = PAX*I_ERI_D2z_S_F3x_S_M1_vrr+WPX*I_ERI_D2z_S_F3x_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3y_S_F3x_S_M1_vrr = PAY*I_ERI_D2y_S_F3x_S_M1_vrr+WPY*I_ERI_D2y_S_F3x_S_M2_vrr+2*oned2z*I_ERI_Py_S_F3x_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F3x_S_M2_vrr;
      Double I_ERI_F2yz_S_F3x_S_M1_vrr = PAZ*I_ERI_D2y_S_F3x_S_M1_vrr+WPZ*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_F3z_S_F3x_S_M1_vrr = PAZ*I_ERI_D2z_S_F3x_S_M1_vrr+WPZ*I_ERI_D2z_S_F3x_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F3x_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F3x_S_M2_vrr;
      Double I_ERI_F3x_S_F2xy_S_M1_vrr = PAX*I_ERI_D2x_S_F2xy_S_M1_vrr+WPX*I_ERI_D2x_S_F2xy_S_M2_vrr+2*oned2z*I_ERI_Px_S_F2xy_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_F2xy_S_F2xy_S_M1_vrr = PAY*I_ERI_D2x_S_F2xy_S_M1_vrr+WPY*I_ERI_D2x_S_F2xy_S_M2_vrr+oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xz_S_F2xy_S_M1_vrr = PAZ*I_ERI_D2x_S_F2xy_S_M1_vrr+WPZ*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2xy_S_M1_vrr = PAX*I_ERI_D2y_S_F2xy_S_M1_vrr+WPX*I_ERI_D2y_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2xy_S_M1_vrr = PAX*I_ERI_D2z_S_F2xy_S_M1_vrr+WPX*I_ERI_D2z_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxy_S_M2_vrr;
      Double I_ERI_F3y_S_F2xy_S_M1_vrr = PAY*I_ERI_D2y_S_F2xy_S_M1_vrr+WPY*I_ERI_D2y_S_F2xy_S_M2_vrr+2*oned2z*I_ERI_Py_S_F2xy_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F2xy_S_M2_vrr+oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_F2yz_S_F2xy_S_M1_vrr = PAZ*I_ERI_D2y_S_F2xy_S_M1_vrr+WPZ*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_F3z_S_F2xy_S_M1_vrr = PAZ*I_ERI_D2z_S_F2xy_S_M1_vrr+WPZ*I_ERI_D2z_S_F2xy_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F2xy_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xy_S_M2_vrr;
      Double I_ERI_F3x_S_F2xz_S_M1_vrr = PAX*I_ERI_D2x_S_F2xz_S_M1_vrr+WPX*I_ERI_D2x_S_F2xz_S_M2_vrr+2*oned2z*I_ERI_Px_S_F2xz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_F2xy_S_F2xz_S_M1_vrr = PAY*I_ERI_D2x_S_F2xz_S_M1_vrr+WPY*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_F2xz_S_F2xz_S_M1_vrr = PAZ*I_ERI_D2x_S_F2xz_S_M1_vrr+WPZ*I_ERI_D2x_S_F2xz_S_M2_vrr+oned2k*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2xz_S_M1_vrr = PAX*I_ERI_D2y_S_F2xz_S_M1_vrr+WPX*I_ERI_D2y_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2xz_S_M1_vrr = PAX*I_ERI_D2z_S_F2xz_S_M1_vrr+WPX*I_ERI_D2z_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M2_vrr;
      Double I_ERI_F3y_S_F2xz_S_M1_vrr = PAY*I_ERI_D2y_S_F2xz_S_M1_vrr+WPY*I_ERI_D2y_S_F2xz_S_M2_vrr+2*oned2z*I_ERI_Py_S_F2xz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F2xz_S_M2_vrr;
      Double I_ERI_F2yz_S_F2xz_S_M1_vrr = PAZ*I_ERI_D2y_S_F2xz_S_M1_vrr+WPZ*I_ERI_D2y_S_F2xz_S_M2_vrr+oned2k*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_F3z_S_F2xz_S_M1_vrr = PAZ*I_ERI_D2z_S_F2xz_S_M1_vrr+WPZ*I_ERI_D2z_S_F2xz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F2xz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xz_S_M2_vrr+oned2k*I_ERI_D2z_S_D2x_S_M2_vrr;
      Double I_ERI_F3x_S_Fx2y_S_M1_vrr = PAX*I_ERI_D2x_S_Fx2y_S_M1_vrr+WPX*I_ERI_D2x_S_Fx2y_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fx2y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2y_S_M2_vrr+oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_M1_vrr = PAY*I_ERI_D2x_S_Fx2y_S_M1_vrr+WPY*I_ERI_D2x_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_F2xz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_D2x_S_Fx2y_S_M1_vrr+WPZ*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fx2y_S_M1_vrr = PAX*I_ERI_D2y_S_Fx2y_S_M1_vrr+WPX*I_ERI_D2y_S_Fx2y_S_M2_vrr+oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fx2y_S_M1_vrr = PAX*I_ERI_D2z_S_Fx2y_S_M1_vrr+WPX*I_ERI_D2z_S_Fx2y_S_M2_vrr+oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3y_S_Fx2y_S_M1_vrr = PAY*I_ERI_D2y_S_Fx2y_S_M1_vrr+WPY*I_ERI_D2y_S_Fx2y_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fx2y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_F2yz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_D2y_S_Fx2y_S_M1_vrr+WPZ*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_D2z_S_Fx2y_S_M1_vrr+WPZ*I_ERI_D2z_S_Fx2y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fx2y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3x_S_Fxyz_S_M1_vrr = PAX*I_ERI_D2x_S_Fxyz_S_M1_vrr+WPX*I_ERI_D2x_S_Fxyz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fxyz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_F2xy_S_Fxyz_S_M1_vrr = PAY*I_ERI_D2x_S_Fxyz_S_M1_vrr+WPY*I_ERI_D2x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_F2xz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_D2x_S_Fxyz_S_M1_vrr+WPZ*I_ERI_D2x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fxyz_S_M1_vrr = PAX*I_ERI_D2y_S_Fxyz_S_M1_vrr+WPX*I_ERI_D2y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fxyz_S_M1_vrr = PAX*I_ERI_D2z_S_Fxyz_S_M1_vrr+WPX*I_ERI_D2z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2z_S_Dyz_S_M2_vrr;
      Double I_ERI_F3y_S_Fxyz_S_M1_vrr = PAY*I_ERI_D2y_S_Fxyz_S_M1_vrr+WPY*I_ERI_D2y_S_Fxyz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fxyz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_F2yz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_D2y_S_Fxyz_S_M1_vrr+WPZ*I_ERI_D2y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_F3z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_D2z_S_Fxyz_S_M1_vrr+WPZ*I_ERI_D2z_S_Fxyz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fxyz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fxyz_S_M2_vrr+oned2k*I_ERI_D2z_S_Dxy_S_M2_vrr;
      Double I_ERI_F3x_S_Fx2z_S_M1_vrr = PAX*I_ERI_D2x_S_Fx2z_S_M1_vrr+WPX*I_ERI_D2x_S_Fx2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fx2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2z_S_M2_vrr+oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_M1_vrr = PAY*I_ERI_D2x_S_Fx2z_S_M1_vrr+WPY*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2xz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_D2x_S_Fx2z_S_M1_vrr+WPZ*I_ERI_D2x_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fx2z_S_M1_vrr = PAX*I_ERI_D2y_S_Fx2z_S_M1_vrr+WPX*I_ERI_D2y_S_Fx2z_S_M2_vrr+oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fx2z_S_M1_vrr = PAX*I_ERI_D2z_S_Fx2z_S_M1_vrr+WPX*I_ERI_D2z_S_Fx2z_S_M2_vrr+oned2k*I_ERI_D2z_S_D2z_S_M2_vrr;
      Double I_ERI_F3y_S_Fx2z_S_M1_vrr = PAY*I_ERI_D2y_S_Fx2z_S_M1_vrr+WPY*I_ERI_D2y_S_Fx2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fx2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2yz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_D2y_S_Fx2z_S_M1_vrr+WPZ*I_ERI_D2y_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M2_vrr;
      Double I_ERI_F3z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_D2z_S_Fx2z_S_M1_vrr+WPZ*I_ERI_D2z_S_Fx2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fx2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M2_vrr;
      Double I_ERI_F3x_S_F3y_S_M1_vrr = PAX*I_ERI_D2x_S_F3y_S_M1_vrr+WPX*I_ERI_D2x_S_F3y_S_M2_vrr+2*oned2z*I_ERI_Px_S_F3y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F3y_S_M2_vrr;
      Double I_ERI_F2xy_S_F3y_S_M1_vrr = PAY*I_ERI_D2x_S_F3y_S_M1_vrr+WPY*I_ERI_D2x_S_F3y_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_F2xz_S_F3y_S_M1_vrr = PAZ*I_ERI_D2x_S_F3y_S_M1_vrr+WPZ*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3y_S_M1_vrr = PAX*I_ERI_D2y_S_F3y_S_M1_vrr+WPX*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3y_S_M1_vrr = PAX*I_ERI_D2z_S_F3y_S_M1_vrr+WPX*I_ERI_D2z_S_F3y_S_M2_vrr;
      Double I_ERI_F3y_S_F3y_S_M1_vrr = PAY*I_ERI_D2y_S_F3y_S_M1_vrr+WPY*I_ERI_D2y_S_F3y_S_M2_vrr+2*oned2z*I_ERI_Py_S_F3y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F3y_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_F2yz_S_F3y_S_M1_vrr = PAZ*I_ERI_D2y_S_F3y_S_M1_vrr+WPZ*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_F3z_S_F3y_S_M1_vrr = PAZ*I_ERI_D2z_S_F3y_S_M1_vrr+WPZ*I_ERI_D2z_S_F3y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F3y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F3y_S_M2_vrr;
      Double I_ERI_F3x_S_F2yz_S_M1_vrr = PAX*I_ERI_D2x_S_F2yz_S_M1_vrr+WPX*I_ERI_D2x_S_F2yz_S_M2_vrr+2*oned2z*I_ERI_Px_S_F2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F2yz_S_M2_vrr;
      Double I_ERI_F2xy_S_F2yz_S_M1_vrr = PAY*I_ERI_D2x_S_F2yz_S_M1_vrr+WPY*I_ERI_D2x_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_F2xz_S_F2yz_S_M1_vrr = PAZ*I_ERI_D2x_S_F2yz_S_M1_vrr+WPZ*I_ERI_D2x_S_F2yz_S_M2_vrr+oned2k*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_F2yz_S_M1_vrr = PAX*I_ERI_D2y_S_F2yz_S_M1_vrr+WPX*I_ERI_D2y_S_F2yz_S_M2_vrr;
      Double I_ERI_Fx2z_S_F2yz_S_M1_vrr = PAX*I_ERI_D2z_S_F2yz_S_M1_vrr+WPX*I_ERI_D2z_S_F2yz_S_M2_vrr;
      Double I_ERI_F3y_S_F2yz_S_M1_vrr = PAY*I_ERI_D2y_S_F2yz_S_M1_vrr+WPY*I_ERI_D2y_S_F2yz_S_M2_vrr+2*oned2z*I_ERI_Py_S_F2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_F2yz_S_F2yz_S_M1_vrr = PAZ*I_ERI_D2y_S_F2yz_S_M1_vrr+WPZ*I_ERI_D2y_S_F2yz_S_M2_vrr+oned2k*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_F3z_S_F2yz_S_M1_vrr = PAZ*I_ERI_D2z_S_F2yz_S_M1_vrr+WPZ*I_ERI_D2z_S_F2yz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F2yz_S_M2_vrr+oned2k*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3x_S_Fy2z_S_M1_vrr = PAX*I_ERI_D2x_S_Fy2z_S_M1_vrr+WPX*I_ERI_D2x_S_Fy2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Fy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Fy2z_S_M2_vrr;
      Double I_ERI_F2xy_S_Fy2z_S_M1_vrr = PAY*I_ERI_D2x_S_Fy2z_S_M1_vrr+WPY*I_ERI_D2x_S_Fy2z_S_M2_vrr+oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_F2xz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_D2x_S_Fy2z_S_M1_vrr+WPZ*I_ERI_D2x_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M2_vrr;
      Double I_ERI_Fx2y_S_Fy2z_S_M1_vrr = PAX*I_ERI_D2y_S_Fy2z_S_M1_vrr+WPX*I_ERI_D2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Fy2z_S_M1_vrr = PAX*I_ERI_D2z_S_Fy2z_S_M1_vrr+WPX*I_ERI_D2z_S_Fy2z_S_M2_vrr;
      Double I_ERI_F3y_S_Fy2z_S_M1_vrr = PAY*I_ERI_D2y_S_Fy2z_S_M1_vrr+WPY*I_ERI_D2y_S_Fy2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Fy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Fy2z_S_M2_vrr+oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_F2yz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_D2y_S_Fy2z_S_M1_vrr+WPZ*I_ERI_D2y_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_F3z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_D2z_S_Fy2z_S_M1_vrr+WPZ*I_ERI_D2z_S_Fy2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Fy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M2_vrr;
      Double I_ERI_F3x_S_F3z_S_M1_vrr = PAX*I_ERI_D2x_S_F3z_S_M1_vrr+WPX*I_ERI_D2x_S_F3z_S_M2_vrr+2*oned2z*I_ERI_Px_S_F3z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_F3z_S_M2_vrr;
      Double I_ERI_F2xy_S_F3z_S_M1_vrr = PAY*I_ERI_D2x_S_F3z_S_M1_vrr+WPY*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_F2xz_S_F3z_S_M1_vrr = PAZ*I_ERI_D2x_S_F3z_S_M1_vrr+WPZ*I_ERI_D2x_S_F3z_S_M2_vrr+3*oned2k*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2y_S_F3z_S_M1_vrr = PAX*I_ERI_D2y_S_F3z_S_M1_vrr+WPX*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_Fx2z_S_F3z_S_M1_vrr = PAX*I_ERI_D2z_S_F3z_S_M1_vrr+WPX*I_ERI_D2z_S_F3z_S_M2_vrr;
      Double I_ERI_F3y_S_F3z_S_M1_vrr = PAY*I_ERI_D2y_S_F3z_S_M1_vrr+WPY*I_ERI_D2y_S_F3z_S_M2_vrr+2*oned2z*I_ERI_Py_S_F3z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_F3z_S_M2_vrr;
      Double I_ERI_F2yz_S_F3z_S_M1_vrr = PAZ*I_ERI_D2y_S_F3z_S_M1_vrr+WPZ*I_ERI_D2y_S_F3z_S_M2_vrr+3*oned2k*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_F3z_S_F3z_S_M1_vrr = PAZ*I_ERI_D2z_S_F3z_S_M1_vrr+WPZ*I_ERI_D2z_S_F3z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_F3z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_F3z_S_M2_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 45 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_M1_vrr = PAX*I_ERI_Px_S_G4x_S_M1_vrr+WPX*I_ERI_Px_S_G4x_S_M2_vrr+oned2z*I_ERI_S_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M2_vrr+4*oned2k*I_ERI_Px_S_F3x_S_M2_vrr;
      Double I_ERI_D2y_S_G4x_S_M1_vrr = PAY*I_ERI_Py_S_G4x_S_M1_vrr+WPY*I_ERI_Py_S_G4x_S_M2_vrr+oned2z*I_ERI_S_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_D2z_S_G4x_S_M1_vrr = PAZ*I_ERI_Pz_S_G4x_S_M1_vrr+WPZ*I_ERI_Pz_S_G4x_S_M2_vrr+oned2z*I_ERI_S_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_D2x_S_G3xy_S_M1_vrr = PAX*I_ERI_Px_S_G3xy_S_M1_vrr+WPX*I_ERI_Px_S_G3xy_S_M2_vrr+oned2z*I_ERI_S_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_Px_S_F2xy_S_M2_vrr;
      Double I_ERI_D2y_S_G3xy_S_M1_vrr = PAY*I_ERI_Py_S_G3xy_S_M1_vrr+WPY*I_ERI_Py_S_G3xy_S_M2_vrr+oned2z*I_ERI_S_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M2_vrr+oned2k*I_ERI_Py_S_F3x_S_M2_vrr;
      Double I_ERI_D2z_S_G3xy_S_M1_vrr = PAZ*I_ERI_Pz_S_G3xy_S_M1_vrr+WPZ*I_ERI_Pz_S_G3xy_S_M2_vrr+oned2z*I_ERI_S_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M2_vrr;
      Double I_ERI_D2x_S_G3xz_S_M1_vrr = PAX*I_ERI_Px_S_G3xz_S_M1_vrr+WPX*I_ERI_Px_S_G3xz_S_M2_vrr+oned2z*I_ERI_S_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_Px_S_F2xz_S_M2_vrr;
      Double I_ERI_D2y_S_G3xz_S_M1_vrr = PAY*I_ERI_Py_S_G3xz_S_M1_vrr+WPY*I_ERI_Py_S_G3xz_S_M2_vrr+oned2z*I_ERI_S_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M2_vrr;
      Double I_ERI_D2z_S_G3xz_S_M1_vrr = PAZ*I_ERI_Pz_S_G3xz_S_M1_vrr+WPZ*I_ERI_Pz_S_G3xz_S_M2_vrr+oned2z*I_ERI_S_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M2_vrr+oned2k*I_ERI_Pz_S_F3x_S_M2_vrr;
      Double I_ERI_D2x_S_G2x2y_S_M1_vrr = PAX*I_ERI_Px_S_G2x2y_S_M1_vrr+WPX*I_ERI_Px_S_G2x2y_S_M2_vrr+oned2z*I_ERI_S_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2y_S_G2x2y_S_M1_vrr = PAY*I_ERI_Py_S_G2x2y_S_M1_vrr+WPY*I_ERI_Py_S_G2x2y_S_M2_vrr+oned2z*I_ERI_S_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M2_vrr;
      Double I_ERI_D2z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_Pz_S_G2x2y_S_M1_vrr+WPZ*I_ERI_Pz_S_G2x2y_S_M2_vrr+oned2z*I_ERI_S_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M2_vrr;
      Double I_ERI_D2x_S_G2xyz_S_M1_vrr = PAX*I_ERI_Px_S_G2xyz_S_M1_vrr+WPX*I_ERI_Px_S_G2xyz_S_M2_vrr+oned2z*I_ERI_S_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M2_vrr;
      Double I_ERI_D2y_S_G2xyz_S_M1_vrr = PAY*I_ERI_Py_S_G2xyz_S_M1_vrr+WPY*I_ERI_Py_S_G2xyz_S_M2_vrr+oned2z*I_ERI_S_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_Py_S_F2xz_S_M2_vrr;
      Double I_ERI_D2z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_Pz_S_G2xyz_S_M1_vrr+WPZ*I_ERI_Pz_S_G2xyz_S_M2_vrr+oned2z*I_ERI_S_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_Pz_S_F2xy_S_M2_vrr;
      Double I_ERI_D2x_S_G2x2z_S_M1_vrr = PAX*I_ERI_Px_S_G2x2z_S_M1_vrr+WPX*I_ERI_Px_S_G2x2z_S_M2_vrr+oned2z*I_ERI_S_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2y_S_G2x2z_S_M1_vrr = PAY*I_ERI_Py_S_G2x2z_S_M1_vrr+WPY*I_ERI_Py_S_G2x2z_S_M2_vrr+oned2z*I_ERI_S_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M2_vrr;
      Double I_ERI_D2z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_Pz_S_G2x2z_S_M1_vrr+WPZ*I_ERI_Pz_S_G2x2z_S_M2_vrr+oned2z*I_ERI_S_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M2_vrr;
      Double I_ERI_D2x_S_Gx3y_S_M1_vrr = PAX*I_ERI_Px_S_Gx3y_S_M1_vrr+WPX*I_ERI_Px_S_Gx3y_S_M2_vrr+oned2z*I_ERI_S_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M2_vrr+oned2k*I_ERI_Px_S_F3y_S_M2_vrr;
      Double I_ERI_D2y_S_Gx3y_S_M1_vrr = PAY*I_ERI_Py_S_Gx3y_S_M1_vrr+WPY*I_ERI_Py_S_Gx3y_S_M2_vrr+oned2z*I_ERI_S_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_Py_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_Pz_S_Gx3y_S_M1_vrr+WPZ*I_ERI_Pz_S_Gx3y_S_M2_vrr+oned2z*I_ERI_S_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M2_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_M1_vrr = PAX*I_ERI_Px_S_Gx2yz_S_M1_vrr+WPX*I_ERI_Px_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_Px_S_F2yz_S_M2_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_M1_vrr = PAY*I_ERI_Py_S_Gx2yz_S_M1_vrr+WPY*I_ERI_Py_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_Py_S_Fxyz_S_M2_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_Pz_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_Pz_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_Pz_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_M1_vrr = PAX*I_ERI_Px_S_Gxy2z_S_M1_vrr+WPX*I_ERI_Px_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Px_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_M1_vrr = PAY*I_ERI_Py_S_Gxy2z_S_M1_vrr+WPY*I_ERI_Py_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Py_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Fxyz_S_M2_vrr;
      Double I_ERI_D2x_S_Gx3z_S_M1_vrr = PAX*I_ERI_Px_S_Gx3z_S_M1_vrr+WPX*I_ERI_Px_S_Gx3z_S_M2_vrr+oned2z*I_ERI_S_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M2_vrr+oned2k*I_ERI_Px_S_F3z_S_M2_vrr;
      Double I_ERI_D2y_S_Gx3z_S_M1_vrr = PAY*I_ERI_Py_S_Gx3z_S_M1_vrr+WPY*I_ERI_Py_S_Gx3z_S_M2_vrr+oned2z*I_ERI_S_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M2_vrr;
      Double I_ERI_D2z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_Pz_S_Gx3z_S_M1_vrr+WPZ*I_ERI_Pz_S_Gx3z_S_M2_vrr+oned2z*I_ERI_S_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2x_S_G4y_S_M1_vrr = PAX*I_ERI_Px_S_G4y_S_M1_vrr+WPX*I_ERI_Px_S_G4y_S_M2_vrr+oned2z*I_ERI_S_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_D2y_S_G4y_S_M1_vrr = PAY*I_ERI_Py_S_G4y_S_M1_vrr+WPY*I_ERI_Py_S_G4y_S_M2_vrr+oned2z*I_ERI_S_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M2_vrr+4*oned2k*I_ERI_Py_S_F3y_S_M2_vrr;
      Double I_ERI_D2z_S_G4y_S_M1_vrr = PAZ*I_ERI_Pz_S_G4y_S_M1_vrr+WPZ*I_ERI_Pz_S_G4y_S_M2_vrr+oned2z*I_ERI_S_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_D2x_S_G3yz_S_M1_vrr = PAX*I_ERI_Px_S_G3yz_S_M1_vrr+WPX*I_ERI_Px_S_G3yz_S_M2_vrr+oned2z*I_ERI_S_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M2_vrr;
      Double I_ERI_D2y_S_G3yz_S_M1_vrr = PAY*I_ERI_Py_S_G3yz_S_M1_vrr+WPY*I_ERI_Py_S_G3yz_S_M2_vrr+oned2z*I_ERI_S_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_Py_S_F2yz_S_M2_vrr;
      Double I_ERI_D2z_S_G3yz_S_M1_vrr = PAZ*I_ERI_Pz_S_G3yz_S_M1_vrr+WPZ*I_ERI_Pz_S_G3yz_S_M2_vrr+oned2z*I_ERI_S_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M2_vrr+oned2k*I_ERI_Pz_S_F3y_S_M2_vrr;
      Double I_ERI_D2x_S_G2y2z_S_M1_vrr = PAX*I_ERI_Px_S_G2y2z_S_M1_vrr+WPX*I_ERI_Px_S_G2y2z_S_M2_vrr+oned2z*I_ERI_S_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M2_vrr;
      Double I_ERI_D2y_S_G2y2z_S_M1_vrr = PAY*I_ERI_Py_S_G2y2z_S_M1_vrr+WPY*I_ERI_Py_S_G2y2z_S_M2_vrr+oned2z*I_ERI_S_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_Py_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_Pz_S_G2y2z_S_M1_vrr+WPZ*I_ERI_Pz_S_G2y2z_S_M2_vrr+oned2z*I_ERI_S_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M2_vrr;
      Double I_ERI_D2x_S_Gy3z_S_M1_vrr = PAX*I_ERI_Px_S_Gy3z_S_M1_vrr+WPX*I_ERI_Px_S_Gy3z_S_M2_vrr+oned2z*I_ERI_S_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M2_vrr;
      Double I_ERI_D2y_S_Gy3z_S_M1_vrr = PAY*I_ERI_Py_S_Gy3z_S_M1_vrr+WPY*I_ERI_Py_S_Gy3z_S_M2_vrr+oned2z*I_ERI_S_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M2_vrr+oned2k*I_ERI_Py_S_F3z_S_M2_vrr;
      Double I_ERI_D2z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_Pz_S_Gy3z_S_M1_vrr+WPZ*I_ERI_Pz_S_Gy3z_S_M2_vrr+oned2z*I_ERI_S_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2x_S_G4z_S_M1_vrr = PAX*I_ERI_Px_S_G4z_S_M1_vrr+WPX*I_ERI_Px_S_G4z_S_M2_vrr+oned2z*I_ERI_S_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M2_vrr;
      Double I_ERI_D2y_S_G4z_S_M1_vrr = PAY*I_ERI_Py_S_G4z_S_M1_vrr+WPY*I_ERI_Py_S_G4z_S_M2_vrr+oned2z*I_ERI_S_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M2_vrr;
      Double I_ERI_D2z_S_G4z_S_M1_vrr = PAZ*I_ERI_Pz_S_G4z_S_M1_vrr+WPZ*I_ERI_Pz_S_G4z_S_M2_vrr+oned2z*I_ERI_S_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M2_vrr+4*oned2k*I_ERI_Pz_S_F3z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_G4x_S_M1_vrr = PAX*I_ERI_D2x_S_G4x_S_M1_vrr+WPX*I_ERI_D2x_S_G4x_S_M2_vrr+2*oned2z*I_ERI_Px_S_G4x_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G4x_S_M2_vrr+4*oned2k*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_F2xy_S_G4x_S_M1_vrr = PAY*I_ERI_D2x_S_G4x_S_M1_vrr+WPY*I_ERI_D2x_S_G4x_S_M2_vrr;
      Double I_ERI_F2xz_S_G4x_S_M1_vrr = PAZ*I_ERI_D2x_S_G4x_S_M1_vrr+WPZ*I_ERI_D2x_S_G4x_S_M2_vrr;
      Double I_ERI_Fx2y_S_G4x_S_M1_vrr = PAX*I_ERI_D2y_S_G4x_S_M1_vrr+WPX*I_ERI_D2y_S_G4x_S_M2_vrr+4*oned2k*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_Fx2z_S_G4x_S_M1_vrr = PAX*I_ERI_D2z_S_G4x_S_M1_vrr+WPX*I_ERI_D2z_S_G4x_S_M2_vrr+4*oned2k*I_ERI_D2z_S_F3x_S_M2_vrr;
      Double I_ERI_F3y_S_G4x_S_M1_vrr = PAY*I_ERI_D2y_S_G4x_S_M1_vrr+WPY*I_ERI_D2y_S_G4x_S_M2_vrr+2*oned2z*I_ERI_Py_S_G4x_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G4x_S_M2_vrr;
      Double I_ERI_F2yz_S_G4x_S_M1_vrr = PAZ*I_ERI_D2y_S_G4x_S_M1_vrr+WPZ*I_ERI_D2y_S_G4x_S_M2_vrr;
      Double I_ERI_F3z_S_G4x_S_M1_vrr = PAZ*I_ERI_D2z_S_G4x_S_M1_vrr+WPZ*I_ERI_D2z_S_G4x_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G4x_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G4x_S_M2_vrr;
      Double I_ERI_F3x_S_G3xy_S_M1_vrr = PAX*I_ERI_D2x_S_G3xy_S_M1_vrr+WPX*I_ERI_D2x_S_G3xy_S_M2_vrr+2*oned2z*I_ERI_Px_S_G3xy_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_F2xy_S_G3xy_S_M1_vrr = PAY*I_ERI_D2x_S_G3xy_S_M1_vrr+WPY*I_ERI_D2x_S_G3xy_S_M2_vrr+oned2k*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_F2xz_S_G3xy_S_M1_vrr = PAZ*I_ERI_D2x_S_G3xy_S_M1_vrr+WPZ*I_ERI_D2x_S_G3xy_S_M2_vrr;
      Double I_ERI_Fx2y_S_G3xy_S_M1_vrr = PAX*I_ERI_D2y_S_G3xy_S_M1_vrr+WPX*I_ERI_D2y_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_Fx2z_S_G3xy_S_M1_vrr = PAX*I_ERI_D2z_S_G3xy_S_M1_vrr+WPX*I_ERI_D2z_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_D2z_S_F2xy_S_M2_vrr;
      Double I_ERI_F3y_S_G3xy_S_M1_vrr = PAY*I_ERI_D2y_S_G3xy_S_M1_vrr+WPY*I_ERI_D2y_S_G3xy_S_M2_vrr+2*oned2z*I_ERI_Py_S_G3xy_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G3xy_S_M2_vrr+oned2k*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_F2yz_S_G3xy_S_M1_vrr = PAZ*I_ERI_D2y_S_G3xy_S_M1_vrr+WPZ*I_ERI_D2y_S_G3xy_S_M2_vrr;
      Double I_ERI_F3z_S_G3xy_S_M1_vrr = PAZ*I_ERI_D2z_S_G3xy_S_M1_vrr+WPZ*I_ERI_D2z_S_G3xy_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G3xy_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xy_S_M2_vrr;
      Double I_ERI_F3x_S_G3xz_S_M1_vrr = PAX*I_ERI_D2x_S_G3xz_S_M1_vrr+WPX*I_ERI_D2x_S_G3xz_S_M2_vrr+2*oned2z*I_ERI_Px_S_G3xz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_F2xy_S_G3xz_S_M1_vrr = PAY*I_ERI_D2x_S_G3xz_S_M1_vrr+WPY*I_ERI_D2x_S_G3xz_S_M2_vrr;
      Double I_ERI_F2xz_S_G3xz_S_M1_vrr = PAZ*I_ERI_D2x_S_G3xz_S_M1_vrr+WPZ*I_ERI_D2x_S_G3xz_S_M2_vrr+oned2k*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_Fx2y_S_G3xz_S_M1_vrr = PAX*I_ERI_D2y_S_G3xz_S_M1_vrr+WPX*I_ERI_D2y_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_D2y_S_F2xz_S_M2_vrr;
      Double I_ERI_Fx2z_S_G3xz_S_M1_vrr = PAX*I_ERI_D2z_S_G3xz_S_M1_vrr+WPX*I_ERI_D2z_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_D2z_S_F2xz_S_M2_vrr;
      Double I_ERI_F3y_S_G3xz_S_M1_vrr = PAY*I_ERI_D2y_S_G3xz_S_M1_vrr+WPY*I_ERI_D2y_S_G3xz_S_M2_vrr+2*oned2z*I_ERI_Py_S_G3xz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G3xz_S_M2_vrr;
      Double I_ERI_F2yz_S_G3xz_S_M1_vrr = PAZ*I_ERI_D2y_S_G3xz_S_M1_vrr+WPZ*I_ERI_D2y_S_G3xz_S_M2_vrr+oned2k*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_F3z_S_G3xz_S_M1_vrr = PAZ*I_ERI_D2z_S_G3xz_S_M1_vrr+WPZ*I_ERI_D2z_S_G3xz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G3xz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xz_S_M2_vrr+oned2k*I_ERI_D2z_S_F3x_S_M2_vrr;
      Double I_ERI_F3x_S_G2x2y_S_M1_vrr = PAX*I_ERI_D2x_S_G2x2y_S_M1_vrr+WPX*I_ERI_D2x_S_G2x2y_S_M2_vrr+2*oned2z*I_ERI_Px_S_G2x2y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_F2xy_S_G2x2y_S_M1_vrr = PAY*I_ERI_D2x_S_G2x2y_S_M1_vrr+WPY*I_ERI_D2x_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_F2xz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_D2x_S_G2x2y_S_M1_vrr+WPZ*I_ERI_D2x_S_G2x2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_G2x2y_S_M1_vrr = PAX*I_ERI_D2y_S_G2x2y_S_M1_vrr+WPX*I_ERI_D2y_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_G2x2y_S_M1_vrr = PAX*I_ERI_D2z_S_G2x2y_S_M1_vrr+WPX*I_ERI_D2z_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3y_S_G2x2y_S_M1_vrr = PAY*I_ERI_D2y_S_G2x2y_S_M1_vrr+WPY*I_ERI_D2y_S_G2x2y_S_M2_vrr+2*oned2z*I_ERI_Py_S_G2x2y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_F2yz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_D2y_S_G2x2y_S_M1_vrr+WPZ*I_ERI_D2y_S_G2x2y_S_M2_vrr;
      Double I_ERI_F3z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_D2z_S_G2x2y_S_M1_vrr+WPZ*I_ERI_D2z_S_G2x2y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G2x2y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2y_S_M2_vrr;
      Double I_ERI_F3x_S_G2xyz_S_M1_vrr = PAX*I_ERI_D2x_S_G2xyz_S_M1_vrr+WPX*I_ERI_D2x_S_G2xyz_S_M2_vrr+2*oned2z*I_ERI_Px_S_G2xyz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M2_vrr;
      Double I_ERI_F2xy_S_G2xyz_S_M1_vrr = PAY*I_ERI_D2x_S_G2xyz_S_M1_vrr+WPY*I_ERI_D2x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_F2xz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_D2x_S_G2xyz_S_M1_vrr+WPZ*I_ERI_D2x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_Fx2y_S_G2xyz_S_M1_vrr = PAX*I_ERI_D2y_S_G2xyz_S_M1_vrr+WPX*I_ERI_D2y_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M2_vrr;
      Double I_ERI_Fx2z_S_G2xyz_S_M1_vrr = PAX*I_ERI_D2z_S_G2xyz_S_M1_vrr+WPX*I_ERI_D2z_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M2_vrr;
      Double I_ERI_F3y_S_G2xyz_S_M1_vrr = PAY*I_ERI_D2y_S_G2xyz_S_M1_vrr+WPY*I_ERI_D2y_S_G2xyz_S_M2_vrr+2*oned2z*I_ERI_Py_S_G2xyz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2y_S_F2xz_S_M2_vrr;
      Double I_ERI_F2yz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_D2y_S_G2xyz_S_M1_vrr+WPZ*I_ERI_D2y_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_F3z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_D2z_S_G2xyz_S_M1_vrr+WPZ*I_ERI_D2z_S_G2xyz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G2xyz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G2xyz_S_M2_vrr+oned2k*I_ERI_D2z_S_F2xy_S_M2_vrr;
      Double I_ERI_F3x_S_G2x2z_S_M1_vrr = PAX*I_ERI_D2x_S_G2x2z_S_M1_vrr+WPX*I_ERI_D2x_S_G2x2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_G2x2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2xy_S_G2x2z_S_M1_vrr = PAY*I_ERI_D2x_S_G2x2z_S_M1_vrr+WPY*I_ERI_D2x_S_G2x2z_S_M2_vrr;
      Double I_ERI_F2xz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_D2x_S_G2x2z_S_M1_vrr+WPZ*I_ERI_D2x_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_Fx2y_S_G2x2z_S_M1_vrr = PAX*I_ERI_D2y_S_G2x2z_S_M1_vrr+WPX*I_ERI_D2y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fx2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_G2x2z_S_M1_vrr = PAX*I_ERI_D2z_S_G2x2z_S_M1_vrr+WPX*I_ERI_D2z_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Fx2z_S_M2_vrr;
      Double I_ERI_F3y_S_G2x2z_S_M1_vrr = PAY*I_ERI_D2y_S_G2x2z_S_M1_vrr+WPY*I_ERI_D2y_S_G2x2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_G2x2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2z_S_M2_vrr;
      Double I_ERI_F2yz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_D2y_S_G2x2z_S_M1_vrr+WPZ*I_ERI_D2y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_F2xz_S_M2_vrr;
      Double I_ERI_F3z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_D2z_S_G2x2z_S_M1_vrr+WPZ*I_ERI_D2z_S_G2x2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G2x2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_F2xz_S_M2_vrr;
      Double I_ERI_F3x_S_Gx3y_S_M1_vrr = PAX*I_ERI_D2x_S_Gx3y_S_M1_vrr+WPX*I_ERI_D2x_S_Gx3y_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gx3y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3y_S_M2_vrr+oned2k*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_F2xy_S_Gx3y_S_M1_vrr = PAY*I_ERI_D2x_S_Gx3y_S_M1_vrr+WPY*I_ERI_D2x_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_F2xz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_D2x_S_Gx3y_S_M1_vrr+WPZ*I_ERI_D2x_S_Gx3y_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gx3y_S_M1_vrr = PAX*I_ERI_D2y_S_Gx3y_S_M1_vrr+WPX*I_ERI_D2y_S_Gx3y_S_M2_vrr+oned2k*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gx3y_S_M1_vrr = PAX*I_ERI_D2z_S_Gx3y_S_M1_vrr+WPX*I_ERI_D2z_S_Gx3y_S_M2_vrr+oned2k*I_ERI_D2z_S_F3y_S_M2_vrr;
      Double I_ERI_F3y_S_Gx3y_S_M1_vrr = PAY*I_ERI_D2y_S_Gx3y_S_M1_vrr+WPY*I_ERI_D2y_S_Gx3y_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gx3y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_F2yz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_D2y_S_Gx3y_S_M1_vrr+WPZ*I_ERI_D2y_S_Gx3y_S_M2_vrr;
      Double I_ERI_F3z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_D2z_S_Gx3y_S_M1_vrr+WPZ*I_ERI_D2z_S_Gx3y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gx3y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3y_S_M2_vrr;
      Double I_ERI_F3x_S_Gx2yz_S_M1_vrr = PAX*I_ERI_D2x_S_Gx2yz_S_M1_vrr+WPX*I_ERI_D2x_S_Gx2yz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gx2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2x_S_F2yz_S_M2_vrr;
      Double I_ERI_F2xy_S_Gx2yz_S_M1_vrr = PAY*I_ERI_D2x_S_Gx2yz_S_M1_vrr+WPY*I_ERI_D2x_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M2_vrr;
      Double I_ERI_F2xz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_D2x_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_D2x_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gx2yz_S_M1_vrr = PAX*I_ERI_D2y_S_Gx2yz_S_M1_vrr+WPX*I_ERI_D2y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2y_S_F2yz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gx2yz_S_M1_vrr = PAX*I_ERI_D2z_S_Gx2yz_S_M1_vrr+WPX*I_ERI_D2z_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2z_S_F2yz_S_M2_vrr;
      Double I_ERI_F3y_S_Gx2yz_S_M1_vrr = PAY*I_ERI_D2y_S_Gx2yz_S_M1_vrr+WPY*I_ERI_D2y_S_Gx2yz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gx2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M2_vrr;
      Double I_ERI_F2yz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_D2y_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_D2y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_D2z_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_D2z_S_Gx2yz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gx2yz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_D2z_S_Fx2y_S_M2_vrr;
      Double I_ERI_F3x_S_Gxy2z_S_M1_vrr = PAX*I_ERI_D2x_S_Gxy2z_S_M1_vrr+WPX*I_ERI_D2x_S_Gxy2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gxy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2x_S_Fy2z_S_M2_vrr;
      Double I_ERI_F2xy_S_Gxy2z_S_M1_vrr = PAY*I_ERI_D2x_S_Gxy2z_S_M1_vrr+WPY*I_ERI_D2x_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2xz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_D2x_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_D2x_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gxy2z_S_M1_vrr = PAX*I_ERI_D2y_S_Gxy2z_S_M1_vrr+WPX*I_ERI_D2y_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gxy2z_S_M1_vrr = PAX*I_ERI_D2z_S_Gxy2z_S_M1_vrr+WPX*I_ERI_D2z_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2z_S_Fy2z_S_M2_vrr;
      Double I_ERI_F3y_S_Gxy2z_S_M1_vrr = PAY*I_ERI_D2y_S_Gxy2z_S_M1_vrr+WPY*I_ERI_D2y_S_Gxy2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gxy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_D2y_S_Fx2z_S_M2_vrr;
      Double I_ERI_F2yz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_D2y_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_D2y_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M2_vrr;
      Double I_ERI_F3z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_D2z_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_D2z_S_Gxy2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gxy2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M2_vrr;
      Double I_ERI_F3x_S_Gx3z_S_M1_vrr = PAX*I_ERI_D2x_S_Gx3z_S_M1_vrr+WPX*I_ERI_D2x_S_Gx3z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gx3z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3z_S_M2_vrr+oned2k*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_F2xy_S_Gx3z_S_M1_vrr = PAY*I_ERI_D2x_S_Gx3z_S_M1_vrr+WPY*I_ERI_D2x_S_Gx3z_S_M2_vrr;
      Double I_ERI_F2xz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_D2x_S_Gx3z_S_M1_vrr+WPZ*I_ERI_D2x_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gx3z_S_M1_vrr = PAX*I_ERI_D2y_S_Gx3z_S_M1_vrr+WPX*I_ERI_D2y_S_Gx3z_S_M2_vrr+oned2k*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gx3z_S_M1_vrr = PAX*I_ERI_D2z_S_Gx3z_S_M1_vrr+WPX*I_ERI_D2z_S_Gx3z_S_M2_vrr+oned2k*I_ERI_D2z_S_F3z_S_M2_vrr;
      Double I_ERI_F3y_S_Gx3z_S_M1_vrr = PAY*I_ERI_D2y_S_Gx3z_S_M1_vrr+WPY*I_ERI_D2y_S_Gx3z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gx3z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3z_S_M2_vrr;
      Double I_ERI_F2yz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_D2y_S_Gx3z_S_M1_vrr+WPZ*I_ERI_D2y_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_D2y_S_Fx2z_S_M2_vrr;
      Double I_ERI_F3z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_D2z_S_Gx3z_S_M1_vrr+WPZ*I_ERI_D2z_S_Gx3z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gx3z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_D2z_S_Fx2z_S_M2_vrr;
      Double I_ERI_F3x_S_G4y_S_M1_vrr = PAX*I_ERI_D2x_S_G4y_S_M1_vrr+WPX*I_ERI_D2x_S_G4y_S_M2_vrr+2*oned2z*I_ERI_Px_S_G4y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G4y_S_M2_vrr;
      Double I_ERI_F2xy_S_G4y_S_M1_vrr = PAY*I_ERI_D2x_S_G4y_S_M1_vrr+WPY*I_ERI_D2x_S_G4y_S_M2_vrr+4*oned2k*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_F2xz_S_G4y_S_M1_vrr = PAZ*I_ERI_D2x_S_G4y_S_M1_vrr+WPZ*I_ERI_D2x_S_G4y_S_M2_vrr;
      Double I_ERI_Fx2y_S_G4y_S_M1_vrr = PAX*I_ERI_D2y_S_G4y_S_M1_vrr+WPX*I_ERI_D2y_S_G4y_S_M2_vrr;
      Double I_ERI_Fx2z_S_G4y_S_M1_vrr = PAX*I_ERI_D2z_S_G4y_S_M1_vrr+WPX*I_ERI_D2z_S_G4y_S_M2_vrr;
      Double I_ERI_F3y_S_G4y_S_M1_vrr = PAY*I_ERI_D2y_S_G4y_S_M1_vrr+WPY*I_ERI_D2y_S_G4y_S_M2_vrr+2*oned2z*I_ERI_Py_S_G4y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G4y_S_M2_vrr+4*oned2k*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_F2yz_S_G4y_S_M1_vrr = PAZ*I_ERI_D2y_S_G4y_S_M1_vrr+WPZ*I_ERI_D2y_S_G4y_S_M2_vrr;
      Double I_ERI_F3z_S_G4y_S_M1_vrr = PAZ*I_ERI_D2z_S_G4y_S_M1_vrr+WPZ*I_ERI_D2z_S_G4y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G4y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G4y_S_M2_vrr;
      Double I_ERI_F3x_S_G3yz_S_M1_vrr = PAX*I_ERI_D2x_S_G3yz_S_M1_vrr+WPX*I_ERI_D2x_S_G3yz_S_M2_vrr+2*oned2z*I_ERI_Px_S_G3yz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G3yz_S_M2_vrr;
      Double I_ERI_F2xy_S_G3yz_S_M1_vrr = PAY*I_ERI_D2x_S_G3yz_S_M1_vrr+WPY*I_ERI_D2x_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_D2x_S_F2yz_S_M2_vrr;
      Double I_ERI_F2xz_S_G3yz_S_M1_vrr = PAZ*I_ERI_D2x_S_G3yz_S_M1_vrr+WPZ*I_ERI_D2x_S_G3yz_S_M2_vrr+oned2k*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_Fx2y_S_G3yz_S_M1_vrr = PAX*I_ERI_D2y_S_G3yz_S_M1_vrr+WPX*I_ERI_D2y_S_G3yz_S_M2_vrr;
      Double I_ERI_Fx2z_S_G3yz_S_M1_vrr = PAX*I_ERI_D2z_S_G3yz_S_M1_vrr+WPX*I_ERI_D2z_S_G3yz_S_M2_vrr;
      Double I_ERI_F3y_S_G3yz_S_M1_vrr = PAY*I_ERI_D2y_S_G3yz_S_M1_vrr+WPY*I_ERI_D2y_S_G3yz_S_M2_vrr+2*oned2z*I_ERI_Py_S_G3yz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_D2y_S_F2yz_S_M2_vrr;
      Double I_ERI_F2yz_S_G3yz_S_M1_vrr = PAZ*I_ERI_D2y_S_G3yz_S_M1_vrr+WPZ*I_ERI_D2y_S_G3yz_S_M2_vrr+oned2k*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_F3z_S_G3yz_S_M1_vrr = PAZ*I_ERI_D2z_S_G3yz_S_M1_vrr+WPZ*I_ERI_D2z_S_G3yz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G3yz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G3yz_S_M2_vrr+oned2k*I_ERI_D2z_S_F3y_S_M2_vrr;
      Double I_ERI_F3x_S_G2y2z_S_M1_vrr = PAX*I_ERI_D2x_S_G2y2z_S_M1_vrr+WPX*I_ERI_D2x_S_G2y2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_G2y2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G2y2z_S_M2_vrr;
      Double I_ERI_F2xy_S_G2y2z_S_M1_vrr = PAY*I_ERI_D2x_S_G2y2z_S_M1_vrr+WPY*I_ERI_D2x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Fy2z_S_M2_vrr;
      Double I_ERI_F2xz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_D2x_S_G2y2z_S_M1_vrr+WPZ*I_ERI_D2x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_F2yz_S_M2_vrr;
      Double I_ERI_Fx2y_S_G2y2z_S_M1_vrr = PAX*I_ERI_D2y_S_G2y2z_S_M1_vrr+WPX*I_ERI_D2y_S_G2y2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_G2y2z_S_M1_vrr = PAX*I_ERI_D2z_S_G2y2z_S_M1_vrr+WPX*I_ERI_D2z_S_G2y2z_S_M2_vrr;
      Double I_ERI_F3y_S_G2y2z_S_M1_vrr = PAY*I_ERI_D2y_S_G2y2z_S_M1_vrr+WPY*I_ERI_D2y_S_G2y2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_G2y2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_F2yz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_D2y_S_G2y2z_S_M1_vrr+WPZ*I_ERI_D2y_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_F2yz_S_M2_vrr;
      Double I_ERI_F3z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_D2z_S_G2y2z_S_M1_vrr+WPZ*I_ERI_D2z_S_G2y2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G2y2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_F2yz_S_M2_vrr;
      Double I_ERI_F3x_S_Gy3z_S_M1_vrr = PAX*I_ERI_D2x_S_Gy3z_S_M1_vrr+WPX*I_ERI_D2x_S_Gy3z_S_M2_vrr+2*oned2z*I_ERI_Px_S_Gy3z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Gy3z_S_M2_vrr;
      Double I_ERI_F2xy_S_Gy3z_S_M1_vrr = PAY*I_ERI_D2x_S_Gy3z_S_M1_vrr+WPY*I_ERI_D2x_S_Gy3z_S_M2_vrr+oned2k*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_F2xz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_D2x_S_Gy3z_S_M1_vrr+WPZ*I_ERI_D2x_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_D2x_S_Fy2z_S_M2_vrr;
      Double I_ERI_Fx2y_S_Gy3z_S_M1_vrr = PAX*I_ERI_D2y_S_Gy3z_S_M1_vrr+WPX*I_ERI_D2y_S_Gy3z_S_M2_vrr;
      Double I_ERI_Fx2z_S_Gy3z_S_M1_vrr = PAX*I_ERI_D2z_S_Gy3z_S_M1_vrr+WPX*I_ERI_D2z_S_Gy3z_S_M2_vrr;
      Double I_ERI_F3y_S_Gy3z_S_M1_vrr = PAY*I_ERI_D2y_S_Gy3z_S_M1_vrr+WPY*I_ERI_D2y_S_Gy3z_S_M2_vrr+2*oned2z*I_ERI_Py_S_Gy3z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Gy3z_S_M2_vrr+oned2k*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_F2yz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_D2y_S_Gy3z_S_M1_vrr+WPZ*I_ERI_D2y_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_D2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_F3z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_D2z_S_Gy3z_S_M1_vrr+WPZ*I_ERI_D2z_S_Gy3z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Gy3z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_D2z_S_Fy2z_S_M2_vrr;
      Double I_ERI_F3x_S_G4z_S_M1_vrr = PAX*I_ERI_D2x_S_G4z_S_M1_vrr+WPX*I_ERI_D2x_S_G4z_S_M2_vrr+2*oned2z*I_ERI_Px_S_G4z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_G4z_S_M2_vrr;
      Double I_ERI_F2xy_S_G4z_S_M1_vrr = PAY*I_ERI_D2x_S_G4z_S_M1_vrr+WPY*I_ERI_D2x_S_G4z_S_M2_vrr;
      Double I_ERI_F2xz_S_G4z_S_M1_vrr = PAZ*I_ERI_D2x_S_G4z_S_M1_vrr+WPZ*I_ERI_D2x_S_G4z_S_M2_vrr+4*oned2k*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_Fx2y_S_G4z_S_M1_vrr = PAX*I_ERI_D2y_S_G4z_S_M1_vrr+WPX*I_ERI_D2y_S_G4z_S_M2_vrr;
      Double I_ERI_Fx2z_S_G4z_S_M1_vrr = PAX*I_ERI_D2z_S_G4z_S_M1_vrr+WPX*I_ERI_D2z_S_G4z_S_M2_vrr;
      Double I_ERI_F3y_S_G4z_S_M1_vrr = PAY*I_ERI_D2y_S_G4z_S_M1_vrr+WPY*I_ERI_D2y_S_G4z_S_M2_vrr+2*oned2z*I_ERI_Py_S_G4z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_G4z_S_M2_vrr;
      Double I_ERI_F2yz_S_G4z_S_M1_vrr = PAZ*I_ERI_D2y_S_G4z_S_M1_vrr+WPZ*I_ERI_D2y_S_G4z_S_M2_vrr+4*oned2k*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_F3z_S_G4z_S_M1_vrr = PAZ*I_ERI_D2z_S_G4z_S_M1_vrr+WPZ*I_ERI_D2z_S_G4z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_G4z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_G4z_S_M2_vrr+4*oned2k*I_ERI_D2z_S_F3z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_F3x_S_M1_vrr = PAX*I_ERI_F3x_S_F3x_S_M1_vrr+WPX*I_ERI_F3x_S_F3x_S_M2_vrr+3*oned2z*I_ERI_D2x_S_F3x_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_F3x_S_M2_vrr+3*oned2k*I_ERI_F3x_S_D2x_S_M2_vrr;
      Double I_ERI_G3xy_S_F3x_S_M1_vrr = PAY*I_ERI_F3x_S_F3x_S_M1_vrr+WPY*I_ERI_F3x_S_F3x_S_M2_vrr;
      Double I_ERI_G3xz_S_F3x_S_M1_vrr = PAZ*I_ERI_F3x_S_F3x_S_M1_vrr+WPZ*I_ERI_F3x_S_F3x_S_M2_vrr;
      Double I_ERI_G2x2y_S_F3x_S_M1_vrr = PAY*I_ERI_F2xy_S_F3x_S_M1_vrr+WPY*I_ERI_F2xy_S_F3x_S_M2_vrr+oned2z*I_ERI_D2x_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_G2x2z_S_F3x_S_M1_vrr = PAZ*I_ERI_F2xz_S_F3x_S_M1_vrr+WPZ*I_ERI_F2xz_S_F3x_S_M2_vrr+oned2z*I_ERI_D2x_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F3x_S_M2_vrr;
      Double I_ERI_Gx3y_S_F3x_S_M1_vrr = PAX*I_ERI_F3y_S_F3x_S_M1_vrr+WPX*I_ERI_F3y_S_F3x_S_M2_vrr+3*oned2k*I_ERI_F3y_S_D2x_S_M2_vrr;
      Double I_ERI_Gx3z_S_F3x_S_M1_vrr = PAX*I_ERI_F3z_S_F3x_S_M1_vrr+WPX*I_ERI_F3z_S_F3x_S_M2_vrr+3*oned2k*I_ERI_F3z_S_D2x_S_M2_vrr;
      Double I_ERI_G4y_S_F3x_S_M1_vrr = PAY*I_ERI_F3y_S_F3x_S_M1_vrr+WPY*I_ERI_F3y_S_F3x_S_M2_vrr+3*oned2z*I_ERI_D2y_S_F3x_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_G3yz_S_F3x_S_M1_vrr = PAZ*I_ERI_F3y_S_F3x_S_M1_vrr+WPZ*I_ERI_F3y_S_F3x_S_M2_vrr;
      Double I_ERI_G2y2z_S_F3x_S_M1_vrr = PAZ*I_ERI_F2yz_S_F3x_S_M1_vrr+WPZ*I_ERI_F2yz_S_F3x_S_M2_vrr+oned2z*I_ERI_D2y_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_F3x_S_M2_vrr;
      Double I_ERI_Gy3z_S_F3x_S_M1_vrr = PAY*I_ERI_F3z_S_F3x_S_M1_vrr+WPY*I_ERI_F3z_S_F3x_S_M2_vrr;
      Double I_ERI_G4z_S_F3x_S_M1_vrr = PAZ*I_ERI_F3z_S_F3x_S_M1_vrr+WPZ*I_ERI_F3z_S_F3x_S_M2_vrr+3*oned2z*I_ERI_D2z_S_F3x_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_F3x_S_M2_vrr;
      Double I_ERI_G4x_S_F2xy_S_M1_vrr = PAX*I_ERI_F3x_S_F2xy_S_M1_vrr+WPX*I_ERI_F3x_S_F2xy_S_M2_vrr+3*oned2z*I_ERI_D2x_S_F2xy_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Dxy_S_M2_vrr;
      Double I_ERI_G3xy_S_F2xy_S_M1_vrr = PAY*I_ERI_F3x_S_F2xy_S_M1_vrr+WPY*I_ERI_F3x_S_F2xy_S_M2_vrr+oned2k*I_ERI_F3x_S_D2x_S_M2_vrr;
      Double I_ERI_G3xz_S_F2xy_S_M1_vrr = PAZ*I_ERI_F3x_S_F2xy_S_M1_vrr+WPZ*I_ERI_F3x_S_F2xy_S_M2_vrr;
      Double I_ERI_G2x2y_S_F2xy_S_M1_vrr = PAY*I_ERI_F2xy_S_F2xy_S_M1_vrr+WPY*I_ERI_F2xy_S_F2xy_S_M2_vrr+oned2z*I_ERI_D2x_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F2xy_S_M2_vrr+oned2k*I_ERI_F2xy_S_D2x_S_M2_vrr;
      Double I_ERI_G2x2z_S_F2xy_S_M1_vrr = PAZ*I_ERI_F2xz_S_F2xy_S_M1_vrr+WPZ*I_ERI_F2xz_S_F2xy_S_M2_vrr+oned2z*I_ERI_D2x_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F2xy_S_M2_vrr;
      Double I_ERI_Gx3y_S_F2xy_S_M1_vrr = PAX*I_ERI_F3y_S_F2xy_S_M1_vrr+WPX*I_ERI_F3y_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Dxy_S_M2_vrr;
      Double I_ERI_Gx3z_S_F2xy_S_M1_vrr = PAX*I_ERI_F3z_S_F2xy_S_M1_vrr+WPX*I_ERI_F3z_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Dxy_S_M2_vrr;
      Double I_ERI_G4y_S_F2xy_S_M1_vrr = PAY*I_ERI_F3y_S_F2xy_S_M1_vrr+WPY*I_ERI_F3y_S_F2xy_S_M2_vrr+3*oned2z*I_ERI_D2y_S_F2xy_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_F2xy_S_M2_vrr+oned2k*I_ERI_F3y_S_D2x_S_M2_vrr;
      Double I_ERI_G3yz_S_F2xy_S_M1_vrr = PAZ*I_ERI_F3y_S_F2xy_S_M1_vrr+WPZ*I_ERI_F3y_S_F2xy_S_M2_vrr;
      Double I_ERI_G2y2z_S_F2xy_S_M1_vrr = PAZ*I_ERI_F2yz_S_F2xy_S_M1_vrr+WPZ*I_ERI_F2yz_S_F2xy_S_M2_vrr+oned2z*I_ERI_D2y_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_F2xy_S_M2_vrr;
      Double I_ERI_Gy3z_S_F2xy_S_M1_vrr = PAY*I_ERI_F3z_S_F2xy_S_M1_vrr+WPY*I_ERI_F3z_S_F2xy_S_M2_vrr+oned2k*I_ERI_F3z_S_D2x_S_M2_vrr;
      Double I_ERI_G4z_S_F2xy_S_M1_vrr = PAZ*I_ERI_F3z_S_F2xy_S_M1_vrr+WPZ*I_ERI_F3z_S_F2xy_S_M2_vrr+3*oned2z*I_ERI_D2z_S_F2xy_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_F2xy_S_M2_vrr;
      Double I_ERI_G4x_S_F2xz_S_M1_vrr = PAX*I_ERI_F3x_S_F2xz_S_M1_vrr+WPX*I_ERI_F3x_S_F2xz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_F2xz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Dxz_S_M2_vrr;
      Double I_ERI_G3xy_S_F2xz_S_M1_vrr = PAY*I_ERI_F3x_S_F2xz_S_M1_vrr+WPY*I_ERI_F3x_S_F2xz_S_M2_vrr;
      Double I_ERI_G3xz_S_F2xz_S_M1_vrr = PAZ*I_ERI_F3x_S_F2xz_S_M1_vrr+WPZ*I_ERI_F3x_S_F2xz_S_M2_vrr+oned2k*I_ERI_F3x_S_D2x_S_M2_vrr;
      Double I_ERI_G2x2y_S_F2xz_S_M1_vrr = PAY*I_ERI_F2xy_S_F2xz_S_M1_vrr+WPY*I_ERI_F2xy_S_F2xz_S_M2_vrr+oned2z*I_ERI_D2x_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F2xz_S_M2_vrr;
      Double I_ERI_G2x2z_S_F2xz_S_M1_vrr = PAZ*I_ERI_F2xz_S_F2xz_S_M1_vrr+WPZ*I_ERI_F2xz_S_F2xz_S_M2_vrr+oned2z*I_ERI_D2x_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F2xz_S_M2_vrr+oned2k*I_ERI_F2xz_S_D2x_S_M2_vrr;
      Double I_ERI_Gx3y_S_F2xz_S_M1_vrr = PAX*I_ERI_F3y_S_F2xz_S_M1_vrr+WPX*I_ERI_F3y_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Dxz_S_M2_vrr;
      Double I_ERI_Gx3z_S_F2xz_S_M1_vrr = PAX*I_ERI_F3z_S_F2xz_S_M1_vrr+WPX*I_ERI_F3z_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Dxz_S_M2_vrr;
      Double I_ERI_G4y_S_F2xz_S_M1_vrr = PAY*I_ERI_F3y_S_F2xz_S_M1_vrr+WPY*I_ERI_F3y_S_F2xz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_F2xz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_F2xz_S_M2_vrr;
      Double I_ERI_G3yz_S_F2xz_S_M1_vrr = PAZ*I_ERI_F3y_S_F2xz_S_M1_vrr+WPZ*I_ERI_F3y_S_F2xz_S_M2_vrr+oned2k*I_ERI_F3y_S_D2x_S_M2_vrr;
      Double I_ERI_G2y2z_S_F2xz_S_M1_vrr = PAZ*I_ERI_F2yz_S_F2xz_S_M1_vrr+WPZ*I_ERI_F2yz_S_F2xz_S_M2_vrr+oned2z*I_ERI_D2y_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_F2xz_S_M2_vrr+oned2k*I_ERI_F2yz_S_D2x_S_M2_vrr;
      Double I_ERI_Gy3z_S_F2xz_S_M1_vrr = PAY*I_ERI_F3z_S_F2xz_S_M1_vrr+WPY*I_ERI_F3z_S_F2xz_S_M2_vrr;
      Double I_ERI_G4z_S_F2xz_S_M1_vrr = PAZ*I_ERI_F3z_S_F2xz_S_M1_vrr+WPZ*I_ERI_F3z_S_F2xz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_F2xz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_F2xz_S_M2_vrr+oned2k*I_ERI_F3z_S_D2x_S_M2_vrr;
      Double I_ERI_G4x_S_Fx2y_S_M1_vrr = PAX*I_ERI_F3x_S_Fx2y_S_M1_vrr+WPX*I_ERI_F3x_S_Fx2y_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Fx2y_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Fx2y_S_M2_vrr+oned2k*I_ERI_F3x_S_D2y_S_M2_vrr;
      Double I_ERI_G3xy_S_Fx2y_S_M1_vrr = PAY*I_ERI_F3x_S_Fx2y_S_M1_vrr+WPY*I_ERI_F3x_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Dxy_S_M2_vrr;
      Double I_ERI_G3xz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_F3x_S_Fx2y_S_M1_vrr+WPZ*I_ERI_F3x_S_Fx2y_S_M2_vrr;
      Double I_ERI_G2x2y_S_Fx2y_S_M1_vrr = PAY*I_ERI_F2xy_S_Fx2y_S_M1_vrr+WPY*I_ERI_F2xy_S_Fx2y_S_M2_vrr+oned2z*I_ERI_D2x_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Dxy_S_M2_vrr;
      Double I_ERI_G2x2z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_F2xz_S_Fx2y_S_M1_vrr+WPZ*I_ERI_F2xz_S_Fx2y_S_M2_vrr+oned2z*I_ERI_D2x_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Fx2y_S_M2_vrr;
      Double I_ERI_Gx3y_S_Fx2y_S_M1_vrr = PAX*I_ERI_F3y_S_Fx2y_S_M1_vrr+WPX*I_ERI_F3y_S_Fx2y_S_M2_vrr+oned2k*I_ERI_F3y_S_D2y_S_M2_vrr;
      Double I_ERI_Gx3z_S_Fx2y_S_M1_vrr = PAX*I_ERI_F3z_S_Fx2y_S_M1_vrr+WPX*I_ERI_F3z_S_Fx2y_S_M2_vrr+oned2k*I_ERI_F3z_S_D2y_S_M2_vrr;
      Double I_ERI_G4y_S_Fx2y_S_M1_vrr = PAY*I_ERI_F3y_S_Fx2y_S_M1_vrr+WPY*I_ERI_F3y_S_Fx2y_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Fx2y_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Dxy_S_M2_vrr;
      Double I_ERI_G3yz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_F3y_S_Fx2y_S_M1_vrr+WPZ*I_ERI_F3y_S_Fx2y_S_M2_vrr;
      Double I_ERI_G2y2z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_F2yz_S_Fx2y_S_M1_vrr+WPZ*I_ERI_F2yz_S_Fx2y_S_M2_vrr+oned2z*I_ERI_D2y_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_Gy3z_S_Fx2y_S_M1_vrr = PAY*I_ERI_F3z_S_Fx2y_S_M1_vrr+WPY*I_ERI_F3z_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Dxy_S_M2_vrr;
      Double I_ERI_G4z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_F3z_S_Fx2y_S_M1_vrr+WPZ*I_ERI_F3z_S_Fx2y_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Fx2y_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Fx2y_S_M2_vrr;
      Double I_ERI_G4x_S_Fxyz_S_M1_vrr = PAX*I_ERI_F3x_S_Fxyz_S_M1_vrr+WPX*I_ERI_F3x_S_Fxyz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Fxyz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F3x_S_Dyz_S_M2_vrr;
      Double I_ERI_G3xy_S_Fxyz_S_M1_vrr = PAY*I_ERI_F3x_S_Fxyz_S_M1_vrr+WPY*I_ERI_F3x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F3x_S_Dxz_S_M2_vrr;
      Double I_ERI_G3xz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_F3x_S_Fxyz_S_M1_vrr+WPZ*I_ERI_F3x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F3x_S_Dxy_S_M2_vrr;
      Double I_ERI_G2x2y_S_Fxyz_S_M1_vrr = PAY*I_ERI_F2xy_S_Fxyz_S_M1_vrr+WPY*I_ERI_F2xy_S_Fxyz_S_M2_vrr+oned2z*I_ERI_D2x_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F2xy_S_Dxz_S_M2_vrr;
      Double I_ERI_G2x2z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_F2xz_S_Fxyz_S_M1_vrr+WPZ*I_ERI_F2xz_S_Fxyz_S_M2_vrr+oned2z*I_ERI_D2x_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F2xz_S_Dxy_S_M2_vrr;
      Double I_ERI_Gx3y_S_Fxyz_S_M1_vrr = PAX*I_ERI_F3y_S_Fxyz_S_M1_vrr+WPX*I_ERI_F3y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F3y_S_Dyz_S_M2_vrr;
      Double I_ERI_Gx3z_S_Fxyz_S_M1_vrr = PAX*I_ERI_F3z_S_Fxyz_S_M1_vrr+WPX*I_ERI_F3z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F3z_S_Dyz_S_M2_vrr;
      Double I_ERI_G4y_S_Fxyz_S_M1_vrr = PAY*I_ERI_F3y_S_Fxyz_S_M1_vrr+WPY*I_ERI_F3y_S_Fxyz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Fxyz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F3y_S_Dxz_S_M2_vrr;
      Double I_ERI_G3yz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_F3y_S_Fxyz_S_M1_vrr+WPZ*I_ERI_F3y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F3y_S_Dxy_S_M2_vrr;
      Double I_ERI_G2y2z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_F2yz_S_Fxyz_S_M1_vrr+WPZ*I_ERI_F2yz_S_Fxyz_S_M2_vrr+oned2z*I_ERI_D2y_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F2yz_S_Dxy_S_M2_vrr;
      Double I_ERI_Gy3z_S_Fxyz_S_M1_vrr = PAY*I_ERI_F3z_S_Fxyz_S_M1_vrr+WPY*I_ERI_F3z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F3z_S_Dxz_S_M2_vrr;
      Double I_ERI_G4z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_F3z_S_Fxyz_S_M1_vrr+WPZ*I_ERI_F3z_S_Fxyz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Fxyz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_F3z_S_Dxy_S_M2_vrr;
      Double I_ERI_G4x_S_Fx2z_S_M1_vrr = PAX*I_ERI_F3x_S_Fx2z_S_M1_vrr+WPX*I_ERI_F3x_S_Fx2z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Fx2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Fx2z_S_M2_vrr+oned2k*I_ERI_F3x_S_D2z_S_M2_vrr;
      Double I_ERI_G3xy_S_Fx2z_S_M1_vrr = PAY*I_ERI_F3x_S_Fx2z_S_M1_vrr+WPY*I_ERI_F3x_S_Fx2z_S_M2_vrr;
      Double I_ERI_G3xz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_F3x_S_Fx2z_S_M1_vrr+WPZ*I_ERI_F3x_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Dxz_S_M2_vrr;
      Double I_ERI_G2x2y_S_Fx2z_S_M1_vrr = PAY*I_ERI_F2xy_S_Fx2z_S_M1_vrr+WPY*I_ERI_F2xy_S_Fx2z_S_M2_vrr+oned2z*I_ERI_D2x_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Fx2z_S_M2_vrr;
      Double I_ERI_G2x2z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_F2xz_S_Fx2z_S_M1_vrr+WPZ*I_ERI_F2xz_S_Fx2z_S_M2_vrr+oned2z*I_ERI_D2x_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_Dxz_S_M2_vrr;
      Double I_ERI_Gx3y_S_Fx2z_S_M1_vrr = PAX*I_ERI_F3y_S_Fx2z_S_M1_vrr+WPX*I_ERI_F3y_S_Fx2z_S_M2_vrr+oned2k*I_ERI_F3y_S_D2z_S_M2_vrr;
      Double I_ERI_Gx3z_S_Fx2z_S_M1_vrr = PAX*I_ERI_F3z_S_Fx2z_S_M1_vrr+WPX*I_ERI_F3z_S_Fx2z_S_M2_vrr+oned2k*I_ERI_F3z_S_D2z_S_M2_vrr;
      Double I_ERI_G4y_S_Fx2z_S_M1_vrr = PAY*I_ERI_F3y_S_Fx2z_S_M1_vrr+WPY*I_ERI_F3y_S_Fx2z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Fx2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Fx2z_S_M2_vrr;
      Double I_ERI_G3yz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_F3y_S_Fx2z_S_M1_vrr+WPZ*I_ERI_F3y_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Dxz_S_M2_vrr;
      Double I_ERI_G2y2z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_F2yz_S_Fx2z_S_M1_vrr+WPZ*I_ERI_F2yz_S_Fx2z_S_M2_vrr+oned2z*I_ERI_D2y_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_Dxz_S_M2_vrr;
      Double I_ERI_Gy3z_S_Fx2z_S_M1_vrr = PAY*I_ERI_F3z_S_Fx2z_S_M1_vrr+WPY*I_ERI_F3z_S_Fx2z_S_M2_vrr;
      Double I_ERI_G4z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_F3z_S_Fx2z_S_M1_vrr+WPZ*I_ERI_F3z_S_Fx2z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Fx2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Dxz_S_M2_vrr;
      Double I_ERI_G4x_S_F3y_S_M1_vrr = PAX*I_ERI_F3x_S_F3y_S_M1_vrr+WPX*I_ERI_F3x_S_F3y_S_M2_vrr+3*oned2z*I_ERI_D2x_S_F3y_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_G3xy_S_F3y_S_M1_vrr = PAY*I_ERI_F3x_S_F3y_S_M1_vrr+WPY*I_ERI_F3x_S_F3y_S_M2_vrr+3*oned2k*I_ERI_F3x_S_D2y_S_M2_vrr;
      Double I_ERI_G3xz_S_F3y_S_M1_vrr = PAZ*I_ERI_F3x_S_F3y_S_M1_vrr+WPZ*I_ERI_F3x_S_F3y_S_M2_vrr;
      Double I_ERI_G2x2y_S_F3y_S_M1_vrr = PAY*I_ERI_F2xy_S_F3y_S_M1_vrr+WPY*I_ERI_F2xy_S_F3y_S_M2_vrr+oned2z*I_ERI_D2x_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F3y_S_M2_vrr+3*oned2k*I_ERI_F2xy_S_D2y_S_M2_vrr;
      Double I_ERI_G2x2z_S_F3y_S_M1_vrr = PAZ*I_ERI_F2xz_S_F3y_S_M1_vrr+WPZ*I_ERI_F2xz_S_F3y_S_M2_vrr+oned2z*I_ERI_D2x_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F3y_S_M2_vrr;
      Double I_ERI_Gx3y_S_F3y_S_M1_vrr = PAX*I_ERI_F3y_S_F3y_S_M1_vrr+WPX*I_ERI_F3y_S_F3y_S_M2_vrr;
      Double I_ERI_Gx3z_S_F3y_S_M1_vrr = PAX*I_ERI_F3z_S_F3y_S_M1_vrr+WPX*I_ERI_F3z_S_F3y_S_M2_vrr;
      Double I_ERI_G4y_S_F3y_S_M1_vrr = PAY*I_ERI_F3y_S_F3y_S_M1_vrr+WPY*I_ERI_F3y_S_F3y_S_M2_vrr+3*oned2z*I_ERI_D2y_S_F3y_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_F3y_S_M2_vrr+3*oned2k*I_ERI_F3y_S_D2y_S_M2_vrr;
      Double I_ERI_G3yz_S_F3y_S_M1_vrr = PAZ*I_ERI_F3y_S_F3y_S_M1_vrr+WPZ*I_ERI_F3y_S_F3y_S_M2_vrr;
      Double I_ERI_G2y2z_S_F3y_S_M1_vrr = PAZ*I_ERI_F2yz_S_F3y_S_M1_vrr+WPZ*I_ERI_F2yz_S_F3y_S_M2_vrr+oned2z*I_ERI_D2y_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_F3y_S_M2_vrr;
      Double I_ERI_Gy3z_S_F3y_S_M1_vrr = PAY*I_ERI_F3z_S_F3y_S_M1_vrr+WPY*I_ERI_F3z_S_F3y_S_M2_vrr+3*oned2k*I_ERI_F3z_S_D2y_S_M2_vrr;
      Double I_ERI_G4z_S_F3y_S_M1_vrr = PAZ*I_ERI_F3z_S_F3y_S_M1_vrr+WPZ*I_ERI_F3z_S_F3y_S_M2_vrr+3*oned2z*I_ERI_D2z_S_F3y_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_F3y_S_M2_vrr;
      Double I_ERI_G4x_S_F2yz_S_M1_vrr = PAX*I_ERI_F3x_S_F2yz_S_M1_vrr+WPX*I_ERI_F3x_S_F2yz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_F2yz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_F2yz_S_M2_vrr;
      Double I_ERI_G3xy_S_F2yz_S_M1_vrr = PAY*I_ERI_F3x_S_F2yz_S_M1_vrr+WPY*I_ERI_F3x_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Dyz_S_M2_vrr;
      Double I_ERI_G3xz_S_F2yz_S_M1_vrr = PAZ*I_ERI_F3x_S_F2yz_S_M1_vrr+WPZ*I_ERI_F3x_S_F2yz_S_M2_vrr+oned2k*I_ERI_F3x_S_D2y_S_M2_vrr;
      Double I_ERI_G2x2y_S_F2yz_S_M1_vrr = PAY*I_ERI_F2xy_S_F2yz_S_M1_vrr+WPY*I_ERI_F2xy_S_F2yz_S_M2_vrr+oned2z*I_ERI_D2x_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Dyz_S_M2_vrr;
      Double I_ERI_G2x2z_S_F2yz_S_M1_vrr = PAZ*I_ERI_F2xz_S_F2yz_S_M1_vrr+WPZ*I_ERI_F2xz_S_F2yz_S_M2_vrr+oned2z*I_ERI_D2x_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F2yz_S_M2_vrr+oned2k*I_ERI_F2xz_S_D2y_S_M2_vrr;
      Double I_ERI_Gx3y_S_F2yz_S_M1_vrr = PAX*I_ERI_F3y_S_F2yz_S_M1_vrr+WPX*I_ERI_F3y_S_F2yz_S_M2_vrr;
      Double I_ERI_Gx3z_S_F2yz_S_M1_vrr = PAX*I_ERI_F3z_S_F2yz_S_M1_vrr+WPX*I_ERI_F3z_S_F2yz_S_M2_vrr;
      Double I_ERI_G4y_S_F2yz_S_M1_vrr = PAY*I_ERI_F3y_S_F2yz_S_M1_vrr+WPY*I_ERI_F3y_S_F2yz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_F2yz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Dyz_S_M2_vrr;
      Double I_ERI_G3yz_S_F2yz_S_M1_vrr = PAZ*I_ERI_F3y_S_F2yz_S_M1_vrr+WPZ*I_ERI_F3y_S_F2yz_S_M2_vrr+oned2k*I_ERI_F3y_S_D2y_S_M2_vrr;
      Double I_ERI_G2y2z_S_F2yz_S_M1_vrr = PAZ*I_ERI_F2yz_S_F2yz_S_M1_vrr+WPZ*I_ERI_F2yz_S_F2yz_S_M2_vrr+oned2z*I_ERI_D2y_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_F2yz_S_M2_vrr+oned2k*I_ERI_F2yz_S_D2y_S_M2_vrr;
      Double I_ERI_Gy3z_S_F2yz_S_M1_vrr = PAY*I_ERI_F3z_S_F2yz_S_M1_vrr+WPY*I_ERI_F3z_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Dyz_S_M2_vrr;
      Double I_ERI_G4z_S_F2yz_S_M1_vrr = PAZ*I_ERI_F3z_S_F2yz_S_M1_vrr+WPZ*I_ERI_F3z_S_F2yz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_F2yz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_F2yz_S_M2_vrr+oned2k*I_ERI_F3z_S_D2y_S_M2_vrr;
      Double I_ERI_G4x_S_Fy2z_S_M1_vrr = PAX*I_ERI_F3x_S_Fy2z_S_M1_vrr+WPX*I_ERI_F3x_S_Fy2z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Fy2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Fy2z_S_M2_vrr;
      Double I_ERI_G3xy_S_Fy2z_S_M1_vrr = PAY*I_ERI_F3x_S_Fy2z_S_M1_vrr+WPY*I_ERI_F3x_S_Fy2z_S_M2_vrr+oned2k*I_ERI_F3x_S_D2z_S_M2_vrr;
      Double I_ERI_G3xz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_F3x_S_Fy2z_S_M1_vrr+WPZ*I_ERI_F3x_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Dyz_S_M2_vrr;
      Double I_ERI_G2x2y_S_Fy2z_S_M1_vrr = PAY*I_ERI_F2xy_S_Fy2z_S_M1_vrr+WPY*I_ERI_F2xy_S_Fy2z_S_M2_vrr+oned2z*I_ERI_D2x_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Fy2z_S_M2_vrr+oned2k*I_ERI_F2xy_S_D2z_S_M2_vrr;
      Double I_ERI_G2x2z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_F2xz_S_Fy2z_S_M1_vrr+WPZ*I_ERI_F2xz_S_Fy2z_S_M2_vrr+oned2z*I_ERI_D2x_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_Dyz_S_M2_vrr;
      Double I_ERI_Gx3y_S_Fy2z_S_M1_vrr = PAX*I_ERI_F3y_S_Fy2z_S_M1_vrr+WPX*I_ERI_F3y_S_Fy2z_S_M2_vrr;
      Double I_ERI_Gx3z_S_Fy2z_S_M1_vrr = PAX*I_ERI_F3z_S_Fy2z_S_M1_vrr+WPX*I_ERI_F3z_S_Fy2z_S_M2_vrr;
      Double I_ERI_G4y_S_Fy2z_S_M1_vrr = PAY*I_ERI_F3y_S_Fy2z_S_M1_vrr+WPY*I_ERI_F3y_S_Fy2z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Fy2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Fy2z_S_M2_vrr+oned2k*I_ERI_F3y_S_D2z_S_M2_vrr;
      Double I_ERI_G3yz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_F3y_S_Fy2z_S_M1_vrr+WPZ*I_ERI_F3y_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Dyz_S_M2_vrr;
      Double I_ERI_G2y2z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_F2yz_S_Fy2z_S_M1_vrr+WPZ*I_ERI_F2yz_S_Fy2z_S_M2_vrr+oned2z*I_ERI_D2y_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_Dyz_S_M2_vrr;
      Double I_ERI_Gy3z_S_Fy2z_S_M1_vrr = PAY*I_ERI_F3z_S_Fy2z_S_M1_vrr+WPY*I_ERI_F3z_S_Fy2z_S_M2_vrr+oned2k*I_ERI_F3z_S_D2z_S_M2_vrr;
      Double I_ERI_G4z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_F3z_S_Fy2z_S_M1_vrr+WPZ*I_ERI_F3z_S_Fy2z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Fy2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Dyz_S_M2_vrr;
      Double I_ERI_G4x_S_F3z_S_M1_vrr = PAX*I_ERI_F3x_S_F3z_S_M1_vrr+WPX*I_ERI_F3x_S_F3z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_F3z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_G3xy_S_F3z_S_M1_vrr = PAY*I_ERI_F3x_S_F3z_S_M1_vrr+WPY*I_ERI_F3x_S_F3z_S_M2_vrr;
      Double I_ERI_G3xz_S_F3z_S_M1_vrr = PAZ*I_ERI_F3x_S_F3z_S_M1_vrr+WPZ*I_ERI_F3x_S_F3z_S_M2_vrr+3*oned2k*I_ERI_F3x_S_D2z_S_M2_vrr;
      Double I_ERI_G2x2y_S_F3z_S_M1_vrr = PAY*I_ERI_F2xy_S_F3z_S_M1_vrr+WPY*I_ERI_F2xy_S_F3z_S_M2_vrr+oned2z*I_ERI_D2x_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F3z_S_M2_vrr;
      Double I_ERI_G2x2z_S_F3z_S_M1_vrr = PAZ*I_ERI_F2xz_S_F3z_S_M1_vrr+WPZ*I_ERI_F2xz_S_F3z_S_M2_vrr+oned2z*I_ERI_D2x_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_F3z_S_M2_vrr+3*oned2k*I_ERI_F2xz_S_D2z_S_M2_vrr;
      Double I_ERI_Gx3y_S_F3z_S_M1_vrr = PAX*I_ERI_F3y_S_F3z_S_M1_vrr+WPX*I_ERI_F3y_S_F3z_S_M2_vrr;
      Double I_ERI_Gx3z_S_F3z_S_M1_vrr = PAX*I_ERI_F3z_S_F3z_S_M1_vrr+WPX*I_ERI_F3z_S_F3z_S_M2_vrr;
      Double I_ERI_G4y_S_F3z_S_M1_vrr = PAY*I_ERI_F3y_S_F3z_S_M1_vrr+WPY*I_ERI_F3y_S_F3z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_F3z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_F3z_S_M2_vrr;
      Double I_ERI_G3yz_S_F3z_S_M1_vrr = PAZ*I_ERI_F3y_S_F3z_S_M1_vrr+WPZ*I_ERI_F3y_S_F3z_S_M2_vrr+3*oned2k*I_ERI_F3y_S_D2z_S_M2_vrr;
      Double I_ERI_G2y2z_S_F3z_S_M1_vrr = PAZ*I_ERI_F2yz_S_F3z_S_M1_vrr+WPZ*I_ERI_F2yz_S_F3z_S_M2_vrr+oned2z*I_ERI_D2y_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_F3z_S_M2_vrr+3*oned2k*I_ERI_F2yz_S_D2z_S_M2_vrr;
      Double I_ERI_Gy3z_S_F3z_S_M1_vrr = PAY*I_ERI_F3z_S_F3z_S_M1_vrr+WPY*I_ERI_F3z_S_F3z_S_M2_vrr;
      Double I_ERI_G4z_S_F3z_S_M1_vrr = PAZ*I_ERI_F3z_S_F3z_S_M1_vrr+WPZ*I_ERI_F3z_S_F3z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_F3z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_F3z_S_M2_vrr+3*oned2k*I_ERI_F3z_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 45 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_G4x_S_M1_vrr = PAX*I_ERI_F3x_S_G4x_S_M1_vrr+WPX*I_ERI_F3x_S_G4x_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G4x_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G4x_S_M2_vrr+4*oned2k*I_ERI_F3x_S_F3x_S_M2_vrr;
      Double I_ERI_G3xy_S_G4x_S_M1_vrr = PAY*I_ERI_F3x_S_G4x_S_M1_vrr+WPY*I_ERI_F3x_S_G4x_S_M2_vrr;
      Double I_ERI_G3xz_S_G4x_S_M1_vrr = PAZ*I_ERI_F3x_S_G4x_S_M1_vrr+WPZ*I_ERI_F3x_S_G4x_S_M2_vrr;
      Double I_ERI_G2x2y_S_G4x_S_M1_vrr = PAY*I_ERI_F2xy_S_G4x_S_M1_vrr+WPY*I_ERI_F2xy_S_G4x_S_M2_vrr+oned2z*I_ERI_D2x_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G4x_S_M2_vrr;
      Double I_ERI_G2x2z_S_G4x_S_M1_vrr = PAZ*I_ERI_F2xz_S_G4x_S_M1_vrr+WPZ*I_ERI_F2xz_S_G4x_S_M2_vrr+oned2z*I_ERI_D2x_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G4x_S_M2_vrr;
      Double I_ERI_Gx3y_S_G4x_S_M1_vrr = PAX*I_ERI_F3y_S_G4x_S_M1_vrr+WPX*I_ERI_F3y_S_G4x_S_M2_vrr+4*oned2k*I_ERI_F3y_S_F3x_S_M2_vrr;
      Double I_ERI_Gx3z_S_G4x_S_M1_vrr = PAX*I_ERI_F3z_S_G4x_S_M1_vrr+WPX*I_ERI_F3z_S_G4x_S_M2_vrr+4*oned2k*I_ERI_F3z_S_F3x_S_M2_vrr;
      Double I_ERI_G4y_S_G4x_S_M1_vrr = PAY*I_ERI_F3y_S_G4x_S_M1_vrr+WPY*I_ERI_F3y_S_G4x_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G4x_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G4x_S_M2_vrr;
      Double I_ERI_G3yz_S_G4x_S_M1_vrr = PAZ*I_ERI_F3y_S_G4x_S_M1_vrr+WPZ*I_ERI_F3y_S_G4x_S_M2_vrr;
      Double I_ERI_G2y2z_S_G4x_S_M1_vrr = PAZ*I_ERI_F2yz_S_G4x_S_M1_vrr+WPZ*I_ERI_F2yz_S_G4x_S_M2_vrr+oned2z*I_ERI_D2y_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G4x_S_M2_vrr;
      Double I_ERI_Gy3z_S_G4x_S_M1_vrr = PAY*I_ERI_F3z_S_G4x_S_M1_vrr+WPY*I_ERI_F3z_S_G4x_S_M2_vrr;
      Double I_ERI_G4z_S_G4x_S_M1_vrr = PAZ*I_ERI_F3z_S_G4x_S_M1_vrr+WPZ*I_ERI_F3z_S_G4x_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G4x_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G4x_S_M2_vrr;
      Double I_ERI_G4x_S_G3xy_S_M1_vrr = PAX*I_ERI_F3x_S_G3xy_S_M1_vrr+WPX*I_ERI_F3x_S_G3xy_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G3xy_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_F3x_S_F2xy_S_M2_vrr;
      Double I_ERI_G3xy_S_G3xy_S_M1_vrr = PAY*I_ERI_F3x_S_G3xy_S_M1_vrr+WPY*I_ERI_F3x_S_G3xy_S_M2_vrr+oned2k*I_ERI_F3x_S_F3x_S_M2_vrr;
      Double I_ERI_G3xz_S_G3xy_S_M1_vrr = PAZ*I_ERI_F3x_S_G3xy_S_M1_vrr+WPZ*I_ERI_F3x_S_G3xy_S_M2_vrr;
      Double I_ERI_G2x2y_S_G3xy_S_M1_vrr = PAY*I_ERI_F2xy_S_G3xy_S_M1_vrr+WPY*I_ERI_F2xy_S_G3xy_S_M2_vrr+oned2z*I_ERI_D2x_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G3xy_S_M2_vrr+oned2k*I_ERI_F2xy_S_F3x_S_M2_vrr;
      Double I_ERI_G2x2z_S_G3xy_S_M1_vrr = PAZ*I_ERI_F2xz_S_G3xy_S_M1_vrr+WPZ*I_ERI_F2xz_S_G3xy_S_M2_vrr+oned2z*I_ERI_D2x_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G3xy_S_M2_vrr;
      Double I_ERI_Gx3y_S_G3xy_S_M1_vrr = PAX*I_ERI_F3y_S_G3xy_S_M1_vrr+WPX*I_ERI_F3y_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_F3y_S_F2xy_S_M2_vrr;
      Double I_ERI_Gx3z_S_G3xy_S_M1_vrr = PAX*I_ERI_F3z_S_G3xy_S_M1_vrr+WPX*I_ERI_F3z_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_F3z_S_F2xy_S_M2_vrr;
      Double I_ERI_G4y_S_G3xy_S_M1_vrr = PAY*I_ERI_F3y_S_G3xy_S_M1_vrr+WPY*I_ERI_F3y_S_G3xy_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G3xy_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G3xy_S_M2_vrr+oned2k*I_ERI_F3y_S_F3x_S_M2_vrr;
      Double I_ERI_G3yz_S_G3xy_S_M1_vrr = PAZ*I_ERI_F3y_S_G3xy_S_M1_vrr+WPZ*I_ERI_F3y_S_G3xy_S_M2_vrr;
      Double I_ERI_G2y2z_S_G3xy_S_M1_vrr = PAZ*I_ERI_F2yz_S_G3xy_S_M1_vrr+WPZ*I_ERI_F2yz_S_G3xy_S_M2_vrr+oned2z*I_ERI_D2y_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G3xy_S_M2_vrr;
      Double I_ERI_Gy3z_S_G3xy_S_M1_vrr = PAY*I_ERI_F3z_S_G3xy_S_M1_vrr+WPY*I_ERI_F3z_S_G3xy_S_M2_vrr+oned2k*I_ERI_F3z_S_F3x_S_M2_vrr;
      Double I_ERI_G4z_S_G3xy_S_M1_vrr = PAZ*I_ERI_F3z_S_G3xy_S_M1_vrr+WPZ*I_ERI_F3z_S_G3xy_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G3xy_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G3xy_S_M2_vrr;
      Double I_ERI_G4x_S_G3xz_S_M1_vrr = PAX*I_ERI_F3x_S_G3xz_S_M1_vrr+WPX*I_ERI_F3x_S_G3xz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G3xz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_F3x_S_F2xz_S_M2_vrr;
      Double I_ERI_G3xy_S_G3xz_S_M1_vrr = PAY*I_ERI_F3x_S_G3xz_S_M1_vrr+WPY*I_ERI_F3x_S_G3xz_S_M2_vrr;
      Double I_ERI_G3xz_S_G3xz_S_M1_vrr = PAZ*I_ERI_F3x_S_G3xz_S_M1_vrr+WPZ*I_ERI_F3x_S_G3xz_S_M2_vrr+oned2k*I_ERI_F3x_S_F3x_S_M2_vrr;
      Double I_ERI_G2x2y_S_G3xz_S_M1_vrr = PAY*I_ERI_F2xy_S_G3xz_S_M1_vrr+WPY*I_ERI_F2xy_S_G3xz_S_M2_vrr+oned2z*I_ERI_D2x_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G3xz_S_M2_vrr;
      Double I_ERI_G2x2z_S_G3xz_S_M1_vrr = PAZ*I_ERI_F2xz_S_G3xz_S_M1_vrr+WPZ*I_ERI_F2xz_S_G3xz_S_M2_vrr+oned2z*I_ERI_D2x_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G3xz_S_M2_vrr+oned2k*I_ERI_F2xz_S_F3x_S_M2_vrr;
      Double I_ERI_Gx3y_S_G3xz_S_M1_vrr = PAX*I_ERI_F3y_S_G3xz_S_M1_vrr+WPX*I_ERI_F3y_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_F3y_S_F2xz_S_M2_vrr;
      Double I_ERI_Gx3z_S_G3xz_S_M1_vrr = PAX*I_ERI_F3z_S_G3xz_S_M1_vrr+WPX*I_ERI_F3z_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_F3z_S_F2xz_S_M2_vrr;
      Double I_ERI_G4y_S_G3xz_S_M1_vrr = PAY*I_ERI_F3y_S_G3xz_S_M1_vrr+WPY*I_ERI_F3y_S_G3xz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G3xz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G3xz_S_M2_vrr;
      Double I_ERI_G3yz_S_G3xz_S_M1_vrr = PAZ*I_ERI_F3y_S_G3xz_S_M1_vrr+WPZ*I_ERI_F3y_S_G3xz_S_M2_vrr+oned2k*I_ERI_F3y_S_F3x_S_M2_vrr;
      Double I_ERI_G2y2z_S_G3xz_S_M1_vrr = PAZ*I_ERI_F2yz_S_G3xz_S_M1_vrr+WPZ*I_ERI_F2yz_S_G3xz_S_M2_vrr+oned2z*I_ERI_D2y_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G3xz_S_M2_vrr+oned2k*I_ERI_F2yz_S_F3x_S_M2_vrr;
      Double I_ERI_Gy3z_S_G3xz_S_M1_vrr = PAY*I_ERI_F3z_S_G3xz_S_M1_vrr+WPY*I_ERI_F3z_S_G3xz_S_M2_vrr;
      Double I_ERI_G4z_S_G3xz_S_M1_vrr = PAZ*I_ERI_F3z_S_G3xz_S_M1_vrr+WPZ*I_ERI_F3z_S_G3xz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G3xz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G3xz_S_M2_vrr+oned2k*I_ERI_F3z_S_F3x_S_M2_vrr;
      Double I_ERI_G4x_S_G2x2y_S_M1_vrr = PAX*I_ERI_F3x_S_G2x2y_S_M1_vrr+WPX*I_ERI_F3x_S_G2x2y_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G2x2y_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Fx2y_S_M2_vrr;
      Double I_ERI_G3xy_S_G2x2y_S_M1_vrr = PAY*I_ERI_F3x_S_G2x2y_S_M1_vrr+WPY*I_ERI_F3x_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_F3x_S_F2xy_S_M2_vrr;
      Double I_ERI_G3xz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_F3x_S_G2x2y_S_M1_vrr+WPZ*I_ERI_F3x_S_G2x2y_S_M2_vrr;
      Double I_ERI_G2x2y_S_G2x2y_S_M1_vrr = PAY*I_ERI_F2xy_S_G2x2y_S_M1_vrr+WPY*I_ERI_F2xy_S_G2x2y_S_M2_vrr+oned2z*I_ERI_D2x_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_F2xy_S_M2_vrr;
      Double I_ERI_G2x2z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_F2xz_S_G2x2y_S_M1_vrr+WPZ*I_ERI_F2xz_S_G2x2y_S_M2_vrr+oned2z*I_ERI_D2x_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G2x2y_S_M2_vrr;
      Double I_ERI_Gx3y_S_G2x2y_S_M1_vrr = PAX*I_ERI_F3y_S_G2x2y_S_M1_vrr+WPX*I_ERI_F3y_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Fx2y_S_M2_vrr;
      Double I_ERI_Gx3z_S_G2x2y_S_M1_vrr = PAX*I_ERI_F3z_S_G2x2y_S_M1_vrr+WPX*I_ERI_F3z_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Fx2y_S_M2_vrr;
      Double I_ERI_G4y_S_G2x2y_S_M1_vrr = PAY*I_ERI_F3y_S_G2x2y_S_M1_vrr+WPY*I_ERI_F3y_S_G2x2y_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G2x2y_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_F3y_S_F2xy_S_M2_vrr;
      Double I_ERI_G3yz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_F3y_S_G2x2y_S_M1_vrr+WPZ*I_ERI_F3y_S_G2x2y_S_M2_vrr;
      Double I_ERI_G2y2z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_F2yz_S_G2x2y_S_M1_vrr+WPZ*I_ERI_F2yz_S_G2x2y_S_M2_vrr+oned2z*I_ERI_D2y_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G2x2y_S_M2_vrr;
      Double I_ERI_Gy3z_S_G2x2y_S_M1_vrr = PAY*I_ERI_F3z_S_G2x2y_S_M1_vrr+WPY*I_ERI_F3z_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_F3z_S_F2xy_S_M2_vrr;
      Double I_ERI_G4z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_F3z_S_G2x2y_S_M1_vrr+WPZ*I_ERI_F3z_S_G2x2y_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G2x2y_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G2x2y_S_M2_vrr;
      Double I_ERI_G4x_S_G2xyz_S_M1_vrr = PAX*I_ERI_F3x_S_G2xyz_S_M1_vrr+WPX*I_ERI_F3x_S_G2xyz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G2xyz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M2_vrr;
      Double I_ERI_G3xy_S_G2xyz_S_M1_vrr = PAY*I_ERI_F3x_S_G2xyz_S_M1_vrr+WPY*I_ERI_F3x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_F3x_S_F2xz_S_M2_vrr;
      Double I_ERI_G3xz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_F3x_S_G2xyz_S_M1_vrr+WPZ*I_ERI_F3x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_F3x_S_F2xy_S_M2_vrr;
      Double I_ERI_G2x2y_S_G2xyz_S_M1_vrr = PAY*I_ERI_F2xy_S_G2xyz_S_M1_vrr+WPY*I_ERI_F2xy_S_G2xyz_S_M2_vrr+oned2z*I_ERI_D2x_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_F2xy_S_F2xz_S_M2_vrr;
      Double I_ERI_G2x2z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_F2xz_S_G2xyz_S_M1_vrr+WPZ*I_ERI_F2xz_S_G2xyz_S_M2_vrr+oned2z*I_ERI_D2x_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_F2xz_S_F2xy_S_M2_vrr;
      Double I_ERI_Gx3y_S_G2xyz_S_M1_vrr = PAX*I_ERI_F3y_S_G2xyz_S_M1_vrr+WPX*I_ERI_F3y_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M2_vrr;
      Double I_ERI_Gx3z_S_G2xyz_S_M1_vrr = PAX*I_ERI_F3z_S_G2xyz_S_M1_vrr+WPX*I_ERI_F3z_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M2_vrr;
      Double I_ERI_G4y_S_G2xyz_S_M1_vrr = PAY*I_ERI_F3y_S_G2xyz_S_M1_vrr+WPY*I_ERI_F3y_S_G2xyz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G2xyz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G2xyz_S_M2_vrr+oned2k*I_ERI_F3y_S_F2xz_S_M2_vrr;
      Double I_ERI_G3yz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_F3y_S_G2xyz_S_M1_vrr+WPZ*I_ERI_F3y_S_G2xyz_S_M2_vrr+oned2k*I_ERI_F3y_S_F2xy_S_M2_vrr;
      Double I_ERI_G2y2z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_F2yz_S_G2xyz_S_M1_vrr+WPZ*I_ERI_F2yz_S_G2xyz_S_M2_vrr+oned2z*I_ERI_D2y_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G2xyz_S_M2_vrr+oned2k*I_ERI_F2yz_S_F2xy_S_M2_vrr;
      Double I_ERI_Gy3z_S_G2xyz_S_M1_vrr = PAY*I_ERI_F3z_S_G2xyz_S_M1_vrr+WPY*I_ERI_F3z_S_G2xyz_S_M2_vrr+oned2k*I_ERI_F3z_S_F2xz_S_M2_vrr;
      Double I_ERI_G4z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_F3z_S_G2xyz_S_M1_vrr+WPZ*I_ERI_F3z_S_G2xyz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G2xyz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G2xyz_S_M2_vrr+oned2k*I_ERI_F3z_S_F2xy_S_M2_vrr;
      Double I_ERI_G4x_S_G2x2z_S_M1_vrr = PAX*I_ERI_F3x_S_G2x2z_S_M1_vrr+WPX*I_ERI_F3x_S_G2x2z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G2x2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Fx2z_S_M2_vrr;
      Double I_ERI_G3xy_S_G2x2z_S_M1_vrr = PAY*I_ERI_F3x_S_G2x2z_S_M1_vrr+WPY*I_ERI_F3x_S_G2x2z_S_M2_vrr;
      Double I_ERI_G3xz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_F3x_S_G2x2z_S_M1_vrr+WPZ*I_ERI_F3x_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_F3x_S_F2xz_S_M2_vrr;
      Double I_ERI_G2x2y_S_G2x2z_S_M1_vrr = PAY*I_ERI_F2xy_S_G2x2z_S_M1_vrr+WPY*I_ERI_F2xy_S_G2x2z_S_M2_vrr+oned2z*I_ERI_D2x_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G2x2z_S_M2_vrr;
      Double I_ERI_G2x2z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_F2xz_S_G2x2z_S_M1_vrr+WPZ*I_ERI_F2xz_S_G2x2z_S_M2_vrr+oned2z*I_ERI_D2x_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_F2xz_S_M2_vrr;
      Double I_ERI_Gx3y_S_G2x2z_S_M1_vrr = PAX*I_ERI_F3y_S_G2x2z_S_M1_vrr+WPX*I_ERI_F3y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Fx2z_S_M2_vrr;
      Double I_ERI_Gx3z_S_G2x2z_S_M1_vrr = PAX*I_ERI_F3z_S_G2x2z_S_M1_vrr+WPX*I_ERI_F3z_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Fx2z_S_M2_vrr;
      Double I_ERI_G4y_S_G2x2z_S_M1_vrr = PAY*I_ERI_F3y_S_G2x2z_S_M1_vrr+WPY*I_ERI_F3y_S_G2x2z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G2x2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G2x2z_S_M2_vrr;
      Double I_ERI_G3yz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_F3y_S_G2x2z_S_M1_vrr+WPZ*I_ERI_F3y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_F3y_S_F2xz_S_M2_vrr;
      Double I_ERI_G2y2z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_F2yz_S_G2x2z_S_M1_vrr+WPZ*I_ERI_F2yz_S_G2x2z_S_M2_vrr+oned2z*I_ERI_D2y_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_F2xz_S_M2_vrr;
      Double I_ERI_Gy3z_S_G2x2z_S_M1_vrr = PAY*I_ERI_F3z_S_G2x2z_S_M1_vrr+WPY*I_ERI_F3z_S_G2x2z_S_M2_vrr;
      Double I_ERI_G4z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_F3z_S_G2x2z_S_M1_vrr+WPZ*I_ERI_F3z_S_G2x2z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G2x2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_F3z_S_F2xz_S_M2_vrr;
      Double I_ERI_G4x_S_Gx3y_S_M1_vrr = PAX*I_ERI_F3x_S_Gx3y_S_M1_vrr+WPX*I_ERI_F3x_S_Gx3y_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Gx3y_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx3y_S_M2_vrr+oned2k*I_ERI_F3x_S_F3y_S_M2_vrr;
      Double I_ERI_G3xy_S_Gx3y_S_M1_vrr = PAY*I_ERI_F3x_S_Gx3y_S_M1_vrr+WPY*I_ERI_F3x_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_F3x_S_Fx2y_S_M2_vrr;
      Double I_ERI_G3xz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_F3x_S_Gx3y_S_M1_vrr+WPZ*I_ERI_F3x_S_Gx3y_S_M2_vrr;
      Double I_ERI_G2x2y_S_Gx3y_S_M1_vrr = PAY*I_ERI_F2xy_S_Gx3y_S_M1_vrr+WPY*I_ERI_F2xy_S_Gx3y_S_M2_vrr+oned2z*I_ERI_D2x_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_F2xy_S_Fx2y_S_M2_vrr;
      Double I_ERI_G2x2z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_F2xz_S_Gx3y_S_M1_vrr+WPZ*I_ERI_F2xz_S_Gx3y_S_M2_vrr+oned2z*I_ERI_D2x_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gx3y_S_M2_vrr;
      Double I_ERI_Gx3y_S_Gx3y_S_M1_vrr = PAX*I_ERI_F3y_S_Gx3y_S_M1_vrr+WPX*I_ERI_F3y_S_Gx3y_S_M2_vrr+oned2k*I_ERI_F3y_S_F3y_S_M2_vrr;
      Double I_ERI_Gx3z_S_Gx3y_S_M1_vrr = PAX*I_ERI_F3z_S_Gx3y_S_M1_vrr+WPX*I_ERI_F3z_S_Gx3y_S_M2_vrr+oned2k*I_ERI_F3z_S_F3y_S_M2_vrr;
      Double I_ERI_G4y_S_Gx3y_S_M1_vrr = PAY*I_ERI_F3y_S_Gx3y_S_M1_vrr+WPY*I_ERI_F3y_S_Gx3y_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Gx3y_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_F3y_S_Fx2y_S_M2_vrr;
      Double I_ERI_G3yz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_F3y_S_Gx3y_S_M1_vrr+WPZ*I_ERI_F3y_S_Gx3y_S_M2_vrr;
      Double I_ERI_G2y2z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_F2yz_S_Gx3y_S_M1_vrr+WPZ*I_ERI_F2yz_S_Gx3y_S_M2_vrr+oned2z*I_ERI_D2y_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Gx3y_S_M2_vrr;
      Double I_ERI_Gy3z_S_Gx3y_S_M1_vrr = PAY*I_ERI_F3z_S_Gx3y_S_M1_vrr+WPY*I_ERI_F3z_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_F3z_S_Fx2y_S_M2_vrr;
      Double I_ERI_G4z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_F3z_S_Gx3y_S_M1_vrr+WPZ*I_ERI_F3z_S_Gx3y_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Gx3y_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx3y_S_M2_vrr;
      Double I_ERI_G4x_S_Gx2yz_S_M1_vrr = PAX*I_ERI_F3x_S_Gx2yz_S_M1_vrr+WPX*I_ERI_F3x_S_Gx2yz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Gx2yz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_F3x_S_F2yz_S_M2_vrr;
      Double I_ERI_G3xy_S_Gx2yz_S_M1_vrr = PAY*I_ERI_F3x_S_Gx2yz_S_M1_vrr+WPY*I_ERI_F3x_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M2_vrr;
      Double I_ERI_G3xz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_F3x_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_F3x_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_F3x_S_Fx2y_S_M2_vrr;
      Double I_ERI_G2x2y_S_Gx2yz_S_M1_vrr = PAY*I_ERI_F2xy_S_Gx2yz_S_M1_vrr+WPY*I_ERI_F2xy_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_D2x_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Fxyz_S_M2_vrr;
      Double I_ERI_G2x2z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_F2xz_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_F2xz_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_D2x_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_F2xz_S_Fx2y_S_M2_vrr;
      Double I_ERI_Gx3y_S_Gx2yz_S_M1_vrr = PAX*I_ERI_F3y_S_Gx2yz_S_M1_vrr+WPX*I_ERI_F3y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_F3y_S_F2yz_S_M2_vrr;
      Double I_ERI_Gx3z_S_Gx2yz_S_M1_vrr = PAX*I_ERI_F3z_S_Gx2yz_S_M1_vrr+WPX*I_ERI_F3z_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_F3z_S_F2yz_S_M2_vrr;
      Double I_ERI_G4y_S_Gx2yz_S_M1_vrr = PAY*I_ERI_F3y_S_Gx2yz_S_M1_vrr+WPY*I_ERI_F3y_S_Gx2yz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Gx2yz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M2_vrr;
      Double I_ERI_G3yz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_F3y_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_F3y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_F3y_S_Fx2y_S_M2_vrr;
      Double I_ERI_G2y2z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_F2yz_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_F2yz_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_D2y_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_F2yz_S_Fx2y_S_M2_vrr;
      Double I_ERI_Gy3z_S_Gx2yz_S_M1_vrr = PAY*I_ERI_F3z_S_Gx2yz_S_M1_vrr+WPY*I_ERI_F3z_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M2_vrr;
      Double I_ERI_G4z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_F3z_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_F3z_S_Gx2yz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Gx2yz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_F3z_S_Fx2y_S_M2_vrr;
      Double I_ERI_G4x_S_Gxy2z_S_M1_vrr = PAX*I_ERI_F3x_S_Gxy2z_S_M1_vrr+WPX*I_ERI_F3x_S_Gxy2z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Gxy2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_F3x_S_Fy2z_S_M2_vrr;
      Double I_ERI_G3xy_S_Gxy2z_S_M1_vrr = PAY*I_ERI_F3x_S_Gxy2z_S_M1_vrr+WPY*I_ERI_F3x_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_F3x_S_Fx2z_S_M2_vrr;
      Double I_ERI_G3xz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_F3x_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_F3x_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M2_vrr;
      Double I_ERI_G2x2y_S_Gxy2z_S_M1_vrr = PAY*I_ERI_F2xy_S_Gxy2z_S_M1_vrr+WPY*I_ERI_F2xy_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_D2x_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_F2xy_S_Fx2z_S_M2_vrr;
      Double I_ERI_G2x2z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_F2xz_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_F2xz_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_D2x_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_Fxyz_S_M2_vrr;
      Double I_ERI_Gx3y_S_Gxy2z_S_M1_vrr = PAX*I_ERI_F3y_S_Gxy2z_S_M1_vrr+WPX*I_ERI_F3y_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_F3y_S_Fy2z_S_M2_vrr;
      Double I_ERI_Gx3z_S_Gxy2z_S_M1_vrr = PAX*I_ERI_F3z_S_Gxy2z_S_M1_vrr+WPX*I_ERI_F3z_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_F3z_S_Fy2z_S_M2_vrr;
      Double I_ERI_G4y_S_Gxy2z_S_M1_vrr = PAY*I_ERI_F3y_S_Gxy2z_S_M1_vrr+WPY*I_ERI_F3y_S_Gxy2z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Gxy2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_F3y_S_Fx2z_S_M2_vrr;
      Double I_ERI_G3yz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_F3y_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_F3y_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M2_vrr;
      Double I_ERI_G2y2z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_F2yz_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_F2yz_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_D2y_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_Fxyz_S_M2_vrr;
      Double I_ERI_Gy3z_S_Gxy2z_S_M1_vrr = PAY*I_ERI_F3z_S_Gxy2z_S_M1_vrr+WPY*I_ERI_F3z_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_F3z_S_Fx2z_S_M2_vrr;
      Double I_ERI_G4z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_F3z_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_F3z_S_Gxy2z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Gxy2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M2_vrr;
      Double I_ERI_G4x_S_Gx3z_S_M1_vrr = PAX*I_ERI_F3x_S_Gx3z_S_M1_vrr+WPX*I_ERI_F3x_S_Gx3z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Gx3z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx3z_S_M2_vrr+oned2k*I_ERI_F3x_S_F3z_S_M2_vrr;
      Double I_ERI_G3xy_S_Gx3z_S_M1_vrr = PAY*I_ERI_F3x_S_Gx3z_S_M1_vrr+WPY*I_ERI_F3x_S_Gx3z_S_M2_vrr;
      Double I_ERI_G3xz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_F3x_S_Gx3z_S_M1_vrr+WPZ*I_ERI_F3x_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_F3x_S_Fx2z_S_M2_vrr;
      Double I_ERI_G2x2y_S_Gx3z_S_M1_vrr = PAY*I_ERI_F2xy_S_Gx3z_S_M1_vrr+WPY*I_ERI_F2xy_S_Gx3z_S_M2_vrr+oned2z*I_ERI_D2x_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gx3z_S_M2_vrr;
      Double I_ERI_G2x2z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_F2xz_S_Gx3z_S_M1_vrr+WPZ*I_ERI_F2xz_S_Gx3z_S_M2_vrr+oned2z*I_ERI_D2x_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_F2xz_S_Fx2z_S_M2_vrr;
      Double I_ERI_Gx3y_S_Gx3z_S_M1_vrr = PAX*I_ERI_F3y_S_Gx3z_S_M1_vrr+WPX*I_ERI_F3y_S_Gx3z_S_M2_vrr+oned2k*I_ERI_F3y_S_F3z_S_M2_vrr;
      Double I_ERI_Gx3z_S_Gx3z_S_M1_vrr = PAX*I_ERI_F3z_S_Gx3z_S_M1_vrr+WPX*I_ERI_F3z_S_Gx3z_S_M2_vrr+oned2k*I_ERI_F3z_S_F3z_S_M2_vrr;
      Double I_ERI_G4y_S_Gx3z_S_M1_vrr = PAY*I_ERI_F3y_S_Gx3z_S_M1_vrr+WPY*I_ERI_F3y_S_Gx3z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Gx3z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx3z_S_M2_vrr;
      Double I_ERI_G3yz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_F3y_S_Gx3z_S_M1_vrr+WPZ*I_ERI_F3y_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_F3y_S_Fx2z_S_M2_vrr;
      Double I_ERI_G2y2z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_F2yz_S_Gx3z_S_M1_vrr+WPZ*I_ERI_F2yz_S_Gx3z_S_M2_vrr+oned2z*I_ERI_D2y_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_F2yz_S_Fx2z_S_M2_vrr;
      Double I_ERI_Gy3z_S_Gx3z_S_M1_vrr = PAY*I_ERI_F3z_S_Gx3z_S_M1_vrr+WPY*I_ERI_F3z_S_Gx3z_S_M2_vrr;
      Double I_ERI_G4z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_F3z_S_Gx3z_S_M1_vrr+WPZ*I_ERI_F3z_S_Gx3z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Gx3z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_F3z_S_Fx2z_S_M2_vrr;
      Double I_ERI_G4x_S_G4y_S_M1_vrr = PAX*I_ERI_F3x_S_G4y_S_M1_vrr+WPX*I_ERI_F3x_S_G4y_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G4y_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G4y_S_M2_vrr;
      Double I_ERI_G3xy_S_G4y_S_M1_vrr = PAY*I_ERI_F3x_S_G4y_S_M1_vrr+WPY*I_ERI_F3x_S_G4y_S_M2_vrr+4*oned2k*I_ERI_F3x_S_F3y_S_M2_vrr;
      Double I_ERI_G3xz_S_G4y_S_M1_vrr = PAZ*I_ERI_F3x_S_G4y_S_M1_vrr+WPZ*I_ERI_F3x_S_G4y_S_M2_vrr;
      Double I_ERI_G2x2y_S_G4y_S_M1_vrr = PAY*I_ERI_F2xy_S_G4y_S_M1_vrr+WPY*I_ERI_F2xy_S_G4y_S_M2_vrr+oned2z*I_ERI_D2x_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G4y_S_M2_vrr+4*oned2k*I_ERI_F2xy_S_F3y_S_M2_vrr;
      Double I_ERI_G2x2z_S_G4y_S_M1_vrr = PAZ*I_ERI_F2xz_S_G4y_S_M1_vrr+WPZ*I_ERI_F2xz_S_G4y_S_M2_vrr+oned2z*I_ERI_D2x_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G4y_S_M2_vrr;
      Double I_ERI_Gx3y_S_G4y_S_M1_vrr = PAX*I_ERI_F3y_S_G4y_S_M1_vrr+WPX*I_ERI_F3y_S_G4y_S_M2_vrr;
      Double I_ERI_Gx3z_S_G4y_S_M1_vrr = PAX*I_ERI_F3z_S_G4y_S_M1_vrr+WPX*I_ERI_F3z_S_G4y_S_M2_vrr;
      Double I_ERI_G4y_S_G4y_S_M1_vrr = PAY*I_ERI_F3y_S_G4y_S_M1_vrr+WPY*I_ERI_F3y_S_G4y_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G4y_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G4y_S_M2_vrr+4*oned2k*I_ERI_F3y_S_F3y_S_M2_vrr;
      Double I_ERI_G3yz_S_G4y_S_M1_vrr = PAZ*I_ERI_F3y_S_G4y_S_M1_vrr+WPZ*I_ERI_F3y_S_G4y_S_M2_vrr;
      Double I_ERI_G2y2z_S_G4y_S_M1_vrr = PAZ*I_ERI_F2yz_S_G4y_S_M1_vrr+WPZ*I_ERI_F2yz_S_G4y_S_M2_vrr+oned2z*I_ERI_D2y_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G4y_S_M2_vrr;
      Double I_ERI_Gy3z_S_G4y_S_M1_vrr = PAY*I_ERI_F3z_S_G4y_S_M1_vrr+WPY*I_ERI_F3z_S_G4y_S_M2_vrr+4*oned2k*I_ERI_F3z_S_F3y_S_M2_vrr;
      Double I_ERI_G4z_S_G4y_S_M1_vrr = PAZ*I_ERI_F3z_S_G4y_S_M1_vrr+WPZ*I_ERI_F3z_S_G4y_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G4y_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G4y_S_M2_vrr;
      Double I_ERI_G4x_S_G3yz_S_M1_vrr = PAX*I_ERI_F3x_S_G3yz_S_M1_vrr+WPX*I_ERI_F3x_S_G3yz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G3yz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G3yz_S_M2_vrr;
      Double I_ERI_G3xy_S_G3yz_S_M1_vrr = PAY*I_ERI_F3x_S_G3yz_S_M1_vrr+WPY*I_ERI_F3x_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_F3x_S_F2yz_S_M2_vrr;
      Double I_ERI_G3xz_S_G3yz_S_M1_vrr = PAZ*I_ERI_F3x_S_G3yz_S_M1_vrr+WPZ*I_ERI_F3x_S_G3yz_S_M2_vrr+oned2k*I_ERI_F3x_S_F3y_S_M2_vrr;
      Double I_ERI_G2x2y_S_G3yz_S_M1_vrr = PAY*I_ERI_F2xy_S_G3yz_S_M1_vrr+WPY*I_ERI_F2xy_S_G3yz_S_M2_vrr+oned2z*I_ERI_D2x_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_F2xy_S_F2yz_S_M2_vrr;
      Double I_ERI_G2x2z_S_G3yz_S_M1_vrr = PAZ*I_ERI_F2xz_S_G3yz_S_M1_vrr+WPZ*I_ERI_F2xz_S_G3yz_S_M2_vrr+oned2z*I_ERI_D2x_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G3yz_S_M2_vrr+oned2k*I_ERI_F2xz_S_F3y_S_M2_vrr;
      Double I_ERI_Gx3y_S_G3yz_S_M1_vrr = PAX*I_ERI_F3y_S_G3yz_S_M1_vrr+WPX*I_ERI_F3y_S_G3yz_S_M2_vrr;
      Double I_ERI_Gx3z_S_G3yz_S_M1_vrr = PAX*I_ERI_F3z_S_G3yz_S_M1_vrr+WPX*I_ERI_F3z_S_G3yz_S_M2_vrr;
      Double I_ERI_G4y_S_G3yz_S_M1_vrr = PAY*I_ERI_F3y_S_G3yz_S_M1_vrr+WPY*I_ERI_F3y_S_G3yz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G3yz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_F3y_S_F2yz_S_M2_vrr;
      Double I_ERI_G3yz_S_G3yz_S_M1_vrr = PAZ*I_ERI_F3y_S_G3yz_S_M1_vrr+WPZ*I_ERI_F3y_S_G3yz_S_M2_vrr+oned2k*I_ERI_F3y_S_F3y_S_M2_vrr;
      Double I_ERI_G2y2z_S_G3yz_S_M1_vrr = PAZ*I_ERI_F2yz_S_G3yz_S_M1_vrr+WPZ*I_ERI_F2yz_S_G3yz_S_M2_vrr+oned2z*I_ERI_D2y_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G3yz_S_M2_vrr+oned2k*I_ERI_F2yz_S_F3y_S_M2_vrr;
      Double I_ERI_Gy3z_S_G3yz_S_M1_vrr = PAY*I_ERI_F3z_S_G3yz_S_M1_vrr+WPY*I_ERI_F3z_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_F3z_S_F2yz_S_M2_vrr;
      Double I_ERI_G4z_S_G3yz_S_M1_vrr = PAZ*I_ERI_F3z_S_G3yz_S_M1_vrr+WPZ*I_ERI_F3z_S_G3yz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G3yz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G3yz_S_M2_vrr+oned2k*I_ERI_F3z_S_F3y_S_M2_vrr;
      Double I_ERI_G4x_S_G2y2z_S_M1_vrr = PAX*I_ERI_F3x_S_G2y2z_S_M1_vrr+WPX*I_ERI_F3x_S_G2y2z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G2y2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G2y2z_S_M2_vrr;
      Double I_ERI_G3xy_S_G2y2z_S_M1_vrr = PAY*I_ERI_F3x_S_G2y2z_S_M1_vrr+WPY*I_ERI_F3x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_F3x_S_Fy2z_S_M2_vrr;
      Double I_ERI_G3xz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_F3x_S_G2y2z_S_M1_vrr+WPZ*I_ERI_F3x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_F3x_S_F2yz_S_M2_vrr;
      Double I_ERI_G2x2y_S_G2y2z_S_M1_vrr = PAY*I_ERI_F2xy_S_G2y2z_S_M1_vrr+WPY*I_ERI_F2xy_S_G2y2z_S_M2_vrr+oned2z*I_ERI_D2x_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Fy2z_S_M2_vrr;
      Double I_ERI_G2x2z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_F2xz_S_G2y2z_S_M1_vrr+WPZ*I_ERI_F2xz_S_G2y2z_S_M2_vrr+oned2z*I_ERI_D2x_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_F2yz_S_M2_vrr;
      Double I_ERI_Gx3y_S_G2y2z_S_M1_vrr = PAX*I_ERI_F3y_S_G2y2z_S_M1_vrr+WPX*I_ERI_F3y_S_G2y2z_S_M2_vrr;
      Double I_ERI_Gx3z_S_G2y2z_S_M1_vrr = PAX*I_ERI_F3z_S_G2y2z_S_M1_vrr+WPX*I_ERI_F3z_S_G2y2z_S_M2_vrr;
      Double I_ERI_G4y_S_G2y2z_S_M1_vrr = PAY*I_ERI_F3y_S_G2y2z_S_M1_vrr+WPY*I_ERI_F3y_S_G2y2z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G2y2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_F3y_S_Fy2z_S_M2_vrr;
      Double I_ERI_G3yz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_F3y_S_G2y2z_S_M1_vrr+WPZ*I_ERI_F3y_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_F3y_S_F2yz_S_M2_vrr;
      Double I_ERI_G2y2z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_F2yz_S_G2y2z_S_M1_vrr+WPZ*I_ERI_F2yz_S_G2y2z_S_M2_vrr+oned2z*I_ERI_D2y_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_F2yz_S_M2_vrr;
      Double I_ERI_Gy3z_S_G2y2z_S_M1_vrr = PAY*I_ERI_F3z_S_G2y2z_S_M1_vrr+WPY*I_ERI_F3z_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_F3z_S_Fy2z_S_M2_vrr;
      Double I_ERI_G4z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_F3z_S_G2y2z_S_M1_vrr+WPZ*I_ERI_F3z_S_G2y2z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G2y2z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_F3z_S_F2yz_S_M2_vrr;
      Double I_ERI_G4x_S_Gy3z_S_M1_vrr = PAX*I_ERI_F3x_S_Gy3z_S_M1_vrr+WPX*I_ERI_F3x_S_Gy3z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Gy3z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Gy3z_S_M2_vrr;
      Double I_ERI_G3xy_S_Gy3z_S_M1_vrr = PAY*I_ERI_F3x_S_Gy3z_S_M1_vrr+WPY*I_ERI_F3x_S_Gy3z_S_M2_vrr+oned2k*I_ERI_F3x_S_F3z_S_M2_vrr;
      Double I_ERI_G3xz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_F3x_S_Gy3z_S_M1_vrr+WPZ*I_ERI_F3x_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_F3x_S_Fy2z_S_M2_vrr;
      Double I_ERI_G2x2y_S_Gy3z_S_M1_vrr = PAY*I_ERI_F2xy_S_Gy3z_S_M1_vrr+WPY*I_ERI_F2xy_S_Gy3z_S_M2_vrr+oned2z*I_ERI_D2x_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gy3z_S_M2_vrr+oned2k*I_ERI_F2xy_S_F3z_S_M2_vrr;
      Double I_ERI_G2x2z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_F2xz_S_Gy3z_S_M1_vrr+WPZ*I_ERI_F2xz_S_Gy3z_S_M2_vrr+oned2z*I_ERI_D2x_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_F2xz_S_Fy2z_S_M2_vrr;
      Double I_ERI_Gx3y_S_Gy3z_S_M1_vrr = PAX*I_ERI_F3y_S_Gy3z_S_M1_vrr+WPX*I_ERI_F3y_S_Gy3z_S_M2_vrr;
      Double I_ERI_Gx3z_S_Gy3z_S_M1_vrr = PAX*I_ERI_F3z_S_Gy3z_S_M1_vrr+WPX*I_ERI_F3z_S_Gy3z_S_M2_vrr;
      Double I_ERI_G4y_S_Gy3z_S_M1_vrr = PAY*I_ERI_F3y_S_Gy3z_S_M1_vrr+WPY*I_ERI_F3y_S_Gy3z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Gy3z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Gy3z_S_M2_vrr+oned2k*I_ERI_F3y_S_F3z_S_M2_vrr;
      Double I_ERI_G3yz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_F3y_S_Gy3z_S_M1_vrr+WPZ*I_ERI_F3y_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_F3y_S_Fy2z_S_M2_vrr;
      Double I_ERI_G2y2z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_F2yz_S_Gy3z_S_M1_vrr+WPZ*I_ERI_F2yz_S_Gy3z_S_M2_vrr+oned2z*I_ERI_D2y_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_F2yz_S_Fy2z_S_M2_vrr;
      Double I_ERI_Gy3z_S_Gy3z_S_M1_vrr = PAY*I_ERI_F3z_S_Gy3z_S_M1_vrr+WPY*I_ERI_F3z_S_Gy3z_S_M2_vrr+oned2k*I_ERI_F3z_S_F3z_S_M2_vrr;
      Double I_ERI_G4z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_F3z_S_Gy3z_S_M1_vrr+WPZ*I_ERI_F3z_S_Gy3z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Gy3z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_F3z_S_Fy2z_S_M2_vrr;
      Double I_ERI_G4x_S_G4z_S_M1_vrr = PAX*I_ERI_F3x_S_G4z_S_M1_vrr+WPX*I_ERI_F3x_S_G4z_S_M2_vrr+3*oned2z*I_ERI_D2x_S_G4z_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_G4z_S_M2_vrr;
      Double I_ERI_G3xy_S_G4z_S_M1_vrr = PAY*I_ERI_F3x_S_G4z_S_M1_vrr+WPY*I_ERI_F3x_S_G4z_S_M2_vrr;
      Double I_ERI_G3xz_S_G4z_S_M1_vrr = PAZ*I_ERI_F3x_S_G4z_S_M1_vrr+WPZ*I_ERI_F3x_S_G4z_S_M2_vrr+4*oned2k*I_ERI_F3x_S_F3z_S_M2_vrr;
      Double I_ERI_G2x2y_S_G4z_S_M1_vrr = PAY*I_ERI_F2xy_S_G4z_S_M1_vrr+WPY*I_ERI_F2xy_S_G4z_S_M2_vrr+oned2z*I_ERI_D2x_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G4z_S_M2_vrr;
      Double I_ERI_G2x2z_S_G4z_S_M1_vrr = PAZ*I_ERI_F2xz_S_G4z_S_M1_vrr+WPZ*I_ERI_F2xz_S_G4z_S_M2_vrr+oned2z*I_ERI_D2x_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_G4z_S_M2_vrr+4*oned2k*I_ERI_F2xz_S_F3z_S_M2_vrr;
      Double I_ERI_Gx3y_S_G4z_S_M1_vrr = PAX*I_ERI_F3y_S_G4z_S_M1_vrr+WPX*I_ERI_F3y_S_G4z_S_M2_vrr;
      Double I_ERI_Gx3z_S_G4z_S_M1_vrr = PAX*I_ERI_F3z_S_G4z_S_M1_vrr+WPX*I_ERI_F3z_S_G4z_S_M2_vrr;
      Double I_ERI_G4y_S_G4z_S_M1_vrr = PAY*I_ERI_F3y_S_G4z_S_M1_vrr+WPY*I_ERI_F3y_S_G4z_S_M2_vrr+3*oned2z*I_ERI_D2y_S_G4z_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_G4z_S_M2_vrr;
      Double I_ERI_G3yz_S_G4z_S_M1_vrr = PAZ*I_ERI_F3y_S_G4z_S_M1_vrr+WPZ*I_ERI_F3y_S_G4z_S_M2_vrr+4*oned2k*I_ERI_F3y_S_F3z_S_M2_vrr;
      Double I_ERI_G2y2z_S_G4z_S_M1_vrr = PAZ*I_ERI_F2yz_S_G4z_S_M1_vrr+WPZ*I_ERI_F2yz_S_G4z_S_M2_vrr+oned2z*I_ERI_D2y_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_G4z_S_M2_vrr+4*oned2k*I_ERI_F2yz_S_F3z_S_M2_vrr;
      Double I_ERI_Gy3z_S_G4z_S_M1_vrr = PAY*I_ERI_F3z_S_G4z_S_M1_vrr+WPY*I_ERI_F3z_S_G4z_S_M2_vrr;
      Double I_ERI_G4z_S_G4z_S_M1_vrr = PAZ*I_ERI_F3z_S_G4z_S_M1_vrr+WPZ*I_ERI_F3z_S_G4z_S_M2_vrr+3*oned2z*I_ERI_D2z_S_G4z_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_G4z_S_M2_vrr+4*oned2k*I_ERI_F3z_S_F3z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 50 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_D_S_M2
       ************************************************************/
      Double I_ERI_H5x_S_F3x_S_M1_vrr = PAX*I_ERI_G4x_S_F3x_S_M1_vrr+WPX*I_ERI_G4x_S_F3x_S_M2_vrr+4*oned2z*I_ERI_F3x_S_F3x_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_F3x_S_M2_vrr+3*oned2k*I_ERI_G4x_S_D2x_S_M2_vrr;
      Double I_ERI_H4xy_S_F3x_S_M1_vrr = PAY*I_ERI_G4x_S_F3x_S_M1_vrr+WPY*I_ERI_G4x_S_F3x_S_M2_vrr;
      Double I_ERI_H4xz_S_F3x_S_M1_vrr = PAZ*I_ERI_G4x_S_F3x_S_M1_vrr+WPZ*I_ERI_G4x_S_F3x_S_M2_vrr;
      Double I_ERI_H3x2y_S_F3x_S_M1_vrr = PAY*I_ERI_G3xy_S_F3x_S_M1_vrr+WPY*I_ERI_G3xy_S_F3x_S_M2_vrr+oned2z*I_ERI_F3x_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F3x_S_M2_vrr;
      Double I_ERI_H3x2z_S_F3x_S_M1_vrr = PAZ*I_ERI_G3xz_S_F3x_S_M1_vrr+WPZ*I_ERI_G3xz_S_F3x_S_M2_vrr+oned2z*I_ERI_F3x_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F3x_S_M2_vrr;
      Double I_ERI_H2x3y_S_F3x_S_M1_vrr = PAX*I_ERI_Gx3y_S_F3x_S_M1_vrr+WPX*I_ERI_Gx3y_S_F3x_S_M2_vrr+oned2z*I_ERI_F3y_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F3x_S_M2_vrr+3*oned2k*I_ERI_Gx3y_S_D2x_S_M2_vrr;
      Double I_ERI_H2x2yz_S_F3x_S_M1_vrr = PAZ*I_ERI_G2x2y_S_F3x_S_M1_vrr+WPZ*I_ERI_G2x2y_S_F3x_S_M2_vrr;
      Double I_ERI_H2x3z_S_F3x_S_M1_vrr = PAX*I_ERI_Gx3z_S_F3x_S_M1_vrr+WPX*I_ERI_Gx3z_S_F3x_S_M2_vrr+oned2z*I_ERI_F3z_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F3x_S_M2_vrr+3*oned2k*I_ERI_Gx3z_S_D2x_S_M2_vrr;
      Double I_ERI_Hx4y_S_F3x_S_M1_vrr = PAX*I_ERI_G4y_S_F3x_S_M1_vrr+WPX*I_ERI_G4y_S_F3x_S_M2_vrr+3*oned2k*I_ERI_G4y_S_D2x_S_M2_vrr;
      Double I_ERI_Hx4z_S_F3x_S_M1_vrr = PAX*I_ERI_G4z_S_F3x_S_M1_vrr+WPX*I_ERI_G4z_S_F3x_S_M2_vrr+3*oned2k*I_ERI_G4z_S_D2x_S_M2_vrr;
      Double I_ERI_H5y_S_F3x_S_M1_vrr = PAY*I_ERI_G4y_S_F3x_S_M1_vrr+WPY*I_ERI_G4y_S_F3x_S_M2_vrr+4*oned2z*I_ERI_F3y_S_F3x_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_F3x_S_M2_vrr;
      Double I_ERI_H4yz_S_F3x_S_M1_vrr = PAZ*I_ERI_G4y_S_F3x_S_M1_vrr+WPZ*I_ERI_G4y_S_F3x_S_M2_vrr;
      Double I_ERI_H3y2z_S_F3x_S_M1_vrr = PAZ*I_ERI_G3yz_S_F3x_S_M1_vrr+WPZ*I_ERI_G3yz_S_F3x_S_M2_vrr+oned2z*I_ERI_F3y_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F3x_S_M2_vrr;
      Double I_ERI_H2y3z_S_F3x_S_M1_vrr = PAY*I_ERI_Gy3z_S_F3x_S_M1_vrr+WPY*I_ERI_Gy3z_S_F3x_S_M2_vrr+oned2z*I_ERI_F3z_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F3x_S_M2_vrr;
      Double I_ERI_Hy4z_S_F3x_S_M1_vrr = PAY*I_ERI_G4z_S_F3x_S_M1_vrr+WPY*I_ERI_G4z_S_F3x_S_M2_vrr;
      Double I_ERI_H5z_S_F3x_S_M1_vrr = PAZ*I_ERI_G4z_S_F3x_S_M1_vrr+WPZ*I_ERI_G4z_S_F3x_S_M2_vrr+4*oned2z*I_ERI_F3z_S_F3x_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_F3x_S_M2_vrr;
      Double I_ERI_H5x_S_F2xy_S_M1_vrr = PAX*I_ERI_G4x_S_F2xy_S_M1_vrr+WPX*I_ERI_G4x_S_F2xy_S_M2_vrr+4*oned2z*I_ERI_F3x_S_F2xy_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Dxy_S_M2_vrr;
      Double I_ERI_H4xy_S_F2xy_S_M1_vrr = PAY*I_ERI_G4x_S_F2xy_S_M1_vrr+WPY*I_ERI_G4x_S_F2xy_S_M2_vrr+oned2k*I_ERI_G4x_S_D2x_S_M2_vrr;
      Double I_ERI_H4xz_S_F2xy_S_M1_vrr = PAZ*I_ERI_G4x_S_F2xy_S_M1_vrr+WPZ*I_ERI_G4x_S_F2xy_S_M2_vrr;
      Double I_ERI_H3x2y_S_F2xy_S_M1_vrr = PAY*I_ERI_G3xy_S_F2xy_S_M1_vrr+WPY*I_ERI_G3xy_S_F2xy_S_M2_vrr+oned2z*I_ERI_F3x_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F2xy_S_M2_vrr+oned2k*I_ERI_G3xy_S_D2x_S_M2_vrr;
      Double I_ERI_H3x2z_S_F2xy_S_M1_vrr = PAZ*I_ERI_G3xz_S_F2xy_S_M1_vrr+WPZ*I_ERI_G3xz_S_F2xy_S_M2_vrr+oned2z*I_ERI_F3x_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F2xy_S_M2_vrr;
      Double I_ERI_H2x3y_S_F2xy_S_M1_vrr = PAX*I_ERI_Gx3y_S_F2xy_S_M1_vrr+WPX*I_ERI_Gx3y_S_F2xy_S_M2_vrr+oned2z*I_ERI_F3y_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_Gx3y_S_Dxy_S_M2_vrr;
      Double I_ERI_H2x2yz_S_F2xy_S_M1_vrr = PAZ*I_ERI_G2x2y_S_F2xy_S_M1_vrr+WPZ*I_ERI_G2x2y_S_F2xy_S_M2_vrr;
      Double I_ERI_H2x3z_S_F2xy_S_M1_vrr = PAX*I_ERI_Gx3z_S_F2xy_S_M1_vrr+WPX*I_ERI_Gx3z_S_F2xy_S_M2_vrr+oned2z*I_ERI_F3z_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_Gx3z_S_Dxy_S_M2_vrr;
      Double I_ERI_Hx4y_S_F2xy_S_M1_vrr = PAX*I_ERI_G4y_S_F2xy_S_M1_vrr+WPX*I_ERI_G4y_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Dxy_S_M2_vrr;
      Double I_ERI_Hx4z_S_F2xy_S_M1_vrr = PAX*I_ERI_G4z_S_F2xy_S_M1_vrr+WPX*I_ERI_G4z_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Dxy_S_M2_vrr;
      Double I_ERI_H5y_S_F2xy_S_M1_vrr = PAY*I_ERI_G4y_S_F2xy_S_M1_vrr+WPY*I_ERI_G4y_S_F2xy_S_M2_vrr+4*oned2z*I_ERI_F3y_S_F2xy_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_F2xy_S_M2_vrr+oned2k*I_ERI_G4y_S_D2x_S_M2_vrr;
      Double I_ERI_H4yz_S_F2xy_S_M1_vrr = PAZ*I_ERI_G4y_S_F2xy_S_M1_vrr+WPZ*I_ERI_G4y_S_F2xy_S_M2_vrr;
      Double I_ERI_H3y2z_S_F2xy_S_M1_vrr = PAZ*I_ERI_G3yz_S_F2xy_S_M1_vrr+WPZ*I_ERI_G3yz_S_F2xy_S_M2_vrr+oned2z*I_ERI_F3y_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F2xy_S_M2_vrr;
      Double I_ERI_H2y3z_S_F2xy_S_M1_vrr = PAY*I_ERI_Gy3z_S_F2xy_S_M1_vrr+WPY*I_ERI_Gy3z_S_F2xy_S_M2_vrr+oned2z*I_ERI_F3z_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F2xy_S_M2_vrr+oned2k*I_ERI_Gy3z_S_D2x_S_M2_vrr;
      Double I_ERI_Hy4z_S_F2xy_S_M1_vrr = PAY*I_ERI_G4z_S_F2xy_S_M1_vrr+WPY*I_ERI_G4z_S_F2xy_S_M2_vrr+oned2k*I_ERI_G4z_S_D2x_S_M2_vrr;
      Double I_ERI_H5z_S_F2xy_S_M1_vrr = PAZ*I_ERI_G4z_S_F2xy_S_M1_vrr+WPZ*I_ERI_G4z_S_F2xy_S_M2_vrr+4*oned2z*I_ERI_F3z_S_F2xy_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_F2xy_S_M2_vrr;
      Double I_ERI_H5x_S_F2xz_S_M1_vrr = PAX*I_ERI_G4x_S_F2xz_S_M1_vrr+WPX*I_ERI_G4x_S_F2xz_S_M2_vrr+4*oned2z*I_ERI_F3x_S_F2xz_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Dxz_S_M2_vrr;
      Double I_ERI_H4xy_S_F2xz_S_M1_vrr = PAY*I_ERI_G4x_S_F2xz_S_M1_vrr+WPY*I_ERI_G4x_S_F2xz_S_M2_vrr;
      Double I_ERI_H4xz_S_F2xz_S_M1_vrr = PAZ*I_ERI_G4x_S_F2xz_S_M1_vrr+WPZ*I_ERI_G4x_S_F2xz_S_M2_vrr+oned2k*I_ERI_G4x_S_D2x_S_M2_vrr;
      Double I_ERI_H3x2y_S_F2xz_S_M1_vrr = PAY*I_ERI_G3xy_S_F2xz_S_M1_vrr+WPY*I_ERI_G3xy_S_F2xz_S_M2_vrr+oned2z*I_ERI_F3x_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F2xz_S_M2_vrr;
      Double I_ERI_H3x2z_S_F2xz_S_M1_vrr = PAZ*I_ERI_G3xz_S_F2xz_S_M1_vrr+WPZ*I_ERI_G3xz_S_F2xz_S_M2_vrr+oned2z*I_ERI_F3x_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F2xz_S_M2_vrr+oned2k*I_ERI_G3xz_S_D2x_S_M2_vrr;
      Double I_ERI_H2x3y_S_F2xz_S_M1_vrr = PAX*I_ERI_Gx3y_S_F2xz_S_M1_vrr+WPX*I_ERI_Gx3y_S_F2xz_S_M2_vrr+oned2z*I_ERI_F3y_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_Gx3y_S_Dxz_S_M2_vrr;
      Double I_ERI_H2x2yz_S_F2xz_S_M1_vrr = PAZ*I_ERI_G2x2y_S_F2xz_S_M1_vrr+WPZ*I_ERI_G2x2y_S_F2xz_S_M2_vrr+oned2k*I_ERI_G2x2y_S_D2x_S_M2_vrr;
      Double I_ERI_H2x3z_S_F2xz_S_M1_vrr = PAX*I_ERI_Gx3z_S_F2xz_S_M1_vrr+WPX*I_ERI_Gx3z_S_F2xz_S_M2_vrr+oned2z*I_ERI_F3z_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_Gx3z_S_Dxz_S_M2_vrr;
      Double I_ERI_Hx4y_S_F2xz_S_M1_vrr = PAX*I_ERI_G4y_S_F2xz_S_M1_vrr+WPX*I_ERI_G4y_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Dxz_S_M2_vrr;
      Double I_ERI_Hx4z_S_F2xz_S_M1_vrr = PAX*I_ERI_G4z_S_F2xz_S_M1_vrr+WPX*I_ERI_G4z_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Dxz_S_M2_vrr;
      Double I_ERI_H5y_S_F2xz_S_M1_vrr = PAY*I_ERI_G4y_S_F2xz_S_M1_vrr+WPY*I_ERI_G4y_S_F2xz_S_M2_vrr+4*oned2z*I_ERI_F3y_S_F2xz_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_F2xz_S_M2_vrr;
      Double I_ERI_H4yz_S_F2xz_S_M1_vrr = PAZ*I_ERI_G4y_S_F2xz_S_M1_vrr+WPZ*I_ERI_G4y_S_F2xz_S_M2_vrr+oned2k*I_ERI_G4y_S_D2x_S_M2_vrr;
      Double I_ERI_H3y2z_S_F2xz_S_M1_vrr = PAZ*I_ERI_G3yz_S_F2xz_S_M1_vrr+WPZ*I_ERI_G3yz_S_F2xz_S_M2_vrr+oned2z*I_ERI_F3y_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F2xz_S_M2_vrr+oned2k*I_ERI_G3yz_S_D2x_S_M2_vrr;
      Double I_ERI_H2y3z_S_F2xz_S_M1_vrr = PAY*I_ERI_Gy3z_S_F2xz_S_M1_vrr+WPY*I_ERI_Gy3z_S_F2xz_S_M2_vrr+oned2z*I_ERI_F3z_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F2xz_S_M2_vrr;
      Double I_ERI_Hy4z_S_F2xz_S_M1_vrr = PAY*I_ERI_G4z_S_F2xz_S_M1_vrr+WPY*I_ERI_G4z_S_F2xz_S_M2_vrr;
      Double I_ERI_H5z_S_F2xz_S_M1_vrr = PAZ*I_ERI_G4z_S_F2xz_S_M1_vrr+WPZ*I_ERI_G4z_S_F2xz_S_M2_vrr+4*oned2z*I_ERI_F3z_S_F2xz_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_F2xz_S_M2_vrr+oned2k*I_ERI_G4z_S_D2x_S_M2_vrr;
      Double I_ERI_H5x_S_Fx2y_S_M1_vrr = PAX*I_ERI_G4x_S_Fx2y_S_M1_vrr+WPX*I_ERI_G4x_S_Fx2y_S_M2_vrr+4*oned2z*I_ERI_F3x_S_Fx2y_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_Fx2y_S_M2_vrr+oned2k*I_ERI_G4x_S_D2y_S_M2_vrr;
      Double I_ERI_H4xy_S_Fx2y_S_M1_vrr = PAY*I_ERI_G4x_S_Fx2y_S_M1_vrr+WPY*I_ERI_G4x_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Dxy_S_M2_vrr;
      Double I_ERI_H4xz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_G4x_S_Fx2y_S_M1_vrr+WPZ*I_ERI_G4x_S_Fx2y_S_M2_vrr;
      Double I_ERI_H3x2y_S_Fx2y_S_M1_vrr = PAY*I_ERI_G3xy_S_Fx2y_S_M1_vrr+WPY*I_ERI_G3xy_S_Fx2y_S_M2_vrr+oned2z*I_ERI_F3x_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_G3xy_S_Dxy_S_M2_vrr;
      Double I_ERI_H3x2z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_G3xz_S_Fx2y_S_M1_vrr+WPZ*I_ERI_G3xz_S_Fx2y_S_M2_vrr+oned2z*I_ERI_F3x_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Fx2y_S_M2_vrr;
      Double I_ERI_H2x3y_S_Fx2y_S_M1_vrr = PAX*I_ERI_Gx3y_S_Fx2y_S_M1_vrr+WPX*I_ERI_Gx3y_S_Fx2y_S_M2_vrr+oned2z*I_ERI_F3y_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Fx2y_S_M2_vrr+oned2k*I_ERI_Gx3y_S_D2y_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_G2x2y_S_Fx2y_S_M1_vrr+WPZ*I_ERI_G2x2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_H2x3z_S_Fx2y_S_M1_vrr = PAX*I_ERI_Gx3z_S_Fx2y_S_M1_vrr+WPX*I_ERI_Gx3z_S_Fx2y_S_M2_vrr+oned2z*I_ERI_F3z_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Fx2y_S_M2_vrr+oned2k*I_ERI_Gx3z_S_D2y_S_M2_vrr;
      Double I_ERI_Hx4y_S_Fx2y_S_M1_vrr = PAX*I_ERI_G4y_S_Fx2y_S_M1_vrr+WPX*I_ERI_G4y_S_Fx2y_S_M2_vrr+oned2k*I_ERI_G4y_S_D2y_S_M2_vrr;
      Double I_ERI_Hx4z_S_Fx2y_S_M1_vrr = PAX*I_ERI_G4z_S_Fx2y_S_M1_vrr+WPX*I_ERI_G4z_S_Fx2y_S_M2_vrr+oned2k*I_ERI_G4z_S_D2y_S_M2_vrr;
      Double I_ERI_H5y_S_Fx2y_S_M1_vrr = PAY*I_ERI_G4y_S_Fx2y_S_M1_vrr+WPY*I_ERI_G4y_S_Fx2y_S_M2_vrr+4*oned2z*I_ERI_F3y_S_Fx2y_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Dxy_S_M2_vrr;
      Double I_ERI_H4yz_S_Fx2y_S_M1_vrr = PAZ*I_ERI_G4y_S_Fx2y_S_M1_vrr+WPZ*I_ERI_G4y_S_Fx2y_S_M2_vrr;
      Double I_ERI_H3y2z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_G3yz_S_Fx2y_S_M1_vrr+WPZ*I_ERI_G3yz_S_Fx2y_S_M2_vrr+oned2z*I_ERI_F3y_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Fx2y_S_M2_vrr;
      Double I_ERI_H2y3z_S_Fx2y_S_M1_vrr = PAY*I_ERI_Gy3z_S_Fx2y_S_M1_vrr+WPY*I_ERI_Gy3z_S_Fx2y_S_M2_vrr+oned2z*I_ERI_F3z_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_Gy3z_S_Dxy_S_M2_vrr;
      Double I_ERI_Hy4z_S_Fx2y_S_M1_vrr = PAY*I_ERI_G4z_S_Fx2y_S_M1_vrr+WPY*I_ERI_G4z_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Dxy_S_M2_vrr;
      Double I_ERI_H5z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_G4z_S_Fx2y_S_M1_vrr+WPZ*I_ERI_G4z_S_Fx2y_S_M2_vrr+4*oned2z*I_ERI_F3z_S_Fx2y_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_Fx2y_S_M2_vrr;
      Double I_ERI_H5x_S_Fxyz_S_M1_vrr = PAX*I_ERI_G4x_S_Fxyz_S_M1_vrr+WPX*I_ERI_G4x_S_Fxyz_S_M2_vrr+4*oned2z*I_ERI_F3x_S_Fxyz_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G4x_S_Dyz_S_M2_vrr;
      Double I_ERI_H4xy_S_Fxyz_S_M1_vrr = PAY*I_ERI_G4x_S_Fxyz_S_M1_vrr+WPY*I_ERI_G4x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G4x_S_Dxz_S_M2_vrr;
      Double I_ERI_H4xz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_G4x_S_Fxyz_S_M1_vrr+WPZ*I_ERI_G4x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G4x_S_Dxy_S_M2_vrr;
      Double I_ERI_H3x2y_S_Fxyz_S_M1_vrr = PAY*I_ERI_G3xy_S_Fxyz_S_M1_vrr+WPY*I_ERI_G3xy_S_Fxyz_S_M2_vrr+oned2z*I_ERI_F3x_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G3xy_S_Dxz_S_M2_vrr;
      Double I_ERI_H3x2z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_G3xz_S_Fxyz_S_M1_vrr+WPZ*I_ERI_G3xz_S_Fxyz_S_M2_vrr+oned2z*I_ERI_F3x_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G3xz_S_Dxy_S_M2_vrr;
      Double I_ERI_H2x3y_S_Fxyz_S_M1_vrr = PAX*I_ERI_Gx3y_S_Fxyz_S_M1_vrr+WPX*I_ERI_Gx3y_S_Fxyz_S_M2_vrr+oned2z*I_ERI_F3y_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Gx3y_S_Dyz_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_G2x2y_S_Fxyz_S_M1_vrr+WPZ*I_ERI_G2x2y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G2x2y_S_Dxy_S_M2_vrr;
      Double I_ERI_H2x3z_S_Fxyz_S_M1_vrr = PAX*I_ERI_Gx3z_S_Fxyz_S_M1_vrr+WPX*I_ERI_Gx3z_S_Fxyz_S_M2_vrr+oned2z*I_ERI_F3z_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Gx3z_S_Dyz_S_M2_vrr;
      Double I_ERI_Hx4y_S_Fxyz_S_M1_vrr = PAX*I_ERI_G4y_S_Fxyz_S_M1_vrr+WPX*I_ERI_G4y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G4y_S_Dyz_S_M2_vrr;
      Double I_ERI_Hx4z_S_Fxyz_S_M1_vrr = PAX*I_ERI_G4z_S_Fxyz_S_M1_vrr+WPX*I_ERI_G4z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G4z_S_Dyz_S_M2_vrr;
      Double I_ERI_H5y_S_Fxyz_S_M1_vrr = PAY*I_ERI_G4y_S_Fxyz_S_M1_vrr+WPY*I_ERI_G4y_S_Fxyz_S_M2_vrr+4*oned2z*I_ERI_F3y_S_Fxyz_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G4y_S_Dxz_S_M2_vrr;
      Double I_ERI_H4yz_S_Fxyz_S_M1_vrr = PAZ*I_ERI_G4y_S_Fxyz_S_M1_vrr+WPZ*I_ERI_G4y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G4y_S_Dxy_S_M2_vrr;
      Double I_ERI_H3y2z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_G3yz_S_Fxyz_S_M1_vrr+WPZ*I_ERI_G3yz_S_Fxyz_S_M2_vrr+oned2z*I_ERI_F3y_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G3yz_S_Dxy_S_M2_vrr;
      Double I_ERI_H2y3z_S_Fxyz_S_M1_vrr = PAY*I_ERI_Gy3z_S_Fxyz_S_M1_vrr+WPY*I_ERI_Gy3z_S_Fxyz_S_M2_vrr+oned2z*I_ERI_F3z_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Gy3z_S_Dxz_S_M2_vrr;
      Double I_ERI_Hy4z_S_Fxyz_S_M1_vrr = PAY*I_ERI_G4z_S_Fxyz_S_M1_vrr+WPY*I_ERI_G4z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G4z_S_Dxz_S_M2_vrr;
      Double I_ERI_H5z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_G4z_S_Fxyz_S_M1_vrr+WPZ*I_ERI_G4z_S_Fxyz_S_M2_vrr+4*oned2z*I_ERI_F3z_S_Fxyz_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_Fxyz_S_M2_vrr+oned2k*I_ERI_G4z_S_Dxy_S_M2_vrr;
      Double I_ERI_H5x_S_Fx2z_S_M1_vrr = PAX*I_ERI_G4x_S_Fx2z_S_M1_vrr+WPX*I_ERI_G4x_S_Fx2z_S_M2_vrr+4*oned2z*I_ERI_F3x_S_Fx2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_Fx2z_S_M2_vrr+oned2k*I_ERI_G4x_S_D2z_S_M2_vrr;
      Double I_ERI_H4xy_S_Fx2z_S_M1_vrr = PAY*I_ERI_G4x_S_Fx2z_S_M1_vrr+WPY*I_ERI_G4x_S_Fx2z_S_M2_vrr;
      Double I_ERI_H4xz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_G4x_S_Fx2z_S_M1_vrr+WPZ*I_ERI_G4x_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Dxz_S_M2_vrr;
      Double I_ERI_H3x2y_S_Fx2z_S_M1_vrr = PAY*I_ERI_G3xy_S_Fx2z_S_M1_vrr+WPY*I_ERI_G3xy_S_Fx2z_S_M2_vrr+oned2z*I_ERI_F3x_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Fx2z_S_M2_vrr;
      Double I_ERI_H3x2z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_G3xz_S_Fx2z_S_M1_vrr+WPZ*I_ERI_G3xz_S_Fx2z_S_M2_vrr+oned2z*I_ERI_F3x_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_G3xz_S_Dxz_S_M2_vrr;
      Double I_ERI_H2x3y_S_Fx2z_S_M1_vrr = PAX*I_ERI_Gx3y_S_Fx2z_S_M1_vrr+WPX*I_ERI_Gx3y_S_Fx2z_S_M2_vrr+oned2z*I_ERI_F3y_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Fx2z_S_M2_vrr+oned2k*I_ERI_Gx3y_S_D2z_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_G2x2y_S_Fx2z_S_M1_vrr+WPZ*I_ERI_G2x2y_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_G2x2y_S_Dxz_S_M2_vrr;
      Double I_ERI_H2x3z_S_Fx2z_S_M1_vrr = PAX*I_ERI_Gx3z_S_Fx2z_S_M1_vrr+WPX*I_ERI_Gx3z_S_Fx2z_S_M2_vrr+oned2z*I_ERI_F3z_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Fx2z_S_M2_vrr+oned2k*I_ERI_Gx3z_S_D2z_S_M2_vrr;
      Double I_ERI_Hx4y_S_Fx2z_S_M1_vrr = PAX*I_ERI_G4y_S_Fx2z_S_M1_vrr+WPX*I_ERI_G4y_S_Fx2z_S_M2_vrr+oned2k*I_ERI_G4y_S_D2z_S_M2_vrr;
      Double I_ERI_Hx4z_S_Fx2z_S_M1_vrr = PAX*I_ERI_G4z_S_Fx2z_S_M1_vrr+WPX*I_ERI_G4z_S_Fx2z_S_M2_vrr+oned2k*I_ERI_G4z_S_D2z_S_M2_vrr;
      Double I_ERI_H5y_S_Fx2z_S_M1_vrr = PAY*I_ERI_G4y_S_Fx2z_S_M1_vrr+WPY*I_ERI_G4y_S_Fx2z_S_M2_vrr+4*oned2z*I_ERI_F3y_S_Fx2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_Fx2z_S_M2_vrr;
      Double I_ERI_H4yz_S_Fx2z_S_M1_vrr = PAZ*I_ERI_G4y_S_Fx2z_S_M1_vrr+WPZ*I_ERI_G4y_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Dxz_S_M2_vrr;
      Double I_ERI_H3y2z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_G3yz_S_Fx2z_S_M1_vrr+WPZ*I_ERI_G3yz_S_Fx2z_S_M2_vrr+oned2z*I_ERI_F3y_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_G3yz_S_Dxz_S_M2_vrr;
      Double I_ERI_H2y3z_S_Fx2z_S_M1_vrr = PAY*I_ERI_Gy3z_S_Fx2z_S_M1_vrr+WPY*I_ERI_Gy3z_S_Fx2z_S_M2_vrr+oned2z*I_ERI_F3z_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Fx2z_S_M2_vrr;
      Double I_ERI_Hy4z_S_Fx2z_S_M1_vrr = PAY*I_ERI_G4z_S_Fx2z_S_M1_vrr+WPY*I_ERI_G4z_S_Fx2z_S_M2_vrr;
      Double I_ERI_H5z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_G4z_S_Fx2z_S_M1_vrr+WPZ*I_ERI_G4z_S_Fx2z_S_M2_vrr+4*oned2z*I_ERI_F3z_S_Fx2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Dxz_S_M2_vrr;
      Double I_ERI_H5x_S_F3y_S_M1_vrr = PAX*I_ERI_G4x_S_F3y_S_M1_vrr+WPX*I_ERI_G4x_S_F3y_S_M2_vrr+4*oned2z*I_ERI_F3x_S_F3y_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_F3y_S_M2_vrr;
      Double I_ERI_H4xy_S_F3y_S_M1_vrr = PAY*I_ERI_G4x_S_F3y_S_M1_vrr+WPY*I_ERI_G4x_S_F3y_S_M2_vrr+3*oned2k*I_ERI_G4x_S_D2y_S_M2_vrr;
      Double I_ERI_H4xz_S_F3y_S_M1_vrr = PAZ*I_ERI_G4x_S_F3y_S_M1_vrr+WPZ*I_ERI_G4x_S_F3y_S_M2_vrr;
      Double I_ERI_H3x2y_S_F3y_S_M1_vrr = PAY*I_ERI_G3xy_S_F3y_S_M1_vrr+WPY*I_ERI_G3xy_S_F3y_S_M2_vrr+oned2z*I_ERI_F3x_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F3y_S_M2_vrr+3*oned2k*I_ERI_G3xy_S_D2y_S_M2_vrr;
      Double I_ERI_H3x2z_S_F3y_S_M1_vrr = PAZ*I_ERI_G3xz_S_F3y_S_M1_vrr+WPZ*I_ERI_G3xz_S_F3y_S_M2_vrr+oned2z*I_ERI_F3x_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F3y_S_M2_vrr;
      Double I_ERI_H2x3y_S_F3y_S_M1_vrr = PAX*I_ERI_Gx3y_S_F3y_S_M1_vrr+WPX*I_ERI_Gx3y_S_F3y_S_M2_vrr+oned2z*I_ERI_F3y_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F3y_S_M2_vrr;
      Double I_ERI_H2x2yz_S_F3y_S_M1_vrr = PAZ*I_ERI_G2x2y_S_F3y_S_M1_vrr+WPZ*I_ERI_G2x2y_S_F3y_S_M2_vrr;
      Double I_ERI_H2x3z_S_F3y_S_M1_vrr = PAX*I_ERI_Gx3z_S_F3y_S_M1_vrr+WPX*I_ERI_Gx3z_S_F3y_S_M2_vrr+oned2z*I_ERI_F3z_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F3y_S_M2_vrr;
      Double I_ERI_Hx4y_S_F3y_S_M1_vrr = PAX*I_ERI_G4y_S_F3y_S_M1_vrr+WPX*I_ERI_G4y_S_F3y_S_M2_vrr;
      Double I_ERI_Hx4z_S_F3y_S_M1_vrr = PAX*I_ERI_G4z_S_F3y_S_M1_vrr+WPX*I_ERI_G4z_S_F3y_S_M2_vrr;
      Double I_ERI_H5y_S_F3y_S_M1_vrr = PAY*I_ERI_G4y_S_F3y_S_M1_vrr+WPY*I_ERI_G4y_S_F3y_S_M2_vrr+4*oned2z*I_ERI_F3y_S_F3y_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_F3y_S_M2_vrr+3*oned2k*I_ERI_G4y_S_D2y_S_M2_vrr;
      Double I_ERI_H4yz_S_F3y_S_M1_vrr = PAZ*I_ERI_G4y_S_F3y_S_M1_vrr+WPZ*I_ERI_G4y_S_F3y_S_M2_vrr;
      Double I_ERI_H3y2z_S_F3y_S_M1_vrr = PAZ*I_ERI_G3yz_S_F3y_S_M1_vrr+WPZ*I_ERI_G3yz_S_F3y_S_M2_vrr+oned2z*I_ERI_F3y_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F3y_S_M2_vrr;
      Double I_ERI_H2y3z_S_F3y_S_M1_vrr = PAY*I_ERI_Gy3z_S_F3y_S_M1_vrr+WPY*I_ERI_Gy3z_S_F3y_S_M2_vrr+oned2z*I_ERI_F3z_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F3y_S_M2_vrr+3*oned2k*I_ERI_Gy3z_S_D2y_S_M2_vrr;
      Double I_ERI_Hy4z_S_F3y_S_M1_vrr = PAY*I_ERI_G4z_S_F3y_S_M1_vrr+WPY*I_ERI_G4z_S_F3y_S_M2_vrr+3*oned2k*I_ERI_G4z_S_D2y_S_M2_vrr;
      Double I_ERI_H5z_S_F3y_S_M1_vrr = PAZ*I_ERI_G4z_S_F3y_S_M1_vrr+WPZ*I_ERI_G4z_S_F3y_S_M2_vrr+4*oned2z*I_ERI_F3z_S_F3y_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_F3y_S_M2_vrr;
      Double I_ERI_H5x_S_F2yz_S_M1_vrr = PAX*I_ERI_G4x_S_F2yz_S_M1_vrr+WPX*I_ERI_G4x_S_F2yz_S_M2_vrr+4*oned2z*I_ERI_F3x_S_F2yz_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_F2yz_S_M2_vrr;
      Double I_ERI_H4xy_S_F2yz_S_M1_vrr = PAY*I_ERI_G4x_S_F2yz_S_M1_vrr+WPY*I_ERI_G4x_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Dyz_S_M2_vrr;
      Double I_ERI_H4xz_S_F2yz_S_M1_vrr = PAZ*I_ERI_G4x_S_F2yz_S_M1_vrr+WPZ*I_ERI_G4x_S_F2yz_S_M2_vrr+oned2k*I_ERI_G4x_S_D2y_S_M2_vrr;
      Double I_ERI_H3x2y_S_F2yz_S_M1_vrr = PAY*I_ERI_G3xy_S_F2yz_S_M1_vrr+WPY*I_ERI_G3xy_S_F2yz_S_M2_vrr+oned2z*I_ERI_F3x_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_G3xy_S_Dyz_S_M2_vrr;
      Double I_ERI_H3x2z_S_F2yz_S_M1_vrr = PAZ*I_ERI_G3xz_S_F2yz_S_M1_vrr+WPZ*I_ERI_G3xz_S_F2yz_S_M2_vrr+oned2z*I_ERI_F3x_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F2yz_S_M2_vrr+oned2k*I_ERI_G3xz_S_D2y_S_M2_vrr;
      Double I_ERI_H2x3y_S_F2yz_S_M1_vrr = PAX*I_ERI_Gx3y_S_F2yz_S_M1_vrr+WPX*I_ERI_Gx3y_S_F2yz_S_M2_vrr+oned2z*I_ERI_F3y_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F2yz_S_M2_vrr;
      Double I_ERI_H2x2yz_S_F2yz_S_M1_vrr = PAZ*I_ERI_G2x2y_S_F2yz_S_M1_vrr+WPZ*I_ERI_G2x2y_S_F2yz_S_M2_vrr+oned2k*I_ERI_G2x2y_S_D2y_S_M2_vrr;
      Double I_ERI_H2x3z_S_F2yz_S_M1_vrr = PAX*I_ERI_Gx3z_S_F2yz_S_M1_vrr+WPX*I_ERI_Gx3z_S_F2yz_S_M2_vrr+oned2z*I_ERI_F3z_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F2yz_S_M2_vrr;
      Double I_ERI_Hx4y_S_F2yz_S_M1_vrr = PAX*I_ERI_G4y_S_F2yz_S_M1_vrr+WPX*I_ERI_G4y_S_F2yz_S_M2_vrr;
      Double I_ERI_Hx4z_S_F2yz_S_M1_vrr = PAX*I_ERI_G4z_S_F2yz_S_M1_vrr+WPX*I_ERI_G4z_S_F2yz_S_M2_vrr;
      Double I_ERI_H5y_S_F2yz_S_M1_vrr = PAY*I_ERI_G4y_S_F2yz_S_M1_vrr+WPY*I_ERI_G4y_S_F2yz_S_M2_vrr+4*oned2z*I_ERI_F3y_S_F2yz_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Dyz_S_M2_vrr;
      Double I_ERI_H4yz_S_F2yz_S_M1_vrr = PAZ*I_ERI_G4y_S_F2yz_S_M1_vrr+WPZ*I_ERI_G4y_S_F2yz_S_M2_vrr+oned2k*I_ERI_G4y_S_D2y_S_M2_vrr;
      Double I_ERI_H3y2z_S_F2yz_S_M1_vrr = PAZ*I_ERI_G3yz_S_F2yz_S_M1_vrr+WPZ*I_ERI_G3yz_S_F2yz_S_M2_vrr+oned2z*I_ERI_F3y_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F2yz_S_M2_vrr+oned2k*I_ERI_G3yz_S_D2y_S_M2_vrr;
      Double I_ERI_H2y3z_S_F2yz_S_M1_vrr = PAY*I_ERI_Gy3z_S_F2yz_S_M1_vrr+WPY*I_ERI_Gy3z_S_F2yz_S_M2_vrr+oned2z*I_ERI_F3z_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_Gy3z_S_Dyz_S_M2_vrr;
      Double I_ERI_Hy4z_S_F2yz_S_M1_vrr = PAY*I_ERI_G4z_S_F2yz_S_M1_vrr+WPY*I_ERI_G4z_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Dyz_S_M2_vrr;
      Double I_ERI_H5z_S_F2yz_S_M1_vrr = PAZ*I_ERI_G4z_S_F2yz_S_M1_vrr+WPZ*I_ERI_G4z_S_F2yz_S_M2_vrr+4*oned2z*I_ERI_F3z_S_F2yz_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_F2yz_S_M2_vrr+oned2k*I_ERI_G4z_S_D2y_S_M2_vrr;
      Double I_ERI_H5x_S_Fy2z_S_M1_vrr = PAX*I_ERI_G4x_S_Fy2z_S_M1_vrr+WPX*I_ERI_G4x_S_Fy2z_S_M2_vrr+4*oned2z*I_ERI_F3x_S_Fy2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_Fy2z_S_M2_vrr;
      Double I_ERI_H4xy_S_Fy2z_S_M1_vrr = PAY*I_ERI_G4x_S_Fy2z_S_M1_vrr+WPY*I_ERI_G4x_S_Fy2z_S_M2_vrr+oned2k*I_ERI_G4x_S_D2z_S_M2_vrr;
      Double I_ERI_H4xz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_G4x_S_Fy2z_S_M1_vrr+WPZ*I_ERI_G4x_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Dyz_S_M2_vrr;
      Double I_ERI_H3x2y_S_Fy2z_S_M1_vrr = PAY*I_ERI_G3xy_S_Fy2z_S_M1_vrr+WPY*I_ERI_G3xy_S_Fy2z_S_M2_vrr+oned2z*I_ERI_F3x_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Fy2z_S_M2_vrr+oned2k*I_ERI_G3xy_S_D2z_S_M2_vrr;
      Double I_ERI_H3x2z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_G3xz_S_Fy2z_S_M1_vrr+WPZ*I_ERI_G3xz_S_Fy2z_S_M2_vrr+oned2z*I_ERI_F3x_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_G3xz_S_Dyz_S_M2_vrr;
      Double I_ERI_H2x3y_S_Fy2z_S_M1_vrr = PAX*I_ERI_Gx3y_S_Fy2z_S_M1_vrr+WPX*I_ERI_Gx3y_S_Fy2z_S_M2_vrr+oned2z*I_ERI_F3y_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Fy2z_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_G2x2y_S_Fy2z_S_M1_vrr+WPZ*I_ERI_G2x2y_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_G2x2y_S_Dyz_S_M2_vrr;
      Double I_ERI_H2x3z_S_Fy2z_S_M1_vrr = PAX*I_ERI_Gx3z_S_Fy2z_S_M1_vrr+WPX*I_ERI_Gx3z_S_Fy2z_S_M2_vrr+oned2z*I_ERI_F3z_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Fy2z_S_M2_vrr;
      Double I_ERI_Hx4y_S_Fy2z_S_M1_vrr = PAX*I_ERI_G4y_S_Fy2z_S_M1_vrr+WPX*I_ERI_G4y_S_Fy2z_S_M2_vrr;
      Double I_ERI_Hx4z_S_Fy2z_S_M1_vrr = PAX*I_ERI_G4z_S_Fy2z_S_M1_vrr+WPX*I_ERI_G4z_S_Fy2z_S_M2_vrr;
      Double I_ERI_H5y_S_Fy2z_S_M1_vrr = PAY*I_ERI_G4y_S_Fy2z_S_M1_vrr+WPY*I_ERI_G4y_S_Fy2z_S_M2_vrr+4*oned2z*I_ERI_F3y_S_Fy2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_Fy2z_S_M2_vrr+oned2k*I_ERI_G4y_S_D2z_S_M2_vrr;
      Double I_ERI_H4yz_S_Fy2z_S_M1_vrr = PAZ*I_ERI_G4y_S_Fy2z_S_M1_vrr+WPZ*I_ERI_G4y_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Dyz_S_M2_vrr;
      Double I_ERI_H3y2z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_G3yz_S_Fy2z_S_M1_vrr+WPZ*I_ERI_G3yz_S_Fy2z_S_M2_vrr+oned2z*I_ERI_F3y_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_G3yz_S_Dyz_S_M2_vrr;
      Double I_ERI_H2y3z_S_Fy2z_S_M1_vrr = PAY*I_ERI_Gy3z_S_Fy2z_S_M1_vrr+WPY*I_ERI_Gy3z_S_Fy2z_S_M2_vrr+oned2z*I_ERI_F3z_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Fy2z_S_M2_vrr+oned2k*I_ERI_Gy3z_S_D2z_S_M2_vrr;
      Double I_ERI_Hy4z_S_Fy2z_S_M1_vrr = PAY*I_ERI_G4z_S_Fy2z_S_M1_vrr+WPY*I_ERI_G4z_S_Fy2z_S_M2_vrr+oned2k*I_ERI_G4z_S_D2z_S_M2_vrr;
      Double I_ERI_H5z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_G4z_S_Fy2z_S_M1_vrr+WPZ*I_ERI_G4z_S_Fy2z_S_M2_vrr+4*oned2z*I_ERI_F3z_S_Fy2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Dyz_S_M2_vrr;
      Double I_ERI_H5x_S_F3z_S_M1_vrr = PAX*I_ERI_G4x_S_F3z_S_M1_vrr+WPX*I_ERI_G4x_S_F3z_S_M2_vrr+4*oned2z*I_ERI_F3x_S_F3z_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_F3z_S_M2_vrr;
      Double I_ERI_H4xy_S_F3z_S_M1_vrr = PAY*I_ERI_G4x_S_F3z_S_M1_vrr+WPY*I_ERI_G4x_S_F3z_S_M2_vrr;
      Double I_ERI_H4xz_S_F3z_S_M1_vrr = PAZ*I_ERI_G4x_S_F3z_S_M1_vrr+WPZ*I_ERI_G4x_S_F3z_S_M2_vrr+3*oned2k*I_ERI_G4x_S_D2z_S_M2_vrr;
      Double I_ERI_H3x2y_S_F3z_S_M1_vrr = PAY*I_ERI_G3xy_S_F3z_S_M1_vrr+WPY*I_ERI_G3xy_S_F3z_S_M2_vrr+oned2z*I_ERI_F3x_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F3z_S_M2_vrr;
      Double I_ERI_H3x2z_S_F3z_S_M1_vrr = PAZ*I_ERI_G3xz_S_F3z_S_M1_vrr+WPZ*I_ERI_G3xz_S_F3z_S_M2_vrr+oned2z*I_ERI_F3x_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_F3z_S_M2_vrr+3*oned2k*I_ERI_G3xz_S_D2z_S_M2_vrr;
      Double I_ERI_H2x3y_S_F3z_S_M1_vrr = PAX*I_ERI_Gx3y_S_F3z_S_M1_vrr+WPX*I_ERI_Gx3y_S_F3z_S_M2_vrr+oned2z*I_ERI_F3y_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F3z_S_M2_vrr;
      Double I_ERI_H2x2yz_S_F3z_S_M1_vrr = PAZ*I_ERI_G2x2y_S_F3z_S_M1_vrr+WPZ*I_ERI_G2x2y_S_F3z_S_M2_vrr+3*oned2k*I_ERI_G2x2y_S_D2z_S_M2_vrr;
      Double I_ERI_H2x3z_S_F3z_S_M1_vrr = PAX*I_ERI_Gx3z_S_F3z_S_M1_vrr+WPX*I_ERI_Gx3z_S_F3z_S_M2_vrr+oned2z*I_ERI_F3z_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F3z_S_M2_vrr;
      Double I_ERI_Hx4y_S_F3z_S_M1_vrr = PAX*I_ERI_G4y_S_F3z_S_M1_vrr+WPX*I_ERI_G4y_S_F3z_S_M2_vrr;
      Double I_ERI_Hx4z_S_F3z_S_M1_vrr = PAX*I_ERI_G4z_S_F3z_S_M1_vrr+WPX*I_ERI_G4z_S_F3z_S_M2_vrr;
      Double I_ERI_H5y_S_F3z_S_M1_vrr = PAY*I_ERI_G4y_S_F3z_S_M1_vrr+WPY*I_ERI_G4y_S_F3z_S_M2_vrr+4*oned2z*I_ERI_F3y_S_F3z_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_F3z_S_M2_vrr;
      Double I_ERI_H4yz_S_F3z_S_M1_vrr = PAZ*I_ERI_G4y_S_F3z_S_M1_vrr+WPZ*I_ERI_G4y_S_F3z_S_M2_vrr+3*oned2k*I_ERI_G4y_S_D2z_S_M2_vrr;
      Double I_ERI_H3y2z_S_F3z_S_M1_vrr = PAZ*I_ERI_G3yz_S_F3z_S_M1_vrr+WPZ*I_ERI_G3yz_S_F3z_S_M2_vrr+oned2z*I_ERI_F3y_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_F3z_S_M2_vrr+3*oned2k*I_ERI_G3yz_S_D2z_S_M2_vrr;
      Double I_ERI_H2y3z_S_F3z_S_M1_vrr = PAY*I_ERI_Gy3z_S_F3z_S_M1_vrr+WPY*I_ERI_Gy3z_S_F3z_S_M2_vrr+oned2z*I_ERI_F3z_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_F3z_S_M2_vrr;
      Double I_ERI_Hy4z_S_F3z_S_M1_vrr = PAY*I_ERI_G4z_S_F3z_S_M1_vrr+WPY*I_ERI_G4z_S_F3z_S_M2_vrr;
      Double I_ERI_H5z_S_F3z_S_M1_vrr = PAZ*I_ERI_G4z_S_F3z_S_M1_vrr+WPZ*I_ERI_G4z_S_F3z_S_M2_vrr+4*oned2z*I_ERI_F3z_S_F3z_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_F3z_S_M2_vrr+3*oned2k*I_ERI_G4z_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 75 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_F_S_M2
       ************************************************************/
      Double I_ERI_H5x_S_G4x_S_M1_vrr = PAX*I_ERI_G4x_S_G4x_S_M1_vrr+WPX*I_ERI_G4x_S_G4x_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G4x_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G4x_S_M2_vrr+4*oned2k*I_ERI_G4x_S_F3x_S_M2_vrr;
      Double I_ERI_H4xy_S_G4x_S_M1_vrr = PAY*I_ERI_G4x_S_G4x_S_M1_vrr+WPY*I_ERI_G4x_S_G4x_S_M2_vrr;
      Double I_ERI_H4xz_S_G4x_S_M1_vrr = PAZ*I_ERI_G4x_S_G4x_S_M1_vrr+WPZ*I_ERI_G4x_S_G4x_S_M2_vrr;
      Double I_ERI_H3x2y_S_G4x_S_M1_vrr = PAY*I_ERI_G3xy_S_G4x_S_M1_vrr+WPY*I_ERI_G3xy_S_G4x_S_M2_vrr+oned2z*I_ERI_F3x_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G4x_S_M2_vrr;
      Double I_ERI_H3x2z_S_G4x_S_M1_vrr = PAZ*I_ERI_G3xz_S_G4x_S_M1_vrr+WPZ*I_ERI_G3xz_S_G4x_S_M2_vrr+oned2z*I_ERI_F3x_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G4x_S_M2_vrr;
      Double I_ERI_H2x3y_S_G4x_S_M1_vrr = PAX*I_ERI_Gx3y_S_G4x_S_M1_vrr+WPX*I_ERI_Gx3y_S_G4x_S_M2_vrr+oned2z*I_ERI_F3y_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G4x_S_M2_vrr+4*oned2k*I_ERI_Gx3y_S_F3x_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G4x_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G4x_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G4x_S_M2_vrr;
      Double I_ERI_H2x3z_S_G4x_S_M1_vrr = PAX*I_ERI_Gx3z_S_G4x_S_M1_vrr+WPX*I_ERI_Gx3z_S_G4x_S_M2_vrr+oned2z*I_ERI_F3z_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G4x_S_M2_vrr+4*oned2k*I_ERI_Gx3z_S_F3x_S_M2_vrr;
      Double I_ERI_Hx4y_S_G4x_S_M1_vrr = PAX*I_ERI_G4y_S_G4x_S_M1_vrr+WPX*I_ERI_G4y_S_G4x_S_M2_vrr+4*oned2k*I_ERI_G4y_S_F3x_S_M2_vrr;
      Double I_ERI_Hx4z_S_G4x_S_M1_vrr = PAX*I_ERI_G4z_S_G4x_S_M1_vrr+WPX*I_ERI_G4z_S_G4x_S_M2_vrr+4*oned2k*I_ERI_G4z_S_F3x_S_M2_vrr;
      Double I_ERI_H5y_S_G4x_S_M1_vrr = PAY*I_ERI_G4y_S_G4x_S_M1_vrr+WPY*I_ERI_G4y_S_G4x_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G4x_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G4x_S_M2_vrr;
      Double I_ERI_H4yz_S_G4x_S_M1_vrr = PAZ*I_ERI_G4y_S_G4x_S_M1_vrr+WPZ*I_ERI_G4y_S_G4x_S_M2_vrr;
      Double I_ERI_H3y2z_S_G4x_S_M1_vrr = PAZ*I_ERI_G3yz_S_G4x_S_M1_vrr+WPZ*I_ERI_G3yz_S_G4x_S_M2_vrr+oned2z*I_ERI_F3y_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G4x_S_M2_vrr;
      Double I_ERI_H2y3z_S_G4x_S_M1_vrr = PAY*I_ERI_Gy3z_S_G4x_S_M1_vrr+WPY*I_ERI_Gy3z_S_G4x_S_M2_vrr+oned2z*I_ERI_F3z_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G4x_S_M2_vrr;
      Double I_ERI_Hy4z_S_G4x_S_M1_vrr = PAY*I_ERI_G4z_S_G4x_S_M1_vrr+WPY*I_ERI_G4z_S_G4x_S_M2_vrr;
      Double I_ERI_H5z_S_G4x_S_M1_vrr = PAZ*I_ERI_G4z_S_G4x_S_M1_vrr+WPZ*I_ERI_G4z_S_G4x_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G4x_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G4x_S_M2_vrr;
      Double I_ERI_H5x_S_G3xy_S_M1_vrr = PAX*I_ERI_G4x_S_G3xy_S_M1_vrr+WPX*I_ERI_G4x_S_G3xy_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G3xy_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_G4x_S_F2xy_S_M2_vrr;
      Double I_ERI_H4xy_S_G3xy_S_M1_vrr = PAY*I_ERI_G4x_S_G3xy_S_M1_vrr+WPY*I_ERI_G4x_S_G3xy_S_M2_vrr+oned2k*I_ERI_G4x_S_F3x_S_M2_vrr;
      Double I_ERI_H4xz_S_G3xy_S_M1_vrr = PAZ*I_ERI_G4x_S_G3xy_S_M1_vrr+WPZ*I_ERI_G4x_S_G3xy_S_M2_vrr;
      Double I_ERI_H3x2y_S_G3xy_S_M1_vrr = PAY*I_ERI_G3xy_S_G3xy_S_M1_vrr+WPY*I_ERI_G3xy_S_G3xy_S_M2_vrr+oned2z*I_ERI_F3x_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G3xy_S_M2_vrr+oned2k*I_ERI_G3xy_S_F3x_S_M2_vrr;
      Double I_ERI_H3x2z_S_G3xy_S_M1_vrr = PAZ*I_ERI_G3xz_S_G3xy_S_M1_vrr+WPZ*I_ERI_G3xz_S_G3xy_S_M2_vrr+oned2z*I_ERI_F3x_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G3xy_S_M2_vrr;
      Double I_ERI_H2x3y_S_G3xy_S_M1_vrr = PAX*I_ERI_Gx3y_S_G3xy_S_M1_vrr+WPX*I_ERI_Gx3y_S_G3xy_S_M2_vrr+oned2z*I_ERI_F3y_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_Gx3y_S_F2xy_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G3xy_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G3xy_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G3xy_S_M2_vrr;
      Double I_ERI_H2x3z_S_G3xy_S_M1_vrr = PAX*I_ERI_Gx3z_S_G3xy_S_M1_vrr+WPX*I_ERI_Gx3z_S_G3xy_S_M2_vrr+oned2z*I_ERI_F3z_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_Gx3z_S_F2xy_S_M2_vrr;
      Double I_ERI_Hx4y_S_G3xy_S_M1_vrr = PAX*I_ERI_G4y_S_G3xy_S_M1_vrr+WPX*I_ERI_G4y_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_G4y_S_F2xy_S_M2_vrr;
      Double I_ERI_Hx4z_S_G3xy_S_M1_vrr = PAX*I_ERI_G4z_S_G3xy_S_M1_vrr+WPX*I_ERI_G4z_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_G4z_S_F2xy_S_M2_vrr;
      Double I_ERI_H5y_S_G3xy_S_M1_vrr = PAY*I_ERI_G4y_S_G3xy_S_M1_vrr+WPY*I_ERI_G4y_S_G3xy_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G3xy_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G3xy_S_M2_vrr+oned2k*I_ERI_G4y_S_F3x_S_M2_vrr;
      Double I_ERI_H4yz_S_G3xy_S_M1_vrr = PAZ*I_ERI_G4y_S_G3xy_S_M1_vrr+WPZ*I_ERI_G4y_S_G3xy_S_M2_vrr;
      Double I_ERI_H3y2z_S_G3xy_S_M1_vrr = PAZ*I_ERI_G3yz_S_G3xy_S_M1_vrr+WPZ*I_ERI_G3yz_S_G3xy_S_M2_vrr+oned2z*I_ERI_F3y_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G3xy_S_M2_vrr;
      Double I_ERI_H2y3z_S_G3xy_S_M1_vrr = PAY*I_ERI_Gy3z_S_G3xy_S_M1_vrr+WPY*I_ERI_Gy3z_S_G3xy_S_M2_vrr+oned2z*I_ERI_F3z_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G3xy_S_M2_vrr+oned2k*I_ERI_Gy3z_S_F3x_S_M2_vrr;
      Double I_ERI_Hy4z_S_G3xy_S_M1_vrr = PAY*I_ERI_G4z_S_G3xy_S_M1_vrr+WPY*I_ERI_G4z_S_G3xy_S_M2_vrr+oned2k*I_ERI_G4z_S_F3x_S_M2_vrr;
      Double I_ERI_H5z_S_G3xy_S_M1_vrr = PAZ*I_ERI_G4z_S_G3xy_S_M1_vrr+WPZ*I_ERI_G4z_S_G3xy_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G3xy_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G3xy_S_M2_vrr;
      Double I_ERI_H5x_S_G3xz_S_M1_vrr = PAX*I_ERI_G4x_S_G3xz_S_M1_vrr+WPX*I_ERI_G4x_S_G3xz_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G3xz_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_G4x_S_F2xz_S_M2_vrr;
      Double I_ERI_H4xy_S_G3xz_S_M1_vrr = PAY*I_ERI_G4x_S_G3xz_S_M1_vrr+WPY*I_ERI_G4x_S_G3xz_S_M2_vrr;
      Double I_ERI_H4xz_S_G3xz_S_M1_vrr = PAZ*I_ERI_G4x_S_G3xz_S_M1_vrr+WPZ*I_ERI_G4x_S_G3xz_S_M2_vrr+oned2k*I_ERI_G4x_S_F3x_S_M2_vrr;
      Double I_ERI_H3x2y_S_G3xz_S_M1_vrr = PAY*I_ERI_G3xy_S_G3xz_S_M1_vrr+WPY*I_ERI_G3xy_S_G3xz_S_M2_vrr+oned2z*I_ERI_F3x_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G3xz_S_M2_vrr;
      Double I_ERI_H3x2z_S_G3xz_S_M1_vrr = PAZ*I_ERI_G3xz_S_G3xz_S_M1_vrr+WPZ*I_ERI_G3xz_S_G3xz_S_M2_vrr+oned2z*I_ERI_F3x_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G3xz_S_M2_vrr+oned2k*I_ERI_G3xz_S_F3x_S_M2_vrr;
      Double I_ERI_H2x3y_S_G3xz_S_M1_vrr = PAX*I_ERI_Gx3y_S_G3xz_S_M1_vrr+WPX*I_ERI_Gx3y_S_G3xz_S_M2_vrr+oned2z*I_ERI_F3y_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_Gx3y_S_F2xz_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G3xz_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G3xz_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G3xz_S_M2_vrr+oned2k*I_ERI_G2x2y_S_F3x_S_M2_vrr;
      Double I_ERI_H2x3z_S_G3xz_S_M1_vrr = PAX*I_ERI_Gx3z_S_G3xz_S_M1_vrr+WPX*I_ERI_Gx3z_S_G3xz_S_M2_vrr+oned2z*I_ERI_F3z_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_Gx3z_S_F2xz_S_M2_vrr;
      Double I_ERI_Hx4y_S_G3xz_S_M1_vrr = PAX*I_ERI_G4y_S_G3xz_S_M1_vrr+WPX*I_ERI_G4y_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_G4y_S_F2xz_S_M2_vrr;
      Double I_ERI_Hx4z_S_G3xz_S_M1_vrr = PAX*I_ERI_G4z_S_G3xz_S_M1_vrr+WPX*I_ERI_G4z_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_G4z_S_F2xz_S_M2_vrr;
      Double I_ERI_H5y_S_G3xz_S_M1_vrr = PAY*I_ERI_G4y_S_G3xz_S_M1_vrr+WPY*I_ERI_G4y_S_G3xz_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G3xz_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G3xz_S_M2_vrr;
      Double I_ERI_H4yz_S_G3xz_S_M1_vrr = PAZ*I_ERI_G4y_S_G3xz_S_M1_vrr+WPZ*I_ERI_G4y_S_G3xz_S_M2_vrr+oned2k*I_ERI_G4y_S_F3x_S_M2_vrr;
      Double I_ERI_H3y2z_S_G3xz_S_M1_vrr = PAZ*I_ERI_G3yz_S_G3xz_S_M1_vrr+WPZ*I_ERI_G3yz_S_G3xz_S_M2_vrr+oned2z*I_ERI_F3y_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G3xz_S_M2_vrr+oned2k*I_ERI_G3yz_S_F3x_S_M2_vrr;
      Double I_ERI_H2y3z_S_G3xz_S_M1_vrr = PAY*I_ERI_Gy3z_S_G3xz_S_M1_vrr+WPY*I_ERI_Gy3z_S_G3xz_S_M2_vrr+oned2z*I_ERI_F3z_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G3xz_S_M2_vrr;
      Double I_ERI_Hy4z_S_G3xz_S_M1_vrr = PAY*I_ERI_G4z_S_G3xz_S_M1_vrr+WPY*I_ERI_G4z_S_G3xz_S_M2_vrr;
      Double I_ERI_H5z_S_G3xz_S_M1_vrr = PAZ*I_ERI_G4z_S_G3xz_S_M1_vrr+WPZ*I_ERI_G4z_S_G3xz_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G3xz_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G3xz_S_M2_vrr+oned2k*I_ERI_G4z_S_F3x_S_M2_vrr;
      Double I_ERI_H5x_S_G2x2y_S_M1_vrr = PAX*I_ERI_G4x_S_G2x2y_S_M1_vrr+WPX*I_ERI_G4x_S_G2x2y_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G2x2y_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Fx2y_S_M2_vrr;
      Double I_ERI_H4xy_S_G2x2y_S_M1_vrr = PAY*I_ERI_G4x_S_G2x2y_S_M1_vrr+WPY*I_ERI_G4x_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_G4x_S_F2xy_S_M2_vrr;
      Double I_ERI_H4xz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_G4x_S_G2x2y_S_M1_vrr+WPZ*I_ERI_G4x_S_G2x2y_S_M2_vrr;
      Double I_ERI_H3x2y_S_G2x2y_S_M1_vrr = PAY*I_ERI_G3xy_S_G2x2y_S_M1_vrr+WPY*I_ERI_G3xy_S_G2x2y_S_M2_vrr+oned2z*I_ERI_F3x_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_G3xy_S_F2xy_S_M2_vrr;
      Double I_ERI_H3x2z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_G3xz_S_G2x2y_S_M1_vrr+WPZ*I_ERI_G3xz_S_G2x2y_S_M2_vrr+oned2z*I_ERI_F3x_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G2x2y_S_M2_vrr;
      Double I_ERI_H2x3y_S_G2x2y_S_M1_vrr = PAX*I_ERI_Gx3y_S_G2x2y_S_M1_vrr+WPX*I_ERI_Gx3y_S_G2x2y_S_M2_vrr+oned2z*I_ERI_F3y_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Gx3y_S_Fx2y_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G2x2y_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G2x2y_S_M2_vrr;
      Double I_ERI_H2x3z_S_G2x2y_S_M1_vrr = PAX*I_ERI_Gx3z_S_G2x2y_S_M1_vrr+WPX*I_ERI_Gx3z_S_G2x2y_S_M2_vrr+oned2z*I_ERI_F3z_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Gx3z_S_Fx2y_S_M2_vrr;
      Double I_ERI_Hx4y_S_G2x2y_S_M1_vrr = PAX*I_ERI_G4y_S_G2x2y_S_M1_vrr+WPX*I_ERI_G4y_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Fx2y_S_M2_vrr;
      Double I_ERI_Hx4z_S_G2x2y_S_M1_vrr = PAX*I_ERI_G4z_S_G2x2y_S_M1_vrr+WPX*I_ERI_G4z_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Fx2y_S_M2_vrr;
      Double I_ERI_H5y_S_G2x2y_S_M1_vrr = PAY*I_ERI_G4y_S_G2x2y_S_M1_vrr+WPY*I_ERI_G4y_S_G2x2y_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G2x2y_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_G4y_S_F2xy_S_M2_vrr;
      Double I_ERI_H4yz_S_G2x2y_S_M1_vrr = PAZ*I_ERI_G4y_S_G2x2y_S_M1_vrr+WPZ*I_ERI_G4y_S_G2x2y_S_M2_vrr;
      Double I_ERI_H3y2z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_G3yz_S_G2x2y_S_M1_vrr+WPZ*I_ERI_G3yz_S_G2x2y_S_M2_vrr+oned2z*I_ERI_F3y_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G2x2y_S_M2_vrr;
      Double I_ERI_H2y3z_S_G2x2y_S_M1_vrr = PAY*I_ERI_Gy3z_S_G2x2y_S_M1_vrr+WPY*I_ERI_Gy3z_S_G2x2y_S_M2_vrr+oned2z*I_ERI_F3z_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Gy3z_S_F2xy_S_M2_vrr;
      Double I_ERI_Hy4z_S_G2x2y_S_M1_vrr = PAY*I_ERI_G4z_S_G2x2y_S_M1_vrr+WPY*I_ERI_G4z_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_G4z_S_F2xy_S_M2_vrr;
      Double I_ERI_H5z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_G4z_S_G2x2y_S_M1_vrr+WPZ*I_ERI_G4z_S_G2x2y_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G2x2y_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G2x2y_S_M2_vrr;
      Double I_ERI_H5x_S_G2xyz_S_M1_vrr = PAX*I_ERI_G4x_S_G2xyz_S_M1_vrr+WPX*I_ERI_G4x_S_G2xyz_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G2xyz_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Fxyz_S_M2_vrr;
      Double I_ERI_H4xy_S_G2xyz_S_M1_vrr = PAY*I_ERI_G4x_S_G2xyz_S_M1_vrr+WPY*I_ERI_G4x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G4x_S_F2xz_S_M2_vrr;
      Double I_ERI_H4xz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_G4x_S_G2xyz_S_M1_vrr+WPZ*I_ERI_G4x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G4x_S_F2xy_S_M2_vrr;
      Double I_ERI_H3x2y_S_G2xyz_S_M1_vrr = PAY*I_ERI_G3xy_S_G2xyz_S_M1_vrr+WPY*I_ERI_G3xy_S_G2xyz_S_M2_vrr+oned2z*I_ERI_F3x_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G3xy_S_F2xz_S_M2_vrr;
      Double I_ERI_H3x2z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_G3xz_S_G2xyz_S_M1_vrr+WPZ*I_ERI_G3xz_S_G2xyz_S_M2_vrr+oned2z*I_ERI_F3x_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G3xz_S_F2xy_S_M2_vrr;
      Double I_ERI_H2x3y_S_G2xyz_S_M1_vrr = PAX*I_ERI_Gx3y_S_G2xyz_S_M1_vrr+WPX*I_ERI_Gx3y_S_G2xyz_S_M2_vrr+oned2z*I_ERI_F3y_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_Gx3y_S_Fxyz_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G2xyz_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G2x2y_S_F2xy_S_M2_vrr;
      Double I_ERI_H2x3z_S_G2xyz_S_M1_vrr = PAX*I_ERI_Gx3z_S_G2xyz_S_M1_vrr+WPX*I_ERI_Gx3z_S_G2xyz_S_M2_vrr+oned2z*I_ERI_F3z_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_Gx3z_S_Fxyz_S_M2_vrr;
      Double I_ERI_Hx4y_S_G2xyz_S_M1_vrr = PAX*I_ERI_G4y_S_G2xyz_S_M1_vrr+WPX*I_ERI_G4y_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Fxyz_S_M2_vrr;
      Double I_ERI_Hx4z_S_G2xyz_S_M1_vrr = PAX*I_ERI_G4z_S_G2xyz_S_M1_vrr+WPX*I_ERI_G4z_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Fxyz_S_M2_vrr;
      Double I_ERI_H5y_S_G2xyz_S_M1_vrr = PAY*I_ERI_G4y_S_G2xyz_S_M1_vrr+WPY*I_ERI_G4y_S_G2xyz_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G2xyz_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G4y_S_F2xz_S_M2_vrr;
      Double I_ERI_H4yz_S_G2xyz_S_M1_vrr = PAZ*I_ERI_G4y_S_G2xyz_S_M1_vrr+WPZ*I_ERI_G4y_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G4y_S_F2xy_S_M2_vrr;
      Double I_ERI_H3y2z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_G3yz_S_G2xyz_S_M1_vrr+WPZ*I_ERI_G3yz_S_G2xyz_S_M2_vrr+oned2z*I_ERI_F3y_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G3yz_S_F2xy_S_M2_vrr;
      Double I_ERI_H2y3z_S_G2xyz_S_M1_vrr = PAY*I_ERI_Gy3z_S_G2xyz_S_M1_vrr+WPY*I_ERI_Gy3z_S_G2xyz_S_M2_vrr+oned2z*I_ERI_F3z_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G2xyz_S_M2_vrr+oned2k*I_ERI_Gy3z_S_F2xz_S_M2_vrr;
      Double I_ERI_Hy4z_S_G2xyz_S_M1_vrr = PAY*I_ERI_G4z_S_G2xyz_S_M1_vrr+WPY*I_ERI_G4z_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G4z_S_F2xz_S_M2_vrr;
      Double I_ERI_H5z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_G4z_S_G2xyz_S_M1_vrr+WPZ*I_ERI_G4z_S_G2xyz_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G2xyz_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G2xyz_S_M2_vrr+oned2k*I_ERI_G4z_S_F2xy_S_M2_vrr;
      Double I_ERI_H5x_S_G2x2z_S_M1_vrr = PAX*I_ERI_G4x_S_G2x2z_S_M1_vrr+WPX*I_ERI_G4x_S_G2x2z_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G2x2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Fx2z_S_M2_vrr;
      Double I_ERI_H4xy_S_G2x2z_S_M1_vrr = PAY*I_ERI_G4x_S_G2x2z_S_M1_vrr+WPY*I_ERI_G4x_S_G2x2z_S_M2_vrr;
      Double I_ERI_H4xz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_G4x_S_G2x2z_S_M1_vrr+WPZ*I_ERI_G4x_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_G4x_S_F2xz_S_M2_vrr;
      Double I_ERI_H3x2y_S_G2x2z_S_M1_vrr = PAY*I_ERI_G3xy_S_G2x2z_S_M1_vrr+WPY*I_ERI_G3xy_S_G2x2z_S_M2_vrr+oned2z*I_ERI_F3x_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G2x2z_S_M2_vrr;
      Double I_ERI_H3x2z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_G3xz_S_G2x2z_S_M1_vrr+WPZ*I_ERI_G3xz_S_G2x2z_S_M2_vrr+oned2z*I_ERI_F3x_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_G3xz_S_F2xz_S_M2_vrr;
      Double I_ERI_H2x3y_S_G2x2z_S_M1_vrr = PAX*I_ERI_Gx3y_S_G2x2z_S_M1_vrr+WPX*I_ERI_Gx3y_S_G2x2z_S_M2_vrr+oned2z*I_ERI_F3y_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_Gx3y_S_Fx2z_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G2x2z_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_G2x2y_S_F2xz_S_M2_vrr;
      Double I_ERI_H2x3z_S_G2x2z_S_M1_vrr = PAX*I_ERI_Gx3z_S_G2x2z_S_M1_vrr+WPX*I_ERI_Gx3z_S_G2x2z_S_M2_vrr+oned2z*I_ERI_F3z_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_Gx3z_S_Fx2z_S_M2_vrr;
      Double I_ERI_Hx4y_S_G2x2z_S_M1_vrr = PAX*I_ERI_G4y_S_G2x2z_S_M1_vrr+WPX*I_ERI_G4y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Fx2z_S_M2_vrr;
      Double I_ERI_Hx4z_S_G2x2z_S_M1_vrr = PAX*I_ERI_G4z_S_G2x2z_S_M1_vrr+WPX*I_ERI_G4z_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Fx2z_S_M2_vrr;
      Double I_ERI_H5y_S_G2x2z_S_M1_vrr = PAY*I_ERI_G4y_S_G2x2z_S_M1_vrr+WPY*I_ERI_G4y_S_G2x2z_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G2x2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G2x2z_S_M2_vrr;
      Double I_ERI_H4yz_S_G2x2z_S_M1_vrr = PAZ*I_ERI_G4y_S_G2x2z_S_M1_vrr+WPZ*I_ERI_G4y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_G4y_S_F2xz_S_M2_vrr;
      Double I_ERI_H3y2z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_G3yz_S_G2x2z_S_M1_vrr+WPZ*I_ERI_G3yz_S_G2x2z_S_M2_vrr+oned2z*I_ERI_F3y_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_G3yz_S_F2xz_S_M2_vrr;
      Double I_ERI_H2y3z_S_G2x2z_S_M1_vrr = PAY*I_ERI_Gy3z_S_G2x2z_S_M1_vrr+WPY*I_ERI_Gy3z_S_G2x2z_S_M2_vrr+oned2z*I_ERI_F3z_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G2x2z_S_M2_vrr;
      Double I_ERI_Hy4z_S_G2x2z_S_M1_vrr = PAY*I_ERI_G4z_S_G2x2z_S_M1_vrr+WPY*I_ERI_G4z_S_G2x2z_S_M2_vrr;
      Double I_ERI_H5z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_G4z_S_G2x2z_S_M1_vrr+WPZ*I_ERI_G4z_S_G2x2z_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G2x2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_G4z_S_F2xz_S_M2_vrr;
      Double I_ERI_H5x_S_Gx3y_S_M1_vrr = PAX*I_ERI_G4x_S_Gx3y_S_M1_vrr+WPX*I_ERI_G4x_S_Gx3y_S_M2_vrr+4*oned2z*I_ERI_F3x_S_Gx3y_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_Gx3y_S_M2_vrr+oned2k*I_ERI_G4x_S_F3y_S_M2_vrr;
      Double I_ERI_H4xy_S_Gx3y_S_M1_vrr = PAY*I_ERI_G4x_S_Gx3y_S_M1_vrr+WPY*I_ERI_G4x_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_G4x_S_Fx2y_S_M2_vrr;
      Double I_ERI_H4xz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_G4x_S_Gx3y_S_M1_vrr+WPZ*I_ERI_G4x_S_Gx3y_S_M2_vrr;
      Double I_ERI_H3x2y_S_Gx3y_S_M1_vrr = PAY*I_ERI_G3xy_S_Gx3y_S_M1_vrr+WPY*I_ERI_G3xy_S_Gx3y_S_M2_vrr+oned2z*I_ERI_F3x_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_G3xy_S_Fx2y_S_M2_vrr;
      Double I_ERI_H3x2z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_G3xz_S_Gx3y_S_M1_vrr+WPZ*I_ERI_G3xz_S_Gx3y_S_M2_vrr+oned2z*I_ERI_F3x_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gx3y_S_M2_vrr;
      Double I_ERI_H2x3y_S_Gx3y_S_M1_vrr = PAX*I_ERI_Gx3y_S_Gx3y_S_M1_vrr+WPX*I_ERI_Gx3y_S_Gx3y_S_M2_vrr+oned2z*I_ERI_F3y_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gx3y_S_M2_vrr+oned2k*I_ERI_Gx3y_S_F3y_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_G2x2y_S_Gx3y_S_M1_vrr+WPZ*I_ERI_G2x2y_S_Gx3y_S_M2_vrr;
      Double I_ERI_H2x3z_S_Gx3y_S_M1_vrr = PAX*I_ERI_Gx3z_S_Gx3y_S_M1_vrr+WPX*I_ERI_Gx3z_S_Gx3y_S_M2_vrr+oned2z*I_ERI_F3z_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gx3y_S_M2_vrr+oned2k*I_ERI_Gx3z_S_F3y_S_M2_vrr;
      Double I_ERI_Hx4y_S_Gx3y_S_M1_vrr = PAX*I_ERI_G4y_S_Gx3y_S_M1_vrr+WPX*I_ERI_G4y_S_Gx3y_S_M2_vrr+oned2k*I_ERI_G4y_S_F3y_S_M2_vrr;
      Double I_ERI_Hx4z_S_Gx3y_S_M1_vrr = PAX*I_ERI_G4z_S_Gx3y_S_M1_vrr+WPX*I_ERI_G4z_S_Gx3y_S_M2_vrr+oned2k*I_ERI_G4z_S_F3y_S_M2_vrr;
      Double I_ERI_H5y_S_Gx3y_S_M1_vrr = PAY*I_ERI_G4y_S_Gx3y_S_M1_vrr+WPY*I_ERI_G4y_S_Gx3y_S_M2_vrr+4*oned2z*I_ERI_F3y_S_Gx3y_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_G4y_S_Fx2y_S_M2_vrr;
      Double I_ERI_H4yz_S_Gx3y_S_M1_vrr = PAZ*I_ERI_G4y_S_Gx3y_S_M1_vrr+WPZ*I_ERI_G4y_S_Gx3y_S_M2_vrr;
      Double I_ERI_H3y2z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_G3yz_S_Gx3y_S_M1_vrr+WPZ*I_ERI_G3yz_S_Gx3y_S_M2_vrr+oned2z*I_ERI_F3y_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gx3y_S_M2_vrr;
      Double I_ERI_H2y3z_S_Gx3y_S_M1_vrr = PAY*I_ERI_Gy3z_S_Gx3y_S_M1_vrr+WPY*I_ERI_Gy3z_S_Gx3y_S_M2_vrr+oned2z*I_ERI_F3z_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_Gy3z_S_Fx2y_S_M2_vrr;
      Double I_ERI_Hy4z_S_Gx3y_S_M1_vrr = PAY*I_ERI_G4z_S_Gx3y_S_M1_vrr+WPY*I_ERI_G4z_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_G4z_S_Fx2y_S_M2_vrr;
      Double I_ERI_H5z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_G4z_S_Gx3y_S_M1_vrr+WPZ*I_ERI_G4z_S_Gx3y_S_M2_vrr+4*oned2z*I_ERI_F3z_S_Gx3y_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_Gx3y_S_M2_vrr;
      Double I_ERI_H5x_S_Gx2yz_S_M1_vrr = PAX*I_ERI_G4x_S_Gx2yz_S_M1_vrr+WPX*I_ERI_G4x_S_Gx2yz_S_M2_vrr+4*oned2z*I_ERI_F3x_S_Gx2yz_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_G4x_S_F2yz_S_M2_vrr;
      Double I_ERI_H4xy_S_Gx2yz_S_M1_vrr = PAY*I_ERI_G4x_S_Gx2yz_S_M1_vrr+WPY*I_ERI_G4x_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Fxyz_S_M2_vrr;
      Double I_ERI_H4xz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_G4x_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_G4x_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_G4x_S_Fx2y_S_M2_vrr;
      Double I_ERI_H3x2y_S_Gx2yz_S_M1_vrr = PAY*I_ERI_G3xy_S_Gx2yz_S_M1_vrr+WPY*I_ERI_G3xy_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_F3x_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_G3xy_S_Fxyz_S_M2_vrr;
      Double I_ERI_H3x2z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_G3xz_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_G3xz_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_F3x_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_G3xz_S_Fx2y_S_M2_vrr;
      Double I_ERI_H2x3y_S_Gx2yz_S_M1_vrr = PAX*I_ERI_Gx3y_S_Gx2yz_S_M1_vrr+WPX*I_ERI_Gx3y_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_F3y_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_Gx3y_S_F2yz_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_G2x2y_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_G2x2y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_G2x2y_S_Fx2y_S_M2_vrr;
      Double I_ERI_H2x3z_S_Gx2yz_S_M1_vrr = PAX*I_ERI_Gx3z_S_Gx2yz_S_M1_vrr+WPX*I_ERI_Gx3z_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_F3z_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_Gx3z_S_F2yz_S_M2_vrr;
      Double I_ERI_Hx4y_S_Gx2yz_S_M1_vrr = PAX*I_ERI_G4y_S_Gx2yz_S_M1_vrr+WPX*I_ERI_G4y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_G4y_S_F2yz_S_M2_vrr;
      Double I_ERI_Hx4z_S_Gx2yz_S_M1_vrr = PAX*I_ERI_G4z_S_Gx2yz_S_M1_vrr+WPX*I_ERI_G4z_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_G4z_S_F2yz_S_M2_vrr;
      Double I_ERI_H5y_S_Gx2yz_S_M1_vrr = PAY*I_ERI_G4y_S_Gx2yz_S_M1_vrr+WPY*I_ERI_G4y_S_Gx2yz_S_M2_vrr+4*oned2z*I_ERI_F3y_S_Gx2yz_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Fxyz_S_M2_vrr;
      Double I_ERI_H4yz_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_G4y_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_G4y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_G4y_S_Fx2y_S_M2_vrr;
      Double I_ERI_H3y2z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_G3yz_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_G3yz_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_F3y_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_G3yz_S_Fx2y_S_M2_vrr;
      Double I_ERI_H2y3z_S_Gx2yz_S_M1_vrr = PAY*I_ERI_Gy3z_S_Gx2yz_S_M1_vrr+WPY*I_ERI_Gy3z_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_F3z_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_Gy3z_S_Fxyz_S_M2_vrr;
      Double I_ERI_Hy4z_S_Gx2yz_S_M1_vrr = PAY*I_ERI_G4z_S_Gx2yz_S_M1_vrr+WPY*I_ERI_G4z_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Fxyz_S_M2_vrr;
      Double I_ERI_H5z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_G4z_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_G4z_S_Gx2yz_S_M2_vrr+4*oned2z*I_ERI_F3z_S_Gx2yz_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_G4z_S_Fx2y_S_M2_vrr;
      Double I_ERI_H5x_S_Gxy2z_S_M1_vrr = PAX*I_ERI_G4x_S_Gxy2z_S_M1_vrr+WPX*I_ERI_G4x_S_Gxy2z_S_M2_vrr+4*oned2z*I_ERI_F3x_S_Gxy2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_G4x_S_Fy2z_S_M2_vrr;
      Double I_ERI_H4xy_S_Gxy2z_S_M1_vrr = PAY*I_ERI_G4x_S_Gxy2z_S_M1_vrr+WPY*I_ERI_G4x_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_G4x_S_Fx2z_S_M2_vrr;
      Double I_ERI_H4xz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_G4x_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_G4x_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Fxyz_S_M2_vrr;
      Double I_ERI_H3x2y_S_Gxy2z_S_M1_vrr = PAY*I_ERI_G3xy_S_Gxy2z_S_M1_vrr+WPY*I_ERI_G3xy_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_F3x_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_G3xy_S_Fx2z_S_M2_vrr;
      Double I_ERI_H3x2z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_G3xz_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_G3xz_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_F3x_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_G3xz_S_Fxyz_S_M2_vrr;
      Double I_ERI_H2x3y_S_Gxy2z_S_M1_vrr = PAX*I_ERI_Gx3y_S_Gxy2z_S_M1_vrr+WPX*I_ERI_Gx3y_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_F3y_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Gx3y_S_Fy2z_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_G2x2y_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_G2x2y_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_G2x2y_S_Fxyz_S_M2_vrr;
      Double I_ERI_H2x3z_S_Gxy2z_S_M1_vrr = PAX*I_ERI_Gx3z_S_Gxy2z_S_M1_vrr+WPX*I_ERI_Gx3z_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_F3z_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Gx3z_S_Fy2z_S_M2_vrr;
      Double I_ERI_Hx4y_S_Gxy2z_S_M1_vrr = PAX*I_ERI_G4y_S_Gxy2z_S_M1_vrr+WPX*I_ERI_G4y_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_G4y_S_Fy2z_S_M2_vrr;
      Double I_ERI_Hx4z_S_Gxy2z_S_M1_vrr = PAX*I_ERI_G4z_S_Gxy2z_S_M1_vrr+WPX*I_ERI_G4z_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_G4z_S_Fy2z_S_M2_vrr;
      Double I_ERI_H5y_S_Gxy2z_S_M1_vrr = PAY*I_ERI_G4y_S_Gxy2z_S_M1_vrr+WPY*I_ERI_G4y_S_Gxy2z_S_M2_vrr+4*oned2z*I_ERI_F3y_S_Gxy2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_G4y_S_Fx2z_S_M2_vrr;
      Double I_ERI_H4yz_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_G4y_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_G4y_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Fxyz_S_M2_vrr;
      Double I_ERI_H3y2z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_G3yz_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_G3yz_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_F3y_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_G3yz_S_Fxyz_S_M2_vrr;
      Double I_ERI_H2y3z_S_Gxy2z_S_M1_vrr = PAY*I_ERI_Gy3z_S_Gxy2z_S_M1_vrr+WPY*I_ERI_Gy3z_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_F3z_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Gy3z_S_Fx2z_S_M2_vrr;
      Double I_ERI_Hy4z_S_Gxy2z_S_M1_vrr = PAY*I_ERI_G4z_S_Gxy2z_S_M1_vrr+WPY*I_ERI_G4z_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_G4z_S_Fx2z_S_M2_vrr;
      Double I_ERI_H5z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_G4z_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_G4z_S_Gxy2z_S_M2_vrr+4*oned2z*I_ERI_F3z_S_Gxy2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Fxyz_S_M2_vrr;
      Double I_ERI_H5x_S_Gx3z_S_M1_vrr = PAX*I_ERI_G4x_S_Gx3z_S_M1_vrr+WPX*I_ERI_G4x_S_Gx3z_S_M2_vrr+4*oned2z*I_ERI_F3x_S_Gx3z_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_Gx3z_S_M2_vrr+oned2k*I_ERI_G4x_S_F3z_S_M2_vrr;
      Double I_ERI_H4xy_S_Gx3z_S_M1_vrr = PAY*I_ERI_G4x_S_Gx3z_S_M1_vrr+WPY*I_ERI_G4x_S_Gx3z_S_M2_vrr;
      Double I_ERI_H4xz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_G4x_S_Gx3z_S_M1_vrr+WPZ*I_ERI_G4x_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_G4x_S_Fx2z_S_M2_vrr;
      Double I_ERI_H3x2y_S_Gx3z_S_M1_vrr = PAY*I_ERI_G3xy_S_Gx3z_S_M1_vrr+WPY*I_ERI_G3xy_S_Gx3z_S_M2_vrr+oned2z*I_ERI_F3x_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gx3z_S_M2_vrr;
      Double I_ERI_H3x2z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_G3xz_S_Gx3z_S_M1_vrr+WPZ*I_ERI_G3xz_S_Gx3z_S_M2_vrr+oned2z*I_ERI_F3x_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_G3xz_S_Fx2z_S_M2_vrr;
      Double I_ERI_H2x3y_S_Gx3z_S_M1_vrr = PAX*I_ERI_Gx3y_S_Gx3z_S_M1_vrr+WPX*I_ERI_Gx3y_S_Gx3z_S_M2_vrr+oned2z*I_ERI_F3y_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gx3z_S_M2_vrr+oned2k*I_ERI_Gx3y_S_F3z_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_G2x2y_S_Gx3z_S_M1_vrr+WPZ*I_ERI_G2x2y_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_G2x2y_S_Fx2z_S_M2_vrr;
      Double I_ERI_H2x3z_S_Gx3z_S_M1_vrr = PAX*I_ERI_Gx3z_S_Gx3z_S_M1_vrr+WPX*I_ERI_Gx3z_S_Gx3z_S_M2_vrr+oned2z*I_ERI_F3z_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gx3z_S_M2_vrr+oned2k*I_ERI_Gx3z_S_F3z_S_M2_vrr;
      Double I_ERI_Hx4y_S_Gx3z_S_M1_vrr = PAX*I_ERI_G4y_S_Gx3z_S_M1_vrr+WPX*I_ERI_G4y_S_Gx3z_S_M2_vrr+oned2k*I_ERI_G4y_S_F3z_S_M2_vrr;
      Double I_ERI_Hx4z_S_Gx3z_S_M1_vrr = PAX*I_ERI_G4z_S_Gx3z_S_M1_vrr+WPX*I_ERI_G4z_S_Gx3z_S_M2_vrr+oned2k*I_ERI_G4z_S_F3z_S_M2_vrr;
      Double I_ERI_H5y_S_Gx3z_S_M1_vrr = PAY*I_ERI_G4y_S_Gx3z_S_M1_vrr+WPY*I_ERI_G4y_S_Gx3z_S_M2_vrr+4*oned2z*I_ERI_F3y_S_Gx3z_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_Gx3z_S_M2_vrr;
      Double I_ERI_H4yz_S_Gx3z_S_M1_vrr = PAZ*I_ERI_G4y_S_Gx3z_S_M1_vrr+WPZ*I_ERI_G4y_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_G4y_S_Fx2z_S_M2_vrr;
      Double I_ERI_H3y2z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_G3yz_S_Gx3z_S_M1_vrr+WPZ*I_ERI_G3yz_S_Gx3z_S_M2_vrr+oned2z*I_ERI_F3y_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_G3yz_S_Fx2z_S_M2_vrr;
      Double I_ERI_H2y3z_S_Gx3z_S_M1_vrr = PAY*I_ERI_Gy3z_S_Gx3z_S_M1_vrr+WPY*I_ERI_Gy3z_S_Gx3z_S_M2_vrr+oned2z*I_ERI_F3z_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gx3z_S_M2_vrr;
      Double I_ERI_Hy4z_S_Gx3z_S_M1_vrr = PAY*I_ERI_G4z_S_Gx3z_S_M1_vrr+WPY*I_ERI_G4z_S_Gx3z_S_M2_vrr;
      Double I_ERI_H5z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_G4z_S_Gx3z_S_M1_vrr+WPZ*I_ERI_G4z_S_Gx3z_S_M2_vrr+4*oned2z*I_ERI_F3z_S_Gx3z_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_G4z_S_Fx2z_S_M2_vrr;
      Double I_ERI_H5x_S_G4y_S_M1_vrr = PAX*I_ERI_G4x_S_G4y_S_M1_vrr+WPX*I_ERI_G4x_S_G4y_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G4y_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G4y_S_M2_vrr;
      Double I_ERI_H4xy_S_G4y_S_M1_vrr = PAY*I_ERI_G4x_S_G4y_S_M1_vrr+WPY*I_ERI_G4x_S_G4y_S_M2_vrr+4*oned2k*I_ERI_G4x_S_F3y_S_M2_vrr;
      Double I_ERI_H4xz_S_G4y_S_M1_vrr = PAZ*I_ERI_G4x_S_G4y_S_M1_vrr+WPZ*I_ERI_G4x_S_G4y_S_M2_vrr;
      Double I_ERI_H3x2y_S_G4y_S_M1_vrr = PAY*I_ERI_G3xy_S_G4y_S_M1_vrr+WPY*I_ERI_G3xy_S_G4y_S_M2_vrr+oned2z*I_ERI_F3x_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G4y_S_M2_vrr+4*oned2k*I_ERI_G3xy_S_F3y_S_M2_vrr;
      Double I_ERI_H3x2z_S_G4y_S_M1_vrr = PAZ*I_ERI_G3xz_S_G4y_S_M1_vrr+WPZ*I_ERI_G3xz_S_G4y_S_M2_vrr+oned2z*I_ERI_F3x_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G4y_S_M2_vrr;
      Double I_ERI_H2x3y_S_G4y_S_M1_vrr = PAX*I_ERI_Gx3y_S_G4y_S_M1_vrr+WPX*I_ERI_Gx3y_S_G4y_S_M2_vrr+oned2z*I_ERI_F3y_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G4y_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G4y_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G4y_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G4y_S_M2_vrr;
      Double I_ERI_H2x3z_S_G4y_S_M1_vrr = PAX*I_ERI_Gx3z_S_G4y_S_M1_vrr+WPX*I_ERI_Gx3z_S_G4y_S_M2_vrr+oned2z*I_ERI_F3z_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G4y_S_M2_vrr;
      Double I_ERI_Hx4y_S_G4y_S_M1_vrr = PAX*I_ERI_G4y_S_G4y_S_M1_vrr+WPX*I_ERI_G4y_S_G4y_S_M2_vrr;
      Double I_ERI_Hx4z_S_G4y_S_M1_vrr = PAX*I_ERI_G4z_S_G4y_S_M1_vrr+WPX*I_ERI_G4z_S_G4y_S_M2_vrr;
      Double I_ERI_H5y_S_G4y_S_M1_vrr = PAY*I_ERI_G4y_S_G4y_S_M1_vrr+WPY*I_ERI_G4y_S_G4y_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G4y_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G4y_S_M2_vrr+4*oned2k*I_ERI_G4y_S_F3y_S_M2_vrr;
      Double I_ERI_H4yz_S_G4y_S_M1_vrr = PAZ*I_ERI_G4y_S_G4y_S_M1_vrr+WPZ*I_ERI_G4y_S_G4y_S_M2_vrr;
      Double I_ERI_H3y2z_S_G4y_S_M1_vrr = PAZ*I_ERI_G3yz_S_G4y_S_M1_vrr+WPZ*I_ERI_G3yz_S_G4y_S_M2_vrr+oned2z*I_ERI_F3y_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G4y_S_M2_vrr;
      Double I_ERI_H2y3z_S_G4y_S_M1_vrr = PAY*I_ERI_Gy3z_S_G4y_S_M1_vrr+WPY*I_ERI_Gy3z_S_G4y_S_M2_vrr+oned2z*I_ERI_F3z_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G4y_S_M2_vrr+4*oned2k*I_ERI_Gy3z_S_F3y_S_M2_vrr;
      Double I_ERI_Hy4z_S_G4y_S_M1_vrr = PAY*I_ERI_G4z_S_G4y_S_M1_vrr+WPY*I_ERI_G4z_S_G4y_S_M2_vrr+4*oned2k*I_ERI_G4z_S_F3y_S_M2_vrr;
      Double I_ERI_H5z_S_G4y_S_M1_vrr = PAZ*I_ERI_G4z_S_G4y_S_M1_vrr+WPZ*I_ERI_G4z_S_G4y_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G4y_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G4y_S_M2_vrr;
      Double I_ERI_H5x_S_G3yz_S_M1_vrr = PAX*I_ERI_G4x_S_G3yz_S_M1_vrr+WPX*I_ERI_G4x_S_G3yz_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G3yz_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G3yz_S_M2_vrr;
      Double I_ERI_H4xy_S_G3yz_S_M1_vrr = PAY*I_ERI_G4x_S_G3yz_S_M1_vrr+WPY*I_ERI_G4x_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_G4x_S_F2yz_S_M2_vrr;
      Double I_ERI_H4xz_S_G3yz_S_M1_vrr = PAZ*I_ERI_G4x_S_G3yz_S_M1_vrr+WPZ*I_ERI_G4x_S_G3yz_S_M2_vrr+oned2k*I_ERI_G4x_S_F3y_S_M2_vrr;
      Double I_ERI_H3x2y_S_G3yz_S_M1_vrr = PAY*I_ERI_G3xy_S_G3yz_S_M1_vrr+WPY*I_ERI_G3xy_S_G3yz_S_M2_vrr+oned2z*I_ERI_F3x_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_G3xy_S_F2yz_S_M2_vrr;
      Double I_ERI_H3x2z_S_G3yz_S_M1_vrr = PAZ*I_ERI_G3xz_S_G3yz_S_M1_vrr+WPZ*I_ERI_G3xz_S_G3yz_S_M2_vrr+oned2z*I_ERI_F3x_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G3yz_S_M2_vrr+oned2k*I_ERI_G3xz_S_F3y_S_M2_vrr;
      Double I_ERI_H2x3y_S_G3yz_S_M1_vrr = PAX*I_ERI_Gx3y_S_G3yz_S_M1_vrr+WPX*I_ERI_Gx3y_S_G3yz_S_M2_vrr+oned2z*I_ERI_F3y_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G3yz_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G3yz_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G3yz_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G3yz_S_M2_vrr+oned2k*I_ERI_G2x2y_S_F3y_S_M2_vrr;
      Double I_ERI_H2x3z_S_G3yz_S_M1_vrr = PAX*I_ERI_Gx3z_S_G3yz_S_M1_vrr+WPX*I_ERI_Gx3z_S_G3yz_S_M2_vrr+oned2z*I_ERI_F3z_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G3yz_S_M2_vrr;
      Double I_ERI_Hx4y_S_G3yz_S_M1_vrr = PAX*I_ERI_G4y_S_G3yz_S_M1_vrr+WPX*I_ERI_G4y_S_G3yz_S_M2_vrr;
      Double I_ERI_Hx4z_S_G3yz_S_M1_vrr = PAX*I_ERI_G4z_S_G3yz_S_M1_vrr+WPX*I_ERI_G4z_S_G3yz_S_M2_vrr;
      Double I_ERI_H5y_S_G3yz_S_M1_vrr = PAY*I_ERI_G4y_S_G3yz_S_M1_vrr+WPY*I_ERI_G4y_S_G3yz_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G3yz_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_G4y_S_F2yz_S_M2_vrr;
      Double I_ERI_H4yz_S_G3yz_S_M1_vrr = PAZ*I_ERI_G4y_S_G3yz_S_M1_vrr+WPZ*I_ERI_G4y_S_G3yz_S_M2_vrr+oned2k*I_ERI_G4y_S_F3y_S_M2_vrr;
      Double I_ERI_H3y2z_S_G3yz_S_M1_vrr = PAZ*I_ERI_G3yz_S_G3yz_S_M1_vrr+WPZ*I_ERI_G3yz_S_G3yz_S_M2_vrr+oned2z*I_ERI_F3y_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G3yz_S_M2_vrr+oned2k*I_ERI_G3yz_S_F3y_S_M2_vrr;
      Double I_ERI_H2y3z_S_G3yz_S_M1_vrr = PAY*I_ERI_Gy3z_S_G3yz_S_M1_vrr+WPY*I_ERI_Gy3z_S_G3yz_S_M2_vrr+oned2z*I_ERI_F3z_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_Gy3z_S_F2yz_S_M2_vrr;
      Double I_ERI_Hy4z_S_G3yz_S_M1_vrr = PAY*I_ERI_G4z_S_G3yz_S_M1_vrr+WPY*I_ERI_G4z_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_G4z_S_F2yz_S_M2_vrr;
      Double I_ERI_H5z_S_G3yz_S_M1_vrr = PAZ*I_ERI_G4z_S_G3yz_S_M1_vrr+WPZ*I_ERI_G4z_S_G3yz_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G3yz_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G3yz_S_M2_vrr+oned2k*I_ERI_G4z_S_F3y_S_M2_vrr;
      Double I_ERI_H5x_S_G2y2z_S_M1_vrr = PAX*I_ERI_G4x_S_G2y2z_S_M1_vrr+WPX*I_ERI_G4x_S_G2y2z_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G2y2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G2y2z_S_M2_vrr;
      Double I_ERI_H4xy_S_G2y2z_S_M1_vrr = PAY*I_ERI_G4x_S_G2y2z_S_M1_vrr+WPY*I_ERI_G4x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G4x_S_Fy2z_S_M2_vrr;
      Double I_ERI_H4xz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_G4x_S_G2y2z_S_M1_vrr+WPZ*I_ERI_G4x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G4x_S_F2yz_S_M2_vrr;
      Double I_ERI_H3x2y_S_G2y2z_S_M1_vrr = PAY*I_ERI_G3xy_S_G2y2z_S_M1_vrr+WPY*I_ERI_G3xy_S_G2y2z_S_M2_vrr+oned2z*I_ERI_F3x_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G3xy_S_Fy2z_S_M2_vrr;
      Double I_ERI_H3x2z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_G3xz_S_G2y2z_S_M1_vrr+WPZ*I_ERI_G3xz_S_G2y2z_S_M2_vrr+oned2z*I_ERI_F3x_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G3xz_S_F2yz_S_M2_vrr;
      Double I_ERI_H2x3y_S_G2y2z_S_M1_vrr = PAX*I_ERI_Gx3y_S_G2y2z_S_M1_vrr+WPX*I_ERI_Gx3y_S_G2y2z_S_M2_vrr+oned2z*I_ERI_F3y_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G2y2z_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G2y2z_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G2x2y_S_F2yz_S_M2_vrr;
      Double I_ERI_H2x3z_S_G2y2z_S_M1_vrr = PAX*I_ERI_Gx3z_S_G2y2z_S_M1_vrr+WPX*I_ERI_Gx3z_S_G2y2z_S_M2_vrr+oned2z*I_ERI_F3z_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G2y2z_S_M2_vrr;
      Double I_ERI_Hx4y_S_G2y2z_S_M1_vrr = PAX*I_ERI_G4y_S_G2y2z_S_M1_vrr+WPX*I_ERI_G4y_S_G2y2z_S_M2_vrr;
      Double I_ERI_Hx4z_S_G2y2z_S_M1_vrr = PAX*I_ERI_G4z_S_G2y2z_S_M1_vrr+WPX*I_ERI_G4z_S_G2y2z_S_M2_vrr;
      Double I_ERI_H5y_S_G2y2z_S_M1_vrr = PAY*I_ERI_G4y_S_G2y2z_S_M1_vrr+WPY*I_ERI_G4y_S_G2y2z_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G2y2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G4y_S_Fy2z_S_M2_vrr;
      Double I_ERI_H4yz_S_G2y2z_S_M1_vrr = PAZ*I_ERI_G4y_S_G2y2z_S_M1_vrr+WPZ*I_ERI_G4y_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G4y_S_F2yz_S_M2_vrr;
      Double I_ERI_H3y2z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_G3yz_S_G2y2z_S_M1_vrr+WPZ*I_ERI_G3yz_S_G2y2z_S_M2_vrr+oned2z*I_ERI_F3y_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G3yz_S_F2yz_S_M2_vrr;
      Double I_ERI_H2y3z_S_G2y2z_S_M1_vrr = PAY*I_ERI_Gy3z_S_G2y2z_S_M1_vrr+WPY*I_ERI_Gy3z_S_G2y2z_S_M2_vrr+oned2z*I_ERI_F3z_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_Gy3z_S_Fy2z_S_M2_vrr;
      Double I_ERI_Hy4z_S_G2y2z_S_M1_vrr = PAY*I_ERI_G4z_S_G2y2z_S_M1_vrr+WPY*I_ERI_G4z_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G4z_S_Fy2z_S_M2_vrr;
      Double I_ERI_H5z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_G4z_S_G2y2z_S_M1_vrr+WPZ*I_ERI_G4z_S_G2y2z_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G2y2z_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_G4z_S_F2yz_S_M2_vrr;
      Double I_ERI_H5x_S_Gy3z_S_M1_vrr = PAX*I_ERI_G4x_S_Gy3z_S_M1_vrr+WPX*I_ERI_G4x_S_Gy3z_S_M2_vrr+4*oned2z*I_ERI_F3x_S_Gy3z_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_Gy3z_S_M2_vrr;
      Double I_ERI_H4xy_S_Gy3z_S_M1_vrr = PAY*I_ERI_G4x_S_Gy3z_S_M1_vrr+WPY*I_ERI_G4x_S_Gy3z_S_M2_vrr+oned2k*I_ERI_G4x_S_F3z_S_M2_vrr;
      Double I_ERI_H4xz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_G4x_S_Gy3z_S_M1_vrr+WPZ*I_ERI_G4x_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_G4x_S_Fy2z_S_M2_vrr;
      Double I_ERI_H3x2y_S_Gy3z_S_M1_vrr = PAY*I_ERI_G3xy_S_Gy3z_S_M1_vrr+WPY*I_ERI_G3xy_S_Gy3z_S_M2_vrr+oned2z*I_ERI_F3x_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gy3z_S_M2_vrr+oned2k*I_ERI_G3xy_S_F3z_S_M2_vrr;
      Double I_ERI_H3x2z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_G3xz_S_Gy3z_S_M1_vrr+WPZ*I_ERI_G3xz_S_Gy3z_S_M2_vrr+oned2z*I_ERI_F3x_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_G3xz_S_Fy2z_S_M2_vrr;
      Double I_ERI_H2x3y_S_Gy3z_S_M1_vrr = PAX*I_ERI_Gx3y_S_Gy3z_S_M1_vrr+WPX*I_ERI_Gx3y_S_Gy3z_S_M2_vrr+oned2z*I_ERI_F3y_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gy3z_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_G2x2y_S_Gy3z_S_M1_vrr+WPZ*I_ERI_G2x2y_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_G2x2y_S_Fy2z_S_M2_vrr;
      Double I_ERI_H2x3z_S_Gy3z_S_M1_vrr = PAX*I_ERI_Gx3z_S_Gy3z_S_M1_vrr+WPX*I_ERI_Gx3z_S_Gy3z_S_M2_vrr+oned2z*I_ERI_F3z_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gy3z_S_M2_vrr;
      Double I_ERI_Hx4y_S_Gy3z_S_M1_vrr = PAX*I_ERI_G4y_S_Gy3z_S_M1_vrr+WPX*I_ERI_G4y_S_Gy3z_S_M2_vrr;
      Double I_ERI_Hx4z_S_Gy3z_S_M1_vrr = PAX*I_ERI_G4z_S_Gy3z_S_M1_vrr+WPX*I_ERI_G4z_S_Gy3z_S_M2_vrr;
      Double I_ERI_H5y_S_Gy3z_S_M1_vrr = PAY*I_ERI_G4y_S_Gy3z_S_M1_vrr+WPY*I_ERI_G4y_S_Gy3z_S_M2_vrr+4*oned2z*I_ERI_F3y_S_Gy3z_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_Gy3z_S_M2_vrr+oned2k*I_ERI_G4y_S_F3z_S_M2_vrr;
      Double I_ERI_H4yz_S_Gy3z_S_M1_vrr = PAZ*I_ERI_G4y_S_Gy3z_S_M1_vrr+WPZ*I_ERI_G4y_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_G4y_S_Fy2z_S_M2_vrr;
      Double I_ERI_H3y2z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_G3yz_S_Gy3z_S_M1_vrr+WPZ*I_ERI_G3yz_S_Gy3z_S_M2_vrr+oned2z*I_ERI_F3y_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_G3yz_S_Fy2z_S_M2_vrr;
      Double I_ERI_H2y3z_S_Gy3z_S_M1_vrr = PAY*I_ERI_Gy3z_S_Gy3z_S_M1_vrr+WPY*I_ERI_Gy3z_S_Gy3z_S_M2_vrr+oned2z*I_ERI_F3z_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_Gy3z_S_M2_vrr+oned2k*I_ERI_Gy3z_S_F3z_S_M2_vrr;
      Double I_ERI_Hy4z_S_Gy3z_S_M1_vrr = PAY*I_ERI_G4z_S_Gy3z_S_M1_vrr+WPY*I_ERI_G4z_S_Gy3z_S_M2_vrr+oned2k*I_ERI_G4z_S_F3z_S_M2_vrr;
      Double I_ERI_H5z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_G4z_S_Gy3z_S_M1_vrr+WPZ*I_ERI_G4z_S_Gy3z_S_M2_vrr+4*oned2z*I_ERI_F3z_S_Gy3z_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_G4z_S_Fy2z_S_M2_vrr;
      Double I_ERI_H5x_S_G4z_S_M1_vrr = PAX*I_ERI_G4x_S_G4z_S_M1_vrr+WPX*I_ERI_G4x_S_G4z_S_M2_vrr+4*oned2z*I_ERI_F3x_S_G4z_S_M1_vrr-4*rhod2zsq*I_ERI_F3x_S_G4z_S_M2_vrr;
      Double I_ERI_H4xy_S_G4z_S_M1_vrr = PAY*I_ERI_G4x_S_G4z_S_M1_vrr+WPY*I_ERI_G4x_S_G4z_S_M2_vrr;
      Double I_ERI_H4xz_S_G4z_S_M1_vrr = PAZ*I_ERI_G4x_S_G4z_S_M1_vrr+WPZ*I_ERI_G4x_S_G4z_S_M2_vrr+4*oned2k*I_ERI_G4x_S_F3z_S_M2_vrr;
      Double I_ERI_H3x2y_S_G4z_S_M1_vrr = PAY*I_ERI_G3xy_S_G4z_S_M1_vrr+WPY*I_ERI_G3xy_S_G4z_S_M2_vrr+oned2z*I_ERI_F3x_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G4z_S_M2_vrr;
      Double I_ERI_H3x2z_S_G4z_S_M1_vrr = PAZ*I_ERI_G3xz_S_G4z_S_M1_vrr+WPZ*I_ERI_G3xz_S_G4z_S_M2_vrr+oned2z*I_ERI_F3x_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_F3x_S_G4z_S_M2_vrr+4*oned2k*I_ERI_G3xz_S_F3z_S_M2_vrr;
      Double I_ERI_H2x3y_S_G4z_S_M1_vrr = PAX*I_ERI_Gx3y_S_G4z_S_M1_vrr+WPX*I_ERI_Gx3y_S_G4z_S_M2_vrr+oned2z*I_ERI_F3y_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G4z_S_M2_vrr;
      Double I_ERI_H2x2yz_S_G4z_S_M1_vrr = PAZ*I_ERI_G2x2y_S_G4z_S_M1_vrr+WPZ*I_ERI_G2x2y_S_G4z_S_M2_vrr+4*oned2k*I_ERI_G2x2y_S_F3z_S_M2_vrr;
      Double I_ERI_H2x3z_S_G4z_S_M1_vrr = PAX*I_ERI_Gx3z_S_G4z_S_M1_vrr+WPX*I_ERI_Gx3z_S_G4z_S_M2_vrr+oned2z*I_ERI_F3z_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G4z_S_M2_vrr;
      Double I_ERI_Hx4y_S_G4z_S_M1_vrr = PAX*I_ERI_G4y_S_G4z_S_M1_vrr+WPX*I_ERI_G4y_S_G4z_S_M2_vrr;
      Double I_ERI_Hx4z_S_G4z_S_M1_vrr = PAX*I_ERI_G4z_S_G4z_S_M1_vrr+WPX*I_ERI_G4z_S_G4z_S_M2_vrr;
      Double I_ERI_H5y_S_G4z_S_M1_vrr = PAY*I_ERI_G4y_S_G4z_S_M1_vrr+WPY*I_ERI_G4y_S_G4z_S_M2_vrr+4*oned2z*I_ERI_F3y_S_G4z_S_M1_vrr-4*rhod2zsq*I_ERI_F3y_S_G4z_S_M2_vrr;
      Double I_ERI_H4yz_S_G4z_S_M1_vrr = PAZ*I_ERI_G4y_S_G4z_S_M1_vrr+WPZ*I_ERI_G4y_S_G4z_S_M2_vrr+4*oned2k*I_ERI_G4y_S_F3z_S_M2_vrr;
      Double I_ERI_H3y2z_S_G4z_S_M1_vrr = PAZ*I_ERI_G3yz_S_G4z_S_M1_vrr+WPZ*I_ERI_G3yz_S_G4z_S_M2_vrr+oned2z*I_ERI_F3y_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_F3y_S_G4z_S_M2_vrr+4*oned2k*I_ERI_G3yz_S_F3z_S_M2_vrr;
      Double I_ERI_H2y3z_S_G4z_S_M1_vrr = PAY*I_ERI_Gy3z_S_G4z_S_M1_vrr+WPY*I_ERI_Gy3z_S_G4z_S_M2_vrr+oned2z*I_ERI_F3z_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_F3z_S_G4z_S_M2_vrr;
      Double I_ERI_Hy4z_S_G4z_S_M1_vrr = PAY*I_ERI_G4z_S_G4z_S_M1_vrr+WPY*I_ERI_G4z_S_G4z_S_M2_vrr;
      Double I_ERI_H5z_S_G4z_S_M1_vrr = PAZ*I_ERI_G4z_S_G4z_S_M1_vrr+WPZ*I_ERI_G4z_S_G4z_S_M2_vrr+4*oned2z*I_ERI_F3z_S_G4z_S_M1_vrr-4*rhod2zsq*I_ERI_F3z_S_G4z_S_M2_vrr+4*oned2k*I_ERI_G4z_S_F3z_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_S_D2x_S_vrr = QCX*I_ERI_S_S_Px_S_vrr+WQX*I_ERI_S_S_Px_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_D2y_S_vrr = QCY*I_ERI_S_S_Py_S_vrr+WQY*I_ERI_S_S_Py_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_D2z_S_vrr = QCZ*I_ERI_S_S_Pz_S_vrr+WQZ*I_ERI_S_S_Pz_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_S_S_F3x_S_vrr = QCX*I_ERI_S_S_D2x_S_vrr+WQX*I_ERI_S_S_D2x_S_M1_vrr+2*oned2e*I_ERI_S_S_Px_S_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_F2xy_S_vrr = QCY*I_ERI_S_S_D2x_S_vrr+WQY*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_F2xz_S_vrr = QCZ*I_ERI_S_S_D2x_S_vrr+WQZ*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_Fx2y_S_vrr = QCX*I_ERI_S_S_D2y_S_vrr+WQX*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_Fx2z_S_vrr = QCX*I_ERI_S_S_D2z_S_vrr+WQX*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_S_F3y_S_vrr = QCY*I_ERI_S_S_D2y_S_vrr+WQY*I_ERI_S_S_D2y_S_M1_vrr+2*oned2e*I_ERI_S_S_Py_S_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_S_F2yz_S_vrr = QCZ*I_ERI_S_S_D2y_S_vrr+WQZ*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_F3z_S_vrr = QCZ*I_ERI_S_S_D2z_S_vrr+WQZ*I_ERI_S_S_D2z_S_M1_vrr+2*oned2e*I_ERI_S_S_Pz_S_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       ************************************************************/
      Double I_ERI_S_S_G4x_S_vrr = QCX*I_ERI_S_S_F3x_S_vrr+WQX*I_ERI_S_S_F3x_S_M1_vrr+3*oned2e*I_ERI_S_S_D2x_S_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_G3xy_S_vrr = QCY*I_ERI_S_S_F3x_S_vrr+WQY*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_S_S_G3xz_S_vrr = QCZ*I_ERI_S_S_F3x_S_vrr+WQZ*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_S_S_G2x2y_S_vrr = QCY*I_ERI_S_S_F2xy_S_vrr+WQY*I_ERI_S_S_F2xy_S_M1_vrr+oned2e*I_ERI_S_S_D2x_S_vrr-rhod2esq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_G2xyz_S_vrr = QCZ*I_ERI_S_S_F2xy_S_vrr+WQZ*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_S_S_G2x2z_S_vrr = QCZ*I_ERI_S_S_F2xz_S_vrr+WQZ*I_ERI_S_S_F2xz_S_M1_vrr+oned2e*I_ERI_S_S_D2x_S_vrr-rhod2esq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_Gx3y_S_vrr = QCX*I_ERI_S_S_F3y_S_vrr+WQX*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_S_S_Gx2yz_S_vrr = QCZ*I_ERI_S_S_Fx2y_S_vrr+WQZ*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_S_S_Gxy2z_S_vrr = QCY*I_ERI_S_S_Fx2z_S_vrr+WQY*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_S_S_Gx3z_S_vrr = QCX*I_ERI_S_S_F3z_S_vrr+WQX*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_S_S_G4y_S_vrr = QCY*I_ERI_S_S_F3y_S_vrr+WQY*I_ERI_S_S_F3y_S_M1_vrr+3*oned2e*I_ERI_S_S_D2y_S_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_G3yz_S_vrr = QCZ*I_ERI_S_S_F3y_S_vrr+WQZ*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_S_S_G2y2z_S_vrr = QCZ*I_ERI_S_S_F2yz_S_vrr+WQZ*I_ERI_S_S_F2yz_S_M1_vrr+oned2e*I_ERI_S_S_D2y_S_vrr-rhod2esq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_Gy3z_S_vrr = QCY*I_ERI_S_S_F3z_S_vrr+WQY*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_S_S_G4z_S_vrr = QCZ*I_ERI_S_S_F3z_S_vrr+WQZ*I_ERI_S_S_F3z_S_M1_vrr+3*oned2e*I_ERI_S_S_D2z_S_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_G_S
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       ************************************************************/
      Double I_ERI_Px_S_G4x_S_vrr = PAX*I_ERI_S_S_G4x_S_vrr+WPX*I_ERI_S_S_G4x_S_M1_vrr+4*oned2k*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Py_S_G4x_S_vrr = PAY*I_ERI_S_S_G4x_S_vrr+WPY*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_Pz_S_G4x_S_vrr = PAZ*I_ERI_S_S_G4x_S_vrr+WPZ*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_Px_S_G3xy_S_vrr = PAX*I_ERI_S_S_G3xy_S_vrr+WPX*I_ERI_S_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_Py_S_G3xy_S_vrr = PAY*I_ERI_S_S_G3xy_S_vrr+WPY*I_ERI_S_S_G3xy_S_M1_vrr+oned2k*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Pz_S_G3xy_S_vrr = PAZ*I_ERI_S_S_G3xy_S_vrr+WPZ*I_ERI_S_S_G3xy_S_M1_vrr;
      Double I_ERI_Px_S_G3xz_S_vrr = PAX*I_ERI_S_S_G3xz_S_vrr+WPX*I_ERI_S_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Py_S_G3xz_S_vrr = PAY*I_ERI_S_S_G3xz_S_vrr+WPY*I_ERI_S_S_G3xz_S_M1_vrr;
      Double I_ERI_Pz_S_G3xz_S_vrr = PAZ*I_ERI_S_S_G3xz_S_vrr+WPZ*I_ERI_S_S_G3xz_S_M1_vrr+oned2k*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Px_S_G2x2y_S_vrr = PAX*I_ERI_S_S_G2x2y_S_vrr+WPX*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_Py_S_G2x2y_S_vrr = PAY*I_ERI_S_S_G2x2y_S_vrr+WPY*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_Pz_S_G2x2y_S_vrr = PAZ*I_ERI_S_S_G2x2y_S_vrr+WPZ*I_ERI_S_S_G2x2y_S_M1_vrr;
      Double I_ERI_Px_S_G2xyz_S_vrr = PAX*I_ERI_S_S_G2xyz_S_vrr+WPX*I_ERI_S_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M1_vrr;
      Double I_ERI_Py_S_G2xyz_S_vrr = PAY*I_ERI_S_S_G2xyz_S_vrr+WPY*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Pz_S_G2xyz_S_vrr = PAZ*I_ERI_S_S_G2xyz_S_vrr+WPZ*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_Px_S_G2x2z_S_vrr = PAX*I_ERI_S_S_G2x2z_S_vrr+WPX*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Py_S_G2x2z_S_vrr = PAY*I_ERI_S_S_G2x2z_S_vrr+WPY*I_ERI_S_S_G2x2z_S_M1_vrr;
      Double I_ERI_Pz_S_G2x2z_S_vrr = PAZ*I_ERI_S_S_G2x2z_S_vrr+WPZ*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Px_S_Gx3y_S_vrr = PAX*I_ERI_S_S_Gx3y_S_vrr+WPX*I_ERI_S_S_Gx3y_S_M1_vrr+oned2k*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Py_S_Gx3y_S_vrr = PAY*I_ERI_S_S_Gx3y_S_vrr+WPY*I_ERI_S_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_Pz_S_Gx3y_S_vrr = PAZ*I_ERI_S_S_Gx3y_S_vrr+WPZ*I_ERI_S_S_Gx3y_S_M1_vrr;
      Double I_ERI_Px_S_Gx2yz_S_vrr = PAX*I_ERI_S_S_Gx2yz_S_vrr+WPX*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Py_S_Gx2yz_S_vrr = PAY*I_ERI_S_S_Gx2yz_S_vrr+WPY*I_ERI_S_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M1_vrr;
      Double I_ERI_Pz_S_Gx2yz_S_vrr = PAZ*I_ERI_S_S_Gx2yz_S_vrr+WPZ*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_Px_S_Gxy2z_S_vrr = PAX*I_ERI_S_S_Gxy2z_S_vrr+WPX*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Py_S_Gxy2z_S_vrr = PAY*I_ERI_S_S_Gxy2z_S_vrr+WPY*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Pz_S_Gxy2z_S_vrr = PAZ*I_ERI_S_S_Gxy2z_S_vrr+WPZ*I_ERI_S_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M1_vrr;
      Double I_ERI_Px_S_Gx3z_S_vrr = PAX*I_ERI_S_S_Gx3z_S_vrr+WPX*I_ERI_S_S_Gx3z_S_M1_vrr+oned2k*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Py_S_Gx3z_S_vrr = PAY*I_ERI_S_S_Gx3z_S_vrr+WPY*I_ERI_S_S_Gx3z_S_M1_vrr;
      Double I_ERI_Pz_S_Gx3z_S_vrr = PAZ*I_ERI_S_S_Gx3z_S_vrr+WPZ*I_ERI_S_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Px_S_G4y_S_vrr = PAX*I_ERI_S_S_G4y_S_vrr+WPX*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_Py_S_G4y_S_vrr = PAY*I_ERI_S_S_G4y_S_vrr+WPY*I_ERI_S_S_G4y_S_M1_vrr+4*oned2k*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Pz_S_G4y_S_vrr = PAZ*I_ERI_S_S_G4y_S_vrr+WPZ*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_Px_S_G3yz_S_vrr = PAX*I_ERI_S_S_G3yz_S_vrr+WPX*I_ERI_S_S_G3yz_S_M1_vrr;
      Double I_ERI_Py_S_G3yz_S_vrr = PAY*I_ERI_S_S_G3yz_S_vrr+WPY*I_ERI_S_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Pz_S_G3yz_S_vrr = PAZ*I_ERI_S_S_G3yz_S_vrr+WPZ*I_ERI_S_S_G3yz_S_M1_vrr+oned2k*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Px_S_G2y2z_S_vrr = PAX*I_ERI_S_S_G2y2z_S_vrr+WPX*I_ERI_S_S_G2y2z_S_M1_vrr;
      Double I_ERI_Py_S_G2y2z_S_vrr = PAY*I_ERI_S_S_G2y2z_S_vrr+WPY*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Pz_S_G2y2z_S_vrr = PAZ*I_ERI_S_S_G2y2z_S_vrr+WPZ*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Px_S_Gy3z_S_vrr = PAX*I_ERI_S_S_Gy3z_S_vrr+WPX*I_ERI_S_S_Gy3z_S_M1_vrr;
      Double I_ERI_Py_S_Gy3z_S_vrr = PAY*I_ERI_S_S_Gy3z_S_vrr+WPY*I_ERI_S_S_Gy3z_S_M1_vrr+oned2k*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Pz_S_Gy3z_S_vrr = PAZ*I_ERI_S_S_Gy3z_S_vrr+WPZ*I_ERI_S_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Px_S_G4z_S_vrr = PAX*I_ERI_S_S_G4z_S_vrr+WPX*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_Py_S_G4z_S_vrr = PAY*I_ERI_S_S_G4z_S_vrr+WPY*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_Pz_S_G4z_S_vrr = PAZ*I_ERI_S_S_G4z_S_vrr+WPZ*I_ERI_S_S_G4z_S_M1_vrr+4*oned2k*I_ERI_S_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 45 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_G_S
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_G_S
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_vrr = PAX*I_ERI_Px_S_G4x_S_vrr+WPX*I_ERI_Px_S_G4x_S_M1_vrr+oned2z*I_ERI_S_S_G4x_S_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M1_vrr+4*oned2k*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_D2y_S_G4x_S_vrr = PAY*I_ERI_Py_S_G4x_S_vrr+WPY*I_ERI_Py_S_G4x_S_M1_vrr+oned2z*I_ERI_S_S_G4x_S_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_D2z_S_G4x_S_vrr = PAZ*I_ERI_Pz_S_G4x_S_vrr+WPZ*I_ERI_Pz_S_G4x_S_M1_vrr+oned2z*I_ERI_S_S_G4x_S_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_D2x_S_G3xy_S_vrr = PAX*I_ERI_Px_S_G3xy_S_vrr+WPX*I_ERI_Px_S_G3xy_S_M1_vrr+oned2z*I_ERI_S_S_G3xy_S_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_Px_S_F2xy_S_M1_vrr;
      Double I_ERI_D2y_S_G3xy_S_vrr = PAY*I_ERI_Py_S_G3xy_S_vrr+WPY*I_ERI_Py_S_G3xy_S_M1_vrr+oned2z*I_ERI_S_S_G3xy_S_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M1_vrr+oned2k*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_D2z_S_G3xy_S_vrr = PAZ*I_ERI_Pz_S_G3xy_S_vrr+WPZ*I_ERI_Pz_S_G3xy_S_M1_vrr+oned2z*I_ERI_S_S_G3xy_S_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M1_vrr;
      Double I_ERI_D2x_S_G3xz_S_vrr = PAX*I_ERI_Px_S_G3xz_S_vrr+WPX*I_ERI_Px_S_G3xz_S_M1_vrr+oned2z*I_ERI_S_S_G3xz_S_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_Px_S_F2xz_S_M1_vrr;
      Double I_ERI_D2y_S_G3xz_S_vrr = PAY*I_ERI_Py_S_G3xz_S_vrr+WPY*I_ERI_Py_S_G3xz_S_M1_vrr+oned2z*I_ERI_S_S_G3xz_S_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M1_vrr;
      Double I_ERI_D2z_S_G3xz_S_vrr = PAZ*I_ERI_Pz_S_G3xz_S_vrr+WPZ*I_ERI_Pz_S_G3xz_S_M1_vrr+oned2z*I_ERI_S_S_G3xz_S_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M1_vrr+oned2k*I_ERI_Pz_S_F3x_S_M1_vrr;
      Double I_ERI_D2x_S_G2x2y_S_vrr = PAX*I_ERI_Px_S_G2x2y_S_vrr+WPX*I_ERI_Px_S_G2x2y_S_M1_vrr+oned2z*I_ERI_S_S_G2x2y_S_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2y_S_G2x2y_S_vrr = PAY*I_ERI_Py_S_G2x2y_S_vrr+WPY*I_ERI_Py_S_G2x2y_S_M1_vrr+oned2z*I_ERI_S_S_G2x2y_S_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M1_vrr;
      Double I_ERI_D2z_S_G2x2y_S_vrr = PAZ*I_ERI_Pz_S_G2x2y_S_vrr+WPZ*I_ERI_Pz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_S_S_G2x2y_S_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M1_vrr;
      Double I_ERI_D2x_S_G2xyz_S_vrr = PAX*I_ERI_Px_S_G2xyz_S_vrr+WPX*I_ERI_Px_S_G2xyz_S_M1_vrr+oned2z*I_ERI_S_S_G2xyz_S_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M1_vrr;
      Double I_ERI_D2y_S_G2xyz_S_vrr = PAY*I_ERI_Py_S_G2xyz_S_vrr+WPY*I_ERI_Py_S_G2xyz_S_M1_vrr+oned2z*I_ERI_S_S_G2xyz_S_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Py_S_F2xz_S_M1_vrr;
      Double I_ERI_D2z_S_G2xyz_S_vrr = PAZ*I_ERI_Pz_S_G2xyz_S_vrr+WPZ*I_ERI_Pz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_S_S_G2xyz_S_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Pz_S_F2xy_S_M1_vrr;
      Double I_ERI_D2x_S_G2x2z_S_vrr = PAX*I_ERI_Px_S_G2x2z_S_vrr+WPX*I_ERI_Px_S_G2x2z_S_M1_vrr+oned2z*I_ERI_S_S_G2x2z_S_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2y_S_G2x2z_S_vrr = PAY*I_ERI_Py_S_G2x2z_S_vrr+WPY*I_ERI_Py_S_G2x2z_S_M1_vrr+oned2z*I_ERI_S_S_G2x2z_S_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M1_vrr;
      Double I_ERI_D2z_S_G2x2z_S_vrr = PAZ*I_ERI_Pz_S_G2x2z_S_vrr+WPZ*I_ERI_Pz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_S_S_G2x2z_S_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M1_vrr;
      Double I_ERI_D2x_S_Gx3y_S_vrr = PAX*I_ERI_Px_S_Gx3y_S_vrr+WPX*I_ERI_Px_S_Gx3y_S_M1_vrr+oned2z*I_ERI_S_S_Gx3y_S_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M1_vrr+oned2k*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_D2y_S_Gx3y_S_vrr = PAY*I_ERI_Py_S_Gx3y_S_vrr+WPY*I_ERI_Py_S_Gx3y_S_M1_vrr+oned2z*I_ERI_S_S_Gx3y_S_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Py_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2z_S_Gx3y_S_vrr = PAZ*I_ERI_Pz_S_Gx3y_S_vrr+WPZ*I_ERI_Pz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_S_S_Gx3y_S_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M1_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_vrr = PAX*I_ERI_Px_S_Gx2yz_S_vrr+WPX*I_ERI_Px_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_S_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Px_S_F2yz_S_M1_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_vrr = PAY*I_ERI_Py_S_Gx2yz_S_vrr+WPY*I_ERI_Py_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_S_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Py_S_Fxyz_S_M1_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_vrr = PAZ*I_ERI_Pz_S_Gx2yz_S_vrr+WPZ*I_ERI_Pz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_S_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Pz_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_vrr = PAX*I_ERI_Px_S_Gxy2z_S_vrr+WPX*I_ERI_Px_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_S_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Px_S_Fy2z_S_M1_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_vrr = PAY*I_ERI_Py_S_Gxy2z_S_vrr+WPY*I_ERI_Py_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_S_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Py_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_vrr = PAZ*I_ERI_Pz_S_Gxy2z_S_vrr+WPZ*I_ERI_Pz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_S_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Fxyz_S_M1_vrr;
      Double I_ERI_D2x_S_Gx3z_S_vrr = PAX*I_ERI_Px_S_Gx3z_S_vrr+WPX*I_ERI_Px_S_Gx3z_S_M1_vrr+oned2z*I_ERI_S_S_Gx3z_S_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M1_vrr+oned2k*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_D2y_S_Gx3z_S_vrr = PAY*I_ERI_Py_S_Gx3z_S_vrr+WPY*I_ERI_Py_S_Gx3z_S_M1_vrr+oned2z*I_ERI_S_S_Gx3z_S_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M1_vrr;
      Double I_ERI_D2z_S_Gx3z_S_vrr = PAZ*I_ERI_Pz_S_Gx3z_S_vrr+WPZ*I_ERI_Pz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_S_S_Gx3z_S_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Pz_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2x_S_G4y_S_vrr = PAX*I_ERI_Px_S_G4y_S_vrr+WPX*I_ERI_Px_S_G4y_S_M1_vrr+oned2z*I_ERI_S_S_G4y_S_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_D2y_S_G4y_S_vrr = PAY*I_ERI_Py_S_G4y_S_vrr+WPY*I_ERI_Py_S_G4y_S_M1_vrr+oned2z*I_ERI_S_S_G4y_S_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Py_S_F3y_S_M1_vrr;
      Double I_ERI_D2z_S_G4y_S_vrr = PAZ*I_ERI_Pz_S_G4y_S_vrr+WPZ*I_ERI_Pz_S_G4y_S_M1_vrr+oned2z*I_ERI_S_S_G4y_S_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_D2x_S_G3yz_S_vrr = PAX*I_ERI_Px_S_G3yz_S_vrr+WPX*I_ERI_Px_S_G3yz_S_M1_vrr+oned2z*I_ERI_S_S_G3yz_S_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M1_vrr;
      Double I_ERI_D2y_S_G3yz_S_vrr = PAY*I_ERI_Py_S_G3yz_S_vrr+WPY*I_ERI_Py_S_G3yz_S_M1_vrr+oned2z*I_ERI_S_S_G3yz_S_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Py_S_F2yz_S_M1_vrr;
      Double I_ERI_D2z_S_G3yz_S_vrr = PAZ*I_ERI_Pz_S_G3yz_S_vrr+WPZ*I_ERI_Pz_S_G3yz_S_M1_vrr+oned2z*I_ERI_S_S_G3yz_S_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M1_vrr+oned2k*I_ERI_Pz_S_F3y_S_M1_vrr;
      Double I_ERI_D2x_S_G2y2z_S_vrr = PAX*I_ERI_Px_S_G2y2z_S_vrr+WPX*I_ERI_Px_S_G2y2z_S_M1_vrr+oned2z*I_ERI_S_S_G2y2z_S_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M1_vrr;
      Double I_ERI_D2y_S_G2y2z_S_vrr = PAY*I_ERI_Py_S_G2y2z_S_vrr+WPY*I_ERI_Py_S_G2y2z_S_M1_vrr+oned2z*I_ERI_S_S_G2y2z_S_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Fy2z_S_M1_vrr;
      Double I_ERI_D2z_S_G2y2z_S_vrr = PAZ*I_ERI_Pz_S_G2y2z_S_vrr+WPZ*I_ERI_Pz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_S_S_G2y2z_S_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M1_vrr;
      Double I_ERI_D2x_S_Gy3z_S_vrr = PAX*I_ERI_Px_S_Gy3z_S_vrr+WPX*I_ERI_Px_S_Gy3z_S_M1_vrr+oned2z*I_ERI_S_S_Gy3z_S_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M1_vrr;
      Double I_ERI_D2y_S_Gy3z_S_vrr = PAY*I_ERI_Py_S_Gy3z_S_vrr+WPY*I_ERI_Py_S_Gy3z_S_M1_vrr+oned2z*I_ERI_S_S_Gy3z_S_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Py_S_F3z_S_M1_vrr;
      Double I_ERI_D2z_S_Gy3z_S_vrr = PAZ*I_ERI_Pz_S_Gy3z_S_vrr+WPZ*I_ERI_Pz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_S_S_Gy3z_S_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Pz_S_Fy2z_S_M1_vrr;
      Double I_ERI_D2x_S_G4z_S_vrr = PAX*I_ERI_Px_S_G4z_S_vrr+WPX*I_ERI_Px_S_G4z_S_M1_vrr+oned2z*I_ERI_S_S_G4z_S_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_D2y_S_G4z_S_vrr = PAY*I_ERI_Py_S_G4z_S_vrr+WPY*I_ERI_Py_S_G4z_S_M1_vrr+oned2z*I_ERI_S_S_G4z_S_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_D2z_S_G4z_S_vrr = PAZ*I_ERI_Pz_S_G4z_S_vrr+WPZ*I_ERI_Pz_S_G4z_S_M1_vrr+oned2z*I_ERI_S_S_G4z_S_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Pz_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_G_S
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_G_S
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_G4x_S_vrr = PAX*I_ERI_D2x_S_G4x_S_vrr+WPX*I_ERI_D2x_S_G4x_S_M1_vrr+2*oned2z*I_ERI_Px_S_G4x_S_vrr-2*rhod2zsq*I_ERI_Px_S_G4x_S_M1_vrr+4*oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xy_S_G4x_S_vrr = PAY*I_ERI_D2x_S_G4x_S_vrr+WPY*I_ERI_D2x_S_G4x_S_M1_vrr;
      Double I_ERI_F2xz_S_G4x_S_vrr = PAZ*I_ERI_D2x_S_G4x_S_vrr+WPZ*I_ERI_D2x_S_G4x_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4x_S_vrr = PAX*I_ERI_D2y_S_G4x_S_vrr+WPX*I_ERI_D2y_S_G4x_S_M1_vrr+4*oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4x_S_vrr = PAX*I_ERI_D2z_S_G4x_S_vrr+WPX*I_ERI_D2z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3y_S_G4x_S_vrr = PAY*I_ERI_D2y_S_G4x_S_vrr+WPY*I_ERI_D2y_S_G4x_S_M1_vrr+2*oned2z*I_ERI_Py_S_G4x_S_vrr-2*rhod2zsq*I_ERI_Py_S_G4x_S_M1_vrr;
      Double I_ERI_F2yz_S_G4x_S_vrr = PAZ*I_ERI_D2y_S_G4x_S_vrr+WPZ*I_ERI_D2y_S_G4x_S_M1_vrr;
      Double I_ERI_F3z_S_G4x_S_vrr = PAZ*I_ERI_D2z_S_G4x_S_vrr+WPZ*I_ERI_D2z_S_G4x_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G4x_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G4x_S_M1_vrr;
      Double I_ERI_F3x_S_G3xy_S_vrr = PAX*I_ERI_D2x_S_G3xy_S_vrr+WPX*I_ERI_D2x_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_Px_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_Px_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_F2xy_S_G3xy_S_vrr = PAY*I_ERI_D2x_S_G3xy_S_vrr+WPY*I_ERI_D2x_S_G3xy_S_M1_vrr+oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xz_S_G3xy_S_vrr = PAZ*I_ERI_D2x_S_G3xy_S_vrr+WPZ*I_ERI_D2x_S_G3xy_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3xy_S_vrr = PAX*I_ERI_D2y_S_G3xy_S_vrr+WPX*I_ERI_D2y_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3xy_S_vrr = PAX*I_ERI_D2z_S_G3xy_S_vrr+WPX*I_ERI_D2z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_F3y_S_G3xy_S_vrr = PAY*I_ERI_D2y_S_G3xy_S_vrr+WPY*I_ERI_D2y_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_Py_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_Py_S_G3xy_S_M1_vrr+oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_F2yz_S_G3xy_S_vrr = PAZ*I_ERI_D2y_S_G3xy_S_vrr+WPZ*I_ERI_D2y_S_G3xy_S_M1_vrr;
      Double I_ERI_F3z_S_G3xy_S_vrr = PAZ*I_ERI_D2z_S_G3xy_S_vrr+WPZ*I_ERI_D2z_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xy_S_M1_vrr;
      Double I_ERI_F3x_S_G3xz_S_vrr = PAX*I_ERI_D2x_S_G3xz_S_vrr+WPX*I_ERI_D2x_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_Px_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_Px_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_F2xy_S_G3xz_S_vrr = PAY*I_ERI_D2x_S_G3xz_S_vrr+WPY*I_ERI_D2x_S_G3xz_S_M1_vrr;
      Double I_ERI_F2xz_S_G3xz_S_vrr = PAZ*I_ERI_D2x_S_G3xz_S_vrr+WPZ*I_ERI_D2x_S_G3xz_S_M1_vrr+oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3xz_S_vrr = PAX*I_ERI_D2y_S_G3xz_S_vrr+WPX*I_ERI_D2y_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3xz_S_vrr = PAX*I_ERI_D2z_S_G3xz_S_vrr+WPX*I_ERI_D2z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3y_S_G3xz_S_vrr = PAY*I_ERI_D2y_S_G3xz_S_vrr+WPY*I_ERI_D2y_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_Py_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_Py_S_G3xz_S_M1_vrr;
      Double I_ERI_F2yz_S_G3xz_S_vrr = PAZ*I_ERI_D2y_S_G3xz_S_vrr+WPZ*I_ERI_D2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_F3z_S_G3xz_S_vrr = PAZ*I_ERI_D2z_S_G3xz_S_vrr+WPZ*I_ERI_D2z_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xz_S_M1_vrr+oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3x_S_G2x2y_S_vrr = PAX*I_ERI_D2x_S_G2x2y_S_vrr+WPX*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2xy_S_G2x2y_S_vrr = PAY*I_ERI_D2x_S_G2x2y_S_vrr+WPY*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_F2xz_S_G2x2y_S_vrr = PAZ*I_ERI_D2x_S_G2x2y_S_vrr+WPZ*I_ERI_D2x_S_G2x2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2x2y_S_vrr = PAX*I_ERI_D2y_S_G2x2y_S_vrr+WPX*I_ERI_D2y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2x2y_S_vrr = PAX*I_ERI_D2z_S_G2x2y_S_vrr+WPX*I_ERI_D2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3y_S_G2x2y_S_vrr = PAY*I_ERI_D2y_S_G2x2y_S_vrr+WPY*I_ERI_D2y_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_F2yz_S_G2x2y_S_vrr = PAZ*I_ERI_D2y_S_G2x2y_S_vrr+WPZ*I_ERI_D2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_F3z_S_G2x2y_S_vrr = PAZ*I_ERI_D2z_S_G2x2y_S_vrr+WPZ*I_ERI_D2z_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2y_S_M1_vrr;
      Double I_ERI_F3x_S_G2xyz_S_vrr = PAX*I_ERI_D2x_S_G2xyz_S_vrr+WPX*I_ERI_D2x_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M1_vrr;
      Double I_ERI_F2xy_S_G2xyz_S_vrr = PAY*I_ERI_D2x_S_G2xyz_S_vrr+WPY*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_F2xz_S_G2xyz_S_vrr = PAZ*I_ERI_D2x_S_G2xyz_S_vrr+WPZ*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2xyz_S_vrr = PAX*I_ERI_D2y_S_G2xyz_S_vrr+WPX*I_ERI_D2y_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2xyz_S_vrr = PAX*I_ERI_D2z_S_G2xyz_S_vrr+WPX*I_ERI_D2z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_F3y_S_G2xyz_S_vrr = PAY*I_ERI_D2y_S_G2xyz_S_vrr+WPY*I_ERI_D2y_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_F2yz_S_G2xyz_S_vrr = PAZ*I_ERI_D2y_S_G2xyz_S_vrr+WPZ*I_ERI_D2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_F3z_S_G2xyz_S_vrr = PAZ*I_ERI_D2z_S_G2xyz_S_vrr+WPZ*I_ERI_D2z_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_F3x_S_G2x2z_S_vrr = PAX*I_ERI_D2x_S_G2x2z_S_vrr+WPX*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xy_S_G2x2z_S_vrr = PAY*I_ERI_D2x_S_G2x2z_S_vrr+WPY*I_ERI_D2x_S_G2x2z_S_M1_vrr;
      Double I_ERI_F2xz_S_G2x2z_S_vrr = PAZ*I_ERI_D2x_S_G2x2z_S_vrr+WPZ*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2x2z_S_vrr = PAX*I_ERI_D2y_S_G2x2z_S_vrr+WPX*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2x2z_S_vrr = PAX*I_ERI_D2z_S_G2x2z_S_vrr+WPX*I_ERI_D2z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3y_S_G2x2z_S_vrr = PAY*I_ERI_D2y_S_G2x2z_S_vrr+WPY*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2z_S_M1_vrr;
      Double I_ERI_F2yz_S_G2x2z_S_vrr = PAZ*I_ERI_D2y_S_G2x2z_S_vrr+WPZ*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_F3z_S_G2x2z_S_vrr = PAZ*I_ERI_D2z_S_G2x2z_S_vrr+WPZ*I_ERI_D2z_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3x_S_Gx3y_S_vrr = PAX*I_ERI_D2x_S_Gx3y_S_vrr+WPX*I_ERI_D2x_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3y_S_M1_vrr+oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx3y_S_vrr = PAY*I_ERI_D2x_S_Gx3y_S_vrr+WPY*I_ERI_D2x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx3y_S_vrr = PAZ*I_ERI_D2x_S_Gx3y_S_vrr+WPZ*I_ERI_D2x_S_Gx3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx3y_S_vrr = PAX*I_ERI_D2y_S_Gx3y_S_vrr+WPX*I_ERI_D2y_S_Gx3y_S_M1_vrr+oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx3y_S_vrr = PAX*I_ERI_D2z_S_Gx3y_S_vrr+WPX*I_ERI_D2z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3y_S_Gx3y_S_vrr = PAY*I_ERI_D2y_S_Gx3y_S_vrr+WPY*I_ERI_D2y_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx3y_S_vrr = PAZ*I_ERI_D2y_S_Gx3y_S_vrr+WPZ*I_ERI_D2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_F3z_S_Gx3y_S_vrr = PAZ*I_ERI_D2z_S_Gx3y_S_vrr+WPZ*I_ERI_D2z_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3y_S_M1_vrr;
      Double I_ERI_F3x_S_Gx2yz_S_vrr = PAX*I_ERI_D2x_S_Gx2yz_S_vrr+WPX*I_ERI_D2x_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx2yz_S_vrr = PAY*I_ERI_D2x_S_Gx2yz_S_vrr+WPY*I_ERI_D2x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx2yz_S_vrr = PAZ*I_ERI_D2x_S_Gx2yz_S_vrr+WPZ*I_ERI_D2x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx2yz_S_vrr = PAX*I_ERI_D2y_S_Gx2yz_S_vrr+WPX*I_ERI_D2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx2yz_S_vrr = PAX*I_ERI_D2z_S_Gx2yz_S_vrr+WPX*I_ERI_D2z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3y_S_Gx2yz_S_vrr = PAY*I_ERI_D2y_S_Gx2yz_S_vrr+WPY*I_ERI_D2y_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx2yz_S_vrr = PAZ*I_ERI_D2y_S_Gx2yz_S_vrr+WPZ*I_ERI_D2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3z_S_Gx2yz_S_vrr = PAZ*I_ERI_D2z_S_Gx2yz_S_vrr+WPZ*I_ERI_D2z_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3x_S_Gxy2z_S_vrr = PAX*I_ERI_D2x_S_Gxy2z_S_vrr+WPX*I_ERI_D2x_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gxy2z_S_vrr = PAY*I_ERI_D2x_S_Gxy2z_S_vrr+WPY*I_ERI_D2x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gxy2z_S_vrr = PAZ*I_ERI_D2x_S_Gxy2z_S_vrr+WPZ*I_ERI_D2x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gxy2z_S_vrr = PAX*I_ERI_D2y_S_Gxy2z_S_vrr+WPX*I_ERI_D2y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gxy2z_S_vrr = PAX*I_ERI_D2z_S_Gxy2z_S_vrr+WPX*I_ERI_D2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3y_S_Gxy2z_S_vrr = PAY*I_ERI_D2y_S_Gxy2z_S_vrr+WPY*I_ERI_D2y_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gxy2z_S_vrr = PAZ*I_ERI_D2y_S_Gxy2z_S_vrr+WPZ*I_ERI_D2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_F3z_S_Gxy2z_S_vrr = PAZ*I_ERI_D2z_S_Gxy2z_S_vrr+WPZ*I_ERI_D2z_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_F3x_S_Gx3z_S_vrr = PAX*I_ERI_D2x_S_Gx3z_S_vrr+WPX*I_ERI_D2x_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3z_S_M1_vrr+oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx3z_S_vrr = PAY*I_ERI_D2x_S_Gx3z_S_vrr+WPY*I_ERI_D2x_S_Gx3z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx3z_S_vrr = PAZ*I_ERI_D2x_S_Gx3z_S_vrr+WPZ*I_ERI_D2x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx3z_S_vrr = PAX*I_ERI_D2y_S_Gx3z_S_vrr+WPX*I_ERI_D2y_S_Gx3z_S_M1_vrr+oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx3z_S_vrr = PAX*I_ERI_D2z_S_Gx3z_S_vrr+WPX*I_ERI_D2z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3y_S_Gx3z_S_vrr = PAY*I_ERI_D2y_S_Gx3z_S_vrr+WPY*I_ERI_D2y_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx3z_S_vrr = PAZ*I_ERI_D2y_S_Gx3z_S_vrr+WPZ*I_ERI_D2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3z_S_Gx3z_S_vrr = PAZ*I_ERI_D2z_S_Gx3z_S_vrr+WPZ*I_ERI_D2z_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3x_S_G4y_S_vrr = PAX*I_ERI_D2x_S_G4y_S_vrr+WPX*I_ERI_D2x_S_G4y_S_M1_vrr+2*oned2z*I_ERI_Px_S_G4y_S_vrr-2*rhod2zsq*I_ERI_Px_S_G4y_S_M1_vrr;
      Double I_ERI_F2xy_S_G4y_S_vrr = PAY*I_ERI_D2x_S_G4y_S_vrr+WPY*I_ERI_D2x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_F2xz_S_G4y_S_vrr = PAZ*I_ERI_D2x_S_G4y_S_vrr+WPZ*I_ERI_D2x_S_G4y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4y_S_vrr = PAX*I_ERI_D2y_S_G4y_S_vrr+WPX*I_ERI_D2y_S_G4y_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4y_S_vrr = PAX*I_ERI_D2z_S_G4y_S_vrr+WPX*I_ERI_D2z_S_G4y_S_M1_vrr;
      Double I_ERI_F3y_S_G4y_S_vrr = PAY*I_ERI_D2y_S_G4y_S_vrr+WPY*I_ERI_D2y_S_G4y_S_M1_vrr+2*oned2z*I_ERI_Py_S_G4y_S_vrr-2*rhod2zsq*I_ERI_Py_S_G4y_S_M1_vrr+4*oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_F2yz_S_G4y_S_vrr = PAZ*I_ERI_D2y_S_G4y_S_vrr+WPZ*I_ERI_D2y_S_G4y_S_M1_vrr;
      Double I_ERI_F3z_S_G4y_S_vrr = PAZ*I_ERI_D2z_S_G4y_S_vrr+WPZ*I_ERI_D2z_S_G4y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G4y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G4y_S_M1_vrr;
      Double I_ERI_F3x_S_G3yz_S_vrr = PAX*I_ERI_D2x_S_G3yz_S_vrr+WPX*I_ERI_D2x_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_Px_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_Px_S_G3yz_S_M1_vrr;
      Double I_ERI_F2xy_S_G3yz_S_vrr = PAY*I_ERI_D2x_S_G3yz_S_vrr+WPY*I_ERI_D2x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_F2xz_S_G3yz_S_vrr = PAZ*I_ERI_D2x_S_G3yz_S_vrr+WPZ*I_ERI_D2x_S_G3yz_S_M1_vrr+oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3yz_S_vrr = PAX*I_ERI_D2y_S_G3yz_S_vrr+WPX*I_ERI_D2y_S_G3yz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3yz_S_vrr = PAX*I_ERI_D2z_S_G3yz_S_vrr+WPX*I_ERI_D2z_S_G3yz_S_M1_vrr;
      Double I_ERI_F3y_S_G3yz_S_vrr = PAY*I_ERI_D2y_S_G3yz_S_vrr+WPY*I_ERI_D2y_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_Py_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_Py_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_F2yz_S_G3yz_S_vrr = PAZ*I_ERI_D2y_S_G3yz_S_vrr+WPZ*I_ERI_D2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_F3z_S_G3yz_S_vrr = PAZ*I_ERI_D2z_S_G3yz_S_vrr+WPZ*I_ERI_D2z_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G3yz_S_M1_vrr+oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3x_S_G2y2z_S_vrr = PAX*I_ERI_D2x_S_G2y2z_S_vrr+WPX*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2y2z_S_M1_vrr;
      Double I_ERI_F2xy_S_G2y2z_S_vrr = PAY*I_ERI_D2x_S_G2y2z_S_vrr+WPY*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2xz_S_G2y2z_S_vrr = PAZ*I_ERI_D2x_S_G2y2z_S_vrr+WPZ*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2y2z_S_vrr = PAX*I_ERI_D2y_S_G2y2z_S_vrr+WPX*I_ERI_D2y_S_G2y2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2y2z_S_vrr = PAX*I_ERI_D2z_S_G2y2z_S_vrr+WPX*I_ERI_D2z_S_G2y2z_S_M1_vrr;
      Double I_ERI_F3y_S_G2y2z_S_vrr = PAY*I_ERI_D2y_S_G2y2z_S_vrr+WPY*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2yz_S_G2y2z_S_vrr = PAZ*I_ERI_D2y_S_G2y2z_S_vrr+WPZ*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_F3z_S_G2y2z_S_vrr = PAZ*I_ERI_D2z_S_G2y2z_S_vrr+WPZ*I_ERI_D2z_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3x_S_Gy3z_S_vrr = PAX*I_ERI_D2x_S_Gy3z_S_vrr+WPX*I_ERI_D2x_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gy3z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gy3z_S_vrr = PAY*I_ERI_D2x_S_Gy3z_S_vrr+WPY*I_ERI_D2x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gy3z_S_vrr = PAZ*I_ERI_D2x_S_Gy3z_S_vrr+WPZ*I_ERI_D2x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gy3z_S_vrr = PAX*I_ERI_D2y_S_Gy3z_S_vrr+WPX*I_ERI_D2y_S_Gy3z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gy3z_S_vrr = PAX*I_ERI_D2z_S_Gy3z_S_vrr+WPX*I_ERI_D2z_S_Gy3z_S_M1_vrr;
      Double I_ERI_F3y_S_Gy3z_S_vrr = PAY*I_ERI_D2y_S_Gy3z_S_vrr+WPY*I_ERI_D2y_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gy3z_S_M1_vrr+oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gy3z_S_vrr = PAZ*I_ERI_D2y_S_Gy3z_S_vrr+WPZ*I_ERI_D2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3z_S_Gy3z_S_vrr = PAZ*I_ERI_D2z_S_Gy3z_S_vrr+WPZ*I_ERI_D2z_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3x_S_G4z_S_vrr = PAX*I_ERI_D2x_S_G4z_S_vrr+WPX*I_ERI_D2x_S_G4z_S_M1_vrr+2*oned2z*I_ERI_Px_S_G4z_S_vrr-2*rhod2zsq*I_ERI_Px_S_G4z_S_M1_vrr;
      Double I_ERI_F2xy_S_G4z_S_vrr = PAY*I_ERI_D2x_S_G4z_S_vrr+WPY*I_ERI_D2x_S_G4z_S_M1_vrr;
      Double I_ERI_F2xz_S_G4z_S_vrr = PAZ*I_ERI_D2x_S_G4z_S_vrr+WPZ*I_ERI_D2x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4z_S_vrr = PAX*I_ERI_D2y_S_G4z_S_vrr+WPX*I_ERI_D2y_S_G4z_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4z_S_vrr = PAX*I_ERI_D2z_S_G4z_S_vrr+WPX*I_ERI_D2z_S_G4z_S_M1_vrr;
      Double I_ERI_F3y_S_G4z_S_vrr = PAY*I_ERI_D2y_S_G4z_S_vrr+WPY*I_ERI_D2y_S_G4z_S_M1_vrr+2*oned2z*I_ERI_Py_S_G4z_S_vrr-2*rhod2zsq*I_ERI_Py_S_G4z_S_M1_vrr;
      Double I_ERI_F2yz_S_G4z_S_vrr = PAZ*I_ERI_D2y_S_G4z_S_vrr+WPZ*I_ERI_D2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_F3z_S_G4z_S_vrr = PAZ*I_ERI_D2z_S_G4z_S_vrr+WPZ*I_ERI_D2z_S_G4z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G4z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G4z_S_M1_vrr+4*oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_G_S
       * RHS shell quartet name: SQ_ERI_F_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_G_S
       * RHS shell quartet name: SQ_ERI_D_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_G4x_S_vrr = PAX*I_ERI_F3x_S_G4x_S_vrr+WPX*I_ERI_F3x_S_G4x_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G4x_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G4x_S_M1_vrr+4*oned2k*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G3xy_S_G4x_S_vrr = PAY*I_ERI_F3x_S_G4x_S_vrr+WPY*I_ERI_F3x_S_G4x_S_M1_vrr;
      Double I_ERI_G3xz_S_G4x_S_vrr = PAZ*I_ERI_F3x_S_G4x_S_vrr+WPZ*I_ERI_F3x_S_G4x_S_M1_vrr;
      Double I_ERI_G2x2y_S_G4x_S_vrr = PAY*I_ERI_F2xy_S_G4x_S_vrr+WPY*I_ERI_F2xy_S_G4x_S_M1_vrr+oned2z*I_ERI_D2x_S_G4x_S_vrr-rhod2zsq*I_ERI_D2x_S_G4x_S_M1_vrr;
      Double I_ERI_G2xyz_S_G4x_S_vrr = PAZ*I_ERI_F2xy_S_G4x_S_vrr+WPZ*I_ERI_F2xy_S_G4x_S_M1_vrr;
      Double I_ERI_G2x2z_S_G4x_S_vrr = PAZ*I_ERI_F2xz_S_G4x_S_vrr+WPZ*I_ERI_F2xz_S_G4x_S_M1_vrr+oned2z*I_ERI_D2x_S_G4x_S_vrr-rhod2zsq*I_ERI_D2x_S_G4x_S_M1_vrr;
      Double I_ERI_Gx3y_S_G4x_S_vrr = PAX*I_ERI_F3y_S_G4x_S_vrr+WPX*I_ERI_F3y_S_G4x_S_M1_vrr+4*oned2k*I_ERI_F3y_S_F3x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G4x_S_vrr = PAZ*I_ERI_Fx2y_S_G4x_S_vrr+WPZ*I_ERI_Fx2y_S_G4x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G4x_S_vrr = PAY*I_ERI_Fx2z_S_G4x_S_vrr+WPY*I_ERI_Fx2z_S_G4x_S_M1_vrr;
      Double I_ERI_Gx3z_S_G4x_S_vrr = PAX*I_ERI_F3z_S_G4x_S_vrr+WPX*I_ERI_F3z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_F3z_S_F3x_S_M1_vrr;
      Double I_ERI_G4y_S_G4x_S_vrr = PAY*I_ERI_F3y_S_G4x_S_vrr+WPY*I_ERI_F3y_S_G4x_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G4x_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G4x_S_M1_vrr;
      Double I_ERI_G3yz_S_G4x_S_vrr = PAZ*I_ERI_F3y_S_G4x_S_vrr+WPZ*I_ERI_F3y_S_G4x_S_M1_vrr;
      Double I_ERI_G2y2z_S_G4x_S_vrr = PAZ*I_ERI_F2yz_S_G4x_S_vrr+WPZ*I_ERI_F2yz_S_G4x_S_M1_vrr+oned2z*I_ERI_D2y_S_G4x_S_vrr-rhod2zsq*I_ERI_D2y_S_G4x_S_M1_vrr;
      Double I_ERI_Gy3z_S_G4x_S_vrr = PAY*I_ERI_F3z_S_G4x_S_vrr+WPY*I_ERI_F3z_S_G4x_S_M1_vrr;
      Double I_ERI_G4z_S_G4x_S_vrr = PAZ*I_ERI_F3z_S_G4x_S_vrr+WPZ*I_ERI_F3z_S_G4x_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G4x_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G4x_S_M1_vrr;
      Double I_ERI_G4x_S_G3xy_S_vrr = PAX*I_ERI_F3x_S_G3xy_S_vrr+WPX*I_ERI_F3x_S_G3xy_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G3xy_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_F3x_S_F2xy_S_M1_vrr;
      Double I_ERI_G3xy_S_G3xy_S_vrr = PAY*I_ERI_F3x_S_G3xy_S_vrr+WPY*I_ERI_F3x_S_G3xy_S_M1_vrr+oned2k*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G3xz_S_G3xy_S_vrr = PAZ*I_ERI_F3x_S_G3xy_S_vrr+WPZ*I_ERI_F3x_S_G3xy_S_M1_vrr;
      Double I_ERI_G2x2y_S_G3xy_S_vrr = PAY*I_ERI_F2xy_S_G3xy_S_vrr+WPY*I_ERI_F2xy_S_G3xy_S_M1_vrr+oned2z*I_ERI_D2x_S_G3xy_S_vrr-rhod2zsq*I_ERI_D2x_S_G3xy_S_M1_vrr+oned2k*I_ERI_F2xy_S_F3x_S_M1_vrr;
      Double I_ERI_G2xyz_S_G3xy_S_vrr = PAZ*I_ERI_F2xy_S_G3xy_S_vrr+WPZ*I_ERI_F2xy_S_G3xy_S_M1_vrr;
      Double I_ERI_G2x2z_S_G3xy_S_vrr = PAZ*I_ERI_F2xz_S_G3xy_S_vrr+WPZ*I_ERI_F2xz_S_G3xy_S_M1_vrr+oned2z*I_ERI_D2x_S_G3xy_S_vrr-rhod2zsq*I_ERI_D2x_S_G3xy_S_M1_vrr;
      Double I_ERI_Gx3y_S_G3xy_S_vrr = PAX*I_ERI_F3y_S_G3xy_S_vrr+WPX*I_ERI_F3y_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_F3y_S_F2xy_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G3xy_S_vrr = PAZ*I_ERI_Fx2y_S_G3xy_S_vrr+WPZ*I_ERI_Fx2y_S_G3xy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G3xy_S_vrr = PAY*I_ERI_Fx2z_S_G3xy_S_vrr+WPY*I_ERI_Fx2z_S_G3xy_S_M1_vrr+oned2k*I_ERI_Fx2z_S_F3x_S_M1_vrr;
      Double I_ERI_Gx3z_S_G3xy_S_vrr = PAX*I_ERI_F3z_S_G3xy_S_vrr+WPX*I_ERI_F3z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_F3z_S_F2xy_S_M1_vrr;
      Double I_ERI_G4y_S_G3xy_S_vrr = PAY*I_ERI_F3y_S_G3xy_S_vrr+WPY*I_ERI_F3y_S_G3xy_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G3xy_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G3xy_S_M1_vrr+oned2k*I_ERI_F3y_S_F3x_S_M1_vrr;
      Double I_ERI_G3yz_S_G3xy_S_vrr = PAZ*I_ERI_F3y_S_G3xy_S_vrr+WPZ*I_ERI_F3y_S_G3xy_S_M1_vrr;
      Double I_ERI_G2y2z_S_G3xy_S_vrr = PAZ*I_ERI_F2yz_S_G3xy_S_vrr+WPZ*I_ERI_F2yz_S_G3xy_S_M1_vrr+oned2z*I_ERI_D2y_S_G3xy_S_vrr-rhod2zsq*I_ERI_D2y_S_G3xy_S_M1_vrr;
      Double I_ERI_Gy3z_S_G3xy_S_vrr = PAY*I_ERI_F3z_S_G3xy_S_vrr+WPY*I_ERI_F3z_S_G3xy_S_M1_vrr+oned2k*I_ERI_F3z_S_F3x_S_M1_vrr;
      Double I_ERI_G4z_S_G3xy_S_vrr = PAZ*I_ERI_F3z_S_G3xy_S_vrr+WPZ*I_ERI_F3z_S_G3xy_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G3xy_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G3xy_S_M1_vrr;
      Double I_ERI_G4x_S_G3xz_S_vrr = PAX*I_ERI_F3x_S_G3xz_S_vrr+WPX*I_ERI_F3x_S_G3xz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G3xz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_F3x_S_F2xz_S_M1_vrr;
      Double I_ERI_G3xy_S_G3xz_S_vrr = PAY*I_ERI_F3x_S_G3xz_S_vrr+WPY*I_ERI_F3x_S_G3xz_S_M1_vrr;
      Double I_ERI_G3xz_S_G3xz_S_vrr = PAZ*I_ERI_F3x_S_G3xz_S_vrr+WPZ*I_ERI_F3x_S_G3xz_S_M1_vrr+oned2k*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G2x2y_S_G3xz_S_vrr = PAY*I_ERI_F2xy_S_G3xz_S_vrr+WPY*I_ERI_F2xy_S_G3xz_S_M1_vrr+oned2z*I_ERI_D2x_S_G3xz_S_vrr-rhod2zsq*I_ERI_D2x_S_G3xz_S_M1_vrr;
      Double I_ERI_G2xyz_S_G3xz_S_vrr = PAZ*I_ERI_F2xy_S_G3xz_S_vrr+WPZ*I_ERI_F2xy_S_G3xz_S_M1_vrr+oned2k*I_ERI_F2xy_S_F3x_S_M1_vrr;
      Double I_ERI_G2x2z_S_G3xz_S_vrr = PAZ*I_ERI_F2xz_S_G3xz_S_vrr+WPZ*I_ERI_F2xz_S_G3xz_S_M1_vrr+oned2z*I_ERI_D2x_S_G3xz_S_vrr-rhod2zsq*I_ERI_D2x_S_G3xz_S_M1_vrr+oned2k*I_ERI_F2xz_S_F3x_S_M1_vrr;
      Double I_ERI_Gx3y_S_G3xz_S_vrr = PAX*I_ERI_F3y_S_G3xz_S_vrr+WPX*I_ERI_F3y_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_F3y_S_F2xz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G3xz_S_vrr = PAZ*I_ERI_Fx2y_S_G3xz_S_vrr+WPZ*I_ERI_Fx2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_F3x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G3xz_S_vrr = PAY*I_ERI_Fx2z_S_G3xz_S_vrr+WPY*I_ERI_Fx2z_S_G3xz_S_M1_vrr;
      Double I_ERI_Gx3z_S_G3xz_S_vrr = PAX*I_ERI_F3z_S_G3xz_S_vrr+WPX*I_ERI_F3z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_F3z_S_F2xz_S_M1_vrr;
      Double I_ERI_G4y_S_G3xz_S_vrr = PAY*I_ERI_F3y_S_G3xz_S_vrr+WPY*I_ERI_F3y_S_G3xz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G3xz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G3xz_S_M1_vrr;
      Double I_ERI_G3yz_S_G3xz_S_vrr = PAZ*I_ERI_F3y_S_G3xz_S_vrr+WPZ*I_ERI_F3y_S_G3xz_S_M1_vrr+oned2k*I_ERI_F3y_S_F3x_S_M1_vrr;
      Double I_ERI_G2y2z_S_G3xz_S_vrr = PAZ*I_ERI_F2yz_S_G3xz_S_vrr+WPZ*I_ERI_F2yz_S_G3xz_S_M1_vrr+oned2z*I_ERI_D2y_S_G3xz_S_vrr-rhod2zsq*I_ERI_D2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_F2yz_S_F3x_S_M1_vrr;
      Double I_ERI_Gy3z_S_G3xz_S_vrr = PAY*I_ERI_F3z_S_G3xz_S_vrr+WPY*I_ERI_F3z_S_G3xz_S_M1_vrr;
      Double I_ERI_G4z_S_G3xz_S_vrr = PAZ*I_ERI_F3z_S_G3xz_S_vrr+WPZ*I_ERI_F3z_S_G3xz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G3xz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G3xz_S_M1_vrr+oned2k*I_ERI_F3z_S_F3x_S_M1_vrr;
      Double I_ERI_G4x_S_G2x2y_S_vrr = PAX*I_ERI_F3x_S_G2x2y_S_vrr+WPX*I_ERI_F3x_S_G2x2y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G2x2y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fx2y_S_M1_vrr;
      Double I_ERI_G3xy_S_G2x2y_S_vrr = PAY*I_ERI_F3x_S_G2x2y_S_vrr+WPY*I_ERI_F3x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3x_S_F2xy_S_M1_vrr;
      Double I_ERI_G3xz_S_G2x2y_S_vrr = PAZ*I_ERI_F3x_S_G2x2y_S_vrr+WPZ*I_ERI_F3x_S_G2x2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_G2x2y_S_vrr = PAY*I_ERI_F2xy_S_G2x2y_S_vrr+WPY*I_ERI_F2xy_S_G2x2y_S_M1_vrr+oned2z*I_ERI_D2x_S_G2x2y_S_vrr-rhod2zsq*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_F2xy_S_M1_vrr;
      Double I_ERI_G2xyz_S_G2x2y_S_vrr = PAZ*I_ERI_F2xy_S_G2x2y_S_vrr+WPZ*I_ERI_F2xy_S_G2x2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_G2x2y_S_vrr = PAZ*I_ERI_F2xz_S_G2x2y_S_vrr+WPZ*I_ERI_F2xz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_D2x_S_G2x2y_S_vrr-rhod2zsq*I_ERI_D2x_S_G2x2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_G2x2y_S_vrr = PAX*I_ERI_F3y_S_G2x2y_S_vrr+WPX*I_ERI_F3y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G2x2y_S_vrr = PAZ*I_ERI_Fx2y_S_G2x2y_S_vrr+WPZ*I_ERI_Fx2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G2x2y_S_vrr = PAY*I_ERI_Fx2z_S_G2x2y_S_vrr+WPY*I_ERI_Fx2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_F2xy_S_M1_vrr;
      Double I_ERI_Gx3z_S_G2x2y_S_vrr = PAX*I_ERI_F3z_S_G2x2y_S_vrr+WPX*I_ERI_F3z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_G4y_S_G2x2y_S_vrr = PAY*I_ERI_F3y_S_G2x2y_S_vrr+WPY*I_ERI_F3y_S_G2x2y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G2x2y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3y_S_F2xy_S_M1_vrr;
      Double I_ERI_G3yz_S_G2x2y_S_vrr = PAZ*I_ERI_F3y_S_G2x2y_S_vrr+WPZ*I_ERI_F3y_S_G2x2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_G2x2y_S_vrr = PAZ*I_ERI_F2yz_S_G2x2y_S_vrr+WPZ*I_ERI_F2yz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_D2y_S_G2x2y_S_vrr-rhod2zsq*I_ERI_D2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_G2x2y_S_vrr = PAY*I_ERI_F3z_S_G2x2y_S_vrr+WPY*I_ERI_F3z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_F3z_S_F2xy_S_M1_vrr;
      Double I_ERI_G4z_S_G2x2y_S_vrr = PAZ*I_ERI_F3z_S_G2x2y_S_vrr+WPZ*I_ERI_F3z_S_G2x2y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G2x2y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G2x2y_S_M1_vrr;
      Double I_ERI_G4x_S_G2xyz_S_vrr = PAX*I_ERI_F3x_S_G2xyz_S_vrr+WPX*I_ERI_F3x_S_G2xyz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G2xyz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M1_vrr;
      Double I_ERI_G3xy_S_G2xyz_S_vrr = PAY*I_ERI_F3x_S_G2xyz_S_vrr+WPY*I_ERI_F3x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3x_S_F2xz_S_M1_vrr;
      Double I_ERI_G3xz_S_G2xyz_S_vrr = PAZ*I_ERI_F3x_S_G2xyz_S_vrr+WPZ*I_ERI_F3x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3x_S_F2xy_S_M1_vrr;
      Double I_ERI_G2x2y_S_G2xyz_S_vrr = PAY*I_ERI_F2xy_S_G2xyz_S_vrr+WPY*I_ERI_F2xy_S_G2xyz_S_M1_vrr+oned2z*I_ERI_D2x_S_G2xyz_S_vrr-rhod2zsq*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_F2xz_S_M1_vrr;
      Double I_ERI_G2xyz_S_G2xyz_S_vrr = PAZ*I_ERI_F2xy_S_G2xyz_S_vrr+WPZ*I_ERI_F2xy_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_F2xy_S_M1_vrr;
      Double I_ERI_G2x2z_S_G2xyz_S_vrr = PAZ*I_ERI_F2xz_S_G2xyz_S_vrr+WPZ*I_ERI_F2xz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_D2x_S_G2xyz_S_vrr-rhod2zsq*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F2xz_S_F2xy_S_M1_vrr;
      Double I_ERI_Gx3y_S_G2xyz_S_vrr = PAX*I_ERI_F3y_S_G2xyz_S_vrr+WPX*I_ERI_F3y_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G2xyz_S_vrr = PAZ*I_ERI_Fx2y_S_G2xyz_S_vrr+WPZ*I_ERI_Fx2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G2xyz_S_vrr = PAY*I_ERI_Fx2z_S_G2xyz_S_vrr+WPY*I_ERI_Fx2z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Fx2z_S_F2xz_S_M1_vrr;
      Double I_ERI_Gx3z_S_G2xyz_S_vrr = PAX*I_ERI_F3z_S_G2xyz_S_vrr+WPX*I_ERI_F3z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_G4y_S_G2xyz_S_vrr = PAY*I_ERI_F3y_S_G2xyz_S_vrr+WPY*I_ERI_F3y_S_G2xyz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G2xyz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3y_S_F2xz_S_M1_vrr;
      Double I_ERI_G3yz_S_G2xyz_S_vrr = PAZ*I_ERI_F3y_S_G2xyz_S_vrr+WPZ*I_ERI_F3y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3y_S_F2xy_S_M1_vrr;
      Double I_ERI_G2y2z_S_G2xyz_S_vrr = PAZ*I_ERI_F2yz_S_G2xyz_S_vrr+WPZ*I_ERI_F2yz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_D2y_S_G2xyz_S_vrr-rhod2zsq*I_ERI_D2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F2yz_S_F2xy_S_M1_vrr;
      Double I_ERI_Gy3z_S_G2xyz_S_vrr = PAY*I_ERI_F3z_S_G2xyz_S_vrr+WPY*I_ERI_F3z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3z_S_F2xz_S_M1_vrr;
      Double I_ERI_G4z_S_G2xyz_S_vrr = PAZ*I_ERI_F3z_S_G2xyz_S_vrr+WPZ*I_ERI_F3z_S_G2xyz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G2xyz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_F3z_S_F2xy_S_M1_vrr;
      Double I_ERI_G4x_S_G2x2z_S_vrr = PAX*I_ERI_F3x_S_G2x2z_S_vrr+WPX*I_ERI_F3x_S_G2x2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G2x2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3xy_S_G2x2z_S_vrr = PAY*I_ERI_F3x_S_G2x2z_S_vrr+WPY*I_ERI_F3x_S_G2x2z_S_M1_vrr;
      Double I_ERI_G3xz_S_G2x2z_S_vrr = PAZ*I_ERI_F3x_S_G2x2z_S_vrr+WPZ*I_ERI_F3x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_F2xz_S_M1_vrr;
      Double I_ERI_G2x2y_S_G2x2z_S_vrr = PAY*I_ERI_F2xy_S_G2x2z_S_vrr+WPY*I_ERI_F2xy_S_G2x2z_S_M1_vrr+oned2z*I_ERI_D2x_S_G2x2z_S_vrr-rhod2zsq*I_ERI_D2x_S_G2x2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_G2x2z_S_vrr = PAZ*I_ERI_F2xy_S_G2x2z_S_vrr+WPZ*I_ERI_F2xy_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_F2xz_S_M1_vrr;
      Double I_ERI_G2x2z_S_G2x2z_S_vrr = PAZ*I_ERI_F2xz_S_G2x2z_S_vrr+WPZ*I_ERI_F2xz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_D2x_S_G2x2z_S_vrr-rhod2zsq*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_F2xz_S_M1_vrr;
      Double I_ERI_Gx3y_S_G2x2z_S_vrr = PAX*I_ERI_F3y_S_G2x2z_S_vrr+WPX*I_ERI_F3y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G2x2z_S_vrr = PAZ*I_ERI_Fx2y_S_G2x2z_S_vrr+WPZ*I_ERI_Fx2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_F2xz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G2x2z_S_vrr = PAY*I_ERI_Fx2z_S_G2x2z_S_vrr+WPY*I_ERI_Fx2z_S_G2x2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_G2x2z_S_vrr = PAX*I_ERI_F3z_S_G2x2z_S_vrr+WPX*I_ERI_F3z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_G4y_S_G2x2z_S_vrr = PAY*I_ERI_F3y_S_G2x2z_S_vrr+WPY*I_ERI_F3y_S_G2x2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G2x2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G2x2z_S_M1_vrr;
      Double I_ERI_G3yz_S_G2x2z_S_vrr = PAZ*I_ERI_F3y_S_G2x2z_S_vrr+WPZ*I_ERI_F3y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_F2xz_S_M1_vrr;
      Double I_ERI_G2y2z_S_G2x2z_S_vrr = PAZ*I_ERI_F2yz_S_G2x2z_S_vrr+WPZ*I_ERI_F2yz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_D2y_S_G2x2z_S_vrr-rhod2zsq*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_F2xz_S_M1_vrr;
      Double I_ERI_Gy3z_S_G2x2z_S_vrr = PAY*I_ERI_F3z_S_G2x2z_S_vrr+WPY*I_ERI_F3z_S_G2x2z_S_M1_vrr;
      Double I_ERI_G4z_S_G2x2z_S_vrr = PAZ*I_ERI_F3z_S_G2x2z_S_vrr+WPZ*I_ERI_F3z_S_G2x2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G2x2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_F2xz_S_M1_vrr;
      Double I_ERI_G4x_S_Gx3y_S_vrr = PAX*I_ERI_F3x_S_Gx3y_S_vrr+WPX*I_ERI_F3x_S_Gx3y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gx3y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx3y_S_M1_vrr+oned2k*I_ERI_F3x_S_F3y_S_M1_vrr;
      Double I_ERI_G3xy_S_Gx3y_S_vrr = PAY*I_ERI_F3x_S_Gx3y_S_vrr+WPY*I_ERI_F3x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_F3x_S_Fx2y_S_M1_vrr;
      Double I_ERI_G3xz_S_Gx3y_S_vrr = PAZ*I_ERI_F3x_S_Gx3y_S_vrr+WPZ*I_ERI_F3x_S_Gx3y_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gx3y_S_vrr = PAY*I_ERI_F2xy_S_Gx3y_S_vrr+WPY*I_ERI_F2xy_S_Gx3y_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx3y_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gx3y_S_vrr = PAZ*I_ERI_F2xy_S_Gx3y_S_vrr+WPZ*I_ERI_F2xy_S_Gx3y_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gx3y_S_vrr = PAZ*I_ERI_F2xz_S_Gx3y_S_vrr+WPZ*I_ERI_F2xz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx3y_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx3y_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gx3y_S_vrr = PAX*I_ERI_F3y_S_Gx3y_S_vrr+WPX*I_ERI_F3y_S_Gx3y_S_M1_vrr+oned2k*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gx3y_S_vrr = PAZ*I_ERI_Fx2y_S_Gx3y_S_vrr+WPZ*I_ERI_Fx2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gx3y_S_vrr = PAY*I_ERI_Fx2z_S_Gx3y_S_vrr+WPY*I_ERI_Fx2z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gx3y_S_vrr = PAX*I_ERI_F3z_S_Gx3y_S_vrr+WPX*I_ERI_F3z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_F3z_S_F3y_S_M1_vrr;
      Double I_ERI_G4y_S_Gx3y_S_vrr = PAY*I_ERI_F3y_S_Gx3y_S_vrr+WPY*I_ERI_F3y_S_Gx3y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gx3y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_F3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_G3yz_S_Gx3y_S_vrr = PAZ*I_ERI_F3y_S_Gx3y_S_vrr+WPZ*I_ERI_F3y_S_Gx3y_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gx3y_S_vrr = PAZ*I_ERI_F2yz_S_Gx3y_S_vrr+WPZ*I_ERI_F2yz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_D2y_S_Gx3y_S_vrr-rhod2zsq*I_ERI_D2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gx3y_S_vrr = PAY*I_ERI_F3z_S_Gx3y_S_vrr+WPY*I_ERI_F3z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_F3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_G4z_S_Gx3y_S_vrr = PAZ*I_ERI_F3z_S_Gx3y_S_vrr+WPZ*I_ERI_F3z_S_Gx3y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gx3y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx3y_S_M1_vrr;
      Double I_ERI_G4x_S_Gx2yz_S_vrr = PAX*I_ERI_F3x_S_Gx2yz_S_vrr+WPX*I_ERI_F3x_S_Gx2yz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gx2yz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3x_S_F2yz_S_M1_vrr;
      Double I_ERI_G3xy_S_Gx2yz_S_vrr = PAY*I_ERI_F3x_S_Gx2yz_S_vrr+WPY*I_ERI_F3x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M1_vrr;
      Double I_ERI_G3xz_S_Gx2yz_S_vrr = PAZ*I_ERI_F3x_S_Gx2yz_S_vrr+WPZ*I_ERI_F3x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3x_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gx2yz_S_vrr = PAY*I_ERI_F2xy_S_Gx2yz_S_vrr+WPY*I_ERI_F2xy_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Fxyz_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gx2yz_S_vrr = PAZ*I_ERI_F2xy_S_Gx2yz_S_vrr+WPZ*I_ERI_F2xy_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F2xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gx2yz_S_vrr = PAZ*I_ERI_F2xz_S_Gx2yz_S_vrr+WPZ*I_ERI_F2xz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F2xz_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gx2yz_S_vrr = PAX*I_ERI_F3y_S_Gx2yz_S_vrr+WPX*I_ERI_F3y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3y_S_F2yz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gx2yz_S_vrr = PAZ*I_ERI_Fx2y_S_Gx2yz_S_vrr+WPZ*I_ERI_Fx2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gx2yz_S_vrr = PAY*I_ERI_Fx2z_S_Gx2yz_S_vrr+WPY*I_ERI_Fx2z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gx2yz_S_vrr = PAX*I_ERI_F3z_S_Gx2yz_S_vrr+WPX*I_ERI_F3z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3z_S_F2yz_S_M1_vrr;
      Double I_ERI_G4y_S_Gx2yz_S_vrr = PAY*I_ERI_F3y_S_Gx2yz_S_vrr+WPY*I_ERI_F3y_S_Gx2yz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gx2yz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M1_vrr;
      Double I_ERI_G3yz_S_Gx2yz_S_vrr = PAZ*I_ERI_F3y_S_Gx2yz_S_vrr+WPZ*I_ERI_F3y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gx2yz_S_vrr = PAZ*I_ERI_F2yz_S_Gx2yz_S_vrr+WPZ*I_ERI_F2yz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_D2y_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_D2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F2yz_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gx2yz_S_vrr = PAY*I_ERI_F3z_S_Gx2yz_S_vrr+WPY*I_ERI_F3z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_G4z_S_Gx2yz_S_vrr = PAZ*I_ERI_F3z_S_Gx2yz_S_vrr+WPZ*I_ERI_F3z_S_Gx2yz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gx2yz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_F3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_G4x_S_Gxy2z_S_vrr = PAX*I_ERI_F3x_S_Gxy2z_S_vrr+WPX*I_ERI_F3x_S_Gxy2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gxy2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3x_S_Fy2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Gxy2z_S_vrr = PAY*I_ERI_F3x_S_Gxy2z_S_vrr+WPY*I_ERI_F3x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Gxy2z_S_vrr = PAZ*I_ERI_F3x_S_Gxy2z_S_vrr+WPZ*I_ERI_F3x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fxyz_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gxy2z_S_vrr = PAY*I_ERI_F2xy_S_Gxy2z_S_vrr+WPY*I_ERI_F2xy_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F2xy_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gxy2z_S_vrr = PAZ*I_ERI_F2xy_S_Gxy2z_S_vrr+WPZ*I_ERI_F2xy_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Fxyz_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gxy2z_S_vrr = PAZ*I_ERI_F2xz_S_Gxy2z_S_vrr+WPZ*I_ERI_F2xz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gxy2z_S_vrr = PAX*I_ERI_F3y_S_Gxy2z_S_vrr+WPX*I_ERI_F3y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gxy2z_S_vrr = PAZ*I_ERI_Fx2y_S_Gxy2z_S_vrr+WPZ*I_ERI_Fx2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gxy2z_S_vrr = PAY*I_ERI_Fx2z_S_Gxy2z_S_vrr+WPY*I_ERI_Fx2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Fx2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gxy2z_S_vrr = PAX*I_ERI_F3z_S_Gxy2z_S_vrr+WPX*I_ERI_F3z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_G4y_S_Gxy2z_S_vrr = PAY*I_ERI_F3y_S_Gxy2z_S_vrr+WPY*I_ERI_F3y_S_Gxy2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gxy2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3y_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Gxy2z_S_vrr = PAZ*I_ERI_F3y_S_Gxy2z_S_vrr+WPZ*I_ERI_F3y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fxyz_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gxy2z_S_vrr = PAZ*I_ERI_F2yz_S_Gxy2z_S_vrr+WPZ*I_ERI_F2yz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_D2y_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_D2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Fxyz_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gxy2z_S_vrr = PAY*I_ERI_F3z_S_Gxy2z_S_vrr+WPY*I_ERI_F3z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_F3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_G4z_S_Gxy2z_S_vrr = PAZ*I_ERI_F3z_S_Gxy2z_S_vrr+WPZ*I_ERI_F3z_S_Gxy2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gxy2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_G4x_S_Gx3z_S_vrr = PAX*I_ERI_F3x_S_Gx3z_S_vrr+WPX*I_ERI_F3x_S_Gx3z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gx3z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gx3z_S_M1_vrr+oned2k*I_ERI_F3x_S_F3z_S_M1_vrr;
      Double I_ERI_G3xy_S_Gx3z_S_vrr = PAY*I_ERI_F3x_S_Gx3z_S_vrr+WPY*I_ERI_F3x_S_Gx3z_S_M1_vrr;
      Double I_ERI_G3xz_S_Gx3z_S_vrr = PAZ*I_ERI_F3x_S_Gx3z_S_vrr+WPZ*I_ERI_F3x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F3x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gx3z_S_vrr = PAY*I_ERI_F2xy_S_Gx3z_S_vrr+WPY*I_ERI_F2xy_S_Gx3z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx3z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx3z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gx3z_S_vrr = PAZ*I_ERI_F2xy_S_Gx3z_S_vrr+WPZ*I_ERI_F2xy_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gx3z_S_vrr = PAZ*I_ERI_F2xz_S_Gx3z_S_vrr+WPZ*I_ERI_F2xz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gx3z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gx3z_S_vrr = PAX*I_ERI_F3y_S_Gx3z_S_vrr+WPX*I_ERI_F3y_S_Gx3z_S_M1_vrr+oned2k*I_ERI_F3y_S_F3z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gx3z_S_vrr = PAZ*I_ERI_Fx2y_S_Gx3z_S_vrr+WPZ*I_ERI_Fx2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gx3z_S_vrr = PAY*I_ERI_Fx2z_S_Gx3z_S_vrr+WPY*I_ERI_Fx2z_S_Gx3z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gx3z_S_vrr = PAX*I_ERI_F3z_S_Gx3z_S_vrr+WPX*I_ERI_F3z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_G4y_S_Gx3z_S_vrr = PAY*I_ERI_F3y_S_Gx3z_S_vrr+WPY*I_ERI_F3y_S_Gx3z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gx3z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gx3z_S_M1_vrr;
      Double I_ERI_G3yz_S_Gx3z_S_vrr = PAZ*I_ERI_F3y_S_Gx3z_S_vrr+WPZ*I_ERI_F3y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F3y_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gx3z_S_vrr = PAZ*I_ERI_F2yz_S_Gx3z_S_vrr+WPZ*I_ERI_F2yz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_D2y_S_Gx3z_S_vrr-rhod2zsq*I_ERI_D2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gx3z_S_vrr = PAY*I_ERI_F3z_S_Gx3z_S_vrr+WPY*I_ERI_F3z_S_Gx3z_S_M1_vrr;
      Double I_ERI_G4z_S_Gx3z_S_vrr = PAZ*I_ERI_F3z_S_Gx3z_S_vrr+WPZ*I_ERI_F3z_S_Gx3z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gx3z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_F3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_G4x_S_G4y_S_vrr = PAX*I_ERI_F3x_S_G4y_S_vrr+WPX*I_ERI_F3x_S_G4y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G4y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G4y_S_M1_vrr;
      Double I_ERI_G3xy_S_G4y_S_vrr = PAY*I_ERI_F3x_S_G4y_S_vrr+WPY*I_ERI_F3x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_F3x_S_F3y_S_M1_vrr;
      Double I_ERI_G3xz_S_G4y_S_vrr = PAZ*I_ERI_F3x_S_G4y_S_vrr+WPZ*I_ERI_F3x_S_G4y_S_M1_vrr;
      Double I_ERI_G2x2y_S_G4y_S_vrr = PAY*I_ERI_F2xy_S_G4y_S_vrr+WPY*I_ERI_F2xy_S_G4y_S_M1_vrr+oned2z*I_ERI_D2x_S_G4y_S_vrr-rhod2zsq*I_ERI_D2x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_F2xy_S_F3y_S_M1_vrr;
      Double I_ERI_G2xyz_S_G4y_S_vrr = PAZ*I_ERI_F2xy_S_G4y_S_vrr+WPZ*I_ERI_F2xy_S_G4y_S_M1_vrr;
      Double I_ERI_G2x2z_S_G4y_S_vrr = PAZ*I_ERI_F2xz_S_G4y_S_vrr+WPZ*I_ERI_F2xz_S_G4y_S_M1_vrr+oned2z*I_ERI_D2x_S_G4y_S_vrr-rhod2zsq*I_ERI_D2x_S_G4y_S_M1_vrr;
      Double I_ERI_Gx3y_S_G4y_S_vrr = PAX*I_ERI_F3y_S_G4y_S_vrr+WPX*I_ERI_F3y_S_G4y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G4y_S_vrr = PAZ*I_ERI_Fx2y_S_G4y_S_vrr+WPZ*I_ERI_Fx2y_S_G4y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G4y_S_vrr = PAY*I_ERI_Fx2z_S_G4y_S_vrr+WPY*I_ERI_Fx2z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Fx2z_S_F3y_S_M1_vrr;
      Double I_ERI_Gx3z_S_G4y_S_vrr = PAX*I_ERI_F3z_S_G4y_S_vrr+WPX*I_ERI_F3z_S_G4y_S_M1_vrr;
      Double I_ERI_G4y_S_G4y_S_vrr = PAY*I_ERI_F3y_S_G4y_S_vrr+WPY*I_ERI_F3y_S_G4y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G4y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G4y_S_M1_vrr+4*oned2k*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_G3yz_S_G4y_S_vrr = PAZ*I_ERI_F3y_S_G4y_S_vrr+WPZ*I_ERI_F3y_S_G4y_S_M1_vrr;
      Double I_ERI_G2y2z_S_G4y_S_vrr = PAZ*I_ERI_F2yz_S_G4y_S_vrr+WPZ*I_ERI_F2yz_S_G4y_S_M1_vrr+oned2z*I_ERI_D2y_S_G4y_S_vrr-rhod2zsq*I_ERI_D2y_S_G4y_S_M1_vrr;
      Double I_ERI_Gy3z_S_G4y_S_vrr = PAY*I_ERI_F3z_S_G4y_S_vrr+WPY*I_ERI_F3z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_F3z_S_F3y_S_M1_vrr;
      Double I_ERI_G4z_S_G4y_S_vrr = PAZ*I_ERI_F3z_S_G4y_S_vrr+WPZ*I_ERI_F3z_S_G4y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G4y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G4y_S_M1_vrr;
      Double I_ERI_G4x_S_G3yz_S_vrr = PAX*I_ERI_F3x_S_G3yz_S_vrr+WPX*I_ERI_F3x_S_G3yz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G3yz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G3yz_S_M1_vrr;
      Double I_ERI_G3xy_S_G3yz_S_vrr = PAY*I_ERI_F3x_S_G3yz_S_vrr+WPY*I_ERI_F3x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_F3x_S_F2yz_S_M1_vrr;
      Double I_ERI_G3xz_S_G3yz_S_vrr = PAZ*I_ERI_F3x_S_G3yz_S_vrr+WPZ*I_ERI_F3x_S_G3yz_S_M1_vrr+oned2k*I_ERI_F3x_S_F3y_S_M1_vrr;
      Double I_ERI_G2x2y_S_G3yz_S_vrr = PAY*I_ERI_F2xy_S_G3yz_S_vrr+WPY*I_ERI_F2xy_S_G3yz_S_M1_vrr+oned2z*I_ERI_D2x_S_G3yz_S_vrr-rhod2zsq*I_ERI_D2x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_F2yz_S_M1_vrr;
      Double I_ERI_G2xyz_S_G3yz_S_vrr = PAZ*I_ERI_F2xy_S_G3yz_S_vrr+WPZ*I_ERI_F2xy_S_G3yz_S_M1_vrr+oned2k*I_ERI_F2xy_S_F3y_S_M1_vrr;
      Double I_ERI_G2x2z_S_G3yz_S_vrr = PAZ*I_ERI_F2xz_S_G3yz_S_vrr+WPZ*I_ERI_F2xz_S_G3yz_S_M1_vrr+oned2z*I_ERI_D2x_S_G3yz_S_vrr-rhod2zsq*I_ERI_D2x_S_G3yz_S_M1_vrr+oned2k*I_ERI_F2xz_S_F3y_S_M1_vrr;
      Double I_ERI_Gx3y_S_G3yz_S_vrr = PAX*I_ERI_F3y_S_G3yz_S_vrr+WPX*I_ERI_F3y_S_G3yz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G3yz_S_vrr = PAZ*I_ERI_Fx2y_S_G3yz_S_vrr+WPZ*I_ERI_Fx2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_F3y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G3yz_S_vrr = PAY*I_ERI_Fx2z_S_G3yz_S_vrr+WPY*I_ERI_Fx2z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_F2yz_S_M1_vrr;
      Double I_ERI_Gx3z_S_G3yz_S_vrr = PAX*I_ERI_F3z_S_G3yz_S_vrr+WPX*I_ERI_F3z_S_G3yz_S_M1_vrr;
      Double I_ERI_G4y_S_G3yz_S_vrr = PAY*I_ERI_F3y_S_G3yz_S_vrr+WPY*I_ERI_F3y_S_G3yz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G3yz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_F3y_S_F2yz_S_M1_vrr;
      Double I_ERI_G3yz_S_G3yz_S_vrr = PAZ*I_ERI_F3y_S_G3yz_S_vrr+WPZ*I_ERI_F3y_S_G3yz_S_M1_vrr+oned2k*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_G2y2z_S_G3yz_S_vrr = PAZ*I_ERI_F2yz_S_G3yz_S_vrr+WPZ*I_ERI_F2yz_S_G3yz_S_M1_vrr+oned2z*I_ERI_D2y_S_G3yz_S_vrr-rhod2zsq*I_ERI_D2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_F2yz_S_F3y_S_M1_vrr;
      Double I_ERI_Gy3z_S_G3yz_S_vrr = PAY*I_ERI_F3z_S_G3yz_S_vrr+WPY*I_ERI_F3z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_F3z_S_F2yz_S_M1_vrr;
      Double I_ERI_G4z_S_G3yz_S_vrr = PAZ*I_ERI_F3z_S_G3yz_S_vrr+WPZ*I_ERI_F3z_S_G3yz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G3yz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G3yz_S_M1_vrr+oned2k*I_ERI_F3z_S_F3y_S_M1_vrr;
      Double I_ERI_G4x_S_G2y2z_S_vrr = PAX*I_ERI_F3x_S_G2y2z_S_vrr+WPX*I_ERI_F3x_S_G2y2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G2y2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G2y2z_S_M1_vrr;
      Double I_ERI_G3xy_S_G2y2z_S_vrr = PAY*I_ERI_F3x_S_G2y2z_S_vrr+WPY*I_ERI_F3x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Fy2z_S_M1_vrr;
      Double I_ERI_G3xz_S_G2y2z_S_vrr = PAZ*I_ERI_F3x_S_G2y2z_S_vrr+WPZ*I_ERI_F3x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_F2yz_S_M1_vrr;
      Double I_ERI_G2x2y_S_G2y2z_S_vrr = PAY*I_ERI_F2xy_S_G2y2z_S_vrr+WPY*I_ERI_F2xy_S_G2y2z_S_M1_vrr+oned2z*I_ERI_D2x_S_G2y2z_S_vrr-rhod2zsq*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Fy2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_G2y2z_S_vrr = PAZ*I_ERI_F2xy_S_G2y2z_S_vrr+WPZ*I_ERI_F2xy_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_F2yz_S_M1_vrr;
      Double I_ERI_G2x2z_S_G2y2z_S_vrr = PAZ*I_ERI_F2xz_S_G2y2z_S_vrr+WPZ*I_ERI_F2xz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_D2x_S_G2y2z_S_vrr-rhod2zsq*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_F2yz_S_M1_vrr;
      Double I_ERI_Gx3y_S_G2y2z_S_vrr = PAX*I_ERI_F3y_S_G2y2z_S_vrr+WPX*I_ERI_F3y_S_G2y2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G2y2z_S_vrr = PAZ*I_ERI_Fx2y_S_G2y2z_S_vrr+WPZ*I_ERI_Fx2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_F2yz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G2y2z_S_vrr = PAY*I_ERI_Fx2z_S_G2y2z_S_vrr+WPY*I_ERI_Fx2z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_G2y2z_S_vrr = PAX*I_ERI_F3z_S_G2y2z_S_vrr+WPX*I_ERI_F3z_S_G2y2z_S_M1_vrr;
      Double I_ERI_G4y_S_G2y2z_S_vrr = PAY*I_ERI_F3y_S_G2y2z_S_vrr+WPY*I_ERI_F3y_S_G2y2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G2y2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_G3yz_S_G2y2z_S_vrr = PAZ*I_ERI_F3y_S_G2y2z_S_vrr+WPZ*I_ERI_F3y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_F2yz_S_M1_vrr;
      Double I_ERI_G2y2z_S_G2y2z_S_vrr = PAZ*I_ERI_F2yz_S_G2y2z_S_vrr+WPZ*I_ERI_F2yz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_D2y_S_G2y2z_S_vrr-rhod2zsq*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_F2yz_S_M1_vrr;
      Double I_ERI_Gy3z_S_G2y2z_S_vrr = PAY*I_ERI_F3z_S_G2y2z_S_vrr+WPY*I_ERI_F3z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_G4z_S_G2y2z_S_vrr = PAZ*I_ERI_F3z_S_G2y2z_S_vrr+WPZ*I_ERI_F3z_S_G2y2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G2y2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_F2yz_S_M1_vrr;
      Double I_ERI_G4x_S_Gy3z_S_vrr = PAX*I_ERI_F3x_S_Gy3z_S_vrr+WPX*I_ERI_F3x_S_Gy3z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Gy3z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Gy3z_S_M1_vrr;
      Double I_ERI_G3xy_S_Gy3z_S_vrr = PAY*I_ERI_F3x_S_Gy3z_S_vrr+WPY*I_ERI_F3x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_F3x_S_F3z_S_M1_vrr;
      Double I_ERI_G3xz_S_Gy3z_S_vrr = PAZ*I_ERI_F3x_S_Gy3z_S_vrr+WPZ*I_ERI_F3x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F3x_S_Fy2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_Gy3z_S_vrr = PAY*I_ERI_F2xy_S_Gy3z_S_vrr+WPY*I_ERI_F2xy_S_Gy3z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gy3z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_F2xy_S_F3z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Gy3z_S_vrr = PAZ*I_ERI_F2xy_S_Gy3z_S_vrr+WPZ*I_ERI_F2xy_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_Fy2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_Gy3z_S_vrr = PAZ*I_ERI_F2xz_S_Gy3z_S_vrr+WPZ*I_ERI_F2xz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_D2x_S_Gy3z_S_vrr-rhod2zsq*I_ERI_D2x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_Gy3z_S_vrr = PAX*I_ERI_F3y_S_Gy3z_S_vrr+WPX*I_ERI_F3y_S_Gy3z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Gy3z_S_vrr = PAZ*I_ERI_Fx2y_S_Gy3z_S_vrr+WPZ*I_ERI_Fx2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Gy3z_S_vrr = PAY*I_ERI_Fx2z_S_Gy3z_S_vrr+WPY*I_ERI_Fx2z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Fx2z_S_F3z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Gy3z_S_vrr = PAX*I_ERI_F3z_S_Gy3z_S_vrr+WPX*I_ERI_F3z_S_Gy3z_S_M1_vrr;
      Double I_ERI_G4y_S_Gy3z_S_vrr = PAY*I_ERI_F3y_S_Gy3z_S_vrr+WPY*I_ERI_F3y_S_Gy3z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Gy3z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Gy3z_S_M1_vrr+oned2k*I_ERI_F3y_S_F3z_S_M1_vrr;
      Double I_ERI_G3yz_S_Gy3z_S_vrr = PAZ*I_ERI_F3y_S_Gy3z_S_vrr+WPZ*I_ERI_F3y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_Gy3z_S_vrr = PAZ*I_ERI_F2yz_S_Gy3z_S_vrr+WPZ*I_ERI_F2yz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_D2y_S_Gy3z_S_vrr-rhod2zsq*I_ERI_D2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_Gy3z_S_vrr = PAY*I_ERI_F3z_S_Gy3z_S_vrr+WPY*I_ERI_F3z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_G4z_S_Gy3z_S_vrr = PAZ*I_ERI_F3z_S_Gy3z_S_vrr+WPZ*I_ERI_F3z_S_Gy3z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Gy3z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_F3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_G4x_S_G4z_S_vrr = PAX*I_ERI_F3x_S_G4z_S_vrr+WPX*I_ERI_F3x_S_G4z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_G4z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_G4z_S_M1_vrr;
      Double I_ERI_G3xy_S_G4z_S_vrr = PAY*I_ERI_F3x_S_G4z_S_vrr+WPY*I_ERI_F3x_S_G4z_S_M1_vrr;
      Double I_ERI_G3xz_S_G4z_S_vrr = PAZ*I_ERI_F3x_S_G4z_S_vrr+WPZ*I_ERI_F3x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F3x_S_F3z_S_M1_vrr;
      Double I_ERI_G2x2y_S_G4z_S_vrr = PAY*I_ERI_F2xy_S_G4z_S_vrr+WPY*I_ERI_F2xy_S_G4z_S_M1_vrr+oned2z*I_ERI_D2x_S_G4z_S_vrr-rhod2zsq*I_ERI_D2x_S_G4z_S_M1_vrr;
      Double I_ERI_G2xyz_S_G4z_S_vrr = PAZ*I_ERI_F2xy_S_G4z_S_vrr+WPZ*I_ERI_F2xy_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F2xy_S_F3z_S_M1_vrr;
      Double I_ERI_G2x2z_S_G4z_S_vrr = PAZ*I_ERI_F2xz_S_G4z_S_vrr+WPZ*I_ERI_F2xz_S_G4z_S_M1_vrr+oned2z*I_ERI_D2x_S_G4z_S_vrr-rhod2zsq*I_ERI_D2x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F2xz_S_F3z_S_M1_vrr;
      Double I_ERI_Gx3y_S_G4z_S_vrr = PAX*I_ERI_F3y_S_G4z_S_vrr+WPX*I_ERI_F3y_S_G4z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_G4z_S_vrr = PAZ*I_ERI_Fx2y_S_G4z_S_vrr+WPZ*I_ERI_Fx2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Fx2y_S_F3z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_G4z_S_vrr = PAY*I_ERI_Fx2z_S_G4z_S_vrr+WPY*I_ERI_Fx2z_S_G4z_S_M1_vrr;
      Double I_ERI_Gx3z_S_G4z_S_vrr = PAX*I_ERI_F3z_S_G4z_S_vrr+WPX*I_ERI_F3z_S_G4z_S_M1_vrr;
      Double I_ERI_G4y_S_G4z_S_vrr = PAY*I_ERI_F3y_S_G4z_S_vrr+WPY*I_ERI_F3y_S_G4z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_G4z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_G4z_S_M1_vrr;
      Double I_ERI_G3yz_S_G4z_S_vrr = PAZ*I_ERI_F3y_S_G4z_S_vrr+WPZ*I_ERI_F3y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F3y_S_F3z_S_M1_vrr;
      Double I_ERI_G2y2z_S_G4z_S_vrr = PAZ*I_ERI_F2yz_S_G4z_S_vrr+WPZ*I_ERI_F2yz_S_G4z_S_M1_vrr+oned2z*I_ERI_D2y_S_G4z_S_vrr-rhod2zsq*I_ERI_D2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F2yz_S_F3z_S_M1_vrr;
      Double I_ERI_Gy3z_S_G4z_S_vrr = PAY*I_ERI_F3z_S_G4z_S_vrr+WPY*I_ERI_F3z_S_G4z_S_M1_vrr;
      Double I_ERI_G4z_S_G4z_S_vrr = PAZ*I_ERI_F3z_S_G4z_S_vrr+WPZ*I_ERI_F3z_S_G4z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_G4z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_G4z_S_M1_vrr+4*oned2k*I_ERI_F3z_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_G_S
       * RHS shell quartet name: SQ_ERI_G_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_G_S
       * RHS shell quartet name: SQ_ERI_F_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_F_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_G4x_S_vrr = PAX*I_ERI_G4x_S_G4x_S_vrr+WPX*I_ERI_G4x_S_G4x_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G4x_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G4x_S_M1_vrr+4*oned2k*I_ERI_G4x_S_F3x_S_M1_vrr;
      Double I_ERI_H4xy_S_G4x_S_vrr = PAY*I_ERI_G4x_S_G4x_S_vrr+WPY*I_ERI_G4x_S_G4x_S_M1_vrr;
      Double I_ERI_H4xz_S_G4x_S_vrr = PAZ*I_ERI_G4x_S_G4x_S_vrr+WPZ*I_ERI_G4x_S_G4x_S_M1_vrr;
      Double I_ERI_H3x2y_S_G4x_S_vrr = PAY*I_ERI_G3xy_S_G4x_S_vrr+WPY*I_ERI_G3xy_S_G4x_S_M1_vrr+oned2z*I_ERI_F3x_S_G4x_S_vrr-rhod2zsq*I_ERI_F3x_S_G4x_S_M1_vrr;
      Double I_ERI_H3xyz_S_G4x_S_vrr = PAZ*I_ERI_G3xy_S_G4x_S_vrr+WPZ*I_ERI_G3xy_S_G4x_S_M1_vrr;
      Double I_ERI_H3x2z_S_G4x_S_vrr = PAZ*I_ERI_G3xz_S_G4x_S_vrr+WPZ*I_ERI_G3xz_S_G4x_S_M1_vrr+oned2z*I_ERI_F3x_S_G4x_S_vrr-rhod2zsq*I_ERI_F3x_S_G4x_S_M1_vrr;
      Double I_ERI_H2x3y_S_G4x_S_vrr = PAX*I_ERI_Gx3y_S_G4x_S_vrr+WPX*I_ERI_Gx3y_S_G4x_S_M1_vrr+oned2z*I_ERI_F3y_S_G4x_S_vrr-rhod2zsq*I_ERI_F3y_S_G4x_S_M1_vrr+4*oned2k*I_ERI_Gx3y_S_F3x_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G4x_S_vrr = PAZ*I_ERI_G2x2y_S_G4x_S_vrr+WPZ*I_ERI_G2x2y_S_G4x_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G4x_S_vrr = PAY*I_ERI_G2x2z_S_G4x_S_vrr+WPY*I_ERI_G2x2z_S_G4x_S_M1_vrr;
      Double I_ERI_H2x3z_S_G4x_S_vrr = PAX*I_ERI_Gx3z_S_G4x_S_vrr+WPX*I_ERI_Gx3z_S_G4x_S_M1_vrr+oned2z*I_ERI_F3z_S_G4x_S_vrr-rhod2zsq*I_ERI_F3z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_Gx3z_S_F3x_S_M1_vrr;
      Double I_ERI_Hx4y_S_G4x_S_vrr = PAX*I_ERI_G4y_S_G4x_S_vrr+WPX*I_ERI_G4y_S_G4x_S_M1_vrr+4*oned2k*I_ERI_G4y_S_F3x_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G4x_S_vrr = PAZ*I_ERI_Gx3y_S_G4x_S_vrr+WPZ*I_ERI_Gx3y_S_G4x_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G4x_S_vrr = PAX*I_ERI_G2y2z_S_G4x_S_vrr+WPX*I_ERI_G2y2z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_G2y2z_S_F3x_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G4x_S_vrr = PAY*I_ERI_Gx3z_S_G4x_S_vrr+WPY*I_ERI_Gx3z_S_G4x_S_M1_vrr;
      Double I_ERI_Hx4z_S_G4x_S_vrr = PAX*I_ERI_G4z_S_G4x_S_vrr+WPX*I_ERI_G4z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_G4z_S_F3x_S_M1_vrr;
      Double I_ERI_H5y_S_G4x_S_vrr = PAY*I_ERI_G4y_S_G4x_S_vrr+WPY*I_ERI_G4y_S_G4x_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G4x_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G4x_S_M1_vrr;
      Double I_ERI_H4yz_S_G4x_S_vrr = PAZ*I_ERI_G4y_S_G4x_S_vrr+WPZ*I_ERI_G4y_S_G4x_S_M1_vrr;
      Double I_ERI_H3y2z_S_G4x_S_vrr = PAZ*I_ERI_G3yz_S_G4x_S_vrr+WPZ*I_ERI_G3yz_S_G4x_S_M1_vrr+oned2z*I_ERI_F3y_S_G4x_S_vrr-rhod2zsq*I_ERI_F3y_S_G4x_S_M1_vrr;
      Double I_ERI_H2y3z_S_G4x_S_vrr = PAY*I_ERI_Gy3z_S_G4x_S_vrr+WPY*I_ERI_Gy3z_S_G4x_S_M1_vrr+oned2z*I_ERI_F3z_S_G4x_S_vrr-rhod2zsq*I_ERI_F3z_S_G4x_S_M1_vrr;
      Double I_ERI_Hy4z_S_G4x_S_vrr = PAY*I_ERI_G4z_S_G4x_S_vrr+WPY*I_ERI_G4z_S_G4x_S_M1_vrr;
      Double I_ERI_H5z_S_G4x_S_vrr = PAZ*I_ERI_G4z_S_G4x_S_vrr+WPZ*I_ERI_G4z_S_G4x_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G4x_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G4x_S_M1_vrr;
      Double I_ERI_H5x_S_G3xy_S_vrr = PAX*I_ERI_G4x_S_G3xy_S_vrr+WPX*I_ERI_G4x_S_G3xy_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G3xy_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_G4x_S_F2xy_S_M1_vrr;
      Double I_ERI_H4xy_S_G3xy_S_vrr = PAY*I_ERI_G4x_S_G3xy_S_vrr+WPY*I_ERI_G4x_S_G3xy_S_M1_vrr+oned2k*I_ERI_G4x_S_F3x_S_M1_vrr;
      Double I_ERI_H4xz_S_G3xy_S_vrr = PAZ*I_ERI_G4x_S_G3xy_S_vrr+WPZ*I_ERI_G4x_S_G3xy_S_M1_vrr;
      Double I_ERI_H3x2y_S_G3xy_S_vrr = PAY*I_ERI_G3xy_S_G3xy_S_vrr+WPY*I_ERI_G3xy_S_G3xy_S_M1_vrr+oned2z*I_ERI_F3x_S_G3xy_S_vrr-rhod2zsq*I_ERI_F3x_S_G3xy_S_M1_vrr+oned2k*I_ERI_G3xy_S_F3x_S_M1_vrr;
      Double I_ERI_H3xyz_S_G3xy_S_vrr = PAZ*I_ERI_G3xy_S_G3xy_S_vrr+WPZ*I_ERI_G3xy_S_G3xy_S_M1_vrr;
      Double I_ERI_H3x2z_S_G3xy_S_vrr = PAZ*I_ERI_G3xz_S_G3xy_S_vrr+WPZ*I_ERI_G3xz_S_G3xy_S_M1_vrr+oned2z*I_ERI_F3x_S_G3xy_S_vrr-rhod2zsq*I_ERI_F3x_S_G3xy_S_M1_vrr;
      Double I_ERI_H2x3y_S_G3xy_S_vrr = PAX*I_ERI_Gx3y_S_G3xy_S_vrr+WPX*I_ERI_Gx3y_S_G3xy_S_M1_vrr+oned2z*I_ERI_F3y_S_G3xy_S_vrr-rhod2zsq*I_ERI_F3y_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_Gx3y_S_F2xy_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G3xy_S_vrr = PAZ*I_ERI_G2x2y_S_G3xy_S_vrr+WPZ*I_ERI_G2x2y_S_G3xy_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G3xy_S_vrr = PAY*I_ERI_G2x2z_S_G3xy_S_vrr+WPY*I_ERI_G2x2z_S_G3xy_S_M1_vrr+oned2k*I_ERI_G2x2z_S_F3x_S_M1_vrr;
      Double I_ERI_H2x3z_S_G3xy_S_vrr = PAX*I_ERI_Gx3z_S_G3xy_S_vrr+WPX*I_ERI_Gx3z_S_G3xy_S_M1_vrr+oned2z*I_ERI_F3z_S_G3xy_S_vrr-rhod2zsq*I_ERI_F3z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_Gx3z_S_F2xy_S_M1_vrr;
      Double I_ERI_Hx4y_S_G3xy_S_vrr = PAX*I_ERI_G4y_S_G3xy_S_vrr+WPX*I_ERI_G4y_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_G4y_S_F2xy_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G3xy_S_vrr = PAZ*I_ERI_Gx3y_S_G3xy_S_vrr+WPZ*I_ERI_Gx3y_S_G3xy_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G3xy_S_vrr = PAX*I_ERI_G2y2z_S_G3xy_S_vrr+WPX*I_ERI_G2y2z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_F2xy_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G3xy_S_vrr = PAY*I_ERI_Gx3z_S_G3xy_S_vrr+WPY*I_ERI_Gx3z_S_G3xy_S_M1_vrr+oned2k*I_ERI_Gx3z_S_F3x_S_M1_vrr;
      Double I_ERI_Hx4z_S_G3xy_S_vrr = PAX*I_ERI_G4z_S_G3xy_S_vrr+WPX*I_ERI_G4z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_G4z_S_F2xy_S_M1_vrr;
      Double I_ERI_H5y_S_G3xy_S_vrr = PAY*I_ERI_G4y_S_G3xy_S_vrr+WPY*I_ERI_G4y_S_G3xy_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G3xy_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G3xy_S_M1_vrr+oned2k*I_ERI_G4y_S_F3x_S_M1_vrr;
      Double I_ERI_H4yz_S_G3xy_S_vrr = PAZ*I_ERI_G4y_S_G3xy_S_vrr+WPZ*I_ERI_G4y_S_G3xy_S_M1_vrr;
      Double I_ERI_H3y2z_S_G3xy_S_vrr = PAZ*I_ERI_G3yz_S_G3xy_S_vrr+WPZ*I_ERI_G3yz_S_G3xy_S_M1_vrr+oned2z*I_ERI_F3y_S_G3xy_S_vrr-rhod2zsq*I_ERI_F3y_S_G3xy_S_M1_vrr;
      Double I_ERI_H2y3z_S_G3xy_S_vrr = PAY*I_ERI_Gy3z_S_G3xy_S_vrr+WPY*I_ERI_Gy3z_S_G3xy_S_M1_vrr+oned2z*I_ERI_F3z_S_G3xy_S_vrr-rhod2zsq*I_ERI_F3z_S_G3xy_S_M1_vrr+oned2k*I_ERI_Gy3z_S_F3x_S_M1_vrr;
      Double I_ERI_Hy4z_S_G3xy_S_vrr = PAY*I_ERI_G4z_S_G3xy_S_vrr+WPY*I_ERI_G4z_S_G3xy_S_M1_vrr+oned2k*I_ERI_G4z_S_F3x_S_M1_vrr;
      Double I_ERI_H5z_S_G3xy_S_vrr = PAZ*I_ERI_G4z_S_G3xy_S_vrr+WPZ*I_ERI_G4z_S_G3xy_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G3xy_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G3xy_S_M1_vrr;
      Double I_ERI_H5x_S_G3xz_S_vrr = PAX*I_ERI_G4x_S_G3xz_S_vrr+WPX*I_ERI_G4x_S_G3xz_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G3xz_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_G4x_S_F2xz_S_M1_vrr;
      Double I_ERI_H4xy_S_G3xz_S_vrr = PAY*I_ERI_G4x_S_G3xz_S_vrr+WPY*I_ERI_G4x_S_G3xz_S_M1_vrr;
      Double I_ERI_H4xz_S_G3xz_S_vrr = PAZ*I_ERI_G4x_S_G3xz_S_vrr+WPZ*I_ERI_G4x_S_G3xz_S_M1_vrr+oned2k*I_ERI_G4x_S_F3x_S_M1_vrr;
      Double I_ERI_H3x2y_S_G3xz_S_vrr = PAY*I_ERI_G3xy_S_G3xz_S_vrr+WPY*I_ERI_G3xy_S_G3xz_S_M1_vrr+oned2z*I_ERI_F3x_S_G3xz_S_vrr-rhod2zsq*I_ERI_F3x_S_G3xz_S_M1_vrr;
      Double I_ERI_H3xyz_S_G3xz_S_vrr = PAZ*I_ERI_G3xy_S_G3xz_S_vrr+WPZ*I_ERI_G3xy_S_G3xz_S_M1_vrr+oned2k*I_ERI_G3xy_S_F3x_S_M1_vrr;
      Double I_ERI_H3x2z_S_G3xz_S_vrr = PAZ*I_ERI_G3xz_S_G3xz_S_vrr+WPZ*I_ERI_G3xz_S_G3xz_S_M1_vrr+oned2z*I_ERI_F3x_S_G3xz_S_vrr-rhod2zsq*I_ERI_F3x_S_G3xz_S_M1_vrr+oned2k*I_ERI_G3xz_S_F3x_S_M1_vrr;
      Double I_ERI_H2x3y_S_G3xz_S_vrr = PAX*I_ERI_Gx3y_S_G3xz_S_vrr+WPX*I_ERI_Gx3y_S_G3xz_S_M1_vrr+oned2z*I_ERI_F3y_S_G3xz_S_vrr-rhod2zsq*I_ERI_F3y_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_Gx3y_S_F2xz_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G3xz_S_vrr = PAZ*I_ERI_G2x2y_S_G3xz_S_vrr+WPZ*I_ERI_G2x2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_G2x2y_S_F3x_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G3xz_S_vrr = PAY*I_ERI_G2x2z_S_G3xz_S_vrr+WPY*I_ERI_G2x2z_S_G3xz_S_M1_vrr;
      Double I_ERI_H2x3z_S_G3xz_S_vrr = PAX*I_ERI_Gx3z_S_G3xz_S_vrr+WPX*I_ERI_Gx3z_S_G3xz_S_M1_vrr+oned2z*I_ERI_F3z_S_G3xz_S_vrr-rhod2zsq*I_ERI_F3z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_Gx3z_S_F2xz_S_M1_vrr;
      Double I_ERI_Hx4y_S_G3xz_S_vrr = PAX*I_ERI_G4y_S_G3xz_S_vrr+WPX*I_ERI_G4y_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_G4y_S_F2xz_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G3xz_S_vrr = PAZ*I_ERI_Gx3y_S_G3xz_S_vrr+WPZ*I_ERI_Gx3y_S_G3xz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_F3x_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G3xz_S_vrr = PAX*I_ERI_G2y2z_S_G3xz_S_vrr+WPX*I_ERI_G2y2z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_F2xz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G3xz_S_vrr = PAY*I_ERI_Gx3z_S_G3xz_S_vrr+WPY*I_ERI_Gx3z_S_G3xz_S_M1_vrr;
      Double I_ERI_Hx4z_S_G3xz_S_vrr = PAX*I_ERI_G4z_S_G3xz_S_vrr+WPX*I_ERI_G4z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_G4z_S_F2xz_S_M1_vrr;
      Double I_ERI_H5y_S_G3xz_S_vrr = PAY*I_ERI_G4y_S_G3xz_S_vrr+WPY*I_ERI_G4y_S_G3xz_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G3xz_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G3xz_S_M1_vrr;
      Double I_ERI_H4yz_S_G3xz_S_vrr = PAZ*I_ERI_G4y_S_G3xz_S_vrr+WPZ*I_ERI_G4y_S_G3xz_S_M1_vrr+oned2k*I_ERI_G4y_S_F3x_S_M1_vrr;
      Double I_ERI_H3y2z_S_G3xz_S_vrr = PAZ*I_ERI_G3yz_S_G3xz_S_vrr+WPZ*I_ERI_G3yz_S_G3xz_S_M1_vrr+oned2z*I_ERI_F3y_S_G3xz_S_vrr-rhod2zsq*I_ERI_F3y_S_G3xz_S_M1_vrr+oned2k*I_ERI_G3yz_S_F3x_S_M1_vrr;
      Double I_ERI_H2y3z_S_G3xz_S_vrr = PAY*I_ERI_Gy3z_S_G3xz_S_vrr+WPY*I_ERI_Gy3z_S_G3xz_S_M1_vrr+oned2z*I_ERI_F3z_S_G3xz_S_vrr-rhod2zsq*I_ERI_F3z_S_G3xz_S_M1_vrr;
      Double I_ERI_Hy4z_S_G3xz_S_vrr = PAY*I_ERI_G4z_S_G3xz_S_vrr+WPY*I_ERI_G4z_S_G3xz_S_M1_vrr;
      Double I_ERI_H5z_S_G3xz_S_vrr = PAZ*I_ERI_G4z_S_G3xz_S_vrr+WPZ*I_ERI_G4z_S_G3xz_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G3xz_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G3xz_S_M1_vrr+oned2k*I_ERI_G4z_S_F3x_S_M1_vrr;
      Double I_ERI_H5x_S_G2x2y_S_vrr = PAX*I_ERI_G4x_S_G2x2y_S_vrr+WPX*I_ERI_G4x_S_G2x2y_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G2x2y_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_G4x_S_Fx2y_S_M1_vrr;
      Double I_ERI_H4xy_S_G2x2y_S_vrr = PAY*I_ERI_G4x_S_G2x2y_S_vrr+WPY*I_ERI_G4x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_G4x_S_F2xy_S_M1_vrr;
      Double I_ERI_H4xz_S_G2x2y_S_vrr = PAZ*I_ERI_G4x_S_G2x2y_S_vrr+WPZ*I_ERI_G4x_S_G2x2y_S_M1_vrr;
      Double I_ERI_H3x2y_S_G2x2y_S_vrr = PAY*I_ERI_G3xy_S_G2x2y_S_vrr+WPY*I_ERI_G3xy_S_G2x2y_S_M1_vrr+oned2z*I_ERI_F3x_S_G2x2y_S_vrr-rhod2zsq*I_ERI_F3x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_F2xy_S_M1_vrr;
      Double I_ERI_H3xyz_S_G2x2y_S_vrr = PAZ*I_ERI_G3xy_S_G2x2y_S_vrr+WPZ*I_ERI_G3xy_S_G2x2y_S_M1_vrr;
      Double I_ERI_H3x2z_S_G2x2y_S_vrr = PAZ*I_ERI_G3xz_S_G2x2y_S_vrr+WPZ*I_ERI_G3xz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_F3x_S_G2x2y_S_vrr-rhod2zsq*I_ERI_F3x_S_G2x2y_S_M1_vrr;
      Double I_ERI_H2x3y_S_G2x2y_S_vrr = PAX*I_ERI_Gx3y_S_G2x2y_S_vrr+WPX*I_ERI_Gx3y_S_G2x2y_S_M1_vrr+oned2z*I_ERI_F3y_S_G2x2y_S_vrr-rhod2zsq*I_ERI_F3y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G2x2y_S_vrr = PAZ*I_ERI_G2x2y_S_G2x2y_S_vrr+WPZ*I_ERI_G2x2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G2x2y_S_vrr = PAY*I_ERI_G2x2z_S_G2x2y_S_vrr+WPY*I_ERI_G2x2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_G2x2z_S_F2xy_S_M1_vrr;
      Double I_ERI_H2x3z_S_G2x2y_S_vrr = PAX*I_ERI_Gx3z_S_G2x2y_S_vrr+WPX*I_ERI_Gx3z_S_G2x2y_S_M1_vrr+oned2z*I_ERI_F3z_S_G2x2y_S_vrr-rhod2zsq*I_ERI_F3z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Hx4y_S_G2x2y_S_vrr = PAX*I_ERI_G4y_S_G2x2y_S_vrr+WPX*I_ERI_G4y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_G4y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G2x2y_S_vrr = PAZ*I_ERI_Gx3y_S_G2x2y_S_vrr+WPZ*I_ERI_Gx3y_S_G2x2y_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G2x2y_S_vrr = PAX*I_ERI_G2y2z_S_G2x2y_S_vrr+WPX*I_ERI_G2y2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_G2y2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G2x2y_S_vrr = PAY*I_ERI_Gx3z_S_G2x2y_S_vrr+WPY*I_ERI_Gx3z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_F2xy_S_M1_vrr;
      Double I_ERI_Hx4z_S_G2x2y_S_vrr = PAX*I_ERI_G4z_S_G2x2y_S_vrr+WPX*I_ERI_G4z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_G4z_S_Fx2y_S_M1_vrr;
      Double I_ERI_H5y_S_G2x2y_S_vrr = PAY*I_ERI_G4y_S_G2x2y_S_vrr+WPY*I_ERI_G4y_S_G2x2y_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G2x2y_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_G4y_S_F2xy_S_M1_vrr;
      Double I_ERI_H4yz_S_G2x2y_S_vrr = PAZ*I_ERI_G4y_S_G2x2y_S_vrr+WPZ*I_ERI_G4y_S_G2x2y_S_M1_vrr;
      Double I_ERI_H3y2z_S_G2x2y_S_vrr = PAZ*I_ERI_G3yz_S_G2x2y_S_vrr+WPZ*I_ERI_G3yz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_F3y_S_G2x2y_S_vrr-rhod2zsq*I_ERI_F3y_S_G2x2y_S_M1_vrr;
      Double I_ERI_H2y3z_S_G2x2y_S_vrr = PAY*I_ERI_Gy3z_S_G2x2y_S_vrr+WPY*I_ERI_Gy3z_S_G2x2y_S_M1_vrr+oned2z*I_ERI_F3z_S_G2x2y_S_vrr-rhod2zsq*I_ERI_F3z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Gy3z_S_F2xy_S_M1_vrr;
      Double I_ERI_Hy4z_S_G2x2y_S_vrr = PAY*I_ERI_G4z_S_G2x2y_S_vrr+WPY*I_ERI_G4z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_G4z_S_F2xy_S_M1_vrr;
      Double I_ERI_H5z_S_G2x2y_S_vrr = PAZ*I_ERI_G4z_S_G2x2y_S_vrr+WPZ*I_ERI_G4z_S_G2x2y_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G2x2y_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G2x2y_S_M1_vrr;
      Double I_ERI_H5x_S_G2xyz_S_vrr = PAX*I_ERI_G4x_S_G2xyz_S_vrr+WPX*I_ERI_G4x_S_G2xyz_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G2xyz_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_G4x_S_Fxyz_S_M1_vrr;
      Double I_ERI_H4xy_S_G2xyz_S_vrr = PAY*I_ERI_G4x_S_G2xyz_S_vrr+WPY*I_ERI_G4x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G4x_S_F2xz_S_M1_vrr;
      Double I_ERI_H4xz_S_G2xyz_S_vrr = PAZ*I_ERI_G4x_S_G2xyz_S_vrr+WPZ*I_ERI_G4x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G4x_S_F2xy_S_M1_vrr;
      Double I_ERI_H3x2y_S_G2xyz_S_vrr = PAY*I_ERI_G3xy_S_G2xyz_S_vrr+WPY*I_ERI_G3xy_S_G2xyz_S_M1_vrr+oned2z*I_ERI_F3x_S_G2xyz_S_vrr-rhod2zsq*I_ERI_F3x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G3xy_S_F2xz_S_M1_vrr;
      Double I_ERI_H3xyz_S_G2xyz_S_vrr = PAZ*I_ERI_G3xy_S_G2xyz_S_vrr+WPZ*I_ERI_G3xy_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G3xy_S_F2xy_S_M1_vrr;
      Double I_ERI_H3x2z_S_G2xyz_S_vrr = PAZ*I_ERI_G3xz_S_G2xyz_S_vrr+WPZ*I_ERI_G3xz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_F3x_S_G2xyz_S_vrr-rhod2zsq*I_ERI_F3x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G3xz_S_F2xy_S_M1_vrr;
      Double I_ERI_H2x3y_S_G2xyz_S_vrr = PAX*I_ERI_Gx3y_S_G2xyz_S_vrr+WPX*I_ERI_Gx3y_S_G2xyz_S_M1_vrr+oned2z*I_ERI_F3y_S_G2xyz_S_vrr-rhod2zsq*I_ERI_F3y_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_Fxyz_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G2xyz_S_vrr = PAZ*I_ERI_G2x2y_S_G2xyz_S_vrr+WPZ*I_ERI_G2x2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G2x2y_S_F2xy_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G2xyz_S_vrr = PAY*I_ERI_G2x2z_S_G2xyz_S_vrr+WPY*I_ERI_G2x2z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G2x2z_S_F2xz_S_M1_vrr;
      Double I_ERI_H2x3z_S_G2xyz_S_vrr = PAX*I_ERI_Gx3z_S_G2xyz_S_vrr+WPX*I_ERI_Gx3z_S_G2xyz_S_M1_vrr+oned2z*I_ERI_F3z_S_G2xyz_S_vrr-rhod2zsq*I_ERI_F3z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Hx4y_S_G2xyz_S_vrr = PAX*I_ERI_G4y_S_G2xyz_S_vrr+WPX*I_ERI_G4y_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_G4y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G2xyz_S_vrr = PAZ*I_ERI_Gx3y_S_G2xyz_S_vrr+WPZ*I_ERI_Gx3y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_F2xy_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G2xyz_S_vrr = PAX*I_ERI_G2y2z_S_G2xyz_S_vrr+WPX*I_ERI_G2y2z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_G2y2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G2xyz_S_vrr = PAY*I_ERI_Gx3z_S_G2xyz_S_vrr+WPY*I_ERI_Gx3z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Gx3z_S_F2xz_S_M1_vrr;
      Double I_ERI_Hx4z_S_G2xyz_S_vrr = PAX*I_ERI_G4z_S_G2xyz_S_vrr+WPX*I_ERI_G4z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_G4z_S_Fxyz_S_M1_vrr;
      Double I_ERI_H5y_S_G2xyz_S_vrr = PAY*I_ERI_G4y_S_G2xyz_S_vrr+WPY*I_ERI_G4y_S_G2xyz_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G2xyz_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G4y_S_F2xz_S_M1_vrr;
      Double I_ERI_H4yz_S_G2xyz_S_vrr = PAZ*I_ERI_G4y_S_G2xyz_S_vrr+WPZ*I_ERI_G4y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G4y_S_F2xy_S_M1_vrr;
      Double I_ERI_H3y2z_S_G2xyz_S_vrr = PAZ*I_ERI_G3yz_S_G2xyz_S_vrr+WPZ*I_ERI_G3yz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_F3y_S_G2xyz_S_vrr-rhod2zsq*I_ERI_F3y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G3yz_S_F2xy_S_M1_vrr;
      Double I_ERI_H2y3z_S_G2xyz_S_vrr = PAY*I_ERI_Gy3z_S_G2xyz_S_vrr+WPY*I_ERI_Gy3z_S_G2xyz_S_M1_vrr+oned2z*I_ERI_F3z_S_G2xyz_S_vrr-rhod2zsq*I_ERI_F3z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Gy3z_S_F2xz_S_M1_vrr;
      Double I_ERI_Hy4z_S_G2xyz_S_vrr = PAY*I_ERI_G4z_S_G2xyz_S_vrr+WPY*I_ERI_G4z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G4z_S_F2xz_S_M1_vrr;
      Double I_ERI_H5z_S_G2xyz_S_vrr = PAZ*I_ERI_G4z_S_G2xyz_S_vrr+WPZ*I_ERI_G4z_S_G2xyz_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G2xyz_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_G4z_S_F2xy_S_M1_vrr;
      Double I_ERI_H5x_S_G2x2z_S_vrr = PAX*I_ERI_G4x_S_G2x2z_S_vrr+WPX*I_ERI_G4x_S_G2x2z_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G2x2z_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G4x_S_Fx2z_S_M1_vrr;
      Double I_ERI_H4xy_S_G2x2z_S_vrr = PAY*I_ERI_G4x_S_G2x2z_S_vrr+WPY*I_ERI_G4x_S_G2x2z_S_M1_vrr;
      Double I_ERI_H4xz_S_G2x2z_S_vrr = PAZ*I_ERI_G4x_S_G2x2z_S_vrr+WPZ*I_ERI_G4x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G4x_S_F2xz_S_M1_vrr;
      Double I_ERI_H3x2y_S_G2x2z_S_vrr = PAY*I_ERI_G3xy_S_G2x2z_S_vrr+WPY*I_ERI_G3xy_S_G2x2z_S_M1_vrr+oned2z*I_ERI_F3x_S_G2x2z_S_vrr-rhod2zsq*I_ERI_F3x_S_G2x2z_S_M1_vrr;
      Double I_ERI_H3xyz_S_G2x2z_S_vrr = PAZ*I_ERI_G3xy_S_G2x2z_S_vrr+WPZ*I_ERI_G3xy_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_F2xz_S_M1_vrr;
      Double I_ERI_H3x2z_S_G2x2z_S_vrr = PAZ*I_ERI_G3xz_S_G2x2z_S_vrr+WPZ*I_ERI_G3xz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_F3x_S_G2x2z_S_vrr-rhod2zsq*I_ERI_F3x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_F2xz_S_M1_vrr;
      Double I_ERI_H2x3y_S_G2x2z_S_vrr = PAX*I_ERI_Gx3y_S_G2x2z_S_vrr+WPX*I_ERI_Gx3y_S_G2x2z_S_M1_vrr+oned2z*I_ERI_F3y_S_G2x2z_S_vrr-rhod2zsq*I_ERI_F3y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_Fx2z_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G2x2z_S_vrr = PAZ*I_ERI_G2x2y_S_G2x2z_S_vrr+WPZ*I_ERI_G2x2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G2x2y_S_F2xz_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G2x2z_S_vrr = PAY*I_ERI_G2x2z_S_G2x2z_S_vrr+WPY*I_ERI_G2x2z_S_G2x2z_S_M1_vrr;
      Double I_ERI_H2x3z_S_G2x2z_S_vrr = PAX*I_ERI_Gx3z_S_G2x2z_S_vrr+WPX*I_ERI_Gx3z_S_G2x2z_S_M1_vrr+oned2z*I_ERI_F3z_S_G2x2z_S_vrr-rhod2zsq*I_ERI_F3z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Hx4y_S_G2x2z_S_vrr = PAX*I_ERI_G4y_S_G2x2z_S_vrr+WPX*I_ERI_G4y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G4y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G2x2z_S_vrr = PAZ*I_ERI_Gx3y_S_G2x2z_S_vrr+WPZ*I_ERI_Gx3y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_F2xz_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G2x2z_S_vrr = PAX*I_ERI_G2y2z_S_G2x2z_S_vrr+WPX*I_ERI_G2y2z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G2y2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G2x2z_S_vrr = PAY*I_ERI_Gx3z_S_G2x2z_S_vrr+WPY*I_ERI_Gx3z_S_G2x2z_S_M1_vrr;
      Double I_ERI_Hx4z_S_G2x2z_S_vrr = PAX*I_ERI_G4z_S_G2x2z_S_vrr+WPX*I_ERI_G4z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G4z_S_Fx2z_S_M1_vrr;
      Double I_ERI_H5y_S_G2x2z_S_vrr = PAY*I_ERI_G4y_S_G2x2z_S_vrr+WPY*I_ERI_G4y_S_G2x2z_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G2x2z_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G2x2z_S_M1_vrr;
      Double I_ERI_H4yz_S_G2x2z_S_vrr = PAZ*I_ERI_G4y_S_G2x2z_S_vrr+WPZ*I_ERI_G4y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G4y_S_F2xz_S_M1_vrr;
      Double I_ERI_H3y2z_S_G2x2z_S_vrr = PAZ*I_ERI_G3yz_S_G2x2z_S_vrr+WPZ*I_ERI_G3yz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_F3y_S_G2x2z_S_vrr-rhod2zsq*I_ERI_F3y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_F2xz_S_M1_vrr;
      Double I_ERI_H2y3z_S_G2x2z_S_vrr = PAY*I_ERI_Gy3z_S_G2x2z_S_vrr+WPY*I_ERI_Gy3z_S_G2x2z_S_M1_vrr+oned2z*I_ERI_F3z_S_G2x2z_S_vrr-rhod2zsq*I_ERI_F3z_S_G2x2z_S_M1_vrr;
      Double I_ERI_Hy4z_S_G2x2z_S_vrr = PAY*I_ERI_G4z_S_G2x2z_S_vrr+WPY*I_ERI_G4z_S_G2x2z_S_M1_vrr;
      Double I_ERI_H5z_S_G2x2z_S_vrr = PAZ*I_ERI_G4z_S_G2x2z_S_vrr+WPZ*I_ERI_G4z_S_G2x2z_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G2x2z_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_G4z_S_F2xz_S_M1_vrr;
      Double I_ERI_H5x_S_Gx3y_S_vrr = PAX*I_ERI_G4x_S_Gx3y_S_vrr+WPX*I_ERI_G4x_S_Gx3y_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Gx3y_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Gx3y_S_M1_vrr+oned2k*I_ERI_G4x_S_F3y_S_M1_vrr;
      Double I_ERI_H4xy_S_Gx3y_S_vrr = PAY*I_ERI_G4x_S_Gx3y_S_vrr+WPY*I_ERI_G4x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_G4x_S_Fx2y_S_M1_vrr;
      Double I_ERI_H4xz_S_Gx3y_S_vrr = PAZ*I_ERI_G4x_S_Gx3y_S_vrr+WPZ*I_ERI_G4x_S_Gx3y_S_M1_vrr;
      Double I_ERI_H3x2y_S_Gx3y_S_vrr = PAY*I_ERI_G3xy_S_Gx3y_S_vrr+WPY*I_ERI_G3xy_S_Gx3y_S_M1_vrr+oned2z*I_ERI_F3x_S_Gx3y_S_vrr-rhod2zsq*I_ERI_F3x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_G3xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_H3xyz_S_Gx3y_S_vrr = PAZ*I_ERI_G3xy_S_Gx3y_S_vrr+WPZ*I_ERI_G3xy_S_Gx3y_S_M1_vrr;
      Double I_ERI_H3x2z_S_Gx3y_S_vrr = PAZ*I_ERI_G3xz_S_Gx3y_S_vrr+WPZ*I_ERI_G3xz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_F3x_S_Gx3y_S_vrr-rhod2zsq*I_ERI_F3x_S_Gx3y_S_M1_vrr;
      Double I_ERI_H2x3y_S_Gx3y_S_vrr = PAX*I_ERI_Gx3y_S_Gx3y_S_vrr+WPX*I_ERI_Gx3y_S_Gx3y_S_M1_vrr+oned2z*I_ERI_F3y_S_Gx3y_S_vrr-rhod2zsq*I_ERI_F3y_S_Gx3y_S_M1_vrr+oned2k*I_ERI_Gx3y_S_F3y_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Gx3y_S_vrr = PAZ*I_ERI_G2x2y_S_Gx3y_S_vrr+WPZ*I_ERI_G2x2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Gx3y_S_vrr = PAY*I_ERI_G2x2z_S_Gx3y_S_vrr+WPY*I_ERI_G2x2z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_H2x3z_S_Gx3y_S_vrr = PAX*I_ERI_Gx3z_S_Gx3y_S_vrr+WPX*I_ERI_Gx3z_S_Gx3y_S_M1_vrr+oned2z*I_ERI_F3z_S_Gx3y_S_vrr-rhod2zsq*I_ERI_F3z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_Gx3z_S_F3y_S_M1_vrr;
      Double I_ERI_Hx4y_S_Gx3y_S_vrr = PAX*I_ERI_G4y_S_Gx3y_S_vrr+WPX*I_ERI_G4y_S_Gx3y_S_M1_vrr+oned2k*I_ERI_G4y_S_F3y_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Gx3y_S_vrr = PAZ*I_ERI_Gx3y_S_Gx3y_S_vrr+WPZ*I_ERI_Gx3y_S_Gx3y_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Gx3y_S_vrr = PAX*I_ERI_G2y2z_S_Gx3y_S_vrr+WPX*I_ERI_G2y2z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_G2y2z_S_F3y_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Gx3y_S_vrr = PAY*I_ERI_Gx3z_S_Gx3y_S_vrr+WPY*I_ERI_Gx3z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Gx3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Hx4z_S_Gx3y_S_vrr = PAX*I_ERI_G4z_S_Gx3y_S_vrr+WPX*I_ERI_G4z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_G4z_S_F3y_S_M1_vrr;
      Double I_ERI_H5y_S_Gx3y_S_vrr = PAY*I_ERI_G4y_S_Gx3y_S_vrr+WPY*I_ERI_G4y_S_Gx3y_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Gx3y_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_G4y_S_Fx2y_S_M1_vrr;
      Double I_ERI_H4yz_S_Gx3y_S_vrr = PAZ*I_ERI_G4y_S_Gx3y_S_vrr+WPZ*I_ERI_G4y_S_Gx3y_S_M1_vrr;
      Double I_ERI_H3y2z_S_Gx3y_S_vrr = PAZ*I_ERI_G3yz_S_Gx3y_S_vrr+WPZ*I_ERI_G3yz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_F3y_S_Gx3y_S_vrr-rhod2zsq*I_ERI_F3y_S_Gx3y_S_M1_vrr;
      Double I_ERI_H2y3z_S_Gx3y_S_vrr = PAY*I_ERI_Gy3z_S_Gx3y_S_vrr+WPY*I_ERI_Gy3z_S_Gx3y_S_M1_vrr+oned2z*I_ERI_F3z_S_Gx3y_S_vrr-rhod2zsq*I_ERI_F3z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Gy3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Hy4z_S_Gx3y_S_vrr = PAY*I_ERI_G4z_S_Gx3y_S_vrr+WPY*I_ERI_G4z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_G4z_S_Fx2y_S_M1_vrr;
      Double I_ERI_H5z_S_Gx3y_S_vrr = PAZ*I_ERI_G4z_S_Gx3y_S_vrr+WPZ*I_ERI_G4z_S_Gx3y_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Gx3y_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Gx3y_S_M1_vrr;
      Double I_ERI_H5x_S_Gx2yz_S_vrr = PAX*I_ERI_G4x_S_Gx2yz_S_vrr+WPX*I_ERI_G4x_S_Gx2yz_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Gx2yz_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G4x_S_F2yz_S_M1_vrr;
      Double I_ERI_H4xy_S_Gx2yz_S_vrr = PAY*I_ERI_G4x_S_Gx2yz_S_vrr+WPY*I_ERI_G4x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_G4x_S_Fxyz_S_M1_vrr;
      Double I_ERI_H4xz_S_Gx2yz_S_vrr = PAZ*I_ERI_G4x_S_Gx2yz_S_vrr+WPZ*I_ERI_G4x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G4x_S_Fx2y_S_M1_vrr;
      Double I_ERI_H3x2y_S_Gx2yz_S_vrr = PAY*I_ERI_G3xy_S_Gx2yz_S_vrr+WPY*I_ERI_G3xy_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_F3x_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_F3x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_Fxyz_S_M1_vrr;
      Double I_ERI_H3xyz_S_Gx2yz_S_vrr = PAZ*I_ERI_G3xy_S_Gx2yz_S_vrr+WPZ*I_ERI_G3xy_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G3xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_H3x2z_S_Gx2yz_S_vrr = PAZ*I_ERI_G3xz_S_Gx2yz_S_vrr+WPZ*I_ERI_G3xz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_F3x_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_F3x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G3xz_S_Fx2y_S_M1_vrr;
      Double I_ERI_H2x3y_S_Gx2yz_S_vrr = PAX*I_ERI_Gx3y_S_Gx2yz_S_vrr+WPX*I_ERI_Gx3y_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_F3y_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_F3y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_F2yz_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Gx2yz_S_vrr = PAZ*I_ERI_G2x2y_S_Gx2yz_S_vrr+WPZ*I_ERI_G2x2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G2x2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Gx2yz_S_vrr = PAY*I_ERI_G2x2z_S_Gx2yz_S_vrr+WPY*I_ERI_G2x2z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_G2x2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_H2x3z_S_Gx2yz_S_vrr = PAX*I_ERI_Gx3z_S_Gx2yz_S_vrr+WPX*I_ERI_Gx3z_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_F3z_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_F3z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Gx3z_S_F2yz_S_M1_vrr;
      Double I_ERI_Hx4y_S_Gx2yz_S_vrr = PAX*I_ERI_G4y_S_Gx2yz_S_vrr+WPX*I_ERI_G4y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G4y_S_F2yz_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Gx2yz_S_vrr = PAZ*I_ERI_Gx3y_S_Gx2yz_S_vrr+WPZ*I_ERI_Gx3y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Gx2yz_S_vrr = PAX*I_ERI_G2y2z_S_Gx2yz_S_vrr+WPX*I_ERI_G2y2z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G2y2z_S_F2yz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Gx2yz_S_vrr = PAY*I_ERI_Gx3z_S_Gx2yz_S_vrr+WPY*I_ERI_Gx3z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Hx4z_S_Gx2yz_S_vrr = PAX*I_ERI_G4z_S_Gx2yz_S_vrr+WPX*I_ERI_G4z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G4z_S_F2yz_S_M1_vrr;
      Double I_ERI_H5y_S_Gx2yz_S_vrr = PAY*I_ERI_G4y_S_Gx2yz_S_vrr+WPY*I_ERI_G4y_S_Gx2yz_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Gx2yz_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_G4y_S_Fxyz_S_M1_vrr;
      Double I_ERI_H4yz_S_Gx2yz_S_vrr = PAZ*I_ERI_G4y_S_Gx2yz_S_vrr+WPZ*I_ERI_G4y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G4y_S_Fx2y_S_M1_vrr;
      Double I_ERI_H3y2z_S_Gx2yz_S_vrr = PAZ*I_ERI_G3yz_S_Gx2yz_S_vrr+WPZ*I_ERI_G3yz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_F3y_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_F3y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G3yz_S_Fx2y_S_M1_vrr;
      Double I_ERI_H2y3z_S_Gx2yz_S_vrr = PAY*I_ERI_Gy3z_S_Gx2yz_S_vrr+WPY*I_ERI_Gy3z_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_F3z_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_F3z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Gy3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Hy4z_S_Gx2yz_S_vrr = PAY*I_ERI_G4z_S_Gx2yz_S_vrr+WPY*I_ERI_G4z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_G4z_S_Fxyz_S_M1_vrr;
      Double I_ERI_H5z_S_Gx2yz_S_vrr = PAZ*I_ERI_G4z_S_Gx2yz_S_vrr+WPZ*I_ERI_G4z_S_Gx2yz_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Gx2yz_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_G4z_S_Fx2y_S_M1_vrr;
      Double I_ERI_H5x_S_Gxy2z_S_vrr = PAX*I_ERI_G4x_S_Gxy2z_S_vrr+WPX*I_ERI_G4x_S_Gxy2z_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Gxy2z_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_G4x_S_Fy2z_S_M1_vrr;
      Double I_ERI_H4xy_S_Gxy2z_S_vrr = PAY*I_ERI_G4x_S_Gxy2z_S_vrr+WPY*I_ERI_G4x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_G4x_S_Fx2z_S_M1_vrr;
      Double I_ERI_H4xz_S_Gxy2z_S_vrr = PAZ*I_ERI_G4x_S_Gxy2z_S_vrr+WPZ*I_ERI_G4x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_G4x_S_Fxyz_S_M1_vrr;
      Double I_ERI_H3x2y_S_Gxy2z_S_vrr = PAY*I_ERI_G3xy_S_Gxy2z_S_vrr+WPY*I_ERI_G3xy_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_F3x_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_F3x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_G3xy_S_Fx2z_S_M1_vrr;
      Double I_ERI_H3xyz_S_Gxy2z_S_vrr = PAZ*I_ERI_G3xy_S_Gxy2z_S_vrr+WPZ*I_ERI_G3xy_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_Fxyz_S_M1_vrr;
      Double I_ERI_H3x2z_S_Gxy2z_S_vrr = PAZ*I_ERI_G3xz_S_Gxy2z_S_vrr+WPZ*I_ERI_G3xz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_F3x_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_F3x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_Fxyz_S_M1_vrr;
      Double I_ERI_H2x3y_S_Gxy2z_S_vrr = PAX*I_ERI_Gx3y_S_Gxy2z_S_vrr+WPX*I_ERI_Gx3y_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_F3y_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_F3y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Gxy2z_S_vrr = PAZ*I_ERI_G2x2y_S_Gxy2z_S_vrr+WPZ*I_ERI_G2x2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_G2x2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Gxy2z_S_vrr = PAY*I_ERI_G2x2z_S_Gxy2z_S_vrr+WPY*I_ERI_G2x2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_G2x2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_H2x3z_S_Gxy2z_S_vrr = PAX*I_ERI_Gx3z_S_Gxy2z_S_vrr+WPX*I_ERI_Gx3z_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_F3z_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_F3z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Hx4y_S_Gxy2z_S_vrr = PAX*I_ERI_G4y_S_Gxy2z_S_vrr+WPX*I_ERI_G4y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_G4y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Gxy2z_S_vrr = PAZ*I_ERI_Gx3y_S_Gxy2z_S_vrr+WPZ*I_ERI_Gx3y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Gxy2z_S_vrr = PAX*I_ERI_G2y2z_S_Gxy2z_S_vrr+WPX*I_ERI_G2y2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_G2y2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Gxy2z_S_vrr = PAY*I_ERI_Gx3z_S_Gxy2z_S_vrr+WPY*I_ERI_Gx3z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Hx4z_S_Gxy2z_S_vrr = PAX*I_ERI_G4z_S_Gxy2z_S_vrr+WPX*I_ERI_G4z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_G4z_S_Fy2z_S_M1_vrr;
      Double I_ERI_H5y_S_Gxy2z_S_vrr = PAY*I_ERI_G4y_S_Gxy2z_S_vrr+WPY*I_ERI_G4y_S_Gxy2z_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Gxy2z_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_G4y_S_Fx2z_S_M1_vrr;
      Double I_ERI_H4yz_S_Gxy2z_S_vrr = PAZ*I_ERI_G4y_S_Gxy2z_S_vrr+WPZ*I_ERI_G4y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_G4y_S_Fxyz_S_M1_vrr;
      Double I_ERI_H3y2z_S_Gxy2z_S_vrr = PAZ*I_ERI_G3yz_S_Gxy2z_S_vrr+WPZ*I_ERI_G3yz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_F3y_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_F3y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_Fxyz_S_M1_vrr;
      Double I_ERI_H2y3z_S_Gxy2z_S_vrr = PAY*I_ERI_Gy3z_S_Gxy2z_S_vrr+WPY*I_ERI_Gy3z_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_F3z_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_F3z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Gy3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Hy4z_S_Gxy2z_S_vrr = PAY*I_ERI_G4z_S_Gxy2z_S_vrr+WPY*I_ERI_G4z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_G4z_S_Fx2z_S_M1_vrr;
      Double I_ERI_H5z_S_Gxy2z_S_vrr = PAZ*I_ERI_G4z_S_Gxy2z_S_vrr+WPZ*I_ERI_G4z_S_Gxy2z_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Gxy2z_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_G4z_S_Fxyz_S_M1_vrr;
      Double I_ERI_H5x_S_Gx3z_S_vrr = PAX*I_ERI_G4x_S_Gx3z_S_vrr+WPX*I_ERI_G4x_S_Gx3z_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Gx3z_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Gx3z_S_M1_vrr+oned2k*I_ERI_G4x_S_F3z_S_M1_vrr;
      Double I_ERI_H4xy_S_Gx3z_S_vrr = PAY*I_ERI_G4x_S_Gx3z_S_vrr+WPY*I_ERI_G4x_S_Gx3z_S_M1_vrr;
      Double I_ERI_H4xz_S_Gx3z_S_vrr = PAZ*I_ERI_G4x_S_Gx3z_S_vrr+WPZ*I_ERI_G4x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_G4x_S_Fx2z_S_M1_vrr;
      Double I_ERI_H3x2y_S_Gx3z_S_vrr = PAY*I_ERI_G3xy_S_Gx3z_S_vrr+WPY*I_ERI_G3xy_S_Gx3z_S_M1_vrr+oned2z*I_ERI_F3x_S_Gx3z_S_vrr-rhod2zsq*I_ERI_F3x_S_Gx3z_S_M1_vrr;
      Double I_ERI_H3xyz_S_Gx3z_S_vrr = PAZ*I_ERI_G3xy_S_Gx3z_S_vrr+WPZ*I_ERI_G3xy_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_G3xy_S_Fx2z_S_M1_vrr;
      Double I_ERI_H3x2z_S_Gx3z_S_vrr = PAZ*I_ERI_G3xz_S_Gx3z_S_vrr+WPZ*I_ERI_G3xz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_F3x_S_Gx3z_S_vrr-rhod2zsq*I_ERI_F3x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_G3xz_S_Fx2z_S_M1_vrr;
      Double I_ERI_H2x3y_S_Gx3z_S_vrr = PAX*I_ERI_Gx3y_S_Gx3z_S_vrr+WPX*I_ERI_Gx3y_S_Gx3z_S_M1_vrr+oned2z*I_ERI_F3y_S_Gx3z_S_vrr-rhod2zsq*I_ERI_F3y_S_Gx3z_S_M1_vrr+oned2k*I_ERI_Gx3y_S_F3z_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Gx3z_S_vrr = PAZ*I_ERI_G2x2y_S_Gx3z_S_vrr+WPZ*I_ERI_G2x2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Gx3z_S_vrr = PAY*I_ERI_G2x2z_S_Gx3z_S_vrr+WPY*I_ERI_G2x2z_S_Gx3z_S_M1_vrr;
      Double I_ERI_H2x3z_S_Gx3z_S_vrr = PAX*I_ERI_Gx3z_S_Gx3z_S_vrr+WPX*I_ERI_Gx3z_S_Gx3z_S_M1_vrr+oned2z*I_ERI_F3z_S_Gx3z_S_vrr-rhod2zsq*I_ERI_F3z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_Gx3z_S_F3z_S_M1_vrr;
      Double I_ERI_Hx4y_S_Gx3z_S_vrr = PAX*I_ERI_G4y_S_Gx3z_S_vrr+WPX*I_ERI_G4y_S_Gx3z_S_M1_vrr+oned2k*I_ERI_G4y_S_F3z_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Gx3z_S_vrr = PAZ*I_ERI_Gx3y_S_Gx3z_S_vrr+WPZ*I_ERI_Gx3y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Gx3y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Gx3z_S_vrr = PAX*I_ERI_G2y2z_S_Gx3z_S_vrr+WPX*I_ERI_G2y2z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_G2y2z_S_F3z_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Gx3z_S_vrr = PAY*I_ERI_Gx3z_S_Gx3z_S_vrr+WPY*I_ERI_Gx3z_S_Gx3z_S_M1_vrr;
      Double I_ERI_Hx4z_S_Gx3z_S_vrr = PAX*I_ERI_G4z_S_Gx3z_S_vrr+WPX*I_ERI_G4z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_G4z_S_F3z_S_M1_vrr;
      Double I_ERI_H5y_S_Gx3z_S_vrr = PAY*I_ERI_G4y_S_Gx3z_S_vrr+WPY*I_ERI_G4y_S_Gx3z_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Gx3z_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Gx3z_S_M1_vrr;
      Double I_ERI_H4yz_S_Gx3z_S_vrr = PAZ*I_ERI_G4y_S_Gx3z_S_vrr+WPZ*I_ERI_G4y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_G4y_S_Fx2z_S_M1_vrr;
      Double I_ERI_H3y2z_S_Gx3z_S_vrr = PAZ*I_ERI_G3yz_S_Gx3z_S_vrr+WPZ*I_ERI_G3yz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_F3y_S_Gx3z_S_vrr-rhod2zsq*I_ERI_F3y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_G3yz_S_Fx2z_S_M1_vrr;
      Double I_ERI_H2y3z_S_Gx3z_S_vrr = PAY*I_ERI_Gy3z_S_Gx3z_S_vrr+WPY*I_ERI_Gy3z_S_Gx3z_S_M1_vrr+oned2z*I_ERI_F3z_S_Gx3z_S_vrr-rhod2zsq*I_ERI_F3z_S_Gx3z_S_M1_vrr;
      Double I_ERI_Hy4z_S_Gx3z_S_vrr = PAY*I_ERI_G4z_S_Gx3z_S_vrr+WPY*I_ERI_G4z_S_Gx3z_S_M1_vrr;
      Double I_ERI_H5z_S_Gx3z_S_vrr = PAZ*I_ERI_G4z_S_Gx3z_S_vrr+WPZ*I_ERI_G4z_S_Gx3z_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Gx3z_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_G4z_S_Fx2z_S_M1_vrr;
      Double I_ERI_H5x_S_G4y_S_vrr = PAX*I_ERI_G4x_S_G4y_S_vrr+WPX*I_ERI_G4x_S_G4y_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G4y_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G4y_S_M1_vrr;
      Double I_ERI_H4xy_S_G4y_S_vrr = PAY*I_ERI_G4x_S_G4y_S_vrr+WPY*I_ERI_G4x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_G4x_S_F3y_S_M1_vrr;
      Double I_ERI_H4xz_S_G4y_S_vrr = PAZ*I_ERI_G4x_S_G4y_S_vrr+WPZ*I_ERI_G4x_S_G4y_S_M1_vrr;
      Double I_ERI_H3x2y_S_G4y_S_vrr = PAY*I_ERI_G3xy_S_G4y_S_vrr+WPY*I_ERI_G3xy_S_G4y_S_M1_vrr+oned2z*I_ERI_F3x_S_G4y_S_vrr-rhod2zsq*I_ERI_F3x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_G3xy_S_F3y_S_M1_vrr;
      Double I_ERI_H3xyz_S_G4y_S_vrr = PAZ*I_ERI_G3xy_S_G4y_S_vrr+WPZ*I_ERI_G3xy_S_G4y_S_M1_vrr;
      Double I_ERI_H3x2z_S_G4y_S_vrr = PAZ*I_ERI_G3xz_S_G4y_S_vrr+WPZ*I_ERI_G3xz_S_G4y_S_M1_vrr+oned2z*I_ERI_F3x_S_G4y_S_vrr-rhod2zsq*I_ERI_F3x_S_G4y_S_M1_vrr;
      Double I_ERI_H2x3y_S_G4y_S_vrr = PAX*I_ERI_Gx3y_S_G4y_S_vrr+WPX*I_ERI_Gx3y_S_G4y_S_M1_vrr+oned2z*I_ERI_F3y_S_G4y_S_vrr-rhod2zsq*I_ERI_F3y_S_G4y_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G4y_S_vrr = PAZ*I_ERI_G2x2y_S_G4y_S_vrr+WPZ*I_ERI_G2x2y_S_G4y_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G4y_S_vrr = PAY*I_ERI_G2x2z_S_G4y_S_vrr+WPY*I_ERI_G2x2z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_G2x2z_S_F3y_S_M1_vrr;
      Double I_ERI_H2x3z_S_G4y_S_vrr = PAX*I_ERI_Gx3z_S_G4y_S_vrr+WPX*I_ERI_Gx3z_S_G4y_S_M1_vrr+oned2z*I_ERI_F3z_S_G4y_S_vrr-rhod2zsq*I_ERI_F3z_S_G4y_S_M1_vrr;
      Double I_ERI_Hx4y_S_G4y_S_vrr = PAX*I_ERI_G4y_S_G4y_S_vrr+WPX*I_ERI_G4y_S_G4y_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G4y_S_vrr = PAZ*I_ERI_Gx3y_S_G4y_S_vrr+WPZ*I_ERI_Gx3y_S_G4y_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G4y_S_vrr = PAX*I_ERI_G2y2z_S_G4y_S_vrr+WPX*I_ERI_G2y2z_S_G4y_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G4y_S_vrr = PAY*I_ERI_Gx3z_S_G4y_S_vrr+WPY*I_ERI_Gx3z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Gx3z_S_F3y_S_M1_vrr;
      Double I_ERI_Hx4z_S_G4y_S_vrr = PAX*I_ERI_G4z_S_G4y_S_vrr+WPX*I_ERI_G4z_S_G4y_S_M1_vrr;
      Double I_ERI_H5y_S_G4y_S_vrr = PAY*I_ERI_G4y_S_G4y_S_vrr+WPY*I_ERI_G4y_S_G4y_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G4y_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G4y_S_M1_vrr+4*oned2k*I_ERI_G4y_S_F3y_S_M1_vrr;
      Double I_ERI_H4yz_S_G4y_S_vrr = PAZ*I_ERI_G4y_S_G4y_S_vrr+WPZ*I_ERI_G4y_S_G4y_S_M1_vrr;
      Double I_ERI_H3y2z_S_G4y_S_vrr = PAZ*I_ERI_G3yz_S_G4y_S_vrr+WPZ*I_ERI_G3yz_S_G4y_S_M1_vrr+oned2z*I_ERI_F3y_S_G4y_S_vrr-rhod2zsq*I_ERI_F3y_S_G4y_S_M1_vrr;
      Double I_ERI_H2y3z_S_G4y_S_vrr = PAY*I_ERI_Gy3z_S_G4y_S_vrr+WPY*I_ERI_Gy3z_S_G4y_S_M1_vrr+oned2z*I_ERI_F3z_S_G4y_S_vrr-rhod2zsq*I_ERI_F3z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Gy3z_S_F3y_S_M1_vrr;
      Double I_ERI_Hy4z_S_G4y_S_vrr = PAY*I_ERI_G4z_S_G4y_S_vrr+WPY*I_ERI_G4z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_G4z_S_F3y_S_M1_vrr;
      Double I_ERI_H5z_S_G4y_S_vrr = PAZ*I_ERI_G4z_S_G4y_S_vrr+WPZ*I_ERI_G4z_S_G4y_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G4y_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G4y_S_M1_vrr;
      Double I_ERI_H5x_S_G3yz_S_vrr = PAX*I_ERI_G4x_S_G3yz_S_vrr+WPX*I_ERI_G4x_S_G3yz_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G3yz_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G3yz_S_M1_vrr;
      Double I_ERI_H4xy_S_G3yz_S_vrr = PAY*I_ERI_G4x_S_G3yz_S_vrr+WPY*I_ERI_G4x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_G4x_S_F2yz_S_M1_vrr;
      Double I_ERI_H4xz_S_G3yz_S_vrr = PAZ*I_ERI_G4x_S_G3yz_S_vrr+WPZ*I_ERI_G4x_S_G3yz_S_M1_vrr+oned2k*I_ERI_G4x_S_F3y_S_M1_vrr;
      Double I_ERI_H3x2y_S_G3yz_S_vrr = PAY*I_ERI_G3xy_S_G3yz_S_vrr+WPY*I_ERI_G3xy_S_G3yz_S_M1_vrr+oned2z*I_ERI_F3x_S_G3yz_S_vrr-rhod2zsq*I_ERI_F3x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_G3xy_S_F2yz_S_M1_vrr;
      Double I_ERI_H3xyz_S_G3yz_S_vrr = PAZ*I_ERI_G3xy_S_G3yz_S_vrr+WPZ*I_ERI_G3xy_S_G3yz_S_M1_vrr+oned2k*I_ERI_G3xy_S_F3y_S_M1_vrr;
      Double I_ERI_H3x2z_S_G3yz_S_vrr = PAZ*I_ERI_G3xz_S_G3yz_S_vrr+WPZ*I_ERI_G3xz_S_G3yz_S_M1_vrr+oned2z*I_ERI_F3x_S_G3yz_S_vrr-rhod2zsq*I_ERI_F3x_S_G3yz_S_M1_vrr+oned2k*I_ERI_G3xz_S_F3y_S_M1_vrr;
      Double I_ERI_H2x3y_S_G3yz_S_vrr = PAX*I_ERI_Gx3y_S_G3yz_S_vrr+WPX*I_ERI_Gx3y_S_G3yz_S_M1_vrr+oned2z*I_ERI_F3y_S_G3yz_S_vrr-rhod2zsq*I_ERI_F3y_S_G3yz_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G3yz_S_vrr = PAZ*I_ERI_G2x2y_S_G3yz_S_vrr+WPZ*I_ERI_G2x2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_G2x2y_S_F3y_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G3yz_S_vrr = PAY*I_ERI_G2x2z_S_G3yz_S_vrr+WPY*I_ERI_G2x2z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_F2yz_S_M1_vrr;
      Double I_ERI_H2x3z_S_G3yz_S_vrr = PAX*I_ERI_Gx3z_S_G3yz_S_vrr+WPX*I_ERI_Gx3z_S_G3yz_S_M1_vrr+oned2z*I_ERI_F3z_S_G3yz_S_vrr-rhod2zsq*I_ERI_F3z_S_G3yz_S_M1_vrr;
      Double I_ERI_Hx4y_S_G3yz_S_vrr = PAX*I_ERI_G4y_S_G3yz_S_vrr+WPX*I_ERI_G4y_S_G3yz_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G3yz_S_vrr = PAZ*I_ERI_Gx3y_S_G3yz_S_vrr+WPZ*I_ERI_Gx3y_S_G3yz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_F3y_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G3yz_S_vrr = PAX*I_ERI_G2y2z_S_G3yz_S_vrr+WPX*I_ERI_G2y2z_S_G3yz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G3yz_S_vrr = PAY*I_ERI_Gx3z_S_G3yz_S_vrr+WPY*I_ERI_Gx3z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Gx3z_S_F2yz_S_M1_vrr;
      Double I_ERI_Hx4z_S_G3yz_S_vrr = PAX*I_ERI_G4z_S_G3yz_S_vrr+WPX*I_ERI_G4z_S_G3yz_S_M1_vrr;
      Double I_ERI_H5y_S_G3yz_S_vrr = PAY*I_ERI_G4y_S_G3yz_S_vrr+WPY*I_ERI_G4y_S_G3yz_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G3yz_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_G4y_S_F2yz_S_M1_vrr;
      Double I_ERI_H4yz_S_G3yz_S_vrr = PAZ*I_ERI_G4y_S_G3yz_S_vrr+WPZ*I_ERI_G4y_S_G3yz_S_M1_vrr+oned2k*I_ERI_G4y_S_F3y_S_M1_vrr;
      Double I_ERI_H3y2z_S_G3yz_S_vrr = PAZ*I_ERI_G3yz_S_G3yz_S_vrr+WPZ*I_ERI_G3yz_S_G3yz_S_M1_vrr+oned2z*I_ERI_F3y_S_G3yz_S_vrr-rhod2zsq*I_ERI_F3y_S_G3yz_S_M1_vrr+oned2k*I_ERI_G3yz_S_F3y_S_M1_vrr;
      Double I_ERI_H2y3z_S_G3yz_S_vrr = PAY*I_ERI_Gy3z_S_G3yz_S_vrr+WPY*I_ERI_Gy3z_S_G3yz_S_M1_vrr+oned2z*I_ERI_F3z_S_G3yz_S_vrr-rhod2zsq*I_ERI_F3z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Gy3z_S_F2yz_S_M1_vrr;
      Double I_ERI_Hy4z_S_G3yz_S_vrr = PAY*I_ERI_G4z_S_G3yz_S_vrr+WPY*I_ERI_G4z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_G4z_S_F2yz_S_M1_vrr;
      Double I_ERI_H5z_S_G3yz_S_vrr = PAZ*I_ERI_G4z_S_G3yz_S_vrr+WPZ*I_ERI_G4z_S_G3yz_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G3yz_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G3yz_S_M1_vrr+oned2k*I_ERI_G4z_S_F3y_S_M1_vrr;
      Double I_ERI_H5x_S_G2y2z_S_vrr = PAX*I_ERI_G4x_S_G2y2z_S_vrr+WPX*I_ERI_G4x_S_G2y2z_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G2y2z_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G2y2z_S_M1_vrr;
      Double I_ERI_H4xy_S_G2y2z_S_vrr = PAY*I_ERI_G4x_S_G2y2z_S_vrr+WPY*I_ERI_G4x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G4x_S_Fy2z_S_M1_vrr;
      Double I_ERI_H4xz_S_G2y2z_S_vrr = PAZ*I_ERI_G4x_S_G2y2z_S_vrr+WPZ*I_ERI_G4x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G4x_S_F2yz_S_M1_vrr;
      Double I_ERI_H3x2y_S_G2y2z_S_vrr = PAY*I_ERI_G3xy_S_G2y2z_S_vrr+WPY*I_ERI_G3xy_S_G2y2z_S_M1_vrr+oned2z*I_ERI_F3x_S_G2y2z_S_vrr-rhod2zsq*I_ERI_F3x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_Fy2z_S_M1_vrr;
      Double I_ERI_H3xyz_S_G2y2z_S_vrr = PAZ*I_ERI_G3xy_S_G2y2z_S_vrr+WPZ*I_ERI_G3xy_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_F2yz_S_M1_vrr;
      Double I_ERI_H3x2z_S_G2y2z_S_vrr = PAZ*I_ERI_G3xz_S_G2y2z_S_vrr+WPZ*I_ERI_G3xz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_F3x_S_G2y2z_S_vrr-rhod2zsq*I_ERI_F3x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_F2yz_S_M1_vrr;
      Double I_ERI_H2x3y_S_G2y2z_S_vrr = PAX*I_ERI_Gx3y_S_G2y2z_S_vrr+WPX*I_ERI_Gx3y_S_G2y2z_S_M1_vrr+oned2z*I_ERI_F3y_S_G2y2z_S_vrr-rhod2zsq*I_ERI_F3y_S_G2y2z_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G2y2z_S_vrr = PAZ*I_ERI_G2x2y_S_G2y2z_S_vrr+WPZ*I_ERI_G2x2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G2x2y_S_F2yz_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G2y2z_S_vrr = PAY*I_ERI_G2x2z_S_G2y2z_S_vrr+WPY*I_ERI_G2x2z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G2x2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_H2x3z_S_G2y2z_S_vrr = PAX*I_ERI_Gx3z_S_G2y2z_S_vrr+WPX*I_ERI_Gx3z_S_G2y2z_S_M1_vrr+oned2z*I_ERI_F3z_S_G2y2z_S_vrr-rhod2zsq*I_ERI_F3z_S_G2y2z_S_M1_vrr;
      Double I_ERI_Hx4y_S_G2y2z_S_vrr = PAX*I_ERI_G4y_S_G2y2z_S_vrr+WPX*I_ERI_G4y_S_G2y2z_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G2y2z_S_vrr = PAZ*I_ERI_Gx3y_S_G2y2z_S_vrr+WPZ*I_ERI_Gx3y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_F2yz_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G2y2z_S_vrr = PAX*I_ERI_G2y2z_S_G2y2z_S_vrr+WPX*I_ERI_G2y2z_S_G2y2z_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G2y2z_S_vrr = PAY*I_ERI_Gx3z_S_G2y2z_S_vrr+WPY*I_ERI_Gx3z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Hx4z_S_G2y2z_S_vrr = PAX*I_ERI_G4z_S_G2y2z_S_vrr+WPX*I_ERI_G4z_S_G2y2z_S_M1_vrr;
      Double I_ERI_H5y_S_G2y2z_S_vrr = PAY*I_ERI_G4y_S_G2y2z_S_vrr+WPY*I_ERI_G4y_S_G2y2z_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G2y2z_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G4y_S_Fy2z_S_M1_vrr;
      Double I_ERI_H4yz_S_G2y2z_S_vrr = PAZ*I_ERI_G4y_S_G2y2z_S_vrr+WPZ*I_ERI_G4y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G4y_S_F2yz_S_M1_vrr;
      Double I_ERI_H3y2z_S_G2y2z_S_vrr = PAZ*I_ERI_G3yz_S_G2y2z_S_vrr+WPZ*I_ERI_G3yz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_F3y_S_G2y2z_S_vrr-rhod2zsq*I_ERI_F3y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_F2yz_S_M1_vrr;
      Double I_ERI_H2y3z_S_G2y2z_S_vrr = PAY*I_ERI_Gy3z_S_G2y2z_S_vrr+WPY*I_ERI_Gy3z_S_G2y2z_S_M1_vrr+oned2z*I_ERI_F3z_S_G2y2z_S_vrr-rhod2zsq*I_ERI_F3z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Gy3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Hy4z_S_G2y2z_S_vrr = PAY*I_ERI_G4z_S_G2y2z_S_vrr+WPY*I_ERI_G4z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G4z_S_Fy2z_S_M1_vrr;
      Double I_ERI_H5z_S_G2y2z_S_vrr = PAZ*I_ERI_G4z_S_G2y2z_S_vrr+WPZ*I_ERI_G4z_S_G2y2z_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G2y2z_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_G4z_S_F2yz_S_M1_vrr;
      Double I_ERI_H5x_S_Gy3z_S_vrr = PAX*I_ERI_G4x_S_Gy3z_S_vrr+WPX*I_ERI_G4x_S_Gy3z_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Gy3z_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Gy3z_S_M1_vrr;
      Double I_ERI_H4xy_S_Gy3z_S_vrr = PAY*I_ERI_G4x_S_Gy3z_S_vrr+WPY*I_ERI_G4x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_G4x_S_F3z_S_M1_vrr;
      Double I_ERI_H4xz_S_Gy3z_S_vrr = PAZ*I_ERI_G4x_S_Gy3z_S_vrr+WPZ*I_ERI_G4x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_G4x_S_Fy2z_S_M1_vrr;
      Double I_ERI_H3x2y_S_Gy3z_S_vrr = PAY*I_ERI_G3xy_S_Gy3z_S_vrr+WPY*I_ERI_G3xy_S_Gy3z_S_M1_vrr+oned2z*I_ERI_F3x_S_Gy3z_S_vrr-rhod2zsq*I_ERI_F3x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_G3xy_S_F3z_S_M1_vrr;
      Double I_ERI_H3xyz_S_Gy3z_S_vrr = PAZ*I_ERI_G3xy_S_Gy3z_S_vrr+WPZ*I_ERI_G3xy_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_G3xy_S_Fy2z_S_M1_vrr;
      Double I_ERI_H3x2z_S_Gy3z_S_vrr = PAZ*I_ERI_G3xz_S_Gy3z_S_vrr+WPZ*I_ERI_G3xz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_F3x_S_Gy3z_S_vrr-rhod2zsq*I_ERI_F3x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_G3xz_S_Fy2z_S_M1_vrr;
      Double I_ERI_H2x3y_S_Gy3z_S_vrr = PAX*I_ERI_Gx3y_S_Gy3z_S_vrr+WPX*I_ERI_Gx3y_S_Gy3z_S_M1_vrr+oned2z*I_ERI_F3y_S_Gy3z_S_vrr-rhod2zsq*I_ERI_F3y_S_Gy3z_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Gy3z_S_vrr = PAZ*I_ERI_G2x2y_S_Gy3z_S_vrr+WPZ*I_ERI_G2x2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Gy3z_S_vrr = PAY*I_ERI_G2x2z_S_Gy3z_S_vrr+WPY*I_ERI_G2x2z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_G2x2z_S_F3z_S_M1_vrr;
      Double I_ERI_H2x3z_S_Gy3z_S_vrr = PAX*I_ERI_Gx3z_S_Gy3z_S_vrr+WPX*I_ERI_Gx3z_S_Gy3z_S_M1_vrr+oned2z*I_ERI_F3z_S_Gy3z_S_vrr-rhod2zsq*I_ERI_F3z_S_Gy3z_S_M1_vrr;
      Double I_ERI_Hx4y_S_Gy3z_S_vrr = PAX*I_ERI_G4y_S_Gy3z_S_vrr+WPX*I_ERI_G4y_S_Gy3z_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Gy3z_S_vrr = PAZ*I_ERI_Gx3y_S_Gy3z_S_vrr+WPZ*I_ERI_Gx3y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Gx3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Gy3z_S_vrr = PAX*I_ERI_G2y2z_S_Gy3z_S_vrr+WPX*I_ERI_G2y2z_S_Gy3z_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Gy3z_S_vrr = PAY*I_ERI_Gx3z_S_Gy3z_S_vrr+WPY*I_ERI_Gx3z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Gx3z_S_F3z_S_M1_vrr;
      Double I_ERI_Hx4z_S_Gy3z_S_vrr = PAX*I_ERI_G4z_S_Gy3z_S_vrr+WPX*I_ERI_G4z_S_Gy3z_S_M1_vrr;
      Double I_ERI_H5y_S_Gy3z_S_vrr = PAY*I_ERI_G4y_S_Gy3z_S_vrr+WPY*I_ERI_G4y_S_Gy3z_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Gy3z_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Gy3z_S_M1_vrr+oned2k*I_ERI_G4y_S_F3z_S_M1_vrr;
      Double I_ERI_H4yz_S_Gy3z_S_vrr = PAZ*I_ERI_G4y_S_Gy3z_S_vrr+WPZ*I_ERI_G4y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_G4y_S_Fy2z_S_M1_vrr;
      Double I_ERI_H3y2z_S_Gy3z_S_vrr = PAZ*I_ERI_G3yz_S_Gy3z_S_vrr+WPZ*I_ERI_G3yz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_F3y_S_Gy3z_S_vrr-rhod2zsq*I_ERI_F3y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_G3yz_S_Fy2z_S_M1_vrr;
      Double I_ERI_H2y3z_S_Gy3z_S_vrr = PAY*I_ERI_Gy3z_S_Gy3z_S_vrr+WPY*I_ERI_Gy3z_S_Gy3z_S_M1_vrr+oned2z*I_ERI_F3z_S_Gy3z_S_vrr-rhod2zsq*I_ERI_F3z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Gy3z_S_F3z_S_M1_vrr;
      Double I_ERI_Hy4z_S_Gy3z_S_vrr = PAY*I_ERI_G4z_S_Gy3z_S_vrr+WPY*I_ERI_G4z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_G4z_S_F3z_S_M1_vrr;
      Double I_ERI_H5z_S_Gy3z_S_vrr = PAZ*I_ERI_G4z_S_Gy3z_S_vrr+WPZ*I_ERI_G4z_S_Gy3z_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Gy3z_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_G4z_S_Fy2z_S_M1_vrr;
      Double I_ERI_H5x_S_G4z_S_vrr = PAX*I_ERI_G4x_S_G4z_S_vrr+WPX*I_ERI_G4x_S_G4z_S_M1_vrr+4*oned2z*I_ERI_F3x_S_G4z_S_vrr-4*rhod2zsq*I_ERI_F3x_S_G4z_S_M1_vrr;
      Double I_ERI_H4xy_S_G4z_S_vrr = PAY*I_ERI_G4x_S_G4z_S_vrr+WPY*I_ERI_G4x_S_G4z_S_M1_vrr;
      Double I_ERI_H4xz_S_G4z_S_vrr = PAZ*I_ERI_G4x_S_G4z_S_vrr+WPZ*I_ERI_G4x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_G4x_S_F3z_S_M1_vrr;
      Double I_ERI_H3x2y_S_G4z_S_vrr = PAY*I_ERI_G3xy_S_G4z_S_vrr+WPY*I_ERI_G3xy_S_G4z_S_M1_vrr+oned2z*I_ERI_F3x_S_G4z_S_vrr-rhod2zsq*I_ERI_F3x_S_G4z_S_M1_vrr;
      Double I_ERI_H3xyz_S_G4z_S_vrr = PAZ*I_ERI_G3xy_S_G4z_S_vrr+WPZ*I_ERI_G3xy_S_G4z_S_M1_vrr+4*oned2k*I_ERI_G3xy_S_F3z_S_M1_vrr;
      Double I_ERI_H3x2z_S_G4z_S_vrr = PAZ*I_ERI_G3xz_S_G4z_S_vrr+WPZ*I_ERI_G3xz_S_G4z_S_M1_vrr+oned2z*I_ERI_F3x_S_G4z_S_vrr-rhod2zsq*I_ERI_F3x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_G3xz_S_F3z_S_M1_vrr;
      Double I_ERI_H2x3y_S_G4z_S_vrr = PAX*I_ERI_Gx3y_S_G4z_S_vrr+WPX*I_ERI_Gx3y_S_G4z_S_M1_vrr+oned2z*I_ERI_F3y_S_G4z_S_vrr-rhod2zsq*I_ERI_F3y_S_G4z_S_M1_vrr;
      Double I_ERI_H2x2yz_S_G4z_S_vrr = PAZ*I_ERI_G2x2y_S_G4z_S_vrr+WPZ*I_ERI_G2x2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_G2x2y_S_F3z_S_M1_vrr;
      Double I_ERI_H2xy2z_S_G4z_S_vrr = PAY*I_ERI_G2x2z_S_G4z_S_vrr+WPY*I_ERI_G2x2z_S_G4z_S_M1_vrr;
      Double I_ERI_H2x3z_S_G4z_S_vrr = PAX*I_ERI_Gx3z_S_G4z_S_vrr+WPX*I_ERI_Gx3z_S_G4z_S_M1_vrr+oned2z*I_ERI_F3z_S_G4z_S_vrr-rhod2zsq*I_ERI_F3z_S_G4z_S_M1_vrr;
      Double I_ERI_Hx4y_S_G4z_S_vrr = PAX*I_ERI_G4y_S_G4z_S_vrr+WPX*I_ERI_G4y_S_G4z_S_M1_vrr;
      Double I_ERI_Hx3yz_S_G4z_S_vrr = PAZ*I_ERI_Gx3y_S_G4z_S_vrr+WPZ*I_ERI_Gx3y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Gx3y_S_F3z_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_G4z_S_vrr = PAX*I_ERI_G2y2z_S_G4z_S_vrr+WPX*I_ERI_G2y2z_S_G4z_S_M1_vrr;
      Double I_ERI_Hxy3z_S_G4z_S_vrr = PAY*I_ERI_Gx3z_S_G4z_S_vrr+WPY*I_ERI_Gx3z_S_G4z_S_M1_vrr;
      Double I_ERI_Hx4z_S_G4z_S_vrr = PAX*I_ERI_G4z_S_G4z_S_vrr+WPX*I_ERI_G4z_S_G4z_S_M1_vrr;
      Double I_ERI_H5y_S_G4z_S_vrr = PAY*I_ERI_G4y_S_G4z_S_vrr+WPY*I_ERI_G4y_S_G4z_S_M1_vrr+4*oned2z*I_ERI_F3y_S_G4z_S_vrr-4*rhod2zsq*I_ERI_F3y_S_G4z_S_M1_vrr;
      Double I_ERI_H4yz_S_G4z_S_vrr = PAZ*I_ERI_G4y_S_G4z_S_vrr+WPZ*I_ERI_G4y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_G4y_S_F3z_S_M1_vrr;
      Double I_ERI_H3y2z_S_G4z_S_vrr = PAZ*I_ERI_G3yz_S_G4z_S_vrr+WPZ*I_ERI_G3yz_S_G4z_S_M1_vrr+oned2z*I_ERI_F3y_S_G4z_S_vrr-rhod2zsq*I_ERI_F3y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_G3yz_S_F3z_S_M1_vrr;
      Double I_ERI_H2y3z_S_G4z_S_vrr = PAY*I_ERI_Gy3z_S_G4z_S_vrr+WPY*I_ERI_Gy3z_S_G4z_S_M1_vrr+oned2z*I_ERI_F3z_S_G4z_S_vrr-rhod2zsq*I_ERI_F3z_S_G4z_S_M1_vrr;
      Double I_ERI_Hy4z_S_G4z_S_vrr = PAY*I_ERI_G4z_S_G4z_S_vrr+WPY*I_ERI_G4z_S_G4z_S_M1_vrr;
      Double I_ERI_H5z_S_G4z_S_vrr = PAZ*I_ERI_G4z_S_G4z_S_vrr+WPZ*I_ERI_G4z_S_G4z_S_M1_vrr+4*oned2z*I_ERI_F3z_S_G4z_S_vrr-4*rhod2zsq*I_ERI_F3z_S_G4z_S_M1_vrr+4*oned2k*I_ERI_G4z_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_G_S
       * RHS shell quartet name: SQ_ERI_H_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_G_S
       * RHS shell quartet name: SQ_ERI_G_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_F_S_M1
       ************************************************************/
      Double I_ERI_I6x_S_G4x_S_vrr = PAX*I_ERI_H5x_S_G4x_S_vrr+WPX*I_ERI_H5x_S_G4x_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G4x_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G4x_S_M1_vrr+4*oned2k*I_ERI_H5x_S_F3x_S_M1_vrr;
      Double I_ERI_I5xy_S_G4x_S_vrr = PAY*I_ERI_H5x_S_G4x_S_vrr+WPY*I_ERI_H5x_S_G4x_S_M1_vrr;
      Double I_ERI_I5xz_S_G4x_S_vrr = PAZ*I_ERI_H5x_S_G4x_S_vrr+WPZ*I_ERI_H5x_S_G4x_S_M1_vrr;
      Double I_ERI_I4x2y_S_G4x_S_vrr = PAY*I_ERI_H4xy_S_G4x_S_vrr+WPY*I_ERI_H4xy_S_G4x_S_M1_vrr+oned2z*I_ERI_G4x_S_G4x_S_vrr-rhod2zsq*I_ERI_G4x_S_G4x_S_M1_vrr;
      Double I_ERI_I4xyz_S_G4x_S_vrr = PAZ*I_ERI_H4xy_S_G4x_S_vrr+WPZ*I_ERI_H4xy_S_G4x_S_M1_vrr;
      Double I_ERI_I4x2z_S_G4x_S_vrr = PAZ*I_ERI_H4xz_S_G4x_S_vrr+WPZ*I_ERI_H4xz_S_G4x_S_M1_vrr+oned2z*I_ERI_G4x_S_G4x_S_vrr-rhod2zsq*I_ERI_G4x_S_G4x_S_M1_vrr;
      Double I_ERI_I3x3y_S_G4x_S_vrr = PAY*I_ERI_H3x2y_S_G4x_S_vrr+WPY*I_ERI_H3x2y_S_G4x_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G4x_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G4x_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G4x_S_vrr = PAZ*I_ERI_H3x2y_S_G4x_S_vrr+WPZ*I_ERI_H3x2y_S_G4x_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G4x_S_vrr = PAY*I_ERI_H3x2z_S_G4x_S_vrr+WPY*I_ERI_H3x2z_S_G4x_S_M1_vrr;
      Double I_ERI_I3x3z_S_G4x_S_vrr = PAZ*I_ERI_H3x2z_S_G4x_S_vrr+WPZ*I_ERI_H3x2z_S_G4x_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G4x_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G4x_S_M1_vrr;
      Double I_ERI_I2x4y_S_G4x_S_vrr = PAX*I_ERI_Hx4y_S_G4x_S_vrr+WPX*I_ERI_Hx4y_S_G4x_S_M1_vrr+oned2z*I_ERI_G4y_S_G4x_S_vrr-rhod2zsq*I_ERI_G4y_S_G4x_S_M1_vrr+4*oned2k*I_ERI_Hx4y_S_F3x_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G4x_S_vrr = PAZ*I_ERI_H2x3y_S_G4x_S_vrr+WPZ*I_ERI_H2x3y_S_G4x_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G4x_S_vrr = PAZ*I_ERI_H2x2yz_S_G4x_S_vrr+WPZ*I_ERI_H2x2yz_S_G4x_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G4x_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G4x_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G4x_S_vrr = PAY*I_ERI_H2x3z_S_G4x_S_vrr+WPY*I_ERI_H2x3z_S_G4x_S_M1_vrr;
      Double I_ERI_I2x4z_S_G4x_S_vrr = PAX*I_ERI_Hx4z_S_G4x_S_vrr+WPX*I_ERI_Hx4z_S_G4x_S_M1_vrr+oned2z*I_ERI_G4z_S_G4x_S_vrr-rhod2zsq*I_ERI_G4z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_Hx4z_S_F3x_S_M1_vrr;
      Double I_ERI_Ix5y_S_G4x_S_vrr = PAX*I_ERI_H5y_S_G4x_S_vrr+WPX*I_ERI_H5y_S_G4x_S_M1_vrr+4*oned2k*I_ERI_H5y_S_F3x_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G4x_S_vrr = PAZ*I_ERI_Hx4y_S_G4x_S_vrr+WPZ*I_ERI_Hx4y_S_G4x_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G4x_S_vrr = PAX*I_ERI_H3y2z_S_G4x_S_vrr+WPX*I_ERI_H3y2z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_H3y2z_S_F3x_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G4x_S_vrr = PAX*I_ERI_H2y3z_S_G4x_S_vrr+WPX*I_ERI_H2y3z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_H2y3z_S_F3x_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G4x_S_vrr = PAY*I_ERI_Hx4z_S_G4x_S_vrr+WPY*I_ERI_Hx4z_S_G4x_S_M1_vrr;
      Double I_ERI_Ix5z_S_G4x_S_vrr = PAX*I_ERI_H5z_S_G4x_S_vrr+WPX*I_ERI_H5z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_H5z_S_F3x_S_M1_vrr;
      Double I_ERI_I6y_S_G4x_S_vrr = PAY*I_ERI_H5y_S_G4x_S_vrr+WPY*I_ERI_H5y_S_G4x_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G4x_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G4x_S_M1_vrr;
      Double I_ERI_I5yz_S_G4x_S_vrr = PAZ*I_ERI_H5y_S_G4x_S_vrr+WPZ*I_ERI_H5y_S_G4x_S_M1_vrr;
      Double I_ERI_I4y2z_S_G4x_S_vrr = PAZ*I_ERI_H4yz_S_G4x_S_vrr+WPZ*I_ERI_H4yz_S_G4x_S_M1_vrr+oned2z*I_ERI_G4y_S_G4x_S_vrr-rhod2zsq*I_ERI_G4y_S_G4x_S_M1_vrr;
      Double I_ERI_I3y3z_S_G4x_S_vrr = PAZ*I_ERI_H3y2z_S_G4x_S_vrr+WPZ*I_ERI_H3y2z_S_G4x_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G4x_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G4x_S_M1_vrr;
      Double I_ERI_I2y4z_S_G4x_S_vrr = PAY*I_ERI_Hy4z_S_G4x_S_vrr+WPY*I_ERI_Hy4z_S_G4x_S_M1_vrr+oned2z*I_ERI_G4z_S_G4x_S_vrr-rhod2zsq*I_ERI_G4z_S_G4x_S_M1_vrr;
      Double I_ERI_Iy5z_S_G4x_S_vrr = PAY*I_ERI_H5z_S_G4x_S_vrr+WPY*I_ERI_H5z_S_G4x_S_M1_vrr;
      Double I_ERI_I6z_S_G4x_S_vrr = PAZ*I_ERI_H5z_S_G4x_S_vrr+WPZ*I_ERI_H5z_S_G4x_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G4x_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G4x_S_M1_vrr;
      Double I_ERI_I6x_S_G3xy_S_vrr = PAX*I_ERI_H5x_S_G3xy_S_vrr+WPX*I_ERI_H5x_S_G3xy_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G3xy_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_H5x_S_F2xy_S_M1_vrr;
      Double I_ERI_I5xy_S_G3xy_S_vrr = PAY*I_ERI_H5x_S_G3xy_S_vrr+WPY*I_ERI_H5x_S_G3xy_S_M1_vrr+oned2k*I_ERI_H5x_S_F3x_S_M1_vrr;
      Double I_ERI_I5xz_S_G3xy_S_vrr = PAZ*I_ERI_H5x_S_G3xy_S_vrr+WPZ*I_ERI_H5x_S_G3xy_S_M1_vrr;
      Double I_ERI_I4x2y_S_G3xy_S_vrr = PAY*I_ERI_H4xy_S_G3xy_S_vrr+WPY*I_ERI_H4xy_S_G3xy_S_M1_vrr+oned2z*I_ERI_G4x_S_G3xy_S_vrr-rhod2zsq*I_ERI_G4x_S_G3xy_S_M1_vrr+oned2k*I_ERI_H4xy_S_F3x_S_M1_vrr;
      Double I_ERI_I4xyz_S_G3xy_S_vrr = PAZ*I_ERI_H4xy_S_G3xy_S_vrr+WPZ*I_ERI_H4xy_S_G3xy_S_M1_vrr;
      Double I_ERI_I4x2z_S_G3xy_S_vrr = PAZ*I_ERI_H4xz_S_G3xy_S_vrr+WPZ*I_ERI_H4xz_S_G3xy_S_M1_vrr+oned2z*I_ERI_G4x_S_G3xy_S_vrr-rhod2zsq*I_ERI_G4x_S_G3xy_S_M1_vrr;
      Double I_ERI_I3x3y_S_G3xy_S_vrr = PAY*I_ERI_H3x2y_S_G3xy_S_vrr+WPY*I_ERI_H3x2y_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G3xy_S_M1_vrr+oned2k*I_ERI_H3x2y_S_F3x_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G3xy_S_vrr = PAZ*I_ERI_H3x2y_S_G3xy_S_vrr+WPZ*I_ERI_H3x2y_S_G3xy_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G3xy_S_vrr = PAY*I_ERI_H3x2z_S_G3xy_S_vrr+WPY*I_ERI_H3x2z_S_G3xy_S_M1_vrr+oned2k*I_ERI_H3x2z_S_F3x_S_M1_vrr;
      Double I_ERI_I3x3z_S_G3xy_S_vrr = PAZ*I_ERI_H3x2z_S_G3xy_S_vrr+WPZ*I_ERI_H3x2z_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G3xy_S_M1_vrr;
      Double I_ERI_I2x4y_S_G3xy_S_vrr = PAX*I_ERI_Hx4y_S_G3xy_S_vrr+WPX*I_ERI_Hx4y_S_G3xy_S_M1_vrr+oned2z*I_ERI_G4y_S_G3xy_S_vrr-rhod2zsq*I_ERI_G4y_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_Hx4y_S_F2xy_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G3xy_S_vrr = PAZ*I_ERI_H2x3y_S_G3xy_S_vrr+WPZ*I_ERI_H2x3y_S_G3xy_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G3xy_S_vrr = PAZ*I_ERI_H2x2yz_S_G3xy_S_vrr+WPZ*I_ERI_H2x2yz_S_G3xy_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G3xy_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G3xy_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G3xy_S_vrr = PAY*I_ERI_H2x3z_S_G3xy_S_vrr+WPY*I_ERI_H2x3z_S_G3xy_S_M1_vrr+oned2k*I_ERI_H2x3z_S_F3x_S_M1_vrr;
      Double I_ERI_I2x4z_S_G3xy_S_vrr = PAX*I_ERI_Hx4z_S_G3xy_S_vrr+WPX*I_ERI_Hx4z_S_G3xy_S_M1_vrr+oned2z*I_ERI_G4z_S_G3xy_S_vrr-rhod2zsq*I_ERI_G4z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_Hx4z_S_F2xy_S_M1_vrr;
      Double I_ERI_Ix5y_S_G3xy_S_vrr = PAX*I_ERI_H5y_S_G3xy_S_vrr+WPX*I_ERI_H5y_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_H5y_S_F2xy_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G3xy_S_vrr = PAZ*I_ERI_Hx4y_S_G3xy_S_vrr+WPZ*I_ERI_Hx4y_S_G3xy_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G3xy_S_vrr = PAX*I_ERI_H3y2z_S_G3xy_S_vrr+WPX*I_ERI_H3y2z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_H3y2z_S_F2xy_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G3xy_S_vrr = PAX*I_ERI_H2y3z_S_G3xy_S_vrr+WPX*I_ERI_H2y3z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_H2y3z_S_F2xy_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G3xy_S_vrr = PAY*I_ERI_Hx4z_S_G3xy_S_vrr+WPY*I_ERI_Hx4z_S_G3xy_S_M1_vrr+oned2k*I_ERI_Hx4z_S_F3x_S_M1_vrr;
      Double I_ERI_Ix5z_S_G3xy_S_vrr = PAX*I_ERI_H5z_S_G3xy_S_vrr+WPX*I_ERI_H5z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_H5z_S_F2xy_S_M1_vrr;
      Double I_ERI_I6y_S_G3xy_S_vrr = PAY*I_ERI_H5y_S_G3xy_S_vrr+WPY*I_ERI_H5y_S_G3xy_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G3xy_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G3xy_S_M1_vrr+oned2k*I_ERI_H5y_S_F3x_S_M1_vrr;
      Double I_ERI_I5yz_S_G3xy_S_vrr = PAZ*I_ERI_H5y_S_G3xy_S_vrr+WPZ*I_ERI_H5y_S_G3xy_S_M1_vrr;
      Double I_ERI_I4y2z_S_G3xy_S_vrr = PAZ*I_ERI_H4yz_S_G3xy_S_vrr+WPZ*I_ERI_H4yz_S_G3xy_S_M1_vrr+oned2z*I_ERI_G4y_S_G3xy_S_vrr-rhod2zsq*I_ERI_G4y_S_G3xy_S_M1_vrr;
      Double I_ERI_I3y3z_S_G3xy_S_vrr = PAZ*I_ERI_H3y2z_S_G3xy_S_vrr+WPZ*I_ERI_H3y2z_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G3xy_S_M1_vrr;
      Double I_ERI_I2y4z_S_G3xy_S_vrr = PAY*I_ERI_Hy4z_S_G3xy_S_vrr+WPY*I_ERI_Hy4z_S_G3xy_S_M1_vrr+oned2z*I_ERI_G4z_S_G3xy_S_vrr-rhod2zsq*I_ERI_G4z_S_G3xy_S_M1_vrr+oned2k*I_ERI_Hy4z_S_F3x_S_M1_vrr;
      Double I_ERI_Iy5z_S_G3xy_S_vrr = PAY*I_ERI_H5z_S_G3xy_S_vrr+WPY*I_ERI_H5z_S_G3xy_S_M1_vrr+oned2k*I_ERI_H5z_S_F3x_S_M1_vrr;
      Double I_ERI_I6z_S_G3xy_S_vrr = PAZ*I_ERI_H5z_S_G3xy_S_vrr+WPZ*I_ERI_H5z_S_G3xy_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G3xy_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G3xy_S_M1_vrr;
      Double I_ERI_I6x_S_G3xz_S_vrr = PAX*I_ERI_H5x_S_G3xz_S_vrr+WPX*I_ERI_H5x_S_G3xz_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G3xz_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_H5x_S_F2xz_S_M1_vrr;
      Double I_ERI_I5xy_S_G3xz_S_vrr = PAY*I_ERI_H5x_S_G3xz_S_vrr+WPY*I_ERI_H5x_S_G3xz_S_M1_vrr;
      Double I_ERI_I5xz_S_G3xz_S_vrr = PAZ*I_ERI_H5x_S_G3xz_S_vrr+WPZ*I_ERI_H5x_S_G3xz_S_M1_vrr+oned2k*I_ERI_H5x_S_F3x_S_M1_vrr;
      Double I_ERI_I4x2y_S_G3xz_S_vrr = PAY*I_ERI_H4xy_S_G3xz_S_vrr+WPY*I_ERI_H4xy_S_G3xz_S_M1_vrr+oned2z*I_ERI_G4x_S_G3xz_S_vrr-rhod2zsq*I_ERI_G4x_S_G3xz_S_M1_vrr;
      Double I_ERI_I4xyz_S_G3xz_S_vrr = PAZ*I_ERI_H4xy_S_G3xz_S_vrr+WPZ*I_ERI_H4xy_S_G3xz_S_M1_vrr+oned2k*I_ERI_H4xy_S_F3x_S_M1_vrr;
      Double I_ERI_I4x2z_S_G3xz_S_vrr = PAZ*I_ERI_H4xz_S_G3xz_S_vrr+WPZ*I_ERI_H4xz_S_G3xz_S_M1_vrr+oned2z*I_ERI_G4x_S_G3xz_S_vrr-rhod2zsq*I_ERI_G4x_S_G3xz_S_M1_vrr+oned2k*I_ERI_H4xz_S_F3x_S_M1_vrr;
      Double I_ERI_I3x3y_S_G3xz_S_vrr = PAY*I_ERI_H3x2y_S_G3xz_S_vrr+WPY*I_ERI_H3x2y_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G3xz_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G3xz_S_vrr = PAZ*I_ERI_H3x2y_S_G3xz_S_vrr+WPZ*I_ERI_H3x2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_H3x2y_S_F3x_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G3xz_S_vrr = PAY*I_ERI_H3x2z_S_G3xz_S_vrr+WPY*I_ERI_H3x2z_S_G3xz_S_M1_vrr;
      Double I_ERI_I3x3z_S_G3xz_S_vrr = PAZ*I_ERI_H3x2z_S_G3xz_S_vrr+WPZ*I_ERI_H3x2z_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G3xz_S_M1_vrr+oned2k*I_ERI_H3x2z_S_F3x_S_M1_vrr;
      Double I_ERI_I2x4y_S_G3xz_S_vrr = PAX*I_ERI_Hx4y_S_G3xz_S_vrr+WPX*I_ERI_Hx4y_S_G3xz_S_M1_vrr+oned2z*I_ERI_G4y_S_G3xz_S_vrr-rhod2zsq*I_ERI_G4y_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_Hx4y_S_F2xz_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G3xz_S_vrr = PAZ*I_ERI_H2x3y_S_G3xz_S_vrr+WPZ*I_ERI_H2x3y_S_G3xz_S_M1_vrr+oned2k*I_ERI_H2x3y_S_F3x_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G3xz_S_vrr = PAZ*I_ERI_H2x2yz_S_G3xz_S_vrr+WPZ*I_ERI_H2x2yz_S_G3xz_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G3xz_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_H2x2yz_S_F3x_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G3xz_S_vrr = PAY*I_ERI_H2x3z_S_G3xz_S_vrr+WPY*I_ERI_H2x3z_S_G3xz_S_M1_vrr;
      Double I_ERI_I2x4z_S_G3xz_S_vrr = PAX*I_ERI_Hx4z_S_G3xz_S_vrr+WPX*I_ERI_Hx4z_S_G3xz_S_M1_vrr+oned2z*I_ERI_G4z_S_G3xz_S_vrr-rhod2zsq*I_ERI_G4z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_Hx4z_S_F2xz_S_M1_vrr;
      Double I_ERI_Ix5y_S_G3xz_S_vrr = PAX*I_ERI_H5y_S_G3xz_S_vrr+WPX*I_ERI_H5y_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_H5y_S_F2xz_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G3xz_S_vrr = PAZ*I_ERI_Hx4y_S_G3xz_S_vrr+WPZ*I_ERI_Hx4y_S_G3xz_S_M1_vrr+oned2k*I_ERI_Hx4y_S_F3x_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G3xz_S_vrr = PAX*I_ERI_H3y2z_S_G3xz_S_vrr+WPX*I_ERI_H3y2z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_H3y2z_S_F2xz_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G3xz_S_vrr = PAX*I_ERI_H2y3z_S_G3xz_S_vrr+WPX*I_ERI_H2y3z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_H2y3z_S_F2xz_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G3xz_S_vrr = PAY*I_ERI_Hx4z_S_G3xz_S_vrr+WPY*I_ERI_Hx4z_S_G3xz_S_M1_vrr;
      Double I_ERI_Ix5z_S_G3xz_S_vrr = PAX*I_ERI_H5z_S_G3xz_S_vrr+WPX*I_ERI_H5z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_H5z_S_F2xz_S_M1_vrr;
      Double I_ERI_I6y_S_G3xz_S_vrr = PAY*I_ERI_H5y_S_G3xz_S_vrr+WPY*I_ERI_H5y_S_G3xz_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G3xz_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G3xz_S_M1_vrr;
      Double I_ERI_I5yz_S_G3xz_S_vrr = PAZ*I_ERI_H5y_S_G3xz_S_vrr+WPZ*I_ERI_H5y_S_G3xz_S_M1_vrr+oned2k*I_ERI_H5y_S_F3x_S_M1_vrr;
      Double I_ERI_I4y2z_S_G3xz_S_vrr = PAZ*I_ERI_H4yz_S_G3xz_S_vrr+WPZ*I_ERI_H4yz_S_G3xz_S_M1_vrr+oned2z*I_ERI_G4y_S_G3xz_S_vrr-rhod2zsq*I_ERI_G4y_S_G3xz_S_M1_vrr+oned2k*I_ERI_H4yz_S_F3x_S_M1_vrr;
      Double I_ERI_I3y3z_S_G3xz_S_vrr = PAZ*I_ERI_H3y2z_S_G3xz_S_vrr+WPZ*I_ERI_H3y2z_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G3xz_S_M1_vrr+oned2k*I_ERI_H3y2z_S_F3x_S_M1_vrr;
      Double I_ERI_I2y4z_S_G3xz_S_vrr = PAY*I_ERI_Hy4z_S_G3xz_S_vrr+WPY*I_ERI_Hy4z_S_G3xz_S_M1_vrr+oned2z*I_ERI_G4z_S_G3xz_S_vrr-rhod2zsq*I_ERI_G4z_S_G3xz_S_M1_vrr;
      Double I_ERI_Iy5z_S_G3xz_S_vrr = PAY*I_ERI_H5z_S_G3xz_S_vrr+WPY*I_ERI_H5z_S_G3xz_S_M1_vrr;
      Double I_ERI_I6z_S_G3xz_S_vrr = PAZ*I_ERI_H5z_S_G3xz_S_vrr+WPZ*I_ERI_H5z_S_G3xz_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G3xz_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G3xz_S_M1_vrr+oned2k*I_ERI_H5z_S_F3x_S_M1_vrr;
      Double I_ERI_I6x_S_G2x2y_S_vrr = PAX*I_ERI_H5x_S_G2x2y_S_vrr+WPX*I_ERI_H5x_S_G2x2y_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G2x2y_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H5x_S_Fx2y_S_M1_vrr;
      Double I_ERI_I5xy_S_G2x2y_S_vrr = PAY*I_ERI_H5x_S_G2x2y_S_vrr+WPY*I_ERI_H5x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H5x_S_F2xy_S_M1_vrr;
      Double I_ERI_I5xz_S_G2x2y_S_vrr = PAZ*I_ERI_H5x_S_G2x2y_S_vrr+WPZ*I_ERI_H5x_S_G2x2y_S_M1_vrr;
      Double I_ERI_I4x2y_S_G2x2y_S_vrr = PAY*I_ERI_H4xy_S_G2x2y_S_vrr+WPY*I_ERI_H4xy_S_G2x2y_S_M1_vrr+oned2z*I_ERI_G4x_S_G2x2y_S_vrr-rhod2zsq*I_ERI_G4x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H4xy_S_F2xy_S_M1_vrr;
      Double I_ERI_I4xyz_S_G2x2y_S_vrr = PAZ*I_ERI_H4xy_S_G2x2y_S_vrr+WPZ*I_ERI_H4xy_S_G2x2y_S_M1_vrr;
      Double I_ERI_I4x2z_S_G2x2y_S_vrr = PAZ*I_ERI_H4xz_S_G2x2y_S_vrr+WPZ*I_ERI_H4xz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_G4x_S_G2x2y_S_vrr-rhod2zsq*I_ERI_G4x_S_G2x2y_S_M1_vrr;
      Double I_ERI_I3x3y_S_G2x2y_S_vrr = PAY*I_ERI_H3x2y_S_G2x2y_S_vrr+WPY*I_ERI_H3x2y_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H3x2y_S_F2xy_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G2x2y_S_vrr = PAZ*I_ERI_H3x2y_S_G2x2y_S_vrr+WPZ*I_ERI_H3x2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G2x2y_S_vrr = PAY*I_ERI_H3x2z_S_G2x2y_S_vrr+WPY*I_ERI_H3x2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H3x2z_S_F2xy_S_M1_vrr;
      Double I_ERI_I3x3z_S_G2x2y_S_vrr = PAZ*I_ERI_H3x2z_S_G2x2y_S_vrr+WPZ*I_ERI_H3x2z_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G2x2y_S_M1_vrr;
      Double I_ERI_I2x4y_S_G2x2y_S_vrr = PAX*I_ERI_Hx4y_S_G2x2y_S_vrr+WPX*I_ERI_Hx4y_S_G2x2y_S_M1_vrr+oned2z*I_ERI_G4y_S_G2x2y_S_vrr-rhod2zsq*I_ERI_G4y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Hx4y_S_Fx2y_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G2x2y_S_vrr = PAZ*I_ERI_H2x3y_S_G2x2y_S_vrr+WPZ*I_ERI_H2x3y_S_G2x2y_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G2x2y_S_vrr = PAZ*I_ERI_H2x2yz_S_G2x2y_S_vrr+WPZ*I_ERI_H2x2yz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G2x2y_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G2x2y_S_vrr = PAY*I_ERI_H2x3z_S_G2x2y_S_vrr+WPY*I_ERI_H2x3z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H2x3z_S_F2xy_S_M1_vrr;
      Double I_ERI_I2x4z_S_G2x2y_S_vrr = PAX*I_ERI_Hx4z_S_G2x2y_S_vrr+WPX*I_ERI_Hx4z_S_G2x2y_S_M1_vrr+oned2z*I_ERI_G4z_S_G2x2y_S_vrr-rhod2zsq*I_ERI_G4z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Hx4z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Ix5y_S_G2x2y_S_vrr = PAX*I_ERI_H5y_S_G2x2y_S_vrr+WPX*I_ERI_H5y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H5y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G2x2y_S_vrr = PAZ*I_ERI_Hx4y_S_G2x2y_S_vrr+WPZ*I_ERI_Hx4y_S_G2x2y_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G2x2y_S_vrr = PAX*I_ERI_H3y2z_S_G2x2y_S_vrr+WPX*I_ERI_H3y2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H3y2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G2x2y_S_vrr = PAX*I_ERI_H2y3z_S_G2x2y_S_vrr+WPX*I_ERI_H2y3z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H2y3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G2x2y_S_vrr = PAY*I_ERI_Hx4z_S_G2x2y_S_vrr+WPY*I_ERI_Hx4z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Hx4z_S_F2xy_S_M1_vrr;
      Double I_ERI_Ix5z_S_G2x2y_S_vrr = PAX*I_ERI_H5z_S_G2x2y_S_vrr+WPX*I_ERI_H5z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H5z_S_Fx2y_S_M1_vrr;
      Double I_ERI_I6y_S_G2x2y_S_vrr = PAY*I_ERI_H5y_S_G2x2y_S_vrr+WPY*I_ERI_H5y_S_G2x2y_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G2x2y_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H5y_S_F2xy_S_M1_vrr;
      Double I_ERI_I5yz_S_G2x2y_S_vrr = PAZ*I_ERI_H5y_S_G2x2y_S_vrr+WPZ*I_ERI_H5y_S_G2x2y_S_M1_vrr;
      Double I_ERI_I4y2z_S_G2x2y_S_vrr = PAZ*I_ERI_H4yz_S_G2x2y_S_vrr+WPZ*I_ERI_H4yz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_G4y_S_G2x2y_S_vrr-rhod2zsq*I_ERI_G4y_S_G2x2y_S_M1_vrr;
      Double I_ERI_I3y3z_S_G2x2y_S_vrr = PAZ*I_ERI_H3y2z_S_G2x2y_S_vrr+WPZ*I_ERI_H3y2z_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G2x2y_S_M1_vrr;
      Double I_ERI_I2y4z_S_G2x2y_S_vrr = PAY*I_ERI_Hy4z_S_G2x2y_S_vrr+WPY*I_ERI_Hy4z_S_G2x2y_S_M1_vrr+oned2z*I_ERI_G4z_S_G2x2y_S_vrr-rhod2zsq*I_ERI_G4z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Hy4z_S_F2xy_S_M1_vrr;
      Double I_ERI_Iy5z_S_G2x2y_S_vrr = PAY*I_ERI_H5z_S_G2x2y_S_vrr+WPY*I_ERI_H5z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_H5z_S_F2xy_S_M1_vrr;
      Double I_ERI_I6z_S_G2x2y_S_vrr = PAZ*I_ERI_H5z_S_G2x2y_S_vrr+WPZ*I_ERI_H5z_S_G2x2y_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G2x2y_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G2x2y_S_M1_vrr;
      Double I_ERI_I6x_S_G2xyz_S_vrr = PAX*I_ERI_H5x_S_G2xyz_S_vrr+WPX*I_ERI_H5x_S_G2xyz_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G2xyz_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_H5x_S_Fxyz_S_M1_vrr;
      Double I_ERI_I5xy_S_G2xyz_S_vrr = PAY*I_ERI_H5x_S_G2xyz_S_vrr+WPY*I_ERI_H5x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H5x_S_F2xz_S_M1_vrr;
      Double I_ERI_I5xz_S_G2xyz_S_vrr = PAZ*I_ERI_H5x_S_G2xyz_S_vrr+WPZ*I_ERI_H5x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H5x_S_F2xy_S_M1_vrr;
      Double I_ERI_I4x2y_S_G2xyz_S_vrr = PAY*I_ERI_H4xy_S_G2xyz_S_vrr+WPY*I_ERI_H4xy_S_G2xyz_S_M1_vrr+oned2z*I_ERI_G4x_S_G2xyz_S_vrr-rhod2zsq*I_ERI_G4x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H4xy_S_F2xz_S_M1_vrr;
      Double I_ERI_I4xyz_S_G2xyz_S_vrr = PAZ*I_ERI_H4xy_S_G2xyz_S_vrr+WPZ*I_ERI_H4xy_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H4xy_S_F2xy_S_M1_vrr;
      Double I_ERI_I4x2z_S_G2xyz_S_vrr = PAZ*I_ERI_H4xz_S_G2xyz_S_vrr+WPZ*I_ERI_H4xz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_G4x_S_G2xyz_S_vrr-rhod2zsq*I_ERI_G4x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H4xz_S_F2xy_S_M1_vrr;
      Double I_ERI_I3x3y_S_G2xyz_S_vrr = PAY*I_ERI_H3x2y_S_G2xyz_S_vrr+WPY*I_ERI_H3x2y_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H3x2y_S_F2xz_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G2xyz_S_vrr = PAZ*I_ERI_H3x2y_S_G2xyz_S_vrr+WPZ*I_ERI_H3x2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H3x2y_S_F2xy_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G2xyz_S_vrr = PAY*I_ERI_H3x2z_S_G2xyz_S_vrr+WPY*I_ERI_H3x2z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H3x2z_S_F2xz_S_M1_vrr;
      Double I_ERI_I3x3z_S_G2xyz_S_vrr = PAZ*I_ERI_H3x2z_S_G2xyz_S_vrr+WPZ*I_ERI_H3x2z_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H3x2z_S_F2xy_S_M1_vrr;
      Double I_ERI_I2x4y_S_G2xyz_S_vrr = PAX*I_ERI_Hx4y_S_G2xyz_S_vrr+WPX*I_ERI_Hx4y_S_G2xyz_S_M1_vrr+oned2z*I_ERI_G4y_S_G2xyz_S_vrr-rhod2zsq*I_ERI_G4y_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_Hx4y_S_Fxyz_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G2xyz_S_vrr = PAZ*I_ERI_H2x3y_S_G2xyz_S_vrr+WPZ*I_ERI_H2x3y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H2x3y_S_F2xy_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G2xyz_S_vrr = PAZ*I_ERI_H2x2yz_S_G2xyz_S_vrr+WPZ*I_ERI_H2x2yz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G2xyz_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H2x2yz_S_F2xy_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G2xyz_S_vrr = PAY*I_ERI_H2x3z_S_G2xyz_S_vrr+WPY*I_ERI_H2x3z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H2x3z_S_F2xz_S_M1_vrr;
      Double I_ERI_I2x4z_S_G2xyz_S_vrr = PAX*I_ERI_Hx4z_S_G2xyz_S_vrr+WPX*I_ERI_Hx4z_S_G2xyz_S_M1_vrr+oned2z*I_ERI_G4z_S_G2xyz_S_vrr-rhod2zsq*I_ERI_G4z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_Hx4z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Ix5y_S_G2xyz_S_vrr = PAX*I_ERI_H5y_S_G2xyz_S_vrr+WPX*I_ERI_H5y_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_H5y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G2xyz_S_vrr = PAZ*I_ERI_Hx4y_S_G2xyz_S_vrr+WPZ*I_ERI_Hx4y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Hx4y_S_F2xy_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G2xyz_S_vrr = PAX*I_ERI_H3y2z_S_G2xyz_S_vrr+WPX*I_ERI_H3y2z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_H3y2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G2xyz_S_vrr = PAX*I_ERI_H2y3z_S_G2xyz_S_vrr+WPX*I_ERI_H2y3z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_H2y3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G2xyz_S_vrr = PAY*I_ERI_Hx4z_S_G2xyz_S_vrr+WPY*I_ERI_Hx4z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Hx4z_S_F2xz_S_M1_vrr;
      Double I_ERI_Ix5z_S_G2xyz_S_vrr = PAX*I_ERI_H5z_S_G2xyz_S_vrr+WPX*I_ERI_H5z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_H5z_S_Fxyz_S_M1_vrr;
      Double I_ERI_I6y_S_G2xyz_S_vrr = PAY*I_ERI_H5y_S_G2xyz_S_vrr+WPY*I_ERI_H5y_S_G2xyz_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G2xyz_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H5y_S_F2xz_S_M1_vrr;
      Double I_ERI_I5yz_S_G2xyz_S_vrr = PAZ*I_ERI_H5y_S_G2xyz_S_vrr+WPZ*I_ERI_H5y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H5y_S_F2xy_S_M1_vrr;
      Double I_ERI_I4y2z_S_G2xyz_S_vrr = PAZ*I_ERI_H4yz_S_G2xyz_S_vrr+WPZ*I_ERI_H4yz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_G4y_S_G2xyz_S_vrr-rhod2zsq*I_ERI_G4y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H4yz_S_F2xy_S_M1_vrr;
      Double I_ERI_I3y3z_S_G2xyz_S_vrr = PAZ*I_ERI_H3y2z_S_G2xyz_S_vrr+WPZ*I_ERI_H3y2z_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H3y2z_S_F2xy_S_M1_vrr;
      Double I_ERI_I2y4z_S_G2xyz_S_vrr = PAY*I_ERI_Hy4z_S_G2xyz_S_vrr+WPY*I_ERI_Hy4z_S_G2xyz_S_M1_vrr+oned2z*I_ERI_G4z_S_G2xyz_S_vrr-rhod2zsq*I_ERI_G4z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Hy4z_S_F2xz_S_M1_vrr;
      Double I_ERI_Iy5z_S_G2xyz_S_vrr = PAY*I_ERI_H5z_S_G2xyz_S_vrr+WPY*I_ERI_H5z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H5z_S_F2xz_S_M1_vrr;
      Double I_ERI_I6z_S_G2xyz_S_vrr = PAZ*I_ERI_H5z_S_G2xyz_S_vrr+WPZ*I_ERI_H5z_S_G2xyz_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G2xyz_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_H5z_S_F2xy_S_M1_vrr;
      Double I_ERI_I6x_S_G2x2z_S_vrr = PAX*I_ERI_H5x_S_G2x2z_S_vrr+WPX*I_ERI_H5x_S_G2x2z_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G2x2z_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H5x_S_Fx2z_S_M1_vrr;
      Double I_ERI_I5xy_S_G2x2z_S_vrr = PAY*I_ERI_H5x_S_G2x2z_S_vrr+WPY*I_ERI_H5x_S_G2x2z_S_M1_vrr;
      Double I_ERI_I5xz_S_G2x2z_S_vrr = PAZ*I_ERI_H5x_S_G2x2z_S_vrr+WPZ*I_ERI_H5x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H5x_S_F2xz_S_M1_vrr;
      Double I_ERI_I4x2y_S_G2x2z_S_vrr = PAY*I_ERI_H4xy_S_G2x2z_S_vrr+WPY*I_ERI_H4xy_S_G2x2z_S_M1_vrr+oned2z*I_ERI_G4x_S_G2x2z_S_vrr-rhod2zsq*I_ERI_G4x_S_G2x2z_S_M1_vrr;
      Double I_ERI_I4xyz_S_G2x2z_S_vrr = PAZ*I_ERI_H4xy_S_G2x2z_S_vrr+WPZ*I_ERI_H4xy_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H4xy_S_F2xz_S_M1_vrr;
      Double I_ERI_I4x2z_S_G2x2z_S_vrr = PAZ*I_ERI_H4xz_S_G2x2z_S_vrr+WPZ*I_ERI_H4xz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_G4x_S_G2x2z_S_vrr-rhod2zsq*I_ERI_G4x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H4xz_S_F2xz_S_M1_vrr;
      Double I_ERI_I3x3y_S_G2x2z_S_vrr = PAY*I_ERI_H3x2y_S_G2x2z_S_vrr+WPY*I_ERI_H3x2y_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G2x2z_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G2x2z_S_vrr = PAZ*I_ERI_H3x2y_S_G2x2z_S_vrr+WPZ*I_ERI_H3x2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H3x2y_S_F2xz_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G2x2z_S_vrr = PAY*I_ERI_H3x2z_S_G2x2z_S_vrr+WPY*I_ERI_H3x2z_S_G2x2z_S_M1_vrr;
      Double I_ERI_I3x3z_S_G2x2z_S_vrr = PAZ*I_ERI_H3x2z_S_G2x2z_S_vrr+WPZ*I_ERI_H3x2z_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H3x2z_S_F2xz_S_M1_vrr;
      Double I_ERI_I2x4y_S_G2x2z_S_vrr = PAX*I_ERI_Hx4y_S_G2x2z_S_vrr+WPX*I_ERI_Hx4y_S_G2x2z_S_M1_vrr+oned2z*I_ERI_G4y_S_G2x2z_S_vrr-rhod2zsq*I_ERI_G4y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Hx4y_S_Fx2z_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G2x2z_S_vrr = PAZ*I_ERI_H2x3y_S_G2x2z_S_vrr+WPZ*I_ERI_H2x3y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H2x3y_S_F2xz_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G2x2z_S_vrr = PAZ*I_ERI_H2x2yz_S_G2x2z_S_vrr+WPZ*I_ERI_H2x2yz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G2x2z_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H2x2yz_S_F2xz_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G2x2z_S_vrr = PAY*I_ERI_H2x3z_S_G2x2z_S_vrr+WPY*I_ERI_H2x3z_S_G2x2z_S_M1_vrr;
      Double I_ERI_I2x4z_S_G2x2z_S_vrr = PAX*I_ERI_Hx4z_S_G2x2z_S_vrr+WPX*I_ERI_Hx4z_S_G2x2z_S_M1_vrr+oned2z*I_ERI_G4z_S_G2x2z_S_vrr-rhod2zsq*I_ERI_G4z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Hx4z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Ix5y_S_G2x2z_S_vrr = PAX*I_ERI_H5y_S_G2x2z_S_vrr+WPX*I_ERI_H5y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H5y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G2x2z_S_vrr = PAZ*I_ERI_Hx4y_S_G2x2z_S_vrr+WPZ*I_ERI_Hx4y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Hx4y_S_F2xz_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G2x2z_S_vrr = PAX*I_ERI_H3y2z_S_G2x2z_S_vrr+WPX*I_ERI_H3y2z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H3y2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G2x2z_S_vrr = PAX*I_ERI_H2y3z_S_G2x2z_S_vrr+WPX*I_ERI_H2y3z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H2y3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G2x2z_S_vrr = PAY*I_ERI_Hx4z_S_G2x2z_S_vrr+WPY*I_ERI_Hx4z_S_G2x2z_S_M1_vrr;
      Double I_ERI_Ix5z_S_G2x2z_S_vrr = PAX*I_ERI_H5z_S_G2x2z_S_vrr+WPX*I_ERI_H5z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H5z_S_Fx2z_S_M1_vrr;
      Double I_ERI_I6y_S_G2x2z_S_vrr = PAY*I_ERI_H5y_S_G2x2z_S_vrr+WPY*I_ERI_H5y_S_G2x2z_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G2x2z_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G2x2z_S_M1_vrr;
      Double I_ERI_I5yz_S_G2x2z_S_vrr = PAZ*I_ERI_H5y_S_G2x2z_S_vrr+WPZ*I_ERI_H5y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H5y_S_F2xz_S_M1_vrr;
      Double I_ERI_I4y2z_S_G2x2z_S_vrr = PAZ*I_ERI_H4yz_S_G2x2z_S_vrr+WPZ*I_ERI_H4yz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_G4y_S_G2x2z_S_vrr-rhod2zsq*I_ERI_G4y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H4yz_S_F2xz_S_M1_vrr;
      Double I_ERI_I3y3z_S_G2x2z_S_vrr = PAZ*I_ERI_H3y2z_S_G2x2z_S_vrr+WPZ*I_ERI_H3y2z_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H3y2z_S_F2xz_S_M1_vrr;
      Double I_ERI_I2y4z_S_G2x2z_S_vrr = PAY*I_ERI_Hy4z_S_G2x2z_S_vrr+WPY*I_ERI_Hy4z_S_G2x2z_S_M1_vrr+oned2z*I_ERI_G4z_S_G2x2z_S_vrr-rhod2zsq*I_ERI_G4z_S_G2x2z_S_M1_vrr;
      Double I_ERI_Iy5z_S_G2x2z_S_vrr = PAY*I_ERI_H5z_S_G2x2z_S_vrr+WPY*I_ERI_H5z_S_G2x2z_S_M1_vrr;
      Double I_ERI_I6z_S_G2x2z_S_vrr = PAZ*I_ERI_H5z_S_G2x2z_S_vrr+WPZ*I_ERI_H5z_S_G2x2z_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G2x2z_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_H5z_S_F2xz_S_M1_vrr;
      Double I_ERI_I6x_S_Gx3y_S_vrr = PAX*I_ERI_H5x_S_Gx3y_S_vrr+WPX*I_ERI_H5x_S_Gx3y_S_M1_vrr+5*oned2z*I_ERI_G4x_S_Gx3y_S_vrr-5*rhod2zsq*I_ERI_G4x_S_Gx3y_S_M1_vrr+oned2k*I_ERI_H5x_S_F3y_S_M1_vrr;
      Double I_ERI_I5xy_S_Gx3y_S_vrr = PAY*I_ERI_H5x_S_Gx3y_S_vrr+WPY*I_ERI_H5x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_H5x_S_Fx2y_S_M1_vrr;
      Double I_ERI_I5xz_S_Gx3y_S_vrr = PAZ*I_ERI_H5x_S_Gx3y_S_vrr+WPZ*I_ERI_H5x_S_Gx3y_S_M1_vrr;
      Double I_ERI_I4x2y_S_Gx3y_S_vrr = PAY*I_ERI_H4xy_S_Gx3y_S_vrr+WPY*I_ERI_H4xy_S_Gx3y_S_M1_vrr+oned2z*I_ERI_G4x_S_Gx3y_S_vrr-rhod2zsq*I_ERI_G4x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_H4xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_I4xyz_S_Gx3y_S_vrr = PAZ*I_ERI_H4xy_S_Gx3y_S_vrr+WPZ*I_ERI_H4xy_S_Gx3y_S_M1_vrr;
      Double I_ERI_I4x2z_S_Gx3y_S_vrr = PAZ*I_ERI_H4xz_S_Gx3y_S_vrr+WPZ*I_ERI_H4xz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_G4x_S_Gx3y_S_vrr-rhod2zsq*I_ERI_G4x_S_Gx3y_S_M1_vrr;
      Double I_ERI_I3x3y_S_Gx3y_S_vrr = PAY*I_ERI_H3x2y_S_Gx3y_S_vrr+WPY*I_ERI_H3x2y_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_H3x2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Gx3y_S_vrr = PAZ*I_ERI_H3x2y_S_Gx3y_S_vrr+WPZ*I_ERI_H3x2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Gx3y_S_vrr = PAY*I_ERI_H3x2z_S_Gx3y_S_vrr+WPY*I_ERI_H3x2z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_H3x2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_I3x3z_S_Gx3y_S_vrr = PAZ*I_ERI_H3x2z_S_Gx3y_S_vrr+WPZ*I_ERI_H3x2z_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_Gx3y_S_M1_vrr;
      Double I_ERI_I2x4y_S_Gx3y_S_vrr = PAX*I_ERI_Hx4y_S_Gx3y_S_vrr+WPX*I_ERI_Hx4y_S_Gx3y_S_M1_vrr+oned2z*I_ERI_G4y_S_Gx3y_S_vrr-rhod2zsq*I_ERI_G4y_S_Gx3y_S_M1_vrr+oned2k*I_ERI_Hx4y_S_F3y_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Gx3y_S_vrr = PAZ*I_ERI_H2x3y_S_Gx3y_S_vrr+WPZ*I_ERI_H2x3y_S_Gx3y_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Gx3y_S_vrr = PAZ*I_ERI_H2x2yz_S_Gx3y_S_vrr+WPZ*I_ERI_H2x2yz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_G2x2y_S_Gx3y_S_vrr-rhod2zsq*I_ERI_G2x2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Gx3y_S_vrr = PAY*I_ERI_H2x3z_S_Gx3y_S_vrr+WPY*I_ERI_H2x3z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_H2x3z_S_Fx2y_S_M1_vrr;
      Double I_ERI_I2x4z_S_Gx3y_S_vrr = PAX*I_ERI_Hx4z_S_Gx3y_S_vrr+WPX*I_ERI_Hx4z_S_Gx3y_S_M1_vrr+oned2z*I_ERI_G4z_S_Gx3y_S_vrr-rhod2zsq*I_ERI_G4z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_Hx4z_S_F3y_S_M1_vrr;
      Double I_ERI_Ix5y_S_Gx3y_S_vrr = PAX*I_ERI_H5y_S_Gx3y_S_vrr+WPX*I_ERI_H5y_S_Gx3y_S_M1_vrr+oned2k*I_ERI_H5y_S_F3y_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Gx3y_S_vrr = PAZ*I_ERI_Hx4y_S_Gx3y_S_vrr+WPZ*I_ERI_Hx4y_S_Gx3y_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Gx3y_S_vrr = PAX*I_ERI_H3y2z_S_Gx3y_S_vrr+WPX*I_ERI_H3y2z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_H3y2z_S_F3y_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Gx3y_S_vrr = PAX*I_ERI_H2y3z_S_Gx3y_S_vrr+WPX*I_ERI_H2y3z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_H2y3z_S_F3y_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Gx3y_S_vrr = PAY*I_ERI_Hx4z_S_Gx3y_S_vrr+WPY*I_ERI_Hx4z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Hx4z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Ix5z_S_Gx3y_S_vrr = PAX*I_ERI_H5z_S_Gx3y_S_vrr+WPX*I_ERI_H5z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_H5z_S_F3y_S_M1_vrr;
      Double I_ERI_I6y_S_Gx3y_S_vrr = PAY*I_ERI_H5y_S_Gx3y_S_vrr+WPY*I_ERI_H5y_S_Gx3y_S_M1_vrr+5*oned2z*I_ERI_G4y_S_Gx3y_S_vrr-5*rhod2zsq*I_ERI_G4y_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_H5y_S_Fx2y_S_M1_vrr;
      Double I_ERI_I5yz_S_Gx3y_S_vrr = PAZ*I_ERI_H5y_S_Gx3y_S_vrr+WPZ*I_ERI_H5y_S_Gx3y_S_M1_vrr;
      Double I_ERI_I4y2z_S_Gx3y_S_vrr = PAZ*I_ERI_H4yz_S_Gx3y_S_vrr+WPZ*I_ERI_H4yz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_G4y_S_Gx3y_S_vrr-rhod2zsq*I_ERI_G4y_S_Gx3y_S_M1_vrr;
      Double I_ERI_I3y3z_S_Gx3y_S_vrr = PAZ*I_ERI_H3y2z_S_Gx3y_S_vrr+WPZ*I_ERI_H3y2z_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_Gx3y_S_M1_vrr;
      Double I_ERI_I2y4z_S_Gx3y_S_vrr = PAY*I_ERI_Hy4z_S_Gx3y_S_vrr+WPY*I_ERI_Hy4z_S_Gx3y_S_M1_vrr+oned2z*I_ERI_G4z_S_Gx3y_S_vrr-rhod2zsq*I_ERI_G4z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Hy4z_S_Fx2y_S_M1_vrr;
      Double I_ERI_Iy5z_S_Gx3y_S_vrr = PAY*I_ERI_H5z_S_Gx3y_S_vrr+WPY*I_ERI_H5z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_H5z_S_Fx2y_S_M1_vrr;
      Double I_ERI_I6z_S_Gx3y_S_vrr = PAZ*I_ERI_H5z_S_Gx3y_S_vrr+WPZ*I_ERI_H5z_S_Gx3y_S_M1_vrr+5*oned2z*I_ERI_G4z_S_Gx3y_S_vrr-5*rhod2zsq*I_ERI_G4z_S_Gx3y_S_M1_vrr;
      Double I_ERI_I6x_S_Gx2yz_S_vrr = PAX*I_ERI_H5x_S_Gx2yz_S_vrr+WPX*I_ERI_H5x_S_Gx2yz_S_M1_vrr+5*oned2z*I_ERI_G4x_S_Gx2yz_S_vrr-5*rhod2zsq*I_ERI_G4x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H5x_S_F2yz_S_M1_vrr;
      Double I_ERI_I5xy_S_Gx2yz_S_vrr = PAY*I_ERI_H5x_S_Gx2yz_S_vrr+WPY*I_ERI_H5x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_H5x_S_Fxyz_S_M1_vrr;
      Double I_ERI_I5xz_S_Gx2yz_S_vrr = PAZ*I_ERI_H5x_S_Gx2yz_S_vrr+WPZ*I_ERI_H5x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H5x_S_Fx2y_S_M1_vrr;
      Double I_ERI_I4x2y_S_Gx2yz_S_vrr = PAY*I_ERI_H4xy_S_Gx2yz_S_vrr+WPY*I_ERI_H4xy_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_G4x_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_G4x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_H4xy_S_Fxyz_S_M1_vrr;
      Double I_ERI_I4xyz_S_Gx2yz_S_vrr = PAZ*I_ERI_H4xy_S_Gx2yz_S_vrr+WPZ*I_ERI_H4xy_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H4xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_I4x2z_S_Gx2yz_S_vrr = PAZ*I_ERI_H4xz_S_Gx2yz_S_vrr+WPZ*I_ERI_H4xz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_G4x_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_G4x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H4xz_S_Fx2y_S_M1_vrr;
      Double I_ERI_I3x3y_S_Gx2yz_S_vrr = PAY*I_ERI_H3x2y_S_Gx2yz_S_vrr+WPY*I_ERI_H3x2y_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_H3x2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Gx2yz_S_vrr = PAZ*I_ERI_H3x2y_S_Gx2yz_S_vrr+WPZ*I_ERI_H3x2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H3x2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Gx2yz_S_vrr = PAY*I_ERI_H3x2z_S_Gx2yz_S_vrr+WPY*I_ERI_H3x2z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_H3x2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_I3x3z_S_Gx2yz_S_vrr = PAZ*I_ERI_H3x2z_S_Gx2yz_S_vrr+WPZ*I_ERI_H3x2z_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H3x2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_I2x4y_S_Gx2yz_S_vrr = PAX*I_ERI_Hx4y_S_Gx2yz_S_vrr+WPX*I_ERI_Hx4y_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_G4y_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_G4y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Hx4y_S_F2yz_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Gx2yz_S_vrr = PAZ*I_ERI_H2x3y_S_Gx2yz_S_vrr+WPZ*I_ERI_H2x3y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H2x3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Gx2yz_S_vrr = PAZ*I_ERI_H2x2yz_S_Gx2yz_S_vrr+WPZ*I_ERI_H2x2yz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_G2x2y_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_G2x2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H2x2yz_S_Fx2y_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Gx2yz_S_vrr = PAY*I_ERI_H2x3z_S_Gx2yz_S_vrr+WPY*I_ERI_H2x3z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_H2x3z_S_Fxyz_S_M1_vrr;
      Double I_ERI_I2x4z_S_Gx2yz_S_vrr = PAX*I_ERI_Hx4z_S_Gx2yz_S_vrr+WPX*I_ERI_Hx4z_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_G4z_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_G4z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Hx4z_S_F2yz_S_M1_vrr;
      Double I_ERI_Ix5y_S_Gx2yz_S_vrr = PAX*I_ERI_H5y_S_Gx2yz_S_vrr+WPX*I_ERI_H5y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H5y_S_F2yz_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Gx2yz_S_vrr = PAZ*I_ERI_Hx4y_S_Gx2yz_S_vrr+WPZ*I_ERI_Hx4y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Hx4y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Gx2yz_S_vrr = PAX*I_ERI_H3y2z_S_Gx2yz_S_vrr+WPX*I_ERI_H3y2z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H3y2z_S_F2yz_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Gx2yz_S_vrr = PAX*I_ERI_H2y3z_S_Gx2yz_S_vrr+WPX*I_ERI_H2y3z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H2y3z_S_F2yz_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Gx2yz_S_vrr = PAY*I_ERI_Hx4z_S_Gx2yz_S_vrr+WPY*I_ERI_Hx4z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Hx4z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Ix5z_S_Gx2yz_S_vrr = PAX*I_ERI_H5z_S_Gx2yz_S_vrr+WPX*I_ERI_H5z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H5z_S_F2yz_S_M1_vrr;
      Double I_ERI_I6y_S_Gx2yz_S_vrr = PAY*I_ERI_H5y_S_Gx2yz_S_vrr+WPY*I_ERI_H5y_S_Gx2yz_S_M1_vrr+5*oned2z*I_ERI_G4y_S_Gx2yz_S_vrr-5*rhod2zsq*I_ERI_G4y_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_H5y_S_Fxyz_S_M1_vrr;
      Double I_ERI_I5yz_S_Gx2yz_S_vrr = PAZ*I_ERI_H5y_S_Gx2yz_S_vrr+WPZ*I_ERI_H5y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H5y_S_Fx2y_S_M1_vrr;
      Double I_ERI_I4y2z_S_Gx2yz_S_vrr = PAZ*I_ERI_H4yz_S_Gx2yz_S_vrr+WPZ*I_ERI_H4yz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_G4y_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_G4y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H4yz_S_Fx2y_S_M1_vrr;
      Double I_ERI_I3y3z_S_Gx2yz_S_vrr = PAZ*I_ERI_H3y2z_S_Gx2yz_S_vrr+WPZ*I_ERI_H3y2z_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H3y2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_I2y4z_S_Gx2yz_S_vrr = PAY*I_ERI_Hy4z_S_Gx2yz_S_vrr+WPY*I_ERI_Hy4z_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_G4z_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_G4z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Hy4z_S_Fxyz_S_M1_vrr;
      Double I_ERI_Iy5z_S_Gx2yz_S_vrr = PAY*I_ERI_H5z_S_Gx2yz_S_vrr+WPY*I_ERI_H5z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_H5z_S_Fxyz_S_M1_vrr;
      Double I_ERI_I6z_S_Gx2yz_S_vrr = PAZ*I_ERI_H5z_S_Gx2yz_S_vrr+WPZ*I_ERI_H5z_S_Gx2yz_S_M1_vrr+5*oned2z*I_ERI_G4z_S_Gx2yz_S_vrr-5*rhod2zsq*I_ERI_G4z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_H5z_S_Fx2y_S_M1_vrr;
      Double I_ERI_I6x_S_Gxy2z_S_vrr = PAX*I_ERI_H5x_S_Gxy2z_S_vrr+WPX*I_ERI_H5x_S_Gxy2z_S_M1_vrr+5*oned2z*I_ERI_G4x_S_Gxy2z_S_vrr-5*rhod2zsq*I_ERI_G4x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H5x_S_Fy2z_S_M1_vrr;
      Double I_ERI_I5xy_S_Gxy2z_S_vrr = PAY*I_ERI_H5x_S_Gxy2z_S_vrr+WPY*I_ERI_H5x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H5x_S_Fx2z_S_M1_vrr;
      Double I_ERI_I5xz_S_Gxy2z_S_vrr = PAZ*I_ERI_H5x_S_Gxy2z_S_vrr+WPZ*I_ERI_H5x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H5x_S_Fxyz_S_M1_vrr;
      Double I_ERI_I4x2y_S_Gxy2z_S_vrr = PAY*I_ERI_H4xy_S_Gxy2z_S_vrr+WPY*I_ERI_H4xy_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_G4x_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_G4x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H4xy_S_Fx2z_S_M1_vrr;
      Double I_ERI_I4xyz_S_Gxy2z_S_vrr = PAZ*I_ERI_H4xy_S_Gxy2z_S_vrr+WPZ*I_ERI_H4xy_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H4xy_S_Fxyz_S_M1_vrr;
      Double I_ERI_I4x2z_S_Gxy2z_S_vrr = PAZ*I_ERI_H4xz_S_Gxy2z_S_vrr+WPZ*I_ERI_H4xz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_G4x_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_G4x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H4xz_S_Fxyz_S_M1_vrr;
      Double I_ERI_I3x3y_S_Gxy2z_S_vrr = PAY*I_ERI_H3x2y_S_Gxy2z_S_vrr+WPY*I_ERI_H3x2y_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H3x2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Gxy2z_S_vrr = PAZ*I_ERI_H3x2y_S_Gxy2z_S_vrr+WPZ*I_ERI_H3x2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H3x2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Gxy2z_S_vrr = PAY*I_ERI_H3x2z_S_Gxy2z_S_vrr+WPY*I_ERI_H3x2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H3x2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_I3x3z_S_Gxy2z_S_vrr = PAZ*I_ERI_H3x2z_S_Gxy2z_S_vrr+WPZ*I_ERI_H3x2z_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H3x2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_I2x4y_S_Gxy2z_S_vrr = PAX*I_ERI_Hx4y_S_Gxy2z_S_vrr+WPX*I_ERI_Hx4y_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_G4y_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_G4y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Hx4y_S_Fy2z_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Gxy2z_S_vrr = PAZ*I_ERI_H2x3y_S_Gxy2z_S_vrr+WPZ*I_ERI_H2x3y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H2x3y_S_Fxyz_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Gxy2z_S_vrr = PAZ*I_ERI_H2x2yz_S_Gxy2z_S_vrr+WPZ*I_ERI_H2x2yz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_G2x2y_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_G2x2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H2x2yz_S_Fxyz_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Gxy2z_S_vrr = PAY*I_ERI_H2x3z_S_Gxy2z_S_vrr+WPY*I_ERI_H2x3z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H2x3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_I2x4z_S_Gxy2z_S_vrr = PAX*I_ERI_Hx4z_S_Gxy2z_S_vrr+WPX*I_ERI_Hx4z_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_G4z_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_G4z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Hx4z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Ix5y_S_Gxy2z_S_vrr = PAX*I_ERI_H5y_S_Gxy2z_S_vrr+WPX*I_ERI_H5y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H5y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Gxy2z_S_vrr = PAZ*I_ERI_Hx4y_S_Gxy2z_S_vrr+WPZ*I_ERI_Hx4y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Hx4y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Gxy2z_S_vrr = PAX*I_ERI_H3y2z_S_Gxy2z_S_vrr+WPX*I_ERI_H3y2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H3y2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Gxy2z_S_vrr = PAX*I_ERI_H2y3z_S_Gxy2z_S_vrr+WPX*I_ERI_H2y3z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H2y3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Gxy2z_S_vrr = PAY*I_ERI_Hx4z_S_Gxy2z_S_vrr+WPY*I_ERI_Hx4z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Hx4z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Ix5z_S_Gxy2z_S_vrr = PAX*I_ERI_H5z_S_Gxy2z_S_vrr+WPX*I_ERI_H5z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H5z_S_Fy2z_S_M1_vrr;
      Double I_ERI_I6y_S_Gxy2z_S_vrr = PAY*I_ERI_H5y_S_Gxy2z_S_vrr+WPY*I_ERI_H5y_S_Gxy2z_S_M1_vrr+5*oned2z*I_ERI_G4y_S_Gxy2z_S_vrr-5*rhod2zsq*I_ERI_G4y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H5y_S_Fx2z_S_M1_vrr;
      Double I_ERI_I5yz_S_Gxy2z_S_vrr = PAZ*I_ERI_H5y_S_Gxy2z_S_vrr+WPZ*I_ERI_H5y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H5y_S_Fxyz_S_M1_vrr;
      Double I_ERI_I4y2z_S_Gxy2z_S_vrr = PAZ*I_ERI_H4yz_S_Gxy2z_S_vrr+WPZ*I_ERI_H4yz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_G4y_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_G4y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H4yz_S_Fxyz_S_M1_vrr;
      Double I_ERI_I3y3z_S_Gxy2z_S_vrr = PAZ*I_ERI_H3y2z_S_Gxy2z_S_vrr+WPZ*I_ERI_H3y2z_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H3y2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_I2y4z_S_Gxy2z_S_vrr = PAY*I_ERI_Hy4z_S_Gxy2z_S_vrr+WPY*I_ERI_Hy4z_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_G4z_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_G4z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Hy4z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Iy5z_S_Gxy2z_S_vrr = PAY*I_ERI_H5z_S_Gxy2z_S_vrr+WPY*I_ERI_H5z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_H5z_S_Fx2z_S_M1_vrr;
      Double I_ERI_I6z_S_Gxy2z_S_vrr = PAZ*I_ERI_H5z_S_Gxy2z_S_vrr+WPZ*I_ERI_H5z_S_Gxy2z_S_M1_vrr+5*oned2z*I_ERI_G4z_S_Gxy2z_S_vrr-5*rhod2zsq*I_ERI_G4z_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_H5z_S_Fxyz_S_M1_vrr;
      Double I_ERI_I6x_S_Gx3z_S_vrr = PAX*I_ERI_H5x_S_Gx3z_S_vrr+WPX*I_ERI_H5x_S_Gx3z_S_M1_vrr+5*oned2z*I_ERI_G4x_S_Gx3z_S_vrr-5*rhod2zsq*I_ERI_G4x_S_Gx3z_S_M1_vrr+oned2k*I_ERI_H5x_S_F3z_S_M1_vrr;
      Double I_ERI_I5xy_S_Gx3z_S_vrr = PAY*I_ERI_H5x_S_Gx3z_S_vrr+WPY*I_ERI_H5x_S_Gx3z_S_M1_vrr;
      Double I_ERI_I5xz_S_Gx3z_S_vrr = PAZ*I_ERI_H5x_S_Gx3z_S_vrr+WPZ*I_ERI_H5x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H5x_S_Fx2z_S_M1_vrr;
      Double I_ERI_I4x2y_S_Gx3z_S_vrr = PAY*I_ERI_H4xy_S_Gx3z_S_vrr+WPY*I_ERI_H4xy_S_Gx3z_S_M1_vrr+oned2z*I_ERI_G4x_S_Gx3z_S_vrr-rhod2zsq*I_ERI_G4x_S_Gx3z_S_M1_vrr;
      Double I_ERI_I4xyz_S_Gx3z_S_vrr = PAZ*I_ERI_H4xy_S_Gx3z_S_vrr+WPZ*I_ERI_H4xy_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H4xy_S_Fx2z_S_M1_vrr;
      Double I_ERI_I4x2z_S_Gx3z_S_vrr = PAZ*I_ERI_H4xz_S_Gx3z_S_vrr+WPZ*I_ERI_H4xz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_G4x_S_Gx3z_S_vrr-rhod2zsq*I_ERI_G4x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H4xz_S_Fx2z_S_M1_vrr;
      Double I_ERI_I3x3y_S_Gx3z_S_vrr = PAY*I_ERI_H3x2y_S_Gx3z_S_vrr+WPY*I_ERI_H3x2y_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_Gx3z_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Gx3z_S_vrr = PAZ*I_ERI_H3x2y_S_Gx3z_S_vrr+WPZ*I_ERI_H3x2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H3x2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Gx3z_S_vrr = PAY*I_ERI_H3x2z_S_Gx3z_S_vrr+WPY*I_ERI_H3x2z_S_Gx3z_S_M1_vrr;
      Double I_ERI_I3x3z_S_Gx3z_S_vrr = PAZ*I_ERI_H3x2z_S_Gx3z_S_vrr+WPZ*I_ERI_H3x2z_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H3x2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_I2x4y_S_Gx3z_S_vrr = PAX*I_ERI_Hx4y_S_Gx3z_S_vrr+WPX*I_ERI_Hx4y_S_Gx3z_S_M1_vrr+oned2z*I_ERI_G4y_S_Gx3z_S_vrr-rhod2zsq*I_ERI_G4y_S_Gx3z_S_M1_vrr+oned2k*I_ERI_Hx4y_S_F3z_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Gx3z_S_vrr = PAZ*I_ERI_H2x3y_S_Gx3z_S_vrr+WPZ*I_ERI_H2x3y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H2x3y_S_Fx2z_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Gx3z_S_vrr = PAZ*I_ERI_H2x2yz_S_Gx3z_S_vrr+WPZ*I_ERI_H2x2yz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_G2x2y_S_Gx3z_S_vrr-rhod2zsq*I_ERI_G2x2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H2x2yz_S_Fx2z_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Gx3z_S_vrr = PAY*I_ERI_H2x3z_S_Gx3z_S_vrr+WPY*I_ERI_H2x3z_S_Gx3z_S_M1_vrr;
      Double I_ERI_I2x4z_S_Gx3z_S_vrr = PAX*I_ERI_Hx4z_S_Gx3z_S_vrr+WPX*I_ERI_Hx4z_S_Gx3z_S_M1_vrr+oned2z*I_ERI_G4z_S_Gx3z_S_vrr-rhod2zsq*I_ERI_G4z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_Hx4z_S_F3z_S_M1_vrr;
      Double I_ERI_Ix5y_S_Gx3z_S_vrr = PAX*I_ERI_H5y_S_Gx3z_S_vrr+WPX*I_ERI_H5y_S_Gx3z_S_M1_vrr+oned2k*I_ERI_H5y_S_F3z_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Gx3z_S_vrr = PAZ*I_ERI_Hx4y_S_Gx3z_S_vrr+WPZ*I_ERI_Hx4y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Hx4y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Gx3z_S_vrr = PAX*I_ERI_H3y2z_S_Gx3z_S_vrr+WPX*I_ERI_H3y2z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_H3y2z_S_F3z_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Gx3z_S_vrr = PAX*I_ERI_H2y3z_S_Gx3z_S_vrr+WPX*I_ERI_H2y3z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_H2y3z_S_F3z_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Gx3z_S_vrr = PAY*I_ERI_Hx4z_S_Gx3z_S_vrr+WPY*I_ERI_Hx4z_S_Gx3z_S_M1_vrr;
      Double I_ERI_Ix5z_S_Gx3z_S_vrr = PAX*I_ERI_H5z_S_Gx3z_S_vrr+WPX*I_ERI_H5z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_H5z_S_F3z_S_M1_vrr;
      Double I_ERI_I6y_S_Gx3z_S_vrr = PAY*I_ERI_H5y_S_Gx3z_S_vrr+WPY*I_ERI_H5y_S_Gx3z_S_M1_vrr+5*oned2z*I_ERI_G4y_S_Gx3z_S_vrr-5*rhod2zsq*I_ERI_G4y_S_Gx3z_S_M1_vrr;
      Double I_ERI_I5yz_S_Gx3z_S_vrr = PAZ*I_ERI_H5y_S_Gx3z_S_vrr+WPZ*I_ERI_H5y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H5y_S_Fx2z_S_M1_vrr;
      Double I_ERI_I4y2z_S_Gx3z_S_vrr = PAZ*I_ERI_H4yz_S_Gx3z_S_vrr+WPZ*I_ERI_H4yz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_G4y_S_Gx3z_S_vrr-rhod2zsq*I_ERI_G4y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H4yz_S_Fx2z_S_M1_vrr;
      Double I_ERI_I3y3z_S_Gx3z_S_vrr = PAZ*I_ERI_H3y2z_S_Gx3z_S_vrr+WPZ*I_ERI_H3y2z_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H3y2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_I2y4z_S_Gx3z_S_vrr = PAY*I_ERI_Hy4z_S_Gx3z_S_vrr+WPY*I_ERI_Hy4z_S_Gx3z_S_M1_vrr+oned2z*I_ERI_G4z_S_Gx3z_S_vrr-rhod2zsq*I_ERI_G4z_S_Gx3z_S_M1_vrr;
      Double I_ERI_Iy5z_S_Gx3z_S_vrr = PAY*I_ERI_H5z_S_Gx3z_S_vrr+WPY*I_ERI_H5z_S_Gx3z_S_M1_vrr;
      Double I_ERI_I6z_S_Gx3z_S_vrr = PAZ*I_ERI_H5z_S_Gx3z_S_vrr+WPZ*I_ERI_H5z_S_Gx3z_S_M1_vrr+5*oned2z*I_ERI_G4z_S_Gx3z_S_vrr-5*rhod2zsq*I_ERI_G4z_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_H5z_S_Fx2z_S_M1_vrr;
      Double I_ERI_I6x_S_G4y_S_vrr = PAX*I_ERI_H5x_S_G4y_S_vrr+WPX*I_ERI_H5x_S_G4y_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G4y_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G4y_S_M1_vrr;
      Double I_ERI_I5xy_S_G4y_S_vrr = PAY*I_ERI_H5x_S_G4y_S_vrr+WPY*I_ERI_H5x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_H5x_S_F3y_S_M1_vrr;
      Double I_ERI_I5xz_S_G4y_S_vrr = PAZ*I_ERI_H5x_S_G4y_S_vrr+WPZ*I_ERI_H5x_S_G4y_S_M1_vrr;
      Double I_ERI_I4x2y_S_G4y_S_vrr = PAY*I_ERI_H4xy_S_G4y_S_vrr+WPY*I_ERI_H4xy_S_G4y_S_M1_vrr+oned2z*I_ERI_G4x_S_G4y_S_vrr-rhod2zsq*I_ERI_G4x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_H4xy_S_F3y_S_M1_vrr;
      Double I_ERI_I4xyz_S_G4y_S_vrr = PAZ*I_ERI_H4xy_S_G4y_S_vrr+WPZ*I_ERI_H4xy_S_G4y_S_M1_vrr;
      Double I_ERI_I4x2z_S_G4y_S_vrr = PAZ*I_ERI_H4xz_S_G4y_S_vrr+WPZ*I_ERI_H4xz_S_G4y_S_M1_vrr+oned2z*I_ERI_G4x_S_G4y_S_vrr-rhod2zsq*I_ERI_G4x_S_G4y_S_M1_vrr;
      Double I_ERI_I3x3y_S_G4y_S_vrr = PAY*I_ERI_H3x2y_S_G4y_S_vrr+WPY*I_ERI_H3x2y_S_G4y_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G4y_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G4y_S_M1_vrr+4*oned2k*I_ERI_H3x2y_S_F3y_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G4y_S_vrr = PAZ*I_ERI_H3x2y_S_G4y_S_vrr+WPZ*I_ERI_H3x2y_S_G4y_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G4y_S_vrr = PAY*I_ERI_H3x2z_S_G4y_S_vrr+WPY*I_ERI_H3x2z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_H3x2z_S_F3y_S_M1_vrr;
      Double I_ERI_I3x3z_S_G4y_S_vrr = PAZ*I_ERI_H3x2z_S_G4y_S_vrr+WPZ*I_ERI_H3x2z_S_G4y_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G4y_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G4y_S_M1_vrr;
      Double I_ERI_I2x4y_S_G4y_S_vrr = PAX*I_ERI_Hx4y_S_G4y_S_vrr+WPX*I_ERI_Hx4y_S_G4y_S_M1_vrr+oned2z*I_ERI_G4y_S_G4y_S_vrr-rhod2zsq*I_ERI_G4y_S_G4y_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G4y_S_vrr = PAZ*I_ERI_H2x3y_S_G4y_S_vrr+WPZ*I_ERI_H2x3y_S_G4y_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G4y_S_vrr = PAZ*I_ERI_H2x2yz_S_G4y_S_vrr+WPZ*I_ERI_H2x2yz_S_G4y_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G4y_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G4y_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G4y_S_vrr = PAY*I_ERI_H2x3z_S_G4y_S_vrr+WPY*I_ERI_H2x3z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_H2x3z_S_F3y_S_M1_vrr;
      Double I_ERI_I2x4z_S_G4y_S_vrr = PAX*I_ERI_Hx4z_S_G4y_S_vrr+WPX*I_ERI_Hx4z_S_G4y_S_M1_vrr+oned2z*I_ERI_G4z_S_G4y_S_vrr-rhod2zsq*I_ERI_G4z_S_G4y_S_M1_vrr;
      Double I_ERI_Ix5y_S_G4y_S_vrr = PAX*I_ERI_H5y_S_G4y_S_vrr+WPX*I_ERI_H5y_S_G4y_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G4y_S_vrr = PAZ*I_ERI_Hx4y_S_G4y_S_vrr+WPZ*I_ERI_Hx4y_S_G4y_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G4y_S_vrr = PAX*I_ERI_H3y2z_S_G4y_S_vrr+WPX*I_ERI_H3y2z_S_G4y_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G4y_S_vrr = PAX*I_ERI_H2y3z_S_G4y_S_vrr+WPX*I_ERI_H2y3z_S_G4y_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G4y_S_vrr = PAY*I_ERI_Hx4z_S_G4y_S_vrr+WPY*I_ERI_Hx4z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Hx4z_S_F3y_S_M1_vrr;
      Double I_ERI_Ix5z_S_G4y_S_vrr = PAX*I_ERI_H5z_S_G4y_S_vrr+WPX*I_ERI_H5z_S_G4y_S_M1_vrr;
      Double I_ERI_I6y_S_G4y_S_vrr = PAY*I_ERI_H5y_S_G4y_S_vrr+WPY*I_ERI_H5y_S_G4y_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G4y_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G4y_S_M1_vrr+4*oned2k*I_ERI_H5y_S_F3y_S_M1_vrr;
      Double I_ERI_I5yz_S_G4y_S_vrr = PAZ*I_ERI_H5y_S_G4y_S_vrr+WPZ*I_ERI_H5y_S_G4y_S_M1_vrr;
      Double I_ERI_I4y2z_S_G4y_S_vrr = PAZ*I_ERI_H4yz_S_G4y_S_vrr+WPZ*I_ERI_H4yz_S_G4y_S_M1_vrr+oned2z*I_ERI_G4y_S_G4y_S_vrr-rhod2zsq*I_ERI_G4y_S_G4y_S_M1_vrr;
      Double I_ERI_I3y3z_S_G4y_S_vrr = PAZ*I_ERI_H3y2z_S_G4y_S_vrr+WPZ*I_ERI_H3y2z_S_G4y_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G4y_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G4y_S_M1_vrr;
      Double I_ERI_I2y4z_S_G4y_S_vrr = PAY*I_ERI_Hy4z_S_G4y_S_vrr+WPY*I_ERI_Hy4z_S_G4y_S_M1_vrr+oned2z*I_ERI_G4z_S_G4y_S_vrr-rhod2zsq*I_ERI_G4z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Hy4z_S_F3y_S_M1_vrr;
      Double I_ERI_Iy5z_S_G4y_S_vrr = PAY*I_ERI_H5z_S_G4y_S_vrr+WPY*I_ERI_H5z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_H5z_S_F3y_S_M1_vrr;
      Double I_ERI_I6z_S_G4y_S_vrr = PAZ*I_ERI_H5z_S_G4y_S_vrr+WPZ*I_ERI_H5z_S_G4y_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G4y_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G4y_S_M1_vrr;
      Double I_ERI_I6x_S_G3yz_S_vrr = PAX*I_ERI_H5x_S_G3yz_S_vrr+WPX*I_ERI_H5x_S_G3yz_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G3yz_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G3yz_S_M1_vrr;
      Double I_ERI_I5xy_S_G3yz_S_vrr = PAY*I_ERI_H5x_S_G3yz_S_vrr+WPY*I_ERI_H5x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_H5x_S_F2yz_S_M1_vrr;
      Double I_ERI_I5xz_S_G3yz_S_vrr = PAZ*I_ERI_H5x_S_G3yz_S_vrr+WPZ*I_ERI_H5x_S_G3yz_S_M1_vrr+oned2k*I_ERI_H5x_S_F3y_S_M1_vrr;
      Double I_ERI_I4x2y_S_G3yz_S_vrr = PAY*I_ERI_H4xy_S_G3yz_S_vrr+WPY*I_ERI_H4xy_S_G3yz_S_M1_vrr+oned2z*I_ERI_G4x_S_G3yz_S_vrr-rhod2zsq*I_ERI_G4x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_H4xy_S_F2yz_S_M1_vrr;
      Double I_ERI_I4xyz_S_G3yz_S_vrr = PAZ*I_ERI_H4xy_S_G3yz_S_vrr+WPZ*I_ERI_H4xy_S_G3yz_S_M1_vrr+oned2k*I_ERI_H4xy_S_F3y_S_M1_vrr;
      Double I_ERI_I4x2z_S_G3yz_S_vrr = PAZ*I_ERI_H4xz_S_G3yz_S_vrr+WPZ*I_ERI_H4xz_S_G3yz_S_M1_vrr+oned2z*I_ERI_G4x_S_G3yz_S_vrr-rhod2zsq*I_ERI_G4x_S_G3yz_S_M1_vrr+oned2k*I_ERI_H4xz_S_F3y_S_M1_vrr;
      Double I_ERI_I3x3y_S_G3yz_S_vrr = PAY*I_ERI_H3x2y_S_G3yz_S_vrr+WPY*I_ERI_H3x2y_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_H3x2y_S_F2yz_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G3yz_S_vrr = PAZ*I_ERI_H3x2y_S_G3yz_S_vrr+WPZ*I_ERI_H3x2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_H3x2y_S_F3y_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G3yz_S_vrr = PAY*I_ERI_H3x2z_S_G3yz_S_vrr+WPY*I_ERI_H3x2z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_H3x2z_S_F2yz_S_M1_vrr;
      Double I_ERI_I3x3z_S_G3yz_S_vrr = PAZ*I_ERI_H3x2z_S_G3yz_S_vrr+WPZ*I_ERI_H3x2z_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G3yz_S_M1_vrr+oned2k*I_ERI_H3x2z_S_F3y_S_M1_vrr;
      Double I_ERI_I2x4y_S_G3yz_S_vrr = PAX*I_ERI_Hx4y_S_G3yz_S_vrr+WPX*I_ERI_Hx4y_S_G3yz_S_M1_vrr+oned2z*I_ERI_G4y_S_G3yz_S_vrr-rhod2zsq*I_ERI_G4y_S_G3yz_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G3yz_S_vrr = PAZ*I_ERI_H2x3y_S_G3yz_S_vrr+WPZ*I_ERI_H2x3y_S_G3yz_S_M1_vrr+oned2k*I_ERI_H2x3y_S_F3y_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G3yz_S_vrr = PAZ*I_ERI_H2x2yz_S_G3yz_S_vrr+WPZ*I_ERI_H2x2yz_S_G3yz_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G3yz_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_H2x2yz_S_F3y_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G3yz_S_vrr = PAY*I_ERI_H2x3z_S_G3yz_S_vrr+WPY*I_ERI_H2x3z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_H2x3z_S_F2yz_S_M1_vrr;
      Double I_ERI_I2x4z_S_G3yz_S_vrr = PAX*I_ERI_Hx4z_S_G3yz_S_vrr+WPX*I_ERI_Hx4z_S_G3yz_S_M1_vrr+oned2z*I_ERI_G4z_S_G3yz_S_vrr-rhod2zsq*I_ERI_G4z_S_G3yz_S_M1_vrr;
      Double I_ERI_Ix5y_S_G3yz_S_vrr = PAX*I_ERI_H5y_S_G3yz_S_vrr+WPX*I_ERI_H5y_S_G3yz_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G3yz_S_vrr = PAZ*I_ERI_Hx4y_S_G3yz_S_vrr+WPZ*I_ERI_Hx4y_S_G3yz_S_M1_vrr+oned2k*I_ERI_Hx4y_S_F3y_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G3yz_S_vrr = PAX*I_ERI_H3y2z_S_G3yz_S_vrr+WPX*I_ERI_H3y2z_S_G3yz_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G3yz_S_vrr = PAX*I_ERI_H2y3z_S_G3yz_S_vrr+WPX*I_ERI_H2y3z_S_G3yz_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G3yz_S_vrr = PAY*I_ERI_Hx4z_S_G3yz_S_vrr+WPY*I_ERI_Hx4z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Hx4z_S_F2yz_S_M1_vrr;
      Double I_ERI_Ix5z_S_G3yz_S_vrr = PAX*I_ERI_H5z_S_G3yz_S_vrr+WPX*I_ERI_H5z_S_G3yz_S_M1_vrr;
      Double I_ERI_I6y_S_G3yz_S_vrr = PAY*I_ERI_H5y_S_G3yz_S_vrr+WPY*I_ERI_H5y_S_G3yz_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G3yz_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_H5y_S_F2yz_S_M1_vrr;
      Double I_ERI_I5yz_S_G3yz_S_vrr = PAZ*I_ERI_H5y_S_G3yz_S_vrr+WPZ*I_ERI_H5y_S_G3yz_S_M1_vrr+oned2k*I_ERI_H5y_S_F3y_S_M1_vrr;
      Double I_ERI_I4y2z_S_G3yz_S_vrr = PAZ*I_ERI_H4yz_S_G3yz_S_vrr+WPZ*I_ERI_H4yz_S_G3yz_S_M1_vrr+oned2z*I_ERI_G4y_S_G3yz_S_vrr-rhod2zsq*I_ERI_G4y_S_G3yz_S_M1_vrr+oned2k*I_ERI_H4yz_S_F3y_S_M1_vrr;
      Double I_ERI_I3y3z_S_G3yz_S_vrr = PAZ*I_ERI_H3y2z_S_G3yz_S_vrr+WPZ*I_ERI_H3y2z_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G3yz_S_M1_vrr+oned2k*I_ERI_H3y2z_S_F3y_S_M1_vrr;
      Double I_ERI_I2y4z_S_G3yz_S_vrr = PAY*I_ERI_Hy4z_S_G3yz_S_vrr+WPY*I_ERI_Hy4z_S_G3yz_S_M1_vrr+oned2z*I_ERI_G4z_S_G3yz_S_vrr-rhod2zsq*I_ERI_G4z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Hy4z_S_F2yz_S_M1_vrr;
      Double I_ERI_Iy5z_S_G3yz_S_vrr = PAY*I_ERI_H5z_S_G3yz_S_vrr+WPY*I_ERI_H5z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_H5z_S_F2yz_S_M1_vrr;
      Double I_ERI_I6z_S_G3yz_S_vrr = PAZ*I_ERI_H5z_S_G3yz_S_vrr+WPZ*I_ERI_H5z_S_G3yz_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G3yz_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G3yz_S_M1_vrr+oned2k*I_ERI_H5z_S_F3y_S_M1_vrr;
      Double I_ERI_I6x_S_G2y2z_S_vrr = PAX*I_ERI_H5x_S_G2y2z_S_vrr+WPX*I_ERI_H5x_S_G2y2z_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G2y2z_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G2y2z_S_M1_vrr;
      Double I_ERI_I5xy_S_G2y2z_S_vrr = PAY*I_ERI_H5x_S_G2y2z_S_vrr+WPY*I_ERI_H5x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H5x_S_Fy2z_S_M1_vrr;
      Double I_ERI_I5xz_S_G2y2z_S_vrr = PAZ*I_ERI_H5x_S_G2y2z_S_vrr+WPZ*I_ERI_H5x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H5x_S_F2yz_S_M1_vrr;
      Double I_ERI_I4x2y_S_G2y2z_S_vrr = PAY*I_ERI_H4xy_S_G2y2z_S_vrr+WPY*I_ERI_H4xy_S_G2y2z_S_M1_vrr+oned2z*I_ERI_G4x_S_G2y2z_S_vrr-rhod2zsq*I_ERI_G4x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H4xy_S_Fy2z_S_M1_vrr;
      Double I_ERI_I4xyz_S_G2y2z_S_vrr = PAZ*I_ERI_H4xy_S_G2y2z_S_vrr+WPZ*I_ERI_H4xy_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H4xy_S_F2yz_S_M1_vrr;
      Double I_ERI_I4x2z_S_G2y2z_S_vrr = PAZ*I_ERI_H4xz_S_G2y2z_S_vrr+WPZ*I_ERI_H4xz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_G4x_S_G2y2z_S_vrr-rhod2zsq*I_ERI_G4x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H4xz_S_F2yz_S_M1_vrr;
      Double I_ERI_I3x3y_S_G2y2z_S_vrr = PAY*I_ERI_H3x2y_S_G2y2z_S_vrr+WPY*I_ERI_H3x2y_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H3x2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G2y2z_S_vrr = PAZ*I_ERI_H3x2y_S_G2y2z_S_vrr+WPZ*I_ERI_H3x2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H3x2y_S_F2yz_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G2y2z_S_vrr = PAY*I_ERI_H3x2z_S_G2y2z_S_vrr+WPY*I_ERI_H3x2z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H3x2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_I3x3z_S_G2y2z_S_vrr = PAZ*I_ERI_H3x2z_S_G2y2z_S_vrr+WPZ*I_ERI_H3x2z_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H3x2z_S_F2yz_S_M1_vrr;
      Double I_ERI_I2x4y_S_G2y2z_S_vrr = PAX*I_ERI_Hx4y_S_G2y2z_S_vrr+WPX*I_ERI_Hx4y_S_G2y2z_S_M1_vrr+oned2z*I_ERI_G4y_S_G2y2z_S_vrr-rhod2zsq*I_ERI_G4y_S_G2y2z_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G2y2z_S_vrr = PAZ*I_ERI_H2x3y_S_G2y2z_S_vrr+WPZ*I_ERI_H2x3y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H2x3y_S_F2yz_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G2y2z_S_vrr = PAZ*I_ERI_H2x2yz_S_G2y2z_S_vrr+WPZ*I_ERI_H2x2yz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G2y2z_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H2x2yz_S_F2yz_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G2y2z_S_vrr = PAY*I_ERI_H2x3z_S_G2y2z_S_vrr+WPY*I_ERI_H2x3z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H2x3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_I2x4z_S_G2y2z_S_vrr = PAX*I_ERI_Hx4z_S_G2y2z_S_vrr+WPX*I_ERI_Hx4z_S_G2y2z_S_M1_vrr+oned2z*I_ERI_G4z_S_G2y2z_S_vrr-rhod2zsq*I_ERI_G4z_S_G2y2z_S_M1_vrr;
      Double I_ERI_Ix5y_S_G2y2z_S_vrr = PAX*I_ERI_H5y_S_G2y2z_S_vrr+WPX*I_ERI_H5y_S_G2y2z_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G2y2z_S_vrr = PAZ*I_ERI_Hx4y_S_G2y2z_S_vrr+WPZ*I_ERI_Hx4y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Hx4y_S_F2yz_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G2y2z_S_vrr = PAX*I_ERI_H3y2z_S_G2y2z_S_vrr+WPX*I_ERI_H3y2z_S_G2y2z_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G2y2z_S_vrr = PAX*I_ERI_H2y3z_S_G2y2z_S_vrr+WPX*I_ERI_H2y3z_S_G2y2z_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G2y2z_S_vrr = PAY*I_ERI_Hx4z_S_G2y2z_S_vrr+WPY*I_ERI_Hx4z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Hx4z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Ix5z_S_G2y2z_S_vrr = PAX*I_ERI_H5z_S_G2y2z_S_vrr+WPX*I_ERI_H5z_S_G2y2z_S_M1_vrr;
      Double I_ERI_I6y_S_G2y2z_S_vrr = PAY*I_ERI_H5y_S_G2y2z_S_vrr+WPY*I_ERI_H5y_S_G2y2z_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G2y2z_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H5y_S_Fy2z_S_M1_vrr;
      Double I_ERI_I5yz_S_G2y2z_S_vrr = PAZ*I_ERI_H5y_S_G2y2z_S_vrr+WPZ*I_ERI_H5y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H5y_S_F2yz_S_M1_vrr;
      Double I_ERI_I4y2z_S_G2y2z_S_vrr = PAZ*I_ERI_H4yz_S_G2y2z_S_vrr+WPZ*I_ERI_H4yz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_G4y_S_G2y2z_S_vrr-rhod2zsq*I_ERI_G4y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H4yz_S_F2yz_S_M1_vrr;
      Double I_ERI_I3y3z_S_G2y2z_S_vrr = PAZ*I_ERI_H3y2z_S_G2y2z_S_vrr+WPZ*I_ERI_H3y2z_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H3y2z_S_F2yz_S_M1_vrr;
      Double I_ERI_I2y4z_S_G2y2z_S_vrr = PAY*I_ERI_Hy4z_S_G2y2z_S_vrr+WPY*I_ERI_Hy4z_S_G2y2z_S_M1_vrr+oned2z*I_ERI_G4z_S_G2y2z_S_vrr-rhod2zsq*I_ERI_G4z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Hy4z_S_Fy2z_S_M1_vrr;
      Double I_ERI_Iy5z_S_G2y2z_S_vrr = PAY*I_ERI_H5z_S_G2y2z_S_vrr+WPY*I_ERI_H5z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H5z_S_Fy2z_S_M1_vrr;
      Double I_ERI_I6z_S_G2y2z_S_vrr = PAZ*I_ERI_H5z_S_G2y2z_S_vrr+WPZ*I_ERI_H5z_S_G2y2z_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G2y2z_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_H5z_S_F2yz_S_M1_vrr;
      Double I_ERI_I6x_S_Gy3z_S_vrr = PAX*I_ERI_H5x_S_Gy3z_S_vrr+WPX*I_ERI_H5x_S_Gy3z_S_M1_vrr+5*oned2z*I_ERI_G4x_S_Gy3z_S_vrr-5*rhod2zsq*I_ERI_G4x_S_Gy3z_S_M1_vrr;
      Double I_ERI_I5xy_S_Gy3z_S_vrr = PAY*I_ERI_H5x_S_Gy3z_S_vrr+WPY*I_ERI_H5x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_H5x_S_F3z_S_M1_vrr;
      Double I_ERI_I5xz_S_Gy3z_S_vrr = PAZ*I_ERI_H5x_S_Gy3z_S_vrr+WPZ*I_ERI_H5x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H5x_S_Fy2z_S_M1_vrr;
      Double I_ERI_I4x2y_S_Gy3z_S_vrr = PAY*I_ERI_H4xy_S_Gy3z_S_vrr+WPY*I_ERI_H4xy_S_Gy3z_S_M1_vrr+oned2z*I_ERI_G4x_S_Gy3z_S_vrr-rhod2zsq*I_ERI_G4x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_H4xy_S_F3z_S_M1_vrr;
      Double I_ERI_I4xyz_S_Gy3z_S_vrr = PAZ*I_ERI_H4xy_S_Gy3z_S_vrr+WPZ*I_ERI_H4xy_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H4xy_S_Fy2z_S_M1_vrr;
      Double I_ERI_I4x2z_S_Gy3z_S_vrr = PAZ*I_ERI_H4xz_S_Gy3z_S_vrr+WPZ*I_ERI_H4xz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_G4x_S_Gy3z_S_vrr-rhod2zsq*I_ERI_G4x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H4xz_S_Fy2z_S_M1_vrr;
      Double I_ERI_I3x3y_S_Gy3z_S_vrr = PAY*I_ERI_H3x2y_S_Gy3z_S_vrr+WPY*I_ERI_H3x2y_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_Gy3z_S_M1_vrr+oned2k*I_ERI_H3x2y_S_F3z_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Gy3z_S_vrr = PAZ*I_ERI_H3x2y_S_Gy3z_S_vrr+WPZ*I_ERI_H3x2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H3x2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Gy3z_S_vrr = PAY*I_ERI_H3x2z_S_Gy3z_S_vrr+WPY*I_ERI_H3x2z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_H3x2z_S_F3z_S_M1_vrr;
      Double I_ERI_I3x3z_S_Gy3z_S_vrr = PAZ*I_ERI_H3x2z_S_Gy3z_S_vrr+WPZ*I_ERI_H3x2z_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H3x2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_I2x4y_S_Gy3z_S_vrr = PAX*I_ERI_Hx4y_S_Gy3z_S_vrr+WPX*I_ERI_Hx4y_S_Gy3z_S_M1_vrr+oned2z*I_ERI_G4y_S_Gy3z_S_vrr-rhod2zsq*I_ERI_G4y_S_Gy3z_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Gy3z_S_vrr = PAZ*I_ERI_H2x3y_S_Gy3z_S_vrr+WPZ*I_ERI_H2x3y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H2x3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Gy3z_S_vrr = PAZ*I_ERI_H2x2yz_S_Gy3z_S_vrr+WPZ*I_ERI_H2x2yz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_G2x2y_S_Gy3z_S_vrr-rhod2zsq*I_ERI_G2x2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H2x2yz_S_Fy2z_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Gy3z_S_vrr = PAY*I_ERI_H2x3z_S_Gy3z_S_vrr+WPY*I_ERI_H2x3z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_H2x3z_S_F3z_S_M1_vrr;
      Double I_ERI_I2x4z_S_Gy3z_S_vrr = PAX*I_ERI_Hx4z_S_Gy3z_S_vrr+WPX*I_ERI_Hx4z_S_Gy3z_S_M1_vrr+oned2z*I_ERI_G4z_S_Gy3z_S_vrr-rhod2zsq*I_ERI_G4z_S_Gy3z_S_M1_vrr;
      Double I_ERI_Ix5y_S_Gy3z_S_vrr = PAX*I_ERI_H5y_S_Gy3z_S_vrr+WPX*I_ERI_H5y_S_Gy3z_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Gy3z_S_vrr = PAZ*I_ERI_Hx4y_S_Gy3z_S_vrr+WPZ*I_ERI_Hx4y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Hx4y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Gy3z_S_vrr = PAX*I_ERI_H3y2z_S_Gy3z_S_vrr+WPX*I_ERI_H3y2z_S_Gy3z_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Gy3z_S_vrr = PAX*I_ERI_H2y3z_S_Gy3z_S_vrr+WPX*I_ERI_H2y3z_S_Gy3z_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Gy3z_S_vrr = PAY*I_ERI_Hx4z_S_Gy3z_S_vrr+WPY*I_ERI_Hx4z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Hx4z_S_F3z_S_M1_vrr;
      Double I_ERI_Ix5z_S_Gy3z_S_vrr = PAX*I_ERI_H5z_S_Gy3z_S_vrr+WPX*I_ERI_H5z_S_Gy3z_S_M1_vrr;
      Double I_ERI_I6y_S_Gy3z_S_vrr = PAY*I_ERI_H5y_S_Gy3z_S_vrr+WPY*I_ERI_H5y_S_Gy3z_S_M1_vrr+5*oned2z*I_ERI_G4y_S_Gy3z_S_vrr-5*rhod2zsq*I_ERI_G4y_S_Gy3z_S_M1_vrr+oned2k*I_ERI_H5y_S_F3z_S_M1_vrr;
      Double I_ERI_I5yz_S_Gy3z_S_vrr = PAZ*I_ERI_H5y_S_Gy3z_S_vrr+WPZ*I_ERI_H5y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H5y_S_Fy2z_S_M1_vrr;
      Double I_ERI_I4y2z_S_Gy3z_S_vrr = PAZ*I_ERI_H4yz_S_Gy3z_S_vrr+WPZ*I_ERI_H4yz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_G4y_S_Gy3z_S_vrr-rhod2zsq*I_ERI_G4y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H4yz_S_Fy2z_S_M1_vrr;
      Double I_ERI_I3y3z_S_Gy3z_S_vrr = PAZ*I_ERI_H3y2z_S_Gy3z_S_vrr+WPZ*I_ERI_H3y2z_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H3y2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_I2y4z_S_Gy3z_S_vrr = PAY*I_ERI_Hy4z_S_Gy3z_S_vrr+WPY*I_ERI_Hy4z_S_Gy3z_S_M1_vrr+oned2z*I_ERI_G4z_S_Gy3z_S_vrr-rhod2zsq*I_ERI_G4z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Hy4z_S_F3z_S_M1_vrr;
      Double I_ERI_Iy5z_S_Gy3z_S_vrr = PAY*I_ERI_H5z_S_Gy3z_S_vrr+WPY*I_ERI_H5z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_H5z_S_F3z_S_M1_vrr;
      Double I_ERI_I6z_S_Gy3z_S_vrr = PAZ*I_ERI_H5z_S_Gy3z_S_vrr+WPZ*I_ERI_H5z_S_Gy3z_S_M1_vrr+5*oned2z*I_ERI_G4z_S_Gy3z_S_vrr-5*rhod2zsq*I_ERI_G4z_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_H5z_S_Fy2z_S_M1_vrr;
      Double I_ERI_I6x_S_G4z_S_vrr = PAX*I_ERI_H5x_S_G4z_S_vrr+WPX*I_ERI_H5x_S_G4z_S_M1_vrr+5*oned2z*I_ERI_G4x_S_G4z_S_vrr-5*rhod2zsq*I_ERI_G4x_S_G4z_S_M1_vrr;
      Double I_ERI_I5xy_S_G4z_S_vrr = PAY*I_ERI_H5x_S_G4z_S_vrr+WPY*I_ERI_H5x_S_G4z_S_M1_vrr;
      Double I_ERI_I5xz_S_G4z_S_vrr = PAZ*I_ERI_H5x_S_G4z_S_vrr+WPZ*I_ERI_H5x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H5x_S_F3z_S_M1_vrr;
      Double I_ERI_I4x2y_S_G4z_S_vrr = PAY*I_ERI_H4xy_S_G4z_S_vrr+WPY*I_ERI_H4xy_S_G4z_S_M1_vrr+oned2z*I_ERI_G4x_S_G4z_S_vrr-rhod2zsq*I_ERI_G4x_S_G4z_S_M1_vrr;
      Double I_ERI_I4xyz_S_G4z_S_vrr = PAZ*I_ERI_H4xy_S_G4z_S_vrr+WPZ*I_ERI_H4xy_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H4xy_S_F3z_S_M1_vrr;
      Double I_ERI_I4x2z_S_G4z_S_vrr = PAZ*I_ERI_H4xz_S_G4z_S_vrr+WPZ*I_ERI_H4xz_S_G4z_S_M1_vrr+oned2z*I_ERI_G4x_S_G4z_S_vrr-rhod2zsq*I_ERI_G4x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H4xz_S_F3z_S_M1_vrr;
      Double I_ERI_I3x3y_S_G4z_S_vrr = PAY*I_ERI_H3x2y_S_G4z_S_vrr+WPY*I_ERI_H3x2y_S_G4z_S_M1_vrr+2*oned2z*I_ERI_G3xy_S_G4z_S_vrr-2*rhod2zsq*I_ERI_G3xy_S_G4z_S_M1_vrr;
      Double I_ERI_I3x2yz_S_G4z_S_vrr = PAZ*I_ERI_H3x2y_S_G4z_S_vrr+WPZ*I_ERI_H3x2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H3x2y_S_F3z_S_M1_vrr;
      Double I_ERI_I3xy2z_S_G4z_S_vrr = PAY*I_ERI_H3x2z_S_G4z_S_vrr+WPY*I_ERI_H3x2z_S_G4z_S_M1_vrr;
      Double I_ERI_I3x3z_S_G4z_S_vrr = PAZ*I_ERI_H3x2z_S_G4z_S_vrr+WPZ*I_ERI_H3x2z_S_G4z_S_M1_vrr+2*oned2z*I_ERI_G3xz_S_G4z_S_vrr-2*rhod2zsq*I_ERI_G3xz_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H3x2z_S_F3z_S_M1_vrr;
      Double I_ERI_I2x4y_S_G4z_S_vrr = PAX*I_ERI_Hx4y_S_G4z_S_vrr+WPX*I_ERI_Hx4y_S_G4z_S_M1_vrr+oned2z*I_ERI_G4y_S_G4z_S_vrr-rhod2zsq*I_ERI_G4y_S_G4z_S_M1_vrr;
      Double I_ERI_I2x3yz_S_G4z_S_vrr = PAZ*I_ERI_H2x3y_S_G4z_S_vrr+WPZ*I_ERI_H2x3y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H2x3y_S_F3z_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_G4z_S_vrr = PAZ*I_ERI_H2x2yz_S_G4z_S_vrr+WPZ*I_ERI_H2x2yz_S_G4z_S_M1_vrr+oned2z*I_ERI_G2x2y_S_G4z_S_vrr-rhod2zsq*I_ERI_G2x2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H2x2yz_S_F3z_S_M1_vrr;
      Double I_ERI_I2xy3z_S_G4z_S_vrr = PAY*I_ERI_H2x3z_S_G4z_S_vrr+WPY*I_ERI_H2x3z_S_G4z_S_M1_vrr;
      Double I_ERI_I2x4z_S_G4z_S_vrr = PAX*I_ERI_Hx4z_S_G4z_S_vrr+WPX*I_ERI_Hx4z_S_G4z_S_M1_vrr+oned2z*I_ERI_G4z_S_G4z_S_vrr-rhod2zsq*I_ERI_G4z_S_G4z_S_M1_vrr;
      Double I_ERI_Ix5y_S_G4z_S_vrr = PAX*I_ERI_H5y_S_G4z_S_vrr+WPX*I_ERI_H5y_S_G4z_S_M1_vrr;
      Double I_ERI_Ix4yz_S_G4z_S_vrr = PAZ*I_ERI_Hx4y_S_G4z_S_vrr+WPZ*I_ERI_Hx4y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Hx4y_S_F3z_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_G4z_S_vrr = PAX*I_ERI_H3y2z_S_G4z_S_vrr+WPX*I_ERI_H3y2z_S_G4z_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_G4z_S_vrr = PAX*I_ERI_H2y3z_S_G4z_S_vrr+WPX*I_ERI_H2y3z_S_G4z_S_M1_vrr;
      Double I_ERI_Ixy4z_S_G4z_S_vrr = PAY*I_ERI_Hx4z_S_G4z_S_vrr+WPY*I_ERI_Hx4z_S_G4z_S_M1_vrr;
      Double I_ERI_Ix5z_S_G4z_S_vrr = PAX*I_ERI_H5z_S_G4z_S_vrr+WPX*I_ERI_H5z_S_G4z_S_M1_vrr;
      Double I_ERI_I6y_S_G4z_S_vrr = PAY*I_ERI_H5y_S_G4z_S_vrr+WPY*I_ERI_H5y_S_G4z_S_M1_vrr+5*oned2z*I_ERI_G4y_S_G4z_S_vrr-5*rhod2zsq*I_ERI_G4y_S_G4z_S_M1_vrr;
      Double I_ERI_I5yz_S_G4z_S_vrr = PAZ*I_ERI_H5y_S_G4z_S_vrr+WPZ*I_ERI_H5y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H5y_S_F3z_S_M1_vrr;
      Double I_ERI_I4y2z_S_G4z_S_vrr = PAZ*I_ERI_H4yz_S_G4z_S_vrr+WPZ*I_ERI_H4yz_S_G4z_S_M1_vrr+oned2z*I_ERI_G4y_S_G4z_S_vrr-rhod2zsq*I_ERI_G4y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H4yz_S_F3z_S_M1_vrr;
      Double I_ERI_I3y3z_S_G4z_S_vrr = PAZ*I_ERI_H3y2z_S_G4z_S_vrr+WPZ*I_ERI_H3y2z_S_G4z_S_M1_vrr+2*oned2z*I_ERI_G3yz_S_G4z_S_vrr-2*rhod2zsq*I_ERI_G3yz_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H3y2z_S_F3z_S_M1_vrr;
      Double I_ERI_I2y4z_S_G4z_S_vrr = PAY*I_ERI_Hy4z_S_G4z_S_vrr+WPY*I_ERI_Hy4z_S_G4z_S_M1_vrr+oned2z*I_ERI_G4z_S_G4z_S_vrr-rhod2zsq*I_ERI_G4z_S_G4z_S_M1_vrr;
      Double I_ERI_Iy5z_S_G4z_S_vrr = PAY*I_ERI_H5z_S_G4z_S_vrr+WPY*I_ERI_H5z_S_G4z_S_M1_vrr;
      Double I_ERI_I6z_S_G4z_S_vrr = PAZ*I_ERI_H5z_S_G4z_S_vrr+WPZ*I_ERI_H5z_S_G4z_S_M1_vrr+5*oned2z*I_ERI_G4z_S_G4z_S_vrr-5*rhod2zsq*I_ERI_G4z_S_G4z_S_M1_vrr+4*oned2k*I_ERI_H5z_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_G_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_I6x_S_G4x_S += I_ERI_I6x_S_G4x_S_vrr;
      I_ERI_I5xy_S_G4x_S += I_ERI_I5xy_S_G4x_S_vrr;
      I_ERI_I5xz_S_G4x_S += I_ERI_I5xz_S_G4x_S_vrr;
      I_ERI_I4x2y_S_G4x_S += I_ERI_I4x2y_S_G4x_S_vrr;
      I_ERI_I4xyz_S_G4x_S += I_ERI_I4xyz_S_G4x_S_vrr;
      I_ERI_I4x2z_S_G4x_S += I_ERI_I4x2z_S_G4x_S_vrr;
      I_ERI_I3x3y_S_G4x_S += I_ERI_I3x3y_S_G4x_S_vrr;
      I_ERI_I3x2yz_S_G4x_S += I_ERI_I3x2yz_S_G4x_S_vrr;
      I_ERI_I3xy2z_S_G4x_S += I_ERI_I3xy2z_S_G4x_S_vrr;
      I_ERI_I3x3z_S_G4x_S += I_ERI_I3x3z_S_G4x_S_vrr;
      I_ERI_I2x4y_S_G4x_S += I_ERI_I2x4y_S_G4x_S_vrr;
      I_ERI_I2x3yz_S_G4x_S += I_ERI_I2x3yz_S_G4x_S_vrr;
      I_ERI_I2x2y2z_S_G4x_S += I_ERI_I2x2y2z_S_G4x_S_vrr;
      I_ERI_I2xy3z_S_G4x_S += I_ERI_I2xy3z_S_G4x_S_vrr;
      I_ERI_I2x4z_S_G4x_S += I_ERI_I2x4z_S_G4x_S_vrr;
      I_ERI_Ix5y_S_G4x_S += I_ERI_Ix5y_S_G4x_S_vrr;
      I_ERI_Ix4yz_S_G4x_S += I_ERI_Ix4yz_S_G4x_S_vrr;
      I_ERI_Ix3y2z_S_G4x_S += I_ERI_Ix3y2z_S_G4x_S_vrr;
      I_ERI_Ix2y3z_S_G4x_S += I_ERI_Ix2y3z_S_G4x_S_vrr;
      I_ERI_Ixy4z_S_G4x_S += I_ERI_Ixy4z_S_G4x_S_vrr;
      I_ERI_Ix5z_S_G4x_S += I_ERI_Ix5z_S_G4x_S_vrr;
      I_ERI_I6y_S_G4x_S += I_ERI_I6y_S_G4x_S_vrr;
      I_ERI_I5yz_S_G4x_S += I_ERI_I5yz_S_G4x_S_vrr;
      I_ERI_I4y2z_S_G4x_S += I_ERI_I4y2z_S_G4x_S_vrr;
      I_ERI_I3y3z_S_G4x_S += I_ERI_I3y3z_S_G4x_S_vrr;
      I_ERI_I2y4z_S_G4x_S += I_ERI_I2y4z_S_G4x_S_vrr;
      I_ERI_Iy5z_S_G4x_S += I_ERI_Iy5z_S_G4x_S_vrr;
      I_ERI_I6z_S_G4x_S += I_ERI_I6z_S_G4x_S_vrr;
      I_ERI_I6x_S_G3xy_S += I_ERI_I6x_S_G3xy_S_vrr;
      I_ERI_I5xy_S_G3xy_S += I_ERI_I5xy_S_G3xy_S_vrr;
      I_ERI_I5xz_S_G3xy_S += I_ERI_I5xz_S_G3xy_S_vrr;
      I_ERI_I4x2y_S_G3xy_S += I_ERI_I4x2y_S_G3xy_S_vrr;
      I_ERI_I4xyz_S_G3xy_S += I_ERI_I4xyz_S_G3xy_S_vrr;
      I_ERI_I4x2z_S_G3xy_S += I_ERI_I4x2z_S_G3xy_S_vrr;
      I_ERI_I3x3y_S_G3xy_S += I_ERI_I3x3y_S_G3xy_S_vrr;
      I_ERI_I3x2yz_S_G3xy_S += I_ERI_I3x2yz_S_G3xy_S_vrr;
      I_ERI_I3xy2z_S_G3xy_S += I_ERI_I3xy2z_S_G3xy_S_vrr;
      I_ERI_I3x3z_S_G3xy_S += I_ERI_I3x3z_S_G3xy_S_vrr;
      I_ERI_I2x4y_S_G3xy_S += I_ERI_I2x4y_S_G3xy_S_vrr;
      I_ERI_I2x3yz_S_G3xy_S += I_ERI_I2x3yz_S_G3xy_S_vrr;
      I_ERI_I2x2y2z_S_G3xy_S += I_ERI_I2x2y2z_S_G3xy_S_vrr;
      I_ERI_I2xy3z_S_G3xy_S += I_ERI_I2xy3z_S_G3xy_S_vrr;
      I_ERI_I2x4z_S_G3xy_S += I_ERI_I2x4z_S_G3xy_S_vrr;
      I_ERI_Ix5y_S_G3xy_S += I_ERI_Ix5y_S_G3xy_S_vrr;
      I_ERI_Ix4yz_S_G3xy_S += I_ERI_Ix4yz_S_G3xy_S_vrr;
      I_ERI_Ix3y2z_S_G3xy_S += I_ERI_Ix3y2z_S_G3xy_S_vrr;
      I_ERI_Ix2y3z_S_G3xy_S += I_ERI_Ix2y3z_S_G3xy_S_vrr;
      I_ERI_Ixy4z_S_G3xy_S += I_ERI_Ixy4z_S_G3xy_S_vrr;
      I_ERI_Ix5z_S_G3xy_S += I_ERI_Ix5z_S_G3xy_S_vrr;
      I_ERI_I6y_S_G3xy_S += I_ERI_I6y_S_G3xy_S_vrr;
      I_ERI_I5yz_S_G3xy_S += I_ERI_I5yz_S_G3xy_S_vrr;
      I_ERI_I4y2z_S_G3xy_S += I_ERI_I4y2z_S_G3xy_S_vrr;
      I_ERI_I3y3z_S_G3xy_S += I_ERI_I3y3z_S_G3xy_S_vrr;
      I_ERI_I2y4z_S_G3xy_S += I_ERI_I2y4z_S_G3xy_S_vrr;
      I_ERI_Iy5z_S_G3xy_S += I_ERI_Iy5z_S_G3xy_S_vrr;
      I_ERI_I6z_S_G3xy_S += I_ERI_I6z_S_G3xy_S_vrr;
      I_ERI_I6x_S_G3xz_S += I_ERI_I6x_S_G3xz_S_vrr;
      I_ERI_I5xy_S_G3xz_S += I_ERI_I5xy_S_G3xz_S_vrr;
      I_ERI_I5xz_S_G3xz_S += I_ERI_I5xz_S_G3xz_S_vrr;
      I_ERI_I4x2y_S_G3xz_S += I_ERI_I4x2y_S_G3xz_S_vrr;
      I_ERI_I4xyz_S_G3xz_S += I_ERI_I4xyz_S_G3xz_S_vrr;
      I_ERI_I4x2z_S_G3xz_S += I_ERI_I4x2z_S_G3xz_S_vrr;
      I_ERI_I3x3y_S_G3xz_S += I_ERI_I3x3y_S_G3xz_S_vrr;
      I_ERI_I3x2yz_S_G3xz_S += I_ERI_I3x2yz_S_G3xz_S_vrr;
      I_ERI_I3xy2z_S_G3xz_S += I_ERI_I3xy2z_S_G3xz_S_vrr;
      I_ERI_I3x3z_S_G3xz_S += I_ERI_I3x3z_S_G3xz_S_vrr;
      I_ERI_I2x4y_S_G3xz_S += I_ERI_I2x4y_S_G3xz_S_vrr;
      I_ERI_I2x3yz_S_G3xz_S += I_ERI_I2x3yz_S_G3xz_S_vrr;
      I_ERI_I2x2y2z_S_G3xz_S += I_ERI_I2x2y2z_S_G3xz_S_vrr;
      I_ERI_I2xy3z_S_G3xz_S += I_ERI_I2xy3z_S_G3xz_S_vrr;
      I_ERI_I2x4z_S_G3xz_S += I_ERI_I2x4z_S_G3xz_S_vrr;
      I_ERI_Ix5y_S_G3xz_S += I_ERI_Ix5y_S_G3xz_S_vrr;
      I_ERI_Ix4yz_S_G3xz_S += I_ERI_Ix4yz_S_G3xz_S_vrr;
      I_ERI_Ix3y2z_S_G3xz_S += I_ERI_Ix3y2z_S_G3xz_S_vrr;
      I_ERI_Ix2y3z_S_G3xz_S += I_ERI_Ix2y3z_S_G3xz_S_vrr;
      I_ERI_Ixy4z_S_G3xz_S += I_ERI_Ixy4z_S_G3xz_S_vrr;
      I_ERI_Ix5z_S_G3xz_S += I_ERI_Ix5z_S_G3xz_S_vrr;
      I_ERI_I6y_S_G3xz_S += I_ERI_I6y_S_G3xz_S_vrr;
      I_ERI_I5yz_S_G3xz_S += I_ERI_I5yz_S_G3xz_S_vrr;
      I_ERI_I4y2z_S_G3xz_S += I_ERI_I4y2z_S_G3xz_S_vrr;
      I_ERI_I3y3z_S_G3xz_S += I_ERI_I3y3z_S_G3xz_S_vrr;
      I_ERI_I2y4z_S_G3xz_S += I_ERI_I2y4z_S_G3xz_S_vrr;
      I_ERI_Iy5z_S_G3xz_S += I_ERI_Iy5z_S_G3xz_S_vrr;
      I_ERI_I6z_S_G3xz_S += I_ERI_I6z_S_G3xz_S_vrr;
      I_ERI_I6x_S_G2x2y_S += I_ERI_I6x_S_G2x2y_S_vrr;
      I_ERI_I5xy_S_G2x2y_S += I_ERI_I5xy_S_G2x2y_S_vrr;
      I_ERI_I5xz_S_G2x2y_S += I_ERI_I5xz_S_G2x2y_S_vrr;
      I_ERI_I4x2y_S_G2x2y_S += I_ERI_I4x2y_S_G2x2y_S_vrr;
      I_ERI_I4xyz_S_G2x2y_S += I_ERI_I4xyz_S_G2x2y_S_vrr;
      I_ERI_I4x2z_S_G2x2y_S += I_ERI_I4x2z_S_G2x2y_S_vrr;
      I_ERI_I3x3y_S_G2x2y_S += I_ERI_I3x3y_S_G2x2y_S_vrr;
      I_ERI_I3x2yz_S_G2x2y_S += I_ERI_I3x2yz_S_G2x2y_S_vrr;
      I_ERI_I3xy2z_S_G2x2y_S += I_ERI_I3xy2z_S_G2x2y_S_vrr;
      I_ERI_I3x3z_S_G2x2y_S += I_ERI_I3x3z_S_G2x2y_S_vrr;
      I_ERI_I2x4y_S_G2x2y_S += I_ERI_I2x4y_S_G2x2y_S_vrr;
      I_ERI_I2x3yz_S_G2x2y_S += I_ERI_I2x3yz_S_G2x2y_S_vrr;
      I_ERI_I2x2y2z_S_G2x2y_S += I_ERI_I2x2y2z_S_G2x2y_S_vrr;
      I_ERI_I2xy3z_S_G2x2y_S += I_ERI_I2xy3z_S_G2x2y_S_vrr;
      I_ERI_I2x4z_S_G2x2y_S += I_ERI_I2x4z_S_G2x2y_S_vrr;
      I_ERI_Ix5y_S_G2x2y_S += I_ERI_Ix5y_S_G2x2y_S_vrr;
      I_ERI_Ix4yz_S_G2x2y_S += I_ERI_Ix4yz_S_G2x2y_S_vrr;
      I_ERI_Ix3y2z_S_G2x2y_S += I_ERI_Ix3y2z_S_G2x2y_S_vrr;
      I_ERI_Ix2y3z_S_G2x2y_S += I_ERI_Ix2y3z_S_G2x2y_S_vrr;
      I_ERI_Ixy4z_S_G2x2y_S += I_ERI_Ixy4z_S_G2x2y_S_vrr;
      I_ERI_Ix5z_S_G2x2y_S += I_ERI_Ix5z_S_G2x2y_S_vrr;
      I_ERI_I6y_S_G2x2y_S += I_ERI_I6y_S_G2x2y_S_vrr;
      I_ERI_I5yz_S_G2x2y_S += I_ERI_I5yz_S_G2x2y_S_vrr;
      I_ERI_I4y2z_S_G2x2y_S += I_ERI_I4y2z_S_G2x2y_S_vrr;
      I_ERI_I3y3z_S_G2x2y_S += I_ERI_I3y3z_S_G2x2y_S_vrr;
      I_ERI_I2y4z_S_G2x2y_S += I_ERI_I2y4z_S_G2x2y_S_vrr;
      I_ERI_Iy5z_S_G2x2y_S += I_ERI_Iy5z_S_G2x2y_S_vrr;
      I_ERI_I6z_S_G2x2y_S += I_ERI_I6z_S_G2x2y_S_vrr;
      I_ERI_I6x_S_G2xyz_S += I_ERI_I6x_S_G2xyz_S_vrr;
      I_ERI_I5xy_S_G2xyz_S += I_ERI_I5xy_S_G2xyz_S_vrr;
      I_ERI_I5xz_S_G2xyz_S += I_ERI_I5xz_S_G2xyz_S_vrr;
      I_ERI_I4x2y_S_G2xyz_S += I_ERI_I4x2y_S_G2xyz_S_vrr;
      I_ERI_I4xyz_S_G2xyz_S += I_ERI_I4xyz_S_G2xyz_S_vrr;
      I_ERI_I4x2z_S_G2xyz_S += I_ERI_I4x2z_S_G2xyz_S_vrr;
      I_ERI_I3x3y_S_G2xyz_S += I_ERI_I3x3y_S_G2xyz_S_vrr;
      I_ERI_I3x2yz_S_G2xyz_S += I_ERI_I3x2yz_S_G2xyz_S_vrr;
      I_ERI_I3xy2z_S_G2xyz_S += I_ERI_I3xy2z_S_G2xyz_S_vrr;
      I_ERI_I3x3z_S_G2xyz_S += I_ERI_I3x3z_S_G2xyz_S_vrr;
      I_ERI_I2x4y_S_G2xyz_S += I_ERI_I2x4y_S_G2xyz_S_vrr;
      I_ERI_I2x3yz_S_G2xyz_S += I_ERI_I2x3yz_S_G2xyz_S_vrr;
      I_ERI_I2x2y2z_S_G2xyz_S += I_ERI_I2x2y2z_S_G2xyz_S_vrr;
      I_ERI_I2xy3z_S_G2xyz_S += I_ERI_I2xy3z_S_G2xyz_S_vrr;
      I_ERI_I2x4z_S_G2xyz_S += I_ERI_I2x4z_S_G2xyz_S_vrr;
      I_ERI_Ix5y_S_G2xyz_S += I_ERI_Ix5y_S_G2xyz_S_vrr;
      I_ERI_Ix4yz_S_G2xyz_S += I_ERI_Ix4yz_S_G2xyz_S_vrr;
      I_ERI_Ix3y2z_S_G2xyz_S += I_ERI_Ix3y2z_S_G2xyz_S_vrr;
      I_ERI_Ix2y3z_S_G2xyz_S += I_ERI_Ix2y3z_S_G2xyz_S_vrr;
      I_ERI_Ixy4z_S_G2xyz_S += I_ERI_Ixy4z_S_G2xyz_S_vrr;
      I_ERI_Ix5z_S_G2xyz_S += I_ERI_Ix5z_S_G2xyz_S_vrr;
      I_ERI_I6y_S_G2xyz_S += I_ERI_I6y_S_G2xyz_S_vrr;
      I_ERI_I5yz_S_G2xyz_S += I_ERI_I5yz_S_G2xyz_S_vrr;
      I_ERI_I4y2z_S_G2xyz_S += I_ERI_I4y2z_S_G2xyz_S_vrr;
      I_ERI_I3y3z_S_G2xyz_S += I_ERI_I3y3z_S_G2xyz_S_vrr;
      I_ERI_I2y4z_S_G2xyz_S += I_ERI_I2y4z_S_G2xyz_S_vrr;
      I_ERI_Iy5z_S_G2xyz_S += I_ERI_Iy5z_S_G2xyz_S_vrr;
      I_ERI_I6z_S_G2xyz_S += I_ERI_I6z_S_G2xyz_S_vrr;
      I_ERI_I6x_S_G2x2z_S += I_ERI_I6x_S_G2x2z_S_vrr;
      I_ERI_I5xy_S_G2x2z_S += I_ERI_I5xy_S_G2x2z_S_vrr;
      I_ERI_I5xz_S_G2x2z_S += I_ERI_I5xz_S_G2x2z_S_vrr;
      I_ERI_I4x2y_S_G2x2z_S += I_ERI_I4x2y_S_G2x2z_S_vrr;
      I_ERI_I4xyz_S_G2x2z_S += I_ERI_I4xyz_S_G2x2z_S_vrr;
      I_ERI_I4x2z_S_G2x2z_S += I_ERI_I4x2z_S_G2x2z_S_vrr;
      I_ERI_I3x3y_S_G2x2z_S += I_ERI_I3x3y_S_G2x2z_S_vrr;
      I_ERI_I3x2yz_S_G2x2z_S += I_ERI_I3x2yz_S_G2x2z_S_vrr;
      I_ERI_I3xy2z_S_G2x2z_S += I_ERI_I3xy2z_S_G2x2z_S_vrr;
      I_ERI_I3x3z_S_G2x2z_S += I_ERI_I3x3z_S_G2x2z_S_vrr;
      I_ERI_I2x4y_S_G2x2z_S += I_ERI_I2x4y_S_G2x2z_S_vrr;
      I_ERI_I2x3yz_S_G2x2z_S += I_ERI_I2x3yz_S_G2x2z_S_vrr;
      I_ERI_I2x2y2z_S_G2x2z_S += I_ERI_I2x2y2z_S_G2x2z_S_vrr;
      I_ERI_I2xy3z_S_G2x2z_S += I_ERI_I2xy3z_S_G2x2z_S_vrr;
      I_ERI_I2x4z_S_G2x2z_S += I_ERI_I2x4z_S_G2x2z_S_vrr;
      I_ERI_Ix5y_S_G2x2z_S += I_ERI_Ix5y_S_G2x2z_S_vrr;
      I_ERI_Ix4yz_S_G2x2z_S += I_ERI_Ix4yz_S_G2x2z_S_vrr;
      I_ERI_Ix3y2z_S_G2x2z_S += I_ERI_Ix3y2z_S_G2x2z_S_vrr;
      I_ERI_Ix2y3z_S_G2x2z_S += I_ERI_Ix2y3z_S_G2x2z_S_vrr;
      I_ERI_Ixy4z_S_G2x2z_S += I_ERI_Ixy4z_S_G2x2z_S_vrr;
      I_ERI_Ix5z_S_G2x2z_S += I_ERI_Ix5z_S_G2x2z_S_vrr;
      I_ERI_I6y_S_G2x2z_S += I_ERI_I6y_S_G2x2z_S_vrr;
      I_ERI_I5yz_S_G2x2z_S += I_ERI_I5yz_S_G2x2z_S_vrr;
      I_ERI_I4y2z_S_G2x2z_S += I_ERI_I4y2z_S_G2x2z_S_vrr;
      I_ERI_I3y3z_S_G2x2z_S += I_ERI_I3y3z_S_G2x2z_S_vrr;
      I_ERI_I2y4z_S_G2x2z_S += I_ERI_I2y4z_S_G2x2z_S_vrr;
      I_ERI_Iy5z_S_G2x2z_S += I_ERI_Iy5z_S_G2x2z_S_vrr;
      I_ERI_I6z_S_G2x2z_S += I_ERI_I6z_S_G2x2z_S_vrr;
      I_ERI_I6x_S_Gx3y_S += I_ERI_I6x_S_Gx3y_S_vrr;
      I_ERI_I5xy_S_Gx3y_S += I_ERI_I5xy_S_Gx3y_S_vrr;
      I_ERI_I5xz_S_Gx3y_S += I_ERI_I5xz_S_Gx3y_S_vrr;
      I_ERI_I4x2y_S_Gx3y_S += I_ERI_I4x2y_S_Gx3y_S_vrr;
      I_ERI_I4xyz_S_Gx3y_S += I_ERI_I4xyz_S_Gx3y_S_vrr;
      I_ERI_I4x2z_S_Gx3y_S += I_ERI_I4x2z_S_Gx3y_S_vrr;
      I_ERI_I3x3y_S_Gx3y_S += I_ERI_I3x3y_S_Gx3y_S_vrr;
      I_ERI_I3x2yz_S_Gx3y_S += I_ERI_I3x2yz_S_Gx3y_S_vrr;
      I_ERI_I3xy2z_S_Gx3y_S += I_ERI_I3xy2z_S_Gx3y_S_vrr;
      I_ERI_I3x3z_S_Gx3y_S += I_ERI_I3x3z_S_Gx3y_S_vrr;
      I_ERI_I2x4y_S_Gx3y_S += I_ERI_I2x4y_S_Gx3y_S_vrr;
      I_ERI_I2x3yz_S_Gx3y_S += I_ERI_I2x3yz_S_Gx3y_S_vrr;
      I_ERI_I2x2y2z_S_Gx3y_S += I_ERI_I2x2y2z_S_Gx3y_S_vrr;
      I_ERI_I2xy3z_S_Gx3y_S += I_ERI_I2xy3z_S_Gx3y_S_vrr;
      I_ERI_I2x4z_S_Gx3y_S += I_ERI_I2x4z_S_Gx3y_S_vrr;
      I_ERI_Ix5y_S_Gx3y_S += I_ERI_Ix5y_S_Gx3y_S_vrr;
      I_ERI_Ix4yz_S_Gx3y_S += I_ERI_Ix4yz_S_Gx3y_S_vrr;
      I_ERI_Ix3y2z_S_Gx3y_S += I_ERI_Ix3y2z_S_Gx3y_S_vrr;
      I_ERI_Ix2y3z_S_Gx3y_S += I_ERI_Ix2y3z_S_Gx3y_S_vrr;
      I_ERI_Ixy4z_S_Gx3y_S += I_ERI_Ixy4z_S_Gx3y_S_vrr;
      I_ERI_Ix5z_S_Gx3y_S += I_ERI_Ix5z_S_Gx3y_S_vrr;
      I_ERI_I6y_S_Gx3y_S += I_ERI_I6y_S_Gx3y_S_vrr;
      I_ERI_I5yz_S_Gx3y_S += I_ERI_I5yz_S_Gx3y_S_vrr;
      I_ERI_I4y2z_S_Gx3y_S += I_ERI_I4y2z_S_Gx3y_S_vrr;
      I_ERI_I3y3z_S_Gx3y_S += I_ERI_I3y3z_S_Gx3y_S_vrr;
      I_ERI_I2y4z_S_Gx3y_S += I_ERI_I2y4z_S_Gx3y_S_vrr;
      I_ERI_Iy5z_S_Gx3y_S += I_ERI_Iy5z_S_Gx3y_S_vrr;
      I_ERI_I6z_S_Gx3y_S += I_ERI_I6z_S_Gx3y_S_vrr;
      I_ERI_I6x_S_Gx2yz_S += I_ERI_I6x_S_Gx2yz_S_vrr;
      I_ERI_I5xy_S_Gx2yz_S += I_ERI_I5xy_S_Gx2yz_S_vrr;
      I_ERI_I5xz_S_Gx2yz_S += I_ERI_I5xz_S_Gx2yz_S_vrr;
      I_ERI_I4x2y_S_Gx2yz_S += I_ERI_I4x2y_S_Gx2yz_S_vrr;
      I_ERI_I4xyz_S_Gx2yz_S += I_ERI_I4xyz_S_Gx2yz_S_vrr;
      I_ERI_I4x2z_S_Gx2yz_S += I_ERI_I4x2z_S_Gx2yz_S_vrr;
      I_ERI_I3x3y_S_Gx2yz_S += I_ERI_I3x3y_S_Gx2yz_S_vrr;
      I_ERI_I3x2yz_S_Gx2yz_S += I_ERI_I3x2yz_S_Gx2yz_S_vrr;
      I_ERI_I3xy2z_S_Gx2yz_S += I_ERI_I3xy2z_S_Gx2yz_S_vrr;
      I_ERI_I3x3z_S_Gx2yz_S += I_ERI_I3x3z_S_Gx2yz_S_vrr;
      I_ERI_I2x4y_S_Gx2yz_S += I_ERI_I2x4y_S_Gx2yz_S_vrr;
      I_ERI_I2x3yz_S_Gx2yz_S += I_ERI_I2x3yz_S_Gx2yz_S_vrr;
      I_ERI_I2x2y2z_S_Gx2yz_S += I_ERI_I2x2y2z_S_Gx2yz_S_vrr;
      I_ERI_I2xy3z_S_Gx2yz_S += I_ERI_I2xy3z_S_Gx2yz_S_vrr;
      I_ERI_I2x4z_S_Gx2yz_S += I_ERI_I2x4z_S_Gx2yz_S_vrr;
      I_ERI_Ix5y_S_Gx2yz_S += I_ERI_Ix5y_S_Gx2yz_S_vrr;
      I_ERI_Ix4yz_S_Gx2yz_S += I_ERI_Ix4yz_S_Gx2yz_S_vrr;
      I_ERI_Ix3y2z_S_Gx2yz_S += I_ERI_Ix3y2z_S_Gx2yz_S_vrr;
      I_ERI_Ix2y3z_S_Gx2yz_S += I_ERI_Ix2y3z_S_Gx2yz_S_vrr;
      I_ERI_Ixy4z_S_Gx2yz_S += I_ERI_Ixy4z_S_Gx2yz_S_vrr;
      I_ERI_Ix5z_S_Gx2yz_S += I_ERI_Ix5z_S_Gx2yz_S_vrr;
      I_ERI_I6y_S_Gx2yz_S += I_ERI_I6y_S_Gx2yz_S_vrr;
      I_ERI_I5yz_S_Gx2yz_S += I_ERI_I5yz_S_Gx2yz_S_vrr;
      I_ERI_I4y2z_S_Gx2yz_S += I_ERI_I4y2z_S_Gx2yz_S_vrr;
      I_ERI_I3y3z_S_Gx2yz_S += I_ERI_I3y3z_S_Gx2yz_S_vrr;
      I_ERI_I2y4z_S_Gx2yz_S += I_ERI_I2y4z_S_Gx2yz_S_vrr;
      I_ERI_Iy5z_S_Gx2yz_S += I_ERI_Iy5z_S_Gx2yz_S_vrr;
      I_ERI_I6z_S_Gx2yz_S += I_ERI_I6z_S_Gx2yz_S_vrr;
      I_ERI_I6x_S_Gxy2z_S += I_ERI_I6x_S_Gxy2z_S_vrr;
      I_ERI_I5xy_S_Gxy2z_S += I_ERI_I5xy_S_Gxy2z_S_vrr;
      I_ERI_I5xz_S_Gxy2z_S += I_ERI_I5xz_S_Gxy2z_S_vrr;
      I_ERI_I4x2y_S_Gxy2z_S += I_ERI_I4x2y_S_Gxy2z_S_vrr;
      I_ERI_I4xyz_S_Gxy2z_S += I_ERI_I4xyz_S_Gxy2z_S_vrr;
      I_ERI_I4x2z_S_Gxy2z_S += I_ERI_I4x2z_S_Gxy2z_S_vrr;
      I_ERI_I3x3y_S_Gxy2z_S += I_ERI_I3x3y_S_Gxy2z_S_vrr;
      I_ERI_I3x2yz_S_Gxy2z_S += I_ERI_I3x2yz_S_Gxy2z_S_vrr;
      I_ERI_I3xy2z_S_Gxy2z_S += I_ERI_I3xy2z_S_Gxy2z_S_vrr;
      I_ERI_I3x3z_S_Gxy2z_S += I_ERI_I3x3z_S_Gxy2z_S_vrr;
      I_ERI_I2x4y_S_Gxy2z_S += I_ERI_I2x4y_S_Gxy2z_S_vrr;
      I_ERI_I2x3yz_S_Gxy2z_S += I_ERI_I2x3yz_S_Gxy2z_S_vrr;
      I_ERI_I2x2y2z_S_Gxy2z_S += I_ERI_I2x2y2z_S_Gxy2z_S_vrr;
      I_ERI_I2xy3z_S_Gxy2z_S += I_ERI_I2xy3z_S_Gxy2z_S_vrr;
      I_ERI_I2x4z_S_Gxy2z_S += I_ERI_I2x4z_S_Gxy2z_S_vrr;
      I_ERI_Ix5y_S_Gxy2z_S += I_ERI_Ix5y_S_Gxy2z_S_vrr;
      I_ERI_Ix4yz_S_Gxy2z_S += I_ERI_Ix4yz_S_Gxy2z_S_vrr;
      I_ERI_Ix3y2z_S_Gxy2z_S += I_ERI_Ix3y2z_S_Gxy2z_S_vrr;
      I_ERI_Ix2y3z_S_Gxy2z_S += I_ERI_Ix2y3z_S_Gxy2z_S_vrr;
      I_ERI_Ixy4z_S_Gxy2z_S += I_ERI_Ixy4z_S_Gxy2z_S_vrr;
      I_ERI_Ix5z_S_Gxy2z_S += I_ERI_Ix5z_S_Gxy2z_S_vrr;
      I_ERI_I6y_S_Gxy2z_S += I_ERI_I6y_S_Gxy2z_S_vrr;
      I_ERI_I5yz_S_Gxy2z_S += I_ERI_I5yz_S_Gxy2z_S_vrr;
      I_ERI_I4y2z_S_Gxy2z_S += I_ERI_I4y2z_S_Gxy2z_S_vrr;
      I_ERI_I3y3z_S_Gxy2z_S += I_ERI_I3y3z_S_Gxy2z_S_vrr;
      I_ERI_I2y4z_S_Gxy2z_S += I_ERI_I2y4z_S_Gxy2z_S_vrr;
      I_ERI_Iy5z_S_Gxy2z_S += I_ERI_Iy5z_S_Gxy2z_S_vrr;
      I_ERI_I6z_S_Gxy2z_S += I_ERI_I6z_S_Gxy2z_S_vrr;
      I_ERI_I6x_S_Gx3z_S += I_ERI_I6x_S_Gx3z_S_vrr;
      I_ERI_I5xy_S_Gx3z_S += I_ERI_I5xy_S_Gx3z_S_vrr;
      I_ERI_I5xz_S_Gx3z_S += I_ERI_I5xz_S_Gx3z_S_vrr;
      I_ERI_I4x2y_S_Gx3z_S += I_ERI_I4x2y_S_Gx3z_S_vrr;
      I_ERI_I4xyz_S_Gx3z_S += I_ERI_I4xyz_S_Gx3z_S_vrr;
      I_ERI_I4x2z_S_Gx3z_S += I_ERI_I4x2z_S_Gx3z_S_vrr;
      I_ERI_I3x3y_S_Gx3z_S += I_ERI_I3x3y_S_Gx3z_S_vrr;
      I_ERI_I3x2yz_S_Gx3z_S += I_ERI_I3x2yz_S_Gx3z_S_vrr;
      I_ERI_I3xy2z_S_Gx3z_S += I_ERI_I3xy2z_S_Gx3z_S_vrr;
      I_ERI_I3x3z_S_Gx3z_S += I_ERI_I3x3z_S_Gx3z_S_vrr;
      I_ERI_I2x4y_S_Gx3z_S += I_ERI_I2x4y_S_Gx3z_S_vrr;
      I_ERI_I2x3yz_S_Gx3z_S += I_ERI_I2x3yz_S_Gx3z_S_vrr;
      I_ERI_I2x2y2z_S_Gx3z_S += I_ERI_I2x2y2z_S_Gx3z_S_vrr;
      I_ERI_I2xy3z_S_Gx3z_S += I_ERI_I2xy3z_S_Gx3z_S_vrr;
      I_ERI_I2x4z_S_Gx3z_S += I_ERI_I2x4z_S_Gx3z_S_vrr;
      I_ERI_Ix5y_S_Gx3z_S += I_ERI_Ix5y_S_Gx3z_S_vrr;
      I_ERI_Ix4yz_S_Gx3z_S += I_ERI_Ix4yz_S_Gx3z_S_vrr;
      I_ERI_Ix3y2z_S_Gx3z_S += I_ERI_Ix3y2z_S_Gx3z_S_vrr;
      I_ERI_Ix2y3z_S_Gx3z_S += I_ERI_Ix2y3z_S_Gx3z_S_vrr;
      I_ERI_Ixy4z_S_Gx3z_S += I_ERI_Ixy4z_S_Gx3z_S_vrr;
      I_ERI_Ix5z_S_Gx3z_S += I_ERI_Ix5z_S_Gx3z_S_vrr;
      I_ERI_I6y_S_Gx3z_S += I_ERI_I6y_S_Gx3z_S_vrr;
      I_ERI_I5yz_S_Gx3z_S += I_ERI_I5yz_S_Gx3z_S_vrr;
      I_ERI_I4y2z_S_Gx3z_S += I_ERI_I4y2z_S_Gx3z_S_vrr;
      I_ERI_I3y3z_S_Gx3z_S += I_ERI_I3y3z_S_Gx3z_S_vrr;
      I_ERI_I2y4z_S_Gx3z_S += I_ERI_I2y4z_S_Gx3z_S_vrr;
      I_ERI_Iy5z_S_Gx3z_S += I_ERI_Iy5z_S_Gx3z_S_vrr;
      I_ERI_I6z_S_Gx3z_S += I_ERI_I6z_S_Gx3z_S_vrr;
      I_ERI_I6x_S_G4y_S += I_ERI_I6x_S_G4y_S_vrr;
      I_ERI_I5xy_S_G4y_S += I_ERI_I5xy_S_G4y_S_vrr;
      I_ERI_I5xz_S_G4y_S += I_ERI_I5xz_S_G4y_S_vrr;
      I_ERI_I4x2y_S_G4y_S += I_ERI_I4x2y_S_G4y_S_vrr;
      I_ERI_I4xyz_S_G4y_S += I_ERI_I4xyz_S_G4y_S_vrr;
      I_ERI_I4x2z_S_G4y_S += I_ERI_I4x2z_S_G4y_S_vrr;
      I_ERI_I3x3y_S_G4y_S += I_ERI_I3x3y_S_G4y_S_vrr;
      I_ERI_I3x2yz_S_G4y_S += I_ERI_I3x2yz_S_G4y_S_vrr;
      I_ERI_I3xy2z_S_G4y_S += I_ERI_I3xy2z_S_G4y_S_vrr;
      I_ERI_I3x3z_S_G4y_S += I_ERI_I3x3z_S_G4y_S_vrr;
      I_ERI_I2x4y_S_G4y_S += I_ERI_I2x4y_S_G4y_S_vrr;
      I_ERI_I2x3yz_S_G4y_S += I_ERI_I2x3yz_S_G4y_S_vrr;
      I_ERI_I2x2y2z_S_G4y_S += I_ERI_I2x2y2z_S_G4y_S_vrr;
      I_ERI_I2xy3z_S_G4y_S += I_ERI_I2xy3z_S_G4y_S_vrr;
      I_ERI_I2x4z_S_G4y_S += I_ERI_I2x4z_S_G4y_S_vrr;
      I_ERI_Ix5y_S_G4y_S += I_ERI_Ix5y_S_G4y_S_vrr;
      I_ERI_Ix4yz_S_G4y_S += I_ERI_Ix4yz_S_G4y_S_vrr;
      I_ERI_Ix3y2z_S_G4y_S += I_ERI_Ix3y2z_S_G4y_S_vrr;
      I_ERI_Ix2y3z_S_G4y_S += I_ERI_Ix2y3z_S_G4y_S_vrr;
      I_ERI_Ixy4z_S_G4y_S += I_ERI_Ixy4z_S_G4y_S_vrr;
      I_ERI_Ix5z_S_G4y_S += I_ERI_Ix5z_S_G4y_S_vrr;
      I_ERI_I6y_S_G4y_S += I_ERI_I6y_S_G4y_S_vrr;
      I_ERI_I5yz_S_G4y_S += I_ERI_I5yz_S_G4y_S_vrr;
      I_ERI_I4y2z_S_G4y_S += I_ERI_I4y2z_S_G4y_S_vrr;
      I_ERI_I3y3z_S_G4y_S += I_ERI_I3y3z_S_G4y_S_vrr;
      I_ERI_I2y4z_S_G4y_S += I_ERI_I2y4z_S_G4y_S_vrr;
      I_ERI_Iy5z_S_G4y_S += I_ERI_Iy5z_S_G4y_S_vrr;
      I_ERI_I6z_S_G4y_S += I_ERI_I6z_S_G4y_S_vrr;
      I_ERI_I6x_S_G3yz_S += I_ERI_I6x_S_G3yz_S_vrr;
      I_ERI_I5xy_S_G3yz_S += I_ERI_I5xy_S_G3yz_S_vrr;
      I_ERI_I5xz_S_G3yz_S += I_ERI_I5xz_S_G3yz_S_vrr;
      I_ERI_I4x2y_S_G3yz_S += I_ERI_I4x2y_S_G3yz_S_vrr;
      I_ERI_I4xyz_S_G3yz_S += I_ERI_I4xyz_S_G3yz_S_vrr;
      I_ERI_I4x2z_S_G3yz_S += I_ERI_I4x2z_S_G3yz_S_vrr;
      I_ERI_I3x3y_S_G3yz_S += I_ERI_I3x3y_S_G3yz_S_vrr;
      I_ERI_I3x2yz_S_G3yz_S += I_ERI_I3x2yz_S_G3yz_S_vrr;
      I_ERI_I3xy2z_S_G3yz_S += I_ERI_I3xy2z_S_G3yz_S_vrr;
      I_ERI_I3x3z_S_G3yz_S += I_ERI_I3x3z_S_G3yz_S_vrr;
      I_ERI_I2x4y_S_G3yz_S += I_ERI_I2x4y_S_G3yz_S_vrr;
      I_ERI_I2x3yz_S_G3yz_S += I_ERI_I2x3yz_S_G3yz_S_vrr;
      I_ERI_I2x2y2z_S_G3yz_S += I_ERI_I2x2y2z_S_G3yz_S_vrr;
      I_ERI_I2xy3z_S_G3yz_S += I_ERI_I2xy3z_S_G3yz_S_vrr;
      I_ERI_I2x4z_S_G3yz_S += I_ERI_I2x4z_S_G3yz_S_vrr;
      I_ERI_Ix5y_S_G3yz_S += I_ERI_Ix5y_S_G3yz_S_vrr;
      I_ERI_Ix4yz_S_G3yz_S += I_ERI_Ix4yz_S_G3yz_S_vrr;
      I_ERI_Ix3y2z_S_G3yz_S += I_ERI_Ix3y2z_S_G3yz_S_vrr;
      I_ERI_Ix2y3z_S_G3yz_S += I_ERI_Ix2y3z_S_G3yz_S_vrr;
      I_ERI_Ixy4z_S_G3yz_S += I_ERI_Ixy4z_S_G3yz_S_vrr;
      I_ERI_Ix5z_S_G3yz_S += I_ERI_Ix5z_S_G3yz_S_vrr;
      I_ERI_I6y_S_G3yz_S += I_ERI_I6y_S_G3yz_S_vrr;
      I_ERI_I5yz_S_G3yz_S += I_ERI_I5yz_S_G3yz_S_vrr;
      I_ERI_I4y2z_S_G3yz_S += I_ERI_I4y2z_S_G3yz_S_vrr;
      I_ERI_I3y3z_S_G3yz_S += I_ERI_I3y3z_S_G3yz_S_vrr;
      I_ERI_I2y4z_S_G3yz_S += I_ERI_I2y4z_S_G3yz_S_vrr;
      I_ERI_Iy5z_S_G3yz_S += I_ERI_Iy5z_S_G3yz_S_vrr;
      I_ERI_I6z_S_G3yz_S += I_ERI_I6z_S_G3yz_S_vrr;
      I_ERI_I6x_S_G2y2z_S += I_ERI_I6x_S_G2y2z_S_vrr;
      I_ERI_I5xy_S_G2y2z_S += I_ERI_I5xy_S_G2y2z_S_vrr;
      I_ERI_I5xz_S_G2y2z_S += I_ERI_I5xz_S_G2y2z_S_vrr;
      I_ERI_I4x2y_S_G2y2z_S += I_ERI_I4x2y_S_G2y2z_S_vrr;
      I_ERI_I4xyz_S_G2y2z_S += I_ERI_I4xyz_S_G2y2z_S_vrr;
      I_ERI_I4x2z_S_G2y2z_S += I_ERI_I4x2z_S_G2y2z_S_vrr;
      I_ERI_I3x3y_S_G2y2z_S += I_ERI_I3x3y_S_G2y2z_S_vrr;
      I_ERI_I3x2yz_S_G2y2z_S += I_ERI_I3x2yz_S_G2y2z_S_vrr;
      I_ERI_I3xy2z_S_G2y2z_S += I_ERI_I3xy2z_S_G2y2z_S_vrr;
      I_ERI_I3x3z_S_G2y2z_S += I_ERI_I3x3z_S_G2y2z_S_vrr;
      I_ERI_I2x4y_S_G2y2z_S += I_ERI_I2x4y_S_G2y2z_S_vrr;
      I_ERI_I2x3yz_S_G2y2z_S += I_ERI_I2x3yz_S_G2y2z_S_vrr;
      I_ERI_I2x2y2z_S_G2y2z_S += I_ERI_I2x2y2z_S_G2y2z_S_vrr;
      I_ERI_I2xy3z_S_G2y2z_S += I_ERI_I2xy3z_S_G2y2z_S_vrr;
      I_ERI_I2x4z_S_G2y2z_S += I_ERI_I2x4z_S_G2y2z_S_vrr;
      I_ERI_Ix5y_S_G2y2z_S += I_ERI_Ix5y_S_G2y2z_S_vrr;
      I_ERI_Ix4yz_S_G2y2z_S += I_ERI_Ix4yz_S_G2y2z_S_vrr;
      I_ERI_Ix3y2z_S_G2y2z_S += I_ERI_Ix3y2z_S_G2y2z_S_vrr;
      I_ERI_Ix2y3z_S_G2y2z_S += I_ERI_Ix2y3z_S_G2y2z_S_vrr;
      I_ERI_Ixy4z_S_G2y2z_S += I_ERI_Ixy4z_S_G2y2z_S_vrr;
      I_ERI_Ix5z_S_G2y2z_S += I_ERI_Ix5z_S_G2y2z_S_vrr;
      I_ERI_I6y_S_G2y2z_S += I_ERI_I6y_S_G2y2z_S_vrr;
      I_ERI_I5yz_S_G2y2z_S += I_ERI_I5yz_S_G2y2z_S_vrr;
      I_ERI_I4y2z_S_G2y2z_S += I_ERI_I4y2z_S_G2y2z_S_vrr;
      I_ERI_I3y3z_S_G2y2z_S += I_ERI_I3y3z_S_G2y2z_S_vrr;
      I_ERI_I2y4z_S_G2y2z_S += I_ERI_I2y4z_S_G2y2z_S_vrr;
      I_ERI_Iy5z_S_G2y2z_S += I_ERI_Iy5z_S_G2y2z_S_vrr;
      I_ERI_I6z_S_G2y2z_S += I_ERI_I6z_S_G2y2z_S_vrr;
      I_ERI_I6x_S_Gy3z_S += I_ERI_I6x_S_Gy3z_S_vrr;
      I_ERI_I5xy_S_Gy3z_S += I_ERI_I5xy_S_Gy3z_S_vrr;
      I_ERI_I5xz_S_Gy3z_S += I_ERI_I5xz_S_Gy3z_S_vrr;
      I_ERI_I4x2y_S_Gy3z_S += I_ERI_I4x2y_S_Gy3z_S_vrr;
      I_ERI_I4xyz_S_Gy3z_S += I_ERI_I4xyz_S_Gy3z_S_vrr;
      I_ERI_I4x2z_S_Gy3z_S += I_ERI_I4x2z_S_Gy3z_S_vrr;
      I_ERI_I3x3y_S_Gy3z_S += I_ERI_I3x3y_S_Gy3z_S_vrr;
      I_ERI_I3x2yz_S_Gy3z_S += I_ERI_I3x2yz_S_Gy3z_S_vrr;
      I_ERI_I3xy2z_S_Gy3z_S += I_ERI_I3xy2z_S_Gy3z_S_vrr;
      I_ERI_I3x3z_S_Gy3z_S += I_ERI_I3x3z_S_Gy3z_S_vrr;
      I_ERI_I2x4y_S_Gy3z_S += I_ERI_I2x4y_S_Gy3z_S_vrr;
      I_ERI_I2x3yz_S_Gy3z_S += I_ERI_I2x3yz_S_Gy3z_S_vrr;
      I_ERI_I2x2y2z_S_Gy3z_S += I_ERI_I2x2y2z_S_Gy3z_S_vrr;
      I_ERI_I2xy3z_S_Gy3z_S += I_ERI_I2xy3z_S_Gy3z_S_vrr;
      I_ERI_I2x4z_S_Gy3z_S += I_ERI_I2x4z_S_Gy3z_S_vrr;
      I_ERI_Ix5y_S_Gy3z_S += I_ERI_Ix5y_S_Gy3z_S_vrr;
      I_ERI_Ix4yz_S_Gy3z_S += I_ERI_Ix4yz_S_Gy3z_S_vrr;
      I_ERI_Ix3y2z_S_Gy3z_S += I_ERI_Ix3y2z_S_Gy3z_S_vrr;
      I_ERI_Ix2y3z_S_Gy3z_S += I_ERI_Ix2y3z_S_Gy3z_S_vrr;
      I_ERI_Ixy4z_S_Gy3z_S += I_ERI_Ixy4z_S_Gy3z_S_vrr;
      I_ERI_Ix5z_S_Gy3z_S += I_ERI_Ix5z_S_Gy3z_S_vrr;
      I_ERI_I6y_S_Gy3z_S += I_ERI_I6y_S_Gy3z_S_vrr;
      I_ERI_I5yz_S_Gy3z_S += I_ERI_I5yz_S_Gy3z_S_vrr;
      I_ERI_I4y2z_S_Gy3z_S += I_ERI_I4y2z_S_Gy3z_S_vrr;
      I_ERI_I3y3z_S_Gy3z_S += I_ERI_I3y3z_S_Gy3z_S_vrr;
      I_ERI_I2y4z_S_Gy3z_S += I_ERI_I2y4z_S_Gy3z_S_vrr;
      I_ERI_Iy5z_S_Gy3z_S += I_ERI_Iy5z_S_Gy3z_S_vrr;
      I_ERI_I6z_S_Gy3z_S += I_ERI_I6z_S_Gy3z_S_vrr;
      I_ERI_I6x_S_G4z_S += I_ERI_I6x_S_G4z_S_vrr;
      I_ERI_I5xy_S_G4z_S += I_ERI_I5xy_S_G4z_S_vrr;
      I_ERI_I5xz_S_G4z_S += I_ERI_I5xz_S_G4z_S_vrr;
      I_ERI_I4x2y_S_G4z_S += I_ERI_I4x2y_S_G4z_S_vrr;
      I_ERI_I4xyz_S_G4z_S += I_ERI_I4xyz_S_G4z_S_vrr;
      I_ERI_I4x2z_S_G4z_S += I_ERI_I4x2z_S_G4z_S_vrr;
      I_ERI_I3x3y_S_G4z_S += I_ERI_I3x3y_S_G4z_S_vrr;
      I_ERI_I3x2yz_S_G4z_S += I_ERI_I3x2yz_S_G4z_S_vrr;
      I_ERI_I3xy2z_S_G4z_S += I_ERI_I3xy2z_S_G4z_S_vrr;
      I_ERI_I3x3z_S_G4z_S += I_ERI_I3x3z_S_G4z_S_vrr;
      I_ERI_I2x4y_S_G4z_S += I_ERI_I2x4y_S_G4z_S_vrr;
      I_ERI_I2x3yz_S_G4z_S += I_ERI_I2x3yz_S_G4z_S_vrr;
      I_ERI_I2x2y2z_S_G4z_S += I_ERI_I2x2y2z_S_G4z_S_vrr;
      I_ERI_I2xy3z_S_G4z_S += I_ERI_I2xy3z_S_G4z_S_vrr;
      I_ERI_I2x4z_S_G4z_S += I_ERI_I2x4z_S_G4z_S_vrr;
      I_ERI_Ix5y_S_G4z_S += I_ERI_Ix5y_S_G4z_S_vrr;
      I_ERI_Ix4yz_S_G4z_S += I_ERI_Ix4yz_S_G4z_S_vrr;
      I_ERI_Ix3y2z_S_G4z_S += I_ERI_Ix3y2z_S_G4z_S_vrr;
      I_ERI_Ix2y3z_S_G4z_S += I_ERI_Ix2y3z_S_G4z_S_vrr;
      I_ERI_Ixy4z_S_G4z_S += I_ERI_Ixy4z_S_G4z_S_vrr;
      I_ERI_Ix5z_S_G4z_S += I_ERI_Ix5z_S_G4z_S_vrr;
      I_ERI_I6y_S_G4z_S += I_ERI_I6y_S_G4z_S_vrr;
      I_ERI_I5yz_S_G4z_S += I_ERI_I5yz_S_G4z_S_vrr;
      I_ERI_I4y2z_S_G4z_S += I_ERI_I4y2z_S_G4z_S_vrr;
      I_ERI_I3y3z_S_G4z_S += I_ERI_I3y3z_S_G4z_S_vrr;
      I_ERI_I2y4z_S_G4z_S += I_ERI_I2y4z_S_G4z_S_vrr;
      I_ERI_Iy5z_S_G4z_S += I_ERI_Iy5z_S_G4z_S_vrr;
      I_ERI_I6z_S_G4z_S += I_ERI_I6z_S_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_G_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_H5x_S_G4x_S += I_ERI_H5x_S_G4x_S_vrr;
      I_ERI_H4xy_S_G4x_S += I_ERI_H4xy_S_G4x_S_vrr;
      I_ERI_H4xz_S_G4x_S += I_ERI_H4xz_S_G4x_S_vrr;
      I_ERI_H3x2y_S_G4x_S += I_ERI_H3x2y_S_G4x_S_vrr;
      I_ERI_H3xyz_S_G4x_S += I_ERI_H3xyz_S_G4x_S_vrr;
      I_ERI_H3x2z_S_G4x_S += I_ERI_H3x2z_S_G4x_S_vrr;
      I_ERI_H2x3y_S_G4x_S += I_ERI_H2x3y_S_G4x_S_vrr;
      I_ERI_H2x2yz_S_G4x_S += I_ERI_H2x2yz_S_G4x_S_vrr;
      I_ERI_H2xy2z_S_G4x_S += I_ERI_H2xy2z_S_G4x_S_vrr;
      I_ERI_H2x3z_S_G4x_S += I_ERI_H2x3z_S_G4x_S_vrr;
      I_ERI_Hx4y_S_G4x_S += I_ERI_Hx4y_S_G4x_S_vrr;
      I_ERI_Hx3yz_S_G4x_S += I_ERI_Hx3yz_S_G4x_S_vrr;
      I_ERI_Hx2y2z_S_G4x_S += I_ERI_Hx2y2z_S_G4x_S_vrr;
      I_ERI_Hxy3z_S_G4x_S += I_ERI_Hxy3z_S_G4x_S_vrr;
      I_ERI_Hx4z_S_G4x_S += I_ERI_Hx4z_S_G4x_S_vrr;
      I_ERI_H5y_S_G4x_S += I_ERI_H5y_S_G4x_S_vrr;
      I_ERI_H4yz_S_G4x_S += I_ERI_H4yz_S_G4x_S_vrr;
      I_ERI_H3y2z_S_G4x_S += I_ERI_H3y2z_S_G4x_S_vrr;
      I_ERI_H2y3z_S_G4x_S += I_ERI_H2y3z_S_G4x_S_vrr;
      I_ERI_Hy4z_S_G4x_S += I_ERI_Hy4z_S_G4x_S_vrr;
      I_ERI_H5z_S_G4x_S += I_ERI_H5z_S_G4x_S_vrr;
      I_ERI_H5x_S_G3xy_S += I_ERI_H5x_S_G3xy_S_vrr;
      I_ERI_H4xy_S_G3xy_S += I_ERI_H4xy_S_G3xy_S_vrr;
      I_ERI_H4xz_S_G3xy_S += I_ERI_H4xz_S_G3xy_S_vrr;
      I_ERI_H3x2y_S_G3xy_S += I_ERI_H3x2y_S_G3xy_S_vrr;
      I_ERI_H3xyz_S_G3xy_S += I_ERI_H3xyz_S_G3xy_S_vrr;
      I_ERI_H3x2z_S_G3xy_S += I_ERI_H3x2z_S_G3xy_S_vrr;
      I_ERI_H2x3y_S_G3xy_S += I_ERI_H2x3y_S_G3xy_S_vrr;
      I_ERI_H2x2yz_S_G3xy_S += I_ERI_H2x2yz_S_G3xy_S_vrr;
      I_ERI_H2xy2z_S_G3xy_S += I_ERI_H2xy2z_S_G3xy_S_vrr;
      I_ERI_H2x3z_S_G3xy_S += I_ERI_H2x3z_S_G3xy_S_vrr;
      I_ERI_Hx4y_S_G3xy_S += I_ERI_Hx4y_S_G3xy_S_vrr;
      I_ERI_Hx3yz_S_G3xy_S += I_ERI_Hx3yz_S_G3xy_S_vrr;
      I_ERI_Hx2y2z_S_G3xy_S += I_ERI_Hx2y2z_S_G3xy_S_vrr;
      I_ERI_Hxy3z_S_G3xy_S += I_ERI_Hxy3z_S_G3xy_S_vrr;
      I_ERI_Hx4z_S_G3xy_S += I_ERI_Hx4z_S_G3xy_S_vrr;
      I_ERI_H5y_S_G3xy_S += I_ERI_H5y_S_G3xy_S_vrr;
      I_ERI_H4yz_S_G3xy_S += I_ERI_H4yz_S_G3xy_S_vrr;
      I_ERI_H3y2z_S_G3xy_S += I_ERI_H3y2z_S_G3xy_S_vrr;
      I_ERI_H2y3z_S_G3xy_S += I_ERI_H2y3z_S_G3xy_S_vrr;
      I_ERI_Hy4z_S_G3xy_S += I_ERI_Hy4z_S_G3xy_S_vrr;
      I_ERI_H5z_S_G3xy_S += I_ERI_H5z_S_G3xy_S_vrr;
      I_ERI_H5x_S_G3xz_S += I_ERI_H5x_S_G3xz_S_vrr;
      I_ERI_H4xy_S_G3xz_S += I_ERI_H4xy_S_G3xz_S_vrr;
      I_ERI_H4xz_S_G3xz_S += I_ERI_H4xz_S_G3xz_S_vrr;
      I_ERI_H3x2y_S_G3xz_S += I_ERI_H3x2y_S_G3xz_S_vrr;
      I_ERI_H3xyz_S_G3xz_S += I_ERI_H3xyz_S_G3xz_S_vrr;
      I_ERI_H3x2z_S_G3xz_S += I_ERI_H3x2z_S_G3xz_S_vrr;
      I_ERI_H2x3y_S_G3xz_S += I_ERI_H2x3y_S_G3xz_S_vrr;
      I_ERI_H2x2yz_S_G3xz_S += I_ERI_H2x2yz_S_G3xz_S_vrr;
      I_ERI_H2xy2z_S_G3xz_S += I_ERI_H2xy2z_S_G3xz_S_vrr;
      I_ERI_H2x3z_S_G3xz_S += I_ERI_H2x3z_S_G3xz_S_vrr;
      I_ERI_Hx4y_S_G3xz_S += I_ERI_Hx4y_S_G3xz_S_vrr;
      I_ERI_Hx3yz_S_G3xz_S += I_ERI_Hx3yz_S_G3xz_S_vrr;
      I_ERI_Hx2y2z_S_G3xz_S += I_ERI_Hx2y2z_S_G3xz_S_vrr;
      I_ERI_Hxy3z_S_G3xz_S += I_ERI_Hxy3z_S_G3xz_S_vrr;
      I_ERI_Hx4z_S_G3xz_S += I_ERI_Hx4z_S_G3xz_S_vrr;
      I_ERI_H5y_S_G3xz_S += I_ERI_H5y_S_G3xz_S_vrr;
      I_ERI_H4yz_S_G3xz_S += I_ERI_H4yz_S_G3xz_S_vrr;
      I_ERI_H3y2z_S_G3xz_S += I_ERI_H3y2z_S_G3xz_S_vrr;
      I_ERI_H2y3z_S_G3xz_S += I_ERI_H2y3z_S_G3xz_S_vrr;
      I_ERI_Hy4z_S_G3xz_S += I_ERI_Hy4z_S_G3xz_S_vrr;
      I_ERI_H5z_S_G3xz_S += I_ERI_H5z_S_G3xz_S_vrr;
      I_ERI_H5x_S_G2x2y_S += I_ERI_H5x_S_G2x2y_S_vrr;
      I_ERI_H4xy_S_G2x2y_S += I_ERI_H4xy_S_G2x2y_S_vrr;
      I_ERI_H4xz_S_G2x2y_S += I_ERI_H4xz_S_G2x2y_S_vrr;
      I_ERI_H3x2y_S_G2x2y_S += I_ERI_H3x2y_S_G2x2y_S_vrr;
      I_ERI_H3xyz_S_G2x2y_S += I_ERI_H3xyz_S_G2x2y_S_vrr;
      I_ERI_H3x2z_S_G2x2y_S += I_ERI_H3x2z_S_G2x2y_S_vrr;
      I_ERI_H2x3y_S_G2x2y_S += I_ERI_H2x3y_S_G2x2y_S_vrr;
      I_ERI_H2x2yz_S_G2x2y_S += I_ERI_H2x2yz_S_G2x2y_S_vrr;
      I_ERI_H2xy2z_S_G2x2y_S += I_ERI_H2xy2z_S_G2x2y_S_vrr;
      I_ERI_H2x3z_S_G2x2y_S += I_ERI_H2x3z_S_G2x2y_S_vrr;
      I_ERI_Hx4y_S_G2x2y_S += I_ERI_Hx4y_S_G2x2y_S_vrr;
      I_ERI_Hx3yz_S_G2x2y_S += I_ERI_Hx3yz_S_G2x2y_S_vrr;
      I_ERI_Hx2y2z_S_G2x2y_S += I_ERI_Hx2y2z_S_G2x2y_S_vrr;
      I_ERI_Hxy3z_S_G2x2y_S += I_ERI_Hxy3z_S_G2x2y_S_vrr;
      I_ERI_Hx4z_S_G2x2y_S += I_ERI_Hx4z_S_G2x2y_S_vrr;
      I_ERI_H5y_S_G2x2y_S += I_ERI_H5y_S_G2x2y_S_vrr;
      I_ERI_H4yz_S_G2x2y_S += I_ERI_H4yz_S_G2x2y_S_vrr;
      I_ERI_H3y2z_S_G2x2y_S += I_ERI_H3y2z_S_G2x2y_S_vrr;
      I_ERI_H2y3z_S_G2x2y_S += I_ERI_H2y3z_S_G2x2y_S_vrr;
      I_ERI_Hy4z_S_G2x2y_S += I_ERI_Hy4z_S_G2x2y_S_vrr;
      I_ERI_H5z_S_G2x2y_S += I_ERI_H5z_S_G2x2y_S_vrr;
      I_ERI_H5x_S_G2xyz_S += I_ERI_H5x_S_G2xyz_S_vrr;
      I_ERI_H4xy_S_G2xyz_S += I_ERI_H4xy_S_G2xyz_S_vrr;
      I_ERI_H4xz_S_G2xyz_S += I_ERI_H4xz_S_G2xyz_S_vrr;
      I_ERI_H3x2y_S_G2xyz_S += I_ERI_H3x2y_S_G2xyz_S_vrr;
      I_ERI_H3xyz_S_G2xyz_S += I_ERI_H3xyz_S_G2xyz_S_vrr;
      I_ERI_H3x2z_S_G2xyz_S += I_ERI_H3x2z_S_G2xyz_S_vrr;
      I_ERI_H2x3y_S_G2xyz_S += I_ERI_H2x3y_S_G2xyz_S_vrr;
      I_ERI_H2x2yz_S_G2xyz_S += I_ERI_H2x2yz_S_G2xyz_S_vrr;
      I_ERI_H2xy2z_S_G2xyz_S += I_ERI_H2xy2z_S_G2xyz_S_vrr;
      I_ERI_H2x3z_S_G2xyz_S += I_ERI_H2x3z_S_G2xyz_S_vrr;
      I_ERI_Hx4y_S_G2xyz_S += I_ERI_Hx4y_S_G2xyz_S_vrr;
      I_ERI_Hx3yz_S_G2xyz_S += I_ERI_Hx3yz_S_G2xyz_S_vrr;
      I_ERI_Hx2y2z_S_G2xyz_S += I_ERI_Hx2y2z_S_G2xyz_S_vrr;
      I_ERI_Hxy3z_S_G2xyz_S += I_ERI_Hxy3z_S_G2xyz_S_vrr;
      I_ERI_Hx4z_S_G2xyz_S += I_ERI_Hx4z_S_G2xyz_S_vrr;
      I_ERI_H5y_S_G2xyz_S += I_ERI_H5y_S_G2xyz_S_vrr;
      I_ERI_H4yz_S_G2xyz_S += I_ERI_H4yz_S_G2xyz_S_vrr;
      I_ERI_H3y2z_S_G2xyz_S += I_ERI_H3y2z_S_G2xyz_S_vrr;
      I_ERI_H2y3z_S_G2xyz_S += I_ERI_H2y3z_S_G2xyz_S_vrr;
      I_ERI_Hy4z_S_G2xyz_S += I_ERI_Hy4z_S_G2xyz_S_vrr;
      I_ERI_H5z_S_G2xyz_S += I_ERI_H5z_S_G2xyz_S_vrr;
      I_ERI_H5x_S_G2x2z_S += I_ERI_H5x_S_G2x2z_S_vrr;
      I_ERI_H4xy_S_G2x2z_S += I_ERI_H4xy_S_G2x2z_S_vrr;
      I_ERI_H4xz_S_G2x2z_S += I_ERI_H4xz_S_G2x2z_S_vrr;
      I_ERI_H3x2y_S_G2x2z_S += I_ERI_H3x2y_S_G2x2z_S_vrr;
      I_ERI_H3xyz_S_G2x2z_S += I_ERI_H3xyz_S_G2x2z_S_vrr;
      I_ERI_H3x2z_S_G2x2z_S += I_ERI_H3x2z_S_G2x2z_S_vrr;
      I_ERI_H2x3y_S_G2x2z_S += I_ERI_H2x3y_S_G2x2z_S_vrr;
      I_ERI_H2x2yz_S_G2x2z_S += I_ERI_H2x2yz_S_G2x2z_S_vrr;
      I_ERI_H2xy2z_S_G2x2z_S += I_ERI_H2xy2z_S_G2x2z_S_vrr;
      I_ERI_H2x3z_S_G2x2z_S += I_ERI_H2x3z_S_G2x2z_S_vrr;
      I_ERI_Hx4y_S_G2x2z_S += I_ERI_Hx4y_S_G2x2z_S_vrr;
      I_ERI_Hx3yz_S_G2x2z_S += I_ERI_Hx3yz_S_G2x2z_S_vrr;
      I_ERI_Hx2y2z_S_G2x2z_S += I_ERI_Hx2y2z_S_G2x2z_S_vrr;
      I_ERI_Hxy3z_S_G2x2z_S += I_ERI_Hxy3z_S_G2x2z_S_vrr;
      I_ERI_Hx4z_S_G2x2z_S += I_ERI_Hx4z_S_G2x2z_S_vrr;
      I_ERI_H5y_S_G2x2z_S += I_ERI_H5y_S_G2x2z_S_vrr;
      I_ERI_H4yz_S_G2x2z_S += I_ERI_H4yz_S_G2x2z_S_vrr;
      I_ERI_H3y2z_S_G2x2z_S += I_ERI_H3y2z_S_G2x2z_S_vrr;
      I_ERI_H2y3z_S_G2x2z_S += I_ERI_H2y3z_S_G2x2z_S_vrr;
      I_ERI_Hy4z_S_G2x2z_S += I_ERI_Hy4z_S_G2x2z_S_vrr;
      I_ERI_H5z_S_G2x2z_S += I_ERI_H5z_S_G2x2z_S_vrr;
      I_ERI_H5x_S_Gx3y_S += I_ERI_H5x_S_Gx3y_S_vrr;
      I_ERI_H4xy_S_Gx3y_S += I_ERI_H4xy_S_Gx3y_S_vrr;
      I_ERI_H4xz_S_Gx3y_S += I_ERI_H4xz_S_Gx3y_S_vrr;
      I_ERI_H3x2y_S_Gx3y_S += I_ERI_H3x2y_S_Gx3y_S_vrr;
      I_ERI_H3xyz_S_Gx3y_S += I_ERI_H3xyz_S_Gx3y_S_vrr;
      I_ERI_H3x2z_S_Gx3y_S += I_ERI_H3x2z_S_Gx3y_S_vrr;
      I_ERI_H2x3y_S_Gx3y_S += I_ERI_H2x3y_S_Gx3y_S_vrr;
      I_ERI_H2x2yz_S_Gx3y_S += I_ERI_H2x2yz_S_Gx3y_S_vrr;
      I_ERI_H2xy2z_S_Gx3y_S += I_ERI_H2xy2z_S_Gx3y_S_vrr;
      I_ERI_H2x3z_S_Gx3y_S += I_ERI_H2x3z_S_Gx3y_S_vrr;
      I_ERI_Hx4y_S_Gx3y_S += I_ERI_Hx4y_S_Gx3y_S_vrr;
      I_ERI_Hx3yz_S_Gx3y_S += I_ERI_Hx3yz_S_Gx3y_S_vrr;
      I_ERI_Hx2y2z_S_Gx3y_S += I_ERI_Hx2y2z_S_Gx3y_S_vrr;
      I_ERI_Hxy3z_S_Gx3y_S += I_ERI_Hxy3z_S_Gx3y_S_vrr;
      I_ERI_Hx4z_S_Gx3y_S += I_ERI_Hx4z_S_Gx3y_S_vrr;
      I_ERI_H5y_S_Gx3y_S += I_ERI_H5y_S_Gx3y_S_vrr;
      I_ERI_H4yz_S_Gx3y_S += I_ERI_H4yz_S_Gx3y_S_vrr;
      I_ERI_H3y2z_S_Gx3y_S += I_ERI_H3y2z_S_Gx3y_S_vrr;
      I_ERI_H2y3z_S_Gx3y_S += I_ERI_H2y3z_S_Gx3y_S_vrr;
      I_ERI_Hy4z_S_Gx3y_S += I_ERI_Hy4z_S_Gx3y_S_vrr;
      I_ERI_H5z_S_Gx3y_S += I_ERI_H5z_S_Gx3y_S_vrr;
      I_ERI_H5x_S_Gx2yz_S += I_ERI_H5x_S_Gx2yz_S_vrr;
      I_ERI_H4xy_S_Gx2yz_S += I_ERI_H4xy_S_Gx2yz_S_vrr;
      I_ERI_H4xz_S_Gx2yz_S += I_ERI_H4xz_S_Gx2yz_S_vrr;
      I_ERI_H3x2y_S_Gx2yz_S += I_ERI_H3x2y_S_Gx2yz_S_vrr;
      I_ERI_H3xyz_S_Gx2yz_S += I_ERI_H3xyz_S_Gx2yz_S_vrr;
      I_ERI_H3x2z_S_Gx2yz_S += I_ERI_H3x2z_S_Gx2yz_S_vrr;
      I_ERI_H2x3y_S_Gx2yz_S += I_ERI_H2x3y_S_Gx2yz_S_vrr;
      I_ERI_H2x2yz_S_Gx2yz_S += I_ERI_H2x2yz_S_Gx2yz_S_vrr;
      I_ERI_H2xy2z_S_Gx2yz_S += I_ERI_H2xy2z_S_Gx2yz_S_vrr;
      I_ERI_H2x3z_S_Gx2yz_S += I_ERI_H2x3z_S_Gx2yz_S_vrr;
      I_ERI_Hx4y_S_Gx2yz_S += I_ERI_Hx4y_S_Gx2yz_S_vrr;
      I_ERI_Hx3yz_S_Gx2yz_S += I_ERI_Hx3yz_S_Gx2yz_S_vrr;
      I_ERI_Hx2y2z_S_Gx2yz_S += I_ERI_Hx2y2z_S_Gx2yz_S_vrr;
      I_ERI_Hxy3z_S_Gx2yz_S += I_ERI_Hxy3z_S_Gx2yz_S_vrr;
      I_ERI_Hx4z_S_Gx2yz_S += I_ERI_Hx4z_S_Gx2yz_S_vrr;
      I_ERI_H5y_S_Gx2yz_S += I_ERI_H5y_S_Gx2yz_S_vrr;
      I_ERI_H4yz_S_Gx2yz_S += I_ERI_H4yz_S_Gx2yz_S_vrr;
      I_ERI_H3y2z_S_Gx2yz_S += I_ERI_H3y2z_S_Gx2yz_S_vrr;
      I_ERI_H2y3z_S_Gx2yz_S += I_ERI_H2y3z_S_Gx2yz_S_vrr;
      I_ERI_Hy4z_S_Gx2yz_S += I_ERI_Hy4z_S_Gx2yz_S_vrr;
      I_ERI_H5z_S_Gx2yz_S += I_ERI_H5z_S_Gx2yz_S_vrr;
      I_ERI_H5x_S_Gxy2z_S += I_ERI_H5x_S_Gxy2z_S_vrr;
      I_ERI_H4xy_S_Gxy2z_S += I_ERI_H4xy_S_Gxy2z_S_vrr;
      I_ERI_H4xz_S_Gxy2z_S += I_ERI_H4xz_S_Gxy2z_S_vrr;
      I_ERI_H3x2y_S_Gxy2z_S += I_ERI_H3x2y_S_Gxy2z_S_vrr;
      I_ERI_H3xyz_S_Gxy2z_S += I_ERI_H3xyz_S_Gxy2z_S_vrr;
      I_ERI_H3x2z_S_Gxy2z_S += I_ERI_H3x2z_S_Gxy2z_S_vrr;
      I_ERI_H2x3y_S_Gxy2z_S += I_ERI_H2x3y_S_Gxy2z_S_vrr;
      I_ERI_H2x2yz_S_Gxy2z_S += I_ERI_H2x2yz_S_Gxy2z_S_vrr;
      I_ERI_H2xy2z_S_Gxy2z_S += I_ERI_H2xy2z_S_Gxy2z_S_vrr;
      I_ERI_H2x3z_S_Gxy2z_S += I_ERI_H2x3z_S_Gxy2z_S_vrr;
      I_ERI_Hx4y_S_Gxy2z_S += I_ERI_Hx4y_S_Gxy2z_S_vrr;
      I_ERI_Hx3yz_S_Gxy2z_S += I_ERI_Hx3yz_S_Gxy2z_S_vrr;
      I_ERI_Hx2y2z_S_Gxy2z_S += I_ERI_Hx2y2z_S_Gxy2z_S_vrr;
      I_ERI_Hxy3z_S_Gxy2z_S += I_ERI_Hxy3z_S_Gxy2z_S_vrr;
      I_ERI_Hx4z_S_Gxy2z_S += I_ERI_Hx4z_S_Gxy2z_S_vrr;
      I_ERI_H5y_S_Gxy2z_S += I_ERI_H5y_S_Gxy2z_S_vrr;
      I_ERI_H4yz_S_Gxy2z_S += I_ERI_H4yz_S_Gxy2z_S_vrr;
      I_ERI_H3y2z_S_Gxy2z_S += I_ERI_H3y2z_S_Gxy2z_S_vrr;
      I_ERI_H2y3z_S_Gxy2z_S += I_ERI_H2y3z_S_Gxy2z_S_vrr;
      I_ERI_Hy4z_S_Gxy2z_S += I_ERI_Hy4z_S_Gxy2z_S_vrr;
      I_ERI_H5z_S_Gxy2z_S += I_ERI_H5z_S_Gxy2z_S_vrr;
      I_ERI_H5x_S_Gx3z_S += I_ERI_H5x_S_Gx3z_S_vrr;
      I_ERI_H4xy_S_Gx3z_S += I_ERI_H4xy_S_Gx3z_S_vrr;
      I_ERI_H4xz_S_Gx3z_S += I_ERI_H4xz_S_Gx3z_S_vrr;
      I_ERI_H3x2y_S_Gx3z_S += I_ERI_H3x2y_S_Gx3z_S_vrr;
      I_ERI_H3xyz_S_Gx3z_S += I_ERI_H3xyz_S_Gx3z_S_vrr;
      I_ERI_H3x2z_S_Gx3z_S += I_ERI_H3x2z_S_Gx3z_S_vrr;
      I_ERI_H2x3y_S_Gx3z_S += I_ERI_H2x3y_S_Gx3z_S_vrr;
      I_ERI_H2x2yz_S_Gx3z_S += I_ERI_H2x2yz_S_Gx3z_S_vrr;
      I_ERI_H2xy2z_S_Gx3z_S += I_ERI_H2xy2z_S_Gx3z_S_vrr;
      I_ERI_H2x3z_S_Gx3z_S += I_ERI_H2x3z_S_Gx3z_S_vrr;
      I_ERI_Hx4y_S_Gx3z_S += I_ERI_Hx4y_S_Gx3z_S_vrr;
      I_ERI_Hx3yz_S_Gx3z_S += I_ERI_Hx3yz_S_Gx3z_S_vrr;
      I_ERI_Hx2y2z_S_Gx3z_S += I_ERI_Hx2y2z_S_Gx3z_S_vrr;
      I_ERI_Hxy3z_S_Gx3z_S += I_ERI_Hxy3z_S_Gx3z_S_vrr;
      I_ERI_Hx4z_S_Gx3z_S += I_ERI_Hx4z_S_Gx3z_S_vrr;
      I_ERI_H5y_S_Gx3z_S += I_ERI_H5y_S_Gx3z_S_vrr;
      I_ERI_H4yz_S_Gx3z_S += I_ERI_H4yz_S_Gx3z_S_vrr;
      I_ERI_H3y2z_S_Gx3z_S += I_ERI_H3y2z_S_Gx3z_S_vrr;
      I_ERI_H2y3z_S_Gx3z_S += I_ERI_H2y3z_S_Gx3z_S_vrr;
      I_ERI_Hy4z_S_Gx3z_S += I_ERI_Hy4z_S_Gx3z_S_vrr;
      I_ERI_H5z_S_Gx3z_S += I_ERI_H5z_S_Gx3z_S_vrr;
      I_ERI_H5x_S_G4y_S += I_ERI_H5x_S_G4y_S_vrr;
      I_ERI_H4xy_S_G4y_S += I_ERI_H4xy_S_G4y_S_vrr;
      I_ERI_H4xz_S_G4y_S += I_ERI_H4xz_S_G4y_S_vrr;
      I_ERI_H3x2y_S_G4y_S += I_ERI_H3x2y_S_G4y_S_vrr;
      I_ERI_H3xyz_S_G4y_S += I_ERI_H3xyz_S_G4y_S_vrr;
      I_ERI_H3x2z_S_G4y_S += I_ERI_H3x2z_S_G4y_S_vrr;
      I_ERI_H2x3y_S_G4y_S += I_ERI_H2x3y_S_G4y_S_vrr;
      I_ERI_H2x2yz_S_G4y_S += I_ERI_H2x2yz_S_G4y_S_vrr;
      I_ERI_H2xy2z_S_G4y_S += I_ERI_H2xy2z_S_G4y_S_vrr;
      I_ERI_H2x3z_S_G4y_S += I_ERI_H2x3z_S_G4y_S_vrr;
      I_ERI_Hx4y_S_G4y_S += I_ERI_Hx4y_S_G4y_S_vrr;
      I_ERI_Hx3yz_S_G4y_S += I_ERI_Hx3yz_S_G4y_S_vrr;
      I_ERI_Hx2y2z_S_G4y_S += I_ERI_Hx2y2z_S_G4y_S_vrr;
      I_ERI_Hxy3z_S_G4y_S += I_ERI_Hxy3z_S_G4y_S_vrr;
      I_ERI_Hx4z_S_G4y_S += I_ERI_Hx4z_S_G4y_S_vrr;
      I_ERI_H5y_S_G4y_S += I_ERI_H5y_S_G4y_S_vrr;
      I_ERI_H4yz_S_G4y_S += I_ERI_H4yz_S_G4y_S_vrr;
      I_ERI_H3y2z_S_G4y_S += I_ERI_H3y2z_S_G4y_S_vrr;
      I_ERI_H2y3z_S_G4y_S += I_ERI_H2y3z_S_G4y_S_vrr;
      I_ERI_Hy4z_S_G4y_S += I_ERI_Hy4z_S_G4y_S_vrr;
      I_ERI_H5z_S_G4y_S += I_ERI_H5z_S_G4y_S_vrr;
      I_ERI_H5x_S_G3yz_S += I_ERI_H5x_S_G3yz_S_vrr;
      I_ERI_H4xy_S_G3yz_S += I_ERI_H4xy_S_G3yz_S_vrr;
      I_ERI_H4xz_S_G3yz_S += I_ERI_H4xz_S_G3yz_S_vrr;
      I_ERI_H3x2y_S_G3yz_S += I_ERI_H3x2y_S_G3yz_S_vrr;
      I_ERI_H3xyz_S_G3yz_S += I_ERI_H3xyz_S_G3yz_S_vrr;
      I_ERI_H3x2z_S_G3yz_S += I_ERI_H3x2z_S_G3yz_S_vrr;
      I_ERI_H2x3y_S_G3yz_S += I_ERI_H2x3y_S_G3yz_S_vrr;
      I_ERI_H2x2yz_S_G3yz_S += I_ERI_H2x2yz_S_G3yz_S_vrr;
      I_ERI_H2xy2z_S_G3yz_S += I_ERI_H2xy2z_S_G3yz_S_vrr;
      I_ERI_H2x3z_S_G3yz_S += I_ERI_H2x3z_S_G3yz_S_vrr;
      I_ERI_Hx4y_S_G3yz_S += I_ERI_Hx4y_S_G3yz_S_vrr;
      I_ERI_Hx3yz_S_G3yz_S += I_ERI_Hx3yz_S_G3yz_S_vrr;
      I_ERI_Hx2y2z_S_G3yz_S += I_ERI_Hx2y2z_S_G3yz_S_vrr;
      I_ERI_Hxy3z_S_G3yz_S += I_ERI_Hxy3z_S_G3yz_S_vrr;
      I_ERI_Hx4z_S_G3yz_S += I_ERI_Hx4z_S_G3yz_S_vrr;
      I_ERI_H5y_S_G3yz_S += I_ERI_H5y_S_G3yz_S_vrr;
      I_ERI_H4yz_S_G3yz_S += I_ERI_H4yz_S_G3yz_S_vrr;
      I_ERI_H3y2z_S_G3yz_S += I_ERI_H3y2z_S_G3yz_S_vrr;
      I_ERI_H2y3z_S_G3yz_S += I_ERI_H2y3z_S_G3yz_S_vrr;
      I_ERI_Hy4z_S_G3yz_S += I_ERI_Hy4z_S_G3yz_S_vrr;
      I_ERI_H5z_S_G3yz_S += I_ERI_H5z_S_G3yz_S_vrr;
      I_ERI_H5x_S_G2y2z_S += I_ERI_H5x_S_G2y2z_S_vrr;
      I_ERI_H4xy_S_G2y2z_S += I_ERI_H4xy_S_G2y2z_S_vrr;
      I_ERI_H4xz_S_G2y2z_S += I_ERI_H4xz_S_G2y2z_S_vrr;
      I_ERI_H3x2y_S_G2y2z_S += I_ERI_H3x2y_S_G2y2z_S_vrr;
      I_ERI_H3xyz_S_G2y2z_S += I_ERI_H3xyz_S_G2y2z_S_vrr;
      I_ERI_H3x2z_S_G2y2z_S += I_ERI_H3x2z_S_G2y2z_S_vrr;
      I_ERI_H2x3y_S_G2y2z_S += I_ERI_H2x3y_S_G2y2z_S_vrr;
      I_ERI_H2x2yz_S_G2y2z_S += I_ERI_H2x2yz_S_G2y2z_S_vrr;
      I_ERI_H2xy2z_S_G2y2z_S += I_ERI_H2xy2z_S_G2y2z_S_vrr;
      I_ERI_H2x3z_S_G2y2z_S += I_ERI_H2x3z_S_G2y2z_S_vrr;
      I_ERI_Hx4y_S_G2y2z_S += I_ERI_Hx4y_S_G2y2z_S_vrr;
      I_ERI_Hx3yz_S_G2y2z_S += I_ERI_Hx3yz_S_G2y2z_S_vrr;
      I_ERI_Hx2y2z_S_G2y2z_S += I_ERI_Hx2y2z_S_G2y2z_S_vrr;
      I_ERI_Hxy3z_S_G2y2z_S += I_ERI_Hxy3z_S_G2y2z_S_vrr;
      I_ERI_Hx4z_S_G2y2z_S += I_ERI_Hx4z_S_G2y2z_S_vrr;
      I_ERI_H5y_S_G2y2z_S += I_ERI_H5y_S_G2y2z_S_vrr;
      I_ERI_H4yz_S_G2y2z_S += I_ERI_H4yz_S_G2y2z_S_vrr;
      I_ERI_H3y2z_S_G2y2z_S += I_ERI_H3y2z_S_G2y2z_S_vrr;
      I_ERI_H2y3z_S_G2y2z_S += I_ERI_H2y3z_S_G2y2z_S_vrr;
      I_ERI_Hy4z_S_G2y2z_S += I_ERI_Hy4z_S_G2y2z_S_vrr;
      I_ERI_H5z_S_G2y2z_S += I_ERI_H5z_S_G2y2z_S_vrr;
      I_ERI_H5x_S_Gy3z_S += I_ERI_H5x_S_Gy3z_S_vrr;
      I_ERI_H4xy_S_Gy3z_S += I_ERI_H4xy_S_Gy3z_S_vrr;
      I_ERI_H4xz_S_Gy3z_S += I_ERI_H4xz_S_Gy3z_S_vrr;
      I_ERI_H3x2y_S_Gy3z_S += I_ERI_H3x2y_S_Gy3z_S_vrr;
      I_ERI_H3xyz_S_Gy3z_S += I_ERI_H3xyz_S_Gy3z_S_vrr;
      I_ERI_H3x2z_S_Gy3z_S += I_ERI_H3x2z_S_Gy3z_S_vrr;
      I_ERI_H2x3y_S_Gy3z_S += I_ERI_H2x3y_S_Gy3z_S_vrr;
      I_ERI_H2x2yz_S_Gy3z_S += I_ERI_H2x2yz_S_Gy3z_S_vrr;
      I_ERI_H2xy2z_S_Gy3z_S += I_ERI_H2xy2z_S_Gy3z_S_vrr;
      I_ERI_H2x3z_S_Gy3z_S += I_ERI_H2x3z_S_Gy3z_S_vrr;
      I_ERI_Hx4y_S_Gy3z_S += I_ERI_Hx4y_S_Gy3z_S_vrr;
      I_ERI_Hx3yz_S_Gy3z_S += I_ERI_Hx3yz_S_Gy3z_S_vrr;
      I_ERI_Hx2y2z_S_Gy3z_S += I_ERI_Hx2y2z_S_Gy3z_S_vrr;
      I_ERI_Hxy3z_S_Gy3z_S += I_ERI_Hxy3z_S_Gy3z_S_vrr;
      I_ERI_Hx4z_S_Gy3z_S += I_ERI_Hx4z_S_Gy3z_S_vrr;
      I_ERI_H5y_S_Gy3z_S += I_ERI_H5y_S_Gy3z_S_vrr;
      I_ERI_H4yz_S_Gy3z_S += I_ERI_H4yz_S_Gy3z_S_vrr;
      I_ERI_H3y2z_S_Gy3z_S += I_ERI_H3y2z_S_Gy3z_S_vrr;
      I_ERI_H2y3z_S_Gy3z_S += I_ERI_H2y3z_S_Gy3z_S_vrr;
      I_ERI_Hy4z_S_Gy3z_S += I_ERI_Hy4z_S_Gy3z_S_vrr;
      I_ERI_H5z_S_Gy3z_S += I_ERI_H5z_S_Gy3z_S_vrr;
      I_ERI_H5x_S_G4z_S += I_ERI_H5x_S_G4z_S_vrr;
      I_ERI_H4xy_S_G4z_S += I_ERI_H4xy_S_G4z_S_vrr;
      I_ERI_H4xz_S_G4z_S += I_ERI_H4xz_S_G4z_S_vrr;
      I_ERI_H3x2y_S_G4z_S += I_ERI_H3x2y_S_G4z_S_vrr;
      I_ERI_H3xyz_S_G4z_S += I_ERI_H3xyz_S_G4z_S_vrr;
      I_ERI_H3x2z_S_G4z_S += I_ERI_H3x2z_S_G4z_S_vrr;
      I_ERI_H2x3y_S_G4z_S += I_ERI_H2x3y_S_G4z_S_vrr;
      I_ERI_H2x2yz_S_G4z_S += I_ERI_H2x2yz_S_G4z_S_vrr;
      I_ERI_H2xy2z_S_G4z_S += I_ERI_H2xy2z_S_G4z_S_vrr;
      I_ERI_H2x3z_S_G4z_S += I_ERI_H2x3z_S_G4z_S_vrr;
      I_ERI_Hx4y_S_G4z_S += I_ERI_Hx4y_S_G4z_S_vrr;
      I_ERI_Hx3yz_S_G4z_S += I_ERI_Hx3yz_S_G4z_S_vrr;
      I_ERI_Hx2y2z_S_G4z_S += I_ERI_Hx2y2z_S_G4z_S_vrr;
      I_ERI_Hxy3z_S_G4z_S += I_ERI_Hxy3z_S_G4z_S_vrr;
      I_ERI_Hx4z_S_G4z_S += I_ERI_Hx4z_S_G4z_S_vrr;
      I_ERI_H5y_S_G4z_S += I_ERI_H5y_S_G4z_S_vrr;
      I_ERI_H4yz_S_G4z_S += I_ERI_H4yz_S_G4z_S_vrr;
      I_ERI_H3y2z_S_G4z_S += I_ERI_H3y2z_S_G4z_S_vrr;
      I_ERI_H2y3z_S_G4z_S += I_ERI_H2y3z_S_G4z_S_vrr;
      I_ERI_Hy4z_S_G4z_S += I_ERI_Hy4z_S_G4z_S_vrr;
      I_ERI_H5z_S_G4z_S += I_ERI_H5z_S_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_G_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_G4x_S_G4x_S += I_ERI_G4x_S_G4x_S_vrr;
      I_ERI_G3xy_S_G4x_S += I_ERI_G3xy_S_G4x_S_vrr;
      I_ERI_G3xz_S_G4x_S += I_ERI_G3xz_S_G4x_S_vrr;
      I_ERI_G2x2y_S_G4x_S += I_ERI_G2x2y_S_G4x_S_vrr;
      I_ERI_G2xyz_S_G4x_S += I_ERI_G2xyz_S_G4x_S_vrr;
      I_ERI_G2x2z_S_G4x_S += I_ERI_G2x2z_S_G4x_S_vrr;
      I_ERI_Gx3y_S_G4x_S += I_ERI_Gx3y_S_G4x_S_vrr;
      I_ERI_Gx2yz_S_G4x_S += I_ERI_Gx2yz_S_G4x_S_vrr;
      I_ERI_Gxy2z_S_G4x_S += I_ERI_Gxy2z_S_G4x_S_vrr;
      I_ERI_Gx3z_S_G4x_S += I_ERI_Gx3z_S_G4x_S_vrr;
      I_ERI_G4y_S_G4x_S += I_ERI_G4y_S_G4x_S_vrr;
      I_ERI_G3yz_S_G4x_S += I_ERI_G3yz_S_G4x_S_vrr;
      I_ERI_G2y2z_S_G4x_S += I_ERI_G2y2z_S_G4x_S_vrr;
      I_ERI_Gy3z_S_G4x_S += I_ERI_Gy3z_S_G4x_S_vrr;
      I_ERI_G4z_S_G4x_S += I_ERI_G4z_S_G4x_S_vrr;
      I_ERI_G4x_S_G3xy_S += I_ERI_G4x_S_G3xy_S_vrr;
      I_ERI_G3xy_S_G3xy_S += I_ERI_G3xy_S_G3xy_S_vrr;
      I_ERI_G3xz_S_G3xy_S += I_ERI_G3xz_S_G3xy_S_vrr;
      I_ERI_G2x2y_S_G3xy_S += I_ERI_G2x2y_S_G3xy_S_vrr;
      I_ERI_G2xyz_S_G3xy_S += I_ERI_G2xyz_S_G3xy_S_vrr;
      I_ERI_G2x2z_S_G3xy_S += I_ERI_G2x2z_S_G3xy_S_vrr;
      I_ERI_Gx3y_S_G3xy_S += I_ERI_Gx3y_S_G3xy_S_vrr;
      I_ERI_Gx2yz_S_G3xy_S += I_ERI_Gx2yz_S_G3xy_S_vrr;
      I_ERI_Gxy2z_S_G3xy_S += I_ERI_Gxy2z_S_G3xy_S_vrr;
      I_ERI_Gx3z_S_G3xy_S += I_ERI_Gx3z_S_G3xy_S_vrr;
      I_ERI_G4y_S_G3xy_S += I_ERI_G4y_S_G3xy_S_vrr;
      I_ERI_G3yz_S_G3xy_S += I_ERI_G3yz_S_G3xy_S_vrr;
      I_ERI_G2y2z_S_G3xy_S += I_ERI_G2y2z_S_G3xy_S_vrr;
      I_ERI_Gy3z_S_G3xy_S += I_ERI_Gy3z_S_G3xy_S_vrr;
      I_ERI_G4z_S_G3xy_S += I_ERI_G4z_S_G3xy_S_vrr;
      I_ERI_G4x_S_G3xz_S += I_ERI_G4x_S_G3xz_S_vrr;
      I_ERI_G3xy_S_G3xz_S += I_ERI_G3xy_S_G3xz_S_vrr;
      I_ERI_G3xz_S_G3xz_S += I_ERI_G3xz_S_G3xz_S_vrr;
      I_ERI_G2x2y_S_G3xz_S += I_ERI_G2x2y_S_G3xz_S_vrr;
      I_ERI_G2xyz_S_G3xz_S += I_ERI_G2xyz_S_G3xz_S_vrr;
      I_ERI_G2x2z_S_G3xz_S += I_ERI_G2x2z_S_G3xz_S_vrr;
      I_ERI_Gx3y_S_G3xz_S += I_ERI_Gx3y_S_G3xz_S_vrr;
      I_ERI_Gx2yz_S_G3xz_S += I_ERI_Gx2yz_S_G3xz_S_vrr;
      I_ERI_Gxy2z_S_G3xz_S += I_ERI_Gxy2z_S_G3xz_S_vrr;
      I_ERI_Gx3z_S_G3xz_S += I_ERI_Gx3z_S_G3xz_S_vrr;
      I_ERI_G4y_S_G3xz_S += I_ERI_G4y_S_G3xz_S_vrr;
      I_ERI_G3yz_S_G3xz_S += I_ERI_G3yz_S_G3xz_S_vrr;
      I_ERI_G2y2z_S_G3xz_S += I_ERI_G2y2z_S_G3xz_S_vrr;
      I_ERI_Gy3z_S_G3xz_S += I_ERI_Gy3z_S_G3xz_S_vrr;
      I_ERI_G4z_S_G3xz_S += I_ERI_G4z_S_G3xz_S_vrr;
      I_ERI_G4x_S_G2x2y_S += I_ERI_G4x_S_G2x2y_S_vrr;
      I_ERI_G3xy_S_G2x2y_S += I_ERI_G3xy_S_G2x2y_S_vrr;
      I_ERI_G3xz_S_G2x2y_S += I_ERI_G3xz_S_G2x2y_S_vrr;
      I_ERI_G2x2y_S_G2x2y_S += I_ERI_G2x2y_S_G2x2y_S_vrr;
      I_ERI_G2xyz_S_G2x2y_S += I_ERI_G2xyz_S_G2x2y_S_vrr;
      I_ERI_G2x2z_S_G2x2y_S += I_ERI_G2x2z_S_G2x2y_S_vrr;
      I_ERI_Gx3y_S_G2x2y_S += I_ERI_Gx3y_S_G2x2y_S_vrr;
      I_ERI_Gx2yz_S_G2x2y_S += I_ERI_Gx2yz_S_G2x2y_S_vrr;
      I_ERI_Gxy2z_S_G2x2y_S += I_ERI_Gxy2z_S_G2x2y_S_vrr;
      I_ERI_Gx3z_S_G2x2y_S += I_ERI_Gx3z_S_G2x2y_S_vrr;
      I_ERI_G4y_S_G2x2y_S += I_ERI_G4y_S_G2x2y_S_vrr;
      I_ERI_G3yz_S_G2x2y_S += I_ERI_G3yz_S_G2x2y_S_vrr;
      I_ERI_G2y2z_S_G2x2y_S += I_ERI_G2y2z_S_G2x2y_S_vrr;
      I_ERI_Gy3z_S_G2x2y_S += I_ERI_Gy3z_S_G2x2y_S_vrr;
      I_ERI_G4z_S_G2x2y_S += I_ERI_G4z_S_G2x2y_S_vrr;
      I_ERI_G4x_S_G2xyz_S += I_ERI_G4x_S_G2xyz_S_vrr;
      I_ERI_G3xy_S_G2xyz_S += I_ERI_G3xy_S_G2xyz_S_vrr;
      I_ERI_G3xz_S_G2xyz_S += I_ERI_G3xz_S_G2xyz_S_vrr;
      I_ERI_G2x2y_S_G2xyz_S += I_ERI_G2x2y_S_G2xyz_S_vrr;
      I_ERI_G2xyz_S_G2xyz_S += I_ERI_G2xyz_S_G2xyz_S_vrr;
      I_ERI_G2x2z_S_G2xyz_S += I_ERI_G2x2z_S_G2xyz_S_vrr;
      I_ERI_Gx3y_S_G2xyz_S += I_ERI_Gx3y_S_G2xyz_S_vrr;
      I_ERI_Gx2yz_S_G2xyz_S += I_ERI_Gx2yz_S_G2xyz_S_vrr;
      I_ERI_Gxy2z_S_G2xyz_S += I_ERI_Gxy2z_S_G2xyz_S_vrr;
      I_ERI_Gx3z_S_G2xyz_S += I_ERI_Gx3z_S_G2xyz_S_vrr;
      I_ERI_G4y_S_G2xyz_S += I_ERI_G4y_S_G2xyz_S_vrr;
      I_ERI_G3yz_S_G2xyz_S += I_ERI_G3yz_S_G2xyz_S_vrr;
      I_ERI_G2y2z_S_G2xyz_S += I_ERI_G2y2z_S_G2xyz_S_vrr;
      I_ERI_Gy3z_S_G2xyz_S += I_ERI_Gy3z_S_G2xyz_S_vrr;
      I_ERI_G4z_S_G2xyz_S += I_ERI_G4z_S_G2xyz_S_vrr;
      I_ERI_G4x_S_G2x2z_S += I_ERI_G4x_S_G2x2z_S_vrr;
      I_ERI_G3xy_S_G2x2z_S += I_ERI_G3xy_S_G2x2z_S_vrr;
      I_ERI_G3xz_S_G2x2z_S += I_ERI_G3xz_S_G2x2z_S_vrr;
      I_ERI_G2x2y_S_G2x2z_S += I_ERI_G2x2y_S_G2x2z_S_vrr;
      I_ERI_G2xyz_S_G2x2z_S += I_ERI_G2xyz_S_G2x2z_S_vrr;
      I_ERI_G2x2z_S_G2x2z_S += I_ERI_G2x2z_S_G2x2z_S_vrr;
      I_ERI_Gx3y_S_G2x2z_S += I_ERI_Gx3y_S_G2x2z_S_vrr;
      I_ERI_Gx2yz_S_G2x2z_S += I_ERI_Gx2yz_S_G2x2z_S_vrr;
      I_ERI_Gxy2z_S_G2x2z_S += I_ERI_Gxy2z_S_G2x2z_S_vrr;
      I_ERI_Gx3z_S_G2x2z_S += I_ERI_Gx3z_S_G2x2z_S_vrr;
      I_ERI_G4y_S_G2x2z_S += I_ERI_G4y_S_G2x2z_S_vrr;
      I_ERI_G3yz_S_G2x2z_S += I_ERI_G3yz_S_G2x2z_S_vrr;
      I_ERI_G2y2z_S_G2x2z_S += I_ERI_G2y2z_S_G2x2z_S_vrr;
      I_ERI_Gy3z_S_G2x2z_S += I_ERI_Gy3z_S_G2x2z_S_vrr;
      I_ERI_G4z_S_G2x2z_S += I_ERI_G4z_S_G2x2z_S_vrr;
      I_ERI_G4x_S_Gx3y_S += I_ERI_G4x_S_Gx3y_S_vrr;
      I_ERI_G3xy_S_Gx3y_S += I_ERI_G3xy_S_Gx3y_S_vrr;
      I_ERI_G3xz_S_Gx3y_S += I_ERI_G3xz_S_Gx3y_S_vrr;
      I_ERI_G2x2y_S_Gx3y_S += I_ERI_G2x2y_S_Gx3y_S_vrr;
      I_ERI_G2xyz_S_Gx3y_S += I_ERI_G2xyz_S_Gx3y_S_vrr;
      I_ERI_G2x2z_S_Gx3y_S += I_ERI_G2x2z_S_Gx3y_S_vrr;
      I_ERI_Gx3y_S_Gx3y_S += I_ERI_Gx3y_S_Gx3y_S_vrr;
      I_ERI_Gx2yz_S_Gx3y_S += I_ERI_Gx2yz_S_Gx3y_S_vrr;
      I_ERI_Gxy2z_S_Gx3y_S += I_ERI_Gxy2z_S_Gx3y_S_vrr;
      I_ERI_Gx3z_S_Gx3y_S += I_ERI_Gx3z_S_Gx3y_S_vrr;
      I_ERI_G4y_S_Gx3y_S += I_ERI_G4y_S_Gx3y_S_vrr;
      I_ERI_G3yz_S_Gx3y_S += I_ERI_G3yz_S_Gx3y_S_vrr;
      I_ERI_G2y2z_S_Gx3y_S += I_ERI_G2y2z_S_Gx3y_S_vrr;
      I_ERI_Gy3z_S_Gx3y_S += I_ERI_Gy3z_S_Gx3y_S_vrr;
      I_ERI_G4z_S_Gx3y_S += I_ERI_G4z_S_Gx3y_S_vrr;
      I_ERI_G4x_S_Gx2yz_S += I_ERI_G4x_S_Gx2yz_S_vrr;
      I_ERI_G3xy_S_Gx2yz_S += I_ERI_G3xy_S_Gx2yz_S_vrr;
      I_ERI_G3xz_S_Gx2yz_S += I_ERI_G3xz_S_Gx2yz_S_vrr;
      I_ERI_G2x2y_S_Gx2yz_S += I_ERI_G2x2y_S_Gx2yz_S_vrr;
      I_ERI_G2xyz_S_Gx2yz_S += I_ERI_G2xyz_S_Gx2yz_S_vrr;
      I_ERI_G2x2z_S_Gx2yz_S += I_ERI_G2x2z_S_Gx2yz_S_vrr;
      I_ERI_Gx3y_S_Gx2yz_S += I_ERI_Gx3y_S_Gx2yz_S_vrr;
      I_ERI_Gx2yz_S_Gx2yz_S += I_ERI_Gx2yz_S_Gx2yz_S_vrr;
      I_ERI_Gxy2z_S_Gx2yz_S += I_ERI_Gxy2z_S_Gx2yz_S_vrr;
      I_ERI_Gx3z_S_Gx2yz_S += I_ERI_Gx3z_S_Gx2yz_S_vrr;
      I_ERI_G4y_S_Gx2yz_S += I_ERI_G4y_S_Gx2yz_S_vrr;
      I_ERI_G3yz_S_Gx2yz_S += I_ERI_G3yz_S_Gx2yz_S_vrr;
      I_ERI_G2y2z_S_Gx2yz_S += I_ERI_G2y2z_S_Gx2yz_S_vrr;
      I_ERI_Gy3z_S_Gx2yz_S += I_ERI_Gy3z_S_Gx2yz_S_vrr;
      I_ERI_G4z_S_Gx2yz_S += I_ERI_G4z_S_Gx2yz_S_vrr;
      I_ERI_G4x_S_Gxy2z_S += I_ERI_G4x_S_Gxy2z_S_vrr;
      I_ERI_G3xy_S_Gxy2z_S += I_ERI_G3xy_S_Gxy2z_S_vrr;
      I_ERI_G3xz_S_Gxy2z_S += I_ERI_G3xz_S_Gxy2z_S_vrr;
      I_ERI_G2x2y_S_Gxy2z_S += I_ERI_G2x2y_S_Gxy2z_S_vrr;
      I_ERI_G2xyz_S_Gxy2z_S += I_ERI_G2xyz_S_Gxy2z_S_vrr;
      I_ERI_G2x2z_S_Gxy2z_S += I_ERI_G2x2z_S_Gxy2z_S_vrr;
      I_ERI_Gx3y_S_Gxy2z_S += I_ERI_Gx3y_S_Gxy2z_S_vrr;
      I_ERI_Gx2yz_S_Gxy2z_S += I_ERI_Gx2yz_S_Gxy2z_S_vrr;
      I_ERI_Gxy2z_S_Gxy2z_S += I_ERI_Gxy2z_S_Gxy2z_S_vrr;
      I_ERI_Gx3z_S_Gxy2z_S += I_ERI_Gx3z_S_Gxy2z_S_vrr;
      I_ERI_G4y_S_Gxy2z_S += I_ERI_G4y_S_Gxy2z_S_vrr;
      I_ERI_G3yz_S_Gxy2z_S += I_ERI_G3yz_S_Gxy2z_S_vrr;
      I_ERI_G2y2z_S_Gxy2z_S += I_ERI_G2y2z_S_Gxy2z_S_vrr;
      I_ERI_Gy3z_S_Gxy2z_S += I_ERI_Gy3z_S_Gxy2z_S_vrr;
      I_ERI_G4z_S_Gxy2z_S += I_ERI_G4z_S_Gxy2z_S_vrr;
      I_ERI_G4x_S_Gx3z_S += I_ERI_G4x_S_Gx3z_S_vrr;
      I_ERI_G3xy_S_Gx3z_S += I_ERI_G3xy_S_Gx3z_S_vrr;
      I_ERI_G3xz_S_Gx3z_S += I_ERI_G3xz_S_Gx3z_S_vrr;
      I_ERI_G2x2y_S_Gx3z_S += I_ERI_G2x2y_S_Gx3z_S_vrr;
      I_ERI_G2xyz_S_Gx3z_S += I_ERI_G2xyz_S_Gx3z_S_vrr;
      I_ERI_G2x2z_S_Gx3z_S += I_ERI_G2x2z_S_Gx3z_S_vrr;
      I_ERI_Gx3y_S_Gx3z_S += I_ERI_Gx3y_S_Gx3z_S_vrr;
      I_ERI_Gx2yz_S_Gx3z_S += I_ERI_Gx2yz_S_Gx3z_S_vrr;
      I_ERI_Gxy2z_S_Gx3z_S += I_ERI_Gxy2z_S_Gx3z_S_vrr;
      I_ERI_Gx3z_S_Gx3z_S += I_ERI_Gx3z_S_Gx3z_S_vrr;
      I_ERI_G4y_S_Gx3z_S += I_ERI_G4y_S_Gx3z_S_vrr;
      I_ERI_G3yz_S_Gx3z_S += I_ERI_G3yz_S_Gx3z_S_vrr;
      I_ERI_G2y2z_S_Gx3z_S += I_ERI_G2y2z_S_Gx3z_S_vrr;
      I_ERI_Gy3z_S_Gx3z_S += I_ERI_Gy3z_S_Gx3z_S_vrr;
      I_ERI_G4z_S_Gx3z_S += I_ERI_G4z_S_Gx3z_S_vrr;
      I_ERI_G4x_S_G4y_S += I_ERI_G4x_S_G4y_S_vrr;
      I_ERI_G3xy_S_G4y_S += I_ERI_G3xy_S_G4y_S_vrr;
      I_ERI_G3xz_S_G4y_S += I_ERI_G3xz_S_G4y_S_vrr;
      I_ERI_G2x2y_S_G4y_S += I_ERI_G2x2y_S_G4y_S_vrr;
      I_ERI_G2xyz_S_G4y_S += I_ERI_G2xyz_S_G4y_S_vrr;
      I_ERI_G2x2z_S_G4y_S += I_ERI_G2x2z_S_G4y_S_vrr;
      I_ERI_Gx3y_S_G4y_S += I_ERI_Gx3y_S_G4y_S_vrr;
      I_ERI_Gx2yz_S_G4y_S += I_ERI_Gx2yz_S_G4y_S_vrr;
      I_ERI_Gxy2z_S_G4y_S += I_ERI_Gxy2z_S_G4y_S_vrr;
      I_ERI_Gx3z_S_G4y_S += I_ERI_Gx3z_S_G4y_S_vrr;
      I_ERI_G4y_S_G4y_S += I_ERI_G4y_S_G4y_S_vrr;
      I_ERI_G3yz_S_G4y_S += I_ERI_G3yz_S_G4y_S_vrr;
      I_ERI_G2y2z_S_G4y_S += I_ERI_G2y2z_S_G4y_S_vrr;
      I_ERI_Gy3z_S_G4y_S += I_ERI_Gy3z_S_G4y_S_vrr;
      I_ERI_G4z_S_G4y_S += I_ERI_G4z_S_G4y_S_vrr;
      I_ERI_G4x_S_G3yz_S += I_ERI_G4x_S_G3yz_S_vrr;
      I_ERI_G3xy_S_G3yz_S += I_ERI_G3xy_S_G3yz_S_vrr;
      I_ERI_G3xz_S_G3yz_S += I_ERI_G3xz_S_G3yz_S_vrr;
      I_ERI_G2x2y_S_G3yz_S += I_ERI_G2x2y_S_G3yz_S_vrr;
      I_ERI_G2xyz_S_G3yz_S += I_ERI_G2xyz_S_G3yz_S_vrr;
      I_ERI_G2x2z_S_G3yz_S += I_ERI_G2x2z_S_G3yz_S_vrr;
      I_ERI_Gx3y_S_G3yz_S += I_ERI_Gx3y_S_G3yz_S_vrr;
      I_ERI_Gx2yz_S_G3yz_S += I_ERI_Gx2yz_S_G3yz_S_vrr;
      I_ERI_Gxy2z_S_G3yz_S += I_ERI_Gxy2z_S_G3yz_S_vrr;
      I_ERI_Gx3z_S_G3yz_S += I_ERI_Gx3z_S_G3yz_S_vrr;
      I_ERI_G4y_S_G3yz_S += I_ERI_G4y_S_G3yz_S_vrr;
      I_ERI_G3yz_S_G3yz_S += I_ERI_G3yz_S_G3yz_S_vrr;
      I_ERI_G2y2z_S_G3yz_S += I_ERI_G2y2z_S_G3yz_S_vrr;
      I_ERI_Gy3z_S_G3yz_S += I_ERI_Gy3z_S_G3yz_S_vrr;
      I_ERI_G4z_S_G3yz_S += I_ERI_G4z_S_G3yz_S_vrr;
      I_ERI_G4x_S_G2y2z_S += I_ERI_G4x_S_G2y2z_S_vrr;
      I_ERI_G3xy_S_G2y2z_S += I_ERI_G3xy_S_G2y2z_S_vrr;
      I_ERI_G3xz_S_G2y2z_S += I_ERI_G3xz_S_G2y2z_S_vrr;
      I_ERI_G2x2y_S_G2y2z_S += I_ERI_G2x2y_S_G2y2z_S_vrr;
      I_ERI_G2xyz_S_G2y2z_S += I_ERI_G2xyz_S_G2y2z_S_vrr;
      I_ERI_G2x2z_S_G2y2z_S += I_ERI_G2x2z_S_G2y2z_S_vrr;
      I_ERI_Gx3y_S_G2y2z_S += I_ERI_Gx3y_S_G2y2z_S_vrr;
      I_ERI_Gx2yz_S_G2y2z_S += I_ERI_Gx2yz_S_G2y2z_S_vrr;
      I_ERI_Gxy2z_S_G2y2z_S += I_ERI_Gxy2z_S_G2y2z_S_vrr;
      I_ERI_Gx3z_S_G2y2z_S += I_ERI_Gx3z_S_G2y2z_S_vrr;
      I_ERI_G4y_S_G2y2z_S += I_ERI_G4y_S_G2y2z_S_vrr;
      I_ERI_G3yz_S_G2y2z_S += I_ERI_G3yz_S_G2y2z_S_vrr;
      I_ERI_G2y2z_S_G2y2z_S += I_ERI_G2y2z_S_G2y2z_S_vrr;
      I_ERI_Gy3z_S_G2y2z_S += I_ERI_Gy3z_S_G2y2z_S_vrr;
      I_ERI_G4z_S_G2y2z_S += I_ERI_G4z_S_G2y2z_S_vrr;
      I_ERI_G4x_S_Gy3z_S += I_ERI_G4x_S_Gy3z_S_vrr;
      I_ERI_G3xy_S_Gy3z_S += I_ERI_G3xy_S_Gy3z_S_vrr;
      I_ERI_G3xz_S_Gy3z_S += I_ERI_G3xz_S_Gy3z_S_vrr;
      I_ERI_G2x2y_S_Gy3z_S += I_ERI_G2x2y_S_Gy3z_S_vrr;
      I_ERI_G2xyz_S_Gy3z_S += I_ERI_G2xyz_S_Gy3z_S_vrr;
      I_ERI_G2x2z_S_Gy3z_S += I_ERI_G2x2z_S_Gy3z_S_vrr;
      I_ERI_Gx3y_S_Gy3z_S += I_ERI_Gx3y_S_Gy3z_S_vrr;
      I_ERI_Gx2yz_S_Gy3z_S += I_ERI_Gx2yz_S_Gy3z_S_vrr;
      I_ERI_Gxy2z_S_Gy3z_S += I_ERI_Gxy2z_S_Gy3z_S_vrr;
      I_ERI_Gx3z_S_Gy3z_S += I_ERI_Gx3z_S_Gy3z_S_vrr;
      I_ERI_G4y_S_Gy3z_S += I_ERI_G4y_S_Gy3z_S_vrr;
      I_ERI_G3yz_S_Gy3z_S += I_ERI_G3yz_S_Gy3z_S_vrr;
      I_ERI_G2y2z_S_Gy3z_S += I_ERI_G2y2z_S_Gy3z_S_vrr;
      I_ERI_Gy3z_S_Gy3z_S += I_ERI_Gy3z_S_Gy3z_S_vrr;
      I_ERI_G4z_S_Gy3z_S += I_ERI_G4z_S_Gy3z_S_vrr;
      I_ERI_G4x_S_G4z_S += I_ERI_G4x_S_G4z_S_vrr;
      I_ERI_G3xy_S_G4z_S += I_ERI_G3xy_S_G4z_S_vrr;
      I_ERI_G3xz_S_G4z_S += I_ERI_G3xz_S_G4z_S_vrr;
      I_ERI_G2x2y_S_G4z_S += I_ERI_G2x2y_S_G4z_S_vrr;
      I_ERI_G2xyz_S_G4z_S += I_ERI_G2xyz_S_G4z_S_vrr;
      I_ERI_G2x2z_S_G4z_S += I_ERI_G2x2z_S_G4z_S_vrr;
      I_ERI_Gx3y_S_G4z_S += I_ERI_Gx3y_S_G4z_S_vrr;
      I_ERI_Gx2yz_S_G4z_S += I_ERI_Gx2yz_S_G4z_S_vrr;
      I_ERI_Gxy2z_S_G4z_S += I_ERI_Gxy2z_S_G4z_S_vrr;
      I_ERI_Gx3z_S_G4z_S += I_ERI_Gx3z_S_G4z_S_vrr;
      I_ERI_G4y_S_G4z_S += I_ERI_G4y_S_G4z_S_vrr;
      I_ERI_G3yz_S_G4z_S += I_ERI_G3yz_S_G4z_S_vrr;
      I_ERI_G2y2z_S_G4z_S += I_ERI_G2y2z_S_G4z_S_vrr;
      I_ERI_Gy3z_S_G4z_S += I_ERI_Gy3z_S_G4z_S_vrr;
      I_ERI_G4z_S_G4z_S += I_ERI_G4z_S_G4z_S_vrr;
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
   * shell quartet name: SQ_ERI_G_P_G_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_G_S
   * RHS shell quartet name: SQ_ERI_G_S_G_S
   ************************************************************/
  Double I_ERI_G4x_Px_G4x_S = I_ERI_H5x_S_G4x_S+ABX*I_ERI_G4x_S_G4x_S;
  Double I_ERI_G3xy_Px_G4x_S = I_ERI_H4xy_S_G4x_S+ABX*I_ERI_G3xy_S_G4x_S;
  Double I_ERI_G3xz_Px_G4x_S = I_ERI_H4xz_S_G4x_S+ABX*I_ERI_G3xz_S_G4x_S;
  Double I_ERI_G2x2y_Px_G4x_S = I_ERI_H3x2y_S_G4x_S+ABX*I_ERI_G2x2y_S_G4x_S;
  Double I_ERI_G2xyz_Px_G4x_S = I_ERI_H3xyz_S_G4x_S+ABX*I_ERI_G2xyz_S_G4x_S;
  Double I_ERI_G2x2z_Px_G4x_S = I_ERI_H3x2z_S_G4x_S+ABX*I_ERI_G2x2z_S_G4x_S;
  Double I_ERI_Gx3y_Px_G4x_S = I_ERI_H2x3y_S_G4x_S+ABX*I_ERI_Gx3y_S_G4x_S;
  Double I_ERI_Gx2yz_Px_G4x_S = I_ERI_H2x2yz_S_G4x_S+ABX*I_ERI_Gx2yz_S_G4x_S;
  Double I_ERI_Gxy2z_Px_G4x_S = I_ERI_H2xy2z_S_G4x_S+ABX*I_ERI_Gxy2z_S_G4x_S;
  Double I_ERI_Gx3z_Px_G4x_S = I_ERI_H2x3z_S_G4x_S+ABX*I_ERI_Gx3z_S_G4x_S;
  Double I_ERI_G4y_Px_G4x_S = I_ERI_Hx4y_S_G4x_S+ABX*I_ERI_G4y_S_G4x_S;
  Double I_ERI_G3yz_Px_G4x_S = I_ERI_Hx3yz_S_G4x_S+ABX*I_ERI_G3yz_S_G4x_S;
  Double I_ERI_G2y2z_Px_G4x_S = I_ERI_Hx2y2z_S_G4x_S+ABX*I_ERI_G2y2z_S_G4x_S;
  Double I_ERI_Gy3z_Px_G4x_S = I_ERI_Hxy3z_S_G4x_S+ABX*I_ERI_Gy3z_S_G4x_S;
  Double I_ERI_G4z_Px_G4x_S = I_ERI_Hx4z_S_G4x_S+ABX*I_ERI_G4z_S_G4x_S;
  Double I_ERI_G4x_Py_G4x_S = I_ERI_H4xy_S_G4x_S+ABY*I_ERI_G4x_S_G4x_S;
  Double I_ERI_G3xy_Py_G4x_S = I_ERI_H3x2y_S_G4x_S+ABY*I_ERI_G3xy_S_G4x_S;
  Double I_ERI_G3xz_Py_G4x_S = I_ERI_H3xyz_S_G4x_S+ABY*I_ERI_G3xz_S_G4x_S;
  Double I_ERI_G2x2y_Py_G4x_S = I_ERI_H2x3y_S_G4x_S+ABY*I_ERI_G2x2y_S_G4x_S;
  Double I_ERI_G2xyz_Py_G4x_S = I_ERI_H2x2yz_S_G4x_S+ABY*I_ERI_G2xyz_S_G4x_S;
  Double I_ERI_G2x2z_Py_G4x_S = I_ERI_H2xy2z_S_G4x_S+ABY*I_ERI_G2x2z_S_G4x_S;
  Double I_ERI_Gx3y_Py_G4x_S = I_ERI_Hx4y_S_G4x_S+ABY*I_ERI_Gx3y_S_G4x_S;
  Double I_ERI_Gx2yz_Py_G4x_S = I_ERI_Hx3yz_S_G4x_S+ABY*I_ERI_Gx2yz_S_G4x_S;
  Double I_ERI_Gxy2z_Py_G4x_S = I_ERI_Hx2y2z_S_G4x_S+ABY*I_ERI_Gxy2z_S_G4x_S;
  Double I_ERI_Gx3z_Py_G4x_S = I_ERI_Hxy3z_S_G4x_S+ABY*I_ERI_Gx3z_S_G4x_S;
  Double I_ERI_G4y_Py_G4x_S = I_ERI_H5y_S_G4x_S+ABY*I_ERI_G4y_S_G4x_S;
  Double I_ERI_G3yz_Py_G4x_S = I_ERI_H4yz_S_G4x_S+ABY*I_ERI_G3yz_S_G4x_S;
  Double I_ERI_G2y2z_Py_G4x_S = I_ERI_H3y2z_S_G4x_S+ABY*I_ERI_G2y2z_S_G4x_S;
  Double I_ERI_Gy3z_Py_G4x_S = I_ERI_H2y3z_S_G4x_S+ABY*I_ERI_Gy3z_S_G4x_S;
  Double I_ERI_G4z_Py_G4x_S = I_ERI_Hy4z_S_G4x_S+ABY*I_ERI_G4z_S_G4x_S;
  Double I_ERI_G4x_Pz_G4x_S = I_ERI_H4xz_S_G4x_S+ABZ*I_ERI_G4x_S_G4x_S;
  Double I_ERI_G3xy_Pz_G4x_S = I_ERI_H3xyz_S_G4x_S+ABZ*I_ERI_G3xy_S_G4x_S;
  Double I_ERI_G3xz_Pz_G4x_S = I_ERI_H3x2z_S_G4x_S+ABZ*I_ERI_G3xz_S_G4x_S;
  Double I_ERI_G2x2y_Pz_G4x_S = I_ERI_H2x2yz_S_G4x_S+ABZ*I_ERI_G2x2y_S_G4x_S;
  Double I_ERI_G2xyz_Pz_G4x_S = I_ERI_H2xy2z_S_G4x_S+ABZ*I_ERI_G2xyz_S_G4x_S;
  Double I_ERI_G2x2z_Pz_G4x_S = I_ERI_H2x3z_S_G4x_S+ABZ*I_ERI_G2x2z_S_G4x_S;
  Double I_ERI_Gx3y_Pz_G4x_S = I_ERI_Hx3yz_S_G4x_S+ABZ*I_ERI_Gx3y_S_G4x_S;
  Double I_ERI_Gx2yz_Pz_G4x_S = I_ERI_Hx2y2z_S_G4x_S+ABZ*I_ERI_Gx2yz_S_G4x_S;
  Double I_ERI_Gxy2z_Pz_G4x_S = I_ERI_Hxy3z_S_G4x_S+ABZ*I_ERI_Gxy2z_S_G4x_S;
  Double I_ERI_Gx3z_Pz_G4x_S = I_ERI_Hx4z_S_G4x_S+ABZ*I_ERI_Gx3z_S_G4x_S;
  Double I_ERI_G4y_Pz_G4x_S = I_ERI_H4yz_S_G4x_S+ABZ*I_ERI_G4y_S_G4x_S;
  Double I_ERI_G3yz_Pz_G4x_S = I_ERI_H3y2z_S_G4x_S+ABZ*I_ERI_G3yz_S_G4x_S;
  Double I_ERI_G2y2z_Pz_G4x_S = I_ERI_H2y3z_S_G4x_S+ABZ*I_ERI_G2y2z_S_G4x_S;
  Double I_ERI_Gy3z_Pz_G4x_S = I_ERI_Hy4z_S_G4x_S+ABZ*I_ERI_Gy3z_S_G4x_S;
  Double I_ERI_G4z_Pz_G4x_S = I_ERI_H5z_S_G4x_S+ABZ*I_ERI_G4z_S_G4x_S;
  Double I_ERI_G4x_Px_G3xy_S = I_ERI_H5x_S_G3xy_S+ABX*I_ERI_G4x_S_G3xy_S;
  Double I_ERI_G3xy_Px_G3xy_S = I_ERI_H4xy_S_G3xy_S+ABX*I_ERI_G3xy_S_G3xy_S;
  Double I_ERI_G3xz_Px_G3xy_S = I_ERI_H4xz_S_G3xy_S+ABX*I_ERI_G3xz_S_G3xy_S;
  Double I_ERI_G2x2y_Px_G3xy_S = I_ERI_H3x2y_S_G3xy_S+ABX*I_ERI_G2x2y_S_G3xy_S;
  Double I_ERI_G2xyz_Px_G3xy_S = I_ERI_H3xyz_S_G3xy_S+ABX*I_ERI_G2xyz_S_G3xy_S;
  Double I_ERI_G2x2z_Px_G3xy_S = I_ERI_H3x2z_S_G3xy_S+ABX*I_ERI_G2x2z_S_G3xy_S;
  Double I_ERI_Gx3y_Px_G3xy_S = I_ERI_H2x3y_S_G3xy_S+ABX*I_ERI_Gx3y_S_G3xy_S;
  Double I_ERI_Gx2yz_Px_G3xy_S = I_ERI_H2x2yz_S_G3xy_S+ABX*I_ERI_Gx2yz_S_G3xy_S;
  Double I_ERI_Gxy2z_Px_G3xy_S = I_ERI_H2xy2z_S_G3xy_S+ABX*I_ERI_Gxy2z_S_G3xy_S;
  Double I_ERI_Gx3z_Px_G3xy_S = I_ERI_H2x3z_S_G3xy_S+ABX*I_ERI_Gx3z_S_G3xy_S;
  Double I_ERI_G4y_Px_G3xy_S = I_ERI_Hx4y_S_G3xy_S+ABX*I_ERI_G4y_S_G3xy_S;
  Double I_ERI_G3yz_Px_G3xy_S = I_ERI_Hx3yz_S_G3xy_S+ABX*I_ERI_G3yz_S_G3xy_S;
  Double I_ERI_G2y2z_Px_G3xy_S = I_ERI_Hx2y2z_S_G3xy_S+ABX*I_ERI_G2y2z_S_G3xy_S;
  Double I_ERI_Gy3z_Px_G3xy_S = I_ERI_Hxy3z_S_G3xy_S+ABX*I_ERI_Gy3z_S_G3xy_S;
  Double I_ERI_G4z_Px_G3xy_S = I_ERI_Hx4z_S_G3xy_S+ABX*I_ERI_G4z_S_G3xy_S;
  Double I_ERI_G4x_Py_G3xy_S = I_ERI_H4xy_S_G3xy_S+ABY*I_ERI_G4x_S_G3xy_S;
  Double I_ERI_G3xy_Py_G3xy_S = I_ERI_H3x2y_S_G3xy_S+ABY*I_ERI_G3xy_S_G3xy_S;
  Double I_ERI_G3xz_Py_G3xy_S = I_ERI_H3xyz_S_G3xy_S+ABY*I_ERI_G3xz_S_G3xy_S;
  Double I_ERI_G2x2y_Py_G3xy_S = I_ERI_H2x3y_S_G3xy_S+ABY*I_ERI_G2x2y_S_G3xy_S;
  Double I_ERI_G2xyz_Py_G3xy_S = I_ERI_H2x2yz_S_G3xy_S+ABY*I_ERI_G2xyz_S_G3xy_S;
  Double I_ERI_G2x2z_Py_G3xy_S = I_ERI_H2xy2z_S_G3xy_S+ABY*I_ERI_G2x2z_S_G3xy_S;
  Double I_ERI_Gx3y_Py_G3xy_S = I_ERI_Hx4y_S_G3xy_S+ABY*I_ERI_Gx3y_S_G3xy_S;
  Double I_ERI_Gx2yz_Py_G3xy_S = I_ERI_Hx3yz_S_G3xy_S+ABY*I_ERI_Gx2yz_S_G3xy_S;
  Double I_ERI_Gxy2z_Py_G3xy_S = I_ERI_Hx2y2z_S_G3xy_S+ABY*I_ERI_Gxy2z_S_G3xy_S;
  Double I_ERI_Gx3z_Py_G3xy_S = I_ERI_Hxy3z_S_G3xy_S+ABY*I_ERI_Gx3z_S_G3xy_S;
  Double I_ERI_G4y_Py_G3xy_S = I_ERI_H5y_S_G3xy_S+ABY*I_ERI_G4y_S_G3xy_S;
  Double I_ERI_G3yz_Py_G3xy_S = I_ERI_H4yz_S_G3xy_S+ABY*I_ERI_G3yz_S_G3xy_S;
  Double I_ERI_G2y2z_Py_G3xy_S = I_ERI_H3y2z_S_G3xy_S+ABY*I_ERI_G2y2z_S_G3xy_S;
  Double I_ERI_Gy3z_Py_G3xy_S = I_ERI_H2y3z_S_G3xy_S+ABY*I_ERI_Gy3z_S_G3xy_S;
  Double I_ERI_G4z_Py_G3xy_S = I_ERI_Hy4z_S_G3xy_S+ABY*I_ERI_G4z_S_G3xy_S;
  Double I_ERI_G4x_Pz_G3xy_S = I_ERI_H4xz_S_G3xy_S+ABZ*I_ERI_G4x_S_G3xy_S;
  Double I_ERI_G3xy_Pz_G3xy_S = I_ERI_H3xyz_S_G3xy_S+ABZ*I_ERI_G3xy_S_G3xy_S;
  Double I_ERI_G3xz_Pz_G3xy_S = I_ERI_H3x2z_S_G3xy_S+ABZ*I_ERI_G3xz_S_G3xy_S;
  Double I_ERI_G2x2y_Pz_G3xy_S = I_ERI_H2x2yz_S_G3xy_S+ABZ*I_ERI_G2x2y_S_G3xy_S;
  Double I_ERI_G2xyz_Pz_G3xy_S = I_ERI_H2xy2z_S_G3xy_S+ABZ*I_ERI_G2xyz_S_G3xy_S;
  Double I_ERI_G2x2z_Pz_G3xy_S = I_ERI_H2x3z_S_G3xy_S+ABZ*I_ERI_G2x2z_S_G3xy_S;
  Double I_ERI_Gx3y_Pz_G3xy_S = I_ERI_Hx3yz_S_G3xy_S+ABZ*I_ERI_Gx3y_S_G3xy_S;
  Double I_ERI_Gx2yz_Pz_G3xy_S = I_ERI_Hx2y2z_S_G3xy_S+ABZ*I_ERI_Gx2yz_S_G3xy_S;
  Double I_ERI_Gxy2z_Pz_G3xy_S = I_ERI_Hxy3z_S_G3xy_S+ABZ*I_ERI_Gxy2z_S_G3xy_S;
  Double I_ERI_Gx3z_Pz_G3xy_S = I_ERI_Hx4z_S_G3xy_S+ABZ*I_ERI_Gx3z_S_G3xy_S;
  Double I_ERI_G4y_Pz_G3xy_S = I_ERI_H4yz_S_G3xy_S+ABZ*I_ERI_G4y_S_G3xy_S;
  Double I_ERI_G3yz_Pz_G3xy_S = I_ERI_H3y2z_S_G3xy_S+ABZ*I_ERI_G3yz_S_G3xy_S;
  Double I_ERI_G2y2z_Pz_G3xy_S = I_ERI_H2y3z_S_G3xy_S+ABZ*I_ERI_G2y2z_S_G3xy_S;
  Double I_ERI_Gy3z_Pz_G3xy_S = I_ERI_Hy4z_S_G3xy_S+ABZ*I_ERI_Gy3z_S_G3xy_S;
  Double I_ERI_G4z_Pz_G3xy_S = I_ERI_H5z_S_G3xy_S+ABZ*I_ERI_G4z_S_G3xy_S;
  Double I_ERI_G4x_Px_G3xz_S = I_ERI_H5x_S_G3xz_S+ABX*I_ERI_G4x_S_G3xz_S;
  Double I_ERI_G3xy_Px_G3xz_S = I_ERI_H4xy_S_G3xz_S+ABX*I_ERI_G3xy_S_G3xz_S;
  Double I_ERI_G3xz_Px_G3xz_S = I_ERI_H4xz_S_G3xz_S+ABX*I_ERI_G3xz_S_G3xz_S;
  Double I_ERI_G2x2y_Px_G3xz_S = I_ERI_H3x2y_S_G3xz_S+ABX*I_ERI_G2x2y_S_G3xz_S;
  Double I_ERI_G2xyz_Px_G3xz_S = I_ERI_H3xyz_S_G3xz_S+ABX*I_ERI_G2xyz_S_G3xz_S;
  Double I_ERI_G2x2z_Px_G3xz_S = I_ERI_H3x2z_S_G3xz_S+ABX*I_ERI_G2x2z_S_G3xz_S;
  Double I_ERI_Gx3y_Px_G3xz_S = I_ERI_H2x3y_S_G3xz_S+ABX*I_ERI_Gx3y_S_G3xz_S;
  Double I_ERI_Gx2yz_Px_G3xz_S = I_ERI_H2x2yz_S_G3xz_S+ABX*I_ERI_Gx2yz_S_G3xz_S;
  Double I_ERI_Gxy2z_Px_G3xz_S = I_ERI_H2xy2z_S_G3xz_S+ABX*I_ERI_Gxy2z_S_G3xz_S;
  Double I_ERI_Gx3z_Px_G3xz_S = I_ERI_H2x3z_S_G3xz_S+ABX*I_ERI_Gx3z_S_G3xz_S;
  Double I_ERI_G4y_Px_G3xz_S = I_ERI_Hx4y_S_G3xz_S+ABX*I_ERI_G4y_S_G3xz_S;
  Double I_ERI_G3yz_Px_G3xz_S = I_ERI_Hx3yz_S_G3xz_S+ABX*I_ERI_G3yz_S_G3xz_S;
  Double I_ERI_G2y2z_Px_G3xz_S = I_ERI_Hx2y2z_S_G3xz_S+ABX*I_ERI_G2y2z_S_G3xz_S;
  Double I_ERI_Gy3z_Px_G3xz_S = I_ERI_Hxy3z_S_G3xz_S+ABX*I_ERI_Gy3z_S_G3xz_S;
  Double I_ERI_G4z_Px_G3xz_S = I_ERI_Hx4z_S_G3xz_S+ABX*I_ERI_G4z_S_G3xz_S;
  Double I_ERI_G4x_Py_G3xz_S = I_ERI_H4xy_S_G3xz_S+ABY*I_ERI_G4x_S_G3xz_S;
  Double I_ERI_G3xy_Py_G3xz_S = I_ERI_H3x2y_S_G3xz_S+ABY*I_ERI_G3xy_S_G3xz_S;
  Double I_ERI_G3xz_Py_G3xz_S = I_ERI_H3xyz_S_G3xz_S+ABY*I_ERI_G3xz_S_G3xz_S;
  Double I_ERI_G2x2y_Py_G3xz_S = I_ERI_H2x3y_S_G3xz_S+ABY*I_ERI_G2x2y_S_G3xz_S;
  Double I_ERI_G2xyz_Py_G3xz_S = I_ERI_H2x2yz_S_G3xz_S+ABY*I_ERI_G2xyz_S_G3xz_S;
  Double I_ERI_G2x2z_Py_G3xz_S = I_ERI_H2xy2z_S_G3xz_S+ABY*I_ERI_G2x2z_S_G3xz_S;
  Double I_ERI_Gx3y_Py_G3xz_S = I_ERI_Hx4y_S_G3xz_S+ABY*I_ERI_Gx3y_S_G3xz_S;
  Double I_ERI_Gx2yz_Py_G3xz_S = I_ERI_Hx3yz_S_G3xz_S+ABY*I_ERI_Gx2yz_S_G3xz_S;
  Double I_ERI_Gxy2z_Py_G3xz_S = I_ERI_Hx2y2z_S_G3xz_S+ABY*I_ERI_Gxy2z_S_G3xz_S;
  Double I_ERI_Gx3z_Py_G3xz_S = I_ERI_Hxy3z_S_G3xz_S+ABY*I_ERI_Gx3z_S_G3xz_S;
  Double I_ERI_G4y_Py_G3xz_S = I_ERI_H5y_S_G3xz_S+ABY*I_ERI_G4y_S_G3xz_S;
  Double I_ERI_G3yz_Py_G3xz_S = I_ERI_H4yz_S_G3xz_S+ABY*I_ERI_G3yz_S_G3xz_S;
  Double I_ERI_G2y2z_Py_G3xz_S = I_ERI_H3y2z_S_G3xz_S+ABY*I_ERI_G2y2z_S_G3xz_S;
  Double I_ERI_Gy3z_Py_G3xz_S = I_ERI_H2y3z_S_G3xz_S+ABY*I_ERI_Gy3z_S_G3xz_S;
  Double I_ERI_G4z_Py_G3xz_S = I_ERI_Hy4z_S_G3xz_S+ABY*I_ERI_G4z_S_G3xz_S;
  Double I_ERI_G4x_Pz_G3xz_S = I_ERI_H4xz_S_G3xz_S+ABZ*I_ERI_G4x_S_G3xz_S;
  Double I_ERI_G3xy_Pz_G3xz_S = I_ERI_H3xyz_S_G3xz_S+ABZ*I_ERI_G3xy_S_G3xz_S;
  Double I_ERI_G3xz_Pz_G3xz_S = I_ERI_H3x2z_S_G3xz_S+ABZ*I_ERI_G3xz_S_G3xz_S;
  Double I_ERI_G2x2y_Pz_G3xz_S = I_ERI_H2x2yz_S_G3xz_S+ABZ*I_ERI_G2x2y_S_G3xz_S;
  Double I_ERI_G2xyz_Pz_G3xz_S = I_ERI_H2xy2z_S_G3xz_S+ABZ*I_ERI_G2xyz_S_G3xz_S;
  Double I_ERI_G2x2z_Pz_G3xz_S = I_ERI_H2x3z_S_G3xz_S+ABZ*I_ERI_G2x2z_S_G3xz_S;
  Double I_ERI_Gx3y_Pz_G3xz_S = I_ERI_Hx3yz_S_G3xz_S+ABZ*I_ERI_Gx3y_S_G3xz_S;
  Double I_ERI_Gx2yz_Pz_G3xz_S = I_ERI_Hx2y2z_S_G3xz_S+ABZ*I_ERI_Gx2yz_S_G3xz_S;
  Double I_ERI_Gxy2z_Pz_G3xz_S = I_ERI_Hxy3z_S_G3xz_S+ABZ*I_ERI_Gxy2z_S_G3xz_S;
  Double I_ERI_Gx3z_Pz_G3xz_S = I_ERI_Hx4z_S_G3xz_S+ABZ*I_ERI_Gx3z_S_G3xz_S;
  Double I_ERI_G4y_Pz_G3xz_S = I_ERI_H4yz_S_G3xz_S+ABZ*I_ERI_G4y_S_G3xz_S;
  Double I_ERI_G3yz_Pz_G3xz_S = I_ERI_H3y2z_S_G3xz_S+ABZ*I_ERI_G3yz_S_G3xz_S;
  Double I_ERI_G2y2z_Pz_G3xz_S = I_ERI_H2y3z_S_G3xz_S+ABZ*I_ERI_G2y2z_S_G3xz_S;
  Double I_ERI_Gy3z_Pz_G3xz_S = I_ERI_Hy4z_S_G3xz_S+ABZ*I_ERI_Gy3z_S_G3xz_S;
  Double I_ERI_G4z_Pz_G3xz_S = I_ERI_H5z_S_G3xz_S+ABZ*I_ERI_G4z_S_G3xz_S;
  Double I_ERI_G4x_Px_G2x2y_S = I_ERI_H5x_S_G2x2y_S+ABX*I_ERI_G4x_S_G2x2y_S;
  Double I_ERI_G3xy_Px_G2x2y_S = I_ERI_H4xy_S_G2x2y_S+ABX*I_ERI_G3xy_S_G2x2y_S;
  Double I_ERI_G3xz_Px_G2x2y_S = I_ERI_H4xz_S_G2x2y_S+ABX*I_ERI_G3xz_S_G2x2y_S;
  Double I_ERI_G2x2y_Px_G2x2y_S = I_ERI_H3x2y_S_G2x2y_S+ABX*I_ERI_G2x2y_S_G2x2y_S;
  Double I_ERI_G2xyz_Px_G2x2y_S = I_ERI_H3xyz_S_G2x2y_S+ABX*I_ERI_G2xyz_S_G2x2y_S;
  Double I_ERI_G2x2z_Px_G2x2y_S = I_ERI_H3x2z_S_G2x2y_S+ABX*I_ERI_G2x2z_S_G2x2y_S;
  Double I_ERI_Gx3y_Px_G2x2y_S = I_ERI_H2x3y_S_G2x2y_S+ABX*I_ERI_Gx3y_S_G2x2y_S;
  Double I_ERI_Gx2yz_Px_G2x2y_S = I_ERI_H2x2yz_S_G2x2y_S+ABX*I_ERI_Gx2yz_S_G2x2y_S;
  Double I_ERI_Gxy2z_Px_G2x2y_S = I_ERI_H2xy2z_S_G2x2y_S+ABX*I_ERI_Gxy2z_S_G2x2y_S;
  Double I_ERI_Gx3z_Px_G2x2y_S = I_ERI_H2x3z_S_G2x2y_S+ABX*I_ERI_Gx3z_S_G2x2y_S;
  Double I_ERI_G4y_Px_G2x2y_S = I_ERI_Hx4y_S_G2x2y_S+ABX*I_ERI_G4y_S_G2x2y_S;
  Double I_ERI_G3yz_Px_G2x2y_S = I_ERI_Hx3yz_S_G2x2y_S+ABX*I_ERI_G3yz_S_G2x2y_S;
  Double I_ERI_G2y2z_Px_G2x2y_S = I_ERI_Hx2y2z_S_G2x2y_S+ABX*I_ERI_G2y2z_S_G2x2y_S;
  Double I_ERI_Gy3z_Px_G2x2y_S = I_ERI_Hxy3z_S_G2x2y_S+ABX*I_ERI_Gy3z_S_G2x2y_S;
  Double I_ERI_G4z_Px_G2x2y_S = I_ERI_Hx4z_S_G2x2y_S+ABX*I_ERI_G4z_S_G2x2y_S;
  Double I_ERI_G4x_Py_G2x2y_S = I_ERI_H4xy_S_G2x2y_S+ABY*I_ERI_G4x_S_G2x2y_S;
  Double I_ERI_G3xy_Py_G2x2y_S = I_ERI_H3x2y_S_G2x2y_S+ABY*I_ERI_G3xy_S_G2x2y_S;
  Double I_ERI_G3xz_Py_G2x2y_S = I_ERI_H3xyz_S_G2x2y_S+ABY*I_ERI_G3xz_S_G2x2y_S;
  Double I_ERI_G2x2y_Py_G2x2y_S = I_ERI_H2x3y_S_G2x2y_S+ABY*I_ERI_G2x2y_S_G2x2y_S;
  Double I_ERI_G2xyz_Py_G2x2y_S = I_ERI_H2x2yz_S_G2x2y_S+ABY*I_ERI_G2xyz_S_G2x2y_S;
  Double I_ERI_G2x2z_Py_G2x2y_S = I_ERI_H2xy2z_S_G2x2y_S+ABY*I_ERI_G2x2z_S_G2x2y_S;
  Double I_ERI_Gx3y_Py_G2x2y_S = I_ERI_Hx4y_S_G2x2y_S+ABY*I_ERI_Gx3y_S_G2x2y_S;
  Double I_ERI_Gx2yz_Py_G2x2y_S = I_ERI_Hx3yz_S_G2x2y_S+ABY*I_ERI_Gx2yz_S_G2x2y_S;
  Double I_ERI_Gxy2z_Py_G2x2y_S = I_ERI_Hx2y2z_S_G2x2y_S+ABY*I_ERI_Gxy2z_S_G2x2y_S;
  Double I_ERI_Gx3z_Py_G2x2y_S = I_ERI_Hxy3z_S_G2x2y_S+ABY*I_ERI_Gx3z_S_G2x2y_S;
  Double I_ERI_G4y_Py_G2x2y_S = I_ERI_H5y_S_G2x2y_S+ABY*I_ERI_G4y_S_G2x2y_S;
  Double I_ERI_G3yz_Py_G2x2y_S = I_ERI_H4yz_S_G2x2y_S+ABY*I_ERI_G3yz_S_G2x2y_S;
  Double I_ERI_G2y2z_Py_G2x2y_S = I_ERI_H3y2z_S_G2x2y_S+ABY*I_ERI_G2y2z_S_G2x2y_S;
  Double I_ERI_Gy3z_Py_G2x2y_S = I_ERI_H2y3z_S_G2x2y_S+ABY*I_ERI_Gy3z_S_G2x2y_S;
  Double I_ERI_G4z_Py_G2x2y_S = I_ERI_Hy4z_S_G2x2y_S+ABY*I_ERI_G4z_S_G2x2y_S;
  Double I_ERI_G4x_Pz_G2x2y_S = I_ERI_H4xz_S_G2x2y_S+ABZ*I_ERI_G4x_S_G2x2y_S;
  Double I_ERI_G3xy_Pz_G2x2y_S = I_ERI_H3xyz_S_G2x2y_S+ABZ*I_ERI_G3xy_S_G2x2y_S;
  Double I_ERI_G3xz_Pz_G2x2y_S = I_ERI_H3x2z_S_G2x2y_S+ABZ*I_ERI_G3xz_S_G2x2y_S;
  Double I_ERI_G2x2y_Pz_G2x2y_S = I_ERI_H2x2yz_S_G2x2y_S+ABZ*I_ERI_G2x2y_S_G2x2y_S;
  Double I_ERI_G2xyz_Pz_G2x2y_S = I_ERI_H2xy2z_S_G2x2y_S+ABZ*I_ERI_G2xyz_S_G2x2y_S;
  Double I_ERI_G2x2z_Pz_G2x2y_S = I_ERI_H2x3z_S_G2x2y_S+ABZ*I_ERI_G2x2z_S_G2x2y_S;
  Double I_ERI_Gx3y_Pz_G2x2y_S = I_ERI_Hx3yz_S_G2x2y_S+ABZ*I_ERI_Gx3y_S_G2x2y_S;
  Double I_ERI_Gx2yz_Pz_G2x2y_S = I_ERI_Hx2y2z_S_G2x2y_S+ABZ*I_ERI_Gx2yz_S_G2x2y_S;
  Double I_ERI_Gxy2z_Pz_G2x2y_S = I_ERI_Hxy3z_S_G2x2y_S+ABZ*I_ERI_Gxy2z_S_G2x2y_S;
  Double I_ERI_Gx3z_Pz_G2x2y_S = I_ERI_Hx4z_S_G2x2y_S+ABZ*I_ERI_Gx3z_S_G2x2y_S;
  Double I_ERI_G4y_Pz_G2x2y_S = I_ERI_H4yz_S_G2x2y_S+ABZ*I_ERI_G4y_S_G2x2y_S;
  Double I_ERI_G3yz_Pz_G2x2y_S = I_ERI_H3y2z_S_G2x2y_S+ABZ*I_ERI_G3yz_S_G2x2y_S;
  Double I_ERI_G2y2z_Pz_G2x2y_S = I_ERI_H2y3z_S_G2x2y_S+ABZ*I_ERI_G2y2z_S_G2x2y_S;
  Double I_ERI_Gy3z_Pz_G2x2y_S = I_ERI_Hy4z_S_G2x2y_S+ABZ*I_ERI_Gy3z_S_G2x2y_S;
  Double I_ERI_G4z_Pz_G2x2y_S = I_ERI_H5z_S_G2x2y_S+ABZ*I_ERI_G4z_S_G2x2y_S;
  Double I_ERI_G4x_Px_G2xyz_S = I_ERI_H5x_S_G2xyz_S+ABX*I_ERI_G4x_S_G2xyz_S;
  Double I_ERI_G3xy_Px_G2xyz_S = I_ERI_H4xy_S_G2xyz_S+ABX*I_ERI_G3xy_S_G2xyz_S;
  Double I_ERI_G3xz_Px_G2xyz_S = I_ERI_H4xz_S_G2xyz_S+ABX*I_ERI_G3xz_S_G2xyz_S;
  Double I_ERI_G2x2y_Px_G2xyz_S = I_ERI_H3x2y_S_G2xyz_S+ABX*I_ERI_G2x2y_S_G2xyz_S;
  Double I_ERI_G2xyz_Px_G2xyz_S = I_ERI_H3xyz_S_G2xyz_S+ABX*I_ERI_G2xyz_S_G2xyz_S;
  Double I_ERI_G2x2z_Px_G2xyz_S = I_ERI_H3x2z_S_G2xyz_S+ABX*I_ERI_G2x2z_S_G2xyz_S;
  Double I_ERI_Gx3y_Px_G2xyz_S = I_ERI_H2x3y_S_G2xyz_S+ABX*I_ERI_Gx3y_S_G2xyz_S;
  Double I_ERI_Gx2yz_Px_G2xyz_S = I_ERI_H2x2yz_S_G2xyz_S+ABX*I_ERI_Gx2yz_S_G2xyz_S;
  Double I_ERI_Gxy2z_Px_G2xyz_S = I_ERI_H2xy2z_S_G2xyz_S+ABX*I_ERI_Gxy2z_S_G2xyz_S;
  Double I_ERI_Gx3z_Px_G2xyz_S = I_ERI_H2x3z_S_G2xyz_S+ABX*I_ERI_Gx3z_S_G2xyz_S;
  Double I_ERI_G4y_Px_G2xyz_S = I_ERI_Hx4y_S_G2xyz_S+ABX*I_ERI_G4y_S_G2xyz_S;
  Double I_ERI_G3yz_Px_G2xyz_S = I_ERI_Hx3yz_S_G2xyz_S+ABX*I_ERI_G3yz_S_G2xyz_S;
  Double I_ERI_G2y2z_Px_G2xyz_S = I_ERI_Hx2y2z_S_G2xyz_S+ABX*I_ERI_G2y2z_S_G2xyz_S;
  Double I_ERI_Gy3z_Px_G2xyz_S = I_ERI_Hxy3z_S_G2xyz_S+ABX*I_ERI_Gy3z_S_G2xyz_S;
  Double I_ERI_G4z_Px_G2xyz_S = I_ERI_Hx4z_S_G2xyz_S+ABX*I_ERI_G4z_S_G2xyz_S;
  Double I_ERI_G4x_Py_G2xyz_S = I_ERI_H4xy_S_G2xyz_S+ABY*I_ERI_G4x_S_G2xyz_S;
  Double I_ERI_G3xy_Py_G2xyz_S = I_ERI_H3x2y_S_G2xyz_S+ABY*I_ERI_G3xy_S_G2xyz_S;
  Double I_ERI_G3xz_Py_G2xyz_S = I_ERI_H3xyz_S_G2xyz_S+ABY*I_ERI_G3xz_S_G2xyz_S;
  Double I_ERI_G2x2y_Py_G2xyz_S = I_ERI_H2x3y_S_G2xyz_S+ABY*I_ERI_G2x2y_S_G2xyz_S;
  Double I_ERI_G2xyz_Py_G2xyz_S = I_ERI_H2x2yz_S_G2xyz_S+ABY*I_ERI_G2xyz_S_G2xyz_S;
  Double I_ERI_G2x2z_Py_G2xyz_S = I_ERI_H2xy2z_S_G2xyz_S+ABY*I_ERI_G2x2z_S_G2xyz_S;
  Double I_ERI_Gx3y_Py_G2xyz_S = I_ERI_Hx4y_S_G2xyz_S+ABY*I_ERI_Gx3y_S_G2xyz_S;
  Double I_ERI_Gx2yz_Py_G2xyz_S = I_ERI_Hx3yz_S_G2xyz_S+ABY*I_ERI_Gx2yz_S_G2xyz_S;
  Double I_ERI_Gxy2z_Py_G2xyz_S = I_ERI_Hx2y2z_S_G2xyz_S+ABY*I_ERI_Gxy2z_S_G2xyz_S;
  Double I_ERI_Gx3z_Py_G2xyz_S = I_ERI_Hxy3z_S_G2xyz_S+ABY*I_ERI_Gx3z_S_G2xyz_S;
  Double I_ERI_G4y_Py_G2xyz_S = I_ERI_H5y_S_G2xyz_S+ABY*I_ERI_G4y_S_G2xyz_S;
  Double I_ERI_G3yz_Py_G2xyz_S = I_ERI_H4yz_S_G2xyz_S+ABY*I_ERI_G3yz_S_G2xyz_S;
  Double I_ERI_G2y2z_Py_G2xyz_S = I_ERI_H3y2z_S_G2xyz_S+ABY*I_ERI_G2y2z_S_G2xyz_S;
  Double I_ERI_Gy3z_Py_G2xyz_S = I_ERI_H2y3z_S_G2xyz_S+ABY*I_ERI_Gy3z_S_G2xyz_S;
  Double I_ERI_G4z_Py_G2xyz_S = I_ERI_Hy4z_S_G2xyz_S+ABY*I_ERI_G4z_S_G2xyz_S;
  Double I_ERI_G4x_Pz_G2xyz_S = I_ERI_H4xz_S_G2xyz_S+ABZ*I_ERI_G4x_S_G2xyz_S;
  Double I_ERI_G3xy_Pz_G2xyz_S = I_ERI_H3xyz_S_G2xyz_S+ABZ*I_ERI_G3xy_S_G2xyz_S;
  Double I_ERI_G3xz_Pz_G2xyz_S = I_ERI_H3x2z_S_G2xyz_S+ABZ*I_ERI_G3xz_S_G2xyz_S;
  Double I_ERI_G2x2y_Pz_G2xyz_S = I_ERI_H2x2yz_S_G2xyz_S+ABZ*I_ERI_G2x2y_S_G2xyz_S;
  Double I_ERI_G2xyz_Pz_G2xyz_S = I_ERI_H2xy2z_S_G2xyz_S+ABZ*I_ERI_G2xyz_S_G2xyz_S;
  Double I_ERI_G2x2z_Pz_G2xyz_S = I_ERI_H2x3z_S_G2xyz_S+ABZ*I_ERI_G2x2z_S_G2xyz_S;
  Double I_ERI_Gx3y_Pz_G2xyz_S = I_ERI_Hx3yz_S_G2xyz_S+ABZ*I_ERI_Gx3y_S_G2xyz_S;
  Double I_ERI_Gx2yz_Pz_G2xyz_S = I_ERI_Hx2y2z_S_G2xyz_S+ABZ*I_ERI_Gx2yz_S_G2xyz_S;
  Double I_ERI_Gxy2z_Pz_G2xyz_S = I_ERI_Hxy3z_S_G2xyz_S+ABZ*I_ERI_Gxy2z_S_G2xyz_S;
  Double I_ERI_Gx3z_Pz_G2xyz_S = I_ERI_Hx4z_S_G2xyz_S+ABZ*I_ERI_Gx3z_S_G2xyz_S;
  Double I_ERI_G4y_Pz_G2xyz_S = I_ERI_H4yz_S_G2xyz_S+ABZ*I_ERI_G4y_S_G2xyz_S;
  Double I_ERI_G3yz_Pz_G2xyz_S = I_ERI_H3y2z_S_G2xyz_S+ABZ*I_ERI_G3yz_S_G2xyz_S;
  Double I_ERI_G2y2z_Pz_G2xyz_S = I_ERI_H2y3z_S_G2xyz_S+ABZ*I_ERI_G2y2z_S_G2xyz_S;
  Double I_ERI_Gy3z_Pz_G2xyz_S = I_ERI_Hy4z_S_G2xyz_S+ABZ*I_ERI_Gy3z_S_G2xyz_S;
  Double I_ERI_G4z_Pz_G2xyz_S = I_ERI_H5z_S_G2xyz_S+ABZ*I_ERI_G4z_S_G2xyz_S;
  Double I_ERI_G4x_Px_G2x2z_S = I_ERI_H5x_S_G2x2z_S+ABX*I_ERI_G4x_S_G2x2z_S;
  Double I_ERI_G3xy_Px_G2x2z_S = I_ERI_H4xy_S_G2x2z_S+ABX*I_ERI_G3xy_S_G2x2z_S;
  Double I_ERI_G3xz_Px_G2x2z_S = I_ERI_H4xz_S_G2x2z_S+ABX*I_ERI_G3xz_S_G2x2z_S;
  Double I_ERI_G2x2y_Px_G2x2z_S = I_ERI_H3x2y_S_G2x2z_S+ABX*I_ERI_G2x2y_S_G2x2z_S;
  Double I_ERI_G2xyz_Px_G2x2z_S = I_ERI_H3xyz_S_G2x2z_S+ABX*I_ERI_G2xyz_S_G2x2z_S;
  Double I_ERI_G2x2z_Px_G2x2z_S = I_ERI_H3x2z_S_G2x2z_S+ABX*I_ERI_G2x2z_S_G2x2z_S;
  Double I_ERI_Gx3y_Px_G2x2z_S = I_ERI_H2x3y_S_G2x2z_S+ABX*I_ERI_Gx3y_S_G2x2z_S;
  Double I_ERI_Gx2yz_Px_G2x2z_S = I_ERI_H2x2yz_S_G2x2z_S+ABX*I_ERI_Gx2yz_S_G2x2z_S;
  Double I_ERI_Gxy2z_Px_G2x2z_S = I_ERI_H2xy2z_S_G2x2z_S+ABX*I_ERI_Gxy2z_S_G2x2z_S;
  Double I_ERI_Gx3z_Px_G2x2z_S = I_ERI_H2x3z_S_G2x2z_S+ABX*I_ERI_Gx3z_S_G2x2z_S;
  Double I_ERI_G4y_Px_G2x2z_S = I_ERI_Hx4y_S_G2x2z_S+ABX*I_ERI_G4y_S_G2x2z_S;
  Double I_ERI_G3yz_Px_G2x2z_S = I_ERI_Hx3yz_S_G2x2z_S+ABX*I_ERI_G3yz_S_G2x2z_S;
  Double I_ERI_G2y2z_Px_G2x2z_S = I_ERI_Hx2y2z_S_G2x2z_S+ABX*I_ERI_G2y2z_S_G2x2z_S;
  Double I_ERI_Gy3z_Px_G2x2z_S = I_ERI_Hxy3z_S_G2x2z_S+ABX*I_ERI_Gy3z_S_G2x2z_S;
  Double I_ERI_G4z_Px_G2x2z_S = I_ERI_Hx4z_S_G2x2z_S+ABX*I_ERI_G4z_S_G2x2z_S;
  Double I_ERI_G4x_Py_G2x2z_S = I_ERI_H4xy_S_G2x2z_S+ABY*I_ERI_G4x_S_G2x2z_S;
  Double I_ERI_G3xy_Py_G2x2z_S = I_ERI_H3x2y_S_G2x2z_S+ABY*I_ERI_G3xy_S_G2x2z_S;
  Double I_ERI_G3xz_Py_G2x2z_S = I_ERI_H3xyz_S_G2x2z_S+ABY*I_ERI_G3xz_S_G2x2z_S;
  Double I_ERI_G2x2y_Py_G2x2z_S = I_ERI_H2x3y_S_G2x2z_S+ABY*I_ERI_G2x2y_S_G2x2z_S;
  Double I_ERI_G2xyz_Py_G2x2z_S = I_ERI_H2x2yz_S_G2x2z_S+ABY*I_ERI_G2xyz_S_G2x2z_S;
  Double I_ERI_G2x2z_Py_G2x2z_S = I_ERI_H2xy2z_S_G2x2z_S+ABY*I_ERI_G2x2z_S_G2x2z_S;
  Double I_ERI_Gx3y_Py_G2x2z_S = I_ERI_Hx4y_S_G2x2z_S+ABY*I_ERI_Gx3y_S_G2x2z_S;
  Double I_ERI_Gx2yz_Py_G2x2z_S = I_ERI_Hx3yz_S_G2x2z_S+ABY*I_ERI_Gx2yz_S_G2x2z_S;
  Double I_ERI_Gxy2z_Py_G2x2z_S = I_ERI_Hx2y2z_S_G2x2z_S+ABY*I_ERI_Gxy2z_S_G2x2z_S;
  Double I_ERI_Gx3z_Py_G2x2z_S = I_ERI_Hxy3z_S_G2x2z_S+ABY*I_ERI_Gx3z_S_G2x2z_S;
  Double I_ERI_G4y_Py_G2x2z_S = I_ERI_H5y_S_G2x2z_S+ABY*I_ERI_G4y_S_G2x2z_S;
  Double I_ERI_G3yz_Py_G2x2z_S = I_ERI_H4yz_S_G2x2z_S+ABY*I_ERI_G3yz_S_G2x2z_S;
  Double I_ERI_G2y2z_Py_G2x2z_S = I_ERI_H3y2z_S_G2x2z_S+ABY*I_ERI_G2y2z_S_G2x2z_S;
  Double I_ERI_Gy3z_Py_G2x2z_S = I_ERI_H2y3z_S_G2x2z_S+ABY*I_ERI_Gy3z_S_G2x2z_S;
  Double I_ERI_G4z_Py_G2x2z_S = I_ERI_Hy4z_S_G2x2z_S+ABY*I_ERI_G4z_S_G2x2z_S;
  Double I_ERI_G4x_Pz_G2x2z_S = I_ERI_H4xz_S_G2x2z_S+ABZ*I_ERI_G4x_S_G2x2z_S;
  Double I_ERI_G3xy_Pz_G2x2z_S = I_ERI_H3xyz_S_G2x2z_S+ABZ*I_ERI_G3xy_S_G2x2z_S;
  Double I_ERI_G3xz_Pz_G2x2z_S = I_ERI_H3x2z_S_G2x2z_S+ABZ*I_ERI_G3xz_S_G2x2z_S;
  Double I_ERI_G2x2y_Pz_G2x2z_S = I_ERI_H2x2yz_S_G2x2z_S+ABZ*I_ERI_G2x2y_S_G2x2z_S;
  Double I_ERI_G2xyz_Pz_G2x2z_S = I_ERI_H2xy2z_S_G2x2z_S+ABZ*I_ERI_G2xyz_S_G2x2z_S;
  Double I_ERI_G2x2z_Pz_G2x2z_S = I_ERI_H2x3z_S_G2x2z_S+ABZ*I_ERI_G2x2z_S_G2x2z_S;
  Double I_ERI_Gx3y_Pz_G2x2z_S = I_ERI_Hx3yz_S_G2x2z_S+ABZ*I_ERI_Gx3y_S_G2x2z_S;
  Double I_ERI_Gx2yz_Pz_G2x2z_S = I_ERI_Hx2y2z_S_G2x2z_S+ABZ*I_ERI_Gx2yz_S_G2x2z_S;
  Double I_ERI_Gxy2z_Pz_G2x2z_S = I_ERI_Hxy3z_S_G2x2z_S+ABZ*I_ERI_Gxy2z_S_G2x2z_S;
  Double I_ERI_Gx3z_Pz_G2x2z_S = I_ERI_Hx4z_S_G2x2z_S+ABZ*I_ERI_Gx3z_S_G2x2z_S;
  Double I_ERI_G4y_Pz_G2x2z_S = I_ERI_H4yz_S_G2x2z_S+ABZ*I_ERI_G4y_S_G2x2z_S;
  Double I_ERI_G3yz_Pz_G2x2z_S = I_ERI_H3y2z_S_G2x2z_S+ABZ*I_ERI_G3yz_S_G2x2z_S;
  Double I_ERI_G2y2z_Pz_G2x2z_S = I_ERI_H2y3z_S_G2x2z_S+ABZ*I_ERI_G2y2z_S_G2x2z_S;
  Double I_ERI_Gy3z_Pz_G2x2z_S = I_ERI_Hy4z_S_G2x2z_S+ABZ*I_ERI_Gy3z_S_G2x2z_S;
  Double I_ERI_G4z_Pz_G2x2z_S = I_ERI_H5z_S_G2x2z_S+ABZ*I_ERI_G4z_S_G2x2z_S;
  Double I_ERI_G4x_Px_Gx3y_S = I_ERI_H5x_S_Gx3y_S+ABX*I_ERI_G4x_S_Gx3y_S;
  Double I_ERI_G3xy_Px_Gx3y_S = I_ERI_H4xy_S_Gx3y_S+ABX*I_ERI_G3xy_S_Gx3y_S;
  Double I_ERI_G3xz_Px_Gx3y_S = I_ERI_H4xz_S_Gx3y_S+ABX*I_ERI_G3xz_S_Gx3y_S;
  Double I_ERI_G2x2y_Px_Gx3y_S = I_ERI_H3x2y_S_Gx3y_S+ABX*I_ERI_G2x2y_S_Gx3y_S;
  Double I_ERI_G2xyz_Px_Gx3y_S = I_ERI_H3xyz_S_Gx3y_S+ABX*I_ERI_G2xyz_S_Gx3y_S;
  Double I_ERI_G2x2z_Px_Gx3y_S = I_ERI_H3x2z_S_Gx3y_S+ABX*I_ERI_G2x2z_S_Gx3y_S;
  Double I_ERI_Gx3y_Px_Gx3y_S = I_ERI_H2x3y_S_Gx3y_S+ABX*I_ERI_Gx3y_S_Gx3y_S;
  Double I_ERI_Gx2yz_Px_Gx3y_S = I_ERI_H2x2yz_S_Gx3y_S+ABX*I_ERI_Gx2yz_S_Gx3y_S;
  Double I_ERI_Gxy2z_Px_Gx3y_S = I_ERI_H2xy2z_S_Gx3y_S+ABX*I_ERI_Gxy2z_S_Gx3y_S;
  Double I_ERI_Gx3z_Px_Gx3y_S = I_ERI_H2x3z_S_Gx3y_S+ABX*I_ERI_Gx3z_S_Gx3y_S;
  Double I_ERI_G4y_Px_Gx3y_S = I_ERI_Hx4y_S_Gx3y_S+ABX*I_ERI_G4y_S_Gx3y_S;
  Double I_ERI_G3yz_Px_Gx3y_S = I_ERI_Hx3yz_S_Gx3y_S+ABX*I_ERI_G3yz_S_Gx3y_S;
  Double I_ERI_G2y2z_Px_Gx3y_S = I_ERI_Hx2y2z_S_Gx3y_S+ABX*I_ERI_G2y2z_S_Gx3y_S;
  Double I_ERI_Gy3z_Px_Gx3y_S = I_ERI_Hxy3z_S_Gx3y_S+ABX*I_ERI_Gy3z_S_Gx3y_S;
  Double I_ERI_G4z_Px_Gx3y_S = I_ERI_Hx4z_S_Gx3y_S+ABX*I_ERI_G4z_S_Gx3y_S;
  Double I_ERI_G4x_Py_Gx3y_S = I_ERI_H4xy_S_Gx3y_S+ABY*I_ERI_G4x_S_Gx3y_S;
  Double I_ERI_G3xy_Py_Gx3y_S = I_ERI_H3x2y_S_Gx3y_S+ABY*I_ERI_G3xy_S_Gx3y_S;
  Double I_ERI_G3xz_Py_Gx3y_S = I_ERI_H3xyz_S_Gx3y_S+ABY*I_ERI_G3xz_S_Gx3y_S;
  Double I_ERI_G2x2y_Py_Gx3y_S = I_ERI_H2x3y_S_Gx3y_S+ABY*I_ERI_G2x2y_S_Gx3y_S;
  Double I_ERI_G2xyz_Py_Gx3y_S = I_ERI_H2x2yz_S_Gx3y_S+ABY*I_ERI_G2xyz_S_Gx3y_S;
  Double I_ERI_G2x2z_Py_Gx3y_S = I_ERI_H2xy2z_S_Gx3y_S+ABY*I_ERI_G2x2z_S_Gx3y_S;
  Double I_ERI_Gx3y_Py_Gx3y_S = I_ERI_Hx4y_S_Gx3y_S+ABY*I_ERI_Gx3y_S_Gx3y_S;
  Double I_ERI_Gx2yz_Py_Gx3y_S = I_ERI_Hx3yz_S_Gx3y_S+ABY*I_ERI_Gx2yz_S_Gx3y_S;
  Double I_ERI_Gxy2z_Py_Gx3y_S = I_ERI_Hx2y2z_S_Gx3y_S+ABY*I_ERI_Gxy2z_S_Gx3y_S;
  Double I_ERI_Gx3z_Py_Gx3y_S = I_ERI_Hxy3z_S_Gx3y_S+ABY*I_ERI_Gx3z_S_Gx3y_S;
  Double I_ERI_G4y_Py_Gx3y_S = I_ERI_H5y_S_Gx3y_S+ABY*I_ERI_G4y_S_Gx3y_S;
  Double I_ERI_G3yz_Py_Gx3y_S = I_ERI_H4yz_S_Gx3y_S+ABY*I_ERI_G3yz_S_Gx3y_S;
  Double I_ERI_G2y2z_Py_Gx3y_S = I_ERI_H3y2z_S_Gx3y_S+ABY*I_ERI_G2y2z_S_Gx3y_S;
  Double I_ERI_Gy3z_Py_Gx3y_S = I_ERI_H2y3z_S_Gx3y_S+ABY*I_ERI_Gy3z_S_Gx3y_S;
  Double I_ERI_G4z_Py_Gx3y_S = I_ERI_Hy4z_S_Gx3y_S+ABY*I_ERI_G4z_S_Gx3y_S;
  Double I_ERI_G4x_Pz_Gx3y_S = I_ERI_H4xz_S_Gx3y_S+ABZ*I_ERI_G4x_S_Gx3y_S;
  Double I_ERI_G3xy_Pz_Gx3y_S = I_ERI_H3xyz_S_Gx3y_S+ABZ*I_ERI_G3xy_S_Gx3y_S;
  Double I_ERI_G3xz_Pz_Gx3y_S = I_ERI_H3x2z_S_Gx3y_S+ABZ*I_ERI_G3xz_S_Gx3y_S;
  Double I_ERI_G2x2y_Pz_Gx3y_S = I_ERI_H2x2yz_S_Gx3y_S+ABZ*I_ERI_G2x2y_S_Gx3y_S;
  Double I_ERI_G2xyz_Pz_Gx3y_S = I_ERI_H2xy2z_S_Gx3y_S+ABZ*I_ERI_G2xyz_S_Gx3y_S;
  Double I_ERI_G2x2z_Pz_Gx3y_S = I_ERI_H2x3z_S_Gx3y_S+ABZ*I_ERI_G2x2z_S_Gx3y_S;
  Double I_ERI_Gx3y_Pz_Gx3y_S = I_ERI_Hx3yz_S_Gx3y_S+ABZ*I_ERI_Gx3y_S_Gx3y_S;
  Double I_ERI_Gx2yz_Pz_Gx3y_S = I_ERI_Hx2y2z_S_Gx3y_S+ABZ*I_ERI_Gx2yz_S_Gx3y_S;
  Double I_ERI_Gxy2z_Pz_Gx3y_S = I_ERI_Hxy3z_S_Gx3y_S+ABZ*I_ERI_Gxy2z_S_Gx3y_S;
  Double I_ERI_Gx3z_Pz_Gx3y_S = I_ERI_Hx4z_S_Gx3y_S+ABZ*I_ERI_Gx3z_S_Gx3y_S;
  Double I_ERI_G4y_Pz_Gx3y_S = I_ERI_H4yz_S_Gx3y_S+ABZ*I_ERI_G4y_S_Gx3y_S;
  Double I_ERI_G3yz_Pz_Gx3y_S = I_ERI_H3y2z_S_Gx3y_S+ABZ*I_ERI_G3yz_S_Gx3y_S;
  Double I_ERI_G2y2z_Pz_Gx3y_S = I_ERI_H2y3z_S_Gx3y_S+ABZ*I_ERI_G2y2z_S_Gx3y_S;
  Double I_ERI_Gy3z_Pz_Gx3y_S = I_ERI_Hy4z_S_Gx3y_S+ABZ*I_ERI_Gy3z_S_Gx3y_S;
  Double I_ERI_G4z_Pz_Gx3y_S = I_ERI_H5z_S_Gx3y_S+ABZ*I_ERI_G4z_S_Gx3y_S;
  Double I_ERI_G4x_Px_Gx2yz_S = I_ERI_H5x_S_Gx2yz_S+ABX*I_ERI_G4x_S_Gx2yz_S;
  Double I_ERI_G3xy_Px_Gx2yz_S = I_ERI_H4xy_S_Gx2yz_S+ABX*I_ERI_G3xy_S_Gx2yz_S;
  Double I_ERI_G3xz_Px_Gx2yz_S = I_ERI_H4xz_S_Gx2yz_S+ABX*I_ERI_G3xz_S_Gx2yz_S;
  Double I_ERI_G2x2y_Px_Gx2yz_S = I_ERI_H3x2y_S_Gx2yz_S+ABX*I_ERI_G2x2y_S_Gx2yz_S;
  Double I_ERI_G2xyz_Px_Gx2yz_S = I_ERI_H3xyz_S_Gx2yz_S+ABX*I_ERI_G2xyz_S_Gx2yz_S;
  Double I_ERI_G2x2z_Px_Gx2yz_S = I_ERI_H3x2z_S_Gx2yz_S+ABX*I_ERI_G2x2z_S_Gx2yz_S;
  Double I_ERI_Gx3y_Px_Gx2yz_S = I_ERI_H2x3y_S_Gx2yz_S+ABX*I_ERI_Gx3y_S_Gx2yz_S;
  Double I_ERI_Gx2yz_Px_Gx2yz_S = I_ERI_H2x2yz_S_Gx2yz_S+ABX*I_ERI_Gx2yz_S_Gx2yz_S;
  Double I_ERI_Gxy2z_Px_Gx2yz_S = I_ERI_H2xy2z_S_Gx2yz_S+ABX*I_ERI_Gxy2z_S_Gx2yz_S;
  Double I_ERI_Gx3z_Px_Gx2yz_S = I_ERI_H2x3z_S_Gx2yz_S+ABX*I_ERI_Gx3z_S_Gx2yz_S;
  Double I_ERI_G4y_Px_Gx2yz_S = I_ERI_Hx4y_S_Gx2yz_S+ABX*I_ERI_G4y_S_Gx2yz_S;
  Double I_ERI_G3yz_Px_Gx2yz_S = I_ERI_Hx3yz_S_Gx2yz_S+ABX*I_ERI_G3yz_S_Gx2yz_S;
  Double I_ERI_G2y2z_Px_Gx2yz_S = I_ERI_Hx2y2z_S_Gx2yz_S+ABX*I_ERI_G2y2z_S_Gx2yz_S;
  Double I_ERI_Gy3z_Px_Gx2yz_S = I_ERI_Hxy3z_S_Gx2yz_S+ABX*I_ERI_Gy3z_S_Gx2yz_S;
  Double I_ERI_G4z_Px_Gx2yz_S = I_ERI_Hx4z_S_Gx2yz_S+ABX*I_ERI_G4z_S_Gx2yz_S;
  Double I_ERI_G4x_Py_Gx2yz_S = I_ERI_H4xy_S_Gx2yz_S+ABY*I_ERI_G4x_S_Gx2yz_S;
  Double I_ERI_G3xy_Py_Gx2yz_S = I_ERI_H3x2y_S_Gx2yz_S+ABY*I_ERI_G3xy_S_Gx2yz_S;
  Double I_ERI_G3xz_Py_Gx2yz_S = I_ERI_H3xyz_S_Gx2yz_S+ABY*I_ERI_G3xz_S_Gx2yz_S;
  Double I_ERI_G2x2y_Py_Gx2yz_S = I_ERI_H2x3y_S_Gx2yz_S+ABY*I_ERI_G2x2y_S_Gx2yz_S;
  Double I_ERI_G2xyz_Py_Gx2yz_S = I_ERI_H2x2yz_S_Gx2yz_S+ABY*I_ERI_G2xyz_S_Gx2yz_S;
  Double I_ERI_G2x2z_Py_Gx2yz_S = I_ERI_H2xy2z_S_Gx2yz_S+ABY*I_ERI_G2x2z_S_Gx2yz_S;
  Double I_ERI_Gx3y_Py_Gx2yz_S = I_ERI_Hx4y_S_Gx2yz_S+ABY*I_ERI_Gx3y_S_Gx2yz_S;
  Double I_ERI_Gx2yz_Py_Gx2yz_S = I_ERI_Hx3yz_S_Gx2yz_S+ABY*I_ERI_Gx2yz_S_Gx2yz_S;
  Double I_ERI_Gxy2z_Py_Gx2yz_S = I_ERI_Hx2y2z_S_Gx2yz_S+ABY*I_ERI_Gxy2z_S_Gx2yz_S;
  Double I_ERI_Gx3z_Py_Gx2yz_S = I_ERI_Hxy3z_S_Gx2yz_S+ABY*I_ERI_Gx3z_S_Gx2yz_S;
  Double I_ERI_G4y_Py_Gx2yz_S = I_ERI_H5y_S_Gx2yz_S+ABY*I_ERI_G4y_S_Gx2yz_S;
  Double I_ERI_G3yz_Py_Gx2yz_S = I_ERI_H4yz_S_Gx2yz_S+ABY*I_ERI_G3yz_S_Gx2yz_S;
  Double I_ERI_G2y2z_Py_Gx2yz_S = I_ERI_H3y2z_S_Gx2yz_S+ABY*I_ERI_G2y2z_S_Gx2yz_S;
  Double I_ERI_Gy3z_Py_Gx2yz_S = I_ERI_H2y3z_S_Gx2yz_S+ABY*I_ERI_Gy3z_S_Gx2yz_S;
  Double I_ERI_G4z_Py_Gx2yz_S = I_ERI_Hy4z_S_Gx2yz_S+ABY*I_ERI_G4z_S_Gx2yz_S;
  Double I_ERI_G4x_Pz_Gx2yz_S = I_ERI_H4xz_S_Gx2yz_S+ABZ*I_ERI_G4x_S_Gx2yz_S;
  Double I_ERI_G3xy_Pz_Gx2yz_S = I_ERI_H3xyz_S_Gx2yz_S+ABZ*I_ERI_G3xy_S_Gx2yz_S;
  Double I_ERI_G3xz_Pz_Gx2yz_S = I_ERI_H3x2z_S_Gx2yz_S+ABZ*I_ERI_G3xz_S_Gx2yz_S;
  Double I_ERI_G2x2y_Pz_Gx2yz_S = I_ERI_H2x2yz_S_Gx2yz_S+ABZ*I_ERI_G2x2y_S_Gx2yz_S;
  Double I_ERI_G2xyz_Pz_Gx2yz_S = I_ERI_H2xy2z_S_Gx2yz_S+ABZ*I_ERI_G2xyz_S_Gx2yz_S;
  Double I_ERI_G2x2z_Pz_Gx2yz_S = I_ERI_H2x3z_S_Gx2yz_S+ABZ*I_ERI_G2x2z_S_Gx2yz_S;
  Double I_ERI_Gx3y_Pz_Gx2yz_S = I_ERI_Hx3yz_S_Gx2yz_S+ABZ*I_ERI_Gx3y_S_Gx2yz_S;
  Double I_ERI_Gx2yz_Pz_Gx2yz_S = I_ERI_Hx2y2z_S_Gx2yz_S+ABZ*I_ERI_Gx2yz_S_Gx2yz_S;
  Double I_ERI_Gxy2z_Pz_Gx2yz_S = I_ERI_Hxy3z_S_Gx2yz_S+ABZ*I_ERI_Gxy2z_S_Gx2yz_S;
  Double I_ERI_Gx3z_Pz_Gx2yz_S = I_ERI_Hx4z_S_Gx2yz_S+ABZ*I_ERI_Gx3z_S_Gx2yz_S;
  Double I_ERI_G4y_Pz_Gx2yz_S = I_ERI_H4yz_S_Gx2yz_S+ABZ*I_ERI_G4y_S_Gx2yz_S;
  Double I_ERI_G3yz_Pz_Gx2yz_S = I_ERI_H3y2z_S_Gx2yz_S+ABZ*I_ERI_G3yz_S_Gx2yz_S;
  Double I_ERI_G2y2z_Pz_Gx2yz_S = I_ERI_H2y3z_S_Gx2yz_S+ABZ*I_ERI_G2y2z_S_Gx2yz_S;
  Double I_ERI_Gy3z_Pz_Gx2yz_S = I_ERI_Hy4z_S_Gx2yz_S+ABZ*I_ERI_Gy3z_S_Gx2yz_S;
  Double I_ERI_G4z_Pz_Gx2yz_S = I_ERI_H5z_S_Gx2yz_S+ABZ*I_ERI_G4z_S_Gx2yz_S;
  Double I_ERI_G4x_Px_Gxy2z_S = I_ERI_H5x_S_Gxy2z_S+ABX*I_ERI_G4x_S_Gxy2z_S;
  Double I_ERI_G3xy_Px_Gxy2z_S = I_ERI_H4xy_S_Gxy2z_S+ABX*I_ERI_G3xy_S_Gxy2z_S;
  Double I_ERI_G3xz_Px_Gxy2z_S = I_ERI_H4xz_S_Gxy2z_S+ABX*I_ERI_G3xz_S_Gxy2z_S;
  Double I_ERI_G2x2y_Px_Gxy2z_S = I_ERI_H3x2y_S_Gxy2z_S+ABX*I_ERI_G2x2y_S_Gxy2z_S;
  Double I_ERI_G2xyz_Px_Gxy2z_S = I_ERI_H3xyz_S_Gxy2z_S+ABX*I_ERI_G2xyz_S_Gxy2z_S;
  Double I_ERI_G2x2z_Px_Gxy2z_S = I_ERI_H3x2z_S_Gxy2z_S+ABX*I_ERI_G2x2z_S_Gxy2z_S;
  Double I_ERI_Gx3y_Px_Gxy2z_S = I_ERI_H2x3y_S_Gxy2z_S+ABX*I_ERI_Gx3y_S_Gxy2z_S;
  Double I_ERI_Gx2yz_Px_Gxy2z_S = I_ERI_H2x2yz_S_Gxy2z_S+ABX*I_ERI_Gx2yz_S_Gxy2z_S;
  Double I_ERI_Gxy2z_Px_Gxy2z_S = I_ERI_H2xy2z_S_Gxy2z_S+ABX*I_ERI_Gxy2z_S_Gxy2z_S;
  Double I_ERI_Gx3z_Px_Gxy2z_S = I_ERI_H2x3z_S_Gxy2z_S+ABX*I_ERI_Gx3z_S_Gxy2z_S;
  Double I_ERI_G4y_Px_Gxy2z_S = I_ERI_Hx4y_S_Gxy2z_S+ABX*I_ERI_G4y_S_Gxy2z_S;
  Double I_ERI_G3yz_Px_Gxy2z_S = I_ERI_Hx3yz_S_Gxy2z_S+ABX*I_ERI_G3yz_S_Gxy2z_S;
  Double I_ERI_G2y2z_Px_Gxy2z_S = I_ERI_Hx2y2z_S_Gxy2z_S+ABX*I_ERI_G2y2z_S_Gxy2z_S;
  Double I_ERI_Gy3z_Px_Gxy2z_S = I_ERI_Hxy3z_S_Gxy2z_S+ABX*I_ERI_Gy3z_S_Gxy2z_S;
  Double I_ERI_G4z_Px_Gxy2z_S = I_ERI_Hx4z_S_Gxy2z_S+ABX*I_ERI_G4z_S_Gxy2z_S;
  Double I_ERI_G4x_Py_Gxy2z_S = I_ERI_H4xy_S_Gxy2z_S+ABY*I_ERI_G4x_S_Gxy2z_S;
  Double I_ERI_G3xy_Py_Gxy2z_S = I_ERI_H3x2y_S_Gxy2z_S+ABY*I_ERI_G3xy_S_Gxy2z_S;
  Double I_ERI_G3xz_Py_Gxy2z_S = I_ERI_H3xyz_S_Gxy2z_S+ABY*I_ERI_G3xz_S_Gxy2z_S;
  Double I_ERI_G2x2y_Py_Gxy2z_S = I_ERI_H2x3y_S_Gxy2z_S+ABY*I_ERI_G2x2y_S_Gxy2z_S;
  Double I_ERI_G2xyz_Py_Gxy2z_S = I_ERI_H2x2yz_S_Gxy2z_S+ABY*I_ERI_G2xyz_S_Gxy2z_S;
  Double I_ERI_G2x2z_Py_Gxy2z_S = I_ERI_H2xy2z_S_Gxy2z_S+ABY*I_ERI_G2x2z_S_Gxy2z_S;
  Double I_ERI_Gx3y_Py_Gxy2z_S = I_ERI_Hx4y_S_Gxy2z_S+ABY*I_ERI_Gx3y_S_Gxy2z_S;
  Double I_ERI_Gx2yz_Py_Gxy2z_S = I_ERI_Hx3yz_S_Gxy2z_S+ABY*I_ERI_Gx2yz_S_Gxy2z_S;
  Double I_ERI_Gxy2z_Py_Gxy2z_S = I_ERI_Hx2y2z_S_Gxy2z_S+ABY*I_ERI_Gxy2z_S_Gxy2z_S;
  Double I_ERI_Gx3z_Py_Gxy2z_S = I_ERI_Hxy3z_S_Gxy2z_S+ABY*I_ERI_Gx3z_S_Gxy2z_S;
  Double I_ERI_G4y_Py_Gxy2z_S = I_ERI_H5y_S_Gxy2z_S+ABY*I_ERI_G4y_S_Gxy2z_S;
  Double I_ERI_G3yz_Py_Gxy2z_S = I_ERI_H4yz_S_Gxy2z_S+ABY*I_ERI_G3yz_S_Gxy2z_S;
  Double I_ERI_G2y2z_Py_Gxy2z_S = I_ERI_H3y2z_S_Gxy2z_S+ABY*I_ERI_G2y2z_S_Gxy2z_S;
  Double I_ERI_Gy3z_Py_Gxy2z_S = I_ERI_H2y3z_S_Gxy2z_S+ABY*I_ERI_Gy3z_S_Gxy2z_S;
  Double I_ERI_G4z_Py_Gxy2z_S = I_ERI_Hy4z_S_Gxy2z_S+ABY*I_ERI_G4z_S_Gxy2z_S;
  Double I_ERI_G4x_Pz_Gxy2z_S = I_ERI_H4xz_S_Gxy2z_S+ABZ*I_ERI_G4x_S_Gxy2z_S;
  Double I_ERI_G3xy_Pz_Gxy2z_S = I_ERI_H3xyz_S_Gxy2z_S+ABZ*I_ERI_G3xy_S_Gxy2z_S;
  Double I_ERI_G3xz_Pz_Gxy2z_S = I_ERI_H3x2z_S_Gxy2z_S+ABZ*I_ERI_G3xz_S_Gxy2z_S;
  Double I_ERI_G2x2y_Pz_Gxy2z_S = I_ERI_H2x2yz_S_Gxy2z_S+ABZ*I_ERI_G2x2y_S_Gxy2z_S;
  Double I_ERI_G2xyz_Pz_Gxy2z_S = I_ERI_H2xy2z_S_Gxy2z_S+ABZ*I_ERI_G2xyz_S_Gxy2z_S;
  Double I_ERI_G2x2z_Pz_Gxy2z_S = I_ERI_H2x3z_S_Gxy2z_S+ABZ*I_ERI_G2x2z_S_Gxy2z_S;
  Double I_ERI_Gx3y_Pz_Gxy2z_S = I_ERI_Hx3yz_S_Gxy2z_S+ABZ*I_ERI_Gx3y_S_Gxy2z_S;
  Double I_ERI_Gx2yz_Pz_Gxy2z_S = I_ERI_Hx2y2z_S_Gxy2z_S+ABZ*I_ERI_Gx2yz_S_Gxy2z_S;
  Double I_ERI_Gxy2z_Pz_Gxy2z_S = I_ERI_Hxy3z_S_Gxy2z_S+ABZ*I_ERI_Gxy2z_S_Gxy2z_S;
  Double I_ERI_Gx3z_Pz_Gxy2z_S = I_ERI_Hx4z_S_Gxy2z_S+ABZ*I_ERI_Gx3z_S_Gxy2z_S;
  Double I_ERI_G4y_Pz_Gxy2z_S = I_ERI_H4yz_S_Gxy2z_S+ABZ*I_ERI_G4y_S_Gxy2z_S;
  Double I_ERI_G3yz_Pz_Gxy2z_S = I_ERI_H3y2z_S_Gxy2z_S+ABZ*I_ERI_G3yz_S_Gxy2z_S;
  Double I_ERI_G2y2z_Pz_Gxy2z_S = I_ERI_H2y3z_S_Gxy2z_S+ABZ*I_ERI_G2y2z_S_Gxy2z_S;
  Double I_ERI_Gy3z_Pz_Gxy2z_S = I_ERI_Hy4z_S_Gxy2z_S+ABZ*I_ERI_Gy3z_S_Gxy2z_S;
  Double I_ERI_G4z_Pz_Gxy2z_S = I_ERI_H5z_S_Gxy2z_S+ABZ*I_ERI_G4z_S_Gxy2z_S;
  Double I_ERI_G4x_Px_Gx3z_S = I_ERI_H5x_S_Gx3z_S+ABX*I_ERI_G4x_S_Gx3z_S;
  Double I_ERI_G3xy_Px_Gx3z_S = I_ERI_H4xy_S_Gx3z_S+ABX*I_ERI_G3xy_S_Gx3z_S;
  Double I_ERI_G3xz_Px_Gx3z_S = I_ERI_H4xz_S_Gx3z_S+ABX*I_ERI_G3xz_S_Gx3z_S;
  Double I_ERI_G2x2y_Px_Gx3z_S = I_ERI_H3x2y_S_Gx3z_S+ABX*I_ERI_G2x2y_S_Gx3z_S;
  Double I_ERI_G2xyz_Px_Gx3z_S = I_ERI_H3xyz_S_Gx3z_S+ABX*I_ERI_G2xyz_S_Gx3z_S;
  Double I_ERI_G2x2z_Px_Gx3z_S = I_ERI_H3x2z_S_Gx3z_S+ABX*I_ERI_G2x2z_S_Gx3z_S;
  Double I_ERI_Gx3y_Px_Gx3z_S = I_ERI_H2x3y_S_Gx3z_S+ABX*I_ERI_Gx3y_S_Gx3z_S;
  Double I_ERI_Gx2yz_Px_Gx3z_S = I_ERI_H2x2yz_S_Gx3z_S+ABX*I_ERI_Gx2yz_S_Gx3z_S;
  Double I_ERI_Gxy2z_Px_Gx3z_S = I_ERI_H2xy2z_S_Gx3z_S+ABX*I_ERI_Gxy2z_S_Gx3z_S;
  Double I_ERI_Gx3z_Px_Gx3z_S = I_ERI_H2x3z_S_Gx3z_S+ABX*I_ERI_Gx3z_S_Gx3z_S;
  Double I_ERI_G4y_Px_Gx3z_S = I_ERI_Hx4y_S_Gx3z_S+ABX*I_ERI_G4y_S_Gx3z_S;
  Double I_ERI_G3yz_Px_Gx3z_S = I_ERI_Hx3yz_S_Gx3z_S+ABX*I_ERI_G3yz_S_Gx3z_S;
  Double I_ERI_G2y2z_Px_Gx3z_S = I_ERI_Hx2y2z_S_Gx3z_S+ABX*I_ERI_G2y2z_S_Gx3z_S;
  Double I_ERI_Gy3z_Px_Gx3z_S = I_ERI_Hxy3z_S_Gx3z_S+ABX*I_ERI_Gy3z_S_Gx3z_S;
  Double I_ERI_G4z_Px_Gx3z_S = I_ERI_Hx4z_S_Gx3z_S+ABX*I_ERI_G4z_S_Gx3z_S;
  Double I_ERI_G4x_Py_Gx3z_S = I_ERI_H4xy_S_Gx3z_S+ABY*I_ERI_G4x_S_Gx3z_S;
  Double I_ERI_G3xy_Py_Gx3z_S = I_ERI_H3x2y_S_Gx3z_S+ABY*I_ERI_G3xy_S_Gx3z_S;
  Double I_ERI_G3xz_Py_Gx3z_S = I_ERI_H3xyz_S_Gx3z_S+ABY*I_ERI_G3xz_S_Gx3z_S;
  Double I_ERI_G2x2y_Py_Gx3z_S = I_ERI_H2x3y_S_Gx3z_S+ABY*I_ERI_G2x2y_S_Gx3z_S;
  Double I_ERI_G2xyz_Py_Gx3z_S = I_ERI_H2x2yz_S_Gx3z_S+ABY*I_ERI_G2xyz_S_Gx3z_S;
  Double I_ERI_G2x2z_Py_Gx3z_S = I_ERI_H2xy2z_S_Gx3z_S+ABY*I_ERI_G2x2z_S_Gx3z_S;
  Double I_ERI_Gx3y_Py_Gx3z_S = I_ERI_Hx4y_S_Gx3z_S+ABY*I_ERI_Gx3y_S_Gx3z_S;
  Double I_ERI_Gx2yz_Py_Gx3z_S = I_ERI_Hx3yz_S_Gx3z_S+ABY*I_ERI_Gx2yz_S_Gx3z_S;
  Double I_ERI_Gxy2z_Py_Gx3z_S = I_ERI_Hx2y2z_S_Gx3z_S+ABY*I_ERI_Gxy2z_S_Gx3z_S;
  Double I_ERI_Gx3z_Py_Gx3z_S = I_ERI_Hxy3z_S_Gx3z_S+ABY*I_ERI_Gx3z_S_Gx3z_S;
  Double I_ERI_G4y_Py_Gx3z_S = I_ERI_H5y_S_Gx3z_S+ABY*I_ERI_G4y_S_Gx3z_S;
  Double I_ERI_G3yz_Py_Gx3z_S = I_ERI_H4yz_S_Gx3z_S+ABY*I_ERI_G3yz_S_Gx3z_S;
  Double I_ERI_G2y2z_Py_Gx3z_S = I_ERI_H3y2z_S_Gx3z_S+ABY*I_ERI_G2y2z_S_Gx3z_S;
  Double I_ERI_Gy3z_Py_Gx3z_S = I_ERI_H2y3z_S_Gx3z_S+ABY*I_ERI_Gy3z_S_Gx3z_S;
  Double I_ERI_G4z_Py_Gx3z_S = I_ERI_Hy4z_S_Gx3z_S+ABY*I_ERI_G4z_S_Gx3z_S;
  Double I_ERI_G4x_Pz_Gx3z_S = I_ERI_H4xz_S_Gx3z_S+ABZ*I_ERI_G4x_S_Gx3z_S;
  Double I_ERI_G3xy_Pz_Gx3z_S = I_ERI_H3xyz_S_Gx3z_S+ABZ*I_ERI_G3xy_S_Gx3z_S;
  Double I_ERI_G3xz_Pz_Gx3z_S = I_ERI_H3x2z_S_Gx3z_S+ABZ*I_ERI_G3xz_S_Gx3z_S;
  Double I_ERI_G2x2y_Pz_Gx3z_S = I_ERI_H2x2yz_S_Gx3z_S+ABZ*I_ERI_G2x2y_S_Gx3z_S;
  Double I_ERI_G2xyz_Pz_Gx3z_S = I_ERI_H2xy2z_S_Gx3z_S+ABZ*I_ERI_G2xyz_S_Gx3z_S;
  Double I_ERI_G2x2z_Pz_Gx3z_S = I_ERI_H2x3z_S_Gx3z_S+ABZ*I_ERI_G2x2z_S_Gx3z_S;
  Double I_ERI_Gx3y_Pz_Gx3z_S = I_ERI_Hx3yz_S_Gx3z_S+ABZ*I_ERI_Gx3y_S_Gx3z_S;
  Double I_ERI_Gx2yz_Pz_Gx3z_S = I_ERI_Hx2y2z_S_Gx3z_S+ABZ*I_ERI_Gx2yz_S_Gx3z_S;
  Double I_ERI_Gxy2z_Pz_Gx3z_S = I_ERI_Hxy3z_S_Gx3z_S+ABZ*I_ERI_Gxy2z_S_Gx3z_S;
  Double I_ERI_Gx3z_Pz_Gx3z_S = I_ERI_Hx4z_S_Gx3z_S+ABZ*I_ERI_Gx3z_S_Gx3z_S;
  Double I_ERI_G4y_Pz_Gx3z_S = I_ERI_H4yz_S_Gx3z_S+ABZ*I_ERI_G4y_S_Gx3z_S;
  Double I_ERI_G3yz_Pz_Gx3z_S = I_ERI_H3y2z_S_Gx3z_S+ABZ*I_ERI_G3yz_S_Gx3z_S;
  Double I_ERI_G2y2z_Pz_Gx3z_S = I_ERI_H2y3z_S_Gx3z_S+ABZ*I_ERI_G2y2z_S_Gx3z_S;
  Double I_ERI_Gy3z_Pz_Gx3z_S = I_ERI_Hy4z_S_Gx3z_S+ABZ*I_ERI_Gy3z_S_Gx3z_S;
  Double I_ERI_G4z_Pz_Gx3z_S = I_ERI_H5z_S_Gx3z_S+ABZ*I_ERI_G4z_S_Gx3z_S;
  Double I_ERI_G4x_Px_G4y_S = I_ERI_H5x_S_G4y_S+ABX*I_ERI_G4x_S_G4y_S;
  Double I_ERI_G3xy_Px_G4y_S = I_ERI_H4xy_S_G4y_S+ABX*I_ERI_G3xy_S_G4y_S;
  Double I_ERI_G3xz_Px_G4y_S = I_ERI_H4xz_S_G4y_S+ABX*I_ERI_G3xz_S_G4y_S;
  Double I_ERI_G2x2y_Px_G4y_S = I_ERI_H3x2y_S_G4y_S+ABX*I_ERI_G2x2y_S_G4y_S;
  Double I_ERI_G2xyz_Px_G4y_S = I_ERI_H3xyz_S_G4y_S+ABX*I_ERI_G2xyz_S_G4y_S;
  Double I_ERI_G2x2z_Px_G4y_S = I_ERI_H3x2z_S_G4y_S+ABX*I_ERI_G2x2z_S_G4y_S;
  Double I_ERI_Gx3y_Px_G4y_S = I_ERI_H2x3y_S_G4y_S+ABX*I_ERI_Gx3y_S_G4y_S;
  Double I_ERI_Gx2yz_Px_G4y_S = I_ERI_H2x2yz_S_G4y_S+ABX*I_ERI_Gx2yz_S_G4y_S;
  Double I_ERI_Gxy2z_Px_G4y_S = I_ERI_H2xy2z_S_G4y_S+ABX*I_ERI_Gxy2z_S_G4y_S;
  Double I_ERI_Gx3z_Px_G4y_S = I_ERI_H2x3z_S_G4y_S+ABX*I_ERI_Gx3z_S_G4y_S;
  Double I_ERI_G4y_Px_G4y_S = I_ERI_Hx4y_S_G4y_S+ABX*I_ERI_G4y_S_G4y_S;
  Double I_ERI_G3yz_Px_G4y_S = I_ERI_Hx3yz_S_G4y_S+ABX*I_ERI_G3yz_S_G4y_S;
  Double I_ERI_G2y2z_Px_G4y_S = I_ERI_Hx2y2z_S_G4y_S+ABX*I_ERI_G2y2z_S_G4y_S;
  Double I_ERI_Gy3z_Px_G4y_S = I_ERI_Hxy3z_S_G4y_S+ABX*I_ERI_Gy3z_S_G4y_S;
  Double I_ERI_G4z_Px_G4y_S = I_ERI_Hx4z_S_G4y_S+ABX*I_ERI_G4z_S_G4y_S;
  Double I_ERI_G4x_Py_G4y_S = I_ERI_H4xy_S_G4y_S+ABY*I_ERI_G4x_S_G4y_S;
  Double I_ERI_G3xy_Py_G4y_S = I_ERI_H3x2y_S_G4y_S+ABY*I_ERI_G3xy_S_G4y_S;
  Double I_ERI_G3xz_Py_G4y_S = I_ERI_H3xyz_S_G4y_S+ABY*I_ERI_G3xz_S_G4y_S;
  Double I_ERI_G2x2y_Py_G4y_S = I_ERI_H2x3y_S_G4y_S+ABY*I_ERI_G2x2y_S_G4y_S;
  Double I_ERI_G2xyz_Py_G4y_S = I_ERI_H2x2yz_S_G4y_S+ABY*I_ERI_G2xyz_S_G4y_S;
  Double I_ERI_G2x2z_Py_G4y_S = I_ERI_H2xy2z_S_G4y_S+ABY*I_ERI_G2x2z_S_G4y_S;
  Double I_ERI_Gx3y_Py_G4y_S = I_ERI_Hx4y_S_G4y_S+ABY*I_ERI_Gx3y_S_G4y_S;
  Double I_ERI_Gx2yz_Py_G4y_S = I_ERI_Hx3yz_S_G4y_S+ABY*I_ERI_Gx2yz_S_G4y_S;
  Double I_ERI_Gxy2z_Py_G4y_S = I_ERI_Hx2y2z_S_G4y_S+ABY*I_ERI_Gxy2z_S_G4y_S;
  Double I_ERI_Gx3z_Py_G4y_S = I_ERI_Hxy3z_S_G4y_S+ABY*I_ERI_Gx3z_S_G4y_S;
  Double I_ERI_G4y_Py_G4y_S = I_ERI_H5y_S_G4y_S+ABY*I_ERI_G4y_S_G4y_S;
  Double I_ERI_G3yz_Py_G4y_S = I_ERI_H4yz_S_G4y_S+ABY*I_ERI_G3yz_S_G4y_S;
  Double I_ERI_G2y2z_Py_G4y_S = I_ERI_H3y2z_S_G4y_S+ABY*I_ERI_G2y2z_S_G4y_S;
  Double I_ERI_Gy3z_Py_G4y_S = I_ERI_H2y3z_S_G4y_S+ABY*I_ERI_Gy3z_S_G4y_S;
  Double I_ERI_G4z_Py_G4y_S = I_ERI_Hy4z_S_G4y_S+ABY*I_ERI_G4z_S_G4y_S;
  Double I_ERI_G4x_Pz_G4y_S = I_ERI_H4xz_S_G4y_S+ABZ*I_ERI_G4x_S_G4y_S;
  Double I_ERI_G3xy_Pz_G4y_S = I_ERI_H3xyz_S_G4y_S+ABZ*I_ERI_G3xy_S_G4y_S;
  Double I_ERI_G3xz_Pz_G4y_S = I_ERI_H3x2z_S_G4y_S+ABZ*I_ERI_G3xz_S_G4y_S;
  Double I_ERI_G2x2y_Pz_G4y_S = I_ERI_H2x2yz_S_G4y_S+ABZ*I_ERI_G2x2y_S_G4y_S;
  Double I_ERI_G2xyz_Pz_G4y_S = I_ERI_H2xy2z_S_G4y_S+ABZ*I_ERI_G2xyz_S_G4y_S;
  Double I_ERI_G2x2z_Pz_G4y_S = I_ERI_H2x3z_S_G4y_S+ABZ*I_ERI_G2x2z_S_G4y_S;
  Double I_ERI_Gx3y_Pz_G4y_S = I_ERI_Hx3yz_S_G4y_S+ABZ*I_ERI_Gx3y_S_G4y_S;
  Double I_ERI_Gx2yz_Pz_G4y_S = I_ERI_Hx2y2z_S_G4y_S+ABZ*I_ERI_Gx2yz_S_G4y_S;
  Double I_ERI_Gxy2z_Pz_G4y_S = I_ERI_Hxy3z_S_G4y_S+ABZ*I_ERI_Gxy2z_S_G4y_S;
  Double I_ERI_Gx3z_Pz_G4y_S = I_ERI_Hx4z_S_G4y_S+ABZ*I_ERI_Gx3z_S_G4y_S;
  Double I_ERI_G4y_Pz_G4y_S = I_ERI_H4yz_S_G4y_S+ABZ*I_ERI_G4y_S_G4y_S;
  Double I_ERI_G3yz_Pz_G4y_S = I_ERI_H3y2z_S_G4y_S+ABZ*I_ERI_G3yz_S_G4y_S;
  Double I_ERI_G2y2z_Pz_G4y_S = I_ERI_H2y3z_S_G4y_S+ABZ*I_ERI_G2y2z_S_G4y_S;
  Double I_ERI_Gy3z_Pz_G4y_S = I_ERI_Hy4z_S_G4y_S+ABZ*I_ERI_Gy3z_S_G4y_S;
  Double I_ERI_G4z_Pz_G4y_S = I_ERI_H5z_S_G4y_S+ABZ*I_ERI_G4z_S_G4y_S;
  Double I_ERI_G4x_Px_G3yz_S = I_ERI_H5x_S_G3yz_S+ABX*I_ERI_G4x_S_G3yz_S;
  Double I_ERI_G3xy_Px_G3yz_S = I_ERI_H4xy_S_G3yz_S+ABX*I_ERI_G3xy_S_G3yz_S;
  Double I_ERI_G3xz_Px_G3yz_S = I_ERI_H4xz_S_G3yz_S+ABX*I_ERI_G3xz_S_G3yz_S;
  Double I_ERI_G2x2y_Px_G3yz_S = I_ERI_H3x2y_S_G3yz_S+ABX*I_ERI_G2x2y_S_G3yz_S;
  Double I_ERI_G2xyz_Px_G3yz_S = I_ERI_H3xyz_S_G3yz_S+ABX*I_ERI_G2xyz_S_G3yz_S;
  Double I_ERI_G2x2z_Px_G3yz_S = I_ERI_H3x2z_S_G3yz_S+ABX*I_ERI_G2x2z_S_G3yz_S;
  Double I_ERI_Gx3y_Px_G3yz_S = I_ERI_H2x3y_S_G3yz_S+ABX*I_ERI_Gx3y_S_G3yz_S;
  Double I_ERI_Gx2yz_Px_G3yz_S = I_ERI_H2x2yz_S_G3yz_S+ABX*I_ERI_Gx2yz_S_G3yz_S;
  Double I_ERI_Gxy2z_Px_G3yz_S = I_ERI_H2xy2z_S_G3yz_S+ABX*I_ERI_Gxy2z_S_G3yz_S;
  Double I_ERI_Gx3z_Px_G3yz_S = I_ERI_H2x3z_S_G3yz_S+ABX*I_ERI_Gx3z_S_G3yz_S;
  Double I_ERI_G4y_Px_G3yz_S = I_ERI_Hx4y_S_G3yz_S+ABX*I_ERI_G4y_S_G3yz_S;
  Double I_ERI_G3yz_Px_G3yz_S = I_ERI_Hx3yz_S_G3yz_S+ABX*I_ERI_G3yz_S_G3yz_S;
  Double I_ERI_G2y2z_Px_G3yz_S = I_ERI_Hx2y2z_S_G3yz_S+ABX*I_ERI_G2y2z_S_G3yz_S;
  Double I_ERI_Gy3z_Px_G3yz_S = I_ERI_Hxy3z_S_G3yz_S+ABX*I_ERI_Gy3z_S_G3yz_S;
  Double I_ERI_G4z_Px_G3yz_S = I_ERI_Hx4z_S_G3yz_S+ABX*I_ERI_G4z_S_G3yz_S;
  Double I_ERI_G4x_Py_G3yz_S = I_ERI_H4xy_S_G3yz_S+ABY*I_ERI_G4x_S_G3yz_S;
  Double I_ERI_G3xy_Py_G3yz_S = I_ERI_H3x2y_S_G3yz_S+ABY*I_ERI_G3xy_S_G3yz_S;
  Double I_ERI_G3xz_Py_G3yz_S = I_ERI_H3xyz_S_G3yz_S+ABY*I_ERI_G3xz_S_G3yz_S;
  Double I_ERI_G2x2y_Py_G3yz_S = I_ERI_H2x3y_S_G3yz_S+ABY*I_ERI_G2x2y_S_G3yz_S;
  Double I_ERI_G2xyz_Py_G3yz_S = I_ERI_H2x2yz_S_G3yz_S+ABY*I_ERI_G2xyz_S_G3yz_S;
  Double I_ERI_G2x2z_Py_G3yz_S = I_ERI_H2xy2z_S_G3yz_S+ABY*I_ERI_G2x2z_S_G3yz_S;
  Double I_ERI_Gx3y_Py_G3yz_S = I_ERI_Hx4y_S_G3yz_S+ABY*I_ERI_Gx3y_S_G3yz_S;
  Double I_ERI_Gx2yz_Py_G3yz_S = I_ERI_Hx3yz_S_G3yz_S+ABY*I_ERI_Gx2yz_S_G3yz_S;
  Double I_ERI_Gxy2z_Py_G3yz_S = I_ERI_Hx2y2z_S_G3yz_S+ABY*I_ERI_Gxy2z_S_G3yz_S;
  Double I_ERI_Gx3z_Py_G3yz_S = I_ERI_Hxy3z_S_G3yz_S+ABY*I_ERI_Gx3z_S_G3yz_S;
  Double I_ERI_G4y_Py_G3yz_S = I_ERI_H5y_S_G3yz_S+ABY*I_ERI_G4y_S_G3yz_S;
  Double I_ERI_G3yz_Py_G3yz_S = I_ERI_H4yz_S_G3yz_S+ABY*I_ERI_G3yz_S_G3yz_S;
  Double I_ERI_G2y2z_Py_G3yz_S = I_ERI_H3y2z_S_G3yz_S+ABY*I_ERI_G2y2z_S_G3yz_S;
  Double I_ERI_Gy3z_Py_G3yz_S = I_ERI_H2y3z_S_G3yz_S+ABY*I_ERI_Gy3z_S_G3yz_S;
  Double I_ERI_G4z_Py_G3yz_S = I_ERI_Hy4z_S_G3yz_S+ABY*I_ERI_G4z_S_G3yz_S;
  Double I_ERI_G4x_Pz_G3yz_S = I_ERI_H4xz_S_G3yz_S+ABZ*I_ERI_G4x_S_G3yz_S;
  Double I_ERI_G3xy_Pz_G3yz_S = I_ERI_H3xyz_S_G3yz_S+ABZ*I_ERI_G3xy_S_G3yz_S;
  Double I_ERI_G3xz_Pz_G3yz_S = I_ERI_H3x2z_S_G3yz_S+ABZ*I_ERI_G3xz_S_G3yz_S;
  Double I_ERI_G2x2y_Pz_G3yz_S = I_ERI_H2x2yz_S_G3yz_S+ABZ*I_ERI_G2x2y_S_G3yz_S;
  Double I_ERI_G2xyz_Pz_G3yz_S = I_ERI_H2xy2z_S_G3yz_S+ABZ*I_ERI_G2xyz_S_G3yz_S;
  Double I_ERI_G2x2z_Pz_G3yz_S = I_ERI_H2x3z_S_G3yz_S+ABZ*I_ERI_G2x2z_S_G3yz_S;
  Double I_ERI_Gx3y_Pz_G3yz_S = I_ERI_Hx3yz_S_G3yz_S+ABZ*I_ERI_Gx3y_S_G3yz_S;
  Double I_ERI_Gx2yz_Pz_G3yz_S = I_ERI_Hx2y2z_S_G3yz_S+ABZ*I_ERI_Gx2yz_S_G3yz_S;
  Double I_ERI_Gxy2z_Pz_G3yz_S = I_ERI_Hxy3z_S_G3yz_S+ABZ*I_ERI_Gxy2z_S_G3yz_S;
  Double I_ERI_Gx3z_Pz_G3yz_S = I_ERI_Hx4z_S_G3yz_S+ABZ*I_ERI_Gx3z_S_G3yz_S;
  Double I_ERI_G4y_Pz_G3yz_S = I_ERI_H4yz_S_G3yz_S+ABZ*I_ERI_G4y_S_G3yz_S;
  Double I_ERI_G3yz_Pz_G3yz_S = I_ERI_H3y2z_S_G3yz_S+ABZ*I_ERI_G3yz_S_G3yz_S;
  Double I_ERI_G2y2z_Pz_G3yz_S = I_ERI_H2y3z_S_G3yz_S+ABZ*I_ERI_G2y2z_S_G3yz_S;
  Double I_ERI_Gy3z_Pz_G3yz_S = I_ERI_Hy4z_S_G3yz_S+ABZ*I_ERI_Gy3z_S_G3yz_S;
  Double I_ERI_G4z_Pz_G3yz_S = I_ERI_H5z_S_G3yz_S+ABZ*I_ERI_G4z_S_G3yz_S;
  Double I_ERI_G4x_Px_G2y2z_S = I_ERI_H5x_S_G2y2z_S+ABX*I_ERI_G4x_S_G2y2z_S;
  Double I_ERI_G3xy_Px_G2y2z_S = I_ERI_H4xy_S_G2y2z_S+ABX*I_ERI_G3xy_S_G2y2z_S;
  Double I_ERI_G3xz_Px_G2y2z_S = I_ERI_H4xz_S_G2y2z_S+ABX*I_ERI_G3xz_S_G2y2z_S;
  Double I_ERI_G2x2y_Px_G2y2z_S = I_ERI_H3x2y_S_G2y2z_S+ABX*I_ERI_G2x2y_S_G2y2z_S;
  Double I_ERI_G2xyz_Px_G2y2z_S = I_ERI_H3xyz_S_G2y2z_S+ABX*I_ERI_G2xyz_S_G2y2z_S;
  Double I_ERI_G2x2z_Px_G2y2z_S = I_ERI_H3x2z_S_G2y2z_S+ABX*I_ERI_G2x2z_S_G2y2z_S;
  Double I_ERI_Gx3y_Px_G2y2z_S = I_ERI_H2x3y_S_G2y2z_S+ABX*I_ERI_Gx3y_S_G2y2z_S;
  Double I_ERI_Gx2yz_Px_G2y2z_S = I_ERI_H2x2yz_S_G2y2z_S+ABX*I_ERI_Gx2yz_S_G2y2z_S;
  Double I_ERI_Gxy2z_Px_G2y2z_S = I_ERI_H2xy2z_S_G2y2z_S+ABX*I_ERI_Gxy2z_S_G2y2z_S;
  Double I_ERI_Gx3z_Px_G2y2z_S = I_ERI_H2x3z_S_G2y2z_S+ABX*I_ERI_Gx3z_S_G2y2z_S;
  Double I_ERI_G4y_Px_G2y2z_S = I_ERI_Hx4y_S_G2y2z_S+ABX*I_ERI_G4y_S_G2y2z_S;
  Double I_ERI_G3yz_Px_G2y2z_S = I_ERI_Hx3yz_S_G2y2z_S+ABX*I_ERI_G3yz_S_G2y2z_S;
  Double I_ERI_G2y2z_Px_G2y2z_S = I_ERI_Hx2y2z_S_G2y2z_S+ABX*I_ERI_G2y2z_S_G2y2z_S;
  Double I_ERI_Gy3z_Px_G2y2z_S = I_ERI_Hxy3z_S_G2y2z_S+ABX*I_ERI_Gy3z_S_G2y2z_S;
  Double I_ERI_G4z_Px_G2y2z_S = I_ERI_Hx4z_S_G2y2z_S+ABX*I_ERI_G4z_S_G2y2z_S;
  Double I_ERI_G4x_Py_G2y2z_S = I_ERI_H4xy_S_G2y2z_S+ABY*I_ERI_G4x_S_G2y2z_S;
  Double I_ERI_G3xy_Py_G2y2z_S = I_ERI_H3x2y_S_G2y2z_S+ABY*I_ERI_G3xy_S_G2y2z_S;
  Double I_ERI_G3xz_Py_G2y2z_S = I_ERI_H3xyz_S_G2y2z_S+ABY*I_ERI_G3xz_S_G2y2z_S;
  Double I_ERI_G2x2y_Py_G2y2z_S = I_ERI_H2x3y_S_G2y2z_S+ABY*I_ERI_G2x2y_S_G2y2z_S;
  Double I_ERI_G2xyz_Py_G2y2z_S = I_ERI_H2x2yz_S_G2y2z_S+ABY*I_ERI_G2xyz_S_G2y2z_S;
  Double I_ERI_G2x2z_Py_G2y2z_S = I_ERI_H2xy2z_S_G2y2z_S+ABY*I_ERI_G2x2z_S_G2y2z_S;
  Double I_ERI_Gx3y_Py_G2y2z_S = I_ERI_Hx4y_S_G2y2z_S+ABY*I_ERI_Gx3y_S_G2y2z_S;
  Double I_ERI_Gx2yz_Py_G2y2z_S = I_ERI_Hx3yz_S_G2y2z_S+ABY*I_ERI_Gx2yz_S_G2y2z_S;
  Double I_ERI_Gxy2z_Py_G2y2z_S = I_ERI_Hx2y2z_S_G2y2z_S+ABY*I_ERI_Gxy2z_S_G2y2z_S;
  Double I_ERI_Gx3z_Py_G2y2z_S = I_ERI_Hxy3z_S_G2y2z_S+ABY*I_ERI_Gx3z_S_G2y2z_S;
  Double I_ERI_G4y_Py_G2y2z_S = I_ERI_H5y_S_G2y2z_S+ABY*I_ERI_G4y_S_G2y2z_S;
  Double I_ERI_G3yz_Py_G2y2z_S = I_ERI_H4yz_S_G2y2z_S+ABY*I_ERI_G3yz_S_G2y2z_S;
  Double I_ERI_G2y2z_Py_G2y2z_S = I_ERI_H3y2z_S_G2y2z_S+ABY*I_ERI_G2y2z_S_G2y2z_S;
  Double I_ERI_Gy3z_Py_G2y2z_S = I_ERI_H2y3z_S_G2y2z_S+ABY*I_ERI_Gy3z_S_G2y2z_S;
  Double I_ERI_G4z_Py_G2y2z_S = I_ERI_Hy4z_S_G2y2z_S+ABY*I_ERI_G4z_S_G2y2z_S;
  Double I_ERI_G4x_Pz_G2y2z_S = I_ERI_H4xz_S_G2y2z_S+ABZ*I_ERI_G4x_S_G2y2z_S;
  Double I_ERI_G3xy_Pz_G2y2z_S = I_ERI_H3xyz_S_G2y2z_S+ABZ*I_ERI_G3xy_S_G2y2z_S;
  Double I_ERI_G3xz_Pz_G2y2z_S = I_ERI_H3x2z_S_G2y2z_S+ABZ*I_ERI_G3xz_S_G2y2z_S;
  Double I_ERI_G2x2y_Pz_G2y2z_S = I_ERI_H2x2yz_S_G2y2z_S+ABZ*I_ERI_G2x2y_S_G2y2z_S;
  Double I_ERI_G2xyz_Pz_G2y2z_S = I_ERI_H2xy2z_S_G2y2z_S+ABZ*I_ERI_G2xyz_S_G2y2z_S;
  Double I_ERI_G2x2z_Pz_G2y2z_S = I_ERI_H2x3z_S_G2y2z_S+ABZ*I_ERI_G2x2z_S_G2y2z_S;
  Double I_ERI_Gx3y_Pz_G2y2z_S = I_ERI_Hx3yz_S_G2y2z_S+ABZ*I_ERI_Gx3y_S_G2y2z_S;
  Double I_ERI_Gx2yz_Pz_G2y2z_S = I_ERI_Hx2y2z_S_G2y2z_S+ABZ*I_ERI_Gx2yz_S_G2y2z_S;
  Double I_ERI_Gxy2z_Pz_G2y2z_S = I_ERI_Hxy3z_S_G2y2z_S+ABZ*I_ERI_Gxy2z_S_G2y2z_S;
  Double I_ERI_Gx3z_Pz_G2y2z_S = I_ERI_Hx4z_S_G2y2z_S+ABZ*I_ERI_Gx3z_S_G2y2z_S;
  Double I_ERI_G4y_Pz_G2y2z_S = I_ERI_H4yz_S_G2y2z_S+ABZ*I_ERI_G4y_S_G2y2z_S;
  Double I_ERI_G3yz_Pz_G2y2z_S = I_ERI_H3y2z_S_G2y2z_S+ABZ*I_ERI_G3yz_S_G2y2z_S;
  Double I_ERI_G2y2z_Pz_G2y2z_S = I_ERI_H2y3z_S_G2y2z_S+ABZ*I_ERI_G2y2z_S_G2y2z_S;
  Double I_ERI_Gy3z_Pz_G2y2z_S = I_ERI_Hy4z_S_G2y2z_S+ABZ*I_ERI_Gy3z_S_G2y2z_S;
  Double I_ERI_G4z_Pz_G2y2z_S = I_ERI_H5z_S_G2y2z_S+ABZ*I_ERI_G4z_S_G2y2z_S;
  Double I_ERI_G4x_Px_Gy3z_S = I_ERI_H5x_S_Gy3z_S+ABX*I_ERI_G4x_S_Gy3z_S;
  Double I_ERI_G3xy_Px_Gy3z_S = I_ERI_H4xy_S_Gy3z_S+ABX*I_ERI_G3xy_S_Gy3z_S;
  Double I_ERI_G3xz_Px_Gy3z_S = I_ERI_H4xz_S_Gy3z_S+ABX*I_ERI_G3xz_S_Gy3z_S;
  Double I_ERI_G2x2y_Px_Gy3z_S = I_ERI_H3x2y_S_Gy3z_S+ABX*I_ERI_G2x2y_S_Gy3z_S;
  Double I_ERI_G2xyz_Px_Gy3z_S = I_ERI_H3xyz_S_Gy3z_S+ABX*I_ERI_G2xyz_S_Gy3z_S;
  Double I_ERI_G2x2z_Px_Gy3z_S = I_ERI_H3x2z_S_Gy3z_S+ABX*I_ERI_G2x2z_S_Gy3z_S;
  Double I_ERI_Gx3y_Px_Gy3z_S = I_ERI_H2x3y_S_Gy3z_S+ABX*I_ERI_Gx3y_S_Gy3z_S;
  Double I_ERI_Gx2yz_Px_Gy3z_S = I_ERI_H2x2yz_S_Gy3z_S+ABX*I_ERI_Gx2yz_S_Gy3z_S;
  Double I_ERI_Gxy2z_Px_Gy3z_S = I_ERI_H2xy2z_S_Gy3z_S+ABX*I_ERI_Gxy2z_S_Gy3z_S;
  Double I_ERI_Gx3z_Px_Gy3z_S = I_ERI_H2x3z_S_Gy3z_S+ABX*I_ERI_Gx3z_S_Gy3z_S;
  Double I_ERI_G4y_Px_Gy3z_S = I_ERI_Hx4y_S_Gy3z_S+ABX*I_ERI_G4y_S_Gy3z_S;
  Double I_ERI_G3yz_Px_Gy3z_S = I_ERI_Hx3yz_S_Gy3z_S+ABX*I_ERI_G3yz_S_Gy3z_S;
  Double I_ERI_G2y2z_Px_Gy3z_S = I_ERI_Hx2y2z_S_Gy3z_S+ABX*I_ERI_G2y2z_S_Gy3z_S;
  Double I_ERI_Gy3z_Px_Gy3z_S = I_ERI_Hxy3z_S_Gy3z_S+ABX*I_ERI_Gy3z_S_Gy3z_S;
  Double I_ERI_G4z_Px_Gy3z_S = I_ERI_Hx4z_S_Gy3z_S+ABX*I_ERI_G4z_S_Gy3z_S;
  Double I_ERI_G4x_Py_Gy3z_S = I_ERI_H4xy_S_Gy3z_S+ABY*I_ERI_G4x_S_Gy3z_S;
  Double I_ERI_G3xy_Py_Gy3z_S = I_ERI_H3x2y_S_Gy3z_S+ABY*I_ERI_G3xy_S_Gy3z_S;
  Double I_ERI_G3xz_Py_Gy3z_S = I_ERI_H3xyz_S_Gy3z_S+ABY*I_ERI_G3xz_S_Gy3z_S;
  Double I_ERI_G2x2y_Py_Gy3z_S = I_ERI_H2x3y_S_Gy3z_S+ABY*I_ERI_G2x2y_S_Gy3z_S;
  Double I_ERI_G2xyz_Py_Gy3z_S = I_ERI_H2x2yz_S_Gy3z_S+ABY*I_ERI_G2xyz_S_Gy3z_S;
  Double I_ERI_G2x2z_Py_Gy3z_S = I_ERI_H2xy2z_S_Gy3z_S+ABY*I_ERI_G2x2z_S_Gy3z_S;
  Double I_ERI_Gx3y_Py_Gy3z_S = I_ERI_Hx4y_S_Gy3z_S+ABY*I_ERI_Gx3y_S_Gy3z_S;
  Double I_ERI_Gx2yz_Py_Gy3z_S = I_ERI_Hx3yz_S_Gy3z_S+ABY*I_ERI_Gx2yz_S_Gy3z_S;
  Double I_ERI_Gxy2z_Py_Gy3z_S = I_ERI_Hx2y2z_S_Gy3z_S+ABY*I_ERI_Gxy2z_S_Gy3z_S;
  Double I_ERI_Gx3z_Py_Gy3z_S = I_ERI_Hxy3z_S_Gy3z_S+ABY*I_ERI_Gx3z_S_Gy3z_S;
  Double I_ERI_G4y_Py_Gy3z_S = I_ERI_H5y_S_Gy3z_S+ABY*I_ERI_G4y_S_Gy3z_S;
  Double I_ERI_G3yz_Py_Gy3z_S = I_ERI_H4yz_S_Gy3z_S+ABY*I_ERI_G3yz_S_Gy3z_S;
  Double I_ERI_G2y2z_Py_Gy3z_S = I_ERI_H3y2z_S_Gy3z_S+ABY*I_ERI_G2y2z_S_Gy3z_S;
  Double I_ERI_Gy3z_Py_Gy3z_S = I_ERI_H2y3z_S_Gy3z_S+ABY*I_ERI_Gy3z_S_Gy3z_S;
  Double I_ERI_G4z_Py_Gy3z_S = I_ERI_Hy4z_S_Gy3z_S+ABY*I_ERI_G4z_S_Gy3z_S;
  Double I_ERI_G4x_Pz_Gy3z_S = I_ERI_H4xz_S_Gy3z_S+ABZ*I_ERI_G4x_S_Gy3z_S;
  Double I_ERI_G3xy_Pz_Gy3z_S = I_ERI_H3xyz_S_Gy3z_S+ABZ*I_ERI_G3xy_S_Gy3z_S;
  Double I_ERI_G3xz_Pz_Gy3z_S = I_ERI_H3x2z_S_Gy3z_S+ABZ*I_ERI_G3xz_S_Gy3z_S;
  Double I_ERI_G2x2y_Pz_Gy3z_S = I_ERI_H2x2yz_S_Gy3z_S+ABZ*I_ERI_G2x2y_S_Gy3z_S;
  Double I_ERI_G2xyz_Pz_Gy3z_S = I_ERI_H2xy2z_S_Gy3z_S+ABZ*I_ERI_G2xyz_S_Gy3z_S;
  Double I_ERI_G2x2z_Pz_Gy3z_S = I_ERI_H2x3z_S_Gy3z_S+ABZ*I_ERI_G2x2z_S_Gy3z_S;
  Double I_ERI_Gx3y_Pz_Gy3z_S = I_ERI_Hx3yz_S_Gy3z_S+ABZ*I_ERI_Gx3y_S_Gy3z_S;
  Double I_ERI_Gx2yz_Pz_Gy3z_S = I_ERI_Hx2y2z_S_Gy3z_S+ABZ*I_ERI_Gx2yz_S_Gy3z_S;
  Double I_ERI_Gxy2z_Pz_Gy3z_S = I_ERI_Hxy3z_S_Gy3z_S+ABZ*I_ERI_Gxy2z_S_Gy3z_S;
  Double I_ERI_Gx3z_Pz_Gy3z_S = I_ERI_Hx4z_S_Gy3z_S+ABZ*I_ERI_Gx3z_S_Gy3z_S;
  Double I_ERI_G4y_Pz_Gy3z_S = I_ERI_H4yz_S_Gy3z_S+ABZ*I_ERI_G4y_S_Gy3z_S;
  Double I_ERI_G3yz_Pz_Gy3z_S = I_ERI_H3y2z_S_Gy3z_S+ABZ*I_ERI_G3yz_S_Gy3z_S;
  Double I_ERI_G2y2z_Pz_Gy3z_S = I_ERI_H2y3z_S_Gy3z_S+ABZ*I_ERI_G2y2z_S_Gy3z_S;
  Double I_ERI_Gy3z_Pz_Gy3z_S = I_ERI_Hy4z_S_Gy3z_S+ABZ*I_ERI_Gy3z_S_Gy3z_S;
  Double I_ERI_G4z_Pz_Gy3z_S = I_ERI_H5z_S_Gy3z_S+ABZ*I_ERI_G4z_S_Gy3z_S;
  Double I_ERI_G4x_Px_G4z_S = I_ERI_H5x_S_G4z_S+ABX*I_ERI_G4x_S_G4z_S;
  Double I_ERI_G3xy_Px_G4z_S = I_ERI_H4xy_S_G4z_S+ABX*I_ERI_G3xy_S_G4z_S;
  Double I_ERI_G3xz_Px_G4z_S = I_ERI_H4xz_S_G4z_S+ABX*I_ERI_G3xz_S_G4z_S;
  Double I_ERI_G2x2y_Px_G4z_S = I_ERI_H3x2y_S_G4z_S+ABX*I_ERI_G2x2y_S_G4z_S;
  Double I_ERI_G2xyz_Px_G4z_S = I_ERI_H3xyz_S_G4z_S+ABX*I_ERI_G2xyz_S_G4z_S;
  Double I_ERI_G2x2z_Px_G4z_S = I_ERI_H3x2z_S_G4z_S+ABX*I_ERI_G2x2z_S_G4z_S;
  Double I_ERI_Gx3y_Px_G4z_S = I_ERI_H2x3y_S_G4z_S+ABX*I_ERI_Gx3y_S_G4z_S;
  Double I_ERI_Gx2yz_Px_G4z_S = I_ERI_H2x2yz_S_G4z_S+ABX*I_ERI_Gx2yz_S_G4z_S;
  Double I_ERI_Gxy2z_Px_G4z_S = I_ERI_H2xy2z_S_G4z_S+ABX*I_ERI_Gxy2z_S_G4z_S;
  Double I_ERI_Gx3z_Px_G4z_S = I_ERI_H2x3z_S_G4z_S+ABX*I_ERI_Gx3z_S_G4z_S;
  Double I_ERI_G4y_Px_G4z_S = I_ERI_Hx4y_S_G4z_S+ABX*I_ERI_G4y_S_G4z_S;
  Double I_ERI_G3yz_Px_G4z_S = I_ERI_Hx3yz_S_G4z_S+ABX*I_ERI_G3yz_S_G4z_S;
  Double I_ERI_G2y2z_Px_G4z_S = I_ERI_Hx2y2z_S_G4z_S+ABX*I_ERI_G2y2z_S_G4z_S;
  Double I_ERI_Gy3z_Px_G4z_S = I_ERI_Hxy3z_S_G4z_S+ABX*I_ERI_Gy3z_S_G4z_S;
  Double I_ERI_G4z_Px_G4z_S = I_ERI_Hx4z_S_G4z_S+ABX*I_ERI_G4z_S_G4z_S;
  Double I_ERI_G4x_Py_G4z_S = I_ERI_H4xy_S_G4z_S+ABY*I_ERI_G4x_S_G4z_S;
  Double I_ERI_G3xy_Py_G4z_S = I_ERI_H3x2y_S_G4z_S+ABY*I_ERI_G3xy_S_G4z_S;
  Double I_ERI_G3xz_Py_G4z_S = I_ERI_H3xyz_S_G4z_S+ABY*I_ERI_G3xz_S_G4z_S;
  Double I_ERI_G2x2y_Py_G4z_S = I_ERI_H2x3y_S_G4z_S+ABY*I_ERI_G2x2y_S_G4z_S;
  Double I_ERI_G2xyz_Py_G4z_S = I_ERI_H2x2yz_S_G4z_S+ABY*I_ERI_G2xyz_S_G4z_S;
  Double I_ERI_G2x2z_Py_G4z_S = I_ERI_H2xy2z_S_G4z_S+ABY*I_ERI_G2x2z_S_G4z_S;
  Double I_ERI_Gx3y_Py_G4z_S = I_ERI_Hx4y_S_G4z_S+ABY*I_ERI_Gx3y_S_G4z_S;
  Double I_ERI_Gx2yz_Py_G4z_S = I_ERI_Hx3yz_S_G4z_S+ABY*I_ERI_Gx2yz_S_G4z_S;
  Double I_ERI_Gxy2z_Py_G4z_S = I_ERI_Hx2y2z_S_G4z_S+ABY*I_ERI_Gxy2z_S_G4z_S;
  Double I_ERI_Gx3z_Py_G4z_S = I_ERI_Hxy3z_S_G4z_S+ABY*I_ERI_Gx3z_S_G4z_S;
  Double I_ERI_G4y_Py_G4z_S = I_ERI_H5y_S_G4z_S+ABY*I_ERI_G4y_S_G4z_S;
  Double I_ERI_G3yz_Py_G4z_S = I_ERI_H4yz_S_G4z_S+ABY*I_ERI_G3yz_S_G4z_S;
  Double I_ERI_G2y2z_Py_G4z_S = I_ERI_H3y2z_S_G4z_S+ABY*I_ERI_G2y2z_S_G4z_S;
  Double I_ERI_Gy3z_Py_G4z_S = I_ERI_H2y3z_S_G4z_S+ABY*I_ERI_Gy3z_S_G4z_S;
  Double I_ERI_G4z_Py_G4z_S = I_ERI_Hy4z_S_G4z_S+ABY*I_ERI_G4z_S_G4z_S;
  Double I_ERI_G4x_Pz_G4z_S = I_ERI_H4xz_S_G4z_S+ABZ*I_ERI_G4x_S_G4z_S;
  Double I_ERI_G3xy_Pz_G4z_S = I_ERI_H3xyz_S_G4z_S+ABZ*I_ERI_G3xy_S_G4z_S;
  Double I_ERI_G3xz_Pz_G4z_S = I_ERI_H3x2z_S_G4z_S+ABZ*I_ERI_G3xz_S_G4z_S;
  Double I_ERI_G2x2y_Pz_G4z_S = I_ERI_H2x2yz_S_G4z_S+ABZ*I_ERI_G2x2y_S_G4z_S;
  Double I_ERI_G2xyz_Pz_G4z_S = I_ERI_H2xy2z_S_G4z_S+ABZ*I_ERI_G2xyz_S_G4z_S;
  Double I_ERI_G2x2z_Pz_G4z_S = I_ERI_H2x3z_S_G4z_S+ABZ*I_ERI_G2x2z_S_G4z_S;
  Double I_ERI_Gx3y_Pz_G4z_S = I_ERI_Hx3yz_S_G4z_S+ABZ*I_ERI_Gx3y_S_G4z_S;
  Double I_ERI_Gx2yz_Pz_G4z_S = I_ERI_Hx2y2z_S_G4z_S+ABZ*I_ERI_Gx2yz_S_G4z_S;
  Double I_ERI_Gxy2z_Pz_G4z_S = I_ERI_Hxy3z_S_G4z_S+ABZ*I_ERI_Gxy2z_S_G4z_S;
  Double I_ERI_Gx3z_Pz_G4z_S = I_ERI_Hx4z_S_G4z_S+ABZ*I_ERI_Gx3z_S_G4z_S;
  Double I_ERI_G4y_Pz_G4z_S = I_ERI_H4yz_S_G4z_S+ABZ*I_ERI_G4y_S_G4z_S;
  Double I_ERI_G3yz_Pz_G4z_S = I_ERI_H3y2z_S_G4z_S+ABZ*I_ERI_G3yz_S_G4z_S;
  Double I_ERI_G2y2z_Pz_G4z_S = I_ERI_H2y3z_S_G4z_S+ABZ*I_ERI_G2y2z_S_G4z_S;
  Double I_ERI_Gy3z_Pz_G4z_S = I_ERI_Hy4z_S_G4z_S+ABZ*I_ERI_Gy3z_S_G4z_S;
  Double I_ERI_G4z_Pz_G4z_S = I_ERI_H5z_S_G4z_S+ABZ*I_ERI_G4z_S_G4z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_H_P_G_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 105 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_G_S
   * RHS shell quartet name: SQ_ERI_H_S_G_S
   ************************************************************/
  Double I_ERI_H5x_Px_G4x_S = I_ERI_I6x_S_G4x_S+ABX*I_ERI_H5x_S_G4x_S;
  Double I_ERI_H4xy_Px_G4x_S = I_ERI_I5xy_S_G4x_S+ABX*I_ERI_H4xy_S_G4x_S;
  Double I_ERI_H4xz_Px_G4x_S = I_ERI_I5xz_S_G4x_S+ABX*I_ERI_H4xz_S_G4x_S;
  Double I_ERI_H3x2y_Px_G4x_S = I_ERI_I4x2y_S_G4x_S+ABX*I_ERI_H3x2y_S_G4x_S;
  Double I_ERI_H3xyz_Px_G4x_S = I_ERI_I4xyz_S_G4x_S+ABX*I_ERI_H3xyz_S_G4x_S;
  Double I_ERI_H3x2z_Px_G4x_S = I_ERI_I4x2z_S_G4x_S+ABX*I_ERI_H3x2z_S_G4x_S;
  Double I_ERI_H2x3y_Px_G4x_S = I_ERI_I3x3y_S_G4x_S+ABX*I_ERI_H2x3y_S_G4x_S;
  Double I_ERI_H2x2yz_Px_G4x_S = I_ERI_I3x2yz_S_G4x_S+ABX*I_ERI_H2x2yz_S_G4x_S;
  Double I_ERI_H2xy2z_Px_G4x_S = I_ERI_I3xy2z_S_G4x_S+ABX*I_ERI_H2xy2z_S_G4x_S;
  Double I_ERI_H2x3z_Px_G4x_S = I_ERI_I3x3z_S_G4x_S+ABX*I_ERI_H2x3z_S_G4x_S;
  Double I_ERI_Hx4y_Px_G4x_S = I_ERI_I2x4y_S_G4x_S+ABX*I_ERI_Hx4y_S_G4x_S;
  Double I_ERI_Hx3yz_Px_G4x_S = I_ERI_I2x3yz_S_G4x_S+ABX*I_ERI_Hx3yz_S_G4x_S;
  Double I_ERI_Hx2y2z_Px_G4x_S = I_ERI_I2x2y2z_S_G4x_S+ABX*I_ERI_Hx2y2z_S_G4x_S;
  Double I_ERI_Hxy3z_Px_G4x_S = I_ERI_I2xy3z_S_G4x_S+ABX*I_ERI_Hxy3z_S_G4x_S;
  Double I_ERI_Hx4z_Px_G4x_S = I_ERI_I2x4z_S_G4x_S+ABX*I_ERI_Hx4z_S_G4x_S;
  Double I_ERI_H5y_Px_G4x_S = I_ERI_Ix5y_S_G4x_S+ABX*I_ERI_H5y_S_G4x_S;
  Double I_ERI_H4yz_Px_G4x_S = I_ERI_Ix4yz_S_G4x_S+ABX*I_ERI_H4yz_S_G4x_S;
  Double I_ERI_H3y2z_Px_G4x_S = I_ERI_Ix3y2z_S_G4x_S+ABX*I_ERI_H3y2z_S_G4x_S;
  Double I_ERI_H2y3z_Px_G4x_S = I_ERI_Ix2y3z_S_G4x_S+ABX*I_ERI_H2y3z_S_G4x_S;
  Double I_ERI_Hy4z_Px_G4x_S = I_ERI_Ixy4z_S_G4x_S+ABX*I_ERI_Hy4z_S_G4x_S;
  Double I_ERI_H5z_Px_G4x_S = I_ERI_Ix5z_S_G4x_S+ABX*I_ERI_H5z_S_G4x_S;
  Double I_ERI_H4xy_Py_G4x_S = I_ERI_I4x2y_S_G4x_S+ABY*I_ERI_H4xy_S_G4x_S;
  Double I_ERI_H4xz_Py_G4x_S = I_ERI_I4xyz_S_G4x_S+ABY*I_ERI_H4xz_S_G4x_S;
  Double I_ERI_H3x2y_Py_G4x_S = I_ERI_I3x3y_S_G4x_S+ABY*I_ERI_H3x2y_S_G4x_S;
  Double I_ERI_H3xyz_Py_G4x_S = I_ERI_I3x2yz_S_G4x_S+ABY*I_ERI_H3xyz_S_G4x_S;
  Double I_ERI_H3x2z_Py_G4x_S = I_ERI_I3xy2z_S_G4x_S+ABY*I_ERI_H3x2z_S_G4x_S;
  Double I_ERI_H2x3y_Py_G4x_S = I_ERI_I2x4y_S_G4x_S+ABY*I_ERI_H2x3y_S_G4x_S;
  Double I_ERI_H2x2yz_Py_G4x_S = I_ERI_I2x3yz_S_G4x_S+ABY*I_ERI_H2x2yz_S_G4x_S;
  Double I_ERI_H2xy2z_Py_G4x_S = I_ERI_I2x2y2z_S_G4x_S+ABY*I_ERI_H2xy2z_S_G4x_S;
  Double I_ERI_H2x3z_Py_G4x_S = I_ERI_I2xy3z_S_G4x_S+ABY*I_ERI_H2x3z_S_G4x_S;
  Double I_ERI_Hx4y_Py_G4x_S = I_ERI_Ix5y_S_G4x_S+ABY*I_ERI_Hx4y_S_G4x_S;
  Double I_ERI_Hx3yz_Py_G4x_S = I_ERI_Ix4yz_S_G4x_S+ABY*I_ERI_Hx3yz_S_G4x_S;
  Double I_ERI_Hx2y2z_Py_G4x_S = I_ERI_Ix3y2z_S_G4x_S+ABY*I_ERI_Hx2y2z_S_G4x_S;
  Double I_ERI_Hxy3z_Py_G4x_S = I_ERI_Ix2y3z_S_G4x_S+ABY*I_ERI_Hxy3z_S_G4x_S;
  Double I_ERI_Hx4z_Py_G4x_S = I_ERI_Ixy4z_S_G4x_S+ABY*I_ERI_Hx4z_S_G4x_S;
  Double I_ERI_H5y_Py_G4x_S = I_ERI_I6y_S_G4x_S+ABY*I_ERI_H5y_S_G4x_S;
  Double I_ERI_H4yz_Py_G4x_S = I_ERI_I5yz_S_G4x_S+ABY*I_ERI_H4yz_S_G4x_S;
  Double I_ERI_H3y2z_Py_G4x_S = I_ERI_I4y2z_S_G4x_S+ABY*I_ERI_H3y2z_S_G4x_S;
  Double I_ERI_H2y3z_Py_G4x_S = I_ERI_I3y3z_S_G4x_S+ABY*I_ERI_H2y3z_S_G4x_S;
  Double I_ERI_Hy4z_Py_G4x_S = I_ERI_I2y4z_S_G4x_S+ABY*I_ERI_Hy4z_S_G4x_S;
  Double I_ERI_H5z_Py_G4x_S = I_ERI_Iy5z_S_G4x_S+ABY*I_ERI_H5z_S_G4x_S;
  Double I_ERI_H4xz_Pz_G4x_S = I_ERI_I4x2z_S_G4x_S+ABZ*I_ERI_H4xz_S_G4x_S;
  Double I_ERI_H3xyz_Pz_G4x_S = I_ERI_I3xy2z_S_G4x_S+ABZ*I_ERI_H3xyz_S_G4x_S;
  Double I_ERI_H3x2z_Pz_G4x_S = I_ERI_I3x3z_S_G4x_S+ABZ*I_ERI_H3x2z_S_G4x_S;
  Double I_ERI_H2x2yz_Pz_G4x_S = I_ERI_I2x2y2z_S_G4x_S+ABZ*I_ERI_H2x2yz_S_G4x_S;
  Double I_ERI_H2xy2z_Pz_G4x_S = I_ERI_I2xy3z_S_G4x_S+ABZ*I_ERI_H2xy2z_S_G4x_S;
  Double I_ERI_H2x3z_Pz_G4x_S = I_ERI_I2x4z_S_G4x_S+ABZ*I_ERI_H2x3z_S_G4x_S;
  Double I_ERI_Hx3yz_Pz_G4x_S = I_ERI_Ix3y2z_S_G4x_S+ABZ*I_ERI_Hx3yz_S_G4x_S;
  Double I_ERI_Hx2y2z_Pz_G4x_S = I_ERI_Ix2y3z_S_G4x_S+ABZ*I_ERI_Hx2y2z_S_G4x_S;
  Double I_ERI_Hxy3z_Pz_G4x_S = I_ERI_Ixy4z_S_G4x_S+ABZ*I_ERI_Hxy3z_S_G4x_S;
  Double I_ERI_Hx4z_Pz_G4x_S = I_ERI_Ix5z_S_G4x_S+ABZ*I_ERI_Hx4z_S_G4x_S;
  Double I_ERI_H4yz_Pz_G4x_S = I_ERI_I4y2z_S_G4x_S+ABZ*I_ERI_H4yz_S_G4x_S;
  Double I_ERI_H3y2z_Pz_G4x_S = I_ERI_I3y3z_S_G4x_S+ABZ*I_ERI_H3y2z_S_G4x_S;
  Double I_ERI_H2y3z_Pz_G4x_S = I_ERI_I2y4z_S_G4x_S+ABZ*I_ERI_H2y3z_S_G4x_S;
  Double I_ERI_Hy4z_Pz_G4x_S = I_ERI_Iy5z_S_G4x_S+ABZ*I_ERI_Hy4z_S_G4x_S;
  Double I_ERI_H5z_Pz_G4x_S = I_ERI_I6z_S_G4x_S+ABZ*I_ERI_H5z_S_G4x_S;
  Double I_ERI_H5x_Px_G3xy_S = I_ERI_I6x_S_G3xy_S+ABX*I_ERI_H5x_S_G3xy_S;
  Double I_ERI_H4xy_Px_G3xy_S = I_ERI_I5xy_S_G3xy_S+ABX*I_ERI_H4xy_S_G3xy_S;
  Double I_ERI_H4xz_Px_G3xy_S = I_ERI_I5xz_S_G3xy_S+ABX*I_ERI_H4xz_S_G3xy_S;
  Double I_ERI_H3x2y_Px_G3xy_S = I_ERI_I4x2y_S_G3xy_S+ABX*I_ERI_H3x2y_S_G3xy_S;
  Double I_ERI_H3xyz_Px_G3xy_S = I_ERI_I4xyz_S_G3xy_S+ABX*I_ERI_H3xyz_S_G3xy_S;
  Double I_ERI_H3x2z_Px_G3xy_S = I_ERI_I4x2z_S_G3xy_S+ABX*I_ERI_H3x2z_S_G3xy_S;
  Double I_ERI_H2x3y_Px_G3xy_S = I_ERI_I3x3y_S_G3xy_S+ABX*I_ERI_H2x3y_S_G3xy_S;
  Double I_ERI_H2x2yz_Px_G3xy_S = I_ERI_I3x2yz_S_G3xy_S+ABX*I_ERI_H2x2yz_S_G3xy_S;
  Double I_ERI_H2xy2z_Px_G3xy_S = I_ERI_I3xy2z_S_G3xy_S+ABX*I_ERI_H2xy2z_S_G3xy_S;
  Double I_ERI_H2x3z_Px_G3xy_S = I_ERI_I3x3z_S_G3xy_S+ABX*I_ERI_H2x3z_S_G3xy_S;
  Double I_ERI_Hx4y_Px_G3xy_S = I_ERI_I2x4y_S_G3xy_S+ABX*I_ERI_Hx4y_S_G3xy_S;
  Double I_ERI_Hx3yz_Px_G3xy_S = I_ERI_I2x3yz_S_G3xy_S+ABX*I_ERI_Hx3yz_S_G3xy_S;
  Double I_ERI_Hx2y2z_Px_G3xy_S = I_ERI_I2x2y2z_S_G3xy_S+ABX*I_ERI_Hx2y2z_S_G3xy_S;
  Double I_ERI_Hxy3z_Px_G3xy_S = I_ERI_I2xy3z_S_G3xy_S+ABX*I_ERI_Hxy3z_S_G3xy_S;
  Double I_ERI_Hx4z_Px_G3xy_S = I_ERI_I2x4z_S_G3xy_S+ABX*I_ERI_Hx4z_S_G3xy_S;
  Double I_ERI_H5y_Px_G3xy_S = I_ERI_Ix5y_S_G3xy_S+ABX*I_ERI_H5y_S_G3xy_S;
  Double I_ERI_H4yz_Px_G3xy_S = I_ERI_Ix4yz_S_G3xy_S+ABX*I_ERI_H4yz_S_G3xy_S;
  Double I_ERI_H3y2z_Px_G3xy_S = I_ERI_Ix3y2z_S_G3xy_S+ABX*I_ERI_H3y2z_S_G3xy_S;
  Double I_ERI_H2y3z_Px_G3xy_S = I_ERI_Ix2y3z_S_G3xy_S+ABX*I_ERI_H2y3z_S_G3xy_S;
  Double I_ERI_Hy4z_Px_G3xy_S = I_ERI_Ixy4z_S_G3xy_S+ABX*I_ERI_Hy4z_S_G3xy_S;
  Double I_ERI_H5z_Px_G3xy_S = I_ERI_Ix5z_S_G3xy_S+ABX*I_ERI_H5z_S_G3xy_S;
  Double I_ERI_H4xy_Py_G3xy_S = I_ERI_I4x2y_S_G3xy_S+ABY*I_ERI_H4xy_S_G3xy_S;
  Double I_ERI_H4xz_Py_G3xy_S = I_ERI_I4xyz_S_G3xy_S+ABY*I_ERI_H4xz_S_G3xy_S;
  Double I_ERI_H3x2y_Py_G3xy_S = I_ERI_I3x3y_S_G3xy_S+ABY*I_ERI_H3x2y_S_G3xy_S;
  Double I_ERI_H3xyz_Py_G3xy_S = I_ERI_I3x2yz_S_G3xy_S+ABY*I_ERI_H3xyz_S_G3xy_S;
  Double I_ERI_H3x2z_Py_G3xy_S = I_ERI_I3xy2z_S_G3xy_S+ABY*I_ERI_H3x2z_S_G3xy_S;
  Double I_ERI_H2x3y_Py_G3xy_S = I_ERI_I2x4y_S_G3xy_S+ABY*I_ERI_H2x3y_S_G3xy_S;
  Double I_ERI_H2x2yz_Py_G3xy_S = I_ERI_I2x3yz_S_G3xy_S+ABY*I_ERI_H2x2yz_S_G3xy_S;
  Double I_ERI_H2xy2z_Py_G3xy_S = I_ERI_I2x2y2z_S_G3xy_S+ABY*I_ERI_H2xy2z_S_G3xy_S;
  Double I_ERI_H2x3z_Py_G3xy_S = I_ERI_I2xy3z_S_G3xy_S+ABY*I_ERI_H2x3z_S_G3xy_S;
  Double I_ERI_Hx4y_Py_G3xy_S = I_ERI_Ix5y_S_G3xy_S+ABY*I_ERI_Hx4y_S_G3xy_S;
  Double I_ERI_Hx3yz_Py_G3xy_S = I_ERI_Ix4yz_S_G3xy_S+ABY*I_ERI_Hx3yz_S_G3xy_S;
  Double I_ERI_Hx2y2z_Py_G3xy_S = I_ERI_Ix3y2z_S_G3xy_S+ABY*I_ERI_Hx2y2z_S_G3xy_S;
  Double I_ERI_Hxy3z_Py_G3xy_S = I_ERI_Ix2y3z_S_G3xy_S+ABY*I_ERI_Hxy3z_S_G3xy_S;
  Double I_ERI_Hx4z_Py_G3xy_S = I_ERI_Ixy4z_S_G3xy_S+ABY*I_ERI_Hx4z_S_G3xy_S;
  Double I_ERI_H5y_Py_G3xy_S = I_ERI_I6y_S_G3xy_S+ABY*I_ERI_H5y_S_G3xy_S;
  Double I_ERI_H4yz_Py_G3xy_S = I_ERI_I5yz_S_G3xy_S+ABY*I_ERI_H4yz_S_G3xy_S;
  Double I_ERI_H3y2z_Py_G3xy_S = I_ERI_I4y2z_S_G3xy_S+ABY*I_ERI_H3y2z_S_G3xy_S;
  Double I_ERI_H2y3z_Py_G3xy_S = I_ERI_I3y3z_S_G3xy_S+ABY*I_ERI_H2y3z_S_G3xy_S;
  Double I_ERI_Hy4z_Py_G3xy_S = I_ERI_I2y4z_S_G3xy_S+ABY*I_ERI_Hy4z_S_G3xy_S;
  Double I_ERI_H5z_Py_G3xy_S = I_ERI_Iy5z_S_G3xy_S+ABY*I_ERI_H5z_S_G3xy_S;
  Double I_ERI_H4xz_Pz_G3xy_S = I_ERI_I4x2z_S_G3xy_S+ABZ*I_ERI_H4xz_S_G3xy_S;
  Double I_ERI_H3xyz_Pz_G3xy_S = I_ERI_I3xy2z_S_G3xy_S+ABZ*I_ERI_H3xyz_S_G3xy_S;
  Double I_ERI_H3x2z_Pz_G3xy_S = I_ERI_I3x3z_S_G3xy_S+ABZ*I_ERI_H3x2z_S_G3xy_S;
  Double I_ERI_H2x2yz_Pz_G3xy_S = I_ERI_I2x2y2z_S_G3xy_S+ABZ*I_ERI_H2x2yz_S_G3xy_S;
  Double I_ERI_H2xy2z_Pz_G3xy_S = I_ERI_I2xy3z_S_G3xy_S+ABZ*I_ERI_H2xy2z_S_G3xy_S;
  Double I_ERI_H2x3z_Pz_G3xy_S = I_ERI_I2x4z_S_G3xy_S+ABZ*I_ERI_H2x3z_S_G3xy_S;
  Double I_ERI_Hx3yz_Pz_G3xy_S = I_ERI_Ix3y2z_S_G3xy_S+ABZ*I_ERI_Hx3yz_S_G3xy_S;
  Double I_ERI_Hx2y2z_Pz_G3xy_S = I_ERI_Ix2y3z_S_G3xy_S+ABZ*I_ERI_Hx2y2z_S_G3xy_S;
  Double I_ERI_Hxy3z_Pz_G3xy_S = I_ERI_Ixy4z_S_G3xy_S+ABZ*I_ERI_Hxy3z_S_G3xy_S;
  Double I_ERI_Hx4z_Pz_G3xy_S = I_ERI_Ix5z_S_G3xy_S+ABZ*I_ERI_Hx4z_S_G3xy_S;
  Double I_ERI_H4yz_Pz_G3xy_S = I_ERI_I4y2z_S_G3xy_S+ABZ*I_ERI_H4yz_S_G3xy_S;
  Double I_ERI_H3y2z_Pz_G3xy_S = I_ERI_I3y3z_S_G3xy_S+ABZ*I_ERI_H3y2z_S_G3xy_S;
  Double I_ERI_H2y3z_Pz_G3xy_S = I_ERI_I2y4z_S_G3xy_S+ABZ*I_ERI_H2y3z_S_G3xy_S;
  Double I_ERI_Hy4z_Pz_G3xy_S = I_ERI_Iy5z_S_G3xy_S+ABZ*I_ERI_Hy4z_S_G3xy_S;
  Double I_ERI_H5z_Pz_G3xy_S = I_ERI_I6z_S_G3xy_S+ABZ*I_ERI_H5z_S_G3xy_S;
  Double I_ERI_H5x_Px_G3xz_S = I_ERI_I6x_S_G3xz_S+ABX*I_ERI_H5x_S_G3xz_S;
  Double I_ERI_H4xy_Px_G3xz_S = I_ERI_I5xy_S_G3xz_S+ABX*I_ERI_H4xy_S_G3xz_S;
  Double I_ERI_H4xz_Px_G3xz_S = I_ERI_I5xz_S_G3xz_S+ABX*I_ERI_H4xz_S_G3xz_S;
  Double I_ERI_H3x2y_Px_G3xz_S = I_ERI_I4x2y_S_G3xz_S+ABX*I_ERI_H3x2y_S_G3xz_S;
  Double I_ERI_H3xyz_Px_G3xz_S = I_ERI_I4xyz_S_G3xz_S+ABX*I_ERI_H3xyz_S_G3xz_S;
  Double I_ERI_H3x2z_Px_G3xz_S = I_ERI_I4x2z_S_G3xz_S+ABX*I_ERI_H3x2z_S_G3xz_S;
  Double I_ERI_H2x3y_Px_G3xz_S = I_ERI_I3x3y_S_G3xz_S+ABX*I_ERI_H2x3y_S_G3xz_S;
  Double I_ERI_H2x2yz_Px_G3xz_S = I_ERI_I3x2yz_S_G3xz_S+ABX*I_ERI_H2x2yz_S_G3xz_S;
  Double I_ERI_H2xy2z_Px_G3xz_S = I_ERI_I3xy2z_S_G3xz_S+ABX*I_ERI_H2xy2z_S_G3xz_S;
  Double I_ERI_H2x3z_Px_G3xz_S = I_ERI_I3x3z_S_G3xz_S+ABX*I_ERI_H2x3z_S_G3xz_S;
  Double I_ERI_Hx4y_Px_G3xz_S = I_ERI_I2x4y_S_G3xz_S+ABX*I_ERI_Hx4y_S_G3xz_S;
  Double I_ERI_Hx3yz_Px_G3xz_S = I_ERI_I2x3yz_S_G3xz_S+ABX*I_ERI_Hx3yz_S_G3xz_S;
  Double I_ERI_Hx2y2z_Px_G3xz_S = I_ERI_I2x2y2z_S_G3xz_S+ABX*I_ERI_Hx2y2z_S_G3xz_S;
  Double I_ERI_Hxy3z_Px_G3xz_S = I_ERI_I2xy3z_S_G3xz_S+ABX*I_ERI_Hxy3z_S_G3xz_S;
  Double I_ERI_Hx4z_Px_G3xz_S = I_ERI_I2x4z_S_G3xz_S+ABX*I_ERI_Hx4z_S_G3xz_S;
  Double I_ERI_H5y_Px_G3xz_S = I_ERI_Ix5y_S_G3xz_S+ABX*I_ERI_H5y_S_G3xz_S;
  Double I_ERI_H4yz_Px_G3xz_S = I_ERI_Ix4yz_S_G3xz_S+ABX*I_ERI_H4yz_S_G3xz_S;
  Double I_ERI_H3y2z_Px_G3xz_S = I_ERI_Ix3y2z_S_G3xz_S+ABX*I_ERI_H3y2z_S_G3xz_S;
  Double I_ERI_H2y3z_Px_G3xz_S = I_ERI_Ix2y3z_S_G3xz_S+ABX*I_ERI_H2y3z_S_G3xz_S;
  Double I_ERI_Hy4z_Px_G3xz_S = I_ERI_Ixy4z_S_G3xz_S+ABX*I_ERI_Hy4z_S_G3xz_S;
  Double I_ERI_H5z_Px_G3xz_S = I_ERI_Ix5z_S_G3xz_S+ABX*I_ERI_H5z_S_G3xz_S;
  Double I_ERI_H4xy_Py_G3xz_S = I_ERI_I4x2y_S_G3xz_S+ABY*I_ERI_H4xy_S_G3xz_S;
  Double I_ERI_H4xz_Py_G3xz_S = I_ERI_I4xyz_S_G3xz_S+ABY*I_ERI_H4xz_S_G3xz_S;
  Double I_ERI_H3x2y_Py_G3xz_S = I_ERI_I3x3y_S_G3xz_S+ABY*I_ERI_H3x2y_S_G3xz_S;
  Double I_ERI_H3xyz_Py_G3xz_S = I_ERI_I3x2yz_S_G3xz_S+ABY*I_ERI_H3xyz_S_G3xz_S;
  Double I_ERI_H3x2z_Py_G3xz_S = I_ERI_I3xy2z_S_G3xz_S+ABY*I_ERI_H3x2z_S_G3xz_S;
  Double I_ERI_H2x3y_Py_G3xz_S = I_ERI_I2x4y_S_G3xz_S+ABY*I_ERI_H2x3y_S_G3xz_S;
  Double I_ERI_H2x2yz_Py_G3xz_S = I_ERI_I2x3yz_S_G3xz_S+ABY*I_ERI_H2x2yz_S_G3xz_S;
  Double I_ERI_H2xy2z_Py_G3xz_S = I_ERI_I2x2y2z_S_G3xz_S+ABY*I_ERI_H2xy2z_S_G3xz_S;
  Double I_ERI_H2x3z_Py_G3xz_S = I_ERI_I2xy3z_S_G3xz_S+ABY*I_ERI_H2x3z_S_G3xz_S;
  Double I_ERI_Hx4y_Py_G3xz_S = I_ERI_Ix5y_S_G3xz_S+ABY*I_ERI_Hx4y_S_G3xz_S;
  Double I_ERI_Hx3yz_Py_G3xz_S = I_ERI_Ix4yz_S_G3xz_S+ABY*I_ERI_Hx3yz_S_G3xz_S;
  Double I_ERI_Hx2y2z_Py_G3xz_S = I_ERI_Ix3y2z_S_G3xz_S+ABY*I_ERI_Hx2y2z_S_G3xz_S;
  Double I_ERI_Hxy3z_Py_G3xz_S = I_ERI_Ix2y3z_S_G3xz_S+ABY*I_ERI_Hxy3z_S_G3xz_S;
  Double I_ERI_Hx4z_Py_G3xz_S = I_ERI_Ixy4z_S_G3xz_S+ABY*I_ERI_Hx4z_S_G3xz_S;
  Double I_ERI_H5y_Py_G3xz_S = I_ERI_I6y_S_G3xz_S+ABY*I_ERI_H5y_S_G3xz_S;
  Double I_ERI_H4yz_Py_G3xz_S = I_ERI_I5yz_S_G3xz_S+ABY*I_ERI_H4yz_S_G3xz_S;
  Double I_ERI_H3y2z_Py_G3xz_S = I_ERI_I4y2z_S_G3xz_S+ABY*I_ERI_H3y2z_S_G3xz_S;
  Double I_ERI_H2y3z_Py_G3xz_S = I_ERI_I3y3z_S_G3xz_S+ABY*I_ERI_H2y3z_S_G3xz_S;
  Double I_ERI_Hy4z_Py_G3xz_S = I_ERI_I2y4z_S_G3xz_S+ABY*I_ERI_Hy4z_S_G3xz_S;
  Double I_ERI_H5z_Py_G3xz_S = I_ERI_Iy5z_S_G3xz_S+ABY*I_ERI_H5z_S_G3xz_S;
  Double I_ERI_H4xz_Pz_G3xz_S = I_ERI_I4x2z_S_G3xz_S+ABZ*I_ERI_H4xz_S_G3xz_S;
  Double I_ERI_H3xyz_Pz_G3xz_S = I_ERI_I3xy2z_S_G3xz_S+ABZ*I_ERI_H3xyz_S_G3xz_S;
  Double I_ERI_H3x2z_Pz_G3xz_S = I_ERI_I3x3z_S_G3xz_S+ABZ*I_ERI_H3x2z_S_G3xz_S;
  Double I_ERI_H2x2yz_Pz_G3xz_S = I_ERI_I2x2y2z_S_G3xz_S+ABZ*I_ERI_H2x2yz_S_G3xz_S;
  Double I_ERI_H2xy2z_Pz_G3xz_S = I_ERI_I2xy3z_S_G3xz_S+ABZ*I_ERI_H2xy2z_S_G3xz_S;
  Double I_ERI_H2x3z_Pz_G3xz_S = I_ERI_I2x4z_S_G3xz_S+ABZ*I_ERI_H2x3z_S_G3xz_S;
  Double I_ERI_Hx3yz_Pz_G3xz_S = I_ERI_Ix3y2z_S_G3xz_S+ABZ*I_ERI_Hx3yz_S_G3xz_S;
  Double I_ERI_Hx2y2z_Pz_G3xz_S = I_ERI_Ix2y3z_S_G3xz_S+ABZ*I_ERI_Hx2y2z_S_G3xz_S;
  Double I_ERI_Hxy3z_Pz_G3xz_S = I_ERI_Ixy4z_S_G3xz_S+ABZ*I_ERI_Hxy3z_S_G3xz_S;
  Double I_ERI_Hx4z_Pz_G3xz_S = I_ERI_Ix5z_S_G3xz_S+ABZ*I_ERI_Hx4z_S_G3xz_S;
  Double I_ERI_H4yz_Pz_G3xz_S = I_ERI_I4y2z_S_G3xz_S+ABZ*I_ERI_H4yz_S_G3xz_S;
  Double I_ERI_H3y2z_Pz_G3xz_S = I_ERI_I3y3z_S_G3xz_S+ABZ*I_ERI_H3y2z_S_G3xz_S;
  Double I_ERI_H2y3z_Pz_G3xz_S = I_ERI_I2y4z_S_G3xz_S+ABZ*I_ERI_H2y3z_S_G3xz_S;
  Double I_ERI_Hy4z_Pz_G3xz_S = I_ERI_Iy5z_S_G3xz_S+ABZ*I_ERI_Hy4z_S_G3xz_S;
  Double I_ERI_H5z_Pz_G3xz_S = I_ERI_I6z_S_G3xz_S+ABZ*I_ERI_H5z_S_G3xz_S;
  Double I_ERI_H5x_Px_G2x2y_S = I_ERI_I6x_S_G2x2y_S+ABX*I_ERI_H5x_S_G2x2y_S;
  Double I_ERI_H4xy_Px_G2x2y_S = I_ERI_I5xy_S_G2x2y_S+ABX*I_ERI_H4xy_S_G2x2y_S;
  Double I_ERI_H4xz_Px_G2x2y_S = I_ERI_I5xz_S_G2x2y_S+ABX*I_ERI_H4xz_S_G2x2y_S;
  Double I_ERI_H3x2y_Px_G2x2y_S = I_ERI_I4x2y_S_G2x2y_S+ABX*I_ERI_H3x2y_S_G2x2y_S;
  Double I_ERI_H3xyz_Px_G2x2y_S = I_ERI_I4xyz_S_G2x2y_S+ABX*I_ERI_H3xyz_S_G2x2y_S;
  Double I_ERI_H3x2z_Px_G2x2y_S = I_ERI_I4x2z_S_G2x2y_S+ABX*I_ERI_H3x2z_S_G2x2y_S;
  Double I_ERI_H2x3y_Px_G2x2y_S = I_ERI_I3x3y_S_G2x2y_S+ABX*I_ERI_H2x3y_S_G2x2y_S;
  Double I_ERI_H2x2yz_Px_G2x2y_S = I_ERI_I3x2yz_S_G2x2y_S+ABX*I_ERI_H2x2yz_S_G2x2y_S;
  Double I_ERI_H2xy2z_Px_G2x2y_S = I_ERI_I3xy2z_S_G2x2y_S+ABX*I_ERI_H2xy2z_S_G2x2y_S;
  Double I_ERI_H2x3z_Px_G2x2y_S = I_ERI_I3x3z_S_G2x2y_S+ABX*I_ERI_H2x3z_S_G2x2y_S;
  Double I_ERI_Hx4y_Px_G2x2y_S = I_ERI_I2x4y_S_G2x2y_S+ABX*I_ERI_Hx4y_S_G2x2y_S;
  Double I_ERI_Hx3yz_Px_G2x2y_S = I_ERI_I2x3yz_S_G2x2y_S+ABX*I_ERI_Hx3yz_S_G2x2y_S;
  Double I_ERI_Hx2y2z_Px_G2x2y_S = I_ERI_I2x2y2z_S_G2x2y_S+ABX*I_ERI_Hx2y2z_S_G2x2y_S;
  Double I_ERI_Hxy3z_Px_G2x2y_S = I_ERI_I2xy3z_S_G2x2y_S+ABX*I_ERI_Hxy3z_S_G2x2y_S;
  Double I_ERI_Hx4z_Px_G2x2y_S = I_ERI_I2x4z_S_G2x2y_S+ABX*I_ERI_Hx4z_S_G2x2y_S;
  Double I_ERI_H5y_Px_G2x2y_S = I_ERI_Ix5y_S_G2x2y_S+ABX*I_ERI_H5y_S_G2x2y_S;
  Double I_ERI_H4yz_Px_G2x2y_S = I_ERI_Ix4yz_S_G2x2y_S+ABX*I_ERI_H4yz_S_G2x2y_S;
  Double I_ERI_H3y2z_Px_G2x2y_S = I_ERI_Ix3y2z_S_G2x2y_S+ABX*I_ERI_H3y2z_S_G2x2y_S;
  Double I_ERI_H2y3z_Px_G2x2y_S = I_ERI_Ix2y3z_S_G2x2y_S+ABX*I_ERI_H2y3z_S_G2x2y_S;
  Double I_ERI_Hy4z_Px_G2x2y_S = I_ERI_Ixy4z_S_G2x2y_S+ABX*I_ERI_Hy4z_S_G2x2y_S;
  Double I_ERI_H5z_Px_G2x2y_S = I_ERI_Ix5z_S_G2x2y_S+ABX*I_ERI_H5z_S_G2x2y_S;
  Double I_ERI_H4xy_Py_G2x2y_S = I_ERI_I4x2y_S_G2x2y_S+ABY*I_ERI_H4xy_S_G2x2y_S;
  Double I_ERI_H4xz_Py_G2x2y_S = I_ERI_I4xyz_S_G2x2y_S+ABY*I_ERI_H4xz_S_G2x2y_S;
  Double I_ERI_H3x2y_Py_G2x2y_S = I_ERI_I3x3y_S_G2x2y_S+ABY*I_ERI_H3x2y_S_G2x2y_S;
  Double I_ERI_H3xyz_Py_G2x2y_S = I_ERI_I3x2yz_S_G2x2y_S+ABY*I_ERI_H3xyz_S_G2x2y_S;
  Double I_ERI_H3x2z_Py_G2x2y_S = I_ERI_I3xy2z_S_G2x2y_S+ABY*I_ERI_H3x2z_S_G2x2y_S;
  Double I_ERI_H2x3y_Py_G2x2y_S = I_ERI_I2x4y_S_G2x2y_S+ABY*I_ERI_H2x3y_S_G2x2y_S;
  Double I_ERI_H2x2yz_Py_G2x2y_S = I_ERI_I2x3yz_S_G2x2y_S+ABY*I_ERI_H2x2yz_S_G2x2y_S;
  Double I_ERI_H2xy2z_Py_G2x2y_S = I_ERI_I2x2y2z_S_G2x2y_S+ABY*I_ERI_H2xy2z_S_G2x2y_S;
  Double I_ERI_H2x3z_Py_G2x2y_S = I_ERI_I2xy3z_S_G2x2y_S+ABY*I_ERI_H2x3z_S_G2x2y_S;
  Double I_ERI_Hx4y_Py_G2x2y_S = I_ERI_Ix5y_S_G2x2y_S+ABY*I_ERI_Hx4y_S_G2x2y_S;
  Double I_ERI_Hx3yz_Py_G2x2y_S = I_ERI_Ix4yz_S_G2x2y_S+ABY*I_ERI_Hx3yz_S_G2x2y_S;
  Double I_ERI_Hx2y2z_Py_G2x2y_S = I_ERI_Ix3y2z_S_G2x2y_S+ABY*I_ERI_Hx2y2z_S_G2x2y_S;
  Double I_ERI_Hxy3z_Py_G2x2y_S = I_ERI_Ix2y3z_S_G2x2y_S+ABY*I_ERI_Hxy3z_S_G2x2y_S;
  Double I_ERI_Hx4z_Py_G2x2y_S = I_ERI_Ixy4z_S_G2x2y_S+ABY*I_ERI_Hx4z_S_G2x2y_S;
  Double I_ERI_H5y_Py_G2x2y_S = I_ERI_I6y_S_G2x2y_S+ABY*I_ERI_H5y_S_G2x2y_S;
  Double I_ERI_H4yz_Py_G2x2y_S = I_ERI_I5yz_S_G2x2y_S+ABY*I_ERI_H4yz_S_G2x2y_S;
  Double I_ERI_H3y2z_Py_G2x2y_S = I_ERI_I4y2z_S_G2x2y_S+ABY*I_ERI_H3y2z_S_G2x2y_S;
  Double I_ERI_H2y3z_Py_G2x2y_S = I_ERI_I3y3z_S_G2x2y_S+ABY*I_ERI_H2y3z_S_G2x2y_S;
  Double I_ERI_Hy4z_Py_G2x2y_S = I_ERI_I2y4z_S_G2x2y_S+ABY*I_ERI_Hy4z_S_G2x2y_S;
  Double I_ERI_H5z_Py_G2x2y_S = I_ERI_Iy5z_S_G2x2y_S+ABY*I_ERI_H5z_S_G2x2y_S;
  Double I_ERI_H4xz_Pz_G2x2y_S = I_ERI_I4x2z_S_G2x2y_S+ABZ*I_ERI_H4xz_S_G2x2y_S;
  Double I_ERI_H3xyz_Pz_G2x2y_S = I_ERI_I3xy2z_S_G2x2y_S+ABZ*I_ERI_H3xyz_S_G2x2y_S;
  Double I_ERI_H3x2z_Pz_G2x2y_S = I_ERI_I3x3z_S_G2x2y_S+ABZ*I_ERI_H3x2z_S_G2x2y_S;
  Double I_ERI_H2x2yz_Pz_G2x2y_S = I_ERI_I2x2y2z_S_G2x2y_S+ABZ*I_ERI_H2x2yz_S_G2x2y_S;
  Double I_ERI_H2xy2z_Pz_G2x2y_S = I_ERI_I2xy3z_S_G2x2y_S+ABZ*I_ERI_H2xy2z_S_G2x2y_S;
  Double I_ERI_H2x3z_Pz_G2x2y_S = I_ERI_I2x4z_S_G2x2y_S+ABZ*I_ERI_H2x3z_S_G2x2y_S;
  Double I_ERI_Hx3yz_Pz_G2x2y_S = I_ERI_Ix3y2z_S_G2x2y_S+ABZ*I_ERI_Hx3yz_S_G2x2y_S;
  Double I_ERI_Hx2y2z_Pz_G2x2y_S = I_ERI_Ix2y3z_S_G2x2y_S+ABZ*I_ERI_Hx2y2z_S_G2x2y_S;
  Double I_ERI_Hxy3z_Pz_G2x2y_S = I_ERI_Ixy4z_S_G2x2y_S+ABZ*I_ERI_Hxy3z_S_G2x2y_S;
  Double I_ERI_Hx4z_Pz_G2x2y_S = I_ERI_Ix5z_S_G2x2y_S+ABZ*I_ERI_Hx4z_S_G2x2y_S;
  Double I_ERI_H4yz_Pz_G2x2y_S = I_ERI_I4y2z_S_G2x2y_S+ABZ*I_ERI_H4yz_S_G2x2y_S;
  Double I_ERI_H3y2z_Pz_G2x2y_S = I_ERI_I3y3z_S_G2x2y_S+ABZ*I_ERI_H3y2z_S_G2x2y_S;
  Double I_ERI_H2y3z_Pz_G2x2y_S = I_ERI_I2y4z_S_G2x2y_S+ABZ*I_ERI_H2y3z_S_G2x2y_S;
  Double I_ERI_Hy4z_Pz_G2x2y_S = I_ERI_Iy5z_S_G2x2y_S+ABZ*I_ERI_Hy4z_S_G2x2y_S;
  Double I_ERI_H5z_Pz_G2x2y_S = I_ERI_I6z_S_G2x2y_S+ABZ*I_ERI_H5z_S_G2x2y_S;
  Double I_ERI_H5x_Px_G2xyz_S = I_ERI_I6x_S_G2xyz_S+ABX*I_ERI_H5x_S_G2xyz_S;
  Double I_ERI_H4xy_Px_G2xyz_S = I_ERI_I5xy_S_G2xyz_S+ABX*I_ERI_H4xy_S_G2xyz_S;
  Double I_ERI_H4xz_Px_G2xyz_S = I_ERI_I5xz_S_G2xyz_S+ABX*I_ERI_H4xz_S_G2xyz_S;
  Double I_ERI_H3x2y_Px_G2xyz_S = I_ERI_I4x2y_S_G2xyz_S+ABX*I_ERI_H3x2y_S_G2xyz_S;
  Double I_ERI_H3xyz_Px_G2xyz_S = I_ERI_I4xyz_S_G2xyz_S+ABX*I_ERI_H3xyz_S_G2xyz_S;
  Double I_ERI_H3x2z_Px_G2xyz_S = I_ERI_I4x2z_S_G2xyz_S+ABX*I_ERI_H3x2z_S_G2xyz_S;
  Double I_ERI_H2x3y_Px_G2xyz_S = I_ERI_I3x3y_S_G2xyz_S+ABX*I_ERI_H2x3y_S_G2xyz_S;
  Double I_ERI_H2x2yz_Px_G2xyz_S = I_ERI_I3x2yz_S_G2xyz_S+ABX*I_ERI_H2x2yz_S_G2xyz_S;
  Double I_ERI_H2xy2z_Px_G2xyz_S = I_ERI_I3xy2z_S_G2xyz_S+ABX*I_ERI_H2xy2z_S_G2xyz_S;
  Double I_ERI_H2x3z_Px_G2xyz_S = I_ERI_I3x3z_S_G2xyz_S+ABX*I_ERI_H2x3z_S_G2xyz_S;
  Double I_ERI_Hx4y_Px_G2xyz_S = I_ERI_I2x4y_S_G2xyz_S+ABX*I_ERI_Hx4y_S_G2xyz_S;
  Double I_ERI_Hx3yz_Px_G2xyz_S = I_ERI_I2x3yz_S_G2xyz_S+ABX*I_ERI_Hx3yz_S_G2xyz_S;
  Double I_ERI_Hx2y2z_Px_G2xyz_S = I_ERI_I2x2y2z_S_G2xyz_S+ABX*I_ERI_Hx2y2z_S_G2xyz_S;
  Double I_ERI_Hxy3z_Px_G2xyz_S = I_ERI_I2xy3z_S_G2xyz_S+ABX*I_ERI_Hxy3z_S_G2xyz_S;
  Double I_ERI_Hx4z_Px_G2xyz_S = I_ERI_I2x4z_S_G2xyz_S+ABX*I_ERI_Hx4z_S_G2xyz_S;
  Double I_ERI_H5y_Px_G2xyz_S = I_ERI_Ix5y_S_G2xyz_S+ABX*I_ERI_H5y_S_G2xyz_S;
  Double I_ERI_H4yz_Px_G2xyz_S = I_ERI_Ix4yz_S_G2xyz_S+ABX*I_ERI_H4yz_S_G2xyz_S;
  Double I_ERI_H3y2z_Px_G2xyz_S = I_ERI_Ix3y2z_S_G2xyz_S+ABX*I_ERI_H3y2z_S_G2xyz_S;
  Double I_ERI_H2y3z_Px_G2xyz_S = I_ERI_Ix2y3z_S_G2xyz_S+ABX*I_ERI_H2y3z_S_G2xyz_S;
  Double I_ERI_Hy4z_Px_G2xyz_S = I_ERI_Ixy4z_S_G2xyz_S+ABX*I_ERI_Hy4z_S_G2xyz_S;
  Double I_ERI_H5z_Px_G2xyz_S = I_ERI_Ix5z_S_G2xyz_S+ABX*I_ERI_H5z_S_G2xyz_S;
  Double I_ERI_H4xy_Py_G2xyz_S = I_ERI_I4x2y_S_G2xyz_S+ABY*I_ERI_H4xy_S_G2xyz_S;
  Double I_ERI_H4xz_Py_G2xyz_S = I_ERI_I4xyz_S_G2xyz_S+ABY*I_ERI_H4xz_S_G2xyz_S;
  Double I_ERI_H3x2y_Py_G2xyz_S = I_ERI_I3x3y_S_G2xyz_S+ABY*I_ERI_H3x2y_S_G2xyz_S;
  Double I_ERI_H3xyz_Py_G2xyz_S = I_ERI_I3x2yz_S_G2xyz_S+ABY*I_ERI_H3xyz_S_G2xyz_S;
  Double I_ERI_H3x2z_Py_G2xyz_S = I_ERI_I3xy2z_S_G2xyz_S+ABY*I_ERI_H3x2z_S_G2xyz_S;
  Double I_ERI_H2x3y_Py_G2xyz_S = I_ERI_I2x4y_S_G2xyz_S+ABY*I_ERI_H2x3y_S_G2xyz_S;
  Double I_ERI_H2x2yz_Py_G2xyz_S = I_ERI_I2x3yz_S_G2xyz_S+ABY*I_ERI_H2x2yz_S_G2xyz_S;
  Double I_ERI_H2xy2z_Py_G2xyz_S = I_ERI_I2x2y2z_S_G2xyz_S+ABY*I_ERI_H2xy2z_S_G2xyz_S;
  Double I_ERI_H2x3z_Py_G2xyz_S = I_ERI_I2xy3z_S_G2xyz_S+ABY*I_ERI_H2x3z_S_G2xyz_S;
  Double I_ERI_Hx4y_Py_G2xyz_S = I_ERI_Ix5y_S_G2xyz_S+ABY*I_ERI_Hx4y_S_G2xyz_S;
  Double I_ERI_Hx3yz_Py_G2xyz_S = I_ERI_Ix4yz_S_G2xyz_S+ABY*I_ERI_Hx3yz_S_G2xyz_S;
  Double I_ERI_Hx2y2z_Py_G2xyz_S = I_ERI_Ix3y2z_S_G2xyz_S+ABY*I_ERI_Hx2y2z_S_G2xyz_S;
  Double I_ERI_Hxy3z_Py_G2xyz_S = I_ERI_Ix2y3z_S_G2xyz_S+ABY*I_ERI_Hxy3z_S_G2xyz_S;
  Double I_ERI_Hx4z_Py_G2xyz_S = I_ERI_Ixy4z_S_G2xyz_S+ABY*I_ERI_Hx4z_S_G2xyz_S;
  Double I_ERI_H5y_Py_G2xyz_S = I_ERI_I6y_S_G2xyz_S+ABY*I_ERI_H5y_S_G2xyz_S;
  Double I_ERI_H4yz_Py_G2xyz_S = I_ERI_I5yz_S_G2xyz_S+ABY*I_ERI_H4yz_S_G2xyz_S;
  Double I_ERI_H3y2z_Py_G2xyz_S = I_ERI_I4y2z_S_G2xyz_S+ABY*I_ERI_H3y2z_S_G2xyz_S;
  Double I_ERI_H2y3z_Py_G2xyz_S = I_ERI_I3y3z_S_G2xyz_S+ABY*I_ERI_H2y3z_S_G2xyz_S;
  Double I_ERI_Hy4z_Py_G2xyz_S = I_ERI_I2y4z_S_G2xyz_S+ABY*I_ERI_Hy4z_S_G2xyz_S;
  Double I_ERI_H5z_Py_G2xyz_S = I_ERI_Iy5z_S_G2xyz_S+ABY*I_ERI_H5z_S_G2xyz_S;
  Double I_ERI_H4xz_Pz_G2xyz_S = I_ERI_I4x2z_S_G2xyz_S+ABZ*I_ERI_H4xz_S_G2xyz_S;
  Double I_ERI_H3xyz_Pz_G2xyz_S = I_ERI_I3xy2z_S_G2xyz_S+ABZ*I_ERI_H3xyz_S_G2xyz_S;
  Double I_ERI_H3x2z_Pz_G2xyz_S = I_ERI_I3x3z_S_G2xyz_S+ABZ*I_ERI_H3x2z_S_G2xyz_S;
  Double I_ERI_H2x2yz_Pz_G2xyz_S = I_ERI_I2x2y2z_S_G2xyz_S+ABZ*I_ERI_H2x2yz_S_G2xyz_S;
  Double I_ERI_H2xy2z_Pz_G2xyz_S = I_ERI_I2xy3z_S_G2xyz_S+ABZ*I_ERI_H2xy2z_S_G2xyz_S;
  Double I_ERI_H2x3z_Pz_G2xyz_S = I_ERI_I2x4z_S_G2xyz_S+ABZ*I_ERI_H2x3z_S_G2xyz_S;
  Double I_ERI_Hx3yz_Pz_G2xyz_S = I_ERI_Ix3y2z_S_G2xyz_S+ABZ*I_ERI_Hx3yz_S_G2xyz_S;
  Double I_ERI_Hx2y2z_Pz_G2xyz_S = I_ERI_Ix2y3z_S_G2xyz_S+ABZ*I_ERI_Hx2y2z_S_G2xyz_S;
  Double I_ERI_Hxy3z_Pz_G2xyz_S = I_ERI_Ixy4z_S_G2xyz_S+ABZ*I_ERI_Hxy3z_S_G2xyz_S;
  Double I_ERI_Hx4z_Pz_G2xyz_S = I_ERI_Ix5z_S_G2xyz_S+ABZ*I_ERI_Hx4z_S_G2xyz_S;
  Double I_ERI_H4yz_Pz_G2xyz_S = I_ERI_I4y2z_S_G2xyz_S+ABZ*I_ERI_H4yz_S_G2xyz_S;
  Double I_ERI_H3y2z_Pz_G2xyz_S = I_ERI_I3y3z_S_G2xyz_S+ABZ*I_ERI_H3y2z_S_G2xyz_S;
  Double I_ERI_H2y3z_Pz_G2xyz_S = I_ERI_I2y4z_S_G2xyz_S+ABZ*I_ERI_H2y3z_S_G2xyz_S;
  Double I_ERI_Hy4z_Pz_G2xyz_S = I_ERI_Iy5z_S_G2xyz_S+ABZ*I_ERI_Hy4z_S_G2xyz_S;
  Double I_ERI_H5z_Pz_G2xyz_S = I_ERI_I6z_S_G2xyz_S+ABZ*I_ERI_H5z_S_G2xyz_S;
  Double I_ERI_H5x_Px_G2x2z_S = I_ERI_I6x_S_G2x2z_S+ABX*I_ERI_H5x_S_G2x2z_S;
  Double I_ERI_H4xy_Px_G2x2z_S = I_ERI_I5xy_S_G2x2z_S+ABX*I_ERI_H4xy_S_G2x2z_S;
  Double I_ERI_H4xz_Px_G2x2z_S = I_ERI_I5xz_S_G2x2z_S+ABX*I_ERI_H4xz_S_G2x2z_S;
  Double I_ERI_H3x2y_Px_G2x2z_S = I_ERI_I4x2y_S_G2x2z_S+ABX*I_ERI_H3x2y_S_G2x2z_S;
  Double I_ERI_H3xyz_Px_G2x2z_S = I_ERI_I4xyz_S_G2x2z_S+ABX*I_ERI_H3xyz_S_G2x2z_S;
  Double I_ERI_H3x2z_Px_G2x2z_S = I_ERI_I4x2z_S_G2x2z_S+ABX*I_ERI_H3x2z_S_G2x2z_S;
  Double I_ERI_H2x3y_Px_G2x2z_S = I_ERI_I3x3y_S_G2x2z_S+ABX*I_ERI_H2x3y_S_G2x2z_S;
  Double I_ERI_H2x2yz_Px_G2x2z_S = I_ERI_I3x2yz_S_G2x2z_S+ABX*I_ERI_H2x2yz_S_G2x2z_S;
  Double I_ERI_H2xy2z_Px_G2x2z_S = I_ERI_I3xy2z_S_G2x2z_S+ABX*I_ERI_H2xy2z_S_G2x2z_S;
  Double I_ERI_H2x3z_Px_G2x2z_S = I_ERI_I3x3z_S_G2x2z_S+ABX*I_ERI_H2x3z_S_G2x2z_S;
  Double I_ERI_Hx4y_Px_G2x2z_S = I_ERI_I2x4y_S_G2x2z_S+ABX*I_ERI_Hx4y_S_G2x2z_S;
  Double I_ERI_Hx3yz_Px_G2x2z_S = I_ERI_I2x3yz_S_G2x2z_S+ABX*I_ERI_Hx3yz_S_G2x2z_S;
  Double I_ERI_Hx2y2z_Px_G2x2z_S = I_ERI_I2x2y2z_S_G2x2z_S+ABX*I_ERI_Hx2y2z_S_G2x2z_S;
  Double I_ERI_Hxy3z_Px_G2x2z_S = I_ERI_I2xy3z_S_G2x2z_S+ABX*I_ERI_Hxy3z_S_G2x2z_S;
  Double I_ERI_Hx4z_Px_G2x2z_S = I_ERI_I2x4z_S_G2x2z_S+ABX*I_ERI_Hx4z_S_G2x2z_S;
  Double I_ERI_H5y_Px_G2x2z_S = I_ERI_Ix5y_S_G2x2z_S+ABX*I_ERI_H5y_S_G2x2z_S;
  Double I_ERI_H4yz_Px_G2x2z_S = I_ERI_Ix4yz_S_G2x2z_S+ABX*I_ERI_H4yz_S_G2x2z_S;
  Double I_ERI_H3y2z_Px_G2x2z_S = I_ERI_Ix3y2z_S_G2x2z_S+ABX*I_ERI_H3y2z_S_G2x2z_S;
  Double I_ERI_H2y3z_Px_G2x2z_S = I_ERI_Ix2y3z_S_G2x2z_S+ABX*I_ERI_H2y3z_S_G2x2z_S;
  Double I_ERI_Hy4z_Px_G2x2z_S = I_ERI_Ixy4z_S_G2x2z_S+ABX*I_ERI_Hy4z_S_G2x2z_S;
  Double I_ERI_H5z_Px_G2x2z_S = I_ERI_Ix5z_S_G2x2z_S+ABX*I_ERI_H5z_S_G2x2z_S;
  Double I_ERI_H4xy_Py_G2x2z_S = I_ERI_I4x2y_S_G2x2z_S+ABY*I_ERI_H4xy_S_G2x2z_S;
  Double I_ERI_H4xz_Py_G2x2z_S = I_ERI_I4xyz_S_G2x2z_S+ABY*I_ERI_H4xz_S_G2x2z_S;
  Double I_ERI_H3x2y_Py_G2x2z_S = I_ERI_I3x3y_S_G2x2z_S+ABY*I_ERI_H3x2y_S_G2x2z_S;
  Double I_ERI_H3xyz_Py_G2x2z_S = I_ERI_I3x2yz_S_G2x2z_S+ABY*I_ERI_H3xyz_S_G2x2z_S;
  Double I_ERI_H3x2z_Py_G2x2z_S = I_ERI_I3xy2z_S_G2x2z_S+ABY*I_ERI_H3x2z_S_G2x2z_S;
  Double I_ERI_H2x3y_Py_G2x2z_S = I_ERI_I2x4y_S_G2x2z_S+ABY*I_ERI_H2x3y_S_G2x2z_S;
  Double I_ERI_H2x2yz_Py_G2x2z_S = I_ERI_I2x3yz_S_G2x2z_S+ABY*I_ERI_H2x2yz_S_G2x2z_S;
  Double I_ERI_H2xy2z_Py_G2x2z_S = I_ERI_I2x2y2z_S_G2x2z_S+ABY*I_ERI_H2xy2z_S_G2x2z_S;
  Double I_ERI_H2x3z_Py_G2x2z_S = I_ERI_I2xy3z_S_G2x2z_S+ABY*I_ERI_H2x3z_S_G2x2z_S;
  Double I_ERI_Hx4y_Py_G2x2z_S = I_ERI_Ix5y_S_G2x2z_S+ABY*I_ERI_Hx4y_S_G2x2z_S;
  Double I_ERI_Hx3yz_Py_G2x2z_S = I_ERI_Ix4yz_S_G2x2z_S+ABY*I_ERI_Hx3yz_S_G2x2z_S;
  Double I_ERI_Hx2y2z_Py_G2x2z_S = I_ERI_Ix3y2z_S_G2x2z_S+ABY*I_ERI_Hx2y2z_S_G2x2z_S;
  Double I_ERI_Hxy3z_Py_G2x2z_S = I_ERI_Ix2y3z_S_G2x2z_S+ABY*I_ERI_Hxy3z_S_G2x2z_S;
  Double I_ERI_Hx4z_Py_G2x2z_S = I_ERI_Ixy4z_S_G2x2z_S+ABY*I_ERI_Hx4z_S_G2x2z_S;
  Double I_ERI_H5y_Py_G2x2z_S = I_ERI_I6y_S_G2x2z_S+ABY*I_ERI_H5y_S_G2x2z_S;
  Double I_ERI_H4yz_Py_G2x2z_S = I_ERI_I5yz_S_G2x2z_S+ABY*I_ERI_H4yz_S_G2x2z_S;
  Double I_ERI_H3y2z_Py_G2x2z_S = I_ERI_I4y2z_S_G2x2z_S+ABY*I_ERI_H3y2z_S_G2x2z_S;
  Double I_ERI_H2y3z_Py_G2x2z_S = I_ERI_I3y3z_S_G2x2z_S+ABY*I_ERI_H2y3z_S_G2x2z_S;
  Double I_ERI_Hy4z_Py_G2x2z_S = I_ERI_I2y4z_S_G2x2z_S+ABY*I_ERI_Hy4z_S_G2x2z_S;
  Double I_ERI_H5z_Py_G2x2z_S = I_ERI_Iy5z_S_G2x2z_S+ABY*I_ERI_H5z_S_G2x2z_S;
  Double I_ERI_H4xz_Pz_G2x2z_S = I_ERI_I4x2z_S_G2x2z_S+ABZ*I_ERI_H4xz_S_G2x2z_S;
  Double I_ERI_H3xyz_Pz_G2x2z_S = I_ERI_I3xy2z_S_G2x2z_S+ABZ*I_ERI_H3xyz_S_G2x2z_S;
  Double I_ERI_H3x2z_Pz_G2x2z_S = I_ERI_I3x3z_S_G2x2z_S+ABZ*I_ERI_H3x2z_S_G2x2z_S;
  Double I_ERI_H2x2yz_Pz_G2x2z_S = I_ERI_I2x2y2z_S_G2x2z_S+ABZ*I_ERI_H2x2yz_S_G2x2z_S;
  Double I_ERI_H2xy2z_Pz_G2x2z_S = I_ERI_I2xy3z_S_G2x2z_S+ABZ*I_ERI_H2xy2z_S_G2x2z_S;
  Double I_ERI_H2x3z_Pz_G2x2z_S = I_ERI_I2x4z_S_G2x2z_S+ABZ*I_ERI_H2x3z_S_G2x2z_S;
  Double I_ERI_Hx3yz_Pz_G2x2z_S = I_ERI_Ix3y2z_S_G2x2z_S+ABZ*I_ERI_Hx3yz_S_G2x2z_S;
  Double I_ERI_Hx2y2z_Pz_G2x2z_S = I_ERI_Ix2y3z_S_G2x2z_S+ABZ*I_ERI_Hx2y2z_S_G2x2z_S;
  Double I_ERI_Hxy3z_Pz_G2x2z_S = I_ERI_Ixy4z_S_G2x2z_S+ABZ*I_ERI_Hxy3z_S_G2x2z_S;
  Double I_ERI_Hx4z_Pz_G2x2z_S = I_ERI_Ix5z_S_G2x2z_S+ABZ*I_ERI_Hx4z_S_G2x2z_S;
  Double I_ERI_H4yz_Pz_G2x2z_S = I_ERI_I4y2z_S_G2x2z_S+ABZ*I_ERI_H4yz_S_G2x2z_S;
  Double I_ERI_H3y2z_Pz_G2x2z_S = I_ERI_I3y3z_S_G2x2z_S+ABZ*I_ERI_H3y2z_S_G2x2z_S;
  Double I_ERI_H2y3z_Pz_G2x2z_S = I_ERI_I2y4z_S_G2x2z_S+ABZ*I_ERI_H2y3z_S_G2x2z_S;
  Double I_ERI_Hy4z_Pz_G2x2z_S = I_ERI_Iy5z_S_G2x2z_S+ABZ*I_ERI_Hy4z_S_G2x2z_S;
  Double I_ERI_H5z_Pz_G2x2z_S = I_ERI_I6z_S_G2x2z_S+ABZ*I_ERI_H5z_S_G2x2z_S;
  Double I_ERI_H5x_Px_Gx3y_S = I_ERI_I6x_S_Gx3y_S+ABX*I_ERI_H5x_S_Gx3y_S;
  Double I_ERI_H4xy_Px_Gx3y_S = I_ERI_I5xy_S_Gx3y_S+ABX*I_ERI_H4xy_S_Gx3y_S;
  Double I_ERI_H4xz_Px_Gx3y_S = I_ERI_I5xz_S_Gx3y_S+ABX*I_ERI_H4xz_S_Gx3y_S;
  Double I_ERI_H3x2y_Px_Gx3y_S = I_ERI_I4x2y_S_Gx3y_S+ABX*I_ERI_H3x2y_S_Gx3y_S;
  Double I_ERI_H3xyz_Px_Gx3y_S = I_ERI_I4xyz_S_Gx3y_S+ABX*I_ERI_H3xyz_S_Gx3y_S;
  Double I_ERI_H3x2z_Px_Gx3y_S = I_ERI_I4x2z_S_Gx3y_S+ABX*I_ERI_H3x2z_S_Gx3y_S;
  Double I_ERI_H2x3y_Px_Gx3y_S = I_ERI_I3x3y_S_Gx3y_S+ABX*I_ERI_H2x3y_S_Gx3y_S;
  Double I_ERI_H2x2yz_Px_Gx3y_S = I_ERI_I3x2yz_S_Gx3y_S+ABX*I_ERI_H2x2yz_S_Gx3y_S;
  Double I_ERI_H2xy2z_Px_Gx3y_S = I_ERI_I3xy2z_S_Gx3y_S+ABX*I_ERI_H2xy2z_S_Gx3y_S;
  Double I_ERI_H2x3z_Px_Gx3y_S = I_ERI_I3x3z_S_Gx3y_S+ABX*I_ERI_H2x3z_S_Gx3y_S;
  Double I_ERI_Hx4y_Px_Gx3y_S = I_ERI_I2x4y_S_Gx3y_S+ABX*I_ERI_Hx4y_S_Gx3y_S;
  Double I_ERI_Hx3yz_Px_Gx3y_S = I_ERI_I2x3yz_S_Gx3y_S+ABX*I_ERI_Hx3yz_S_Gx3y_S;
  Double I_ERI_Hx2y2z_Px_Gx3y_S = I_ERI_I2x2y2z_S_Gx3y_S+ABX*I_ERI_Hx2y2z_S_Gx3y_S;
  Double I_ERI_Hxy3z_Px_Gx3y_S = I_ERI_I2xy3z_S_Gx3y_S+ABX*I_ERI_Hxy3z_S_Gx3y_S;
  Double I_ERI_Hx4z_Px_Gx3y_S = I_ERI_I2x4z_S_Gx3y_S+ABX*I_ERI_Hx4z_S_Gx3y_S;
  Double I_ERI_H5y_Px_Gx3y_S = I_ERI_Ix5y_S_Gx3y_S+ABX*I_ERI_H5y_S_Gx3y_S;
  Double I_ERI_H4yz_Px_Gx3y_S = I_ERI_Ix4yz_S_Gx3y_S+ABX*I_ERI_H4yz_S_Gx3y_S;
  Double I_ERI_H3y2z_Px_Gx3y_S = I_ERI_Ix3y2z_S_Gx3y_S+ABX*I_ERI_H3y2z_S_Gx3y_S;
  Double I_ERI_H2y3z_Px_Gx3y_S = I_ERI_Ix2y3z_S_Gx3y_S+ABX*I_ERI_H2y3z_S_Gx3y_S;
  Double I_ERI_Hy4z_Px_Gx3y_S = I_ERI_Ixy4z_S_Gx3y_S+ABX*I_ERI_Hy4z_S_Gx3y_S;
  Double I_ERI_H5z_Px_Gx3y_S = I_ERI_Ix5z_S_Gx3y_S+ABX*I_ERI_H5z_S_Gx3y_S;
  Double I_ERI_H4xy_Py_Gx3y_S = I_ERI_I4x2y_S_Gx3y_S+ABY*I_ERI_H4xy_S_Gx3y_S;
  Double I_ERI_H4xz_Py_Gx3y_S = I_ERI_I4xyz_S_Gx3y_S+ABY*I_ERI_H4xz_S_Gx3y_S;
  Double I_ERI_H3x2y_Py_Gx3y_S = I_ERI_I3x3y_S_Gx3y_S+ABY*I_ERI_H3x2y_S_Gx3y_S;
  Double I_ERI_H3xyz_Py_Gx3y_S = I_ERI_I3x2yz_S_Gx3y_S+ABY*I_ERI_H3xyz_S_Gx3y_S;
  Double I_ERI_H3x2z_Py_Gx3y_S = I_ERI_I3xy2z_S_Gx3y_S+ABY*I_ERI_H3x2z_S_Gx3y_S;
  Double I_ERI_H2x3y_Py_Gx3y_S = I_ERI_I2x4y_S_Gx3y_S+ABY*I_ERI_H2x3y_S_Gx3y_S;
  Double I_ERI_H2x2yz_Py_Gx3y_S = I_ERI_I2x3yz_S_Gx3y_S+ABY*I_ERI_H2x2yz_S_Gx3y_S;
  Double I_ERI_H2xy2z_Py_Gx3y_S = I_ERI_I2x2y2z_S_Gx3y_S+ABY*I_ERI_H2xy2z_S_Gx3y_S;
  Double I_ERI_H2x3z_Py_Gx3y_S = I_ERI_I2xy3z_S_Gx3y_S+ABY*I_ERI_H2x3z_S_Gx3y_S;
  Double I_ERI_Hx4y_Py_Gx3y_S = I_ERI_Ix5y_S_Gx3y_S+ABY*I_ERI_Hx4y_S_Gx3y_S;
  Double I_ERI_Hx3yz_Py_Gx3y_S = I_ERI_Ix4yz_S_Gx3y_S+ABY*I_ERI_Hx3yz_S_Gx3y_S;
  Double I_ERI_Hx2y2z_Py_Gx3y_S = I_ERI_Ix3y2z_S_Gx3y_S+ABY*I_ERI_Hx2y2z_S_Gx3y_S;
  Double I_ERI_Hxy3z_Py_Gx3y_S = I_ERI_Ix2y3z_S_Gx3y_S+ABY*I_ERI_Hxy3z_S_Gx3y_S;
  Double I_ERI_Hx4z_Py_Gx3y_S = I_ERI_Ixy4z_S_Gx3y_S+ABY*I_ERI_Hx4z_S_Gx3y_S;
  Double I_ERI_H5y_Py_Gx3y_S = I_ERI_I6y_S_Gx3y_S+ABY*I_ERI_H5y_S_Gx3y_S;
  Double I_ERI_H4yz_Py_Gx3y_S = I_ERI_I5yz_S_Gx3y_S+ABY*I_ERI_H4yz_S_Gx3y_S;
  Double I_ERI_H3y2z_Py_Gx3y_S = I_ERI_I4y2z_S_Gx3y_S+ABY*I_ERI_H3y2z_S_Gx3y_S;
  Double I_ERI_H2y3z_Py_Gx3y_S = I_ERI_I3y3z_S_Gx3y_S+ABY*I_ERI_H2y3z_S_Gx3y_S;
  Double I_ERI_Hy4z_Py_Gx3y_S = I_ERI_I2y4z_S_Gx3y_S+ABY*I_ERI_Hy4z_S_Gx3y_S;
  Double I_ERI_H5z_Py_Gx3y_S = I_ERI_Iy5z_S_Gx3y_S+ABY*I_ERI_H5z_S_Gx3y_S;
  Double I_ERI_H4xz_Pz_Gx3y_S = I_ERI_I4x2z_S_Gx3y_S+ABZ*I_ERI_H4xz_S_Gx3y_S;
  Double I_ERI_H3xyz_Pz_Gx3y_S = I_ERI_I3xy2z_S_Gx3y_S+ABZ*I_ERI_H3xyz_S_Gx3y_S;
  Double I_ERI_H3x2z_Pz_Gx3y_S = I_ERI_I3x3z_S_Gx3y_S+ABZ*I_ERI_H3x2z_S_Gx3y_S;
  Double I_ERI_H2x2yz_Pz_Gx3y_S = I_ERI_I2x2y2z_S_Gx3y_S+ABZ*I_ERI_H2x2yz_S_Gx3y_S;
  Double I_ERI_H2xy2z_Pz_Gx3y_S = I_ERI_I2xy3z_S_Gx3y_S+ABZ*I_ERI_H2xy2z_S_Gx3y_S;
  Double I_ERI_H2x3z_Pz_Gx3y_S = I_ERI_I2x4z_S_Gx3y_S+ABZ*I_ERI_H2x3z_S_Gx3y_S;
  Double I_ERI_Hx3yz_Pz_Gx3y_S = I_ERI_Ix3y2z_S_Gx3y_S+ABZ*I_ERI_Hx3yz_S_Gx3y_S;
  Double I_ERI_Hx2y2z_Pz_Gx3y_S = I_ERI_Ix2y3z_S_Gx3y_S+ABZ*I_ERI_Hx2y2z_S_Gx3y_S;
  Double I_ERI_Hxy3z_Pz_Gx3y_S = I_ERI_Ixy4z_S_Gx3y_S+ABZ*I_ERI_Hxy3z_S_Gx3y_S;
  Double I_ERI_Hx4z_Pz_Gx3y_S = I_ERI_Ix5z_S_Gx3y_S+ABZ*I_ERI_Hx4z_S_Gx3y_S;
  Double I_ERI_H4yz_Pz_Gx3y_S = I_ERI_I4y2z_S_Gx3y_S+ABZ*I_ERI_H4yz_S_Gx3y_S;
  Double I_ERI_H3y2z_Pz_Gx3y_S = I_ERI_I3y3z_S_Gx3y_S+ABZ*I_ERI_H3y2z_S_Gx3y_S;
  Double I_ERI_H2y3z_Pz_Gx3y_S = I_ERI_I2y4z_S_Gx3y_S+ABZ*I_ERI_H2y3z_S_Gx3y_S;
  Double I_ERI_Hy4z_Pz_Gx3y_S = I_ERI_Iy5z_S_Gx3y_S+ABZ*I_ERI_Hy4z_S_Gx3y_S;
  Double I_ERI_H5z_Pz_Gx3y_S = I_ERI_I6z_S_Gx3y_S+ABZ*I_ERI_H5z_S_Gx3y_S;
  Double I_ERI_H5x_Px_Gx2yz_S = I_ERI_I6x_S_Gx2yz_S+ABX*I_ERI_H5x_S_Gx2yz_S;
  Double I_ERI_H4xy_Px_Gx2yz_S = I_ERI_I5xy_S_Gx2yz_S+ABX*I_ERI_H4xy_S_Gx2yz_S;
  Double I_ERI_H4xz_Px_Gx2yz_S = I_ERI_I5xz_S_Gx2yz_S+ABX*I_ERI_H4xz_S_Gx2yz_S;
  Double I_ERI_H3x2y_Px_Gx2yz_S = I_ERI_I4x2y_S_Gx2yz_S+ABX*I_ERI_H3x2y_S_Gx2yz_S;
  Double I_ERI_H3xyz_Px_Gx2yz_S = I_ERI_I4xyz_S_Gx2yz_S+ABX*I_ERI_H3xyz_S_Gx2yz_S;
  Double I_ERI_H3x2z_Px_Gx2yz_S = I_ERI_I4x2z_S_Gx2yz_S+ABX*I_ERI_H3x2z_S_Gx2yz_S;
  Double I_ERI_H2x3y_Px_Gx2yz_S = I_ERI_I3x3y_S_Gx2yz_S+ABX*I_ERI_H2x3y_S_Gx2yz_S;
  Double I_ERI_H2x2yz_Px_Gx2yz_S = I_ERI_I3x2yz_S_Gx2yz_S+ABX*I_ERI_H2x2yz_S_Gx2yz_S;
  Double I_ERI_H2xy2z_Px_Gx2yz_S = I_ERI_I3xy2z_S_Gx2yz_S+ABX*I_ERI_H2xy2z_S_Gx2yz_S;
  Double I_ERI_H2x3z_Px_Gx2yz_S = I_ERI_I3x3z_S_Gx2yz_S+ABX*I_ERI_H2x3z_S_Gx2yz_S;
  Double I_ERI_Hx4y_Px_Gx2yz_S = I_ERI_I2x4y_S_Gx2yz_S+ABX*I_ERI_Hx4y_S_Gx2yz_S;
  Double I_ERI_Hx3yz_Px_Gx2yz_S = I_ERI_I2x3yz_S_Gx2yz_S+ABX*I_ERI_Hx3yz_S_Gx2yz_S;
  Double I_ERI_Hx2y2z_Px_Gx2yz_S = I_ERI_I2x2y2z_S_Gx2yz_S+ABX*I_ERI_Hx2y2z_S_Gx2yz_S;
  Double I_ERI_Hxy3z_Px_Gx2yz_S = I_ERI_I2xy3z_S_Gx2yz_S+ABX*I_ERI_Hxy3z_S_Gx2yz_S;
  Double I_ERI_Hx4z_Px_Gx2yz_S = I_ERI_I2x4z_S_Gx2yz_S+ABX*I_ERI_Hx4z_S_Gx2yz_S;
  Double I_ERI_H5y_Px_Gx2yz_S = I_ERI_Ix5y_S_Gx2yz_S+ABX*I_ERI_H5y_S_Gx2yz_S;
  Double I_ERI_H4yz_Px_Gx2yz_S = I_ERI_Ix4yz_S_Gx2yz_S+ABX*I_ERI_H4yz_S_Gx2yz_S;
  Double I_ERI_H3y2z_Px_Gx2yz_S = I_ERI_Ix3y2z_S_Gx2yz_S+ABX*I_ERI_H3y2z_S_Gx2yz_S;
  Double I_ERI_H2y3z_Px_Gx2yz_S = I_ERI_Ix2y3z_S_Gx2yz_S+ABX*I_ERI_H2y3z_S_Gx2yz_S;
  Double I_ERI_Hy4z_Px_Gx2yz_S = I_ERI_Ixy4z_S_Gx2yz_S+ABX*I_ERI_Hy4z_S_Gx2yz_S;
  Double I_ERI_H5z_Px_Gx2yz_S = I_ERI_Ix5z_S_Gx2yz_S+ABX*I_ERI_H5z_S_Gx2yz_S;
  Double I_ERI_H4xy_Py_Gx2yz_S = I_ERI_I4x2y_S_Gx2yz_S+ABY*I_ERI_H4xy_S_Gx2yz_S;
  Double I_ERI_H4xz_Py_Gx2yz_S = I_ERI_I4xyz_S_Gx2yz_S+ABY*I_ERI_H4xz_S_Gx2yz_S;
  Double I_ERI_H3x2y_Py_Gx2yz_S = I_ERI_I3x3y_S_Gx2yz_S+ABY*I_ERI_H3x2y_S_Gx2yz_S;
  Double I_ERI_H3xyz_Py_Gx2yz_S = I_ERI_I3x2yz_S_Gx2yz_S+ABY*I_ERI_H3xyz_S_Gx2yz_S;
  Double I_ERI_H3x2z_Py_Gx2yz_S = I_ERI_I3xy2z_S_Gx2yz_S+ABY*I_ERI_H3x2z_S_Gx2yz_S;
  Double I_ERI_H2x3y_Py_Gx2yz_S = I_ERI_I2x4y_S_Gx2yz_S+ABY*I_ERI_H2x3y_S_Gx2yz_S;
  Double I_ERI_H2x2yz_Py_Gx2yz_S = I_ERI_I2x3yz_S_Gx2yz_S+ABY*I_ERI_H2x2yz_S_Gx2yz_S;
  Double I_ERI_H2xy2z_Py_Gx2yz_S = I_ERI_I2x2y2z_S_Gx2yz_S+ABY*I_ERI_H2xy2z_S_Gx2yz_S;
  Double I_ERI_H2x3z_Py_Gx2yz_S = I_ERI_I2xy3z_S_Gx2yz_S+ABY*I_ERI_H2x3z_S_Gx2yz_S;
  Double I_ERI_Hx4y_Py_Gx2yz_S = I_ERI_Ix5y_S_Gx2yz_S+ABY*I_ERI_Hx4y_S_Gx2yz_S;
  Double I_ERI_Hx3yz_Py_Gx2yz_S = I_ERI_Ix4yz_S_Gx2yz_S+ABY*I_ERI_Hx3yz_S_Gx2yz_S;
  Double I_ERI_Hx2y2z_Py_Gx2yz_S = I_ERI_Ix3y2z_S_Gx2yz_S+ABY*I_ERI_Hx2y2z_S_Gx2yz_S;
  Double I_ERI_Hxy3z_Py_Gx2yz_S = I_ERI_Ix2y3z_S_Gx2yz_S+ABY*I_ERI_Hxy3z_S_Gx2yz_S;
  Double I_ERI_Hx4z_Py_Gx2yz_S = I_ERI_Ixy4z_S_Gx2yz_S+ABY*I_ERI_Hx4z_S_Gx2yz_S;
  Double I_ERI_H5y_Py_Gx2yz_S = I_ERI_I6y_S_Gx2yz_S+ABY*I_ERI_H5y_S_Gx2yz_S;
  Double I_ERI_H4yz_Py_Gx2yz_S = I_ERI_I5yz_S_Gx2yz_S+ABY*I_ERI_H4yz_S_Gx2yz_S;
  Double I_ERI_H3y2z_Py_Gx2yz_S = I_ERI_I4y2z_S_Gx2yz_S+ABY*I_ERI_H3y2z_S_Gx2yz_S;
  Double I_ERI_H2y3z_Py_Gx2yz_S = I_ERI_I3y3z_S_Gx2yz_S+ABY*I_ERI_H2y3z_S_Gx2yz_S;
  Double I_ERI_Hy4z_Py_Gx2yz_S = I_ERI_I2y4z_S_Gx2yz_S+ABY*I_ERI_Hy4z_S_Gx2yz_S;
  Double I_ERI_H5z_Py_Gx2yz_S = I_ERI_Iy5z_S_Gx2yz_S+ABY*I_ERI_H5z_S_Gx2yz_S;
  Double I_ERI_H4xz_Pz_Gx2yz_S = I_ERI_I4x2z_S_Gx2yz_S+ABZ*I_ERI_H4xz_S_Gx2yz_S;
  Double I_ERI_H3xyz_Pz_Gx2yz_S = I_ERI_I3xy2z_S_Gx2yz_S+ABZ*I_ERI_H3xyz_S_Gx2yz_S;
  Double I_ERI_H3x2z_Pz_Gx2yz_S = I_ERI_I3x3z_S_Gx2yz_S+ABZ*I_ERI_H3x2z_S_Gx2yz_S;
  Double I_ERI_H2x2yz_Pz_Gx2yz_S = I_ERI_I2x2y2z_S_Gx2yz_S+ABZ*I_ERI_H2x2yz_S_Gx2yz_S;
  Double I_ERI_H2xy2z_Pz_Gx2yz_S = I_ERI_I2xy3z_S_Gx2yz_S+ABZ*I_ERI_H2xy2z_S_Gx2yz_S;
  Double I_ERI_H2x3z_Pz_Gx2yz_S = I_ERI_I2x4z_S_Gx2yz_S+ABZ*I_ERI_H2x3z_S_Gx2yz_S;
  Double I_ERI_Hx3yz_Pz_Gx2yz_S = I_ERI_Ix3y2z_S_Gx2yz_S+ABZ*I_ERI_Hx3yz_S_Gx2yz_S;
  Double I_ERI_Hx2y2z_Pz_Gx2yz_S = I_ERI_Ix2y3z_S_Gx2yz_S+ABZ*I_ERI_Hx2y2z_S_Gx2yz_S;
  Double I_ERI_Hxy3z_Pz_Gx2yz_S = I_ERI_Ixy4z_S_Gx2yz_S+ABZ*I_ERI_Hxy3z_S_Gx2yz_S;
  Double I_ERI_Hx4z_Pz_Gx2yz_S = I_ERI_Ix5z_S_Gx2yz_S+ABZ*I_ERI_Hx4z_S_Gx2yz_S;
  Double I_ERI_H4yz_Pz_Gx2yz_S = I_ERI_I4y2z_S_Gx2yz_S+ABZ*I_ERI_H4yz_S_Gx2yz_S;
  Double I_ERI_H3y2z_Pz_Gx2yz_S = I_ERI_I3y3z_S_Gx2yz_S+ABZ*I_ERI_H3y2z_S_Gx2yz_S;
  Double I_ERI_H2y3z_Pz_Gx2yz_S = I_ERI_I2y4z_S_Gx2yz_S+ABZ*I_ERI_H2y3z_S_Gx2yz_S;
  Double I_ERI_Hy4z_Pz_Gx2yz_S = I_ERI_Iy5z_S_Gx2yz_S+ABZ*I_ERI_Hy4z_S_Gx2yz_S;
  Double I_ERI_H5z_Pz_Gx2yz_S = I_ERI_I6z_S_Gx2yz_S+ABZ*I_ERI_H5z_S_Gx2yz_S;
  Double I_ERI_H5x_Px_Gxy2z_S = I_ERI_I6x_S_Gxy2z_S+ABX*I_ERI_H5x_S_Gxy2z_S;
  Double I_ERI_H4xy_Px_Gxy2z_S = I_ERI_I5xy_S_Gxy2z_S+ABX*I_ERI_H4xy_S_Gxy2z_S;
  Double I_ERI_H4xz_Px_Gxy2z_S = I_ERI_I5xz_S_Gxy2z_S+ABX*I_ERI_H4xz_S_Gxy2z_S;
  Double I_ERI_H3x2y_Px_Gxy2z_S = I_ERI_I4x2y_S_Gxy2z_S+ABX*I_ERI_H3x2y_S_Gxy2z_S;
  Double I_ERI_H3xyz_Px_Gxy2z_S = I_ERI_I4xyz_S_Gxy2z_S+ABX*I_ERI_H3xyz_S_Gxy2z_S;
  Double I_ERI_H3x2z_Px_Gxy2z_S = I_ERI_I4x2z_S_Gxy2z_S+ABX*I_ERI_H3x2z_S_Gxy2z_S;
  Double I_ERI_H2x3y_Px_Gxy2z_S = I_ERI_I3x3y_S_Gxy2z_S+ABX*I_ERI_H2x3y_S_Gxy2z_S;
  Double I_ERI_H2x2yz_Px_Gxy2z_S = I_ERI_I3x2yz_S_Gxy2z_S+ABX*I_ERI_H2x2yz_S_Gxy2z_S;
  Double I_ERI_H2xy2z_Px_Gxy2z_S = I_ERI_I3xy2z_S_Gxy2z_S+ABX*I_ERI_H2xy2z_S_Gxy2z_S;
  Double I_ERI_H2x3z_Px_Gxy2z_S = I_ERI_I3x3z_S_Gxy2z_S+ABX*I_ERI_H2x3z_S_Gxy2z_S;
  Double I_ERI_Hx4y_Px_Gxy2z_S = I_ERI_I2x4y_S_Gxy2z_S+ABX*I_ERI_Hx4y_S_Gxy2z_S;
  Double I_ERI_Hx3yz_Px_Gxy2z_S = I_ERI_I2x3yz_S_Gxy2z_S+ABX*I_ERI_Hx3yz_S_Gxy2z_S;
  Double I_ERI_Hx2y2z_Px_Gxy2z_S = I_ERI_I2x2y2z_S_Gxy2z_S+ABX*I_ERI_Hx2y2z_S_Gxy2z_S;
  Double I_ERI_Hxy3z_Px_Gxy2z_S = I_ERI_I2xy3z_S_Gxy2z_S+ABX*I_ERI_Hxy3z_S_Gxy2z_S;
  Double I_ERI_Hx4z_Px_Gxy2z_S = I_ERI_I2x4z_S_Gxy2z_S+ABX*I_ERI_Hx4z_S_Gxy2z_S;
  Double I_ERI_H5y_Px_Gxy2z_S = I_ERI_Ix5y_S_Gxy2z_S+ABX*I_ERI_H5y_S_Gxy2z_S;
  Double I_ERI_H4yz_Px_Gxy2z_S = I_ERI_Ix4yz_S_Gxy2z_S+ABX*I_ERI_H4yz_S_Gxy2z_S;
  Double I_ERI_H3y2z_Px_Gxy2z_S = I_ERI_Ix3y2z_S_Gxy2z_S+ABX*I_ERI_H3y2z_S_Gxy2z_S;
  Double I_ERI_H2y3z_Px_Gxy2z_S = I_ERI_Ix2y3z_S_Gxy2z_S+ABX*I_ERI_H2y3z_S_Gxy2z_S;
  Double I_ERI_Hy4z_Px_Gxy2z_S = I_ERI_Ixy4z_S_Gxy2z_S+ABX*I_ERI_Hy4z_S_Gxy2z_S;
  Double I_ERI_H5z_Px_Gxy2z_S = I_ERI_Ix5z_S_Gxy2z_S+ABX*I_ERI_H5z_S_Gxy2z_S;
  Double I_ERI_H4xy_Py_Gxy2z_S = I_ERI_I4x2y_S_Gxy2z_S+ABY*I_ERI_H4xy_S_Gxy2z_S;
  Double I_ERI_H4xz_Py_Gxy2z_S = I_ERI_I4xyz_S_Gxy2z_S+ABY*I_ERI_H4xz_S_Gxy2z_S;
  Double I_ERI_H3x2y_Py_Gxy2z_S = I_ERI_I3x3y_S_Gxy2z_S+ABY*I_ERI_H3x2y_S_Gxy2z_S;
  Double I_ERI_H3xyz_Py_Gxy2z_S = I_ERI_I3x2yz_S_Gxy2z_S+ABY*I_ERI_H3xyz_S_Gxy2z_S;
  Double I_ERI_H3x2z_Py_Gxy2z_S = I_ERI_I3xy2z_S_Gxy2z_S+ABY*I_ERI_H3x2z_S_Gxy2z_S;
  Double I_ERI_H2x3y_Py_Gxy2z_S = I_ERI_I2x4y_S_Gxy2z_S+ABY*I_ERI_H2x3y_S_Gxy2z_S;
  Double I_ERI_H2x2yz_Py_Gxy2z_S = I_ERI_I2x3yz_S_Gxy2z_S+ABY*I_ERI_H2x2yz_S_Gxy2z_S;
  Double I_ERI_H2xy2z_Py_Gxy2z_S = I_ERI_I2x2y2z_S_Gxy2z_S+ABY*I_ERI_H2xy2z_S_Gxy2z_S;
  Double I_ERI_H2x3z_Py_Gxy2z_S = I_ERI_I2xy3z_S_Gxy2z_S+ABY*I_ERI_H2x3z_S_Gxy2z_S;
  Double I_ERI_Hx4y_Py_Gxy2z_S = I_ERI_Ix5y_S_Gxy2z_S+ABY*I_ERI_Hx4y_S_Gxy2z_S;
  Double I_ERI_Hx3yz_Py_Gxy2z_S = I_ERI_Ix4yz_S_Gxy2z_S+ABY*I_ERI_Hx3yz_S_Gxy2z_S;
  Double I_ERI_Hx2y2z_Py_Gxy2z_S = I_ERI_Ix3y2z_S_Gxy2z_S+ABY*I_ERI_Hx2y2z_S_Gxy2z_S;
  Double I_ERI_Hxy3z_Py_Gxy2z_S = I_ERI_Ix2y3z_S_Gxy2z_S+ABY*I_ERI_Hxy3z_S_Gxy2z_S;
  Double I_ERI_Hx4z_Py_Gxy2z_S = I_ERI_Ixy4z_S_Gxy2z_S+ABY*I_ERI_Hx4z_S_Gxy2z_S;
  Double I_ERI_H5y_Py_Gxy2z_S = I_ERI_I6y_S_Gxy2z_S+ABY*I_ERI_H5y_S_Gxy2z_S;
  Double I_ERI_H4yz_Py_Gxy2z_S = I_ERI_I5yz_S_Gxy2z_S+ABY*I_ERI_H4yz_S_Gxy2z_S;
  Double I_ERI_H3y2z_Py_Gxy2z_S = I_ERI_I4y2z_S_Gxy2z_S+ABY*I_ERI_H3y2z_S_Gxy2z_S;
  Double I_ERI_H2y3z_Py_Gxy2z_S = I_ERI_I3y3z_S_Gxy2z_S+ABY*I_ERI_H2y3z_S_Gxy2z_S;
  Double I_ERI_Hy4z_Py_Gxy2z_S = I_ERI_I2y4z_S_Gxy2z_S+ABY*I_ERI_Hy4z_S_Gxy2z_S;
  Double I_ERI_H5z_Py_Gxy2z_S = I_ERI_Iy5z_S_Gxy2z_S+ABY*I_ERI_H5z_S_Gxy2z_S;
  Double I_ERI_H4xz_Pz_Gxy2z_S = I_ERI_I4x2z_S_Gxy2z_S+ABZ*I_ERI_H4xz_S_Gxy2z_S;
  Double I_ERI_H3xyz_Pz_Gxy2z_S = I_ERI_I3xy2z_S_Gxy2z_S+ABZ*I_ERI_H3xyz_S_Gxy2z_S;
  Double I_ERI_H3x2z_Pz_Gxy2z_S = I_ERI_I3x3z_S_Gxy2z_S+ABZ*I_ERI_H3x2z_S_Gxy2z_S;
  Double I_ERI_H2x2yz_Pz_Gxy2z_S = I_ERI_I2x2y2z_S_Gxy2z_S+ABZ*I_ERI_H2x2yz_S_Gxy2z_S;
  Double I_ERI_H2xy2z_Pz_Gxy2z_S = I_ERI_I2xy3z_S_Gxy2z_S+ABZ*I_ERI_H2xy2z_S_Gxy2z_S;
  Double I_ERI_H2x3z_Pz_Gxy2z_S = I_ERI_I2x4z_S_Gxy2z_S+ABZ*I_ERI_H2x3z_S_Gxy2z_S;
  Double I_ERI_Hx3yz_Pz_Gxy2z_S = I_ERI_Ix3y2z_S_Gxy2z_S+ABZ*I_ERI_Hx3yz_S_Gxy2z_S;
  Double I_ERI_Hx2y2z_Pz_Gxy2z_S = I_ERI_Ix2y3z_S_Gxy2z_S+ABZ*I_ERI_Hx2y2z_S_Gxy2z_S;
  Double I_ERI_Hxy3z_Pz_Gxy2z_S = I_ERI_Ixy4z_S_Gxy2z_S+ABZ*I_ERI_Hxy3z_S_Gxy2z_S;
  Double I_ERI_Hx4z_Pz_Gxy2z_S = I_ERI_Ix5z_S_Gxy2z_S+ABZ*I_ERI_Hx4z_S_Gxy2z_S;
  Double I_ERI_H4yz_Pz_Gxy2z_S = I_ERI_I4y2z_S_Gxy2z_S+ABZ*I_ERI_H4yz_S_Gxy2z_S;
  Double I_ERI_H3y2z_Pz_Gxy2z_S = I_ERI_I3y3z_S_Gxy2z_S+ABZ*I_ERI_H3y2z_S_Gxy2z_S;
  Double I_ERI_H2y3z_Pz_Gxy2z_S = I_ERI_I2y4z_S_Gxy2z_S+ABZ*I_ERI_H2y3z_S_Gxy2z_S;
  Double I_ERI_Hy4z_Pz_Gxy2z_S = I_ERI_Iy5z_S_Gxy2z_S+ABZ*I_ERI_Hy4z_S_Gxy2z_S;
  Double I_ERI_H5z_Pz_Gxy2z_S = I_ERI_I6z_S_Gxy2z_S+ABZ*I_ERI_H5z_S_Gxy2z_S;
  Double I_ERI_H5x_Px_Gx3z_S = I_ERI_I6x_S_Gx3z_S+ABX*I_ERI_H5x_S_Gx3z_S;
  Double I_ERI_H4xy_Px_Gx3z_S = I_ERI_I5xy_S_Gx3z_S+ABX*I_ERI_H4xy_S_Gx3z_S;
  Double I_ERI_H4xz_Px_Gx3z_S = I_ERI_I5xz_S_Gx3z_S+ABX*I_ERI_H4xz_S_Gx3z_S;
  Double I_ERI_H3x2y_Px_Gx3z_S = I_ERI_I4x2y_S_Gx3z_S+ABX*I_ERI_H3x2y_S_Gx3z_S;
  Double I_ERI_H3xyz_Px_Gx3z_S = I_ERI_I4xyz_S_Gx3z_S+ABX*I_ERI_H3xyz_S_Gx3z_S;
  Double I_ERI_H3x2z_Px_Gx3z_S = I_ERI_I4x2z_S_Gx3z_S+ABX*I_ERI_H3x2z_S_Gx3z_S;
  Double I_ERI_H2x3y_Px_Gx3z_S = I_ERI_I3x3y_S_Gx3z_S+ABX*I_ERI_H2x3y_S_Gx3z_S;
  Double I_ERI_H2x2yz_Px_Gx3z_S = I_ERI_I3x2yz_S_Gx3z_S+ABX*I_ERI_H2x2yz_S_Gx3z_S;
  Double I_ERI_H2xy2z_Px_Gx3z_S = I_ERI_I3xy2z_S_Gx3z_S+ABX*I_ERI_H2xy2z_S_Gx3z_S;
  Double I_ERI_H2x3z_Px_Gx3z_S = I_ERI_I3x3z_S_Gx3z_S+ABX*I_ERI_H2x3z_S_Gx3z_S;
  Double I_ERI_Hx4y_Px_Gx3z_S = I_ERI_I2x4y_S_Gx3z_S+ABX*I_ERI_Hx4y_S_Gx3z_S;
  Double I_ERI_Hx3yz_Px_Gx3z_S = I_ERI_I2x3yz_S_Gx3z_S+ABX*I_ERI_Hx3yz_S_Gx3z_S;
  Double I_ERI_Hx2y2z_Px_Gx3z_S = I_ERI_I2x2y2z_S_Gx3z_S+ABX*I_ERI_Hx2y2z_S_Gx3z_S;
  Double I_ERI_Hxy3z_Px_Gx3z_S = I_ERI_I2xy3z_S_Gx3z_S+ABX*I_ERI_Hxy3z_S_Gx3z_S;
  Double I_ERI_Hx4z_Px_Gx3z_S = I_ERI_I2x4z_S_Gx3z_S+ABX*I_ERI_Hx4z_S_Gx3z_S;
  Double I_ERI_H5y_Px_Gx3z_S = I_ERI_Ix5y_S_Gx3z_S+ABX*I_ERI_H5y_S_Gx3z_S;
  Double I_ERI_H4yz_Px_Gx3z_S = I_ERI_Ix4yz_S_Gx3z_S+ABX*I_ERI_H4yz_S_Gx3z_S;
  Double I_ERI_H3y2z_Px_Gx3z_S = I_ERI_Ix3y2z_S_Gx3z_S+ABX*I_ERI_H3y2z_S_Gx3z_S;
  Double I_ERI_H2y3z_Px_Gx3z_S = I_ERI_Ix2y3z_S_Gx3z_S+ABX*I_ERI_H2y3z_S_Gx3z_S;
  Double I_ERI_Hy4z_Px_Gx3z_S = I_ERI_Ixy4z_S_Gx3z_S+ABX*I_ERI_Hy4z_S_Gx3z_S;
  Double I_ERI_H5z_Px_Gx3z_S = I_ERI_Ix5z_S_Gx3z_S+ABX*I_ERI_H5z_S_Gx3z_S;
  Double I_ERI_H4xy_Py_Gx3z_S = I_ERI_I4x2y_S_Gx3z_S+ABY*I_ERI_H4xy_S_Gx3z_S;
  Double I_ERI_H4xz_Py_Gx3z_S = I_ERI_I4xyz_S_Gx3z_S+ABY*I_ERI_H4xz_S_Gx3z_S;
  Double I_ERI_H3x2y_Py_Gx3z_S = I_ERI_I3x3y_S_Gx3z_S+ABY*I_ERI_H3x2y_S_Gx3z_S;
  Double I_ERI_H3xyz_Py_Gx3z_S = I_ERI_I3x2yz_S_Gx3z_S+ABY*I_ERI_H3xyz_S_Gx3z_S;
  Double I_ERI_H3x2z_Py_Gx3z_S = I_ERI_I3xy2z_S_Gx3z_S+ABY*I_ERI_H3x2z_S_Gx3z_S;
  Double I_ERI_H2x3y_Py_Gx3z_S = I_ERI_I2x4y_S_Gx3z_S+ABY*I_ERI_H2x3y_S_Gx3z_S;
  Double I_ERI_H2x2yz_Py_Gx3z_S = I_ERI_I2x3yz_S_Gx3z_S+ABY*I_ERI_H2x2yz_S_Gx3z_S;
  Double I_ERI_H2xy2z_Py_Gx3z_S = I_ERI_I2x2y2z_S_Gx3z_S+ABY*I_ERI_H2xy2z_S_Gx3z_S;
  Double I_ERI_H2x3z_Py_Gx3z_S = I_ERI_I2xy3z_S_Gx3z_S+ABY*I_ERI_H2x3z_S_Gx3z_S;
  Double I_ERI_Hx4y_Py_Gx3z_S = I_ERI_Ix5y_S_Gx3z_S+ABY*I_ERI_Hx4y_S_Gx3z_S;
  Double I_ERI_Hx3yz_Py_Gx3z_S = I_ERI_Ix4yz_S_Gx3z_S+ABY*I_ERI_Hx3yz_S_Gx3z_S;
  Double I_ERI_Hx2y2z_Py_Gx3z_S = I_ERI_Ix3y2z_S_Gx3z_S+ABY*I_ERI_Hx2y2z_S_Gx3z_S;
  Double I_ERI_Hxy3z_Py_Gx3z_S = I_ERI_Ix2y3z_S_Gx3z_S+ABY*I_ERI_Hxy3z_S_Gx3z_S;
  Double I_ERI_Hx4z_Py_Gx3z_S = I_ERI_Ixy4z_S_Gx3z_S+ABY*I_ERI_Hx4z_S_Gx3z_S;
  Double I_ERI_H5y_Py_Gx3z_S = I_ERI_I6y_S_Gx3z_S+ABY*I_ERI_H5y_S_Gx3z_S;
  Double I_ERI_H4yz_Py_Gx3z_S = I_ERI_I5yz_S_Gx3z_S+ABY*I_ERI_H4yz_S_Gx3z_S;
  Double I_ERI_H3y2z_Py_Gx3z_S = I_ERI_I4y2z_S_Gx3z_S+ABY*I_ERI_H3y2z_S_Gx3z_S;
  Double I_ERI_H2y3z_Py_Gx3z_S = I_ERI_I3y3z_S_Gx3z_S+ABY*I_ERI_H2y3z_S_Gx3z_S;
  Double I_ERI_Hy4z_Py_Gx3z_S = I_ERI_I2y4z_S_Gx3z_S+ABY*I_ERI_Hy4z_S_Gx3z_S;
  Double I_ERI_H5z_Py_Gx3z_S = I_ERI_Iy5z_S_Gx3z_S+ABY*I_ERI_H5z_S_Gx3z_S;
  Double I_ERI_H4xz_Pz_Gx3z_S = I_ERI_I4x2z_S_Gx3z_S+ABZ*I_ERI_H4xz_S_Gx3z_S;
  Double I_ERI_H3xyz_Pz_Gx3z_S = I_ERI_I3xy2z_S_Gx3z_S+ABZ*I_ERI_H3xyz_S_Gx3z_S;
  Double I_ERI_H3x2z_Pz_Gx3z_S = I_ERI_I3x3z_S_Gx3z_S+ABZ*I_ERI_H3x2z_S_Gx3z_S;
  Double I_ERI_H2x2yz_Pz_Gx3z_S = I_ERI_I2x2y2z_S_Gx3z_S+ABZ*I_ERI_H2x2yz_S_Gx3z_S;
  Double I_ERI_H2xy2z_Pz_Gx3z_S = I_ERI_I2xy3z_S_Gx3z_S+ABZ*I_ERI_H2xy2z_S_Gx3z_S;
  Double I_ERI_H2x3z_Pz_Gx3z_S = I_ERI_I2x4z_S_Gx3z_S+ABZ*I_ERI_H2x3z_S_Gx3z_S;
  Double I_ERI_Hx3yz_Pz_Gx3z_S = I_ERI_Ix3y2z_S_Gx3z_S+ABZ*I_ERI_Hx3yz_S_Gx3z_S;
  Double I_ERI_Hx2y2z_Pz_Gx3z_S = I_ERI_Ix2y3z_S_Gx3z_S+ABZ*I_ERI_Hx2y2z_S_Gx3z_S;
  Double I_ERI_Hxy3z_Pz_Gx3z_S = I_ERI_Ixy4z_S_Gx3z_S+ABZ*I_ERI_Hxy3z_S_Gx3z_S;
  Double I_ERI_Hx4z_Pz_Gx3z_S = I_ERI_Ix5z_S_Gx3z_S+ABZ*I_ERI_Hx4z_S_Gx3z_S;
  Double I_ERI_H4yz_Pz_Gx3z_S = I_ERI_I4y2z_S_Gx3z_S+ABZ*I_ERI_H4yz_S_Gx3z_S;
  Double I_ERI_H3y2z_Pz_Gx3z_S = I_ERI_I3y3z_S_Gx3z_S+ABZ*I_ERI_H3y2z_S_Gx3z_S;
  Double I_ERI_H2y3z_Pz_Gx3z_S = I_ERI_I2y4z_S_Gx3z_S+ABZ*I_ERI_H2y3z_S_Gx3z_S;
  Double I_ERI_Hy4z_Pz_Gx3z_S = I_ERI_Iy5z_S_Gx3z_S+ABZ*I_ERI_Hy4z_S_Gx3z_S;
  Double I_ERI_H5z_Pz_Gx3z_S = I_ERI_I6z_S_Gx3z_S+ABZ*I_ERI_H5z_S_Gx3z_S;
  Double I_ERI_H5x_Px_G4y_S = I_ERI_I6x_S_G4y_S+ABX*I_ERI_H5x_S_G4y_S;
  Double I_ERI_H4xy_Px_G4y_S = I_ERI_I5xy_S_G4y_S+ABX*I_ERI_H4xy_S_G4y_S;
  Double I_ERI_H4xz_Px_G4y_S = I_ERI_I5xz_S_G4y_S+ABX*I_ERI_H4xz_S_G4y_S;
  Double I_ERI_H3x2y_Px_G4y_S = I_ERI_I4x2y_S_G4y_S+ABX*I_ERI_H3x2y_S_G4y_S;
  Double I_ERI_H3xyz_Px_G4y_S = I_ERI_I4xyz_S_G4y_S+ABX*I_ERI_H3xyz_S_G4y_S;
  Double I_ERI_H3x2z_Px_G4y_S = I_ERI_I4x2z_S_G4y_S+ABX*I_ERI_H3x2z_S_G4y_S;
  Double I_ERI_H2x3y_Px_G4y_S = I_ERI_I3x3y_S_G4y_S+ABX*I_ERI_H2x3y_S_G4y_S;
  Double I_ERI_H2x2yz_Px_G4y_S = I_ERI_I3x2yz_S_G4y_S+ABX*I_ERI_H2x2yz_S_G4y_S;
  Double I_ERI_H2xy2z_Px_G4y_S = I_ERI_I3xy2z_S_G4y_S+ABX*I_ERI_H2xy2z_S_G4y_S;
  Double I_ERI_H2x3z_Px_G4y_S = I_ERI_I3x3z_S_G4y_S+ABX*I_ERI_H2x3z_S_G4y_S;
  Double I_ERI_Hx4y_Px_G4y_S = I_ERI_I2x4y_S_G4y_S+ABX*I_ERI_Hx4y_S_G4y_S;
  Double I_ERI_Hx3yz_Px_G4y_S = I_ERI_I2x3yz_S_G4y_S+ABX*I_ERI_Hx3yz_S_G4y_S;
  Double I_ERI_Hx2y2z_Px_G4y_S = I_ERI_I2x2y2z_S_G4y_S+ABX*I_ERI_Hx2y2z_S_G4y_S;
  Double I_ERI_Hxy3z_Px_G4y_S = I_ERI_I2xy3z_S_G4y_S+ABX*I_ERI_Hxy3z_S_G4y_S;
  Double I_ERI_Hx4z_Px_G4y_S = I_ERI_I2x4z_S_G4y_S+ABX*I_ERI_Hx4z_S_G4y_S;
  Double I_ERI_H5y_Px_G4y_S = I_ERI_Ix5y_S_G4y_S+ABX*I_ERI_H5y_S_G4y_S;
  Double I_ERI_H4yz_Px_G4y_S = I_ERI_Ix4yz_S_G4y_S+ABX*I_ERI_H4yz_S_G4y_S;
  Double I_ERI_H3y2z_Px_G4y_S = I_ERI_Ix3y2z_S_G4y_S+ABX*I_ERI_H3y2z_S_G4y_S;
  Double I_ERI_H2y3z_Px_G4y_S = I_ERI_Ix2y3z_S_G4y_S+ABX*I_ERI_H2y3z_S_G4y_S;
  Double I_ERI_Hy4z_Px_G4y_S = I_ERI_Ixy4z_S_G4y_S+ABX*I_ERI_Hy4z_S_G4y_S;
  Double I_ERI_H5z_Px_G4y_S = I_ERI_Ix5z_S_G4y_S+ABX*I_ERI_H5z_S_G4y_S;
  Double I_ERI_H4xy_Py_G4y_S = I_ERI_I4x2y_S_G4y_S+ABY*I_ERI_H4xy_S_G4y_S;
  Double I_ERI_H4xz_Py_G4y_S = I_ERI_I4xyz_S_G4y_S+ABY*I_ERI_H4xz_S_G4y_S;
  Double I_ERI_H3x2y_Py_G4y_S = I_ERI_I3x3y_S_G4y_S+ABY*I_ERI_H3x2y_S_G4y_S;
  Double I_ERI_H3xyz_Py_G4y_S = I_ERI_I3x2yz_S_G4y_S+ABY*I_ERI_H3xyz_S_G4y_S;
  Double I_ERI_H3x2z_Py_G4y_S = I_ERI_I3xy2z_S_G4y_S+ABY*I_ERI_H3x2z_S_G4y_S;
  Double I_ERI_H2x3y_Py_G4y_S = I_ERI_I2x4y_S_G4y_S+ABY*I_ERI_H2x3y_S_G4y_S;
  Double I_ERI_H2x2yz_Py_G4y_S = I_ERI_I2x3yz_S_G4y_S+ABY*I_ERI_H2x2yz_S_G4y_S;
  Double I_ERI_H2xy2z_Py_G4y_S = I_ERI_I2x2y2z_S_G4y_S+ABY*I_ERI_H2xy2z_S_G4y_S;
  Double I_ERI_H2x3z_Py_G4y_S = I_ERI_I2xy3z_S_G4y_S+ABY*I_ERI_H2x3z_S_G4y_S;
  Double I_ERI_Hx4y_Py_G4y_S = I_ERI_Ix5y_S_G4y_S+ABY*I_ERI_Hx4y_S_G4y_S;
  Double I_ERI_Hx3yz_Py_G4y_S = I_ERI_Ix4yz_S_G4y_S+ABY*I_ERI_Hx3yz_S_G4y_S;
  Double I_ERI_Hx2y2z_Py_G4y_S = I_ERI_Ix3y2z_S_G4y_S+ABY*I_ERI_Hx2y2z_S_G4y_S;
  Double I_ERI_Hxy3z_Py_G4y_S = I_ERI_Ix2y3z_S_G4y_S+ABY*I_ERI_Hxy3z_S_G4y_S;
  Double I_ERI_Hx4z_Py_G4y_S = I_ERI_Ixy4z_S_G4y_S+ABY*I_ERI_Hx4z_S_G4y_S;
  Double I_ERI_H5y_Py_G4y_S = I_ERI_I6y_S_G4y_S+ABY*I_ERI_H5y_S_G4y_S;
  Double I_ERI_H4yz_Py_G4y_S = I_ERI_I5yz_S_G4y_S+ABY*I_ERI_H4yz_S_G4y_S;
  Double I_ERI_H3y2z_Py_G4y_S = I_ERI_I4y2z_S_G4y_S+ABY*I_ERI_H3y2z_S_G4y_S;
  Double I_ERI_H2y3z_Py_G4y_S = I_ERI_I3y3z_S_G4y_S+ABY*I_ERI_H2y3z_S_G4y_S;
  Double I_ERI_Hy4z_Py_G4y_S = I_ERI_I2y4z_S_G4y_S+ABY*I_ERI_Hy4z_S_G4y_S;
  Double I_ERI_H5z_Py_G4y_S = I_ERI_Iy5z_S_G4y_S+ABY*I_ERI_H5z_S_G4y_S;
  Double I_ERI_H4xz_Pz_G4y_S = I_ERI_I4x2z_S_G4y_S+ABZ*I_ERI_H4xz_S_G4y_S;
  Double I_ERI_H3xyz_Pz_G4y_S = I_ERI_I3xy2z_S_G4y_S+ABZ*I_ERI_H3xyz_S_G4y_S;
  Double I_ERI_H3x2z_Pz_G4y_S = I_ERI_I3x3z_S_G4y_S+ABZ*I_ERI_H3x2z_S_G4y_S;
  Double I_ERI_H2x2yz_Pz_G4y_S = I_ERI_I2x2y2z_S_G4y_S+ABZ*I_ERI_H2x2yz_S_G4y_S;
  Double I_ERI_H2xy2z_Pz_G4y_S = I_ERI_I2xy3z_S_G4y_S+ABZ*I_ERI_H2xy2z_S_G4y_S;
  Double I_ERI_H2x3z_Pz_G4y_S = I_ERI_I2x4z_S_G4y_S+ABZ*I_ERI_H2x3z_S_G4y_S;
  Double I_ERI_Hx3yz_Pz_G4y_S = I_ERI_Ix3y2z_S_G4y_S+ABZ*I_ERI_Hx3yz_S_G4y_S;
  Double I_ERI_Hx2y2z_Pz_G4y_S = I_ERI_Ix2y3z_S_G4y_S+ABZ*I_ERI_Hx2y2z_S_G4y_S;
  Double I_ERI_Hxy3z_Pz_G4y_S = I_ERI_Ixy4z_S_G4y_S+ABZ*I_ERI_Hxy3z_S_G4y_S;
  Double I_ERI_Hx4z_Pz_G4y_S = I_ERI_Ix5z_S_G4y_S+ABZ*I_ERI_Hx4z_S_G4y_S;
  Double I_ERI_H4yz_Pz_G4y_S = I_ERI_I4y2z_S_G4y_S+ABZ*I_ERI_H4yz_S_G4y_S;
  Double I_ERI_H3y2z_Pz_G4y_S = I_ERI_I3y3z_S_G4y_S+ABZ*I_ERI_H3y2z_S_G4y_S;
  Double I_ERI_H2y3z_Pz_G4y_S = I_ERI_I2y4z_S_G4y_S+ABZ*I_ERI_H2y3z_S_G4y_S;
  Double I_ERI_Hy4z_Pz_G4y_S = I_ERI_Iy5z_S_G4y_S+ABZ*I_ERI_Hy4z_S_G4y_S;
  Double I_ERI_H5z_Pz_G4y_S = I_ERI_I6z_S_G4y_S+ABZ*I_ERI_H5z_S_G4y_S;
  Double I_ERI_H5x_Px_G3yz_S = I_ERI_I6x_S_G3yz_S+ABX*I_ERI_H5x_S_G3yz_S;
  Double I_ERI_H4xy_Px_G3yz_S = I_ERI_I5xy_S_G3yz_S+ABX*I_ERI_H4xy_S_G3yz_S;
  Double I_ERI_H4xz_Px_G3yz_S = I_ERI_I5xz_S_G3yz_S+ABX*I_ERI_H4xz_S_G3yz_S;
  Double I_ERI_H3x2y_Px_G3yz_S = I_ERI_I4x2y_S_G3yz_S+ABX*I_ERI_H3x2y_S_G3yz_S;
  Double I_ERI_H3xyz_Px_G3yz_S = I_ERI_I4xyz_S_G3yz_S+ABX*I_ERI_H3xyz_S_G3yz_S;
  Double I_ERI_H3x2z_Px_G3yz_S = I_ERI_I4x2z_S_G3yz_S+ABX*I_ERI_H3x2z_S_G3yz_S;
  Double I_ERI_H2x3y_Px_G3yz_S = I_ERI_I3x3y_S_G3yz_S+ABX*I_ERI_H2x3y_S_G3yz_S;
  Double I_ERI_H2x2yz_Px_G3yz_S = I_ERI_I3x2yz_S_G3yz_S+ABX*I_ERI_H2x2yz_S_G3yz_S;
  Double I_ERI_H2xy2z_Px_G3yz_S = I_ERI_I3xy2z_S_G3yz_S+ABX*I_ERI_H2xy2z_S_G3yz_S;
  Double I_ERI_H2x3z_Px_G3yz_S = I_ERI_I3x3z_S_G3yz_S+ABX*I_ERI_H2x3z_S_G3yz_S;
  Double I_ERI_Hx4y_Px_G3yz_S = I_ERI_I2x4y_S_G3yz_S+ABX*I_ERI_Hx4y_S_G3yz_S;
  Double I_ERI_Hx3yz_Px_G3yz_S = I_ERI_I2x3yz_S_G3yz_S+ABX*I_ERI_Hx3yz_S_G3yz_S;
  Double I_ERI_Hx2y2z_Px_G3yz_S = I_ERI_I2x2y2z_S_G3yz_S+ABX*I_ERI_Hx2y2z_S_G3yz_S;
  Double I_ERI_Hxy3z_Px_G3yz_S = I_ERI_I2xy3z_S_G3yz_S+ABX*I_ERI_Hxy3z_S_G3yz_S;
  Double I_ERI_Hx4z_Px_G3yz_S = I_ERI_I2x4z_S_G3yz_S+ABX*I_ERI_Hx4z_S_G3yz_S;
  Double I_ERI_H5y_Px_G3yz_S = I_ERI_Ix5y_S_G3yz_S+ABX*I_ERI_H5y_S_G3yz_S;
  Double I_ERI_H4yz_Px_G3yz_S = I_ERI_Ix4yz_S_G3yz_S+ABX*I_ERI_H4yz_S_G3yz_S;
  Double I_ERI_H3y2z_Px_G3yz_S = I_ERI_Ix3y2z_S_G3yz_S+ABX*I_ERI_H3y2z_S_G3yz_S;
  Double I_ERI_H2y3z_Px_G3yz_S = I_ERI_Ix2y3z_S_G3yz_S+ABX*I_ERI_H2y3z_S_G3yz_S;
  Double I_ERI_Hy4z_Px_G3yz_S = I_ERI_Ixy4z_S_G3yz_S+ABX*I_ERI_Hy4z_S_G3yz_S;
  Double I_ERI_H5z_Px_G3yz_S = I_ERI_Ix5z_S_G3yz_S+ABX*I_ERI_H5z_S_G3yz_S;
  Double I_ERI_H4xy_Py_G3yz_S = I_ERI_I4x2y_S_G3yz_S+ABY*I_ERI_H4xy_S_G3yz_S;
  Double I_ERI_H4xz_Py_G3yz_S = I_ERI_I4xyz_S_G3yz_S+ABY*I_ERI_H4xz_S_G3yz_S;
  Double I_ERI_H3x2y_Py_G3yz_S = I_ERI_I3x3y_S_G3yz_S+ABY*I_ERI_H3x2y_S_G3yz_S;
  Double I_ERI_H3xyz_Py_G3yz_S = I_ERI_I3x2yz_S_G3yz_S+ABY*I_ERI_H3xyz_S_G3yz_S;
  Double I_ERI_H3x2z_Py_G3yz_S = I_ERI_I3xy2z_S_G3yz_S+ABY*I_ERI_H3x2z_S_G3yz_S;
  Double I_ERI_H2x3y_Py_G3yz_S = I_ERI_I2x4y_S_G3yz_S+ABY*I_ERI_H2x3y_S_G3yz_S;
  Double I_ERI_H2x2yz_Py_G3yz_S = I_ERI_I2x3yz_S_G3yz_S+ABY*I_ERI_H2x2yz_S_G3yz_S;
  Double I_ERI_H2xy2z_Py_G3yz_S = I_ERI_I2x2y2z_S_G3yz_S+ABY*I_ERI_H2xy2z_S_G3yz_S;
  Double I_ERI_H2x3z_Py_G3yz_S = I_ERI_I2xy3z_S_G3yz_S+ABY*I_ERI_H2x3z_S_G3yz_S;
  Double I_ERI_Hx4y_Py_G3yz_S = I_ERI_Ix5y_S_G3yz_S+ABY*I_ERI_Hx4y_S_G3yz_S;
  Double I_ERI_Hx3yz_Py_G3yz_S = I_ERI_Ix4yz_S_G3yz_S+ABY*I_ERI_Hx3yz_S_G3yz_S;
  Double I_ERI_Hx2y2z_Py_G3yz_S = I_ERI_Ix3y2z_S_G3yz_S+ABY*I_ERI_Hx2y2z_S_G3yz_S;
  Double I_ERI_Hxy3z_Py_G3yz_S = I_ERI_Ix2y3z_S_G3yz_S+ABY*I_ERI_Hxy3z_S_G3yz_S;
  Double I_ERI_Hx4z_Py_G3yz_S = I_ERI_Ixy4z_S_G3yz_S+ABY*I_ERI_Hx4z_S_G3yz_S;
  Double I_ERI_H5y_Py_G3yz_S = I_ERI_I6y_S_G3yz_S+ABY*I_ERI_H5y_S_G3yz_S;
  Double I_ERI_H4yz_Py_G3yz_S = I_ERI_I5yz_S_G3yz_S+ABY*I_ERI_H4yz_S_G3yz_S;
  Double I_ERI_H3y2z_Py_G3yz_S = I_ERI_I4y2z_S_G3yz_S+ABY*I_ERI_H3y2z_S_G3yz_S;
  Double I_ERI_H2y3z_Py_G3yz_S = I_ERI_I3y3z_S_G3yz_S+ABY*I_ERI_H2y3z_S_G3yz_S;
  Double I_ERI_Hy4z_Py_G3yz_S = I_ERI_I2y4z_S_G3yz_S+ABY*I_ERI_Hy4z_S_G3yz_S;
  Double I_ERI_H5z_Py_G3yz_S = I_ERI_Iy5z_S_G3yz_S+ABY*I_ERI_H5z_S_G3yz_S;
  Double I_ERI_H4xz_Pz_G3yz_S = I_ERI_I4x2z_S_G3yz_S+ABZ*I_ERI_H4xz_S_G3yz_S;
  Double I_ERI_H3xyz_Pz_G3yz_S = I_ERI_I3xy2z_S_G3yz_S+ABZ*I_ERI_H3xyz_S_G3yz_S;
  Double I_ERI_H3x2z_Pz_G3yz_S = I_ERI_I3x3z_S_G3yz_S+ABZ*I_ERI_H3x2z_S_G3yz_S;
  Double I_ERI_H2x2yz_Pz_G3yz_S = I_ERI_I2x2y2z_S_G3yz_S+ABZ*I_ERI_H2x2yz_S_G3yz_S;
  Double I_ERI_H2xy2z_Pz_G3yz_S = I_ERI_I2xy3z_S_G3yz_S+ABZ*I_ERI_H2xy2z_S_G3yz_S;
  Double I_ERI_H2x3z_Pz_G3yz_S = I_ERI_I2x4z_S_G3yz_S+ABZ*I_ERI_H2x3z_S_G3yz_S;
  Double I_ERI_Hx3yz_Pz_G3yz_S = I_ERI_Ix3y2z_S_G3yz_S+ABZ*I_ERI_Hx3yz_S_G3yz_S;
  Double I_ERI_Hx2y2z_Pz_G3yz_S = I_ERI_Ix2y3z_S_G3yz_S+ABZ*I_ERI_Hx2y2z_S_G3yz_S;
  Double I_ERI_Hxy3z_Pz_G3yz_S = I_ERI_Ixy4z_S_G3yz_S+ABZ*I_ERI_Hxy3z_S_G3yz_S;
  Double I_ERI_Hx4z_Pz_G3yz_S = I_ERI_Ix5z_S_G3yz_S+ABZ*I_ERI_Hx4z_S_G3yz_S;
  Double I_ERI_H4yz_Pz_G3yz_S = I_ERI_I4y2z_S_G3yz_S+ABZ*I_ERI_H4yz_S_G3yz_S;
  Double I_ERI_H3y2z_Pz_G3yz_S = I_ERI_I3y3z_S_G3yz_S+ABZ*I_ERI_H3y2z_S_G3yz_S;
  Double I_ERI_H2y3z_Pz_G3yz_S = I_ERI_I2y4z_S_G3yz_S+ABZ*I_ERI_H2y3z_S_G3yz_S;
  Double I_ERI_Hy4z_Pz_G3yz_S = I_ERI_Iy5z_S_G3yz_S+ABZ*I_ERI_Hy4z_S_G3yz_S;
  Double I_ERI_H5z_Pz_G3yz_S = I_ERI_I6z_S_G3yz_S+ABZ*I_ERI_H5z_S_G3yz_S;
  Double I_ERI_H5x_Px_G2y2z_S = I_ERI_I6x_S_G2y2z_S+ABX*I_ERI_H5x_S_G2y2z_S;
  Double I_ERI_H4xy_Px_G2y2z_S = I_ERI_I5xy_S_G2y2z_S+ABX*I_ERI_H4xy_S_G2y2z_S;
  Double I_ERI_H4xz_Px_G2y2z_S = I_ERI_I5xz_S_G2y2z_S+ABX*I_ERI_H4xz_S_G2y2z_S;
  Double I_ERI_H3x2y_Px_G2y2z_S = I_ERI_I4x2y_S_G2y2z_S+ABX*I_ERI_H3x2y_S_G2y2z_S;
  Double I_ERI_H3xyz_Px_G2y2z_S = I_ERI_I4xyz_S_G2y2z_S+ABX*I_ERI_H3xyz_S_G2y2z_S;
  Double I_ERI_H3x2z_Px_G2y2z_S = I_ERI_I4x2z_S_G2y2z_S+ABX*I_ERI_H3x2z_S_G2y2z_S;
  Double I_ERI_H2x3y_Px_G2y2z_S = I_ERI_I3x3y_S_G2y2z_S+ABX*I_ERI_H2x3y_S_G2y2z_S;
  Double I_ERI_H2x2yz_Px_G2y2z_S = I_ERI_I3x2yz_S_G2y2z_S+ABX*I_ERI_H2x2yz_S_G2y2z_S;
  Double I_ERI_H2xy2z_Px_G2y2z_S = I_ERI_I3xy2z_S_G2y2z_S+ABX*I_ERI_H2xy2z_S_G2y2z_S;
  Double I_ERI_H2x3z_Px_G2y2z_S = I_ERI_I3x3z_S_G2y2z_S+ABX*I_ERI_H2x3z_S_G2y2z_S;
  Double I_ERI_Hx4y_Px_G2y2z_S = I_ERI_I2x4y_S_G2y2z_S+ABX*I_ERI_Hx4y_S_G2y2z_S;
  Double I_ERI_Hx3yz_Px_G2y2z_S = I_ERI_I2x3yz_S_G2y2z_S+ABX*I_ERI_Hx3yz_S_G2y2z_S;
  Double I_ERI_Hx2y2z_Px_G2y2z_S = I_ERI_I2x2y2z_S_G2y2z_S+ABX*I_ERI_Hx2y2z_S_G2y2z_S;
  Double I_ERI_Hxy3z_Px_G2y2z_S = I_ERI_I2xy3z_S_G2y2z_S+ABX*I_ERI_Hxy3z_S_G2y2z_S;
  Double I_ERI_Hx4z_Px_G2y2z_S = I_ERI_I2x4z_S_G2y2z_S+ABX*I_ERI_Hx4z_S_G2y2z_S;
  Double I_ERI_H5y_Px_G2y2z_S = I_ERI_Ix5y_S_G2y2z_S+ABX*I_ERI_H5y_S_G2y2z_S;
  Double I_ERI_H4yz_Px_G2y2z_S = I_ERI_Ix4yz_S_G2y2z_S+ABX*I_ERI_H4yz_S_G2y2z_S;
  Double I_ERI_H3y2z_Px_G2y2z_S = I_ERI_Ix3y2z_S_G2y2z_S+ABX*I_ERI_H3y2z_S_G2y2z_S;
  Double I_ERI_H2y3z_Px_G2y2z_S = I_ERI_Ix2y3z_S_G2y2z_S+ABX*I_ERI_H2y3z_S_G2y2z_S;
  Double I_ERI_Hy4z_Px_G2y2z_S = I_ERI_Ixy4z_S_G2y2z_S+ABX*I_ERI_Hy4z_S_G2y2z_S;
  Double I_ERI_H5z_Px_G2y2z_S = I_ERI_Ix5z_S_G2y2z_S+ABX*I_ERI_H5z_S_G2y2z_S;
  Double I_ERI_H4xy_Py_G2y2z_S = I_ERI_I4x2y_S_G2y2z_S+ABY*I_ERI_H4xy_S_G2y2z_S;
  Double I_ERI_H4xz_Py_G2y2z_S = I_ERI_I4xyz_S_G2y2z_S+ABY*I_ERI_H4xz_S_G2y2z_S;
  Double I_ERI_H3x2y_Py_G2y2z_S = I_ERI_I3x3y_S_G2y2z_S+ABY*I_ERI_H3x2y_S_G2y2z_S;
  Double I_ERI_H3xyz_Py_G2y2z_S = I_ERI_I3x2yz_S_G2y2z_S+ABY*I_ERI_H3xyz_S_G2y2z_S;
  Double I_ERI_H3x2z_Py_G2y2z_S = I_ERI_I3xy2z_S_G2y2z_S+ABY*I_ERI_H3x2z_S_G2y2z_S;
  Double I_ERI_H2x3y_Py_G2y2z_S = I_ERI_I2x4y_S_G2y2z_S+ABY*I_ERI_H2x3y_S_G2y2z_S;
  Double I_ERI_H2x2yz_Py_G2y2z_S = I_ERI_I2x3yz_S_G2y2z_S+ABY*I_ERI_H2x2yz_S_G2y2z_S;
  Double I_ERI_H2xy2z_Py_G2y2z_S = I_ERI_I2x2y2z_S_G2y2z_S+ABY*I_ERI_H2xy2z_S_G2y2z_S;
  Double I_ERI_H2x3z_Py_G2y2z_S = I_ERI_I2xy3z_S_G2y2z_S+ABY*I_ERI_H2x3z_S_G2y2z_S;
  Double I_ERI_Hx4y_Py_G2y2z_S = I_ERI_Ix5y_S_G2y2z_S+ABY*I_ERI_Hx4y_S_G2y2z_S;
  Double I_ERI_Hx3yz_Py_G2y2z_S = I_ERI_Ix4yz_S_G2y2z_S+ABY*I_ERI_Hx3yz_S_G2y2z_S;
  Double I_ERI_Hx2y2z_Py_G2y2z_S = I_ERI_Ix3y2z_S_G2y2z_S+ABY*I_ERI_Hx2y2z_S_G2y2z_S;
  Double I_ERI_Hxy3z_Py_G2y2z_S = I_ERI_Ix2y3z_S_G2y2z_S+ABY*I_ERI_Hxy3z_S_G2y2z_S;
  Double I_ERI_Hx4z_Py_G2y2z_S = I_ERI_Ixy4z_S_G2y2z_S+ABY*I_ERI_Hx4z_S_G2y2z_S;
  Double I_ERI_H5y_Py_G2y2z_S = I_ERI_I6y_S_G2y2z_S+ABY*I_ERI_H5y_S_G2y2z_S;
  Double I_ERI_H4yz_Py_G2y2z_S = I_ERI_I5yz_S_G2y2z_S+ABY*I_ERI_H4yz_S_G2y2z_S;
  Double I_ERI_H3y2z_Py_G2y2z_S = I_ERI_I4y2z_S_G2y2z_S+ABY*I_ERI_H3y2z_S_G2y2z_S;
  Double I_ERI_H2y3z_Py_G2y2z_S = I_ERI_I3y3z_S_G2y2z_S+ABY*I_ERI_H2y3z_S_G2y2z_S;
  Double I_ERI_Hy4z_Py_G2y2z_S = I_ERI_I2y4z_S_G2y2z_S+ABY*I_ERI_Hy4z_S_G2y2z_S;
  Double I_ERI_H5z_Py_G2y2z_S = I_ERI_Iy5z_S_G2y2z_S+ABY*I_ERI_H5z_S_G2y2z_S;
  Double I_ERI_H4xz_Pz_G2y2z_S = I_ERI_I4x2z_S_G2y2z_S+ABZ*I_ERI_H4xz_S_G2y2z_S;
  Double I_ERI_H3xyz_Pz_G2y2z_S = I_ERI_I3xy2z_S_G2y2z_S+ABZ*I_ERI_H3xyz_S_G2y2z_S;
  Double I_ERI_H3x2z_Pz_G2y2z_S = I_ERI_I3x3z_S_G2y2z_S+ABZ*I_ERI_H3x2z_S_G2y2z_S;
  Double I_ERI_H2x2yz_Pz_G2y2z_S = I_ERI_I2x2y2z_S_G2y2z_S+ABZ*I_ERI_H2x2yz_S_G2y2z_S;
  Double I_ERI_H2xy2z_Pz_G2y2z_S = I_ERI_I2xy3z_S_G2y2z_S+ABZ*I_ERI_H2xy2z_S_G2y2z_S;
  Double I_ERI_H2x3z_Pz_G2y2z_S = I_ERI_I2x4z_S_G2y2z_S+ABZ*I_ERI_H2x3z_S_G2y2z_S;
  Double I_ERI_Hx3yz_Pz_G2y2z_S = I_ERI_Ix3y2z_S_G2y2z_S+ABZ*I_ERI_Hx3yz_S_G2y2z_S;
  Double I_ERI_Hx2y2z_Pz_G2y2z_S = I_ERI_Ix2y3z_S_G2y2z_S+ABZ*I_ERI_Hx2y2z_S_G2y2z_S;
  Double I_ERI_Hxy3z_Pz_G2y2z_S = I_ERI_Ixy4z_S_G2y2z_S+ABZ*I_ERI_Hxy3z_S_G2y2z_S;
  Double I_ERI_Hx4z_Pz_G2y2z_S = I_ERI_Ix5z_S_G2y2z_S+ABZ*I_ERI_Hx4z_S_G2y2z_S;
  Double I_ERI_H4yz_Pz_G2y2z_S = I_ERI_I4y2z_S_G2y2z_S+ABZ*I_ERI_H4yz_S_G2y2z_S;
  Double I_ERI_H3y2z_Pz_G2y2z_S = I_ERI_I3y3z_S_G2y2z_S+ABZ*I_ERI_H3y2z_S_G2y2z_S;
  Double I_ERI_H2y3z_Pz_G2y2z_S = I_ERI_I2y4z_S_G2y2z_S+ABZ*I_ERI_H2y3z_S_G2y2z_S;
  Double I_ERI_Hy4z_Pz_G2y2z_S = I_ERI_Iy5z_S_G2y2z_S+ABZ*I_ERI_Hy4z_S_G2y2z_S;
  Double I_ERI_H5z_Pz_G2y2z_S = I_ERI_I6z_S_G2y2z_S+ABZ*I_ERI_H5z_S_G2y2z_S;
  Double I_ERI_H5x_Px_Gy3z_S = I_ERI_I6x_S_Gy3z_S+ABX*I_ERI_H5x_S_Gy3z_S;
  Double I_ERI_H4xy_Px_Gy3z_S = I_ERI_I5xy_S_Gy3z_S+ABX*I_ERI_H4xy_S_Gy3z_S;
  Double I_ERI_H4xz_Px_Gy3z_S = I_ERI_I5xz_S_Gy3z_S+ABX*I_ERI_H4xz_S_Gy3z_S;
  Double I_ERI_H3x2y_Px_Gy3z_S = I_ERI_I4x2y_S_Gy3z_S+ABX*I_ERI_H3x2y_S_Gy3z_S;
  Double I_ERI_H3xyz_Px_Gy3z_S = I_ERI_I4xyz_S_Gy3z_S+ABX*I_ERI_H3xyz_S_Gy3z_S;
  Double I_ERI_H3x2z_Px_Gy3z_S = I_ERI_I4x2z_S_Gy3z_S+ABX*I_ERI_H3x2z_S_Gy3z_S;
  Double I_ERI_H2x3y_Px_Gy3z_S = I_ERI_I3x3y_S_Gy3z_S+ABX*I_ERI_H2x3y_S_Gy3z_S;
  Double I_ERI_H2x2yz_Px_Gy3z_S = I_ERI_I3x2yz_S_Gy3z_S+ABX*I_ERI_H2x2yz_S_Gy3z_S;
  Double I_ERI_H2xy2z_Px_Gy3z_S = I_ERI_I3xy2z_S_Gy3z_S+ABX*I_ERI_H2xy2z_S_Gy3z_S;
  Double I_ERI_H2x3z_Px_Gy3z_S = I_ERI_I3x3z_S_Gy3z_S+ABX*I_ERI_H2x3z_S_Gy3z_S;
  Double I_ERI_Hx4y_Px_Gy3z_S = I_ERI_I2x4y_S_Gy3z_S+ABX*I_ERI_Hx4y_S_Gy3z_S;
  Double I_ERI_Hx3yz_Px_Gy3z_S = I_ERI_I2x3yz_S_Gy3z_S+ABX*I_ERI_Hx3yz_S_Gy3z_S;
  Double I_ERI_Hx2y2z_Px_Gy3z_S = I_ERI_I2x2y2z_S_Gy3z_S+ABX*I_ERI_Hx2y2z_S_Gy3z_S;
  Double I_ERI_Hxy3z_Px_Gy3z_S = I_ERI_I2xy3z_S_Gy3z_S+ABX*I_ERI_Hxy3z_S_Gy3z_S;
  Double I_ERI_Hx4z_Px_Gy3z_S = I_ERI_I2x4z_S_Gy3z_S+ABX*I_ERI_Hx4z_S_Gy3z_S;
  Double I_ERI_H5y_Px_Gy3z_S = I_ERI_Ix5y_S_Gy3z_S+ABX*I_ERI_H5y_S_Gy3z_S;
  Double I_ERI_H4yz_Px_Gy3z_S = I_ERI_Ix4yz_S_Gy3z_S+ABX*I_ERI_H4yz_S_Gy3z_S;
  Double I_ERI_H3y2z_Px_Gy3z_S = I_ERI_Ix3y2z_S_Gy3z_S+ABX*I_ERI_H3y2z_S_Gy3z_S;
  Double I_ERI_H2y3z_Px_Gy3z_S = I_ERI_Ix2y3z_S_Gy3z_S+ABX*I_ERI_H2y3z_S_Gy3z_S;
  Double I_ERI_Hy4z_Px_Gy3z_S = I_ERI_Ixy4z_S_Gy3z_S+ABX*I_ERI_Hy4z_S_Gy3z_S;
  Double I_ERI_H5z_Px_Gy3z_S = I_ERI_Ix5z_S_Gy3z_S+ABX*I_ERI_H5z_S_Gy3z_S;
  Double I_ERI_H4xy_Py_Gy3z_S = I_ERI_I4x2y_S_Gy3z_S+ABY*I_ERI_H4xy_S_Gy3z_S;
  Double I_ERI_H4xz_Py_Gy3z_S = I_ERI_I4xyz_S_Gy3z_S+ABY*I_ERI_H4xz_S_Gy3z_S;
  Double I_ERI_H3x2y_Py_Gy3z_S = I_ERI_I3x3y_S_Gy3z_S+ABY*I_ERI_H3x2y_S_Gy3z_S;
  Double I_ERI_H3xyz_Py_Gy3z_S = I_ERI_I3x2yz_S_Gy3z_S+ABY*I_ERI_H3xyz_S_Gy3z_S;
  Double I_ERI_H3x2z_Py_Gy3z_S = I_ERI_I3xy2z_S_Gy3z_S+ABY*I_ERI_H3x2z_S_Gy3z_S;
  Double I_ERI_H2x3y_Py_Gy3z_S = I_ERI_I2x4y_S_Gy3z_S+ABY*I_ERI_H2x3y_S_Gy3z_S;
  Double I_ERI_H2x2yz_Py_Gy3z_S = I_ERI_I2x3yz_S_Gy3z_S+ABY*I_ERI_H2x2yz_S_Gy3z_S;
  Double I_ERI_H2xy2z_Py_Gy3z_S = I_ERI_I2x2y2z_S_Gy3z_S+ABY*I_ERI_H2xy2z_S_Gy3z_S;
  Double I_ERI_H2x3z_Py_Gy3z_S = I_ERI_I2xy3z_S_Gy3z_S+ABY*I_ERI_H2x3z_S_Gy3z_S;
  Double I_ERI_Hx4y_Py_Gy3z_S = I_ERI_Ix5y_S_Gy3z_S+ABY*I_ERI_Hx4y_S_Gy3z_S;
  Double I_ERI_Hx3yz_Py_Gy3z_S = I_ERI_Ix4yz_S_Gy3z_S+ABY*I_ERI_Hx3yz_S_Gy3z_S;
  Double I_ERI_Hx2y2z_Py_Gy3z_S = I_ERI_Ix3y2z_S_Gy3z_S+ABY*I_ERI_Hx2y2z_S_Gy3z_S;
  Double I_ERI_Hxy3z_Py_Gy3z_S = I_ERI_Ix2y3z_S_Gy3z_S+ABY*I_ERI_Hxy3z_S_Gy3z_S;
  Double I_ERI_Hx4z_Py_Gy3z_S = I_ERI_Ixy4z_S_Gy3z_S+ABY*I_ERI_Hx4z_S_Gy3z_S;
  Double I_ERI_H5y_Py_Gy3z_S = I_ERI_I6y_S_Gy3z_S+ABY*I_ERI_H5y_S_Gy3z_S;
  Double I_ERI_H4yz_Py_Gy3z_S = I_ERI_I5yz_S_Gy3z_S+ABY*I_ERI_H4yz_S_Gy3z_S;
  Double I_ERI_H3y2z_Py_Gy3z_S = I_ERI_I4y2z_S_Gy3z_S+ABY*I_ERI_H3y2z_S_Gy3z_S;
  Double I_ERI_H2y3z_Py_Gy3z_S = I_ERI_I3y3z_S_Gy3z_S+ABY*I_ERI_H2y3z_S_Gy3z_S;
  Double I_ERI_Hy4z_Py_Gy3z_S = I_ERI_I2y4z_S_Gy3z_S+ABY*I_ERI_Hy4z_S_Gy3z_S;
  Double I_ERI_H5z_Py_Gy3z_S = I_ERI_Iy5z_S_Gy3z_S+ABY*I_ERI_H5z_S_Gy3z_S;
  Double I_ERI_H4xz_Pz_Gy3z_S = I_ERI_I4x2z_S_Gy3z_S+ABZ*I_ERI_H4xz_S_Gy3z_S;
  Double I_ERI_H3xyz_Pz_Gy3z_S = I_ERI_I3xy2z_S_Gy3z_S+ABZ*I_ERI_H3xyz_S_Gy3z_S;
  Double I_ERI_H3x2z_Pz_Gy3z_S = I_ERI_I3x3z_S_Gy3z_S+ABZ*I_ERI_H3x2z_S_Gy3z_S;
  Double I_ERI_H2x2yz_Pz_Gy3z_S = I_ERI_I2x2y2z_S_Gy3z_S+ABZ*I_ERI_H2x2yz_S_Gy3z_S;
  Double I_ERI_H2xy2z_Pz_Gy3z_S = I_ERI_I2xy3z_S_Gy3z_S+ABZ*I_ERI_H2xy2z_S_Gy3z_S;
  Double I_ERI_H2x3z_Pz_Gy3z_S = I_ERI_I2x4z_S_Gy3z_S+ABZ*I_ERI_H2x3z_S_Gy3z_S;
  Double I_ERI_Hx3yz_Pz_Gy3z_S = I_ERI_Ix3y2z_S_Gy3z_S+ABZ*I_ERI_Hx3yz_S_Gy3z_S;
  Double I_ERI_Hx2y2z_Pz_Gy3z_S = I_ERI_Ix2y3z_S_Gy3z_S+ABZ*I_ERI_Hx2y2z_S_Gy3z_S;
  Double I_ERI_Hxy3z_Pz_Gy3z_S = I_ERI_Ixy4z_S_Gy3z_S+ABZ*I_ERI_Hxy3z_S_Gy3z_S;
  Double I_ERI_Hx4z_Pz_Gy3z_S = I_ERI_Ix5z_S_Gy3z_S+ABZ*I_ERI_Hx4z_S_Gy3z_S;
  Double I_ERI_H4yz_Pz_Gy3z_S = I_ERI_I4y2z_S_Gy3z_S+ABZ*I_ERI_H4yz_S_Gy3z_S;
  Double I_ERI_H3y2z_Pz_Gy3z_S = I_ERI_I3y3z_S_Gy3z_S+ABZ*I_ERI_H3y2z_S_Gy3z_S;
  Double I_ERI_H2y3z_Pz_Gy3z_S = I_ERI_I2y4z_S_Gy3z_S+ABZ*I_ERI_H2y3z_S_Gy3z_S;
  Double I_ERI_Hy4z_Pz_Gy3z_S = I_ERI_Iy5z_S_Gy3z_S+ABZ*I_ERI_Hy4z_S_Gy3z_S;
  Double I_ERI_H5z_Pz_Gy3z_S = I_ERI_I6z_S_Gy3z_S+ABZ*I_ERI_H5z_S_Gy3z_S;
  Double I_ERI_H5x_Px_G4z_S = I_ERI_I6x_S_G4z_S+ABX*I_ERI_H5x_S_G4z_S;
  Double I_ERI_H4xy_Px_G4z_S = I_ERI_I5xy_S_G4z_S+ABX*I_ERI_H4xy_S_G4z_S;
  Double I_ERI_H4xz_Px_G4z_S = I_ERI_I5xz_S_G4z_S+ABX*I_ERI_H4xz_S_G4z_S;
  Double I_ERI_H3x2y_Px_G4z_S = I_ERI_I4x2y_S_G4z_S+ABX*I_ERI_H3x2y_S_G4z_S;
  Double I_ERI_H3xyz_Px_G4z_S = I_ERI_I4xyz_S_G4z_S+ABX*I_ERI_H3xyz_S_G4z_S;
  Double I_ERI_H3x2z_Px_G4z_S = I_ERI_I4x2z_S_G4z_S+ABX*I_ERI_H3x2z_S_G4z_S;
  Double I_ERI_H2x3y_Px_G4z_S = I_ERI_I3x3y_S_G4z_S+ABX*I_ERI_H2x3y_S_G4z_S;
  Double I_ERI_H2x2yz_Px_G4z_S = I_ERI_I3x2yz_S_G4z_S+ABX*I_ERI_H2x2yz_S_G4z_S;
  Double I_ERI_H2xy2z_Px_G4z_S = I_ERI_I3xy2z_S_G4z_S+ABX*I_ERI_H2xy2z_S_G4z_S;
  Double I_ERI_H2x3z_Px_G4z_S = I_ERI_I3x3z_S_G4z_S+ABX*I_ERI_H2x3z_S_G4z_S;
  Double I_ERI_Hx4y_Px_G4z_S = I_ERI_I2x4y_S_G4z_S+ABX*I_ERI_Hx4y_S_G4z_S;
  Double I_ERI_Hx3yz_Px_G4z_S = I_ERI_I2x3yz_S_G4z_S+ABX*I_ERI_Hx3yz_S_G4z_S;
  Double I_ERI_Hx2y2z_Px_G4z_S = I_ERI_I2x2y2z_S_G4z_S+ABX*I_ERI_Hx2y2z_S_G4z_S;
  Double I_ERI_Hxy3z_Px_G4z_S = I_ERI_I2xy3z_S_G4z_S+ABX*I_ERI_Hxy3z_S_G4z_S;
  Double I_ERI_Hx4z_Px_G4z_S = I_ERI_I2x4z_S_G4z_S+ABX*I_ERI_Hx4z_S_G4z_S;
  Double I_ERI_H5y_Px_G4z_S = I_ERI_Ix5y_S_G4z_S+ABX*I_ERI_H5y_S_G4z_S;
  Double I_ERI_H4yz_Px_G4z_S = I_ERI_Ix4yz_S_G4z_S+ABX*I_ERI_H4yz_S_G4z_S;
  Double I_ERI_H3y2z_Px_G4z_S = I_ERI_Ix3y2z_S_G4z_S+ABX*I_ERI_H3y2z_S_G4z_S;
  Double I_ERI_H2y3z_Px_G4z_S = I_ERI_Ix2y3z_S_G4z_S+ABX*I_ERI_H2y3z_S_G4z_S;
  Double I_ERI_Hy4z_Px_G4z_S = I_ERI_Ixy4z_S_G4z_S+ABX*I_ERI_Hy4z_S_G4z_S;
  Double I_ERI_H5z_Px_G4z_S = I_ERI_Ix5z_S_G4z_S+ABX*I_ERI_H5z_S_G4z_S;
  Double I_ERI_H4xy_Py_G4z_S = I_ERI_I4x2y_S_G4z_S+ABY*I_ERI_H4xy_S_G4z_S;
  Double I_ERI_H4xz_Py_G4z_S = I_ERI_I4xyz_S_G4z_S+ABY*I_ERI_H4xz_S_G4z_S;
  Double I_ERI_H3x2y_Py_G4z_S = I_ERI_I3x3y_S_G4z_S+ABY*I_ERI_H3x2y_S_G4z_S;
  Double I_ERI_H3xyz_Py_G4z_S = I_ERI_I3x2yz_S_G4z_S+ABY*I_ERI_H3xyz_S_G4z_S;
  Double I_ERI_H3x2z_Py_G4z_S = I_ERI_I3xy2z_S_G4z_S+ABY*I_ERI_H3x2z_S_G4z_S;
  Double I_ERI_H2x3y_Py_G4z_S = I_ERI_I2x4y_S_G4z_S+ABY*I_ERI_H2x3y_S_G4z_S;
  Double I_ERI_H2x2yz_Py_G4z_S = I_ERI_I2x3yz_S_G4z_S+ABY*I_ERI_H2x2yz_S_G4z_S;
  Double I_ERI_H2xy2z_Py_G4z_S = I_ERI_I2x2y2z_S_G4z_S+ABY*I_ERI_H2xy2z_S_G4z_S;
  Double I_ERI_H2x3z_Py_G4z_S = I_ERI_I2xy3z_S_G4z_S+ABY*I_ERI_H2x3z_S_G4z_S;
  Double I_ERI_Hx4y_Py_G4z_S = I_ERI_Ix5y_S_G4z_S+ABY*I_ERI_Hx4y_S_G4z_S;
  Double I_ERI_Hx3yz_Py_G4z_S = I_ERI_Ix4yz_S_G4z_S+ABY*I_ERI_Hx3yz_S_G4z_S;
  Double I_ERI_Hx2y2z_Py_G4z_S = I_ERI_Ix3y2z_S_G4z_S+ABY*I_ERI_Hx2y2z_S_G4z_S;
  Double I_ERI_Hxy3z_Py_G4z_S = I_ERI_Ix2y3z_S_G4z_S+ABY*I_ERI_Hxy3z_S_G4z_S;
  Double I_ERI_Hx4z_Py_G4z_S = I_ERI_Ixy4z_S_G4z_S+ABY*I_ERI_Hx4z_S_G4z_S;
  Double I_ERI_H5y_Py_G4z_S = I_ERI_I6y_S_G4z_S+ABY*I_ERI_H5y_S_G4z_S;
  Double I_ERI_H4yz_Py_G4z_S = I_ERI_I5yz_S_G4z_S+ABY*I_ERI_H4yz_S_G4z_S;
  Double I_ERI_H3y2z_Py_G4z_S = I_ERI_I4y2z_S_G4z_S+ABY*I_ERI_H3y2z_S_G4z_S;
  Double I_ERI_H2y3z_Py_G4z_S = I_ERI_I3y3z_S_G4z_S+ABY*I_ERI_H2y3z_S_G4z_S;
  Double I_ERI_Hy4z_Py_G4z_S = I_ERI_I2y4z_S_G4z_S+ABY*I_ERI_Hy4z_S_G4z_S;
  Double I_ERI_H5z_Py_G4z_S = I_ERI_Iy5z_S_G4z_S+ABY*I_ERI_H5z_S_G4z_S;
  Double I_ERI_H4xz_Pz_G4z_S = I_ERI_I4x2z_S_G4z_S+ABZ*I_ERI_H4xz_S_G4z_S;
  Double I_ERI_H3xyz_Pz_G4z_S = I_ERI_I3xy2z_S_G4z_S+ABZ*I_ERI_H3xyz_S_G4z_S;
  Double I_ERI_H3x2z_Pz_G4z_S = I_ERI_I3x3z_S_G4z_S+ABZ*I_ERI_H3x2z_S_G4z_S;
  Double I_ERI_H2x2yz_Pz_G4z_S = I_ERI_I2x2y2z_S_G4z_S+ABZ*I_ERI_H2x2yz_S_G4z_S;
  Double I_ERI_H2xy2z_Pz_G4z_S = I_ERI_I2xy3z_S_G4z_S+ABZ*I_ERI_H2xy2z_S_G4z_S;
  Double I_ERI_H2x3z_Pz_G4z_S = I_ERI_I2x4z_S_G4z_S+ABZ*I_ERI_H2x3z_S_G4z_S;
  Double I_ERI_Hx3yz_Pz_G4z_S = I_ERI_Ix3y2z_S_G4z_S+ABZ*I_ERI_Hx3yz_S_G4z_S;
  Double I_ERI_Hx2y2z_Pz_G4z_S = I_ERI_Ix2y3z_S_G4z_S+ABZ*I_ERI_Hx2y2z_S_G4z_S;
  Double I_ERI_Hxy3z_Pz_G4z_S = I_ERI_Ixy4z_S_G4z_S+ABZ*I_ERI_Hxy3z_S_G4z_S;
  Double I_ERI_Hx4z_Pz_G4z_S = I_ERI_Ix5z_S_G4z_S+ABZ*I_ERI_Hx4z_S_G4z_S;
  Double I_ERI_H4yz_Pz_G4z_S = I_ERI_I4y2z_S_G4z_S+ABZ*I_ERI_H4yz_S_G4z_S;
  Double I_ERI_H3y2z_Pz_G4z_S = I_ERI_I3y3z_S_G4z_S+ABZ*I_ERI_H3y2z_S_G4z_S;
  Double I_ERI_H2y3z_Pz_G4z_S = I_ERI_I2y4z_S_G4z_S+ABZ*I_ERI_H2y3z_S_G4z_S;
  Double I_ERI_Hy4z_Pz_G4z_S = I_ERI_Iy5z_S_G4z_S+ABZ*I_ERI_Hy4z_S_G4z_S;
  Double I_ERI_H5z_Pz_G4z_S = I_ERI_I6z_S_G4z_S+ABZ*I_ERI_H5z_S_G4z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_G_D_G_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_P_G_S
   * RHS shell quartet name: SQ_ERI_G_P_G_S
   ************************************************************/
  abcd[0] = I_ERI_H5x_Px_G4x_S+ABX*I_ERI_G4x_Px_G4x_S;
  abcd[1] = I_ERI_H4xy_Px_G4x_S+ABX*I_ERI_G3xy_Px_G4x_S;
  abcd[2] = I_ERI_H4xz_Px_G4x_S+ABX*I_ERI_G3xz_Px_G4x_S;
  abcd[3] = I_ERI_H3x2y_Px_G4x_S+ABX*I_ERI_G2x2y_Px_G4x_S;
  abcd[4] = I_ERI_H3xyz_Px_G4x_S+ABX*I_ERI_G2xyz_Px_G4x_S;
  abcd[5] = I_ERI_H3x2z_Px_G4x_S+ABX*I_ERI_G2x2z_Px_G4x_S;
  abcd[6] = I_ERI_H2x3y_Px_G4x_S+ABX*I_ERI_Gx3y_Px_G4x_S;
  abcd[7] = I_ERI_H2x2yz_Px_G4x_S+ABX*I_ERI_Gx2yz_Px_G4x_S;
  abcd[8] = I_ERI_H2xy2z_Px_G4x_S+ABX*I_ERI_Gxy2z_Px_G4x_S;
  abcd[9] = I_ERI_H2x3z_Px_G4x_S+ABX*I_ERI_Gx3z_Px_G4x_S;
  abcd[10] = I_ERI_Hx4y_Px_G4x_S+ABX*I_ERI_G4y_Px_G4x_S;
  abcd[11] = I_ERI_Hx3yz_Px_G4x_S+ABX*I_ERI_G3yz_Px_G4x_S;
  abcd[12] = I_ERI_Hx2y2z_Px_G4x_S+ABX*I_ERI_G2y2z_Px_G4x_S;
  abcd[13] = I_ERI_Hxy3z_Px_G4x_S+ABX*I_ERI_Gy3z_Px_G4x_S;
  abcd[14] = I_ERI_Hx4z_Px_G4x_S+ABX*I_ERI_G4z_Px_G4x_S;
  abcd[15] = I_ERI_H4xy_Px_G4x_S+ABY*I_ERI_G4x_Px_G4x_S;
  abcd[16] = I_ERI_H3x2y_Px_G4x_S+ABY*I_ERI_G3xy_Px_G4x_S;
  abcd[17] = I_ERI_H3xyz_Px_G4x_S+ABY*I_ERI_G3xz_Px_G4x_S;
  abcd[18] = I_ERI_H2x3y_Px_G4x_S+ABY*I_ERI_G2x2y_Px_G4x_S;
  abcd[19] = I_ERI_H2x2yz_Px_G4x_S+ABY*I_ERI_G2xyz_Px_G4x_S;
  abcd[20] = I_ERI_H2xy2z_Px_G4x_S+ABY*I_ERI_G2x2z_Px_G4x_S;
  abcd[21] = I_ERI_Hx4y_Px_G4x_S+ABY*I_ERI_Gx3y_Px_G4x_S;
  abcd[22] = I_ERI_Hx3yz_Px_G4x_S+ABY*I_ERI_Gx2yz_Px_G4x_S;
  abcd[23] = I_ERI_Hx2y2z_Px_G4x_S+ABY*I_ERI_Gxy2z_Px_G4x_S;
  abcd[24] = I_ERI_Hxy3z_Px_G4x_S+ABY*I_ERI_Gx3z_Px_G4x_S;
  abcd[25] = I_ERI_H5y_Px_G4x_S+ABY*I_ERI_G4y_Px_G4x_S;
  abcd[26] = I_ERI_H4yz_Px_G4x_S+ABY*I_ERI_G3yz_Px_G4x_S;
  abcd[27] = I_ERI_H3y2z_Px_G4x_S+ABY*I_ERI_G2y2z_Px_G4x_S;
  abcd[28] = I_ERI_H2y3z_Px_G4x_S+ABY*I_ERI_Gy3z_Px_G4x_S;
  abcd[29] = I_ERI_Hy4z_Px_G4x_S+ABY*I_ERI_G4z_Px_G4x_S;
  abcd[30] = I_ERI_H4xz_Px_G4x_S+ABZ*I_ERI_G4x_Px_G4x_S;
  abcd[31] = I_ERI_H3xyz_Px_G4x_S+ABZ*I_ERI_G3xy_Px_G4x_S;
  abcd[32] = I_ERI_H3x2z_Px_G4x_S+ABZ*I_ERI_G3xz_Px_G4x_S;
  abcd[33] = I_ERI_H2x2yz_Px_G4x_S+ABZ*I_ERI_G2x2y_Px_G4x_S;
  abcd[34] = I_ERI_H2xy2z_Px_G4x_S+ABZ*I_ERI_G2xyz_Px_G4x_S;
  abcd[35] = I_ERI_H2x3z_Px_G4x_S+ABZ*I_ERI_G2x2z_Px_G4x_S;
  abcd[36] = I_ERI_Hx3yz_Px_G4x_S+ABZ*I_ERI_Gx3y_Px_G4x_S;
  abcd[37] = I_ERI_Hx2y2z_Px_G4x_S+ABZ*I_ERI_Gx2yz_Px_G4x_S;
  abcd[38] = I_ERI_Hxy3z_Px_G4x_S+ABZ*I_ERI_Gxy2z_Px_G4x_S;
  abcd[39] = I_ERI_Hx4z_Px_G4x_S+ABZ*I_ERI_Gx3z_Px_G4x_S;
  abcd[40] = I_ERI_H4yz_Px_G4x_S+ABZ*I_ERI_G4y_Px_G4x_S;
  abcd[41] = I_ERI_H3y2z_Px_G4x_S+ABZ*I_ERI_G3yz_Px_G4x_S;
  abcd[42] = I_ERI_H2y3z_Px_G4x_S+ABZ*I_ERI_G2y2z_Px_G4x_S;
  abcd[43] = I_ERI_Hy4z_Px_G4x_S+ABZ*I_ERI_Gy3z_Px_G4x_S;
  abcd[44] = I_ERI_H5z_Px_G4x_S+ABZ*I_ERI_G4z_Px_G4x_S;
  abcd[45] = I_ERI_H4xy_Py_G4x_S+ABY*I_ERI_G4x_Py_G4x_S;
  abcd[46] = I_ERI_H3x2y_Py_G4x_S+ABY*I_ERI_G3xy_Py_G4x_S;
  abcd[47] = I_ERI_H3xyz_Py_G4x_S+ABY*I_ERI_G3xz_Py_G4x_S;
  abcd[48] = I_ERI_H2x3y_Py_G4x_S+ABY*I_ERI_G2x2y_Py_G4x_S;
  abcd[49] = I_ERI_H2x2yz_Py_G4x_S+ABY*I_ERI_G2xyz_Py_G4x_S;
  abcd[50] = I_ERI_H2xy2z_Py_G4x_S+ABY*I_ERI_G2x2z_Py_G4x_S;
  abcd[51] = I_ERI_Hx4y_Py_G4x_S+ABY*I_ERI_Gx3y_Py_G4x_S;
  abcd[52] = I_ERI_Hx3yz_Py_G4x_S+ABY*I_ERI_Gx2yz_Py_G4x_S;
  abcd[53] = I_ERI_Hx2y2z_Py_G4x_S+ABY*I_ERI_Gxy2z_Py_G4x_S;
  abcd[54] = I_ERI_Hxy3z_Py_G4x_S+ABY*I_ERI_Gx3z_Py_G4x_S;
  abcd[55] = I_ERI_H5y_Py_G4x_S+ABY*I_ERI_G4y_Py_G4x_S;
  abcd[56] = I_ERI_H4yz_Py_G4x_S+ABY*I_ERI_G3yz_Py_G4x_S;
  abcd[57] = I_ERI_H3y2z_Py_G4x_S+ABY*I_ERI_G2y2z_Py_G4x_S;
  abcd[58] = I_ERI_H2y3z_Py_G4x_S+ABY*I_ERI_Gy3z_Py_G4x_S;
  abcd[59] = I_ERI_Hy4z_Py_G4x_S+ABY*I_ERI_G4z_Py_G4x_S;
  abcd[60] = I_ERI_H4xz_Py_G4x_S+ABZ*I_ERI_G4x_Py_G4x_S;
  abcd[61] = I_ERI_H3xyz_Py_G4x_S+ABZ*I_ERI_G3xy_Py_G4x_S;
  abcd[62] = I_ERI_H3x2z_Py_G4x_S+ABZ*I_ERI_G3xz_Py_G4x_S;
  abcd[63] = I_ERI_H2x2yz_Py_G4x_S+ABZ*I_ERI_G2x2y_Py_G4x_S;
  abcd[64] = I_ERI_H2xy2z_Py_G4x_S+ABZ*I_ERI_G2xyz_Py_G4x_S;
  abcd[65] = I_ERI_H2x3z_Py_G4x_S+ABZ*I_ERI_G2x2z_Py_G4x_S;
  abcd[66] = I_ERI_Hx3yz_Py_G4x_S+ABZ*I_ERI_Gx3y_Py_G4x_S;
  abcd[67] = I_ERI_Hx2y2z_Py_G4x_S+ABZ*I_ERI_Gx2yz_Py_G4x_S;
  abcd[68] = I_ERI_Hxy3z_Py_G4x_S+ABZ*I_ERI_Gxy2z_Py_G4x_S;
  abcd[69] = I_ERI_Hx4z_Py_G4x_S+ABZ*I_ERI_Gx3z_Py_G4x_S;
  abcd[70] = I_ERI_H4yz_Py_G4x_S+ABZ*I_ERI_G4y_Py_G4x_S;
  abcd[71] = I_ERI_H3y2z_Py_G4x_S+ABZ*I_ERI_G3yz_Py_G4x_S;
  abcd[72] = I_ERI_H2y3z_Py_G4x_S+ABZ*I_ERI_G2y2z_Py_G4x_S;
  abcd[73] = I_ERI_Hy4z_Py_G4x_S+ABZ*I_ERI_Gy3z_Py_G4x_S;
  abcd[74] = I_ERI_H5z_Py_G4x_S+ABZ*I_ERI_G4z_Py_G4x_S;
  abcd[75] = I_ERI_H4xz_Pz_G4x_S+ABZ*I_ERI_G4x_Pz_G4x_S;
  abcd[76] = I_ERI_H3xyz_Pz_G4x_S+ABZ*I_ERI_G3xy_Pz_G4x_S;
  abcd[77] = I_ERI_H3x2z_Pz_G4x_S+ABZ*I_ERI_G3xz_Pz_G4x_S;
  abcd[78] = I_ERI_H2x2yz_Pz_G4x_S+ABZ*I_ERI_G2x2y_Pz_G4x_S;
  abcd[79] = I_ERI_H2xy2z_Pz_G4x_S+ABZ*I_ERI_G2xyz_Pz_G4x_S;
  abcd[80] = I_ERI_H2x3z_Pz_G4x_S+ABZ*I_ERI_G2x2z_Pz_G4x_S;
  abcd[81] = I_ERI_Hx3yz_Pz_G4x_S+ABZ*I_ERI_Gx3y_Pz_G4x_S;
  abcd[82] = I_ERI_Hx2y2z_Pz_G4x_S+ABZ*I_ERI_Gx2yz_Pz_G4x_S;
  abcd[83] = I_ERI_Hxy3z_Pz_G4x_S+ABZ*I_ERI_Gxy2z_Pz_G4x_S;
  abcd[84] = I_ERI_Hx4z_Pz_G4x_S+ABZ*I_ERI_Gx3z_Pz_G4x_S;
  abcd[85] = I_ERI_H4yz_Pz_G4x_S+ABZ*I_ERI_G4y_Pz_G4x_S;
  abcd[86] = I_ERI_H3y2z_Pz_G4x_S+ABZ*I_ERI_G3yz_Pz_G4x_S;
  abcd[87] = I_ERI_H2y3z_Pz_G4x_S+ABZ*I_ERI_G2y2z_Pz_G4x_S;
  abcd[88] = I_ERI_Hy4z_Pz_G4x_S+ABZ*I_ERI_Gy3z_Pz_G4x_S;
  abcd[89] = I_ERI_H5z_Pz_G4x_S+ABZ*I_ERI_G4z_Pz_G4x_S;
  abcd[90] = I_ERI_H5x_Px_G3xy_S+ABX*I_ERI_G4x_Px_G3xy_S;
  abcd[91] = I_ERI_H4xy_Px_G3xy_S+ABX*I_ERI_G3xy_Px_G3xy_S;
  abcd[92] = I_ERI_H4xz_Px_G3xy_S+ABX*I_ERI_G3xz_Px_G3xy_S;
  abcd[93] = I_ERI_H3x2y_Px_G3xy_S+ABX*I_ERI_G2x2y_Px_G3xy_S;
  abcd[94] = I_ERI_H3xyz_Px_G3xy_S+ABX*I_ERI_G2xyz_Px_G3xy_S;
  abcd[95] = I_ERI_H3x2z_Px_G3xy_S+ABX*I_ERI_G2x2z_Px_G3xy_S;
  abcd[96] = I_ERI_H2x3y_Px_G3xy_S+ABX*I_ERI_Gx3y_Px_G3xy_S;
  abcd[97] = I_ERI_H2x2yz_Px_G3xy_S+ABX*I_ERI_Gx2yz_Px_G3xy_S;
  abcd[98] = I_ERI_H2xy2z_Px_G3xy_S+ABX*I_ERI_Gxy2z_Px_G3xy_S;
  abcd[99] = I_ERI_H2x3z_Px_G3xy_S+ABX*I_ERI_Gx3z_Px_G3xy_S;
  abcd[100] = I_ERI_Hx4y_Px_G3xy_S+ABX*I_ERI_G4y_Px_G3xy_S;
  abcd[101] = I_ERI_Hx3yz_Px_G3xy_S+ABX*I_ERI_G3yz_Px_G3xy_S;
  abcd[102] = I_ERI_Hx2y2z_Px_G3xy_S+ABX*I_ERI_G2y2z_Px_G3xy_S;
  abcd[103] = I_ERI_Hxy3z_Px_G3xy_S+ABX*I_ERI_Gy3z_Px_G3xy_S;
  abcd[104] = I_ERI_Hx4z_Px_G3xy_S+ABX*I_ERI_G4z_Px_G3xy_S;
  abcd[105] = I_ERI_H4xy_Px_G3xy_S+ABY*I_ERI_G4x_Px_G3xy_S;
  abcd[106] = I_ERI_H3x2y_Px_G3xy_S+ABY*I_ERI_G3xy_Px_G3xy_S;
  abcd[107] = I_ERI_H3xyz_Px_G3xy_S+ABY*I_ERI_G3xz_Px_G3xy_S;
  abcd[108] = I_ERI_H2x3y_Px_G3xy_S+ABY*I_ERI_G2x2y_Px_G3xy_S;
  abcd[109] = I_ERI_H2x2yz_Px_G3xy_S+ABY*I_ERI_G2xyz_Px_G3xy_S;
  abcd[110] = I_ERI_H2xy2z_Px_G3xy_S+ABY*I_ERI_G2x2z_Px_G3xy_S;
  abcd[111] = I_ERI_Hx4y_Px_G3xy_S+ABY*I_ERI_Gx3y_Px_G3xy_S;
  abcd[112] = I_ERI_Hx3yz_Px_G3xy_S+ABY*I_ERI_Gx2yz_Px_G3xy_S;
  abcd[113] = I_ERI_Hx2y2z_Px_G3xy_S+ABY*I_ERI_Gxy2z_Px_G3xy_S;
  abcd[114] = I_ERI_Hxy3z_Px_G3xy_S+ABY*I_ERI_Gx3z_Px_G3xy_S;
  abcd[115] = I_ERI_H5y_Px_G3xy_S+ABY*I_ERI_G4y_Px_G3xy_S;
  abcd[116] = I_ERI_H4yz_Px_G3xy_S+ABY*I_ERI_G3yz_Px_G3xy_S;
  abcd[117] = I_ERI_H3y2z_Px_G3xy_S+ABY*I_ERI_G2y2z_Px_G3xy_S;
  abcd[118] = I_ERI_H2y3z_Px_G3xy_S+ABY*I_ERI_Gy3z_Px_G3xy_S;
  abcd[119] = I_ERI_Hy4z_Px_G3xy_S+ABY*I_ERI_G4z_Px_G3xy_S;
  abcd[120] = I_ERI_H4xz_Px_G3xy_S+ABZ*I_ERI_G4x_Px_G3xy_S;
  abcd[121] = I_ERI_H3xyz_Px_G3xy_S+ABZ*I_ERI_G3xy_Px_G3xy_S;
  abcd[122] = I_ERI_H3x2z_Px_G3xy_S+ABZ*I_ERI_G3xz_Px_G3xy_S;
  abcd[123] = I_ERI_H2x2yz_Px_G3xy_S+ABZ*I_ERI_G2x2y_Px_G3xy_S;
  abcd[124] = I_ERI_H2xy2z_Px_G3xy_S+ABZ*I_ERI_G2xyz_Px_G3xy_S;
  abcd[125] = I_ERI_H2x3z_Px_G3xy_S+ABZ*I_ERI_G2x2z_Px_G3xy_S;
  abcd[126] = I_ERI_Hx3yz_Px_G3xy_S+ABZ*I_ERI_Gx3y_Px_G3xy_S;
  abcd[127] = I_ERI_Hx2y2z_Px_G3xy_S+ABZ*I_ERI_Gx2yz_Px_G3xy_S;
  abcd[128] = I_ERI_Hxy3z_Px_G3xy_S+ABZ*I_ERI_Gxy2z_Px_G3xy_S;
  abcd[129] = I_ERI_Hx4z_Px_G3xy_S+ABZ*I_ERI_Gx3z_Px_G3xy_S;
  abcd[130] = I_ERI_H4yz_Px_G3xy_S+ABZ*I_ERI_G4y_Px_G3xy_S;
  abcd[131] = I_ERI_H3y2z_Px_G3xy_S+ABZ*I_ERI_G3yz_Px_G3xy_S;
  abcd[132] = I_ERI_H2y3z_Px_G3xy_S+ABZ*I_ERI_G2y2z_Px_G3xy_S;
  abcd[133] = I_ERI_Hy4z_Px_G3xy_S+ABZ*I_ERI_Gy3z_Px_G3xy_S;
  abcd[134] = I_ERI_H5z_Px_G3xy_S+ABZ*I_ERI_G4z_Px_G3xy_S;
  abcd[135] = I_ERI_H4xy_Py_G3xy_S+ABY*I_ERI_G4x_Py_G3xy_S;
  abcd[136] = I_ERI_H3x2y_Py_G3xy_S+ABY*I_ERI_G3xy_Py_G3xy_S;
  abcd[137] = I_ERI_H3xyz_Py_G3xy_S+ABY*I_ERI_G3xz_Py_G3xy_S;
  abcd[138] = I_ERI_H2x3y_Py_G3xy_S+ABY*I_ERI_G2x2y_Py_G3xy_S;
  abcd[139] = I_ERI_H2x2yz_Py_G3xy_S+ABY*I_ERI_G2xyz_Py_G3xy_S;
  abcd[140] = I_ERI_H2xy2z_Py_G3xy_S+ABY*I_ERI_G2x2z_Py_G3xy_S;
  abcd[141] = I_ERI_Hx4y_Py_G3xy_S+ABY*I_ERI_Gx3y_Py_G3xy_S;
  abcd[142] = I_ERI_Hx3yz_Py_G3xy_S+ABY*I_ERI_Gx2yz_Py_G3xy_S;
  abcd[143] = I_ERI_Hx2y2z_Py_G3xy_S+ABY*I_ERI_Gxy2z_Py_G3xy_S;
  abcd[144] = I_ERI_Hxy3z_Py_G3xy_S+ABY*I_ERI_Gx3z_Py_G3xy_S;
  abcd[145] = I_ERI_H5y_Py_G3xy_S+ABY*I_ERI_G4y_Py_G3xy_S;
  abcd[146] = I_ERI_H4yz_Py_G3xy_S+ABY*I_ERI_G3yz_Py_G3xy_S;
  abcd[147] = I_ERI_H3y2z_Py_G3xy_S+ABY*I_ERI_G2y2z_Py_G3xy_S;
  abcd[148] = I_ERI_H2y3z_Py_G3xy_S+ABY*I_ERI_Gy3z_Py_G3xy_S;
  abcd[149] = I_ERI_Hy4z_Py_G3xy_S+ABY*I_ERI_G4z_Py_G3xy_S;
  abcd[150] = I_ERI_H4xz_Py_G3xy_S+ABZ*I_ERI_G4x_Py_G3xy_S;
  abcd[151] = I_ERI_H3xyz_Py_G3xy_S+ABZ*I_ERI_G3xy_Py_G3xy_S;
  abcd[152] = I_ERI_H3x2z_Py_G3xy_S+ABZ*I_ERI_G3xz_Py_G3xy_S;
  abcd[153] = I_ERI_H2x2yz_Py_G3xy_S+ABZ*I_ERI_G2x2y_Py_G3xy_S;
  abcd[154] = I_ERI_H2xy2z_Py_G3xy_S+ABZ*I_ERI_G2xyz_Py_G3xy_S;
  abcd[155] = I_ERI_H2x3z_Py_G3xy_S+ABZ*I_ERI_G2x2z_Py_G3xy_S;
  abcd[156] = I_ERI_Hx3yz_Py_G3xy_S+ABZ*I_ERI_Gx3y_Py_G3xy_S;
  abcd[157] = I_ERI_Hx2y2z_Py_G3xy_S+ABZ*I_ERI_Gx2yz_Py_G3xy_S;
  abcd[158] = I_ERI_Hxy3z_Py_G3xy_S+ABZ*I_ERI_Gxy2z_Py_G3xy_S;
  abcd[159] = I_ERI_Hx4z_Py_G3xy_S+ABZ*I_ERI_Gx3z_Py_G3xy_S;
  abcd[160] = I_ERI_H4yz_Py_G3xy_S+ABZ*I_ERI_G4y_Py_G3xy_S;
  abcd[161] = I_ERI_H3y2z_Py_G3xy_S+ABZ*I_ERI_G3yz_Py_G3xy_S;
  abcd[162] = I_ERI_H2y3z_Py_G3xy_S+ABZ*I_ERI_G2y2z_Py_G3xy_S;
  abcd[163] = I_ERI_Hy4z_Py_G3xy_S+ABZ*I_ERI_Gy3z_Py_G3xy_S;
  abcd[164] = I_ERI_H5z_Py_G3xy_S+ABZ*I_ERI_G4z_Py_G3xy_S;
  abcd[165] = I_ERI_H4xz_Pz_G3xy_S+ABZ*I_ERI_G4x_Pz_G3xy_S;
  abcd[166] = I_ERI_H3xyz_Pz_G3xy_S+ABZ*I_ERI_G3xy_Pz_G3xy_S;
  abcd[167] = I_ERI_H3x2z_Pz_G3xy_S+ABZ*I_ERI_G3xz_Pz_G3xy_S;
  abcd[168] = I_ERI_H2x2yz_Pz_G3xy_S+ABZ*I_ERI_G2x2y_Pz_G3xy_S;
  abcd[169] = I_ERI_H2xy2z_Pz_G3xy_S+ABZ*I_ERI_G2xyz_Pz_G3xy_S;
  abcd[170] = I_ERI_H2x3z_Pz_G3xy_S+ABZ*I_ERI_G2x2z_Pz_G3xy_S;
  abcd[171] = I_ERI_Hx3yz_Pz_G3xy_S+ABZ*I_ERI_Gx3y_Pz_G3xy_S;
  abcd[172] = I_ERI_Hx2y2z_Pz_G3xy_S+ABZ*I_ERI_Gx2yz_Pz_G3xy_S;
  abcd[173] = I_ERI_Hxy3z_Pz_G3xy_S+ABZ*I_ERI_Gxy2z_Pz_G3xy_S;
  abcd[174] = I_ERI_Hx4z_Pz_G3xy_S+ABZ*I_ERI_Gx3z_Pz_G3xy_S;
  abcd[175] = I_ERI_H4yz_Pz_G3xy_S+ABZ*I_ERI_G4y_Pz_G3xy_S;
  abcd[176] = I_ERI_H3y2z_Pz_G3xy_S+ABZ*I_ERI_G3yz_Pz_G3xy_S;
  abcd[177] = I_ERI_H2y3z_Pz_G3xy_S+ABZ*I_ERI_G2y2z_Pz_G3xy_S;
  abcd[178] = I_ERI_Hy4z_Pz_G3xy_S+ABZ*I_ERI_Gy3z_Pz_G3xy_S;
  abcd[179] = I_ERI_H5z_Pz_G3xy_S+ABZ*I_ERI_G4z_Pz_G3xy_S;
  abcd[180] = I_ERI_H5x_Px_G3xz_S+ABX*I_ERI_G4x_Px_G3xz_S;
  abcd[181] = I_ERI_H4xy_Px_G3xz_S+ABX*I_ERI_G3xy_Px_G3xz_S;
  abcd[182] = I_ERI_H4xz_Px_G3xz_S+ABX*I_ERI_G3xz_Px_G3xz_S;
  abcd[183] = I_ERI_H3x2y_Px_G3xz_S+ABX*I_ERI_G2x2y_Px_G3xz_S;
  abcd[184] = I_ERI_H3xyz_Px_G3xz_S+ABX*I_ERI_G2xyz_Px_G3xz_S;
  abcd[185] = I_ERI_H3x2z_Px_G3xz_S+ABX*I_ERI_G2x2z_Px_G3xz_S;
  abcd[186] = I_ERI_H2x3y_Px_G3xz_S+ABX*I_ERI_Gx3y_Px_G3xz_S;
  abcd[187] = I_ERI_H2x2yz_Px_G3xz_S+ABX*I_ERI_Gx2yz_Px_G3xz_S;
  abcd[188] = I_ERI_H2xy2z_Px_G3xz_S+ABX*I_ERI_Gxy2z_Px_G3xz_S;
  abcd[189] = I_ERI_H2x3z_Px_G3xz_S+ABX*I_ERI_Gx3z_Px_G3xz_S;
  abcd[190] = I_ERI_Hx4y_Px_G3xz_S+ABX*I_ERI_G4y_Px_G3xz_S;
  abcd[191] = I_ERI_Hx3yz_Px_G3xz_S+ABX*I_ERI_G3yz_Px_G3xz_S;
  abcd[192] = I_ERI_Hx2y2z_Px_G3xz_S+ABX*I_ERI_G2y2z_Px_G3xz_S;
  abcd[193] = I_ERI_Hxy3z_Px_G3xz_S+ABX*I_ERI_Gy3z_Px_G3xz_S;
  abcd[194] = I_ERI_Hx4z_Px_G3xz_S+ABX*I_ERI_G4z_Px_G3xz_S;
  abcd[195] = I_ERI_H4xy_Px_G3xz_S+ABY*I_ERI_G4x_Px_G3xz_S;
  abcd[196] = I_ERI_H3x2y_Px_G3xz_S+ABY*I_ERI_G3xy_Px_G3xz_S;
  abcd[197] = I_ERI_H3xyz_Px_G3xz_S+ABY*I_ERI_G3xz_Px_G3xz_S;
  abcd[198] = I_ERI_H2x3y_Px_G3xz_S+ABY*I_ERI_G2x2y_Px_G3xz_S;
  abcd[199] = I_ERI_H2x2yz_Px_G3xz_S+ABY*I_ERI_G2xyz_Px_G3xz_S;
  abcd[200] = I_ERI_H2xy2z_Px_G3xz_S+ABY*I_ERI_G2x2z_Px_G3xz_S;
  abcd[201] = I_ERI_Hx4y_Px_G3xz_S+ABY*I_ERI_Gx3y_Px_G3xz_S;
  abcd[202] = I_ERI_Hx3yz_Px_G3xz_S+ABY*I_ERI_Gx2yz_Px_G3xz_S;
  abcd[203] = I_ERI_Hx2y2z_Px_G3xz_S+ABY*I_ERI_Gxy2z_Px_G3xz_S;
  abcd[204] = I_ERI_Hxy3z_Px_G3xz_S+ABY*I_ERI_Gx3z_Px_G3xz_S;
  abcd[205] = I_ERI_H5y_Px_G3xz_S+ABY*I_ERI_G4y_Px_G3xz_S;
  abcd[206] = I_ERI_H4yz_Px_G3xz_S+ABY*I_ERI_G3yz_Px_G3xz_S;
  abcd[207] = I_ERI_H3y2z_Px_G3xz_S+ABY*I_ERI_G2y2z_Px_G3xz_S;
  abcd[208] = I_ERI_H2y3z_Px_G3xz_S+ABY*I_ERI_Gy3z_Px_G3xz_S;
  abcd[209] = I_ERI_Hy4z_Px_G3xz_S+ABY*I_ERI_G4z_Px_G3xz_S;
  abcd[210] = I_ERI_H4xz_Px_G3xz_S+ABZ*I_ERI_G4x_Px_G3xz_S;
  abcd[211] = I_ERI_H3xyz_Px_G3xz_S+ABZ*I_ERI_G3xy_Px_G3xz_S;
  abcd[212] = I_ERI_H3x2z_Px_G3xz_S+ABZ*I_ERI_G3xz_Px_G3xz_S;
  abcd[213] = I_ERI_H2x2yz_Px_G3xz_S+ABZ*I_ERI_G2x2y_Px_G3xz_S;
  abcd[214] = I_ERI_H2xy2z_Px_G3xz_S+ABZ*I_ERI_G2xyz_Px_G3xz_S;
  abcd[215] = I_ERI_H2x3z_Px_G3xz_S+ABZ*I_ERI_G2x2z_Px_G3xz_S;
  abcd[216] = I_ERI_Hx3yz_Px_G3xz_S+ABZ*I_ERI_Gx3y_Px_G3xz_S;
  abcd[217] = I_ERI_Hx2y2z_Px_G3xz_S+ABZ*I_ERI_Gx2yz_Px_G3xz_S;
  abcd[218] = I_ERI_Hxy3z_Px_G3xz_S+ABZ*I_ERI_Gxy2z_Px_G3xz_S;
  abcd[219] = I_ERI_Hx4z_Px_G3xz_S+ABZ*I_ERI_Gx3z_Px_G3xz_S;
  abcd[220] = I_ERI_H4yz_Px_G3xz_S+ABZ*I_ERI_G4y_Px_G3xz_S;
  abcd[221] = I_ERI_H3y2z_Px_G3xz_S+ABZ*I_ERI_G3yz_Px_G3xz_S;
  abcd[222] = I_ERI_H2y3z_Px_G3xz_S+ABZ*I_ERI_G2y2z_Px_G3xz_S;
  abcd[223] = I_ERI_Hy4z_Px_G3xz_S+ABZ*I_ERI_Gy3z_Px_G3xz_S;
  abcd[224] = I_ERI_H5z_Px_G3xz_S+ABZ*I_ERI_G4z_Px_G3xz_S;
  abcd[225] = I_ERI_H4xy_Py_G3xz_S+ABY*I_ERI_G4x_Py_G3xz_S;
  abcd[226] = I_ERI_H3x2y_Py_G3xz_S+ABY*I_ERI_G3xy_Py_G3xz_S;
  abcd[227] = I_ERI_H3xyz_Py_G3xz_S+ABY*I_ERI_G3xz_Py_G3xz_S;
  abcd[228] = I_ERI_H2x3y_Py_G3xz_S+ABY*I_ERI_G2x2y_Py_G3xz_S;
  abcd[229] = I_ERI_H2x2yz_Py_G3xz_S+ABY*I_ERI_G2xyz_Py_G3xz_S;
  abcd[230] = I_ERI_H2xy2z_Py_G3xz_S+ABY*I_ERI_G2x2z_Py_G3xz_S;
  abcd[231] = I_ERI_Hx4y_Py_G3xz_S+ABY*I_ERI_Gx3y_Py_G3xz_S;
  abcd[232] = I_ERI_Hx3yz_Py_G3xz_S+ABY*I_ERI_Gx2yz_Py_G3xz_S;
  abcd[233] = I_ERI_Hx2y2z_Py_G3xz_S+ABY*I_ERI_Gxy2z_Py_G3xz_S;
  abcd[234] = I_ERI_Hxy3z_Py_G3xz_S+ABY*I_ERI_Gx3z_Py_G3xz_S;
  abcd[235] = I_ERI_H5y_Py_G3xz_S+ABY*I_ERI_G4y_Py_G3xz_S;
  abcd[236] = I_ERI_H4yz_Py_G3xz_S+ABY*I_ERI_G3yz_Py_G3xz_S;
  abcd[237] = I_ERI_H3y2z_Py_G3xz_S+ABY*I_ERI_G2y2z_Py_G3xz_S;
  abcd[238] = I_ERI_H2y3z_Py_G3xz_S+ABY*I_ERI_Gy3z_Py_G3xz_S;
  abcd[239] = I_ERI_Hy4z_Py_G3xz_S+ABY*I_ERI_G4z_Py_G3xz_S;
  abcd[240] = I_ERI_H4xz_Py_G3xz_S+ABZ*I_ERI_G4x_Py_G3xz_S;
  abcd[241] = I_ERI_H3xyz_Py_G3xz_S+ABZ*I_ERI_G3xy_Py_G3xz_S;
  abcd[242] = I_ERI_H3x2z_Py_G3xz_S+ABZ*I_ERI_G3xz_Py_G3xz_S;
  abcd[243] = I_ERI_H2x2yz_Py_G3xz_S+ABZ*I_ERI_G2x2y_Py_G3xz_S;
  abcd[244] = I_ERI_H2xy2z_Py_G3xz_S+ABZ*I_ERI_G2xyz_Py_G3xz_S;
  abcd[245] = I_ERI_H2x3z_Py_G3xz_S+ABZ*I_ERI_G2x2z_Py_G3xz_S;
  abcd[246] = I_ERI_Hx3yz_Py_G3xz_S+ABZ*I_ERI_Gx3y_Py_G3xz_S;
  abcd[247] = I_ERI_Hx2y2z_Py_G3xz_S+ABZ*I_ERI_Gx2yz_Py_G3xz_S;
  abcd[248] = I_ERI_Hxy3z_Py_G3xz_S+ABZ*I_ERI_Gxy2z_Py_G3xz_S;
  abcd[249] = I_ERI_Hx4z_Py_G3xz_S+ABZ*I_ERI_Gx3z_Py_G3xz_S;
  abcd[250] = I_ERI_H4yz_Py_G3xz_S+ABZ*I_ERI_G4y_Py_G3xz_S;
  abcd[251] = I_ERI_H3y2z_Py_G3xz_S+ABZ*I_ERI_G3yz_Py_G3xz_S;
  abcd[252] = I_ERI_H2y3z_Py_G3xz_S+ABZ*I_ERI_G2y2z_Py_G3xz_S;
  abcd[253] = I_ERI_Hy4z_Py_G3xz_S+ABZ*I_ERI_Gy3z_Py_G3xz_S;
  abcd[254] = I_ERI_H5z_Py_G3xz_S+ABZ*I_ERI_G4z_Py_G3xz_S;
  abcd[255] = I_ERI_H4xz_Pz_G3xz_S+ABZ*I_ERI_G4x_Pz_G3xz_S;
  abcd[256] = I_ERI_H3xyz_Pz_G3xz_S+ABZ*I_ERI_G3xy_Pz_G3xz_S;
  abcd[257] = I_ERI_H3x2z_Pz_G3xz_S+ABZ*I_ERI_G3xz_Pz_G3xz_S;
  abcd[258] = I_ERI_H2x2yz_Pz_G3xz_S+ABZ*I_ERI_G2x2y_Pz_G3xz_S;
  abcd[259] = I_ERI_H2xy2z_Pz_G3xz_S+ABZ*I_ERI_G2xyz_Pz_G3xz_S;
  abcd[260] = I_ERI_H2x3z_Pz_G3xz_S+ABZ*I_ERI_G2x2z_Pz_G3xz_S;
  abcd[261] = I_ERI_Hx3yz_Pz_G3xz_S+ABZ*I_ERI_Gx3y_Pz_G3xz_S;
  abcd[262] = I_ERI_Hx2y2z_Pz_G3xz_S+ABZ*I_ERI_Gx2yz_Pz_G3xz_S;
  abcd[263] = I_ERI_Hxy3z_Pz_G3xz_S+ABZ*I_ERI_Gxy2z_Pz_G3xz_S;
  abcd[264] = I_ERI_Hx4z_Pz_G3xz_S+ABZ*I_ERI_Gx3z_Pz_G3xz_S;
  abcd[265] = I_ERI_H4yz_Pz_G3xz_S+ABZ*I_ERI_G4y_Pz_G3xz_S;
  abcd[266] = I_ERI_H3y2z_Pz_G3xz_S+ABZ*I_ERI_G3yz_Pz_G3xz_S;
  abcd[267] = I_ERI_H2y3z_Pz_G3xz_S+ABZ*I_ERI_G2y2z_Pz_G3xz_S;
  abcd[268] = I_ERI_Hy4z_Pz_G3xz_S+ABZ*I_ERI_Gy3z_Pz_G3xz_S;
  abcd[269] = I_ERI_H5z_Pz_G3xz_S+ABZ*I_ERI_G4z_Pz_G3xz_S;
  abcd[270] = I_ERI_H5x_Px_G2x2y_S+ABX*I_ERI_G4x_Px_G2x2y_S;
  abcd[271] = I_ERI_H4xy_Px_G2x2y_S+ABX*I_ERI_G3xy_Px_G2x2y_S;
  abcd[272] = I_ERI_H4xz_Px_G2x2y_S+ABX*I_ERI_G3xz_Px_G2x2y_S;
  abcd[273] = I_ERI_H3x2y_Px_G2x2y_S+ABX*I_ERI_G2x2y_Px_G2x2y_S;
  abcd[274] = I_ERI_H3xyz_Px_G2x2y_S+ABX*I_ERI_G2xyz_Px_G2x2y_S;
  abcd[275] = I_ERI_H3x2z_Px_G2x2y_S+ABX*I_ERI_G2x2z_Px_G2x2y_S;
  abcd[276] = I_ERI_H2x3y_Px_G2x2y_S+ABX*I_ERI_Gx3y_Px_G2x2y_S;
  abcd[277] = I_ERI_H2x2yz_Px_G2x2y_S+ABX*I_ERI_Gx2yz_Px_G2x2y_S;
  abcd[278] = I_ERI_H2xy2z_Px_G2x2y_S+ABX*I_ERI_Gxy2z_Px_G2x2y_S;
  abcd[279] = I_ERI_H2x3z_Px_G2x2y_S+ABX*I_ERI_Gx3z_Px_G2x2y_S;
  abcd[280] = I_ERI_Hx4y_Px_G2x2y_S+ABX*I_ERI_G4y_Px_G2x2y_S;
  abcd[281] = I_ERI_Hx3yz_Px_G2x2y_S+ABX*I_ERI_G3yz_Px_G2x2y_S;
  abcd[282] = I_ERI_Hx2y2z_Px_G2x2y_S+ABX*I_ERI_G2y2z_Px_G2x2y_S;
  abcd[283] = I_ERI_Hxy3z_Px_G2x2y_S+ABX*I_ERI_Gy3z_Px_G2x2y_S;
  abcd[284] = I_ERI_Hx4z_Px_G2x2y_S+ABX*I_ERI_G4z_Px_G2x2y_S;
  abcd[285] = I_ERI_H4xy_Px_G2x2y_S+ABY*I_ERI_G4x_Px_G2x2y_S;
  abcd[286] = I_ERI_H3x2y_Px_G2x2y_S+ABY*I_ERI_G3xy_Px_G2x2y_S;
  abcd[287] = I_ERI_H3xyz_Px_G2x2y_S+ABY*I_ERI_G3xz_Px_G2x2y_S;
  abcd[288] = I_ERI_H2x3y_Px_G2x2y_S+ABY*I_ERI_G2x2y_Px_G2x2y_S;
  abcd[289] = I_ERI_H2x2yz_Px_G2x2y_S+ABY*I_ERI_G2xyz_Px_G2x2y_S;
  abcd[290] = I_ERI_H2xy2z_Px_G2x2y_S+ABY*I_ERI_G2x2z_Px_G2x2y_S;
  abcd[291] = I_ERI_Hx4y_Px_G2x2y_S+ABY*I_ERI_Gx3y_Px_G2x2y_S;
  abcd[292] = I_ERI_Hx3yz_Px_G2x2y_S+ABY*I_ERI_Gx2yz_Px_G2x2y_S;
  abcd[293] = I_ERI_Hx2y2z_Px_G2x2y_S+ABY*I_ERI_Gxy2z_Px_G2x2y_S;
  abcd[294] = I_ERI_Hxy3z_Px_G2x2y_S+ABY*I_ERI_Gx3z_Px_G2x2y_S;
  abcd[295] = I_ERI_H5y_Px_G2x2y_S+ABY*I_ERI_G4y_Px_G2x2y_S;
  abcd[296] = I_ERI_H4yz_Px_G2x2y_S+ABY*I_ERI_G3yz_Px_G2x2y_S;
  abcd[297] = I_ERI_H3y2z_Px_G2x2y_S+ABY*I_ERI_G2y2z_Px_G2x2y_S;
  abcd[298] = I_ERI_H2y3z_Px_G2x2y_S+ABY*I_ERI_Gy3z_Px_G2x2y_S;
  abcd[299] = I_ERI_Hy4z_Px_G2x2y_S+ABY*I_ERI_G4z_Px_G2x2y_S;
  abcd[300] = I_ERI_H4xz_Px_G2x2y_S+ABZ*I_ERI_G4x_Px_G2x2y_S;
  abcd[301] = I_ERI_H3xyz_Px_G2x2y_S+ABZ*I_ERI_G3xy_Px_G2x2y_S;
  abcd[302] = I_ERI_H3x2z_Px_G2x2y_S+ABZ*I_ERI_G3xz_Px_G2x2y_S;
  abcd[303] = I_ERI_H2x2yz_Px_G2x2y_S+ABZ*I_ERI_G2x2y_Px_G2x2y_S;
  abcd[304] = I_ERI_H2xy2z_Px_G2x2y_S+ABZ*I_ERI_G2xyz_Px_G2x2y_S;
  abcd[305] = I_ERI_H2x3z_Px_G2x2y_S+ABZ*I_ERI_G2x2z_Px_G2x2y_S;
  abcd[306] = I_ERI_Hx3yz_Px_G2x2y_S+ABZ*I_ERI_Gx3y_Px_G2x2y_S;
  abcd[307] = I_ERI_Hx2y2z_Px_G2x2y_S+ABZ*I_ERI_Gx2yz_Px_G2x2y_S;
  abcd[308] = I_ERI_Hxy3z_Px_G2x2y_S+ABZ*I_ERI_Gxy2z_Px_G2x2y_S;
  abcd[309] = I_ERI_Hx4z_Px_G2x2y_S+ABZ*I_ERI_Gx3z_Px_G2x2y_S;
  abcd[310] = I_ERI_H4yz_Px_G2x2y_S+ABZ*I_ERI_G4y_Px_G2x2y_S;
  abcd[311] = I_ERI_H3y2z_Px_G2x2y_S+ABZ*I_ERI_G3yz_Px_G2x2y_S;
  abcd[312] = I_ERI_H2y3z_Px_G2x2y_S+ABZ*I_ERI_G2y2z_Px_G2x2y_S;
  abcd[313] = I_ERI_Hy4z_Px_G2x2y_S+ABZ*I_ERI_Gy3z_Px_G2x2y_S;
  abcd[314] = I_ERI_H5z_Px_G2x2y_S+ABZ*I_ERI_G4z_Px_G2x2y_S;
  abcd[315] = I_ERI_H4xy_Py_G2x2y_S+ABY*I_ERI_G4x_Py_G2x2y_S;
  abcd[316] = I_ERI_H3x2y_Py_G2x2y_S+ABY*I_ERI_G3xy_Py_G2x2y_S;
  abcd[317] = I_ERI_H3xyz_Py_G2x2y_S+ABY*I_ERI_G3xz_Py_G2x2y_S;
  abcd[318] = I_ERI_H2x3y_Py_G2x2y_S+ABY*I_ERI_G2x2y_Py_G2x2y_S;
  abcd[319] = I_ERI_H2x2yz_Py_G2x2y_S+ABY*I_ERI_G2xyz_Py_G2x2y_S;
  abcd[320] = I_ERI_H2xy2z_Py_G2x2y_S+ABY*I_ERI_G2x2z_Py_G2x2y_S;
  abcd[321] = I_ERI_Hx4y_Py_G2x2y_S+ABY*I_ERI_Gx3y_Py_G2x2y_S;
  abcd[322] = I_ERI_Hx3yz_Py_G2x2y_S+ABY*I_ERI_Gx2yz_Py_G2x2y_S;
  abcd[323] = I_ERI_Hx2y2z_Py_G2x2y_S+ABY*I_ERI_Gxy2z_Py_G2x2y_S;
  abcd[324] = I_ERI_Hxy3z_Py_G2x2y_S+ABY*I_ERI_Gx3z_Py_G2x2y_S;
  abcd[325] = I_ERI_H5y_Py_G2x2y_S+ABY*I_ERI_G4y_Py_G2x2y_S;
  abcd[326] = I_ERI_H4yz_Py_G2x2y_S+ABY*I_ERI_G3yz_Py_G2x2y_S;
  abcd[327] = I_ERI_H3y2z_Py_G2x2y_S+ABY*I_ERI_G2y2z_Py_G2x2y_S;
  abcd[328] = I_ERI_H2y3z_Py_G2x2y_S+ABY*I_ERI_Gy3z_Py_G2x2y_S;
  abcd[329] = I_ERI_Hy4z_Py_G2x2y_S+ABY*I_ERI_G4z_Py_G2x2y_S;
  abcd[330] = I_ERI_H4xz_Py_G2x2y_S+ABZ*I_ERI_G4x_Py_G2x2y_S;
  abcd[331] = I_ERI_H3xyz_Py_G2x2y_S+ABZ*I_ERI_G3xy_Py_G2x2y_S;
  abcd[332] = I_ERI_H3x2z_Py_G2x2y_S+ABZ*I_ERI_G3xz_Py_G2x2y_S;
  abcd[333] = I_ERI_H2x2yz_Py_G2x2y_S+ABZ*I_ERI_G2x2y_Py_G2x2y_S;
  abcd[334] = I_ERI_H2xy2z_Py_G2x2y_S+ABZ*I_ERI_G2xyz_Py_G2x2y_S;
  abcd[335] = I_ERI_H2x3z_Py_G2x2y_S+ABZ*I_ERI_G2x2z_Py_G2x2y_S;
  abcd[336] = I_ERI_Hx3yz_Py_G2x2y_S+ABZ*I_ERI_Gx3y_Py_G2x2y_S;
  abcd[337] = I_ERI_Hx2y2z_Py_G2x2y_S+ABZ*I_ERI_Gx2yz_Py_G2x2y_S;
  abcd[338] = I_ERI_Hxy3z_Py_G2x2y_S+ABZ*I_ERI_Gxy2z_Py_G2x2y_S;
  abcd[339] = I_ERI_Hx4z_Py_G2x2y_S+ABZ*I_ERI_Gx3z_Py_G2x2y_S;
  abcd[340] = I_ERI_H4yz_Py_G2x2y_S+ABZ*I_ERI_G4y_Py_G2x2y_S;
  abcd[341] = I_ERI_H3y2z_Py_G2x2y_S+ABZ*I_ERI_G3yz_Py_G2x2y_S;
  abcd[342] = I_ERI_H2y3z_Py_G2x2y_S+ABZ*I_ERI_G2y2z_Py_G2x2y_S;
  abcd[343] = I_ERI_Hy4z_Py_G2x2y_S+ABZ*I_ERI_Gy3z_Py_G2x2y_S;
  abcd[344] = I_ERI_H5z_Py_G2x2y_S+ABZ*I_ERI_G4z_Py_G2x2y_S;
  abcd[345] = I_ERI_H4xz_Pz_G2x2y_S+ABZ*I_ERI_G4x_Pz_G2x2y_S;
  abcd[346] = I_ERI_H3xyz_Pz_G2x2y_S+ABZ*I_ERI_G3xy_Pz_G2x2y_S;
  abcd[347] = I_ERI_H3x2z_Pz_G2x2y_S+ABZ*I_ERI_G3xz_Pz_G2x2y_S;
  abcd[348] = I_ERI_H2x2yz_Pz_G2x2y_S+ABZ*I_ERI_G2x2y_Pz_G2x2y_S;
  abcd[349] = I_ERI_H2xy2z_Pz_G2x2y_S+ABZ*I_ERI_G2xyz_Pz_G2x2y_S;
  abcd[350] = I_ERI_H2x3z_Pz_G2x2y_S+ABZ*I_ERI_G2x2z_Pz_G2x2y_S;
  abcd[351] = I_ERI_Hx3yz_Pz_G2x2y_S+ABZ*I_ERI_Gx3y_Pz_G2x2y_S;
  abcd[352] = I_ERI_Hx2y2z_Pz_G2x2y_S+ABZ*I_ERI_Gx2yz_Pz_G2x2y_S;
  abcd[353] = I_ERI_Hxy3z_Pz_G2x2y_S+ABZ*I_ERI_Gxy2z_Pz_G2x2y_S;
  abcd[354] = I_ERI_Hx4z_Pz_G2x2y_S+ABZ*I_ERI_Gx3z_Pz_G2x2y_S;
  abcd[355] = I_ERI_H4yz_Pz_G2x2y_S+ABZ*I_ERI_G4y_Pz_G2x2y_S;
  abcd[356] = I_ERI_H3y2z_Pz_G2x2y_S+ABZ*I_ERI_G3yz_Pz_G2x2y_S;
  abcd[357] = I_ERI_H2y3z_Pz_G2x2y_S+ABZ*I_ERI_G2y2z_Pz_G2x2y_S;
  abcd[358] = I_ERI_Hy4z_Pz_G2x2y_S+ABZ*I_ERI_Gy3z_Pz_G2x2y_S;
  abcd[359] = I_ERI_H5z_Pz_G2x2y_S+ABZ*I_ERI_G4z_Pz_G2x2y_S;
  abcd[360] = I_ERI_H5x_Px_G2xyz_S+ABX*I_ERI_G4x_Px_G2xyz_S;
  abcd[361] = I_ERI_H4xy_Px_G2xyz_S+ABX*I_ERI_G3xy_Px_G2xyz_S;
  abcd[362] = I_ERI_H4xz_Px_G2xyz_S+ABX*I_ERI_G3xz_Px_G2xyz_S;
  abcd[363] = I_ERI_H3x2y_Px_G2xyz_S+ABX*I_ERI_G2x2y_Px_G2xyz_S;
  abcd[364] = I_ERI_H3xyz_Px_G2xyz_S+ABX*I_ERI_G2xyz_Px_G2xyz_S;
  abcd[365] = I_ERI_H3x2z_Px_G2xyz_S+ABX*I_ERI_G2x2z_Px_G2xyz_S;
  abcd[366] = I_ERI_H2x3y_Px_G2xyz_S+ABX*I_ERI_Gx3y_Px_G2xyz_S;
  abcd[367] = I_ERI_H2x2yz_Px_G2xyz_S+ABX*I_ERI_Gx2yz_Px_G2xyz_S;
  abcd[368] = I_ERI_H2xy2z_Px_G2xyz_S+ABX*I_ERI_Gxy2z_Px_G2xyz_S;
  abcd[369] = I_ERI_H2x3z_Px_G2xyz_S+ABX*I_ERI_Gx3z_Px_G2xyz_S;
  abcd[370] = I_ERI_Hx4y_Px_G2xyz_S+ABX*I_ERI_G4y_Px_G2xyz_S;
  abcd[371] = I_ERI_Hx3yz_Px_G2xyz_S+ABX*I_ERI_G3yz_Px_G2xyz_S;
  abcd[372] = I_ERI_Hx2y2z_Px_G2xyz_S+ABX*I_ERI_G2y2z_Px_G2xyz_S;
  abcd[373] = I_ERI_Hxy3z_Px_G2xyz_S+ABX*I_ERI_Gy3z_Px_G2xyz_S;
  abcd[374] = I_ERI_Hx4z_Px_G2xyz_S+ABX*I_ERI_G4z_Px_G2xyz_S;
  abcd[375] = I_ERI_H4xy_Px_G2xyz_S+ABY*I_ERI_G4x_Px_G2xyz_S;
  abcd[376] = I_ERI_H3x2y_Px_G2xyz_S+ABY*I_ERI_G3xy_Px_G2xyz_S;
  abcd[377] = I_ERI_H3xyz_Px_G2xyz_S+ABY*I_ERI_G3xz_Px_G2xyz_S;
  abcd[378] = I_ERI_H2x3y_Px_G2xyz_S+ABY*I_ERI_G2x2y_Px_G2xyz_S;
  abcd[379] = I_ERI_H2x2yz_Px_G2xyz_S+ABY*I_ERI_G2xyz_Px_G2xyz_S;
  abcd[380] = I_ERI_H2xy2z_Px_G2xyz_S+ABY*I_ERI_G2x2z_Px_G2xyz_S;
  abcd[381] = I_ERI_Hx4y_Px_G2xyz_S+ABY*I_ERI_Gx3y_Px_G2xyz_S;
  abcd[382] = I_ERI_Hx3yz_Px_G2xyz_S+ABY*I_ERI_Gx2yz_Px_G2xyz_S;
  abcd[383] = I_ERI_Hx2y2z_Px_G2xyz_S+ABY*I_ERI_Gxy2z_Px_G2xyz_S;
  abcd[384] = I_ERI_Hxy3z_Px_G2xyz_S+ABY*I_ERI_Gx3z_Px_G2xyz_S;
  abcd[385] = I_ERI_H5y_Px_G2xyz_S+ABY*I_ERI_G4y_Px_G2xyz_S;
  abcd[386] = I_ERI_H4yz_Px_G2xyz_S+ABY*I_ERI_G3yz_Px_G2xyz_S;
  abcd[387] = I_ERI_H3y2z_Px_G2xyz_S+ABY*I_ERI_G2y2z_Px_G2xyz_S;
  abcd[388] = I_ERI_H2y3z_Px_G2xyz_S+ABY*I_ERI_Gy3z_Px_G2xyz_S;
  abcd[389] = I_ERI_Hy4z_Px_G2xyz_S+ABY*I_ERI_G4z_Px_G2xyz_S;
  abcd[390] = I_ERI_H4xz_Px_G2xyz_S+ABZ*I_ERI_G4x_Px_G2xyz_S;
  abcd[391] = I_ERI_H3xyz_Px_G2xyz_S+ABZ*I_ERI_G3xy_Px_G2xyz_S;
  abcd[392] = I_ERI_H3x2z_Px_G2xyz_S+ABZ*I_ERI_G3xz_Px_G2xyz_S;
  abcd[393] = I_ERI_H2x2yz_Px_G2xyz_S+ABZ*I_ERI_G2x2y_Px_G2xyz_S;
  abcd[394] = I_ERI_H2xy2z_Px_G2xyz_S+ABZ*I_ERI_G2xyz_Px_G2xyz_S;
  abcd[395] = I_ERI_H2x3z_Px_G2xyz_S+ABZ*I_ERI_G2x2z_Px_G2xyz_S;
  abcd[396] = I_ERI_Hx3yz_Px_G2xyz_S+ABZ*I_ERI_Gx3y_Px_G2xyz_S;
  abcd[397] = I_ERI_Hx2y2z_Px_G2xyz_S+ABZ*I_ERI_Gx2yz_Px_G2xyz_S;
  abcd[398] = I_ERI_Hxy3z_Px_G2xyz_S+ABZ*I_ERI_Gxy2z_Px_G2xyz_S;
  abcd[399] = I_ERI_Hx4z_Px_G2xyz_S+ABZ*I_ERI_Gx3z_Px_G2xyz_S;
  abcd[400] = I_ERI_H4yz_Px_G2xyz_S+ABZ*I_ERI_G4y_Px_G2xyz_S;
  abcd[401] = I_ERI_H3y2z_Px_G2xyz_S+ABZ*I_ERI_G3yz_Px_G2xyz_S;
  abcd[402] = I_ERI_H2y3z_Px_G2xyz_S+ABZ*I_ERI_G2y2z_Px_G2xyz_S;
  abcd[403] = I_ERI_Hy4z_Px_G2xyz_S+ABZ*I_ERI_Gy3z_Px_G2xyz_S;
  abcd[404] = I_ERI_H5z_Px_G2xyz_S+ABZ*I_ERI_G4z_Px_G2xyz_S;
  abcd[405] = I_ERI_H4xy_Py_G2xyz_S+ABY*I_ERI_G4x_Py_G2xyz_S;
  abcd[406] = I_ERI_H3x2y_Py_G2xyz_S+ABY*I_ERI_G3xy_Py_G2xyz_S;
  abcd[407] = I_ERI_H3xyz_Py_G2xyz_S+ABY*I_ERI_G3xz_Py_G2xyz_S;
  abcd[408] = I_ERI_H2x3y_Py_G2xyz_S+ABY*I_ERI_G2x2y_Py_G2xyz_S;
  abcd[409] = I_ERI_H2x2yz_Py_G2xyz_S+ABY*I_ERI_G2xyz_Py_G2xyz_S;
  abcd[410] = I_ERI_H2xy2z_Py_G2xyz_S+ABY*I_ERI_G2x2z_Py_G2xyz_S;
  abcd[411] = I_ERI_Hx4y_Py_G2xyz_S+ABY*I_ERI_Gx3y_Py_G2xyz_S;
  abcd[412] = I_ERI_Hx3yz_Py_G2xyz_S+ABY*I_ERI_Gx2yz_Py_G2xyz_S;
  abcd[413] = I_ERI_Hx2y2z_Py_G2xyz_S+ABY*I_ERI_Gxy2z_Py_G2xyz_S;
  abcd[414] = I_ERI_Hxy3z_Py_G2xyz_S+ABY*I_ERI_Gx3z_Py_G2xyz_S;
  abcd[415] = I_ERI_H5y_Py_G2xyz_S+ABY*I_ERI_G4y_Py_G2xyz_S;
  abcd[416] = I_ERI_H4yz_Py_G2xyz_S+ABY*I_ERI_G3yz_Py_G2xyz_S;
  abcd[417] = I_ERI_H3y2z_Py_G2xyz_S+ABY*I_ERI_G2y2z_Py_G2xyz_S;
  abcd[418] = I_ERI_H2y3z_Py_G2xyz_S+ABY*I_ERI_Gy3z_Py_G2xyz_S;
  abcd[419] = I_ERI_Hy4z_Py_G2xyz_S+ABY*I_ERI_G4z_Py_G2xyz_S;
  abcd[420] = I_ERI_H4xz_Py_G2xyz_S+ABZ*I_ERI_G4x_Py_G2xyz_S;
  abcd[421] = I_ERI_H3xyz_Py_G2xyz_S+ABZ*I_ERI_G3xy_Py_G2xyz_S;
  abcd[422] = I_ERI_H3x2z_Py_G2xyz_S+ABZ*I_ERI_G3xz_Py_G2xyz_S;
  abcd[423] = I_ERI_H2x2yz_Py_G2xyz_S+ABZ*I_ERI_G2x2y_Py_G2xyz_S;
  abcd[424] = I_ERI_H2xy2z_Py_G2xyz_S+ABZ*I_ERI_G2xyz_Py_G2xyz_S;
  abcd[425] = I_ERI_H2x3z_Py_G2xyz_S+ABZ*I_ERI_G2x2z_Py_G2xyz_S;
  abcd[426] = I_ERI_Hx3yz_Py_G2xyz_S+ABZ*I_ERI_Gx3y_Py_G2xyz_S;
  abcd[427] = I_ERI_Hx2y2z_Py_G2xyz_S+ABZ*I_ERI_Gx2yz_Py_G2xyz_S;
  abcd[428] = I_ERI_Hxy3z_Py_G2xyz_S+ABZ*I_ERI_Gxy2z_Py_G2xyz_S;
  abcd[429] = I_ERI_Hx4z_Py_G2xyz_S+ABZ*I_ERI_Gx3z_Py_G2xyz_S;
  abcd[430] = I_ERI_H4yz_Py_G2xyz_S+ABZ*I_ERI_G4y_Py_G2xyz_S;
  abcd[431] = I_ERI_H3y2z_Py_G2xyz_S+ABZ*I_ERI_G3yz_Py_G2xyz_S;
  abcd[432] = I_ERI_H2y3z_Py_G2xyz_S+ABZ*I_ERI_G2y2z_Py_G2xyz_S;
  abcd[433] = I_ERI_Hy4z_Py_G2xyz_S+ABZ*I_ERI_Gy3z_Py_G2xyz_S;
  abcd[434] = I_ERI_H5z_Py_G2xyz_S+ABZ*I_ERI_G4z_Py_G2xyz_S;
  abcd[435] = I_ERI_H4xz_Pz_G2xyz_S+ABZ*I_ERI_G4x_Pz_G2xyz_S;
  abcd[436] = I_ERI_H3xyz_Pz_G2xyz_S+ABZ*I_ERI_G3xy_Pz_G2xyz_S;
  abcd[437] = I_ERI_H3x2z_Pz_G2xyz_S+ABZ*I_ERI_G3xz_Pz_G2xyz_S;
  abcd[438] = I_ERI_H2x2yz_Pz_G2xyz_S+ABZ*I_ERI_G2x2y_Pz_G2xyz_S;
  abcd[439] = I_ERI_H2xy2z_Pz_G2xyz_S+ABZ*I_ERI_G2xyz_Pz_G2xyz_S;
  abcd[440] = I_ERI_H2x3z_Pz_G2xyz_S+ABZ*I_ERI_G2x2z_Pz_G2xyz_S;
  abcd[441] = I_ERI_Hx3yz_Pz_G2xyz_S+ABZ*I_ERI_Gx3y_Pz_G2xyz_S;
  abcd[442] = I_ERI_Hx2y2z_Pz_G2xyz_S+ABZ*I_ERI_Gx2yz_Pz_G2xyz_S;
  abcd[443] = I_ERI_Hxy3z_Pz_G2xyz_S+ABZ*I_ERI_Gxy2z_Pz_G2xyz_S;
  abcd[444] = I_ERI_Hx4z_Pz_G2xyz_S+ABZ*I_ERI_Gx3z_Pz_G2xyz_S;
  abcd[445] = I_ERI_H4yz_Pz_G2xyz_S+ABZ*I_ERI_G4y_Pz_G2xyz_S;
  abcd[446] = I_ERI_H3y2z_Pz_G2xyz_S+ABZ*I_ERI_G3yz_Pz_G2xyz_S;
  abcd[447] = I_ERI_H2y3z_Pz_G2xyz_S+ABZ*I_ERI_G2y2z_Pz_G2xyz_S;
  abcd[448] = I_ERI_Hy4z_Pz_G2xyz_S+ABZ*I_ERI_Gy3z_Pz_G2xyz_S;
  abcd[449] = I_ERI_H5z_Pz_G2xyz_S+ABZ*I_ERI_G4z_Pz_G2xyz_S;
  abcd[450] = I_ERI_H5x_Px_G2x2z_S+ABX*I_ERI_G4x_Px_G2x2z_S;
  abcd[451] = I_ERI_H4xy_Px_G2x2z_S+ABX*I_ERI_G3xy_Px_G2x2z_S;
  abcd[452] = I_ERI_H4xz_Px_G2x2z_S+ABX*I_ERI_G3xz_Px_G2x2z_S;
  abcd[453] = I_ERI_H3x2y_Px_G2x2z_S+ABX*I_ERI_G2x2y_Px_G2x2z_S;
  abcd[454] = I_ERI_H3xyz_Px_G2x2z_S+ABX*I_ERI_G2xyz_Px_G2x2z_S;
  abcd[455] = I_ERI_H3x2z_Px_G2x2z_S+ABX*I_ERI_G2x2z_Px_G2x2z_S;
  abcd[456] = I_ERI_H2x3y_Px_G2x2z_S+ABX*I_ERI_Gx3y_Px_G2x2z_S;
  abcd[457] = I_ERI_H2x2yz_Px_G2x2z_S+ABX*I_ERI_Gx2yz_Px_G2x2z_S;
  abcd[458] = I_ERI_H2xy2z_Px_G2x2z_S+ABX*I_ERI_Gxy2z_Px_G2x2z_S;
  abcd[459] = I_ERI_H2x3z_Px_G2x2z_S+ABX*I_ERI_Gx3z_Px_G2x2z_S;
  abcd[460] = I_ERI_Hx4y_Px_G2x2z_S+ABX*I_ERI_G4y_Px_G2x2z_S;
  abcd[461] = I_ERI_Hx3yz_Px_G2x2z_S+ABX*I_ERI_G3yz_Px_G2x2z_S;
  abcd[462] = I_ERI_Hx2y2z_Px_G2x2z_S+ABX*I_ERI_G2y2z_Px_G2x2z_S;
  abcd[463] = I_ERI_Hxy3z_Px_G2x2z_S+ABX*I_ERI_Gy3z_Px_G2x2z_S;
  abcd[464] = I_ERI_Hx4z_Px_G2x2z_S+ABX*I_ERI_G4z_Px_G2x2z_S;
  abcd[465] = I_ERI_H4xy_Px_G2x2z_S+ABY*I_ERI_G4x_Px_G2x2z_S;
  abcd[466] = I_ERI_H3x2y_Px_G2x2z_S+ABY*I_ERI_G3xy_Px_G2x2z_S;
  abcd[467] = I_ERI_H3xyz_Px_G2x2z_S+ABY*I_ERI_G3xz_Px_G2x2z_S;
  abcd[468] = I_ERI_H2x3y_Px_G2x2z_S+ABY*I_ERI_G2x2y_Px_G2x2z_S;
  abcd[469] = I_ERI_H2x2yz_Px_G2x2z_S+ABY*I_ERI_G2xyz_Px_G2x2z_S;
  abcd[470] = I_ERI_H2xy2z_Px_G2x2z_S+ABY*I_ERI_G2x2z_Px_G2x2z_S;
  abcd[471] = I_ERI_Hx4y_Px_G2x2z_S+ABY*I_ERI_Gx3y_Px_G2x2z_S;
  abcd[472] = I_ERI_Hx3yz_Px_G2x2z_S+ABY*I_ERI_Gx2yz_Px_G2x2z_S;
  abcd[473] = I_ERI_Hx2y2z_Px_G2x2z_S+ABY*I_ERI_Gxy2z_Px_G2x2z_S;
  abcd[474] = I_ERI_Hxy3z_Px_G2x2z_S+ABY*I_ERI_Gx3z_Px_G2x2z_S;
  abcd[475] = I_ERI_H5y_Px_G2x2z_S+ABY*I_ERI_G4y_Px_G2x2z_S;
  abcd[476] = I_ERI_H4yz_Px_G2x2z_S+ABY*I_ERI_G3yz_Px_G2x2z_S;
  abcd[477] = I_ERI_H3y2z_Px_G2x2z_S+ABY*I_ERI_G2y2z_Px_G2x2z_S;
  abcd[478] = I_ERI_H2y3z_Px_G2x2z_S+ABY*I_ERI_Gy3z_Px_G2x2z_S;
  abcd[479] = I_ERI_Hy4z_Px_G2x2z_S+ABY*I_ERI_G4z_Px_G2x2z_S;
  abcd[480] = I_ERI_H4xz_Px_G2x2z_S+ABZ*I_ERI_G4x_Px_G2x2z_S;
  abcd[481] = I_ERI_H3xyz_Px_G2x2z_S+ABZ*I_ERI_G3xy_Px_G2x2z_S;
  abcd[482] = I_ERI_H3x2z_Px_G2x2z_S+ABZ*I_ERI_G3xz_Px_G2x2z_S;
  abcd[483] = I_ERI_H2x2yz_Px_G2x2z_S+ABZ*I_ERI_G2x2y_Px_G2x2z_S;
  abcd[484] = I_ERI_H2xy2z_Px_G2x2z_S+ABZ*I_ERI_G2xyz_Px_G2x2z_S;
  abcd[485] = I_ERI_H2x3z_Px_G2x2z_S+ABZ*I_ERI_G2x2z_Px_G2x2z_S;
  abcd[486] = I_ERI_Hx3yz_Px_G2x2z_S+ABZ*I_ERI_Gx3y_Px_G2x2z_S;
  abcd[487] = I_ERI_Hx2y2z_Px_G2x2z_S+ABZ*I_ERI_Gx2yz_Px_G2x2z_S;
  abcd[488] = I_ERI_Hxy3z_Px_G2x2z_S+ABZ*I_ERI_Gxy2z_Px_G2x2z_S;
  abcd[489] = I_ERI_Hx4z_Px_G2x2z_S+ABZ*I_ERI_Gx3z_Px_G2x2z_S;
  abcd[490] = I_ERI_H4yz_Px_G2x2z_S+ABZ*I_ERI_G4y_Px_G2x2z_S;
  abcd[491] = I_ERI_H3y2z_Px_G2x2z_S+ABZ*I_ERI_G3yz_Px_G2x2z_S;
  abcd[492] = I_ERI_H2y3z_Px_G2x2z_S+ABZ*I_ERI_G2y2z_Px_G2x2z_S;
  abcd[493] = I_ERI_Hy4z_Px_G2x2z_S+ABZ*I_ERI_Gy3z_Px_G2x2z_S;
  abcd[494] = I_ERI_H5z_Px_G2x2z_S+ABZ*I_ERI_G4z_Px_G2x2z_S;
  abcd[495] = I_ERI_H4xy_Py_G2x2z_S+ABY*I_ERI_G4x_Py_G2x2z_S;
  abcd[496] = I_ERI_H3x2y_Py_G2x2z_S+ABY*I_ERI_G3xy_Py_G2x2z_S;
  abcd[497] = I_ERI_H3xyz_Py_G2x2z_S+ABY*I_ERI_G3xz_Py_G2x2z_S;
  abcd[498] = I_ERI_H2x3y_Py_G2x2z_S+ABY*I_ERI_G2x2y_Py_G2x2z_S;
  abcd[499] = I_ERI_H2x2yz_Py_G2x2z_S+ABY*I_ERI_G2xyz_Py_G2x2z_S;
  abcd[500] = I_ERI_H2xy2z_Py_G2x2z_S+ABY*I_ERI_G2x2z_Py_G2x2z_S;
  abcd[501] = I_ERI_Hx4y_Py_G2x2z_S+ABY*I_ERI_Gx3y_Py_G2x2z_S;
  abcd[502] = I_ERI_Hx3yz_Py_G2x2z_S+ABY*I_ERI_Gx2yz_Py_G2x2z_S;
  abcd[503] = I_ERI_Hx2y2z_Py_G2x2z_S+ABY*I_ERI_Gxy2z_Py_G2x2z_S;
  abcd[504] = I_ERI_Hxy3z_Py_G2x2z_S+ABY*I_ERI_Gx3z_Py_G2x2z_S;
  abcd[505] = I_ERI_H5y_Py_G2x2z_S+ABY*I_ERI_G4y_Py_G2x2z_S;
  abcd[506] = I_ERI_H4yz_Py_G2x2z_S+ABY*I_ERI_G3yz_Py_G2x2z_S;
  abcd[507] = I_ERI_H3y2z_Py_G2x2z_S+ABY*I_ERI_G2y2z_Py_G2x2z_S;
  abcd[508] = I_ERI_H2y3z_Py_G2x2z_S+ABY*I_ERI_Gy3z_Py_G2x2z_S;
  abcd[509] = I_ERI_Hy4z_Py_G2x2z_S+ABY*I_ERI_G4z_Py_G2x2z_S;
  abcd[510] = I_ERI_H4xz_Py_G2x2z_S+ABZ*I_ERI_G4x_Py_G2x2z_S;
  abcd[511] = I_ERI_H3xyz_Py_G2x2z_S+ABZ*I_ERI_G3xy_Py_G2x2z_S;
  abcd[512] = I_ERI_H3x2z_Py_G2x2z_S+ABZ*I_ERI_G3xz_Py_G2x2z_S;
  abcd[513] = I_ERI_H2x2yz_Py_G2x2z_S+ABZ*I_ERI_G2x2y_Py_G2x2z_S;
  abcd[514] = I_ERI_H2xy2z_Py_G2x2z_S+ABZ*I_ERI_G2xyz_Py_G2x2z_S;
  abcd[515] = I_ERI_H2x3z_Py_G2x2z_S+ABZ*I_ERI_G2x2z_Py_G2x2z_S;
  abcd[516] = I_ERI_Hx3yz_Py_G2x2z_S+ABZ*I_ERI_Gx3y_Py_G2x2z_S;
  abcd[517] = I_ERI_Hx2y2z_Py_G2x2z_S+ABZ*I_ERI_Gx2yz_Py_G2x2z_S;
  abcd[518] = I_ERI_Hxy3z_Py_G2x2z_S+ABZ*I_ERI_Gxy2z_Py_G2x2z_S;
  abcd[519] = I_ERI_Hx4z_Py_G2x2z_S+ABZ*I_ERI_Gx3z_Py_G2x2z_S;
  abcd[520] = I_ERI_H4yz_Py_G2x2z_S+ABZ*I_ERI_G4y_Py_G2x2z_S;
  abcd[521] = I_ERI_H3y2z_Py_G2x2z_S+ABZ*I_ERI_G3yz_Py_G2x2z_S;
  abcd[522] = I_ERI_H2y3z_Py_G2x2z_S+ABZ*I_ERI_G2y2z_Py_G2x2z_S;
  abcd[523] = I_ERI_Hy4z_Py_G2x2z_S+ABZ*I_ERI_Gy3z_Py_G2x2z_S;
  abcd[524] = I_ERI_H5z_Py_G2x2z_S+ABZ*I_ERI_G4z_Py_G2x2z_S;
  abcd[525] = I_ERI_H4xz_Pz_G2x2z_S+ABZ*I_ERI_G4x_Pz_G2x2z_S;
  abcd[526] = I_ERI_H3xyz_Pz_G2x2z_S+ABZ*I_ERI_G3xy_Pz_G2x2z_S;
  abcd[527] = I_ERI_H3x2z_Pz_G2x2z_S+ABZ*I_ERI_G3xz_Pz_G2x2z_S;
  abcd[528] = I_ERI_H2x2yz_Pz_G2x2z_S+ABZ*I_ERI_G2x2y_Pz_G2x2z_S;
  abcd[529] = I_ERI_H2xy2z_Pz_G2x2z_S+ABZ*I_ERI_G2xyz_Pz_G2x2z_S;
  abcd[530] = I_ERI_H2x3z_Pz_G2x2z_S+ABZ*I_ERI_G2x2z_Pz_G2x2z_S;
  abcd[531] = I_ERI_Hx3yz_Pz_G2x2z_S+ABZ*I_ERI_Gx3y_Pz_G2x2z_S;
  abcd[532] = I_ERI_Hx2y2z_Pz_G2x2z_S+ABZ*I_ERI_Gx2yz_Pz_G2x2z_S;
  abcd[533] = I_ERI_Hxy3z_Pz_G2x2z_S+ABZ*I_ERI_Gxy2z_Pz_G2x2z_S;
  abcd[534] = I_ERI_Hx4z_Pz_G2x2z_S+ABZ*I_ERI_Gx3z_Pz_G2x2z_S;
  abcd[535] = I_ERI_H4yz_Pz_G2x2z_S+ABZ*I_ERI_G4y_Pz_G2x2z_S;
  abcd[536] = I_ERI_H3y2z_Pz_G2x2z_S+ABZ*I_ERI_G3yz_Pz_G2x2z_S;
  abcd[537] = I_ERI_H2y3z_Pz_G2x2z_S+ABZ*I_ERI_G2y2z_Pz_G2x2z_S;
  abcd[538] = I_ERI_Hy4z_Pz_G2x2z_S+ABZ*I_ERI_Gy3z_Pz_G2x2z_S;
  abcd[539] = I_ERI_H5z_Pz_G2x2z_S+ABZ*I_ERI_G4z_Pz_G2x2z_S;
  abcd[540] = I_ERI_H5x_Px_Gx3y_S+ABX*I_ERI_G4x_Px_Gx3y_S;
  abcd[541] = I_ERI_H4xy_Px_Gx3y_S+ABX*I_ERI_G3xy_Px_Gx3y_S;
  abcd[542] = I_ERI_H4xz_Px_Gx3y_S+ABX*I_ERI_G3xz_Px_Gx3y_S;
  abcd[543] = I_ERI_H3x2y_Px_Gx3y_S+ABX*I_ERI_G2x2y_Px_Gx3y_S;
  abcd[544] = I_ERI_H3xyz_Px_Gx3y_S+ABX*I_ERI_G2xyz_Px_Gx3y_S;
  abcd[545] = I_ERI_H3x2z_Px_Gx3y_S+ABX*I_ERI_G2x2z_Px_Gx3y_S;
  abcd[546] = I_ERI_H2x3y_Px_Gx3y_S+ABX*I_ERI_Gx3y_Px_Gx3y_S;
  abcd[547] = I_ERI_H2x2yz_Px_Gx3y_S+ABX*I_ERI_Gx2yz_Px_Gx3y_S;
  abcd[548] = I_ERI_H2xy2z_Px_Gx3y_S+ABX*I_ERI_Gxy2z_Px_Gx3y_S;
  abcd[549] = I_ERI_H2x3z_Px_Gx3y_S+ABX*I_ERI_Gx3z_Px_Gx3y_S;
  abcd[550] = I_ERI_Hx4y_Px_Gx3y_S+ABX*I_ERI_G4y_Px_Gx3y_S;
  abcd[551] = I_ERI_Hx3yz_Px_Gx3y_S+ABX*I_ERI_G3yz_Px_Gx3y_S;
  abcd[552] = I_ERI_Hx2y2z_Px_Gx3y_S+ABX*I_ERI_G2y2z_Px_Gx3y_S;
  abcd[553] = I_ERI_Hxy3z_Px_Gx3y_S+ABX*I_ERI_Gy3z_Px_Gx3y_S;
  abcd[554] = I_ERI_Hx4z_Px_Gx3y_S+ABX*I_ERI_G4z_Px_Gx3y_S;
  abcd[555] = I_ERI_H4xy_Px_Gx3y_S+ABY*I_ERI_G4x_Px_Gx3y_S;
  abcd[556] = I_ERI_H3x2y_Px_Gx3y_S+ABY*I_ERI_G3xy_Px_Gx3y_S;
  abcd[557] = I_ERI_H3xyz_Px_Gx3y_S+ABY*I_ERI_G3xz_Px_Gx3y_S;
  abcd[558] = I_ERI_H2x3y_Px_Gx3y_S+ABY*I_ERI_G2x2y_Px_Gx3y_S;
  abcd[559] = I_ERI_H2x2yz_Px_Gx3y_S+ABY*I_ERI_G2xyz_Px_Gx3y_S;
  abcd[560] = I_ERI_H2xy2z_Px_Gx3y_S+ABY*I_ERI_G2x2z_Px_Gx3y_S;
  abcd[561] = I_ERI_Hx4y_Px_Gx3y_S+ABY*I_ERI_Gx3y_Px_Gx3y_S;
  abcd[562] = I_ERI_Hx3yz_Px_Gx3y_S+ABY*I_ERI_Gx2yz_Px_Gx3y_S;
  abcd[563] = I_ERI_Hx2y2z_Px_Gx3y_S+ABY*I_ERI_Gxy2z_Px_Gx3y_S;
  abcd[564] = I_ERI_Hxy3z_Px_Gx3y_S+ABY*I_ERI_Gx3z_Px_Gx3y_S;
  abcd[565] = I_ERI_H5y_Px_Gx3y_S+ABY*I_ERI_G4y_Px_Gx3y_S;
  abcd[566] = I_ERI_H4yz_Px_Gx3y_S+ABY*I_ERI_G3yz_Px_Gx3y_S;
  abcd[567] = I_ERI_H3y2z_Px_Gx3y_S+ABY*I_ERI_G2y2z_Px_Gx3y_S;
  abcd[568] = I_ERI_H2y3z_Px_Gx3y_S+ABY*I_ERI_Gy3z_Px_Gx3y_S;
  abcd[569] = I_ERI_Hy4z_Px_Gx3y_S+ABY*I_ERI_G4z_Px_Gx3y_S;
  abcd[570] = I_ERI_H4xz_Px_Gx3y_S+ABZ*I_ERI_G4x_Px_Gx3y_S;
  abcd[571] = I_ERI_H3xyz_Px_Gx3y_S+ABZ*I_ERI_G3xy_Px_Gx3y_S;
  abcd[572] = I_ERI_H3x2z_Px_Gx3y_S+ABZ*I_ERI_G3xz_Px_Gx3y_S;
  abcd[573] = I_ERI_H2x2yz_Px_Gx3y_S+ABZ*I_ERI_G2x2y_Px_Gx3y_S;
  abcd[574] = I_ERI_H2xy2z_Px_Gx3y_S+ABZ*I_ERI_G2xyz_Px_Gx3y_S;
  abcd[575] = I_ERI_H2x3z_Px_Gx3y_S+ABZ*I_ERI_G2x2z_Px_Gx3y_S;
  abcd[576] = I_ERI_Hx3yz_Px_Gx3y_S+ABZ*I_ERI_Gx3y_Px_Gx3y_S;
  abcd[577] = I_ERI_Hx2y2z_Px_Gx3y_S+ABZ*I_ERI_Gx2yz_Px_Gx3y_S;
  abcd[578] = I_ERI_Hxy3z_Px_Gx3y_S+ABZ*I_ERI_Gxy2z_Px_Gx3y_S;
  abcd[579] = I_ERI_Hx4z_Px_Gx3y_S+ABZ*I_ERI_Gx3z_Px_Gx3y_S;
  abcd[580] = I_ERI_H4yz_Px_Gx3y_S+ABZ*I_ERI_G4y_Px_Gx3y_S;
  abcd[581] = I_ERI_H3y2z_Px_Gx3y_S+ABZ*I_ERI_G3yz_Px_Gx3y_S;
  abcd[582] = I_ERI_H2y3z_Px_Gx3y_S+ABZ*I_ERI_G2y2z_Px_Gx3y_S;
  abcd[583] = I_ERI_Hy4z_Px_Gx3y_S+ABZ*I_ERI_Gy3z_Px_Gx3y_S;
  abcd[584] = I_ERI_H5z_Px_Gx3y_S+ABZ*I_ERI_G4z_Px_Gx3y_S;
  abcd[585] = I_ERI_H4xy_Py_Gx3y_S+ABY*I_ERI_G4x_Py_Gx3y_S;
  abcd[586] = I_ERI_H3x2y_Py_Gx3y_S+ABY*I_ERI_G3xy_Py_Gx3y_S;
  abcd[587] = I_ERI_H3xyz_Py_Gx3y_S+ABY*I_ERI_G3xz_Py_Gx3y_S;
  abcd[588] = I_ERI_H2x3y_Py_Gx3y_S+ABY*I_ERI_G2x2y_Py_Gx3y_S;
  abcd[589] = I_ERI_H2x2yz_Py_Gx3y_S+ABY*I_ERI_G2xyz_Py_Gx3y_S;
  abcd[590] = I_ERI_H2xy2z_Py_Gx3y_S+ABY*I_ERI_G2x2z_Py_Gx3y_S;
  abcd[591] = I_ERI_Hx4y_Py_Gx3y_S+ABY*I_ERI_Gx3y_Py_Gx3y_S;
  abcd[592] = I_ERI_Hx3yz_Py_Gx3y_S+ABY*I_ERI_Gx2yz_Py_Gx3y_S;
  abcd[593] = I_ERI_Hx2y2z_Py_Gx3y_S+ABY*I_ERI_Gxy2z_Py_Gx3y_S;
  abcd[594] = I_ERI_Hxy3z_Py_Gx3y_S+ABY*I_ERI_Gx3z_Py_Gx3y_S;
  abcd[595] = I_ERI_H5y_Py_Gx3y_S+ABY*I_ERI_G4y_Py_Gx3y_S;
  abcd[596] = I_ERI_H4yz_Py_Gx3y_S+ABY*I_ERI_G3yz_Py_Gx3y_S;
  abcd[597] = I_ERI_H3y2z_Py_Gx3y_S+ABY*I_ERI_G2y2z_Py_Gx3y_S;
  abcd[598] = I_ERI_H2y3z_Py_Gx3y_S+ABY*I_ERI_Gy3z_Py_Gx3y_S;
  abcd[599] = I_ERI_Hy4z_Py_Gx3y_S+ABY*I_ERI_G4z_Py_Gx3y_S;
  abcd[600] = I_ERI_H4xz_Py_Gx3y_S+ABZ*I_ERI_G4x_Py_Gx3y_S;
  abcd[601] = I_ERI_H3xyz_Py_Gx3y_S+ABZ*I_ERI_G3xy_Py_Gx3y_S;
  abcd[602] = I_ERI_H3x2z_Py_Gx3y_S+ABZ*I_ERI_G3xz_Py_Gx3y_S;
  abcd[603] = I_ERI_H2x2yz_Py_Gx3y_S+ABZ*I_ERI_G2x2y_Py_Gx3y_S;
  abcd[604] = I_ERI_H2xy2z_Py_Gx3y_S+ABZ*I_ERI_G2xyz_Py_Gx3y_S;
  abcd[605] = I_ERI_H2x3z_Py_Gx3y_S+ABZ*I_ERI_G2x2z_Py_Gx3y_S;
  abcd[606] = I_ERI_Hx3yz_Py_Gx3y_S+ABZ*I_ERI_Gx3y_Py_Gx3y_S;
  abcd[607] = I_ERI_Hx2y2z_Py_Gx3y_S+ABZ*I_ERI_Gx2yz_Py_Gx3y_S;
  abcd[608] = I_ERI_Hxy3z_Py_Gx3y_S+ABZ*I_ERI_Gxy2z_Py_Gx3y_S;
  abcd[609] = I_ERI_Hx4z_Py_Gx3y_S+ABZ*I_ERI_Gx3z_Py_Gx3y_S;
  abcd[610] = I_ERI_H4yz_Py_Gx3y_S+ABZ*I_ERI_G4y_Py_Gx3y_S;
  abcd[611] = I_ERI_H3y2z_Py_Gx3y_S+ABZ*I_ERI_G3yz_Py_Gx3y_S;
  abcd[612] = I_ERI_H2y3z_Py_Gx3y_S+ABZ*I_ERI_G2y2z_Py_Gx3y_S;
  abcd[613] = I_ERI_Hy4z_Py_Gx3y_S+ABZ*I_ERI_Gy3z_Py_Gx3y_S;
  abcd[614] = I_ERI_H5z_Py_Gx3y_S+ABZ*I_ERI_G4z_Py_Gx3y_S;
  abcd[615] = I_ERI_H4xz_Pz_Gx3y_S+ABZ*I_ERI_G4x_Pz_Gx3y_S;
  abcd[616] = I_ERI_H3xyz_Pz_Gx3y_S+ABZ*I_ERI_G3xy_Pz_Gx3y_S;
  abcd[617] = I_ERI_H3x2z_Pz_Gx3y_S+ABZ*I_ERI_G3xz_Pz_Gx3y_S;
  abcd[618] = I_ERI_H2x2yz_Pz_Gx3y_S+ABZ*I_ERI_G2x2y_Pz_Gx3y_S;
  abcd[619] = I_ERI_H2xy2z_Pz_Gx3y_S+ABZ*I_ERI_G2xyz_Pz_Gx3y_S;
  abcd[620] = I_ERI_H2x3z_Pz_Gx3y_S+ABZ*I_ERI_G2x2z_Pz_Gx3y_S;
  abcd[621] = I_ERI_Hx3yz_Pz_Gx3y_S+ABZ*I_ERI_Gx3y_Pz_Gx3y_S;
  abcd[622] = I_ERI_Hx2y2z_Pz_Gx3y_S+ABZ*I_ERI_Gx2yz_Pz_Gx3y_S;
  abcd[623] = I_ERI_Hxy3z_Pz_Gx3y_S+ABZ*I_ERI_Gxy2z_Pz_Gx3y_S;
  abcd[624] = I_ERI_Hx4z_Pz_Gx3y_S+ABZ*I_ERI_Gx3z_Pz_Gx3y_S;
  abcd[625] = I_ERI_H4yz_Pz_Gx3y_S+ABZ*I_ERI_G4y_Pz_Gx3y_S;
  abcd[626] = I_ERI_H3y2z_Pz_Gx3y_S+ABZ*I_ERI_G3yz_Pz_Gx3y_S;
  abcd[627] = I_ERI_H2y3z_Pz_Gx3y_S+ABZ*I_ERI_G2y2z_Pz_Gx3y_S;
  abcd[628] = I_ERI_Hy4z_Pz_Gx3y_S+ABZ*I_ERI_Gy3z_Pz_Gx3y_S;
  abcd[629] = I_ERI_H5z_Pz_Gx3y_S+ABZ*I_ERI_G4z_Pz_Gx3y_S;
  abcd[630] = I_ERI_H5x_Px_Gx2yz_S+ABX*I_ERI_G4x_Px_Gx2yz_S;
  abcd[631] = I_ERI_H4xy_Px_Gx2yz_S+ABX*I_ERI_G3xy_Px_Gx2yz_S;
  abcd[632] = I_ERI_H4xz_Px_Gx2yz_S+ABX*I_ERI_G3xz_Px_Gx2yz_S;
  abcd[633] = I_ERI_H3x2y_Px_Gx2yz_S+ABX*I_ERI_G2x2y_Px_Gx2yz_S;
  abcd[634] = I_ERI_H3xyz_Px_Gx2yz_S+ABX*I_ERI_G2xyz_Px_Gx2yz_S;
  abcd[635] = I_ERI_H3x2z_Px_Gx2yz_S+ABX*I_ERI_G2x2z_Px_Gx2yz_S;
  abcd[636] = I_ERI_H2x3y_Px_Gx2yz_S+ABX*I_ERI_Gx3y_Px_Gx2yz_S;
  abcd[637] = I_ERI_H2x2yz_Px_Gx2yz_S+ABX*I_ERI_Gx2yz_Px_Gx2yz_S;
  abcd[638] = I_ERI_H2xy2z_Px_Gx2yz_S+ABX*I_ERI_Gxy2z_Px_Gx2yz_S;
  abcd[639] = I_ERI_H2x3z_Px_Gx2yz_S+ABX*I_ERI_Gx3z_Px_Gx2yz_S;
  abcd[640] = I_ERI_Hx4y_Px_Gx2yz_S+ABX*I_ERI_G4y_Px_Gx2yz_S;
  abcd[641] = I_ERI_Hx3yz_Px_Gx2yz_S+ABX*I_ERI_G3yz_Px_Gx2yz_S;
  abcd[642] = I_ERI_Hx2y2z_Px_Gx2yz_S+ABX*I_ERI_G2y2z_Px_Gx2yz_S;
  abcd[643] = I_ERI_Hxy3z_Px_Gx2yz_S+ABX*I_ERI_Gy3z_Px_Gx2yz_S;
  abcd[644] = I_ERI_Hx4z_Px_Gx2yz_S+ABX*I_ERI_G4z_Px_Gx2yz_S;
  abcd[645] = I_ERI_H4xy_Px_Gx2yz_S+ABY*I_ERI_G4x_Px_Gx2yz_S;
  abcd[646] = I_ERI_H3x2y_Px_Gx2yz_S+ABY*I_ERI_G3xy_Px_Gx2yz_S;
  abcd[647] = I_ERI_H3xyz_Px_Gx2yz_S+ABY*I_ERI_G3xz_Px_Gx2yz_S;
  abcd[648] = I_ERI_H2x3y_Px_Gx2yz_S+ABY*I_ERI_G2x2y_Px_Gx2yz_S;
  abcd[649] = I_ERI_H2x2yz_Px_Gx2yz_S+ABY*I_ERI_G2xyz_Px_Gx2yz_S;
  abcd[650] = I_ERI_H2xy2z_Px_Gx2yz_S+ABY*I_ERI_G2x2z_Px_Gx2yz_S;
  abcd[651] = I_ERI_Hx4y_Px_Gx2yz_S+ABY*I_ERI_Gx3y_Px_Gx2yz_S;
  abcd[652] = I_ERI_Hx3yz_Px_Gx2yz_S+ABY*I_ERI_Gx2yz_Px_Gx2yz_S;
  abcd[653] = I_ERI_Hx2y2z_Px_Gx2yz_S+ABY*I_ERI_Gxy2z_Px_Gx2yz_S;
  abcd[654] = I_ERI_Hxy3z_Px_Gx2yz_S+ABY*I_ERI_Gx3z_Px_Gx2yz_S;
  abcd[655] = I_ERI_H5y_Px_Gx2yz_S+ABY*I_ERI_G4y_Px_Gx2yz_S;
  abcd[656] = I_ERI_H4yz_Px_Gx2yz_S+ABY*I_ERI_G3yz_Px_Gx2yz_S;
  abcd[657] = I_ERI_H3y2z_Px_Gx2yz_S+ABY*I_ERI_G2y2z_Px_Gx2yz_S;
  abcd[658] = I_ERI_H2y3z_Px_Gx2yz_S+ABY*I_ERI_Gy3z_Px_Gx2yz_S;
  abcd[659] = I_ERI_Hy4z_Px_Gx2yz_S+ABY*I_ERI_G4z_Px_Gx2yz_S;
  abcd[660] = I_ERI_H4xz_Px_Gx2yz_S+ABZ*I_ERI_G4x_Px_Gx2yz_S;
  abcd[661] = I_ERI_H3xyz_Px_Gx2yz_S+ABZ*I_ERI_G3xy_Px_Gx2yz_S;
  abcd[662] = I_ERI_H3x2z_Px_Gx2yz_S+ABZ*I_ERI_G3xz_Px_Gx2yz_S;
  abcd[663] = I_ERI_H2x2yz_Px_Gx2yz_S+ABZ*I_ERI_G2x2y_Px_Gx2yz_S;
  abcd[664] = I_ERI_H2xy2z_Px_Gx2yz_S+ABZ*I_ERI_G2xyz_Px_Gx2yz_S;
  abcd[665] = I_ERI_H2x3z_Px_Gx2yz_S+ABZ*I_ERI_G2x2z_Px_Gx2yz_S;
  abcd[666] = I_ERI_Hx3yz_Px_Gx2yz_S+ABZ*I_ERI_Gx3y_Px_Gx2yz_S;
  abcd[667] = I_ERI_Hx2y2z_Px_Gx2yz_S+ABZ*I_ERI_Gx2yz_Px_Gx2yz_S;
  abcd[668] = I_ERI_Hxy3z_Px_Gx2yz_S+ABZ*I_ERI_Gxy2z_Px_Gx2yz_S;
  abcd[669] = I_ERI_Hx4z_Px_Gx2yz_S+ABZ*I_ERI_Gx3z_Px_Gx2yz_S;
  abcd[670] = I_ERI_H4yz_Px_Gx2yz_S+ABZ*I_ERI_G4y_Px_Gx2yz_S;
  abcd[671] = I_ERI_H3y2z_Px_Gx2yz_S+ABZ*I_ERI_G3yz_Px_Gx2yz_S;
  abcd[672] = I_ERI_H2y3z_Px_Gx2yz_S+ABZ*I_ERI_G2y2z_Px_Gx2yz_S;
  abcd[673] = I_ERI_Hy4z_Px_Gx2yz_S+ABZ*I_ERI_Gy3z_Px_Gx2yz_S;
  abcd[674] = I_ERI_H5z_Px_Gx2yz_S+ABZ*I_ERI_G4z_Px_Gx2yz_S;
  abcd[675] = I_ERI_H4xy_Py_Gx2yz_S+ABY*I_ERI_G4x_Py_Gx2yz_S;
  abcd[676] = I_ERI_H3x2y_Py_Gx2yz_S+ABY*I_ERI_G3xy_Py_Gx2yz_S;
  abcd[677] = I_ERI_H3xyz_Py_Gx2yz_S+ABY*I_ERI_G3xz_Py_Gx2yz_S;
  abcd[678] = I_ERI_H2x3y_Py_Gx2yz_S+ABY*I_ERI_G2x2y_Py_Gx2yz_S;
  abcd[679] = I_ERI_H2x2yz_Py_Gx2yz_S+ABY*I_ERI_G2xyz_Py_Gx2yz_S;
  abcd[680] = I_ERI_H2xy2z_Py_Gx2yz_S+ABY*I_ERI_G2x2z_Py_Gx2yz_S;
  abcd[681] = I_ERI_Hx4y_Py_Gx2yz_S+ABY*I_ERI_Gx3y_Py_Gx2yz_S;
  abcd[682] = I_ERI_Hx3yz_Py_Gx2yz_S+ABY*I_ERI_Gx2yz_Py_Gx2yz_S;
  abcd[683] = I_ERI_Hx2y2z_Py_Gx2yz_S+ABY*I_ERI_Gxy2z_Py_Gx2yz_S;
  abcd[684] = I_ERI_Hxy3z_Py_Gx2yz_S+ABY*I_ERI_Gx3z_Py_Gx2yz_S;
  abcd[685] = I_ERI_H5y_Py_Gx2yz_S+ABY*I_ERI_G4y_Py_Gx2yz_S;
  abcd[686] = I_ERI_H4yz_Py_Gx2yz_S+ABY*I_ERI_G3yz_Py_Gx2yz_S;
  abcd[687] = I_ERI_H3y2z_Py_Gx2yz_S+ABY*I_ERI_G2y2z_Py_Gx2yz_S;
  abcd[688] = I_ERI_H2y3z_Py_Gx2yz_S+ABY*I_ERI_Gy3z_Py_Gx2yz_S;
  abcd[689] = I_ERI_Hy4z_Py_Gx2yz_S+ABY*I_ERI_G4z_Py_Gx2yz_S;
  abcd[690] = I_ERI_H4xz_Py_Gx2yz_S+ABZ*I_ERI_G4x_Py_Gx2yz_S;
  abcd[691] = I_ERI_H3xyz_Py_Gx2yz_S+ABZ*I_ERI_G3xy_Py_Gx2yz_S;
  abcd[692] = I_ERI_H3x2z_Py_Gx2yz_S+ABZ*I_ERI_G3xz_Py_Gx2yz_S;
  abcd[693] = I_ERI_H2x2yz_Py_Gx2yz_S+ABZ*I_ERI_G2x2y_Py_Gx2yz_S;
  abcd[694] = I_ERI_H2xy2z_Py_Gx2yz_S+ABZ*I_ERI_G2xyz_Py_Gx2yz_S;
  abcd[695] = I_ERI_H2x3z_Py_Gx2yz_S+ABZ*I_ERI_G2x2z_Py_Gx2yz_S;
  abcd[696] = I_ERI_Hx3yz_Py_Gx2yz_S+ABZ*I_ERI_Gx3y_Py_Gx2yz_S;
  abcd[697] = I_ERI_Hx2y2z_Py_Gx2yz_S+ABZ*I_ERI_Gx2yz_Py_Gx2yz_S;
  abcd[698] = I_ERI_Hxy3z_Py_Gx2yz_S+ABZ*I_ERI_Gxy2z_Py_Gx2yz_S;
  abcd[699] = I_ERI_Hx4z_Py_Gx2yz_S+ABZ*I_ERI_Gx3z_Py_Gx2yz_S;
  abcd[700] = I_ERI_H4yz_Py_Gx2yz_S+ABZ*I_ERI_G4y_Py_Gx2yz_S;
  abcd[701] = I_ERI_H3y2z_Py_Gx2yz_S+ABZ*I_ERI_G3yz_Py_Gx2yz_S;
  abcd[702] = I_ERI_H2y3z_Py_Gx2yz_S+ABZ*I_ERI_G2y2z_Py_Gx2yz_S;
  abcd[703] = I_ERI_Hy4z_Py_Gx2yz_S+ABZ*I_ERI_Gy3z_Py_Gx2yz_S;
  abcd[704] = I_ERI_H5z_Py_Gx2yz_S+ABZ*I_ERI_G4z_Py_Gx2yz_S;
  abcd[705] = I_ERI_H4xz_Pz_Gx2yz_S+ABZ*I_ERI_G4x_Pz_Gx2yz_S;
  abcd[706] = I_ERI_H3xyz_Pz_Gx2yz_S+ABZ*I_ERI_G3xy_Pz_Gx2yz_S;
  abcd[707] = I_ERI_H3x2z_Pz_Gx2yz_S+ABZ*I_ERI_G3xz_Pz_Gx2yz_S;
  abcd[708] = I_ERI_H2x2yz_Pz_Gx2yz_S+ABZ*I_ERI_G2x2y_Pz_Gx2yz_S;
  abcd[709] = I_ERI_H2xy2z_Pz_Gx2yz_S+ABZ*I_ERI_G2xyz_Pz_Gx2yz_S;
  abcd[710] = I_ERI_H2x3z_Pz_Gx2yz_S+ABZ*I_ERI_G2x2z_Pz_Gx2yz_S;
  abcd[711] = I_ERI_Hx3yz_Pz_Gx2yz_S+ABZ*I_ERI_Gx3y_Pz_Gx2yz_S;
  abcd[712] = I_ERI_Hx2y2z_Pz_Gx2yz_S+ABZ*I_ERI_Gx2yz_Pz_Gx2yz_S;
  abcd[713] = I_ERI_Hxy3z_Pz_Gx2yz_S+ABZ*I_ERI_Gxy2z_Pz_Gx2yz_S;
  abcd[714] = I_ERI_Hx4z_Pz_Gx2yz_S+ABZ*I_ERI_Gx3z_Pz_Gx2yz_S;
  abcd[715] = I_ERI_H4yz_Pz_Gx2yz_S+ABZ*I_ERI_G4y_Pz_Gx2yz_S;
  abcd[716] = I_ERI_H3y2z_Pz_Gx2yz_S+ABZ*I_ERI_G3yz_Pz_Gx2yz_S;
  abcd[717] = I_ERI_H2y3z_Pz_Gx2yz_S+ABZ*I_ERI_G2y2z_Pz_Gx2yz_S;
  abcd[718] = I_ERI_Hy4z_Pz_Gx2yz_S+ABZ*I_ERI_Gy3z_Pz_Gx2yz_S;
  abcd[719] = I_ERI_H5z_Pz_Gx2yz_S+ABZ*I_ERI_G4z_Pz_Gx2yz_S;
  abcd[720] = I_ERI_H5x_Px_Gxy2z_S+ABX*I_ERI_G4x_Px_Gxy2z_S;
  abcd[721] = I_ERI_H4xy_Px_Gxy2z_S+ABX*I_ERI_G3xy_Px_Gxy2z_S;
  abcd[722] = I_ERI_H4xz_Px_Gxy2z_S+ABX*I_ERI_G3xz_Px_Gxy2z_S;
  abcd[723] = I_ERI_H3x2y_Px_Gxy2z_S+ABX*I_ERI_G2x2y_Px_Gxy2z_S;
  abcd[724] = I_ERI_H3xyz_Px_Gxy2z_S+ABX*I_ERI_G2xyz_Px_Gxy2z_S;
  abcd[725] = I_ERI_H3x2z_Px_Gxy2z_S+ABX*I_ERI_G2x2z_Px_Gxy2z_S;
  abcd[726] = I_ERI_H2x3y_Px_Gxy2z_S+ABX*I_ERI_Gx3y_Px_Gxy2z_S;
  abcd[727] = I_ERI_H2x2yz_Px_Gxy2z_S+ABX*I_ERI_Gx2yz_Px_Gxy2z_S;
  abcd[728] = I_ERI_H2xy2z_Px_Gxy2z_S+ABX*I_ERI_Gxy2z_Px_Gxy2z_S;
  abcd[729] = I_ERI_H2x3z_Px_Gxy2z_S+ABX*I_ERI_Gx3z_Px_Gxy2z_S;
  abcd[730] = I_ERI_Hx4y_Px_Gxy2z_S+ABX*I_ERI_G4y_Px_Gxy2z_S;
  abcd[731] = I_ERI_Hx3yz_Px_Gxy2z_S+ABX*I_ERI_G3yz_Px_Gxy2z_S;
  abcd[732] = I_ERI_Hx2y2z_Px_Gxy2z_S+ABX*I_ERI_G2y2z_Px_Gxy2z_S;
  abcd[733] = I_ERI_Hxy3z_Px_Gxy2z_S+ABX*I_ERI_Gy3z_Px_Gxy2z_S;
  abcd[734] = I_ERI_Hx4z_Px_Gxy2z_S+ABX*I_ERI_G4z_Px_Gxy2z_S;
  abcd[735] = I_ERI_H4xy_Px_Gxy2z_S+ABY*I_ERI_G4x_Px_Gxy2z_S;
  abcd[736] = I_ERI_H3x2y_Px_Gxy2z_S+ABY*I_ERI_G3xy_Px_Gxy2z_S;
  abcd[737] = I_ERI_H3xyz_Px_Gxy2z_S+ABY*I_ERI_G3xz_Px_Gxy2z_S;
  abcd[738] = I_ERI_H2x3y_Px_Gxy2z_S+ABY*I_ERI_G2x2y_Px_Gxy2z_S;
  abcd[739] = I_ERI_H2x2yz_Px_Gxy2z_S+ABY*I_ERI_G2xyz_Px_Gxy2z_S;
  abcd[740] = I_ERI_H2xy2z_Px_Gxy2z_S+ABY*I_ERI_G2x2z_Px_Gxy2z_S;
  abcd[741] = I_ERI_Hx4y_Px_Gxy2z_S+ABY*I_ERI_Gx3y_Px_Gxy2z_S;
  abcd[742] = I_ERI_Hx3yz_Px_Gxy2z_S+ABY*I_ERI_Gx2yz_Px_Gxy2z_S;
  abcd[743] = I_ERI_Hx2y2z_Px_Gxy2z_S+ABY*I_ERI_Gxy2z_Px_Gxy2z_S;
  abcd[744] = I_ERI_Hxy3z_Px_Gxy2z_S+ABY*I_ERI_Gx3z_Px_Gxy2z_S;
  abcd[745] = I_ERI_H5y_Px_Gxy2z_S+ABY*I_ERI_G4y_Px_Gxy2z_S;
  abcd[746] = I_ERI_H4yz_Px_Gxy2z_S+ABY*I_ERI_G3yz_Px_Gxy2z_S;
  abcd[747] = I_ERI_H3y2z_Px_Gxy2z_S+ABY*I_ERI_G2y2z_Px_Gxy2z_S;
  abcd[748] = I_ERI_H2y3z_Px_Gxy2z_S+ABY*I_ERI_Gy3z_Px_Gxy2z_S;
  abcd[749] = I_ERI_Hy4z_Px_Gxy2z_S+ABY*I_ERI_G4z_Px_Gxy2z_S;
  abcd[750] = I_ERI_H4xz_Px_Gxy2z_S+ABZ*I_ERI_G4x_Px_Gxy2z_S;
  abcd[751] = I_ERI_H3xyz_Px_Gxy2z_S+ABZ*I_ERI_G3xy_Px_Gxy2z_S;
  abcd[752] = I_ERI_H3x2z_Px_Gxy2z_S+ABZ*I_ERI_G3xz_Px_Gxy2z_S;
  abcd[753] = I_ERI_H2x2yz_Px_Gxy2z_S+ABZ*I_ERI_G2x2y_Px_Gxy2z_S;
  abcd[754] = I_ERI_H2xy2z_Px_Gxy2z_S+ABZ*I_ERI_G2xyz_Px_Gxy2z_S;
  abcd[755] = I_ERI_H2x3z_Px_Gxy2z_S+ABZ*I_ERI_G2x2z_Px_Gxy2z_S;
  abcd[756] = I_ERI_Hx3yz_Px_Gxy2z_S+ABZ*I_ERI_Gx3y_Px_Gxy2z_S;
  abcd[757] = I_ERI_Hx2y2z_Px_Gxy2z_S+ABZ*I_ERI_Gx2yz_Px_Gxy2z_S;
  abcd[758] = I_ERI_Hxy3z_Px_Gxy2z_S+ABZ*I_ERI_Gxy2z_Px_Gxy2z_S;
  abcd[759] = I_ERI_Hx4z_Px_Gxy2z_S+ABZ*I_ERI_Gx3z_Px_Gxy2z_S;
  abcd[760] = I_ERI_H4yz_Px_Gxy2z_S+ABZ*I_ERI_G4y_Px_Gxy2z_S;
  abcd[761] = I_ERI_H3y2z_Px_Gxy2z_S+ABZ*I_ERI_G3yz_Px_Gxy2z_S;
  abcd[762] = I_ERI_H2y3z_Px_Gxy2z_S+ABZ*I_ERI_G2y2z_Px_Gxy2z_S;
  abcd[763] = I_ERI_Hy4z_Px_Gxy2z_S+ABZ*I_ERI_Gy3z_Px_Gxy2z_S;
  abcd[764] = I_ERI_H5z_Px_Gxy2z_S+ABZ*I_ERI_G4z_Px_Gxy2z_S;
  abcd[765] = I_ERI_H4xy_Py_Gxy2z_S+ABY*I_ERI_G4x_Py_Gxy2z_S;
  abcd[766] = I_ERI_H3x2y_Py_Gxy2z_S+ABY*I_ERI_G3xy_Py_Gxy2z_S;
  abcd[767] = I_ERI_H3xyz_Py_Gxy2z_S+ABY*I_ERI_G3xz_Py_Gxy2z_S;
  abcd[768] = I_ERI_H2x3y_Py_Gxy2z_S+ABY*I_ERI_G2x2y_Py_Gxy2z_S;
  abcd[769] = I_ERI_H2x2yz_Py_Gxy2z_S+ABY*I_ERI_G2xyz_Py_Gxy2z_S;
  abcd[770] = I_ERI_H2xy2z_Py_Gxy2z_S+ABY*I_ERI_G2x2z_Py_Gxy2z_S;
  abcd[771] = I_ERI_Hx4y_Py_Gxy2z_S+ABY*I_ERI_Gx3y_Py_Gxy2z_S;
  abcd[772] = I_ERI_Hx3yz_Py_Gxy2z_S+ABY*I_ERI_Gx2yz_Py_Gxy2z_S;
  abcd[773] = I_ERI_Hx2y2z_Py_Gxy2z_S+ABY*I_ERI_Gxy2z_Py_Gxy2z_S;
  abcd[774] = I_ERI_Hxy3z_Py_Gxy2z_S+ABY*I_ERI_Gx3z_Py_Gxy2z_S;
  abcd[775] = I_ERI_H5y_Py_Gxy2z_S+ABY*I_ERI_G4y_Py_Gxy2z_S;
  abcd[776] = I_ERI_H4yz_Py_Gxy2z_S+ABY*I_ERI_G3yz_Py_Gxy2z_S;
  abcd[777] = I_ERI_H3y2z_Py_Gxy2z_S+ABY*I_ERI_G2y2z_Py_Gxy2z_S;
  abcd[778] = I_ERI_H2y3z_Py_Gxy2z_S+ABY*I_ERI_Gy3z_Py_Gxy2z_S;
  abcd[779] = I_ERI_Hy4z_Py_Gxy2z_S+ABY*I_ERI_G4z_Py_Gxy2z_S;
  abcd[780] = I_ERI_H4xz_Py_Gxy2z_S+ABZ*I_ERI_G4x_Py_Gxy2z_S;
  abcd[781] = I_ERI_H3xyz_Py_Gxy2z_S+ABZ*I_ERI_G3xy_Py_Gxy2z_S;
  abcd[782] = I_ERI_H3x2z_Py_Gxy2z_S+ABZ*I_ERI_G3xz_Py_Gxy2z_S;
  abcd[783] = I_ERI_H2x2yz_Py_Gxy2z_S+ABZ*I_ERI_G2x2y_Py_Gxy2z_S;
  abcd[784] = I_ERI_H2xy2z_Py_Gxy2z_S+ABZ*I_ERI_G2xyz_Py_Gxy2z_S;
  abcd[785] = I_ERI_H2x3z_Py_Gxy2z_S+ABZ*I_ERI_G2x2z_Py_Gxy2z_S;
  abcd[786] = I_ERI_Hx3yz_Py_Gxy2z_S+ABZ*I_ERI_Gx3y_Py_Gxy2z_S;
  abcd[787] = I_ERI_Hx2y2z_Py_Gxy2z_S+ABZ*I_ERI_Gx2yz_Py_Gxy2z_S;
  abcd[788] = I_ERI_Hxy3z_Py_Gxy2z_S+ABZ*I_ERI_Gxy2z_Py_Gxy2z_S;
  abcd[789] = I_ERI_Hx4z_Py_Gxy2z_S+ABZ*I_ERI_Gx3z_Py_Gxy2z_S;
  abcd[790] = I_ERI_H4yz_Py_Gxy2z_S+ABZ*I_ERI_G4y_Py_Gxy2z_S;
  abcd[791] = I_ERI_H3y2z_Py_Gxy2z_S+ABZ*I_ERI_G3yz_Py_Gxy2z_S;
  abcd[792] = I_ERI_H2y3z_Py_Gxy2z_S+ABZ*I_ERI_G2y2z_Py_Gxy2z_S;
  abcd[793] = I_ERI_Hy4z_Py_Gxy2z_S+ABZ*I_ERI_Gy3z_Py_Gxy2z_S;
  abcd[794] = I_ERI_H5z_Py_Gxy2z_S+ABZ*I_ERI_G4z_Py_Gxy2z_S;
  abcd[795] = I_ERI_H4xz_Pz_Gxy2z_S+ABZ*I_ERI_G4x_Pz_Gxy2z_S;
  abcd[796] = I_ERI_H3xyz_Pz_Gxy2z_S+ABZ*I_ERI_G3xy_Pz_Gxy2z_S;
  abcd[797] = I_ERI_H3x2z_Pz_Gxy2z_S+ABZ*I_ERI_G3xz_Pz_Gxy2z_S;
  abcd[798] = I_ERI_H2x2yz_Pz_Gxy2z_S+ABZ*I_ERI_G2x2y_Pz_Gxy2z_S;
  abcd[799] = I_ERI_H2xy2z_Pz_Gxy2z_S+ABZ*I_ERI_G2xyz_Pz_Gxy2z_S;
  abcd[800] = I_ERI_H2x3z_Pz_Gxy2z_S+ABZ*I_ERI_G2x2z_Pz_Gxy2z_S;
  abcd[801] = I_ERI_Hx3yz_Pz_Gxy2z_S+ABZ*I_ERI_Gx3y_Pz_Gxy2z_S;
  abcd[802] = I_ERI_Hx2y2z_Pz_Gxy2z_S+ABZ*I_ERI_Gx2yz_Pz_Gxy2z_S;
  abcd[803] = I_ERI_Hxy3z_Pz_Gxy2z_S+ABZ*I_ERI_Gxy2z_Pz_Gxy2z_S;
  abcd[804] = I_ERI_Hx4z_Pz_Gxy2z_S+ABZ*I_ERI_Gx3z_Pz_Gxy2z_S;
  abcd[805] = I_ERI_H4yz_Pz_Gxy2z_S+ABZ*I_ERI_G4y_Pz_Gxy2z_S;
  abcd[806] = I_ERI_H3y2z_Pz_Gxy2z_S+ABZ*I_ERI_G3yz_Pz_Gxy2z_S;
  abcd[807] = I_ERI_H2y3z_Pz_Gxy2z_S+ABZ*I_ERI_G2y2z_Pz_Gxy2z_S;
  abcd[808] = I_ERI_Hy4z_Pz_Gxy2z_S+ABZ*I_ERI_Gy3z_Pz_Gxy2z_S;
  abcd[809] = I_ERI_H5z_Pz_Gxy2z_S+ABZ*I_ERI_G4z_Pz_Gxy2z_S;
  abcd[810] = I_ERI_H5x_Px_Gx3z_S+ABX*I_ERI_G4x_Px_Gx3z_S;
  abcd[811] = I_ERI_H4xy_Px_Gx3z_S+ABX*I_ERI_G3xy_Px_Gx3z_S;
  abcd[812] = I_ERI_H4xz_Px_Gx3z_S+ABX*I_ERI_G3xz_Px_Gx3z_S;
  abcd[813] = I_ERI_H3x2y_Px_Gx3z_S+ABX*I_ERI_G2x2y_Px_Gx3z_S;
  abcd[814] = I_ERI_H3xyz_Px_Gx3z_S+ABX*I_ERI_G2xyz_Px_Gx3z_S;
  abcd[815] = I_ERI_H3x2z_Px_Gx3z_S+ABX*I_ERI_G2x2z_Px_Gx3z_S;
  abcd[816] = I_ERI_H2x3y_Px_Gx3z_S+ABX*I_ERI_Gx3y_Px_Gx3z_S;
  abcd[817] = I_ERI_H2x2yz_Px_Gx3z_S+ABX*I_ERI_Gx2yz_Px_Gx3z_S;
  abcd[818] = I_ERI_H2xy2z_Px_Gx3z_S+ABX*I_ERI_Gxy2z_Px_Gx3z_S;
  abcd[819] = I_ERI_H2x3z_Px_Gx3z_S+ABX*I_ERI_Gx3z_Px_Gx3z_S;
  abcd[820] = I_ERI_Hx4y_Px_Gx3z_S+ABX*I_ERI_G4y_Px_Gx3z_S;
  abcd[821] = I_ERI_Hx3yz_Px_Gx3z_S+ABX*I_ERI_G3yz_Px_Gx3z_S;
  abcd[822] = I_ERI_Hx2y2z_Px_Gx3z_S+ABX*I_ERI_G2y2z_Px_Gx3z_S;
  abcd[823] = I_ERI_Hxy3z_Px_Gx3z_S+ABX*I_ERI_Gy3z_Px_Gx3z_S;
  abcd[824] = I_ERI_Hx4z_Px_Gx3z_S+ABX*I_ERI_G4z_Px_Gx3z_S;
  abcd[825] = I_ERI_H4xy_Px_Gx3z_S+ABY*I_ERI_G4x_Px_Gx3z_S;
  abcd[826] = I_ERI_H3x2y_Px_Gx3z_S+ABY*I_ERI_G3xy_Px_Gx3z_S;
  abcd[827] = I_ERI_H3xyz_Px_Gx3z_S+ABY*I_ERI_G3xz_Px_Gx3z_S;
  abcd[828] = I_ERI_H2x3y_Px_Gx3z_S+ABY*I_ERI_G2x2y_Px_Gx3z_S;
  abcd[829] = I_ERI_H2x2yz_Px_Gx3z_S+ABY*I_ERI_G2xyz_Px_Gx3z_S;
  abcd[830] = I_ERI_H2xy2z_Px_Gx3z_S+ABY*I_ERI_G2x2z_Px_Gx3z_S;
  abcd[831] = I_ERI_Hx4y_Px_Gx3z_S+ABY*I_ERI_Gx3y_Px_Gx3z_S;
  abcd[832] = I_ERI_Hx3yz_Px_Gx3z_S+ABY*I_ERI_Gx2yz_Px_Gx3z_S;
  abcd[833] = I_ERI_Hx2y2z_Px_Gx3z_S+ABY*I_ERI_Gxy2z_Px_Gx3z_S;
  abcd[834] = I_ERI_Hxy3z_Px_Gx3z_S+ABY*I_ERI_Gx3z_Px_Gx3z_S;
  abcd[835] = I_ERI_H5y_Px_Gx3z_S+ABY*I_ERI_G4y_Px_Gx3z_S;
  abcd[836] = I_ERI_H4yz_Px_Gx3z_S+ABY*I_ERI_G3yz_Px_Gx3z_S;
  abcd[837] = I_ERI_H3y2z_Px_Gx3z_S+ABY*I_ERI_G2y2z_Px_Gx3z_S;
  abcd[838] = I_ERI_H2y3z_Px_Gx3z_S+ABY*I_ERI_Gy3z_Px_Gx3z_S;
  abcd[839] = I_ERI_Hy4z_Px_Gx3z_S+ABY*I_ERI_G4z_Px_Gx3z_S;
  abcd[840] = I_ERI_H4xz_Px_Gx3z_S+ABZ*I_ERI_G4x_Px_Gx3z_S;
  abcd[841] = I_ERI_H3xyz_Px_Gx3z_S+ABZ*I_ERI_G3xy_Px_Gx3z_S;
  abcd[842] = I_ERI_H3x2z_Px_Gx3z_S+ABZ*I_ERI_G3xz_Px_Gx3z_S;
  abcd[843] = I_ERI_H2x2yz_Px_Gx3z_S+ABZ*I_ERI_G2x2y_Px_Gx3z_S;
  abcd[844] = I_ERI_H2xy2z_Px_Gx3z_S+ABZ*I_ERI_G2xyz_Px_Gx3z_S;
  abcd[845] = I_ERI_H2x3z_Px_Gx3z_S+ABZ*I_ERI_G2x2z_Px_Gx3z_S;
  abcd[846] = I_ERI_Hx3yz_Px_Gx3z_S+ABZ*I_ERI_Gx3y_Px_Gx3z_S;
  abcd[847] = I_ERI_Hx2y2z_Px_Gx3z_S+ABZ*I_ERI_Gx2yz_Px_Gx3z_S;
  abcd[848] = I_ERI_Hxy3z_Px_Gx3z_S+ABZ*I_ERI_Gxy2z_Px_Gx3z_S;
  abcd[849] = I_ERI_Hx4z_Px_Gx3z_S+ABZ*I_ERI_Gx3z_Px_Gx3z_S;
  abcd[850] = I_ERI_H4yz_Px_Gx3z_S+ABZ*I_ERI_G4y_Px_Gx3z_S;
  abcd[851] = I_ERI_H3y2z_Px_Gx3z_S+ABZ*I_ERI_G3yz_Px_Gx3z_S;
  abcd[852] = I_ERI_H2y3z_Px_Gx3z_S+ABZ*I_ERI_G2y2z_Px_Gx3z_S;
  abcd[853] = I_ERI_Hy4z_Px_Gx3z_S+ABZ*I_ERI_Gy3z_Px_Gx3z_S;
  abcd[854] = I_ERI_H5z_Px_Gx3z_S+ABZ*I_ERI_G4z_Px_Gx3z_S;
  abcd[855] = I_ERI_H4xy_Py_Gx3z_S+ABY*I_ERI_G4x_Py_Gx3z_S;
  abcd[856] = I_ERI_H3x2y_Py_Gx3z_S+ABY*I_ERI_G3xy_Py_Gx3z_S;
  abcd[857] = I_ERI_H3xyz_Py_Gx3z_S+ABY*I_ERI_G3xz_Py_Gx3z_S;
  abcd[858] = I_ERI_H2x3y_Py_Gx3z_S+ABY*I_ERI_G2x2y_Py_Gx3z_S;
  abcd[859] = I_ERI_H2x2yz_Py_Gx3z_S+ABY*I_ERI_G2xyz_Py_Gx3z_S;
  abcd[860] = I_ERI_H2xy2z_Py_Gx3z_S+ABY*I_ERI_G2x2z_Py_Gx3z_S;
  abcd[861] = I_ERI_Hx4y_Py_Gx3z_S+ABY*I_ERI_Gx3y_Py_Gx3z_S;
  abcd[862] = I_ERI_Hx3yz_Py_Gx3z_S+ABY*I_ERI_Gx2yz_Py_Gx3z_S;
  abcd[863] = I_ERI_Hx2y2z_Py_Gx3z_S+ABY*I_ERI_Gxy2z_Py_Gx3z_S;
  abcd[864] = I_ERI_Hxy3z_Py_Gx3z_S+ABY*I_ERI_Gx3z_Py_Gx3z_S;
  abcd[865] = I_ERI_H5y_Py_Gx3z_S+ABY*I_ERI_G4y_Py_Gx3z_S;
  abcd[866] = I_ERI_H4yz_Py_Gx3z_S+ABY*I_ERI_G3yz_Py_Gx3z_S;
  abcd[867] = I_ERI_H3y2z_Py_Gx3z_S+ABY*I_ERI_G2y2z_Py_Gx3z_S;
  abcd[868] = I_ERI_H2y3z_Py_Gx3z_S+ABY*I_ERI_Gy3z_Py_Gx3z_S;
  abcd[869] = I_ERI_Hy4z_Py_Gx3z_S+ABY*I_ERI_G4z_Py_Gx3z_S;
  abcd[870] = I_ERI_H4xz_Py_Gx3z_S+ABZ*I_ERI_G4x_Py_Gx3z_S;
  abcd[871] = I_ERI_H3xyz_Py_Gx3z_S+ABZ*I_ERI_G3xy_Py_Gx3z_S;
  abcd[872] = I_ERI_H3x2z_Py_Gx3z_S+ABZ*I_ERI_G3xz_Py_Gx3z_S;
  abcd[873] = I_ERI_H2x2yz_Py_Gx3z_S+ABZ*I_ERI_G2x2y_Py_Gx3z_S;
  abcd[874] = I_ERI_H2xy2z_Py_Gx3z_S+ABZ*I_ERI_G2xyz_Py_Gx3z_S;
  abcd[875] = I_ERI_H2x3z_Py_Gx3z_S+ABZ*I_ERI_G2x2z_Py_Gx3z_S;
  abcd[876] = I_ERI_Hx3yz_Py_Gx3z_S+ABZ*I_ERI_Gx3y_Py_Gx3z_S;
  abcd[877] = I_ERI_Hx2y2z_Py_Gx3z_S+ABZ*I_ERI_Gx2yz_Py_Gx3z_S;
  abcd[878] = I_ERI_Hxy3z_Py_Gx3z_S+ABZ*I_ERI_Gxy2z_Py_Gx3z_S;
  abcd[879] = I_ERI_Hx4z_Py_Gx3z_S+ABZ*I_ERI_Gx3z_Py_Gx3z_S;
  abcd[880] = I_ERI_H4yz_Py_Gx3z_S+ABZ*I_ERI_G4y_Py_Gx3z_S;
  abcd[881] = I_ERI_H3y2z_Py_Gx3z_S+ABZ*I_ERI_G3yz_Py_Gx3z_S;
  abcd[882] = I_ERI_H2y3z_Py_Gx3z_S+ABZ*I_ERI_G2y2z_Py_Gx3z_S;
  abcd[883] = I_ERI_Hy4z_Py_Gx3z_S+ABZ*I_ERI_Gy3z_Py_Gx3z_S;
  abcd[884] = I_ERI_H5z_Py_Gx3z_S+ABZ*I_ERI_G4z_Py_Gx3z_S;
  abcd[885] = I_ERI_H4xz_Pz_Gx3z_S+ABZ*I_ERI_G4x_Pz_Gx3z_S;
  abcd[886] = I_ERI_H3xyz_Pz_Gx3z_S+ABZ*I_ERI_G3xy_Pz_Gx3z_S;
  abcd[887] = I_ERI_H3x2z_Pz_Gx3z_S+ABZ*I_ERI_G3xz_Pz_Gx3z_S;
  abcd[888] = I_ERI_H2x2yz_Pz_Gx3z_S+ABZ*I_ERI_G2x2y_Pz_Gx3z_S;
  abcd[889] = I_ERI_H2xy2z_Pz_Gx3z_S+ABZ*I_ERI_G2xyz_Pz_Gx3z_S;
  abcd[890] = I_ERI_H2x3z_Pz_Gx3z_S+ABZ*I_ERI_G2x2z_Pz_Gx3z_S;
  abcd[891] = I_ERI_Hx3yz_Pz_Gx3z_S+ABZ*I_ERI_Gx3y_Pz_Gx3z_S;
  abcd[892] = I_ERI_Hx2y2z_Pz_Gx3z_S+ABZ*I_ERI_Gx2yz_Pz_Gx3z_S;
  abcd[893] = I_ERI_Hxy3z_Pz_Gx3z_S+ABZ*I_ERI_Gxy2z_Pz_Gx3z_S;
  abcd[894] = I_ERI_Hx4z_Pz_Gx3z_S+ABZ*I_ERI_Gx3z_Pz_Gx3z_S;
  abcd[895] = I_ERI_H4yz_Pz_Gx3z_S+ABZ*I_ERI_G4y_Pz_Gx3z_S;
  abcd[896] = I_ERI_H3y2z_Pz_Gx3z_S+ABZ*I_ERI_G3yz_Pz_Gx3z_S;
  abcd[897] = I_ERI_H2y3z_Pz_Gx3z_S+ABZ*I_ERI_G2y2z_Pz_Gx3z_S;
  abcd[898] = I_ERI_Hy4z_Pz_Gx3z_S+ABZ*I_ERI_Gy3z_Pz_Gx3z_S;
  abcd[899] = I_ERI_H5z_Pz_Gx3z_S+ABZ*I_ERI_G4z_Pz_Gx3z_S;
  abcd[900] = I_ERI_H5x_Px_G4y_S+ABX*I_ERI_G4x_Px_G4y_S;
  abcd[901] = I_ERI_H4xy_Px_G4y_S+ABX*I_ERI_G3xy_Px_G4y_S;
  abcd[902] = I_ERI_H4xz_Px_G4y_S+ABX*I_ERI_G3xz_Px_G4y_S;
  abcd[903] = I_ERI_H3x2y_Px_G4y_S+ABX*I_ERI_G2x2y_Px_G4y_S;
  abcd[904] = I_ERI_H3xyz_Px_G4y_S+ABX*I_ERI_G2xyz_Px_G4y_S;
  abcd[905] = I_ERI_H3x2z_Px_G4y_S+ABX*I_ERI_G2x2z_Px_G4y_S;
  abcd[906] = I_ERI_H2x3y_Px_G4y_S+ABX*I_ERI_Gx3y_Px_G4y_S;
  abcd[907] = I_ERI_H2x2yz_Px_G4y_S+ABX*I_ERI_Gx2yz_Px_G4y_S;
  abcd[908] = I_ERI_H2xy2z_Px_G4y_S+ABX*I_ERI_Gxy2z_Px_G4y_S;
  abcd[909] = I_ERI_H2x3z_Px_G4y_S+ABX*I_ERI_Gx3z_Px_G4y_S;
  abcd[910] = I_ERI_Hx4y_Px_G4y_S+ABX*I_ERI_G4y_Px_G4y_S;
  abcd[911] = I_ERI_Hx3yz_Px_G4y_S+ABX*I_ERI_G3yz_Px_G4y_S;
  abcd[912] = I_ERI_Hx2y2z_Px_G4y_S+ABX*I_ERI_G2y2z_Px_G4y_S;
  abcd[913] = I_ERI_Hxy3z_Px_G4y_S+ABX*I_ERI_Gy3z_Px_G4y_S;
  abcd[914] = I_ERI_Hx4z_Px_G4y_S+ABX*I_ERI_G4z_Px_G4y_S;
  abcd[915] = I_ERI_H4xy_Px_G4y_S+ABY*I_ERI_G4x_Px_G4y_S;
  abcd[916] = I_ERI_H3x2y_Px_G4y_S+ABY*I_ERI_G3xy_Px_G4y_S;
  abcd[917] = I_ERI_H3xyz_Px_G4y_S+ABY*I_ERI_G3xz_Px_G4y_S;
  abcd[918] = I_ERI_H2x3y_Px_G4y_S+ABY*I_ERI_G2x2y_Px_G4y_S;
  abcd[919] = I_ERI_H2x2yz_Px_G4y_S+ABY*I_ERI_G2xyz_Px_G4y_S;
  abcd[920] = I_ERI_H2xy2z_Px_G4y_S+ABY*I_ERI_G2x2z_Px_G4y_S;
  abcd[921] = I_ERI_Hx4y_Px_G4y_S+ABY*I_ERI_Gx3y_Px_G4y_S;
  abcd[922] = I_ERI_Hx3yz_Px_G4y_S+ABY*I_ERI_Gx2yz_Px_G4y_S;
  abcd[923] = I_ERI_Hx2y2z_Px_G4y_S+ABY*I_ERI_Gxy2z_Px_G4y_S;
  abcd[924] = I_ERI_Hxy3z_Px_G4y_S+ABY*I_ERI_Gx3z_Px_G4y_S;
  abcd[925] = I_ERI_H5y_Px_G4y_S+ABY*I_ERI_G4y_Px_G4y_S;
  abcd[926] = I_ERI_H4yz_Px_G4y_S+ABY*I_ERI_G3yz_Px_G4y_S;
  abcd[927] = I_ERI_H3y2z_Px_G4y_S+ABY*I_ERI_G2y2z_Px_G4y_S;
  abcd[928] = I_ERI_H2y3z_Px_G4y_S+ABY*I_ERI_Gy3z_Px_G4y_S;
  abcd[929] = I_ERI_Hy4z_Px_G4y_S+ABY*I_ERI_G4z_Px_G4y_S;
  abcd[930] = I_ERI_H4xz_Px_G4y_S+ABZ*I_ERI_G4x_Px_G4y_S;
  abcd[931] = I_ERI_H3xyz_Px_G4y_S+ABZ*I_ERI_G3xy_Px_G4y_S;
  abcd[932] = I_ERI_H3x2z_Px_G4y_S+ABZ*I_ERI_G3xz_Px_G4y_S;
  abcd[933] = I_ERI_H2x2yz_Px_G4y_S+ABZ*I_ERI_G2x2y_Px_G4y_S;
  abcd[934] = I_ERI_H2xy2z_Px_G4y_S+ABZ*I_ERI_G2xyz_Px_G4y_S;
  abcd[935] = I_ERI_H2x3z_Px_G4y_S+ABZ*I_ERI_G2x2z_Px_G4y_S;
  abcd[936] = I_ERI_Hx3yz_Px_G4y_S+ABZ*I_ERI_Gx3y_Px_G4y_S;
  abcd[937] = I_ERI_Hx2y2z_Px_G4y_S+ABZ*I_ERI_Gx2yz_Px_G4y_S;
  abcd[938] = I_ERI_Hxy3z_Px_G4y_S+ABZ*I_ERI_Gxy2z_Px_G4y_S;
  abcd[939] = I_ERI_Hx4z_Px_G4y_S+ABZ*I_ERI_Gx3z_Px_G4y_S;
  abcd[940] = I_ERI_H4yz_Px_G4y_S+ABZ*I_ERI_G4y_Px_G4y_S;
  abcd[941] = I_ERI_H3y2z_Px_G4y_S+ABZ*I_ERI_G3yz_Px_G4y_S;
  abcd[942] = I_ERI_H2y3z_Px_G4y_S+ABZ*I_ERI_G2y2z_Px_G4y_S;
  abcd[943] = I_ERI_Hy4z_Px_G4y_S+ABZ*I_ERI_Gy3z_Px_G4y_S;
  abcd[944] = I_ERI_H5z_Px_G4y_S+ABZ*I_ERI_G4z_Px_G4y_S;
  abcd[945] = I_ERI_H4xy_Py_G4y_S+ABY*I_ERI_G4x_Py_G4y_S;
  abcd[946] = I_ERI_H3x2y_Py_G4y_S+ABY*I_ERI_G3xy_Py_G4y_S;
  abcd[947] = I_ERI_H3xyz_Py_G4y_S+ABY*I_ERI_G3xz_Py_G4y_S;
  abcd[948] = I_ERI_H2x3y_Py_G4y_S+ABY*I_ERI_G2x2y_Py_G4y_S;
  abcd[949] = I_ERI_H2x2yz_Py_G4y_S+ABY*I_ERI_G2xyz_Py_G4y_S;
  abcd[950] = I_ERI_H2xy2z_Py_G4y_S+ABY*I_ERI_G2x2z_Py_G4y_S;
  abcd[951] = I_ERI_Hx4y_Py_G4y_S+ABY*I_ERI_Gx3y_Py_G4y_S;
  abcd[952] = I_ERI_Hx3yz_Py_G4y_S+ABY*I_ERI_Gx2yz_Py_G4y_S;
  abcd[953] = I_ERI_Hx2y2z_Py_G4y_S+ABY*I_ERI_Gxy2z_Py_G4y_S;
  abcd[954] = I_ERI_Hxy3z_Py_G4y_S+ABY*I_ERI_Gx3z_Py_G4y_S;
  abcd[955] = I_ERI_H5y_Py_G4y_S+ABY*I_ERI_G4y_Py_G4y_S;
  abcd[956] = I_ERI_H4yz_Py_G4y_S+ABY*I_ERI_G3yz_Py_G4y_S;
  abcd[957] = I_ERI_H3y2z_Py_G4y_S+ABY*I_ERI_G2y2z_Py_G4y_S;
  abcd[958] = I_ERI_H2y3z_Py_G4y_S+ABY*I_ERI_Gy3z_Py_G4y_S;
  abcd[959] = I_ERI_Hy4z_Py_G4y_S+ABY*I_ERI_G4z_Py_G4y_S;
  abcd[960] = I_ERI_H4xz_Py_G4y_S+ABZ*I_ERI_G4x_Py_G4y_S;
  abcd[961] = I_ERI_H3xyz_Py_G4y_S+ABZ*I_ERI_G3xy_Py_G4y_S;
  abcd[962] = I_ERI_H3x2z_Py_G4y_S+ABZ*I_ERI_G3xz_Py_G4y_S;
  abcd[963] = I_ERI_H2x2yz_Py_G4y_S+ABZ*I_ERI_G2x2y_Py_G4y_S;
  abcd[964] = I_ERI_H2xy2z_Py_G4y_S+ABZ*I_ERI_G2xyz_Py_G4y_S;
  abcd[965] = I_ERI_H2x3z_Py_G4y_S+ABZ*I_ERI_G2x2z_Py_G4y_S;
  abcd[966] = I_ERI_Hx3yz_Py_G4y_S+ABZ*I_ERI_Gx3y_Py_G4y_S;
  abcd[967] = I_ERI_Hx2y2z_Py_G4y_S+ABZ*I_ERI_Gx2yz_Py_G4y_S;
  abcd[968] = I_ERI_Hxy3z_Py_G4y_S+ABZ*I_ERI_Gxy2z_Py_G4y_S;
  abcd[969] = I_ERI_Hx4z_Py_G4y_S+ABZ*I_ERI_Gx3z_Py_G4y_S;
  abcd[970] = I_ERI_H4yz_Py_G4y_S+ABZ*I_ERI_G4y_Py_G4y_S;
  abcd[971] = I_ERI_H3y2z_Py_G4y_S+ABZ*I_ERI_G3yz_Py_G4y_S;
  abcd[972] = I_ERI_H2y3z_Py_G4y_S+ABZ*I_ERI_G2y2z_Py_G4y_S;
  abcd[973] = I_ERI_Hy4z_Py_G4y_S+ABZ*I_ERI_Gy3z_Py_G4y_S;
  abcd[974] = I_ERI_H5z_Py_G4y_S+ABZ*I_ERI_G4z_Py_G4y_S;
  abcd[975] = I_ERI_H4xz_Pz_G4y_S+ABZ*I_ERI_G4x_Pz_G4y_S;
  abcd[976] = I_ERI_H3xyz_Pz_G4y_S+ABZ*I_ERI_G3xy_Pz_G4y_S;
  abcd[977] = I_ERI_H3x2z_Pz_G4y_S+ABZ*I_ERI_G3xz_Pz_G4y_S;
  abcd[978] = I_ERI_H2x2yz_Pz_G4y_S+ABZ*I_ERI_G2x2y_Pz_G4y_S;
  abcd[979] = I_ERI_H2xy2z_Pz_G4y_S+ABZ*I_ERI_G2xyz_Pz_G4y_S;
  abcd[980] = I_ERI_H2x3z_Pz_G4y_S+ABZ*I_ERI_G2x2z_Pz_G4y_S;
  abcd[981] = I_ERI_Hx3yz_Pz_G4y_S+ABZ*I_ERI_Gx3y_Pz_G4y_S;
  abcd[982] = I_ERI_Hx2y2z_Pz_G4y_S+ABZ*I_ERI_Gx2yz_Pz_G4y_S;
  abcd[983] = I_ERI_Hxy3z_Pz_G4y_S+ABZ*I_ERI_Gxy2z_Pz_G4y_S;
  abcd[984] = I_ERI_Hx4z_Pz_G4y_S+ABZ*I_ERI_Gx3z_Pz_G4y_S;
  abcd[985] = I_ERI_H4yz_Pz_G4y_S+ABZ*I_ERI_G4y_Pz_G4y_S;
  abcd[986] = I_ERI_H3y2z_Pz_G4y_S+ABZ*I_ERI_G3yz_Pz_G4y_S;
  abcd[987] = I_ERI_H2y3z_Pz_G4y_S+ABZ*I_ERI_G2y2z_Pz_G4y_S;
  abcd[988] = I_ERI_Hy4z_Pz_G4y_S+ABZ*I_ERI_Gy3z_Pz_G4y_S;
  abcd[989] = I_ERI_H5z_Pz_G4y_S+ABZ*I_ERI_G4z_Pz_G4y_S;
  abcd[990] = I_ERI_H5x_Px_G3yz_S+ABX*I_ERI_G4x_Px_G3yz_S;
  abcd[991] = I_ERI_H4xy_Px_G3yz_S+ABX*I_ERI_G3xy_Px_G3yz_S;
  abcd[992] = I_ERI_H4xz_Px_G3yz_S+ABX*I_ERI_G3xz_Px_G3yz_S;
  abcd[993] = I_ERI_H3x2y_Px_G3yz_S+ABX*I_ERI_G2x2y_Px_G3yz_S;
  abcd[994] = I_ERI_H3xyz_Px_G3yz_S+ABX*I_ERI_G2xyz_Px_G3yz_S;
  abcd[995] = I_ERI_H3x2z_Px_G3yz_S+ABX*I_ERI_G2x2z_Px_G3yz_S;
  abcd[996] = I_ERI_H2x3y_Px_G3yz_S+ABX*I_ERI_Gx3y_Px_G3yz_S;
  abcd[997] = I_ERI_H2x2yz_Px_G3yz_S+ABX*I_ERI_Gx2yz_Px_G3yz_S;
  abcd[998] = I_ERI_H2xy2z_Px_G3yz_S+ABX*I_ERI_Gxy2z_Px_G3yz_S;
  abcd[999] = I_ERI_H2x3z_Px_G3yz_S+ABX*I_ERI_Gx3z_Px_G3yz_S;
  abcd[1000] = I_ERI_Hx4y_Px_G3yz_S+ABX*I_ERI_G4y_Px_G3yz_S;
  abcd[1001] = I_ERI_Hx3yz_Px_G3yz_S+ABX*I_ERI_G3yz_Px_G3yz_S;
  abcd[1002] = I_ERI_Hx2y2z_Px_G3yz_S+ABX*I_ERI_G2y2z_Px_G3yz_S;
  abcd[1003] = I_ERI_Hxy3z_Px_G3yz_S+ABX*I_ERI_Gy3z_Px_G3yz_S;
  abcd[1004] = I_ERI_Hx4z_Px_G3yz_S+ABX*I_ERI_G4z_Px_G3yz_S;
  abcd[1005] = I_ERI_H4xy_Px_G3yz_S+ABY*I_ERI_G4x_Px_G3yz_S;
  abcd[1006] = I_ERI_H3x2y_Px_G3yz_S+ABY*I_ERI_G3xy_Px_G3yz_S;
  abcd[1007] = I_ERI_H3xyz_Px_G3yz_S+ABY*I_ERI_G3xz_Px_G3yz_S;
  abcd[1008] = I_ERI_H2x3y_Px_G3yz_S+ABY*I_ERI_G2x2y_Px_G3yz_S;
  abcd[1009] = I_ERI_H2x2yz_Px_G3yz_S+ABY*I_ERI_G2xyz_Px_G3yz_S;
  abcd[1010] = I_ERI_H2xy2z_Px_G3yz_S+ABY*I_ERI_G2x2z_Px_G3yz_S;
  abcd[1011] = I_ERI_Hx4y_Px_G3yz_S+ABY*I_ERI_Gx3y_Px_G3yz_S;
  abcd[1012] = I_ERI_Hx3yz_Px_G3yz_S+ABY*I_ERI_Gx2yz_Px_G3yz_S;
  abcd[1013] = I_ERI_Hx2y2z_Px_G3yz_S+ABY*I_ERI_Gxy2z_Px_G3yz_S;
  abcd[1014] = I_ERI_Hxy3z_Px_G3yz_S+ABY*I_ERI_Gx3z_Px_G3yz_S;
  abcd[1015] = I_ERI_H5y_Px_G3yz_S+ABY*I_ERI_G4y_Px_G3yz_S;
  abcd[1016] = I_ERI_H4yz_Px_G3yz_S+ABY*I_ERI_G3yz_Px_G3yz_S;
  abcd[1017] = I_ERI_H3y2z_Px_G3yz_S+ABY*I_ERI_G2y2z_Px_G3yz_S;
  abcd[1018] = I_ERI_H2y3z_Px_G3yz_S+ABY*I_ERI_Gy3z_Px_G3yz_S;
  abcd[1019] = I_ERI_Hy4z_Px_G3yz_S+ABY*I_ERI_G4z_Px_G3yz_S;
  abcd[1020] = I_ERI_H4xz_Px_G3yz_S+ABZ*I_ERI_G4x_Px_G3yz_S;
  abcd[1021] = I_ERI_H3xyz_Px_G3yz_S+ABZ*I_ERI_G3xy_Px_G3yz_S;
  abcd[1022] = I_ERI_H3x2z_Px_G3yz_S+ABZ*I_ERI_G3xz_Px_G3yz_S;
  abcd[1023] = I_ERI_H2x2yz_Px_G3yz_S+ABZ*I_ERI_G2x2y_Px_G3yz_S;
  abcd[1024] = I_ERI_H2xy2z_Px_G3yz_S+ABZ*I_ERI_G2xyz_Px_G3yz_S;
  abcd[1025] = I_ERI_H2x3z_Px_G3yz_S+ABZ*I_ERI_G2x2z_Px_G3yz_S;
  abcd[1026] = I_ERI_Hx3yz_Px_G3yz_S+ABZ*I_ERI_Gx3y_Px_G3yz_S;
  abcd[1027] = I_ERI_Hx2y2z_Px_G3yz_S+ABZ*I_ERI_Gx2yz_Px_G3yz_S;
  abcd[1028] = I_ERI_Hxy3z_Px_G3yz_S+ABZ*I_ERI_Gxy2z_Px_G3yz_S;
  abcd[1029] = I_ERI_Hx4z_Px_G3yz_S+ABZ*I_ERI_Gx3z_Px_G3yz_S;
  abcd[1030] = I_ERI_H4yz_Px_G3yz_S+ABZ*I_ERI_G4y_Px_G3yz_S;
  abcd[1031] = I_ERI_H3y2z_Px_G3yz_S+ABZ*I_ERI_G3yz_Px_G3yz_S;
  abcd[1032] = I_ERI_H2y3z_Px_G3yz_S+ABZ*I_ERI_G2y2z_Px_G3yz_S;
  abcd[1033] = I_ERI_Hy4z_Px_G3yz_S+ABZ*I_ERI_Gy3z_Px_G3yz_S;
  abcd[1034] = I_ERI_H5z_Px_G3yz_S+ABZ*I_ERI_G4z_Px_G3yz_S;
  abcd[1035] = I_ERI_H4xy_Py_G3yz_S+ABY*I_ERI_G4x_Py_G3yz_S;
  abcd[1036] = I_ERI_H3x2y_Py_G3yz_S+ABY*I_ERI_G3xy_Py_G3yz_S;
  abcd[1037] = I_ERI_H3xyz_Py_G3yz_S+ABY*I_ERI_G3xz_Py_G3yz_S;
  abcd[1038] = I_ERI_H2x3y_Py_G3yz_S+ABY*I_ERI_G2x2y_Py_G3yz_S;
  abcd[1039] = I_ERI_H2x2yz_Py_G3yz_S+ABY*I_ERI_G2xyz_Py_G3yz_S;
  abcd[1040] = I_ERI_H2xy2z_Py_G3yz_S+ABY*I_ERI_G2x2z_Py_G3yz_S;
  abcd[1041] = I_ERI_Hx4y_Py_G3yz_S+ABY*I_ERI_Gx3y_Py_G3yz_S;
  abcd[1042] = I_ERI_Hx3yz_Py_G3yz_S+ABY*I_ERI_Gx2yz_Py_G3yz_S;
  abcd[1043] = I_ERI_Hx2y2z_Py_G3yz_S+ABY*I_ERI_Gxy2z_Py_G3yz_S;
  abcd[1044] = I_ERI_Hxy3z_Py_G3yz_S+ABY*I_ERI_Gx3z_Py_G3yz_S;
  abcd[1045] = I_ERI_H5y_Py_G3yz_S+ABY*I_ERI_G4y_Py_G3yz_S;
  abcd[1046] = I_ERI_H4yz_Py_G3yz_S+ABY*I_ERI_G3yz_Py_G3yz_S;
  abcd[1047] = I_ERI_H3y2z_Py_G3yz_S+ABY*I_ERI_G2y2z_Py_G3yz_S;
  abcd[1048] = I_ERI_H2y3z_Py_G3yz_S+ABY*I_ERI_Gy3z_Py_G3yz_S;
  abcd[1049] = I_ERI_Hy4z_Py_G3yz_S+ABY*I_ERI_G4z_Py_G3yz_S;
  abcd[1050] = I_ERI_H4xz_Py_G3yz_S+ABZ*I_ERI_G4x_Py_G3yz_S;
  abcd[1051] = I_ERI_H3xyz_Py_G3yz_S+ABZ*I_ERI_G3xy_Py_G3yz_S;
  abcd[1052] = I_ERI_H3x2z_Py_G3yz_S+ABZ*I_ERI_G3xz_Py_G3yz_S;
  abcd[1053] = I_ERI_H2x2yz_Py_G3yz_S+ABZ*I_ERI_G2x2y_Py_G3yz_S;
  abcd[1054] = I_ERI_H2xy2z_Py_G3yz_S+ABZ*I_ERI_G2xyz_Py_G3yz_S;
  abcd[1055] = I_ERI_H2x3z_Py_G3yz_S+ABZ*I_ERI_G2x2z_Py_G3yz_S;
  abcd[1056] = I_ERI_Hx3yz_Py_G3yz_S+ABZ*I_ERI_Gx3y_Py_G3yz_S;
  abcd[1057] = I_ERI_Hx2y2z_Py_G3yz_S+ABZ*I_ERI_Gx2yz_Py_G3yz_S;
  abcd[1058] = I_ERI_Hxy3z_Py_G3yz_S+ABZ*I_ERI_Gxy2z_Py_G3yz_S;
  abcd[1059] = I_ERI_Hx4z_Py_G3yz_S+ABZ*I_ERI_Gx3z_Py_G3yz_S;
  abcd[1060] = I_ERI_H4yz_Py_G3yz_S+ABZ*I_ERI_G4y_Py_G3yz_S;
  abcd[1061] = I_ERI_H3y2z_Py_G3yz_S+ABZ*I_ERI_G3yz_Py_G3yz_S;
  abcd[1062] = I_ERI_H2y3z_Py_G3yz_S+ABZ*I_ERI_G2y2z_Py_G3yz_S;
  abcd[1063] = I_ERI_Hy4z_Py_G3yz_S+ABZ*I_ERI_Gy3z_Py_G3yz_S;
  abcd[1064] = I_ERI_H5z_Py_G3yz_S+ABZ*I_ERI_G4z_Py_G3yz_S;
  abcd[1065] = I_ERI_H4xz_Pz_G3yz_S+ABZ*I_ERI_G4x_Pz_G3yz_S;
  abcd[1066] = I_ERI_H3xyz_Pz_G3yz_S+ABZ*I_ERI_G3xy_Pz_G3yz_S;
  abcd[1067] = I_ERI_H3x2z_Pz_G3yz_S+ABZ*I_ERI_G3xz_Pz_G3yz_S;
  abcd[1068] = I_ERI_H2x2yz_Pz_G3yz_S+ABZ*I_ERI_G2x2y_Pz_G3yz_S;
  abcd[1069] = I_ERI_H2xy2z_Pz_G3yz_S+ABZ*I_ERI_G2xyz_Pz_G3yz_S;
  abcd[1070] = I_ERI_H2x3z_Pz_G3yz_S+ABZ*I_ERI_G2x2z_Pz_G3yz_S;
  abcd[1071] = I_ERI_Hx3yz_Pz_G3yz_S+ABZ*I_ERI_Gx3y_Pz_G3yz_S;
  abcd[1072] = I_ERI_Hx2y2z_Pz_G3yz_S+ABZ*I_ERI_Gx2yz_Pz_G3yz_S;
  abcd[1073] = I_ERI_Hxy3z_Pz_G3yz_S+ABZ*I_ERI_Gxy2z_Pz_G3yz_S;
  abcd[1074] = I_ERI_Hx4z_Pz_G3yz_S+ABZ*I_ERI_Gx3z_Pz_G3yz_S;
  abcd[1075] = I_ERI_H4yz_Pz_G3yz_S+ABZ*I_ERI_G4y_Pz_G3yz_S;
  abcd[1076] = I_ERI_H3y2z_Pz_G3yz_S+ABZ*I_ERI_G3yz_Pz_G3yz_S;
  abcd[1077] = I_ERI_H2y3z_Pz_G3yz_S+ABZ*I_ERI_G2y2z_Pz_G3yz_S;
  abcd[1078] = I_ERI_Hy4z_Pz_G3yz_S+ABZ*I_ERI_Gy3z_Pz_G3yz_S;
  abcd[1079] = I_ERI_H5z_Pz_G3yz_S+ABZ*I_ERI_G4z_Pz_G3yz_S;
  abcd[1080] = I_ERI_H5x_Px_G2y2z_S+ABX*I_ERI_G4x_Px_G2y2z_S;
  abcd[1081] = I_ERI_H4xy_Px_G2y2z_S+ABX*I_ERI_G3xy_Px_G2y2z_S;
  abcd[1082] = I_ERI_H4xz_Px_G2y2z_S+ABX*I_ERI_G3xz_Px_G2y2z_S;
  abcd[1083] = I_ERI_H3x2y_Px_G2y2z_S+ABX*I_ERI_G2x2y_Px_G2y2z_S;
  abcd[1084] = I_ERI_H3xyz_Px_G2y2z_S+ABX*I_ERI_G2xyz_Px_G2y2z_S;
  abcd[1085] = I_ERI_H3x2z_Px_G2y2z_S+ABX*I_ERI_G2x2z_Px_G2y2z_S;
  abcd[1086] = I_ERI_H2x3y_Px_G2y2z_S+ABX*I_ERI_Gx3y_Px_G2y2z_S;
  abcd[1087] = I_ERI_H2x2yz_Px_G2y2z_S+ABX*I_ERI_Gx2yz_Px_G2y2z_S;
  abcd[1088] = I_ERI_H2xy2z_Px_G2y2z_S+ABX*I_ERI_Gxy2z_Px_G2y2z_S;
  abcd[1089] = I_ERI_H2x3z_Px_G2y2z_S+ABX*I_ERI_Gx3z_Px_G2y2z_S;
  abcd[1090] = I_ERI_Hx4y_Px_G2y2z_S+ABX*I_ERI_G4y_Px_G2y2z_S;
  abcd[1091] = I_ERI_Hx3yz_Px_G2y2z_S+ABX*I_ERI_G3yz_Px_G2y2z_S;
  abcd[1092] = I_ERI_Hx2y2z_Px_G2y2z_S+ABX*I_ERI_G2y2z_Px_G2y2z_S;
  abcd[1093] = I_ERI_Hxy3z_Px_G2y2z_S+ABX*I_ERI_Gy3z_Px_G2y2z_S;
  abcd[1094] = I_ERI_Hx4z_Px_G2y2z_S+ABX*I_ERI_G4z_Px_G2y2z_S;
  abcd[1095] = I_ERI_H4xy_Px_G2y2z_S+ABY*I_ERI_G4x_Px_G2y2z_S;
  abcd[1096] = I_ERI_H3x2y_Px_G2y2z_S+ABY*I_ERI_G3xy_Px_G2y2z_S;
  abcd[1097] = I_ERI_H3xyz_Px_G2y2z_S+ABY*I_ERI_G3xz_Px_G2y2z_S;
  abcd[1098] = I_ERI_H2x3y_Px_G2y2z_S+ABY*I_ERI_G2x2y_Px_G2y2z_S;
  abcd[1099] = I_ERI_H2x2yz_Px_G2y2z_S+ABY*I_ERI_G2xyz_Px_G2y2z_S;
  abcd[1100] = I_ERI_H2xy2z_Px_G2y2z_S+ABY*I_ERI_G2x2z_Px_G2y2z_S;
  abcd[1101] = I_ERI_Hx4y_Px_G2y2z_S+ABY*I_ERI_Gx3y_Px_G2y2z_S;
  abcd[1102] = I_ERI_Hx3yz_Px_G2y2z_S+ABY*I_ERI_Gx2yz_Px_G2y2z_S;
  abcd[1103] = I_ERI_Hx2y2z_Px_G2y2z_S+ABY*I_ERI_Gxy2z_Px_G2y2z_S;
  abcd[1104] = I_ERI_Hxy3z_Px_G2y2z_S+ABY*I_ERI_Gx3z_Px_G2y2z_S;
  abcd[1105] = I_ERI_H5y_Px_G2y2z_S+ABY*I_ERI_G4y_Px_G2y2z_S;
  abcd[1106] = I_ERI_H4yz_Px_G2y2z_S+ABY*I_ERI_G3yz_Px_G2y2z_S;
  abcd[1107] = I_ERI_H3y2z_Px_G2y2z_S+ABY*I_ERI_G2y2z_Px_G2y2z_S;
  abcd[1108] = I_ERI_H2y3z_Px_G2y2z_S+ABY*I_ERI_Gy3z_Px_G2y2z_S;
  abcd[1109] = I_ERI_Hy4z_Px_G2y2z_S+ABY*I_ERI_G4z_Px_G2y2z_S;
  abcd[1110] = I_ERI_H4xz_Px_G2y2z_S+ABZ*I_ERI_G4x_Px_G2y2z_S;
  abcd[1111] = I_ERI_H3xyz_Px_G2y2z_S+ABZ*I_ERI_G3xy_Px_G2y2z_S;
  abcd[1112] = I_ERI_H3x2z_Px_G2y2z_S+ABZ*I_ERI_G3xz_Px_G2y2z_S;
  abcd[1113] = I_ERI_H2x2yz_Px_G2y2z_S+ABZ*I_ERI_G2x2y_Px_G2y2z_S;
  abcd[1114] = I_ERI_H2xy2z_Px_G2y2z_S+ABZ*I_ERI_G2xyz_Px_G2y2z_S;
  abcd[1115] = I_ERI_H2x3z_Px_G2y2z_S+ABZ*I_ERI_G2x2z_Px_G2y2z_S;
  abcd[1116] = I_ERI_Hx3yz_Px_G2y2z_S+ABZ*I_ERI_Gx3y_Px_G2y2z_S;
  abcd[1117] = I_ERI_Hx2y2z_Px_G2y2z_S+ABZ*I_ERI_Gx2yz_Px_G2y2z_S;
  abcd[1118] = I_ERI_Hxy3z_Px_G2y2z_S+ABZ*I_ERI_Gxy2z_Px_G2y2z_S;
  abcd[1119] = I_ERI_Hx4z_Px_G2y2z_S+ABZ*I_ERI_Gx3z_Px_G2y2z_S;
  abcd[1120] = I_ERI_H4yz_Px_G2y2z_S+ABZ*I_ERI_G4y_Px_G2y2z_S;
  abcd[1121] = I_ERI_H3y2z_Px_G2y2z_S+ABZ*I_ERI_G3yz_Px_G2y2z_S;
  abcd[1122] = I_ERI_H2y3z_Px_G2y2z_S+ABZ*I_ERI_G2y2z_Px_G2y2z_S;
  abcd[1123] = I_ERI_Hy4z_Px_G2y2z_S+ABZ*I_ERI_Gy3z_Px_G2y2z_S;
  abcd[1124] = I_ERI_H5z_Px_G2y2z_S+ABZ*I_ERI_G4z_Px_G2y2z_S;
  abcd[1125] = I_ERI_H4xy_Py_G2y2z_S+ABY*I_ERI_G4x_Py_G2y2z_S;
  abcd[1126] = I_ERI_H3x2y_Py_G2y2z_S+ABY*I_ERI_G3xy_Py_G2y2z_S;
  abcd[1127] = I_ERI_H3xyz_Py_G2y2z_S+ABY*I_ERI_G3xz_Py_G2y2z_S;
  abcd[1128] = I_ERI_H2x3y_Py_G2y2z_S+ABY*I_ERI_G2x2y_Py_G2y2z_S;
  abcd[1129] = I_ERI_H2x2yz_Py_G2y2z_S+ABY*I_ERI_G2xyz_Py_G2y2z_S;
  abcd[1130] = I_ERI_H2xy2z_Py_G2y2z_S+ABY*I_ERI_G2x2z_Py_G2y2z_S;
  abcd[1131] = I_ERI_Hx4y_Py_G2y2z_S+ABY*I_ERI_Gx3y_Py_G2y2z_S;
  abcd[1132] = I_ERI_Hx3yz_Py_G2y2z_S+ABY*I_ERI_Gx2yz_Py_G2y2z_S;
  abcd[1133] = I_ERI_Hx2y2z_Py_G2y2z_S+ABY*I_ERI_Gxy2z_Py_G2y2z_S;
  abcd[1134] = I_ERI_Hxy3z_Py_G2y2z_S+ABY*I_ERI_Gx3z_Py_G2y2z_S;
  abcd[1135] = I_ERI_H5y_Py_G2y2z_S+ABY*I_ERI_G4y_Py_G2y2z_S;
  abcd[1136] = I_ERI_H4yz_Py_G2y2z_S+ABY*I_ERI_G3yz_Py_G2y2z_S;
  abcd[1137] = I_ERI_H3y2z_Py_G2y2z_S+ABY*I_ERI_G2y2z_Py_G2y2z_S;
  abcd[1138] = I_ERI_H2y3z_Py_G2y2z_S+ABY*I_ERI_Gy3z_Py_G2y2z_S;
  abcd[1139] = I_ERI_Hy4z_Py_G2y2z_S+ABY*I_ERI_G4z_Py_G2y2z_S;
  abcd[1140] = I_ERI_H4xz_Py_G2y2z_S+ABZ*I_ERI_G4x_Py_G2y2z_S;
  abcd[1141] = I_ERI_H3xyz_Py_G2y2z_S+ABZ*I_ERI_G3xy_Py_G2y2z_S;
  abcd[1142] = I_ERI_H3x2z_Py_G2y2z_S+ABZ*I_ERI_G3xz_Py_G2y2z_S;
  abcd[1143] = I_ERI_H2x2yz_Py_G2y2z_S+ABZ*I_ERI_G2x2y_Py_G2y2z_S;
  abcd[1144] = I_ERI_H2xy2z_Py_G2y2z_S+ABZ*I_ERI_G2xyz_Py_G2y2z_S;
  abcd[1145] = I_ERI_H2x3z_Py_G2y2z_S+ABZ*I_ERI_G2x2z_Py_G2y2z_S;
  abcd[1146] = I_ERI_Hx3yz_Py_G2y2z_S+ABZ*I_ERI_Gx3y_Py_G2y2z_S;
  abcd[1147] = I_ERI_Hx2y2z_Py_G2y2z_S+ABZ*I_ERI_Gx2yz_Py_G2y2z_S;
  abcd[1148] = I_ERI_Hxy3z_Py_G2y2z_S+ABZ*I_ERI_Gxy2z_Py_G2y2z_S;
  abcd[1149] = I_ERI_Hx4z_Py_G2y2z_S+ABZ*I_ERI_Gx3z_Py_G2y2z_S;
  abcd[1150] = I_ERI_H4yz_Py_G2y2z_S+ABZ*I_ERI_G4y_Py_G2y2z_S;
  abcd[1151] = I_ERI_H3y2z_Py_G2y2z_S+ABZ*I_ERI_G3yz_Py_G2y2z_S;
  abcd[1152] = I_ERI_H2y3z_Py_G2y2z_S+ABZ*I_ERI_G2y2z_Py_G2y2z_S;
  abcd[1153] = I_ERI_Hy4z_Py_G2y2z_S+ABZ*I_ERI_Gy3z_Py_G2y2z_S;
  abcd[1154] = I_ERI_H5z_Py_G2y2z_S+ABZ*I_ERI_G4z_Py_G2y2z_S;
  abcd[1155] = I_ERI_H4xz_Pz_G2y2z_S+ABZ*I_ERI_G4x_Pz_G2y2z_S;
  abcd[1156] = I_ERI_H3xyz_Pz_G2y2z_S+ABZ*I_ERI_G3xy_Pz_G2y2z_S;
  abcd[1157] = I_ERI_H3x2z_Pz_G2y2z_S+ABZ*I_ERI_G3xz_Pz_G2y2z_S;
  abcd[1158] = I_ERI_H2x2yz_Pz_G2y2z_S+ABZ*I_ERI_G2x2y_Pz_G2y2z_S;
  abcd[1159] = I_ERI_H2xy2z_Pz_G2y2z_S+ABZ*I_ERI_G2xyz_Pz_G2y2z_S;
  abcd[1160] = I_ERI_H2x3z_Pz_G2y2z_S+ABZ*I_ERI_G2x2z_Pz_G2y2z_S;
  abcd[1161] = I_ERI_Hx3yz_Pz_G2y2z_S+ABZ*I_ERI_Gx3y_Pz_G2y2z_S;
  abcd[1162] = I_ERI_Hx2y2z_Pz_G2y2z_S+ABZ*I_ERI_Gx2yz_Pz_G2y2z_S;
  abcd[1163] = I_ERI_Hxy3z_Pz_G2y2z_S+ABZ*I_ERI_Gxy2z_Pz_G2y2z_S;
  abcd[1164] = I_ERI_Hx4z_Pz_G2y2z_S+ABZ*I_ERI_Gx3z_Pz_G2y2z_S;
  abcd[1165] = I_ERI_H4yz_Pz_G2y2z_S+ABZ*I_ERI_G4y_Pz_G2y2z_S;
  abcd[1166] = I_ERI_H3y2z_Pz_G2y2z_S+ABZ*I_ERI_G3yz_Pz_G2y2z_S;
  abcd[1167] = I_ERI_H2y3z_Pz_G2y2z_S+ABZ*I_ERI_G2y2z_Pz_G2y2z_S;
  abcd[1168] = I_ERI_Hy4z_Pz_G2y2z_S+ABZ*I_ERI_Gy3z_Pz_G2y2z_S;
  abcd[1169] = I_ERI_H5z_Pz_G2y2z_S+ABZ*I_ERI_G4z_Pz_G2y2z_S;
  abcd[1170] = I_ERI_H5x_Px_Gy3z_S+ABX*I_ERI_G4x_Px_Gy3z_S;
  abcd[1171] = I_ERI_H4xy_Px_Gy3z_S+ABX*I_ERI_G3xy_Px_Gy3z_S;
  abcd[1172] = I_ERI_H4xz_Px_Gy3z_S+ABX*I_ERI_G3xz_Px_Gy3z_S;
  abcd[1173] = I_ERI_H3x2y_Px_Gy3z_S+ABX*I_ERI_G2x2y_Px_Gy3z_S;
  abcd[1174] = I_ERI_H3xyz_Px_Gy3z_S+ABX*I_ERI_G2xyz_Px_Gy3z_S;
  abcd[1175] = I_ERI_H3x2z_Px_Gy3z_S+ABX*I_ERI_G2x2z_Px_Gy3z_S;
  abcd[1176] = I_ERI_H2x3y_Px_Gy3z_S+ABX*I_ERI_Gx3y_Px_Gy3z_S;
  abcd[1177] = I_ERI_H2x2yz_Px_Gy3z_S+ABX*I_ERI_Gx2yz_Px_Gy3z_S;
  abcd[1178] = I_ERI_H2xy2z_Px_Gy3z_S+ABX*I_ERI_Gxy2z_Px_Gy3z_S;
  abcd[1179] = I_ERI_H2x3z_Px_Gy3z_S+ABX*I_ERI_Gx3z_Px_Gy3z_S;
  abcd[1180] = I_ERI_Hx4y_Px_Gy3z_S+ABX*I_ERI_G4y_Px_Gy3z_S;
  abcd[1181] = I_ERI_Hx3yz_Px_Gy3z_S+ABX*I_ERI_G3yz_Px_Gy3z_S;
  abcd[1182] = I_ERI_Hx2y2z_Px_Gy3z_S+ABX*I_ERI_G2y2z_Px_Gy3z_S;
  abcd[1183] = I_ERI_Hxy3z_Px_Gy3z_S+ABX*I_ERI_Gy3z_Px_Gy3z_S;
  abcd[1184] = I_ERI_Hx4z_Px_Gy3z_S+ABX*I_ERI_G4z_Px_Gy3z_S;
  abcd[1185] = I_ERI_H4xy_Px_Gy3z_S+ABY*I_ERI_G4x_Px_Gy3z_S;
  abcd[1186] = I_ERI_H3x2y_Px_Gy3z_S+ABY*I_ERI_G3xy_Px_Gy3z_S;
  abcd[1187] = I_ERI_H3xyz_Px_Gy3z_S+ABY*I_ERI_G3xz_Px_Gy3z_S;
  abcd[1188] = I_ERI_H2x3y_Px_Gy3z_S+ABY*I_ERI_G2x2y_Px_Gy3z_S;
  abcd[1189] = I_ERI_H2x2yz_Px_Gy3z_S+ABY*I_ERI_G2xyz_Px_Gy3z_S;
  abcd[1190] = I_ERI_H2xy2z_Px_Gy3z_S+ABY*I_ERI_G2x2z_Px_Gy3z_S;
  abcd[1191] = I_ERI_Hx4y_Px_Gy3z_S+ABY*I_ERI_Gx3y_Px_Gy3z_S;
  abcd[1192] = I_ERI_Hx3yz_Px_Gy3z_S+ABY*I_ERI_Gx2yz_Px_Gy3z_S;
  abcd[1193] = I_ERI_Hx2y2z_Px_Gy3z_S+ABY*I_ERI_Gxy2z_Px_Gy3z_S;
  abcd[1194] = I_ERI_Hxy3z_Px_Gy3z_S+ABY*I_ERI_Gx3z_Px_Gy3z_S;
  abcd[1195] = I_ERI_H5y_Px_Gy3z_S+ABY*I_ERI_G4y_Px_Gy3z_S;
  abcd[1196] = I_ERI_H4yz_Px_Gy3z_S+ABY*I_ERI_G3yz_Px_Gy3z_S;
  abcd[1197] = I_ERI_H3y2z_Px_Gy3z_S+ABY*I_ERI_G2y2z_Px_Gy3z_S;
  abcd[1198] = I_ERI_H2y3z_Px_Gy3z_S+ABY*I_ERI_Gy3z_Px_Gy3z_S;
  abcd[1199] = I_ERI_Hy4z_Px_Gy3z_S+ABY*I_ERI_G4z_Px_Gy3z_S;
  abcd[1200] = I_ERI_H4xz_Px_Gy3z_S+ABZ*I_ERI_G4x_Px_Gy3z_S;
  abcd[1201] = I_ERI_H3xyz_Px_Gy3z_S+ABZ*I_ERI_G3xy_Px_Gy3z_S;
  abcd[1202] = I_ERI_H3x2z_Px_Gy3z_S+ABZ*I_ERI_G3xz_Px_Gy3z_S;
  abcd[1203] = I_ERI_H2x2yz_Px_Gy3z_S+ABZ*I_ERI_G2x2y_Px_Gy3z_S;
  abcd[1204] = I_ERI_H2xy2z_Px_Gy3z_S+ABZ*I_ERI_G2xyz_Px_Gy3z_S;
  abcd[1205] = I_ERI_H2x3z_Px_Gy3z_S+ABZ*I_ERI_G2x2z_Px_Gy3z_S;
  abcd[1206] = I_ERI_Hx3yz_Px_Gy3z_S+ABZ*I_ERI_Gx3y_Px_Gy3z_S;
  abcd[1207] = I_ERI_Hx2y2z_Px_Gy3z_S+ABZ*I_ERI_Gx2yz_Px_Gy3z_S;
  abcd[1208] = I_ERI_Hxy3z_Px_Gy3z_S+ABZ*I_ERI_Gxy2z_Px_Gy3z_S;
  abcd[1209] = I_ERI_Hx4z_Px_Gy3z_S+ABZ*I_ERI_Gx3z_Px_Gy3z_S;
  abcd[1210] = I_ERI_H4yz_Px_Gy3z_S+ABZ*I_ERI_G4y_Px_Gy3z_S;
  abcd[1211] = I_ERI_H3y2z_Px_Gy3z_S+ABZ*I_ERI_G3yz_Px_Gy3z_S;
  abcd[1212] = I_ERI_H2y3z_Px_Gy3z_S+ABZ*I_ERI_G2y2z_Px_Gy3z_S;
  abcd[1213] = I_ERI_Hy4z_Px_Gy3z_S+ABZ*I_ERI_Gy3z_Px_Gy3z_S;
  abcd[1214] = I_ERI_H5z_Px_Gy3z_S+ABZ*I_ERI_G4z_Px_Gy3z_S;
  abcd[1215] = I_ERI_H4xy_Py_Gy3z_S+ABY*I_ERI_G4x_Py_Gy3z_S;
  abcd[1216] = I_ERI_H3x2y_Py_Gy3z_S+ABY*I_ERI_G3xy_Py_Gy3z_S;
  abcd[1217] = I_ERI_H3xyz_Py_Gy3z_S+ABY*I_ERI_G3xz_Py_Gy3z_S;
  abcd[1218] = I_ERI_H2x3y_Py_Gy3z_S+ABY*I_ERI_G2x2y_Py_Gy3z_S;
  abcd[1219] = I_ERI_H2x2yz_Py_Gy3z_S+ABY*I_ERI_G2xyz_Py_Gy3z_S;
  abcd[1220] = I_ERI_H2xy2z_Py_Gy3z_S+ABY*I_ERI_G2x2z_Py_Gy3z_S;
  abcd[1221] = I_ERI_Hx4y_Py_Gy3z_S+ABY*I_ERI_Gx3y_Py_Gy3z_S;
  abcd[1222] = I_ERI_Hx3yz_Py_Gy3z_S+ABY*I_ERI_Gx2yz_Py_Gy3z_S;
  abcd[1223] = I_ERI_Hx2y2z_Py_Gy3z_S+ABY*I_ERI_Gxy2z_Py_Gy3z_S;
  abcd[1224] = I_ERI_Hxy3z_Py_Gy3z_S+ABY*I_ERI_Gx3z_Py_Gy3z_S;
  abcd[1225] = I_ERI_H5y_Py_Gy3z_S+ABY*I_ERI_G4y_Py_Gy3z_S;
  abcd[1226] = I_ERI_H4yz_Py_Gy3z_S+ABY*I_ERI_G3yz_Py_Gy3z_S;
  abcd[1227] = I_ERI_H3y2z_Py_Gy3z_S+ABY*I_ERI_G2y2z_Py_Gy3z_S;
  abcd[1228] = I_ERI_H2y3z_Py_Gy3z_S+ABY*I_ERI_Gy3z_Py_Gy3z_S;
  abcd[1229] = I_ERI_Hy4z_Py_Gy3z_S+ABY*I_ERI_G4z_Py_Gy3z_S;
  abcd[1230] = I_ERI_H4xz_Py_Gy3z_S+ABZ*I_ERI_G4x_Py_Gy3z_S;
  abcd[1231] = I_ERI_H3xyz_Py_Gy3z_S+ABZ*I_ERI_G3xy_Py_Gy3z_S;
  abcd[1232] = I_ERI_H3x2z_Py_Gy3z_S+ABZ*I_ERI_G3xz_Py_Gy3z_S;
  abcd[1233] = I_ERI_H2x2yz_Py_Gy3z_S+ABZ*I_ERI_G2x2y_Py_Gy3z_S;
  abcd[1234] = I_ERI_H2xy2z_Py_Gy3z_S+ABZ*I_ERI_G2xyz_Py_Gy3z_S;
  abcd[1235] = I_ERI_H2x3z_Py_Gy3z_S+ABZ*I_ERI_G2x2z_Py_Gy3z_S;
  abcd[1236] = I_ERI_Hx3yz_Py_Gy3z_S+ABZ*I_ERI_Gx3y_Py_Gy3z_S;
  abcd[1237] = I_ERI_Hx2y2z_Py_Gy3z_S+ABZ*I_ERI_Gx2yz_Py_Gy3z_S;
  abcd[1238] = I_ERI_Hxy3z_Py_Gy3z_S+ABZ*I_ERI_Gxy2z_Py_Gy3z_S;
  abcd[1239] = I_ERI_Hx4z_Py_Gy3z_S+ABZ*I_ERI_Gx3z_Py_Gy3z_S;
  abcd[1240] = I_ERI_H4yz_Py_Gy3z_S+ABZ*I_ERI_G4y_Py_Gy3z_S;
  abcd[1241] = I_ERI_H3y2z_Py_Gy3z_S+ABZ*I_ERI_G3yz_Py_Gy3z_S;
  abcd[1242] = I_ERI_H2y3z_Py_Gy3z_S+ABZ*I_ERI_G2y2z_Py_Gy3z_S;
  abcd[1243] = I_ERI_Hy4z_Py_Gy3z_S+ABZ*I_ERI_Gy3z_Py_Gy3z_S;
  abcd[1244] = I_ERI_H5z_Py_Gy3z_S+ABZ*I_ERI_G4z_Py_Gy3z_S;
  abcd[1245] = I_ERI_H4xz_Pz_Gy3z_S+ABZ*I_ERI_G4x_Pz_Gy3z_S;
  abcd[1246] = I_ERI_H3xyz_Pz_Gy3z_S+ABZ*I_ERI_G3xy_Pz_Gy3z_S;
  abcd[1247] = I_ERI_H3x2z_Pz_Gy3z_S+ABZ*I_ERI_G3xz_Pz_Gy3z_S;
  abcd[1248] = I_ERI_H2x2yz_Pz_Gy3z_S+ABZ*I_ERI_G2x2y_Pz_Gy3z_S;
  abcd[1249] = I_ERI_H2xy2z_Pz_Gy3z_S+ABZ*I_ERI_G2xyz_Pz_Gy3z_S;
  abcd[1250] = I_ERI_H2x3z_Pz_Gy3z_S+ABZ*I_ERI_G2x2z_Pz_Gy3z_S;
  abcd[1251] = I_ERI_Hx3yz_Pz_Gy3z_S+ABZ*I_ERI_Gx3y_Pz_Gy3z_S;
  abcd[1252] = I_ERI_Hx2y2z_Pz_Gy3z_S+ABZ*I_ERI_Gx2yz_Pz_Gy3z_S;
  abcd[1253] = I_ERI_Hxy3z_Pz_Gy3z_S+ABZ*I_ERI_Gxy2z_Pz_Gy3z_S;
  abcd[1254] = I_ERI_Hx4z_Pz_Gy3z_S+ABZ*I_ERI_Gx3z_Pz_Gy3z_S;
  abcd[1255] = I_ERI_H4yz_Pz_Gy3z_S+ABZ*I_ERI_G4y_Pz_Gy3z_S;
  abcd[1256] = I_ERI_H3y2z_Pz_Gy3z_S+ABZ*I_ERI_G3yz_Pz_Gy3z_S;
  abcd[1257] = I_ERI_H2y3z_Pz_Gy3z_S+ABZ*I_ERI_G2y2z_Pz_Gy3z_S;
  abcd[1258] = I_ERI_Hy4z_Pz_Gy3z_S+ABZ*I_ERI_Gy3z_Pz_Gy3z_S;
  abcd[1259] = I_ERI_H5z_Pz_Gy3z_S+ABZ*I_ERI_G4z_Pz_Gy3z_S;
  abcd[1260] = I_ERI_H5x_Px_G4z_S+ABX*I_ERI_G4x_Px_G4z_S;
  abcd[1261] = I_ERI_H4xy_Px_G4z_S+ABX*I_ERI_G3xy_Px_G4z_S;
  abcd[1262] = I_ERI_H4xz_Px_G4z_S+ABX*I_ERI_G3xz_Px_G4z_S;
  abcd[1263] = I_ERI_H3x2y_Px_G4z_S+ABX*I_ERI_G2x2y_Px_G4z_S;
  abcd[1264] = I_ERI_H3xyz_Px_G4z_S+ABX*I_ERI_G2xyz_Px_G4z_S;
  abcd[1265] = I_ERI_H3x2z_Px_G4z_S+ABX*I_ERI_G2x2z_Px_G4z_S;
  abcd[1266] = I_ERI_H2x3y_Px_G4z_S+ABX*I_ERI_Gx3y_Px_G4z_S;
  abcd[1267] = I_ERI_H2x2yz_Px_G4z_S+ABX*I_ERI_Gx2yz_Px_G4z_S;
  abcd[1268] = I_ERI_H2xy2z_Px_G4z_S+ABX*I_ERI_Gxy2z_Px_G4z_S;
  abcd[1269] = I_ERI_H2x3z_Px_G4z_S+ABX*I_ERI_Gx3z_Px_G4z_S;
  abcd[1270] = I_ERI_Hx4y_Px_G4z_S+ABX*I_ERI_G4y_Px_G4z_S;
  abcd[1271] = I_ERI_Hx3yz_Px_G4z_S+ABX*I_ERI_G3yz_Px_G4z_S;
  abcd[1272] = I_ERI_Hx2y2z_Px_G4z_S+ABX*I_ERI_G2y2z_Px_G4z_S;
  abcd[1273] = I_ERI_Hxy3z_Px_G4z_S+ABX*I_ERI_Gy3z_Px_G4z_S;
  abcd[1274] = I_ERI_Hx4z_Px_G4z_S+ABX*I_ERI_G4z_Px_G4z_S;
  abcd[1275] = I_ERI_H4xy_Px_G4z_S+ABY*I_ERI_G4x_Px_G4z_S;
  abcd[1276] = I_ERI_H3x2y_Px_G4z_S+ABY*I_ERI_G3xy_Px_G4z_S;
  abcd[1277] = I_ERI_H3xyz_Px_G4z_S+ABY*I_ERI_G3xz_Px_G4z_S;
  abcd[1278] = I_ERI_H2x3y_Px_G4z_S+ABY*I_ERI_G2x2y_Px_G4z_S;
  abcd[1279] = I_ERI_H2x2yz_Px_G4z_S+ABY*I_ERI_G2xyz_Px_G4z_S;
  abcd[1280] = I_ERI_H2xy2z_Px_G4z_S+ABY*I_ERI_G2x2z_Px_G4z_S;
  abcd[1281] = I_ERI_Hx4y_Px_G4z_S+ABY*I_ERI_Gx3y_Px_G4z_S;
  abcd[1282] = I_ERI_Hx3yz_Px_G4z_S+ABY*I_ERI_Gx2yz_Px_G4z_S;
  abcd[1283] = I_ERI_Hx2y2z_Px_G4z_S+ABY*I_ERI_Gxy2z_Px_G4z_S;
  abcd[1284] = I_ERI_Hxy3z_Px_G4z_S+ABY*I_ERI_Gx3z_Px_G4z_S;
  abcd[1285] = I_ERI_H5y_Px_G4z_S+ABY*I_ERI_G4y_Px_G4z_S;
  abcd[1286] = I_ERI_H4yz_Px_G4z_S+ABY*I_ERI_G3yz_Px_G4z_S;
  abcd[1287] = I_ERI_H3y2z_Px_G4z_S+ABY*I_ERI_G2y2z_Px_G4z_S;
  abcd[1288] = I_ERI_H2y3z_Px_G4z_S+ABY*I_ERI_Gy3z_Px_G4z_S;
  abcd[1289] = I_ERI_Hy4z_Px_G4z_S+ABY*I_ERI_G4z_Px_G4z_S;
  abcd[1290] = I_ERI_H4xz_Px_G4z_S+ABZ*I_ERI_G4x_Px_G4z_S;
  abcd[1291] = I_ERI_H3xyz_Px_G4z_S+ABZ*I_ERI_G3xy_Px_G4z_S;
  abcd[1292] = I_ERI_H3x2z_Px_G4z_S+ABZ*I_ERI_G3xz_Px_G4z_S;
  abcd[1293] = I_ERI_H2x2yz_Px_G4z_S+ABZ*I_ERI_G2x2y_Px_G4z_S;
  abcd[1294] = I_ERI_H2xy2z_Px_G4z_S+ABZ*I_ERI_G2xyz_Px_G4z_S;
  abcd[1295] = I_ERI_H2x3z_Px_G4z_S+ABZ*I_ERI_G2x2z_Px_G4z_S;
  abcd[1296] = I_ERI_Hx3yz_Px_G4z_S+ABZ*I_ERI_Gx3y_Px_G4z_S;
  abcd[1297] = I_ERI_Hx2y2z_Px_G4z_S+ABZ*I_ERI_Gx2yz_Px_G4z_S;
  abcd[1298] = I_ERI_Hxy3z_Px_G4z_S+ABZ*I_ERI_Gxy2z_Px_G4z_S;
  abcd[1299] = I_ERI_Hx4z_Px_G4z_S+ABZ*I_ERI_Gx3z_Px_G4z_S;
  abcd[1300] = I_ERI_H4yz_Px_G4z_S+ABZ*I_ERI_G4y_Px_G4z_S;
  abcd[1301] = I_ERI_H3y2z_Px_G4z_S+ABZ*I_ERI_G3yz_Px_G4z_S;
  abcd[1302] = I_ERI_H2y3z_Px_G4z_S+ABZ*I_ERI_G2y2z_Px_G4z_S;
  abcd[1303] = I_ERI_Hy4z_Px_G4z_S+ABZ*I_ERI_Gy3z_Px_G4z_S;
  abcd[1304] = I_ERI_H5z_Px_G4z_S+ABZ*I_ERI_G4z_Px_G4z_S;
  abcd[1305] = I_ERI_H4xy_Py_G4z_S+ABY*I_ERI_G4x_Py_G4z_S;
  abcd[1306] = I_ERI_H3x2y_Py_G4z_S+ABY*I_ERI_G3xy_Py_G4z_S;
  abcd[1307] = I_ERI_H3xyz_Py_G4z_S+ABY*I_ERI_G3xz_Py_G4z_S;
  abcd[1308] = I_ERI_H2x3y_Py_G4z_S+ABY*I_ERI_G2x2y_Py_G4z_S;
  abcd[1309] = I_ERI_H2x2yz_Py_G4z_S+ABY*I_ERI_G2xyz_Py_G4z_S;
  abcd[1310] = I_ERI_H2xy2z_Py_G4z_S+ABY*I_ERI_G2x2z_Py_G4z_S;
  abcd[1311] = I_ERI_Hx4y_Py_G4z_S+ABY*I_ERI_Gx3y_Py_G4z_S;
  abcd[1312] = I_ERI_Hx3yz_Py_G4z_S+ABY*I_ERI_Gx2yz_Py_G4z_S;
  abcd[1313] = I_ERI_Hx2y2z_Py_G4z_S+ABY*I_ERI_Gxy2z_Py_G4z_S;
  abcd[1314] = I_ERI_Hxy3z_Py_G4z_S+ABY*I_ERI_Gx3z_Py_G4z_S;
  abcd[1315] = I_ERI_H5y_Py_G4z_S+ABY*I_ERI_G4y_Py_G4z_S;
  abcd[1316] = I_ERI_H4yz_Py_G4z_S+ABY*I_ERI_G3yz_Py_G4z_S;
  abcd[1317] = I_ERI_H3y2z_Py_G4z_S+ABY*I_ERI_G2y2z_Py_G4z_S;
  abcd[1318] = I_ERI_H2y3z_Py_G4z_S+ABY*I_ERI_Gy3z_Py_G4z_S;
  abcd[1319] = I_ERI_Hy4z_Py_G4z_S+ABY*I_ERI_G4z_Py_G4z_S;
  abcd[1320] = I_ERI_H4xz_Py_G4z_S+ABZ*I_ERI_G4x_Py_G4z_S;
  abcd[1321] = I_ERI_H3xyz_Py_G4z_S+ABZ*I_ERI_G3xy_Py_G4z_S;
  abcd[1322] = I_ERI_H3x2z_Py_G4z_S+ABZ*I_ERI_G3xz_Py_G4z_S;
  abcd[1323] = I_ERI_H2x2yz_Py_G4z_S+ABZ*I_ERI_G2x2y_Py_G4z_S;
  abcd[1324] = I_ERI_H2xy2z_Py_G4z_S+ABZ*I_ERI_G2xyz_Py_G4z_S;
  abcd[1325] = I_ERI_H2x3z_Py_G4z_S+ABZ*I_ERI_G2x2z_Py_G4z_S;
  abcd[1326] = I_ERI_Hx3yz_Py_G4z_S+ABZ*I_ERI_Gx3y_Py_G4z_S;
  abcd[1327] = I_ERI_Hx2y2z_Py_G4z_S+ABZ*I_ERI_Gx2yz_Py_G4z_S;
  abcd[1328] = I_ERI_Hxy3z_Py_G4z_S+ABZ*I_ERI_Gxy2z_Py_G4z_S;
  abcd[1329] = I_ERI_Hx4z_Py_G4z_S+ABZ*I_ERI_Gx3z_Py_G4z_S;
  abcd[1330] = I_ERI_H4yz_Py_G4z_S+ABZ*I_ERI_G4y_Py_G4z_S;
  abcd[1331] = I_ERI_H3y2z_Py_G4z_S+ABZ*I_ERI_G3yz_Py_G4z_S;
  abcd[1332] = I_ERI_H2y3z_Py_G4z_S+ABZ*I_ERI_G2y2z_Py_G4z_S;
  abcd[1333] = I_ERI_Hy4z_Py_G4z_S+ABZ*I_ERI_Gy3z_Py_G4z_S;
  abcd[1334] = I_ERI_H5z_Py_G4z_S+ABZ*I_ERI_G4z_Py_G4z_S;
  abcd[1335] = I_ERI_H4xz_Pz_G4z_S+ABZ*I_ERI_G4x_Pz_G4z_S;
  abcd[1336] = I_ERI_H3xyz_Pz_G4z_S+ABZ*I_ERI_G3xy_Pz_G4z_S;
  abcd[1337] = I_ERI_H3x2z_Pz_G4z_S+ABZ*I_ERI_G3xz_Pz_G4z_S;
  abcd[1338] = I_ERI_H2x2yz_Pz_G4z_S+ABZ*I_ERI_G2x2y_Pz_G4z_S;
  abcd[1339] = I_ERI_H2xy2z_Pz_G4z_S+ABZ*I_ERI_G2xyz_Pz_G4z_S;
  abcd[1340] = I_ERI_H2x3z_Pz_G4z_S+ABZ*I_ERI_G2x2z_Pz_G4z_S;
  abcd[1341] = I_ERI_Hx3yz_Pz_G4z_S+ABZ*I_ERI_Gx3y_Pz_G4z_S;
  abcd[1342] = I_ERI_Hx2y2z_Pz_G4z_S+ABZ*I_ERI_Gx2yz_Pz_G4z_S;
  abcd[1343] = I_ERI_Hxy3z_Pz_G4z_S+ABZ*I_ERI_Gxy2z_Pz_G4z_S;
  abcd[1344] = I_ERI_Hx4z_Pz_G4z_S+ABZ*I_ERI_Gx3z_Pz_G4z_S;
  abcd[1345] = I_ERI_H4yz_Pz_G4z_S+ABZ*I_ERI_G4y_Pz_G4z_S;
  abcd[1346] = I_ERI_H3y2z_Pz_G4z_S+ABZ*I_ERI_G3yz_Pz_G4z_S;
  abcd[1347] = I_ERI_H2y3z_Pz_G4z_S+ABZ*I_ERI_G2y2z_Pz_G4z_S;
  abcd[1348] = I_ERI_Hy4z_Pz_G4z_S+ABZ*I_ERI_Gy3z_Pz_G4z_S;
  abcd[1349] = I_ERI_H5z_Pz_G4z_S+ABZ*I_ERI_G4z_Pz_G4z_S;
}
