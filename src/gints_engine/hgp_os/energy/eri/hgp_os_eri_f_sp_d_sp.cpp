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

void hgp_os_eri_f_sp_d_sp(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_G4x_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C2001003 = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C2001003 = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C2001003 = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C2001003 = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C2001003 = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C2001003 = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1002000003 = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1002000003 = 0.0E0;
  Double I_ERI_G4x_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_G4x_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4y_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_G4z_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1002001003 = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1002001003 = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1002001003 = 0.0E0;

  // initialize the significance check for VRR part 
  // this will determine that whether we skip the following part 
  bool isSignificant = false;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double onedz = iexp[ip2];
    Double ic2   = icoe[ip2];
    Double ic2_1 = icoe[ip2+1*inp2];
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
      Double jc2_1 = jcoe[jp2+1*jnp2];
      Double fket  = jfac[jp2];
      Double pref      = fbra*fket;
      Double prefactor = pref;

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
       * shell quartet name: SQ_ERI_F_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M4_vrr = PAX*I_ERI_D2x_S_S_S_M4_vrr+WPX*I_ERI_D2x_S_S_S_M5_vrr+2*oned2z*I_ERI_Px_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_F2xy_S_S_S_M4_vrr = PAY*I_ERI_D2x_S_S_S_M4_vrr+WPY*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_F2xz_S_S_S_M4_vrr = PAZ*I_ERI_D2x_S_S_S_M4_vrr+WPZ*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_Fx2y_S_S_S_M4_vrr = PAX*I_ERI_D2y_S_S_S_M4_vrr+WPX*I_ERI_D2y_S_S_S_M5_vrr;
      Double I_ERI_Fx2z_S_S_S_M4_vrr = PAX*I_ERI_D2z_S_S_S_M4_vrr+WPX*I_ERI_D2z_S_S_S_M5_vrr;
      Double I_ERI_F3y_S_S_S_M4_vrr = PAY*I_ERI_D2y_S_S_S_M4_vrr+WPY*I_ERI_D2y_S_S_S_M5_vrr+2*oned2z*I_ERI_Py_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_F2yz_S_S_S_M4_vrr = PAZ*I_ERI_D2y_S_S_S_M4_vrr+WPZ*I_ERI_D2y_S_S_S_M5_vrr;
      Double I_ERI_F3z_S_S_S_M4_vrr = PAZ*I_ERI_D2z_S_S_S_M4_vrr+WPZ*I_ERI_D2z_S_S_S_M5_vrr+2*oned2z*I_ERI_Pz_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_G_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M3_vrr = PAX*I_ERI_F3x_S_S_S_M3_vrr+WPX*I_ERI_F3x_S_S_S_M4_vrr+3*oned2z*I_ERI_D2x_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G3xy_S_S_S_M3_vrr = PAY*I_ERI_F3x_S_S_S_M3_vrr+WPY*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G3xz_S_S_S_M3_vrr = PAZ*I_ERI_F3x_S_S_S_M3_vrr+WPZ*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G2x2y_S_S_S_M3_vrr = PAY*I_ERI_F2xy_S_S_S_M3_vrr+WPY*I_ERI_F2xy_S_S_S_M4_vrr+oned2z*I_ERI_D2x_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G2xyz_S_S_S_M3_vrr = PAZ*I_ERI_F2xy_S_S_S_M3_vrr+WPZ*I_ERI_F2xy_S_S_S_M4_vrr;
      Double I_ERI_G2x2z_S_S_S_M3_vrr = PAZ*I_ERI_F2xz_S_S_S_M3_vrr+WPZ*I_ERI_F2xz_S_S_S_M4_vrr+oned2z*I_ERI_D2x_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Gx3y_S_S_S_M3_vrr = PAX*I_ERI_F3y_S_S_S_M3_vrr+WPX*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_Gx2yz_S_S_S_M3_vrr = PAZ*I_ERI_Fx2y_S_S_S_M3_vrr+WPZ*I_ERI_Fx2y_S_S_S_M4_vrr;
      Double I_ERI_Gxy2z_S_S_S_M3_vrr = PAY*I_ERI_Fx2z_S_S_S_M3_vrr+WPY*I_ERI_Fx2z_S_S_S_M4_vrr;
      Double I_ERI_Gx3z_S_S_S_M3_vrr = PAX*I_ERI_F3z_S_S_S_M3_vrr+WPX*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4y_S_S_S_M3_vrr = PAY*I_ERI_F3y_S_S_S_M3_vrr+WPY*I_ERI_F3y_S_S_S_M4_vrr+3*oned2z*I_ERI_D2y_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_G3yz_S_S_S_M3_vrr = PAZ*I_ERI_F3y_S_S_S_M3_vrr+WPZ*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_G2y2z_S_S_S_M3_vrr = PAZ*I_ERI_F2yz_S_S_S_M3_vrr+WPZ*I_ERI_F2yz_S_S_S_M4_vrr+oned2z*I_ERI_D2y_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Gy3z_S_S_S_M3_vrr = PAY*I_ERI_F3z_S_S_S_M3_vrr+WPY*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4z_S_S_S_M3_vrr = PAZ*I_ERI_F3z_S_S_S_M3_vrr+WPZ*I_ERI_F3z_S_S_S_M4_vrr+3*oned2z*I_ERI_D2z_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_G_S_P_S_M2
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       ************************************************************/
      Double I_ERI_G4x_S_Px_S_M2_vrr = QCX*I_ERI_G4x_S_S_S_M2_vrr+WQX*I_ERI_G4x_S_S_S_M3_vrr+4*oned2k*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G3xy_S_Px_S_M2_vrr = QCX*I_ERI_G3xy_S_S_S_M2_vrr+WQX*I_ERI_G3xy_S_S_S_M3_vrr+3*oned2k*I_ERI_F2xy_S_S_S_M3_vrr;
      Double I_ERI_G3xz_S_Px_S_M2_vrr = QCX*I_ERI_G3xz_S_S_S_M2_vrr+WQX*I_ERI_G3xz_S_S_S_M3_vrr+3*oned2k*I_ERI_F2xz_S_S_S_M3_vrr;
      Double I_ERI_G2x2y_S_Px_S_M2_vrr = QCX*I_ERI_G2x2y_S_S_S_M2_vrr+WQX*I_ERI_G2x2y_S_S_S_M3_vrr+2*oned2k*I_ERI_Fx2y_S_S_S_M3_vrr;
      Double I_ERI_G2xyz_S_Px_S_M2_vrr = QCX*I_ERI_G2xyz_S_S_S_M2_vrr+WQX*I_ERI_G2xyz_S_S_S_M3_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M3_vrr;
      Double I_ERI_G2x2z_S_Px_S_M2_vrr = QCX*I_ERI_G2x2z_S_S_S_M2_vrr+WQX*I_ERI_G2x2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Fx2z_S_S_S_M3_vrr;
      Double I_ERI_Gx3y_S_Px_S_M2_vrr = QCX*I_ERI_Gx3y_S_S_S_M2_vrr+WQX*I_ERI_Gx3y_S_S_S_M3_vrr+oned2k*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_Gx2yz_S_Px_S_M2_vrr = QCX*I_ERI_Gx2yz_S_S_S_M2_vrr+WQX*I_ERI_Gx2yz_S_S_S_M3_vrr+oned2k*I_ERI_F2yz_S_S_S_M3_vrr;
      Double I_ERI_Gxy2z_S_Px_S_M2_vrr = QCX*I_ERI_Gxy2z_S_S_S_M2_vrr+WQX*I_ERI_Gxy2z_S_S_S_M3_vrr+oned2k*I_ERI_Fy2z_S_S_S_M3_vrr;
      Double I_ERI_Gx3z_S_Px_S_M2_vrr = QCX*I_ERI_Gx3z_S_S_S_M2_vrr+WQX*I_ERI_Gx3z_S_S_S_M3_vrr+oned2k*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4y_S_Px_S_M2_vrr = QCX*I_ERI_G4y_S_S_S_M2_vrr+WQX*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_G3yz_S_Px_S_M2_vrr = QCX*I_ERI_G3yz_S_S_S_M2_vrr+WQX*I_ERI_G3yz_S_S_S_M3_vrr;
      Double I_ERI_G2y2z_S_Px_S_M2_vrr = QCX*I_ERI_G2y2z_S_S_S_M2_vrr+WQX*I_ERI_G2y2z_S_S_S_M3_vrr;
      Double I_ERI_Gy3z_S_Px_S_M2_vrr = QCX*I_ERI_Gy3z_S_S_S_M2_vrr+WQX*I_ERI_Gy3z_S_S_S_M3_vrr;
      Double I_ERI_G4z_S_Px_S_M2_vrr = QCX*I_ERI_G4z_S_S_S_M2_vrr+WQX*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_G4x_S_Py_S_M2_vrr = QCY*I_ERI_G4x_S_S_S_M2_vrr+WQY*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_G3xy_S_Py_S_M2_vrr = QCY*I_ERI_G3xy_S_S_S_M2_vrr+WQY*I_ERI_G3xy_S_S_S_M3_vrr+oned2k*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G3xz_S_Py_S_M2_vrr = QCY*I_ERI_G3xz_S_S_S_M2_vrr+WQY*I_ERI_G3xz_S_S_S_M3_vrr;
      Double I_ERI_G2x2y_S_Py_S_M2_vrr = QCY*I_ERI_G2x2y_S_S_S_M2_vrr+WQY*I_ERI_G2x2y_S_S_S_M3_vrr+2*oned2k*I_ERI_F2xy_S_S_S_M3_vrr;
      Double I_ERI_G2xyz_S_Py_S_M2_vrr = QCY*I_ERI_G2xyz_S_S_S_M2_vrr+WQY*I_ERI_G2xyz_S_S_S_M3_vrr+oned2k*I_ERI_F2xz_S_S_S_M3_vrr;
      Double I_ERI_G2x2z_S_Py_S_M2_vrr = QCY*I_ERI_G2x2z_S_S_S_M2_vrr+WQY*I_ERI_G2x2z_S_S_S_M3_vrr;
      Double I_ERI_Gx3y_S_Py_S_M2_vrr = QCY*I_ERI_Gx3y_S_S_S_M2_vrr+WQY*I_ERI_Gx3y_S_S_S_M3_vrr+3*oned2k*I_ERI_Fx2y_S_S_S_M3_vrr;
      Double I_ERI_Gx2yz_S_Py_S_M2_vrr = QCY*I_ERI_Gx2yz_S_S_S_M2_vrr+WQY*I_ERI_Gx2yz_S_S_S_M3_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M3_vrr;
      Double I_ERI_Gxy2z_S_Py_S_M2_vrr = QCY*I_ERI_Gxy2z_S_S_S_M2_vrr+WQY*I_ERI_Gxy2z_S_S_S_M3_vrr+oned2k*I_ERI_Fx2z_S_S_S_M3_vrr;
      Double I_ERI_Gx3z_S_Py_S_M2_vrr = QCY*I_ERI_Gx3z_S_S_S_M2_vrr+WQY*I_ERI_Gx3z_S_S_S_M3_vrr;
      Double I_ERI_G4y_S_Py_S_M2_vrr = QCY*I_ERI_G4y_S_S_S_M2_vrr+WQY*I_ERI_G4y_S_S_S_M3_vrr+4*oned2k*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_G3yz_S_Py_S_M2_vrr = QCY*I_ERI_G3yz_S_S_S_M2_vrr+WQY*I_ERI_G3yz_S_S_S_M3_vrr+3*oned2k*I_ERI_F2yz_S_S_S_M3_vrr;
      Double I_ERI_G2y2z_S_Py_S_M2_vrr = QCY*I_ERI_G2y2z_S_S_S_M2_vrr+WQY*I_ERI_G2y2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Fy2z_S_S_S_M3_vrr;
      Double I_ERI_Gy3z_S_Py_S_M2_vrr = QCY*I_ERI_Gy3z_S_S_S_M2_vrr+WQY*I_ERI_Gy3z_S_S_S_M3_vrr+oned2k*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_G4z_S_Py_S_M2_vrr = QCY*I_ERI_G4z_S_S_S_M2_vrr+WQY*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_G4x_S_Pz_S_M2_vrr = QCZ*I_ERI_G4x_S_S_S_M2_vrr+WQZ*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_G3xy_S_Pz_S_M2_vrr = QCZ*I_ERI_G3xy_S_S_S_M2_vrr+WQZ*I_ERI_G3xy_S_S_S_M3_vrr;
      Double I_ERI_G3xz_S_Pz_S_M2_vrr = QCZ*I_ERI_G3xz_S_S_S_M2_vrr+WQZ*I_ERI_G3xz_S_S_S_M3_vrr+oned2k*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_G2x2y_S_Pz_S_M2_vrr = QCZ*I_ERI_G2x2y_S_S_S_M2_vrr+WQZ*I_ERI_G2x2y_S_S_S_M3_vrr;
      Double I_ERI_G2xyz_S_Pz_S_M2_vrr = QCZ*I_ERI_G2xyz_S_S_S_M2_vrr+WQZ*I_ERI_G2xyz_S_S_S_M3_vrr+oned2k*I_ERI_F2xy_S_S_S_M3_vrr;
      Double I_ERI_G2x2z_S_Pz_S_M2_vrr = QCZ*I_ERI_G2x2z_S_S_S_M2_vrr+WQZ*I_ERI_G2x2z_S_S_S_M3_vrr+2*oned2k*I_ERI_F2xz_S_S_S_M3_vrr;
      Double I_ERI_Gx3y_S_Pz_S_M2_vrr = QCZ*I_ERI_Gx3y_S_S_S_M2_vrr+WQZ*I_ERI_Gx3y_S_S_S_M3_vrr;
      Double I_ERI_Gx2yz_S_Pz_S_M2_vrr = QCZ*I_ERI_Gx2yz_S_S_S_M2_vrr+WQZ*I_ERI_Gx2yz_S_S_S_M3_vrr+oned2k*I_ERI_Fx2y_S_S_S_M3_vrr;
      Double I_ERI_Gxy2z_S_Pz_S_M2_vrr = QCZ*I_ERI_Gxy2z_S_S_S_M2_vrr+WQZ*I_ERI_Gxy2z_S_S_S_M3_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M3_vrr;
      Double I_ERI_Gx3z_S_Pz_S_M2_vrr = QCZ*I_ERI_Gx3z_S_S_S_M2_vrr+WQZ*I_ERI_Gx3z_S_S_S_M3_vrr+3*oned2k*I_ERI_Fx2z_S_S_S_M3_vrr;
      Double I_ERI_G4y_S_Pz_S_M2_vrr = QCZ*I_ERI_G4y_S_S_S_M2_vrr+WQZ*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_G3yz_S_Pz_S_M2_vrr = QCZ*I_ERI_G3yz_S_S_S_M2_vrr+WQZ*I_ERI_G3yz_S_S_S_M3_vrr+oned2k*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_G2y2z_S_Pz_S_M2_vrr = QCZ*I_ERI_G2y2z_S_S_S_M2_vrr+WQZ*I_ERI_G2y2z_S_S_S_M3_vrr+2*oned2k*I_ERI_F2yz_S_S_S_M3_vrr;
      Double I_ERI_Gy3z_S_Pz_S_M2_vrr = QCZ*I_ERI_Gy3z_S_S_S_M2_vrr+WQZ*I_ERI_Gy3z_S_S_S_M3_vrr+3*oned2k*I_ERI_Fy2z_S_S_S_M3_vrr;
      Double I_ERI_G4z_S_Pz_S_M2_vrr = QCZ*I_ERI_G4z_S_S_S_M2_vrr+WQZ*I_ERI_G4z_S_S_S_M3_vrr+4*oned2k*I_ERI_F3z_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 12 integrals are omitted 
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
      Double I_ERI_D2x_S_D2y_S_M1_vrr = QCY*I_ERI_D2x_S_Py_S_M1_vrr+WQY*I_ERI_D2x_S_Py_S_M2_vrr+oned2e*I_ERI_D2x_S_S_S_M1_vrr-rhod2esq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_D2y_S_M1_vrr = QCY*I_ERI_Dxy_S_Py_S_M1_vrr+WQY*I_ERI_Dxy_S_Py_S_M2_vrr+oned2e*I_ERI_Dxy_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Dxz_S_D2y_S_M1_vrr = QCY*I_ERI_Dxz_S_Py_S_M1_vrr+WQY*I_ERI_Dxz_S_Py_S_M2_vrr+oned2e*I_ERI_Dxz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_D2y_S_M1_vrr = QCY*I_ERI_D2y_S_Py_S_M1_vrr+WQY*I_ERI_D2y_S_Py_S_M2_vrr+oned2e*I_ERI_D2y_S_S_S_M1_vrr-rhod2esq*I_ERI_D2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_Dyz_S_D2y_S_M1_vrr = QCY*I_ERI_Dyz_S_Py_S_M1_vrr+WQY*I_ERI_Dyz_S_Py_S_M2_vrr+oned2e*I_ERI_Dyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_D2z_S_D2y_S_M1_vrr = QCY*I_ERI_D2z_S_Py_S_M1_vrr+WQY*I_ERI_D2z_S_Py_S_M2_vrr+oned2e*I_ERI_D2z_S_S_S_M1_vrr-rhod2esq*I_ERI_D2z_S_S_S_M2_vrr;
      Double I_ERI_D2x_S_D2z_S_M1_vrr = QCZ*I_ERI_D2x_S_Pz_S_M1_vrr+WQZ*I_ERI_D2x_S_Pz_S_M2_vrr+oned2e*I_ERI_D2x_S_S_S_M1_vrr-rhod2esq*I_ERI_D2x_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_D2z_S_M1_vrr = QCZ*I_ERI_Dxy_S_Pz_S_M1_vrr+WQZ*I_ERI_Dxy_S_Pz_S_M2_vrr+oned2e*I_ERI_Dxy_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxy_S_S_S_M2_vrr;
      Double I_ERI_Dxz_S_D2z_S_M1_vrr = QCZ*I_ERI_Dxz_S_Pz_S_M1_vrr+WQZ*I_ERI_Dxz_S_Pz_S_M2_vrr+oned2e*I_ERI_Dxz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dxz_S_S_S_M2_vrr+oned2k*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_D2y_S_D2z_S_M1_vrr = QCZ*I_ERI_D2y_S_Pz_S_M1_vrr+WQZ*I_ERI_D2y_S_Pz_S_M2_vrr+oned2e*I_ERI_D2y_S_S_S_M1_vrr-rhod2esq*I_ERI_D2y_S_S_S_M2_vrr;
      Double I_ERI_Dyz_S_D2z_S_M1_vrr = QCZ*I_ERI_Dyz_S_Pz_S_M1_vrr+WQZ*I_ERI_Dyz_S_Pz_S_M2_vrr+oned2e*I_ERI_Dyz_S_S_S_M1_vrr-rhod2esq*I_ERI_Dyz_S_S_S_M2_vrr+oned2k*I_ERI_Py_S_Pz_S_M2_vrr;
      Double I_ERI_D2z_S_D2z_S_M1_vrr = QCZ*I_ERI_D2z_S_Pz_S_M1_vrr+WQZ*I_ERI_D2z_S_Pz_S_M2_vrr+oned2e*I_ERI_D2z_S_S_S_M1_vrr-rhod2esq*I_ERI_D2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_F_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 20 integrals are omitted 
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
       * shell quartet name: SQ_ERI_G_S_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_Px_S_M1_vrr = PAX*I_ERI_F3x_S_Px_S_M1_vrr+WPX*I_ERI_F3x_S_Px_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Px_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Px_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_Px_S_M1_vrr = PAY*I_ERI_F3x_S_Px_S_M1_vrr+WPY*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G3xz_S_Px_S_M1_vrr = PAZ*I_ERI_F3x_S_Px_S_M1_vrr+WPZ*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G2x2y_S_Px_S_M1_vrr = PAY*I_ERI_F2xy_S_Px_S_M1_vrr+WPY*I_ERI_F2xy_S_Px_S_M2_vrr+oned2z*I_ERI_D2x_S_Px_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_G2xyz_S_Px_S_M1_vrr = PAZ*I_ERI_F2xy_S_Px_S_M1_vrr+WPZ*I_ERI_F2xy_S_Px_S_M2_vrr;
      Double I_ERI_G2x2z_S_Px_S_M1_vrr = PAZ*I_ERI_F2xz_S_Px_S_M1_vrr+WPZ*I_ERI_F2xz_S_Px_S_M2_vrr+oned2z*I_ERI_D2x_S_Px_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_Gx3y_S_Px_S_M1_vrr = PAX*I_ERI_F3y_S_Px_S_M1_vrr+WPX*I_ERI_F3y_S_Px_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Px_S_M1_vrr = PAZ*I_ERI_Fx2y_S_Px_S_M1_vrr+WPZ*I_ERI_Fx2y_S_Px_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Px_S_M1_vrr = PAY*I_ERI_Fx2z_S_Px_S_M1_vrr+WPY*I_ERI_Fx2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3z_S_Px_S_M1_vrr = PAX*I_ERI_F3z_S_Px_S_M1_vrr+WPX*I_ERI_F3z_S_Px_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_Px_S_M1_vrr = PAY*I_ERI_F3y_S_Px_S_M1_vrr+WPY*I_ERI_F3y_S_Px_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Px_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_G3yz_S_Px_S_M1_vrr = PAZ*I_ERI_F3y_S_Px_S_M1_vrr+WPZ*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_G2y2z_S_Px_S_M1_vrr = PAZ*I_ERI_F2yz_S_Px_S_M1_vrr+WPZ*I_ERI_F2yz_S_Px_S_M2_vrr+oned2z*I_ERI_D2y_S_Px_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_Gy3z_S_Px_S_M1_vrr = PAY*I_ERI_F3z_S_Px_S_M1_vrr+WPY*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4z_S_Px_S_M1_vrr = PAZ*I_ERI_F3z_S_Px_S_M1_vrr+WPZ*I_ERI_F3z_S_Px_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Px_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_G4x_S_Py_S_M1_vrr = PAX*I_ERI_F3x_S_Py_S_M1_vrr+WPX*I_ERI_F3x_S_Py_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Py_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_G3xy_S_Py_S_M1_vrr = PAY*I_ERI_F3x_S_Py_S_M1_vrr+WPY*I_ERI_F3x_S_Py_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_Py_S_M1_vrr = PAZ*I_ERI_F3x_S_Py_S_M1_vrr+WPZ*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_G2x2y_S_Py_S_M1_vrr = PAY*I_ERI_F2xy_S_Py_S_M1_vrr+WPY*I_ERI_F2xy_S_Py_S_M2_vrr+oned2z*I_ERI_D2x_S_Py_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Py_S_M2_vrr+oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_Py_S_M1_vrr = PAZ*I_ERI_F2xy_S_Py_S_M1_vrr+WPZ*I_ERI_F2xy_S_Py_S_M2_vrr;
      Double I_ERI_G2x2z_S_Py_S_M1_vrr = PAZ*I_ERI_F2xz_S_Py_S_M1_vrr+WPZ*I_ERI_F2xz_S_Py_S_M2_vrr+oned2z*I_ERI_D2x_S_Py_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Gx3y_S_Py_S_M1_vrr = PAX*I_ERI_F3y_S_Py_S_M1_vrr+WPX*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Py_S_M1_vrr = PAZ*I_ERI_Fx2y_S_Py_S_M1_vrr+WPZ*I_ERI_Fx2y_S_Py_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Py_S_M1_vrr = PAY*I_ERI_Fx2z_S_Py_S_M1_vrr+WPY*I_ERI_Fx2z_S_Py_S_M2_vrr+oned2k*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_Py_S_M1_vrr = PAX*I_ERI_F3z_S_Py_S_M1_vrr+WPX*I_ERI_F3z_S_Py_S_M2_vrr;
      Double I_ERI_G4y_S_Py_S_M1_vrr = PAY*I_ERI_F3y_S_Py_S_M1_vrr+WPY*I_ERI_F3y_S_Py_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Py_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Py_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_Py_S_M1_vrr = PAZ*I_ERI_F3y_S_Py_S_M1_vrr+WPZ*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_G2y2z_S_Py_S_M1_vrr = PAZ*I_ERI_F2yz_S_Py_S_M1_vrr+WPZ*I_ERI_F2yz_S_Py_S_M2_vrr+oned2z*I_ERI_D2y_S_Py_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_Gy3z_S_Py_S_M1_vrr = PAY*I_ERI_F3z_S_Py_S_M1_vrr+WPY*I_ERI_F3z_S_Py_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_Py_S_M1_vrr = PAZ*I_ERI_F3z_S_Py_S_M1_vrr+WPZ*I_ERI_F3z_S_Py_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Py_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_G4x_S_Pz_S_M1_vrr = PAX*I_ERI_F3x_S_Pz_S_M1_vrr+WPX*I_ERI_F3x_S_Pz_S_M2_vrr+3*oned2z*I_ERI_D2x_S_Pz_S_M1_vrr-3*rhod2zsq*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_G3xy_S_Pz_S_M1_vrr = PAY*I_ERI_F3x_S_Pz_S_M1_vrr+WPY*I_ERI_F3x_S_Pz_S_M2_vrr;
      Double I_ERI_G3xz_S_Pz_S_M1_vrr = PAZ*I_ERI_F3x_S_Pz_S_M1_vrr+WPZ*I_ERI_F3x_S_Pz_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_Pz_S_M1_vrr = PAY*I_ERI_F2xy_S_Pz_S_M1_vrr+WPY*I_ERI_F2xy_S_Pz_S_M2_vrr+oned2z*I_ERI_D2x_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_G2xyz_S_Pz_S_M1_vrr = PAZ*I_ERI_F2xy_S_Pz_S_M1_vrr+WPZ*I_ERI_F2xy_S_Pz_S_M2_vrr+oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_Pz_S_M1_vrr = PAZ*I_ERI_F2xz_S_Pz_S_M1_vrr+WPZ*I_ERI_F2xz_S_Pz_S_M2_vrr+oned2z*I_ERI_D2x_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_D2x_S_Pz_S_M2_vrr+oned2k*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_Pz_S_M1_vrr = PAX*I_ERI_F3y_S_Pz_S_M1_vrr+WPX*I_ERI_F3y_S_Pz_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Pz_S_M1_vrr = PAZ*I_ERI_Fx2y_S_Pz_S_M1_vrr+WPZ*I_ERI_Fx2y_S_Pz_S_M2_vrr+oned2k*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Pz_S_M1_vrr = PAY*I_ERI_Fx2z_S_Pz_S_M1_vrr+WPY*I_ERI_Fx2z_S_Pz_S_M2_vrr;
      Double I_ERI_Gx3z_S_Pz_S_M1_vrr = PAX*I_ERI_F3z_S_Pz_S_M1_vrr+WPX*I_ERI_F3z_S_Pz_S_M2_vrr;
      Double I_ERI_G4y_S_Pz_S_M1_vrr = PAY*I_ERI_F3y_S_Pz_S_M1_vrr+WPY*I_ERI_F3y_S_Pz_S_M2_vrr+3*oned2z*I_ERI_D2y_S_Pz_S_M1_vrr-3*rhod2zsq*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_G3yz_S_Pz_S_M1_vrr = PAZ*I_ERI_F3y_S_Pz_S_M1_vrr+WPZ*I_ERI_F3y_S_Pz_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_Pz_S_M1_vrr = PAZ*I_ERI_F2yz_S_Pz_S_M1_vrr+WPZ*I_ERI_F2yz_S_Pz_S_M2_vrr+oned2z*I_ERI_D2y_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_D2y_S_Pz_S_M2_vrr+oned2k*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_Pz_S_M1_vrr = PAY*I_ERI_F3z_S_Pz_S_M1_vrr+WPY*I_ERI_F3z_S_Pz_S_M2_vrr;
      Double I_ERI_G4z_S_Pz_S_M1_vrr = PAZ*I_ERI_F3z_S_Pz_S_M1_vrr+WPZ*I_ERI_F3z_S_Pz_S_M2_vrr+3*oned2z*I_ERI_D2z_S_Pz_S_M1_vrr-3*rhod2zsq*I_ERI_D2z_S_Pz_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_D2x_S_M1_vrr = QCX*I_ERI_G4x_S_Px_S_M1_vrr+WQX*I_ERI_G4x_S_Px_S_M2_vrr+oned2e*I_ERI_G4x_S_S_S_M1_vrr-rhod2esq*I_ERI_G4x_S_S_S_M2_vrr+4*oned2k*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G3xy_S_D2x_S_M1_vrr = QCX*I_ERI_G3xy_S_Px_S_M1_vrr+WQX*I_ERI_G3xy_S_Px_S_M2_vrr+oned2e*I_ERI_G3xy_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xy_S_S_S_M2_vrr+3*oned2k*I_ERI_F2xy_S_Px_S_M2_vrr;
      Double I_ERI_G3xz_S_D2x_S_M1_vrr = QCX*I_ERI_G3xz_S_Px_S_M1_vrr+WQX*I_ERI_G3xz_S_Px_S_M2_vrr+oned2e*I_ERI_G3xz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xz_S_S_S_M2_vrr+3*oned2k*I_ERI_F2xz_S_Px_S_M2_vrr;
      Double I_ERI_G2x2y_S_D2x_S_M1_vrr = QCX*I_ERI_G2x2y_S_Px_S_M1_vrr+WQX*I_ERI_G2x2y_S_Px_S_M2_vrr+oned2e*I_ERI_G2x2y_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Fx2y_S_Px_S_M2_vrr;
      Double I_ERI_G2xyz_S_D2x_S_M1_vrr = QCX*I_ERI_G2xyz_S_Px_S_M1_vrr+WQX*I_ERI_G2xyz_S_Px_S_M2_vrr+oned2e*I_ERI_G2xyz_S_S_S_M1_vrr-rhod2esq*I_ERI_G2xyz_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_Px_S_M2_vrr;
      Double I_ERI_G2x2z_S_D2x_S_M1_vrr = QCX*I_ERI_G2x2z_S_Px_S_M1_vrr+WQX*I_ERI_G2x2z_S_Px_S_M2_vrr+oned2e*I_ERI_G2x2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fx2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3y_S_D2x_S_M1_vrr = QCX*I_ERI_Gx3y_S_Px_S_M1_vrr+WQX*I_ERI_Gx3y_S_Px_S_M2_vrr+oned2e*I_ERI_Gx3y_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3y_S_S_S_M2_vrr+oned2k*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_Gx2yz_S_D2x_S_M1_vrr = QCX*I_ERI_Gx2yz_S_Px_S_M1_vrr+WQX*I_ERI_Gx2yz_S_Px_S_M2_vrr+oned2e*I_ERI_Gx2yz_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx2yz_S_S_S_M2_vrr+oned2k*I_ERI_F2yz_S_Px_S_M2_vrr;
      Double I_ERI_Gxy2z_S_D2x_S_M1_vrr = QCX*I_ERI_Gxy2z_S_Px_S_M1_vrr+WQX*I_ERI_Gxy2z_S_Px_S_M2_vrr+oned2e*I_ERI_Gxy2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gxy2z_S_S_S_M2_vrr+oned2k*I_ERI_Fy2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3z_S_D2x_S_M1_vrr = QCX*I_ERI_Gx3z_S_Px_S_M1_vrr+WQX*I_ERI_Gx3z_S_Px_S_M2_vrr+oned2e*I_ERI_Gx3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3z_S_S_S_M2_vrr+oned2k*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4y_S_D2x_S_M1_vrr = QCX*I_ERI_G4y_S_Px_S_M1_vrr+WQX*I_ERI_G4y_S_Px_S_M2_vrr+oned2e*I_ERI_G4y_S_S_S_M1_vrr-rhod2esq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_D2x_S_M1_vrr = QCX*I_ERI_G3yz_S_Px_S_M1_vrr+WQX*I_ERI_G3yz_S_Px_S_M2_vrr+oned2e*I_ERI_G3yz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_D2x_S_M1_vrr = QCX*I_ERI_G2y2z_S_Px_S_M1_vrr+WQX*I_ERI_G2y2z_S_Px_S_M2_vrr+oned2e*I_ERI_G2y2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_D2x_S_M1_vrr = QCX*I_ERI_Gy3z_S_Px_S_M1_vrr+WQX*I_ERI_Gy3z_S_Px_S_M2_vrr+oned2e*I_ERI_Gy3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_D2x_S_M1_vrr = QCX*I_ERI_G4z_S_Px_S_M1_vrr+WQX*I_ERI_G4z_S_Px_S_M2_vrr+oned2e*I_ERI_G4z_S_S_S_M1_vrr-rhod2esq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_G4x_S_Dxy_S_M1_vrr = QCY*I_ERI_G4x_S_Px_S_M1_vrr+WQY*I_ERI_G4x_S_Px_S_M2_vrr;
      Double I_ERI_G3xy_S_Dxy_S_M1_vrr = QCY*I_ERI_G3xy_S_Px_S_M1_vrr+WQY*I_ERI_G3xy_S_Px_S_M2_vrr+oned2k*I_ERI_F3x_S_Px_S_M2_vrr;
      Double I_ERI_G3xz_S_Dxy_S_M1_vrr = QCY*I_ERI_G3xz_S_Px_S_M1_vrr+WQY*I_ERI_G3xz_S_Px_S_M2_vrr;
      Double I_ERI_G2x2y_S_Dxy_S_M1_vrr = QCY*I_ERI_G2x2y_S_Px_S_M1_vrr+WQY*I_ERI_G2x2y_S_Px_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Px_S_M2_vrr;
      Double I_ERI_G2xyz_S_Dxy_S_M1_vrr = QCY*I_ERI_G2xyz_S_Px_S_M1_vrr+WQY*I_ERI_G2xyz_S_Px_S_M2_vrr+oned2k*I_ERI_F2xz_S_Px_S_M2_vrr;
      Double I_ERI_G2x2z_S_Dxy_S_M1_vrr = QCY*I_ERI_G2x2z_S_Px_S_M1_vrr+WQY*I_ERI_G2x2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3y_S_Dxy_S_M1_vrr = QCY*I_ERI_Gx3y_S_Px_S_M1_vrr+WQY*I_ERI_Gx3y_S_Px_S_M2_vrr+3*oned2k*I_ERI_Fx2y_S_Px_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Dxy_S_M1_vrr = QCY*I_ERI_Gx2yz_S_Px_S_M1_vrr+WQY*I_ERI_Gx2yz_S_Px_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_Px_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Dxy_S_M1_vrr = QCY*I_ERI_Gxy2z_S_Px_S_M1_vrr+WQY*I_ERI_Gxy2z_S_Px_S_M2_vrr+oned2k*I_ERI_Fx2z_S_Px_S_M2_vrr;
      Double I_ERI_Gx3z_S_Dxy_S_M1_vrr = QCY*I_ERI_Gx3z_S_Px_S_M1_vrr+WQY*I_ERI_Gx3z_S_Px_S_M2_vrr;
      Double I_ERI_G4y_S_Dxy_S_M1_vrr = QCY*I_ERI_G4y_S_Px_S_M1_vrr+WQY*I_ERI_G4y_S_Px_S_M2_vrr+4*oned2k*I_ERI_F3y_S_Px_S_M2_vrr;
      Double I_ERI_G3yz_S_Dxy_S_M1_vrr = QCY*I_ERI_G3yz_S_Px_S_M1_vrr+WQY*I_ERI_G3yz_S_Px_S_M2_vrr+3*oned2k*I_ERI_F2yz_S_Px_S_M2_vrr;
      Double I_ERI_G2y2z_S_Dxy_S_M1_vrr = QCY*I_ERI_G2y2z_S_Px_S_M1_vrr+WQY*I_ERI_G2y2z_S_Px_S_M2_vrr+2*oned2k*I_ERI_Fy2z_S_Px_S_M2_vrr;
      Double I_ERI_Gy3z_S_Dxy_S_M1_vrr = QCY*I_ERI_Gy3z_S_Px_S_M1_vrr+WQY*I_ERI_Gy3z_S_Px_S_M2_vrr+oned2k*I_ERI_F3z_S_Px_S_M2_vrr;
      Double I_ERI_G4z_S_Dxy_S_M1_vrr = QCY*I_ERI_G4z_S_Px_S_M1_vrr+WQY*I_ERI_G4z_S_Px_S_M2_vrr;
      Double I_ERI_G4x_S_D2y_S_M1_vrr = QCY*I_ERI_G4x_S_Py_S_M1_vrr+WQY*I_ERI_G4x_S_Py_S_M2_vrr+oned2e*I_ERI_G4x_S_S_S_M1_vrr-rhod2esq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_D2y_S_M1_vrr = QCY*I_ERI_G3xy_S_Py_S_M1_vrr+WQY*I_ERI_G3xy_S_Py_S_M2_vrr+oned2e*I_ERI_G3xy_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xy_S_S_S_M2_vrr+oned2k*I_ERI_F3x_S_Py_S_M2_vrr;
      Double I_ERI_G3xz_S_D2y_S_M1_vrr = QCY*I_ERI_G3xz_S_Py_S_M1_vrr+WQY*I_ERI_G3xz_S_Py_S_M2_vrr+oned2e*I_ERI_G3xz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_D2y_S_M1_vrr = QCY*I_ERI_G2x2y_S_Py_S_M1_vrr+WQY*I_ERI_G2x2y_S_Py_S_M2_vrr+oned2e*I_ERI_G2x2y_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_Py_S_M2_vrr;
      Double I_ERI_G2xyz_S_D2y_S_M1_vrr = QCY*I_ERI_G2xyz_S_Py_S_M1_vrr+WQY*I_ERI_G2xyz_S_Py_S_M2_vrr+oned2e*I_ERI_G2xyz_S_S_S_M1_vrr-rhod2esq*I_ERI_G2xyz_S_S_S_M2_vrr+oned2k*I_ERI_F2xz_S_Py_S_M2_vrr;
      Double I_ERI_G2x2z_S_D2y_S_M1_vrr = QCY*I_ERI_G2x2z_S_Py_S_M1_vrr+WQY*I_ERI_G2x2z_S_Py_S_M2_vrr+oned2e*I_ERI_G2x2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_D2y_S_M1_vrr = QCY*I_ERI_Gx3y_S_Py_S_M1_vrr+WQY*I_ERI_Gx3y_S_Py_S_M2_vrr+oned2e*I_ERI_Gx3y_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3y_S_S_S_M2_vrr+3*oned2k*I_ERI_Fx2y_S_Py_S_M2_vrr;
      Double I_ERI_Gx2yz_S_D2y_S_M1_vrr = QCY*I_ERI_Gx2yz_S_Py_S_M1_vrr+WQY*I_ERI_Gx2yz_S_Py_S_M2_vrr+oned2e*I_ERI_Gx2yz_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_Py_S_M2_vrr;
      Double I_ERI_Gxy2z_S_D2y_S_M1_vrr = QCY*I_ERI_Gxy2z_S_Py_S_M1_vrr+WQY*I_ERI_Gxy2z_S_Py_S_M2_vrr+oned2e*I_ERI_Gxy2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gxy2z_S_S_S_M2_vrr+oned2k*I_ERI_Fx2z_S_Py_S_M2_vrr;
      Double I_ERI_Gx3z_S_D2y_S_M1_vrr = QCY*I_ERI_Gx3z_S_Py_S_M1_vrr+WQY*I_ERI_Gx3z_S_Py_S_M2_vrr+oned2e*I_ERI_Gx3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_D2y_S_M1_vrr = QCY*I_ERI_G4y_S_Py_S_M1_vrr+WQY*I_ERI_G4y_S_Py_S_M2_vrr+oned2e*I_ERI_G4y_S_S_S_M1_vrr-rhod2esq*I_ERI_G4y_S_S_S_M2_vrr+4*oned2k*I_ERI_F3y_S_Py_S_M2_vrr;
      Double I_ERI_G3yz_S_D2y_S_M1_vrr = QCY*I_ERI_G3yz_S_Py_S_M1_vrr+WQY*I_ERI_G3yz_S_Py_S_M2_vrr+oned2e*I_ERI_G3yz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3yz_S_S_S_M2_vrr+3*oned2k*I_ERI_F2yz_S_Py_S_M2_vrr;
      Double I_ERI_G2y2z_S_D2y_S_M1_vrr = QCY*I_ERI_G2y2z_S_Py_S_M1_vrr+WQY*I_ERI_G2y2z_S_Py_S_M2_vrr+oned2e*I_ERI_G2y2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fy2z_S_Py_S_M2_vrr;
      Double I_ERI_Gy3z_S_D2y_S_M1_vrr = QCY*I_ERI_Gy3z_S_Py_S_M1_vrr+WQY*I_ERI_Gy3z_S_Py_S_M2_vrr+oned2e*I_ERI_Gy3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gy3z_S_S_S_M2_vrr+oned2k*I_ERI_F3z_S_Py_S_M2_vrr;
      Double I_ERI_G4z_S_D2y_S_M1_vrr = QCY*I_ERI_G4z_S_Py_S_M1_vrr+WQY*I_ERI_G4z_S_Py_S_M2_vrr+oned2e*I_ERI_G4z_S_S_S_M1_vrr-rhod2esq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_G4x_S_D2z_S_M1_vrr = QCZ*I_ERI_G4x_S_Pz_S_M1_vrr+WQZ*I_ERI_G4x_S_Pz_S_M2_vrr+oned2e*I_ERI_G4x_S_S_S_M1_vrr-rhod2esq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_D2z_S_M1_vrr = QCZ*I_ERI_G3xy_S_Pz_S_M1_vrr+WQZ*I_ERI_G3xy_S_Pz_S_M2_vrr+oned2e*I_ERI_G3xy_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_D2z_S_M1_vrr = QCZ*I_ERI_G3xz_S_Pz_S_M1_vrr+WQZ*I_ERI_G3xz_S_Pz_S_M2_vrr+oned2e*I_ERI_G3xz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3xz_S_S_S_M2_vrr+oned2k*I_ERI_F3x_S_Pz_S_M2_vrr;
      Double I_ERI_G2x2y_S_D2z_S_M1_vrr = QCZ*I_ERI_G2x2y_S_Pz_S_M1_vrr+WQZ*I_ERI_G2x2y_S_Pz_S_M2_vrr+oned2e*I_ERI_G2x2y_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_D2z_S_M1_vrr = QCZ*I_ERI_G2xyz_S_Pz_S_M1_vrr+WQZ*I_ERI_G2xyz_S_Pz_S_M2_vrr+oned2e*I_ERI_G2xyz_S_S_S_M1_vrr-rhod2esq*I_ERI_G2xyz_S_S_S_M2_vrr+oned2k*I_ERI_F2xy_S_Pz_S_M2_vrr;
      Double I_ERI_G2x2z_S_D2z_S_M1_vrr = QCZ*I_ERI_G2x2z_S_Pz_S_M1_vrr+WQZ*I_ERI_G2x2z_S_Pz_S_M2_vrr+oned2e*I_ERI_G2x2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_Pz_S_M2_vrr;
      Double I_ERI_Gx3y_S_D2z_S_M1_vrr = QCZ*I_ERI_Gx3y_S_Pz_S_M1_vrr+WQZ*I_ERI_Gx3y_S_Pz_S_M2_vrr+oned2e*I_ERI_Gx3y_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_D2z_S_M1_vrr = QCZ*I_ERI_Gx2yz_S_Pz_S_M1_vrr+WQZ*I_ERI_Gx2yz_S_Pz_S_M2_vrr+oned2e*I_ERI_Gx2yz_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx2yz_S_S_S_M2_vrr+oned2k*I_ERI_Fx2y_S_Pz_S_M2_vrr;
      Double I_ERI_Gxy2z_S_D2z_S_M1_vrr = QCZ*I_ERI_Gxy2z_S_Pz_S_M1_vrr+WQZ*I_ERI_Gxy2z_S_Pz_S_M2_vrr+oned2e*I_ERI_Gxy2z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gxy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_Pz_S_M2_vrr;
      Double I_ERI_Gx3z_S_D2z_S_M1_vrr = QCZ*I_ERI_Gx3z_S_Pz_S_M1_vrr+WQZ*I_ERI_Gx3z_S_Pz_S_M2_vrr+oned2e*I_ERI_Gx3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gx3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Fx2z_S_Pz_S_M2_vrr;
      Double I_ERI_G4y_S_D2z_S_M1_vrr = QCZ*I_ERI_G4y_S_Pz_S_M1_vrr+WQZ*I_ERI_G4y_S_Pz_S_M2_vrr+oned2e*I_ERI_G4y_S_S_S_M1_vrr-rhod2esq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_D2z_S_M1_vrr = QCZ*I_ERI_G3yz_S_Pz_S_M1_vrr+WQZ*I_ERI_G3yz_S_Pz_S_M2_vrr+oned2e*I_ERI_G3yz_S_S_S_M1_vrr-rhod2esq*I_ERI_G3yz_S_S_S_M2_vrr+oned2k*I_ERI_F3y_S_Pz_S_M2_vrr;
      Double I_ERI_G2y2z_S_D2z_S_M1_vrr = QCZ*I_ERI_G2y2z_S_Pz_S_M1_vrr+WQZ*I_ERI_G2y2z_S_Pz_S_M2_vrr+oned2e*I_ERI_G2y2z_S_S_S_M1_vrr-rhod2esq*I_ERI_G2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_Pz_S_M2_vrr;
      Double I_ERI_Gy3z_S_D2z_S_M1_vrr = QCZ*I_ERI_Gy3z_S_Pz_S_M1_vrr+WQZ*I_ERI_Gy3z_S_Pz_S_M2_vrr+oned2e*I_ERI_Gy3z_S_S_S_M1_vrr-rhod2esq*I_ERI_Gy3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Fy2z_S_Pz_S_M2_vrr;
      Double I_ERI_G4z_S_D2z_S_M1_vrr = QCZ*I_ERI_G4z_S_Pz_S_M1_vrr+WQZ*I_ERI_G4z_S_Pz_S_M2_vrr+oned2e*I_ERI_G4z_S_S_S_M1_vrr-rhod2esq*I_ERI_G4z_S_S_S_M2_vrr+4*oned2k*I_ERI_F3z_S_Pz_S_M2_vrr;

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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_S_S_vrr = PAX*I_ERI_Px_S_S_S_vrr+WPX*I_ERI_Px_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_Dxy_S_S_S_vrr = PAY*I_ERI_Px_S_S_S_vrr+WPY*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_S_S_vrr = PAY*I_ERI_Py_S_S_S_vrr+WPY*I_ERI_Py_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
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
       * shell quartet name: SQ_ERI_G_S_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_D_S
       * RHS shell quartet name: SQ_ERI_G_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_P_S
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_F3x_S_vrr = QCX*I_ERI_G4x_S_D2x_S_vrr+WQX*I_ERI_G4x_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G4x_S_Px_S_vrr-2*rhod2esq*I_ERI_G4x_S_Px_S_M1_vrr+4*oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xy_S_F3x_S_vrr = QCX*I_ERI_G3xy_S_D2x_S_vrr+WQX*I_ERI_G3xy_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G3xy_S_Px_S_vrr-2*rhod2esq*I_ERI_G3xy_S_Px_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G3xz_S_F3x_S_vrr = QCX*I_ERI_G3xz_S_D2x_S_vrr+WQX*I_ERI_G3xz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G3xz_S_Px_S_vrr-2*rhod2esq*I_ERI_G3xz_S_Px_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3x_S_vrr = QCX*I_ERI_G2x2y_S_D2x_S_vrr+WQX*I_ERI_G2x2y_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G2x2y_S_Px_S_vrr-2*rhod2esq*I_ERI_G2x2y_S_Px_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3x_S_vrr = QCX*I_ERI_G2xyz_S_D2x_S_vrr+WQX*I_ERI_G2xyz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G2xyz_S_Px_S_vrr-2*rhod2esq*I_ERI_G2xyz_S_Px_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3x_S_vrr = QCX*I_ERI_G2x2z_S_D2x_S_vrr+WQX*I_ERI_G2x2z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G2x2z_S_Px_S_vrr-2*rhod2esq*I_ERI_G2x2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3x_S_vrr = QCX*I_ERI_Gx3y_S_D2x_S_vrr+WQX*I_ERI_Gx3y_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gx3y_S_Px_S_vrr-2*rhod2esq*I_ERI_Gx3y_S_Px_S_M1_vrr+oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3x_S_vrr = QCX*I_ERI_Gx2yz_S_D2x_S_vrr+WQX*I_ERI_Gx2yz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gx2yz_S_Px_S_vrr-2*rhod2esq*I_ERI_Gx2yz_S_Px_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3x_S_vrr = QCX*I_ERI_Gxy2z_S_D2x_S_vrr+WQX*I_ERI_Gxy2z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gxy2z_S_Px_S_vrr-2*rhod2esq*I_ERI_Gxy2z_S_Px_S_M1_vrr+oned2k*I_ERI_Fy2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3x_S_vrr = QCX*I_ERI_Gx3z_S_D2x_S_vrr+WQX*I_ERI_Gx3z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gx3z_S_Px_S_vrr-2*rhod2esq*I_ERI_Gx3z_S_Px_S_M1_vrr+oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4y_S_F3x_S_vrr = QCX*I_ERI_G4y_S_D2x_S_vrr+WQX*I_ERI_G4y_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G4y_S_Px_S_vrr-2*rhod2esq*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_G3yz_S_F3x_S_vrr = QCX*I_ERI_G3yz_S_D2x_S_vrr+WQX*I_ERI_G3yz_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G3yz_S_Px_S_vrr-2*rhod2esq*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3x_S_vrr = QCX*I_ERI_G2y2z_S_D2x_S_vrr+WQX*I_ERI_G2y2z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G2y2z_S_Px_S_vrr-2*rhod2esq*I_ERI_G2y2z_S_Px_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3x_S_vrr = QCX*I_ERI_Gy3z_S_D2x_S_vrr+WQX*I_ERI_Gy3z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_Gy3z_S_Px_S_vrr-2*rhod2esq*I_ERI_Gy3z_S_Px_S_M1_vrr;
      Double I_ERI_G4z_S_F3x_S_vrr = QCX*I_ERI_G4z_S_D2x_S_vrr+WQX*I_ERI_G4z_S_D2x_S_M1_vrr+2*oned2e*I_ERI_G4z_S_Px_S_vrr-2*rhod2esq*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_G4x_S_F2xy_S_vrr = QCY*I_ERI_G4x_S_D2x_S_vrr+WQY*I_ERI_G4x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xy_S_F2xy_S_vrr = QCY*I_ERI_G3xy_S_D2x_S_vrr+WQY*I_ERI_G3xy_S_D2x_S_M1_vrr+oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xz_S_F2xy_S_vrr = QCY*I_ERI_G3xz_S_D2x_S_vrr+WQY*I_ERI_G3xz_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2xy_S_vrr = QCY*I_ERI_G2x2y_S_D2x_S_vrr+WQY*I_ERI_G2x2y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2xy_S_vrr = QCY*I_ERI_G2xyz_S_D2x_S_vrr+WQY*I_ERI_G2xyz_S_D2x_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2xy_S_vrr = QCY*I_ERI_G2x2z_S_D2x_S_vrr+WQY*I_ERI_G2x2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2xy_S_vrr = QCY*I_ERI_Gx3y_S_D2x_S_vrr+WQY*I_ERI_Gx3y_S_D2x_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2xy_S_vrr = QCY*I_ERI_Gx2yz_S_D2x_S_vrr+WQY*I_ERI_Gx2yz_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2xy_S_vrr = QCY*I_ERI_Gxy2z_S_D2x_S_vrr+WQY*I_ERI_Gxy2z_S_D2x_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2xy_S_vrr = QCY*I_ERI_Gx3z_S_D2x_S_vrr+WQY*I_ERI_Gx3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4y_S_F2xy_S_vrr = QCY*I_ERI_G4y_S_D2x_S_vrr+WQY*I_ERI_G4y_S_D2x_S_M1_vrr+4*oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_G3yz_S_F2xy_S_vrr = QCY*I_ERI_G3yz_S_D2x_S_vrr+WQY*I_ERI_G3yz_S_D2x_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_D2x_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2xy_S_vrr = QCY*I_ERI_G2y2z_S_D2x_S_vrr+WQY*I_ERI_G2y2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2xy_S_vrr = QCY*I_ERI_Gy3z_S_D2x_S_vrr+WQY*I_ERI_Gy3z_S_D2x_S_M1_vrr+oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4z_S_F2xy_S_vrr = QCY*I_ERI_G4z_S_D2x_S_vrr+WQY*I_ERI_G4z_S_D2x_S_M1_vrr;
      Double I_ERI_G4x_S_F2xz_S_vrr = QCZ*I_ERI_G4x_S_D2x_S_vrr+WQZ*I_ERI_G4x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xy_S_F2xz_S_vrr = QCZ*I_ERI_G3xy_S_D2x_S_vrr+WQZ*I_ERI_G3xy_S_D2x_S_M1_vrr;
      Double I_ERI_G3xz_S_F2xz_S_vrr = QCZ*I_ERI_G3xz_S_D2x_S_vrr+WQZ*I_ERI_G3xz_S_D2x_S_M1_vrr+oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2xz_S_vrr = QCZ*I_ERI_G2x2y_S_D2x_S_vrr+WQZ*I_ERI_G2x2y_S_D2x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2xz_S_vrr = QCZ*I_ERI_G2xyz_S_D2x_S_vrr+WQZ*I_ERI_G2xyz_S_D2x_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2xz_S_vrr = QCZ*I_ERI_G2x2z_S_D2x_S_vrr+WQZ*I_ERI_G2x2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2xz_S_vrr = QCZ*I_ERI_Gx3y_S_D2x_S_vrr+WQZ*I_ERI_Gx3y_S_D2x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2xz_S_vrr = QCZ*I_ERI_Gx2yz_S_D2x_S_vrr+WQZ*I_ERI_Gx2yz_S_D2x_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2xz_S_vrr = QCZ*I_ERI_Gxy2z_S_D2x_S_vrr+WQZ*I_ERI_Gxy2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2xz_S_vrr = QCZ*I_ERI_Gx3z_S_D2x_S_vrr+WQZ*I_ERI_Gx3z_S_D2x_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_G4y_S_F2xz_S_vrr = QCZ*I_ERI_G4y_S_D2x_S_vrr+WQZ*I_ERI_G4y_S_D2x_S_M1_vrr;
      Double I_ERI_G3yz_S_F2xz_S_vrr = QCZ*I_ERI_G3yz_S_D2x_S_vrr+WQZ*I_ERI_G3yz_S_D2x_S_M1_vrr+oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2xz_S_vrr = QCZ*I_ERI_G2y2z_S_D2x_S_vrr+WQZ*I_ERI_G2y2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_D2x_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2xz_S_vrr = QCZ*I_ERI_Gy3z_S_D2x_S_vrr+WQZ*I_ERI_Gy3z_S_D2x_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_D2x_S_M1_vrr;
      Double I_ERI_G4z_S_F2xz_S_vrr = QCZ*I_ERI_G4z_S_D2x_S_vrr+WQZ*I_ERI_G4z_S_D2x_S_M1_vrr+4*oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4x_S_Fx2y_S_vrr = QCX*I_ERI_G4x_S_D2y_S_vrr+WQX*I_ERI_G4x_S_D2y_S_M1_vrr+4*oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xy_S_Fx2y_S_vrr = QCX*I_ERI_G3xy_S_D2y_S_vrr+WQX*I_ERI_G3xy_S_D2y_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G3xz_S_Fx2y_S_vrr = QCX*I_ERI_G3xz_S_D2y_S_vrr+WQX*I_ERI_G3xz_S_D2y_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fx2y_S_vrr = QCX*I_ERI_G2x2y_S_D2y_S_vrr+WQX*I_ERI_G2x2y_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fx2y_S_vrr = QCX*I_ERI_G2xyz_S_D2y_S_vrr+WQX*I_ERI_G2xyz_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fx2y_S_vrr = QCX*I_ERI_G2x2z_S_D2y_S_vrr+WQX*I_ERI_G2x2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fx2y_S_vrr = QCX*I_ERI_Gx3y_S_D2y_S_vrr+WQX*I_ERI_Gx3y_S_D2y_S_M1_vrr+oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fx2y_S_vrr = QCX*I_ERI_Gx2yz_S_D2y_S_vrr+WQX*I_ERI_Gx2yz_S_D2y_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fx2y_S_vrr = QCX*I_ERI_Gxy2z_S_D2y_S_vrr+WQX*I_ERI_Gxy2z_S_D2y_S_M1_vrr+oned2k*I_ERI_Fy2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fx2y_S_vrr = QCX*I_ERI_Gx3z_S_D2y_S_vrr+WQX*I_ERI_Gx3z_S_D2y_S_M1_vrr+oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4y_S_Fx2y_S_vrr = QCX*I_ERI_G4y_S_D2y_S_vrr+WQX*I_ERI_G4y_S_D2y_S_M1_vrr;
      Double I_ERI_G3yz_S_Fx2y_S_vrr = QCX*I_ERI_G3yz_S_D2y_S_vrr+WQX*I_ERI_G3yz_S_D2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fx2y_S_vrr = QCX*I_ERI_G2y2z_S_D2y_S_vrr+WQX*I_ERI_G2y2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fx2y_S_vrr = QCX*I_ERI_Gy3z_S_D2y_S_vrr+WQX*I_ERI_Gy3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4z_S_Fx2y_S_vrr = QCX*I_ERI_G4z_S_D2y_S_vrr+WQX*I_ERI_G4z_S_D2y_S_M1_vrr;
      Double I_ERI_G4x_S_Fxyz_S_vrr = QCZ*I_ERI_G4x_S_Dxy_S_vrr+WQZ*I_ERI_G4x_S_Dxy_S_M1_vrr;
      Double I_ERI_G3xy_S_Fxyz_S_vrr = QCZ*I_ERI_G3xy_S_Dxy_S_vrr+WQZ*I_ERI_G3xy_S_Dxy_S_M1_vrr;
      Double I_ERI_G3xz_S_Fxyz_S_vrr = QCZ*I_ERI_G3xz_S_Dxy_S_vrr+WQZ*I_ERI_G3xz_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fxyz_S_vrr = QCZ*I_ERI_G2x2y_S_Dxy_S_vrr+WQZ*I_ERI_G2x2y_S_Dxy_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fxyz_S_vrr = QCZ*I_ERI_G2xyz_S_Dxy_S_vrr+WQZ*I_ERI_G2xyz_S_Dxy_S_M1_vrr+oned2k*I_ERI_F2xy_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fxyz_S_vrr = QCZ*I_ERI_G2x2z_S_Dxy_S_vrr+WQZ*I_ERI_G2x2z_S_Dxy_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fxyz_S_vrr = QCZ*I_ERI_Gx3y_S_Dxy_S_vrr+WQZ*I_ERI_Gx3y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fxyz_S_vrr = QCZ*I_ERI_Gx2yz_S_Dxy_S_vrr+WQZ*I_ERI_Gx2yz_S_Dxy_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fxyz_S_vrr = QCZ*I_ERI_Gxy2z_S_Dxy_S_vrr+WQZ*I_ERI_Gxy2z_S_Dxy_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fxyz_S_vrr = QCZ*I_ERI_Gx3z_S_Dxy_S_vrr+WQZ*I_ERI_Gx3z_S_Dxy_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4y_S_Fxyz_S_vrr = QCZ*I_ERI_G4y_S_Dxy_S_vrr+WQZ*I_ERI_G4y_S_Dxy_S_M1_vrr;
      Double I_ERI_G3yz_S_Fxyz_S_vrr = QCZ*I_ERI_G3yz_S_Dxy_S_vrr+WQZ*I_ERI_G3yz_S_Dxy_S_M1_vrr+oned2k*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fxyz_S_vrr = QCZ*I_ERI_G2y2z_S_Dxy_S_vrr+WQZ*I_ERI_G2y2z_S_Dxy_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fxyz_S_vrr = QCZ*I_ERI_Gy3z_S_Dxy_S_vrr+WQZ*I_ERI_Gy3z_S_Dxy_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4z_S_Fxyz_S_vrr = QCZ*I_ERI_G4z_S_Dxy_S_vrr+WQZ*I_ERI_G4z_S_Dxy_S_M1_vrr+4*oned2k*I_ERI_F3z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4x_S_Fx2z_S_vrr = QCX*I_ERI_G4x_S_D2z_S_vrr+WQX*I_ERI_G4x_S_D2z_S_M1_vrr+4*oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Fx2z_S_vrr = QCX*I_ERI_G3xy_S_D2z_S_vrr+WQX*I_ERI_G3xy_S_D2z_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Fx2z_S_vrr = QCX*I_ERI_G3xz_S_D2z_S_vrr+WQX*I_ERI_G3xz_S_D2z_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fx2z_S_vrr = QCX*I_ERI_G2x2y_S_D2z_S_vrr+WQX*I_ERI_G2x2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_D2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fx2z_S_vrr = QCX*I_ERI_G2xyz_S_D2z_S_vrr+WQX*I_ERI_G2xyz_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fx2z_S_vrr = QCX*I_ERI_G2x2z_S_D2z_S_vrr+WQX*I_ERI_G2x2z_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fx2z_S_vrr = QCX*I_ERI_Gx3y_S_D2z_S_vrr+WQX*I_ERI_Gx3y_S_D2z_S_M1_vrr+oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fx2z_S_vrr = QCX*I_ERI_Gx2yz_S_D2z_S_vrr+WQX*I_ERI_Gx2yz_S_D2z_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fx2z_S_vrr = QCX*I_ERI_Gxy2z_S_D2z_S_vrr+WQX*I_ERI_Gxy2z_S_D2z_S_M1_vrr+oned2k*I_ERI_Fy2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fx2z_S_vrr = QCX*I_ERI_Gx3z_S_D2z_S_vrr+WQX*I_ERI_Gx3z_S_D2z_S_M1_vrr+oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4y_S_Fx2z_S_vrr = QCX*I_ERI_G4y_S_D2z_S_vrr+WQX*I_ERI_G4y_S_D2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Fx2z_S_vrr = QCX*I_ERI_G3yz_S_D2z_S_vrr+WQX*I_ERI_G3yz_S_D2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fx2z_S_vrr = QCX*I_ERI_G2y2z_S_D2z_S_vrr+WQX*I_ERI_G2y2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fx2z_S_vrr = QCX*I_ERI_Gy3z_S_D2z_S_vrr+WQX*I_ERI_Gy3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4z_S_Fx2z_S_vrr = QCX*I_ERI_G4z_S_D2z_S_vrr+WQX*I_ERI_G4z_S_D2z_S_M1_vrr;
      Double I_ERI_G4x_S_F3y_S_vrr = QCY*I_ERI_G4x_S_D2y_S_vrr+WQY*I_ERI_G4x_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G4x_S_Py_S_vrr-2*rhod2esq*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_G3xy_S_F3y_S_vrr = QCY*I_ERI_G3xy_S_D2y_S_vrr+WQY*I_ERI_G3xy_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G3xy_S_Py_S_vrr-2*rhod2esq*I_ERI_G3xy_S_Py_S_M1_vrr+oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xz_S_F3y_S_vrr = QCY*I_ERI_G3xz_S_D2y_S_vrr+WQY*I_ERI_G3xz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G3xz_S_Py_S_vrr-2*rhod2esq*I_ERI_G3xz_S_Py_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3y_S_vrr = QCY*I_ERI_G2x2y_S_D2y_S_vrr+WQY*I_ERI_G2x2y_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G2x2y_S_Py_S_vrr-2*rhod2esq*I_ERI_G2x2y_S_Py_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3y_S_vrr = QCY*I_ERI_G2xyz_S_D2y_S_vrr+WQY*I_ERI_G2xyz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G2xyz_S_Py_S_vrr-2*rhod2esq*I_ERI_G2xyz_S_Py_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3y_S_vrr = QCY*I_ERI_G2x2z_S_D2y_S_vrr+WQY*I_ERI_G2x2z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G2x2z_S_Py_S_vrr-2*rhod2esq*I_ERI_G2x2z_S_Py_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3y_S_vrr = QCY*I_ERI_Gx3y_S_D2y_S_vrr+WQY*I_ERI_Gx3y_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gx3y_S_Py_S_vrr-2*rhod2esq*I_ERI_Gx3y_S_Py_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3y_S_vrr = QCY*I_ERI_Gx2yz_S_D2y_S_vrr+WQY*I_ERI_Gx2yz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gx2yz_S_Py_S_vrr-2*rhod2esq*I_ERI_Gx2yz_S_Py_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3y_S_vrr = QCY*I_ERI_Gxy2z_S_D2y_S_vrr+WQY*I_ERI_Gxy2z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gxy2z_S_Py_S_vrr-2*rhod2esq*I_ERI_Gxy2z_S_Py_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3y_S_vrr = QCY*I_ERI_Gx3z_S_D2y_S_vrr+WQY*I_ERI_Gx3z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gx3z_S_Py_S_vrr-2*rhod2esq*I_ERI_Gx3z_S_Py_S_M1_vrr;
      Double I_ERI_G4y_S_F3y_S_vrr = QCY*I_ERI_G4y_S_D2y_S_vrr+WQY*I_ERI_G4y_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G4y_S_Py_S_vrr-2*rhod2esq*I_ERI_G4y_S_Py_S_M1_vrr+4*oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_G3yz_S_F3y_S_vrr = QCY*I_ERI_G3yz_S_D2y_S_vrr+WQY*I_ERI_G3yz_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G3yz_S_Py_S_vrr-2*rhod2esq*I_ERI_G3yz_S_Py_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_D2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3y_S_vrr = QCY*I_ERI_G2y2z_S_D2y_S_vrr+WQY*I_ERI_G2y2z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G2y2z_S_Py_S_vrr-2*rhod2esq*I_ERI_G2y2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3y_S_vrr = QCY*I_ERI_Gy3z_S_D2y_S_vrr+WQY*I_ERI_Gy3z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_Gy3z_S_Py_S_vrr-2*rhod2esq*I_ERI_Gy3z_S_Py_S_M1_vrr+oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4z_S_F3y_S_vrr = QCY*I_ERI_G4z_S_D2y_S_vrr+WQY*I_ERI_G4z_S_D2y_S_M1_vrr+2*oned2e*I_ERI_G4z_S_Py_S_vrr-2*rhod2esq*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_G4x_S_F2yz_S_vrr = QCZ*I_ERI_G4x_S_D2y_S_vrr+WQZ*I_ERI_G4x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xy_S_F2yz_S_vrr = QCZ*I_ERI_G3xy_S_D2y_S_vrr+WQZ*I_ERI_G3xy_S_D2y_S_M1_vrr;
      Double I_ERI_G3xz_S_F2yz_S_vrr = QCZ*I_ERI_G3xz_S_D2y_S_vrr+WQZ*I_ERI_G3xz_S_D2y_S_M1_vrr+oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2yz_S_vrr = QCZ*I_ERI_G2x2y_S_D2y_S_vrr+WQZ*I_ERI_G2x2y_S_D2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2yz_S_vrr = QCZ*I_ERI_G2xyz_S_D2y_S_vrr+WQZ*I_ERI_G2xyz_S_D2y_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2yz_S_vrr = QCZ*I_ERI_G2x2z_S_D2y_S_vrr+WQZ*I_ERI_G2x2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2yz_S_vrr = QCZ*I_ERI_Gx3y_S_D2y_S_vrr+WQZ*I_ERI_Gx3y_S_D2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2yz_S_vrr = QCZ*I_ERI_Gx2yz_S_D2y_S_vrr+WQZ*I_ERI_Gx2yz_S_D2y_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2yz_S_vrr = QCZ*I_ERI_Gxy2z_S_D2y_S_vrr+WQZ*I_ERI_Gxy2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2yz_S_vrr = QCZ*I_ERI_Gx3z_S_D2y_S_vrr+WQZ*I_ERI_Gx3z_S_D2y_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_D2y_S_M1_vrr;
      Double I_ERI_G4y_S_F2yz_S_vrr = QCZ*I_ERI_G4y_S_D2y_S_vrr+WQZ*I_ERI_G4y_S_D2y_S_M1_vrr;
      Double I_ERI_G3yz_S_F2yz_S_vrr = QCZ*I_ERI_G3yz_S_D2y_S_vrr+WQZ*I_ERI_G3yz_S_D2y_S_M1_vrr+oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2yz_S_vrr = QCZ*I_ERI_G2y2z_S_D2y_S_vrr+WQZ*I_ERI_G2y2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_D2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2yz_S_vrr = QCZ*I_ERI_Gy3z_S_D2y_S_vrr+WQZ*I_ERI_Gy3z_S_D2y_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_D2y_S_M1_vrr;
      Double I_ERI_G4z_S_F2yz_S_vrr = QCZ*I_ERI_G4z_S_D2y_S_vrr+WQZ*I_ERI_G4z_S_D2y_S_M1_vrr+4*oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4x_S_Fy2z_S_vrr = QCY*I_ERI_G4x_S_D2z_S_vrr+WQY*I_ERI_G4x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Fy2z_S_vrr = QCY*I_ERI_G3xy_S_D2z_S_vrr+WQY*I_ERI_G3xy_S_D2z_S_M1_vrr+oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Fy2z_S_vrr = QCY*I_ERI_G3xz_S_D2z_S_vrr+WQY*I_ERI_G3xz_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fy2z_S_vrr = QCY*I_ERI_G2x2y_S_D2z_S_vrr+WQY*I_ERI_G2x2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fy2z_S_vrr = QCY*I_ERI_G2xyz_S_D2z_S_vrr+WQY*I_ERI_G2xyz_S_D2z_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fy2z_S_vrr = QCY*I_ERI_G2x2z_S_D2z_S_vrr+WQY*I_ERI_G2x2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fy2z_S_vrr = QCY*I_ERI_Gx3y_S_D2z_S_vrr+WQY*I_ERI_Gx3y_S_D2z_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_D2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fy2z_S_vrr = QCY*I_ERI_Gx2yz_S_D2z_S_vrr+WQY*I_ERI_Gx2yz_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fy2z_S_vrr = QCY*I_ERI_Gxy2z_S_D2z_S_vrr+WQY*I_ERI_Gxy2z_S_D2z_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fy2z_S_vrr = QCY*I_ERI_Gx3z_S_D2z_S_vrr+WQY*I_ERI_Gx3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4y_S_Fy2z_S_vrr = QCY*I_ERI_G4y_S_D2z_S_vrr+WQY*I_ERI_G4y_S_D2z_S_M1_vrr+4*oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Fy2z_S_vrr = QCY*I_ERI_G3yz_S_D2z_S_vrr+WQY*I_ERI_G3yz_S_D2z_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_D2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fy2z_S_vrr = QCY*I_ERI_G2y2z_S_D2z_S_vrr+WQY*I_ERI_G2y2z_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Fy2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fy2z_S_vrr = QCY*I_ERI_Gy3z_S_D2z_S_vrr+WQY*I_ERI_Gy3z_S_D2z_S_M1_vrr+oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4z_S_Fy2z_S_vrr = QCY*I_ERI_G4z_S_D2z_S_vrr+WQY*I_ERI_G4z_S_D2z_S_M1_vrr;
      Double I_ERI_G4x_S_F3z_S_vrr = QCZ*I_ERI_G4x_S_D2z_S_vrr+WQZ*I_ERI_G4x_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G4x_S_Pz_S_vrr-2*rhod2esq*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_G3xy_S_F3z_S_vrr = QCZ*I_ERI_G3xy_S_D2z_S_vrr+WQZ*I_ERI_G3xy_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G3xy_S_Pz_S_vrr-2*rhod2esq*I_ERI_G3xy_S_Pz_S_M1_vrr;
      Double I_ERI_G3xz_S_F3z_S_vrr = QCZ*I_ERI_G3xz_S_D2z_S_vrr+WQZ*I_ERI_G3xz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G3xz_S_Pz_S_vrr-2*rhod2esq*I_ERI_G3xz_S_Pz_S_M1_vrr+oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3z_S_vrr = QCZ*I_ERI_G2x2y_S_D2z_S_vrr+WQZ*I_ERI_G2x2y_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G2x2y_S_Pz_S_vrr-2*rhod2esq*I_ERI_G2x2y_S_Pz_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3z_S_vrr = QCZ*I_ERI_G2xyz_S_D2z_S_vrr+WQZ*I_ERI_G2xyz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G2xyz_S_Pz_S_vrr-2*rhod2esq*I_ERI_G2xyz_S_Pz_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3z_S_vrr = QCZ*I_ERI_G2x2z_S_D2z_S_vrr+WQZ*I_ERI_G2x2z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G2x2z_S_Pz_S_vrr-2*rhod2esq*I_ERI_G2x2z_S_Pz_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3z_S_vrr = QCZ*I_ERI_Gx3y_S_D2z_S_vrr+WQZ*I_ERI_Gx3y_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gx3y_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gx3y_S_Pz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3z_S_vrr = QCZ*I_ERI_Gx2yz_S_D2z_S_vrr+WQZ*I_ERI_Gx2yz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gx2yz_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gx2yz_S_Pz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3z_S_vrr = QCZ*I_ERI_Gxy2z_S_D2z_S_vrr+WQZ*I_ERI_Gxy2z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gxy2z_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gxy2z_S_Pz_S_M1_vrr+2*oned2k*I_ERI_Fxyz_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3z_S_vrr = QCZ*I_ERI_Gx3z_S_D2z_S_vrr+WQZ*I_ERI_Gx3z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gx3z_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gx3z_S_Pz_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_G4y_S_F3z_S_vrr = QCZ*I_ERI_G4y_S_D2z_S_vrr+WQZ*I_ERI_G4y_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G4y_S_Pz_S_vrr-2*rhod2esq*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_G3yz_S_F3z_S_vrr = QCZ*I_ERI_G3yz_S_D2z_S_vrr+WQZ*I_ERI_G3yz_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G3yz_S_Pz_S_vrr-2*rhod2esq*I_ERI_G3yz_S_Pz_S_M1_vrr+oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3z_S_vrr = QCZ*I_ERI_G2y2z_S_D2z_S_vrr+WQZ*I_ERI_G2y2z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G2y2z_S_Pz_S_vrr-2*rhod2esq*I_ERI_G2y2z_S_Pz_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_D2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3z_S_vrr = QCZ*I_ERI_Gy3z_S_D2z_S_vrr+WQZ*I_ERI_Gy3z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_Gy3z_S_Pz_S_vrr-2*rhod2esq*I_ERI_Gy3z_S_Pz_S_M1_vrr+3*oned2k*I_ERI_Fy2z_S_D2z_S_M1_vrr;
      Double I_ERI_G4z_S_F3z_S_vrr = QCZ*I_ERI_G4z_S_D2z_S_vrr+WQZ*I_ERI_G4z_S_D2z_S_M1_vrr+2*oned2e*I_ERI_G4z_S_Pz_S_vrr-2*rhod2esq*I_ERI_G4z_S_Pz_S_M1_vrr+4*oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C2000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C2000003_coefs = ic2*jc2;
      abcd[0] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3x_S_D2x_S_vrr;
      abcd[1] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      abcd[2] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      abcd[3] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      abcd[4] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      abcd[5] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      abcd[6] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3y_S_D2x_S_vrr;
      abcd[7] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      abcd[8] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      abcd[9] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3z_S_D2x_S_vrr;
      abcd[40] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      abcd[41] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      abcd[42] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      abcd[43] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      abcd[44] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      abcd[45] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      abcd[46] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      abcd[47] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      abcd[48] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      abcd[49] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      abcd[80] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      abcd[81] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      abcd[82] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      abcd[83] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      abcd[84] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      abcd[85] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      abcd[86] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      abcd[87] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      abcd[88] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      abcd[89] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      abcd[120] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3x_S_D2y_S_vrr;
      abcd[121] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      abcd[122] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      abcd[123] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      abcd[124] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      abcd[125] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      abcd[126] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3y_S_D2y_S_vrr;
      abcd[127] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      abcd[128] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      abcd[129] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3z_S_D2y_S_vrr;
      abcd[160] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      abcd[161] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      abcd[162] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      abcd[163] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      abcd[164] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      abcd[165] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      abcd[166] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      abcd[167] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      abcd[168] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      abcd[169] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      abcd[200] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3x_S_D2z_S_vrr;
      abcd[201] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      abcd[202] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      abcd[203] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      abcd[204] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      abcd[205] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      abcd[206] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3y_S_D2z_S_vrr;
      abcd[207] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      abcd[208] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      abcd[209] += SQ_ERI_F_S_D_S_C2000003_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_C2001003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_C2001003_coefs = ic2_1*jc2;
      I_ERI_G4x_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_C2001003 += SQ_ERI_G_S_D_S_C2001003_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C2001003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C2001003_coefs = ic2_1*jc2;
      I_ERI_F3x_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C2001003 += SQ_ERI_F_S_D_S_C2001003_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_C1002000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_C1002000003_coefs = ic2*jc2_1;
      I_ERI_F3x_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_C1002000003 += SQ_ERI_F_S_F_S_C1002000003_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1002000003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1002000003_coefs = ic2*jc2_1;
      I_ERI_F3x_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1002000003 += SQ_ERI_F_S_D_S_C1002000003_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_F_S_C1002001003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_F_S_C1002001003_coefs = ic2_1*jc2_1;
      I_ERI_G4x_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_F3x_S_vrr;
      I_ERI_G3xy_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_F3x_S_vrr;
      I_ERI_G3xz_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_F3x_S_vrr;
      I_ERI_G2x2y_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_F3x_S_vrr;
      I_ERI_G2xyz_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_F3x_S_vrr;
      I_ERI_G2x2z_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_F3x_S_vrr;
      I_ERI_Gx3y_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_F3x_S_vrr;
      I_ERI_Gx2yz_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_F3x_S_vrr;
      I_ERI_Gxy2z_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_F3x_S_vrr;
      I_ERI_Gx3z_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_F3x_S_vrr;
      I_ERI_G4y_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_F3x_S_vrr;
      I_ERI_G3yz_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_F3x_S_vrr;
      I_ERI_G2y2z_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_F3x_S_vrr;
      I_ERI_Gy3z_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_F3x_S_vrr;
      I_ERI_G4z_S_F3x_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_F3x_S_vrr;
      I_ERI_G4x_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_F2xy_S_vrr;
      I_ERI_G3xy_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_F2xy_S_vrr;
      I_ERI_G3xz_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_F2xy_S_vrr;
      I_ERI_G2x2y_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_F2xy_S_vrr;
      I_ERI_G2xyz_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_F2xy_S_vrr;
      I_ERI_G2x2z_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_F2xy_S_vrr;
      I_ERI_Gx3y_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_F2xy_S_vrr;
      I_ERI_Gx2yz_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_F2xy_S_vrr;
      I_ERI_Gxy2z_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_F2xy_S_vrr;
      I_ERI_Gx3z_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_F2xy_S_vrr;
      I_ERI_G4y_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_F2xy_S_vrr;
      I_ERI_G3yz_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_F2xy_S_vrr;
      I_ERI_G2y2z_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_F2xy_S_vrr;
      I_ERI_Gy3z_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_F2xy_S_vrr;
      I_ERI_G4z_S_F2xy_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_F2xy_S_vrr;
      I_ERI_G4x_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_F2xz_S_vrr;
      I_ERI_G3xy_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_F2xz_S_vrr;
      I_ERI_G3xz_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_F2xz_S_vrr;
      I_ERI_G2x2y_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_F2xz_S_vrr;
      I_ERI_G2xyz_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_F2xz_S_vrr;
      I_ERI_G2x2z_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_F2xz_S_vrr;
      I_ERI_Gx3y_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_F2xz_S_vrr;
      I_ERI_Gx2yz_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_F2xz_S_vrr;
      I_ERI_Gxy2z_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_F2xz_S_vrr;
      I_ERI_Gx3z_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_F2xz_S_vrr;
      I_ERI_G4y_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_F2xz_S_vrr;
      I_ERI_G3yz_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_F2xz_S_vrr;
      I_ERI_G2y2z_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_F2xz_S_vrr;
      I_ERI_Gy3z_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_F2xz_S_vrr;
      I_ERI_G4z_S_F2xz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_F2xz_S_vrr;
      I_ERI_G4x_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_Fx2y_S_vrr;
      I_ERI_G3xy_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_Fx2y_S_vrr;
      I_ERI_G3xz_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_Fx2y_S_vrr;
      I_ERI_G2x2y_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_Fx2y_S_vrr;
      I_ERI_G2xyz_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_Fx2y_S_vrr;
      I_ERI_G2x2z_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_Fx2y_S_vrr;
      I_ERI_Gx3y_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_Fx2y_S_vrr;
      I_ERI_Gx2yz_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_Fx2y_S_vrr;
      I_ERI_Gxy2z_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_Fx2y_S_vrr;
      I_ERI_Gx3z_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_Fx2y_S_vrr;
      I_ERI_G4y_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_Fx2y_S_vrr;
      I_ERI_G3yz_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_Fx2y_S_vrr;
      I_ERI_G2y2z_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_Fx2y_S_vrr;
      I_ERI_Gy3z_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_Fx2y_S_vrr;
      I_ERI_G4z_S_Fx2y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_Fx2y_S_vrr;
      I_ERI_G4x_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_Fxyz_S_vrr;
      I_ERI_G3xy_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_Fxyz_S_vrr;
      I_ERI_G3xz_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_Fxyz_S_vrr;
      I_ERI_G2x2y_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_Fxyz_S_vrr;
      I_ERI_G2xyz_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_Fxyz_S_vrr;
      I_ERI_G2x2z_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_Fxyz_S_vrr;
      I_ERI_Gx3y_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_Fxyz_S_vrr;
      I_ERI_Gx2yz_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_Fxyz_S_vrr;
      I_ERI_Gxy2z_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_Fxyz_S_vrr;
      I_ERI_Gx3z_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_Fxyz_S_vrr;
      I_ERI_G4y_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_Fxyz_S_vrr;
      I_ERI_G3yz_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_Fxyz_S_vrr;
      I_ERI_G2y2z_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_Fxyz_S_vrr;
      I_ERI_Gy3z_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_Fxyz_S_vrr;
      I_ERI_G4z_S_Fxyz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_Fxyz_S_vrr;
      I_ERI_G4x_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_Fx2z_S_vrr;
      I_ERI_G3xy_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_Fx2z_S_vrr;
      I_ERI_G3xz_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_Fx2z_S_vrr;
      I_ERI_G2x2y_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_Fx2z_S_vrr;
      I_ERI_G2xyz_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_Fx2z_S_vrr;
      I_ERI_G2x2z_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_Fx2z_S_vrr;
      I_ERI_Gx3y_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_Fx2z_S_vrr;
      I_ERI_Gx2yz_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_Fx2z_S_vrr;
      I_ERI_Gxy2z_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_Fx2z_S_vrr;
      I_ERI_Gx3z_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_Fx2z_S_vrr;
      I_ERI_G4y_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_Fx2z_S_vrr;
      I_ERI_G3yz_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_Fx2z_S_vrr;
      I_ERI_G2y2z_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_Fx2z_S_vrr;
      I_ERI_Gy3z_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_Fx2z_S_vrr;
      I_ERI_G4z_S_Fx2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_Fx2z_S_vrr;
      I_ERI_G4x_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_F3y_S_vrr;
      I_ERI_G3xy_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_F3y_S_vrr;
      I_ERI_G3xz_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_F3y_S_vrr;
      I_ERI_G2x2y_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_F3y_S_vrr;
      I_ERI_G2xyz_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_F3y_S_vrr;
      I_ERI_G2x2z_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_F3y_S_vrr;
      I_ERI_Gx3y_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_F3y_S_vrr;
      I_ERI_Gx2yz_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_F3y_S_vrr;
      I_ERI_Gxy2z_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_F3y_S_vrr;
      I_ERI_Gx3z_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_F3y_S_vrr;
      I_ERI_G4y_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_F3y_S_vrr;
      I_ERI_G3yz_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_F3y_S_vrr;
      I_ERI_G2y2z_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_F3y_S_vrr;
      I_ERI_Gy3z_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_F3y_S_vrr;
      I_ERI_G4z_S_F3y_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_F3y_S_vrr;
      I_ERI_G4x_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_F2yz_S_vrr;
      I_ERI_G3xy_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_F2yz_S_vrr;
      I_ERI_G3xz_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_F2yz_S_vrr;
      I_ERI_G2x2y_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_F2yz_S_vrr;
      I_ERI_G2xyz_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_F2yz_S_vrr;
      I_ERI_G2x2z_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_F2yz_S_vrr;
      I_ERI_Gx3y_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_F2yz_S_vrr;
      I_ERI_Gx2yz_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_F2yz_S_vrr;
      I_ERI_Gxy2z_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_F2yz_S_vrr;
      I_ERI_Gx3z_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_F2yz_S_vrr;
      I_ERI_G4y_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_F2yz_S_vrr;
      I_ERI_G3yz_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_F2yz_S_vrr;
      I_ERI_G2y2z_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_F2yz_S_vrr;
      I_ERI_Gy3z_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_F2yz_S_vrr;
      I_ERI_G4z_S_F2yz_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_F2yz_S_vrr;
      I_ERI_G4x_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_Fy2z_S_vrr;
      I_ERI_G3xy_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_Fy2z_S_vrr;
      I_ERI_G3xz_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_Fy2z_S_vrr;
      I_ERI_G2x2y_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_Fy2z_S_vrr;
      I_ERI_G2xyz_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_Fy2z_S_vrr;
      I_ERI_G2x2z_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_Fy2z_S_vrr;
      I_ERI_Gx3y_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_Fy2z_S_vrr;
      I_ERI_Gx2yz_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_Fy2z_S_vrr;
      I_ERI_Gxy2z_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_Fy2z_S_vrr;
      I_ERI_Gx3z_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_Fy2z_S_vrr;
      I_ERI_G4y_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_Fy2z_S_vrr;
      I_ERI_G3yz_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_Fy2z_S_vrr;
      I_ERI_G2y2z_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_Fy2z_S_vrr;
      I_ERI_Gy3z_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_Fy2z_S_vrr;
      I_ERI_G4z_S_Fy2z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_Fy2z_S_vrr;
      I_ERI_G4x_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4x_S_F3z_S_vrr;
      I_ERI_G3xy_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xy_S_F3z_S_vrr;
      I_ERI_G3xz_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3xz_S_F3z_S_vrr;
      I_ERI_G2x2y_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2y_S_F3z_S_vrr;
      I_ERI_G2xyz_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2xyz_S_F3z_S_vrr;
      I_ERI_G2x2z_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2x2z_S_F3z_S_vrr;
      I_ERI_Gx3y_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3y_S_F3z_S_vrr;
      I_ERI_Gx2yz_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx2yz_S_F3z_S_vrr;
      I_ERI_Gxy2z_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gxy2z_S_F3z_S_vrr;
      I_ERI_Gx3z_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gx3z_S_F3z_S_vrr;
      I_ERI_G4y_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4y_S_F3z_S_vrr;
      I_ERI_G3yz_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G3yz_S_F3z_S_vrr;
      I_ERI_G2y2z_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G2y2z_S_F3z_S_vrr;
      I_ERI_Gy3z_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_Gy3z_S_F3z_S_vrr;
      I_ERI_G4z_S_F3z_S_C1002001003 += SQ_ERI_G_S_F_S_C1002001003_coefs*I_ERI_G4z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_D_S_C1002001003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_G_S_D_S_C1002001003_coefs = ic2_1*jc2_1;
      I_ERI_G4x_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S_C1002001003 += SQ_ERI_G_S_D_S_C1002001003_coefs*I_ERI_G4z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_C1002001003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_C1002001003_coefs = ic2_1*jc2_1;
      I_ERI_F3x_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_C1002001003 += SQ_ERI_F_S_F_S_C1002001003_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1002001003
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1002001003_coefs = ic2_1*jc2_1;
      I_ERI_F3x_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1002001003 += SQ_ERI_F_S_D_S_C1002001003_coefs*I_ERI_F3z_S_D2z_S_vrr;
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
  Double CDX = C[0] - D[0];
  Double CDY = C[1] - D[1];
  Double CDZ = C[2] - D[2];

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_P_C1002000003
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_C1002000003
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1002000003
   ************************************************************/
  abcd[240] = I_ERI_F3x_S_F3x_S_C1002000003+CDX*I_ERI_F3x_S_D2x_S_C1002000003;
  abcd[241] = I_ERI_F2xy_S_F3x_S_C1002000003+CDX*I_ERI_F2xy_S_D2x_S_C1002000003;
  abcd[242] = I_ERI_F2xz_S_F3x_S_C1002000003+CDX*I_ERI_F2xz_S_D2x_S_C1002000003;
  abcd[243] = I_ERI_Fx2y_S_F3x_S_C1002000003+CDX*I_ERI_Fx2y_S_D2x_S_C1002000003;
  abcd[244] = I_ERI_Fxyz_S_F3x_S_C1002000003+CDX*I_ERI_Fxyz_S_D2x_S_C1002000003;
  abcd[245] = I_ERI_Fx2z_S_F3x_S_C1002000003+CDX*I_ERI_Fx2z_S_D2x_S_C1002000003;
  abcd[246] = I_ERI_F3y_S_F3x_S_C1002000003+CDX*I_ERI_F3y_S_D2x_S_C1002000003;
  abcd[247] = I_ERI_F2yz_S_F3x_S_C1002000003+CDX*I_ERI_F2yz_S_D2x_S_C1002000003;
  abcd[248] = I_ERI_Fy2z_S_F3x_S_C1002000003+CDX*I_ERI_Fy2z_S_D2x_S_C1002000003;
  abcd[249] = I_ERI_F3z_S_F3x_S_C1002000003+CDX*I_ERI_F3z_S_D2x_S_C1002000003;
  abcd[280] = I_ERI_F3x_S_F2xy_S_C1002000003+CDX*I_ERI_F3x_S_Dxy_S_C1002000003;
  abcd[281] = I_ERI_F2xy_S_F2xy_S_C1002000003+CDX*I_ERI_F2xy_S_Dxy_S_C1002000003;
  abcd[282] = I_ERI_F2xz_S_F2xy_S_C1002000003+CDX*I_ERI_F2xz_S_Dxy_S_C1002000003;
  abcd[283] = I_ERI_Fx2y_S_F2xy_S_C1002000003+CDX*I_ERI_Fx2y_S_Dxy_S_C1002000003;
  abcd[284] = I_ERI_Fxyz_S_F2xy_S_C1002000003+CDX*I_ERI_Fxyz_S_Dxy_S_C1002000003;
  abcd[285] = I_ERI_Fx2z_S_F2xy_S_C1002000003+CDX*I_ERI_Fx2z_S_Dxy_S_C1002000003;
  abcd[286] = I_ERI_F3y_S_F2xy_S_C1002000003+CDX*I_ERI_F3y_S_Dxy_S_C1002000003;
  abcd[287] = I_ERI_F2yz_S_F2xy_S_C1002000003+CDX*I_ERI_F2yz_S_Dxy_S_C1002000003;
  abcd[288] = I_ERI_Fy2z_S_F2xy_S_C1002000003+CDX*I_ERI_Fy2z_S_Dxy_S_C1002000003;
  abcd[289] = I_ERI_F3z_S_F2xy_S_C1002000003+CDX*I_ERI_F3z_S_Dxy_S_C1002000003;
  abcd[320] = I_ERI_F3x_S_F2xz_S_C1002000003+CDX*I_ERI_F3x_S_Dxz_S_C1002000003;
  abcd[321] = I_ERI_F2xy_S_F2xz_S_C1002000003+CDX*I_ERI_F2xy_S_Dxz_S_C1002000003;
  abcd[322] = I_ERI_F2xz_S_F2xz_S_C1002000003+CDX*I_ERI_F2xz_S_Dxz_S_C1002000003;
  abcd[323] = I_ERI_Fx2y_S_F2xz_S_C1002000003+CDX*I_ERI_Fx2y_S_Dxz_S_C1002000003;
  abcd[324] = I_ERI_Fxyz_S_F2xz_S_C1002000003+CDX*I_ERI_Fxyz_S_Dxz_S_C1002000003;
  abcd[325] = I_ERI_Fx2z_S_F2xz_S_C1002000003+CDX*I_ERI_Fx2z_S_Dxz_S_C1002000003;
  abcd[326] = I_ERI_F3y_S_F2xz_S_C1002000003+CDX*I_ERI_F3y_S_Dxz_S_C1002000003;
  abcd[327] = I_ERI_F2yz_S_F2xz_S_C1002000003+CDX*I_ERI_F2yz_S_Dxz_S_C1002000003;
  abcd[328] = I_ERI_Fy2z_S_F2xz_S_C1002000003+CDX*I_ERI_Fy2z_S_Dxz_S_C1002000003;
  abcd[329] = I_ERI_F3z_S_F2xz_S_C1002000003+CDX*I_ERI_F3z_S_Dxz_S_C1002000003;
  abcd[360] = I_ERI_F3x_S_Fx2y_S_C1002000003+CDX*I_ERI_F3x_S_D2y_S_C1002000003;
  abcd[361] = I_ERI_F2xy_S_Fx2y_S_C1002000003+CDX*I_ERI_F2xy_S_D2y_S_C1002000003;
  abcd[362] = I_ERI_F2xz_S_Fx2y_S_C1002000003+CDX*I_ERI_F2xz_S_D2y_S_C1002000003;
  abcd[363] = I_ERI_Fx2y_S_Fx2y_S_C1002000003+CDX*I_ERI_Fx2y_S_D2y_S_C1002000003;
  abcd[364] = I_ERI_Fxyz_S_Fx2y_S_C1002000003+CDX*I_ERI_Fxyz_S_D2y_S_C1002000003;
  abcd[365] = I_ERI_Fx2z_S_Fx2y_S_C1002000003+CDX*I_ERI_Fx2z_S_D2y_S_C1002000003;
  abcd[366] = I_ERI_F3y_S_Fx2y_S_C1002000003+CDX*I_ERI_F3y_S_D2y_S_C1002000003;
  abcd[367] = I_ERI_F2yz_S_Fx2y_S_C1002000003+CDX*I_ERI_F2yz_S_D2y_S_C1002000003;
  abcd[368] = I_ERI_Fy2z_S_Fx2y_S_C1002000003+CDX*I_ERI_Fy2z_S_D2y_S_C1002000003;
  abcd[369] = I_ERI_F3z_S_Fx2y_S_C1002000003+CDX*I_ERI_F3z_S_D2y_S_C1002000003;
  abcd[400] = I_ERI_F3x_S_Fxyz_S_C1002000003+CDX*I_ERI_F3x_S_Dyz_S_C1002000003;
  abcd[401] = I_ERI_F2xy_S_Fxyz_S_C1002000003+CDX*I_ERI_F2xy_S_Dyz_S_C1002000003;
  abcd[402] = I_ERI_F2xz_S_Fxyz_S_C1002000003+CDX*I_ERI_F2xz_S_Dyz_S_C1002000003;
  abcd[403] = I_ERI_Fx2y_S_Fxyz_S_C1002000003+CDX*I_ERI_Fx2y_S_Dyz_S_C1002000003;
  abcd[404] = I_ERI_Fxyz_S_Fxyz_S_C1002000003+CDX*I_ERI_Fxyz_S_Dyz_S_C1002000003;
  abcd[405] = I_ERI_Fx2z_S_Fxyz_S_C1002000003+CDX*I_ERI_Fx2z_S_Dyz_S_C1002000003;
  abcd[406] = I_ERI_F3y_S_Fxyz_S_C1002000003+CDX*I_ERI_F3y_S_Dyz_S_C1002000003;
  abcd[407] = I_ERI_F2yz_S_Fxyz_S_C1002000003+CDX*I_ERI_F2yz_S_Dyz_S_C1002000003;
  abcd[408] = I_ERI_Fy2z_S_Fxyz_S_C1002000003+CDX*I_ERI_Fy2z_S_Dyz_S_C1002000003;
  abcd[409] = I_ERI_F3z_S_Fxyz_S_C1002000003+CDX*I_ERI_F3z_S_Dyz_S_C1002000003;
  abcd[440] = I_ERI_F3x_S_Fx2z_S_C1002000003+CDX*I_ERI_F3x_S_D2z_S_C1002000003;
  abcd[441] = I_ERI_F2xy_S_Fx2z_S_C1002000003+CDX*I_ERI_F2xy_S_D2z_S_C1002000003;
  abcd[442] = I_ERI_F2xz_S_Fx2z_S_C1002000003+CDX*I_ERI_F2xz_S_D2z_S_C1002000003;
  abcd[443] = I_ERI_Fx2y_S_Fx2z_S_C1002000003+CDX*I_ERI_Fx2y_S_D2z_S_C1002000003;
  abcd[444] = I_ERI_Fxyz_S_Fx2z_S_C1002000003+CDX*I_ERI_Fxyz_S_D2z_S_C1002000003;
  abcd[445] = I_ERI_Fx2z_S_Fx2z_S_C1002000003+CDX*I_ERI_Fx2z_S_D2z_S_C1002000003;
  abcd[446] = I_ERI_F3y_S_Fx2z_S_C1002000003+CDX*I_ERI_F3y_S_D2z_S_C1002000003;
  abcd[447] = I_ERI_F2yz_S_Fx2z_S_C1002000003+CDX*I_ERI_F2yz_S_D2z_S_C1002000003;
  abcd[448] = I_ERI_Fy2z_S_Fx2z_S_C1002000003+CDX*I_ERI_Fy2z_S_D2z_S_C1002000003;
  abcd[449] = I_ERI_F3z_S_Fx2z_S_C1002000003+CDX*I_ERI_F3z_S_D2z_S_C1002000003;
  abcd[480] = I_ERI_F3x_S_F2xy_S_C1002000003+CDY*I_ERI_F3x_S_D2x_S_C1002000003;
  abcd[481] = I_ERI_F2xy_S_F2xy_S_C1002000003+CDY*I_ERI_F2xy_S_D2x_S_C1002000003;
  abcd[482] = I_ERI_F2xz_S_F2xy_S_C1002000003+CDY*I_ERI_F2xz_S_D2x_S_C1002000003;
  abcd[483] = I_ERI_Fx2y_S_F2xy_S_C1002000003+CDY*I_ERI_Fx2y_S_D2x_S_C1002000003;
  abcd[484] = I_ERI_Fxyz_S_F2xy_S_C1002000003+CDY*I_ERI_Fxyz_S_D2x_S_C1002000003;
  abcd[485] = I_ERI_Fx2z_S_F2xy_S_C1002000003+CDY*I_ERI_Fx2z_S_D2x_S_C1002000003;
  abcd[486] = I_ERI_F3y_S_F2xy_S_C1002000003+CDY*I_ERI_F3y_S_D2x_S_C1002000003;
  abcd[487] = I_ERI_F2yz_S_F2xy_S_C1002000003+CDY*I_ERI_F2yz_S_D2x_S_C1002000003;
  abcd[488] = I_ERI_Fy2z_S_F2xy_S_C1002000003+CDY*I_ERI_Fy2z_S_D2x_S_C1002000003;
  abcd[489] = I_ERI_F3z_S_F2xy_S_C1002000003+CDY*I_ERI_F3z_S_D2x_S_C1002000003;
  abcd[520] = I_ERI_F3x_S_Fx2y_S_C1002000003+CDY*I_ERI_F3x_S_Dxy_S_C1002000003;
  abcd[521] = I_ERI_F2xy_S_Fx2y_S_C1002000003+CDY*I_ERI_F2xy_S_Dxy_S_C1002000003;
  abcd[522] = I_ERI_F2xz_S_Fx2y_S_C1002000003+CDY*I_ERI_F2xz_S_Dxy_S_C1002000003;
  abcd[523] = I_ERI_Fx2y_S_Fx2y_S_C1002000003+CDY*I_ERI_Fx2y_S_Dxy_S_C1002000003;
  abcd[524] = I_ERI_Fxyz_S_Fx2y_S_C1002000003+CDY*I_ERI_Fxyz_S_Dxy_S_C1002000003;
  abcd[525] = I_ERI_Fx2z_S_Fx2y_S_C1002000003+CDY*I_ERI_Fx2z_S_Dxy_S_C1002000003;
  abcd[526] = I_ERI_F3y_S_Fx2y_S_C1002000003+CDY*I_ERI_F3y_S_Dxy_S_C1002000003;
  abcd[527] = I_ERI_F2yz_S_Fx2y_S_C1002000003+CDY*I_ERI_F2yz_S_Dxy_S_C1002000003;
  abcd[528] = I_ERI_Fy2z_S_Fx2y_S_C1002000003+CDY*I_ERI_Fy2z_S_Dxy_S_C1002000003;
  abcd[529] = I_ERI_F3z_S_Fx2y_S_C1002000003+CDY*I_ERI_F3z_S_Dxy_S_C1002000003;
  abcd[560] = I_ERI_F3x_S_Fxyz_S_C1002000003+CDY*I_ERI_F3x_S_Dxz_S_C1002000003;
  abcd[561] = I_ERI_F2xy_S_Fxyz_S_C1002000003+CDY*I_ERI_F2xy_S_Dxz_S_C1002000003;
  abcd[562] = I_ERI_F2xz_S_Fxyz_S_C1002000003+CDY*I_ERI_F2xz_S_Dxz_S_C1002000003;
  abcd[563] = I_ERI_Fx2y_S_Fxyz_S_C1002000003+CDY*I_ERI_Fx2y_S_Dxz_S_C1002000003;
  abcd[564] = I_ERI_Fxyz_S_Fxyz_S_C1002000003+CDY*I_ERI_Fxyz_S_Dxz_S_C1002000003;
  abcd[565] = I_ERI_Fx2z_S_Fxyz_S_C1002000003+CDY*I_ERI_Fx2z_S_Dxz_S_C1002000003;
  abcd[566] = I_ERI_F3y_S_Fxyz_S_C1002000003+CDY*I_ERI_F3y_S_Dxz_S_C1002000003;
  abcd[567] = I_ERI_F2yz_S_Fxyz_S_C1002000003+CDY*I_ERI_F2yz_S_Dxz_S_C1002000003;
  abcd[568] = I_ERI_Fy2z_S_Fxyz_S_C1002000003+CDY*I_ERI_Fy2z_S_Dxz_S_C1002000003;
  abcd[569] = I_ERI_F3z_S_Fxyz_S_C1002000003+CDY*I_ERI_F3z_S_Dxz_S_C1002000003;
  abcd[600] = I_ERI_F3x_S_F3y_S_C1002000003+CDY*I_ERI_F3x_S_D2y_S_C1002000003;
  abcd[601] = I_ERI_F2xy_S_F3y_S_C1002000003+CDY*I_ERI_F2xy_S_D2y_S_C1002000003;
  abcd[602] = I_ERI_F2xz_S_F3y_S_C1002000003+CDY*I_ERI_F2xz_S_D2y_S_C1002000003;
  abcd[603] = I_ERI_Fx2y_S_F3y_S_C1002000003+CDY*I_ERI_Fx2y_S_D2y_S_C1002000003;
  abcd[604] = I_ERI_Fxyz_S_F3y_S_C1002000003+CDY*I_ERI_Fxyz_S_D2y_S_C1002000003;
  abcd[605] = I_ERI_Fx2z_S_F3y_S_C1002000003+CDY*I_ERI_Fx2z_S_D2y_S_C1002000003;
  abcd[606] = I_ERI_F3y_S_F3y_S_C1002000003+CDY*I_ERI_F3y_S_D2y_S_C1002000003;
  abcd[607] = I_ERI_F2yz_S_F3y_S_C1002000003+CDY*I_ERI_F2yz_S_D2y_S_C1002000003;
  abcd[608] = I_ERI_Fy2z_S_F3y_S_C1002000003+CDY*I_ERI_Fy2z_S_D2y_S_C1002000003;
  abcd[609] = I_ERI_F3z_S_F3y_S_C1002000003+CDY*I_ERI_F3z_S_D2y_S_C1002000003;
  abcd[640] = I_ERI_F3x_S_F2yz_S_C1002000003+CDY*I_ERI_F3x_S_Dyz_S_C1002000003;
  abcd[641] = I_ERI_F2xy_S_F2yz_S_C1002000003+CDY*I_ERI_F2xy_S_Dyz_S_C1002000003;
  abcd[642] = I_ERI_F2xz_S_F2yz_S_C1002000003+CDY*I_ERI_F2xz_S_Dyz_S_C1002000003;
  abcd[643] = I_ERI_Fx2y_S_F2yz_S_C1002000003+CDY*I_ERI_Fx2y_S_Dyz_S_C1002000003;
  abcd[644] = I_ERI_Fxyz_S_F2yz_S_C1002000003+CDY*I_ERI_Fxyz_S_Dyz_S_C1002000003;
  abcd[645] = I_ERI_Fx2z_S_F2yz_S_C1002000003+CDY*I_ERI_Fx2z_S_Dyz_S_C1002000003;
  abcd[646] = I_ERI_F3y_S_F2yz_S_C1002000003+CDY*I_ERI_F3y_S_Dyz_S_C1002000003;
  abcd[647] = I_ERI_F2yz_S_F2yz_S_C1002000003+CDY*I_ERI_F2yz_S_Dyz_S_C1002000003;
  abcd[648] = I_ERI_Fy2z_S_F2yz_S_C1002000003+CDY*I_ERI_Fy2z_S_Dyz_S_C1002000003;
  abcd[649] = I_ERI_F3z_S_F2yz_S_C1002000003+CDY*I_ERI_F3z_S_Dyz_S_C1002000003;
  abcd[680] = I_ERI_F3x_S_Fy2z_S_C1002000003+CDY*I_ERI_F3x_S_D2z_S_C1002000003;
  abcd[681] = I_ERI_F2xy_S_Fy2z_S_C1002000003+CDY*I_ERI_F2xy_S_D2z_S_C1002000003;
  abcd[682] = I_ERI_F2xz_S_Fy2z_S_C1002000003+CDY*I_ERI_F2xz_S_D2z_S_C1002000003;
  abcd[683] = I_ERI_Fx2y_S_Fy2z_S_C1002000003+CDY*I_ERI_Fx2y_S_D2z_S_C1002000003;
  abcd[684] = I_ERI_Fxyz_S_Fy2z_S_C1002000003+CDY*I_ERI_Fxyz_S_D2z_S_C1002000003;
  abcd[685] = I_ERI_Fx2z_S_Fy2z_S_C1002000003+CDY*I_ERI_Fx2z_S_D2z_S_C1002000003;
  abcd[686] = I_ERI_F3y_S_Fy2z_S_C1002000003+CDY*I_ERI_F3y_S_D2z_S_C1002000003;
  abcd[687] = I_ERI_F2yz_S_Fy2z_S_C1002000003+CDY*I_ERI_F2yz_S_D2z_S_C1002000003;
  abcd[688] = I_ERI_Fy2z_S_Fy2z_S_C1002000003+CDY*I_ERI_Fy2z_S_D2z_S_C1002000003;
  abcd[689] = I_ERI_F3z_S_Fy2z_S_C1002000003+CDY*I_ERI_F3z_S_D2z_S_C1002000003;
  abcd[720] = I_ERI_F3x_S_F2xz_S_C1002000003+CDZ*I_ERI_F3x_S_D2x_S_C1002000003;
  abcd[721] = I_ERI_F2xy_S_F2xz_S_C1002000003+CDZ*I_ERI_F2xy_S_D2x_S_C1002000003;
  abcd[722] = I_ERI_F2xz_S_F2xz_S_C1002000003+CDZ*I_ERI_F2xz_S_D2x_S_C1002000003;
  abcd[723] = I_ERI_Fx2y_S_F2xz_S_C1002000003+CDZ*I_ERI_Fx2y_S_D2x_S_C1002000003;
  abcd[724] = I_ERI_Fxyz_S_F2xz_S_C1002000003+CDZ*I_ERI_Fxyz_S_D2x_S_C1002000003;
  abcd[725] = I_ERI_Fx2z_S_F2xz_S_C1002000003+CDZ*I_ERI_Fx2z_S_D2x_S_C1002000003;
  abcd[726] = I_ERI_F3y_S_F2xz_S_C1002000003+CDZ*I_ERI_F3y_S_D2x_S_C1002000003;
  abcd[727] = I_ERI_F2yz_S_F2xz_S_C1002000003+CDZ*I_ERI_F2yz_S_D2x_S_C1002000003;
  abcd[728] = I_ERI_Fy2z_S_F2xz_S_C1002000003+CDZ*I_ERI_Fy2z_S_D2x_S_C1002000003;
  abcd[729] = I_ERI_F3z_S_F2xz_S_C1002000003+CDZ*I_ERI_F3z_S_D2x_S_C1002000003;
  abcd[760] = I_ERI_F3x_S_Fxyz_S_C1002000003+CDZ*I_ERI_F3x_S_Dxy_S_C1002000003;
  abcd[761] = I_ERI_F2xy_S_Fxyz_S_C1002000003+CDZ*I_ERI_F2xy_S_Dxy_S_C1002000003;
  abcd[762] = I_ERI_F2xz_S_Fxyz_S_C1002000003+CDZ*I_ERI_F2xz_S_Dxy_S_C1002000003;
  abcd[763] = I_ERI_Fx2y_S_Fxyz_S_C1002000003+CDZ*I_ERI_Fx2y_S_Dxy_S_C1002000003;
  abcd[764] = I_ERI_Fxyz_S_Fxyz_S_C1002000003+CDZ*I_ERI_Fxyz_S_Dxy_S_C1002000003;
  abcd[765] = I_ERI_Fx2z_S_Fxyz_S_C1002000003+CDZ*I_ERI_Fx2z_S_Dxy_S_C1002000003;
  abcd[766] = I_ERI_F3y_S_Fxyz_S_C1002000003+CDZ*I_ERI_F3y_S_Dxy_S_C1002000003;
  abcd[767] = I_ERI_F2yz_S_Fxyz_S_C1002000003+CDZ*I_ERI_F2yz_S_Dxy_S_C1002000003;
  abcd[768] = I_ERI_Fy2z_S_Fxyz_S_C1002000003+CDZ*I_ERI_Fy2z_S_Dxy_S_C1002000003;
  abcd[769] = I_ERI_F3z_S_Fxyz_S_C1002000003+CDZ*I_ERI_F3z_S_Dxy_S_C1002000003;
  abcd[800] = I_ERI_F3x_S_Fx2z_S_C1002000003+CDZ*I_ERI_F3x_S_Dxz_S_C1002000003;
  abcd[801] = I_ERI_F2xy_S_Fx2z_S_C1002000003+CDZ*I_ERI_F2xy_S_Dxz_S_C1002000003;
  abcd[802] = I_ERI_F2xz_S_Fx2z_S_C1002000003+CDZ*I_ERI_F2xz_S_Dxz_S_C1002000003;
  abcd[803] = I_ERI_Fx2y_S_Fx2z_S_C1002000003+CDZ*I_ERI_Fx2y_S_Dxz_S_C1002000003;
  abcd[804] = I_ERI_Fxyz_S_Fx2z_S_C1002000003+CDZ*I_ERI_Fxyz_S_Dxz_S_C1002000003;
  abcd[805] = I_ERI_Fx2z_S_Fx2z_S_C1002000003+CDZ*I_ERI_Fx2z_S_Dxz_S_C1002000003;
  abcd[806] = I_ERI_F3y_S_Fx2z_S_C1002000003+CDZ*I_ERI_F3y_S_Dxz_S_C1002000003;
  abcd[807] = I_ERI_F2yz_S_Fx2z_S_C1002000003+CDZ*I_ERI_F2yz_S_Dxz_S_C1002000003;
  abcd[808] = I_ERI_Fy2z_S_Fx2z_S_C1002000003+CDZ*I_ERI_Fy2z_S_Dxz_S_C1002000003;
  abcd[809] = I_ERI_F3z_S_Fx2z_S_C1002000003+CDZ*I_ERI_F3z_S_Dxz_S_C1002000003;
  abcd[840] = I_ERI_F3x_S_F2yz_S_C1002000003+CDZ*I_ERI_F3x_S_D2y_S_C1002000003;
  abcd[841] = I_ERI_F2xy_S_F2yz_S_C1002000003+CDZ*I_ERI_F2xy_S_D2y_S_C1002000003;
  abcd[842] = I_ERI_F2xz_S_F2yz_S_C1002000003+CDZ*I_ERI_F2xz_S_D2y_S_C1002000003;
  abcd[843] = I_ERI_Fx2y_S_F2yz_S_C1002000003+CDZ*I_ERI_Fx2y_S_D2y_S_C1002000003;
  abcd[844] = I_ERI_Fxyz_S_F2yz_S_C1002000003+CDZ*I_ERI_Fxyz_S_D2y_S_C1002000003;
  abcd[845] = I_ERI_Fx2z_S_F2yz_S_C1002000003+CDZ*I_ERI_Fx2z_S_D2y_S_C1002000003;
  abcd[846] = I_ERI_F3y_S_F2yz_S_C1002000003+CDZ*I_ERI_F3y_S_D2y_S_C1002000003;
  abcd[847] = I_ERI_F2yz_S_F2yz_S_C1002000003+CDZ*I_ERI_F2yz_S_D2y_S_C1002000003;
  abcd[848] = I_ERI_Fy2z_S_F2yz_S_C1002000003+CDZ*I_ERI_Fy2z_S_D2y_S_C1002000003;
  abcd[849] = I_ERI_F3z_S_F2yz_S_C1002000003+CDZ*I_ERI_F3z_S_D2y_S_C1002000003;
  abcd[880] = I_ERI_F3x_S_Fy2z_S_C1002000003+CDZ*I_ERI_F3x_S_Dyz_S_C1002000003;
  abcd[881] = I_ERI_F2xy_S_Fy2z_S_C1002000003+CDZ*I_ERI_F2xy_S_Dyz_S_C1002000003;
  abcd[882] = I_ERI_F2xz_S_Fy2z_S_C1002000003+CDZ*I_ERI_F2xz_S_Dyz_S_C1002000003;
  abcd[883] = I_ERI_Fx2y_S_Fy2z_S_C1002000003+CDZ*I_ERI_Fx2y_S_Dyz_S_C1002000003;
  abcd[884] = I_ERI_Fxyz_S_Fy2z_S_C1002000003+CDZ*I_ERI_Fxyz_S_Dyz_S_C1002000003;
  abcd[885] = I_ERI_Fx2z_S_Fy2z_S_C1002000003+CDZ*I_ERI_Fx2z_S_Dyz_S_C1002000003;
  abcd[886] = I_ERI_F3y_S_Fy2z_S_C1002000003+CDZ*I_ERI_F3y_S_Dyz_S_C1002000003;
  abcd[887] = I_ERI_F2yz_S_Fy2z_S_C1002000003+CDZ*I_ERI_F2yz_S_Dyz_S_C1002000003;
  abcd[888] = I_ERI_Fy2z_S_Fy2z_S_C1002000003+CDZ*I_ERI_Fy2z_S_Dyz_S_C1002000003;
  abcd[889] = I_ERI_F3z_S_Fy2z_S_C1002000003+CDZ*I_ERI_F3z_S_Dyz_S_C1002000003;
  abcd[920] = I_ERI_F3x_S_F3z_S_C1002000003+CDZ*I_ERI_F3x_S_D2z_S_C1002000003;
  abcd[921] = I_ERI_F2xy_S_F3z_S_C1002000003+CDZ*I_ERI_F2xy_S_D2z_S_C1002000003;
  abcd[922] = I_ERI_F2xz_S_F3z_S_C1002000003+CDZ*I_ERI_F2xz_S_D2z_S_C1002000003;
  abcd[923] = I_ERI_Fx2y_S_F3z_S_C1002000003+CDZ*I_ERI_Fx2y_S_D2z_S_C1002000003;
  abcd[924] = I_ERI_Fxyz_S_F3z_S_C1002000003+CDZ*I_ERI_Fxyz_S_D2z_S_C1002000003;
  abcd[925] = I_ERI_Fx2z_S_F3z_S_C1002000003+CDZ*I_ERI_Fx2z_S_D2z_S_C1002000003;
  abcd[926] = I_ERI_F3y_S_F3z_S_C1002000003+CDZ*I_ERI_F3y_S_D2z_S_C1002000003;
  abcd[927] = I_ERI_F2yz_S_F3z_S_C1002000003+CDZ*I_ERI_F2yz_S_D2z_S_C1002000003;
  abcd[928] = I_ERI_Fy2z_S_F3z_S_C1002000003+CDZ*I_ERI_Fy2z_S_D2z_S_C1002000003;
  abcd[929] = I_ERI_F3z_S_F3z_S_C1002000003+CDZ*I_ERI_F3z_S_D2z_S_C1002000003;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_P_C1002001003
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_C1002001003
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1002001003
   ************************************************************/
  Double I_ERI_F3x_S_D2x_Px_C1002001003 = I_ERI_F3x_S_F3x_S_C1002001003+CDX*I_ERI_F3x_S_D2x_S_C1002001003;
  Double I_ERI_F2xy_S_D2x_Px_C1002001003 = I_ERI_F2xy_S_F3x_S_C1002001003+CDX*I_ERI_F2xy_S_D2x_S_C1002001003;
  Double I_ERI_F2xz_S_D2x_Px_C1002001003 = I_ERI_F2xz_S_F3x_S_C1002001003+CDX*I_ERI_F2xz_S_D2x_S_C1002001003;
  Double I_ERI_Fx2y_S_D2x_Px_C1002001003 = I_ERI_Fx2y_S_F3x_S_C1002001003+CDX*I_ERI_Fx2y_S_D2x_S_C1002001003;
  Double I_ERI_Fxyz_S_D2x_Px_C1002001003 = I_ERI_Fxyz_S_F3x_S_C1002001003+CDX*I_ERI_Fxyz_S_D2x_S_C1002001003;
  Double I_ERI_Fx2z_S_D2x_Px_C1002001003 = I_ERI_Fx2z_S_F3x_S_C1002001003+CDX*I_ERI_Fx2z_S_D2x_S_C1002001003;
  Double I_ERI_F3y_S_D2x_Px_C1002001003 = I_ERI_F3y_S_F3x_S_C1002001003+CDX*I_ERI_F3y_S_D2x_S_C1002001003;
  Double I_ERI_F2yz_S_D2x_Px_C1002001003 = I_ERI_F2yz_S_F3x_S_C1002001003+CDX*I_ERI_F2yz_S_D2x_S_C1002001003;
  Double I_ERI_Fy2z_S_D2x_Px_C1002001003 = I_ERI_Fy2z_S_F3x_S_C1002001003+CDX*I_ERI_Fy2z_S_D2x_S_C1002001003;
  Double I_ERI_F3z_S_D2x_Px_C1002001003 = I_ERI_F3z_S_F3x_S_C1002001003+CDX*I_ERI_F3z_S_D2x_S_C1002001003;
  Double I_ERI_F3x_S_Dxy_Px_C1002001003 = I_ERI_F3x_S_F2xy_S_C1002001003+CDX*I_ERI_F3x_S_Dxy_S_C1002001003;
  Double I_ERI_F2xy_S_Dxy_Px_C1002001003 = I_ERI_F2xy_S_F2xy_S_C1002001003+CDX*I_ERI_F2xy_S_Dxy_S_C1002001003;
  Double I_ERI_F2xz_S_Dxy_Px_C1002001003 = I_ERI_F2xz_S_F2xy_S_C1002001003+CDX*I_ERI_F2xz_S_Dxy_S_C1002001003;
  Double I_ERI_Fx2y_S_Dxy_Px_C1002001003 = I_ERI_Fx2y_S_F2xy_S_C1002001003+CDX*I_ERI_Fx2y_S_Dxy_S_C1002001003;
  Double I_ERI_Fxyz_S_Dxy_Px_C1002001003 = I_ERI_Fxyz_S_F2xy_S_C1002001003+CDX*I_ERI_Fxyz_S_Dxy_S_C1002001003;
  Double I_ERI_Fx2z_S_Dxy_Px_C1002001003 = I_ERI_Fx2z_S_F2xy_S_C1002001003+CDX*I_ERI_Fx2z_S_Dxy_S_C1002001003;
  Double I_ERI_F3y_S_Dxy_Px_C1002001003 = I_ERI_F3y_S_F2xy_S_C1002001003+CDX*I_ERI_F3y_S_Dxy_S_C1002001003;
  Double I_ERI_F2yz_S_Dxy_Px_C1002001003 = I_ERI_F2yz_S_F2xy_S_C1002001003+CDX*I_ERI_F2yz_S_Dxy_S_C1002001003;
  Double I_ERI_Fy2z_S_Dxy_Px_C1002001003 = I_ERI_Fy2z_S_F2xy_S_C1002001003+CDX*I_ERI_Fy2z_S_Dxy_S_C1002001003;
  Double I_ERI_F3z_S_Dxy_Px_C1002001003 = I_ERI_F3z_S_F2xy_S_C1002001003+CDX*I_ERI_F3z_S_Dxy_S_C1002001003;
  Double I_ERI_F3x_S_Dxz_Px_C1002001003 = I_ERI_F3x_S_F2xz_S_C1002001003+CDX*I_ERI_F3x_S_Dxz_S_C1002001003;
  Double I_ERI_F2xy_S_Dxz_Px_C1002001003 = I_ERI_F2xy_S_F2xz_S_C1002001003+CDX*I_ERI_F2xy_S_Dxz_S_C1002001003;
  Double I_ERI_F2xz_S_Dxz_Px_C1002001003 = I_ERI_F2xz_S_F2xz_S_C1002001003+CDX*I_ERI_F2xz_S_Dxz_S_C1002001003;
  Double I_ERI_Fx2y_S_Dxz_Px_C1002001003 = I_ERI_Fx2y_S_F2xz_S_C1002001003+CDX*I_ERI_Fx2y_S_Dxz_S_C1002001003;
  Double I_ERI_Fxyz_S_Dxz_Px_C1002001003 = I_ERI_Fxyz_S_F2xz_S_C1002001003+CDX*I_ERI_Fxyz_S_Dxz_S_C1002001003;
  Double I_ERI_Fx2z_S_Dxz_Px_C1002001003 = I_ERI_Fx2z_S_F2xz_S_C1002001003+CDX*I_ERI_Fx2z_S_Dxz_S_C1002001003;
  Double I_ERI_F3y_S_Dxz_Px_C1002001003 = I_ERI_F3y_S_F2xz_S_C1002001003+CDX*I_ERI_F3y_S_Dxz_S_C1002001003;
  Double I_ERI_F2yz_S_Dxz_Px_C1002001003 = I_ERI_F2yz_S_F2xz_S_C1002001003+CDX*I_ERI_F2yz_S_Dxz_S_C1002001003;
  Double I_ERI_Fy2z_S_Dxz_Px_C1002001003 = I_ERI_Fy2z_S_F2xz_S_C1002001003+CDX*I_ERI_Fy2z_S_Dxz_S_C1002001003;
  Double I_ERI_F3z_S_Dxz_Px_C1002001003 = I_ERI_F3z_S_F2xz_S_C1002001003+CDX*I_ERI_F3z_S_Dxz_S_C1002001003;
  Double I_ERI_F3x_S_D2y_Px_C1002001003 = I_ERI_F3x_S_Fx2y_S_C1002001003+CDX*I_ERI_F3x_S_D2y_S_C1002001003;
  Double I_ERI_F2xy_S_D2y_Px_C1002001003 = I_ERI_F2xy_S_Fx2y_S_C1002001003+CDX*I_ERI_F2xy_S_D2y_S_C1002001003;
  Double I_ERI_F2xz_S_D2y_Px_C1002001003 = I_ERI_F2xz_S_Fx2y_S_C1002001003+CDX*I_ERI_F2xz_S_D2y_S_C1002001003;
  Double I_ERI_Fx2y_S_D2y_Px_C1002001003 = I_ERI_Fx2y_S_Fx2y_S_C1002001003+CDX*I_ERI_Fx2y_S_D2y_S_C1002001003;
  Double I_ERI_Fxyz_S_D2y_Px_C1002001003 = I_ERI_Fxyz_S_Fx2y_S_C1002001003+CDX*I_ERI_Fxyz_S_D2y_S_C1002001003;
  Double I_ERI_Fx2z_S_D2y_Px_C1002001003 = I_ERI_Fx2z_S_Fx2y_S_C1002001003+CDX*I_ERI_Fx2z_S_D2y_S_C1002001003;
  Double I_ERI_F3y_S_D2y_Px_C1002001003 = I_ERI_F3y_S_Fx2y_S_C1002001003+CDX*I_ERI_F3y_S_D2y_S_C1002001003;
  Double I_ERI_F2yz_S_D2y_Px_C1002001003 = I_ERI_F2yz_S_Fx2y_S_C1002001003+CDX*I_ERI_F2yz_S_D2y_S_C1002001003;
  Double I_ERI_Fy2z_S_D2y_Px_C1002001003 = I_ERI_Fy2z_S_Fx2y_S_C1002001003+CDX*I_ERI_Fy2z_S_D2y_S_C1002001003;
  Double I_ERI_F3z_S_D2y_Px_C1002001003 = I_ERI_F3z_S_Fx2y_S_C1002001003+CDX*I_ERI_F3z_S_D2y_S_C1002001003;
  Double I_ERI_F3x_S_Dyz_Px_C1002001003 = I_ERI_F3x_S_Fxyz_S_C1002001003+CDX*I_ERI_F3x_S_Dyz_S_C1002001003;
  Double I_ERI_F2xy_S_Dyz_Px_C1002001003 = I_ERI_F2xy_S_Fxyz_S_C1002001003+CDX*I_ERI_F2xy_S_Dyz_S_C1002001003;
  Double I_ERI_F2xz_S_Dyz_Px_C1002001003 = I_ERI_F2xz_S_Fxyz_S_C1002001003+CDX*I_ERI_F2xz_S_Dyz_S_C1002001003;
  Double I_ERI_Fx2y_S_Dyz_Px_C1002001003 = I_ERI_Fx2y_S_Fxyz_S_C1002001003+CDX*I_ERI_Fx2y_S_Dyz_S_C1002001003;
  Double I_ERI_Fxyz_S_Dyz_Px_C1002001003 = I_ERI_Fxyz_S_Fxyz_S_C1002001003+CDX*I_ERI_Fxyz_S_Dyz_S_C1002001003;
  Double I_ERI_Fx2z_S_Dyz_Px_C1002001003 = I_ERI_Fx2z_S_Fxyz_S_C1002001003+CDX*I_ERI_Fx2z_S_Dyz_S_C1002001003;
  Double I_ERI_F3y_S_Dyz_Px_C1002001003 = I_ERI_F3y_S_Fxyz_S_C1002001003+CDX*I_ERI_F3y_S_Dyz_S_C1002001003;
  Double I_ERI_F2yz_S_Dyz_Px_C1002001003 = I_ERI_F2yz_S_Fxyz_S_C1002001003+CDX*I_ERI_F2yz_S_Dyz_S_C1002001003;
  Double I_ERI_Fy2z_S_Dyz_Px_C1002001003 = I_ERI_Fy2z_S_Fxyz_S_C1002001003+CDX*I_ERI_Fy2z_S_Dyz_S_C1002001003;
  Double I_ERI_F3z_S_Dyz_Px_C1002001003 = I_ERI_F3z_S_Fxyz_S_C1002001003+CDX*I_ERI_F3z_S_Dyz_S_C1002001003;
  Double I_ERI_F3x_S_D2z_Px_C1002001003 = I_ERI_F3x_S_Fx2z_S_C1002001003+CDX*I_ERI_F3x_S_D2z_S_C1002001003;
  Double I_ERI_F2xy_S_D2z_Px_C1002001003 = I_ERI_F2xy_S_Fx2z_S_C1002001003+CDX*I_ERI_F2xy_S_D2z_S_C1002001003;
  Double I_ERI_F2xz_S_D2z_Px_C1002001003 = I_ERI_F2xz_S_Fx2z_S_C1002001003+CDX*I_ERI_F2xz_S_D2z_S_C1002001003;
  Double I_ERI_Fx2y_S_D2z_Px_C1002001003 = I_ERI_Fx2y_S_Fx2z_S_C1002001003+CDX*I_ERI_Fx2y_S_D2z_S_C1002001003;
  Double I_ERI_Fxyz_S_D2z_Px_C1002001003 = I_ERI_Fxyz_S_Fx2z_S_C1002001003+CDX*I_ERI_Fxyz_S_D2z_S_C1002001003;
  Double I_ERI_Fx2z_S_D2z_Px_C1002001003 = I_ERI_Fx2z_S_Fx2z_S_C1002001003+CDX*I_ERI_Fx2z_S_D2z_S_C1002001003;
  Double I_ERI_F3y_S_D2z_Px_C1002001003 = I_ERI_F3y_S_Fx2z_S_C1002001003+CDX*I_ERI_F3y_S_D2z_S_C1002001003;
  Double I_ERI_F2yz_S_D2z_Px_C1002001003 = I_ERI_F2yz_S_Fx2z_S_C1002001003+CDX*I_ERI_F2yz_S_D2z_S_C1002001003;
  Double I_ERI_Fy2z_S_D2z_Px_C1002001003 = I_ERI_Fy2z_S_Fx2z_S_C1002001003+CDX*I_ERI_Fy2z_S_D2z_S_C1002001003;
  Double I_ERI_F3z_S_D2z_Px_C1002001003 = I_ERI_F3z_S_Fx2z_S_C1002001003+CDX*I_ERI_F3z_S_D2z_S_C1002001003;
  Double I_ERI_F3x_S_D2x_Py_C1002001003 = I_ERI_F3x_S_F2xy_S_C1002001003+CDY*I_ERI_F3x_S_D2x_S_C1002001003;
  Double I_ERI_F2xy_S_D2x_Py_C1002001003 = I_ERI_F2xy_S_F2xy_S_C1002001003+CDY*I_ERI_F2xy_S_D2x_S_C1002001003;
  Double I_ERI_F2xz_S_D2x_Py_C1002001003 = I_ERI_F2xz_S_F2xy_S_C1002001003+CDY*I_ERI_F2xz_S_D2x_S_C1002001003;
  Double I_ERI_Fx2y_S_D2x_Py_C1002001003 = I_ERI_Fx2y_S_F2xy_S_C1002001003+CDY*I_ERI_Fx2y_S_D2x_S_C1002001003;
  Double I_ERI_Fxyz_S_D2x_Py_C1002001003 = I_ERI_Fxyz_S_F2xy_S_C1002001003+CDY*I_ERI_Fxyz_S_D2x_S_C1002001003;
  Double I_ERI_Fx2z_S_D2x_Py_C1002001003 = I_ERI_Fx2z_S_F2xy_S_C1002001003+CDY*I_ERI_Fx2z_S_D2x_S_C1002001003;
  Double I_ERI_F3y_S_D2x_Py_C1002001003 = I_ERI_F3y_S_F2xy_S_C1002001003+CDY*I_ERI_F3y_S_D2x_S_C1002001003;
  Double I_ERI_F2yz_S_D2x_Py_C1002001003 = I_ERI_F2yz_S_F2xy_S_C1002001003+CDY*I_ERI_F2yz_S_D2x_S_C1002001003;
  Double I_ERI_Fy2z_S_D2x_Py_C1002001003 = I_ERI_Fy2z_S_F2xy_S_C1002001003+CDY*I_ERI_Fy2z_S_D2x_S_C1002001003;
  Double I_ERI_F3z_S_D2x_Py_C1002001003 = I_ERI_F3z_S_F2xy_S_C1002001003+CDY*I_ERI_F3z_S_D2x_S_C1002001003;
  Double I_ERI_F3x_S_Dxy_Py_C1002001003 = I_ERI_F3x_S_Fx2y_S_C1002001003+CDY*I_ERI_F3x_S_Dxy_S_C1002001003;
  Double I_ERI_F2xy_S_Dxy_Py_C1002001003 = I_ERI_F2xy_S_Fx2y_S_C1002001003+CDY*I_ERI_F2xy_S_Dxy_S_C1002001003;
  Double I_ERI_F2xz_S_Dxy_Py_C1002001003 = I_ERI_F2xz_S_Fx2y_S_C1002001003+CDY*I_ERI_F2xz_S_Dxy_S_C1002001003;
  Double I_ERI_Fx2y_S_Dxy_Py_C1002001003 = I_ERI_Fx2y_S_Fx2y_S_C1002001003+CDY*I_ERI_Fx2y_S_Dxy_S_C1002001003;
  Double I_ERI_Fxyz_S_Dxy_Py_C1002001003 = I_ERI_Fxyz_S_Fx2y_S_C1002001003+CDY*I_ERI_Fxyz_S_Dxy_S_C1002001003;
  Double I_ERI_Fx2z_S_Dxy_Py_C1002001003 = I_ERI_Fx2z_S_Fx2y_S_C1002001003+CDY*I_ERI_Fx2z_S_Dxy_S_C1002001003;
  Double I_ERI_F3y_S_Dxy_Py_C1002001003 = I_ERI_F3y_S_Fx2y_S_C1002001003+CDY*I_ERI_F3y_S_Dxy_S_C1002001003;
  Double I_ERI_F2yz_S_Dxy_Py_C1002001003 = I_ERI_F2yz_S_Fx2y_S_C1002001003+CDY*I_ERI_F2yz_S_Dxy_S_C1002001003;
  Double I_ERI_Fy2z_S_Dxy_Py_C1002001003 = I_ERI_Fy2z_S_Fx2y_S_C1002001003+CDY*I_ERI_Fy2z_S_Dxy_S_C1002001003;
  Double I_ERI_F3z_S_Dxy_Py_C1002001003 = I_ERI_F3z_S_Fx2y_S_C1002001003+CDY*I_ERI_F3z_S_Dxy_S_C1002001003;
  Double I_ERI_F3x_S_Dxz_Py_C1002001003 = I_ERI_F3x_S_Fxyz_S_C1002001003+CDY*I_ERI_F3x_S_Dxz_S_C1002001003;
  Double I_ERI_F2xy_S_Dxz_Py_C1002001003 = I_ERI_F2xy_S_Fxyz_S_C1002001003+CDY*I_ERI_F2xy_S_Dxz_S_C1002001003;
  Double I_ERI_F2xz_S_Dxz_Py_C1002001003 = I_ERI_F2xz_S_Fxyz_S_C1002001003+CDY*I_ERI_F2xz_S_Dxz_S_C1002001003;
  Double I_ERI_Fx2y_S_Dxz_Py_C1002001003 = I_ERI_Fx2y_S_Fxyz_S_C1002001003+CDY*I_ERI_Fx2y_S_Dxz_S_C1002001003;
  Double I_ERI_Fxyz_S_Dxz_Py_C1002001003 = I_ERI_Fxyz_S_Fxyz_S_C1002001003+CDY*I_ERI_Fxyz_S_Dxz_S_C1002001003;
  Double I_ERI_Fx2z_S_Dxz_Py_C1002001003 = I_ERI_Fx2z_S_Fxyz_S_C1002001003+CDY*I_ERI_Fx2z_S_Dxz_S_C1002001003;
  Double I_ERI_F3y_S_Dxz_Py_C1002001003 = I_ERI_F3y_S_Fxyz_S_C1002001003+CDY*I_ERI_F3y_S_Dxz_S_C1002001003;
  Double I_ERI_F2yz_S_Dxz_Py_C1002001003 = I_ERI_F2yz_S_Fxyz_S_C1002001003+CDY*I_ERI_F2yz_S_Dxz_S_C1002001003;
  Double I_ERI_Fy2z_S_Dxz_Py_C1002001003 = I_ERI_Fy2z_S_Fxyz_S_C1002001003+CDY*I_ERI_Fy2z_S_Dxz_S_C1002001003;
  Double I_ERI_F3z_S_Dxz_Py_C1002001003 = I_ERI_F3z_S_Fxyz_S_C1002001003+CDY*I_ERI_F3z_S_Dxz_S_C1002001003;
  Double I_ERI_F3x_S_D2y_Py_C1002001003 = I_ERI_F3x_S_F3y_S_C1002001003+CDY*I_ERI_F3x_S_D2y_S_C1002001003;
  Double I_ERI_F2xy_S_D2y_Py_C1002001003 = I_ERI_F2xy_S_F3y_S_C1002001003+CDY*I_ERI_F2xy_S_D2y_S_C1002001003;
  Double I_ERI_F2xz_S_D2y_Py_C1002001003 = I_ERI_F2xz_S_F3y_S_C1002001003+CDY*I_ERI_F2xz_S_D2y_S_C1002001003;
  Double I_ERI_Fx2y_S_D2y_Py_C1002001003 = I_ERI_Fx2y_S_F3y_S_C1002001003+CDY*I_ERI_Fx2y_S_D2y_S_C1002001003;
  Double I_ERI_Fxyz_S_D2y_Py_C1002001003 = I_ERI_Fxyz_S_F3y_S_C1002001003+CDY*I_ERI_Fxyz_S_D2y_S_C1002001003;
  Double I_ERI_Fx2z_S_D2y_Py_C1002001003 = I_ERI_Fx2z_S_F3y_S_C1002001003+CDY*I_ERI_Fx2z_S_D2y_S_C1002001003;
  Double I_ERI_F3y_S_D2y_Py_C1002001003 = I_ERI_F3y_S_F3y_S_C1002001003+CDY*I_ERI_F3y_S_D2y_S_C1002001003;
  Double I_ERI_F2yz_S_D2y_Py_C1002001003 = I_ERI_F2yz_S_F3y_S_C1002001003+CDY*I_ERI_F2yz_S_D2y_S_C1002001003;
  Double I_ERI_Fy2z_S_D2y_Py_C1002001003 = I_ERI_Fy2z_S_F3y_S_C1002001003+CDY*I_ERI_Fy2z_S_D2y_S_C1002001003;
  Double I_ERI_F3z_S_D2y_Py_C1002001003 = I_ERI_F3z_S_F3y_S_C1002001003+CDY*I_ERI_F3z_S_D2y_S_C1002001003;
  Double I_ERI_F3x_S_Dyz_Py_C1002001003 = I_ERI_F3x_S_F2yz_S_C1002001003+CDY*I_ERI_F3x_S_Dyz_S_C1002001003;
  Double I_ERI_F2xy_S_Dyz_Py_C1002001003 = I_ERI_F2xy_S_F2yz_S_C1002001003+CDY*I_ERI_F2xy_S_Dyz_S_C1002001003;
  Double I_ERI_F2xz_S_Dyz_Py_C1002001003 = I_ERI_F2xz_S_F2yz_S_C1002001003+CDY*I_ERI_F2xz_S_Dyz_S_C1002001003;
  Double I_ERI_Fx2y_S_Dyz_Py_C1002001003 = I_ERI_Fx2y_S_F2yz_S_C1002001003+CDY*I_ERI_Fx2y_S_Dyz_S_C1002001003;
  Double I_ERI_Fxyz_S_Dyz_Py_C1002001003 = I_ERI_Fxyz_S_F2yz_S_C1002001003+CDY*I_ERI_Fxyz_S_Dyz_S_C1002001003;
  Double I_ERI_Fx2z_S_Dyz_Py_C1002001003 = I_ERI_Fx2z_S_F2yz_S_C1002001003+CDY*I_ERI_Fx2z_S_Dyz_S_C1002001003;
  Double I_ERI_F3y_S_Dyz_Py_C1002001003 = I_ERI_F3y_S_F2yz_S_C1002001003+CDY*I_ERI_F3y_S_Dyz_S_C1002001003;
  Double I_ERI_F2yz_S_Dyz_Py_C1002001003 = I_ERI_F2yz_S_F2yz_S_C1002001003+CDY*I_ERI_F2yz_S_Dyz_S_C1002001003;
  Double I_ERI_Fy2z_S_Dyz_Py_C1002001003 = I_ERI_Fy2z_S_F2yz_S_C1002001003+CDY*I_ERI_Fy2z_S_Dyz_S_C1002001003;
  Double I_ERI_F3z_S_Dyz_Py_C1002001003 = I_ERI_F3z_S_F2yz_S_C1002001003+CDY*I_ERI_F3z_S_Dyz_S_C1002001003;
  Double I_ERI_F3x_S_D2z_Py_C1002001003 = I_ERI_F3x_S_Fy2z_S_C1002001003+CDY*I_ERI_F3x_S_D2z_S_C1002001003;
  Double I_ERI_F2xy_S_D2z_Py_C1002001003 = I_ERI_F2xy_S_Fy2z_S_C1002001003+CDY*I_ERI_F2xy_S_D2z_S_C1002001003;
  Double I_ERI_F2xz_S_D2z_Py_C1002001003 = I_ERI_F2xz_S_Fy2z_S_C1002001003+CDY*I_ERI_F2xz_S_D2z_S_C1002001003;
  Double I_ERI_Fx2y_S_D2z_Py_C1002001003 = I_ERI_Fx2y_S_Fy2z_S_C1002001003+CDY*I_ERI_Fx2y_S_D2z_S_C1002001003;
  Double I_ERI_Fxyz_S_D2z_Py_C1002001003 = I_ERI_Fxyz_S_Fy2z_S_C1002001003+CDY*I_ERI_Fxyz_S_D2z_S_C1002001003;
  Double I_ERI_Fx2z_S_D2z_Py_C1002001003 = I_ERI_Fx2z_S_Fy2z_S_C1002001003+CDY*I_ERI_Fx2z_S_D2z_S_C1002001003;
  Double I_ERI_F3y_S_D2z_Py_C1002001003 = I_ERI_F3y_S_Fy2z_S_C1002001003+CDY*I_ERI_F3y_S_D2z_S_C1002001003;
  Double I_ERI_F2yz_S_D2z_Py_C1002001003 = I_ERI_F2yz_S_Fy2z_S_C1002001003+CDY*I_ERI_F2yz_S_D2z_S_C1002001003;
  Double I_ERI_Fy2z_S_D2z_Py_C1002001003 = I_ERI_Fy2z_S_Fy2z_S_C1002001003+CDY*I_ERI_Fy2z_S_D2z_S_C1002001003;
  Double I_ERI_F3z_S_D2z_Py_C1002001003 = I_ERI_F3z_S_Fy2z_S_C1002001003+CDY*I_ERI_F3z_S_D2z_S_C1002001003;
  Double I_ERI_F3x_S_D2x_Pz_C1002001003 = I_ERI_F3x_S_F2xz_S_C1002001003+CDZ*I_ERI_F3x_S_D2x_S_C1002001003;
  Double I_ERI_F2xy_S_D2x_Pz_C1002001003 = I_ERI_F2xy_S_F2xz_S_C1002001003+CDZ*I_ERI_F2xy_S_D2x_S_C1002001003;
  Double I_ERI_F2xz_S_D2x_Pz_C1002001003 = I_ERI_F2xz_S_F2xz_S_C1002001003+CDZ*I_ERI_F2xz_S_D2x_S_C1002001003;
  Double I_ERI_Fx2y_S_D2x_Pz_C1002001003 = I_ERI_Fx2y_S_F2xz_S_C1002001003+CDZ*I_ERI_Fx2y_S_D2x_S_C1002001003;
  Double I_ERI_Fxyz_S_D2x_Pz_C1002001003 = I_ERI_Fxyz_S_F2xz_S_C1002001003+CDZ*I_ERI_Fxyz_S_D2x_S_C1002001003;
  Double I_ERI_Fx2z_S_D2x_Pz_C1002001003 = I_ERI_Fx2z_S_F2xz_S_C1002001003+CDZ*I_ERI_Fx2z_S_D2x_S_C1002001003;
  Double I_ERI_F3y_S_D2x_Pz_C1002001003 = I_ERI_F3y_S_F2xz_S_C1002001003+CDZ*I_ERI_F3y_S_D2x_S_C1002001003;
  Double I_ERI_F2yz_S_D2x_Pz_C1002001003 = I_ERI_F2yz_S_F2xz_S_C1002001003+CDZ*I_ERI_F2yz_S_D2x_S_C1002001003;
  Double I_ERI_Fy2z_S_D2x_Pz_C1002001003 = I_ERI_Fy2z_S_F2xz_S_C1002001003+CDZ*I_ERI_Fy2z_S_D2x_S_C1002001003;
  Double I_ERI_F3z_S_D2x_Pz_C1002001003 = I_ERI_F3z_S_F2xz_S_C1002001003+CDZ*I_ERI_F3z_S_D2x_S_C1002001003;
  Double I_ERI_F3x_S_Dxy_Pz_C1002001003 = I_ERI_F3x_S_Fxyz_S_C1002001003+CDZ*I_ERI_F3x_S_Dxy_S_C1002001003;
  Double I_ERI_F2xy_S_Dxy_Pz_C1002001003 = I_ERI_F2xy_S_Fxyz_S_C1002001003+CDZ*I_ERI_F2xy_S_Dxy_S_C1002001003;
  Double I_ERI_F2xz_S_Dxy_Pz_C1002001003 = I_ERI_F2xz_S_Fxyz_S_C1002001003+CDZ*I_ERI_F2xz_S_Dxy_S_C1002001003;
  Double I_ERI_Fx2y_S_Dxy_Pz_C1002001003 = I_ERI_Fx2y_S_Fxyz_S_C1002001003+CDZ*I_ERI_Fx2y_S_Dxy_S_C1002001003;
  Double I_ERI_Fxyz_S_Dxy_Pz_C1002001003 = I_ERI_Fxyz_S_Fxyz_S_C1002001003+CDZ*I_ERI_Fxyz_S_Dxy_S_C1002001003;
  Double I_ERI_Fx2z_S_Dxy_Pz_C1002001003 = I_ERI_Fx2z_S_Fxyz_S_C1002001003+CDZ*I_ERI_Fx2z_S_Dxy_S_C1002001003;
  Double I_ERI_F3y_S_Dxy_Pz_C1002001003 = I_ERI_F3y_S_Fxyz_S_C1002001003+CDZ*I_ERI_F3y_S_Dxy_S_C1002001003;
  Double I_ERI_F2yz_S_Dxy_Pz_C1002001003 = I_ERI_F2yz_S_Fxyz_S_C1002001003+CDZ*I_ERI_F2yz_S_Dxy_S_C1002001003;
  Double I_ERI_Fy2z_S_Dxy_Pz_C1002001003 = I_ERI_Fy2z_S_Fxyz_S_C1002001003+CDZ*I_ERI_Fy2z_S_Dxy_S_C1002001003;
  Double I_ERI_F3z_S_Dxy_Pz_C1002001003 = I_ERI_F3z_S_Fxyz_S_C1002001003+CDZ*I_ERI_F3z_S_Dxy_S_C1002001003;
  Double I_ERI_F3x_S_Dxz_Pz_C1002001003 = I_ERI_F3x_S_Fx2z_S_C1002001003+CDZ*I_ERI_F3x_S_Dxz_S_C1002001003;
  Double I_ERI_F2xy_S_Dxz_Pz_C1002001003 = I_ERI_F2xy_S_Fx2z_S_C1002001003+CDZ*I_ERI_F2xy_S_Dxz_S_C1002001003;
  Double I_ERI_F2xz_S_Dxz_Pz_C1002001003 = I_ERI_F2xz_S_Fx2z_S_C1002001003+CDZ*I_ERI_F2xz_S_Dxz_S_C1002001003;
  Double I_ERI_Fx2y_S_Dxz_Pz_C1002001003 = I_ERI_Fx2y_S_Fx2z_S_C1002001003+CDZ*I_ERI_Fx2y_S_Dxz_S_C1002001003;
  Double I_ERI_Fxyz_S_Dxz_Pz_C1002001003 = I_ERI_Fxyz_S_Fx2z_S_C1002001003+CDZ*I_ERI_Fxyz_S_Dxz_S_C1002001003;
  Double I_ERI_Fx2z_S_Dxz_Pz_C1002001003 = I_ERI_Fx2z_S_Fx2z_S_C1002001003+CDZ*I_ERI_Fx2z_S_Dxz_S_C1002001003;
  Double I_ERI_F3y_S_Dxz_Pz_C1002001003 = I_ERI_F3y_S_Fx2z_S_C1002001003+CDZ*I_ERI_F3y_S_Dxz_S_C1002001003;
  Double I_ERI_F2yz_S_Dxz_Pz_C1002001003 = I_ERI_F2yz_S_Fx2z_S_C1002001003+CDZ*I_ERI_F2yz_S_Dxz_S_C1002001003;
  Double I_ERI_Fy2z_S_Dxz_Pz_C1002001003 = I_ERI_Fy2z_S_Fx2z_S_C1002001003+CDZ*I_ERI_Fy2z_S_Dxz_S_C1002001003;
  Double I_ERI_F3z_S_Dxz_Pz_C1002001003 = I_ERI_F3z_S_Fx2z_S_C1002001003+CDZ*I_ERI_F3z_S_Dxz_S_C1002001003;
  Double I_ERI_F3x_S_D2y_Pz_C1002001003 = I_ERI_F3x_S_F2yz_S_C1002001003+CDZ*I_ERI_F3x_S_D2y_S_C1002001003;
  Double I_ERI_F2xy_S_D2y_Pz_C1002001003 = I_ERI_F2xy_S_F2yz_S_C1002001003+CDZ*I_ERI_F2xy_S_D2y_S_C1002001003;
  Double I_ERI_F2xz_S_D2y_Pz_C1002001003 = I_ERI_F2xz_S_F2yz_S_C1002001003+CDZ*I_ERI_F2xz_S_D2y_S_C1002001003;
  Double I_ERI_Fx2y_S_D2y_Pz_C1002001003 = I_ERI_Fx2y_S_F2yz_S_C1002001003+CDZ*I_ERI_Fx2y_S_D2y_S_C1002001003;
  Double I_ERI_Fxyz_S_D2y_Pz_C1002001003 = I_ERI_Fxyz_S_F2yz_S_C1002001003+CDZ*I_ERI_Fxyz_S_D2y_S_C1002001003;
  Double I_ERI_Fx2z_S_D2y_Pz_C1002001003 = I_ERI_Fx2z_S_F2yz_S_C1002001003+CDZ*I_ERI_Fx2z_S_D2y_S_C1002001003;
  Double I_ERI_F3y_S_D2y_Pz_C1002001003 = I_ERI_F3y_S_F2yz_S_C1002001003+CDZ*I_ERI_F3y_S_D2y_S_C1002001003;
  Double I_ERI_F2yz_S_D2y_Pz_C1002001003 = I_ERI_F2yz_S_F2yz_S_C1002001003+CDZ*I_ERI_F2yz_S_D2y_S_C1002001003;
  Double I_ERI_Fy2z_S_D2y_Pz_C1002001003 = I_ERI_Fy2z_S_F2yz_S_C1002001003+CDZ*I_ERI_Fy2z_S_D2y_S_C1002001003;
  Double I_ERI_F3z_S_D2y_Pz_C1002001003 = I_ERI_F3z_S_F2yz_S_C1002001003+CDZ*I_ERI_F3z_S_D2y_S_C1002001003;
  Double I_ERI_F3x_S_Dyz_Pz_C1002001003 = I_ERI_F3x_S_Fy2z_S_C1002001003+CDZ*I_ERI_F3x_S_Dyz_S_C1002001003;
  Double I_ERI_F2xy_S_Dyz_Pz_C1002001003 = I_ERI_F2xy_S_Fy2z_S_C1002001003+CDZ*I_ERI_F2xy_S_Dyz_S_C1002001003;
  Double I_ERI_F2xz_S_Dyz_Pz_C1002001003 = I_ERI_F2xz_S_Fy2z_S_C1002001003+CDZ*I_ERI_F2xz_S_Dyz_S_C1002001003;
  Double I_ERI_Fx2y_S_Dyz_Pz_C1002001003 = I_ERI_Fx2y_S_Fy2z_S_C1002001003+CDZ*I_ERI_Fx2y_S_Dyz_S_C1002001003;
  Double I_ERI_Fxyz_S_Dyz_Pz_C1002001003 = I_ERI_Fxyz_S_Fy2z_S_C1002001003+CDZ*I_ERI_Fxyz_S_Dyz_S_C1002001003;
  Double I_ERI_Fx2z_S_Dyz_Pz_C1002001003 = I_ERI_Fx2z_S_Fy2z_S_C1002001003+CDZ*I_ERI_Fx2z_S_Dyz_S_C1002001003;
  Double I_ERI_F3y_S_Dyz_Pz_C1002001003 = I_ERI_F3y_S_Fy2z_S_C1002001003+CDZ*I_ERI_F3y_S_Dyz_S_C1002001003;
  Double I_ERI_F2yz_S_Dyz_Pz_C1002001003 = I_ERI_F2yz_S_Fy2z_S_C1002001003+CDZ*I_ERI_F2yz_S_Dyz_S_C1002001003;
  Double I_ERI_Fy2z_S_Dyz_Pz_C1002001003 = I_ERI_Fy2z_S_Fy2z_S_C1002001003+CDZ*I_ERI_Fy2z_S_Dyz_S_C1002001003;
  Double I_ERI_F3z_S_Dyz_Pz_C1002001003 = I_ERI_F3z_S_Fy2z_S_C1002001003+CDZ*I_ERI_F3z_S_Dyz_S_C1002001003;
  Double I_ERI_F3x_S_D2z_Pz_C1002001003 = I_ERI_F3x_S_F3z_S_C1002001003+CDZ*I_ERI_F3x_S_D2z_S_C1002001003;
  Double I_ERI_F2xy_S_D2z_Pz_C1002001003 = I_ERI_F2xy_S_F3z_S_C1002001003+CDZ*I_ERI_F2xy_S_D2z_S_C1002001003;
  Double I_ERI_F2xz_S_D2z_Pz_C1002001003 = I_ERI_F2xz_S_F3z_S_C1002001003+CDZ*I_ERI_F2xz_S_D2z_S_C1002001003;
  Double I_ERI_Fx2y_S_D2z_Pz_C1002001003 = I_ERI_Fx2y_S_F3z_S_C1002001003+CDZ*I_ERI_Fx2y_S_D2z_S_C1002001003;
  Double I_ERI_Fxyz_S_D2z_Pz_C1002001003 = I_ERI_Fxyz_S_F3z_S_C1002001003+CDZ*I_ERI_Fxyz_S_D2z_S_C1002001003;
  Double I_ERI_Fx2z_S_D2z_Pz_C1002001003 = I_ERI_Fx2z_S_F3z_S_C1002001003+CDZ*I_ERI_Fx2z_S_D2z_S_C1002001003;
  Double I_ERI_F3y_S_D2z_Pz_C1002001003 = I_ERI_F3y_S_F3z_S_C1002001003+CDZ*I_ERI_F3y_S_D2z_S_C1002001003;
  Double I_ERI_F2yz_S_D2z_Pz_C1002001003 = I_ERI_F2yz_S_F3z_S_C1002001003+CDZ*I_ERI_F2yz_S_D2z_S_C1002001003;
  Double I_ERI_Fy2z_S_D2z_Pz_C1002001003 = I_ERI_Fy2z_S_F3z_S_C1002001003+CDZ*I_ERI_Fy2z_S_D2z_S_C1002001003;
  Double I_ERI_F3z_S_D2z_Pz_C1002001003 = I_ERI_F3z_S_F3z_S_C1002001003+CDZ*I_ERI_F3z_S_D2z_S_C1002001003;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_D_P_C1002001003
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_F_S_C1002001003
   * RHS shell quartet name: SQ_ERI_G_S_D_S_C1002001003
   ************************************************************/
  Double I_ERI_G4x_S_D2x_Px_C1002001003 = I_ERI_G4x_S_F3x_S_C1002001003+CDX*I_ERI_G4x_S_D2x_S_C1002001003;
  Double I_ERI_G3xy_S_D2x_Px_C1002001003 = I_ERI_G3xy_S_F3x_S_C1002001003+CDX*I_ERI_G3xy_S_D2x_S_C1002001003;
  Double I_ERI_G3xz_S_D2x_Px_C1002001003 = I_ERI_G3xz_S_F3x_S_C1002001003+CDX*I_ERI_G3xz_S_D2x_S_C1002001003;
  Double I_ERI_G2x2y_S_D2x_Px_C1002001003 = I_ERI_G2x2y_S_F3x_S_C1002001003+CDX*I_ERI_G2x2y_S_D2x_S_C1002001003;
  Double I_ERI_G2xyz_S_D2x_Px_C1002001003 = I_ERI_G2xyz_S_F3x_S_C1002001003+CDX*I_ERI_G2xyz_S_D2x_S_C1002001003;
  Double I_ERI_G2x2z_S_D2x_Px_C1002001003 = I_ERI_G2x2z_S_F3x_S_C1002001003+CDX*I_ERI_G2x2z_S_D2x_S_C1002001003;
  Double I_ERI_Gx3y_S_D2x_Px_C1002001003 = I_ERI_Gx3y_S_F3x_S_C1002001003+CDX*I_ERI_Gx3y_S_D2x_S_C1002001003;
  Double I_ERI_Gx2yz_S_D2x_Px_C1002001003 = I_ERI_Gx2yz_S_F3x_S_C1002001003+CDX*I_ERI_Gx2yz_S_D2x_S_C1002001003;
  Double I_ERI_Gxy2z_S_D2x_Px_C1002001003 = I_ERI_Gxy2z_S_F3x_S_C1002001003+CDX*I_ERI_Gxy2z_S_D2x_S_C1002001003;
  Double I_ERI_Gx3z_S_D2x_Px_C1002001003 = I_ERI_Gx3z_S_F3x_S_C1002001003+CDX*I_ERI_Gx3z_S_D2x_S_C1002001003;
  Double I_ERI_G4y_S_D2x_Px_C1002001003 = I_ERI_G4y_S_F3x_S_C1002001003+CDX*I_ERI_G4y_S_D2x_S_C1002001003;
  Double I_ERI_G3yz_S_D2x_Px_C1002001003 = I_ERI_G3yz_S_F3x_S_C1002001003+CDX*I_ERI_G3yz_S_D2x_S_C1002001003;
  Double I_ERI_G2y2z_S_D2x_Px_C1002001003 = I_ERI_G2y2z_S_F3x_S_C1002001003+CDX*I_ERI_G2y2z_S_D2x_S_C1002001003;
  Double I_ERI_Gy3z_S_D2x_Px_C1002001003 = I_ERI_Gy3z_S_F3x_S_C1002001003+CDX*I_ERI_Gy3z_S_D2x_S_C1002001003;
  Double I_ERI_G4z_S_D2x_Px_C1002001003 = I_ERI_G4z_S_F3x_S_C1002001003+CDX*I_ERI_G4z_S_D2x_S_C1002001003;
  Double I_ERI_G4x_S_Dxy_Px_C1002001003 = I_ERI_G4x_S_F2xy_S_C1002001003+CDX*I_ERI_G4x_S_Dxy_S_C1002001003;
  Double I_ERI_G3xy_S_Dxy_Px_C1002001003 = I_ERI_G3xy_S_F2xy_S_C1002001003+CDX*I_ERI_G3xy_S_Dxy_S_C1002001003;
  Double I_ERI_G3xz_S_Dxy_Px_C1002001003 = I_ERI_G3xz_S_F2xy_S_C1002001003+CDX*I_ERI_G3xz_S_Dxy_S_C1002001003;
  Double I_ERI_G2x2y_S_Dxy_Px_C1002001003 = I_ERI_G2x2y_S_F2xy_S_C1002001003+CDX*I_ERI_G2x2y_S_Dxy_S_C1002001003;
  Double I_ERI_G2xyz_S_Dxy_Px_C1002001003 = I_ERI_G2xyz_S_F2xy_S_C1002001003+CDX*I_ERI_G2xyz_S_Dxy_S_C1002001003;
  Double I_ERI_G2x2z_S_Dxy_Px_C1002001003 = I_ERI_G2x2z_S_F2xy_S_C1002001003+CDX*I_ERI_G2x2z_S_Dxy_S_C1002001003;
  Double I_ERI_Gx3y_S_Dxy_Px_C1002001003 = I_ERI_Gx3y_S_F2xy_S_C1002001003+CDX*I_ERI_Gx3y_S_Dxy_S_C1002001003;
  Double I_ERI_Gx2yz_S_Dxy_Px_C1002001003 = I_ERI_Gx2yz_S_F2xy_S_C1002001003+CDX*I_ERI_Gx2yz_S_Dxy_S_C1002001003;
  Double I_ERI_Gxy2z_S_Dxy_Px_C1002001003 = I_ERI_Gxy2z_S_F2xy_S_C1002001003+CDX*I_ERI_Gxy2z_S_Dxy_S_C1002001003;
  Double I_ERI_Gx3z_S_Dxy_Px_C1002001003 = I_ERI_Gx3z_S_F2xy_S_C1002001003+CDX*I_ERI_Gx3z_S_Dxy_S_C1002001003;
  Double I_ERI_G4y_S_Dxy_Px_C1002001003 = I_ERI_G4y_S_F2xy_S_C1002001003+CDX*I_ERI_G4y_S_Dxy_S_C1002001003;
  Double I_ERI_G3yz_S_Dxy_Px_C1002001003 = I_ERI_G3yz_S_F2xy_S_C1002001003+CDX*I_ERI_G3yz_S_Dxy_S_C1002001003;
  Double I_ERI_G2y2z_S_Dxy_Px_C1002001003 = I_ERI_G2y2z_S_F2xy_S_C1002001003+CDX*I_ERI_G2y2z_S_Dxy_S_C1002001003;
  Double I_ERI_Gy3z_S_Dxy_Px_C1002001003 = I_ERI_Gy3z_S_F2xy_S_C1002001003+CDX*I_ERI_Gy3z_S_Dxy_S_C1002001003;
  Double I_ERI_G4z_S_Dxy_Px_C1002001003 = I_ERI_G4z_S_F2xy_S_C1002001003+CDX*I_ERI_G4z_S_Dxy_S_C1002001003;
  Double I_ERI_G4x_S_Dxz_Px_C1002001003 = I_ERI_G4x_S_F2xz_S_C1002001003+CDX*I_ERI_G4x_S_Dxz_S_C1002001003;
  Double I_ERI_G3xy_S_Dxz_Px_C1002001003 = I_ERI_G3xy_S_F2xz_S_C1002001003+CDX*I_ERI_G3xy_S_Dxz_S_C1002001003;
  Double I_ERI_G3xz_S_Dxz_Px_C1002001003 = I_ERI_G3xz_S_F2xz_S_C1002001003+CDX*I_ERI_G3xz_S_Dxz_S_C1002001003;
  Double I_ERI_G2x2y_S_Dxz_Px_C1002001003 = I_ERI_G2x2y_S_F2xz_S_C1002001003+CDX*I_ERI_G2x2y_S_Dxz_S_C1002001003;
  Double I_ERI_G2xyz_S_Dxz_Px_C1002001003 = I_ERI_G2xyz_S_F2xz_S_C1002001003+CDX*I_ERI_G2xyz_S_Dxz_S_C1002001003;
  Double I_ERI_G2x2z_S_Dxz_Px_C1002001003 = I_ERI_G2x2z_S_F2xz_S_C1002001003+CDX*I_ERI_G2x2z_S_Dxz_S_C1002001003;
  Double I_ERI_Gx3y_S_Dxz_Px_C1002001003 = I_ERI_Gx3y_S_F2xz_S_C1002001003+CDX*I_ERI_Gx3y_S_Dxz_S_C1002001003;
  Double I_ERI_Gx2yz_S_Dxz_Px_C1002001003 = I_ERI_Gx2yz_S_F2xz_S_C1002001003+CDX*I_ERI_Gx2yz_S_Dxz_S_C1002001003;
  Double I_ERI_Gxy2z_S_Dxz_Px_C1002001003 = I_ERI_Gxy2z_S_F2xz_S_C1002001003+CDX*I_ERI_Gxy2z_S_Dxz_S_C1002001003;
  Double I_ERI_Gx3z_S_Dxz_Px_C1002001003 = I_ERI_Gx3z_S_F2xz_S_C1002001003+CDX*I_ERI_Gx3z_S_Dxz_S_C1002001003;
  Double I_ERI_G4y_S_Dxz_Px_C1002001003 = I_ERI_G4y_S_F2xz_S_C1002001003+CDX*I_ERI_G4y_S_Dxz_S_C1002001003;
  Double I_ERI_G3yz_S_Dxz_Px_C1002001003 = I_ERI_G3yz_S_F2xz_S_C1002001003+CDX*I_ERI_G3yz_S_Dxz_S_C1002001003;
  Double I_ERI_G2y2z_S_Dxz_Px_C1002001003 = I_ERI_G2y2z_S_F2xz_S_C1002001003+CDX*I_ERI_G2y2z_S_Dxz_S_C1002001003;
  Double I_ERI_Gy3z_S_Dxz_Px_C1002001003 = I_ERI_Gy3z_S_F2xz_S_C1002001003+CDX*I_ERI_Gy3z_S_Dxz_S_C1002001003;
  Double I_ERI_G4z_S_Dxz_Px_C1002001003 = I_ERI_G4z_S_F2xz_S_C1002001003+CDX*I_ERI_G4z_S_Dxz_S_C1002001003;
  Double I_ERI_G4x_S_D2y_Px_C1002001003 = I_ERI_G4x_S_Fx2y_S_C1002001003+CDX*I_ERI_G4x_S_D2y_S_C1002001003;
  Double I_ERI_G3xy_S_D2y_Px_C1002001003 = I_ERI_G3xy_S_Fx2y_S_C1002001003+CDX*I_ERI_G3xy_S_D2y_S_C1002001003;
  Double I_ERI_G3xz_S_D2y_Px_C1002001003 = I_ERI_G3xz_S_Fx2y_S_C1002001003+CDX*I_ERI_G3xz_S_D2y_S_C1002001003;
  Double I_ERI_G2x2y_S_D2y_Px_C1002001003 = I_ERI_G2x2y_S_Fx2y_S_C1002001003+CDX*I_ERI_G2x2y_S_D2y_S_C1002001003;
  Double I_ERI_G2xyz_S_D2y_Px_C1002001003 = I_ERI_G2xyz_S_Fx2y_S_C1002001003+CDX*I_ERI_G2xyz_S_D2y_S_C1002001003;
  Double I_ERI_G2x2z_S_D2y_Px_C1002001003 = I_ERI_G2x2z_S_Fx2y_S_C1002001003+CDX*I_ERI_G2x2z_S_D2y_S_C1002001003;
  Double I_ERI_Gx3y_S_D2y_Px_C1002001003 = I_ERI_Gx3y_S_Fx2y_S_C1002001003+CDX*I_ERI_Gx3y_S_D2y_S_C1002001003;
  Double I_ERI_Gx2yz_S_D2y_Px_C1002001003 = I_ERI_Gx2yz_S_Fx2y_S_C1002001003+CDX*I_ERI_Gx2yz_S_D2y_S_C1002001003;
  Double I_ERI_Gxy2z_S_D2y_Px_C1002001003 = I_ERI_Gxy2z_S_Fx2y_S_C1002001003+CDX*I_ERI_Gxy2z_S_D2y_S_C1002001003;
  Double I_ERI_Gx3z_S_D2y_Px_C1002001003 = I_ERI_Gx3z_S_Fx2y_S_C1002001003+CDX*I_ERI_Gx3z_S_D2y_S_C1002001003;
  Double I_ERI_G4y_S_D2y_Px_C1002001003 = I_ERI_G4y_S_Fx2y_S_C1002001003+CDX*I_ERI_G4y_S_D2y_S_C1002001003;
  Double I_ERI_G3yz_S_D2y_Px_C1002001003 = I_ERI_G3yz_S_Fx2y_S_C1002001003+CDX*I_ERI_G3yz_S_D2y_S_C1002001003;
  Double I_ERI_G2y2z_S_D2y_Px_C1002001003 = I_ERI_G2y2z_S_Fx2y_S_C1002001003+CDX*I_ERI_G2y2z_S_D2y_S_C1002001003;
  Double I_ERI_Gy3z_S_D2y_Px_C1002001003 = I_ERI_Gy3z_S_Fx2y_S_C1002001003+CDX*I_ERI_Gy3z_S_D2y_S_C1002001003;
  Double I_ERI_G4z_S_D2y_Px_C1002001003 = I_ERI_G4z_S_Fx2y_S_C1002001003+CDX*I_ERI_G4z_S_D2y_S_C1002001003;
  Double I_ERI_G4x_S_Dyz_Px_C1002001003 = I_ERI_G4x_S_Fxyz_S_C1002001003+CDX*I_ERI_G4x_S_Dyz_S_C1002001003;
  Double I_ERI_G3xy_S_Dyz_Px_C1002001003 = I_ERI_G3xy_S_Fxyz_S_C1002001003+CDX*I_ERI_G3xy_S_Dyz_S_C1002001003;
  Double I_ERI_G3xz_S_Dyz_Px_C1002001003 = I_ERI_G3xz_S_Fxyz_S_C1002001003+CDX*I_ERI_G3xz_S_Dyz_S_C1002001003;
  Double I_ERI_G2x2y_S_Dyz_Px_C1002001003 = I_ERI_G2x2y_S_Fxyz_S_C1002001003+CDX*I_ERI_G2x2y_S_Dyz_S_C1002001003;
  Double I_ERI_G2xyz_S_Dyz_Px_C1002001003 = I_ERI_G2xyz_S_Fxyz_S_C1002001003+CDX*I_ERI_G2xyz_S_Dyz_S_C1002001003;
  Double I_ERI_G2x2z_S_Dyz_Px_C1002001003 = I_ERI_G2x2z_S_Fxyz_S_C1002001003+CDX*I_ERI_G2x2z_S_Dyz_S_C1002001003;
  Double I_ERI_Gx3y_S_Dyz_Px_C1002001003 = I_ERI_Gx3y_S_Fxyz_S_C1002001003+CDX*I_ERI_Gx3y_S_Dyz_S_C1002001003;
  Double I_ERI_Gx2yz_S_Dyz_Px_C1002001003 = I_ERI_Gx2yz_S_Fxyz_S_C1002001003+CDX*I_ERI_Gx2yz_S_Dyz_S_C1002001003;
  Double I_ERI_Gxy2z_S_Dyz_Px_C1002001003 = I_ERI_Gxy2z_S_Fxyz_S_C1002001003+CDX*I_ERI_Gxy2z_S_Dyz_S_C1002001003;
  Double I_ERI_Gx3z_S_Dyz_Px_C1002001003 = I_ERI_Gx3z_S_Fxyz_S_C1002001003+CDX*I_ERI_Gx3z_S_Dyz_S_C1002001003;
  Double I_ERI_G4y_S_Dyz_Px_C1002001003 = I_ERI_G4y_S_Fxyz_S_C1002001003+CDX*I_ERI_G4y_S_Dyz_S_C1002001003;
  Double I_ERI_G3yz_S_Dyz_Px_C1002001003 = I_ERI_G3yz_S_Fxyz_S_C1002001003+CDX*I_ERI_G3yz_S_Dyz_S_C1002001003;
  Double I_ERI_G2y2z_S_Dyz_Px_C1002001003 = I_ERI_G2y2z_S_Fxyz_S_C1002001003+CDX*I_ERI_G2y2z_S_Dyz_S_C1002001003;
  Double I_ERI_Gy3z_S_Dyz_Px_C1002001003 = I_ERI_Gy3z_S_Fxyz_S_C1002001003+CDX*I_ERI_Gy3z_S_Dyz_S_C1002001003;
  Double I_ERI_G4z_S_Dyz_Px_C1002001003 = I_ERI_G4z_S_Fxyz_S_C1002001003+CDX*I_ERI_G4z_S_Dyz_S_C1002001003;
  Double I_ERI_G4x_S_D2z_Px_C1002001003 = I_ERI_G4x_S_Fx2z_S_C1002001003+CDX*I_ERI_G4x_S_D2z_S_C1002001003;
  Double I_ERI_G3xy_S_D2z_Px_C1002001003 = I_ERI_G3xy_S_Fx2z_S_C1002001003+CDX*I_ERI_G3xy_S_D2z_S_C1002001003;
  Double I_ERI_G3xz_S_D2z_Px_C1002001003 = I_ERI_G3xz_S_Fx2z_S_C1002001003+CDX*I_ERI_G3xz_S_D2z_S_C1002001003;
  Double I_ERI_G2x2y_S_D2z_Px_C1002001003 = I_ERI_G2x2y_S_Fx2z_S_C1002001003+CDX*I_ERI_G2x2y_S_D2z_S_C1002001003;
  Double I_ERI_G2xyz_S_D2z_Px_C1002001003 = I_ERI_G2xyz_S_Fx2z_S_C1002001003+CDX*I_ERI_G2xyz_S_D2z_S_C1002001003;
  Double I_ERI_G2x2z_S_D2z_Px_C1002001003 = I_ERI_G2x2z_S_Fx2z_S_C1002001003+CDX*I_ERI_G2x2z_S_D2z_S_C1002001003;
  Double I_ERI_Gx3y_S_D2z_Px_C1002001003 = I_ERI_Gx3y_S_Fx2z_S_C1002001003+CDX*I_ERI_Gx3y_S_D2z_S_C1002001003;
  Double I_ERI_Gx2yz_S_D2z_Px_C1002001003 = I_ERI_Gx2yz_S_Fx2z_S_C1002001003+CDX*I_ERI_Gx2yz_S_D2z_S_C1002001003;
  Double I_ERI_Gxy2z_S_D2z_Px_C1002001003 = I_ERI_Gxy2z_S_Fx2z_S_C1002001003+CDX*I_ERI_Gxy2z_S_D2z_S_C1002001003;
  Double I_ERI_Gx3z_S_D2z_Px_C1002001003 = I_ERI_Gx3z_S_Fx2z_S_C1002001003+CDX*I_ERI_Gx3z_S_D2z_S_C1002001003;
  Double I_ERI_G4y_S_D2z_Px_C1002001003 = I_ERI_G4y_S_Fx2z_S_C1002001003+CDX*I_ERI_G4y_S_D2z_S_C1002001003;
  Double I_ERI_G3yz_S_D2z_Px_C1002001003 = I_ERI_G3yz_S_Fx2z_S_C1002001003+CDX*I_ERI_G3yz_S_D2z_S_C1002001003;
  Double I_ERI_G2y2z_S_D2z_Px_C1002001003 = I_ERI_G2y2z_S_Fx2z_S_C1002001003+CDX*I_ERI_G2y2z_S_D2z_S_C1002001003;
  Double I_ERI_Gy3z_S_D2z_Px_C1002001003 = I_ERI_Gy3z_S_Fx2z_S_C1002001003+CDX*I_ERI_Gy3z_S_D2z_S_C1002001003;
  Double I_ERI_G4z_S_D2z_Px_C1002001003 = I_ERI_G4z_S_Fx2z_S_C1002001003+CDX*I_ERI_G4z_S_D2z_S_C1002001003;
  Double I_ERI_G4x_S_D2x_Py_C1002001003 = I_ERI_G4x_S_F2xy_S_C1002001003+CDY*I_ERI_G4x_S_D2x_S_C1002001003;
  Double I_ERI_G3xy_S_D2x_Py_C1002001003 = I_ERI_G3xy_S_F2xy_S_C1002001003+CDY*I_ERI_G3xy_S_D2x_S_C1002001003;
  Double I_ERI_G3xz_S_D2x_Py_C1002001003 = I_ERI_G3xz_S_F2xy_S_C1002001003+CDY*I_ERI_G3xz_S_D2x_S_C1002001003;
  Double I_ERI_G2x2y_S_D2x_Py_C1002001003 = I_ERI_G2x2y_S_F2xy_S_C1002001003+CDY*I_ERI_G2x2y_S_D2x_S_C1002001003;
  Double I_ERI_G2xyz_S_D2x_Py_C1002001003 = I_ERI_G2xyz_S_F2xy_S_C1002001003+CDY*I_ERI_G2xyz_S_D2x_S_C1002001003;
  Double I_ERI_G2x2z_S_D2x_Py_C1002001003 = I_ERI_G2x2z_S_F2xy_S_C1002001003+CDY*I_ERI_G2x2z_S_D2x_S_C1002001003;
  Double I_ERI_Gx3y_S_D2x_Py_C1002001003 = I_ERI_Gx3y_S_F2xy_S_C1002001003+CDY*I_ERI_Gx3y_S_D2x_S_C1002001003;
  Double I_ERI_Gx2yz_S_D2x_Py_C1002001003 = I_ERI_Gx2yz_S_F2xy_S_C1002001003+CDY*I_ERI_Gx2yz_S_D2x_S_C1002001003;
  Double I_ERI_Gxy2z_S_D2x_Py_C1002001003 = I_ERI_Gxy2z_S_F2xy_S_C1002001003+CDY*I_ERI_Gxy2z_S_D2x_S_C1002001003;
  Double I_ERI_Gx3z_S_D2x_Py_C1002001003 = I_ERI_Gx3z_S_F2xy_S_C1002001003+CDY*I_ERI_Gx3z_S_D2x_S_C1002001003;
  Double I_ERI_G4y_S_D2x_Py_C1002001003 = I_ERI_G4y_S_F2xy_S_C1002001003+CDY*I_ERI_G4y_S_D2x_S_C1002001003;
  Double I_ERI_G3yz_S_D2x_Py_C1002001003 = I_ERI_G3yz_S_F2xy_S_C1002001003+CDY*I_ERI_G3yz_S_D2x_S_C1002001003;
  Double I_ERI_G2y2z_S_D2x_Py_C1002001003 = I_ERI_G2y2z_S_F2xy_S_C1002001003+CDY*I_ERI_G2y2z_S_D2x_S_C1002001003;
  Double I_ERI_Gy3z_S_D2x_Py_C1002001003 = I_ERI_Gy3z_S_F2xy_S_C1002001003+CDY*I_ERI_Gy3z_S_D2x_S_C1002001003;
  Double I_ERI_G4z_S_D2x_Py_C1002001003 = I_ERI_G4z_S_F2xy_S_C1002001003+CDY*I_ERI_G4z_S_D2x_S_C1002001003;
  Double I_ERI_G4x_S_Dxy_Py_C1002001003 = I_ERI_G4x_S_Fx2y_S_C1002001003+CDY*I_ERI_G4x_S_Dxy_S_C1002001003;
  Double I_ERI_G3xy_S_Dxy_Py_C1002001003 = I_ERI_G3xy_S_Fx2y_S_C1002001003+CDY*I_ERI_G3xy_S_Dxy_S_C1002001003;
  Double I_ERI_G3xz_S_Dxy_Py_C1002001003 = I_ERI_G3xz_S_Fx2y_S_C1002001003+CDY*I_ERI_G3xz_S_Dxy_S_C1002001003;
  Double I_ERI_G2x2y_S_Dxy_Py_C1002001003 = I_ERI_G2x2y_S_Fx2y_S_C1002001003+CDY*I_ERI_G2x2y_S_Dxy_S_C1002001003;
  Double I_ERI_G2xyz_S_Dxy_Py_C1002001003 = I_ERI_G2xyz_S_Fx2y_S_C1002001003+CDY*I_ERI_G2xyz_S_Dxy_S_C1002001003;
  Double I_ERI_G2x2z_S_Dxy_Py_C1002001003 = I_ERI_G2x2z_S_Fx2y_S_C1002001003+CDY*I_ERI_G2x2z_S_Dxy_S_C1002001003;
  Double I_ERI_Gx3y_S_Dxy_Py_C1002001003 = I_ERI_Gx3y_S_Fx2y_S_C1002001003+CDY*I_ERI_Gx3y_S_Dxy_S_C1002001003;
  Double I_ERI_Gx2yz_S_Dxy_Py_C1002001003 = I_ERI_Gx2yz_S_Fx2y_S_C1002001003+CDY*I_ERI_Gx2yz_S_Dxy_S_C1002001003;
  Double I_ERI_Gxy2z_S_Dxy_Py_C1002001003 = I_ERI_Gxy2z_S_Fx2y_S_C1002001003+CDY*I_ERI_Gxy2z_S_Dxy_S_C1002001003;
  Double I_ERI_Gx3z_S_Dxy_Py_C1002001003 = I_ERI_Gx3z_S_Fx2y_S_C1002001003+CDY*I_ERI_Gx3z_S_Dxy_S_C1002001003;
  Double I_ERI_G4y_S_Dxy_Py_C1002001003 = I_ERI_G4y_S_Fx2y_S_C1002001003+CDY*I_ERI_G4y_S_Dxy_S_C1002001003;
  Double I_ERI_G3yz_S_Dxy_Py_C1002001003 = I_ERI_G3yz_S_Fx2y_S_C1002001003+CDY*I_ERI_G3yz_S_Dxy_S_C1002001003;
  Double I_ERI_G2y2z_S_Dxy_Py_C1002001003 = I_ERI_G2y2z_S_Fx2y_S_C1002001003+CDY*I_ERI_G2y2z_S_Dxy_S_C1002001003;
  Double I_ERI_Gy3z_S_Dxy_Py_C1002001003 = I_ERI_Gy3z_S_Fx2y_S_C1002001003+CDY*I_ERI_Gy3z_S_Dxy_S_C1002001003;
  Double I_ERI_G4z_S_Dxy_Py_C1002001003 = I_ERI_G4z_S_Fx2y_S_C1002001003+CDY*I_ERI_G4z_S_Dxy_S_C1002001003;
  Double I_ERI_G4x_S_Dxz_Py_C1002001003 = I_ERI_G4x_S_Fxyz_S_C1002001003+CDY*I_ERI_G4x_S_Dxz_S_C1002001003;
  Double I_ERI_G3xy_S_Dxz_Py_C1002001003 = I_ERI_G3xy_S_Fxyz_S_C1002001003+CDY*I_ERI_G3xy_S_Dxz_S_C1002001003;
  Double I_ERI_G3xz_S_Dxz_Py_C1002001003 = I_ERI_G3xz_S_Fxyz_S_C1002001003+CDY*I_ERI_G3xz_S_Dxz_S_C1002001003;
  Double I_ERI_G2x2y_S_Dxz_Py_C1002001003 = I_ERI_G2x2y_S_Fxyz_S_C1002001003+CDY*I_ERI_G2x2y_S_Dxz_S_C1002001003;
  Double I_ERI_G2xyz_S_Dxz_Py_C1002001003 = I_ERI_G2xyz_S_Fxyz_S_C1002001003+CDY*I_ERI_G2xyz_S_Dxz_S_C1002001003;
  Double I_ERI_G2x2z_S_Dxz_Py_C1002001003 = I_ERI_G2x2z_S_Fxyz_S_C1002001003+CDY*I_ERI_G2x2z_S_Dxz_S_C1002001003;
  Double I_ERI_Gx3y_S_Dxz_Py_C1002001003 = I_ERI_Gx3y_S_Fxyz_S_C1002001003+CDY*I_ERI_Gx3y_S_Dxz_S_C1002001003;
  Double I_ERI_Gx2yz_S_Dxz_Py_C1002001003 = I_ERI_Gx2yz_S_Fxyz_S_C1002001003+CDY*I_ERI_Gx2yz_S_Dxz_S_C1002001003;
  Double I_ERI_Gxy2z_S_Dxz_Py_C1002001003 = I_ERI_Gxy2z_S_Fxyz_S_C1002001003+CDY*I_ERI_Gxy2z_S_Dxz_S_C1002001003;
  Double I_ERI_Gx3z_S_Dxz_Py_C1002001003 = I_ERI_Gx3z_S_Fxyz_S_C1002001003+CDY*I_ERI_Gx3z_S_Dxz_S_C1002001003;
  Double I_ERI_G4y_S_Dxz_Py_C1002001003 = I_ERI_G4y_S_Fxyz_S_C1002001003+CDY*I_ERI_G4y_S_Dxz_S_C1002001003;
  Double I_ERI_G3yz_S_Dxz_Py_C1002001003 = I_ERI_G3yz_S_Fxyz_S_C1002001003+CDY*I_ERI_G3yz_S_Dxz_S_C1002001003;
  Double I_ERI_G2y2z_S_Dxz_Py_C1002001003 = I_ERI_G2y2z_S_Fxyz_S_C1002001003+CDY*I_ERI_G2y2z_S_Dxz_S_C1002001003;
  Double I_ERI_Gy3z_S_Dxz_Py_C1002001003 = I_ERI_Gy3z_S_Fxyz_S_C1002001003+CDY*I_ERI_Gy3z_S_Dxz_S_C1002001003;
  Double I_ERI_G4z_S_Dxz_Py_C1002001003 = I_ERI_G4z_S_Fxyz_S_C1002001003+CDY*I_ERI_G4z_S_Dxz_S_C1002001003;
  Double I_ERI_G4x_S_D2y_Py_C1002001003 = I_ERI_G4x_S_F3y_S_C1002001003+CDY*I_ERI_G4x_S_D2y_S_C1002001003;
  Double I_ERI_G3xy_S_D2y_Py_C1002001003 = I_ERI_G3xy_S_F3y_S_C1002001003+CDY*I_ERI_G3xy_S_D2y_S_C1002001003;
  Double I_ERI_G3xz_S_D2y_Py_C1002001003 = I_ERI_G3xz_S_F3y_S_C1002001003+CDY*I_ERI_G3xz_S_D2y_S_C1002001003;
  Double I_ERI_G2x2y_S_D2y_Py_C1002001003 = I_ERI_G2x2y_S_F3y_S_C1002001003+CDY*I_ERI_G2x2y_S_D2y_S_C1002001003;
  Double I_ERI_G2xyz_S_D2y_Py_C1002001003 = I_ERI_G2xyz_S_F3y_S_C1002001003+CDY*I_ERI_G2xyz_S_D2y_S_C1002001003;
  Double I_ERI_G2x2z_S_D2y_Py_C1002001003 = I_ERI_G2x2z_S_F3y_S_C1002001003+CDY*I_ERI_G2x2z_S_D2y_S_C1002001003;
  Double I_ERI_Gx3y_S_D2y_Py_C1002001003 = I_ERI_Gx3y_S_F3y_S_C1002001003+CDY*I_ERI_Gx3y_S_D2y_S_C1002001003;
  Double I_ERI_Gx2yz_S_D2y_Py_C1002001003 = I_ERI_Gx2yz_S_F3y_S_C1002001003+CDY*I_ERI_Gx2yz_S_D2y_S_C1002001003;
  Double I_ERI_Gxy2z_S_D2y_Py_C1002001003 = I_ERI_Gxy2z_S_F3y_S_C1002001003+CDY*I_ERI_Gxy2z_S_D2y_S_C1002001003;
  Double I_ERI_Gx3z_S_D2y_Py_C1002001003 = I_ERI_Gx3z_S_F3y_S_C1002001003+CDY*I_ERI_Gx3z_S_D2y_S_C1002001003;
  Double I_ERI_G4y_S_D2y_Py_C1002001003 = I_ERI_G4y_S_F3y_S_C1002001003+CDY*I_ERI_G4y_S_D2y_S_C1002001003;
  Double I_ERI_G3yz_S_D2y_Py_C1002001003 = I_ERI_G3yz_S_F3y_S_C1002001003+CDY*I_ERI_G3yz_S_D2y_S_C1002001003;
  Double I_ERI_G2y2z_S_D2y_Py_C1002001003 = I_ERI_G2y2z_S_F3y_S_C1002001003+CDY*I_ERI_G2y2z_S_D2y_S_C1002001003;
  Double I_ERI_Gy3z_S_D2y_Py_C1002001003 = I_ERI_Gy3z_S_F3y_S_C1002001003+CDY*I_ERI_Gy3z_S_D2y_S_C1002001003;
  Double I_ERI_G4z_S_D2y_Py_C1002001003 = I_ERI_G4z_S_F3y_S_C1002001003+CDY*I_ERI_G4z_S_D2y_S_C1002001003;
  Double I_ERI_G4x_S_Dyz_Py_C1002001003 = I_ERI_G4x_S_F2yz_S_C1002001003+CDY*I_ERI_G4x_S_Dyz_S_C1002001003;
  Double I_ERI_G3xy_S_Dyz_Py_C1002001003 = I_ERI_G3xy_S_F2yz_S_C1002001003+CDY*I_ERI_G3xy_S_Dyz_S_C1002001003;
  Double I_ERI_G3xz_S_Dyz_Py_C1002001003 = I_ERI_G3xz_S_F2yz_S_C1002001003+CDY*I_ERI_G3xz_S_Dyz_S_C1002001003;
  Double I_ERI_G2x2y_S_Dyz_Py_C1002001003 = I_ERI_G2x2y_S_F2yz_S_C1002001003+CDY*I_ERI_G2x2y_S_Dyz_S_C1002001003;
  Double I_ERI_G2xyz_S_Dyz_Py_C1002001003 = I_ERI_G2xyz_S_F2yz_S_C1002001003+CDY*I_ERI_G2xyz_S_Dyz_S_C1002001003;
  Double I_ERI_G2x2z_S_Dyz_Py_C1002001003 = I_ERI_G2x2z_S_F2yz_S_C1002001003+CDY*I_ERI_G2x2z_S_Dyz_S_C1002001003;
  Double I_ERI_Gx3y_S_Dyz_Py_C1002001003 = I_ERI_Gx3y_S_F2yz_S_C1002001003+CDY*I_ERI_Gx3y_S_Dyz_S_C1002001003;
  Double I_ERI_Gx2yz_S_Dyz_Py_C1002001003 = I_ERI_Gx2yz_S_F2yz_S_C1002001003+CDY*I_ERI_Gx2yz_S_Dyz_S_C1002001003;
  Double I_ERI_Gxy2z_S_Dyz_Py_C1002001003 = I_ERI_Gxy2z_S_F2yz_S_C1002001003+CDY*I_ERI_Gxy2z_S_Dyz_S_C1002001003;
  Double I_ERI_Gx3z_S_Dyz_Py_C1002001003 = I_ERI_Gx3z_S_F2yz_S_C1002001003+CDY*I_ERI_Gx3z_S_Dyz_S_C1002001003;
  Double I_ERI_G4y_S_Dyz_Py_C1002001003 = I_ERI_G4y_S_F2yz_S_C1002001003+CDY*I_ERI_G4y_S_Dyz_S_C1002001003;
  Double I_ERI_G3yz_S_Dyz_Py_C1002001003 = I_ERI_G3yz_S_F2yz_S_C1002001003+CDY*I_ERI_G3yz_S_Dyz_S_C1002001003;
  Double I_ERI_G2y2z_S_Dyz_Py_C1002001003 = I_ERI_G2y2z_S_F2yz_S_C1002001003+CDY*I_ERI_G2y2z_S_Dyz_S_C1002001003;
  Double I_ERI_Gy3z_S_Dyz_Py_C1002001003 = I_ERI_Gy3z_S_F2yz_S_C1002001003+CDY*I_ERI_Gy3z_S_Dyz_S_C1002001003;
  Double I_ERI_G4z_S_Dyz_Py_C1002001003 = I_ERI_G4z_S_F2yz_S_C1002001003+CDY*I_ERI_G4z_S_Dyz_S_C1002001003;
  Double I_ERI_G4x_S_D2z_Py_C1002001003 = I_ERI_G4x_S_Fy2z_S_C1002001003+CDY*I_ERI_G4x_S_D2z_S_C1002001003;
  Double I_ERI_G3xy_S_D2z_Py_C1002001003 = I_ERI_G3xy_S_Fy2z_S_C1002001003+CDY*I_ERI_G3xy_S_D2z_S_C1002001003;
  Double I_ERI_G3xz_S_D2z_Py_C1002001003 = I_ERI_G3xz_S_Fy2z_S_C1002001003+CDY*I_ERI_G3xz_S_D2z_S_C1002001003;
  Double I_ERI_G2x2y_S_D2z_Py_C1002001003 = I_ERI_G2x2y_S_Fy2z_S_C1002001003+CDY*I_ERI_G2x2y_S_D2z_S_C1002001003;
  Double I_ERI_G2xyz_S_D2z_Py_C1002001003 = I_ERI_G2xyz_S_Fy2z_S_C1002001003+CDY*I_ERI_G2xyz_S_D2z_S_C1002001003;
  Double I_ERI_G2x2z_S_D2z_Py_C1002001003 = I_ERI_G2x2z_S_Fy2z_S_C1002001003+CDY*I_ERI_G2x2z_S_D2z_S_C1002001003;
  Double I_ERI_Gx3y_S_D2z_Py_C1002001003 = I_ERI_Gx3y_S_Fy2z_S_C1002001003+CDY*I_ERI_Gx3y_S_D2z_S_C1002001003;
  Double I_ERI_Gx2yz_S_D2z_Py_C1002001003 = I_ERI_Gx2yz_S_Fy2z_S_C1002001003+CDY*I_ERI_Gx2yz_S_D2z_S_C1002001003;
  Double I_ERI_Gxy2z_S_D2z_Py_C1002001003 = I_ERI_Gxy2z_S_Fy2z_S_C1002001003+CDY*I_ERI_Gxy2z_S_D2z_S_C1002001003;
  Double I_ERI_Gx3z_S_D2z_Py_C1002001003 = I_ERI_Gx3z_S_Fy2z_S_C1002001003+CDY*I_ERI_Gx3z_S_D2z_S_C1002001003;
  Double I_ERI_G4y_S_D2z_Py_C1002001003 = I_ERI_G4y_S_Fy2z_S_C1002001003+CDY*I_ERI_G4y_S_D2z_S_C1002001003;
  Double I_ERI_G3yz_S_D2z_Py_C1002001003 = I_ERI_G3yz_S_Fy2z_S_C1002001003+CDY*I_ERI_G3yz_S_D2z_S_C1002001003;
  Double I_ERI_G2y2z_S_D2z_Py_C1002001003 = I_ERI_G2y2z_S_Fy2z_S_C1002001003+CDY*I_ERI_G2y2z_S_D2z_S_C1002001003;
  Double I_ERI_Gy3z_S_D2z_Py_C1002001003 = I_ERI_Gy3z_S_Fy2z_S_C1002001003+CDY*I_ERI_Gy3z_S_D2z_S_C1002001003;
  Double I_ERI_G4z_S_D2z_Py_C1002001003 = I_ERI_G4z_S_Fy2z_S_C1002001003+CDY*I_ERI_G4z_S_D2z_S_C1002001003;
  Double I_ERI_G4x_S_D2x_Pz_C1002001003 = I_ERI_G4x_S_F2xz_S_C1002001003+CDZ*I_ERI_G4x_S_D2x_S_C1002001003;
  Double I_ERI_G3xy_S_D2x_Pz_C1002001003 = I_ERI_G3xy_S_F2xz_S_C1002001003+CDZ*I_ERI_G3xy_S_D2x_S_C1002001003;
  Double I_ERI_G3xz_S_D2x_Pz_C1002001003 = I_ERI_G3xz_S_F2xz_S_C1002001003+CDZ*I_ERI_G3xz_S_D2x_S_C1002001003;
  Double I_ERI_G2x2y_S_D2x_Pz_C1002001003 = I_ERI_G2x2y_S_F2xz_S_C1002001003+CDZ*I_ERI_G2x2y_S_D2x_S_C1002001003;
  Double I_ERI_G2xyz_S_D2x_Pz_C1002001003 = I_ERI_G2xyz_S_F2xz_S_C1002001003+CDZ*I_ERI_G2xyz_S_D2x_S_C1002001003;
  Double I_ERI_G2x2z_S_D2x_Pz_C1002001003 = I_ERI_G2x2z_S_F2xz_S_C1002001003+CDZ*I_ERI_G2x2z_S_D2x_S_C1002001003;
  Double I_ERI_Gx3y_S_D2x_Pz_C1002001003 = I_ERI_Gx3y_S_F2xz_S_C1002001003+CDZ*I_ERI_Gx3y_S_D2x_S_C1002001003;
  Double I_ERI_Gx2yz_S_D2x_Pz_C1002001003 = I_ERI_Gx2yz_S_F2xz_S_C1002001003+CDZ*I_ERI_Gx2yz_S_D2x_S_C1002001003;
  Double I_ERI_Gxy2z_S_D2x_Pz_C1002001003 = I_ERI_Gxy2z_S_F2xz_S_C1002001003+CDZ*I_ERI_Gxy2z_S_D2x_S_C1002001003;
  Double I_ERI_Gx3z_S_D2x_Pz_C1002001003 = I_ERI_Gx3z_S_F2xz_S_C1002001003+CDZ*I_ERI_Gx3z_S_D2x_S_C1002001003;
  Double I_ERI_G4y_S_D2x_Pz_C1002001003 = I_ERI_G4y_S_F2xz_S_C1002001003+CDZ*I_ERI_G4y_S_D2x_S_C1002001003;
  Double I_ERI_G3yz_S_D2x_Pz_C1002001003 = I_ERI_G3yz_S_F2xz_S_C1002001003+CDZ*I_ERI_G3yz_S_D2x_S_C1002001003;
  Double I_ERI_G2y2z_S_D2x_Pz_C1002001003 = I_ERI_G2y2z_S_F2xz_S_C1002001003+CDZ*I_ERI_G2y2z_S_D2x_S_C1002001003;
  Double I_ERI_Gy3z_S_D2x_Pz_C1002001003 = I_ERI_Gy3z_S_F2xz_S_C1002001003+CDZ*I_ERI_Gy3z_S_D2x_S_C1002001003;
  Double I_ERI_G4z_S_D2x_Pz_C1002001003 = I_ERI_G4z_S_F2xz_S_C1002001003+CDZ*I_ERI_G4z_S_D2x_S_C1002001003;
  Double I_ERI_G4x_S_Dxy_Pz_C1002001003 = I_ERI_G4x_S_Fxyz_S_C1002001003+CDZ*I_ERI_G4x_S_Dxy_S_C1002001003;
  Double I_ERI_G3xy_S_Dxy_Pz_C1002001003 = I_ERI_G3xy_S_Fxyz_S_C1002001003+CDZ*I_ERI_G3xy_S_Dxy_S_C1002001003;
  Double I_ERI_G3xz_S_Dxy_Pz_C1002001003 = I_ERI_G3xz_S_Fxyz_S_C1002001003+CDZ*I_ERI_G3xz_S_Dxy_S_C1002001003;
  Double I_ERI_G2x2y_S_Dxy_Pz_C1002001003 = I_ERI_G2x2y_S_Fxyz_S_C1002001003+CDZ*I_ERI_G2x2y_S_Dxy_S_C1002001003;
  Double I_ERI_G2xyz_S_Dxy_Pz_C1002001003 = I_ERI_G2xyz_S_Fxyz_S_C1002001003+CDZ*I_ERI_G2xyz_S_Dxy_S_C1002001003;
  Double I_ERI_G2x2z_S_Dxy_Pz_C1002001003 = I_ERI_G2x2z_S_Fxyz_S_C1002001003+CDZ*I_ERI_G2x2z_S_Dxy_S_C1002001003;
  Double I_ERI_Gx3y_S_Dxy_Pz_C1002001003 = I_ERI_Gx3y_S_Fxyz_S_C1002001003+CDZ*I_ERI_Gx3y_S_Dxy_S_C1002001003;
  Double I_ERI_Gx2yz_S_Dxy_Pz_C1002001003 = I_ERI_Gx2yz_S_Fxyz_S_C1002001003+CDZ*I_ERI_Gx2yz_S_Dxy_S_C1002001003;
  Double I_ERI_Gxy2z_S_Dxy_Pz_C1002001003 = I_ERI_Gxy2z_S_Fxyz_S_C1002001003+CDZ*I_ERI_Gxy2z_S_Dxy_S_C1002001003;
  Double I_ERI_Gx3z_S_Dxy_Pz_C1002001003 = I_ERI_Gx3z_S_Fxyz_S_C1002001003+CDZ*I_ERI_Gx3z_S_Dxy_S_C1002001003;
  Double I_ERI_G4y_S_Dxy_Pz_C1002001003 = I_ERI_G4y_S_Fxyz_S_C1002001003+CDZ*I_ERI_G4y_S_Dxy_S_C1002001003;
  Double I_ERI_G3yz_S_Dxy_Pz_C1002001003 = I_ERI_G3yz_S_Fxyz_S_C1002001003+CDZ*I_ERI_G3yz_S_Dxy_S_C1002001003;
  Double I_ERI_G2y2z_S_Dxy_Pz_C1002001003 = I_ERI_G2y2z_S_Fxyz_S_C1002001003+CDZ*I_ERI_G2y2z_S_Dxy_S_C1002001003;
  Double I_ERI_Gy3z_S_Dxy_Pz_C1002001003 = I_ERI_Gy3z_S_Fxyz_S_C1002001003+CDZ*I_ERI_Gy3z_S_Dxy_S_C1002001003;
  Double I_ERI_G4z_S_Dxy_Pz_C1002001003 = I_ERI_G4z_S_Fxyz_S_C1002001003+CDZ*I_ERI_G4z_S_Dxy_S_C1002001003;
  Double I_ERI_G4x_S_Dxz_Pz_C1002001003 = I_ERI_G4x_S_Fx2z_S_C1002001003+CDZ*I_ERI_G4x_S_Dxz_S_C1002001003;
  Double I_ERI_G3xy_S_Dxz_Pz_C1002001003 = I_ERI_G3xy_S_Fx2z_S_C1002001003+CDZ*I_ERI_G3xy_S_Dxz_S_C1002001003;
  Double I_ERI_G3xz_S_Dxz_Pz_C1002001003 = I_ERI_G3xz_S_Fx2z_S_C1002001003+CDZ*I_ERI_G3xz_S_Dxz_S_C1002001003;
  Double I_ERI_G2x2y_S_Dxz_Pz_C1002001003 = I_ERI_G2x2y_S_Fx2z_S_C1002001003+CDZ*I_ERI_G2x2y_S_Dxz_S_C1002001003;
  Double I_ERI_G2xyz_S_Dxz_Pz_C1002001003 = I_ERI_G2xyz_S_Fx2z_S_C1002001003+CDZ*I_ERI_G2xyz_S_Dxz_S_C1002001003;
  Double I_ERI_G2x2z_S_Dxz_Pz_C1002001003 = I_ERI_G2x2z_S_Fx2z_S_C1002001003+CDZ*I_ERI_G2x2z_S_Dxz_S_C1002001003;
  Double I_ERI_Gx3y_S_Dxz_Pz_C1002001003 = I_ERI_Gx3y_S_Fx2z_S_C1002001003+CDZ*I_ERI_Gx3y_S_Dxz_S_C1002001003;
  Double I_ERI_Gx2yz_S_Dxz_Pz_C1002001003 = I_ERI_Gx2yz_S_Fx2z_S_C1002001003+CDZ*I_ERI_Gx2yz_S_Dxz_S_C1002001003;
  Double I_ERI_Gxy2z_S_Dxz_Pz_C1002001003 = I_ERI_Gxy2z_S_Fx2z_S_C1002001003+CDZ*I_ERI_Gxy2z_S_Dxz_S_C1002001003;
  Double I_ERI_Gx3z_S_Dxz_Pz_C1002001003 = I_ERI_Gx3z_S_Fx2z_S_C1002001003+CDZ*I_ERI_Gx3z_S_Dxz_S_C1002001003;
  Double I_ERI_G4y_S_Dxz_Pz_C1002001003 = I_ERI_G4y_S_Fx2z_S_C1002001003+CDZ*I_ERI_G4y_S_Dxz_S_C1002001003;
  Double I_ERI_G3yz_S_Dxz_Pz_C1002001003 = I_ERI_G3yz_S_Fx2z_S_C1002001003+CDZ*I_ERI_G3yz_S_Dxz_S_C1002001003;
  Double I_ERI_G2y2z_S_Dxz_Pz_C1002001003 = I_ERI_G2y2z_S_Fx2z_S_C1002001003+CDZ*I_ERI_G2y2z_S_Dxz_S_C1002001003;
  Double I_ERI_Gy3z_S_Dxz_Pz_C1002001003 = I_ERI_Gy3z_S_Fx2z_S_C1002001003+CDZ*I_ERI_Gy3z_S_Dxz_S_C1002001003;
  Double I_ERI_G4z_S_Dxz_Pz_C1002001003 = I_ERI_G4z_S_Fx2z_S_C1002001003+CDZ*I_ERI_G4z_S_Dxz_S_C1002001003;
  Double I_ERI_G4x_S_D2y_Pz_C1002001003 = I_ERI_G4x_S_F2yz_S_C1002001003+CDZ*I_ERI_G4x_S_D2y_S_C1002001003;
  Double I_ERI_G3xy_S_D2y_Pz_C1002001003 = I_ERI_G3xy_S_F2yz_S_C1002001003+CDZ*I_ERI_G3xy_S_D2y_S_C1002001003;
  Double I_ERI_G3xz_S_D2y_Pz_C1002001003 = I_ERI_G3xz_S_F2yz_S_C1002001003+CDZ*I_ERI_G3xz_S_D2y_S_C1002001003;
  Double I_ERI_G2x2y_S_D2y_Pz_C1002001003 = I_ERI_G2x2y_S_F2yz_S_C1002001003+CDZ*I_ERI_G2x2y_S_D2y_S_C1002001003;
  Double I_ERI_G2xyz_S_D2y_Pz_C1002001003 = I_ERI_G2xyz_S_F2yz_S_C1002001003+CDZ*I_ERI_G2xyz_S_D2y_S_C1002001003;
  Double I_ERI_G2x2z_S_D2y_Pz_C1002001003 = I_ERI_G2x2z_S_F2yz_S_C1002001003+CDZ*I_ERI_G2x2z_S_D2y_S_C1002001003;
  Double I_ERI_Gx3y_S_D2y_Pz_C1002001003 = I_ERI_Gx3y_S_F2yz_S_C1002001003+CDZ*I_ERI_Gx3y_S_D2y_S_C1002001003;
  Double I_ERI_Gx2yz_S_D2y_Pz_C1002001003 = I_ERI_Gx2yz_S_F2yz_S_C1002001003+CDZ*I_ERI_Gx2yz_S_D2y_S_C1002001003;
  Double I_ERI_Gxy2z_S_D2y_Pz_C1002001003 = I_ERI_Gxy2z_S_F2yz_S_C1002001003+CDZ*I_ERI_Gxy2z_S_D2y_S_C1002001003;
  Double I_ERI_Gx3z_S_D2y_Pz_C1002001003 = I_ERI_Gx3z_S_F2yz_S_C1002001003+CDZ*I_ERI_Gx3z_S_D2y_S_C1002001003;
  Double I_ERI_G4y_S_D2y_Pz_C1002001003 = I_ERI_G4y_S_F2yz_S_C1002001003+CDZ*I_ERI_G4y_S_D2y_S_C1002001003;
  Double I_ERI_G3yz_S_D2y_Pz_C1002001003 = I_ERI_G3yz_S_F2yz_S_C1002001003+CDZ*I_ERI_G3yz_S_D2y_S_C1002001003;
  Double I_ERI_G2y2z_S_D2y_Pz_C1002001003 = I_ERI_G2y2z_S_F2yz_S_C1002001003+CDZ*I_ERI_G2y2z_S_D2y_S_C1002001003;
  Double I_ERI_Gy3z_S_D2y_Pz_C1002001003 = I_ERI_Gy3z_S_F2yz_S_C1002001003+CDZ*I_ERI_Gy3z_S_D2y_S_C1002001003;
  Double I_ERI_G4z_S_D2y_Pz_C1002001003 = I_ERI_G4z_S_F2yz_S_C1002001003+CDZ*I_ERI_G4z_S_D2y_S_C1002001003;
  Double I_ERI_G4x_S_Dyz_Pz_C1002001003 = I_ERI_G4x_S_Fy2z_S_C1002001003+CDZ*I_ERI_G4x_S_Dyz_S_C1002001003;
  Double I_ERI_G3xy_S_Dyz_Pz_C1002001003 = I_ERI_G3xy_S_Fy2z_S_C1002001003+CDZ*I_ERI_G3xy_S_Dyz_S_C1002001003;
  Double I_ERI_G3xz_S_Dyz_Pz_C1002001003 = I_ERI_G3xz_S_Fy2z_S_C1002001003+CDZ*I_ERI_G3xz_S_Dyz_S_C1002001003;
  Double I_ERI_G2x2y_S_Dyz_Pz_C1002001003 = I_ERI_G2x2y_S_Fy2z_S_C1002001003+CDZ*I_ERI_G2x2y_S_Dyz_S_C1002001003;
  Double I_ERI_G2xyz_S_Dyz_Pz_C1002001003 = I_ERI_G2xyz_S_Fy2z_S_C1002001003+CDZ*I_ERI_G2xyz_S_Dyz_S_C1002001003;
  Double I_ERI_G2x2z_S_Dyz_Pz_C1002001003 = I_ERI_G2x2z_S_Fy2z_S_C1002001003+CDZ*I_ERI_G2x2z_S_Dyz_S_C1002001003;
  Double I_ERI_Gx3y_S_Dyz_Pz_C1002001003 = I_ERI_Gx3y_S_Fy2z_S_C1002001003+CDZ*I_ERI_Gx3y_S_Dyz_S_C1002001003;
  Double I_ERI_Gx2yz_S_Dyz_Pz_C1002001003 = I_ERI_Gx2yz_S_Fy2z_S_C1002001003+CDZ*I_ERI_Gx2yz_S_Dyz_S_C1002001003;
  Double I_ERI_Gxy2z_S_Dyz_Pz_C1002001003 = I_ERI_Gxy2z_S_Fy2z_S_C1002001003+CDZ*I_ERI_Gxy2z_S_Dyz_S_C1002001003;
  Double I_ERI_Gx3z_S_Dyz_Pz_C1002001003 = I_ERI_Gx3z_S_Fy2z_S_C1002001003+CDZ*I_ERI_Gx3z_S_Dyz_S_C1002001003;
  Double I_ERI_G4y_S_Dyz_Pz_C1002001003 = I_ERI_G4y_S_Fy2z_S_C1002001003+CDZ*I_ERI_G4y_S_Dyz_S_C1002001003;
  Double I_ERI_G3yz_S_Dyz_Pz_C1002001003 = I_ERI_G3yz_S_Fy2z_S_C1002001003+CDZ*I_ERI_G3yz_S_Dyz_S_C1002001003;
  Double I_ERI_G2y2z_S_Dyz_Pz_C1002001003 = I_ERI_G2y2z_S_Fy2z_S_C1002001003+CDZ*I_ERI_G2y2z_S_Dyz_S_C1002001003;
  Double I_ERI_Gy3z_S_Dyz_Pz_C1002001003 = I_ERI_Gy3z_S_Fy2z_S_C1002001003+CDZ*I_ERI_Gy3z_S_Dyz_S_C1002001003;
  Double I_ERI_G4z_S_Dyz_Pz_C1002001003 = I_ERI_G4z_S_Fy2z_S_C1002001003+CDZ*I_ERI_G4z_S_Dyz_S_C1002001003;
  Double I_ERI_G4x_S_D2z_Pz_C1002001003 = I_ERI_G4x_S_F3z_S_C1002001003+CDZ*I_ERI_G4x_S_D2z_S_C1002001003;
  Double I_ERI_G3xy_S_D2z_Pz_C1002001003 = I_ERI_G3xy_S_F3z_S_C1002001003+CDZ*I_ERI_G3xy_S_D2z_S_C1002001003;
  Double I_ERI_G3xz_S_D2z_Pz_C1002001003 = I_ERI_G3xz_S_F3z_S_C1002001003+CDZ*I_ERI_G3xz_S_D2z_S_C1002001003;
  Double I_ERI_G2x2y_S_D2z_Pz_C1002001003 = I_ERI_G2x2y_S_F3z_S_C1002001003+CDZ*I_ERI_G2x2y_S_D2z_S_C1002001003;
  Double I_ERI_G2xyz_S_D2z_Pz_C1002001003 = I_ERI_G2xyz_S_F3z_S_C1002001003+CDZ*I_ERI_G2xyz_S_D2z_S_C1002001003;
  Double I_ERI_G2x2z_S_D2z_Pz_C1002001003 = I_ERI_G2x2z_S_F3z_S_C1002001003+CDZ*I_ERI_G2x2z_S_D2z_S_C1002001003;
  Double I_ERI_Gx3y_S_D2z_Pz_C1002001003 = I_ERI_Gx3y_S_F3z_S_C1002001003+CDZ*I_ERI_Gx3y_S_D2z_S_C1002001003;
  Double I_ERI_Gx2yz_S_D2z_Pz_C1002001003 = I_ERI_Gx2yz_S_F3z_S_C1002001003+CDZ*I_ERI_Gx2yz_S_D2z_S_C1002001003;
  Double I_ERI_Gxy2z_S_D2z_Pz_C1002001003 = I_ERI_Gxy2z_S_F3z_S_C1002001003+CDZ*I_ERI_Gxy2z_S_D2z_S_C1002001003;
  Double I_ERI_Gx3z_S_D2z_Pz_C1002001003 = I_ERI_Gx3z_S_F3z_S_C1002001003+CDZ*I_ERI_Gx3z_S_D2z_S_C1002001003;
  Double I_ERI_G4y_S_D2z_Pz_C1002001003 = I_ERI_G4y_S_F3z_S_C1002001003+CDZ*I_ERI_G4y_S_D2z_S_C1002001003;
  Double I_ERI_G3yz_S_D2z_Pz_C1002001003 = I_ERI_G3yz_S_F3z_S_C1002001003+CDZ*I_ERI_G3yz_S_D2z_S_C1002001003;
  Double I_ERI_G2y2z_S_D2z_Pz_C1002001003 = I_ERI_G2y2z_S_F3z_S_C1002001003+CDZ*I_ERI_G2y2z_S_D2z_S_C1002001003;
  Double I_ERI_Gy3z_S_D2z_Pz_C1002001003 = I_ERI_Gy3z_S_F3z_S_C1002001003+CDZ*I_ERI_Gy3z_S_D2z_S_C1002001003;
  Double I_ERI_G4z_S_D2z_Pz_C1002001003 = I_ERI_G4z_S_F3z_S_C1002001003+CDZ*I_ERI_G4z_S_D2z_S_C1002001003;

  /************************************************************
   * declare the HRR2 result shell quartets in array form
   ************************************************************/

  /************************************************************
   * initilize the HRR steps : build the AB/CD variables
   ************************************************************/
  Double ABX = A[0] - B[0];
  Double ABY = A[1] - B[1];
  Double ABZ = A[2] - B[2];

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_S_C2001003
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S_C2001003
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C2001003
   ************************************************************/
  abcd[10] = I_ERI_G4x_S_D2x_S_C2001003+ABX*I_ERI_F3x_S_D2x_S_C2001003;
  abcd[11] = I_ERI_G3xy_S_D2x_S_C2001003+ABX*I_ERI_F2xy_S_D2x_S_C2001003;
  abcd[12] = I_ERI_G3xz_S_D2x_S_C2001003+ABX*I_ERI_F2xz_S_D2x_S_C2001003;
  abcd[13] = I_ERI_G2x2y_S_D2x_S_C2001003+ABX*I_ERI_Fx2y_S_D2x_S_C2001003;
  abcd[14] = I_ERI_G2xyz_S_D2x_S_C2001003+ABX*I_ERI_Fxyz_S_D2x_S_C2001003;
  abcd[15] = I_ERI_G2x2z_S_D2x_S_C2001003+ABX*I_ERI_Fx2z_S_D2x_S_C2001003;
  abcd[16] = I_ERI_Gx3y_S_D2x_S_C2001003+ABX*I_ERI_F3y_S_D2x_S_C2001003;
  abcd[17] = I_ERI_Gx2yz_S_D2x_S_C2001003+ABX*I_ERI_F2yz_S_D2x_S_C2001003;
  abcd[18] = I_ERI_Gxy2z_S_D2x_S_C2001003+ABX*I_ERI_Fy2z_S_D2x_S_C2001003;
  abcd[19] = I_ERI_Gx3z_S_D2x_S_C2001003+ABX*I_ERI_F3z_S_D2x_S_C2001003;
  abcd[20] = I_ERI_G3xy_S_D2x_S_C2001003+ABY*I_ERI_F3x_S_D2x_S_C2001003;
  abcd[21] = I_ERI_G2x2y_S_D2x_S_C2001003+ABY*I_ERI_F2xy_S_D2x_S_C2001003;
  abcd[22] = I_ERI_G2xyz_S_D2x_S_C2001003+ABY*I_ERI_F2xz_S_D2x_S_C2001003;
  abcd[23] = I_ERI_Gx3y_S_D2x_S_C2001003+ABY*I_ERI_Fx2y_S_D2x_S_C2001003;
  abcd[24] = I_ERI_Gx2yz_S_D2x_S_C2001003+ABY*I_ERI_Fxyz_S_D2x_S_C2001003;
  abcd[25] = I_ERI_Gxy2z_S_D2x_S_C2001003+ABY*I_ERI_Fx2z_S_D2x_S_C2001003;
  abcd[26] = I_ERI_G4y_S_D2x_S_C2001003+ABY*I_ERI_F3y_S_D2x_S_C2001003;
  abcd[27] = I_ERI_G3yz_S_D2x_S_C2001003+ABY*I_ERI_F2yz_S_D2x_S_C2001003;
  abcd[28] = I_ERI_G2y2z_S_D2x_S_C2001003+ABY*I_ERI_Fy2z_S_D2x_S_C2001003;
  abcd[29] = I_ERI_Gy3z_S_D2x_S_C2001003+ABY*I_ERI_F3z_S_D2x_S_C2001003;
  abcd[30] = I_ERI_G3xz_S_D2x_S_C2001003+ABZ*I_ERI_F3x_S_D2x_S_C2001003;
  abcd[31] = I_ERI_G2xyz_S_D2x_S_C2001003+ABZ*I_ERI_F2xy_S_D2x_S_C2001003;
  abcd[32] = I_ERI_G2x2z_S_D2x_S_C2001003+ABZ*I_ERI_F2xz_S_D2x_S_C2001003;
  abcd[33] = I_ERI_Gx2yz_S_D2x_S_C2001003+ABZ*I_ERI_Fx2y_S_D2x_S_C2001003;
  abcd[34] = I_ERI_Gxy2z_S_D2x_S_C2001003+ABZ*I_ERI_Fxyz_S_D2x_S_C2001003;
  abcd[35] = I_ERI_Gx3z_S_D2x_S_C2001003+ABZ*I_ERI_Fx2z_S_D2x_S_C2001003;
  abcd[36] = I_ERI_G3yz_S_D2x_S_C2001003+ABZ*I_ERI_F3y_S_D2x_S_C2001003;
  abcd[37] = I_ERI_G2y2z_S_D2x_S_C2001003+ABZ*I_ERI_F2yz_S_D2x_S_C2001003;
  abcd[38] = I_ERI_Gy3z_S_D2x_S_C2001003+ABZ*I_ERI_Fy2z_S_D2x_S_C2001003;
  abcd[39] = I_ERI_G4z_S_D2x_S_C2001003+ABZ*I_ERI_F3z_S_D2x_S_C2001003;
  abcd[50] = I_ERI_G4x_S_Dxy_S_C2001003+ABX*I_ERI_F3x_S_Dxy_S_C2001003;
  abcd[51] = I_ERI_G3xy_S_Dxy_S_C2001003+ABX*I_ERI_F2xy_S_Dxy_S_C2001003;
  abcd[52] = I_ERI_G3xz_S_Dxy_S_C2001003+ABX*I_ERI_F2xz_S_Dxy_S_C2001003;
  abcd[53] = I_ERI_G2x2y_S_Dxy_S_C2001003+ABX*I_ERI_Fx2y_S_Dxy_S_C2001003;
  abcd[54] = I_ERI_G2xyz_S_Dxy_S_C2001003+ABX*I_ERI_Fxyz_S_Dxy_S_C2001003;
  abcd[55] = I_ERI_G2x2z_S_Dxy_S_C2001003+ABX*I_ERI_Fx2z_S_Dxy_S_C2001003;
  abcd[56] = I_ERI_Gx3y_S_Dxy_S_C2001003+ABX*I_ERI_F3y_S_Dxy_S_C2001003;
  abcd[57] = I_ERI_Gx2yz_S_Dxy_S_C2001003+ABX*I_ERI_F2yz_S_Dxy_S_C2001003;
  abcd[58] = I_ERI_Gxy2z_S_Dxy_S_C2001003+ABX*I_ERI_Fy2z_S_Dxy_S_C2001003;
  abcd[59] = I_ERI_Gx3z_S_Dxy_S_C2001003+ABX*I_ERI_F3z_S_Dxy_S_C2001003;
  abcd[60] = I_ERI_G3xy_S_Dxy_S_C2001003+ABY*I_ERI_F3x_S_Dxy_S_C2001003;
  abcd[61] = I_ERI_G2x2y_S_Dxy_S_C2001003+ABY*I_ERI_F2xy_S_Dxy_S_C2001003;
  abcd[62] = I_ERI_G2xyz_S_Dxy_S_C2001003+ABY*I_ERI_F2xz_S_Dxy_S_C2001003;
  abcd[63] = I_ERI_Gx3y_S_Dxy_S_C2001003+ABY*I_ERI_Fx2y_S_Dxy_S_C2001003;
  abcd[64] = I_ERI_Gx2yz_S_Dxy_S_C2001003+ABY*I_ERI_Fxyz_S_Dxy_S_C2001003;
  abcd[65] = I_ERI_Gxy2z_S_Dxy_S_C2001003+ABY*I_ERI_Fx2z_S_Dxy_S_C2001003;
  abcd[66] = I_ERI_G4y_S_Dxy_S_C2001003+ABY*I_ERI_F3y_S_Dxy_S_C2001003;
  abcd[67] = I_ERI_G3yz_S_Dxy_S_C2001003+ABY*I_ERI_F2yz_S_Dxy_S_C2001003;
  abcd[68] = I_ERI_G2y2z_S_Dxy_S_C2001003+ABY*I_ERI_Fy2z_S_Dxy_S_C2001003;
  abcd[69] = I_ERI_Gy3z_S_Dxy_S_C2001003+ABY*I_ERI_F3z_S_Dxy_S_C2001003;
  abcd[70] = I_ERI_G3xz_S_Dxy_S_C2001003+ABZ*I_ERI_F3x_S_Dxy_S_C2001003;
  abcd[71] = I_ERI_G2xyz_S_Dxy_S_C2001003+ABZ*I_ERI_F2xy_S_Dxy_S_C2001003;
  abcd[72] = I_ERI_G2x2z_S_Dxy_S_C2001003+ABZ*I_ERI_F2xz_S_Dxy_S_C2001003;
  abcd[73] = I_ERI_Gx2yz_S_Dxy_S_C2001003+ABZ*I_ERI_Fx2y_S_Dxy_S_C2001003;
  abcd[74] = I_ERI_Gxy2z_S_Dxy_S_C2001003+ABZ*I_ERI_Fxyz_S_Dxy_S_C2001003;
  abcd[75] = I_ERI_Gx3z_S_Dxy_S_C2001003+ABZ*I_ERI_Fx2z_S_Dxy_S_C2001003;
  abcd[76] = I_ERI_G3yz_S_Dxy_S_C2001003+ABZ*I_ERI_F3y_S_Dxy_S_C2001003;
  abcd[77] = I_ERI_G2y2z_S_Dxy_S_C2001003+ABZ*I_ERI_F2yz_S_Dxy_S_C2001003;
  abcd[78] = I_ERI_Gy3z_S_Dxy_S_C2001003+ABZ*I_ERI_Fy2z_S_Dxy_S_C2001003;
  abcd[79] = I_ERI_G4z_S_Dxy_S_C2001003+ABZ*I_ERI_F3z_S_Dxy_S_C2001003;
  abcd[90] = I_ERI_G4x_S_Dxz_S_C2001003+ABX*I_ERI_F3x_S_Dxz_S_C2001003;
  abcd[91] = I_ERI_G3xy_S_Dxz_S_C2001003+ABX*I_ERI_F2xy_S_Dxz_S_C2001003;
  abcd[92] = I_ERI_G3xz_S_Dxz_S_C2001003+ABX*I_ERI_F2xz_S_Dxz_S_C2001003;
  abcd[93] = I_ERI_G2x2y_S_Dxz_S_C2001003+ABX*I_ERI_Fx2y_S_Dxz_S_C2001003;
  abcd[94] = I_ERI_G2xyz_S_Dxz_S_C2001003+ABX*I_ERI_Fxyz_S_Dxz_S_C2001003;
  abcd[95] = I_ERI_G2x2z_S_Dxz_S_C2001003+ABX*I_ERI_Fx2z_S_Dxz_S_C2001003;
  abcd[96] = I_ERI_Gx3y_S_Dxz_S_C2001003+ABX*I_ERI_F3y_S_Dxz_S_C2001003;
  abcd[97] = I_ERI_Gx2yz_S_Dxz_S_C2001003+ABX*I_ERI_F2yz_S_Dxz_S_C2001003;
  abcd[98] = I_ERI_Gxy2z_S_Dxz_S_C2001003+ABX*I_ERI_Fy2z_S_Dxz_S_C2001003;
  abcd[99] = I_ERI_Gx3z_S_Dxz_S_C2001003+ABX*I_ERI_F3z_S_Dxz_S_C2001003;
  abcd[100] = I_ERI_G3xy_S_Dxz_S_C2001003+ABY*I_ERI_F3x_S_Dxz_S_C2001003;
  abcd[101] = I_ERI_G2x2y_S_Dxz_S_C2001003+ABY*I_ERI_F2xy_S_Dxz_S_C2001003;
  abcd[102] = I_ERI_G2xyz_S_Dxz_S_C2001003+ABY*I_ERI_F2xz_S_Dxz_S_C2001003;
  abcd[103] = I_ERI_Gx3y_S_Dxz_S_C2001003+ABY*I_ERI_Fx2y_S_Dxz_S_C2001003;
  abcd[104] = I_ERI_Gx2yz_S_Dxz_S_C2001003+ABY*I_ERI_Fxyz_S_Dxz_S_C2001003;
  abcd[105] = I_ERI_Gxy2z_S_Dxz_S_C2001003+ABY*I_ERI_Fx2z_S_Dxz_S_C2001003;
  abcd[106] = I_ERI_G4y_S_Dxz_S_C2001003+ABY*I_ERI_F3y_S_Dxz_S_C2001003;
  abcd[107] = I_ERI_G3yz_S_Dxz_S_C2001003+ABY*I_ERI_F2yz_S_Dxz_S_C2001003;
  abcd[108] = I_ERI_G2y2z_S_Dxz_S_C2001003+ABY*I_ERI_Fy2z_S_Dxz_S_C2001003;
  abcd[109] = I_ERI_Gy3z_S_Dxz_S_C2001003+ABY*I_ERI_F3z_S_Dxz_S_C2001003;
  abcd[110] = I_ERI_G3xz_S_Dxz_S_C2001003+ABZ*I_ERI_F3x_S_Dxz_S_C2001003;
  abcd[111] = I_ERI_G2xyz_S_Dxz_S_C2001003+ABZ*I_ERI_F2xy_S_Dxz_S_C2001003;
  abcd[112] = I_ERI_G2x2z_S_Dxz_S_C2001003+ABZ*I_ERI_F2xz_S_Dxz_S_C2001003;
  abcd[113] = I_ERI_Gx2yz_S_Dxz_S_C2001003+ABZ*I_ERI_Fx2y_S_Dxz_S_C2001003;
  abcd[114] = I_ERI_Gxy2z_S_Dxz_S_C2001003+ABZ*I_ERI_Fxyz_S_Dxz_S_C2001003;
  abcd[115] = I_ERI_Gx3z_S_Dxz_S_C2001003+ABZ*I_ERI_Fx2z_S_Dxz_S_C2001003;
  abcd[116] = I_ERI_G3yz_S_Dxz_S_C2001003+ABZ*I_ERI_F3y_S_Dxz_S_C2001003;
  abcd[117] = I_ERI_G2y2z_S_Dxz_S_C2001003+ABZ*I_ERI_F2yz_S_Dxz_S_C2001003;
  abcd[118] = I_ERI_Gy3z_S_Dxz_S_C2001003+ABZ*I_ERI_Fy2z_S_Dxz_S_C2001003;
  abcd[119] = I_ERI_G4z_S_Dxz_S_C2001003+ABZ*I_ERI_F3z_S_Dxz_S_C2001003;
  abcd[130] = I_ERI_G4x_S_D2y_S_C2001003+ABX*I_ERI_F3x_S_D2y_S_C2001003;
  abcd[131] = I_ERI_G3xy_S_D2y_S_C2001003+ABX*I_ERI_F2xy_S_D2y_S_C2001003;
  abcd[132] = I_ERI_G3xz_S_D2y_S_C2001003+ABX*I_ERI_F2xz_S_D2y_S_C2001003;
  abcd[133] = I_ERI_G2x2y_S_D2y_S_C2001003+ABX*I_ERI_Fx2y_S_D2y_S_C2001003;
  abcd[134] = I_ERI_G2xyz_S_D2y_S_C2001003+ABX*I_ERI_Fxyz_S_D2y_S_C2001003;
  abcd[135] = I_ERI_G2x2z_S_D2y_S_C2001003+ABX*I_ERI_Fx2z_S_D2y_S_C2001003;
  abcd[136] = I_ERI_Gx3y_S_D2y_S_C2001003+ABX*I_ERI_F3y_S_D2y_S_C2001003;
  abcd[137] = I_ERI_Gx2yz_S_D2y_S_C2001003+ABX*I_ERI_F2yz_S_D2y_S_C2001003;
  abcd[138] = I_ERI_Gxy2z_S_D2y_S_C2001003+ABX*I_ERI_Fy2z_S_D2y_S_C2001003;
  abcd[139] = I_ERI_Gx3z_S_D2y_S_C2001003+ABX*I_ERI_F3z_S_D2y_S_C2001003;
  abcd[140] = I_ERI_G3xy_S_D2y_S_C2001003+ABY*I_ERI_F3x_S_D2y_S_C2001003;
  abcd[141] = I_ERI_G2x2y_S_D2y_S_C2001003+ABY*I_ERI_F2xy_S_D2y_S_C2001003;
  abcd[142] = I_ERI_G2xyz_S_D2y_S_C2001003+ABY*I_ERI_F2xz_S_D2y_S_C2001003;
  abcd[143] = I_ERI_Gx3y_S_D2y_S_C2001003+ABY*I_ERI_Fx2y_S_D2y_S_C2001003;
  abcd[144] = I_ERI_Gx2yz_S_D2y_S_C2001003+ABY*I_ERI_Fxyz_S_D2y_S_C2001003;
  abcd[145] = I_ERI_Gxy2z_S_D2y_S_C2001003+ABY*I_ERI_Fx2z_S_D2y_S_C2001003;
  abcd[146] = I_ERI_G4y_S_D2y_S_C2001003+ABY*I_ERI_F3y_S_D2y_S_C2001003;
  abcd[147] = I_ERI_G3yz_S_D2y_S_C2001003+ABY*I_ERI_F2yz_S_D2y_S_C2001003;
  abcd[148] = I_ERI_G2y2z_S_D2y_S_C2001003+ABY*I_ERI_Fy2z_S_D2y_S_C2001003;
  abcd[149] = I_ERI_Gy3z_S_D2y_S_C2001003+ABY*I_ERI_F3z_S_D2y_S_C2001003;
  abcd[150] = I_ERI_G3xz_S_D2y_S_C2001003+ABZ*I_ERI_F3x_S_D2y_S_C2001003;
  abcd[151] = I_ERI_G2xyz_S_D2y_S_C2001003+ABZ*I_ERI_F2xy_S_D2y_S_C2001003;
  abcd[152] = I_ERI_G2x2z_S_D2y_S_C2001003+ABZ*I_ERI_F2xz_S_D2y_S_C2001003;
  abcd[153] = I_ERI_Gx2yz_S_D2y_S_C2001003+ABZ*I_ERI_Fx2y_S_D2y_S_C2001003;
  abcd[154] = I_ERI_Gxy2z_S_D2y_S_C2001003+ABZ*I_ERI_Fxyz_S_D2y_S_C2001003;
  abcd[155] = I_ERI_Gx3z_S_D2y_S_C2001003+ABZ*I_ERI_Fx2z_S_D2y_S_C2001003;
  abcd[156] = I_ERI_G3yz_S_D2y_S_C2001003+ABZ*I_ERI_F3y_S_D2y_S_C2001003;
  abcd[157] = I_ERI_G2y2z_S_D2y_S_C2001003+ABZ*I_ERI_F2yz_S_D2y_S_C2001003;
  abcd[158] = I_ERI_Gy3z_S_D2y_S_C2001003+ABZ*I_ERI_Fy2z_S_D2y_S_C2001003;
  abcd[159] = I_ERI_G4z_S_D2y_S_C2001003+ABZ*I_ERI_F3z_S_D2y_S_C2001003;
  abcd[170] = I_ERI_G4x_S_Dyz_S_C2001003+ABX*I_ERI_F3x_S_Dyz_S_C2001003;
  abcd[171] = I_ERI_G3xy_S_Dyz_S_C2001003+ABX*I_ERI_F2xy_S_Dyz_S_C2001003;
  abcd[172] = I_ERI_G3xz_S_Dyz_S_C2001003+ABX*I_ERI_F2xz_S_Dyz_S_C2001003;
  abcd[173] = I_ERI_G2x2y_S_Dyz_S_C2001003+ABX*I_ERI_Fx2y_S_Dyz_S_C2001003;
  abcd[174] = I_ERI_G2xyz_S_Dyz_S_C2001003+ABX*I_ERI_Fxyz_S_Dyz_S_C2001003;
  abcd[175] = I_ERI_G2x2z_S_Dyz_S_C2001003+ABX*I_ERI_Fx2z_S_Dyz_S_C2001003;
  abcd[176] = I_ERI_Gx3y_S_Dyz_S_C2001003+ABX*I_ERI_F3y_S_Dyz_S_C2001003;
  abcd[177] = I_ERI_Gx2yz_S_Dyz_S_C2001003+ABX*I_ERI_F2yz_S_Dyz_S_C2001003;
  abcd[178] = I_ERI_Gxy2z_S_Dyz_S_C2001003+ABX*I_ERI_Fy2z_S_Dyz_S_C2001003;
  abcd[179] = I_ERI_Gx3z_S_Dyz_S_C2001003+ABX*I_ERI_F3z_S_Dyz_S_C2001003;
  abcd[180] = I_ERI_G3xy_S_Dyz_S_C2001003+ABY*I_ERI_F3x_S_Dyz_S_C2001003;
  abcd[181] = I_ERI_G2x2y_S_Dyz_S_C2001003+ABY*I_ERI_F2xy_S_Dyz_S_C2001003;
  abcd[182] = I_ERI_G2xyz_S_Dyz_S_C2001003+ABY*I_ERI_F2xz_S_Dyz_S_C2001003;
  abcd[183] = I_ERI_Gx3y_S_Dyz_S_C2001003+ABY*I_ERI_Fx2y_S_Dyz_S_C2001003;
  abcd[184] = I_ERI_Gx2yz_S_Dyz_S_C2001003+ABY*I_ERI_Fxyz_S_Dyz_S_C2001003;
  abcd[185] = I_ERI_Gxy2z_S_Dyz_S_C2001003+ABY*I_ERI_Fx2z_S_Dyz_S_C2001003;
  abcd[186] = I_ERI_G4y_S_Dyz_S_C2001003+ABY*I_ERI_F3y_S_Dyz_S_C2001003;
  abcd[187] = I_ERI_G3yz_S_Dyz_S_C2001003+ABY*I_ERI_F2yz_S_Dyz_S_C2001003;
  abcd[188] = I_ERI_G2y2z_S_Dyz_S_C2001003+ABY*I_ERI_Fy2z_S_Dyz_S_C2001003;
  abcd[189] = I_ERI_Gy3z_S_Dyz_S_C2001003+ABY*I_ERI_F3z_S_Dyz_S_C2001003;
  abcd[190] = I_ERI_G3xz_S_Dyz_S_C2001003+ABZ*I_ERI_F3x_S_Dyz_S_C2001003;
  abcd[191] = I_ERI_G2xyz_S_Dyz_S_C2001003+ABZ*I_ERI_F2xy_S_Dyz_S_C2001003;
  abcd[192] = I_ERI_G2x2z_S_Dyz_S_C2001003+ABZ*I_ERI_F2xz_S_Dyz_S_C2001003;
  abcd[193] = I_ERI_Gx2yz_S_Dyz_S_C2001003+ABZ*I_ERI_Fx2y_S_Dyz_S_C2001003;
  abcd[194] = I_ERI_Gxy2z_S_Dyz_S_C2001003+ABZ*I_ERI_Fxyz_S_Dyz_S_C2001003;
  abcd[195] = I_ERI_Gx3z_S_Dyz_S_C2001003+ABZ*I_ERI_Fx2z_S_Dyz_S_C2001003;
  abcd[196] = I_ERI_G3yz_S_Dyz_S_C2001003+ABZ*I_ERI_F3y_S_Dyz_S_C2001003;
  abcd[197] = I_ERI_G2y2z_S_Dyz_S_C2001003+ABZ*I_ERI_F2yz_S_Dyz_S_C2001003;
  abcd[198] = I_ERI_Gy3z_S_Dyz_S_C2001003+ABZ*I_ERI_Fy2z_S_Dyz_S_C2001003;
  abcd[199] = I_ERI_G4z_S_Dyz_S_C2001003+ABZ*I_ERI_F3z_S_Dyz_S_C2001003;
  abcd[210] = I_ERI_G4x_S_D2z_S_C2001003+ABX*I_ERI_F3x_S_D2z_S_C2001003;
  abcd[211] = I_ERI_G3xy_S_D2z_S_C2001003+ABX*I_ERI_F2xy_S_D2z_S_C2001003;
  abcd[212] = I_ERI_G3xz_S_D2z_S_C2001003+ABX*I_ERI_F2xz_S_D2z_S_C2001003;
  abcd[213] = I_ERI_G2x2y_S_D2z_S_C2001003+ABX*I_ERI_Fx2y_S_D2z_S_C2001003;
  abcd[214] = I_ERI_G2xyz_S_D2z_S_C2001003+ABX*I_ERI_Fxyz_S_D2z_S_C2001003;
  abcd[215] = I_ERI_G2x2z_S_D2z_S_C2001003+ABX*I_ERI_Fx2z_S_D2z_S_C2001003;
  abcd[216] = I_ERI_Gx3y_S_D2z_S_C2001003+ABX*I_ERI_F3y_S_D2z_S_C2001003;
  abcd[217] = I_ERI_Gx2yz_S_D2z_S_C2001003+ABX*I_ERI_F2yz_S_D2z_S_C2001003;
  abcd[218] = I_ERI_Gxy2z_S_D2z_S_C2001003+ABX*I_ERI_Fy2z_S_D2z_S_C2001003;
  abcd[219] = I_ERI_Gx3z_S_D2z_S_C2001003+ABX*I_ERI_F3z_S_D2z_S_C2001003;
  abcd[220] = I_ERI_G3xy_S_D2z_S_C2001003+ABY*I_ERI_F3x_S_D2z_S_C2001003;
  abcd[221] = I_ERI_G2x2y_S_D2z_S_C2001003+ABY*I_ERI_F2xy_S_D2z_S_C2001003;
  abcd[222] = I_ERI_G2xyz_S_D2z_S_C2001003+ABY*I_ERI_F2xz_S_D2z_S_C2001003;
  abcd[223] = I_ERI_Gx3y_S_D2z_S_C2001003+ABY*I_ERI_Fx2y_S_D2z_S_C2001003;
  abcd[224] = I_ERI_Gx2yz_S_D2z_S_C2001003+ABY*I_ERI_Fxyz_S_D2z_S_C2001003;
  abcd[225] = I_ERI_Gxy2z_S_D2z_S_C2001003+ABY*I_ERI_Fx2z_S_D2z_S_C2001003;
  abcd[226] = I_ERI_G4y_S_D2z_S_C2001003+ABY*I_ERI_F3y_S_D2z_S_C2001003;
  abcd[227] = I_ERI_G3yz_S_D2z_S_C2001003+ABY*I_ERI_F2yz_S_D2z_S_C2001003;
  abcd[228] = I_ERI_G2y2z_S_D2z_S_C2001003+ABY*I_ERI_Fy2z_S_D2z_S_C2001003;
  abcd[229] = I_ERI_Gy3z_S_D2z_S_C2001003+ABY*I_ERI_F3z_S_D2z_S_C2001003;
  abcd[230] = I_ERI_G3xz_S_D2z_S_C2001003+ABZ*I_ERI_F3x_S_D2z_S_C2001003;
  abcd[231] = I_ERI_G2xyz_S_D2z_S_C2001003+ABZ*I_ERI_F2xy_S_D2z_S_C2001003;
  abcd[232] = I_ERI_G2x2z_S_D2z_S_C2001003+ABZ*I_ERI_F2xz_S_D2z_S_C2001003;
  abcd[233] = I_ERI_Gx2yz_S_D2z_S_C2001003+ABZ*I_ERI_Fx2y_S_D2z_S_C2001003;
  abcd[234] = I_ERI_Gxy2z_S_D2z_S_C2001003+ABZ*I_ERI_Fxyz_S_D2z_S_C2001003;
  abcd[235] = I_ERI_Gx3z_S_D2z_S_C2001003+ABZ*I_ERI_Fx2z_S_D2z_S_C2001003;
  abcd[236] = I_ERI_G3yz_S_D2z_S_C2001003+ABZ*I_ERI_F3y_S_D2z_S_C2001003;
  abcd[237] = I_ERI_G2y2z_S_D2z_S_C2001003+ABZ*I_ERI_F2yz_S_D2z_S_C2001003;
  abcd[238] = I_ERI_Gy3z_S_D2z_S_C2001003+ABZ*I_ERI_Fy2z_S_D2z_S_C2001003;
  abcd[239] = I_ERI_G4z_S_D2z_S_C2001003+ABZ*I_ERI_F3z_S_D2z_S_C2001003;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_D_P_C1002001003
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_P_C1002001003
   * RHS shell quartet name: SQ_ERI_F_S_D_P_C1002001003
   ************************************************************/
  abcd[250] = I_ERI_G4x_S_D2x_Px_C1002001003+ABX*I_ERI_F3x_S_D2x_Px_C1002001003;
  abcd[251] = I_ERI_G3xy_S_D2x_Px_C1002001003+ABX*I_ERI_F2xy_S_D2x_Px_C1002001003;
  abcd[252] = I_ERI_G3xz_S_D2x_Px_C1002001003+ABX*I_ERI_F2xz_S_D2x_Px_C1002001003;
  abcd[253] = I_ERI_G2x2y_S_D2x_Px_C1002001003+ABX*I_ERI_Fx2y_S_D2x_Px_C1002001003;
  abcd[254] = I_ERI_G2xyz_S_D2x_Px_C1002001003+ABX*I_ERI_Fxyz_S_D2x_Px_C1002001003;
  abcd[255] = I_ERI_G2x2z_S_D2x_Px_C1002001003+ABX*I_ERI_Fx2z_S_D2x_Px_C1002001003;
  abcd[256] = I_ERI_Gx3y_S_D2x_Px_C1002001003+ABX*I_ERI_F3y_S_D2x_Px_C1002001003;
  abcd[257] = I_ERI_Gx2yz_S_D2x_Px_C1002001003+ABX*I_ERI_F2yz_S_D2x_Px_C1002001003;
  abcd[258] = I_ERI_Gxy2z_S_D2x_Px_C1002001003+ABX*I_ERI_Fy2z_S_D2x_Px_C1002001003;
  abcd[259] = I_ERI_Gx3z_S_D2x_Px_C1002001003+ABX*I_ERI_F3z_S_D2x_Px_C1002001003;
  abcd[260] = I_ERI_G3xy_S_D2x_Px_C1002001003+ABY*I_ERI_F3x_S_D2x_Px_C1002001003;
  abcd[261] = I_ERI_G2x2y_S_D2x_Px_C1002001003+ABY*I_ERI_F2xy_S_D2x_Px_C1002001003;
  abcd[262] = I_ERI_G2xyz_S_D2x_Px_C1002001003+ABY*I_ERI_F2xz_S_D2x_Px_C1002001003;
  abcd[263] = I_ERI_Gx3y_S_D2x_Px_C1002001003+ABY*I_ERI_Fx2y_S_D2x_Px_C1002001003;
  abcd[264] = I_ERI_Gx2yz_S_D2x_Px_C1002001003+ABY*I_ERI_Fxyz_S_D2x_Px_C1002001003;
  abcd[265] = I_ERI_Gxy2z_S_D2x_Px_C1002001003+ABY*I_ERI_Fx2z_S_D2x_Px_C1002001003;
  abcd[266] = I_ERI_G4y_S_D2x_Px_C1002001003+ABY*I_ERI_F3y_S_D2x_Px_C1002001003;
  abcd[267] = I_ERI_G3yz_S_D2x_Px_C1002001003+ABY*I_ERI_F2yz_S_D2x_Px_C1002001003;
  abcd[268] = I_ERI_G2y2z_S_D2x_Px_C1002001003+ABY*I_ERI_Fy2z_S_D2x_Px_C1002001003;
  abcd[269] = I_ERI_Gy3z_S_D2x_Px_C1002001003+ABY*I_ERI_F3z_S_D2x_Px_C1002001003;
  abcd[270] = I_ERI_G3xz_S_D2x_Px_C1002001003+ABZ*I_ERI_F3x_S_D2x_Px_C1002001003;
  abcd[271] = I_ERI_G2xyz_S_D2x_Px_C1002001003+ABZ*I_ERI_F2xy_S_D2x_Px_C1002001003;
  abcd[272] = I_ERI_G2x2z_S_D2x_Px_C1002001003+ABZ*I_ERI_F2xz_S_D2x_Px_C1002001003;
  abcd[273] = I_ERI_Gx2yz_S_D2x_Px_C1002001003+ABZ*I_ERI_Fx2y_S_D2x_Px_C1002001003;
  abcd[274] = I_ERI_Gxy2z_S_D2x_Px_C1002001003+ABZ*I_ERI_Fxyz_S_D2x_Px_C1002001003;
  abcd[275] = I_ERI_Gx3z_S_D2x_Px_C1002001003+ABZ*I_ERI_Fx2z_S_D2x_Px_C1002001003;
  abcd[276] = I_ERI_G3yz_S_D2x_Px_C1002001003+ABZ*I_ERI_F3y_S_D2x_Px_C1002001003;
  abcd[277] = I_ERI_G2y2z_S_D2x_Px_C1002001003+ABZ*I_ERI_F2yz_S_D2x_Px_C1002001003;
  abcd[278] = I_ERI_Gy3z_S_D2x_Px_C1002001003+ABZ*I_ERI_Fy2z_S_D2x_Px_C1002001003;
  abcd[279] = I_ERI_G4z_S_D2x_Px_C1002001003+ABZ*I_ERI_F3z_S_D2x_Px_C1002001003;
  abcd[290] = I_ERI_G4x_S_Dxy_Px_C1002001003+ABX*I_ERI_F3x_S_Dxy_Px_C1002001003;
  abcd[291] = I_ERI_G3xy_S_Dxy_Px_C1002001003+ABX*I_ERI_F2xy_S_Dxy_Px_C1002001003;
  abcd[292] = I_ERI_G3xz_S_Dxy_Px_C1002001003+ABX*I_ERI_F2xz_S_Dxy_Px_C1002001003;
  abcd[293] = I_ERI_G2x2y_S_Dxy_Px_C1002001003+ABX*I_ERI_Fx2y_S_Dxy_Px_C1002001003;
  abcd[294] = I_ERI_G2xyz_S_Dxy_Px_C1002001003+ABX*I_ERI_Fxyz_S_Dxy_Px_C1002001003;
  abcd[295] = I_ERI_G2x2z_S_Dxy_Px_C1002001003+ABX*I_ERI_Fx2z_S_Dxy_Px_C1002001003;
  abcd[296] = I_ERI_Gx3y_S_Dxy_Px_C1002001003+ABX*I_ERI_F3y_S_Dxy_Px_C1002001003;
  abcd[297] = I_ERI_Gx2yz_S_Dxy_Px_C1002001003+ABX*I_ERI_F2yz_S_Dxy_Px_C1002001003;
  abcd[298] = I_ERI_Gxy2z_S_Dxy_Px_C1002001003+ABX*I_ERI_Fy2z_S_Dxy_Px_C1002001003;
  abcd[299] = I_ERI_Gx3z_S_Dxy_Px_C1002001003+ABX*I_ERI_F3z_S_Dxy_Px_C1002001003;
  abcd[300] = I_ERI_G3xy_S_Dxy_Px_C1002001003+ABY*I_ERI_F3x_S_Dxy_Px_C1002001003;
  abcd[301] = I_ERI_G2x2y_S_Dxy_Px_C1002001003+ABY*I_ERI_F2xy_S_Dxy_Px_C1002001003;
  abcd[302] = I_ERI_G2xyz_S_Dxy_Px_C1002001003+ABY*I_ERI_F2xz_S_Dxy_Px_C1002001003;
  abcd[303] = I_ERI_Gx3y_S_Dxy_Px_C1002001003+ABY*I_ERI_Fx2y_S_Dxy_Px_C1002001003;
  abcd[304] = I_ERI_Gx2yz_S_Dxy_Px_C1002001003+ABY*I_ERI_Fxyz_S_Dxy_Px_C1002001003;
  abcd[305] = I_ERI_Gxy2z_S_Dxy_Px_C1002001003+ABY*I_ERI_Fx2z_S_Dxy_Px_C1002001003;
  abcd[306] = I_ERI_G4y_S_Dxy_Px_C1002001003+ABY*I_ERI_F3y_S_Dxy_Px_C1002001003;
  abcd[307] = I_ERI_G3yz_S_Dxy_Px_C1002001003+ABY*I_ERI_F2yz_S_Dxy_Px_C1002001003;
  abcd[308] = I_ERI_G2y2z_S_Dxy_Px_C1002001003+ABY*I_ERI_Fy2z_S_Dxy_Px_C1002001003;
  abcd[309] = I_ERI_Gy3z_S_Dxy_Px_C1002001003+ABY*I_ERI_F3z_S_Dxy_Px_C1002001003;
  abcd[310] = I_ERI_G3xz_S_Dxy_Px_C1002001003+ABZ*I_ERI_F3x_S_Dxy_Px_C1002001003;
  abcd[311] = I_ERI_G2xyz_S_Dxy_Px_C1002001003+ABZ*I_ERI_F2xy_S_Dxy_Px_C1002001003;
  abcd[312] = I_ERI_G2x2z_S_Dxy_Px_C1002001003+ABZ*I_ERI_F2xz_S_Dxy_Px_C1002001003;
  abcd[313] = I_ERI_Gx2yz_S_Dxy_Px_C1002001003+ABZ*I_ERI_Fx2y_S_Dxy_Px_C1002001003;
  abcd[314] = I_ERI_Gxy2z_S_Dxy_Px_C1002001003+ABZ*I_ERI_Fxyz_S_Dxy_Px_C1002001003;
  abcd[315] = I_ERI_Gx3z_S_Dxy_Px_C1002001003+ABZ*I_ERI_Fx2z_S_Dxy_Px_C1002001003;
  abcd[316] = I_ERI_G3yz_S_Dxy_Px_C1002001003+ABZ*I_ERI_F3y_S_Dxy_Px_C1002001003;
  abcd[317] = I_ERI_G2y2z_S_Dxy_Px_C1002001003+ABZ*I_ERI_F2yz_S_Dxy_Px_C1002001003;
  abcd[318] = I_ERI_Gy3z_S_Dxy_Px_C1002001003+ABZ*I_ERI_Fy2z_S_Dxy_Px_C1002001003;
  abcd[319] = I_ERI_G4z_S_Dxy_Px_C1002001003+ABZ*I_ERI_F3z_S_Dxy_Px_C1002001003;
  abcd[330] = I_ERI_G4x_S_Dxz_Px_C1002001003+ABX*I_ERI_F3x_S_Dxz_Px_C1002001003;
  abcd[331] = I_ERI_G3xy_S_Dxz_Px_C1002001003+ABX*I_ERI_F2xy_S_Dxz_Px_C1002001003;
  abcd[332] = I_ERI_G3xz_S_Dxz_Px_C1002001003+ABX*I_ERI_F2xz_S_Dxz_Px_C1002001003;
  abcd[333] = I_ERI_G2x2y_S_Dxz_Px_C1002001003+ABX*I_ERI_Fx2y_S_Dxz_Px_C1002001003;
  abcd[334] = I_ERI_G2xyz_S_Dxz_Px_C1002001003+ABX*I_ERI_Fxyz_S_Dxz_Px_C1002001003;
  abcd[335] = I_ERI_G2x2z_S_Dxz_Px_C1002001003+ABX*I_ERI_Fx2z_S_Dxz_Px_C1002001003;
  abcd[336] = I_ERI_Gx3y_S_Dxz_Px_C1002001003+ABX*I_ERI_F3y_S_Dxz_Px_C1002001003;
  abcd[337] = I_ERI_Gx2yz_S_Dxz_Px_C1002001003+ABX*I_ERI_F2yz_S_Dxz_Px_C1002001003;
  abcd[338] = I_ERI_Gxy2z_S_Dxz_Px_C1002001003+ABX*I_ERI_Fy2z_S_Dxz_Px_C1002001003;
  abcd[339] = I_ERI_Gx3z_S_Dxz_Px_C1002001003+ABX*I_ERI_F3z_S_Dxz_Px_C1002001003;
  abcd[340] = I_ERI_G3xy_S_Dxz_Px_C1002001003+ABY*I_ERI_F3x_S_Dxz_Px_C1002001003;
  abcd[341] = I_ERI_G2x2y_S_Dxz_Px_C1002001003+ABY*I_ERI_F2xy_S_Dxz_Px_C1002001003;
  abcd[342] = I_ERI_G2xyz_S_Dxz_Px_C1002001003+ABY*I_ERI_F2xz_S_Dxz_Px_C1002001003;
  abcd[343] = I_ERI_Gx3y_S_Dxz_Px_C1002001003+ABY*I_ERI_Fx2y_S_Dxz_Px_C1002001003;
  abcd[344] = I_ERI_Gx2yz_S_Dxz_Px_C1002001003+ABY*I_ERI_Fxyz_S_Dxz_Px_C1002001003;
  abcd[345] = I_ERI_Gxy2z_S_Dxz_Px_C1002001003+ABY*I_ERI_Fx2z_S_Dxz_Px_C1002001003;
  abcd[346] = I_ERI_G4y_S_Dxz_Px_C1002001003+ABY*I_ERI_F3y_S_Dxz_Px_C1002001003;
  abcd[347] = I_ERI_G3yz_S_Dxz_Px_C1002001003+ABY*I_ERI_F2yz_S_Dxz_Px_C1002001003;
  abcd[348] = I_ERI_G2y2z_S_Dxz_Px_C1002001003+ABY*I_ERI_Fy2z_S_Dxz_Px_C1002001003;
  abcd[349] = I_ERI_Gy3z_S_Dxz_Px_C1002001003+ABY*I_ERI_F3z_S_Dxz_Px_C1002001003;
  abcd[350] = I_ERI_G3xz_S_Dxz_Px_C1002001003+ABZ*I_ERI_F3x_S_Dxz_Px_C1002001003;
  abcd[351] = I_ERI_G2xyz_S_Dxz_Px_C1002001003+ABZ*I_ERI_F2xy_S_Dxz_Px_C1002001003;
  abcd[352] = I_ERI_G2x2z_S_Dxz_Px_C1002001003+ABZ*I_ERI_F2xz_S_Dxz_Px_C1002001003;
  abcd[353] = I_ERI_Gx2yz_S_Dxz_Px_C1002001003+ABZ*I_ERI_Fx2y_S_Dxz_Px_C1002001003;
  abcd[354] = I_ERI_Gxy2z_S_Dxz_Px_C1002001003+ABZ*I_ERI_Fxyz_S_Dxz_Px_C1002001003;
  abcd[355] = I_ERI_Gx3z_S_Dxz_Px_C1002001003+ABZ*I_ERI_Fx2z_S_Dxz_Px_C1002001003;
  abcd[356] = I_ERI_G3yz_S_Dxz_Px_C1002001003+ABZ*I_ERI_F3y_S_Dxz_Px_C1002001003;
  abcd[357] = I_ERI_G2y2z_S_Dxz_Px_C1002001003+ABZ*I_ERI_F2yz_S_Dxz_Px_C1002001003;
  abcd[358] = I_ERI_Gy3z_S_Dxz_Px_C1002001003+ABZ*I_ERI_Fy2z_S_Dxz_Px_C1002001003;
  abcd[359] = I_ERI_G4z_S_Dxz_Px_C1002001003+ABZ*I_ERI_F3z_S_Dxz_Px_C1002001003;
  abcd[370] = I_ERI_G4x_S_D2y_Px_C1002001003+ABX*I_ERI_F3x_S_D2y_Px_C1002001003;
  abcd[371] = I_ERI_G3xy_S_D2y_Px_C1002001003+ABX*I_ERI_F2xy_S_D2y_Px_C1002001003;
  abcd[372] = I_ERI_G3xz_S_D2y_Px_C1002001003+ABX*I_ERI_F2xz_S_D2y_Px_C1002001003;
  abcd[373] = I_ERI_G2x2y_S_D2y_Px_C1002001003+ABX*I_ERI_Fx2y_S_D2y_Px_C1002001003;
  abcd[374] = I_ERI_G2xyz_S_D2y_Px_C1002001003+ABX*I_ERI_Fxyz_S_D2y_Px_C1002001003;
  abcd[375] = I_ERI_G2x2z_S_D2y_Px_C1002001003+ABX*I_ERI_Fx2z_S_D2y_Px_C1002001003;
  abcd[376] = I_ERI_Gx3y_S_D2y_Px_C1002001003+ABX*I_ERI_F3y_S_D2y_Px_C1002001003;
  abcd[377] = I_ERI_Gx2yz_S_D2y_Px_C1002001003+ABX*I_ERI_F2yz_S_D2y_Px_C1002001003;
  abcd[378] = I_ERI_Gxy2z_S_D2y_Px_C1002001003+ABX*I_ERI_Fy2z_S_D2y_Px_C1002001003;
  abcd[379] = I_ERI_Gx3z_S_D2y_Px_C1002001003+ABX*I_ERI_F3z_S_D2y_Px_C1002001003;
  abcd[380] = I_ERI_G3xy_S_D2y_Px_C1002001003+ABY*I_ERI_F3x_S_D2y_Px_C1002001003;
  abcd[381] = I_ERI_G2x2y_S_D2y_Px_C1002001003+ABY*I_ERI_F2xy_S_D2y_Px_C1002001003;
  abcd[382] = I_ERI_G2xyz_S_D2y_Px_C1002001003+ABY*I_ERI_F2xz_S_D2y_Px_C1002001003;
  abcd[383] = I_ERI_Gx3y_S_D2y_Px_C1002001003+ABY*I_ERI_Fx2y_S_D2y_Px_C1002001003;
  abcd[384] = I_ERI_Gx2yz_S_D2y_Px_C1002001003+ABY*I_ERI_Fxyz_S_D2y_Px_C1002001003;
  abcd[385] = I_ERI_Gxy2z_S_D2y_Px_C1002001003+ABY*I_ERI_Fx2z_S_D2y_Px_C1002001003;
  abcd[386] = I_ERI_G4y_S_D2y_Px_C1002001003+ABY*I_ERI_F3y_S_D2y_Px_C1002001003;
  abcd[387] = I_ERI_G3yz_S_D2y_Px_C1002001003+ABY*I_ERI_F2yz_S_D2y_Px_C1002001003;
  abcd[388] = I_ERI_G2y2z_S_D2y_Px_C1002001003+ABY*I_ERI_Fy2z_S_D2y_Px_C1002001003;
  abcd[389] = I_ERI_Gy3z_S_D2y_Px_C1002001003+ABY*I_ERI_F3z_S_D2y_Px_C1002001003;
  abcd[390] = I_ERI_G3xz_S_D2y_Px_C1002001003+ABZ*I_ERI_F3x_S_D2y_Px_C1002001003;
  abcd[391] = I_ERI_G2xyz_S_D2y_Px_C1002001003+ABZ*I_ERI_F2xy_S_D2y_Px_C1002001003;
  abcd[392] = I_ERI_G2x2z_S_D2y_Px_C1002001003+ABZ*I_ERI_F2xz_S_D2y_Px_C1002001003;
  abcd[393] = I_ERI_Gx2yz_S_D2y_Px_C1002001003+ABZ*I_ERI_Fx2y_S_D2y_Px_C1002001003;
  abcd[394] = I_ERI_Gxy2z_S_D2y_Px_C1002001003+ABZ*I_ERI_Fxyz_S_D2y_Px_C1002001003;
  abcd[395] = I_ERI_Gx3z_S_D2y_Px_C1002001003+ABZ*I_ERI_Fx2z_S_D2y_Px_C1002001003;
  abcd[396] = I_ERI_G3yz_S_D2y_Px_C1002001003+ABZ*I_ERI_F3y_S_D2y_Px_C1002001003;
  abcd[397] = I_ERI_G2y2z_S_D2y_Px_C1002001003+ABZ*I_ERI_F2yz_S_D2y_Px_C1002001003;
  abcd[398] = I_ERI_Gy3z_S_D2y_Px_C1002001003+ABZ*I_ERI_Fy2z_S_D2y_Px_C1002001003;
  abcd[399] = I_ERI_G4z_S_D2y_Px_C1002001003+ABZ*I_ERI_F3z_S_D2y_Px_C1002001003;
  abcd[410] = I_ERI_G4x_S_Dyz_Px_C1002001003+ABX*I_ERI_F3x_S_Dyz_Px_C1002001003;
  abcd[411] = I_ERI_G3xy_S_Dyz_Px_C1002001003+ABX*I_ERI_F2xy_S_Dyz_Px_C1002001003;
  abcd[412] = I_ERI_G3xz_S_Dyz_Px_C1002001003+ABX*I_ERI_F2xz_S_Dyz_Px_C1002001003;
  abcd[413] = I_ERI_G2x2y_S_Dyz_Px_C1002001003+ABX*I_ERI_Fx2y_S_Dyz_Px_C1002001003;
  abcd[414] = I_ERI_G2xyz_S_Dyz_Px_C1002001003+ABX*I_ERI_Fxyz_S_Dyz_Px_C1002001003;
  abcd[415] = I_ERI_G2x2z_S_Dyz_Px_C1002001003+ABX*I_ERI_Fx2z_S_Dyz_Px_C1002001003;
  abcd[416] = I_ERI_Gx3y_S_Dyz_Px_C1002001003+ABX*I_ERI_F3y_S_Dyz_Px_C1002001003;
  abcd[417] = I_ERI_Gx2yz_S_Dyz_Px_C1002001003+ABX*I_ERI_F2yz_S_Dyz_Px_C1002001003;
  abcd[418] = I_ERI_Gxy2z_S_Dyz_Px_C1002001003+ABX*I_ERI_Fy2z_S_Dyz_Px_C1002001003;
  abcd[419] = I_ERI_Gx3z_S_Dyz_Px_C1002001003+ABX*I_ERI_F3z_S_Dyz_Px_C1002001003;
  abcd[420] = I_ERI_G3xy_S_Dyz_Px_C1002001003+ABY*I_ERI_F3x_S_Dyz_Px_C1002001003;
  abcd[421] = I_ERI_G2x2y_S_Dyz_Px_C1002001003+ABY*I_ERI_F2xy_S_Dyz_Px_C1002001003;
  abcd[422] = I_ERI_G2xyz_S_Dyz_Px_C1002001003+ABY*I_ERI_F2xz_S_Dyz_Px_C1002001003;
  abcd[423] = I_ERI_Gx3y_S_Dyz_Px_C1002001003+ABY*I_ERI_Fx2y_S_Dyz_Px_C1002001003;
  abcd[424] = I_ERI_Gx2yz_S_Dyz_Px_C1002001003+ABY*I_ERI_Fxyz_S_Dyz_Px_C1002001003;
  abcd[425] = I_ERI_Gxy2z_S_Dyz_Px_C1002001003+ABY*I_ERI_Fx2z_S_Dyz_Px_C1002001003;
  abcd[426] = I_ERI_G4y_S_Dyz_Px_C1002001003+ABY*I_ERI_F3y_S_Dyz_Px_C1002001003;
  abcd[427] = I_ERI_G3yz_S_Dyz_Px_C1002001003+ABY*I_ERI_F2yz_S_Dyz_Px_C1002001003;
  abcd[428] = I_ERI_G2y2z_S_Dyz_Px_C1002001003+ABY*I_ERI_Fy2z_S_Dyz_Px_C1002001003;
  abcd[429] = I_ERI_Gy3z_S_Dyz_Px_C1002001003+ABY*I_ERI_F3z_S_Dyz_Px_C1002001003;
  abcd[430] = I_ERI_G3xz_S_Dyz_Px_C1002001003+ABZ*I_ERI_F3x_S_Dyz_Px_C1002001003;
  abcd[431] = I_ERI_G2xyz_S_Dyz_Px_C1002001003+ABZ*I_ERI_F2xy_S_Dyz_Px_C1002001003;
  abcd[432] = I_ERI_G2x2z_S_Dyz_Px_C1002001003+ABZ*I_ERI_F2xz_S_Dyz_Px_C1002001003;
  abcd[433] = I_ERI_Gx2yz_S_Dyz_Px_C1002001003+ABZ*I_ERI_Fx2y_S_Dyz_Px_C1002001003;
  abcd[434] = I_ERI_Gxy2z_S_Dyz_Px_C1002001003+ABZ*I_ERI_Fxyz_S_Dyz_Px_C1002001003;
  abcd[435] = I_ERI_Gx3z_S_Dyz_Px_C1002001003+ABZ*I_ERI_Fx2z_S_Dyz_Px_C1002001003;
  abcd[436] = I_ERI_G3yz_S_Dyz_Px_C1002001003+ABZ*I_ERI_F3y_S_Dyz_Px_C1002001003;
  abcd[437] = I_ERI_G2y2z_S_Dyz_Px_C1002001003+ABZ*I_ERI_F2yz_S_Dyz_Px_C1002001003;
  abcd[438] = I_ERI_Gy3z_S_Dyz_Px_C1002001003+ABZ*I_ERI_Fy2z_S_Dyz_Px_C1002001003;
  abcd[439] = I_ERI_G4z_S_Dyz_Px_C1002001003+ABZ*I_ERI_F3z_S_Dyz_Px_C1002001003;
  abcd[450] = I_ERI_G4x_S_D2z_Px_C1002001003+ABX*I_ERI_F3x_S_D2z_Px_C1002001003;
  abcd[451] = I_ERI_G3xy_S_D2z_Px_C1002001003+ABX*I_ERI_F2xy_S_D2z_Px_C1002001003;
  abcd[452] = I_ERI_G3xz_S_D2z_Px_C1002001003+ABX*I_ERI_F2xz_S_D2z_Px_C1002001003;
  abcd[453] = I_ERI_G2x2y_S_D2z_Px_C1002001003+ABX*I_ERI_Fx2y_S_D2z_Px_C1002001003;
  abcd[454] = I_ERI_G2xyz_S_D2z_Px_C1002001003+ABX*I_ERI_Fxyz_S_D2z_Px_C1002001003;
  abcd[455] = I_ERI_G2x2z_S_D2z_Px_C1002001003+ABX*I_ERI_Fx2z_S_D2z_Px_C1002001003;
  abcd[456] = I_ERI_Gx3y_S_D2z_Px_C1002001003+ABX*I_ERI_F3y_S_D2z_Px_C1002001003;
  abcd[457] = I_ERI_Gx2yz_S_D2z_Px_C1002001003+ABX*I_ERI_F2yz_S_D2z_Px_C1002001003;
  abcd[458] = I_ERI_Gxy2z_S_D2z_Px_C1002001003+ABX*I_ERI_Fy2z_S_D2z_Px_C1002001003;
  abcd[459] = I_ERI_Gx3z_S_D2z_Px_C1002001003+ABX*I_ERI_F3z_S_D2z_Px_C1002001003;
  abcd[460] = I_ERI_G3xy_S_D2z_Px_C1002001003+ABY*I_ERI_F3x_S_D2z_Px_C1002001003;
  abcd[461] = I_ERI_G2x2y_S_D2z_Px_C1002001003+ABY*I_ERI_F2xy_S_D2z_Px_C1002001003;
  abcd[462] = I_ERI_G2xyz_S_D2z_Px_C1002001003+ABY*I_ERI_F2xz_S_D2z_Px_C1002001003;
  abcd[463] = I_ERI_Gx3y_S_D2z_Px_C1002001003+ABY*I_ERI_Fx2y_S_D2z_Px_C1002001003;
  abcd[464] = I_ERI_Gx2yz_S_D2z_Px_C1002001003+ABY*I_ERI_Fxyz_S_D2z_Px_C1002001003;
  abcd[465] = I_ERI_Gxy2z_S_D2z_Px_C1002001003+ABY*I_ERI_Fx2z_S_D2z_Px_C1002001003;
  abcd[466] = I_ERI_G4y_S_D2z_Px_C1002001003+ABY*I_ERI_F3y_S_D2z_Px_C1002001003;
  abcd[467] = I_ERI_G3yz_S_D2z_Px_C1002001003+ABY*I_ERI_F2yz_S_D2z_Px_C1002001003;
  abcd[468] = I_ERI_G2y2z_S_D2z_Px_C1002001003+ABY*I_ERI_Fy2z_S_D2z_Px_C1002001003;
  abcd[469] = I_ERI_Gy3z_S_D2z_Px_C1002001003+ABY*I_ERI_F3z_S_D2z_Px_C1002001003;
  abcd[470] = I_ERI_G3xz_S_D2z_Px_C1002001003+ABZ*I_ERI_F3x_S_D2z_Px_C1002001003;
  abcd[471] = I_ERI_G2xyz_S_D2z_Px_C1002001003+ABZ*I_ERI_F2xy_S_D2z_Px_C1002001003;
  abcd[472] = I_ERI_G2x2z_S_D2z_Px_C1002001003+ABZ*I_ERI_F2xz_S_D2z_Px_C1002001003;
  abcd[473] = I_ERI_Gx2yz_S_D2z_Px_C1002001003+ABZ*I_ERI_Fx2y_S_D2z_Px_C1002001003;
  abcd[474] = I_ERI_Gxy2z_S_D2z_Px_C1002001003+ABZ*I_ERI_Fxyz_S_D2z_Px_C1002001003;
  abcd[475] = I_ERI_Gx3z_S_D2z_Px_C1002001003+ABZ*I_ERI_Fx2z_S_D2z_Px_C1002001003;
  abcd[476] = I_ERI_G3yz_S_D2z_Px_C1002001003+ABZ*I_ERI_F3y_S_D2z_Px_C1002001003;
  abcd[477] = I_ERI_G2y2z_S_D2z_Px_C1002001003+ABZ*I_ERI_F2yz_S_D2z_Px_C1002001003;
  abcd[478] = I_ERI_Gy3z_S_D2z_Px_C1002001003+ABZ*I_ERI_Fy2z_S_D2z_Px_C1002001003;
  abcd[479] = I_ERI_G4z_S_D2z_Px_C1002001003+ABZ*I_ERI_F3z_S_D2z_Px_C1002001003;
  abcd[490] = I_ERI_G4x_S_D2x_Py_C1002001003+ABX*I_ERI_F3x_S_D2x_Py_C1002001003;
  abcd[491] = I_ERI_G3xy_S_D2x_Py_C1002001003+ABX*I_ERI_F2xy_S_D2x_Py_C1002001003;
  abcd[492] = I_ERI_G3xz_S_D2x_Py_C1002001003+ABX*I_ERI_F2xz_S_D2x_Py_C1002001003;
  abcd[493] = I_ERI_G2x2y_S_D2x_Py_C1002001003+ABX*I_ERI_Fx2y_S_D2x_Py_C1002001003;
  abcd[494] = I_ERI_G2xyz_S_D2x_Py_C1002001003+ABX*I_ERI_Fxyz_S_D2x_Py_C1002001003;
  abcd[495] = I_ERI_G2x2z_S_D2x_Py_C1002001003+ABX*I_ERI_Fx2z_S_D2x_Py_C1002001003;
  abcd[496] = I_ERI_Gx3y_S_D2x_Py_C1002001003+ABX*I_ERI_F3y_S_D2x_Py_C1002001003;
  abcd[497] = I_ERI_Gx2yz_S_D2x_Py_C1002001003+ABX*I_ERI_F2yz_S_D2x_Py_C1002001003;
  abcd[498] = I_ERI_Gxy2z_S_D2x_Py_C1002001003+ABX*I_ERI_Fy2z_S_D2x_Py_C1002001003;
  abcd[499] = I_ERI_Gx3z_S_D2x_Py_C1002001003+ABX*I_ERI_F3z_S_D2x_Py_C1002001003;
  abcd[500] = I_ERI_G3xy_S_D2x_Py_C1002001003+ABY*I_ERI_F3x_S_D2x_Py_C1002001003;
  abcd[501] = I_ERI_G2x2y_S_D2x_Py_C1002001003+ABY*I_ERI_F2xy_S_D2x_Py_C1002001003;
  abcd[502] = I_ERI_G2xyz_S_D2x_Py_C1002001003+ABY*I_ERI_F2xz_S_D2x_Py_C1002001003;
  abcd[503] = I_ERI_Gx3y_S_D2x_Py_C1002001003+ABY*I_ERI_Fx2y_S_D2x_Py_C1002001003;
  abcd[504] = I_ERI_Gx2yz_S_D2x_Py_C1002001003+ABY*I_ERI_Fxyz_S_D2x_Py_C1002001003;
  abcd[505] = I_ERI_Gxy2z_S_D2x_Py_C1002001003+ABY*I_ERI_Fx2z_S_D2x_Py_C1002001003;
  abcd[506] = I_ERI_G4y_S_D2x_Py_C1002001003+ABY*I_ERI_F3y_S_D2x_Py_C1002001003;
  abcd[507] = I_ERI_G3yz_S_D2x_Py_C1002001003+ABY*I_ERI_F2yz_S_D2x_Py_C1002001003;
  abcd[508] = I_ERI_G2y2z_S_D2x_Py_C1002001003+ABY*I_ERI_Fy2z_S_D2x_Py_C1002001003;
  abcd[509] = I_ERI_Gy3z_S_D2x_Py_C1002001003+ABY*I_ERI_F3z_S_D2x_Py_C1002001003;
  abcd[510] = I_ERI_G3xz_S_D2x_Py_C1002001003+ABZ*I_ERI_F3x_S_D2x_Py_C1002001003;
  abcd[511] = I_ERI_G2xyz_S_D2x_Py_C1002001003+ABZ*I_ERI_F2xy_S_D2x_Py_C1002001003;
  abcd[512] = I_ERI_G2x2z_S_D2x_Py_C1002001003+ABZ*I_ERI_F2xz_S_D2x_Py_C1002001003;
  abcd[513] = I_ERI_Gx2yz_S_D2x_Py_C1002001003+ABZ*I_ERI_Fx2y_S_D2x_Py_C1002001003;
  abcd[514] = I_ERI_Gxy2z_S_D2x_Py_C1002001003+ABZ*I_ERI_Fxyz_S_D2x_Py_C1002001003;
  abcd[515] = I_ERI_Gx3z_S_D2x_Py_C1002001003+ABZ*I_ERI_Fx2z_S_D2x_Py_C1002001003;
  abcd[516] = I_ERI_G3yz_S_D2x_Py_C1002001003+ABZ*I_ERI_F3y_S_D2x_Py_C1002001003;
  abcd[517] = I_ERI_G2y2z_S_D2x_Py_C1002001003+ABZ*I_ERI_F2yz_S_D2x_Py_C1002001003;
  abcd[518] = I_ERI_Gy3z_S_D2x_Py_C1002001003+ABZ*I_ERI_Fy2z_S_D2x_Py_C1002001003;
  abcd[519] = I_ERI_G4z_S_D2x_Py_C1002001003+ABZ*I_ERI_F3z_S_D2x_Py_C1002001003;
  abcd[530] = I_ERI_G4x_S_Dxy_Py_C1002001003+ABX*I_ERI_F3x_S_Dxy_Py_C1002001003;
  abcd[531] = I_ERI_G3xy_S_Dxy_Py_C1002001003+ABX*I_ERI_F2xy_S_Dxy_Py_C1002001003;
  abcd[532] = I_ERI_G3xz_S_Dxy_Py_C1002001003+ABX*I_ERI_F2xz_S_Dxy_Py_C1002001003;
  abcd[533] = I_ERI_G2x2y_S_Dxy_Py_C1002001003+ABX*I_ERI_Fx2y_S_Dxy_Py_C1002001003;
  abcd[534] = I_ERI_G2xyz_S_Dxy_Py_C1002001003+ABX*I_ERI_Fxyz_S_Dxy_Py_C1002001003;
  abcd[535] = I_ERI_G2x2z_S_Dxy_Py_C1002001003+ABX*I_ERI_Fx2z_S_Dxy_Py_C1002001003;
  abcd[536] = I_ERI_Gx3y_S_Dxy_Py_C1002001003+ABX*I_ERI_F3y_S_Dxy_Py_C1002001003;
  abcd[537] = I_ERI_Gx2yz_S_Dxy_Py_C1002001003+ABX*I_ERI_F2yz_S_Dxy_Py_C1002001003;
  abcd[538] = I_ERI_Gxy2z_S_Dxy_Py_C1002001003+ABX*I_ERI_Fy2z_S_Dxy_Py_C1002001003;
  abcd[539] = I_ERI_Gx3z_S_Dxy_Py_C1002001003+ABX*I_ERI_F3z_S_Dxy_Py_C1002001003;
  abcd[540] = I_ERI_G3xy_S_Dxy_Py_C1002001003+ABY*I_ERI_F3x_S_Dxy_Py_C1002001003;
  abcd[541] = I_ERI_G2x2y_S_Dxy_Py_C1002001003+ABY*I_ERI_F2xy_S_Dxy_Py_C1002001003;
  abcd[542] = I_ERI_G2xyz_S_Dxy_Py_C1002001003+ABY*I_ERI_F2xz_S_Dxy_Py_C1002001003;
  abcd[543] = I_ERI_Gx3y_S_Dxy_Py_C1002001003+ABY*I_ERI_Fx2y_S_Dxy_Py_C1002001003;
  abcd[544] = I_ERI_Gx2yz_S_Dxy_Py_C1002001003+ABY*I_ERI_Fxyz_S_Dxy_Py_C1002001003;
  abcd[545] = I_ERI_Gxy2z_S_Dxy_Py_C1002001003+ABY*I_ERI_Fx2z_S_Dxy_Py_C1002001003;
  abcd[546] = I_ERI_G4y_S_Dxy_Py_C1002001003+ABY*I_ERI_F3y_S_Dxy_Py_C1002001003;
  abcd[547] = I_ERI_G3yz_S_Dxy_Py_C1002001003+ABY*I_ERI_F2yz_S_Dxy_Py_C1002001003;
  abcd[548] = I_ERI_G2y2z_S_Dxy_Py_C1002001003+ABY*I_ERI_Fy2z_S_Dxy_Py_C1002001003;
  abcd[549] = I_ERI_Gy3z_S_Dxy_Py_C1002001003+ABY*I_ERI_F3z_S_Dxy_Py_C1002001003;
  abcd[550] = I_ERI_G3xz_S_Dxy_Py_C1002001003+ABZ*I_ERI_F3x_S_Dxy_Py_C1002001003;
  abcd[551] = I_ERI_G2xyz_S_Dxy_Py_C1002001003+ABZ*I_ERI_F2xy_S_Dxy_Py_C1002001003;
  abcd[552] = I_ERI_G2x2z_S_Dxy_Py_C1002001003+ABZ*I_ERI_F2xz_S_Dxy_Py_C1002001003;
  abcd[553] = I_ERI_Gx2yz_S_Dxy_Py_C1002001003+ABZ*I_ERI_Fx2y_S_Dxy_Py_C1002001003;
  abcd[554] = I_ERI_Gxy2z_S_Dxy_Py_C1002001003+ABZ*I_ERI_Fxyz_S_Dxy_Py_C1002001003;
  abcd[555] = I_ERI_Gx3z_S_Dxy_Py_C1002001003+ABZ*I_ERI_Fx2z_S_Dxy_Py_C1002001003;
  abcd[556] = I_ERI_G3yz_S_Dxy_Py_C1002001003+ABZ*I_ERI_F3y_S_Dxy_Py_C1002001003;
  abcd[557] = I_ERI_G2y2z_S_Dxy_Py_C1002001003+ABZ*I_ERI_F2yz_S_Dxy_Py_C1002001003;
  abcd[558] = I_ERI_Gy3z_S_Dxy_Py_C1002001003+ABZ*I_ERI_Fy2z_S_Dxy_Py_C1002001003;
  abcd[559] = I_ERI_G4z_S_Dxy_Py_C1002001003+ABZ*I_ERI_F3z_S_Dxy_Py_C1002001003;
  abcd[570] = I_ERI_G4x_S_Dxz_Py_C1002001003+ABX*I_ERI_F3x_S_Dxz_Py_C1002001003;
  abcd[571] = I_ERI_G3xy_S_Dxz_Py_C1002001003+ABX*I_ERI_F2xy_S_Dxz_Py_C1002001003;
  abcd[572] = I_ERI_G3xz_S_Dxz_Py_C1002001003+ABX*I_ERI_F2xz_S_Dxz_Py_C1002001003;
  abcd[573] = I_ERI_G2x2y_S_Dxz_Py_C1002001003+ABX*I_ERI_Fx2y_S_Dxz_Py_C1002001003;
  abcd[574] = I_ERI_G2xyz_S_Dxz_Py_C1002001003+ABX*I_ERI_Fxyz_S_Dxz_Py_C1002001003;
  abcd[575] = I_ERI_G2x2z_S_Dxz_Py_C1002001003+ABX*I_ERI_Fx2z_S_Dxz_Py_C1002001003;
  abcd[576] = I_ERI_Gx3y_S_Dxz_Py_C1002001003+ABX*I_ERI_F3y_S_Dxz_Py_C1002001003;
  abcd[577] = I_ERI_Gx2yz_S_Dxz_Py_C1002001003+ABX*I_ERI_F2yz_S_Dxz_Py_C1002001003;
  abcd[578] = I_ERI_Gxy2z_S_Dxz_Py_C1002001003+ABX*I_ERI_Fy2z_S_Dxz_Py_C1002001003;
  abcd[579] = I_ERI_Gx3z_S_Dxz_Py_C1002001003+ABX*I_ERI_F3z_S_Dxz_Py_C1002001003;
  abcd[580] = I_ERI_G3xy_S_Dxz_Py_C1002001003+ABY*I_ERI_F3x_S_Dxz_Py_C1002001003;
  abcd[581] = I_ERI_G2x2y_S_Dxz_Py_C1002001003+ABY*I_ERI_F2xy_S_Dxz_Py_C1002001003;
  abcd[582] = I_ERI_G2xyz_S_Dxz_Py_C1002001003+ABY*I_ERI_F2xz_S_Dxz_Py_C1002001003;
  abcd[583] = I_ERI_Gx3y_S_Dxz_Py_C1002001003+ABY*I_ERI_Fx2y_S_Dxz_Py_C1002001003;
  abcd[584] = I_ERI_Gx2yz_S_Dxz_Py_C1002001003+ABY*I_ERI_Fxyz_S_Dxz_Py_C1002001003;
  abcd[585] = I_ERI_Gxy2z_S_Dxz_Py_C1002001003+ABY*I_ERI_Fx2z_S_Dxz_Py_C1002001003;
  abcd[586] = I_ERI_G4y_S_Dxz_Py_C1002001003+ABY*I_ERI_F3y_S_Dxz_Py_C1002001003;
  abcd[587] = I_ERI_G3yz_S_Dxz_Py_C1002001003+ABY*I_ERI_F2yz_S_Dxz_Py_C1002001003;
  abcd[588] = I_ERI_G2y2z_S_Dxz_Py_C1002001003+ABY*I_ERI_Fy2z_S_Dxz_Py_C1002001003;
  abcd[589] = I_ERI_Gy3z_S_Dxz_Py_C1002001003+ABY*I_ERI_F3z_S_Dxz_Py_C1002001003;
  abcd[590] = I_ERI_G3xz_S_Dxz_Py_C1002001003+ABZ*I_ERI_F3x_S_Dxz_Py_C1002001003;
  abcd[591] = I_ERI_G2xyz_S_Dxz_Py_C1002001003+ABZ*I_ERI_F2xy_S_Dxz_Py_C1002001003;
  abcd[592] = I_ERI_G2x2z_S_Dxz_Py_C1002001003+ABZ*I_ERI_F2xz_S_Dxz_Py_C1002001003;
  abcd[593] = I_ERI_Gx2yz_S_Dxz_Py_C1002001003+ABZ*I_ERI_Fx2y_S_Dxz_Py_C1002001003;
  abcd[594] = I_ERI_Gxy2z_S_Dxz_Py_C1002001003+ABZ*I_ERI_Fxyz_S_Dxz_Py_C1002001003;
  abcd[595] = I_ERI_Gx3z_S_Dxz_Py_C1002001003+ABZ*I_ERI_Fx2z_S_Dxz_Py_C1002001003;
  abcd[596] = I_ERI_G3yz_S_Dxz_Py_C1002001003+ABZ*I_ERI_F3y_S_Dxz_Py_C1002001003;
  abcd[597] = I_ERI_G2y2z_S_Dxz_Py_C1002001003+ABZ*I_ERI_F2yz_S_Dxz_Py_C1002001003;
  abcd[598] = I_ERI_Gy3z_S_Dxz_Py_C1002001003+ABZ*I_ERI_Fy2z_S_Dxz_Py_C1002001003;
  abcd[599] = I_ERI_G4z_S_Dxz_Py_C1002001003+ABZ*I_ERI_F3z_S_Dxz_Py_C1002001003;
  abcd[610] = I_ERI_G4x_S_D2y_Py_C1002001003+ABX*I_ERI_F3x_S_D2y_Py_C1002001003;
  abcd[611] = I_ERI_G3xy_S_D2y_Py_C1002001003+ABX*I_ERI_F2xy_S_D2y_Py_C1002001003;
  abcd[612] = I_ERI_G3xz_S_D2y_Py_C1002001003+ABX*I_ERI_F2xz_S_D2y_Py_C1002001003;
  abcd[613] = I_ERI_G2x2y_S_D2y_Py_C1002001003+ABX*I_ERI_Fx2y_S_D2y_Py_C1002001003;
  abcd[614] = I_ERI_G2xyz_S_D2y_Py_C1002001003+ABX*I_ERI_Fxyz_S_D2y_Py_C1002001003;
  abcd[615] = I_ERI_G2x2z_S_D2y_Py_C1002001003+ABX*I_ERI_Fx2z_S_D2y_Py_C1002001003;
  abcd[616] = I_ERI_Gx3y_S_D2y_Py_C1002001003+ABX*I_ERI_F3y_S_D2y_Py_C1002001003;
  abcd[617] = I_ERI_Gx2yz_S_D2y_Py_C1002001003+ABX*I_ERI_F2yz_S_D2y_Py_C1002001003;
  abcd[618] = I_ERI_Gxy2z_S_D2y_Py_C1002001003+ABX*I_ERI_Fy2z_S_D2y_Py_C1002001003;
  abcd[619] = I_ERI_Gx3z_S_D2y_Py_C1002001003+ABX*I_ERI_F3z_S_D2y_Py_C1002001003;
  abcd[620] = I_ERI_G3xy_S_D2y_Py_C1002001003+ABY*I_ERI_F3x_S_D2y_Py_C1002001003;
  abcd[621] = I_ERI_G2x2y_S_D2y_Py_C1002001003+ABY*I_ERI_F2xy_S_D2y_Py_C1002001003;
  abcd[622] = I_ERI_G2xyz_S_D2y_Py_C1002001003+ABY*I_ERI_F2xz_S_D2y_Py_C1002001003;
  abcd[623] = I_ERI_Gx3y_S_D2y_Py_C1002001003+ABY*I_ERI_Fx2y_S_D2y_Py_C1002001003;
  abcd[624] = I_ERI_Gx2yz_S_D2y_Py_C1002001003+ABY*I_ERI_Fxyz_S_D2y_Py_C1002001003;
  abcd[625] = I_ERI_Gxy2z_S_D2y_Py_C1002001003+ABY*I_ERI_Fx2z_S_D2y_Py_C1002001003;
  abcd[626] = I_ERI_G4y_S_D2y_Py_C1002001003+ABY*I_ERI_F3y_S_D2y_Py_C1002001003;
  abcd[627] = I_ERI_G3yz_S_D2y_Py_C1002001003+ABY*I_ERI_F2yz_S_D2y_Py_C1002001003;
  abcd[628] = I_ERI_G2y2z_S_D2y_Py_C1002001003+ABY*I_ERI_Fy2z_S_D2y_Py_C1002001003;
  abcd[629] = I_ERI_Gy3z_S_D2y_Py_C1002001003+ABY*I_ERI_F3z_S_D2y_Py_C1002001003;
  abcd[630] = I_ERI_G3xz_S_D2y_Py_C1002001003+ABZ*I_ERI_F3x_S_D2y_Py_C1002001003;
  abcd[631] = I_ERI_G2xyz_S_D2y_Py_C1002001003+ABZ*I_ERI_F2xy_S_D2y_Py_C1002001003;
  abcd[632] = I_ERI_G2x2z_S_D2y_Py_C1002001003+ABZ*I_ERI_F2xz_S_D2y_Py_C1002001003;
  abcd[633] = I_ERI_Gx2yz_S_D2y_Py_C1002001003+ABZ*I_ERI_Fx2y_S_D2y_Py_C1002001003;
  abcd[634] = I_ERI_Gxy2z_S_D2y_Py_C1002001003+ABZ*I_ERI_Fxyz_S_D2y_Py_C1002001003;
  abcd[635] = I_ERI_Gx3z_S_D2y_Py_C1002001003+ABZ*I_ERI_Fx2z_S_D2y_Py_C1002001003;
  abcd[636] = I_ERI_G3yz_S_D2y_Py_C1002001003+ABZ*I_ERI_F3y_S_D2y_Py_C1002001003;
  abcd[637] = I_ERI_G2y2z_S_D2y_Py_C1002001003+ABZ*I_ERI_F2yz_S_D2y_Py_C1002001003;
  abcd[638] = I_ERI_Gy3z_S_D2y_Py_C1002001003+ABZ*I_ERI_Fy2z_S_D2y_Py_C1002001003;
  abcd[639] = I_ERI_G4z_S_D2y_Py_C1002001003+ABZ*I_ERI_F3z_S_D2y_Py_C1002001003;
  abcd[650] = I_ERI_G4x_S_Dyz_Py_C1002001003+ABX*I_ERI_F3x_S_Dyz_Py_C1002001003;
  abcd[651] = I_ERI_G3xy_S_Dyz_Py_C1002001003+ABX*I_ERI_F2xy_S_Dyz_Py_C1002001003;
  abcd[652] = I_ERI_G3xz_S_Dyz_Py_C1002001003+ABX*I_ERI_F2xz_S_Dyz_Py_C1002001003;
  abcd[653] = I_ERI_G2x2y_S_Dyz_Py_C1002001003+ABX*I_ERI_Fx2y_S_Dyz_Py_C1002001003;
  abcd[654] = I_ERI_G2xyz_S_Dyz_Py_C1002001003+ABX*I_ERI_Fxyz_S_Dyz_Py_C1002001003;
  abcd[655] = I_ERI_G2x2z_S_Dyz_Py_C1002001003+ABX*I_ERI_Fx2z_S_Dyz_Py_C1002001003;
  abcd[656] = I_ERI_Gx3y_S_Dyz_Py_C1002001003+ABX*I_ERI_F3y_S_Dyz_Py_C1002001003;
  abcd[657] = I_ERI_Gx2yz_S_Dyz_Py_C1002001003+ABX*I_ERI_F2yz_S_Dyz_Py_C1002001003;
  abcd[658] = I_ERI_Gxy2z_S_Dyz_Py_C1002001003+ABX*I_ERI_Fy2z_S_Dyz_Py_C1002001003;
  abcd[659] = I_ERI_Gx3z_S_Dyz_Py_C1002001003+ABX*I_ERI_F3z_S_Dyz_Py_C1002001003;
  abcd[660] = I_ERI_G3xy_S_Dyz_Py_C1002001003+ABY*I_ERI_F3x_S_Dyz_Py_C1002001003;
  abcd[661] = I_ERI_G2x2y_S_Dyz_Py_C1002001003+ABY*I_ERI_F2xy_S_Dyz_Py_C1002001003;
  abcd[662] = I_ERI_G2xyz_S_Dyz_Py_C1002001003+ABY*I_ERI_F2xz_S_Dyz_Py_C1002001003;
  abcd[663] = I_ERI_Gx3y_S_Dyz_Py_C1002001003+ABY*I_ERI_Fx2y_S_Dyz_Py_C1002001003;
  abcd[664] = I_ERI_Gx2yz_S_Dyz_Py_C1002001003+ABY*I_ERI_Fxyz_S_Dyz_Py_C1002001003;
  abcd[665] = I_ERI_Gxy2z_S_Dyz_Py_C1002001003+ABY*I_ERI_Fx2z_S_Dyz_Py_C1002001003;
  abcd[666] = I_ERI_G4y_S_Dyz_Py_C1002001003+ABY*I_ERI_F3y_S_Dyz_Py_C1002001003;
  abcd[667] = I_ERI_G3yz_S_Dyz_Py_C1002001003+ABY*I_ERI_F2yz_S_Dyz_Py_C1002001003;
  abcd[668] = I_ERI_G2y2z_S_Dyz_Py_C1002001003+ABY*I_ERI_Fy2z_S_Dyz_Py_C1002001003;
  abcd[669] = I_ERI_Gy3z_S_Dyz_Py_C1002001003+ABY*I_ERI_F3z_S_Dyz_Py_C1002001003;
  abcd[670] = I_ERI_G3xz_S_Dyz_Py_C1002001003+ABZ*I_ERI_F3x_S_Dyz_Py_C1002001003;
  abcd[671] = I_ERI_G2xyz_S_Dyz_Py_C1002001003+ABZ*I_ERI_F2xy_S_Dyz_Py_C1002001003;
  abcd[672] = I_ERI_G2x2z_S_Dyz_Py_C1002001003+ABZ*I_ERI_F2xz_S_Dyz_Py_C1002001003;
  abcd[673] = I_ERI_Gx2yz_S_Dyz_Py_C1002001003+ABZ*I_ERI_Fx2y_S_Dyz_Py_C1002001003;
  abcd[674] = I_ERI_Gxy2z_S_Dyz_Py_C1002001003+ABZ*I_ERI_Fxyz_S_Dyz_Py_C1002001003;
  abcd[675] = I_ERI_Gx3z_S_Dyz_Py_C1002001003+ABZ*I_ERI_Fx2z_S_Dyz_Py_C1002001003;
  abcd[676] = I_ERI_G3yz_S_Dyz_Py_C1002001003+ABZ*I_ERI_F3y_S_Dyz_Py_C1002001003;
  abcd[677] = I_ERI_G2y2z_S_Dyz_Py_C1002001003+ABZ*I_ERI_F2yz_S_Dyz_Py_C1002001003;
  abcd[678] = I_ERI_Gy3z_S_Dyz_Py_C1002001003+ABZ*I_ERI_Fy2z_S_Dyz_Py_C1002001003;
  abcd[679] = I_ERI_G4z_S_Dyz_Py_C1002001003+ABZ*I_ERI_F3z_S_Dyz_Py_C1002001003;
  abcd[690] = I_ERI_G4x_S_D2z_Py_C1002001003+ABX*I_ERI_F3x_S_D2z_Py_C1002001003;
  abcd[691] = I_ERI_G3xy_S_D2z_Py_C1002001003+ABX*I_ERI_F2xy_S_D2z_Py_C1002001003;
  abcd[692] = I_ERI_G3xz_S_D2z_Py_C1002001003+ABX*I_ERI_F2xz_S_D2z_Py_C1002001003;
  abcd[693] = I_ERI_G2x2y_S_D2z_Py_C1002001003+ABX*I_ERI_Fx2y_S_D2z_Py_C1002001003;
  abcd[694] = I_ERI_G2xyz_S_D2z_Py_C1002001003+ABX*I_ERI_Fxyz_S_D2z_Py_C1002001003;
  abcd[695] = I_ERI_G2x2z_S_D2z_Py_C1002001003+ABX*I_ERI_Fx2z_S_D2z_Py_C1002001003;
  abcd[696] = I_ERI_Gx3y_S_D2z_Py_C1002001003+ABX*I_ERI_F3y_S_D2z_Py_C1002001003;
  abcd[697] = I_ERI_Gx2yz_S_D2z_Py_C1002001003+ABX*I_ERI_F2yz_S_D2z_Py_C1002001003;
  abcd[698] = I_ERI_Gxy2z_S_D2z_Py_C1002001003+ABX*I_ERI_Fy2z_S_D2z_Py_C1002001003;
  abcd[699] = I_ERI_Gx3z_S_D2z_Py_C1002001003+ABX*I_ERI_F3z_S_D2z_Py_C1002001003;
  abcd[700] = I_ERI_G3xy_S_D2z_Py_C1002001003+ABY*I_ERI_F3x_S_D2z_Py_C1002001003;
  abcd[701] = I_ERI_G2x2y_S_D2z_Py_C1002001003+ABY*I_ERI_F2xy_S_D2z_Py_C1002001003;
  abcd[702] = I_ERI_G2xyz_S_D2z_Py_C1002001003+ABY*I_ERI_F2xz_S_D2z_Py_C1002001003;
  abcd[703] = I_ERI_Gx3y_S_D2z_Py_C1002001003+ABY*I_ERI_Fx2y_S_D2z_Py_C1002001003;
  abcd[704] = I_ERI_Gx2yz_S_D2z_Py_C1002001003+ABY*I_ERI_Fxyz_S_D2z_Py_C1002001003;
  abcd[705] = I_ERI_Gxy2z_S_D2z_Py_C1002001003+ABY*I_ERI_Fx2z_S_D2z_Py_C1002001003;
  abcd[706] = I_ERI_G4y_S_D2z_Py_C1002001003+ABY*I_ERI_F3y_S_D2z_Py_C1002001003;
  abcd[707] = I_ERI_G3yz_S_D2z_Py_C1002001003+ABY*I_ERI_F2yz_S_D2z_Py_C1002001003;
  abcd[708] = I_ERI_G2y2z_S_D2z_Py_C1002001003+ABY*I_ERI_Fy2z_S_D2z_Py_C1002001003;
  abcd[709] = I_ERI_Gy3z_S_D2z_Py_C1002001003+ABY*I_ERI_F3z_S_D2z_Py_C1002001003;
  abcd[710] = I_ERI_G3xz_S_D2z_Py_C1002001003+ABZ*I_ERI_F3x_S_D2z_Py_C1002001003;
  abcd[711] = I_ERI_G2xyz_S_D2z_Py_C1002001003+ABZ*I_ERI_F2xy_S_D2z_Py_C1002001003;
  abcd[712] = I_ERI_G2x2z_S_D2z_Py_C1002001003+ABZ*I_ERI_F2xz_S_D2z_Py_C1002001003;
  abcd[713] = I_ERI_Gx2yz_S_D2z_Py_C1002001003+ABZ*I_ERI_Fx2y_S_D2z_Py_C1002001003;
  abcd[714] = I_ERI_Gxy2z_S_D2z_Py_C1002001003+ABZ*I_ERI_Fxyz_S_D2z_Py_C1002001003;
  abcd[715] = I_ERI_Gx3z_S_D2z_Py_C1002001003+ABZ*I_ERI_Fx2z_S_D2z_Py_C1002001003;
  abcd[716] = I_ERI_G3yz_S_D2z_Py_C1002001003+ABZ*I_ERI_F3y_S_D2z_Py_C1002001003;
  abcd[717] = I_ERI_G2y2z_S_D2z_Py_C1002001003+ABZ*I_ERI_F2yz_S_D2z_Py_C1002001003;
  abcd[718] = I_ERI_Gy3z_S_D2z_Py_C1002001003+ABZ*I_ERI_Fy2z_S_D2z_Py_C1002001003;
  abcd[719] = I_ERI_G4z_S_D2z_Py_C1002001003+ABZ*I_ERI_F3z_S_D2z_Py_C1002001003;
  abcd[730] = I_ERI_G4x_S_D2x_Pz_C1002001003+ABX*I_ERI_F3x_S_D2x_Pz_C1002001003;
  abcd[731] = I_ERI_G3xy_S_D2x_Pz_C1002001003+ABX*I_ERI_F2xy_S_D2x_Pz_C1002001003;
  abcd[732] = I_ERI_G3xz_S_D2x_Pz_C1002001003+ABX*I_ERI_F2xz_S_D2x_Pz_C1002001003;
  abcd[733] = I_ERI_G2x2y_S_D2x_Pz_C1002001003+ABX*I_ERI_Fx2y_S_D2x_Pz_C1002001003;
  abcd[734] = I_ERI_G2xyz_S_D2x_Pz_C1002001003+ABX*I_ERI_Fxyz_S_D2x_Pz_C1002001003;
  abcd[735] = I_ERI_G2x2z_S_D2x_Pz_C1002001003+ABX*I_ERI_Fx2z_S_D2x_Pz_C1002001003;
  abcd[736] = I_ERI_Gx3y_S_D2x_Pz_C1002001003+ABX*I_ERI_F3y_S_D2x_Pz_C1002001003;
  abcd[737] = I_ERI_Gx2yz_S_D2x_Pz_C1002001003+ABX*I_ERI_F2yz_S_D2x_Pz_C1002001003;
  abcd[738] = I_ERI_Gxy2z_S_D2x_Pz_C1002001003+ABX*I_ERI_Fy2z_S_D2x_Pz_C1002001003;
  abcd[739] = I_ERI_Gx3z_S_D2x_Pz_C1002001003+ABX*I_ERI_F3z_S_D2x_Pz_C1002001003;
  abcd[740] = I_ERI_G3xy_S_D2x_Pz_C1002001003+ABY*I_ERI_F3x_S_D2x_Pz_C1002001003;
  abcd[741] = I_ERI_G2x2y_S_D2x_Pz_C1002001003+ABY*I_ERI_F2xy_S_D2x_Pz_C1002001003;
  abcd[742] = I_ERI_G2xyz_S_D2x_Pz_C1002001003+ABY*I_ERI_F2xz_S_D2x_Pz_C1002001003;
  abcd[743] = I_ERI_Gx3y_S_D2x_Pz_C1002001003+ABY*I_ERI_Fx2y_S_D2x_Pz_C1002001003;
  abcd[744] = I_ERI_Gx2yz_S_D2x_Pz_C1002001003+ABY*I_ERI_Fxyz_S_D2x_Pz_C1002001003;
  abcd[745] = I_ERI_Gxy2z_S_D2x_Pz_C1002001003+ABY*I_ERI_Fx2z_S_D2x_Pz_C1002001003;
  abcd[746] = I_ERI_G4y_S_D2x_Pz_C1002001003+ABY*I_ERI_F3y_S_D2x_Pz_C1002001003;
  abcd[747] = I_ERI_G3yz_S_D2x_Pz_C1002001003+ABY*I_ERI_F2yz_S_D2x_Pz_C1002001003;
  abcd[748] = I_ERI_G2y2z_S_D2x_Pz_C1002001003+ABY*I_ERI_Fy2z_S_D2x_Pz_C1002001003;
  abcd[749] = I_ERI_Gy3z_S_D2x_Pz_C1002001003+ABY*I_ERI_F3z_S_D2x_Pz_C1002001003;
  abcd[750] = I_ERI_G3xz_S_D2x_Pz_C1002001003+ABZ*I_ERI_F3x_S_D2x_Pz_C1002001003;
  abcd[751] = I_ERI_G2xyz_S_D2x_Pz_C1002001003+ABZ*I_ERI_F2xy_S_D2x_Pz_C1002001003;
  abcd[752] = I_ERI_G2x2z_S_D2x_Pz_C1002001003+ABZ*I_ERI_F2xz_S_D2x_Pz_C1002001003;
  abcd[753] = I_ERI_Gx2yz_S_D2x_Pz_C1002001003+ABZ*I_ERI_Fx2y_S_D2x_Pz_C1002001003;
  abcd[754] = I_ERI_Gxy2z_S_D2x_Pz_C1002001003+ABZ*I_ERI_Fxyz_S_D2x_Pz_C1002001003;
  abcd[755] = I_ERI_Gx3z_S_D2x_Pz_C1002001003+ABZ*I_ERI_Fx2z_S_D2x_Pz_C1002001003;
  abcd[756] = I_ERI_G3yz_S_D2x_Pz_C1002001003+ABZ*I_ERI_F3y_S_D2x_Pz_C1002001003;
  abcd[757] = I_ERI_G2y2z_S_D2x_Pz_C1002001003+ABZ*I_ERI_F2yz_S_D2x_Pz_C1002001003;
  abcd[758] = I_ERI_Gy3z_S_D2x_Pz_C1002001003+ABZ*I_ERI_Fy2z_S_D2x_Pz_C1002001003;
  abcd[759] = I_ERI_G4z_S_D2x_Pz_C1002001003+ABZ*I_ERI_F3z_S_D2x_Pz_C1002001003;
  abcd[770] = I_ERI_G4x_S_Dxy_Pz_C1002001003+ABX*I_ERI_F3x_S_Dxy_Pz_C1002001003;
  abcd[771] = I_ERI_G3xy_S_Dxy_Pz_C1002001003+ABX*I_ERI_F2xy_S_Dxy_Pz_C1002001003;
  abcd[772] = I_ERI_G3xz_S_Dxy_Pz_C1002001003+ABX*I_ERI_F2xz_S_Dxy_Pz_C1002001003;
  abcd[773] = I_ERI_G2x2y_S_Dxy_Pz_C1002001003+ABX*I_ERI_Fx2y_S_Dxy_Pz_C1002001003;
  abcd[774] = I_ERI_G2xyz_S_Dxy_Pz_C1002001003+ABX*I_ERI_Fxyz_S_Dxy_Pz_C1002001003;
  abcd[775] = I_ERI_G2x2z_S_Dxy_Pz_C1002001003+ABX*I_ERI_Fx2z_S_Dxy_Pz_C1002001003;
  abcd[776] = I_ERI_Gx3y_S_Dxy_Pz_C1002001003+ABX*I_ERI_F3y_S_Dxy_Pz_C1002001003;
  abcd[777] = I_ERI_Gx2yz_S_Dxy_Pz_C1002001003+ABX*I_ERI_F2yz_S_Dxy_Pz_C1002001003;
  abcd[778] = I_ERI_Gxy2z_S_Dxy_Pz_C1002001003+ABX*I_ERI_Fy2z_S_Dxy_Pz_C1002001003;
  abcd[779] = I_ERI_Gx3z_S_Dxy_Pz_C1002001003+ABX*I_ERI_F3z_S_Dxy_Pz_C1002001003;
  abcd[780] = I_ERI_G3xy_S_Dxy_Pz_C1002001003+ABY*I_ERI_F3x_S_Dxy_Pz_C1002001003;
  abcd[781] = I_ERI_G2x2y_S_Dxy_Pz_C1002001003+ABY*I_ERI_F2xy_S_Dxy_Pz_C1002001003;
  abcd[782] = I_ERI_G2xyz_S_Dxy_Pz_C1002001003+ABY*I_ERI_F2xz_S_Dxy_Pz_C1002001003;
  abcd[783] = I_ERI_Gx3y_S_Dxy_Pz_C1002001003+ABY*I_ERI_Fx2y_S_Dxy_Pz_C1002001003;
  abcd[784] = I_ERI_Gx2yz_S_Dxy_Pz_C1002001003+ABY*I_ERI_Fxyz_S_Dxy_Pz_C1002001003;
  abcd[785] = I_ERI_Gxy2z_S_Dxy_Pz_C1002001003+ABY*I_ERI_Fx2z_S_Dxy_Pz_C1002001003;
  abcd[786] = I_ERI_G4y_S_Dxy_Pz_C1002001003+ABY*I_ERI_F3y_S_Dxy_Pz_C1002001003;
  abcd[787] = I_ERI_G3yz_S_Dxy_Pz_C1002001003+ABY*I_ERI_F2yz_S_Dxy_Pz_C1002001003;
  abcd[788] = I_ERI_G2y2z_S_Dxy_Pz_C1002001003+ABY*I_ERI_Fy2z_S_Dxy_Pz_C1002001003;
  abcd[789] = I_ERI_Gy3z_S_Dxy_Pz_C1002001003+ABY*I_ERI_F3z_S_Dxy_Pz_C1002001003;
  abcd[790] = I_ERI_G3xz_S_Dxy_Pz_C1002001003+ABZ*I_ERI_F3x_S_Dxy_Pz_C1002001003;
  abcd[791] = I_ERI_G2xyz_S_Dxy_Pz_C1002001003+ABZ*I_ERI_F2xy_S_Dxy_Pz_C1002001003;
  abcd[792] = I_ERI_G2x2z_S_Dxy_Pz_C1002001003+ABZ*I_ERI_F2xz_S_Dxy_Pz_C1002001003;
  abcd[793] = I_ERI_Gx2yz_S_Dxy_Pz_C1002001003+ABZ*I_ERI_Fx2y_S_Dxy_Pz_C1002001003;
  abcd[794] = I_ERI_Gxy2z_S_Dxy_Pz_C1002001003+ABZ*I_ERI_Fxyz_S_Dxy_Pz_C1002001003;
  abcd[795] = I_ERI_Gx3z_S_Dxy_Pz_C1002001003+ABZ*I_ERI_Fx2z_S_Dxy_Pz_C1002001003;
  abcd[796] = I_ERI_G3yz_S_Dxy_Pz_C1002001003+ABZ*I_ERI_F3y_S_Dxy_Pz_C1002001003;
  abcd[797] = I_ERI_G2y2z_S_Dxy_Pz_C1002001003+ABZ*I_ERI_F2yz_S_Dxy_Pz_C1002001003;
  abcd[798] = I_ERI_Gy3z_S_Dxy_Pz_C1002001003+ABZ*I_ERI_Fy2z_S_Dxy_Pz_C1002001003;
  abcd[799] = I_ERI_G4z_S_Dxy_Pz_C1002001003+ABZ*I_ERI_F3z_S_Dxy_Pz_C1002001003;
  abcd[810] = I_ERI_G4x_S_Dxz_Pz_C1002001003+ABX*I_ERI_F3x_S_Dxz_Pz_C1002001003;
  abcd[811] = I_ERI_G3xy_S_Dxz_Pz_C1002001003+ABX*I_ERI_F2xy_S_Dxz_Pz_C1002001003;
  abcd[812] = I_ERI_G3xz_S_Dxz_Pz_C1002001003+ABX*I_ERI_F2xz_S_Dxz_Pz_C1002001003;
  abcd[813] = I_ERI_G2x2y_S_Dxz_Pz_C1002001003+ABX*I_ERI_Fx2y_S_Dxz_Pz_C1002001003;
  abcd[814] = I_ERI_G2xyz_S_Dxz_Pz_C1002001003+ABX*I_ERI_Fxyz_S_Dxz_Pz_C1002001003;
  abcd[815] = I_ERI_G2x2z_S_Dxz_Pz_C1002001003+ABX*I_ERI_Fx2z_S_Dxz_Pz_C1002001003;
  abcd[816] = I_ERI_Gx3y_S_Dxz_Pz_C1002001003+ABX*I_ERI_F3y_S_Dxz_Pz_C1002001003;
  abcd[817] = I_ERI_Gx2yz_S_Dxz_Pz_C1002001003+ABX*I_ERI_F2yz_S_Dxz_Pz_C1002001003;
  abcd[818] = I_ERI_Gxy2z_S_Dxz_Pz_C1002001003+ABX*I_ERI_Fy2z_S_Dxz_Pz_C1002001003;
  abcd[819] = I_ERI_Gx3z_S_Dxz_Pz_C1002001003+ABX*I_ERI_F3z_S_Dxz_Pz_C1002001003;
  abcd[820] = I_ERI_G3xy_S_Dxz_Pz_C1002001003+ABY*I_ERI_F3x_S_Dxz_Pz_C1002001003;
  abcd[821] = I_ERI_G2x2y_S_Dxz_Pz_C1002001003+ABY*I_ERI_F2xy_S_Dxz_Pz_C1002001003;
  abcd[822] = I_ERI_G2xyz_S_Dxz_Pz_C1002001003+ABY*I_ERI_F2xz_S_Dxz_Pz_C1002001003;
  abcd[823] = I_ERI_Gx3y_S_Dxz_Pz_C1002001003+ABY*I_ERI_Fx2y_S_Dxz_Pz_C1002001003;
  abcd[824] = I_ERI_Gx2yz_S_Dxz_Pz_C1002001003+ABY*I_ERI_Fxyz_S_Dxz_Pz_C1002001003;
  abcd[825] = I_ERI_Gxy2z_S_Dxz_Pz_C1002001003+ABY*I_ERI_Fx2z_S_Dxz_Pz_C1002001003;
  abcd[826] = I_ERI_G4y_S_Dxz_Pz_C1002001003+ABY*I_ERI_F3y_S_Dxz_Pz_C1002001003;
  abcd[827] = I_ERI_G3yz_S_Dxz_Pz_C1002001003+ABY*I_ERI_F2yz_S_Dxz_Pz_C1002001003;
  abcd[828] = I_ERI_G2y2z_S_Dxz_Pz_C1002001003+ABY*I_ERI_Fy2z_S_Dxz_Pz_C1002001003;
  abcd[829] = I_ERI_Gy3z_S_Dxz_Pz_C1002001003+ABY*I_ERI_F3z_S_Dxz_Pz_C1002001003;
  abcd[830] = I_ERI_G3xz_S_Dxz_Pz_C1002001003+ABZ*I_ERI_F3x_S_Dxz_Pz_C1002001003;
  abcd[831] = I_ERI_G2xyz_S_Dxz_Pz_C1002001003+ABZ*I_ERI_F2xy_S_Dxz_Pz_C1002001003;
  abcd[832] = I_ERI_G2x2z_S_Dxz_Pz_C1002001003+ABZ*I_ERI_F2xz_S_Dxz_Pz_C1002001003;
  abcd[833] = I_ERI_Gx2yz_S_Dxz_Pz_C1002001003+ABZ*I_ERI_Fx2y_S_Dxz_Pz_C1002001003;
  abcd[834] = I_ERI_Gxy2z_S_Dxz_Pz_C1002001003+ABZ*I_ERI_Fxyz_S_Dxz_Pz_C1002001003;
  abcd[835] = I_ERI_Gx3z_S_Dxz_Pz_C1002001003+ABZ*I_ERI_Fx2z_S_Dxz_Pz_C1002001003;
  abcd[836] = I_ERI_G3yz_S_Dxz_Pz_C1002001003+ABZ*I_ERI_F3y_S_Dxz_Pz_C1002001003;
  abcd[837] = I_ERI_G2y2z_S_Dxz_Pz_C1002001003+ABZ*I_ERI_F2yz_S_Dxz_Pz_C1002001003;
  abcd[838] = I_ERI_Gy3z_S_Dxz_Pz_C1002001003+ABZ*I_ERI_Fy2z_S_Dxz_Pz_C1002001003;
  abcd[839] = I_ERI_G4z_S_Dxz_Pz_C1002001003+ABZ*I_ERI_F3z_S_Dxz_Pz_C1002001003;
  abcd[850] = I_ERI_G4x_S_D2y_Pz_C1002001003+ABX*I_ERI_F3x_S_D2y_Pz_C1002001003;
  abcd[851] = I_ERI_G3xy_S_D2y_Pz_C1002001003+ABX*I_ERI_F2xy_S_D2y_Pz_C1002001003;
  abcd[852] = I_ERI_G3xz_S_D2y_Pz_C1002001003+ABX*I_ERI_F2xz_S_D2y_Pz_C1002001003;
  abcd[853] = I_ERI_G2x2y_S_D2y_Pz_C1002001003+ABX*I_ERI_Fx2y_S_D2y_Pz_C1002001003;
  abcd[854] = I_ERI_G2xyz_S_D2y_Pz_C1002001003+ABX*I_ERI_Fxyz_S_D2y_Pz_C1002001003;
  abcd[855] = I_ERI_G2x2z_S_D2y_Pz_C1002001003+ABX*I_ERI_Fx2z_S_D2y_Pz_C1002001003;
  abcd[856] = I_ERI_Gx3y_S_D2y_Pz_C1002001003+ABX*I_ERI_F3y_S_D2y_Pz_C1002001003;
  abcd[857] = I_ERI_Gx2yz_S_D2y_Pz_C1002001003+ABX*I_ERI_F2yz_S_D2y_Pz_C1002001003;
  abcd[858] = I_ERI_Gxy2z_S_D2y_Pz_C1002001003+ABX*I_ERI_Fy2z_S_D2y_Pz_C1002001003;
  abcd[859] = I_ERI_Gx3z_S_D2y_Pz_C1002001003+ABX*I_ERI_F3z_S_D2y_Pz_C1002001003;
  abcd[860] = I_ERI_G3xy_S_D2y_Pz_C1002001003+ABY*I_ERI_F3x_S_D2y_Pz_C1002001003;
  abcd[861] = I_ERI_G2x2y_S_D2y_Pz_C1002001003+ABY*I_ERI_F2xy_S_D2y_Pz_C1002001003;
  abcd[862] = I_ERI_G2xyz_S_D2y_Pz_C1002001003+ABY*I_ERI_F2xz_S_D2y_Pz_C1002001003;
  abcd[863] = I_ERI_Gx3y_S_D2y_Pz_C1002001003+ABY*I_ERI_Fx2y_S_D2y_Pz_C1002001003;
  abcd[864] = I_ERI_Gx2yz_S_D2y_Pz_C1002001003+ABY*I_ERI_Fxyz_S_D2y_Pz_C1002001003;
  abcd[865] = I_ERI_Gxy2z_S_D2y_Pz_C1002001003+ABY*I_ERI_Fx2z_S_D2y_Pz_C1002001003;
  abcd[866] = I_ERI_G4y_S_D2y_Pz_C1002001003+ABY*I_ERI_F3y_S_D2y_Pz_C1002001003;
  abcd[867] = I_ERI_G3yz_S_D2y_Pz_C1002001003+ABY*I_ERI_F2yz_S_D2y_Pz_C1002001003;
  abcd[868] = I_ERI_G2y2z_S_D2y_Pz_C1002001003+ABY*I_ERI_Fy2z_S_D2y_Pz_C1002001003;
  abcd[869] = I_ERI_Gy3z_S_D2y_Pz_C1002001003+ABY*I_ERI_F3z_S_D2y_Pz_C1002001003;
  abcd[870] = I_ERI_G3xz_S_D2y_Pz_C1002001003+ABZ*I_ERI_F3x_S_D2y_Pz_C1002001003;
  abcd[871] = I_ERI_G2xyz_S_D2y_Pz_C1002001003+ABZ*I_ERI_F2xy_S_D2y_Pz_C1002001003;
  abcd[872] = I_ERI_G2x2z_S_D2y_Pz_C1002001003+ABZ*I_ERI_F2xz_S_D2y_Pz_C1002001003;
  abcd[873] = I_ERI_Gx2yz_S_D2y_Pz_C1002001003+ABZ*I_ERI_Fx2y_S_D2y_Pz_C1002001003;
  abcd[874] = I_ERI_Gxy2z_S_D2y_Pz_C1002001003+ABZ*I_ERI_Fxyz_S_D2y_Pz_C1002001003;
  abcd[875] = I_ERI_Gx3z_S_D2y_Pz_C1002001003+ABZ*I_ERI_Fx2z_S_D2y_Pz_C1002001003;
  abcd[876] = I_ERI_G3yz_S_D2y_Pz_C1002001003+ABZ*I_ERI_F3y_S_D2y_Pz_C1002001003;
  abcd[877] = I_ERI_G2y2z_S_D2y_Pz_C1002001003+ABZ*I_ERI_F2yz_S_D2y_Pz_C1002001003;
  abcd[878] = I_ERI_Gy3z_S_D2y_Pz_C1002001003+ABZ*I_ERI_Fy2z_S_D2y_Pz_C1002001003;
  abcd[879] = I_ERI_G4z_S_D2y_Pz_C1002001003+ABZ*I_ERI_F3z_S_D2y_Pz_C1002001003;
  abcd[890] = I_ERI_G4x_S_Dyz_Pz_C1002001003+ABX*I_ERI_F3x_S_Dyz_Pz_C1002001003;
  abcd[891] = I_ERI_G3xy_S_Dyz_Pz_C1002001003+ABX*I_ERI_F2xy_S_Dyz_Pz_C1002001003;
  abcd[892] = I_ERI_G3xz_S_Dyz_Pz_C1002001003+ABX*I_ERI_F2xz_S_Dyz_Pz_C1002001003;
  abcd[893] = I_ERI_G2x2y_S_Dyz_Pz_C1002001003+ABX*I_ERI_Fx2y_S_Dyz_Pz_C1002001003;
  abcd[894] = I_ERI_G2xyz_S_Dyz_Pz_C1002001003+ABX*I_ERI_Fxyz_S_Dyz_Pz_C1002001003;
  abcd[895] = I_ERI_G2x2z_S_Dyz_Pz_C1002001003+ABX*I_ERI_Fx2z_S_Dyz_Pz_C1002001003;
  abcd[896] = I_ERI_Gx3y_S_Dyz_Pz_C1002001003+ABX*I_ERI_F3y_S_Dyz_Pz_C1002001003;
  abcd[897] = I_ERI_Gx2yz_S_Dyz_Pz_C1002001003+ABX*I_ERI_F2yz_S_Dyz_Pz_C1002001003;
  abcd[898] = I_ERI_Gxy2z_S_Dyz_Pz_C1002001003+ABX*I_ERI_Fy2z_S_Dyz_Pz_C1002001003;
  abcd[899] = I_ERI_Gx3z_S_Dyz_Pz_C1002001003+ABX*I_ERI_F3z_S_Dyz_Pz_C1002001003;
  abcd[900] = I_ERI_G3xy_S_Dyz_Pz_C1002001003+ABY*I_ERI_F3x_S_Dyz_Pz_C1002001003;
  abcd[901] = I_ERI_G2x2y_S_Dyz_Pz_C1002001003+ABY*I_ERI_F2xy_S_Dyz_Pz_C1002001003;
  abcd[902] = I_ERI_G2xyz_S_Dyz_Pz_C1002001003+ABY*I_ERI_F2xz_S_Dyz_Pz_C1002001003;
  abcd[903] = I_ERI_Gx3y_S_Dyz_Pz_C1002001003+ABY*I_ERI_Fx2y_S_Dyz_Pz_C1002001003;
  abcd[904] = I_ERI_Gx2yz_S_Dyz_Pz_C1002001003+ABY*I_ERI_Fxyz_S_Dyz_Pz_C1002001003;
  abcd[905] = I_ERI_Gxy2z_S_Dyz_Pz_C1002001003+ABY*I_ERI_Fx2z_S_Dyz_Pz_C1002001003;
  abcd[906] = I_ERI_G4y_S_Dyz_Pz_C1002001003+ABY*I_ERI_F3y_S_Dyz_Pz_C1002001003;
  abcd[907] = I_ERI_G3yz_S_Dyz_Pz_C1002001003+ABY*I_ERI_F2yz_S_Dyz_Pz_C1002001003;
  abcd[908] = I_ERI_G2y2z_S_Dyz_Pz_C1002001003+ABY*I_ERI_Fy2z_S_Dyz_Pz_C1002001003;
  abcd[909] = I_ERI_Gy3z_S_Dyz_Pz_C1002001003+ABY*I_ERI_F3z_S_Dyz_Pz_C1002001003;
  abcd[910] = I_ERI_G3xz_S_Dyz_Pz_C1002001003+ABZ*I_ERI_F3x_S_Dyz_Pz_C1002001003;
  abcd[911] = I_ERI_G2xyz_S_Dyz_Pz_C1002001003+ABZ*I_ERI_F2xy_S_Dyz_Pz_C1002001003;
  abcd[912] = I_ERI_G2x2z_S_Dyz_Pz_C1002001003+ABZ*I_ERI_F2xz_S_Dyz_Pz_C1002001003;
  abcd[913] = I_ERI_Gx2yz_S_Dyz_Pz_C1002001003+ABZ*I_ERI_Fx2y_S_Dyz_Pz_C1002001003;
  abcd[914] = I_ERI_Gxy2z_S_Dyz_Pz_C1002001003+ABZ*I_ERI_Fxyz_S_Dyz_Pz_C1002001003;
  abcd[915] = I_ERI_Gx3z_S_Dyz_Pz_C1002001003+ABZ*I_ERI_Fx2z_S_Dyz_Pz_C1002001003;
  abcd[916] = I_ERI_G3yz_S_Dyz_Pz_C1002001003+ABZ*I_ERI_F3y_S_Dyz_Pz_C1002001003;
  abcd[917] = I_ERI_G2y2z_S_Dyz_Pz_C1002001003+ABZ*I_ERI_F2yz_S_Dyz_Pz_C1002001003;
  abcd[918] = I_ERI_Gy3z_S_Dyz_Pz_C1002001003+ABZ*I_ERI_Fy2z_S_Dyz_Pz_C1002001003;
  abcd[919] = I_ERI_G4z_S_Dyz_Pz_C1002001003+ABZ*I_ERI_F3z_S_Dyz_Pz_C1002001003;
  abcd[930] = I_ERI_G4x_S_D2z_Pz_C1002001003+ABX*I_ERI_F3x_S_D2z_Pz_C1002001003;
  abcd[931] = I_ERI_G3xy_S_D2z_Pz_C1002001003+ABX*I_ERI_F2xy_S_D2z_Pz_C1002001003;
  abcd[932] = I_ERI_G3xz_S_D2z_Pz_C1002001003+ABX*I_ERI_F2xz_S_D2z_Pz_C1002001003;
  abcd[933] = I_ERI_G2x2y_S_D2z_Pz_C1002001003+ABX*I_ERI_Fx2y_S_D2z_Pz_C1002001003;
  abcd[934] = I_ERI_G2xyz_S_D2z_Pz_C1002001003+ABX*I_ERI_Fxyz_S_D2z_Pz_C1002001003;
  abcd[935] = I_ERI_G2x2z_S_D2z_Pz_C1002001003+ABX*I_ERI_Fx2z_S_D2z_Pz_C1002001003;
  abcd[936] = I_ERI_Gx3y_S_D2z_Pz_C1002001003+ABX*I_ERI_F3y_S_D2z_Pz_C1002001003;
  abcd[937] = I_ERI_Gx2yz_S_D2z_Pz_C1002001003+ABX*I_ERI_F2yz_S_D2z_Pz_C1002001003;
  abcd[938] = I_ERI_Gxy2z_S_D2z_Pz_C1002001003+ABX*I_ERI_Fy2z_S_D2z_Pz_C1002001003;
  abcd[939] = I_ERI_Gx3z_S_D2z_Pz_C1002001003+ABX*I_ERI_F3z_S_D2z_Pz_C1002001003;
  abcd[940] = I_ERI_G3xy_S_D2z_Pz_C1002001003+ABY*I_ERI_F3x_S_D2z_Pz_C1002001003;
  abcd[941] = I_ERI_G2x2y_S_D2z_Pz_C1002001003+ABY*I_ERI_F2xy_S_D2z_Pz_C1002001003;
  abcd[942] = I_ERI_G2xyz_S_D2z_Pz_C1002001003+ABY*I_ERI_F2xz_S_D2z_Pz_C1002001003;
  abcd[943] = I_ERI_Gx3y_S_D2z_Pz_C1002001003+ABY*I_ERI_Fx2y_S_D2z_Pz_C1002001003;
  abcd[944] = I_ERI_Gx2yz_S_D2z_Pz_C1002001003+ABY*I_ERI_Fxyz_S_D2z_Pz_C1002001003;
  abcd[945] = I_ERI_Gxy2z_S_D2z_Pz_C1002001003+ABY*I_ERI_Fx2z_S_D2z_Pz_C1002001003;
  abcd[946] = I_ERI_G4y_S_D2z_Pz_C1002001003+ABY*I_ERI_F3y_S_D2z_Pz_C1002001003;
  abcd[947] = I_ERI_G3yz_S_D2z_Pz_C1002001003+ABY*I_ERI_F2yz_S_D2z_Pz_C1002001003;
  abcd[948] = I_ERI_G2y2z_S_D2z_Pz_C1002001003+ABY*I_ERI_Fy2z_S_D2z_Pz_C1002001003;
  abcd[949] = I_ERI_Gy3z_S_D2z_Pz_C1002001003+ABY*I_ERI_F3z_S_D2z_Pz_C1002001003;
  abcd[950] = I_ERI_G3xz_S_D2z_Pz_C1002001003+ABZ*I_ERI_F3x_S_D2z_Pz_C1002001003;
  abcd[951] = I_ERI_G2xyz_S_D2z_Pz_C1002001003+ABZ*I_ERI_F2xy_S_D2z_Pz_C1002001003;
  abcd[952] = I_ERI_G2x2z_S_D2z_Pz_C1002001003+ABZ*I_ERI_F2xz_S_D2z_Pz_C1002001003;
  abcd[953] = I_ERI_Gx2yz_S_D2z_Pz_C1002001003+ABZ*I_ERI_Fx2y_S_D2z_Pz_C1002001003;
  abcd[954] = I_ERI_Gxy2z_S_D2z_Pz_C1002001003+ABZ*I_ERI_Fxyz_S_D2z_Pz_C1002001003;
  abcd[955] = I_ERI_Gx3z_S_D2z_Pz_C1002001003+ABZ*I_ERI_Fx2z_S_D2z_Pz_C1002001003;
  abcd[956] = I_ERI_G3yz_S_D2z_Pz_C1002001003+ABZ*I_ERI_F3y_S_D2z_Pz_C1002001003;
  abcd[957] = I_ERI_G2y2z_S_D2z_Pz_C1002001003+ABZ*I_ERI_F2yz_S_D2z_Pz_C1002001003;
  abcd[958] = I_ERI_Gy3z_S_D2z_Pz_C1002001003+ABZ*I_ERI_Fy2z_S_D2z_Pz_C1002001003;
  abcd[959] = I_ERI_G4z_S_D2z_Pz_C1002001003+ABZ*I_ERI_F3z_S_D2z_Pz_C1002001003;
}
