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

void hgp_os_eri_d_d_f_p(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
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
  Double I_ERI_G4x_S_F3x_S = 0.0E0;
  Double I_ERI_G3xy_S_F3x_S = 0.0E0;
  Double I_ERI_G3xz_S_F3x_S = 0.0E0;
  Double I_ERI_G2x2y_S_F3x_S = 0.0E0;
  Double I_ERI_G2xyz_S_F3x_S = 0.0E0;
  Double I_ERI_G2x2z_S_F3x_S = 0.0E0;
  Double I_ERI_Gx3y_S_F3x_S = 0.0E0;
  Double I_ERI_Gx2yz_S_F3x_S = 0.0E0;
  Double I_ERI_Gxy2z_S_F3x_S = 0.0E0;
  Double I_ERI_Gx3z_S_F3x_S = 0.0E0;
  Double I_ERI_G4y_S_F3x_S = 0.0E0;
  Double I_ERI_G3yz_S_F3x_S = 0.0E0;
  Double I_ERI_G2y2z_S_F3x_S = 0.0E0;
  Double I_ERI_Gy3z_S_F3x_S = 0.0E0;
  Double I_ERI_G4z_S_F3x_S = 0.0E0;
  Double I_ERI_G4x_S_F2xy_S = 0.0E0;
  Double I_ERI_G3xy_S_F2xy_S = 0.0E0;
  Double I_ERI_G3xz_S_F2xy_S = 0.0E0;
  Double I_ERI_G2x2y_S_F2xy_S = 0.0E0;
  Double I_ERI_G2xyz_S_F2xy_S = 0.0E0;
  Double I_ERI_G2x2z_S_F2xy_S = 0.0E0;
  Double I_ERI_Gx3y_S_F2xy_S = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xy_S = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xy_S = 0.0E0;
  Double I_ERI_Gx3z_S_F2xy_S = 0.0E0;
  Double I_ERI_G4y_S_F2xy_S = 0.0E0;
  Double I_ERI_G3yz_S_F2xy_S = 0.0E0;
  Double I_ERI_G2y2z_S_F2xy_S = 0.0E0;
  Double I_ERI_Gy3z_S_F2xy_S = 0.0E0;
  Double I_ERI_G4z_S_F2xy_S = 0.0E0;
  Double I_ERI_G4x_S_F2xz_S = 0.0E0;
  Double I_ERI_G3xy_S_F2xz_S = 0.0E0;
  Double I_ERI_G3xz_S_F2xz_S = 0.0E0;
  Double I_ERI_G2x2y_S_F2xz_S = 0.0E0;
  Double I_ERI_G2xyz_S_F2xz_S = 0.0E0;
  Double I_ERI_G2x2z_S_F2xz_S = 0.0E0;
  Double I_ERI_Gx3y_S_F2xz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_F2xz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_F2xz_S = 0.0E0;
  Double I_ERI_Gx3z_S_F2xz_S = 0.0E0;
  Double I_ERI_G4y_S_F2xz_S = 0.0E0;
  Double I_ERI_G3yz_S_F2xz_S = 0.0E0;
  Double I_ERI_G2y2z_S_F2xz_S = 0.0E0;
  Double I_ERI_Gy3z_S_F2xz_S = 0.0E0;
  Double I_ERI_G4z_S_F2xz_S = 0.0E0;
  Double I_ERI_G4x_S_Fx2y_S = 0.0E0;
  Double I_ERI_G3xy_S_Fx2y_S = 0.0E0;
  Double I_ERI_G3xz_S_Fx2y_S = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2y_S = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2y_S = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2y_S = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2y_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2y_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2y_S = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2y_S = 0.0E0;
  Double I_ERI_G4y_S_Fx2y_S = 0.0E0;
  Double I_ERI_G3yz_S_Fx2y_S = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2y_S = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2y_S = 0.0E0;
  Double I_ERI_G4z_S_Fx2y_S = 0.0E0;
  Double I_ERI_G4x_S_Fxyz_S = 0.0E0;
  Double I_ERI_G3xy_S_Fxyz_S = 0.0E0;
  Double I_ERI_G3xz_S_Fxyz_S = 0.0E0;
  Double I_ERI_G2x2y_S_Fxyz_S = 0.0E0;
  Double I_ERI_G2xyz_S_Fxyz_S = 0.0E0;
  Double I_ERI_G2x2z_S_Fxyz_S = 0.0E0;
  Double I_ERI_Gx3y_S_Fxyz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Fxyz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Fxyz_S = 0.0E0;
  Double I_ERI_Gx3z_S_Fxyz_S = 0.0E0;
  Double I_ERI_G4y_S_Fxyz_S = 0.0E0;
  Double I_ERI_G3yz_S_Fxyz_S = 0.0E0;
  Double I_ERI_G2y2z_S_Fxyz_S = 0.0E0;
  Double I_ERI_Gy3z_S_Fxyz_S = 0.0E0;
  Double I_ERI_G4z_S_Fxyz_S = 0.0E0;
  Double I_ERI_G4x_S_Fx2z_S = 0.0E0;
  Double I_ERI_G3xy_S_Fx2z_S = 0.0E0;
  Double I_ERI_G3xz_S_Fx2z_S = 0.0E0;
  Double I_ERI_G2x2y_S_Fx2z_S = 0.0E0;
  Double I_ERI_G2xyz_S_Fx2z_S = 0.0E0;
  Double I_ERI_G2x2z_S_Fx2z_S = 0.0E0;
  Double I_ERI_Gx3y_S_Fx2z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Fx2z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Fx2z_S = 0.0E0;
  Double I_ERI_Gx3z_S_Fx2z_S = 0.0E0;
  Double I_ERI_G4y_S_Fx2z_S = 0.0E0;
  Double I_ERI_G3yz_S_Fx2z_S = 0.0E0;
  Double I_ERI_G2y2z_S_Fx2z_S = 0.0E0;
  Double I_ERI_Gy3z_S_Fx2z_S = 0.0E0;
  Double I_ERI_G4z_S_Fx2z_S = 0.0E0;
  Double I_ERI_G4x_S_F3y_S = 0.0E0;
  Double I_ERI_G3xy_S_F3y_S = 0.0E0;
  Double I_ERI_G3xz_S_F3y_S = 0.0E0;
  Double I_ERI_G2x2y_S_F3y_S = 0.0E0;
  Double I_ERI_G2xyz_S_F3y_S = 0.0E0;
  Double I_ERI_G2x2z_S_F3y_S = 0.0E0;
  Double I_ERI_Gx3y_S_F3y_S = 0.0E0;
  Double I_ERI_Gx2yz_S_F3y_S = 0.0E0;
  Double I_ERI_Gxy2z_S_F3y_S = 0.0E0;
  Double I_ERI_Gx3z_S_F3y_S = 0.0E0;
  Double I_ERI_G4y_S_F3y_S = 0.0E0;
  Double I_ERI_G3yz_S_F3y_S = 0.0E0;
  Double I_ERI_G2y2z_S_F3y_S = 0.0E0;
  Double I_ERI_Gy3z_S_F3y_S = 0.0E0;
  Double I_ERI_G4z_S_F3y_S = 0.0E0;
  Double I_ERI_G4x_S_F2yz_S = 0.0E0;
  Double I_ERI_G3xy_S_F2yz_S = 0.0E0;
  Double I_ERI_G3xz_S_F2yz_S = 0.0E0;
  Double I_ERI_G2x2y_S_F2yz_S = 0.0E0;
  Double I_ERI_G2xyz_S_F2yz_S = 0.0E0;
  Double I_ERI_G2x2z_S_F2yz_S = 0.0E0;
  Double I_ERI_Gx3y_S_F2yz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_F2yz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_F2yz_S = 0.0E0;
  Double I_ERI_Gx3z_S_F2yz_S = 0.0E0;
  Double I_ERI_G4y_S_F2yz_S = 0.0E0;
  Double I_ERI_G3yz_S_F2yz_S = 0.0E0;
  Double I_ERI_G2y2z_S_F2yz_S = 0.0E0;
  Double I_ERI_Gy3z_S_F2yz_S = 0.0E0;
  Double I_ERI_G4z_S_F2yz_S = 0.0E0;
  Double I_ERI_G4x_S_Fy2z_S = 0.0E0;
  Double I_ERI_G3xy_S_Fy2z_S = 0.0E0;
  Double I_ERI_G3xz_S_Fy2z_S = 0.0E0;
  Double I_ERI_G2x2y_S_Fy2z_S = 0.0E0;
  Double I_ERI_G2xyz_S_Fy2z_S = 0.0E0;
  Double I_ERI_G2x2z_S_Fy2z_S = 0.0E0;
  Double I_ERI_Gx3y_S_Fy2z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Fy2z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Fy2z_S = 0.0E0;
  Double I_ERI_Gx3z_S_Fy2z_S = 0.0E0;
  Double I_ERI_G4y_S_Fy2z_S = 0.0E0;
  Double I_ERI_G3yz_S_Fy2z_S = 0.0E0;
  Double I_ERI_G2y2z_S_Fy2z_S = 0.0E0;
  Double I_ERI_Gy3z_S_Fy2z_S = 0.0E0;
  Double I_ERI_G4z_S_Fy2z_S = 0.0E0;
  Double I_ERI_G4x_S_F3z_S = 0.0E0;
  Double I_ERI_G3xy_S_F3z_S = 0.0E0;
  Double I_ERI_G3xz_S_F3z_S = 0.0E0;
  Double I_ERI_G2x2y_S_F3z_S = 0.0E0;
  Double I_ERI_G2xyz_S_F3z_S = 0.0E0;
  Double I_ERI_G2x2z_S_F3z_S = 0.0E0;
  Double I_ERI_Gx3y_S_F3z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_F3z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_F3z_S = 0.0E0;
  Double I_ERI_Gx3z_S_F3z_S = 0.0E0;
  Double I_ERI_G4y_S_F3z_S = 0.0E0;
  Double I_ERI_G3yz_S_F3z_S = 0.0E0;
  Double I_ERI_G2y2z_S_F3z_S = 0.0E0;
  Double I_ERI_Gy3z_S_F3z_S = 0.0E0;
  Double I_ERI_G4z_S_F3z_S = 0.0E0;
  Double I_ERI_F3x_S_G4x_S = 0.0E0;
  Double I_ERI_F2xy_S_G4x_S = 0.0E0;
  Double I_ERI_F2xz_S_G4x_S = 0.0E0;
  Double I_ERI_Fx2y_S_G4x_S = 0.0E0;
  Double I_ERI_Fxyz_S_G4x_S = 0.0E0;
  Double I_ERI_Fx2z_S_G4x_S = 0.0E0;
  Double I_ERI_F3y_S_G4x_S = 0.0E0;
  Double I_ERI_F2yz_S_G4x_S = 0.0E0;
  Double I_ERI_Fy2z_S_G4x_S = 0.0E0;
  Double I_ERI_F3z_S_G4x_S = 0.0E0;
  Double I_ERI_F3x_S_G3xy_S = 0.0E0;
  Double I_ERI_F2xy_S_G3xy_S = 0.0E0;
  Double I_ERI_F2xz_S_G3xy_S = 0.0E0;
  Double I_ERI_Fx2y_S_G3xy_S = 0.0E0;
  Double I_ERI_Fxyz_S_G3xy_S = 0.0E0;
  Double I_ERI_Fx2z_S_G3xy_S = 0.0E0;
  Double I_ERI_F3y_S_G3xy_S = 0.0E0;
  Double I_ERI_F2yz_S_G3xy_S = 0.0E0;
  Double I_ERI_Fy2z_S_G3xy_S = 0.0E0;
  Double I_ERI_F3z_S_G3xy_S = 0.0E0;
  Double I_ERI_F3x_S_G3xz_S = 0.0E0;
  Double I_ERI_F2xy_S_G3xz_S = 0.0E0;
  Double I_ERI_F2xz_S_G3xz_S = 0.0E0;
  Double I_ERI_Fx2y_S_G3xz_S = 0.0E0;
  Double I_ERI_Fxyz_S_G3xz_S = 0.0E0;
  Double I_ERI_Fx2z_S_G3xz_S = 0.0E0;
  Double I_ERI_F3y_S_G3xz_S = 0.0E0;
  Double I_ERI_F2yz_S_G3xz_S = 0.0E0;
  Double I_ERI_Fy2z_S_G3xz_S = 0.0E0;
  Double I_ERI_F3z_S_G3xz_S = 0.0E0;
  Double I_ERI_F3x_S_G2x2y_S = 0.0E0;
  Double I_ERI_F2xy_S_G2x2y_S = 0.0E0;
  Double I_ERI_F2xz_S_G2x2y_S = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2y_S = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2y_S = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_F3y_S_G2x2y_S = 0.0E0;
  Double I_ERI_F2yz_S_G2x2y_S = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_F3z_S_G2x2y_S = 0.0E0;
  Double I_ERI_F3x_S_G2xyz_S = 0.0E0;
  Double I_ERI_F2xy_S_G2xyz_S = 0.0E0;
  Double I_ERI_F2xz_S_G2xyz_S = 0.0E0;
  Double I_ERI_Fx2y_S_G2xyz_S = 0.0E0;
  Double I_ERI_Fxyz_S_G2xyz_S = 0.0E0;
  Double I_ERI_Fx2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_F3y_S_G2xyz_S = 0.0E0;
  Double I_ERI_F2yz_S_G2xyz_S = 0.0E0;
  Double I_ERI_Fy2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_F3z_S_G2xyz_S = 0.0E0;
  Double I_ERI_F3x_S_G2x2z_S = 0.0E0;
  Double I_ERI_F2xy_S_G2x2z_S = 0.0E0;
  Double I_ERI_F2xz_S_G2x2z_S = 0.0E0;
  Double I_ERI_Fx2y_S_G2x2z_S = 0.0E0;
  Double I_ERI_Fxyz_S_G2x2z_S = 0.0E0;
  Double I_ERI_Fx2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_F3y_S_G2x2z_S = 0.0E0;
  Double I_ERI_F2yz_S_G2x2z_S = 0.0E0;
  Double I_ERI_Fy2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_F3z_S_G2x2z_S = 0.0E0;
  Double I_ERI_F3x_S_Gx3y_S = 0.0E0;
  Double I_ERI_F2xy_S_Gx3y_S = 0.0E0;
  Double I_ERI_F2xz_S_Gx3y_S = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3y_S = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3y_S = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_F3y_S_Gx3y_S = 0.0E0;
  Double I_ERI_F2yz_S_Gx3y_S = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_F3z_S_Gx3y_S = 0.0E0;
  Double I_ERI_F3x_S_Gx2yz_S = 0.0E0;
  Double I_ERI_F2xy_S_Gx2yz_S = 0.0E0;
  Double I_ERI_F2xz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Fx2y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Fxyz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Fx2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_F3y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_F2yz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Fy2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_F3z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_F3x_S_Gxy2z_S = 0.0E0;
  Double I_ERI_F2xy_S_Gxy2z_S = 0.0E0;
  Double I_ERI_F2xz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Fx2y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Fxyz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Fx2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_F3y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_F2yz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Fy2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_F3z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_F3x_S_Gx3z_S = 0.0E0;
  Double I_ERI_F2xy_S_Gx3z_S = 0.0E0;
  Double I_ERI_F2xz_S_Gx3z_S = 0.0E0;
  Double I_ERI_Fx2y_S_Gx3z_S = 0.0E0;
  Double I_ERI_Fxyz_S_Gx3z_S = 0.0E0;
  Double I_ERI_Fx2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_F3y_S_Gx3z_S = 0.0E0;
  Double I_ERI_F2yz_S_Gx3z_S = 0.0E0;
  Double I_ERI_Fy2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_F3z_S_Gx3z_S = 0.0E0;
  Double I_ERI_F3x_S_G4y_S = 0.0E0;
  Double I_ERI_F2xy_S_G4y_S = 0.0E0;
  Double I_ERI_F2xz_S_G4y_S = 0.0E0;
  Double I_ERI_Fx2y_S_G4y_S = 0.0E0;
  Double I_ERI_Fxyz_S_G4y_S = 0.0E0;
  Double I_ERI_Fx2z_S_G4y_S = 0.0E0;
  Double I_ERI_F3y_S_G4y_S = 0.0E0;
  Double I_ERI_F2yz_S_G4y_S = 0.0E0;
  Double I_ERI_Fy2z_S_G4y_S = 0.0E0;
  Double I_ERI_F3z_S_G4y_S = 0.0E0;
  Double I_ERI_F3x_S_G3yz_S = 0.0E0;
  Double I_ERI_F2xy_S_G3yz_S = 0.0E0;
  Double I_ERI_F2xz_S_G3yz_S = 0.0E0;
  Double I_ERI_Fx2y_S_G3yz_S = 0.0E0;
  Double I_ERI_Fxyz_S_G3yz_S = 0.0E0;
  Double I_ERI_Fx2z_S_G3yz_S = 0.0E0;
  Double I_ERI_F3y_S_G3yz_S = 0.0E0;
  Double I_ERI_F2yz_S_G3yz_S = 0.0E0;
  Double I_ERI_Fy2z_S_G3yz_S = 0.0E0;
  Double I_ERI_F3z_S_G3yz_S = 0.0E0;
  Double I_ERI_F3x_S_G2y2z_S = 0.0E0;
  Double I_ERI_F2xy_S_G2y2z_S = 0.0E0;
  Double I_ERI_F2xz_S_G2y2z_S = 0.0E0;
  Double I_ERI_Fx2y_S_G2y2z_S = 0.0E0;
  Double I_ERI_Fxyz_S_G2y2z_S = 0.0E0;
  Double I_ERI_Fx2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_F3y_S_G2y2z_S = 0.0E0;
  Double I_ERI_F2yz_S_G2y2z_S = 0.0E0;
  Double I_ERI_Fy2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_F3z_S_G2y2z_S = 0.0E0;
  Double I_ERI_F3x_S_Gy3z_S = 0.0E0;
  Double I_ERI_F2xy_S_Gy3z_S = 0.0E0;
  Double I_ERI_F2xz_S_Gy3z_S = 0.0E0;
  Double I_ERI_Fx2y_S_Gy3z_S = 0.0E0;
  Double I_ERI_Fxyz_S_Gy3z_S = 0.0E0;
  Double I_ERI_Fx2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_F3y_S_Gy3z_S = 0.0E0;
  Double I_ERI_F2yz_S_Gy3z_S = 0.0E0;
  Double I_ERI_Fy2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_F3z_S_Gy3z_S = 0.0E0;
  Double I_ERI_F3x_S_G4z_S = 0.0E0;
  Double I_ERI_F2xy_S_G4z_S = 0.0E0;
  Double I_ERI_F2xz_S_G4z_S = 0.0E0;
  Double I_ERI_Fx2y_S_G4z_S = 0.0E0;
  Double I_ERI_Fxyz_S_G4z_S = 0.0E0;
  Double I_ERI_Fx2z_S_G4z_S = 0.0E0;
  Double I_ERI_F3y_S_G4z_S = 0.0E0;
  Double I_ERI_F2yz_S_G4z_S = 0.0E0;
  Double I_ERI_Fy2z_S_G4z_S = 0.0E0;
  Double I_ERI_F3z_S_G4z_S = 0.0E0;
  Double I_ERI_F3x_S_F3x_S = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S = 0.0E0;
  Double I_ERI_F3y_S_F3x_S = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S = 0.0E0;
  Double I_ERI_F3z_S_F3x_S = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S = 0.0E0;
  Double I_ERI_F3x_S_F3y_S = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S = 0.0E0;
  Double I_ERI_F3y_S_F3y_S = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S = 0.0E0;
  Double I_ERI_F3z_S_F3y_S = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S = 0.0E0;
  Double I_ERI_F3x_S_F3z_S = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S = 0.0E0;
  Double I_ERI_F3y_S_F3z_S = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S = 0.0E0;
  Double I_ERI_F3z_S_F3z_S = 0.0E0;
  Double I_ERI_D2x_S_G4x_S = 0.0E0;
  Double I_ERI_Dxy_S_G4x_S = 0.0E0;
  Double I_ERI_Dxz_S_G4x_S = 0.0E0;
  Double I_ERI_D2y_S_G4x_S = 0.0E0;
  Double I_ERI_Dyz_S_G4x_S = 0.0E0;
  Double I_ERI_D2z_S_G4x_S = 0.0E0;
  Double I_ERI_D2x_S_G3xy_S = 0.0E0;
  Double I_ERI_Dxy_S_G3xy_S = 0.0E0;
  Double I_ERI_Dxz_S_G3xy_S = 0.0E0;
  Double I_ERI_D2y_S_G3xy_S = 0.0E0;
  Double I_ERI_Dyz_S_G3xy_S = 0.0E0;
  Double I_ERI_D2z_S_G3xy_S = 0.0E0;
  Double I_ERI_D2x_S_G3xz_S = 0.0E0;
  Double I_ERI_Dxy_S_G3xz_S = 0.0E0;
  Double I_ERI_Dxz_S_G3xz_S = 0.0E0;
  Double I_ERI_D2y_S_G3xz_S = 0.0E0;
  Double I_ERI_Dyz_S_G3xz_S = 0.0E0;
  Double I_ERI_D2z_S_G3xz_S = 0.0E0;
  Double I_ERI_D2x_S_G2x2y_S = 0.0E0;
  Double I_ERI_Dxy_S_G2x2y_S = 0.0E0;
  Double I_ERI_Dxz_S_G2x2y_S = 0.0E0;
  Double I_ERI_D2y_S_G2x2y_S = 0.0E0;
  Double I_ERI_Dyz_S_G2x2y_S = 0.0E0;
  Double I_ERI_D2z_S_G2x2y_S = 0.0E0;
  Double I_ERI_D2x_S_G2xyz_S = 0.0E0;
  Double I_ERI_Dxy_S_G2xyz_S = 0.0E0;
  Double I_ERI_Dxz_S_G2xyz_S = 0.0E0;
  Double I_ERI_D2y_S_G2xyz_S = 0.0E0;
  Double I_ERI_Dyz_S_G2xyz_S = 0.0E0;
  Double I_ERI_D2z_S_G2xyz_S = 0.0E0;
  Double I_ERI_D2x_S_G2x2z_S = 0.0E0;
  Double I_ERI_Dxy_S_G2x2z_S = 0.0E0;
  Double I_ERI_Dxz_S_G2x2z_S = 0.0E0;
  Double I_ERI_D2y_S_G2x2z_S = 0.0E0;
  Double I_ERI_Dyz_S_G2x2z_S = 0.0E0;
  Double I_ERI_D2z_S_G2x2z_S = 0.0E0;
  Double I_ERI_D2x_S_Gx3y_S = 0.0E0;
  Double I_ERI_Dxy_S_Gx3y_S = 0.0E0;
  Double I_ERI_Dxz_S_Gx3y_S = 0.0E0;
  Double I_ERI_D2y_S_Gx3y_S = 0.0E0;
  Double I_ERI_Dyz_S_Gx3y_S = 0.0E0;
  Double I_ERI_D2z_S_Gx3y_S = 0.0E0;
  Double I_ERI_D2x_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Dxy_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Dxz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_D2y_S_Gx2yz_S = 0.0E0;
  Double I_ERI_Dyz_S_Gx2yz_S = 0.0E0;
  Double I_ERI_D2z_S_Gx2yz_S = 0.0E0;
  Double I_ERI_D2x_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Dxy_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Dxz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_D2y_S_Gxy2z_S = 0.0E0;
  Double I_ERI_Dyz_S_Gxy2z_S = 0.0E0;
  Double I_ERI_D2z_S_Gxy2z_S = 0.0E0;
  Double I_ERI_D2x_S_Gx3z_S = 0.0E0;
  Double I_ERI_Dxy_S_Gx3z_S = 0.0E0;
  Double I_ERI_Dxz_S_Gx3z_S = 0.0E0;
  Double I_ERI_D2y_S_Gx3z_S = 0.0E0;
  Double I_ERI_Dyz_S_Gx3z_S = 0.0E0;
  Double I_ERI_D2z_S_Gx3z_S = 0.0E0;
  Double I_ERI_D2x_S_G4y_S = 0.0E0;
  Double I_ERI_Dxy_S_G4y_S = 0.0E0;
  Double I_ERI_Dxz_S_G4y_S = 0.0E0;
  Double I_ERI_D2y_S_G4y_S = 0.0E0;
  Double I_ERI_Dyz_S_G4y_S = 0.0E0;
  Double I_ERI_D2z_S_G4y_S = 0.0E0;
  Double I_ERI_D2x_S_G3yz_S = 0.0E0;
  Double I_ERI_Dxy_S_G3yz_S = 0.0E0;
  Double I_ERI_Dxz_S_G3yz_S = 0.0E0;
  Double I_ERI_D2y_S_G3yz_S = 0.0E0;
  Double I_ERI_Dyz_S_G3yz_S = 0.0E0;
  Double I_ERI_D2z_S_G3yz_S = 0.0E0;
  Double I_ERI_D2x_S_G2y2z_S = 0.0E0;
  Double I_ERI_Dxy_S_G2y2z_S = 0.0E0;
  Double I_ERI_Dxz_S_G2y2z_S = 0.0E0;
  Double I_ERI_D2y_S_G2y2z_S = 0.0E0;
  Double I_ERI_Dyz_S_G2y2z_S = 0.0E0;
  Double I_ERI_D2z_S_G2y2z_S = 0.0E0;
  Double I_ERI_D2x_S_Gy3z_S = 0.0E0;
  Double I_ERI_Dxy_S_Gy3z_S = 0.0E0;
  Double I_ERI_Dxz_S_Gy3z_S = 0.0E0;
  Double I_ERI_D2y_S_Gy3z_S = 0.0E0;
  Double I_ERI_Dyz_S_Gy3z_S = 0.0E0;
  Double I_ERI_D2z_S_Gy3z_S = 0.0E0;
  Double I_ERI_D2x_S_G4z_S = 0.0E0;
  Double I_ERI_Dxy_S_G4z_S = 0.0E0;
  Double I_ERI_Dxz_S_G4z_S = 0.0E0;
  Double I_ERI_D2y_S_G4z_S = 0.0E0;
  Double I_ERI_Dyz_S_G4z_S = 0.0E0;
  Double I_ERI_D2z_S_G4z_S = 0.0E0;
  Double I_ERI_D2x_S_F3x_S = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S = 0.0E0;
  Double I_ERI_D2y_S_F3x_S = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S = 0.0E0;
  Double I_ERI_D2z_S_F3x_S = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S = 0.0E0;
  Double I_ERI_D2x_S_F3y_S = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S = 0.0E0;
  Double I_ERI_D2y_S_F3y_S = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S = 0.0E0;
  Double I_ERI_D2z_S_F3y_S = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S = 0.0E0;
  Double I_ERI_D2x_S_F3z_S = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S = 0.0E0;
  Double I_ERI_D2y_S_F3z_S = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S = 0.0E0;
  Double I_ERI_D2z_S_F3z_S = 0.0E0;

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
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER51;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER49*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER47*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER45*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER43*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER25*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER23*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER21*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = 1.0E0+u2*ONEOVER19*I_ERI_S_S_S_S_M8_vrr;
        I_ERI_S_S_S_S_M8_vrr = ONEOVER17*I_ERI_S_S_S_S_M8_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M8_vrr  = f*I_ERI_S_S_S_S_M8_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
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
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
        I_ERI_S_S_S_S_M5_vrr = I_ERI_S_S_S_S_M5_vrr*erfPref_11;
        I_ERI_S_S_S_S_M6_vrr = I_ERI_S_S_S_S_M6_vrr*erfPref_13;
        I_ERI_S_S_S_S_M7_vrr = I_ERI_S_S_S_S_M7_vrr*erfPref_15;
        I_ERI_S_S_S_S_M8_vrr = I_ERI_S_S_S_S_M8_vrr*erfPref_17;
      }

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
       * shell quartet name: SQ_ERI_S_S_D_S_M6
       * expanding position: KET1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M6_vrr = QCX*I_ERI_S_S_Px_S_M6_vrr+WQX*I_ERI_S_S_Px_S_M7_vrr+oned2e*I_ERI_S_S_S_S_M6_vrr-rhod2esq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_D2y_S_M6_vrr = QCY*I_ERI_S_S_Py_S_M6_vrr+WQY*I_ERI_S_S_Py_S_M7_vrr+oned2e*I_ERI_S_S_S_S_M6_vrr-rhod2esq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_S_S_D2z_S_M6_vrr = QCZ*I_ERI_S_S_Pz_S_M6_vrr+WQZ*I_ERI_S_S_Pz_S_M7_vrr+oned2e*I_ERI_S_S_S_S_M6_vrr-rhod2esq*I_ERI_S_S_S_S_M7_vrr;

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
       * shell quartet name: SQ_ERI_S_S_D_S_M5
       * expanding position: KET1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M5_vrr = QCX*I_ERI_S_S_Px_S_M5_vrr+WQX*I_ERI_S_S_Px_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_Dxy_S_M5_vrr = QCY*I_ERI_S_S_Px_S_M5_vrr+WQY*I_ERI_S_S_Px_S_M6_vrr;
      Double I_ERI_S_S_D2y_S_M5_vrr = QCY*I_ERI_S_S_Py_S_M5_vrr+WQY*I_ERI_S_S_Py_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;
      Double I_ERI_S_S_D2z_S_M5_vrr = QCZ*I_ERI_S_S_Pz_S_M5_vrr+WQZ*I_ERI_S_S_Pz_S_M6_vrr+oned2e*I_ERI_S_S_S_S_M5_vrr-rhod2esq*I_ERI_S_S_S_S_M6_vrr;

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
       * shell quartet name: SQ_ERI_S_S_F_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       ************************************************************/
      Double I_ERI_S_S_F3x_S_M4_vrr = QCX*I_ERI_S_S_D2x_S_M4_vrr+WQX*I_ERI_S_S_D2x_S_M5_vrr+2*oned2e*I_ERI_S_S_Px_S_M4_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_F2xy_S_M4_vrr = QCY*I_ERI_S_S_D2x_S_M4_vrr+WQY*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_F2xz_S_M4_vrr = QCZ*I_ERI_S_S_D2x_S_M4_vrr+WQZ*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_Fx2y_S_M4_vrr = QCX*I_ERI_S_S_D2y_S_M4_vrr+WQX*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_Fxyz_S_M4_vrr = QCZ*I_ERI_S_S_Dxy_S_M4_vrr+WQZ*I_ERI_S_S_Dxy_S_M5_vrr;
      Double I_ERI_S_S_Fx2z_S_M4_vrr = QCX*I_ERI_S_S_D2z_S_M4_vrr+WQX*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_S_S_F3y_S_M4_vrr = QCY*I_ERI_S_S_D2y_S_M4_vrr+WQY*I_ERI_S_S_D2y_S_M5_vrr+2*oned2e*I_ERI_S_S_Py_S_M4_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M5_vrr;
      Double I_ERI_S_S_F2yz_S_M4_vrr = QCZ*I_ERI_S_S_D2y_S_M4_vrr+WQZ*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_Fy2z_S_M4_vrr = QCY*I_ERI_S_S_D2z_S_M4_vrr+WQY*I_ERI_S_S_D2z_S_M5_vrr;
      Double I_ERI_S_S_F3z_S_M4_vrr = QCZ*I_ERI_S_S_D2z_S_M4_vrr+WQZ*I_ERI_S_S_D2z_S_M5_vrr+2*oned2e*I_ERI_S_S_Pz_S_M4_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S_M4
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M5
       ************************************************************/
      Double I_ERI_S_S_G4x_S_M4_vrr = QCX*I_ERI_S_S_F3x_S_M4_vrr+WQX*I_ERI_S_S_F3x_S_M5_vrr+3*oned2e*I_ERI_S_S_D2x_S_M4_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_G3xy_S_M4_vrr = QCY*I_ERI_S_S_F3x_S_M4_vrr+WQY*I_ERI_S_S_F3x_S_M5_vrr;
      Double I_ERI_S_S_G3xz_S_M4_vrr = QCZ*I_ERI_S_S_F3x_S_M4_vrr+WQZ*I_ERI_S_S_F3x_S_M5_vrr;
      Double I_ERI_S_S_G2x2y_S_M4_vrr = QCY*I_ERI_S_S_F2xy_S_M4_vrr+WQY*I_ERI_S_S_F2xy_S_M5_vrr+oned2e*I_ERI_S_S_D2x_S_M4_vrr-rhod2esq*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_G2xyz_S_M4_vrr = QCZ*I_ERI_S_S_F2xy_S_M4_vrr+WQZ*I_ERI_S_S_F2xy_S_M5_vrr;
      Double I_ERI_S_S_G2x2z_S_M4_vrr = QCZ*I_ERI_S_S_F2xz_S_M4_vrr+WQZ*I_ERI_S_S_F2xz_S_M5_vrr+oned2e*I_ERI_S_S_D2x_S_M4_vrr-rhod2esq*I_ERI_S_S_D2x_S_M5_vrr;
      Double I_ERI_S_S_Gx3y_S_M4_vrr = QCX*I_ERI_S_S_F3y_S_M4_vrr+WQX*I_ERI_S_S_F3y_S_M5_vrr;
      Double I_ERI_S_S_Gx2yz_S_M4_vrr = QCZ*I_ERI_S_S_Fx2y_S_M4_vrr+WQZ*I_ERI_S_S_Fx2y_S_M5_vrr;
      Double I_ERI_S_S_Gxy2z_S_M4_vrr = QCY*I_ERI_S_S_Fx2z_S_M4_vrr+WQY*I_ERI_S_S_Fx2z_S_M5_vrr;
      Double I_ERI_S_S_Gx3z_S_M4_vrr = QCX*I_ERI_S_S_F3z_S_M4_vrr+WQX*I_ERI_S_S_F3z_S_M5_vrr;
      Double I_ERI_S_S_G4y_S_M4_vrr = QCY*I_ERI_S_S_F3y_S_M4_vrr+WQY*I_ERI_S_S_F3y_S_M5_vrr+3*oned2e*I_ERI_S_S_D2y_S_M4_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_G3yz_S_M4_vrr = QCZ*I_ERI_S_S_F3y_S_M4_vrr+WQZ*I_ERI_S_S_F3y_S_M5_vrr;
      Double I_ERI_S_S_G2y2z_S_M4_vrr = QCZ*I_ERI_S_S_F2yz_S_M4_vrr+WQZ*I_ERI_S_S_F2yz_S_M5_vrr+oned2e*I_ERI_S_S_D2y_S_M4_vrr-rhod2esq*I_ERI_S_S_D2y_S_M5_vrr;
      Double I_ERI_S_S_Gy3z_S_M4_vrr = QCY*I_ERI_S_S_F3z_S_M4_vrr+WQY*I_ERI_S_S_F3z_S_M5_vrr;
      Double I_ERI_S_S_G4z_S_M4_vrr = QCZ*I_ERI_S_S_F3z_S_M4_vrr+WQZ*I_ERI_S_S_F3z_S_M5_vrr+3*oned2e*I_ERI_S_S_D2z_S_M4_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_P_S_F_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_M3_vrr = PAX*I_ERI_S_S_F3x_S_M3_vrr+WPX*I_ERI_S_S_F3x_S_M4_vrr+3*oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Py_S_F3x_S_M3_vrr = PAY*I_ERI_S_S_F3x_S_M3_vrr+WPY*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Pz_S_F3x_S_M3_vrr = PAZ*I_ERI_S_S_F3x_S_M3_vrr+WPZ*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Px_S_F2xy_S_M3_vrr = PAX*I_ERI_S_S_F2xy_S_M3_vrr+WPX*I_ERI_S_S_F2xy_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Py_S_F2xy_S_M3_vrr = PAY*I_ERI_S_S_F2xy_S_M3_vrr+WPY*I_ERI_S_S_F2xy_S_M4_vrr+oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Pz_S_F2xy_S_M3_vrr = PAZ*I_ERI_S_S_F2xy_S_M3_vrr+WPZ*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_Px_S_F2xz_S_M3_vrr = PAX*I_ERI_S_S_F2xz_S_M3_vrr+WPX*I_ERI_S_S_F2xz_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Py_S_F2xz_S_M3_vrr = PAY*I_ERI_S_S_F2xz_S_M3_vrr+WPY*I_ERI_S_S_F2xz_S_M4_vrr;
      Double I_ERI_Pz_S_F2xz_S_M3_vrr = PAZ*I_ERI_S_S_F2xz_S_M3_vrr+WPZ*I_ERI_S_S_F2xz_S_M4_vrr+oned2k*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_Px_S_Fx2y_S_M3_vrr = PAX*I_ERI_S_S_Fx2y_S_M3_vrr+WPX*I_ERI_S_S_Fx2y_S_M4_vrr+oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Py_S_Fx2y_S_M3_vrr = PAY*I_ERI_S_S_Fx2y_S_M3_vrr+WPY*I_ERI_S_S_Fx2y_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Pz_S_Fx2y_S_M3_vrr = PAZ*I_ERI_S_S_Fx2y_S_M3_vrr+WPZ*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_Px_S_Fxyz_S_M3_vrr = PAX*I_ERI_S_S_Fxyz_S_M3_vrr+WPX*I_ERI_S_S_Fxyz_S_M4_vrr+oned2k*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Py_S_Fxyz_S_M3_vrr = PAY*I_ERI_S_S_Fxyz_S_M3_vrr+WPY*I_ERI_S_S_Fxyz_S_M4_vrr+oned2k*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Pz_S_Fxyz_S_M3_vrr = PAZ*I_ERI_S_S_Fxyz_S_M3_vrr+WPZ*I_ERI_S_S_Fxyz_S_M4_vrr+oned2k*I_ERI_S_S_Dxy_S_M4_vrr;
      Double I_ERI_Px_S_Fx2z_S_M3_vrr = PAX*I_ERI_S_S_Fx2z_S_M3_vrr+WPX*I_ERI_S_S_Fx2z_S_M4_vrr+oned2k*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Py_S_Fx2z_S_M3_vrr = PAY*I_ERI_S_S_Fx2z_S_M3_vrr+WPY*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_Pz_S_Fx2z_S_M3_vrr = PAZ*I_ERI_S_S_Fx2z_S_M3_vrr+WPZ*I_ERI_S_S_Fx2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M4_vrr;
      Double I_ERI_Px_S_F3y_S_M3_vrr = PAX*I_ERI_S_S_F3y_S_M3_vrr+WPX*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Py_S_F3y_S_M3_vrr = PAY*I_ERI_S_S_F3y_S_M3_vrr+WPY*I_ERI_S_S_F3y_S_M4_vrr+3*oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Pz_S_F3y_S_M3_vrr = PAZ*I_ERI_S_S_F3y_S_M3_vrr+WPZ*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Px_S_F2yz_S_M3_vrr = PAX*I_ERI_S_S_F2yz_S_M3_vrr+WPX*I_ERI_S_S_F2yz_S_M4_vrr;
      Double I_ERI_Py_S_F2yz_S_M3_vrr = PAY*I_ERI_S_S_F2yz_S_M3_vrr+WPY*I_ERI_S_S_F2yz_S_M4_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Pz_S_F2yz_S_M3_vrr = PAZ*I_ERI_S_S_F2yz_S_M3_vrr+WPZ*I_ERI_S_S_F2yz_S_M4_vrr+oned2k*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_Px_S_Fy2z_S_M3_vrr = PAX*I_ERI_S_S_Fy2z_S_M3_vrr+WPX*I_ERI_S_S_Fy2z_S_M4_vrr;
      Double I_ERI_Py_S_Fy2z_S_M3_vrr = PAY*I_ERI_S_S_Fy2z_S_M3_vrr+WPY*I_ERI_S_S_Fy2z_S_M4_vrr+oned2k*I_ERI_S_S_D2z_S_M4_vrr;
      Double I_ERI_Pz_S_Fy2z_S_M3_vrr = PAZ*I_ERI_S_S_Fy2z_S_M3_vrr+WPZ*I_ERI_S_S_Fy2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M4_vrr;
      Double I_ERI_Px_S_F3z_S_M3_vrr = PAX*I_ERI_S_S_F3z_S_M3_vrr+WPX*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Py_S_F3z_S_M3_vrr = PAY*I_ERI_S_S_F3z_S_M3_vrr+WPY*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Pz_S_F3z_S_M3_vrr = PAZ*I_ERI_S_S_F3z_S_M3_vrr+WPZ*I_ERI_S_S_F3z_S_M4_vrr+3*oned2k*I_ERI_S_S_D2z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_G_S_M3
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M4
       ************************************************************/
      Double I_ERI_S_S_G4x_S_M3_vrr = QCX*I_ERI_S_S_F3x_S_M3_vrr+WQX*I_ERI_S_S_F3x_S_M4_vrr+3*oned2e*I_ERI_S_S_D2x_S_M3_vrr-3*rhod2esq*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_G3xy_S_M3_vrr = QCY*I_ERI_S_S_F3x_S_M3_vrr+WQY*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_S_S_G3xz_S_M3_vrr = QCZ*I_ERI_S_S_F3x_S_M3_vrr+WQZ*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_S_S_G2x2y_S_M3_vrr = QCY*I_ERI_S_S_F2xy_S_M3_vrr+WQY*I_ERI_S_S_F2xy_S_M4_vrr+oned2e*I_ERI_S_S_D2x_S_M3_vrr-rhod2esq*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_G2xyz_S_M3_vrr = QCZ*I_ERI_S_S_F2xy_S_M3_vrr+WQZ*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_S_S_G2x2z_S_M3_vrr = QCZ*I_ERI_S_S_F2xz_S_M3_vrr+WQZ*I_ERI_S_S_F2xz_S_M4_vrr+oned2e*I_ERI_S_S_D2x_S_M3_vrr-rhod2esq*I_ERI_S_S_D2x_S_M4_vrr;
      Double I_ERI_S_S_Gx3y_S_M3_vrr = QCX*I_ERI_S_S_F3y_S_M3_vrr+WQX*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_S_S_Gx2yz_S_M3_vrr = QCZ*I_ERI_S_S_Fx2y_S_M3_vrr+WQZ*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_S_S_Gxy2z_S_M3_vrr = QCY*I_ERI_S_S_Fx2z_S_M3_vrr+WQY*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_S_S_Gx3z_S_M3_vrr = QCX*I_ERI_S_S_F3z_S_M3_vrr+WQX*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_S_S_G4y_S_M3_vrr = QCY*I_ERI_S_S_F3y_S_M3_vrr+WQY*I_ERI_S_S_F3y_S_M4_vrr+3*oned2e*I_ERI_S_S_D2y_S_M3_vrr-3*rhod2esq*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_G3yz_S_M3_vrr = QCZ*I_ERI_S_S_F3y_S_M3_vrr+WQZ*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_S_S_G2y2z_S_M3_vrr = QCZ*I_ERI_S_S_F2yz_S_M3_vrr+WQZ*I_ERI_S_S_F2yz_S_M4_vrr+oned2e*I_ERI_S_S_D2y_S_M3_vrr-rhod2esq*I_ERI_S_S_D2y_S_M4_vrr;
      Double I_ERI_S_S_Gy3z_S_M3_vrr = QCY*I_ERI_S_S_F3z_S_M3_vrr+WQY*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_S_S_G4z_S_M3_vrr = QCZ*I_ERI_S_S_F3z_S_M3_vrr+WQZ*I_ERI_S_S_F3z_S_M4_vrr+3*oned2e*I_ERI_S_S_D2z_S_M3_vrr-3*rhod2esq*I_ERI_S_S_D2z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_G_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M4
       ************************************************************/
      Double I_ERI_Px_S_G4x_S_M3_vrr = PAX*I_ERI_S_S_G4x_S_M3_vrr+WPX*I_ERI_S_S_G4x_S_M4_vrr+4*oned2k*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Py_S_G4x_S_M3_vrr = PAY*I_ERI_S_S_G4x_S_M3_vrr+WPY*I_ERI_S_S_G4x_S_M4_vrr;
      Double I_ERI_Pz_S_G4x_S_M3_vrr = PAZ*I_ERI_S_S_G4x_S_M3_vrr+WPZ*I_ERI_S_S_G4x_S_M4_vrr;
      Double I_ERI_Px_S_G3xy_S_M3_vrr = PAX*I_ERI_S_S_G3xy_S_M3_vrr+WPX*I_ERI_S_S_G3xy_S_M4_vrr+3*oned2k*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_Py_S_G3xy_S_M3_vrr = PAY*I_ERI_S_S_G3xy_S_M3_vrr+WPY*I_ERI_S_S_G3xy_S_M4_vrr+oned2k*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Pz_S_G3xy_S_M3_vrr = PAZ*I_ERI_S_S_G3xy_S_M3_vrr+WPZ*I_ERI_S_S_G3xy_S_M4_vrr;
      Double I_ERI_Px_S_G3xz_S_M3_vrr = PAX*I_ERI_S_S_G3xz_S_M3_vrr+WPX*I_ERI_S_S_G3xz_S_M4_vrr+3*oned2k*I_ERI_S_S_F2xz_S_M4_vrr;
      Double I_ERI_Py_S_G3xz_S_M3_vrr = PAY*I_ERI_S_S_G3xz_S_M3_vrr+WPY*I_ERI_S_S_G3xz_S_M4_vrr;
      Double I_ERI_Pz_S_G3xz_S_M3_vrr = PAZ*I_ERI_S_S_G3xz_S_M3_vrr+WPZ*I_ERI_S_S_G3xz_S_M4_vrr+oned2k*I_ERI_S_S_F3x_S_M4_vrr;
      Double I_ERI_Px_S_G2x2y_S_M3_vrr = PAX*I_ERI_S_S_G2x2y_S_M3_vrr+WPX*I_ERI_S_S_G2x2y_S_M4_vrr+2*oned2k*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_Py_S_G2x2y_S_M3_vrr = PAY*I_ERI_S_S_G2x2y_S_M3_vrr+WPY*I_ERI_S_S_G2x2y_S_M4_vrr+2*oned2k*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_Pz_S_G2x2y_S_M3_vrr = PAZ*I_ERI_S_S_G2x2y_S_M3_vrr+WPZ*I_ERI_S_S_G2x2y_S_M4_vrr;
      Double I_ERI_Px_S_G2xyz_S_M3_vrr = PAX*I_ERI_S_S_G2xyz_S_M3_vrr+WPX*I_ERI_S_S_G2xyz_S_M4_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M4_vrr;
      Double I_ERI_Py_S_G2xyz_S_M3_vrr = PAY*I_ERI_S_S_G2xyz_S_M3_vrr+WPY*I_ERI_S_S_G2xyz_S_M4_vrr+oned2k*I_ERI_S_S_F2xz_S_M4_vrr;
      Double I_ERI_Pz_S_G2xyz_S_M3_vrr = PAZ*I_ERI_S_S_G2xyz_S_M3_vrr+WPZ*I_ERI_S_S_G2xyz_S_M4_vrr+oned2k*I_ERI_S_S_F2xy_S_M4_vrr;
      Double I_ERI_Px_S_G2x2z_S_M3_vrr = PAX*I_ERI_S_S_G2x2z_S_M3_vrr+WPX*I_ERI_S_S_G2x2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_Py_S_G2x2z_S_M3_vrr = PAY*I_ERI_S_S_G2x2z_S_M3_vrr+WPY*I_ERI_S_S_G2x2z_S_M4_vrr;
      Double I_ERI_Pz_S_G2x2z_S_M3_vrr = PAZ*I_ERI_S_S_G2x2z_S_M3_vrr+WPZ*I_ERI_S_S_G2x2z_S_M4_vrr+2*oned2k*I_ERI_S_S_F2xz_S_M4_vrr;
      Double I_ERI_Px_S_Gx3y_S_M3_vrr = PAX*I_ERI_S_S_Gx3y_S_M3_vrr+WPX*I_ERI_S_S_Gx3y_S_M4_vrr+oned2k*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Py_S_Gx3y_S_M3_vrr = PAY*I_ERI_S_S_Gx3y_S_M3_vrr+WPY*I_ERI_S_S_Gx3y_S_M4_vrr+3*oned2k*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_Pz_S_Gx3y_S_M3_vrr = PAZ*I_ERI_S_S_Gx3y_S_M3_vrr+WPZ*I_ERI_S_S_Gx3y_S_M4_vrr;
      Double I_ERI_Px_S_Gx2yz_S_M3_vrr = PAX*I_ERI_S_S_Gx2yz_S_M3_vrr+WPX*I_ERI_S_S_Gx2yz_S_M4_vrr+oned2k*I_ERI_S_S_F2yz_S_M4_vrr;
      Double I_ERI_Py_S_Gx2yz_S_M3_vrr = PAY*I_ERI_S_S_Gx2yz_S_M3_vrr+WPY*I_ERI_S_S_Gx2yz_S_M4_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M4_vrr;
      Double I_ERI_Pz_S_Gx2yz_S_M3_vrr = PAZ*I_ERI_S_S_Gx2yz_S_M3_vrr+WPZ*I_ERI_S_S_Gx2yz_S_M4_vrr+oned2k*I_ERI_S_S_Fx2y_S_M4_vrr;
      Double I_ERI_Px_S_Gxy2z_S_M3_vrr = PAX*I_ERI_S_S_Gxy2z_S_M3_vrr+WPX*I_ERI_S_S_Gxy2z_S_M4_vrr+oned2k*I_ERI_S_S_Fy2z_S_M4_vrr;
      Double I_ERI_Py_S_Gxy2z_S_M3_vrr = PAY*I_ERI_S_S_Gxy2z_S_M3_vrr+WPY*I_ERI_S_S_Gxy2z_S_M4_vrr+oned2k*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_Pz_S_Gxy2z_S_M3_vrr = PAZ*I_ERI_S_S_Gxy2z_S_M3_vrr+WPZ*I_ERI_S_S_Gxy2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M4_vrr;
      Double I_ERI_Px_S_Gx3z_S_M3_vrr = PAX*I_ERI_S_S_Gx3z_S_M3_vrr+WPX*I_ERI_S_S_Gx3z_S_M4_vrr+oned2k*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Py_S_Gx3z_S_M3_vrr = PAY*I_ERI_S_S_Gx3z_S_M3_vrr+WPY*I_ERI_S_S_Gx3z_S_M4_vrr;
      Double I_ERI_Pz_S_Gx3z_S_M3_vrr = PAZ*I_ERI_S_S_Gx3z_S_M3_vrr+WPZ*I_ERI_S_S_Gx3z_S_M4_vrr+3*oned2k*I_ERI_S_S_Fx2z_S_M4_vrr;
      Double I_ERI_Px_S_G4y_S_M3_vrr = PAX*I_ERI_S_S_G4y_S_M3_vrr+WPX*I_ERI_S_S_G4y_S_M4_vrr;
      Double I_ERI_Py_S_G4y_S_M3_vrr = PAY*I_ERI_S_S_G4y_S_M3_vrr+WPY*I_ERI_S_S_G4y_S_M4_vrr+4*oned2k*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Pz_S_G4y_S_M3_vrr = PAZ*I_ERI_S_S_G4y_S_M3_vrr+WPZ*I_ERI_S_S_G4y_S_M4_vrr;
      Double I_ERI_Px_S_G3yz_S_M3_vrr = PAX*I_ERI_S_S_G3yz_S_M3_vrr+WPX*I_ERI_S_S_G3yz_S_M4_vrr;
      Double I_ERI_Py_S_G3yz_S_M3_vrr = PAY*I_ERI_S_S_G3yz_S_M3_vrr+WPY*I_ERI_S_S_G3yz_S_M4_vrr+3*oned2k*I_ERI_S_S_F2yz_S_M4_vrr;
      Double I_ERI_Pz_S_G3yz_S_M3_vrr = PAZ*I_ERI_S_S_G3yz_S_M3_vrr+WPZ*I_ERI_S_S_G3yz_S_M4_vrr+oned2k*I_ERI_S_S_F3y_S_M4_vrr;
      Double I_ERI_Px_S_G2y2z_S_M3_vrr = PAX*I_ERI_S_S_G2y2z_S_M3_vrr+WPX*I_ERI_S_S_G2y2z_S_M4_vrr;
      Double I_ERI_Py_S_G2y2z_S_M3_vrr = PAY*I_ERI_S_S_G2y2z_S_M3_vrr+WPY*I_ERI_S_S_G2y2z_S_M4_vrr+2*oned2k*I_ERI_S_S_Fy2z_S_M4_vrr;
      Double I_ERI_Pz_S_G2y2z_S_M3_vrr = PAZ*I_ERI_S_S_G2y2z_S_M3_vrr+WPZ*I_ERI_S_S_G2y2z_S_M4_vrr+2*oned2k*I_ERI_S_S_F2yz_S_M4_vrr;
      Double I_ERI_Px_S_Gy3z_S_M3_vrr = PAX*I_ERI_S_S_Gy3z_S_M3_vrr+WPX*I_ERI_S_S_Gy3z_S_M4_vrr;
      Double I_ERI_Py_S_Gy3z_S_M3_vrr = PAY*I_ERI_S_S_Gy3z_S_M3_vrr+WPY*I_ERI_S_S_Gy3z_S_M4_vrr+oned2k*I_ERI_S_S_F3z_S_M4_vrr;
      Double I_ERI_Pz_S_Gy3z_S_M3_vrr = PAZ*I_ERI_S_S_Gy3z_S_M3_vrr+WPZ*I_ERI_S_S_Gy3z_S_M4_vrr+3*oned2k*I_ERI_S_S_Fy2z_S_M4_vrr;
      Double I_ERI_Px_S_G4z_S_M3_vrr = PAX*I_ERI_S_S_G4z_S_M3_vrr+WPX*I_ERI_S_S_G4z_S_M4_vrr;
      Double I_ERI_Py_S_G4z_S_M3_vrr = PAY*I_ERI_S_S_G4z_S_M3_vrr+WPY*I_ERI_S_S_G4z_S_M4_vrr;
      Double I_ERI_Pz_S_G4z_S_M3_vrr = PAZ*I_ERI_S_S_G4z_S_M3_vrr+WPZ*I_ERI_S_S_G4z_S_M4_vrr+4*oned2k*I_ERI_S_S_F3z_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_D_S_P_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M2_vrr = PAX*I_ERI_Px_S_Px_S_M2_vrr+WPX*I_ERI_Px_S_Px_S_M3_vrr+oned2z*I_ERI_S_S_Px_S_M2_vrr-rhod2zsq*I_ERI_S_S_Px_S_M3_vrr+oned2k*I_ERI_Px_S_S_S_M3_vrr;
      Double I_ERI_D2y_S_Px_S_M2_vrr = PAY*I_ERI_Py_S_Px_S_M2_vrr+WPY*I_ERI_Py_S_Px_S_M3_vrr+oned2z*I_ERI_S_S_Px_S_M2_vrr-rhod2zsq*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_D2z_S_Px_S_M2_vrr = PAZ*I_ERI_Pz_S_Px_S_M2_vrr+WPZ*I_ERI_Pz_S_Px_S_M3_vrr+oned2z*I_ERI_S_S_Px_S_M2_vrr-rhod2zsq*I_ERI_S_S_Px_S_M3_vrr;
      Double I_ERI_D2x_S_Py_S_M2_vrr = PAX*I_ERI_Px_S_Py_S_M2_vrr+WPX*I_ERI_Px_S_Py_S_M3_vrr+oned2z*I_ERI_S_S_Py_S_M2_vrr-rhod2zsq*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_D2y_S_Py_S_M2_vrr = PAY*I_ERI_Py_S_Py_S_M2_vrr+WPY*I_ERI_Py_S_Py_S_M3_vrr+oned2z*I_ERI_S_S_Py_S_M2_vrr-rhod2zsq*I_ERI_S_S_Py_S_M3_vrr+oned2k*I_ERI_Py_S_S_S_M3_vrr;
      Double I_ERI_D2z_S_Py_S_M2_vrr = PAZ*I_ERI_Pz_S_Py_S_M2_vrr+WPZ*I_ERI_Pz_S_Py_S_M3_vrr+oned2z*I_ERI_S_S_Py_S_M2_vrr-rhod2zsq*I_ERI_S_S_Py_S_M3_vrr;
      Double I_ERI_D2x_S_Pz_S_M2_vrr = PAX*I_ERI_Px_S_Pz_S_M2_vrr+WPX*I_ERI_Px_S_Pz_S_M3_vrr+oned2z*I_ERI_S_S_Pz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_D2y_S_Pz_S_M2_vrr = PAY*I_ERI_Py_S_Pz_S_M2_vrr+WPY*I_ERI_Py_S_Pz_S_M3_vrr+oned2z*I_ERI_S_S_Pz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M3_vrr;
      Double I_ERI_D2z_S_Pz_S_M2_vrr = PAZ*I_ERI_Pz_S_Pz_S_M2_vrr+WPZ*I_ERI_Pz_S_Pz_S_M3_vrr+oned2z*I_ERI_S_S_Pz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M3_vrr+oned2k*I_ERI_Pz_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_P_S_G_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M3
       ************************************************************/
      Double I_ERI_Px_S_G4x_S_M2_vrr = PAX*I_ERI_S_S_G4x_S_M2_vrr+WPX*I_ERI_S_S_G4x_S_M3_vrr+4*oned2k*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Py_S_G4x_S_M2_vrr = PAY*I_ERI_S_S_G4x_S_M2_vrr+WPY*I_ERI_S_S_G4x_S_M3_vrr;
      Double I_ERI_Pz_S_G4x_S_M2_vrr = PAZ*I_ERI_S_S_G4x_S_M2_vrr+WPZ*I_ERI_S_S_G4x_S_M3_vrr;
      Double I_ERI_Px_S_G3xy_S_M2_vrr = PAX*I_ERI_S_S_G3xy_S_M2_vrr+WPX*I_ERI_S_S_G3xy_S_M3_vrr+3*oned2k*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Py_S_G3xy_S_M2_vrr = PAY*I_ERI_S_S_G3xy_S_M2_vrr+WPY*I_ERI_S_S_G3xy_S_M3_vrr+oned2k*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Pz_S_G3xy_S_M2_vrr = PAZ*I_ERI_S_S_G3xy_S_M2_vrr+WPZ*I_ERI_S_S_G3xy_S_M3_vrr;
      Double I_ERI_Px_S_G3xz_S_M2_vrr = PAX*I_ERI_S_S_G3xz_S_M2_vrr+WPX*I_ERI_S_S_G3xz_S_M3_vrr+3*oned2k*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Py_S_G3xz_S_M2_vrr = PAY*I_ERI_S_S_G3xz_S_M2_vrr+WPY*I_ERI_S_S_G3xz_S_M3_vrr;
      Double I_ERI_Pz_S_G3xz_S_M2_vrr = PAZ*I_ERI_S_S_G3xz_S_M2_vrr+WPZ*I_ERI_S_S_G3xz_S_M3_vrr+oned2k*I_ERI_S_S_F3x_S_M3_vrr;
      Double I_ERI_Px_S_G2x2y_S_M2_vrr = PAX*I_ERI_S_S_G2x2y_S_M2_vrr+WPX*I_ERI_S_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Py_S_G2x2y_S_M2_vrr = PAY*I_ERI_S_S_G2x2y_S_M2_vrr+WPY*I_ERI_S_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Pz_S_G2x2y_S_M2_vrr = PAZ*I_ERI_S_S_G2x2y_S_M2_vrr+WPZ*I_ERI_S_S_G2x2y_S_M3_vrr;
      Double I_ERI_Px_S_G2xyz_S_M2_vrr = PAX*I_ERI_S_S_G2xyz_S_M2_vrr+WPX*I_ERI_S_S_G2xyz_S_M3_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M3_vrr;
      Double I_ERI_Py_S_G2xyz_S_M2_vrr = PAY*I_ERI_S_S_G2xyz_S_M2_vrr+WPY*I_ERI_S_S_G2xyz_S_M3_vrr+oned2k*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Pz_S_G2xyz_S_M2_vrr = PAZ*I_ERI_S_S_G2xyz_S_M2_vrr+WPZ*I_ERI_S_S_G2xyz_S_M3_vrr+oned2k*I_ERI_S_S_F2xy_S_M3_vrr;
      Double I_ERI_Px_S_G2x2z_S_M2_vrr = PAX*I_ERI_S_S_G2x2z_S_M2_vrr+WPX*I_ERI_S_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Py_S_G2x2z_S_M2_vrr = PAY*I_ERI_S_S_G2x2z_S_M2_vrr+WPY*I_ERI_S_S_G2x2z_S_M3_vrr;
      Double I_ERI_Pz_S_G2x2z_S_M2_vrr = PAZ*I_ERI_S_S_G2x2z_S_M2_vrr+WPZ*I_ERI_S_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_S_S_F2xz_S_M3_vrr;
      Double I_ERI_Px_S_Gx3y_S_M2_vrr = PAX*I_ERI_S_S_Gx3y_S_M2_vrr+WPX*I_ERI_S_S_Gx3y_S_M3_vrr+oned2k*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Py_S_Gx3y_S_M2_vrr = PAY*I_ERI_S_S_Gx3y_S_M2_vrr+WPY*I_ERI_S_S_Gx3y_S_M3_vrr+3*oned2k*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Pz_S_Gx3y_S_M2_vrr = PAZ*I_ERI_S_S_Gx3y_S_M2_vrr+WPZ*I_ERI_S_S_Gx3y_S_M3_vrr;
      Double I_ERI_Px_S_Gx2yz_S_M2_vrr = PAX*I_ERI_S_S_Gx2yz_S_M2_vrr+WPX*I_ERI_S_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Py_S_Gx2yz_S_M2_vrr = PAY*I_ERI_S_S_Gx2yz_S_M2_vrr+WPY*I_ERI_S_S_Gx2yz_S_M3_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M3_vrr;
      Double I_ERI_Pz_S_Gx2yz_S_M2_vrr = PAZ*I_ERI_S_S_Gx2yz_S_M2_vrr+WPZ*I_ERI_S_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_S_S_Fx2y_S_M3_vrr;
      Double I_ERI_Px_S_Gxy2z_S_M2_vrr = PAX*I_ERI_S_S_Gxy2z_S_M2_vrr+WPX*I_ERI_S_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_Py_S_Gxy2z_S_M2_vrr = PAY*I_ERI_S_S_Gxy2z_S_M2_vrr+WPY*I_ERI_S_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Pz_S_Gxy2z_S_M2_vrr = PAZ*I_ERI_S_S_Gxy2z_S_M2_vrr+WPZ*I_ERI_S_S_Gxy2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Fxyz_S_M3_vrr;
      Double I_ERI_Px_S_Gx3z_S_M2_vrr = PAX*I_ERI_S_S_Gx3z_S_M2_vrr+WPX*I_ERI_S_S_Gx3z_S_M3_vrr+oned2k*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Py_S_Gx3z_S_M2_vrr = PAY*I_ERI_S_S_Gx3z_S_M2_vrr+WPY*I_ERI_S_S_Gx3z_S_M3_vrr;
      Double I_ERI_Pz_S_Gx3z_S_M2_vrr = PAZ*I_ERI_S_S_Gx3z_S_M2_vrr+WPZ*I_ERI_S_S_Gx3z_S_M3_vrr+3*oned2k*I_ERI_S_S_Fx2z_S_M3_vrr;
      Double I_ERI_Px_S_G4y_S_M2_vrr = PAX*I_ERI_S_S_G4y_S_M2_vrr+WPX*I_ERI_S_S_G4y_S_M3_vrr;
      Double I_ERI_Py_S_G4y_S_M2_vrr = PAY*I_ERI_S_S_G4y_S_M2_vrr+WPY*I_ERI_S_S_G4y_S_M3_vrr+4*oned2k*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Pz_S_G4y_S_M2_vrr = PAZ*I_ERI_S_S_G4y_S_M2_vrr+WPZ*I_ERI_S_S_G4y_S_M3_vrr;
      Double I_ERI_Px_S_G3yz_S_M2_vrr = PAX*I_ERI_S_S_G3yz_S_M2_vrr+WPX*I_ERI_S_S_G3yz_S_M3_vrr;
      Double I_ERI_Py_S_G3yz_S_M2_vrr = PAY*I_ERI_S_S_G3yz_S_M2_vrr+WPY*I_ERI_S_S_G3yz_S_M3_vrr+3*oned2k*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Pz_S_G3yz_S_M2_vrr = PAZ*I_ERI_S_S_G3yz_S_M2_vrr+WPZ*I_ERI_S_S_G3yz_S_M3_vrr+oned2k*I_ERI_S_S_F3y_S_M3_vrr;
      Double I_ERI_Px_S_G2y2z_S_M2_vrr = PAX*I_ERI_S_S_G2y2z_S_M2_vrr+WPX*I_ERI_S_S_G2y2z_S_M3_vrr;
      Double I_ERI_Py_S_G2y2z_S_M2_vrr = PAY*I_ERI_S_S_G2y2z_S_M2_vrr+WPY*I_ERI_S_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_Pz_S_G2y2z_S_M2_vrr = PAZ*I_ERI_S_S_G2y2z_S_M2_vrr+WPZ*I_ERI_S_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_S_S_F2yz_S_M3_vrr;
      Double I_ERI_Px_S_Gy3z_S_M2_vrr = PAX*I_ERI_S_S_Gy3z_S_M2_vrr+WPX*I_ERI_S_S_Gy3z_S_M3_vrr;
      Double I_ERI_Py_S_Gy3z_S_M2_vrr = PAY*I_ERI_S_S_Gy3z_S_M2_vrr+WPY*I_ERI_S_S_Gy3z_S_M3_vrr+oned2k*I_ERI_S_S_F3z_S_M3_vrr;
      Double I_ERI_Pz_S_Gy3z_S_M2_vrr = PAZ*I_ERI_S_S_Gy3z_S_M2_vrr+WPZ*I_ERI_S_S_Gy3z_S_M3_vrr+3*oned2k*I_ERI_S_S_Fy2z_S_M3_vrr;
      Double I_ERI_Px_S_G4z_S_M2_vrr = PAX*I_ERI_S_S_G4z_S_M2_vrr+WPX*I_ERI_S_S_G4z_S_M3_vrr;
      Double I_ERI_Py_S_G4z_S_M2_vrr = PAY*I_ERI_S_S_G4z_S_M2_vrr+WPY*I_ERI_S_S_G4z_S_M3_vrr;
      Double I_ERI_Pz_S_G4z_S_M2_vrr = PAZ*I_ERI_S_S_G4z_S_M2_vrr+WPZ*I_ERI_S_S_G4z_S_M3_vrr+4*oned2k*I_ERI_S_S_F3z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 45 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M3
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_M2_vrr = PAX*I_ERI_Px_S_G4x_S_M2_vrr+WPX*I_ERI_Px_S_G4x_S_M3_vrr+oned2z*I_ERI_S_S_G4x_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M3_vrr+4*oned2k*I_ERI_Px_S_F3x_S_M3_vrr;
      Double I_ERI_D2y_S_G4x_S_M2_vrr = PAY*I_ERI_Py_S_G4x_S_M2_vrr+WPY*I_ERI_Py_S_G4x_S_M3_vrr+oned2z*I_ERI_S_S_G4x_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M3_vrr;
      Double I_ERI_D2z_S_G4x_S_M2_vrr = PAZ*I_ERI_Pz_S_G4x_S_M2_vrr+WPZ*I_ERI_Pz_S_G4x_S_M3_vrr+oned2z*I_ERI_S_S_G4x_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M3_vrr;
      Double I_ERI_D2x_S_G3xy_S_M2_vrr = PAX*I_ERI_Px_S_G3xy_S_M2_vrr+WPX*I_ERI_Px_S_G3xy_S_M3_vrr+oned2z*I_ERI_S_S_G3xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M3_vrr+3*oned2k*I_ERI_Px_S_F2xy_S_M3_vrr;
      Double I_ERI_D2y_S_G3xy_S_M2_vrr = PAY*I_ERI_Py_S_G3xy_S_M2_vrr+WPY*I_ERI_Py_S_G3xy_S_M3_vrr+oned2z*I_ERI_S_S_G3xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M3_vrr+oned2k*I_ERI_Py_S_F3x_S_M3_vrr;
      Double I_ERI_D2z_S_G3xy_S_M2_vrr = PAZ*I_ERI_Pz_S_G3xy_S_M2_vrr+WPZ*I_ERI_Pz_S_G3xy_S_M3_vrr+oned2z*I_ERI_S_S_G3xy_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M3_vrr;
      Double I_ERI_D2x_S_G3xz_S_M2_vrr = PAX*I_ERI_Px_S_G3xz_S_M2_vrr+WPX*I_ERI_Px_S_G3xz_S_M3_vrr+oned2z*I_ERI_S_S_G3xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M3_vrr+3*oned2k*I_ERI_Px_S_F2xz_S_M3_vrr;
      Double I_ERI_D2y_S_G3xz_S_M2_vrr = PAY*I_ERI_Py_S_G3xz_S_M2_vrr+WPY*I_ERI_Py_S_G3xz_S_M3_vrr+oned2z*I_ERI_S_S_G3xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M3_vrr;
      Double I_ERI_D2z_S_G3xz_S_M2_vrr = PAZ*I_ERI_Pz_S_G3xz_S_M2_vrr+WPZ*I_ERI_Pz_S_G3xz_S_M3_vrr+oned2z*I_ERI_S_S_G3xz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M3_vrr+oned2k*I_ERI_Pz_S_F3x_S_M3_vrr;
      Double I_ERI_D2x_S_G2x2y_S_M2_vrr = PAX*I_ERI_Px_S_G2x2y_S_M2_vrr+WPX*I_ERI_Px_S_G2x2y_S_M3_vrr+oned2z*I_ERI_S_S_G2x2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_Px_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2y_S_G2x2y_S_M2_vrr = PAY*I_ERI_Py_S_G2x2y_S_M2_vrr+WPY*I_ERI_Py_S_G2x2y_S_M3_vrr+oned2z*I_ERI_S_S_G2x2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M3_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M3_vrr;
      Double I_ERI_D2z_S_G2x2y_S_M2_vrr = PAZ*I_ERI_Pz_S_G2x2y_S_M2_vrr+WPZ*I_ERI_Pz_S_G2x2y_S_M3_vrr+oned2z*I_ERI_S_S_G2x2y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M3_vrr;
      Double I_ERI_D2x_S_G2xyz_S_M2_vrr = PAX*I_ERI_Px_S_G2xyz_S_M2_vrr+WPX*I_ERI_Px_S_G2xyz_S_M3_vrr+oned2z*I_ERI_S_S_G2xyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M3_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M3_vrr;
      Double I_ERI_D2y_S_G2xyz_S_M2_vrr = PAY*I_ERI_Py_S_G2xyz_S_M2_vrr+WPY*I_ERI_Py_S_G2xyz_S_M3_vrr+oned2z*I_ERI_S_S_G2xyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M3_vrr+oned2k*I_ERI_Py_S_F2xz_S_M3_vrr;
      Double I_ERI_D2z_S_G2xyz_S_M2_vrr = PAZ*I_ERI_Pz_S_G2xyz_S_M2_vrr+WPZ*I_ERI_Pz_S_G2xyz_S_M3_vrr+oned2z*I_ERI_S_S_G2xyz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M3_vrr+oned2k*I_ERI_Pz_S_F2xy_S_M3_vrr;
      Double I_ERI_D2x_S_G2x2z_S_M2_vrr = PAX*I_ERI_Px_S_G2x2z_S_M2_vrr+WPX*I_ERI_Px_S_G2x2z_S_M3_vrr+oned2z*I_ERI_S_S_G2x2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_Px_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2y_S_G2x2z_S_M2_vrr = PAY*I_ERI_Py_S_G2x2z_S_M2_vrr+WPY*I_ERI_Py_S_G2x2z_S_M3_vrr+oned2z*I_ERI_S_S_G2x2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M3_vrr;
      Double I_ERI_D2z_S_G2x2z_S_M2_vrr = PAZ*I_ERI_Pz_S_G2x2z_S_M2_vrr+WPZ*I_ERI_Pz_S_G2x2z_S_M3_vrr+oned2z*I_ERI_S_S_G2x2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M3_vrr;
      Double I_ERI_D2x_S_Gx3y_S_M2_vrr = PAX*I_ERI_Px_S_Gx3y_S_M2_vrr+WPX*I_ERI_Px_S_Gx3y_S_M3_vrr+oned2z*I_ERI_S_S_Gx3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M3_vrr+oned2k*I_ERI_Px_S_F3y_S_M3_vrr;
      Double I_ERI_D2y_S_Gx3y_S_M2_vrr = PAY*I_ERI_Py_S_Gx3y_S_M2_vrr+WPY*I_ERI_Py_S_Gx3y_S_M3_vrr+oned2z*I_ERI_S_S_Gx3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M3_vrr+3*oned2k*I_ERI_Py_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2z_S_Gx3y_S_M2_vrr = PAZ*I_ERI_Pz_S_Gx3y_S_M2_vrr+WPZ*I_ERI_Pz_S_Gx3y_S_M3_vrr+oned2z*I_ERI_S_S_Gx3y_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M3_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_M2_vrr = PAX*I_ERI_Px_S_Gx2yz_S_M2_vrr+WPX*I_ERI_Px_S_Gx2yz_S_M3_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_Px_S_F2yz_S_M3_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_M2_vrr = PAY*I_ERI_Py_S_Gx2yz_S_M2_vrr+WPY*I_ERI_Py_S_Gx2yz_S_M3_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M3_vrr+2*oned2k*I_ERI_Py_S_Fxyz_S_M3_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_M2_vrr = PAZ*I_ERI_Pz_S_Gx2yz_S_M2_vrr+WPZ*I_ERI_Pz_S_Gx2yz_S_M3_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M3_vrr+oned2k*I_ERI_Pz_S_Fx2y_S_M3_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_M2_vrr = PAX*I_ERI_Px_S_Gxy2z_S_M2_vrr+WPX*I_ERI_Px_S_Gxy2z_S_M3_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_Px_S_Fy2z_S_M3_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_M2_vrr = PAY*I_ERI_Py_S_Gxy2z_S_M2_vrr+WPY*I_ERI_Py_S_Gxy2z_S_M3_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M3_vrr+oned2k*I_ERI_Py_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_M2_vrr = PAZ*I_ERI_Pz_S_Gxy2z_S_M2_vrr+WPZ*I_ERI_Pz_S_Gxy2z_S_M3_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_Fxyz_S_M3_vrr;
      Double I_ERI_D2x_S_Gx3z_S_M2_vrr = PAX*I_ERI_Px_S_Gx3z_S_M2_vrr+WPX*I_ERI_Px_S_Gx3z_S_M3_vrr+oned2z*I_ERI_S_S_Gx3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M3_vrr+oned2k*I_ERI_Px_S_F3z_S_M3_vrr;
      Double I_ERI_D2y_S_Gx3z_S_M2_vrr = PAY*I_ERI_Py_S_Gx3z_S_M2_vrr+WPY*I_ERI_Py_S_Gx3z_S_M3_vrr+oned2z*I_ERI_S_S_Gx3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M3_vrr;
      Double I_ERI_D2z_S_Gx3z_S_M2_vrr = PAZ*I_ERI_Pz_S_Gx3z_S_M2_vrr+WPZ*I_ERI_Pz_S_Gx3z_S_M3_vrr+oned2z*I_ERI_S_S_Gx3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M3_vrr+3*oned2k*I_ERI_Pz_S_Fx2z_S_M3_vrr;
      Double I_ERI_D2x_S_G4y_S_M2_vrr = PAX*I_ERI_Px_S_G4y_S_M2_vrr+WPX*I_ERI_Px_S_G4y_S_M3_vrr+oned2z*I_ERI_S_S_G4y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M3_vrr;
      Double I_ERI_D2y_S_G4y_S_M2_vrr = PAY*I_ERI_Py_S_G4y_S_M2_vrr+WPY*I_ERI_Py_S_G4y_S_M3_vrr+oned2z*I_ERI_S_S_G4y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M3_vrr+4*oned2k*I_ERI_Py_S_F3y_S_M3_vrr;
      Double I_ERI_D2z_S_G4y_S_M2_vrr = PAZ*I_ERI_Pz_S_G4y_S_M2_vrr+WPZ*I_ERI_Pz_S_G4y_S_M3_vrr+oned2z*I_ERI_S_S_G4y_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M3_vrr;
      Double I_ERI_D2x_S_G3yz_S_M2_vrr = PAX*I_ERI_Px_S_G3yz_S_M2_vrr+WPX*I_ERI_Px_S_G3yz_S_M3_vrr+oned2z*I_ERI_S_S_G3yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M3_vrr;
      Double I_ERI_D2y_S_G3yz_S_M2_vrr = PAY*I_ERI_Py_S_G3yz_S_M2_vrr+WPY*I_ERI_Py_S_G3yz_S_M3_vrr+oned2z*I_ERI_S_S_G3yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M3_vrr+3*oned2k*I_ERI_Py_S_F2yz_S_M3_vrr;
      Double I_ERI_D2z_S_G3yz_S_M2_vrr = PAZ*I_ERI_Pz_S_G3yz_S_M2_vrr+WPZ*I_ERI_Pz_S_G3yz_S_M3_vrr+oned2z*I_ERI_S_S_G3yz_S_M2_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M3_vrr+oned2k*I_ERI_Pz_S_F3y_S_M3_vrr;
      Double I_ERI_D2x_S_G2y2z_S_M2_vrr = PAX*I_ERI_Px_S_G2y2z_S_M2_vrr+WPX*I_ERI_Px_S_G2y2z_S_M3_vrr+oned2z*I_ERI_S_S_G2y2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M3_vrr;
      Double I_ERI_D2y_S_G2y2z_S_M2_vrr = PAY*I_ERI_Py_S_G2y2z_S_M2_vrr+WPY*I_ERI_Py_S_G2y2z_S_M3_vrr+oned2z*I_ERI_S_S_G2y2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_Py_S_Fy2z_S_M3_vrr;
      Double I_ERI_D2z_S_G2y2z_S_M2_vrr = PAZ*I_ERI_Pz_S_G2y2z_S_M2_vrr+WPZ*I_ERI_Pz_S_G2y2z_S_M3_vrr+oned2z*I_ERI_S_S_G2y2z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M3_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M3_vrr;
      Double I_ERI_D2x_S_Gy3z_S_M2_vrr = PAX*I_ERI_Px_S_Gy3z_S_M2_vrr+WPX*I_ERI_Px_S_Gy3z_S_M3_vrr+oned2z*I_ERI_S_S_Gy3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M3_vrr;
      Double I_ERI_D2y_S_Gy3z_S_M2_vrr = PAY*I_ERI_Py_S_Gy3z_S_M2_vrr+WPY*I_ERI_Py_S_Gy3z_S_M3_vrr+oned2z*I_ERI_S_S_Gy3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M3_vrr+oned2k*I_ERI_Py_S_F3z_S_M3_vrr;
      Double I_ERI_D2z_S_Gy3z_S_M2_vrr = PAZ*I_ERI_Pz_S_Gy3z_S_M2_vrr+WPZ*I_ERI_Pz_S_Gy3z_S_M3_vrr+oned2z*I_ERI_S_S_Gy3z_S_M2_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M3_vrr+3*oned2k*I_ERI_Pz_S_Fy2z_S_M3_vrr;
      Double I_ERI_D2x_S_G4z_S_M2_vrr = PAX*I_ERI_Px_S_G4z_S_M2_vrr+WPX*I_ERI_Px_S_G4z_S_M3_vrr+oned2z*I_ERI_S_S_G4z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M3_vrr;
      Double I_ERI_D2y_S_G4z_S_M2_vrr = PAY*I_ERI_Py_S_G4z_S_M2_vrr+WPY*I_ERI_Py_S_G4z_S_M3_vrr+oned2z*I_ERI_S_S_G4z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M3_vrr;
      Double I_ERI_D2z_S_G4z_S_M2_vrr = PAZ*I_ERI_Pz_S_G4z_S_M2_vrr+WPZ*I_ERI_Pz_S_G4z_S_M3_vrr+oned2z*I_ERI_S_S_G4z_S_M2_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M3_vrr+4*oned2k*I_ERI_Pz_S_F3z_S_M3_vrr;

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
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M1_vrr = QCX*I_ERI_S_S_Px_S_M1_vrr+WQX*I_ERI_S_S_Px_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Dxy_S_M1_vrr = QCY*I_ERI_S_S_Px_S_M1_vrr+WQY*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_S_Dxz_S_M1_vrr = QCZ*I_ERI_S_S_Px_S_M1_vrr+WQZ*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_S_S_D2y_S_M1_vrr = QCY*I_ERI_S_S_Py_S_M1_vrr+WQY*I_ERI_S_S_Py_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_S_S_Dyz_S_M1_vrr = QCZ*I_ERI_S_S_Py_S_M1_vrr+WQZ*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_S_S_D2z_S_M1_vrr = QCZ*I_ERI_S_S_Pz_S_M1_vrr+WQZ*I_ERI_S_S_Pz_S_M2_vrr+oned2e*I_ERI_S_S_S_S_M1_vrr-rhod2esq*I_ERI_S_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_M1_vrr = PAX*I_ERI_S_S_D2x_S_M1_vrr+WPX*I_ERI_S_S_D2x_S_M2_vrr+2*oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Py_S_D2x_S_M1_vrr = PAY*I_ERI_S_S_D2x_S_M1_vrr+WPY*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Pz_S_D2x_S_M1_vrr = PAZ*I_ERI_S_S_D2x_S_M1_vrr+WPZ*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_Px_S_Dxy_S_M1_vrr = PAX*I_ERI_S_S_Dxy_S_M1_vrr+WPX*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Py_S_Dxy_S_M1_vrr = PAY*I_ERI_S_S_Dxy_S_M1_vrr+WPY*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Pz_S_Dxy_S_M1_vrr = PAZ*I_ERI_S_S_Dxy_S_M1_vrr+WPZ*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_Px_S_Dxz_S_M1_vrr = PAX*I_ERI_S_S_Dxz_S_M1_vrr+WPX*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Py_S_Dxz_S_M1_vrr = PAY*I_ERI_S_S_Dxz_S_M1_vrr+WPY*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_Pz_S_Dxz_S_M1_vrr = PAZ*I_ERI_S_S_Dxz_S_M1_vrr+WPZ*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_Px_S_D2y_S_M1_vrr = PAX*I_ERI_S_S_D2y_S_M1_vrr+WPX*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Py_S_D2y_S_M1_vrr = PAY*I_ERI_S_S_D2y_S_M1_vrr+WPY*I_ERI_S_S_D2y_S_M2_vrr+2*oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Pz_S_D2y_S_M1_vrr = PAZ*I_ERI_S_S_D2y_S_M1_vrr+WPZ*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Px_S_Dyz_S_M1_vrr = PAX*I_ERI_S_S_Dyz_S_M1_vrr+WPX*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Py_S_Dyz_S_M1_vrr = PAY*I_ERI_S_S_Dyz_S_M1_vrr+WPY*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Pz_S_Dyz_S_M1_vrr = PAZ*I_ERI_S_S_Dyz_S_M1_vrr+WPZ*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Px_S_D2z_S_M1_vrr = PAX*I_ERI_S_S_D2z_S_M1_vrr+WPX*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Py_S_D2z_S_M1_vrr = PAY*I_ERI_S_S_D2z_S_M1_vrr+WPY*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Pz_S_D2z_S_M1_vrr = PAZ*I_ERI_S_S_D2z_S_M1_vrr+WPZ*I_ERI_S_S_D2z_S_M2_vrr+2*oned2k*I_ERI_S_S_Pz_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 12 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_M1_vrr = PAX*I_ERI_Px_S_D2x_S_M1_vrr+WPX*I_ERI_Px_S_D2x_S_M2_vrr+oned2z*I_ERI_S_S_D2x_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M2_vrr+2*oned2k*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_Dxy_S_D2x_S_M1_vrr = PAY*I_ERI_Px_S_D2x_S_M1_vrr+WPY*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_D2y_S_D2x_S_M1_vrr = PAY*I_ERI_Py_S_D2x_S_M1_vrr+WPY*I_ERI_Py_S_D2x_S_M2_vrr+oned2z*I_ERI_S_S_D2x_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_D2z_S_D2x_S_M1_vrr = PAZ*I_ERI_Pz_S_D2x_S_M1_vrr+WPZ*I_ERI_Pz_S_D2x_S_M2_vrr+oned2z*I_ERI_S_S_D2x_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M2_vrr;
      Double I_ERI_D2x_S_Dxy_S_M1_vrr = PAX*I_ERI_Px_S_Dxy_S_M1_vrr+WPX*I_ERI_Px_S_Dxy_S_M2_vrr+oned2z*I_ERI_S_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_Dxy_S_Dxy_S_M1_vrr = PAY*I_ERI_Px_S_Dxy_S_M1_vrr+WPY*I_ERI_Px_S_Dxy_S_M2_vrr+oned2k*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_D2y_S_Dxy_S_M1_vrr = PAY*I_ERI_Py_S_Dxy_S_M1_vrr+WPY*I_ERI_Py_S_Dxy_S_M2_vrr+oned2z*I_ERI_S_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M2_vrr+oned2k*I_ERI_Py_S_Px_S_M2_vrr;
      Double I_ERI_D2z_S_Dxy_S_M1_vrr = PAZ*I_ERI_Pz_S_Dxy_S_M1_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M2_vrr+oned2z*I_ERI_S_S_Dxy_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M2_vrr;
      Double I_ERI_D2x_S_Dxz_S_M1_vrr = PAX*I_ERI_Px_S_Dxz_S_M1_vrr+WPX*I_ERI_Px_S_Dxz_S_M2_vrr+oned2z*I_ERI_S_S_Dxz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_Dxy_S_Dxz_S_M1_vrr = PAY*I_ERI_Px_S_Dxz_S_M1_vrr+WPY*I_ERI_Px_S_Dxz_S_M2_vrr;
      Double I_ERI_D2y_S_Dxz_S_M1_vrr = PAY*I_ERI_Py_S_Dxz_S_M1_vrr+WPY*I_ERI_Py_S_Dxz_S_M2_vrr+oned2z*I_ERI_S_S_Dxz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M2_vrr;
      Double I_ERI_D2z_S_Dxz_S_M1_vrr = PAZ*I_ERI_Pz_S_Dxz_S_M1_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M2_vrr+oned2z*I_ERI_S_S_Dxz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M2_vrr+oned2k*I_ERI_Pz_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_D2y_S_M1_vrr = PAX*I_ERI_Px_S_D2y_S_M1_vrr+WPX*I_ERI_Px_S_D2y_S_M2_vrr+oned2z*I_ERI_S_S_D2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_Dxy_S_D2y_S_M1_vrr = PAY*I_ERI_Px_S_D2y_S_M1_vrr+WPY*I_ERI_Px_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Px_S_Py_S_M2_vrr;
      Double I_ERI_D2y_S_D2y_S_M1_vrr = PAY*I_ERI_Py_S_D2y_S_M1_vrr+WPY*I_ERI_Py_S_D2y_S_M2_vrr+oned2z*I_ERI_S_S_D2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M2_vrr+2*oned2k*I_ERI_Py_S_Py_S_M2_vrr;
      Double I_ERI_D2z_S_D2y_S_M1_vrr = PAZ*I_ERI_Pz_S_D2y_S_M1_vrr+WPZ*I_ERI_Pz_S_D2y_S_M2_vrr+oned2z*I_ERI_S_S_D2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M2_vrr;
      Double I_ERI_D2x_S_Dyz_S_M1_vrr = PAX*I_ERI_Px_S_Dyz_S_M1_vrr+WPX*I_ERI_Px_S_Dyz_S_M2_vrr+oned2z*I_ERI_S_S_Dyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M2_vrr;
      Double I_ERI_Dxy_S_Dyz_S_M1_vrr = PAY*I_ERI_Px_S_Dyz_S_M1_vrr+WPY*I_ERI_Px_S_Dyz_S_M2_vrr+oned2k*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_D2y_S_Dyz_S_M1_vrr = PAY*I_ERI_Py_S_Dyz_S_M1_vrr+WPY*I_ERI_Py_S_Dyz_S_M2_vrr+oned2z*I_ERI_S_S_Dyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_Py_S_Pz_S_M2_vrr;
      Double I_ERI_D2z_S_Dyz_S_M1_vrr = PAZ*I_ERI_Pz_S_Dyz_S_M1_vrr+WPZ*I_ERI_Pz_S_Dyz_S_M2_vrr+oned2z*I_ERI_S_S_Dyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M2_vrr+oned2k*I_ERI_Pz_S_Py_S_M2_vrr;
      Double I_ERI_D2x_S_D2z_S_M1_vrr = PAX*I_ERI_Px_S_D2z_S_M1_vrr+WPX*I_ERI_Px_S_D2z_S_M2_vrr+oned2z*I_ERI_S_S_D2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_Dxy_S_D2z_S_M1_vrr = PAY*I_ERI_Px_S_D2z_S_M1_vrr+WPY*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_D2y_S_D2z_S_M1_vrr = PAY*I_ERI_Py_S_D2z_S_M1_vrr+WPY*I_ERI_Py_S_D2z_S_M2_vrr+oned2z*I_ERI_S_S_D2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M2_vrr;
      Double I_ERI_D2z_S_D2z_S_M1_vrr = PAZ*I_ERI_Pz_S_D2z_S_M1_vrr+WPZ*I_ERI_Pz_S_D2z_S_M2_vrr+oned2z*I_ERI_S_S_D2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M2_vrr;

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
       * totally 20 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_M1_vrr = PAX*I_ERI_Px_S_F3x_S_M1_vrr+WPX*I_ERI_Px_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_Dxy_S_F3x_S_M1_vrr = PAY*I_ERI_Px_S_F3x_S_M1_vrr+WPY*I_ERI_Px_S_F3x_S_M2_vrr;
      Double I_ERI_D2y_S_F3x_S_M1_vrr = PAY*I_ERI_Py_S_F3x_S_M1_vrr+WPY*I_ERI_Py_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_D2z_S_F3x_S_M1_vrr = PAZ*I_ERI_Pz_S_F3x_S_M1_vrr+WPZ*I_ERI_Pz_S_F3x_S_M2_vrr+oned2z*I_ERI_S_S_F3x_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M2_vrr;
      Double I_ERI_D2x_S_F2xy_S_M1_vrr = PAX*I_ERI_Px_S_F2xy_S_M1_vrr+WPX*I_ERI_Px_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M2_vrr;
      Double I_ERI_Dxy_S_F2xy_S_M1_vrr = PAY*I_ERI_Px_S_F2xy_S_M1_vrr+WPY*I_ERI_Px_S_F2xy_S_M2_vrr+oned2k*I_ERI_Px_S_D2x_S_M2_vrr;
      Double I_ERI_D2y_S_F2xy_S_M1_vrr = PAY*I_ERI_Py_S_F2xy_S_M1_vrr+WPY*I_ERI_Py_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr+oned2k*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_D2z_S_F2xy_S_M1_vrr = PAZ*I_ERI_Pz_S_F2xy_S_M1_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M2_vrr+oned2z*I_ERI_S_S_F2xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M2_vrr;
      Double I_ERI_D2x_S_F2xz_S_M1_vrr = PAX*I_ERI_Px_S_F2xz_S_M1_vrr+WPX*I_ERI_Px_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M2_vrr;
      Double I_ERI_Dxy_S_F2xz_S_M1_vrr = PAY*I_ERI_Px_S_F2xz_S_M1_vrr+WPY*I_ERI_Px_S_F2xz_S_M2_vrr;
      Double I_ERI_D2y_S_F2xz_S_M1_vrr = PAY*I_ERI_Py_S_F2xz_S_M1_vrr+WPY*I_ERI_Py_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr;
      Double I_ERI_D2z_S_F2xz_S_M1_vrr = PAZ*I_ERI_Pz_S_F2xz_S_M1_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M2_vrr+oned2z*I_ERI_S_S_F2xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M2_vrr+oned2k*I_ERI_Pz_S_D2x_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2y_S_M1_vrr = PAX*I_ERI_Px_S_Fx2y_S_M1_vrr+WPX*I_ERI_Px_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr+oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_Dxy_S_Fx2y_S_M1_vrr = PAY*I_ERI_Px_S_Fx2y_S_M1_vrr+WPY*I_ERI_Px_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2y_S_M1_vrr = PAY*I_ERI_Py_S_Fx2y_S_M1_vrr+WPY*I_ERI_Py_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2y_S_M1_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_M1_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M2_vrr+oned2z*I_ERI_S_S_Fx2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2x_S_Fxyz_S_M1_vrr = PAX*I_ERI_Px_S_Fxyz_S_M1_vrr+WPX*I_ERI_Px_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Px_S_Dyz_S_M2_vrr;
      Double I_ERI_Dxy_S_Fxyz_S_M1_vrr = PAY*I_ERI_Px_S_Fxyz_S_M1_vrr+WPY*I_ERI_Px_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Px_S_Dxz_S_M2_vrr;
      Double I_ERI_D2y_S_Fxyz_S_M1_vrr = PAY*I_ERI_Py_S_Fxyz_S_M1_vrr+WPY*I_ERI_Py_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Py_S_Dxz_S_M2_vrr;
      Double I_ERI_D2z_S_Fxyz_S_M1_vrr = PAZ*I_ERI_Pz_S_Fxyz_S_M1_vrr+WPZ*I_ERI_Pz_S_Fxyz_S_M2_vrr+oned2z*I_ERI_S_S_Fxyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M2_vrr+oned2k*I_ERI_Pz_S_Dxy_S_M2_vrr;
      Double I_ERI_D2x_S_Fx2z_S_M1_vrr = PAX*I_ERI_Px_S_Fx2z_S_M1_vrr+WPX*I_ERI_Px_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr+oned2k*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_Dxy_S_Fx2z_S_M1_vrr = PAY*I_ERI_Px_S_Fx2z_S_M1_vrr+WPY*I_ERI_Px_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2y_S_Fx2z_S_M1_vrr = PAY*I_ERI_Py_S_Fx2z_S_M1_vrr+WPY*I_ERI_Py_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2z_S_Fx2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M2_vrr+oned2z*I_ERI_S_S_Fx2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M2_vrr;
      Double I_ERI_D2x_S_F3y_S_M1_vrr = PAX*I_ERI_Px_S_F3y_S_M1_vrr+WPX*I_ERI_Px_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_Dxy_S_F3y_S_M1_vrr = PAY*I_ERI_Px_S_F3y_S_M1_vrr+WPY*I_ERI_Px_S_F3y_S_M2_vrr+3*oned2k*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_D2y_S_F3y_S_M1_vrr = PAY*I_ERI_Py_S_F3y_S_M1_vrr+WPY*I_ERI_Py_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M2_vrr;
      Double I_ERI_D2z_S_F3y_S_M1_vrr = PAZ*I_ERI_Pz_S_F3y_S_M1_vrr+WPZ*I_ERI_Pz_S_F3y_S_M2_vrr+oned2z*I_ERI_S_S_F3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M2_vrr;
      Double I_ERI_D2x_S_F2yz_S_M1_vrr = PAX*I_ERI_Px_S_F2yz_S_M1_vrr+WPX*I_ERI_Px_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr;
      Double I_ERI_Dxy_S_F2yz_S_M1_vrr = PAY*I_ERI_Px_S_F2yz_S_M1_vrr+WPY*I_ERI_Px_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Dyz_S_M2_vrr;
      Double I_ERI_D2y_S_F2yz_S_M1_vrr = PAY*I_ERI_Py_S_F2yz_S_M1_vrr+WPY*I_ERI_Py_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M2_vrr;
      Double I_ERI_D2z_S_F2yz_S_M1_vrr = PAZ*I_ERI_Pz_S_F2yz_S_M1_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M2_vrr+oned2z*I_ERI_S_S_F2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M2_vrr+oned2k*I_ERI_Pz_S_D2y_S_M2_vrr;
      Double I_ERI_D2x_S_Fy2z_S_M1_vrr = PAX*I_ERI_Px_S_Fy2z_S_M1_vrr+WPX*I_ERI_Px_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr;
      Double I_ERI_Dxy_S_Fy2z_S_M1_vrr = PAY*I_ERI_Px_S_Fy2z_S_M1_vrr+WPY*I_ERI_Px_S_Fy2z_S_M2_vrr+oned2k*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_D2y_S_Fy2z_S_M1_vrr = PAY*I_ERI_Py_S_Fy2z_S_M1_vrr+WPY*I_ERI_Py_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr+oned2k*I_ERI_Py_S_D2z_S_M2_vrr;
      Double I_ERI_D2z_S_Fy2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Fy2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Fy2z_S_M2_vrr+oned2z*I_ERI_S_S_Fy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Dyz_S_M2_vrr;
      Double I_ERI_D2x_S_F3z_S_M1_vrr = PAX*I_ERI_Px_S_F3z_S_M1_vrr+WPX*I_ERI_Px_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_Dxy_S_F3z_S_M1_vrr = PAY*I_ERI_Px_S_F3z_S_M1_vrr+WPY*I_ERI_Px_S_F3z_S_M2_vrr;
      Double I_ERI_D2y_S_F3z_S_M1_vrr = PAY*I_ERI_Py_S_F3z_S_M1_vrr+WPY*I_ERI_Py_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr;
      Double I_ERI_D2z_S_F3z_S_M1_vrr = PAZ*I_ERI_Pz_S_F3z_S_M1_vrr+WPZ*I_ERI_Pz_S_F3z_S_M2_vrr+oned2z*I_ERI_S_S_F3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 12 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M2
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M2
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_M1_vrr = PAX*I_ERI_D2x_S_D2x_S_M1_vrr+WPX*I_ERI_D2x_S_D2x_S_M2_vrr+2*oned2z*I_ERI_Px_S_D2x_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_D2x_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_F2xy_S_D2x_S_M1_vrr = PAY*I_ERI_D2x_S_D2x_S_M1_vrr+WPY*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_F2xz_S_D2x_S_M1_vrr = PAZ*I_ERI_D2x_S_D2x_S_M1_vrr+WPZ*I_ERI_D2x_S_D2x_S_M2_vrr;
      Double I_ERI_Fx2y_S_D2x_S_M1_vrr = PAX*I_ERI_D2y_S_D2x_S_M1_vrr+WPX*I_ERI_D2y_S_D2x_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_Fx2z_S_D2x_S_M1_vrr = PAX*I_ERI_D2z_S_D2x_S_M1_vrr+WPX*I_ERI_D2z_S_D2x_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_F3y_S_D2x_S_M1_vrr = PAY*I_ERI_D2y_S_D2x_S_M1_vrr+WPY*I_ERI_D2y_S_D2x_S_M2_vrr+2*oned2z*I_ERI_Py_S_D2x_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_D2x_S_M2_vrr;
      Double I_ERI_F2yz_S_D2x_S_M1_vrr = PAZ*I_ERI_D2y_S_D2x_S_M1_vrr+WPZ*I_ERI_D2y_S_D2x_S_M2_vrr;
      Double I_ERI_F3z_S_D2x_S_M1_vrr = PAZ*I_ERI_D2z_S_D2x_S_M1_vrr+WPZ*I_ERI_D2z_S_D2x_S_M2_vrr+2*oned2z*I_ERI_Pz_S_D2x_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_D2x_S_M2_vrr;
      Double I_ERI_F3x_S_Dxy_S_M1_vrr = PAX*I_ERI_D2x_S_Dxy_S_M1_vrr+WPX*I_ERI_D2x_S_Dxy_S_M2_vrr+2*oned2z*I_ERI_Px_S_Dxy_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_F2xy_S_Dxy_S_M1_vrr = PAY*I_ERI_D2x_S_Dxy_S_M1_vrr+WPY*I_ERI_D2x_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_F2xz_S_Dxy_S_M1_vrr = PAZ*I_ERI_D2x_S_Dxy_S_M1_vrr+WPZ*I_ERI_D2x_S_Dxy_S_M2_vrr;
      Double I_ERI_Fx2y_S_Dxy_S_M1_vrr = PAX*I_ERI_D2y_S_Dxy_S_M1_vrr+WPX*I_ERI_D2y_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_Fx2z_S_Dxy_S_M1_vrr = PAX*I_ERI_D2z_S_Dxy_S_M1_vrr+WPX*I_ERI_D2z_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_F3y_S_Dxy_S_M1_vrr = PAY*I_ERI_D2y_S_Dxy_S_M1_vrr+WPY*I_ERI_D2y_S_Dxy_S_M2_vrr+2*oned2z*I_ERI_Py_S_Dxy_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Dxy_S_M2_vrr+oned2k*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_F2yz_S_Dxy_S_M1_vrr = PAZ*I_ERI_D2y_S_Dxy_S_M1_vrr+WPZ*I_ERI_D2y_S_Dxy_S_M2_vrr;
      Double I_ERI_F3z_S_Dxy_S_M1_vrr = PAZ*I_ERI_D2z_S_Dxy_S_M1_vrr+WPZ*I_ERI_D2z_S_Dxy_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Dxy_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxy_S_M2_vrr;
      Double I_ERI_F3x_S_Dxz_S_M1_vrr = PAX*I_ERI_D2x_S_Dxz_S_M1_vrr+WPX*I_ERI_D2x_S_Dxz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Dxz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_F2xy_S_Dxz_S_M1_vrr = PAY*I_ERI_D2x_S_Dxz_S_M1_vrr+WPY*I_ERI_D2x_S_Dxz_S_M2_vrr;
      Double I_ERI_F2xz_S_Dxz_S_M1_vrr = PAZ*I_ERI_D2x_S_Dxz_S_M1_vrr+WPZ*I_ERI_D2x_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2x_S_Px_S_M2_vrr;
      Double I_ERI_Fx2y_S_Dxz_S_M1_vrr = PAX*I_ERI_D2y_S_Dxz_S_M1_vrr+WPX*I_ERI_D2y_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Dxz_S_M1_vrr = PAX*I_ERI_D2z_S_Dxz_S_M1_vrr+WPX*I_ERI_D2z_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2z_S_Pz_S_M2_vrr;
      Double I_ERI_F3y_S_Dxz_S_M1_vrr = PAY*I_ERI_D2y_S_Dxz_S_M1_vrr+WPY*I_ERI_D2y_S_Dxz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Dxz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Dxz_S_M2_vrr;
      Double I_ERI_F2yz_S_Dxz_S_M1_vrr = PAZ*I_ERI_D2y_S_Dxz_S_M1_vrr+WPZ*I_ERI_D2y_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2y_S_Px_S_M2_vrr;
      Double I_ERI_F3z_S_Dxz_S_M1_vrr = PAZ*I_ERI_D2z_S_Dxz_S_M1_vrr+WPZ*I_ERI_D2z_S_Dxz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Dxz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxz_S_M2_vrr+oned2k*I_ERI_D2z_S_Px_S_M2_vrr;
      Double I_ERI_F3x_S_D2y_S_M1_vrr = PAX*I_ERI_D2x_S_D2y_S_M1_vrr+WPX*I_ERI_D2x_S_D2y_S_M2_vrr+2*oned2z*I_ERI_Px_S_D2y_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_D2y_S_M2_vrr;
      Double I_ERI_F2xy_S_D2y_S_M1_vrr = PAY*I_ERI_D2x_S_D2y_S_M1_vrr+WPY*I_ERI_D2x_S_D2y_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_F2xz_S_D2y_S_M1_vrr = PAZ*I_ERI_D2x_S_D2y_S_M1_vrr+WPZ*I_ERI_D2x_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2y_S_D2y_S_M1_vrr = PAX*I_ERI_D2y_S_D2y_S_M1_vrr+WPX*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_Fx2z_S_D2y_S_M1_vrr = PAX*I_ERI_D2z_S_D2y_S_M1_vrr+WPX*I_ERI_D2z_S_D2y_S_M2_vrr;
      Double I_ERI_F3y_S_D2y_S_M1_vrr = PAY*I_ERI_D2y_S_D2y_S_M1_vrr+WPY*I_ERI_D2y_S_D2y_S_M2_vrr+2*oned2z*I_ERI_Py_S_D2y_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_D2y_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_F2yz_S_D2y_S_M1_vrr = PAZ*I_ERI_D2y_S_D2y_S_M1_vrr+WPZ*I_ERI_D2y_S_D2y_S_M2_vrr;
      Double I_ERI_F3z_S_D2y_S_M1_vrr = PAZ*I_ERI_D2z_S_D2y_S_M1_vrr+WPZ*I_ERI_D2z_S_D2y_S_M2_vrr+2*oned2z*I_ERI_Pz_S_D2y_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_D2y_S_M2_vrr;
      Double I_ERI_F3x_S_Dyz_S_M1_vrr = PAX*I_ERI_D2x_S_Dyz_S_M1_vrr+WPX*I_ERI_D2x_S_Dyz_S_M2_vrr+2*oned2z*I_ERI_Px_S_Dyz_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_Dyz_S_M2_vrr;
      Double I_ERI_F2xy_S_Dyz_S_M1_vrr = PAY*I_ERI_D2x_S_Dyz_S_M1_vrr+WPY*I_ERI_D2x_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_F2xz_S_Dyz_S_M1_vrr = PAZ*I_ERI_D2x_S_Dyz_S_M1_vrr+WPZ*I_ERI_D2x_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2x_S_Py_S_M2_vrr;
      Double I_ERI_Fx2y_S_Dyz_S_M1_vrr = PAX*I_ERI_D2y_S_Dyz_S_M1_vrr+WPX*I_ERI_D2y_S_Dyz_S_M2_vrr;
      Double I_ERI_Fx2z_S_Dyz_S_M1_vrr = PAX*I_ERI_D2z_S_Dyz_S_M1_vrr+WPX*I_ERI_D2z_S_Dyz_S_M2_vrr;
      Double I_ERI_F3y_S_Dyz_S_M1_vrr = PAY*I_ERI_D2y_S_Dyz_S_M1_vrr+WPY*I_ERI_D2y_S_Dyz_S_M2_vrr+2*oned2z*I_ERI_Py_S_Dyz_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_F2yz_S_Dyz_S_M1_vrr = PAZ*I_ERI_D2y_S_Dyz_S_M1_vrr+WPZ*I_ERI_D2y_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2y_S_Py_S_M2_vrr;
      Double I_ERI_F3z_S_Dyz_S_M1_vrr = PAZ*I_ERI_D2z_S_Dyz_S_M1_vrr+WPZ*I_ERI_D2z_S_Dyz_S_M2_vrr+2*oned2z*I_ERI_Pz_S_Dyz_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_Dyz_S_M2_vrr+oned2k*I_ERI_D2z_S_Py_S_M2_vrr;
      Double I_ERI_F3x_S_D2z_S_M1_vrr = PAX*I_ERI_D2x_S_D2z_S_M1_vrr+WPX*I_ERI_D2x_S_D2z_S_M2_vrr+2*oned2z*I_ERI_Px_S_D2z_S_M1_vrr-2*rhod2zsq*I_ERI_Px_S_D2z_S_M2_vrr;
      Double I_ERI_F2xy_S_D2z_S_M1_vrr = PAY*I_ERI_D2x_S_D2z_S_M1_vrr+WPY*I_ERI_D2x_S_D2z_S_M2_vrr;
      Double I_ERI_F2xz_S_D2z_S_M1_vrr = PAZ*I_ERI_D2x_S_D2z_S_M1_vrr+WPZ*I_ERI_D2x_S_D2z_S_M2_vrr+2*oned2k*I_ERI_D2x_S_Pz_S_M2_vrr;
      Double I_ERI_Fx2y_S_D2z_S_M1_vrr = PAX*I_ERI_D2y_S_D2z_S_M1_vrr+WPX*I_ERI_D2y_S_D2z_S_M2_vrr;
      Double I_ERI_Fx2z_S_D2z_S_M1_vrr = PAX*I_ERI_D2z_S_D2z_S_M1_vrr+WPX*I_ERI_D2z_S_D2z_S_M2_vrr;
      Double I_ERI_F3y_S_D2z_S_M1_vrr = PAY*I_ERI_D2y_S_D2z_S_M1_vrr+WPY*I_ERI_D2y_S_D2z_S_M2_vrr+2*oned2z*I_ERI_Py_S_D2z_S_M1_vrr-2*rhod2zsq*I_ERI_Py_S_D2z_S_M2_vrr;
      Double I_ERI_F2yz_S_D2z_S_M1_vrr = PAZ*I_ERI_D2y_S_D2z_S_M1_vrr+WPZ*I_ERI_D2y_S_D2z_S_M2_vrr+2*oned2k*I_ERI_D2y_S_Pz_S_M2_vrr;
      Double I_ERI_F3z_S_D2z_S_M1_vrr = PAZ*I_ERI_D2z_S_D2z_S_M1_vrr+WPZ*I_ERI_D2z_S_D2z_S_M2_vrr+2*oned2z*I_ERI_Pz_S_D2z_S_M1_vrr-2*rhod2zsq*I_ERI_Pz_S_D2z_S_M2_vrr+2*oned2k*I_ERI_D2z_S_Pz_S_M2_vrr;

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
       * totally 30 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_M1_vrr = PAX*I_ERI_Px_S_G4x_S_M1_vrr+WPX*I_ERI_Px_S_G4x_S_M2_vrr+oned2z*I_ERI_S_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M2_vrr+4*oned2k*I_ERI_Px_S_F3x_S_M2_vrr;
      Double I_ERI_Dxy_S_G4x_S_M1_vrr = PAY*I_ERI_Px_S_G4x_S_M1_vrr+WPY*I_ERI_Px_S_G4x_S_M2_vrr;
      Double I_ERI_D2y_S_G4x_S_M1_vrr = PAY*I_ERI_Py_S_G4x_S_M1_vrr+WPY*I_ERI_Py_S_G4x_S_M2_vrr+oned2z*I_ERI_S_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_D2z_S_G4x_S_M1_vrr = PAZ*I_ERI_Pz_S_G4x_S_M1_vrr+WPZ*I_ERI_Pz_S_G4x_S_M2_vrr+oned2z*I_ERI_S_S_G4x_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M2_vrr;
      Double I_ERI_D2x_S_G3xy_S_M1_vrr = PAX*I_ERI_Px_S_G3xy_S_M1_vrr+WPX*I_ERI_Px_S_G3xy_S_M2_vrr+oned2z*I_ERI_S_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M2_vrr+3*oned2k*I_ERI_Px_S_F2xy_S_M2_vrr;
      Double I_ERI_Dxy_S_G3xy_S_M1_vrr = PAY*I_ERI_Px_S_G3xy_S_M1_vrr+WPY*I_ERI_Px_S_G3xy_S_M2_vrr+oned2k*I_ERI_Px_S_F3x_S_M2_vrr;
      Double I_ERI_D2y_S_G3xy_S_M1_vrr = PAY*I_ERI_Py_S_G3xy_S_M1_vrr+WPY*I_ERI_Py_S_G3xy_S_M2_vrr+oned2z*I_ERI_S_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M2_vrr+oned2k*I_ERI_Py_S_F3x_S_M2_vrr;
      Double I_ERI_D2z_S_G3xy_S_M1_vrr = PAZ*I_ERI_Pz_S_G3xy_S_M1_vrr+WPZ*I_ERI_Pz_S_G3xy_S_M2_vrr+oned2z*I_ERI_S_S_G3xy_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M2_vrr;
      Double I_ERI_D2x_S_G3xz_S_M1_vrr = PAX*I_ERI_Px_S_G3xz_S_M1_vrr+WPX*I_ERI_Px_S_G3xz_S_M2_vrr+oned2z*I_ERI_S_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M2_vrr+3*oned2k*I_ERI_Px_S_F2xz_S_M2_vrr;
      Double I_ERI_Dxy_S_G3xz_S_M1_vrr = PAY*I_ERI_Px_S_G3xz_S_M1_vrr+WPY*I_ERI_Px_S_G3xz_S_M2_vrr;
      Double I_ERI_D2y_S_G3xz_S_M1_vrr = PAY*I_ERI_Py_S_G3xz_S_M1_vrr+WPY*I_ERI_Py_S_G3xz_S_M2_vrr+oned2z*I_ERI_S_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M2_vrr;
      Double I_ERI_D2z_S_G3xz_S_M1_vrr = PAZ*I_ERI_Pz_S_G3xz_S_M1_vrr+WPZ*I_ERI_Pz_S_G3xz_S_M2_vrr+oned2z*I_ERI_S_S_G3xz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M2_vrr+oned2k*I_ERI_Pz_S_F3x_S_M2_vrr;
      Double I_ERI_D2x_S_G2x2y_S_M1_vrr = PAX*I_ERI_Px_S_G2x2y_S_M1_vrr+WPX*I_ERI_Px_S_G2x2y_S_M2_vrr+oned2z*I_ERI_S_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fx2y_S_M2_vrr;
      Double I_ERI_Dxy_S_G2x2y_S_M1_vrr = PAY*I_ERI_Px_S_G2x2y_S_M1_vrr+WPY*I_ERI_Px_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Px_S_F2xy_S_M2_vrr;
      Double I_ERI_D2y_S_G2x2y_S_M1_vrr = PAY*I_ERI_Py_S_G2x2y_S_M1_vrr+WPY*I_ERI_Py_S_G2x2y_S_M2_vrr+oned2z*I_ERI_S_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M2_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M2_vrr;
      Double I_ERI_D2z_S_G2x2y_S_M1_vrr = PAZ*I_ERI_Pz_S_G2x2y_S_M1_vrr+WPZ*I_ERI_Pz_S_G2x2y_S_M2_vrr+oned2z*I_ERI_S_S_G2x2y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M2_vrr;
      Double I_ERI_D2x_S_G2xyz_S_M1_vrr = PAX*I_ERI_Px_S_G2xyz_S_M1_vrr+WPX*I_ERI_Px_S_G2xyz_S_M2_vrr+oned2z*I_ERI_S_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M2_vrr;
      Double I_ERI_Dxy_S_G2xyz_S_M1_vrr = PAY*I_ERI_Px_S_G2xyz_S_M1_vrr+WPY*I_ERI_Px_S_G2xyz_S_M2_vrr+oned2k*I_ERI_Px_S_F2xz_S_M2_vrr;
      Double I_ERI_D2y_S_G2xyz_S_M1_vrr = PAY*I_ERI_Py_S_G2xyz_S_M1_vrr+WPY*I_ERI_Py_S_G2xyz_S_M2_vrr+oned2z*I_ERI_S_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_Py_S_F2xz_S_M2_vrr;
      Double I_ERI_D2z_S_G2xyz_S_M1_vrr = PAZ*I_ERI_Pz_S_G2xyz_S_M1_vrr+WPZ*I_ERI_Pz_S_G2xyz_S_M2_vrr+oned2z*I_ERI_S_S_G2xyz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M2_vrr+oned2k*I_ERI_Pz_S_F2xy_S_M2_vrr;
      Double I_ERI_D2x_S_G2x2z_S_M1_vrr = PAX*I_ERI_Px_S_G2x2z_S_M1_vrr+WPX*I_ERI_Px_S_G2x2z_S_M2_vrr+oned2z*I_ERI_S_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fx2z_S_M2_vrr;
      Double I_ERI_Dxy_S_G2x2z_S_M1_vrr = PAY*I_ERI_Px_S_G2x2z_S_M1_vrr+WPY*I_ERI_Px_S_G2x2z_S_M2_vrr;
      Double I_ERI_D2y_S_G2x2z_S_M1_vrr = PAY*I_ERI_Py_S_G2x2z_S_M1_vrr+WPY*I_ERI_Py_S_G2x2z_S_M2_vrr+oned2z*I_ERI_S_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M2_vrr;
      Double I_ERI_D2z_S_G2x2z_S_M1_vrr = PAZ*I_ERI_Pz_S_G2x2z_S_M1_vrr+WPZ*I_ERI_Pz_S_G2x2z_S_M2_vrr+oned2z*I_ERI_S_S_G2x2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M2_vrr;
      Double I_ERI_D2x_S_Gx3y_S_M1_vrr = PAX*I_ERI_Px_S_Gx3y_S_M1_vrr+WPX*I_ERI_Px_S_Gx3y_S_M2_vrr+oned2z*I_ERI_S_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M2_vrr+oned2k*I_ERI_Px_S_F3y_S_M2_vrr;
      Double I_ERI_Dxy_S_Gx3y_S_M1_vrr = PAY*I_ERI_Px_S_Gx3y_S_M1_vrr+WPY*I_ERI_Px_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_Px_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2y_S_Gx3y_S_M1_vrr = PAY*I_ERI_Py_S_Gx3y_S_M1_vrr+WPY*I_ERI_Py_S_Gx3y_S_M2_vrr+oned2z*I_ERI_S_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M2_vrr+3*oned2k*I_ERI_Py_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2z_S_Gx3y_S_M1_vrr = PAZ*I_ERI_Pz_S_Gx3y_S_M1_vrr+WPZ*I_ERI_Pz_S_Gx3y_S_M2_vrr+oned2z*I_ERI_S_S_Gx3y_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M2_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_M1_vrr = PAX*I_ERI_Px_S_Gx2yz_S_M1_vrr+WPX*I_ERI_Px_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_Px_S_F2yz_S_M2_vrr;
      Double I_ERI_Dxy_S_Gx2yz_S_M1_vrr = PAY*I_ERI_Px_S_Gx2yz_S_M1_vrr+WPY*I_ERI_Px_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M2_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_M1_vrr = PAY*I_ERI_Py_S_Gx2yz_S_M1_vrr+WPY*I_ERI_Py_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M2_vrr+2*oned2k*I_ERI_Py_S_Fxyz_S_M2_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_M1_vrr = PAZ*I_ERI_Pz_S_Gx2yz_S_M1_vrr+WPZ*I_ERI_Pz_S_Gx2yz_S_M2_vrr+oned2z*I_ERI_S_S_Gx2yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M2_vrr+oned2k*I_ERI_Pz_S_Fx2y_S_M2_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_M1_vrr = PAX*I_ERI_Px_S_Gxy2z_S_M1_vrr+WPX*I_ERI_Px_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Px_S_Fy2z_S_M2_vrr;
      Double I_ERI_Dxy_S_Gxy2z_S_M1_vrr = PAY*I_ERI_Px_S_Gxy2z_S_M1_vrr+WPY*I_ERI_Px_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Px_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_M1_vrr = PAY*I_ERI_Py_S_Gxy2z_S_M1_vrr+WPY*I_ERI_Py_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M2_vrr+oned2k*I_ERI_Py_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_M1_vrr = PAZ*I_ERI_Pz_S_Gxy2z_S_M1_vrr+WPZ*I_ERI_Pz_S_Gxy2z_S_M2_vrr+oned2z*I_ERI_S_S_Gxy2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_Fxyz_S_M2_vrr;
      Double I_ERI_D2x_S_Gx3z_S_M1_vrr = PAX*I_ERI_Px_S_Gx3z_S_M1_vrr+WPX*I_ERI_Px_S_Gx3z_S_M2_vrr+oned2z*I_ERI_S_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M2_vrr+oned2k*I_ERI_Px_S_F3z_S_M2_vrr;
      Double I_ERI_Dxy_S_Gx3z_S_M1_vrr = PAY*I_ERI_Px_S_Gx3z_S_M1_vrr+WPY*I_ERI_Px_S_Gx3z_S_M2_vrr;
      Double I_ERI_D2y_S_Gx3z_S_M1_vrr = PAY*I_ERI_Py_S_Gx3z_S_M1_vrr+WPY*I_ERI_Py_S_Gx3z_S_M2_vrr+oned2z*I_ERI_S_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M2_vrr;
      Double I_ERI_D2z_S_Gx3z_S_M1_vrr = PAZ*I_ERI_Pz_S_Gx3z_S_M1_vrr+WPZ*I_ERI_Pz_S_Gx3z_S_M2_vrr+oned2z*I_ERI_S_S_Gx3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_Fx2z_S_M2_vrr;
      Double I_ERI_D2x_S_G4y_S_M1_vrr = PAX*I_ERI_Px_S_G4y_S_M1_vrr+WPX*I_ERI_Px_S_G4y_S_M2_vrr+oned2z*I_ERI_S_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_Dxy_S_G4y_S_M1_vrr = PAY*I_ERI_Px_S_G4y_S_M1_vrr+WPY*I_ERI_Px_S_G4y_S_M2_vrr+4*oned2k*I_ERI_Px_S_F3y_S_M2_vrr;
      Double I_ERI_D2y_S_G4y_S_M1_vrr = PAY*I_ERI_Py_S_G4y_S_M1_vrr+WPY*I_ERI_Py_S_G4y_S_M2_vrr+oned2z*I_ERI_S_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M2_vrr+4*oned2k*I_ERI_Py_S_F3y_S_M2_vrr;
      Double I_ERI_D2z_S_G4y_S_M1_vrr = PAZ*I_ERI_Pz_S_G4y_S_M1_vrr+WPZ*I_ERI_Pz_S_G4y_S_M2_vrr+oned2z*I_ERI_S_S_G4y_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M2_vrr;
      Double I_ERI_D2x_S_G3yz_S_M1_vrr = PAX*I_ERI_Px_S_G3yz_S_M1_vrr+WPX*I_ERI_Px_S_G3yz_S_M2_vrr+oned2z*I_ERI_S_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M2_vrr;
      Double I_ERI_Dxy_S_G3yz_S_M1_vrr = PAY*I_ERI_Px_S_G3yz_S_M1_vrr+WPY*I_ERI_Px_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_Px_S_F2yz_S_M2_vrr;
      Double I_ERI_D2y_S_G3yz_S_M1_vrr = PAY*I_ERI_Py_S_G3yz_S_M1_vrr+WPY*I_ERI_Py_S_G3yz_S_M2_vrr+oned2z*I_ERI_S_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M2_vrr+3*oned2k*I_ERI_Py_S_F2yz_S_M2_vrr;
      Double I_ERI_D2z_S_G3yz_S_M1_vrr = PAZ*I_ERI_Pz_S_G3yz_S_M1_vrr+WPZ*I_ERI_Pz_S_G3yz_S_M2_vrr+oned2z*I_ERI_S_S_G3yz_S_M1_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M2_vrr+oned2k*I_ERI_Pz_S_F3y_S_M2_vrr;
      Double I_ERI_D2x_S_G2y2z_S_M1_vrr = PAX*I_ERI_Px_S_G2y2z_S_M1_vrr+WPX*I_ERI_Px_S_G2y2z_S_M2_vrr+oned2z*I_ERI_S_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M2_vrr;
      Double I_ERI_Dxy_S_G2y2z_S_M1_vrr = PAY*I_ERI_Px_S_G2y2z_S_M1_vrr+WPY*I_ERI_Px_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_Px_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2y_S_G2y2z_S_M1_vrr = PAY*I_ERI_Py_S_G2y2z_S_M1_vrr+WPY*I_ERI_Py_S_G2y2z_S_M2_vrr+oned2z*I_ERI_S_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_Py_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2z_S_G2y2z_S_M1_vrr = PAZ*I_ERI_Pz_S_G2y2z_S_M1_vrr+WPZ*I_ERI_Pz_S_G2y2z_S_M2_vrr+oned2z*I_ERI_S_S_G2y2z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M2_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M2_vrr;
      Double I_ERI_D2x_S_Gy3z_S_M1_vrr = PAX*I_ERI_Px_S_Gy3z_S_M1_vrr+WPX*I_ERI_Px_S_Gy3z_S_M2_vrr+oned2z*I_ERI_S_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M2_vrr;
      Double I_ERI_Dxy_S_Gy3z_S_M1_vrr = PAY*I_ERI_Px_S_Gy3z_S_M1_vrr+WPY*I_ERI_Px_S_Gy3z_S_M2_vrr+oned2k*I_ERI_Px_S_F3z_S_M2_vrr;
      Double I_ERI_D2y_S_Gy3z_S_M1_vrr = PAY*I_ERI_Py_S_Gy3z_S_M1_vrr+WPY*I_ERI_Py_S_Gy3z_S_M2_vrr+oned2z*I_ERI_S_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M2_vrr+oned2k*I_ERI_Py_S_F3z_S_M2_vrr;
      Double I_ERI_D2z_S_Gy3z_S_M1_vrr = PAZ*I_ERI_Pz_S_Gy3z_S_M1_vrr+WPZ*I_ERI_Pz_S_Gy3z_S_M2_vrr+oned2z*I_ERI_S_S_Gy3z_S_M1_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M2_vrr+3*oned2k*I_ERI_Pz_S_Fy2z_S_M2_vrr;
      Double I_ERI_D2x_S_G4z_S_M1_vrr = PAX*I_ERI_Px_S_G4z_S_M1_vrr+WPX*I_ERI_Px_S_G4z_S_M2_vrr+oned2z*I_ERI_S_S_G4z_S_M1_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M2_vrr;
      Double I_ERI_Dxy_S_G4z_S_M1_vrr = PAY*I_ERI_Px_S_G4z_S_M1_vrr+WPY*I_ERI_Px_S_G4z_S_M2_vrr;
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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_S_D2x_S_vrr = QCX*I_ERI_S_S_Px_S_vrr+WQX*I_ERI_S_S_Px_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Dxy_S_vrr = QCY*I_ERI_S_S_Px_S_vrr+WQY*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_D2y_S_vrr = QCY*I_ERI_S_S_Py_S_vrr+WQY*I_ERI_S_S_Py_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_D2z_S_vrr = QCZ*I_ERI_S_S_Pz_S_vrr+WQZ*I_ERI_S_S_Pz_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_S_S_F_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_S_S_F3x_S_vrr = QCX*I_ERI_S_S_D2x_S_vrr+WQX*I_ERI_S_S_D2x_S_M1_vrr+2*oned2e*I_ERI_S_S_Px_S_vrr-2*rhod2esq*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_F2xy_S_vrr = QCY*I_ERI_S_S_D2x_S_vrr+WQY*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_F2xz_S_vrr = QCZ*I_ERI_S_S_D2x_S_vrr+WQZ*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_S_S_Fx2y_S_vrr = QCX*I_ERI_S_S_D2y_S_vrr+WQX*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_Fxyz_S_vrr = QCZ*I_ERI_S_S_Dxy_S_vrr+WQZ*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_S_S_Fx2z_S_vrr = QCX*I_ERI_S_S_D2z_S_vrr+WQX*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_S_F3y_S_vrr = QCY*I_ERI_S_S_D2y_S_vrr+WQY*I_ERI_S_S_D2y_S_M1_vrr+2*oned2e*I_ERI_S_S_Py_S_vrr-2*rhod2esq*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_S_F2yz_S_vrr = QCZ*I_ERI_S_S_D2y_S_vrr+WQZ*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_S_S_Fy2z_S_vrr = QCY*I_ERI_S_S_D2z_S_vrr+WQY*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_S_S_F3z_S_vrr = QCZ*I_ERI_S_S_D2z_S_vrr+WQZ*I_ERI_S_S_D2z_S_M1_vrr+2*oned2e*I_ERI_S_S_Pz_S_vrr-2*rhod2esq*I_ERI_S_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_F_S
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       ************************************************************/
      Double I_ERI_Px_S_F3x_S_vrr = PAX*I_ERI_S_S_F3x_S_vrr+WPX*I_ERI_S_S_F3x_S_M1_vrr+3*oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Py_S_F3x_S_vrr = PAY*I_ERI_S_S_F3x_S_vrr+WPY*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Pz_S_F3x_S_vrr = PAZ*I_ERI_S_S_F3x_S_vrr+WPZ*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Px_S_F2xy_S_vrr = PAX*I_ERI_S_S_F2xy_S_vrr+WPX*I_ERI_S_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_Py_S_F2xy_S_vrr = PAY*I_ERI_S_S_F2xy_S_vrr+WPY*I_ERI_S_S_F2xy_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Pz_S_F2xy_S_vrr = PAZ*I_ERI_S_S_F2xy_S_vrr+WPZ*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_Px_S_F2xz_S_vrr = PAX*I_ERI_S_S_F2xz_S_vrr+WPX*I_ERI_S_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Py_S_F2xz_S_vrr = PAY*I_ERI_S_S_F2xz_S_vrr+WPY*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Pz_S_F2xz_S_vrr = PAZ*I_ERI_S_S_F2xz_S_vrr+WPZ*I_ERI_S_S_F2xz_S_M1_vrr+oned2k*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Px_S_Fx2y_S_vrr = PAX*I_ERI_S_S_Fx2y_S_vrr+WPX*I_ERI_S_S_Fx2y_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Py_S_Fx2y_S_vrr = PAY*I_ERI_S_S_Fx2y_S_vrr+WPY*I_ERI_S_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_Pz_S_Fx2y_S_vrr = PAZ*I_ERI_S_S_Fx2y_S_vrr+WPZ*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_Px_S_Fxyz_S_vrr = PAX*I_ERI_S_S_Fxyz_S_vrr+WPX*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Py_S_Fxyz_S_vrr = PAY*I_ERI_S_S_Fxyz_S_vrr+WPY*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Pz_S_Fxyz_S_vrr = PAZ*I_ERI_S_S_Fxyz_S_vrr+WPZ*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_Px_S_Fx2z_S_vrr = PAX*I_ERI_S_S_Fx2z_S_vrr+WPX*I_ERI_S_S_Fx2z_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Py_S_Fx2z_S_vrr = PAY*I_ERI_S_S_Fx2z_S_vrr+WPY*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Pz_S_Fx2z_S_vrr = PAZ*I_ERI_S_S_Fx2z_S_vrr+WPZ*I_ERI_S_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Px_S_F3y_S_vrr = PAX*I_ERI_S_S_F3y_S_vrr+WPX*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Py_S_F3y_S_vrr = PAY*I_ERI_S_S_F3y_S_vrr+WPY*I_ERI_S_S_F3y_S_M1_vrr+3*oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Pz_S_F3y_S_vrr = PAZ*I_ERI_S_S_F3y_S_vrr+WPZ*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Px_S_F2yz_S_vrr = PAX*I_ERI_S_S_F2yz_S_vrr+WPX*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Py_S_F2yz_S_vrr = PAY*I_ERI_S_S_F2yz_S_vrr+WPY*I_ERI_S_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Pz_S_F2yz_S_vrr = PAZ*I_ERI_S_S_F2yz_S_vrr+WPZ*I_ERI_S_S_F2yz_S_M1_vrr+oned2k*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Px_S_Fy2z_S_vrr = PAX*I_ERI_S_S_Fy2z_S_vrr+WPX*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Py_S_Fy2z_S_vrr = PAY*I_ERI_S_S_Fy2z_S_vrr+WPY*I_ERI_S_S_Fy2z_S_M1_vrr+oned2k*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Pz_S_Fy2z_S_vrr = PAZ*I_ERI_S_S_Fy2z_S_vrr+WPZ*I_ERI_S_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Px_S_F3z_S_vrr = PAX*I_ERI_S_S_F3z_S_vrr+WPX*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Py_S_F3z_S_vrr = PAY*I_ERI_S_S_F3z_S_vrr+WPY*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Pz_S_F3z_S_vrr = PAZ*I_ERI_S_S_F3z_S_vrr+WPZ*I_ERI_S_S_F3z_S_M1_vrr+3*oned2k*I_ERI_S_S_D2z_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_D_S_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_F_S
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_F_S
       * RHS shell quartet name: SQ_ERI_S_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_F3x_S_vrr = PAX*I_ERI_Px_S_F3x_S_vrr+WPX*I_ERI_Px_S_F3x_S_M1_vrr+oned2z*I_ERI_S_S_F3x_S_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M1_vrr+3*oned2k*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Dxy_S_F3x_S_vrr = PAY*I_ERI_Px_S_F3x_S_vrr+WPY*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_Dxz_S_F3x_S_vrr = PAZ*I_ERI_Px_S_F3x_S_vrr+WPZ*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_D2y_S_F3x_S_vrr = PAY*I_ERI_Py_S_F3x_S_vrr+WPY*I_ERI_Py_S_F3x_S_M1_vrr+oned2z*I_ERI_S_S_F3x_S_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_Dyz_S_F3x_S_vrr = PAZ*I_ERI_Py_S_F3x_S_vrr+WPZ*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_D2z_S_F3x_S_vrr = PAZ*I_ERI_Pz_S_F3x_S_vrr+WPZ*I_ERI_Pz_S_F3x_S_M1_vrr+oned2z*I_ERI_S_S_F3x_S_vrr-rhod2zsq*I_ERI_S_S_F3x_S_M1_vrr;
      Double I_ERI_D2x_S_F2xy_S_vrr = PAX*I_ERI_Px_S_F2xy_S_vrr+WPX*I_ERI_Px_S_F2xy_S_M1_vrr+oned2z*I_ERI_S_S_F2xy_S_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_Dxy_S_F2xy_S_vrr = PAY*I_ERI_Px_S_F2xy_S_vrr+WPY*I_ERI_Px_S_F2xy_S_M1_vrr+oned2k*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Dxz_S_F2xy_S_vrr = PAZ*I_ERI_Px_S_F2xy_S_vrr+WPZ*I_ERI_Px_S_F2xy_S_M1_vrr;
      Double I_ERI_D2y_S_F2xy_S_vrr = PAY*I_ERI_Py_S_F2xy_S_vrr+WPY*I_ERI_Py_S_F2xy_S_M1_vrr+oned2z*I_ERI_S_S_F2xy_S_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M1_vrr+oned2k*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_Dyz_S_F2xy_S_vrr = PAZ*I_ERI_Py_S_F2xy_S_vrr+WPZ*I_ERI_Py_S_F2xy_S_M1_vrr;
      Double I_ERI_D2z_S_F2xy_S_vrr = PAZ*I_ERI_Pz_S_F2xy_S_vrr+WPZ*I_ERI_Pz_S_F2xy_S_M1_vrr+oned2z*I_ERI_S_S_F2xy_S_vrr-rhod2zsq*I_ERI_S_S_F2xy_S_M1_vrr;
      Double I_ERI_D2x_S_F2xz_S_vrr = PAX*I_ERI_Px_S_F2xz_S_vrr+WPX*I_ERI_Px_S_F2xz_S_M1_vrr+oned2z*I_ERI_S_S_F2xz_S_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M1_vrr;
      Double I_ERI_Dxy_S_F2xz_S_vrr = PAY*I_ERI_Px_S_F2xz_S_vrr+WPY*I_ERI_Px_S_F2xz_S_M1_vrr;
      Double I_ERI_Dxz_S_F2xz_S_vrr = PAZ*I_ERI_Px_S_F2xz_S_vrr+WPZ*I_ERI_Px_S_F2xz_S_M1_vrr+oned2k*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_D2y_S_F2xz_S_vrr = PAY*I_ERI_Py_S_F2xz_S_vrr+WPY*I_ERI_Py_S_F2xz_S_M1_vrr+oned2z*I_ERI_S_S_F2xz_S_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M1_vrr;
      Double I_ERI_Dyz_S_F2xz_S_vrr = PAZ*I_ERI_Py_S_F2xz_S_vrr+WPZ*I_ERI_Py_S_F2xz_S_M1_vrr+oned2k*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_D2z_S_F2xz_S_vrr = PAZ*I_ERI_Pz_S_F2xz_S_vrr+WPZ*I_ERI_Pz_S_F2xz_S_M1_vrr+oned2z*I_ERI_S_S_F2xz_S_vrr-rhod2zsq*I_ERI_S_S_F2xz_S_M1_vrr+oned2k*I_ERI_Pz_S_D2x_S_M1_vrr;
      Double I_ERI_D2x_S_Fx2y_S_vrr = PAX*I_ERI_Px_S_Fx2y_S_vrr+WPX*I_ERI_Px_S_Fx2y_S_M1_vrr+oned2z*I_ERI_S_S_Fx2y_S_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M1_vrr+oned2k*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_Dxy_S_Fx2y_S_vrr = PAY*I_ERI_Px_S_Fx2y_S_vrr+WPY*I_ERI_Px_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_Dxz_S_Fx2y_S_vrr = PAZ*I_ERI_Px_S_Fx2y_S_vrr+WPZ*I_ERI_Px_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2y_S_Fx2y_S_vrr = PAY*I_ERI_Py_S_Fx2y_S_vrr+WPY*I_ERI_Py_S_Fx2y_S_M1_vrr+oned2z*I_ERI_S_S_Fx2y_S_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_Py_S_Dxy_S_M1_vrr;
      Double I_ERI_Dyz_S_Fx2y_S_vrr = PAZ*I_ERI_Py_S_Fx2y_S_vrr+WPZ*I_ERI_Py_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2z_S_Fx2y_S_vrr = PAZ*I_ERI_Pz_S_Fx2y_S_vrr+WPZ*I_ERI_Pz_S_Fx2y_S_M1_vrr+oned2z*I_ERI_S_S_Fx2y_S_vrr-rhod2zsq*I_ERI_S_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2x_S_Fxyz_S_vrr = PAX*I_ERI_Px_S_Fxyz_S_vrr+WPX*I_ERI_Px_S_Fxyz_S_M1_vrr+oned2z*I_ERI_S_S_Fxyz_S_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Px_S_Dyz_S_M1_vrr;
      Double I_ERI_Dxy_S_Fxyz_S_vrr = PAY*I_ERI_Px_S_Fxyz_S_vrr+WPY*I_ERI_Px_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Px_S_Dxz_S_M1_vrr;
      Double I_ERI_Dxz_S_Fxyz_S_vrr = PAZ*I_ERI_Px_S_Fxyz_S_vrr+WPZ*I_ERI_Px_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_D2y_S_Fxyz_S_vrr = PAY*I_ERI_Py_S_Fxyz_S_vrr+WPY*I_ERI_Py_S_Fxyz_S_M1_vrr+oned2z*I_ERI_S_S_Fxyz_S_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Py_S_Dxz_S_M1_vrr;
      Double I_ERI_Dyz_S_Fxyz_S_vrr = PAZ*I_ERI_Py_S_Fxyz_S_vrr+WPZ*I_ERI_Py_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Py_S_Dxy_S_M1_vrr;
      Double I_ERI_D2z_S_Fxyz_S_vrr = PAZ*I_ERI_Pz_S_Fxyz_S_vrr+WPZ*I_ERI_Pz_S_Fxyz_S_M1_vrr+oned2z*I_ERI_S_S_Fxyz_S_vrr-rhod2zsq*I_ERI_S_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Pz_S_Dxy_S_M1_vrr;
      Double I_ERI_D2x_S_Fx2z_S_vrr = PAX*I_ERI_Px_S_Fx2z_S_vrr+WPX*I_ERI_Px_S_Fx2z_S_M1_vrr+oned2z*I_ERI_S_S_Fx2z_S_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M1_vrr+oned2k*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Dxy_S_Fx2z_S_vrr = PAY*I_ERI_Px_S_Fx2z_S_vrr+WPY*I_ERI_Px_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dxz_S_Fx2z_S_vrr = PAZ*I_ERI_Px_S_Fx2z_S_vrr+WPZ*I_ERI_Px_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dxz_S_M1_vrr;
      Double I_ERI_D2y_S_Fx2z_S_vrr = PAY*I_ERI_Py_S_Fx2z_S_vrr+WPY*I_ERI_Py_S_Fx2z_S_M1_vrr+oned2z*I_ERI_S_S_Fx2z_S_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dyz_S_Fx2z_S_vrr = PAZ*I_ERI_Py_S_Fx2z_S_vrr+WPZ*I_ERI_Py_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Dxz_S_M1_vrr;
      Double I_ERI_D2z_S_Fx2z_S_vrr = PAZ*I_ERI_Pz_S_Fx2z_S_vrr+WPZ*I_ERI_Pz_S_Fx2z_S_M1_vrr+oned2z*I_ERI_S_S_Fx2z_S_vrr-rhod2zsq*I_ERI_S_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Dxz_S_M1_vrr;
      Double I_ERI_D2x_S_F3y_S_vrr = PAX*I_ERI_Px_S_F3y_S_vrr+WPX*I_ERI_Px_S_F3y_S_M1_vrr+oned2z*I_ERI_S_S_F3y_S_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_Dxy_S_F3y_S_vrr = PAY*I_ERI_Px_S_F3y_S_vrr+WPY*I_ERI_Px_S_F3y_S_M1_vrr+3*oned2k*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_Dxz_S_F3y_S_vrr = PAZ*I_ERI_Px_S_F3y_S_vrr+WPZ*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_D2y_S_F3y_S_vrr = PAY*I_ERI_Py_S_F3y_S_vrr+WPY*I_ERI_Py_S_F3y_S_M1_vrr+oned2z*I_ERI_S_S_F3y_S_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M1_vrr+3*oned2k*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_Dyz_S_F3y_S_vrr = PAZ*I_ERI_Py_S_F3y_S_vrr+WPZ*I_ERI_Py_S_F3y_S_M1_vrr;
      Double I_ERI_D2z_S_F3y_S_vrr = PAZ*I_ERI_Pz_S_F3y_S_vrr+WPZ*I_ERI_Pz_S_F3y_S_M1_vrr+oned2z*I_ERI_S_S_F3y_S_vrr-rhod2zsq*I_ERI_S_S_F3y_S_M1_vrr;
      Double I_ERI_D2x_S_F2yz_S_vrr = PAX*I_ERI_Px_S_F2yz_S_vrr+WPX*I_ERI_Px_S_F2yz_S_M1_vrr+oned2z*I_ERI_S_S_F2yz_S_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M1_vrr;
      Double I_ERI_Dxy_S_F2yz_S_vrr = PAY*I_ERI_Px_S_F2yz_S_vrr+WPY*I_ERI_Px_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dyz_S_M1_vrr;
      Double I_ERI_Dxz_S_F2yz_S_vrr = PAZ*I_ERI_Px_S_F2yz_S_vrr+WPZ*I_ERI_Px_S_F2yz_S_M1_vrr+oned2k*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_D2y_S_F2yz_S_vrr = PAY*I_ERI_Py_S_F2yz_S_vrr+WPY*I_ERI_Py_S_F2yz_S_M1_vrr+oned2z*I_ERI_S_S_F2yz_S_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M1_vrr;
      Double I_ERI_Dyz_S_F2yz_S_vrr = PAZ*I_ERI_Py_S_F2yz_S_vrr+WPZ*I_ERI_Py_S_F2yz_S_M1_vrr+oned2k*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_D2z_S_F2yz_S_vrr = PAZ*I_ERI_Pz_S_F2yz_S_vrr+WPZ*I_ERI_Pz_S_F2yz_S_M1_vrr+oned2z*I_ERI_S_S_F2yz_S_vrr-rhod2zsq*I_ERI_S_S_F2yz_S_M1_vrr+oned2k*I_ERI_Pz_S_D2y_S_M1_vrr;
      Double I_ERI_D2x_S_Fy2z_S_vrr = PAX*I_ERI_Px_S_Fy2z_S_vrr+WPX*I_ERI_Px_S_Fy2z_S_M1_vrr+oned2z*I_ERI_S_S_Fy2z_S_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M1_vrr;
      Double I_ERI_Dxy_S_Fy2z_S_vrr = PAY*I_ERI_Px_S_Fy2z_S_vrr+WPY*I_ERI_Px_S_Fy2z_S_M1_vrr+oned2k*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Dxz_S_Fy2z_S_vrr = PAZ*I_ERI_Px_S_Fy2z_S_vrr+WPZ*I_ERI_Px_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Dyz_S_M1_vrr;
      Double I_ERI_D2y_S_Fy2z_S_vrr = PAY*I_ERI_Py_S_Fy2z_S_vrr+WPY*I_ERI_Py_S_Fy2z_S_M1_vrr+oned2z*I_ERI_S_S_Fy2z_S_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M1_vrr+oned2k*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_Dyz_S_Fy2z_S_vrr = PAZ*I_ERI_Py_S_Fy2z_S_vrr+WPZ*I_ERI_Py_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Dyz_S_M1_vrr;
      Double I_ERI_D2z_S_Fy2z_S_vrr = PAZ*I_ERI_Pz_S_Fy2z_S_vrr+WPZ*I_ERI_Pz_S_Fy2z_S_M1_vrr+oned2z*I_ERI_S_S_Fy2z_S_vrr-rhod2zsq*I_ERI_S_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Dyz_S_M1_vrr;
      Double I_ERI_D2x_S_F3z_S_vrr = PAX*I_ERI_Px_S_F3z_S_vrr+WPX*I_ERI_Px_S_F3z_S_M1_vrr+oned2z*I_ERI_S_S_F3z_S_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Dxy_S_F3z_S_vrr = PAY*I_ERI_Px_S_F3z_S_vrr+WPY*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_Dxz_S_F3z_S_vrr = PAZ*I_ERI_Px_S_F3z_S_vrr+WPZ*I_ERI_Px_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_D2y_S_F3z_S_vrr = PAY*I_ERI_Py_S_F3z_S_vrr+WPY*I_ERI_Py_S_F3z_S_M1_vrr+oned2z*I_ERI_S_S_F3z_S_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M1_vrr;
      Double I_ERI_Dyz_S_F3z_S_vrr = PAZ*I_ERI_Py_S_F3z_S_vrr+WPZ*I_ERI_Py_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_D2z_S_F3z_S_vrr = PAZ*I_ERI_Pz_S_F3z_S_vrr+WPZ*I_ERI_Pz_S_F3z_S_M1_vrr+oned2z*I_ERI_S_S_F3z_S_vrr-rhod2zsq*I_ERI_S_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Pz_S_D2z_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_F_S_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_F_S
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_F3x_S_vrr = PAX*I_ERI_D2x_S_F3x_S_vrr+WPX*I_ERI_D2x_S_F3x_S_M1_vrr+2*oned2z*I_ERI_Px_S_F3x_S_vrr-2*rhod2zsq*I_ERI_Px_S_F3x_S_M1_vrr+3*oned2k*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xy_S_F3x_S_vrr = PAY*I_ERI_D2x_S_F3x_S_vrr+WPY*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xz_S_F3x_S_vrr = PAZ*I_ERI_D2x_S_F3x_S_vrr+WPZ*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2y_S_F3x_S_vrr = PAX*I_ERI_D2y_S_F3x_S_vrr+WPX*I_ERI_D2y_S_F3x_S_M1_vrr+3*oned2k*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fxyz_S_F3x_S_vrr = PAZ*I_ERI_Dxy_S_F3x_S_vrr+WPZ*I_ERI_Dxy_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2z_S_F3x_S_vrr = PAX*I_ERI_D2z_S_F3x_S_vrr+WPX*I_ERI_D2z_S_F3x_S_M1_vrr+3*oned2k*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3y_S_F3x_S_vrr = PAY*I_ERI_D2y_S_F3x_S_vrr+WPY*I_ERI_D2y_S_F3x_S_M1_vrr+2*oned2z*I_ERI_Py_S_F3x_S_vrr-2*rhod2zsq*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_F2yz_S_F3x_S_vrr = PAZ*I_ERI_D2y_S_F3x_S_vrr+WPZ*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Fy2z_S_F3x_S_vrr = PAY*I_ERI_D2z_S_F3x_S_vrr+WPY*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3z_S_F3x_S_vrr = PAZ*I_ERI_D2z_S_F3x_S_vrr+WPZ*I_ERI_D2z_S_F3x_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F3x_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F3x_S_M1_vrr;
      Double I_ERI_F3x_S_F2xy_S_vrr = PAX*I_ERI_D2x_S_F2xy_S_vrr+WPX*I_ERI_D2x_S_F2xy_S_M1_vrr+2*oned2z*I_ERI_Px_S_F2xy_S_vrr-2*rhod2zsq*I_ERI_Px_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_F2xy_S_F2xy_S_vrr = PAY*I_ERI_D2x_S_F2xy_S_vrr+WPY*I_ERI_D2x_S_F2xy_S_M1_vrr+oned2k*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xz_S_F2xy_S_vrr = PAZ*I_ERI_D2x_S_F2xy_S_vrr+WPZ*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2y_S_F2xy_S_vrr = PAX*I_ERI_D2y_S_F2xy_S_vrr+WPX*I_ERI_D2y_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Fxyz_S_F2xy_S_vrr = PAZ*I_ERI_Dxy_S_F2xy_S_vrr+WPZ*I_ERI_Dxy_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2z_S_F2xy_S_vrr = PAX*I_ERI_D2z_S_F2xy_S_vrr+WPX*I_ERI_D2z_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dxy_S_M1_vrr;
      Double I_ERI_F3y_S_F2xy_S_vrr = PAY*I_ERI_D2y_S_F2xy_S_vrr+WPY*I_ERI_D2y_S_F2xy_S_M1_vrr+2*oned2z*I_ERI_Py_S_F2xy_S_vrr-2*rhod2zsq*I_ERI_Py_S_F2xy_S_M1_vrr+oned2k*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_F2yz_S_F2xy_S_vrr = PAZ*I_ERI_D2y_S_F2xy_S_vrr+WPZ*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Fy2z_S_F2xy_S_vrr = PAY*I_ERI_D2z_S_F2xy_S_vrr+WPY*I_ERI_D2z_S_F2xy_S_M1_vrr+oned2k*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3z_S_F2xy_S_vrr = PAZ*I_ERI_D2z_S_F2xy_S_vrr+WPZ*I_ERI_D2z_S_F2xy_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F2xy_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xy_S_M1_vrr;
      Double I_ERI_F3x_S_F2xz_S_vrr = PAX*I_ERI_D2x_S_F2xz_S_vrr+WPX*I_ERI_D2x_S_F2xz_S_M1_vrr+2*oned2z*I_ERI_Px_S_F2xz_S_vrr-2*rhod2zsq*I_ERI_Px_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M1_vrr;
      Double I_ERI_F2xy_S_F2xz_S_vrr = PAY*I_ERI_D2x_S_F2xz_S_vrr+WPY*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_F2xz_S_F2xz_S_vrr = PAZ*I_ERI_D2x_S_F2xz_S_vrr+WPZ*I_ERI_D2x_S_F2xz_S_M1_vrr+oned2k*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2y_S_F2xz_S_vrr = PAX*I_ERI_D2y_S_F2xz_S_vrr+WPX*I_ERI_D2y_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M1_vrr;
      Double I_ERI_Fxyz_S_F2xz_S_vrr = PAZ*I_ERI_Dxy_S_F2xz_S_vrr+WPZ*I_ERI_Dxy_S_F2xz_S_M1_vrr+oned2k*I_ERI_Dxy_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2z_S_F2xz_S_vrr = PAX*I_ERI_D2z_S_F2xz_S_vrr+WPX*I_ERI_D2z_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M1_vrr;
      Double I_ERI_F3y_S_F2xz_S_vrr = PAY*I_ERI_D2y_S_F2xz_S_vrr+WPY*I_ERI_D2y_S_F2xz_S_M1_vrr+2*oned2z*I_ERI_Py_S_F2xz_S_vrr-2*rhod2zsq*I_ERI_Py_S_F2xz_S_M1_vrr;
      Double I_ERI_F2yz_S_F2xz_S_vrr = PAZ*I_ERI_D2y_S_F2xz_S_vrr+WPZ*I_ERI_D2y_S_F2xz_S_M1_vrr+oned2k*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fy2z_S_F2xz_S_vrr = PAY*I_ERI_D2z_S_F2xz_S_vrr+WPY*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3z_S_F2xz_S_vrr = PAZ*I_ERI_D2z_S_F2xz_S_vrr+WPZ*I_ERI_D2z_S_F2xz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F2xz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F2xz_S_M1_vrr+oned2k*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3x_S_Fx2y_S_vrr = PAX*I_ERI_D2x_S_Fx2y_S_vrr+WPX*I_ERI_D2x_S_Fx2y_S_M1_vrr+2*oned2z*I_ERI_Px_S_Fx2y_S_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2y_S_M1_vrr+oned2k*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_F2xy_S_Fx2y_S_vrr = PAY*I_ERI_D2x_S_Fx2y_S_vrr+WPY*I_ERI_D2x_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_F2xz_S_Fx2y_S_vrr = PAZ*I_ERI_D2x_S_Fx2y_S_vrr+WPZ*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fx2y_S_vrr = PAX*I_ERI_D2y_S_Fx2y_S_vrr+WPX*I_ERI_D2y_S_Fx2y_S_M1_vrr+oned2k*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fx2y_S_vrr = PAZ*I_ERI_Dxy_S_Fx2y_S_vrr+WPZ*I_ERI_Dxy_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fx2y_S_vrr = PAX*I_ERI_D2z_S_Fx2y_S_vrr+WPX*I_ERI_D2z_S_Fx2y_S_M1_vrr+oned2k*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3y_S_Fx2y_S_vrr = PAY*I_ERI_D2y_S_Fx2y_S_vrr+WPY*I_ERI_D2y_S_Fx2y_S_M1_vrr+2*oned2z*I_ERI_Py_S_Fx2y_S_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_F2yz_S_Fx2y_S_vrr = PAZ*I_ERI_D2y_S_Fx2y_S_vrr+WPZ*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fx2y_S_vrr = PAY*I_ERI_D2z_S_Fx2y_S_vrr+WPY*I_ERI_D2z_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dxy_S_M1_vrr;
      Double I_ERI_F3z_S_Fx2y_S_vrr = PAZ*I_ERI_D2z_S_Fx2y_S_vrr+WPZ*I_ERI_D2z_S_Fx2y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Fx2y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3x_S_Fxyz_S_vrr = PAX*I_ERI_D2x_S_Fxyz_S_vrr+WPX*I_ERI_D2x_S_Fxyz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Fxyz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Dyz_S_M1_vrr;
      Double I_ERI_F2xy_S_Fxyz_S_vrr = PAY*I_ERI_D2x_S_Fxyz_S_vrr+WPY*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Dxz_S_M1_vrr;
      Double I_ERI_F2xz_S_Fxyz_S_vrr = PAZ*I_ERI_D2x_S_Fxyz_S_vrr+WPZ*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fxyz_S_vrr = PAX*I_ERI_D2y_S_Fxyz_S_vrr+WPX*I_ERI_D2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Dyz_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fxyz_S_vrr = PAZ*I_ERI_Dxy_S_Fxyz_S_vrr+WPZ*I_ERI_Dxy_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Dxy_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fxyz_S_vrr = PAX*I_ERI_D2z_S_Fxyz_S_vrr+WPX*I_ERI_D2z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Dyz_S_M1_vrr;
      Double I_ERI_F3y_S_Fxyz_S_vrr = PAY*I_ERI_D2y_S_Fxyz_S_vrr+WPY*I_ERI_D2y_S_Fxyz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Fxyz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Dxz_S_M1_vrr;
      Double I_ERI_F2yz_S_Fxyz_S_vrr = PAZ*I_ERI_D2y_S_Fxyz_S_vrr+WPZ*I_ERI_D2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fxyz_S_vrr = PAY*I_ERI_D2z_S_Fxyz_S_vrr+WPY*I_ERI_D2z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Dxz_S_M1_vrr;
      Double I_ERI_F3z_S_Fxyz_S_vrr = PAZ*I_ERI_D2z_S_Fxyz_S_vrr+WPZ*I_ERI_D2z_S_Fxyz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Fxyz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Fxyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Dxy_S_M1_vrr;
      Double I_ERI_F3x_S_Fx2z_S_vrr = PAX*I_ERI_D2x_S_Fx2z_S_vrr+WPX*I_ERI_D2x_S_Fx2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Fx2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Fx2z_S_M1_vrr+oned2k*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Fx2z_S_vrr = PAY*I_ERI_D2x_S_Fx2z_S_vrr+WPY*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Fx2z_S_vrr = PAZ*I_ERI_D2x_S_Fx2z_S_vrr+WPZ*I_ERI_D2x_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dxz_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fx2z_S_vrr = PAX*I_ERI_D2y_S_Fx2z_S_vrr+WPX*I_ERI_D2y_S_Fx2z_S_M1_vrr+oned2k*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fx2z_S_vrr = PAZ*I_ERI_Dxy_S_Fx2z_S_vrr+WPZ*I_ERI_Dxy_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Dxz_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fx2z_S_vrr = PAX*I_ERI_D2z_S_Fx2z_S_vrr+WPX*I_ERI_D2z_S_Fx2z_S_M1_vrr+oned2k*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3y_S_Fx2z_S_vrr = PAY*I_ERI_D2y_S_Fx2z_S_vrr+WPY*I_ERI_D2y_S_Fx2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Fx2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Fx2z_S_vrr = PAZ*I_ERI_D2y_S_Fx2z_S_vrr+WPZ*I_ERI_D2y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dxz_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fx2z_S_vrr = PAY*I_ERI_D2z_S_Fx2z_S_vrr+WPY*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3z_S_Fx2z_S_vrr = PAZ*I_ERI_D2z_S_Fx2z_S_vrr+WPZ*I_ERI_D2z_S_Fx2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Fx2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dxz_S_M1_vrr;
      Double I_ERI_F3x_S_F3y_S_vrr = PAX*I_ERI_D2x_S_F3y_S_vrr+WPX*I_ERI_D2x_S_F3y_S_M1_vrr+2*oned2z*I_ERI_Px_S_F3y_S_vrr-2*rhod2zsq*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_F2xy_S_F3y_S_vrr = PAY*I_ERI_D2x_S_F3y_S_vrr+WPY*I_ERI_D2x_S_F3y_S_M1_vrr+3*oned2k*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_F2xz_S_F3y_S_vrr = PAZ*I_ERI_D2x_S_F3y_S_vrr+WPZ*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_F3y_S_vrr = PAX*I_ERI_D2y_S_F3y_S_vrr+WPX*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fxyz_S_F3y_S_vrr = PAZ*I_ERI_Dxy_S_F3y_S_vrr+WPZ*I_ERI_Dxy_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2z_S_F3y_S_vrr = PAX*I_ERI_D2z_S_F3y_S_vrr+WPX*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3y_S_F3y_S_vrr = PAY*I_ERI_D2y_S_F3y_S_vrr+WPY*I_ERI_D2y_S_F3y_S_M1_vrr+2*oned2z*I_ERI_Py_S_F3y_S_vrr-2*rhod2zsq*I_ERI_Py_S_F3y_S_M1_vrr+3*oned2k*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_F2yz_S_F3y_S_vrr = PAZ*I_ERI_D2y_S_F3y_S_vrr+WPZ*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fy2z_S_F3y_S_vrr = PAY*I_ERI_D2z_S_F3y_S_vrr+WPY*I_ERI_D2z_S_F3y_S_M1_vrr+3*oned2k*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3z_S_F3y_S_vrr = PAZ*I_ERI_D2z_S_F3y_S_vrr+WPZ*I_ERI_D2z_S_F3y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F3y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F3y_S_M1_vrr;
      Double I_ERI_F3x_S_F2yz_S_vrr = PAX*I_ERI_D2x_S_F2yz_S_vrr+WPX*I_ERI_D2x_S_F2yz_S_M1_vrr+2*oned2z*I_ERI_Px_S_F2yz_S_vrr-2*rhod2zsq*I_ERI_Px_S_F2yz_S_M1_vrr;
      Double I_ERI_F2xy_S_F2yz_S_vrr = PAY*I_ERI_D2x_S_F2yz_S_vrr+WPY*I_ERI_D2x_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M1_vrr;
      Double I_ERI_F2xz_S_F2yz_S_vrr = PAZ*I_ERI_D2x_S_F2yz_S_vrr+WPZ*I_ERI_D2x_S_F2yz_S_M1_vrr+oned2k*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_F2yz_S_vrr = PAX*I_ERI_D2y_S_F2yz_S_vrr+WPX*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_Fxyz_S_F2yz_S_vrr = PAZ*I_ERI_Dxy_S_F2yz_S_vrr+WPZ*I_ERI_Dxy_S_F2yz_S_M1_vrr+oned2k*I_ERI_Dxy_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_F2yz_S_vrr = PAX*I_ERI_D2z_S_F2yz_S_vrr+WPX*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3y_S_F2yz_S_vrr = PAY*I_ERI_D2y_S_F2yz_S_vrr+WPY*I_ERI_D2y_S_F2yz_S_M1_vrr+2*oned2z*I_ERI_Py_S_F2yz_S_vrr-2*rhod2zsq*I_ERI_Py_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M1_vrr;
      Double I_ERI_F2yz_S_F2yz_S_vrr = PAZ*I_ERI_D2y_S_F2yz_S_vrr+WPZ*I_ERI_D2y_S_F2yz_S_M1_vrr+oned2k*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_F2yz_S_vrr = PAY*I_ERI_D2z_S_F2yz_S_vrr+WPY*I_ERI_D2z_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M1_vrr;
      Double I_ERI_F3z_S_F2yz_S_vrr = PAZ*I_ERI_D2z_S_F2yz_S_vrr+WPZ*I_ERI_D2z_S_F2yz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F2yz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F2yz_S_M1_vrr+oned2k*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3x_S_Fy2z_S_vrr = PAX*I_ERI_D2x_S_Fy2z_S_vrr+WPX*I_ERI_D2x_S_Fy2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Fy2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Fy2z_S_vrr = PAY*I_ERI_D2x_S_Fy2z_S_vrr+WPY*I_ERI_D2x_S_Fy2z_S_M1_vrr+oned2k*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Fy2z_S_vrr = PAZ*I_ERI_D2x_S_Fy2z_S_vrr+WPZ*I_ERI_D2x_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Dyz_S_M1_vrr;
      Double I_ERI_Fx2y_S_Fy2z_S_vrr = PAX*I_ERI_D2y_S_Fy2z_S_vrr+WPX*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Fy2z_S_vrr = PAZ*I_ERI_Dxy_S_Fy2z_S_vrr+WPZ*I_ERI_Dxy_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Dyz_S_M1_vrr;
      Double I_ERI_Fx2z_S_Fy2z_S_vrr = PAX*I_ERI_D2z_S_Fy2z_S_vrr+WPX*I_ERI_D2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3y_S_Fy2z_S_vrr = PAY*I_ERI_D2y_S_Fy2z_S_vrr+WPY*I_ERI_D2y_S_Fy2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Fy2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Fy2z_S_M1_vrr+oned2k*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Fy2z_S_vrr = PAZ*I_ERI_D2y_S_Fy2z_S_vrr+WPZ*I_ERI_D2y_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Dyz_S_M1_vrr;
      Double I_ERI_Fy2z_S_Fy2z_S_vrr = PAY*I_ERI_D2z_S_Fy2z_S_vrr+WPY*I_ERI_D2z_S_Fy2z_S_M1_vrr+oned2k*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3z_S_Fy2z_S_vrr = PAZ*I_ERI_D2z_S_Fy2z_S_vrr+WPZ*I_ERI_D2z_S_Fy2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Fy2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Dyz_S_M1_vrr;
      Double I_ERI_F3x_S_F3z_S_vrr = PAX*I_ERI_D2x_S_F3z_S_vrr+WPX*I_ERI_D2x_S_F3z_S_M1_vrr+2*oned2z*I_ERI_Px_S_F3z_S_vrr-2*rhod2zsq*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_F2xy_S_F3z_S_vrr = PAY*I_ERI_D2x_S_F3z_S_vrr+WPY*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xz_S_F3z_S_vrr = PAZ*I_ERI_D2x_S_F3z_S_vrr+WPZ*I_ERI_D2x_S_F3z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_F3z_S_vrr = PAX*I_ERI_D2y_S_F3z_S_vrr+WPX*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_Fxyz_S_F3z_S_vrr = PAZ*I_ERI_Dxy_S_F3z_S_vrr+WPZ*I_ERI_Dxy_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Dxy_S_D2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_F3z_S_vrr = PAX*I_ERI_D2z_S_F3z_S_vrr+WPX*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3y_S_F3z_S_vrr = PAY*I_ERI_D2y_S_F3z_S_vrr+WPY*I_ERI_D2y_S_F3z_S_M1_vrr+2*oned2z*I_ERI_Py_S_F3z_S_vrr-2*rhod2zsq*I_ERI_Py_S_F3z_S_M1_vrr;
      Double I_ERI_F2yz_S_F3z_S_vrr = PAZ*I_ERI_D2y_S_F3z_S_vrr+WPZ*I_ERI_D2y_S_F3z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Fy2z_S_F3z_S_vrr = PAY*I_ERI_D2z_S_F3z_S_vrr+WPY*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3z_S_F3z_S_vrr = PAZ*I_ERI_D2z_S_F3z_S_vrr+WPZ*I_ERI_D2z_S_F3z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_F3z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_F3z_S_M1_vrr+3*oned2k*I_ERI_D2z_S_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_G_S
       * RHS shell quartet name: SQ_ERI_P_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_G_S
       * RHS shell quartet name: SQ_ERI_S_S_G_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_F_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_G4x_S_vrr = PAX*I_ERI_Px_S_G4x_S_vrr+WPX*I_ERI_Px_S_G4x_S_M1_vrr+oned2z*I_ERI_S_S_G4x_S_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M1_vrr+4*oned2k*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_Dxy_S_G4x_S_vrr = PAY*I_ERI_Px_S_G4x_S_vrr+WPY*I_ERI_Px_S_G4x_S_M1_vrr;
      Double I_ERI_Dxz_S_G4x_S_vrr = PAZ*I_ERI_Px_S_G4x_S_vrr+WPZ*I_ERI_Px_S_G4x_S_M1_vrr;
      Double I_ERI_D2y_S_G4x_S_vrr = PAY*I_ERI_Py_S_G4x_S_vrr+WPY*I_ERI_Py_S_G4x_S_M1_vrr+oned2z*I_ERI_S_S_G4x_S_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_Dyz_S_G4x_S_vrr = PAZ*I_ERI_Py_S_G4x_S_vrr+WPZ*I_ERI_Py_S_G4x_S_M1_vrr;
      Double I_ERI_D2z_S_G4x_S_vrr = PAZ*I_ERI_Pz_S_G4x_S_vrr+WPZ*I_ERI_Pz_S_G4x_S_M1_vrr+oned2z*I_ERI_S_S_G4x_S_vrr-rhod2zsq*I_ERI_S_S_G4x_S_M1_vrr;
      Double I_ERI_D2x_S_G3xy_S_vrr = PAX*I_ERI_Px_S_G3xy_S_vrr+WPX*I_ERI_Px_S_G3xy_S_M1_vrr+oned2z*I_ERI_S_S_G3xy_S_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_Px_S_F2xy_S_M1_vrr;
      Double I_ERI_Dxy_S_G3xy_S_vrr = PAY*I_ERI_Px_S_G3xy_S_vrr+WPY*I_ERI_Px_S_G3xy_S_M1_vrr+oned2k*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_Dxz_S_G3xy_S_vrr = PAZ*I_ERI_Px_S_G3xy_S_vrr+WPZ*I_ERI_Px_S_G3xy_S_M1_vrr;
      Double I_ERI_D2y_S_G3xy_S_vrr = PAY*I_ERI_Py_S_G3xy_S_vrr+WPY*I_ERI_Py_S_G3xy_S_M1_vrr+oned2z*I_ERI_S_S_G3xy_S_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M1_vrr+oned2k*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_Dyz_S_G3xy_S_vrr = PAZ*I_ERI_Py_S_G3xy_S_vrr+WPZ*I_ERI_Py_S_G3xy_S_M1_vrr;
      Double I_ERI_D2z_S_G3xy_S_vrr = PAZ*I_ERI_Pz_S_G3xy_S_vrr+WPZ*I_ERI_Pz_S_G3xy_S_M1_vrr+oned2z*I_ERI_S_S_G3xy_S_vrr-rhod2zsq*I_ERI_S_S_G3xy_S_M1_vrr;
      Double I_ERI_D2x_S_G3xz_S_vrr = PAX*I_ERI_Px_S_G3xz_S_vrr+WPX*I_ERI_Px_S_G3xz_S_M1_vrr+oned2z*I_ERI_S_S_G3xz_S_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_Px_S_F2xz_S_M1_vrr;
      Double I_ERI_Dxy_S_G3xz_S_vrr = PAY*I_ERI_Px_S_G3xz_S_vrr+WPY*I_ERI_Px_S_G3xz_S_M1_vrr;
      Double I_ERI_Dxz_S_G3xz_S_vrr = PAZ*I_ERI_Px_S_G3xz_S_vrr+WPZ*I_ERI_Px_S_G3xz_S_M1_vrr+oned2k*I_ERI_Px_S_F3x_S_M1_vrr;
      Double I_ERI_D2y_S_G3xz_S_vrr = PAY*I_ERI_Py_S_G3xz_S_vrr+WPY*I_ERI_Py_S_G3xz_S_M1_vrr+oned2z*I_ERI_S_S_G3xz_S_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M1_vrr;
      Double I_ERI_Dyz_S_G3xz_S_vrr = PAZ*I_ERI_Py_S_G3xz_S_vrr+WPZ*I_ERI_Py_S_G3xz_S_M1_vrr+oned2k*I_ERI_Py_S_F3x_S_M1_vrr;
      Double I_ERI_D2z_S_G3xz_S_vrr = PAZ*I_ERI_Pz_S_G3xz_S_vrr+WPZ*I_ERI_Pz_S_G3xz_S_M1_vrr+oned2z*I_ERI_S_S_G3xz_S_vrr-rhod2zsq*I_ERI_S_S_G3xz_S_M1_vrr+oned2k*I_ERI_Pz_S_F3x_S_M1_vrr;
      Double I_ERI_D2x_S_G2x2y_S_vrr = PAX*I_ERI_Px_S_G2x2y_S_vrr+WPX*I_ERI_Px_S_G2x2y_S_M1_vrr+oned2z*I_ERI_S_S_G2x2y_S_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fx2y_S_M1_vrr;
      Double I_ERI_Dxy_S_G2x2y_S_vrr = PAY*I_ERI_Px_S_G2x2y_S_vrr+WPY*I_ERI_Px_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_F2xy_S_M1_vrr;
      Double I_ERI_Dxz_S_G2x2y_S_vrr = PAZ*I_ERI_Px_S_G2x2y_S_vrr+WPZ*I_ERI_Px_S_G2x2y_S_M1_vrr;
      Double I_ERI_D2y_S_G2x2y_S_vrr = PAY*I_ERI_Py_S_G2x2y_S_vrr+WPY*I_ERI_Py_S_G2x2y_S_M1_vrr+oned2z*I_ERI_S_S_G2x2y_S_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_Py_S_F2xy_S_M1_vrr;
      Double I_ERI_Dyz_S_G2x2y_S_vrr = PAZ*I_ERI_Py_S_G2x2y_S_vrr+WPZ*I_ERI_Py_S_G2x2y_S_M1_vrr;
      Double I_ERI_D2z_S_G2x2y_S_vrr = PAZ*I_ERI_Pz_S_G2x2y_S_vrr+WPZ*I_ERI_Pz_S_G2x2y_S_M1_vrr+oned2z*I_ERI_S_S_G2x2y_S_vrr-rhod2zsq*I_ERI_S_S_G2x2y_S_M1_vrr;
      Double I_ERI_D2x_S_G2xyz_S_vrr = PAX*I_ERI_Px_S_G2xyz_S_vrr+WPX*I_ERI_Px_S_G2xyz_S_M1_vrr+oned2z*I_ERI_S_S_G2xyz_S_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M1_vrr;
      Double I_ERI_Dxy_S_G2xyz_S_vrr = PAY*I_ERI_Px_S_G2xyz_S_vrr+WPY*I_ERI_Px_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Px_S_F2xz_S_M1_vrr;
      Double I_ERI_Dxz_S_G2xyz_S_vrr = PAZ*I_ERI_Px_S_G2xyz_S_vrr+WPZ*I_ERI_Px_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Px_S_F2xy_S_M1_vrr;
      Double I_ERI_D2y_S_G2xyz_S_vrr = PAY*I_ERI_Py_S_G2xyz_S_vrr+WPY*I_ERI_Py_S_G2xyz_S_M1_vrr+oned2z*I_ERI_S_S_G2xyz_S_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Py_S_F2xz_S_M1_vrr;
      Double I_ERI_Dyz_S_G2xyz_S_vrr = PAZ*I_ERI_Py_S_G2xyz_S_vrr+WPZ*I_ERI_Py_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Py_S_F2xy_S_M1_vrr;
      Double I_ERI_D2z_S_G2xyz_S_vrr = PAZ*I_ERI_Pz_S_G2xyz_S_vrr+WPZ*I_ERI_Pz_S_G2xyz_S_M1_vrr+oned2z*I_ERI_S_S_G2xyz_S_vrr-rhod2zsq*I_ERI_S_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Pz_S_F2xy_S_M1_vrr;
      Double I_ERI_D2x_S_G2x2z_S_vrr = PAX*I_ERI_Px_S_G2x2z_S_vrr+WPX*I_ERI_Px_S_G2x2z_S_M1_vrr+oned2z*I_ERI_S_S_G2x2z_S_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dxy_S_G2x2z_S_vrr = PAY*I_ERI_Px_S_G2x2z_S_vrr+WPY*I_ERI_Px_S_G2x2z_S_M1_vrr;
      Double I_ERI_Dxz_S_G2x2z_S_vrr = PAZ*I_ERI_Px_S_G2x2z_S_vrr+WPZ*I_ERI_Px_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_F2xz_S_M1_vrr;
      Double I_ERI_D2y_S_G2x2z_S_vrr = PAY*I_ERI_Py_S_G2x2z_S_vrr+WPY*I_ERI_Py_S_G2x2z_S_M1_vrr+oned2z*I_ERI_S_S_G2x2z_S_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M1_vrr;
      Double I_ERI_Dyz_S_G2x2z_S_vrr = PAZ*I_ERI_Py_S_G2x2z_S_vrr+WPZ*I_ERI_Py_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_F2xz_S_M1_vrr;
      Double I_ERI_D2z_S_G2x2z_S_vrr = PAZ*I_ERI_Pz_S_G2x2z_S_vrr+WPZ*I_ERI_Pz_S_G2x2z_S_M1_vrr+oned2z*I_ERI_S_S_G2x2z_S_vrr-rhod2zsq*I_ERI_S_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F2xz_S_M1_vrr;
      Double I_ERI_D2x_S_Gx3y_S_vrr = PAX*I_ERI_Px_S_Gx3y_S_vrr+WPX*I_ERI_Px_S_Gx3y_S_M1_vrr+oned2z*I_ERI_S_S_Gx3y_S_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M1_vrr+oned2k*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_Dxy_S_Gx3y_S_vrr = PAY*I_ERI_Px_S_Gx3y_S_vrr+WPY*I_ERI_Px_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Px_S_Fx2y_S_M1_vrr;
      Double I_ERI_Dxz_S_Gx3y_S_vrr = PAZ*I_ERI_Px_S_Gx3y_S_vrr+WPZ*I_ERI_Px_S_Gx3y_S_M1_vrr;
      Double I_ERI_D2y_S_Gx3y_S_vrr = PAY*I_ERI_Py_S_Gx3y_S_vrr+WPY*I_ERI_Py_S_Gx3y_S_M1_vrr+oned2z*I_ERI_S_S_Gx3y_S_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_Py_S_Fx2y_S_M1_vrr;
      Double I_ERI_Dyz_S_Gx3y_S_vrr = PAZ*I_ERI_Py_S_Gx3y_S_vrr+WPZ*I_ERI_Py_S_Gx3y_S_M1_vrr;
      Double I_ERI_D2z_S_Gx3y_S_vrr = PAZ*I_ERI_Pz_S_Gx3y_S_vrr+WPZ*I_ERI_Pz_S_Gx3y_S_M1_vrr+oned2z*I_ERI_S_S_Gx3y_S_vrr-rhod2zsq*I_ERI_S_S_Gx3y_S_M1_vrr;
      Double I_ERI_D2x_S_Gx2yz_S_vrr = PAX*I_ERI_Px_S_Gx2yz_S_vrr+WPX*I_ERI_Px_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_S_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Px_S_F2yz_S_M1_vrr;
      Double I_ERI_Dxy_S_Gx2yz_S_vrr = PAY*I_ERI_Px_S_Gx2yz_S_vrr+WPY*I_ERI_Px_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M1_vrr;
      Double I_ERI_Dxz_S_Gx2yz_S_vrr = PAZ*I_ERI_Px_S_Gx2yz_S_vrr+WPZ*I_ERI_Px_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Px_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2y_S_Gx2yz_S_vrr = PAY*I_ERI_Py_S_Gx2yz_S_vrr+WPY*I_ERI_Py_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_S_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_Py_S_Fxyz_S_M1_vrr;
      Double I_ERI_Dyz_S_Gx2yz_S_vrr = PAZ*I_ERI_Py_S_Gx2yz_S_vrr+WPZ*I_ERI_Py_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Py_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2z_S_Gx2yz_S_vrr = PAZ*I_ERI_Pz_S_Gx2yz_S_vrr+WPZ*I_ERI_Pz_S_Gx2yz_S_M1_vrr+oned2z*I_ERI_S_S_Gx2yz_S_vrr-rhod2zsq*I_ERI_S_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Pz_S_Fx2y_S_M1_vrr;
      Double I_ERI_D2x_S_Gxy2z_S_vrr = PAX*I_ERI_Px_S_Gxy2z_S_vrr+WPX*I_ERI_Px_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_S_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Px_S_Fy2z_S_M1_vrr;
      Double I_ERI_Dxy_S_Gxy2z_S_vrr = PAY*I_ERI_Px_S_Gxy2z_S_vrr+WPY*I_ERI_Px_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Px_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dxz_S_Gxy2z_S_vrr = PAZ*I_ERI_Px_S_Gxy2z_S_vrr+WPZ*I_ERI_Px_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fxyz_S_M1_vrr;
      Double I_ERI_D2y_S_Gxy2z_S_vrr = PAY*I_ERI_Py_S_Gxy2z_S_vrr+WPY*I_ERI_Py_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_S_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_Py_S_Fx2z_S_M1_vrr;
      Double I_ERI_Dyz_S_Gxy2z_S_vrr = PAZ*I_ERI_Py_S_Gxy2z_S_vrr+WPZ*I_ERI_Py_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Fxyz_S_M1_vrr;
      Double I_ERI_D2z_S_Gxy2z_S_vrr = PAZ*I_ERI_Pz_S_Gxy2z_S_vrr+WPZ*I_ERI_Pz_S_Gxy2z_S_M1_vrr+oned2z*I_ERI_S_S_Gxy2z_S_vrr-rhod2zsq*I_ERI_S_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Fxyz_S_M1_vrr;
      Double I_ERI_D2x_S_Gx3z_S_vrr = PAX*I_ERI_Px_S_Gx3z_S_vrr+WPX*I_ERI_Px_S_Gx3z_S_M1_vrr+oned2z*I_ERI_S_S_Gx3z_S_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M1_vrr+oned2k*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_Dxy_S_Gx3z_S_vrr = PAY*I_ERI_Px_S_Gx3z_S_vrr+WPY*I_ERI_Px_S_Gx3z_S_M1_vrr;
      Double I_ERI_Dxz_S_Gx3z_S_vrr = PAZ*I_ERI_Px_S_Gx3z_S_vrr+WPZ*I_ERI_Px_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Px_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2y_S_Gx3z_S_vrr = PAY*I_ERI_Py_S_Gx3z_S_vrr+WPY*I_ERI_Py_S_Gx3z_S_M1_vrr+oned2z*I_ERI_S_S_Gx3z_S_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M1_vrr;
      Double I_ERI_Dyz_S_Gx3z_S_vrr = PAZ*I_ERI_Py_S_Gx3z_S_vrr+WPZ*I_ERI_Py_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Py_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2z_S_Gx3z_S_vrr = PAZ*I_ERI_Pz_S_Gx3z_S_vrr+WPZ*I_ERI_Pz_S_Gx3z_S_M1_vrr+oned2z*I_ERI_S_S_Gx3z_S_vrr-rhod2zsq*I_ERI_S_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Pz_S_Fx2z_S_M1_vrr;
      Double I_ERI_D2x_S_G4y_S_vrr = PAX*I_ERI_Px_S_G4y_S_vrr+WPX*I_ERI_Px_S_G4y_S_M1_vrr+oned2z*I_ERI_S_S_G4y_S_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_Dxy_S_G4y_S_vrr = PAY*I_ERI_Px_S_G4y_S_vrr+WPY*I_ERI_Px_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_Dxz_S_G4y_S_vrr = PAZ*I_ERI_Px_S_G4y_S_vrr+WPZ*I_ERI_Px_S_G4y_S_M1_vrr;
      Double I_ERI_D2y_S_G4y_S_vrr = PAY*I_ERI_Py_S_G4y_S_vrr+WPY*I_ERI_Py_S_G4y_S_M1_vrr+oned2z*I_ERI_S_S_G4y_S_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M1_vrr+4*oned2k*I_ERI_Py_S_F3y_S_M1_vrr;
      Double I_ERI_Dyz_S_G4y_S_vrr = PAZ*I_ERI_Py_S_G4y_S_vrr+WPZ*I_ERI_Py_S_G4y_S_M1_vrr;
      Double I_ERI_D2z_S_G4y_S_vrr = PAZ*I_ERI_Pz_S_G4y_S_vrr+WPZ*I_ERI_Pz_S_G4y_S_M1_vrr+oned2z*I_ERI_S_S_G4y_S_vrr-rhod2zsq*I_ERI_S_S_G4y_S_M1_vrr;
      Double I_ERI_D2x_S_G3yz_S_vrr = PAX*I_ERI_Px_S_G3yz_S_vrr+WPX*I_ERI_Px_S_G3yz_S_M1_vrr+oned2z*I_ERI_S_S_G3yz_S_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M1_vrr;
      Double I_ERI_Dxy_S_G3yz_S_vrr = PAY*I_ERI_Px_S_G3yz_S_vrr+WPY*I_ERI_Px_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Px_S_F2yz_S_M1_vrr;
      Double I_ERI_Dxz_S_G3yz_S_vrr = PAZ*I_ERI_Px_S_G3yz_S_vrr+WPZ*I_ERI_Px_S_G3yz_S_M1_vrr+oned2k*I_ERI_Px_S_F3y_S_M1_vrr;
      Double I_ERI_D2y_S_G3yz_S_vrr = PAY*I_ERI_Py_S_G3yz_S_vrr+WPY*I_ERI_Py_S_G3yz_S_M1_vrr+oned2z*I_ERI_S_S_G3yz_S_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_Py_S_F2yz_S_M1_vrr;
      Double I_ERI_Dyz_S_G3yz_S_vrr = PAZ*I_ERI_Py_S_G3yz_S_vrr+WPZ*I_ERI_Py_S_G3yz_S_M1_vrr+oned2k*I_ERI_Py_S_F3y_S_M1_vrr;
      Double I_ERI_D2z_S_G3yz_S_vrr = PAZ*I_ERI_Pz_S_G3yz_S_vrr+WPZ*I_ERI_Pz_S_G3yz_S_M1_vrr+oned2z*I_ERI_S_S_G3yz_S_vrr-rhod2zsq*I_ERI_S_S_G3yz_S_M1_vrr+oned2k*I_ERI_Pz_S_F3y_S_M1_vrr;
      Double I_ERI_D2x_S_G2y2z_S_vrr = PAX*I_ERI_Px_S_G2y2z_S_vrr+WPX*I_ERI_Px_S_G2y2z_S_M1_vrr+oned2z*I_ERI_S_S_G2y2z_S_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M1_vrr;
      Double I_ERI_Dxy_S_G2y2z_S_vrr = PAY*I_ERI_Px_S_G2y2z_S_vrr+WPY*I_ERI_Px_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Fy2z_S_M1_vrr;
      Double I_ERI_Dxz_S_G2y2z_S_vrr = PAZ*I_ERI_Px_S_G2y2z_S_vrr+WPZ*I_ERI_Px_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_F2yz_S_M1_vrr;
      Double I_ERI_D2y_S_G2y2z_S_vrr = PAY*I_ERI_Py_S_G2y2z_S_vrr+WPY*I_ERI_Py_S_G2y2z_S_M1_vrr+oned2z*I_ERI_S_S_G2y2z_S_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Fy2z_S_M1_vrr;
      Double I_ERI_Dyz_S_G2y2z_S_vrr = PAZ*I_ERI_Py_S_G2y2z_S_vrr+WPZ*I_ERI_Py_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_F2yz_S_M1_vrr;
      Double I_ERI_D2z_S_G2y2z_S_vrr = PAZ*I_ERI_Pz_S_G2y2z_S_vrr+WPZ*I_ERI_Pz_S_G2y2z_S_M1_vrr+oned2z*I_ERI_S_S_G2y2z_S_vrr-rhod2zsq*I_ERI_S_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_F2yz_S_M1_vrr;
      Double I_ERI_D2x_S_Gy3z_S_vrr = PAX*I_ERI_Px_S_Gy3z_S_vrr+WPX*I_ERI_Px_S_Gy3z_S_M1_vrr+oned2z*I_ERI_S_S_Gy3z_S_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M1_vrr;
      Double I_ERI_Dxy_S_Gy3z_S_vrr = PAY*I_ERI_Px_S_Gy3z_S_vrr+WPY*I_ERI_Px_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_Dxz_S_Gy3z_S_vrr = PAZ*I_ERI_Px_S_Gy3z_S_vrr+WPZ*I_ERI_Px_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Px_S_Fy2z_S_M1_vrr;
      Double I_ERI_D2y_S_Gy3z_S_vrr = PAY*I_ERI_Py_S_Gy3z_S_vrr+WPY*I_ERI_Py_S_Gy3z_S_M1_vrr+oned2z*I_ERI_S_S_Gy3z_S_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M1_vrr+oned2k*I_ERI_Py_S_F3z_S_M1_vrr;
      Double I_ERI_Dyz_S_Gy3z_S_vrr = PAZ*I_ERI_Py_S_Gy3z_S_vrr+WPZ*I_ERI_Py_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Py_S_Fy2z_S_M1_vrr;
      Double I_ERI_D2z_S_Gy3z_S_vrr = PAZ*I_ERI_Pz_S_Gy3z_S_vrr+WPZ*I_ERI_Pz_S_Gy3z_S_M1_vrr+oned2z*I_ERI_S_S_Gy3z_S_vrr-rhod2zsq*I_ERI_S_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Pz_S_Fy2z_S_M1_vrr;
      Double I_ERI_D2x_S_G4z_S_vrr = PAX*I_ERI_Px_S_G4z_S_vrr+WPX*I_ERI_Px_S_G4z_S_M1_vrr+oned2z*I_ERI_S_S_G4z_S_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_Dxy_S_G4z_S_vrr = PAY*I_ERI_Px_S_G4z_S_vrr+WPY*I_ERI_Px_S_G4z_S_M1_vrr;
      Double I_ERI_Dxz_S_G4z_S_vrr = PAZ*I_ERI_Px_S_G4z_S_vrr+WPZ*I_ERI_Px_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Px_S_F3z_S_M1_vrr;
      Double I_ERI_D2y_S_G4z_S_vrr = PAY*I_ERI_Py_S_G4z_S_vrr+WPY*I_ERI_Py_S_G4z_S_M1_vrr+oned2z*I_ERI_S_S_G4z_S_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M1_vrr;
      Double I_ERI_Dyz_S_G4z_S_vrr = PAZ*I_ERI_Py_S_G4z_S_vrr+WPZ*I_ERI_Py_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Py_S_F3z_S_M1_vrr;
      Double I_ERI_D2z_S_G4z_S_vrr = PAZ*I_ERI_Pz_S_G4z_S_vrr+WPZ*I_ERI_Pz_S_G4z_S_M1_vrr+oned2z*I_ERI_S_S_G4z_S_vrr-rhod2zsq*I_ERI_S_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Pz_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
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
      Double I_ERI_Fxyz_S_G4x_S_vrr = PAZ*I_ERI_Dxy_S_G4x_S_vrr+WPZ*I_ERI_Dxy_S_G4x_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4x_S_vrr = PAX*I_ERI_D2z_S_G4x_S_vrr+WPX*I_ERI_D2z_S_G4x_S_M1_vrr+4*oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3y_S_G4x_S_vrr = PAY*I_ERI_D2y_S_G4x_S_vrr+WPY*I_ERI_D2y_S_G4x_S_M1_vrr+2*oned2z*I_ERI_Py_S_G4x_S_vrr-2*rhod2zsq*I_ERI_Py_S_G4x_S_M1_vrr;
      Double I_ERI_F2yz_S_G4x_S_vrr = PAZ*I_ERI_D2y_S_G4x_S_vrr+WPZ*I_ERI_D2y_S_G4x_S_M1_vrr;
      Double I_ERI_Fy2z_S_G4x_S_vrr = PAY*I_ERI_D2z_S_G4x_S_vrr+WPY*I_ERI_D2z_S_G4x_S_M1_vrr;
      Double I_ERI_F3z_S_G4x_S_vrr = PAZ*I_ERI_D2z_S_G4x_S_vrr+WPZ*I_ERI_D2z_S_G4x_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G4x_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G4x_S_M1_vrr;
      Double I_ERI_F3x_S_G3xy_S_vrr = PAX*I_ERI_D2x_S_G3xy_S_vrr+WPX*I_ERI_D2x_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_Px_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_Px_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_F2xy_S_G3xy_S_vrr = PAY*I_ERI_D2x_S_G3xy_S_vrr+WPY*I_ERI_D2x_S_G3xy_S_M1_vrr+oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_F2xz_S_G3xy_S_vrr = PAZ*I_ERI_D2x_S_G3xy_S_vrr+WPZ*I_ERI_D2x_S_G3xy_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3xy_S_vrr = PAX*I_ERI_D2y_S_G3xy_S_vrr+WPX*I_ERI_D2y_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Fxyz_S_G3xy_S_vrr = PAZ*I_ERI_Dxy_S_G3xy_S_vrr+WPZ*I_ERI_Dxy_S_G3xy_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3xy_S_vrr = PAX*I_ERI_D2z_S_G3xy_S_vrr+WPX*I_ERI_D2z_S_G3xy_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_F3y_S_G3xy_S_vrr = PAY*I_ERI_D2y_S_G3xy_S_vrr+WPY*I_ERI_D2y_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_Py_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_Py_S_G3xy_S_M1_vrr+oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_F2yz_S_G3xy_S_vrr = PAZ*I_ERI_D2y_S_G3xy_S_vrr+WPZ*I_ERI_D2y_S_G3xy_S_M1_vrr;
      Double I_ERI_Fy2z_S_G3xy_S_vrr = PAY*I_ERI_D2z_S_G3xy_S_vrr+WPY*I_ERI_D2z_S_G3xy_S_M1_vrr+oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3z_S_G3xy_S_vrr = PAZ*I_ERI_D2z_S_G3xy_S_vrr+WPZ*I_ERI_D2z_S_G3xy_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G3xy_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xy_S_M1_vrr;
      Double I_ERI_F3x_S_G3xz_S_vrr = PAX*I_ERI_D2x_S_G3xz_S_vrr+WPX*I_ERI_D2x_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_Px_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_Px_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_F2xy_S_G3xz_S_vrr = PAY*I_ERI_D2x_S_G3xz_S_vrr+WPY*I_ERI_D2x_S_G3xz_S_M1_vrr;
      Double I_ERI_F2xz_S_G3xz_S_vrr = PAZ*I_ERI_D2x_S_G3xz_S_vrr+WPZ*I_ERI_D2x_S_G3xz_S_M1_vrr+oned2k*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3xz_S_vrr = PAX*I_ERI_D2y_S_G3xz_S_vrr+WPX*I_ERI_D2y_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_Fxyz_S_G3xz_S_vrr = PAZ*I_ERI_Dxy_S_G3xz_S_vrr+WPZ*I_ERI_Dxy_S_G3xz_S_M1_vrr+oned2k*I_ERI_Dxy_S_F3x_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3xz_S_vrr = PAX*I_ERI_D2z_S_G3xz_S_vrr+WPX*I_ERI_D2z_S_G3xz_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3y_S_G3xz_S_vrr = PAY*I_ERI_D2y_S_G3xz_S_vrr+WPY*I_ERI_D2y_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_Py_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_Py_S_G3xz_S_M1_vrr;
      Double I_ERI_F2yz_S_G3xz_S_vrr = PAZ*I_ERI_D2y_S_G3xz_S_vrr+WPZ*I_ERI_D2y_S_G3xz_S_M1_vrr+oned2k*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Fy2z_S_G3xz_S_vrr = PAY*I_ERI_D2z_S_G3xz_S_vrr+WPY*I_ERI_D2z_S_G3xz_S_M1_vrr;
      Double I_ERI_F3z_S_G3xz_S_vrr = PAZ*I_ERI_D2z_S_G3xz_S_vrr+WPZ*I_ERI_D2z_S_G3xz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G3xz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G3xz_S_M1_vrr+oned2k*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_F3x_S_G2x2y_S_vrr = PAX*I_ERI_D2x_S_G2x2y_S_vrr+WPX*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2xy_S_G2x2y_S_vrr = PAY*I_ERI_D2x_S_G2x2y_S_vrr+WPY*I_ERI_D2x_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_F2xz_S_G2x2y_S_vrr = PAZ*I_ERI_D2x_S_G2x2y_S_vrr+WPZ*I_ERI_D2x_S_G2x2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2x2y_S_vrr = PAX*I_ERI_D2y_S_G2x2y_S_vrr+WPX*I_ERI_D2y_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_G2x2y_S_vrr = PAZ*I_ERI_Dxy_S_G2x2y_S_vrr+WPZ*I_ERI_Dxy_S_G2x2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2x2y_S_vrr = PAX*I_ERI_D2z_S_G2x2y_S_vrr+WPX*I_ERI_D2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3y_S_G2x2y_S_vrr = PAY*I_ERI_D2y_S_G2x2y_S_vrr+WPY*I_ERI_D2y_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_F2yz_S_G2x2y_S_vrr = PAZ*I_ERI_D2y_S_G2x2y_S_vrr+WPZ*I_ERI_D2y_S_G2x2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_G2x2y_S_vrr = PAY*I_ERI_D2z_S_G2x2y_S_vrr+WPY*I_ERI_D2z_S_G2x2y_S_M1_vrr+2*oned2k*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_F3z_S_G2x2y_S_vrr = PAZ*I_ERI_D2z_S_G2x2y_S_vrr+WPZ*I_ERI_D2z_S_G2x2y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2x2y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2y_S_M1_vrr;
      Double I_ERI_F3x_S_G2xyz_S_vrr = PAX*I_ERI_D2x_S_G2xyz_S_vrr+WPX*I_ERI_D2x_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M1_vrr;
      Double I_ERI_F2xy_S_G2xyz_S_vrr = PAY*I_ERI_D2x_S_G2xyz_S_vrr+WPY*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_F2xz_S_G2xyz_S_vrr = PAZ*I_ERI_D2x_S_G2xyz_S_vrr+WPZ*I_ERI_D2x_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2xyz_S_vrr = PAX*I_ERI_D2y_S_G2xyz_S_vrr+WPX*I_ERI_D2y_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Fxyz_S_G2xyz_S_vrr = PAZ*I_ERI_Dxy_S_G2xyz_S_vrr+WPZ*I_ERI_Dxy_S_G2xyz_S_M1_vrr+oned2k*I_ERI_Dxy_S_F2xy_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2xyz_S_vrr = PAX*I_ERI_D2z_S_G2xyz_S_vrr+WPX*I_ERI_D2z_S_G2xyz_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_F3y_S_G2xyz_S_vrr = PAY*I_ERI_D2y_S_G2xyz_S_vrr+WPY*I_ERI_D2y_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_F2yz_S_G2xyz_S_vrr = PAZ*I_ERI_D2y_S_G2xyz_S_vrr+WPZ*I_ERI_D2y_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Fy2z_S_G2xyz_S_vrr = PAY*I_ERI_D2z_S_G2xyz_S_vrr+WPY*I_ERI_D2z_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3z_S_G2xyz_S_vrr = PAZ*I_ERI_D2z_S_G2xyz_S_vrr+WPZ*I_ERI_D2z_S_G2xyz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2xyz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2xyz_S_M1_vrr+oned2k*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_F3x_S_G2x2z_S_vrr = PAX*I_ERI_D2x_S_G2x2z_S_vrr+WPX*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xy_S_G2x2z_S_vrr = PAY*I_ERI_D2x_S_G2x2z_S_vrr+WPY*I_ERI_D2x_S_G2x2z_S_M1_vrr;
      Double I_ERI_F2xz_S_G2x2z_S_vrr = PAZ*I_ERI_D2x_S_G2x2z_S_vrr+WPZ*I_ERI_D2x_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2x2z_S_vrr = PAX*I_ERI_D2y_S_G2x2z_S_vrr+WPX*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_G2x2z_S_vrr = PAZ*I_ERI_Dxy_S_G2x2z_S_vrr+WPZ*I_ERI_Dxy_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_F2xz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2x2z_S_vrr = PAX*I_ERI_D2z_S_G2x2z_S_vrr+WPX*I_ERI_D2z_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3y_S_G2x2z_S_vrr = PAY*I_ERI_D2y_S_G2x2z_S_vrr+WPY*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2x2z_S_M1_vrr;
      Double I_ERI_F2yz_S_G2x2z_S_vrr = PAZ*I_ERI_D2y_S_G2x2z_S_vrr+WPZ*I_ERI_D2y_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_Fy2z_S_G2x2z_S_vrr = PAY*I_ERI_D2z_S_G2x2z_S_vrr+WPY*I_ERI_D2z_S_G2x2z_S_M1_vrr;
      Double I_ERI_F3z_S_G2x2z_S_vrr = PAZ*I_ERI_D2z_S_G2x2z_S_vrr+WPZ*I_ERI_D2z_S_G2x2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2x2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2x2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_F2xz_S_M1_vrr;
      Double I_ERI_F3x_S_Gx3y_S_vrr = PAX*I_ERI_D2x_S_Gx3y_S_vrr+WPX*I_ERI_D2x_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3y_S_M1_vrr+oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx3y_S_vrr = PAY*I_ERI_D2x_S_Gx3y_S_vrr+WPY*I_ERI_D2x_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx3y_S_vrr = PAZ*I_ERI_D2x_S_Gx3y_S_vrr+WPZ*I_ERI_D2x_S_Gx3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx3y_S_vrr = PAX*I_ERI_D2y_S_Gx3y_S_vrr+WPX*I_ERI_D2y_S_Gx3y_S_M1_vrr+oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gx3y_S_vrr = PAZ*I_ERI_Dxy_S_Gx3y_S_vrr+WPZ*I_ERI_Dxy_S_Gx3y_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx3y_S_vrr = PAX*I_ERI_D2z_S_Gx3y_S_vrr+WPX*I_ERI_D2z_S_Gx3y_S_M1_vrr+oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3y_S_Gx3y_S_vrr = PAY*I_ERI_D2y_S_Gx3y_S_vrr+WPY*I_ERI_D2y_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx3y_S_vrr = PAZ*I_ERI_D2y_S_Gx3y_S_vrr+WPZ*I_ERI_D2y_S_Gx3y_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gx3y_S_vrr = PAY*I_ERI_D2z_S_Gx3y_S_vrr+WPY*I_ERI_D2z_S_Gx3y_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3z_S_Gx3y_S_vrr = PAZ*I_ERI_D2z_S_Gx3y_S_vrr+WPZ*I_ERI_D2z_S_Gx3y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gx3y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3y_S_M1_vrr;
      Double I_ERI_F3x_S_Gx2yz_S_vrr = PAX*I_ERI_D2x_S_Gx2yz_S_vrr+WPX*I_ERI_D2x_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx2yz_S_vrr = PAY*I_ERI_D2x_S_Gx2yz_S_vrr+WPY*I_ERI_D2x_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx2yz_S_vrr = PAZ*I_ERI_D2x_S_Gx2yz_S_vrr+WPZ*I_ERI_D2x_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx2yz_S_vrr = PAX*I_ERI_D2y_S_Gx2yz_S_vrr+WPX*I_ERI_D2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gx2yz_S_vrr = PAZ*I_ERI_Dxy_S_Gx2yz_S_vrr+WPZ*I_ERI_Dxy_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_Dxy_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx2yz_S_vrr = PAX*I_ERI_D2z_S_Gx2yz_S_vrr+WPX*I_ERI_D2z_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3y_S_Gx2yz_S_vrr = PAY*I_ERI_D2y_S_Gx2yz_S_vrr+WPY*I_ERI_D2y_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx2yz_S_vrr = PAZ*I_ERI_D2y_S_Gx2yz_S_vrr+WPZ*I_ERI_D2y_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gx2yz_S_vrr = PAY*I_ERI_D2z_S_Gx2yz_S_vrr+WPY*I_ERI_D2z_S_Gx2yz_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_F3z_S_Gx2yz_S_vrr = PAZ*I_ERI_D2z_S_Gx2yz_S_vrr+WPZ*I_ERI_D2z_S_Gx2yz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gx2yz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx2yz_S_M1_vrr+oned2k*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_F3x_S_Gxy2z_S_vrr = PAX*I_ERI_D2x_S_Gxy2z_S_vrr+WPX*I_ERI_D2x_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gxy2z_S_vrr = PAY*I_ERI_D2x_S_Gxy2z_S_vrr+WPY*I_ERI_D2x_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gxy2z_S_vrr = PAZ*I_ERI_D2x_S_Gxy2z_S_vrr+WPZ*I_ERI_D2x_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fxyz_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gxy2z_S_vrr = PAX*I_ERI_D2y_S_Gxy2z_S_vrr+WPX*I_ERI_D2y_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gxy2z_S_vrr = PAZ*I_ERI_Dxy_S_Gxy2z_S_vrr+WPZ*I_ERI_Dxy_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Fxyz_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gxy2z_S_vrr = PAX*I_ERI_D2z_S_Gxy2z_S_vrr+WPX*I_ERI_D2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3y_S_Gxy2z_S_vrr = PAY*I_ERI_D2y_S_Gxy2z_S_vrr+WPY*I_ERI_D2y_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gxy2z_S_vrr = PAZ*I_ERI_D2y_S_Gxy2z_S_vrr+WPZ*I_ERI_D2y_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fxyz_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gxy2z_S_vrr = PAY*I_ERI_D2z_S_Gxy2z_S_vrr+WPY*I_ERI_D2z_S_Gxy2z_S_M1_vrr+oned2k*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3z_S_Gxy2z_S_vrr = PAZ*I_ERI_D2z_S_Gxy2z_S_vrr+WPZ*I_ERI_D2z_S_Gxy2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gxy2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gxy2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fxyz_S_M1_vrr;
      Double I_ERI_F3x_S_Gx3z_S_vrr = PAX*I_ERI_D2x_S_Gx3z_S_vrr+WPX*I_ERI_D2x_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gx3z_S_M1_vrr+oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gx3z_S_vrr = PAY*I_ERI_D2x_S_Gx3z_S_vrr+WPY*I_ERI_D2x_S_Gx3z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gx3z_S_vrr = PAZ*I_ERI_D2x_S_Gx3z_S_vrr+WPZ*I_ERI_D2x_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gx3z_S_vrr = PAX*I_ERI_D2y_S_Gx3z_S_vrr+WPX*I_ERI_D2y_S_Gx3z_S_M1_vrr+oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gx3z_S_vrr = PAZ*I_ERI_Dxy_S_Gx3z_S_vrr+WPZ*I_ERI_Dxy_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_Dxy_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gx3z_S_vrr = PAX*I_ERI_D2z_S_Gx3z_S_vrr+WPX*I_ERI_D2z_S_Gx3z_S_M1_vrr+oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3y_S_Gx3z_S_vrr = PAY*I_ERI_D2y_S_Gx3z_S_vrr+WPY*I_ERI_D2y_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gx3z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gx3z_S_vrr = PAZ*I_ERI_D2y_S_Gx3z_S_vrr+WPZ*I_ERI_D2y_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gx3z_S_vrr = PAY*I_ERI_D2z_S_Gx3z_S_vrr+WPY*I_ERI_D2z_S_Gx3z_S_M1_vrr;
      Double I_ERI_F3z_S_Gx3z_S_vrr = PAZ*I_ERI_D2z_S_Gx3z_S_vrr+WPZ*I_ERI_D2z_S_Gx3z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gx3z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gx3z_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_F3x_S_G4y_S_vrr = PAX*I_ERI_D2x_S_G4y_S_vrr+WPX*I_ERI_D2x_S_G4y_S_M1_vrr+2*oned2z*I_ERI_Px_S_G4y_S_vrr-2*rhod2zsq*I_ERI_Px_S_G4y_S_M1_vrr;
      Double I_ERI_F2xy_S_G4y_S_vrr = PAY*I_ERI_D2x_S_G4y_S_vrr+WPY*I_ERI_D2x_S_G4y_S_M1_vrr+4*oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_F2xz_S_G4y_S_vrr = PAZ*I_ERI_D2x_S_G4y_S_vrr+WPZ*I_ERI_D2x_S_G4y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4y_S_vrr = PAX*I_ERI_D2y_S_G4y_S_vrr+WPX*I_ERI_D2y_S_G4y_S_M1_vrr;
      Double I_ERI_Fxyz_S_G4y_S_vrr = PAZ*I_ERI_Dxy_S_G4y_S_vrr+WPZ*I_ERI_Dxy_S_G4y_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4y_S_vrr = PAX*I_ERI_D2z_S_G4y_S_vrr+WPX*I_ERI_D2z_S_G4y_S_M1_vrr;
      Double I_ERI_F3y_S_G4y_S_vrr = PAY*I_ERI_D2y_S_G4y_S_vrr+WPY*I_ERI_D2y_S_G4y_S_M1_vrr+2*oned2z*I_ERI_Py_S_G4y_S_vrr-2*rhod2zsq*I_ERI_Py_S_G4y_S_M1_vrr+4*oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_F2yz_S_G4y_S_vrr = PAZ*I_ERI_D2y_S_G4y_S_vrr+WPZ*I_ERI_D2y_S_G4y_S_M1_vrr;
      Double I_ERI_Fy2z_S_G4y_S_vrr = PAY*I_ERI_D2z_S_G4y_S_vrr+WPY*I_ERI_D2z_S_G4y_S_M1_vrr+4*oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3z_S_G4y_S_vrr = PAZ*I_ERI_D2z_S_G4y_S_vrr+WPZ*I_ERI_D2z_S_G4y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G4y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G4y_S_M1_vrr;
      Double I_ERI_F3x_S_G3yz_S_vrr = PAX*I_ERI_D2x_S_G3yz_S_vrr+WPX*I_ERI_D2x_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_Px_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_Px_S_G3yz_S_M1_vrr;
      Double I_ERI_F2xy_S_G3yz_S_vrr = PAY*I_ERI_D2x_S_G3yz_S_vrr+WPY*I_ERI_D2x_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_F2xz_S_G3yz_S_vrr = PAZ*I_ERI_D2x_S_G3yz_S_vrr+WPZ*I_ERI_D2x_S_G3yz_S_M1_vrr+oned2k*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2y_S_G3yz_S_vrr = PAX*I_ERI_D2y_S_G3yz_S_vrr+WPX*I_ERI_D2y_S_G3yz_S_M1_vrr;
      Double I_ERI_Fxyz_S_G3yz_S_vrr = PAZ*I_ERI_Dxy_S_G3yz_S_vrr+WPZ*I_ERI_Dxy_S_G3yz_S_M1_vrr+oned2k*I_ERI_Dxy_S_F3y_S_M1_vrr;
      Double I_ERI_Fx2z_S_G3yz_S_vrr = PAX*I_ERI_D2z_S_G3yz_S_vrr+WPX*I_ERI_D2z_S_G3yz_S_M1_vrr;
      Double I_ERI_F3y_S_G3yz_S_vrr = PAY*I_ERI_D2y_S_G3yz_S_vrr+WPY*I_ERI_D2y_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_Py_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_Py_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_F2yz_S_G3yz_S_vrr = PAZ*I_ERI_D2y_S_G3yz_S_vrr+WPZ*I_ERI_D2y_S_G3yz_S_M1_vrr+oned2k*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Fy2z_S_G3yz_S_vrr = PAY*I_ERI_D2z_S_G3yz_S_vrr+WPY*I_ERI_D2z_S_G3yz_S_M1_vrr+3*oned2k*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3z_S_G3yz_S_vrr = PAZ*I_ERI_D2z_S_G3yz_S_vrr+WPZ*I_ERI_D2z_S_G3yz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G3yz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G3yz_S_M1_vrr+oned2k*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_F3x_S_G2y2z_S_vrr = PAX*I_ERI_D2x_S_G2y2z_S_vrr+WPX*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_G2y2z_S_M1_vrr;
      Double I_ERI_F2xy_S_G2y2z_S_vrr = PAY*I_ERI_D2x_S_G2y2z_S_vrr+WPY*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2xz_S_G2y2z_S_vrr = PAZ*I_ERI_D2x_S_G2y2z_S_vrr+WPZ*I_ERI_D2x_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_Fx2y_S_G2y2z_S_vrr = PAX*I_ERI_D2y_S_G2y2z_S_vrr+WPX*I_ERI_D2y_S_G2y2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_G2y2z_S_vrr = PAZ*I_ERI_Dxy_S_G2y2z_S_vrr+WPZ*I_ERI_Dxy_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_F2yz_S_M1_vrr;
      Double I_ERI_Fx2z_S_G2y2z_S_vrr = PAX*I_ERI_D2z_S_G2y2z_S_vrr+WPX*I_ERI_D2z_S_G2y2z_S_M1_vrr;
      Double I_ERI_F3y_S_G2y2z_S_vrr = PAY*I_ERI_D2y_S_G2y2z_S_vrr+WPY*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_F2yz_S_G2y2z_S_vrr = PAZ*I_ERI_D2y_S_G2y2z_S_vrr+WPZ*I_ERI_D2y_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_F2yz_S_M1_vrr;
      Double I_ERI_Fy2z_S_G2y2z_S_vrr = PAY*I_ERI_D2z_S_G2y2z_S_vrr+WPY*I_ERI_D2z_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3z_S_G2y2z_S_vrr = PAZ*I_ERI_D2z_S_G2y2z_S_vrr+WPZ*I_ERI_D2z_S_G2y2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G2y2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G2y2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_F2yz_S_M1_vrr;
      Double I_ERI_F3x_S_Gy3z_S_vrr = PAX*I_ERI_D2x_S_Gy3z_S_vrr+WPX*I_ERI_D2x_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_Px_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_Px_S_Gy3z_S_M1_vrr;
      Double I_ERI_F2xy_S_Gy3z_S_vrr = PAY*I_ERI_D2x_S_Gy3z_S_vrr+WPY*I_ERI_D2x_S_Gy3z_S_M1_vrr+oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_F2xz_S_Gy3z_S_vrr = PAZ*I_ERI_D2x_S_Gy3z_S_vrr+WPZ*I_ERI_D2x_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fx2y_S_Gy3z_S_vrr = PAX*I_ERI_D2y_S_Gy3z_S_vrr+WPX*I_ERI_D2y_S_Gy3z_S_M1_vrr;
      Double I_ERI_Fxyz_S_Gy3z_S_vrr = PAZ*I_ERI_Dxy_S_Gy3z_S_vrr+WPZ*I_ERI_Dxy_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_Dxy_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fx2z_S_Gy3z_S_vrr = PAX*I_ERI_D2z_S_Gy3z_S_vrr+WPX*I_ERI_D2z_S_Gy3z_S_M1_vrr;
      Double I_ERI_F3y_S_Gy3z_S_vrr = PAY*I_ERI_D2y_S_Gy3z_S_vrr+WPY*I_ERI_D2y_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_Py_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_Py_S_Gy3z_S_M1_vrr+oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_F2yz_S_Gy3z_S_vrr = PAZ*I_ERI_D2y_S_Gy3z_S_vrr+WPZ*I_ERI_D2y_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_D2y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Fy2z_S_Gy3z_S_vrr = PAY*I_ERI_D2z_S_Gy3z_S_vrr+WPY*I_ERI_D2z_S_Gy3z_S_M1_vrr+oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;
      Double I_ERI_F3z_S_Gy3z_S_vrr = PAZ*I_ERI_D2z_S_Gy3z_S_vrr+WPZ*I_ERI_D2z_S_Gy3z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Gy3z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Gy3z_S_M1_vrr+3*oned2k*I_ERI_D2z_S_Fy2z_S_M1_vrr;
      Double I_ERI_F3x_S_G4z_S_vrr = PAX*I_ERI_D2x_S_G4z_S_vrr+WPX*I_ERI_D2x_S_G4z_S_M1_vrr+2*oned2z*I_ERI_Px_S_G4z_S_vrr-2*rhod2zsq*I_ERI_Px_S_G4z_S_M1_vrr;
      Double I_ERI_F2xy_S_G4z_S_vrr = PAY*I_ERI_D2x_S_G4z_S_vrr+WPY*I_ERI_D2x_S_G4z_S_M1_vrr;
      Double I_ERI_F2xz_S_G4z_S_vrr = PAZ*I_ERI_D2x_S_G4z_S_vrr+WPZ*I_ERI_D2x_S_G4z_S_M1_vrr+4*oned2k*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2y_S_G4z_S_vrr = PAX*I_ERI_D2y_S_G4z_S_vrr+WPX*I_ERI_D2y_S_G4z_S_M1_vrr;
      Double I_ERI_Fxyz_S_G4z_S_vrr = PAZ*I_ERI_Dxy_S_G4z_S_vrr+WPZ*I_ERI_Dxy_S_G4z_S_M1_vrr+4*oned2k*I_ERI_Dxy_S_F3z_S_M1_vrr;
      Double I_ERI_Fx2z_S_G4z_S_vrr = PAX*I_ERI_D2z_S_G4z_S_vrr+WPX*I_ERI_D2z_S_G4z_S_M1_vrr;
      Double I_ERI_F3y_S_G4z_S_vrr = PAY*I_ERI_D2y_S_G4z_S_vrr+WPY*I_ERI_D2y_S_G4z_S_M1_vrr+2*oned2z*I_ERI_Py_S_G4z_S_vrr-2*rhod2zsq*I_ERI_Py_S_G4z_S_M1_vrr;
      Double I_ERI_F2yz_S_G4z_S_vrr = PAZ*I_ERI_D2y_S_G4z_S_vrr+WPZ*I_ERI_D2y_S_G4z_S_M1_vrr+4*oned2k*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_Fy2z_S_G4z_S_vrr = PAY*I_ERI_D2z_S_G4z_S_vrr+WPY*I_ERI_D2z_S_G4z_S_M1_vrr;
      Double I_ERI_F3z_S_G4z_S_vrr = PAZ*I_ERI_D2z_S_G4z_S_vrr+WPZ*I_ERI_D2z_S_G4z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_G4z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_G4z_S_M1_vrr+4*oned2k*I_ERI_D2z_S_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_F_S
       * RHS shell quartet name: SQ_ERI_F_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_F_S
       * RHS shell quartet name: SQ_ERI_D_S_F_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_D_S_M1
       ************************************************************/
      Double I_ERI_G4x_S_F3x_S_vrr = PAX*I_ERI_F3x_S_F3x_S_vrr+WPX*I_ERI_F3x_S_F3x_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F3x_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F3x_S_M1_vrr+3*oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xy_S_F3x_S_vrr = PAY*I_ERI_F3x_S_F3x_S_vrr+WPY*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G3xz_S_F3x_S_vrr = PAZ*I_ERI_F3x_S_F3x_S_vrr+WPZ*I_ERI_F3x_S_F3x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3x_S_vrr = PAY*I_ERI_F2xy_S_F3x_S_vrr+WPY*I_ERI_F2xy_S_F3x_S_M1_vrr+oned2z*I_ERI_D2x_S_F3x_S_vrr-rhod2zsq*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3x_S_vrr = PAZ*I_ERI_F2xy_S_F3x_S_vrr+WPZ*I_ERI_F2xy_S_F3x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3x_S_vrr = PAZ*I_ERI_F2xz_S_F3x_S_vrr+WPZ*I_ERI_F2xz_S_F3x_S_M1_vrr+oned2z*I_ERI_D2x_S_F3x_S_vrr-rhod2zsq*I_ERI_D2x_S_F3x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3x_S_vrr = PAX*I_ERI_F3y_S_F3x_S_vrr+WPX*I_ERI_F3y_S_F3x_S_M1_vrr+3*oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3x_S_vrr = PAZ*I_ERI_Fx2y_S_F3x_S_vrr+WPZ*I_ERI_Fx2y_S_F3x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3x_S_vrr = PAY*I_ERI_Fx2z_S_F3x_S_vrr+WPY*I_ERI_Fx2z_S_F3x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3x_S_vrr = PAX*I_ERI_F3z_S_F3x_S_vrr+WPX*I_ERI_F3z_S_F3x_S_M1_vrr+3*oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4y_S_F3x_S_vrr = PAY*I_ERI_F3y_S_F3x_S_vrr+WPY*I_ERI_F3y_S_F3x_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F3x_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_G3yz_S_F3x_S_vrr = PAZ*I_ERI_F3y_S_F3x_S_vrr+WPZ*I_ERI_F3y_S_F3x_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3x_S_vrr = PAZ*I_ERI_F2yz_S_F3x_S_vrr+WPZ*I_ERI_F2yz_S_F3x_S_M1_vrr+oned2z*I_ERI_D2y_S_F3x_S_vrr-rhod2zsq*I_ERI_D2y_S_F3x_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3x_S_vrr = PAY*I_ERI_F3z_S_F3x_S_vrr+WPY*I_ERI_F3z_S_F3x_S_M1_vrr;
      Double I_ERI_G4z_S_F3x_S_vrr = PAZ*I_ERI_F3z_S_F3x_S_vrr+WPZ*I_ERI_F3z_S_F3x_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F3x_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F3x_S_M1_vrr;
      Double I_ERI_G4x_S_F2xy_S_vrr = PAX*I_ERI_F3x_S_F2xy_S_vrr+WPX*I_ERI_F3x_S_F2xy_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F2xy_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_G3xy_S_F2xy_S_vrr = PAY*I_ERI_F3x_S_F2xy_S_vrr+WPY*I_ERI_F3x_S_F2xy_S_M1_vrr+oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G3xz_S_F2xy_S_vrr = PAZ*I_ERI_F3x_S_F2xy_S_vrr+WPZ*I_ERI_F3x_S_F2xy_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2xy_S_vrr = PAY*I_ERI_F2xy_S_F2xy_S_vrr+WPY*I_ERI_F2xy_S_F2xy_S_M1_vrr+oned2z*I_ERI_D2x_S_F2xy_S_vrr-rhod2zsq*I_ERI_D2x_S_F2xy_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2xy_S_vrr = PAZ*I_ERI_F2xy_S_F2xy_S_vrr+WPZ*I_ERI_F2xy_S_F2xy_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2xy_S_vrr = PAZ*I_ERI_F2xz_S_F2xy_S_vrr+WPZ*I_ERI_F2xz_S_F2xy_S_M1_vrr+oned2z*I_ERI_D2x_S_F2xy_S_vrr-rhod2zsq*I_ERI_D2x_S_F2xy_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2xy_S_vrr = PAX*I_ERI_F3y_S_F2xy_S_vrr+WPX*I_ERI_F3y_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2xy_S_vrr = PAZ*I_ERI_Fx2y_S_F2xy_S_vrr+WPZ*I_ERI_Fx2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2xy_S_vrr = PAY*I_ERI_Fx2z_S_F2xy_S_vrr+WPY*I_ERI_Fx2z_S_F2xy_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2xy_S_vrr = PAX*I_ERI_F3z_S_F2xy_S_vrr+WPX*I_ERI_F3z_S_F2xy_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4y_S_F2xy_S_vrr = PAY*I_ERI_F3y_S_F2xy_S_vrr+WPY*I_ERI_F3y_S_F2xy_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F2xy_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F2xy_S_M1_vrr+oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_G3yz_S_F2xy_S_vrr = PAZ*I_ERI_F3y_S_F2xy_S_vrr+WPZ*I_ERI_F3y_S_F2xy_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2xy_S_vrr = PAZ*I_ERI_F2yz_S_F2xy_S_vrr+WPZ*I_ERI_F2yz_S_F2xy_S_M1_vrr+oned2z*I_ERI_D2y_S_F2xy_S_vrr-rhod2zsq*I_ERI_D2y_S_F2xy_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2xy_S_vrr = PAY*I_ERI_F3z_S_F2xy_S_vrr+WPY*I_ERI_F3z_S_F2xy_S_M1_vrr+oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4z_S_F2xy_S_vrr = PAZ*I_ERI_F3z_S_F2xy_S_vrr+WPZ*I_ERI_F3z_S_F2xy_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F2xy_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F2xy_S_M1_vrr;
      Double I_ERI_G4x_S_F2xz_S_vrr = PAX*I_ERI_F3x_S_F2xz_S_vrr+WPX*I_ERI_F3x_S_F2xz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F2xz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dxz_S_M1_vrr;
      Double I_ERI_G3xy_S_F2xz_S_vrr = PAY*I_ERI_F3x_S_F2xz_S_vrr+WPY*I_ERI_F3x_S_F2xz_S_M1_vrr;
      Double I_ERI_G3xz_S_F2xz_S_vrr = PAZ*I_ERI_F3x_S_F2xz_S_vrr+WPZ*I_ERI_F3x_S_F2xz_S_M1_vrr+oned2k*I_ERI_F3x_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2xz_S_vrr = PAY*I_ERI_F2xy_S_F2xz_S_vrr+WPY*I_ERI_F2xy_S_F2xz_S_M1_vrr+oned2z*I_ERI_D2x_S_F2xz_S_vrr-rhod2zsq*I_ERI_D2x_S_F2xz_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2xz_S_vrr = PAZ*I_ERI_F2xy_S_F2xz_S_vrr+WPZ*I_ERI_F2xy_S_F2xz_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2x_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2xz_S_vrr = PAZ*I_ERI_F2xz_S_F2xz_S_vrr+WPZ*I_ERI_F2xz_S_F2xz_S_M1_vrr+oned2z*I_ERI_D2x_S_F2xz_S_vrr-rhod2zsq*I_ERI_D2x_S_F2xz_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2x_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2xz_S_vrr = PAX*I_ERI_F3y_S_F2xz_S_vrr+WPX*I_ERI_F3y_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dxz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2xz_S_vrr = PAZ*I_ERI_Fx2y_S_F2xz_S_vrr+WPZ*I_ERI_Fx2y_S_F2xz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2x_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2xz_S_vrr = PAY*I_ERI_Fx2z_S_F2xz_S_vrr+WPY*I_ERI_Fx2z_S_F2xz_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2xz_S_vrr = PAX*I_ERI_F3z_S_F2xz_S_vrr+WPX*I_ERI_F3z_S_F2xz_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dxz_S_M1_vrr;
      Double I_ERI_G4y_S_F2xz_S_vrr = PAY*I_ERI_F3y_S_F2xz_S_vrr+WPY*I_ERI_F3y_S_F2xz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F2xz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F2xz_S_M1_vrr;
      Double I_ERI_G3yz_S_F2xz_S_vrr = PAZ*I_ERI_F3y_S_F2xz_S_vrr+WPZ*I_ERI_F3y_S_F2xz_S_M1_vrr+oned2k*I_ERI_F3y_S_D2x_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2xz_S_vrr = PAZ*I_ERI_F2yz_S_F2xz_S_vrr+WPZ*I_ERI_F2yz_S_F2xz_S_M1_vrr+oned2z*I_ERI_D2y_S_F2xz_S_vrr-rhod2zsq*I_ERI_D2y_S_F2xz_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2x_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2xz_S_vrr = PAY*I_ERI_F3z_S_F2xz_S_vrr+WPY*I_ERI_F3z_S_F2xz_S_M1_vrr;
      Double I_ERI_G4z_S_F2xz_S_vrr = PAZ*I_ERI_F3z_S_F2xz_S_vrr+WPZ*I_ERI_F3z_S_F2xz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F2xz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F2xz_S_M1_vrr+oned2k*I_ERI_F3z_S_D2x_S_M1_vrr;
      Double I_ERI_G4x_S_Fx2y_S_vrr = PAX*I_ERI_F3x_S_Fx2y_S_vrr+WPX*I_ERI_F3x_S_Fx2y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Fx2y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Fx2y_S_M1_vrr+oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xy_S_Fx2y_S_vrr = PAY*I_ERI_F3x_S_Fx2y_S_vrr+WPY*I_ERI_F3x_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_G3xz_S_Fx2y_S_vrr = PAZ*I_ERI_F3x_S_Fx2y_S_vrr+WPZ*I_ERI_F3x_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fx2y_S_vrr = PAY*I_ERI_F2xy_S_Fx2y_S_vrr+WPY*I_ERI_F2xy_S_Fx2y_S_M1_vrr+oned2z*I_ERI_D2x_S_Fx2y_S_vrr-rhod2zsq*I_ERI_D2x_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Dxy_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fx2y_S_vrr = PAZ*I_ERI_F2xy_S_Fx2y_S_vrr+WPZ*I_ERI_F2xy_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fx2y_S_vrr = PAZ*I_ERI_F2xz_S_Fx2y_S_vrr+WPZ*I_ERI_F2xz_S_Fx2y_S_M1_vrr+oned2z*I_ERI_D2x_S_Fx2y_S_vrr-rhod2zsq*I_ERI_D2x_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fx2y_S_vrr = PAX*I_ERI_F3y_S_Fx2y_S_vrr+WPX*I_ERI_F3y_S_Fx2y_S_M1_vrr+oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fx2y_S_vrr = PAZ*I_ERI_Fx2y_S_Fx2y_S_vrr+WPZ*I_ERI_Fx2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fx2y_S_vrr = PAY*I_ERI_Fx2z_S_Fx2y_S_vrr+WPY*I_ERI_Fx2z_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fx2y_S_vrr = PAX*I_ERI_F3z_S_Fx2y_S_vrr+WPX*I_ERI_F3z_S_Fx2y_S_M1_vrr+oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4y_S_Fx2y_S_vrr = PAY*I_ERI_F3y_S_Fx2y_S_vrr+WPY*I_ERI_F3y_S_Fx2y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Fx2y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_G3yz_S_Fx2y_S_vrr = PAZ*I_ERI_F3y_S_Fx2y_S_vrr+WPZ*I_ERI_F3y_S_Fx2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fx2y_S_vrr = PAZ*I_ERI_F2yz_S_Fx2y_S_vrr+WPZ*I_ERI_F2yz_S_Fx2y_S_M1_vrr+oned2z*I_ERI_D2y_S_Fx2y_S_vrr-rhod2zsq*I_ERI_D2y_S_Fx2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fx2y_S_vrr = PAY*I_ERI_F3z_S_Fx2y_S_vrr+WPY*I_ERI_F3z_S_Fx2y_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4z_S_Fx2y_S_vrr = PAZ*I_ERI_F3z_S_Fx2y_S_vrr+WPZ*I_ERI_F3z_S_Fx2y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Fx2y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Fx2y_S_M1_vrr;
      Double I_ERI_G4x_S_Fxyz_S_vrr = PAX*I_ERI_F3x_S_Fxyz_S_vrr+WPX*I_ERI_F3x_S_Fxyz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Fxyz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3x_S_Dyz_S_M1_vrr;
      Double I_ERI_G3xy_S_Fxyz_S_vrr = PAY*I_ERI_F3x_S_Fxyz_S_vrr+WPY*I_ERI_F3x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3x_S_Dxz_S_M1_vrr;
      Double I_ERI_G3xz_S_Fxyz_S_vrr = PAZ*I_ERI_F3x_S_Fxyz_S_vrr+WPZ*I_ERI_F3x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3x_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fxyz_S_vrr = PAY*I_ERI_F2xy_S_Fxyz_S_vrr+WPY*I_ERI_F2xy_S_Fxyz_S_M1_vrr+oned2z*I_ERI_D2x_S_Fxyz_S_vrr-rhod2zsq*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_Dxz_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fxyz_S_vrr = PAZ*I_ERI_F2xy_S_Fxyz_S_vrr+WPZ*I_ERI_F2xy_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F2xy_S_Dxy_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fxyz_S_vrr = PAZ*I_ERI_F2xz_S_Fxyz_S_vrr+WPZ*I_ERI_F2xz_S_Fxyz_S_M1_vrr+oned2z*I_ERI_D2x_S_Fxyz_S_vrr-rhod2zsq*I_ERI_D2x_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F2xz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fxyz_S_vrr = PAX*I_ERI_F3y_S_Fxyz_S_vrr+WPX*I_ERI_F3y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3y_S_Dyz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fxyz_S_vrr = PAZ*I_ERI_Fx2y_S_Fxyz_S_vrr+WPZ*I_ERI_Fx2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fxyz_S_vrr = PAY*I_ERI_Fx2z_S_Fxyz_S_vrr+WPY*I_ERI_Fx2z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_Fx2z_S_Dxz_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fxyz_S_vrr = PAX*I_ERI_F3z_S_Fxyz_S_vrr+WPX*I_ERI_F3z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3z_S_Dyz_S_M1_vrr;
      Double I_ERI_G4y_S_Fxyz_S_vrr = PAY*I_ERI_F3y_S_Fxyz_S_vrr+WPY*I_ERI_F3y_S_Fxyz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Fxyz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3y_S_Dxz_S_M1_vrr;
      Double I_ERI_G3yz_S_Fxyz_S_vrr = PAZ*I_ERI_F3y_S_Fxyz_S_vrr+WPZ*I_ERI_F3y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3y_S_Dxy_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fxyz_S_vrr = PAZ*I_ERI_F2yz_S_Fxyz_S_vrr+WPZ*I_ERI_F2yz_S_Fxyz_S_M1_vrr+oned2z*I_ERI_D2y_S_Fxyz_S_vrr-rhod2zsq*I_ERI_D2y_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F2yz_S_Dxy_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fxyz_S_vrr = PAY*I_ERI_F3z_S_Fxyz_S_vrr+WPY*I_ERI_F3z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3z_S_Dxz_S_M1_vrr;
      Double I_ERI_G4z_S_Fxyz_S_vrr = PAZ*I_ERI_F3z_S_Fxyz_S_vrr+WPZ*I_ERI_F3z_S_Fxyz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Fxyz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Fxyz_S_M1_vrr+oned2k*I_ERI_F3z_S_Dxy_S_M1_vrr;
      Double I_ERI_G4x_S_Fx2z_S_vrr = PAX*I_ERI_F3x_S_Fx2z_S_vrr+WPX*I_ERI_F3x_S_Fx2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Fx2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Fx2z_S_M1_vrr+oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Fx2z_S_vrr = PAY*I_ERI_F3x_S_Fx2z_S_vrr+WPY*I_ERI_F3x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Fx2z_S_vrr = PAZ*I_ERI_F3x_S_Fx2z_S_vrr+WPZ*I_ERI_F3x_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dxz_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fx2z_S_vrr = PAY*I_ERI_F2xy_S_Fx2z_S_vrr+WPY*I_ERI_F2xy_S_Fx2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Fx2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Fx2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fx2z_S_vrr = PAZ*I_ERI_F2xy_S_Fx2z_S_vrr+WPZ*I_ERI_F2xy_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Dxz_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fx2z_S_vrr = PAZ*I_ERI_F2xz_S_Fx2z_S_vrr+WPZ*I_ERI_F2xz_S_Fx2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Fx2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Dxz_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fx2z_S_vrr = PAX*I_ERI_F3y_S_Fx2z_S_vrr+WPX*I_ERI_F3y_S_Fx2z_S_M1_vrr+oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fx2z_S_vrr = PAZ*I_ERI_Fx2y_S_Fx2z_S_vrr+WPZ*I_ERI_Fx2y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_Dxz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fx2z_S_vrr = PAY*I_ERI_Fx2z_S_Fx2z_S_vrr+WPY*I_ERI_Fx2z_S_Fx2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fx2z_S_vrr = PAX*I_ERI_F3z_S_Fx2z_S_vrr+WPX*I_ERI_F3z_S_Fx2z_S_M1_vrr+oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4y_S_Fx2z_S_vrr = PAY*I_ERI_F3y_S_Fx2z_S_vrr+WPY*I_ERI_F3y_S_Fx2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Fx2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Fx2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Fx2z_S_vrr = PAZ*I_ERI_F3y_S_Fx2z_S_vrr+WPZ*I_ERI_F3y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dxz_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fx2z_S_vrr = PAZ*I_ERI_F2yz_S_Fx2z_S_vrr+WPZ*I_ERI_F2yz_S_Fx2z_S_M1_vrr+oned2z*I_ERI_D2y_S_Fx2z_S_vrr-rhod2zsq*I_ERI_D2y_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Dxz_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fx2z_S_vrr = PAY*I_ERI_F3z_S_Fx2z_S_vrr+WPY*I_ERI_F3z_S_Fx2z_S_M1_vrr;
      Double I_ERI_G4z_S_Fx2z_S_vrr = PAZ*I_ERI_F3z_S_Fx2z_S_vrr+WPZ*I_ERI_F3z_S_Fx2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Fx2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Fx2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dxz_S_M1_vrr;
      Double I_ERI_G4x_S_F3y_S_vrr = PAX*I_ERI_F3x_S_F3y_S_vrr+WPX*I_ERI_F3x_S_F3y_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F3y_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_G3xy_S_F3y_S_vrr = PAY*I_ERI_F3x_S_F3y_S_vrr+WPY*I_ERI_F3x_S_F3y_S_M1_vrr+3*oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G3xz_S_F3y_S_vrr = PAZ*I_ERI_F3x_S_F3y_S_vrr+WPZ*I_ERI_F3x_S_F3y_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3y_S_vrr = PAY*I_ERI_F2xy_S_F3y_S_vrr+WPY*I_ERI_F2xy_S_F3y_S_M1_vrr+oned2z*I_ERI_D2x_S_F3y_S_vrr-rhod2zsq*I_ERI_D2x_S_F3y_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3y_S_vrr = PAZ*I_ERI_F2xy_S_F3y_S_vrr+WPZ*I_ERI_F2xy_S_F3y_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3y_S_vrr = PAZ*I_ERI_F2xz_S_F3y_S_vrr+WPZ*I_ERI_F2xz_S_F3y_S_M1_vrr+oned2z*I_ERI_D2x_S_F3y_S_vrr-rhod2zsq*I_ERI_D2x_S_F3y_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3y_S_vrr = PAX*I_ERI_F3y_S_F3y_S_vrr+WPX*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3y_S_vrr = PAZ*I_ERI_Fx2y_S_F3y_S_vrr+WPZ*I_ERI_Fx2y_S_F3y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3y_S_vrr = PAY*I_ERI_Fx2z_S_F3y_S_vrr+WPY*I_ERI_Fx2z_S_F3y_S_M1_vrr+3*oned2k*I_ERI_Fx2z_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3y_S_vrr = PAX*I_ERI_F3z_S_F3y_S_vrr+WPX*I_ERI_F3z_S_F3y_S_M1_vrr;
      Double I_ERI_G4y_S_F3y_S_vrr = PAY*I_ERI_F3y_S_F3y_S_vrr+WPY*I_ERI_F3y_S_F3y_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F3y_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F3y_S_M1_vrr+3*oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_G3yz_S_F3y_S_vrr = PAZ*I_ERI_F3y_S_F3y_S_vrr+WPZ*I_ERI_F3y_S_F3y_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3y_S_vrr = PAZ*I_ERI_F2yz_S_F3y_S_vrr+WPZ*I_ERI_F2yz_S_F3y_S_M1_vrr+oned2z*I_ERI_D2y_S_F3y_S_vrr-rhod2zsq*I_ERI_D2y_S_F3y_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3y_S_vrr = PAY*I_ERI_F3z_S_F3y_S_vrr+WPY*I_ERI_F3z_S_F3y_S_M1_vrr+3*oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4z_S_F3y_S_vrr = PAZ*I_ERI_F3z_S_F3y_S_vrr+WPZ*I_ERI_F3z_S_F3y_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F3y_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F3y_S_M1_vrr;
      Double I_ERI_G4x_S_F2yz_S_vrr = PAX*I_ERI_F3x_S_F2yz_S_vrr+WPX*I_ERI_F3x_S_F2yz_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F2yz_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F2yz_S_M1_vrr;
      Double I_ERI_G3xy_S_F2yz_S_vrr = PAY*I_ERI_F3x_S_F2yz_S_vrr+WPY*I_ERI_F3x_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dyz_S_M1_vrr;
      Double I_ERI_G3xz_S_F2yz_S_vrr = PAZ*I_ERI_F3x_S_F2yz_S_vrr+WPZ*I_ERI_F3x_S_F2yz_S_M1_vrr+oned2k*I_ERI_F3x_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2y_S_F2yz_S_vrr = PAY*I_ERI_F2xy_S_F2yz_S_vrr+WPY*I_ERI_F2xy_S_F2yz_S_M1_vrr+oned2z*I_ERI_D2x_S_F2yz_S_vrr-rhod2zsq*I_ERI_D2x_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Dyz_S_M1_vrr;
      Double I_ERI_G2xyz_S_F2yz_S_vrr = PAZ*I_ERI_F2xy_S_F2yz_S_vrr+WPZ*I_ERI_F2xy_S_F2yz_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2y_S_M1_vrr;
      Double I_ERI_G2x2z_S_F2yz_S_vrr = PAZ*I_ERI_F2xz_S_F2yz_S_vrr+WPZ*I_ERI_F2xz_S_F2yz_S_M1_vrr+oned2z*I_ERI_D2x_S_F2yz_S_vrr-rhod2zsq*I_ERI_D2x_S_F2yz_S_M1_vrr+oned2k*I_ERI_F2xz_S_D2y_S_M1_vrr;
      Double I_ERI_Gx3y_S_F2yz_S_vrr = PAX*I_ERI_F3y_S_F2yz_S_vrr+WPX*I_ERI_F3y_S_F2yz_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F2yz_S_vrr = PAZ*I_ERI_Fx2y_S_F2yz_S_vrr+WPZ*I_ERI_Fx2y_S_F2yz_S_M1_vrr+oned2k*I_ERI_Fx2y_S_D2y_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F2yz_S_vrr = PAY*I_ERI_Fx2z_S_F2yz_S_vrr+WPY*I_ERI_Fx2z_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_Fx2z_S_Dyz_S_M1_vrr;
      Double I_ERI_Gx3z_S_F2yz_S_vrr = PAX*I_ERI_F3z_S_F2yz_S_vrr+WPX*I_ERI_F3z_S_F2yz_S_M1_vrr;
      Double I_ERI_G4y_S_F2yz_S_vrr = PAY*I_ERI_F3y_S_F2yz_S_vrr+WPY*I_ERI_F3y_S_F2yz_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F2yz_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dyz_S_M1_vrr;
      Double I_ERI_G3yz_S_F2yz_S_vrr = PAZ*I_ERI_F3y_S_F2yz_S_vrr+WPZ*I_ERI_F3y_S_F2yz_S_M1_vrr+oned2k*I_ERI_F3y_S_D2y_S_M1_vrr;
      Double I_ERI_G2y2z_S_F2yz_S_vrr = PAZ*I_ERI_F2yz_S_F2yz_S_vrr+WPZ*I_ERI_F2yz_S_F2yz_S_M1_vrr+oned2z*I_ERI_D2y_S_F2yz_S_vrr-rhod2zsq*I_ERI_D2y_S_F2yz_S_M1_vrr+oned2k*I_ERI_F2yz_S_D2y_S_M1_vrr;
      Double I_ERI_Gy3z_S_F2yz_S_vrr = PAY*I_ERI_F3z_S_F2yz_S_vrr+WPY*I_ERI_F3z_S_F2yz_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dyz_S_M1_vrr;
      Double I_ERI_G4z_S_F2yz_S_vrr = PAZ*I_ERI_F3z_S_F2yz_S_vrr+WPZ*I_ERI_F3z_S_F2yz_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F2yz_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F2yz_S_M1_vrr+oned2k*I_ERI_F3z_S_D2y_S_M1_vrr;
      Double I_ERI_G4x_S_Fy2z_S_vrr = PAX*I_ERI_F3x_S_Fy2z_S_vrr+WPX*I_ERI_F3x_S_Fy2z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_Fy2z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_Fy2z_S_M1_vrr;
      Double I_ERI_G3xy_S_Fy2z_S_vrr = PAY*I_ERI_F3x_S_Fy2z_S_vrr+WPY*I_ERI_F3x_S_Fy2z_S_M1_vrr+oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G3xz_S_Fy2z_S_vrr = PAZ*I_ERI_F3x_S_Fy2z_S_vrr+WPZ*I_ERI_F3x_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F3x_S_Dyz_S_M1_vrr;
      Double I_ERI_G2x2y_S_Fy2z_S_vrr = PAY*I_ERI_F2xy_S_Fy2z_S_vrr+WPY*I_ERI_F2xy_S_Fy2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Fy2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Fy2z_S_M1_vrr+oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G2xyz_S_Fy2z_S_vrr = PAZ*I_ERI_F2xy_S_Fy2z_S_vrr+WPZ*I_ERI_F2xy_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F2xy_S_Dyz_S_M1_vrr;
      Double I_ERI_G2x2z_S_Fy2z_S_vrr = PAZ*I_ERI_F2xz_S_Fy2z_S_vrr+WPZ*I_ERI_F2xz_S_Fy2z_S_M1_vrr+oned2z*I_ERI_D2x_S_Fy2z_S_vrr-rhod2zsq*I_ERI_D2x_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F2xz_S_Dyz_S_M1_vrr;
      Double I_ERI_Gx3y_S_Fy2z_S_vrr = PAX*I_ERI_F3y_S_Fy2z_S_vrr+WPX*I_ERI_F3y_S_Fy2z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_Fy2z_S_vrr = PAZ*I_ERI_Fx2y_S_Fy2z_S_vrr+WPZ*I_ERI_Fx2y_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_Fx2y_S_Dyz_S_M1_vrr;
      Double I_ERI_Gxy2z_S_Fy2z_S_vrr = PAY*I_ERI_Fx2z_S_Fy2z_S_vrr+WPY*I_ERI_Fx2z_S_Fy2z_S_M1_vrr+oned2k*I_ERI_Fx2z_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3z_S_Fy2z_S_vrr = PAX*I_ERI_F3z_S_Fy2z_S_vrr+WPX*I_ERI_F3z_S_Fy2z_S_M1_vrr;
      Double I_ERI_G4y_S_Fy2z_S_vrr = PAY*I_ERI_F3y_S_Fy2z_S_vrr+WPY*I_ERI_F3y_S_Fy2z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_Fy2z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_Fy2z_S_M1_vrr+oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_G3yz_S_Fy2z_S_vrr = PAZ*I_ERI_F3y_S_Fy2z_S_vrr+WPZ*I_ERI_F3y_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F3y_S_Dyz_S_M1_vrr;
      Double I_ERI_G2y2z_S_Fy2z_S_vrr = PAZ*I_ERI_F2yz_S_Fy2z_S_vrr+WPZ*I_ERI_F2yz_S_Fy2z_S_M1_vrr+oned2z*I_ERI_D2y_S_Fy2z_S_vrr-rhod2zsq*I_ERI_D2y_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F2yz_S_Dyz_S_M1_vrr;
      Double I_ERI_Gy3z_S_Fy2z_S_vrr = PAY*I_ERI_F3z_S_Fy2z_S_vrr+WPY*I_ERI_F3z_S_Fy2z_S_M1_vrr+oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;
      Double I_ERI_G4z_S_Fy2z_S_vrr = PAZ*I_ERI_F3z_S_Fy2z_S_vrr+WPZ*I_ERI_F3z_S_Fy2z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_Fy2z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_Fy2z_S_M1_vrr+2*oned2k*I_ERI_F3z_S_Dyz_S_M1_vrr;
      Double I_ERI_G4x_S_F3z_S_vrr = PAX*I_ERI_F3x_S_F3z_S_vrr+WPX*I_ERI_F3x_S_F3z_S_M1_vrr+3*oned2z*I_ERI_D2x_S_F3z_S_vrr-3*rhod2zsq*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_G3xy_S_F3z_S_vrr = PAY*I_ERI_F3x_S_F3z_S_vrr+WPY*I_ERI_F3x_S_F3z_S_M1_vrr;
      Double I_ERI_G3xz_S_F3z_S_vrr = PAZ*I_ERI_F3x_S_F3z_S_vrr+WPZ*I_ERI_F3x_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F3x_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2y_S_F3z_S_vrr = PAY*I_ERI_F2xy_S_F3z_S_vrr+WPY*I_ERI_F2xy_S_F3z_S_M1_vrr+oned2z*I_ERI_D2x_S_F3z_S_vrr-rhod2zsq*I_ERI_D2x_S_F3z_S_M1_vrr;
      Double I_ERI_G2xyz_S_F3z_S_vrr = PAZ*I_ERI_F2xy_S_F3z_S_vrr+WPZ*I_ERI_F2xy_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F2xy_S_D2z_S_M1_vrr;
      Double I_ERI_G2x2z_S_F3z_S_vrr = PAZ*I_ERI_F2xz_S_F3z_S_vrr+WPZ*I_ERI_F2xz_S_F3z_S_M1_vrr+oned2z*I_ERI_D2x_S_F3z_S_vrr-rhod2zsq*I_ERI_D2x_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F2xz_S_D2z_S_M1_vrr;
      Double I_ERI_Gx3y_S_F3z_S_vrr = PAX*I_ERI_F3y_S_F3z_S_vrr+WPX*I_ERI_F3y_S_F3z_S_M1_vrr;
      Double I_ERI_Gx2yz_S_F3z_S_vrr = PAZ*I_ERI_Fx2y_S_F3z_S_vrr+WPZ*I_ERI_Fx2y_S_F3z_S_M1_vrr+3*oned2k*I_ERI_Fx2y_S_D2z_S_M1_vrr;
      Double I_ERI_Gxy2z_S_F3z_S_vrr = PAY*I_ERI_Fx2z_S_F3z_S_vrr+WPY*I_ERI_Fx2z_S_F3z_S_M1_vrr;
      Double I_ERI_Gx3z_S_F3z_S_vrr = PAX*I_ERI_F3z_S_F3z_S_vrr+WPX*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_G4y_S_F3z_S_vrr = PAY*I_ERI_F3y_S_F3z_S_vrr+WPY*I_ERI_F3y_S_F3z_S_M1_vrr+3*oned2z*I_ERI_D2y_S_F3z_S_vrr-3*rhod2zsq*I_ERI_D2y_S_F3z_S_M1_vrr;
      Double I_ERI_G3yz_S_F3z_S_vrr = PAZ*I_ERI_F3y_S_F3z_S_vrr+WPZ*I_ERI_F3y_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F3y_S_D2z_S_M1_vrr;
      Double I_ERI_G2y2z_S_F3z_S_vrr = PAZ*I_ERI_F2yz_S_F3z_S_vrr+WPZ*I_ERI_F2yz_S_F3z_S_M1_vrr+oned2z*I_ERI_D2y_S_F3z_S_vrr-rhod2zsq*I_ERI_D2y_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F2yz_S_D2z_S_M1_vrr;
      Double I_ERI_Gy3z_S_F3z_S_vrr = PAY*I_ERI_F3z_S_F3z_S_vrr+WPY*I_ERI_F3z_S_F3z_S_M1_vrr;
      Double I_ERI_G4z_S_F3z_S_vrr = PAZ*I_ERI_F3z_S_F3z_S_vrr+WPZ*I_ERI_F3z_S_F3z_S_M1_vrr+3*oned2z*I_ERI_D2z_S_F3z_S_vrr-3*rhod2zsq*I_ERI_D2z_S_F3z_S_M1_vrr+3*oned2k*I_ERI_F3z_S_D2z_S_M1_vrr;

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

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_F_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_G4x_S_F3x_S += I_ERI_G4x_S_F3x_S_vrr;
      I_ERI_G3xy_S_F3x_S += I_ERI_G3xy_S_F3x_S_vrr;
      I_ERI_G3xz_S_F3x_S += I_ERI_G3xz_S_F3x_S_vrr;
      I_ERI_G2x2y_S_F3x_S += I_ERI_G2x2y_S_F3x_S_vrr;
      I_ERI_G2xyz_S_F3x_S += I_ERI_G2xyz_S_F3x_S_vrr;
      I_ERI_G2x2z_S_F3x_S += I_ERI_G2x2z_S_F3x_S_vrr;
      I_ERI_Gx3y_S_F3x_S += I_ERI_Gx3y_S_F3x_S_vrr;
      I_ERI_Gx2yz_S_F3x_S += I_ERI_Gx2yz_S_F3x_S_vrr;
      I_ERI_Gxy2z_S_F3x_S += I_ERI_Gxy2z_S_F3x_S_vrr;
      I_ERI_Gx3z_S_F3x_S += I_ERI_Gx3z_S_F3x_S_vrr;
      I_ERI_G4y_S_F3x_S += I_ERI_G4y_S_F3x_S_vrr;
      I_ERI_G3yz_S_F3x_S += I_ERI_G3yz_S_F3x_S_vrr;
      I_ERI_G2y2z_S_F3x_S += I_ERI_G2y2z_S_F3x_S_vrr;
      I_ERI_Gy3z_S_F3x_S += I_ERI_Gy3z_S_F3x_S_vrr;
      I_ERI_G4z_S_F3x_S += I_ERI_G4z_S_F3x_S_vrr;
      I_ERI_G4x_S_F2xy_S += I_ERI_G4x_S_F2xy_S_vrr;
      I_ERI_G3xy_S_F2xy_S += I_ERI_G3xy_S_F2xy_S_vrr;
      I_ERI_G3xz_S_F2xy_S += I_ERI_G3xz_S_F2xy_S_vrr;
      I_ERI_G2x2y_S_F2xy_S += I_ERI_G2x2y_S_F2xy_S_vrr;
      I_ERI_G2xyz_S_F2xy_S += I_ERI_G2xyz_S_F2xy_S_vrr;
      I_ERI_G2x2z_S_F2xy_S += I_ERI_G2x2z_S_F2xy_S_vrr;
      I_ERI_Gx3y_S_F2xy_S += I_ERI_Gx3y_S_F2xy_S_vrr;
      I_ERI_Gx2yz_S_F2xy_S += I_ERI_Gx2yz_S_F2xy_S_vrr;
      I_ERI_Gxy2z_S_F2xy_S += I_ERI_Gxy2z_S_F2xy_S_vrr;
      I_ERI_Gx3z_S_F2xy_S += I_ERI_Gx3z_S_F2xy_S_vrr;
      I_ERI_G4y_S_F2xy_S += I_ERI_G4y_S_F2xy_S_vrr;
      I_ERI_G3yz_S_F2xy_S += I_ERI_G3yz_S_F2xy_S_vrr;
      I_ERI_G2y2z_S_F2xy_S += I_ERI_G2y2z_S_F2xy_S_vrr;
      I_ERI_Gy3z_S_F2xy_S += I_ERI_Gy3z_S_F2xy_S_vrr;
      I_ERI_G4z_S_F2xy_S += I_ERI_G4z_S_F2xy_S_vrr;
      I_ERI_G4x_S_F2xz_S += I_ERI_G4x_S_F2xz_S_vrr;
      I_ERI_G3xy_S_F2xz_S += I_ERI_G3xy_S_F2xz_S_vrr;
      I_ERI_G3xz_S_F2xz_S += I_ERI_G3xz_S_F2xz_S_vrr;
      I_ERI_G2x2y_S_F2xz_S += I_ERI_G2x2y_S_F2xz_S_vrr;
      I_ERI_G2xyz_S_F2xz_S += I_ERI_G2xyz_S_F2xz_S_vrr;
      I_ERI_G2x2z_S_F2xz_S += I_ERI_G2x2z_S_F2xz_S_vrr;
      I_ERI_Gx3y_S_F2xz_S += I_ERI_Gx3y_S_F2xz_S_vrr;
      I_ERI_Gx2yz_S_F2xz_S += I_ERI_Gx2yz_S_F2xz_S_vrr;
      I_ERI_Gxy2z_S_F2xz_S += I_ERI_Gxy2z_S_F2xz_S_vrr;
      I_ERI_Gx3z_S_F2xz_S += I_ERI_Gx3z_S_F2xz_S_vrr;
      I_ERI_G4y_S_F2xz_S += I_ERI_G4y_S_F2xz_S_vrr;
      I_ERI_G3yz_S_F2xz_S += I_ERI_G3yz_S_F2xz_S_vrr;
      I_ERI_G2y2z_S_F2xz_S += I_ERI_G2y2z_S_F2xz_S_vrr;
      I_ERI_Gy3z_S_F2xz_S += I_ERI_Gy3z_S_F2xz_S_vrr;
      I_ERI_G4z_S_F2xz_S += I_ERI_G4z_S_F2xz_S_vrr;
      I_ERI_G4x_S_Fx2y_S += I_ERI_G4x_S_Fx2y_S_vrr;
      I_ERI_G3xy_S_Fx2y_S += I_ERI_G3xy_S_Fx2y_S_vrr;
      I_ERI_G3xz_S_Fx2y_S += I_ERI_G3xz_S_Fx2y_S_vrr;
      I_ERI_G2x2y_S_Fx2y_S += I_ERI_G2x2y_S_Fx2y_S_vrr;
      I_ERI_G2xyz_S_Fx2y_S += I_ERI_G2xyz_S_Fx2y_S_vrr;
      I_ERI_G2x2z_S_Fx2y_S += I_ERI_G2x2z_S_Fx2y_S_vrr;
      I_ERI_Gx3y_S_Fx2y_S += I_ERI_Gx3y_S_Fx2y_S_vrr;
      I_ERI_Gx2yz_S_Fx2y_S += I_ERI_Gx2yz_S_Fx2y_S_vrr;
      I_ERI_Gxy2z_S_Fx2y_S += I_ERI_Gxy2z_S_Fx2y_S_vrr;
      I_ERI_Gx3z_S_Fx2y_S += I_ERI_Gx3z_S_Fx2y_S_vrr;
      I_ERI_G4y_S_Fx2y_S += I_ERI_G4y_S_Fx2y_S_vrr;
      I_ERI_G3yz_S_Fx2y_S += I_ERI_G3yz_S_Fx2y_S_vrr;
      I_ERI_G2y2z_S_Fx2y_S += I_ERI_G2y2z_S_Fx2y_S_vrr;
      I_ERI_Gy3z_S_Fx2y_S += I_ERI_Gy3z_S_Fx2y_S_vrr;
      I_ERI_G4z_S_Fx2y_S += I_ERI_G4z_S_Fx2y_S_vrr;
      I_ERI_G4x_S_Fxyz_S += I_ERI_G4x_S_Fxyz_S_vrr;
      I_ERI_G3xy_S_Fxyz_S += I_ERI_G3xy_S_Fxyz_S_vrr;
      I_ERI_G3xz_S_Fxyz_S += I_ERI_G3xz_S_Fxyz_S_vrr;
      I_ERI_G2x2y_S_Fxyz_S += I_ERI_G2x2y_S_Fxyz_S_vrr;
      I_ERI_G2xyz_S_Fxyz_S += I_ERI_G2xyz_S_Fxyz_S_vrr;
      I_ERI_G2x2z_S_Fxyz_S += I_ERI_G2x2z_S_Fxyz_S_vrr;
      I_ERI_Gx3y_S_Fxyz_S += I_ERI_Gx3y_S_Fxyz_S_vrr;
      I_ERI_Gx2yz_S_Fxyz_S += I_ERI_Gx2yz_S_Fxyz_S_vrr;
      I_ERI_Gxy2z_S_Fxyz_S += I_ERI_Gxy2z_S_Fxyz_S_vrr;
      I_ERI_Gx3z_S_Fxyz_S += I_ERI_Gx3z_S_Fxyz_S_vrr;
      I_ERI_G4y_S_Fxyz_S += I_ERI_G4y_S_Fxyz_S_vrr;
      I_ERI_G3yz_S_Fxyz_S += I_ERI_G3yz_S_Fxyz_S_vrr;
      I_ERI_G2y2z_S_Fxyz_S += I_ERI_G2y2z_S_Fxyz_S_vrr;
      I_ERI_Gy3z_S_Fxyz_S += I_ERI_Gy3z_S_Fxyz_S_vrr;
      I_ERI_G4z_S_Fxyz_S += I_ERI_G4z_S_Fxyz_S_vrr;
      I_ERI_G4x_S_Fx2z_S += I_ERI_G4x_S_Fx2z_S_vrr;
      I_ERI_G3xy_S_Fx2z_S += I_ERI_G3xy_S_Fx2z_S_vrr;
      I_ERI_G3xz_S_Fx2z_S += I_ERI_G3xz_S_Fx2z_S_vrr;
      I_ERI_G2x2y_S_Fx2z_S += I_ERI_G2x2y_S_Fx2z_S_vrr;
      I_ERI_G2xyz_S_Fx2z_S += I_ERI_G2xyz_S_Fx2z_S_vrr;
      I_ERI_G2x2z_S_Fx2z_S += I_ERI_G2x2z_S_Fx2z_S_vrr;
      I_ERI_Gx3y_S_Fx2z_S += I_ERI_Gx3y_S_Fx2z_S_vrr;
      I_ERI_Gx2yz_S_Fx2z_S += I_ERI_Gx2yz_S_Fx2z_S_vrr;
      I_ERI_Gxy2z_S_Fx2z_S += I_ERI_Gxy2z_S_Fx2z_S_vrr;
      I_ERI_Gx3z_S_Fx2z_S += I_ERI_Gx3z_S_Fx2z_S_vrr;
      I_ERI_G4y_S_Fx2z_S += I_ERI_G4y_S_Fx2z_S_vrr;
      I_ERI_G3yz_S_Fx2z_S += I_ERI_G3yz_S_Fx2z_S_vrr;
      I_ERI_G2y2z_S_Fx2z_S += I_ERI_G2y2z_S_Fx2z_S_vrr;
      I_ERI_Gy3z_S_Fx2z_S += I_ERI_Gy3z_S_Fx2z_S_vrr;
      I_ERI_G4z_S_Fx2z_S += I_ERI_G4z_S_Fx2z_S_vrr;
      I_ERI_G4x_S_F3y_S += I_ERI_G4x_S_F3y_S_vrr;
      I_ERI_G3xy_S_F3y_S += I_ERI_G3xy_S_F3y_S_vrr;
      I_ERI_G3xz_S_F3y_S += I_ERI_G3xz_S_F3y_S_vrr;
      I_ERI_G2x2y_S_F3y_S += I_ERI_G2x2y_S_F3y_S_vrr;
      I_ERI_G2xyz_S_F3y_S += I_ERI_G2xyz_S_F3y_S_vrr;
      I_ERI_G2x2z_S_F3y_S += I_ERI_G2x2z_S_F3y_S_vrr;
      I_ERI_Gx3y_S_F3y_S += I_ERI_Gx3y_S_F3y_S_vrr;
      I_ERI_Gx2yz_S_F3y_S += I_ERI_Gx2yz_S_F3y_S_vrr;
      I_ERI_Gxy2z_S_F3y_S += I_ERI_Gxy2z_S_F3y_S_vrr;
      I_ERI_Gx3z_S_F3y_S += I_ERI_Gx3z_S_F3y_S_vrr;
      I_ERI_G4y_S_F3y_S += I_ERI_G4y_S_F3y_S_vrr;
      I_ERI_G3yz_S_F3y_S += I_ERI_G3yz_S_F3y_S_vrr;
      I_ERI_G2y2z_S_F3y_S += I_ERI_G2y2z_S_F3y_S_vrr;
      I_ERI_Gy3z_S_F3y_S += I_ERI_Gy3z_S_F3y_S_vrr;
      I_ERI_G4z_S_F3y_S += I_ERI_G4z_S_F3y_S_vrr;
      I_ERI_G4x_S_F2yz_S += I_ERI_G4x_S_F2yz_S_vrr;
      I_ERI_G3xy_S_F2yz_S += I_ERI_G3xy_S_F2yz_S_vrr;
      I_ERI_G3xz_S_F2yz_S += I_ERI_G3xz_S_F2yz_S_vrr;
      I_ERI_G2x2y_S_F2yz_S += I_ERI_G2x2y_S_F2yz_S_vrr;
      I_ERI_G2xyz_S_F2yz_S += I_ERI_G2xyz_S_F2yz_S_vrr;
      I_ERI_G2x2z_S_F2yz_S += I_ERI_G2x2z_S_F2yz_S_vrr;
      I_ERI_Gx3y_S_F2yz_S += I_ERI_Gx3y_S_F2yz_S_vrr;
      I_ERI_Gx2yz_S_F2yz_S += I_ERI_Gx2yz_S_F2yz_S_vrr;
      I_ERI_Gxy2z_S_F2yz_S += I_ERI_Gxy2z_S_F2yz_S_vrr;
      I_ERI_Gx3z_S_F2yz_S += I_ERI_Gx3z_S_F2yz_S_vrr;
      I_ERI_G4y_S_F2yz_S += I_ERI_G4y_S_F2yz_S_vrr;
      I_ERI_G3yz_S_F2yz_S += I_ERI_G3yz_S_F2yz_S_vrr;
      I_ERI_G2y2z_S_F2yz_S += I_ERI_G2y2z_S_F2yz_S_vrr;
      I_ERI_Gy3z_S_F2yz_S += I_ERI_Gy3z_S_F2yz_S_vrr;
      I_ERI_G4z_S_F2yz_S += I_ERI_G4z_S_F2yz_S_vrr;
      I_ERI_G4x_S_Fy2z_S += I_ERI_G4x_S_Fy2z_S_vrr;
      I_ERI_G3xy_S_Fy2z_S += I_ERI_G3xy_S_Fy2z_S_vrr;
      I_ERI_G3xz_S_Fy2z_S += I_ERI_G3xz_S_Fy2z_S_vrr;
      I_ERI_G2x2y_S_Fy2z_S += I_ERI_G2x2y_S_Fy2z_S_vrr;
      I_ERI_G2xyz_S_Fy2z_S += I_ERI_G2xyz_S_Fy2z_S_vrr;
      I_ERI_G2x2z_S_Fy2z_S += I_ERI_G2x2z_S_Fy2z_S_vrr;
      I_ERI_Gx3y_S_Fy2z_S += I_ERI_Gx3y_S_Fy2z_S_vrr;
      I_ERI_Gx2yz_S_Fy2z_S += I_ERI_Gx2yz_S_Fy2z_S_vrr;
      I_ERI_Gxy2z_S_Fy2z_S += I_ERI_Gxy2z_S_Fy2z_S_vrr;
      I_ERI_Gx3z_S_Fy2z_S += I_ERI_Gx3z_S_Fy2z_S_vrr;
      I_ERI_G4y_S_Fy2z_S += I_ERI_G4y_S_Fy2z_S_vrr;
      I_ERI_G3yz_S_Fy2z_S += I_ERI_G3yz_S_Fy2z_S_vrr;
      I_ERI_G2y2z_S_Fy2z_S += I_ERI_G2y2z_S_Fy2z_S_vrr;
      I_ERI_Gy3z_S_Fy2z_S += I_ERI_Gy3z_S_Fy2z_S_vrr;
      I_ERI_G4z_S_Fy2z_S += I_ERI_G4z_S_Fy2z_S_vrr;
      I_ERI_G4x_S_F3z_S += I_ERI_G4x_S_F3z_S_vrr;
      I_ERI_G3xy_S_F3z_S += I_ERI_G3xy_S_F3z_S_vrr;
      I_ERI_G3xz_S_F3z_S += I_ERI_G3xz_S_F3z_S_vrr;
      I_ERI_G2x2y_S_F3z_S += I_ERI_G2x2y_S_F3z_S_vrr;
      I_ERI_G2xyz_S_F3z_S += I_ERI_G2xyz_S_F3z_S_vrr;
      I_ERI_G2x2z_S_F3z_S += I_ERI_G2x2z_S_F3z_S_vrr;
      I_ERI_Gx3y_S_F3z_S += I_ERI_Gx3y_S_F3z_S_vrr;
      I_ERI_Gx2yz_S_F3z_S += I_ERI_Gx2yz_S_F3z_S_vrr;
      I_ERI_Gxy2z_S_F3z_S += I_ERI_Gxy2z_S_F3z_S_vrr;
      I_ERI_Gx3z_S_F3z_S += I_ERI_Gx3z_S_F3z_S_vrr;
      I_ERI_G4y_S_F3z_S += I_ERI_G4y_S_F3z_S_vrr;
      I_ERI_G3yz_S_F3z_S += I_ERI_G3yz_S_F3z_S_vrr;
      I_ERI_G2y2z_S_F3z_S += I_ERI_G2y2z_S_F3z_S_vrr;
      I_ERI_Gy3z_S_F3z_S += I_ERI_Gy3z_S_F3z_S_vrr;
      I_ERI_G4z_S_F3z_S += I_ERI_G4z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_G_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_F3x_S_G4x_S += I_ERI_F3x_S_G4x_S_vrr;
      I_ERI_F2xy_S_G4x_S += I_ERI_F2xy_S_G4x_S_vrr;
      I_ERI_F2xz_S_G4x_S += I_ERI_F2xz_S_G4x_S_vrr;
      I_ERI_Fx2y_S_G4x_S += I_ERI_Fx2y_S_G4x_S_vrr;
      I_ERI_Fxyz_S_G4x_S += I_ERI_Fxyz_S_G4x_S_vrr;
      I_ERI_Fx2z_S_G4x_S += I_ERI_Fx2z_S_G4x_S_vrr;
      I_ERI_F3y_S_G4x_S += I_ERI_F3y_S_G4x_S_vrr;
      I_ERI_F2yz_S_G4x_S += I_ERI_F2yz_S_G4x_S_vrr;
      I_ERI_Fy2z_S_G4x_S += I_ERI_Fy2z_S_G4x_S_vrr;
      I_ERI_F3z_S_G4x_S += I_ERI_F3z_S_G4x_S_vrr;
      I_ERI_F3x_S_G3xy_S += I_ERI_F3x_S_G3xy_S_vrr;
      I_ERI_F2xy_S_G3xy_S += I_ERI_F2xy_S_G3xy_S_vrr;
      I_ERI_F2xz_S_G3xy_S += I_ERI_F2xz_S_G3xy_S_vrr;
      I_ERI_Fx2y_S_G3xy_S += I_ERI_Fx2y_S_G3xy_S_vrr;
      I_ERI_Fxyz_S_G3xy_S += I_ERI_Fxyz_S_G3xy_S_vrr;
      I_ERI_Fx2z_S_G3xy_S += I_ERI_Fx2z_S_G3xy_S_vrr;
      I_ERI_F3y_S_G3xy_S += I_ERI_F3y_S_G3xy_S_vrr;
      I_ERI_F2yz_S_G3xy_S += I_ERI_F2yz_S_G3xy_S_vrr;
      I_ERI_Fy2z_S_G3xy_S += I_ERI_Fy2z_S_G3xy_S_vrr;
      I_ERI_F3z_S_G3xy_S += I_ERI_F3z_S_G3xy_S_vrr;
      I_ERI_F3x_S_G3xz_S += I_ERI_F3x_S_G3xz_S_vrr;
      I_ERI_F2xy_S_G3xz_S += I_ERI_F2xy_S_G3xz_S_vrr;
      I_ERI_F2xz_S_G3xz_S += I_ERI_F2xz_S_G3xz_S_vrr;
      I_ERI_Fx2y_S_G3xz_S += I_ERI_Fx2y_S_G3xz_S_vrr;
      I_ERI_Fxyz_S_G3xz_S += I_ERI_Fxyz_S_G3xz_S_vrr;
      I_ERI_Fx2z_S_G3xz_S += I_ERI_Fx2z_S_G3xz_S_vrr;
      I_ERI_F3y_S_G3xz_S += I_ERI_F3y_S_G3xz_S_vrr;
      I_ERI_F2yz_S_G3xz_S += I_ERI_F2yz_S_G3xz_S_vrr;
      I_ERI_Fy2z_S_G3xz_S += I_ERI_Fy2z_S_G3xz_S_vrr;
      I_ERI_F3z_S_G3xz_S += I_ERI_F3z_S_G3xz_S_vrr;
      I_ERI_F3x_S_G2x2y_S += I_ERI_F3x_S_G2x2y_S_vrr;
      I_ERI_F2xy_S_G2x2y_S += I_ERI_F2xy_S_G2x2y_S_vrr;
      I_ERI_F2xz_S_G2x2y_S += I_ERI_F2xz_S_G2x2y_S_vrr;
      I_ERI_Fx2y_S_G2x2y_S += I_ERI_Fx2y_S_G2x2y_S_vrr;
      I_ERI_Fxyz_S_G2x2y_S += I_ERI_Fxyz_S_G2x2y_S_vrr;
      I_ERI_Fx2z_S_G2x2y_S += I_ERI_Fx2z_S_G2x2y_S_vrr;
      I_ERI_F3y_S_G2x2y_S += I_ERI_F3y_S_G2x2y_S_vrr;
      I_ERI_F2yz_S_G2x2y_S += I_ERI_F2yz_S_G2x2y_S_vrr;
      I_ERI_Fy2z_S_G2x2y_S += I_ERI_Fy2z_S_G2x2y_S_vrr;
      I_ERI_F3z_S_G2x2y_S += I_ERI_F3z_S_G2x2y_S_vrr;
      I_ERI_F3x_S_G2xyz_S += I_ERI_F3x_S_G2xyz_S_vrr;
      I_ERI_F2xy_S_G2xyz_S += I_ERI_F2xy_S_G2xyz_S_vrr;
      I_ERI_F2xz_S_G2xyz_S += I_ERI_F2xz_S_G2xyz_S_vrr;
      I_ERI_Fx2y_S_G2xyz_S += I_ERI_Fx2y_S_G2xyz_S_vrr;
      I_ERI_Fxyz_S_G2xyz_S += I_ERI_Fxyz_S_G2xyz_S_vrr;
      I_ERI_Fx2z_S_G2xyz_S += I_ERI_Fx2z_S_G2xyz_S_vrr;
      I_ERI_F3y_S_G2xyz_S += I_ERI_F3y_S_G2xyz_S_vrr;
      I_ERI_F2yz_S_G2xyz_S += I_ERI_F2yz_S_G2xyz_S_vrr;
      I_ERI_Fy2z_S_G2xyz_S += I_ERI_Fy2z_S_G2xyz_S_vrr;
      I_ERI_F3z_S_G2xyz_S += I_ERI_F3z_S_G2xyz_S_vrr;
      I_ERI_F3x_S_G2x2z_S += I_ERI_F3x_S_G2x2z_S_vrr;
      I_ERI_F2xy_S_G2x2z_S += I_ERI_F2xy_S_G2x2z_S_vrr;
      I_ERI_F2xz_S_G2x2z_S += I_ERI_F2xz_S_G2x2z_S_vrr;
      I_ERI_Fx2y_S_G2x2z_S += I_ERI_Fx2y_S_G2x2z_S_vrr;
      I_ERI_Fxyz_S_G2x2z_S += I_ERI_Fxyz_S_G2x2z_S_vrr;
      I_ERI_Fx2z_S_G2x2z_S += I_ERI_Fx2z_S_G2x2z_S_vrr;
      I_ERI_F3y_S_G2x2z_S += I_ERI_F3y_S_G2x2z_S_vrr;
      I_ERI_F2yz_S_G2x2z_S += I_ERI_F2yz_S_G2x2z_S_vrr;
      I_ERI_Fy2z_S_G2x2z_S += I_ERI_Fy2z_S_G2x2z_S_vrr;
      I_ERI_F3z_S_G2x2z_S += I_ERI_F3z_S_G2x2z_S_vrr;
      I_ERI_F3x_S_Gx3y_S += I_ERI_F3x_S_Gx3y_S_vrr;
      I_ERI_F2xy_S_Gx3y_S += I_ERI_F2xy_S_Gx3y_S_vrr;
      I_ERI_F2xz_S_Gx3y_S += I_ERI_F2xz_S_Gx3y_S_vrr;
      I_ERI_Fx2y_S_Gx3y_S += I_ERI_Fx2y_S_Gx3y_S_vrr;
      I_ERI_Fxyz_S_Gx3y_S += I_ERI_Fxyz_S_Gx3y_S_vrr;
      I_ERI_Fx2z_S_Gx3y_S += I_ERI_Fx2z_S_Gx3y_S_vrr;
      I_ERI_F3y_S_Gx3y_S += I_ERI_F3y_S_Gx3y_S_vrr;
      I_ERI_F2yz_S_Gx3y_S += I_ERI_F2yz_S_Gx3y_S_vrr;
      I_ERI_Fy2z_S_Gx3y_S += I_ERI_Fy2z_S_Gx3y_S_vrr;
      I_ERI_F3z_S_Gx3y_S += I_ERI_F3z_S_Gx3y_S_vrr;
      I_ERI_F3x_S_Gx2yz_S += I_ERI_F3x_S_Gx2yz_S_vrr;
      I_ERI_F2xy_S_Gx2yz_S += I_ERI_F2xy_S_Gx2yz_S_vrr;
      I_ERI_F2xz_S_Gx2yz_S += I_ERI_F2xz_S_Gx2yz_S_vrr;
      I_ERI_Fx2y_S_Gx2yz_S += I_ERI_Fx2y_S_Gx2yz_S_vrr;
      I_ERI_Fxyz_S_Gx2yz_S += I_ERI_Fxyz_S_Gx2yz_S_vrr;
      I_ERI_Fx2z_S_Gx2yz_S += I_ERI_Fx2z_S_Gx2yz_S_vrr;
      I_ERI_F3y_S_Gx2yz_S += I_ERI_F3y_S_Gx2yz_S_vrr;
      I_ERI_F2yz_S_Gx2yz_S += I_ERI_F2yz_S_Gx2yz_S_vrr;
      I_ERI_Fy2z_S_Gx2yz_S += I_ERI_Fy2z_S_Gx2yz_S_vrr;
      I_ERI_F3z_S_Gx2yz_S += I_ERI_F3z_S_Gx2yz_S_vrr;
      I_ERI_F3x_S_Gxy2z_S += I_ERI_F3x_S_Gxy2z_S_vrr;
      I_ERI_F2xy_S_Gxy2z_S += I_ERI_F2xy_S_Gxy2z_S_vrr;
      I_ERI_F2xz_S_Gxy2z_S += I_ERI_F2xz_S_Gxy2z_S_vrr;
      I_ERI_Fx2y_S_Gxy2z_S += I_ERI_Fx2y_S_Gxy2z_S_vrr;
      I_ERI_Fxyz_S_Gxy2z_S += I_ERI_Fxyz_S_Gxy2z_S_vrr;
      I_ERI_Fx2z_S_Gxy2z_S += I_ERI_Fx2z_S_Gxy2z_S_vrr;
      I_ERI_F3y_S_Gxy2z_S += I_ERI_F3y_S_Gxy2z_S_vrr;
      I_ERI_F2yz_S_Gxy2z_S += I_ERI_F2yz_S_Gxy2z_S_vrr;
      I_ERI_Fy2z_S_Gxy2z_S += I_ERI_Fy2z_S_Gxy2z_S_vrr;
      I_ERI_F3z_S_Gxy2z_S += I_ERI_F3z_S_Gxy2z_S_vrr;
      I_ERI_F3x_S_Gx3z_S += I_ERI_F3x_S_Gx3z_S_vrr;
      I_ERI_F2xy_S_Gx3z_S += I_ERI_F2xy_S_Gx3z_S_vrr;
      I_ERI_F2xz_S_Gx3z_S += I_ERI_F2xz_S_Gx3z_S_vrr;
      I_ERI_Fx2y_S_Gx3z_S += I_ERI_Fx2y_S_Gx3z_S_vrr;
      I_ERI_Fxyz_S_Gx3z_S += I_ERI_Fxyz_S_Gx3z_S_vrr;
      I_ERI_Fx2z_S_Gx3z_S += I_ERI_Fx2z_S_Gx3z_S_vrr;
      I_ERI_F3y_S_Gx3z_S += I_ERI_F3y_S_Gx3z_S_vrr;
      I_ERI_F2yz_S_Gx3z_S += I_ERI_F2yz_S_Gx3z_S_vrr;
      I_ERI_Fy2z_S_Gx3z_S += I_ERI_Fy2z_S_Gx3z_S_vrr;
      I_ERI_F3z_S_Gx3z_S += I_ERI_F3z_S_Gx3z_S_vrr;
      I_ERI_F3x_S_G4y_S += I_ERI_F3x_S_G4y_S_vrr;
      I_ERI_F2xy_S_G4y_S += I_ERI_F2xy_S_G4y_S_vrr;
      I_ERI_F2xz_S_G4y_S += I_ERI_F2xz_S_G4y_S_vrr;
      I_ERI_Fx2y_S_G4y_S += I_ERI_Fx2y_S_G4y_S_vrr;
      I_ERI_Fxyz_S_G4y_S += I_ERI_Fxyz_S_G4y_S_vrr;
      I_ERI_Fx2z_S_G4y_S += I_ERI_Fx2z_S_G4y_S_vrr;
      I_ERI_F3y_S_G4y_S += I_ERI_F3y_S_G4y_S_vrr;
      I_ERI_F2yz_S_G4y_S += I_ERI_F2yz_S_G4y_S_vrr;
      I_ERI_Fy2z_S_G4y_S += I_ERI_Fy2z_S_G4y_S_vrr;
      I_ERI_F3z_S_G4y_S += I_ERI_F3z_S_G4y_S_vrr;
      I_ERI_F3x_S_G3yz_S += I_ERI_F3x_S_G3yz_S_vrr;
      I_ERI_F2xy_S_G3yz_S += I_ERI_F2xy_S_G3yz_S_vrr;
      I_ERI_F2xz_S_G3yz_S += I_ERI_F2xz_S_G3yz_S_vrr;
      I_ERI_Fx2y_S_G3yz_S += I_ERI_Fx2y_S_G3yz_S_vrr;
      I_ERI_Fxyz_S_G3yz_S += I_ERI_Fxyz_S_G3yz_S_vrr;
      I_ERI_Fx2z_S_G3yz_S += I_ERI_Fx2z_S_G3yz_S_vrr;
      I_ERI_F3y_S_G3yz_S += I_ERI_F3y_S_G3yz_S_vrr;
      I_ERI_F2yz_S_G3yz_S += I_ERI_F2yz_S_G3yz_S_vrr;
      I_ERI_Fy2z_S_G3yz_S += I_ERI_Fy2z_S_G3yz_S_vrr;
      I_ERI_F3z_S_G3yz_S += I_ERI_F3z_S_G3yz_S_vrr;
      I_ERI_F3x_S_G2y2z_S += I_ERI_F3x_S_G2y2z_S_vrr;
      I_ERI_F2xy_S_G2y2z_S += I_ERI_F2xy_S_G2y2z_S_vrr;
      I_ERI_F2xz_S_G2y2z_S += I_ERI_F2xz_S_G2y2z_S_vrr;
      I_ERI_Fx2y_S_G2y2z_S += I_ERI_Fx2y_S_G2y2z_S_vrr;
      I_ERI_Fxyz_S_G2y2z_S += I_ERI_Fxyz_S_G2y2z_S_vrr;
      I_ERI_Fx2z_S_G2y2z_S += I_ERI_Fx2z_S_G2y2z_S_vrr;
      I_ERI_F3y_S_G2y2z_S += I_ERI_F3y_S_G2y2z_S_vrr;
      I_ERI_F2yz_S_G2y2z_S += I_ERI_F2yz_S_G2y2z_S_vrr;
      I_ERI_Fy2z_S_G2y2z_S += I_ERI_Fy2z_S_G2y2z_S_vrr;
      I_ERI_F3z_S_G2y2z_S += I_ERI_F3z_S_G2y2z_S_vrr;
      I_ERI_F3x_S_Gy3z_S += I_ERI_F3x_S_Gy3z_S_vrr;
      I_ERI_F2xy_S_Gy3z_S += I_ERI_F2xy_S_Gy3z_S_vrr;
      I_ERI_F2xz_S_Gy3z_S += I_ERI_F2xz_S_Gy3z_S_vrr;
      I_ERI_Fx2y_S_Gy3z_S += I_ERI_Fx2y_S_Gy3z_S_vrr;
      I_ERI_Fxyz_S_Gy3z_S += I_ERI_Fxyz_S_Gy3z_S_vrr;
      I_ERI_Fx2z_S_Gy3z_S += I_ERI_Fx2z_S_Gy3z_S_vrr;
      I_ERI_F3y_S_Gy3z_S += I_ERI_F3y_S_Gy3z_S_vrr;
      I_ERI_F2yz_S_Gy3z_S += I_ERI_F2yz_S_Gy3z_S_vrr;
      I_ERI_Fy2z_S_Gy3z_S += I_ERI_Fy2z_S_Gy3z_S_vrr;
      I_ERI_F3z_S_Gy3z_S += I_ERI_F3z_S_Gy3z_S_vrr;
      I_ERI_F3x_S_G4z_S += I_ERI_F3x_S_G4z_S_vrr;
      I_ERI_F2xy_S_G4z_S += I_ERI_F2xy_S_G4z_S_vrr;
      I_ERI_F2xz_S_G4z_S += I_ERI_F2xz_S_G4z_S_vrr;
      I_ERI_Fx2y_S_G4z_S += I_ERI_Fx2y_S_G4z_S_vrr;
      I_ERI_Fxyz_S_G4z_S += I_ERI_Fxyz_S_G4z_S_vrr;
      I_ERI_Fx2z_S_G4z_S += I_ERI_Fx2z_S_G4z_S_vrr;
      I_ERI_F3y_S_G4z_S += I_ERI_F3y_S_G4z_S_vrr;
      I_ERI_F2yz_S_G4z_S += I_ERI_F2yz_S_G4z_S_vrr;
      I_ERI_Fy2z_S_G4z_S += I_ERI_Fy2z_S_G4z_S_vrr;
      I_ERI_F3z_S_G4z_S += I_ERI_F3z_S_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_F3x_S_F3x_S += I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S += I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S += I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S += I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S += I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S += I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S += I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S += I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S += I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S += I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S += I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S += I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S += I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S += I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S += I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S += I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S += I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S += I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S += I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S += I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S += I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S += I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S += I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S += I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S += I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S += I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S += I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S += I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S += I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S += I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S += I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S += I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S += I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S += I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S += I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S += I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S += I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S += I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S += I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S += I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S += I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S += I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S += I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S += I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S += I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S += I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S += I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S += I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S += I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S += I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S += I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S += I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S += I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S += I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S += I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S += I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S += I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S += I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S += I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S += I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S += I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S += I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S += I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S += I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S += I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S += I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S += I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S += I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S += I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S += I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S += I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S += I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S += I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S += I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S += I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S += I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S += I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S += I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S += I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S += I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S += I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S += I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S += I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S += I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S += I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S += I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S += I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S += I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S += I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S += I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S += I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S += I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S += I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S += I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S += I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S += I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S += I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S += I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S += I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S += I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_G_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_D2x_S_G4x_S += I_ERI_D2x_S_G4x_S_vrr;
      I_ERI_Dxy_S_G4x_S += I_ERI_Dxy_S_G4x_S_vrr;
      I_ERI_Dxz_S_G4x_S += I_ERI_Dxz_S_G4x_S_vrr;
      I_ERI_D2y_S_G4x_S += I_ERI_D2y_S_G4x_S_vrr;
      I_ERI_Dyz_S_G4x_S += I_ERI_Dyz_S_G4x_S_vrr;
      I_ERI_D2z_S_G4x_S += I_ERI_D2z_S_G4x_S_vrr;
      I_ERI_D2x_S_G3xy_S += I_ERI_D2x_S_G3xy_S_vrr;
      I_ERI_Dxy_S_G3xy_S += I_ERI_Dxy_S_G3xy_S_vrr;
      I_ERI_Dxz_S_G3xy_S += I_ERI_Dxz_S_G3xy_S_vrr;
      I_ERI_D2y_S_G3xy_S += I_ERI_D2y_S_G3xy_S_vrr;
      I_ERI_Dyz_S_G3xy_S += I_ERI_Dyz_S_G3xy_S_vrr;
      I_ERI_D2z_S_G3xy_S += I_ERI_D2z_S_G3xy_S_vrr;
      I_ERI_D2x_S_G3xz_S += I_ERI_D2x_S_G3xz_S_vrr;
      I_ERI_Dxy_S_G3xz_S += I_ERI_Dxy_S_G3xz_S_vrr;
      I_ERI_Dxz_S_G3xz_S += I_ERI_Dxz_S_G3xz_S_vrr;
      I_ERI_D2y_S_G3xz_S += I_ERI_D2y_S_G3xz_S_vrr;
      I_ERI_Dyz_S_G3xz_S += I_ERI_Dyz_S_G3xz_S_vrr;
      I_ERI_D2z_S_G3xz_S += I_ERI_D2z_S_G3xz_S_vrr;
      I_ERI_D2x_S_G2x2y_S += I_ERI_D2x_S_G2x2y_S_vrr;
      I_ERI_Dxy_S_G2x2y_S += I_ERI_Dxy_S_G2x2y_S_vrr;
      I_ERI_Dxz_S_G2x2y_S += I_ERI_Dxz_S_G2x2y_S_vrr;
      I_ERI_D2y_S_G2x2y_S += I_ERI_D2y_S_G2x2y_S_vrr;
      I_ERI_Dyz_S_G2x2y_S += I_ERI_Dyz_S_G2x2y_S_vrr;
      I_ERI_D2z_S_G2x2y_S += I_ERI_D2z_S_G2x2y_S_vrr;
      I_ERI_D2x_S_G2xyz_S += I_ERI_D2x_S_G2xyz_S_vrr;
      I_ERI_Dxy_S_G2xyz_S += I_ERI_Dxy_S_G2xyz_S_vrr;
      I_ERI_Dxz_S_G2xyz_S += I_ERI_Dxz_S_G2xyz_S_vrr;
      I_ERI_D2y_S_G2xyz_S += I_ERI_D2y_S_G2xyz_S_vrr;
      I_ERI_Dyz_S_G2xyz_S += I_ERI_Dyz_S_G2xyz_S_vrr;
      I_ERI_D2z_S_G2xyz_S += I_ERI_D2z_S_G2xyz_S_vrr;
      I_ERI_D2x_S_G2x2z_S += I_ERI_D2x_S_G2x2z_S_vrr;
      I_ERI_Dxy_S_G2x2z_S += I_ERI_Dxy_S_G2x2z_S_vrr;
      I_ERI_Dxz_S_G2x2z_S += I_ERI_Dxz_S_G2x2z_S_vrr;
      I_ERI_D2y_S_G2x2z_S += I_ERI_D2y_S_G2x2z_S_vrr;
      I_ERI_Dyz_S_G2x2z_S += I_ERI_Dyz_S_G2x2z_S_vrr;
      I_ERI_D2z_S_G2x2z_S += I_ERI_D2z_S_G2x2z_S_vrr;
      I_ERI_D2x_S_Gx3y_S += I_ERI_D2x_S_Gx3y_S_vrr;
      I_ERI_Dxy_S_Gx3y_S += I_ERI_Dxy_S_Gx3y_S_vrr;
      I_ERI_Dxz_S_Gx3y_S += I_ERI_Dxz_S_Gx3y_S_vrr;
      I_ERI_D2y_S_Gx3y_S += I_ERI_D2y_S_Gx3y_S_vrr;
      I_ERI_Dyz_S_Gx3y_S += I_ERI_Dyz_S_Gx3y_S_vrr;
      I_ERI_D2z_S_Gx3y_S += I_ERI_D2z_S_Gx3y_S_vrr;
      I_ERI_D2x_S_Gx2yz_S += I_ERI_D2x_S_Gx2yz_S_vrr;
      I_ERI_Dxy_S_Gx2yz_S += I_ERI_Dxy_S_Gx2yz_S_vrr;
      I_ERI_Dxz_S_Gx2yz_S += I_ERI_Dxz_S_Gx2yz_S_vrr;
      I_ERI_D2y_S_Gx2yz_S += I_ERI_D2y_S_Gx2yz_S_vrr;
      I_ERI_Dyz_S_Gx2yz_S += I_ERI_Dyz_S_Gx2yz_S_vrr;
      I_ERI_D2z_S_Gx2yz_S += I_ERI_D2z_S_Gx2yz_S_vrr;
      I_ERI_D2x_S_Gxy2z_S += I_ERI_D2x_S_Gxy2z_S_vrr;
      I_ERI_Dxy_S_Gxy2z_S += I_ERI_Dxy_S_Gxy2z_S_vrr;
      I_ERI_Dxz_S_Gxy2z_S += I_ERI_Dxz_S_Gxy2z_S_vrr;
      I_ERI_D2y_S_Gxy2z_S += I_ERI_D2y_S_Gxy2z_S_vrr;
      I_ERI_Dyz_S_Gxy2z_S += I_ERI_Dyz_S_Gxy2z_S_vrr;
      I_ERI_D2z_S_Gxy2z_S += I_ERI_D2z_S_Gxy2z_S_vrr;
      I_ERI_D2x_S_Gx3z_S += I_ERI_D2x_S_Gx3z_S_vrr;
      I_ERI_Dxy_S_Gx3z_S += I_ERI_Dxy_S_Gx3z_S_vrr;
      I_ERI_Dxz_S_Gx3z_S += I_ERI_Dxz_S_Gx3z_S_vrr;
      I_ERI_D2y_S_Gx3z_S += I_ERI_D2y_S_Gx3z_S_vrr;
      I_ERI_Dyz_S_Gx3z_S += I_ERI_Dyz_S_Gx3z_S_vrr;
      I_ERI_D2z_S_Gx3z_S += I_ERI_D2z_S_Gx3z_S_vrr;
      I_ERI_D2x_S_G4y_S += I_ERI_D2x_S_G4y_S_vrr;
      I_ERI_Dxy_S_G4y_S += I_ERI_Dxy_S_G4y_S_vrr;
      I_ERI_Dxz_S_G4y_S += I_ERI_Dxz_S_G4y_S_vrr;
      I_ERI_D2y_S_G4y_S += I_ERI_D2y_S_G4y_S_vrr;
      I_ERI_Dyz_S_G4y_S += I_ERI_Dyz_S_G4y_S_vrr;
      I_ERI_D2z_S_G4y_S += I_ERI_D2z_S_G4y_S_vrr;
      I_ERI_D2x_S_G3yz_S += I_ERI_D2x_S_G3yz_S_vrr;
      I_ERI_Dxy_S_G3yz_S += I_ERI_Dxy_S_G3yz_S_vrr;
      I_ERI_Dxz_S_G3yz_S += I_ERI_Dxz_S_G3yz_S_vrr;
      I_ERI_D2y_S_G3yz_S += I_ERI_D2y_S_G3yz_S_vrr;
      I_ERI_Dyz_S_G3yz_S += I_ERI_Dyz_S_G3yz_S_vrr;
      I_ERI_D2z_S_G3yz_S += I_ERI_D2z_S_G3yz_S_vrr;
      I_ERI_D2x_S_G2y2z_S += I_ERI_D2x_S_G2y2z_S_vrr;
      I_ERI_Dxy_S_G2y2z_S += I_ERI_Dxy_S_G2y2z_S_vrr;
      I_ERI_Dxz_S_G2y2z_S += I_ERI_Dxz_S_G2y2z_S_vrr;
      I_ERI_D2y_S_G2y2z_S += I_ERI_D2y_S_G2y2z_S_vrr;
      I_ERI_Dyz_S_G2y2z_S += I_ERI_Dyz_S_G2y2z_S_vrr;
      I_ERI_D2z_S_G2y2z_S += I_ERI_D2z_S_G2y2z_S_vrr;
      I_ERI_D2x_S_Gy3z_S += I_ERI_D2x_S_Gy3z_S_vrr;
      I_ERI_Dxy_S_Gy3z_S += I_ERI_Dxy_S_Gy3z_S_vrr;
      I_ERI_Dxz_S_Gy3z_S += I_ERI_Dxz_S_Gy3z_S_vrr;
      I_ERI_D2y_S_Gy3z_S += I_ERI_D2y_S_Gy3z_S_vrr;
      I_ERI_Dyz_S_Gy3z_S += I_ERI_Dyz_S_Gy3z_S_vrr;
      I_ERI_D2z_S_Gy3z_S += I_ERI_D2z_S_Gy3z_S_vrr;
      I_ERI_D2x_S_G4z_S += I_ERI_D2x_S_G4z_S_vrr;
      I_ERI_Dxy_S_G4z_S += I_ERI_Dxy_S_G4z_S_vrr;
      I_ERI_Dxz_S_G4z_S += I_ERI_Dxz_S_G4z_S_vrr;
      I_ERI_D2y_S_G4z_S += I_ERI_D2y_S_G4z_S_vrr;
      I_ERI_Dyz_S_G4z_S += I_ERI_Dyz_S_G4z_S_vrr;
      I_ERI_D2z_S_G4z_S += I_ERI_D2z_S_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_D2x_S_F3x_S += I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S += I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S += I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S += I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S += I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S += I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S += I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S += I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S += I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S += I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S += I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S += I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S += I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S += I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S += I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S += I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S += I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S += I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S += I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S += I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S += I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S += I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S += I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S += I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S += I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S += I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S += I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S += I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S += I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S += I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S += I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S += I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S += I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S += I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S += I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S += I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S += I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S += I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S += I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S += I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S += I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S += I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S += I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S += I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S += I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S += I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S += I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S += I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S += I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S += I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S += I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S += I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S += I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S += I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S += I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S += I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S += I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S += I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S += I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S += I_ERI_D2z_S_F3z_S_vrr;
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
   * shell quartet name: SQ_ERI_D_S_F_P
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_G_S
   * RHS shell quartet name: SQ_ERI_D_S_F_S
   ************************************************************/
  Double I_ERI_D2x_S_F3x_Px = I_ERI_D2x_S_G4x_S+CDX*I_ERI_D2x_S_F3x_S;
  Double I_ERI_Dxy_S_F3x_Px = I_ERI_Dxy_S_G4x_S+CDX*I_ERI_Dxy_S_F3x_S;
  Double I_ERI_Dxz_S_F3x_Px = I_ERI_Dxz_S_G4x_S+CDX*I_ERI_Dxz_S_F3x_S;
  Double I_ERI_D2y_S_F3x_Px = I_ERI_D2y_S_G4x_S+CDX*I_ERI_D2y_S_F3x_S;
  Double I_ERI_Dyz_S_F3x_Px = I_ERI_Dyz_S_G4x_S+CDX*I_ERI_Dyz_S_F3x_S;
  Double I_ERI_D2z_S_F3x_Px = I_ERI_D2z_S_G4x_S+CDX*I_ERI_D2z_S_F3x_S;
  Double I_ERI_D2x_S_F2xy_Px = I_ERI_D2x_S_G3xy_S+CDX*I_ERI_D2x_S_F2xy_S;
  Double I_ERI_Dxy_S_F2xy_Px = I_ERI_Dxy_S_G3xy_S+CDX*I_ERI_Dxy_S_F2xy_S;
  Double I_ERI_Dxz_S_F2xy_Px = I_ERI_Dxz_S_G3xy_S+CDX*I_ERI_Dxz_S_F2xy_S;
  Double I_ERI_D2y_S_F2xy_Px = I_ERI_D2y_S_G3xy_S+CDX*I_ERI_D2y_S_F2xy_S;
  Double I_ERI_Dyz_S_F2xy_Px = I_ERI_Dyz_S_G3xy_S+CDX*I_ERI_Dyz_S_F2xy_S;
  Double I_ERI_D2z_S_F2xy_Px = I_ERI_D2z_S_G3xy_S+CDX*I_ERI_D2z_S_F2xy_S;
  Double I_ERI_D2x_S_F2xz_Px = I_ERI_D2x_S_G3xz_S+CDX*I_ERI_D2x_S_F2xz_S;
  Double I_ERI_Dxy_S_F2xz_Px = I_ERI_Dxy_S_G3xz_S+CDX*I_ERI_Dxy_S_F2xz_S;
  Double I_ERI_Dxz_S_F2xz_Px = I_ERI_Dxz_S_G3xz_S+CDX*I_ERI_Dxz_S_F2xz_S;
  Double I_ERI_D2y_S_F2xz_Px = I_ERI_D2y_S_G3xz_S+CDX*I_ERI_D2y_S_F2xz_S;
  Double I_ERI_Dyz_S_F2xz_Px = I_ERI_Dyz_S_G3xz_S+CDX*I_ERI_Dyz_S_F2xz_S;
  Double I_ERI_D2z_S_F2xz_Px = I_ERI_D2z_S_G3xz_S+CDX*I_ERI_D2z_S_F2xz_S;
  Double I_ERI_D2x_S_Fx2y_Px = I_ERI_D2x_S_G2x2y_S+CDX*I_ERI_D2x_S_Fx2y_S;
  Double I_ERI_Dxy_S_Fx2y_Px = I_ERI_Dxy_S_G2x2y_S+CDX*I_ERI_Dxy_S_Fx2y_S;
  Double I_ERI_Dxz_S_Fx2y_Px = I_ERI_Dxz_S_G2x2y_S+CDX*I_ERI_Dxz_S_Fx2y_S;
  Double I_ERI_D2y_S_Fx2y_Px = I_ERI_D2y_S_G2x2y_S+CDX*I_ERI_D2y_S_Fx2y_S;
  Double I_ERI_Dyz_S_Fx2y_Px = I_ERI_Dyz_S_G2x2y_S+CDX*I_ERI_Dyz_S_Fx2y_S;
  Double I_ERI_D2z_S_Fx2y_Px = I_ERI_D2z_S_G2x2y_S+CDX*I_ERI_D2z_S_Fx2y_S;
  Double I_ERI_D2x_S_Fxyz_Px = I_ERI_D2x_S_G2xyz_S+CDX*I_ERI_D2x_S_Fxyz_S;
  Double I_ERI_Dxy_S_Fxyz_Px = I_ERI_Dxy_S_G2xyz_S+CDX*I_ERI_Dxy_S_Fxyz_S;
  Double I_ERI_Dxz_S_Fxyz_Px = I_ERI_Dxz_S_G2xyz_S+CDX*I_ERI_Dxz_S_Fxyz_S;
  Double I_ERI_D2y_S_Fxyz_Px = I_ERI_D2y_S_G2xyz_S+CDX*I_ERI_D2y_S_Fxyz_S;
  Double I_ERI_Dyz_S_Fxyz_Px = I_ERI_Dyz_S_G2xyz_S+CDX*I_ERI_Dyz_S_Fxyz_S;
  Double I_ERI_D2z_S_Fxyz_Px = I_ERI_D2z_S_G2xyz_S+CDX*I_ERI_D2z_S_Fxyz_S;
  Double I_ERI_D2x_S_Fx2z_Px = I_ERI_D2x_S_G2x2z_S+CDX*I_ERI_D2x_S_Fx2z_S;
  Double I_ERI_Dxy_S_Fx2z_Px = I_ERI_Dxy_S_G2x2z_S+CDX*I_ERI_Dxy_S_Fx2z_S;
  Double I_ERI_Dxz_S_Fx2z_Px = I_ERI_Dxz_S_G2x2z_S+CDX*I_ERI_Dxz_S_Fx2z_S;
  Double I_ERI_D2y_S_Fx2z_Px = I_ERI_D2y_S_G2x2z_S+CDX*I_ERI_D2y_S_Fx2z_S;
  Double I_ERI_Dyz_S_Fx2z_Px = I_ERI_Dyz_S_G2x2z_S+CDX*I_ERI_Dyz_S_Fx2z_S;
  Double I_ERI_D2z_S_Fx2z_Px = I_ERI_D2z_S_G2x2z_S+CDX*I_ERI_D2z_S_Fx2z_S;
  Double I_ERI_D2x_S_F3y_Px = I_ERI_D2x_S_Gx3y_S+CDX*I_ERI_D2x_S_F3y_S;
  Double I_ERI_Dxy_S_F3y_Px = I_ERI_Dxy_S_Gx3y_S+CDX*I_ERI_Dxy_S_F3y_S;
  Double I_ERI_Dxz_S_F3y_Px = I_ERI_Dxz_S_Gx3y_S+CDX*I_ERI_Dxz_S_F3y_S;
  Double I_ERI_D2y_S_F3y_Px = I_ERI_D2y_S_Gx3y_S+CDX*I_ERI_D2y_S_F3y_S;
  Double I_ERI_Dyz_S_F3y_Px = I_ERI_Dyz_S_Gx3y_S+CDX*I_ERI_Dyz_S_F3y_S;
  Double I_ERI_D2z_S_F3y_Px = I_ERI_D2z_S_Gx3y_S+CDX*I_ERI_D2z_S_F3y_S;
  Double I_ERI_D2x_S_F2yz_Px = I_ERI_D2x_S_Gx2yz_S+CDX*I_ERI_D2x_S_F2yz_S;
  Double I_ERI_Dxy_S_F2yz_Px = I_ERI_Dxy_S_Gx2yz_S+CDX*I_ERI_Dxy_S_F2yz_S;
  Double I_ERI_Dxz_S_F2yz_Px = I_ERI_Dxz_S_Gx2yz_S+CDX*I_ERI_Dxz_S_F2yz_S;
  Double I_ERI_D2y_S_F2yz_Px = I_ERI_D2y_S_Gx2yz_S+CDX*I_ERI_D2y_S_F2yz_S;
  Double I_ERI_Dyz_S_F2yz_Px = I_ERI_Dyz_S_Gx2yz_S+CDX*I_ERI_Dyz_S_F2yz_S;
  Double I_ERI_D2z_S_F2yz_Px = I_ERI_D2z_S_Gx2yz_S+CDX*I_ERI_D2z_S_F2yz_S;
  Double I_ERI_D2x_S_Fy2z_Px = I_ERI_D2x_S_Gxy2z_S+CDX*I_ERI_D2x_S_Fy2z_S;
  Double I_ERI_Dxy_S_Fy2z_Px = I_ERI_Dxy_S_Gxy2z_S+CDX*I_ERI_Dxy_S_Fy2z_S;
  Double I_ERI_Dxz_S_Fy2z_Px = I_ERI_Dxz_S_Gxy2z_S+CDX*I_ERI_Dxz_S_Fy2z_S;
  Double I_ERI_D2y_S_Fy2z_Px = I_ERI_D2y_S_Gxy2z_S+CDX*I_ERI_D2y_S_Fy2z_S;
  Double I_ERI_Dyz_S_Fy2z_Px = I_ERI_Dyz_S_Gxy2z_S+CDX*I_ERI_Dyz_S_Fy2z_S;
  Double I_ERI_D2z_S_Fy2z_Px = I_ERI_D2z_S_Gxy2z_S+CDX*I_ERI_D2z_S_Fy2z_S;
  Double I_ERI_D2x_S_F3z_Px = I_ERI_D2x_S_Gx3z_S+CDX*I_ERI_D2x_S_F3z_S;
  Double I_ERI_Dxy_S_F3z_Px = I_ERI_Dxy_S_Gx3z_S+CDX*I_ERI_Dxy_S_F3z_S;
  Double I_ERI_Dxz_S_F3z_Px = I_ERI_Dxz_S_Gx3z_S+CDX*I_ERI_Dxz_S_F3z_S;
  Double I_ERI_D2y_S_F3z_Px = I_ERI_D2y_S_Gx3z_S+CDX*I_ERI_D2y_S_F3z_S;
  Double I_ERI_Dyz_S_F3z_Px = I_ERI_Dyz_S_Gx3z_S+CDX*I_ERI_Dyz_S_F3z_S;
  Double I_ERI_D2z_S_F3z_Px = I_ERI_D2z_S_Gx3z_S+CDX*I_ERI_D2z_S_F3z_S;
  Double I_ERI_D2x_S_F3x_Py = I_ERI_D2x_S_G3xy_S+CDY*I_ERI_D2x_S_F3x_S;
  Double I_ERI_Dxy_S_F3x_Py = I_ERI_Dxy_S_G3xy_S+CDY*I_ERI_Dxy_S_F3x_S;
  Double I_ERI_Dxz_S_F3x_Py = I_ERI_Dxz_S_G3xy_S+CDY*I_ERI_Dxz_S_F3x_S;
  Double I_ERI_D2y_S_F3x_Py = I_ERI_D2y_S_G3xy_S+CDY*I_ERI_D2y_S_F3x_S;
  Double I_ERI_Dyz_S_F3x_Py = I_ERI_Dyz_S_G3xy_S+CDY*I_ERI_Dyz_S_F3x_S;
  Double I_ERI_D2z_S_F3x_Py = I_ERI_D2z_S_G3xy_S+CDY*I_ERI_D2z_S_F3x_S;
  Double I_ERI_D2x_S_F2xy_Py = I_ERI_D2x_S_G2x2y_S+CDY*I_ERI_D2x_S_F2xy_S;
  Double I_ERI_Dxy_S_F2xy_Py = I_ERI_Dxy_S_G2x2y_S+CDY*I_ERI_Dxy_S_F2xy_S;
  Double I_ERI_Dxz_S_F2xy_Py = I_ERI_Dxz_S_G2x2y_S+CDY*I_ERI_Dxz_S_F2xy_S;
  Double I_ERI_D2y_S_F2xy_Py = I_ERI_D2y_S_G2x2y_S+CDY*I_ERI_D2y_S_F2xy_S;
  Double I_ERI_Dyz_S_F2xy_Py = I_ERI_Dyz_S_G2x2y_S+CDY*I_ERI_Dyz_S_F2xy_S;
  Double I_ERI_D2z_S_F2xy_Py = I_ERI_D2z_S_G2x2y_S+CDY*I_ERI_D2z_S_F2xy_S;
  Double I_ERI_D2x_S_F2xz_Py = I_ERI_D2x_S_G2xyz_S+CDY*I_ERI_D2x_S_F2xz_S;
  Double I_ERI_Dxy_S_F2xz_Py = I_ERI_Dxy_S_G2xyz_S+CDY*I_ERI_Dxy_S_F2xz_S;
  Double I_ERI_Dxz_S_F2xz_Py = I_ERI_Dxz_S_G2xyz_S+CDY*I_ERI_Dxz_S_F2xz_S;
  Double I_ERI_D2y_S_F2xz_Py = I_ERI_D2y_S_G2xyz_S+CDY*I_ERI_D2y_S_F2xz_S;
  Double I_ERI_Dyz_S_F2xz_Py = I_ERI_Dyz_S_G2xyz_S+CDY*I_ERI_Dyz_S_F2xz_S;
  Double I_ERI_D2z_S_F2xz_Py = I_ERI_D2z_S_G2xyz_S+CDY*I_ERI_D2z_S_F2xz_S;
  Double I_ERI_D2x_S_Fx2y_Py = I_ERI_D2x_S_Gx3y_S+CDY*I_ERI_D2x_S_Fx2y_S;
  Double I_ERI_Dxy_S_Fx2y_Py = I_ERI_Dxy_S_Gx3y_S+CDY*I_ERI_Dxy_S_Fx2y_S;
  Double I_ERI_Dxz_S_Fx2y_Py = I_ERI_Dxz_S_Gx3y_S+CDY*I_ERI_Dxz_S_Fx2y_S;
  Double I_ERI_D2y_S_Fx2y_Py = I_ERI_D2y_S_Gx3y_S+CDY*I_ERI_D2y_S_Fx2y_S;
  Double I_ERI_Dyz_S_Fx2y_Py = I_ERI_Dyz_S_Gx3y_S+CDY*I_ERI_Dyz_S_Fx2y_S;
  Double I_ERI_D2z_S_Fx2y_Py = I_ERI_D2z_S_Gx3y_S+CDY*I_ERI_D2z_S_Fx2y_S;
  Double I_ERI_D2x_S_Fxyz_Py = I_ERI_D2x_S_Gx2yz_S+CDY*I_ERI_D2x_S_Fxyz_S;
  Double I_ERI_Dxy_S_Fxyz_Py = I_ERI_Dxy_S_Gx2yz_S+CDY*I_ERI_Dxy_S_Fxyz_S;
  Double I_ERI_Dxz_S_Fxyz_Py = I_ERI_Dxz_S_Gx2yz_S+CDY*I_ERI_Dxz_S_Fxyz_S;
  Double I_ERI_D2y_S_Fxyz_Py = I_ERI_D2y_S_Gx2yz_S+CDY*I_ERI_D2y_S_Fxyz_S;
  Double I_ERI_Dyz_S_Fxyz_Py = I_ERI_Dyz_S_Gx2yz_S+CDY*I_ERI_Dyz_S_Fxyz_S;
  Double I_ERI_D2z_S_Fxyz_Py = I_ERI_D2z_S_Gx2yz_S+CDY*I_ERI_D2z_S_Fxyz_S;
  Double I_ERI_D2x_S_Fx2z_Py = I_ERI_D2x_S_Gxy2z_S+CDY*I_ERI_D2x_S_Fx2z_S;
  Double I_ERI_Dxy_S_Fx2z_Py = I_ERI_Dxy_S_Gxy2z_S+CDY*I_ERI_Dxy_S_Fx2z_S;
  Double I_ERI_Dxz_S_Fx2z_Py = I_ERI_Dxz_S_Gxy2z_S+CDY*I_ERI_Dxz_S_Fx2z_S;
  Double I_ERI_D2y_S_Fx2z_Py = I_ERI_D2y_S_Gxy2z_S+CDY*I_ERI_D2y_S_Fx2z_S;
  Double I_ERI_Dyz_S_Fx2z_Py = I_ERI_Dyz_S_Gxy2z_S+CDY*I_ERI_Dyz_S_Fx2z_S;
  Double I_ERI_D2z_S_Fx2z_Py = I_ERI_D2z_S_Gxy2z_S+CDY*I_ERI_D2z_S_Fx2z_S;
  Double I_ERI_D2x_S_F3y_Py = I_ERI_D2x_S_G4y_S+CDY*I_ERI_D2x_S_F3y_S;
  Double I_ERI_Dxy_S_F3y_Py = I_ERI_Dxy_S_G4y_S+CDY*I_ERI_Dxy_S_F3y_S;
  Double I_ERI_Dxz_S_F3y_Py = I_ERI_Dxz_S_G4y_S+CDY*I_ERI_Dxz_S_F3y_S;
  Double I_ERI_D2y_S_F3y_Py = I_ERI_D2y_S_G4y_S+CDY*I_ERI_D2y_S_F3y_S;
  Double I_ERI_Dyz_S_F3y_Py = I_ERI_Dyz_S_G4y_S+CDY*I_ERI_Dyz_S_F3y_S;
  Double I_ERI_D2z_S_F3y_Py = I_ERI_D2z_S_G4y_S+CDY*I_ERI_D2z_S_F3y_S;
  Double I_ERI_D2x_S_F2yz_Py = I_ERI_D2x_S_G3yz_S+CDY*I_ERI_D2x_S_F2yz_S;
  Double I_ERI_Dxy_S_F2yz_Py = I_ERI_Dxy_S_G3yz_S+CDY*I_ERI_Dxy_S_F2yz_S;
  Double I_ERI_Dxz_S_F2yz_Py = I_ERI_Dxz_S_G3yz_S+CDY*I_ERI_Dxz_S_F2yz_S;
  Double I_ERI_D2y_S_F2yz_Py = I_ERI_D2y_S_G3yz_S+CDY*I_ERI_D2y_S_F2yz_S;
  Double I_ERI_Dyz_S_F2yz_Py = I_ERI_Dyz_S_G3yz_S+CDY*I_ERI_Dyz_S_F2yz_S;
  Double I_ERI_D2z_S_F2yz_Py = I_ERI_D2z_S_G3yz_S+CDY*I_ERI_D2z_S_F2yz_S;
  Double I_ERI_D2x_S_Fy2z_Py = I_ERI_D2x_S_G2y2z_S+CDY*I_ERI_D2x_S_Fy2z_S;
  Double I_ERI_Dxy_S_Fy2z_Py = I_ERI_Dxy_S_G2y2z_S+CDY*I_ERI_Dxy_S_Fy2z_S;
  Double I_ERI_Dxz_S_Fy2z_Py = I_ERI_Dxz_S_G2y2z_S+CDY*I_ERI_Dxz_S_Fy2z_S;
  Double I_ERI_D2y_S_Fy2z_Py = I_ERI_D2y_S_G2y2z_S+CDY*I_ERI_D2y_S_Fy2z_S;
  Double I_ERI_Dyz_S_Fy2z_Py = I_ERI_Dyz_S_G2y2z_S+CDY*I_ERI_Dyz_S_Fy2z_S;
  Double I_ERI_D2z_S_Fy2z_Py = I_ERI_D2z_S_G2y2z_S+CDY*I_ERI_D2z_S_Fy2z_S;
  Double I_ERI_D2x_S_F3z_Py = I_ERI_D2x_S_Gy3z_S+CDY*I_ERI_D2x_S_F3z_S;
  Double I_ERI_Dxy_S_F3z_Py = I_ERI_Dxy_S_Gy3z_S+CDY*I_ERI_Dxy_S_F3z_S;
  Double I_ERI_Dxz_S_F3z_Py = I_ERI_Dxz_S_Gy3z_S+CDY*I_ERI_Dxz_S_F3z_S;
  Double I_ERI_D2y_S_F3z_Py = I_ERI_D2y_S_Gy3z_S+CDY*I_ERI_D2y_S_F3z_S;
  Double I_ERI_Dyz_S_F3z_Py = I_ERI_Dyz_S_Gy3z_S+CDY*I_ERI_Dyz_S_F3z_S;
  Double I_ERI_D2z_S_F3z_Py = I_ERI_D2z_S_Gy3z_S+CDY*I_ERI_D2z_S_F3z_S;
  Double I_ERI_D2x_S_F3x_Pz = I_ERI_D2x_S_G3xz_S+CDZ*I_ERI_D2x_S_F3x_S;
  Double I_ERI_Dxy_S_F3x_Pz = I_ERI_Dxy_S_G3xz_S+CDZ*I_ERI_Dxy_S_F3x_S;
  Double I_ERI_Dxz_S_F3x_Pz = I_ERI_Dxz_S_G3xz_S+CDZ*I_ERI_Dxz_S_F3x_S;
  Double I_ERI_D2y_S_F3x_Pz = I_ERI_D2y_S_G3xz_S+CDZ*I_ERI_D2y_S_F3x_S;
  Double I_ERI_Dyz_S_F3x_Pz = I_ERI_Dyz_S_G3xz_S+CDZ*I_ERI_Dyz_S_F3x_S;
  Double I_ERI_D2z_S_F3x_Pz = I_ERI_D2z_S_G3xz_S+CDZ*I_ERI_D2z_S_F3x_S;
  Double I_ERI_D2x_S_F2xy_Pz = I_ERI_D2x_S_G2xyz_S+CDZ*I_ERI_D2x_S_F2xy_S;
  Double I_ERI_Dxy_S_F2xy_Pz = I_ERI_Dxy_S_G2xyz_S+CDZ*I_ERI_Dxy_S_F2xy_S;
  Double I_ERI_Dxz_S_F2xy_Pz = I_ERI_Dxz_S_G2xyz_S+CDZ*I_ERI_Dxz_S_F2xy_S;
  Double I_ERI_D2y_S_F2xy_Pz = I_ERI_D2y_S_G2xyz_S+CDZ*I_ERI_D2y_S_F2xy_S;
  Double I_ERI_Dyz_S_F2xy_Pz = I_ERI_Dyz_S_G2xyz_S+CDZ*I_ERI_Dyz_S_F2xy_S;
  Double I_ERI_D2z_S_F2xy_Pz = I_ERI_D2z_S_G2xyz_S+CDZ*I_ERI_D2z_S_F2xy_S;
  Double I_ERI_D2x_S_F2xz_Pz = I_ERI_D2x_S_G2x2z_S+CDZ*I_ERI_D2x_S_F2xz_S;
  Double I_ERI_Dxy_S_F2xz_Pz = I_ERI_Dxy_S_G2x2z_S+CDZ*I_ERI_Dxy_S_F2xz_S;
  Double I_ERI_Dxz_S_F2xz_Pz = I_ERI_Dxz_S_G2x2z_S+CDZ*I_ERI_Dxz_S_F2xz_S;
  Double I_ERI_D2y_S_F2xz_Pz = I_ERI_D2y_S_G2x2z_S+CDZ*I_ERI_D2y_S_F2xz_S;
  Double I_ERI_Dyz_S_F2xz_Pz = I_ERI_Dyz_S_G2x2z_S+CDZ*I_ERI_Dyz_S_F2xz_S;
  Double I_ERI_D2z_S_F2xz_Pz = I_ERI_D2z_S_G2x2z_S+CDZ*I_ERI_D2z_S_F2xz_S;
  Double I_ERI_D2x_S_Fx2y_Pz = I_ERI_D2x_S_Gx2yz_S+CDZ*I_ERI_D2x_S_Fx2y_S;
  Double I_ERI_Dxy_S_Fx2y_Pz = I_ERI_Dxy_S_Gx2yz_S+CDZ*I_ERI_Dxy_S_Fx2y_S;
  Double I_ERI_Dxz_S_Fx2y_Pz = I_ERI_Dxz_S_Gx2yz_S+CDZ*I_ERI_Dxz_S_Fx2y_S;
  Double I_ERI_D2y_S_Fx2y_Pz = I_ERI_D2y_S_Gx2yz_S+CDZ*I_ERI_D2y_S_Fx2y_S;
  Double I_ERI_Dyz_S_Fx2y_Pz = I_ERI_Dyz_S_Gx2yz_S+CDZ*I_ERI_Dyz_S_Fx2y_S;
  Double I_ERI_D2z_S_Fx2y_Pz = I_ERI_D2z_S_Gx2yz_S+CDZ*I_ERI_D2z_S_Fx2y_S;
  Double I_ERI_D2x_S_Fxyz_Pz = I_ERI_D2x_S_Gxy2z_S+CDZ*I_ERI_D2x_S_Fxyz_S;
  Double I_ERI_Dxy_S_Fxyz_Pz = I_ERI_Dxy_S_Gxy2z_S+CDZ*I_ERI_Dxy_S_Fxyz_S;
  Double I_ERI_Dxz_S_Fxyz_Pz = I_ERI_Dxz_S_Gxy2z_S+CDZ*I_ERI_Dxz_S_Fxyz_S;
  Double I_ERI_D2y_S_Fxyz_Pz = I_ERI_D2y_S_Gxy2z_S+CDZ*I_ERI_D2y_S_Fxyz_S;
  Double I_ERI_Dyz_S_Fxyz_Pz = I_ERI_Dyz_S_Gxy2z_S+CDZ*I_ERI_Dyz_S_Fxyz_S;
  Double I_ERI_D2z_S_Fxyz_Pz = I_ERI_D2z_S_Gxy2z_S+CDZ*I_ERI_D2z_S_Fxyz_S;
  Double I_ERI_D2x_S_Fx2z_Pz = I_ERI_D2x_S_Gx3z_S+CDZ*I_ERI_D2x_S_Fx2z_S;
  Double I_ERI_Dxy_S_Fx2z_Pz = I_ERI_Dxy_S_Gx3z_S+CDZ*I_ERI_Dxy_S_Fx2z_S;
  Double I_ERI_Dxz_S_Fx2z_Pz = I_ERI_Dxz_S_Gx3z_S+CDZ*I_ERI_Dxz_S_Fx2z_S;
  Double I_ERI_D2y_S_Fx2z_Pz = I_ERI_D2y_S_Gx3z_S+CDZ*I_ERI_D2y_S_Fx2z_S;
  Double I_ERI_Dyz_S_Fx2z_Pz = I_ERI_Dyz_S_Gx3z_S+CDZ*I_ERI_Dyz_S_Fx2z_S;
  Double I_ERI_D2z_S_Fx2z_Pz = I_ERI_D2z_S_Gx3z_S+CDZ*I_ERI_D2z_S_Fx2z_S;
  Double I_ERI_D2x_S_F3y_Pz = I_ERI_D2x_S_G3yz_S+CDZ*I_ERI_D2x_S_F3y_S;
  Double I_ERI_Dxy_S_F3y_Pz = I_ERI_Dxy_S_G3yz_S+CDZ*I_ERI_Dxy_S_F3y_S;
  Double I_ERI_Dxz_S_F3y_Pz = I_ERI_Dxz_S_G3yz_S+CDZ*I_ERI_Dxz_S_F3y_S;
  Double I_ERI_D2y_S_F3y_Pz = I_ERI_D2y_S_G3yz_S+CDZ*I_ERI_D2y_S_F3y_S;
  Double I_ERI_Dyz_S_F3y_Pz = I_ERI_Dyz_S_G3yz_S+CDZ*I_ERI_Dyz_S_F3y_S;
  Double I_ERI_D2z_S_F3y_Pz = I_ERI_D2z_S_G3yz_S+CDZ*I_ERI_D2z_S_F3y_S;
  Double I_ERI_D2x_S_F2yz_Pz = I_ERI_D2x_S_G2y2z_S+CDZ*I_ERI_D2x_S_F2yz_S;
  Double I_ERI_Dxy_S_F2yz_Pz = I_ERI_Dxy_S_G2y2z_S+CDZ*I_ERI_Dxy_S_F2yz_S;
  Double I_ERI_Dxz_S_F2yz_Pz = I_ERI_Dxz_S_G2y2z_S+CDZ*I_ERI_Dxz_S_F2yz_S;
  Double I_ERI_D2y_S_F2yz_Pz = I_ERI_D2y_S_G2y2z_S+CDZ*I_ERI_D2y_S_F2yz_S;
  Double I_ERI_Dyz_S_F2yz_Pz = I_ERI_Dyz_S_G2y2z_S+CDZ*I_ERI_Dyz_S_F2yz_S;
  Double I_ERI_D2z_S_F2yz_Pz = I_ERI_D2z_S_G2y2z_S+CDZ*I_ERI_D2z_S_F2yz_S;
  Double I_ERI_D2x_S_Fy2z_Pz = I_ERI_D2x_S_Gy3z_S+CDZ*I_ERI_D2x_S_Fy2z_S;
  Double I_ERI_Dxy_S_Fy2z_Pz = I_ERI_Dxy_S_Gy3z_S+CDZ*I_ERI_Dxy_S_Fy2z_S;
  Double I_ERI_Dxz_S_Fy2z_Pz = I_ERI_Dxz_S_Gy3z_S+CDZ*I_ERI_Dxz_S_Fy2z_S;
  Double I_ERI_D2y_S_Fy2z_Pz = I_ERI_D2y_S_Gy3z_S+CDZ*I_ERI_D2y_S_Fy2z_S;
  Double I_ERI_Dyz_S_Fy2z_Pz = I_ERI_Dyz_S_Gy3z_S+CDZ*I_ERI_Dyz_S_Fy2z_S;
  Double I_ERI_D2z_S_Fy2z_Pz = I_ERI_D2z_S_Gy3z_S+CDZ*I_ERI_D2z_S_Fy2z_S;
  Double I_ERI_D2x_S_F3z_Pz = I_ERI_D2x_S_G4z_S+CDZ*I_ERI_D2x_S_F3z_S;
  Double I_ERI_Dxy_S_F3z_Pz = I_ERI_Dxy_S_G4z_S+CDZ*I_ERI_Dxy_S_F3z_S;
  Double I_ERI_Dxz_S_F3z_Pz = I_ERI_Dxz_S_G4z_S+CDZ*I_ERI_Dxz_S_F3z_S;
  Double I_ERI_D2y_S_F3z_Pz = I_ERI_D2y_S_G4z_S+CDZ*I_ERI_D2y_S_F3z_S;
  Double I_ERI_Dyz_S_F3z_Pz = I_ERI_Dyz_S_G4z_S+CDZ*I_ERI_Dyz_S_F3z_S;
  Double I_ERI_D2z_S_F3z_Pz = I_ERI_D2z_S_G4z_S+CDZ*I_ERI_D2z_S_F3z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_F_P
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_G_S
   * RHS shell quartet name: SQ_ERI_F_S_F_S
   ************************************************************/
  Double I_ERI_F3x_S_F3x_Px = I_ERI_F3x_S_G4x_S+CDX*I_ERI_F3x_S_F3x_S;
  Double I_ERI_F2xy_S_F3x_Px = I_ERI_F2xy_S_G4x_S+CDX*I_ERI_F2xy_S_F3x_S;
  Double I_ERI_F2xz_S_F3x_Px = I_ERI_F2xz_S_G4x_S+CDX*I_ERI_F2xz_S_F3x_S;
  Double I_ERI_Fx2y_S_F3x_Px = I_ERI_Fx2y_S_G4x_S+CDX*I_ERI_Fx2y_S_F3x_S;
  Double I_ERI_Fxyz_S_F3x_Px = I_ERI_Fxyz_S_G4x_S+CDX*I_ERI_Fxyz_S_F3x_S;
  Double I_ERI_Fx2z_S_F3x_Px = I_ERI_Fx2z_S_G4x_S+CDX*I_ERI_Fx2z_S_F3x_S;
  Double I_ERI_F3y_S_F3x_Px = I_ERI_F3y_S_G4x_S+CDX*I_ERI_F3y_S_F3x_S;
  Double I_ERI_F2yz_S_F3x_Px = I_ERI_F2yz_S_G4x_S+CDX*I_ERI_F2yz_S_F3x_S;
  Double I_ERI_Fy2z_S_F3x_Px = I_ERI_Fy2z_S_G4x_S+CDX*I_ERI_Fy2z_S_F3x_S;
  Double I_ERI_F3z_S_F3x_Px = I_ERI_F3z_S_G4x_S+CDX*I_ERI_F3z_S_F3x_S;
  Double I_ERI_F3x_S_F2xy_Px = I_ERI_F3x_S_G3xy_S+CDX*I_ERI_F3x_S_F2xy_S;
  Double I_ERI_F2xy_S_F2xy_Px = I_ERI_F2xy_S_G3xy_S+CDX*I_ERI_F2xy_S_F2xy_S;
  Double I_ERI_F2xz_S_F2xy_Px = I_ERI_F2xz_S_G3xy_S+CDX*I_ERI_F2xz_S_F2xy_S;
  Double I_ERI_Fx2y_S_F2xy_Px = I_ERI_Fx2y_S_G3xy_S+CDX*I_ERI_Fx2y_S_F2xy_S;
  Double I_ERI_Fxyz_S_F2xy_Px = I_ERI_Fxyz_S_G3xy_S+CDX*I_ERI_Fxyz_S_F2xy_S;
  Double I_ERI_Fx2z_S_F2xy_Px = I_ERI_Fx2z_S_G3xy_S+CDX*I_ERI_Fx2z_S_F2xy_S;
  Double I_ERI_F3y_S_F2xy_Px = I_ERI_F3y_S_G3xy_S+CDX*I_ERI_F3y_S_F2xy_S;
  Double I_ERI_F2yz_S_F2xy_Px = I_ERI_F2yz_S_G3xy_S+CDX*I_ERI_F2yz_S_F2xy_S;
  Double I_ERI_Fy2z_S_F2xy_Px = I_ERI_Fy2z_S_G3xy_S+CDX*I_ERI_Fy2z_S_F2xy_S;
  Double I_ERI_F3z_S_F2xy_Px = I_ERI_F3z_S_G3xy_S+CDX*I_ERI_F3z_S_F2xy_S;
  Double I_ERI_F3x_S_F2xz_Px = I_ERI_F3x_S_G3xz_S+CDX*I_ERI_F3x_S_F2xz_S;
  Double I_ERI_F2xy_S_F2xz_Px = I_ERI_F2xy_S_G3xz_S+CDX*I_ERI_F2xy_S_F2xz_S;
  Double I_ERI_F2xz_S_F2xz_Px = I_ERI_F2xz_S_G3xz_S+CDX*I_ERI_F2xz_S_F2xz_S;
  Double I_ERI_Fx2y_S_F2xz_Px = I_ERI_Fx2y_S_G3xz_S+CDX*I_ERI_Fx2y_S_F2xz_S;
  Double I_ERI_Fxyz_S_F2xz_Px = I_ERI_Fxyz_S_G3xz_S+CDX*I_ERI_Fxyz_S_F2xz_S;
  Double I_ERI_Fx2z_S_F2xz_Px = I_ERI_Fx2z_S_G3xz_S+CDX*I_ERI_Fx2z_S_F2xz_S;
  Double I_ERI_F3y_S_F2xz_Px = I_ERI_F3y_S_G3xz_S+CDX*I_ERI_F3y_S_F2xz_S;
  Double I_ERI_F2yz_S_F2xz_Px = I_ERI_F2yz_S_G3xz_S+CDX*I_ERI_F2yz_S_F2xz_S;
  Double I_ERI_Fy2z_S_F2xz_Px = I_ERI_Fy2z_S_G3xz_S+CDX*I_ERI_Fy2z_S_F2xz_S;
  Double I_ERI_F3z_S_F2xz_Px = I_ERI_F3z_S_G3xz_S+CDX*I_ERI_F3z_S_F2xz_S;
  Double I_ERI_F3x_S_Fx2y_Px = I_ERI_F3x_S_G2x2y_S+CDX*I_ERI_F3x_S_Fx2y_S;
  Double I_ERI_F2xy_S_Fx2y_Px = I_ERI_F2xy_S_G2x2y_S+CDX*I_ERI_F2xy_S_Fx2y_S;
  Double I_ERI_F2xz_S_Fx2y_Px = I_ERI_F2xz_S_G2x2y_S+CDX*I_ERI_F2xz_S_Fx2y_S;
  Double I_ERI_Fx2y_S_Fx2y_Px = I_ERI_Fx2y_S_G2x2y_S+CDX*I_ERI_Fx2y_S_Fx2y_S;
  Double I_ERI_Fxyz_S_Fx2y_Px = I_ERI_Fxyz_S_G2x2y_S+CDX*I_ERI_Fxyz_S_Fx2y_S;
  Double I_ERI_Fx2z_S_Fx2y_Px = I_ERI_Fx2z_S_G2x2y_S+CDX*I_ERI_Fx2z_S_Fx2y_S;
  Double I_ERI_F3y_S_Fx2y_Px = I_ERI_F3y_S_G2x2y_S+CDX*I_ERI_F3y_S_Fx2y_S;
  Double I_ERI_F2yz_S_Fx2y_Px = I_ERI_F2yz_S_G2x2y_S+CDX*I_ERI_F2yz_S_Fx2y_S;
  Double I_ERI_Fy2z_S_Fx2y_Px = I_ERI_Fy2z_S_G2x2y_S+CDX*I_ERI_Fy2z_S_Fx2y_S;
  Double I_ERI_F3z_S_Fx2y_Px = I_ERI_F3z_S_G2x2y_S+CDX*I_ERI_F3z_S_Fx2y_S;
  Double I_ERI_F3x_S_Fxyz_Px = I_ERI_F3x_S_G2xyz_S+CDX*I_ERI_F3x_S_Fxyz_S;
  Double I_ERI_F2xy_S_Fxyz_Px = I_ERI_F2xy_S_G2xyz_S+CDX*I_ERI_F2xy_S_Fxyz_S;
  Double I_ERI_F2xz_S_Fxyz_Px = I_ERI_F2xz_S_G2xyz_S+CDX*I_ERI_F2xz_S_Fxyz_S;
  Double I_ERI_Fx2y_S_Fxyz_Px = I_ERI_Fx2y_S_G2xyz_S+CDX*I_ERI_Fx2y_S_Fxyz_S;
  Double I_ERI_Fxyz_S_Fxyz_Px = I_ERI_Fxyz_S_G2xyz_S+CDX*I_ERI_Fxyz_S_Fxyz_S;
  Double I_ERI_Fx2z_S_Fxyz_Px = I_ERI_Fx2z_S_G2xyz_S+CDX*I_ERI_Fx2z_S_Fxyz_S;
  Double I_ERI_F3y_S_Fxyz_Px = I_ERI_F3y_S_G2xyz_S+CDX*I_ERI_F3y_S_Fxyz_S;
  Double I_ERI_F2yz_S_Fxyz_Px = I_ERI_F2yz_S_G2xyz_S+CDX*I_ERI_F2yz_S_Fxyz_S;
  Double I_ERI_Fy2z_S_Fxyz_Px = I_ERI_Fy2z_S_G2xyz_S+CDX*I_ERI_Fy2z_S_Fxyz_S;
  Double I_ERI_F3z_S_Fxyz_Px = I_ERI_F3z_S_G2xyz_S+CDX*I_ERI_F3z_S_Fxyz_S;
  Double I_ERI_F3x_S_Fx2z_Px = I_ERI_F3x_S_G2x2z_S+CDX*I_ERI_F3x_S_Fx2z_S;
  Double I_ERI_F2xy_S_Fx2z_Px = I_ERI_F2xy_S_G2x2z_S+CDX*I_ERI_F2xy_S_Fx2z_S;
  Double I_ERI_F2xz_S_Fx2z_Px = I_ERI_F2xz_S_G2x2z_S+CDX*I_ERI_F2xz_S_Fx2z_S;
  Double I_ERI_Fx2y_S_Fx2z_Px = I_ERI_Fx2y_S_G2x2z_S+CDX*I_ERI_Fx2y_S_Fx2z_S;
  Double I_ERI_Fxyz_S_Fx2z_Px = I_ERI_Fxyz_S_G2x2z_S+CDX*I_ERI_Fxyz_S_Fx2z_S;
  Double I_ERI_Fx2z_S_Fx2z_Px = I_ERI_Fx2z_S_G2x2z_S+CDX*I_ERI_Fx2z_S_Fx2z_S;
  Double I_ERI_F3y_S_Fx2z_Px = I_ERI_F3y_S_G2x2z_S+CDX*I_ERI_F3y_S_Fx2z_S;
  Double I_ERI_F2yz_S_Fx2z_Px = I_ERI_F2yz_S_G2x2z_S+CDX*I_ERI_F2yz_S_Fx2z_S;
  Double I_ERI_Fy2z_S_Fx2z_Px = I_ERI_Fy2z_S_G2x2z_S+CDX*I_ERI_Fy2z_S_Fx2z_S;
  Double I_ERI_F3z_S_Fx2z_Px = I_ERI_F3z_S_G2x2z_S+CDX*I_ERI_F3z_S_Fx2z_S;
  Double I_ERI_F3x_S_F3y_Px = I_ERI_F3x_S_Gx3y_S+CDX*I_ERI_F3x_S_F3y_S;
  Double I_ERI_F2xy_S_F3y_Px = I_ERI_F2xy_S_Gx3y_S+CDX*I_ERI_F2xy_S_F3y_S;
  Double I_ERI_F2xz_S_F3y_Px = I_ERI_F2xz_S_Gx3y_S+CDX*I_ERI_F2xz_S_F3y_S;
  Double I_ERI_Fx2y_S_F3y_Px = I_ERI_Fx2y_S_Gx3y_S+CDX*I_ERI_Fx2y_S_F3y_S;
  Double I_ERI_Fxyz_S_F3y_Px = I_ERI_Fxyz_S_Gx3y_S+CDX*I_ERI_Fxyz_S_F3y_S;
  Double I_ERI_Fx2z_S_F3y_Px = I_ERI_Fx2z_S_Gx3y_S+CDX*I_ERI_Fx2z_S_F3y_S;
  Double I_ERI_F3y_S_F3y_Px = I_ERI_F3y_S_Gx3y_S+CDX*I_ERI_F3y_S_F3y_S;
  Double I_ERI_F2yz_S_F3y_Px = I_ERI_F2yz_S_Gx3y_S+CDX*I_ERI_F2yz_S_F3y_S;
  Double I_ERI_Fy2z_S_F3y_Px = I_ERI_Fy2z_S_Gx3y_S+CDX*I_ERI_Fy2z_S_F3y_S;
  Double I_ERI_F3z_S_F3y_Px = I_ERI_F3z_S_Gx3y_S+CDX*I_ERI_F3z_S_F3y_S;
  Double I_ERI_F3x_S_F2yz_Px = I_ERI_F3x_S_Gx2yz_S+CDX*I_ERI_F3x_S_F2yz_S;
  Double I_ERI_F2xy_S_F2yz_Px = I_ERI_F2xy_S_Gx2yz_S+CDX*I_ERI_F2xy_S_F2yz_S;
  Double I_ERI_F2xz_S_F2yz_Px = I_ERI_F2xz_S_Gx2yz_S+CDX*I_ERI_F2xz_S_F2yz_S;
  Double I_ERI_Fx2y_S_F2yz_Px = I_ERI_Fx2y_S_Gx2yz_S+CDX*I_ERI_Fx2y_S_F2yz_S;
  Double I_ERI_Fxyz_S_F2yz_Px = I_ERI_Fxyz_S_Gx2yz_S+CDX*I_ERI_Fxyz_S_F2yz_S;
  Double I_ERI_Fx2z_S_F2yz_Px = I_ERI_Fx2z_S_Gx2yz_S+CDX*I_ERI_Fx2z_S_F2yz_S;
  Double I_ERI_F3y_S_F2yz_Px = I_ERI_F3y_S_Gx2yz_S+CDX*I_ERI_F3y_S_F2yz_S;
  Double I_ERI_F2yz_S_F2yz_Px = I_ERI_F2yz_S_Gx2yz_S+CDX*I_ERI_F2yz_S_F2yz_S;
  Double I_ERI_Fy2z_S_F2yz_Px = I_ERI_Fy2z_S_Gx2yz_S+CDX*I_ERI_Fy2z_S_F2yz_S;
  Double I_ERI_F3z_S_F2yz_Px = I_ERI_F3z_S_Gx2yz_S+CDX*I_ERI_F3z_S_F2yz_S;
  Double I_ERI_F3x_S_Fy2z_Px = I_ERI_F3x_S_Gxy2z_S+CDX*I_ERI_F3x_S_Fy2z_S;
  Double I_ERI_F2xy_S_Fy2z_Px = I_ERI_F2xy_S_Gxy2z_S+CDX*I_ERI_F2xy_S_Fy2z_S;
  Double I_ERI_F2xz_S_Fy2z_Px = I_ERI_F2xz_S_Gxy2z_S+CDX*I_ERI_F2xz_S_Fy2z_S;
  Double I_ERI_Fx2y_S_Fy2z_Px = I_ERI_Fx2y_S_Gxy2z_S+CDX*I_ERI_Fx2y_S_Fy2z_S;
  Double I_ERI_Fxyz_S_Fy2z_Px = I_ERI_Fxyz_S_Gxy2z_S+CDX*I_ERI_Fxyz_S_Fy2z_S;
  Double I_ERI_Fx2z_S_Fy2z_Px = I_ERI_Fx2z_S_Gxy2z_S+CDX*I_ERI_Fx2z_S_Fy2z_S;
  Double I_ERI_F3y_S_Fy2z_Px = I_ERI_F3y_S_Gxy2z_S+CDX*I_ERI_F3y_S_Fy2z_S;
  Double I_ERI_F2yz_S_Fy2z_Px = I_ERI_F2yz_S_Gxy2z_S+CDX*I_ERI_F2yz_S_Fy2z_S;
  Double I_ERI_Fy2z_S_Fy2z_Px = I_ERI_Fy2z_S_Gxy2z_S+CDX*I_ERI_Fy2z_S_Fy2z_S;
  Double I_ERI_F3z_S_Fy2z_Px = I_ERI_F3z_S_Gxy2z_S+CDX*I_ERI_F3z_S_Fy2z_S;
  Double I_ERI_F3x_S_F3z_Px = I_ERI_F3x_S_Gx3z_S+CDX*I_ERI_F3x_S_F3z_S;
  Double I_ERI_F2xy_S_F3z_Px = I_ERI_F2xy_S_Gx3z_S+CDX*I_ERI_F2xy_S_F3z_S;
  Double I_ERI_F2xz_S_F3z_Px = I_ERI_F2xz_S_Gx3z_S+CDX*I_ERI_F2xz_S_F3z_S;
  Double I_ERI_Fx2y_S_F3z_Px = I_ERI_Fx2y_S_Gx3z_S+CDX*I_ERI_Fx2y_S_F3z_S;
  Double I_ERI_Fxyz_S_F3z_Px = I_ERI_Fxyz_S_Gx3z_S+CDX*I_ERI_Fxyz_S_F3z_S;
  Double I_ERI_Fx2z_S_F3z_Px = I_ERI_Fx2z_S_Gx3z_S+CDX*I_ERI_Fx2z_S_F3z_S;
  Double I_ERI_F3y_S_F3z_Px = I_ERI_F3y_S_Gx3z_S+CDX*I_ERI_F3y_S_F3z_S;
  Double I_ERI_F2yz_S_F3z_Px = I_ERI_F2yz_S_Gx3z_S+CDX*I_ERI_F2yz_S_F3z_S;
  Double I_ERI_Fy2z_S_F3z_Px = I_ERI_Fy2z_S_Gx3z_S+CDX*I_ERI_Fy2z_S_F3z_S;
  Double I_ERI_F3z_S_F3z_Px = I_ERI_F3z_S_Gx3z_S+CDX*I_ERI_F3z_S_F3z_S;
  Double I_ERI_F3x_S_F3x_Py = I_ERI_F3x_S_G3xy_S+CDY*I_ERI_F3x_S_F3x_S;
  Double I_ERI_F2xy_S_F3x_Py = I_ERI_F2xy_S_G3xy_S+CDY*I_ERI_F2xy_S_F3x_S;
  Double I_ERI_F2xz_S_F3x_Py = I_ERI_F2xz_S_G3xy_S+CDY*I_ERI_F2xz_S_F3x_S;
  Double I_ERI_Fx2y_S_F3x_Py = I_ERI_Fx2y_S_G3xy_S+CDY*I_ERI_Fx2y_S_F3x_S;
  Double I_ERI_Fxyz_S_F3x_Py = I_ERI_Fxyz_S_G3xy_S+CDY*I_ERI_Fxyz_S_F3x_S;
  Double I_ERI_Fx2z_S_F3x_Py = I_ERI_Fx2z_S_G3xy_S+CDY*I_ERI_Fx2z_S_F3x_S;
  Double I_ERI_F3y_S_F3x_Py = I_ERI_F3y_S_G3xy_S+CDY*I_ERI_F3y_S_F3x_S;
  Double I_ERI_F2yz_S_F3x_Py = I_ERI_F2yz_S_G3xy_S+CDY*I_ERI_F2yz_S_F3x_S;
  Double I_ERI_Fy2z_S_F3x_Py = I_ERI_Fy2z_S_G3xy_S+CDY*I_ERI_Fy2z_S_F3x_S;
  Double I_ERI_F3z_S_F3x_Py = I_ERI_F3z_S_G3xy_S+CDY*I_ERI_F3z_S_F3x_S;
  Double I_ERI_F3x_S_F2xy_Py = I_ERI_F3x_S_G2x2y_S+CDY*I_ERI_F3x_S_F2xy_S;
  Double I_ERI_F2xy_S_F2xy_Py = I_ERI_F2xy_S_G2x2y_S+CDY*I_ERI_F2xy_S_F2xy_S;
  Double I_ERI_F2xz_S_F2xy_Py = I_ERI_F2xz_S_G2x2y_S+CDY*I_ERI_F2xz_S_F2xy_S;
  Double I_ERI_Fx2y_S_F2xy_Py = I_ERI_Fx2y_S_G2x2y_S+CDY*I_ERI_Fx2y_S_F2xy_S;
  Double I_ERI_Fxyz_S_F2xy_Py = I_ERI_Fxyz_S_G2x2y_S+CDY*I_ERI_Fxyz_S_F2xy_S;
  Double I_ERI_Fx2z_S_F2xy_Py = I_ERI_Fx2z_S_G2x2y_S+CDY*I_ERI_Fx2z_S_F2xy_S;
  Double I_ERI_F3y_S_F2xy_Py = I_ERI_F3y_S_G2x2y_S+CDY*I_ERI_F3y_S_F2xy_S;
  Double I_ERI_F2yz_S_F2xy_Py = I_ERI_F2yz_S_G2x2y_S+CDY*I_ERI_F2yz_S_F2xy_S;
  Double I_ERI_Fy2z_S_F2xy_Py = I_ERI_Fy2z_S_G2x2y_S+CDY*I_ERI_Fy2z_S_F2xy_S;
  Double I_ERI_F3z_S_F2xy_Py = I_ERI_F3z_S_G2x2y_S+CDY*I_ERI_F3z_S_F2xy_S;
  Double I_ERI_F3x_S_F2xz_Py = I_ERI_F3x_S_G2xyz_S+CDY*I_ERI_F3x_S_F2xz_S;
  Double I_ERI_F2xy_S_F2xz_Py = I_ERI_F2xy_S_G2xyz_S+CDY*I_ERI_F2xy_S_F2xz_S;
  Double I_ERI_F2xz_S_F2xz_Py = I_ERI_F2xz_S_G2xyz_S+CDY*I_ERI_F2xz_S_F2xz_S;
  Double I_ERI_Fx2y_S_F2xz_Py = I_ERI_Fx2y_S_G2xyz_S+CDY*I_ERI_Fx2y_S_F2xz_S;
  Double I_ERI_Fxyz_S_F2xz_Py = I_ERI_Fxyz_S_G2xyz_S+CDY*I_ERI_Fxyz_S_F2xz_S;
  Double I_ERI_Fx2z_S_F2xz_Py = I_ERI_Fx2z_S_G2xyz_S+CDY*I_ERI_Fx2z_S_F2xz_S;
  Double I_ERI_F3y_S_F2xz_Py = I_ERI_F3y_S_G2xyz_S+CDY*I_ERI_F3y_S_F2xz_S;
  Double I_ERI_F2yz_S_F2xz_Py = I_ERI_F2yz_S_G2xyz_S+CDY*I_ERI_F2yz_S_F2xz_S;
  Double I_ERI_Fy2z_S_F2xz_Py = I_ERI_Fy2z_S_G2xyz_S+CDY*I_ERI_Fy2z_S_F2xz_S;
  Double I_ERI_F3z_S_F2xz_Py = I_ERI_F3z_S_G2xyz_S+CDY*I_ERI_F3z_S_F2xz_S;
  Double I_ERI_F3x_S_Fx2y_Py = I_ERI_F3x_S_Gx3y_S+CDY*I_ERI_F3x_S_Fx2y_S;
  Double I_ERI_F2xy_S_Fx2y_Py = I_ERI_F2xy_S_Gx3y_S+CDY*I_ERI_F2xy_S_Fx2y_S;
  Double I_ERI_F2xz_S_Fx2y_Py = I_ERI_F2xz_S_Gx3y_S+CDY*I_ERI_F2xz_S_Fx2y_S;
  Double I_ERI_Fx2y_S_Fx2y_Py = I_ERI_Fx2y_S_Gx3y_S+CDY*I_ERI_Fx2y_S_Fx2y_S;
  Double I_ERI_Fxyz_S_Fx2y_Py = I_ERI_Fxyz_S_Gx3y_S+CDY*I_ERI_Fxyz_S_Fx2y_S;
  Double I_ERI_Fx2z_S_Fx2y_Py = I_ERI_Fx2z_S_Gx3y_S+CDY*I_ERI_Fx2z_S_Fx2y_S;
  Double I_ERI_F3y_S_Fx2y_Py = I_ERI_F3y_S_Gx3y_S+CDY*I_ERI_F3y_S_Fx2y_S;
  Double I_ERI_F2yz_S_Fx2y_Py = I_ERI_F2yz_S_Gx3y_S+CDY*I_ERI_F2yz_S_Fx2y_S;
  Double I_ERI_Fy2z_S_Fx2y_Py = I_ERI_Fy2z_S_Gx3y_S+CDY*I_ERI_Fy2z_S_Fx2y_S;
  Double I_ERI_F3z_S_Fx2y_Py = I_ERI_F3z_S_Gx3y_S+CDY*I_ERI_F3z_S_Fx2y_S;
  Double I_ERI_F3x_S_Fxyz_Py = I_ERI_F3x_S_Gx2yz_S+CDY*I_ERI_F3x_S_Fxyz_S;
  Double I_ERI_F2xy_S_Fxyz_Py = I_ERI_F2xy_S_Gx2yz_S+CDY*I_ERI_F2xy_S_Fxyz_S;
  Double I_ERI_F2xz_S_Fxyz_Py = I_ERI_F2xz_S_Gx2yz_S+CDY*I_ERI_F2xz_S_Fxyz_S;
  Double I_ERI_Fx2y_S_Fxyz_Py = I_ERI_Fx2y_S_Gx2yz_S+CDY*I_ERI_Fx2y_S_Fxyz_S;
  Double I_ERI_Fxyz_S_Fxyz_Py = I_ERI_Fxyz_S_Gx2yz_S+CDY*I_ERI_Fxyz_S_Fxyz_S;
  Double I_ERI_Fx2z_S_Fxyz_Py = I_ERI_Fx2z_S_Gx2yz_S+CDY*I_ERI_Fx2z_S_Fxyz_S;
  Double I_ERI_F3y_S_Fxyz_Py = I_ERI_F3y_S_Gx2yz_S+CDY*I_ERI_F3y_S_Fxyz_S;
  Double I_ERI_F2yz_S_Fxyz_Py = I_ERI_F2yz_S_Gx2yz_S+CDY*I_ERI_F2yz_S_Fxyz_S;
  Double I_ERI_Fy2z_S_Fxyz_Py = I_ERI_Fy2z_S_Gx2yz_S+CDY*I_ERI_Fy2z_S_Fxyz_S;
  Double I_ERI_F3z_S_Fxyz_Py = I_ERI_F3z_S_Gx2yz_S+CDY*I_ERI_F3z_S_Fxyz_S;
  Double I_ERI_F3x_S_Fx2z_Py = I_ERI_F3x_S_Gxy2z_S+CDY*I_ERI_F3x_S_Fx2z_S;
  Double I_ERI_F2xy_S_Fx2z_Py = I_ERI_F2xy_S_Gxy2z_S+CDY*I_ERI_F2xy_S_Fx2z_S;
  Double I_ERI_F2xz_S_Fx2z_Py = I_ERI_F2xz_S_Gxy2z_S+CDY*I_ERI_F2xz_S_Fx2z_S;
  Double I_ERI_Fx2y_S_Fx2z_Py = I_ERI_Fx2y_S_Gxy2z_S+CDY*I_ERI_Fx2y_S_Fx2z_S;
  Double I_ERI_Fxyz_S_Fx2z_Py = I_ERI_Fxyz_S_Gxy2z_S+CDY*I_ERI_Fxyz_S_Fx2z_S;
  Double I_ERI_Fx2z_S_Fx2z_Py = I_ERI_Fx2z_S_Gxy2z_S+CDY*I_ERI_Fx2z_S_Fx2z_S;
  Double I_ERI_F3y_S_Fx2z_Py = I_ERI_F3y_S_Gxy2z_S+CDY*I_ERI_F3y_S_Fx2z_S;
  Double I_ERI_F2yz_S_Fx2z_Py = I_ERI_F2yz_S_Gxy2z_S+CDY*I_ERI_F2yz_S_Fx2z_S;
  Double I_ERI_Fy2z_S_Fx2z_Py = I_ERI_Fy2z_S_Gxy2z_S+CDY*I_ERI_Fy2z_S_Fx2z_S;
  Double I_ERI_F3z_S_Fx2z_Py = I_ERI_F3z_S_Gxy2z_S+CDY*I_ERI_F3z_S_Fx2z_S;
  Double I_ERI_F3x_S_F3y_Py = I_ERI_F3x_S_G4y_S+CDY*I_ERI_F3x_S_F3y_S;
  Double I_ERI_F2xy_S_F3y_Py = I_ERI_F2xy_S_G4y_S+CDY*I_ERI_F2xy_S_F3y_S;
  Double I_ERI_F2xz_S_F3y_Py = I_ERI_F2xz_S_G4y_S+CDY*I_ERI_F2xz_S_F3y_S;
  Double I_ERI_Fx2y_S_F3y_Py = I_ERI_Fx2y_S_G4y_S+CDY*I_ERI_Fx2y_S_F3y_S;
  Double I_ERI_Fxyz_S_F3y_Py = I_ERI_Fxyz_S_G4y_S+CDY*I_ERI_Fxyz_S_F3y_S;
  Double I_ERI_Fx2z_S_F3y_Py = I_ERI_Fx2z_S_G4y_S+CDY*I_ERI_Fx2z_S_F3y_S;
  Double I_ERI_F3y_S_F3y_Py = I_ERI_F3y_S_G4y_S+CDY*I_ERI_F3y_S_F3y_S;
  Double I_ERI_F2yz_S_F3y_Py = I_ERI_F2yz_S_G4y_S+CDY*I_ERI_F2yz_S_F3y_S;
  Double I_ERI_Fy2z_S_F3y_Py = I_ERI_Fy2z_S_G4y_S+CDY*I_ERI_Fy2z_S_F3y_S;
  Double I_ERI_F3z_S_F3y_Py = I_ERI_F3z_S_G4y_S+CDY*I_ERI_F3z_S_F3y_S;
  Double I_ERI_F3x_S_F2yz_Py = I_ERI_F3x_S_G3yz_S+CDY*I_ERI_F3x_S_F2yz_S;
  Double I_ERI_F2xy_S_F2yz_Py = I_ERI_F2xy_S_G3yz_S+CDY*I_ERI_F2xy_S_F2yz_S;
  Double I_ERI_F2xz_S_F2yz_Py = I_ERI_F2xz_S_G3yz_S+CDY*I_ERI_F2xz_S_F2yz_S;
  Double I_ERI_Fx2y_S_F2yz_Py = I_ERI_Fx2y_S_G3yz_S+CDY*I_ERI_Fx2y_S_F2yz_S;
  Double I_ERI_Fxyz_S_F2yz_Py = I_ERI_Fxyz_S_G3yz_S+CDY*I_ERI_Fxyz_S_F2yz_S;
  Double I_ERI_Fx2z_S_F2yz_Py = I_ERI_Fx2z_S_G3yz_S+CDY*I_ERI_Fx2z_S_F2yz_S;
  Double I_ERI_F3y_S_F2yz_Py = I_ERI_F3y_S_G3yz_S+CDY*I_ERI_F3y_S_F2yz_S;
  Double I_ERI_F2yz_S_F2yz_Py = I_ERI_F2yz_S_G3yz_S+CDY*I_ERI_F2yz_S_F2yz_S;
  Double I_ERI_Fy2z_S_F2yz_Py = I_ERI_Fy2z_S_G3yz_S+CDY*I_ERI_Fy2z_S_F2yz_S;
  Double I_ERI_F3z_S_F2yz_Py = I_ERI_F3z_S_G3yz_S+CDY*I_ERI_F3z_S_F2yz_S;
  Double I_ERI_F3x_S_Fy2z_Py = I_ERI_F3x_S_G2y2z_S+CDY*I_ERI_F3x_S_Fy2z_S;
  Double I_ERI_F2xy_S_Fy2z_Py = I_ERI_F2xy_S_G2y2z_S+CDY*I_ERI_F2xy_S_Fy2z_S;
  Double I_ERI_F2xz_S_Fy2z_Py = I_ERI_F2xz_S_G2y2z_S+CDY*I_ERI_F2xz_S_Fy2z_S;
  Double I_ERI_Fx2y_S_Fy2z_Py = I_ERI_Fx2y_S_G2y2z_S+CDY*I_ERI_Fx2y_S_Fy2z_S;
  Double I_ERI_Fxyz_S_Fy2z_Py = I_ERI_Fxyz_S_G2y2z_S+CDY*I_ERI_Fxyz_S_Fy2z_S;
  Double I_ERI_Fx2z_S_Fy2z_Py = I_ERI_Fx2z_S_G2y2z_S+CDY*I_ERI_Fx2z_S_Fy2z_S;
  Double I_ERI_F3y_S_Fy2z_Py = I_ERI_F3y_S_G2y2z_S+CDY*I_ERI_F3y_S_Fy2z_S;
  Double I_ERI_F2yz_S_Fy2z_Py = I_ERI_F2yz_S_G2y2z_S+CDY*I_ERI_F2yz_S_Fy2z_S;
  Double I_ERI_Fy2z_S_Fy2z_Py = I_ERI_Fy2z_S_G2y2z_S+CDY*I_ERI_Fy2z_S_Fy2z_S;
  Double I_ERI_F3z_S_Fy2z_Py = I_ERI_F3z_S_G2y2z_S+CDY*I_ERI_F3z_S_Fy2z_S;
  Double I_ERI_F3x_S_F3z_Py = I_ERI_F3x_S_Gy3z_S+CDY*I_ERI_F3x_S_F3z_S;
  Double I_ERI_F2xy_S_F3z_Py = I_ERI_F2xy_S_Gy3z_S+CDY*I_ERI_F2xy_S_F3z_S;
  Double I_ERI_F2xz_S_F3z_Py = I_ERI_F2xz_S_Gy3z_S+CDY*I_ERI_F2xz_S_F3z_S;
  Double I_ERI_Fx2y_S_F3z_Py = I_ERI_Fx2y_S_Gy3z_S+CDY*I_ERI_Fx2y_S_F3z_S;
  Double I_ERI_Fxyz_S_F3z_Py = I_ERI_Fxyz_S_Gy3z_S+CDY*I_ERI_Fxyz_S_F3z_S;
  Double I_ERI_Fx2z_S_F3z_Py = I_ERI_Fx2z_S_Gy3z_S+CDY*I_ERI_Fx2z_S_F3z_S;
  Double I_ERI_F3y_S_F3z_Py = I_ERI_F3y_S_Gy3z_S+CDY*I_ERI_F3y_S_F3z_S;
  Double I_ERI_F2yz_S_F3z_Py = I_ERI_F2yz_S_Gy3z_S+CDY*I_ERI_F2yz_S_F3z_S;
  Double I_ERI_Fy2z_S_F3z_Py = I_ERI_Fy2z_S_Gy3z_S+CDY*I_ERI_Fy2z_S_F3z_S;
  Double I_ERI_F3z_S_F3z_Py = I_ERI_F3z_S_Gy3z_S+CDY*I_ERI_F3z_S_F3z_S;
  Double I_ERI_F3x_S_F3x_Pz = I_ERI_F3x_S_G3xz_S+CDZ*I_ERI_F3x_S_F3x_S;
  Double I_ERI_F2xy_S_F3x_Pz = I_ERI_F2xy_S_G3xz_S+CDZ*I_ERI_F2xy_S_F3x_S;
  Double I_ERI_F2xz_S_F3x_Pz = I_ERI_F2xz_S_G3xz_S+CDZ*I_ERI_F2xz_S_F3x_S;
  Double I_ERI_Fx2y_S_F3x_Pz = I_ERI_Fx2y_S_G3xz_S+CDZ*I_ERI_Fx2y_S_F3x_S;
  Double I_ERI_Fxyz_S_F3x_Pz = I_ERI_Fxyz_S_G3xz_S+CDZ*I_ERI_Fxyz_S_F3x_S;
  Double I_ERI_Fx2z_S_F3x_Pz = I_ERI_Fx2z_S_G3xz_S+CDZ*I_ERI_Fx2z_S_F3x_S;
  Double I_ERI_F3y_S_F3x_Pz = I_ERI_F3y_S_G3xz_S+CDZ*I_ERI_F3y_S_F3x_S;
  Double I_ERI_F2yz_S_F3x_Pz = I_ERI_F2yz_S_G3xz_S+CDZ*I_ERI_F2yz_S_F3x_S;
  Double I_ERI_Fy2z_S_F3x_Pz = I_ERI_Fy2z_S_G3xz_S+CDZ*I_ERI_Fy2z_S_F3x_S;
  Double I_ERI_F3z_S_F3x_Pz = I_ERI_F3z_S_G3xz_S+CDZ*I_ERI_F3z_S_F3x_S;
  Double I_ERI_F3x_S_F2xy_Pz = I_ERI_F3x_S_G2xyz_S+CDZ*I_ERI_F3x_S_F2xy_S;
  Double I_ERI_F2xy_S_F2xy_Pz = I_ERI_F2xy_S_G2xyz_S+CDZ*I_ERI_F2xy_S_F2xy_S;
  Double I_ERI_F2xz_S_F2xy_Pz = I_ERI_F2xz_S_G2xyz_S+CDZ*I_ERI_F2xz_S_F2xy_S;
  Double I_ERI_Fx2y_S_F2xy_Pz = I_ERI_Fx2y_S_G2xyz_S+CDZ*I_ERI_Fx2y_S_F2xy_S;
  Double I_ERI_Fxyz_S_F2xy_Pz = I_ERI_Fxyz_S_G2xyz_S+CDZ*I_ERI_Fxyz_S_F2xy_S;
  Double I_ERI_Fx2z_S_F2xy_Pz = I_ERI_Fx2z_S_G2xyz_S+CDZ*I_ERI_Fx2z_S_F2xy_S;
  Double I_ERI_F3y_S_F2xy_Pz = I_ERI_F3y_S_G2xyz_S+CDZ*I_ERI_F3y_S_F2xy_S;
  Double I_ERI_F2yz_S_F2xy_Pz = I_ERI_F2yz_S_G2xyz_S+CDZ*I_ERI_F2yz_S_F2xy_S;
  Double I_ERI_Fy2z_S_F2xy_Pz = I_ERI_Fy2z_S_G2xyz_S+CDZ*I_ERI_Fy2z_S_F2xy_S;
  Double I_ERI_F3z_S_F2xy_Pz = I_ERI_F3z_S_G2xyz_S+CDZ*I_ERI_F3z_S_F2xy_S;
  Double I_ERI_F3x_S_F2xz_Pz = I_ERI_F3x_S_G2x2z_S+CDZ*I_ERI_F3x_S_F2xz_S;
  Double I_ERI_F2xy_S_F2xz_Pz = I_ERI_F2xy_S_G2x2z_S+CDZ*I_ERI_F2xy_S_F2xz_S;
  Double I_ERI_F2xz_S_F2xz_Pz = I_ERI_F2xz_S_G2x2z_S+CDZ*I_ERI_F2xz_S_F2xz_S;
  Double I_ERI_Fx2y_S_F2xz_Pz = I_ERI_Fx2y_S_G2x2z_S+CDZ*I_ERI_Fx2y_S_F2xz_S;
  Double I_ERI_Fxyz_S_F2xz_Pz = I_ERI_Fxyz_S_G2x2z_S+CDZ*I_ERI_Fxyz_S_F2xz_S;
  Double I_ERI_Fx2z_S_F2xz_Pz = I_ERI_Fx2z_S_G2x2z_S+CDZ*I_ERI_Fx2z_S_F2xz_S;
  Double I_ERI_F3y_S_F2xz_Pz = I_ERI_F3y_S_G2x2z_S+CDZ*I_ERI_F3y_S_F2xz_S;
  Double I_ERI_F2yz_S_F2xz_Pz = I_ERI_F2yz_S_G2x2z_S+CDZ*I_ERI_F2yz_S_F2xz_S;
  Double I_ERI_Fy2z_S_F2xz_Pz = I_ERI_Fy2z_S_G2x2z_S+CDZ*I_ERI_Fy2z_S_F2xz_S;
  Double I_ERI_F3z_S_F2xz_Pz = I_ERI_F3z_S_G2x2z_S+CDZ*I_ERI_F3z_S_F2xz_S;
  Double I_ERI_F3x_S_Fx2y_Pz = I_ERI_F3x_S_Gx2yz_S+CDZ*I_ERI_F3x_S_Fx2y_S;
  Double I_ERI_F2xy_S_Fx2y_Pz = I_ERI_F2xy_S_Gx2yz_S+CDZ*I_ERI_F2xy_S_Fx2y_S;
  Double I_ERI_F2xz_S_Fx2y_Pz = I_ERI_F2xz_S_Gx2yz_S+CDZ*I_ERI_F2xz_S_Fx2y_S;
  Double I_ERI_Fx2y_S_Fx2y_Pz = I_ERI_Fx2y_S_Gx2yz_S+CDZ*I_ERI_Fx2y_S_Fx2y_S;
  Double I_ERI_Fxyz_S_Fx2y_Pz = I_ERI_Fxyz_S_Gx2yz_S+CDZ*I_ERI_Fxyz_S_Fx2y_S;
  Double I_ERI_Fx2z_S_Fx2y_Pz = I_ERI_Fx2z_S_Gx2yz_S+CDZ*I_ERI_Fx2z_S_Fx2y_S;
  Double I_ERI_F3y_S_Fx2y_Pz = I_ERI_F3y_S_Gx2yz_S+CDZ*I_ERI_F3y_S_Fx2y_S;
  Double I_ERI_F2yz_S_Fx2y_Pz = I_ERI_F2yz_S_Gx2yz_S+CDZ*I_ERI_F2yz_S_Fx2y_S;
  Double I_ERI_Fy2z_S_Fx2y_Pz = I_ERI_Fy2z_S_Gx2yz_S+CDZ*I_ERI_Fy2z_S_Fx2y_S;
  Double I_ERI_F3z_S_Fx2y_Pz = I_ERI_F3z_S_Gx2yz_S+CDZ*I_ERI_F3z_S_Fx2y_S;
  Double I_ERI_F3x_S_Fxyz_Pz = I_ERI_F3x_S_Gxy2z_S+CDZ*I_ERI_F3x_S_Fxyz_S;
  Double I_ERI_F2xy_S_Fxyz_Pz = I_ERI_F2xy_S_Gxy2z_S+CDZ*I_ERI_F2xy_S_Fxyz_S;
  Double I_ERI_F2xz_S_Fxyz_Pz = I_ERI_F2xz_S_Gxy2z_S+CDZ*I_ERI_F2xz_S_Fxyz_S;
  Double I_ERI_Fx2y_S_Fxyz_Pz = I_ERI_Fx2y_S_Gxy2z_S+CDZ*I_ERI_Fx2y_S_Fxyz_S;
  Double I_ERI_Fxyz_S_Fxyz_Pz = I_ERI_Fxyz_S_Gxy2z_S+CDZ*I_ERI_Fxyz_S_Fxyz_S;
  Double I_ERI_Fx2z_S_Fxyz_Pz = I_ERI_Fx2z_S_Gxy2z_S+CDZ*I_ERI_Fx2z_S_Fxyz_S;
  Double I_ERI_F3y_S_Fxyz_Pz = I_ERI_F3y_S_Gxy2z_S+CDZ*I_ERI_F3y_S_Fxyz_S;
  Double I_ERI_F2yz_S_Fxyz_Pz = I_ERI_F2yz_S_Gxy2z_S+CDZ*I_ERI_F2yz_S_Fxyz_S;
  Double I_ERI_Fy2z_S_Fxyz_Pz = I_ERI_Fy2z_S_Gxy2z_S+CDZ*I_ERI_Fy2z_S_Fxyz_S;
  Double I_ERI_F3z_S_Fxyz_Pz = I_ERI_F3z_S_Gxy2z_S+CDZ*I_ERI_F3z_S_Fxyz_S;
  Double I_ERI_F3x_S_Fx2z_Pz = I_ERI_F3x_S_Gx3z_S+CDZ*I_ERI_F3x_S_Fx2z_S;
  Double I_ERI_F2xy_S_Fx2z_Pz = I_ERI_F2xy_S_Gx3z_S+CDZ*I_ERI_F2xy_S_Fx2z_S;
  Double I_ERI_F2xz_S_Fx2z_Pz = I_ERI_F2xz_S_Gx3z_S+CDZ*I_ERI_F2xz_S_Fx2z_S;
  Double I_ERI_Fx2y_S_Fx2z_Pz = I_ERI_Fx2y_S_Gx3z_S+CDZ*I_ERI_Fx2y_S_Fx2z_S;
  Double I_ERI_Fxyz_S_Fx2z_Pz = I_ERI_Fxyz_S_Gx3z_S+CDZ*I_ERI_Fxyz_S_Fx2z_S;
  Double I_ERI_Fx2z_S_Fx2z_Pz = I_ERI_Fx2z_S_Gx3z_S+CDZ*I_ERI_Fx2z_S_Fx2z_S;
  Double I_ERI_F3y_S_Fx2z_Pz = I_ERI_F3y_S_Gx3z_S+CDZ*I_ERI_F3y_S_Fx2z_S;
  Double I_ERI_F2yz_S_Fx2z_Pz = I_ERI_F2yz_S_Gx3z_S+CDZ*I_ERI_F2yz_S_Fx2z_S;
  Double I_ERI_Fy2z_S_Fx2z_Pz = I_ERI_Fy2z_S_Gx3z_S+CDZ*I_ERI_Fy2z_S_Fx2z_S;
  Double I_ERI_F3z_S_Fx2z_Pz = I_ERI_F3z_S_Gx3z_S+CDZ*I_ERI_F3z_S_Fx2z_S;
  Double I_ERI_F3x_S_F3y_Pz = I_ERI_F3x_S_G3yz_S+CDZ*I_ERI_F3x_S_F3y_S;
  Double I_ERI_F2xy_S_F3y_Pz = I_ERI_F2xy_S_G3yz_S+CDZ*I_ERI_F2xy_S_F3y_S;
  Double I_ERI_F2xz_S_F3y_Pz = I_ERI_F2xz_S_G3yz_S+CDZ*I_ERI_F2xz_S_F3y_S;
  Double I_ERI_Fx2y_S_F3y_Pz = I_ERI_Fx2y_S_G3yz_S+CDZ*I_ERI_Fx2y_S_F3y_S;
  Double I_ERI_Fxyz_S_F3y_Pz = I_ERI_Fxyz_S_G3yz_S+CDZ*I_ERI_Fxyz_S_F3y_S;
  Double I_ERI_Fx2z_S_F3y_Pz = I_ERI_Fx2z_S_G3yz_S+CDZ*I_ERI_Fx2z_S_F3y_S;
  Double I_ERI_F3y_S_F3y_Pz = I_ERI_F3y_S_G3yz_S+CDZ*I_ERI_F3y_S_F3y_S;
  Double I_ERI_F2yz_S_F3y_Pz = I_ERI_F2yz_S_G3yz_S+CDZ*I_ERI_F2yz_S_F3y_S;
  Double I_ERI_Fy2z_S_F3y_Pz = I_ERI_Fy2z_S_G3yz_S+CDZ*I_ERI_Fy2z_S_F3y_S;
  Double I_ERI_F3z_S_F3y_Pz = I_ERI_F3z_S_G3yz_S+CDZ*I_ERI_F3z_S_F3y_S;
  Double I_ERI_F3x_S_F2yz_Pz = I_ERI_F3x_S_G2y2z_S+CDZ*I_ERI_F3x_S_F2yz_S;
  Double I_ERI_F2xy_S_F2yz_Pz = I_ERI_F2xy_S_G2y2z_S+CDZ*I_ERI_F2xy_S_F2yz_S;
  Double I_ERI_F2xz_S_F2yz_Pz = I_ERI_F2xz_S_G2y2z_S+CDZ*I_ERI_F2xz_S_F2yz_S;
  Double I_ERI_Fx2y_S_F2yz_Pz = I_ERI_Fx2y_S_G2y2z_S+CDZ*I_ERI_Fx2y_S_F2yz_S;
  Double I_ERI_Fxyz_S_F2yz_Pz = I_ERI_Fxyz_S_G2y2z_S+CDZ*I_ERI_Fxyz_S_F2yz_S;
  Double I_ERI_Fx2z_S_F2yz_Pz = I_ERI_Fx2z_S_G2y2z_S+CDZ*I_ERI_Fx2z_S_F2yz_S;
  Double I_ERI_F3y_S_F2yz_Pz = I_ERI_F3y_S_G2y2z_S+CDZ*I_ERI_F3y_S_F2yz_S;
  Double I_ERI_F2yz_S_F2yz_Pz = I_ERI_F2yz_S_G2y2z_S+CDZ*I_ERI_F2yz_S_F2yz_S;
  Double I_ERI_Fy2z_S_F2yz_Pz = I_ERI_Fy2z_S_G2y2z_S+CDZ*I_ERI_Fy2z_S_F2yz_S;
  Double I_ERI_F3z_S_F2yz_Pz = I_ERI_F3z_S_G2y2z_S+CDZ*I_ERI_F3z_S_F2yz_S;
  Double I_ERI_F3x_S_Fy2z_Pz = I_ERI_F3x_S_Gy3z_S+CDZ*I_ERI_F3x_S_Fy2z_S;
  Double I_ERI_F2xy_S_Fy2z_Pz = I_ERI_F2xy_S_Gy3z_S+CDZ*I_ERI_F2xy_S_Fy2z_S;
  Double I_ERI_F2xz_S_Fy2z_Pz = I_ERI_F2xz_S_Gy3z_S+CDZ*I_ERI_F2xz_S_Fy2z_S;
  Double I_ERI_Fx2y_S_Fy2z_Pz = I_ERI_Fx2y_S_Gy3z_S+CDZ*I_ERI_Fx2y_S_Fy2z_S;
  Double I_ERI_Fxyz_S_Fy2z_Pz = I_ERI_Fxyz_S_Gy3z_S+CDZ*I_ERI_Fxyz_S_Fy2z_S;
  Double I_ERI_Fx2z_S_Fy2z_Pz = I_ERI_Fx2z_S_Gy3z_S+CDZ*I_ERI_Fx2z_S_Fy2z_S;
  Double I_ERI_F3y_S_Fy2z_Pz = I_ERI_F3y_S_Gy3z_S+CDZ*I_ERI_F3y_S_Fy2z_S;
  Double I_ERI_F2yz_S_Fy2z_Pz = I_ERI_F2yz_S_Gy3z_S+CDZ*I_ERI_F2yz_S_Fy2z_S;
  Double I_ERI_Fy2z_S_Fy2z_Pz = I_ERI_Fy2z_S_Gy3z_S+CDZ*I_ERI_Fy2z_S_Fy2z_S;
  Double I_ERI_F3z_S_Fy2z_Pz = I_ERI_F3z_S_Gy3z_S+CDZ*I_ERI_F3z_S_Fy2z_S;
  Double I_ERI_F3x_S_F3z_Pz = I_ERI_F3x_S_G4z_S+CDZ*I_ERI_F3x_S_F3z_S;
  Double I_ERI_F2xy_S_F3z_Pz = I_ERI_F2xy_S_G4z_S+CDZ*I_ERI_F2xy_S_F3z_S;
  Double I_ERI_F2xz_S_F3z_Pz = I_ERI_F2xz_S_G4z_S+CDZ*I_ERI_F2xz_S_F3z_S;
  Double I_ERI_Fx2y_S_F3z_Pz = I_ERI_Fx2y_S_G4z_S+CDZ*I_ERI_Fx2y_S_F3z_S;
  Double I_ERI_Fxyz_S_F3z_Pz = I_ERI_Fxyz_S_G4z_S+CDZ*I_ERI_Fxyz_S_F3z_S;
  Double I_ERI_Fx2z_S_F3z_Pz = I_ERI_Fx2z_S_G4z_S+CDZ*I_ERI_Fx2z_S_F3z_S;
  Double I_ERI_F3y_S_F3z_Pz = I_ERI_F3y_S_G4z_S+CDZ*I_ERI_F3y_S_F3z_S;
  Double I_ERI_F2yz_S_F3z_Pz = I_ERI_F2yz_S_G4z_S+CDZ*I_ERI_F2yz_S_F3z_S;
  Double I_ERI_Fy2z_S_F3z_Pz = I_ERI_Fy2z_S_G4z_S+CDZ*I_ERI_Fy2z_S_F3z_S;
  Double I_ERI_F3z_S_F3z_Pz = I_ERI_F3z_S_G4z_S+CDZ*I_ERI_F3z_S_F3z_S;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_F_P
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_G_S
   * RHS shell quartet name: SQ_ERI_G_S_F_S
   ************************************************************/
  Double I_ERI_G4x_S_F3x_Px = I_ERI_G4x_S_G4x_S+CDX*I_ERI_G4x_S_F3x_S;
  Double I_ERI_G3xy_S_F3x_Px = I_ERI_G3xy_S_G4x_S+CDX*I_ERI_G3xy_S_F3x_S;
  Double I_ERI_G3xz_S_F3x_Px = I_ERI_G3xz_S_G4x_S+CDX*I_ERI_G3xz_S_F3x_S;
  Double I_ERI_G2x2y_S_F3x_Px = I_ERI_G2x2y_S_G4x_S+CDX*I_ERI_G2x2y_S_F3x_S;
  Double I_ERI_G2xyz_S_F3x_Px = I_ERI_G2xyz_S_G4x_S+CDX*I_ERI_G2xyz_S_F3x_S;
  Double I_ERI_G2x2z_S_F3x_Px = I_ERI_G2x2z_S_G4x_S+CDX*I_ERI_G2x2z_S_F3x_S;
  Double I_ERI_Gx3y_S_F3x_Px = I_ERI_Gx3y_S_G4x_S+CDX*I_ERI_Gx3y_S_F3x_S;
  Double I_ERI_Gx2yz_S_F3x_Px = I_ERI_Gx2yz_S_G4x_S+CDX*I_ERI_Gx2yz_S_F3x_S;
  Double I_ERI_Gxy2z_S_F3x_Px = I_ERI_Gxy2z_S_G4x_S+CDX*I_ERI_Gxy2z_S_F3x_S;
  Double I_ERI_Gx3z_S_F3x_Px = I_ERI_Gx3z_S_G4x_S+CDX*I_ERI_Gx3z_S_F3x_S;
  Double I_ERI_G4y_S_F3x_Px = I_ERI_G4y_S_G4x_S+CDX*I_ERI_G4y_S_F3x_S;
  Double I_ERI_G3yz_S_F3x_Px = I_ERI_G3yz_S_G4x_S+CDX*I_ERI_G3yz_S_F3x_S;
  Double I_ERI_G2y2z_S_F3x_Px = I_ERI_G2y2z_S_G4x_S+CDX*I_ERI_G2y2z_S_F3x_S;
  Double I_ERI_Gy3z_S_F3x_Px = I_ERI_Gy3z_S_G4x_S+CDX*I_ERI_Gy3z_S_F3x_S;
  Double I_ERI_G4z_S_F3x_Px = I_ERI_G4z_S_G4x_S+CDX*I_ERI_G4z_S_F3x_S;
  Double I_ERI_G4x_S_F2xy_Px = I_ERI_G4x_S_G3xy_S+CDX*I_ERI_G4x_S_F2xy_S;
  Double I_ERI_G3xy_S_F2xy_Px = I_ERI_G3xy_S_G3xy_S+CDX*I_ERI_G3xy_S_F2xy_S;
  Double I_ERI_G3xz_S_F2xy_Px = I_ERI_G3xz_S_G3xy_S+CDX*I_ERI_G3xz_S_F2xy_S;
  Double I_ERI_G2x2y_S_F2xy_Px = I_ERI_G2x2y_S_G3xy_S+CDX*I_ERI_G2x2y_S_F2xy_S;
  Double I_ERI_G2xyz_S_F2xy_Px = I_ERI_G2xyz_S_G3xy_S+CDX*I_ERI_G2xyz_S_F2xy_S;
  Double I_ERI_G2x2z_S_F2xy_Px = I_ERI_G2x2z_S_G3xy_S+CDX*I_ERI_G2x2z_S_F2xy_S;
  Double I_ERI_Gx3y_S_F2xy_Px = I_ERI_Gx3y_S_G3xy_S+CDX*I_ERI_Gx3y_S_F2xy_S;
  Double I_ERI_Gx2yz_S_F2xy_Px = I_ERI_Gx2yz_S_G3xy_S+CDX*I_ERI_Gx2yz_S_F2xy_S;
  Double I_ERI_Gxy2z_S_F2xy_Px = I_ERI_Gxy2z_S_G3xy_S+CDX*I_ERI_Gxy2z_S_F2xy_S;
  Double I_ERI_Gx3z_S_F2xy_Px = I_ERI_Gx3z_S_G3xy_S+CDX*I_ERI_Gx3z_S_F2xy_S;
  Double I_ERI_G4y_S_F2xy_Px = I_ERI_G4y_S_G3xy_S+CDX*I_ERI_G4y_S_F2xy_S;
  Double I_ERI_G3yz_S_F2xy_Px = I_ERI_G3yz_S_G3xy_S+CDX*I_ERI_G3yz_S_F2xy_S;
  Double I_ERI_G2y2z_S_F2xy_Px = I_ERI_G2y2z_S_G3xy_S+CDX*I_ERI_G2y2z_S_F2xy_S;
  Double I_ERI_Gy3z_S_F2xy_Px = I_ERI_Gy3z_S_G3xy_S+CDX*I_ERI_Gy3z_S_F2xy_S;
  Double I_ERI_G4z_S_F2xy_Px = I_ERI_G4z_S_G3xy_S+CDX*I_ERI_G4z_S_F2xy_S;
  Double I_ERI_G4x_S_F2xz_Px = I_ERI_G4x_S_G3xz_S+CDX*I_ERI_G4x_S_F2xz_S;
  Double I_ERI_G3xy_S_F2xz_Px = I_ERI_G3xy_S_G3xz_S+CDX*I_ERI_G3xy_S_F2xz_S;
  Double I_ERI_G3xz_S_F2xz_Px = I_ERI_G3xz_S_G3xz_S+CDX*I_ERI_G3xz_S_F2xz_S;
  Double I_ERI_G2x2y_S_F2xz_Px = I_ERI_G2x2y_S_G3xz_S+CDX*I_ERI_G2x2y_S_F2xz_S;
  Double I_ERI_G2xyz_S_F2xz_Px = I_ERI_G2xyz_S_G3xz_S+CDX*I_ERI_G2xyz_S_F2xz_S;
  Double I_ERI_G2x2z_S_F2xz_Px = I_ERI_G2x2z_S_G3xz_S+CDX*I_ERI_G2x2z_S_F2xz_S;
  Double I_ERI_Gx3y_S_F2xz_Px = I_ERI_Gx3y_S_G3xz_S+CDX*I_ERI_Gx3y_S_F2xz_S;
  Double I_ERI_Gx2yz_S_F2xz_Px = I_ERI_Gx2yz_S_G3xz_S+CDX*I_ERI_Gx2yz_S_F2xz_S;
  Double I_ERI_Gxy2z_S_F2xz_Px = I_ERI_Gxy2z_S_G3xz_S+CDX*I_ERI_Gxy2z_S_F2xz_S;
  Double I_ERI_Gx3z_S_F2xz_Px = I_ERI_Gx3z_S_G3xz_S+CDX*I_ERI_Gx3z_S_F2xz_S;
  Double I_ERI_G4y_S_F2xz_Px = I_ERI_G4y_S_G3xz_S+CDX*I_ERI_G4y_S_F2xz_S;
  Double I_ERI_G3yz_S_F2xz_Px = I_ERI_G3yz_S_G3xz_S+CDX*I_ERI_G3yz_S_F2xz_S;
  Double I_ERI_G2y2z_S_F2xz_Px = I_ERI_G2y2z_S_G3xz_S+CDX*I_ERI_G2y2z_S_F2xz_S;
  Double I_ERI_Gy3z_S_F2xz_Px = I_ERI_Gy3z_S_G3xz_S+CDX*I_ERI_Gy3z_S_F2xz_S;
  Double I_ERI_G4z_S_F2xz_Px = I_ERI_G4z_S_G3xz_S+CDX*I_ERI_G4z_S_F2xz_S;
  Double I_ERI_G4x_S_Fx2y_Px = I_ERI_G4x_S_G2x2y_S+CDX*I_ERI_G4x_S_Fx2y_S;
  Double I_ERI_G3xy_S_Fx2y_Px = I_ERI_G3xy_S_G2x2y_S+CDX*I_ERI_G3xy_S_Fx2y_S;
  Double I_ERI_G3xz_S_Fx2y_Px = I_ERI_G3xz_S_G2x2y_S+CDX*I_ERI_G3xz_S_Fx2y_S;
  Double I_ERI_G2x2y_S_Fx2y_Px = I_ERI_G2x2y_S_G2x2y_S+CDX*I_ERI_G2x2y_S_Fx2y_S;
  Double I_ERI_G2xyz_S_Fx2y_Px = I_ERI_G2xyz_S_G2x2y_S+CDX*I_ERI_G2xyz_S_Fx2y_S;
  Double I_ERI_G2x2z_S_Fx2y_Px = I_ERI_G2x2z_S_G2x2y_S+CDX*I_ERI_G2x2z_S_Fx2y_S;
  Double I_ERI_Gx3y_S_Fx2y_Px = I_ERI_Gx3y_S_G2x2y_S+CDX*I_ERI_Gx3y_S_Fx2y_S;
  Double I_ERI_Gx2yz_S_Fx2y_Px = I_ERI_Gx2yz_S_G2x2y_S+CDX*I_ERI_Gx2yz_S_Fx2y_S;
  Double I_ERI_Gxy2z_S_Fx2y_Px = I_ERI_Gxy2z_S_G2x2y_S+CDX*I_ERI_Gxy2z_S_Fx2y_S;
  Double I_ERI_Gx3z_S_Fx2y_Px = I_ERI_Gx3z_S_G2x2y_S+CDX*I_ERI_Gx3z_S_Fx2y_S;
  Double I_ERI_G4y_S_Fx2y_Px = I_ERI_G4y_S_G2x2y_S+CDX*I_ERI_G4y_S_Fx2y_S;
  Double I_ERI_G3yz_S_Fx2y_Px = I_ERI_G3yz_S_G2x2y_S+CDX*I_ERI_G3yz_S_Fx2y_S;
  Double I_ERI_G2y2z_S_Fx2y_Px = I_ERI_G2y2z_S_G2x2y_S+CDX*I_ERI_G2y2z_S_Fx2y_S;
  Double I_ERI_Gy3z_S_Fx2y_Px = I_ERI_Gy3z_S_G2x2y_S+CDX*I_ERI_Gy3z_S_Fx2y_S;
  Double I_ERI_G4z_S_Fx2y_Px = I_ERI_G4z_S_G2x2y_S+CDX*I_ERI_G4z_S_Fx2y_S;
  Double I_ERI_G4x_S_Fxyz_Px = I_ERI_G4x_S_G2xyz_S+CDX*I_ERI_G4x_S_Fxyz_S;
  Double I_ERI_G3xy_S_Fxyz_Px = I_ERI_G3xy_S_G2xyz_S+CDX*I_ERI_G3xy_S_Fxyz_S;
  Double I_ERI_G3xz_S_Fxyz_Px = I_ERI_G3xz_S_G2xyz_S+CDX*I_ERI_G3xz_S_Fxyz_S;
  Double I_ERI_G2x2y_S_Fxyz_Px = I_ERI_G2x2y_S_G2xyz_S+CDX*I_ERI_G2x2y_S_Fxyz_S;
  Double I_ERI_G2xyz_S_Fxyz_Px = I_ERI_G2xyz_S_G2xyz_S+CDX*I_ERI_G2xyz_S_Fxyz_S;
  Double I_ERI_G2x2z_S_Fxyz_Px = I_ERI_G2x2z_S_G2xyz_S+CDX*I_ERI_G2x2z_S_Fxyz_S;
  Double I_ERI_Gx3y_S_Fxyz_Px = I_ERI_Gx3y_S_G2xyz_S+CDX*I_ERI_Gx3y_S_Fxyz_S;
  Double I_ERI_Gx2yz_S_Fxyz_Px = I_ERI_Gx2yz_S_G2xyz_S+CDX*I_ERI_Gx2yz_S_Fxyz_S;
  Double I_ERI_Gxy2z_S_Fxyz_Px = I_ERI_Gxy2z_S_G2xyz_S+CDX*I_ERI_Gxy2z_S_Fxyz_S;
  Double I_ERI_Gx3z_S_Fxyz_Px = I_ERI_Gx3z_S_G2xyz_S+CDX*I_ERI_Gx3z_S_Fxyz_S;
  Double I_ERI_G4y_S_Fxyz_Px = I_ERI_G4y_S_G2xyz_S+CDX*I_ERI_G4y_S_Fxyz_S;
  Double I_ERI_G3yz_S_Fxyz_Px = I_ERI_G3yz_S_G2xyz_S+CDX*I_ERI_G3yz_S_Fxyz_S;
  Double I_ERI_G2y2z_S_Fxyz_Px = I_ERI_G2y2z_S_G2xyz_S+CDX*I_ERI_G2y2z_S_Fxyz_S;
  Double I_ERI_Gy3z_S_Fxyz_Px = I_ERI_Gy3z_S_G2xyz_S+CDX*I_ERI_Gy3z_S_Fxyz_S;
  Double I_ERI_G4z_S_Fxyz_Px = I_ERI_G4z_S_G2xyz_S+CDX*I_ERI_G4z_S_Fxyz_S;
  Double I_ERI_G4x_S_Fx2z_Px = I_ERI_G4x_S_G2x2z_S+CDX*I_ERI_G4x_S_Fx2z_S;
  Double I_ERI_G3xy_S_Fx2z_Px = I_ERI_G3xy_S_G2x2z_S+CDX*I_ERI_G3xy_S_Fx2z_S;
  Double I_ERI_G3xz_S_Fx2z_Px = I_ERI_G3xz_S_G2x2z_S+CDX*I_ERI_G3xz_S_Fx2z_S;
  Double I_ERI_G2x2y_S_Fx2z_Px = I_ERI_G2x2y_S_G2x2z_S+CDX*I_ERI_G2x2y_S_Fx2z_S;
  Double I_ERI_G2xyz_S_Fx2z_Px = I_ERI_G2xyz_S_G2x2z_S+CDX*I_ERI_G2xyz_S_Fx2z_S;
  Double I_ERI_G2x2z_S_Fx2z_Px = I_ERI_G2x2z_S_G2x2z_S+CDX*I_ERI_G2x2z_S_Fx2z_S;
  Double I_ERI_Gx3y_S_Fx2z_Px = I_ERI_Gx3y_S_G2x2z_S+CDX*I_ERI_Gx3y_S_Fx2z_S;
  Double I_ERI_Gx2yz_S_Fx2z_Px = I_ERI_Gx2yz_S_G2x2z_S+CDX*I_ERI_Gx2yz_S_Fx2z_S;
  Double I_ERI_Gxy2z_S_Fx2z_Px = I_ERI_Gxy2z_S_G2x2z_S+CDX*I_ERI_Gxy2z_S_Fx2z_S;
  Double I_ERI_Gx3z_S_Fx2z_Px = I_ERI_Gx3z_S_G2x2z_S+CDX*I_ERI_Gx3z_S_Fx2z_S;
  Double I_ERI_G4y_S_Fx2z_Px = I_ERI_G4y_S_G2x2z_S+CDX*I_ERI_G4y_S_Fx2z_S;
  Double I_ERI_G3yz_S_Fx2z_Px = I_ERI_G3yz_S_G2x2z_S+CDX*I_ERI_G3yz_S_Fx2z_S;
  Double I_ERI_G2y2z_S_Fx2z_Px = I_ERI_G2y2z_S_G2x2z_S+CDX*I_ERI_G2y2z_S_Fx2z_S;
  Double I_ERI_Gy3z_S_Fx2z_Px = I_ERI_Gy3z_S_G2x2z_S+CDX*I_ERI_Gy3z_S_Fx2z_S;
  Double I_ERI_G4z_S_Fx2z_Px = I_ERI_G4z_S_G2x2z_S+CDX*I_ERI_G4z_S_Fx2z_S;
  Double I_ERI_G4x_S_F3y_Px = I_ERI_G4x_S_Gx3y_S+CDX*I_ERI_G4x_S_F3y_S;
  Double I_ERI_G3xy_S_F3y_Px = I_ERI_G3xy_S_Gx3y_S+CDX*I_ERI_G3xy_S_F3y_S;
  Double I_ERI_G3xz_S_F3y_Px = I_ERI_G3xz_S_Gx3y_S+CDX*I_ERI_G3xz_S_F3y_S;
  Double I_ERI_G2x2y_S_F3y_Px = I_ERI_G2x2y_S_Gx3y_S+CDX*I_ERI_G2x2y_S_F3y_S;
  Double I_ERI_G2xyz_S_F3y_Px = I_ERI_G2xyz_S_Gx3y_S+CDX*I_ERI_G2xyz_S_F3y_S;
  Double I_ERI_G2x2z_S_F3y_Px = I_ERI_G2x2z_S_Gx3y_S+CDX*I_ERI_G2x2z_S_F3y_S;
  Double I_ERI_Gx3y_S_F3y_Px = I_ERI_Gx3y_S_Gx3y_S+CDX*I_ERI_Gx3y_S_F3y_S;
  Double I_ERI_Gx2yz_S_F3y_Px = I_ERI_Gx2yz_S_Gx3y_S+CDX*I_ERI_Gx2yz_S_F3y_S;
  Double I_ERI_Gxy2z_S_F3y_Px = I_ERI_Gxy2z_S_Gx3y_S+CDX*I_ERI_Gxy2z_S_F3y_S;
  Double I_ERI_Gx3z_S_F3y_Px = I_ERI_Gx3z_S_Gx3y_S+CDX*I_ERI_Gx3z_S_F3y_S;
  Double I_ERI_G4y_S_F3y_Px = I_ERI_G4y_S_Gx3y_S+CDX*I_ERI_G4y_S_F3y_S;
  Double I_ERI_G3yz_S_F3y_Px = I_ERI_G3yz_S_Gx3y_S+CDX*I_ERI_G3yz_S_F3y_S;
  Double I_ERI_G2y2z_S_F3y_Px = I_ERI_G2y2z_S_Gx3y_S+CDX*I_ERI_G2y2z_S_F3y_S;
  Double I_ERI_Gy3z_S_F3y_Px = I_ERI_Gy3z_S_Gx3y_S+CDX*I_ERI_Gy3z_S_F3y_S;
  Double I_ERI_G4z_S_F3y_Px = I_ERI_G4z_S_Gx3y_S+CDX*I_ERI_G4z_S_F3y_S;
  Double I_ERI_G4x_S_F2yz_Px = I_ERI_G4x_S_Gx2yz_S+CDX*I_ERI_G4x_S_F2yz_S;
  Double I_ERI_G3xy_S_F2yz_Px = I_ERI_G3xy_S_Gx2yz_S+CDX*I_ERI_G3xy_S_F2yz_S;
  Double I_ERI_G3xz_S_F2yz_Px = I_ERI_G3xz_S_Gx2yz_S+CDX*I_ERI_G3xz_S_F2yz_S;
  Double I_ERI_G2x2y_S_F2yz_Px = I_ERI_G2x2y_S_Gx2yz_S+CDX*I_ERI_G2x2y_S_F2yz_S;
  Double I_ERI_G2xyz_S_F2yz_Px = I_ERI_G2xyz_S_Gx2yz_S+CDX*I_ERI_G2xyz_S_F2yz_S;
  Double I_ERI_G2x2z_S_F2yz_Px = I_ERI_G2x2z_S_Gx2yz_S+CDX*I_ERI_G2x2z_S_F2yz_S;
  Double I_ERI_Gx3y_S_F2yz_Px = I_ERI_Gx3y_S_Gx2yz_S+CDX*I_ERI_Gx3y_S_F2yz_S;
  Double I_ERI_Gx2yz_S_F2yz_Px = I_ERI_Gx2yz_S_Gx2yz_S+CDX*I_ERI_Gx2yz_S_F2yz_S;
  Double I_ERI_Gxy2z_S_F2yz_Px = I_ERI_Gxy2z_S_Gx2yz_S+CDX*I_ERI_Gxy2z_S_F2yz_S;
  Double I_ERI_Gx3z_S_F2yz_Px = I_ERI_Gx3z_S_Gx2yz_S+CDX*I_ERI_Gx3z_S_F2yz_S;
  Double I_ERI_G4y_S_F2yz_Px = I_ERI_G4y_S_Gx2yz_S+CDX*I_ERI_G4y_S_F2yz_S;
  Double I_ERI_G3yz_S_F2yz_Px = I_ERI_G3yz_S_Gx2yz_S+CDX*I_ERI_G3yz_S_F2yz_S;
  Double I_ERI_G2y2z_S_F2yz_Px = I_ERI_G2y2z_S_Gx2yz_S+CDX*I_ERI_G2y2z_S_F2yz_S;
  Double I_ERI_Gy3z_S_F2yz_Px = I_ERI_Gy3z_S_Gx2yz_S+CDX*I_ERI_Gy3z_S_F2yz_S;
  Double I_ERI_G4z_S_F2yz_Px = I_ERI_G4z_S_Gx2yz_S+CDX*I_ERI_G4z_S_F2yz_S;
  Double I_ERI_G4x_S_Fy2z_Px = I_ERI_G4x_S_Gxy2z_S+CDX*I_ERI_G4x_S_Fy2z_S;
  Double I_ERI_G3xy_S_Fy2z_Px = I_ERI_G3xy_S_Gxy2z_S+CDX*I_ERI_G3xy_S_Fy2z_S;
  Double I_ERI_G3xz_S_Fy2z_Px = I_ERI_G3xz_S_Gxy2z_S+CDX*I_ERI_G3xz_S_Fy2z_S;
  Double I_ERI_G2x2y_S_Fy2z_Px = I_ERI_G2x2y_S_Gxy2z_S+CDX*I_ERI_G2x2y_S_Fy2z_S;
  Double I_ERI_G2xyz_S_Fy2z_Px = I_ERI_G2xyz_S_Gxy2z_S+CDX*I_ERI_G2xyz_S_Fy2z_S;
  Double I_ERI_G2x2z_S_Fy2z_Px = I_ERI_G2x2z_S_Gxy2z_S+CDX*I_ERI_G2x2z_S_Fy2z_S;
  Double I_ERI_Gx3y_S_Fy2z_Px = I_ERI_Gx3y_S_Gxy2z_S+CDX*I_ERI_Gx3y_S_Fy2z_S;
  Double I_ERI_Gx2yz_S_Fy2z_Px = I_ERI_Gx2yz_S_Gxy2z_S+CDX*I_ERI_Gx2yz_S_Fy2z_S;
  Double I_ERI_Gxy2z_S_Fy2z_Px = I_ERI_Gxy2z_S_Gxy2z_S+CDX*I_ERI_Gxy2z_S_Fy2z_S;
  Double I_ERI_Gx3z_S_Fy2z_Px = I_ERI_Gx3z_S_Gxy2z_S+CDX*I_ERI_Gx3z_S_Fy2z_S;
  Double I_ERI_G4y_S_Fy2z_Px = I_ERI_G4y_S_Gxy2z_S+CDX*I_ERI_G4y_S_Fy2z_S;
  Double I_ERI_G3yz_S_Fy2z_Px = I_ERI_G3yz_S_Gxy2z_S+CDX*I_ERI_G3yz_S_Fy2z_S;
  Double I_ERI_G2y2z_S_Fy2z_Px = I_ERI_G2y2z_S_Gxy2z_S+CDX*I_ERI_G2y2z_S_Fy2z_S;
  Double I_ERI_Gy3z_S_Fy2z_Px = I_ERI_Gy3z_S_Gxy2z_S+CDX*I_ERI_Gy3z_S_Fy2z_S;
  Double I_ERI_G4z_S_Fy2z_Px = I_ERI_G4z_S_Gxy2z_S+CDX*I_ERI_G4z_S_Fy2z_S;
  Double I_ERI_G4x_S_F3z_Px = I_ERI_G4x_S_Gx3z_S+CDX*I_ERI_G4x_S_F3z_S;
  Double I_ERI_G3xy_S_F3z_Px = I_ERI_G3xy_S_Gx3z_S+CDX*I_ERI_G3xy_S_F3z_S;
  Double I_ERI_G3xz_S_F3z_Px = I_ERI_G3xz_S_Gx3z_S+CDX*I_ERI_G3xz_S_F3z_S;
  Double I_ERI_G2x2y_S_F3z_Px = I_ERI_G2x2y_S_Gx3z_S+CDX*I_ERI_G2x2y_S_F3z_S;
  Double I_ERI_G2xyz_S_F3z_Px = I_ERI_G2xyz_S_Gx3z_S+CDX*I_ERI_G2xyz_S_F3z_S;
  Double I_ERI_G2x2z_S_F3z_Px = I_ERI_G2x2z_S_Gx3z_S+CDX*I_ERI_G2x2z_S_F3z_S;
  Double I_ERI_Gx3y_S_F3z_Px = I_ERI_Gx3y_S_Gx3z_S+CDX*I_ERI_Gx3y_S_F3z_S;
  Double I_ERI_Gx2yz_S_F3z_Px = I_ERI_Gx2yz_S_Gx3z_S+CDX*I_ERI_Gx2yz_S_F3z_S;
  Double I_ERI_Gxy2z_S_F3z_Px = I_ERI_Gxy2z_S_Gx3z_S+CDX*I_ERI_Gxy2z_S_F3z_S;
  Double I_ERI_Gx3z_S_F3z_Px = I_ERI_Gx3z_S_Gx3z_S+CDX*I_ERI_Gx3z_S_F3z_S;
  Double I_ERI_G4y_S_F3z_Px = I_ERI_G4y_S_Gx3z_S+CDX*I_ERI_G4y_S_F3z_S;
  Double I_ERI_G3yz_S_F3z_Px = I_ERI_G3yz_S_Gx3z_S+CDX*I_ERI_G3yz_S_F3z_S;
  Double I_ERI_G2y2z_S_F3z_Px = I_ERI_G2y2z_S_Gx3z_S+CDX*I_ERI_G2y2z_S_F3z_S;
  Double I_ERI_Gy3z_S_F3z_Px = I_ERI_Gy3z_S_Gx3z_S+CDX*I_ERI_Gy3z_S_F3z_S;
  Double I_ERI_G4z_S_F3z_Px = I_ERI_G4z_S_Gx3z_S+CDX*I_ERI_G4z_S_F3z_S;
  Double I_ERI_G4x_S_F3x_Py = I_ERI_G4x_S_G3xy_S+CDY*I_ERI_G4x_S_F3x_S;
  Double I_ERI_G3xy_S_F3x_Py = I_ERI_G3xy_S_G3xy_S+CDY*I_ERI_G3xy_S_F3x_S;
  Double I_ERI_G3xz_S_F3x_Py = I_ERI_G3xz_S_G3xy_S+CDY*I_ERI_G3xz_S_F3x_S;
  Double I_ERI_G2x2y_S_F3x_Py = I_ERI_G2x2y_S_G3xy_S+CDY*I_ERI_G2x2y_S_F3x_S;
  Double I_ERI_G2xyz_S_F3x_Py = I_ERI_G2xyz_S_G3xy_S+CDY*I_ERI_G2xyz_S_F3x_S;
  Double I_ERI_G2x2z_S_F3x_Py = I_ERI_G2x2z_S_G3xy_S+CDY*I_ERI_G2x2z_S_F3x_S;
  Double I_ERI_Gx3y_S_F3x_Py = I_ERI_Gx3y_S_G3xy_S+CDY*I_ERI_Gx3y_S_F3x_S;
  Double I_ERI_Gx2yz_S_F3x_Py = I_ERI_Gx2yz_S_G3xy_S+CDY*I_ERI_Gx2yz_S_F3x_S;
  Double I_ERI_Gxy2z_S_F3x_Py = I_ERI_Gxy2z_S_G3xy_S+CDY*I_ERI_Gxy2z_S_F3x_S;
  Double I_ERI_Gx3z_S_F3x_Py = I_ERI_Gx3z_S_G3xy_S+CDY*I_ERI_Gx3z_S_F3x_S;
  Double I_ERI_G4y_S_F3x_Py = I_ERI_G4y_S_G3xy_S+CDY*I_ERI_G4y_S_F3x_S;
  Double I_ERI_G3yz_S_F3x_Py = I_ERI_G3yz_S_G3xy_S+CDY*I_ERI_G3yz_S_F3x_S;
  Double I_ERI_G2y2z_S_F3x_Py = I_ERI_G2y2z_S_G3xy_S+CDY*I_ERI_G2y2z_S_F3x_S;
  Double I_ERI_Gy3z_S_F3x_Py = I_ERI_Gy3z_S_G3xy_S+CDY*I_ERI_Gy3z_S_F3x_S;
  Double I_ERI_G4z_S_F3x_Py = I_ERI_G4z_S_G3xy_S+CDY*I_ERI_G4z_S_F3x_S;
  Double I_ERI_G4x_S_F2xy_Py = I_ERI_G4x_S_G2x2y_S+CDY*I_ERI_G4x_S_F2xy_S;
  Double I_ERI_G3xy_S_F2xy_Py = I_ERI_G3xy_S_G2x2y_S+CDY*I_ERI_G3xy_S_F2xy_S;
  Double I_ERI_G3xz_S_F2xy_Py = I_ERI_G3xz_S_G2x2y_S+CDY*I_ERI_G3xz_S_F2xy_S;
  Double I_ERI_G2x2y_S_F2xy_Py = I_ERI_G2x2y_S_G2x2y_S+CDY*I_ERI_G2x2y_S_F2xy_S;
  Double I_ERI_G2xyz_S_F2xy_Py = I_ERI_G2xyz_S_G2x2y_S+CDY*I_ERI_G2xyz_S_F2xy_S;
  Double I_ERI_G2x2z_S_F2xy_Py = I_ERI_G2x2z_S_G2x2y_S+CDY*I_ERI_G2x2z_S_F2xy_S;
  Double I_ERI_Gx3y_S_F2xy_Py = I_ERI_Gx3y_S_G2x2y_S+CDY*I_ERI_Gx3y_S_F2xy_S;
  Double I_ERI_Gx2yz_S_F2xy_Py = I_ERI_Gx2yz_S_G2x2y_S+CDY*I_ERI_Gx2yz_S_F2xy_S;
  Double I_ERI_Gxy2z_S_F2xy_Py = I_ERI_Gxy2z_S_G2x2y_S+CDY*I_ERI_Gxy2z_S_F2xy_S;
  Double I_ERI_Gx3z_S_F2xy_Py = I_ERI_Gx3z_S_G2x2y_S+CDY*I_ERI_Gx3z_S_F2xy_S;
  Double I_ERI_G4y_S_F2xy_Py = I_ERI_G4y_S_G2x2y_S+CDY*I_ERI_G4y_S_F2xy_S;
  Double I_ERI_G3yz_S_F2xy_Py = I_ERI_G3yz_S_G2x2y_S+CDY*I_ERI_G3yz_S_F2xy_S;
  Double I_ERI_G2y2z_S_F2xy_Py = I_ERI_G2y2z_S_G2x2y_S+CDY*I_ERI_G2y2z_S_F2xy_S;
  Double I_ERI_Gy3z_S_F2xy_Py = I_ERI_Gy3z_S_G2x2y_S+CDY*I_ERI_Gy3z_S_F2xy_S;
  Double I_ERI_G4z_S_F2xy_Py = I_ERI_G4z_S_G2x2y_S+CDY*I_ERI_G4z_S_F2xy_S;
  Double I_ERI_G4x_S_F2xz_Py = I_ERI_G4x_S_G2xyz_S+CDY*I_ERI_G4x_S_F2xz_S;
  Double I_ERI_G3xy_S_F2xz_Py = I_ERI_G3xy_S_G2xyz_S+CDY*I_ERI_G3xy_S_F2xz_S;
  Double I_ERI_G3xz_S_F2xz_Py = I_ERI_G3xz_S_G2xyz_S+CDY*I_ERI_G3xz_S_F2xz_S;
  Double I_ERI_G2x2y_S_F2xz_Py = I_ERI_G2x2y_S_G2xyz_S+CDY*I_ERI_G2x2y_S_F2xz_S;
  Double I_ERI_G2xyz_S_F2xz_Py = I_ERI_G2xyz_S_G2xyz_S+CDY*I_ERI_G2xyz_S_F2xz_S;
  Double I_ERI_G2x2z_S_F2xz_Py = I_ERI_G2x2z_S_G2xyz_S+CDY*I_ERI_G2x2z_S_F2xz_S;
  Double I_ERI_Gx3y_S_F2xz_Py = I_ERI_Gx3y_S_G2xyz_S+CDY*I_ERI_Gx3y_S_F2xz_S;
  Double I_ERI_Gx2yz_S_F2xz_Py = I_ERI_Gx2yz_S_G2xyz_S+CDY*I_ERI_Gx2yz_S_F2xz_S;
  Double I_ERI_Gxy2z_S_F2xz_Py = I_ERI_Gxy2z_S_G2xyz_S+CDY*I_ERI_Gxy2z_S_F2xz_S;
  Double I_ERI_Gx3z_S_F2xz_Py = I_ERI_Gx3z_S_G2xyz_S+CDY*I_ERI_Gx3z_S_F2xz_S;
  Double I_ERI_G4y_S_F2xz_Py = I_ERI_G4y_S_G2xyz_S+CDY*I_ERI_G4y_S_F2xz_S;
  Double I_ERI_G3yz_S_F2xz_Py = I_ERI_G3yz_S_G2xyz_S+CDY*I_ERI_G3yz_S_F2xz_S;
  Double I_ERI_G2y2z_S_F2xz_Py = I_ERI_G2y2z_S_G2xyz_S+CDY*I_ERI_G2y2z_S_F2xz_S;
  Double I_ERI_Gy3z_S_F2xz_Py = I_ERI_Gy3z_S_G2xyz_S+CDY*I_ERI_Gy3z_S_F2xz_S;
  Double I_ERI_G4z_S_F2xz_Py = I_ERI_G4z_S_G2xyz_S+CDY*I_ERI_G4z_S_F2xz_S;
  Double I_ERI_G4x_S_Fx2y_Py = I_ERI_G4x_S_Gx3y_S+CDY*I_ERI_G4x_S_Fx2y_S;
  Double I_ERI_G3xy_S_Fx2y_Py = I_ERI_G3xy_S_Gx3y_S+CDY*I_ERI_G3xy_S_Fx2y_S;
  Double I_ERI_G3xz_S_Fx2y_Py = I_ERI_G3xz_S_Gx3y_S+CDY*I_ERI_G3xz_S_Fx2y_S;
  Double I_ERI_G2x2y_S_Fx2y_Py = I_ERI_G2x2y_S_Gx3y_S+CDY*I_ERI_G2x2y_S_Fx2y_S;
  Double I_ERI_G2xyz_S_Fx2y_Py = I_ERI_G2xyz_S_Gx3y_S+CDY*I_ERI_G2xyz_S_Fx2y_S;
  Double I_ERI_G2x2z_S_Fx2y_Py = I_ERI_G2x2z_S_Gx3y_S+CDY*I_ERI_G2x2z_S_Fx2y_S;
  Double I_ERI_Gx3y_S_Fx2y_Py = I_ERI_Gx3y_S_Gx3y_S+CDY*I_ERI_Gx3y_S_Fx2y_S;
  Double I_ERI_Gx2yz_S_Fx2y_Py = I_ERI_Gx2yz_S_Gx3y_S+CDY*I_ERI_Gx2yz_S_Fx2y_S;
  Double I_ERI_Gxy2z_S_Fx2y_Py = I_ERI_Gxy2z_S_Gx3y_S+CDY*I_ERI_Gxy2z_S_Fx2y_S;
  Double I_ERI_Gx3z_S_Fx2y_Py = I_ERI_Gx3z_S_Gx3y_S+CDY*I_ERI_Gx3z_S_Fx2y_S;
  Double I_ERI_G4y_S_Fx2y_Py = I_ERI_G4y_S_Gx3y_S+CDY*I_ERI_G4y_S_Fx2y_S;
  Double I_ERI_G3yz_S_Fx2y_Py = I_ERI_G3yz_S_Gx3y_S+CDY*I_ERI_G3yz_S_Fx2y_S;
  Double I_ERI_G2y2z_S_Fx2y_Py = I_ERI_G2y2z_S_Gx3y_S+CDY*I_ERI_G2y2z_S_Fx2y_S;
  Double I_ERI_Gy3z_S_Fx2y_Py = I_ERI_Gy3z_S_Gx3y_S+CDY*I_ERI_Gy3z_S_Fx2y_S;
  Double I_ERI_G4z_S_Fx2y_Py = I_ERI_G4z_S_Gx3y_S+CDY*I_ERI_G4z_S_Fx2y_S;
  Double I_ERI_G4x_S_Fxyz_Py = I_ERI_G4x_S_Gx2yz_S+CDY*I_ERI_G4x_S_Fxyz_S;
  Double I_ERI_G3xy_S_Fxyz_Py = I_ERI_G3xy_S_Gx2yz_S+CDY*I_ERI_G3xy_S_Fxyz_S;
  Double I_ERI_G3xz_S_Fxyz_Py = I_ERI_G3xz_S_Gx2yz_S+CDY*I_ERI_G3xz_S_Fxyz_S;
  Double I_ERI_G2x2y_S_Fxyz_Py = I_ERI_G2x2y_S_Gx2yz_S+CDY*I_ERI_G2x2y_S_Fxyz_S;
  Double I_ERI_G2xyz_S_Fxyz_Py = I_ERI_G2xyz_S_Gx2yz_S+CDY*I_ERI_G2xyz_S_Fxyz_S;
  Double I_ERI_G2x2z_S_Fxyz_Py = I_ERI_G2x2z_S_Gx2yz_S+CDY*I_ERI_G2x2z_S_Fxyz_S;
  Double I_ERI_Gx3y_S_Fxyz_Py = I_ERI_Gx3y_S_Gx2yz_S+CDY*I_ERI_Gx3y_S_Fxyz_S;
  Double I_ERI_Gx2yz_S_Fxyz_Py = I_ERI_Gx2yz_S_Gx2yz_S+CDY*I_ERI_Gx2yz_S_Fxyz_S;
  Double I_ERI_Gxy2z_S_Fxyz_Py = I_ERI_Gxy2z_S_Gx2yz_S+CDY*I_ERI_Gxy2z_S_Fxyz_S;
  Double I_ERI_Gx3z_S_Fxyz_Py = I_ERI_Gx3z_S_Gx2yz_S+CDY*I_ERI_Gx3z_S_Fxyz_S;
  Double I_ERI_G4y_S_Fxyz_Py = I_ERI_G4y_S_Gx2yz_S+CDY*I_ERI_G4y_S_Fxyz_S;
  Double I_ERI_G3yz_S_Fxyz_Py = I_ERI_G3yz_S_Gx2yz_S+CDY*I_ERI_G3yz_S_Fxyz_S;
  Double I_ERI_G2y2z_S_Fxyz_Py = I_ERI_G2y2z_S_Gx2yz_S+CDY*I_ERI_G2y2z_S_Fxyz_S;
  Double I_ERI_Gy3z_S_Fxyz_Py = I_ERI_Gy3z_S_Gx2yz_S+CDY*I_ERI_Gy3z_S_Fxyz_S;
  Double I_ERI_G4z_S_Fxyz_Py = I_ERI_G4z_S_Gx2yz_S+CDY*I_ERI_G4z_S_Fxyz_S;
  Double I_ERI_G4x_S_Fx2z_Py = I_ERI_G4x_S_Gxy2z_S+CDY*I_ERI_G4x_S_Fx2z_S;
  Double I_ERI_G3xy_S_Fx2z_Py = I_ERI_G3xy_S_Gxy2z_S+CDY*I_ERI_G3xy_S_Fx2z_S;
  Double I_ERI_G3xz_S_Fx2z_Py = I_ERI_G3xz_S_Gxy2z_S+CDY*I_ERI_G3xz_S_Fx2z_S;
  Double I_ERI_G2x2y_S_Fx2z_Py = I_ERI_G2x2y_S_Gxy2z_S+CDY*I_ERI_G2x2y_S_Fx2z_S;
  Double I_ERI_G2xyz_S_Fx2z_Py = I_ERI_G2xyz_S_Gxy2z_S+CDY*I_ERI_G2xyz_S_Fx2z_S;
  Double I_ERI_G2x2z_S_Fx2z_Py = I_ERI_G2x2z_S_Gxy2z_S+CDY*I_ERI_G2x2z_S_Fx2z_S;
  Double I_ERI_Gx3y_S_Fx2z_Py = I_ERI_Gx3y_S_Gxy2z_S+CDY*I_ERI_Gx3y_S_Fx2z_S;
  Double I_ERI_Gx2yz_S_Fx2z_Py = I_ERI_Gx2yz_S_Gxy2z_S+CDY*I_ERI_Gx2yz_S_Fx2z_S;
  Double I_ERI_Gxy2z_S_Fx2z_Py = I_ERI_Gxy2z_S_Gxy2z_S+CDY*I_ERI_Gxy2z_S_Fx2z_S;
  Double I_ERI_Gx3z_S_Fx2z_Py = I_ERI_Gx3z_S_Gxy2z_S+CDY*I_ERI_Gx3z_S_Fx2z_S;
  Double I_ERI_G4y_S_Fx2z_Py = I_ERI_G4y_S_Gxy2z_S+CDY*I_ERI_G4y_S_Fx2z_S;
  Double I_ERI_G3yz_S_Fx2z_Py = I_ERI_G3yz_S_Gxy2z_S+CDY*I_ERI_G3yz_S_Fx2z_S;
  Double I_ERI_G2y2z_S_Fx2z_Py = I_ERI_G2y2z_S_Gxy2z_S+CDY*I_ERI_G2y2z_S_Fx2z_S;
  Double I_ERI_Gy3z_S_Fx2z_Py = I_ERI_Gy3z_S_Gxy2z_S+CDY*I_ERI_Gy3z_S_Fx2z_S;
  Double I_ERI_G4z_S_Fx2z_Py = I_ERI_G4z_S_Gxy2z_S+CDY*I_ERI_G4z_S_Fx2z_S;
  Double I_ERI_G4x_S_F3y_Py = I_ERI_G4x_S_G4y_S+CDY*I_ERI_G4x_S_F3y_S;
  Double I_ERI_G3xy_S_F3y_Py = I_ERI_G3xy_S_G4y_S+CDY*I_ERI_G3xy_S_F3y_S;
  Double I_ERI_G3xz_S_F3y_Py = I_ERI_G3xz_S_G4y_S+CDY*I_ERI_G3xz_S_F3y_S;
  Double I_ERI_G2x2y_S_F3y_Py = I_ERI_G2x2y_S_G4y_S+CDY*I_ERI_G2x2y_S_F3y_S;
  Double I_ERI_G2xyz_S_F3y_Py = I_ERI_G2xyz_S_G4y_S+CDY*I_ERI_G2xyz_S_F3y_S;
  Double I_ERI_G2x2z_S_F3y_Py = I_ERI_G2x2z_S_G4y_S+CDY*I_ERI_G2x2z_S_F3y_S;
  Double I_ERI_Gx3y_S_F3y_Py = I_ERI_Gx3y_S_G4y_S+CDY*I_ERI_Gx3y_S_F3y_S;
  Double I_ERI_Gx2yz_S_F3y_Py = I_ERI_Gx2yz_S_G4y_S+CDY*I_ERI_Gx2yz_S_F3y_S;
  Double I_ERI_Gxy2z_S_F3y_Py = I_ERI_Gxy2z_S_G4y_S+CDY*I_ERI_Gxy2z_S_F3y_S;
  Double I_ERI_Gx3z_S_F3y_Py = I_ERI_Gx3z_S_G4y_S+CDY*I_ERI_Gx3z_S_F3y_S;
  Double I_ERI_G4y_S_F3y_Py = I_ERI_G4y_S_G4y_S+CDY*I_ERI_G4y_S_F3y_S;
  Double I_ERI_G3yz_S_F3y_Py = I_ERI_G3yz_S_G4y_S+CDY*I_ERI_G3yz_S_F3y_S;
  Double I_ERI_G2y2z_S_F3y_Py = I_ERI_G2y2z_S_G4y_S+CDY*I_ERI_G2y2z_S_F3y_S;
  Double I_ERI_Gy3z_S_F3y_Py = I_ERI_Gy3z_S_G4y_S+CDY*I_ERI_Gy3z_S_F3y_S;
  Double I_ERI_G4z_S_F3y_Py = I_ERI_G4z_S_G4y_S+CDY*I_ERI_G4z_S_F3y_S;
  Double I_ERI_G4x_S_F2yz_Py = I_ERI_G4x_S_G3yz_S+CDY*I_ERI_G4x_S_F2yz_S;
  Double I_ERI_G3xy_S_F2yz_Py = I_ERI_G3xy_S_G3yz_S+CDY*I_ERI_G3xy_S_F2yz_S;
  Double I_ERI_G3xz_S_F2yz_Py = I_ERI_G3xz_S_G3yz_S+CDY*I_ERI_G3xz_S_F2yz_S;
  Double I_ERI_G2x2y_S_F2yz_Py = I_ERI_G2x2y_S_G3yz_S+CDY*I_ERI_G2x2y_S_F2yz_S;
  Double I_ERI_G2xyz_S_F2yz_Py = I_ERI_G2xyz_S_G3yz_S+CDY*I_ERI_G2xyz_S_F2yz_S;
  Double I_ERI_G2x2z_S_F2yz_Py = I_ERI_G2x2z_S_G3yz_S+CDY*I_ERI_G2x2z_S_F2yz_S;
  Double I_ERI_Gx3y_S_F2yz_Py = I_ERI_Gx3y_S_G3yz_S+CDY*I_ERI_Gx3y_S_F2yz_S;
  Double I_ERI_Gx2yz_S_F2yz_Py = I_ERI_Gx2yz_S_G3yz_S+CDY*I_ERI_Gx2yz_S_F2yz_S;
  Double I_ERI_Gxy2z_S_F2yz_Py = I_ERI_Gxy2z_S_G3yz_S+CDY*I_ERI_Gxy2z_S_F2yz_S;
  Double I_ERI_Gx3z_S_F2yz_Py = I_ERI_Gx3z_S_G3yz_S+CDY*I_ERI_Gx3z_S_F2yz_S;
  Double I_ERI_G4y_S_F2yz_Py = I_ERI_G4y_S_G3yz_S+CDY*I_ERI_G4y_S_F2yz_S;
  Double I_ERI_G3yz_S_F2yz_Py = I_ERI_G3yz_S_G3yz_S+CDY*I_ERI_G3yz_S_F2yz_S;
  Double I_ERI_G2y2z_S_F2yz_Py = I_ERI_G2y2z_S_G3yz_S+CDY*I_ERI_G2y2z_S_F2yz_S;
  Double I_ERI_Gy3z_S_F2yz_Py = I_ERI_Gy3z_S_G3yz_S+CDY*I_ERI_Gy3z_S_F2yz_S;
  Double I_ERI_G4z_S_F2yz_Py = I_ERI_G4z_S_G3yz_S+CDY*I_ERI_G4z_S_F2yz_S;
  Double I_ERI_G4x_S_Fy2z_Py = I_ERI_G4x_S_G2y2z_S+CDY*I_ERI_G4x_S_Fy2z_S;
  Double I_ERI_G3xy_S_Fy2z_Py = I_ERI_G3xy_S_G2y2z_S+CDY*I_ERI_G3xy_S_Fy2z_S;
  Double I_ERI_G3xz_S_Fy2z_Py = I_ERI_G3xz_S_G2y2z_S+CDY*I_ERI_G3xz_S_Fy2z_S;
  Double I_ERI_G2x2y_S_Fy2z_Py = I_ERI_G2x2y_S_G2y2z_S+CDY*I_ERI_G2x2y_S_Fy2z_S;
  Double I_ERI_G2xyz_S_Fy2z_Py = I_ERI_G2xyz_S_G2y2z_S+CDY*I_ERI_G2xyz_S_Fy2z_S;
  Double I_ERI_G2x2z_S_Fy2z_Py = I_ERI_G2x2z_S_G2y2z_S+CDY*I_ERI_G2x2z_S_Fy2z_S;
  Double I_ERI_Gx3y_S_Fy2z_Py = I_ERI_Gx3y_S_G2y2z_S+CDY*I_ERI_Gx3y_S_Fy2z_S;
  Double I_ERI_Gx2yz_S_Fy2z_Py = I_ERI_Gx2yz_S_G2y2z_S+CDY*I_ERI_Gx2yz_S_Fy2z_S;
  Double I_ERI_Gxy2z_S_Fy2z_Py = I_ERI_Gxy2z_S_G2y2z_S+CDY*I_ERI_Gxy2z_S_Fy2z_S;
  Double I_ERI_Gx3z_S_Fy2z_Py = I_ERI_Gx3z_S_G2y2z_S+CDY*I_ERI_Gx3z_S_Fy2z_S;
  Double I_ERI_G4y_S_Fy2z_Py = I_ERI_G4y_S_G2y2z_S+CDY*I_ERI_G4y_S_Fy2z_S;
  Double I_ERI_G3yz_S_Fy2z_Py = I_ERI_G3yz_S_G2y2z_S+CDY*I_ERI_G3yz_S_Fy2z_S;
  Double I_ERI_G2y2z_S_Fy2z_Py = I_ERI_G2y2z_S_G2y2z_S+CDY*I_ERI_G2y2z_S_Fy2z_S;
  Double I_ERI_Gy3z_S_Fy2z_Py = I_ERI_Gy3z_S_G2y2z_S+CDY*I_ERI_Gy3z_S_Fy2z_S;
  Double I_ERI_G4z_S_Fy2z_Py = I_ERI_G4z_S_G2y2z_S+CDY*I_ERI_G4z_S_Fy2z_S;
  Double I_ERI_G4x_S_F3z_Py = I_ERI_G4x_S_Gy3z_S+CDY*I_ERI_G4x_S_F3z_S;
  Double I_ERI_G3xy_S_F3z_Py = I_ERI_G3xy_S_Gy3z_S+CDY*I_ERI_G3xy_S_F3z_S;
  Double I_ERI_G3xz_S_F3z_Py = I_ERI_G3xz_S_Gy3z_S+CDY*I_ERI_G3xz_S_F3z_S;
  Double I_ERI_G2x2y_S_F3z_Py = I_ERI_G2x2y_S_Gy3z_S+CDY*I_ERI_G2x2y_S_F3z_S;
  Double I_ERI_G2xyz_S_F3z_Py = I_ERI_G2xyz_S_Gy3z_S+CDY*I_ERI_G2xyz_S_F3z_S;
  Double I_ERI_G2x2z_S_F3z_Py = I_ERI_G2x2z_S_Gy3z_S+CDY*I_ERI_G2x2z_S_F3z_S;
  Double I_ERI_Gx3y_S_F3z_Py = I_ERI_Gx3y_S_Gy3z_S+CDY*I_ERI_Gx3y_S_F3z_S;
  Double I_ERI_Gx2yz_S_F3z_Py = I_ERI_Gx2yz_S_Gy3z_S+CDY*I_ERI_Gx2yz_S_F3z_S;
  Double I_ERI_Gxy2z_S_F3z_Py = I_ERI_Gxy2z_S_Gy3z_S+CDY*I_ERI_Gxy2z_S_F3z_S;
  Double I_ERI_Gx3z_S_F3z_Py = I_ERI_Gx3z_S_Gy3z_S+CDY*I_ERI_Gx3z_S_F3z_S;
  Double I_ERI_G4y_S_F3z_Py = I_ERI_G4y_S_Gy3z_S+CDY*I_ERI_G4y_S_F3z_S;
  Double I_ERI_G3yz_S_F3z_Py = I_ERI_G3yz_S_Gy3z_S+CDY*I_ERI_G3yz_S_F3z_S;
  Double I_ERI_G2y2z_S_F3z_Py = I_ERI_G2y2z_S_Gy3z_S+CDY*I_ERI_G2y2z_S_F3z_S;
  Double I_ERI_Gy3z_S_F3z_Py = I_ERI_Gy3z_S_Gy3z_S+CDY*I_ERI_Gy3z_S_F3z_S;
  Double I_ERI_G4z_S_F3z_Py = I_ERI_G4z_S_Gy3z_S+CDY*I_ERI_G4z_S_F3z_S;
  Double I_ERI_G4x_S_F3x_Pz = I_ERI_G4x_S_G3xz_S+CDZ*I_ERI_G4x_S_F3x_S;
  Double I_ERI_G3xy_S_F3x_Pz = I_ERI_G3xy_S_G3xz_S+CDZ*I_ERI_G3xy_S_F3x_S;
  Double I_ERI_G3xz_S_F3x_Pz = I_ERI_G3xz_S_G3xz_S+CDZ*I_ERI_G3xz_S_F3x_S;
  Double I_ERI_G2x2y_S_F3x_Pz = I_ERI_G2x2y_S_G3xz_S+CDZ*I_ERI_G2x2y_S_F3x_S;
  Double I_ERI_G2xyz_S_F3x_Pz = I_ERI_G2xyz_S_G3xz_S+CDZ*I_ERI_G2xyz_S_F3x_S;
  Double I_ERI_G2x2z_S_F3x_Pz = I_ERI_G2x2z_S_G3xz_S+CDZ*I_ERI_G2x2z_S_F3x_S;
  Double I_ERI_Gx3y_S_F3x_Pz = I_ERI_Gx3y_S_G3xz_S+CDZ*I_ERI_Gx3y_S_F3x_S;
  Double I_ERI_Gx2yz_S_F3x_Pz = I_ERI_Gx2yz_S_G3xz_S+CDZ*I_ERI_Gx2yz_S_F3x_S;
  Double I_ERI_Gxy2z_S_F3x_Pz = I_ERI_Gxy2z_S_G3xz_S+CDZ*I_ERI_Gxy2z_S_F3x_S;
  Double I_ERI_Gx3z_S_F3x_Pz = I_ERI_Gx3z_S_G3xz_S+CDZ*I_ERI_Gx3z_S_F3x_S;
  Double I_ERI_G4y_S_F3x_Pz = I_ERI_G4y_S_G3xz_S+CDZ*I_ERI_G4y_S_F3x_S;
  Double I_ERI_G3yz_S_F3x_Pz = I_ERI_G3yz_S_G3xz_S+CDZ*I_ERI_G3yz_S_F3x_S;
  Double I_ERI_G2y2z_S_F3x_Pz = I_ERI_G2y2z_S_G3xz_S+CDZ*I_ERI_G2y2z_S_F3x_S;
  Double I_ERI_Gy3z_S_F3x_Pz = I_ERI_Gy3z_S_G3xz_S+CDZ*I_ERI_Gy3z_S_F3x_S;
  Double I_ERI_G4z_S_F3x_Pz = I_ERI_G4z_S_G3xz_S+CDZ*I_ERI_G4z_S_F3x_S;
  Double I_ERI_G4x_S_F2xy_Pz = I_ERI_G4x_S_G2xyz_S+CDZ*I_ERI_G4x_S_F2xy_S;
  Double I_ERI_G3xy_S_F2xy_Pz = I_ERI_G3xy_S_G2xyz_S+CDZ*I_ERI_G3xy_S_F2xy_S;
  Double I_ERI_G3xz_S_F2xy_Pz = I_ERI_G3xz_S_G2xyz_S+CDZ*I_ERI_G3xz_S_F2xy_S;
  Double I_ERI_G2x2y_S_F2xy_Pz = I_ERI_G2x2y_S_G2xyz_S+CDZ*I_ERI_G2x2y_S_F2xy_S;
  Double I_ERI_G2xyz_S_F2xy_Pz = I_ERI_G2xyz_S_G2xyz_S+CDZ*I_ERI_G2xyz_S_F2xy_S;
  Double I_ERI_G2x2z_S_F2xy_Pz = I_ERI_G2x2z_S_G2xyz_S+CDZ*I_ERI_G2x2z_S_F2xy_S;
  Double I_ERI_Gx3y_S_F2xy_Pz = I_ERI_Gx3y_S_G2xyz_S+CDZ*I_ERI_Gx3y_S_F2xy_S;
  Double I_ERI_Gx2yz_S_F2xy_Pz = I_ERI_Gx2yz_S_G2xyz_S+CDZ*I_ERI_Gx2yz_S_F2xy_S;
  Double I_ERI_Gxy2z_S_F2xy_Pz = I_ERI_Gxy2z_S_G2xyz_S+CDZ*I_ERI_Gxy2z_S_F2xy_S;
  Double I_ERI_Gx3z_S_F2xy_Pz = I_ERI_Gx3z_S_G2xyz_S+CDZ*I_ERI_Gx3z_S_F2xy_S;
  Double I_ERI_G4y_S_F2xy_Pz = I_ERI_G4y_S_G2xyz_S+CDZ*I_ERI_G4y_S_F2xy_S;
  Double I_ERI_G3yz_S_F2xy_Pz = I_ERI_G3yz_S_G2xyz_S+CDZ*I_ERI_G3yz_S_F2xy_S;
  Double I_ERI_G2y2z_S_F2xy_Pz = I_ERI_G2y2z_S_G2xyz_S+CDZ*I_ERI_G2y2z_S_F2xy_S;
  Double I_ERI_Gy3z_S_F2xy_Pz = I_ERI_Gy3z_S_G2xyz_S+CDZ*I_ERI_Gy3z_S_F2xy_S;
  Double I_ERI_G4z_S_F2xy_Pz = I_ERI_G4z_S_G2xyz_S+CDZ*I_ERI_G4z_S_F2xy_S;
  Double I_ERI_G4x_S_F2xz_Pz = I_ERI_G4x_S_G2x2z_S+CDZ*I_ERI_G4x_S_F2xz_S;
  Double I_ERI_G3xy_S_F2xz_Pz = I_ERI_G3xy_S_G2x2z_S+CDZ*I_ERI_G3xy_S_F2xz_S;
  Double I_ERI_G3xz_S_F2xz_Pz = I_ERI_G3xz_S_G2x2z_S+CDZ*I_ERI_G3xz_S_F2xz_S;
  Double I_ERI_G2x2y_S_F2xz_Pz = I_ERI_G2x2y_S_G2x2z_S+CDZ*I_ERI_G2x2y_S_F2xz_S;
  Double I_ERI_G2xyz_S_F2xz_Pz = I_ERI_G2xyz_S_G2x2z_S+CDZ*I_ERI_G2xyz_S_F2xz_S;
  Double I_ERI_G2x2z_S_F2xz_Pz = I_ERI_G2x2z_S_G2x2z_S+CDZ*I_ERI_G2x2z_S_F2xz_S;
  Double I_ERI_Gx3y_S_F2xz_Pz = I_ERI_Gx3y_S_G2x2z_S+CDZ*I_ERI_Gx3y_S_F2xz_S;
  Double I_ERI_Gx2yz_S_F2xz_Pz = I_ERI_Gx2yz_S_G2x2z_S+CDZ*I_ERI_Gx2yz_S_F2xz_S;
  Double I_ERI_Gxy2z_S_F2xz_Pz = I_ERI_Gxy2z_S_G2x2z_S+CDZ*I_ERI_Gxy2z_S_F2xz_S;
  Double I_ERI_Gx3z_S_F2xz_Pz = I_ERI_Gx3z_S_G2x2z_S+CDZ*I_ERI_Gx3z_S_F2xz_S;
  Double I_ERI_G4y_S_F2xz_Pz = I_ERI_G4y_S_G2x2z_S+CDZ*I_ERI_G4y_S_F2xz_S;
  Double I_ERI_G3yz_S_F2xz_Pz = I_ERI_G3yz_S_G2x2z_S+CDZ*I_ERI_G3yz_S_F2xz_S;
  Double I_ERI_G2y2z_S_F2xz_Pz = I_ERI_G2y2z_S_G2x2z_S+CDZ*I_ERI_G2y2z_S_F2xz_S;
  Double I_ERI_Gy3z_S_F2xz_Pz = I_ERI_Gy3z_S_G2x2z_S+CDZ*I_ERI_Gy3z_S_F2xz_S;
  Double I_ERI_G4z_S_F2xz_Pz = I_ERI_G4z_S_G2x2z_S+CDZ*I_ERI_G4z_S_F2xz_S;
  Double I_ERI_G4x_S_Fx2y_Pz = I_ERI_G4x_S_Gx2yz_S+CDZ*I_ERI_G4x_S_Fx2y_S;
  Double I_ERI_G3xy_S_Fx2y_Pz = I_ERI_G3xy_S_Gx2yz_S+CDZ*I_ERI_G3xy_S_Fx2y_S;
  Double I_ERI_G3xz_S_Fx2y_Pz = I_ERI_G3xz_S_Gx2yz_S+CDZ*I_ERI_G3xz_S_Fx2y_S;
  Double I_ERI_G2x2y_S_Fx2y_Pz = I_ERI_G2x2y_S_Gx2yz_S+CDZ*I_ERI_G2x2y_S_Fx2y_S;
  Double I_ERI_G2xyz_S_Fx2y_Pz = I_ERI_G2xyz_S_Gx2yz_S+CDZ*I_ERI_G2xyz_S_Fx2y_S;
  Double I_ERI_G2x2z_S_Fx2y_Pz = I_ERI_G2x2z_S_Gx2yz_S+CDZ*I_ERI_G2x2z_S_Fx2y_S;
  Double I_ERI_Gx3y_S_Fx2y_Pz = I_ERI_Gx3y_S_Gx2yz_S+CDZ*I_ERI_Gx3y_S_Fx2y_S;
  Double I_ERI_Gx2yz_S_Fx2y_Pz = I_ERI_Gx2yz_S_Gx2yz_S+CDZ*I_ERI_Gx2yz_S_Fx2y_S;
  Double I_ERI_Gxy2z_S_Fx2y_Pz = I_ERI_Gxy2z_S_Gx2yz_S+CDZ*I_ERI_Gxy2z_S_Fx2y_S;
  Double I_ERI_Gx3z_S_Fx2y_Pz = I_ERI_Gx3z_S_Gx2yz_S+CDZ*I_ERI_Gx3z_S_Fx2y_S;
  Double I_ERI_G4y_S_Fx2y_Pz = I_ERI_G4y_S_Gx2yz_S+CDZ*I_ERI_G4y_S_Fx2y_S;
  Double I_ERI_G3yz_S_Fx2y_Pz = I_ERI_G3yz_S_Gx2yz_S+CDZ*I_ERI_G3yz_S_Fx2y_S;
  Double I_ERI_G2y2z_S_Fx2y_Pz = I_ERI_G2y2z_S_Gx2yz_S+CDZ*I_ERI_G2y2z_S_Fx2y_S;
  Double I_ERI_Gy3z_S_Fx2y_Pz = I_ERI_Gy3z_S_Gx2yz_S+CDZ*I_ERI_Gy3z_S_Fx2y_S;
  Double I_ERI_G4z_S_Fx2y_Pz = I_ERI_G4z_S_Gx2yz_S+CDZ*I_ERI_G4z_S_Fx2y_S;
  Double I_ERI_G4x_S_Fxyz_Pz = I_ERI_G4x_S_Gxy2z_S+CDZ*I_ERI_G4x_S_Fxyz_S;
  Double I_ERI_G3xy_S_Fxyz_Pz = I_ERI_G3xy_S_Gxy2z_S+CDZ*I_ERI_G3xy_S_Fxyz_S;
  Double I_ERI_G3xz_S_Fxyz_Pz = I_ERI_G3xz_S_Gxy2z_S+CDZ*I_ERI_G3xz_S_Fxyz_S;
  Double I_ERI_G2x2y_S_Fxyz_Pz = I_ERI_G2x2y_S_Gxy2z_S+CDZ*I_ERI_G2x2y_S_Fxyz_S;
  Double I_ERI_G2xyz_S_Fxyz_Pz = I_ERI_G2xyz_S_Gxy2z_S+CDZ*I_ERI_G2xyz_S_Fxyz_S;
  Double I_ERI_G2x2z_S_Fxyz_Pz = I_ERI_G2x2z_S_Gxy2z_S+CDZ*I_ERI_G2x2z_S_Fxyz_S;
  Double I_ERI_Gx3y_S_Fxyz_Pz = I_ERI_Gx3y_S_Gxy2z_S+CDZ*I_ERI_Gx3y_S_Fxyz_S;
  Double I_ERI_Gx2yz_S_Fxyz_Pz = I_ERI_Gx2yz_S_Gxy2z_S+CDZ*I_ERI_Gx2yz_S_Fxyz_S;
  Double I_ERI_Gxy2z_S_Fxyz_Pz = I_ERI_Gxy2z_S_Gxy2z_S+CDZ*I_ERI_Gxy2z_S_Fxyz_S;
  Double I_ERI_Gx3z_S_Fxyz_Pz = I_ERI_Gx3z_S_Gxy2z_S+CDZ*I_ERI_Gx3z_S_Fxyz_S;
  Double I_ERI_G4y_S_Fxyz_Pz = I_ERI_G4y_S_Gxy2z_S+CDZ*I_ERI_G4y_S_Fxyz_S;
  Double I_ERI_G3yz_S_Fxyz_Pz = I_ERI_G3yz_S_Gxy2z_S+CDZ*I_ERI_G3yz_S_Fxyz_S;
  Double I_ERI_G2y2z_S_Fxyz_Pz = I_ERI_G2y2z_S_Gxy2z_S+CDZ*I_ERI_G2y2z_S_Fxyz_S;
  Double I_ERI_Gy3z_S_Fxyz_Pz = I_ERI_Gy3z_S_Gxy2z_S+CDZ*I_ERI_Gy3z_S_Fxyz_S;
  Double I_ERI_G4z_S_Fxyz_Pz = I_ERI_G4z_S_Gxy2z_S+CDZ*I_ERI_G4z_S_Fxyz_S;
  Double I_ERI_G4x_S_Fx2z_Pz = I_ERI_G4x_S_Gx3z_S+CDZ*I_ERI_G4x_S_Fx2z_S;
  Double I_ERI_G3xy_S_Fx2z_Pz = I_ERI_G3xy_S_Gx3z_S+CDZ*I_ERI_G3xy_S_Fx2z_S;
  Double I_ERI_G3xz_S_Fx2z_Pz = I_ERI_G3xz_S_Gx3z_S+CDZ*I_ERI_G3xz_S_Fx2z_S;
  Double I_ERI_G2x2y_S_Fx2z_Pz = I_ERI_G2x2y_S_Gx3z_S+CDZ*I_ERI_G2x2y_S_Fx2z_S;
  Double I_ERI_G2xyz_S_Fx2z_Pz = I_ERI_G2xyz_S_Gx3z_S+CDZ*I_ERI_G2xyz_S_Fx2z_S;
  Double I_ERI_G2x2z_S_Fx2z_Pz = I_ERI_G2x2z_S_Gx3z_S+CDZ*I_ERI_G2x2z_S_Fx2z_S;
  Double I_ERI_Gx3y_S_Fx2z_Pz = I_ERI_Gx3y_S_Gx3z_S+CDZ*I_ERI_Gx3y_S_Fx2z_S;
  Double I_ERI_Gx2yz_S_Fx2z_Pz = I_ERI_Gx2yz_S_Gx3z_S+CDZ*I_ERI_Gx2yz_S_Fx2z_S;
  Double I_ERI_Gxy2z_S_Fx2z_Pz = I_ERI_Gxy2z_S_Gx3z_S+CDZ*I_ERI_Gxy2z_S_Fx2z_S;
  Double I_ERI_Gx3z_S_Fx2z_Pz = I_ERI_Gx3z_S_Gx3z_S+CDZ*I_ERI_Gx3z_S_Fx2z_S;
  Double I_ERI_G4y_S_Fx2z_Pz = I_ERI_G4y_S_Gx3z_S+CDZ*I_ERI_G4y_S_Fx2z_S;
  Double I_ERI_G3yz_S_Fx2z_Pz = I_ERI_G3yz_S_Gx3z_S+CDZ*I_ERI_G3yz_S_Fx2z_S;
  Double I_ERI_G2y2z_S_Fx2z_Pz = I_ERI_G2y2z_S_Gx3z_S+CDZ*I_ERI_G2y2z_S_Fx2z_S;
  Double I_ERI_Gy3z_S_Fx2z_Pz = I_ERI_Gy3z_S_Gx3z_S+CDZ*I_ERI_Gy3z_S_Fx2z_S;
  Double I_ERI_G4z_S_Fx2z_Pz = I_ERI_G4z_S_Gx3z_S+CDZ*I_ERI_G4z_S_Fx2z_S;
  Double I_ERI_G4x_S_F3y_Pz = I_ERI_G4x_S_G3yz_S+CDZ*I_ERI_G4x_S_F3y_S;
  Double I_ERI_G3xy_S_F3y_Pz = I_ERI_G3xy_S_G3yz_S+CDZ*I_ERI_G3xy_S_F3y_S;
  Double I_ERI_G3xz_S_F3y_Pz = I_ERI_G3xz_S_G3yz_S+CDZ*I_ERI_G3xz_S_F3y_S;
  Double I_ERI_G2x2y_S_F3y_Pz = I_ERI_G2x2y_S_G3yz_S+CDZ*I_ERI_G2x2y_S_F3y_S;
  Double I_ERI_G2xyz_S_F3y_Pz = I_ERI_G2xyz_S_G3yz_S+CDZ*I_ERI_G2xyz_S_F3y_S;
  Double I_ERI_G2x2z_S_F3y_Pz = I_ERI_G2x2z_S_G3yz_S+CDZ*I_ERI_G2x2z_S_F3y_S;
  Double I_ERI_Gx3y_S_F3y_Pz = I_ERI_Gx3y_S_G3yz_S+CDZ*I_ERI_Gx3y_S_F3y_S;
  Double I_ERI_Gx2yz_S_F3y_Pz = I_ERI_Gx2yz_S_G3yz_S+CDZ*I_ERI_Gx2yz_S_F3y_S;
  Double I_ERI_Gxy2z_S_F3y_Pz = I_ERI_Gxy2z_S_G3yz_S+CDZ*I_ERI_Gxy2z_S_F3y_S;
  Double I_ERI_Gx3z_S_F3y_Pz = I_ERI_Gx3z_S_G3yz_S+CDZ*I_ERI_Gx3z_S_F3y_S;
  Double I_ERI_G4y_S_F3y_Pz = I_ERI_G4y_S_G3yz_S+CDZ*I_ERI_G4y_S_F3y_S;
  Double I_ERI_G3yz_S_F3y_Pz = I_ERI_G3yz_S_G3yz_S+CDZ*I_ERI_G3yz_S_F3y_S;
  Double I_ERI_G2y2z_S_F3y_Pz = I_ERI_G2y2z_S_G3yz_S+CDZ*I_ERI_G2y2z_S_F3y_S;
  Double I_ERI_Gy3z_S_F3y_Pz = I_ERI_Gy3z_S_G3yz_S+CDZ*I_ERI_Gy3z_S_F3y_S;
  Double I_ERI_G4z_S_F3y_Pz = I_ERI_G4z_S_G3yz_S+CDZ*I_ERI_G4z_S_F3y_S;
  Double I_ERI_G4x_S_F2yz_Pz = I_ERI_G4x_S_G2y2z_S+CDZ*I_ERI_G4x_S_F2yz_S;
  Double I_ERI_G3xy_S_F2yz_Pz = I_ERI_G3xy_S_G2y2z_S+CDZ*I_ERI_G3xy_S_F2yz_S;
  Double I_ERI_G3xz_S_F2yz_Pz = I_ERI_G3xz_S_G2y2z_S+CDZ*I_ERI_G3xz_S_F2yz_S;
  Double I_ERI_G2x2y_S_F2yz_Pz = I_ERI_G2x2y_S_G2y2z_S+CDZ*I_ERI_G2x2y_S_F2yz_S;
  Double I_ERI_G2xyz_S_F2yz_Pz = I_ERI_G2xyz_S_G2y2z_S+CDZ*I_ERI_G2xyz_S_F2yz_S;
  Double I_ERI_G2x2z_S_F2yz_Pz = I_ERI_G2x2z_S_G2y2z_S+CDZ*I_ERI_G2x2z_S_F2yz_S;
  Double I_ERI_Gx3y_S_F2yz_Pz = I_ERI_Gx3y_S_G2y2z_S+CDZ*I_ERI_Gx3y_S_F2yz_S;
  Double I_ERI_Gx2yz_S_F2yz_Pz = I_ERI_Gx2yz_S_G2y2z_S+CDZ*I_ERI_Gx2yz_S_F2yz_S;
  Double I_ERI_Gxy2z_S_F2yz_Pz = I_ERI_Gxy2z_S_G2y2z_S+CDZ*I_ERI_Gxy2z_S_F2yz_S;
  Double I_ERI_Gx3z_S_F2yz_Pz = I_ERI_Gx3z_S_G2y2z_S+CDZ*I_ERI_Gx3z_S_F2yz_S;
  Double I_ERI_G4y_S_F2yz_Pz = I_ERI_G4y_S_G2y2z_S+CDZ*I_ERI_G4y_S_F2yz_S;
  Double I_ERI_G3yz_S_F2yz_Pz = I_ERI_G3yz_S_G2y2z_S+CDZ*I_ERI_G3yz_S_F2yz_S;
  Double I_ERI_G2y2z_S_F2yz_Pz = I_ERI_G2y2z_S_G2y2z_S+CDZ*I_ERI_G2y2z_S_F2yz_S;
  Double I_ERI_Gy3z_S_F2yz_Pz = I_ERI_Gy3z_S_G2y2z_S+CDZ*I_ERI_Gy3z_S_F2yz_S;
  Double I_ERI_G4z_S_F2yz_Pz = I_ERI_G4z_S_G2y2z_S+CDZ*I_ERI_G4z_S_F2yz_S;
  Double I_ERI_G4x_S_Fy2z_Pz = I_ERI_G4x_S_Gy3z_S+CDZ*I_ERI_G4x_S_Fy2z_S;
  Double I_ERI_G3xy_S_Fy2z_Pz = I_ERI_G3xy_S_Gy3z_S+CDZ*I_ERI_G3xy_S_Fy2z_S;
  Double I_ERI_G3xz_S_Fy2z_Pz = I_ERI_G3xz_S_Gy3z_S+CDZ*I_ERI_G3xz_S_Fy2z_S;
  Double I_ERI_G2x2y_S_Fy2z_Pz = I_ERI_G2x2y_S_Gy3z_S+CDZ*I_ERI_G2x2y_S_Fy2z_S;
  Double I_ERI_G2xyz_S_Fy2z_Pz = I_ERI_G2xyz_S_Gy3z_S+CDZ*I_ERI_G2xyz_S_Fy2z_S;
  Double I_ERI_G2x2z_S_Fy2z_Pz = I_ERI_G2x2z_S_Gy3z_S+CDZ*I_ERI_G2x2z_S_Fy2z_S;
  Double I_ERI_Gx3y_S_Fy2z_Pz = I_ERI_Gx3y_S_Gy3z_S+CDZ*I_ERI_Gx3y_S_Fy2z_S;
  Double I_ERI_Gx2yz_S_Fy2z_Pz = I_ERI_Gx2yz_S_Gy3z_S+CDZ*I_ERI_Gx2yz_S_Fy2z_S;
  Double I_ERI_Gxy2z_S_Fy2z_Pz = I_ERI_Gxy2z_S_Gy3z_S+CDZ*I_ERI_Gxy2z_S_Fy2z_S;
  Double I_ERI_Gx3z_S_Fy2z_Pz = I_ERI_Gx3z_S_Gy3z_S+CDZ*I_ERI_Gx3z_S_Fy2z_S;
  Double I_ERI_G4y_S_Fy2z_Pz = I_ERI_G4y_S_Gy3z_S+CDZ*I_ERI_G4y_S_Fy2z_S;
  Double I_ERI_G3yz_S_Fy2z_Pz = I_ERI_G3yz_S_Gy3z_S+CDZ*I_ERI_G3yz_S_Fy2z_S;
  Double I_ERI_G2y2z_S_Fy2z_Pz = I_ERI_G2y2z_S_Gy3z_S+CDZ*I_ERI_G2y2z_S_Fy2z_S;
  Double I_ERI_Gy3z_S_Fy2z_Pz = I_ERI_Gy3z_S_Gy3z_S+CDZ*I_ERI_Gy3z_S_Fy2z_S;
  Double I_ERI_G4z_S_Fy2z_Pz = I_ERI_G4z_S_Gy3z_S+CDZ*I_ERI_G4z_S_Fy2z_S;
  Double I_ERI_G4x_S_F3z_Pz = I_ERI_G4x_S_G4z_S+CDZ*I_ERI_G4x_S_F3z_S;
  Double I_ERI_G3xy_S_F3z_Pz = I_ERI_G3xy_S_G4z_S+CDZ*I_ERI_G3xy_S_F3z_S;
  Double I_ERI_G3xz_S_F3z_Pz = I_ERI_G3xz_S_G4z_S+CDZ*I_ERI_G3xz_S_F3z_S;
  Double I_ERI_G2x2y_S_F3z_Pz = I_ERI_G2x2y_S_G4z_S+CDZ*I_ERI_G2x2y_S_F3z_S;
  Double I_ERI_G2xyz_S_F3z_Pz = I_ERI_G2xyz_S_G4z_S+CDZ*I_ERI_G2xyz_S_F3z_S;
  Double I_ERI_G2x2z_S_F3z_Pz = I_ERI_G2x2z_S_G4z_S+CDZ*I_ERI_G2x2z_S_F3z_S;
  Double I_ERI_Gx3y_S_F3z_Pz = I_ERI_Gx3y_S_G4z_S+CDZ*I_ERI_Gx3y_S_F3z_S;
  Double I_ERI_Gx2yz_S_F3z_Pz = I_ERI_Gx2yz_S_G4z_S+CDZ*I_ERI_Gx2yz_S_F3z_S;
  Double I_ERI_Gxy2z_S_F3z_Pz = I_ERI_Gxy2z_S_G4z_S+CDZ*I_ERI_Gxy2z_S_F3z_S;
  Double I_ERI_Gx3z_S_F3z_Pz = I_ERI_Gx3z_S_G4z_S+CDZ*I_ERI_Gx3z_S_F3z_S;
  Double I_ERI_G4y_S_F3z_Pz = I_ERI_G4y_S_G4z_S+CDZ*I_ERI_G4y_S_F3z_S;
  Double I_ERI_G3yz_S_F3z_Pz = I_ERI_G3yz_S_G4z_S+CDZ*I_ERI_G3yz_S_F3z_S;
  Double I_ERI_G2y2z_S_F3z_Pz = I_ERI_G2y2z_S_G4z_S+CDZ*I_ERI_G2y2z_S_F3z_S;
  Double I_ERI_Gy3z_S_F3z_Pz = I_ERI_Gy3z_S_G4z_S+CDZ*I_ERI_Gy3z_S_F3z_S;
  Double I_ERI_G4z_S_F3z_Pz = I_ERI_G4z_S_G4z_S+CDZ*I_ERI_G4z_S_F3z_S;

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
   * shell quartet name: SQ_ERI_D_P_F_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_P
   * RHS shell quartet name: SQ_ERI_D_S_F_P
   ************************************************************/
  Double I_ERI_D2x_Px_F3x_Px = I_ERI_F3x_S_F3x_Px+ABX*I_ERI_D2x_S_F3x_Px;
  Double I_ERI_Dxy_Px_F3x_Px = I_ERI_F2xy_S_F3x_Px+ABX*I_ERI_Dxy_S_F3x_Px;
  Double I_ERI_Dxz_Px_F3x_Px = I_ERI_F2xz_S_F3x_Px+ABX*I_ERI_Dxz_S_F3x_Px;
  Double I_ERI_D2y_Px_F3x_Px = I_ERI_Fx2y_S_F3x_Px+ABX*I_ERI_D2y_S_F3x_Px;
  Double I_ERI_Dyz_Px_F3x_Px = I_ERI_Fxyz_S_F3x_Px+ABX*I_ERI_Dyz_S_F3x_Px;
  Double I_ERI_D2z_Px_F3x_Px = I_ERI_Fx2z_S_F3x_Px+ABX*I_ERI_D2z_S_F3x_Px;
  Double I_ERI_D2x_Py_F3x_Px = I_ERI_F2xy_S_F3x_Px+ABY*I_ERI_D2x_S_F3x_Px;
  Double I_ERI_Dxy_Py_F3x_Px = I_ERI_Fx2y_S_F3x_Px+ABY*I_ERI_Dxy_S_F3x_Px;
  Double I_ERI_Dxz_Py_F3x_Px = I_ERI_Fxyz_S_F3x_Px+ABY*I_ERI_Dxz_S_F3x_Px;
  Double I_ERI_D2y_Py_F3x_Px = I_ERI_F3y_S_F3x_Px+ABY*I_ERI_D2y_S_F3x_Px;
  Double I_ERI_Dyz_Py_F3x_Px = I_ERI_F2yz_S_F3x_Px+ABY*I_ERI_Dyz_S_F3x_Px;
  Double I_ERI_D2z_Py_F3x_Px = I_ERI_Fy2z_S_F3x_Px+ABY*I_ERI_D2z_S_F3x_Px;
  Double I_ERI_D2x_Pz_F3x_Px = I_ERI_F2xz_S_F3x_Px+ABZ*I_ERI_D2x_S_F3x_Px;
  Double I_ERI_Dxy_Pz_F3x_Px = I_ERI_Fxyz_S_F3x_Px+ABZ*I_ERI_Dxy_S_F3x_Px;
  Double I_ERI_Dxz_Pz_F3x_Px = I_ERI_Fx2z_S_F3x_Px+ABZ*I_ERI_Dxz_S_F3x_Px;
  Double I_ERI_D2y_Pz_F3x_Px = I_ERI_F2yz_S_F3x_Px+ABZ*I_ERI_D2y_S_F3x_Px;
  Double I_ERI_Dyz_Pz_F3x_Px = I_ERI_Fy2z_S_F3x_Px+ABZ*I_ERI_Dyz_S_F3x_Px;
  Double I_ERI_D2z_Pz_F3x_Px = I_ERI_F3z_S_F3x_Px+ABZ*I_ERI_D2z_S_F3x_Px;
  Double I_ERI_D2x_Px_F2xy_Px = I_ERI_F3x_S_F2xy_Px+ABX*I_ERI_D2x_S_F2xy_Px;
  Double I_ERI_Dxy_Px_F2xy_Px = I_ERI_F2xy_S_F2xy_Px+ABX*I_ERI_Dxy_S_F2xy_Px;
  Double I_ERI_Dxz_Px_F2xy_Px = I_ERI_F2xz_S_F2xy_Px+ABX*I_ERI_Dxz_S_F2xy_Px;
  Double I_ERI_D2y_Px_F2xy_Px = I_ERI_Fx2y_S_F2xy_Px+ABX*I_ERI_D2y_S_F2xy_Px;
  Double I_ERI_Dyz_Px_F2xy_Px = I_ERI_Fxyz_S_F2xy_Px+ABX*I_ERI_Dyz_S_F2xy_Px;
  Double I_ERI_D2z_Px_F2xy_Px = I_ERI_Fx2z_S_F2xy_Px+ABX*I_ERI_D2z_S_F2xy_Px;
  Double I_ERI_D2x_Py_F2xy_Px = I_ERI_F2xy_S_F2xy_Px+ABY*I_ERI_D2x_S_F2xy_Px;
  Double I_ERI_Dxy_Py_F2xy_Px = I_ERI_Fx2y_S_F2xy_Px+ABY*I_ERI_Dxy_S_F2xy_Px;
  Double I_ERI_Dxz_Py_F2xy_Px = I_ERI_Fxyz_S_F2xy_Px+ABY*I_ERI_Dxz_S_F2xy_Px;
  Double I_ERI_D2y_Py_F2xy_Px = I_ERI_F3y_S_F2xy_Px+ABY*I_ERI_D2y_S_F2xy_Px;
  Double I_ERI_Dyz_Py_F2xy_Px = I_ERI_F2yz_S_F2xy_Px+ABY*I_ERI_Dyz_S_F2xy_Px;
  Double I_ERI_D2z_Py_F2xy_Px = I_ERI_Fy2z_S_F2xy_Px+ABY*I_ERI_D2z_S_F2xy_Px;
  Double I_ERI_D2x_Pz_F2xy_Px = I_ERI_F2xz_S_F2xy_Px+ABZ*I_ERI_D2x_S_F2xy_Px;
  Double I_ERI_Dxy_Pz_F2xy_Px = I_ERI_Fxyz_S_F2xy_Px+ABZ*I_ERI_Dxy_S_F2xy_Px;
  Double I_ERI_Dxz_Pz_F2xy_Px = I_ERI_Fx2z_S_F2xy_Px+ABZ*I_ERI_Dxz_S_F2xy_Px;
  Double I_ERI_D2y_Pz_F2xy_Px = I_ERI_F2yz_S_F2xy_Px+ABZ*I_ERI_D2y_S_F2xy_Px;
  Double I_ERI_Dyz_Pz_F2xy_Px = I_ERI_Fy2z_S_F2xy_Px+ABZ*I_ERI_Dyz_S_F2xy_Px;
  Double I_ERI_D2z_Pz_F2xy_Px = I_ERI_F3z_S_F2xy_Px+ABZ*I_ERI_D2z_S_F2xy_Px;
  Double I_ERI_D2x_Px_F2xz_Px = I_ERI_F3x_S_F2xz_Px+ABX*I_ERI_D2x_S_F2xz_Px;
  Double I_ERI_Dxy_Px_F2xz_Px = I_ERI_F2xy_S_F2xz_Px+ABX*I_ERI_Dxy_S_F2xz_Px;
  Double I_ERI_Dxz_Px_F2xz_Px = I_ERI_F2xz_S_F2xz_Px+ABX*I_ERI_Dxz_S_F2xz_Px;
  Double I_ERI_D2y_Px_F2xz_Px = I_ERI_Fx2y_S_F2xz_Px+ABX*I_ERI_D2y_S_F2xz_Px;
  Double I_ERI_Dyz_Px_F2xz_Px = I_ERI_Fxyz_S_F2xz_Px+ABX*I_ERI_Dyz_S_F2xz_Px;
  Double I_ERI_D2z_Px_F2xz_Px = I_ERI_Fx2z_S_F2xz_Px+ABX*I_ERI_D2z_S_F2xz_Px;
  Double I_ERI_D2x_Py_F2xz_Px = I_ERI_F2xy_S_F2xz_Px+ABY*I_ERI_D2x_S_F2xz_Px;
  Double I_ERI_Dxy_Py_F2xz_Px = I_ERI_Fx2y_S_F2xz_Px+ABY*I_ERI_Dxy_S_F2xz_Px;
  Double I_ERI_Dxz_Py_F2xz_Px = I_ERI_Fxyz_S_F2xz_Px+ABY*I_ERI_Dxz_S_F2xz_Px;
  Double I_ERI_D2y_Py_F2xz_Px = I_ERI_F3y_S_F2xz_Px+ABY*I_ERI_D2y_S_F2xz_Px;
  Double I_ERI_Dyz_Py_F2xz_Px = I_ERI_F2yz_S_F2xz_Px+ABY*I_ERI_Dyz_S_F2xz_Px;
  Double I_ERI_D2z_Py_F2xz_Px = I_ERI_Fy2z_S_F2xz_Px+ABY*I_ERI_D2z_S_F2xz_Px;
  Double I_ERI_D2x_Pz_F2xz_Px = I_ERI_F2xz_S_F2xz_Px+ABZ*I_ERI_D2x_S_F2xz_Px;
  Double I_ERI_Dxy_Pz_F2xz_Px = I_ERI_Fxyz_S_F2xz_Px+ABZ*I_ERI_Dxy_S_F2xz_Px;
  Double I_ERI_Dxz_Pz_F2xz_Px = I_ERI_Fx2z_S_F2xz_Px+ABZ*I_ERI_Dxz_S_F2xz_Px;
  Double I_ERI_D2y_Pz_F2xz_Px = I_ERI_F2yz_S_F2xz_Px+ABZ*I_ERI_D2y_S_F2xz_Px;
  Double I_ERI_Dyz_Pz_F2xz_Px = I_ERI_Fy2z_S_F2xz_Px+ABZ*I_ERI_Dyz_S_F2xz_Px;
  Double I_ERI_D2z_Pz_F2xz_Px = I_ERI_F3z_S_F2xz_Px+ABZ*I_ERI_D2z_S_F2xz_Px;
  Double I_ERI_D2x_Px_Fx2y_Px = I_ERI_F3x_S_Fx2y_Px+ABX*I_ERI_D2x_S_Fx2y_Px;
  Double I_ERI_Dxy_Px_Fx2y_Px = I_ERI_F2xy_S_Fx2y_Px+ABX*I_ERI_Dxy_S_Fx2y_Px;
  Double I_ERI_Dxz_Px_Fx2y_Px = I_ERI_F2xz_S_Fx2y_Px+ABX*I_ERI_Dxz_S_Fx2y_Px;
  Double I_ERI_D2y_Px_Fx2y_Px = I_ERI_Fx2y_S_Fx2y_Px+ABX*I_ERI_D2y_S_Fx2y_Px;
  Double I_ERI_Dyz_Px_Fx2y_Px = I_ERI_Fxyz_S_Fx2y_Px+ABX*I_ERI_Dyz_S_Fx2y_Px;
  Double I_ERI_D2z_Px_Fx2y_Px = I_ERI_Fx2z_S_Fx2y_Px+ABX*I_ERI_D2z_S_Fx2y_Px;
  Double I_ERI_D2x_Py_Fx2y_Px = I_ERI_F2xy_S_Fx2y_Px+ABY*I_ERI_D2x_S_Fx2y_Px;
  Double I_ERI_Dxy_Py_Fx2y_Px = I_ERI_Fx2y_S_Fx2y_Px+ABY*I_ERI_Dxy_S_Fx2y_Px;
  Double I_ERI_Dxz_Py_Fx2y_Px = I_ERI_Fxyz_S_Fx2y_Px+ABY*I_ERI_Dxz_S_Fx2y_Px;
  Double I_ERI_D2y_Py_Fx2y_Px = I_ERI_F3y_S_Fx2y_Px+ABY*I_ERI_D2y_S_Fx2y_Px;
  Double I_ERI_Dyz_Py_Fx2y_Px = I_ERI_F2yz_S_Fx2y_Px+ABY*I_ERI_Dyz_S_Fx2y_Px;
  Double I_ERI_D2z_Py_Fx2y_Px = I_ERI_Fy2z_S_Fx2y_Px+ABY*I_ERI_D2z_S_Fx2y_Px;
  Double I_ERI_D2x_Pz_Fx2y_Px = I_ERI_F2xz_S_Fx2y_Px+ABZ*I_ERI_D2x_S_Fx2y_Px;
  Double I_ERI_Dxy_Pz_Fx2y_Px = I_ERI_Fxyz_S_Fx2y_Px+ABZ*I_ERI_Dxy_S_Fx2y_Px;
  Double I_ERI_Dxz_Pz_Fx2y_Px = I_ERI_Fx2z_S_Fx2y_Px+ABZ*I_ERI_Dxz_S_Fx2y_Px;
  Double I_ERI_D2y_Pz_Fx2y_Px = I_ERI_F2yz_S_Fx2y_Px+ABZ*I_ERI_D2y_S_Fx2y_Px;
  Double I_ERI_Dyz_Pz_Fx2y_Px = I_ERI_Fy2z_S_Fx2y_Px+ABZ*I_ERI_Dyz_S_Fx2y_Px;
  Double I_ERI_D2z_Pz_Fx2y_Px = I_ERI_F3z_S_Fx2y_Px+ABZ*I_ERI_D2z_S_Fx2y_Px;
  Double I_ERI_D2x_Px_Fxyz_Px = I_ERI_F3x_S_Fxyz_Px+ABX*I_ERI_D2x_S_Fxyz_Px;
  Double I_ERI_Dxy_Px_Fxyz_Px = I_ERI_F2xy_S_Fxyz_Px+ABX*I_ERI_Dxy_S_Fxyz_Px;
  Double I_ERI_Dxz_Px_Fxyz_Px = I_ERI_F2xz_S_Fxyz_Px+ABX*I_ERI_Dxz_S_Fxyz_Px;
  Double I_ERI_D2y_Px_Fxyz_Px = I_ERI_Fx2y_S_Fxyz_Px+ABX*I_ERI_D2y_S_Fxyz_Px;
  Double I_ERI_Dyz_Px_Fxyz_Px = I_ERI_Fxyz_S_Fxyz_Px+ABX*I_ERI_Dyz_S_Fxyz_Px;
  Double I_ERI_D2z_Px_Fxyz_Px = I_ERI_Fx2z_S_Fxyz_Px+ABX*I_ERI_D2z_S_Fxyz_Px;
  Double I_ERI_D2x_Py_Fxyz_Px = I_ERI_F2xy_S_Fxyz_Px+ABY*I_ERI_D2x_S_Fxyz_Px;
  Double I_ERI_Dxy_Py_Fxyz_Px = I_ERI_Fx2y_S_Fxyz_Px+ABY*I_ERI_Dxy_S_Fxyz_Px;
  Double I_ERI_Dxz_Py_Fxyz_Px = I_ERI_Fxyz_S_Fxyz_Px+ABY*I_ERI_Dxz_S_Fxyz_Px;
  Double I_ERI_D2y_Py_Fxyz_Px = I_ERI_F3y_S_Fxyz_Px+ABY*I_ERI_D2y_S_Fxyz_Px;
  Double I_ERI_Dyz_Py_Fxyz_Px = I_ERI_F2yz_S_Fxyz_Px+ABY*I_ERI_Dyz_S_Fxyz_Px;
  Double I_ERI_D2z_Py_Fxyz_Px = I_ERI_Fy2z_S_Fxyz_Px+ABY*I_ERI_D2z_S_Fxyz_Px;
  Double I_ERI_D2x_Pz_Fxyz_Px = I_ERI_F2xz_S_Fxyz_Px+ABZ*I_ERI_D2x_S_Fxyz_Px;
  Double I_ERI_Dxy_Pz_Fxyz_Px = I_ERI_Fxyz_S_Fxyz_Px+ABZ*I_ERI_Dxy_S_Fxyz_Px;
  Double I_ERI_Dxz_Pz_Fxyz_Px = I_ERI_Fx2z_S_Fxyz_Px+ABZ*I_ERI_Dxz_S_Fxyz_Px;
  Double I_ERI_D2y_Pz_Fxyz_Px = I_ERI_F2yz_S_Fxyz_Px+ABZ*I_ERI_D2y_S_Fxyz_Px;
  Double I_ERI_Dyz_Pz_Fxyz_Px = I_ERI_Fy2z_S_Fxyz_Px+ABZ*I_ERI_Dyz_S_Fxyz_Px;
  Double I_ERI_D2z_Pz_Fxyz_Px = I_ERI_F3z_S_Fxyz_Px+ABZ*I_ERI_D2z_S_Fxyz_Px;
  Double I_ERI_D2x_Px_Fx2z_Px = I_ERI_F3x_S_Fx2z_Px+ABX*I_ERI_D2x_S_Fx2z_Px;
  Double I_ERI_Dxy_Px_Fx2z_Px = I_ERI_F2xy_S_Fx2z_Px+ABX*I_ERI_Dxy_S_Fx2z_Px;
  Double I_ERI_Dxz_Px_Fx2z_Px = I_ERI_F2xz_S_Fx2z_Px+ABX*I_ERI_Dxz_S_Fx2z_Px;
  Double I_ERI_D2y_Px_Fx2z_Px = I_ERI_Fx2y_S_Fx2z_Px+ABX*I_ERI_D2y_S_Fx2z_Px;
  Double I_ERI_Dyz_Px_Fx2z_Px = I_ERI_Fxyz_S_Fx2z_Px+ABX*I_ERI_Dyz_S_Fx2z_Px;
  Double I_ERI_D2z_Px_Fx2z_Px = I_ERI_Fx2z_S_Fx2z_Px+ABX*I_ERI_D2z_S_Fx2z_Px;
  Double I_ERI_D2x_Py_Fx2z_Px = I_ERI_F2xy_S_Fx2z_Px+ABY*I_ERI_D2x_S_Fx2z_Px;
  Double I_ERI_Dxy_Py_Fx2z_Px = I_ERI_Fx2y_S_Fx2z_Px+ABY*I_ERI_Dxy_S_Fx2z_Px;
  Double I_ERI_Dxz_Py_Fx2z_Px = I_ERI_Fxyz_S_Fx2z_Px+ABY*I_ERI_Dxz_S_Fx2z_Px;
  Double I_ERI_D2y_Py_Fx2z_Px = I_ERI_F3y_S_Fx2z_Px+ABY*I_ERI_D2y_S_Fx2z_Px;
  Double I_ERI_Dyz_Py_Fx2z_Px = I_ERI_F2yz_S_Fx2z_Px+ABY*I_ERI_Dyz_S_Fx2z_Px;
  Double I_ERI_D2z_Py_Fx2z_Px = I_ERI_Fy2z_S_Fx2z_Px+ABY*I_ERI_D2z_S_Fx2z_Px;
  Double I_ERI_D2x_Pz_Fx2z_Px = I_ERI_F2xz_S_Fx2z_Px+ABZ*I_ERI_D2x_S_Fx2z_Px;
  Double I_ERI_Dxy_Pz_Fx2z_Px = I_ERI_Fxyz_S_Fx2z_Px+ABZ*I_ERI_Dxy_S_Fx2z_Px;
  Double I_ERI_Dxz_Pz_Fx2z_Px = I_ERI_Fx2z_S_Fx2z_Px+ABZ*I_ERI_Dxz_S_Fx2z_Px;
  Double I_ERI_D2y_Pz_Fx2z_Px = I_ERI_F2yz_S_Fx2z_Px+ABZ*I_ERI_D2y_S_Fx2z_Px;
  Double I_ERI_Dyz_Pz_Fx2z_Px = I_ERI_Fy2z_S_Fx2z_Px+ABZ*I_ERI_Dyz_S_Fx2z_Px;
  Double I_ERI_D2z_Pz_Fx2z_Px = I_ERI_F3z_S_Fx2z_Px+ABZ*I_ERI_D2z_S_Fx2z_Px;
  Double I_ERI_D2x_Px_F3y_Px = I_ERI_F3x_S_F3y_Px+ABX*I_ERI_D2x_S_F3y_Px;
  Double I_ERI_Dxy_Px_F3y_Px = I_ERI_F2xy_S_F3y_Px+ABX*I_ERI_Dxy_S_F3y_Px;
  Double I_ERI_Dxz_Px_F3y_Px = I_ERI_F2xz_S_F3y_Px+ABX*I_ERI_Dxz_S_F3y_Px;
  Double I_ERI_D2y_Px_F3y_Px = I_ERI_Fx2y_S_F3y_Px+ABX*I_ERI_D2y_S_F3y_Px;
  Double I_ERI_Dyz_Px_F3y_Px = I_ERI_Fxyz_S_F3y_Px+ABX*I_ERI_Dyz_S_F3y_Px;
  Double I_ERI_D2z_Px_F3y_Px = I_ERI_Fx2z_S_F3y_Px+ABX*I_ERI_D2z_S_F3y_Px;
  Double I_ERI_D2x_Py_F3y_Px = I_ERI_F2xy_S_F3y_Px+ABY*I_ERI_D2x_S_F3y_Px;
  Double I_ERI_Dxy_Py_F3y_Px = I_ERI_Fx2y_S_F3y_Px+ABY*I_ERI_Dxy_S_F3y_Px;
  Double I_ERI_Dxz_Py_F3y_Px = I_ERI_Fxyz_S_F3y_Px+ABY*I_ERI_Dxz_S_F3y_Px;
  Double I_ERI_D2y_Py_F3y_Px = I_ERI_F3y_S_F3y_Px+ABY*I_ERI_D2y_S_F3y_Px;
  Double I_ERI_Dyz_Py_F3y_Px = I_ERI_F2yz_S_F3y_Px+ABY*I_ERI_Dyz_S_F3y_Px;
  Double I_ERI_D2z_Py_F3y_Px = I_ERI_Fy2z_S_F3y_Px+ABY*I_ERI_D2z_S_F3y_Px;
  Double I_ERI_D2x_Pz_F3y_Px = I_ERI_F2xz_S_F3y_Px+ABZ*I_ERI_D2x_S_F3y_Px;
  Double I_ERI_Dxy_Pz_F3y_Px = I_ERI_Fxyz_S_F3y_Px+ABZ*I_ERI_Dxy_S_F3y_Px;
  Double I_ERI_Dxz_Pz_F3y_Px = I_ERI_Fx2z_S_F3y_Px+ABZ*I_ERI_Dxz_S_F3y_Px;
  Double I_ERI_D2y_Pz_F3y_Px = I_ERI_F2yz_S_F3y_Px+ABZ*I_ERI_D2y_S_F3y_Px;
  Double I_ERI_Dyz_Pz_F3y_Px = I_ERI_Fy2z_S_F3y_Px+ABZ*I_ERI_Dyz_S_F3y_Px;
  Double I_ERI_D2z_Pz_F3y_Px = I_ERI_F3z_S_F3y_Px+ABZ*I_ERI_D2z_S_F3y_Px;
  Double I_ERI_D2x_Px_F2yz_Px = I_ERI_F3x_S_F2yz_Px+ABX*I_ERI_D2x_S_F2yz_Px;
  Double I_ERI_Dxy_Px_F2yz_Px = I_ERI_F2xy_S_F2yz_Px+ABX*I_ERI_Dxy_S_F2yz_Px;
  Double I_ERI_Dxz_Px_F2yz_Px = I_ERI_F2xz_S_F2yz_Px+ABX*I_ERI_Dxz_S_F2yz_Px;
  Double I_ERI_D2y_Px_F2yz_Px = I_ERI_Fx2y_S_F2yz_Px+ABX*I_ERI_D2y_S_F2yz_Px;
  Double I_ERI_Dyz_Px_F2yz_Px = I_ERI_Fxyz_S_F2yz_Px+ABX*I_ERI_Dyz_S_F2yz_Px;
  Double I_ERI_D2z_Px_F2yz_Px = I_ERI_Fx2z_S_F2yz_Px+ABX*I_ERI_D2z_S_F2yz_Px;
  Double I_ERI_D2x_Py_F2yz_Px = I_ERI_F2xy_S_F2yz_Px+ABY*I_ERI_D2x_S_F2yz_Px;
  Double I_ERI_Dxy_Py_F2yz_Px = I_ERI_Fx2y_S_F2yz_Px+ABY*I_ERI_Dxy_S_F2yz_Px;
  Double I_ERI_Dxz_Py_F2yz_Px = I_ERI_Fxyz_S_F2yz_Px+ABY*I_ERI_Dxz_S_F2yz_Px;
  Double I_ERI_D2y_Py_F2yz_Px = I_ERI_F3y_S_F2yz_Px+ABY*I_ERI_D2y_S_F2yz_Px;
  Double I_ERI_Dyz_Py_F2yz_Px = I_ERI_F2yz_S_F2yz_Px+ABY*I_ERI_Dyz_S_F2yz_Px;
  Double I_ERI_D2z_Py_F2yz_Px = I_ERI_Fy2z_S_F2yz_Px+ABY*I_ERI_D2z_S_F2yz_Px;
  Double I_ERI_D2x_Pz_F2yz_Px = I_ERI_F2xz_S_F2yz_Px+ABZ*I_ERI_D2x_S_F2yz_Px;
  Double I_ERI_Dxy_Pz_F2yz_Px = I_ERI_Fxyz_S_F2yz_Px+ABZ*I_ERI_Dxy_S_F2yz_Px;
  Double I_ERI_Dxz_Pz_F2yz_Px = I_ERI_Fx2z_S_F2yz_Px+ABZ*I_ERI_Dxz_S_F2yz_Px;
  Double I_ERI_D2y_Pz_F2yz_Px = I_ERI_F2yz_S_F2yz_Px+ABZ*I_ERI_D2y_S_F2yz_Px;
  Double I_ERI_Dyz_Pz_F2yz_Px = I_ERI_Fy2z_S_F2yz_Px+ABZ*I_ERI_Dyz_S_F2yz_Px;
  Double I_ERI_D2z_Pz_F2yz_Px = I_ERI_F3z_S_F2yz_Px+ABZ*I_ERI_D2z_S_F2yz_Px;
  Double I_ERI_D2x_Px_Fy2z_Px = I_ERI_F3x_S_Fy2z_Px+ABX*I_ERI_D2x_S_Fy2z_Px;
  Double I_ERI_Dxy_Px_Fy2z_Px = I_ERI_F2xy_S_Fy2z_Px+ABX*I_ERI_Dxy_S_Fy2z_Px;
  Double I_ERI_Dxz_Px_Fy2z_Px = I_ERI_F2xz_S_Fy2z_Px+ABX*I_ERI_Dxz_S_Fy2z_Px;
  Double I_ERI_D2y_Px_Fy2z_Px = I_ERI_Fx2y_S_Fy2z_Px+ABX*I_ERI_D2y_S_Fy2z_Px;
  Double I_ERI_Dyz_Px_Fy2z_Px = I_ERI_Fxyz_S_Fy2z_Px+ABX*I_ERI_Dyz_S_Fy2z_Px;
  Double I_ERI_D2z_Px_Fy2z_Px = I_ERI_Fx2z_S_Fy2z_Px+ABX*I_ERI_D2z_S_Fy2z_Px;
  Double I_ERI_D2x_Py_Fy2z_Px = I_ERI_F2xy_S_Fy2z_Px+ABY*I_ERI_D2x_S_Fy2z_Px;
  Double I_ERI_Dxy_Py_Fy2z_Px = I_ERI_Fx2y_S_Fy2z_Px+ABY*I_ERI_Dxy_S_Fy2z_Px;
  Double I_ERI_Dxz_Py_Fy2z_Px = I_ERI_Fxyz_S_Fy2z_Px+ABY*I_ERI_Dxz_S_Fy2z_Px;
  Double I_ERI_D2y_Py_Fy2z_Px = I_ERI_F3y_S_Fy2z_Px+ABY*I_ERI_D2y_S_Fy2z_Px;
  Double I_ERI_Dyz_Py_Fy2z_Px = I_ERI_F2yz_S_Fy2z_Px+ABY*I_ERI_Dyz_S_Fy2z_Px;
  Double I_ERI_D2z_Py_Fy2z_Px = I_ERI_Fy2z_S_Fy2z_Px+ABY*I_ERI_D2z_S_Fy2z_Px;
  Double I_ERI_D2x_Pz_Fy2z_Px = I_ERI_F2xz_S_Fy2z_Px+ABZ*I_ERI_D2x_S_Fy2z_Px;
  Double I_ERI_Dxy_Pz_Fy2z_Px = I_ERI_Fxyz_S_Fy2z_Px+ABZ*I_ERI_Dxy_S_Fy2z_Px;
  Double I_ERI_Dxz_Pz_Fy2z_Px = I_ERI_Fx2z_S_Fy2z_Px+ABZ*I_ERI_Dxz_S_Fy2z_Px;
  Double I_ERI_D2y_Pz_Fy2z_Px = I_ERI_F2yz_S_Fy2z_Px+ABZ*I_ERI_D2y_S_Fy2z_Px;
  Double I_ERI_Dyz_Pz_Fy2z_Px = I_ERI_Fy2z_S_Fy2z_Px+ABZ*I_ERI_Dyz_S_Fy2z_Px;
  Double I_ERI_D2z_Pz_Fy2z_Px = I_ERI_F3z_S_Fy2z_Px+ABZ*I_ERI_D2z_S_Fy2z_Px;
  Double I_ERI_D2x_Px_F3z_Px = I_ERI_F3x_S_F3z_Px+ABX*I_ERI_D2x_S_F3z_Px;
  Double I_ERI_Dxy_Px_F3z_Px = I_ERI_F2xy_S_F3z_Px+ABX*I_ERI_Dxy_S_F3z_Px;
  Double I_ERI_Dxz_Px_F3z_Px = I_ERI_F2xz_S_F3z_Px+ABX*I_ERI_Dxz_S_F3z_Px;
  Double I_ERI_D2y_Px_F3z_Px = I_ERI_Fx2y_S_F3z_Px+ABX*I_ERI_D2y_S_F3z_Px;
  Double I_ERI_Dyz_Px_F3z_Px = I_ERI_Fxyz_S_F3z_Px+ABX*I_ERI_Dyz_S_F3z_Px;
  Double I_ERI_D2z_Px_F3z_Px = I_ERI_Fx2z_S_F3z_Px+ABX*I_ERI_D2z_S_F3z_Px;
  Double I_ERI_D2x_Py_F3z_Px = I_ERI_F2xy_S_F3z_Px+ABY*I_ERI_D2x_S_F3z_Px;
  Double I_ERI_Dxy_Py_F3z_Px = I_ERI_Fx2y_S_F3z_Px+ABY*I_ERI_Dxy_S_F3z_Px;
  Double I_ERI_Dxz_Py_F3z_Px = I_ERI_Fxyz_S_F3z_Px+ABY*I_ERI_Dxz_S_F3z_Px;
  Double I_ERI_D2y_Py_F3z_Px = I_ERI_F3y_S_F3z_Px+ABY*I_ERI_D2y_S_F3z_Px;
  Double I_ERI_Dyz_Py_F3z_Px = I_ERI_F2yz_S_F3z_Px+ABY*I_ERI_Dyz_S_F3z_Px;
  Double I_ERI_D2z_Py_F3z_Px = I_ERI_Fy2z_S_F3z_Px+ABY*I_ERI_D2z_S_F3z_Px;
  Double I_ERI_D2x_Pz_F3z_Px = I_ERI_F2xz_S_F3z_Px+ABZ*I_ERI_D2x_S_F3z_Px;
  Double I_ERI_Dxy_Pz_F3z_Px = I_ERI_Fxyz_S_F3z_Px+ABZ*I_ERI_Dxy_S_F3z_Px;
  Double I_ERI_Dxz_Pz_F3z_Px = I_ERI_Fx2z_S_F3z_Px+ABZ*I_ERI_Dxz_S_F3z_Px;
  Double I_ERI_D2y_Pz_F3z_Px = I_ERI_F2yz_S_F3z_Px+ABZ*I_ERI_D2y_S_F3z_Px;
  Double I_ERI_Dyz_Pz_F3z_Px = I_ERI_Fy2z_S_F3z_Px+ABZ*I_ERI_Dyz_S_F3z_Px;
  Double I_ERI_D2z_Pz_F3z_Px = I_ERI_F3z_S_F3z_Px+ABZ*I_ERI_D2z_S_F3z_Px;
  Double I_ERI_D2x_Px_F3x_Py = I_ERI_F3x_S_F3x_Py+ABX*I_ERI_D2x_S_F3x_Py;
  Double I_ERI_Dxy_Px_F3x_Py = I_ERI_F2xy_S_F3x_Py+ABX*I_ERI_Dxy_S_F3x_Py;
  Double I_ERI_Dxz_Px_F3x_Py = I_ERI_F2xz_S_F3x_Py+ABX*I_ERI_Dxz_S_F3x_Py;
  Double I_ERI_D2y_Px_F3x_Py = I_ERI_Fx2y_S_F3x_Py+ABX*I_ERI_D2y_S_F3x_Py;
  Double I_ERI_Dyz_Px_F3x_Py = I_ERI_Fxyz_S_F3x_Py+ABX*I_ERI_Dyz_S_F3x_Py;
  Double I_ERI_D2z_Px_F3x_Py = I_ERI_Fx2z_S_F3x_Py+ABX*I_ERI_D2z_S_F3x_Py;
  Double I_ERI_D2x_Py_F3x_Py = I_ERI_F2xy_S_F3x_Py+ABY*I_ERI_D2x_S_F3x_Py;
  Double I_ERI_Dxy_Py_F3x_Py = I_ERI_Fx2y_S_F3x_Py+ABY*I_ERI_Dxy_S_F3x_Py;
  Double I_ERI_Dxz_Py_F3x_Py = I_ERI_Fxyz_S_F3x_Py+ABY*I_ERI_Dxz_S_F3x_Py;
  Double I_ERI_D2y_Py_F3x_Py = I_ERI_F3y_S_F3x_Py+ABY*I_ERI_D2y_S_F3x_Py;
  Double I_ERI_Dyz_Py_F3x_Py = I_ERI_F2yz_S_F3x_Py+ABY*I_ERI_Dyz_S_F3x_Py;
  Double I_ERI_D2z_Py_F3x_Py = I_ERI_Fy2z_S_F3x_Py+ABY*I_ERI_D2z_S_F3x_Py;
  Double I_ERI_D2x_Pz_F3x_Py = I_ERI_F2xz_S_F3x_Py+ABZ*I_ERI_D2x_S_F3x_Py;
  Double I_ERI_Dxy_Pz_F3x_Py = I_ERI_Fxyz_S_F3x_Py+ABZ*I_ERI_Dxy_S_F3x_Py;
  Double I_ERI_Dxz_Pz_F3x_Py = I_ERI_Fx2z_S_F3x_Py+ABZ*I_ERI_Dxz_S_F3x_Py;
  Double I_ERI_D2y_Pz_F3x_Py = I_ERI_F2yz_S_F3x_Py+ABZ*I_ERI_D2y_S_F3x_Py;
  Double I_ERI_Dyz_Pz_F3x_Py = I_ERI_Fy2z_S_F3x_Py+ABZ*I_ERI_Dyz_S_F3x_Py;
  Double I_ERI_D2z_Pz_F3x_Py = I_ERI_F3z_S_F3x_Py+ABZ*I_ERI_D2z_S_F3x_Py;
  Double I_ERI_D2x_Px_F2xy_Py = I_ERI_F3x_S_F2xy_Py+ABX*I_ERI_D2x_S_F2xy_Py;
  Double I_ERI_Dxy_Px_F2xy_Py = I_ERI_F2xy_S_F2xy_Py+ABX*I_ERI_Dxy_S_F2xy_Py;
  Double I_ERI_Dxz_Px_F2xy_Py = I_ERI_F2xz_S_F2xy_Py+ABX*I_ERI_Dxz_S_F2xy_Py;
  Double I_ERI_D2y_Px_F2xy_Py = I_ERI_Fx2y_S_F2xy_Py+ABX*I_ERI_D2y_S_F2xy_Py;
  Double I_ERI_Dyz_Px_F2xy_Py = I_ERI_Fxyz_S_F2xy_Py+ABX*I_ERI_Dyz_S_F2xy_Py;
  Double I_ERI_D2z_Px_F2xy_Py = I_ERI_Fx2z_S_F2xy_Py+ABX*I_ERI_D2z_S_F2xy_Py;
  Double I_ERI_D2x_Py_F2xy_Py = I_ERI_F2xy_S_F2xy_Py+ABY*I_ERI_D2x_S_F2xy_Py;
  Double I_ERI_Dxy_Py_F2xy_Py = I_ERI_Fx2y_S_F2xy_Py+ABY*I_ERI_Dxy_S_F2xy_Py;
  Double I_ERI_Dxz_Py_F2xy_Py = I_ERI_Fxyz_S_F2xy_Py+ABY*I_ERI_Dxz_S_F2xy_Py;
  Double I_ERI_D2y_Py_F2xy_Py = I_ERI_F3y_S_F2xy_Py+ABY*I_ERI_D2y_S_F2xy_Py;
  Double I_ERI_Dyz_Py_F2xy_Py = I_ERI_F2yz_S_F2xy_Py+ABY*I_ERI_Dyz_S_F2xy_Py;
  Double I_ERI_D2z_Py_F2xy_Py = I_ERI_Fy2z_S_F2xy_Py+ABY*I_ERI_D2z_S_F2xy_Py;
  Double I_ERI_D2x_Pz_F2xy_Py = I_ERI_F2xz_S_F2xy_Py+ABZ*I_ERI_D2x_S_F2xy_Py;
  Double I_ERI_Dxy_Pz_F2xy_Py = I_ERI_Fxyz_S_F2xy_Py+ABZ*I_ERI_Dxy_S_F2xy_Py;
  Double I_ERI_Dxz_Pz_F2xy_Py = I_ERI_Fx2z_S_F2xy_Py+ABZ*I_ERI_Dxz_S_F2xy_Py;
  Double I_ERI_D2y_Pz_F2xy_Py = I_ERI_F2yz_S_F2xy_Py+ABZ*I_ERI_D2y_S_F2xy_Py;
  Double I_ERI_Dyz_Pz_F2xy_Py = I_ERI_Fy2z_S_F2xy_Py+ABZ*I_ERI_Dyz_S_F2xy_Py;
  Double I_ERI_D2z_Pz_F2xy_Py = I_ERI_F3z_S_F2xy_Py+ABZ*I_ERI_D2z_S_F2xy_Py;
  Double I_ERI_D2x_Px_F2xz_Py = I_ERI_F3x_S_F2xz_Py+ABX*I_ERI_D2x_S_F2xz_Py;
  Double I_ERI_Dxy_Px_F2xz_Py = I_ERI_F2xy_S_F2xz_Py+ABX*I_ERI_Dxy_S_F2xz_Py;
  Double I_ERI_Dxz_Px_F2xz_Py = I_ERI_F2xz_S_F2xz_Py+ABX*I_ERI_Dxz_S_F2xz_Py;
  Double I_ERI_D2y_Px_F2xz_Py = I_ERI_Fx2y_S_F2xz_Py+ABX*I_ERI_D2y_S_F2xz_Py;
  Double I_ERI_Dyz_Px_F2xz_Py = I_ERI_Fxyz_S_F2xz_Py+ABX*I_ERI_Dyz_S_F2xz_Py;
  Double I_ERI_D2z_Px_F2xz_Py = I_ERI_Fx2z_S_F2xz_Py+ABX*I_ERI_D2z_S_F2xz_Py;
  Double I_ERI_D2x_Py_F2xz_Py = I_ERI_F2xy_S_F2xz_Py+ABY*I_ERI_D2x_S_F2xz_Py;
  Double I_ERI_Dxy_Py_F2xz_Py = I_ERI_Fx2y_S_F2xz_Py+ABY*I_ERI_Dxy_S_F2xz_Py;
  Double I_ERI_Dxz_Py_F2xz_Py = I_ERI_Fxyz_S_F2xz_Py+ABY*I_ERI_Dxz_S_F2xz_Py;
  Double I_ERI_D2y_Py_F2xz_Py = I_ERI_F3y_S_F2xz_Py+ABY*I_ERI_D2y_S_F2xz_Py;
  Double I_ERI_Dyz_Py_F2xz_Py = I_ERI_F2yz_S_F2xz_Py+ABY*I_ERI_Dyz_S_F2xz_Py;
  Double I_ERI_D2z_Py_F2xz_Py = I_ERI_Fy2z_S_F2xz_Py+ABY*I_ERI_D2z_S_F2xz_Py;
  Double I_ERI_D2x_Pz_F2xz_Py = I_ERI_F2xz_S_F2xz_Py+ABZ*I_ERI_D2x_S_F2xz_Py;
  Double I_ERI_Dxy_Pz_F2xz_Py = I_ERI_Fxyz_S_F2xz_Py+ABZ*I_ERI_Dxy_S_F2xz_Py;
  Double I_ERI_Dxz_Pz_F2xz_Py = I_ERI_Fx2z_S_F2xz_Py+ABZ*I_ERI_Dxz_S_F2xz_Py;
  Double I_ERI_D2y_Pz_F2xz_Py = I_ERI_F2yz_S_F2xz_Py+ABZ*I_ERI_D2y_S_F2xz_Py;
  Double I_ERI_Dyz_Pz_F2xz_Py = I_ERI_Fy2z_S_F2xz_Py+ABZ*I_ERI_Dyz_S_F2xz_Py;
  Double I_ERI_D2z_Pz_F2xz_Py = I_ERI_F3z_S_F2xz_Py+ABZ*I_ERI_D2z_S_F2xz_Py;
  Double I_ERI_D2x_Px_Fx2y_Py = I_ERI_F3x_S_Fx2y_Py+ABX*I_ERI_D2x_S_Fx2y_Py;
  Double I_ERI_Dxy_Px_Fx2y_Py = I_ERI_F2xy_S_Fx2y_Py+ABX*I_ERI_Dxy_S_Fx2y_Py;
  Double I_ERI_Dxz_Px_Fx2y_Py = I_ERI_F2xz_S_Fx2y_Py+ABX*I_ERI_Dxz_S_Fx2y_Py;
  Double I_ERI_D2y_Px_Fx2y_Py = I_ERI_Fx2y_S_Fx2y_Py+ABX*I_ERI_D2y_S_Fx2y_Py;
  Double I_ERI_Dyz_Px_Fx2y_Py = I_ERI_Fxyz_S_Fx2y_Py+ABX*I_ERI_Dyz_S_Fx2y_Py;
  Double I_ERI_D2z_Px_Fx2y_Py = I_ERI_Fx2z_S_Fx2y_Py+ABX*I_ERI_D2z_S_Fx2y_Py;
  Double I_ERI_D2x_Py_Fx2y_Py = I_ERI_F2xy_S_Fx2y_Py+ABY*I_ERI_D2x_S_Fx2y_Py;
  Double I_ERI_Dxy_Py_Fx2y_Py = I_ERI_Fx2y_S_Fx2y_Py+ABY*I_ERI_Dxy_S_Fx2y_Py;
  Double I_ERI_Dxz_Py_Fx2y_Py = I_ERI_Fxyz_S_Fx2y_Py+ABY*I_ERI_Dxz_S_Fx2y_Py;
  Double I_ERI_D2y_Py_Fx2y_Py = I_ERI_F3y_S_Fx2y_Py+ABY*I_ERI_D2y_S_Fx2y_Py;
  Double I_ERI_Dyz_Py_Fx2y_Py = I_ERI_F2yz_S_Fx2y_Py+ABY*I_ERI_Dyz_S_Fx2y_Py;
  Double I_ERI_D2z_Py_Fx2y_Py = I_ERI_Fy2z_S_Fx2y_Py+ABY*I_ERI_D2z_S_Fx2y_Py;
  Double I_ERI_D2x_Pz_Fx2y_Py = I_ERI_F2xz_S_Fx2y_Py+ABZ*I_ERI_D2x_S_Fx2y_Py;
  Double I_ERI_Dxy_Pz_Fx2y_Py = I_ERI_Fxyz_S_Fx2y_Py+ABZ*I_ERI_Dxy_S_Fx2y_Py;
  Double I_ERI_Dxz_Pz_Fx2y_Py = I_ERI_Fx2z_S_Fx2y_Py+ABZ*I_ERI_Dxz_S_Fx2y_Py;
  Double I_ERI_D2y_Pz_Fx2y_Py = I_ERI_F2yz_S_Fx2y_Py+ABZ*I_ERI_D2y_S_Fx2y_Py;
  Double I_ERI_Dyz_Pz_Fx2y_Py = I_ERI_Fy2z_S_Fx2y_Py+ABZ*I_ERI_Dyz_S_Fx2y_Py;
  Double I_ERI_D2z_Pz_Fx2y_Py = I_ERI_F3z_S_Fx2y_Py+ABZ*I_ERI_D2z_S_Fx2y_Py;
  Double I_ERI_D2x_Px_Fxyz_Py = I_ERI_F3x_S_Fxyz_Py+ABX*I_ERI_D2x_S_Fxyz_Py;
  Double I_ERI_Dxy_Px_Fxyz_Py = I_ERI_F2xy_S_Fxyz_Py+ABX*I_ERI_Dxy_S_Fxyz_Py;
  Double I_ERI_Dxz_Px_Fxyz_Py = I_ERI_F2xz_S_Fxyz_Py+ABX*I_ERI_Dxz_S_Fxyz_Py;
  Double I_ERI_D2y_Px_Fxyz_Py = I_ERI_Fx2y_S_Fxyz_Py+ABX*I_ERI_D2y_S_Fxyz_Py;
  Double I_ERI_Dyz_Px_Fxyz_Py = I_ERI_Fxyz_S_Fxyz_Py+ABX*I_ERI_Dyz_S_Fxyz_Py;
  Double I_ERI_D2z_Px_Fxyz_Py = I_ERI_Fx2z_S_Fxyz_Py+ABX*I_ERI_D2z_S_Fxyz_Py;
  Double I_ERI_D2x_Py_Fxyz_Py = I_ERI_F2xy_S_Fxyz_Py+ABY*I_ERI_D2x_S_Fxyz_Py;
  Double I_ERI_Dxy_Py_Fxyz_Py = I_ERI_Fx2y_S_Fxyz_Py+ABY*I_ERI_Dxy_S_Fxyz_Py;
  Double I_ERI_Dxz_Py_Fxyz_Py = I_ERI_Fxyz_S_Fxyz_Py+ABY*I_ERI_Dxz_S_Fxyz_Py;
  Double I_ERI_D2y_Py_Fxyz_Py = I_ERI_F3y_S_Fxyz_Py+ABY*I_ERI_D2y_S_Fxyz_Py;
  Double I_ERI_Dyz_Py_Fxyz_Py = I_ERI_F2yz_S_Fxyz_Py+ABY*I_ERI_Dyz_S_Fxyz_Py;
  Double I_ERI_D2z_Py_Fxyz_Py = I_ERI_Fy2z_S_Fxyz_Py+ABY*I_ERI_D2z_S_Fxyz_Py;
  Double I_ERI_D2x_Pz_Fxyz_Py = I_ERI_F2xz_S_Fxyz_Py+ABZ*I_ERI_D2x_S_Fxyz_Py;
  Double I_ERI_Dxy_Pz_Fxyz_Py = I_ERI_Fxyz_S_Fxyz_Py+ABZ*I_ERI_Dxy_S_Fxyz_Py;
  Double I_ERI_Dxz_Pz_Fxyz_Py = I_ERI_Fx2z_S_Fxyz_Py+ABZ*I_ERI_Dxz_S_Fxyz_Py;
  Double I_ERI_D2y_Pz_Fxyz_Py = I_ERI_F2yz_S_Fxyz_Py+ABZ*I_ERI_D2y_S_Fxyz_Py;
  Double I_ERI_Dyz_Pz_Fxyz_Py = I_ERI_Fy2z_S_Fxyz_Py+ABZ*I_ERI_Dyz_S_Fxyz_Py;
  Double I_ERI_D2z_Pz_Fxyz_Py = I_ERI_F3z_S_Fxyz_Py+ABZ*I_ERI_D2z_S_Fxyz_Py;
  Double I_ERI_D2x_Px_Fx2z_Py = I_ERI_F3x_S_Fx2z_Py+ABX*I_ERI_D2x_S_Fx2z_Py;
  Double I_ERI_Dxy_Px_Fx2z_Py = I_ERI_F2xy_S_Fx2z_Py+ABX*I_ERI_Dxy_S_Fx2z_Py;
  Double I_ERI_Dxz_Px_Fx2z_Py = I_ERI_F2xz_S_Fx2z_Py+ABX*I_ERI_Dxz_S_Fx2z_Py;
  Double I_ERI_D2y_Px_Fx2z_Py = I_ERI_Fx2y_S_Fx2z_Py+ABX*I_ERI_D2y_S_Fx2z_Py;
  Double I_ERI_Dyz_Px_Fx2z_Py = I_ERI_Fxyz_S_Fx2z_Py+ABX*I_ERI_Dyz_S_Fx2z_Py;
  Double I_ERI_D2z_Px_Fx2z_Py = I_ERI_Fx2z_S_Fx2z_Py+ABX*I_ERI_D2z_S_Fx2z_Py;
  Double I_ERI_D2x_Py_Fx2z_Py = I_ERI_F2xy_S_Fx2z_Py+ABY*I_ERI_D2x_S_Fx2z_Py;
  Double I_ERI_Dxy_Py_Fx2z_Py = I_ERI_Fx2y_S_Fx2z_Py+ABY*I_ERI_Dxy_S_Fx2z_Py;
  Double I_ERI_Dxz_Py_Fx2z_Py = I_ERI_Fxyz_S_Fx2z_Py+ABY*I_ERI_Dxz_S_Fx2z_Py;
  Double I_ERI_D2y_Py_Fx2z_Py = I_ERI_F3y_S_Fx2z_Py+ABY*I_ERI_D2y_S_Fx2z_Py;
  Double I_ERI_Dyz_Py_Fx2z_Py = I_ERI_F2yz_S_Fx2z_Py+ABY*I_ERI_Dyz_S_Fx2z_Py;
  Double I_ERI_D2z_Py_Fx2z_Py = I_ERI_Fy2z_S_Fx2z_Py+ABY*I_ERI_D2z_S_Fx2z_Py;
  Double I_ERI_D2x_Pz_Fx2z_Py = I_ERI_F2xz_S_Fx2z_Py+ABZ*I_ERI_D2x_S_Fx2z_Py;
  Double I_ERI_Dxy_Pz_Fx2z_Py = I_ERI_Fxyz_S_Fx2z_Py+ABZ*I_ERI_Dxy_S_Fx2z_Py;
  Double I_ERI_Dxz_Pz_Fx2z_Py = I_ERI_Fx2z_S_Fx2z_Py+ABZ*I_ERI_Dxz_S_Fx2z_Py;
  Double I_ERI_D2y_Pz_Fx2z_Py = I_ERI_F2yz_S_Fx2z_Py+ABZ*I_ERI_D2y_S_Fx2z_Py;
  Double I_ERI_Dyz_Pz_Fx2z_Py = I_ERI_Fy2z_S_Fx2z_Py+ABZ*I_ERI_Dyz_S_Fx2z_Py;
  Double I_ERI_D2z_Pz_Fx2z_Py = I_ERI_F3z_S_Fx2z_Py+ABZ*I_ERI_D2z_S_Fx2z_Py;
  Double I_ERI_D2x_Px_F3y_Py = I_ERI_F3x_S_F3y_Py+ABX*I_ERI_D2x_S_F3y_Py;
  Double I_ERI_Dxy_Px_F3y_Py = I_ERI_F2xy_S_F3y_Py+ABX*I_ERI_Dxy_S_F3y_Py;
  Double I_ERI_Dxz_Px_F3y_Py = I_ERI_F2xz_S_F3y_Py+ABX*I_ERI_Dxz_S_F3y_Py;
  Double I_ERI_D2y_Px_F3y_Py = I_ERI_Fx2y_S_F3y_Py+ABX*I_ERI_D2y_S_F3y_Py;
  Double I_ERI_Dyz_Px_F3y_Py = I_ERI_Fxyz_S_F3y_Py+ABX*I_ERI_Dyz_S_F3y_Py;
  Double I_ERI_D2z_Px_F3y_Py = I_ERI_Fx2z_S_F3y_Py+ABX*I_ERI_D2z_S_F3y_Py;
  Double I_ERI_D2x_Py_F3y_Py = I_ERI_F2xy_S_F3y_Py+ABY*I_ERI_D2x_S_F3y_Py;
  Double I_ERI_Dxy_Py_F3y_Py = I_ERI_Fx2y_S_F3y_Py+ABY*I_ERI_Dxy_S_F3y_Py;
  Double I_ERI_Dxz_Py_F3y_Py = I_ERI_Fxyz_S_F3y_Py+ABY*I_ERI_Dxz_S_F3y_Py;
  Double I_ERI_D2y_Py_F3y_Py = I_ERI_F3y_S_F3y_Py+ABY*I_ERI_D2y_S_F3y_Py;
  Double I_ERI_Dyz_Py_F3y_Py = I_ERI_F2yz_S_F3y_Py+ABY*I_ERI_Dyz_S_F3y_Py;
  Double I_ERI_D2z_Py_F3y_Py = I_ERI_Fy2z_S_F3y_Py+ABY*I_ERI_D2z_S_F3y_Py;
  Double I_ERI_D2x_Pz_F3y_Py = I_ERI_F2xz_S_F3y_Py+ABZ*I_ERI_D2x_S_F3y_Py;
  Double I_ERI_Dxy_Pz_F3y_Py = I_ERI_Fxyz_S_F3y_Py+ABZ*I_ERI_Dxy_S_F3y_Py;
  Double I_ERI_Dxz_Pz_F3y_Py = I_ERI_Fx2z_S_F3y_Py+ABZ*I_ERI_Dxz_S_F3y_Py;
  Double I_ERI_D2y_Pz_F3y_Py = I_ERI_F2yz_S_F3y_Py+ABZ*I_ERI_D2y_S_F3y_Py;
  Double I_ERI_Dyz_Pz_F3y_Py = I_ERI_Fy2z_S_F3y_Py+ABZ*I_ERI_Dyz_S_F3y_Py;
  Double I_ERI_D2z_Pz_F3y_Py = I_ERI_F3z_S_F3y_Py+ABZ*I_ERI_D2z_S_F3y_Py;
  Double I_ERI_D2x_Px_F2yz_Py = I_ERI_F3x_S_F2yz_Py+ABX*I_ERI_D2x_S_F2yz_Py;
  Double I_ERI_Dxy_Px_F2yz_Py = I_ERI_F2xy_S_F2yz_Py+ABX*I_ERI_Dxy_S_F2yz_Py;
  Double I_ERI_Dxz_Px_F2yz_Py = I_ERI_F2xz_S_F2yz_Py+ABX*I_ERI_Dxz_S_F2yz_Py;
  Double I_ERI_D2y_Px_F2yz_Py = I_ERI_Fx2y_S_F2yz_Py+ABX*I_ERI_D2y_S_F2yz_Py;
  Double I_ERI_Dyz_Px_F2yz_Py = I_ERI_Fxyz_S_F2yz_Py+ABX*I_ERI_Dyz_S_F2yz_Py;
  Double I_ERI_D2z_Px_F2yz_Py = I_ERI_Fx2z_S_F2yz_Py+ABX*I_ERI_D2z_S_F2yz_Py;
  Double I_ERI_D2x_Py_F2yz_Py = I_ERI_F2xy_S_F2yz_Py+ABY*I_ERI_D2x_S_F2yz_Py;
  Double I_ERI_Dxy_Py_F2yz_Py = I_ERI_Fx2y_S_F2yz_Py+ABY*I_ERI_Dxy_S_F2yz_Py;
  Double I_ERI_Dxz_Py_F2yz_Py = I_ERI_Fxyz_S_F2yz_Py+ABY*I_ERI_Dxz_S_F2yz_Py;
  Double I_ERI_D2y_Py_F2yz_Py = I_ERI_F3y_S_F2yz_Py+ABY*I_ERI_D2y_S_F2yz_Py;
  Double I_ERI_Dyz_Py_F2yz_Py = I_ERI_F2yz_S_F2yz_Py+ABY*I_ERI_Dyz_S_F2yz_Py;
  Double I_ERI_D2z_Py_F2yz_Py = I_ERI_Fy2z_S_F2yz_Py+ABY*I_ERI_D2z_S_F2yz_Py;
  Double I_ERI_D2x_Pz_F2yz_Py = I_ERI_F2xz_S_F2yz_Py+ABZ*I_ERI_D2x_S_F2yz_Py;
  Double I_ERI_Dxy_Pz_F2yz_Py = I_ERI_Fxyz_S_F2yz_Py+ABZ*I_ERI_Dxy_S_F2yz_Py;
  Double I_ERI_Dxz_Pz_F2yz_Py = I_ERI_Fx2z_S_F2yz_Py+ABZ*I_ERI_Dxz_S_F2yz_Py;
  Double I_ERI_D2y_Pz_F2yz_Py = I_ERI_F2yz_S_F2yz_Py+ABZ*I_ERI_D2y_S_F2yz_Py;
  Double I_ERI_Dyz_Pz_F2yz_Py = I_ERI_Fy2z_S_F2yz_Py+ABZ*I_ERI_Dyz_S_F2yz_Py;
  Double I_ERI_D2z_Pz_F2yz_Py = I_ERI_F3z_S_F2yz_Py+ABZ*I_ERI_D2z_S_F2yz_Py;
  Double I_ERI_D2x_Px_Fy2z_Py = I_ERI_F3x_S_Fy2z_Py+ABX*I_ERI_D2x_S_Fy2z_Py;
  Double I_ERI_Dxy_Px_Fy2z_Py = I_ERI_F2xy_S_Fy2z_Py+ABX*I_ERI_Dxy_S_Fy2z_Py;
  Double I_ERI_Dxz_Px_Fy2z_Py = I_ERI_F2xz_S_Fy2z_Py+ABX*I_ERI_Dxz_S_Fy2z_Py;
  Double I_ERI_D2y_Px_Fy2z_Py = I_ERI_Fx2y_S_Fy2z_Py+ABX*I_ERI_D2y_S_Fy2z_Py;
  Double I_ERI_Dyz_Px_Fy2z_Py = I_ERI_Fxyz_S_Fy2z_Py+ABX*I_ERI_Dyz_S_Fy2z_Py;
  Double I_ERI_D2z_Px_Fy2z_Py = I_ERI_Fx2z_S_Fy2z_Py+ABX*I_ERI_D2z_S_Fy2z_Py;
  Double I_ERI_D2x_Py_Fy2z_Py = I_ERI_F2xy_S_Fy2z_Py+ABY*I_ERI_D2x_S_Fy2z_Py;
  Double I_ERI_Dxy_Py_Fy2z_Py = I_ERI_Fx2y_S_Fy2z_Py+ABY*I_ERI_Dxy_S_Fy2z_Py;
  Double I_ERI_Dxz_Py_Fy2z_Py = I_ERI_Fxyz_S_Fy2z_Py+ABY*I_ERI_Dxz_S_Fy2z_Py;
  Double I_ERI_D2y_Py_Fy2z_Py = I_ERI_F3y_S_Fy2z_Py+ABY*I_ERI_D2y_S_Fy2z_Py;
  Double I_ERI_Dyz_Py_Fy2z_Py = I_ERI_F2yz_S_Fy2z_Py+ABY*I_ERI_Dyz_S_Fy2z_Py;
  Double I_ERI_D2z_Py_Fy2z_Py = I_ERI_Fy2z_S_Fy2z_Py+ABY*I_ERI_D2z_S_Fy2z_Py;
  Double I_ERI_D2x_Pz_Fy2z_Py = I_ERI_F2xz_S_Fy2z_Py+ABZ*I_ERI_D2x_S_Fy2z_Py;
  Double I_ERI_Dxy_Pz_Fy2z_Py = I_ERI_Fxyz_S_Fy2z_Py+ABZ*I_ERI_Dxy_S_Fy2z_Py;
  Double I_ERI_Dxz_Pz_Fy2z_Py = I_ERI_Fx2z_S_Fy2z_Py+ABZ*I_ERI_Dxz_S_Fy2z_Py;
  Double I_ERI_D2y_Pz_Fy2z_Py = I_ERI_F2yz_S_Fy2z_Py+ABZ*I_ERI_D2y_S_Fy2z_Py;
  Double I_ERI_Dyz_Pz_Fy2z_Py = I_ERI_Fy2z_S_Fy2z_Py+ABZ*I_ERI_Dyz_S_Fy2z_Py;
  Double I_ERI_D2z_Pz_Fy2z_Py = I_ERI_F3z_S_Fy2z_Py+ABZ*I_ERI_D2z_S_Fy2z_Py;
  Double I_ERI_D2x_Px_F3z_Py = I_ERI_F3x_S_F3z_Py+ABX*I_ERI_D2x_S_F3z_Py;
  Double I_ERI_Dxy_Px_F3z_Py = I_ERI_F2xy_S_F3z_Py+ABX*I_ERI_Dxy_S_F3z_Py;
  Double I_ERI_Dxz_Px_F3z_Py = I_ERI_F2xz_S_F3z_Py+ABX*I_ERI_Dxz_S_F3z_Py;
  Double I_ERI_D2y_Px_F3z_Py = I_ERI_Fx2y_S_F3z_Py+ABX*I_ERI_D2y_S_F3z_Py;
  Double I_ERI_Dyz_Px_F3z_Py = I_ERI_Fxyz_S_F3z_Py+ABX*I_ERI_Dyz_S_F3z_Py;
  Double I_ERI_D2z_Px_F3z_Py = I_ERI_Fx2z_S_F3z_Py+ABX*I_ERI_D2z_S_F3z_Py;
  Double I_ERI_D2x_Py_F3z_Py = I_ERI_F2xy_S_F3z_Py+ABY*I_ERI_D2x_S_F3z_Py;
  Double I_ERI_Dxy_Py_F3z_Py = I_ERI_Fx2y_S_F3z_Py+ABY*I_ERI_Dxy_S_F3z_Py;
  Double I_ERI_Dxz_Py_F3z_Py = I_ERI_Fxyz_S_F3z_Py+ABY*I_ERI_Dxz_S_F3z_Py;
  Double I_ERI_D2y_Py_F3z_Py = I_ERI_F3y_S_F3z_Py+ABY*I_ERI_D2y_S_F3z_Py;
  Double I_ERI_Dyz_Py_F3z_Py = I_ERI_F2yz_S_F3z_Py+ABY*I_ERI_Dyz_S_F3z_Py;
  Double I_ERI_D2z_Py_F3z_Py = I_ERI_Fy2z_S_F3z_Py+ABY*I_ERI_D2z_S_F3z_Py;
  Double I_ERI_D2x_Pz_F3z_Py = I_ERI_F2xz_S_F3z_Py+ABZ*I_ERI_D2x_S_F3z_Py;
  Double I_ERI_Dxy_Pz_F3z_Py = I_ERI_Fxyz_S_F3z_Py+ABZ*I_ERI_Dxy_S_F3z_Py;
  Double I_ERI_Dxz_Pz_F3z_Py = I_ERI_Fx2z_S_F3z_Py+ABZ*I_ERI_Dxz_S_F3z_Py;
  Double I_ERI_D2y_Pz_F3z_Py = I_ERI_F2yz_S_F3z_Py+ABZ*I_ERI_D2y_S_F3z_Py;
  Double I_ERI_Dyz_Pz_F3z_Py = I_ERI_Fy2z_S_F3z_Py+ABZ*I_ERI_Dyz_S_F3z_Py;
  Double I_ERI_D2z_Pz_F3z_Py = I_ERI_F3z_S_F3z_Py+ABZ*I_ERI_D2z_S_F3z_Py;
  Double I_ERI_D2x_Px_F3x_Pz = I_ERI_F3x_S_F3x_Pz+ABX*I_ERI_D2x_S_F3x_Pz;
  Double I_ERI_Dxy_Px_F3x_Pz = I_ERI_F2xy_S_F3x_Pz+ABX*I_ERI_Dxy_S_F3x_Pz;
  Double I_ERI_Dxz_Px_F3x_Pz = I_ERI_F2xz_S_F3x_Pz+ABX*I_ERI_Dxz_S_F3x_Pz;
  Double I_ERI_D2y_Px_F3x_Pz = I_ERI_Fx2y_S_F3x_Pz+ABX*I_ERI_D2y_S_F3x_Pz;
  Double I_ERI_Dyz_Px_F3x_Pz = I_ERI_Fxyz_S_F3x_Pz+ABX*I_ERI_Dyz_S_F3x_Pz;
  Double I_ERI_D2z_Px_F3x_Pz = I_ERI_Fx2z_S_F3x_Pz+ABX*I_ERI_D2z_S_F3x_Pz;
  Double I_ERI_D2x_Py_F3x_Pz = I_ERI_F2xy_S_F3x_Pz+ABY*I_ERI_D2x_S_F3x_Pz;
  Double I_ERI_Dxy_Py_F3x_Pz = I_ERI_Fx2y_S_F3x_Pz+ABY*I_ERI_Dxy_S_F3x_Pz;
  Double I_ERI_Dxz_Py_F3x_Pz = I_ERI_Fxyz_S_F3x_Pz+ABY*I_ERI_Dxz_S_F3x_Pz;
  Double I_ERI_D2y_Py_F3x_Pz = I_ERI_F3y_S_F3x_Pz+ABY*I_ERI_D2y_S_F3x_Pz;
  Double I_ERI_Dyz_Py_F3x_Pz = I_ERI_F2yz_S_F3x_Pz+ABY*I_ERI_Dyz_S_F3x_Pz;
  Double I_ERI_D2z_Py_F3x_Pz = I_ERI_Fy2z_S_F3x_Pz+ABY*I_ERI_D2z_S_F3x_Pz;
  Double I_ERI_D2x_Pz_F3x_Pz = I_ERI_F2xz_S_F3x_Pz+ABZ*I_ERI_D2x_S_F3x_Pz;
  Double I_ERI_Dxy_Pz_F3x_Pz = I_ERI_Fxyz_S_F3x_Pz+ABZ*I_ERI_Dxy_S_F3x_Pz;
  Double I_ERI_Dxz_Pz_F3x_Pz = I_ERI_Fx2z_S_F3x_Pz+ABZ*I_ERI_Dxz_S_F3x_Pz;
  Double I_ERI_D2y_Pz_F3x_Pz = I_ERI_F2yz_S_F3x_Pz+ABZ*I_ERI_D2y_S_F3x_Pz;
  Double I_ERI_Dyz_Pz_F3x_Pz = I_ERI_Fy2z_S_F3x_Pz+ABZ*I_ERI_Dyz_S_F3x_Pz;
  Double I_ERI_D2z_Pz_F3x_Pz = I_ERI_F3z_S_F3x_Pz+ABZ*I_ERI_D2z_S_F3x_Pz;
  Double I_ERI_D2x_Px_F2xy_Pz = I_ERI_F3x_S_F2xy_Pz+ABX*I_ERI_D2x_S_F2xy_Pz;
  Double I_ERI_Dxy_Px_F2xy_Pz = I_ERI_F2xy_S_F2xy_Pz+ABX*I_ERI_Dxy_S_F2xy_Pz;
  Double I_ERI_Dxz_Px_F2xy_Pz = I_ERI_F2xz_S_F2xy_Pz+ABX*I_ERI_Dxz_S_F2xy_Pz;
  Double I_ERI_D2y_Px_F2xy_Pz = I_ERI_Fx2y_S_F2xy_Pz+ABX*I_ERI_D2y_S_F2xy_Pz;
  Double I_ERI_Dyz_Px_F2xy_Pz = I_ERI_Fxyz_S_F2xy_Pz+ABX*I_ERI_Dyz_S_F2xy_Pz;
  Double I_ERI_D2z_Px_F2xy_Pz = I_ERI_Fx2z_S_F2xy_Pz+ABX*I_ERI_D2z_S_F2xy_Pz;
  Double I_ERI_D2x_Py_F2xy_Pz = I_ERI_F2xy_S_F2xy_Pz+ABY*I_ERI_D2x_S_F2xy_Pz;
  Double I_ERI_Dxy_Py_F2xy_Pz = I_ERI_Fx2y_S_F2xy_Pz+ABY*I_ERI_Dxy_S_F2xy_Pz;
  Double I_ERI_Dxz_Py_F2xy_Pz = I_ERI_Fxyz_S_F2xy_Pz+ABY*I_ERI_Dxz_S_F2xy_Pz;
  Double I_ERI_D2y_Py_F2xy_Pz = I_ERI_F3y_S_F2xy_Pz+ABY*I_ERI_D2y_S_F2xy_Pz;
  Double I_ERI_Dyz_Py_F2xy_Pz = I_ERI_F2yz_S_F2xy_Pz+ABY*I_ERI_Dyz_S_F2xy_Pz;
  Double I_ERI_D2z_Py_F2xy_Pz = I_ERI_Fy2z_S_F2xy_Pz+ABY*I_ERI_D2z_S_F2xy_Pz;
  Double I_ERI_D2x_Pz_F2xy_Pz = I_ERI_F2xz_S_F2xy_Pz+ABZ*I_ERI_D2x_S_F2xy_Pz;
  Double I_ERI_Dxy_Pz_F2xy_Pz = I_ERI_Fxyz_S_F2xy_Pz+ABZ*I_ERI_Dxy_S_F2xy_Pz;
  Double I_ERI_Dxz_Pz_F2xy_Pz = I_ERI_Fx2z_S_F2xy_Pz+ABZ*I_ERI_Dxz_S_F2xy_Pz;
  Double I_ERI_D2y_Pz_F2xy_Pz = I_ERI_F2yz_S_F2xy_Pz+ABZ*I_ERI_D2y_S_F2xy_Pz;
  Double I_ERI_Dyz_Pz_F2xy_Pz = I_ERI_Fy2z_S_F2xy_Pz+ABZ*I_ERI_Dyz_S_F2xy_Pz;
  Double I_ERI_D2z_Pz_F2xy_Pz = I_ERI_F3z_S_F2xy_Pz+ABZ*I_ERI_D2z_S_F2xy_Pz;
  Double I_ERI_D2x_Px_F2xz_Pz = I_ERI_F3x_S_F2xz_Pz+ABX*I_ERI_D2x_S_F2xz_Pz;
  Double I_ERI_Dxy_Px_F2xz_Pz = I_ERI_F2xy_S_F2xz_Pz+ABX*I_ERI_Dxy_S_F2xz_Pz;
  Double I_ERI_Dxz_Px_F2xz_Pz = I_ERI_F2xz_S_F2xz_Pz+ABX*I_ERI_Dxz_S_F2xz_Pz;
  Double I_ERI_D2y_Px_F2xz_Pz = I_ERI_Fx2y_S_F2xz_Pz+ABX*I_ERI_D2y_S_F2xz_Pz;
  Double I_ERI_Dyz_Px_F2xz_Pz = I_ERI_Fxyz_S_F2xz_Pz+ABX*I_ERI_Dyz_S_F2xz_Pz;
  Double I_ERI_D2z_Px_F2xz_Pz = I_ERI_Fx2z_S_F2xz_Pz+ABX*I_ERI_D2z_S_F2xz_Pz;
  Double I_ERI_D2x_Py_F2xz_Pz = I_ERI_F2xy_S_F2xz_Pz+ABY*I_ERI_D2x_S_F2xz_Pz;
  Double I_ERI_Dxy_Py_F2xz_Pz = I_ERI_Fx2y_S_F2xz_Pz+ABY*I_ERI_Dxy_S_F2xz_Pz;
  Double I_ERI_Dxz_Py_F2xz_Pz = I_ERI_Fxyz_S_F2xz_Pz+ABY*I_ERI_Dxz_S_F2xz_Pz;
  Double I_ERI_D2y_Py_F2xz_Pz = I_ERI_F3y_S_F2xz_Pz+ABY*I_ERI_D2y_S_F2xz_Pz;
  Double I_ERI_Dyz_Py_F2xz_Pz = I_ERI_F2yz_S_F2xz_Pz+ABY*I_ERI_Dyz_S_F2xz_Pz;
  Double I_ERI_D2z_Py_F2xz_Pz = I_ERI_Fy2z_S_F2xz_Pz+ABY*I_ERI_D2z_S_F2xz_Pz;
  Double I_ERI_D2x_Pz_F2xz_Pz = I_ERI_F2xz_S_F2xz_Pz+ABZ*I_ERI_D2x_S_F2xz_Pz;
  Double I_ERI_Dxy_Pz_F2xz_Pz = I_ERI_Fxyz_S_F2xz_Pz+ABZ*I_ERI_Dxy_S_F2xz_Pz;
  Double I_ERI_Dxz_Pz_F2xz_Pz = I_ERI_Fx2z_S_F2xz_Pz+ABZ*I_ERI_Dxz_S_F2xz_Pz;
  Double I_ERI_D2y_Pz_F2xz_Pz = I_ERI_F2yz_S_F2xz_Pz+ABZ*I_ERI_D2y_S_F2xz_Pz;
  Double I_ERI_Dyz_Pz_F2xz_Pz = I_ERI_Fy2z_S_F2xz_Pz+ABZ*I_ERI_Dyz_S_F2xz_Pz;
  Double I_ERI_D2z_Pz_F2xz_Pz = I_ERI_F3z_S_F2xz_Pz+ABZ*I_ERI_D2z_S_F2xz_Pz;
  Double I_ERI_D2x_Px_Fx2y_Pz = I_ERI_F3x_S_Fx2y_Pz+ABX*I_ERI_D2x_S_Fx2y_Pz;
  Double I_ERI_Dxy_Px_Fx2y_Pz = I_ERI_F2xy_S_Fx2y_Pz+ABX*I_ERI_Dxy_S_Fx2y_Pz;
  Double I_ERI_Dxz_Px_Fx2y_Pz = I_ERI_F2xz_S_Fx2y_Pz+ABX*I_ERI_Dxz_S_Fx2y_Pz;
  Double I_ERI_D2y_Px_Fx2y_Pz = I_ERI_Fx2y_S_Fx2y_Pz+ABX*I_ERI_D2y_S_Fx2y_Pz;
  Double I_ERI_Dyz_Px_Fx2y_Pz = I_ERI_Fxyz_S_Fx2y_Pz+ABX*I_ERI_Dyz_S_Fx2y_Pz;
  Double I_ERI_D2z_Px_Fx2y_Pz = I_ERI_Fx2z_S_Fx2y_Pz+ABX*I_ERI_D2z_S_Fx2y_Pz;
  Double I_ERI_D2x_Py_Fx2y_Pz = I_ERI_F2xy_S_Fx2y_Pz+ABY*I_ERI_D2x_S_Fx2y_Pz;
  Double I_ERI_Dxy_Py_Fx2y_Pz = I_ERI_Fx2y_S_Fx2y_Pz+ABY*I_ERI_Dxy_S_Fx2y_Pz;
  Double I_ERI_Dxz_Py_Fx2y_Pz = I_ERI_Fxyz_S_Fx2y_Pz+ABY*I_ERI_Dxz_S_Fx2y_Pz;
  Double I_ERI_D2y_Py_Fx2y_Pz = I_ERI_F3y_S_Fx2y_Pz+ABY*I_ERI_D2y_S_Fx2y_Pz;
  Double I_ERI_Dyz_Py_Fx2y_Pz = I_ERI_F2yz_S_Fx2y_Pz+ABY*I_ERI_Dyz_S_Fx2y_Pz;
  Double I_ERI_D2z_Py_Fx2y_Pz = I_ERI_Fy2z_S_Fx2y_Pz+ABY*I_ERI_D2z_S_Fx2y_Pz;
  Double I_ERI_D2x_Pz_Fx2y_Pz = I_ERI_F2xz_S_Fx2y_Pz+ABZ*I_ERI_D2x_S_Fx2y_Pz;
  Double I_ERI_Dxy_Pz_Fx2y_Pz = I_ERI_Fxyz_S_Fx2y_Pz+ABZ*I_ERI_Dxy_S_Fx2y_Pz;
  Double I_ERI_Dxz_Pz_Fx2y_Pz = I_ERI_Fx2z_S_Fx2y_Pz+ABZ*I_ERI_Dxz_S_Fx2y_Pz;
  Double I_ERI_D2y_Pz_Fx2y_Pz = I_ERI_F2yz_S_Fx2y_Pz+ABZ*I_ERI_D2y_S_Fx2y_Pz;
  Double I_ERI_Dyz_Pz_Fx2y_Pz = I_ERI_Fy2z_S_Fx2y_Pz+ABZ*I_ERI_Dyz_S_Fx2y_Pz;
  Double I_ERI_D2z_Pz_Fx2y_Pz = I_ERI_F3z_S_Fx2y_Pz+ABZ*I_ERI_D2z_S_Fx2y_Pz;
  Double I_ERI_D2x_Px_Fxyz_Pz = I_ERI_F3x_S_Fxyz_Pz+ABX*I_ERI_D2x_S_Fxyz_Pz;
  Double I_ERI_Dxy_Px_Fxyz_Pz = I_ERI_F2xy_S_Fxyz_Pz+ABX*I_ERI_Dxy_S_Fxyz_Pz;
  Double I_ERI_Dxz_Px_Fxyz_Pz = I_ERI_F2xz_S_Fxyz_Pz+ABX*I_ERI_Dxz_S_Fxyz_Pz;
  Double I_ERI_D2y_Px_Fxyz_Pz = I_ERI_Fx2y_S_Fxyz_Pz+ABX*I_ERI_D2y_S_Fxyz_Pz;
  Double I_ERI_Dyz_Px_Fxyz_Pz = I_ERI_Fxyz_S_Fxyz_Pz+ABX*I_ERI_Dyz_S_Fxyz_Pz;
  Double I_ERI_D2z_Px_Fxyz_Pz = I_ERI_Fx2z_S_Fxyz_Pz+ABX*I_ERI_D2z_S_Fxyz_Pz;
  Double I_ERI_D2x_Py_Fxyz_Pz = I_ERI_F2xy_S_Fxyz_Pz+ABY*I_ERI_D2x_S_Fxyz_Pz;
  Double I_ERI_Dxy_Py_Fxyz_Pz = I_ERI_Fx2y_S_Fxyz_Pz+ABY*I_ERI_Dxy_S_Fxyz_Pz;
  Double I_ERI_Dxz_Py_Fxyz_Pz = I_ERI_Fxyz_S_Fxyz_Pz+ABY*I_ERI_Dxz_S_Fxyz_Pz;
  Double I_ERI_D2y_Py_Fxyz_Pz = I_ERI_F3y_S_Fxyz_Pz+ABY*I_ERI_D2y_S_Fxyz_Pz;
  Double I_ERI_Dyz_Py_Fxyz_Pz = I_ERI_F2yz_S_Fxyz_Pz+ABY*I_ERI_Dyz_S_Fxyz_Pz;
  Double I_ERI_D2z_Py_Fxyz_Pz = I_ERI_Fy2z_S_Fxyz_Pz+ABY*I_ERI_D2z_S_Fxyz_Pz;
  Double I_ERI_D2x_Pz_Fxyz_Pz = I_ERI_F2xz_S_Fxyz_Pz+ABZ*I_ERI_D2x_S_Fxyz_Pz;
  Double I_ERI_Dxy_Pz_Fxyz_Pz = I_ERI_Fxyz_S_Fxyz_Pz+ABZ*I_ERI_Dxy_S_Fxyz_Pz;
  Double I_ERI_Dxz_Pz_Fxyz_Pz = I_ERI_Fx2z_S_Fxyz_Pz+ABZ*I_ERI_Dxz_S_Fxyz_Pz;
  Double I_ERI_D2y_Pz_Fxyz_Pz = I_ERI_F2yz_S_Fxyz_Pz+ABZ*I_ERI_D2y_S_Fxyz_Pz;
  Double I_ERI_Dyz_Pz_Fxyz_Pz = I_ERI_Fy2z_S_Fxyz_Pz+ABZ*I_ERI_Dyz_S_Fxyz_Pz;
  Double I_ERI_D2z_Pz_Fxyz_Pz = I_ERI_F3z_S_Fxyz_Pz+ABZ*I_ERI_D2z_S_Fxyz_Pz;
  Double I_ERI_D2x_Px_Fx2z_Pz = I_ERI_F3x_S_Fx2z_Pz+ABX*I_ERI_D2x_S_Fx2z_Pz;
  Double I_ERI_Dxy_Px_Fx2z_Pz = I_ERI_F2xy_S_Fx2z_Pz+ABX*I_ERI_Dxy_S_Fx2z_Pz;
  Double I_ERI_Dxz_Px_Fx2z_Pz = I_ERI_F2xz_S_Fx2z_Pz+ABX*I_ERI_Dxz_S_Fx2z_Pz;
  Double I_ERI_D2y_Px_Fx2z_Pz = I_ERI_Fx2y_S_Fx2z_Pz+ABX*I_ERI_D2y_S_Fx2z_Pz;
  Double I_ERI_Dyz_Px_Fx2z_Pz = I_ERI_Fxyz_S_Fx2z_Pz+ABX*I_ERI_Dyz_S_Fx2z_Pz;
  Double I_ERI_D2z_Px_Fx2z_Pz = I_ERI_Fx2z_S_Fx2z_Pz+ABX*I_ERI_D2z_S_Fx2z_Pz;
  Double I_ERI_D2x_Py_Fx2z_Pz = I_ERI_F2xy_S_Fx2z_Pz+ABY*I_ERI_D2x_S_Fx2z_Pz;
  Double I_ERI_Dxy_Py_Fx2z_Pz = I_ERI_Fx2y_S_Fx2z_Pz+ABY*I_ERI_Dxy_S_Fx2z_Pz;
  Double I_ERI_Dxz_Py_Fx2z_Pz = I_ERI_Fxyz_S_Fx2z_Pz+ABY*I_ERI_Dxz_S_Fx2z_Pz;
  Double I_ERI_D2y_Py_Fx2z_Pz = I_ERI_F3y_S_Fx2z_Pz+ABY*I_ERI_D2y_S_Fx2z_Pz;
  Double I_ERI_Dyz_Py_Fx2z_Pz = I_ERI_F2yz_S_Fx2z_Pz+ABY*I_ERI_Dyz_S_Fx2z_Pz;
  Double I_ERI_D2z_Py_Fx2z_Pz = I_ERI_Fy2z_S_Fx2z_Pz+ABY*I_ERI_D2z_S_Fx2z_Pz;
  Double I_ERI_D2x_Pz_Fx2z_Pz = I_ERI_F2xz_S_Fx2z_Pz+ABZ*I_ERI_D2x_S_Fx2z_Pz;
  Double I_ERI_Dxy_Pz_Fx2z_Pz = I_ERI_Fxyz_S_Fx2z_Pz+ABZ*I_ERI_Dxy_S_Fx2z_Pz;
  Double I_ERI_Dxz_Pz_Fx2z_Pz = I_ERI_Fx2z_S_Fx2z_Pz+ABZ*I_ERI_Dxz_S_Fx2z_Pz;
  Double I_ERI_D2y_Pz_Fx2z_Pz = I_ERI_F2yz_S_Fx2z_Pz+ABZ*I_ERI_D2y_S_Fx2z_Pz;
  Double I_ERI_Dyz_Pz_Fx2z_Pz = I_ERI_Fy2z_S_Fx2z_Pz+ABZ*I_ERI_Dyz_S_Fx2z_Pz;
  Double I_ERI_D2z_Pz_Fx2z_Pz = I_ERI_F3z_S_Fx2z_Pz+ABZ*I_ERI_D2z_S_Fx2z_Pz;
  Double I_ERI_D2x_Px_F3y_Pz = I_ERI_F3x_S_F3y_Pz+ABX*I_ERI_D2x_S_F3y_Pz;
  Double I_ERI_Dxy_Px_F3y_Pz = I_ERI_F2xy_S_F3y_Pz+ABX*I_ERI_Dxy_S_F3y_Pz;
  Double I_ERI_Dxz_Px_F3y_Pz = I_ERI_F2xz_S_F3y_Pz+ABX*I_ERI_Dxz_S_F3y_Pz;
  Double I_ERI_D2y_Px_F3y_Pz = I_ERI_Fx2y_S_F3y_Pz+ABX*I_ERI_D2y_S_F3y_Pz;
  Double I_ERI_Dyz_Px_F3y_Pz = I_ERI_Fxyz_S_F3y_Pz+ABX*I_ERI_Dyz_S_F3y_Pz;
  Double I_ERI_D2z_Px_F3y_Pz = I_ERI_Fx2z_S_F3y_Pz+ABX*I_ERI_D2z_S_F3y_Pz;
  Double I_ERI_D2x_Py_F3y_Pz = I_ERI_F2xy_S_F3y_Pz+ABY*I_ERI_D2x_S_F3y_Pz;
  Double I_ERI_Dxy_Py_F3y_Pz = I_ERI_Fx2y_S_F3y_Pz+ABY*I_ERI_Dxy_S_F3y_Pz;
  Double I_ERI_Dxz_Py_F3y_Pz = I_ERI_Fxyz_S_F3y_Pz+ABY*I_ERI_Dxz_S_F3y_Pz;
  Double I_ERI_D2y_Py_F3y_Pz = I_ERI_F3y_S_F3y_Pz+ABY*I_ERI_D2y_S_F3y_Pz;
  Double I_ERI_Dyz_Py_F3y_Pz = I_ERI_F2yz_S_F3y_Pz+ABY*I_ERI_Dyz_S_F3y_Pz;
  Double I_ERI_D2z_Py_F3y_Pz = I_ERI_Fy2z_S_F3y_Pz+ABY*I_ERI_D2z_S_F3y_Pz;
  Double I_ERI_D2x_Pz_F3y_Pz = I_ERI_F2xz_S_F3y_Pz+ABZ*I_ERI_D2x_S_F3y_Pz;
  Double I_ERI_Dxy_Pz_F3y_Pz = I_ERI_Fxyz_S_F3y_Pz+ABZ*I_ERI_Dxy_S_F3y_Pz;
  Double I_ERI_Dxz_Pz_F3y_Pz = I_ERI_Fx2z_S_F3y_Pz+ABZ*I_ERI_Dxz_S_F3y_Pz;
  Double I_ERI_D2y_Pz_F3y_Pz = I_ERI_F2yz_S_F3y_Pz+ABZ*I_ERI_D2y_S_F3y_Pz;
  Double I_ERI_Dyz_Pz_F3y_Pz = I_ERI_Fy2z_S_F3y_Pz+ABZ*I_ERI_Dyz_S_F3y_Pz;
  Double I_ERI_D2z_Pz_F3y_Pz = I_ERI_F3z_S_F3y_Pz+ABZ*I_ERI_D2z_S_F3y_Pz;
  Double I_ERI_D2x_Px_F2yz_Pz = I_ERI_F3x_S_F2yz_Pz+ABX*I_ERI_D2x_S_F2yz_Pz;
  Double I_ERI_Dxy_Px_F2yz_Pz = I_ERI_F2xy_S_F2yz_Pz+ABX*I_ERI_Dxy_S_F2yz_Pz;
  Double I_ERI_Dxz_Px_F2yz_Pz = I_ERI_F2xz_S_F2yz_Pz+ABX*I_ERI_Dxz_S_F2yz_Pz;
  Double I_ERI_D2y_Px_F2yz_Pz = I_ERI_Fx2y_S_F2yz_Pz+ABX*I_ERI_D2y_S_F2yz_Pz;
  Double I_ERI_Dyz_Px_F2yz_Pz = I_ERI_Fxyz_S_F2yz_Pz+ABX*I_ERI_Dyz_S_F2yz_Pz;
  Double I_ERI_D2z_Px_F2yz_Pz = I_ERI_Fx2z_S_F2yz_Pz+ABX*I_ERI_D2z_S_F2yz_Pz;
  Double I_ERI_D2x_Py_F2yz_Pz = I_ERI_F2xy_S_F2yz_Pz+ABY*I_ERI_D2x_S_F2yz_Pz;
  Double I_ERI_Dxy_Py_F2yz_Pz = I_ERI_Fx2y_S_F2yz_Pz+ABY*I_ERI_Dxy_S_F2yz_Pz;
  Double I_ERI_Dxz_Py_F2yz_Pz = I_ERI_Fxyz_S_F2yz_Pz+ABY*I_ERI_Dxz_S_F2yz_Pz;
  Double I_ERI_D2y_Py_F2yz_Pz = I_ERI_F3y_S_F2yz_Pz+ABY*I_ERI_D2y_S_F2yz_Pz;
  Double I_ERI_Dyz_Py_F2yz_Pz = I_ERI_F2yz_S_F2yz_Pz+ABY*I_ERI_Dyz_S_F2yz_Pz;
  Double I_ERI_D2z_Py_F2yz_Pz = I_ERI_Fy2z_S_F2yz_Pz+ABY*I_ERI_D2z_S_F2yz_Pz;
  Double I_ERI_D2x_Pz_F2yz_Pz = I_ERI_F2xz_S_F2yz_Pz+ABZ*I_ERI_D2x_S_F2yz_Pz;
  Double I_ERI_Dxy_Pz_F2yz_Pz = I_ERI_Fxyz_S_F2yz_Pz+ABZ*I_ERI_Dxy_S_F2yz_Pz;
  Double I_ERI_Dxz_Pz_F2yz_Pz = I_ERI_Fx2z_S_F2yz_Pz+ABZ*I_ERI_Dxz_S_F2yz_Pz;
  Double I_ERI_D2y_Pz_F2yz_Pz = I_ERI_F2yz_S_F2yz_Pz+ABZ*I_ERI_D2y_S_F2yz_Pz;
  Double I_ERI_Dyz_Pz_F2yz_Pz = I_ERI_Fy2z_S_F2yz_Pz+ABZ*I_ERI_Dyz_S_F2yz_Pz;
  Double I_ERI_D2z_Pz_F2yz_Pz = I_ERI_F3z_S_F2yz_Pz+ABZ*I_ERI_D2z_S_F2yz_Pz;
  Double I_ERI_D2x_Px_Fy2z_Pz = I_ERI_F3x_S_Fy2z_Pz+ABX*I_ERI_D2x_S_Fy2z_Pz;
  Double I_ERI_Dxy_Px_Fy2z_Pz = I_ERI_F2xy_S_Fy2z_Pz+ABX*I_ERI_Dxy_S_Fy2z_Pz;
  Double I_ERI_Dxz_Px_Fy2z_Pz = I_ERI_F2xz_S_Fy2z_Pz+ABX*I_ERI_Dxz_S_Fy2z_Pz;
  Double I_ERI_D2y_Px_Fy2z_Pz = I_ERI_Fx2y_S_Fy2z_Pz+ABX*I_ERI_D2y_S_Fy2z_Pz;
  Double I_ERI_Dyz_Px_Fy2z_Pz = I_ERI_Fxyz_S_Fy2z_Pz+ABX*I_ERI_Dyz_S_Fy2z_Pz;
  Double I_ERI_D2z_Px_Fy2z_Pz = I_ERI_Fx2z_S_Fy2z_Pz+ABX*I_ERI_D2z_S_Fy2z_Pz;
  Double I_ERI_D2x_Py_Fy2z_Pz = I_ERI_F2xy_S_Fy2z_Pz+ABY*I_ERI_D2x_S_Fy2z_Pz;
  Double I_ERI_Dxy_Py_Fy2z_Pz = I_ERI_Fx2y_S_Fy2z_Pz+ABY*I_ERI_Dxy_S_Fy2z_Pz;
  Double I_ERI_Dxz_Py_Fy2z_Pz = I_ERI_Fxyz_S_Fy2z_Pz+ABY*I_ERI_Dxz_S_Fy2z_Pz;
  Double I_ERI_D2y_Py_Fy2z_Pz = I_ERI_F3y_S_Fy2z_Pz+ABY*I_ERI_D2y_S_Fy2z_Pz;
  Double I_ERI_Dyz_Py_Fy2z_Pz = I_ERI_F2yz_S_Fy2z_Pz+ABY*I_ERI_Dyz_S_Fy2z_Pz;
  Double I_ERI_D2z_Py_Fy2z_Pz = I_ERI_Fy2z_S_Fy2z_Pz+ABY*I_ERI_D2z_S_Fy2z_Pz;
  Double I_ERI_D2x_Pz_Fy2z_Pz = I_ERI_F2xz_S_Fy2z_Pz+ABZ*I_ERI_D2x_S_Fy2z_Pz;
  Double I_ERI_Dxy_Pz_Fy2z_Pz = I_ERI_Fxyz_S_Fy2z_Pz+ABZ*I_ERI_Dxy_S_Fy2z_Pz;
  Double I_ERI_Dxz_Pz_Fy2z_Pz = I_ERI_Fx2z_S_Fy2z_Pz+ABZ*I_ERI_Dxz_S_Fy2z_Pz;
  Double I_ERI_D2y_Pz_Fy2z_Pz = I_ERI_F2yz_S_Fy2z_Pz+ABZ*I_ERI_D2y_S_Fy2z_Pz;
  Double I_ERI_Dyz_Pz_Fy2z_Pz = I_ERI_Fy2z_S_Fy2z_Pz+ABZ*I_ERI_Dyz_S_Fy2z_Pz;
  Double I_ERI_D2z_Pz_Fy2z_Pz = I_ERI_F3z_S_Fy2z_Pz+ABZ*I_ERI_D2z_S_Fy2z_Pz;
  Double I_ERI_D2x_Px_F3z_Pz = I_ERI_F3x_S_F3z_Pz+ABX*I_ERI_D2x_S_F3z_Pz;
  Double I_ERI_Dxy_Px_F3z_Pz = I_ERI_F2xy_S_F3z_Pz+ABX*I_ERI_Dxy_S_F3z_Pz;
  Double I_ERI_Dxz_Px_F3z_Pz = I_ERI_F2xz_S_F3z_Pz+ABX*I_ERI_Dxz_S_F3z_Pz;
  Double I_ERI_D2y_Px_F3z_Pz = I_ERI_Fx2y_S_F3z_Pz+ABX*I_ERI_D2y_S_F3z_Pz;
  Double I_ERI_Dyz_Px_F3z_Pz = I_ERI_Fxyz_S_F3z_Pz+ABX*I_ERI_Dyz_S_F3z_Pz;
  Double I_ERI_D2z_Px_F3z_Pz = I_ERI_Fx2z_S_F3z_Pz+ABX*I_ERI_D2z_S_F3z_Pz;
  Double I_ERI_D2x_Py_F3z_Pz = I_ERI_F2xy_S_F3z_Pz+ABY*I_ERI_D2x_S_F3z_Pz;
  Double I_ERI_Dxy_Py_F3z_Pz = I_ERI_Fx2y_S_F3z_Pz+ABY*I_ERI_Dxy_S_F3z_Pz;
  Double I_ERI_Dxz_Py_F3z_Pz = I_ERI_Fxyz_S_F3z_Pz+ABY*I_ERI_Dxz_S_F3z_Pz;
  Double I_ERI_D2y_Py_F3z_Pz = I_ERI_F3y_S_F3z_Pz+ABY*I_ERI_D2y_S_F3z_Pz;
  Double I_ERI_Dyz_Py_F3z_Pz = I_ERI_F2yz_S_F3z_Pz+ABY*I_ERI_Dyz_S_F3z_Pz;
  Double I_ERI_D2z_Py_F3z_Pz = I_ERI_Fy2z_S_F3z_Pz+ABY*I_ERI_D2z_S_F3z_Pz;
  Double I_ERI_D2x_Pz_F3z_Pz = I_ERI_F2xz_S_F3z_Pz+ABZ*I_ERI_D2x_S_F3z_Pz;
  Double I_ERI_Dxy_Pz_F3z_Pz = I_ERI_Fxyz_S_F3z_Pz+ABZ*I_ERI_Dxy_S_F3z_Pz;
  Double I_ERI_Dxz_Pz_F3z_Pz = I_ERI_Fx2z_S_F3z_Pz+ABZ*I_ERI_Dxz_S_F3z_Pz;
  Double I_ERI_D2y_Pz_F3z_Pz = I_ERI_F2yz_S_F3z_Pz+ABZ*I_ERI_D2y_S_F3z_Pz;
  Double I_ERI_Dyz_Pz_F3z_Pz = I_ERI_Fy2z_S_F3z_Pz+ABZ*I_ERI_Dyz_S_F3z_Pz;
  Double I_ERI_D2z_Pz_F3z_Pz = I_ERI_F3z_S_F3z_Pz+ABZ*I_ERI_D2z_S_F3z_Pz;

  /************************************************************
   * shell quartet name: SQ_ERI_F_P_F_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 150 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_F_P
   * RHS shell quartet name: SQ_ERI_F_S_F_P
   ************************************************************/
  Double I_ERI_F3x_Px_F3x_Px = I_ERI_G4x_S_F3x_Px+ABX*I_ERI_F3x_S_F3x_Px;
  Double I_ERI_F2xy_Px_F3x_Px = I_ERI_G3xy_S_F3x_Px+ABX*I_ERI_F2xy_S_F3x_Px;
  Double I_ERI_F2xz_Px_F3x_Px = I_ERI_G3xz_S_F3x_Px+ABX*I_ERI_F2xz_S_F3x_Px;
  Double I_ERI_Fx2y_Px_F3x_Px = I_ERI_G2x2y_S_F3x_Px+ABX*I_ERI_Fx2y_S_F3x_Px;
  Double I_ERI_Fxyz_Px_F3x_Px = I_ERI_G2xyz_S_F3x_Px+ABX*I_ERI_Fxyz_S_F3x_Px;
  Double I_ERI_Fx2z_Px_F3x_Px = I_ERI_G2x2z_S_F3x_Px+ABX*I_ERI_Fx2z_S_F3x_Px;
  Double I_ERI_F3y_Px_F3x_Px = I_ERI_Gx3y_S_F3x_Px+ABX*I_ERI_F3y_S_F3x_Px;
  Double I_ERI_F2yz_Px_F3x_Px = I_ERI_Gx2yz_S_F3x_Px+ABX*I_ERI_F2yz_S_F3x_Px;
  Double I_ERI_Fy2z_Px_F3x_Px = I_ERI_Gxy2z_S_F3x_Px+ABX*I_ERI_Fy2z_S_F3x_Px;
  Double I_ERI_F3z_Px_F3x_Px = I_ERI_Gx3z_S_F3x_Px+ABX*I_ERI_F3z_S_F3x_Px;
  Double I_ERI_F2xy_Py_F3x_Px = I_ERI_G2x2y_S_F3x_Px+ABY*I_ERI_F2xy_S_F3x_Px;
  Double I_ERI_F2xz_Py_F3x_Px = I_ERI_G2xyz_S_F3x_Px+ABY*I_ERI_F2xz_S_F3x_Px;
  Double I_ERI_Fx2y_Py_F3x_Px = I_ERI_Gx3y_S_F3x_Px+ABY*I_ERI_Fx2y_S_F3x_Px;
  Double I_ERI_Fxyz_Py_F3x_Px = I_ERI_Gx2yz_S_F3x_Px+ABY*I_ERI_Fxyz_S_F3x_Px;
  Double I_ERI_Fx2z_Py_F3x_Px = I_ERI_Gxy2z_S_F3x_Px+ABY*I_ERI_Fx2z_S_F3x_Px;
  Double I_ERI_F3y_Py_F3x_Px = I_ERI_G4y_S_F3x_Px+ABY*I_ERI_F3y_S_F3x_Px;
  Double I_ERI_F2yz_Py_F3x_Px = I_ERI_G3yz_S_F3x_Px+ABY*I_ERI_F2yz_S_F3x_Px;
  Double I_ERI_Fy2z_Py_F3x_Px = I_ERI_G2y2z_S_F3x_Px+ABY*I_ERI_Fy2z_S_F3x_Px;
  Double I_ERI_F3z_Py_F3x_Px = I_ERI_Gy3z_S_F3x_Px+ABY*I_ERI_F3z_S_F3x_Px;
  Double I_ERI_F2xz_Pz_F3x_Px = I_ERI_G2x2z_S_F3x_Px+ABZ*I_ERI_F2xz_S_F3x_Px;
  Double I_ERI_Fxyz_Pz_F3x_Px = I_ERI_Gxy2z_S_F3x_Px+ABZ*I_ERI_Fxyz_S_F3x_Px;
  Double I_ERI_Fx2z_Pz_F3x_Px = I_ERI_Gx3z_S_F3x_Px+ABZ*I_ERI_Fx2z_S_F3x_Px;
  Double I_ERI_F2yz_Pz_F3x_Px = I_ERI_G2y2z_S_F3x_Px+ABZ*I_ERI_F2yz_S_F3x_Px;
  Double I_ERI_Fy2z_Pz_F3x_Px = I_ERI_Gy3z_S_F3x_Px+ABZ*I_ERI_Fy2z_S_F3x_Px;
  Double I_ERI_F3z_Pz_F3x_Px = I_ERI_G4z_S_F3x_Px+ABZ*I_ERI_F3z_S_F3x_Px;
  Double I_ERI_F3x_Px_F2xy_Px = I_ERI_G4x_S_F2xy_Px+ABX*I_ERI_F3x_S_F2xy_Px;
  Double I_ERI_F2xy_Px_F2xy_Px = I_ERI_G3xy_S_F2xy_Px+ABX*I_ERI_F2xy_S_F2xy_Px;
  Double I_ERI_F2xz_Px_F2xy_Px = I_ERI_G3xz_S_F2xy_Px+ABX*I_ERI_F2xz_S_F2xy_Px;
  Double I_ERI_Fx2y_Px_F2xy_Px = I_ERI_G2x2y_S_F2xy_Px+ABX*I_ERI_Fx2y_S_F2xy_Px;
  Double I_ERI_Fxyz_Px_F2xy_Px = I_ERI_G2xyz_S_F2xy_Px+ABX*I_ERI_Fxyz_S_F2xy_Px;
  Double I_ERI_Fx2z_Px_F2xy_Px = I_ERI_G2x2z_S_F2xy_Px+ABX*I_ERI_Fx2z_S_F2xy_Px;
  Double I_ERI_F3y_Px_F2xy_Px = I_ERI_Gx3y_S_F2xy_Px+ABX*I_ERI_F3y_S_F2xy_Px;
  Double I_ERI_F2yz_Px_F2xy_Px = I_ERI_Gx2yz_S_F2xy_Px+ABX*I_ERI_F2yz_S_F2xy_Px;
  Double I_ERI_Fy2z_Px_F2xy_Px = I_ERI_Gxy2z_S_F2xy_Px+ABX*I_ERI_Fy2z_S_F2xy_Px;
  Double I_ERI_F3z_Px_F2xy_Px = I_ERI_Gx3z_S_F2xy_Px+ABX*I_ERI_F3z_S_F2xy_Px;
  Double I_ERI_F2xy_Py_F2xy_Px = I_ERI_G2x2y_S_F2xy_Px+ABY*I_ERI_F2xy_S_F2xy_Px;
  Double I_ERI_F2xz_Py_F2xy_Px = I_ERI_G2xyz_S_F2xy_Px+ABY*I_ERI_F2xz_S_F2xy_Px;
  Double I_ERI_Fx2y_Py_F2xy_Px = I_ERI_Gx3y_S_F2xy_Px+ABY*I_ERI_Fx2y_S_F2xy_Px;
  Double I_ERI_Fxyz_Py_F2xy_Px = I_ERI_Gx2yz_S_F2xy_Px+ABY*I_ERI_Fxyz_S_F2xy_Px;
  Double I_ERI_Fx2z_Py_F2xy_Px = I_ERI_Gxy2z_S_F2xy_Px+ABY*I_ERI_Fx2z_S_F2xy_Px;
  Double I_ERI_F3y_Py_F2xy_Px = I_ERI_G4y_S_F2xy_Px+ABY*I_ERI_F3y_S_F2xy_Px;
  Double I_ERI_F2yz_Py_F2xy_Px = I_ERI_G3yz_S_F2xy_Px+ABY*I_ERI_F2yz_S_F2xy_Px;
  Double I_ERI_Fy2z_Py_F2xy_Px = I_ERI_G2y2z_S_F2xy_Px+ABY*I_ERI_Fy2z_S_F2xy_Px;
  Double I_ERI_F3z_Py_F2xy_Px = I_ERI_Gy3z_S_F2xy_Px+ABY*I_ERI_F3z_S_F2xy_Px;
  Double I_ERI_F2xz_Pz_F2xy_Px = I_ERI_G2x2z_S_F2xy_Px+ABZ*I_ERI_F2xz_S_F2xy_Px;
  Double I_ERI_Fxyz_Pz_F2xy_Px = I_ERI_Gxy2z_S_F2xy_Px+ABZ*I_ERI_Fxyz_S_F2xy_Px;
  Double I_ERI_Fx2z_Pz_F2xy_Px = I_ERI_Gx3z_S_F2xy_Px+ABZ*I_ERI_Fx2z_S_F2xy_Px;
  Double I_ERI_F2yz_Pz_F2xy_Px = I_ERI_G2y2z_S_F2xy_Px+ABZ*I_ERI_F2yz_S_F2xy_Px;
  Double I_ERI_Fy2z_Pz_F2xy_Px = I_ERI_Gy3z_S_F2xy_Px+ABZ*I_ERI_Fy2z_S_F2xy_Px;
  Double I_ERI_F3z_Pz_F2xy_Px = I_ERI_G4z_S_F2xy_Px+ABZ*I_ERI_F3z_S_F2xy_Px;
  Double I_ERI_F3x_Px_F2xz_Px = I_ERI_G4x_S_F2xz_Px+ABX*I_ERI_F3x_S_F2xz_Px;
  Double I_ERI_F2xy_Px_F2xz_Px = I_ERI_G3xy_S_F2xz_Px+ABX*I_ERI_F2xy_S_F2xz_Px;
  Double I_ERI_F2xz_Px_F2xz_Px = I_ERI_G3xz_S_F2xz_Px+ABX*I_ERI_F2xz_S_F2xz_Px;
  Double I_ERI_Fx2y_Px_F2xz_Px = I_ERI_G2x2y_S_F2xz_Px+ABX*I_ERI_Fx2y_S_F2xz_Px;
  Double I_ERI_Fxyz_Px_F2xz_Px = I_ERI_G2xyz_S_F2xz_Px+ABX*I_ERI_Fxyz_S_F2xz_Px;
  Double I_ERI_Fx2z_Px_F2xz_Px = I_ERI_G2x2z_S_F2xz_Px+ABX*I_ERI_Fx2z_S_F2xz_Px;
  Double I_ERI_F3y_Px_F2xz_Px = I_ERI_Gx3y_S_F2xz_Px+ABX*I_ERI_F3y_S_F2xz_Px;
  Double I_ERI_F2yz_Px_F2xz_Px = I_ERI_Gx2yz_S_F2xz_Px+ABX*I_ERI_F2yz_S_F2xz_Px;
  Double I_ERI_Fy2z_Px_F2xz_Px = I_ERI_Gxy2z_S_F2xz_Px+ABX*I_ERI_Fy2z_S_F2xz_Px;
  Double I_ERI_F3z_Px_F2xz_Px = I_ERI_Gx3z_S_F2xz_Px+ABX*I_ERI_F3z_S_F2xz_Px;
  Double I_ERI_F2xy_Py_F2xz_Px = I_ERI_G2x2y_S_F2xz_Px+ABY*I_ERI_F2xy_S_F2xz_Px;
  Double I_ERI_F2xz_Py_F2xz_Px = I_ERI_G2xyz_S_F2xz_Px+ABY*I_ERI_F2xz_S_F2xz_Px;
  Double I_ERI_Fx2y_Py_F2xz_Px = I_ERI_Gx3y_S_F2xz_Px+ABY*I_ERI_Fx2y_S_F2xz_Px;
  Double I_ERI_Fxyz_Py_F2xz_Px = I_ERI_Gx2yz_S_F2xz_Px+ABY*I_ERI_Fxyz_S_F2xz_Px;
  Double I_ERI_Fx2z_Py_F2xz_Px = I_ERI_Gxy2z_S_F2xz_Px+ABY*I_ERI_Fx2z_S_F2xz_Px;
  Double I_ERI_F3y_Py_F2xz_Px = I_ERI_G4y_S_F2xz_Px+ABY*I_ERI_F3y_S_F2xz_Px;
  Double I_ERI_F2yz_Py_F2xz_Px = I_ERI_G3yz_S_F2xz_Px+ABY*I_ERI_F2yz_S_F2xz_Px;
  Double I_ERI_Fy2z_Py_F2xz_Px = I_ERI_G2y2z_S_F2xz_Px+ABY*I_ERI_Fy2z_S_F2xz_Px;
  Double I_ERI_F3z_Py_F2xz_Px = I_ERI_Gy3z_S_F2xz_Px+ABY*I_ERI_F3z_S_F2xz_Px;
  Double I_ERI_F2xz_Pz_F2xz_Px = I_ERI_G2x2z_S_F2xz_Px+ABZ*I_ERI_F2xz_S_F2xz_Px;
  Double I_ERI_Fxyz_Pz_F2xz_Px = I_ERI_Gxy2z_S_F2xz_Px+ABZ*I_ERI_Fxyz_S_F2xz_Px;
  Double I_ERI_Fx2z_Pz_F2xz_Px = I_ERI_Gx3z_S_F2xz_Px+ABZ*I_ERI_Fx2z_S_F2xz_Px;
  Double I_ERI_F2yz_Pz_F2xz_Px = I_ERI_G2y2z_S_F2xz_Px+ABZ*I_ERI_F2yz_S_F2xz_Px;
  Double I_ERI_Fy2z_Pz_F2xz_Px = I_ERI_Gy3z_S_F2xz_Px+ABZ*I_ERI_Fy2z_S_F2xz_Px;
  Double I_ERI_F3z_Pz_F2xz_Px = I_ERI_G4z_S_F2xz_Px+ABZ*I_ERI_F3z_S_F2xz_Px;
  Double I_ERI_F3x_Px_Fx2y_Px = I_ERI_G4x_S_Fx2y_Px+ABX*I_ERI_F3x_S_Fx2y_Px;
  Double I_ERI_F2xy_Px_Fx2y_Px = I_ERI_G3xy_S_Fx2y_Px+ABX*I_ERI_F2xy_S_Fx2y_Px;
  Double I_ERI_F2xz_Px_Fx2y_Px = I_ERI_G3xz_S_Fx2y_Px+ABX*I_ERI_F2xz_S_Fx2y_Px;
  Double I_ERI_Fx2y_Px_Fx2y_Px = I_ERI_G2x2y_S_Fx2y_Px+ABX*I_ERI_Fx2y_S_Fx2y_Px;
  Double I_ERI_Fxyz_Px_Fx2y_Px = I_ERI_G2xyz_S_Fx2y_Px+ABX*I_ERI_Fxyz_S_Fx2y_Px;
  Double I_ERI_Fx2z_Px_Fx2y_Px = I_ERI_G2x2z_S_Fx2y_Px+ABX*I_ERI_Fx2z_S_Fx2y_Px;
  Double I_ERI_F3y_Px_Fx2y_Px = I_ERI_Gx3y_S_Fx2y_Px+ABX*I_ERI_F3y_S_Fx2y_Px;
  Double I_ERI_F2yz_Px_Fx2y_Px = I_ERI_Gx2yz_S_Fx2y_Px+ABX*I_ERI_F2yz_S_Fx2y_Px;
  Double I_ERI_Fy2z_Px_Fx2y_Px = I_ERI_Gxy2z_S_Fx2y_Px+ABX*I_ERI_Fy2z_S_Fx2y_Px;
  Double I_ERI_F3z_Px_Fx2y_Px = I_ERI_Gx3z_S_Fx2y_Px+ABX*I_ERI_F3z_S_Fx2y_Px;
  Double I_ERI_F2xy_Py_Fx2y_Px = I_ERI_G2x2y_S_Fx2y_Px+ABY*I_ERI_F2xy_S_Fx2y_Px;
  Double I_ERI_F2xz_Py_Fx2y_Px = I_ERI_G2xyz_S_Fx2y_Px+ABY*I_ERI_F2xz_S_Fx2y_Px;
  Double I_ERI_Fx2y_Py_Fx2y_Px = I_ERI_Gx3y_S_Fx2y_Px+ABY*I_ERI_Fx2y_S_Fx2y_Px;
  Double I_ERI_Fxyz_Py_Fx2y_Px = I_ERI_Gx2yz_S_Fx2y_Px+ABY*I_ERI_Fxyz_S_Fx2y_Px;
  Double I_ERI_Fx2z_Py_Fx2y_Px = I_ERI_Gxy2z_S_Fx2y_Px+ABY*I_ERI_Fx2z_S_Fx2y_Px;
  Double I_ERI_F3y_Py_Fx2y_Px = I_ERI_G4y_S_Fx2y_Px+ABY*I_ERI_F3y_S_Fx2y_Px;
  Double I_ERI_F2yz_Py_Fx2y_Px = I_ERI_G3yz_S_Fx2y_Px+ABY*I_ERI_F2yz_S_Fx2y_Px;
  Double I_ERI_Fy2z_Py_Fx2y_Px = I_ERI_G2y2z_S_Fx2y_Px+ABY*I_ERI_Fy2z_S_Fx2y_Px;
  Double I_ERI_F3z_Py_Fx2y_Px = I_ERI_Gy3z_S_Fx2y_Px+ABY*I_ERI_F3z_S_Fx2y_Px;
  Double I_ERI_F2xz_Pz_Fx2y_Px = I_ERI_G2x2z_S_Fx2y_Px+ABZ*I_ERI_F2xz_S_Fx2y_Px;
  Double I_ERI_Fxyz_Pz_Fx2y_Px = I_ERI_Gxy2z_S_Fx2y_Px+ABZ*I_ERI_Fxyz_S_Fx2y_Px;
  Double I_ERI_Fx2z_Pz_Fx2y_Px = I_ERI_Gx3z_S_Fx2y_Px+ABZ*I_ERI_Fx2z_S_Fx2y_Px;
  Double I_ERI_F2yz_Pz_Fx2y_Px = I_ERI_G2y2z_S_Fx2y_Px+ABZ*I_ERI_F2yz_S_Fx2y_Px;
  Double I_ERI_Fy2z_Pz_Fx2y_Px = I_ERI_Gy3z_S_Fx2y_Px+ABZ*I_ERI_Fy2z_S_Fx2y_Px;
  Double I_ERI_F3z_Pz_Fx2y_Px = I_ERI_G4z_S_Fx2y_Px+ABZ*I_ERI_F3z_S_Fx2y_Px;
  Double I_ERI_F3x_Px_Fxyz_Px = I_ERI_G4x_S_Fxyz_Px+ABX*I_ERI_F3x_S_Fxyz_Px;
  Double I_ERI_F2xy_Px_Fxyz_Px = I_ERI_G3xy_S_Fxyz_Px+ABX*I_ERI_F2xy_S_Fxyz_Px;
  Double I_ERI_F2xz_Px_Fxyz_Px = I_ERI_G3xz_S_Fxyz_Px+ABX*I_ERI_F2xz_S_Fxyz_Px;
  Double I_ERI_Fx2y_Px_Fxyz_Px = I_ERI_G2x2y_S_Fxyz_Px+ABX*I_ERI_Fx2y_S_Fxyz_Px;
  Double I_ERI_Fxyz_Px_Fxyz_Px = I_ERI_G2xyz_S_Fxyz_Px+ABX*I_ERI_Fxyz_S_Fxyz_Px;
  Double I_ERI_Fx2z_Px_Fxyz_Px = I_ERI_G2x2z_S_Fxyz_Px+ABX*I_ERI_Fx2z_S_Fxyz_Px;
  Double I_ERI_F3y_Px_Fxyz_Px = I_ERI_Gx3y_S_Fxyz_Px+ABX*I_ERI_F3y_S_Fxyz_Px;
  Double I_ERI_F2yz_Px_Fxyz_Px = I_ERI_Gx2yz_S_Fxyz_Px+ABX*I_ERI_F2yz_S_Fxyz_Px;
  Double I_ERI_Fy2z_Px_Fxyz_Px = I_ERI_Gxy2z_S_Fxyz_Px+ABX*I_ERI_Fy2z_S_Fxyz_Px;
  Double I_ERI_F3z_Px_Fxyz_Px = I_ERI_Gx3z_S_Fxyz_Px+ABX*I_ERI_F3z_S_Fxyz_Px;
  Double I_ERI_F2xy_Py_Fxyz_Px = I_ERI_G2x2y_S_Fxyz_Px+ABY*I_ERI_F2xy_S_Fxyz_Px;
  Double I_ERI_F2xz_Py_Fxyz_Px = I_ERI_G2xyz_S_Fxyz_Px+ABY*I_ERI_F2xz_S_Fxyz_Px;
  Double I_ERI_Fx2y_Py_Fxyz_Px = I_ERI_Gx3y_S_Fxyz_Px+ABY*I_ERI_Fx2y_S_Fxyz_Px;
  Double I_ERI_Fxyz_Py_Fxyz_Px = I_ERI_Gx2yz_S_Fxyz_Px+ABY*I_ERI_Fxyz_S_Fxyz_Px;
  Double I_ERI_Fx2z_Py_Fxyz_Px = I_ERI_Gxy2z_S_Fxyz_Px+ABY*I_ERI_Fx2z_S_Fxyz_Px;
  Double I_ERI_F3y_Py_Fxyz_Px = I_ERI_G4y_S_Fxyz_Px+ABY*I_ERI_F3y_S_Fxyz_Px;
  Double I_ERI_F2yz_Py_Fxyz_Px = I_ERI_G3yz_S_Fxyz_Px+ABY*I_ERI_F2yz_S_Fxyz_Px;
  Double I_ERI_Fy2z_Py_Fxyz_Px = I_ERI_G2y2z_S_Fxyz_Px+ABY*I_ERI_Fy2z_S_Fxyz_Px;
  Double I_ERI_F3z_Py_Fxyz_Px = I_ERI_Gy3z_S_Fxyz_Px+ABY*I_ERI_F3z_S_Fxyz_Px;
  Double I_ERI_F2xz_Pz_Fxyz_Px = I_ERI_G2x2z_S_Fxyz_Px+ABZ*I_ERI_F2xz_S_Fxyz_Px;
  Double I_ERI_Fxyz_Pz_Fxyz_Px = I_ERI_Gxy2z_S_Fxyz_Px+ABZ*I_ERI_Fxyz_S_Fxyz_Px;
  Double I_ERI_Fx2z_Pz_Fxyz_Px = I_ERI_Gx3z_S_Fxyz_Px+ABZ*I_ERI_Fx2z_S_Fxyz_Px;
  Double I_ERI_F2yz_Pz_Fxyz_Px = I_ERI_G2y2z_S_Fxyz_Px+ABZ*I_ERI_F2yz_S_Fxyz_Px;
  Double I_ERI_Fy2z_Pz_Fxyz_Px = I_ERI_Gy3z_S_Fxyz_Px+ABZ*I_ERI_Fy2z_S_Fxyz_Px;
  Double I_ERI_F3z_Pz_Fxyz_Px = I_ERI_G4z_S_Fxyz_Px+ABZ*I_ERI_F3z_S_Fxyz_Px;
  Double I_ERI_F3x_Px_Fx2z_Px = I_ERI_G4x_S_Fx2z_Px+ABX*I_ERI_F3x_S_Fx2z_Px;
  Double I_ERI_F2xy_Px_Fx2z_Px = I_ERI_G3xy_S_Fx2z_Px+ABX*I_ERI_F2xy_S_Fx2z_Px;
  Double I_ERI_F2xz_Px_Fx2z_Px = I_ERI_G3xz_S_Fx2z_Px+ABX*I_ERI_F2xz_S_Fx2z_Px;
  Double I_ERI_Fx2y_Px_Fx2z_Px = I_ERI_G2x2y_S_Fx2z_Px+ABX*I_ERI_Fx2y_S_Fx2z_Px;
  Double I_ERI_Fxyz_Px_Fx2z_Px = I_ERI_G2xyz_S_Fx2z_Px+ABX*I_ERI_Fxyz_S_Fx2z_Px;
  Double I_ERI_Fx2z_Px_Fx2z_Px = I_ERI_G2x2z_S_Fx2z_Px+ABX*I_ERI_Fx2z_S_Fx2z_Px;
  Double I_ERI_F3y_Px_Fx2z_Px = I_ERI_Gx3y_S_Fx2z_Px+ABX*I_ERI_F3y_S_Fx2z_Px;
  Double I_ERI_F2yz_Px_Fx2z_Px = I_ERI_Gx2yz_S_Fx2z_Px+ABX*I_ERI_F2yz_S_Fx2z_Px;
  Double I_ERI_Fy2z_Px_Fx2z_Px = I_ERI_Gxy2z_S_Fx2z_Px+ABX*I_ERI_Fy2z_S_Fx2z_Px;
  Double I_ERI_F3z_Px_Fx2z_Px = I_ERI_Gx3z_S_Fx2z_Px+ABX*I_ERI_F3z_S_Fx2z_Px;
  Double I_ERI_F2xy_Py_Fx2z_Px = I_ERI_G2x2y_S_Fx2z_Px+ABY*I_ERI_F2xy_S_Fx2z_Px;
  Double I_ERI_F2xz_Py_Fx2z_Px = I_ERI_G2xyz_S_Fx2z_Px+ABY*I_ERI_F2xz_S_Fx2z_Px;
  Double I_ERI_Fx2y_Py_Fx2z_Px = I_ERI_Gx3y_S_Fx2z_Px+ABY*I_ERI_Fx2y_S_Fx2z_Px;
  Double I_ERI_Fxyz_Py_Fx2z_Px = I_ERI_Gx2yz_S_Fx2z_Px+ABY*I_ERI_Fxyz_S_Fx2z_Px;
  Double I_ERI_Fx2z_Py_Fx2z_Px = I_ERI_Gxy2z_S_Fx2z_Px+ABY*I_ERI_Fx2z_S_Fx2z_Px;
  Double I_ERI_F3y_Py_Fx2z_Px = I_ERI_G4y_S_Fx2z_Px+ABY*I_ERI_F3y_S_Fx2z_Px;
  Double I_ERI_F2yz_Py_Fx2z_Px = I_ERI_G3yz_S_Fx2z_Px+ABY*I_ERI_F2yz_S_Fx2z_Px;
  Double I_ERI_Fy2z_Py_Fx2z_Px = I_ERI_G2y2z_S_Fx2z_Px+ABY*I_ERI_Fy2z_S_Fx2z_Px;
  Double I_ERI_F3z_Py_Fx2z_Px = I_ERI_Gy3z_S_Fx2z_Px+ABY*I_ERI_F3z_S_Fx2z_Px;
  Double I_ERI_F2xz_Pz_Fx2z_Px = I_ERI_G2x2z_S_Fx2z_Px+ABZ*I_ERI_F2xz_S_Fx2z_Px;
  Double I_ERI_Fxyz_Pz_Fx2z_Px = I_ERI_Gxy2z_S_Fx2z_Px+ABZ*I_ERI_Fxyz_S_Fx2z_Px;
  Double I_ERI_Fx2z_Pz_Fx2z_Px = I_ERI_Gx3z_S_Fx2z_Px+ABZ*I_ERI_Fx2z_S_Fx2z_Px;
  Double I_ERI_F2yz_Pz_Fx2z_Px = I_ERI_G2y2z_S_Fx2z_Px+ABZ*I_ERI_F2yz_S_Fx2z_Px;
  Double I_ERI_Fy2z_Pz_Fx2z_Px = I_ERI_Gy3z_S_Fx2z_Px+ABZ*I_ERI_Fy2z_S_Fx2z_Px;
  Double I_ERI_F3z_Pz_Fx2z_Px = I_ERI_G4z_S_Fx2z_Px+ABZ*I_ERI_F3z_S_Fx2z_Px;
  Double I_ERI_F3x_Px_F3y_Px = I_ERI_G4x_S_F3y_Px+ABX*I_ERI_F3x_S_F3y_Px;
  Double I_ERI_F2xy_Px_F3y_Px = I_ERI_G3xy_S_F3y_Px+ABX*I_ERI_F2xy_S_F3y_Px;
  Double I_ERI_F2xz_Px_F3y_Px = I_ERI_G3xz_S_F3y_Px+ABX*I_ERI_F2xz_S_F3y_Px;
  Double I_ERI_Fx2y_Px_F3y_Px = I_ERI_G2x2y_S_F3y_Px+ABX*I_ERI_Fx2y_S_F3y_Px;
  Double I_ERI_Fxyz_Px_F3y_Px = I_ERI_G2xyz_S_F3y_Px+ABX*I_ERI_Fxyz_S_F3y_Px;
  Double I_ERI_Fx2z_Px_F3y_Px = I_ERI_G2x2z_S_F3y_Px+ABX*I_ERI_Fx2z_S_F3y_Px;
  Double I_ERI_F3y_Px_F3y_Px = I_ERI_Gx3y_S_F3y_Px+ABX*I_ERI_F3y_S_F3y_Px;
  Double I_ERI_F2yz_Px_F3y_Px = I_ERI_Gx2yz_S_F3y_Px+ABX*I_ERI_F2yz_S_F3y_Px;
  Double I_ERI_Fy2z_Px_F3y_Px = I_ERI_Gxy2z_S_F3y_Px+ABX*I_ERI_Fy2z_S_F3y_Px;
  Double I_ERI_F3z_Px_F3y_Px = I_ERI_Gx3z_S_F3y_Px+ABX*I_ERI_F3z_S_F3y_Px;
  Double I_ERI_F2xy_Py_F3y_Px = I_ERI_G2x2y_S_F3y_Px+ABY*I_ERI_F2xy_S_F3y_Px;
  Double I_ERI_F2xz_Py_F3y_Px = I_ERI_G2xyz_S_F3y_Px+ABY*I_ERI_F2xz_S_F3y_Px;
  Double I_ERI_Fx2y_Py_F3y_Px = I_ERI_Gx3y_S_F3y_Px+ABY*I_ERI_Fx2y_S_F3y_Px;
  Double I_ERI_Fxyz_Py_F3y_Px = I_ERI_Gx2yz_S_F3y_Px+ABY*I_ERI_Fxyz_S_F3y_Px;
  Double I_ERI_Fx2z_Py_F3y_Px = I_ERI_Gxy2z_S_F3y_Px+ABY*I_ERI_Fx2z_S_F3y_Px;
  Double I_ERI_F3y_Py_F3y_Px = I_ERI_G4y_S_F3y_Px+ABY*I_ERI_F3y_S_F3y_Px;
  Double I_ERI_F2yz_Py_F3y_Px = I_ERI_G3yz_S_F3y_Px+ABY*I_ERI_F2yz_S_F3y_Px;
  Double I_ERI_Fy2z_Py_F3y_Px = I_ERI_G2y2z_S_F3y_Px+ABY*I_ERI_Fy2z_S_F3y_Px;
  Double I_ERI_F3z_Py_F3y_Px = I_ERI_Gy3z_S_F3y_Px+ABY*I_ERI_F3z_S_F3y_Px;
  Double I_ERI_F2xz_Pz_F3y_Px = I_ERI_G2x2z_S_F3y_Px+ABZ*I_ERI_F2xz_S_F3y_Px;
  Double I_ERI_Fxyz_Pz_F3y_Px = I_ERI_Gxy2z_S_F3y_Px+ABZ*I_ERI_Fxyz_S_F3y_Px;
  Double I_ERI_Fx2z_Pz_F3y_Px = I_ERI_Gx3z_S_F3y_Px+ABZ*I_ERI_Fx2z_S_F3y_Px;
  Double I_ERI_F2yz_Pz_F3y_Px = I_ERI_G2y2z_S_F3y_Px+ABZ*I_ERI_F2yz_S_F3y_Px;
  Double I_ERI_Fy2z_Pz_F3y_Px = I_ERI_Gy3z_S_F3y_Px+ABZ*I_ERI_Fy2z_S_F3y_Px;
  Double I_ERI_F3z_Pz_F3y_Px = I_ERI_G4z_S_F3y_Px+ABZ*I_ERI_F3z_S_F3y_Px;
  Double I_ERI_F3x_Px_F2yz_Px = I_ERI_G4x_S_F2yz_Px+ABX*I_ERI_F3x_S_F2yz_Px;
  Double I_ERI_F2xy_Px_F2yz_Px = I_ERI_G3xy_S_F2yz_Px+ABX*I_ERI_F2xy_S_F2yz_Px;
  Double I_ERI_F2xz_Px_F2yz_Px = I_ERI_G3xz_S_F2yz_Px+ABX*I_ERI_F2xz_S_F2yz_Px;
  Double I_ERI_Fx2y_Px_F2yz_Px = I_ERI_G2x2y_S_F2yz_Px+ABX*I_ERI_Fx2y_S_F2yz_Px;
  Double I_ERI_Fxyz_Px_F2yz_Px = I_ERI_G2xyz_S_F2yz_Px+ABX*I_ERI_Fxyz_S_F2yz_Px;
  Double I_ERI_Fx2z_Px_F2yz_Px = I_ERI_G2x2z_S_F2yz_Px+ABX*I_ERI_Fx2z_S_F2yz_Px;
  Double I_ERI_F3y_Px_F2yz_Px = I_ERI_Gx3y_S_F2yz_Px+ABX*I_ERI_F3y_S_F2yz_Px;
  Double I_ERI_F2yz_Px_F2yz_Px = I_ERI_Gx2yz_S_F2yz_Px+ABX*I_ERI_F2yz_S_F2yz_Px;
  Double I_ERI_Fy2z_Px_F2yz_Px = I_ERI_Gxy2z_S_F2yz_Px+ABX*I_ERI_Fy2z_S_F2yz_Px;
  Double I_ERI_F3z_Px_F2yz_Px = I_ERI_Gx3z_S_F2yz_Px+ABX*I_ERI_F3z_S_F2yz_Px;
  Double I_ERI_F2xy_Py_F2yz_Px = I_ERI_G2x2y_S_F2yz_Px+ABY*I_ERI_F2xy_S_F2yz_Px;
  Double I_ERI_F2xz_Py_F2yz_Px = I_ERI_G2xyz_S_F2yz_Px+ABY*I_ERI_F2xz_S_F2yz_Px;
  Double I_ERI_Fx2y_Py_F2yz_Px = I_ERI_Gx3y_S_F2yz_Px+ABY*I_ERI_Fx2y_S_F2yz_Px;
  Double I_ERI_Fxyz_Py_F2yz_Px = I_ERI_Gx2yz_S_F2yz_Px+ABY*I_ERI_Fxyz_S_F2yz_Px;
  Double I_ERI_Fx2z_Py_F2yz_Px = I_ERI_Gxy2z_S_F2yz_Px+ABY*I_ERI_Fx2z_S_F2yz_Px;
  Double I_ERI_F3y_Py_F2yz_Px = I_ERI_G4y_S_F2yz_Px+ABY*I_ERI_F3y_S_F2yz_Px;
  Double I_ERI_F2yz_Py_F2yz_Px = I_ERI_G3yz_S_F2yz_Px+ABY*I_ERI_F2yz_S_F2yz_Px;
  Double I_ERI_Fy2z_Py_F2yz_Px = I_ERI_G2y2z_S_F2yz_Px+ABY*I_ERI_Fy2z_S_F2yz_Px;
  Double I_ERI_F3z_Py_F2yz_Px = I_ERI_Gy3z_S_F2yz_Px+ABY*I_ERI_F3z_S_F2yz_Px;
  Double I_ERI_F2xz_Pz_F2yz_Px = I_ERI_G2x2z_S_F2yz_Px+ABZ*I_ERI_F2xz_S_F2yz_Px;
  Double I_ERI_Fxyz_Pz_F2yz_Px = I_ERI_Gxy2z_S_F2yz_Px+ABZ*I_ERI_Fxyz_S_F2yz_Px;
  Double I_ERI_Fx2z_Pz_F2yz_Px = I_ERI_Gx3z_S_F2yz_Px+ABZ*I_ERI_Fx2z_S_F2yz_Px;
  Double I_ERI_F2yz_Pz_F2yz_Px = I_ERI_G2y2z_S_F2yz_Px+ABZ*I_ERI_F2yz_S_F2yz_Px;
  Double I_ERI_Fy2z_Pz_F2yz_Px = I_ERI_Gy3z_S_F2yz_Px+ABZ*I_ERI_Fy2z_S_F2yz_Px;
  Double I_ERI_F3z_Pz_F2yz_Px = I_ERI_G4z_S_F2yz_Px+ABZ*I_ERI_F3z_S_F2yz_Px;
  Double I_ERI_F3x_Px_Fy2z_Px = I_ERI_G4x_S_Fy2z_Px+ABX*I_ERI_F3x_S_Fy2z_Px;
  Double I_ERI_F2xy_Px_Fy2z_Px = I_ERI_G3xy_S_Fy2z_Px+ABX*I_ERI_F2xy_S_Fy2z_Px;
  Double I_ERI_F2xz_Px_Fy2z_Px = I_ERI_G3xz_S_Fy2z_Px+ABX*I_ERI_F2xz_S_Fy2z_Px;
  Double I_ERI_Fx2y_Px_Fy2z_Px = I_ERI_G2x2y_S_Fy2z_Px+ABX*I_ERI_Fx2y_S_Fy2z_Px;
  Double I_ERI_Fxyz_Px_Fy2z_Px = I_ERI_G2xyz_S_Fy2z_Px+ABX*I_ERI_Fxyz_S_Fy2z_Px;
  Double I_ERI_Fx2z_Px_Fy2z_Px = I_ERI_G2x2z_S_Fy2z_Px+ABX*I_ERI_Fx2z_S_Fy2z_Px;
  Double I_ERI_F3y_Px_Fy2z_Px = I_ERI_Gx3y_S_Fy2z_Px+ABX*I_ERI_F3y_S_Fy2z_Px;
  Double I_ERI_F2yz_Px_Fy2z_Px = I_ERI_Gx2yz_S_Fy2z_Px+ABX*I_ERI_F2yz_S_Fy2z_Px;
  Double I_ERI_Fy2z_Px_Fy2z_Px = I_ERI_Gxy2z_S_Fy2z_Px+ABX*I_ERI_Fy2z_S_Fy2z_Px;
  Double I_ERI_F3z_Px_Fy2z_Px = I_ERI_Gx3z_S_Fy2z_Px+ABX*I_ERI_F3z_S_Fy2z_Px;
  Double I_ERI_F2xy_Py_Fy2z_Px = I_ERI_G2x2y_S_Fy2z_Px+ABY*I_ERI_F2xy_S_Fy2z_Px;
  Double I_ERI_F2xz_Py_Fy2z_Px = I_ERI_G2xyz_S_Fy2z_Px+ABY*I_ERI_F2xz_S_Fy2z_Px;
  Double I_ERI_Fx2y_Py_Fy2z_Px = I_ERI_Gx3y_S_Fy2z_Px+ABY*I_ERI_Fx2y_S_Fy2z_Px;
  Double I_ERI_Fxyz_Py_Fy2z_Px = I_ERI_Gx2yz_S_Fy2z_Px+ABY*I_ERI_Fxyz_S_Fy2z_Px;
  Double I_ERI_Fx2z_Py_Fy2z_Px = I_ERI_Gxy2z_S_Fy2z_Px+ABY*I_ERI_Fx2z_S_Fy2z_Px;
  Double I_ERI_F3y_Py_Fy2z_Px = I_ERI_G4y_S_Fy2z_Px+ABY*I_ERI_F3y_S_Fy2z_Px;
  Double I_ERI_F2yz_Py_Fy2z_Px = I_ERI_G3yz_S_Fy2z_Px+ABY*I_ERI_F2yz_S_Fy2z_Px;
  Double I_ERI_Fy2z_Py_Fy2z_Px = I_ERI_G2y2z_S_Fy2z_Px+ABY*I_ERI_Fy2z_S_Fy2z_Px;
  Double I_ERI_F3z_Py_Fy2z_Px = I_ERI_Gy3z_S_Fy2z_Px+ABY*I_ERI_F3z_S_Fy2z_Px;
  Double I_ERI_F2xz_Pz_Fy2z_Px = I_ERI_G2x2z_S_Fy2z_Px+ABZ*I_ERI_F2xz_S_Fy2z_Px;
  Double I_ERI_Fxyz_Pz_Fy2z_Px = I_ERI_Gxy2z_S_Fy2z_Px+ABZ*I_ERI_Fxyz_S_Fy2z_Px;
  Double I_ERI_Fx2z_Pz_Fy2z_Px = I_ERI_Gx3z_S_Fy2z_Px+ABZ*I_ERI_Fx2z_S_Fy2z_Px;
  Double I_ERI_F2yz_Pz_Fy2z_Px = I_ERI_G2y2z_S_Fy2z_Px+ABZ*I_ERI_F2yz_S_Fy2z_Px;
  Double I_ERI_Fy2z_Pz_Fy2z_Px = I_ERI_Gy3z_S_Fy2z_Px+ABZ*I_ERI_Fy2z_S_Fy2z_Px;
  Double I_ERI_F3z_Pz_Fy2z_Px = I_ERI_G4z_S_Fy2z_Px+ABZ*I_ERI_F3z_S_Fy2z_Px;
  Double I_ERI_F3x_Px_F3z_Px = I_ERI_G4x_S_F3z_Px+ABX*I_ERI_F3x_S_F3z_Px;
  Double I_ERI_F2xy_Px_F3z_Px = I_ERI_G3xy_S_F3z_Px+ABX*I_ERI_F2xy_S_F3z_Px;
  Double I_ERI_F2xz_Px_F3z_Px = I_ERI_G3xz_S_F3z_Px+ABX*I_ERI_F2xz_S_F3z_Px;
  Double I_ERI_Fx2y_Px_F3z_Px = I_ERI_G2x2y_S_F3z_Px+ABX*I_ERI_Fx2y_S_F3z_Px;
  Double I_ERI_Fxyz_Px_F3z_Px = I_ERI_G2xyz_S_F3z_Px+ABX*I_ERI_Fxyz_S_F3z_Px;
  Double I_ERI_Fx2z_Px_F3z_Px = I_ERI_G2x2z_S_F3z_Px+ABX*I_ERI_Fx2z_S_F3z_Px;
  Double I_ERI_F3y_Px_F3z_Px = I_ERI_Gx3y_S_F3z_Px+ABX*I_ERI_F3y_S_F3z_Px;
  Double I_ERI_F2yz_Px_F3z_Px = I_ERI_Gx2yz_S_F3z_Px+ABX*I_ERI_F2yz_S_F3z_Px;
  Double I_ERI_Fy2z_Px_F3z_Px = I_ERI_Gxy2z_S_F3z_Px+ABX*I_ERI_Fy2z_S_F3z_Px;
  Double I_ERI_F3z_Px_F3z_Px = I_ERI_Gx3z_S_F3z_Px+ABX*I_ERI_F3z_S_F3z_Px;
  Double I_ERI_F2xy_Py_F3z_Px = I_ERI_G2x2y_S_F3z_Px+ABY*I_ERI_F2xy_S_F3z_Px;
  Double I_ERI_F2xz_Py_F3z_Px = I_ERI_G2xyz_S_F3z_Px+ABY*I_ERI_F2xz_S_F3z_Px;
  Double I_ERI_Fx2y_Py_F3z_Px = I_ERI_Gx3y_S_F3z_Px+ABY*I_ERI_Fx2y_S_F3z_Px;
  Double I_ERI_Fxyz_Py_F3z_Px = I_ERI_Gx2yz_S_F3z_Px+ABY*I_ERI_Fxyz_S_F3z_Px;
  Double I_ERI_Fx2z_Py_F3z_Px = I_ERI_Gxy2z_S_F3z_Px+ABY*I_ERI_Fx2z_S_F3z_Px;
  Double I_ERI_F3y_Py_F3z_Px = I_ERI_G4y_S_F3z_Px+ABY*I_ERI_F3y_S_F3z_Px;
  Double I_ERI_F2yz_Py_F3z_Px = I_ERI_G3yz_S_F3z_Px+ABY*I_ERI_F2yz_S_F3z_Px;
  Double I_ERI_Fy2z_Py_F3z_Px = I_ERI_G2y2z_S_F3z_Px+ABY*I_ERI_Fy2z_S_F3z_Px;
  Double I_ERI_F3z_Py_F3z_Px = I_ERI_Gy3z_S_F3z_Px+ABY*I_ERI_F3z_S_F3z_Px;
  Double I_ERI_F2xz_Pz_F3z_Px = I_ERI_G2x2z_S_F3z_Px+ABZ*I_ERI_F2xz_S_F3z_Px;
  Double I_ERI_Fxyz_Pz_F3z_Px = I_ERI_Gxy2z_S_F3z_Px+ABZ*I_ERI_Fxyz_S_F3z_Px;
  Double I_ERI_Fx2z_Pz_F3z_Px = I_ERI_Gx3z_S_F3z_Px+ABZ*I_ERI_Fx2z_S_F3z_Px;
  Double I_ERI_F2yz_Pz_F3z_Px = I_ERI_G2y2z_S_F3z_Px+ABZ*I_ERI_F2yz_S_F3z_Px;
  Double I_ERI_Fy2z_Pz_F3z_Px = I_ERI_Gy3z_S_F3z_Px+ABZ*I_ERI_Fy2z_S_F3z_Px;
  Double I_ERI_F3z_Pz_F3z_Px = I_ERI_G4z_S_F3z_Px+ABZ*I_ERI_F3z_S_F3z_Px;
  Double I_ERI_F3x_Px_F3x_Py = I_ERI_G4x_S_F3x_Py+ABX*I_ERI_F3x_S_F3x_Py;
  Double I_ERI_F2xy_Px_F3x_Py = I_ERI_G3xy_S_F3x_Py+ABX*I_ERI_F2xy_S_F3x_Py;
  Double I_ERI_F2xz_Px_F3x_Py = I_ERI_G3xz_S_F3x_Py+ABX*I_ERI_F2xz_S_F3x_Py;
  Double I_ERI_Fx2y_Px_F3x_Py = I_ERI_G2x2y_S_F3x_Py+ABX*I_ERI_Fx2y_S_F3x_Py;
  Double I_ERI_Fxyz_Px_F3x_Py = I_ERI_G2xyz_S_F3x_Py+ABX*I_ERI_Fxyz_S_F3x_Py;
  Double I_ERI_Fx2z_Px_F3x_Py = I_ERI_G2x2z_S_F3x_Py+ABX*I_ERI_Fx2z_S_F3x_Py;
  Double I_ERI_F3y_Px_F3x_Py = I_ERI_Gx3y_S_F3x_Py+ABX*I_ERI_F3y_S_F3x_Py;
  Double I_ERI_F2yz_Px_F3x_Py = I_ERI_Gx2yz_S_F3x_Py+ABX*I_ERI_F2yz_S_F3x_Py;
  Double I_ERI_Fy2z_Px_F3x_Py = I_ERI_Gxy2z_S_F3x_Py+ABX*I_ERI_Fy2z_S_F3x_Py;
  Double I_ERI_F3z_Px_F3x_Py = I_ERI_Gx3z_S_F3x_Py+ABX*I_ERI_F3z_S_F3x_Py;
  Double I_ERI_F2xy_Py_F3x_Py = I_ERI_G2x2y_S_F3x_Py+ABY*I_ERI_F2xy_S_F3x_Py;
  Double I_ERI_F2xz_Py_F3x_Py = I_ERI_G2xyz_S_F3x_Py+ABY*I_ERI_F2xz_S_F3x_Py;
  Double I_ERI_Fx2y_Py_F3x_Py = I_ERI_Gx3y_S_F3x_Py+ABY*I_ERI_Fx2y_S_F3x_Py;
  Double I_ERI_Fxyz_Py_F3x_Py = I_ERI_Gx2yz_S_F3x_Py+ABY*I_ERI_Fxyz_S_F3x_Py;
  Double I_ERI_Fx2z_Py_F3x_Py = I_ERI_Gxy2z_S_F3x_Py+ABY*I_ERI_Fx2z_S_F3x_Py;
  Double I_ERI_F3y_Py_F3x_Py = I_ERI_G4y_S_F3x_Py+ABY*I_ERI_F3y_S_F3x_Py;
  Double I_ERI_F2yz_Py_F3x_Py = I_ERI_G3yz_S_F3x_Py+ABY*I_ERI_F2yz_S_F3x_Py;
  Double I_ERI_Fy2z_Py_F3x_Py = I_ERI_G2y2z_S_F3x_Py+ABY*I_ERI_Fy2z_S_F3x_Py;
  Double I_ERI_F3z_Py_F3x_Py = I_ERI_Gy3z_S_F3x_Py+ABY*I_ERI_F3z_S_F3x_Py;
  Double I_ERI_F2xz_Pz_F3x_Py = I_ERI_G2x2z_S_F3x_Py+ABZ*I_ERI_F2xz_S_F3x_Py;
  Double I_ERI_Fxyz_Pz_F3x_Py = I_ERI_Gxy2z_S_F3x_Py+ABZ*I_ERI_Fxyz_S_F3x_Py;
  Double I_ERI_Fx2z_Pz_F3x_Py = I_ERI_Gx3z_S_F3x_Py+ABZ*I_ERI_Fx2z_S_F3x_Py;
  Double I_ERI_F2yz_Pz_F3x_Py = I_ERI_G2y2z_S_F3x_Py+ABZ*I_ERI_F2yz_S_F3x_Py;
  Double I_ERI_Fy2z_Pz_F3x_Py = I_ERI_Gy3z_S_F3x_Py+ABZ*I_ERI_Fy2z_S_F3x_Py;
  Double I_ERI_F3z_Pz_F3x_Py = I_ERI_G4z_S_F3x_Py+ABZ*I_ERI_F3z_S_F3x_Py;
  Double I_ERI_F3x_Px_F2xy_Py = I_ERI_G4x_S_F2xy_Py+ABX*I_ERI_F3x_S_F2xy_Py;
  Double I_ERI_F2xy_Px_F2xy_Py = I_ERI_G3xy_S_F2xy_Py+ABX*I_ERI_F2xy_S_F2xy_Py;
  Double I_ERI_F2xz_Px_F2xy_Py = I_ERI_G3xz_S_F2xy_Py+ABX*I_ERI_F2xz_S_F2xy_Py;
  Double I_ERI_Fx2y_Px_F2xy_Py = I_ERI_G2x2y_S_F2xy_Py+ABX*I_ERI_Fx2y_S_F2xy_Py;
  Double I_ERI_Fxyz_Px_F2xy_Py = I_ERI_G2xyz_S_F2xy_Py+ABX*I_ERI_Fxyz_S_F2xy_Py;
  Double I_ERI_Fx2z_Px_F2xy_Py = I_ERI_G2x2z_S_F2xy_Py+ABX*I_ERI_Fx2z_S_F2xy_Py;
  Double I_ERI_F3y_Px_F2xy_Py = I_ERI_Gx3y_S_F2xy_Py+ABX*I_ERI_F3y_S_F2xy_Py;
  Double I_ERI_F2yz_Px_F2xy_Py = I_ERI_Gx2yz_S_F2xy_Py+ABX*I_ERI_F2yz_S_F2xy_Py;
  Double I_ERI_Fy2z_Px_F2xy_Py = I_ERI_Gxy2z_S_F2xy_Py+ABX*I_ERI_Fy2z_S_F2xy_Py;
  Double I_ERI_F3z_Px_F2xy_Py = I_ERI_Gx3z_S_F2xy_Py+ABX*I_ERI_F3z_S_F2xy_Py;
  Double I_ERI_F2xy_Py_F2xy_Py = I_ERI_G2x2y_S_F2xy_Py+ABY*I_ERI_F2xy_S_F2xy_Py;
  Double I_ERI_F2xz_Py_F2xy_Py = I_ERI_G2xyz_S_F2xy_Py+ABY*I_ERI_F2xz_S_F2xy_Py;
  Double I_ERI_Fx2y_Py_F2xy_Py = I_ERI_Gx3y_S_F2xy_Py+ABY*I_ERI_Fx2y_S_F2xy_Py;
  Double I_ERI_Fxyz_Py_F2xy_Py = I_ERI_Gx2yz_S_F2xy_Py+ABY*I_ERI_Fxyz_S_F2xy_Py;
  Double I_ERI_Fx2z_Py_F2xy_Py = I_ERI_Gxy2z_S_F2xy_Py+ABY*I_ERI_Fx2z_S_F2xy_Py;
  Double I_ERI_F3y_Py_F2xy_Py = I_ERI_G4y_S_F2xy_Py+ABY*I_ERI_F3y_S_F2xy_Py;
  Double I_ERI_F2yz_Py_F2xy_Py = I_ERI_G3yz_S_F2xy_Py+ABY*I_ERI_F2yz_S_F2xy_Py;
  Double I_ERI_Fy2z_Py_F2xy_Py = I_ERI_G2y2z_S_F2xy_Py+ABY*I_ERI_Fy2z_S_F2xy_Py;
  Double I_ERI_F3z_Py_F2xy_Py = I_ERI_Gy3z_S_F2xy_Py+ABY*I_ERI_F3z_S_F2xy_Py;
  Double I_ERI_F2xz_Pz_F2xy_Py = I_ERI_G2x2z_S_F2xy_Py+ABZ*I_ERI_F2xz_S_F2xy_Py;
  Double I_ERI_Fxyz_Pz_F2xy_Py = I_ERI_Gxy2z_S_F2xy_Py+ABZ*I_ERI_Fxyz_S_F2xy_Py;
  Double I_ERI_Fx2z_Pz_F2xy_Py = I_ERI_Gx3z_S_F2xy_Py+ABZ*I_ERI_Fx2z_S_F2xy_Py;
  Double I_ERI_F2yz_Pz_F2xy_Py = I_ERI_G2y2z_S_F2xy_Py+ABZ*I_ERI_F2yz_S_F2xy_Py;
  Double I_ERI_Fy2z_Pz_F2xy_Py = I_ERI_Gy3z_S_F2xy_Py+ABZ*I_ERI_Fy2z_S_F2xy_Py;
  Double I_ERI_F3z_Pz_F2xy_Py = I_ERI_G4z_S_F2xy_Py+ABZ*I_ERI_F3z_S_F2xy_Py;
  Double I_ERI_F3x_Px_F2xz_Py = I_ERI_G4x_S_F2xz_Py+ABX*I_ERI_F3x_S_F2xz_Py;
  Double I_ERI_F2xy_Px_F2xz_Py = I_ERI_G3xy_S_F2xz_Py+ABX*I_ERI_F2xy_S_F2xz_Py;
  Double I_ERI_F2xz_Px_F2xz_Py = I_ERI_G3xz_S_F2xz_Py+ABX*I_ERI_F2xz_S_F2xz_Py;
  Double I_ERI_Fx2y_Px_F2xz_Py = I_ERI_G2x2y_S_F2xz_Py+ABX*I_ERI_Fx2y_S_F2xz_Py;
  Double I_ERI_Fxyz_Px_F2xz_Py = I_ERI_G2xyz_S_F2xz_Py+ABX*I_ERI_Fxyz_S_F2xz_Py;
  Double I_ERI_Fx2z_Px_F2xz_Py = I_ERI_G2x2z_S_F2xz_Py+ABX*I_ERI_Fx2z_S_F2xz_Py;
  Double I_ERI_F3y_Px_F2xz_Py = I_ERI_Gx3y_S_F2xz_Py+ABX*I_ERI_F3y_S_F2xz_Py;
  Double I_ERI_F2yz_Px_F2xz_Py = I_ERI_Gx2yz_S_F2xz_Py+ABX*I_ERI_F2yz_S_F2xz_Py;
  Double I_ERI_Fy2z_Px_F2xz_Py = I_ERI_Gxy2z_S_F2xz_Py+ABX*I_ERI_Fy2z_S_F2xz_Py;
  Double I_ERI_F3z_Px_F2xz_Py = I_ERI_Gx3z_S_F2xz_Py+ABX*I_ERI_F3z_S_F2xz_Py;
  Double I_ERI_F2xy_Py_F2xz_Py = I_ERI_G2x2y_S_F2xz_Py+ABY*I_ERI_F2xy_S_F2xz_Py;
  Double I_ERI_F2xz_Py_F2xz_Py = I_ERI_G2xyz_S_F2xz_Py+ABY*I_ERI_F2xz_S_F2xz_Py;
  Double I_ERI_Fx2y_Py_F2xz_Py = I_ERI_Gx3y_S_F2xz_Py+ABY*I_ERI_Fx2y_S_F2xz_Py;
  Double I_ERI_Fxyz_Py_F2xz_Py = I_ERI_Gx2yz_S_F2xz_Py+ABY*I_ERI_Fxyz_S_F2xz_Py;
  Double I_ERI_Fx2z_Py_F2xz_Py = I_ERI_Gxy2z_S_F2xz_Py+ABY*I_ERI_Fx2z_S_F2xz_Py;
  Double I_ERI_F3y_Py_F2xz_Py = I_ERI_G4y_S_F2xz_Py+ABY*I_ERI_F3y_S_F2xz_Py;
  Double I_ERI_F2yz_Py_F2xz_Py = I_ERI_G3yz_S_F2xz_Py+ABY*I_ERI_F2yz_S_F2xz_Py;
  Double I_ERI_Fy2z_Py_F2xz_Py = I_ERI_G2y2z_S_F2xz_Py+ABY*I_ERI_Fy2z_S_F2xz_Py;
  Double I_ERI_F3z_Py_F2xz_Py = I_ERI_Gy3z_S_F2xz_Py+ABY*I_ERI_F3z_S_F2xz_Py;
  Double I_ERI_F2xz_Pz_F2xz_Py = I_ERI_G2x2z_S_F2xz_Py+ABZ*I_ERI_F2xz_S_F2xz_Py;
  Double I_ERI_Fxyz_Pz_F2xz_Py = I_ERI_Gxy2z_S_F2xz_Py+ABZ*I_ERI_Fxyz_S_F2xz_Py;
  Double I_ERI_Fx2z_Pz_F2xz_Py = I_ERI_Gx3z_S_F2xz_Py+ABZ*I_ERI_Fx2z_S_F2xz_Py;
  Double I_ERI_F2yz_Pz_F2xz_Py = I_ERI_G2y2z_S_F2xz_Py+ABZ*I_ERI_F2yz_S_F2xz_Py;
  Double I_ERI_Fy2z_Pz_F2xz_Py = I_ERI_Gy3z_S_F2xz_Py+ABZ*I_ERI_Fy2z_S_F2xz_Py;
  Double I_ERI_F3z_Pz_F2xz_Py = I_ERI_G4z_S_F2xz_Py+ABZ*I_ERI_F3z_S_F2xz_Py;
  Double I_ERI_F3x_Px_Fx2y_Py = I_ERI_G4x_S_Fx2y_Py+ABX*I_ERI_F3x_S_Fx2y_Py;
  Double I_ERI_F2xy_Px_Fx2y_Py = I_ERI_G3xy_S_Fx2y_Py+ABX*I_ERI_F2xy_S_Fx2y_Py;
  Double I_ERI_F2xz_Px_Fx2y_Py = I_ERI_G3xz_S_Fx2y_Py+ABX*I_ERI_F2xz_S_Fx2y_Py;
  Double I_ERI_Fx2y_Px_Fx2y_Py = I_ERI_G2x2y_S_Fx2y_Py+ABX*I_ERI_Fx2y_S_Fx2y_Py;
  Double I_ERI_Fxyz_Px_Fx2y_Py = I_ERI_G2xyz_S_Fx2y_Py+ABX*I_ERI_Fxyz_S_Fx2y_Py;
  Double I_ERI_Fx2z_Px_Fx2y_Py = I_ERI_G2x2z_S_Fx2y_Py+ABX*I_ERI_Fx2z_S_Fx2y_Py;
  Double I_ERI_F3y_Px_Fx2y_Py = I_ERI_Gx3y_S_Fx2y_Py+ABX*I_ERI_F3y_S_Fx2y_Py;
  Double I_ERI_F2yz_Px_Fx2y_Py = I_ERI_Gx2yz_S_Fx2y_Py+ABX*I_ERI_F2yz_S_Fx2y_Py;
  Double I_ERI_Fy2z_Px_Fx2y_Py = I_ERI_Gxy2z_S_Fx2y_Py+ABX*I_ERI_Fy2z_S_Fx2y_Py;
  Double I_ERI_F3z_Px_Fx2y_Py = I_ERI_Gx3z_S_Fx2y_Py+ABX*I_ERI_F3z_S_Fx2y_Py;
  Double I_ERI_F2xy_Py_Fx2y_Py = I_ERI_G2x2y_S_Fx2y_Py+ABY*I_ERI_F2xy_S_Fx2y_Py;
  Double I_ERI_F2xz_Py_Fx2y_Py = I_ERI_G2xyz_S_Fx2y_Py+ABY*I_ERI_F2xz_S_Fx2y_Py;
  Double I_ERI_Fx2y_Py_Fx2y_Py = I_ERI_Gx3y_S_Fx2y_Py+ABY*I_ERI_Fx2y_S_Fx2y_Py;
  Double I_ERI_Fxyz_Py_Fx2y_Py = I_ERI_Gx2yz_S_Fx2y_Py+ABY*I_ERI_Fxyz_S_Fx2y_Py;
  Double I_ERI_Fx2z_Py_Fx2y_Py = I_ERI_Gxy2z_S_Fx2y_Py+ABY*I_ERI_Fx2z_S_Fx2y_Py;
  Double I_ERI_F3y_Py_Fx2y_Py = I_ERI_G4y_S_Fx2y_Py+ABY*I_ERI_F3y_S_Fx2y_Py;
  Double I_ERI_F2yz_Py_Fx2y_Py = I_ERI_G3yz_S_Fx2y_Py+ABY*I_ERI_F2yz_S_Fx2y_Py;
  Double I_ERI_Fy2z_Py_Fx2y_Py = I_ERI_G2y2z_S_Fx2y_Py+ABY*I_ERI_Fy2z_S_Fx2y_Py;
  Double I_ERI_F3z_Py_Fx2y_Py = I_ERI_Gy3z_S_Fx2y_Py+ABY*I_ERI_F3z_S_Fx2y_Py;
  Double I_ERI_F2xz_Pz_Fx2y_Py = I_ERI_G2x2z_S_Fx2y_Py+ABZ*I_ERI_F2xz_S_Fx2y_Py;
  Double I_ERI_Fxyz_Pz_Fx2y_Py = I_ERI_Gxy2z_S_Fx2y_Py+ABZ*I_ERI_Fxyz_S_Fx2y_Py;
  Double I_ERI_Fx2z_Pz_Fx2y_Py = I_ERI_Gx3z_S_Fx2y_Py+ABZ*I_ERI_Fx2z_S_Fx2y_Py;
  Double I_ERI_F2yz_Pz_Fx2y_Py = I_ERI_G2y2z_S_Fx2y_Py+ABZ*I_ERI_F2yz_S_Fx2y_Py;
  Double I_ERI_Fy2z_Pz_Fx2y_Py = I_ERI_Gy3z_S_Fx2y_Py+ABZ*I_ERI_Fy2z_S_Fx2y_Py;
  Double I_ERI_F3z_Pz_Fx2y_Py = I_ERI_G4z_S_Fx2y_Py+ABZ*I_ERI_F3z_S_Fx2y_Py;
  Double I_ERI_F3x_Px_Fxyz_Py = I_ERI_G4x_S_Fxyz_Py+ABX*I_ERI_F3x_S_Fxyz_Py;
  Double I_ERI_F2xy_Px_Fxyz_Py = I_ERI_G3xy_S_Fxyz_Py+ABX*I_ERI_F2xy_S_Fxyz_Py;
  Double I_ERI_F2xz_Px_Fxyz_Py = I_ERI_G3xz_S_Fxyz_Py+ABX*I_ERI_F2xz_S_Fxyz_Py;
  Double I_ERI_Fx2y_Px_Fxyz_Py = I_ERI_G2x2y_S_Fxyz_Py+ABX*I_ERI_Fx2y_S_Fxyz_Py;
  Double I_ERI_Fxyz_Px_Fxyz_Py = I_ERI_G2xyz_S_Fxyz_Py+ABX*I_ERI_Fxyz_S_Fxyz_Py;
  Double I_ERI_Fx2z_Px_Fxyz_Py = I_ERI_G2x2z_S_Fxyz_Py+ABX*I_ERI_Fx2z_S_Fxyz_Py;
  Double I_ERI_F3y_Px_Fxyz_Py = I_ERI_Gx3y_S_Fxyz_Py+ABX*I_ERI_F3y_S_Fxyz_Py;
  Double I_ERI_F2yz_Px_Fxyz_Py = I_ERI_Gx2yz_S_Fxyz_Py+ABX*I_ERI_F2yz_S_Fxyz_Py;
  Double I_ERI_Fy2z_Px_Fxyz_Py = I_ERI_Gxy2z_S_Fxyz_Py+ABX*I_ERI_Fy2z_S_Fxyz_Py;
  Double I_ERI_F3z_Px_Fxyz_Py = I_ERI_Gx3z_S_Fxyz_Py+ABX*I_ERI_F3z_S_Fxyz_Py;
  Double I_ERI_F2xy_Py_Fxyz_Py = I_ERI_G2x2y_S_Fxyz_Py+ABY*I_ERI_F2xy_S_Fxyz_Py;
  Double I_ERI_F2xz_Py_Fxyz_Py = I_ERI_G2xyz_S_Fxyz_Py+ABY*I_ERI_F2xz_S_Fxyz_Py;
  Double I_ERI_Fx2y_Py_Fxyz_Py = I_ERI_Gx3y_S_Fxyz_Py+ABY*I_ERI_Fx2y_S_Fxyz_Py;
  Double I_ERI_Fxyz_Py_Fxyz_Py = I_ERI_Gx2yz_S_Fxyz_Py+ABY*I_ERI_Fxyz_S_Fxyz_Py;
  Double I_ERI_Fx2z_Py_Fxyz_Py = I_ERI_Gxy2z_S_Fxyz_Py+ABY*I_ERI_Fx2z_S_Fxyz_Py;
  Double I_ERI_F3y_Py_Fxyz_Py = I_ERI_G4y_S_Fxyz_Py+ABY*I_ERI_F3y_S_Fxyz_Py;
  Double I_ERI_F2yz_Py_Fxyz_Py = I_ERI_G3yz_S_Fxyz_Py+ABY*I_ERI_F2yz_S_Fxyz_Py;
  Double I_ERI_Fy2z_Py_Fxyz_Py = I_ERI_G2y2z_S_Fxyz_Py+ABY*I_ERI_Fy2z_S_Fxyz_Py;
  Double I_ERI_F3z_Py_Fxyz_Py = I_ERI_Gy3z_S_Fxyz_Py+ABY*I_ERI_F3z_S_Fxyz_Py;
  Double I_ERI_F2xz_Pz_Fxyz_Py = I_ERI_G2x2z_S_Fxyz_Py+ABZ*I_ERI_F2xz_S_Fxyz_Py;
  Double I_ERI_Fxyz_Pz_Fxyz_Py = I_ERI_Gxy2z_S_Fxyz_Py+ABZ*I_ERI_Fxyz_S_Fxyz_Py;
  Double I_ERI_Fx2z_Pz_Fxyz_Py = I_ERI_Gx3z_S_Fxyz_Py+ABZ*I_ERI_Fx2z_S_Fxyz_Py;
  Double I_ERI_F2yz_Pz_Fxyz_Py = I_ERI_G2y2z_S_Fxyz_Py+ABZ*I_ERI_F2yz_S_Fxyz_Py;
  Double I_ERI_Fy2z_Pz_Fxyz_Py = I_ERI_Gy3z_S_Fxyz_Py+ABZ*I_ERI_Fy2z_S_Fxyz_Py;
  Double I_ERI_F3z_Pz_Fxyz_Py = I_ERI_G4z_S_Fxyz_Py+ABZ*I_ERI_F3z_S_Fxyz_Py;
  Double I_ERI_F3x_Px_Fx2z_Py = I_ERI_G4x_S_Fx2z_Py+ABX*I_ERI_F3x_S_Fx2z_Py;
  Double I_ERI_F2xy_Px_Fx2z_Py = I_ERI_G3xy_S_Fx2z_Py+ABX*I_ERI_F2xy_S_Fx2z_Py;
  Double I_ERI_F2xz_Px_Fx2z_Py = I_ERI_G3xz_S_Fx2z_Py+ABX*I_ERI_F2xz_S_Fx2z_Py;
  Double I_ERI_Fx2y_Px_Fx2z_Py = I_ERI_G2x2y_S_Fx2z_Py+ABX*I_ERI_Fx2y_S_Fx2z_Py;
  Double I_ERI_Fxyz_Px_Fx2z_Py = I_ERI_G2xyz_S_Fx2z_Py+ABX*I_ERI_Fxyz_S_Fx2z_Py;
  Double I_ERI_Fx2z_Px_Fx2z_Py = I_ERI_G2x2z_S_Fx2z_Py+ABX*I_ERI_Fx2z_S_Fx2z_Py;
  Double I_ERI_F3y_Px_Fx2z_Py = I_ERI_Gx3y_S_Fx2z_Py+ABX*I_ERI_F3y_S_Fx2z_Py;
  Double I_ERI_F2yz_Px_Fx2z_Py = I_ERI_Gx2yz_S_Fx2z_Py+ABX*I_ERI_F2yz_S_Fx2z_Py;
  Double I_ERI_Fy2z_Px_Fx2z_Py = I_ERI_Gxy2z_S_Fx2z_Py+ABX*I_ERI_Fy2z_S_Fx2z_Py;
  Double I_ERI_F3z_Px_Fx2z_Py = I_ERI_Gx3z_S_Fx2z_Py+ABX*I_ERI_F3z_S_Fx2z_Py;
  Double I_ERI_F2xy_Py_Fx2z_Py = I_ERI_G2x2y_S_Fx2z_Py+ABY*I_ERI_F2xy_S_Fx2z_Py;
  Double I_ERI_F2xz_Py_Fx2z_Py = I_ERI_G2xyz_S_Fx2z_Py+ABY*I_ERI_F2xz_S_Fx2z_Py;
  Double I_ERI_Fx2y_Py_Fx2z_Py = I_ERI_Gx3y_S_Fx2z_Py+ABY*I_ERI_Fx2y_S_Fx2z_Py;
  Double I_ERI_Fxyz_Py_Fx2z_Py = I_ERI_Gx2yz_S_Fx2z_Py+ABY*I_ERI_Fxyz_S_Fx2z_Py;
  Double I_ERI_Fx2z_Py_Fx2z_Py = I_ERI_Gxy2z_S_Fx2z_Py+ABY*I_ERI_Fx2z_S_Fx2z_Py;
  Double I_ERI_F3y_Py_Fx2z_Py = I_ERI_G4y_S_Fx2z_Py+ABY*I_ERI_F3y_S_Fx2z_Py;
  Double I_ERI_F2yz_Py_Fx2z_Py = I_ERI_G3yz_S_Fx2z_Py+ABY*I_ERI_F2yz_S_Fx2z_Py;
  Double I_ERI_Fy2z_Py_Fx2z_Py = I_ERI_G2y2z_S_Fx2z_Py+ABY*I_ERI_Fy2z_S_Fx2z_Py;
  Double I_ERI_F3z_Py_Fx2z_Py = I_ERI_Gy3z_S_Fx2z_Py+ABY*I_ERI_F3z_S_Fx2z_Py;
  Double I_ERI_F2xz_Pz_Fx2z_Py = I_ERI_G2x2z_S_Fx2z_Py+ABZ*I_ERI_F2xz_S_Fx2z_Py;
  Double I_ERI_Fxyz_Pz_Fx2z_Py = I_ERI_Gxy2z_S_Fx2z_Py+ABZ*I_ERI_Fxyz_S_Fx2z_Py;
  Double I_ERI_Fx2z_Pz_Fx2z_Py = I_ERI_Gx3z_S_Fx2z_Py+ABZ*I_ERI_Fx2z_S_Fx2z_Py;
  Double I_ERI_F2yz_Pz_Fx2z_Py = I_ERI_G2y2z_S_Fx2z_Py+ABZ*I_ERI_F2yz_S_Fx2z_Py;
  Double I_ERI_Fy2z_Pz_Fx2z_Py = I_ERI_Gy3z_S_Fx2z_Py+ABZ*I_ERI_Fy2z_S_Fx2z_Py;
  Double I_ERI_F3z_Pz_Fx2z_Py = I_ERI_G4z_S_Fx2z_Py+ABZ*I_ERI_F3z_S_Fx2z_Py;
  Double I_ERI_F3x_Px_F3y_Py = I_ERI_G4x_S_F3y_Py+ABX*I_ERI_F3x_S_F3y_Py;
  Double I_ERI_F2xy_Px_F3y_Py = I_ERI_G3xy_S_F3y_Py+ABX*I_ERI_F2xy_S_F3y_Py;
  Double I_ERI_F2xz_Px_F3y_Py = I_ERI_G3xz_S_F3y_Py+ABX*I_ERI_F2xz_S_F3y_Py;
  Double I_ERI_Fx2y_Px_F3y_Py = I_ERI_G2x2y_S_F3y_Py+ABX*I_ERI_Fx2y_S_F3y_Py;
  Double I_ERI_Fxyz_Px_F3y_Py = I_ERI_G2xyz_S_F3y_Py+ABX*I_ERI_Fxyz_S_F3y_Py;
  Double I_ERI_Fx2z_Px_F3y_Py = I_ERI_G2x2z_S_F3y_Py+ABX*I_ERI_Fx2z_S_F3y_Py;
  Double I_ERI_F3y_Px_F3y_Py = I_ERI_Gx3y_S_F3y_Py+ABX*I_ERI_F3y_S_F3y_Py;
  Double I_ERI_F2yz_Px_F3y_Py = I_ERI_Gx2yz_S_F3y_Py+ABX*I_ERI_F2yz_S_F3y_Py;
  Double I_ERI_Fy2z_Px_F3y_Py = I_ERI_Gxy2z_S_F3y_Py+ABX*I_ERI_Fy2z_S_F3y_Py;
  Double I_ERI_F3z_Px_F3y_Py = I_ERI_Gx3z_S_F3y_Py+ABX*I_ERI_F3z_S_F3y_Py;
  Double I_ERI_F2xy_Py_F3y_Py = I_ERI_G2x2y_S_F3y_Py+ABY*I_ERI_F2xy_S_F3y_Py;
  Double I_ERI_F2xz_Py_F3y_Py = I_ERI_G2xyz_S_F3y_Py+ABY*I_ERI_F2xz_S_F3y_Py;
  Double I_ERI_Fx2y_Py_F3y_Py = I_ERI_Gx3y_S_F3y_Py+ABY*I_ERI_Fx2y_S_F3y_Py;
  Double I_ERI_Fxyz_Py_F3y_Py = I_ERI_Gx2yz_S_F3y_Py+ABY*I_ERI_Fxyz_S_F3y_Py;
  Double I_ERI_Fx2z_Py_F3y_Py = I_ERI_Gxy2z_S_F3y_Py+ABY*I_ERI_Fx2z_S_F3y_Py;
  Double I_ERI_F3y_Py_F3y_Py = I_ERI_G4y_S_F3y_Py+ABY*I_ERI_F3y_S_F3y_Py;
  Double I_ERI_F2yz_Py_F3y_Py = I_ERI_G3yz_S_F3y_Py+ABY*I_ERI_F2yz_S_F3y_Py;
  Double I_ERI_Fy2z_Py_F3y_Py = I_ERI_G2y2z_S_F3y_Py+ABY*I_ERI_Fy2z_S_F3y_Py;
  Double I_ERI_F3z_Py_F3y_Py = I_ERI_Gy3z_S_F3y_Py+ABY*I_ERI_F3z_S_F3y_Py;
  Double I_ERI_F2xz_Pz_F3y_Py = I_ERI_G2x2z_S_F3y_Py+ABZ*I_ERI_F2xz_S_F3y_Py;
  Double I_ERI_Fxyz_Pz_F3y_Py = I_ERI_Gxy2z_S_F3y_Py+ABZ*I_ERI_Fxyz_S_F3y_Py;
  Double I_ERI_Fx2z_Pz_F3y_Py = I_ERI_Gx3z_S_F3y_Py+ABZ*I_ERI_Fx2z_S_F3y_Py;
  Double I_ERI_F2yz_Pz_F3y_Py = I_ERI_G2y2z_S_F3y_Py+ABZ*I_ERI_F2yz_S_F3y_Py;
  Double I_ERI_Fy2z_Pz_F3y_Py = I_ERI_Gy3z_S_F3y_Py+ABZ*I_ERI_Fy2z_S_F3y_Py;
  Double I_ERI_F3z_Pz_F3y_Py = I_ERI_G4z_S_F3y_Py+ABZ*I_ERI_F3z_S_F3y_Py;
  Double I_ERI_F3x_Px_F2yz_Py = I_ERI_G4x_S_F2yz_Py+ABX*I_ERI_F3x_S_F2yz_Py;
  Double I_ERI_F2xy_Px_F2yz_Py = I_ERI_G3xy_S_F2yz_Py+ABX*I_ERI_F2xy_S_F2yz_Py;
  Double I_ERI_F2xz_Px_F2yz_Py = I_ERI_G3xz_S_F2yz_Py+ABX*I_ERI_F2xz_S_F2yz_Py;
  Double I_ERI_Fx2y_Px_F2yz_Py = I_ERI_G2x2y_S_F2yz_Py+ABX*I_ERI_Fx2y_S_F2yz_Py;
  Double I_ERI_Fxyz_Px_F2yz_Py = I_ERI_G2xyz_S_F2yz_Py+ABX*I_ERI_Fxyz_S_F2yz_Py;
  Double I_ERI_Fx2z_Px_F2yz_Py = I_ERI_G2x2z_S_F2yz_Py+ABX*I_ERI_Fx2z_S_F2yz_Py;
  Double I_ERI_F3y_Px_F2yz_Py = I_ERI_Gx3y_S_F2yz_Py+ABX*I_ERI_F3y_S_F2yz_Py;
  Double I_ERI_F2yz_Px_F2yz_Py = I_ERI_Gx2yz_S_F2yz_Py+ABX*I_ERI_F2yz_S_F2yz_Py;
  Double I_ERI_Fy2z_Px_F2yz_Py = I_ERI_Gxy2z_S_F2yz_Py+ABX*I_ERI_Fy2z_S_F2yz_Py;
  Double I_ERI_F3z_Px_F2yz_Py = I_ERI_Gx3z_S_F2yz_Py+ABX*I_ERI_F3z_S_F2yz_Py;
  Double I_ERI_F2xy_Py_F2yz_Py = I_ERI_G2x2y_S_F2yz_Py+ABY*I_ERI_F2xy_S_F2yz_Py;
  Double I_ERI_F2xz_Py_F2yz_Py = I_ERI_G2xyz_S_F2yz_Py+ABY*I_ERI_F2xz_S_F2yz_Py;
  Double I_ERI_Fx2y_Py_F2yz_Py = I_ERI_Gx3y_S_F2yz_Py+ABY*I_ERI_Fx2y_S_F2yz_Py;
  Double I_ERI_Fxyz_Py_F2yz_Py = I_ERI_Gx2yz_S_F2yz_Py+ABY*I_ERI_Fxyz_S_F2yz_Py;
  Double I_ERI_Fx2z_Py_F2yz_Py = I_ERI_Gxy2z_S_F2yz_Py+ABY*I_ERI_Fx2z_S_F2yz_Py;
  Double I_ERI_F3y_Py_F2yz_Py = I_ERI_G4y_S_F2yz_Py+ABY*I_ERI_F3y_S_F2yz_Py;
  Double I_ERI_F2yz_Py_F2yz_Py = I_ERI_G3yz_S_F2yz_Py+ABY*I_ERI_F2yz_S_F2yz_Py;
  Double I_ERI_Fy2z_Py_F2yz_Py = I_ERI_G2y2z_S_F2yz_Py+ABY*I_ERI_Fy2z_S_F2yz_Py;
  Double I_ERI_F3z_Py_F2yz_Py = I_ERI_Gy3z_S_F2yz_Py+ABY*I_ERI_F3z_S_F2yz_Py;
  Double I_ERI_F2xz_Pz_F2yz_Py = I_ERI_G2x2z_S_F2yz_Py+ABZ*I_ERI_F2xz_S_F2yz_Py;
  Double I_ERI_Fxyz_Pz_F2yz_Py = I_ERI_Gxy2z_S_F2yz_Py+ABZ*I_ERI_Fxyz_S_F2yz_Py;
  Double I_ERI_Fx2z_Pz_F2yz_Py = I_ERI_Gx3z_S_F2yz_Py+ABZ*I_ERI_Fx2z_S_F2yz_Py;
  Double I_ERI_F2yz_Pz_F2yz_Py = I_ERI_G2y2z_S_F2yz_Py+ABZ*I_ERI_F2yz_S_F2yz_Py;
  Double I_ERI_Fy2z_Pz_F2yz_Py = I_ERI_Gy3z_S_F2yz_Py+ABZ*I_ERI_Fy2z_S_F2yz_Py;
  Double I_ERI_F3z_Pz_F2yz_Py = I_ERI_G4z_S_F2yz_Py+ABZ*I_ERI_F3z_S_F2yz_Py;
  Double I_ERI_F3x_Px_Fy2z_Py = I_ERI_G4x_S_Fy2z_Py+ABX*I_ERI_F3x_S_Fy2z_Py;
  Double I_ERI_F2xy_Px_Fy2z_Py = I_ERI_G3xy_S_Fy2z_Py+ABX*I_ERI_F2xy_S_Fy2z_Py;
  Double I_ERI_F2xz_Px_Fy2z_Py = I_ERI_G3xz_S_Fy2z_Py+ABX*I_ERI_F2xz_S_Fy2z_Py;
  Double I_ERI_Fx2y_Px_Fy2z_Py = I_ERI_G2x2y_S_Fy2z_Py+ABX*I_ERI_Fx2y_S_Fy2z_Py;
  Double I_ERI_Fxyz_Px_Fy2z_Py = I_ERI_G2xyz_S_Fy2z_Py+ABX*I_ERI_Fxyz_S_Fy2z_Py;
  Double I_ERI_Fx2z_Px_Fy2z_Py = I_ERI_G2x2z_S_Fy2z_Py+ABX*I_ERI_Fx2z_S_Fy2z_Py;
  Double I_ERI_F3y_Px_Fy2z_Py = I_ERI_Gx3y_S_Fy2z_Py+ABX*I_ERI_F3y_S_Fy2z_Py;
  Double I_ERI_F2yz_Px_Fy2z_Py = I_ERI_Gx2yz_S_Fy2z_Py+ABX*I_ERI_F2yz_S_Fy2z_Py;
  Double I_ERI_Fy2z_Px_Fy2z_Py = I_ERI_Gxy2z_S_Fy2z_Py+ABX*I_ERI_Fy2z_S_Fy2z_Py;
  Double I_ERI_F3z_Px_Fy2z_Py = I_ERI_Gx3z_S_Fy2z_Py+ABX*I_ERI_F3z_S_Fy2z_Py;
  Double I_ERI_F2xy_Py_Fy2z_Py = I_ERI_G2x2y_S_Fy2z_Py+ABY*I_ERI_F2xy_S_Fy2z_Py;
  Double I_ERI_F2xz_Py_Fy2z_Py = I_ERI_G2xyz_S_Fy2z_Py+ABY*I_ERI_F2xz_S_Fy2z_Py;
  Double I_ERI_Fx2y_Py_Fy2z_Py = I_ERI_Gx3y_S_Fy2z_Py+ABY*I_ERI_Fx2y_S_Fy2z_Py;
  Double I_ERI_Fxyz_Py_Fy2z_Py = I_ERI_Gx2yz_S_Fy2z_Py+ABY*I_ERI_Fxyz_S_Fy2z_Py;
  Double I_ERI_Fx2z_Py_Fy2z_Py = I_ERI_Gxy2z_S_Fy2z_Py+ABY*I_ERI_Fx2z_S_Fy2z_Py;
  Double I_ERI_F3y_Py_Fy2z_Py = I_ERI_G4y_S_Fy2z_Py+ABY*I_ERI_F3y_S_Fy2z_Py;
  Double I_ERI_F2yz_Py_Fy2z_Py = I_ERI_G3yz_S_Fy2z_Py+ABY*I_ERI_F2yz_S_Fy2z_Py;
  Double I_ERI_Fy2z_Py_Fy2z_Py = I_ERI_G2y2z_S_Fy2z_Py+ABY*I_ERI_Fy2z_S_Fy2z_Py;
  Double I_ERI_F3z_Py_Fy2z_Py = I_ERI_Gy3z_S_Fy2z_Py+ABY*I_ERI_F3z_S_Fy2z_Py;
  Double I_ERI_F2xz_Pz_Fy2z_Py = I_ERI_G2x2z_S_Fy2z_Py+ABZ*I_ERI_F2xz_S_Fy2z_Py;
  Double I_ERI_Fxyz_Pz_Fy2z_Py = I_ERI_Gxy2z_S_Fy2z_Py+ABZ*I_ERI_Fxyz_S_Fy2z_Py;
  Double I_ERI_Fx2z_Pz_Fy2z_Py = I_ERI_Gx3z_S_Fy2z_Py+ABZ*I_ERI_Fx2z_S_Fy2z_Py;
  Double I_ERI_F2yz_Pz_Fy2z_Py = I_ERI_G2y2z_S_Fy2z_Py+ABZ*I_ERI_F2yz_S_Fy2z_Py;
  Double I_ERI_Fy2z_Pz_Fy2z_Py = I_ERI_Gy3z_S_Fy2z_Py+ABZ*I_ERI_Fy2z_S_Fy2z_Py;
  Double I_ERI_F3z_Pz_Fy2z_Py = I_ERI_G4z_S_Fy2z_Py+ABZ*I_ERI_F3z_S_Fy2z_Py;
  Double I_ERI_F3x_Px_F3z_Py = I_ERI_G4x_S_F3z_Py+ABX*I_ERI_F3x_S_F3z_Py;
  Double I_ERI_F2xy_Px_F3z_Py = I_ERI_G3xy_S_F3z_Py+ABX*I_ERI_F2xy_S_F3z_Py;
  Double I_ERI_F2xz_Px_F3z_Py = I_ERI_G3xz_S_F3z_Py+ABX*I_ERI_F2xz_S_F3z_Py;
  Double I_ERI_Fx2y_Px_F3z_Py = I_ERI_G2x2y_S_F3z_Py+ABX*I_ERI_Fx2y_S_F3z_Py;
  Double I_ERI_Fxyz_Px_F3z_Py = I_ERI_G2xyz_S_F3z_Py+ABX*I_ERI_Fxyz_S_F3z_Py;
  Double I_ERI_Fx2z_Px_F3z_Py = I_ERI_G2x2z_S_F3z_Py+ABX*I_ERI_Fx2z_S_F3z_Py;
  Double I_ERI_F3y_Px_F3z_Py = I_ERI_Gx3y_S_F3z_Py+ABX*I_ERI_F3y_S_F3z_Py;
  Double I_ERI_F2yz_Px_F3z_Py = I_ERI_Gx2yz_S_F3z_Py+ABX*I_ERI_F2yz_S_F3z_Py;
  Double I_ERI_Fy2z_Px_F3z_Py = I_ERI_Gxy2z_S_F3z_Py+ABX*I_ERI_Fy2z_S_F3z_Py;
  Double I_ERI_F3z_Px_F3z_Py = I_ERI_Gx3z_S_F3z_Py+ABX*I_ERI_F3z_S_F3z_Py;
  Double I_ERI_F2xy_Py_F3z_Py = I_ERI_G2x2y_S_F3z_Py+ABY*I_ERI_F2xy_S_F3z_Py;
  Double I_ERI_F2xz_Py_F3z_Py = I_ERI_G2xyz_S_F3z_Py+ABY*I_ERI_F2xz_S_F3z_Py;
  Double I_ERI_Fx2y_Py_F3z_Py = I_ERI_Gx3y_S_F3z_Py+ABY*I_ERI_Fx2y_S_F3z_Py;
  Double I_ERI_Fxyz_Py_F3z_Py = I_ERI_Gx2yz_S_F3z_Py+ABY*I_ERI_Fxyz_S_F3z_Py;
  Double I_ERI_Fx2z_Py_F3z_Py = I_ERI_Gxy2z_S_F3z_Py+ABY*I_ERI_Fx2z_S_F3z_Py;
  Double I_ERI_F3y_Py_F3z_Py = I_ERI_G4y_S_F3z_Py+ABY*I_ERI_F3y_S_F3z_Py;
  Double I_ERI_F2yz_Py_F3z_Py = I_ERI_G3yz_S_F3z_Py+ABY*I_ERI_F2yz_S_F3z_Py;
  Double I_ERI_Fy2z_Py_F3z_Py = I_ERI_G2y2z_S_F3z_Py+ABY*I_ERI_Fy2z_S_F3z_Py;
  Double I_ERI_F3z_Py_F3z_Py = I_ERI_Gy3z_S_F3z_Py+ABY*I_ERI_F3z_S_F3z_Py;
  Double I_ERI_F2xz_Pz_F3z_Py = I_ERI_G2x2z_S_F3z_Py+ABZ*I_ERI_F2xz_S_F3z_Py;
  Double I_ERI_Fxyz_Pz_F3z_Py = I_ERI_Gxy2z_S_F3z_Py+ABZ*I_ERI_Fxyz_S_F3z_Py;
  Double I_ERI_Fx2z_Pz_F3z_Py = I_ERI_Gx3z_S_F3z_Py+ABZ*I_ERI_Fx2z_S_F3z_Py;
  Double I_ERI_F2yz_Pz_F3z_Py = I_ERI_G2y2z_S_F3z_Py+ABZ*I_ERI_F2yz_S_F3z_Py;
  Double I_ERI_Fy2z_Pz_F3z_Py = I_ERI_Gy3z_S_F3z_Py+ABZ*I_ERI_Fy2z_S_F3z_Py;
  Double I_ERI_F3z_Pz_F3z_Py = I_ERI_G4z_S_F3z_Py+ABZ*I_ERI_F3z_S_F3z_Py;
  Double I_ERI_F3x_Px_F3x_Pz = I_ERI_G4x_S_F3x_Pz+ABX*I_ERI_F3x_S_F3x_Pz;
  Double I_ERI_F2xy_Px_F3x_Pz = I_ERI_G3xy_S_F3x_Pz+ABX*I_ERI_F2xy_S_F3x_Pz;
  Double I_ERI_F2xz_Px_F3x_Pz = I_ERI_G3xz_S_F3x_Pz+ABX*I_ERI_F2xz_S_F3x_Pz;
  Double I_ERI_Fx2y_Px_F3x_Pz = I_ERI_G2x2y_S_F3x_Pz+ABX*I_ERI_Fx2y_S_F3x_Pz;
  Double I_ERI_Fxyz_Px_F3x_Pz = I_ERI_G2xyz_S_F3x_Pz+ABX*I_ERI_Fxyz_S_F3x_Pz;
  Double I_ERI_Fx2z_Px_F3x_Pz = I_ERI_G2x2z_S_F3x_Pz+ABX*I_ERI_Fx2z_S_F3x_Pz;
  Double I_ERI_F3y_Px_F3x_Pz = I_ERI_Gx3y_S_F3x_Pz+ABX*I_ERI_F3y_S_F3x_Pz;
  Double I_ERI_F2yz_Px_F3x_Pz = I_ERI_Gx2yz_S_F3x_Pz+ABX*I_ERI_F2yz_S_F3x_Pz;
  Double I_ERI_Fy2z_Px_F3x_Pz = I_ERI_Gxy2z_S_F3x_Pz+ABX*I_ERI_Fy2z_S_F3x_Pz;
  Double I_ERI_F3z_Px_F3x_Pz = I_ERI_Gx3z_S_F3x_Pz+ABX*I_ERI_F3z_S_F3x_Pz;
  Double I_ERI_F2xy_Py_F3x_Pz = I_ERI_G2x2y_S_F3x_Pz+ABY*I_ERI_F2xy_S_F3x_Pz;
  Double I_ERI_F2xz_Py_F3x_Pz = I_ERI_G2xyz_S_F3x_Pz+ABY*I_ERI_F2xz_S_F3x_Pz;
  Double I_ERI_Fx2y_Py_F3x_Pz = I_ERI_Gx3y_S_F3x_Pz+ABY*I_ERI_Fx2y_S_F3x_Pz;
  Double I_ERI_Fxyz_Py_F3x_Pz = I_ERI_Gx2yz_S_F3x_Pz+ABY*I_ERI_Fxyz_S_F3x_Pz;
  Double I_ERI_Fx2z_Py_F3x_Pz = I_ERI_Gxy2z_S_F3x_Pz+ABY*I_ERI_Fx2z_S_F3x_Pz;
  Double I_ERI_F3y_Py_F3x_Pz = I_ERI_G4y_S_F3x_Pz+ABY*I_ERI_F3y_S_F3x_Pz;
  Double I_ERI_F2yz_Py_F3x_Pz = I_ERI_G3yz_S_F3x_Pz+ABY*I_ERI_F2yz_S_F3x_Pz;
  Double I_ERI_Fy2z_Py_F3x_Pz = I_ERI_G2y2z_S_F3x_Pz+ABY*I_ERI_Fy2z_S_F3x_Pz;
  Double I_ERI_F3z_Py_F3x_Pz = I_ERI_Gy3z_S_F3x_Pz+ABY*I_ERI_F3z_S_F3x_Pz;
  Double I_ERI_F2xz_Pz_F3x_Pz = I_ERI_G2x2z_S_F3x_Pz+ABZ*I_ERI_F2xz_S_F3x_Pz;
  Double I_ERI_Fxyz_Pz_F3x_Pz = I_ERI_Gxy2z_S_F3x_Pz+ABZ*I_ERI_Fxyz_S_F3x_Pz;
  Double I_ERI_Fx2z_Pz_F3x_Pz = I_ERI_Gx3z_S_F3x_Pz+ABZ*I_ERI_Fx2z_S_F3x_Pz;
  Double I_ERI_F2yz_Pz_F3x_Pz = I_ERI_G2y2z_S_F3x_Pz+ABZ*I_ERI_F2yz_S_F3x_Pz;
  Double I_ERI_Fy2z_Pz_F3x_Pz = I_ERI_Gy3z_S_F3x_Pz+ABZ*I_ERI_Fy2z_S_F3x_Pz;
  Double I_ERI_F3z_Pz_F3x_Pz = I_ERI_G4z_S_F3x_Pz+ABZ*I_ERI_F3z_S_F3x_Pz;
  Double I_ERI_F3x_Px_F2xy_Pz = I_ERI_G4x_S_F2xy_Pz+ABX*I_ERI_F3x_S_F2xy_Pz;
  Double I_ERI_F2xy_Px_F2xy_Pz = I_ERI_G3xy_S_F2xy_Pz+ABX*I_ERI_F2xy_S_F2xy_Pz;
  Double I_ERI_F2xz_Px_F2xy_Pz = I_ERI_G3xz_S_F2xy_Pz+ABX*I_ERI_F2xz_S_F2xy_Pz;
  Double I_ERI_Fx2y_Px_F2xy_Pz = I_ERI_G2x2y_S_F2xy_Pz+ABX*I_ERI_Fx2y_S_F2xy_Pz;
  Double I_ERI_Fxyz_Px_F2xy_Pz = I_ERI_G2xyz_S_F2xy_Pz+ABX*I_ERI_Fxyz_S_F2xy_Pz;
  Double I_ERI_Fx2z_Px_F2xy_Pz = I_ERI_G2x2z_S_F2xy_Pz+ABX*I_ERI_Fx2z_S_F2xy_Pz;
  Double I_ERI_F3y_Px_F2xy_Pz = I_ERI_Gx3y_S_F2xy_Pz+ABX*I_ERI_F3y_S_F2xy_Pz;
  Double I_ERI_F2yz_Px_F2xy_Pz = I_ERI_Gx2yz_S_F2xy_Pz+ABX*I_ERI_F2yz_S_F2xy_Pz;
  Double I_ERI_Fy2z_Px_F2xy_Pz = I_ERI_Gxy2z_S_F2xy_Pz+ABX*I_ERI_Fy2z_S_F2xy_Pz;
  Double I_ERI_F3z_Px_F2xy_Pz = I_ERI_Gx3z_S_F2xy_Pz+ABX*I_ERI_F3z_S_F2xy_Pz;
  Double I_ERI_F2xy_Py_F2xy_Pz = I_ERI_G2x2y_S_F2xy_Pz+ABY*I_ERI_F2xy_S_F2xy_Pz;
  Double I_ERI_F2xz_Py_F2xy_Pz = I_ERI_G2xyz_S_F2xy_Pz+ABY*I_ERI_F2xz_S_F2xy_Pz;
  Double I_ERI_Fx2y_Py_F2xy_Pz = I_ERI_Gx3y_S_F2xy_Pz+ABY*I_ERI_Fx2y_S_F2xy_Pz;
  Double I_ERI_Fxyz_Py_F2xy_Pz = I_ERI_Gx2yz_S_F2xy_Pz+ABY*I_ERI_Fxyz_S_F2xy_Pz;
  Double I_ERI_Fx2z_Py_F2xy_Pz = I_ERI_Gxy2z_S_F2xy_Pz+ABY*I_ERI_Fx2z_S_F2xy_Pz;
  Double I_ERI_F3y_Py_F2xy_Pz = I_ERI_G4y_S_F2xy_Pz+ABY*I_ERI_F3y_S_F2xy_Pz;
  Double I_ERI_F2yz_Py_F2xy_Pz = I_ERI_G3yz_S_F2xy_Pz+ABY*I_ERI_F2yz_S_F2xy_Pz;
  Double I_ERI_Fy2z_Py_F2xy_Pz = I_ERI_G2y2z_S_F2xy_Pz+ABY*I_ERI_Fy2z_S_F2xy_Pz;
  Double I_ERI_F3z_Py_F2xy_Pz = I_ERI_Gy3z_S_F2xy_Pz+ABY*I_ERI_F3z_S_F2xy_Pz;
  Double I_ERI_F2xz_Pz_F2xy_Pz = I_ERI_G2x2z_S_F2xy_Pz+ABZ*I_ERI_F2xz_S_F2xy_Pz;
  Double I_ERI_Fxyz_Pz_F2xy_Pz = I_ERI_Gxy2z_S_F2xy_Pz+ABZ*I_ERI_Fxyz_S_F2xy_Pz;
  Double I_ERI_Fx2z_Pz_F2xy_Pz = I_ERI_Gx3z_S_F2xy_Pz+ABZ*I_ERI_Fx2z_S_F2xy_Pz;
  Double I_ERI_F2yz_Pz_F2xy_Pz = I_ERI_G2y2z_S_F2xy_Pz+ABZ*I_ERI_F2yz_S_F2xy_Pz;
  Double I_ERI_Fy2z_Pz_F2xy_Pz = I_ERI_Gy3z_S_F2xy_Pz+ABZ*I_ERI_Fy2z_S_F2xy_Pz;
  Double I_ERI_F3z_Pz_F2xy_Pz = I_ERI_G4z_S_F2xy_Pz+ABZ*I_ERI_F3z_S_F2xy_Pz;
  Double I_ERI_F3x_Px_F2xz_Pz = I_ERI_G4x_S_F2xz_Pz+ABX*I_ERI_F3x_S_F2xz_Pz;
  Double I_ERI_F2xy_Px_F2xz_Pz = I_ERI_G3xy_S_F2xz_Pz+ABX*I_ERI_F2xy_S_F2xz_Pz;
  Double I_ERI_F2xz_Px_F2xz_Pz = I_ERI_G3xz_S_F2xz_Pz+ABX*I_ERI_F2xz_S_F2xz_Pz;
  Double I_ERI_Fx2y_Px_F2xz_Pz = I_ERI_G2x2y_S_F2xz_Pz+ABX*I_ERI_Fx2y_S_F2xz_Pz;
  Double I_ERI_Fxyz_Px_F2xz_Pz = I_ERI_G2xyz_S_F2xz_Pz+ABX*I_ERI_Fxyz_S_F2xz_Pz;
  Double I_ERI_Fx2z_Px_F2xz_Pz = I_ERI_G2x2z_S_F2xz_Pz+ABX*I_ERI_Fx2z_S_F2xz_Pz;
  Double I_ERI_F3y_Px_F2xz_Pz = I_ERI_Gx3y_S_F2xz_Pz+ABX*I_ERI_F3y_S_F2xz_Pz;
  Double I_ERI_F2yz_Px_F2xz_Pz = I_ERI_Gx2yz_S_F2xz_Pz+ABX*I_ERI_F2yz_S_F2xz_Pz;
  Double I_ERI_Fy2z_Px_F2xz_Pz = I_ERI_Gxy2z_S_F2xz_Pz+ABX*I_ERI_Fy2z_S_F2xz_Pz;
  Double I_ERI_F3z_Px_F2xz_Pz = I_ERI_Gx3z_S_F2xz_Pz+ABX*I_ERI_F3z_S_F2xz_Pz;
  Double I_ERI_F2xy_Py_F2xz_Pz = I_ERI_G2x2y_S_F2xz_Pz+ABY*I_ERI_F2xy_S_F2xz_Pz;
  Double I_ERI_F2xz_Py_F2xz_Pz = I_ERI_G2xyz_S_F2xz_Pz+ABY*I_ERI_F2xz_S_F2xz_Pz;
  Double I_ERI_Fx2y_Py_F2xz_Pz = I_ERI_Gx3y_S_F2xz_Pz+ABY*I_ERI_Fx2y_S_F2xz_Pz;
  Double I_ERI_Fxyz_Py_F2xz_Pz = I_ERI_Gx2yz_S_F2xz_Pz+ABY*I_ERI_Fxyz_S_F2xz_Pz;
  Double I_ERI_Fx2z_Py_F2xz_Pz = I_ERI_Gxy2z_S_F2xz_Pz+ABY*I_ERI_Fx2z_S_F2xz_Pz;
  Double I_ERI_F3y_Py_F2xz_Pz = I_ERI_G4y_S_F2xz_Pz+ABY*I_ERI_F3y_S_F2xz_Pz;
  Double I_ERI_F2yz_Py_F2xz_Pz = I_ERI_G3yz_S_F2xz_Pz+ABY*I_ERI_F2yz_S_F2xz_Pz;
  Double I_ERI_Fy2z_Py_F2xz_Pz = I_ERI_G2y2z_S_F2xz_Pz+ABY*I_ERI_Fy2z_S_F2xz_Pz;
  Double I_ERI_F3z_Py_F2xz_Pz = I_ERI_Gy3z_S_F2xz_Pz+ABY*I_ERI_F3z_S_F2xz_Pz;
  Double I_ERI_F2xz_Pz_F2xz_Pz = I_ERI_G2x2z_S_F2xz_Pz+ABZ*I_ERI_F2xz_S_F2xz_Pz;
  Double I_ERI_Fxyz_Pz_F2xz_Pz = I_ERI_Gxy2z_S_F2xz_Pz+ABZ*I_ERI_Fxyz_S_F2xz_Pz;
  Double I_ERI_Fx2z_Pz_F2xz_Pz = I_ERI_Gx3z_S_F2xz_Pz+ABZ*I_ERI_Fx2z_S_F2xz_Pz;
  Double I_ERI_F2yz_Pz_F2xz_Pz = I_ERI_G2y2z_S_F2xz_Pz+ABZ*I_ERI_F2yz_S_F2xz_Pz;
  Double I_ERI_Fy2z_Pz_F2xz_Pz = I_ERI_Gy3z_S_F2xz_Pz+ABZ*I_ERI_Fy2z_S_F2xz_Pz;
  Double I_ERI_F3z_Pz_F2xz_Pz = I_ERI_G4z_S_F2xz_Pz+ABZ*I_ERI_F3z_S_F2xz_Pz;
  Double I_ERI_F3x_Px_Fx2y_Pz = I_ERI_G4x_S_Fx2y_Pz+ABX*I_ERI_F3x_S_Fx2y_Pz;
  Double I_ERI_F2xy_Px_Fx2y_Pz = I_ERI_G3xy_S_Fx2y_Pz+ABX*I_ERI_F2xy_S_Fx2y_Pz;
  Double I_ERI_F2xz_Px_Fx2y_Pz = I_ERI_G3xz_S_Fx2y_Pz+ABX*I_ERI_F2xz_S_Fx2y_Pz;
  Double I_ERI_Fx2y_Px_Fx2y_Pz = I_ERI_G2x2y_S_Fx2y_Pz+ABX*I_ERI_Fx2y_S_Fx2y_Pz;
  Double I_ERI_Fxyz_Px_Fx2y_Pz = I_ERI_G2xyz_S_Fx2y_Pz+ABX*I_ERI_Fxyz_S_Fx2y_Pz;
  Double I_ERI_Fx2z_Px_Fx2y_Pz = I_ERI_G2x2z_S_Fx2y_Pz+ABX*I_ERI_Fx2z_S_Fx2y_Pz;
  Double I_ERI_F3y_Px_Fx2y_Pz = I_ERI_Gx3y_S_Fx2y_Pz+ABX*I_ERI_F3y_S_Fx2y_Pz;
  Double I_ERI_F2yz_Px_Fx2y_Pz = I_ERI_Gx2yz_S_Fx2y_Pz+ABX*I_ERI_F2yz_S_Fx2y_Pz;
  Double I_ERI_Fy2z_Px_Fx2y_Pz = I_ERI_Gxy2z_S_Fx2y_Pz+ABX*I_ERI_Fy2z_S_Fx2y_Pz;
  Double I_ERI_F3z_Px_Fx2y_Pz = I_ERI_Gx3z_S_Fx2y_Pz+ABX*I_ERI_F3z_S_Fx2y_Pz;
  Double I_ERI_F2xy_Py_Fx2y_Pz = I_ERI_G2x2y_S_Fx2y_Pz+ABY*I_ERI_F2xy_S_Fx2y_Pz;
  Double I_ERI_F2xz_Py_Fx2y_Pz = I_ERI_G2xyz_S_Fx2y_Pz+ABY*I_ERI_F2xz_S_Fx2y_Pz;
  Double I_ERI_Fx2y_Py_Fx2y_Pz = I_ERI_Gx3y_S_Fx2y_Pz+ABY*I_ERI_Fx2y_S_Fx2y_Pz;
  Double I_ERI_Fxyz_Py_Fx2y_Pz = I_ERI_Gx2yz_S_Fx2y_Pz+ABY*I_ERI_Fxyz_S_Fx2y_Pz;
  Double I_ERI_Fx2z_Py_Fx2y_Pz = I_ERI_Gxy2z_S_Fx2y_Pz+ABY*I_ERI_Fx2z_S_Fx2y_Pz;
  Double I_ERI_F3y_Py_Fx2y_Pz = I_ERI_G4y_S_Fx2y_Pz+ABY*I_ERI_F3y_S_Fx2y_Pz;
  Double I_ERI_F2yz_Py_Fx2y_Pz = I_ERI_G3yz_S_Fx2y_Pz+ABY*I_ERI_F2yz_S_Fx2y_Pz;
  Double I_ERI_Fy2z_Py_Fx2y_Pz = I_ERI_G2y2z_S_Fx2y_Pz+ABY*I_ERI_Fy2z_S_Fx2y_Pz;
  Double I_ERI_F3z_Py_Fx2y_Pz = I_ERI_Gy3z_S_Fx2y_Pz+ABY*I_ERI_F3z_S_Fx2y_Pz;
  Double I_ERI_F2xz_Pz_Fx2y_Pz = I_ERI_G2x2z_S_Fx2y_Pz+ABZ*I_ERI_F2xz_S_Fx2y_Pz;
  Double I_ERI_Fxyz_Pz_Fx2y_Pz = I_ERI_Gxy2z_S_Fx2y_Pz+ABZ*I_ERI_Fxyz_S_Fx2y_Pz;
  Double I_ERI_Fx2z_Pz_Fx2y_Pz = I_ERI_Gx3z_S_Fx2y_Pz+ABZ*I_ERI_Fx2z_S_Fx2y_Pz;
  Double I_ERI_F2yz_Pz_Fx2y_Pz = I_ERI_G2y2z_S_Fx2y_Pz+ABZ*I_ERI_F2yz_S_Fx2y_Pz;
  Double I_ERI_Fy2z_Pz_Fx2y_Pz = I_ERI_Gy3z_S_Fx2y_Pz+ABZ*I_ERI_Fy2z_S_Fx2y_Pz;
  Double I_ERI_F3z_Pz_Fx2y_Pz = I_ERI_G4z_S_Fx2y_Pz+ABZ*I_ERI_F3z_S_Fx2y_Pz;
  Double I_ERI_F3x_Px_Fxyz_Pz = I_ERI_G4x_S_Fxyz_Pz+ABX*I_ERI_F3x_S_Fxyz_Pz;
  Double I_ERI_F2xy_Px_Fxyz_Pz = I_ERI_G3xy_S_Fxyz_Pz+ABX*I_ERI_F2xy_S_Fxyz_Pz;
  Double I_ERI_F2xz_Px_Fxyz_Pz = I_ERI_G3xz_S_Fxyz_Pz+ABX*I_ERI_F2xz_S_Fxyz_Pz;
  Double I_ERI_Fx2y_Px_Fxyz_Pz = I_ERI_G2x2y_S_Fxyz_Pz+ABX*I_ERI_Fx2y_S_Fxyz_Pz;
  Double I_ERI_Fxyz_Px_Fxyz_Pz = I_ERI_G2xyz_S_Fxyz_Pz+ABX*I_ERI_Fxyz_S_Fxyz_Pz;
  Double I_ERI_Fx2z_Px_Fxyz_Pz = I_ERI_G2x2z_S_Fxyz_Pz+ABX*I_ERI_Fx2z_S_Fxyz_Pz;
  Double I_ERI_F3y_Px_Fxyz_Pz = I_ERI_Gx3y_S_Fxyz_Pz+ABX*I_ERI_F3y_S_Fxyz_Pz;
  Double I_ERI_F2yz_Px_Fxyz_Pz = I_ERI_Gx2yz_S_Fxyz_Pz+ABX*I_ERI_F2yz_S_Fxyz_Pz;
  Double I_ERI_Fy2z_Px_Fxyz_Pz = I_ERI_Gxy2z_S_Fxyz_Pz+ABX*I_ERI_Fy2z_S_Fxyz_Pz;
  Double I_ERI_F3z_Px_Fxyz_Pz = I_ERI_Gx3z_S_Fxyz_Pz+ABX*I_ERI_F3z_S_Fxyz_Pz;
  Double I_ERI_F2xy_Py_Fxyz_Pz = I_ERI_G2x2y_S_Fxyz_Pz+ABY*I_ERI_F2xy_S_Fxyz_Pz;
  Double I_ERI_F2xz_Py_Fxyz_Pz = I_ERI_G2xyz_S_Fxyz_Pz+ABY*I_ERI_F2xz_S_Fxyz_Pz;
  Double I_ERI_Fx2y_Py_Fxyz_Pz = I_ERI_Gx3y_S_Fxyz_Pz+ABY*I_ERI_Fx2y_S_Fxyz_Pz;
  Double I_ERI_Fxyz_Py_Fxyz_Pz = I_ERI_Gx2yz_S_Fxyz_Pz+ABY*I_ERI_Fxyz_S_Fxyz_Pz;
  Double I_ERI_Fx2z_Py_Fxyz_Pz = I_ERI_Gxy2z_S_Fxyz_Pz+ABY*I_ERI_Fx2z_S_Fxyz_Pz;
  Double I_ERI_F3y_Py_Fxyz_Pz = I_ERI_G4y_S_Fxyz_Pz+ABY*I_ERI_F3y_S_Fxyz_Pz;
  Double I_ERI_F2yz_Py_Fxyz_Pz = I_ERI_G3yz_S_Fxyz_Pz+ABY*I_ERI_F2yz_S_Fxyz_Pz;
  Double I_ERI_Fy2z_Py_Fxyz_Pz = I_ERI_G2y2z_S_Fxyz_Pz+ABY*I_ERI_Fy2z_S_Fxyz_Pz;
  Double I_ERI_F3z_Py_Fxyz_Pz = I_ERI_Gy3z_S_Fxyz_Pz+ABY*I_ERI_F3z_S_Fxyz_Pz;
  Double I_ERI_F2xz_Pz_Fxyz_Pz = I_ERI_G2x2z_S_Fxyz_Pz+ABZ*I_ERI_F2xz_S_Fxyz_Pz;
  Double I_ERI_Fxyz_Pz_Fxyz_Pz = I_ERI_Gxy2z_S_Fxyz_Pz+ABZ*I_ERI_Fxyz_S_Fxyz_Pz;
  Double I_ERI_Fx2z_Pz_Fxyz_Pz = I_ERI_Gx3z_S_Fxyz_Pz+ABZ*I_ERI_Fx2z_S_Fxyz_Pz;
  Double I_ERI_F2yz_Pz_Fxyz_Pz = I_ERI_G2y2z_S_Fxyz_Pz+ABZ*I_ERI_F2yz_S_Fxyz_Pz;
  Double I_ERI_Fy2z_Pz_Fxyz_Pz = I_ERI_Gy3z_S_Fxyz_Pz+ABZ*I_ERI_Fy2z_S_Fxyz_Pz;
  Double I_ERI_F3z_Pz_Fxyz_Pz = I_ERI_G4z_S_Fxyz_Pz+ABZ*I_ERI_F3z_S_Fxyz_Pz;
  Double I_ERI_F3x_Px_Fx2z_Pz = I_ERI_G4x_S_Fx2z_Pz+ABX*I_ERI_F3x_S_Fx2z_Pz;
  Double I_ERI_F2xy_Px_Fx2z_Pz = I_ERI_G3xy_S_Fx2z_Pz+ABX*I_ERI_F2xy_S_Fx2z_Pz;
  Double I_ERI_F2xz_Px_Fx2z_Pz = I_ERI_G3xz_S_Fx2z_Pz+ABX*I_ERI_F2xz_S_Fx2z_Pz;
  Double I_ERI_Fx2y_Px_Fx2z_Pz = I_ERI_G2x2y_S_Fx2z_Pz+ABX*I_ERI_Fx2y_S_Fx2z_Pz;
  Double I_ERI_Fxyz_Px_Fx2z_Pz = I_ERI_G2xyz_S_Fx2z_Pz+ABX*I_ERI_Fxyz_S_Fx2z_Pz;
  Double I_ERI_Fx2z_Px_Fx2z_Pz = I_ERI_G2x2z_S_Fx2z_Pz+ABX*I_ERI_Fx2z_S_Fx2z_Pz;
  Double I_ERI_F3y_Px_Fx2z_Pz = I_ERI_Gx3y_S_Fx2z_Pz+ABX*I_ERI_F3y_S_Fx2z_Pz;
  Double I_ERI_F2yz_Px_Fx2z_Pz = I_ERI_Gx2yz_S_Fx2z_Pz+ABX*I_ERI_F2yz_S_Fx2z_Pz;
  Double I_ERI_Fy2z_Px_Fx2z_Pz = I_ERI_Gxy2z_S_Fx2z_Pz+ABX*I_ERI_Fy2z_S_Fx2z_Pz;
  Double I_ERI_F3z_Px_Fx2z_Pz = I_ERI_Gx3z_S_Fx2z_Pz+ABX*I_ERI_F3z_S_Fx2z_Pz;
  Double I_ERI_F2xy_Py_Fx2z_Pz = I_ERI_G2x2y_S_Fx2z_Pz+ABY*I_ERI_F2xy_S_Fx2z_Pz;
  Double I_ERI_F2xz_Py_Fx2z_Pz = I_ERI_G2xyz_S_Fx2z_Pz+ABY*I_ERI_F2xz_S_Fx2z_Pz;
  Double I_ERI_Fx2y_Py_Fx2z_Pz = I_ERI_Gx3y_S_Fx2z_Pz+ABY*I_ERI_Fx2y_S_Fx2z_Pz;
  Double I_ERI_Fxyz_Py_Fx2z_Pz = I_ERI_Gx2yz_S_Fx2z_Pz+ABY*I_ERI_Fxyz_S_Fx2z_Pz;
  Double I_ERI_Fx2z_Py_Fx2z_Pz = I_ERI_Gxy2z_S_Fx2z_Pz+ABY*I_ERI_Fx2z_S_Fx2z_Pz;
  Double I_ERI_F3y_Py_Fx2z_Pz = I_ERI_G4y_S_Fx2z_Pz+ABY*I_ERI_F3y_S_Fx2z_Pz;
  Double I_ERI_F2yz_Py_Fx2z_Pz = I_ERI_G3yz_S_Fx2z_Pz+ABY*I_ERI_F2yz_S_Fx2z_Pz;
  Double I_ERI_Fy2z_Py_Fx2z_Pz = I_ERI_G2y2z_S_Fx2z_Pz+ABY*I_ERI_Fy2z_S_Fx2z_Pz;
  Double I_ERI_F3z_Py_Fx2z_Pz = I_ERI_Gy3z_S_Fx2z_Pz+ABY*I_ERI_F3z_S_Fx2z_Pz;
  Double I_ERI_F2xz_Pz_Fx2z_Pz = I_ERI_G2x2z_S_Fx2z_Pz+ABZ*I_ERI_F2xz_S_Fx2z_Pz;
  Double I_ERI_Fxyz_Pz_Fx2z_Pz = I_ERI_Gxy2z_S_Fx2z_Pz+ABZ*I_ERI_Fxyz_S_Fx2z_Pz;
  Double I_ERI_Fx2z_Pz_Fx2z_Pz = I_ERI_Gx3z_S_Fx2z_Pz+ABZ*I_ERI_Fx2z_S_Fx2z_Pz;
  Double I_ERI_F2yz_Pz_Fx2z_Pz = I_ERI_G2y2z_S_Fx2z_Pz+ABZ*I_ERI_F2yz_S_Fx2z_Pz;
  Double I_ERI_Fy2z_Pz_Fx2z_Pz = I_ERI_Gy3z_S_Fx2z_Pz+ABZ*I_ERI_Fy2z_S_Fx2z_Pz;
  Double I_ERI_F3z_Pz_Fx2z_Pz = I_ERI_G4z_S_Fx2z_Pz+ABZ*I_ERI_F3z_S_Fx2z_Pz;
  Double I_ERI_F3x_Px_F3y_Pz = I_ERI_G4x_S_F3y_Pz+ABX*I_ERI_F3x_S_F3y_Pz;
  Double I_ERI_F2xy_Px_F3y_Pz = I_ERI_G3xy_S_F3y_Pz+ABX*I_ERI_F2xy_S_F3y_Pz;
  Double I_ERI_F2xz_Px_F3y_Pz = I_ERI_G3xz_S_F3y_Pz+ABX*I_ERI_F2xz_S_F3y_Pz;
  Double I_ERI_Fx2y_Px_F3y_Pz = I_ERI_G2x2y_S_F3y_Pz+ABX*I_ERI_Fx2y_S_F3y_Pz;
  Double I_ERI_Fxyz_Px_F3y_Pz = I_ERI_G2xyz_S_F3y_Pz+ABX*I_ERI_Fxyz_S_F3y_Pz;
  Double I_ERI_Fx2z_Px_F3y_Pz = I_ERI_G2x2z_S_F3y_Pz+ABX*I_ERI_Fx2z_S_F3y_Pz;
  Double I_ERI_F3y_Px_F3y_Pz = I_ERI_Gx3y_S_F3y_Pz+ABX*I_ERI_F3y_S_F3y_Pz;
  Double I_ERI_F2yz_Px_F3y_Pz = I_ERI_Gx2yz_S_F3y_Pz+ABX*I_ERI_F2yz_S_F3y_Pz;
  Double I_ERI_Fy2z_Px_F3y_Pz = I_ERI_Gxy2z_S_F3y_Pz+ABX*I_ERI_Fy2z_S_F3y_Pz;
  Double I_ERI_F3z_Px_F3y_Pz = I_ERI_Gx3z_S_F3y_Pz+ABX*I_ERI_F3z_S_F3y_Pz;
  Double I_ERI_F2xy_Py_F3y_Pz = I_ERI_G2x2y_S_F3y_Pz+ABY*I_ERI_F2xy_S_F3y_Pz;
  Double I_ERI_F2xz_Py_F3y_Pz = I_ERI_G2xyz_S_F3y_Pz+ABY*I_ERI_F2xz_S_F3y_Pz;
  Double I_ERI_Fx2y_Py_F3y_Pz = I_ERI_Gx3y_S_F3y_Pz+ABY*I_ERI_Fx2y_S_F3y_Pz;
  Double I_ERI_Fxyz_Py_F3y_Pz = I_ERI_Gx2yz_S_F3y_Pz+ABY*I_ERI_Fxyz_S_F3y_Pz;
  Double I_ERI_Fx2z_Py_F3y_Pz = I_ERI_Gxy2z_S_F3y_Pz+ABY*I_ERI_Fx2z_S_F3y_Pz;
  Double I_ERI_F3y_Py_F3y_Pz = I_ERI_G4y_S_F3y_Pz+ABY*I_ERI_F3y_S_F3y_Pz;
  Double I_ERI_F2yz_Py_F3y_Pz = I_ERI_G3yz_S_F3y_Pz+ABY*I_ERI_F2yz_S_F3y_Pz;
  Double I_ERI_Fy2z_Py_F3y_Pz = I_ERI_G2y2z_S_F3y_Pz+ABY*I_ERI_Fy2z_S_F3y_Pz;
  Double I_ERI_F3z_Py_F3y_Pz = I_ERI_Gy3z_S_F3y_Pz+ABY*I_ERI_F3z_S_F3y_Pz;
  Double I_ERI_F2xz_Pz_F3y_Pz = I_ERI_G2x2z_S_F3y_Pz+ABZ*I_ERI_F2xz_S_F3y_Pz;
  Double I_ERI_Fxyz_Pz_F3y_Pz = I_ERI_Gxy2z_S_F3y_Pz+ABZ*I_ERI_Fxyz_S_F3y_Pz;
  Double I_ERI_Fx2z_Pz_F3y_Pz = I_ERI_Gx3z_S_F3y_Pz+ABZ*I_ERI_Fx2z_S_F3y_Pz;
  Double I_ERI_F2yz_Pz_F3y_Pz = I_ERI_G2y2z_S_F3y_Pz+ABZ*I_ERI_F2yz_S_F3y_Pz;
  Double I_ERI_Fy2z_Pz_F3y_Pz = I_ERI_Gy3z_S_F3y_Pz+ABZ*I_ERI_Fy2z_S_F3y_Pz;
  Double I_ERI_F3z_Pz_F3y_Pz = I_ERI_G4z_S_F3y_Pz+ABZ*I_ERI_F3z_S_F3y_Pz;
  Double I_ERI_F3x_Px_F2yz_Pz = I_ERI_G4x_S_F2yz_Pz+ABX*I_ERI_F3x_S_F2yz_Pz;
  Double I_ERI_F2xy_Px_F2yz_Pz = I_ERI_G3xy_S_F2yz_Pz+ABX*I_ERI_F2xy_S_F2yz_Pz;
  Double I_ERI_F2xz_Px_F2yz_Pz = I_ERI_G3xz_S_F2yz_Pz+ABX*I_ERI_F2xz_S_F2yz_Pz;
  Double I_ERI_Fx2y_Px_F2yz_Pz = I_ERI_G2x2y_S_F2yz_Pz+ABX*I_ERI_Fx2y_S_F2yz_Pz;
  Double I_ERI_Fxyz_Px_F2yz_Pz = I_ERI_G2xyz_S_F2yz_Pz+ABX*I_ERI_Fxyz_S_F2yz_Pz;
  Double I_ERI_Fx2z_Px_F2yz_Pz = I_ERI_G2x2z_S_F2yz_Pz+ABX*I_ERI_Fx2z_S_F2yz_Pz;
  Double I_ERI_F3y_Px_F2yz_Pz = I_ERI_Gx3y_S_F2yz_Pz+ABX*I_ERI_F3y_S_F2yz_Pz;
  Double I_ERI_F2yz_Px_F2yz_Pz = I_ERI_Gx2yz_S_F2yz_Pz+ABX*I_ERI_F2yz_S_F2yz_Pz;
  Double I_ERI_Fy2z_Px_F2yz_Pz = I_ERI_Gxy2z_S_F2yz_Pz+ABX*I_ERI_Fy2z_S_F2yz_Pz;
  Double I_ERI_F3z_Px_F2yz_Pz = I_ERI_Gx3z_S_F2yz_Pz+ABX*I_ERI_F3z_S_F2yz_Pz;
  Double I_ERI_F2xy_Py_F2yz_Pz = I_ERI_G2x2y_S_F2yz_Pz+ABY*I_ERI_F2xy_S_F2yz_Pz;
  Double I_ERI_F2xz_Py_F2yz_Pz = I_ERI_G2xyz_S_F2yz_Pz+ABY*I_ERI_F2xz_S_F2yz_Pz;
  Double I_ERI_Fx2y_Py_F2yz_Pz = I_ERI_Gx3y_S_F2yz_Pz+ABY*I_ERI_Fx2y_S_F2yz_Pz;
  Double I_ERI_Fxyz_Py_F2yz_Pz = I_ERI_Gx2yz_S_F2yz_Pz+ABY*I_ERI_Fxyz_S_F2yz_Pz;
  Double I_ERI_Fx2z_Py_F2yz_Pz = I_ERI_Gxy2z_S_F2yz_Pz+ABY*I_ERI_Fx2z_S_F2yz_Pz;
  Double I_ERI_F3y_Py_F2yz_Pz = I_ERI_G4y_S_F2yz_Pz+ABY*I_ERI_F3y_S_F2yz_Pz;
  Double I_ERI_F2yz_Py_F2yz_Pz = I_ERI_G3yz_S_F2yz_Pz+ABY*I_ERI_F2yz_S_F2yz_Pz;
  Double I_ERI_Fy2z_Py_F2yz_Pz = I_ERI_G2y2z_S_F2yz_Pz+ABY*I_ERI_Fy2z_S_F2yz_Pz;
  Double I_ERI_F3z_Py_F2yz_Pz = I_ERI_Gy3z_S_F2yz_Pz+ABY*I_ERI_F3z_S_F2yz_Pz;
  Double I_ERI_F2xz_Pz_F2yz_Pz = I_ERI_G2x2z_S_F2yz_Pz+ABZ*I_ERI_F2xz_S_F2yz_Pz;
  Double I_ERI_Fxyz_Pz_F2yz_Pz = I_ERI_Gxy2z_S_F2yz_Pz+ABZ*I_ERI_Fxyz_S_F2yz_Pz;
  Double I_ERI_Fx2z_Pz_F2yz_Pz = I_ERI_Gx3z_S_F2yz_Pz+ABZ*I_ERI_Fx2z_S_F2yz_Pz;
  Double I_ERI_F2yz_Pz_F2yz_Pz = I_ERI_G2y2z_S_F2yz_Pz+ABZ*I_ERI_F2yz_S_F2yz_Pz;
  Double I_ERI_Fy2z_Pz_F2yz_Pz = I_ERI_Gy3z_S_F2yz_Pz+ABZ*I_ERI_Fy2z_S_F2yz_Pz;
  Double I_ERI_F3z_Pz_F2yz_Pz = I_ERI_G4z_S_F2yz_Pz+ABZ*I_ERI_F3z_S_F2yz_Pz;
  Double I_ERI_F3x_Px_Fy2z_Pz = I_ERI_G4x_S_Fy2z_Pz+ABX*I_ERI_F3x_S_Fy2z_Pz;
  Double I_ERI_F2xy_Px_Fy2z_Pz = I_ERI_G3xy_S_Fy2z_Pz+ABX*I_ERI_F2xy_S_Fy2z_Pz;
  Double I_ERI_F2xz_Px_Fy2z_Pz = I_ERI_G3xz_S_Fy2z_Pz+ABX*I_ERI_F2xz_S_Fy2z_Pz;
  Double I_ERI_Fx2y_Px_Fy2z_Pz = I_ERI_G2x2y_S_Fy2z_Pz+ABX*I_ERI_Fx2y_S_Fy2z_Pz;
  Double I_ERI_Fxyz_Px_Fy2z_Pz = I_ERI_G2xyz_S_Fy2z_Pz+ABX*I_ERI_Fxyz_S_Fy2z_Pz;
  Double I_ERI_Fx2z_Px_Fy2z_Pz = I_ERI_G2x2z_S_Fy2z_Pz+ABX*I_ERI_Fx2z_S_Fy2z_Pz;
  Double I_ERI_F3y_Px_Fy2z_Pz = I_ERI_Gx3y_S_Fy2z_Pz+ABX*I_ERI_F3y_S_Fy2z_Pz;
  Double I_ERI_F2yz_Px_Fy2z_Pz = I_ERI_Gx2yz_S_Fy2z_Pz+ABX*I_ERI_F2yz_S_Fy2z_Pz;
  Double I_ERI_Fy2z_Px_Fy2z_Pz = I_ERI_Gxy2z_S_Fy2z_Pz+ABX*I_ERI_Fy2z_S_Fy2z_Pz;
  Double I_ERI_F3z_Px_Fy2z_Pz = I_ERI_Gx3z_S_Fy2z_Pz+ABX*I_ERI_F3z_S_Fy2z_Pz;
  Double I_ERI_F2xy_Py_Fy2z_Pz = I_ERI_G2x2y_S_Fy2z_Pz+ABY*I_ERI_F2xy_S_Fy2z_Pz;
  Double I_ERI_F2xz_Py_Fy2z_Pz = I_ERI_G2xyz_S_Fy2z_Pz+ABY*I_ERI_F2xz_S_Fy2z_Pz;
  Double I_ERI_Fx2y_Py_Fy2z_Pz = I_ERI_Gx3y_S_Fy2z_Pz+ABY*I_ERI_Fx2y_S_Fy2z_Pz;
  Double I_ERI_Fxyz_Py_Fy2z_Pz = I_ERI_Gx2yz_S_Fy2z_Pz+ABY*I_ERI_Fxyz_S_Fy2z_Pz;
  Double I_ERI_Fx2z_Py_Fy2z_Pz = I_ERI_Gxy2z_S_Fy2z_Pz+ABY*I_ERI_Fx2z_S_Fy2z_Pz;
  Double I_ERI_F3y_Py_Fy2z_Pz = I_ERI_G4y_S_Fy2z_Pz+ABY*I_ERI_F3y_S_Fy2z_Pz;
  Double I_ERI_F2yz_Py_Fy2z_Pz = I_ERI_G3yz_S_Fy2z_Pz+ABY*I_ERI_F2yz_S_Fy2z_Pz;
  Double I_ERI_Fy2z_Py_Fy2z_Pz = I_ERI_G2y2z_S_Fy2z_Pz+ABY*I_ERI_Fy2z_S_Fy2z_Pz;
  Double I_ERI_F3z_Py_Fy2z_Pz = I_ERI_Gy3z_S_Fy2z_Pz+ABY*I_ERI_F3z_S_Fy2z_Pz;
  Double I_ERI_F2xz_Pz_Fy2z_Pz = I_ERI_G2x2z_S_Fy2z_Pz+ABZ*I_ERI_F2xz_S_Fy2z_Pz;
  Double I_ERI_Fxyz_Pz_Fy2z_Pz = I_ERI_Gxy2z_S_Fy2z_Pz+ABZ*I_ERI_Fxyz_S_Fy2z_Pz;
  Double I_ERI_Fx2z_Pz_Fy2z_Pz = I_ERI_Gx3z_S_Fy2z_Pz+ABZ*I_ERI_Fx2z_S_Fy2z_Pz;
  Double I_ERI_F2yz_Pz_Fy2z_Pz = I_ERI_G2y2z_S_Fy2z_Pz+ABZ*I_ERI_F2yz_S_Fy2z_Pz;
  Double I_ERI_Fy2z_Pz_Fy2z_Pz = I_ERI_Gy3z_S_Fy2z_Pz+ABZ*I_ERI_Fy2z_S_Fy2z_Pz;
  Double I_ERI_F3z_Pz_Fy2z_Pz = I_ERI_G4z_S_Fy2z_Pz+ABZ*I_ERI_F3z_S_Fy2z_Pz;
  Double I_ERI_F3x_Px_F3z_Pz = I_ERI_G4x_S_F3z_Pz+ABX*I_ERI_F3x_S_F3z_Pz;
  Double I_ERI_F2xy_Px_F3z_Pz = I_ERI_G3xy_S_F3z_Pz+ABX*I_ERI_F2xy_S_F3z_Pz;
  Double I_ERI_F2xz_Px_F3z_Pz = I_ERI_G3xz_S_F3z_Pz+ABX*I_ERI_F2xz_S_F3z_Pz;
  Double I_ERI_Fx2y_Px_F3z_Pz = I_ERI_G2x2y_S_F3z_Pz+ABX*I_ERI_Fx2y_S_F3z_Pz;
  Double I_ERI_Fxyz_Px_F3z_Pz = I_ERI_G2xyz_S_F3z_Pz+ABX*I_ERI_Fxyz_S_F3z_Pz;
  Double I_ERI_Fx2z_Px_F3z_Pz = I_ERI_G2x2z_S_F3z_Pz+ABX*I_ERI_Fx2z_S_F3z_Pz;
  Double I_ERI_F3y_Px_F3z_Pz = I_ERI_Gx3y_S_F3z_Pz+ABX*I_ERI_F3y_S_F3z_Pz;
  Double I_ERI_F2yz_Px_F3z_Pz = I_ERI_Gx2yz_S_F3z_Pz+ABX*I_ERI_F2yz_S_F3z_Pz;
  Double I_ERI_Fy2z_Px_F3z_Pz = I_ERI_Gxy2z_S_F3z_Pz+ABX*I_ERI_Fy2z_S_F3z_Pz;
  Double I_ERI_F3z_Px_F3z_Pz = I_ERI_Gx3z_S_F3z_Pz+ABX*I_ERI_F3z_S_F3z_Pz;
  Double I_ERI_F2xy_Py_F3z_Pz = I_ERI_G2x2y_S_F3z_Pz+ABY*I_ERI_F2xy_S_F3z_Pz;
  Double I_ERI_F2xz_Py_F3z_Pz = I_ERI_G2xyz_S_F3z_Pz+ABY*I_ERI_F2xz_S_F3z_Pz;
  Double I_ERI_Fx2y_Py_F3z_Pz = I_ERI_Gx3y_S_F3z_Pz+ABY*I_ERI_Fx2y_S_F3z_Pz;
  Double I_ERI_Fxyz_Py_F3z_Pz = I_ERI_Gx2yz_S_F3z_Pz+ABY*I_ERI_Fxyz_S_F3z_Pz;
  Double I_ERI_Fx2z_Py_F3z_Pz = I_ERI_Gxy2z_S_F3z_Pz+ABY*I_ERI_Fx2z_S_F3z_Pz;
  Double I_ERI_F3y_Py_F3z_Pz = I_ERI_G4y_S_F3z_Pz+ABY*I_ERI_F3y_S_F3z_Pz;
  Double I_ERI_F2yz_Py_F3z_Pz = I_ERI_G3yz_S_F3z_Pz+ABY*I_ERI_F2yz_S_F3z_Pz;
  Double I_ERI_Fy2z_Py_F3z_Pz = I_ERI_G2y2z_S_F3z_Pz+ABY*I_ERI_Fy2z_S_F3z_Pz;
  Double I_ERI_F3z_Py_F3z_Pz = I_ERI_Gy3z_S_F3z_Pz+ABY*I_ERI_F3z_S_F3z_Pz;
  Double I_ERI_F2xz_Pz_F3z_Pz = I_ERI_G2x2z_S_F3z_Pz+ABZ*I_ERI_F2xz_S_F3z_Pz;
  Double I_ERI_Fxyz_Pz_F3z_Pz = I_ERI_Gxy2z_S_F3z_Pz+ABZ*I_ERI_Fxyz_S_F3z_Pz;
  Double I_ERI_Fx2z_Pz_F3z_Pz = I_ERI_Gx3z_S_F3z_Pz+ABZ*I_ERI_Fx2z_S_F3z_Pz;
  Double I_ERI_F2yz_Pz_F3z_Pz = I_ERI_G2y2z_S_F3z_Pz+ABZ*I_ERI_F2yz_S_F3z_Pz;
  Double I_ERI_Fy2z_Pz_F3z_Pz = I_ERI_Gy3z_S_F3z_Pz+ABZ*I_ERI_Fy2z_S_F3z_Pz;
  Double I_ERI_F3z_Pz_F3z_Pz = I_ERI_G4z_S_F3z_Pz+ABZ*I_ERI_F3z_S_F3z_Pz;

  /************************************************************
   * shell quartet name: SQ_ERI_D_D_F_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_P_F_P
   * RHS shell quartet name: SQ_ERI_D_P_F_P
   ************************************************************/
  abcd[0] = I_ERI_F3x_Px_F3x_Px+ABX*I_ERI_D2x_Px_F3x_Px;
  abcd[1] = I_ERI_F2xy_Px_F3x_Px+ABX*I_ERI_Dxy_Px_F3x_Px;
  abcd[2] = I_ERI_F2xz_Px_F3x_Px+ABX*I_ERI_Dxz_Px_F3x_Px;
  abcd[3] = I_ERI_Fx2y_Px_F3x_Px+ABX*I_ERI_D2y_Px_F3x_Px;
  abcd[4] = I_ERI_Fxyz_Px_F3x_Px+ABX*I_ERI_Dyz_Px_F3x_Px;
  abcd[5] = I_ERI_Fx2z_Px_F3x_Px+ABX*I_ERI_D2z_Px_F3x_Px;
  abcd[6] = I_ERI_F2xy_Px_F3x_Px+ABY*I_ERI_D2x_Px_F3x_Px;
  abcd[7] = I_ERI_Fx2y_Px_F3x_Px+ABY*I_ERI_Dxy_Px_F3x_Px;
  abcd[8] = I_ERI_Fxyz_Px_F3x_Px+ABY*I_ERI_Dxz_Px_F3x_Px;
  abcd[9] = I_ERI_F3y_Px_F3x_Px+ABY*I_ERI_D2y_Px_F3x_Px;
  abcd[10] = I_ERI_F2yz_Px_F3x_Px+ABY*I_ERI_Dyz_Px_F3x_Px;
  abcd[11] = I_ERI_Fy2z_Px_F3x_Px+ABY*I_ERI_D2z_Px_F3x_Px;
  abcd[12] = I_ERI_F2xz_Px_F3x_Px+ABZ*I_ERI_D2x_Px_F3x_Px;
  abcd[13] = I_ERI_Fxyz_Px_F3x_Px+ABZ*I_ERI_Dxy_Px_F3x_Px;
  abcd[14] = I_ERI_Fx2z_Px_F3x_Px+ABZ*I_ERI_Dxz_Px_F3x_Px;
  abcd[15] = I_ERI_F2yz_Px_F3x_Px+ABZ*I_ERI_D2y_Px_F3x_Px;
  abcd[16] = I_ERI_Fy2z_Px_F3x_Px+ABZ*I_ERI_Dyz_Px_F3x_Px;
  abcd[17] = I_ERI_F3z_Px_F3x_Px+ABZ*I_ERI_D2z_Px_F3x_Px;
  abcd[18] = I_ERI_F2xy_Py_F3x_Px+ABY*I_ERI_D2x_Py_F3x_Px;
  abcd[19] = I_ERI_Fx2y_Py_F3x_Px+ABY*I_ERI_Dxy_Py_F3x_Px;
  abcd[20] = I_ERI_Fxyz_Py_F3x_Px+ABY*I_ERI_Dxz_Py_F3x_Px;
  abcd[21] = I_ERI_F3y_Py_F3x_Px+ABY*I_ERI_D2y_Py_F3x_Px;
  abcd[22] = I_ERI_F2yz_Py_F3x_Px+ABY*I_ERI_Dyz_Py_F3x_Px;
  abcd[23] = I_ERI_Fy2z_Py_F3x_Px+ABY*I_ERI_D2z_Py_F3x_Px;
  abcd[24] = I_ERI_F2xz_Py_F3x_Px+ABZ*I_ERI_D2x_Py_F3x_Px;
  abcd[25] = I_ERI_Fxyz_Py_F3x_Px+ABZ*I_ERI_Dxy_Py_F3x_Px;
  abcd[26] = I_ERI_Fx2z_Py_F3x_Px+ABZ*I_ERI_Dxz_Py_F3x_Px;
  abcd[27] = I_ERI_F2yz_Py_F3x_Px+ABZ*I_ERI_D2y_Py_F3x_Px;
  abcd[28] = I_ERI_Fy2z_Py_F3x_Px+ABZ*I_ERI_Dyz_Py_F3x_Px;
  abcd[29] = I_ERI_F3z_Py_F3x_Px+ABZ*I_ERI_D2z_Py_F3x_Px;
  abcd[30] = I_ERI_F2xz_Pz_F3x_Px+ABZ*I_ERI_D2x_Pz_F3x_Px;
  abcd[31] = I_ERI_Fxyz_Pz_F3x_Px+ABZ*I_ERI_Dxy_Pz_F3x_Px;
  abcd[32] = I_ERI_Fx2z_Pz_F3x_Px+ABZ*I_ERI_Dxz_Pz_F3x_Px;
  abcd[33] = I_ERI_F2yz_Pz_F3x_Px+ABZ*I_ERI_D2y_Pz_F3x_Px;
  abcd[34] = I_ERI_Fy2z_Pz_F3x_Px+ABZ*I_ERI_Dyz_Pz_F3x_Px;
  abcd[35] = I_ERI_F3z_Pz_F3x_Px+ABZ*I_ERI_D2z_Pz_F3x_Px;
  abcd[36] = I_ERI_F3x_Px_F2xy_Px+ABX*I_ERI_D2x_Px_F2xy_Px;
  abcd[37] = I_ERI_F2xy_Px_F2xy_Px+ABX*I_ERI_Dxy_Px_F2xy_Px;
  abcd[38] = I_ERI_F2xz_Px_F2xy_Px+ABX*I_ERI_Dxz_Px_F2xy_Px;
  abcd[39] = I_ERI_Fx2y_Px_F2xy_Px+ABX*I_ERI_D2y_Px_F2xy_Px;
  abcd[40] = I_ERI_Fxyz_Px_F2xy_Px+ABX*I_ERI_Dyz_Px_F2xy_Px;
  abcd[41] = I_ERI_Fx2z_Px_F2xy_Px+ABX*I_ERI_D2z_Px_F2xy_Px;
  abcd[42] = I_ERI_F2xy_Px_F2xy_Px+ABY*I_ERI_D2x_Px_F2xy_Px;
  abcd[43] = I_ERI_Fx2y_Px_F2xy_Px+ABY*I_ERI_Dxy_Px_F2xy_Px;
  abcd[44] = I_ERI_Fxyz_Px_F2xy_Px+ABY*I_ERI_Dxz_Px_F2xy_Px;
  abcd[45] = I_ERI_F3y_Px_F2xy_Px+ABY*I_ERI_D2y_Px_F2xy_Px;
  abcd[46] = I_ERI_F2yz_Px_F2xy_Px+ABY*I_ERI_Dyz_Px_F2xy_Px;
  abcd[47] = I_ERI_Fy2z_Px_F2xy_Px+ABY*I_ERI_D2z_Px_F2xy_Px;
  abcd[48] = I_ERI_F2xz_Px_F2xy_Px+ABZ*I_ERI_D2x_Px_F2xy_Px;
  abcd[49] = I_ERI_Fxyz_Px_F2xy_Px+ABZ*I_ERI_Dxy_Px_F2xy_Px;
  abcd[50] = I_ERI_Fx2z_Px_F2xy_Px+ABZ*I_ERI_Dxz_Px_F2xy_Px;
  abcd[51] = I_ERI_F2yz_Px_F2xy_Px+ABZ*I_ERI_D2y_Px_F2xy_Px;
  abcd[52] = I_ERI_Fy2z_Px_F2xy_Px+ABZ*I_ERI_Dyz_Px_F2xy_Px;
  abcd[53] = I_ERI_F3z_Px_F2xy_Px+ABZ*I_ERI_D2z_Px_F2xy_Px;
  abcd[54] = I_ERI_F2xy_Py_F2xy_Px+ABY*I_ERI_D2x_Py_F2xy_Px;
  abcd[55] = I_ERI_Fx2y_Py_F2xy_Px+ABY*I_ERI_Dxy_Py_F2xy_Px;
  abcd[56] = I_ERI_Fxyz_Py_F2xy_Px+ABY*I_ERI_Dxz_Py_F2xy_Px;
  abcd[57] = I_ERI_F3y_Py_F2xy_Px+ABY*I_ERI_D2y_Py_F2xy_Px;
  abcd[58] = I_ERI_F2yz_Py_F2xy_Px+ABY*I_ERI_Dyz_Py_F2xy_Px;
  abcd[59] = I_ERI_Fy2z_Py_F2xy_Px+ABY*I_ERI_D2z_Py_F2xy_Px;
  abcd[60] = I_ERI_F2xz_Py_F2xy_Px+ABZ*I_ERI_D2x_Py_F2xy_Px;
  abcd[61] = I_ERI_Fxyz_Py_F2xy_Px+ABZ*I_ERI_Dxy_Py_F2xy_Px;
  abcd[62] = I_ERI_Fx2z_Py_F2xy_Px+ABZ*I_ERI_Dxz_Py_F2xy_Px;
  abcd[63] = I_ERI_F2yz_Py_F2xy_Px+ABZ*I_ERI_D2y_Py_F2xy_Px;
  abcd[64] = I_ERI_Fy2z_Py_F2xy_Px+ABZ*I_ERI_Dyz_Py_F2xy_Px;
  abcd[65] = I_ERI_F3z_Py_F2xy_Px+ABZ*I_ERI_D2z_Py_F2xy_Px;
  abcd[66] = I_ERI_F2xz_Pz_F2xy_Px+ABZ*I_ERI_D2x_Pz_F2xy_Px;
  abcd[67] = I_ERI_Fxyz_Pz_F2xy_Px+ABZ*I_ERI_Dxy_Pz_F2xy_Px;
  abcd[68] = I_ERI_Fx2z_Pz_F2xy_Px+ABZ*I_ERI_Dxz_Pz_F2xy_Px;
  abcd[69] = I_ERI_F2yz_Pz_F2xy_Px+ABZ*I_ERI_D2y_Pz_F2xy_Px;
  abcd[70] = I_ERI_Fy2z_Pz_F2xy_Px+ABZ*I_ERI_Dyz_Pz_F2xy_Px;
  abcd[71] = I_ERI_F3z_Pz_F2xy_Px+ABZ*I_ERI_D2z_Pz_F2xy_Px;
  abcd[72] = I_ERI_F3x_Px_F2xz_Px+ABX*I_ERI_D2x_Px_F2xz_Px;
  abcd[73] = I_ERI_F2xy_Px_F2xz_Px+ABX*I_ERI_Dxy_Px_F2xz_Px;
  abcd[74] = I_ERI_F2xz_Px_F2xz_Px+ABX*I_ERI_Dxz_Px_F2xz_Px;
  abcd[75] = I_ERI_Fx2y_Px_F2xz_Px+ABX*I_ERI_D2y_Px_F2xz_Px;
  abcd[76] = I_ERI_Fxyz_Px_F2xz_Px+ABX*I_ERI_Dyz_Px_F2xz_Px;
  abcd[77] = I_ERI_Fx2z_Px_F2xz_Px+ABX*I_ERI_D2z_Px_F2xz_Px;
  abcd[78] = I_ERI_F2xy_Px_F2xz_Px+ABY*I_ERI_D2x_Px_F2xz_Px;
  abcd[79] = I_ERI_Fx2y_Px_F2xz_Px+ABY*I_ERI_Dxy_Px_F2xz_Px;
  abcd[80] = I_ERI_Fxyz_Px_F2xz_Px+ABY*I_ERI_Dxz_Px_F2xz_Px;
  abcd[81] = I_ERI_F3y_Px_F2xz_Px+ABY*I_ERI_D2y_Px_F2xz_Px;
  abcd[82] = I_ERI_F2yz_Px_F2xz_Px+ABY*I_ERI_Dyz_Px_F2xz_Px;
  abcd[83] = I_ERI_Fy2z_Px_F2xz_Px+ABY*I_ERI_D2z_Px_F2xz_Px;
  abcd[84] = I_ERI_F2xz_Px_F2xz_Px+ABZ*I_ERI_D2x_Px_F2xz_Px;
  abcd[85] = I_ERI_Fxyz_Px_F2xz_Px+ABZ*I_ERI_Dxy_Px_F2xz_Px;
  abcd[86] = I_ERI_Fx2z_Px_F2xz_Px+ABZ*I_ERI_Dxz_Px_F2xz_Px;
  abcd[87] = I_ERI_F2yz_Px_F2xz_Px+ABZ*I_ERI_D2y_Px_F2xz_Px;
  abcd[88] = I_ERI_Fy2z_Px_F2xz_Px+ABZ*I_ERI_Dyz_Px_F2xz_Px;
  abcd[89] = I_ERI_F3z_Px_F2xz_Px+ABZ*I_ERI_D2z_Px_F2xz_Px;
  abcd[90] = I_ERI_F2xy_Py_F2xz_Px+ABY*I_ERI_D2x_Py_F2xz_Px;
  abcd[91] = I_ERI_Fx2y_Py_F2xz_Px+ABY*I_ERI_Dxy_Py_F2xz_Px;
  abcd[92] = I_ERI_Fxyz_Py_F2xz_Px+ABY*I_ERI_Dxz_Py_F2xz_Px;
  abcd[93] = I_ERI_F3y_Py_F2xz_Px+ABY*I_ERI_D2y_Py_F2xz_Px;
  abcd[94] = I_ERI_F2yz_Py_F2xz_Px+ABY*I_ERI_Dyz_Py_F2xz_Px;
  abcd[95] = I_ERI_Fy2z_Py_F2xz_Px+ABY*I_ERI_D2z_Py_F2xz_Px;
  abcd[96] = I_ERI_F2xz_Py_F2xz_Px+ABZ*I_ERI_D2x_Py_F2xz_Px;
  abcd[97] = I_ERI_Fxyz_Py_F2xz_Px+ABZ*I_ERI_Dxy_Py_F2xz_Px;
  abcd[98] = I_ERI_Fx2z_Py_F2xz_Px+ABZ*I_ERI_Dxz_Py_F2xz_Px;
  abcd[99] = I_ERI_F2yz_Py_F2xz_Px+ABZ*I_ERI_D2y_Py_F2xz_Px;
  abcd[100] = I_ERI_Fy2z_Py_F2xz_Px+ABZ*I_ERI_Dyz_Py_F2xz_Px;
  abcd[101] = I_ERI_F3z_Py_F2xz_Px+ABZ*I_ERI_D2z_Py_F2xz_Px;
  abcd[102] = I_ERI_F2xz_Pz_F2xz_Px+ABZ*I_ERI_D2x_Pz_F2xz_Px;
  abcd[103] = I_ERI_Fxyz_Pz_F2xz_Px+ABZ*I_ERI_Dxy_Pz_F2xz_Px;
  abcd[104] = I_ERI_Fx2z_Pz_F2xz_Px+ABZ*I_ERI_Dxz_Pz_F2xz_Px;
  abcd[105] = I_ERI_F2yz_Pz_F2xz_Px+ABZ*I_ERI_D2y_Pz_F2xz_Px;
  abcd[106] = I_ERI_Fy2z_Pz_F2xz_Px+ABZ*I_ERI_Dyz_Pz_F2xz_Px;
  abcd[107] = I_ERI_F3z_Pz_F2xz_Px+ABZ*I_ERI_D2z_Pz_F2xz_Px;
  abcd[108] = I_ERI_F3x_Px_Fx2y_Px+ABX*I_ERI_D2x_Px_Fx2y_Px;
  abcd[109] = I_ERI_F2xy_Px_Fx2y_Px+ABX*I_ERI_Dxy_Px_Fx2y_Px;
  abcd[110] = I_ERI_F2xz_Px_Fx2y_Px+ABX*I_ERI_Dxz_Px_Fx2y_Px;
  abcd[111] = I_ERI_Fx2y_Px_Fx2y_Px+ABX*I_ERI_D2y_Px_Fx2y_Px;
  abcd[112] = I_ERI_Fxyz_Px_Fx2y_Px+ABX*I_ERI_Dyz_Px_Fx2y_Px;
  abcd[113] = I_ERI_Fx2z_Px_Fx2y_Px+ABX*I_ERI_D2z_Px_Fx2y_Px;
  abcd[114] = I_ERI_F2xy_Px_Fx2y_Px+ABY*I_ERI_D2x_Px_Fx2y_Px;
  abcd[115] = I_ERI_Fx2y_Px_Fx2y_Px+ABY*I_ERI_Dxy_Px_Fx2y_Px;
  abcd[116] = I_ERI_Fxyz_Px_Fx2y_Px+ABY*I_ERI_Dxz_Px_Fx2y_Px;
  abcd[117] = I_ERI_F3y_Px_Fx2y_Px+ABY*I_ERI_D2y_Px_Fx2y_Px;
  abcd[118] = I_ERI_F2yz_Px_Fx2y_Px+ABY*I_ERI_Dyz_Px_Fx2y_Px;
  abcd[119] = I_ERI_Fy2z_Px_Fx2y_Px+ABY*I_ERI_D2z_Px_Fx2y_Px;
  abcd[120] = I_ERI_F2xz_Px_Fx2y_Px+ABZ*I_ERI_D2x_Px_Fx2y_Px;
  abcd[121] = I_ERI_Fxyz_Px_Fx2y_Px+ABZ*I_ERI_Dxy_Px_Fx2y_Px;
  abcd[122] = I_ERI_Fx2z_Px_Fx2y_Px+ABZ*I_ERI_Dxz_Px_Fx2y_Px;
  abcd[123] = I_ERI_F2yz_Px_Fx2y_Px+ABZ*I_ERI_D2y_Px_Fx2y_Px;
  abcd[124] = I_ERI_Fy2z_Px_Fx2y_Px+ABZ*I_ERI_Dyz_Px_Fx2y_Px;
  abcd[125] = I_ERI_F3z_Px_Fx2y_Px+ABZ*I_ERI_D2z_Px_Fx2y_Px;
  abcd[126] = I_ERI_F2xy_Py_Fx2y_Px+ABY*I_ERI_D2x_Py_Fx2y_Px;
  abcd[127] = I_ERI_Fx2y_Py_Fx2y_Px+ABY*I_ERI_Dxy_Py_Fx2y_Px;
  abcd[128] = I_ERI_Fxyz_Py_Fx2y_Px+ABY*I_ERI_Dxz_Py_Fx2y_Px;
  abcd[129] = I_ERI_F3y_Py_Fx2y_Px+ABY*I_ERI_D2y_Py_Fx2y_Px;
  abcd[130] = I_ERI_F2yz_Py_Fx2y_Px+ABY*I_ERI_Dyz_Py_Fx2y_Px;
  abcd[131] = I_ERI_Fy2z_Py_Fx2y_Px+ABY*I_ERI_D2z_Py_Fx2y_Px;
  abcd[132] = I_ERI_F2xz_Py_Fx2y_Px+ABZ*I_ERI_D2x_Py_Fx2y_Px;
  abcd[133] = I_ERI_Fxyz_Py_Fx2y_Px+ABZ*I_ERI_Dxy_Py_Fx2y_Px;
  abcd[134] = I_ERI_Fx2z_Py_Fx2y_Px+ABZ*I_ERI_Dxz_Py_Fx2y_Px;
  abcd[135] = I_ERI_F2yz_Py_Fx2y_Px+ABZ*I_ERI_D2y_Py_Fx2y_Px;
  abcd[136] = I_ERI_Fy2z_Py_Fx2y_Px+ABZ*I_ERI_Dyz_Py_Fx2y_Px;
  abcd[137] = I_ERI_F3z_Py_Fx2y_Px+ABZ*I_ERI_D2z_Py_Fx2y_Px;
  abcd[138] = I_ERI_F2xz_Pz_Fx2y_Px+ABZ*I_ERI_D2x_Pz_Fx2y_Px;
  abcd[139] = I_ERI_Fxyz_Pz_Fx2y_Px+ABZ*I_ERI_Dxy_Pz_Fx2y_Px;
  abcd[140] = I_ERI_Fx2z_Pz_Fx2y_Px+ABZ*I_ERI_Dxz_Pz_Fx2y_Px;
  abcd[141] = I_ERI_F2yz_Pz_Fx2y_Px+ABZ*I_ERI_D2y_Pz_Fx2y_Px;
  abcd[142] = I_ERI_Fy2z_Pz_Fx2y_Px+ABZ*I_ERI_Dyz_Pz_Fx2y_Px;
  abcd[143] = I_ERI_F3z_Pz_Fx2y_Px+ABZ*I_ERI_D2z_Pz_Fx2y_Px;
  abcd[144] = I_ERI_F3x_Px_Fxyz_Px+ABX*I_ERI_D2x_Px_Fxyz_Px;
  abcd[145] = I_ERI_F2xy_Px_Fxyz_Px+ABX*I_ERI_Dxy_Px_Fxyz_Px;
  abcd[146] = I_ERI_F2xz_Px_Fxyz_Px+ABX*I_ERI_Dxz_Px_Fxyz_Px;
  abcd[147] = I_ERI_Fx2y_Px_Fxyz_Px+ABX*I_ERI_D2y_Px_Fxyz_Px;
  abcd[148] = I_ERI_Fxyz_Px_Fxyz_Px+ABX*I_ERI_Dyz_Px_Fxyz_Px;
  abcd[149] = I_ERI_Fx2z_Px_Fxyz_Px+ABX*I_ERI_D2z_Px_Fxyz_Px;
  abcd[150] = I_ERI_F2xy_Px_Fxyz_Px+ABY*I_ERI_D2x_Px_Fxyz_Px;
  abcd[151] = I_ERI_Fx2y_Px_Fxyz_Px+ABY*I_ERI_Dxy_Px_Fxyz_Px;
  abcd[152] = I_ERI_Fxyz_Px_Fxyz_Px+ABY*I_ERI_Dxz_Px_Fxyz_Px;
  abcd[153] = I_ERI_F3y_Px_Fxyz_Px+ABY*I_ERI_D2y_Px_Fxyz_Px;
  abcd[154] = I_ERI_F2yz_Px_Fxyz_Px+ABY*I_ERI_Dyz_Px_Fxyz_Px;
  abcd[155] = I_ERI_Fy2z_Px_Fxyz_Px+ABY*I_ERI_D2z_Px_Fxyz_Px;
  abcd[156] = I_ERI_F2xz_Px_Fxyz_Px+ABZ*I_ERI_D2x_Px_Fxyz_Px;
  abcd[157] = I_ERI_Fxyz_Px_Fxyz_Px+ABZ*I_ERI_Dxy_Px_Fxyz_Px;
  abcd[158] = I_ERI_Fx2z_Px_Fxyz_Px+ABZ*I_ERI_Dxz_Px_Fxyz_Px;
  abcd[159] = I_ERI_F2yz_Px_Fxyz_Px+ABZ*I_ERI_D2y_Px_Fxyz_Px;
  abcd[160] = I_ERI_Fy2z_Px_Fxyz_Px+ABZ*I_ERI_Dyz_Px_Fxyz_Px;
  abcd[161] = I_ERI_F3z_Px_Fxyz_Px+ABZ*I_ERI_D2z_Px_Fxyz_Px;
  abcd[162] = I_ERI_F2xy_Py_Fxyz_Px+ABY*I_ERI_D2x_Py_Fxyz_Px;
  abcd[163] = I_ERI_Fx2y_Py_Fxyz_Px+ABY*I_ERI_Dxy_Py_Fxyz_Px;
  abcd[164] = I_ERI_Fxyz_Py_Fxyz_Px+ABY*I_ERI_Dxz_Py_Fxyz_Px;
  abcd[165] = I_ERI_F3y_Py_Fxyz_Px+ABY*I_ERI_D2y_Py_Fxyz_Px;
  abcd[166] = I_ERI_F2yz_Py_Fxyz_Px+ABY*I_ERI_Dyz_Py_Fxyz_Px;
  abcd[167] = I_ERI_Fy2z_Py_Fxyz_Px+ABY*I_ERI_D2z_Py_Fxyz_Px;
  abcd[168] = I_ERI_F2xz_Py_Fxyz_Px+ABZ*I_ERI_D2x_Py_Fxyz_Px;
  abcd[169] = I_ERI_Fxyz_Py_Fxyz_Px+ABZ*I_ERI_Dxy_Py_Fxyz_Px;
  abcd[170] = I_ERI_Fx2z_Py_Fxyz_Px+ABZ*I_ERI_Dxz_Py_Fxyz_Px;
  abcd[171] = I_ERI_F2yz_Py_Fxyz_Px+ABZ*I_ERI_D2y_Py_Fxyz_Px;
  abcd[172] = I_ERI_Fy2z_Py_Fxyz_Px+ABZ*I_ERI_Dyz_Py_Fxyz_Px;
  abcd[173] = I_ERI_F3z_Py_Fxyz_Px+ABZ*I_ERI_D2z_Py_Fxyz_Px;
  abcd[174] = I_ERI_F2xz_Pz_Fxyz_Px+ABZ*I_ERI_D2x_Pz_Fxyz_Px;
  abcd[175] = I_ERI_Fxyz_Pz_Fxyz_Px+ABZ*I_ERI_Dxy_Pz_Fxyz_Px;
  abcd[176] = I_ERI_Fx2z_Pz_Fxyz_Px+ABZ*I_ERI_Dxz_Pz_Fxyz_Px;
  abcd[177] = I_ERI_F2yz_Pz_Fxyz_Px+ABZ*I_ERI_D2y_Pz_Fxyz_Px;
  abcd[178] = I_ERI_Fy2z_Pz_Fxyz_Px+ABZ*I_ERI_Dyz_Pz_Fxyz_Px;
  abcd[179] = I_ERI_F3z_Pz_Fxyz_Px+ABZ*I_ERI_D2z_Pz_Fxyz_Px;
  abcd[180] = I_ERI_F3x_Px_Fx2z_Px+ABX*I_ERI_D2x_Px_Fx2z_Px;
  abcd[181] = I_ERI_F2xy_Px_Fx2z_Px+ABX*I_ERI_Dxy_Px_Fx2z_Px;
  abcd[182] = I_ERI_F2xz_Px_Fx2z_Px+ABX*I_ERI_Dxz_Px_Fx2z_Px;
  abcd[183] = I_ERI_Fx2y_Px_Fx2z_Px+ABX*I_ERI_D2y_Px_Fx2z_Px;
  abcd[184] = I_ERI_Fxyz_Px_Fx2z_Px+ABX*I_ERI_Dyz_Px_Fx2z_Px;
  abcd[185] = I_ERI_Fx2z_Px_Fx2z_Px+ABX*I_ERI_D2z_Px_Fx2z_Px;
  abcd[186] = I_ERI_F2xy_Px_Fx2z_Px+ABY*I_ERI_D2x_Px_Fx2z_Px;
  abcd[187] = I_ERI_Fx2y_Px_Fx2z_Px+ABY*I_ERI_Dxy_Px_Fx2z_Px;
  abcd[188] = I_ERI_Fxyz_Px_Fx2z_Px+ABY*I_ERI_Dxz_Px_Fx2z_Px;
  abcd[189] = I_ERI_F3y_Px_Fx2z_Px+ABY*I_ERI_D2y_Px_Fx2z_Px;
  abcd[190] = I_ERI_F2yz_Px_Fx2z_Px+ABY*I_ERI_Dyz_Px_Fx2z_Px;
  abcd[191] = I_ERI_Fy2z_Px_Fx2z_Px+ABY*I_ERI_D2z_Px_Fx2z_Px;
  abcd[192] = I_ERI_F2xz_Px_Fx2z_Px+ABZ*I_ERI_D2x_Px_Fx2z_Px;
  abcd[193] = I_ERI_Fxyz_Px_Fx2z_Px+ABZ*I_ERI_Dxy_Px_Fx2z_Px;
  abcd[194] = I_ERI_Fx2z_Px_Fx2z_Px+ABZ*I_ERI_Dxz_Px_Fx2z_Px;
  abcd[195] = I_ERI_F2yz_Px_Fx2z_Px+ABZ*I_ERI_D2y_Px_Fx2z_Px;
  abcd[196] = I_ERI_Fy2z_Px_Fx2z_Px+ABZ*I_ERI_Dyz_Px_Fx2z_Px;
  abcd[197] = I_ERI_F3z_Px_Fx2z_Px+ABZ*I_ERI_D2z_Px_Fx2z_Px;
  abcd[198] = I_ERI_F2xy_Py_Fx2z_Px+ABY*I_ERI_D2x_Py_Fx2z_Px;
  abcd[199] = I_ERI_Fx2y_Py_Fx2z_Px+ABY*I_ERI_Dxy_Py_Fx2z_Px;
  abcd[200] = I_ERI_Fxyz_Py_Fx2z_Px+ABY*I_ERI_Dxz_Py_Fx2z_Px;
  abcd[201] = I_ERI_F3y_Py_Fx2z_Px+ABY*I_ERI_D2y_Py_Fx2z_Px;
  abcd[202] = I_ERI_F2yz_Py_Fx2z_Px+ABY*I_ERI_Dyz_Py_Fx2z_Px;
  abcd[203] = I_ERI_Fy2z_Py_Fx2z_Px+ABY*I_ERI_D2z_Py_Fx2z_Px;
  abcd[204] = I_ERI_F2xz_Py_Fx2z_Px+ABZ*I_ERI_D2x_Py_Fx2z_Px;
  abcd[205] = I_ERI_Fxyz_Py_Fx2z_Px+ABZ*I_ERI_Dxy_Py_Fx2z_Px;
  abcd[206] = I_ERI_Fx2z_Py_Fx2z_Px+ABZ*I_ERI_Dxz_Py_Fx2z_Px;
  abcd[207] = I_ERI_F2yz_Py_Fx2z_Px+ABZ*I_ERI_D2y_Py_Fx2z_Px;
  abcd[208] = I_ERI_Fy2z_Py_Fx2z_Px+ABZ*I_ERI_Dyz_Py_Fx2z_Px;
  abcd[209] = I_ERI_F3z_Py_Fx2z_Px+ABZ*I_ERI_D2z_Py_Fx2z_Px;
  abcd[210] = I_ERI_F2xz_Pz_Fx2z_Px+ABZ*I_ERI_D2x_Pz_Fx2z_Px;
  abcd[211] = I_ERI_Fxyz_Pz_Fx2z_Px+ABZ*I_ERI_Dxy_Pz_Fx2z_Px;
  abcd[212] = I_ERI_Fx2z_Pz_Fx2z_Px+ABZ*I_ERI_Dxz_Pz_Fx2z_Px;
  abcd[213] = I_ERI_F2yz_Pz_Fx2z_Px+ABZ*I_ERI_D2y_Pz_Fx2z_Px;
  abcd[214] = I_ERI_Fy2z_Pz_Fx2z_Px+ABZ*I_ERI_Dyz_Pz_Fx2z_Px;
  abcd[215] = I_ERI_F3z_Pz_Fx2z_Px+ABZ*I_ERI_D2z_Pz_Fx2z_Px;
  abcd[216] = I_ERI_F3x_Px_F3y_Px+ABX*I_ERI_D2x_Px_F3y_Px;
  abcd[217] = I_ERI_F2xy_Px_F3y_Px+ABX*I_ERI_Dxy_Px_F3y_Px;
  abcd[218] = I_ERI_F2xz_Px_F3y_Px+ABX*I_ERI_Dxz_Px_F3y_Px;
  abcd[219] = I_ERI_Fx2y_Px_F3y_Px+ABX*I_ERI_D2y_Px_F3y_Px;
  abcd[220] = I_ERI_Fxyz_Px_F3y_Px+ABX*I_ERI_Dyz_Px_F3y_Px;
  abcd[221] = I_ERI_Fx2z_Px_F3y_Px+ABX*I_ERI_D2z_Px_F3y_Px;
  abcd[222] = I_ERI_F2xy_Px_F3y_Px+ABY*I_ERI_D2x_Px_F3y_Px;
  abcd[223] = I_ERI_Fx2y_Px_F3y_Px+ABY*I_ERI_Dxy_Px_F3y_Px;
  abcd[224] = I_ERI_Fxyz_Px_F3y_Px+ABY*I_ERI_Dxz_Px_F3y_Px;
  abcd[225] = I_ERI_F3y_Px_F3y_Px+ABY*I_ERI_D2y_Px_F3y_Px;
  abcd[226] = I_ERI_F2yz_Px_F3y_Px+ABY*I_ERI_Dyz_Px_F3y_Px;
  abcd[227] = I_ERI_Fy2z_Px_F3y_Px+ABY*I_ERI_D2z_Px_F3y_Px;
  abcd[228] = I_ERI_F2xz_Px_F3y_Px+ABZ*I_ERI_D2x_Px_F3y_Px;
  abcd[229] = I_ERI_Fxyz_Px_F3y_Px+ABZ*I_ERI_Dxy_Px_F3y_Px;
  abcd[230] = I_ERI_Fx2z_Px_F3y_Px+ABZ*I_ERI_Dxz_Px_F3y_Px;
  abcd[231] = I_ERI_F2yz_Px_F3y_Px+ABZ*I_ERI_D2y_Px_F3y_Px;
  abcd[232] = I_ERI_Fy2z_Px_F3y_Px+ABZ*I_ERI_Dyz_Px_F3y_Px;
  abcd[233] = I_ERI_F3z_Px_F3y_Px+ABZ*I_ERI_D2z_Px_F3y_Px;
  abcd[234] = I_ERI_F2xy_Py_F3y_Px+ABY*I_ERI_D2x_Py_F3y_Px;
  abcd[235] = I_ERI_Fx2y_Py_F3y_Px+ABY*I_ERI_Dxy_Py_F3y_Px;
  abcd[236] = I_ERI_Fxyz_Py_F3y_Px+ABY*I_ERI_Dxz_Py_F3y_Px;
  abcd[237] = I_ERI_F3y_Py_F3y_Px+ABY*I_ERI_D2y_Py_F3y_Px;
  abcd[238] = I_ERI_F2yz_Py_F3y_Px+ABY*I_ERI_Dyz_Py_F3y_Px;
  abcd[239] = I_ERI_Fy2z_Py_F3y_Px+ABY*I_ERI_D2z_Py_F3y_Px;
  abcd[240] = I_ERI_F2xz_Py_F3y_Px+ABZ*I_ERI_D2x_Py_F3y_Px;
  abcd[241] = I_ERI_Fxyz_Py_F3y_Px+ABZ*I_ERI_Dxy_Py_F3y_Px;
  abcd[242] = I_ERI_Fx2z_Py_F3y_Px+ABZ*I_ERI_Dxz_Py_F3y_Px;
  abcd[243] = I_ERI_F2yz_Py_F3y_Px+ABZ*I_ERI_D2y_Py_F3y_Px;
  abcd[244] = I_ERI_Fy2z_Py_F3y_Px+ABZ*I_ERI_Dyz_Py_F3y_Px;
  abcd[245] = I_ERI_F3z_Py_F3y_Px+ABZ*I_ERI_D2z_Py_F3y_Px;
  abcd[246] = I_ERI_F2xz_Pz_F3y_Px+ABZ*I_ERI_D2x_Pz_F3y_Px;
  abcd[247] = I_ERI_Fxyz_Pz_F3y_Px+ABZ*I_ERI_Dxy_Pz_F3y_Px;
  abcd[248] = I_ERI_Fx2z_Pz_F3y_Px+ABZ*I_ERI_Dxz_Pz_F3y_Px;
  abcd[249] = I_ERI_F2yz_Pz_F3y_Px+ABZ*I_ERI_D2y_Pz_F3y_Px;
  abcd[250] = I_ERI_Fy2z_Pz_F3y_Px+ABZ*I_ERI_Dyz_Pz_F3y_Px;
  abcd[251] = I_ERI_F3z_Pz_F3y_Px+ABZ*I_ERI_D2z_Pz_F3y_Px;
  abcd[252] = I_ERI_F3x_Px_F2yz_Px+ABX*I_ERI_D2x_Px_F2yz_Px;
  abcd[253] = I_ERI_F2xy_Px_F2yz_Px+ABX*I_ERI_Dxy_Px_F2yz_Px;
  abcd[254] = I_ERI_F2xz_Px_F2yz_Px+ABX*I_ERI_Dxz_Px_F2yz_Px;
  abcd[255] = I_ERI_Fx2y_Px_F2yz_Px+ABX*I_ERI_D2y_Px_F2yz_Px;
  abcd[256] = I_ERI_Fxyz_Px_F2yz_Px+ABX*I_ERI_Dyz_Px_F2yz_Px;
  abcd[257] = I_ERI_Fx2z_Px_F2yz_Px+ABX*I_ERI_D2z_Px_F2yz_Px;
  abcd[258] = I_ERI_F2xy_Px_F2yz_Px+ABY*I_ERI_D2x_Px_F2yz_Px;
  abcd[259] = I_ERI_Fx2y_Px_F2yz_Px+ABY*I_ERI_Dxy_Px_F2yz_Px;
  abcd[260] = I_ERI_Fxyz_Px_F2yz_Px+ABY*I_ERI_Dxz_Px_F2yz_Px;
  abcd[261] = I_ERI_F3y_Px_F2yz_Px+ABY*I_ERI_D2y_Px_F2yz_Px;
  abcd[262] = I_ERI_F2yz_Px_F2yz_Px+ABY*I_ERI_Dyz_Px_F2yz_Px;
  abcd[263] = I_ERI_Fy2z_Px_F2yz_Px+ABY*I_ERI_D2z_Px_F2yz_Px;
  abcd[264] = I_ERI_F2xz_Px_F2yz_Px+ABZ*I_ERI_D2x_Px_F2yz_Px;
  abcd[265] = I_ERI_Fxyz_Px_F2yz_Px+ABZ*I_ERI_Dxy_Px_F2yz_Px;
  abcd[266] = I_ERI_Fx2z_Px_F2yz_Px+ABZ*I_ERI_Dxz_Px_F2yz_Px;
  abcd[267] = I_ERI_F2yz_Px_F2yz_Px+ABZ*I_ERI_D2y_Px_F2yz_Px;
  abcd[268] = I_ERI_Fy2z_Px_F2yz_Px+ABZ*I_ERI_Dyz_Px_F2yz_Px;
  abcd[269] = I_ERI_F3z_Px_F2yz_Px+ABZ*I_ERI_D2z_Px_F2yz_Px;
  abcd[270] = I_ERI_F2xy_Py_F2yz_Px+ABY*I_ERI_D2x_Py_F2yz_Px;
  abcd[271] = I_ERI_Fx2y_Py_F2yz_Px+ABY*I_ERI_Dxy_Py_F2yz_Px;
  abcd[272] = I_ERI_Fxyz_Py_F2yz_Px+ABY*I_ERI_Dxz_Py_F2yz_Px;
  abcd[273] = I_ERI_F3y_Py_F2yz_Px+ABY*I_ERI_D2y_Py_F2yz_Px;
  abcd[274] = I_ERI_F2yz_Py_F2yz_Px+ABY*I_ERI_Dyz_Py_F2yz_Px;
  abcd[275] = I_ERI_Fy2z_Py_F2yz_Px+ABY*I_ERI_D2z_Py_F2yz_Px;
  abcd[276] = I_ERI_F2xz_Py_F2yz_Px+ABZ*I_ERI_D2x_Py_F2yz_Px;
  abcd[277] = I_ERI_Fxyz_Py_F2yz_Px+ABZ*I_ERI_Dxy_Py_F2yz_Px;
  abcd[278] = I_ERI_Fx2z_Py_F2yz_Px+ABZ*I_ERI_Dxz_Py_F2yz_Px;
  abcd[279] = I_ERI_F2yz_Py_F2yz_Px+ABZ*I_ERI_D2y_Py_F2yz_Px;
  abcd[280] = I_ERI_Fy2z_Py_F2yz_Px+ABZ*I_ERI_Dyz_Py_F2yz_Px;
  abcd[281] = I_ERI_F3z_Py_F2yz_Px+ABZ*I_ERI_D2z_Py_F2yz_Px;
  abcd[282] = I_ERI_F2xz_Pz_F2yz_Px+ABZ*I_ERI_D2x_Pz_F2yz_Px;
  abcd[283] = I_ERI_Fxyz_Pz_F2yz_Px+ABZ*I_ERI_Dxy_Pz_F2yz_Px;
  abcd[284] = I_ERI_Fx2z_Pz_F2yz_Px+ABZ*I_ERI_Dxz_Pz_F2yz_Px;
  abcd[285] = I_ERI_F2yz_Pz_F2yz_Px+ABZ*I_ERI_D2y_Pz_F2yz_Px;
  abcd[286] = I_ERI_Fy2z_Pz_F2yz_Px+ABZ*I_ERI_Dyz_Pz_F2yz_Px;
  abcd[287] = I_ERI_F3z_Pz_F2yz_Px+ABZ*I_ERI_D2z_Pz_F2yz_Px;
  abcd[288] = I_ERI_F3x_Px_Fy2z_Px+ABX*I_ERI_D2x_Px_Fy2z_Px;
  abcd[289] = I_ERI_F2xy_Px_Fy2z_Px+ABX*I_ERI_Dxy_Px_Fy2z_Px;
  abcd[290] = I_ERI_F2xz_Px_Fy2z_Px+ABX*I_ERI_Dxz_Px_Fy2z_Px;
  abcd[291] = I_ERI_Fx2y_Px_Fy2z_Px+ABX*I_ERI_D2y_Px_Fy2z_Px;
  abcd[292] = I_ERI_Fxyz_Px_Fy2z_Px+ABX*I_ERI_Dyz_Px_Fy2z_Px;
  abcd[293] = I_ERI_Fx2z_Px_Fy2z_Px+ABX*I_ERI_D2z_Px_Fy2z_Px;
  abcd[294] = I_ERI_F2xy_Px_Fy2z_Px+ABY*I_ERI_D2x_Px_Fy2z_Px;
  abcd[295] = I_ERI_Fx2y_Px_Fy2z_Px+ABY*I_ERI_Dxy_Px_Fy2z_Px;
  abcd[296] = I_ERI_Fxyz_Px_Fy2z_Px+ABY*I_ERI_Dxz_Px_Fy2z_Px;
  abcd[297] = I_ERI_F3y_Px_Fy2z_Px+ABY*I_ERI_D2y_Px_Fy2z_Px;
  abcd[298] = I_ERI_F2yz_Px_Fy2z_Px+ABY*I_ERI_Dyz_Px_Fy2z_Px;
  abcd[299] = I_ERI_Fy2z_Px_Fy2z_Px+ABY*I_ERI_D2z_Px_Fy2z_Px;
  abcd[300] = I_ERI_F2xz_Px_Fy2z_Px+ABZ*I_ERI_D2x_Px_Fy2z_Px;
  abcd[301] = I_ERI_Fxyz_Px_Fy2z_Px+ABZ*I_ERI_Dxy_Px_Fy2z_Px;
  abcd[302] = I_ERI_Fx2z_Px_Fy2z_Px+ABZ*I_ERI_Dxz_Px_Fy2z_Px;
  abcd[303] = I_ERI_F2yz_Px_Fy2z_Px+ABZ*I_ERI_D2y_Px_Fy2z_Px;
  abcd[304] = I_ERI_Fy2z_Px_Fy2z_Px+ABZ*I_ERI_Dyz_Px_Fy2z_Px;
  abcd[305] = I_ERI_F3z_Px_Fy2z_Px+ABZ*I_ERI_D2z_Px_Fy2z_Px;
  abcd[306] = I_ERI_F2xy_Py_Fy2z_Px+ABY*I_ERI_D2x_Py_Fy2z_Px;
  abcd[307] = I_ERI_Fx2y_Py_Fy2z_Px+ABY*I_ERI_Dxy_Py_Fy2z_Px;
  abcd[308] = I_ERI_Fxyz_Py_Fy2z_Px+ABY*I_ERI_Dxz_Py_Fy2z_Px;
  abcd[309] = I_ERI_F3y_Py_Fy2z_Px+ABY*I_ERI_D2y_Py_Fy2z_Px;
  abcd[310] = I_ERI_F2yz_Py_Fy2z_Px+ABY*I_ERI_Dyz_Py_Fy2z_Px;
  abcd[311] = I_ERI_Fy2z_Py_Fy2z_Px+ABY*I_ERI_D2z_Py_Fy2z_Px;
  abcd[312] = I_ERI_F2xz_Py_Fy2z_Px+ABZ*I_ERI_D2x_Py_Fy2z_Px;
  abcd[313] = I_ERI_Fxyz_Py_Fy2z_Px+ABZ*I_ERI_Dxy_Py_Fy2z_Px;
  abcd[314] = I_ERI_Fx2z_Py_Fy2z_Px+ABZ*I_ERI_Dxz_Py_Fy2z_Px;
  abcd[315] = I_ERI_F2yz_Py_Fy2z_Px+ABZ*I_ERI_D2y_Py_Fy2z_Px;
  abcd[316] = I_ERI_Fy2z_Py_Fy2z_Px+ABZ*I_ERI_Dyz_Py_Fy2z_Px;
  abcd[317] = I_ERI_F3z_Py_Fy2z_Px+ABZ*I_ERI_D2z_Py_Fy2z_Px;
  abcd[318] = I_ERI_F2xz_Pz_Fy2z_Px+ABZ*I_ERI_D2x_Pz_Fy2z_Px;
  abcd[319] = I_ERI_Fxyz_Pz_Fy2z_Px+ABZ*I_ERI_Dxy_Pz_Fy2z_Px;
  abcd[320] = I_ERI_Fx2z_Pz_Fy2z_Px+ABZ*I_ERI_Dxz_Pz_Fy2z_Px;
  abcd[321] = I_ERI_F2yz_Pz_Fy2z_Px+ABZ*I_ERI_D2y_Pz_Fy2z_Px;
  abcd[322] = I_ERI_Fy2z_Pz_Fy2z_Px+ABZ*I_ERI_Dyz_Pz_Fy2z_Px;
  abcd[323] = I_ERI_F3z_Pz_Fy2z_Px+ABZ*I_ERI_D2z_Pz_Fy2z_Px;
  abcd[324] = I_ERI_F3x_Px_F3z_Px+ABX*I_ERI_D2x_Px_F3z_Px;
  abcd[325] = I_ERI_F2xy_Px_F3z_Px+ABX*I_ERI_Dxy_Px_F3z_Px;
  abcd[326] = I_ERI_F2xz_Px_F3z_Px+ABX*I_ERI_Dxz_Px_F3z_Px;
  abcd[327] = I_ERI_Fx2y_Px_F3z_Px+ABX*I_ERI_D2y_Px_F3z_Px;
  abcd[328] = I_ERI_Fxyz_Px_F3z_Px+ABX*I_ERI_Dyz_Px_F3z_Px;
  abcd[329] = I_ERI_Fx2z_Px_F3z_Px+ABX*I_ERI_D2z_Px_F3z_Px;
  abcd[330] = I_ERI_F2xy_Px_F3z_Px+ABY*I_ERI_D2x_Px_F3z_Px;
  abcd[331] = I_ERI_Fx2y_Px_F3z_Px+ABY*I_ERI_Dxy_Px_F3z_Px;
  abcd[332] = I_ERI_Fxyz_Px_F3z_Px+ABY*I_ERI_Dxz_Px_F3z_Px;
  abcd[333] = I_ERI_F3y_Px_F3z_Px+ABY*I_ERI_D2y_Px_F3z_Px;
  abcd[334] = I_ERI_F2yz_Px_F3z_Px+ABY*I_ERI_Dyz_Px_F3z_Px;
  abcd[335] = I_ERI_Fy2z_Px_F3z_Px+ABY*I_ERI_D2z_Px_F3z_Px;
  abcd[336] = I_ERI_F2xz_Px_F3z_Px+ABZ*I_ERI_D2x_Px_F3z_Px;
  abcd[337] = I_ERI_Fxyz_Px_F3z_Px+ABZ*I_ERI_Dxy_Px_F3z_Px;
  abcd[338] = I_ERI_Fx2z_Px_F3z_Px+ABZ*I_ERI_Dxz_Px_F3z_Px;
  abcd[339] = I_ERI_F2yz_Px_F3z_Px+ABZ*I_ERI_D2y_Px_F3z_Px;
  abcd[340] = I_ERI_Fy2z_Px_F3z_Px+ABZ*I_ERI_Dyz_Px_F3z_Px;
  abcd[341] = I_ERI_F3z_Px_F3z_Px+ABZ*I_ERI_D2z_Px_F3z_Px;
  abcd[342] = I_ERI_F2xy_Py_F3z_Px+ABY*I_ERI_D2x_Py_F3z_Px;
  abcd[343] = I_ERI_Fx2y_Py_F3z_Px+ABY*I_ERI_Dxy_Py_F3z_Px;
  abcd[344] = I_ERI_Fxyz_Py_F3z_Px+ABY*I_ERI_Dxz_Py_F3z_Px;
  abcd[345] = I_ERI_F3y_Py_F3z_Px+ABY*I_ERI_D2y_Py_F3z_Px;
  abcd[346] = I_ERI_F2yz_Py_F3z_Px+ABY*I_ERI_Dyz_Py_F3z_Px;
  abcd[347] = I_ERI_Fy2z_Py_F3z_Px+ABY*I_ERI_D2z_Py_F3z_Px;
  abcd[348] = I_ERI_F2xz_Py_F3z_Px+ABZ*I_ERI_D2x_Py_F3z_Px;
  abcd[349] = I_ERI_Fxyz_Py_F3z_Px+ABZ*I_ERI_Dxy_Py_F3z_Px;
  abcd[350] = I_ERI_Fx2z_Py_F3z_Px+ABZ*I_ERI_Dxz_Py_F3z_Px;
  abcd[351] = I_ERI_F2yz_Py_F3z_Px+ABZ*I_ERI_D2y_Py_F3z_Px;
  abcd[352] = I_ERI_Fy2z_Py_F3z_Px+ABZ*I_ERI_Dyz_Py_F3z_Px;
  abcd[353] = I_ERI_F3z_Py_F3z_Px+ABZ*I_ERI_D2z_Py_F3z_Px;
  abcd[354] = I_ERI_F2xz_Pz_F3z_Px+ABZ*I_ERI_D2x_Pz_F3z_Px;
  abcd[355] = I_ERI_Fxyz_Pz_F3z_Px+ABZ*I_ERI_Dxy_Pz_F3z_Px;
  abcd[356] = I_ERI_Fx2z_Pz_F3z_Px+ABZ*I_ERI_Dxz_Pz_F3z_Px;
  abcd[357] = I_ERI_F2yz_Pz_F3z_Px+ABZ*I_ERI_D2y_Pz_F3z_Px;
  abcd[358] = I_ERI_Fy2z_Pz_F3z_Px+ABZ*I_ERI_Dyz_Pz_F3z_Px;
  abcd[359] = I_ERI_F3z_Pz_F3z_Px+ABZ*I_ERI_D2z_Pz_F3z_Px;
  abcd[360] = I_ERI_F3x_Px_F3x_Py+ABX*I_ERI_D2x_Px_F3x_Py;
  abcd[361] = I_ERI_F2xy_Px_F3x_Py+ABX*I_ERI_Dxy_Px_F3x_Py;
  abcd[362] = I_ERI_F2xz_Px_F3x_Py+ABX*I_ERI_Dxz_Px_F3x_Py;
  abcd[363] = I_ERI_Fx2y_Px_F3x_Py+ABX*I_ERI_D2y_Px_F3x_Py;
  abcd[364] = I_ERI_Fxyz_Px_F3x_Py+ABX*I_ERI_Dyz_Px_F3x_Py;
  abcd[365] = I_ERI_Fx2z_Px_F3x_Py+ABX*I_ERI_D2z_Px_F3x_Py;
  abcd[366] = I_ERI_F2xy_Px_F3x_Py+ABY*I_ERI_D2x_Px_F3x_Py;
  abcd[367] = I_ERI_Fx2y_Px_F3x_Py+ABY*I_ERI_Dxy_Px_F3x_Py;
  abcd[368] = I_ERI_Fxyz_Px_F3x_Py+ABY*I_ERI_Dxz_Px_F3x_Py;
  abcd[369] = I_ERI_F3y_Px_F3x_Py+ABY*I_ERI_D2y_Px_F3x_Py;
  abcd[370] = I_ERI_F2yz_Px_F3x_Py+ABY*I_ERI_Dyz_Px_F3x_Py;
  abcd[371] = I_ERI_Fy2z_Px_F3x_Py+ABY*I_ERI_D2z_Px_F3x_Py;
  abcd[372] = I_ERI_F2xz_Px_F3x_Py+ABZ*I_ERI_D2x_Px_F3x_Py;
  abcd[373] = I_ERI_Fxyz_Px_F3x_Py+ABZ*I_ERI_Dxy_Px_F3x_Py;
  abcd[374] = I_ERI_Fx2z_Px_F3x_Py+ABZ*I_ERI_Dxz_Px_F3x_Py;
  abcd[375] = I_ERI_F2yz_Px_F3x_Py+ABZ*I_ERI_D2y_Px_F3x_Py;
  abcd[376] = I_ERI_Fy2z_Px_F3x_Py+ABZ*I_ERI_Dyz_Px_F3x_Py;
  abcd[377] = I_ERI_F3z_Px_F3x_Py+ABZ*I_ERI_D2z_Px_F3x_Py;
  abcd[378] = I_ERI_F2xy_Py_F3x_Py+ABY*I_ERI_D2x_Py_F3x_Py;
  abcd[379] = I_ERI_Fx2y_Py_F3x_Py+ABY*I_ERI_Dxy_Py_F3x_Py;
  abcd[380] = I_ERI_Fxyz_Py_F3x_Py+ABY*I_ERI_Dxz_Py_F3x_Py;
  abcd[381] = I_ERI_F3y_Py_F3x_Py+ABY*I_ERI_D2y_Py_F3x_Py;
  abcd[382] = I_ERI_F2yz_Py_F3x_Py+ABY*I_ERI_Dyz_Py_F3x_Py;
  abcd[383] = I_ERI_Fy2z_Py_F3x_Py+ABY*I_ERI_D2z_Py_F3x_Py;
  abcd[384] = I_ERI_F2xz_Py_F3x_Py+ABZ*I_ERI_D2x_Py_F3x_Py;
  abcd[385] = I_ERI_Fxyz_Py_F3x_Py+ABZ*I_ERI_Dxy_Py_F3x_Py;
  abcd[386] = I_ERI_Fx2z_Py_F3x_Py+ABZ*I_ERI_Dxz_Py_F3x_Py;
  abcd[387] = I_ERI_F2yz_Py_F3x_Py+ABZ*I_ERI_D2y_Py_F3x_Py;
  abcd[388] = I_ERI_Fy2z_Py_F3x_Py+ABZ*I_ERI_Dyz_Py_F3x_Py;
  abcd[389] = I_ERI_F3z_Py_F3x_Py+ABZ*I_ERI_D2z_Py_F3x_Py;
  abcd[390] = I_ERI_F2xz_Pz_F3x_Py+ABZ*I_ERI_D2x_Pz_F3x_Py;
  abcd[391] = I_ERI_Fxyz_Pz_F3x_Py+ABZ*I_ERI_Dxy_Pz_F3x_Py;
  abcd[392] = I_ERI_Fx2z_Pz_F3x_Py+ABZ*I_ERI_Dxz_Pz_F3x_Py;
  abcd[393] = I_ERI_F2yz_Pz_F3x_Py+ABZ*I_ERI_D2y_Pz_F3x_Py;
  abcd[394] = I_ERI_Fy2z_Pz_F3x_Py+ABZ*I_ERI_Dyz_Pz_F3x_Py;
  abcd[395] = I_ERI_F3z_Pz_F3x_Py+ABZ*I_ERI_D2z_Pz_F3x_Py;
  abcd[396] = I_ERI_F3x_Px_F2xy_Py+ABX*I_ERI_D2x_Px_F2xy_Py;
  abcd[397] = I_ERI_F2xy_Px_F2xy_Py+ABX*I_ERI_Dxy_Px_F2xy_Py;
  abcd[398] = I_ERI_F2xz_Px_F2xy_Py+ABX*I_ERI_Dxz_Px_F2xy_Py;
  abcd[399] = I_ERI_Fx2y_Px_F2xy_Py+ABX*I_ERI_D2y_Px_F2xy_Py;
  abcd[400] = I_ERI_Fxyz_Px_F2xy_Py+ABX*I_ERI_Dyz_Px_F2xy_Py;
  abcd[401] = I_ERI_Fx2z_Px_F2xy_Py+ABX*I_ERI_D2z_Px_F2xy_Py;
  abcd[402] = I_ERI_F2xy_Px_F2xy_Py+ABY*I_ERI_D2x_Px_F2xy_Py;
  abcd[403] = I_ERI_Fx2y_Px_F2xy_Py+ABY*I_ERI_Dxy_Px_F2xy_Py;
  abcd[404] = I_ERI_Fxyz_Px_F2xy_Py+ABY*I_ERI_Dxz_Px_F2xy_Py;
  abcd[405] = I_ERI_F3y_Px_F2xy_Py+ABY*I_ERI_D2y_Px_F2xy_Py;
  abcd[406] = I_ERI_F2yz_Px_F2xy_Py+ABY*I_ERI_Dyz_Px_F2xy_Py;
  abcd[407] = I_ERI_Fy2z_Px_F2xy_Py+ABY*I_ERI_D2z_Px_F2xy_Py;
  abcd[408] = I_ERI_F2xz_Px_F2xy_Py+ABZ*I_ERI_D2x_Px_F2xy_Py;
  abcd[409] = I_ERI_Fxyz_Px_F2xy_Py+ABZ*I_ERI_Dxy_Px_F2xy_Py;
  abcd[410] = I_ERI_Fx2z_Px_F2xy_Py+ABZ*I_ERI_Dxz_Px_F2xy_Py;
  abcd[411] = I_ERI_F2yz_Px_F2xy_Py+ABZ*I_ERI_D2y_Px_F2xy_Py;
  abcd[412] = I_ERI_Fy2z_Px_F2xy_Py+ABZ*I_ERI_Dyz_Px_F2xy_Py;
  abcd[413] = I_ERI_F3z_Px_F2xy_Py+ABZ*I_ERI_D2z_Px_F2xy_Py;
  abcd[414] = I_ERI_F2xy_Py_F2xy_Py+ABY*I_ERI_D2x_Py_F2xy_Py;
  abcd[415] = I_ERI_Fx2y_Py_F2xy_Py+ABY*I_ERI_Dxy_Py_F2xy_Py;
  abcd[416] = I_ERI_Fxyz_Py_F2xy_Py+ABY*I_ERI_Dxz_Py_F2xy_Py;
  abcd[417] = I_ERI_F3y_Py_F2xy_Py+ABY*I_ERI_D2y_Py_F2xy_Py;
  abcd[418] = I_ERI_F2yz_Py_F2xy_Py+ABY*I_ERI_Dyz_Py_F2xy_Py;
  abcd[419] = I_ERI_Fy2z_Py_F2xy_Py+ABY*I_ERI_D2z_Py_F2xy_Py;
  abcd[420] = I_ERI_F2xz_Py_F2xy_Py+ABZ*I_ERI_D2x_Py_F2xy_Py;
  abcd[421] = I_ERI_Fxyz_Py_F2xy_Py+ABZ*I_ERI_Dxy_Py_F2xy_Py;
  abcd[422] = I_ERI_Fx2z_Py_F2xy_Py+ABZ*I_ERI_Dxz_Py_F2xy_Py;
  abcd[423] = I_ERI_F2yz_Py_F2xy_Py+ABZ*I_ERI_D2y_Py_F2xy_Py;
  abcd[424] = I_ERI_Fy2z_Py_F2xy_Py+ABZ*I_ERI_Dyz_Py_F2xy_Py;
  abcd[425] = I_ERI_F3z_Py_F2xy_Py+ABZ*I_ERI_D2z_Py_F2xy_Py;
  abcd[426] = I_ERI_F2xz_Pz_F2xy_Py+ABZ*I_ERI_D2x_Pz_F2xy_Py;
  abcd[427] = I_ERI_Fxyz_Pz_F2xy_Py+ABZ*I_ERI_Dxy_Pz_F2xy_Py;
  abcd[428] = I_ERI_Fx2z_Pz_F2xy_Py+ABZ*I_ERI_Dxz_Pz_F2xy_Py;
  abcd[429] = I_ERI_F2yz_Pz_F2xy_Py+ABZ*I_ERI_D2y_Pz_F2xy_Py;
  abcd[430] = I_ERI_Fy2z_Pz_F2xy_Py+ABZ*I_ERI_Dyz_Pz_F2xy_Py;
  abcd[431] = I_ERI_F3z_Pz_F2xy_Py+ABZ*I_ERI_D2z_Pz_F2xy_Py;
  abcd[432] = I_ERI_F3x_Px_F2xz_Py+ABX*I_ERI_D2x_Px_F2xz_Py;
  abcd[433] = I_ERI_F2xy_Px_F2xz_Py+ABX*I_ERI_Dxy_Px_F2xz_Py;
  abcd[434] = I_ERI_F2xz_Px_F2xz_Py+ABX*I_ERI_Dxz_Px_F2xz_Py;
  abcd[435] = I_ERI_Fx2y_Px_F2xz_Py+ABX*I_ERI_D2y_Px_F2xz_Py;
  abcd[436] = I_ERI_Fxyz_Px_F2xz_Py+ABX*I_ERI_Dyz_Px_F2xz_Py;
  abcd[437] = I_ERI_Fx2z_Px_F2xz_Py+ABX*I_ERI_D2z_Px_F2xz_Py;
  abcd[438] = I_ERI_F2xy_Px_F2xz_Py+ABY*I_ERI_D2x_Px_F2xz_Py;
  abcd[439] = I_ERI_Fx2y_Px_F2xz_Py+ABY*I_ERI_Dxy_Px_F2xz_Py;
  abcd[440] = I_ERI_Fxyz_Px_F2xz_Py+ABY*I_ERI_Dxz_Px_F2xz_Py;
  abcd[441] = I_ERI_F3y_Px_F2xz_Py+ABY*I_ERI_D2y_Px_F2xz_Py;
  abcd[442] = I_ERI_F2yz_Px_F2xz_Py+ABY*I_ERI_Dyz_Px_F2xz_Py;
  abcd[443] = I_ERI_Fy2z_Px_F2xz_Py+ABY*I_ERI_D2z_Px_F2xz_Py;
  abcd[444] = I_ERI_F2xz_Px_F2xz_Py+ABZ*I_ERI_D2x_Px_F2xz_Py;
  abcd[445] = I_ERI_Fxyz_Px_F2xz_Py+ABZ*I_ERI_Dxy_Px_F2xz_Py;
  abcd[446] = I_ERI_Fx2z_Px_F2xz_Py+ABZ*I_ERI_Dxz_Px_F2xz_Py;
  abcd[447] = I_ERI_F2yz_Px_F2xz_Py+ABZ*I_ERI_D2y_Px_F2xz_Py;
  abcd[448] = I_ERI_Fy2z_Px_F2xz_Py+ABZ*I_ERI_Dyz_Px_F2xz_Py;
  abcd[449] = I_ERI_F3z_Px_F2xz_Py+ABZ*I_ERI_D2z_Px_F2xz_Py;
  abcd[450] = I_ERI_F2xy_Py_F2xz_Py+ABY*I_ERI_D2x_Py_F2xz_Py;
  abcd[451] = I_ERI_Fx2y_Py_F2xz_Py+ABY*I_ERI_Dxy_Py_F2xz_Py;
  abcd[452] = I_ERI_Fxyz_Py_F2xz_Py+ABY*I_ERI_Dxz_Py_F2xz_Py;
  abcd[453] = I_ERI_F3y_Py_F2xz_Py+ABY*I_ERI_D2y_Py_F2xz_Py;
  abcd[454] = I_ERI_F2yz_Py_F2xz_Py+ABY*I_ERI_Dyz_Py_F2xz_Py;
  abcd[455] = I_ERI_Fy2z_Py_F2xz_Py+ABY*I_ERI_D2z_Py_F2xz_Py;
  abcd[456] = I_ERI_F2xz_Py_F2xz_Py+ABZ*I_ERI_D2x_Py_F2xz_Py;
  abcd[457] = I_ERI_Fxyz_Py_F2xz_Py+ABZ*I_ERI_Dxy_Py_F2xz_Py;
  abcd[458] = I_ERI_Fx2z_Py_F2xz_Py+ABZ*I_ERI_Dxz_Py_F2xz_Py;
  abcd[459] = I_ERI_F2yz_Py_F2xz_Py+ABZ*I_ERI_D2y_Py_F2xz_Py;
  abcd[460] = I_ERI_Fy2z_Py_F2xz_Py+ABZ*I_ERI_Dyz_Py_F2xz_Py;
  abcd[461] = I_ERI_F3z_Py_F2xz_Py+ABZ*I_ERI_D2z_Py_F2xz_Py;
  abcd[462] = I_ERI_F2xz_Pz_F2xz_Py+ABZ*I_ERI_D2x_Pz_F2xz_Py;
  abcd[463] = I_ERI_Fxyz_Pz_F2xz_Py+ABZ*I_ERI_Dxy_Pz_F2xz_Py;
  abcd[464] = I_ERI_Fx2z_Pz_F2xz_Py+ABZ*I_ERI_Dxz_Pz_F2xz_Py;
  abcd[465] = I_ERI_F2yz_Pz_F2xz_Py+ABZ*I_ERI_D2y_Pz_F2xz_Py;
  abcd[466] = I_ERI_Fy2z_Pz_F2xz_Py+ABZ*I_ERI_Dyz_Pz_F2xz_Py;
  abcd[467] = I_ERI_F3z_Pz_F2xz_Py+ABZ*I_ERI_D2z_Pz_F2xz_Py;
  abcd[468] = I_ERI_F3x_Px_Fx2y_Py+ABX*I_ERI_D2x_Px_Fx2y_Py;
  abcd[469] = I_ERI_F2xy_Px_Fx2y_Py+ABX*I_ERI_Dxy_Px_Fx2y_Py;
  abcd[470] = I_ERI_F2xz_Px_Fx2y_Py+ABX*I_ERI_Dxz_Px_Fx2y_Py;
  abcd[471] = I_ERI_Fx2y_Px_Fx2y_Py+ABX*I_ERI_D2y_Px_Fx2y_Py;
  abcd[472] = I_ERI_Fxyz_Px_Fx2y_Py+ABX*I_ERI_Dyz_Px_Fx2y_Py;
  abcd[473] = I_ERI_Fx2z_Px_Fx2y_Py+ABX*I_ERI_D2z_Px_Fx2y_Py;
  abcd[474] = I_ERI_F2xy_Px_Fx2y_Py+ABY*I_ERI_D2x_Px_Fx2y_Py;
  abcd[475] = I_ERI_Fx2y_Px_Fx2y_Py+ABY*I_ERI_Dxy_Px_Fx2y_Py;
  abcd[476] = I_ERI_Fxyz_Px_Fx2y_Py+ABY*I_ERI_Dxz_Px_Fx2y_Py;
  abcd[477] = I_ERI_F3y_Px_Fx2y_Py+ABY*I_ERI_D2y_Px_Fx2y_Py;
  abcd[478] = I_ERI_F2yz_Px_Fx2y_Py+ABY*I_ERI_Dyz_Px_Fx2y_Py;
  abcd[479] = I_ERI_Fy2z_Px_Fx2y_Py+ABY*I_ERI_D2z_Px_Fx2y_Py;
  abcd[480] = I_ERI_F2xz_Px_Fx2y_Py+ABZ*I_ERI_D2x_Px_Fx2y_Py;
  abcd[481] = I_ERI_Fxyz_Px_Fx2y_Py+ABZ*I_ERI_Dxy_Px_Fx2y_Py;
  abcd[482] = I_ERI_Fx2z_Px_Fx2y_Py+ABZ*I_ERI_Dxz_Px_Fx2y_Py;
  abcd[483] = I_ERI_F2yz_Px_Fx2y_Py+ABZ*I_ERI_D2y_Px_Fx2y_Py;
  abcd[484] = I_ERI_Fy2z_Px_Fx2y_Py+ABZ*I_ERI_Dyz_Px_Fx2y_Py;
  abcd[485] = I_ERI_F3z_Px_Fx2y_Py+ABZ*I_ERI_D2z_Px_Fx2y_Py;
  abcd[486] = I_ERI_F2xy_Py_Fx2y_Py+ABY*I_ERI_D2x_Py_Fx2y_Py;
  abcd[487] = I_ERI_Fx2y_Py_Fx2y_Py+ABY*I_ERI_Dxy_Py_Fx2y_Py;
  abcd[488] = I_ERI_Fxyz_Py_Fx2y_Py+ABY*I_ERI_Dxz_Py_Fx2y_Py;
  abcd[489] = I_ERI_F3y_Py_Fx2y_Py+ABY*I_ERI_D2y_Py_Fx2y_Py;
  abcd[490] = I_ERI_F2yz_Py_Fx2y_Py+ABY*I_ERI_Dyz_Py_Fx2y_Py;
  abcd[491] = I_ERI_Fy2z_Py_Fx2y_Py+ABY*I_ERI_D2z_Py_Fx2y_Py;
  abcd[492] = I_ERI_F2xz_Py_Fx2y_Py+ABZ*I_ERI_D2x_Py_Fx2y_Py;
  abcd[493] = I_ERI_Fxyz_Py_Fx2y_Py+ABZ*I_ERI_Dxy_Py_Fx2y_Py;
  abcd[494] = I_ERI_Fx2z_Py_Fx2y_Py+ABZ*I_ERI_Dxz_Py_Fx2y_Py;
  abcd[495] = I_ERI_F2yz_Py_Fx2y_Py+ABZ*I_ERI_D2y_Py_Fx2y_Py;
  abcd[496] = I_ERI_Fy2z_Py_Fx2y_Py+ABZ*I_ERI_Dyz_Py_Fx2y_Py;
  abcd[497] = I_ERI_F3z_Py_Fx2y_Py+ABZ*I_ERI_D2z_Py_Fx2y_Py;
  abcd[498] = I_ERI_F2xz_Pz_Fx2y_Py+ABZ*I_ERI_D2x_Pz_Fx2y_Py;
  abcd[499] = I_ERI_Fxyz_Pz_Fx2y_Py+ABZ*I_ERI_Dxy_Pz_Fx2y_Py;
  abcd[500] = I_ERI_Fx2z_Pz_Fx2y_Py+ABZ*I_ERI_Dxz_Pz_Fx2y_Py;
  abcd[501] = I_ERI_F2yz_Pz_Fx2y_Py+ABZ*I_ERI_D2y_Pz_Fx2y_Py;
  abcd[502] = I_ERI_Fy2z_Pz_Fx2y_Py+ABZ*I_ERI_Dyz_Pz_Fx2y_Py;
  abcd[503] = I_ERI_F3z_Pz_Fx2y_Py+ABZ*I_ERI_D2z_Pz_Fx2y_Py;
  abcd[504] = I_ERI_F3x_Px_Fxyz_Py+ABX*I_ERI_D2x_Px_Fxyz_Py;
  abcd[505] = I_ERI_F2xy_Px_Fxyz_Py+ABX*I_ERI_Dxy_Px_Fxyz_Py;
  abcd[506] = I_ERI_F2xz_Px_Fxyz_Py+ABX*I_ERI_Dxz_Px_Fxyz_Py;
  abcd[507] = I_ERI_Fx2y_Px_Fxyz_Py+ABX*I_ERI_D2y_Px_Fxyz_Py;
  abcd[508] = I_ERI_Fxyz_Px_Fxyz_Py+ABX*I_ERI_Dyz_Px_Fxyz_Py;
  abcd[509] = I_ERI_Fx2z_Px_Fxyz_Py+ABX*I_ERI_D2z_Px_Fxyz_Py;
  abcd[510] = I_ERI_F2xy_Px_Fxyz_Py+ABY*I_ERI_D2x_Px_Fxyz_Py;
  abcd[511] = I_ERI_Fx2y_Px_Fxyz_Py+ABY*I_ERI_Dxy_Px_Fxyz_Py;
  abcd[512] = I_ERI_Fxyz_Px_Fxyz_Py+ABY*I_ERI_Dxz_Px_Fxyz_Py;
  abcd[513] = I_ERI_F3y_Px_Fxyz_Py+ABY*I_ERI_D2y_Px_Fxyz_Py;
  abcd[514] = I_ERI_F2yz_Px_Fxyz_Py+ABY*I_ERI_Dyz_Px_Fxyz_Py;
  abcd[515] = I_ERI_Fy2z_Px_Fxyz_Py+ABY*I_ERI_D2z_Px_Fxyz_Py;
  abcd[516] = I_ERI_F2xz_Px_Fxyz_Py+ABZ*I_ERI_D2x_Px_Fxyz_Py;
  abcd[517] = I_ERI_Fxyz_Px_Fxyz_Py+ABZ*I_ERI_Dxy_Px_Fxyz_Py;
  abcd[518] = I_ERI_Fx2z_Px_Fxyz_Py+ABZ*I_ERI_Dxz_Px_Fxyz_Py;
  abcd[519] = I_ERI_F2yz_Px_Fxyz_Py+ABZ*I_ERI_D2y_Px_Fxyz_Py;
  abcd[520] = I_ERI_Fy2z_Px_Fxyz_Py+ABZ*I_ERI_Dyz_Px_Fxyz_Py;
  abcd[521] = I_ERI_F3z_Px_Fxyz_Py+ABZ*I_ERI_D2z_Px_Fxyz_Py;
  abcd[522] = I_ERI_F2xy_Py_Fxyz_Py+ABY*I_ERI_D2x_Py_Fxyz_Py;
  abcd[523] = I_ERI_Fx2y_Py_Fxyz_Py+ABY*I_ERI_Dxy_Py_Fxyz_Py;
  abcd[524] = I_ERI_Fxyz_Py_Fxyz_Py+ABY*I_ERI_Dxz_Py_Fxyz_Py;
  abcd[525] = I_ERI_F3y_Py_Fxyz_Py+ABY*I_ERI_D2y_Py_Fxyz_Py;
  abcd[526] = I_ERI_F2yz_Py_Fxyz_Py+ABY*I_ERI_Dyz_Py_Fxyz_Py;
  abcd[527] = I_ERI_Fy2z_Py_Fxyz_Py+ABY*I_ERI_D2z_Py_Fxyz_Py;
  abcd[528] = I_ERI_F2xz_Py_Fxyz_Py+ABZ*I_ERI_D2x_Py_Fxyz_Py;
  abcd[529] = I_ERI_Fxyz_Py_Fxyz_Py+ABZ*I_ERI_Dxy_Py_Fxyz_Py;
  abcd[530] = I_ERI_Fx2z_Py_Fxyz_Py+ABZ*I_ERI_Dxz_Py_Fxyz_Py;
  abcd[531] = I_ERI_F2yz_Py_Fxyz_Py+ABZ*I_ERI_D2y_Py_Fxyz_Py;
  abcd[532] = I_ERI_Fy2z_Py_Fxyz_Py+ABZ*I_ERI_Dyz_Py_Fxyz_Py;
  abcd[533] = I_ERI_F3z_Py_Fxyz_Py+ABZ*I_ERI_D2z_Py_Fxyz_Py;
  abcd[534] = I_ERI_F2xz_Pz_Fxyz_Py+ABZ*I_ERI_D2x_Pz_Fxyz_Py;
  abcd[535] = I_ERI_Fxyz_Pz_Fxyz_Py+ABZ*I_ERI_Dxy_Pz_Fxyz_Py;
  abcd[536] = I_ERI_Fx2z_Pz_Fxyz_Py+ABZ*I_ERI_Dxz_Pz_Fxyz_Py;
  abcd[537] = I_ERI_F2yz_Pz_Fxyz_Py+ABZ*I_ERI_D2y_Pz_Fxyz_Py;
  abcd[538] = I_ERI_Fy2z_Pz_Fxyz_Py+ABZ*I_ERI_Dyz_Pz_Fxyz_Py;
  abcd[539] = I_ERI_F3z_Pz_Fxyz_Py+ABZ*I_ERI_D2z_Pz_Fxyz_Py;
  abcd[540] = I_ERI_F3x_Px_Fx2z_Py+ABX*I_ERI_D2x_Px_Fx2z_Py;
  abcd[541] = I_ERI_F2xy_Px_Fx2z_Py+ABX*I_ERI_Dxy_Px_Fx2z_Py;
  abcd[542] = I_ERI_F2xz_Px_Fx2z_Py+ABX*I_ERI_Dxz_Px_Fx2z_Py;
  abcd[543] = I_ERI_Fx2y_Px_Fx2z_Py+ABX*I_ERI_D2y_Px_Fx2z_Py;
  abcd[544] = I_ERI_Fxyz_Px_Fx2z_Py+ABX*I_ERI_Dyz_Px_Fx2z_Py;
  abcd[545] = I_ERI_Fx2z_Px_Fx2z_Py+ABX*I_ERI_D2z_Px_Fx2z_Py;
  abcd[546] = I_ERI_F2xy_Px_Fx2z_Py+ABY*I_ERI_D2x_Px_Fx2z_Py;
  abcd[547] = I_ERI_Fx2y_Px_Fx2z_Py+ABY*I_ERI_Dxy_Px_Fx2z_Py;
  abcd[548] = I_ERI_Fxyz_Px_Fx2z_Py+ABY*I_ERI_Dxz_Px_Fx2z_Py;
  abcd[549] = I_ERI_F3y_Px_Fx2z_Py+ABY*I_ERI_D2y_Px_Fx2z_Py;
  abcd[550] = I_ERI_F2yz_Px_Fx2z_Py+ABY*I_ERI_Dyz_Px_Fx2z_Py;
  abcd[551] = I_ERI_Fy2z_Px_Fx2z_Py+ABY*I_ERI_D2z_Px_Fx2z_Py;
  abcd[552] = I_ERI_F2xz_Px_Fx2z_Py+ABZ*I_ERI_D2x_Px_Fx2z_Py;
  abcd[553] = I_ERI_Fxyz_Px_Fx2z_Py+ABZ*I_ERI_Dxy_Px_Fx2z_Py;
  abcd[554] = I_ERI_Fx2z_Px_Fx2z_Py+ABZ*I_ERI_Dxz_Px_Fx2z_Py;
  abcd[555] = I_ERI_F2yz_Px_Fx2z_Py+ABZ*I_ERI_D2y_Px_Fx2z_Py;
  abcd[556] = I_ERI_Fy2z_Px_Fx2z_Py+ABZ*I_ERI_Dyz_Px_Fx2z_Py;
  abcd[557] = I_ERI_F3z_Px_Fx2z_Py+ABZ*I_ERI_D2z_Px_Fx2z_Py;
  abcd[558] = I_ERI_F2xy_Py_Fx2z_Py+ABY*I_ERI_D2x_Py_Fx2z_Py;
  abcd[559] = I_ERI_Fx2y_Py_Fx2z_Py+ABY*I_ERI_Dxy_Py_Fx2z_Py;
  abcd[560] = I_ERI_Fxyz_Py_Fx2z_Py+ABY*I_ERI_Dxz_Py_Fx2z_Py;
  abcd[561] = I_ERI_F3y_Py_Fx2z_Py+ABY*I_ERI_D2y_Py_Fx2z_Py;
  abcd[562] = I_ERI_F2yz_Py_Fx2z_Py+ABY*I_ERI_Dyz_Py_Fx2z_Py;
  abcd[563] = I_ERI_Fy2z_Py_Fx2z_Py+ABY*I_ERI_D2z_Py_Fx2z_Py;
  abcd[564] = I_ERI_F2xz_Py_Fx2z_Py+ABZ*I_ERI_D2x_Py_Fx2z_Py;
  abcd[565] = I_ERI_Fxyz_Py_Fx2z_Py+ABZ*I_ERI_Dxy_Py_Fx2z_Py;
  abcd[566] = I_ERI_Fx2z_Py_Fx2z_Py+ABZ*I_ERI_Dxz_Py_Fx2z_Py;
  abcd[567] = I_ERI_F2yz_Py_Fx2z_Py+ABZ*I_ERI_D2y_Py_Fx2z_Py;
  abcd[568] = I_ERI_Fy2z_Py_Fx2z_Py+ABZ*I_ERI_Dyz_Py_Fx2z_Py;
  abcd[569] = I_ERI_F3z_Py_Fx2z_Py+ABZ*I_ERI_D2z_Py_Fx2z_Py;
  abcd[570] = I_ERI_F2xz_Pz_Fx2z_Py+ABZ*I_ERI_D2x_Pz_Fx2z_Py;
  abcd[571] = I_ERI_Fxyz_Pz_Fx2z_Py+ABZ*I_ERI_Dxy_Pz_Fx2z_Py;
  abcd[572] = I_ERI_Fx2z_Pz_Fx2z_Py+ABZ*I_ERI_Dxz_Pz_Fx2z_Py;
  abcd[573] = I_ERI_F2yz_Pz_Fx2z_Py+ABZ*I_ERI_D2y_Pz_Fx2z_Py;
  abcd[574] = I_ERI_Fy2z_Pz_Fx2z_Py+ABZ*I_ERI_Dyz_Pz_Fx2z_Py;
  abcd[575] = I_ERI_F3z_Pz_Fx2z_Py+ABZ*I_ERI_D2z_Pz_Fx2z_Py;
  abcd[576] = I_ERI_F3x_Px_F3y_Py+ABX*I_ERI_D2x_Px_F3y_Py;
  abcd[577] = I_ERI_F2xy_Px_F3y_Py+ABX*I_ERI_Dxy_Px_F3y_Py;
  abcd[578] = I_ERI_F2xz_Px_F3y_Py+ABX*I_ERI_Dxz_Px_F3y_Py;
  abcd[579] = I_ERI_Fx2y_Px_F3y_Py+ABX*I_ERI_D2y_Px_F3y_Py;
  abcd[580] = I_ERI_Fxyz_Px_F3y_Py+ABX*I_ERI_Dyz_Px_F3y_Py;
  abcd[581] = I_ERI_Fx2z_Px_F3y_Py+ABX*I_ERI_D2z_Px_F3y_Py;
  abcd[582] = I_ERI_F2xy_Px_F3y_Py+ABY*I_ERI_D2x_Px_F3y_Py;
  abcd[583] = I_ERI_Fx2y_Px_F3y_Py+ABY*I_ERI_Dxy_Px_F3y_Py;
  abcd[584] = I_ERI_Fxyz_Px_F3y_Py+ABY*I_ERI_Dxz_Px_F3y_Py;
  abcd[585] = I_ERI_F3y_Px_F3y_Py+ABY*I_ERI_D2y_Px_F3y_Py;
  abcd[586] = I_ERI_F2yz_Px_F3y_Py+ABY*I_ERI_Dyz_Px_F3y_Py;
  abcd[587] = I_ERI_Fy2z_Px_F3y_Py+ABY*I_ERI_D2z_Px_F3y_Py;
  abcd[588] = I_ERI_F2xz_Px_F3y_Py+ABZ*I_ERI_D2x_Px_F3y_Py;
  abcd[589] = I_ERI_Fxyz_Px_F3y_Py+ABZ*I_ERI_Dxy_Px_F3y_Py;
  abcd[590] = I_ERI_Fx2z_Px_F3y_Py+ABZ*I_ERI_Dxz_Px_F3y_Py;
  abcd[591] = I_ERI_F2yz_Px_F3y_Py+ABZ*I_ERI_D2y_Px_F3y_Py;
  abcd[592] = I_ERI_Fy2z_Px_F3y_Py+ABZ*I_ERI_Dyz_Px_F3y_Py;
  abcd[593] = I_ERI_F3z_Px_F3y_Py+ABZ*I_ERI_D2z_Px_F3y_Py;
  abcd[594] = I_ERI_F2xy_Py_F3y_Py+ABY*I_ERI_D2x_Py_F3y_Py;
  abcd[595] = I_ERI_Fx2y_Py_F3y_Py+ABY*I_ERI_Dxy_Py_F3y_Py;
  abcd[596] = I_ERI_Fxyz_Py_F3y_Py+ABY*I_ERI_Dxz_Py_F3y_Py;
  abcd[597] = I_ERI_F3y_Py_F3y_Py+ABY*I_ERI_D2y_Py_F3y_Py;
  abcd[598] = I_ERI_F2yz_Py_F3y_Py+ABY*I_ERI_Dyz_Py_F3y_Py;
  abcd[599] = I_ERI_Fy2z_Py_F3y_Py+ABY*I_ERI_D2z_Py_F3y_Py;
  abcd[600] = I_ERI_F2xz_Py_F3y_Py+ABZ*I_ERI_D2x_Py_F3y_Py;
  abcd[601] = I_ERI_Fxyz_Py_F3y_Py+ABZ*I_ERI_Dxy_Py_F3y_Py;
  abcd[602] = I_ERI_Fx2z_Py_F3y_Py+ABZ*I_ERI_Dxz_Py_F3y_Py;
  abcd[603] = I_ERI_F2yz_Py_F3y_Py+ABZ*I_ERI_D2y_Py_F3y_Py;
  abcd[604] = I_ERI_Fy2z_Py_F3y_Py+ABZ*I_ERI_Dyz_Py_F3y_Py;
  abcd[605] = I_ERI_F3z_Py_F3y_Py+ABZ*I_ERI_D2z_Py_F3y_Py;
  abcd[606] = I_ERI_F2xz_Pz_F3y_Py+ABZ*I_ERI_D2x_Pz_F3y_Py;
  abcd[607] = I_ERI_Fxyz_Pz_F3y_Py+ABZ*I_ERI_Dxy_Pz_F3y_Py;
  abcd[608] = I_ERI_Fx2z_Pz_F3y_Py+ABZ*I_ERI_Dxz_Pz_F3y_Py;
  abcd[609] = I_ERI_F2yz_Pz_F3y_Py+ABZ*I_ERI_D2y_Pz_F3y_Py;
  abcd[610] = I_ERI_Fy2z_Pz_F3y_Py+ABZ*I_ERI_Dyz_Pz_F3y_Py;
  abcd[611] = I_ERI_F3z_Pz_F3y_Py+ABZ*I_ERI_D2z_Pz_F3y_Py;
  abcd[612] = I_ERI_F3x_Px_F2yz_Py+ABX*I_ERI_D2x_Px_F2yz_Py;
  abcd[613] = I_ERI_F2xy_Px_F2yz_Py+ABX*I_ERI_Dxy_Px_F2yz_Py;
  abcd[614] = I_ERI_F2xz_Px_F2yz_Py+ABX*I_ERI_Dxz_Px_F2yz_Py;
  abcd[615] = I_ERI_Fx2y_Px_F2yz_Py+ABX*I_ERI_D2y_Px_F2yz_Py;
  abcd[616] = I_ERI_Fxyz_Px_F2yz_Py+ABX*I_ERI_Dyz_Px_F2yz_Py;
  abcd[617] = I_ERI_Fx2z_Px_F2yz_Py+ABX*I_ERI_D2z_Px_F2yz_Py;
  abcd[618] = I_ERI_F2xy_Px_F2yz_Py+ABY*I_ERI_D2x_Px_F2yz_Py;
  abcd[619] = I_ERI_Fx2y_Px_F2yz_Py+ABY*I_ERI_Dxy_Px_F2yz_Py;
  abcd[620] = I_ERI_Fxyz_Px_F2yz_Py+ABY*I_ERI_Dxz_Px_F2yz_Py;
  abcd[621] = I_ERI_F3y_Px_F2yz_Py+ABY*I_ERI_D2y_Px_F2yz_Py;
  abcd[622] = I_ERI_F2yz_Px_F2yz_Py+ABY*I_ERI_Dyz_Px_F2yz_Py;
  abcd[623] = I_ERI_Fy2z_Px_F2yz_Py+ABY*I_ERI_D2z_Px_F2yz_Py;
  abcd[624] = I_ERI_F2xz_Px_F2yz_Py+ABZ*I_ERI_D2x_Px_F2yz_Py;
  abcd[625] = I_ERI_Fxyz_Px_F2yz_Py+ABZ*I_ERI_Dxy_Px_F2yz_Py;
  abcd[626] = I_ERI_Fx2z_Px_F2yz_Py+ABZ*I_ERI_Dxz_Px_F2yz_Py;
  abcd[627] = I_ERI_F2yz_Px_F2yz_Py+ABZ*I_ERI_D2y_Px_F2yz_Py;
  abcd[628] = I_ERI_Fy2z_Px_F2yz_Py+ABZ*I_ERI_Dyz_Px_F2yz_Py;
  abcd[629] = I_ERI_F3z_Px_F2yz_Py+ABZ*I_ERI_D2z_Px_F2yz_Py;
  abcd[630] = I_ERI_F2xy_Py_F2yz_Py+ABY*I_ERI_D2x_Py_F2yz_Py;
  abcd[631] = I_ERI_Fx2y_Py_F2yz_Py+ABY*I_ERI_Dxy_Py_F2yz_Py;
  abcd[632] = I_ERI_Fxyz_Py_F2yz_Py+ABY*I_ERI_Dxz_Py_F2yz_Py;
  abcd[633] = I_ERI_F3y_Py_F2yz_Py+ABY*I_ERI_D2y_Py_F2yz_Py;
  abcd[634] = I_ERI_F2yz_Py_F2yz_Py+ABY*I_ERI_Dyz_Py_F2yz_Py;
  abcd[635] = I_ERI_Fy2z_Py_F2yz_Py+ABY*I_ERI_D2z_Py_F2yz_Py;
  abcd[636] = I_ERI_F2xz_Py_F2yz_Py+ABZ*I_ERI_D2x_Py_F2yz_Py;
  abcd[637] = I_ERI_Fxyz_Py_F2yz_Py+ABZ*I_ERI_Dxy_Py_F2yz_Py;
  abcd[638] = I_ERI_Fx2z_Py_F2yz_Py+ABZ*I_ERI_Dxz_Py_F2yz_Py;
  abcd[639] = I_ERI_F2yz_Py_F2yz_Py+ABZ*I_ERI_D2y_Py_F2yz_Py;
  abcd[640] = I_ERI_Fy2z_Py_F2yz_Py+ABZ*I_ERI_Dyz_Py_F2yz_Py;
  abcd[641] = I_ERI_F3z_Py_F2yz_Py+ABZ*I_ERI_D2z_Py_F2yz_Py;
  abcd[642] = I_ERI_F2xz_Pz_F2yz_Py+ABZ*I_ERI_D2x_Pz_F2yz_Py;
  abcd[643] = I_ERI_Fxyz_Pz_F2yz_Py+ABZ*I_ERI_Dxy_Pz_F2yz_Py;
  abcd[644] = I_ERI_Fx2z_Pz_F2yz_Py+ABZ*I_ERI_Dxz_Pz_F2yz_Py;
  abcd[645] = I_ERI_F2yz_Pz_F2yz_Py+ABZ*I_ERI_D2y_Pz_F2yz_Py;
  abcd[646] = I_ERI_Fy2z_Pz_F2yz_Py+ABZ*I_ERI_Dyz_Pz_F2yz_Py;
  abcd[647] = I_ERI_F3z_Pz_F2yz_Py+ABZ*I_ERI_D2z_Pz_F2yz_Py;
  abcd[648] = I_ERI_F3x_Px_Fy2z_Py+ABX*I_ERI_D2x_Px_Fy2z_Py;
  abcd[649] = I_ERI_F2xy_Px_Fy2z_Py+ABX*I_ERI_Dxy_Px_Fy2z_Py;
  abcd[650] = I_ERI_F2xz_Px_Fy2z_Py+ABX*I_ERI_Dxz_Px_Fy2z_Py;
  abcd[651] = I_ERI_Fx2y_Px_Fy2z_Py+ABX*I_ERI_D2y_Px_Fy2z_Py;
  abcd[652] = I_ERI_Fxyz_Px_Fy2z_Py+ABX*I_ERI_Dyz_Px_Fy2z_Py;
  abcd[653] = I_ERI_Fx2z_Px_Fy2z_Py+ABX*I_ERI_D2z_Px_Fy2z_Py;
  abcd[654] = I_ERI_F2xy_Px_Fy2z_Py+ABY*I_ERI_D2x_Px_Fy2z_Py;
  abcd[655] = I_ERI_Fx2y_Px_Fy2z_Py+ABY*I_ERI_Dxy_Px_Fy2z_Py;
  abcd[656] = I_ERI_Fxyz_Px_Fy2z_Py+ABY*I_ERI_Dxz_Px_Fy2z_Py;
  abcd[657] = I_ERI_F3y_Px_Fy2z_Py+ABY*I_ERI_D2y_Px_Fy2z_Py;
  abcd[658] = I_ERI_F2yz_Px_Fy2z_Py+ABY*I_ERI_Dyz_Px_Fy2z_Py;
  abcd[659] = I_ERI_Fy2z_Px_Fy2z_Py+ABY*I_ERI_D2z_Px_Fy2z_Py;
  abcd[660] = I_ERI_F2xz_Px_Fy2z_Py+ABZ*I_ERI_D2x_Px_Fy2z_Py;
  abcd[661] = I_ERI_Fxyz_Px_Fy2z_Py+ABZ*I_ERI_Dxy_Px_Fy2z_Py;
  abcd[662] = I_ERI_Fx2z_Px_Fy2z_Py+ABZ*I_ERI_Dxz_Px_Fy2z_Py;
  abcd[663] = I_ERI_F2yz_Px_Fy2z_Py+ABZ*I_ERI_D2y_Px_Fy2z_Py;
  abcd[664] = I_ERI_Fy2z_Px_Fy2z_Py+ABZ*I_ERI_Dyz_Px_Fy2z_Py;
  abcd[665] = I_ERI_F3z_Px_Fy2z_Py+ABZ*I_ERI_D2z_Px_Fy2z_Py;
  abcd[666] = I_ERI_F2xy_Py_Fy2z_Py+ABY*I_ERI_D2x_Py_Fy2z_Py;
  abcd[667] = I_ERI_Fx2y_Py_Fy2z_Py+ABY*I_ERI_Dxy_Py_Fy2z_Py;
  abcd[668] = I_ERI_Fxyz_Py_Fy2z_Py+ABY*I_ERI_Dxz_Py_Fy2z_Py;
  abcd[669] = I_ERI_F3y_Py_Fy2z_Py+ABY*I_ERI_D2y_Py_Fy2z_Py;
  abcd[670] = I_ERI_F2yz_Py_Fy2z_Py+ABY*I_ERI_Dyz_Py_Fy2z_Py;
  abcd[671] = I_ERI_Fy2z_Py_Fy2z_Py+ABY*I_ERI_D2z_Py_Fy2z_Py;
  abcd[672] = I_ERI_F2xz_Py_Fy2z_Py+ABZ*I_ERI_D2x_Py_Fy2z_Py;
  abcd[673] = I_ERI_Fxyz_Py_Fy2z_Py+ABZ*I_ERI_Dxy_Py_Fy2z_Py;
  abcd[674] = I_ERI_Fx2z_Py_Fy2z_Py+ABZ*I_ERI_Dxz_Py_Fy2z_Py;
  abcd[675] = I_ERI_F2yz_Py_Fy2z_Py+ABZ*I_ERI_D2y_Py_Fy2z_Py;
  abcd[676] = I_ERI_Fy2z_Py_Fy2z_Py+ABZ*I_ERI_Dyz_Py_Fy2z_Py;
  abcd[677] = I_ERI_F3z_Py_Fy2z_Py+ABZ*I_ERI_D2z_Py_Fy2z_Py;
  abcd[678] = I_ERI_F2xz_Pz_Fy2z_Py+ABZ*I_ERI_D2x_Pz_Fy2z_Py;
  abcd[679] = I_ERI_Fxyz_Pz_Fy2z_Py+ABZ*I_ERI_Dxy_Pz_Fy2z_Py;
  abcd[680] = I_ERI_Fx2z_Pz_Fy2z_Py+ABZ*I_ERI_Dxz_Pz_Fy2z_Py;
  abcd[681] = I_ERI_F2yz_Pz_Fy2z_Py+ABZ*I_ERI_D2y_Pz_Fy2z_Py;
  abcd[682] = I_ERI_Fy2z_Pz_Fy2z_Py+ABZ*I_ERI_Dyz_Pz_Fy2z_Py;
  abcd[683] = I_ERI_F3z_Pz_Fy2z_Py+ABZ*I_ERI_D2z_Pz_Fy2z_Py;
  abcd[684] = I_ERI_F3x_Px_F3z_Py+ABX*I_ERI_D2x_Px_F3z_Py;
  abcd[685] = I_ERI_F2xy_Px_F3z_Py+ABX*I_ERI_Dxy_Px_F3z_Py;
  abcd[686] = I_ERI_F2xz_Px_F3z_Py+ABX*I_ERI_Dxz_Px_F3z_Py;
  abcd[687] = I_ERI_Fx2y_Px_F3z_Py+ABX*I_ERI_D2y_Px_F3z_Py;
  abcd[688] = I_ERI_Fxyz_Px_F3z_Py+ABX*I_ERI_Dyz_Px_F3z_Py;
  abcd[689] = I_ERI_Fx2z_Px_F3z_Py+ABX*I_ERI_D2z_Px_F3z_Py;
  abcd[690] = I_ERI_F2xy_Px_F3z_Py+ABY*I_ERI_D2x_Px_F3z_Py;
  abcd[691] = I_ERI_Fx2y_Px_F3z_Py+ABY*I_ERI_Dxy_Px_F3z_Py;
  abcd[692] = I_ERI_Fxyz_Px_F3z_Py+ABY*I_ERI_Dxz_Px_F3z_Py;
  abcd[693] = I_ERI_F3y_Px_F3z_Py+ABY*I_ERI_D2y_Px_F3z_Py;
  abcd[694] = I_ERI_F2yz_Px_F3z_Py+ABY*I_ERI_Dyz_Px_F3z_Py;
  abcd[695] = I_ERI_Fy2z_Px_F3z_Py+ABY*I_ERI_D2z_Px_F3z_Py;
  abcd[696] = I_ERI_F2xz_Px_F3z_Py+ABZ*I_ERI_D2x_Px_F3z_Py;
  abcd[697] = I_ERI_Fxyz_Px_F3z_Py+ABZ*I_ERI_Dxy_Px_F3z_Py;
  abcd[698] = I_ERI_Fx2z_Px_F3z_Py+ABZ*I_ERI_Dxz_Px_F3z_Py;
  abcd[699] = I_ERI_F2yz_Px_F3z_Py+ABZ*I_ERI_D2y_Px_F3z_Py;
  abcd[700] = I_ERI_Fy2z_Px_F3z_Py+ABZ*I_ERI_Dyz_Px_F3z_Py;
  abcd[701] = I_ERI_F3z_Px_F3z_Py+ABZ*I_ERI_D2z_Px_F3z_Py;
  abcd[702] = I_ERI_F2xy_Py_F3z_Py+ABY*I_ERI_D2x_Py_F3z_Py;
  abcd[703] = I_ERI_Fx2y_Py_F3z_Py+ABY*I_ERI_Dxy_Py_F3z_Py;
  abcd[704] = I_ERI_Fxyz_Py_F3z_Py+ABY*I_ERI_Dxz_Py_F3z_Py;
  abcd[705] = I_ERI_F3y_Py_F3z_Py+ABY*I_ERI_D2y_Py_F3z_Py;
  abcd[706] = I_ERI_F2yz_Py_F3z_Py+ABY*I_ERI_Dyz_Py_F3z_Py;
  abcd[707] = I_ERI_Fy2z_Py_F3z_Py+ABY*I_ERI_D2z_Py_F3z_Py;
  abcd[708] = I_ERI_F2xz_Py_F3z_Py+ABZ*I_ERI_D2x_Py_F3z_Py;
  abcd[709] = I_ERI_Fxyz_Py_F3z_Py+ABZ*I_ERI_Dxy_Py_F3z_Py;
  abcd[710] = I_ERI_Fx2z_Py_F3z_Py+ABZ*I_ERI_Dxz_Py_F3z_Py;
  abcd[711] = I_ERI_F2yz_Py_F3z_Py+ABZ*I_ERI_D2y_Py_F3z_Py;
  abcd[712] = I_ERI_Fy2z_Py_F3z_Py+ABZ*I_ERI_Dyz_Py_F3z_Py;
  abcd[713] = I_ERI_F3z_Py_F3z_Py+ABZ*I_ERI_D2z_Py_F3z_Py;
  abcd[714] = I_ERI_F2xz_Pz_F3z_Py+ABZ*I_ERI_D2x_Pz_F3z_Py;
  abcd[715] = I_ERI_Fxyz_Pz_F3z_Py+ABZ*I_ERI_Dxy_Pz_F3z_Py;
  abcd[716] = I_ERI_Fx2z_Pz_F3z_Py+ABZ*I_ERI_Dxz_Pz_F3z_Py;
  abcd[717] = I_ERI_F2yz_Pz_F3z_Py+ABZ*I_ERI_D2y_Pz_F3z_Py;
  abcd[718] = I_ERI_Fy2z_Pz_F3z_Py+ABZ*I_ERI_Dyz_Pz_F3z_Py;
  abcd[719] = I_ERI_F3z_Pz_F3z_Py+ABZ*I_ERI_D2z_Pz_F3z_Py;
  abcd[720] = I_ERI_F3x_Px_F3x_Pz+ABX*I_ERI_D2x_Px_F3x_Pz;
  abcd[721] = I_ERI_F2xy_Px_F3x_Pz+ABX*I_ERI_Dxy_Px_F3x_Pz;
  abcd[722] = I_ERI_F2xz_Px_F3x_Pz+ABX*I_ERI_Dxz_Px_F3x_Pz;
  abcd[723] = I_ERI_Fx2y_Px_F3x_Pz+ABX*I_ERI_D2y_Px_F3x_Pz;
  abcd[724] = I_ERI_Fxyz_Px_F3x_Pz+ABX*I_ERI_Dyz_Px_F3x_Pz;
  abcd[725] = I_ERI_Fx2z_Px_F3x_Pz+ABX*I_ERI_D2z_Px_F3x_Pz;
  abcd[726] = I_ERI_F2xy_Px_F3x_Pz+ABY*I_ERI_D2x_Px_F3x_Pz;
  abcd[727] = I_ERI_Fx2y_Px_F3x_Pz+ABY*I_ERI_Dxy_Px_F3x_Pz;
  abcd[728] = I_ERI_Fxyz_Px_F3x_Pz+ABY*I_ERI_Dxz_Px_F3x_Pz;
  abcd[729] = I_ERI_F3y_Px_F3x_Pz+ABY*I_ERI_D2y_Px_F3x_Pz;
  abcd[730] = I_ERI_F2yz_Px_F3x_Pz+ABY*I_ERI_Dyz_Px_F3x_Pz;
  abcd[731] = I_ERI_Fy2z_Px_F3x_Pz+ABY*I_ERI_D2z_Px_F3x_Pz;
  abcd[732] = I_ERI_F2xz_Px_F3x_Pz+ABZ*I_ERI_D2x_Px_F3x_Pz;
  abcd[733] = I_ERI_Fxyz_Px_F3x_Pz+ABZ*I_ERI_Dxy_Px_F3x_Pz;
  abcd[734] = I_ERI_Fx2z_Px_F3x_Pz+ABZ*I_ERI_Dxz_Px_F3x_Pz;
  abcd[735] = I_ERI_F2yz_Px_F3x_Pz+ABZ*I_ERI_D2y_Px_F3x_Pz;
  abcd[736] = I_ERI_Fy2z_Px_F3x_Pz+ABZ*I_ERI_Dyz_Px_F3x_Pz;
  abcd[737] = I_ERI_F3z_Px_F3x_Pz+ABZ*I_ERI_D2z_Px_F3x_Pz;
  abcd[738] = I_ERI_F2xy_Py_F3x_Pz+ABY*I_ERI_D2x_Py_F3x_Pz;
  abcd[739] = I_ERI_Fx2y_Py_F3x_Pz+ABY*I_ERI_Dxy_Py_F3x_Pz;
  abcd[740] = I_ERI_Fxyz_Py_F3x_Pz+ABY*I_ERI_Dxz_Py_F3x_Pz;
  abcd[741] = I_ERI_F3y_Py_F3x_Pz+ABY*I_ERI_D2y_Py_F3x_Pz;
  abcd[742] = I_ERI_F2yz_Py_F3x_Pz+ABY*I_ERI_Dyz_Py_F3x_Pz;
  abcd[743] = I_ERI_Fy2z_Py_F3x_Pz+ABY*I_ERI_D2z_Py_F3x_Pz;
  abcd[744] = I_ERI_F2xz_Py_F3x_Pz+ABZ*I_ERI_D2x_Py_F3x_Pz;
  abcd[745] = I_ERI_Fxyz_Py_F3x_Pz+ABZ*I_ERI_Dxy_Py_F3x_Pz;
  abcd[746] = I_ERI_Fx2z_Py_F3x_Pz+ABZ*I_ERI_Dxz_Py_F3x_Pz;
  abcd[747] = I_ERI_F2yz_Py_F3x_Pz+ABZ*I_ERI_D2y_Py_F3x_Pz;
  abcd[748] = I_ERI_Fy2z_Py_F3x_Pz+ABZ*I_ERI_Dyz_Py_F3x_Pz;
  abcd[749] = I_ERI_F3z_Py_F3x_Pz+ABZ*I_ERI_D2z_Py_F3x_Pz;
  abcd[750] = I_ERI_F2xz_Pz_F3x_Pz+ABZ*I_ERI_D2x_Pz_F3x_Pz;
  abcd[751] = I_ERI_Fxyz_Pz_F3x_Pz+ABZ*I_ERI_Dxy_Pz_F3x_Pz;
  abcd[752] = I_ERI_Fx2z_Pz_F3x_Pz+ABZ*I_ERI_Dxz_Pz_F3x_Pz;
  abcd[753] = I_ERI_F2yz_Pz_F3x_Pz+ABZ*I_ERI_D2y_Pz_F3x_Pz;
  abcd[754] = I_ERI_Fy2z_Pz_F3x_Pz+ABZ*I_ERI_Dyz_Pz_F3x_Pz;
  abcd[755] = I_ERI_F3z_Pz_F3x_Pz+ABZ*I_ERI_D2z_Pz_F3x_Pz;
  abcd[756] = I_ERI_F3x_Px_F2xy_Pz+ABX*I_ERI_D2x_Px_F2xy_Pz;
  abcd[757] = I_ERI_F2xy_Px_F2xy_Pz+ABX*I_ERI_Dxy_Px_F2xy_Pz;
  abcd[758] = I_ERI_F2xz_Px_F2xy_Pz+ABX*I_ERI_Dxz_Px_F2xy_Pz;
  abcd[759] = I_ERI_Fx2y_Px_F2xy_Pz+ABX*I_ERI_D2y_Px_F2xy_Pz;
  abcd[760] = I_ERI_Fxyz_Px_F2xy_Pz+ABX*I_ERI_Dyz_Px_F2xy_Pz;
  abcd[761] = I_ERI_Fx2z_Px_F2xy_Pz+ABX*I_ERI_D2z_Px_F2xy_Pz;
  abcd[762] = I_ERI_F2xy_Px_F2xy_Pz+ABY*I_ERI_D2x_Px_F2xy_Pz;
  abcd[763] = I_ERI_Fx2y_Px_F2xy_Pz+ABY*I_ERI_Dxy_Px_F2xy_Pz;
  abcd[764] = I_ERI_Fxyz_Px_F2xy_Pz+ABY*I_ERI_Dxz_Px_F2xy_Pz;
  abcd[765] = I_ERI_F3y_Px_F2xy_Pz+ABY*I_ERI_D2y_Px_F2xy_Pz;
  abcd[766] = I_ERI_F2yz_Px_F2xy_Pz+ABY*I_ERI_Dyz_Px_F2xy_Pz;
  abcd[767] = I_ERI_Fy2z_Px_F2xy_Pz+ABY*I_ERI_D2z_Px_F2xy_Pz;
  abcd[768] = I_ERI_F2xz_Px_F2xy_Pz+ABZ*I_ERI_D2x_Px_F2xy_Pz;
  abcd[769] = I_ERI_Fxyz_Px_F2xy_Pz+ABZ*I_ERI_Dxy_Px_F2xy_Pz;
  abcd[770] = I_ERI_Fx2z_Px_F2xy_Pz+ABZ*I_ERI_Dxz_Px_F2xy_Pz;
  abcd[771] = I_ERI_F2yz_Px_F2xy_Pz+ABZ*I_ERI_D2y_Px_F2xy_Pz;
  abcd[772] = I_ERI_Fy2z_Px_F2xy_Pz+ABZ*I_ERI_Dyz_Px_F2xy_Pz;
  abcd[773] = I_ERI_F3z_Px_F2xy_Pz+ABZ*I_ERI_D2z_Px_F2xy_Pz;
  abcd[774] = I_ERI_F2xy_Py_F2xy_Pz+ABY*I_ERI_D2x_Py_F2xy_Pz;
  abcd[775] = I_ERI_Fx2y_Py_F2xy_Pz+ABY*I_ERI_Dxy_Py_F2xy_Pz;
  abcd[776] = I_ERI_Fxyz_Py_F2xy_Pz+ABY*I_ERI_Dxz_Py_F2xy_Pz;
  abcd[777] = I_ERI_F3y_Py_F2xy_Pz+ABY*I_ERI_D2y_Py_F2xy_Pz;
  abcd[778] = I_ERI_F2yz_Py_F2xy_Pz+ABY*I_ERI_Dyz_Py_F2xy_Pz;
  abcd[779] = I_ERI_Fy2z_Py_F2xy_Pz+ABY*I_ERI_D2z_Py_F2xy_Pz;
  abcd[780] = I_ERI_F2xz_Py_F2xy_Pz+ABZ*I_ERI_D2x_Py_F2xy_Pz;
  abcd[781] = I_ERI_Fxyz_Py_F2xy_Pz+ABZ*I_ERI_Dxy_Py_F2xy_Pz;
  abcd[782] = I_ERI_Fx2z_Py_F2xy_Pz+ABZ*I_ERI_Dxz_Py_F2xy_Pz;
  abcd[783] = I_ERI_F2yz_Py_F2xy_Pz+ABZ*I_ERI_D2y_Py_F2xy_Pz;
  abcd[784] = I_ERI_Fy2z_Py_F2xy_Pz+ABZ*I_ERI_Dyz_Py_F2xy_Pz;
  abcd[785] = I_ERI_F3z_Py_F2xy_Pz+ABZ*I_ERI_D2z_Py_F2xy_Pz;
  abcd[786] = I_ERI_F2xz_Pz_F2xy_Pz+ABZ*I_ERI_D2x_Pz_F2xy_Pz;
  abcd[787] = I_ERI_Fxyz_Pz_F2xy_Pz+ABZ*I_ERI_Dxy_Pz_F2xy_Pz;
  abcd[788] = I_ERI_Fx2z_Pz_F2xy_Pz+ABZ*I_ERI_Dxz_Pz_F2xy_Pz;
  abcd[789] = I_ERI_F2yz_Pz_F2xy_Pz+ABZ*I_ERI_D2y_Pz_F2xy_Pz;
  abcd[790] = I_ERI_Fy2z_Pz_F2xy_Pz+ABZ*I_ERI_Dyz_Pz_F2xy_Pz;
  abcd[791] = I_ERI_F3z_Pz_F2xy_Pz+ABZ*I_ERI_D2z_Pz_F2xy_Pz;
  abcd[792] = I_ERI_F3x_Px_F2xz_Pz+ABX*I_ERI_D2x_Px_F2xz_Pz;
  abcd[793] = I_ERI_F2xy_Px_F2xz_Pz+ABX*I_ERI_Dxy_Px_F2xz_Pz;
  abcd[794] = I_ERI_F2xz_Px_F2xz_Pz+ABX*I_ERI_Dxz_Px_F2xz_Pz;
  abcd[795] = I_ERI_Fx2y_Px_F2xz_Pz+ABX*I_ERI_D2y_Px_F2xz_Pz;
  abcd[796] = I_ERI_Fxyz_Px_F2xz_Pz+ABX*I_ERI_Dyz_Px_F2xz_Pz;
  abcd[797] = I_ERI_Fx2z_Px_F2xz_Pz+ABX*I_ERI_D2z_Px_F2xz_Pz;
  abcd[798] = I_ERI_F2xy_Px_F2xz_Pz+ABY*I_ERI_D2x_Px_F2xz_Pz;
  abcd[799] = I_ERI_Fx2y_Px_F2xz_Pz+ABY*I_ERI_Dxy_Px_F2xz_Pz;
  abcd[800] = I_ERI_Fxyz_Px_F2xz_Pz+ABY*I_ERI_Dxz_Px_F2xz_Pz;
  abcd[801] = I_ERI_F3y_Px_F2xz_Pz+ABY*I_ERI_D2y_Px_F2xz_Pz;
  abcd[802] = I_ERI_F2yz_Px_F2xz_Pz+ABY*I_ERI_Dyz_Px_F2xz_Pz;
  abcd[803] = I_ERI_Fy2z_Px_F2xz_Pz+ABY*I_ERI_D2z_Px_F2xz_Pz;
  abcd[804] = I_ERI_F2xz_Px_F2xz_Pz+ABZ*I_ERI_D2x_Px_F2xz_Pz;
  abcd[805] = I_ERI_Fxyz_Px_F2xz_Pz+ABZ*I_ERI_Dxy_Px_F2xz_Pz;
  abcd[806] = I_ERI_Fx2z_Px_F2xz_Pz+ABZ*I_ERI_Dxz_Px_F2xz_Pz;
  abcd[807] = I_ERI_F2yz_Px_F2xz_Pz+ABZ*I_ERI_D2y_Px_F2xz_Pz;
  abcd[808] = I_ERI_Fy2z_Px_F2xz_Pz+ABZ*I_ERI_Dyz_Px_F2xz_Pz;
  abcd[809] = I_ERI_F3z_Px_F2xz_Pz+ABZ*I_ERI_D2z_Px_F2xz_Pz;
  abcd[810] = I_ERI_F2xy_Py_F2xz_Pz+ABY*I_ERI_D2x_Py_F2xz_Pz;
  abcd[811] = I_ERI_Fx2y_Py_F2xz_Pz+ABY*I_ERI_Dxy_Py_F2xz_Pz;
  abcd[812] = I_ERI_Fxyz_Py_F2xz_Pz+ABY*I_ERI_Dxz_Py_F2xz_Pz;
  abcd[813] = I_ERI_F3y_Py_F2xz_Pz+ABY*I_ERI_D2y_Py_F2xz_Pz;
  abcd[814] = I_ERI_F2yz_Py_F2xz_Pz+ABY*I_ERI_Dyz_Py_F2xz_Pz;
  abcd[815] = I_ERI_Fy2z_Py_F2xz_Pz+ABY*I_ERI_D2z_Py_F2xz_Pz;
  abcd[816] = I_ERI_F2xz_Py_F2xz_Pz+ABZ*I_ERI_D2x_Py_F2xz_Pz;
  abcd[817] = I_ERI_Fxyz_Py_F2xz_Pz+ABZ*I_ERI_Dxy_Py_F2xz_Pz;
  abcd[818] = I_ERI_Fx2z_Py_F2xz_Pz+ABZ*I_ERI_Dxz_Py_F2xz_Pz;
  abcd[819] = I_ERI_F2yz_Py_F2xz_Pz+ABZ*I_ERI_D2y_Py_F2xz_Pz;
  abcd[820] = I_ERI_Fy2z_Py_F2xz_Pz+ABZ*I_ERI_Dyz_Py_F2xz_Pz;
  abcd[821] = I_ERI_F3z_Py_F2xz_Pz+ABZ*I_ERI_D2z_Py_F2xz_Pz;
  abcd[822] = I_ERI_F2xz_Pz_F2xz_Pz+ABZ*I_ERI_D2x_Pz_F2xz_Pz;
  abcd[823] = I_ERI_Fxyz_Pz_F2xz_Pz+ABZ*I_ERI_Dxy_Pz_F2xz_Pz;
  abcd[824] = I_ERI_Fx2z_Pz_F2xz_Pz+ABZ*I_ERI_Dxz_Pz_F2xz_Pz;
  abcd[825] = I_ERI_F2yz_Pz_F2xz_Pz+ABZ*I_ERI_D2y_Pz_F2xz_Pz;
  abcd[826] = I_ERI_Fy2z_Pz_F2xz_Pz+ABZ*I_ERI_Dyz_Pz_F2xz_Pz;
  abcd[827] = I_ERI_F3z_Pz_F2xz_Pz+ABZ*I_ERI_D2z_Pz_F2xz_Pz;
  abcd[828] = I_ERI_F3x_Px_Fx2y_Pz+ABX*I_ERI_D2x_Px_Fx2y_Pz;
  abcd[829] = I_ERI_F2xy_Px_Fx2y_Pz+ABX*I_ERI_Dxy_Px_Fx2y_Pz;
  abcd[830] = I_ERI_F2xz_Px_Fx2y_Pz+ABX*I_ERI_Dxz_Px_Fx2y_Pz;
  abcd[831] = I_ERI_Fx2y_Px_Fx2y_Pz+ABX*I_ERI_D2y_Px_Fx2y_Pz;
  abcd[832] = I_ERI_Fxyz_Px_Fx2y_Pz+ABX*I_ERI_Dyz_Px_Fx2y_Pz;
  abcd[833] = I_ERI_Fx2z_Px_Fx2y_Pz+ABX*I_ERI_D2z_Px_Fx2y_Pz;
  abcd[834] = I_ERI_F2xy_Px_Fx2y_Pz+ABY*I_ERI_D2x_Px_Fx2y_Pz;
  abcd[835] = I_ERI_Fx2y_Px_Fx2y_Pz+ABY*I_ERI_Dxy_Px_Fx2y_Pz;
  abcd[836] = I_ERI_Fxyz_Px_Fx2y_Pz+ABY*I_ERI_Dxz_Px_Fx2y_Pz;
  abcd[837] = I_ERI_F3y_Px_Fx2y_Pz+ABY*I_ERI_D2y_Px_Fx2y_Pz;
  abcd[838] = I_ERI_F2yz_Px_Fx2y_Pz+ABY*I_ERI_Dyz_Px_Fx2y_Pz;
  abcd[839] = I_ERI_Fy2z_Px_Fx2y_Pz+ABY*I_ERI_D2z_Px_Fx2y_Pz;
  abcd[840] = I_ERI_F2xz_Px_Fx2y_Pz+ABZ*I_ERI_D2x_Px_Fx2y_Pz;
  abcd[841] = I_ERI_Fxyz_Px_Fx2y_Pz+ABZ*I_ERI_Dxy_Px_Fx2y_Pz;
  abcd[842] = I_ERI_Fx2z_Px_Fx2y_Pz+ABZ*I_ERI_Dxz_Px_Fx2y_Pz;
  abcd[843] = I_ERI_F2yz_Px_Fx2y_Pz+ABZ*I_ERI_D2y_Px_Fx2y_Pz;
  abcd[844] = I_ERI_Fy2z_Px_Fx2y_Pz+ABZ*I_ERI_Dyz_Px_Fx2y_Pz;
  abcd[845] = I_ERI_F3z_Px_Fx2y_Pz+ABZ*I_ERI_D2z_Px_Fx2y_Pz;
  abcd[846] = I_ERI_F2xy_Py_Fx2y_Pz+ABY*I_ERI_D2x_Py_Fx2y_Pz;
  abcd[847] = I_ERI_Fx2y_Py_Fx2y_Pz+ABY*I_ERI_Dxy_Py_Fx2y_Pz;
  abcd[848] = I_ERI_Fxyz_Py_Fx2y_Pz+ABY*I_ERI_Dxz_Py_Fx2y_Pz;
  abcd[849] = I_ERI_F3y_Py_Fx2y_Pz+ABY*I_ERI_D2y_Py_Fx2y_Pz;
  abcd[850] = I_ERI_F2yz_Py_Fx2y_Pz+ABY*I_ERI_Dyz_Py_Fx2y_Pz;
  abcd[851] = I_ERI_Fy2z_Py_Fx2y_Pz+ABY*I_ERI_D2z_Py_Fx2y_Pz;
  abcd[852] = I_ERI_F2xz_Py_Fx2y_Pz+ABZ*I_ERI_D2x_Py_Fx2y_Pz;
  abcd[853] = I_ERI_Fxyz_Py_Fx2y_Pz+ABZ*I_ERI_Dxy_Py_Fx2y_Pz;
  abcd[854] = I_ERI_Fx2z_Py_Fx2y_Pz+ABZ*I_ERI_Dxz_Py_Fx2y_Pz;
  abcd[855] = I_ERI_F2yz_Py_Fx2y_Pz+ABZ*I_ERI_D2y_Py_Fx2y_Pz;
  abcd[856] = I_ERI_Fy2z_Py_Fx2y_Pz+ABZ*I_ERI_Dyz_Py_Fx2y_Pz;
  abcd[857] = I_ERI_F3z_Py_Fx2y_Pz+ABZ*I_ERI_D2z_Py_Fx2y_Pz;
  abcd[858] = I_ERI_F2xz_Pz_Fx2y_Pz+ABZ*I_ERI_D2x_Pz_Fx2y_Pz;
  abcd[859] = I_ERI_Fxyz_Pz_Fx2y_Pz+ABZ*I_ERI_Dxy_Pz_Fx2y_Pz;
  abcd[860] = I_ERI_Fx2z_Pz_Fx2y_Pz+ABZ*I_ERI_Dxz_Pz_Fx2y_Pz;
  abcd[861] = I_ERI_F2yz_Pz_Fx2y_Pz+ABZ*I_ERI_D2y_Pz_Fx2y_Pz;
  abcd[862] = I_ERI_Fy2z_Pz_Fx2y_Pz+ABZ*I_ERI_Dyz_Pz_Fx2y_Pz;
  abcd[863] = I_ERI_F3z_Pz_Fx2y_Pz+ABZ*I_ERI_D2z_Pz_Fx2y_Pz;
  abcd[864] = I_ERI_F3x_Px_Fxyz_Pz+ABX*I_ERI_D2x_Px_Fxyz_Pz;
  abcd[865] = I_ERI_F2xy_Px_Fxyz_Pz+ABX*I_ERI_Dxy_Px_Fxyz_Pz;
  abcd[866] = I_ERI_F2xz_Px_Fxyz_Pz+ABX*I_ERI_Dxz_Px_Fxyz_Pz;
  abcd[867] = I_ERI_Fx2y_Px_Fxyz_Pz+ABX*I_ERI_D2y_Px_Fxyz_Pz;
  abcd[868] = I_ERI_Fxyz_Px_Fxyz_Pz+ABX*I_ERI_Dyz_Px_Fxyz_Pz;
  abcd[869] = I_ERI_Fx2z_Px_Fxyz_Pz+ABX*I_ERI_D2z_Px_Fxyz_Pz;
  abcd[870] = I_ERI_F2xy_Px_Fxyz_Pz+ABY*I_ERI_D2x_Px_Fxyz_Pz;
  abcd[871] = I_ERI_Fx2y_Px_Fxyz_Pz+ABY*I_ERI_Dxy_Px_Fxyz_Pz;
  abcd[872] = I_ERI_Fxyz_Px_Fxyz_Pz+ABY*I_ERI_Dxz_Px_Fxyz_Pz;
  abcd[873] = I_ERI_F3y_Px_Fxyz_Pz+ABY*I_ERI_D2y_Px_Fxyz_Pz;
  abcd[874] = I_ERI_F2yz_Px_Fxyz_Pz+ABY*I_ERI_Dyz_Px_Fxyz_Pz;
  abcd[875] = I_ERI_Fy2z_Px_Fxyz_Pz+ABY*I_ERI_D2z_Px_Fxyz_Pz;
  abcd[876] = I_ERI_F2xz_Px_Fxyz_Pz+ABZ*I_ERI_D2x_Px_Fxyz_Pz;
  abcd[877] = I_ERI_Fxyz_Px_Fxyz_Pz+ABZ*I_ERI_Dxy_Px_Fxyz_Pz;
  abcd[878] = I_ERI_Fx2z_Px_Fxyz_Pz+ABZ*I_ERI_Dxz_Px_Fxyz_Pz;
  abcd[879] = I_ERI_F2yz_Px_Fxyz_Pz+ABZ*I_ERI_D2y_Px_Fxyz_Pz;
  abcd[880] = I_ERI_Fy2z_Px_Fxyz_Pz+ABZ*I_ERI_Dyz_Px_Fxyz_Pz;
  abcd[881] = I_ERI_F3z_Px_Fxyz_Pz+ABZ*I_ERI_D2z_Px_Fxyz_Pz;
  abcd[882] = I_ERI_F2xy_Py_Fxyz_Pz+ABY*I_ERI_D2x_Py_Fxyz_Pz;
  abcd[883] = I_ERI_Fx2y_Py_Fxyz_Pz+ABY*I_ERI_Dxy_Py_Fxyz_Pz;
  abcd[884] = I_ERI_Fxyz_Py_Fxyz_Pz+ABY*I_ERI_Dxz_Py_Fxyz_Pz;
  abcd[885] = I_ERI_F3y_Py_Fxyz_Pz+ABY*I_ERI_D2y_Py_Fxyz_Pz;
  abcd[886] = I_ERI_F2yz_Py_Fxyz_Pz+ABY*I_ERI_Dyz_Py_Fxyz_Pz;
  abcd[887] = I_ERI_Fy2z_Py_Fxyz_Pz+ABY*I_ERI_D2z_Py_Fxyz_Pz;
  abcd[888] = I_ERI_F2xz_Py_Fxyz_Pz+ABZ*I_ERI_D2x_Py_Fxyz_Pz;
  abcd[889] = I_ERI_Fxyz_Py_Fxyz_Pz+ABZ*I_ERI_Dxy_Py_Fxyz_Pz;
  abcd[890] = I_ERI_Fx2z_Py_Fxyz_Pz+ABZ*I_ERI_Dxz_Py_Fxyz_Pz;
  abcd[891] = I_ERI_F2yz_Py_Fxyz_Pz+ABZ*I_ERI_D2y_Py_Fxyz_Pz;
  abcd[892] = I_ERI_Fy2z_Py_Fxyz_Pz+ABZ*I_ERI_Dyz_Py_Fxyz_Pz;
  abcd[893] = I_ERI_F3z_Py_Fxyz_Pz+ABZ*I_ERI_D2z_Py_Fxyz_Pz;
  abcd[894] = I_ERI_F2xz_Pz_Fxyz_Pz+ABZ*I_ERI_D2x_Pz_Fxyz_Pz;
  abcd[895] = I_ERI_Fxyz_Pz_Fxyz_Pz+ABZ*I_ERI_Dxy_Pz_Fxyz_Pz;
  abcd[896] = I_ERI_Fx2z_Pz_Fxyz_Pz+ABZ*I_ERI_Dxz_Pz_Fxyz_Pz;
  abcd[897] = I_ERI_F2yz_Pz_Fxyz_Pz+ABZ*I_ERI_D2y_Pz_Fxyz_Pz;
  abcd[898] = I_ERI_Fy2z_Pz_Fxyz_Pz+ABZ*I_ERI_Dyz_Pz_Fxyz_Pz;
  abcd[899] = I_ERI_F3z_Pz_Fxyz_Pz+ABZ*I_ERI_D2z_Pz_Fxyz_Pz;
  abcd[900] = I_ERI_F3x_Px_Fx2z_Pz+ABX*I_ERI_D2x_Px_Fx2z_Pz;
  abcd[901] = I_ERI_F2xy_Px_Fx2z_Pz+ABX*I_ERI_Dxy_Px_Fx2z_Pz;
  abcd[902] = I_ERI_F2xz_Px_Fx2z_Pz+ABX*I_ERI_Dxz_Px_Fx2z_Pz;
  abcd[903] = I_ERI_Fx2y_Px_Fx2z_Pz+ABX*I_ERI_D2y_Px_Fx2z_Pz;
  abcd[904] = I_ERI_Fxyz_Px_Fx2z_Pz+ABX*I_ERI_Dyz_Px_Fx2z_Pz;
  abcd[905] = I_ERI_Fx2z_Px_Fx2z_Pz+ABX*I_ERI_D2z_Px_Fx2z_Pz;
  abcd[906] = I_ERI_F2xy_Px_Fx2z_Pz+ABY*I_ERI_D2x_Px_Fx2z_Pz;
  abcd[907] = I_ERI_Fx2y_Px_Fx2z_Pz+ABY*I_ERI_Dxy_Px_Fx2z_Pz;
  abcd[908] = I_ERI_Fxyz_Px_Fx2z_Pz+ABY*I_ERI_Dxz_Px_Fx2z_Pz;
  abcd[909] = I_ERI_F3y_Px_Fx2z_Pz+ABY*I_ERI_D2y_Px_Fx2z_Pz;
  abcd[910] = I_ERI_F2yz_Px_Fx2z_Pz+ABY*I_ERI_Dyz_Px_Fx2z_Pz;
  abcd[911] = I_ERI_Fy2z_Px_Fx2z_Pz+ABY*I_ERI_D2z_Px_Fx2z_Pz;
  abcd[912] = I_ERI_F2xz_Px_Fx2z_Pz+ABZ*I_ERI_D2x_Px_Fx2z_Pz;
  abcd[913] = I_ERI_Fxyz_Px_Fx2z_Pz+ABZ*I_ERI_Dxy_Px_Fx2z_Pz;
  abcd[914] = I_ERI_Fx2z_Px_Fx2z_Pz+ABZ*I_ERI_Dxz_Px_Fx2z_Pz;
  abcd[915] = I_ERI_F2yz_Px_Fx2z_Pz+ABZ*I_ERI_D2y_Px_Fx2z_Pz;
  abcd[916] = I_ERI_Fy2z_Px_Fx2z_Pz+ABZ*I_ERI_Dyz_Px_Fx2z_Pz;
  abcd[917] = I_ERI_F3z_Px_Fx2z_Pz+ABZ*I_ERI_D2z_Px_Fx2z_Pz;
  abcd[918] = I_ERI_F2xy_Py_Fx2z_Pz+ABY*I_ERI_D2x_Py_Fx2z_Pz;
  abcd[919] = I_ERI_Fx2y_Py_Fx2z_Pz+ABY*I_ERI_Dxy_Py_Fx2z_Pz;
  abcd[920] = I_ERI_Fxyz_Py_Fx2z_Pz+ABY*I_ERI_Dxz_Py_Fx2z_Pz;
  abcd[921] = I_ERI_F3y_Py_Fx2z_Pz+ABY*I_ERI_D2y_Py_Fx2z_Pz;
  abcd[922] = I_ERI_F2yz_Py_Fx2z_Pz+ABY*I_ERI_Dyz_Py_Fx2z_Pz;
  abcd[923] = I_ERI_Fy2z_Py_Fx2z_Pz+ABY*I_ERI_D2z_Py_Fx2z_Pz;
  abcd[924] = I_ERI_F2xz_Py_Fx2z_Pz+ABZ*I_ERI_D2x_Py_Fx2z_Pz;
  abcd[925] = I_ERI_Fxyz_Py_Fx2z_Pz+ABZ*I_ERI_Dxy_Py_Fx2z_Pz;
  abcd[926] = I_ERI_Fx2z_Py_Fx2z_Pz+ABZ*I_ERI_Dxz_Py_Fx2z_Pz;
  abcd[927] = I_ERI_F2yz_Py_Fx2z_Pz+ABZ*I_ERI_D2y_Py_Fx2z_Pz;
  abcd[928] = I_ERI_Fy2z_Py_Fx2z_Pz+ABZ*I_ERI_Dyz_Py_Fx2z_Pz;
  abcd[929] = I_ERI_F3z_Py_Fx2z_Pz+ABZ*I_ERI_D2z_Py_Fx2z_Pz;
  abcd[930] = I_ERI_F2xz_Pz_Fx2z_Pz+ABZ*I_ERI_D2x_Pz_Fx2z_Pz;
  abcd[931] = I_ERI_Fxyz_Pz_Fx2z_Pz+ABZ*I_ERI_Dxy_Pz_Fx2z_Pz;
  abcd[932] = I_ERI_Fx2z_Pz_Fx2z_Pz+ABZ*I_ERI_Dxz_Pz_Fx2z_Pz;
  abcd[933] = I_ERI_F2yz_Pz_Fx2z_Pz+ABZ*I_ERI_D2y_Pz_Fx2z_Pz;
  abcd[934] = I_ERI_Fy2z_Pz_Fx2z_Pz+ABZ*I_ERI_Dyz_Pz_Fx2z_Pz;
  abcd[935] = I_ERI_F3z_Pz_Fx2z_Pz+ABZ*I_ERI_D2z_Pz_Fx2z_Pz;
  abcd[936] = I_ERI_F3x_Px_F3y_Pz+ABX*I_ERI_D2x_Px_F3y_Pz;
  abcd[937] = I_ERI_F2xy_Px_F3y_Pz+ABX*I_ERI_Dxy_Px_F3y_Pz;
  abcd[938] = I_ERI_F2xz_Px_F3y_Pz+ABX*I_ERI_Dxz_Px_F3y_Pz;
  abcd[939] = I_ERI_Fx2y_Px_F3y_Pz+ABX*I_ERI_D2y_Px_F3y_Pz;
  abcd[940] = I_ERI_Fxyz_Px_F3y_Pz+ABX*I_ERI_Dyz_Px_F3y_Pz;
  abcd[941] = I_ERI_Fx2z_Px_F3y_Pz+ABX*I_ERI_D2z_Px_F3y_Pz;
  abcd[942] = I_ERI_F2xy_Px_F3y_Pz+ABY*I_ERI_D2x_Px_F3y_Pz;
  abcd[943] = I_ERI_Fx2y_Px_F3y_Pz+ABY*I_ERI_Dxy_Px_F3y_Pz;
  abcd[944] = I_ERI_Fxyz_Px_F3y_Pz+ABY*I_ERI_Dxz_Px_F3y_Pz;
  abcd[945] = I_ERI_F3y_Px_F3y_Pz+ABY*I_ERI_D2y_Px_F3y_Pz;
  abcd[946] = I_ERI_F2yz_Px_F3y_Pz+ABY*I_ERI_Dyz_Px_F3y_Pz;
  abcd[947] = I_ERI_Fy2z_Px_F3y_Pz+ABY*I_ERI_D2z_Px_F3y_Pz;
  abcd[948] = I_ERI_F2xz_Px_F3y_Pz+ABZ*I_ERI_D2x_Px_F3y_Pz;
  abcd[949] = I_ERI_Fxyz_Px_F3y_Pz+ABZ*I_ERI_Dxy_Px_F3y_Pz;
  abcd[950] = I_ERI_Fx2z_Px_F3y_Pz+ABZ*I_ERI_Dxz_Px_F3y_Pz;
  abcd[951] = I_ERI_F2yz_Px_F3y_Pz+ABZ*I_ERI_D2y_Px_F3y_Pz;
  abcd[952] = I_ERI_Fy2z_Px_F3y_Pz+ABZ*I_ERI_Dyz_Px_F3y_Pz;
  abcd[953] = I_ERI_F3z_Px_F3y_Pz+ABZ*I_ERI_D2z_Px_F3y_Pz;
  abcd[954] = I_ERI_F2xy_Py_F3y_Pz+ABY*I_ERI_D2x_Py_F3y_Pz;
  abcd[955] = I_ERI_Fx2y_Py_F3y_Pz+ABY*I_ERI_Dxy_Py_F3y_Pz;
  abcd[956] = I_ERI_Fxyz_Py_F3y_Pz+ABY*I_ERI_Dxz_Py_F3y_Pz;
  abcd[957] = I_ERI_F3y_Py_F3y_Pz+ABY*I_ERI_D2y_Py_F3y_Pz;
  abcd[958] = I_ERI_F2yz_Py_F3y_Pz+ABY*I_ERI_Dyz_Py_F3y_Pz;
  abcd[959] = I_ERI_Fy2z_Py_F3y_Pz+ABY*I_ERI_D2z_Py_F3y_Pz;
  abcd[960] = I_ERI_F2xz_Py_F3y_Pz+ABZ*I_ERI_D2x_Py_F3y_Pz;
  abcd[961] = I_ERI_Fxyz_Py_F3y_Pz+ABZ*I_ERI_Dxy_Py_F3y_Pz;
  abcd[962] = I_ERI_Fx2z_Py_F3y_Pz+ABZ*I_ERI_Dxz_Py_F3y_Pz;
  abcd[963] = I_ERI_F2yz_Py_F3y_Pz+ABZ*I_ERI_D2y_Py_F3y_Pz;
  abcd[964] = I_ERI_Fy2z_Py_F3y_Pz+ABZ*I_ERI_Dyz_Py_F3y_Pz;
  abcd[965] = I_ERI_F3z_Py_F3y_Pz+ABZ*I_ERI_D2z_Py_F3y_Pz;
  abcd[966] = I_ERI_F2xz_Pz_F3y_Pz+ABZ*I_ERI_D2x_Pz_F3y_Pz;
  abcd[967] = I_ERI_Fxyz_Pz_F3y_Pz+ABZ*I_ERI_Dxy_Pz_F3y_Pz;
  abcd[968] = I_ERI_Fx2z_Pz_F3y_Pz+ABZ*I_ERI_Dxz_Pz_F3y_Pz;
  abcd[969] = I_ERI_F2yz_Pz_F3y_Pz+ABZ*I_ERI_D2y_Pz_F3y_Pz;
  abcd[970] = I_ERI_Fy2z_Pz_F3y_Pz+ABZ*I_ERI_Dyz_Pz_F3y_Pz;
  abcd[971] = I_ERI_F3z_Pz_F3y_Pz+ABZ*I_ERI_D2z_Pz_F3y_Pz;
  abcd[972] = I_ERI_F3x_Px_F2yz_Pz+ABX*I_ERI_D2x_Px_F2yz_Pz;
  abcd[973] = I_ERI_F2xy_Px_F2yz_Pz+ABX*I_ERI_Dxy_Px_F2yz_Pz;
  abcd[974] = I_ERI_F2xz_Px_F2yz_Pz+ABX*I_ERI_Dxz_Px_F2yz_Pz;
  abcd[975] = I_ERI_Fx2y_Px_F2yz_Pz+ABX*I_ERI_D2y_Px_F2yz_Pz;
  abcd[976] = I_ERI_Fxyz_Px_F2yz_Pz+ABX*I_ERI_Dyz_Px_F2yz_Pz;
  abcd[977] = I_ERI_Fx2z_Px_F2yz_Pz+ABX*I_ERI_D2z_Px_F2yz_Pz;
  abcd[978] = I_ERI_F2xy_Px_F2yz_Pz+ABY*I_ERI_D2x_Px_F2yz_Pz;
  abcd[979] = I_ERI_Fx2y_Px_F2yz_Pz+ABY*I_ERI_Dxy_Px_F2yz_Pz;
  abcd[980] = I_ERI_Fxyz_Px_F2yz_Pz+ABY*I_ERI_Dxz_Px_F2yz_Pz;
  abcd[981] = I_ERI_F3y_Px_F2yz_Pz+ABY*I_ERI_D2y_Px_F2yz_Pz;
  abcd[982] = I_ERI_F2yz_Px_F2yz_Pz+ABY*I_ERI_Dyz_Px_F2yz_Pz;
  abcd[983] = I_ERI_Fy2z_Px_F2yz_Pz+ABY*I_ERI_D2z_Px_F2yz_Pz;
  abcd[984] = I_ERI_F2xz_Px_F2yz_Pz+ABZ*I_ERI_D2x_Px_F2yz_Pz;
  abcd[985] = I_ERI_Fxyz_Px_F2yz_Pz+ABZ*I_ERI_Dxy_Px_F2yz_Pz;
  abcd[986] = I_ERI_Fx2z_Px_F2yz_Pz+ABZ*I_ERI_Dxz_Px_F2yz_Pz;
  abcd[987] = I_ERI_F2yz_Px_F2yz_Pz+ABZ*I_ERI_D2y_Px_F2yz_Pz;
  abcd[988] = I_ERI_Fy2z_Px_F2yz_Pz+ABZ*I_ERI_Dyz_Px_F2yz_Pz;
  abcd[989] = I_ERI_F3z_Px_F2yz_Pz+ABZ*I_ERI_D2z_Px_F2yz_Pz;
  abcd[990] = I_ERI_F2xy_Py_F2yz_Pz+ABY*I_ERI_D2x_Py_F2yz_Pz;
  abcd[991] = I_ERI_Fx2y_Py_F2yz_Pz+ABY*I_ERI_Dxy_Py_F2yz_Pz;
  abcd[992] = I_ERI_Fxyz_Py_F2yz_Pz+ABY*I_ERI_Dxz_Py_F2yz_Pz;
  abcd[993] = I_ERI_F3y_Py_F2yz_Pz+ABY*I_ERI_D2y_Py_F2yz_Pz;
  abcd[994] = I_ERI_F2yz_Py_F2yz_Pz+ABY*I_ERI_Dyz_Py_F2yz_Pz;
  abcd[995] = I_ERI_Fy2z_Py_F2yz_Pz+ABY*I_ERI_D2z_Py_F2yz_Pz;
  abcd[996] = I_ERI_F2xz_Py_F2yz_Pz+ABZ*I_ERI_D2x_Py_F2yz_Pz;
  abcd[997] = I_ERI_Fxyz_Py_F2yz_Pz+ABZ*I_ERI_Dxy_Py_F2yz_Pz;
  abcd[998] = I_ERI_Fx2z_Py_F2yz_Pz+ABZ*I_ERI_Dxz_Py_F2yz_Pz;
  abcd[999] = I_ERI_F2yz_Py_F2yz_Pz+ABZ*I_ERI_D2y_Py_F2yz_Pz;
  abcd[1000] = I_ERI_Fy2z_Py_F2yz_Pz+ABZ*I_ERI_Dyz_Py_F2yz_Pz;
  abcd[1001] = I_ERI_F3z_Py_F2yz_Pz+ABZ*I_ERI_D2z_Py_F2yz_Pz;
  abcd[1002] = I_ERI_F2xz_Pz_F2yz_Pz+ABZ*I_ERI_D2x_Pz_F2yz_Pz;
  abcd[1003] = I_ERI_Fxyz_Pz_F2yz_Pz+ABZ*I_ERI_Dxy_Pz_F2yz_Pz;
  abcd[1004] = I_ERI_Fx2z_Pz_F2yz_Pz+ABZ*I_ERI_Dxz_Pz_F2yz_Pz;
  abcd[1005] = I_ERI_F2yz_Pz_F2yz_Pz+ABZ*I_ERI_D2y_Pz_F2yz_Pz;
  abcd[1006] = I_ERI_Fy2z_Pz_F2yz_Pz+ABZ*I_ERI_Dyz_Pz_F2yz_Pz;
  abcd[1007] = I_ERI_F3z_Pz_F2yz_Pz+ABZ*I_ERI_D2z_Pz_F2yz_Pz;
  abcd[1008] = I_ERI_F3x_Px_Fy2z_Pz+ABX*I_ERI_D2x_Px_Fy2z_Pz;
  abcd[1009] = I_ERI_F2xy_Px_Fy2z_Pz+ABX*I_ERI_Dxy_Px_Fy2z_Pz;
  abcd[1010] = I_ERI_F2xz_Px_Fy2z_Pz+ABX*I_ERI_Dxz_Px_Fy2z_Pz;
  abcd[1011] = I_ERI_Fx2y_Px_Fy2z_Pz+ABX*I_ERI_D2y_Px_Fy2z_Pz;
  abcd[1012] = I_ERI_Fxyz_Px_Fy2z_Pz+ABX*I_ERI_Dyz_Px_Fy2z_Pz;
  abcd[1013] = I_ERI_Fx2z_Px_Fy2z_Pz+ABX*I_ERI_D2z_Px_Fy2z_Pz;
  abcd[1014] = I_ERI_F2xy_Px_Fy2z_Pz+ABY*I_ERI_D2x_Px_Fy2z_Pz;
  abcd[1015] = I_ERI_Fx2y_Px_Fy2z_Pz+ABY*I_ERI_Dxy_Px_Fy2z_Pz;
  abcd[1016] = I_ERI_Fxyz_Px_Fy2z_Pz+ABY*I_ERI_Dxz_Px_Fy2z_Pz;
  abcd[1017] = I_ERI_F3y_Px_Fy2z_Pz+ABY*I_ERI_D2y_Px_Fy2z_Pz;
  abcd[1018] = I_ERI_F2yz_Px_Fy2z_Pz+ABY*I_ERI_Dyz_Px_Fy2z_Pz;
  abcd[1019] = I_ERI_Fy2z_Px_Fy2z_Pz+ABY*I_ERI_D2z_Px_Fy2z_Pz;
  abcd[1020] = I_ERI_F2xz_Px_Fy2z_Pz+ABZ*I_ERI_D2x_Px_Fy2z_Pz;
  abcd[1021] = I_ERI_Fxyz_Px_Fy2z_Pz+ABZ*I_ERI_Dxy_Px_Fy2z_Pz;
  abcd[1022] = I_ERI_Fx2z_Px_Fy2z_Pz+ABZ*I_ERI_Dxz_Px_Fy2z_Pz;
  abcd[1023] = I_ERI_F2yz_Px_Fy2z_Pz+ABZ*I_ERI_D2y_Px_Fy2z_Pz;
  abcd[1024] = I_ERI_Fy2z_Px_Fy2z_Pz+ABZ*I_ERI_Dyz_Px_Fy2z_Pz;
  abcd[1025] = I_ERI_F3z_Px_Fy2z_Pz+ABZ*I_ERI_D2z_Px_Fy2z_Pz;
  abcd[1026] = I_ERI_F2xy_Py_Fy2z_Pz+ABY*I_ERI_D2x_Py_Fy2z_Pz;
  abcd[1027] = I_ERI_Fx2y_Py_Fy2z_Pz+ABY*I_ERI_Dxy_Py_Fy2z_Pz;
  abcd[1028] = I_ERI_Fxyz_Py_Fy2z_Pz+ABY*I_ERI_Dxz_Py_Fy2z_Pz;
  abcd[1029] = I_ERI_F3y_Py_Fy2z_Pz+ABY*I_ERI_D2y_Py_Fy2z_Pz;
  abcd[1030] = I_ERI_F2yz_Py_Fy2z_Pz+ABY*I_ERI_Dyz_Py_Fy2z_Pz;
  abcd[1031] = I_ERI_Fy2z_Py_Fy2z_Pz+ABY*I_ERI_D2z_Py_Fy2z_Pz;
  abcd[1032] = I_ERI_F2xz_Py_Fy2z_Pz+ABZ*I_ERI_D2x_Py_Fy2z_Pz;
  abcd[1033] = I_ERI_Fxyz_Py_Fy2z_Pz+ABZ*I_ERI_Dxy_Py_Fy2z_Pz;
  abcd[1034] = I_ERI_Fx2z_Py_Fy2z_Pz+ABZ*I_ERI_Dxz_Py_Fy2z_Pz;
  abcd[1035] = I_ERI_F2yz_Py_Fy2z_Pz+ABZ*I_ERI_D2y_Py_Fy2z_Pz;
  abcd[1036] = I_ERI_Fy2z_Py_Fy2z_Pz+ABZ*I_ERI_Dyz_Py_Fy2z_Pz;
  abcd[1037] = I_ERI_F3z_Py_Fy2z_Pz+ABZ*I_ERI_D2z_Py_Fy2z_Pz;
  abcd[1038] = I_ERI_F2xz_Pz_Fy2z_Pz+ABZ*I_ERI_D2x_Pz_Fy2z_Pz;
  abcd[1039] = I_ERI_Fxyz_Pz_Fy2z_Pz+ABZ*I_ERI_Dxy_Pz_Fy2z_Pz;
  abcd[1040] = I_ERI_Fx2z_Pz_Fy2z_Pz+ABZ*I_ERI_Dxz_Pz_Fy2z_Pz;
  abcd[1041] = I_ERI_F2yz_Pz_Fy2z_Pz+ABZ*I_ERI_D2y_Pz_Fy2z_Pz;
  abcd[1042] = I_ERI_Fy2z_Pz_Fy2z_Pz+ABZ*I_ERI_Dyz_Pz_Fy2z_Pz;
  abcd[1043] = I_ERI_F3z_Pz_Fy2z_Pz+ABZ*I_ERI_D2z_Pz_Fy2z_Pz;
  abcd[1044] = I_ERI_F3x_Px_F3z_Pz+ABX*I_ERI_D2x_Px_F3z_Pz;
  abcd[1045] = I_ERI_F2xy_Px_F3z_Pz+ABX*I_ERI_Dxy_Px_F3z_Pz;
  abcd[1046] = I_ERI_F2xz_Px_F3z_Pz+ABX*I_ERI_Dxz_Px_F3z_Pz;
  abcd[1047] = I_ERI_Fx2y_Px_F3z_Pz+ABX*I_ERI_D2y_Px_F3z_Pz;
  abcd[1048] = I_ERI_Fxyz_Px_F3z_Pz+ABX*I_ERI_Dyz_Px_F3z_Pz;
  abcd[1049] = I_ERI_Fx2z_Px_F3z_Pz+ABX*I_ERI_D2z_Px_F3z_Pz;
  abcd[1050] = I_ERI_F2xy_Px_F3z_Pz+ABY*I_ERI_D2x_Px_F3z_Pz;
  abcd[1051] = I_ERI_Fx2y_Px_F3z_Pz+ABY*I_ERI_Dxy_Px_F3z_Pz;
  abcd[1052] = I_ERI_Fxyz_Px_F3z_Pz+ABY*I_ERI_Dxz_Px_F3z_Pz;
  abcd[1053] = I_ERI_F3y_Px_F3z_Pz+ABY*I_ERI_D2y_Px_F3z_Pz;
  abcd[1054] = I_ERI_F2yz_Px_F3z_Pz+ABY*I_ERI_Dyz_Px_F3z_Pz;
  abcd[1055] = I_ERI_Fy2z_Px_F3z_Pz+ABY*I_ERI_D2z_Px_F3z_Pz;
  abcd[1056] = I_ERI_F2xz_Px_F3z_Pz+ABZ*I_ERI_D2x_Px_F3z_Pz;
  abcd[1057] = I_ERI_Fxyz_Px_F3z_Pz+ABZ*I_ERI_Dxy_Px_F3z_Pz;
  abcd[1058] = I_ERI_Fx2z_Px_F3z_Pz+ABZ*I_ERI_Dxz_Px_F3z_Pz;
  abcd[1059] = I_ERI_F2yz_Px_F3z_Pz+ABZ*I_ERI_D2y_Px_F3z_Pz;
  abcd[1060] = I_ERI_Fy2z_Px_F3z_Pz+ABZ*I_ERI_Dyz_Px_F3z_Pz;
  abcd[1061] = I_ERI_F3z_Px_F3z_Pz+ABZ*I_ERI_D2z_Px_F3z_Pz;
  abcd[1062] = I_ERI_F2xy_Py_F3z_Pz+ABY*I_ERI_D2x_Py_F3z_Pz;
  abcd[1063] = I_ERI_Fx2y_Py_F3z_Pz+ABY*I_ERI_Dxy_Py_F3z_Pz;
  abcd[1064] = I_ERI_Fxyz_Py_F3z_Pz+ABY*I_ERI_Dxz_Py_F3z_Pz;
  abcd[1065] = I_ERI_F3y_Py_F3z_Pz+ABY*I_ERI_D2y_Py_F3z_Pz;
  abcd[1066] = I_ERI_F2yz_Py_F3z_Pz+ABY*I_ERI_Dyz_Py_F3z_Pz;
  abcd[1067] = I_ERI_Fy2z_Py_F3z_Pz+ABY*I_ERI_D2z_Py_F3z_Pz;
  abcd[1068] = I_ERI_F2xz_Py_F3z_Pz+ABZ*I_ERI_D2x_Py_F3z_Pz;
  abcd[1069] = I_ERI_Fxyz_Py_F3z_Pz+ABZ*I_ERI_Dxy_Py_F3z_Pz;
  abcd[1070] = I_ERI_Fx2z_Py_F3z_Pz+ABZ*I_ERI_Dxz_Py_F3z_Pz;
  abcd[1071] = I_ERI_F2yz_Py_F3z_Pz+ABZ*I_ERI_D2y_Py_F3z_Pz;
  abcd[1072] = I_ERI_Fy2z_Py_F3z_Pz+ABZ*I_ERI_Dyz_Py_F3z_Pz;
  abcd[1073] = I_ERI_F3z_Py_F3z_Pz+ABZ*I_ERI_D2z_Py_F3z_Pz;
  abcd[1074] = I_ERI_F2xz_Pz_F3z_Pz+ABZ*I_ERI_D2x_Pz_F3z_Pz;
  abcd[1075] = I_ERI_Fxyz_Pz_F3z_Pz+ABZ*I_ERI_Dxy_Pz_F3z_Pz;
  abcd[1076] = I_ERI_Fx2z_Pz_F3z_Pz+ABZ*I_ERI_Dxz_Pz_F3z_Pz;
  abcd[1077] = I_ERI_F2yz_Pz_F3z_Pz+ABZ*I_ERI_D2y_Pz_F3z_Pz;
  abcd[1078] = I_ERI_Fy2z_Pz_F3z_Pz+ABZ*I_ERI_Dyz_Pz_F3z_Pz;
  abcd[1079] = I_ERI_F3z_Pz_F3z_Pz+ABZ*I_ERI_D2z_Pz_F3z_Pz;
}
