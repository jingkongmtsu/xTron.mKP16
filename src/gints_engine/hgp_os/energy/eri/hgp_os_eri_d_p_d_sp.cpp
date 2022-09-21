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

void hgp_os_eri_d_p_d_sp(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_F3x_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C2001002 = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C2001002 = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C2001002 = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C2001002 = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C2001002 = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C2001002 = 0.0E0;
  Double I_ERI_F3x_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_F3x_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3y_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_F3z_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_F3x_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_F2xy_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_F2xz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_Fx2y_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_Fxyz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_Fx2z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_F3y_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_F2yz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_Fy2z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_F3z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_D2x_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_Dxy_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_Dxz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_D2y_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_Dyz_S_C1002001002 = 0.0E0;
  Double I_ERI_D2x_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxy_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dxz_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2y_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_Dyz_S_D2z_S_C1002001002 = 0.0E0;
  Double I_ERI_D2z_S_D2z_S_C1002001002 = 0.0E0;

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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M5
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M5
       ************************************************************/
      Double I_ERI_S_S_D2x_S_M4_vrr = QCX*I_ERI_S_S_Px_S_M4_vrr+WQX*I_ERI_S_S_Px_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_Dxy_S_M4_vrr = QCY*I_ERI_S_S_Px_S_M4_vrr+WQY*I_ERI_S_S_Px_S_M5_vrr;
      Double I_ERI_S_S_D2y_S_M4_vrr = QCY*I_ERI_S_S_Py_S_M4_vrr+WQY*I_ERI_S_S_Py_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;
      Double I_ERI_S_S_D2z_S_M4_vrr = QCZ*I_ERI_S_S_Pz_S_M4_vrr+WQZ*I_ERI_S_S_Pz_S_M5_vrr+oned2e*I_ERI_S_S_S_S_M4_vrr-rhod2esq*I_ERI_S_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_D_S_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M2
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_Px_S_M1_vrr = PAX*I_ERI_Px_S_Px_S_M1_vrr+WPX*I_ERI_Px_S_Px_S_M2_vrr+oned2z*I_ERI_S_S_Px_S_M1_vrr-rhod2zsq*I_ERI_S_S_Px_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_Px_S_M1_vrr = PAY*I_ERI_Px_S_Px_S_M1_vrr+WPY*I_ERI_Px_S_Px_S_M2_vrr;
      Double I_ERI_D2y_S_Px_S_M1_vrr = PAY*I_ERI_Py_S_Px_S_M1_vrr+WPY*I_ERI_Py_S_Px_S_M2_vrr+oned2z*I_ERI_S_S_Px_S_M1_vrr-rhod2zsq*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_D2z_S_Px_S_M1_vrr = PAZ*I_ERI_Pz_S_Px_S_M1_vrr+WPZ*I_ERI_Pz_S_Px_S_M2_vrr+oned2z*I_ERI_S_S_Px_S_M1_vrr-rhod2zsq*I_ERI_S_S_Px_S_M2_vrr;
      Double I_ERI_D2x_S_Py_S_M1_vrr = PAX*I_ERI_Px_S_Py_S_M1_vrr+WPX*I_ERI_Px_S_Py_S_M2_vrr+oned2z*I_ERI_S_S_Py_S_M1_vrr-rhod2zsq*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_Dxy_S_Py_S_M1_vrr = PAY*I_ERI_Px_S_Py_S_M1_vrr+WPY*I_ERI_Px_S_Py_S_M2_vrr+oned2k*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_Py_S_M1_vrr = PAY*I_ERI_Py_S_Py_S_M1_vrr+WPY*I_ERI_Py_S_Py_S_M2_vrr+oned2z*I_ERI_S_S_Py_S_M1_vrr-rhod2zsq*I_ERI_S_S_Py_S_M2_vrr+oned2k*I_ERI_Py_S_S_S_M2_vrr;
      Double I_ERI_D2z_S_Py_S_M1_vrr = PAZ*I_ERI_Pz_S_Py_S_M1_vrr+WPZ*I_ERI_Pz_S_Py_S_M2_vrr+oned2z*I_ERI_S_S_Py_S_M1_vrr-rhod2zsq*I_ERI_S_S_Py_S_M2_vrr;
      Double I_ERI_D2x_S_Pz_S_M1_vrr = PAX*I_ERI_Px_S_Pz_S_M1_vrr+WPX*I_ERI_Px_S_Pz_S_M2_vrr+oned2z*I_ERI_S_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_Dxy_S_Pz_S_M1_vrr = PAY*I_ERI_Px_S_Pz_S_M1_vrr+WPY*I_ERI_Px_S_Pz_S_M2_vrr;
      Double I_ERI_D2y_S_Pz_S_M1_vrr = PAY*I_ERI_Py_S_Pz_S_M1_vrr+WPY*I_ERI_Py_S_Pz_S_M2_vrr+oned2z*I_ERI_S_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M2_vrr;
      Double I_ERI_D2z_S_Pz_S_M1_vrr = PAZ*I_ERI_Pz_S_Pz_S_M1_vrr+WPZ*I_ERI_Pz_S_Pz_S_M2_vrr+oned2z*I_ERI_S_S_Pz_S_M1_vrr-rhod2zsq*I_ERI_S_S_Pz_S_M2_vrr+oned2k*I_ERI_Pz_S_S_S_M2_vrr;

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
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_P_S
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_S_S_D2x_S_vrr = QCX*I_ERI_S_S_Px_S_vrr+WQX*I_ERI_S_S_Px_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Dxy_S_vrr = QCY*I_ERI_S_S_Px_S_vrr+WQY*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_Dxz_S_vrr = QCZ*I_ERI_S_S_Px_S_vrr+WQZ*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_S_S_D2y_S_vrr = QCY*I_ERI_S_S_Py_S_vrr+WQY*I_ERI_S_S_Py_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_S_S_Dyz_S_vrr = QCZ*I_ERI_S_S_Py_S_vrr+WQZ*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_S_S_D2z_S_vrr = QCZ*I_ERI_S_S_Pz_S_vrr+WQZ*I_ERI_S_S_Pz_S_M1_vrr+oned2e*I_ERI_S_S_S_S_vrr-rhod2esq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_P_S_M1
       ************************************************************/
      Double I_ERI_Px_S_D2x_S_vrr = PAX*I_ERI_S_S_D2x_S_vrr+WPX*I_ERI_S_S_D2x_S_M1_vrr+2*oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Py_S_D2x_S_vrr = PAY*I_ERI_S_S_D2x_S_vrr+WPY*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Pz_S_D2x_S_vrr = PAZ*I_ERI_S_S_D2x_S_vrr+WPZ*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Px_S_Dxy_S_vrr = PAX*I_ERI_S_S_Dxy_S_vrr+WPX*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Py_S_Dxy_S_vrr = PAY*I_ERI_S_S_Dxy_S_vrr+WPY*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Pz_S_Dxy_S_vrr = PAZ*I_ERI_S_S_Dxy_S_vrr+WPZ*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_Px_S_Dxz_S_vrr = PAX*I_ERI_S_S_Dxz_S_vrr+WPX*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Py_S_Dxz_S_vrr = PAY*I_ERI_S_S_Dxz_S_vrr+WPY*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Pz_S_Dxz_S_vrr = PAZ*I_ERI_S_S_Dxz_S_vrr+WPZ*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_S_S_Px_S_M1_vrr;
      Double I_ERI_Px_S_D2y_S_vrr = PAX*I_ERI_S_S_D2y_S_vrr+WPX*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Py_S_D2y_S_vrr = PAY*I_ERI_S_S_D2y_S_vrr+WPY*I_ERI_S_S_D2y_S_M1_vrr+2*oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Pz_S_D2y_S_vrr = PAZ*I_ERI_S_S_D2y_S_vrr+WPZ*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Px_S_Dyz_S_vrr = PAX*I_ERI_S_S_Dyz_S_vrr+WPX*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Py_S_Dyz_S_vrr = PAY*I_ERI_S_S_Dyz_S_vrr+WPY*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_S_S_Pz_S_M1_vrr;
      Double I_ERI_Pz_S_Dyz_S_vrr = PAZ*I_ERI_S_S_Dyz_S_vrr+WPZ*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_S_S_Py_S_M1_vrr;
      Double I_ERI_Px_S_D2z_S_vrr = PAX*I_ERI_S_S_D2z_S_vrr+WPX*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Py_S_D2z_S_vrr = PAY*I_ERI_S_S_D2z_S_vrr+WPY*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Pz_S_D2z_S_vrr = PAZ*I_ERI_S_S_D2z_S_vrr+WPZ*I_ERI_S_S_D2z_S_M1_vrr+2*oned2k*I_ERI_S_S_Pz_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_D_S_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_D_S
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_D_S
       * RHS shell quartet name: SQ_ERI_S_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_P_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_D2x_S_vrr = PAX*I_ERI_Px_S_D2x_S_vrr+WPX*I_ERI_Px_S_D2x_S_M1_vrr+oned2z*I_ERI_S_S_D2x_S_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M1_vrr+2*oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Dxy_S_D2x_S_vrr = PAY*I_ERI_Px_S_D2x_S_vrr+WPY*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_Dxz_S_D2x_S_vrr = PAZ*I_ERI_Px_S_D2x_S_vrr+WPZ*I_ERI_Px_S_D2x_S_M1_vrr;
      Double I_ERI_D2y_S_D2x_S_vrr = PAY*I_ERI_Py_S_D2x_S_vrr+WPY*I_ERI_Py_S_D2x_S_M1_vrr+oned2z*I_ERI_S_S_D2x_S_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_Dyz_S_D2x_S_vrr = PAZ*I_ERI_Py_S_D2x_S_vrr+WPZ*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_D2z_S_D2x_S_vrr = PAZ*I_ERI_Pz_S_D2x_S_vrr+WPZ*I_ERI_Pz_S_D2x_S_M1_vrr+oned2z*I_ERI_S_S_D2x_S_vrr-rhod2zsq*I_ERI_S_S_D2x_S_M1_vrr;
      Double I_ERI_D2x_S_Dxy_S_vrr = PAX*I_ERI_Px_S_Dxy_S_vrr+WPX*I_ERI_Px_S_Dxy_S_M1_vrr+oned2z*I_ERI_S_S_Dxy_S_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_Dxy_S_Dxy_S_vrr = PAY*I_ERI_Px_S_Dxy_S_vrr+WPY*I_ERI_Px_S_Dxy_S_M1_vrr+oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_Dxz_S_Dxy_S_vrr = PAZ*I_ERI_Px_S_Dxy_S_vrr+WPZ*I_ERI_Px_S_Dxy_S_M1_vrr;
      Double I_ERI_D2y_S_Dxy_S_vrr = PAY*I_ERI_Py_S_Dxy_S_vrr+WPY*I_ERI_Py_S_Dxy_S_M1_vrr+oned2z*I_ERI_S_S_Dxy_S_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M1_vrr+oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_Dyz_S_Dxy_S_vrr = PAZ*I_ERI_Py_S_Dxy_S_vrr+WPZ*I_ERI_Py_S_Dxy_S_M1_vrr;
      Double I_ERI_D2z_S_Dxy_S_vrr = PAZ*I_ERI_Pz_S_Dxy_S_vrr+WPZ*I_ERI_Pz_S_Dxy_S_M1_vrr+oned2z*I_ERI_S_S_Dxy_S_vrr-rhod2zsq*I_ERI_S_S_Dxy_S_M1_vrr;
      Double I_ERI_D2x_S_Dxz_S_vrr = PAX*I_ERI_Px_S_Dxz_S_vrr+WPX*I_ERI_Px_S_Dxz_S_M1_vrr+oned2z*I_ERI_S_S_Dxz_S_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_Dxy_S_Dxz_S_vrr = PAY*I_ERI_Px_S_Dxz_S_vrr+WPY*I_ERI_Px_S_Dxz_S_M1_vrr;
      Double I_ERI_Dxz_S_Dxz_S_vrr = PAZ*I_ERI_Px_S_Dxz_S_vrr+WPZ*I_ERI_Px_S_Dxz_S_M1_vrr+oned2k*I_ERI_Px_S_Px_S_M1_vrr;
      Double I_ERI_D2y_S_Dxz_S_vrr = PAY*I_ERI_Py_S_Dxz_S_vrr+WPY*I_ERI_Py_S_Dxz_S_M1_vrr+oned2z*I_ERI_S_S_Dxz_S_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M1_vrr;
      Double I_ERI_Dyz_S_Dxz_S_vrr = PAZ*I_ERI_Py_S_Dxz_S_vrr+WPZ*I_ERI_Py_S_Dxz_S_M1_vrr+oned2k*I_ERI_Py_S_Px_S_M1_vrr;
      Double I_ERI_D2z_S_Dxz_S_vrr = PAZ*I_ERI_Pz_S_Dxz_S_vrr+WPZ*I_ERI_Pz_S_Dxz_S_M1_vrr+oned2z*I_ERI_S_S_Dxz_S_vrr-rhod2zsq*I_ERI_S_S_Dxz_S_M1_vrr+oned2k*I_ERI_Pz_S_Px_S_M1_vrr;
      Double I_ERI_D2x_S_D2y_S_vrr = PAX*I_ERI_Px_S_D2y_S_vrr+WPX*I_ERI_Px_S_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_Dxy_S_D2y_S_vrr = PAY*I_ERI_Px_S_D2y_S_vrr+WPY*I_ERI_Px_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_Dxz_S_D2y_S_vrr = PAZ*I_ERI_Px_S_D2y_S_vrr+WPZ*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_D2y_S_D2y_S_vrr = PAY*I_ERI_Py_S_D2y_S_vrr+WPY*I_ERI_Py_S_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr+2*oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_Dyz_S_D2y_S_vrr = PAZ*I_ERI_Py_S_D2y_S_vrr+WPZ*I_ERI_Py_S_D2y_S_M1_vrr;
      Double I_ERI_D2z_S_D2y_S_vrr = PAZ*I_ERI_Pz_S_D2y_S_vrr+WPZ*I_ERI_Pz_S_D2y_S_M1_vrr+oned2z*I_ERI_S_S_D2y_S_vrr-rhod2zsq*I_ERI_S_S_D2y_S_M1_vrr;
      Double I_ERI_D2x_S_Dyz_S_vrr = PAX*I_ERI_Px_S_Dyz_S_vrr+WPX*I_ERI_Px_S_Dyz_S_M1_vrr+oned2z*I_ERI_S_S_Dyz_S_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M1_vrr;
      Double I_ERI_Dxy_S_Dyz_S_vrr = PAY*I_ERI_Px_S_Dyz_S_vrr+WPY*I_ERI_Px_S_Dyz_S_M1_vrr+oned2k*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_Dxz_S_Dyz_S_vrr = PAZ*I_ERI_Px_S_Dyz_S_vrr+WPZ*I_ERI_Px_S_Dyz_S_M1_vrr+oned2k*I_ERI_Px_S_Py_S_M1_vrr;
      Double I_ERI_D2y_S_Dyz_S_vrr = PAY*I_ERI_Py_S_Dyz_S_vrr+WPY*I_ERI_Py_S_Dyz_S_M1_vrr+oned2z*I_ERI_S_S_Dyz_S_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_Dyz_S_Dyz_S_vrr = PAZ*I_ERI_Py_S_Dyz_S_vrr+WPZ*I_ERI_Py_S_Dyz_S_M1_vrr+oned2k*I_ERI_Py_S_Py_S_M1_vrr;
      Double I_ERI_D2z_S_Dyz_S_vrr = PAZ*I_ERI_Pz_S_Dyz_S_vrr+WPZ*I_ERI_Pz_S_Dyz_S_M1_vrr+oned2z*I_ERI_S_S_Dyz_S_vrr-rhod2zsq*I_ERI_S_S_Dyz_S_M1_vrr+oned2k*I_ERI_Pz_S_Py_S_M1_vrr;
      Double I_ERI_D2x_S_D2z_S_vrr = PAX*I_ERI_Px_S_D2z_S_vrr+WPX*I_ERI_Px_S_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Dxy_S_D2z_S_vrr = PAY*I_ERI_Px_S_D2z_S_vrr+WPY*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_Dxz_S_D2z_S_vrr = PAZ*I_ERI_Px_S_D2z_S_vrr+WPZ*I_ERI_Px_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Px_S_Pz_S_M1_vrr;
      Double I_ERI_D2y_S_D2z_S_vrr = PAY*I_ERI_Py_S_D2z_S_vrr+WPY*I_ERI_Py_S_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr;
      Double I_ERI_Dyz_S_D2z_S_vrr = PAZ*I_ERI_Py_S_D2z_S_vrr+WPZ*I_ERI_Py_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Py_S_Pz_S_M1_vrr;
      Double I_ERI_D2z_S_D2z_S_vrr = PAZ*I_ERI_Pz_S_D2z_S_vrr+WPZ*I_ERI_Pz_S_D2z_S_M1_vrr+oned2z*I_ERI_S_S_D2z_S_vrr-rhod2zsq*I_ERI_S_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Pz_S_Pz_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_F_S_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_D_S
       * RHS shell quartet name: SQ_ERI_D_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_D_S
       * RHS shell quartet name: SQ_ERI_P_S_D_S_M1
       * RHS shell quartet name: SQ_ERI_D_S_P_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_D2x_S_vrr = PAX*I_ERI_D2x_S_D2x_S_vrr+WPX*I_ERI_D2x_S_D2x_S_M1_vrr+2*oned2z*I_ERI_Px_S_D2x_S_vrr-2*rhod2zsq*I_ERI_Px_S_D2x_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_F2xy_S_D2x_S_vrr = PAY*I_ERI_D2x_S_D2x_S_vrr+WPY*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_F2xz_S_D2x_S_vrr = PAZ*I_ERI_D2x_S_D2x_S_vrr+WPZ*I_ERI_D2x_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2y_S_D2x_S_vrr = PAX*I_ERI_D2y_S_D2x_S_vrr+WPX*I_ERI_D2y_S_D2x_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Fxyz_S_D2x_S_vrr = PAZ*I_ERI_Dxy_S_D2x_S_vrr+WPZ*I_ERI_Dxy_S_D2x_S_M1_vrr;
      Double I_ERI_Fx2z_S_D2x_S_vrr = PAX*I_ERI_D2z_S_D2x_S_vrr+WPX*I_ERI_D2z_S_D2x_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3y_S_D2x_S_vrr = PAY*I_ERI_D2y_S_D2x_S_vrr+WPY*I_ERI_D2y_S_D2x_S_M1_vrr+2*oned2z*I_ERI_Py_S_D2x_S_vrr-2*rhod2zsq*I_ERI_Py_S_D2x_S_M1_vrr;
      Double I_ERI_F2yz_S_D2x_S_vrr = PAZ*I_ERI_D2y_S_D2x_S_vrr+WPZ*I_ERI_D2y_S_D2x_S_M1_vrr;
      Double I_ERI_Fy2z_S_D2x_S_vrr = PAY*I_ERI_D2z_S_D2x_S_vrr+WPY*I_ERI_D2z_S_D2x_S_M1_vrr;
      Double I_ERI_F3z_S_D2x_S_vrr = PAZ*I_ERI_D2z_S_D2x_S_vrr+WPZ*I_ERI_D2z_S_D2x_S_M1_vrr+2*oned2z*I_ERI_Pz_S_D2x_S_vrr-2*rhod2zsq*I_ERI_Pz_S_D2x_S_M1_vrr;
      Double I_ERI_F3x_S_Dxy_S_vrr = PAX*I_ERI_D2x_S_Dxy_S_vrr+WPX*I_ERI_D2x_S_Dxy_S_M1_vrr+2*oned2z*I_ERI_Px_S_Dxy_S_vrr-2*rhod2zsq*I_ERI_Px_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_F2xy_S_Dxy_S_vrr = PAY*I_ERI_D2x_S_Dxy_S_vrr+WPY*I_ERI_D2x_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_F2xz_S_Dxy_S_vrr = PAZ*I_ERI_D2x_S_Dxy_S_vrr+WPZ*I_ERI_D2x_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2y_S_Dxy_S_vrr = PAX*I_ERI_D2y_S_Dxy_S_vrr+WPX*I_ERI_D2y_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Fxyz_S_Dxy_S_vrr = PAZ*I_ERI_Dxy_S_Dxy_S_vrr+WPZ*I_ERI_Dxy_S_Dxy_S_M1_vrr;
      Double I_ERI_Fx2z_S_Dxy_S_vrr = PAX*I_ERI_D2z_S_Dxy_S_vrr+WPX*I_ERI_D2z_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3y_S_Dxy_S_vrr = PAY*I_ERI_D2y_S_Dxy_S_vrr+WPY*I_ERI_D2y_S_Dxy_S_M1_vrr+2*oned2z*I_ERI_Py_S_Dxy_S_vrr-2*rhod2zsq*I_ERI_Py_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_F2yz_S_Dxy_S_vrr = PAZ*I_ERI_D2y_S_Dxy_S_vrr+WPZ*I_ERI_D2y_S_Dxy_S_M1_vrr;
      Double I_ERI_Fy2z_S_Dxy_S_vrr = PAY*I_ERI_D2z_S_Dxy_S_vrr+WPY*I_ERI_D2z_S_Dxy_S_M1_vrr+oned2k*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3z_S_Dxy_S_vrr = PAZ*I_ERI_D2z_S_Dxy_S_vrr+WPZ*I_ERI_D2z_S_Dxy_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Dxy_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxy_S_M1_vrr;
      Double I_ERI_F3x_S_Dxz_S_vrr = PAX*I_ERI_D2x_S_Dxz_S_vrr+WPX*I_ERI_D2x_S_Dxz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Dxz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_F2xy_S_Dxz_S_vrr = PAY*I_ERI_D2x_S_Dxz_S_vrr+WPY*I_ERI_D2x_S_Dxz_S_M1_vrr;
      Double I_ERI_F2xz_S_Dxz_S_vrr = PAZ*I_ERI_D2x_S_Dxz_S_vrr+WPZ*I_ERI_D2x_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2x_S_Px_S_M1_vrr;
      Double I_ERI_Fx2y_S_Dxz_S_vrr = PAX*I_ERI_D2y_S_Dxz_S_vrr+WPX*I_ERI_D2y_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_Fxyz_S_Dxz_S_vrr = PAZ*I_ERI_Dxy_S_Dxz_S_vrr+WPZ*I_ERI_Dxy_S_Dxz_S_M1_vrr+oned2k*I_ERI_Dxy_S_Px_S_M1_vrr;
      Double I_ERI_Fx2z_S_Dxz_S_vrr = PAX*I_ERI_D2z_S_Dxz_S_vrr+WPX*I_ERI_D2z_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2z_S_Pz_S_M1_vrr;
      Double I_ERI_F3y_S_Dxz_S_vrr = PAY*I_ERI_D2y_S_Dxz_S_vrr+WPY*I_ERI_D2y_S_Dxz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Dxz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Dxz_S_M1_vrr;
      Double I_ERI_F2yz_S_Dxz_S_vrr = PAZ*I_ERI_D2y_S_Dxz_S_vrr+WPZ*I_ERI_D2y_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2y_S_Px_S_M1_vrr;
      Double I_ERI_Fy2z_S_Dxz_S_vrr = PAY*I_ERI_D2z_S_Dxz_S_vrr+WPY*I_ERI_D2z_S_Dxz_S_M1_vrr;
      Double I_ERI_F3z_S_Dxz_S_vrr = PAZ*I_ERI_D2z_S_Dxz_S_vrr+WPZ*I_ERI_D2z_S_Dxz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Dxz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Dxz_S_M1_vrr+oned2k*I_ERI_D2z_S_Px_S_M1_vrr;
      Double I_ERI_F3x_S_D2y_S_vrr = PAX*I_ERI_D2x_S_D2y_S_vrr+WPX*I_ERI_D2x_S_D2y_S_M1_vrr+2*oned2z*I_ERI_Px_S_D2y_S_vrr-2*rhod2zsq*I_ERI_Px_S_D2y_S_M1_vrr;
      Double I_ERI_F2xy_S_D2y_S_vrr = PAY*I_ERI_D2x_S_D2y_S_vrr+WPY*I_ERI_D2x_S_D2y_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_F2xz_S_D2y_S_vrr = PAZ*I_ERI_D2x_S_D2y_S_vrr+WPZ*I_ERI_D2x_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2y_S_D2y_S_vrr = PAX*I_ERI_D2y_S_D2y_S_vrr+WPX*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fxyz_S_D2y_S_vrr = PAZ*I_ERI_Dxy_S_D2y_S_vrr+WPZ*I_ERI_Dxy_S_D2y_S_M1_vrr;
      Double I_ERI_Fx2z_S_D2y_S_vrr = PAX*I_ERI_D2z_S_D2y_S_vrr+WPX*I_ERI_D2z_S_D2y_S_M1_vrr;
      Double I_ERI_F3y_S_D2y_S_vrr = PAY*I_ERI_D2y_S_D2y_S_vrr+WPY*I_ERI_D2y_S_D2y_S_M1_vrr+2*oned2z*I_ERI_Py_S_D2y_S_vrr-2*rhod2zsq*I_ERI_Py_S_D2y_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_F2yz_S_D2y_S_vrr = PAZ*I_ERI_D2y_S_D2y_S_vrr+WPZ*I_ERI_D2y_S_D2y_S_M1_vrr;
      Double I_ERI_Fy2z_S_D2y_S_vrr = PAY*I_ERI_D2z_S_D2y_S_vrr+WPY*I_ERI_D2z_S_D2y_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3z_S_D2y_S_vrr = PAZ*I_ERI_D2z_S_D2y_S_vrr+WPZ*I_ERI_D2z_S_D2y_S_M1_vrr+2*oned2z*I_ERI_Pz_S_D2y_S_vrr-2*rhod2zsq*I_ERI_Pz_S_D2y_S_M1_vrr;
      Double I_ERI_F3x_S_Dyz_S_vrr = PAX*I_ERI_D2x_S_Dyz_S_vrr+WPX*I_ERI_D2x_S_Dyz_S_M1_vrr+2*oned2z*I_ERI_Px_S_Dyz_S_vrr-2*rhod2zsq*I_ERI_Px_S_Dyz_S_M1_vrr;
      Double I_ERI_F2xy_S_Dyz_S_vrr = PAY*I_ERI_D2x_S_Dyz_S_vrr+WPY*I_ERI_D2x_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_F2xz_S_Dyz_S_vrr = PAZ*I_ERI_D2x_S_Dyz_S_vrr+WPZ*I_ERI_D2x_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2x_S_Py_S_M1_vrr;
      Double I_ERI_Fx2y_S_Dyz_S_vrr = PAX*I_ERI_D2y_S_Dyz_S_vrr+WPX*I_ERI_D2y_S_Dyz_S_M1_vrr;
      Double I_ERI_Fxyz_S_Dyz_S_vrr = PAZ*I_ERI_Dxy_S_Dyz_S_vrr+WPZ*I_ERI_Dxy_S_Dyz_S_M1_vrr+oned2k*I_ERI_Dxy_S_Py_S_M1_vrr;
      Double I_ERI_Fx2z_S_Dyz_S_vrr = PAX*I_ERI_D2z_S_Dyz_S_vrr+WPX*I_ERI_D2z_S_Dyz_S_M1_vrr;
      Double I_ERI_F3y_S_Dyz_S_vrr = PAY*I_ERI_D2y_S_Dyz_S_vrr+WPY*I_ERI_D2y_S_Dyz_S_M1_vrr+2*oned2z*I_ERI_Py_S_Dyz_S_vrr-2*rhod2zsq*I_ERI_Py_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_F2yz_S_Dyz_S_vrr = PAZ*I_ERI_D2y_S_Dyz_S_vrr+WPZ*I_ERI_D2y_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2y_S_Py_S_M1_vrr;
      Double I_ERI_Fy2z_S_Dyz_S_vrr = PAY*I_ERI_D2z_S_Dyz_S_vrr+WPY*I_ERI_D2z_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Pz_S_M1_vrr;
      Double I_ERI_F3z_S_Dyz_S_vrr = PAZ*I_ERI_D2z_S_Dyz_S_vrr+WPZ*I_ERI_D2z_S_Dyz_S_M1_vrr+2*oned2z*I_ERI_Pz_S_Dyz_S_vrr-2*rhod2zsq*I_ERI_Pz_S_Dyz_S_M1_vrr+oned2k*I_ERI_D2z_S_Py_S_M1_vrr;
      Double I_ERI_F3x_S_D2z_S_vrr = PAX*I_ERI_D2x_S_D2z_S_vrr+WPX*I_ERI_D2x_S_D2z_S_M1_vrr+2*oned2z*I_ERI_Px_S_D2z_S_vrr-2*rhod2zsq*I_ERI_Px_S_D2z_S_M1_vrr;
      Double I_ERI_F2xy_S_D2z_S_vrr = PAY*I_ERI_D2x_S_D2z_S_vrr+WPY*I_ERI_D2x_S_D2z_S_M1_vrr;
      Double I_ERI_F2xz_S_D2z_S_vrr = PAZ*I_ERI_D2x_S_D2z_S_vrr+WPZ*I_ERI_D2x_S_D2z_S_M1_vrr+2*oned2k*I_ERI_D2x_S_Pz_S_M1_vrr;
      Double I_ERI_Fx2y_S_D2z_S_vrr = PAX*I_ERI_D2y_S_D2z_S_vrr+WPX*I_ERI_D2y_S_D2z_S_M1_vrr;
      Double I_ERI_Fxyz_S_D2z_S_vrr = PAZ*I_ERI_Dxy_S_D2z_S_vrr+WPZ*I_ERI_Dxy_S_D2z_S_M1_vrr+2*oned2k*I_ERI_Dxy_S_Pz_S_M1_vrr;
      Double I_ERI_Fx2z_S_D2z_S_vrr = PAX*I_ERI_D2z_S_D2z_S_vrr+WPX*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3y_S_D2z_S_vrr = PAY*I_ERI_D2y_S_D2z_S_vrr+WPY*I_ERI_D2y_S_D2z_S_M1_vrr+2*oned2z*I_ERI_Py_S_D2z_S_vrr-2*rhod2zsq*I_ERI_Py_S_D2z_S_M1_vrr;
      Double I_ERI_F2yz_S_D2z_S_vrr = PAZ*I_ERI_D2y_S_D2z_S_vrr+WPZ*I_ERI_D2y_S_D2z_S_M1_vrr+2*oned2k*I_ERI_D2y_S_Pz_S_M1_vrr;
      Double I_ERI_Fy2z_S_D2z_S_vrr = PAY*I_ERI_D2z_S_D2z_S_vrr+WPY*I_ERI_D2z_S_D2z_S_M1_vrr;
      Double I_ERI_F3z_S_D2z_S_vrr = PAZ*I_ERI_D2z_S_D2z_S_vrr+WPZ*I_ERI_D2z_S_D2z_S_M1_vrr+2*oned2z*I_ERI_Pz_S_D2z_S_vrr-2*rhod2zsq*I_ERI_Pz_S_D2z_S_M1_vrr+2*oned2k*I_ERI_D2z_S_Pz_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_F_S_D_S_C2001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C2001002_coefs = ic2*jc2;
      I_ERI_F3x_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C2001002 += SQ_ERI_F_S_D_S_C2001002_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C2001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C2001002_coefs = ic2*jc2;
      I_ERI_D2x_S_D2x_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C2001002 += SQ_ERI_D_S_D_S_C2001002_coefs*I_ERI_D2z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_F_S_C1002001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_F_S_C1002001002_coefs = ic2*jc2_1;
      I_ERI_F3x_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_F3x_S_vrr;
      I_ERI_F2xy_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_F3x_S_vrr;
      I_ERI_F2xz_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_F3x_S_vrr;
      I_ERI_Fx2y_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_F3x_S_vrr;
      I_ERI_Fxyz_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_F3x_S_vrr;
      I_ERI_Fx2z_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_F3x_S_vrr;
      I_ERI_F3y_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_F3x_S_vrr;
      I_ERI_F2yz_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_F3x_S_vrr;
      I_ERI_Fy2z_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_F3x_S_vrr;
      I_ERI_F3z_S_F3x_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_F3x_S_vrr;
      I_ERI_F3x_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_F2xy_S_vrr;
      I_ERI_F2xy_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_F2xy_S_vrr;
      I_ERI_F2xz_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_F2xy_S_vrr;
      I_ERI_Fx2y_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_F2xy_S_vrr;
      I_ERI_Fxyz_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_F2xy_S_vrr;
      I_ERI_Fx2z_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_F2xy_S_vrr;
      I_ERI_F3y_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_F2xy_S_vrr;
      I_ERI_F2yz_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_F2xy_S_vrr;
      I_ERI_Fy2z_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_F2xy_S_vrr;
      I_ERI_F3z_S_F2xy_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_F2xy_S_vrr;
      I_ERI_F3x_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_F2xz_S_vrr;
      I_ERI_F2xy_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_F2xz_S_vrr;
      I_ERI_F2xz_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_F2xz_S_vrr;
      I_ERI_Fx2y_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_F2xz_S_vrr;
      I_ERI_Fxyz_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_F2xz_S_vrr;
      I_ERI_Fx2z_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_F2xz_S_vrr;
      I_ERI_F3y_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_F2xz_S_vrr;
      I_ERI_F2yz_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_F2xz_S_vrr;
      I_ERI_Fy2z_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_F2xz_S_vrr;
      I_ERI_F3z_S_F2xz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_F2xz_S_vrr;
      I_ERI_F3x_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_Fx2y_S_vrr;
      I_ERI_F2xy_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_Fx2y_S_vrr;
      I_ERI_F2xz_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_Fx2y_S_vrr;
      I_ERI_Fx2y_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_Fx2y_S_vrr;
      I_ERI_Fxyz_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_Fx2y_S_vrr;
      I_ERI_Fx2z_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_Fx2y_S_vrr;
      I_ERI_F3y_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_Fx2y_S_vrr;
      I_ERI_F2yz_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_Fx2y_S_vrr;
      I_ERI_Fy2z_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_Fx2y_S_vrr;
      I_ERI_F3z_S_Fx2y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_Fx2y_S_vrr;
      I_ERI_F3x_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_Fxyz_S_vrr;
      I_ERI_F2xy_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_Fxyz_S_vrr;
      I_ERI_F2xz_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_Fxyz_S_vrr;
      I_ERI_Fx2y_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_Fxyz_S_vrr;
      I_ERI_Fxyz_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_Fxyz_S_vrr;
      I_ERI_Fx2z_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_Fxyz_S_vrr;
      I_ERI_F3y_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_Fxyz_S_vrr;
      I_ERI_F2yz_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_Fxyz_S_vrr;
      I_ERI_Fy2z_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_Fxyz_S_vrr;
      I_ERI_F3z_S_Fxyz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_Fxyz_S_vrr;
      I_ERI_F3x_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_Fx2z_S_vrr;
      I_ERI_F2xy_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_Fx2z_S_vrr;
      I_ERI_F2xz_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_Fx2z_S_vrr;
      I_ERI_Fx2y_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_Fx2z_S_vrr;
      I_ERI_Fxyz_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_Fx2z_S_vrr;
      I_ERI_Fx2z_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_Fx2z_S_vrr;
      I_ERI_F3y_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_Fx2z_S_vrr;
      I_ERI_F2yz_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_Fx2z_S_vrr;
      I_ERI_Fy2z_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_Fx2z_S_vrr;
      I_ERI_F3z_S_Fx2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_Fx2z_S_vrr;
      I_ERI_F3x_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_F3y_S_vrr;
      I_ERI_F2xy_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_F3y_S_vrr;
      I_ERI_F2xz_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_F3y_S_vrr;
      I_ERI_Fx2y_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_F3y_S_vrr;
      I_ERI_Fxyz_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_F3y_S_vrr;
      I_ERI_Fx2z_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_F3y_S_vrr;
      I_ERI_F3y_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_F3y_S_vrr;
      I_ERI_F2yz_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_F3y_S_vrr;
      I_ERI_Fy2z_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_F3y_S_vrr;
      I_ERI_F3z_S_F3y_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_F3y_S_vrr;
      I_ERI_F3x_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_F2yz_S_vrr;
      I_ERI_F2xy_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_F2yz_S_vrr;
      I_ERI_F2xz_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_F2yz_S_vrr;
      I_ERI_Fx2y_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_F2yz_S_vrr;
      I_ERI_Fxyz_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_F2yz_S_vrr;
      I_ERI_Fx2z_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_F2yz_S_vrr;
      I_ERI_F3y_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_F2yz_S_vrr;
      I_ERI_F2yz_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_F2yz_S_vrr;
      I_ERI_Fy2z_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_F2yz_S_vrr;
      I_ERI_F3z_S_F2yz_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_F2yz_S_vrr;
      I_ERI_F3x_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_Fy2z_S_vrr;
      I_ERI_F2xy_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_Fy2z_S_vrr;
      I_ERI_F2xz_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_Fy2z_S_vrr;
      I_ERI_Fx2y_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_Fy2z_S_vrr;
      I_ERI_Fxyz_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_Fy2z_S_vrr;
      I_ERI_Fx2z_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_Fy2z_S_vrr;
      I_ERI_F3y_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_Fy2z_S_vrr;
      I_ERI_F2yz_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_Fy2z_S_vrr;
      I_ERI_Fy2z_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_Fy2z_S_vrr;
      I_ERI_F3z_S_Fy2z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_Fy2z_S_vrr;
      I_ERI_F3x_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3x_S_F3z_S_vrr;
      I_ERI_F2xy_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xy_S_F3z_S_vrr;
      I_ERI_F2xz_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2xz_S_F3z_S_vrr;
      I_ERI_Fx2y_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2y_S_F3z_S_vrr;
      I_ERI_Fxyz_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fxyz_S_F3z_S_vrr;
      I_ERI_Fx2z_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fx2z_S_F3z_S_vrr;
      I_ERI_F3y_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3y_S_F3z_S_vrr;
      I_ERI_F2yz_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F2yz_S_F3z_S_vrr;
      I_ERI_Fy2z_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_Fy2z_S_F3z_S_vrr;
      I_ERI_F3z_S_F3z_S_C1002001002 += SQ_ERI_F_S_F_S_C1002001002_coefs*I_ERI_F3z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S_C1002001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_F_S_D_S_C1002001002_coefs = ic2*jc2_1;
      I_ERI_F3x_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S_C1002001002 += SQ_ERI_F_S_D_S_C1002001002_coefs*I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_F_S_C1002001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_F_S_C1002001002_coefs = ic2*jc2_1;
      I_ERI_D2x_S_F3x_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_F3x_S_vrr;
      I_ERI_Dxy_S_F3x_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_F3x_S_vrr;
      I_ERI_Dxz_S_F3x_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_F3x_S_vrr;
      I_ERI_D2y_S_F3x_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_F3x_S_vrr;
      I_ERI_Dyz_S_F3x_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_F3x_S_vrr;
      I_ERI_D2z_S_F3x_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_F3x_S_vrr;
      I_ERI_D2x_S_F2xy_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_F2xy_S_vrr;
      I_ERI_Dxy_S_F2xy_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_F2xy_S_vrr;
      I_ERI_Dxz_S_F2xy_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_F2xy_S_vrr;
      I_ERI_D2y_S_F2xy_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_F2xy_S_vrr;
      I_ERI_Dyz_S_F2xy_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_F2xy_S_vrr;
      I_ERI_D2z_S_F2xy_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_F2xy_S_vrr;
      I_ERI_D2x_S_F2xz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_F2xz_S_vrr;
      I_ERI_Dxy_S_F2xz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_F2xz_S_vrr;
      I_ERI_Dxz_S_F2xz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_F2xz_S_vrr;
      I_ERI_D2y_S_F2xz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_F2xz_S_vrr;
      I_ERI_Dyz_S_F2xz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_F2xz_S_vrr;
      I_ERI_D2z_S_F2xz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_F2xz_S_vrr;
      I_ERI_D2x_S_Fx2y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_Fx2y_S_vrr;
      I_ERI_Dxy_S_Fx2y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_Fx2y_S_vrr;
      I_ERI_Dxz_S_Fx2y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_Fx2y_S_vrr;
      I_ERI_D2y_S_Fx2y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_Fx2y_S_vrr;
      I_ERI_Dyz_S_Fx2y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_Fx2y_S_vrr;
      I_ERI_D2z_S_Fx2y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_Fx2y_S_vrr;
      I_ERI_D2x_S_Fxyz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_Fxyz_S_vrr;
      I_ERI_Dxy_S_Fxyz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_Fxyz_S_vrr;
      I_ERI_Dxz_S_Fxyz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_Fxyz_S_vrr;
      I_ERI_D2y_S_Fxyz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_Fxyz_S_vrr;
      I_ERI_Dyz_S_Fxyz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_Fxyz_S_vrr;
      I_ERI_D2z_S_Fxyz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_Fxyz_S_vrr;
      I_ERI_D2x_S_Fx2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_Fx2z_S_vrr;
      I_ERI_Dxy_S_Fx2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_Fx2z_S_vrr;
      I_ERI_Dxz_S_Fx2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_Fx2z_S_vrr;
      I_ERI_D2y_S_Fx2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_Fx2z_S_vrr;
      I_ERI_Dyz_S_Fx2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_Fx2z_S_vrr;
      I_ERI_D2z_S_Fx2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_Fx2z_S_vrr;
      I_ERI_D2x_S_F3y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_F3y_S_vrr;
      I_ERI_Dxy_S_F3y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_F3y_S_vrr;
      I_ERI_Dxz_S_F3y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_F3y_S_vrr;
      I_ERI_D2y_S_F3y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_F3y_S_vrr;
      I_ERI_Dyz_S_F3y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_F3y_S_vrr;
      I_ERI_D2z_S_F3y_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_F3y_S_vrr;
      I_ERI_D2x_S_F2yz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_F2yz_S_vrr;
      I_ERI_Dxy_S_F2yz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_F2yz_S_vrr;
      I_ERI_Dxz_S_F2yz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_F2yz_S_vrr;
      I_ERI_D2y_S_F2yz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_F2yz_S_vrr;
      I_ERI_Dyz_S_F2yz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_F2yz_S_vrr;
      I_ERI_D2z_S_F2yz_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_F2yz_S_vrr;
      I_ERI_D2x_S_Fy2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_Fy2z_S_vrr;
      I_ERI_Dxy_S_Fy2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_Fy2z_S_vrr;
      I_ERI_Dxz_S_Fy2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_Fy2z_S_vrr;
      I_ERI_D2y_S_Fy2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_Fy2z_S_vrr;
      I_ERI_Dyz_S_Fy2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_Fy2z_S_vrr;
      I_ERI_D2z_S_Fy2z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_Fy2z_S_vrr;
      I_ERI_D2x_S_F3z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2x_S_F3z_S_vrr;
      I_ERI_Dxy_S_F3z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxy_S_F3z_S_vrr;
      I_ERI_Dxz_S_F3z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dxz_S_F3z_S_vrr;
      I_ERI_D2y_S_F3z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2y_S_F3z_S_vrr;
      I_ERI_Dyz_S_F3z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_Dyz_S_F3z_S_vrr;
      I_ERI_D2z_S_F3z_S_C1002001002 += SQ_ERI_D_S_F_S_C1002001002_coefs*I_ERI_D2z_S_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_D_S_C1002001002
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_ERI_D_S_D_S_C1002001002_coefs = ic2*jc2_1;
      I_ERI_D2x_S_D2x_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2x_S_D2x_S_vrr;
      I_ERI_Dxy_S_D2x_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxy_S_D2x_S_vrr;
      I_ERI_Dxz_S_D2x_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxz_S_D2x_S_vrr;
      I_ERI_D2y_S_D2x_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2y_S_D2x_S_vrr;
      I_ERI_Dyz_S_D2x_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dyz_S_D2x_S_vrr;
      I_ERI_D2z_S_D2x_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2z_S_D2x_S_vrr;
      I_ERI_D2x_S_Dxy_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2x_S_Dxy_S_vrr;
      I_ERI_Dxy_S_Dxy_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxy_S_Dxy_S_vrr;
      I_ERI_Dxz_S_Dxy_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxz_S_Dxy_S_vrr;
      I_ERI_D2y_S_Dxy_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2y_S_Dxy_S_vrr;
      I_ERI_Dyz_S_Dxy_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dyz_S_Dxy_S_vrr;
      I_ERI_D2z_S_Dxy_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2z_S_Dxy_S_vrr;
      I_ERI_D2x_S_Dxz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2x_S_Dxz_S_vrr;
      I_ERI_Dxy_S_Dxz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxy_S_Dxz_S_vrr;
      I_ERI_Dxz_S_Dxz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxz_S_Dxz_S_vrr;
      I_ERI_D2y_S_Dxz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2y_S_Dxz_S_vrr;
      I_ERI_Dyz_S_Dxz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dyz_S_Dxz_S_vrr;
      I_ERI_D2z_S_Dxz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2z_S_Dxz_S_vrr;
      I_ERI_D2x_S_D2y_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2x_S_D2y_S_vrr;
      I_ERI_Dxy_S_D2y_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxy_S_D2y_S_vrr;
      I_ERI_Dxz_S_D2y_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxz_S_D2y_S_vrr;
      I_ERI_D2y_S_D2y_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2y_S_D2y_S_vrr;
      I_ERI_Dyz_S_D2y_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dyz_S_D2y_S_vrr;
      I_ERI_D2z_S_D2y_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2z_S_D2y_S_vrr;
      I_ERI_D2x_S_Dyz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2x_S_Dyz_S_vrr;
      I_ERI_Dxy_S_Dyz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxy_S_Dyz_S_vrr;
      I_ERI_Dxz_S_Dyz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxz_S_Dyz_S_vrr;
      I_ERI_D2y_S_Dyz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2y_S_Dyz_S_vrr;
      I_ERI_Dyz_S_Dyz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dyz_S_Dyz_S_vrr;
      I_ERI_D2z_S_Dyz_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2z_S_Dyz_S_vrr;
      I_ERI_D2x_S_D2z_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2x_S_D2z_S_vrr;
      I_ERI_Dxy_S_D2z_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxy_S_D2z_S_vrr;
      I_ERI_Dxz_S_D2z_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dxz_S_D2z_S_vrr;
      I_ERI_D2y_S_D2z_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2y_S_D2z_S_vrr;
      I_ERI_Dyz_S_D2z_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_Dyz_S_D2z_S_vrr;
      I_ERI_D2z_S_D2z_S_C1002001002 += SQ_ERI_D_S_D_S_C1002001002_coefs*I_ERI_D2z_S_D2z_S_vrr;
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
   * shell quartet name: SQ_ERI_D_S_D_P_C1002001002
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_D_S_F_S_C1002001002
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C1002001002
   ************************************************************/
  Double I_ERI_D2x_S_D2x_Px_C1002001002 = I_ERI_D2x_S_F3x_S_C1002001002+CDX*I_ERI_D2x_S_D2x_S_C1002001002;
  Double I_ERI_Dxy_S_D2x_Px_C1002001002 = I_ERI_Dxy_S_F3x_S_C1002001002+CDX*I_ERI_Dxy_S_D2x_S_C1002001002;
  Double I_ERI_Dxz_S_D2x_Px_C1002001002 = I_ERI_Dxz_S_F3x_S_C1002001002+CDX*I_ERI_Dxz_S_D2x_S_C1002001002;
  Double I_ERI_D2y_S_D2x_Px_C1002001002 = I_ERI_D2y_S_F3x_S_C1002001002+CDX*I_ERI_D2y_S_D2x_S_C1002001002;
  Double I_ERI_Dyz_S_D2x_Px_C1002001002 = I_ERI_Dyz_S_F3x_S_C1002001002+CDX*I_ERI_Dyz_S_D2x_S_C1002001002;
  Double I_ERI_D2z_S_D2x_Px_C1002001002 = I_ERI_D2z_S_F3x_S_C1002001002+CDX*I_ERI_D2z_S_D2x_S_C1002001002;
  Double I_ERI_D2x_S_Dxy_Px_C1002001002 = I_ERI_D2x_S_F2xy_S_C1002001002+CDX*I_ERI_D2x_S_Dxy_S_C1002001002;
  Double I_ERI_Dxy_S_Dxy_Px_C1002001002 = I_ERI_Dxy_S_F2xy_S_C1002001002+CDX*I_ERI_Dxy_S_Dxy_S_C1002001002;
  Double I_ERI_Dxz_S_Dxy_Px_C1002001002 = I_ERI_Dxz_S_F2xy_S_C1002001002+CDX*I_ERI_Dxz_S_Dxy_S_C1002001002;
  Double I_ERI_D2y_S_Dxy_Px_C1002001002 = I_ERI_D2y_S_F2xy_S_C1002001002+CDX*I_ERI_D2y_S_Dxy_S_C1002001002;
  Double I_ERI_Dyz_S_Dxy_Px_C1002001002 = I_ERI_Dyz_S_F2xy_S_C1002001002+CDX*I_ERI_Dyz_S_Dxy_S_C1002001002;
  Double I_ERI_D2z_S_Dxy_Px_C1002001002 = I_ERI_D2z_S_F2xy_S_C1002001002+CDX*I_ERI_D2z_S_Dxy_S_C1002001002;
  Double I_ERI_D2x_S_Dxz_Px_C1002001002 = I_ERI_D2x_S_F2xz_S_C1002001002+CDX*I_ERI_D2x_S_Dxz_S_C1002001002;
  Double I_ERI_Dxy_S_Dxz_Px_C1002001002 = I_ERI_Dxy_S_F2xz_S_C1002001002+CDX*I_ERI_Dxy_S_Dxz_S_C1002001002;
  Double I_ERI_Dxz_S_Dxz_Px_C1002001002 = I_ERI_Dxz_S_F2xz_S_C1002001002+CDX*I_ERI_Dxz_S_Dxz_S_C1002001002;
  Double I_ERI_D2y_S_Dxz_Px_C1002001002 = I_ERI_D2y_S_F2xz_S_C1002001002+CDX*I_ERI_D2y_S_Dxz_S_C1002001002;
  Double I_ERI_Dyz_S_Dxz_Px_C1002001002 = I_ERI_Dyz_S_F2xz_S_C1002001002+CDX*I_ERI_Dyz_S_Dxz_S_C1002001002;
  Double I_ERI_D2z_S_Dxz_Px_C1002001002 = I_ERI_D2z_S_F2xz_S_C1002001002+CDX*I_ERI_D2z_S_Dxz_S_C1002001002;
  Double I_ERI_D2x_S_D2y_Px_C1002001002 = I_ERI_D2x_S_Fx2y_S_C1002001002+CDX*I_ERI_D2x_S_D2y_S_C1002001002;
  Double I_ERI_Dxy_S_D2y_Px_C1002001002 = I_ERI_Dxy_S_Fx2y_S_C1002001002+CDX*I_ERI_Dxy_S_D2y_S_C1002001002;
  Double I_ERI_Dxz_S_D2y_Px_C1002001002 = I_ERI_Dxz_S_Fx2y_S_C1002001002+CDX*I_ERI_Dxz_S_D2y_S_C1002001002;
  Double I_ERI_D2y_S_D2y_Px_C1002001002 = I_ERI_D2y_S_Fx2y_S_C1002001002+CDX*I_ERI_D2y_S_D2y_S_C1002001002;
  Double I_ERI_Dyz_S_D2y_Px_C1002001002 = I_ERI_Dyz_S_Fx2y_S_C1002001002+CDX*I_ERI_Dyz_S_D2y_S_C1002001002;
  Double I_ERI_D2z_S_D2y_Px_C1002001002 = I_ERI_D2z_S_Fx2y_S_C1002001002+CDX*I_ERI_D2z_S_D2y_S_C1002001002;
  Double I_ERI_D2x_S_Dyz_Px_C1002001002 = I_ERI_D2x_S_Fxyz_S_C1002001002+CDX*I_ERI_D2x_S_Dyz_S_C1002001002;
  Double I_ERI_Dxy_S_Dyz_Px_C1002001002 = I_ERI_Dxy_S_Fxyz_S_C1002001002+CDX*I_ERI_Dxy_S_Dyz_S_C1002001002;
  Double I_ERI_Dxz_S_Dyz_Px_C1002001002 = I_ERI_Dxz_S_Fxyz_S_C1002001002+CDX*I_ERI_Dxz_S_Dyz_S_C1002001002;
  Double I_ERI_D2y_S_Dyz_Px_C1002001002 = I_ERI_D2y_S_Fxyz_S_C1002001002+CDX*I_ERI_D2y_S_Dyz_S_C1002001002;
  Double I_ERI_Dyz_S_Dyz_Px_C1002001002 = I_ERI_Dyz_S_Fxyz_S_C1002001002+CDX*I_ERI_Dyz_S_Dyz_S_C1002001002;
  Double I_ERI_D2z_S_Dyz_Px_C1002001002 = I_ERI_D2z_S_Fxyz_S_C1002001002+CDX*I_ERI_D2z_S_Dyz_S_C1002001002;
  Double I_ERI_D2x_S_D2z_Px_C1002001002 = I_ERI_D2x_S_Fx2z_S_C1002001002+CDX*I_ERI_D2x_S_D2z_S_C1002001002;
  Double I_ERI_Dxy_S_D2z_Px_C1002001002 = I_ERI_Dxy_S_Fx2z_S_C1002001002+CDX*I_ERI_Dxy_S_D2z_S_C1002001002;
  Double I_ERI_Dxz_S_D2z_Px_C1002001002 = I_ERI_Dxz_S_Fx2z_S_C1002001002+CDX*I_ERI_Dxz_S_D2z_S_C1002001002;
  Double I_ERI_D2y_S_D2z_Px_C1002001002 = I_ERI_D2y_S_Fx2z_S_C1002001002+CDX*I_ERI_D2y_S_D2z_S_C1002001002;
  Double I_ERI_Dyz_S_D2z_Px_C1002001002 = I_ERI_Dyz_S_Fx2z_S_C1002001002+CDX*I_ERI_Dyz_S_D2z_S_C1002001002;
  Double I_ERI_D2z_S_D2z_Px_C1002001002 = I_ERI_D2z_S_Fx2z_S_C1002001002+CDX*I_ERI_D2z_S_D2z_S_C1002001002;
  Double I_ERI_D2x_S_D2x_Py_C1002001002 = I_ERI_D2x_S_F2xy_S_C1002001002+CDY*I_ERI_D2x_S_D2x_S_C1002001002;
  Double I_ERI_Dxy_S_D2x_Py_C1002001002 = I_ERI_Dxy_S_F2xy_S_C1002001002+CDY*I_ERI_Dxy_S_D2x_S_C1002001002;
  Double I_ERI_Dxz_S_D2x_Py_C1002001002 = I_ERI_Dxz_S_F2xy_S_C1002001002+CDY*I_ERI_Dxz_S_D2x_S_C1002001002;
  Double I_ERI_D2y_S_D2x_Py_C1002001002 = I_ERI_D2y_S_F2xy_S_C1002001002+CDY*I_ERI_D2y_S_D2x_S_C1002001002;
  Double I_ERI_Dyz_S_D2x_Py_C1002001002 = I_ERI_Dyz_S_F2xy_S_C1002001002+CDY*I_ERI_Dyz_S_D2x_S_C1002001002;
  Double I_ERI_D2z_S_D2x_Py_C1002001002 = I_ERI_D2z_S_F2xy_S_C1002001002+CDY*I_ERI_D2z_S_D2x_S_C1002001002;
  Double I_ERI_D2x_S_Dxy_Py_C1002001002 = I_ERI_D2x_S_Fx2y_S_C1002001002+CDY*I_ERI_D2x_S_Dxy_S_C1002001002;
  Double I_ERI_Dxy_S_Dxy_Py_C1002001002 = I_ERI_Dxy_S_Fx2y_S_C1002001002+CDY*I_ERI_Dxy_S_Dxy_S_C1002001002;
  Double I_ERI_Dxz_S_Dxy_Py_C1002001002 = I_ERI_Dxz_S_Fx2y_S_C1002001002+CDY*I_ERI_Dxz_S_Dxy_S_C1002001002;
  Double I_ERI_D2y_S_Dxy_Py_C1002001002 = I_ERI_D2y_S_Fx2y_S_C1002001002+CDY*I_ERI_D2y_S_Dxy_S_C1002001002;
  Double I_ERI_Dyz_S_Dxy_Py_C1002001002 = I_ERI_Dyz_S_Fx2y_S_C1002001002+CDY*I_ERI_Dyz_S_Dxy_S_C1002001002;
  Double I_ERI_D2z_S_Dxy_Py_C1002001002 = I_ERI_D2z_S_Fx2y_S_C1002001002+CDY*I_ERI_D2z_S_Dxy_S_C1002001002;
  Double I_ERI_D2x_S_Dxz_Py_C1002001002 = I_ERI_D2x_S_Fxyz_S_C1002001002+CDY*I_ERI_D2x_S_Dxz_S_C1002001002;
  Double I_ERI_Dxy_S_Dxz_Py_C1002001002 = I_ERI_Dxy_S_Fxyz_S_C1002001002+CDY*I_ERI_Dxy_S_Dxz_S_C1002001002;
  Double I_ERI_Dxz_S_Dxz_Py_C1002001002 = I_ERI_Dxz_S_Fxyz_S_C1002001002+CDY*I_ERI_Dxz_S_Dxz_S_C1002001002;
  Double I_ERI_D2y_S_Dxz_Py_C1002001002 = I_ERI_D2y_S_Fxyz_S_C1002001002+CDY*I_ERI_D2y_S_Dxz_S_C1002001002;
  Double I_ERI_Dyz_S_Dxz_Py_C1002001002 = I_ERI_Dyz_S_Fxyz_S_C1002001002+CDY*I_ERI_Dyz_S_Dxz_S_C1002001002;
  Double I_ERI_D2z_S_Dxz_Py_C1002001002 = I_ERI_D2z_S_Fxyz_S_C1002001002+CDY*I_ERI_D2z_S_Dxz_S_C1002001002;
  Double I_ERI_D2x_S_D2y_Py_C1002001002 = I_ERI_D2x_S_F3y_S_C1002001002+CDY*I_ERI_D2x_S_D2y_S_C1002001002;
  Double I_ERI_Dxy_S_D2y_Py_C1002001002 = I_ERI_Dxy_S_F3y_S_C1002001002+CDY*I_ERI_Dxy_S_D2y_S_C1002001002;
  Double I_ERI_Dxz_S_D2y_Py_C1002001002 = I_ERI_Dxz_S_F3y_S_C1002001002+CDY*I_ERI_Dxz_S_D2y_S_C1002001002;
  Double I_ERI_D2y_S_D2y_Py_C1002001002 = I_ERI_D2y_S_F3y_S_C1002001002+CDY*I_ERI_D2y_S_D2y_S_C1002001002;
  Double I_ERI_Dyz_S_D2y_Py_C1002001002 = I_ERI_Dyz_S_F3y_S_C1002001002+CDY*I_ERI_Dyz_S_D2y_S_C1002001002;
  Double I_ERI_D2z_S_D2y_Py_C1002001002 = I_ERI_D2z_S_F3y_S_C1002001002+CDY*I_ERI_D2z_S_D2y_S_C1002001002;
  Double I_ERI_D2x_S_Dyz_Py_C1002001002 = I_ERI_D2x_S_F2yz_S_C1002001002+CDY*I_ERI_D2x_S_Dyz_S_C1002001002;
  Double I_ERI_Dxy_S_Dyz_Py_C1002001002 = I_ERI_Dxy_S_F2yz_S_C1002001002+CDY*I_ERI_Dxy_S_Dyz_S_C1002001002;
  Double I_ERI_Dxz_S_Dyz_Py_C1002001002 = I_ERI_Dxz_S_F2yz_S_C1002001002+CDY*I_ERI_Dxz_S_Dyz_S_C1002001002;
  Double I_ERI_D2y_S_Dyz_Py_C1002001002 = I_ERI_D2y_S_F2yz_S_C1002001002+CDY*I_ERI_D2y_S_Dyz_S_C1002001002;
  Double I_ERI_Dyz_S_Dyz_Py_C1002001002 = I_ERI_Dyz_S_F2yz_S_C1002001002+CDY*I_ERI_Dyz_S_Dyz_S_C1002001002;
  Double I_ERI_D2z_S_Dyz_Py_C1002001002 = I_ERI_D2z_S_F2yz_S_C1002001002+CDY*I_ERI_D2z_S_Dyz_S_C1002001002;
  Double I_ERI_D2x_S_D2z_Py_C1002001002 = I_ERI_D2x_S_Fy2z_S_C1002001002+CDY*I_ERI_D2x_S_D2z_S_C1002001002;
  Double I_ERI_Dxy_S_D2z_Py_C1002001002 = I_ERI_Dxy_S_Fy2z_S_C1002001002+CDY*I_ERI_Dxy_S_D2z_S_C1002001002;
  Double I_ERI_Dxz_S_D2z_Py_C1002001002 = I_ERI_Dxz_S_Fy2z_S_C1002001002+CDY*I_ERI_Dxz_S_D2z_S_C1002001002;
  Double I_ERI_D2y_S_D2z_Py_C1002001002 = I_ERI_D2y_S_Fy2z_S_C1002001002+CDY*I_ERI_D2y_S_D2z_S_C1002001002;
  Double I_ERI_Dyz_S_D2z_Py_C1002001002 = I_ERI_Dyz_S_Fy2z_S_C1002001002+CDY*I_ERI_Dyz_S_D2z_S_C1002001002;
  Double I_ERI_D2z_S_D2z_Py_C1002001002 = I_ERI_D2z_S_Fy2z_S_C1002001002+CDY*I_ERI_D2z_S_D2z_S_C1002001002;
  Double I_ERI_D2x_S_D2x_Pz_C1002001002 = I_ERI_D2x_S_F2xz_S_C1002001002+CDZ*I_ERI_D2x_S_D2x_S_C1002001002;
  Double I_ERI_Dxy_S_D2x_Pz_C1002001002 = I_ERI_Dxy_S_F2xz_S_C1002001002+CDZ*I_ERI_Dxy_S_D2x_S_C1002001002;
  Double I_ERI_Dxz_S_D2x_Pz_C1002001002 = I_ERI_Dxz_S_F2xz_S_C1002001002+CDZ*I_ERI_Dxz_S_D2x_S_C1002001002;
  Double I_ERI_D2y_S_D2x_Pz_C1002001002 = I_ERI_D2y_S_F2xz_S_C1002001002+CDZ*I_ERI_D2y_S_D2x_S_C1002001002;
  Double I_ERI_Dyz_S_D2x_Pz_C1002001002 = I_ERI_Dyz_S_F2xz_S_C1002001002+CDZ*I_ERI_Dyz_S_D2x_S_C1002001002;
  Double I_ERI_D2z_S_D2x_Pz_C1002001002 = I_ERI_D2z_S_F2xz_S_C1002001002+CDZ*I_ERI_D2z_S_D2x_S_C1002001002;
  Double I_ERI_D2x_S_Dxy_Pz_C1002001002 = I_ERI_D2x_S_Fxyz_S_C1002001002+CDZ*I_ERI_D2x_S_Dxy_S_C1002001002;
  Double I_ERI_Dxy_S_Dxy_Pz_C1002001002 = I_ERI_Dxy_S_Fxyz_S_C1002001002+CDZ*I_ERI_Dxy_S_Dxy_S_C1002001002;
  Double I_ERI_Dxz_S_Dxy_Pz_C1002001002 = I_ERI_Dxz_S_Fxyz_S_C1002001002+CDZ*I_ERI_Dxz_S_Dxy_S_C1002001002;
  Double I_ERI_D2y_S_Dxy_Pz_C1002001002 = I_ERI_D2y_S_Fxyz_S_C1002001002+CDZ*I_ERI_D2y_S_Dxy_S_C1002001002;
  Double I_ERI_Dyz_S_Dxy_Pz_C1002001002 = I_ERI_Dyz_S_Fxyz_S_C1002001002+CDZ*I_ERI_Dyz_S_Dxy_S_C1002001002;
  Double I_ERI_D2z_S_Dxy_Pz_C1002001002 = I_ERI_D2z_S_Fxyz_S_C1002001002+CDZ*I_ERI_D2z_S_Dxy_S_C1002001002;
  Double I_ERI_D2x_S_Dxz_Pz_C1002001002 = I_ERI_D2x_S_Fx2z_S_C1002001002+CDZ*I_ERI_D2x_S_Dxz_S_C1002001002;
  Double I_ERI_Dxy_S_Dxz_Pz_C1002001002 = I_ERI_Dxy_S_Fx2z_S_C1002001002+CDZ*I_ERI_Dxy_S_Dxz_S_C1002001002;
  Double I_ERI_Dxz_S_Dxz_Pz_C1002001002 = I_ERI_Dxz_S_Fx2z_S_C1002001002+CDZ*I_ERI_Dxz_S_Dxz_S_C1002001002;
  Double I_ERI_D2y_S_Dxz_Pz_C1002001002 = I_ERI_D2y_S_Fx2z_S_C1002001002+CDZ*I_ERI_D2y_S_Dxz_S_C1002001002;
  Double I_ERI_Dyz_S_Dxz_Pz_C1002001002 = I_ERI_Dyz_S_Fx2z_S_C1002001002+CDZ*I_ERI_Dyz_S_Dxz_S_C1002001002;
  Double I_ERI_D2z_S_Dxz_Pz_C1002001002 = I_ERI_D2z_S_Fx2z_S_C1002001002+CDZ*I_ERI_D2z_S_Dxz_S_C1002001002;
  Double I_ERI_D2x_S_D2y_Pz_C1002001002 = I_ERI_D2x_S_F2yz_S_C1002001002+CDZ*I_ERI_D2x_S_D2y_S_C1002001002;
  Double I_ERI_Dxy_S_D2y_Pz_C1002001002 = I_ERI_Dxy_S_F2yz_S_C1002001002+CDZ*I_ERI_Dxy_S_D2y_S_C1002001002;
  Double I_ERI_Dxz_S_D2y_Pz_C1002001002 = I_ERI_Dxz_S_F2yz_S_C1002001002+CDZ*I_ERI_Dxz_S_D2y_S_C1002001002;
  Double I_ERI_D2y_S_D2y_Pz_C1002001002 = I_ERI_D2y_S_F2yz_S_C1002001002+CDZ*I_ERI_D2y_S_D2y_S_C1002001002;
  Double I_ERI_Dyz_S_D2y_Pz_C1002001002 = I_ERI_Dyz_S_F2yz_S_C1002001002+CDZ*I_ERI_Dyz_S_D2y_S_C1002001002;
  Double I_ERI_D2z_S_D2y_Pz_C1002001002 = I_ERI_D2z_S_F2yz_S_C1002001002+CDZ*I_ERI_D2z_S_D2y_S_C1002001002;
  Double I_ERI_D2x_S_Dyz_Pz_C1002001002 = I_ERI_D2x_S_Fy2z_S_C1002001002+CDZ*I_ERI_D2x_S_Dyz_S_C1002001002;
  Double I_ERI_Dxy_S_Dyz_Pz_C1002001002 = I_ERI_Dxy_S_Fy2z_S_C1002001002+CDZ*I_ERI_Dxy_S_Dyz_S_C1002001002;
  Double I_ERI_Dxz_S_Dyz_Pz_C1002001002 = I_ERI_Dxz_S_Fy2z_S_C1002001002+CDZ*I_ERI_Dxz_S_Dyz_S_C1002001002;
  Double I_ERI_D2y_S_Dyz_Pz_C1002001002 = I_ERI_D2y_S_Fy2z_S_C1002001002+CDZ*I_ERI_D2y_S_Dyz_S_C1002001002;
  Double I_ERI_Dyz_S_Dyz_Pz_C1002001002 = I_ERI_Dyz_S_Fy2z_S_C1002001002+CDZ*I_ERI_Dyz_S_Dyz_S_C1002001002;
  Double I_ERI_D2z_S_Dyz_Pz_C1002001002 = I_ERI_D2z_S_Fy2z_S_C1002001002+CDZ*I_ERI_D2z_S_Dyz_S_C1002001002;
  Double I_ERI_D2x_S_D2z_Pz_C1002001002 = I_ERI_D2x_S_F3z_S_C1002001002+CDZ*I_ERI_D2x_S_D2z_S_C1002001002;
  Double I_ERI_Dxy_S_D2z_Pz_C1002001002 = I_ERI_Dxy_S_F3z_S_C1002001002+CDZ*I_ERI_Dxy_S_D2z_S_C1002001002;
  Double I_ERI_Dxz_S_D2z_Pz_C1002001002 = I_ERI_Dxz_S_F3z_S_C1002001002+CDZ*I_ERI_Dxz_S_D2z_S_C1002001002;
  Double I_ERI_D2y_S_D2z_Pz_C1002001002 = I_ERI_D2y_S_F3z_S_C1002001002+CDZ*I_ERI_D2y_S_D2z_S_C1002001002;
  Double I_ERI_Dyz_S_D2z_Pz_C1002001002 = I_ERI_Dyz_S_F3z_S_C1002001002+CDZ*I_ERI_Dyz_S_D2z_S_C1002001002;
  Double I_ERI_D2z_S_D2z_Pz_C1002001002 = I_ERI_D2z_S_F3z_S_C1002001002+CDZ*I_ERI_D2z_S_D2z_S_C1002001002;

  /************************************************************
   * shell quartet name: SQ_ERI_F_S_D_P_C1002001002
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_F_S_C1002001002
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C1002001002
   ************************************************************/
  Double I_ERI_F3x_S_D2x_Px_C1002001002 = I_ERI_F3x_S_F3x_S_C1002001002+CDX*I_ERI_F3x_S_D2x_S_C1002001002;
  Double I_ERI_F2xy_S_D2x_Px_C1002001002 = I_ERI_F2xy_S_F3x_S_C1002001002+CDX*I_ERI_F2xy_S_D2x_S_C1002001002;
  Double I_ERI_F2xz_S_D2x_Px_C1002001002 = I_ERI_F2xz_S_F3x_S_C1002001002+CDX*I_ERI_F2xz_S_D2x_S_C1002001002;
  Double I_ERI_Fx2y_S_D2x_Px_C1002001002 = I_ERI_Fx2y_S_F3x_S_C1002001002+CDX*I_ERI_Fx2y_S_D2x_S_C1002001002;
  Double I_ERI_Fxyz_S_D2x_Px_C1002001002 = I_ERI_Fxyz_S_F3x_S_C1002001002+CDX*I_ERI_Fxyz_S_D2x_S_C1002001002;
  Double I_ERI_Fx2z_S_D2x_Px_C1002001002 = I_ERI_Fx2z_S_F3x_S_C1002001002+CDX*I_ERI_Fx2z_S_D2x_S_C1002001002;
  Double I_ERI_F3y_S_D2x_Px_C1002001002 = I_ERI_F3y_S_F3x_S_C1002001002+CDX*I_ERI_F3y_S_D2x_S_C1002001002;
  Double I_ERI_F2yz_S_D2x_Px_C1002001002 = I_ERI_F2yz_S_F3x_S_C1002001002+CDX*I_ERI_F2yz_S_D2x_S_C1002001002;
  Double I_ERI_Fy2z_S_D2x_Px_C1002001002 = I_ERI_Fy2z_S_F3x_S_C1002001002+CDX*I_ERI_Fy2z_S_D2x_S_C1002001002;
  Double I_ERI_F3z_S_D2x_Px_C1002001002 = I_ERI_F3z_S_F3x_S_C1002001002+CDX*I_ERI_F3z_S_D2x_S_C1002001002;
  Double I_ERI_F3x_S_Dxy_Px_C1002001002 = I_ERI_F3x_S_F2xy_S_C1002001002+CDX*I_ERI_F3x_S_Dxy_S_C1002001002;
  Double I_ERI_F2xy_S_Dxy_Px_C1002001002 = I_ERI_F2xy_S_F2xy_S_C1002001002+CDX*I_ERI_F2xy_S_Dxy_S_C1002001002;
  Double I_ERI_F2xz_S_Dxy_Px_C1002001002 = I_ERI_F2xz_S_F2xy_S_C1002001002+CDX*I_ERI_F2xz_S_Dxy_S_C1002001002;
  Double I_ERI_Fx2y_S_Dxy_Px_C1002001002 = I_ERI_Fx2y_S_F2xy_S_C1002001002+CDX*I_ERI_Fx2y_S_Dxy_S_C1002001002;
  Double I_ERI_Fxyz_S_Dxy_Px_C1002001002 = I_ERI_Fxyz_S_F2xy_S_C1002001002+CDX*I_ERI_Fxyz_S_Dxy_S_C1002001002;
  Double I_ERI_Fx2z_S_Dxy_Px_C1002001002 = I_ERI_Fx2z_S_F2xy_S_C1002001002+CDX*I_ERI_Fx2z_S_Dxy_S_C1002001002;
  Double I_ERI_F3y_S_Dxy_Px_C1002001002 = I_ERI_F3y_S_F2xy_S_C1002001002+CDX*I_ERI_F3y_S_Dxy_S_C1002001002;
  Double I_ERI_F2yz_S_Dxy_Px_C1002001002 = I_ERI_F2yz_S_F2xy_S_C1002001002+CDX*I_ERI_F2yz_S_Dxy_S_C1002001002;
  Double I_ERI_Fy2z_S_Dxy_Px_C1002001002 = I_ERI_Fy2z_S_F2xy_S_C1002001002+CDX*I_ERI_Fy2z_S_Dxy_S_C1002001002;
  Double I_ERI_F3z_S_Dxy_Px_C1002001002 = I_ERI_F3z_S_F2xy_S_C1002001002+CDX*I_ERI_F3z_S_Dxy_S_C1002001002;
  Double I_ERI_F3x_S_Dxz_Px_C1002001002 = I_ERI_F3x_S_F2xz_S_C1002001002+CDX*I_ERI_F3x_S_Dxz_S_C1002001002;
  Double I_ERI_F2xy_S_Dxz_Px_C1002001002 = I_ERI_F2xy_S_F2xz_S_C1002001002+CDX*I_ERI_F2xy_S_Dxz_S_C1002001002;
  Double I_ERI_F2xz_S_Dxz_Px_C1002001002 = I_ERI_F2xz_S_F2xz_S_C1002001002+CDX*I_ERI_F2xz_S_Dxz_S_C1002001002;
  Double I_ERI_Fx2y_S_Dxz_Px_C1002001002 = I_ERI_Fx2y_S_F2xz_S_C1002001002+CDX*I_ERI_Fx2y_S_Dxz_S_C1002001002;
  Double I_ERI_Fxyz_S_Dxz_Px_C1002001002 = I_ERI_Fxyz_S_F2xz_S_C1002001002+CDX*I_ERI_Fxyz_S_Dxz_S_C1002001002;
  Double I_ERI_Fx2z_S_Dxz_Px_C1002001002 = I_ERI_Fx2z_S_F2xz_S_C1002001002+CDX*I_ERI_Fx2z_S_Dxz_S_C1002001002;
  Double I_ERI_F3y_S_Dxz_Px_C1002001002 = I_ERI_F3y_S_F2xz_S_C1002001002+CDX*I_ERI_F3y_S_Dxz_S_C1002001002;
  Double I_ERI_F2yz_S_Dxz_Px_C1002001002 = I_ERI_F2yz_S_F2xz_S_C1002001002+CDX*I_ERI_F2yz_S_Dxz_S_C1002001002;
  Double I_ERI_Fy2z_S_Dxz_Px_C1002001002 = I_ERI_Fy2z_S_F2xz_S_C1002001002+CDX*I_ERI_Fy2z_S_Dxz_S_C1002001002;
  Double I_ERI_F3z_S_Dxz_Px_C1002001002 = I_ERI_F3z_S_F2xz_S_C1002001002+CDX*I_ERI_F3z_S_Dxz_S_C1002001002;
  Double I_ERI_F3x_S_D2y_Px_C1002001002 = I_ERI_F3x_S_Fx2y_S_C1002001002+CDX*I_ERI_F3x_S_D2y_S_C1002001002;
  Double I_ERI_F2xy_S_D2y_Px_C1002001002 = I_ERI_F2xy_S_Fx2y_S_C1002001002+CDX*I_ERI_F2xy_S_D2y_S_C1002001002;
  Double I_ERI_F2xz_S_D2y_Px_C1002001002 = I_ERI_F2xz_S_Fx2y_S_C1002001002+CDX*I_ERI_F2xz_S_D2y_S_C1002001002;
  Double I_ERI_Fx2y_S_D2y_Px_C1002001002 = I_ERI_Fx2y_S_Fx2y_S_C1002001002+CDX*I_ERI_Fx2y_S_D2y_S_C1002001002;
  Double I_ERI_Fxyz_S_D2y_Px_C1002001002 = I_ERI_Fxyz_S_Fx2y_S_C1002001002+CDX*I_ERI_Fxyz_S_D2y_S_C1002001002;
  Double I_ERI_Fx2z_S_D2y_Px_C1002001002 = I_ERI_Fx2z_S_Fx2y_S_C1002001002+CDX*I_ERI_Fx2z_S_D2y_S_C1002001002;
  Double I_ERI_F3y_S_D2y_Px_C1002001002 = I_ERI_F3y_S_Fx2y_S_C1002001002+CDX*I_ERI_F3y_S_D2y_S_C1002001002;
  Double I_ERI_F2yz_S_D2y_Px_C1002001002 = I_ERI_F2yz_S_Fx2y_S_C1002001002+CDX*I_ERI_F2yz_S_D2y_S_C1002001002;
  Double I_ERI_Fy2z_S_D2y_Px_C1002001002 = I_ERI_Fy2z_S_Fx2y_S_C1002001002+CDX*I_ERI_Fy2z_S_D2y_S_C1002001002;
  Double I_ERI_F3z_S_D2y_Px_C1002001002 = I_ERI_F3z_S_Fx2y_S_C1002001002+CDX*I_ERI_F3z_S_D2y_S_C1002001002;
  Double I_ERI_F3x_S_Dyz_Px_C1002001002 = I_ERI_F3x_S_Fxyz_S_C1002001002+CDX*I_ERI_F3x_S_Dyz_S_C1002001002;
  Double I_ERI_F2xy_S_Dyz_Px_C1002001002 = I_ERI_F2xy_S_Fxyz_S_C1002001002+CDX*I_ERI_F2xy_S_Dyz_S_C1002001002;
  Double I_ERI_F2xz_S_Dyz_Px_C1002001002 = I_ERI_F2xz_S_Fxyz_S_C1002001002+CDX*I_ERI_F2xz_S_Dyz_S_C1002001002;
  Double I_ERI_Fx2y_S_Dyz_Px_C1002001002 = I_ERI_Fx2y_S_Fxyz_S_C1002001002+CDX*I_ERI_Fx2y_S_Dyz_S_C1002001002;
  Double I_ERI_Fxyz_S_Dyz_Px_C1002001002 = I_ERI_Fxyz_S_Fxyz_S_C1002001002+CDX*I_ERI_Fxyz_S_Dyz_S_C1002001002;
  Double I_ERI_Fx2z_S_Dyz_Px_C1002001002 = I_ERI_Fx2z_S_Fxyz_S_C1002001002+CDX*I_ERI_Fx2z_S_Dyz_S_C1002001002;
  Double I_ERI_F3y_S_Dyz_Px_C1002001002 = I_ERI_F3y_S_Fxyz_S_C1002001002+CDX*I_ERI_F3y_S_Dyz_S_C1002001002;
  Double I_ERI_F2yz_S_Dyz_Px_C1002001002 = I_ERI_F2yz_S_Fxyz_S_C1002001002+CDX*I_ERI_F2yz_S_Dyz_S_C1002001002;
  Double I_ERI_Fy2z_S_Dyz_Px_C1002001002 = I_ERI_Fy2z_S_Fxyz_S_C1002001002+CDX*I_ERI_Fy2z_S_Dyz_S_C1002001002;
  Double I_ERI_F3z_S_Dyz_Px_C1002001002 = I_ERI_F3z_S_Fxyz_S_C1002001002+CDX*I_ERI_F3z_S_Dyz_S_C1002001002;
  Double I_ERI_F3x_S_D2z_Px_C1002001002 = I_ERI_F3x_S_Fx2z_S_C1002001002+CDX*I_ERI_F3x_S_D2z_S_C1002001002;
  Double I_ERI_F2xy_S_D2z_Px_C1002001002 = I_ERI_F2xy_S_Fx2z_S_C1002001002+CDX*I_ERI_F2xy_S_D2z_S_C1002001002;
  Double I_ERI_F2xz_S_D2z_Px_C1002001002 = I_ERI_F2xz_S_Fx2z_S_C1002001002+CDX*I_ERI_F2xz_S_D2z_S_C1002001002;
  Double I_ERI_Fx2y_S_D2z_Px_C1002001002 = I_ERI_Fx2y_S_Fx2z_S_C1002001002+CDX*I_ERI_Fx2y_S_D2z_S_C1002001002;
  Double I_ERI_Fxyz_S_D2z_Px_C1002001002 = I_ERI_Fxyz_S_Fx2z_S_C1002001002+CDX*I_ERI_Fxyz_S_D2z_S_C1002001002;
  Double I_ERI_Fx2z_S_D2z_Px_C1002001002 = I_ERI_Fx2z_S_Fx2z_S_C1002001002+CDX*I_ERI_Fx2z_S_D2z_S_C1002001002;
  Double I_ERI_F3y_S_D2z_Px_C1002001002 = I_ERI_F3y_S_Fx2z_S_C1002001002+CDX*I_ERI_F3y_S_D2z_S_C1002001002;
  Double I_ERI_F2yz_S_D2z_Px_C1002001002 = I_ERI_F2yz_S_Fx2z_S_C1002001002+CDX*I_ERI_F2yz_S_D2z_S_C1002001002;
  Double I_ERI_Fy2z_S_D2z_Px_C1002001002 = I_ERI_Fy2z_S_Fx2z_S_C1002001002+CDX*I_ERI_Fy2z_S_D2z_S_C1002001002;
  Double I_ERI_F3z_S_D2z_Px_C1002001002 = I_ERI_F3z_S_Fx2z_S_C1002001002+CDX*I_ERI_F3z_S_D2z_S_C1002001002;
  Double I_ERI_F3x_S_D2x_Py_C1002001002 = I_ERI_F3x_S_F2xy_S_C1002001002+CDY*I_ERI_F3x_S_D2x_S_C1002001002;
  Double I_ERI_F2xy_S_D2x_Py_C1002001002 = I_ERI_F2xy_S_F2xy_S_C1002001002+CDY*I_ERI_F2xy_S_D2x_S_C1002001002;
  Double I_ERI_F2xz_S_D2x_Py_C1002001002 = I_ERI_F2xz_S_F2xy_S_C1002001002+CDY*I_ERI_F2xz_S_D2x_S_C1002001002;
  Double I_ERI_Fx2y_S_D2x_Py_C1002001002 = I_ERI_Fx2y_S_F2xy_S_C1002001002+CDY*I_ERI_Fx2y_S_D2x_S_C1002001002;
  Double I_ERI_Fxyz_S_D2x_Py_C1002001002 = I_ERI_Fxyz_S_F2xy_S_C1002001002+CDY*I_ERI_Fxyz_S_D2x_S_C1002001002;
  Double I_ERI_Fx2z_S_D2x_Py_C1002001002 = I_ERI_Fx2z_S_F2xy_S_C1002001002+CDY*I_ERI_Fx2z_S_D2x_S_C1002001002;
  Double I_ERI_F3y_S_D2x_Py_C1002001002 = I_ERI_F3y_S_F2xy_S_C1002001002+CDY*I_ERI_F3y_S_D2x_S_C1002001002;
  Double I_ERI_F2yz_S_D2x_Py_C1002001002 = I_ERI_F2yz_S_F2xy_S_C1002001002+CDY*I_ERI_F2yz_S_D2x_S_C1002001002;
  Double I_ERI_Fy2z_S_D2x_Py_C1002001002 = I_ERI_Fy2z_S_F2xy_S_C1002001002+CDY*I_ERI_Fy2z_S_D2x_S_C1002001002;
  Double I_ERI_F3z_S_D2x_Py_C1002001002 = I_ERI_F3z_S_F2xy_S_C1002001002+CDY*I_ERI_F3z_S_D2x_S_C1002001002;
  Double I_ERI_F3x_S_Dxy_Py_C1002001002 = I_ERI_F3x_S_Fx2y_S_C1002001002+CDY*I_ERI_F3x_S_Dxy_S_C1002001002;
  Double I_ERI_F2xy_S_Dxy_Py_C1002001002 = I_ERI_F2xy_S_Fx2y_S_C1002001002+CDY*I_ERI_F2xy_S_Dxy_S_C1002001002;
  Double I_ERI_F2xz_S_Dxy_Py_C1002001002 = I_ERI_F2xz_S_Fx2y_S_C1002001002+CDY*I_ERI_F2xz_S_Dxy_S_C1002001002;
  Double I_ERI_Fx2y_S_Dxy_Py_C1002001002 = I_ERI_Fx2y_S_Fx2y_S_C1002001002+CDY*I_ERI_Fx2y_S_Dxy_S_C1002001002;
  Double I_ERI_Fxyz_S_Dxy_Py_C1002001002 = I_ERI_Fxyz_S_Fx2y_S_C1002001002+CDY*I_ERI_Fxyz_S_Dxy_S_C1002001002;
  Double I_ERI_Fx2z_S_Dxy_Py_C1002001002 = I_ERI_Fx2z_S_Fx2y_S_C1002001002+CDY*I_ERI_Fx2z_S_Dxy_S_C1002001002;
  Double I_ERI_F3y_S_Dxy_Py_C1002001002 = I_ERI_F3y_S_Fx2y_S_C1002001002+CDY*I_ERI_F3y_S_Dxy_S_C1002001002;
  Double I_ERI_F2yz_S_Dxy_Py_C1002001002 = I_ERI_F2yz_S_Fx2y_S_C1002001002+CDY*I_ERI_F2yz_S_Dxy_S_C1002001002;
  Double I_ERI_Fy2z_S_Dxy_Py_C1002001002 = I_ERI_Fy2z_S_Fx2y_S_C1002001002+CDY*I_ERI_Fy2z_S_Dxy_S_C1002001002;
  Double I_ERI_F3z_S_Dxy_Py_C1002001002 = I_ERI_F3z_S_Fx2y_S_C1002001002+CDY*I_ERI_F3z_S_Dxy_S_C1002001002;
  Double I_ERI_F3x_S_Dxz_Py_C1002001002 = I_ERI_F3x_S_Fxyz_S_C1002001002+CDY*I_ERI_F3x_S_Dxz_S_C1002001002;
  Double I_ERI_F2xy_S_Dxz_Py_C1002001002 = I_ERI_F2xy_S_Fxyz_S_C1002001002+CDY*I_ERI_F2xy_S_Dxz_S_C1002001002;
  Double I_ERI_F2xz_S_Dxz_Py_C1002001002 = I_ERI_F2xz_S_Fxyz_S_C1002001002+CDY*I_ERI_F2xz_S_Dxz_S_C1002001002;
  Double I_ERI_Fx2y_S_Dxz_Py_C1002001002 = I_ERI_Fx2y_S_Fxyz_S_C1002001002+CDY*I_ERI_Fx2y_S_Dxz_S_C1002001002;
  Double I_ERI_Fxyz_S_Dxz_Py_C1002001002 = I_ERI_Fxyz_S_Fxyz_S_C1002001002+CDY*I_ERI_Fxyz_S_Dxz_S_C1002001002;
  Double I_ERI_Fx2z_S_Dxz_Py_C1002001002 = I_ERI_Fx2z_S_Fxyz_S_C1002001002+CDY*I_ERI_Fx2z_S_Dxz_S_C1002001002;
  Double I_ERI_F3y_S_Dxz_Py_C1002001002 = I_ERI_F3y_S_Fxyz_S_C1002001002+CDY*I_ERI_F3y_S_Dxz_S_C1002001002;
  Double I_ERI_F2yz_S_Dxz_Py_C1002001002 = I_ERI_F2yz_S_Fxyz_S_C1002001002+CDY*I_ERI_F2yz_S_Dxz_S_C1002001002;
  Double I_ERI_Fy2z_S_Dxz_Py_C1002001002 = I_ERI_Fy2z_S_Fxyz_S_C1002001002+CDY*I_ERI_Fy2z_S_Dxz_S_C1002001002;
  Double I_ERI_F3z_S_Dxz_Py_C1002001002 = I_ERI_F3z_S_Fxyz_S_C1002001002+CDY*I_ERI_F3z_S_Dxz_S_C1002001002;
  Double I_ERI_F3x_S_D2y_Py_C1002001002 = I_ERI_F3x_S_F3y_S_C1002001002+CDY*I_ERI_F3x_S_D2y_S_C1002001002;
  Double I_ERI_F2xy_S_D2y_Py_C1002001002 = I_ERI_F2xy_S_F3y_S_C1002001002+CDY*I_ERI_F2xy_S_D2y_S_C1002001002;
  Double I_ERI_F2xz_S_D2y_Py_C1002001002 = I_ERI_F2xz_S_F3y_S_C1002001002+CDY*I_ERI_F2xz_S_D2y_S_C1002001002;
  Double I_ERI_Fx2y_S_D2y_Py_C1002001002 = I_ERI_Fx2y_S_F3y_S_C1002001002+CDY*I_ERI_Fx2y_S_D2y_S_C1002001002;
  Double I_ERI_Fxyz_S_D2y_Py_C1002001002 = I_ERI_Fxyz_S_F3y_S_C1002001002+CDY*I_ERI_Fxyz_S_D2y_S_C1002001002;
  Double I_ERI_Fx2z_S_D2y_Py_C1002001002 = I_ERI_Fx2z_S_F3y_S_C1002001002+CDY*I_ERI_Fx2z_S_D2y_S_C1002001002;
  Double I_ERI_F3y_S_D2y_Py_C1002001002 = I_ERI_F3y_S_F3y_S_C1002001002+CDY*I_ERI_F3y_S_D2y_S_C1002001002;
  Double I_ERI_F2yz_S_D2y_Py_C1002001002 = I_ERI_F2yz_S_F3y_S_C1002001002+CDY*I_ERI_F2yz_S_D2y_S_C1002001002;
  Double I_ERI_Fy2z_S_D2y_Py_C1002001002 = I_ERI_Fy2z_S_F3y_S_C1002001002+CDY*I_ERI_Fy2z_S_D2y_S_C1002001002;
  Double I_ERI_F3z_S_D2y_Py_C1002001002 = I_ERI_F3z_S_F3y_S_C1002001002+CDY*I_ERI_F3z_S_D2y_S_C1002001002;
  Double I_ERI_F3x_S_Dyz_Py_C1002001002 = I_ERI_F3x_S_F2yz_S_C1002001002+CDY*I_ERI_F3x_S_Dyz_S_C1002001002;
  Double I_ERI_F2xy_S_Dyz_Py_C1002001002 = I_ERI_F2xy_S_F2yz_S_C1002001002+CDY*I_ERI_F2xy_S_Dyz_S_C1002001002;
  Double I_ERI_F2xz_S_Dyz_Py_C1002001002 = I_ERI_F2xz_S_F2yz_S_C1002001002+CDY*I_ERI_F2xz_S_Dyz_S_C1002001002;
  Double I_ERI_Fx2y_S_Dyz_Py_C1002001002 = I_ERI_Fx2y_S_F2yz_S_C1002001002+CDY*I_ERI_Fx2y_S_Dyz_S_C1002001002;
  Double I_ERI_Fxyz_S_Dyz_Py_C1002001002 = I_ERI_Fxyz_S_F2yz_S_C1002001002+CDY*I_ERI_Fxyz_S_Dyz_S_C1002001002;
  Double I_ERI_Fx2z_S_Dyz_Py_C1002001002 = I_ERI_Fx2z_S_F2yz_S_C1002001002+CDY*I_ERI_Fx2z_S_Dyz_S_C1002001002;
  Double I_ERI_F3y_S_Dyz_Py_C1002001002 = I_ERI_F3y_S_F2yz_S_C1002001002+CDY*I_ERI_F3y_S_Dyz_S_C1002001002;
  Double I_ERI_F2yz_S_Dyz_Py_C1002001002 = I_ERI_F2yz_S_F2yz_S_C1002001002+CDY*I_ERI_F2yz_S_Dyz_S_C1002001002;
  Double I_ERI_Fy2z_S_Dyz_Py_C1002001002 = I_ERI_Fy2z_S_F2yz_S_C1002001002+CDY*I_ERI_Fy2z_S_Dyz_S_C1002001002;
  Double I_ERI_F3z_S_Dyz_Py_C1002001002 = I_ERI_F3z_S_F2yz_S_C1002001002+CDY*I_ERI_F3z_S_Dyz_S_C1002001002;
  Double I_ERI_F3x_S_D2z_Py_C1002001002 = I_ERI_F3x_S_Fy2z_S_C1002001002+CDY*I_ERI_F3x_S_D2z_S_C1002001002;
  Double I_ERI_F2xy_S_D2z_Py_C1002001002 = I_ERI_F2xy_S_Fy2z_S_C1002001002+CDY*I_ERI_F2xy_S_D2z_S_C1002001002;
  Double I_ERI_F2xz_S_D2z_Py_C1002001002 = I_ERI_F2xz_S_Fy2z_S_C1002001002+CDY*I_ERI_F2xz_S_D2z_S_C1002001002;
  Double I_ERI_Fx2y_S_D2z_Py_C1002001002 = I_ERI_Fx2y_S_Fy2z_S_C1002001002+CDY*I_ERI_Fx2y_S_D2z_S_C1002001002;
  Double I_ERI_Fxyz_S_D2z_Py_C1002001002 = I_ERI_Fxyz_S_Fy2z_S_C1002001002+CDY*I_ERI_Fxyz_S_D2z_S_C1002001002;
  Double I_ERI_Fx2z_S_D2z_Py_C1002001002 = I_ERI_Fx2z_S_Fy2z_S_C1002001002+CDY*I_ERI_Fx2z_S_D2z_S_C1002001002;
  Double I_ERI_F3y_S_D2z_Py_C1002001002 = I_ERI_F3y_S_Fy2z_S_C1002001002+CDY*I_ERI_F3y_S_D2z_S_C1002001002;
  Double I_ERI_F2yz_S_D2z_Py_C1002001002 = I_ERI_F2yz_S_Fy2z_S_C1002001002+CDY*I_ERI_F2yz_S_D2z_S_C1002001002;
  Double I_ERI_Fy2z_S_D2z_Py_C1002001002 = I_ERI_Fy2z_S_Fy2z_S_C1002001002+CDY*I_ERI_Fy2z_S_D2z_S_C1002001002;
  Double I_ERI_F3z_S_D2z_Py_C1002001002 = I_ERI_F3z_S_Fy2z_S_C1002001002+CDY*I_ERI_F3z_S_D2z_S_C1002001002;
  Double I_ERI_F3x_S_D2x_Pz_C1002001002 = I_ERI_F3x_S_F2xz_S_C1002001002+CDZ*I_ERI_F3x_S_D2x_S_C1002001002;
  Double I_ERI_F2xy_S_D2x_Pz_C1002001002 = I_ERI_F2xy_S_F2xz_S_C1002001002+CDZ*I_ERI_F2xy_S_D2x_S_C1002001002;
  Double I_ERI_F2xz_S_D2x_Pz_C1002001002 = I_ERI_F2xz_S_F2xz_S_C1002001002+CDZ*I_ERI_F2xz_S_D2x_S_C1002001002;
  Double I_ERI_Fx2y_S_D2x_Pz_C1002001002 = I_ERI_Fx2y_S_F2xz_S_C1002001002+CDZ*I_ERI_Fx2y_S_D2x_S_C1002001002;
  Double I_ERI_Fxyz_S_D2x_Pz_C1002001002 = I_ERI_Fxyz_S_F2xz_S_C1002001002+CDZ*I_ERI_Fxyz_S_D2x_S_C1002001002;
  Double I_ERI_Fx2z_S_D2x_Pz_C1002001002 = I_ERI_Fx2z_S_F2xz_S_C1002001002+CDZ*I_ERI_Fx2z_S_D2x_S_C1002001002;
  Double I_ERI_F3y_S_D2x_Pz_C1002001002 = I_ERI_F3y_S_F2xz_S_C1002001002+CDZ*I_ERI_F3y_S_D2x_S_C1002001002;
  Double I_ERI_F2yz_S_D2x_Pz_C1002001002 = I_ERI_F2yz_S_F2xz_S_C1002001002+CDZ*I_ERI_F2yz_S_D2x_S_C1002001002;
  Double I_ERI_Fy2z_S_D2x_Pz_C1002001002 = I_ERI_Fy2z_S_F2xz_S_C1002001002+CDZ*I_ERI_Fy2z_S_D2x_S_C1002001002;
  Double I_ERI_F3z_S_D2x_Pz_C1002001002 = I_ERI_F3z_S_F2xz_S_C1002001002+CDZ*I_ERI_F3z_S_D2x_S_C1002001002;
  Double I_ERI_F3x_S_Dxy_Pz_C1002001002 = I_ERI_F3x_S_Fxyz_S_C1002001002+CDZ*I_ERI_F3x_S_Dxy_S_C1002001002;
  Double I_ERI_F2xy_S_Dxy_Pz_C1002001002 = I_ERI_F2xy_S_Fxyz_S_C1002001002+CDZ*I_ERI_F2xy_S_Dxy_S_C1002001002;
  Double I_ERI_F2xz_S_Dxy_Pz_C1002001002 = I_ERI_F2xz_S_Fxyz_S_C1002001002+CDZ*I_ERI_F2xz_S_Dxy_S_C1002001002;
  Double I_ERI_Fx2y_S_Dxy_Pz_C1002001002 = I_ERI_Fx2y_S_Fxyz_S_C1002001002+CDZ*I_ERI_Fx2y_S_Dxy_S_C1002001002;
  Double I_ERI_Fxyz_S_Dxy_Pz_C1002001002 = I_ERI_Fxyz_S_Fxyz_S_C1002001002+CDZ*I_ERI_Fxyz_S_Dxy_S_C1002001002;
  Double I_ERI_Fx2z_S_Dxy_Pz_C1002001002 = I_ERI_Fx2z_S_Fxyz_S_C1002001002+CDZ*I_ERI_Fx2z_S_Dxy_S_C1002001002;
  Double I_ERI_F3y_S_Dxy_Pz_C1002001002 = I_ERI_F3y_S_Fxyz_S_C1002001002+CDZ*I_ERI_F3y_S_Dxy_S_C1002001002;
  Double I_ERI_F2yz_S_Dxy_Pz_C1002001002 = I_ERI_F2yz_S_Fxyz_S_C1002001002+CDZ*I_ERI_F2yz_S_Dxy_S_C1002001002;
  Double I_ERI_Fy2z_S_Dxy_Pz_C1002001002 = I_ERI_Fy2z_S_Fxyz_S_C1002001002+CDZ*I_ERI_Fy2z_S_Dxy_S_C1002001002;
  Double I_ERI_F3z_S_Dxy_Pz_C1002001002 = I_ERI_F3z_S_Fxyz_S_C1002001002+CDZ*I_ERI_F3z_S_Dxy_S_C1002001002;
  Double I_ERI_F3x_S_Dxz_Pz_C1002001002 = I_ERI_F3x_S_Fx2z_S_C1002001002+CDZ*I_ERI_F3x_S_Dxz_S_C1002001002;
  Double I_ERI_F2xy_S_Dxz_Pz_C1002001002 = I_ERI_F2xy_S_Fx2z_S_C1002001002+CDZ*I_ERI_F2xy_S_Dxz_S_C1002001002;
  Double I_ERI_F2xz_S_Dxz_Pz_C1002001002 = I_ERI_F2xz_S_Fx2z_S_C1002001002+CDZ*I_ERI_F2xz_S_Dxz_S_C1002001002;
  Double I_ERI_Fx2y_S_Dxz_Pz_C1002001002 = I_ERI_Fx2y_S_Fx2z_S_C1002001002+CDZ*I_ERI_Fx2y_S_Dxz_S_C1002001002;
  Double I_ERI_Fxyz_S_Dxz_Pz_C1002001002 = I_ERI_Fxyz_S_Fx2z_S_C1002001002+CDZ*I_ERI_Fxyz_S_Dxz_S_C1002001002;
  Double I_ERI_Fx2z_S_Dxz_Pz_C1002001002 = I_ERI_Fx2z_S_Fx2z_S_C1002001002+CDZ*I_ERI_Fx2z_S_Dxz_S_C1002001002;
  Double I_ERI_F3y_S_Dxz_Pz_C1002001002 = I_ERI_F3y_S_Fx2z_S_C1002001002+CDZ*I_ERI_F3y_S_Dxz_S_C1002001002;
  Double I_ERI_F2yz_S_Dxz_Pz_C1002001002 = I_ERI_F2yz_S_Fx2z_S_C1002001002+CDZ*I_ERI_F2yz_S_Dxz_S_C1002001002;
  Double I_ERI_Fy2z_S_Dxz_Pz_C1002001002 = I_ERI_Fy2z_S_Fx2z_S_C1002001002+CDZ*I_ERI_Fy2z_S_Dxz_S_C1002001002;
  Double I_ERI_F3z_S_Dxz_Pz_C1002001002 = I_ERI_F3z_S_Fx2z_S_C1002001002+CDZ*I_ERI_F3z_S_Dxz_S_C1002001002;
  Double I_ERI_F3x_S_D2y_Pz_C1002001002 = I_ERI_F3x_S_F2yz_S_C1002001002+CDZ*I_ERI_F3x_S_D2y_S_C1002001002;
  Double I_ERI_F2xy_S_D2y_Pz_C1002001002 = I_ERI_F2xy_S_F2yz_S_C1002001002+CDZ*I_ERI_F2xy_S_D2y_S_C1002001002;
  Double I_ERI_F2xz_S_D2y_Pz_C1002001002 = I_ERI_F2xz_S_F2yz_S_C1002001002+CDZ*I_ERI_F2xz_S_D2y_S_C1002001002;
  Double I_ERI_Fx2y_S_D2y_Pz_C1002001002 = I_ERI_Fx2y_S_F2yz_S_C1002001002+CDZ*I_ERI_Fx2y_S_D2y_S_C1002001002;
  Double I_ERI_Fxyz_S_D2y_Pz_C1002001002 = I_ERI_Fxyz_S_F2yz_S_C1002001002+CDZ*I_ERI_Fxyz_S_D2y_S_C1002001002;
  Double I_ERI_Fx2z_S_D2y_Pz_C1002001002 = I_ERI_Fx2z_S_F2yz_S_C1002001002+CDZ*I_ERI_Fx2z_S_D2y_S_C1002001002;
  Double I_ERI_F3y_S_D2y_Pz_C1002001002 = I_ERI_F3y_S_F2yz_S_C1002001002+CDZ*I_ERI_F3y_S_D2y_S_C1002001002;
  Double I_ERI_F2yz_S_D2y_Pz_C1002001002 = I_ERI_F2yz_S_F2yz_S_C1002001002+CDZ*I_ERI_F2yz_S_D2y_S_C1002001002;
  Double I_ERI_Fy2z_S_D2y_Pz_C1002001002 = I_ERI_Fy2z_S_F2yz_S_C1002001002+CDZ*I_ERI_Fy2z_S_D2y_S_C1002001002;
  Double I_ERI_F3z_S_D2y_Pz_C1002001002 = I_ERI_F3z_S_F2yz_S_C1002001002+CDZ*I_ERI_F3z_S_D2y_S_C1002001002;
  Double I_ERI_F3x_S_Dyz_Pz_C1002001002 = I_ERI_F3x_S_Fy2z_S_C1002001002+CDZ*I_ERI_F3x_S_Dyz_S_C1002001002;
  Double I_ERI_F2xy_S_Dyz_Pz_C1002001002 = I_ERI_F2xy_S_Fy2z_S_C1002001002+CDZ*I_ERI_F2xy_S_Dyz_S_C1002001002;
  Double I_ERI_F2xz_S_Dyz_Pz_C1002001002 = I_ERI_F2xz_S_Fy2z_S_C1002001002+CDZ*I_ERI_F2xz_S_Dyz_S_C1002001002;
  Double I_ERI_Fx2y_S_Dyz_Pz_C1002001002 = I_ERI_Fx2y_S_Fy2z_S_C1002001002+CDZ*I_ERI_Fx2y_S_Dyz_S_C1002001002;
  Double I_ERI_Fxyz_S_Dyz_Pz_C1002001002 = I_ERI_Fxyz_S_Fy2z_S_C1002001002+CDZ*I_ERI_Fxyz_S_Dyz_S_C1002001002;
  Double I_ERI_Fx2z_S_Dyz_Pz_C1002001002 = I_ERI_Fx2z_S_Fy2z_S_C1002001002+CDZ*I_ERI_Fx2z_S_Dyz_S_C1002001002;
  Double I_ERI_F3y_S_Dyz_Pz_C1002001002 = I_ERI_F3y_S_Fy2z_S_C1002001002+CDZ*I_ERI_F3y_S_Dyz_S_C1002001002;
  Double I_ERI_F2yz_S_Dyz_Pz_C1002001002 = I_ERI_F2yz_S_Fy2z_S_C1002001002+CDZ*I_ERI_F2yz_S_Dyz_S_C1002001002;
  Double I_ERI_Fy2z_S_Dyz_Pz_C1002001002 = I_ERI_Fy2z_S_Fy2z_S_C1002001002+CDZ*I_ERI_Fy2z_S_Dyz_S_C1002001002;
  Double I_ERI_F3z_S_Dyz_Pz_C1002001002 = I_ERI_F3z_S_Fy2z_S_C1002001002+CDZ*I_ERI_F3z_S_Dyz_S_C1002001002;
  Double I_ERI_F3x_S_D2z_Pz_C1002001002 = I_ERI_F3x_S_F3z_S_C1002001002+CDZ*I_ERI_F3x_S_D2z_S_C1002001002;
  Double I_ERI_F2xy_S_D2z_Pz_C1002001002 = I_ERI_F2xy_S_F3z_S_C1002001002+CDZ*I_ERI_F2xy_S_D2z_S_C1002001002;
  Double I_ERI_F2xz_S_D2z_Pz_C1002001002 = I_ERI_F2xz_S_F3z_S_C1002001002+CDZ*I_ERI_F2xz_S_D2z_S_C1002001002;
  Double I_ERI_Fx2y_S_D2z_Pz_C1002001002 = I_ERI_Fx2y_S_F3z_S_C1002001002+CDZ*I_ERI_Fx2y_S_D2z_S_C1002001002;
  Double I_ERI_Fxyz_S_D2z_Pz_C1002001002 = I_ERI_Fxyz_S_F3z_S_C1002001002+CDZ*I_ERI_Fxyz_S_D2z_S_C1002001002;
  Double I_ERI_Fx2z_S_D2z_Pz_C1002001002 = I_ERI_Fx2z_S_F3z_S_C1002001002+CDZ*I_ERI_Fx2z_S_D2z_S_C1002001002;
  Double I_ERI_F3y_S_D2z_Pz_C1002001002 = I_ERI_F3y_S_F3z_S_C1002001002+CDZ*I_ERI_F3y_S_D2z_S_C1002001002;
  Double I_ERI_F2yz_S_D2z_Pz_C1002001002 = I_ERI_F2yz_S_F3z_S_C1002001002+CDZ*I_ERI_F2yz_S_D2z_S_C1002001002;
  Double I_ERI_Fy2z_S_D2z_Pz_C1002001002 = I_ERI_Fy2z_S_F3z_S_C1002001002+CDZ*I_ERI_Fy2z_S_D2z_S_C1002001002;
  Double I_ERI_F3z_S_D2z_Pz_C1002001002 = I_ERI_F3z_S_F3z_S_C1002001002+CDZ*I_ERI_F3z_S_D2z_S_C1002001002;

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
   * shell quartet name: SQ_ERI_D_P_D_S_C2001002
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S_C2001002
   * RHS shell quartet name: SQ_ERI_D_S_D_S_C2001002
   ************************************************************/
  abcd[0] = I_ERI_F3x_S_D2x_S_C2001002+ABX*I_ERI_D2x_S_D2x_S_C2001002;
  abcd[1] = I_ERI_F2xy_S_D2x_S_C2001002+ABX*I_ERI_Dxy_S_D2x_S_C2001002;
  abcd[2] = I_ERI_F2xz_S_D2x_S_C2001002+ABX*I_ERI_Dxz_S_D2x_S_C2001002;
  abcd[3] = I_ERI_Fx2y_S_D2x_S_C2001002+ABX*I_ERI_D2y_S_D2x_S_C2001002;
  abcd[4] = I_ERI_Fxyz_S_D2x_S_C2001002+ABX*I_ERI_Dyz_S_D2x_S_C2001002;
  abcd[5] = I_ERI_Fx2z_S_D2x_S_C2001002+ABX*I_ERI_D2z_S_D2x_S_C2001002;
  abcd[6] = I_ERI_F2xy_S_D2x_S_C2001002+ABY*I_ERI_D2x_S_D2x_S_C2001002;
  abcd[7] = I_ERI_Fx2y_S_D2x_S_C2001002+ABY*I_ERI_Dxy_S_D2x_S_C2001002;
  abcd[8] = I_ERI_Fxyz_S_D2x_S_C2001002+ABY*I_ERI_Dxz_S_D2x_S_C2001002;
  abcd[9] = I_ERI_F3y_S_D2x_S_C2001002+ABY*I_ERI_D2y_S_D2x_S_C2001002;
  abcd[10] = I_ERI_F2yz_S_D2x_S_C2001002+ABY*I_ERI_Dyz_S_D2x_S_C2001002;
  abcd[11] = I_ERI_Fy2z_S_D2x_S_C2001002+ABY*I_ERI_D2z_S_D2x_S_C2001002;
  abcd[12] = I_ERI_F2xz_S_D2x_S_C2001002+ABZ*I_ERI_D2x_S_D2x_S_C2001002;
  abcd[13] = I_ERI_Fxyz_S_D2x_S_C2001002+ABZ*I_ERI_Dxy_S_D2x_S_C2001002;
  abcd[14] = I_ERI_Fx2z_S_D2x_S_C2001002+ABZ*I_ERI_Dxz_S_D2x_S_C2001002;
  abcd[15] = I_ERI_F2yz_S_D2x_S_C2001002+ABZ*I_ERI_D2y_S_D2x_S_C2001002;
  abcd[16] = I_ERI_Fy2z_S_D2x_S_C2001002+ABZ*I_ERI_Dyz_S_D2x_S_C2001002;
  abcd[17] = I_ERI_F3z_S_D2x_S_C2001002+ABZ*I_ERI_D2z_S_D2x_S_C2001002;
  abcd[18] = I_ERI_F3x_S_Dxy_S_C2001002+ABX*I_ERI_D2x_S_Dxy_S_C2001002;
  abcd[19] = I_ERI_F2xy_S_Dxy_S_C2001002+ABX*I_ERI_Dxy_S_Dxy_S_C2001002;
  abcd[20] = I_ERI_F2xz_S_Dxy_S_C2001002+ABX*I_ERI_Dxz_S_Dxy_S_C2001002;
  abcd[21] = I_ERI_Fx2y_S_Dxy_S_C2001002+ABX*I_ERI_D2y_S_Dxy_S_C2001002;
  abcd[22] = I_ERI_Fxyz_S_Dxy_S_C2001002+ABX*I_ERI_Dyz_S_Dxy_S_C2001002;
  abcd[23] = I_ERI_Fx2z_S_Dxy_S_C2001002+ABX*I_ERI_D2z_S_Dxy_S_C2001002;
  abcd[24] = I_ERI_F2xy_S_Dxy_S_C2001002+ABY*I_ERI_D2x_S_Dxy_S_C2001002;
  abcd[25] = I_ERI_Fx2y_S_Dxy_S_C2001002+ABY*I_ERI_Dxy_S_Dxy_S_C2001002;
  abcd[26] = I_ERI_Fxyz_S_Dxy_S_C2001002+ABY*I_ERI_Dxz_S_Dxy_S_C2001002;
  abcd[27] = I_ERI_F3y_S_Dxy_S_C2001002+ABY*I_ERI_D2y_S_Dxy_S_C2001002;
  abcd[28] = I_ERI_F2yz_S_Dxy_S_C2001002+ABY*I_ERI_Dyz_S_Dxy_S_C2001002;
  abcd[29] = I_ERI_Fy2z_S_Dxy_S_C2001002+ABY*I_ERI_D2z_S_Dxy_S_C2001002;
  abcd[30] = I_ERI_F2xz_S_Dxy_S_C2001002+ABZ*I_ERI_D2x_S_Dxy_S_C2001002;
  abcd[31] = I_ERI_Fxyz_S_Dxy_S_C2001002+ABZ*I_ERI_Dxy_S_Dxy_S_C2001002;
  abcd[32] = I_ERI_Fx2z_S_Dxy_S_C2001002+ABZ*I_ERI_Dxz_S_Dxy_S_C2001002;
  abcd[33] = I_ERI_F2yz_S_Dxy_S_C2001002+ABZ*I_ERI_D2y_S_Dxy_S_C2001002;
  abcd[34] = I_ERI_Fy2z_S_Dxy_S_C2001002+ABZ*I_ERI_Dyz_S_Dxy_S_C2001002;
  abcd[35] = I_ERI_F3z_S_Dxy_S_C2001002+ABZ*I_ERI_D2z_S_Dxy_S_C2001002;
  abcd[36] = I_ERI_F3x_S_Dxz_S_C2001002+ABX*I_ERI_D2x_S_Dxz_S_C2001002;
  abcd[37] = I_ERI_F2xy_S_Dxz_S_C2001002+ABX*I_ERI_Dxy_S_Dxz_S_C2001002;
  abcd[38] = I_ERI_F2xz_S_Dxz_S_C2001002+ABX*I_ERI_Dxz_S_Dxz_S_C2001002;
  abcd[39] = I_ERI_Fx2y_S_Dxz_S_C2001002+ABX*I_ERI_D2y_S_Dxz_S_C2001002;
  abcd[40] = I_ERI_Fxyz_S_Dxz_S_C2001002+ABX*I_ERI_Dyz_S_Dxz_S_C2001002;
  abcd[41] = I_ERI_Fx2z_S_Dxz_S_C2001002+ABX*I_ERI_D2z_S_Dxz_S_C2001002;
  abcd[42] = I_ERI_F2xy_S_Dxz_S_C2001002+ABY*I_ERI_D2x_S_Dxz_S_C2001002;
  abcd[43] = I_ERI_Fx2y_S_Dxz_S_C2001002+ABY*I_ERI_Dxy_S_Dxz_S_C2001002;
  abcd[44] = I_ERI_Fxyz_S_Dxz_S_C2001002+ABY*I_ERI_Dxz_S_Dxz_S_C2001002;
  abcd[45] = I_ERI_F3y_S_Dxz_S_C2001002+ABY*I_ERI_D2y_S_Dxz_S_C2001002;
  abcd[46] = I_ERI_F2yz_S_Dxz_S_C2001002+ABY*I_ERI_Dyz_S_Dxz_S_C2001002;
  abcd[47] = I_ERI_Fy2z_S_Dxz_S_C2001002+ABY*I_ERI_D2z_S_Dxz_S_C2001002;
  abcd[48] = I_ERI_F2xz_S_Dxz_S_C2001002+ABZ*I_ERI_D2x_S_Dxz_S_C2001002;
  abcd[49] = I_ERI_Fxyz_S_Dxz_S_C2001002+ABZ*I_ERI_Dxy_S_Dxz_S_C2001002;
  abcd[50] = I_ERI_Fx2z_S_Dxz_S_C2001002+ABZ*I_ERI_Dxz_S_Dxz_S_C2001002;
  abcd[51] = I_ERI_F2yz_S_Dxz_S_C2001002+ABZ*I_ERI_D2y_S_Dxz_S_C2001002;
  abcd[52] = I_ERI_Fy2z_S_Dxz_S_C2001002+ABZ*I_ERI_Dyz_S_Dxz_S_C2001002;
  abcd[53] = I_ERI_F3z_S_Dxz_S_C2001002+ABZ*I_ERI_D2z_S_Dxz_S_C2001002;
  abcd[54] = I_ERI_F3x_S_D2y_S_C2001002+ABX*I_ERI_D2x_S_D2y_S_C2001002;
  abcd[55] = I_ERI_F2xy_S_D2y_S_C2001002+ABX*I_ERI_Dxy_S_D2y_S_C2001002;
  abcd[56] = I_ERI_F2xz_S_D2y_S_C2001002+ABX*I_ERI_Dxz_S_D2y_S_C2001002;
  abcd[57] = I_ERI_Fx2y_S_D2y_S_C2001002+ABX*I_ERI_D2y_S_D2y_S_C2001002;
  abcd[58] = I_ERI_Fxyz_S_D2y_S_C2001002+ABX*I_ERI_Dyz_S_D2y_S_C2001002;
  abcd[59] = I_ERI_Fx2z_S_D2y_S_C2001002+ABX*I_ERI_D2z_S_D2y_S_C2001002;
  abcd[60] = I_ERI_F2xy_S_D2y_S_C2001002+ABY*I_ERI_D2x_S_D2y_S_C2001002;
  abcd[61] = I_ERI_Fx2y_S_D2y_S_C2001002+ABY*I_ERI_Dxy_S_D2y_S_C2001002;
  abcd[62] = I_ERI_Fxyz_S_D2y_S_C2001002+ABY*I_ERI_Dxz_S_D2y_S_C2001002;
  abcd[63] = I_ERI_F3y_S_D2y_S_C2001002+ABY*I_ERI_D2y_S_D2y_S_C2001002;
  abcd[64] = I_ERI_F2yz_S_D2y_S_C2001002+ABY*I_ERI_Dyz_S_D2y_S_C2001002;
  abcd[65] = I_ERI_Fy2z_S_D2y_S_C2001002+ABY*I_ERI_D2z_S_D2y_S_C2001002;
  abcd[66] = I_ERI_F2xz_S_D2y_S_C2001002+ABZ*I_ERI_D2x_S_D2y_S_C2001002;
  abcd[67] = I_ERI_Fxyz_S_D2y_S_C2001002+ABZ*I_ERI_Dxy_S_D2y_S_C2001002;
  abcd[68] = I_ERI_Fx2z_S_D2y_S_C2001002+ABZ*I_ERI_Dxz_S_D2y_S_C2001002;
  abcd[69] = I_ERI_F2yz_S_D2y_S_C2001002+ABZ*I_ERI_D2y_S_D2y_S_C2001002;
  abcd[70] = I_ERI_Fy2z_S_D2y_S_C2001002+ABZ*I_ERI_Dyz_S_D2y_S_C2001002;
  abcd[71] = I_ERI_F3z_S_D2y_S_C2001002+ABZ*I_ERI_D2z_S_D2y_S_C2001002;
  abcd[72] = I_ERI_F3x_S_Dyz_S_C2001002+ABX*I_ERI_D2x_S_Dyz_S_C2001002;
  abcd[73] = I_ERI_F2xy_S_Dyz_S_C2001002+ABX*I_ERI_Dxy_S_Dyz_S_C2001002;
  abcd[74] = I_ERI_F2xz_S_Dyz_S_C2001002+ABX*I_ERI_Dxz_S_Dyz_S_C2001002;
  abcd[75] = I_ERI_Fx2y_S_Dyz_S_C2001002+ABX*I_ERI_D2y_S_Dyz_S_C2001002;
  abcd[76] = I_ERI_Fxyz_S_Dyz_S_C2001002+ABX*I_ERI_Dyz_S_Dyz_S_C2001002;
  abcd[77] = I_ERI_Fx2z_S_Dyz_S_C2001002+ABX*I_ERI_D2z_S_Dyz_S_C2001002;
  abcd[78] = I_ERI_F2xy_S_Dyz_S_C2001002+ABY*I_ERI_D2x_S_Dyz_S_C2001002;
  abcd[79] = I_ERI_Fx2y_S_Dyz_S_C2001002+ABY*I_ERI_Dxy_S_Dyz_S_C2001002;
  abcd[80] = I_ERI_Fxyz_S_Dyz_S_C2001002+ABY*I_ERI_Dxz_S_Dyz_S_C2001002;
  abcd[81] = I_ERI_F3y_S_Dyz_S_C2001002+ABY*I_ERI_D2y_S_Dyz_S_C2001002;
  abcd[82] = I_ERI_F2yz_S_Dyz_S_C2001002+ABY*I_ERI_Dyz_S_Dyz_S_C2001002;
  abcd[83] = I_ERI_Fy2z_S_Dyz_S_C2001002+ABY*I_ERI_D2z_S_Dyz_S_C2001002;
  abcd[84] = I_ERI_F2xz_S_Dyz_S_C2001002+ABZ*I_ERI_D2x_S_Dyz_S_C2001002;
  abcd[85] = I_ERI_Fxyz_S_Dyz_S_C2001002+ABZ*I_ERI_Dxy_S_Dyz_S_C2001002;
  abcd[86] = I_ERI_Fx2z_S_Dyz_S_C2001002+ABZ*I_ERI_Dxz_S_Dyz_S_C2001002;
  abcd[87] = I_ERI_F2yz_S_Dyz_S_C2001002+ABZ*I_ERI_D2y_S_Dyz_S_C2001002;
  abcd[88] = I_ERI_Fy2z_S_Dyz_S_C2001002+ABZ*I_ERI_Dyz_S_Dyz_S_C2001002;
  abcd[89] = I_ERI_F3z_S_Dyz_S_C2001002+ABZ*I_ERI_D2z_S_Dyz_S_C2001002;
  abcd[90] = I_ERI_F3x_S_D2z_S_C2001002+ABX*I_ERI_D2x_S_D2z_S_C2001002;
  abcd[91] = I_ERI_F2xy_S_D2z_S_C2001002+ABX*I_ERI_Dxy_S_D2z_S_C2001002;
  abcd[92] = I_ERI_F2xz_S_D2z_S_C2001002+ABX*I_ERI_Dxz_S_D2z_S_C2001002;
  abcd[93] = I_ERI_Fx2y_S_D2z_S_C2001002+ABX*I_ERI_D2y_S_D2z_S_C2001002;
  abcd[94] = I_ERI_Fxyz_S_D2z_S_C2001002+ABX*I_ERI_Dyz_S_D2z_S_C2001002;
  abcd[95] = I_ERI_Fx2z_S_D2z_S_C2001002+ABX*I_ERI_D2z_S_D2z_S_C2001002;
  abcd[96] = I_ERI_F2xy_S_D2z_S_C2001002+ABY*I_ERI_D2x_S_D2z_S_C2001002;
  abcd[97] = I_ERI_Fx2y_S_D2z_S_C2001002+ABY*I_ERI_Dxy_S_D2z_S_C2001002;
  abcd[98] = I_ERI_Fxyz_S_D2z_S_C2001002+ABY*I_ERI_Dxz_S_D2z_S_C2001002;
  abcd[99] = I_ERI_F3y_S_D2z_S_C2001002+ABY*I_ERI_D2y_S_D2z_S_C2001002;
  abcd[100] = I_ERI_F2yz_S_D2z_S_C2001002+ABY*I_ERI_Dyz_S_D2z_S_C2001002;
  abcd[101] = I_ERI_Fy2z_S_D2z_S_C2001002+ABY*I_ERI_D2z_S_D2z_S_C2001002;
  abcd[102] = I_ERI_F2xz_S_D2z_S_C2001002+ABZ*I_ERI_D2x_S_D2z_S_C2001002;
  abcd[103] = I_ERI_Fxyz_S_D2z_S_C2001002+ABZ*I_ERI_Dxy_S_D2z_S_C2001002;
  abcd[104] = I_ERI_Fx2z_S_D2z_S_C2001002+ABZ*I_ERI_Dxz_S_D2z_S_C2001002;
  abcd[105] = I_ERI_F2yz_S_D2z_S_C2001002+ABZ*I_ERI_D2y_S_D2z_S_C2001002;
  abcd[106] = I_ERI_Fy2z_S_D2z_S_C2001002+ABZ*I_ERI_Dyz_S_D2z_S_C2001002;
  abcd[107] = I_ERI_F3z_S_D2z_S_C2001002+ABZ*I_ERI_D2z_S_D2z_S_C2001002;

  /************************************************************
   * shell quartet name: SQ_ERI_D_P_D_P_C1002001002
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_P_C1002001002
   * RHS shell quartet name: SQ_ERI_D_S_D_P_C1002001002
   ************************************************************/
  abcd[108] = I_ERI_F3x_S_D2x_Px_C1002001002+ABX*I_ERI_D2x_S_D2x_Px_C1002001002;
  abcd[109] = I_ERI_F2xy_S_D2x_Px_C1002001002+ABX*I_ERI_Dxy_S_D2x_Px_C1002001002;
  abcd[110] = I_ERI_F2xz_S_D2x_Px_C1002001002+ABX*I_ERI_Dxz_S_D2x_Px_C1002001002;
  abcd[111] = I_ERI_Fx2y_S_D2x_Px_C1002001002+ABX*I_ERI_D2y_S_D2x_Px_C1002001002;
  abcd[112] = I_ERI_Fxyz_S_D2x_Px_C1002001002+ABX*I_ERI_Dyz_S_D2x_Px_C1002001002;
  abcd[113] = I_ERI_Fx2z_S_D2x_Px_C1002001002+ABX*I_ERI_D2z_S_D2x_Px_C1002001002;
  abcd[114] = I_ERI_F2xy_S_D2x_Px_C1002001002+ABY*I_ERI_D2x_S_D2x_Px_C1002001002;
  abcd[115] = I_ERI_Fx2y_S_D2x_Px_C1002001002+ABY*I_ERI_Dxy_S_D2x_Px_C1002001002;
  abcd[116] = I_ERI_Fxyz_S_D2x_Px_C1002001002+ABY*I_ERI_Dxz_S_D2x_Px_C1002001002;
  abcd[117] = I_ERI_F3y_S_D2x_Px_C1002001002+ABY*I_ERI_D2y_S_D2x_Px_C1002001002;
  abcd[118] = I_ERI_F2yz_S_D2x_Px_C1002001002+ABY*I_ERI_Dyz_S_D2x_Px_C1002001002;
  abcd[119] = I_ERI_Fy2z_S_D2x_Px_C1002001002+ABY*I_ERI_D2z_S_D2x_Px_C1002001002;
  abcd[120] = I_ERI_F2xz_S_D2x_Px_C1002001002+ABZ*I_ERI_D2x_S_D2x_Px_C1002001002;
  abcd[121] = I_ERI_Fxyz_S_D2x_Px_C1002001002+ABZ*I_ERI_Dxy_S_D2x_Px_C1002001002;
  abcd[122] = I_ERI_Fx2z_S_D2x_Px_C1002001002+ABZ*I_ERI_Dxz_S_D2x_Px_C1002001002;
  abcd[123] = I_ERI_F2yz_S_D2x_Px_C1002001002+ABZ*I_ERI_D2y_S_D2x_Px_C1002001002;
  abcd[124] = I_ERI_Fy2z_S_D2x_Px_C1002001002+ABZ*I_ERI_Dyz_S_D2x_Px_C1002001002;
  abcd[125] = I_ERI_F3z_S_D2x_Px_C1002001002+ABZ*I_ERI_D2z_S_D2x_Px_C1002001002;
  abcd[126] = I_ERI_F3x_S_Dxy_Px_C1002001002+ABX*I_ERI_D2x_S_Dxy_Px_C1002001002;
  abcd[127] = I_ERI_F2xy_S_Dxy_Px_C1002001002+ABX*I_ERI_Dxy_S_Dxy_Px_C1002001002;
  abcd[128] = I_ERI_F2xz_S_Dxy_Px_C1002001002+ABX*I_ERI_Dxz_S_Dxy_Px_C1002001002;
  abcd[129] = I_ERI_Fx2y_S_Dxy_Px_C1002001002+ABX*I_ERI_D2y_S_Dxy_Px_C1002001002;
  abcd[130] = I_ERI_Fxyz_S_Dxy_Px_C1002001002+ABX*I_ERI_Dyz_S_Dxy_Px_C1002001002;
  abcd[131] = I_ERI_Fx2z_S_Dxy_Px_C1002001002+ABX*I_ERI_D2z_S_Dxy_Px_C1002001002;
  abcd[132] = I_ERI_F2xy_S_Dxy_Px_C1002001002+ABY*I_ERI_D2x_S_Dxy_Px_C1002001002;
  abcd[133] = I_ERI_Fx2y_S_Dxy_Px_C1002001002+ABY*I_ERI_Dxy_S_Dxy_Px_C1002001002;
  abcd[134] = I_ERI_Fxyz_S_Dxy_Px_C1002001002+ABY*I_ERI_Dxz_S_Dxy_Px_C1002001002;
  abcd[135] = I_ERI_F3y_S_Dxy_Px_C1002001002+ABY*I_ERI_D2y_S_Dxy_Px_C1002001002;
  abcd[136] = I_ERI_F2yz_S_Dxy_Px_C1002001002+ABY*I_ERI_Dyz_S_Dxy_Px_C1002001002;
  abcd[137] = I_ERI_Fy2z_S_Dxy_Px_C1002001002+ABY*I_ERI_D2z_S_Dxy_Px_C1002001002;
  abcd[138] = I_ERI_F2xz_S_Dxy_Px_C1002001002+ABZ*I_ERI_D2x_S_Dxy_Px_C1002001002;
  abcd[139] = I_ERI_Fxyz_S_Dxy_Px_C1002001002+ABZ*I_ERI_Dxy_S_Dxy_Px_C1002001002;
  abcd[140] = I_ERI_Fx2z_S_Dxy_Px_C1002001002+ABZ*I_ERI_Dxz_S_Dxy_Px_C1002001002;
  abcd[141] = I_ERI_F2yz_S_Dxy_Px_C1002001002+ABZ*I_ERI_D2y_S_Dxy_Px_C1002001002;
  abcd[142] = I_ERI_Fy2z_S_Dxy_Px_C1002001002+ABZ*I_ERI_Dyz_S_Dxy_Px_C1002001002;
  abcd[143] = I_ERI_F3z_S_Dxy_Px_C1002001002+ABZ*I_ERI_D2z_S_Dxy_Px_C1002001002;
  abcd[144] = I_ERI_F3x_S_Dxz_Px_C1002001002+ABX*I_ERI_D2x_S_Dxz_Px_C1002001002;
  abcd[145] = I_ERI_F2xy_S_Dxz_Px_C1002001002+ABX*I_ERI_Dxy_S_Dxz_Px_C1002001002;
  abcd[146] = I_ERI_F2xz_S_Dxz_Px_C1002001002+ABX*I_ERI_Dxz_S_Dxz_Px_C1002001002;
  abcd[147] = I_ERI_Fx2y_S_Dxz_Px_C1002001002+ABX*I_ERI_D2y_S_Dxz_Px_C1002001002;
  abcd[148] = I_ERI_Fxyz_S_Dxz_Px_C1002001002+ABX*I_ERI_Dyz_S_Dxz_Px_C1002001002;
  abcd[149] = I_ERI_Fx2z_S_Dxz_Px_C1002001002+ABX*I_ERI_D2z_S_Dxz_Px_C1002001002;
  abcd[150] = I_ERI_F2xy_S_Dxz_Px_C1002001002+ABY*I_ERI_D2x_S_Dxz_Px_C1002001002;
  abcd[151] = I_ERI_Fx2y_S_Dxz_Px_C1002001002+ABY*I_ERI_Dxy_S_Dxz_Px_C1002001002;
  abcd[152] = I_ERI_Fxyz_S_Dxz_Px_C1002001002+ABY*I_ERI_Dxz_S_Dxz_Px_C1002001002;
  abcd[153] = I_ERI_F3y_S_Dxz_Px_C1002001002+ABY*I_ERI_D2y_S_Dxz_Px_C1002001002;
  abcd[154] = I_ERI_F2yz_S_Dxz_Px_C1002001002+ABY*I_ERI_Dyz_S_Dxz_Px_C1002001002;
  abcd[155] = I_ERI_Fy2z_S_Dxz_Px_C1002001002+ABY*I_ERI_D2z_S_Dxz_Px_C1002001002;
  abcd[156] = I_ERI_F2xz_S_Dxz_Px_C1002001002+ABZ*I_ERI_D2x_S_Dxz_Px_C1002001002;
  abcd[157] = I_ERI_Fxyz_S_Dxz_Px_C1002001002+ABZ*I_ERI_Dxy_S_Dxz_Px_C1002001002;
  abcd[158] = I_ERI_Fx2z_S_Dxz_Px_C1002001002+ABZ*I_ERI_Dxz_S_Dxz_Px_C1002001002;
  abcd[159] = I_ERI_F2yz_S_Dxz_Px_C1002001002+ABZ*I_ERI_D2y_S_Dxz_Px_C1002001002;
  abcd[160] = I_ERI_Fy2z_S_Dxz_Px_C1002001002+ABZ*I_ERI_Dyz_S_Dxz_Px_C1002001002;
  abcd[161] = I_ERI_F3z_S_Dxz_Px_C1002001002+ABZ*I_ERI_D2z_S_Dxz_Px_C1002001002;
  abcd[162] = I_ERI_F3x_S_D2y_Px_C1002001002+ABX*I_ERI_D2x_S_D2y_Px_C1002001002;
  abcd[163] = I_ERI_F2xy_S_D2y_Px_C1002001002+ABX*I_ERI_Dxy_S_D2y_Px_C1002001002;
  abcd[164] = I_ERI_F2xz_S_D2y_Px_C1002001002+ABX*I_ERI_Dxz_S_D2y_Px_C1002001002;
  abcd[165] = I_ERI_Fx2y_S_D2y_Px_C1002001002+ABX*I_ERI_D2y_S_D2y_Px_C1002001002;
  abcd[166] = I_ERI_Fxyz_S_D2y_Px_C1002001002+ABX*I_ERI_Dyz_S_D2y_Px_C1002001002;
  abcd[167] = I_ERI_Fx2z_S_D2y_Px_C1002001002+ABX*I_ERI_D2z_S_D2y_Px_C1002001002;
  abcd[168] = I_ERI_F2xy_S_D2y_Px_C1002001002+ABY*I_ERI_D2x_S_D2y_Px_C1002001002;
  abcd[169] = I_ERI_Fx2y_S_D2y_Px_C1002001002+ABY*I_ERI_Dxy_S_D2y_Px_C1002001002;
  abcd[170] = I_ERI_Fxyz_S_D2y_Px_C1002001002+ABY*I_ERI_Dxz_S_D2y_Px_C1002001002;
  abcd[171] = I_ERI_F3y_S_D2y_Px_C1002001002+ABY*I_ERI_D2y_S_D2y_Px_C1002001002;
  abcd[172] = I_ERI_F2yz_S_D2y_Px_C1002001002+ABY*I_ERI_Dyz_S_D2y_Px_C1002001002;
  abcd[173] = I_ERI_Fy2z_S_D2y_Px_C1002001002+ABY*I_ERI_D2z_S_D2y_Px_C1002001002;
  abcd[174] = I_ERI_F2xz_S_D2y_Px_C1002001002+ABZ*I_ERI_D2x_S_D2y_Px_C1002001002;
  abcd[175] = I_ERI_Fxyz_S_D2y_Px_C1002001002+ABZ*I_ERI_Dxy_S_D2y_Px_C1002001002;
  abcd[176] = I_ERI_Fx2z_S_D2y_Px_C1002001002+ABZ*I_ERI_Dxz_S_D2y_Px_C1002001002;
  abcd[177] = I_ERI_F2yz_S_D2y_Px_C1002001002+ABZ*I_ERI_D2y_S_D2y_Px_C1002001002;
  abcd[178] = I_ERI_Fy2z_S_D2y_Px_C1002001002+ABZ*I_ERI_Dyz_S_D2y_Px_C1002001002;
  abcd[179] = I_ERI_F3z_S_D2y_Px_C1002001002+ABZ*I_ERI_D2z_S_D2y_Px_C1002001002;
  abcd[180] = I_ERI_F3x_S_Dyz_Px_C1002001002+ABX*I_ERI_D2x_S_Dyz_Px_C1002001002;
  abcd[181] = I_ERI_F2xy_S_Dyz_Px_C1002001002+ABX*I_ERI_Dxy_S_Dyz_Px_C1002001002;
  abcd[182] = I_ERI_F2xz_S_Dyz_Px_C1002001002+ABX*I_ERI_Dxz_S_Dyz_Px_C1002001002;
  abcd[183] = I_ERI_Fx2y_S_Dyz_Px_C1002001002+ABX*I_ERI_D2y_S_Dyz_Px_C1002001002;
  abcd[184] = I_ERI_Fxyz_S_Dyz_Px_C1002001002+ABX*I_ERI_Dyz_S_Dyz_Px_C1002001002;
  abcd[185] = I_ERI_Fx2z_S_Dyz_Px_C1002001002+ABX*I_ERI_D2z_S_Dyz_Px_C1002001002;
  abcd[186] = I_ERI_F2xy_S_Dyz_Px_C1002001002+ABY*I_ERI_D2x_S_Dyz_Px_C1002001002;
  abcd[187] = I_ERI_Fx2y_S_Dyz_Px_C1002001002+ABY*I_ERI_Dxy_S_Dyz_Px_C1002001002;
  abcd[188] = I_ERI_Fxyz_S_Dyz_Px_C1002001002+ABY*I_ERI_Dxz_S_Dyz_Px_C1002001002;
  abcd[189] = I_ERI_F3y_S_Dyz_Px_C1002001002+ABY*I_ERI_D2y_S_Dyz_Px_C1002001002;
  abcd[190] = I_ERI_F2yz_S_Dyz_Px_C1002001002+ABY*I_ERI_Dyz_S_Dyz_Px_C1002001002;
  abcd[191] = I_ERI_Fy2z_S_Dyz_Px_C1002001002+ABY*I_ERI_D2z_S_Dyz_Px_C1002001002;
  abcd[192] = I_ERI_F2xz_S_Dyz_Px_C1002001002+ABZ*I_ERI_D2x_S_Dyz_Px_C1002001002;
  abcd[193] = I_ERI_Fxyz_S_Dyz_Px_C1002001002+ABZ*I_ERI_Dxy_S_Dyz_Px_C1002001002;
  abcd[194] = I_ERI_Fx2z_S_Dyz_Px_C1002001002+ABZ*I_ERI_Dxz_S_Dyz_Px_C1002001002;
  abcd[195] = I_ERI_F2yz_S_Dyz_Px_C1002001002+ABZ*I_ERI_D2y_S_Dyz_Px_C1002001002;
  abcd[196] = I_ERI_Fy2z_S_Dyz_Px_C1002001002+ABZ*I_ERI_Dyz_S_Dyz_Px_C1002001002;
  abcd[197] = I_ERI_F3z_S_Dyz_Px_C1002001002+ABZ*I_ERI_D2z_S_Dyz_Px_C1002001002;
  abcd[198] = I_ERI_F3x_S_D2z_Px_C1002001002+ABX*I_ERI_D2x_S_D2z_Px_C1002001002;
  abcd[199] = I_ERI_F2xy_S_D2z_Px_C1002001002+ABX*I_ERI_Dxy_S_D2z_Px_C1002001002;
  abcd[200] = I_ERI_F2xz_S_D2z_Px_C1002001002+ABX*I_ERI_Dxz_S_D2z_Px_C1002001002;
  abcd[201] = I_ERI_Fx2y_S_D2z_Px_C1002001002+ABX*I_ERI_D2y_S_D2z_Px_C1002001002;
  abcd[202] = I_ERI_Fxyz_S_D2z_Px_C1002001002+ABX*I_ERI_Dyz_S_D2z_Px_C1002001002;
  abcd[203] = I_ERI_Fx2z_S_D2z_Px_C1002001002+ABX*I_ERI_D2z_S_D2z_Px_C1002001002;
  abcd[204] = I_ERI_F2xy_S_D2z_Px_C1002001002+ABY*I_ERI_D2x_S_D2z_Px_C1002001002;
  abcd[205] = I_ERI_Fx2y_S_D2z_Px_C1002001002+ABY*I_ERI_Dxy_S_D2z_Px_C1002001002;
  abcd[206] = I_ERI_Fxyz_S_D2z_Px_C1002001002+ABY*I_ERI_Dxz_S_D2z_Px_C1002001002;
  abcd[207] = I_ERI_F3y_S_D2z_Px_C1002001002+ABY*I_ERI_D2y_S_D2z_Px_C1002001002;
  abcd[208] = I_ERI_F2yz_S_D2z_Px_C1002001002+ABY*I_ERI_Dyz_S_D2z_Px_C1002001002;
  abcd[209] = I_ERI_Fy2z_S_D2z_Px_C1002001002+ABY*I_ERI_D2z_S_D2z_Px_C1002001002;
  abcd[210] = I_ERI_F2xz_S_D2z_Px_C1002001002+ABZ*I_ERI_D2x_S_D2z_Px_C1002001002;
  abcd[211] = I_ERI_Fxyz_S_D2z_Px_C1002001002+ABZ*I_ERI_Dxy_S_D2z_Px_C1002001002;
  abcd[212] = I_ERI_Fx2z_S_D2z_Px_C1002001002+ABZ*I_ERI_Dxz_S_D2z_Px_C1002001002;
  abcd[213] = I_ERI_F2yz_S_D2z_Px_C1002001002+ABZ*I_ERI_D2y_S_D2z_Px_C1002001002;
  abcd[214] = I_ERI_Fy2z_S_D2z_Px_C1002001002+ABZ*I_ERI_Dyz_S_D2z_Px_C1002001002;
  abcd[215] = I_ERI_F3z_S_D2z_Px_C1002001002+ABZ*I_ERI_D2z_S_D2z_Px_C1002001002;
  abcd[216] = I_ERI_F3x_S_D2x_Py_C1002001002+ABX*I_ERI_D2x_S_D2x_Py_C1002001002;
  abcd[217] = I_ERI_F2xy_S_D2x_Py_C1002001002+ABX*I_ERI_Dxy_S_D2x_Py_C1002001002;
  abcd[218] = I_ERI_F2xz_S_D2x_Py_C1002001002+ABX*I_ERI_Dxz_S_D2x_Py_C1002001002;
  abcd[219] = I_ERI_Fx2y_S_D2x_Py_C1002001002+ABX*I_ERI_D2y_S_D2x_Py_C1002001002;
  abcd[220] = I_ERI_Fxyz_S_D2x_Py_C1002001002+ABX*I_ERI_Dyz_S_D2x_Py_C1002001002;
  abcd[221] = I_ERI_Fx2z_S_D2x_Py_C1002001002+ABX*I_ERI_D2z_S_D2x_Py_C1002001002;
  abcd[222] = I_ERI_F2xy_S_D2x_Py_C1002001002+ABY*I_ERI_D2x_S_D2x_Py_C1002001002;
  abcd[223] = I_ERI_Fx2y_S_D2x_Py_C1002001002+ABY*I_ERI_Dxy_S_D2x_Py_C1002001002;
  abcd[224] = I_ERI_Fxyz_S_D2x_Py_C1002001002+ABY*I_ERI_Dxz_S_D2x_Py_C1002001002;
  abcd[225] = I_ERI_F3y_S_D2x_Py_C1002001002+ABY*I_ERI_D2y_S_D2x_Py_C1002001002;
  abcd[226] = I_ERI_F2yz_S_D2x_Py_C1002001002+ABY*I_ERI_Dyz_S_D2x_Py_C1002001002;
  abcd[227] = I_ERI_Fy2z_S_D2x_Py_C1002001002+ABY*I_ERI_D2z_S_D2x_Py_C1002001002;
  abcd[228] = I_ERI_F2xz_S_D2x_Py_C1002001002+ABZ*I_ERI_D2x_S_D2x_Py_C1002001002;
  abcd[229] = I_ERI_Fxyz_S_D2x_Py_C1002001002+ABZ*I_ERI_Dxy_S_D2x_Py_C1002001002;
  abcd[230] = I_ERI_Fx2z_S_D2x_Py_C1002001002+ABZ*I_ERI_Dxz_S_D2x_Py_C1002001002;
  abcd[231] = I_ERI_F2yz_S_D2x_Py_C1002001002+ABZ*I_ERI_D2y_S_D2x_Py_C1002001002;
  abcd[232] = I_ERI_Fy2z_S_D2x_Py_C1002001002+ABZ*I_ERI_Dyz_S_D2x_Py_C1002001002;
  abcd[233] = I_ERI_F3z_S_D2x_Py_C1002001002+ABZ*I_ERI_D2z_S_D2x_Py_C1002001002;
  abcd[234] = I_ERI_F3x_S_Dxy_Py_C1002001002+ABX*I_ERI_D2x_S_Dxy_Py_C1002001002;
  abcd[235] = I_ERI_F2xy_S_Dxy_Py_C1002001002+ABX*I_ERI_Dxy_S_Dxy_Py_C1002001002;
  abcd[236] = I_ERI_F2xz_S_Dxy_Py_C1002001002+ABX*I_ERI_Dxz_S_Dxy_Py_C1002001002;
  abcd[237] = I_ERI_Fx2y_S_Dxy_Py_C1002001002+ABX*I_ERI_D2y_S_Dxy_Py_C1002001002;
  abcd[238] = I_ERI_Fxyz_S_Dxy_Py_C1002001002+ABX*I_ERI_Dyz_S_Dxy_Py_C1002001002;
  abcd[239] = I_ERI_Fx2z_S_Dxy_Py_C1002001002+ABX*I_ERI_D2z_S_Dxy_Py_C1002001002;
  abcd[240] = I_ERI_F2xy_S_Dxy_Py_C1002001002+ABY*I_ERI_D2x_S_Dxy_Py_C1002001002;
  abcd[241] = I_ERI_Fx2y_S_Dxy_Py_C1002001002+ABY*I_ERI_Dxy_S_Dxy_Py_C1002001002;
  abcd[242] = I_ERI_Fxyz_S_Dxy_Py_C1002001002+ABY*I_ERI_Dxz_S_Dxy_Py_C1002001002;
  abcd[243] = I_ERI_F3y_S_Dxy_Py_C1002001002+ABY*I_ERI_D2y_S_Dxy_Py_C1002001002;
  abcd[244] = I_ERI_F2yz_S_Dxy_Py_C1002001002+ABY*I_ERI_Dyz_S_Dxy_Py_C1002001002;
  abcd[245] = I_ERI_Fy2z_S_Dxy_Py_C1002001002+ABY*I_ERI_D2z_S_Dxy_Py_C1002001002;
  abcd[246] = I_ERI_F2xz_S_Dxy_Py_C1002001002+ABZ*I_ERI_D2x_S_Dxy_Py_C1002001002;
  abcd[247] = I_ERI_Fxyz_S_Dxy_Py_C1002001002+ABZ*I_ERI_Dxy_S_Dxy_Py_C1002001002;
  abcd[248] = I_ERI_Fx2z_S_Dxy_Py_C1002001002+ABZ*I_ERI_Dxz_S_Dxy_Py_C1002001002;
  abcd[249] = I_ERI_F2yz_S_Dxy_Py_C1002001002+ABZ*I_ERI_D2y_S_Dxy_Py_C1002001002;
  abcd[250] = I_ERI_Fy2z_S_Dxy_Py_C1002001002+ABZ*I_ERI_Dyz_S_Dxy_Py_C1002001002;
  abcd[251] = I_ERI_F3z_S_Dxy_Py_C1002001002+ABZ*I_ERI_D2z_S_Dxy_Py_C1002001002;
  abcd[252] = I_ERI_F3x_S_Dxz_Py_C1002001002+ABX*I_ERI_D2x_S_Dxz_Py_C1002001002;
  abcd[253] = I_ERI_F2xy_S_Dxz_Py_C1002001002+ABX*I_ERI_Dxy_S_Dxz_Py_C1002001002;
  abcd[254] = I_ERI_F2xz_S_Dxz_Py_C1002001002+ABX*I_ERI_Dxz_S_Dxz_Py_C1002001002;
  abcd[255] = I_ERI_Fx2y_S_Dxz_Py_C1002001002+ABX*I_ERI_D2y_S_Dxz_Py_C1002001002;
  abcd[256] = I_ERI_Fxyz_S_Dxz_Py_C1002001002+ABX*I_ERI_Dyz_S_Dxz_Py_C1002001002;
  abcd[257] = I_ERI_Fx2z_S_Dxz_Py_C1002001002+ABX*I_ERI_D2z_S_Dxz_Py_C1002001002;
  abcd[258] = I_ERI_F2xy_S_Dxz_Py_C1002001002+ABY*I_ERI_D2x_S_Dxz_Py_C1002001002;
  abcd[259] = I_ERI_Fx2y_S_Dxz_Py_C1002001002+ABY*I_ERI_Dxy_S_Dxz_Py_C1002001002;
  abcd[260] = I_ERI_Fxyz_S_Dxz_Py_C1002001002+ABY*I_ERI_Dxz_S_Dxz_Py_C1002001002;
  abcd[261] = I_ERI_F3y_S_Dxz_Py_C1002001002+ABY*I_ERI_D2y_S_Dxz_Py_C1002001002;
  abcd[262] = I_ERI_F2yz_S_Dxz_Py_C1002001002+ABY*I_ERI_Dyz_S_Dxz_Py_C1002001002;
  abcd[263] = I_ERI_Fy2z_S_Dxz_Py_C1002001002+ABY*I_ERI_D2z_S_Dxz_Py_C1002001002;
  abcd[264] = I_ERI_F2xz_S_Dxz_Py_C1002001002+ABZ*I_ERI_D2x_S_Dxz_Py_C1002001002;
  abcd[265] = I_ERI_Fxyz_S_Dxz_Py_C1002001002+ABZ*I_ERI_Dxy_S_Dxz_Py_C1002001002;
  abcd[266] = I_ERI_Fx2z_S_Dxz_Py_C1002001002+ABZ*I_ERI_Dxz_S_Dxz_Py_C1002001002;
  abcd[267] = I_ERI_F2yz_S_Dxz_Py_C1002001002+ABZ*I_ERI_D2y_S_Dxz_Py_C1002001002;
  abcd[268] = I_ERI_Fy2z_S_Dxz_Py_C1002001002+ABZ*I_ERI_Dyz_S_Dxz_Py_C1002001002;
  abcd[269] = I_ERI_F3z_S_Dxz_Py_C1002001002+ABZ*I_ERI_D2z_S_Dxz_Py_C1002001002;
  abcd[270] = I_ERI_F3x_S_D2y_Py_C1002001002+ABX*I_ERI_D2x_S_D2y_Py_C1002001002;
  abcd[271] = I_ERI_F2xy_S_D2y_Py_C1002001002+ABX*I_ERI_Dxy_S_D2y_Py_C1002001002;
  abcd[272] = I_ERI_F2xz_S_D2y_Py_C1002001002+ABX*I_ERI_Dxz_S_D2y_Py_C1002001002;
  abcd[273] = I_ERI_Fx2y_S_D2y_Py_C1002001002+ABX*I_ERI_D2y_S_D2y_Py_C1002001002;
  abcd[274] = I_ERI_Fxyz_S_D2y_Py_C1002001002+ABX*I_ERI_Dyz_S_D2y_Py_C1002001002;
  abcd[275] = I_ERI_Fx2z_S_D2y_Py_C1002001002+ABX*I_ERI_D2z_S_D2y_Py_C1002001002;
  abcd[276] = I_ERI_F2xy_S_D2y_Py_C1002001002+ABY*I_ERI_D2x_S_D2y_Py_C1002001002;
  abcd[277] = I_ERI_Fx2y_S_D2y_Py_C1002001002+ABY*I_ERI_Dxy_S_D2y_Py_C1002001002;
  abcd[278] = I_ERI_Fxyz_S_D2y_Py_C1002001002+ABY*I_ERI_Dxz_S_D2y_Py_C1002001002;
  abcd[279] = I_ERI_F3y_S_D2y_Py_C1002001002+ABY*I_ERI_D2y_S_D2y_Py_C1002001002;
  abcd[280] = I_ERI_F2yz_S_D2y_Py_C1002001002+ABY*I_ERI_Dyz_S_D2y_Py_C1002001002;
  abcd[281] = I_ERI_Fy2z_S_D2y_Py_C1002001002+ABY*I_ERI_D2z_S_D2y_Py_C1002001002;
  abcd[282] = I_ERI_F2xz_S_D2y_Py_C1002001002+ABZ*I_ERI_D2x_S_D2y_Py_C1002001002;
  abcd[283] = I_ERI_Fxyz_S_D2y_Py_C1002001002+ABZ*I_ERI_Dxy_S_D2y_Py_C1002001002;
  abcd[284] = I_ERI_Fx2z_S_D2y_Py_C1002001002+ABZ*I_ERI_Dxz_S_D2y_Py_C1002001002;
  abcd[285] = I_ERI_F2yz_S_D2y_Py_C1002001002+ABZ*I_ERI_D2y_S_D2y_Py_C1002001002;
  abcd[286] = I_ERI_Fy2z_S_D2y_Py_C1002001002+ABZ*I_ERI_Dyz_S_D2y_Py_C1002001002;
  abcd[287] = I_ERI_F3z_S_D2y_Py_C1002001002+ABZ*I_ERI_D2z_S_D2y_Py_C1002001002;
  abcd[288] = I_ERI_F3x_S_Dyz_Py_C1002001002+ABX*I_ERI_D2x_S_Dyz_Py_C1002001002;
  abcd[289] = I_ERI_F2xy_S_Dyz_Py_C1002001002+ABX*I_ERI_Dxy_S_Dyz_Py_C1002001002;
  abcd[290] = I_ERI_F2xz_S_Dyz_Py_C1002001002+ABX*I_ERI_Dxz_S_Dyz_Py_C1002001002;
  abcd[291] = I_ERI_Fx2y_S_Dyz_Py_C1002001002+ABX*I_ERI_D2y_S_Dyz_Py_C1002001002;
  abcd[292] = I_ERI_Fxyz_S_Dyz_Py_C1002001002+ABX*I_ERI_Dyz_S_Dyz_Py_C1002001002;
  abcd[293] = I_ERI_Fx2z_S_Dyz_Py_C1002001002+ABX*I_ERI_D2z_S_Dyz_Py_C1002001002;
  abcd[294] = I_ERI_F2xy_S_Dyz_Py_C1002001002+ABY*I_ERI_D2x_S_Dyz_Py_C1002001002;
  abcd[295] = I_ERI_Fx2y_S_Dyz_Py_C1002001002+ABY*I_ERI_Dxy_S_Dyz_Py_C1002001002;
  abcd[296] = I_ERI_Fxyz_S_Dyz_Py_C1002001002+ABY*I_ERI_Dxz_S_Dyz_Py_C1002001002;
  abcd[297] = I_ERI_F3y_S_Dyz_Py_C1002001002+ABY*I_ERI_D2y_S_Dyz_Py_C1002001002;
  abcd[298] = I_ERI_F2yz_S_Dyz_Py_C1002001002+ABY*I_ERI_Dyz_S_Dyz_Py_C1002001002;
  abcd[299] = I_ERI_Fy2z_S_Dyz_Py_C1002001002+ABY*I_ERI_D2z_S_Dyz_Py_C1002001002;
  abcd[300] = I_ERI_F2xz_S_Dyz_Py_C1002001002+ABZ*I_ERI_D2x_S_Dyz_Py_C1002001002;
  abcd[301] = I_ERI_Fxyz_S_Dyz_Py_C1002001002+ABZ*I_ERI_Dxy_S_Dyz_Py_C1002001002;
  abcd[302] = I_ERI_Fx2z_S_Dyz_Py_C1002001002+ABZ*I_ERI_Dxz_S_Dyz_Py_C1002001002;
  abcd[303] = I_ERI_F2yz_S_Dyz_Py_C1002001002+ABZ*I_ERI_D2y_S_Dyz_Py_C1002001002;
  abcd[304] = I_ERI_Fy2z_S_Dyz_Py_C1002001002+ABZ*I_ERI_Dyz_S_Dyz_Py_C1002001002;
  abcd[305] = I_ERI_F3z_S_Dyz_Py_C1002001002+ABZ*I_ERI_D2z_S_Dyz_Py_C1002001002;
  abcd[306] = I_ERI_F3x_S_D2z_Py_C1002001002+ABX*I_ERI_D2x_S_D2z_Py_C1002001002;
  abcd[307] = I_ERI_F2xy_S_D2z_Py_C1002001002+ABX*I_ERI_Dxy_S_D2z_Py_C1002001002;
  abcd[308] = I_ERI_F2xz_S_D2z_Py_C1002001002+ABX*I_ERI_Dxz_S_D2z_Py_C1002001002;
  abcd[309] = I_ERI_Fx2y_S_D2z_Py_C1002001002+ABX*I_ERI_D2y_S_D2z_Py_C1002001002;
  abcd[310] = I_ERI_Fxyz_S_D2z_Py_C1002001002+ABX*I_ERI_Dyz_S_D2z_Py_C1002001002;
  abcd[311] = I_ERI_Fx2z_S_D2z_Py_C1002001002+ABX*I_ERI_D2z_S_D2z_Py_C1002001002;
  abcd[312] = I_ERI_F2xy_S_D2z_Py_C1002001002+ABY*I_ERI_D2x_S_D2z_Py_C1002001002;
  abcd[313] = I_ERI_Fx2y_S_D2z_Py_C1002001002+ABY*I_ERI_Dxy_S_D2z_Py_C1002001002;
  abcd[314] = I_ERI_Fxyz_S_D2z_Py_C1002001002+ABY*I_ERI_Dxz_S_D2z_Py_C1002001002;
  abcd[315] = I_ERI_F3y_S_D2z_Py_C1002001002+ABY*I_ERI_D2y_S_D2z_Py_C1002001002;
  abcd[316] = I_ERI_F2yz_S_D2z_Py_C1002001002+ABY*I_ERI_Dyz_S_D2z_Py_C1002001002;
  abcd[317] = I_ERI_Fy2z_S_D2z_Py_C1002001002+ABY*I_ERI_D2z_S_D2z_Py_C1002001002;
  abcd[318] = I_ERI_F2xz_S_D2z_Py_C1002001002+ABZ*I_ERI_D2x_S_D2z_Py_C1002001002;
  abcd[319] = I_ERI_Fxyz_S_D2z_Py_C1002001002+ABZ*I_ERI_Dxy_S_D2z_Py_C1002001002;
  abcd[320] = I_ERI_Fx2z_S_D2z_Py_C1002001002+ABZ*I_ERI_Dxz_S_D2z_Py_C1002001002;
  abcd[321] = I_ERI_F2yz_S_D2z_Py_C1002001002+ABZ*I_ERI_D2y_S_D2z_Py_C1002001002;
  abcd[322] = I_ERI_Fy2z_S_D2z_Py_C1002001002+ABZ*I_ERI_Dyz_S_D2z_Py_C1002001002;
  abcd[323] = I_ERI_F3z_S_D2z_Py_C1002001002+ABZ*I_ERI_D2z_S_D2z_Py_C1002001002;
  abcd[324] = I_ERI_F3x_S_D2x_Pz_C1002001002+ABX*I_ERI_D2x_S_D2x_Pz_C1002001002;
  abcd[325] = I_ERI_F2xy_S_D2x_Pz_C1002001002+ABX*I_ERI_Dxy_S_D2x_Pz_C1002001002;
  abcd[326] = I_ERI_F2xz_S_D2x_Pz_C1002001002+ABX*I_ERI_Dxz_S_D2x_Pz_C1002001002;
  abcd[327] = I_ERI_Fx2y_S_D2x_Pz_C1002001002+ABX*I_ERI_D2y_S_D2x_Pz_C1002001002;
  abcd[328] = I_ERI_Fxyz_S_D2x_Pz_C1002001002+ABX*I_ERI_Dyz_S_D2x_Pz_C1002001002;
  abcd[329] = I_ERI_Fx2z_S_D2x_Pz_C1002001002+ABX*I_ERI_D2z_S_D2x_Pz_C1002001002;
  abcd[330] = I_ERI_F2xy_S_D2x_Pz_C1002001002+ABY*I_ERI_D2x_S_D2x_Pz_C1002001002;
  abcd[331] = I_ERI_Fx2y_S_D2x_Pz_C1002001002+ABY*I_ERI_Dxy_S_D2x_Pz_C1002001002;
  abcd[332] = I_ERI_Fxyz_S_D2x_Pz_C1002001002+ABY*I_ERI_Dxz_S_D2x_Pz_C1002001002;
  abcd[333] = I_ERI_F3y_S_D2x_Pz_C1002001002+ABY*I_ERI_D2y_S_D2x_Pz_C1002001002;
  abcd[334] = I_ERI_F2yz_S_D2x_Pz_C1002001002+ABY*I_ERI_Dyz_S_D2x_Pz_C1002001002;
  abcd[335] = I_ERI_Fy2z_S_D2x_Pz_C1002001002+ABY*I_ERI_D2z_S_D2x_Pz_C1002001002;
  abcd[336] = I_ERI_F2xz_S_D2x_Pz_C1002001002+ABZ*I_ERI_D2x_S_D2x_Pz_C1002001002;
  abcd[337] = I_ERI_Fxyz_S_D2x_Pz_C1002001002+ABZ*I_ERI_Dxy_S_D2x_Pz_C1002001002;
  abcd[338] = I_ERI_Fx2z_S_D2x_Pz_C1002001002+ABZ*I_ERI_Dxz_S_D2x_Pz_C1002001002;
  abcd[339] = I_ERI_F2yz_S_D2x_Pz_C1002001002+ABZ*I_ERI_D2y_S_D2x_Pz_C1002001002;
  abcd[340] = I_ERI_Fy2z_S_D2x_Pz_C1002001002+ABZ*I_ERI_Dyz_S_D2x_Pz_C1002001002;
  abcd[341] = I_ERI_F3z_S_D2x_Pz_C1002001002+ABZ*I_ERI_D2z_S_D2x_Pz_C1002001002;
  abcd[342] = I_ERI_F3x_S_Dxy_Pz_C1002001002+ABX*I_ERI_D2x_S_Dxy_Pz_C1002001002;
  abcd[343] = I_ERI_F2xy_S_Dxy_Pz_C1002001002+ABX*I_ERI_Dxy_S_Dxy_Pz_C1002001002;
  abcd[344] = I_ERI_F2xz_S_Dxy_Pz_C1002001002+ABX*I_ERI_Dxz_S_Dxy_Pz_C1002001002;
  abcd[345] = I_ERI_Fx2y_S_Dxy_Pz_C1002001002+ABX*I_ERI_D2y_S_Dxy_Pz_C1002001002;
  abcd[346] = I_ERI_Fxyz_S_Dxy_Pz_C1002001002+ABX*I_ERI_Dyz_S_Dxy_Pz_C1002001002;
  abcd[347] = I_ERI_Fx2z_S_Dxy_Pz_C1002001002+ABX*I_ERI_D2z_S_Dxy_Pz_C1002001002;
  abcd[348] = I_ERI_F2xy_S_Dxy_Pz_C1002001002+ABY*I_ERI_D2x_S_Dxy_Pz_C1002001002;
  abcd[349] = I_ERI_Fx2y_S_Dxy_Pz_C1002001002+ABY*I_ERI_Dxy_S_Dxy_Pz_C1002001002;
  abcd[350] = I_ERI_Fxyz_S_Dxy_Pz_C1002001002+ABY*I_ERI_Dxz_S_Dxy_Pz_C1002001002;
  abcd[351] = I_ERI_F3y_S_Dxy_Pz_C1002001002+ABY*I_ERI_D2y_S_Dxy_Pz_C1002001002;
  abcd[352] = I_ERI_F2yz_S_Dxy_Pz_C1002001002+ABY*I_ERI_Dyz_S_Dxy_Pz_C1002001002;
  abcd[353] = I_ERI_Fy2z_S_Dxy_Pz_C1002001002+ABY*I_ERI_D2z_S_Dxy_Pz_C1002001002;
  abcd[354] = I_ERI_F2xz_S_Dxy_Pz_C1002001002+ABZ*I_ERI_D2x_S_Dxy_Pz_C1002001002;
  abcd[355] = I_ERI_Fxyz_S_Dxy_Pz_C1002001002+ABZ*I_ERI_Dxy_S_Dxy_Pz_C1002001002;
  abcd[356] = I_ERI_Fx2z_S_Dxy_Pz_C1002001002+ABZ*I_ERI_Dxz_S_Dxy_Pz_C1002001002;
  abcd[357] = I_ERI_F2yz_S_Dxy_Pz_C1002001002+ABZ*I_ERI_D2y_S_Dxy_Pz_C1002001002;
  abcd[358] = I_ERI_Fy2z_S_Dxy_Pz_C1002001002+ABZ*I_ERI_Dyz_S_Dxy_Pz_C1002001002;
  abcd[359] = I_ERI_F3z_S_Dxy_Pz_C1002001002+ABZ*I_ERI_D2z_S_Dxy_Pz_C1002001002;
  abcd[360] = I_ERI_F3x_S_Dxz_Pz_C1002001002+ABX*I_ERI_D2x_S_Dxz_Pz_C1002001002;
  abcd[361] = I_ERI_F2xy_S_Dxz_Pz_C1002001002+ABX*I_ERI_Dxy_S_Dxz_Pz_C1002001002;
  abcd[362] = I_ERI_F2xz_S_Dxz_Pz_C1002001002+ABX*I_ERI_Dxz_S_Dxz_Pz_C1002001002;
  abcd[363] = I_ERI_Fx2y_S_Dxz_Pz_C1002001002+ABX*I_ERI_D2y_S_Dxz_Pz_C1002001002;
  abcd[364] = I_ERI_Fxyz_S_Dxz_Pz_C1002001002+ABX*I_ERI_Dyz_S_Dxz_Pz_C1002001002;
  abcd[365] = I_ERI_Fx2z_S_Dxz_Pz_C1002001002+ABX*I_ERI_D2z_S_Dxz_Pz_C1002001002;
  abcd[366] = I_ERI_F2xy_S_Dxz_Pz_C1002001002+ABY*I_ERI_D2x_S_Dxz_Pz_C1002001002;
  abcd[367] = I_ERI_Fx2y_S_Dxz_Pz_C1002001002+ABY*I_ERI_Dxy_S_Dxz_Pz_C1002001002;
  abcd[368] = I_ERI_Fxyz_S_Dxz_Pz_C1002001002+ABY*I_ERI_Dxz_S_Dxz_Pz_C1002001002;
  abcd[369] = I_ERI_F3y_S_Dxz_Pz_C1002001002+ABY*I_ERI_D2y_S_Dxz_Pz_C1002001002;
  abcd[370] = I_ERI_F2yz_S_Dxz_Pz_C1002001002+ABY*I_ERI_Dyz_S_Dxz_Pz_C1002001002;
  abcd[371] = I_ERI_Fy2z_S_Dxz_Pz_C1002001002+ABY*I_ERI_D2z_S_Dxz_Pz_C1002001002;
  abcd[372] = I_ERI_F2xz_S_Dxz_Pz_C1002001002+ABZ*I_ERI_D2x_S_Dxz_Pz_C1002001002;
  abcd[373] = I_ERI_Fxyz_S_Dxz_Pz_C1002001002+ABZ*I_ERI_Dxy_S_Dxz_Pz_C1002001002;
  abcd[374] = I_ERI_Fx2z_S_Dxz_Pz_C1002001002+ABZ*I_ERI_Dxz_S_Dxz_Pz_C1002001002;
  abcd[375] = I_ERI_F2yz_S_Dxz_Pz_C1002001002+ABZ*I_ERI_D2y_S_Dxz_Pz_C1002001002;
  abcd[376] = I_ERI_Fy2z_S_Dxz_Pz_C1002001002+ABZ*I_ERI_Dyz_S_Dxz_Pz_C1002001002;
  abcd[377] = I_ERI_F3z_S_Dxz_Pz_C1002001002+ABZ*I_ERI_D2z_S_Dxz_Pz_C1002001002;
  abcd[378] = I_ERI_F3x_S_D2y_Pz_C1002001002+ABX*I_ERI_D2x_S_D2y_Pz_C1002001002;
  abcd[379] = I_ERI_F2xy_S_D2y_Pz_C1002001002+ABX*I_ERI_Dxy_S_D2y_Pz_C1002001002;
  abcd[380] = I_ERI_F2xz_S_D2y_Pz_C1002001002+ABX*I_ERI_Dxz_S_D2y_Pz_C1002001002;
  abcd[381] = I_ERI_Fx2y_S_D2y_Pz_C1002001002+ABX*I_ERI_D2y_S_D2y_Pz_C1002001002;
  abcd[382] = I_ERI_Fxyz_S_D2y_Pz_C1002001002+ABX*I_ERI_Dyz_S_D2y_Pz_C1002001002;
  abcd[383] = I_ERI_Fx2z_S_D2y_Pz_C1002001002+ABX*I_ERI_D2z_S_D2y_Pz_C1002001002;
  abcd[384] = I_ERI_F2xy_S_D2y_Pz_C1002001002+ABY*I_ERI_D2x_S_D2y_Pz_C1002001002;
  abcd[385] = I_ERI_Fx2y_S_D2y_Pz_C1002001002+ABY*I_ERI_Dxy_S_D2y_Pz_C1002001002;
  abcd[386] = I_ERI_Fxyz_S_D2y_Pz_C1002001002+ABY*I_ERI_Dxz_S_D2y_Pz_C1002001002;
  abcd[387] = I_ERI_F3y_S_D2y_Pz_C1002001002+ABY*I_ERI_D2y_S_D2y_Pz_C1002001002;
  abcd[388] = I_ERI_F2yz_S_D2y_Pz_C1002001002+ABY*I_ERI_Dyz_S_D2y_Pz_C1002001002;
  abcd[389] = I_ERI_Fy2z_S_D2y_Pz_C1002001002+ABY*I_ERI_D2z_S_D2y_Pz_C1002001002;
  abcd[390] = I_ERI_F2xz_S_D2y_Pz_C1002001002+ABZ*I_ERI_D2x_S_D2y_Pz_C1002001002;
  abcd[391] = I_ERI_Fxyz_S_D2y_Pz_C1002001002+ABZ*I_ERI_Dxy_S_D2y_Pz_C1002001002;
  abcd[392] = I_ERI_Fx2z_S_D2y_Pz_C1002001002+ABZ*I_ERI_Dxz_S_D2y_Pz_C1002001002;
  abcd[393] = I_ERI_F2yz_S_D2y_Pz_C1002001002+ABZ*I_ERI_D2y_S_D2y_Pz_C1002001002;
  abcd[394] = I_ERI_Fy2z_S_D2y_Pz_C1002001002+ABZ*I_ERI_Dyz_S_D2y_Pz_C1002001002;
  abcd[395] = I_ERI_F3z_S_D2y_Pz_C1002001002+ABZ*I_ERI_D2z_S_D2y_Pz_C1002001002;
  abcd[396] = I_ERI_F3x_S_Dyz_Pz_C1002001002+ABX*I_ERI_D2x_S_Dyz_Pz_C1002001002;
  abcd[397] = I_ERI_F2xy_S_Dyz_Pz_C1002001002+ABX*I_ERI_Dxy_S_Dyz_Pz_C1002001002;
  abcd[398] = I_ERI_F2xz_S_Dyz_Pz_C1002001002+ABX*I_ERI_Dxz_S_Dyz_Pz_C1002001002;
  abcd[399] = I_ERI_Fx2y_S_Dyz_Pz_C1002001002+ABX*I_ERI_D2y_S_Dyz_Pz_C1002001002;
  abcd[400] = I_ERI_Fxyz_S_Dyz_Pz_C1002001002+ABX*I_ERI_Dyz_S_Dyz_Pz_C1002001002;
  abcd[401] = I_ERI_Fx2z_S_Dyz_Pz_C1002001002+ABX*I_ERI_D2z_S_Dyz_Pz_C1002001002;
  abcd[402] = I_ERI_F2xy_S_Dyz_Pz_C1002001002+ABY*I_ERI_D2x_S_Dyz_Pz_C1002001002;
  abcd[403] = I_ERI_Fx2y_S_Dyz_Pz_C1002001002+ABY*I_ERI_Dxy_S_Dyz_Pz_C1002001002;
  abcd[404] = I_ERI_Fxyz_S_Dyz_Pz_C1002001002+ABY*I_ERI_Dxz_S_Dyz_Pz_C1002001002;
  abcd[405] = I_ERI_F3y_S_Dyz_Pz_C1002001002+ABY*I_ERI_D2y_S_Dyz_Pz_C1002001002;
  abcd[406] = I_ERI_F2yz_S_Dyz_Pz_C1002001002+ABY*I_ERI_Dyz_S_Dyz_Pz_C1002001002;
  abcd[407] = I_ERI_Fy2z_S_Dyz_Pz_C1002001002+ABY*I_ERI_D2z_S_Dyz_Pz_C1002001002;
  abcd[408] = I_ERI_F2xz_S_Dyz_Pz_C1002001002+ABZ*I_ERI_D2x_S_Dyz_Pz_C1002001002;
  abcd[409] = I_ERI_Fxyz_S_Dyz_Pz_C1002001002+ABZ*I_ERI_Dxy_S_Dyz_Pz_C1002001002;
  abcd[410] = I_ERI_Fx2z_S_Dyz_Pz_C1002001002+ABZ*I_ERI_Dxz_S_Dyz_Pz_C1002001002;
  abcd[411] = I_ERI_F2yz_S_Dyz_Pz_C1002001002+ABZ*I_ERI_D2y_S_Dyz_Pz_C1002001002;
  abcd[412] = I_ERI_Fy2z_S_Dyz_Pz_C1002001002+ABZ*I_ERI_Dyz_S_Dyz_Pz_C1002001002;
  abcd[413] = I_ERI_F3z_S_Dyz_Pz_C1002001002+ABZ*I_ERI_D2z_S_Dyz_Pz_C1002001002;
  abcd[414] = I_ERI_F3x_S_D2z_Pz_C1002001002+ABX*I_ERI_D2x_S_D2z_Pz_C1002001002;
  abcd[415] = I_ERI_F2xy_S_D2z_Pz_C1002001002+ABX*I_ERI_Dxy_S_D2z_Pz_C1002001002;
  abcd[416] = I_ERI_F2xz_S_D2z_Pz_C1002001002+ABX*I_ERI_Dxz_S_D2z_Pz_C1002001002;
  abcd[417] = I_ERI_Fx2y_S_D2z_Pz_C1002001002+ABX*I_ERI_D2y_S_D2z_Pz_C1002001002;
  abcd[418] = I_ERI_Fxyz_S_D2z_Pz_C1002001002+ABX*I_ERI_Dyz_S_D2z_Pz_C1002001002;
  abcd[419] = I_ERI_Fx2z_S_D2z_Pz_C1002001002+ABX*I_ERI_D2z_S_D2z_Pz_C1002001002;
  abcd[420] = I_ERI_F2xy_S_D2z_Pz_C1002001002+ABY*I_ERI_D2x_S_D2z_Pz_C1002001002;
  abcd[421] = I_ERI_Fx2y_S_D2z_Pz_C1002001002+ABY*I_ERI_Dxy_S_D2z_Pz_C1002001002;
  abcd[422] = I_ERI_Fxyz_S_D2z_Pz_C1002001002+ABY*I_ERI_Dxz_S_D2z_Pz_C1002001002;
  abcd[423] = I_ERI_F3y_S_D2z_Pz_C1002001002+ABY*I_ERI_D2y_S_D2z_Pz_C1002001002;
  abcd[424] = I_ERI_F2yz_S_D2z_Pz_C1002001002+ABY*I_ERI_Dyz_S_D2z_Pz_C1002001002;
  abcd[425] = I_ERI_Fy2z_S_D2z_Pz_C1002001002+ABY*I_ERI_D2z_S_D2z_Pz_C1002001002;
  abcd[426] = I_ERI_F2xz_S_D2z_Pz_C1002001002+ABZ*I_ERI_D2x_S_D2z_Pz_C1002001002;
  abcd[427] = I_ERI_Fxyz_S_D2z_Pz_C1002001002+ABZ*I_ERI_Dxy_S_D2z_Pz_C1002001002;
  abcd[428] = I_ERI_Fx2z_S_D2z_Pz_C1002001002+ABZ*I_ERI_Dxz_S_D2z_Pz_C1002001002;
  abcd[429] = I_ERI_F2yz_S_D2z_Pz_C1002001002+ABZ*I_ERI_D2y_S_D2z_Pz_C1002001002;
  abcd[430] = I_ERI_Fy2z_S_D2z_Pz_C1002001002+ABZ*I_ERI_Dyz_S_D2z_Pz_C1002001002;
  abcd[431] = I_ERI_F3z_S_D2z_Pz_C1002001002+ABZ*I_ERI_D2z_S_D2z_Pz_C1002001002;
}
