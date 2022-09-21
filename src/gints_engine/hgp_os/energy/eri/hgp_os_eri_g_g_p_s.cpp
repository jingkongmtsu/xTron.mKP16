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

void hgp_os_eri_g_g_p_s(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_L8x_S_Px_S = 0.0E0;
  Double I_ERI_L7xy_S_Px_S = 0.0E0;
  Double I_ERI_L7xz_S_Px_S = 0.0E0;
  Double I_ERI_L6x2y_S_Px_S = 0.0E0;
  Double I_ERI_L6xyz_S_Px_S = 0.0E0;
  Double I_ERI_L6x2z_S_Px_S = 0.0E0;
  Double I_ERI_L5x3y_S_Px_S = 0.0E0;
  Double I_ERI_L5x2yz_S_Px_S = 0.0E0;
  Double I_ERI_L5xy2z_S_Px_S = 0.0E0;
  Double I_ERI_L5x3z_S_Px_S = 0.0E0;
  Double I_ERI_L4x4y_S_Px_S = 0.0E0;
  Double I_ERI_L4x3yz_S_Px_S = 0.0E0;
  Double I_ERI_L4x2y2z_S_Px_S = 0.0E0;
  Double I_ERI_L4xy3z_S_Px_S = 0.0E0;
  Double I_ERI_L4x4z_S_Px_S = 0.0E0;
  Double I_ERI_L3x5y_S_Px_S = 0.0E0;
  Double I_ERI_L3x4yz_S_Px_S = 0.0E0;
  Double I_ERI_L3x3y2z_S_Px_S = 0.0E0;
  Double I_ERI_L3x2y3z_S_Px_S = 0.0E0;
  Double I_ERI_L3xy4z_S_Px_S = 0.0E0;
  Double I_ERI_L3x5z_S_Px_S = 0.0E0;
  Double I_ERI_L2x6y_S_Px_S = 0.0E0;
  Double I_ERI_L2x5yz_S_Px_S = 0.0E0;
  Double I_ERI_L2x4y2z_S_Px_S = 0.0E0;
  Double I_ERI_L2x3y3z_S_Px_S = 0.0E0;
  Double I_ERI_L2x2y4z_S_Px_S = 0.0E0;
  Double I_ERI_L2xy5z_S_Px_S = 0.0E0;
  Double I_ERI_L2x6z_S_Px_S = 0.0E0;
  Double I_ERI_Lx7y_S_Px_S = 0.0E0;
  Double I_ERI_Lx6yz_S_Px_S = 0.0E0;
  Double I_ERI_Lx5y2z_S_Px_S = 0.0E0;
  Double I_ERI_Lx4y3z_S_Px_S = 0.0E0;
  Double I_ERI_Lx3y4z_S_Px_S = 0.0E0;
  Double I_ERI_Lx2y5z_S_Px_S = 0.0E0;
  Double I_ERI_Lxy6z_S_Px_S = 0.0E0;
  Double I_ERI_Lx7z_S_Px_S = 0.0E0;
  Double I_ERI_L8y_S_Px_S = 0.0E0;
  Double I_ERI_L7yz_S_Px_S = 0.0E0;
  Double I_ERI_L6y2z_S_Px_S = 0.0E0;
  Double I_ERI_L5y3z_S_Px_S = 0.0E0;
  Double I_ERI_L4y4z_S_Px_S = 0.0E0;
  Double I_ERI_L3y5z_S_Px_S = 0.0E0;
  Double I_ERI_L2y6z_S_Px_S = 0.0E0;
  Double I_ERI_Ly7z_S_Px_S = 0.0E0;
  Double I_ERI_L8z_S_Px_S = 0.0E0;
  Double I_ERI_L8x_S_Py_S = 0.0E0;
  Double I_ERI_L7xy_S_Py_S = 0.0E0;
  Double I_ERI_L7xz_S_Py_S = 0.0E0;
  Double I_ERI_L6x2y_S_Py_S = 0.0E0;
  Double I_ERI_L6xyz_S_Py_S = 0.0E0;
  Double I_ERI_L6x2z_S_Py_S = 0.0E0;
  Double I_ERI_L5x3y_S_Py_S = 0.0E0;
  Double I_ERI_L5x2yz_S_Py_S = 0.0E0;
  Double I_ERI_L5xy2z_S_Py_S = 0.0E0;
  Double I_ERI_L5x3z_S_Py_S = 0.0E0;
  Double I_ERI_L4x4y_S_Py_S = 0.0E0;
  Double I_ERI_L4x3yz_S_Py_S = 0.0E0;
  Double I_ERI_L4x2y2z_S_Py_S = 0.0E0;
  Double I_ERI_L4xy3z_S_Py_S = 0.0E0;
  Double I_ERI_L4x4z_S_Py_S = 0.0E0;
  Double I_ERI_L3x5y_S_Py_S = 0.0E0;
  Double I_ERI_L3x4yz_S_Py_S = 0.0E0;
  Double I_ERI_L3x3y2z_S_Py_S = 0.0E0;
  Double I_ERI_L3x2y3z_S_Py_S = 0.0E0;
  Double I_ERI_L3xy4z_S_Py_S = 0.0E0;
  Double I_ERI_L3x5z_S_Py_S = 0.0E0;
  Double I_ERI_L2x6y_S_Py_S = 0.0E0;
  Double I_ERI_L2x5yz_S_Py_S = 0.0E0;
  Double I_ERI_L2x4y2z_S_Py_S = 0.0E0;
  Double I_ERI_L2x3y3z_S_Py_S = 0.0E0;
  Double I_ERI_L2x2y4z_S_Py_S = 0.0E0;
  Double I_ERI_L2xy5z_S_Py_S = 0.0E0;
  Double I_ERI_L2x6z_S_Py_S = 0.0E0;
  Double I_ERI_Lx7y_S_Py_S = 0.0E0;
  Double I_ERI_Lx6yz_S_Py_S = 0.0E0;
  Double I_ERI_Lx5y2z_S_Py_S = 0.0E0;
  Double I_ERI_Lx4y3z_S_Py_S = 0.0E0;
  Double I_ERI_Lx3y4z_S_Py_S = 0.0E0;
  Double I_ERI_Lx2y5z_S_Py_S = 0.0E0;
  Double I_ERI_Lxy6z_S_Py_S = 0.0E0;
  Double I_ERI_Lx7z_S_Py_S = 0.0E0;
  Double I_ERI_L8y_S_Py_S = 0.0E0;
  Double I_ERI_L7yz_S_Py_S = 0.0E0;
  Double I_ERI_L6y2z_S_Py_S = 0.0E0;
  Double I_ERI_L5y3z_S_Py_S = 0.0E0;
  Double I_ERI_L4y4z_S_Py_S = 0.0E0;
  Double I_ERI_L3y5z_S_Py_S = 0.0E0;
  Double I_ERI_L2y6z_S_Py_S = 0.0E0;
  Double I_ERI_Ly7z_S_Py_S = 0.0E0;
  Double I_ERI_L8z_S_Py_S = 0.0E0;
  Double I_ERI_L8x_S_Pz_S = 0.0E0;
  Double I_ERI_L7xy_S_Pz_S = 0.0E0;
  Double I_ERI_L7xz_S_Pz_S = 0.0E0;
  Double I_ERI_L6x2y_S_Pz_S = 0.0E0;
  Double I_ERI_L6xyz_S_Pz_S = 0.0E0;
  Double I_ERI_L6x2z_S_Pz_S = 0.0E0;
  Double I_ERI_L5x3y_S_Pz_S = 0.0E0;
  Double I_ERI_L5x2yz_S_Pz_S = 0.0E0;
  Double I_ERI_L5xy2z_S_Pz_S = 0.0E0;
  Double I_ERI_L5x3z_S_Pz_S = 0.0E0;
  Double I_ERI_L4x4y_S_Pz_S = 0.0E0;
  Double I_ERI_L4x3yz_S_Pz_S = 0.0E0;
  Double I_ERI_L4x2y2z_S_Pz_S = 0.0E0;
  Double I_ERI_L4xy3z_S_Pz_S = 0.0E0;
  Double I_ERI_L4x4z_S_Pz_S = 0.0E0;
  Double I_ERI_L3x5y_S_Pz_S = 0.0E0;
  Double I_ERI_L3x4yz_S_Pz_S = 0.0E0;
  Double I_ERI_L3x3y2z_S_Pz_S = 0.0E0;
  Double I_ERI_L3x2y3z_S_Pz_S = 0.0E0;
  Double I_ERI_L3xy4z_S_Pz_S = 0.0E0;
  Double I_ERI_L3x5z_S_Pz_S = 0.0E0;
  Double I_ERI_L2x6y_S_Pz_S = 0.0E0;
  Double I_ERI_L2x5yz_S_Pz_S = 0.0E0;
  Double I_ERI_L2x4y2z_S_Pz_S = 0.0E0;
  Double I_ERI_L2x3y3z_S_Pz_S = 0.0E0;
  Double I_ERI_L2x2y4z_S_Pz_S = 0.0E0;
  Double I_ERI_L2xy5z_S_Pz_S = 0.0E0;
  Double I_ERI_L2x6z_S_Pz_S = 0.0E0;
  Double I_ERI_Lx7y_S_Pz_S = 0.0E0;
  Double I_ERI_Lx6yz_S_Pz_S = 0.0E0;
  Double I_ERI_Lx5y2z_S_Pz_S = 0.0E0;
  Double I_ERI_Lx4y3z_S_Pz_S = 0.0E0;
  Double I_ERI_Lx3y4z_S_Pz_S = 0.0E0;
  Double I_ERI_Lx2y5z_S_Pz_S = 0.0E0;
  Double I_ERI_Lxy6z_S_Pz_S = 0.0E0;
  Double I_ERI_Lx7z_S_Pz_S = 0.0E0;
  Double I_ERI_L8y_S_Pz_S = 0.0E0;
  Double I_ERI_L7yz_S_Pz_S = 0.0E0;
  Double I_ERI_L6y2z_S_Pz_S = 0.0E0;
  Double I_ERI_L5y3z_S_Pz_S = 0.0E0;
  Double I_ERI_L4y4z_S_Pz_S = 0.0E0;
  Double I_ERI_L3y5z_S_Pz_S = 0.0E0;
  Double I_ERI_L2y6z_S_Pz_S = 0.0E0;
  Double I_ERI_Ly7z_S_Pz_S = 0.0E0;
  Double I_ERI_L8z_S_Pz_S = 0.0E0;
  Double I_ERI_K7x_S_Px_S = 0.0E0;
  Double I_ERI_K6xy_S_Px_S = 0.0E0;
  Double I_ERI_K6xz_S_Px_S = 0.0E0;
  Double I_ERI_K5x2y_S_Px_S = 0.0E0;
  Double I_ERI_K5xyz_S_Px_S = 0.0E0;
  Double I_ERI_K5x2z_S_Px_S = 0.0E0;
  Double I_ERI_K4x3y_S_Px_S = 0.0E0;
  Double I_ERI_K4x2yz_S_Px_S = 0.0E0;
  Double I_ERI_K4xy2z_S_Px_S = 0.0E0;
  Double I_ERI_K4x3z_S_Px_S = 0.0E0;
  Double I_ERI_K3x4y_S_Px_S = 0.0E0;
  Double I_ERI_K3x3yz_S_Px_S = 0.0E0;
  Double I_ERI_K3x2y2z_S_Px_S = 0.0E0;
  Double I_ERI_K3xy3z_S_Px_S = 0.0E0;
  Double I_ERI_K3x4z_S_Px_S = 0.0E0;
  Double I_ERI_K2x5y_S_Px_S = 0.0E0;
  Double I_ERI_K2x4yz_S_Px_S = 0.0E0;
  Double I_ERI_K2x3y2z_S_Px_S = 0.0E0;
  Double I_ERI_K2x2y3z_S_Px_S = 0.0E0;
  Double I_ERI_K2xy4z_S_Px_S = 0.0E0;
  Double I_ERI_K2x5z_S_Px_S = 0.0E0;
  Double I_ERI_Kx6y_S_Px_S = 0.0E0;
  Double I_ERI_Kx5yz_S_Px_S = 0.0E0;
  Double I_ERI_Kx4y2z_S_Px_S = 0.0E0;
  Double I_ERI_Kx3y3z_S_Px_S = 0.0E0;
  Double I_ERI_Kx2y4z_S_Px_S = 0.0E0;
  Double I_ERI_Kxy5z_S_Px_S = 0.0E0;
  Double I_ERI_Kx6z_S_Px_S = 0.0E0;
  Double I_ERI_K7y_S_Px_S = 0.0E0;
  Double I_ERI_K6yz_S_Px_S = 0.0E0;
  Double I_ERI_K5y2z_S_Px_S = 0.0E0;
  Double I_ERI_K4y3z_S_Px_S = 0.0E0;
  Double I_ERI_K3y4z_S_Px_S = 0.0E0;
  Double I_ERI_K2y5z_S_Px_S = 0.0E0;
  Double I_ERI_Ky6z_S_Px_S = 0.0E0;
  Double I_ERI_K7z_S_Px_S = 0.0E0;
  Double I_ERI_K7x_S_Py_S = 0.0E0;
  Double I_ERI_K6xy_S_Py_S = 0.0E0;
  Double I_ERI_K6xz_S_Py_S = 0.0E0;
  Double I_ERI_K5x2y_S_Py_S = 0.0E0;
  Double I_ERI_K5xyz_S_Py_S = 0.0E0;
  Double I_ERI_K5x2z_S_Py_S = 0.0E0;
  Double I_ERI_K4x3y_S_Py_S = 0.0E0;
  Double I_ERI_K4x2yz_S_Py_S = 0.0E0;
  Double I_ERI_K4xy2z_S_Py_S = 0.0E0;
  Double I_ERI_K4x3z_S_Py_S = 0.0E0;
  Double I_ERI_K3x4y_S_Py_S = 0.0E0;
  Double I_ERI_K3x3yz_S_Py_S = 0.0E0;
  Double I_ERI_K3x2y2z_S_Py_S = 0.0E0;
  Double I_ERI_K3xy3z_S_Py_S = 0.0E0;
  Double I_ERI_K3x4z_S_Py_S = 0.0E0;
  Double I_ERI_K2x5y_S_Py_S = 0.0E0;
  Double I_ERI_K2x4yz_S_Py_S = 0.0E0;
  Double I_ERI_K2x3y2z_S_Py_S = 0.0E0;
  Double I_ERI_K2x2y3z_S_Py_S = 0.0E0;
  Double I_ERI_K2xy4z_S_Py_S = 0.0E0;
  Double I_ERI_K2x5z_S_Py_S = 0.0E0;
  Double I_ERI_Kx6y_S_Py_S = 0.0E0;
  Double I_ERI_Kx5yz_S_Py_S = 0.0E0;
  Double I_ERI_Kx4y2z_S_Py_S = 0.0E0;
  Double I_ERI_Kx3y3z_S_Py_S = 0.0E0;
  Double I_ERI_Kx2y4z_S_Py_S = 0.0E0;
  Double I_ERI_Kxy5z_S_Py_S = 0.0E0;
  Double I_ERI_Kx6z_S_Py_S = 0.0E0;
  Double I_ERI_K7y_S_Py_S = 0.0E0;
  Double I_ERI_K6yz_S_Py_S = 0.0E0;
  Double I_ERI_K5y2z_S_Py_S = 0.0E0;
  Double I_ERI_K4y3z_S_Py_S = 0.0E0;
  Double I_ERI_K3y4z_S_Py_S = 0.0E0;
  Double I_ERI_K2y5z_S_Py_S = 0.0E0;
  Double I_ERI_Ky6z_S_Py_S = 0.0E0;
  Double I_ERI_K7z_S_Py_S = 0.0E0;
  Double I_ERI_K7x_S_Pz_S = 0.0E0;
  Double I_ERI_K6xy_S_Pz_S = 0.0E0;
  Double I_ERI_K6xz_S_Pz_S = 0.0E0;
  Double I_ERI_K5x2y_S_Pz_S = 0.0E0;
  Double I_ERI_K5xyz_S_Pz_S = 0.0E0;
  Double I_ERI_K5x2z_S_Pz_S = 0.0E0;
  Double I_ERI_K4x3y_S_Pz_S = 0.0E0;
  Double I_ERI_K4x2yz_S_Pz_S = 0.0E0;
  Double I_ERI_K4xy2z_S_Pz_S = 0.0E0;
  Double I_ERI_K4x3z_S_Pz_S = 0.0E0;
  Double I_ERI_K3x4y_S_Pz_S = 0.0E0;
  Double I_ERI_K3x3yz_S_Pz_S = 0.0E0;
  Double I_ERI_K3x2y2z_S_Pz_S = 0.0E0;
  Double I_ERI_K3xy3z_S_Pz_S = 0.0E0;
  Double I_ERI_K3x4z_S_Pz_S = 0.0E0;
  Double I_ERI_K2x5y_S_Pz_S = 0.0E0;
  Double I_ERI_K2x4yz_S_Pz_S = 0.0E0;
  Double I_ERI_K2x3y2z_S_Pz_S = 0.0E0;
  Double I_ERI_K2x2y3z_S_Pz_S = 0.0E0;
  Double I_ERI_K2xy4z_S_Pz_S = 0.0E0;
  Double I_ERI_K2x5z_S_Pz_S = 0.0E0;
  Double I_ERI_Kx6y_S_Pz_S = 0.0E0;
  Double I_ERI_Kx5yz_S_Pz_S = 0.0E0;
  Double I_ERI_Kx4y2z_S_Pz_S = 0.0E0;
  Double I_ERI_Kx3y3z_S_Pz_S = 0.0E0;
  Double I_ERI_Kx2y4z_S_Pz_S = 0.0E0;
  Double I_ERI_Kxy5z_S_Pz_S = 0.0E0;
  Double I_ERI_Kx6z_S_Pz_S = 0.0E0;
  Double I_ERI_K7y_S_Pz_S = 0.0E0;
  Double I_ERI_K6yz_S_Pz_S = 0.0E0;
  Double I_ERI_K5y2z_S_Pz_S = 0.0E0;
  Double I_ERI_K4y3z_S_Pz_S = 0.0E0;
  Double I_ERI_K3y4z_S_Pz_S = 0.0E0;
  Double I_ERI_K2y5z_S_Pz_S = 0.0E0;
  Double I_ERI_Ky6z_S_Pz_S = 0.0E0;
  Double I_ERI_K7z_S_Pz_S = 0.0E0;
  Double I_ERI_I6x_S_Px_S = 0.0E0;
  Double I_ERI_I5xy_S_Px_S = 0.0E0;
  Double I_ERI_I5xz_S_Px_S = 0.0E0;
  Double I_ERI_I4x2y_S_Px_S = 0.0E0;
  Double I_ERI_I4xyz_S_Px_S = 0.0E0;
  Double I_ERI_I4x2z_S_Px_S = 0.0E0;
  Double I_ERI_I3x3y_S_Px_S = 0.0E0;
  Double I_ERI_I3x2yz_S_Px_S = 0.0E0;
  Double I_ERI_I3xy2z_S_Px_S = 0.0E0;
  Double I_ERI_I3x3z_S_Px_S = 0.0E0;
  Double I_ERI_I2x4y_S_Px_S = 0.0E0;
  Double I_ERI_I2x3yz_S_Px_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_Px_S = 0.0E0;
  Double I_ERI_I2xy3z_S_Px_S = 0.0E0;
  Double I_ERI_I2x4z_S_Px_S = 0.0E0;
  Double I_ERI_Ix5y_S_Px_S = 0.0E0;
  Double I_ERI_Ix4yz_S_Px_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_Px_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_Px_S = 0.0E0;
  Double I_ERI_Ixy4z_S_Px_S = 0.0E0;
  Double I_ERI_Ix5z_S_Px_S = 0.0E0;
  Double I_ERI_I6y_S_Px_S = 0.0E0;
  Double I_ERI_I5yz_S_Px_S = 0.0E0;
  Double I_ERI_I4y2z_S_Px_S = 0.0E0;
  Double I_ERI_I3y3z_S_Px_S = 0.0E0;
  Double I_ERI_I2y4z_S_Px_S = 0.0E0;
  Double I_ERI_Iy5z_S_Px_S = 0.0E0;
  Double I_ERI_I6z_S_Px_S = 0.0E0;
  Double I_ERI_I6x_S_Py_S = 0.0E0;
  Double I_ERI_I5xy_S_Py_S = 0.0E0;
  Double I_ERI_I5xz_S_Py_S = 0.0E0;
  Double I_ERI_I4x2y_S_Py_S = 0.0E0;
  Double I_ERI_I4xyz_S_Py_S = 0.0E0;
  Double I_ERI_I4x2z_S_Py_S = 0.0E0;
  Double I_ERI_I3x3y_S_Py_S = 0.0E0;
  Double I_ERI_I3x2yz_S_Py_S = 0.0E0;
  Double I_ERI_I3xy2z_S_Py_S = 0.0E0;
  Double I_ERI_I3x3z_S_Py_S = 0.0E0;
  Double I_ERI_I2x4y_S_Py_S = 0.0E0;
  Double I_ERI_I2x3yz_S_Py_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_Py_S = 0.0E0;
  Double I_ERI_I2xy3z_S_Py_S = 0.0E0;
  Double I_ERI_I2x4z_S_Py_S = 0.0E0;
  Double I_ERI_Ix5y_S_Py_S = 0.0E0;
  Double I_ERI_Ix4yz_S_Py_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_Py_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_Py_S = 0.0E0;
  Double I_ERI_Ixy4z_S_Py_S = 0.0E0;
  Double I_ERI_Ix5z_S_Py_S = 0.0E0;
  Double I_ERI_I6y_S_Py_S = 0.0E0;
  Double I_ERI_I5yz_S_Py_S = 0.0E0;
  Double I_ERI_I4y2z_S_Py_S = 0.0E0;
  Double I_ERI_I3y3z_S_Py_S = 0.0E0;
  Double I_ERI_I2y4z_S_Py_S = 0.0E0;
  Double I_ERI_Iy5z_S_Py_S = 0.0E0;
  Double I_ERI_I6z_S_Py_S = 0.0E0;
  Double I_ERI_I6x_S_Pz_S = 0.0E0;
  Double I_ERI_I5xy_S_Pz_S = 0.0E0;
  Double I_ERI_I5xz_S_Pz_S = 0.0E0;
  Double I_ERI_I4x2y_S_Pz_S = 0.0E0;
  Double I_ERI_I4xyz_S_Pz_S = 0.0E0;
  Double I_ERI_I4x2z_S_Pz_S = 0.0E0;
  Double I_ERI_I3x3y_S_Pz_S = 0.0E0;
  Double I_ERI_I3x2yz_S_Pz_S = 0.0E0;
  Double I_ERI_I3xy2z_S_Pz_S = 0.0E0;
  Double I_ERI_I3x3z_S_Pz_S = 0.0E0;
  Double I_ERI_I2x4y_S_Pz_S = 0.0E0;
  Double I_ERI_I2x3yz_S_Pz_S = 0.0E0;
  Double I_ERI_I2x2y2z_S_Pz_S = 0.0E0;
  Double I_ERI_I2xy3z_S_Pz_S = 0.0E0;
  Double I_ERI_I2x4z_S_Pz_S = 0.0E0;
  Double I_ERI_Ix5y_S_Pz_S = 0.0E0;
  Double I_ERI_Ix4yz_S_Pz_S = 0.0E0;
  Double I_ERI_Ix3y2z_S_Pz_S = 0.0E0;
  Double I_ERI_Ix2y3z_S_Pz_S = 0.0E0;
  Double I_ERI_Ixy4z_S_Pz_S = 0.0E0;
  Double I_ERI_Ix5z_S_Pz_S = 0.0E0;
  Double I_ERI_I6y_S_Pz_S = 0.0E0;
  Double I_ERI_I5yz_S_Pz_S = 0.0E0;
  Double I_ERI_I4y2z_S_Pz_S = 0.0E0;
  Double I_ERI_I3y3z_S_Pz_S = 0.0E0;
  Double I_ERI_I2y4z_S_Pz_S = 0.0E0;
  Double I_ERI_Iy5z_S_Pz_S = 0.0E0;
  Double I_ERI_I6z_S_Pz_S = 0.0E0;
  Double I_ERI_H5x_S_Px_S = 0.0E0;
  Double I_ERI_H4xy_S_Px_S = 0.0E0;
  Double I_ERI_H4xz_S_Px_S = 0.0E0;
  Double I_ERI_H3x2y_S_Px_S = 0.0E0;
  Double I_ERI_H3xyz_S_Px_S = 0.0E0;
  Double I_ERI_H3x2z_S_Px_S = 0.0E0;
  Double I_ERI_H2x3y_S_Px_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Px_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Px_S = 0.0E0;
  Double I_ERI_H2x3z_S_Px_S = 0.0E0;
  Double I_ERI_Hx4y_S_Px_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Px_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Px_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Px_S = 0.0E0;
  Double I_ERI_Hx4z_S_Px_S = 0.0E0;
  Double I_ERI_H5y_S_Px_S = 0.0E0;
  Double I_ERI_H4yz_S_Px_S = 0.0E0;
  Double I_ERI_H3y2z_S_Px_S = 0.0E0;
  Double I_ERI_H2y3z_S_Px_S = 0.0E0;
  Double I_ERI_Hy4z_S_Px_S = 0.0E0;
  Double I_ERI_H5z_S_Px_S = 0.0E0;
  Double I_ERI_H5x_S_Py_S = 0.0E0;
  Double I_ERI_H4xy_S_Py_S = 0.0E0;
  Double I_ERI_H4xz_S_Py_S = 0.0E0;
  Double I_ERI_H3x2y_S_Py_S = 0.0E0;
  Double I_ERI_H3xyz_S_Py_S = 0.0E0;
  Double I_ERI_H3x2z_S_Py_S = 0.0E0;
  Double I_ERI_H2x3y_S_Py_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Py_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Py_S = 0.0E0;
  Double I_ERI_H2x3z_S_Py_S = 0.0E0;
  Double I_ERI_Hx4y_S_Py_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Py_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Py_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Py_S = 0.0E0;
  Double I_ERI_Hx4z_S_Py_S = 0.0E0;
  Double I_ERI_H5y_S_Py_S = 0.0E0;
  Double I_ERI_H4yz_S_Py_S = 0.0E0;
  Double I_ERI_H3y2z_S_Py_S = 0.0E0;
  Double I_ERI_H2y3z_S_Py_S = 0.0E0;
  Double I_ERI_Hy4z_S_Py_S = 0.0E0;
  Double I_ERI_H5z_S_Py_S = 0.0E0;
  Double I_ERI_H5x_S_Pz_S = 0.0E0;
  Double I_ERI_H4xy_S_Pz_S = 0.0E0;
  Double I_ERI_H4xz_S_Pz_S = 0.0E0;
  Double I_ERI_H3x2y_S_Pz_S = 0.0E0;
  Double I_ERI_H3xyz_S_Pz_S = 0.0E0;
  Double I_ERI_H3x2z_S_Pz_S = 0.0E0;
  Double I_ERI_H2x3y_S_Pz_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Pz_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Pz_S = 0.0E0;
  Double I_ERI_H2x3z_S_Pz_S = 0.0E0;
  Double I_ERI_Hx4y_S_Pz_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Pz_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Pz_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Pz_S = 0.0E0;
  Double I_ERI_Hx4z_S_Pz_S = 0.0E0;
  Double I_ERI_H5y_S_Pz_S = 0.0E0;
  Double I_ERI_H4yz_S_Pz_S = 0.0E0;
  Double I_ERI_H3y2z_S_Pz_S = 0.0E0;
  Double I_ERI_H2y3z_S_Pz_S = 0.0E0;
  Double I_ERI_Hy4z_S_Pz_S = 0.0E0;
  Double I_ERI_H5z_S_Pz_S = 0.0E0;
  Double I_ERI_G4x_S_Px_S = 0.0E0;
  Double I_ERI_G3xy_S_Px_S = 0.0E0;
  Double I_ERI_G3xz_S_Px_S = 0.0E0;
  Double I_ERI_G2x2y_S_Px_S = 0.0E0;
  Double I_ERI_G2xyz_S_Px_S = 0.0E0;
  Double I_ERI_G2x2z_S_Px_S = 0.0E0;
  Double I_ERI_Gx3y_S_Px_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Px_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Px_S = 0.0E0;
  Double I_ERI_Gx3z_S_Px_S = 0.0E0;
  Double I_ERI_G4y_S_Px_S = 0.0E0;
  Double I_ERI_G3yz_S_Px_S = 0.0E0;
  Double I_ERI_G2y2z_S_Px_S = 0.0E0;
  Double I_ERI_Gy3z_S_Px_S = 0.0E0;
  Double I_ERI_G4z_S_Px_S = 0.0E0;
  Double I_ERI_G4x_S_Py_S = 0.0E0;
  Double I_ERI_G3xy_S_Py_S = 0.0E0;
  Double I_ERI_G3xz_S_Py_S = 0.0E0;
  Double I_ERI_G2x2y_S_Py_S = 0.0E0;
  Double I_ERI_G2xyz_S_Py_S = 0.0E0;
  Double I_ERI_G2x2z_S_Py_S = 0.0E0;
  Double I_ERI_Gx3y_S_Py_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Py_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Py_S = 0.0E0;
  Double I_ERI_Gx3z_S_Py_S = 0.0E0;
  Double I_ERI_G4y_S_Py_S = 0.0E0;
  Double I_ERI_G3yz_S_Py_S = 0.0E0;
  Double I_ERI_G2y2z_S_Py_S = 0.0E0;
  Double I_ERI_Gy3z_S_Py_S = 0.0E0;
  Double I_ERI_G4z_S_Py_S = 0.0E0;
  Double I_ERI_G4x_S_Pz_S = 0.0E0;
  Double I_ERI_G3xy_S_Pz_S = 0.0E0;
  Double I_ERI_G3xz_S_Pz_S = 0.0E0;
  Double I_ERI_G2x2y_S_Pz_S = 0.0E0;
  Double I_ERI_G2xyz_S_Pz_S = 0.0E0;
  Double I_ERI_G2x2z_S_Pz_S = 0.0E0;
  Double I_ERI_Gx3y_S_Pz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Pz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Pz_S = 0.0E0;
  Double I_ERI_Gx3z_S_Pz_S = 0.0E0;
  Double I_ERI_G4y_S_Pz_S = 0.0E0;
  Double I_ERI_G3yz_S_Pz_S = 0.0E0;
  Double I_ERI_G2y2z_S_Pz_S = 0.0E0;
  Double I_ERI_Gy3z_S_Pz_S = 0.0E0;
  Double I_ERI_G4z_S_Pz_S = 0.0E0;

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
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER53;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER51*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER49*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER47*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER45*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER43*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER41*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER39*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER37*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER35*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER33*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER31*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER29*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER27*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER25*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER23*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = 1.0E0+u2*ONEOVER21*I_ERI_S_S_S_S_M9_vrr;
        I_ERI_S_S_S_S_M9_vrr = ONEOVER19*I_ERI_S_S_S_S_M9_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ERI_S_S_S_S_M9_vrr  = f*I_ERI_S_S_S_S_M9_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
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
        I_ERI_S_S_S_S_M1_vrr = I_ERI_S_S_S_S_M1_vrr*erfPref_3;
        I_ERI_S_S_S_S_M2_vrr = I_ERI_S_S_S_S_M2_vrr*erfPref_5;
        I_ERI_S_S_S_S_M3_vrr = I_ERI_S_S_S_S_M3_vrr*erfPref_7;
        I_ERI_S_S_S_S_M4_vrr = I_ERI_S_S_S_S_M4_vrr*erfPref_9;
        I_ERI_S_S_S_S_M5_vrr = I_ERI_S_S_S_S_M5_vrr*erfPref_11;
        I_ERI_S_S_S_S_M6_vrr = I_ERI_S_S_S_S_M6_vrr*erfPref_13;
        I_ERI_S_S_S_S_M7_vrr = I_ERI_S_S_S_S_M7_vrr*erfPref_15;
        I_ERI_S_S_S_S_M8_vrr = I_ERI_S_S_S_S_M8_vrr*erfPref_17;
        I_ERI_S_S_S_S_M9_vrr = I_ERI_S_S_S_S_M9_vrr*erfPref_19;
      }

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M8
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M8
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M9
       ************************************************************/
      Double I_ERI_Px_S_S_S_M8_vrr = PAX*I_ERI_S_S_S_S_M8_vrr+WPX*I_ERI_S_S_S_S_M9_vrr;
      Double I_ERI_Py_S_S_S_M8_vrr = PAY*I_ERI_S_S_S_S_M8_vrr+WPY*I_ERI_S_S_S_S_M9_vrr;
      Double I_ERI_Pz_S_S_S_M8_vrr = PAZ*I_ERI_S_S_S_S_M8_vrr+WPZ*I_ERI_S_S_S_S_M9_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_P_S_S_S_M7
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M8
       ************************************************************/
      Double I_ERI_Px_S_S_S_M7_vrr = PAX*I_ERI_S_S_S_S_M7_vrr+WPX*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_Py_S_S_S_M7_vrr = PAY*I_ERI_S_S_S_S_M7_vrr+WPY*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_Pz_S_S_S_M7_vrr = PAZ*I_ERI_S_S_S_S_M7_vrr+WPZ*I_ERI_S_S_S_S_M8_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_D_S_S_S_M7
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M8
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M8
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M7_vrr = PAX*I_ERI_Px_S_S_S_M7_vrr+WPX*I_ERI_Px_S_S_S_M8_vrr+oned2z*I_ERI_S_S_S_S_M7_vrr-rhod2zsq*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_D2y_S_S_S_M7_vrr = PAY*I_ERI_Py_S_S_S_M7_vrr+WPY*I_ERI_Py_S_S_S_M8_vrr+oned2z*I_ERI_S_S_S_S_M7_vrr-rhod2zsq*I_ERI_S_S_S_S_M8_vrr;
      Double I_ERI_D2z_S_S_S_M7_vrr = PAZ*I_ERI_Pz_S_S_S_M7_vrr+WPZ*I_ERI_Pz_S_S_S_M8_vrr+oned2z*I_ERI_S_S_S_S_M7_vrr-rhod2zsq*I_ERI_S_S_S_S_M8_vrr;

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
       * shell quartet name: SQ_ERI_D_S_S_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M7
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M6_vrr = PAX*I_ERI_Px_S_S_S_M6_vrr+WPX*I_ERI_Px_S_S_S_M7_vrr+oned2z*I_ERI_S_S_S_S_M6_vrr-rhod2zsq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_D2y_S_S_S_M6_vrr = PAY*I_ERI_Py_S_S_S_M6_vrr+WPY*I_ERI_Py_S_S_S_M7_vrr+oned2z*I_ERI_S_S_S_S_M6_vrr-rhod2zsq*I_ERI_S_S_S_S_M7_vrr;
      Double I_ERI_D2z_S_S_S_M6_vrr = PAZ*I_ERI_Pz_S_S_S_M6_vrr+WPZ*I_ERI_Pz_S_S_S_M7_vrr+oned2z*I_ERI_S_S_S_S_M6_vrr-rhod2zsq*I_ERI_S_S_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M7
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M7
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M6_vrr = PAX*I_ERI_D2x_S_S_S_M6_vrr+WPX*I_ERI_D2x_S_S_S_M7_vrr+2*oned2z*I_ERI_Px_S_S_S_M6_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M7_vrr;
      Double I_ERI_F3y_S_S_S_M6_vrr = PAY*I_ERI_D2y_S_S_S_M6_vrr+WPY*I_ERI_D2y_S_S_S_M7_vrr+2*oned2z*I_ERI_Py_S_S_S_M6_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M7_vrr;
      Double I_ERI_F3z_S_S_S_M6_vrr = PAZ*I_ERI_D2z_S_S_S_M6_vrr+WPZ*I_ERI_D2z_S_S_S_M7_vrr+2*oned2z*I_ERI_Pz_S_S_S_M6_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M7_vrr;

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
       * shell quartet name: SQ_ERI_F_S_S_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M6
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M5_vrr = PAX*I_ERI_D2x_S_S_S_M5_vrr+WPX*I_ERI_D2x_S_S_S_M6_vrr+2*oned2z*I_ERI_Px_S_S_S_M5_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M6_vrr;
      Double I_ERI_F3y_S_S_S_M5_vrr = PAY*I_ERI_D2y_S_S_S_M5_vrr+WPY*I_ERI_D2y_S_S_S_M6_vrr+2*oned2z*I_ERI_Py_S_S_S_M5_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M6_vrr;
      Double I_ERI_F3z_S_S_S_M5_vrr = PAZ*I_ERI_D2z_S_S_S_M5_vrr+WPZ*I_ERI_D2z_S_S_S_M6_vrr+2*oned2z*I_ERI_Pz_S_S_S_M5_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M6
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M6
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M5_vrr = PAX*I_ERI_F3x_S_S_S_M5_vrr+WPX*I_ERI_F3x_S_S_S_M6_vrr+3*oned2z*I_ERI_D2x_S_S_S_M5_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M6_vrr;
      Double I_ERI_G3xy_S_S_S_M5_vrr = PAY*I_ERI_F3x_S_S_S_M5_vrr+WPY*I_ERI_F3x_S_S_S_M6_vrr;
      Double I_ERI_G3xz_S_S_S_M5_vrr = PAZ*I_ERI_F3x_S_S_S_M5_vrr+WPZ*I_ERI_F3x_S_S_S_M6_vrr;
      Double I_ERI_G4y_S_S_S_M5_vrr = PAY*I_ERI_F3y_S_S_S_M5_vrr+WPY*I_ERI_F3y_S_S_S_M6_vrr+3*oned2z*I_ERI_D2y_S_S_S_M5_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M6_vrr;
      Double I_ERI_G3yz_S_S_S_M5_vrr = PAZ*I_ERI_F3y_S_S_S_M5_vrr+WPZ*I_ERI_F3y_S_S_S_M6_vrr;
      Double I_ERI_G4z_S_S_S_M5_vrr = PAZ*I_ERI_F3z_S_S_S_M5_vrr+WPZ*I_ERI_F3z_S_S_S_M6_vrr+3*oned2z*I_ERI_D2z_S_S_S_M5_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M6_vrr;

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
       * shell quartet name: SQ_ERI_F_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M4_vrr = PAX*I_ERI_D2x_S_S_S_M4_vrr+WPX*I_ERI_D2x_S_S_S_M5_vrr+2*oned2z*I_ERI_Px_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_F2xy_S_S_S_M4_vrr = PAY*I_ERI_D2x_S_S_S_M4_vrr+WPY*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_F3y_S_S_S_M4_vrr = PAY*I_ERI_D2y_S_S_S_M4_vrr+WPY*I_ERI_D2y_S_S_S_M5_vrr+2*oned2z*I_ERI_Py_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_F3z_S_S_S_M4_vrr = PAZ*I_ERI_D2z_S_S_S_M4_vrr+WPZ*I_ERI_D2z_S_S_S_M5_vrr+2*oned2z*I_ERI_Pz_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M4_vrr = PAX*I_ERI_F3x_S_S_S_M4_vrr+WPX*I_ERI_F3x_S_S_S_M5_vrr+3*oned2z*I_ERI_D2x_S_S_S_M4_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_G3xy_S_S_S_M4_vrr = PAY*I_ERI_F3x_S_S_S_M4_vrr+WPY*I_ERI_F3x_S_S_S_M5_vrr;
      Double I_ERI_G3xz_S_S_S_M4_vrr = PAZ*I_ERI_F3x_S_S_S_M4_vrr+WPZ*I_ERI_F3x_S_S_S_M5_vrr;
      Double I_ERI_Gx3y_S_S_S_M4_vrr = PAX*I_ERI_F3y_S_S_S_M4_vrr+WPX*I_ERI_F3y_S_S_S_M5_vrr;
      Double I_ERI_Gx3z_S_S_S_M4_vrr = PAX*I_ERI_F3z_S_S_S_M4_vrr+WPX*I_ERI_F3z_S_S_S_M5_vrr;
      Double I_ERI_G4y_S_S_S_M4_vrr = PAY*I_ERI_F3y_S_S_S_M4_vrr+WPY*I_ERI_F3y_S_S_S_M5_vrr+3*oned2z*I_ERI_D2y_S_S_S_M4_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M5_vrr;
      Double I_ERI_G3yz_S_S_S_M4_vrr = PAZ*I_ERI_F3y_S_S_S_M4_vrr+WPZ*I_ERI_F3y_S_S_S_M5_vrr;
      Double I_ERI_G4z_S_S_S_M4_vrr = PAZ*I_ERI_F3z_S_S_S_M4_vrr+WPZ*I_ERI_F3z_S_S_S_M5_vrr+3*oned2z*I_ERI_D2z_S_S_S_M4_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M5
       ************************************************************/
      Double I_ERI_H5x_S_S_S_M4_vrr = PAX*I_ERI_G4x_S_S_S_M4_vrr+WPX*I_ERI_G4x_S_S_S_M5_vrr+4*oned2z*I_ERI_F3x_S_S_S_M4_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M5_vrr;
      Double I_ERI_H4xy_S_S_S_M4_vrr = PAY*I_ERI_G4x_S_S_S_M4_vrr+WPY*I_ERI_G4x_S_S_S_M5_vrr;
      Double I_ERI_H4xz_S_S_S_M4_vrr = PAZ*I_ERI_G4x_S_S_S_M4_vrr+WPZ*I_ERI_G4x_S_S_S_M5_vrr;
      Double I_ERI_H3x2y_S_S_S_M4_vrr = PAY*I_ERI_G3xy_S_S_S_M4_vrr+WPY*I_ERI_G3xy_S_S_S_M5_vrr+oned2z*I_ERI_F3x_S_S_S_M4_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M5_vrr;
      Double I_ERI_H3x2z_S_S_S_M4_vrr = PAZ*I_ERI_G3xz_S_S_S_M4_vrr+WPZ*I_ERI_G3xz_S_S_S_M5_vrr+oned2z*I_ERI_F3x_S_S_S_M4_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M5_vrr;
      Double I_ERI_Hx4y_S_S_S_M4_vrr = PAX*I_ERI_G4y_S_S_S_M4_vrr+WPX*I_ERI_G4y_S_S_S_M5_vrr;
      Double I_ERI_Hx4z_S_S_S_M4_vrr = PAX*I_ERI_G4z_S_S_S_M4_vrr+WPX*I_ERI_G4z_S_S_S_M5_vrr;
      Double I_ERI_H5y_S_S_S_M4_vrr = PAY*I_ERI_G4y_S_S_S_M4_vrr+WPY*I_ERI_G4y_S_S_S_M5_vrr+4*oned2z*I_ERI_F3y_S_S_S_M4_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M5_vrr;
      Double I_ERI_H4yz_S_S_S_M4_vrr = PAZ*I_ERI_G4y_S_S_S_M4_vrr+WPZ*I_ERI_G4y_S_S_S_M5_vrr;
      Double I_ERI_H3y2z_S_S_S_M4_vrr = PAZ*I_ERI_G3yz_S_S_S_M4_vrr+WPZ*I_ERI_G3yz_S_S_S_M5_vrr+oned2z*I_ERI_F3y_S_S_S_M4_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M5_vrr;
      Double I_ERI_Hy4z_S_S_S_M4_vrr = PAY*I_ERI_G4z_S_S_S_M4_vrr+WPY*I_ERI_G4z_S_S_S_M5_vrr;
      Double I_ERI_H5z_S_S_S_M4_vrr = PAZ*I_ERI_G4z_S_S_S_M4_vrr+WPZ*I_ERI_G4z_S_S_S_M5_vrr+4*oned2z*I_ERI_F3z_S_S_S_M4_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M5_vrr;

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
       * shell quartet name: SQ_ERI_G_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 5 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M3_vrr = PAX*I_ERI_F3x_S_S_S_M3_vrr+WPX*I_ERI_F3x_S_S_S_M4_vrr+3*oned2z*I_ERI_D2x_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G3xy_S_S_S_M3_vrr = PAY*I_ERI_F3x_S_S_S_M3_vrr+WPY*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G3xz_S_S_S_M3_vrr = PAZ*I_ERI_F3x_S_S_S_M3_vrr+WPZ*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G2x2y_S_S_S_M3_vrr = PAY*I_ERI_F2xy_S_S_S_M3_vrr+WPY*I_ERI_F2xy_S_S_S_M4_vrr+oned2z*I_ERI_D2x_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Gx3y_S_S_S_M3_vrr = PAX*I_ERI_F3y_S_S_S_M3_vrr+WPX*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_Gx3z_S_S_S_M3_vrr = PAX*I_ERI_F3z_S_S_S_M3_vrr+WPX*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4y_S_S_S_M3_vrr = PAY*I_ERI_F3y_S_S_S_M3_vrr+WPY*I_ERI_F3y_S_S_S_M4_vrr+3*oned2z*I_ERI_D2y_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_G3yz_S_S_S_M3_vrr = PAZ*I_ERI_F3y_S_S_S_M3_vrr+WPZ*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_Gy3z_S_S_S_M3_vrr = PAY*I_ERI_F3z_S_S_S_M3_vrr+WPY*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4z_S_S_S_M3_vrr = PAZ*I_ERI_F3z_S_S_S_M3_vrr+WPZ*I_ERI_F3z_S_S_S_M4_vrr+3*oned2z*I_ERI_D2z_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       ************************************************************/
      Double I_ERI_H5x_S_S_S_M3_vrr = PAX*I_ERI_G4x_S_S_S_M3_vrr+WPX*I_ERI_G4x_S_S_S_M4_vrr+4*oned2z*I_ERI_F3x_S_S_S_M3_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_H4xy_S_S_S_M3_vrr = PAY*I_ERI_G4x_S_S_S_M3_vrr+WPY*I_ERI_G4x_S_S_S_M4_vrr;
      Double I_ERI_H4xz_S_S_S_M3_vrr = PAZ*I_ERI_G4x_S_S_S_M3_vrr+WPZ*I_ERI_G4x_S_S_S_M4_vrr;
      Double I_ERI_H3x2y_S_S_S_M3_vrr = PAY*I_ERI_G3xy_S_S_S_M3_vrr+WPY*I_ERI_G3xy_S_S_S_M4_vrr+oned2z*I_ERI_F3x_S_S_S_M3_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_H3x2z_S_S_S_M3_vrr = PAZ*I_ERI_G3xz_S_S_S_M3_vrr+WPZ*I_ERI_G3xz_S_S_S_M4_vrr+oned2z*I_ERI_F3x_S_S_S_M3_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_H2x3y_S_S_S_M3_vrr = PAX*I_ERI_Gx3y_S_S_S_M3_vrr+WPX*I_ERI_Gx3y_S_S_S_M4_vrr+oned2z*I_ERI_F3y_S_S_S_M3_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_H2x3z_S_S_S_M3_vrr = PAX*I_ERI_Gx3z_S_S_S_M3_vrr+WPX*I_ERI_Gx3z_S_S_S_M4_vrr+oned2z*I_ERI_F3z_S_S_S_M3_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_Hx4y_S_S_S_M3_vrr = PAX*I_ERI_G4y_S_S_S_M3_vrr+WPX*I_ERI_G4y_S_S_S_M4_vrr;
      Double I_ERI_Hx4z_S_S_S_M3_vrr = PAX*I_ERI_G4z_S_S_S_M3_vrr+WPX*I_ERI_G4z_S_S_S_M4_vrr;
      Double I_ERI_H5y_S_S_S_M3_vrr = PAY*I_ERI_G4y_S_S_S_M3_vrr+WPY*I_ERI_G4y_S_S_S_M4_vrr+4*oned2z*I_ERI_F3y_S_S_S_M3_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_H4yz_S_S_S_M3_vrr = PAZ*I_ERI_G4y_S_S_S_M3_vrr+WPZ*I_ERI_G4y_S_S_S_M4_vrr;
      Double I_ERI_H3y2z_S_S_S_M3_vrr = PAZ*I_ERI_G3yz_S_S_S_M3_vrr+WPZ*I_ERI_G3yz_S_S_S_M4_vrr+oned2z*I_ERI_F3y_S_S_S_M3_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_Hy4z_S_S_S_M3_vrr = PAY*I_ERI_G4z_S_S_S_M3_vrr+WPY*I_ERI_G4z_S_S_S_M4_vrr;
      Double I_ERI_H5z_S_S_S_M3_vrr = PAZ*I_ERI_G4z_S_S_S_M3_vrr+WPZ*I_ERI_G4z_S_S_S_M4_vrr+4*oned2z*I_ERI_F3z_S_S_S_M3_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 10 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M4
       ************************************************************/
      Double I_ERI_I6x_S_S_S_M3_vrr = PAX*I_ERI_H5x_S_S_S_M3_vrr+WPX*I_ERI_H5x_S_S_S_M4_vrr+5*oned2z*I_ERI_G4x_S_S_S_M3_vrr-5*rhod2zsq*I_ERI_G4x_S_S_S_M4_vrr;
      Double I_ERI_I5xy_S_S_S_M3_vrr = PAY*I_ERI_H5x_S_S_S_M3_vrr+WPY*I_ERI_H5x_S_S_S_M4_vrr;
      Double I_ERI_I5xz_S_S_S_M3_vrr = PAZ*I_ERI_H5x_S_S_S_M3_vrr+WPZ*I_ERI_H5x_S_S_S_M4_vrr;
      Double I_ERI_I4x2y_S_S_S_M3_vrr = PAY*I_ERI_H4xy_S_S_S_M3_vrr+WPY*I_ERI_H4xy_S_S_S_M4_vrr+oned2z*I_ERI_G4x_S_S_S_M3_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M4_vrr;
      Double I_ERI_I4x2z_S_S_S_M3_vrr = PAZ*I_ERI_H4xz_S_S_S_M3_vrr+WPZ*I_ERI_H4xz_S_S_S_M4_vrr+oned2z*I_ERI_G4x_S_S_S_M3_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M4_vrr;
      Double I_ERI_I3x3y_S_S_S_M3_vrr = PAY*I_ERI_H3x2y_S_S_S_M3_vrr+WPY*I_ERI_H3x2y_S_S_S_M4_vrr+2*oned2z*I_ERI_G3xy_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_G3xy_S_S_S_M4_vrr;
      Double I_ERI_I3x3z_S_S_S_M3_vrr = PAZ*I_ERI_H3x2z_S_S_S_M3_vrr+WPZ*I_ERI_H3x2z_S_S_S_M4_vrr+2*oned2z*I_ERI_G3xz_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_G3xz_S_S_S_M4_vrr;
      Double I_ERI_I2x4y_S_S_S_M3_vrr = PAX*I_ERI_Hx4y_S_S_S_M3_vrr+WPX*I_ERI_Hx4y_S_S_S_M4_vrr+oned2z*I_ERI_G4y_S_S_S_M3_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M4_vrr;
      Double I_ERI_I2x4z_S_S_S_M3_vrr = PAX*I_ERI_Hx4z_S_S_S_M3_vrr+WPX*I_ERI_Hx4z_S_S_S_M4_vrr+oned2z*I_ERI_G4z_S_S_S_M3_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M4_vrr;
      Double I_ERI_Ix5y_S_S_S_M3_vrr = PAX*I_ERI_H5y_S_S_S_M3_vrr+WPX*I_ERI_H5y_S_S_S_M4_vrr;
      Double I_ERI_Ix5z_S_S_S_M3_vrr = PAX*I_ERI_H5z_S_S_S_M3_vrr+WPX*I_ERI_H5z_S_S_S_M4_vrr;
      Double I_ERI_I6y_S_S_S_M3_vrr = PAY*I_ERI_H5y_S_S_S_M3_vrr+WPY*I_ERI_H5y_S_S_S_M4_vrr+5*oned2z*I_ERI_G4y_S_S_S_M3_vrr-5*rhod2zsq*I_ERI_G4y_S_S_S_M4_vrr;
      Double I_ERI_I5yz_S_S_S_M3_vrr = PAZ*I_ERI_H5y_S_S_S_M3_vrr+WPZ*I_ERI_H5y_S_S_S_M4_vrr;
      Double I_ERI_I4y2z_S_S_S_M3_vrr = PAZ*I_ERI_H4yz_S_S_S_M3_vrr+WPZ*I_ERI_H4yz_S_S_S_M4_vrr+oned2z*I_ERI_G4y_S_S_S_M3_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M4_vrr;
      Double I_ERI_I3y3z_S_S_S_M3_vrr = PAZ*I_ERI_H3y2z_S_S_S_M3_vrr+WPZ*I_ERI_H3y2z_S_S_S_M4_vrr+2*oned2z*I_ERI_G3yz_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_G3yz_S_S_S_M4_vrr;
      Double I_ERI_I2y4z_S_S_S_M3_vrr = PAY*I_ERI_Hy4z_S_S_S_M3_vrr+WPY*I_ERI_Hy4z_S_S_S_M4_vrr+oned2z*I_ERI_G4z_S_S_S_M3_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M4_vrr;
      Double I_ERI_Iy5z_S_S_S_M3_vrr = PAY*I_ERI_H5z_S_S_S_M3_vrr+WPY*I_ERI_H5z_S_S_S_M4_vrr;
      Double I_ERI_I6z_S_S_S_M3_vrr = PAZ*I_ERI_H5z_S_S_S_M3_vrr+WPZ*I_ERI_H5z_S_S_S_M4_vrr+5*oned2z*I_ERI_G4z_S_S_S_M3_vrr-5*rhod2zsq*I_ERI_G4z_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_H_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 5 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       ************************************************************/
      Double I_ERI_H5x_S_S_S_M2_vrr = PAX*I_ERI_G4x_S_S_S_M2_vrr+WPX*I_ERI_G4x_S_S_S_M3_vrr+4*oned2z*I_ERI_F3x_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H4xy_S_S_S_M2_vrr = PAY*I_ERI_G4x_S_S_S_M2_vrr+WPY*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_H4xz_S_S_S_M2_vrr = PAZ*I_ERI_G4x_S_S_S_M2_vrr+WPZ*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_H3x2y_S_S_S_M2_vrr = PAY*I_ERI_G3xy_S_S_S_M2_vrr+WPY*I_ERI_G3xy_S_S_S_M3_vrr+oned2z*I_ERI_F3x_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H3x2z_S_S_S_M2_vrr = PAZ*I_ERI_G3xz_S_S_S_M2_vrr+WPZ*I_ERI_G3xz_S_S_S_M3_vrr+oned2z*I_ERI_F3x_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H2x3y_S_S_S_M2_vrr = PAX*I_ERI_Gx3y_S_S_S_M2_vrr+WPX*I_ERI_Gx3y_S_S_S_M3_vrr+oned2z*I_ERI_F3y_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H2x2yz_S_S_S_M2_vrr = PAZ*I_ERI_G2x2y_S_S_S_M2_vrr+WPZ*I_ERI_G2x2y_S_S_S_M3_vrr;
      Double I_ERI_H2x3z_S_S_S_M2_vrr = PAX*I_ERI_Gx3z_S_S_S_M2_vrr+WPX*I_ERI_Gx3z_S_S_S_M3_vrr+oned2z*I_ERI_F3z_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_Hx4y_S_S_S_M2_vrr = PAX*I_ERI_G4y_S_S_S_M2_vrr+WPX*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_Hx4z_S_S_S_M2_vrr = PAX*I_ERI_G4z_S_S_S_M2_vrr+WPX*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_H5y_S_S_S_M2_vrr = PAY*I_ERI_G4y_S_S_S_M2_vrr+WPY*I_ERI_G4y_S_S_S_M3_vrr+4*oned2z*I_ERI_F3y_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H4yz_S_S_S_M2_vrr = PAZ*I_ERI_G4y_S_S_S_M2_vrr+WPZ*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_H3y2z_S_S_S_M2_vrr = PAZ*I_ERI_G3yz_S_S_S_M2_vrr+WPZ*I_ERI_G3yz_S_S_S_M3_vrr+oned2z*I_ERI_F3y_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H2y3z_S_S_S_M2_vrr = PAY*I_ERI_Gy3z_S_S_S_M2_vrr+WPY*I_ERI_Gy3z_S_S_S_M3_vrr+oned2z*I_ERI_F3z_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_Hy4z_S_S_S_M2_vrr = PAY*I_ERI_G4z_S_S_S_M2_vrr+WPY*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_H5z_S_S_S_M2_vrr = PAZ*I_ERI_G4z_S_S_S_M2_vrr+WPZ*I_ERI_G4z_S_S_S_M3_vrr+4*oned2z*I_ERI_F3z_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       ************************************************************/
      Double I_ERI_I6x_S_S_S_M2_vrr = PAX*I_ERI_H5x_S_S_S_M2_vrr+WPX*I_ERI_H5x_S_S_S_M3_vrr+5*oned2z*I_ERI_G4x_S_S_S_M2_vrr-5*rhod2zsq*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_I5xy_S_S_S_M2_vrr = PAY*I_ERI_H5x_S_S_S_M2_vrr+WPY*I_ERI_H5x_S_S_S_M3_vrr;
      Double I_ERI_I5xz_S_S_S_M2_vrr = PAZ*I_ERI_H5x_S_S_S_M2_vrr+WPZ*I_ERI_H5x_S_S_S_M3_vrr;
      Double I_ERI_I4x2y_S_S_S_M2_vrr = PAY*I_ERI_H4xy_S_S_S_M2_vrr+WPY*I_ERI_H4xy_S_S_S_M3_vrr+oned2z*I_ERI_G4x_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_I4x2z_S_S_S_M2_vrr = PAZ*I_ERI_H4xz_S_S_S_M2_vrr+WPZ*I_ERI_H4xz_S_S_S_M3_vrr+oned2z*I_ERI_G4x_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_I3x3y_S_S_S_M2_vrr = PAY*I_ERI_H3x2y_S_S_S_M2_vrr+WPY*I_ERI_H3x2y_S_S_S_M3_vrr+2*oned2z*I_ERI_G3xy_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_G3xy_S_S_S_M3_vrr;
      Double I_ERI_I3x2yz_S_S_S_M2_vrr = PAZ*I_ERI_H3x2y_S_S_S_M2_vrr+WPZ*I_ERI_H3x2y_S_S_S_M3_vrr;
      Double I_ERI_I3x3z_S_S_S_M2_vrr = PAZ*I_ERI_H3x2z_S_S_S_M2_vrr+WPZ*I_ERI_H3x2z_S_S_S_M3_vrr+2*oned2z*I_ERI_G3xz_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_G3xz_S_S_S_M3_vrr;
      Double I_ERI_I2x4y_S_S_S_M2_vrr = PAX*I_ERI_Hx4y_S_S_S_M2_vrr+WPX*I_ERI_Hx4y_S_S_S_M3_vrr+oned2z*I_ERI_G4y_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_I2x3yz_S_S_S_M2_vrr = PAZ*I_ERI_H2x3y_S_S_S_M2_vrr+WPZ*I_ERI_H2x3y_S_S_S_M3_vrr;
      Double I_ERI_I2xy3z_S_S_S_M2_vrr = PAY*I_ERI_H2x3z_S_S_S_M2_vrr+WPY*I_ERI_H2x3z_S_S_S_M3_vrr;
      Double I_ERI_I2x4z_S_S_S_M2_vrr = PAX*I_ERI_Hx4z_S_S_S_M2_vrr+WPX*I_ERI_Hx4z_S_S_S_M3_vrr+oned2z*I_ERI_G4z_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_Ix5y_S_S_S_M2_vrr = PAX*I_ERI_H5y_S_S_S_M2_vrr+WPX*I_ERI_H5y_S_S_S_M3_vrr;
      Double I_ERI_Ix5z_S_S_S_M2_vrr = PAX*I_ERI_H5z_S_S_S_M2_vrr+WPX*I_ERI_H5z_S_S_S_M3_vrr;
      Double I_ERI_I6y_S_S_S_M2_vrr = PAY*I_ERI_H5y_S_S_S_M2_vrr+WPY*I_ERI_H5y_S_S_S_M3_vrr+5*oned2z*I_ERI_G4y_S_S_S_M2_vrr-5*rhod2zsq*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_I5yz_S_S_S_M2_vrr = PAZ*I_ERI_H5y_S_S_S_M2_vrr+WPZ*I_ERI_H5y_S_S_S_M3_vrr;
      Double I_ERI_I4y2z_S_S_S_M2_vrr = PAZ*I_ERI_H4yz_S_S_S_M2_vrr+WPZ*I_ERI_H4yz_S_S_S_M3_vrr+oned2z*I_ERI_G4y_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_I3y3z_S_S_S_M2_vrr = PAZ*I_ERI_H3y2z_S_S_S_M2_vrr+WPZ*I_ERI_H3y2z_S_S_S_M3_vrr+2*oned2z*I_ERI_G3yz_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_G3yz_S_S_S_M3_vrr;
      Double I_ERI_I2y4z_S_S_S_M2_vrr = PAY*I_ERI_Hy4z_S_S_S_M2_vrr+WPY*I_ERI_Hy4z_S_S_S_M3_vrr+oned2z*I_ERI_G4z_S_S_S_M2_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_Iy5z_S_S_S_M2_vrr = PAY*I_ERI_H5z_S_S_S_M2_vrr+WPY*I_ERI_H5z_S_S_S_M3_vrr;
      Double I_ERI_I6z_S_S_S_M2_vrr = PAZ*I_ERI_H5z_S_S_S_M2_vrr+WPZ*I_ERI_H5z_S_S_S_M3_vrr+5*oned2z*I_ERI_G4z_S_S_S_M2_vrr-5*rhod2zsq*I_ERI_G4z_S_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 9 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M3
       ************************************************************/
      Double I_ERI_K7x_S_S_S_M2_vrr = PAX*I_ERI_I6x_S_S_S_M2_vrr+WPX*I_ERI_I6x_S_S_S_M3_vrr+6*oned2z*I_ERI_H5x_S_S_S_M2_vrr-6*rhod2zsq*I_ERI_H5x_S_S_S_M3_vrr;
      Double I_ERI_K6xy_S_S_S_M2_vrr = PAY*I_ERI_I6x_S_S_S_M2_vrr+WPY*I_ERI_I6x_S_S_S_M3_vrr;
      Double I_ERI_K6xz_S_S_S_M2_vrr = PAZ*I_ERI_I6x_S_S_S_M2_vrr+WPZ*I_ERI_I6x_S_S_S_M3_vrr;
      Double I_ERI_K5x2y_S_S_S_M2_vrr = PAY*I_ERI_I5xy_S_S_S_M2_vrr+WPY*I_ERI_I5xy_S_S_S_M3_vrr+oned2z*I_ERI_H5x_S_S_S_M2_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M3_vrr;
      Double I_ERI_K5x2z_S_S_S_M2_vrr = PAZ*I_ERI_I5xz_S_S_S_M2_vrr+WPZ*I_ERI_I5xz_S_S_S_M3_vrr+oned2z*I_ERI_H5x_S_S_S_M2_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M3_vrr;
      Double I_ERI_K4x3y_S_S_S_M2_vrr = PAY*I_ERI_I4x2y_S_S_S_M2_vrr+WPY*I_ERI_I4x2y_S_S_S_M3_vrr+2*oned2z*I_ERI_H4xy_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_H4xy_S_S_S_M3_vrr;
      Double I_ERI_K4x2yz_S_S_S_M2_vrr = PAZ*I_ERI_I4x2y_S_S_S_M2_vrr+WPZ*I_ERI_I4x2y_S_S_S_M3_vrr;
      Double I_ERI_K4x3z_S_S_S_M2_vrr = PAZ*I_ERI_I4x2z_S_S_S_M2_vrr+WPZ*I_ERI_I4x2z_S_S_S_M3_vrr+2*oned2z*I_ERI_H4xz_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_H4xz_S_S_S_M3_vrr;
      Double I_ERI_K3x4y_S_S_S_M2_vrr = PAX*I_ERI_I2x4y_S_S_S_M2_vrr+WPX*I_ERI_I2x4y_S_S_S_M3_vrr+2*oned2z*I_ERI_Hx4y_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Hx4y_S_S_S_M3_vrr;
      Double I_ERI_K3x3yz_S_S_S_M2_vrr = PAZ*I_ERI_I3x3y_S_S_S_M2_vrr+WPZ*I_ERI_I3x3y_S_S_S_M3_vrr;
      Double I_ERI_K3xy3z_S_S_S_M2_vrr = PAY*I_ERI_I3x3z_S_S_S_M2_vrr+WPY*I_ERI_I3x3z_S_S_S_M3_vrr;
      Double I_ERI_K3x4z_S_S_S_M2_vrr = PAX*I_ERI_I2x4z_S_S_S_M2_vrr+WPX*I_ERI_I2x4z_S_S_S_M3_vrr+2*oned2z*I_ERI_Hx4z_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Hx4z_S_S_S_M3_vrr;
      Double I_ERI_K2x5y_S_S_S_M2_vrr = PAX*I_ERI_Ix5y_S_S_S_M2_vrr+WPX*I_ERI_Ix5y_S_S_S_M3_vrr+oned2z*I_ERI_H5y_S_S_S_M2_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M3_vrr;
      Double I_ERI_K2x4yz_S_S_S_M2_vrr = PAZ*I_ERI_I2x4y_S_S_S_M2_vrr+WPZ*I_ERI_I2x4y_S_S_S_M3_vrr;
      Double I_ERI_K2xy4z_S_S_S_M2_vrr = PAY*I_ERI_I2x4z_S_S_S_M2_vrr+WPY*I_ERI_I2x4z_S_S_S_M3_vrr;
      Double I_ERI_K2x5z_S_S_S_M2_vrr = PAX*I_ERI_Ix5z_S_S_S_M2_vrr+WPX*I_ERI_Ix5z_S_S_S_M3_vrr+oned2z*I_ERI_H5z_S_S_S_M2_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M3_vrr;
      Double I_ERI_Kx6y_S_S_S_M2_vrr = PAX*I_ERI_I6y_S_S_S_M2_vrr+WPX*I_ERI_I6y_S_S_S_M3_vrr;
      Double I_ERI_Kx3y3z_S_S_S_M2_vrr = PAX*I_ERI_I3y3z_S_S_S_M2_vrr+WPX*I_ERI_I3y3z_S_S_S_M3_vrr;
      Double I_ERI_Kx6z_S_S_S_M2_vrr = PAX*I_ERI_I6z_S_S_S_M2_vrr+WPX*I_ERI_I6z_S_S_S_M3_vrr;
      Double I_ERI_K7y_S_S_S_M2_vrr = PAY*I_ERI_I6y_S_S_S_M2_vrr+WPY*I_ERI_I6y_S_S_S_M3_vrr+6*oned2z*I_ERI_H5y_S_S_S_M2_vrr-6*rhod2zsq*I_ERI_H5y_S_S_S_M3_vrr;
      Double I_ERI_K6yz_S_S_S_M2_vrr = PAZ*I_ERI_I6y_S_S_S_M2_vrr+WPZ*I_ERI_I6y_S_S_S_M3_vrr;
      Double I_ERI_K5y2z_S_S_S_M2_vrr = PAZ*I_ERI_I5yz_S_S_S_M2_vrr+WPZ*I_ERI_I5yz_S_S_S_M3_vrr+oned2z*I_ERI_H5y_S_S_S_M2_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M3_vrr;
      Double I_ERI_K4y3z_S_S_S_M2_vrr = PAZ*I_ERI_I4y2z_S_S_S_M2_vrr+WPZ*I_ERI_I4y2z_S_S_S_M3_vrr+2*oned2z*I_ERI_H4yz_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_H4yz_S_S_S_M3_vrr;
      Double I_ERI_K3y4z_S_S_S_M2_vrr = PAY*I_ERI_I2y4z_S_S_S_M2_vrr+WPY*I_ERI_I2y4z_S_S_S_M3_vrr+2*oned2z*I_ERI_Hy4z_S_S_S_M2_vrr-2*rhod2zsq*I_ERI_Hy4z_S_S_S_M3_vrr;
      Double I_ERI_K2y5z_S_S_S_M2_vrr = PAY*I_ERI_Iy5z_S_S_S_M2_vrr+WPY*I_ERI_Iy5z_S_S_S_M3_vrr+oned2z*I_ERI_H5z_S_S_S_M2_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M3_vrr;
      Double I_ERI_Ky6z_S_S_S_M2_vrr = PAY*I_ERI_I6z_S_S_S_M2_vrr+WPY*I_ERI_I6z_S_S_S_M3_vrr;
      Double I_ERI_K7z_S_S_S_M2_vrr = PAZ*I_ERI_I6z_S_S_S_M2_vrr+WPZ*I_ERI_I6z_S_S_S_M3_vrr+6*oned2z*I_ERI_H5z_S_S_S_M2_vrr-6*rhod2zsq*I_ERI_H5z_S_S_S_M3_vrr;

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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M2
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M1_vrr = PAX*I_ERI_Px_S_S_S_M1_vrr+WPX*I_ERI_Px_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
      Double I_ERI_Dxy_S_S_S_M1_vrr = PAY*I_ERI_Px_S_S_S_M1_vrr+WPY*I_ERI_Px_S_S_S_M2_vrr;
      Double I_ERI_D2y_S_S_S_M1_vrr = PAY*I_ERI_Py_S_S_S_M1_vrr+WPY*I_ERI_Py_S_S_S_M2_vrr+oned2z*I_ERI_S_S_S_S_M1_vrr-rhod2zsq*I_ERI_S_S_S_S_M2_vrr;
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
       * shell quartet name: SQ_ERI_I_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       ************************************************************/
      Double I_ERI_I6x_S_S_S_M1_vrr = PAX*I_ERI_H5x_S_S_S_M1_vrr+WPX*I_ERI_H5x_S_S_S_M2_vrr+5*oned2z*I_ERI_G4x_S_S_S_M1_vrr-5*rhod2zsq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_I5xy_S_S_S_M1_vrr = PAY*I_ERI_H5x_S_S_S_M1_vrr+WPY*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_I5xz_S_S_S_M1_vrr = PAZ*I_ERI_H5x_S_S_S_M1_vrr+WPZ*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_I4x2y_S_S_S_M1_vrr = PAY*I_ERI_H4xy_S_S_S_M1_vrr+WPY*I_ERI_H4xy_S_S_S_M2_vrr+oned2z*I_ERI_G4x_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_I4xyz_S_S_S_M1_vrr = PAZ*I_ERI_H4xy_S_S_S_M1_vrr+WPZ*I_ERI_H4xy_S_S_S_M2_vrr;
      Double I_ERI_I4x2z_S_S_S_M1_vrr = PAZ*I_ERI_H4xz_S_S_S_M1_vrr+WPZ*I_ERI_H4xz_S_S_S_M2_vrr+oned2z*I_ERI_G4x_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_I3x3y_S_S_S_M1_vrr = PAY*I_ERI_H3x2y_S_S_S_M1_vrr+WPY*I_ERI_H3x2y_S_S_S_M2_vrr+2*oned2z*I_ERI_G3xy_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_I3x2yz_S_S_S_M1_vrr = PAZ*I_ERI_H3x2y_S_S_S_M1_vrr+WPZ*I_ERI_H3x2y_S_S_S_M2_vrr;
      Double I_ERI_I3xy2z_S_S_S_M1_vrr = PAY*I_ERI_H3x2z_S_S_S_M1_vrr+WPY*I_ERI_H3x2z_S_S_S_M2_vrr;
      Double I_ERI_I3x3z_S_S_S_M1_vrr = PAZ*I_ERI_H3x2z_S_S_S_M1_vrr+WPZ*I_ERI_H3x2z_S_S_S_M2_vrr+2*oned2z*I_ERI_G3xz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_I2x4y_S_S_S_M1_vrr = PAX*I_ERI_Hx4y_S_S_S_M1_vrr+WPX*I_ERI_Hx4y_S_S_S_M2_vrr+oned2z*I_ERI_G4y_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_I2x3yz_S_S_S_M1_vrr = PAZ*I_ERI_H2x3y_S_S_S_M1_vrr+WPZ*I_ERI_H2x3y_S_S_S_M2_vrr;
      Double I_ERI_I2x2y2z_S_S_S_M1_vrr = PAZ*I_ERI_H2x2yz_S_S_S_M1_vrr+WPZ*I_ERI_H2x2yz_S_S_S_M2_vrr+oned2z*I_ERI_G2x2y_S_S_S_M1_vrr-rhod2zsq*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_I2xy3z_S_S_S_M1_vrr = PAY*I_ERI_H2x3z_S_S_S_M1_vrr+WPY*I_ERI_H2x3z_S_S_S_M2_vrr;
      Double I_ERI_I2x4z_S_S_S_M1_vrr = PAX*I_ERI_Hx4z_S_S_S_M1_vrr+WPX*I_ERI_Hx4z_S_S_S_M2_vrr+oned2z*I_ERI_G4z_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_Ix5y_S_S_S_M1_vrr = PAX*I_ERI_H5y_S_S_S_M1_vrr+WPX*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_Ix4yz_S_S_S_M1_vrr = PAZ*I_ERI_Hx4y_S_S_S_M1_vrr+WPZ*I_ERI_Hx4y_S_S_S_M2_vrr;
      Double I_ERI_Ix3y2z_S_S_S_M1_vrr = PAX*I_ERI_H3y2z_S_S_S_M1_vrr+WPX*I_ERI_H3y2z_S_S_S_M2_vrr;
      Double I_ERI_Ix2y3z_S_S_S_M1_vrr = PAX*I_ERI_H2y3z_S_S_S_M1_vrr+WPX*I_ERI_H2y3z_S_S_S_M2_vrr;
      Double I_ERI_Ixy4z_S_S_S_M1_vrr = PAY*I_ERI_Hx4z_S_S_S_M1_vrr+WPY*I_ERI_Hx4z_S_S_S_M2_vrr;
      Double I_ERI_Ix5z_S_S_S_M1_vrr = PAX*I_ERI_H5z_S_S_S_M1_vrr+WPX*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_I6y_S_S_S_M1_vrr = PAY*I_ERI_H5y_S_S_S_M1_vrr+WPY*I_ERI_H5y_S_S_S_M2_vrr+5*oned2z*I_ERI_G4y_S_S_S_M1_vrr-5*rhod2zsq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_I5yz_S_S_S_M1_vrr = PAZ*I_ERI_H5y_S_S_S_M1_vrr+WPZ*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_I4y2z_S_S_S_M1_vrr = PAZ*I_ERI_H4yz_S_S_S_M1_vrr+WPZ*I_ERI_H4yz_S_S_S_M2_vrr+oned2z*I_ERI_G4y_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_I3y3z_S_S_S_M1_vrr = PAZ*I_ERI_H3y2z_S_S_S_M1_vrr+WPZ*I_ERI_H3y2z_S_S_S_M2_vrr+2*oned2z*I_ERI_G3yz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_I2y4z_S_S_S_M1_vrr = PAY*I_ERI_Hy4z_S_S_S_M1_vrr+WPY*I_ERI_Hy4z_S_S_S_M2_vrr+oned2z*I_ERI_G4z_S_S_S_M1_vrr-rhod2zsq*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_Iy5z_S_S_S_M1_vrr = PAY*I_ERI_H5z_S_S_S_M1_vrr+WPY*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_I6z_S_S_S_M1_vrr = PAZ*I_ERI_H5z_S_S_S_M1_vrr+WPZ*I_ERI_H5z_S_S_S_M2_vrr+5*oned2z*I_ERI_G4z_S_S_S_M1_vrr-5*rhod2zsq*I_ERI_G4z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       ************************************************************/
      Double I_ERI_K7x_S_S_S_M1_vrr = PAX*I_ERI_I6x_S_S_S_M1_vrr+WPX*I_ERI_I6x_S_S_S_M2_vrr+6*oned2z*I_ERI_H5x_S_S_S_M1_vrr-6*rhod2zsq*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_K6xy_S_S_S_M1_vrr = PAY*I_ERI_I6x_S_S_S_M1_vrr+WPY*I_ERI_I6x_S_S_S_M2_vrr;
      Double I_ERI_K6xz_S_S_S_M1_vrr = PAZ*I_ERI_I6x_S_S_S_M1_vrr+WPZ*I_ERI_I6x_S_S_S_M2_vrr;
      Double I_ERI_K5x2y_S_S_S_M1_vrr = PAY*I_ERI_I5xy_S_S_S_M1_vrr+WPY*I_ERI_I5xy_S_S_S_M2_vrr+oned2z*I_ERI_H5x_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_K5xyz_S_S_S_M1_vrr = PAZ*I_ERI_I5xy_S_S_S_M1_vrr+WPZ*I_ERI_I5xy_S_S_S_M2_vrr;
      Double I_ERI_K5x2z_S_S_S_M1_vrr = PAZ*I_ERI_I5xz_S_S_S_M1_vrr+WPZ*I_ERI_I5xz_S_S_S_M2_vrr+oned2z*I_ERI_H5x_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_K4x3y_S_S_S_M1_vrr = PAY*I_ERI_I4x2y_S_S_S_M1_vrr+WPY*I_ERI_I4x2y_S_S_S_M2_vrr+2*oned2z*I_ERI_H4xy_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_H4xy_S_S_S_M2_vrr;
      Double I_ERI_K4x2yz_S_S_S_M1_vrr = PAZ*I_ERI_I4x2y_S_S_S_M1_vrr+WPZ*I_ERI_I4x2y_S_S_S_M2_vrr;
      Double I_ERI_K4xy2z_S_S_S_M1_vrr = PAY*I_ERI_I4x2z_S_S_S_M1_vrr+WPY*I_ERI_I4x2z_S_S_S_M2_vrr;
      Double I_ERI_K4x3z_S_S_S_M1_vrr = PAZ*I_ERI_I4x2z_S_S_S_M1_vrr+WPZ*I_ERI_I4x2z_S_S_S_M2_vrr+2*oned2z*I_ERI_H4xz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_H4xz_S_S_S_M2_vrr;
      Double I_ERI_K3x4y_S_S_S_M1_vrr = PAX*I_ERI_I2x4y_S_S_S_M1_vrr+WPX*I_ERI_I2x4y_S_S_S_M2_vrr+2*oned2z*I_ERI_Hx4y_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Hx4y_S_S_S_M2_vrr;
      Double I_ERI_K3x3yz_S_S_S_M1_vrr = PAZ*I_ERI_I3x3y_S_S_S_M1_vrr+WPZ*I_ERI_I3x3y_S_S_S_M2_vrr;
      Double I_ERI_K3x2y2z_S_S_S_M1_vrr = PAZ*I_ERI_I3x2yz_S_S_S_M1_vrr+WPZ*I_ERI_I3x2yz_S_S_S_M2_vrr+oned2z*I_ERI_H3x2y_S_S_S_M1_vrr-rhod2zsq*I_ERI_H3x2y_S_S_S_M2_vrr;
      Double I_ERI_K3xy3z_S_S_S_M1_vrr = PAY*I_ERI_I3x3z_S_S_S_M1_vrr+WPY*I_ERI_I3x3z_S_S_S_M2_vrr;
      Double I_ERI_K3x4z_S_S_S_M1_vrr = PAX*I_ERI_I2x4z_S_S_S_M1_vrr+WPX*I_ERI_I2x4z_S_S_S_M2_vrr+2*oned2z*I_ERI_Hx4z_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Hx4z_S_S_S_M2_vrr;
      Double I_ERI_K2x5y_S_S_S_M1_vrr = PAX*I_ERI_Ix5y_S_S_S_M1_vrr+WPX*I_ERI_Ix5y_S_S_S_M2_vrr+oned2z*I_ERI_H5y_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_K2x4yz_S_S_S_M1_vrr = PAZ*I_ERI_I2x4y_S_S_S_M1_vrr+WPZ*I_ERI_I2x4y_S_S_S_M2_vrr;
      Double I_ERI_K2x3y2z_S_S_S_M1_vrr = PAZ*I_ERI_I2x3yz_S_S_S_M1_vrr+WPZ*I_ERI_I2x3yz_S_S_S_M2_vrr+oned2z*I_ERI_H2x3y_S_S_S_M1_vrr-rhod2zsq*I_ERI_H2x3y_S_S_S_M2_vrr;
      Double I_ERI_K2x2y3z_S_S_S_M1_vrr = PAY*I_ERI_I2xy3z_S_S_S_M1_vrr+WPY*I_ERI_I2xy3z_S_S_S_M2_vrr+oned2z*I_ERI_H2x3z_S_S_S_M1_vrr-rhod2zsq*I_ERI_H2x3z_S_S_S_M2_vrr;
      Double I_ERI_K2xy4z_S_S_S_M1_vrr = PAY*I_ERI_I2x4z_S_S_S_M1_vrr+WPY*I_ERI_I2x4z_S_S_S_M2_vrr;
      Double I_ERI_K2x5z_S_S_S_M1_vrr = PAX*I_ERI_Ix5z_S_S_S_M1_vrr+WPX*I_ERI_Ix5z_S_S_S_M2_vrr+oned2z*I_ERI_H5z_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_Kx6y_S_S_S_M1_vrr = PAX*I_ERI_I6y_S_S_S_M1_vrr+WPX*I_ERI_I6y_S_S_S_M2_vrr;
      Double I_ERI_Kx5yz_S_S_S_M1_vrr = PAZ*I_ERI_Ix5y_S_S_S_M1_vrr+WPZ*I_ERI_Ix5y_S_S_S_M2_vrr;
      Double I_ERI_Kx4y2z_S_S_S_M1_vrr = PAX*I_ERI_I4y2z_S_S_S_M1_vrr+WPX*I_ERI_I4y2z_S_S_S_M2_vrr;
      Double I_ERI_Kx3y3z_S_S_S_M1_vrr = PAX*I_ERI_I3y3z_S_S_S_M1_vrr+WPX*I_ERI_I3y3z_S_S_S_M2_vrr;
      Double I_ERI_Kx2y4z_S_S_S_M1_vrr = PAX*I_ERI_I2y4z_S_S_S_M1_vrr+WPX*I_ERI_I2y4z_S_S_S_M2_vrr;
      Double I_ERI_Kxy5z_S_S_S_M1_vrr = PAY*I_ERI_Ix5z_S_S_S_M1_vrr+WPY*I_ERI_Ix5z_S_S_S_M2_vrr;
      Double I_ERI_Kx6z_S_S_S_M1_vrr = PAX*I_ERI_I6z_S_S_S_M1_vrr+WPX*I_ERI_I6z_S_S_S_M2_vrr;
      Double I_ERI_K7y_S_S_S_M1_vrr = PAY*I_ERI_I6y_S_S_S_M1_vrr+WPY*I_ERI_I6y_S_S_S_M2_vrr+6*oned2z*I_ERI_H5y_S_S_S_M1_vrr-6*rhod2zsq*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_K6yz_S_S_S_M1_vrr = PAZ*I_ERI_I6y_S_S_S_M1_vrr+WPZ*I_ERI_I6y_S_S_S_M2_vrr;
      Double I_ERI_K5y2z_S_S_S_M1_vrr = PAZ*I_ERI_I5yz_S_S_S_M1_vrr+WPZ*I_ERI_I5yz_S_S_S_M2_vrr+oned2z*I_ERI_H5y_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_K4y3z_S_S_S_M1_vrr = PAZ*I_ERI_I4y2z_S_S_S_M1_vrr+WPZ*I_ERI_I4y2z_S_S_S_M2_vrr+2*oned2z*I_ERI_H4yz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_H4yz_S_S_S_M2_vrr;
      Double I_ERI_K3y4z_S_S_S_M1_vrr = PAY*I_ERI_I2y4z_S_S_S_M1_vrr+WPY*I_ERI_I2y4z_S_S_S_M2_vrr+2*oned2z*I_ERI_Hy4z_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Hy4z_S_S_S_M2_vrr;
      Double I_ERI_K2y5z_S_S_S_M1_vrr = PAY*I_ERI_Iy5z_S_S_S_M1_vrr+WPY*I_ERI_Iy5z_S_S_S_M2_vrr+oned2z*I_ERI_H5z_S_S_S_M1_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_Ky6z_S_S_S_M1_vrr = PAY*I_ERI_I6z_S_S_S_M1_vrr+WPY*I_ERI_I6z_S_S_S_M2_vrr;
      Double I_ERI_K7z_S_S_S_M1_vrr = PAZ*I_ERI_I6z_S_S_S_M1_vrr+WPZ*I_ERI_I6z_S_S_S_M2_vrr+6*oned2z*I_ERI_H5z_S_S_S_M1_vrr-6*rhod2zsq*I_ERI_H5z_S_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_L_S_S_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_K_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_K_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M2
       ************************************************************/
      Double I_ERI_L8x_S_S_S_M1_vrr = PAX*I_ERI_K7x_S_S_S_M1_vrr+WPX*I_ERI_K7x_S_S_S_M2_vrr+7*oned2z*I_ERI_I6x_S_S_S_M1_vrr-7*rhod2zsq*I_ERI_I6x_S_S_S_M2_vrr;
      Double I_ERI_L7xy_S_S_S_M1_vrr = PAY*I_ERI_K7x_S_S_S_M1_vrr+WPY*I_ERI_K7x_S_S_S_M2_vrr;
      Double I_ERI_L7xz_S_S_S_M1_vrr = PAZ*I_ERI_K7x_S_S_S_M1_vrr+WPZ*I_ERI_K7x_S_S_S_M2_vrr;
      Double I_ERI_L6x2y_S_S_S_M1_vrr = PAY*I_ERI_K6xy_S_S_S_M1_vrr+WPY*I_ERI_K6xy_S_S_S_M2_vrr+oned2z*I_ERI_I6x_S_S_S_M1_vrr-rhod2zsq*I_ERI_I6x_S_S_S_M2_vrr;
      Double I_ERI_L6xyz_S_S_S_M1_vrr = PAZ*I_ERI_K6xy_S_S_S_M1_vrr+WPZ*I_ERI_K6xy_S_S_S_M2_vrr;
      Double I_ERI_L6x2z_S_S_S_M1_vrr = PAZ*I_ERI_K6xz_S_S_S_M1_vrr+WPZ*I_ERI_K6xz_S_S_S_M2_vrr+oned2z*I_ERI_I6x_S_S_S_M1_vrr-rhod2zsq*I_ERI_I6x_S_S_S_M2_vrr;
      Double I_ERI_L5x3y_S_S_S_M1_vrr = PAY*I_ERI_K5x2y_S_S_S_M1_vrr+WPY*I_ERI_K5x2y_S_S_S_M2_vrr+2*oned2z*I_ERI_I5xy_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_I5xy_S_S_S_M2_vrr;
      Double I_ERI_L5x2yz_S_S_S_M1_vrr = PAZ*I_ERI_K5x2y_S_S_S_M1_vrr+WPZ*I_ERI_K5x2y_S_S_S_M2_vrr;
      Double I_ERI_L5xy2z_S_S_S_M1_vrr = PAY*I_ERI_K5x2z_S_S_S_M1_vrr+WPY*I_ERI_K5x2z_S_S_S_M2_vrr;
      Double I_ERI_L5x3z_S_S_S_M1_vrr = PAZ*I_ERI_K5x2z_S_S_S_M1_vrr+WPZ*I_ERI_K5x2z_S_S_S_M2_vrr+2*oned2z*I_ERI_I5xz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_I5xz_S_S_S_M2_vrr;
      Double I_ERI_L4x4y_S_S_S_M1_vrr = PAY*I_ERI_K4x3y_S_S_S_M1_vrr+WPY*I_ERI_K4x3y_S_S_S_M2_vrr+3*oned2z*I_ERI_I4x2y_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_I4x2y_S_S_S_M2_vrr;
      Double I_ERI_L4x3yz_S_S_S_M1_vrr = PAZ*I_ERI_K4x3y_S_S_S_M1_vrr+WPZ*I_ERI_K4x3y_S_S_S_M2_vrr;
      Double I_ERI_L4x2y2z_S_S_S_M1_vrr = PAZ*I_ERI_K4x2yz_S_S_S_M1_vrr+WPZ*I_ERI_K4x2yz_S_S_S_M2_vrr+oned2z*I_ERI_I4x2y_S_S_S_M1_vrr-rhod2zsq*I_ERI_I4x2y_S_S_S_M2_vrr;
      Double I_ERI_L4xy3z_S_S_S_M1_vrr = PAY*I_ERI_K4x3z_S_S_S_M1_vrr+WPY*I_ERI_K4x3z_S_S_S_M2_vrr;
      Double I_ERI_L4x4z_S_S_S_M1_vrr = PAZ*I_ERI_K4x3z_S_S_S_M1_vrr+WPZ*I_ERI_K4x3z_S_S_S_M2_vrr+3*oned2z*I_ERI_I4x2z_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_I4x2z_S_S_S_M2_vrr;
      Double I_ERI_L3x5y_S_S_S_M1_vrr = PAX*I_ERI_K2x5y_S_S_S_M1_vrr+WPX*I_ERI_K2x5y_S_S_S_M2_vrr+2*oned2z*I_ERI_Ix5y_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Ix5y_S_S_S_M2_vrr;
      Double I_ERI_L3x4yz_S_S_S_M1_vrr = PAZ*I_ERI_K3x4y_S_S_S_M1_vrr+WPZ*I_ERI_K3x4y_S_S_S_M2_vrr;
      Double I_ERI_L3x3y2z_S_S_S_M1_vrr = PAZ*I_ERI_K3x3yz_S_S_S_M1_vrr+WPZ*I_ERI_K3x3yz_S_S_S_M2_vrr+oned2z*I_ERI_I3x3y_S_S_S_M1_vrr-rhod2zsq*I_ERI_I3x3y_S_S_S_M2_vrr;
      Double I_ERI_L3x2y3z_S_S_S_M1_vrr = PAY*I_ERI_K3xy3z_S_S_S_M1_vrr+WPY*I_ERI_K3xy3z_S_S_S_M2_vrr+oned2z*I_ERI_I3x3z_S_S_S_M1_vrr-rhod2zsq*I_ERI_I3x3z_S_S_S_M2_vrr;
      Double I_ERI_L3xy4z_S_S_S_M1_vrr = PAY*I_ERI_K3x4z_S_S_S_M1_vrr+WPY*I_ERI_K3x4z_S_S_S_M2_vrr;
      Double I_ERI_L3x5z_S_S_S_M1_vrr = PAX*I_ERI_K2x5z_S_S_S_M1_vrr+WPX*I_ERI_K2x5z_S_S_S_M2_vrr+2*oned2z*I_ERI_Ix5z_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Ix5z_S_S_S_M2_vrr;
      Double I_ERI_L2x6y_S_S_S_M1_vrr = PAX*I_ERI_Kx6y_S_S_S_M1_vrr+WPX*I_ERI_Kx6y_S_S_S_M2_vrr+oned2z*I_ERI_I6y_S_S_S_M1_vrr-rhod2zsq*I_ERI_I6y_S_S_S_M2_vrr;
      Double I_ERI_L2x5yz_S_S_S_M1_vrr = PAZ*I_ERI_K2x5y_S_S_S_M1_vrr+WPZ*I_ERI_K2x5y_S_S_S_M2_vrr;
      Double I_ERI_L2x4y2z_S_S_S_M1_vrr = PAZ*I_ERI_K2x4yz_S_S_S_M1_vrr+WPZ*I_ERI_K2x4yz_S_S_S_M2_vrr+oned2z*I_ERI_I2x4y_S_S_S_M1_vrr-rhod2zsq*I_ERI_I2x4y_S_S_S_M2_vrr;
      Double I_ERI_L2x3y3z_S_S_S_M1_vrr = PAX*I_ERI_Kx3y3z_S_S_S_M1_vrr+WPX*I_ERI_Kx3y3z_S_S_S_M2_vrr+oned2z*I_ERI_I3y3z_S_S_S_M1_vrr-rhod2zsq*I_ERI_I3y3z_S_S_S_M2_vrr;
      Double I_ERI_L2x2y4z_S_S_S_M1_vrr = PAY*I_ERI_K2xy4z_S_S_S_M1_vrr+WPY*I_ERI_K2xy4z_S_S_S_M2_vrr+oned2z*I_ERI_I2x4z_S_S_S_M1_vrr-rhod2zsq*I_ERI_I2x4z_S_S_S_M2_vrr;
      Double I_ERI_L2xy5z_S_S_S_M1_vrr = PAY*I_ERI_K2x5z_S_S_S_M1_vrr+WPY*I_ERI_K2x5z_S_S_S_M2_vrr;
      Double I_ERI_L2x6z_S_S_S_M1_vrr = PAX*I_ERI_Kx6z_S_S_S_M1_vrr+WPX*I_ERI_Kx6z_S_S_S_M2_vrr+oned2z*I_ERI_I6z_S_S_S_M1_vrr-rhod2zsq*I_ERI_I6z_S_S_S_M2_vrr;
      Double I_ERI_Lx7y_S_S_S_M1_vrr = PAX*I_ERI_K7y_S_S_S_M1_vrr+WPX*I_ERI_K7y_S_S_S_M2_vrr;
      Double I_ERI_Lx6yz_S_S_S_M1_vrr = PAZ*I_ERI_Kx6y_S_S_S_M1_vrr+WPZ*I_ERI_Kx6y_S_S_S_M2_vrr;
      Double I_ERI_Lx5y2z_S_S_S_M1_vrr = PAX*I_ERI_K5y2z_S_S_S_M1_vrr+WPX*I_ERI_K5y2z_S_S_S_M2_vrr;
      Double I_ERI_Lx4y3z_S_S_S_M1_vrr = PAX*I_ERI_K4y3z_S_S_S_M1_vrr+WPX*I_ERI_K4y3z_S_S_S_M2_vrr;
      Double I_ERI_Lx3y4z_S_S_S_M1_vrr = PAX*I_ERI_K3y4z_S_S_S_M1_vrr+WPX*I_ERI_K3y4z_S_S_S_M2_vrr;
      Double I_ERI_Lx2y5z_S_S_S_M1_vrr = PAX*I_ERI_K2y5z_S_S_S_M1_vrr+WPX*I_ERI_K2y5z_S_S_S_M2_vrr;
      Double I_ERI_Lxy6z_S_S_S_M1_vrr = PAY*I_ERI_Kx6z_S_S_S_M1_vrr+WPY*I_ERI_Kx6z_S_S_S_M2_vrr;
      Double I_ERI_Lx7z_S_S_S_M1_vrr = PAX*I_ERI_K7z_S_S_S_M1_vrr+WPX*I_ERI_K7z_S_S_S_M2_vrr;
      Double I_ERI_L8y_S_S_S_M1_vrr = PAY*I_ERI_K7y_S_S_S_M1_vrr+WPY*I_ERI_K7y_S_S_S_M2_vrr+7*oned2z*I_ERI_I6y_S_S_S_M1_vrr-7*rhod2zsq*I_ERI_I6y_S_S_S_M2_vrr;
      Double I_ERI_L7yz_S_S_S_M1_vrr = PAZ*I_ERI_K7y_S_S_S_M1_vrr+WPZ*I_ERI_K7y_S_S_S_M2_vrr;
      Double I_ERI_L6y2z_S_S_S_M1_vrr = PAZ*I_ERI_K6yz_S_S_S_M1_vrr+WPZ*I_ERI_K6yz_S_S_S_M2_vrr+oned2z*I_ERI_I6y_S_S_S_M1_vrr-rhod2zsq*I_ERI_I6y_S_S_S_M2_vrr;
      Double I_ERI_L5y3z_S_S_S_M1_vrr = PAZ*I_ERI_K5y2z_S_S_S_M1_vrr+WPZ*I_ERI_K5y2z_S_S_S_M2_vrr+2*oned2z*I_ERI_I5yz_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_I5yz_S_S_S_M2_vrr;
      Double I_ERI_L4y4z_S_S_S_M1_vrr = PAZ*I_ERI_K4y3z_S_S_S_M1_vrr+WPZ*I_ERI_K4y3z_S_S_S_M2_vrr+3*oned2z*I_ERI_I4y2z_S_S_S_M1_vrr-3*rhod2zsq*I_ERI_I4y2z_S_S_S_M2_vrr;
      Double I_ERI_L3y5z_S_S_S_M1_vrr = PAY*I_ERI_K2y5z_S_S_S_M1_vrr+WPY*I_ERI_K2y5z_S_S_S_M2_vrr+2*oned2z*I_ERI_Iy5z_S_S_S_M1_vrr-2*rhod2zsq*I_ERI_Iy5z_S_S_S_M2_vrr;
      Double I_ERI_L2y6z_S_S_S_M1_vrr = PAY*I_ERI_Ky6z_S_S_S_M1_vrr+WPY*I_ERI_Ky6z_S_S_S_M2_vrr+oned2z*I_ERI_I6z_S_S_S_M1_vrr-rhod2zsq*I_ERI_I6z_S_S_S_M2_vrr;
      Double I_ERI_Ly7z_S_S_S_M1_vrr = PAY*I_ERI_K7z_S_S_S_M1_vrr+WPY*I_ERI_K7z_S_S_S_M2_vrr;
      Double I_ERI_L8z_S_S_S_M1_vrr = PAZ*I_ERI_K7z_S_S_S_M1_vrr+WPZ*I_ERI_K7z_S_S_S_M2_vrr+7*oned2z*I_ERI_I6z_S_S_S_M1_vrr-7*rhod2zsq*I_ERI_I6z_S_S_S_M2_vrr;

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
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_S_S_S_S
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M1
       ************************************************************/
      Double I_ERI_D2x_S_S_S_vrr = PAX*I_ERI_Px_S_S_S_vrr+WPX*I_ERI_Px_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_D2y_S_S_S_vrr = PAY*I_ERI_Py_S_S_S_vrr+WPY*I_ERI_Py_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;
      Double I_ERI_D2z_S_S_S_vrr = PAZ*I_ERI_Pz_S_S_S_vrr+WPZ*I_ERI_Pz_S_S_S_M1_vrr+oned2z*I_ERI_S_S_S_S_vrr-rhod2zsq*I_ERI_S_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_P_S_S_S
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M1
       ************************************************************/
      Double I_ERI_F3x_S_S_S_vrr = PAX*I_ERI_D2x_S_S_S_vrr+WPX*I_ERI_D2x_S_S_S_M1_vrr+2*oned2z*I_ERI_Px_S_S_S_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M1_vrr;
      Double I_ERI_F2xy_S_S_S_vrr = PAY*I_ERI_D2x_S_S_S_vrr+WPY*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_F2xz_S_S_S_vrr = PAZ*I_ERI_D2x_S_S_S_vrr+WPZ*I_ERI_D2x_S_S_S_M1_vrr;
      Double I_ERI_Fx2y_S_S_S_vrr = PAX*I_ERI_D2y_S_S_S_vrr+WPX*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_Fx2z_S_S_S_vrr = PAX*I_ERI_D2z_S_S_S_vrr+WPX*I_ERI_D2z_S_S_S_M1_vrr;
      Double I_ERI_F3y_S_S_S_vrr = PAY*I_ERI_D2y_S_S_S_vrr+WPY*I_ERI_D2y_S_S_S_M1_vrr+2*oned2z*I_ERI_Py_S_S_S_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M1_vrr;
      Double I_ERI_F2yz_S_S_S_vrr = PAZ*I_ERI_D2y_S_S_S_vrr+WPZ*I_ERI_D2y_S_S_S_M1_vrr;
      Double I_ERI_F3z_S_S_S_vrr = PAZ*I_ERI_D2z_S_S_S_vrr+WPZ*I_ERI_D2z_S_S_S_M1_vrr+2*oned2z*I_ERI_Pz_S_S_S_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M1_vrr;

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
       * shell quartet name: SQ_ERI_I_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_I_S_S_S
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       ************************************************************/
      Double I_ERI_I6x_S_Px_S_vrr = QCX*I_ERI_I6x_S_S_S_vrr+WQX*I_ERI_I6x_S_S_S_M1_vrr+6*oned2k*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I5xy_S_Px_S_vrr = QCX*I_ERI_I5xy_S_S_S_vrr+WQX*I_ERI_I5xy_S_S_S_M1_vrr+5*oned2k*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_I5xz_S_Px_S_vrr = QCX*I_ERI_I5xz_S_S_S_vrr+WQX*I_ERI_I5xz_S_S_S_M1_vrr+5*oned2k*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_I4x2y_S_Px_S_vrr = QCX*I_ERI_I4x2y_S_S_S_vrr+WQX*I_ERI_I4x2y_S_S_S_M1_vrr+4*oned2k*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_I4xyz_S_Px_S_vrr = QCX*I_ERI_I4xyz_S_S_S_vrr+WQX*I_ERI_I4xyz_S_S_S_M1_vrr+4*oned2k*I_ERI_H3xyz_S_S_S_M1_vrr;
      Double I_ERI_I4x2z_S_Px_S_vrr = QCX*I_ERI_I4x2z_S_S_S_vrr+WQX*I_ERI_I4x2z_S_S_S_M1_vrr+4*oned2k*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_I3x3y_S_Px_S_vrr = QCX*I_ERI_I3x3y_S_S_S_vrr+WQX*I_ERI_I3x3y_S_S_S_M1_vrr+3*oned2k*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Px_S_vrr = QCX*I_ERI_I3x2yz_S_S_S_vrr+WQX*I_ERI_I3x2yz_S_S_S_M1_vrr+3*oned2k*I_ERI_H2x2yz_S_S_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Px_S_vrr = QCX*I_ERI_I3xy2z_S_S_S_vrr+WQX*I_ERI_I3xy2z_S_S_S_M1_vrr+3*oned2k*I_ERI_H2xy2z_S_S_S_M1_vrr;
      Double I_ERI_I3x3z_S_Px_S_vrr = QCX*I_ERI_I3x3z_S_S_S_vrr+WQX*I_ERI_I3x3z_S_S_S_M1_vrr+3*oned2k*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_I2x4y_S_Px_S_vrr = QCX*I_ERI_I2x4y_S_S_S_vrr+WQX*I_ERI_I2x4y_S_S_S_M1_vrr+2*oned2k*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Px_S_vrr = QCX*I_ERI_I2x3yz_S_S_S_vrr+WQX*I_ERI_I2x3yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Hx3yz_S_S_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Px_S_vrr = QCX*I_ERI_I2x2y2z_S_S_S_vrr+WQX*I_ERI_I2x2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Hx2y2z_S_S_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Px_S_vrr = QCX*I_ERI_I2xy3z_S_S_S_vrr+WQX*I_ERI_I2xy3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Hxy3z_S_S_S_M1_vrr;
      Double I_ERI_I2x4z_S_Px_S_vrr = QCX*I_ERI_I2x4z_S_S_S_vrr+WQX*I_ERI_I2x4z_S_S_S_M1_vrr+2*oned2k*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5y_S_Px_S_vrr = QCX*I_ERI_Ix5y_S_S_S_vrr+WQX*I_ERI_Ix5y_S_S_S_M1_vrr+oned2k*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Px_S_vrr = QCX*I_ERI_Ix4yz_S_S_S_vrr+WQX*I_ERI_Ix4yz_S_S_S_M1_vrr+oned2k*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Px_S_vrr = QCX*I_ERI_Ix3y2z_S_S_S_vrr+WQX*I_ERI_Ix3y2z_S_S_S_M1_vrr+oned2k*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Px_S_vrr = QCX*I_ERI_Ix2y3z_S_S_S_vrr+WQX*I_ERI_Ix2y3z_S_S_S_M1_vrr+oned2k*I_ERI_H2y3z_S_S_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Px_S_vrr = QCX*I_ERI_Ixy4z_S_S_S_vrr+WQX*I_ERI_Ixy4z_S_S_S_M1_vrr+oned2k*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5z_S_Px_S_vrr = QCX*I_ERI_Ix5z_S_S_S_vrr+WQX*I_ERI_Ix5z_S_S_S_M1_vrr+oned2k*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_I6y_S_Px_S_vrr = QCX*I_ERI_I6y_S_S_S_vrr+WQX*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_I5yz_S_Px_S_vrr = QCX*I_ERI_I5yz_S_S_S_vrr+WQX*I_ERI_I5yz_S_S_S_M1_vrr;
      Double I_ERI_I4y2z_S_Px_S_vrr = QCX*I_ERI_I4y2z_S_S_S_vrr+WQX*I_ERI_I4y2z_S_S_S_M1_vrr;
      Double I_ERI_I3y3z_S_Px_S_vrr = QCX*I_ERI_I3y3z_S_S_S_vrr+WQX*I_ERI_I3y3z_S_S_S_M1_vrr;
      Double I_ERI_I2y4z_S_Px_S_vrr = QCX*I_ERI_I2y4z_S_S_S_vrr+WQX*I_ERI_I2y4z_S_S_S_M1_vrr;
      Double I_ERI_Iy5z_S_Px_S_vrr = QCX*I_ERI_Iy5z_S_S_S_vrr+WQX*I_ERI_Iy5z_S_S_S_M1_vrr;
      Double I_ERI_I6z_S_Px_S_vrr = QCX*I_ERI_I6z_S_S_S_vrr+WQX*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_I6x_S_Py_S_vrr = QCY*I_ERI_I6x_S_S_S_vrr+WQY*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_I5xy_S_Py_S_vrr = QCY*I_ERI_I5xy_S_S_S_vrr+WQY*I_ERI_I5xy_S_S_S_M1_vrr+oned2k*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I5xz_S_Py_S_vrr = QCY*I_ERI_I5xz_S_S_S_vrr+WQY*I_ERI_I5xz_S_S_S_M1_vrr;
      Double I_ERI_I4x2y_S_Py_S_vrr = QCY*I_ERI_I4x2y_S_S_S_vrr+WQY*I_ERI_I4x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_I4xyz_S_Py_S_vrr = QCY*I_ERI_I4xyz_S_S_S_vrr+WQY*I_ERI_I4xyz_S_S_S_M1_vrr+oned2k*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_I4x2z_S_Py_S_vrr = QCY*I_ERI_I4x2z_S_S_S_vrr+WQY*I_ERI_I4x2z_S_S_S_M1_vrr;
      Double I_ERI_I3x3y_S_Py_S_vrr = QCY*I_ERI_I3x3y_S_S_S_vrr+WQY*I_ERI_I3x3y_S_S_S_M1_vrr+3*oned2k*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Py_S_vrr = QCY*I_ERI_I3x2yz_S_S_S_vrr+WQY*I_ERI_I3x2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_H3xyz_S_S_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Py_S_vrr = QCY*I_ERI_I3xy2z_S_S_S_vrr+WQY*I_ERI_I3xy2z_S_S_S_M1_vrr+oned2k*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_I3x3z_S_Py_S_vrr = QCY*I_ERI_I3x3z_S_S_S_vrr+WQY*I_ERI_I3x3z_S_S_S_M1_vrr;
      Double I_ERI_I2x4y_S_Py_S_vrr = QCY*I_ERI_I2x4y_S_S_S_vrr+WQY*I_ERI_I2x4y_S_S_S_M1_vrr+4*oned2k*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Py_S_vrr = QCY*I_ERI_I2x3yz_S_S_S_vrr+WQY*I_ERI_I2x3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_H2x2yz_S_S_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Py_S_vrr = QCY*I_ERI_I2x2y2z_S_S_S_vrr+WQY*I_ERI_I2x2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_H2xy2z_S_S_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Py_S_vrr = QCY*I_ERI_I2xy3z_S_S_S_vrr+WQY*I_ERI_I2xy3z_S_S_S_M1_vrr+oned2k*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_I2x4z_S_Py_S_vrr = QCY*I_ERI_I2x4z_S_S_S_vrr+WQY*I_ERI_I2x4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5y_S_Py_S_vrr = QCY*I_ERI_Ix5y_S_S_S_vrr+WQY*I_ERI_Ix5y_S_S_S_M1_vrr+5*oned2k*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Py_S_vrr = QCY*I_ERI_Ix4yz_S_S_S_vrr+WQY*I_ERI_Ix4yz_S_S_S_M1_vrr+4*oned2k*I_ERI_Hx3yz_S_S_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Py_S_vrr = QCY*I_ERI_Ix3y2z_S_S_S_vrr+WQY*I_ERI_Ix3y2z_S_S_S_M1_vrr+3*oned2k*I_ERI_Hx2y2z_S_S_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Py_S_vrr = QCY*I_ERI_Ix2y3z_S_S_S_vrr+WQY*I_ERI_Ix2y3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Hxy3z_S_S_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Py_S_vrr = QCY*I_ERI_Ixy4z_S_S_S_vrr+WQY*I_ERI_Ixy4z_S_S_S_M1_vrr+oned2k*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_Ix5z_S_Py_S_vrr = QCY*I_ERI_Ix5z_S_S_S_vrr+WQY*I_ERI_Ix5z_S_S_S_M1_vrr;
      Double I_ERI_I6y_S_Py_S_vrr = QCY*I_ERI_I6y_S_S_S_vrr+WQY*I_ERI_I6y_S_S_S_M1_vrr+6*oned2k*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_I5yz_S_Py_S_vrr = QCY*I_ERI_I5yz_S_S_S_vrr+WQY*I_ERI_I5yz_S_S_S_M1_vrr+5*oned2k*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_I4y2z_S_Py_S_vrr = QCY*I_ERI_I4y2z_S_S_S_vrr+WQY*I_ERI_I4y2z_S_S_S_M1_vrr+4*oned2k*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_I3y3z_S_Py_S_vrr = QCY*I_ERI_I3y3z_S_S_S_vrr+WQY*I_ERI_I3y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_H2y3z_S_S_S_M1_vrr;
      Double I_ERI_I2y4z_S_Py_S_vrr = QCY*I_ERI_I2y4z_S_S_S_vrr+WQY*I_ERI_I2y4z_S_S_S_M1_vrr+2*oned2k*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_Iy5z_S_Py_S_vrr = QCY*I_ERI_Iy5z_S_S_S_vrr+WQY*I_ERI_Iy5z_S_S_S_M1_vrr+oned2k*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_I6z_S_Py_S_vrr = QCY*I_ERI_I6z_S_S_S_vrr+WQY*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_I6x_S_Pz_S_vrr = QCZ*I_ERI_I6x_S_S_S_vrr+WQZ*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_I5xy_S_Pz_S_vrr = QCZ*I_ERI_I5xy_S_S_S_vrr+WQZ*I_ERI_I5xy_S_S_S_M1_vrr;
      Double I_ERI_I5xz_S_Pz_S_vrr = QCZ*I_ERI_I5xz_S_S_S_vrr+WQZ*I_ERI_I5xz_S_S_S_M1_vrr+oned2k*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_I4x2y_S_Pz_S_vrr = QCZ*I_ERI_I4x2y_S_S_S_vrr+WQZ*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_I4xyz_S_Pz_S_vrr = QCZ*I_ERI_I4xyz_S_S_S_vrr+WQZ*I_ERI_I4xyz_S_S_S_M1_vrr+oned2k*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_I4x2z_S_Pz_S_vrr = QCZ*I_ERI_I4x2z_S_S_S_vrr+WQZ*I_ERI_I4x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_I3x3y_S_Pz_S_vrr = QCZ*I_ERI_I3x3y_S_S_S_vrr+WQZ*I_ERI_I3x3y_S_S_S_M1_vrr;
      Double I_ERI_I3x2yz_S_Pz_S_vrr = QCZ*I_ERI_I3x2yz_S_S_S_vrr+WQZ*I_ERI_I3x2yz_S_S_S_M1_vrr+oned2k*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_I3xy2z_S_Pz_S_vrr = QCZ*I_ERI_I3xy2z_S_S_S_vrr+WQZ*I_ERI_I3xy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_H3xyz_S_S_S_M1_vrr;
      Double I_ERI_I3x3z_S_Pz_S_vrr = QCZ*I_ERI_I3x3z_S_S_S_vrr+WQZ*I_ERI_I3x3z_S_S_S_M1_vrr+3*oned2k*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_I2x4y_S_Pz_S_vrr = QCZ*I_ERI_I2x4y_S_S_S_vrr+WQZ*I_ERI_I2x4y_S_S_S_M1_vrr;
      Double I_ERI_I2x3yz_S_Pz_S_vrr = QCZ*I_ERI_I2x3yz_S_S_S_vrr+WQZ*I_ERI_I2x3yz_S_S_S_M1_vrr+oned2k*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_I2x2y2z_S_Pz_S_vrr = QCZ*I_ERI_I2x2y2z_S_S_S_vrr+WQZ*I_ERI_I2x2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_H2x2yz_S_S_S_M1_vrr;
      Double I_ERI_I2xy3z_S_Pz_S_vrr = QCZ*I_ERI_I2xy3z_S_S_S_vrr+WQZ*I_ERI_I2xy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_H2xy2z_S_S_S_M1_vrr;
      Double I_ERI_I2x4z_S_Pz_S_vrr = QCZ*I_ERI_I2x4z_S_S_S_vrr+WQZ*I_ERI_I2x4z_S_S_S_M1_vrr+4*oned2k*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_Ix5y_S_Pz_S_vrr = QCZ*I_ERI_Ix5y_S_S_S_vrr+WQZ*I_ERI_Ix5y_S_S_S_M1_vrr;
      Double I_ERI_Ix4yz_S_Pz_S_vrr = QCZ*I_ERI_Ix4yz_S_S_S_vrr+WQZ*I_ERI_Ix4yz_S_S_S_M1_vrr+oned2k*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_Ix3y2z_S_Pz_S_vrr = QCZ*I_ERI_Ix3y2z_S_S_S_vrr+WQZ*I_ERI_Ix3y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Hx3yz_S_S_S_M1_vrr;
      Double I_ERI_Ix2y3z_S_Pz_S_vrr = QCZ*I_ERI_Ix2y3z_S_S_S_vrr+WQZ*I_ERI_Ix2y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Hx2y2z_S_S_S_M1_vrr;
      Double I_ERI_Ixy4z_S_Pz_S_vrr = QCZ*I_ERI_Ixy4z_S_S_S_vrr+WQZ*I_ERI_Ixy4z_S_S_S_M1_vrr+4*oned2k*I_ERI_Hxy3z_S_S_S_M1_vrr;
      Double I_ERI_Ix5z_S_Pz_S_vrr = QCZ*I_ERI_Ix5z_S_S_S_vrr+WQZ*I_ERI_Ix5z_S_S_S_M1_vrr+5*oned2k*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_I6y_S_Pz_S_vrr = QCZ*I_ERI_I6y_S_S_S_vrr+WQZ*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_I5yz_S_Pz_S_vrr = QCZ*I_ERI_I5yz_S_S_S_vrr+WQZ*I_ERI_I5yz_S_S_S_M1_vrr+oned2k*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_I4y2z_S_Pz_S_vrr = QCZ*I_ERI_I4y2z_S_S_S_vrr+WQZ*I_ERI_I4y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_I3y3z_S_Pz_S_vrr = QCZ*I_ERI_I3y3z_S_S_S_vrr+WQZ*I_ERI_I3y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_I2y4z_S_Pz_S_vrr = QCZ*I_ERI_I2y4z_S_S_S_vrr+WQZ*I_ERI_I2y4z_S_S_S_M1_vrr+4*oned2k*I_ERI_H2y3z_S_S_S_M1_vrr;
      Double I_ERI_Iy5z_S_Pz_S_vrr = QCZ*I_ERI_Iy5z_S_S_S_vrr+WQZ*I_ERI_Iy5z_S_S_S_M1_vrr+5*oned2k*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_I6z_S_Pz_S_vrr = QCZ*I_ERI_I6z_S_S_S_vrr+WQZ*I_ERI_I6z_S_S_S_M1_vrr+6*oned2k*I_ERI_H5z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_I_S_S_S
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       ************************************************************/
      Double I_ERI_K7x_S_S_S_vrr = PAX*I_ERI_I6x_S_S_S_vrr+WPX*I_ERI_I6x_S_S_S_M1_vrr+6*oned2z*I_ERI_H5x_S_S_S_vrr-6*rhod2zsq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_K6xy_S_S_S_vrr = PAY*I_ERI_I6x_S_S_S_vrr+WPY*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_K6xz_S_S_S_vrr = PAZ*I_ERI_I6x_S_S_S_vrr+WPZ*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_K5x2y_S_S_S_vrr = PAY*I_ERI_I5xy_S_S_S_vrr+WPY*I_ERI_I5xy_S_S_S_M1_vrr+oned2z*I_ERI_H5x_S_S_S_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_K5xyz_S_S_S_vrr = PAZ*I_ERI_I5xy_S_S_S_vrr+WPZ*I_ERI_I5xy_S_S_S_M1_vrr;
      Double I_ERI_K5x2z_S_S_S_vrr = PAZ*I_ERI_I5xz_S_S_S_vrr+WPZ*I_ERI_I5xz_S_S_S_M1_vrr+oned2z*I_ERI_H5x_S_S_S_vrr-rhod2zsq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_K4x3y_S_S_S_vrr = PAY*I_ERI_I4x2y_S_S_S_vrr+WPY*I_ERI_I4x2y_S_S_S_M1_vrr+2*oned2z*I_ERI_H4xy_S_S_S_vrr-2*rhod2zsq*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_K4x2yz_S_S_S_vrr = PAZ*I_ERI_I4x2y_S_S_S_vrr+WPZ*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_K4xy2z_S_S_S_vrr = PAY*I_ERI_I4x2z_S_S_S_vrr+WPY*I_ERI_I4x2z_S_S_S_M1_vrr;
      Double I_ERI_K4x3z_S_S_S_vrr = PAZ*I_ERI_I4x2z_S_S_S_vrr+WPZ*I_ERI_I4x2z_S_S_S_M1_vrr+2*oned2z*I_ERI_H4xz_S_S_S_vrr-2*rhod2zsq*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_K3x4y_S_S_S_vrr = PAX*I_ERI_I2x4y_S_S_S_vrr+WPX*I_ERI_I2x4y_S_S_S_M1_vrr+2*oned2z*I_ERI_Hx4y_S_S_S_vrr-2*rhod2zsq*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_K3x3yz_S_S_S_vrr = PAZ*I_ERI_I3x3y_S_S_S_vrr+WPZ*I_ERI_I3x3y_S_S_S_M1_vrr;
      Double I_ERI_K3x2y2z_S_S_S_vrr = PAZ*I_ERI_I3x2yz_S_S_S_vrr+WPZ*I_ERI_I3x2yz_S_S_S_M1_vrr+oned2z*I_ERI_H3x2y_S_S_S_vrr-rhod2zsq*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_K3xy3z_S_S_S_vrr = PAY*I_ERI_I3x3z_S_S_S_vrr+WPY*I_ERI_I3x3z_S_S_S_M1_vrr;
      Double I_ERI_K3x4z_S_S_S_vrr = PAX*I_ERI_I2x4z_S_S_S_vrr+WPX*I_ERI_I2x4z_S_S_S_M1_vrr+2*oned2z*I_ERI_Hx4z_S_S_S_vrr-2*rhod2zsq*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5y_S_S_S_vrr = PAX*I_ERI_Ix5y_S_S_S_vrr+WPX*I_ERI_Ix5y_S_S_S_M1_vrr+oned2z*I_ERI_H5y_S_S_S_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_K2x4yz_S_S_S_vrr = PAZ*I_ERI_I2x4y_S_S_S_vrr+WPZ*I_ERI_I2x4y_S_S_S_M1_vrr;
      Double I_ERI_K2x3y2z_S_S_S_vrr = PAZ*I_ERI_I2x3yz_S_S_S_vrr+WPZ*I_ERI_I2x3yz_S_S_S_M1_vrr+oned2z*I_ERI_H2x3y_S_S_S_vrr-rhod2zsq*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_K2x2y3z_S_S_S_vrr = PAY*I_ERI_I2xy3z_S_S_S_vrr+WPY*I_ERI_I2xy3z_S_S_S_M1_vrr+oned2z*I_ERI_H2x3z_S_S_S_vrr-rhod2zsq*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_K2xy4z_S_S_S_vrr = PAY*I_ERI_I2x4z_S_S_S_vrr+WPY*I_ERI_I2x4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5z_S_S_S_vrr = PAX*I_ERI_Ix5z_S_S_S_vrr+WPX*I_ERI_Ix5z_S_S_S_M1_vrr+oned2z*I_ERI_H5z_S_S_S_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6y_S_S_S_vrr = PAX*I_ERI_I6y_S_S_S_vrr+WPX*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_Kx5yz_S_S_S_vrr = PAZ*I_ERI_Ix5y_S_S_S_vrr+WPZ*I_ERI_Ix5y_S_S_S_M1_vrr;
      Double I_ERI_Kx4y2z_S_S_S_vrr = PAX*I_ERI_I4y2z_S_S_S_vrr+WPX*I_ERI_I4y2z_S_S_S_M1_vrr;
      Double I_ERI_Kx3y3z_S_S_S_vrr = PAX*I_ERI_I3y3z_S_S_S_vrr+WPX*I_ERI_I3y3z_S_S_S_M1_vrr;
      Double I_ERI_Kx2y4z_S_S_S_vrr = PAX*I_ERI_I2y4z_S_S_S_vrr+WPX*I_ERI_I2y4z_S_S_S_M1_vrr;
      Double I_ERI_Kxy5z_S_S_S_vrr = PAY*I_ERI_Ix5z_S_S_S_vrr+WPY*I_ERI_Ix5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6z_S_S_S_vrr = PAX*I_ERI_I6z_S_S_S_vrr+WPX*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_K7y_S_S_S_vrr = PAY*I_ERI_I6y_S_S_S_vrr+WPY*I_ERI_I6y_S_S_S_M1_vrr+6*oned2z*I_ERI_H5y_S_S_S_vrr-6*rhod2zsq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_K6yz_S_S_S_vrr = PAZ*I_ERI_I6y_S_S_S_vrr+WPZ*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_K5y2z_S_S_S_vrr = PAZ*I_ERI_I5yz_S_S_S_vrr+WPZ*I_ERI_I5yz_S_S_S_M1_vrr+oned2z*I_ERI_H5y_S_S_S_vrr-rhod2zsq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_K4y3z_S_S_S_vrr = PAZ*I_ERI_I4y2z_S_S_S_vrr+WPZ*I_ERI_I4y2z_S_S_S_M1_vrr+2*oned2z*I_ERI_H4yz_S_S_S_vrr-2*rhod2zsq*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_K3y4z_S_S_S_vrr = PAY*I_ERI_I2y4z_S_S_S_vrr+WPY*I_ERI_I2y4z_S_S_S_M1_vrr+2*oned2z*I_ERI_Hy4z_S_S_S_vrr-2*rhod2zsq*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_K2y5z_S_S_S_vrr = PAY*I_ERI_Iy5z_S_S_S_vrr+WPY*I_ERI_Iy5z_S_S_S_M1_vrr+oned2z*I_ERI_H5z_S_S_S_vrr-rhod2zsq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_Ky6z_S_S_S_vrr = PAY*I_ERI_I6z_S_S_S_vrr+WPY*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_K7z_S_S_S_vrr = PAZ*I_ERI_I6z_S_S_S_vrr+WPZ*I_ERI_I6z_S_S_S_M1_vrr+6*oned2z*I_ERI_H5z_S_S_S_vrr-6*rhod2zsq*I_ERI_H5z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_K_S_S_S
       * RHS shell quartet name: SQ_ERI_K_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       ************************************************************/
      Double I_ERI_K7x_S_Px_S_vrr = QCX*I_ERI_K7x_S_S_S_vrr+WQX*I_ERI_K7x_S_S_S_M1_vrr+7*oned2k*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_K6xy_S_Px_S_vrr = QCX*I_ERI_K6xy_S_S_S_vrr+WQX*I_ERI_K6xy_S_S_S_M1_vrr+6*oned2k*I_ERI_I5xy_S_S_S_M1_vrr;
      Double I_ERI_K6xz_S_Px_S_vrr = QCX*I_ERI_K6xz_S_S_S_vrr+WQX*I_ERI_K6xz_S_S_S_M1_vrr+6*oned2k*I_ERI_I5xz_S_S_S_M1_vrr;
      Double I_ERI_K5x2y_S_Px_S_vrr = QCX*I_ERI_K5x2y_S_S_S_vrr+WQX*I_ERI_K5x2y_S_S_S_M1_vrr+5*oned2k*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_K5xyz_S_Px_S_vrr = QCX*I_ERI_K5xyz_S_S_S_vrr+WQX*I_ERI_K5xyz_S_S_S_M1_vrr+5*oned2k*I_ERI_I4xyz_S_S_S_M1_vrr;
      Double I_ERI_K5x2z_S_Px_S_vrr = QCX*I_ERI_K5x2z_S_S_S_vrr+WQX*I_ERI_K5x2z_S_S_S_M1_vrr+5*oned2k*I_ERI_I4x2z_S_S_S_M1_vrr;
      Double I_ERI_K4x3y_S_Px_S_vrr = QCX*I_ERI_K4x3y_S_S_S_vrr+WQX*I_ERI_K4x3y_S_S_S_M1_vrr+4*oned2k*I_ERI_I3x3y_S_S_S_M1_vrr;
      Double I_ERI_K4x2yz_S_Px_S_vrr = QCX*I_ERI_K4x2yz_S_S_S_vrr+WQX*I_ERI_K4x2yz_S_S_S_M1_vrr+4*oned2k*I_ERI_I3x2yz_S_S_S_M1_vrr;
      Double I_ERI_K4xy2z_S_Px_S_vrr = QCX*I_ERI_K4xy2z_S_S_S_vrr+WQX*I_ERI_K4xy2z_S_S_S_M1_vrr+4*oned2k*I_ERI_I3xy2z_S_S_S_M1_vrr;
      Double I_ERI_K4x3z_S_Px_S_vrr = QCX*I_ERI_K4x3z_S_S_S_vrr+WQX*I_ERI_K4x3z_S_S_S_M1_vrr+4*oned2k*I_ERI_I3x3z_S_S_S_M1_vrr;
      Double I_ERI_K3x4y_S_Px_S_vrr = QCX*I_ERI_K3x4y_S_S_S_vrr+WQX*I_ERI_K3x4y_S_S_S_M1_vrr+3*oned2k*I_ERI_I2x4y_S_S_S_M1_vrr;
      Double I_ERI_K3x3yz_S_Px_S_vrr = QCX*I_ERI_K3x3yz_S_S_S_vrr+WQX*I_ERI_K3x3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_I2x3yz_S_S_S_M1_vrr;
      Double I_ERI_K3x2y2z_S_Px_S_vrr = QCX*I_ERI_K3x2y2z_S_S_S_vrr+WQX*I_ERI_K3x2y2z_S_S_S_M1_vrr+3*oned2k*I_ERI_I2x2y2z_S_S_S_M1_vrr;
      Double I_ERI_K3xy3z_S_Px_S_vrr = QCX*I_ERI_K3xy3z_S_S_S_vrr+WQX*I_ERI_K3xy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_I2xy3z_S_S_S_M1_vrr;
      Double I_ERI_K3x4z_S_Px_S_vrr = QCX*I_ERI_K3x4z_S_S_S_vrr+WQX*I_ERI_K3x4z_S_S_S_M1_vrr+3*oned2k*I_ERI_I2x4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5y_S_Px_S_vrr = QCX*I_ERI_K2x5y_S_S_S_vrr+WQX*I_ERI_K2x5y_S_S_S_M1_vrr+2*oned2k*I_ERI_Ix5y_S_S_S_M1_vrr;
      Double I_ERI_K2x4yz_S_Px_S_vrr = QCX*I_ERI_K2x4yz_S_S_S_vrr+WQX*I_ERI_K2x4yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Ix4yz_S_S_S_M1_vrr;
      Double I_ERI_K2x3y2z_S_Px_S_vrr = QCX*I_ERI_K2x3y2z_S_S_S_vrr+WQX*I_ERI_K2x3y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Ix3y2z_S_S_S_M1_vrr;
      Double I_ERI_K2x2y3z_S_Px_S_vrr = QCX*I_ERI_K2x2y3z_S_S_S_vrr+WQX*I_ERI_K2x2y3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Ix2y3z_S_S_S_M1_vrr;
      Double I_ERI_K2xy4z_S_Px_S_vrr = QCX*I_ERI_K2xy4z_S_S_S_vrr+WQX*I_ERI_K2xy4z_S_S_S_M1_vrr+2*oned2k*I_ERI_Ixy4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5z_S_Px_S_vrr = QCX*I_ERI_K2x5z_S_S_S_vrr+WQX*I_ERI_K2x5z_S_S_S_M1_vrr+2*oned2k*I_ERI_Ix5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6y_S_Px_S_vrr = QCX*I_ERI_Kx6y_S_S_S_vrr+WQX*I_ERI_Kx6y_S_S_S_M1_vrr+oned2k*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_Kx5yz_S_Px_S_vrr = QCX*I_ERI_Kx5yz_S_S_S_vrr+WQX*I_ERI_Kx5yz_S_S_S_M1_vrr+oned2k*I_ERI_I5yz_S_S_S_M1_vrr;
      Double I_ERI_Kx4y2z_S_Px_S_vrr = QCX*I_ERI_Kx4y2z_S_S_S_vrr+WQX*I_ERI_Kx4y2z_S_S_S_M1_vrr+oned2k*I_ERI_I4y2z_S_S_S_M1_vrr;
      Double I_ERI_Kx3y3z_S_Px_S_vrr = QCX*I_ERI_Kx3y3z_S_S_S_vrr+WQX*I_ERI_Kx3y3z_S_S_S_M1_vrr+oned2k*I_ERI_I3y3z_S_S_S_M1_vrr;
      Double I_ERI_Kx2y4z_S_Px_S_vrr = QCX*I_ERI_Kx2y4z_S_S_S_vrr+WQX*I_ERI_Kx2y4z_S_S_S_M1_vrr+oned2k*I_ERI_I2y4z_S_S_S_M1_vrr;
      Double I_ERI_Kxy5z_S_Px_S_vrr = QCX*I_ERI_Kxy5z_S_S_S_vrr+WQX*I_ERI_Kxy5z_S_S_S_M1_vrr+oned2k*I_ERI_Iy5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6z_S_Px_S_vrr = QCX*I_ERI_Kx6z_S_S_S_vrr+WQX*I_ERI_Kx6z_S_S_S_M1_vrr+oned2k*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_K7y_S_Px_S_vrr = QCX*I_ERI_K7y_S_S_S_vrr+WQX*I_ERI_K7y_S_S_S_M1_vrr;
      Double I_ERI_K6yz_S_Px_S_vrr = QCX*I_ERI_K6yz_S_S_S_vrr+WQX*I_ERI_K6yz_S_S_S_M1_vrr;
      Double I_ERI_K5y2z_S_Px_S_vrr = QCX*I_ERI_K5y2z_S_S_S_vrr+WQX*I_ERI_K5y2z_S_S_S_M1_vrr;
      Double I_ERI_K4y3z_S_Px_S_vrr = QCX*I_ERI_K4y3z_S_S_S_vrr+WQX*I_ERI_K4y3z_S_S_S_M1_vrr;
      Double I_ERI_K3y4z_S_Px_S_vrr = QCX*I_ERI_K3y4z_S_S_S_vrr+WQX*I_ERI_K3y4z_S_S_S_M1_vrr;
      Double I_ERI_K2y5z_S_Px_S_vrr = QCX*I_ERI_K2y5z_S_S_S_vrr+WQX*I_ERI_K2y5z_S_S_S_M1_vrr;
      Double I_ERI_Ky6z_S_Px_S_vrr = QCX*I_ERI_Ky6z_S_S_S_vrr+WQX*I_ERI_Ky6z_S_S_S_M1_vrr;
      Double I_ERI_K7z_S_Px_S_vrr = QCX*I_ERI_K7z_S_S_S_vrr+WQX*I_ERI_K7z_S_S_S_M1_vrr;
      Double I_ERI_K7x_S_Py_S_vrr = QCY*I_ERI_K7x_S_S_S_vrr+WQY*I_ERI_K7x_S_S_S_M1_vrr;
      Double I_ERI_K6xy_S_Py_S_vrr = QCY*I_ERI_K6xy_S_S_S_vrr+WQY*I_ERI_K6xy_S_S_S_M1_vrr+oned2k*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_K6xz_S_Py_S_vrr = QCY*I_ERI_K6xz_S_S_S_vrr+WQY*I_ERI_K6xz_S_S_S_M1_vrr;
      Double I_ERI_K5x2y_S_Py_S_vrr = QCY*I_ERI_K5x2y_S_S_S_vrr+WQY*I_ERI_K5x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_I5xy_S_S_S_M1_vrr;
      Double I_ERI_K5xyz_S_Py_S_vrr = QCY*I_ERI_K5xyz_S_S_S_vrr+WQY*I_ERI_K5xyz_S_S_S_M1_vrr+oned2k*I_ERI_I5xz_S_S_S_M1_vrr;
      Double I_ERI_K5x2z_S_Py_S_vrr = QCY*I_ERI_K5x2z_S_S_S_vrr+WQY*I_ERI_K5x2z_S_S_S_M1_vrr;
      Double I_ERI_K4x3y_S_Py_S_vrr = QCY*I_ERI_K4x3y_S_S_S_vrr+WQY*I_ERI_K4x3y_S_S_S_M1_vrr+3*oned2k*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_K4x2yz_S_Py_S_vrr = QCY*I_ERI_K4x2yz_S_S_S_vrr+WQY*I_ERI_K4x2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_I4xyz_S_S_S_M1_vrr;
      Double I_ERI_K4xy2z_S_Py_S_vrr = QCY*I_ERI_K4xy2z_S_S_S_vrr+WQY*I_ERI_K4xy2z_S_S_S_M1_vrr+oned2k*I_ERI_I4x2z_S_S_S_M1_vrr;
      Double I_ERI_K4x3z_S_Py_S_vrr = QCY*I_ERI_K4x3z_S_S_S_vrr+WQY*I_ERI_K4x3z_S_S_S_M1_vrr;
      Double I_ERI_K3x4y_S_Py_S_vrr = QCY*I_ERI_K3x4y_S_S_S_vrr+WQY*I_ERI_K3x4y_S_S_S_M1_vrr+4*oned2k*I_ERI_I3x3y_S_S_S_M1_vrr;
      Double I_ERI_K3x3yz_S_Py_S_vrr = QCY*I_ERI_K3x3yz_S_S_S_vrr+WQY*I_ERI_K3x3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_I3x2yz_S_S_S_M1_vrr;
      Double I_ERI_K3x2y2z_S_Py_S_vrr = QCY*I_ERI_K3x2y2z_S_S_S_vrr+WQY*I_ERI_K3x2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_I3xy2z_S_S_S_M1_vrr;
      Double I_ERI_K3xy3z_S_Py_S_vrr = QCY*I_ERI_K3xy3z_S_S_S_vrr+WQY*I_ERI_K3xy3z_S_S_S_M1_vrr+oned2k*I_ERI_I3x3z_S_S_S_M1_vrr;
      Double I_ERI_K3x4z_S_Py_S_vrr = QCY*I_ERI_K3x4z_S_S_S_vrr+WQY*I_ERI_K3x4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5y_S_Py_S_vrr = QCY*I_ERI_K2x5y_S_S_S_vrr+WQY*I_ERI_K2x5y_S_S_S_M1_vrr+5*oned2k*I_ERI_I2x4y_S_S_S_M1_vrr;
      Double I_ERI_K2x4yz_S_Py_S_vrr = QCY*I_ERI_K2x4yz_S_S_S_vrr+WQY*I_ERI_K2x4yz_S_S_S_M1_vrr+4*oned2k*I_ERI_I2x3yz_S_S_S_M1_vrr;
      Double I_ERI_K2x3y2z_S_Py_S_vrr = QCY*I_ERI_K2x3y2z_S_S_S_vrr+WQY*I_ERI_K2x3y2z_S_S_S_M1_vrr+3*oned2k*I_ERI_I2x2y2z_S_S_S_M1_vrr;
      Double I_ERI_K2x2y3z_S_Py_S_vrr = QCY*I_ERI_K2x2y3z_S_S_S_vrr+WQY*I_ERI_K2x2y3z_S_S_S_M1_vrr+2*oned2k*I_ERI_I2xy3z_S_S_S_M1_vrr;
      Double I_ERI_K2xy4z_S_Py_S_vrr = QCY*I_ERI_K2xy4z_S_S_S_vrr+WQY*I_ERI_K2xy4z_S_S_S_M1_vrr+oned2k*I_ERI_I2x4z_S_S_S_M1_vrr;
      Double I_ERI_K2x5z_S_Py_S_vrr = QCY*I_ERI_K2x5z_S_S_S_vrr+WQY*I_ERI_K2x5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6y_S_Py_S_vrr = QCY*I_ERI_Kx6y_S_S_S_vrr+WQY*I_ERI_Kx6y_S_S_S_M1_vrr+6*oned2k*I_ERI_Ix5y_S_S_S_M1_vrr;
      Double I_ERI_Kx5yz_S_Py_S_vrr = QCY*I_ERI_Kx5yz_S_S_S_vrr+WQY*I_ERI_Kx5yz_S_S_S_M1_vrr+5*oned2k*I_ERI_Ix4yz_S_S_S_M1_vrr;
      Double I_ERI_Kx4y2z_S_Py_S_vrr = QCY*I_ERI_Kx4y2z_S_S_S_vrr+WQY*I_ERI_Kx4y2z_S_S_S_M1_vrr+4*oned2k*I_ERI_Ix3y2z_S_S_S_M1_vrr;
      Double I_ERI_Kx3y3z_S_Py_S_vrr = QCY*I_ERI_Kx3y3z_S_S_S_vrr+WQY*I_ERI_Kx3y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Ix2y3z_S_S_S_M1_vrr;
      Double I_ERI_Kx2y4z_S_Py_S_vrr = QCY*I_ERI_Kx2y4z_S_S_S_vrr+WQY*I_ERI_Kx2y4z_S_S_S_M1_vrr+2*oned2k*I_ERI_Ixy4z_S_S_S_M1_vrr;
      Double I_ERI_Kxy5z_S_Py_S_vrr = QCY*I_ERI_Kxy5z_S_S_S_vrr+WQY*I_ERI_Kxy5z_S_S_S_M1_vrr+oned2k*I_ERI_Ix5z_S_S_S_M1_vrr;
      Double I_ERI_Kx6z_S_Py_S_vrr = QCY*I_ERI_Kx6z_S_S_S_vrr+WQY*I_ERI_Kx6z_S_S_S_M1_vrr;
      Double I_ERI_K7y_S_Py_S_vrr = QCY*I_ERI_K7y_S_S_S_vrr+WQY*I_ERI_K7y_S_S_S_M1_vrr+7*oned2k*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_K6yz_S_Py_S_vrr = QCY*I_ERI_K6yz_S_S_S_vrr+WQY*I_ERI_K6yz_S_S_S_M1_vrr+6*oned2k*I_ERI_I5yz_S_S_S_M1_vrr;
      Double I_ERI_K5y2z_S_Py_S_vrr = QCY*I_ERI_K5y2z_S_S_S_vrr+WQY*I_ERI_K5y2z_S_S_S_M1_vrr+5*oned2k*I_ERI_I4y2z_S_S_S_M1_vrr;
      Double I_ERI_K4y3z_S_Py_S_vrr = QCY*I_ERI_K4y3z_S_S_S_vrr+WQY*I_ERI_K4y3z_S_S_S_M1_vrr+4*oned2k*I_ERI_I3y3z_S_S_S_M1_vrr;
      Double I_ERI_K3y4z_S_Py_S_vrr = QCY*I_ERI_K3y4z_S_S_S_vrr+WQY*I_ERI_K3y4z_S_S_S_M1_vrr+3*oned2k*I_ERI_I2y4z_S_S_S_M1_vrr;
      Double I_ERI_K2y5z_S_Py_S_vrr = QCY*I_ERI_K2y5z_S_S_S_vrr+WQY*I_ERI_K2y5z_S_S_S_M1_vrr+2*oned2k*I_ERI_Iy5z_S_S_S_M1_vrr;
      Double I_ERI_Ky6z_S_Py_S_vrr = QCY*I_ERI_Ky6z_S_S_S_vrr+WQY*I_ERI_Ky6z_S_S_S_M1_vrr+oned2k*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_K7z_S_Py_S_vrr = QCY*I_ERI_K7z_S_S_S_vrr+WQY*I_ERI_K7z_S_S_S_M1_vrr;
      Double I_ERI_K7x_S_Pz_S_vrr = QCZ*I_ERI_K7x_S_S_S_vrr+WQZ*I_ERI_K7x_S_S_S_M1_vrr;
      Double I_ERI_K6xy_S_Pz_S_vrr = QCZ*I_ERI_K6xy_S_S_S_vrr+WQZ*I_ERI_K6xy_S_S_S_M1_vrr;
      Double I_ERI_K6xz_S_Pz_S_vrr = QCZ*I_ERI_K6xz_S_S_S_vrr+WQZ*I_ERI_K6xz_S_S_S_M1_vrr+oned2k*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_K5x2y_S_Pz_S_vrr = QCZ*I_ERI_K5x2y_S_S_S_vrr+WQZ*I_ERI_K5x2y_S_S_S_M1_vrr;
      Double I_ERI_K5xyz_S_Pz_S_vrr = QCZ*I_ERI_K5xyz_S_S_S_vrr+WQZ*I_ERI_K5xyz_S_S_S_M1_vrr+oned2k*I_ERI_I5xy_S_S_S_M1_vrr;
      Double I_ERI_K5x2z_S_Pz_S_vrr = QCZ*I_ERI_K5x2z_S_S_S_vrr+WQZ*I_ERI_K5x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_I5xz_S_S_S_M1_vrr;
      Double I_ERI_K4x3y_S_Pz_S_vrr = QCZ*I_ERI_K4x3y_S_S_S_vrr+WQZ*I_ERI_K4x3y_S_S_S_M1_vrr;
      Double I_ERI_K4x2yz_S_Pz_S_vrr = QCZ*I_ERI_K4x2yz_S_S_S_vrr+WQZ*I_ERI_K4x2yz_S_S_S_M1_vrr+oned2k*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_K4xy2z_S_Pz_S_vrr = QCZ*I_ERI_K4xy2z_S_S_S_vrr+WQZ*I_ERI_K4xy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_I4xyz_S_S_S_M1_vrr;
      Double I_ERI_K4x3z_S_Pz_S_vrr = QCZ*I_ERI_K4x3z_S_S_S_vrr+WQZ*I_ERI_K4x3z_S_S_S_M1_vrr+3*oned2k*I_ERI_I4x2z_S_S_S_M1_vrr;
      Double I_ERI_K3x4y_S_Pz_S_vrr = QCZ*I_ERI_K3x4y_S_S_S_vrr+WQZ*I_ERI_K3x4y_S_S_S_M1_vrr;
      Double I_ERI_K3x3yz_S_Pz_S_vrr = QCZ*I_ERI_K3x3yz_S_S_S_vrr+WQZ*I_ERI_K3x3yz_S_S_S_M1_vrr+oned2k*I_ERI_I3x3y_S_S_S_M1_vrr;
      Double I_ERI_K3x2y2z_S_Pz_S_vrr = QCZ*I_ERI_K3x2y2z_S_S_S_vrr+WQZ*I_ERI_K3x2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_I3x2yz_S_S_S_M1_vrr;
      Double I_ERI_K3xy3z_S_Pz_S_vrr = QCZ*I_ERI_K3xy3z_S_S_S_vrr+WQZ*I_ERI_K3xy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_I3xy2z_S_S_S_M1_vrr;
      Double I_ERI_K3x4z_S_Pz_S_vrr = QCZ*I_ERI_K3x4z_S_S_S_vrr+WQZ*I_ERI_K3x4z_S_S_S_M1_vrr+4*oned2k*I_ERI_I3x3z_S_S_S_M1_vrr;
      Double I_ERI_K2x5y_S_Pz_S_vrr = QCZ*I_ERI_K2x5y_S_S_S_vrr+WQZ*I_ERI_K2x5y_S_S_S_M1_vrr;
      Double I_ERI_K2x4yz_S_Pz_S_vrr = QCZ*I_ERI_K2x4yz_S_S_S_vrr+WQZ*I_ERI_K2x4yz_S_S_S_M1_vrr+oned2k*I_ERI_I2x4y_S_S_S_M1_vrr;
      Double I_ERI_K2x3y2z_S_Pz_S_vrr = QCZ*I_ERI_K2x3y2z_S_S_S_vrr+WQZ*I_ERI_K2x3y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_I2x3yz_S_S_S_M1_vrr;
      Double I_ERI_K2x2y3z_S_Pz_S_vrr = QCZ*I_ERI_K2x2y3z_S_S_S_vrr+WQZ*I_ERI_K2x2y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_I2x2y2z_S_S_S_M1_vrr;
      Double I_ERI_K2xy4z_S_Pz_S_vrr = QCZ*I_ERI_K2xy4z_S_S_S_vrr+WQZ*I_ERI_K2xy4z_S_S_S_M1_vrr+4*oned2k*I_ERI_I2xy3z_S_S_S_M1_vrr;
      Double I_ERI_K2x5z_S_Pz_S_vrr = QCZ*I_ERI_K2x5z_S_S_S_vrr+WQZ*I_ERI_K2x5z_S_S_S_M1_vrr+5*oned2k*I_ERI_I2x4z_S_S_S_M1_vrr;
      Double I_ERI_Kx6y_S_Pz_S_vrr = QCZ*I_ERI_Kx6y_S_S_S_vrr+WQZ*I_ERI_Kx6y_S_S_S_M1_vrr;
      Double I_ERI_Kx5yz_S_Pz_S_vrr = QCZ*I_ERI_Kx5yz_S_S_S_vrr+WQZ*I_ERI_Kx5yz_S_S_S_M1_vrr+oned2k*I_ERI_Ix5y_S_S_S_M1_vrr;
      Double I_ERI_Kx4y2z_S_Pz_S_vrr = QCZ*I_ERI_Kx4y2z_S_S_S_vrr+WQZ*I_ERI_Kx4y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Ix4yz_S_S_S_M1_vrr;
      Double I_ERI_Kx3y3z_S_Pz_S_vrr = QCZ*I_ERI_Kx3y3z_S_S_S_vrr+WQZ*I_ERI_Kx3y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Ix3y2z_S_S_S_M1_vrr;
      Double I_ERI_Kx2y4z_S_Pz_S_vrr = QCZ*I_ERI_Kx2y4z_S_S_S_vrr+WQZ*I_ERI_Kx2y4z_S_S_S_M1_vrr+4*oned2k*I_ERI_Ix2y3z_S_S_S_M1_vrr;
      Double I_ERI_Kxy5z_S_Pz_S_vrr = QCZ*I_ERI_Kxy5z_S_S_S_vrr+WQZ*I_ERI_Kxy5z_S_S_S_M1_vrr+5*oned2k*I_ERI_Ixy4z_S_S_S_M1_vrr;
      Double I_ERI_Kx6z_S_Pz_S_vrr = QCZ*I_ERI_Kx6z_S_S_S_vrr+WQZ*I_ERI_Kx6z_S_S_S_M1_vrr+6*oned2k*I_ERI_Ix5z_S_S_S_M1_vrr;
      Double I_ERI_K7y_S_Pz_S_vrr = QCZ*I_ERI_K7y_S_S_S_vrr+WQZ*I_ERI_K7y_S_S_S_M1_vrr;
      Double I_ERI_K6yz_S_Pz_S_vrr = QCZ*I_ERI_K6yz_S_S_S_vrr+WQZ*I_ERI_K6yz_S_S_S_M1_vrr+oned2k*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_K5y2z_S_Pz_S_vrr = QCZ*I_ERI_K5y2z_S_S_S_vrr+WQZ*I_ERI_K5y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_I5yz_S_S_S_M1_vrr;
      Double I_ERI_K4y3z_S_Pz_S_vrr = QCZ*I_ERI_K4y3z_S_S_S_vrr+WQZ*I_ERI_K4y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_I4y2z_S_S_S_M1_vrr;
      Double I_ERI_K3y4z_S_Pz_S_vrr = QCZ*I_ERI_K3y4z_S_S_S_vrr+WQZ*I_ERI_K3y4z_S_S_S_M1_vrr+4*oned2k*I_ERI_I3y3z_S_S_S_M1_vrr;
      Double I_ERI_K2y5z_S_Pz_S_vrr = QCZ*I_ERI_K2y5z_S_S_S_vrr+WQZ*I_ERI_K2y5z_S_S_S_M1_vrr+5*oned2k*I_ERI_I2y4z_S_S_S_M1_vrr;
      Double I_ERI_Ky6z_S_Pz_S_vrr = QCZ*I_ERI_Ky6z_S_S_S_vrr+WQZ*I_ERI_Ky6z_S_S_S_M1_vrr+6*oned2k*I_ERI_Iy5z_S_S_S_M1_vrr;
      Double I_ERI_K7z_S_Pz_S_vrr = QCZ*I_ERI_K7z_S_S_S_vrr+WQZ*I_ERI_K7z_S_S_S_M1_vrr+7*oned2k*I_ERI_I6z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_L_S_S_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_K_S_S_S
       * RHS shell quartet name: SQ_ERI_K_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_I_S_S_S
       * RHS shell quartet name: SQ_ERI_I_S_S_S_M1
       ************************************************************/
      Double I_ERI_L8x_S_S_S_vrr = PAX*I_ERI_K7x_S_S_S_vrr+WPX*I_ERI_K7x_S_S_S_M1_vrr+7*oned2z*I_ERI_I6x_S_S_S_vrr-7*rhod2zsq*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_L7xy_S_S_S_vrr = PAY*I_ERI_K7x_S_S_S_vrr+WPY*I_ERI_K7x_S_S_S_M1_vrr;
      Double I_ERI_L7xz_S_S_S_vrr = PAZ*I_ERI_K7x_S_S_S_vrr+WPZ*I_ERI_K7x_S_S_S_M1_vrr;
      Double I_ERI_L6x2y_S_S_S_vrr = PAY*I_ERI_K6xy_S_S_S_vrr+WPY*I_ERI_K6xy_S_S_S_M1_vrr+oned2z*I_ERI_I6x_S_S_S_vrr-rhod2zsq*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_L6xyz_S_S_S_vrr = PAZ*I_ERI_K6xy_S_S_S_vrr+WPZ*I_ERI_K6xy_S_S_S_M1_vrr;
      Double I_ERI_L6x2z_S_S_S_vrr = PAZ*I_ERI_K6xz_S_S_S_vrr+WPZ*I_ERI_K6xz_S_S_S_M1_vrr+oned2z*I_ERI_I6x_S_S_S_vrr-rhod2zsq*I_ERI_I6x_S_S_S_M1_vrr;
      Double I_ERI_L5x3y_S_S_S_vrr = PAY*I_ERI_K5x2y_S_S_S_vrr+WPY*I_ERI_K5x2y_S_S_S_M1_vrr+2*oned2z*I_ERI_I5xy_S_S_S_vrr-2*rhod2zsq*I_ERI_I5xy_S_S_S_M1_vrr;
      Double I_ERI_L5x2yz_S_S_S_vrr = PAZ*I_ERI_K5x2y_S_S_S_vrr+WPZ*I_ERI_K5x2y_S_S_S_M1_vrr;
      Double I_ERI_L5xy2z_S_S_S_vrr = PAY*I_ERI_K5x2z_S_S_S_vrr+WPY*I_ERI_K5x2z_S_S_S_M1_vrr;
      Double I_ERI_L5x3z_S_S_S_vrr = PAZ*I_ERI_K5x2z_S_S_S_vrr+WPZ*I_ERI_K5x2z_S_S_S_M1_vrr+2*oned2z*I_ERI_I5xz_S_S_S_vrr-2*rhod2zsq*I_ERI_I5xz_S_S_S_M1_vrr;
      Double I_ERI_L4x4y_S_S_S_vrr = PAY*I_ERI_K4x3y_S_S_S_vrr+WPY*I_ERI_K4x3y_S_S_S_M1_vrr+3*oned2z*I_ERI_I4x2y_S_S_S_vrr-3*rhod2zsq*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_L4x3yz_S_S_S_vrr = PAZ*I_ERI_K4x3y_S_S_S_vrr+WPZ*I_ERI_K4x3y_S_S_S_M1_vrr;
      Double I_ERI_L4x2y2z_S_S_S_vrr = PAZ*I_ERI_K4x2yz_S_S_S_vrr+WPZ*I_ERI_K4x2yz_S_S_S_M1_vrr+oned2z*I_ERI_I4x2y_S_S_S_vrr-rhod2zsq*I_ERI_I4x2y_S_S_S_M1_vrr;
      Double I_ERI_L4xy3z_S_S_S_vrr = PAY*I_ERI_K4x3z_S_S_S_vrr+WPY*I_ERI_K4x3z_S_S_S_M1_vrr;
      Double I_ERI_L4x4z_S_S_S_vrr = PAZ*I_ERI_K4x3z_S_S_S_vrr+WPZ*I_ERI_K4x3z_S_S_S_M1_vrr+3*oned2z*I_ERI_I4x2z_S_S_S_vrr-3*rhod2zsq*I_ERI_I4x2z_S_S_S_M1_vrr;
      Double I_ERI_L3x5y_S_S_S_vrr = PAX*I_ERI_K2x5y_S_S_S_vrr+WPX*I_ERI_K2x5y_S_S_S_M1_vrr+2*oned2z*I_ERI_Ix5y_S_S_S_vrr-2*rhod2zsq*I_ERI_Ix5y_S_S_S_M1_vrr;
      Double I_ERI_L3x4yz_S_S_S_vrr = PAZ*I_ERI_K3x4y_S_S_S_vrr+WPZ*I_ERI_K3x4y_S_S_S_M1_vrr;
      Double I_ERI_L3x3y2z_S_S_S_vrr = PAZ*I_ERI_K3x3yz_S_S_S_vrr+WPZ*I_ERI_K3x3yz_S_S_S_M1_vrr+oned2z*I_ERI_I3x3y_S_S_S_vrr-rhod2zsq*I_ERI_I3x3y_S_S_S_M1_vrr;
      Double I_ERI_L3x2y3z_S_S_S_vrr = PAY*I_ERI_K3xy3z_S_S_S_vrr+WPY*I_ERI_K3xy3z_S_S_S_M1_vrr+oned2z*I_ERI_I3x3z_S_S_S_vrr-rhod2zsq*I_ERI_I3x3z_S_S_S_M1_vrr;
      Double I_ERI_L3xy4z_S_S_S_vrr = PAY*I_ERI_K3x4z_S_S_S_vrr+WPY*I_ERI_K3x4z_S_S_S_M1_vrr;
      Double I_ERI_L3x5z_S_S_S_vrr = PAX*I_ERI_K2x5z_S_S_S_vrr+WPX*I_ERI_K2x5z_S_S_S_M1_vrr+2*oned2z*I_ERI_Ix5z_S_S_S_vrr-2*rhod2zsq*I_ERI_Ix5z_S_S_S_M1_vrr;
      Double I_ERI_L2x6y_S_S_S_vrr = PAX*I_ERI_Kx6y_S_S_S_vrr+WPX*I_ERI_Kx6y_S_S_S_M1_vrr+oned2z*I_ERI_I6y_S_S_S_vrr-rhod2zsq*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_L2x5yz_S_S_S_vrr = PAZ*I_ERI_K2x5y_S_S_S_vrr+WPZ*I_ERI_K2x5y_S_S_S_M1_vrr;
      Double I_ERI_L2x4y2z_S_S_S_vrr = PAZ*I_ERI_K2x4yz_S_S_S_vrr+WPZ*I_ERI_K2x4yz_S_S_S_M1_vrr+oned2z*I_ERI_I2x4y_S_S_S_vrr-rhod2zsq*I_ERI_I2x4y_S_S_S_M1_vrr;
      Double I_ERI_L2x3y3z_S_S_S_vrr = PAX*I_ERI_Kx3y3z_S_S_S_vrr+WPX*I_ERI_Kx3y3z_S_S_S_M1_vrr+oned2z*I_ERI_I3y3z_S_S_S_vrr-rhod2zsq*I_ERI_I3y3z_S_S_S_M1_vrr;
      Double I_ERI_L2x2y4z_S_S_S_vrr = PAY*I_ERI_K2xy4z_S_S_S_vrr+WPY*I_ERI_K2xy4z_S_S_S_M1_vrr+oned2z*I_ERI_I2x4z_S_S_S_vrr-rhod2zsq*I_ERI_I2x4z_S_S_S_M1_vrr;
      Double I_ERI_L2xy5z_S_S_S_vrr = PAY*I_ERI_K2x5z_S_S_S_vrr+WPY*I_ERI_K2x5z_S_S_S_M1_vrr;
      Double I_ERI_L2x6z_S_S_S_vrr = PAX*I_ERI_Kx6z_S_S_S_vrr+WPX*I_ERI_Kx6z_S_S_S_M1_vrr+oned2z*I_ERI_I6z_S_S_S_vrr-rhod2zsq*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_Lx7y_S_S_S_vrr = PAX*I_ERI_K7y_S_S_S_vrr+WPX*I_ERI_K7y_S_S_S_M1_vrr;
      Double I_ERI_Lx6yz_S_S_S_vrr = PAZ*I_ERI_Kx6y_S_S_S_vrr+WPZ*I_ERI_Kx6y_S_S_S_M1_vrr;
      Double I_ERI_Lx5y2z_S_S_S_vrr = PAX*I_ERI_K5y2z_S_S_S_vrr+WPX*I_ERI_K5y2z_S_S_S_M1_vrr;
      Double I_ERI_Lx4y3z_S_S_S_vrr = PAX*I_ERI_K4y3z_S_S_S_vrr+WPX*I_ERI_K4y3z_S_S_S_M1_vrr;
      Double I_ERI_Lx3y4z_S_S_S_vrr = PAX*I_ERI_K3y4z_S_S_S_vrr+WPX*I_ERI_K3y4z_S_S_S_M1_vrr;
      Double I_ERI_Lx2y5z_S_S_S_vrr = PAX*I_ERI_K2y5z_S_S_S_vrr+WPX*I_ERI_K2y5z_S_S_S_M1_vrr;
      Double I_ERI_Lxy6z_S_S_S_vrr = PAY*I_ERI_Kx6z_S_S_S_vrr+WPY*I_ERI_Kx6z_S_S_S_M1_vrr;
      Double I_ERI_Lx7z_S_S_S_vrr = PAX*I_ERI_K7z_S_S_S_vrr+WPX*I_ERI_K7z_S_S_S_M1_vrr;
      Double I_ERI_L8y_S_S_S_vrr = PAY*I_ERI_K7y_S_S_S_vrr+WPY*I_ERI_K7y_S_S_S_M1_vrr+7*oned2z*I_ERI_I6y_S_S_S_vrr-7*rhod2zsq*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_L7yz_S_S_S_vrr = PAZ*I_ERI_K7y_S_S_S_vrr+WPZ*I_ERI_K7y_S_S_S_M1_vrr;
      Double I_ERI_L6y2z_S_S_S_vrr = PAZ*I_ERI_K6yz_S_S_S_vrr+WPZ*I_ERI_K6yz_S_S_S_M1_vrr+oned2z*I_ERI_I6y_S_S_S_vrr-rhod2zsq*I_ERI_I6y_S_S_S_M1_vrr;
      Double I_ERI_L5y3z_S_S_S_vrr = PAZ*I_ERI_K5y2z_S_S_S_vrr+WPZ*I_ERI_K5y2z_S_S_S_M1_vrr+2*oned2z*I_ERI_I5yz_S_S_S_vrr-2*rhod2zsq*I_ERI_I5yz_S_S_S_M1_vrr;
      Double I_ERI_L4y4z_S_S_S_vrr = PAZ*I_ERI_K4y3z_S_S_S_vrr+WPZ*I_ERI_K4y3z_S_S_S_M1_vrr+3*oned2z*I_ERI_I4y2z_S_S_S_vrr-3*rhod2zsq*I_ERI_I4y2z_S_S_S_M1_vrr;
      Double I_ERI_L3y5z_S_S_S_vrr = PAY*I_ERI_K2y5z_S_S_S_vrr+WPY*I_ERI_K2y5z_S_S_S_M1_vrr+2*oned2z*I_ERI_Iy5z_S_S_S_vrr-2*rhod2zsq*I_ERI_Iy5z_S_S_S_M1_vrr;
      Double I_ERI_L2y6z_S_S_S_vrr = PAY*I_ERI_Ky6z_S_S_S_vrr+WPY*I_ERI_Ky6z_S_S_S_M1_vrr+oned2z*I_ERI_I6z_S_S_S_vrr-rhod2zsq*I_ERI_I6z_S_S_S_M1_vrr;
      Double I_ERI_Ly7z_S_S_S_vrr = PAY*I_ERI_K7z_S_S_S_vrr+WPY*I_ERI_K7z_S_S_S_M1_vrr;
      Double I_ERI_L8z_S_S_S_vrr = PAZ*I_ERI_K7z_S_S_S_vrr+WPZ*I_ERI_K7z_S_S_S_M1_vrr+7*oned2z*I_ERI_I6z_S_S_S_vrr-7*rhod2zsq*I_ERI_I6z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_L_S_P_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_L_S_S_S
       * RHS shell quartet name: SQ_ERI_L_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_K_S_S_S_M1
       ************************************************************/
      Double I_ERI_L8x_S_Px_S_vrr = QCX*I_ERI_L8x_S_S_S_vrr+WQX*I_ERI_L8x_S_S_S_M1_vrr+8*oned2k*I_ERI_K7x_S_S_S_M1_vrr;
      Double I_ERI_L7xy_S_Px_S_vrr = QCX*I_ERI_L7xy_S_S_S_vrr+WQX*I_ERI_L7xy_S_S_S_M1_vrr+7*oned2k*I_ERI_K6xy_S_S_S_M1_vrr;
      Double I_ERI_L7xz_S_Px_S_vrr = QCX*I_ERI_L7xz_S_S_S_vrr+WQX*I_ERI_L7xz_S_S_S_M1_vrr+7*oned2k*I_ERI_K6xz_S_S_S_M1_vrr;
      Double I_ERI_L6x2y_S_Px_S_vrr = QCX*I_ERI_L6x2y_S_S_S_vrr+WQX*I_ERI_L6x2y_S_S_S_M1_vrr+6*oned2k*I_ERI_K5x2y_S_S_S_M1_vrr;
      Double I_ERI_L6xyz_S_Px_S_vrr = QCX*I_ERI_L6xyz_S_S_S_vrr+WQX*I_ERI_L6xyz_S_S_S_M1_vrr+6*oned2k*I_ERI_K5xyz_S_S_S_M1_vrr;
      Double I_ERI_L6x2z_S_Px_S_vrr = QCX*I_ERI_L6x2z_S_S_S_vrr+WQX*I_ERI_L6x2z_S_S_S_M1_vrr+6*oned2k*I_ERI_K5x2z_S_S_S_M1_vrr;
      Double I_ERI_L5x3y_S_Px_S_vrr = QCX*I_ERI_L5x3y_S_S_S_vrr+WQX*I_ERI_L5x3y_S_S_S_M1_vrr+5*oned2k*I_ERI_K4x3y_S_S_S_M1_vrr;
      Double I_ERI_L5x2yz_S_Px_S_vrr = QCX*I_ERI_L5x2yz_S_S_S_vrr+WQX*I_ERI_L5x2yz_S_S_S_M1_vrr+5*oned2k*I_ERI_K4x2yz_S_S_S_M1_vrr;
      Double I_ERI_L5xy2z_S_Px_S_vrr = QCX*I_ERI_L5xy2z_S_S_S_vrr+WQX*I_ERI_L5xy2z_S_S_S_M1_vrr+5*oned2k*I_ERI_K4xy2z_S_S_S_M1_vrr;
      Double I_ERI_L5x3z_S_Px_S_vrr = QCX*I_ERI_L5x3z_S_S_S_vrr+WQX*I_ERI_L5x3z_S_S_S_M1_vrr+5*oned2k*I_ERI_K4x3z_S_S_S_M1_vrr;
      Double I_ERI_L4x4y_S_Px_S_vrr = QCX*I_ERI_L4x4y_S_S_S_vrr+WQX*I_ERI_L4x4y_S_S_S_M1_vrr+4*oned2k*I_ERI_K3x4y_S_S_S_M1_vrr;
      Double I_ERI_L4x3yz_S_Px_S_vrr = QCX*I_ERI_L4x3yz_S_S_S_vrr+WQX*I_ERI_L4x3yz_S_S_S_M1_vrr+4*oned2k*I_ERI_K3x3yz_S_S_S_M1_vrr;
      Double I_ERI_L4x2y2z_S_Px_S_vrr = QCX*I_ERI_L4x2y2z_S_S_S_vrr+WQX*I_ERI_L4x2y2z_S_S_S_M1_vrr+4*oned2k*I_ERI_K3x2y2z_S_S_S_M1_vrr;
      Double I_ERI_L4xy3z_S_Px_S_vrr = QCX*I_ERI_L4xy3z_S_S_S_vrr+WQX*I_ERI_L4xy3z_S_S_S_M1_vrr+4*oned2k*I_ERI_K3xy3z_S_S_S_M1_vrr;
      Double I_ERI_L4x4z_S_Px_S_vrr = QCX*I_ERI_L4x4z_S_S_S_vrr+WQX*I_ERI_L4x4z_S_S_S_M1_vrr+4*oned2k*I_ERI_K3x4z_S_S_S_M1_vrr;
      Double I_ERI_L3x5y_S_Px_S_vrr = QCX*I_ERI_L3x5y_S_S_S_vrr+WQX*I_ERI_L3x5y_S_S_S_M1_vrr+3*oned2k*I_ERI_K2x5y_S_S_S_M1_vrr;
      Double I_ERI_L3x4yz_S_Px_S_vrr = QCX*I_ERI_L3x4yz_S_S_S_vrr+WQX*I_ERI_L3x4yz_S_S_S_M1_vrr+3*oned2k*I_ERI_K2x4yz_S_S_S_M1_vrr;
      Double I_ERI_L3x3y2z_S_Px_S_vrr = QCX*I_ERI_L3x3y2z_S_S_S_vrr+WQX*I_ERI_L3x3y2z_S_S_S_M1_vrr+3*oned2k*I_ERI_K2x3y2z_S_S_S_M1_vrr;
      Double I_ERI_L3x2y3z_S_Px_S_vrr = QCX*I_ERI_L3x2y3z_S_S_S_vrr+WQX*I_ERI_L3x2y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_K2x2y3z_S_S_S_M1_vrr;
      Double I_ERI_L3xy4z_S_Px_S_vrr = QCX*I_ERI_L3xy4z_S_S_S_vrr+WQX*I_ERI_L3xy4z_S_S_S_M1_vrr+3*oned2k*I_ERI_K2xy4z_S_S_S_M1_vrr;
      Double I_ERI_L3x5z_S_Px_S_vrr = QCX*I_ERI_L3x5z_S_S_S_vrr+WQX*I_ERI_L3x5z_S_S_S_M1_vrr+3*oned2k*I_ERI_K2x5z_S_S_S_M1_vrr;
      Double I_ERI_L2x6y_S_Px_S_vrr = QCX*I_ERI_L2x6y_S_S_S_vrr+WQX*I_ERI_L2x6y_S_S_S_M1_vrr+2*oned2k*I_ERI_Kx6y_S_S_S_M1_vrr;
      Double I_ERI_L2x5yz_S_Px_S_vrr = QCX*I_ERI_L2x5yz_S_S_S_vrr+WQX*I_ERI_L2x5yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Kx5yz_S_S_S_M1_vrr;
      Double I_ERI_L2x4y2z_S_Px_S_vrr = QCX*I_ERI_L2x4y2z_S_S_S_vrr+WQX*I_ERI_L2x4y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Kx4y2z_S_S_S_M1_vrr;
      Double I_ERI_L2x3y3z_S_Px_S_vrr = QCX*I_ERI_L2x3y3z_S_S_S_vrr+WQX*I_ERI_L2x3y3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Kx3y3z_S_S_S_M1_vrr;
      Double I_ERI_L2x2y4z_S_Px_S_vrr = QCX*I_ERI_L2x2y4z_S_S_S_vrr+WQX*I_ERI_L2x2y4z_S_S_S_M1_vrr+2*oned2k*I_ERI_Kx2y4z_S_S_S_M1_vrr;
      Double I_ERI_L2xy5z_S_Px_S_vrr = QCX*I_ERI_L2xy5z_S_S_S_vrr+WQX*I_ERI_L2xy5z_S_S_S_M1_vrr+2*oned2k*I_ERI_Kxy5z_S_S_S_M1_vrr;
      Double I_ERI_L2x6z_S_Px_S_vrr = QCX*I_ERI_L2x6z_S_S_S_vrr+WQX*I_ERI_L2x6z_S_S_S_M1_vrr+2*oned2k*I_ERI_Kx6z_S_S_S_M1_vrr;
      Double I_ERI_Lx7y_S_Px_S_vrr = QCX*I_ERI_Lx7y_S_S_S_vrr+WQX*I_ERI_Lx7y_S_S_S_M1_vrr+oned2k*I_ERI_K7y_S_S_S_M1_vrr;
      Double I_ERI_Lx6yz_S_Px_S_vrr = QCX*I_ERI_Lx6yz_S_S_S_vrr+WQX*I_ERI_Lx6yz_S_S_S_M1_vrr+oned2k*I_ERI_K6yz_S_S_S_M1_vrr;
      Double I_ERI_Lx5y2z_S_Px_S_vrr = QCX*I_ERI_Lx5y2z_S_S_S_vrr+WQX*I_ERI_Lx5y2z_S_S_S_M1_vrr+oned2k*I_ERI_K5y2z_S_S_S_M1_vrr;
      Double I_ERI_Lx4y3z_S_Px_S_vrr = QCX*I_ERI_Lx4y3z_S_S_S_vrr+WQX*I_ERI_Lx4y3z_S_S_S_M1_vrr+oned2k*I_ERI_K4y3z_S_S_S_M1_vrr;
      Double I_ERI_Lx3y4z_S_Px_S_vrr = QCX*I_ERI_Lx3y4z_S_S_S_vrr+WQX*I_ERI_Lx3y4z_S_S_S_M1_vrr+oned2k*I_ERI_K3y4z_S_S_S_M1_vrr;
      Double I_ERI_Lx2y5z_S_Px_S_vrr = QCX*I_ERI_Lx2y5z_S_S_S_vrr+WQX*I_ERI_Lx2y5z_S_S_S_M1_vrr+oned2k*I_ERI_K2y5z_S_S_S_M1_vrr;
      Double I_ERI_Lxy6z_S_Px_S_vrr = QCX*I_ERI_Lxy6z_S_S_S_vrr+WQX*I_ERI_Lxy6z_S_S_S_M1_vrr+oned2k*I_ERI_Ky6z_S_S_S_M1_vrr;
      Double I_ERI_Lx7z_S_Px_S_vrr = QCX*I_ERI_Lx7z_S_S_S_vrr+WQX*I_ERI_Lx7z_S_S_S_M1_vrr+oned2k*I_ERI_K7z_S_S_S_M1_vrr;
      Double I_ERI_L8y_S_Px_S_vrr = QCX*I_ERI_L8y_S_S_S_vrr+WQX*I_ERI_L8y_S_S_S_M1_vrr;
      Double I_ERI_L7yz_S_Px_S_vrr = QCX*I_ERI_L7yz_S_S_S_vrr+WQX*I_ERI_L7yz_S_S_S_M1_vrr;
      Double I_ERI_L6y2z_S_Px_S_vrr = QCX*I_ERI_L6y2z_S_S_S_vrr+WQX*I_ERI_L6y2z_S_S_S_M1_vrr;
      Double I_ERI_L5y3z_S_Px_S_vrr = QCX*I_ERI_L5y3z_S_S_S_vrr+WQX*I_ERI_L5y3z_S_S_S_M1_vrr;
      Double I_ERI_L4y4z_S_Px_S_vrr = QCX*I_ERI_L4y4z_S_S_S_vrr+WQX*I_ERI_L4y4z_S_S_S_M1_vrr;
      Double I_ERI_L3y5z_S_Px_S_vrr = QCX*I_ERI_L3y5z_S_S_S_vrr+WQX*I_ERI_L3y5z_S_S_S_M1_vrr;
      Double I_ERI_L2y6z_S_Px_S_vrr = QCX*I_ERI_L2y6z_S_S_S_vrr+WQX*I_ERI_L2y6z_S_S_S_M1_vrr;
      Double I_ERI_Ly7z_S_Px_S_vrr = QCX*I_ERI_Ly7z_S_S_S_vrr+WQX*I_ERI_Ly7z_S_S_S_M1_vrr;
      Double I_ERI_L8z_S_Px_S_vrr = QCX*I_ERI_L8z_S_S_S_vrr+WQX*I_ERI_L8z_S_S_S_M1_vrr;
      Double I_ERI_L8x_S_Py_S_vrr = QCY*I_ERI_L8x_S_S_S_vrr+WQY*I_ERI_L8x_S_S_S_M1_vrr;
      Double I_ERI_L7xy_S_Py_S_vrr = QCY*I_ERI_L7xy_S_S_S_vrr+WQY*I_ERI_L7xy_S_S_S_M1_vrr+oned2k*I_ERI_K7x_S_S_S_M1_vrr;
      Double I_ERI_L7xz_S_Py_S_vrr = QCY*I_ERI_L7xz_S_S_S_vrr+WQY*I_ERI_L7xz_S_S_S_M1_vrr;
      Double I_ERI_L6x2y_S_Py_S_vrr = QCY*I_ERI_L6x2y_S_S_S_vrr+WQY*I_ERI_L6x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_K6xy_S_S_S_M1_vrr;
      Double I_ERI_L6xyz_S_Py_S_vrr = QCY*I_ERI_L6xyz_S_S_S_vrr+WQY*I_ERI_L6xyz_S_S_S_M1_vrr+oned2k*I_ERI_K6xz_S_S_S_M1_vrr;
      Double I_ERI_L6x2z_S_Py_S_vrr = QCY*I_ERI_L6x2z_S_S_S_vrr+WQY*I_ERI_L6x2z_S_S_S_M1_vrr;
      Double I_ERI_L5x3y_S_Py_S_vrr = QCY*I_ERI_L5x3y_S_S_S_vrr+WQY*I_ERI_L5x3y_S_S_S_M1_vrr+3*oned2k*I_ERI_K5x2y_S_S_S_M1_vrr;
      Double I_ERI_L5x2yz_S_Py_S_vrr = QCY*I_ERI_L5x2yz_S_S_S_vrr+WQY*I_ERI_L5x2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_K5xyz_S_S_S_M1_vrr;
      Double I_ERI_L5xy2z_S_Py_S_vrr = QCY*I_ERI_L5xy2z_S_S_S_vrr+WQY*I_ERI_L5xy2z_S_S_S_M1_vrr+oned2k*I_ERI_K5x2z_S_S_S_M1_vrr;
      Double I_ERI_L5x3z_S_Py_S_vrr = QCY*I_ERI_L5x3z_S_S_S_vrr+WQY*I_ERI_L5x3z_S_S_S_M1_vrr;
      Double I_ERI_L4x4y_S_Py_S_vrr = QCY*I_ERI_L4x4y_S_S_S_vrr+WQY*I_ERI_L4x4y_S_S_S_M1_vrr+4*oned2k*I_ERI_K4x3y_S_S_S_M1_vrr;
      Double I_ERI_L4x3yz_S_Py_S_vrr = QCY*I_ERI_L4x3yz_S_S_S_vrr+WQY*I_ERI_L4x3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_K4x2yz_S_S_S_M1_vrr;
      Double I_ERI_L4x2y2z_S_Py_S_vrr = QCY*I_ERI_L4x2y2z_S_S_S_vrr+WQY*I_ERI_L4x2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_K4xy2z_S_S_S_M1_vrr;
      Double I_ERI_L4xy3z_S_Py_S_vrr = QCY*I_ERI_L4xy3z_S_S_S_vrr+WQY*I_ERI_L4xy3z_S_S_S_M1_vrr+oned2k*I_ERI_K4x3z_S_S_S_M1_vrr;
      Double I_ERI_L4x4z_S_Py_S_vrr = QCY*I_ERI_L4x4z_S_S_S_vrr+WQY*I_ERI_L4x4z_S_S_S_M1_vrr;
      Double I_ERI_L3x5y_S_Py_S_vrr = QCY*I_ERI_L3x5y_S_S_S_vrr+WQY*I_ERI_L3x5y_S_S_S_M1_vrr+5*oned2k*I_ERI_K3x4y_S_S_S_M1_vrr;
      Double I_ERI_L3x4yz_S_Py_S_vrr = QCY*I_ERI_L3x4yz_S_S_S_vrr+WQY*I_ERI_L3x4yz_S_S_S_M1_vrr+4*oned2k*I_ERI_K3x3yz_S_S_S_M1_vrr;
      Double I_ERI_L3x3y2z_S_Py_S_vrr = QCY*I_ERI_L3x3y2z_S_S_S_vrr+WQY*I_ERI_L3x3y2z_S_S_S_M1_vrr+3*oned2k*I_ERI_K3x2y2z_S_S_S_M1_vrr;
      Double I_ERI_L3x2y3z_S_Py_S_vrr = QCY*I_ERI_L3x2y3z_S_S_S_vrr+WQY*I_ERI_L3x2y3z_S_S_S_M1_vrr+2*oned2k*I_ERI_K3xy3z_S_S_S_M1_vrr;
      Double I_ERI_L3xy4z_S_Py_S_vrr = QCY*I_ERI_L3xy4z_S_S_S_vrr+WQY*I_ERI_L3xy4z_S_S_S_M1_vrr+oned2k*I_ERI_K3x4z_S_S_S_M1_vrr;
      Double I_ERI_L3x5z_S_Py_S_vrr = QCY*I_ERI_L3x5z_S_S_S_vrr+WQY*I_ERI_L3x5z_S_S_S_M1_vrr;
      Double I_ERI_L2x6y_S_Py_S_vrr = QCY*I_ERI_L2x6y_S_S_S_vrr+WQY*I_ERI_L2x6y_S_S_S_M1_vrr+6*oned2k*I_ERI_K2x5y_S_S_S_M1_vrr;
      Double I_ERI_L2x5yz_S_Py_S_vrr = QCY*I_ERI_L2x5yz_S_S_S_vrr+WQY*I_ERI_L2x5yz_S_S_S_M1_vrr+5*oned2k*I_ERI_K2x4yz_S_S_S_M1_vrr;
      Double I_ERI_L2x4y2z_S_Py_S_vrr = QCY*I_ERI_L2x4y2z_S_S_S_vrr+WQY*I_ERI_L2x4y2z_S_S_S_M1_vrr+4*oned2k*I_ERI_K2x3y2z_S_S_S_M1_vrr;
      Double I_ERI_L2x3y3z_S_Py_S_vrr = QCY*I_ERI_L2x3y3z_S_S_S_vrr+WQY*I_ERI_L2x3y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_K2x2y3z_S_S_S_M1_vrr;
      Double I_ERI_L2x2y4z_S_Py_S_vrr = QCY*I_ERI_L2x2y4z_S_S_S_vrr+WQY*I_ERI_L2x2y4z_S_S_S_M1_vrr+2*oned2k*I_ERI_K2xy4z_S_S_S_M1_vrr;
      Double I_ERI_L2xy5z_S_Py_S_vrr = QCY*I_ERI_L2xy5z_S_S_S_vrr+WQY*I_ERI_L2xy5z_S_S_S_M1_vrr+oned2k*I_ERI_K2x5z_S_S_S_M1_vrr;
      Double I_ERI_L2x6z_S_Py_S_vrr = QCY*I_ERI_L2x6z_S_S_S_vrr+WQY*I_ERI_L2x6z_S_S_S_M1_vrr;
      Double I_ERI_Lx7y_S_Py_S_vrr = QCY*I_ERI_Lx7y_S_S_S_vrr+WQY*I_ERI_Lx7y_S_S_S_M1_vrr+7*oned2k*I_ERI_Kx6y_S_S_S_M1_vrr;
      Double I_ERI_Lx6yz_S_Py_S_vrr = QCY*I_ERI_Lx6yz_S_S_S_vrr+WQY*I_ERI_Lx6yz_S_S_S_M1_vrr+6*oned2k*I_ERI_Kx5yz_S_S_S_M1_vrr;
      Double I_ERI_Lx5y2z_S_Py_S_vrr = QCY*I_ERI_Lx5y2z_S_S_S_vrr+WQY*I_ERI_Lx5y2z_S_S_S_M1_vrr+5*oned2k*I_ERI_Kx4y2z_S_S_S_M1_vrr;
      Double I_ERI_Lx4y3z_S_Py_S_vrr = QCY*I_ERI_Lx4y3z_S_S_S_vrr+WQY*I_ERI_Lx4y3z_S_S_S_M1_vrr+4*oned2k*I_ERI_Kx3y3z_S_S_S_M1_vrr;
      Double I_ERI_Lx3y4z_S_Py_S_vrr = QCY*I_ERI_Lx3y4z_S_S_S_vrr+WQY*I_ERI_Lx3y4z_S_S_S_M1_vrr+3*oned2k*I_ERI_Kx2y4z_S_S_S_M1_vrr;
      Double I_ERI_Lx2y5z_S_Py_S_vrr = QCY*I_ERI_Lx2y5z_S_S_S_vrr+WQY*I_ERI_Lx2y5z_S_S_S_M1_vrr+2*oned2k*I_ERI_Kxy5z_S_S_S_M1_vrr;
      Double I_ERI_Lxy6z_S_Py_S_vrr = QCY*I_ERI_Lxy6z_S_S_S_vrr+WQY*I_ERI_Lxy6z_S_S_S_M1_vrr+oned2k*I_ERI_Kx6z_S_S_S_M1_vrr;
      Double I_ERI_Lx7z_S_Py_S_vrr = QCY*I_ERI_Lx7z_S_S_S_vrr+WQY*I_ERI_Lx7z_S_S_S_M1_vrr;
      Double I_ERI_L8y_S_Py_S_vrr = QCY*I_ERI_L8y_S_S_S_vrr+WQY*I_ERI_L8y_S_S_S_M1_vrr+8*oned2k*I_ERI_K7y_S_S_S_M1_vrr;
      Double I_ERI_L7yz_S_Py_S_vrr = QCY*I_ERI_L7yz_S_S_S_vrr+WQY*I_ERI_L7yz_S_S_S_M1_vrr+7*oned2k*I_ERI_K6yz_S_S_S_M1_vrr;
      Double I_ERI_L6y2z_S_Py_S_vrr = QCY*I_ERI_L6y2z_S_S_S_vrr+WQY*I_ERI_L6y2z_S_S_S_M1_vrr+6*oned2k*I_ERI_K5y2z_S_S_S_M1_vrr;
      Double I_ERI_L5y3z_S_Py_S_vrr = QCY*I_ERI_L5y3z_S_S_S_vrr+WQY*I_ERI_L5y3z_S_S_S_M1_vrr+5*oned2k*I_ERI_K4y3z_S_S_S_M1_vrr;
      Double I_ERI_L4y4z_S_Py_S_vrr = QCY*I_ERI_L4y4z_S_S_S_vrr+WQY*I_ERI_L4y4z_S_S_S_M1_vrr+4*oned2k*I_ERI_K3y4z_S_S_S_M1_vrr;
      Double I_ERI_L3y5z_S_Py_S_vrr = QCY*I_ERI_L3y5z_S_S_S_vrr+WQY*I_ERI_L3y5z_S_S_S_M1_vrr+3*oned2k*I_ERI_K2y5z_S_S_S_M1_vrr;
      Double I_ERI_L2y6z_S_Py_S_vrr = QCY*I_ERI_L2y6z_S_S_S_vrr+WQY*I_ERI_L2y6z_S_S_S_M1_vrr+2*oned2k*I_ERI_Ky6z_S_S_S_M1_vrr;
      Double I_ERI_Ly7z_S_Py_S_vrr = QCY*I_ERI_Ly7z_S_S_S_vrr+WQY*I_ERI_Ly7z_S_S_S_M1_vrr+oned2k*I_ERI_K7z_S_S_S_M1_vrr;
      Double I_ERI_L8z_S_Py_S_vrr = QCY*I_ERI_L8z_S_S_S_vrr+WQY*I_ERI_L8z_S_S_S_M1_vrr;
      Double I_ERI_L8x_S_Pz_S_vrr = QCZ*I_ERI_L8x_S_S_S_vrr+WQZ*I_ERI_L8x_S_S_S_M1_vrr;
      Double I_ERI_L7xy_S_Pz_S_vrr = QCZ*I_ERI_L7xy_S_S_S_vrr+WQZ*I_ERI_L7xy_S_S_S_M1_vrr;
      Double I_ERI_L7xz_S_Pz_S_vrr = QCZ*I_ERI_L7xz_S_S_S_vrr+WQZ*I_ERI_L7xz_S_S_S_M1_vrr+oned2k*I_ERI_K7x_S_S_S_M1_vrr;
      Double I_ERI_L6x2y_S_Pz_S_vrr = QCZ*I_ERI_L6x2y_S_S_S_vrr+WQZ*I_ERI_L6x2y_S_S_S_M1_vrr;
      Double I_ERI_L6xyz_S_Pz_S_vrr = QCZ*I_ERI_L6xyz_S_S_S_vrr+WQZ*I_ERI_L6xyz_S_S_S_M1_vrr+oned2k*I_ERI_K6xy_S_S_S_M1_vrr;
      Double I_ERI_L6x2z_S_Pz_S_vrr = QCZ*I_ERI_L6x2z_S_S_S_vrr+WQZ*I_ERI_L6x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_K6xz_S_S_S_M1_vrr;
      Double I_ERI_L5x3y_S_Pz_S_vrr = QCZ*I_ERI_L5x3y_S_S_S_vrr+WQZ*I_ERI_L5x3y_S_S_S_M1_vrr;
      Double I_ERI_L5x2yz_S_Pz_S_vrr = QCZ*I_ERI_L5x2yz_S_S_S_vrr+WQZ*I_ERI_L5x2yz_S_S_S_M1_vrr+oned2k*I_ERI_K5x2y_S_S_S_M1_vrr;
      Double I_ERI_L5xy2z_S_Pz_S_vrr = QCZ*I_ERI_L5xy2z_S_S_S_vrr+WQZ*I_ERI_L5xy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_K5xyz_S_S_S_M1_vrr;
      Double I_ERI_L5x3z_S_Pz_S_vrr = QCZ*I_ERI_L5x3z_S_S_S_vrr+WQZ*I_ERI_L5x3z_S_S_S_M1_vrr+3*oned2k*I_ERI_K5x2z_S_S_S_M1_vrr;
      Double I_ERI_L4x4y_S_Pz_S_vrr = QCZ*I_ERI_L4x4y_S_S_S_vrr+WQZ*I_ERI_L4x4y_S_S_S_M1_vrr;
      Double I_ERI_L4x3yz_S_Pz_S_vrr = QCZ*I_ERI_L4x3yz_S_S_S_vrr+WQZ*I_ERI_L4x3yz_S_S_S_M1_vrr+oned2k*I_ERI_K4x3y_S_S_S_M1_vrr;
      Double I_ERI_L4x2y2z_S_Pz_S_vrr = QCZ*I_ERI_L4x2y2z_S_S_S_vrr+WQZ*I_ERI_L4x2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_K4x2yz_S_S_S_M1_vrr;
      Double I_ERI_L4xy3z_S_Pz_S_vrr = QCZ*I_ERI_L4xy3z_S_S_S_vrr+WQZ*I_ERI_L4xy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_K4xy2z_S_S_S_M1_vrr;
      Double I_ERI_L4x4z_S_Pz_S_vrr = QCZ*I_ERI_L4x4z_S_S_S_vrr+WQZ*I_ERI_L4x4z_S_S_S_M1_vrr+4*oned2k*I_ERI_K4x3z_S_S_S_M1_vrr;
      Double I_ERI_L3x5y_S_Pz_S_vrr = QCZ*I_ERI_L3x5y_S_S_S_vrr+WQZ*I_ERI_L3x5y_S_S_S_M1_vrr;
      Double I_ERI_L3x4yz_S_Pz_S_vrr = QCZ*I_ERI_L3x4yz_S_S_S_vrr+WQZ*I_ERI_L3x4yz_S_S_S_M1_vrr+oned2k*I_ERI_K3x4y_S_S_S_M1_vrr;
      Double I_ERI_L3x3y2z_S_Pz_S_vrr = QCZ*I_ERI_L3x3y2z_S_S_S_vrr+WQZ*I_ERI_L3x3y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_K3x3yz_S_S_S_M1_vrr;
      Double I_ERI_L3x2y3z_S_Pz_S_vrr = QCZ*I_ERI_L3x2y3z_S_S_S_vrr+WQZ*I_ERI_L3x2y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_K3x2y2z_S_S_S_M1_vrr;
      Double I_ERI_L3xy4z_S_Pz_S_vrr = QCZ*I_ERI_L3xy4z_S_S_S_vrr+WQZ*I_ERI_L3xy4z_S_S_S_M1_vrr+4*oned2k*I_ERI_K3xy3z_S_S_S_M1_vrr;
      Double I_ERI_L3x5z_S_Pz_S_vrr = QCZ*I_ERI_L3x5z_S_S_S_vrr+WQZ*I_ERI_L3x5z_S_S_S_M1_vrr+5*oned2k*I_ERI_K3x4z_S_S_S_M1_vrr;
      Double I_ERI_L2x6y_S_Pz_S_vrr = QCZ*I_ERI_L2x6y_S_S_S_vrr+WQZ*I_ERI_L2x6y_S_S_S_M1_vrr;
      Double I_ERI_L2x5yz_S_Pz_S_vrr = QCZ*I_ERI_L2x5yz_S_S_S_vrr+WQZ*I_ERI_L2x5yz_S_S_S_M1_vrr+oned2k*I_ERI_K2x5y_S_S_S_M1_vrr;
      Double I_ERI_L2x4y2z_S_Pz_S_vrr = QCZ*I_ERI_L2x4y2z_S_S_S_vrr+WQZ*I_ERI_L2x4y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_K2x4yz_S_S_S_M1_vrr;
      Double I_ERI_L2x3y3z_S_Pz_S_vrr = QCZ*I_ERI_L2x3y3z_S_S_S_vrr+WQZ*I_ERI_L2x3y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_K2x3y2z_S_S_S_M1_vrr;
      Double I_ERI_L2x2y4z_S_Pz_S_vrr = QCZ*I_ERI_L2x2y4z_S_S_S_vrr+WQZ*I_ERI_L2x2y4z_S_S_S_M1_vrr+4*oned2k*I_ERI_K2x2y3z_S_S_S_M1_vrr;
      Double I_ERI_L2xy5z_S_Pz_S_vrr = QCZ*I_ERI_L2xy5z_S_S_S_vrr+WQZ*I_ERI_L2xy5z_S_S_S_M1_vrr+5*oned2k*I_ERI_K2xy4z_S_S_S_M1_vrr;
      Double I_ERI_L2x6z_S_Pz_S_vrr = QCZ*I_ERI_L2x6z_S_S_S_vrr+WQZ*I_ERI_L2x6z_S_S_S_M1_vrr+6*oned2k*I_ERI_K2x5z_S_S_S_M1_vrr;
      Double I_ERI_Lx7y_S_Pz_S_vrr = QCZ*I_ERI_Lx7y_S_S_S_vrr+WQZ*I_ERI_Lx7y_S_S_S_M1_vrr;
      Double I_ERI_Lx6yz_S_Pz_S_vrr = QCZ*I_ERI_Lx6yz_S_S_S_vrr+WQZ*I_ERI_Lx6yz_S_S_S_M1_vrr+oned2k*I_ERI_Kx6y_S_S_S_M1_vrr;
      Double I_ERI_Lx5y2z_S_Pz_S_vrr = QCZ*I_ERI_Lx5y2z_S_S_S_vrr+WQZ*I_ERI_Lx5y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Kx5yz_S_S_S_M1_vrr;
      Double I_ERI_Lx4y3z_S_Pz_S_vrr = QCZ*I_ERI_Lx4y3z_S_S_S_vrr+WQZ*I_ERI_Lx4y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Kx4y2z_S_S_S_M1_vrr;
      Double I_ERI_Lx3y4z_S_Pz_S_vrr = QCZ*I_ERI_Lx3y4z_S_S_S_vrr+WQZ*I_ERI_Lx3y4z_S_S_S_M1_vrr+4*oned2k*I_ERI_Kx3y3z_S_S_S_M1_vrr;
      Double I_ERI_Lx2y5z_S_Pz_S_vrr = QCZ*I_ERI_Lx2y5z_S_S_S_vrr+WQZ*I_ERI_Lx2y5z_S_S_S_M1_vrr+5*oned2k*I_ERI_Kx2y4z_S_S_S_M1_vrr;
      Double I_ERI_Lxy6z_S_Pz_S_vrr = QCZ*I_ERI_Lxy6z_S_S_S_vrr+WQZ*I_ERI_Lxy6z_S_S_S_M1_vrr+6*oned2k*I_ERI_Kxy5z_S_S_S_M1_vrr;
      Double I_ERI_Lx7z_S_Pz_S_vrr = QCZ*I_ERI_Lx7z_S_S_S_vrr+WQZ*I_ERI_Lx7z_S_S_S_M1_vrr+7*oned2k*I_ERI_Kx6z_S_S_S_M1_vrr;
      Double I_ERI_L8y_S_Pz_S_vrr = QCZ*I_ERI_L8y_S_S_S_vrr+WQZ*I_ERI_L8y_S_S_S_M1_vrr;
      Double I_ERI_L7yz_S_Pz_S_vrr = QCZ*I_ERI_L7yz_S_S_S_vrr+WQZ*I_ERI_L7yz_S_S_S_M1_vrr+oned2k*I_ERI_K7y_S_S_S_M1_vrr;
      Double I_ERI_L6y2z_S_Pz_S_vrr = QCZ*I_ERI_L6y2z_S_S_S_vrr+WQZ*I_ERI_L6y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_K6yz_S_S_S_M1_vrr;
      Double I_ERI_L5y3z_S_Pz_S_vrr = QCZ*I_ERI_L5y3z_S_S_S_vrr+WQZ*I_ERI_L5y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_K5y2z_S_S_S_M1_vrr;
      Double I_ERI_L4y4z_S_Pz_S_vrr = QCZ*I_ERI_L4y4z_S_S_S_vrr+WQZ*I_ERI_L4y4z_S_S_S_M1_vrr+4*oned2k*I_ERI_K4y3z_S_S_S_M1_vrr;
      Double I_ERI_L3y5z_S_Pz_S_vrr = QCZ*I_ERI_L3y5z_S_S_S_vrr+WQZ*I_ERI_L3y5z_S_S_S_M1_vrr+5*oned2k*I_ERI_K3y4z_S_S_S_M1_vrr;
      Double I_ERI_L2y6z_S_Pz_S_vrr = QCZ*I_ERI_L2y6z_S_S_S_vrr+WQZ*I_ERI_L2y6z_S_S_S_M1_vrr+6*oned2k*I_ERI_K2y5z_S_S_S_M1_vrr;
      Double I_ERI_Ly7z_S_Pz_S_vrr = QCZ*I_ERI_Ly7z_S_S_S_vrr+WQZ*I_ERI_Ly7z_S_S_S_M1_vrr+7*oned2k*I_ERI_Ky6z_S_S_S_M1_vrr;
      Double I_ERI_L8z_S_Pz_S_vrr = QCZ*I_ERI_L8z_S_S_S_vrr+WQZ*I_ERI_L8z_S_S_S_M1_vrr+8*oned2k*I_ERI_K7z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_L_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_L8x_S_Px_S += I_ERI_L8x_S_Px_S_vrr;
      I_ERI_L7xy_S_Px_S += I_ERI_L7xy_S_Px_S_vrr;
      I_ERI_L7xz_S_Px_S += I_ERI_L7xz_S_Px_S_vrr;
      I_ERI_L6x2y_S_Px_S += I_ERI_L6x2y_S_Px_S_vrr;
      I_ERI_L6xyz_S_Px_S += I_ERI_L6xyz_S_Px_S_vrr;
      I_ERI_L6x2z_S_Px_S += I_ERI_L6x2z_S_Px_S_vrr;
      I_ERI_L5x3y_S_Px_S += I_ERI_L5x3y_S_Px_S_vrr;
      I_ERI_L5x2yz_S_Px_S += I_ERI_L5x2yz_S_Px_S_vrr;
      I_ERI_L5xy2z_S_Px_S += I_ERI_L5xy2z_S_Px_S_vrr;
      I_ERI_L5x3z_S_Px_S += I_ERI_L5x3z_S_Px_S_vrr;
      I_ERI_L4x4y_S_Px_S += I_ERI_L4x4y_S_Px_S_vrr;
      I_ERI_L4x3yz_S_Px_S += I_ERI_L4x3yz_S_Px_S_vrr;
      I_ERI_L4x2y2z_S_Px_S += I_ERI_L4x2y2z_S_Px_S_vrr;
      I_ERI_L4xy3z_S_Px_S += I_ERI_L4xy3z_S_Px_S_vrr;
      I_ERI_L4x4z_S_Px_S += I_ERI_L4x4z_S_Px_S_vrr;
      I_ERI_L3x5y_S_Px_S += I_ERI_L3x5y_S_Px_S_vrr;
      I_ERI_L3x4yz_S_Px_S += I_ERI_L3x4yz_S_Px_S_vrr;
      I_ERI_L3x3y2z_S_Px_S += I_ERI_L3x3y2z_S_Px_S_vrr;
      I_ERI_L3x2y3z_S_Px_S += I_ERI_L3x2y3z_S_Px_S_vrr;
      I_ERI_L3xy4z_S_Px_S += I_ERI_L3xy4z_S_Px_S_vrr;
      I_ERI_L3x5z_S_Px_S += I_ERI_L3x5z_S_Px_S_vrr;
      I_ERI_L2x6y_S_Px_S += I_ERI_L2x6y_S_Px_S_vrr;
      I_ERI_L2x5yz_S_Px_S += I_ERI_L2x5yz_S_Px_S_vrr;
      I_ERI_L2x4y2z_S_Px_S += I_ERI_L2x4y2z_S_Px_S_vrr;
      I_ERI_L2x3y3z_S_Px_S += I_ERI_L2x3y3z_S_Px_S_vrr;
      I_ERI_L2x2y4z_S_Px_S += I_ERI_L2x2y4z_S_Px_S_vrr;
      I_ERI_L2xy5z_S_Px_S += I_ERI_L2xy5z_S_Px_S_vrr;
      I_ERI_L2x6z_S_Px_S += I_ERI_L2x6z_S_Px_S_vrr;
      I_ERI_Lx7y_S_Px_S += I_ERI_Lx7y_S_Px_S_vrr;
      I_ERI_Lx6yz_S_Px_S += I_ERI_Lx6yz_S_Px_S_vrr;
      I_ERI_Lx5y2z_S_Px_S += I_ERI_Lx5y2z_S_Px_S_vrr;
      I_ERI_Lx4y3z_S_Px_S += I_ERI_Lx4y3z_S_Px_S_vrr;
      I_ERI_Lx3y4z_S_Px_S += I_ERI_Lx3y4z_S_Px_S_vrr;
      I_ERI_Lx2y5z_S_Px_S += I_ERI_Lx2y5z_S_Px_S_vrr;
      I_ERI_Lxy6z_S_Px_S += I_ERI_Lxy6z_S_Px_S_vrr;
      I_ERI_Lx7z_S_Px_S += I_ERI_Lx7z_S_Px_S_vrr;
      I_ERI_L8y_S_Px_S += I_ERI_L8y_S_Px_S_vrr;
      I_ERI_L7yz_S_Px_S += I_ERI_L7yz_S_Px_S_vrr;
      I_ERI_L6y2z_S_Px_S += I_ERI_L6y2z_S_Px_S_vrr;
      I_ERI_L5y3z_S_Px_S += I_ERI_L5y3z_S_Px_S_vrr;
      I_ERI_L4y4z_S_Px_S += I_ERI_L4y4z_S_Px_S_vrr;
      I_ERI_L3y5z_S_Px_S += I_ERI_L3y5z_S_Px_S_vrr;
      I_ERI_L2y6z_S_Px_S += I_ERI_L2y6z_S_Px_S_vrr;
      I_ERI_Ly7z_S_Px_S += I_ERI_Ly7z_S_Px_S_vrr;
      I_ERI_L8z_S_Px_S += I_ERI_L8z_S_Px_S_vrr;
      I_ERI_L8x_S_Py_S += I_ERI_L8x_S_Py_S_vrr;
      I_ERI_L7xy_S_Py_S += I_ERI_L7xy_S_Py_S_vrr;
      I_ERI_L7xz_S_Py_S += I_ERI_L7xz_S_Py_S_vrr;
      I_ERI_L6x2y_S_Py_S += I_ERI_L6x2y_S_Py_S_vrr;
      I_ERI_L6xyz_S_Py_S += I_ERI_L6xyz_S_Py_S_vrr;
      I_ERI_L6x2z_S_Py_S += I_ERI_L6x2z_S_Py_S_vrr;
      I_ERI_L5x3y_S_Py_S += I_ERI_L5x3y_S_Py_S_vrr;
      I_ERI_L5x2yz_S_Py_S += I_ERI_L5x2yz_S_Py_S_vrr;
      I_ERI_L5xy2z_S_Py_S += I_ERI_L5xy2z_S_Py_S_vrr;
      I_ERI_L5x3z_S_Py_S += I_ERI_L5x3z_S_Py_S_vrr;
      I_ERI_L4x4y_S_Py_S += I_ERI_L4x4y_S_Py_S_vrr;
      I_ERI_L4x3yz_S_Py_S += I_ERI_L4x3yz_S_Py_S_vrr;
      I_ERI_L4x2y2z_S_Py_S += I_ERI_L4x2y2z_S_Py_S_vrr;
      I_ERI_L4xy3z_S_Py_S += I_ERI_L4xy3z_S_Py_S_vrr;
      I_ERI_L4x4z_S_Py_S += I_ERI_L4x4z_S_Py_S_vrr;
      I_ERI_L3x5y_S_Py_S += I_ERI_L3x5y_S_Py_S_vrr;
      I_ERI_L3x4yz_S_Py_S += I_ERI_L3x4yz_S_Py_S_vrr;
      I_ERI_L3x3y2z_S_Py_S += I_ERI_L3x3y2z_S_Py_S_vrr;
      I_ERI_L3x2y3z_S_Py_S += I_ERI_L3x2y3z_S_Py_S_vrr;
      I_ERI_L3xy4z_S_Py_S += I_ERI_L3xy4z_S_Py_S_vrr;
      I_ERI_L3x5z_S_Py_S += I_ERI_L3x5z_S_Py_S_vrr;
      I_ERI_L2x6y_S_Py_S += I_ERI_L2x6y_S_Py_S_vrr;
      I_ERI_L2x5yz_S_Py_S += I_ERI_L2x5yz_S_Py_S_vrr;
      I_ERI_L2x4y2z_S_Py_S += I_ERI_L2x4y2z_S_Py_S_vrr;
      I_ERI_L2x3y3z_S_Py_S += I_ERI_L2x3y3z_S_Py_S_vrr;
      I_ERI_L2x2y4z_S_Py_S += I_ERI_L2x2y4z_S_Py_S_vrr;
      I_ERI_L2xy5z_S_Py_S += I_ERI_L2xy5z_S_Py_S_vrr;
      I_ERI_L2x6z_S_Py_S += I_ERI_L2x6z_S_Py_S_vrr;
      I_ERI_Lx7y_S_Py_S += I_ERI_Lx7y_S_Py_S_vrr;
      I_ERI_Lx6yz_S_Py_S += I_ERI_Lx6yz_S_Py_S_vrr;
      I_ERI_Lx5y2z_S_Py_S += I_ERI_Lx5y2z_S_Py_S_vrr;
      I_ERI_Lx4y3z_S_Py_S += I_ERI_Lx4y3z_S_Py_S_vrr;
      I_ERI_Lx3y4z_S_Py_S += I_ERI_Lx3y4z_S_Py_S_vrr;
      I_ERI_Lx2y5z_S_Py_S += I_ERI_Lx2y5z_S_Py_S_vrr;
      I_ERI_Lxy6z_S_Py_S += I_ERI_Lxy6z_S_Py_S_vrr;
      I_ERI_Lx7z_S_Py_S += I_ERI_Lx7z_S_Py_S_vrr;
      I_ERI_L8y_S_Py_S += I_ERI_L8y_S_Py_S_vrr;
      I_ERI_L7yz_S_Py_S += I_ERI_L7yz_S_Py_S_vrr;
      I_ERI_L6y2z_S_Py_S += I_ERI_L6y2z_S_Py_S_vrr;
      I_ERI_L5y3z_S_Py_S += I_ERI_L5y3z_S_Py_S_vrr;
      I_ERI_L4y4z_S_Py_S += I_ERI_L4y4z_S_Py_S_vrr;
      I_ERI_L3y5z_S_Py_S += I_ERI_L3y5z_S_Py_S_vrr;
      I_ERI_L2y6z_S_Py_S += I_ERI_L2y6z_S_Py_S_vrr;
      I_ERI_Ly7z_S_Py_S += I_ERI_Ly7z_S_Py_S_vrr;
      I_ERI_L8z_S_Py_S += I_ERI_L8z_S_Py_S_vrr;
      I_ERI_L8x_S_Pz_S += I_ERI_L8x_S_Pz_S_vrr;
      I_ERI_L7xy_S_Pz_S += I_ERI_L7xy_S_Pz_S_vrr;
      I_ERI_L7xz_S_Pz_S += I_ERI_L7xz_S_Pz_S_vrr;
      I_ERI_L6x2y_S_Pz_S += I_ERI_L6x2y_S_Pz_S_vrr;
      I_ERI_L6xyz_S_Pz_S += I_ERI_L6xyz_S_Pz_S_vrr;
      I_ERI_L6x2z_S_Pz_S += I_ERI_L6x2z_S_Pz_S_vrr;
      I_ERI_L5x3y_S_Pz_S += I_ERI_L5x3y_S_Pz_S_vrr;
      I_ERI_L5x2yz_S_Pz_S += I_ERI_L5x2yz_S_Pz_S_vrr;
      I_ERI_L5xy2z_S_Pz_S += I_ERI_L5xy2z_S_Pz_S_vrr;
      I_ERI_L5x3z_S_Pz_S += I_ERI_L5x3z_S_Pz_S_vrr;
      I_ERI_L4x4y_S_Pz_S += I_ERI_L4x4y_S_Pz_S_vrr;
      I_ERI_L4x3yz_S_Pz_S += I_ERI_L4x3yz_S_Pz_S_vrr;
      I_ERI_L4x2y2z_S_Pz_S += I_ERI_L4x2y2z_S_Pz_S_vrr;
      I_ERI_L4xy3z_S_Pz_S += I_ERI_L4xy3z_S_Pz_S_vrr;
      I_ERI_L4x4z_S_Pz_S += I_ERI_L4x4z_S_Pz_S_vrr;
      I_ERI_L3x5y_S_Pz_S += I_ERI_L3x5y_S_Pz_S_vrr;
      I_ERI_L3x4yz_S_Pz_S += I_ERI_L3x4yz_S_Pz_S_vrr;
      I_ERI_L3x3y2z_S_Pz_S += I_ERI_L3x3y2z_S_Pz_S_vrr;
      I_ERI_L3x2y3z_S_Pz_S += I_ERI_L3x2y3z_S_Pz_S_vrr;
      I_ERI_L3xy4z_S_Pz_S += I_ERI_L3xy4z_S_Pz_S_vrr;
      I_ERI_L3x5z_S_Pz_S += I_ERI_L3x5z_S_Pz_S_vrr;
      I_ERI_L2x6y_S_Pz_S += I_ERI_L2x6y_S_Pz_S_vrr;
      I_ERI_L2x5yz_S_Pz_S += I_ERI_L2x5yz_S_Pz_S_vrr;
      I_ERI_L2x4y2z_S_Pz_S += I_ERI_L2x4y2z_S_Pz_S_vrr;
      I_ERI_L2x3y3z_S_Pz_S += I_ERI_L2x3y3z_S_Pz_S_vrr;
      I_ERI_L2x2y4z_S_Pz_S += I_ERI_L2x2y4z_S_Pz_S_vrr;
      I_ERI_L2xy5z_S_Pz_S += I_ERI_L2xy5z_S_Pz_S_vrr;
      I_ERI_L2x6z_S_Pz_S += I_ERI_L2x6z_S_Pz_S_vrr;
      I_ERI_Lx7y_S_Pz_S += I_ERI_Lx7y_S_Pz_S_vrr;
      I_ERI_Lx6yz_S_Pz_S += I_ERI_Lx6yz_S_Pz_S_vrr;
      I_ERI_Lx5y2z_S_Pz_S += I_ERI_Lx5y2z_S_Pz_S_vrr;
      I_ERI_Lx4y3z_S_Pz_S += I_ERI_Lx4y3z_S_Pz_S_vrr;
      I_ERI_Lx3y4z_S_Pz_S += I_ERI_Lx3y4z_S_Pz_S_vrr;
      I_ERI_Lx2y5z_S_Pz_S += I_ERI_Lx2y5z_S_Pz_S_vrr;
      I_ERI_Lxy6z_S_Pz_S += I_ERI_Lxy6z_S_Pz_S_vrr;
      I_ERI_Lx7z_S_Pz_S += I_ERI_Lx7z_S_Pz_S_vrr;
      I_ERI_L8y_S_Pz_S += I_ERI_L8y_S_Pz_S_vrr;
      I_ERI_L7yz_S_Pz_S += I_ERI_L7yz_S_Pz_S_vrr;
      I_ERI_L6y2z_S_Pz_S += I_ERI_L6y2z_S_Pz_S_vrr;
      I_ERI_L5y3z_S_Pz_S += I_ERI_L5y3z_S_Pz_S_vrr;
      I_ERI_L4y4z_S_Pz_S += I_ERI_L4y4z_S_Pz_S_vrr;
      I_ERI_L3y5z_S_Pz_S += I_ERI_L3y5z_S_Pz_S_vrr;
      I_ERI_L2y6z_S_Pz_S += I_ERI_L2y6z_S_Pz_S_vrr;
      I_ERI_Ly7z_S_Pz_S += I_ERI_Ly7z_S_Pz_S_vrr;
      I_ERI_L8z_S_Pz_S += I_ERI_L8z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_K_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_K7x_S_Px_S += I_ERI_K7x_S_Px_S_vrr;
      I_ERI_K6xy_S_Px_S += I_ERI_K6xy_S_Px_S_vrr;
      I_ERI_K6xz_S_Px_S += I_ERI_K6xz_S_Px_S_vrr;
      I_ERI_K5x2y_S_Px_S += I_ERI_K5x2y_S_Px_S_vrr;
      I_ERI_K5xyz_S_Px_S += I_ERI_K5xyz_S_Px_S_vrr;
      I_ERI_K5x2z_S_Px_S += I_ERI_K5x2z_S_Px_S_vrr;
      I_ERI_K4x3y_S_Px_S += I_ERI_K4x3y_S_Px_S_vrr;
      I_ERI_K4x2yz_S_Px_S += I_ERI_K4x2yz_S_Px_S_vrr;
      I_ERI_K4xy2z_S_Px_S += I_ERI_K4xy2z_S_Px_S_vrr;
      I_ERI_K4x3z_S_Px_S += I_ERI_K4x3z_S_Px_S_vrr;
      I_ERI_K3x4y_S_Px_S += I_ERI_K3x4y_S_Px_S_vrr;
      I_ERI_K3x3yz_S_Px_S += I_ERI_K3x3yz_S_Px_S_vrr;
      I_ERI_K3x2y2z_S_Px_S += I_ERI_K3x2y2z_S_Px_S_vrr;
      I_ERI_K3xy3z_S_Px_S += I_ERI_K3xy3z_S_Px_S_vrr;
      I_ERI_K3x4z_S_Px_S += I_ERI_K3x4z_S_Px_S_vrr;
      I_ERI_K2x5y_S_Px_S += I_ERI_K2x5y_S_Px_S_vrr;
      I_ERI_K2x4yz_S_Px_S += I_ERI_K2x4yz_S_Px_S_vrr;
      I_ERI_K2x3y2z_S_Px_S += I_ERI_K2x3y2z_S_Px_S_vrr;
      I_ERI_K2x2y3z_S_Px_S += I_ERI_K2x2y3z_S_Px_S_vrr;
      I_ERI_K2xy4z_S_Px_S += I_ERI_K2xy4z_S_Px_S_vrr;
      I_ERI_K2x5z_S_Px_S += I_ERI_K2x5z_S_Px_S_vrr;
      I_ERI_Kx6y_S_Px_S += I_ERI_Kx6y_S_Px_S_vrr;
      I_ERI_Kx5yz_S_Px_S += I_ERI_Kx5yz_S_Px_S_vrr;
      I_ERI_Kx4y2z_S_Px_S += I_ERI_Kx4y2z_S_Px_S_vrr;
      I_ERI_Kx3y3z_S_Px_S += I_ERI_Kx3y3z_S_Px_S_vrr;
      I_ERI_Kx2y4z_S_Px_S += I_ERI_Kx2y4z_S_Px_S_vrr;
      I_ERI_Kxy5z_S_Px_S += I_ERI_Kxy5z_S_Px_S_vrr;
      I_ERI_Kx6z_S_Px_S += I_ERI_Kx6z_S_Px_S_vrr;
      I_ERI_K7y_S_Px_S += I_ERI_K7y_S_Px_S_vrr;
      I_ERI_K6yz_S_Px_S += I_ERI_K6yz_S_Px_S_vrr;
      I_ERI_K5y2z_S_Px_S += I_ERI_K5y2z_S_Px_S_vrr;
      I_ERI_K4y3z_S_Px_S += I_ERI_K4y3z_S_Px_S_vrr;
      I_ERI_K3y4z_S_Px_S += I_ERI_K3y4z_S_Px_S_vrr;
      I_ERI_K2y5z_S_Px_S += I_ERI_K2y5z_S_Px_S_vrr;
      I_ERI_Ky6z_S_Px_S += I_ERI_Ky6z_S_Px_S_vrr;
      I_ERI_K7z_S_Px_S += I_ERI_K7z_S_Px_S_vrr;
      I_ERI_K7x_S_Py_S += I_ERI_K7x_S_Py_S_vrr;
      I_ERI_K6xy_S_Py_S += I_ERI_K6xy_S_Py_S_vrr;
      I_ERI_K6xz_S_Py_S += I_ERI_K6xz_S_Py_S_vrr;
      I_ERI_K5x2y_S_Py_S += I_ERI_K5x2y_S_Py_S_vrr;
      I_ERI_K5xyz_S_Py_S += I_ERI_K5xyz_S_Py_S_vrr;
      I_ERI_K5x2z_S_Py_S += I_ERI_K5x2z_S_Py_S_vrr;
      I_ERI_K4x3y_S_Py_S += I_ERI_K4x3y_S_Py_S_vrr;
      I_ERI_K4x2yz_S_Py_S += I_ERI_K4x2yz_S_Py_S_vrr;
      I_ERI_K4xy2z_S_Py_S += I_ERI_K4xy2z_S_Py_S_vrr;
      I_ERI_K4x3z_S_Py_S += I_ERI_K4x3z_S_Py_S_vrr;
      I_ERI_K3x4y_S_Py_S += I_ERI_K3x4y_S_Py_S_vrr;
      I_ERI_K3x3yz_S_Py_S += I_ERI_K3x3yz_S_Py_S_vrr;
      I_ERI_K3x2y2z_S_Py_S += I_ERI_K3x2y2z_S_Py_S_vrr;
      I_ERI_K3xy3z_S_Py_S += I_ERI_K3xy3z_S_Py_S_vrr;
      I_ERI_K3x4z_S_Py_S += I_ERI_K3x4z_S_Py_S_vrr;
      I_ERI_K2x5y_S_Py_S += I_ERI_K2x5y_S_Py_S_vrr;
      I_ERI_K2x4yz_S_Py_S += I_ERI_K2x4yz_S_Py_S_vrr;
      I_ERI_K2x3y2z_S_Py_S += I_ERI_K2x3y2z_S_Py_S_vrr;
      I_ERI_K2x2y3z_S_Py_S += I_ERI_K2x2y3z_S_Py_S_vrr;
      I_ERI_K2xy4z_S_Py_S += I_ERI_K2xy4z_S_Py_S_vrr;
      I_ERI_K2x5z_S_Py_S += I_ERI_K2x5z_S_Py_S_vrr;
      I_ERI_Kx6y_S_Py_S += I_ERI_Kx6y_S_Py_S_vrr;
      I_ERI_Kx5yz_S_Py_S += I_ERI_Kx5yz_S_Py_S_vrr;
      I_ERI_Kx4y2z_S_Py_S += I_ERI_Kx4y2z_S_Py_S_vrr;
      I_ERI_Kx3y3z_S_Py_S += I_ERI_Kx3y3z_S_Py_S_vrr;
      I_ERI_Kx2y4z_S_Py_S += I_ERI_Kx2y4z_S_Py_S_vrr;
      I_ERI_Kxy5z_S_Py_S += I_ERI_Kxy5z_S_Py_S_vrr;
      I_ERI_Kx6z_S_Py_S += I_ERI_Kx6z_S_Py_S_vrr;
      I_ERI_K7y_S_Py_S += I_ERI_K7y_S_Py_S_vrr;
      I_ERI_K6yz_S_Py_S += I_ERI_K6yz_S_Py_S_vrr;
      I_ERI_K5y2z_S_Py_S += I_ERI_K5y2z_S_Py_S_vrr;
      I_ERI_K4y3z_S_Py_S += I_ERI_K4y3z_S_Py_S_vrr;
      I_ERI_K3y4z_S_Py_S += I_ERI_K3y4z_S_Py_S_vrr;
      I_ERI_K2y5z_S_Py_S += I_ERI_K2y5z_S_Py_S_vrr;
      I_ERI_Ky6z_S_Py_S += I_ERI_Ky6z_S_Py_S_vrr;
      I_ERI_K7z_S_Py_S += I_ERI_K7z_S_Py_S_vrr;
      I_ERI_K7x_S_Pz_S += I_ERI_K7x_S_Pz_S_vrr;
      I_ERI_K6xy_S_Pz_S += I_ERI_K6xy_S_Pz_S_vrr;
      I_ERI_K6xz_S_Pz_S += I_ERI_K6xz_S_Pz_S_vrr;
      I_ERI_K5x2y_S_Pz_S += I_ERI_K5x2y_S_Pz_S_vrr;
      I_ERI_K5xyz_S_Pz_S += I_ERI_K5xyz_S_Pz_S_vrr;
      I_ERI_K5x2z_S_Pz_S += I_ERI_K5x2z_S_Pz_S_vrr;
      I_ERI_K4x3y_S_Pz_S += I_ERI_K4x3y_S_Pz_S_vrr;
      I_ERI_K4x2yz_S_Pz_S += I_ERI_K4x2yz_S_Pz_S_vrr;
      I_ERI_K4xy2z_S_Pz_S += I_ERI_K4xy2z_S_Pz_S_vrr;
      I_ERI_K4x3z_S_Pz_S += I_ERI_K4x3z_S_Pz_S_vrr;
      I_ERI_K3x4y_S_Pz_S += I_ERI_K3x4y_S_Pz_S_vrr;
      I_ERI_K3x3yz_S_Pz_S += I_ERI_K3x3yz_S_Pz_S_vrr;
      I_ERI_K3x2y2z_S_Pz_S += I_ERI_K3x2y2z_S_Pz_S_vrr;
      I_ERI_K3xy3z_S_Pz_S += I_ERI_K3xy3z_S_Pz_S_vrr;
      I_ERI_K3x4z_S_Pz_S += I_ERI_K3x4z_S_Pz_S_vrr;
      I_ERI_K2x5y_S_Pz_S += I_ERI_K2x5y_S_Pz_S_vrr;
      I_ERI_K2x4yz_S_Pz_S += I_ERI_K2x4yz_S_Pz_S_vrr;
      I_ERI_K2x3y2z_S_Pz_S += I_ERI_K2x3y2z_S_Pz_S_vrr;
      I_ERI_K2x2y3z_S_Pz_S += I_ERI_K2x2y3z_S_Pz_S_vrr;
      I_ERI_K2xy4z_S_Pz_S += I_ERI_K2xy4z_S_Pz_S_vrr;
      I_ERI_K2x5z_S_Pz_S += I_ERI_K2x5z_S_Pz_S_vrr;
      I_ERI_Kx6y_S_Pz_S += I_ERI_Kx6y_S_Pz_S_vrr;
      I_ERI_Kx5yz_S_Pz_S += I_ERI_Kx5yz_S_Pz_S_vrr;
      I_ERI_Kx4y2z_S_Pz_S += I_ERI_Kx4y2z_S_Pz_S_vrr;
      I_ERI_Kx3y3z_S_Pz_S += I_ERI_Kx3y3z_S_Pz_S_vrr;
      I_ERI_Kx2y4z_S_Pz_S += I_ERI_Kx2y4z_S_Pz_S_vrr;
      I_ERI_Kxy5z_S_Pz_S += I_ERI_Kxy5z_S_Pz_S_vrr;
      I_ERI_Kx6z_S_Pz_S += I_ERI_Kx6z_S_Pz_S_vrr;
      I_ERI_K7y_S_Pz_S += I_ERI_K7y_S_Pz_S_vrr;
      I_ERI_K6yz_S_Pz_S += I_ERI_K6yz_S_Pz_S_vrr;
      I_ERI_K5y2z_S_Pz_S += I_ERI_K5y2z_S_Pz_S_vrr;
      I_ERI_K4y3z_S_Pz_S += I_ERI_K4y3z_S_Pz_S_vrr;
      I_ERI_K3y4z_S_Pz_S += I_ERI_K3y4z_S_Pz_S_vrr;
      I_ERI_K2y5z_S_Pz_S += I_ERI_K2y5z_S_Pz_S_vrr;
      I_ERI_Ky6z_S_Pz_S += I_ERI_Ky6z_S_Pz_S_vrr;
      I_ERI_K7z_S_Pz_S += I_ERI_K7z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_I_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_I6x_S_Px_S += I_ERI_I6x_S_Px_S_vrr;
      I_ERI_I5xy_S_Px_S += I_ERI_I5xy_S_Px_S_vrr;
      I_ERI_I5xz_S_Px_S += I_ERI_I5xz_S_Px_S_vrr;
      I_ERI_I4x2y_S_Px_S += I_ERI_I4x2y_S_Px_S_vrr;
      I_ERI_I4xyz_S_Px_S += I_ERI_I4xyz_S_Px_S_vrr;
      I_ERI_I4x2z_S_Px_S += I_ERI_I4x2z_S_Px_S_vrr;
      I_ERI_I3x3y_S_Px_S += I_ERI_I3x3y_S_Px_S_vrr;
      I_ERI_I3x2yz_S_Px_S += I_ERI_I3x2yz_S_Px_S_vrr;
      I_ERI_I3xy2z_S_Px_S += I_ERI_I3xy2z_S_Px_S_vrr;
      I_ERI_I3x3z_S_Px_S += I_ERI_I3x3z_S_Px_S_vrr;
      I_ERI_I2x4y_S_Px_S += I_ERI_I2x4y_S_Px_S_vrr;
      I_ERI_I2x3yz_S_Px_S += I_ERI_I2x3yz_S_Px_S_vrr;
      I_ERI_I2x2y2z_S_Px_S += I_ERI_I2x2y2z_S_Px_S_vrr;
      I_ERI_I2xy3z_S_Px_S += I_ERI_I2xy3z_S_Px_S_vrr;
      I_ERI_I2x4z_S_Px_S += I_ERI_I2x4z_S_Px_S_vrr;
      I_ERI_Ix5y_S_Px_S += I_ERI_Ix5y_S_Px_S_vrr;
      I_ERI_Ix4yz_S_Px_S += I_ERI_Ix4yz_S_Px_S_vrr;
      I_ERI_Ix3y2z_S_Px_S += I_ERI_Ix3y2z_S_Px_S_vrr;
      I_ERI_Ix2y3z_S_Px_S += I_ERI_Ix2y3z_S_Px_S_vrr;
      I_ERI_Ixy4z_S_Px_S += I_ERI_Ixy4z_S_Px_S_vrr;
      I_ERI_Ix5z_S_Px_S += I_ERI_Ix5z_S_Px_S_vrr;
      I_ERI_I6y_S_Px_S += I_ERI_I6y_S_Px_S_vrr;
      I_ERI_I5yz_S_Px_S += I_ERI_I5yz_S_Px_S_vrr;
      I_ERI_I4y2z_S_Px_S += I_ERI_I4y2z_S_Px_S_vrr;
      I_ERI_I3y3z_S_Px_S += I_ERI_I3y3z_S_Px_S_vrr;
      I_ERI_I2y4z_S_Px_S += I_ERI_I2y4z_S_Px_S_vrr;
      I_ERI_Iy5z_S_Px_S += I_ERI_Iy5z_S_Px_S_vrr;
      I_ERI_I6z_S_Px_S += I_ERI_I6z_S_Px_S_vrr;
      I_ERI_I6x_S_Py_S += I_ERI_I6x_S_Py_S_vrr;
      I_ERI_I5xy_S_Py_S += I_ERI_I5xy_S_Py_S_vrr;
      I_ERI_I5xz_S_Py_S += I_ERI_I5xz_S_Py_S_vrr;
      I_ERI_I4x2y_S_Py_S += I_ERI_I4x2y_S_Py_S_vrr;
      I_ERI_I4xyz_S_Py_S += I_ERI_I4xyz_S_Py_S_vrr;
      I_ERI_I4x2z_S_Py_S += I_ERI_I4x2z_S_Py_S_vrr;
      I_ERI_I3x3y_S_Py_S += I_ERI_I3x3y_S_Py_S_vrr;
      I_ERI_I3x2yz_S_Py_S += I_ERI_I3x2yz_S_Py_S_vrr;
      I_ERI_I3xy2z_S_Py_S += I_ERI_I3xy2z_S_Py_S_vrr;
      I_ERI_I3x3z_S_Py_S += I_ERI_I3x3z_S_Py_S_vrr;
      I_ERI_I2x4y_S_Py_S += I_ERI_I2x4y_S_Py_S_vrr;
      I_ERI_I2x3yz_S_Py_S += I_ERI_I2x3yz_S_Py_S_vrr;
      I_ERI_I2x2y2z_S_Py_S += I_ERI_I2x2y2z_S_Py_S_vrr;
      I_ERI_I2xy3z_S_Py_S += I_ERI_I2xy3z_S_Py_S_vrr;
      I_ERI_I2x4z_S_Py_S += I_ERI_I2x4z_S_Py_S_vrr;
      I_ERI_Ix5y_S_Py_S += I_ERI_Ix5y_S_Py_S_vrr;
      I_ERI_Ix4yz_S_Py_S += I_ERI_Ix4yz_S_Py_S_vrr;
      I_ERI_Ix3y2z_S_Py_S += I_ERI_Ix3y2z_S_Py_S_vrr;
      I_ERI_Ix2y3z_S_Py_S += I_ERI_Ix2y3z_S_Py_S_vrr;
      I_ERI_Ixy4z_S_Py_S += I_ERI_Ixy4z_S_Py_S_vrr;
      I_ERI_Ix5z_S_Py_S += I_ERI_Ix5z_S_Py_S_vrr;
      I_ERI_I6y_S_Py_S += I_ERI_I6y_S_Py_S_vrr;
      I_ERI_I5yz_S_Py_S += I_ERI_I5yz_S_Py_S_vrr;
      I_ERI_I4y2z_S_Py_S += I_ERI_I4y2z_S_Py_S_vrr;
      I_ERI_I3y3z_S_Py_S += I_ERI_I3y3z_S_Py_S_vrr;
      I_ERI_I2y4z_S_Py_S += I_ERI_I2y4z_S_Py_S_vrr;
      I_ERI_Iy5z_S_Py_S += I_ERI_Iy5z_S_Py_S_vrr;
      I_ERI_I6z_S_Py_S += I_ERI_I6z_S_Py_S_vrr;
      I_ERI_I6x_S_Pz_S += I_ERI_I6x_S_Pz_S_vrr;
      I_ERI_I5xy_S_Pz_S += I_ERI_I5xy_S_Pz_S_vrr;
      I_ERI_I5xz_S_Pz_S += I_ERI_I5xz_S_Pz_S_vrr;
      I_ERI_I4x2y_S_Pz_S += I_ERI_I4x2y_S_Pz_S_vrr;
      I_ERI_I4xyz_S_Pz_S += I_ERI_I4xyz_S_Pz_S_vrr;
      I_ERI_I4x2z_S_Pz_S += I_ERI_I4x2z_S_Pz_S_vrr;
      I_ERI_I3x3y_S_Pz_S += I_ERI_I3x3y_S_Pz_S_vrr;
      I_ERI_I3x2yz_S_Pz_S += I_ERI_I3x2yz_S_Pz_S_vrr;
      I_ERI_I3xy2z_S_Pz_S += I_ERI_I3xy2z_S_Pz_S_vrr;
      I_ERI_I3x3z_S_Pz_S += I_ERI_I3x3z_S_Pz_S_vrr;
      I_ERI_I2x4y_S_Pz_S += I_ERI_I2x4y_S_Pz_S_vrr;
      I_ERI_I2x3yz_S_Pz_S += I_ERI_I2x3yz_S_Pz_S_vrr;
      I_ERI_I2x2y2z_S_Pz_S += I_ERI_I2x2y2z_S_Pz_S_vrr;
      I_ERI_I2xy3z_S_Pz_S += I_ERI_I2xy3z_S_Pz_S_vrr;
      I_ERI_I2x4z_S_Pz_S += I_ERI_I2x4z_S_Pz_S_vrr;
      I_ERI_Ix5y_S_Pz_S += I_ERI_Ix5y_S_Pz_S_vrr;
      I_ERI_Ix4yz_S_Pz_S += I_ERI_Ix4yz_S_Pz_S_vrr;
      I_ERI_Ix3y2z_S_Pz_S += I_ERI_Ix3y2z_S_Pz_S_vrr;
      I_ERI_Ix2y3z_S_Pz_S += I_ERI_Ix2y3z_S_Pz_S_vrr;
      I_ERI_Ixy4z_S_Pz_S += I_ERI_Ixy4z_S_Pz_S_vrr;
      I_ERI_Ix5z_S_Pz_S += I_ERI_Ix5z_S_Pz_S_vrr;
      I_ERI_I6y_S_Pz_S += I_ERI_I6y_S_Pz_S_vrr;
      I_ERI_I5yz_S_Pz_S += I_ERI_I5yz_S_Pz_S_vrr;
      I_ERI_I4y2z_S_Pz_S += I_ERI_I4y2z_S_Pz_S_vrr;
      I_ERI_I3y3z_S_Pz_S += I_ERI_I3y3z_S_Pz_S_vrr;
      I_ERI_I2y4z_S_Pz_S += I_ERI_I2y4z_S_Pz_S_vrr;
      I_ERI_Iy5z_S_Pz_S += I_ERI_Iy5z_S_Pz_S_vrr;
      I_ERI_I6z_S_Pz_S += I_ERI_I6z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_H5x_S_Px_S += I_ERI_H5x_S_Px_S_vrr;
      I_ERI_H4xy_S_Px_S += I_ERI_H4xy_S_Px_S_vrr;
      I_ERI_H4xz_S_Px_S += I_ERI_H4xz_S_Px_S_vrr;
      I_ERI_H3x2y_S_Px_S += I_ERI_H3x2y_S_Px_S_vrr;
      I_ERI_H3xyz_S_Px_S += I_ERI_H3xyz_S_Px_S_vrr;
      I_ERI_H3x2z_S_Px_S += I_ERI_H3x2z_S_Px_S_vrr;
      I_ERI_H2x3y_S_Px_S += I_ERI_H2x3y_S_Px_S_vrr;
      I_ERI_H2x2yz_S_Px_S += I_ERI_H2x2yz_S_Px_S_vrr;
      I_ERI_H2xy2z_S_Px_S += I_ERI_H2xy2z_S_Px_S_vrr;
      I_ERI_H2x3z_S_Px_S += I_ERI_H2x3z_S_Px_S_vrr;
      I_ERI_Hx4y_S_Px_S += I_ERI_Hx4y_S_Px_S_vrr;
      I_ERI_Hx3yz_S_Px_S += I_ERI_Hx3yz_S_Px_S_vrr;
      I_ERI_Hx2y2z_S_Px_S += I_ERI_Hx2y2z_S_Px_S_vrr;
      I_ERI_Hxy3z_S_Px_S += I_ERI_Hxy3z_S_Px_S_vrr;
      I_ERI_Hx4z_S_Px_S += I_ERI_Hx4z_S_Px_S_vrr;
      I_ERI_H5y_S_Px_S += I_ERI_H5y_S_Px_S_vrr;
      I_ERI_H4yz_S_Px_S += I_ERI_H4yz_S_Px_S_vrr;
      I_ERI_H3y2z_S_Px_S += I_ERI_H3y2z_S_Px_S_vrr;
      I_ERI_H2y3z_S_Px_S += I_ERI_H2y3z_S_Px_S_vrr;
      I_ERI_Hy4z_S_Px_S += I_ERI_Hy4z_S_Px_S_vrr;
      I_ERI_H5z_S_Px_S += I_ERI_H5z_S_Px_S_vrr;
      I_ERI_H5x_S_Py_S += I_ERI_H5x_S_Py_S_vrr;
      I_ERI_H4xy_S_Py_S += I_ERI_H4xy_S_Py_S_vrr;
      I_ERI_H4xz_S_Py_S += I_ERI_H4xz_S_Py_S_vrr;
      I_ERI_H3x2y_S_Py_S += I_ERI_H3x2y_S_Py_S_vrr;
      I_ERI_H3xyz_S_Py_S += I_ERI_H3xyz_S_Py_S_vrr;
      I_ERI_H3x2z_S_Py_S += I_ERI_H3x2z_S_Py_S_vrr;
      I_ERI_H2x3y_S_Py_S += I_ERI_H2x3y_S_Py_S_vrr;
      I_ERI_H2x2yz_S_Py_S += I_ERI_H2x2yz_S_Py_S_vrr;
      I_ERI_H2xy2z_S_Py_S += I_ERI_H2xy2z_S_Py_S_vrr;
      I_ERI_H2x3z_S_Py_S += I_ERI_H2x3z_S_Py_S_vrr;
      I_ERI_Hx4y_S_Py_S += I_ERI_Hx4y_S_Py_S_vrr;
      I_ERI_Hx3yz_S_Py_S += I_ERI_Hx3yz_S_Py_S_vrr;
      I_ERI_Hx2y2z_S_Py_S += I_ERI_Hx2y2z_S_Py_S_vrr;
      I_ERI_Hxy3z_S_Py_S += I_ERI_Hxy3z_S_Py_S_vrr;
      I_ERI_Hx4z_S_Py_S += I_ERI_Hx4z_S_Py_S_vrr;
      I_ERI_H5y_S_Py_S += I_ERI_H5y_S_Py_S_vrr;
      I_ERI_H4yz_S_Py_S += I_ERI_H4yz_S_Py_S_vrr;
      I_ERI_H3y2z_S_Py_S += I_ERI_H3y2z_S_Py_S_vrr;
      I_ERI_H2y3z_S_Py_S += I_ERI_H2y3z_S_Py_S_vrr;
      I_ERI_Hy4z_S_Py_S += I_ERI_Hy4z_S_Py_S_vrr;
      I_ERI_H5z_S_Py_S += I_ERI_H5z_S_Py_S_vrr;
      I_ERI_H5x_S_Pz_S += I_ERI_H5x_S_Pz_S_vrr;
      I_ERI_H4xy_S_Pz_S += I_ERI_H4xy_S_Pz_S_vrr;
      I_ERI_H4xz_S_Pz_S += I_ERI_H4xz_S_Pz_S_vrr;
      I_ERI_H3x2y_S_Pz_S += I_ERI_H3x2y_S_Pz_S_vrr;
      I_ERI_H3xyz_S_Pz_S += I_ERI_H3xyz_S_Pz_S_vrr;
      I_ERI_H3x2z_S_Pz_S += I_ERI_H3x2z_S_Pz_S_vrr;
      I_ERI_H2x3y_S_Pz_S += I_ERI_H2x3y_S_Pz_S_vrr;
      I_ERI_H2x2yz_S_Pz_S += I_ERI_H2x2yz_S_Pz_S_vrr;
      I_ERI_H2xy2z_S_Pz_S += I_ERI_H2xy2z_S_Pz_S_vrr;
      I_ERI_H2x3z_S_Pz_S += I_ERI_H2x3z_S_Pz_S_vrr;
      I_ERI_Hx4y_S_Pz_S += I_ERI_Hx4y_S_Pz_S_vrr;
      I_ERI_Hx3yz_S_Pz_S += I_ERI_Hx3yz_S_Pz_S_vrr;
      I_ERI_Hx2y2z_S_Pz_S += I_ERI_Hx2y2z_S_Pz_S_vrr;
      I_ERI_Hxy3z_S_Pz_S += I_ERI_Hxy3z_S_Pz_S_vrr;
      I_ERI_Hx4z_S_Pz_S += I_ERI_Hx4z_S_Pz_S_vrr;
      I_ERI_H5y_S_Pz_S += I_ERI_H5y_S_Pz_S_vrr;
      I_ERI_H4yz_S_Pz_S += I_ERI_H4yz_S_Pz_S_vrr;
      I_ERI_H3y2z_S_Pz_S += I_ERI_H3y2z_S_Pz_S_vrr;
      I_ERI_H2y3z_S_Pz_S += I_ERI_H2y3z_S_Pz_S_vrr;
      I_ERI_Hy4z_S_Pz_S += I_ERI_Hy4z_S_Pz_S_vrr;
      I_ERI_H5z_S_Pz_S += I_ERI_H5z_S_Pz_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_G4x_S_Px_S += I_ERI_G4x_S_Px_S_vrr;
      I_ERI_G3xy_S_Px_S += I_ERI_G3xy_S_Px_S_vrr;
      I_ERI_G3xz_S_Px_S += I_ERI_G3xz_S_Px_S_vrr;
      I_ERI_G2x2y_S_Px_S += I_ERI_G2x2y_S_Px_S_vrr;
      I_ERI_G2xyz_S_Px_S += I_ERI_G2xyz_S_Px_S_vrr;
      I_ERI_G2x2z_S_Px_S += I_ERI_G2x2z_S_Px_S_vrr;
      I_ERI_Gx3y_S_Px_S += I_ERI_Gx3y_S_Px_S_vrr;
      I_ERI_Gx2yz_S_Px_S += I_ERI_Gx2yz_S_Px_S_vrr;
      I_ERI_Gxy2z_S_Px_S += I_ERI_Gxy2z_S_Px_S_vrr;
      I_ERI_Gx3z_S_Px_S += I_ERI_Gx3z_S_Px_S_vrr;
      I_ERI_G4y_S_Px_S += I_ERI_G4y_S_Px_S_vrr;
      I_ERI_G3yz_S_Px_S += I_ERI_G3yz_S_Px_S_vrr;
      I_ERI_G2y2z_S_Px_S += I_ERI_G2y2z_S_Px_S_vrr;
      I_ERI_Gy3z_S_Px_S += I_ERI_Gy3z_S_Px_S_vrr;
      I_ERI_G4z_S_Px_S += I_ERI_G4z_S_Px_S_vrr;
      I_ERI_G4x_S_Py_S += I_ERI_G4x_S_Py_S_vrr;
      I_ERI_G3xy_S_Py_S += I_ERI_G3xy_S_Py_S_vrr;
      I_ERI_G3xz_S_Py_S += I_ERI_G3xz_S_Py_S_vrr;
      I_ERI_G2x2y_S_Py_S += I_ERI_G2x2y_S_Py_S_vrr;
      I_ERI_G2xyz_S_Py_S += I_ERI_G2xyz_S_Py_S_vrr;
      I_ERI_G2x2z_S_Py_S += I_ERI_G2x2z_S_Py_S_vrr;
      I_ERI_Gx3y_S_Py_S += I_ERI_Gx3y_S_Py_S_vrr;
      I_ERI_Gx2yz_S_Py_S += I_ERI_Gx2yz_S_Py_S_vrr;
      I_ERI_Gxy2z_S_Py_S += I_ERI_Gxy2z_S_Py_S_vrr;
      I_ERI_Gx3z_S_Py_S += I_ERI_Gx3z_S_Py_S_vrr;
      I_ERI_G4y_S_Py_S += I_ERI_G4y_S_Py_S_vrr;
      I_ERI_G3yz_S_Py_S += I_ERI_G3yz_S_Py_S_vrr;
      I_ERI_G2y2z_S_Py_S += I_ERI_G2y2z_S_Py_S_vrr;
      I_ERI_Gy3z_S_Py_S += I_ERI_Gy3z_S_Py_S_vrr;
      I_ERI_G4z_S_Py_S += I_ERI_G4z_S_Py_S_vrr;
      I_ERI_G4x_S_Pz_S += I_ERI_G4x_S_Pz_S_vrr;
      I_ERI_G3xy_S_Pz_S += I_ERI_G3xy_S_Pz_S_vrr;
      I_ERI_G3xz_S_Pz_S += I_ERI_G3xz_S_Pz_S_vrr;
      I_ERI_G2x2y_S_Pz_S += I_ERI_G2x2y_S_Pz_S_vrr;
      I_ERI_G2xyz_S_Pz_S += I_ERI_G2xyz_S_Pz_S_vrr;
      I_ERI_G2x2z_S_Pz_S += I_ERI_G2x2z_S_Pz_S_vrr;
      I_ERI_Gx3y_S_Pz_S += I_ERI_Gx3y_S_Pz_S_vrr;
      I_ERI_Gx2yz_S_Pz_S += I_ERI_Gx2yz_S_Pz_S_vrr;
      I_ERI_Gxy2z_S_Pz_S += I_ERI_Gxy2z_S_Pz_S_vrr;
      I_ERI_Gx3z_S_Pz_S += I_ERI_Gx3z_S_Pz_S_vrr;
      I_ERI_G4y_S_Pz_S += I_ERI_G4y_S_Pz_S_vrr;
      I_ERI_G3yz_S_Pz_S += I_ERI_G3yz_S_Pz_S_vrr;
      I_ERI_G2y2z_S_Pz_S += I_ERI_G2y2z_S_Pz_S_vrr;
      I_ERI_Gy3z_S_Pz_S += I_ERI_Gy3z_S_Pz_S_vrr;
      I_ERI_G4z_S_Pz_S += I_ERI_G4z_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_G_P_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_P_S
   * RHS shell quartet name: SQ_ERI_G_S_P_S
   ************************************************************/
  Double I_ERI_G4x_Px_Px_S = I_ERI_H5x_S_Px_S+ABX*I_ERI_G4x_S_Px_S;
  Double I_ERI_G3xy_Px_Px_S = I_ERI_H4xy_S_Px_S+ABX*I_ERI_G3xy_S_Px_S;
  Double I_ERI_G3xz_Px_Px_S = I_ERI_H4xz_S_Px_S+ABX*I_ERI_G3xz_S_Px_S;
  Double I_ERI_G2x2y_Px_Px_S = I_ERI_H3x2y_S_Px_S+ABX*I_ERI_G2x2y_S_Px_S;
  Double I_ERI_G2xyz_Px_Px_S = I_ERI_H3xyz_S_Px_S+ABX*I_ERI_G2xyz_S_Px_S;
  Double I_ERI_G2x2z_Px_Px_S = I_ERI_H3x2z_S_Px_S+ABX*I_ERI_G2x2z_S_Px_S;
  Double I_ERI_Gx3y_Px_Px_S = I_ERI_H2x3y_S_Px_S+ABX*I_ERI_Gx3y_S_Px_S;
  Double I_ERI_Gx2yz_Px_Px_S = I_ERI_H2x2yz_S_Px_S+ABX*I_ERI_Gx2yz_S_Px_S;
  Double I_ERI_Gxy2z_Px_Px_S = I_ERI_H2xy2z_S_Px_S+ABX*I_ERI_Gxy2z_S_Px_S;
  Double I_ERI_Gx3z_Px_Px_S = I_ERI_H2x3z_S_Px_S+ABX*I_ERI_Gx3z_S_Px_S;
  Double I_ERI_G4y_Px_Px_S = I_ERI_Hx4y_S_Px_S+ABX*I_ERI_G4y_S_Px_S;
  Double I_ERI_G3yz_Px_Px_S = I_ERI_Hx3yz_S_Px_S+ABX*I_ERI_G3yz_S_Px_S;
  Double I_ERI_G2y2z_Px_Px_S = I_ERI_Hx2y2z_S_Px_S+ABX*I_ERI_G2y2z_S_Px_S;
  Double I_ERI_Gy3z_Px_Px_S = I_ERI_Hxy3z_S_Px_S+ABX*I_ERI_Gy3z_S_Px_S;
  Double I_ERI_G4z_Px_Px_S = I_ERI_Hx4z_S_Px_S+ABX*I_ERI_G4z_S_Px_S;
  Double I_ERI_G4x_Py_Px_S = I_ERI_H4xy_S_Px_S+ABY*I_ERI_G4x_S_Px_S;
  Double I_ERI_G3xy_Py_Px_S = I_ERI_H3x2y_S_Px_S+ABY*I_ERI_G3xy_S_Px_S;
  Double I_ERI_G3xz_Py_Px_S = I_ERI_H3xyz_S_Px_S+ABY*I_ERI_G3xz_S_Px_S;
  Double I_ERI_G2x2y_Py_Px_S = I_ERI_H2x3y_S_Px_S+ABY*I_ERI_G2x2y_S_Px_S;
  Double I_ERI_G2xyz_Py_Px_S = I_ERI_H2x2yz_S_Px_S+ABY*I_ERI_G2xyz_S_Px_S;
  Double I_ERI_G2x2z_Py_Px_S = I_ERI_H2xy2z_S_Px_S+ABY*I_ERI_G2x2z_S_Px_S;
  Double I_ERI_Gx3y_Py_Px_S = I_ERI_Hx4y_S_Px_S+ABY*I_ERI_Gx3y_S_Px_S;
  Double I_ERI_Gx2yz_Py_Px_S = I_ERI_Hx3yz_S_Px_S+ABY*I_ERI_Gx2yz_S_Px_S;
  Double I_ERI_Gxy2z_Py_Px_S = I_ERI_Hx2y2z_S_Px_S+ABY*I_ERI_Gxy2z_S_Px_S;
  Double I_ERI_Gx3z_Py_Px_S = I_ERI_Hxy3z_S_Px_S+ABY*I_ERI_Gx3z_S_Px_S;
  Double I_ERI_G4y_Py_Px_S = I_ERI_H5y_S_Px_S+ABY*I_ERI_G4y_S_Px_S;
  Double I_ERI_G3yz_Py_Px_S = I_ERI_H4yz_S_Px_S+ABY*I_ERI_G3yz_S_Px_S;
  Double I_ERI_G2y2z_Py_Px_S = I_ERI_H3y2z_S_Px_S+ABY*I_ERI_G2y2z_S_Px_S;
  Double I_ERI_Gy3z_Py_Px_S = I_ERI_H2y3z_S_Px_S+ABY*I_ERI_Gy3z_S_Px_S;
  Double I_ERI_G4z_Py_Px_S = I_ERI_Hy4z_S_Px_S+ABY*I_ERI_G4z_S_Px_S;
  Double I_ERI_G4x_Pz_Px_S = I_ERI_H4xz_S_Px_S+ABZ*I_ERI_G4x_S_Px_S;
  Double I_ERI_G3xy_Pz_Px_S = I_ERI_H3xyz_S_Px_S+ABZ*I_ERI_G3xy_S_Px_S;
  Double I_ERI_G3xz_Pz_Px_S = I_ERI_H3x2z_S_Px_S+ABZ*I_ERI_G3xz_S_Px_S;
  Double I_ERI_G2x2y_Pz_Px_S = I_ERI_H2x2yz_S_Px_S+ABZ*I_ERI_G2x2y_S_Px_S;
  Double I_ERI_G2xyz_Pz_Px_S = I_ERI_H2xy2z_S_Px_S+ABZ*I_ERI_G2xyz_S_Px_S;
  Double I_ERI_G2x2z_Pz_Px_S = I_ERI_H2x3z_S_Px_S+ABZ*I_ERI_G2x2z_S_Px_S;
  Double I_ERI_Gx3y_Pz_Px_S = I_ERI_Hx3yz_S_Px_S+ABZ*I_ERI_Gx3y_S_Px_S;
  Double I_ERI_Gx2yz_Pz_Px_S = I_ERI_Hx2y2z_S_Px_S+ABZ*I_ERI_Gx2yz_S_Px_S;
  Double I_ERI_Gxy2z_Pz_Px_S = I_ERI_Hxy3z_S_Px_S+ABZ*I_ERI_Gxy2z_S_Px_S;
  Double I_ERI_Gx3z_Pz_Px_S = I_ERI_Hx4z_S_Px_S+ABZ*I_ERI_Gx3z_S_Px_S;
  Double I_ERI_G4y_Pz_Px_S = I_ERI_H4yz_S_Px_S+ABZ*I_ERI_G4y_S_Px_S;
  Double I_ERI_G3yz_Pz_Px_S = I_ERI_H3y2z_S_Px_S+ABZ*I_ERI_G3yz_S_Px_S;
  Double I_ERI_G2y2z_Pz_Px_S = I_ERI_H2y3z_S_Px_S+ABZ*I_ERI_G2y2z_S_Px_S;
  Double I_ERI_Gy3z_Pz_Px_S = I_ERI_Hy4z_S_Px_S+ABZ*I_ERI_Gy3z_S_Px_S;
  Double I_ERI_G4z_Pz_Px_S = I_ERI_H5z_S_Px_S+ABZ*I_ERI_G4z_S_Px_S;
  Double I_ERI_G4x_Px_Py_S = I_ERI_H5x_S_Py_S+ABX*I_ERI_G4x_S_Py_S;
  Double I_ERI_G3xy_Px_Py_S = I_ERI_H4xy_S_Py_S+ABX*I_ERI_G3xy_S_Py_S;
  Double I_ERI_G3xz_Px_Py_S = I_ERI_H4xz_S_Py_S+ABX*I_ERI_G3xz_S_Py_S;
  Double I_ERI_G2x2y_Px_Py_S = I_ERI_H3x2y_S_Py_S+ABX*I_ERI_G2x2y_S_Py_S;
  Double I_ERI_G2xyz_Px_Py_S = I_ERI_H3xyz_S_Py_S+ABX*I_ERI_G2xyz_S_Py_S;
  Double I_ERI_G2x2z_Px_Py_S = I_ERI_H3x2z_S_Py_S+ABX*I_ERI_G2x2z_S_Py_S;
  Double I_ERI_Gx3y_Px_Py_S = I_ERI_H2x3y_S_Py_S+ABX*I_ERI_Gx3y_S_Py_S;
  Double I_ERI_Gx2yz_Px_Py_S = I_ERI_H2x2yz_S_Py_S+ABX*I_ERI_Gx2yz_S_Py_S;
  Double I_ERI_Gxy2z_Px_Py_S = I_ERI_H2xy2z_S_Py_S+ABX*I_ERI_Gxy2z_S_Py_S;
  Double I_ERI_Gx3z_Px_Py_S = I_ERI_H2x3z_S_Py_S+ABX*I_ERI_Gx3z_S_Py_S;
  Double I_ERI_G4y_Px_Py_S = I_ERI_Hx4y_S_Py_S+ABX*I_ERI_G4y_S_Py_S;
  Double I_ERI_G3yz_Px_Py_S = I_ERI_Hx3yz_S_Py_S+ABX*I_ERI_G3yz_S_Py_S;
  Double I_ERI_G2y2z_Px_Py_S = I_ERI_Hx2y2z_S_Py_S+ABX*I_ERI_G2y2z_S_Py_S;
  Double I_ERI_Gy3z_Px_Py_S = I_ERI_Hxy3z_S_Py_S+ABX*I_ERI_Gy3z_S_Py_S;
  Double I_ERI_G4z_Px_Py_S = I_ERI_Hx4z_S_Py_S+ABX*I_ERI_G4z_S_Py_S;
  Double I_ERI_G4x_Py_Py_S = I_ERI_H4xy_S_Py_S+ABY*I_ERI_G4x_S_Py_S;
  Double I_ERI_G3xy_Py_Py_S = I_ERI_H3x2y_S_Py_S+ABY*I_ERI_G3xy_S_Py_S;
  Double I_ERI_G3xz_Py_Py_S = I_ERI_H3xyz_S_Py_S+ABY*I_ERI_G3xz_S_Py_S;
  Double I_ERI_G2x2y_Py_Py_S = I_ERI_H2x3y_S_Py_S+ABY*I_ERI_G2x2y_S_Py_S;
  Double I_ERI_G2xyz_Py_Py_S = I_ERI_H2x2yz_S_Py_S+ABY*I_ERI_G2xyz_S_Py_S;
  Double I_ERI_G2x2z_Py_Py_S = I_ERI_H2xy2z_S_Py_S+ABY*I_ERI_G2x2z_S_Py_S;
  Double I_ERI_Gx3y_Py_Py_S = I_ERI_Hx4y_S_Py_S+ABY*I_ERI_Gx3y_S_Py_S;
  Double I_ERI_Gx2yz_Py_Py_S = I_ERI_Hx3yz_S_Py_S+ABY*I_ERI_Gx2yz_S_Py_S;
  Double I_ERI_Gxy2z_Py_Py_S = I_ERI_Hx2y2z_S_Py_S+ABY*I_ERI_Gxy2z_S_Py_S;
  Double I_ERI_Gx3z_Py_Py_S = I_ERI_Hxy3z_S_Py_S+ABY*I_ERI_Gx3z_S_Py_S;
  Double I_ERI_G4y_Py_Py_S = I_ERI_H5y_S_Py_S+ABY*I_ERI_G4y_S_Py_S;
  Double I_ERI_G3yz_Py_Py_S = I_ERI_H4yz_S_Py_S+ABY*I_ERI_G3yz_S_Py_S;
  Double I_ERI_G2y2z_Py_Py_S = I_ERI_H3y2z_S_Py_S+ABY*I_ERI_G2y2z_S_Py_S;
  Double I_ERI_Gy3z_Py_Py_S = I_ERI_H2y3z_S_Py_S+ABY*I_ERI_Gy3z_S_Py_S;
  Double I_ERI_G4z_Py_Py_S = I_ERI_Hy4z_S_Py_S+ABY*I_ERI_G4z_S_Py_S;
  Double I_ERI_G4x_Pz_Py_S = I_ERI_H4xz_S_Py_S+ABZ*I_ERI_G4x_S_Py_S;
  Double I_ERI_G3xy_Pz_Py_S = I_ERI_H3xyz_S_Py_S+ABZ*I_ERI_G3xy_S_Py_S;
  Double I_ERI_G3xz_Pz_Py_S = I_ERI_H3x2z_S_Py_S+ABZ*I_ERI_G3xz_S_Py_S;
  Double I_ERI_G2x2y_Pz_Py_S = I_ERI_H2x2yz_S_Py_S+ABZ*I_ERI_G2x2y_S_Py_S;
  Double I_ERI_G2xyz_Pz_Py_S = I_ERI_H2xy2z_S_Py_S+ABZ*I_ERI_G2xyz_S_Py_S;
  Double I_ERI_G2x2z_Pz_Py_S = I_ERI_H2x3z_S_Py_S+ABZ*I_ERI_G2x2z_S_Py_S;
  Double I_ERI_Gx3y_Pz_Py_S = I_ERI_Hx3yz_S_Py_S+ABZ*I_ERI_Gx3y_S_Py_S;
  Double I_ERI_Gx2yz_Pz_Py_S = I_ERI_Hx2y2z_S_Py_S+ABZ*I_ERI_Gx2yz_S_Py_S;
  Double I_ERI_Gxy2z_Pz_Py_S = I_ERI_Hxy3z_S_Py_S+ABZ*I_ERI_Gxy2z_S_Py_S;
  Double I_ERI_Gx3z_Pz_Py_S = I_ERI_Hx4z_S_Py_S+ABZ*I_ERI_Gx3z_S_Py_S;
  Double I_ERI_G4y_Pz_Py_S = I_ERI_H4yz_S_Py_S+ABZ*I_ERI_G4y_S_Py_S;
  Double I_ERI_G3yz_Pz_Py_S = I_ERI_H3y2z_S_Py_S+ABZ*I_ERI_G3yz_S_Py_S;
  Double I_ERI_G2y2z_Pz_Py_S = I_ERI_H2y3z_S_Py_S+ABZ*I_ERI_G2y2z_S_Py_S;
  Double I_ERI_Gy3z_Pz_Py_S = I_ERI_Hy4z_S_Py_S+ABZ*I_ERI_Gy3z_S_Py_S;
  Double I_ERI_G4z_Pz_Py_S = I_ERI_H5z_S_Py_S+ABZ*I_ERI_G4z_S_Py_S;
  Double I_ERI_G4x_Px_Pz_S = I_ERI_H5x_S_Pz_S+ABX*I_ERI_G4x_S_Pz_S;
  Double I_ERI_G3xy_Px_Pz_S = I_ERI_H4xy_S_Pz_S+ABX*I_ERI_G3xy_S_Pz_S;
  Double I_ERI_G3xz_Px_Pz_S = I_ERI_H4xz_S_Pz_S+ABX*I_ERI_G3xz_S_Pz_S;
  Double I_ERI_G2x2y_Px_Pz_S = I_ERI_H3x2y_S_Pz_S+ABX*I_ERI_G2x2y_S_Pz_S;
  Double I_ERI_G2xyz_Px_Pz_S = I_ERI_H3xyz_S_Pz_S+ABX*I_ERI_G2xyz_S_Pz_S;
  Double I_ERI_G2x2z_Px_Pz_S = I_ERI_H3x2z_S_Pz_S+ABX*I_ERI_G2x2z_S_Pz_S;
  Double I_ERI_Gx3y_Px_Pz_S = I_ERI_H2x3y_S_Pz_S+ABX*I_ERI_Gx3y_S_Pz_S;
  Double I_ERI_Gx2yz_Px_Pz_S = I_ERI_H2x2yz_S_Pz_S+ABX*I_ERI_Gx2yz_S_Pz_S;
  Double I_ERI_Gxy2z_Px_Pz_S = I_ERI_H2xy2z_S_Pz_S+ABX*I_ERI_Gxy2z_S_Pz_S;
  Double I_ERI_Gx3z_Px_Pz_S = I_ERI_H2x3z_S_Pz_S+ABX*I_ERI_Gx3z_S_Pz_S;
  Double I_ERI_G4y_Px_Pz_S = I_ERI_Hx4y_S_Pz_S+ABX*I_ERI_G4y_S_Pz_S;
  Double I_ERI_G3yz_Px_Pz_S = I_ERI_Hx3yz_S_Pz_S+ABX*I_ERI_G3yz_S_Pz_S;
  Double I_ERI_G2y2z_Px_Pz_S = I_ERI_Hx2y2z_S_Pz_S+ABX*I_ERI_G2y2z_S_Pz_S;
  Double I_ERI_Gy3z_Px_Pz_S = I_ERI_Hxy3z_S_Pz_S+ABX*I_ERI_Gy3z_S_Pz_S;
  Double I_ERI_G4z_Px_Pz_S = I_ERI_Hx4z_S_Pz_S+ABX*I_ERI_G4z_S_Pz_S;
  Double I_ERI_G4x_Py_Pz_S = I_ERI_H4xy_S_Pz_S+ABY*I_ERI_G4x_S_Pz_S;
  Double I_ERI_G3xy_Py_Pz_S = I_ERI_H3x2y_S_Pz_S+ABY*I_ERI_G3xy_S_Pz_S;
  Double I_ERI_G3xz_Py_Pz_S = I_ERI_H3xyz_S_Pz_S+ABY*I_ERI_G3xz_S_Pz_S;
  Double I_ERI_G2x2y_Py_Pz_S = I_ERI_H2x3y_S_Pz_S+ABY*I_ERI_G2x2y_S_Pz_S;
  Double I_ERI_G2xyz_Py_Pz_S = I_ERI_H2x2yz_S_Pz_S+ABY*I_ERI_G2xyz_S_Pz_S;
  Double I_ERI_G2x2z_Py_Pz_S = I_ERI_H2xy2z_S_Pz_S+ABY*I_ERI_G2x2z_S_Pz_S;
  Double I_ERI_Gx3y_Py_Pz_S = I_ERI_Hx4y_S_Pz_S+ABY*I_ERI_Gx3y_S_Pz_S;
  Double I_ERI_Gx2yz_Py_Pz_S = I_ERI_Hx3yz_S_Pz_S+ABY*I_ERI_Gx2yz_S_Pz_S;
  Double I_ERI_Gxy2z_Py_Pz_S = I_ERI_Hx2y2z_S_Pz_S+ABY*I_ERI_Gxy2z_S_Pz_S;
  Double I_ERI_Gx3z_Py_Pz_S = I_ERI_Hxy3z_S_Pz_S+ABY*I_ERI_Gx3z_S_Pz_S;
  Double I_ERI_G4y_Py_Pz_S = I_ERI_H5y_S_Pz_S+ABY*I_ERI_G4y_S_Pz_S;
  Double I_ERI_G3yz_Py_Pz_S = I_ERI_H4yz_S_Pz_S+ABY*I_ERI_G3yz_S_Pz_S;
  Double I_ERI_G2y2z_Py_Pz_S = I_ERI_H3y2z_S_Pz_S+ABY*I_ERI_G2y2z_S_Pz_S;
  Double I_ERI_Gy3z_Py_Pz_S = I_ERI_H2y3z_S_Pz_S+ABY*I_ERI_Gy3z_S_Pz_S;
  Double I_ERI_G4z_Py_Pz_S = I_ERI_Hy4z_S_Pz_S+ABY*I_ERI_G4z_S_Pz_S;
  Double I_ERI_G4x_Pz_Pz_S = I_ERI_H4xz_S_Pz_S+ABZ*I_ERI_G4x_S_Pz_S;
  Double I_ERI_G3xy_Pz_Pz_S = I_ERI_H3xyz_S_Pz_S+ABZ*I_ERI_G3xy_S_Pz_S;
  Double I_ERI_G3xz_Pz_Pz_S = I_ERI_H3x2z_S_Pz_S+ABZ*I_ERI_G3xz_S_Pz_S;
  Double I_ERI_G2x2y_Pz_Pz_S = I_ERI_H2x2yz_S_Pz_S+ABZ*I_ERI_G2x2y_S_Pz_S;
  Double I_ERI_G2xyz_Pz_Pz_S = I_ERI_H2xy2z_S_Pz_S+ABZ*I_ERI_G2xyz_S_Pz_S;
  Double I_ERI_G2x2z_Pz_Pz_S = I_ERI_H2x3z_S_Pz_S+ABZ*I_ERI_G2x2z_S_Pz_S;
  Double I_ERI_Gx3y_Pz_Pz_S = I_ERI_Hx3yz_S_Pz_S+ABZ*I_ERI_Gx3y_S_Pz_S;
  Double I_ERI_Gx2yz_Pz_Pz_S = I_ERI_Hx2y2z_S_Pz_S+ABZ*I_ERI_Gx2yz_S_Pz_S;
  Double I_ERI_Gxy2z_Pz_Pz_S = I_ERI_Hxy3z_S_Pz_S+ABZ*I_ERI_Gxy2z_S_Pz_S;
  Double I_ERI_Gx3z_Pz_Pz_S = I_ERI_Hx4z_S_Pz_S+ABZ*I_ERI_Gx3z_S_Pz_S;
  Double I_ERI_G4y_Pz_Pz_S = I_ERI_H4yz_S_Pz_S+ABZ*I_ERI_G4y_S_Pz_S;
  Double I_ERI_G3yz_Pz_Pz_S = I_ERI_H3y2z_S_Pz_S+ABZ*I_ERI_G3yz_S_Pz_S;
  Double I_ERI_G2y2z_Pz_Pz_S = I_ERI_H2y3z_S_Pz_S+ABZ*I_ERI_G2y2z_S_Pz_S;
  Double I_ERI_Gy3z_Pz_Pz_S = I_ERI_Hy4z_S_Pz_S+ABZ*I_ERI_Gy3z_S_Pz_S;
  Double I_ERI_G4z_Pz_Pz_S = I_ERI_H5z_S_Pz_S+ABZ*I_ERI_G4z_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_H_P_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_S_P_S
   * RHS shell quartet name: SQ_ERI_H_S_P_S
   ************************************************************/
  Double I_ERI_H5x_Px_Px_S = I_ERI_I6x_S_Px_S+ABX*I_ERI_H5x_S_Px_S;
  Double I_ERI_H4xy_Px_Px_S = I_ERI_I5xy_S_Px_S+ABX*I_ERI_H4xy_S_Px_S;
  Double I_ERI_H4xz_Px_Px_S = I_ERI_I5xz_S_Px_S+ABX*I_ERI_H4xz_S_Px_S;
  Double I_ERI_H3x2y_Px_Px_S = I_ERI_I4x2y_S_Px_S+ABX*I_ERI_H3x2y_S_Px_S;
  Double I_ERI_H3xyz_Px_Px_S = I_ERI_I4xyz_S_Px_S+ABX*I_ERI_H3xyz_S_Px_S;
  Double I_ERI_H3x2z_Px_Px_S = I_ERI_I4x2z_S_Px_S+ABX*I_ERI_H3x2z_S_Px_S;
  Double I_ERI_H2x3y_Px_Px_S = I_ERI_I3x3y_S_Px_S+ABX*I_ERI_H2x3y_S_Px_S;
  Double I_ERI_H2x2yz_Px_Px_S = I_ERI_I3x2yz_S_Px_S+ABX*I_ERI_H2x2yz_S_Px_S;
  Double I_ERI_H2xy2z_Px_Px_S = I_ERI_I3xy2z_S_Px_S+ABX*I_ERI_H2xy2z_S_Px_S;
  Double I_ERI_H2x3z_Px_Px_S = I_ERI_I3x3z_S_Px_S+ABX*I_ERI_H2x3z_S_Px_S;
  Double I_ERI_Hx4y_Px_Px_S = I_ERI_I2x4y_S_Px_S+ABX*I_ERI_Hx4y_S_Px_S;
  Double I_ERI_Hx3yz_Px_Px_S = I_ERI_I2x3yz_S_Px_S+ABX*I_ERI_Hx3yz_S_Px_S;
  Double I_ERI_Hx2y2z_Px_Px_S = I_ERI_I2x2y2z_S_Px_S+ABX*I_ERI_Hx2y2z_S_Px_S;
  Double I_ERI_Hxy3z_Px_Px_S = I_ERI_I2xy3z_S_Px_S+ABX*I_ERI_Hxy3z_S_Px_S;
  Double I_ERI_Hx4z_Px_Px_S = I_ERI_I2x4z_S_Px_S+ABX*I_ERI_Hx4z_S_Px_S;
  Double I_ERI_H5y_Px_Px_S = I_ERI_Ix5y_S_Px_S+ABX*I_ERI_H5y_S_Px_S;
  Double I_ERI_H4yz_Px_Px_S = I_ERI_Ix4yz_S_Px_S+ABX*I_ERI_H4yz_S_Px_S;
  Double I_ERI_H3y2z_Px_Px_S = I_ERI_Ix3y2z_S_Px_S+ABX*I_ERI_H3y2z_S_Px_S;
  Double I_ERI_H2y3z_Px_Px_S = I_ERI_Ix2y3z_S_Px_S+ABX*I_ERI_H2y3z_S_Px_S;
  Double I_ERI_Hy4z_Px_Px_S = I_ERI_Ixy4z_S_Px_S+ABX*I_ERI_Hy4z_S_Px_S;
  Double I_ERI_H5z_Px_Px_S = I_ERI_Ix5z_S_Px_S+ABX*I_ERI_H5z_S_Px_S;
  Double I_ERI_H5x_Py_Px_S = I_ERI_I5xy_S_Px_S+ABY*I_ERI_H5x_S_Px_S;
  Double I_ERI_H4xy_Py_Px_S = I_ERI_I4x2y_S_Px_S+ABY*I_ERI_H4xy_S_Px_S;
  Double I_ERI_H4xz_Py_Px_S = I_ERI_I4xyz_S_Px_S+ABY*I_ERI_H4xz_S_Px_S;
  Double I_ERI_H3x2y_Py_Px_S = I_ERI_I3x3y_S_Px_S+ABY*I_ERI_H3x2y_S_Px_S;
  Double I_ERI_H3xyz_Py_Px_S = I_ERI_I3x2yz_S_Px_S+ABY*I_ERI_H3xyz_S_Px_S;
  Double I_ERI_H3x2z_Py_Px_S = I_ERI_I3xy2z_S_Px_S+ABY*I_ERI_H3x2z_S_Px_S;
  Double I_ERI_H2x3y_Py_Px_S = I_ERI_I2x4y_S_Px_S+ABY*I_ERI_H2x3y_S_Px_S;
  Double I_ERI_H2x2yz_Py_Px_S = I_ERI_I2x3yz_S_Px_S+ABY*I_ERI_H2x2yz_S_Px_S;
  Double I_ERI_H2xy2z_Py_Px_S = I_ERI_I2x2y2z_S_Px_S+ABY*I_ERI_H2xy2z_S_Px_S;
  Double I_ERI_H2x3z_Py_Px_S = I_ERI_I2xy3z_S_Px_S+ABY*I_ERI_H2x3z_S_Px_S;
  Double I_ERI_Hx4y_Py_Px_S = I_ERI_Ix5y_S_Px_S+ABY*I_ERI_Hx4y_S_Px_S;
  Double I_ERI_Hx3yz_Py_Px_S = I_ERI_Ix4yz_S_Px_S+ABY*I_ERI_Hx3yz_S_Px_S;
  Double I_ERI_Hx2y2z_Py_Px_S = I_ERI_Ix3y2z_S_Px_S+ABY*I_ERI_Hx2y2z_S_Px_S;
  Double I_ERI_Hxy3z_Py_Px_S = I_ERI_Ix2y3z_S_Px_S+ABY*I_ERI_Hxy3z_S_Px_S;
  Double I_ERI_Hx4z_Py_Px_S = I_ERI_Ixy4z_S_Px_S+ABY*I_ERI_Hx4z_S_Px_S;
  Double I_ERI_H5y_Py_Px_S = I_ERI_I6y_S_Px_S+ABY*I_ERI_H5y_S_Px_S;
  Double I_ERI_H4yz_Py_Px_S = I_ERI_I5yz_S_Px_S+ABY*I_ERI_H4yz_S_Px_S;
  Double I_ERI_H3y2z_Py_Px_S = I_ERI_I4y2z_S_Px_S+ABY*I_ERI_H3y2z_S_Px_S;
  Double I_ERI_H2y3z_Py_Px_S = I_ERI_I3y3z_S_Px_S+ABY*I_ERI_H2y3z_S_Px_S;
  Double I_ERI_Hy4z_Py_Px_S = I_ERI_I2y4z_S_Px_S+ABY*I_ERI_Hy4z_S_Px_S;
  Double I_ERI_H5z_Py_Px_S = I_ERI_Iy5z_S_Px_S+ABY*I_ERI_H5z_S_Px_S;
  Double I_ERI_H5x_Pz_Px_S = I_ERI_I5xz_S_Px_S+ABZ*I_ERI_H5x_S_Px_S;
  Double I_ERI_H4xy_Pz_Px_S = I_ERI_I4xyz_S_Px_S+ABZ*I_ERI_H4xy_S_Px_S;
  Double I_ERI_H4xz_Pz_Px_S = I_ERI_I4x2z_S_Px_S+ABZ*I_ERI_H4xz_S_Px_S;
  Double I_ERI_H3x2y_Pz_Px_S = I_ERI_I3x2yz_S_Px_S+ABZ*I_ERI_H3x2y_S_Px_S;
  Double I_ERI_H3xyz_Pz_Px_S = I_ERI_I3xy2z_S_Px_S+ABZ*I_ERI_H3xyz_S_Px_S;
  Double I_ERI_H3x2z_Pz_Px_S = I_ERI_I3x3z_S_Px_S+ABZ*I_ERI_H3x2z_S_Px_S;
  Double I_ERI_H2x3y_Pz_Px_S = I_ERI_I2x3yz_S_Px_S+ABZ*I_ERI_H2x3y_S_Px_S;
  Double I_ERI_H2x2yz_Pz_Px_S = I_ERI_I2x2y2z_S_Px_S+ABZ*I_ERI_H2x2yz_S_Px_S;
  Double I_ERI_H2xy2z_Pz_Px_S = I_ERI_I2xy3z_S_Px_S+ABZ*I_ERI_H2xy2z_S_Px_S;
  Double I_ERI_H2x3z_Pz_Px_S = I_ERI_I2x4z_S_Px_S+ABZ*I_ERI_H2x3z_S_Px_S;
  Double I_ERI_Hx4y_Pz_Px_S = I_ERI_Ix4yz_S_Px_S+ABZ*I_ERI_Hx4y_S_Px_S;
  Double I_ERI_Hx3yz_Pz_Px_S = I_ERI_Ix3y2z_S_Px_S+ABZ*I_ERI_Hx3yz_S_Px_S;
  Double I_ERI_Hx2y2z_Pz_Px_S = I_ERI_Ix2y3z_S_Px_S+ABZ*I_ERI_Hx2y2z_S_Px_S;
  Double I_ERI_Hxy3z_Pz_Px_S = I_ERI_Ixy4z_S_Px_S+ABZ*I_ERI_Hxy3z_S_Px_S;
  Double I_ERI_Hx4z_Pz_Px_S = I_ERI_Ix5z_S_Px_S+ABZ*I_ERI_Hx4z_S_Px_S;
  Double I_ERI_H5y_Pz_Px_S = I_ERI_I5yz_S_Px_S+ABZ*I_ERI_H5y_S_Px_S;
  Double I_ERI_H4yz_Pz_Px_S = I_ERI_I4y2z_S_Px_S+ABZ*I_ERI_H4yz_S_Px_S;
  Double I_ERI_H3y2z_Pz_Px_S = I_ERI_I3y3z_S_Px_S+ABZ*I_ERI_H3y2z_S_Px_S;
  Double I_ERI_H2y3z_Pz_Px_S = I_ERI_I2y4z_S_Px_S+ABZ*I_ERI_H2y3z_S_Px_S;
  Double I_ERI_Hy4z_Pz_Px_S = I_ERI_Iy5z_S_Px_S+ABZ*I_ERI_Hy4z_S_Px_S;
  Double I_ERI_H5z_Pz_Px_S = I_ERI_I6z_S_Px_S+ABZ*I_ERI_H5z_S_Px_S;
  Double I_ERI_H5x_Px_Py_S = I_ERI_I6x_S_Py_S+ABX*I_ERI_H5x_S_Py_S;
  Double I_ERI_H4xy_Px_Py_S = I_ERI_I5xy_S_Py_S+ABX*I_ERI_H4xy_S_Py_S;
  Double I_ERI_H4xz_Px_Py_S = I_ERI_I5xz_S_Py_S+ABX*I_ERI_H4xz_S_Py_S;
  Double I_ERI_H3x2y_Px_Py_S = I_ERI_I4x2y_S_Py_S+ABX*I_ERI_H3x2y_S_Py_S;
  Double I_ERI_H3xyz_Px_Py_S = I_ERI_I4xyz_S_Py_S+ABX*I_ERI_H3xyz_S_Py_S;
  Double I_ERI_H3x2z_Px_Py_S = I_ERI_I4x2z_S_Py_S+ABX*I_ERI_H3x2z_S_Py_S;
  Double I_ERI_H2x3y_Px_Py_S = I_ERI_I3x3y_S_Py_S+ABX*I_ERI_H2x3y_S_Py_S;
  Double I_ERI_H2x2yz_Px_Py_S = I_ERI_I3x2yz_S_Py_S+ABX*I_ERI_H2x2yz_S_Py_S;
  Double I_ERI_H2xy2z_Px_Py_S = I_ERI_I3xy2z_S_Py_S+ABX*I_ERI_H2xy2z_S_Py_S;
  Double I_ERI_H2x3z_Px_Py_S = I_ERI_I3x3z_S_Py_S+ABX*I_ERI_H2x3z_S_Py_S;
  Double I_ERI_Hx4y_Px_Py_S = I_ERI_I2x4y_S_Py_S+ABX*I_ERI_Hx4y_S_Py_S;
  Double I_ERI_Hx3yz_Px_Py_S = I_ERI_I2x3yz_S_Py_S+ABX*I_ERI_Hx3yz_S_Py_S;
  Double I_ERI_Hx2y2z_Px_Py_S = I_ERI_I2x2y2z_S_Py_S+ABX*I_ERI_Hx2y2z_S_Py_S;
  Double I_ERI_Hxy3z_Px_Py_S = I_ERI_I2xy3z_S_Py_S+ABX*I_ERI_Hxy3z_S_Py_S;
  Double I_ERI_Hx4z_Px_Py_S = I_ERI_I2x4z_S_Py_S+ABX*I_ERI_Hx4z_S_Py_S;
  Double I_ERI_H5y_Px_Py_S = I_ERI_Ix5y_S_Py_S+ABX*I_ERI_H5y_S_Py_S;
  Double I_ERI_H4yz_Px_Py_S = I_ERI_Ix4yz_S_Py_S+ABX*I_ERI_H4yz_S_Py_S;
  Double I_ERI_H3y2z_Px_Py_S = I_ERI_Ix3y2z_S_Py_S+ABX*I_ERI_H3y2z_S_Py_S;
  Double I_ERI_H2y3z_Px_Py_S = I_ERI_Ix2y3z_S_Py_S+ABX*I_ERI_H2y3z_S_Py_S;
  Double I_ERI_Hy4z_Px_Py_S = I_ERI_Ixy4z_S_Py_S+ABX*I_ERI_Hy4z_S_Py_S;
  Double I_ERI_H5z_Px_Py_S = I_ERI_Ix5z_S_Py_S+ABX*I_ERI_H5z_S_Py_S;
  Double I_ERI_H5x_Py_Py_S = I_ERI_I5xy_S_Py_S+ABY*I_ERI_H5x_S_Py_S;
  Double I_ERI_H4xy_Py_Py_S = I_ERI_I4x2y_S_Py_S+ABY*I_ERI_H4xy_S_Py_S;
  Double I_ERI_H4xz_Py_Py_S = I_ERI_I4xyz_S_Py_S+ABY*I_ERI_H4xz_S_Py_S;
  Double I_ERI_H3x2y_Py_Py_S = I_ERI_I3x3y_S_Py_S+ABY*I_ERI_H3x2y_S_Py_S;
  Double I_ERI_H3xyz_Py_Py_S = I_ERI_I3x2yz_S_Py_S+ABY*I_ERI_H3xyz_S_Py_S;
  Double I_ERI_H3x2z_Py_Py_S = I_ERI_I3xy2z_S_Py_S+ABY*I_ERI_H3x2z_S_Py_S;
  Double I_ERI_H2x3y_Py_Py_S = I_ERI_I2x4y_S_Py_S+ABY*I_ERI_H2x3y_S_Py_S;
  Double I_ERI_H2x2yz_Py_Py_S = I_ERI_I2x3yz_S_Py_S+ABY*I_ERI_H2x2yz_S_Py_S;
  Double I_ERI_H2xy2z_Py_Py_S = I_ERI_I2x2y2z_S_Py_S+ABY*I_ERI_H2xy2z_S_Py_S;
  Double I_ERI_H2x3z_Py_Py_S = I_ERI_I2xy3z_S_Py_S+ABY*I_ERI_H2x3z_S_Py_S;
  Double I_ERI_Hx4y_Py_Py_S = I_ERI_Ix5y_S_Py_S+ABY*I_ERI_Hx4y_S_Py_S;
  Double I_ERI_Hx3yz_Py_Py_S = I_ERI_Ix4yz_S_Py_S+ABY*I_ERI_Hx3yz_S_Py_S;
  Double I_ERI_Hx2y2z_Py_Py_S = I_ERI_Ix3y2z_S_Py_S+ABY*I_ERI_Hx2y2z_S_Py_S;
  Double I_ERI_Hxy3z_Py_Py_S = I_ERI_Ix2y3z_S_Py_S+ABY*I_ERI_Hxy3z_S_Py_S;
  Double I_ERI_Hx4z_Py_Py_S = I_ERI_Ixy4z_S_Py_S+ABY*I_ERI_Hx4z_S_Py_S;
  Double I_ERI_H5y_Py_Py_S = I_ERI_I6y_S_Py_S+ABY*I_ERI_H5y_S_Py_S;
  Double I_ERI_H4yz_Py_Py_S = I_ERI_I5yz_S_Py_S+ABY*I_ERI_H4yz_S_Py_S;
  Double I_ERI_H3y2z_Py_Py_S = I_ERI_I4y2z_S_Py_S+ABY*I_ERI_H3y2z_S_Py_S;
  Double I_ERI_H2y3z_Py_Py_S = I_ERI_I3y3z_S_Py_S+ABY*I_ERI_H2y3z_S_Py_S;
  Double I_ERI_Hy4z_Py_Py_S = I_ERI_I2y4z_S_Py_S+ABY*I_ERI_Hy4z_S_Py_S;
  Double I_ERI_H5z_Py_Py_S = I_ERI_Iy5z_S_Py_S+ABY*I_ERI_H5z_S_Py_S;
  Double I_ERI_H5x_Pz_Py_S = I_ERI_I5xz_S_Py_S+ABZ*I_ERI_H5x_S_Py_S;
  Double I_ERI_H4xy_Pz_Py_S = I_ERI_I4xyz_S_Py_S+ABZ*I_ERI_H4xy_S_Py_S;
  Double I_ERI_H4xz_Pz_Py_S = I_ERI_I4x2z_S_Py_S+ABZ*I_ERI_H4xz_S_Py_S;
  Double I_ERI_H3x2y_Pz_Py_S = I_ERI_I3x2yz_S_Py_S+ABZ*I_ERI_H3x2y_S_Py_S;
  Double I_ERI_H3xyz_Pz_Py_S = I_ERI_I3xy2z_S_Py_S+ABZ*I_ERI_H3xyz_S_Py_S;
  Double I_ERI_H3x2z_Pz_Py_S = I_ERI_I3x3z_S_Py_S+ABZ*I_ERI_H3x2z_S_Py_S;
  Double I_ERI_H2x3y_Pz_Py_S = I_ERI_I2x3yz_S_Py_S+ABZ*I_ERI_H2x3y_S_Py_S;
  Double I_ERI_H2x2yz_Pz_Py_S = I_ERI_I2x2y2z_S_Py_S+ABZ*I_ERI_H2x2yz_S_Py_S;
  Double I_ERI_H2xy2z_Pz_Py_S = I_ERI_I2xy3z_S_Py_S+ABZ*I_ERI_H2xy2z_S_Py_S;
  Double I_ERI_H2x3z_Pz_Py_S = I_ERI_I2x4z_S_Py_S+ABZ*I_ERI_H2x3z_S_Py_S;
  Double I_ERI_Hx4y_Pz_Py_S = I_ERI_Ix4yz_S_Py_S+ABZ*I_ERI_Hx4y_S_Py_S;
  Double I_ERI_Hx3yz_Pz_Py_S = I_ERI_Ix3y2z_S_Py_S+ABZ*I_ERI_Hx3yz_S_Py_S;
  Double I_ERI_Hx2y2z_Pz_Py_S = I_ERI_Ix2y3z_S_Py_S+ABZ*I_ERI_Hx2y2z_S_Py_S;
  Double I_ERI_Hxy3z_Pz_Py_S = I_ERI_Ixy4z_S_Py_S+ABZ*I_ERI_Hxy3z_S_Py_S;
  Double I_ERI_Hx4z_Pz_Py_S = I_ERI_Ix5z_S_Py_S+ABZ*I_ERI_Hx4z_S_Py_S;
  Double I_ERI_H5y_Pz_Py_S = I_ERI_I5yz_S_Py_S+ABZ*I_ERI_H5y_S_Py_S;
  Double I_ERI_H4yz_Pz_Py_S = I_ERI_I4y2z_S_Py_S+ABZ*I_ERI_H4yz_S_Py_S;
  Double I_ERI_H3y2z_Pz_Py_S = I_ERI_I3y3z_S_Py_S+ABZ*I_ERI_H3y2z_S_Py_S;
  Double I_ERI_H2y3z_Pz_Py_S = I_ERI_I2y4z_S_Py_S+ABZ*I_ERI_H2y3z_S_Py_S;
  Double I_ERI_Hy4z_Pz_Py_S = I_ERI_Iy5z_S_Py_S+ABZ*I_ERI_Hy4z_S_Py_S;
  Double I_ERI_H5z_Pz_Py_S = I_ERI_I6z_S_Py_S+ABZ*I_ERI_H5z_S_Py_S;
  Double I_ERI_H5x_Px_Pz_S = I_ERI_I6x_S_Pz_S+ABX*I_ERI_H5x_S_Pz_S;
  Double I_ERI_H4xy_Px_Pz_S = I_ERI_I5xy_S_Pz_S+ABX*I_ERI_H4xy_S_Pz_S;
  Double I_ERI_H4xz_Px_Pz_S = I_ERI_I5xz_S_Pz_S+ABX*I_ERI_H4xz_S_Pz_S;
  Double I_ERI_H3x2y_Px_Pz_S = I_ERI_I4x2y_S_Pz_S+ABX*I_ERI_H3x2y_S_Pz_S;
  Double I_ERI_H3xyz_Px_Pz_S = I_ERI_I4xyz_S_Pz_S+ABX*I_ERI_H3xyz_S_Pz_S;
  Double I_ERI_H3x2z_Px_Pz_S = I_ERI_I4x2z_S_Pz_S+ABX*I_ERI_H3x2z_S_Pz_S;
  Double I_ERI_H2x3y_Px_Pz_S = I_ERI_I3x3y_S_Pz_S+ABX*I_ERI_H2x3y_S_Pz_S;
  Double I_ERI_H2x2yz_Px_Pz_S = I_ERI_I3x2yz_S_Pz_S+ABX*I_ERI_H2x2yz_S_Pz_S;
  Double I_ERI_H2xy2z_Px_Pz_S = I_ERI_I3xy2z_S_Pz_S+ABX*I_ERI_H2xy2z_S_Pz_S;
  Double I_ERI_H2x3z_Px_Pz_S = I_ERI_I3x3z_S_Pz_S+ABX*I_ERI_H2x3z_S_Pz_S;
  Double I_ERI_Hx4y_Px_Pz_S = I_ERI_I2x4y_S_Pz_S+ABX*I_ERI_Hx4y_S_Pz_S;
  Double I_ERI_Hx3yz_Px_Pz_S = I_ERI_I2x3yz_S_Pz_S+ABX*I_ERI_Hx3yz_S_Pz_S;
  Double I_ERI_Hx2y2z_Px_Pz_S = I_ERI_I2x2y2z_S_Pz_S+ABX*I_ERI_Hx2y2z_S_Pz_S;
  Double I_ERI_Hxy3z_Px_Pz_S = I_ERI_I2xy3z_S_Pz_S+ABX*I_ERI_Hxy3z_S_Pz_S;
  Double I_ERI_Hx4z_Px_Pz_S = I_ERI_I2x4z_S_Pz_S+ABX*I_ERI_Hx4z_S_Pz_S;
  Double I_ERI_H5y_Px_Pz_S = I_ERI_Ix5y_S_Pz_S+ABX*I_ERI_H5y_S_Pz_S;
  Double I_ERI_H4yz_Px_Pz_S = I_ERI_Ix4yz_S_Pz_S+ABX*I_ERI_H4yz_S_Pz_S;
  Double I_ERI_H3y2z_Px_Pz_S = I_ERI_Ix3y2z_S_Pz_S+ABX*I_ERI_H3y2z_S_Pz_S;
  Double I_ERI_H2y3z_Px_Pz_S = I_ERI_Ix2y3z_S_Pz_S+ABX*I_ERI_H2y3z_S_Pz_S;
  Double I_ERI_Hy4z_Px_Pz_S = I_ERI_Ixy4z_S_Pz_S+ABX*I_ERI_Hy4z_S_Pz_S;
  Double I_ERI_H5z_Px_Pz_S = I_ERI_Ix5z_S_Pz_S+ABX*I_ERI_H5z_S_Pz_S;
  Double I_ERI_H5x_Py_Pz_S = I_ERI_I5xy_S_Pz_S+ABY*I_ERI_H5x_S_Pz_S;
  Double I_ERI_H4xy_Py_Pz_S = I_ERI_I4x2y_S_Pz_S+ABY*I_ERI_H4xy_S_Pz_S;
  Double I_ERI_H4xz_Py_Pz_S = I_ERI_I4xyz_S_Pz_S+ABY*I_ERI_H4xz_S_Pz_S;
  Double I_ERI_H3x2y_Py_Pz_S = I_ERI_I3x3y_S_Pz_S+ABY*I_ERI_H3x2y_S_Pz_S;
  Double I_ERI_H3xyz_Py_Pz_S = I_ERI_I3x2yz_S_Pz_S+ABY*I_ERI_H3xyz_S_Pz_S;
  Double I_ERI_H3x2z_Py_Pz_S = I_ERI_I3xy2z_S_Pz_S+ABY*I_ERI_H3x2z_S_Pz_S;
  Double I_ERI_H2x3y_Py_Pz_S = I_ERI_I2x4y_S_Pz_S+ABY*I_ERI_H2x3y_S_Pz_S;
  Double I_ERI_H2x2yz_Py_Pz_S = I_ERI_I2x3yz_S_Pz_S+ABY*I_ERI_H2x2yz_S_Pz_S;
  Double I_ERI_H2xy2z_Py_Pz_S = I_ERI_I2x2y2z_S_Pz_S+ABY*I_ERI_H2xy2z_S_Pz_S;
  Double I_ERI_H2x3z_Py_Pz_S = I_ERI_I2xy3z_S_Pz_S+ABY*I_ERI_H2x3z_S_Pz_S;
  Double I_ERI_Hx4y_Py_Pz_S = I_ERI_Ix5y_S_Pz_S+ABY*I_ERI_Hx4y_S_Pz_S;
  Double I_ERI_Hx3yz_Py_Pz_S = I_ERI_Ix4yz_S_Pz_S+ABY*I_ERI_Hx3yz_S_Pz_S;
  Double I_ERI_Hx2y2z_Py_Pz_S = I_ERI_Ix3y2z_S_Pz_S+ABY*I_ERI_Hx2y2z_S_Pz_S;
  Double I_ERI_Hxy3z_Py_Pz_S = I_ERI_Ix2y3z_S_Pz_S+ABY*I_ERI_Hxy3z_S_Pz_S;
  Double I_ERI_Hx4z_Py_Pz_S = I_ERI_Ixy4z_S_Pz_S+ABY*I_ERI_Hx4z_S_Pz_S;
  Double I_ERI_H5y_Py_Pz_S = I_ERI_I6y_S_Pz_S+ABY*I_ERI_H5y_S_Pz_S;
  Double I_ERI_H4yz_Py_Pz_S = I_ERI_I5yz_S_Pz_S+ABY*I_ERI_H4yz_S_Pz_S;
  Double I_ERI_H3y2z_Py_Pz_S = I_ERI_I4y2z_S_Pz_S+ABY*I_ERI_H3y2z_S_Pz_S;
  Double I_ERI_H2y3z_Py_Pz_S = I_ERI_I3y3z_S_Pz_S+ABY*I_ERI_H2y3z_S_Pz_S;
  Double I_ERI_Hy4z_Py_Pz_S = I_ERI_I2y4z_S_Pz_S+ABY*I_ERI_Hy4z_S_Pz_S;
  Double I_ERI_H5z_Py_Pz_S = I_ERI_Iy5z_S_Pz_S+ABY*I_ERI_H5z_S_Pz_S;
  Double I_ERI_H5x_Pz_Pz_S = I_ERI_I5xz_S_Pz_S+ABZ*I_ERI_H5x_S_Pz_S;
  Double I_ERI_H4xy_Pz_Pz_S = I_ERI_I4xyz_S_Pz_S+ABZ*I_ERI_H4xy_S_Pz_S;
  Double I_ERI_H4xz_Pz_Pz_S = I_ERI_I4x2z_S_Pz_S+ABZ*I_ERI_H4xz_S_Pz_S;
  Double I_ERI_H3x2y_Pz_Pz_S = I_ERI_I3x2yz_S_Pz_S+ABZ*I_ERI_H3x2y_S_Pz_S;
  Double I_ERI_H3xyz_Pz_Pz_S = I_ERI_I3xy2z_S_Pz_S+ABZ*I_ERI_H3xyz_S_Pz_S;
  Double I_ERI_H3x2z_Pz_Pz_S = I_ERI_I3x3z_S_Pz_S+ABZ*I_ERI_H3x2z_S_Pz_S;
  Double I_ERI_H2x3y_Pz_Pz_S = I_ERI_I2x3yz_S_Pz_S+ABZ*I_ERI_H2x3y_S_Pz_S;
  Double I_ERI_H2x2yz_Pz_Pz_S = I_ERI_I2x2y2z_S_Pz_S+ABZ*I_ERI_H2x2yz_S_Pz_S;
  Double I_ERI_H2xy2z_Pz_Pz_S = I_ERI_I2xy3z_S_Pz_S+ABZ*I_ERI_H2xy2z_S_Pz_S;
  Double I_ERI_H2x3z_Pz_Pz_S = I_ERI_I2x4z_S_Pz_S+ABZ*I_ERI_H2x3z_S_Pz_S;
  Double I_ERI_Hx4y_Pz_Pz_S = I_ERI_Ix4yz_S_Pz_S+ABZ*I_ERI_Hx4y_S_Pz_S;
  Double I_ERI_Hx3yz_Pz_Pz_S = I_ERI_Ix3y2z_S_Pz_S+ABZ*I_ERI_Hx3yz_S_Pz_S;
  Double I_ERI_Hx2y2z_Pz_Pz_S = I_ERI_Ix2y3z_S_Pz_S+ABZ*I_ERI_Hx2y2z_S_Pz_S;
  Double I_ERI_Hxy3z_Pz_Pz_S = I_ERI_Ixy4z_S_Pz_S+ABZ*I_ERI_Hxy3z_S_Pz_S;
  Double I_ERI_Hx4z_Pz_Pz_S = I_ERI_Ix5z_S_Pz_S+ABZ*I_ERI_Hx4z_S_Pz_S;
  Double I_ERI_H5y_Pz_Pz_S = I_ERI_I5yz_S_Pz_S+ABZ*I_ERI_H5y_S_Pz_S;
  Double I_ERI_H4yz_Pz_Pz_S = I_ERI_I4y2z_S_Pz_S+ABZ*I_ERI_H4yz_S_Pz_S;
  Double I_ERI_H3y2z_Pz_Pz_S = I_ERI_I3y3z_S_Pz_S+ABZ*I_ERI_H3y2z_S_Pz_S;
  Double I_ERI_H2y3z_Pz_Pz_S = I_ERI_I2y4z_S_Pz_S+ABZ*I_ERI_H2y3z_S_Pz_S;
  Double I_ERI_Hy4z_Pz_Pz_S = I_ERI_Iy5z_S_Pz_S+ABZ*I_ERI_Hy4z_S_Pz_S;
  Double I_ERI_H5z_Pz_Pz_S = I_ERI_I6z_S_Pz_S+ABZ*I_ERI_H5z_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_G_D_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 135 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_P_P_S
   * RHS shell quartet name: SQ_ERI_G_P_P_S
   ************************************************************/
  Double I_ERI_G4x_D2x_Px_S = I_ERI_H5x_Px_Px_S+ABX*I_ERI_G4x_Px_Px_S;
  Double I_ERI_G3xy_D2x_Px_S = I_ERI_H4xy_Px_Px_S+ABX*I_ERI_G3xy_Px_Px_S;
  Double I_ERI_G3xz_D2x_Px_S = I_ERI_H4xz_Px_Px_S+ABX*I_ERI_G3xz_Px_Px_S;
  Double I_ERI_G2x2y_D2x_Px_S = I_ERI_H3x2y_Px_Px_S+ABX*I_ERI_G2x2y_Px_Px_S;
  Double I_ERI_G2xyz_D2x_Px_S = I_ERI_H3xyz_Px_Px_S+ABX*I_ERI_G2xyz_Px_Px_S;
  Double I_ERI_G2x2z_D2x_Px_S = I_ERI_H3x2z_Px_Px_S+ABX*I_ERI_G2x2z_Px_Px_S;
  Double I_ERI_Gx3y_D2x_Px_S = I_ERI_H2x3y_Px_Px_S+ABX*I_ERI_Gx3y_Px_Px_S;
  Double I_ERI_Gx2yz_D2x_Px_S = I_ERI_H2x2yz_Px_Px_S+ABX*I_ERI_Gx2yz_Px_Px_S;
  Double I_ERI_Gxy2z_D2x_Px_S = I_ERI_H2xy2z_Px_Px_S+ABX*I_ERI_Gxy2z_Px_Px_S;
  Double I_ERI_Gx3z_D2x_Px_S = I_ERI_H2x3z_Px_Px_S+ABX*I_ERI_Gx3z_Px_Px_S;
  Double I_ERI_G4y_D2x_Px_S = I_ERI_Hx4y_Px_Px_S+ABX*I_ERI_G4y_Px_Px_S;
  Double I_ERI_G3yz_D2x_Px_S = I_ERI_Hx3yz_Px_Px_S+ABX*I_ERI_G3yz_Px_Px_S;
  Double I_ERI_G2y2z_D2x_Px_S = I_ERI_Hx2y2z_Px_Px_S+ABX*I_ERI_G2y2z_Px_Px_S;
  Double I_ERI_Gy3z_D2x_Px_S = I_ERI_Hxy3z_Px_Px_S+ABX*I_ERI_Gy3z_Px_Px_S;
  Double I_ERI_G4z_D2x_Px_S = I_ERI_Hx4z_Px_Px_S+ABX*I_ERI_G4z_Px_Px_S;
  Double I_ERI_G4x_D2y_Px_S = I_ERI_H4xy_Py_Px_S+ABY*I_ERI_G4x_Py_Px_S;
  Double I_ERI_G3xy_D2y_Px_S = I_ERI_H3x2y_Py_Px_S+ABY*I_ERI_G3xy_Py_Px_S;
  Double I_ERI_G3xz_D2y_Px_S = I_ERI_H3xyz_Py_Px_S+ABY*I_ERI_G3xz_Py_Px_S;
  Double I_ERI_G2x2y_D2y_Px_S = I_ERI_H2x3y_Py_Px_S+ABY*I_ERI_G2x2y_Py_Px_S;
  Double I_ERI_G2xyz_D2y_Px_S = I_ERI_H2x2yz_Py_Px_S+ABY*I_ERI_G2xyz_Py_Px_S;
  Double I_ERI_G2x2z_D2y_Px_S = I_ERI_H2xy2z_Py_Px_S+ABY*I_ERI_G2x2z_Py_Px_S;
  Double I_ERI_Gx3y_D2y_Px_S = I_ERI_Hx4y_Py_Px_S+ABY*I_ERI_Gx3y_Py_Px_S;
  Double I_ERI_Gx2yz_D2y_Px_S = I_ERI_Hx3yz_Py_Px_S+ABY*I_ERI_Gx2yz_Py_Px_S;
  Double I_ERI_Gxy2z_D2y_Px_S = I_ERI_Hx2y2z_Py_Px_S+ABY*I_ERI_Gxy2z_Py_Px_S;
  Double I_ERI_Gx3z_D2y_Px_S = I_ERI_Hxy3z_Py_Px_S+ABY*I_ERI_Gx3z_Py_Px_S;
  Double I_ERI_G4y_D2y_Px_S = I_ERI_H5y_Py_Px_S+ABY*I_ERI_G4y_Py_Px_S;
  Double I_ERI_G3yz_D2y_Px_S = I_ERI_H4yz_Py_Px_S+ABY*I_ERI_G3yz_Py_Px_S;
  Double I_ERI_G2y2z_D2y_Px_S = I_ERI_H3y2z_Py_Px_S+ABY*I_ERI_G2y2z_Py_Px_S;
  Double I_ERI_Gy3z_D2y_Px_S = I_ERI_H2y3z_Py_Px_S+ABY*I_ERI_Gy3z_Py_Px_S;
  Double I_ERI_G4z_D2y_Px_S = I_ERI_Hy4z_Py_Px_S+ABY*I_ERI_G4z_Py_Px_S;
  Double I_ERI_G4x_D2z_Px_S = I_ERI_H4xz_Pz_Px_S+ABZ*I_ERI_G4x_Pz_Px_S;
  Double I_ERI_G3xy_D2z_Px_S = I_ERI_H3xyz_Pz_Px_S+ABZ*I_ERI_G3xy_Pz_Px_S;
  Double I_ERI_G3xz_D2z_Px_S = I_ERI_H3x2z_Pz_Px_S+ABZ*I_ERI_G3xz_Pz_Px_S;
  Double I_ERI_G2x2y_D2z_Px_S = I_ERI_H2x2yz_Pz_Px_S+ABZ*I_ERI_G2x2y_Pz_Px_S;
  Double I_ERI_G2xyz_D2z_Px_S = I_ERI_H2xy2z_Pz_Px_S+ABZ*I_ERI_G2xyz_Pz_Px_S;
  Double I_ERI_G2x2z_D2z_Px_S = I_ERI_H2x3z_Pz_Px_S+ABZ*I_ERI_G2x2z_Pz_Px_S;
  Double I_ERI_Gx3y_D2z_Px_S = I_ERI_Hx3yz_Pz_Px_S+ABZ*I_ERI_Gx3y_Pz_Px_S;
  Double I_ERI_Gx2yz_D2z_Px_S = I_ERI_Hx2y2z_Pz_Px_S+ABZ*I_ERI_Gx2yz_Pz_Px_S;
  Double I_ERI_Gxy2z_D2z_Px_S = I_ERI_Hxy3z_Pz_Px_S+ABZ*I_ERI_Gxy2z_Pz_Px_S;
  Double I_ERI_Gx3z_D2z_Px_S = I_ERI_Hx4z_Pz_Px_S+ABZ*I_ERI_Gx3z_Pz_Px_S;
  Double I_ERI_G4y_D2z_Px_S = I_ERI_H4yz_Pz_Px_S+ABZ*I_ERI_G4y_Pz_Px_S;
  Double I_ERI_G3yz_D2z_Px_S = I_ERI_H3y2z_Pz_Px_S+ABZ*I_ERI_G3yz_Pz_Px_S;
  Double I_ERI_G2y2z_D2z_Px_S = I_ERI_H2y3z_Pz_Px_S+ABZ*I_ERI_G2y2z_Pz_Px_S;
  Double I_ERI_Gy3z_D2z_Px_S = I_ERI_Hy4z_Pz_Px_S+ABZ*I_ERI_Gy3z_Pz_Px_S;
  Double I_ERI_G4z_D2z_Px_S = I_ERI_H5z_Pz_Px_S+ABZ*I_ERI_G4z_Pz_Px_S;
  Double I_ERI_G4x_D2x_Py_S = I_ERI_H5x_Px_Py_S+ABX*I_ERI_G4x_Px_Py_S;
  Double I_ERI_G3xy_D2x_Py_S = I_ERI_H4xy_Px_Py_S+ABX*I_ERI_G3xy_Px_Py_S;
  Double I_ERI_G3xz_D2x_Py_S = I_ERI_H4xz_Px_Py_S+ABX*I_ERI_G3xz_Px_Py_S;
  Double I_ERI_G2x2y_D2x_Py_S = I_ERI_H3x2y_Px_Py_S+ABX*I_ERI_G2x2y_Px_Py_S;
  Double I_ERI_G2xyz_D2x_Py_S = I_ERI_H3xyz_Px_Py_S+ABX*I_ERI_G2xyz_Px_Py_S;
  Double I_ERI_G2x2z_D2x_Py_S = I_ERI_H3x2z_Px_Py_S+ABX*I_ERI_G2x2z_Px_Py_S;
  Double I_ERI_Gx3y_D2x_Py_S = I_ERI_H2x3y_Px_Py_S+ABX*I_ERI_Gx3y_Px_Py_S;
  Double I_ERI_Gx2yz_D2x_Py_S = I_ERI_H2x2yz_Px_Py_S+ABX*I_ERI_Gx2yz_Px_Py_S;
  Double I_ERI_Gxy2z_D2x_Py_S = I_ERI_H2xy2z_Px_Py_S+ABX*I_ERI_Gxy2z_Px_Py_S;
  Double I_ERI_Gx3z_D2x_Py_S = I_ERI_H2x3z_Px_Py_S+ABX*I_ERI_Gx3z_Px_Py_S;
  Double I_ERI_G4y_D2x_Py_S = I_ERI_Hx4y_Px_Py_S+ABX*I_ERI_G4y_Px_Py_S;
  Double I_ERI_G3yz_D2x_Py_S = I_ERI_Hx3yz_Px_Py_S+ABX*I_ERI_G3yz_Px_Py_S;
  Double I_ERI_G2y2z_D2x_Py_S = I_ERI_Hx2y2z_Px_Py_S+ABX*I_ERI_G2y2z_Px_Py_S;
  Double I_ERI_Gy3z_D2x_Py_S = I_ERI_Hxy3z_Px_Py_S+ABX*I_ERI_Gy3z_Px_Py_S;
  Double I_ERI_G4z_D2x_Py_S = I_ERI_Hx4z_Px_Py_S+ABX*I_ERI_G4z_Px_Py_S;
  Double I_ERI_G4x_D2y_Py_S = I_ERI_H4xy_Py_Py_S+ABY*I_ERI_G4x_Py_Py_S;
  Double I_ERI_G3xy_D2y_Py_S = I_ERI_H3x2y_Py_Py_S+ABY*I_ERI_G3xy_Py_Py_S;
  Double I_ERI_G3xz_D2y_Py_S = I_ERI_H3xyz_Py_Py_S+ABY*I_ERI_G3xz_Py_Py_S;
  Double I_ERI_G2x2y_D2y_Py_S = I_ERI_H2x3y_Py_Py_S+ABY*I_ERI_G2x2y_Py_Py_S;
  Double I_ERI_G2xyz_D2y_Py_S = I_ERI_H2x2yz_Py_Py_S+ABY*I_ERI_G2xyz_Py_Py_S;
  Double I_ERI_G2x2z_D2y_Py_S = I_ERI_H2xy2z_Py_Py_S+ABY*I_ERI_G2x2z_Py_Py_S;
  Double I_ERI_Gx3y_D2y_Py_S = I_ERI_Hx4y_Py_Py_S+ABY*I_ERI_Gx3y_Py_Py_S;
  Double I_ERI_Gx2yz_D2y_Py_S = I_ERI_Hx3yz_Py_Py_S+ABY*I_ERI_Gx2yz_Py_Py_S;
  Double I_ERI_Gxy2z_D2y_Py_S = I_ERI_Hx2y2z_Py_Py_S+ABY*I_ERI_Gxy2z_Py_Py_S;
  Double I_ERI_Gx3z_D2y_Py_S = I_ERI_Hxy3z_Py_Py_S+ABY*I_ERI_Gx3z_Py_Py_S;
  Double I_ERI_G4y_D2y_Py_S = I_ERI_H5y_Py_Py_S+ABY*I_ERI_G4y_Py_Py_S;
  Double I_ERI_G3yz_D2y_Py_S = I_ERI_H4yz_Py_Py_S+ABY*I_ERI_G3yz_Py_Py_S;
  Double I_ERI_G2y2z_D2y_Py_S = I_ERI_H3y2z_Py_Py_S+ABY*I_ERI_G2y2z_Py_Py_S;
  Double I_ERI_Gy3z_D2y_Py_S = I_ERI_H2y3z_Py_Py_S+ABY*I_ERI_Gy3z_Py_Py_S;
  Double I_ERI_G4z_D2y_Py_S = I_ERI_Hy4z_Py_Py_S+ABY*I_ERI_G4z_Py_Py_S;
  Double I_ERI_G4x_D2z_Py_S = I_ERI_H4xz_Pz_Py_S+ABZ*I_ERI_G4x_Pz_Py_S;
  Double I_ERI_G3xy_D2z_Py_S = I_ERI_H3xyz_Pz_Py_S+ABZ*I_ERI_G3xy_Pz_Py_S;
  Double I_ERI_G3xz_D2z_Py_S = I_ERI_H3x2z_Pz_Py_S+ABZ*I_ERI_G3xz_Pz_Py_S;
  Double I_ERI_G2x2y_D2z_Py_S = I_ERI_H2x2yz_Pz_Py_S+ABZ*I_ERI_G2x2y_Pz_Py_S;
  Double I_ERI_G2xyz_D2z_Py_S = I_ERI_H2xy2z_Pz_Py_S+ABZ*I_ERI_G2xyz_Pz_Py_S;
  Double I_ERI_G2x2z_D2z_Py_S = I_ERI_H2x3z_Pz_Py_S+ABZ*I_ERI_G2x2z_Pz_Py_S;
  Double I_ERI_Gx3y_D2z_Py_S = I_ERI_Hx3yz_Pz_Py_S+ABZ*I_ERI_Gx3y_Pz_Py_S;
  Double I_ERI_Gx2yz_D2z_Py_S = I_ERI_Hx2y2z_Pz_Py_S+ABZ*I_ERI_Gx2yz_Pz_Py_S;
  Double I_ERI_Gxy2z_D2z_Py_S = I_ERI_Hxy3z_Pz_Py_S+ABZ*I_ERI_Gxy2z_Pz_Py_S;
  Double I_ERI_Gx3z_D2z_Py_S = I_ERI_Hx4z_Pz_Py_S+ABZ*I_ERI_Gx3z_Pz_Py_S;
  Double I_ERI_G4y_D2z_Py_S = I_ERI_H4yz_Pz_Py_S+ABZ*I_ERI_G4y_Pz_Py_S;
  Double I_ERI_G3yz_D2z_Py_S = I_ERI_H3y2z_Pz_Py_S+ABZ*I_ERI_G3yz_Pz_Py_S;
  Double I_ERI_G2y2z_D2z_Py_S = I_ERI_H2y3z_Pz_Py_S+ABZ*I_ERI_G2y2z_Pz_Py_S;
  Double I_ERI_Gy3z_D2z_Py_S = I_ERI_Hy4z_Pz_Py_S+ABZ*I_ERI_Gy3z_Pz_Py_S;
  Double I_ERI_G4z_D2z_Py_S = I_ERI_H5z_Pz_Py_S+ABZ*I_ERI_G4z_Pz_Py_S;
  Double I_ERI_G4x_D2x_Pz_S = I_ERI_H5x_Px_Pz_S+ABX*I_ERI_G4x_Px_Pz_S;
  Double I_ERI_G3xy_D2x_Pz_S = I_ERI_H4xy_Px_Pz_S+ABX*I_ERI_G3xy_Px_Pz_S;
  Double I_ERI_G3xz_D2x_Pz_S = I_ERI_H4xz_Px_Pz_S+ABX*I_ERI_G3xz_Px_Pz_S;
  Double I_ERI_G2x2y_D2x_Pz_S = I_ERI_H3x2y_Px_Pz_S+ABX*I_ERI_G2x2y_Px_Pz_S;
  Double I_ERI_G2xyz_D2x_Pz_S = I_ERI_H3xyz_Px_Pz_S+ABX*I_ERI_G2xyz_Px_Pz_S;
  Double I_ERI_G2x2z_D2x_Pz_S = I_ERI_H3x2z_Px_Pz_S+ABX*I_ERI_G2x2z_Px_Pz_S;
  Double I_ERI_Gx3y_D2x_Pz_S = I_ERI_H2x3y_Px_Pz_S+ABX*I_ERI_Gx3y_Px_Pz_S;
  Double I_ERI_Gx2yz_D2x_Pz_S = I_ERI_H2x2yz_Px_Pz_S+ABX*I_ERI_Gx2yz_Px_Pz_S;
  Double I_ERI_Gxy2z_D2x_Pz_S = I_ERI_H2xy2z_Px_Pz_S+ABX*I_ERI_Gxy2z_Px_Pz_S;
  Double I_ERI_Gx3z_D2x_Pz_S = I_ERI_H2x3z_Px_Pz_S+ABX*I_ERI_Gx3z_Px_Pz_S;
  Double I_ERI_G4y_D2x_Pz_S = I_ERI_Hx4y_Px_Pz_S+ABX*I_ERI_G4y_Px_Pz_S;
  Double I_ERI_G3yz_D2x_Pz_S = I_ERI_Hx3yz_Px_Pz_S+ABX*I_ERI_G3yz_Px_Pz_S;
  Double I_ERI_G2y2z_D2x_Pz_S = I_ERI_Hx2y2z_Px_Pz_S+ABX*I_ERI_G2y2z_Px_Pz_S;
  Double I_ERI_Gy3z_D2x_Pz_S = I_ERI_Hxy3z_Px_Pz_S+ABX*I_ERI_Gy3z_Px_Pz_S;
  Double I_ERI_G4z_D2x_Pz_S = I_ERI_Hx4z_Px_Pz_S+ABX*I_ERI_G4z_Px_Pz_S;
  Double I_ERI_G4x_D2y_Pz_S = I_ERI_H4xy_Py_Pz_S+ABY*I_ERI_G4x_Py_Pz_S;
  Double I_ERI_G3xy_D2y_Pz_S = I_ERI_H3x2y_Py_Pz_S+ABY*I_ERI_G3xy_Py_Pz_S;
  Double I_ERI_G3xz_D2y_Pz_S = I_ERI_H3xyz_Py_Pz_S+ABY*I_ERI_G3xz_Py_Pz_S;
  Double I_ERI_G2x2y_D2y_Pz_S = I_ERI_H2x3y_Py_Pz_S+ABY*I_ERI_G2x2y_Py_Pz_S;
  Double I_ERI_G2xyz_D2y_Pz_S = I_ERI_H2x2yz_Py_Pz_S+ABY*I_ERI_G2xyz_Py_Pz_S;
  Double I_ERI_G2x2z_D2y_Pz_S = I_ERI_H2xy2z_Py_Pz_S+ABY*I_ERI_G2x2z_Py_Pz_S;
  Double I_ERI_Gx3y_D2y_Pz_S = I_ERI_Hx4y_Py_Pz_S+ABY*I_ERI_Gx3y_Py_Pz_S;
  Double I_ERI_Gx2yz_D2y_Pz_S = I_ERI_Hx3yz_Py_Pz_S+ABY*I_ERI_Gx2yz_Py_Pz_S;
  Double I_ERI_Gxy2z_D2y_Pz_S = I_ERI_Hx2y2z_Py_Pz_S+ABY*I_ERI_Gxy2z_Py_Pz_S;
  Double I_ERI_Gx3z_D2y_Pz_S = I_ERI_Hxy3z_Py_Pz_S+ABY*I_ERI_Gx3z_Py_Pz_S;
  Double I_ERI_G4y_D2y_Pz_S = I_ERI_H5y_Py_Pz_S+ABY*I_ERI_G4y_Py_Pz_S;
  Double I_ERI_G3yz_D2y_Pz_S = I_ERI_H4yz_Py_Pz_S+ABY*I_ERI_G3yz_Py_Pz_S;
  Double I_ERI_G2y2z_D2y_Pz_S = I_ERI_H3y2z_Py_Pz_S+ABY*I_ERI_G2y2z_Py_Pz_S;
  Double I_ERI_Gy3z_D2y_Pz_S = I_ERI_H2y3z_Py_Pz_S+ABY*I_ERI_Gy3z_Py_Pz_S;
  Double I_ERI_G4z_D2y_Pz_S = I_ERI_Hy4z_Py_Pz_S+ABY*I_ERI_G4z_Py_Pz_S;
  Double I_ERI_G4x_D2z_Pz_S = I_ERI_H4xz_Pz_Pz_S+ABZ*I_ERI_G4x_Pz_Pz_S;
  Double I_ERI_G3xy_D2z_Pz_S = I_ERI_H3xyz_Pz_Pz_S+ABZ*I_ERI_G3xy_Pz_Pz_S;
  Double I_ERI_G3xz_D2z_Pz_S = I_ERI_H3x2z_Pz_Pz_S+ABZ*I_ERI_G3xz_Pz_Pz_S;
  Double I_ERI_G2x2y_D2z_Pz_S = I_ERI_H2x2yz_Pz_Pz_S+ABZ*I_ERI_G2x2y_Pz_Pz_S;
  Double I_ERI_G2xyz_D2z_Pz_S = I_ERI_H2xy2z_Pz_Pz_S+ABZ*I_ERI_G2xyz_Pz_Pz_S;
  Double I_ERI_G2x2z_D2z_Pz_S = I_ERI_H2x3z_Pz_Pz_S+ABZ*I_ERI_G2x2z_Pz_Pz_S;
  Double I_ERI_Gx3y_D2z_Pz_S = I_ERI_Hx3yz_Pz_Pz_S+ABZ*I_ERI_Gx3y_Pz_Pz_S;
  Double I_ERI_Gx2yz_D2z_Pz_S = I_ERI_Hx2y2z_Pz_Pz_S+ABZ*I_ERI_Gx2yz_Pz_Pz_S;
  Double I_ERI_Gxy2z_D2z_Pz_S = I_ERI_Hxy3z_Pz_Pz_S+ABZ*I_ERI_Gxy2z_Pz_Pz_S;
  Double I_ERI_Gx3z_D2z_Pz_S = I_ERI_Hx4z_Pz_Pz_S+ABZ*I_ERI_Gx3z_Pz_Pz_S;
  Double I_ERI_G4y_D2z_Pz_S = I_ERI_H4yz_Pz_Pz_S+ABZ*I_ERI_G4y_Pz_Pz_S;
  Double I_ERI_G3yz_D2z_Pz_S = I_ERI_H3y2z_Pz_Pz_S+ABZ*I_ERI_G3yz_Pz_Pz_S;
  Double I_ERI_G2y2z_D2z_Pz_S = I_ERI_H2y3z_Pz_Pz_S+ABZ*I_ERI_G2y2z_Pz_Pz_S;
  Double I_ERI_Gy3z_D2z_Pz_S = I_ERI_Hy4z_Pz_Pz_S+ABZ*I_ERI_Gy3z_Pz_Pz_S;
  Double I_ERI_G4z_D2z_Pz_S = I_ERI_H5z_Pz_Pz_S+ABZ*I_ERI_G4z_Pz_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_I_P_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 9 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_K_S_P_S
   * RHS shell quartet name: SQ_ERI_I_S_P_S
   ************************************************************/
  Double I_ERI_I6x_Px_Px_S = I_ERI_K7x_S_Px_S+ABX*I_ERI_I6x_S_Px_S;
  Double I_ERI_I5xy_Px_Px_S = I_ERI_K6xy_S_Px_S+ABX*I_ERI_I5xy_S_Px_S;
  Double I_ERI_I5xz_Px_Px_S = I_ERI_K6xz_S_Px_S+ABX*I_ERI_I5xz_S_Px_S;
  Double I_ERI_I4x2y_Px_Px_S = I_ERI_K5x2y_S_Px_S+ABX*I_ERI_I4x2y_S_Px_S;
  Double I_ERI_I4xyz_Px_Px_S = I_ERI_K5xyz_S_Px_S+ABX*I_ERI_I4xyz_S_Px_S;
  Double I_ERI_I4x2z_Px_Px_S = I_ERI_K5x2z_S_Px_S+ABX*I_ERI_I4x2z_S_Px_S;
  Double I_ERI_I3x3y_Px_Px_S = I_ERI_K4x3y_S_Px_S+ABX*I_ERI_I3x3y_S_Px_S;
  Double I_ERI_I3x2yz_Px_Px_S = I_ERI_K4x2yz_S_Px_S+ABX*I_ERI_I3x2yz_S_Px_S;
  Double I_ERI_I3xy2z_Px_Px_S = I_ERI_K4xy2z_S_Px_S+ABX*I_ERI_I3xy2z_S_Px_S;
  Double I_ERI_I3x3z_Px_Px_S = I_ERI_K4x3z_S_Px_S+ABX*I_ERI_I3x3z_S_Px_S;
  Double I_ERI_I2x4y_Px_Px_S = I_ERI_K3x4y_S_Px_S+ABX*I_ERI_I2x4y_S_Px_S;
  Double I_ERI_I2x3yz_Px_Px_S = I_ERI_K3x3yz_S_Px_S+ABX*I_ERI_I2x3yz_S_Px_S;
  Double I_ERI_I2x2y2z_Px_Px_S = I_ERI_K3x2y2z_S_Px_S+ABX*I_ERI_I2x2y2z_S_Px_S;
  Double I_ERI_I2xy3z_Px_Px_S = I_ERI_K3xy3z_S_Px_S+ABX*I_ERI_I2xy3z_S_Px_S;
  Double I_ERI_I2x4z_Px_Px_S = I_ERI_K3x4z_S_Px_S+ABX*I_ERI_I2x4z_S_Px_S;
  Double I_ERI_Ix5y_Px_Px_S = I_ERI_K2x5y_S_Px_S+ABX*I_ERI_Ix5y_S_Px_S;
  Double I_ERI_Ix4yz_Px_Px_S = I_ERI_K2x4yz_S_Px_S+ABX*I_ERI_Ix4yz_S_Px_S;
  Double I_ERI_Ix3y2z_Px_Px_S = I_ERI_K2x3y2z_S_Px_S+ABX*I_ERI_Ix3y2z_S_Px_S;
  Double I_ERI_Ix2y3z_Px_Px_S = I_ERI_K2x2y3z_S_Px_S+ABX*I_ERI_Ix2y3z_S_Px_S;
  Double I_ERI_Ixy4z_Px_Px_S = I_ERI_K2xy4z_S_Px_S+ABX*I_ERI_Ixy4z_S_Px_S;
  Double I_ERI_Ix5z_Px_Px_S = I_ERI_K2x5z_S_Px_S+ABX*I_ERI_Ix5z_S_Px_S;
  Double I_ERI_I6y_Px_Px_S = I_ERI_Kx6y_S_Px_S+ABX*I_ERI_I6y_S_Px_S;
  Double I_ERI_I5yz_Px_Px_S = I_ERI_Kx5yz_S_Px_S+ABX*I_ERI_I5yz_S_Px_S;
  Double I_ERI_I4y2z_Px_Px_S = I_ERI_Kx4y2z_S_Px_S+ABX*I_ERI_I4y2z_S_Px_S;
  Double I_ERI_I3y3z_Px_Px_S = I_ERI_Kx3y3z_S_Px_S+ABX*I_ERI_I3y3z_S_Px_S;
  Double I_ERI_I2y4z_Px_Px_S = I_ERI_Kx2y4z_S_Px_S+ABX*I_ERI_I2y4z_S_Px_S;
  Double I_ERI_Iy5z_Px_Px_S = I_ERI_Kxy5z_S_Px_S+ABX*I_ERI_Iy5z_S_Px_S;
  Double I_ERI_I6z_Px_Px_S = I_ERI_Kx6z_S_Px_S+ABX*I_ERI_I6z_S_Px_S;
  Double I_ERI_I5xy_Py_Px_S = I_ERI_K5x2y_S_Px_S+ABY*I_ERI_I5xy_S_Px_S;
  Double I_ERI_I5xz_Py_Px_S = I_ERI_K5xyz_S_Px_S+ABY*I_ERI_I5xz_S_Px_S;
  Double I_ERI_I4x2y_Py_Px_S = I_ERI_K4x3y_S_Px_S+ABY*I_ERI_I4x2y_S_Px_S;
  Double I_ERI_I4xyz_Py_Px_S = I_ERI_K4x2yz_S_Px_S+ABY*I_ERI_I4xyz_S_Px_S;
  Double I_ERI_I4x2z_Py_Px_S = I_ERI_K4xy2z_S_Px_S+ABY*I_ERI_I4x2z_S_Px_S;
  Double I_ERI_I3x3y_Py_Px_S = I_ERI_K3x4y_S_Px_S+ABY*I_ERI_I3x3y_S_Px_S;
  Double I_ERI_I3x2yz_Py_Px_S = I_ERI_K3x3yz_S_Px_S+ABY*I_ERI_I3x2yz_S_Px_S;
  Double I_ERI_I3xy2z_Py_Px_S = I_ERI_K3x2y2z_S_Px_S+ABY*I_ERI_I3xy2z_S_Px_S;
  Double I_ERI_I3x3z_Py_Px_S = I_ERI_K3xy3z_S_Px_S+ABY*I_ERI_I3x3z_S_Px_S;
  Double I_ERI_I2x4y_Py_Px_S = I_ERI_K2x5y_S_Px_S+ABY*I_ERI_I2x4y_S_Px_S;
  Double I_ERI_I2x3yz_Py_Px_S = I_ERI_K2x4yz_S_Px_S+ABY*I_ERI_I2x3yz_S_Px_S;
  Double I_ERI_I2x2y2z_Py_Px_S = I_ERI_K2x3y2z_S_Px_S+ABY*I_ERI_I2x2y2z_S_Px_S;
  Double I_ERI_I2xy3z_Py_Px_S = I_ERI_K2x2y3z_S_Px_S+ABY*I_ERI_I2xy3z_S_Px_S;
  Double I_ERI_I2x4z_Py_Px_S = I_ERI_K2xy4z_S_Px_S+ABY*I_ERI_I2x4z_S_Px_S;
  Double I_ERI_Ix5y_Py_Px_S = I_ERI_Kx6y_S_Px_S+ABY*I_ERI_Ix5y_S_Px_S;
  Double I_ERI_Ix4yz_Py_Px_S = I_ERI_Kx5yz_S_Px_S+ABY*I_ERI_Ix4yz_S_Px_S;
  Double I_ERI_Ix3y2z_Py_Px_S = I_ERI_Kx4y2z_S_Px_S+ABY*I_ERI_Ix3y2z_S_Px_S;
  Double I_ERI_Ix2y3z_Py_Px_S = I_ERI_Kx3y3z_S_Px_S+ABY*I_ERI_Ix2y3z_S_Px_S;
  Double I_ERI_Ixy4z_Py_Px_S = I_ERI_Kx2y4z_S_Px_S+ABY*I_ERI_Ixy4z_S_Px_S;
  Double I_ERI_Ix5z_Py_Px_S = I_ERI_Kxy5z_S_Px_S+ABY*I_ERI_Ix5z_S_Px_S;
  Double I_ERI_I6y_Py_Px_S = I_ERI_K7y_S_Px_S+ABY*I_ERI_I6y_S_Px_S;
  Double I_ERI_I5yz_Py_Px_S = I_ERI_K6yz_S_Px_S+ABY*I_ERI_I5yz_S_Px_S;
  Double I_ERI_I4y2z_Py_Px_S = I_ERI_K5y2z_S_Px_S+ABY*I_ERI_I4y2z_S_Px_S;
  Double I_ERI_I3y3z_Py_Px_S = I_ERI_K4y3z_S_Px_S+ABY*I_ERI_I3y3z_S_Px_S;
  Double I_ERI_I2y4z_Py_Px_S = I_ERI_K3y4z_S_Px_S+ABY*I_ERI_I2y4z_S_Px_S;
  Double I_ERI_Iy5z_Py_Px_S = I_ERI_K2y5z_S_Px_S+ABY*I_ERI_Iy5z_S_Px_S;
  Double I_ERI_I6z_Py_Px_S = I_ERI_Ky6z_S_Px_S+ABY*I_ERI_I6z_S_Px_S;
  Double I_ERI_I5xy_Pz_Px_S = I_ERI_K5xyz_S_Px_S+ABZ*I_ERI_I5xy_S_Px_S;
  Double I_ERI_I5xz_Pz_Px_S = I_ERI_K5x2z_S_Px_S+ABZ*I_ERI_I5xz_S_Px_S;
  Double I_ERI_I4x2y_Pz_Px_S = I_ERI_K4x2yz_S_Px_S+ABZ*I_ERI_I4x2y_S_Px_S;
  Double I_ERI_I4xyz_Pz_Px_S = I_ERI_K4xy2z_S_Px_S+ABZ*I_ERI_I4xyz_S_Px_S;
  Double I_ERI_I4x2z_Pz_Px_S = I_ERI_K4x3z_S_Px_S+ABZ*I_ERI_I4x2z_S_Px_S;
  Double I_ERI_I3x3y_Pz_Px_S = I_ERI_K3x3yz_S_Px_S+ABZ*I_ERI_I3x3y_S_Px_S;
  Double I_ERI_I3x2yz_Pz_Px_S = I_ERI_K3x2y2z_S_Px_S+ABZ*I_ERI_I3x2yz_S_Px_S;
  Double I_ERI_I3xy2z_Pz_Px_S = I_ERI_K3xy3z_S_Px_S+ABZ*I_ERI_I3xy2z_S_Px_S;
  Double I_ERI_I3x3z_Pz_Px_S = I_ERI_K3x4z_S_Px_S+ABZ*I_ERI_I3x3z_S_Px_S;
  Double I_ERI_I2x4y_Pz_Px_S = I_ERI_K2x4yz_S_Px_S+ABZ*I_ERI_I2x4y_S_Px_S;
  Double I_ERI_I2x3yz_Pz_Px_S = I_ERI_K2x3y2z_S_Px_S+ABZ*I_ERI_I2x3yz_S_Px_S;
  Double I_ERI_I2x2y2z_Pz_Px_S = I_ERI_K2x2y3z_S_Px_S+ABZ*I_ERI_I2x2y2z_S_Px_S;
  Double I_ERI_I2xy3z_Pz_Px_S = I_ERI_K2xy4z_S_Px_S+ABZ*I_ERI_I2xy3z_S_Px_S;
  Double I_ERI_I2x4z_Pz_Px_S = I_ERI_K2x5z_S_Px_S+ABZ*I_ERI_I2x4z_S_Px_S;
  Double I_ERI_Ix5y_Pz_Px_S = I_ERI_Kx5yz_S_Px_S+ABZ*I_ERI_Ix5y_S_Px_S;
  Double I_ERI_Ix4yz_Pz_Px_S = I_ERI_Kx4y2z_S_Px_S+ABZ*I_ERI_Ix4yz_S_Px_S;
  Double I_ERI_Ix3y2z_Pz_Px_S = I_ERI_Kx3y3z_S_Px_S+ABZ*I_ERI_Ix3y2z_S_Px_S;
  Double I_ERI_Ix2y3z_Pz_Px_S = I_ERI_Kx2y4z_S_Px_S+ABZ*I_ERI_Ix2y3z_S_Px_S;
  Double I_ERI_Ixy4z_Pz_Px_S = I_ERI_Kxy5z_S_Px_S+ABZ*I_ERI_Ixy4z_S_Px_S;
  Double I_ERI_Ix5z_Pz_Px_S = I_ERI_Kx6z_S_Px_S+ABZ*I_ERI_Ix5z_S_Px_S;
  Double I_ERI_I5yz_Pz_Px_S = I_ERI_K5y2z_S_Px_S+ABZ*I_ERI_I5yz_S_Px_S;
  Double I_ERI_I4y2z_Pz_Px_S = I_ERI_K4y3z_S_Px_S+ABZ*I_ERI_I4y2z_S_Px_S;
  Double I_ERI_I3y3z_Pz_Px_S = I_ERI_K3y4z_S_Px_S+ABZ*I_ERI_I3y3z_S_Px_S;
  Double I_ERI_I2y4z_Pz_Px_S = I_ERI_K2y5z_S_Px_S+ABZ*I_ERI_I2y4z_S_Px_S;
  Double I_ERI_Iy5z_Pz_Px_S = I_ERI_Ky6z_S_Px_S+ABZ*I_ERI_Iy5z_S_Px_S;
  Double I_ERI_I6z_Pz_Px_S = I_ERI_K7z_S_Px_S+ABZ*I_ERI_I6z_S_Px_S;
  Double I_ERI_I6x_Px_Py_S = I_ERI_K7x_S_Py_S+ABX*I_ERI_I6x_S_Py_S;
  Double I_ERI_I5xy_Px_Py_S = I_ERI_K6xy_S_Py_S+ABX*I_ERI_I5xy_S_Py_S;
  Double I_ERI_I5xz_Px_Py_S = I_ERI_K6xz_S_Py_S+ABX*I_ERI_I5xz_S_Py_S;
  Double I_ERI_I4x2y_Px_Py_S = I_ERI_K5x2y_S_Py_S+ABX*I_ERI_I4x2y_S_Py_S;
  Double I_ERI_I4xyz_Px_Py_S = I_ERI_K5xyz_S_Py_S+ABX*I_ERI_I4xyz_S_Py_S;
  Double I_ERI_I4x2z_Px_Py_S = I_ERI_K5x2z_S_Py_S+ABX*I_ERI_I4x2z_S_Py_S;
  Double I_ERI_I3x3y_Px_Py_S = I_ERI_K4x3y_S_Py_S+ABX*I_ERI_I3x3y_S_Py_S;
  Double I_ERI_I3x2yz_Px_Py_S = I_ERI_K4x2yz_S_Py_S+ABX*I_ERI_I3x2yz_S_Py_S;
  Double I_ERI_I3xy2z_Px_Py_S = I_ERI_K4xy2z_S_Py_S+ABX*I_ERI_I3xy2z_S_Py_S;
  Double I_ERI_I3x3z_Px_Py_S = I_ERI_K4x3z_S_Py_S+ABX*I_ERI_I3x3z_S_Py_S;
  Double I_ERI_I2x4y_Px_Py_S = I_ERI_K3x4y_S_Py_S+ABX*I_ERI_I2x4y_S_Py_S;
  Double I_ERI_I2x3yz_Px_Py_S = I_ERI_K3x3yz_S_Py_S+ABX*I_ERI_I2x3yz_S_Py_S;
  Double I_ERI_I2x2y2z_Px_Py_S = I_ERI_K3x2y2z_S_Py_S+ABX*I_ERI_I2x2y2z_S_Py_S;
  Double I_ERI_I2xy3z_Px_Py_S = I_ERI_K3xy3z_S_Py_S+ABX*I_ERI_I2xy3z_S_Py_S;
  Double I_ERI_I2x4z_Px_Py_S = I_ERI_K3x4z_S_Py_S+ABX*I_ERI_I2x4z_S_Py_S;
  Double I_ERI_Ix5y_Px_Py_S = I_ERI_K2x5y_S_Py_S+ABX*I_ERI_Ix5y_S_Py_S;
  Double I_ERI_Ix4yz_Px_Py_S = I_ERI_K2x4yz_S_Py_S+ABX*I_ERI_Ix4yz_S_Py_S;
  Double I_ERI_Ix3y2z_Px_Py_S = I_ERI_K2x3y2z_S_Py_S+ABX*I_ERI_Ix3y2z_S_Py_S;
  Double I_ERI_Ix2y3z_Px_Py_S = I_ERI_K2x2y3z_S_Py_S+ABX*I_ERI_Ix2y3z_S_Py_S;
  Double I_ERI_Ixy4z_Px_Py_S = I_ERI_K2xy4z_S_Py_S+ABX*I_ERI_Ixy4z_S_Py_S;
  Double I_ERI_Ix5z_Px_Py_S = I_ERI_K2x5z_S_Py_S+ABX*I_ERI_Ix5z_S_Py_S;
  Double I_ERI_I6y_Px_Py_S = I_ERI_Kx6y_S_Py_S+ABX*I_ERI_I6y_S_Py_S;
  Double I_ERI_I5yz_Px_Py_S = I_ERI_Kx5yz_S_Py_S+ABX*I_ERI_I5yz_S_Py_S;
  Double I_ERI_I4y2z_Px_Py_S = I_ERI_Kx4y2z_S_Py_S+ABX*I_ERI_I4y2z_S_Py_S;
  Double I_ERI_I3y3z_Px_Py_S = I_ERI_Kx3y3z_S_Py_S+ABX*I_ERI_I3y3z_S_Py_S;
  Double I_ERI_I2y4z_Px_Py_S = I_ERI_Kx2y4z_S_Py_S+ABX*I_ERI_I2y4z_S_Py_S;
  Double I_ERI_Iy5z_Px_Py_S = I_ERI_Kxy5z_S_Py_S+ABX*I_ERI_Iy5z_S_Py_S;
  Double I_ERI_I6z_Px_Py_S = I_ERI_Kx6z_S_Py_S+ABX*I_ERI_I6z_S_Py_S;
  Double I_ERI_I5xy_Py_Py_S = I_ERI_K5x2y_S_Py_S+ABY*I_ERI_I5xy_S_Py_S;
  Double I_ERI_I5xz_Py_Py_S = I_ERI_K5xyz_S_Py_S+ABY*I_ERI_I5xz_S_Py_S;
  Double I_ERI_I4x2y_Py_Py_S = I_ERI_K4x3y_S_Py_S+ABY*I_ERI_I4x2y_S_Py_S;
  Double I_ERI_I4xyz_Py_Py_S = I_ERI_K4x2yz_S_Py_S+ABY*I_ERI_I4xyz_S_Py_S;
  Double I_ERI_I4x2z_Py_Py_S = I_ERI_K4xy2z_S_Py_S+ABY*I_ERI_I4x2z_S_Py_S;
  Double I_ERI_I3x3y_Py_Py_S = I_ERI_K3x4y_S_Py_S+ABY*I_ERI_I3x3y_S_Py_S;
  Double I_ERI_I3x2yz_Py_Py_S = I_ERI_K3x3yz_S_Py_S+ABY*I_ERI_I3x2yz_S_Py_S;
  Double I_ERI_I3xy2z_Py_Py_S = I_ERI_K3x2y2z_S_Py_S+ABY*I_ERI_I3xy2z_S_Py_S;
  Double I_ERI_I3x3z_Py_Py_S = I_ERI_K3xy3z_S_Py_S+ABY*I_ERI_I3x3z_S_Py_S;
  Double I_ERI_I2x4y_Py_Py_S = I_ERI_K2x5y_S_Py_S+ABY*I_ERI_I2x4y_S_Py_S;
  Double I_ERI_I2x3yz_Py_Py_S = I_ERI_K2x4yz_S_Py_S+ABY*I_ERI_I2x3yz_S_Py_S;
  Double I_ERI_I2x2y2z_Py_Py_S = I_ERI_K2x3y2z_S_Py_S+ABY*I_ERI_I2x2y2z_S_Py_S;
  Double I_ERI_I2xy3z_Py_Py_S = I_ERI_K2x2y3z_S_Py_S+ABY*I_ERI_I2xy3z_S_Py_S;
  Double I_ERI_I2x4z_Py_Py_S = I_ERI_K2xy4z_S_Py_S+ABY*I_ERI_I2x4z_S_Py_S;
  Double I_ERI_Ix5y_Py_Py_S = I_ERI_Kx6y_S_Py_S+ABY*I_ERI_Ix5y_S_Py_S;
  Double I_ERI_Ix4yz_Py_Py_S = I_ERI_Kx5yz_S_Py_S+ABY*I_ERI_Ix4yz_S_Py_S;
  Double I_ERI_Ix3y2z_Py_Py_S = I_ERI_Kx4y2z_S_Py_S+ABY*I_ERI_Ix3y2z_S_Py_S;
  Double I_ERI_Ix2y3z_Py_Py_S = I_ERI_Kx3y3z_S_Py_S+ABY*I_ERI_Ix2y3z_S_Py_S;
  Double I_ERI_Ixy4z_Py_Py_S = I_ERI_Kx2y4z_S_Py_S+ABY*I_ERI_Ixy4z_S_Py_S;
  Double I_ERI_Ix5z_Py_Py_S = I_ERI_Kxy5z_S_Py_S+ABY*I_ERI_Ix5z_S_Py_S;
  Double I_ERI_I6y_Py_Py_S = I_ERI_K7y_S_Py_S+ABY*I_ERI_I6y_S_Py_S;
  Double I_ERI_I5yz_Py_Py_S = I_ERI_K6yz_S_Py_S+ABY*I_ERI_I5yz_S_Py_S;
  Double I_ERI_I4y2z_Py_Py_S = I_ERI_K5y2z_S_Py_S+ABY*I_ERI_I4y2z_S_Py_S;
  Double I_ERI_I3y3z_Py_Py_S = I_ERI_K4y3z_S_Py_S+ABY*I_ERI_I3y3z_S_Py_S;
  Double I_ERI_I2y4z_Py_Py_S = I_ERI_K3y4z_S_Py_S+ABY*I_ERI_I2y4z_S_Py_S;
  Double I_ERI_Iy5z_Py_Py_S = I_ERI_K2y5z_S_Py_S+ABY*I_ERI_Iy5z_S_Py_S;
  Double I_ERI_I6z_Py_Py_S = I_ERI_Ky6z_S_Py_S+ABY*I_ERI_I6z_S_Py_S;
  Double I_ERI_I5xy_Pz_Py_S = I_ERI_K5xyz_S_Py_S+ABZ*I_ERI_I5xy_S_Py_S;
  Double I_ERI_I5xz_Pz_Py_S = I_ERI_K5x2z_S_Py_S+ABZ*I_ERI_I5xz_S_Py_S;
  Double I_ERI_I4x2y_Pz_Py_S = I_ERI_K4x2yz_S_Py_S+ABZ*I_ERI_I4x2y_S_Py_S;
  Double I_ERI_I4xyz_Pz_Py_S = I_ERI_K4xy2z_S_Py_S+ABZ*I_ERI_I4xyz_S_Py_S;
  Double I_ERI_I4x2z_Pz_Py_S = I_ERI_K4x3z_S_Py_S+ABZ*I_ERI_I4x2z_S_Py_S;
  Double I_ERI_I3x3y_Pz_Py_S = I_ERI_K3x3yz_S_Py_S+ABZ*I_ERI_I3x3y_S_Py_S;
  Double I_ERI_I3x2yz_Pz_Py_S = I_ERI_K3x2y2z_S_Py_S+ABZ*I_ERI_I3x2yz_S_Py_S;
  Double I_ERI_I3xy2z_Pz_Py_S = I_ERI_K3xy3z_S_Py_S+ABZ*I_ERI_I3xy2z_S_Py_S;
  Double I_ERI_I3x3z_Pz_Py_S = I_ERI_K3x4z_S_Py_S+ABZ*I_ERI_I3x3z_S_Py_S;
  Double I_ERI_I2x4y_Pz_Py_S = I_ERI_K2x4yz_S_Py_S+ABZ*I_ERI_I2x4y_S_Py_S;
  Double I_ERI_I2x3yz_Pz_Py_S = I_ERI_K2x3y2z_S_Py_S+ABZ*I_ERI_I2x3yz_S_Py_S;
  Double I_ERI_I2x2y2z_Pz_Py_S = I_ERI_K2x2y3z_S_Py_S+ABZ*I_ERI_I2x2y2z_S_Py_S;
  Double I_ERI_I2xy3z_Pz_Py_S = I_ERI_K2xy4z_S_Py_S+ABZ*I_ERI_I2xy3z_S_Py_S;
  Double I_ERI_I2x4z_Pz_Py_S = I_ERI_K2x5z_S_Py_S+ABZ*I_ERI_I2x4z_S_Py_S;
  Double I_ERI_Ix5y_Pz_Py_S = I_ERI_Kx5yz_S_Py_S+ABZ*I_ERI_Ix5y_S_Py_S;
  Double I_ERI_Ix4yz_Pz_Py_S = I_ERI_Kx4y2z_S_Py_S+ABZ*I_ERI_Ix4yz_S_Py_S;
  Double I_ERI_Ix3y2z_Pz_Py_S = I_ERI_Kx3y3z_S_Py_S+ABZ*I_ERI_Ix3y2z_S_Py_S;
  Double I_ERI_Ix2y3z_Pz_Py_S = I_ERI_Kx2y4z_S_Py_S+ABZ*I_ERI_Ix2y3z_S_Py_S;
  Double I_ERI_Ixy4z_Pz_Py_S = I_ERI_Kxy5z_S_Py_S+ABZ*I_ERI_Ixy4z_S_Py_S;
  Double I_ERI_Ix5z_Pz_Py_S = I_ERI_Kx6z_S_Py_S+ABZ*I_ERI_Ix5z_S_Py_S;
  Double I_ERI_I5yz_Pz_Py_S = I_ERI_K5y2z_S_Py_S+ABZ*I_ERI_I5yz_S_Py_S;
  Double I_ERI_I4y2z_Pz_Py_S = I_ERI_K4y3z_S_Py_S+ABZ*I_ERI_I4y2z_S_Py_S;
  Double I_ERI_I3y3z_Pz_Py_S = I_ERI_K3y4z_S_Py_S+ABZ*I_ERI_I3y3z_S_Py_S;
  Double I_ERI_I2y4z_Pz_Py_S = I_ERI_K2y5z_S_Py_S+ABZ*I_ERI_I2y4z_S_Py_S;
  Double I_ERI_Iy5z_Pz_Py_S = I_ERI_Ky6z_S_Py_S+ABZ*I_ERI_Iy5z_S_Py_S;
  Double I_ERI_I6z_Pz_Py_S = I_ERI_K7z_S_Py_S+ABZ*I_ERI_I6z_S_Py_S;
  Double I_ERI_I6x_Px_Pz_S = I_ERI_K7x_S_Pz_S+ABX*I_ERI_I6x_S_Pz_S;
  Double I_ERI_I5xy_Px_Pz_S = I_ERI_K6xy_S_Pz_S+ABX*I_ERI_I5xy_S_Pz_S;
  Double I_ERI_I5xz_Px_Pz_S = I_ERI_K6xz_S_Pz_S+ABX*I_ERI_I5xz_S_Pz_S;
  Double I_ERI_I4x2y_Px_Pz_S = I_ERI_K5x2y_S_Pz_S+ABX*I_ERI_I4x2y_S_Pz_S;
  Double I_ERI_I4xyz_Px_Pz_S = I_ERI_K5xyz_S_Pz_S+ABX*I_ERI_I4xyz_S_Pz_S;
  Double I_ERI_I4x2z_Px_Pz_S = I_ERI_K5x2z_S_Pz_S+ABX*I_ERI_I4x2z_S_Pz_S;
  Double I_ERI_I3x3y_Px_Pz_S = I_ERI_K4x3y_S_Pz_S+ABX*I_ERI_I3x3y_S_Pz_S;
  Double I_ERI_I3x2yz_Px_Pz_S = I_ERI_K4x2yz_S_Pz_S+ABX*I_ERI_I3x2yz_S_Pz_S;
  Double I_ERI_I3xy2z_Px_Pz_S = I_ERI_K4xy2z_S_Pz_S+ABX*I_ERI_I3xy2z_S_Pz_S;
  Double I_ERI_I3x3z_Px_Pz_S = I_ERI_K4x3z_S_Pz_S+ABX*I_ERI_I3x3z_S_Pz_S;
  Double I_ERI_I2x4y_Px_Pz_S = I_ERI_K3x4y_S_Pz_S+ABX*I_ERI_I2x4y_S_Pz_S;
  Double I_ERI_I2x3yz_Px_Pz_S = I_ERI_K3x3yz_S_Pz_S+ABX*I_ERI_I2x3yz_S_Pz_S;
  Double I_ERI_I2x2y2z_Px_Pz_S = I_ERI_K3x2y2z_S_Pz_S+ABX*I_ERI_I2x2y2z_S_Pz_S;
  Double I_ERI_I2xy3z_Px_Pz_S = I_ERI_K3xy3z_S_Pz_S+ABX*I_ERI_I2xy3z_S_Pz_S;
  Double I_ERI_I2x4z_Px_Pz_S = I_ERI_K3x4z_S_Pz_S+ABX*I_ERI_I2x4z_S_Pz_S;
  Double I_ERI_Ix5y_Px_Pz_S = I_ERI_K2x5y_S_Pz_S+ABX*I_ERI_Ix5y_S_Pz_S;
  Double I_ERI_Ix4yz_Px_Pz_S = I_ERI_K2x4yz_S_Pz_S+ABX*I_ERI_Ix4yz_S_Pz_S;
  Double I_ERI_Ix3y2z_Px_Pz_S = I_ERI_K2x3y2z_S_Pz_S+ABX*I_ERI_Ix3y2z_S_Pz_S;
  Double I_ERI_Ix2y3z_Px_Pz_S = I_ERI_K2x2y3z_S_Pz_S+ABX*I_ERI_Ix2y3z_S_Pz_S;
  Double I_ERI_Ixy4z_Px_Pz_S = I_ERI_K2xy4z_S_Pz_S+ABX*I_ERI_Ixy4z_S_Pz_S;
  Double I_ERI_Ix5z_Px_Pz_S = I_ERI_K2x5z_S_Pz_S+ABX*I_ERI_Ix5z_S_Pz_S;
  Double I_ERI_I6y_Px_Pz_S = I_ERI_Kx6y_S_Pz_S+ABX*I_ERI_I6y_S_Pz_S;
  Double I_ERI_I5yz_Px_Pz_S = I_ERI_Kx5yz_S_Pz_S+ABX*I_ERI_I5yz_S_Pz_S;
  Double I_ERI_I4y2z_Px_Pz_S = I_ERI_Kx4y2z_S_Pz_S+ABX*I_ERI_I4y2z_S_Pz_S;
  Double I_ERI_I3y3z_Px_Pz_S = I_ERI_Kx3y3z_S_Pz_S+ABX*I_ERI_I3y3z_S_Pz_S;
  Double I_ERI_I2y4z_Px_Pz_S = I_ERI_Kx2y4z_S_Pz_S+ABX*I_ERI_I2y4z_S_Pz_S;
  Double I_ERI_Iy5z_Px_Pz_S = I_ERI_Kxy5z_S_Pz_S+ABX*I_ERI_Iy5z_S_Pz_S;
  Double I_ERI_I6z_Px_Pz_S = I_ERI_Kx6z_S_Pz_S+ABX*I_ERI_I6z_S_Pz_S;
  Double I_ERI_I5xy_Py_Pz_S = I_ERI_K5x2y_S_Pz_S+ABY*I_ERI_I5xy_S_Pz_S;
  Double I_ERI_I5xz_Py_Pz_S = I_ERI_K5xyz_S_Pz_S+ABY*I_ERI_I5xz_S_Pz_S;
  Double I_ERI_I4x2y_Py_Pz_S = I_ERI_K4x3y_S_Pz_S+ABY*I_ERI_I4x2y_S_Pz_S;
  Double I_ERI_I4xyz_Py_Pz_S = I_ERI_K4x2yz_S_Pz_S+ABY*I_ERI_I4xyz_S_Pz_S;
  Double I_ERI_I4x2z_Py_Pz_S = I_ERI_K4xy2z_S_Pz_S+ABY*I_ERI_I4x2z_S_Pz_S;
  Double I_ERI_I3x3y_Py_Pz_S = I_ERI_K3x4y_S_Pz_S+ABY*I_ERI_I3x3y_S_Pz_S;
  Double I_ERI_I3x2yz_Py_Pz_S = I_ERI_K3x3yz_S_Pz_S+ABY*I_ERI_I3x2yz_S_Pz_S;
  Double I_ERI_I3xy2z_Py_Pz_S = I_ERI_K3x2y2z_S_Pz_S+ABY*I_ERI_I3xy2z_S_Pz_S;
  Double I_ERI_I3x3z_Py_Pz_S = I_ERI_K3xy3z_S_Pz_S+ABY*I_ERI_I3x3z_S_Pz_S;
  Double I_ERI_I2x4y_Py_Pz_S = I_ERI_K2x5y_S_Pz_S+ABY*I_ERI_I2x4y_S_Pz_S;
  Double I_ERI_I2x3yz_Py_Pz_S = I_ERI_K2x4yz_S_Pz_S+ABY*I_ERI_I2x3yz_S_Pz_S;
  Double I_ERI_I2x2y2z_Py_Pz_S = I_ERI_K2x3y2z_S_Pz_S+ABY*I_ERI_I2x2y2z_S_Pz_S;
  Double I_ERI_I2xy3z_Py_Pz_S = I_ERI_K2x2y3z_S_Pz_S+ABY*I_ERI_I2xy3z_S_Pz_S;
  Double I_ERI_I2x4z_Py_Pz_S = I_ERI_K2xy4z_S_Pz_S+ABY*I_ERI_I2x4z_S_Pz_S;
  Double I_ERI_Ix5y_Py_Pz_S = I_ERI_Kx6y_S_Pz_S+ABY*I_ERI_Ix5y_S_Pz_S;
  Double I_ERI_Ix4yz_Py_Pz_S = I_ERI_Kx5yz_S_Pz_S+ABY*I_ERI_Ix4yz_S_Pz_S;
  Double I_ERI_Ix3y2z_Py_Pz_S = I_ERI_Kx4y2z_S_Pz_S+ABY*I_ERI_Ix3y2z_S_Pz_S;
  Double I_ERI_Ix2y3z_Py_Pz_S = I_ERI_Kx3y3z_S_Pz_S+ABY*I_ERI_Ix2y3z_S_Pz_S;
  Double I_ERI_Ixy4z_Py_Pz_S = I_ERI_Kx2y4z_S_Pz_S+ABY*I_ERI_Ixy4z_S_Pz_S;
  Double I_ERI_Ix5z_Py_Pz_S = I_ERI_Kxy5z_S_Pz_S+ABY*I_ERI_Ix5z_S_Pz_S;
  Double I_ERI_I6y_Py_Pz_S = I_ERI_K7y_S_Pz_S+ABY*I_ERI_I6y_S_Pz_S;
  Double I_ERI_I5yz_Py_Pz_S = I_ERI_K6yz_S_Pz_S+ABY*I_ERI_I5yz_S_Pz_S;
  Double I_ERI_I4y2z_Py_Pz_S = I_ERI_K5y2z_S_Pz_S+ABY*I_ERI_I4y2z_S_Pz_S;
  Double I_ERI_I3y3z_Py_Pz_S = I_ERI_K4y3z_S_Pz_S+ABY*I_ERI_I3y3z_S_Pz_S;
  Double I_ERI_I2y4z_Py_Pz_S = I_ERI_K3y4z_S_Pz_S+ABY*I_ERI_I2y4z_S_Pz_S;
  Double I_ERI_Iy5z_Py_Pz_S = I_ERI_K2y5z_S_Pz_S+ABY*I_ERI_Iy5z_S_Pz_S;
  Double I_ERI_I6z_Py_Pz_S = I_ERI_Ky6z_S_Pz_S+ABY*I_ERI_I6z_S_Pz_S;
  Double I_ERI_I5xy_Pz_Pz_S = I_ERI_K5xyz_S_Pz_S+ABZ*I_ERI_I5xy_S_Pz_S;
  Double I_ERI_I5xz_Pz_Pz_S = I_ERI_K5x2z_S_Pz_S+ABZ*I_ERI_I5xz_S_Pz_S;
  Double I_ERI_I4x2y_Pz_Pz_S = I_ERI_K4x2yz_S_Pz_S+ABZ*I_ERI_I4x2y_S_Pz_S;
  Double I_ERI_I4xyz_Pz_Pz_S = I_ERI_K4xy2z_S_Pz_S+ABZ*I_ERI_I4xyz_S_Pz_S;
  Double I_ERI_I4x2z_Pz_Pz_S = I_ERI_K4x3z_S_Pz_S+ABZ*I_ERI_I4x2z_S_Pz_S;
  Double I_ERI_I3x3y_Pz_Pz_S = I_ERI_K3x3yz_S_Pz_S+ABZ*I_ERI_I3x3y_S_Pz_S;
  Double I_ERI_I3x2yz_Pz_Pz_S = I_ERI_K3x2y2z_S_Pz_S+ABZ*I_ERI_I3x2yz_S_Pz_S;
  Double I_ERI_I3xy2z_Pz_Pz_S = I_ERI_K3xy3z_S_Pz_S+ABZ*I_ERI_I3xy2z_S_Pz_S;
  Double I_ERI_I3x3z_Pz_Pz_S = I_ERI_K3x4z_S_Pz_S+ABZ*I_ERI_I3x3z_S_Pz_S;
  Double I_ERI_I2x4y_Pz_Pz_S = I_ERI_K2x4yz_S_Pz_S+ABZ*I_ERI_I2x4y_S_Pz_S;
  Double I_ERI_I2x3yz_Pz_Pz_S = I_ERI_K2x3y2z_S_Pz_S+ABZ*I_ERI_I2x3yz_S_Pz_S;
  Double I_ERI_I2x2y2z_Pz_Pz_S = I_ERI_K2x2y3z_S_Pz_S+ABZ*I_ERI_I2x2y2z_S_Pz_S;
  Double I_ERI_I2xy3z_Pz_Pz_S = I_ERI_K2xy4z_S_Pz_S+ABZ*I_ERI_I2xy3z_S_Pz_S;
  Double I_ERI_I2x4z_Pz_Pz_S = I_ERI_K2x5z_S_Pz_S+ABZ*I_ERI_I2x4z_S_Pz_S;
  Double I_ERI_Ix5y_Pz_Pz_S = I_ERI_Kx5yz_S_Pz_S+ABZ*I_ERI_Ix5y_S_Pz_S;
  Double I_ERI_Ix4yz_Pz_Pz_S = I_ERI_Kx4y2z_S_Pz_S+ABZ*I_ERI_Ix4yz_S_Pz_S;
  Double I_ERI_Ix3y2z_Pz_Pz_S = I_ERI_Kx3y3z_S_Pz_S+ABZ*I_ERI_Ix3y2z_S_Pz_S;
  Double I_ERI_Ix2y3z_Pz_Pz_S = I_ERI_Kx2y4z_S_Pz_S+ABZ*I_ERI_Ix2y3z_S_Pz_S;
  Double I_ERI_Ixy4z_Pz_Pz_S = I_ERI_Kxy5z_S_Pz_S+ABZ*I_ERI_Ixy4z_S_Pz_S;
  Double I_ERI_Ix5z_Pz_Pz_S = I_ERI_Kx6z_S_Pz_S+ABZ*I_ERI_Ix5z_S_Pz_S;
  Double I_ERI_I5yz_Pz_Pz_S = I_ERI_K5y2z_S_Pz_S+ABZ*I_ERI_I5yz_S_Pz_S;
  Double I_ERI_I4y2z_Pz_Pz_S = I_ERI_K4y3z_S_Pz_S+ABZ*I_ERI_I4y2z_S_Pz_S;
  Double I_ERI_I3y3z_Pz_Pz_S = I_ERI_K3y4z_S_Pz_S+ABZ*I_ERI_I3y3z_S_Pz_S;
  Double I_ERI_I2y4z_Pz_Pz_S = I_ERI_K2y5z_S_Pz_S+ABZ*I_ERI_I2y4z_S_Pz_S;
  Double I_ERI_Iy5z_Pz_Pz_S = I_ERI_Ky6z_S_Pz_S+ABZ*I_ERI_Iy5z_S_Pz_S;
  Double I_ERI_I6z_Pz_Pz_S = I_ERI_K7z_S_Pz_S+ABZ*I_ERI_I6z_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_H_D_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 189 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_P_P_S
   * RHS shell quartet name: SQ_ERI_H_P_P_S
   ************************************************************/
  Double I_ERI_H5x_D2x_Px_S = I_ERI_I6x_Px_Px_S+ABX*I_ERI_H5x_Px_Px_S;
  Double I_ERI_H4xy_D2x_Px_S = I_ERI_I5xy_Px_Px_S+ABX*I_ERI_H4xy_Px_Px_S;
  Double I_ERI_H4xz_D2x_Px_S = I_ERI_I5xz_Px_Px_S+ABX*I_ERI_H4xz_Px_Px_S;
  Double I_ERI_H3x2y_D2x_Px_S = I_ERI_I4x2y_Px_Px_S+ABX*I_ERI_H3x2y_Px_Px_S;
  Double I_ERI_H3xyz_D2x_Px_S = I_ERI_I4xyz_Px_Px_S+ABX*I_ERI_H3xyz_Px_Px_S;
  Double I_ERI_H3x2z_D2x_Px_S = I_ERI_I4x2z_Px_Px_S+ABX*I_ERI_H3x2z_Px_Px_S;
  Double I_ERI_H2x3y_D2x_Px_S = I_ERI_I3x3y_Px_Px_S+ABX*I_ERI_H2x3y_Px_Px_S;
  Double I_ERI_H2x2yz_D2x_Px_S = I_ERI_I3x2yz_Px_Px_S+ABX*I_ERI_H2x2yz_Px_Px_S;
  Double I_ERI_H2xy2z_D2x_Px_S = I_ERI_I3xy2z_Px_Px_S+ABX*I_ERI_H2xy2z_Px_Px_S;
  Double I_ERI_H2x3z_D2x_Px_S = I_ERI_I3x3z_Px_Px_S+ABX*I_ERI_H2x3z_Px_Px_S;
  Double I_ERI_Hx4y_D2x_Px_S = I_ERI_I2x4y_Px_Px_S+ABX*I_ERI_Hx4y_Px_Px_S;
  Double I_ERI_Hx3yz_D2x_Px_S = I_ERI_I2x3yz_Px_Px_S+ABX*I_ERI_Hx3yz_Px_Px_S;
  Double I_ERI_Hx2y2z_D2x_Px_S = I_ERI_I2x2y2z_Px_Px_S+ABX*I_ERI_Hx2y2z_Px_Px_S;
  Double I_ERI_Hxy3z_D2x_Px_S = I_ERI_I2xy3z_Px_Px_S+ABX*I_ERI_Hxy3z_Px_Px_S;
  Double I_ERI_Hx4z_D2x_Px_S = I_ERI_I2x4z_Px_Px_S+ABX*I_ERI_Hx4z_Px_Px_S;
  Double I_ERI_H5y_D2x_Px_S = I_ERI_Ix5y_Px_Px_S+ABX*I_ERI_H5y_Px_Px_S;
  Double I_ERI_H4yz_D2x_Px_S = I_ERI_Ix4yz_Px_Px_S+ABX*I_ERI_H4yz_Px_Px_S;
  Double I_ERI_H3y2z_D2x_Px_S = I_ERI_Ix3y2z_Px_Px_S+ABX*I_ERI_H3y2z_Px_Px_S;
  Double I_ERI_H2y3z_D2x_Px_S = I_ERI_Ix2y3z_Px_Px_S+ABX*I_ERI_H2y3z_Px_Px_S;
  Double I_ERI_Hy4z_D2x_Px_S = I_ERI_Ixy4z_Px_Px_S+ABX*I_ERI_Hy4z_Px_Px_S;
  Double I_ERI_H5z_D2x_Px_S = I_ERI_Ix5z_Px_Px_S+ABX*I_ERI_H5z_Px_Px_S;
  Double I_ERI_H5x_D2y_Px_S = I_ERI_I5xy_Py_Px_S+ABY*I_ERI_H5x_Py_Px_S;
  Double I_ERI_H4xy_D2y_Px_S = I_ERI_I4x2y_Py_Px_S+ABY*I_ERI_H4xy_Py_Px_S;
  Double I_ERI_H4xz_D2y_Px_S = I_ERI_I4xyz_Py_Px_S+ABY*I_ERI_H4xz_Py_Px_S;
  Double I_ERI_H3x2y_D2y_Px_S = I_ERI_I3x3y_Py_Px_S+ABY*I_ERI_H3x2y_Py_Px_S;
  Double I_ERI_H3xyz_D2y_Px_S = I_ERI_I3x2yz_Py_Px_S+ABY*I_ERI_H3xyz_Py_Px_S;
  Double I_ERI_H3x2z_D2y_Px_S = I_ERI_I3xy2z_Py_Px_S+ABY*I_ERI_H3x2z_Py_Px_S;
  Double I_ERI_H2x3y_D2y_Px_S = I_ERI_I2x4y_Py_Px_S+ABY*I_ERI_H2x3y_Py_Px_S;
  Double I_ERI_H2x2yz_D2y_Px_S = I_ERI_I2x3yz_Py_Px_S+ABY*I_ERI_H2x2yz_Py_Px_S;
  Double I_ERI_H2xy2z_D2y_Px_S = I_ERI_I2x2y2z_Py_Px_S+ABY*I_ERI_H2xy2z_Py_Px_S;
  Double I_ERI_H2x3z_D2y_Px_S = I_ERI_I2xy3z_Py_Px_S+ABY*I_ERI_H2x3z_Py_Px_S;
  Double I_ERI_Hx4y_D2y_Px_S = I_ERI_Ix5y_Py_Px_S+ABY*I_ERI_Hx4y_Py_Px_S;
  Double I_ERI_Hx3yz_D2y_Px_S = I_ERI_Ix4yz_Py_Px_S+ABY*I_ERI_Hx3yz_Py_Px_S;
  Double I_ERI_Hx2y2z_D2y_Px_S = I_ERI_Ix3y2z_Py_Px_S+ABY*I_ERI_Hx2y2z_Py_Px_S;
  Double I_ERI_Hxy3z_D2y_Px_S = I_ERI_Ix2y3z_Py_Px_S+ABY*I_ERI_Hxy3z_Py_Px_S;
  Double I_ERI_Hx4z_D2y_Px_S = I_ERI_Ixy4z_Py_Px_S+ABY*I_ERI_Hx4z_Py_Px_S;
  Double I_ERI_H5y_D2y_Px_S = I_ERI_I6y_Py_Px_S+ABY*I_ERI_H5y_Py_Px_S;
  Double I_ERI_H4yz_D2y_Px_S = I_ERI_I5yz_Py_Px_S+ABY*I_ERI_H4yz_Py_Px_S;
  Double I_ERI_H3y2z_D2y_Px_S = I_ERI_I4y2z_Py_Px_S+ABY*I_ERI_H3y2z_Py_Px_S;
  Double I_ERI_H2y3z_D2y_Px_S = I_ERI_I3y3z_Py_Px_S+ABY*I_ERI_H2y3z_Py_Px_S;
  Double I_ERI_Hy4z_D2y_Px_S = I_ERI_I2y4z_Py_Px_S+ABY*I_ERI_Hy4z_Py_Px_S;
  Double I_ERI_H5z_D2y_Px_S = I_ERI_Iy5z_Py_Px_S+ABY*I_ERI_H5z_Py_Px_S;
  Double I_ERI_H5x_D2z_Px_S = I_ERI_I5xz_Pz_Px_S+ABZ*I_ERI_H5x_Pz_Px_S;
  Double I_ERI_H4xy_D2z_Px_S = I_ERI_I4xyz_Pz_Px_S+ABZ*I_ERI_H4xy_Pz_Px_S;
  Double I_ERI_H4xz_D2z_Px_S = I_ERI_I4x2z_Pz_Px_S+ABZ*I_ERI_H4xz_Pz_Px_S;
  Double I_ERI_H3x2y_D2z_Px_S = I_ERI_I3x2yz_Pz_Px_S+ABZ*I_ERI_H3x2y_Pz_Px_S;
  Double I_ERI_H3xyz_D2z_Px_S = I_ERI_I3xy2z_Pz_Px_S+ABZ*I_ERI_H3xyz_Pz_Px_S;
  Double I_ERI_H3x2z_D2z_Px_S = I_ERI_I3x3z_Pz_Px_S+ABZ*I_ERI_H3x2z_Pz_Px_S;
  Double I_ERI_H2x3y_D2z_Px_S = I_ERI_I2x3yz_Pz_Px_S+ABZ*I_ERI_H2x3y_Pz_Px_S;
  Double I_ERI_H2x2yz_D2z_Px_S = I_ERI_I2x2y2z_Pz_Px_S+ABZ*I_ERI_H2x2yz_Pz_Px_S;
  Double I_ERI_H2xy2z_D2z_Px_S = I_ERI_I2xy3z_Pz_Px_S+ABZ*I_ERI_H2xy2z_Pz_Px_S;
  Double I_ERI_H2x3z_D2z_Px_S = I_ERI_I2x4z_Pz_Px_S+ABZ*I_ERI_H2x3z_Pz_Px_S;
  Double I_ERI_Hx4y_D2z_Px_S = I_ERI_Ix4yz_Pz_Px_S+ABZ*I_ERI_Hx4y_Pz_Px_S;
  Double I_ERI_Hx3yz_D2z_Px_S = I_ERI_Ix3y2z_Pz_Px_S+ABZ*I_ERI_Hx3yz_Pz_Px_S;
  Double I_ERI_Hx2y2z_D2z_Px_S = I_ERI_Ix2y3z_Pz_Px_S+ABZ*I_ERI_Hx2y2z_Pz_Px_S;
  Double I_ERI_Hxy3z_D2z_Px_S = I_ERI_Ixy4z_Pz_Px_S+ABZ*I_ERI_Hxy3z_Pz_Px_S;
  Double I_ERI_Hx4z_D2z_Px_S = I_ERI_Ix5z_Pz_Px_S+ABZ*I_ERI_Hx4z_Pz_Px_S;
  Double I_ERI_H5y_D2z_Px_S = I_ERI_I5yz_Pz_Px_S+ABZ*I_ERI_H5y_Pz_Px_S;
  Double I_ERI_H4yz_D2z_Px_S = I_ERI_I4y2z_Pz_Px_S+ABZ*I_ERI_H4yz_Pz_Px_S;
  Double I_ERI_H3y2z_D2z_Px_S = I_ERI_I3y3z_Pz_Px_S+ABZ*I_ERI_H3y2z_Pz_Px_S;
  Double I_ERI_H2y3z_D2z_Px_S = I_ERI_I2y4z_Pz_Px_S+ABZ*I_ERI_H2y3z_Pz_Px_S;
  Double I_ERI_Hy4z_D2z_Px_S = I_ERI_Iy5z_Pz_Px_S+ABZ*I_ERI_Hy4z_Pz_Px_S;
  Double I_ERI_H5z_D2z_Px_S = I_ERI_I6z_Pz_Px_S+ABZ*I_ERI_H5z_Pz_Px_S;
  Double I_ERI_H5x_D2x_Py_S = I_ERI_I6x_Px_Py_S+ABX*I_ERI_H5x_Px_Py_S;
  Double I_ERI_H4xy_D2x_Py_S = I_ERI_I5xy_Px_Py_S+ABX*I_ERI_H4xy_Px_Py_S;
  Double I_ERI_H4xz_D2x_Py_S = I_ERI_I5xz_Px_Py_S+ABX*I_ERI_H4xz_Px_Py_S;
  Double I_ERI_H3x2y_D2x_Py_S = I_ERI_I4x2y_Px_Py_S+ABX*I_ERI_H3x2y_Px_Py_S;
  Double I_ERI_H3xyz_D2x_Py_S = I_ERI_I4xyz_Px_Py_S+ABX*I_ERI_H3xyz_Px_Py_S;
  Double I_ERI_H3x2z_D2x_Py_S = I_ERI_I4x2z_Px_Py_S+ABX*I_ERI_H3x2z_Px_Py_S;
  Double I_ERI_H2x3y_D2x_Py_S = I_ERI_I3x3y_Px_Py_S+ABX*I_ERI_H2x3y_Px_Py_S;
  Double I_ERI_H2x2yz_D2x_Py_S = I_ERI_I3x2yz_Px_Py_S+ABX*I_ERI_H2x2yz_Px_Py_S;
  Double I_ERI_H2xy2z_D2x_Py_S = I_ERI_I3xy2z_Px_Py_S+ABX*I_ERI_H2xy2z_Px_Py_S;
  Double I_ERI_H2x3z_D2x_Py_S = I_ERI_I3x3z_Px_Py_S+ABX*I_ERI_H2x3z_Px_Py_S;
  Double I_ERI_Hx4y_D2x_Py_S = I_ERI_I2x4y_Px_Py_S+ABX*I_ERI_Hx4y_Px_Py_S;
  Double I_ERI_Hx3yz_D2x_Py_S = I_ERI_I2x3yz_Px_Py_S+ABX*I_ERI_Hx3yz_Px_Py_S;
  Double I_ERI_Hx2y2z_D2x_Py_S = I_ERI_I2x2y2z_Px_Py_S+ABX*I_ERI_Hx2y2z_Px_Py_S;
  Double I_ERI_Hxy3z_D2x_Py_S = I_ERI_I2xy3z_Px_Py_S+ABX*I_ERI_Hxy3z_Px_Py_S;
  Double I_ERI_Hx4z_D2x_Py_S = I_ERI_I2x4z_Px_Py_S+ABX*I_ERI_Hx4z_Px_Py_S;
  Double I_ERI_H5y_D2x_Py_S = I_ERI_Ix5y_Px_Py_S+ABX*I_ERI_H5y_Px_Py_S;
  Double I_ERI_H4yz_D2x_Py_S = I_ERI_Ix4yz_Px_Py_S+ABX*I_ERI_H4yz_Px_Py_S;
  Double I_ERI_H3y2z_D2x_Py_S = I_ERI_Ix3y2z_Px_Py_S+ABX*I_ERI_H3y2z_Px_Py_S;
  Double I_ERI_H2y3z_D2x_Py_S = I_ERI_Ix2y3z_Px_Py_S+ABX*I_ERI_H2y3z_Px_Py_S;
  Double I_ERI_Hy4z_D2x_Py_S = I_ERI_Ixy4z_Px_Py_S+ABX*I_ERI_Hy4z_Px_Py_S;
  Double I_ERI_H5z_D2x_Py_S = I_ERI_Ix5z_Px_Py_S+ABX*I_ERI_H5z_Px_Py_S;
  Double I_ERI_H5x_D2y_Py_S = I_ERI_I5xy_Py_Py_S+ABY*I_ERI_H5x_Py_Py_S;
  Double I_ERI_H4xy_D2y_Py_S = I_ERI_I4x2y_Py_Py_S+ABY*I_ERI_H4xy_Py_Py_S;
  Double I_ERI_H4xz_D2y_Py_S = I_ERI_I4xyz_Py_Py_S+ABY*I_ERI_H4xz_Py_Py_S;
  Double I_ERI_H3x2y_D2y_Py_S = I_ERI_I3x3y_Py_Py_S+ABY*I_ERI_H3x2y_Py_Py_S;
  Double I_ERI_H3xyz_D2y_Py_S = I_ERI_I3x2yz_Py_Py_S+ABY*I_ERI_H3xyz_Py_Py_S;
  Double I_ERI_H3x2z_D2y_Py_S = I_ERI_I3xy2z_Py_Py_S+ABY*I_ERI_H3x2z_Py_Py_S;
  Double I_ERI_H2x3y_D2y_Py_S = I_ERI_I2x4y_Py_Py_S+ABY*I_ERI_H2x3y_Py_Py_S;
  Double I_ERI_H2x2yz_D2y_Py_S = I_ERI_I2x3yz_Py_Py_S+ABY*I_ERI_H2x2yz_Py_Py_S;
  Double I_ERI_H2xy2z_D2y_Py_S = I_ERI_I2x2y2z_Py_Py_S+ABY*I_ERI_H2xy2z_Py_Py_S;
  Double I_ERI_H2x3z_D2y_Py_S = I_ERI_I2xy3z_Py_Py_S+ABY*I_ERI_H2x3z_Py_Py_S;
  Double I_ERI_Hx4y_D2y_Py_S = I_ERI_Ix5y_Py_Py_S+ABY*I_ERI_Hx4y_Py_Py_S;
  Double I_ERI_Hx3yz_D2y_Py_S = I_ERI_Ix4yz_Py_Py_S+ABY*I_ERI_Hx3yz_Py_Py_S;
  Double I_ERI_Hx2y2z_D2y_Py_S = I_ERI_Ix3y2z_Py_Py_S+ABY*I_ERI_Hx2y2z_Py_Py_S;
  Double I_ERI_Hxy3z_D2y_Py_S = I_ERI_Ix2y3z_Py_Py_S+ABY*I_ERI_Hxy3z_Py_Py_S;
  Double I_ERI_Hx4z_D2y_Py_S = I_ERI_Ixy4z_Py_Py_S+ABY*I_ERI_Hx4z_Py_Py_S;
  Double I_ERI_H5y_D2y_Py_S = I_ERI_I6y_Py_Py_S+ABY*I_ERI_H5y_Py_Py_S;
  Double I_ERI_H4yz_D2y_Py_S = I_ERI_I5yz_Py_Py_S+ABY*I_ERI_H4yz_Py_Py_S;
  Double I_ERI_H3y2z_D2y_Py_S = I_ERI_I4y2z_Py_Py_S+ABY*I_ERI_H3y2z_Py_Py_S;
  Double I_ERI_H2y3z_D2y_Py_S = I_ERI_I3y3z_Py_Py_S+ABY*I_ERI_H2y3z_Py_Py_S;
  Double I_ERI_Hy4z_D2y_Py_S = I_ERI_I2y4z_Py_Py_S+ABY*I_ERI_Hy4z_Py_Py_S;
  Double I_ERI_H5z_D2y_Py_S = I_ERI_Iy5z_Py_Py_S+ABY*I_ERI_H5z_Py_Py_S;
  Double I_ERI_H5x_D2z_Py_S = I_ERI_I5xz_Pz_Py_S+ABZ*I_ERI_H5x_Pz_Py_S;
  Double I_ERI_H4xy_D2z_Py_S = I_ERI_I4xyz_Pz_Py_S+ABZ*I_ERI_H4xy_Pz_Py_S;
  Double I_ERI_H4xz_D2z_Py_S = I_ERI_I4x2z_Pz_Py_S+ABZ*I_ERI_H4xz_Pz_Py_S;
  Double I_ERI_H3x2y_D2z_Py_S = I_ERI_I3x2yz_Pz_Py_S+ABZ*I_ERI_H3x2y_Pz_Py_S;
  Double I_ERI_H3xyz_D2z_Py_S = I_ERI_I3xy2z_Pz_Py_S+ABZ*I_ERI_H3xyz_Pz_Py_S;
  Double I_ERI_H3x2z_D2z_Py_S = I_ERI_I3x3z_Pz_Py_S+ABZ*I_ERI_H3x2z_Pz_Py_S;
  Double I_ERI_H2x3y_D2z_Py_S = I_ERI_I2x3yz_Pz_Py_S+ABZ*I_ERI_H2x3y_Pz_Py_S;
  Double I_ERI_H2x2yz_D2z_Py_S = I_ERI_I2x2y2z_Pz_Py_S+ABZ*I_ERI_H2x2yz_Pz_Py_S;
  Double I_ERI_H2xy2z_D2z_Py_S = I_ERI_I2xy3z_Pz_Py_S+ABZ*I_ERI_H2xy2z_Pz_Py_S;
  Double I_ERI_H2x3z_D2z_Py_S = I_ERI_I2x4z_Pz_Py_S+ABZ*I_ERI_H2x3z_Pz_Py_S;
  Double I_ERI_Hx4y_D2z_Py_S = I_ERI_Ix4yz_Pz_Py_S+ABZ*I_ERI_Hx4y_Pz_Py_S;
  Double I_ERI_Hx3yz_D2z_Py_S = I_ERI_Ix3y2z_Pz_Py_S+ABZ*I_ERI_Hx3yz_Pz_Py_S;
  Double I_ERI_Hx2y2z_D2z_Py_S = I_ERI_Ix2y3z_Pz_Py_S+ABZ*I_ERI_Hx2y2z_Pz_Py_S;
  Double I_ERI_Hxy3z_D2z_Py_S = I_ERI_Ixy4z_Pz_Py_S+ABZ*I_ERI_Hxy3z_Pz_Py_S;
  Double I_ERI_Hx4z_D2z_Py_S = I_ERI_Ix5z_Pz_Py_S+ABZ*I_ERI_Hx4z_Pz_Py_S;
  Double I_ERI_H5y_D2z_Py_S = I_ERI_I5yz_Pz_Py_S+ABZ*I_ERI_H5y_Pz_Py_S;
  Double I_ERI_H4yz_D2z_Py_S = I_ERI_I4y2z_Pz_Py_S+ABZ*I_ERI_H4yz_Pz_Py_S;
  Double I_ERI_H3y2z_D2z_Py_S = I_ERI_I3y3z_Pz_Py_S+ABZ*I_ERI_H3y2z_Pz_Py_S;
  Double I_ERI_H2y3z_D2z_Py_S = I_ERI_I2y4z_Pz_Py_S+ABZ*I_ERI_H2y3z_Pz_Py_S;
  Double I_ERI_Hy4z_D2z_Py_S = I_ERI_Iy5z_Pz_Py_S+ABZ*I_ERI_Hy4z_Pz_Py_S;
  Double I_ERI_H5z_D2z_Py_S = I_ERI_I6z_Pz_Py_S+ABZ*I_ERI_H5z_Pz_Py_S;
  Double I_ERI_H5x_D2x_Pz_S = I_ERI_I6x_Px_Pz_S+ABX*I_ERI_H5x_Px_Pz_S;
  Double I_ERI_H4xy_D2x_Pz_S = I_ERI_I5xy_Px_Pz_S+ABX*I_ERI_H4xy_Px_Pz_S;
  Double I_ERI_H4xz_D2x_Pz_S = I_ERI_I5xz_Px_Pz_S+ABX*I_ERI_H4xz_Px_Pz_S;
  Double I_ERI_H3x2y_D2x_Pz_S = I_ERI_I4x2y_Px_Pz_S+ABX*I_ERI_H3x2y_Px_Pz_S;
  Double I_ERI_H3xyz_D2x_Pz_S = I_ERI_I4xyz_Px_Pz_S+ABX*I_ERI_H3xyz_Px_Pz_S;
  Double I_ERI_H3x2z_D2x_Pz_S = I_ERI_I4x2z_Px_Pz_S+ABX*I_ERI_H3x2z_Px_Pz_S;
  Double I_ERI_H2x3y_D2x_Pz_S = I_ERI_I3x3y_Px_Pz_S+ABX*I_ERI_H2x3y_Px_Pz_S;
  Double I_ERI_H2x2yz_D2x_Pz_S = I_ERI_I3x2yz_Px_Pz_S+ABX*I_ERI_H2x2yz_Px_Pz_S;
  Double I_ERI_H2xy2z_D2x_Pz_S = I_ERI_I3xy2z_Px_Pz_S+ABX*I_ERI_H2xy2z_Px_Pz_S;
  Double I_ERI_H2x3z_D2x_Pz_S = I_ERI_I3x3z_Px_Pz_S+ABX*I_ERI_H2x3z_Px_Pz_S;
  Double I_ERI_Hx4y_D2x_Pz_S = I_ERI_I2x4y_Px_Pz_S+ABX*I_ERI_Hx4y_Px_Pz_S;
  Double I_ERI_Hx3yz_D2x_Pz_S = I_ERI_I2x3yz_Px_Pz_S+ABX*I_ERI_Hx3yz_Px_Pz_S;
  Double I_ERI_Hx2y2z_D2x_Pz_S = I_ERI_I2x2y2z_Px_Pz_S+ABX*I_ERI_Hx2y2z_Px_Pz_S;
  Double I_ERI_Hxy3z_D2x_Pz_S = I_ERI_I2xy3z_Px_Pz_S+ABX*I_ERI_Hxy3z_Px_Pz_S;
  Double I_ERI_Hx4z_D2x_Pz_S = I_ERI_I2x4z_Px_Pz_S+ABX*I_ERI_Hx4z_Px_Pz_S;
  Double I_ERI_H5y_D2x_Pz_S = I_ERI_Ix5y_Px_Pz_S+ABX*I_ERI_H5y_Px_Pz_S;
  Double I_ERI_H4yz_D2x_Pz_S = I_ERI_Ix4yz_Px_Pz_S+ABX*I_ERI_H4yz_Px_Pz_S;
  Double I_ERI_H3y2z_D2x_Pz_S = I_ERI_Ix3y2z_Px_Pz_S+ABX*I_ERI_H3y2z_Px_Pz_S;
  Double I_ERI_H2y3z_D2x_Pz_S = I_ERI_Ix2y3z_Px_Pz_S+ABX*I_ERI_H2y3z_Px_Pz_S;
  Double I_ERI_Hy4z_D2x_Pz_S = I_ERI_Ixy4z_Px_Pz_S+ABX*I_ERI_Hy4z_Px_Pz_S;
  Double I_ERI_H5z_D2x_Pz_S = I_ERI_Ix5z_Px_Pz_S+ABX*I_ERI_H5z_Px_Pz_S;
  Double I_ERI_H5x_D2y_Pz_S = I_ERI_I5xy_Py_Pz_S+ABY*I_ERI_H5x_Py_Pz_S;
  Double I_ERI_H4xy_D2y_Pz_S = I_ERI_I4x2y_Py_Pz_S+ABY*I_ERI_H4xy_Py_Pz_S;
  Double I_ERI_H4xz_D2y_Pz_S = I_ERI_I4xyz_Py_Pz_S+ABY*I_ERI_H4xz_Py_Pz_S;
  Double I_ERI_H3x2y_D2y_Pz_S = I_ERI_I3x3y_Py_Pz_S+ABY*I_ERI_H3x2y_Py_Pz_S;
  Double I_ERI_H3xyz_D2y_Pz_S = I_ERI_I3x2yz_Py_Pz_S+ABY*I_ERI_H3xyz_Py_Pz_S;
  Double I_ERI_H3x2z_D2y_Pz_S = I_ERI_I3xy2z_Py_Pz_S+ABY*I_ERI_H3x2z_Py_Pz_S;
  Double I_ERI_H2x3y_D2y_Pz_S = I_ERI_I2x4y_Py_Pz_S+ABY*I_ERI_H2x3y_Py_Pz_S;
  Double I_ERI_H2x2yz_D2y_Pz_S = I_ERI_I2x3yz_Py_Pz_S+ABY*I_ERI_H2x2yz_Py_Pz_S;
  Double I_ERI_H2xy2z_D2y_Pz_S = I_ERI_I2x2y2z_Py_Pz_S+ABY*I_ERI_H2xy2z_Py_Pz_S;
  Double I_ERI_H2x3z_D2y_Pz_S = I_ERI_I2xy3z_Py_Pz_S+ABY*I_ERI_H2x3z_Py_Pz_S;
  Double I_ERI_Hx4y_D2y_Pz_S = I_ERI_Ix5y_Py_Pz_S+ABY*I_ERI_Hx4y_Py_Pz_S;
  Double I_ERI_Hx3yz_D2y_Pz_S = I_ERI_Ix4yz_Py_Pz_S+ABY*I_ERI_Hx3yz_Py_Pz_S;
  Double I_ERI_Hx2y2z_D2y_Pz_S = I_ERI_Ix3y2z_Py_Pz_S+ABY*I_ERI_Hx2y2z_Py_Pz_S;
  Double I_ERI_Hxy3z_D2y_Pz_S = I_ERI_Ix2y3z_Py_Pz_S+ABY*I_ERI_Hxy3z_Py_Pz_S;
  Double I_ERI_Hx4z_D2y_Pz_S = I_ERI_Ixy4z_Py_Pz_S+ABY*I_ERI_Hx4z_Py_Pz_S;
  Double I_ERI_H5y_D2y_Pz_S = I_ERI_I6y_Py_Pz_S+ABY*I_ERI_H5y_Py_Pz_S;
  Double I_ERI_H4yz_D2y_Pz_S = I_ERI_I5yz_Py_Pz_S+ABY*I_ERI_H4yz_Py_Pz_S;
  Double I_ERI_H3y2z_D2y_Pz_S = I_ERI_I4y2z_Py_Pz_S+ABY*I_ERI_H3y2z_Py_Pz_S;
  Double I_ERI_H2y3z_D2y_Pz_S = I_ERI_I3y3z_Py_Pz_S+ABY*I_ERI_H2y3z_Py_Pz_S;
  Double I_ERI_Hy4z_D2y_Pz_S = I_ERI_I2y4z_Py_Pz_S+ABY*I_ERI_Hy4z_Py_Pz_S;
  Double I_ERI_H5z_D2y_Pz_S = I_ERI_Iy5z_Py_Pz_S+ABY*I_ERI_H5z_Py_Pz_S;
  Double I_ERI_H5x_D2z_Pz_S = I_ERI_I5xz_Pz_Pz_S+ABZ*I_ERI_H5x_Pz_Pz_S;
  Double I_ERI_H4xy_D2z_Pz_S = I_ERI_I4xyz_Pz_Pz_S+ABZ*I_ERI_H4xy_Pz_Pz_S;
  Double I_ERI_H4xz_D2z_Pz_S = I_ERI_I4x2z_Pz_Pz_S+ABZ*I_ERI_H4xz_Pz_Pz_S;
  Double I_ERI_H3x2y_D2z_Pz_S = I_ERI_I3x2yz_Pz_Pz_S+ABZ*I_ERI_H3x2y_Pz_Pz_S;
  Double I_ERI_H3xyz_D2z_Pz_S = I_ERI_I3xy2z_Pz_Pz_S+ABZ*I_ERI_H3xyz_Pz_Pz_S;
  Double I_ERI_H3x2z_D2z_Pz_S = I_ERI_I3x3z_Pz_Pz_S+ABZ*I_ERI_H3x2z_Pz_Pz_S;
  Double I_ERI_H2x3y_D2z_Pz_S = I_ERI_I2x3yz_Pz_Pz_S+ABZ*I_ERI_H2x3y_Pz_Pz_S;
  Double I_ERI_H2x2yz_D2z_Pz_S = I_ERI_I2x2y2z_Pz_Pz_S+ABZ*I_ERI_H2x2yz_Pz_Pz_S;
  Double I_ERI_H2xy2z_D2z_Pz_S = I_ERI_I2xy3z_Pz_Pz_S+ABZ*I_ERI_H2xy2z_Pz_Pz_S;
  Double I_ERI_H2x3z_D2z_Pz_S = I_ERI_I2x4z_Pz_Pz_S+ABZ*I_ERI_H2x3z_Pz_Pz_S;
  Double I_ERI_Hx4y_D2z_Pz_S = I_ERI_Ix4yz_Pz_Pz_S+ABZ*I_ERI_Hx4y_Pz_Pz_S;
  Double I_ERI_Hx3yz_D2z_Pz_S = I_ERI_Ix3y2z_Pz_Pz_S+ABZ*I_ERI_Hx3yz_Pz_Pz_S;
  Double I_ERI_Hx2y2z_D2z_Pz_S = I_ERI_Ix2y3z_Pz_Pz_S+ABZ*I_ERI_Hx2y2z_Pz_Pz_S;
  Double I_ERI_Hxy3z_D2z_Pz_S = I_ERI_Ixy4z_Pz_Pz_S+ABZ*I_ERI_Hxy3z_Pz_Pz_S;
  Double I_ERI_Hx4z_D2z_Pz_S = I_ERI_Ix5z_Pz_Pz_S+ABZ*I_ERI_Hx4z_Pz_Pz_S;
  Double I_ERI_H5y_D2z_Pz_S = I_ERI_I5yz_Pz_Pz_S+ABZ*I_ERI_H5y_Pz_Pz_S;
  Double I_ERI_H4yz_D2z_Pz_S = I_ERI_I4y2z_Pz_Pz_S+ABZ*I_ERI_H4yz_Pz_Pz_S;
  Double I_ERI_H3y2z_D2z_Pz_S = I_ERI_I3y3z_Pz_Pz_S+ABZ*I_ERI_H3y2z_Pz_Pz_S;
  Double I_ERI_H2y3z_D2z_Pz_S = I_ERI_I2y4z_Pz_Pz_S+ABZ*I_ERI_H2y3z_Pz_Pz_S;
  Double I_ERI_Hy4z_D2z_Pz_S = I_ERI_Iy5z_Pz_Pz_S+ABZ*I_ERI_Hy4z_Pz_Pz_S;
  Double I_ERI_H5z_D2z_Pz_S = I_ERI_I6z_Pz_Pz_S+ABZ*I_ERI_H5z_Pz_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_G_F_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 90 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_D_P_S
   * RHS shell quartet name: SQ_ERI_G_D_P_S
   ************************************************************/
  Double I_ERI_G4x_F3x_Px_S = I_ERI_H5x_D2x_Px_S+ABX*I_ERI_G4x_D2x_Px_S;
  Double I_ERI_G3xy_F3x_Px_S = I_ERI_H4xy_D2x_Px_S+ABX*I_ERI_G3xy_D2x_Px_S;
  Double I_ERI_G3xz_F3x_Px_S = I_ERI_H4xz_D2x_Px_S+ABX*I_ERI_G3xz_D2x_Px_S;
  Double I_ERI_G2x2y_F3x_Px_S = I_ERI_H3x2y_D2x_Px_S+ABX*I_ERI_G2x2y_D2x_Px_S;
  Double I_ERI_G2xyz_F3x_Px_S = I_ERI_H3xyz_D2x_Px_S+ABX*I_ERI_G2xyz_D2x_Px_S;
  Double I_ERI_G2x2z_F3x_Px_S = I_ERI_H3x2z_D2x_Px_S+ABX*I_ERI_G2x2z_D2x_Px_S;
  Double I_ERI_Gx3y_F3x_Px_S = I_ERI_H2x3y_D2x_Px_S+ABX*I_ERI_Gx3y_D2x_Px_S;
  Double I_ERI_Gx2yz_F3x_Px_S = I_ERI_H2x2yz_D2x_Px_S+ABX*I_ERI_Gx2yz_D2x_Px_S;
  Double I_ERI_Gxy2z_F3x_Px_S = I_ERI_H2xy2z_D2x_Px_S+ABX*I_ERI_Gxy2z_D2x_Px_S;
  Double I_ERI_Gx3z_F3x_Px_S = I_ERI_H2x3z_D2x_Px_S+ABX*I_ERI_Gx3z_D2x_Px_S;
  Double I_ERI_G4y_F3x_Px_S = I_ERI_Hx4y_D2x_Px_S+ABX*I_ERI_G4y_D2x_Px_S;
  Double I_ERI_G3yz_F3x_Px_S = I_ERI_Hx3yz_D2x_Px_S+ABX*I_ERI_G3yz_D2x_Px_S;
  Double I_ERI_G2y2z_F3x_Px_S = I_ERI_Hx2y2z_D2x_Px_S+ABX*I_ERI_G2y2z_D2x_Px_S;
  Double I_ERI_Gy3z_F3x_Px_S = I_ERI_Hxy3z_D2x_Px_S+ABX*I_ERI_Gy3z_D2x_Px_S;
  Double I_ERI_G4z_F3x_Px_S = I_ERI_Hx4z_D2x_Px_S+ABX*I_ERI_G4z_D2x_Px_S;
  Double I_ERI_G4x_F2xy_Px_S = I_ERI_H4xy_D2x_Px_S+ABY*I_ERI_G4x_D2x_Px_S;
  Double I_ERI_G3xy_F2xy_Px_S = I_ERI_H3x2y_D2x_Px_S+ABY*I_ERI_G3xy_D2x_Px_S;
  Double I_ERI_G3xz_F2xy_Px_S = I_ERI_H3xyz_D2x_Px_S+ABY*I_ERI_G3xz_D2x_Px_S;
  Double I_ERI_G2x2y_F2xy_Px_S = I_ERI_H2x3y_D2x_Px_S+ABY*I_ERI_G2x2y_D2x_Px_S;
  Double I_ERI_G2xyz_F2xy_Px_S = I_ERI_H2x2yz_D2x_Px_S+ABY*I_ERI_G2xyz_D2x_Px_S;
  Double I_ERI_G2x2z_F2xy_Px_S = I_ERI_H2xy2z_D2x_Px_S+ABY*I_ERI_G2x2z_D2x_Px_S;
  Double I_ERI_Gx3y_F2xy_Px_S = I_ERI_Hx4y_D2x_Px_S+ABY*I_ERI_Gx3y_D2x_Px_S;
  Double I_ERI_Gx2yz_F2xy_Px_S = I_ERI_Hx3yz_D2x_Px_S+ABY*I_ERI_Gx2yz_D2x_Px_S;
  Double I_ERI_Gxy2z_F2xy_Px_S = I_ERI_Hx2y2z_D2x_Px_S+ABY*I_ERI_Gxy2z_D2x_Px_S;
  Double I_ERI_Gx3z_F2xy_Px_S = I_ERI_Hxy3z_D2x_Px_S+ABY*I_ERI_Gx3z_D2x_Px_S;
  Double I_ERI_G4y_F2xy_Px_S = I_ERI_H5y_D2x_Px_S+ABY*I_ERI_G4y_D2x_Px_S;
  Double I_ERI_G3yz_F2xy_Px_S = I_ERI_H4yz_D2x_Px_S+ABY*I_ERI_G3yz_D2x_Px_S;
  Double I_ERI_G2y2z_F2xy_Px_S = I_ERI_H3y2z_D2x_Px_S+ABY*I_ERI_G2y2z_D2x_Px_S;
  Double I_ERI_Gy3z_F2xy_Px_S = I_ERI_H2y3z_D2x_Px_S+ABY*I_ERI_Gy3z_D2x_Px_S;
  Double I_ERI_G4z_F2xy_Px_S = I_ERI_Hy4z_D2x_Px_S+ABY*I_ERI_G4z_D2x_Px_S;
  Double I_ERI_G4x_F2xz_Px_S = I_ERI_H4xz_D2x_Px_S+ABZ*I_ERI_G4x_D2x_Px_S;
  Double I_ERI_G3xy_F2xz_Px_S = I_ERI_H3xyz_D2x_Px_S+ABZ*I_ERI_G3xy_D2x_Px_S;
  Double I_ERI_G3xz_F2xz_Px_S = I_ERI_H3x2z_D2x_Px_S+ABZ*I_ERI_G3xz_D2x_Px_S;
  Double I_ERI_G2x2y_F2xz_Px_S = I_ERI_H2x2yz_D2x_Px_S+ABZ*I_ERI_G2x2y_D2x_Px_S;
  Double I_ERI_G2xyz_F2xz_Px_S = I_ERI_H2xy2z_D2x_Px_S+ABZ*I_ERI_G2xyz_D2x_Px_S;
  Double I_ERI_G2x2z_F2xz_Px_S = I_ERI_H2x3z_D2x_Px_S+ABZ*I_ERI_G2x2z_D2x_Px_S;
  Double I_ERI_Gx3y_F2xz_Px_S = I_ERI_Hx3yz_D2x_Px_S+ABZ*I_ERI_Gx3y_D2x_Px_S;
  Double I_ERI_Gx2yz_F2xz_Px_S = I_ERI_Hx2y2z_D2x_Px_S+ABZ*I_ERI_Gx2yz_D2x_Px_S;
  Double I_ERI_Gxy2z_F2xz_Px_S = I_ERI_Hxy3z_D2x_Px_S+ABZ*I_ERI_Gxy2z_D2x_Px_S;
  Double I_ERI_Gx3z_F2xz_Px_S = I_ERI_Hx4z_D2x_Px_S+ABZ*I_ERI_Gx3z_D2x_Px_S;
  Double I_ERI_G4y_F2xz_Px_S = I_ERI_H4yz_D2x_Px_S+ABZ*I_ERI_G4y_D2x_Px_S;
  Double I_ERI_G3yz_F2xz_Px_S = I_ERI_H3y2z_D2x_Px_S+ABZ*I_ERI_G3yz_D2x_Px_S;
  Double I_ERI_G2y2z_F2xz_Px_S = I_ERI_H2y3z_D2x_Px_S+ABZ*I_ERI_G2y2z_D2x_Px_S;
  Double I_ERI_Gy3z_F2xz_Px_S = I_ERI_Hy4z_D2x_Px_S+ABZ*I_ERI_Gy3z_D2x_Px_S;
  Double I_ERI_G4z_F2xz_Px_S = I_ERI_H5z_D2x_Px_S+ABZ*I_ERI_G4z_D2x_Px_S;
  Double I_ERI_G4x_Fx2y_Px_S = I_ERI_H5x_D2y_Px_S+ABX*I_ERI_G4x_D2y_Px_S;
  Double I_ERI_G3xy_Fx2y_Px_S = I_ERI_H4xy_D2y_Px_S+ABX*I_ERI_G3xy_D2y_Px_S;
  Double I_ERI_G3xz_Fx2y_Px_S = I_ERI_H4xz_D2y_Px_S+ABX*I_ERI_G3xz_D2y_Px_S;
  Double I_ERI_G2x2y_Fx2y_Px_S = I_ERI_H3x2y_D2y_Px_S+ABX*I_ERI_G2x2y_D2y_Px_S;
  Double I_ERI_G2xyz_Fx2y_Px_S = I_ERI_H3xyz_D2y_Px_S+ABX*I_ERI_G2xyz_D2y_Px_S;
  Double I_ERI_G2x2z_Fx2y_Px_S = I_ERI_H3x2z_D2y_Px_S+ABX*I_ERI_G2x2z_D2y_Px_S;
  Double I_ERI_Gx3y_Fx2y_Px_S = I_ERI_H2x3y_D2y_Px_S+ABX*I_ERI_Gx3y_D2y_Px_S;
  Double I_ERI_Gx2yz_Fx2y_Px_S = I_ERI_H2x2yz_D2y_Px_S+ABX*I_ERI_Gx2yz_D2y_Px_S;
  Double I_ERI_Gxy2z_Fx2y_Px_S = I_ERI_H2xy2z_D2y_Px_S+ABX*I_ERI_Gxy2z_D2y_Px_S;
  Double I_ERI_Gx3z_Fx2y_Px_S = I_ERI_H2x3z_D2y_Px_S+ABX*I_ERI_Gx3z_D2y_Px_S;
  Double I_ERI_G4y_Fx2y_Px_S = I_ERI_Hx4y_D2y_Px_S+ABX*I_ERI_G4y_D2y_Px_S;
  Double I_ERI_G3yz_Fx2y_Px_S = I_ERI_Hx3yz_D2y_Px_S+ABX*I_ERI_G3yz_D2y_Px_S;
  Double I_ERI_G2y2z_Fx2y_Px_S = I_ERI_Hx2y2z_D2y_Px_S+ABX*I_ERI_G2y2z_D2y_Px_S;
  Double I_ERI_Gy3z_Fx2y_Px_S = I_ERI_Hxy3z_D2y_Px_S+ABX*I_ERI_Gy3z_D2y_Px_S;
  Double I_ERI_G4z_Fx2y_Px_S = I_ERI_Hx4z_D2y_Px_S+ABX*I_ERI_G4z_D2y_Px_S;
  Double I_ERI_G4x_Fx2z_Px_S = I_ERI_H5x_D2z_Px_S+ABX*I_ERI_G4x_D2z_Px_S;
  Double I_ERI_G3xy_Fx2z_Px_S = I_ERI_H4xy_D2z_Px_S+ABX*I_ERI_G3xy_D2z_Px_S;
  Double I_ERI_G3xz_Fx2z_Px_S = I_ERI_H4xz_D2z_Px_S+ABX*I_ERI_G3xz_D2z_Px_S;
  Double I_ERI_G2x2y_Fx2z_Px_S = I_ERI_H3x2y_D2z_Px_S+ABX*I_ERI_G2x2y_D2z_Px_S;
  Double I_ERI_G2xyz_Fx2z_Px_S = I_ERI_H3xyz_D2z_Px_S+ABX*I_ERI_G2xyz_D2z_Px_S;
  Double I_ERI_G2x2z_Fx2z_Px_S = I_ERI_H3x2z_D2z_Px_S+ABX*I_ERI_G2x2z_D2z_Px_S;
  Double I_ERI_Gx3y_Fx2z_Px_S = I_ERI_H2x3y_D2z_Px_S+ABX*I_ERI_Gx3y_D2z_Px_S;
  Double I_ERI_Gx2yz_Fx2z_Px_S = I_ERI_H2x2yz_D2z_Px_S+ABX*I_ERI_Gx2yz_D2z_Px_S;
  Double I_ERI_Gxy2z_Fx2z_Px_S = I_ERI_H2xy2z_D2z_Px_S+ABX*I_ERI_Gxy2z_D2z_Px_S;
  Double I_ERI_Gx3z_Fx2z_Px_S = I_ERI_H2x3z_D2z_Px_S+ABX*I_ERI_Gx3z_D2z_Px_S;
  Double I_ERI_G4y_Fx2z_Px_S = I_ERI_Hx4y_D2z_Px_S+ABX*I_ERI_G4y_D2z_Px_S;
  Double I_ERI_G3yz_Fx2z_Px_S = I_ERI_Hx3yz_D2z_Px_S+ABX*I_ERI_G3yz_D2z_Px_S;
  Double I_ERI_G2y2z_Fx2z_Px_S = I_ERI_Hx2y2z_D2z_Px_S+ABX*I_ERI_G2y2z_D2z_Px_S;
  Double I_ERI_Gy3z_Fx2z_Px_S = I_ERI_Hxy3z_D2z_Px_S+ABX*I_ERI_Gy3z_D2z_Px_S;
  Double I_ERI_G4z_Fx2z_Px_S = I_ERI_Hx4z_D2z_Px_S+ABX*I_ERI_G4z_D2z_Px_S;
  Double I_ERI_G4x_F3y_Px_S = I_ERI_H4xy_D2y_Px_S+ABY*I_ERI_G4x_D2y_Px_S;
  Double I_ERI_G3xy_F3y_Px_S = I_ERI_H3x2y_D2y_Px_S+ABY*I_ERI_G3xy_D2y_Px_S;
  Double I_ERI_G3xz_F3y_Px_S = I_ERI_H3xyz_D2y_Px_S+ABY*I_ERI_G3xz_D2y_Px_S;
  Double I_ERI_G2x2y_F3y_Px_S = I_ERI_H2x3y_D2y_Px_S+ABY*I_ERI_G2x2y_D2y_Px_S;
  Double I_ERI_G2xyz_F3y_Px_S = I_ERI_H2x2yz_D2y_Px_S+ABY*I_ERI_G2xyz_D2y_Px_S;
  Double I_ERI_G2x2z_F3y_Px_S = I_ERI_H2xy2z_D2y_Px_S+ABY*I_ERI_G2x2z_D2y_Px_S;
  Double I_ERI_Gx3y_F3y_Px_S = I_ERI_Hx4y_D2y_Px_S+ABY*I_ERI_Gx3y_D2y_Px_S;
  Double I_ERI_Gx2yz_F3y_Px_S = I_ERI_Hx3yz_D2y_Px_S+ABY*I_ERI_Gx2yz_D2y_Px_S;
  Double I_ERI_Gxy2z_F3y_Px_S = I_ERI_Hx2y2z_D2y_Px_S+ABY*I_ERI_Gxy2z_D2y_Px_S;
  Double I_ERI_Gx3z_F3y_Px_S = I_ERI_Hxy3z_D2y_Px_S+ABY*I_ERI_Gx3z_D2y_Px_S;
  Double I_ERI_G4y_F3y_Px_S = I_ERI_H5y_D2y_Px_S+ABY*I_ERI_G4y_D2y_Px_S;
  Double I_ERI_G3yz_F3y_Px_S = I_ERI_H4yz_D2y_Px_S+ABY*I_ERI_G3yz_D2y_Px_S;
  Double I_ERI_G2y2z_F3y_Px_S = I_ERI_H3y2z_D2y_Px_S+ABY*I_ERI_G2y2z_D2y_Px_S;
  Double I_ERI_Gy3z_F3y_Px_S = I_ERI_H2y3z_D2y_Px_S+ABY*I_ERI_Gy3z_D2y_Px_S;
  Double I_ERI_G4z_F3y_Px_S = I_ERI_Hy4z_D2y_Px_S+ABY*I_ERI_G4z_D2y_Px_S;
  Double I_ERI_G4x_F2yz_Px_S = I_ERI_H4xz_D2y_Px_S+ABZ*I_ERI_G4x_D2y_Px_S;
  Double I_ERI_G3xy_F2yz_Px_S = I_ERI_H3xyz_D2y_Px_S+ABZ*I_ERI_G3xy_D2y_Px_S;
  Double I_ERI_G3xz_F2yz_Px_S = I_ERI_H3x2z_D2y_Px_S+ABZ*I_ERI_G3xz_D2y_Px_S;
  Double I_ERI_G2x2y_F2yz_Px_S = I_ERI_H2x2yz_D2y_Px_S+ABZ*I_ERI_G2x2y_D2y_Px_S;
  Double I_ERI_G2xyz_F2yz_Px_S = I_ERI_H2xy2z_D2y_Px_S+ABZ*I_ERI_G2xyz_D2y_Px_S;
  Double I_ERI_G2x2z_F2yz_Px_S = I_ERI_H2x3z_D2y_Px_S+ABZ*I_ERI_G2x2z_D2y_Px_S;
  Double I_ERI_Gx3y_F2yz_Px_S = I_ERI_Hx3yz_D2y_Px_S+ABZ*I_ERI_Gx3y_D2y_Px_S;
  Double I_ERI_Gx2yz_F2yz_Px_S = I_ERI_Hx2y2z_D2y_Px_S+ABZ*I_ERI_Gx2yz_D2y_Px_S;
  Double I_ERI_Gxy2z_F2yz_Px_S = I_ERI_Hxy3z_D2y_Px_S+ABZ*I_ERI_Gxy2z_D2y_Px_S;
  Double I_ERI_Gx3z_F2yz_Px_S = I_ERI_Hx4z_D2y_Px_S+ABZ*I_ERI_Gx3z_D2y_Px_S;
  Double I_ERI_G4y_F2yz_Px_S = I_ERI_H4yz_D2y_Px_S+ABZ*I_ERI_G4y_D2y_Px_S;
  Double I_ERI_G3yz_F2yz_Px_S = I_ERI_H3y2z_D2y_Px_S+ABZ*I_ERI_G3yz_D2y_Px_S;
  Double I_ERI_G2y2z_F2yz_Px_S = I_ERI_H2y3z_D2y_Px_S+ABZ*I_ERI_G2y2z_D2y_Px_S;
  Double I_ERI_Gy3z_F2yz_Px_S = I_ERI_Hy4z_D2y_Px_S+ABZ*I_ERI_Gy3z_D2y_Px_S;
  Double I_ERI_G4z_F2yz_Px_S = I_ERI_H5z_D2y_Px_S+ABZ*I_ERI_G4z_D2y_Px_S;
  Double I_ERI_G4x_F3z_Px_S = I_ERI_H4xz_D2z_Px_S+ABZ*I_ERI_G4x_D2z_Px_S;
  Double I_ERI_G3xy_F3z_Px_S = I_ERI_H3xyz_D2z_Px_S+ABZ*I_ERI_G3xy_D2z_Px_S;
  Double I_ERI_G3xz_F3z_Px_S = I_ERI_H3x2z_D2z_Px_S+ABZ*I_ERI_G3xz_D2z_Px_S;
  Double I_ERI_G2x2y_F3z_Px_S = I_ERI_H2x2yz_D2z_Px_S+ABZ*I_ERI_G2x2y_D2z_Px_S;
  Double I_ERI_G2xyz_F3z_Px_S = I_ERI_H2xy2z_D2z_Px_S+ABZ*I_ERI_G2xyz_D2z_Px_S;
  Double I_ERI_G2x2z_F3z_Px_S = I_ERI_H2x3z_D2z_Px_S+ABZ*I_ERI_G2x2z_D2z_Px_S;
  Double I_ERI_Gx3y_F3z_Px_S = I_ERI_Hx3yz_D2z_Px_S+ABZ*I_ERI_Gx3y_D2z_Px_S;
  Double I_ERI_Gx2yz_F3z_Px_S = I_ERI_Hx2y2z_D2z_Px_S+ABZ*I_ERI_Gx2yz_D2z_Px_S;
  Double I_ERI_Gxy2z_F3z_Px_S = I_ERI_Hxy3z_D2z_Px_S+ABZ*I_ERI_Gxy2z_D2z_Px_S;
  Double I_ERI_Gx3z_F3z_Px_S = I_ERI_Hx4z_D2z_Px_S+ABZ*I_ERI_Gx3z_D2z_Px_S;
  Double I_ERI_G4y_F3z_Px_S = I_ERI_H4yz_D2z_Px_S+ABZ*I_ERI_G4y_D2z_Px_S;
  Double I_ERI_G3yz_F3z_Px_S = I_ERI_H3y2z_D2z_Px_S+ABZ*I_ERI_G3yz_D2z_Px_S;
  Double I_ERI_G2y2z_F3z_Px_S = I_ERI_H2y3z_D2z_Px_S+ABZ*I_ERI_G2y2z_D2z_Px_S;
  Double I_ERI_Gy3z_F3z_Px_S = I_ERI_Hy4z_D2z_Px_S+ABZ*I_ERI_Gy3z_D2z_Px_S;
  Double I_ERI_G4z_F3z_Px_S = I_ERI_H5z_D2z_Px_S+ABZ*I_ERI_G4z_D2z_Px_S;
  Double I_ERI_G4x_F3x_Py_S = I_ERI_H5x_D2x_Py_S+ABX*I_ERI_G4x_D2x_Py_S;
  Double I_ERI_G3xy_F3x_Py_S = I_ERI_H4xy_D2x_Py_S+ABX*I_ERI_G3xy_D2x_Py_S;
  Double I_ERI_G3xz_F3x_Py_S = I_ERI_H4xz_D2x_Py_S+ABX*I_ERI_G3xz_D2x_Py_S;
  Double I_ERI_G2x2y_F3x_Py_S = I_ERI_H3x2y_D2x_Py_S+ABX*I_ERI_G2x2y_D2x_Py_S;
  Double I_ERI_G2xyz_F3x_Py_S = I_ERI_H3xyz_D2x_Py_S+ABX*I_ERI_G2xyz_D2x_Py_S;
  Double I_ERI_G2x2z_F3x_Py_S = I_ERI_H3x2z_D2x_Py_S+ABX*I_ERI_G2x2z_D2x_Py_S;
  Double I_ERI_Gx3y_F3x_Py_S = I_ERI_H2x3y_D2x_Py_S+ABX*I_ERI_Gx3y_D2x_Py_S;
  Double I_ERI_Gx2yz_F3x_Py_S = I_ERI_H2x2yz_D2x_Py_S+ABX*I_ERI_Gx2yz_D2x_Py_S;
  Double I_ERI_Gxy2z_F3x_Py_S = I_ERI_H2xy2z_D2x_Py_S+ABX*I_ERI_Gxy2z_D2x_Py_S;
  Double I_ERI_Gx3z_F3x_Py_S = I_ERI_H2x3z_D2x_Py_S+ABX*I_ERI_Gx3z_D2x_Py_S;
  Double I_ERI_G4y_F3x_Py_S = I_ERI_Hx4y_D2x_Py_S+ABX*I_ERI_G4y_D2x_Py_S;
  Double I_ERI_G3yz_F3x_Py_S = I_ERI_Hx3yz_D2x_Py_S+ABX*I_ERI_G3yz_D2x_Py_S;
  Double I_ERI_G2y2z_F3x_Py_S = I_ERI_Hx2y2z_D2x_Py_S+ABX*I_ERI_G2y2z_D2x_Py_S;
  Double I_ERI_Gy3z_F3x_Py_S = I_ERI_Hxy3z_D2x_Py_S+ABX*I_ERI_Gy3z_D2x_Py_S;
  Double I_ERI_G4z_F3x_Py_S = I_ERI_Hx4z_D2x_Py_S+ABX*I_ERI_G4z_D2x_Py_S;
  Double I_ERI_G4x_F2xy_Py_S = I_ERI_H4xy_D2x_Py_S+ABY*I_ERI_G4x_D2x_Py_S;
  Double I_ERI_G3xy_F2xy_Py_S = I_ERI_H3x2y_D2x_Py_S+ABY*I_ERI_G3xy_D2x_Py_S;
  Double I_ERI_G3xz_F2xy_Py_S = I_ERI_H3xyz_D2x_Py_S+ABY*I_ERI_G3xz_D2x_Py_S;
  Double I_ERI_G2x2y_F2xy_Py_S = I_ERI_H2x3y_D2x_Py_S+ABY*I_ERI_G2x2y_D2x_Py_S;
  Double I_ERI_G2xyz_F2xy_Py_S = I_ERI_H2x2yz_D2x_Py_S+ABY*I_ERI_G2xyz_D2x_Py_S;
  Double I_ERI_G2x2z_F2xy_Py_S = I_ERI_H2xy2z_D2x_Py_S+ABY*I_ERI_G2x2z_D2x_Py_S;
  Double I_ERI_Gx3y_F2xy_Py_S = I_ERI_Hx4y_D2x_Py_S+ABY*I_ERI_Gx3y_D2x_Py_S;
  Double I_ERI_Gx2yz_F2xy_Py_S = I_ERI_Hx3yz_D2x_Py_S+ABY*I_ERI_Gx2yz_D2x_Py_S;
  Double I_ERI_Gxy2z_F2xy_Py_S = I_ERI_Hx2y2z_D2x_Py_S+ABY*I_ERI_Gxy2z_D2x_Py_S;
  Double I_ERI_Gx3z_F2xy_Py_S = I_ERI_Hxy3z_D2x_Py_S+ABY*I_ERI_Gx3z_D2x_Py_S;
  Double I_ERI_G4y_F2xy_Py_S = I_ERI_H5y_D2x_Py_S+ABY*I_ERI_G4y_D2x_Py_S;
  Double I_ERI_G3yz_F2xy_Py_S = I_ERI_H4yz_D2x_Py_S+ABY*I_ERI_G3yz_D2x_Py_S;
  Double I_ERI_G2y2z_F2xy_Py_S = I_ERI_H3y2z_D2x_Py_S+ABY*I_ERI_G2y2z_D2x_Py_S;
  Double I_ERI_Gy3z_F2xy_Py_S = I_ERI_H2y3z_D2x_Py_S+ABY*I_ERI_Gy3z_D2x_Py_S;
  Double I_ERI_G4z_F2xy_Py_S = I_ERI_Hy4z_D2x_Py_S+ABY*I_ERI_G4z_D2x_Py_S;
  Double I_ERI_G4x_F2xz_Py_S = I_ERI_H4xz_D2x_Py_S+ABZ*I_ERI_G4x_D2x_Py_S;
  Double I_ERI_G3xy_F2xz_Py_S = I_ERI_H3xyz_D2x_Py_S+ABZ*I_ERI_G3xy_D2x_Py_S;
  Double I_ERI_G3xz_F2xz_Py_S = I_ERI_H3x2z_D2x_Py_S+ABZ*I_ERI_G3xz_D2x_Py_S;
  Double I_ERI_G2x2y_F2xz_Py_S = I_ERI_H2x2yz_D2x_Py_S+ABZ*I_ERI_G2x2y_D2x_Py_S;
  Double I_ERI_G2xyz_F2xz_Py_S = I_ERI_H2xy2z_D2x_Py_S+ABZ*I_ERI_G2xyz_D2x_Py_S;
  Double I_ERI_G2x2z_F2xz_Py_S = I_ERI_H2x3z_D2x_Py_S+ABZ*I_ERI_G2x2z_D2x_Py_S;
  Double I_ERI_Gx3y_F2xz_Py_S = I_ERI_Hx3yz_D2x_Py_S+ABZ*I_ERI_Gx3y_D2x_Py_S;
  Double I_ERI_Gx2yz_F2xz_Py_S = I_ERI_Hx2y2z_D2x_Py_S+ABZ*I_ERI_Gx2yz_D2x_Py_S;
  Double I_ERI_Gxy2z_F2xz_Py_S = I_ERI_Hxy3z_D2x_Py_S+ABZ*I_ERI_Gxy2z_D2x_Py_S;
  Double I_ERI_Gx3z_F2xz_Py_S = I_ERI_Hx4z_D2x_Py_S+ABZ*I_ERI_Gx3z_D2x_Py_S;
  Double I_ERI_G4y_F2xz_Py_S = I_ERI_H4yz_D2x_Py_S+ABZ*I_ERI_G4y_D2x_Py_S;
  Double I_ERI_G3yz_F2xz_Py_S = I_ERI_H3y2z_D2x_Py_S+ABZ*I_ERI_G3yz_D2x_Py_S;
  Double I_ERI_G2y2z_F2xz_Py_S = I_ERI_H2y3z_D2x_Py_S+ABZ*I_ERI_G2y2z_D2x_Py_S;
  Double I_ERI_Gy3z_F2xz_Py_S = I_ERI_Hy4z_D2x_Py_S+ABZ*I_ERI_Gy3z_D2x_Py_S;
  Double I_ERI_G4z_F2xz_Py_S = I_ERI_H5z_D2x_Py_S+ABZ*I_ERI_G4z_D2x_Py_S;
  Double I_ERI_G4x_Fx2y_Py_S = I_ERI_H5x_D2y_Py_S+ABX*I_ERI_G4x_D2y_Py_S;
  Double I_ERI_G3xy_Fx2y_Py_S = I_ERI_H4xy_D2y_Py_S+ABX*I_ERI_G3xy_D2y_Py_S;
  Double I_ERI_G3xz_Fx2y_Py_S = I_ERI_H4xz_D2y_Py_S+ABX*I_ERI_G3xz_D2y_Py_S;
  Double I_ERI_G2x2y_Fx2y_Py_S = I_ERI_H3x2y_D2y_Py_S+ABX*I_ERI_G2x2y_D2y_Py_S;
  Double I_ERI_G2xyz_Fx2y_Py_S = I_ERI_H3xyz_D2y_Py_S+ABX*I_ERI_G2xyz_D2y_Py_S;
  Double I_ERI_G2x2z_Fx2y_Py_S = I_ERI_H3x2z_D2y_Py_S+ABX*I_ERI_G2x2z_D2y_Py_S;
  Double I_ERI_Gx3y_Fx2y_Py_S = I_ERI_H2x3y_D2y_Py_S+ABX*I_ERI_Gx3y_D2y_Py_S;
  Double I_ERI_Gx2yz_Fx2y_Py_S = I_ERI_H2x2yz_D2y_Py_S+ABX*I_ERI_Gx2yz_D2y_Py_S;
  Double I_ERI_Gxy2z_Fx2y_Py_S = I_ERI_H2xy2z_D2y_Py_S+ABX*I_ERI_Gxy2z_D2y_Py_S;
  Double I_ERI_Gx3z_Fx2y_Py_S = I_ERI_H2x3z_D2y_Py_S+ABX*I_ERI_Gx3z_D2y_Py_S;
  Double I_ERI_G4y_Fx2y_Py_S = I_ERI_Hx4y_D2y_Py_S+ABX*I_ERI_G4y_D2y_Py_S;
  Double I_ERI_G3yz_Fx2y_Py_S = I_ERI_Hx3yz_D2y_Py_S+ABX*I_ERI_G3yz_D2y_Py_S;
  Double I_ERI_G2y2z_Fx2y_Py_S = I_ERI_Hx2y2z_D2y_Py_S+ABX*I_ERI_G2y2z_D2y_Py_S;
  Double I_ERI_Gy3z_Fx2y_Py_S = I_ERI_Hxy3z_D2y_Py_S+ABX*I_ERI_Gy3z_D2y_Py_S;
  Double I_ERI_G4z_Fx2y_Py_S = I_ERI_Hx4z_D2y_Py_S+ABX*I_ERI_G4z_D2y_Py_S;
  Double I_ERI_G4x_Fx2z_Py_S = I_ERI_H5x_D2z_Py_S+ABX*I_ERI_G4x_D2z_Py_S;
  Double I_ERI_G3xy_Fx2z_Py_S = I_ERI_H4xy_D2z_Py_S+ABX*I_ERI_G3xy_D2z_Py_S;
  Double I_ERI_G3xz_Fx2z_Py_S = I_ERI_H4xz_D2z_Py_S+ABX*I_ERI_G3xz_D2z_Py_S;
  Double I_ERI_G2x2y_Fx2z_Py_S = I_ERI_H3x2y_D2z_Py_S+ABX*I_ERI_G2x2y_D2z_Py_S;
  Double I_ERI_G2xyz_Fx2z_Py_S = I_ERI_H3xyz_D2z_Py_S+ABX*I_ERI_G2xyz_D2z_Py_S;
  Double I_ERI_G2x2z_Fx2z_Py_S = I_ERI_H3x2z_D2z_Py_S+ABX*I_ERI_G2x2z_D2z_Py_S;
  Double I_ERI_Gx3y_Fx2z_Py_S = I_ERI_H2x3y_D2z_Py_S+ABX*I_ERI_Gx3y_D2z_Py_S;
  Double I_ERI_Gx2yz_Fx2z_Py_S = I_ERI_H2x2yz_D2z_Py_S+ABX*I_ERI_Gx2yz_D2z_Py_S;
  Double I_ERI_Gxy2z_Fx2z_Py_S = I_ERI_H2xy2z_D2z_Py_S+ABX*I_ERI_Gxy2z_D2z_Py_S;
  Double I_ERI_Gx3z_Fx2z_Py_S = I_ERI_H2x3z_D2z_Py_S+ABX*I_ERI_Gx3z_D2z_Py_S;
  Double I_ERI_G4y_Fx2z_Py_S = I_ERI_Hx4y_D2z_Py_S+ABX*I_ERI_G4y_D2z_Py_S;
  Double I_ERI_G3yz_Fx2z_Py_S = I_ERI_Hx3yz_D2z_Py_S+ABX*I_ERI_G3yz_D2z_Py_S;
  Double I_ERI_G2y2z_Fx2z_Py_S = I_ERI_Hx2y2z_D2z_Py_S+ABX*I_ERI_G2y2z_D2z_Py_S;
  Double I_ERI_Gy3z_Fx2z_Py_S = I_ERI_Hxy3z_D2z_Py_S+ABX*I_ERI_Gy3z_D2z_Py_S;
  Double I_ERI_G4z_Fx2z_Py_S = I_ERI_Hx4z_D2z_Py_S+ABX*I_ERI_G4z_D2z_Py_S;
  Double I_ERI_G4x_F3y_Py_S = I_ERI_H4xy_D2y_Py_S+ABY*I_ERI_G4x_D2y_Py_S;
  Double I_ERI_G3xy_F3y_Py_S = I_ERI_H3x2y_D2y_Py_S+ABY*I_ERI_G3xy_D2y_Py_S;
  Double I_ERI_G3xz_F3y_Py_S = I_ERI_H3xyz_D2y_Py_S+ABY*I_ERI_G3xz_D2y_Py_S;
  Double I_ERI_G2x2y_F3y_Py_S = I_ERI_H2x3y_D2y_Py_S+ABY*I_ERI_G2x2y_D2y_Py_S;
  Double I_ERI_G2xyz_F3y_Py_S = I_ERI_H2x2yz_D2y_Py_S+ABY*I_ERI_G2xyz_D2y_Py_S;
  Double I_ERI_G2x2z_F3y_Py_S = I_ERI_H2xy2z_D2y_Py_S+ABY*I_ERI_G2x2z_D2y_Py_S;
  Double I_ERI_Gx3y_F3y_Py_S = I_ERI_Hx4y_D2y_Py_S+ABY*I_ERI_Gx3y_D2y_Py_S;
  Double I_ERI_Gx2yz_F3y_Py_S = I_ERI_Hx3yz_D2y_Py_S+ABY*I_ERI_Gx2yz_D2y_Py_S;
  Double I_ERI_Gxy2z_F3y_Py_S = I_ERI_Hx2y2z_D2y_Py_S+ABY*I_ERI_Gxy2z_D2y_Py_S;
  Double I_ERI_Gx3z_F3y_Py_S = I_ERI_Hxy3z_D2y_Py_S+ABY*I_ERI_Gx3z_D2y_Py_S;
  Double I_ERI_G4y_F3y_Py_S = I_ERI_H5y_D2y_Py_S+ABY*I_ERI_G4y_D2y_Py_S;
  Double I_ERI_G3yz_F3y_Py_S = I_ERI_H4yz_D2y_Py_S+ABY*I_ERI_G3yz_D2y_Py_S;
  Double I_ERI_G2y2z_F3y_Py_S = I_ERI_H3y2z_D2y_Py_S+ABY*I_ERI_G2y2z_D2y_Py_S;
  Double I_ERI_Gy3z_F3y_Py_S = I_ERI_H2y3z_D2y_Py_S+ABY*I_ERI_Gy3z_D2y_Py_S;
  Double I_ERI_G4z_F3y_Py_S = I_ERI_Hy4z_D2y_Py_S+ABY*I_ERI_G4z_D2y_Py_S;
  Double I_ERI_G4x_F2yz_Py_S = I_ERI_H4xz_D2y_Py_S+ABZ*I_ERI_G4x_D2y_Py_S;
  Double I_ERI_G3xy_F2yz_Py_S = I_ERI_H3xyz_D2y_Py_S+ABZ*I_ERI_G3xy_D2y_Py_S;
  Double I_ERI_G3xz_F2yz_Py_S = I_ERI_H3x2z_D2y_Py_S+ABZ*I_ERI_G3xz_D2y_Py_S;
  Double I_ERI_G2x2y_F2yz_Py_S = I_ERI_H2x2yz_D2y_Py_S+ABZ*I_ERI_G2x2y_D2y_Py_S;
  Double I_ERI_G2xyz_F2yz_Py_S = I_ERI_H2xy2z_D2y_Py_S+ABZ*I_ERI_G2xyz_D2y_Py_S;
  Double I_ERI_G2x2z_F2yz_Py_S = I_ERI_H2x3z_D2y_Py_S+ABZ*I_ERI_G2x2z_D2y_Py_S;
  Double I_ERI_Gx3y_F2yz_Py_S = I_ERI_Hx3yz_D2y_Py_S+ABZ*I_ERI_Gx3y_D2y_Py_S;
  Double I_ERI_Gx2yz_F2yz_Py_S = I_ERI_Hx2y2z_D2y_Py_S+ABZ*I_ERI_Gx2yz_D2y_Py_S;
  Double I_ERI_Gxy2z_F2yz_Py_S = I_ERI_Hxy3z_D2y_Py_S+ABZ*I_ERI_Gxy2z_D2y_Py_S;
  Double I_ERI_Gx3z_F2yz_Py_S = I_ERI_Hx4z_D2y_Py_S+ABZ*I_ERI_Gx3z_D2y_Py_S;
  Double I_ERI_G4y_F2yz_Py_S = I_ERI_H4yz_D2y_Py_S+ABZ*I_ERI_G4y_D2y_Py_S;
  Double I_ERI_G3yz_F2yz_Py_S = I_ERI_H3y2z_D2y_Py_S+ABZ*I_ERI_G3yz_D2y_Py_S;
  Double I_ERI_G2y2z_F2yz_Py_S = I_ERI_H2y3z_D2y_Py_S+ABZ*I_ERI_G2y2z_D2y_Py_S;
  Double I_ERI_Gy3z_F2yz_Py_S = I_ERI_Hy4z_D2y_Py_S+ABZ*I_ERI_Gy3z_D2y_Py_S;
  Double I_ERI_G4z_F2yz_Py_S = I_ERI_H5z_D2y_Py_S+ABZ*I_ERI_G4z_D2y_Py_S;
  Double I_ERI_G4x_F3z_Py_S = I_ERI_H4xz_D2z_Py_S+ABZ*I_ERI_G4x_D2z_Py_S;
  Double I_ERI_G3xy_F3z_Py_S = I_ERI_H3xyz_D2z_Py_S+ABZ*I_ERI_G3xy_D2z_Py_S;
  Double I_ERI_G3xz_F3z_Py_S = I_ERI_H3x2z_D2z_Py_S+ABZ*I_ERI_G3xz_D2z_Py_S;
  Double I_ERI_G2x2y_F3z_Py_S = I_ERI_H2x2yz_D2z_Py_S+ABZ*I_ERI_G2x2y_D2z_Py_S;
  Double I_ERI_G2xyz_F3z_Py_S = I_ERI_H2xy2z_D2z_Py_S+ABZ*I_ERI_G2xyz_D2z_Py_S;
  Double I_ERI_G2x2z_F3z_Py_S = I_ERI_H2x3z_D2z_Py_S+ABZ*I_ERI_G2x2z_D2z_Py_S;
  Double I_ERI_Gx3y_F3z_Py_S = I_ERI_Hx3yz_D2z_Py_S+ABZ*I_ERI_Gx3y_D2z_Py_S;
  Double I_ERI_Gx2yz_F3z_Py_S = I_ERI_Hx2y2z_D2z_Py_S+ABZ*I_ERI_Gx2yz_D2z_Py_S;
  Double I_ERI_Gxy2z_F3z_Py_S = I_ERI_Hxy3z_D2z_Py_S+ABZ*I_ERI_Gxy2z_D2z_Py_S;
  Double I_ERI_Gx3z_F3z_Py_S = I_ERI_Hx4z_D2z_Py_S+ABZ*I_ERI_Gx3z_D2z_Py_S;
  Double I_ERI_G4y_F3z_Py_S = I_ERI_H4yz_D2z_Py_S+ABZ*I_ERI_G4y_D2z_Py_S;
  Double I_ERI_G3yz_F3z_Py_S = I_ERI_H3y2z_D2z_Py_S+ABZ*I_ERI_G3yz_D2z_Py_S;
  Double I_ERI_G2y2z_F3z_Py_S = I_ERI_H2y3z_D2z_Py_S+ABZ*I_ERI_G2y2z_D2z_Py_S;
  Double I_ERI_Gy3z_F3z_Py_S = I_ERI_Hy4z_D2z_Py_S+ABZ*I_ERI_Gy3z_D2z_Py_S;
  Double I_ERI_G4z_F3z_Py_S = I_ERI_H5z_D2z_Py_S+ABZ*I_ERI_G4z_D2z_Py_S;
  Double I_ERI_G4x_F3x_Pz_S = I_ERI_H5x_D2x_Pz_S+ABX*I_ERI_G4x_D2x_Pz_S;
  Double I_ERI_G3xy_F3x_Pz_S = I_ERI_H4xy_D2x_Pz_S+ABX*I_ERI_G3xy_D2x_Pz_S;
  Double I_ERI_G3xz_F3x_Pz_S = I_ERI_H4xz_D2x_Pz_S+ABX*I_ERI_G3xz_D2x_Pz_S;
  Double I_ERI_G2x2y_F3x_Pz_S = I_ERI_H3x2y_D2x_Pz_S+ABX*I_ERI_G2x2y_D2x_Pz_S;
  Double I_ERI_G2xyz_F3x_Pz_S = I_ERI_H3xyz_D2x_Pz_S+ABX*I_ERI_G2xyz_D2x_Pz_S;
  Double I_ERI_G2x2z_F3x_Pz_S = I_ERI_H3x2z_D2x_Pz_S+ABX*I_ERI_G2x2z_D2x_Pz_S;
  Double I_ERI_Gx3y_F3x_Pz_S = I_ERI_H2x3y_D2x_Pz_S+ABX*I_ERI_Gx3y_D2x_Pz_S;
  Double I_ERI_Gx2yz_F3x_Pz_S = I_ERI_H2x2yz_D2x_Pz_S+ABX*I_ERI_Gx2yz_D2x_Pz_S;
  Double I_ERI_Gxy2z_F3x_Pz_S = I_ERI_H2xy2z_D2x_Pz_S+ABX*I_ERI_Gxy2z_D2x_Pz_S;
  Double I_ERI_Gx3z_F3x_Pz_S = I_ERI_H2x3z_D2x_Pz_S+ABX*I_ERI_Gx3z_D2x_Pz_S;
  Double I_ERI_G4y_F3x_Pz_S = I_ERI_Hx4y_D2x_Pz_S+ABX*I_ERI_G4y_D2x_Pz_S;
  Double I_ERI_G3yz_F3x_Pz_S = I_ERI_Hx3yz_D2x_Pz_S+ABX*I_ERI_G3yz_D2x_Pz_S;
  Double I_ERI_G2y2z_F3x_Pz_S = I_ERI_Hx2y2z_D2x_Pz_S+ABX*I_ERI_G2y2z_D2x_Pz_S;
  Double I_ERI_Gy3z_F3x_Pz_S = I_ERI_Hxy3z_D2x_Pz_S+ABX*I_ERI_Gy3z_D2x_Pz_S;
  Double I_ERI_G4z_F3x_Pz_S = I_ERI_Hx4z_D2x_Pz_S+ABX*I_ERI_G4z_D2x_Pz_S;
  Double I_ERI_G4x_F2xy_Pz_S = I_ERI_H4xy_D2x_Pz_S+ABY*I_ERI_G4x_D2x_Pz_S;
  Double I_ERI_G3xy_F2xy_Pz_S = I_ERI_H3x2y_D2x_Pz_S+ABY*I_ERI_G3xy_D2x_Pz_S;
  Double I_ERI_G3xz_F2xy_Pz_S = I_ERI_H3xyz_D2x_Pz_S+ABY*I_ERI_G3xz_D2x_Pz_S;
  Double I_ERI_G2x2y_F2xy_Pz_S = I_ERI_H2x3y_D2x_Pz_S+ABY*I_ERI_G2x2y_D2x_Pz_S;
  Double I_ERI_G2xyz_F2xy_Pz_S = I_ERI_H2x2yz_D2x_Pz_S+ABY*I_ERI_G2xyz_D2x_Pz_S;
  Double I_ERI_G2x2z_F2xy_Pz_S = I_ERI_H2xy2z_D2x_Pz_S+ABY*I_ERI_G2x2z_D2x_Pz_S;
  Double I_ERI_Gx3y_F2xy_Pz_S = I_ERI_Hx4y_D2x_Pz_S+ABY*I_ERI_Gx3y_D2x_Pz_S;
  Double I_ERI_Gx2yz_F2xy_Pz_S = I_ERI_Hx3yz_D2x_Pz_S+ABY*I_ERI_Gx2yz_D2x_Pz_S;
  Double I_ERI_Gxy2z_F2xy_Pz_S = I_ERI_Hx2y2z_D2x_Pz_S+ABY*I_ERI_Gxy2z_D2x_Pz_S;
  Double I_ERI_Gx3z_F2xy_Pz_S = I_ERI_Hxy3z_D2x_Pz_S+ABY*I_ERI_Gx3z_D2x_Pz_S;
  Double I_ERI_G4y_F2xy_Pz_S = I_ERI_H5y_D2x_Pz_S+ABY*I_ERI_G4y_D2x_Pz_S;
  Double I_ERI_G3yz_F2xy_Pz_S = I_ERI_H4yz_D2x_Pz_S+ABY*I_ERI_G3yz_D2x_Pz_S;
  Double I_ERI_G2y2z_F2xy_Pz_S = I_ERI_H3y2z_D2x_Pz_S+ABY*I_ERI_G2y2z_D2x_Pz_S;
  Double I_ERI_Gy3z_F2xy_Pz_S = I_ERI_H2y3z_D2x_Pz_S+ABY*I_ERI_Gy3z_D2x_Pz_S;
  Double I_ERI_G4z_F2xy_Pz_S = I_ERI_Hy4z_D2x_Pz_S+ABY*I_ERI_G4z_D2x_Pz_S;
  Double I_ERI_G4x_F2xz_Pz_S = I_ERI_H4xz_D2x_Pz_S+ABZ*I_ERI_G4x_D2x_Pz_S;
  Double I_ERI_G3xy_F2xz_Pz_S = I_ERI_H3xyz_D2x_Pz_S+ABZ*I_ERI_G3xy_D2x_Pz_S;
  Double I_ERI_G3xz_F2xz_Pz_S = I_ERI_H3x2z_D2x_Pz_S+ABZ*I_ERI_G3xz_D2x_Pz_S;
  Double I_ERI_G2x2y_F2xz_Pz_S = I_ERI_H2x2yz_D2x_Pz_S+ABZ*I_ERI_G2x2y_D2x_Pz_S;
  Double I_ERI_G2xyz_F2xz_Pz_S = I_ERI_H2xy2z_D2x_Pz_S+ABZ*I_ERI_G2xyz_D2x_Pz_S;
  Double I_ERI_G2x2z_F2xz_Pz_S = I_ERI_H2x3z_D2x_Pz_S+ABZ*I_ERI_G2x2z_D2x_Pz_S;
  Double I_ERI_Gx3y_F2xz_Pz_S = I_ERI_Hx3yz_D2x_Pz_S+ABZ*I_ERI_Gx3y_D2x_Pz_S;
  Double I_ERI_Gx2yz_F2xz_Pz_S = I_ERI_Hx2y2z_D2x_Pz_S+ABZ*I_ERI_Gx2yz_D2x_Pz_S;
  Double I_ERI_Gxy2z_F2xz_Pz_S = I_ERI_Hxy3z_D2x_Pz_S+ABZ*I_ERI_Gxy2z_D2x_Pz_S;
  Double I_ERI_Gx3z_F2xz_Pz_S = I_ERI_Hx4z_D2x_Pz_S+ABZ*I_ERI_Gx3z_D2x_Pz_S;
  Double I_ERI_G4y_F2xz_Pz_S = I_ERI_H4yz_D2x_Pz_S+ABZ*I_ERI_G4y_D2x_Pz_S;
  Double I_ERI_G3yz_F2xz_Pz_S = I_ERI_H3y2z_D2x_Pz_S+ABZ*I_ERI_G3yz_D2x_Pz_S;
  Double I_ERI_G2y2z_F2xz_Pz_S = I_ERI_H2y3z_D2x_Pz_S+ABZ*I_ERI_G2y2z_D2x_Pz_S;
  Double I_ERI_Gy3z_F2xz_Pz_S = I_ERI_Hy4z_D2x_Pz_S+ABZ*I_ERI_Gy3z_D2x_Pz_S;
  Double I_ERI_G4z_F2xz_Pz_S = I_ERI_H5z_D2x_Pz_S+ABZ*I_ERI_G4z_D2x_Pz_S;
  Double I_ERI_G4x_Fx2y_Pz_S = I_ERI_H5x_D2y_Pz_S+ABX*I_ERI_G4x_D2y_Pz_S;
  Double I_ERI_G3xy_Fx2y_Pz_S = I_ERI_H4xy_D2y_Pz_S+ABX*I_ERI_G3xy_D2y_Pz_S;
  Double I_ERI_G3xz_Fx2y_Pz_S = I_ERI_H4xz_D2y_Pz_S+ABX*I_ERI_G3xz_D2y_Pz_S;
  Double I_ERI_G2x2y_Fx2y_Pz_S = I_ERI_H3x2y_D2y_Pz_S+ABX*I_ERI_G2x2y_D2y_Pz_S;
  Double I_ERI_G2xyz_Fx2y_Pz_S = I_ERI_H3xyz_D2y_Pz_S+ABX*I_ERI_G2xyz_D2y_Pz_S;
  Double I_ERI_G2x2z_Fx2y_Pz_S = I_ERI_H3x2z_D2y_Pz_S+ABX*I_ERI_G2x2z_D2y_Pz_S;
  Double I_ERI_Gx3y_Fx2y_Pz_S = I_ERI_H2x3y_D2y_Pz_S+ABX*I_ERI_Gx3y_D2y_Pz_S;
  Double I_ERI_Gx2yz_Fx2y_Pz_S = I_ERI_H2x2yz_D2y_Pz_S+ABX*I_ERI_Gx2yz_D2y_Pz_S;
  Double I_ERI_Gxy2z_Fx2y_Pz_S = I_ERI_H2xy2z_D2y_Pz_S+ABX*I_ERI_Gxy2z_D2y_Pz_S;
  Double I_ERI_Gx3z_Fx2y_Pz_S = I_ERI_H2x3z_D2y_Pz_S+ABX*I_ERI_Gx3z_D2y_Pz_S;
  Double I_ERI_G4y_Fx2y_Pz_S = I_ERI_Hx4y_D2y_Pz_S+ABX*I_ERI_G4y_D2y_Pz_S;
  Double I_ERI_G3yz_Fx2y_Pz_S = I_ERI_Hx3yz_D2y_Pz_S+ABX*I_ERI_G3yz_D2y_Pz_S;
  Double I_ERI_G2y2z_Fx2y_Pz_S = I_ERI_Hx2y2z_D2y_Pz_S+ABX*I_ERI_G2y2z_D2y_Pz_S;
  Double I_ERI_Gy3z_Fx2y_Pz_S = I_ERI_Hxy3z_D2y_Pz_S+ABX*I_ERI_Gy3z_D2y_Pz_S;
  Double I_ERI_G4z_Fx2y_Pz_S = I_ERI_Hx4z_D2y_Pz_S+ABX*I_ERI_G4z_D2y_Pz_S;
  Double I_ERI_G4x_Fx2z_Pz_S = I_ERI_H5x_D2z_Pz_S+ABX*I_ERI_G4x_D2z_Pz_S;
  Double I_ERI_G3xy_Fx2z_Pz_S = I_ERI_H4xy_D2z_Pz_S+ABX*I_ERI_G3xy_D2z_Pz_S;
  Double I_ERI_G3xz_Fx2z_Pz_S = I_ERI_H4xz_D2z_Pz_S+ABX*I_ERI_G3xz_D2z_Pz_S;
  Double I_ERI_G2x2y_Fx2z_Pz_S = I_ERI_H3x2y_D2z_Pz_S+ABX*I_ERI_G2x2y_D2z_Pz_S;
  Double I_ERI_G2xyz_Fx2z_Pz_S = I_ERI_H3xyz_D2z_Pz_S+ABX*I_ERI_G2xyz_D2z_Pz_S;
  Double I_ERI_G2x2z_Fx2z_Pz_S = I_ERI_H3x2z_D2z_Pz_S+ABX*I_ERI_G2x2z_D2z_Pz_S;
  Double I_ERI_Gx3y_Fx2z_Pz_S = I_ERI_H2x3y_D2z_Pz_S+ABX*I_ERI_Gx3y_D2z_Pz_S;
  Double I_ERI_Gx2yz_Fx2z_Pz_S = I_ERI_H2x2yz_D2z_Pz_S+ABX*I_ERI_Gx2yz_D2z_Pz_S;
  Double I_ERI_Gxy2z_Fx2z_Pz_S = I_ERI_H2xy2z_D2z_Pz_S+ABX*I_ERI_Gxy2z_D2z_Pz_S;
  Double I_ERI_Gx3z_Fx2z_Pz_S = I_ERI_H2x3z_D2z_Pz_S+ABX*I_ERI_Gx3z_D2z_Pz_S;
  Double I_ERI_G4y_Fx2z_Pz_S = I_ERI_Hx4y_D2z_Pz_S+ABX*I_ERI_G4y_D2z_Pz_S;
  Double I_ERI_G3yz_Fx2z_Pz_S = I_ERI_Hx3yz_D2z_Pz_S+ABX*I_ERI_G3yz_D2z_Pz_S;
  Double I_ERI_G2y2z_Fx2z_Pz_S = I_ERI_Hx2y2z_D2z_Pz_S+ABX*I_ERI_G2y2z_D2z_Pz_S;
  Double I_ERI_Gy3z_Fx2z_Pz_S = I_ERI_Hxy3z_D2z_Pz_S+ABX*I_ERI_Gy3z_D2z_Pz_S;
  Double I_ERI_G4z_Fx2z_Pz_S = I_ERI_Hx4z_D2z_Pz_S+ABX*I_ERI_G4z_D2z_Pz_S;
  Double I_ERI_G4x_F3y_Pz_S = I_ERI_H4xy_D2y_Pz_S+ABY*I_ERI_G4x_D2y_Pz_S;
  Double I_ERI_G3xy_F3y_Pz_S = I_ERI_H3x2y_D2y_Pz_S+ABY*I_ERI_G3xy_D2y_Pz_S;
  Double I_ERI_G3xz_F3y_Pz_S = I_ERI_H3xyz_D2y_Pz_S+ABY*I_ERI_G3xz_D2y_Pz_S;
  Double I_ERI_G2x2y_F3y_Pz_S = I_ERI_H2x3y_D2y_Pz_S+ABY*I_ERI_G2x2y_D2y_Pz_S;
  Double I_ERI_G2xyz_F3y_Pz_S = I_ERI_H2x2yz_D2y_Pz_S+ABY*I_ERI_G2xyz_D2y_Pz_S;
  Double I_ERI_G2x2z_F3y_Pz_S = I_ERI_H2xy2z_D2y_Pz_S+ABY*I_ERI_G2x2z_D2y_Pz_S;
  Double I_ERI_Gx3y_F3y_Pz_S = I_ERI_Hx4y_D2y_Pz_S+ABY*I_ERI_Gx3y_D2y_Pz_S;
  Double I_ERI_Gx2yz_F3y_Pz_S = I_ERI_Hx3yz_D2y_Pz_S+ABY*I_ERI_Gx2yz_D2y_Pz_S;
  Double I_ERI_Gxy2z_F3y_Pz_S = I_ERI_Hx2y2z_D2y_Pz_S+ABY*I_ERI_Gxy2z_D2y_Pz_S;
  Double I_ERI_Gx3z_F3y_Pz_S = I_ERI_Hxy3z_D2y_Pz_S+ABY*I_ERI_Gx3z_D2y_Pz_S;
  Double I_ERI_G4y_F3y_Pz_S = I_ERI_H5y_D2y_Pz_S+ABY*I_ERI_G4y_D2y_Pz_S;
  Double I_ERI_G3yz_F3y_Pz_S = I_ERI_H4yz_D2y_Pz_S+ABY*I_ERI_G3yz_D2y_Pz_S;
  Double I_ERI_G2y2z_F3y_Pz_S = I_ERI_H3y2z_D2y_Pz_S+ABY*I_ERI_G2y2z_D2y_Pz_S;
  Double I_ERI_Gy3z_F3y_Pz_S = I_ERI_H2y3z_D2y_Pz_S+ABY*I_ERI_Gy3z_D2y_Pz_S;
  Double I_ERI_G4z_F3y_Pz_S = I_ERI_Hy4z_D2y_Pz_S+ABY*I_ERI_G4z_D2y_Pz_S;
  Double I_ERI_G4x_F2yz_Pz_S = I_ERI_H4xz_D2y_Pz_S+ABZ*I_ERI_G4x_D2y_Pz_S;
  Double I_ERI_G3xy_F2yz_Pz_S = I_ERI_H3xyz_D2y_Pz_S+ABZ*I_ERI_G3xy_D2y_Pz_S;
  Double I_ERI_G3xz_F2yz_Pz_S = I_ERI_H3x2z_D2y_Pz_S+ABZ*I_ERI_G3xz_D2y_Pz_S;
  Double I_ERI_G2x2y_F2yz_Pz_S = I_ERI_H2x2yz_D2y_Pz_S+ABZ*I_ERI_G2x2y_D2y_Pz_S;
  Double I_ERI_G2xyz_F2yz_Pz_S = I_ERI_H2xy2z_D2y_Pz_S+ABZ*I_ERI_G2xyz_D2y_Pz_S;
  Double I_ERI_G2x2z_F2yz_Pz_S = I_ERI_H2x3z_D2y_Pz_S+ABZ*I_ERI_G2x2z_D2y_Pz_S;
  Double I_ERI_Gx3y_F2yz_Pz_S = I_ERI_Hx3yz_D2y_Pz_S+ABZ*I_ERI_Gx3y_D2y_Pz_S;
  Double I_ERI_Gx2yz_F2yz_Pz_S = I_ERI_Hx2y2z_D2y_Pz_S+ABZ*I_ERI_Gx2yz_D2y_Pz_S;
  Double I_ERI_Gxy2z_F2yz_Pz_S = I_ERI_Hxy3z_D2y_Pz_S+ABZ*I_ERI_Gxy2z_D2y_Pz_S;
  Double I_ERI_Gx3z_F2yz_Pz_S = I_ERI_Hx4z_D2y_Pz_S+ABZ*I_ERI_Gx3z_D2y_Pz_S;
  Double I_ERI_G4y_F2yz_Pz_S = I_ERI_H4yz_D2y_Pz_S+ABZ*I_ERI_G4y_D2y_Pz_S;
  Double I_ERI_G3yz_F2yz_Pz_S = I_ERI_H3y2z_D2y_Pz_S+ABZ*I_ERI_G3yz_D2y_Pz_S;
  Double I_ERI_G2y2z_F2yz_Pz_S = I_ERI_H2y3z_D2y_Pz_S+ABZ*I_ERI_G2y2z_D2y_Pz_S;
  Double I_ERI_Gy3z_F2yz_Pz_S = I_ERI_Hy4z_D2y_Pz_S+ABZ*I_ERI_Gy3z_D2y_Pz_S;
  Double I_ERI_G4z_F2yz_Pz_S = I_ERI_H5z_D2y_Pz_S+ABZ*I_ERI_G4z_D2y_Pz_S;
  Double I_ERI_G4x_F3z_Pz_S = I_ERI_H4xz_D2z_Pz_S+ABZ*I_ERI_G4x_D2z_Pz_S;
  Double I_ERI_G3xy_F3z_Pz_S = I_ERI_H3xyz_D2z_Pz_S+ABZ*I_ERI_G3xy_D2z_Pz_S;
  Double I_ERI_G3xz_F3z_Pz_S = I_ERI_H3x2z_D2z_Pz_S+ABZ*I_ERI_G3xz_D2z_Pz_S;
  Double I_ERI_G2x2y_F3z_Pz_S = I_ERI_H2x2yz_D2z_Pz_S+ABZ*I_ERI_G2x2y_D2z_Pz_S;
  Double I_ERI_G2xyz_F3z_Pz_S = I_ERI_H2xy2z_D2z_Pz_S+ABZ*I_ERI_G2xyz_D2z_Pz_S;
  Double I_ERI_G2x2z_F3z_Pz_S = I_ERI_H2x3z_D2z_Pz_S+ABZ*I_ERI_G2x2z_D2z_Pz_S;
  Double I_ERI_Gx3y_F3z_Pz_S = I_ERI_Hx3yz_D2z_Pz_S+ABZ*I_ERI_Gx3y_D2z_Pz_S;
  Double I_ERI_Gx2yz_F3z_Pz_S = I_ERI_Hx2y2z_D2z_Pz_S+ABZ*I_ERI_Gx2yz_D2z_Pz_S;
  Double I_ERI_Gxy2z_F3z_Pz_S = I_ERI_Hxy3z_D2z_Pz_S+ABZ*I_ERI_Gxy2z_D2z_Pz_S;
  Double I_ERI_Gx3z_F3z_Pz_S = I_ERI_Hx4z_D2z_Pz_S+ABZ*I_ERI_Gx3z_D2z_Pz_S;
  Double I_ERI_G4y_F3z_Pz_S = I_ERI_H4yz_D2z_Pz_S+ABZ*I_ERI_G4y_D2z_Pz_S;
  Double I_ERI_G3yz_F3z_Pz_S = I_ERI_H3y2z_D2z_Pz_S+ABZ*I_ERI_G3yz_D2z_Pz_S;
  Double I_ERI_G2y2z_F3z_Pz_S = I_ERI_H2y3z_D2z_Pz_S+ABZ*I_ERI_G2y2z_D2z_Pz_S;
  Double I_ERI_Gy3z_F3z_Pz_S = I_ERI_Hy4z_D2z_Pz_S+ABZ*I_ERI_Gy3z_D2z_Pz_S;
  Double I_ERI_G4z_F3z_Pz_S = I_ERI_H5z_D2z_Pz_S+ABZ*I_ERI_G4z_D2z_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_K_P_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 81 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_L_S_P_S
   * RHS shell quartet name: SQ_ERI_K_S_P_S
   ************************************************************/
  Double I_ERI_K7x_Px_Px_S = I_ERI_L8x_S_Px_S+ABX*I_ERI_K7x_S_Px_S;
  Double I_ERI_K6xy_Px_Px_S = I_ERI_L7xy_S_Px_S+ABX*I_ERI_K6xy_S_Px_S;
  Double I_ERI_K6xz_Px_Px_S = I_ERI_L7xz_S_Px_S+ABX*I_ERI_K6xz_S_Px_S;
  Double I_ERI_K5x2y_Px_Px_S = I_ERI_L6x2y_S_Px_S+ABX*I_ERI_K5x2y_S_Px_S;
  Double I_ERI_K5xyz_Px_Px_S = I_ERI_L6xyz_S_Px_S+ABX*I_ERI_K5xyz_S_Px_S;
  Double I_ERI_K5x2z_Px_Px_S = I_ERI_L6x2z_S_Px_S+ABX*I_ERI_K5x2z_S_Px_S;
  Double I_ERI_K4x3y_Px_Px_S = I_ERI_L5x3y_S_Px_S+ABX*I_ERI_K4x3y_S_Px_S;
  Double I_ERI_K4x2yz_Px_Px_S = I_ERI_L5x2yz_S_Px_S+ABX*I_ERI_K4x2yz_S_Px_S;
  Double I_ERI_K4xy2z_Px_Px_S = I_ERI_L5xy2z_S_Px_S+ABX*I_ERI_K4xy2z_S_Px_S;
  Double I_ERI_K4x3z_Px_Px_S = I_ERI_L5x3z_S_Px_S+ABX*I_ERI_K4x3z_S_Px_S;
  Double I_ERI_K3x4y_Px_Px_S = I_ERI_L4x4y_S_Px_S+ABX*I_ERI_K3x4y_S_Px_S;
  Double I_ERI_K3x3yz_Px_Px_S = I_ERI_L4x3yz_S_Px_S+ABX*I_ERI_K3x3yz_S_Px_S;
  Double I_ERI_K3x2y2z_Px_Px_S = I_ERI_L4x2y2z_S_Px_S+ABX*I_ERI_K3x2y2z_S_Px_S;
  Double I_ERI_K3xy3z_Px_Px_S = I_ERI_L4xy3z_S_Px_S+ABX*I_ERI_K3xy3z_S_Px_S;
  Double I_ERI_K3x4z_Px_Px_S = I_ERI_L4x4z_S_Px_S+ABX*I_ERI_K3x4z_S_Px_S;
  Double I_ERI_K2x5y_Px_Px_S = I_ERI_L3x5y_S_Px_S+ABX*I_ERI_K2x5y_S_Px_S;
  Double I_ERI_K2x4yz_Px_Px_S = I_ERI_L3x4yz_S_Px_S+ABX*I_ERI_K2x4yz_S_Px_S;
  Double I_ERI_K2x3y2z_Px_Px_S = I_ERI_L3x3y2z_S_Px_S+ABX*I_ERI_K2x3y2z_S_Px_S;
  Double I_ERI_K2x2y3z_Px_Px_S = I_ERI_L3x2y3z_S_Px_S+ABX*I_ERI_K2x2y3z_S_Px_S;
  Double I_ERI_K2xy4z_Px_Px_S = I_ERI_L3xy4z_S_Px_S+ABX*I_ERI_K2xy4z_S_Px_S;
  Double I_ERI_K2x5z_Px_Px_S = I_ERI_L3x5z_S_Px_S+ABX*I_ERI_K2x5z_S_Px_S;
  Double I_ERI_Kx6y_Px_Px_S = I_ERI_L2x6y_S_Px_S+ABX*I_ERI_Kx6y_S_Px_S;
  Double I_ERI_Kx5yz_Px_Px_S = I_ERI_L2x5yz_S_Px_S+ABX*I_ERI_Kx5yz_S_Px_S;
  Double I_ERI_Kx4y2z_Px_Px_S = I_ERI_L2x4y2z_S_Px_S+ABX*I_ERI_Kx4y2z_S_Px_S;
  Double I_ERI_Kx3y3z_Px_Px_S = I_ERI_L2x3y3z_S_Px_S+ABX*I_ERI_Kx3y3z_S_Px_S;
  Double I_ERI_Kx2y4z_Px_Px_S = I_ERI_L2x2y4z_S_Px_S+ABX*I_ERI_Kx2y4z_S_Px_S;
  Double I_ERI_Kxy5z_Px_Px_S = I_ERI_L2xy5z_S_Px_S+ABX*I_ERI_Kxy5z_S_Px_S;
  Double I_ERI_Kx6z_Px_Px_S = I_ERI_L2x6z_S_Px_S+ABX*I_ERI_Kx6z_S_Px_S;
  Double I_ERI_K5x2y_Py_Px_S = I_ERI_L5x3y_S_Px_S+ABY*I_ERI_K5x2y_S_Px_S;
  Double I_ERI_K5xyz_Py_Px_S = I_ERI_L5x2yz_S_Px_S+ABY*I_ERI_K5xyz_S_Px_S;
  Double I_ERI_K4x3y_Py_Px_S = I_ERI_L4x4y_S_Px_S+ABY*I_ERI_K4x3y_S_Px_S;
  Double I_ERI_K4x2yz_Py_Px_S = I_ERI_L4x3yz_S_Px_S+ABY*I_ERI_K4x2yz_S_Px_S;
  Double I_ERI_K4xy2z_Py_Px_S = I_ERI_L4x2y2z_S_Px_S+ABY*I_ERI_K4xy2z_S_Px_S;
  Double I_ERI_K3x4y_Py_Px_S = I_ERI_L3x5y_S_Px_S+ABY*I_ERI_K3x4y_S_Px_S;
  Double I_ERI_K3x3yz_Py_Px_S = I_ERI_L3x4yz_S_Px_S+ABY*I_ERI_K3x3yz_S_Px_S;
  Double I_ERI_K3x2y2z_Py_Px_S = I_ERI_L3x3y2z_S_Px_S+ABY*I_ERI_K3x2y2z_S_Px_S;
  Double I_ERI_K3xy3z_Py_Px_S = I_ERI_L3x2y3z_S_Px_S+ABY*I_ERI_K3xy3z_S_Px_S;
  Double I_ERI_K2x5y_Py_Px_S = I_ERI_L2x6y_S_Px_S+ABY*I_ERI_K2x5y_S_Px_S;
  Double I_ERI_K2x4yz_Py_Px_S = I_ERI_L2x5yz_S_Px_S+ABY*I_ERI_K2x4yz_S_Px_S;
  Double I_ERI_K2x3y2z_Py_Px_S = I_ERI_L2x4y2z_S_Px_S+ABY*I_ERI_K2x3y2z_S_Px_S;
  Double I_ERI_K2x2y3z_Py_Px_S = I_ERI_L2x3y3z_S_Px_S+ABY*I_ERI_K2x2y3z_S_Px_S;
  Double I_ERI_K2xy4z_Py_Px_S = I_ERI_L2x2y4z_S_Px_S+ABY*I_ERI_K2xy4z_S_Px_S;
  Double I_ERI_Kx6y_Py_Px_S = I_ERI_Lx7y_S_Px_S+ABY*I_ERI_Kx6y_S_Px_S;
  Double I_ERI_Kx5yz_Py_Px_S = I_ERI_Lx6yz_S_Px_S+ABY*I_ERI_Kx5yz_S_Px_S;
  Double I_ERI_Kx4y2z_Py_Px_S = I_ERI_Lx5y2z_S_Px_S+ABY*I_ERI_Kx4y2z_S_Px_S;
  Double I_ERI_Kx3y3z_Py_Px_S = I_ERI_Lx4y3z_S_Px_S+ABY*I_ERI_Kx3y3z_S_Px_S;
  Double I_ERI_Kx2y4z_Py_Px_S = I_ERI_Lx3y4z_S_Px_S+ABY*I_ERI_Kx2y4z_S_Px_S;
  Double I_ERI_Kxy5z_Py_Px_S = I_ERI_Lx2y5z_S_Px_S+ABY*I_ERI_Kxy5z_S_Px_S;
  Double I_ERI_K7y_Py_Px_S = I_ERI_L8y_S_Px_S+ABY*I_ERI_K7y_S_Px_S;
  Double I_ERI_K6yz_Py_Px_S = I_ERI_L7yz_S_Px_S+ABY*I_ERI_K6yz_S_Px_S;
  Double I_ERI_K5y2z_Py_Px_S = I_ERI_L6y2z_S_Px_S+ABY*I_ERI_K5y2z_S_Px_S;
  Double I_ERI_K4y3z_Py_Px_S = I_ERI_L5y3z_S_Px_S+ABY*I_ERI_K4y3z_S_Px_S;
  Double I_ERI_K3y4z_Py_Px_S = I_ERI_L4y4z_S_Px_S+ABY*I_ERI_K3y4z_S_Px_S;
  Double I_ERI_K2y5z_Py_Px_S = I_ERI_L3y5z_S_Px_S+ABY*I_ERI_K2y5z_S_Px_S;
  Double I_ERI_Ky6z_Py_Px_S = I_ERI_L2y6z_S_Px_S+ABY*I_ERI_Ky6z_S_Px_S;
  Double I_ERI_K5xyz_Pz_Px_S = I_ERI_L5xy2z_S_Px_S+ABZ*I_ERI_K5xyz_S_Px_S;
  Double I_ERI_K5x2z_Pz_Px_S = I_ERI_L5x3z_S_Px_S+ABZ*I_ERI_K5x2z_S_Px_S;
  Double I_ERI_K4x2yz_Pz_Px_S = I_ERI_L4x2y2z_S_Px_S+ABZ*I_ERI_K4x2yz_S_Px_S;
  Double I_ERI_K4xy2z_Pz_Px_S = I_ERI_L4xy3z_S_Px_S+ABZ*I_ERI_K4xy2z_S_Px_S;
  Double I_ERI_K4x3z_Pz_Px_S = I_ERI_L4x4z_S_Px_S+ABZ*I_ERI_K4x3z_S_Px_S;
  Double I_ERI_K3x3yz_Pz_Px_S = I_ERI_L3x3y2z_S_Px_S+ABZ*I_ERI_K3x3yz_S_Px_S;
  Double I_ERI_K3x2y2z_Pz_Px_S = I_ERI_L3x2y3z_S_Px_S+ABZ*I_ERI_K3x2y2z_S_Px_S;
  Double I_ERI_K3xy3z_Pz_Px_S = I_ERI_L3xy4z_S_Px_S+ABZ*I_ERI_K3xy3z_S_Px_S;
  Double I_ERI_K3x4z_Pz_Px_S = I_ERI_L3x5z_S_Px_S+ABZ*I_ERI_K3x4z_S_Px_S;
  Double I_ERI_K2x4yz_Pz_Px_S = I_ERI_L2x4y2z_S_Px_S+ABZ*I_ERI_K2x4yz_S_Px_S;
  Double I_ERI_K2x3y2z_Pz_Px_S = I_ERI_L2x3y3z_S_Px_S+ABZ*I_ERI_K2x3y2z_S_Px_S;
  Double I_ERI_K2x2y3z_Pz_Px_S = I_ERI_L2x2y4z_S_Px_S+ABZ*I_ERI_K2x2y3z_S_Px_S;
  Double I_ERI_K2xy4z_Pz_Px_S = I_ERI_L2xy5z_S_Px_S+ABZ*I_ERI_K2xy4z_S_Px_S;
  Double I_ERI_K2x5z_Pz_Px_S = I_ERI_L2x6z_S_Px_S+ABZ*I_ERI_K2x5z_S_Px_S;
  Double I_ERI_Kx5yz_Pz_Px_S = I_ERI_Lx5y2z_S_Px_S+ABZ*I_ERI_Kx5yz_S_Px_S;
  Double I_ERI_Kx4y2z_Pz_Px_S = I_ERI_Lx4y3z_S_Px_S+ABZ*I_ERI_Kx4y2z_S_Px_S;
  Double I_ERI_Kx3y3z_Pz_Px_S = I_ERI_Lx3y4z_S_Px_S+ABZ*I_ERI_Kx3y3z_S_Px_S;
  Double I_ERI_Kx2y4z_Pz_Px_S = I_ERI_Lx2y5z_S_Px_S+ABZ*I_ERI_Kx2y4z_S_Px_S;
  Double I_ERI_Kxy5z_Pz_Px_S = I_ERI_Lxy6z_S_Px_S+ABZ*I_ERI_Kxy5z_S_Px_S;
  Double I_ERI_Kx6z_Pz_Px_S = I_ERI_Lx7z_S_Px_S+ABZ*I_ERI_Kx6z_S_Px_S;
  Double I_ERI_K5y2z_Pz_Px_S = I_ERI_L5y3z_S_Px_S+ABZ*I_ERI_K5y2z_S_Px_S;
  Double I_ERI_K4y3z_Pz_Px_S = I_ERI_L4y4z_S_Px_S+ABZ*I_ERI_K4y3z_S_Px_S;
  Double I_ERI_K3y4z_Pz_Px_S = I_ERI_L3y5z_S_Px_S+ABZ*I_ERI_K3y4z_S_Px_S;
  Double I_ERI_K2y5z_Pz_Px_S = I_ERI_L2y6z_S_Px_S+ABZ*I_ERI_K2y5z_S_Px_S;
  Double I_ERI_Ky6z_Pz_Px_S = I_ERI_Ly7z_S_Px_S+ABZ*I_ERI_Ky6z_S_Px_S;
  Double I_ERI_K7z_Pz_Px_S = I_ERI_L8z_S_Px_S+ABZ*I_ERI_K7z_S_Px_S;
  Double I_ERI_K7x_Px_Py_S = I_ERI_L8x_S_Py_S+ABX*I_ERI_K7x_S_Py_S;
  Double I_ERI_K6xy_Px_Py_S = I_ERI_L7xy_S_Py_S+ABX*I_ERI_K6xy_S_Py_S;
  Double I_ERI_K6xz_Px_Py_S = I_ERI_L7xz_S_Py_S+ABX*I_ERI_K6xz_S_Py_S;
  Double I_ERI_K5x2y_Px_Py_S = I_ERI_L6x2y_S_Py_S+ABX*I_ERI_K5x2y_S_Py_S;
  Double I_ERI_K5xyz_Px_Py_S = I_ERI_L6xyz_S_Py_S+ABX*I_ERI_K5xyz_S_Py_S;
  Double I_ERI_K5x2z_Px_Py_S = I_ERI_L6x2z_S_Py_S+ABX*I_ERI_K5x2z_S_Py_S;
  Double I_ERI_K4x3y_Px_Py_S = I_ERI_L5x3y_S_Py_S+ABX*I_ERI_K4x3y_S_Py_S;
  Double I_ERI_K4x2yz_Px_Py_S = I_ERI_L5x2yz_S_Py_S+ABX*I_ERI_K4x2yz_S_Py_S;
  Double I_ERI_K4xy2z_Px_Py_S = I_ERI_L5xy2z_S_Py_S+ABX*I_ERI_K4xy2z_S_Py_S;
  Double I_ERI_K4x3z_Px_Py_S = I_ERI_L5x3z_S_Py_S+ABX*I_ERI_K4x3z_S_Py_S;
  Double I_ERI_K3x4y_Px_Py_S = I_ERI_L4x4y_S_Py_S+ABX*I_ERI_K3x4y_S_Py_S;
  Double I_ERI_K3x3yz_Px_Py_S = I_ERI_L4x3yz_S_Py_S+ABX*I_ERI_K3x3yz_S_Py_S;
  Double I_ERI_K3x2y2z_Px_Py_S = I_ERI_L4x2y2z_S_Py_S+ABX*I_ERI_K3x2y2z_S_Py_S;
  Double I_ERI_K3xy3z_Px_Py_S = I_ERI_L4xy3z_S_Py_S+ABX*I_ERI_K3xy3z_S_Py_S;
  Double I_ERI_K3x4z_Px_Py_S = I_ERI_L4x4z_S_Py_S+ABX*I_ERI_K3x4z_S_Py_S;
  Double I_ERI_K2x5y_Px_Py_S = I_ERI_L3x5y_S_Py_S+ABX*I_ERI_K2x5y_S_Py_S;
  Double I_ERI_K2x4yz_Px_Py_S = I_ERI_L3x4yz_S_Py_S+ABX*I_ERI_K2x4yz_S_Py_S;
  Double I_ERI_K2x3y2z_Px_Py_S = I_ERI_L3x3y2z_S_Py_S+ABX*I_ERI_K2x3y2z_S_Py_S;
  Double I_ERI_K2x2y3z_Px_Py_S = I_ERI_L3x2y3z_S_Py_S+ABX*I_ERI_K2x2y3z_S_Py_S;
  Double I_ERI_K2xy4z_Px_Py_S = I_ERI_L3xy4z_S_Py_S+ABX*I_ERI_K2xy4z_S_Py_S;
  Double I_ERI_K2x5z_Px_Py_S = I_ERI_L3x5z_S_Py_S+ABX*I_ERI_K2x5z_S_Py_S;
  Double I_ERI_Kx6y_Px_Py_S = I_ERI_L2x6y_S_Py_S+ABX*I_ERI_Kx6y_S_Py_S;
  Double I_ERI_Kx5yz_Px_Py_S = I_ERI_L2x5yz_S_Py_S+ABX*I_ERI_Kx5yz_S_Py_S;
  Double I_ERI_Kx4y2z_Px_Py_S = I_ERI_L2x4y2z_S_Py_S+ABX*I_ERI_Kx4y2z_S_Py_S;
  Double I_ERI_Kx3y3z_Px_Py_S = I_ERI_L2x3y3z_S_Py_S+ABX*I_ERI_Kx3y3z_S_Py_S;
  Double I_ERI_Kx2y4z_Px_Py_S = I_ERI_L2x2y4z_S_Py_S+ABX*I_ERI_Kx2y4z_S_Py_S;
  Double I_ERI_Kxy5z_Px_Py_S = I_ERI_L2xy5z_S_Py_S+ABX*I_ERI_Kxy5z_S_Py_S;
  Double I_ERI_Kx6z_Px_Py_S = I_ERI_L2x6z_S_Py_S+ABX*I_ERI_Kx6z_S_Py_S;
  Double I_ERI_K5x2y_Py_Py_S = I_ERI_L5x3y_S_Py_S+ABY*I_ERI_K5x2y_S_Py_S;
  Double I_ERI_K5xyz_Py_Py_S = I_ERI_L5x2yz_S_Py_S+ABY*I_ERI_K5xyz_S_Py_S;
  Double I_ERI_K4x3y_Py_Py_S = I_ERI_L4x4y_S_Py_S+ABY*I_ERI_K4x3y_S_Py_S;
  Double I_ERI_K4x2yz_Py_Py_S = I_ERI_L4x3yz_S_Py_S+ABY*I_ERI_K4x2yz_S_Py_S;
  Double I_ERI_K4xy2z_Py_Py_S = I_ERI_L4x2y2z_S_Py_S+ABY*I_ERI_K4xy2z_S_Py_S;
  Double I_ERI_K3x4y_Py_Py_S = I_ERI_L3x5y_S_Py_S+ABY*I_ERI_K3x4y_S_Py_S;
  Double I_ERI_K3x3yz_Py_Py_S = I_ERI_L3x4yz_S_Py_S+ABY*I_ERI_K3x3yz_S_Py_S;
  Double I_ERI_K3x2y2z_Py_Py_S = I_ERI_L3x3y2z_S_Py_S+ABY*I_ERI_K3x2y2z_S_Py_S;
  Double I_ERI_K3xy3z_Py_Py_S = I_ERI_L3x2y3z_S_Py_S+ABY*I_ERI_K3xy3z_S_Py_S;
  Double I_ERI_K2x5y_Py_Py_S = I_ERI_L2x6y_S_Py_S+ABY*I_ERI_K2x5y_S_Py_S;
  Double I_ERI_K2x4yz_Py_Py_S = I_ERI_L2x5yz_S_Py_S+ABY*I_ERI_K2x4yz_S_Py_S;
  Double I_ERI_K2x3y2z_Py_Py_S = I_ERI_L2x4y2z_S_Py_S+ABY*I_ERI_K2x3y2z_S_Py_S;
  Double I_ERI_K2x2y3z_Py_Py_S = I_ERI_L2x3y3z_S_Py_S+ABY*I_ERI_K2x2y3z_S_Py_S;
  Double I_ERI_K2xy4z_Py_Py_S = I_ERI_L2x2y4z_S_Py_S+ABY*I_ERI_K2xy4z_S_Py_S;
  Double I_ERI_Kx6y_Py_Py_S = I_ERI_Lx7y_S_Py_S+ABY*I_ERI_Kx6y_S_Py_S;
  Double I_ERI_Kx5yz_Py_Py_S = I_ERI_Lx6yz_S_Py_S+ABY*I_ERI_Kx5yz_S_Py_S;
  Double I_ERI_Kx4y2z_Py_Py_S = I_ERI_Lx5y2z_S_Py_S+ABY*I_ERI_Kx4y2z_S_Py_S;
  Double I_ERI_Kx3y3z_Py_Py_S = I_ERI_Lx4y3z_S_Py_S+ABY*I_ERI_Kx3y3z_S_Py_S;
  Double I_ERI_Kx2y4z_Py_Py_S = I_ERI_Lx3y4z_S_Py_S+ABY*I_ERI_Kx2y4z_S_Py_S;
  Double I_ERI_Kxy5z_Py_Py_S = I_ERI_Lx2y5z_S_Py_S+ABY*I_ERI_Kxy5z_S_Py_S;
  Double I_ERI_K7y_Py_Py_S = I_ERI_L8y_S_Py_S+ABY*I_ERI_K7y_S_Py_S;
  Double I_ERI_K6yz_Py_Py_S = I_ERI_L7yz_S_Py_S+ABY*I_ERI_K6yz_S_Py_S;
  Double I_ERI_K5y2z_Py_Py_S = I_ERI_L6y2z_S_Py_S+ABY*I_ERI_K5y2z_S_Py_S;
  Double I_ERI_K4y3z_Py_Py_S = I_ERI_L5y3z_S_Py_S+ABY*I_ERI_K4y3z_S_Py_S;
  Double I_ERI_K3y4z_Py_Py_S = I_ERI_L4y4z_S_Py_S+ABY*I_ERI_K3y4z_S_Py_S;
  Double I_ERI_K2y5z_Py_Py_S = I_ERI_L3y5z_S_Py_S+ABY*I_ERI_K2y5z_S_Py_S;
  Double I_ERI_Ky6z_Py_Py_S = I_ERI_L2y6z_S_Py_S+ABY*I_ERI_Ky6z_S_Py_S;
  Double I_ERI_K5xyz_Pz_Py_S = I_ERI_L5xy2z_S_Py_S+ABZ*I_ERI_K5xyz_S_Py_S;
  Double I_ERI_K5x2z_Pz_Py_S = I_ERI_L5x3z_S_Py_S+ABZ*I_ERI_K5x2z_S_Py_S;
  Double I_ERI_K4x2yz_Pz_Py_S = I_ERI_L4x2y2z_S_Py_S+ABZ*I_ERI_K4x2yz_S_Py_S;
  Double I_ERI_K4xy2z_Pz_Py_S = I_ERI_L4xy3z_S_Py_S+ABZ*I_ERI_K4xy2z_S_Py_S;
  Double I_ERI_K4x3z_Pz_Py_S = I_ERI_L4x4z_S_Py_S+ABZ*I_ERI_K4x3z_S_Py_S;
  Double I_ERI_K3x3yz_Pz_Py_S = I_ERI_L3x3y2z_S_Py_S+ABZ*I_ERI_K3x3yz_S_Py_S;
  Double I_ERI_K3x2y2z_Pz_Py_S = I_ERI_L3x2y3z_S_Py_S+ABZ*I_ERI_K3x2y2z_S_Py_S;
  Double I_ERI_K3xy3z_Pz_Py_S = I_ERI_L3xy4z_S_Py_S+ABZ*I_ERI_K3xy3z_S_Py_S;
  Double I_ERI_K3x4z_Pz_Py_S = I_ERI_L3x5z_S_Py_S+ABZ*I_ERI_K3x4z_S_Py_S;
  Double I_ERI_K2x4yz_Pz_Py_S = I_ERI_L2x4y2z_S_Py_S+ABZ*I_ERI_K2x4yz_S_Py_S;
  Double I_ERI_K2x3y2z_Pz_Py_S = I_ERI_L2x3y3z_S_Py_S+ABZ*I_ERI_K2x3y2z_S_Py_S;
  Double I_ERI_K2x2y3z_Pz_Py_S = I_ERI_L2x2y4z_S_Py_S+ABZ*I_ERI_K2x2y3z_S_Py_S;
  Double I_ERI_K2xy4z_Pz_Py_S = I_ERI_L2xy5z_S_Py_S+ABZ*I_ERI_K2xy4z_S_Py_S;
  Double I_ERI_K2x5z_Pz_Py_S = I_ERI_L2x6z_S_Py_S+ABZ*I_ERI_K2x5z_S_Py_S;
  Double I_ERI_Kx5yz_Pz_Py_S = I_ERI_Lx5y2z_S_Py_S+ABZ*I_ERI_Kx5yz_S_Py_S;
  Double I_ERI_Kx4y2z_Pz_Py_S = I_ERI_Lx4y3z_S_Py_S+ABZ*I_ERI_Kx4y2z_S_Py_S;
  Double I_ERI_Kx3y3z_Pz_Py_S = I_ERI_Lx3y4z_S_Py_S+ABZ*I_ERI_Kx3y3z_S_Py_S;
  Double I_ERI_Kx2y4z_Pz_Py_S = I_ERI_Lx2y5z_S_Py_S+ABZ*I_ERI_Kx2y4z_S_Py_S;
  Double I_ERI_Kxy5z_Pz_Py_S = I_ERI_Lxy6z_S_Py_S+ABZ*I_ERI_Kxy5z_S_Py_S;
  Double I_ERI_Kx6z_Pz_Py_S = I_ERI_Lx7z_S_Py_S+ABZ*I_ERI_Kx6z_S_Py_S;
  Double I_ERI_K5y2z_Pz_Py_S = I_ERI_L5y3z_S_Py_S+ABZ*I_ERI_K5y2z_S_Py_S;
  Double I_ERI_K4y3z_Pz_Py_S = I_ERI_L4y4z_S_Py_S+ABZ*I_ERI_K4y3z_S_Py_S;
  Double I_ERI_K3y4z_Pz_Py_S = I_ERI_L3y5z_S_Py_S+ABZ*I_ERI_K3y4z_S_Py_S;
  Double I_ERI_K2y5z_Pz_Py_S = I_ERI_L2y6z_S_Py_S+ABZ*I_ERI_K2y5z_S_Py_S;
  Double I_ERI_Ky6z_Pz_Py_S = I_ERI_Ly7z_S_Py_S+ABZ*I_ERI_Ky6z_S_Py_S;
  Double I_ERI_K7z_Pz_Py_S = I_ERI_L8z_S_Py_S+ABZ*I_ERI_K7z_S_Py_S;
  Double I_ERI_K7x_Px_Pz_S = I_ERI_L8x_S_Pz_S+ABX*I_ERI_K7x_S_Pz_S;
  Double I_ERI_K6xy_Px_Pz_S = I_ERI_L7xy_S_Pz_S+ABX*I_ERI_K6xy_S_Pz_S;
  Double I_ERI_K6xz_Px_Pz_S = I_ERI_L7xz_S_Pz_S+ABX*I_ERI_K6xz_S_Pz_S;
  Double I_ERI_K5x2y_Px_Pz_S = I_ERI_L6x2y_S_Pz_S+ABX*I_ERI_K5x2y_S_Pz_S;
  Double I_ERI_K5xyz_Px_Pz_S = I_ERI_L6xyz_S_Pz_S+ABX*I_ERI_K5xyz_S_Pz_S;
  Double I_ERI_K5x2z_Px_Pz_S = I_ERI_L6x2z_S_Pz_S+ABX*I_ERI_K5x2z_S_Pz_S;
  Double I_ERI_K4x3y_Px_Pz_S = I_ERI_L5x3y_S_Pz_S+ABX*I_ERI_K4x3y_S_Pz_S;
  Double I_ERI_K4x2yz_Px_Pz_S = I_ERI_L5x2yz_S_Pz_S+ABX*I_ERI_K4x2yz_S_Pz_S;
  Double I_ERI_K4xy2z_Px_Pz_S = I_ERI_L5xy2z_S_Pz_S+ABX*I_ERI_K4xy2z_S_Pz_S;
  Double I_ERI_K4x3z_Px_Pz_S = I_ERI_L5x3z_S_Pz_S+ABX*I_ERI_K4x3z_S_Pz_S;
  Double I_ERI_K3x4y_Px_Pz_S = I_ERI_L4x4y_S_Pz_S+ABX*I_ERI_K3x4y_S_Pz_S;
  Double I_ERI_K3x3yz_Px_Pz_S = I_ERI_L4x3yz_S_Pz_S+ABX*I_ERI_K3x3yz_S_Pz_S;
  Double I_ERI_K3x2y2z_Px_Pz_S = I_ERI_L4x2y2z_S_Pz_S+ABX*I_ERI_K3x2y2z_S_Pz_S;
  Double I_ERI_K3xy3z_Px_Pz_S = I_ERI_L4xy3z_S_Pz_S+ABX*I_ERI_K3xy3z_S_Pz_S;
  Double I_ERI_K3x4z_Px_Pz_S = I_ERI_L4x4z_S_Pz_S+ABX*I_ERI_K3x4z_S_Pz_S;
  Double I_ERI_K2x5y_Px_Pz_S = I_ERI_L3x5y_S_Pz_S+ABX*I_ERI_K2x5y_S_Pz_S;
  Double I_ERI_K2x4yz_Px_Pz_S = I_ERI_L3x4yz_S_Pz_S+ABX*I_ERI_K2x4yz_S_Pz_S;
  Double I_ERI_K2x3y2z_Px_Pz_S = I_ERI_L3x3y2z_S_Pz_S+ABX*I_ERI_K2x3y2z_S_Pz_S;
  Double I_ERI_K2x2y3z_Px_Pz_S = I_ERI_L3x2y3z_S_Pz_S+ABX*I_ERI_K2x2y3z_S_Pz_S;
  Double I_ERI_K2xy4z_Px_Pz_S = I_ERI_L3xy4z_S_Pz_S+ABX*I_ERI_K2xy4z_S_Pz_S;
  Double I_ERI_K2x5z_Px_Pz_S = I_ERI_L3x5z_S_Pz_S+ABX*I_ERI_K2x5z_S_Pz_S;
  Double I_ERI_Kx6y_Px_Pz_S = I_ERI_L2x6y_S_Pz_S+ABX*I_ERI_Kx6y_S_Pz_S;
  Double I_ERI_Kx5yz_Px_Pz_S = I_ERI_L2x5yz_S_Pz_S+ABX*I_ERI_Kx5yz_S_Pz_S;
  Double I_ERI_Kx4y2z_Px_Pz_S = I_ERI_L2x4y2z_S_Pz_S+ABX*I_ERI_Kx4y2z_S_Pz_S;
  Double I_ERI_Kx3y3z_Px_Pz_S = I_ERI_L2x3y3z_S_Pz_S+ABX*I_ERI_Kx3y3z_S_Pz_S;
  Double I_ERI_Kx2y4z_Px_Pz_S = I_ERI_L2x2y4z_S_Pz_S+ABX*I_ERI_Kx2y4z_S_Pz_S;
  Double I_ERI_Kxy5z_Px_Pz_S = I_ERI_L2xy5z_S_Pz_S+ABX*I_ERI_Kxy5z_S_Pz_S;
  Double I_ERI_Kx6z_Px_Pz_S = I_ERI_L2x6z_S_Pz_S+ABX*I_ERI_Kx6z_S_Pz_S;
  Double I_ERI_K5x2y_Py_Pz_S = I_ERI_L5x3y_S_Pz_S+ABY*I_ERI_K5x2y_S_Pz_S;
  Double I_ERI_K5xyz_Py_Pz_S = I_ERI_L5x2yz_S_Pz_S+ABY*I_ERI_K5xyz_S_Pz_S;
  Double I_ERI_K4x3y_Py_Pz_S = I_ERI_L4x4y_S_Pz_S+ABY*I_ERI_K4x3y_S_Pz_S;
  Double I_ERI_K4x2yz_Py_Pz_S = I_ERI_L4x3yz_S_Pz_S+ABY*I_ERI_K4x2yz_S_Pz_S;
  Double I_ERI_K4xy2z_Py_Pz_S = I_ERI_L4x2y2z_S_Pz_S+ABY*I_ERI_K4xy2z_S_Pz_S;
  Double I_ERI_K3x4y_Py_Pz_S = I_ERI_L3x5y_S_Pz_S+ABY*I_ERI_K3x4y_S_Pz_S;
  Double I_ERI_K3x3yz_Py_Pz_S = I_ERI_L3x4yz_S_Pz_S+ABY*I_ERI_K3x3yz_S_Pz_S;
  Double I_ERI_K3x2y2z_Py_Pz_S = I_ERI_L3x3y2z_S_Pz_S+ABY*I_ERI_K3x2y2z_S_Pz_S;
  Double I_ERI_K3xy3z_Py_Pz_S = I_ERI_L3x2y3z_S_Pz_S+ABY*I_ERI_K3xy3z_S_Pz_S;
  Double I_ERI_K2x5y_Py_Pz_S = I_ERI_L2x6y_S_Pz_S+ABY*I_ERI_K2x5y_S_Pz_S;
  Double I_ERI_K2x4yz_Py_Pz_S = I_ERI_L2x5yz_S_Pz_S+ABY*I_ERI_K2x4yz_S_Pz_S;
  Double I_ERI_K2x3y2z_Py_Pz_S = I_ERI_L2x4y2z_S_Pz_S+ABY*I_ERI_K2x3y2z_S_Pz_S;
  Double I_ERI_K2x2y3z_Py_Pz_S = I_ERI_L2x3y3z_S_Pz_S+ABY*I_ERI_K2x2y3z_S_Pz_S;
  Double I_ERI_K2xy4z_Py_Pz_S = I_ERI_L2x2y4z_S_Pz_S+ABY*I_ERI_K2xy4z_S_Pz_S;
  Double I_ERI_Kx6y_Py_Pz_S = I_ERI_Lx7y_S_Pz_S+ABY*I_ERI_Kx6y_S_Pz_S;
  Double I_ERI_Kx5yz_Py_Pz_S = I_ERI_Lx6yz_S_Pz_S+ABY*I_ERI_Kx5yz_S_Pz_S;
  Double I_ERI_Kx4y2z_Py_Pz_S = I_ERI_Lx5y2z_S_Pz_S+ABY*I_ERI_Kx4y2z_S_Pz_S;
  Double I_ERI_Kx3y3z_Py_Pz_S = I_ERI_Lx4y3z_S_Pz_S+ABY*I_ERI_Kx3y3z_S_Pz_S;
  Double I_ERI_Kx2y4z_Py_Pz_S = I_ERI_Lx3y4z_S_Pz_S+ABY*I_ERI_Kx2y4z_S_Pz_S;
  Double I_ERI_Kxy5z_Py_Pz_S = I_ERI_Lx2y5z_S_Pz_S+ABY*I_ERI_Kxy5z_S_Pz_S;
  Double I_ERI_K7y_Py_Pz_S = I_ERI_L8y_S_Pz_S+ABY*I_ERI_K7y_S_Pz_S;
  Double I_ERI_K6yz_Py_Pz_S = I_ERI_L7yz_S_Pz_S+ABY*I_ERI_K6yz_S_Pz_S;
  Double I_ERI_K5y2z_Py_Pz_S = I_ERI_L6y2z_S_Pz_S+ABY*I_ERI_K5y2z_S_Pz_S;
  Double I_ERI_K4y3z_Py_Pz_S = I_ERI_L5y3z_S_Pz_S+ABY*I_ERI_K4y3z_S_Pz_S;
  Double I_ERI_K3y4z_Py_Pz_S = I_ERI_L4y4z_S_Pz_S+ABY*I_ERI_K3y4z_S_Pz_S;
  Double I_ERI_K2y5z_Py_Pz_S = I_ERI_L3y5z_S_Pz_S+ABY*I_ERI_K2y5z_S_Pz_S;
  Double I_ERI_Ky6z_Py_Pz_S = I_ERI_L2y6z_S_Pz_S+ABY*I_ERI_Ky6z_S_Pz_S;
  Double I_ERI_K5xyz_Pz_Pz_S = I_ERI_L5xy2z_S_Pz_S+ABZ*I_ERI_K5xyz_S_Pz_S;
  Double I_ERI_K5x2z_Pz_Pz_S = I_ERI_L5x3z_S_Pz_S+ABZ*I_ERI_K5x2z_S_Pz_S;
  Double I_ERI_K4x2yz_Pz_Pz_S = I_ERI_L4x2y2z_S_Pz_S+ABZ*I_ERI_K4x2yz_S_Pz_S;
  Double I_ERI_K4xy2z_Pz_Pz_S = I_ERI_L4xy3z_S_Pz_S+ABZ*I_ERI_K4xy2z_S_Pz_S;
  Double I_ERI_K4x3z_Pz_Pz_S = I_ERI_L4x4z_S_Pz_S+ABZ*I_ERI_K4x3z_S_Pz_S;
  Double I_ERI_K3x3yz_Pz_Pz_S = I_ERI_L3x3y2z_S_Pz_S+ABZ*I_ERI_K3x3yz_S_Pz_S;
  Double I_ERI_K3x2y2z_Pz_Pz_S = I_ERI_L3x2y3z_S_Pz_S+ABZ*I_ERI_K3x2y2z_S_Pz_S;
  Double I_ERI_K3xy3z_Pz_Pz_S = I_ERI_L3xy4z_S_Pz_S+ABZ*I_ERI_K3xy3z_S_Pz_S;
  Double I_ERI_K3x4z_Pz_Pz_S = I_ERI_L3x5z_S_Pz_S+ABZ*I_ERI_K3x4z_S_Pz_S;
  Double I_ERI_K2x4yz_Pz_Pz_S = I_ERI_L2x4y2z_S_Pz_S+ABZ*I_ERI_K2x4yz_S_Pz_S;
  Double I_ERI_K2x3y2z_Pz_Pz_S = I_ERI_L2x3y3z_S_Pz_S+ABZ*I_ERI_K2x3y2z_S_Pz_S;
  Double I_ERI_K2x2y3z_Pz_Pz_S = I_ERI_L2x2y4z_S_Pz_S+ABZ*I_ERI_K2x2y3z_S_Pz_S;
  Double I_ERI_K2xy4z_Pz_Pz_S = I_ERI_L2xy5z_S_Pz_S+ABZ*I_ERI_K2xy4z_S_Pz_S;
  Double I_ERI_K2x5z_Pz_Pz_S = I_ERI_L2x6z_S_Pz_S+ABZ*I_ERI_K2x5z_S_Pz_S;
  Double I_ERI_Kx5yz_Pz_Pz_S = I_ERI_Lx5y2z_S_Pz_S+ABZ*I_ERI_Kx5yz_S_Pz_S;
  Double I_ERI_Kx4y2z_Pz_Pz_S = I_ERI_Lx4y3z_S_Pz_S+ABZ*I_ERI_Kx4y2z_S_Pz_S;
  Double I_ERI_Kx3y3z_Pz_Pz_S = I_ERI_Lx3y4z_S_Pz_S+ABZ*I_ERI_Kx3y3z_S_Pz_S;
  Double I_ERI_Kx2y4z_Pz_Pz_S = I_ERI_Lx2y5z_S_Pz_S+ABZ*I_ERI_Kx2y4z_S_Pz_S;
  Double I_ERI_Kxy5z_Pz_Pz_S = I_ERI_Lxy6z_S_Pz_S+ABZ*I_ERI_Kxy5z_S_Pz_S;
  Double I_ERI_Kx6z_Pz_Pz_S = I_ERI_Lx7z_S_Pz_S+ABZ*I_ERI_Kx6z_S_Pz_S;
  Double I_ERI_K5y2z_Pz_Pz_S = I_ERI_L5y3z_S_Pz_S+ABZ*I_ERI_K5y2z_S_Pz_S;
  Double I_ERI_K4y3z_Pz_Pz_S = I_ERI_L4y4z_S_Pz_S+ABZ*I_ERI_K4y3z_S_Pz_S;
  Double I_ERI_K3y4z_Pz_Pz_S = I_ERI_L3y5z_S_Pz_S+ABZ*I_ERI_K3y4z_S_Pz_S;
  Double I_ERI_K2y5z_Pz_Pz_S = I_ERI_L2y6z_S_Pz_S+ABZ*I_ERI_K2y5z_S_Pz_S;
  Double I_ERI_Ky6z_Pz_Pz_S = I_ERI_Ly7z_S_Pz_S+ABZ*I_ERI_Ky6z_S_Pz_S;
  Double I_ERI_K7z_Pz_Pz_S = I_ERI_L8z_S_Pz_S+ABZ*I_ERI_K7z_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_I_D_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 261 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_K_P_P_S
   * RHS shell quartet name: SQ_ERI_I_P_P_S
   ************************************************************/
  Double I_ERI_I6x_D2x_Px_S = I_ERI_K7x_Px_Px_S+ABX*I_ERI_I6x_Px_Px_S;
  Double I_ERI_I5xy_D2x_Px_S = I_ERI_K6xy_Px_Px_S+ABX*I_ERI_I5xy_Px_Px_S;
  Double I_ERI_I5xz_D2x_Px_S = I_ERI_K6xz_Px_Px_S+ABX*I_ERI_I5xz_Px_Px_S;
  Double I_ERI_I4x2y_D2x_Px_S = I_ERI_K5x2y_Px_Px_S+ABX*I_ERI_I4x2y_Px_Px_S;
  Double I_ERI_I4xyz_D2x_Px_S = I_ERI_K5xyz_Px_Px_S+ABX*I_ERI_I4xyz_Px_Px_S;
  Double I_ERI_I4x2z_D2x_Px_S = I_ERI_K5x2z_Px_Px_S+ABX*I_ERI_I4x2z_Px_Px_S;
  Double I_ERI_I3x3y_D2x_Px_S = I_ERI_K4x3y_Px_Px_S+ABX*I_ERI_I3x3y_Px_Px_S;
  Double I_ERI_I3x2yz_D2x_Px_S = I_ERI_K4x2yz_Px_Px_S+ABX*I_ERI_I3x2yz_Px_Px_S;
  Double I_ERI_I3xy2z_D2x_Px_S = I_ERI_K4xy2z_Px_Px_S+ABX*I_ERI_I3xy2z_Px_Px_S;
  Double I_ERI_I3x3z_D2x_Px_S = I_ERI_K4x3z_Px_Px_S+ABX*I_ERI_I3x3z_Px_Px_S;
  Double I_ERI_I2x4y_D2x_Px_S = I_ERI_K3x4y_Px_Px_S+ABX*I_ERI_I2x4y_Px_Px_S;
  Double I_ERI_I2x3yz_D2x_Px_S = I_ERI_K3x3yz_Px_Px_S+ABX*I_ERI_I2x3yz_Px_Px_S;
  Double I_ERI_I2x2y2z_D2x_Px_S = I_ERI_K3x2y2z_Px_Px_S+ABX*I_ERI_I2x2y2z_Px_Px_S;
  Double I_ERI_I2xy3z_D2x_Px_S = I_ERI_K3xy3z_Px_Px_S+ABX*I_ERI_I2xy3z_Px_Px_S;
  Double I_ERI_I2x4z_D2x_Px_S = I_ERI_K3x4z_Px_Px_S+ABX*I_ERI_I2x4z_Px_Px_S;
  Double I_ERI_Ix5y_D2x_Px_S = I_ERI_K2x5y_Px_Px_S+ABX*I_ERI_Ix5y_Px_Px_S;
  Double I_ERI_Ix4yz_D2x_Px_S = I_ERI_K2x4yz_Px_Px_S+ABX*I_ERI_Ix4yz_Px_Px_S;
  Double I_ERI_Ix3y2z_D2x_Px_S = I_ERI_K2x3y2z_Px_Px_S+ABX*I_ERI_Ix3y2z_Px_Px_S;
  Double I_ERI_Ix2y3z_D2x_Px_S = I_ERI_K2x2y3z_Px_Px_S+ABX*I_ERI_Ix2y3z_Px_Px_S;
  Double I_ERI_Ixy4z_D2x_Px_S = I_ERI_K2xy4z_Px_Px_S+ABX*I_ERI_Ixy4z_Px_Px_S;
  Double I_ERI_Ix5z_D2x_Px_S = I_ERI_K2x5z_Px_Px_S+ABX*I_ERI_Ix5z_Px_Px_S;
  Double I_ERI_I6y_D2x_Px_S = I_ERI_Kx6y_Px_Px_S+ABX*I_ERI_I6y_Px_Px_S;
  Double I_ERI_I5yz_D2x_Px_S = I_ERI_Kx5yz_Px_Px_S+ABX*I_ERI_I5yz_Px_Px_S;
  Double I_ERI_I4y2z_D2x_Px_S = I_ERI_Kx4y2z_Px_Px_S+ABX*I_ERI_I4y2z_Px_Px_S;
  Double I_ERI_I3y3z_D2x_Px_S = I_ERI_Kx3y3z_Px_Px_S+ABX*I_ERI_I3y3z_Px_Px_S;
  Double I_ERI_I2y4z_D2x_Px_S = I_ERI_Kx2y4z_Px_Px_S+ABX*I_ERI_I2y4z_Px_Px_S;
  Double I_ERI_Iy5z_D2x_Px_S = I_ERI_Kxy5z_Px_Px_S+ABX*I_ERI_Iy5z_Px_Px_S;
  Double I_ERI_I6z_D2x_Px_S = I_ERI_Kx6z_Px_Px_S+ABX*I_ERI_I6z_Px_Px_S;
  Double I_ERI_I5xy_D2y_Px_S = I_ERI_K5x2y_Py_Px_S+ABY*I_ERI_I5xy_Py_Px_S;
  Double I_ERI_I5xz_D2y_Px_S = I_ERI_K5xyz_Py_Px_S+ABY*I_ERI_I5xz_Py_Px_S;
  Double I_ERI_I4x2y_D2y_Px_S = I_ERI_K4x3y_Py_Px_S+ABY*I_ERI_I4x2y_Py_Px_S;
  Double I_ERI_I4xyz_D2y_Px_S = I_ERI_K4x2yz_Py_Px_S+ABY*I_ERI_I4xyz_Py_Px_S;
  Double I_ERI_I4x2z_D2y_Px_S = I_ERI_K4xy2z_Py_Px_S+ABY*I_ERI_I4x2z_Py_Px_S;
  Double I_ERI_I3x3y_D2y_Px_S = I_ERI_K3x4y_Py_Px_S+ABY*I_ERI_I3x3y_Py_Px_S;
  Double I_ERI_I3x2yz_D2y_Px_S = I_ERI_K3x3yz_Py_Px_S+ABY*I_ERI_I3x2yz_Py_Px_S;
  Double I_ERI_I3xy2z_D2y_Px_S = I_ERI_K3x2y2z_Py_Px_S+ABY*I_ERI_I3xy2z_Py_Px_S;
  Double I_ERI_I3x3z_D2y_Px_S = I_ERI_K3xy3z_Py_Px_S+ABY*I_ERI_I3x3z_Py_Px_S;
  Double I_ERI_I2x4y_D2y_Px_S = I_ERI_K2x5y_Py_Px_S+ABY*I_ERI_I2x4y_Py_Px_S;
  Double I_ERI_I2x3yz_D2y_Px_S = I_ERI_K2x4yz_Py_Px_S+ABY*I_ERI_I2x3yz_Py_Px_S;
  Double I_ERI_I2x2y2z_D2y_Px_S = I_ERI_K2x3y2z_Py_Px_S+ABY*I_ERI_I2x2y2z_Py_Px_S;
  Double I_ERI_I2xy3z_D2y_Px_S = I_ERI_K2x2y3z_Py_Px_S+ABY*I_ERI_I2xy3z_Py_Px_S;
  Double I_ERI_I2x4z_D2y_Px_S = I_ERI_K2xy4z_Py_Px_S+ABY*I_ERI_I2x4z_Py_Px_S;
  Double I_ERI_Ix5y_D2y_Px_S = I_ERI_Kx6y_Py_Px_S+ABY*I_ERI_Ix5y_Py_Px_S;
  Double I_ERI_Ix4yz_D2y_Px_S = I_ERI_Kx5yz_Py_Px_S+ABY*I_ERI_Ix4yz_Py_Px_S;
  Double I_ERI_Ix3y2z_D2y_Px_S = I_ERI_Kx4y2z_Py_Px_S+ABY*I_ERI_Ix3y2z_Py_Px_S;
  Double I_ERI_Ix2y3z_D2y_Px_S = I_ERI_Kx3y3z_Py_Px_S+ABY*I_ERI_Ix2y3z_Py_Px_S;
  Double I_ERI_Ixy4z_D2y_Px_S = I_ERI_Kx2y4z_Py_Px_S+ABY*I_ERI_Ixy4z_Py_Px_S;
  Double I_ERI_Ix5z_D2y_Px_S = I_ERI_Kxy5z_Py_Px_S+ABY*I_ERI_Ix5z_Py_Px_S;
  Double I_ERI_I6y_D2y_Px_S = I_ERI_K7y_Py_Px_S+ABY*I_ERI_I6y_Py_Px_S;
  Double I_ERI_I5yz_D2y_Px_S = I_ERI_K6yz_Py_Px_S+ABY*I_ERI_I5yz_Py_Px_S;
  Double I_ERI_I4y2z_D2y_Px_S = I_ERI_K5y2z_Py_Px_S+ABY*I_ERI_I4y2z_Py_Px_S;
  Double I_ERI_I3y3z_D2y_Px_S = I_ERI_K4y3z_Py_Px_S+ABY*I_ERI_I3y3z_Py_Px_S;
  Double I_ERI_I2y4z_D2y_Px_S = I_ERI_K3y4z_Py_Px_S+ABY*I_ERI_I2y4z_Py_Px_S;
  Double I_ERI_Iy5z_D2y_Px_S = I_ERI_K2y5z_Py_Px_S+ABY*I_ERI_Iy5z_Py_Px_S;
  Double I_ERI_I6z_D2y_Px_S = I_ERI_Ky6z_Py_Px_S+ABY*I_ERI_I6z_Py_Px_S;
  Double I_ERI_I5xy_D2z_Px_S = I_ERI_K5xyz_Pz_Px_S+ABZ*I_ERI_I5xy_Pz_Px_S;
  Double I_ERI_I5xz_D2z_Px_S = I_ERI_K5x2z_Pz_Px_S+ABZ*I_ERI_I5xz_Pz_Px_S;
  Double I_ERI_I4x2y_D2z_Px_S = I_ERI_K4x2yz_Pz_Px_S+ABZ*I_ERI_I4x2y_Pz_Px_S;
  Double I_ERI_I4xyz_D2z_Px_S = I_ERI_K4xy2z_Pz_Px_S+ABZ*I_ERI_I4xyz_Pz_Px_S;
  Double I_ERI_I4x2z_D2z_Px_S = I_ERI_K4x3z_Pz_Px_S+ABZ*I_ERI_I4x2z_Pz_Px_S;
  Double I_ERI_I3x3y_D2z_Px_S = I_ERI_K3x3yz_Pz_Px_S+ABZ*I_ERI_I3x3y_Pz_Px_S;
  Double I_ERI_I3x2yz_D2z_Px_S = I_ERI_K3x2y2z_Pz_Px_S+ABZ*I_ERI_I3x2yz_Pz_Px_S;
  Double I_ERI_I3xy2z_D2z_Px_S = I_ERI_K3xy3z_Pz_Px_S+ABZ*I_ERI_I3xy2z_Pz_Px_S;
  Double I_ERI_I3x3z_D2z_Px_S = I_ERI_K3x4z_Pz_Px_S+ABZ*I_ERI_I3x3z_Pz_Px_S;
  Double I_ERI_I2x4y_D2z_Px_S = I_ERI_K2x4yz_Pz_Px_S+ABZ*I_ERI_I2x4y_Pz_Px_S;
  Double I_ERI_I2x3yz_D2z_Px_S = I_ERI_K2x3y2z_Pz_Px_S+ABZ*I_ERI_I2x3yz_Pz_Px_S;
  Double I_ERI_I2x2y2z_D2z_Px_S = I_ERI_K2x2y3z_Pz_Px_S+ABZ*I_ERI_I2x2y2z_Pz_Px_S;
  Double I_ERI_I2xy3z_D2z_Px_S = I_ERI_K2xy4z_Pz_Px_S+ABZ*I_ERI_I2xy3z_Pz_Px_S;
  Double I_ERI_I2x4z_D2z_Px_S = I_ERI_K2x5z_Pz_Px_S+ABZ*I_ERI_I2x4z_Pz_Px_S;
  Double I_ERI_Ix5y_D2z_Px_S = I_ERI_Kx5yz_Pz_Px_S+ABZ*I_ERI_Ix5y_Pz_Px_S;
  Double I_ERI_Ix4yz_D2z_Px_S = I_ERI_Kx4y2z_Pz_Px_S+ABZ*I_ERI_Ix4yz_Pz_Px_S;
  Double I_ERI_Ix3y2z_D2z_Px_S = I_ERI_Kx3y3z_Pz_Px_S+ABZ*I_ERI_Ix3y2z_Pz_Px_S;
  Double I_ERI_Ix2y3z_D2z_Px_S = I_ERI_Kx2y4z_Pz_Px_S+ABZ*I_ERI_Ix2y3z_Pz_Px_S;
  Double I_ERI_Ixy4z_D2z_Px_S = I_ERI_Kxy5z_Pz_Px_S+ABZ*I_ERI_Ixy4z_Pz_Px_S;
  Double I_ERI_Ix5z_D2z_Px_S = I_ERI_Kx6z_Pz_Px_S+ABZ*I_ERI_Ix5z_Pz_Px_S;
  Double I_ERI_I5yz_D2z_Px_S = I_ERI_K5y2z_Pz_Px_S+ABZ*I_ERI_I5yz_Pz_Px_S;
  Double I_ERI_I4y2z_D2z_Px_S = I_ERI_K4y3z_Pz_Px_S+ABZ*I_ERI_I4y2z_Pz_Px_S;
  Double I_ERI_I3y3z_D2z_Px_S = I_ERI_K3y4z_Pz_Px_S+ABZ*I_ERI_I3y3z_Pz_Px_S;
  Double I_ERI_I2y4z_D2z_Px_S = I_ERI_K2y5z_Pz_Px_S+ABZ*I_ERI_I2y4z_Pz_Px_S;
  Double I_ERI_Iy5z_D2z_Px_S = I_ERI_Ky6z_Pz_Px_S+ABZ*I_ERI_Iy5z_Pz_Px_S;
  Double I_ERI_I6z_D2z_Px_S = I_ERI_K7z_Pz_Px_S+ABZ*I_ERI_I6z_Pz_Px_S;
  Double I_ERI_I6x_D2x_Py_S = I_ERI_K7x_Px_Py_S+ABX*I_ERI_I6x_Px_Py_S;
  Double I_ERI_I5xy_D2x_Py_S = I_ERI_K6xy_Px_Py_S+ABX*I_ERI_I5xy_Px_Py_S;
  Double I_ERI_I5xz_D2x_Py_S = I_ERI_K6xz_Px_Py_S+ABX*I_ERI_I5xz_Px_Py_S;
  Double I_ERI_I4x2y_D2x_Py_S = I_ERI_K5x2y_Px_Py_S+ABX*I_ERI_I4x2y_Px_Py_S;
  Double I_ERI_I4xyz_D2x_Py_S = I_ERI_K5xyz_Px_Py_S+ABX*I_ERI_I4xyz_Px_Py_S;
  Double I_ERI_I4x2z_D2x_Py_S = I_ERI_K5x2z_Px_Py_S+ABX*I_ERI_I4x2z_Px_Py_S;
  Double I_ERI_I3x3y_D2x_Py_S = I_ERI_K4x3y_Px_Py_S+ABX*I_ERI_I3x3y_Px_Py_S;
  Double I_ERI_I3x2yz_D2x_Py_S = I_ERI_K4x2yz_Px_Py_S+ABX*I_ERI_I3x2yz_Px_Py_S;
  Double I_ERI_I3xy2z_D2x_Py_S = I_ERI_K4xy2z_Px_Py_S+ABX*I_ERI_I3xy2z_Px_Py_S;
  Double I_ERI_I3x3z_D2x_Py_S = I_ERI_K4x3z_Px_Py_S+ABX*I_ERI_I3x3z_Px_Py_S;
  Double I_ERI_I2x4y_D2x_Py_S = I_ERI_K3x4y_Px_Py_S+ABX*I_ERI_I2x4y_Px_Py_S;
  Double I_ERI_I2x3yz_D2x_Py_S = I_ERI_K3x3yz_Px_Py_S+ABX*I_ERI_I2x3yz_Px_Py_S;
  Double I_ERI_I2x2y2z_D2x_Py_S = I_ERI_K3x2y2z_Px_Py_S+ABX*I_ERI_I2x2y2z_Px_Py_S;
  Double I_ERI_I2xy3z_D2x_Py_S = I_ERI_K3xy3z_Px_Py_S+ABX*I_ERI_I2xy3z_Px_Py_S;
  Double I_ERI_I2x4z_D2x_Py_S = I_ERI_K3x4z_Px_Py_S+ABX*I_ERI_I2x4z_Px_Py_S;
  Double I_ERI_Ix5y_D2x_Py_S = I_ERI_K2x5y_Px_Py_S+ABX*I_ERI_Ix5y_Px_Py_S;
  Double I_ERI_Ix4yz_D2x_Py_S = I_ERI_K2x4yz_Px_Py_S+ABX*I_ERI_Ix4yz_Px_Py_S;
  Double I_ERI_Ix3y2z_D2x_Py_S = I_ERI_K2x3y2z_Px_Py_S+ABX*I_ERI_Ix3y2z_Px_Py_S;
  Double I_ERI_Ix2y3z_D2x_Py_S = I_ERI_K2x2y3z_Px_Py_S+ABX*I_ERI_Ix2y3z_Px_Py_S;
  Double I_ERI_Ixy4z_D2x_Py_S = I_ERI_K2xy4z_Px_Py_S+ABX*I_ERI_Ixy4z_Px_Py_S;
  Double I_ERI_Ix5z_D2x_Py_S = I_ERI_K2x5z_Px_Py_S+ABX*I_ERI_Ix5z_Px_Py_S;
  Double I_ERI_I6y_D2x_Py_S = I_ERI_Kx6y_Px_Py_S+ABX*I_ERI_I6y_Px_Py_S;
  Double I_ERI_I5yz_D2x_Py_S = I_ERI_Kx5yz_Px_Py_S+ABX*I_ERI_I5yz_Px_Py_S;
  Double I_ERI_I4y2z_D2x_Py_S = I_ERI_Kx4y2z_Px_Py_S+ABX*I_ERI_I4y2z_Px_Py_S;
  Double I_ERI_I3y3z_D2x_Py_S = I_ERI_Kx3y3z_Px_Py_S+ABX*I_ERI_I3y3z_Px_Py_S;
  Double I_ERI_I2y4z_D2x_Py_S = I_ERI_Kx2y4z_Px_Py_S+ABX*I_ERI_I2y4z_Px_Py_S;
  Double I_ERI_Iy5z_D2x_Py_S = I_ERI_Kxy5z_Px_Py_S+ABX*I_ERI_Iy5z_Px_Py_S;
  Double I_ERI_I6z_D2x_Py_S = I_ERI_Kx6z_Px_Py_S+ABX*I_ERI_I6z_Px_Py_S;
  Double I_ERI_I5xy_D2y_Py_S = I_ERI_K5x2y_Py_Py_S+ABY*I_ERI_I5xy_Py_Py_S;
  Double I_ERI_I5xz_D2y_Py_S = I_ERI_K5xyz_Py_Py_S+ABY*I_ERI_I5xz_Py_Py_S;
  Double I_ERI_I4x2y_D2y_Py_S = I_ERI_K4x3y_Py_Py_S+ABY*I_ERI_I4x2y_Py_Py_S;
  Double I_ERI_I4xyz_D2y_Py_S = I_ERI_K4x2yz_Py_Py_S+ABY*I_ERI_I4xyz_Py_Py_S;
  Double I_ERI_I4x2z_D2y_Py_S = I_ERI_K4xy2z_Py_Py_S+ABY*I_ERI_I4x2z_Py_Py_S;
  Double I_ERI_I3x3y_D2y_Py_S = I_ERI_K3x4y_Py_Py_S+ABY*I_ERI_I3x3y_Py_Py_S;
  Double I_ERI_I3x2yz_D2y_Py_S = I_ERI_K3x3yz_Py_Py_S+ABY*I_ERI_I3x2yz_Py_Py_S;
  Double I_ERI_I3xy2z_D2y_Py_S = I_ERI_K3x2y2z_Py_Py_S+ABY*I_ERI_I3xy2z_Py_Py_S;
  Double I_ERI_I3x3z_D2y_Py_S = I_ERI_K3xy3z_Py_Py_S+ABY*I_ERI_I3x3z_Py_Py_S;
  Double I_ERI_I2x4y_D2y_Py_S = I_ERI_K2x5y_Py_Py_S+ABY*I_ERI_I2x4y_Py_Py_S;
  Double I_ERI_I2x3yz_D2y_Py_S = I_ERI_K2x4yz_Py_Py_S+ABY*I_ERI_I2x3yz_Py_Py_S;
  Double I_ERI_I2x2y2z_D2y_Py_S = I_ERI_K2x3y2z_Py_Py_S+ABY*I_ERI_I2x2y2z_Py_Py_S;
  Double I_ERI_I2xy3z_D2y_Py_S = I_ERI_K2x2y3z_Py_Py_S+ABY*I_ERI_I2xy3z_Py_Py_S;
  Double I_ERI_I2x4z_D2y_Py_S = I_ERI_K2xy4z_Py_Py_S+ABY*I_ERI_I2x4z_Py_Py_S;
  Double I_ERI_Ix5y_D2y_Py_S = I_ERI_Kx6y_Py_Py_S+ABY*I_ERI_Ix5y_Py_Py_S;
  Double I_ERI_Ix4yz_D2y_Py_S = I_ERI_Kx5yz_Py_Py_S+ABY*I_ERI_Ix4yz_Py_Py_S;
  Double I_ERI_Ix3y2z_D2y_Py_S = I_ERI_Kx4y2z_Py_Py_S+ABY*I_ERI_Ix3y2z_Py_Py_S;
  Double I_ERI_Ix2y3z_D2y_Py_S = I_ERI_Kx3y3z_Py_Py_S+ABY*I_ERI_Ix2y3z_Py_Py_S;
  Double I_ERI_Ixy4z_D2y_Py_S = I_ERI_Kx2y4z_Py_Py_S+ABY*I_ERI_Ixy4z_Py_Py_S;
  Double I_ERI_Ix5z_D2y_Py_S = I_ERI_Kxy5z_Py_Py_S+ABY*I_ERI_Ix5z_Py_Py_S;
  Double I_ERI_I6y_D2y_Py_S = I_ERI_K7y_Py_Py_S+ABY*I_ERI_I6y_Py_Py_S;
  Double I_ERI_I5yz_D2y_Py_S = I_ERI_K6yz_Py_Py_S+ABY*I_ERI_I5yz_Py_Py_S;
  Double I_ERI_I4y2z_D2y_Py_S = I_ERI_K5y2z_Py_Py_S+ABY*I_ERI_I4y2z_Py_Py_S;
  Double I_ERI_I3y3z_D2y_Py_S = I_ERI_K4y3z_Py_Py_S+ABY*I_ERI_I3y3z_Py_Py_S;
  Double I_ERI_I2y4z_D2y_Py_S = I_ERI_K3y4z_Py_Py_S+ABY*I_ERI_I2y4z_Py_Py_S;
  Double I_ERI_Iy5z_D2y_Py_S = I_ERI_K2y5z_Py_Py_S+ABY*I_ERI_Iy5z_Py_Py_S;
  Double I_ERI_I6z_D2y_Py_S = I_ERI_Ky6z_Py_Py_S+ABY*I_ERI_I6z_Py_Py_S;
  Double I_ERI_I5xy_D2z_Py_S = I_ERI_K5xyz_Pz_Py_S+ABZ*I_ERI_I5xy_Pz_Py_S;
  Double I_ERI_I5xz_D2z_Py_S = I_ERI_K5x2z_Pz_Py_S+ABZ*I_ERI_I5xz_Pz_Py_S;
  Double I_ERI_I4x2y_D2z_Py_S = I_ERI_K4x2yz_Pz_Py_S+ABZ*I_ERI_I4x2y_Pz_Py_S;
  Double I_ERI_I4xyz_D2z_Py_S = I_ERI_K4xy2z_Pz_Py_S+ABZ*I_ERI_I4xyz_Pz_Py_S;
  Double I_ERI_I4x2z_D2z_Py_S = I_ERI_K4x3z_Pz_Py_S+ABZ*I_ERI_I4x2z_Pz_Py_S;
  Double I_ERI_I3x3y_D2z_Py_S = I_ERI_K3x3yz_Pz_Py_S+ABZ*I_ERI_I3x3y_Pz_Py_S;
  Double I_ERI_I3x2yz_D2z_Py_S = I_ERI_K3x2y2z_Pz_Py_S+ABZ*I_ERI_I3x2yz_Pz_Py_S;
  Double I_ERI_I3xy2z_D2z_Py_S = I_ERI_K3xy3z_Pz_Py_S+ABZ*I_ERI_I3xy2z_Pz_Py_S;
  Double I_ERI_I3x3z_D2z_Py_S = I_ERI_K3x4z_Pz_Py_S+ABZ*I_ERI_I3x3z_Pz_Py_S;
  Double I_ERI_I2x4y_D2z_Py_S = I_ERI_K2x4yz_Pz_Py_S+ABZ*I_ERI_I2x4y_Pz_Py_S;
  Double I_ERI_I2x3yz_D2z_Py_S = I_ERI_K2x3y2z_Pz_Py_S+ABZ*I_ERI_I2x3yz_Pz_Py_S;
  Double I_ERI_I2x2y2z_D2z_Py_S = I_ERI_K2x2y3z_Pz_Py_S+ABZ*I_ERI_I2x2y2z_Pz_Py_S;
  Double I_ERI_I2xy3z_D2z_Py_S = I_ERI_K2xy4z_Pz_Py_S+ABZ*I_ERI_I2xy3z_Pz_Py_S;
  Double I_ERI_I2x4z_D2z_Py_S = I_ERI_K2x5z_Pz_Py_S+ABZ*I_ERI_I2x4z_Pz_Py_S;
  Double I_ERI_Ix5y_D2z_Py_S = I_ERI_Kx5yz_Pz_Py_S+ABZ*I_ERI_Ix5y_Pz_Py_S;
  Double I_ERI_Ix4yz_D2z_Py_S = I_ERI_Kx4y2z_Pz_Py_S+ABZ*I_ERI_Ix4yz_Pz_Py_S;
  Double I_ERI_Ix3y2z_D2z_Py_S = I_ERI_Kx3y3z_Pz_Py_S+ABZ*I_ERI_Ix3y2z_Pz_Py_S;
  Double I_ERI_Ix2y3z_D2z_Py_S = I_ERI_Kx2y4z_Pz_Py_S+ABZ*I_ERI_Ix2y3z_Pz_Py_S;
  Double I_ERI_Ixy4z_D2z_Py_S = I_ERI_Kxy5z_Pz_Py_S+ABZ*I_ERI_Ixy4z_Pz_Py_S;
  Double I_ERI_Ix5z_D2z_Py_S = I_ERI_Kx6z_Pz_Py_S+ABZ*I_ERI_Ix5z_Pz_Py_S;
  Double I_ERI_I5yz_D2z_Py_S = I_ERI_K5y2z_Pz_Py_S+ABZ*I_ERI_I5yz_Pz_Py_S;
  Double I_ERI_I4y2z_D2z_Py_S = I_ERI_K4y3z_Pz_Py_S+ABZ*I_ERI_I4y2z_Pz_Py_S;
  Double I_ERI_I3y3z_D2z_Py_S = I_ERI_K3y4z_Pz_Py_S+ABZ*I_ERI_I3y3z_Pz_Py_S;
  Double I_ERI_I2y4z_D2z_Py_S = I_ERI_K2y5z_Pz_Py_S+ABZ*I_ERI_I2y4z_Pz_Py_S;
  Double I_ERI_Iy5z_D2z_Py_S = I_ERI_Ky6z_Pz_Py_S+ABZ*I_ERI_Iy5z_Pz_Py_S;
  Double I_ERI_I6z_D2z_Py_S = I_ERI_K7z_Pz_Py_S+ABZ*I_ERI_I6z_Pz_Py_S;
  Double I_ERI_I6x_D2x_Pz_S = I_ERI_K7x_Px_Pz_S+ABX*I_ERI_I6x_Px_Pz_S;
  Double I_ERI_I5xy_D2x_Pz_S = I_ERI_K6xy_Px_Pz_S+ABX*I_ERI_I5xy_Px_Pz_S;
  Double I_ERI_I5xz_D2x_Pz_S = I_ERI_K6xz_Px_Pz_S+ABX*I_ERI_I5xz_Px_Pz_S;
  Double I_ERI_I4x2y_D2x_Pz_S = I_ERI_K5x2y_Px_Pz_S+ABX*I_ERI_I4x2y_Px_Pz_S;
  Double I_ERI_I4xyz_D2x_Pz_S = I_ERI_K5xyz_Px_Pz_S+ABX*I_ERI_I4xyz_Px_Pz_S;
  Double I_ERI_I4x2z_D2x_Pz_S = I_ERI_K5x2z_Px_Pz_S+ABX*I_ERI_I4x2z_Px_Pz_S;
  Double I_ERI_I3x3y_D2x_Pz_S = I_ERI_K4x3y_Px_Pz_S+ABX*I_ERI_I3x3y_Px_Pz_S;
  Double I_ERI_I3x2yz_D2x_Pz_S = I_ERI_K4x2yz_Px_Pz_S+ABX*I_ERI_I3x2yz_Px_Pz_S;
  Double I_ERI_I3xy2z_D2x_Pz_S = I_ERI_K4xy2z_Px_Pz_S+ABX*I_ERI_I3xy2z_Px_Pz_S;
  Double I_ERI_I3x3z_D2x_Pz_S = I_ERI_K4x3z_Px_Pz_S+ABX*I_ERI_I3x3z_Px_Pz_S;
  Double I_ERI_I2x4y_D2x_Pz_S = I_ERI_K3x4y_Px_Pz_S+ABX*I_ERI_I2x4y_Px_Pz_S;
  Double I_ERI_I2x3yz_D2x_Pz_S = I_ERI_K3x3yz_Px_Pz_S+ABX*I_ERI_I2x3yz_Px_Pz_S;
  Double I_ERI_I2x2y2z_D2x_Pz_S = I_ERI_K3x2y2z_Px_Pz_S+ABX*I_ERI_I2x2y2z_Px_Pz_S;
  Double I_ERI_I2xy3z_D2x_Pz_S = I_ERI_K3xy3z_Px_Pz_S+ABX*I_ERI_I2xy3z_Px_Pz_S;
  Double I_ERI_I2x4z_D2x_Pz_S = I_ERI_K3x4z_Px_Pz_S+ABX*I_ERI_I2x4z_Px_Pz_S;
  Double I_ERI_Ix5y_D2x_Pz_S = I_ERI_K2x5y_Px_Pz_S+ABX*I_ERI_Ix5y_Px_Pz_S;
  Double I_ERI_Ix4yz_D2x_Pz_S = I_ERI_K2x4yz_Px_Pz_S+ABX*I_ERI_Ix4yz_Px_Pz_S;
  Double I_ERI_Ix3y2z_D2x_Pz_S = I_ERI_K2x3y2z_Px_Pz_S+ABX*I_ERI_Ix3y2z_Px_Pz_S;
  Double I_ERI_Ix2y3z_D2x_Pz_S = I_ERI_K2x2y3z_Px_Pz_S+ABX*I_ERI_Ix2y3z_Px_Pz_S;
  Double I_ERI_Ixy4z_D2x_Pz_S = I_ERI_K2xy4z_Px_Pz_S+ABX*I_ERI_Ixy4z_Px_Pz_S;
  Double I_ERI_Ix5z_D2x_Pz_S = I_ERI_K2x5z_Px_Pz_S+ABX*I_ERI_Ix5z_Px_Pz_S;
  Double I_ERI_I6y_D2x_Pz_S = I_ERI_Kx6y_Px_Pz_S+ABX*I_ERI_I6y_Px_Pz_S;
  Double I_ERI_I5yz_D2x_Pz_S = I_ERI_Kx5yz_Px_Pz_S+ABX*I_ERI_I5yz_Px_Pz_S;
  Double I_ERI_I4y2z_D2x_Pz_S = I_ERI_Kx4y2z_Px_Pz_S+ABX*I_ERI_I4y2z_Px_Pz_S;
  Double I_ERI_I3y3z_D2x_Pz_S = I_ERI_Kx3y3z_Px_Pz_S+ABX*I_ERI_I3y3z_Px_Pz_S;
  Double I_ERI_I2y4z_D2x_Pz_S = I_ERI_Kx2y4z_Px_Pz_S+ABX*I_ERI_I2y4z_Px_Pz_S;
  Double I_ERI_Iy5z_D2x_Pz_S = I_ERI_Kxy5z_Px_Pz_S+ABX*I_ERI_Iy5z_Px_Pz_S;
  Double I_ERI_I6z_D2x_Pz_S = I_ERI_Kx6z_Px_Pz_S+ABX*I_ERI_I6z_Px_Pz_S;
  Double I_ERI_I5xy_D2y_Pz_S = I_ERI_K5x2y_Py_Pz_S+ABY*I_ERI_I5xy_Py_Pz_S;
  Double I_ERI_I5xz_D2y_Pz_S = I_ERI_K5xyz_Py_Pz_S+ABY*I_ERI_I5xz_Py_Pz_S;
  Double I_ERI_I4x2y_D2y_Pz_S = I_ERI_K4x3y_Py_Pz_S+ABY*I_ERI_I4x2y_Py_Pz_S;
  Double I_ERI_I4xyz_D2y_Pz_S = I_ERI_K4x2yz_Py_Pz_S+ABY*I_ERI_I4xyz_Py_Pz_S;
  Double I_ERI_I4x2z_D2y_Pz_S = I_ERI_K4xy2z_Py_Pz_S+ABY*I_ERI_I4x2z_Py_Pz_S;
  Double I_ERI_I3x3y_D2y_Pz_S = I_ERI_K3x4y_Py_Pz_S+ABY*I_ERI_I3x3y_Py_Pz_S;
  Double I_ERI_I3x2yz_D2y_Pz_S = I_ERI_K3x3yz_Py_Pz_S+ABY*I_ERI_I3x2yz_Py_Pz_S;
  Double I_ERI_I3xy2z_D2y_Pz_S = I_ERI_K3x2y2z_Py_Pz_S+ABY*I_ERI_I3xy2z_Py_Pz_S;
  Double I_ERI_I3x3z_D2y_Pz_S = I_ERI_K3xy3z_Py_Pz_S+ABY*I_ERI_I3x3z_Py_Pz_S;
  Double I_ERI_I2x4y_D2y_Pz_S = I_ERI_K2x5y_Py_Pz_S+ABY*I_ERI_I2x4y_Py_Pz_S;
  Double I_ERI_I2x3yz_D2y_Pz_S = I_ERI_K2x4yz_Py_Pz_S+ABY*I_ERI_I2x3yz_Py_Pz_S;
  Double I_ERI_I2x2y2z_D2y_Pz_S = I_ERI_K2x3y2z_Py_Pz_S+ABY*I_ERI_I2x2y2z_Py_Pz_S;
  Double I_ERI_I2xy3z_D2y_Pz_S = I_ERI_K2x2y3z_Py_Pz_S+ABY*I_ERI_I2xy3z_Py_Pz_S;
  Double I_ERI_I2x4z_D2y_Pz_S = I_ERI_K2xy4z_Py_Pz_S+ABY*I_ERI_I2x4z_Py_Pz_S;
  Double I_ERI_Ix5y_D2y_Pz_S = I_ERI_Kx6y_Py_Pz_S+ABY*I_ERI_Ix5y_Py_Pz_S;
  Double I_ERI_Ix4yz_D2y_Pz_S = I_ERI_Kx5yz_Py_Pz_S+ABY*I_ERI_Ix4yz_Py_Pz_S;
  Double I_ERI_Ix3y2z_D2y_Pz_S = I_ERI_Kx4y2z_Py_Pz_S+ABY*I_ERI_Ix3y2z_Py_Pz_S;
  Double I_ERI_Ix2y3z_D2y_Pz_S = I_ERI_Kx3y3z_Py_Pz_S+ABY*I_ERI_Ix2y3z_Py_Pz_S;
  Double I_ERI_Ixy4z_D2y_Pz_S = I_ERI_Kx2y4z_Py_Pz_S+ABY*I_ERI_Ixy4z_Py_Pz_S;
  Double I_ERI_Ix5z_D2y_Pz_S = I_ERI_Kxy5z_Py_Pz_S+ABY*I_ERI_Ix5z_Py_Pz_S;
  Double I_ERI_I6y_D2y_Pz_S = I_ERI_K7y_Py_Pz_S+ABY*I_ERI_I6y_Py_Pz_S;
  Double I_ERI_I5yz_D2y_Pz_S = I_ERI_K6yz_Py_Pz_S+ABY*I_ERI_I5yz_Py_Pz_S;
  Double I_ERI_I4y2z_D2y_Pz_S = I_ERI_K5y2z_Py_Pz_S+ABY*I_ERI_I4y2z_Py_Pz_S;
  Double I_ERI_I3y3z_D2y_Pz_S = I_ERI_K4y3z_Py_Pz_S+ABY*I_ERI_I3y3z_Py_Pz_S;
  Double I_ERI_I2y4z_D2y_Pz_S = I_ERI_K3y4z_Py_Pz_S+ABY*I_ERI_I2y4z_Py_Pz_S;
  Double I_ERI_Iy5z_D2y_Pz_S = I_ERI_K2y5z_Py_Pz_S+ABY*I_ERI_Iy5z_Py_Pz_S;
  Double I_ERI_I6z_D2y_Pz_S = I_ERI_Ky6z_Py_Pz_S+ABY*I_ERI_I6z_Py_Pz_S;
  Double I_ERI_I5xy_D2z_Pz_S = I_ERI_K5xyz_Pz_Pz_S+ABZ*I_ERI_I5xy_Pz_Pz_S;
  Double I_ERI_I5xz_D2z_Pz_S = I_ERI_K5x2z_Pz_Pz_S+ABZ*I_ERI_I5xz_Pz_Pz_S;
  Double I_ERI_I4x2y_D2z_Pz_S = I_ERI_K4x2yz_Pz_Pz_S+ABZ*I_ERI_I4x2y_Pz_Pz_S;
  Double I_ERI_I4xyz_D2z_Pz_S = I_ERI_K4xy2z_Pz_Pz_S+ABZ*I_ERI_I4xyz_Pz_Pz_S;
  Double I_ERI_I4x2z_D2z_Pz_S = I_ERI_K4x3z_Pz_Pz_S+ABZ*I_ERI_I4x2z_Pz_Pz_S;
  Double I_ERI_I3x3y_D2z_Pz_S = I_ERI_K3x3yz_Pz_Pz_S+ABZ*I_ERI_I3x3y_Pz_Pz_S;
  Double I_ERI_I3x2yz_D2z_Pz_S = I_ERI_K3x2y2z_Pz_Pz_S+ABZ*I_ERI_I3x2yz_Pz_Pz_S;
  Double I_ERI_I3xy2z_D2z_Pz_S = I_ERI_K3xy3z_Pz_Pz_S+ABZ*I_ERI_I3xy2z_Pz_Pz_S;
  Double I_ERI_I3x3z_D2z_Pz_S = I_ERI_K3x4z_Pz_Pz_S+ABZ*I_ERI_I3x3z_Pz_Pz_S;
  Double I_ERI_I2x4y_D2z_Pz_S = I_ERI_K2x4yz_Pz_Pz_S+ABZ*I_ERI_I2x4y_Pz_Pz_S;
  Double I_ERI_I2x3yz_D2z_Pz_S = I_ERI_K2x3y2z_Pz_Pz_S+ABZ*I_ERI_I2x3yz_Pz_Pz_S;
  Double I_ERI_I2x2y2z_D2z_Pz_S = I_ERI_K2x2y3z_Pz_Pz_S+ABZ*I_ERI_I2x2y2z_Pz_Pz_S;
  Double I_ERI_I2xy3z_D2z_Pz_S = I_ERI_K2xy4z_Pz_Pz_S+ABZ*I_ERI_I2xy3z_Pz_Pz_S;
  Double I_ERI_I2x4z_D2z_Pz_S = I_ERI_K2x5z_Pz_Pz_S+ABZ*I_ERI_I2x4z_Pz_Pz_S;
  Double I_ERI_Ix5y_D2z_Pz_S = I_ERI_Kx5yz_Pz_Pz_S+ABZ*I_ERI_Ix5y_Pz_Pz_S;
  Double I_ERI_Ix4yz_D2z_Pz_S = I_ERI_Kx4y2z_Pz_Pz_S+ABZ*I_ERI_Ix4yz_Pz_Pz_S;
  Double I_ERI_Ix3y2z_D2z_Pz_S = I_ERI_Kx3y3z_Pz_Pz_S+ABZ*I_ERI_Ix3y2z_Pz_Pz_S;
  Double I_ERI_Ix2y3z_D2z_Pz_S = I_ERI_Kx2y4z_Pz_Pz_S+ABZ*I_ERI_Ix2y3z_Pz_Pz_S;
  Double I_ERI_Ixy4z_D2z_Pz_S = I_ERI_Kxy5z_Pz_Pz_S+ABZ*I_ERI_Ixy4z_Pz_Pz_S;
  Double I_ERI_Ix5z_D2z_Pz_S = I_ERI_Kx6z_Pz_Pz_S+ABZ*I_ERI_Ix5z_Pz_Pz_S;
  Double I_ERI_I5yz_D2z_Pz_S = I_ERI_K5y2z_Pz_Pz_S+ABZ*I_ERI_I5yz_Pz_Pz_S;
  Double I_ERI_I4y2z_D2z_Pz_S = I_ERI_K4y3z_Pz_Pz_S+ABZ*I_ERI_I4y2z_Pz_Pz_S;
  Double I_ERI_I3y3z_D2z_Pz_S = I_ERI_K3y4z_Pz_Pz_S+ABZ*I_ERI_I3y3z_Pz_Pz_S;
  Double I_ERI_I2y4z_D2z_Pz_S = I_ERI_K2y5z_Pz_Pz_S+ABZ*I_ERI_I2y4z_Pz_Pz_S;
  Double I_ERI_Iy5z_D2z_Pz_S = I_ERI_Ky6z_Pz_Pz_S+ABZ*I_ERI_Iy5z_Pz_Pz_S;
  Double I_ERI_I6z_D2z_Pz_S = I_ERI_K7z_Pz_Pz_S+ABZ*I_ERI_I6z_Pz_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_H_F_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 201 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_I_D_P_S
   * RHS shell quartet name: SQ_ERI_H_D_P_S
   ************************************************************/
  Double I_ERI_H5x_F3x_Px_S = I_ERI_I6x_D2x_Px_S+ABX*I_ERI_H5x_D2x_Px_S;
  Double I_ERI_H4xy_F3x_Px_S = I_ERI_I5xy_D2x_Px_S+ABX*I_ERI_H4xy_D2x_Px_S;
  Double I_ERI_H4xz_F3x_Px_S = I_ERI_I5xz_D2x_Px_S+ABX*I_ERI_H4xz_D2x_Px_S;
  Double I_ERI_H3x2y_F3x_Px_S = I_ERI_I4x2y_D2x_Px_S+ABX*I_ERI_H3x2y_D2x_Px_S;
  Double I_ERI_H3xyz_F3x_Px_S = I_ERI_I4xyz_D2x_Px_S+ABX*I_ERI_H3xyz_D2x_Px_S;
  Double I_ERI_H3x2z_F3x_Px_S = I_ERI_I4x2z_D2x_Px_S+ABX*I_ERI_H3x2z_D2x_Px_S;
  Double I_ERI_H2x3y_F3x_Px_S = I_ERI_I3x3y_D2x_Px_S+ABX*I_ERI_H2x3y_D2x_Px_S;
  Double I_ERI_H2x2yz_F3x_Px_S = I_ERI_I3x2yz_D2x_Px_S+ABX*I_ERI_H2x2yz_D2x_Px_S;
  Double I_ERI_H2xy2z_F3x_Px_S = I_ERI_I3xy2z_D2x_Px_S+ABX*I_ERI_H2xy2z_D2x_Px_S;
  Double I_ERI_H2x3z_F3x_Px_S = I_ERI_I3x3z_D2x_Px_S+ABX*I_ERI_H2x3z_D2x_Px_S;
  Double I_ERI_Hx4y_F3x_Px_S = I_ERI_I2x4y_D2x_Px_S+ABX*I_ERI_Hx4y_D2x_Px_S;
  Double I_ERI_Hx3yz_F3x_Px_S = I_ERI_I2x3yz_D2x_Px_S+ABX*I_ERI_Hx3yz_D2x_Px_S;
  Double I_ERI_Hx2y2z_F3x_Px_S = I_ERI_I2x2y2z_D2x_Px_S+ABX*I_ERI_Hx2y2z_D2x_Px_S;
  Double I_ERI_Hxy3z_F3x_Px_S = I_ERI_I2xy3z_D2x_Px_S+ABX*I_ERI_Hxy3z_D2x_Px_S;
  Double I_ERI_Hx4z_F3x_Px_S = I_ERI_I2x4z_D2x_Px_S+ABX*I_ERI_Hx4z_D2x_Px_S;
  Double I_ERI_H5y_F3x_Px_S = I_ERI_Ix5y_D2x_Px_S+ABX*I_ERI_H5y_D2x_Px_S;
  Double I_ERI_H4yz_F3x_Px_S = I_ERI_Ix4yz_D2x_Px_S+ABX*I_ERI_H4yz_D2x_Px_S;
  Double I_ERI_H3y2z_F3x_Px_S = I_ERI_Ix3y2z_D2x_Px_S+ABX*I_ERI_H3y2z_D2x_Px_S;
  Double I_ERI_H2y3z_F3x_Px_S = I_ERI_Ix2y3z_D2x_Px_S+ABX*I_ERI_H2y3z_D2x_Px_S;
  Double I_ERI_Hy4z_F3x_Px_S = I_ERI_Ixy4z_D2x_Px_S+ABX*I_ERI_Hy4z_D2x_Px_S;
  Double I_ERI_H5z_F3x_Px_S = I_ERI_Ix5z_D2x_Px_S+ABX*I_ERI_H5z_D2x_Px_S;
  Double I_ERI_H4xy_F2xy_Px_S = I_ERI_I4x2y_D2x_Px_S+ABY*I_ERI_H4xy_D2x_Px_S;
  Double I_ERI_H4xz_F2xy_Px_S = I_ERI_I4xyz_D2x_Px_S+ABY*I_ERI_H4xz_D2x_Px_S;
  Double I_ERI_H3x2y_F2xy_Px_S = I_ERI_I3x3y_D2x_Px_S+ABY*I_ERI_H3x2y_D2x_Px_S;
  Double I_ERI_H3xyz_F2xy_Px_S = I_ERI_I3x2yz_D2x_Px_S+ABY*I_ERI_H3xyz_D2x_Px_S;
  Double I_ERI_H3x2z_F2xy_Px_S = I_ERI_I3xy2z_D2x_Px_S+ABY*I_ERI_H3x2z_D2x_Px_S;
  Double I_ERI_H2x3y_F2xy_Px_S = I_ERI_I2x4y_D2x_Px_S+ABY*I_ERI_H2x3y_D2x_Px_S;
  Double I_ERI_H2x2yz_F2xy_Px_S = I_ERI_I2x3yz_D2x_Px_S+ABY*I_ERI_H2x2yz_D2x_Px_S;
  Double I_ERI_H2xy2z_F2xy_Px_S = I_ERI_I2x2y2z_D2x_Px_S+ABY*I_ERI_H2xy2z_D2x_Px_S;
  Double I_ERI_H2x3z_F2xy_Px_S = I_ERI_I2xy3z_D2x_Px_S+ABY*I_ERI_H2x3z_D2x_Px_S;
  Double I_ERI_Hx4y_F2xy_Px_S = I_ERI_Ix5y_D2x_Px_S+ABY*I_ERI_Hx4y_D2x_Px_S;
  Double I_ERI_Hx3yz_F2xy_Px_S = I_ERI_Ix4yz_D2x_Px_S+ABY*I_ERI_Hx3yz_D2x_Px_S;
  Double I_ERI_Hx2y2z_F2xy_Px_S = I_ERI_Ix3y2z_D2x_Px_S+ABY*I_ERI_Hx2y2z_D2x_Px_S;
  Double I_ERI_Hxy3z_F2xy_Px_S = I_ERI_Ix2y3z_D2x_Px_S+ABY*I_ERI_Hxy3z_D2x_Px_S;
  Double I_ERI_Hx4z_F2xy_Px_S = I_ERI_Ixy4z_D2x_Px_S+ABY*I_ERI_Hx4z_D2x_Px_S;
  Double I_ERI_H5y_F2xy_Px_S = I_ERI_I6y_D2x_Px_S+ABY*I_ERI_H5y_D2x_Px_S;
  Double I_ERI_H4yz_F2xy_Px_S = I_ERI_I5yz_D2x_Px_S+ABY*I_ERI_H4yz_D2x_Px_S;
  Double I_ERI_H3y2z_F2xy_Px_S = I_ERI_I4y2z_D2x_Px_S+ABY*I_ERI_H3y2z_D2x_Px_S;
  Double I_ERI_H2y3z_F2xy_Px_S = I_ERI_I3y3z_D2x_Px_S+ABY*I_ERI_H2y3z_D2x_Px_S;
  Double I_ERI_Hy4z_F2xy_Px_S = I_ERI_I2y4z_D2x_Px_S+ABY*I_ERI_Hy4z_D2x_Px_S;
  Double I_ERI_H5z_F2xy_Px_S = I_ERI_Iy5z_D2x_Px_S+ABY*I_ERI_H5z_D2x_Px_S;
  Double I_ERI_H4xz_F2xz_Px_S = I_ERI_I4x2z_D2x_Px_S+ABZ*I_ERI_H4xz_D2x_Px_S;
  Double I_ERI_H3xyz_F2xz_Px_S = I_ERI_I3xy2z_D2x_Px_S+ABZ*I_ERI_H3xyz_D2x_Px_S;
  Double I_ERI_H3x2z_F2xz_Px_S = I_ERI_I3x3z_D2x_Px_S+ABZ*I_ERI_H3x2z_D2x_Px_S;
  Double I_ERI_H2x2yz_F2xz_Px_S = I_ERI_I2x2y2z_D2x_Px_S+ABZ*I_ERI_H2x2yz_D2x_Px_S;
  Double I_ERI_H2xy2z_F2xz_Px_S = I_ERI_I2xy3z_D2x_Px_S+ABZ*I_ERI_H2xy2z_D2x_Px_S;
  Double I_ERI_H2x3z_F2xz_Px_S = I_ERI_I2x4z_D2x_Px_S+ABZ*I_ERI_H2x3z_D2x_Px_S;
  Double I_ERI_Hx3yz_F2xz_Px_S = I_ERI_Ix3y2z_D2x_Px_S+ABZ*I_ERI_Hx3yz_D2x_Px_S;
  Double I_ERI_Hx2y2z_F2xz_Px_S = I_ERI_Ix2y3z_D2x_Px_S+ABZ*I_ERI_Hx2y2z_D2x_Px_S;
  Double I_ERI_Hxy3z_F2xz_Px_S = I_ERI_Ixy4z_D2x_Px_S+ABZ*I_ERI_Hxy3z_D2x_Px_S;
  Double I_ERI_Hx4z_F2xz_Px_S = I_ERI_Ix5z_D2x_Px_S+ABZ*I_ERI_Hx4z_D2x_Px_S;
  Double I_ERI_H4yz_F2xz_Px_S = I_ERI_I4y2z_D2x_Px_S+ABZ*I_ERI_H4yz_D2x_Px_S;
  Double I_ERI_H3y2z_F2xz_Px_S = I_ERI_I3y3z_D2x_Px_S+ABZ*I_ERI_H3y2z_D2x_Px_S;
  Double I_ERI_H2y3z_F2xz_Px_S = I_ERI_I2y4z_D2x_Px_S+ABZ*I_ERI_H2y3z_D2x_Px_S;
  Double I_ERI_Hy4z_F2xz_Px_S = I_ERI_Iy5z_D2x_Px_S+ABZ*I_ERI_Hy4z_D2x_Px_S;
  Double I_ERI_H5z_F2xz_Px_S = I_ERI_I6z_D2x_Px_S+ABZ*I_ERI_H5z_D2x_Px_S;
  Double I_ERI_H4xz_Fx2y_Px_S = I_ERI_I5xz_D2y_Px_S+ABX*I_ERI_H4xz_D2y_Px_S;
  Double I_ERI_H3xyz_Fx2y_Px_S = I_ERI_I4xyz_D2y_Px_S+ABX*I_ERI_H3xyz_D2y_Px_S;
  Double I_ERI_H3x2z_Fx2y_Px_S = I_ERI_I4x2z_D2y_Px_S+ABX*I_ERI_H3x2z_D2y_Px_S;
  Double I_ERI_H2x2yz_Fx2y_Px_S = I_ERI_I3x2yz_D2y_Px_S+ABX*I_ERI_H2x2yz_D2y_Px_S;
  Double I_ERI_H2xy2z_Fx2y_Px_S = I_ERI_I3xy2z_D2y_Px_S+ABX*I_ERI_H2xy2z_D2y_Px_S;
  Double I_ERI_H2x3z_Fx2y_Px_S = I_ERI_I3x3z_D2y_Px_S+ABX*I_ERI_H2x3z_D2y_Px_S;
  Double I_ERI_Hx3yz_Fx2y_Px_S = I_ERI_I2x3yz_D2y_Px_S+ABX*I_ERI_Hx3yz_D2y_Px_S;
  Double I_ERI_Hx2y2z_Fx2y_Px_S = I_ERI_I2x2y2z_D2y_Px_S+ABX*I_ERI_Hx2y2z_D2y_Px_S;
  Double I_ERI_Hxy3z_Fx2y_Px_S = I_ERI_I2xy3z_D2y_Px_S+ABX*I_ERI_Hxy3z_D2y_Px_S;
  Double I_ERI_Hx4z_Fx2y_Px_S = I_ERI_I2x4z_D2y_Px_S+ABX*I_ERI_Hx4z_D2y_Px_S;
  Double I_ERI_H4yz_Fx2y_Px_S = I_ERI_Ix4yz_D2y_Px_S+ABX*I_ERI_H4yz_D2y_Px_S;
  Double I_ERI_H3y2z_Fx2y_Px_S = I_ERI_Ix3y2z_D2y_Px_S+ABX*I_ERI_H3y2z_D2y_Px_S;
  Double I_ERI_H2y3z_Fx2y_Px_S = I_ERI_Ix2y3z_D2y_Px_S+ABX*I_ERI_H2y3z_D2y_Px_S;
  Double I_ERI_Hy4z_Fx2y_Px_S = I_ERI_Ixy4z_D2y_Px_S+ABX*I_ERI_Hy4z_D2y_Px_S;
  Double I_ERI_H5z_Fx2y_Px_S = I_ERI_Ix5z_D2y_Px_S+ABX*I_ERI_H5z_D2y_Px_S;
  Double I_ERI_H4xy_Fx2z_Px_S = I_ERI_I5xy_D2z_Px_S+ABX*I_ERI_H4xy_D2z_Px_S;
  Double I_ERI_H3x2y_Fx2z_Px_S = I_ERI_I4x2y_D2z_Px_S+ABX*I_ERI_H3x2y_D2z_Px_S;
  Double I_ERI_H3xyz_Fx2z_Px_S = I_ERI_I4xyz_D2z_Px_S+ABX*I_ERI_H3xyz_D2z_Px_S;
  Double I_ERI_H2x3y_Fx2z_Px_S = I_ERI_I3x3y_D2z_Px_S+ABX*I_ERI_H2x3y_D2z_Px_S;
  Double I_ERI_H2x2yz_Fx2z_Px_S = I_ERI_I3x2yz_D2z_Px_S+ABX*I_ERI_H2x2yz_D2z_Px_S;
  Double I_ERI_H2xy2z_Fx2z_Px_S = I_ERI_I3xy2z_D2z_Px_S+ABX*I_ERI_H2xy2z_D2z_Px_S;
  Double I_ERI_Hx4y_Fx2z_Px_S = I_ERI_I2x4y_D2z_Px_S+ABX*I_ERI_Hx4y_D2z_Px_S;
  Double I_ERI_Hx3yz_Fx2z_Px_S = I_ERI_I2x3yz_D2z_Px_S+ABX*I_ERI_Hx3yz_D2z_Px_S;
  Double I_ERI_Hx2y2z_Fx2z_Px_S = I_ERI_I2x2y2z_D2z_Px_S+ABX*I_ERI_Hx2y2z_D2z_Px_S;
  Double I_ERI_Hxy3z_Fx2z_Px_S = I_ERI_I2xy3z_D2z_Px_S+ABX*I_ERI_Hxy3z_D2z_Px_S;
  Double I_ERI_H5y_Fx2z_Px_S = I_ERI_Ix5y_D2z_Px_S+ABX*I_ERI_H5y_D2z_Px_S;
  Double I_ERI_H4yz_Fx2z_Px_S = I_ERI_Ix4yz_D2z_Px_S+ABX*I_ERI_H4yz_D2z_Px_S;
  Double I_ERI_H3y2z_Fx2z_Px_S = I_ERI_Ix3y2z_D2z_Px_S+ABX*I_ERI_H3y2z_D2z_Px_S;
  Double I_ERI_H2y3z_Fx2z_Px_S = I_ERI_Ix2y3z_D2z_Px_S+ABX*I_ERI_H2y3z_D2z_Px_S;
  Double I_ERI_Hy4z_Fx2z_Px_S = I_ERI_Ixy4z_D2z_Px_S+ABX*I_ERI_Hy4z_D2z_Px_S;
  Double I_ERI_H5x_F3y_Px_S = I_ERI_I5xy_D2y_Px_S+ABY*I_ERI_H5x_D2y_Px_S;
  Double I_ERI_H4xy_F3y_Px_S = I_ERI_I4x2y_D2y_Px_S+ABY*I_ERI_H4xy_D2y_Px_S;
  Double I_ERI_H4xz_F3y_Px_S = I_ERI_I4xyz_D2y_Px_S+ABY*I_ERI_H4xz_D2y_Px_S;
  Double I_ERI_H3x2y_F3y_Px_S = I_ERI_I3x3y_D2y_Px_S+ABY*I_ERI_H3x2y_D2y_Px_S;
  Double I_ERI_H3xyz_F3y_Px_S = I_ERI_I3x2yz_D2y_Px_S+ABY*I_ERI_H3xyz_D2y_Px_S;
  Double I_ERI_H3x2z_F3y_Px_S = I_ERI_I3xy2z_D2y_Px_S+ABY*I_ERI_H3x2z_D2y_Px_S;
  Double I_ERI_H2x3y_F3y_Px_S = I_ERI_I2x4y_D2y_Px_S+ABY*I_ERI_H2x3y_D2y_Px_S;
  Double I_ERI_H2x2yz_F3y_Px_S = I_ERI_I2x3yz_D2y_Px_S+ABY*I_ERI_H2x2yz_D2y_Px_S;
  Double I_ERI_H2xy2z_F3y_Px_S = I_ERI_I2x2y2z_D2y_Px_S+ABY*I_ERI_H2xy2z_D2y_Px_S;
  Double I_ERI_H2x3z_F3y_Px_S = I_ERI_I2xy3z_D2y_Px_S+ABY*I_ERI_H2x3z_D2y_Px_S;
  Double I_ERI_Hx4y_F3y_Px_S = I_ERI_Ix5y_D2y_Px_S+ABY*I_ERI_Hx4y_D2y_Px_S;
  Double I_ERI_Hx3yz_F3y_Px_S = I_ERI_Ix4yz_D2y_Px_S+ABY*I_ERI_Hx3yz_D2y_Px_S;
  Double I_ERI_Hx2y2z_F3y_Px_S = I_ERI_Ix3y2z_D2y_Px_S+ABY*I_ERI_Hx2y2z_D2y_Px_S;
  Double I_ERI_Hxy3z_F3y_Px_S = I_ERI_Ix2y3z_D2y_Px_S+ABY*I_ERI_Hxy3z_D2y_Px_S;
  Double I_ERI_Hx4z_F3y_Px_S = I_ERI_Ixy4z_D2y_Px_S+ABY*I_ERI_Hx4z_D2y_Px_S;
  Double I_ERI_H5y_F3y_Px_S = I_ERI_I6y_D2y_Px_S+ABY*I_ERI_H5y_D2y_Px_S;
  Double I_ERI_H4yz_F3y_Px_S = I_ERI_I5yz_D2y_Px_S+ABY*I_ERI_H4yz_D2y_Px_S;
  Double I_ERI_H3y2z_F3y_Px_S = I_ERI_I4y2z_D2y_Px_S+ABY*I_ERI_H3y2z_D2y_Px_S;
  Double I_ERI_H2y3z_F3y_Px_S = I_ERI_I3y3z_D2y_Px_S+ABY*I_ERI_H2y3z_D2y_Px_S;
  Double I_ERI_Hy4z_F3y_Px_S = I_ERI_I2y4z_D2y_Px_S+ABY*I_ERI_Hy4z_D2y_Px_S;
  Double I_ERI_H5z_F3y_Px_S = I_ERI_Iy5z_D2y_Px_S+ABY*I_ERI_H5z_D2y_Px_S;
  Double I_ERI_H4xz_F2yz_Px_S = I_ERI_I4x2z_D2y_Px_S+ABZ*I_ERI_H4xz_D2y_Px_S;
  Double I_ERI_H3xyz_F2yz_Px_S = I_ERI_I3xy2z_D2y_Px_S+ABZ*I_ERI_H3xyz_D2y_Px_S;
  Double I_ERI_H3x2z_F2yz_Px_S = I_ERI_I3x3z_D2y_Px_S+ABZ*I_ERI_H3x2z_D2y_Px_S;
  Double I_ERI_H2x2yz_F2yz_Px_S = I_ERI_I2x2y2z_D2y_Px_S+ABZ*I_ERI_H2x2yz_D2y_Px_S;
  Double I_ERI_H2xy2z_F2yz_Px_S = I_ERI_I2xy3z_D2y_Px_S+ABZ*I_ERI_H2xy2z_D2y_Px_S;
  Double I_ERI_H2x3z_F2yz_Px_S = I_ERI_I2x4z_D2y_Px_S+ABZ*I_ERI_H2x3z_D2y_Px_S;
  Double I_ERI_Hx3yz_F2yz_Px_S = I_ERI_Ix3y2z_D2y_Px_S+ABZ*I_ERI_Hx3yz_D2y_Px_S;
  Double I_ERI_Hx2y2z_F2yz_Px_S = I_ERI_Ix2y3z_D2y_Px_S+ABZ*I_ERI_Hx2y2z_D2y_Px_S;
  Double I_ERI_Hxy3z_F2yz_Px_S = I_ERI_Ixy4z_D2y_Px_S+ABZ*I_ERI_Hxy3z_D2y_Px_S;
  Double I_ERI_Hx4z_F2yz_Px_S = I_ERI_Ix5z_D2y_Px_S+ABZ*I_ERI_Hx4z_D2y_Px_S;
  Double I_ERI_H4yz_F2yz_Px_S = I_ERI_I4y2z_D2y_Px_S+ABZ*I_ERI_H4yz_D2y_Px_S;
  Double I_ERI_H3y2z_F2yz_Px_S = I_ERI_I3y3z_D2y_Px_S+ABZ*I_ERI_H3y2z_D2y_Px_S;
  Double I_ERI_H2y3z_F2yz_Px_S = I_ERI_I2y4z_D2y_Px_S+ABZ*I_ERI_H2y3z_D2y_Px_S;
  Double I_ERI_Hy4z_F2yz_Px_S = I_ERI_Iy5z_D2y_Px_S+ABZ*I_ERI_Hy4z_D2y_Px_S;
  Double I_ERI_H5z_F2yz_Px_S = I_ERI_I6z_D2y_Px_S+ABZ*I_ERI_H5z_D2y_Px_S;
  Double I_ERI_H5x_F3z_Px_S = I_ERI_I5xz_D2z_Px_S+ABZ*I_ERI_H5x_D2z_Px_S;
  Double I_ERI_H4xy_F3z_Px_S = I_ERI_I4xyz_D2z_Px_S+ABZ*I_ERI_H4xy_D2z_Px_S;
  Double I_ERI_H4xz_F3z_Px_S = I_ERI_I4x2z_D2z_Px_S+ABZ*I_ERI_H4xz_D2z_Px_S;
  Double I_ERI_H3x2y_F3z_Px_S = I_ERI_I3x2yz_D2z_Px_S+ABZ*I_ERI_H3x2y_D2z_Px_S;
  Double I_ERI_H3xyz_F3z_Px_S = I_ERI_I3xy2z_D2z_Px_S+ABZ*I_ERI_H3xyz_D2z_Px_S;
  Double I_ERI_H3x2z_F3z_Px_S = I_ERI_I3x3z_D2z_Px_S+ABZ*I_ERI_H3x2z_D2z_Px_S;
  Double I_ERI_H2x3y_F3z_Px_S = I_ERI_I2x3yz_D2z_Px_S+ABZ*I_ERI_H2x3y_D2z_Px_S;
  Double I_ERI_H2x2yz_F3z_Px_S = I_ERI_I2x2y2z_D2z_Px_S+ABZ*I_ERI_H2x2yz_D2z_Px_S;
  Double I_ERI_H2xy2z_F3z_Px_S = I_ERI_I2xy3z_D2z_Px_S+ABZ*I_ERI_H2xy2z_D2z_Px_S;
  Double I_ERI_H2x3z_F3z_Px_S = I_ERI_I2x4z_D2z_Px_S+ABZ*I_ERI_H2x3z_D2z_Px_S;
  Double I_ERI_Hx4y_F3z_Px_S = I_ERI_Ix4yz_D2z_Px_S+ABZ*I_ERI_Hx4y_D2z_Px_S;
  Double I_ERI_Hx3yz_F3z_Px_S = I_ERI_Ix3y2z_D2z_Px_S+ABZ*I_ERI_Hx3yz_D2z_Px_S;
  Double I_ERI_Hx2y2z_F3z_Px_S = I_ERI_Ix2y3z_D2z_Px_S+ABZ*I_ERI_Hx2y2z_D2z_Px_S;
  Double I_ERI_Hxy3z_F3z_Px_S = I_ERI_Ixy4z_D2z_Px_S+ABZ*I_ERI_Hxy3z_D2z_Px_S;
  Double I_ERI_Hx4z_F3z_Px_S = I_ERI_Ix5z_D2z_Px_S+ABZ*I_ERI_Hx4z_D2z_Px_S;
  Double I_ERI_H5y_F3z_Px_S = I_ERI_I5yz_D2z_Px_S+ABZ*I_ERI_H5y_D2z_Px_S;
  Double I_ERI_H4yz_F3z_Px_S = I_ERI_I4y2z_D2z_Px_S+ABZ*I_ERI_H4yz_D2z_Px_S;
  Double I_ERI_H3y2z_F3z_Px_S = I_ERI_I3y3z_D2z_Px_S+ABZ*I_ERI_H3y2z_D2z_Px_S;
  Double I_ERI_H2y3z_F3z_Px_S = I_ERI_I2y4z_D2z_Px_S+ABZ*I_ERI_H2y3z_D2z_Px_S;
  Double I_ERI_Hy4z_F3z_Px_S = I_ERI_Iy5z_D2z_Px_S+ABZ*I_ERI_Hy4z_D2z_Px_S;
  Double I_ERI_H5z_F3z_Px_S = I_ERI_I6z_D2z_Px_S+ABZ*I_ERI_H5z_D2z_Px_S;
  Double I_ERI_H5x_F3x_Py_S = I_ERI_I6x_D2x_Py_S+ABX*I_ERI_H5x_D2x_Py_S;
  Double I_ERI_H4xy_F3x_Py_S = I_ERI_I5xy_D2x_Py_S+ABX*I_ERI_H4xy_D2x_Py_S;
  Double I_ERI_H4xz_F3x_Py_S = I_ERI_I5xz_D2x_Py_S+ABX*I_ERI_H4xz_D2x_Py_S;
  Double I_ERI_H3x2y_F3x_Py_S = I_ERI_I4x2y_D2x_Py_S+ABX*I_ERI_H3x2y_D2x_Py_S;
  Double I_ERI_H3xyz_F3x_Py_S = I_ERI_I4xyz_D2x_Py_S+ABX*I_ERI_H3xyz_D2x_Py_S;
  Double I_ERI_H3x2z_F3x_Py_S = I_ERI_I4x2z_D2x_Py_S+ABX*I_ERI_H3x2z_D2x_Py_S;
  Double I_ERI_H2x3y_F3x_Py_S = I_ERI_I3x3y_D2x_Py_S+ABX*I_ERI_H2x3y_D2x_Py_S;
  Double I_ERI_H2x2yz_F3x_Py_S = I_ERI_I3x2yz_D2x_Py_S+ABX*I_ERI_H2x2yz_D2x_Py_S;
  Double I_ERI_H2xy2z_F3x_Py_S = I_ERI_I3xy2z_D2x_Py_S+ABX*I_ERI_H2xy2z_D2x_Py_S;
  Double I_ERI_H2x3z_F3x_Py_S = I_ERI_I3x3z_D2x_Py_S+ABX*I_ERI_H2x3z_D2x_Py_S;
  Double I_ERI_Hx4y_F3x_Py_S = I_ERI_I2x4y_D2x_Py_S+ABX*I_ERI_Hx4y_D2x_Py_S;
  Double I_ERI_Hx3yz_F3x_Py_S = I_ERI_I2x3yz_D2x_Py_S+ABX*I_ERI_Hx3yz_D2x_Py_S;
  Double I_ERI_Hx2y2z_F3x_Py_S = I_ERI_I2x2y2z_D2x_Py_S+ABX*I_ERI_Hx2y2z_D2x_Py_S;
  Double I_ERI_Hxy3z_F3x_Py_S = I_ERI_I2xy3z_D2x_Py_S+ABX*I_ERI_Hxy3z_D2x_Py_S;
  Double I_ERI_Hx4z_F3x_Py_S = I_ERI_I2x4z_D2x_Py_S+ABX*I_ERI_Hx4z_D2x_Py_S;
  Double I_ERI_H5y_F3x_Py_S = I_ERI_Ix5y_D2x_Py_S+ABX*I_ERI_H5y_D2x_Py_S;
  Double I_ERI_H4yz_F3x_Py_S = I_ERI_Ix4yz_D2x_Py_S+ABX*I_ERI_H4yz_D2x_Py_S;
  Double I_ERI_H3y2z_F3x_Py_S = I_ERI_Ix3y2z_D2x_Py_S+ABX*I_ERI_H3y2z_D2x_Py_S;
  Double I_ERI_H2y3z_F3x_Py_S = I_ERI_Ix2y3z_D2x_Py_S+ABX*I_ERI_H2y3z_D2x_Py_S;
  Double I_ERI_Hy4z_F3x_Py_S = I_ERI_Ixy4z_D2x_Py_S+ABX*I_ERI_Hy4z_D2x_Py_S;
  Double I_ERI_H5z_F3x_Py_S = I_ERI_Ix5z_D2x_Py_S+ABX*I_ERI_H5z_D2x_Py_S;
  Double I_ERI_H4xy_F2xy_Py_S = I_ERI_I4x2y_D2x_Py_S+ABY*I_ERI_H4xy_D2x_Py_S;
  Double I_ERI_H4xz_F2xy_Py_S = I_ERI_I4xyz_D2x_Py_S+ABY*I_ERI_H4xz_D2x_Py_S;
  Double I_ERI_H3x2y_F2xy_Py_S = I_ERI_I3x3y_D2x_Py_S+ABY*I_ERI_H3x2y_D2x_Py_S;
  Double I_ERI_H3xyz_F2xy_Py_S = I_ERI_I3x2yz_D2x_Py_S+ABY*I_ERI_H3xyz_D2x_Py_S;
  Double I_ERI_H3x2z_F2xy_Py_S = I_ERI_I3xy2z_D2x_Py_S+ABY*I_ERI_H3x2z_D2x_Py_S;
  Double I_ERI_H2x3y_F2xy_Py_S = I_ERI_I2x4y_D2x_Py_S+ABY*I_ERI_H2x3y_D2x_Py_S;
  Double I_ERI_H2x2yz_F2xy_Py_S = I_ERI_I2x3yz_D2x_Py_S+ABY*I_ERI_H2x2yz_D2x_Py_S;
  Double I_ERI_H2xy2z_F2xy_Py_S = I_ERI_I2x2y2z_D2x_Py_S+ABY*I_ERI_H2xy2z_D2x_Py_S;
  Double I_ERI_H2x3z_F2xy_Py_S = I_ERI_I2xy3z_D2x_Py_S+ABY*I_ERI_H2x3z_D2x_Py_S;
  Double I_ERI_Hx4y_F2xy_Py_S = I_ERI_Ix5y_D2x_Py_S+ABY*I_ERI_Hx4y_D2x_Py_S;
  Double I_ERI_Hx3yz_F2xy_Py_S = I_ERI_Ix4yz_D2x_Py_S+ABY*I_ERI_Hx3yz_D2x_Py_S;
  Double I_ERI_Hx2y2z_F2xy_Py_S = I_ERI_Ix3y2z_D2x_Py_S+ABY*I_ERI_Hx2y2z_D2x_Py_S;
  Double I_ERI_Hxy3z_F2xy_Py_S = I_ERI_Ix2y3z_D2x_Py_S+ABY*I_ERI_Hxy3z_D2x_Py_S;
  Double I_ERI_Hx4z_F2xy_Py_S = I_ERI_Ixy4z_D2x_Py_S+ABY*I_ERI_Hx4z_D2x_Py_S;
  Double I_ERI_H5y_F2xy_Py_S = I_ERI_I6y_D2x_Py_S+ABY*I_ERI_H5y_D2x_Py_S;
  Double I_ERI_H4yz_F2xy_Py_S = I_ERI_I5yz_D2x_Py_S+ABY*I_ERI_H4yz_D2x_Py_S;
  Double I_ERI_H3y2z_F2xy_Py_S = I_ERI_I4y2z_D2x_Py_S+ABY*I_ERI_H3y2z_D2x_Py_S;
  Double I_ERI_H2y3z_F2xy_Py_S = I_ERI_I3y3z_D2x_Py_S+ABY*I_ERI_H2y3z_D2x_Py_S;
  Double I_ERI_Hy4z_F2xy_Py_S = I_ERI_I2y4z_D2x_Py_S+ABY*I_ERI_Hy4z_D2x_Py_S;
  Double I_ERI_H5z_F2xy_Py_S = I_ERI_Iy5z_D2x_Py_S+ABY*I_ERI_H5z_D2x_Py_S;
  Double I_ERI_H4xz_F2xz_Py_S = I_ERI_I4x2z_D2x_Py_S+ABZ*I_ERI_H4xz_D2x_Py_S;
  Double I_ERI_H3xyz_F2xz_Py_S = I_ERI_I3xy2z_D2x_Py_S+ABZ*I_ERI_H3xyz_D2x_Py_S;
  Double I_ERI_H3x2z_F2xz_Py_S = I_ERI_I3x3z_D2x_Py_S+ABZ*I_ERI_H3x2z_D2x_Py_S;
  Double I_ERI_H2x2yz_F2xz_Py_S = I_ERI_I2x2y2z_D2x_Py_S+ABZ*I_ERI_H2x2yz_D2x_Py_S;
  Double I_ERI_H2xy2z_F2xz_Py_S = I_ERI_I2xy3z_D2x_Py_S+ABZ*I_ERI_H2xy2z_D2x_Py_S;
  Double I_ERI_H2x3z_F2xz_Py_S = I_ERI_I2x4z_D2x_Py_S+ABZ*I_ERI_H2x3z_D2x_Py_S;
  Double I_ERI_Hx3yz_F2xz_Py_S = I_ERI_Ix3y2z_D2x_Py_S+ABZ*I_ERI_Hx3yz_D2x_Py_S;
  Double I_ERI_Hx2y2z_F2xz_Py_S = I_ERI_Ix2y3z_D2x_Py_S+ABZ*I_ERI_Hx2y2z_D2x_Py_S;
  Double I_ERI_Hxy3z_F2xz_Py_S = I_ERI_Ixy4z_D2x_Py_S+ABZ*I_ERI_Hxy3z_D2x_Py_S;
  Double I_ERI_Hx4z_F2xz_Py_S = I_ERI_Ix5z_D2x_Py_S+ABZ*I_ERI_Hx4z_D2x_Py_S;
  Double I_ERI_H4yz_F2xz_Py_S = I_ERI_I4y2z_D2x_Py_S+ABZ*I_ERI_H4yz_D2x_Py_S;
  Double I_ERI_H3y2z_F2xz_Py_S = I_ERI_I3y3z_D2x_Py_S+ABZ*I_ERI_H3y2z_D2x_Py_S;
  Double I_ERI_H2y3z_F2xz_Py_S = I_ERI_I2y4z_D2x_Py_S+ABZ*I_ERI_H2y3z_D2x_Py_S;
  Double I_ERI_Hy4z_F2xz_Py_S = I_ERI_Iy5z_D2x_Py_S+ABZ*I_ERI_Hy4z_D2x_Py_S;
  Double I_ERI_H5z_F2xz_Py_S = I_ERI_I6z_D2x_Py_S+ABZ*I_ERI_H5z_D2x_Py_S;
  Double I_ERI_H4xz_Fx2y_Py_S = I_ERI_I5xz_D2y_Py_S+ABX*I_ERI_H4xz_D2y_Py_S;
  Double I_ERI_H3xyz_Fx2y_Py_S = I_ERI_I4xyz_D2y_Py_S+ABX*I_ERI_H3xyz_D2y_Py_S;
  Double I_ERI_H3x2z_Fx2y_Py_S = I_ERI_I4x2z_D2y_Py_S+ABX*I_ERI_H3x2z_D2y_Py_S;
  Double I_ERI_H2x2yz_Fx2y_Py_S = I_ERI_I3x2yz_D2y_Py_S+ABX*I_ERI_H2x2yz_D2y_Py_S;
  Double I_ERI_H2xy2z_Fx2y_Py_S = I_ERI_I3xy2z_D2y_Py_S+ABX*I_ERI_H2xy2z_D2y_Py_S;
  Double I_ERI_H2x3z_Fx2y_Py_S = I_ERI_I3x3z_D2y_Py_S+ABX*I_ERI_H2x3z_D2y_Py_S;
  Double I_ERI_Hx3yz_Fx2y_Py_S = I_ERI_I2x3yz_D2y_Py_S+ABX*I_ERI_Hx3yz_D2y_Py_S;
  Double I_ERI_Hx2y2z_Fx2y_Py_S = I_ERI_I2x2y2z_D2y_Py_S+ABX*I_ERI_Hx2y2z_D2y_Py_S;
  Double I_ERI_Hxy3z_Fx2y_Py_S = I_ERI_I2xy3z_D2y_Py_S+ABX*I_ERI_Hxy3z_D2y_Py_S;
  Double I_ERI_Hx4z_Fx2y_Py_S = I_ERI_I2x4z_D2y_Py_S+ABX*I_ERI_Hx4z_D2y_Py_S;
  Double I_ERI_H4yz_Fx2y_Py_S = I_ERI_Ix4yz_D2y_Py_S+ABX*I_ERI_H4yz_D2y_Py_S;
  Double I_ERI_H3y2z_Fx2y_Py_S = I_ERI_Ix3y2z_D2y_Py_S+ABX*I_ERI_H3y2z_D2y_Py_S;
  Double I_ERI_H2y3z_Fx2y_Py_S = I_ERI_Ix2y3z_D2y_Py_S+ABX*I_ERI_H2y3z_D2y_Py_S;
  Double I_ERI_Hy4z_Fx2y_Py_S = I_ERI_Ixy4z_D2y_Py_S+ABX*I_ERI_Hy4z_D2y_Py_S;
  Double I_ERI_H5z_Fx2y_Py_S = I_ERI_Ix5z_D2y_Py_S+ABX*I_ERI_H5z_D2y_Py_S;
  Double I_ERI_H4xy_Fx2z_Py_S = I_ERI_I5xy_D2z_Py_S+ABX*I_ERI_H4xy_D2z_Py_S;
  Double I_ERI_H3x2y_Fx2z_Py_S = I_ERI_I4x2y_D2z_Py_S+ABX*I_ERI_H3x2y_D2z_Py_S;
  Double I_ERI_H3xyz_Fx2z_Py_S = I_ERI_I4xyz_D2z_Py_S+ABX*I_ERI_H3xyz_D2z_Py_S;
  Double I_ERI_H2x3y_Fx2z_Py_S = I_ERI_I3x3y_D2z_Py_S+ABX*I_ERI_H2x3y_D2z_Py_S;
  Double I_ERI_H2x2yz_Fx2z_Py_S = I_ERI_I3x2yz_D2z_Py_S+ABX*I_ERI_H2x2yz_D2z_Py_S;
  Double I_ERI_H2xy2z_Fx2z_Py_S = I_ERI_I3xy2z_D2z_Py_S+ABX*I_ERI_H2xy2z_D2z_Py_S;
  Double I_ERI_Hx4y_Fx2z_Py_S = I_ERI_I2x4y_D2z_Py_S+ABX*I_ERI_Hx4y_D2z_Py_S;
  Double I_ERI_Hx3yz_Fx2z_Py_S = I_ERI_I2x3yz_D2z_Py_S+ABX*I_ERI_Hx3yz_D2z_Py_S;
  Double I_ERI_Hx2y2z_Fx2z_Py_S = I_ERI_I2x2y2z_D2z_Py_S+ABX*I_ERI_Hx2y2z_D2z_Py_S;
  Double I_ERI_Hxy3z_Fx2z_Py_S = I_ERI_I2xy3z_D2z_Py_S+ABX*I_ERI_Hxy3z_D2z_Py_S;
  Double I_ERI_H5y_Fx2z_Py_S = I_ERI_Ix5y_D2z_Py_S+ABX*I_ERI_H5y_D2z_Py_S;
  Double I_ERI_H4yz_Fx2z_Py_S = I_ERI_Ix4yz_D2z_Py_S+ABX*I_ERI_H4yz_D2z_Py_S;
  Double I_ERI_H3y2z_Fx2z_Py_S = I_ERI_Ix3y2z_D2z_Py_S+ABX*I_ERI_H3y2z_D2z_Py_S;
  Double I_ERI_H2y3z_Fx2z_Py_S = I_ERI_Ix2y3z_D2z_Py_S+ABX*I_ERI_H2y3z_D2z_Py_S;
  Double I_ERI_Hy4z_Fx2z_Py_S = I_ERI_Ixy4z_D2z_Py_S+ABX*I_ERI_Hy4z_D2z_Py_S;
  Double I_ERI_H5x_F3y_Py_S = I_ERI_I5xy_D2y_Py_S+ABY*I_ERI_H5x_D2y_Py_S;
  Double I_ERI_H4xy_F3y_Py_S = I_ERI_I4x2y_D2y_Py_S+ABY*I_ERI_H4xy_D2y_Py_S;
  Double I_ERI_H4xz_F3y_Py_S = I_ERI_I4xyz_D2y_Py_S+ABY*I_ERI_H4xz_D2y_Py_S;
  Double I_ERI_H3x2y_F3y_Py_S = I_ERI_I3x3y_D2y_Py_S+ABY*I_ERI_H3x2y_D2y_Py_S;
  Double I_ERI_H3xyz_F3y_Py_S = I_ERI_I3x2yz_D2y_Py_S+ABY*I_ERI_H3xyz_D2y_Py_S;
  Double I_ERI_H3x2z_F3y_Py_S = I_ERI_I3xy2z_D2y_Py_S+ABY*I_ERI_H3x2z_D2y_Py_S;
  Double I_ERI_H2x3y_F3y_Py_S = I_ERI_I2x4y_D2y_Py_S+ABY*I_ERI_H2x3y_D2y_Py_S;
  Double I_ERI_H2x2yz_F3y_Py_S = I_ERI_I2x3yz_D2y_Py_S+ABY*I_ERI_H2x2yz_D2y_Py_S;
  Double I_ERI_H2xy2z_F3y_Py_S = I_ERI_I2x2y2z_D2y_Py_S+ABY*I_ERI_H2xy2z_D2y_Py_S;
  Double I_ERI_H2x3z_F3y_Py_S = I_ERI_I2xy3z_D2y_Py_S+ABY*I_ERI_H2x3z_D2y_Py_S;
  Double I_ERI_Hx4y_F3y_Py_S = I_ERI_Ix5y_D2y_Py_S+ABY*I_ERI_Hx4y_D2y_Py_S;
  Double I_ERI_Hx3yz_F3y_Py_S = I_ERI_Ix4yz_D2y_Py_S+ABY*I_ERI_Hx3yz_D2y_Py_S;
  Double I_ERI_Hx2y2z_F3y_Py_S = I_ERI_Ix3y2z_D2y_Py_S+ABY*I_ERI_Hx2y2z_D2y_Py_S;
  Double I_ERI_Hxy3z_F3y_Py_S = I_ERI_Ix2y3z_D2y_Py_S+ABY*I_ERI_Hxy3z_D2y_Py_S;
  Double I_ERI_Hx4z_F3y_Py_S = I_ERI_Ixy4z_D2y_Py_S+ABY*I_ERI_Hx4z_D2y_Py_S;
  Double I_ERI_H5y_F3y_Py_S = I_ERI_I6y_D2y_Py_S+ABY*I_ERI_H5y_D2y_Py_S;
  Double I_ERI_H4yz_F3y_Py_S = I_ERI_I5yz_D2y_Py_S+ABY*I_ERI_H4yz_D2y_Py_S;
  Double I_ERI_H3y2z_F3y_Py_S = I_ERI_I4y2z_D2y_Py_S+ABY*I_ERI_H3y2z_D2y_Py_S;
  Double I_ERI_H2y3z_F3y_Py_S = I_ERI_I3y3z_D2y_Py_S+ABY*I_ERI_H2y3z_D2y_Py_S;
  Double I_ERI_Hy4z_F3y_Py_S = I_ERI_I2y4z_D2y_Py_S+ABY*I_ERI_Hy4z_D2y_Py_S;
  Double I_ERI_H5z_F3y_Py_S = I_ERI_Iy5z_D2y_Py_S+ABY*I_ERI_H5z_D2y_Py_S;
  Double I_ERI_H4xz_F2yz_Py_S = I_ERI_I4x2z_D2y_Py_S+ABZ*I_ERI_H4xz_D2y_Py_S;
  Double I_ERI_H3xyz_F2yz_Py_S = I_ERI_I3xy2z_D2y_Py_S+ABZ*I_ERI_H3xyz_D2y_Py_S;
  Double I_ERI_H3x2z_F2yz_Py_S = I_ERI_I3x3z_D2y_Py_S+ABZ*I_ERI_H3x2z_D2y_Py_S;
  Double I_ERI_H2x2yz_F2yz_Py_S = I_ERI_I2x2y2z_D2y_Py_S+ABZ*I_ERI_H2x2yz_D2y_Py_S;
  Double I_ERI_H2xy2z_F2yz_Py_S = I_ERI_I2xy3z_D2y_Py_S+ABZ*I_ERI_H2xy2z_D2y_Py_S;
  Double I_ERI_H2x3z_F2yz_Py_S = I_ERI_I2x4z_D2y_Py_S+ABZ*I_ERI_H2x3z_D2y_Py_S;
  Double I_ERI_Hx3yz_F2yz_Py_S = I_ERI_Ix3y2z_D2y_Py_S+ABZ*I_ERI_Hx3yz_D2y_Py_S;
  Double I_ERI_Hx2y2z_F2yz_Py_S = I_ERI_Ix2y3z_D2y_Py_S+ABZ*I_ERI_Hx2y2z_D2y_Py_S;
  Double I_ERI_Hxy3z_F2yz_Py_S = I_ERI_Ixy4z_D2y_Py_S+ABZ*I_ERI_Hxy3z_D2y_Py_S;
  Double I_ERI_Hx4z_F2yz_Py_S = I_ERI_Ix5z_D2y_Py_S+ABZ*I_ERI_Hx4z_D2y_Py_S;
  Double I_ERI_H4yz_F2yz_Py_S = I_ERI_I4y2z_D2y_Py_S+ABZ*I_ERI_H4yz_D2y_Py_S;
  Double I_ERI_H3y2z_F2yz_Py_S = I_ERI_I3y3z_D2y_Py_S+ABZ*I_ERI_H3y2z_D2y_Py_S;
  Double I_ERI_H2y3z_F2yz_Py_S = I_ERI_I2y4z_D2y_Py_S+ABZ*I_ERI_H2y3z_D2y_Py_S;
  Double I_ERI_Hy4z_F2yz_Py_S = I_ERI_Iy5z_D2y_Py_S+ABZ*I_ERI_Hy4z_D2y_Py_S;
  Double I_ERI_H5z_F2yz_Py_S = I_ERI_I6z_D2y_Py_S+ABZ*I_ERI_H5z_D2y_Py_S;
  Double I_ERI_H5x_F3z_Py_S = I_ERI_I5xz_D2z_Py_S+ABZ*I_ERI_H5x_D2z_Py_S;
  Double I_ERI_H4xy_F3z_Py_S = I_ERI_I4xyz_D2z_Py_S+ABZ*I_ERI_H4xy_D2z_Py_S;
  Double I_ERI_H4xz_F3z_Py_S = I_ERI_I4x2z_D2z_Py_S+ABZ*I_ERI_H4xz_D2z_Py_S;
  Double I_ERI_H3x2y_F3z_Py_S = I_ERI_I3x2yz_D2z_Py_S+ABZ*I_ERI_H3x2y_D2z_Py_S;
  Double I_ERI_H3xyz_F3z_Py_S = I_ERI_I3xy2z_D2z_Py_S+ABZ*I_ERI_H3xyz_D2z_Py_S;
  Double I_ERI_H3x2z_F3z_Py_S = I_ERI_I3x3z_D2z_Py_S+ABZ*I_ERI_H3x2z_D2z_Py_S;
  Double I_ERI_H2x3y_F3z_Py_S = I_ERI_I2x3yz_D2z_Py_S+ABZ*I_ERI_H2x3y_D2z_Py_S;
  Double I_ERI_H2x2yz_F3z_Py_S = I_ERI_I2x2y2z_D2z_Py_S+ABZ*I_ERI_H2x2yz_D2z_Py_S;
  Double I_ERI_H2xy2z_F3z_Py_S = I_ERI_I2xy3z_D2z_Py_S+ABZ*I_ERI_H2xy2z_D2z_Py_S;
  Double I_ERI_H2x3z_F3z_Py_S = I_ERI_I2x4z_D2z_Py_S+ABZ*I_ERI_H2x3z_D2z_Py_S;
  Double I_ERI_Hx4y_F3z_Py_S = I_ERI_Ix4yz_D2z_Py_S+ABZ*I_ERI_Hx4y_D2z_Py_S;
  Double I_ERI_Hx3yz_F3z_Py_S = I_ERI_Ix3y2z_D2z_Py_S+ABZ*I_ERI_Hx3yz_D2z_Py_S;
  Double I_ERI_Hx2y2z_F3z_Py_S = I_ERI_Ix2y3z_D2z_Py_S+ABZ*I_ERI_Hx2y2z_D2z_Py_S;
  Double I_ERI_Hxy3z_F3z_Py_S = I_ERI_Ixy4z_D2z_Py_S+ABZ*I_ERI_Hxy3z_D2z_Py_S;
  Double I_ERI_Hx4z_F3z_Py_S = I_ERI_Ix5z_D2z_Py_S+ABZ*I_ERI_Hx4z_D2z_Py_S;
  Double I_ERI_H5y_F3z_Py_S = I_ERI_I5yz_D2z_Py_S+ABZ*I_ERI_H5y_D2z_Py_S;
  Double I_ERI_H4yz_F3z_Py_S = I_ERI_I4y2z_D2z_Py_S+ABZ*I_ERI_H4yz_D2z_Py_S;
  Double I_ERI_H3y2z_F3z_Py_S = I_ERI_I3y3z_D2z_Py_S+ABZ*I_ERI_H3y2z_D2z_Py_S;
  Double I_ERI_H2y3z_F3z_Py_S = I_ERI_I2y4z_D2z_Py_S+ABZ*I_ERI_H2y3z_D2z_Py_S;
  Double I_ERI_Hy4z_F3z_Py_S = I_ERI_Iy5z_D2z_Py_S+ABZ*I_ERI_Hy4z_D2z_Py_S;
  Double I_ERI_H5z_F3z_Py_S = I_ERI_I6z_D2z_Py_S+ABZ*I_ERI_H5z_D2z_Py_S;
  Double I_ERI_H5x_F3x_Pz_S = I_ERI_I6x_D2x_Pz_S+ABX*I_ERI_H5x_D2x_Pz_S;
  Double I_ERI_H4xy_F3x_Pz_S = I_ERI_I5xy_D2x_Pz_S+ABX*I_ERI_H4xy_D2x_Pz_S;
  Double I_ERI_H4xz_F3x_Pz_S = I_ERI_I5xz_D2x_Pz_S+ABX*I_ERI_H4xz_D2x_Pz_S;
  Double I_ERI_H3x2y_F3x_Pz_S = I_ERI_I4x2y_D2x_Pz_S+ABX*I_ERI_H3x2y_D2x_Pz_S;
  Double I_ERI_H3xyz_F3x_Pz_S = I_ERI_I4xyz_D2x_Pz_S+ABX*I_ERI_H3xyz_D2x_Pz_S;
  Double I_ERI_H3x2z_F3x_Pz_S = I_ERI_I4x2z_D2x_Pz_S+ABX*I_ERI_H3x2z_D2x_Pz_S;
  Double I_ERI_H2x3y_F3x_Pz_S = I_ERI_I3x3y_D2x_Pz_S+ABX*I_ERI_H2x3y_D2x_Pz_S;
  Double I_ERI_H2x2yz_F3x_Pz_S = I_ERI_I3x2yz_D2x_Pz_S+ABX*I_ERI_H2x2yz_D2x_Pz_S;
  Double I_ERI_H2xy2z_F3x_Pz_S = I_ERI_I3xy2z_D2x_Pz_S+ABX*I_ERI_H2xy2z_D2x_Pz_S;
  Double I_ERI_H2x3z_F3x_Pz_S = I_ERI_I3x3z_D2x_Pz_S+ABX*I_ERI_H2x3z_D2x_Pz_S;
  Double I_ERI_Hx4y_F3x_Pz_S = I_ERI_I2x4y_D2x_Pz_S+ABX*I_ERI_Hx4y_D2x_Pz_S;
  Double I_ERI_Hx3yz_F3x_Pz_S = I_ERI_I2x3yz_D2x_Pz_S+ABX*I_ERI_Hx3yz_D2x_Pz_S;
  Double I_ERI_Hx2y2z_F3x_Pz_S = I_ERI_I2x2y2z_D2x_Pz_S+ABX*I_ERI_Hx2y2z_D2x_Pz_S;
  Double I_ERI_Hxy3z_F3x_Pz_S = I_ERI_I2xy3z_D2x_Pz_S+ABX*I_ERI_Hxy3z_D2x_Pz_S;
  Double I_ERI_Hx4z_F3x_Pz_S = I_ERI_I2x4z_D2x_Pz_S+ABX*I_ERI_Hx4z_D2x_Pz_S;
  Double I_ERI_H5y_F3x_Pz_S = I_ERI_Ix5y_D2x_Pz_S+ABX*I_ERI_H5y_D2x_Pz_S;
  Double I_ERI_H4yz_F3x_Pz_S = I_ERI_Ix4yz_D2x_Pz_S+ABX*I_ERI_H4yz_D2x_Pz_S;
  Double I_ERI_H3y2z_F3x_Pz_S = I_ERI_Ix3y2z_D2x_Pz_S+ABX*I_ERI_H3y2z_D2x_Pz_S;
  Double I_ERI_H2y3z_F3x_Pz_S = I_ERI_Ix2y3z_D2x_Pz_S+ABX*I_ERI_H2y3z_D2x_Pz_S;
  Double I_ERI_Hy4z_F3x_Pz_S = I_ERI_Ixy4z_D2x_Pz_S+ABX*I_ERI_Hy4z_D2x_Pz_S;
  Double I_ERI_H5z_F3x_Pz_S = I_ERI_Ix5z_D2x_Pz_S+ABX*I_ERI_H5z_D2x_Pz_S;
  Double I_ERI_H4xy_F2xy_Pz_S = I_ERI_I4x2y_D2x_Pz_S+ABY*I_ERI_H4xy_D2x_Pz_S;
  Double I_ERI_H4xz_F2xy_Pz_S = I_ERI_I4xyz_D2x_Pz_S+ABY*I_ERI_H4xz_D2x_Pz_S;
  Double I_ERI_H3x2y_F2xy_Pz_S = I_ERI_I3x3y_D2x_Pz_S+ABY*I_ERI_H3x2y_D2x_Pz_S;
  Double I_ERI_H3xyz_F2xy_Pz_S = I_ERI_I3x2yz_D2x_Pz_S+ABY*I_ERI_H3xyz_D2x_Pz_S;
  Double I_ERI_H3x2z_F2xy_Pz_S = I_ERI_I3xy2z_D2x_Pz_S+ABY*I_ERI_H3x2z_D2x_Pz_S;
  Double I_ERI_H2x3y_F2xy_Pz_S = I_ERI_I2x4y_D2x_Pz_S+ABY*I_ERI_H2x3y_D2x_Pz_S;
  Double I_ERI_H2x2yz_F2xy_Pz_S = I_ERI_I2x3yz_D2x_Pz_S+ABY*I_ERI_H2x2yz_D2x_Pz_S;
  Double I_ERI_H2xy2z_F2xy_Pz_S = I_ERI_I2x2y2z_D2x_Pz_S+ABY*I_ERI_H2xy2z_D2x_Pz_S;
  Double I_ERI_H2x3z_F2xy_Pz_S = I_ERI_I2xy3z_D2x_Pz_S+ABY*I_ERI_H2x3z_D2x_Pz_S;
  Double I_ERI_Hx4y_F2xy_Pz_S = I_ERI_Ix5y_D2x_Pz_S+ABY*I_ERI_Hx4y_D2x_Pz_S;
  Double I_ERI_Hx3yz_F2xy_Pz_S = I_ERI_Ix4yz_D2x_Pz_S+ABY*I_ERI_Hx3yz_D2x_Pz_S;
  Double I_ERI_Hx2y2z_F2xy_Pz_S = I_ERI_Ix3y2z_D2x_Pz_S+ABY*I_ERI_Hx2y2z_D2x_Pz_S;
  Double I_ERI_Hxy3z_F2xy_Pz_S = I_ERI_Ix2y3z_D2x_Pz_S+ABY*I_ERI_Hxy3z_D2x_Pz_S;
  Double I_ERI_Hx4z_F2xy_Pz_S = I_ERI_Ixy4z_D2x_Pz_S+ABY*I_ERI_Hx4z_D2x_Pz_S;
  Double I_ERI_H5y_F2xy_Pz_S = I_ERI_I6y_D2x_Pz_S+ABY*I_ERI_H5y_D2x_Pz_S;
  Double I_ERI_H4yz_F2xy_Pz_S = I_ERI_I5yz_D2x_Pz_S+ABY*I_ERI_H4yz_D2x_Pz_S;
  Double I_ERI_H3y2z_F2xy_Pz_S = I_ERI_I4y2z_D2x_Pz_S+ABY*I_ERI_H3y2z_D2x_Pz_S;
  Double I_ERI_H2y3z_F2xy_Pz_S = I_ERI_I3y3z_D2x_Pz_S+ABY*I_ERI_H2y3z_D2x_Pz_S;
  Double I_ERI_Hy4z_F2xy_Pz_S = I_ERI_I2y4z_D2x_Pz_S+ABY*I_ERI_Hy4z_D2x_Pz_S;
  Double I_ERI_H5z_F2xy_Pz_S = I_ERI_Iy5z_D2x_Pz_S+ABY*I_ERI_H5z_D2x_Pz_S;
  Double I_ERI_H4xz_F2xz_Pz_S = I_ERI_I4x2z_D2x_Pz_S+ABZ*I_ERI_H4xz_D2x_Pz_S;
  Double I_ERI_H3xyz_F2xz_Pz_S = I_ERI_I3xy2z_D2x_Pz_S+ABZ*I_ERI_H3xyz_D2x_Pz_S;
  Double I_ERI_H3x2z_F2xz_Pz_S = I_ERI_I3x3z_D2x_Pz_S+ABZ*I_ERI_H3x2z_D2x_Pz_S;
  Double I_ERI_H2x2yz_F2xz_Pz_S = I_ERI_I2x2y2z_D2x_Pz_S+ABZ*I_ERI_H2x2yz_D2x_Pz_S;
  Double I_ERI_H2xy2z_F2xz_Pz_S = I_ERI_I2xy3z_D2x_Pz_S+ABZ*I_ERI_H2xy2z_D2x_Pz_S;
  Double I_ERI_H2x3z_F2xz_Pz_S = I_ERI_I2x4z_D2x_Pz_S+ABZ*I_ERI_H2x3z_D2x_Pz_S;
  Double I_ERI_Hx3yz_F2xz_Pz_S = I_ERI_Ix3y2z_D2x_Pz_S+ABZ*I_ERI_Hx3yz_D2x_Pz_S;
  Double I_ERI_Hx2y2z_F2xz_Pz_S = I_ERI_Ix2y3z_D2x_Pz_S+ABZ*I_ERI_Hx2y2z_D2x_Pz_S;
  Double I_ERI_Hxy3z_F2xz_Pz_S = I_ERI_Ixy4z_D2x_Pz_S+ABZ*I_ERI_Hxy3z_D2x_Pz_S;
  Double I_ERI_Hx4z_F2xz_Pz_S = I_ERI_Ix5z_D2x_Pz_S+ABZ*I_ERI_Hx4z_D2x_Pz_S;
  Double I_ERI_H4yz_F2xz_Pz_S = I_ERI_I4y2z_D2x_Pz_S+ABZ*I_ERI_H4yz_D2x_Pz_S;
  Double I_ERI_H3y2z_F2xz_Pz_S = I_ERI_I3y3z_D2x_Pz_S+ABZ*I_ERI_H3y2z_D2x_Pz_S;
  Double I_ERI_H2y3z_F2xz_Pz_S = I_ERI_I2y4z_D2x_Pz_S+ABZ*I_ERI_H2y3z_D2x_Pz_S;
  Double I_ERI_Hy4z_F2xz_Pz_S = I_ERI_Iy5z_D2x_Pz_S+ABZ*I_ERI_Hy4z_D2x_Pz_S;
  Double I_ERI_H5z_F2xz_Pz_S = I_ERI_I6z_D2x_Pz_S+ABZ*I_ERI_H5z_D2x_Pz_S;
  Double I_ERI_H4xz_Fx2y_Pz_S = I_ERI_I5xz_D2y_Pz_S+ABX*I_ERI_H4xz_D2y_Pz_S;
  Double I_ERI_H3xyz_Fx2y_Pz_S = I_ERI_I4xyz_D2y_Pz_S+ABX*I_ERI_H3xyz_D2y_Pz_S;
  Double I_ERI_H3x2z_Fx2y_Pz_S = I_ERI_I4x2z_D2y_Pz_S+ABX*I_ERI_H3x2z_D2y_Pz_S;
  Double I_ERI_H2x2yz_Fx2y_Pz_S = I_ERI_I3x2yz_D2y_Pz_S+ABX*I_ERI_H2x2yz_D2y_Pz_S;
  Double I_ERI_H2xy2z_Fx2y_Pz_S = I_ERI_I3xy2z_D2y_Pz_S+ABX*I_ERI_H2xy2z_D2y_Pz_S;
  Double I_ERI_H2x3z_Fx2y_Pz_S = I_ERI_I3x3z_D2y_Pz_S+ABX*I_ERI_H2x3z_D2y_Pz_S;
  Double I_ERI_Hx3yz_Fx2y_Pz_S = I_ERI_I2x3yz_D2y_Pz_S+ABX*I_ERI_Hx3yz_D2y_Pz_S;
  Double I_ERI_Hx2y2z_Fx2y_Pz_S = I_ERI_I2x2y2z_D2y_Pz_S+ABX*I_ERI_Hx2y2z_D2y_Pz_S;
  Double I_ERI_Hxy3z_Fx2y_Pz_S = I_ERI_I2xy3z_D2y_Pz_S+ABX*I_ERI_Hxy3z_D2y_Pz_S;
  Double I_ERI_Hx4z_Fx2y_Pz_S = I_ERI_I2x4z_D2y_Pz_S+ABX*I_ERI_Hx4z_D2y_Pz_S;
  Double I_ERI_H4yz_Fx2y_Pz_S = I_ERI_Ix4yz_D2y_Pz_S+ABX*I_ERI_H4yz_D2y_Pz_S;
  Double I_ERI_H3y2z_Fx2y_Pz_S = I_ERI_Ix3y2z_D2y_Pz_S+ABX*I_ERI_H3y2z_D2y_Pz_S;
  Double I_ERI_H2y3z_Fx2y_Pz_S = I_ERI_Ix2y3z_D2y_Pz_S+ABX*I_ERI_H2y3z_D2y_Pz_S;
  Double I_ERI_Hy4z_Fx2y_Pz_S = I_ERI_Ixy4z_D2y_Pz_S+ABX*I_ERI_Hy4z_D2y_Pz_S;
  Double I_ERI_H5z_Fx2y_Pz_S = I_ERI_Ix5z_D2y_Pz_S+ABX*I_ERI_H5z_D2y_Pz_S;
  Double I_ERI_H4xy_Fx2z_Pz_S = I_ERI_I5xy_D2z_Pz_S+ABX*I_ERI_H4xy_D2z_Pz_S;
  Double I_ERI_H3x2y_Fx2z_Pz_S = I_ERI_I4x2y_D2z_Pz_S+ABX*I_ERI_H3x2y_D2z_Pz_S;
  Double I_ERI_H3xyz_Fx2z_Pz_S = I_ERI_I4xyz_D2z_Pz_S+ABX*I_ERI_H3xyz_D2z_Pz_S;
  Double I_ERI_H2x3y_Fx2z_Pz_S = I_ERI_I3x3y_D2z_Pz_S+ABX*I_ERI_H2x3y_D2z_Pz_S;
  Double I_ERI_H2x2yz_Fx2z_Pz_S = I_ERI_I3x2yz_D2z_Pz_S+ABX*I_ERI_H2x2yz_D2z_Pz_S;
  Double I_ERI_H2xy2z_Fx2z_Pz_S = I_ERI_I3xy2z_D2z_Pz_S+ABX*I_ERI_H2xy2z_D2z_Pz_S;
  Double I_ERI_Hx4y_Fx2z_Pz_S = I_ERI_I2x4y_D2z_Pz_S+ABX*I_ERI_Hx4y_D2z_Pz_S;
  Double I_ERI_Hx3yz_Fx2z_Pz_S = I_ERI_I2x3yz_D2z_Pz_S+ABX*I_ERI_Hx3yz_D2z_Pz_S;
  Double I_ERI_Hx2y2z_Fx2z_Pz_S = I_ERI_I2x2y2z_D2z_Pz_S+ABX*I_ERI_Hx2y2z_D2z_Pz_S;
  Double I_ERI_Hxy3z_Fx2z_Pz_S = I_ERI_I2xy3z_D2z_Pz_S+ABX*I_ERI_Hxy3z_D2z_Pz_S;
  Double I_ERI_H5y_Fx2z_Pz_S = I_ERI_Ix5y_D2z_Pz_S+ABX*I_ERI_H5y_D2z_Pz_S;
  Double I_ERI_H4yz_Fx2z_Pz_S = I_ERI_Ix4yz_D2z_Pz_S+ABX*I_ERI_H4yz_D2z_Pz_S;
  Double I_ERI_H3y2z_Fx2z_Pz_S = I_ERI_Ix3y2z_D2z_Pz_S+ABX*I_ERI_H3y2z_D2z_Pz_S;
  Double I_ERI_H2y3z_Fx2z_Pz_S = I_ERI_Ix2y3z_D2z_Pz_S+ABX*I_ERI_H2y3z_D2z_Pz_S;
  Double I_ERI_Hy4z_Fx2z_Pz_S = I_ERI_Ixy4z_D2z_Pz_S+ABX*I_ERI_Hy4z_D2z_Pz_S;
  Double I_ERI_H5x_F3y_Pz_S = I_ERI_I5xy_D2y_Pz_S+ABY*I_ERI_H5x_D2y_Pz_S;
  Double I_ERI_H4xy_F3y_Pz_S = I_ERI_I4x2y_D2y_Pz_S+ABY*I_ERI_H4xy_D2y_Pz_S;
  Double I_ERI_H4xz_F3y_Pz_S = I_ERI_I4xyz_D2y_Pz_S+ABY*I_ERI_H4xz_D2y_Pz_S;
  Double I_ERI_H3x2y_F3y_Pz_S = I_ERI_I3x3y_D2y_Pz_S+ABY*I_ERI_H3x2y_D2y_Pz_S;
  Double I_ERI_H3xyz_F3y_Pz_S = I_ERI_I3x2yz_D2y_Pz_S+ABY*I_ERI_H3xyz_D2y_Pz_S;
  Double I_ERI_H3x2z_F3y_Pz_S = I_ERI_I3xy2z_D2y_Pz_S+ABY*I_ERI_H3x2z_D2y_Pz_S;
  Double I_ERI_H2x3y_F3y_Pz_S = I_ERI_I2x4y_D2y_Pz_S+ABY*I_ERI_H2x3y_D2y_Pz_S;
  Double I_ERI_H2x2yz_F3y_Pz_S = I_ERI_I2x3yz_D2y_Pz_S+ABY*I_ERI_H2x2yz_D2y_Pz_S;
  Double I_ERI_H2xy2z_F3y_Pz_S = I_ERI_I2x2y2z_D2y_Pz_S+ABY*I_ERI_H2xy2z_D2y_Pz_S;
  Double I_ERI_H2x3z_F3y_Pz_S = I_ERI_I2xy3z_D2y_Pz_S+ABY*I_ERI_H2x3z_D2y_Pz_S;
  Double I_ERI_Hx4y_F3y_Pz_S = I_ERI_Ix5y_D2y_Pz_S+ABY*I_ERI_Hx4y_D2y_Pz_S;
  Double I_ERI_Hx3yz_F3y_Pz_S = I_ERI_Ix4yz_D2y_Pz_S+ABY*I_ERI_Hx3yz_D2y_Pz_S;
  Double I_ERI_Hx2y2z_F3y_Pz_S = I_ERI_Ix3y2z_D2y_Pz_S+ABY*I_ERI_Hx2y2z_D2y_Pz_S;
  Double I_ERI_Hxy3z_F3y_Pz_S = I_ERI_Ix2y3z_D2y_Pz_S+ABY*I_ERI_Hxy3z_D2y_Pz_S;
  Double I_ERI_Hx4z_F3y_Pz_S = I_ERI_Ixy4z_D2y_Pz_S+ABY*I_ERI_Hx4z_D2y_Pz_S;
  Double I_ERI_H5y_F3y_Pz_S = I_ERI_I6y_D2y_Pz_S+ABY*I_ERI_H5y_D2y_Pz_S;
  Double I_ERI_H4yz_F3y_Pz_S = I_ERI_I5yz_D2y_Pz_S+ABY*I_ERI_H4yz_D2y_Pz_S;
  Double I_ERI_H3y2z_F3y_Pz_S = I_ERI_I4y2z_D2y_Pz_S+ABY*I_ERI_H3y2z_D2y_Pz_S;
  Double I_ERI_H2y3z_F3y_Pz_S = I_ERI_I3y3z_D2y_Pz_S+ABY*I_ERI_H2y3z_D2y_Pz_S;
  Double I_ERI_Hy4z_F3y_Pz_S = I_ERI_I2y4z_D2y_Pz_S+ABY*I_ERI_Hy4z_D2y_Pz_S;
  Double I_ERI_H5z_F3y_Pz_S = I_ERI_Iy5z_D2y_Pz_S+ABY*I_ERI_H5z_D2y_Pz_S;
  Double I_ERI_H4xz_F2yz_Pz_S = I_ERI_I4x2z_D2y_Pz_S+ABZ*I_ERI_H4xz_D2y_Pz_S;
  Double I_ERI_H3xyz_F2yz_Pz_S = I_ERI_I3xy2z_D2y_Pz_S+ABZ*I_ERI_H3xyz_D2y_Pz_S;
  Double I_ERI_H3x2z_F2yz_Pz_S = I_ERI_I3x3z_D2y_Pz_S+ABZ*I_ERI_H3x2z_D2y_Pz_S;
  Double I_ERI_H2x2yz_F2yz_Pz_S = I_ERI_I2x2y2z_D2y_Pz_S+ABZ*I_ERI_H2x2yz_D2y_Pz_S;
  Double I_ERI_H2xy2z_F2yz_Pz_S = I_ERI_I2xy3z_D2y_Pz_S+ABZ*I_ERI_H2xy2z_D2y_Pz_S;
  Double I_ERI_H2x3z_F2yz_Pz_S = I_ERI_I2x4z_D2y_Pz_S+ABZ*I_ERI_H2x3z_D2y_Pz_S;
  Double I_ERI_Hx3yz_F2yz_Pz_S = I_ERI_Ix3y2z_D2y_Pz_S+ABZ*I_ERI_Hx3yz_D2y_Pz_S;
  Double I_ERI_Hx2y2z_F2yz_Pz_S = I_ERI_Ix2y3z_D2y_Pz_S+ABZ*I_ERI_Hx2y2z_D2y_Pz_S;
  Double I_ERI_Hxy3z_F2yz_Pz_S = I_ERI_Ixy4z_D2y_Pz_S+ABZ*I_ERI_Hxy3z_D2y_Pz_S;
  Double I_ERI_Hx4z_F2yz_Pz_S = I_ERI_Ix5z_D2y_Pz_S+ABZ*I_ERI_Hx4z_D2y_Pz_S;
  Double I_ERI_H4yz_F2yz_Pz_S = I_ERI_I4y2z_D2y_Pz_S+ABZ*I_ERI_H4yz_D2y_Pz_S;
  Double I_ERI_H3y2z_F2yz_Pz_S = I_ERI_I3y3z_D2y_Pz_S+ABZ*I_ERI_H3y2z_D2y_Pz_S;
  Double I_ERI_H2y3z_F2yz_Pz_S = I_ERI_I2y4z_D2y_Pz_S+ABZ*I_ERI_H2y3z_D2y_Pz_S;
  Double I_ERI_Hy4z_F2yz_Pz_S = I_ERI_Iy5z_D2y_Pz_S+ABZ*I_ERI_Hy4z_D2y_Pz_S;
  Double I_ERI_H5z_F2yz_Pz_S = I_ERI_I6z_D2y_Pz_S+ABZ*I_ERI_H5z_D2y_Pz_S;
  Double I_ERI_H5x_F3z_Pz_S = I_ERI_I5xz_D2z_Pz_S+ABZ*I_ERI_H5x_D2z_Pz_S;
  Double I_ERI_H4xy_F3z_Pz_S = I_ERI_I4xyz_D2z_Pz_S+ABZ*I_ERI_H4xy_D2z_Pz_S;
  Double I_ERI_H4xz_F3z_Pz_S = I_ERI_I4x2z_D2z_Pz_S+ABZ*I_ERI_H4xz_D2z_Pz_S;
  Double I_ERI_H3x2y_F3z_Pz_S = I_ERI_I3x2yz_D2z_Pz_S+ABZ*I_ERI_H3x2y_D2z_Pz_S;
  Double I_ERI_H3xyz_F3z_Pz_S = I_ERI_I3xy2z_D2z_Pz_S+ABZ*I_ERI_H3xyz_D2z_Pz_S;
  Double I_ERI_H3x2z_F3z_Pz_S = I_ERI_I3x3z_D2z_Pz_S+ABZ*I_ERI_H3x2z_D2z_Pz_S;
  Double I_ERI_H2x3y_F3z_Pz_S = I_ERI_I2x3yz_D2z_Pz_S+ABZ*I_ERI_H2x3y_D2z_Pz_S;
  Double I_ERI_H2x2yz_F3z_Pz_S = I_ERI_I2x2y2z_D2z_Pz_S+ABZ*I_ERI_H2x2yz_D2z_Pz_S;
  Double I_ERI_H2xy2z_F3z_Pz_S = I_ERI_I2xy3z_D2z_Pz_S+ABZ*I_ERI_H2xy2z_D2z_Pz_S;
  Double I_ERI_H2x3z_F3z_Pz_S = I_ERI_I2x4z_D2z_Pz_S+ABZ*I_ERI_H2x3z_D2z_Pz_S;
  Double I_ERI_Hx4y_F3z_Pz_S = I_ERI_Ix4yz_D2z_Pz_S+ABZ*I_ERI_Hx4y_D2z_Pz_S;
  Double I_ERI_Hx3yz_F3z_Pz_S = I_ERI_Ix3y2z_D2z_Pz_S+ABZ*I_ERI_Hx3yz_D2z_Pz_S;
  Double I_ERI_Hx2y2z_F3z_Pz_S = I_ERI_Ix2y3z_D2z_Pz_S+ABZ*I_ERI_Hx2y2z_D2z_Pz_S;
  Double I_ERI_Hxy3z_F3z_Pz_S = I_ERI_Ixy4z_D2z_Pz_S+ABZ*I_ERI_Hxy3z_D2z_Pz_S;
  Double I_ERI_Hx4z_F3z_Pz_S = I_ERI_Ix5z_D2z_Pz_S+ABZ*I_ERI_Hx4z_D2z_Pz_S;
  Double I_ERI_H5y_F3z_Pz_S = I_ERI_I5yz_D2z_Pz_S+ABZ*I_ERI_H5y_D2z_Pz_S;
  Double I_ERI_H4yz_F3z_Pz_S = I_ERI_I4y2z_D2z_Pz_S+ABZ*I_ERI_H4yz_D2z_Pz_S;
  Double I_ERI_H3y2z_F3z_Pz_S = I_ERI_I3y3z_D2z_Pz_S+ABZ*I_ERI_H3y2z_D2z_Pz_S;
  Double I_ERI_H2y3z_F3z_Pz_S = I_ERI_I2y4z_D2z_Pz_S+ABZ*I_ERI_H2y3z_D2z_Pz_S;
  Double I_ERI_Hy4z_F3z_Pz_S = I_ERI_Iy5z_D2z_Pz_S+ABZ*I_ERI_Hy4z_D2z_Pz_S;
  Double I_ERI_H5z_F3z_Pz_S = I_ERI_I6z_D2z_Pz_S+ABZ*I_ERI_H5z_D2z_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_G_G_P_S
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_F_P_S
   * RHS shell quartet name: SQ_ERI_G_F_P_S
   ************************************************************/
  abcd[0] = I_ERI_H5x_F3x_Px_S+ABX*I_ERI_G4x_F3x_Px_S;
  abcd[1] = I_ERI_H4xy_F3x_Px_S+ABX*I_ERI_G3xy_F3x_Px_S;
  abcd[2] = I_ERI_H4xz_F3x_Px_S+ABX*I_ERI_G3xz_F3x_Px_S;
  abcd[3] = I_ERI_H3x2y_F3x_Px_S+ABX*I_ERI_G2x2y_F3x_Px_S;
  abcd[4] = I_ERI_H3xyz_F3x_Px_S+ABX*I_ERI_G2xyz_F3x_Px_S;
  abcd[5] = I_ERI_H3x2z_F3x_Px_S+ABX*I_ERI_G2x2z_F3x_Px_S;
  abcd[6] = I_ERI_H2x3y_F3x_Px_S+ABX*I_ERI_Gx3y_F3x_Px_S;
  abcd[7] = I_ERI_H2x2yz_F3x_Px_S+ABX*I_ERI_Gx2yz_F3x_Px_S;
  abcd[8] = I_ERI_H2xy2z_F3x_Px_S+ABX*I_ERI_Gxy2z_F3x_Px_S;
  abcd[9] = I_ERI_H2x3z_F3x_Px_S+ABX*I_ERI_Gx3z_F3x_Px_S;
  abcd[10] = I_ERI_Hx4y_F3x_Px_S+ABX*I_ERI_G4y_F3x_Px_S;
  abcd[11] = I_ERI_Hx3yz_F3x_Px_S+ABX*I_ERI_G3yz_F3x_Px_S;
  abcd[12] = I_ERI_Hx2y2z_F3x_Px_S+ABX*I_ERI_G2y2z_F3x_Px_S;
  abcd[13] = I_ERI_Hxy3z_F3x_Px_S+ABX*I_ERI_Gy3z_F3x_Px_S;
  abcd[14] = I_ERI_Hx4z_F3x_Px_S+ABX*I_ERI_G4z_F3x_Px_S;
  abcd[15] = I_ERI_H4xy_F3x_Px_S+ABY*I_ERI_G4x_F3x_Px_S;
  abcd[16] = I_ERI_H3x2y_F3x_Px_S+ABY*I_ERI_G3xy_F3x_Px_S;
  abcd[17] = I_ERI_H3xyz_F3x_Px_S+ABY*I_ERI_G3xz_F3x_Px_S;
  abcd[18] = I_ERI_H2x3y_F3x_Px_S+ABY*I_ERI_G2x2y_F3x_Px_S;
  abcd[19] = I_ERI_H2x2yz_F3x_Px_S+ABY*I_ERI_G2xyz_F3x_Px_S;
  abcd[20] = I_ERI_H2xy2z_F3x_Px_S+ABY*I_ERI_G2x2z_F3x_Px_S;
  abcd[21] = I_ERI_Hx4y_F3x_Px_S+ABY*I_ERI_Gx3y_F3x_Px_S;
  abcd[22] = I_ERI_Hx3yz_F3x_Px_S+ABY*I_ERI_Gx2yz_F3x_Px_S;
  abcd[23] = I_ERI_Hx2y2z_F3x_Px_S+ABY*I_ERI_Gxy2z_F3x_Px_S;
  abcd[24] = I_ERI_Hxy3z_F3x_Px_S+ABY*I_ERI_Gx3z_F3x_Px_S;
  abcd[25] = I_ERI_H5y_F3x_Px_S+ABY*I_ERI_G4y_F3x_Px_S;
  abcd[26] = I_ERI_H4yz_F3x_Px_S+ABY*I_ERI_G3yz_F3x_Px_S;
  abcd[27] = I_ERI_H3y2z_F3x_Px_S+ABY*I_ERI_G2y2z_F3x_Px_S;
  abcd[28] = I_ERI_H2y3z_F3x_Px_S+ABY*I_ERI_Gy3z_F3x_Px_S;
  abcd[29] = I_ERI_Hy4z_F3x_Px_S+ABY*I_ERI_G4z_F3x_Px_S;
  abcd[30] = I_ERI_H4xz_F3x_Px_S+ABZ*I_ERI_G4x_F3x_Px_S;
  abcd[31] = I_ERI_H3xyz_F3x_Px_S+ABZ*I_ERI_G3xy_F3x_Px_S;
  abcd[32] = I_ERI_H3x2z_F3x_Px_S+ABZ*I_ERI_G3xz_F3x_Px_S;
  abcd[33] = I_ERI_H2x2yz_F3x_Px_S+ABZ*I_ERI_G2x2y_F3x_Px_S;
  abcd[34] = I_ERI_H2xy2z_F3x_Px_S+ABZ*I_ERI_G2xyz_F3x_Px_S;
  abcd[35] = I_ERI_H2x3z_F3x_Px_S+ABZ*I_ERI_G2x2z_F3x_Px_S;
  abcd[36] = I_ERI_Hx3yz_F3x_Px_S+ABZ*I_ERI_Gx3y_F3x_Px_S;
  abcd[37] = I_ERI_Hx2y2z_F3x_Px_S+ABZ*I_ERI_Gx2yz_F3x_Px_S;
  abcd[38] = I_ERI_Hxy3z_F3x_Px_S+ABZ*I_ERI_Gxy2z_F3x_Px_S;
  abcd[39] = I_ERI_Hx4z_F3x_Px_S+ABZ*I_ERI_Gx3z_F3x_Px_S;
  abcd[40] = I_ERI_H4yz_F3x_Px_S+ABZ*I_ERI_G4y_F3x_Px_S;
  abcd[41] = I_ERI_H3y2z_F3x_Px_S+ABZ*I_ERI_G3yz_F3x_Px_S;
  abcd[42] = I_ERI_H2y3z_F3x_Px_S+ABZ*I_ERI_G2y2z_F3x_Px_S;
  abcd[43] = I_ERI_Hy4z_F3x_Px_S+ABZ*I_ERI_Gy3z_F3x_Px_S;
  abcd[44] = I_ERI_H5z_F3x_Px_S+ABZ*I_ERI_G4z_F3x_Px_S;
  abcd[45] = I_ERI_H4xy_F2xy_Px_S+ABY*I_ERI_G4x_F2xy_Px_S;
  abcd[46] = I_ERI_H3x2y_F2xy_Px_S+ABY*I_ERI_G3xy_F2xy_Px_S;
  abcd[47] = I_ERI_H3xyz_F2xy_Px_S+ABY*I_ERI_G3xz_F2xy_Px_S;
  abcd[48] = I_ERI_H2x3y_F2xy_Px_S+ABY*I_ERI_G2x2y_F2xy_Px_S;
  abcd[49] = I_ERI_H2x2yz_F2xy_Px_S+ABY*I_ERI_G2xyz_F2xy_Px_S;
  abcd[50] = I_ERI_H2xy2z_F2xy_Px_S+ABY*I_ERI_G2x2z_F2xy_Px_S;
  abcd[51] = I_ERI_Hx4y_F2xy_Px_S+ABY*I_ERI_Gx3y_F2xy_Px_S;
  abcd[52] = I_ERI_Hx3yz_F2xy_Px_S+ABY*I_ERI_Gx2yz_F2xy_Px_S;
  abcd[53] = I_ERI_Hx2y2z_F2xy_Px_S+ABY*I_ERI_Gxy2z_F2xy_Px_S;
  abcd[54] = I_ERI_Hxy3z_F2xy_Px_S+ABY*I_ERI_Gx3z_F2xy_Px_S;
  abcd[55] = I_ERI_H5y_F2xy_Px_S+ABY*I_ERI_G4y_F2xy_Px_S;
  abcd[56] = I_ERI_H4yz_F2xy_Px_S+ABY*I_ERI_G3yz_F2xy_Px_S;
  abcd[57] = I_ERI_H3y2z_F2xy_Px_S+ABY*I_ERI_G2y2z_F2xy_Px_S;
  abcd[58] = I_ERI_H2y3z_F2xy_Px_S+ABY*I_ERI_Gy3z_F2xy_Px_S;
  abcd[59] = I_ERI_Hy4z_F2xy_Px_S+ABY*I_ERI_G4z_F2xy_Px_S;
  abcd[60] = I_ERI_H4xz_F2xy_Px_S+ABZ*I_ERI_G4x_F2xy_Px_S;
  abcd[61] = I_ERI_H3xyz_F2xy_Px_S+ABZ*I_ERI_G3xy_F2xy_Px_S;
  abcd[62] = I_ERI_H3x2z_F2xy_Px_S+ABZ*I_ERI_G3xz_F2xy_Px_S;
  abcd[63] = I_ERI_H2x2yz_F2xy_Px_S+ABZ*I_ERI_G2x2y_F2xy_Px_S;
  abcd[64] = I_ERI_H2xy2z_F2xy_Px_S+ABZ*I_ERI_G2xyz_F2xy_Px_S;
  abcd[65] = I_ERI_H2x3z_F2xy_Px_S+ABZ*I_ERI_G2x2z_F2xy_Px_S;
  abcd[66] = I_ERI_Hx3yz_F2xy_Px_S+ABZ*I_ERI_Gx3y_F2xy_Px_S;
  abcd[67] = I_ERI_Hx2y2z_F2xy_Px_S+ABZ*I_ERI_Gx2yz_F2xy_Px_S;
  abcd[68] = I_ERI_Hxy3z_F2xy_Px_S+ABZ*I_ERI_Gxy2z_F2xy_Px_S;
  abcd[69] = I_ERI_Hx4z_F2xy_Px_S+ABZ*I_ERI_Gx3z_F2xy_Px_S;
  abcd[70] = I_ERI_H4yz_F2xy_Px_S+ABZ*I_ERI_G4y_F2xy_Px_S;
  abcd[71] = I_ERI_H3y2z_F2xy_Px_S+ABZ*I_ERI_G3yz_F2xy_Px_S;
  abcd[72] = I_ERI_H2y3z_F2xy_Px_S+ABZ*I_ERI_G2y2z_F2xy_Px_S;
  abcd[73] = I_ERI_Hy4z_F2xy_Px_S+ABZ*I_ERI_Gy3z_F2xy_Px_S;
  abcd[74] = I_ERI_H5z_F2xy_Px_S+ABZ*I_ERI_G4z_F2xy_Px_S;
  abcd[75] = I_ERI_H4xz_F2xz_Px_S+ABZ*I_ERI_G4x_F2xz_Px_S;
  abcd[76] = I_ERI_H3xyz_F2xz_Px_S+ABZ*I_ERI_G3xy_F2xz_Px_S;
  abcd[77] = I_ERI_H3x2z_F2xz_Px_S+ABZ*I_ERI_G3xz_F2xz_Px_S;
  abcd[78] = I_ERI_H2x2yz_F2xz_Px_S+ABZ*I_ERI_G2x2y_F2xz_Px_S;
  abcd[79] = I_ERI_H2xy2z_F2xz_Px_S+ABZ*I_ERI_G2xyz_F2xz_Px_S;
  abcd[80] = I_ERI_H2x3z_F2xz_Px_S+ABZ*I_ERI_G2x2z_F2xz_Px_S;
  abcd[81] = I_ERI_Hx3yz_F2xz_Px_S+ABZ*I_ERI_Gx3y_F2xz_Px_S;
  abcd[82] = I_ERI_Hx2y2z_F2xz_Px_S+ABZ*I_ERI_Gx2yz_F2xz_Px_S;
  abcd[83] = I_ERI_Hxy3z_F2xz_Px_S+ABZ*I_ERI_Gxy2z_F2xz_Px_S;
  abcd[84] = I_ERI_Hx4z_F2xz_Px_S+ABZ*I_ERI_Gx3z_F2xz_Px_S;
  abcd[85] = I_ERI_H4yz_F2xz_Px_S+ABZ*I_ERI_G4y_F2xz_Px_S;
  abcd[86] = I_ERI_H3y2z_F2xz_Px_S+ABZ*I_ERI_G3yz_F2xz_Px_S;
  abcd[87] = I_ERI_H2y3z_F2xz_Px_S+ABZ*I_ERI_G2y2z_F2xz_Px_S;
  abcd[88] = I_ERI_Hy4z_F2xz_Px_S+ABZ*I_ERI_Gy3z_F2xz_Px_S;
  abcd[89] = I_ERI_H5z_F2xz_Px_S+ABZ*I_ERI_G4z_F2xz_Px_S;
  abcd[90] = I_ERI_H5x_F3y_Px_S+ABX*I_ERI_G4x_F3y_Px_S;
  abcd[91] = I_ERI_H4xy_F3y_Px_S+ABX*I_ERI_G3xy_F3y_Px_S;
  abcd[92] = I_ERI_H4xz_F3y_Px_S+ABX*I_ERI_G3xz_F3y_Px_S;
  abcd[93] = I_ERI_H3x2y_F3y_Px_S+ABX*I_ERI_G2x2y_F3y_Px_S;
  abcd[94] = I_ERI_H3xyz_F3y_Px_S+ABX*I_ERI_G2xyz_F3y_Px_S;
  abcd[95] = I_ERI_H3x2z_F3y_Px_S+ABX*I_ERI_G2x2z_F3y_Px_S;
  abcd[96] = I_ERI_H2x3y_F3y_Px_S+ABX*I_ERI_Gx3y_F3y_Px_S;
  abcd[97] = I_ERI_H2x2yz_F3y_Px_S+ABX*I_ERI_Gx2yz_F3y_Px_S;
  abcd[98] = I_ERI_H2xy2z_F3y_Px_S+ABX*I_ERI_Gxy2z_F3y_Px_S;
  abcd[99] = I_ERI_H2x3z_F3y_Px_S+ABX*I_ERI_Gx3z_F3y_Px_S;
  abcd[100] = I_ERI_Hx4y_F3y_Px_S+ABX*I_ERI_G4y_F3y_Px_S;
  abcd[101] = I_ERI_Hx3yz_F3y_Px_S+ABX*I_ERI_G3yz_F3y_Px_S;
  abcd[102] = I_ERI_Hx2y2z_F3y_Px_S+ABX*I_ERI_G2y2z_F3y_Px_S;
  abcd[103] = I_ERI_Hxy3z_F3y_Px_S+ABX*I_ERI_Gy3z_F3y_Px_S;
  abcd[104] = I_ERI_Hx4z_F3y_Px_S+ABX*I_ERI_G4z_F3y_Px_S;
  abcd[105] = I_ERI_H4xz_Fx2y_Px_S+ABZ*I_ERI_G4x_Fx2y_Px_S;
  abcd[106] = I_ERI_H3xyz_Fx2y_Px_S+ABZ*I_ERI_G3xy_Fx2y_Px_S;
  abcd[107] = I_ERI_H3x2z_Fx2y_Px_S+ABZ*I_ERI_G3xz_Fx2y_Px_S;
  abcd[108] = I_ERI_H2x2yz_Fx2y_Px_S+ABZ*I_ERI_G2x2y_Fx2y_Px_S;
  abcd[109] = I_ERI_H2xy2z_Fx2y_Px_S+ABZ*I_ERI_G2xyz_Fx2y_Px_S;
  abcd[110] = I_ERI_H2x3z_Fx2y_Px_S+ABZ*I_ERI_G2x2z_Fx2y_Px_S;
  abcd[111] = I_ERI_Hx3yz_Fx2y_Px_S+ABZ*I_ERI_Gx3y_Fx2y_Px_S;
  abcd[112] = I_ERI_Hx2y2z_Fx2y_Px_S+ABZ*I_ERI_Gx2yz_Fx2y_Px_S;
  abcd[113] = I_ERI_Hxy3z_Fx2y_Px_S+ABZ*I_ERI_Gxy2z_Fx2y_Px_S;
  abcd[114] = I_ERI_Hx4z_Fx2y_Px_S+ABZ*I_ERI_Gx3z_Fx2y_Px_S;
  abcd[115] = I_ERI_H4yz_Fx2y_Px_S+ABZ*I_ERI_G4y_Fx2y_Px_S;
  abcd[116] = I_ERI_H3y2z_Fx2y_Px_S+ABZ*I_ERI_G3yz_Fx2y_Px_S;
  abcd[117] = I_ERI_H2y3z_Fx2y_Px_S+ABZ*I_ERI_G2y2z_Fx2y_Px_S;
  abcd[118] = I_ERI_Hy4z_Fx2y_Px_S+ABZ*I_ERI_Gy3z_Fx2y_Px_S;
  abcd[119] = I_ERI_H5z_Fx2y_Px_S+ABZ*I_ERI_G4z_Fx2y_Px_S;
  abcd[120] = I_ERI_H4xy_Fx2z_Px_S+ABY*I_ERI_G4x_Fx2z_Px_S;
  abcd[121] = I_ERI_H3x2y_Fx2z_Px_S+ABY*I_ERI_G3xy_Fx2z_Px_S;
  abcd[122] = I_ERI_H3xyz_Fx2z_Px_S+ABY*I_ERI_G3xz_Fx2z_Px_S;
  abcd[123] = I_ERI_H2x3y_Fx2z_Px_S+ABY*I_ERI_G2x2y_Fx2z_Px_S;
  abcd[124] = I_ERI_H2x2yz_Fx2z_Px_S+ABY*I_ERI_G2xyz_Fx2z_Px_S;
  abcd[125] = I_ERI_H2xy2z_Fx2z_Px_S+ABY*I_ERI_G2x2z_Fx2z_Px_S;
  abcd[126] = I_ERI_Hx4y_Fx2z_Px_S+ABY*I_ERI_Gx3y_Fx2z_Px_S;
  abcd[127] = I_ERI_Hx3yz_Fx2z_Px_S+ABY*I_ERI_Gx2yz_Fx2z_Px_S;
  abcd[128] = I_ERI_Hx2y2z_Fx2z_Px_S+ABY*I_ERI_Gxy2z_Fx2z_Px_S;
  abcd[129] = I_ERI_Hxy3z_Fx2z_Px_S+ABY*I_ERI_Gx3z_Fx2z_Px_S;
  abcd[130] = I_ERI_H5y_Fx2z_Px_S+ABY*I_ERI_G4y_Fx2z_Px_S;
  abcd[131] = I_ERI_H4yz_Fx2z_Px_S+ABY*I_ERI_G3yz_Fx2z_Px_S;
  abcd[132] = I_ERI_H3y2z_Fx2z_Px_S+ABY*I_ERI_G2y2z_Fx2z_Px_S;
  abcd[133] = I_ERI_H2y3z_Fx2z_Px_S+ABY*I_ERI_Gy3z_Fx2z_Px_S;
  abcd[134] = I_ERI_Hy4z_Fx2z_Px_S+ABY*I_ERI_G4z_Fx2z_Px_S;
  abcd[135] = I_ERI_H5x_F3z_Px_S+ABX*I_ERI_G4x_F3z_Px_S;
  abcd[136] = I_ERI_H4xy_F3z_Px_S+ABX*I_ERI_G3xy_F3z_Px_S;
  abcd[137] = I_ERI_H4xz_F3z_Px_S+ABX*I_ERI_G3xz_F3z_Px_S;
  abcd[138] = I_ERI_H3x2y_F3z_Px_S+ABX*I_ERI_G2x2y_F3z_Px_S;
  abcd[139] = I_ERI_H3xyz_F3z_Px_S+ABX*I_ERI_G2xyz_F3z_Px_S;
  abcd[140] = I_ERI_H3x2z_F3z_Px_S+ABX*I_ERI_G2x2z_F3z_Px_S;
  abcd[141] = I_ERI_H2x3y_F3z_Px_S+ABX*I_ERI_Gx3y_F3z_Px_S;
  abcd[142] = I_ERI_H2x2yz_F3z_Px_S+ABX*I_ERI_Gx2yz_F3z_Px_S;
  abcd[143] = I_ERI_H2xy2z_F3z_Px_S+ABX*I_ERI_Gxy2z_F3z_Px_S;
  abcd[144] = I_ERI_H2x3z_F3z_Px_S+ABX*I_ERI_Gx3z_F3z_Px_S;
  abcd[145] = I_ERI_Hx4y_F3z_Px_S+ABX*I_ERI_G4y_F3z_Px_S;
  abcd[146] = I_ERI_Hx3yz_F3z_Px_S+ABX*I_ERI_G3yz_F3z_Px_S;
  abcd[147] = I_ERI_Hx2y2z_F3z_Px_S+ABX*I_ERI_G2y2z_F3z_Px_S;
  abcd[148] = I_ERI_Hxy3z_F3z_Px_S+ABX*I_ERI_Gy3z_F3z_Px_S;
  abcd[149] = I_ERI_Hx4z_F3z_Px_S+ABX*I_ERI_G4z_F3z_Px_S;
  abcd[150] = I_ERI_H4xy_F3y_Px_S+ABY*I_ERI_G4x_F3y_Px_S;
  abcd[151] = I_ERI_H3x2y_F3y_Px_S+ABY*I_ERI_G3xy_F3y_Px_S;
  abcd[152] = I_ERI_H3xyz_F3y_Px_S+ABY*I_ERI_G3xz_F3y_Px_S;
  abcd[153] = I_ERI_H2x3y_F3y_Px_S+ABY*I_ERI_G2x2y_F3y_Px_S;
  abcd[154] = I_ERI_H2x2yz_F3y_Px_S+ABY*I_ERI_G2xyz_F3y_Px_S;
  abcd[155] = I_ERI_H2xy2z_F3y_Px_S+ABY*I_ERI_G2x2z_F3y_Px_S;
  abcd[156] = I_ERI_Hx4y_F3y_Px_S+ABY*I_ERI_Gx3y_F3y_Px_S;
  abcd[157] = I_ERI_Hx3yz_F3y_Px_S+ABY*I_ERI_Gx2yz_F3y_Px_S;
  abcd[158] = I_ERI_Hx2y2z_F3y_Px_S+ABY*I_ERI_Gxy2z_F3y_Px_S;
  abcd[159] = I_ERI_Hxy3z_F3y_Px_S+ABY*I_ERI_Gx3z_F3y_Px_S;
  abcd[160] = I_ERI_H5y_F3y_Px_S+ABY*I_ERI_G4y_F3y_Px_S;
  abcd[161] = I_ERI_H4yz_F3y_Px_S+ABY*I_ERI_G3yz_F3y_Px_S;
  abcd[162] = I_ERI_H3y2z_F3y_Px_S+ABY*I_ERI_G2y2z_F3y_Px_S;
  abcd[163] = I_ERI_H2y3z_F3y_Px_S+ABY*I_ERI_Gy3z_F3y_Px_S;
  abcd[164] = I_ERI_Hy4z_F3y_Px_S+ABY*I_ERI_G4z_F3y_Px_S;
  abcd[165] = I_ERI_H4xz_F3y_Px_S+ABZ*I_ERI_G4x_F3y_Px_S;
  abcd[166] = I_ERI_H3xyz_F3y_Px_S+ABZ*I_ERI_G3xy_F3y_Px_S;
  abcd[167] = I_ERI_H3x2z_F3y_Px_S+ABZ*I_ERI_G3xz_F3y_Px_S;
  abcd[168] = I_ERI_H2x2yz_F3y_Px_S+ABZ*I_ERI_G2x2y_F3y_Px_S;
  abcd[169] = I_ERI_H2xy2z_F3y_Px_S+ABZ*I_ERI_G2xyz_F3y_Px_S;
  abcd[170] = I_ERI_H2x3z_F3y_Px_S+ABZ*I_ERI_G2x2z_F3y_Px_S;
  abcd[171] = I_ERI_Hx3yz_F3y_Px_S+ABZ*I_ERI_Gx3y_F3y_Px_S;
  abcd[172] = I_ERI_Hx2y2z_F3y_Px_S+ABZ*I_ERI_Gx2yz_F3y_Px_S;
  abcd[173] = I_ERI_Hxy3z_F3y_Px_S+ABZ*I_ERI_Gxy2z_F3y_Px_S;
  abcd[174] = I_ERI_Hx4z_F3y_Px_S+ABZ*I_ERI_Gx3z_F3y_Px_S;
  abcd[175] = I_ERI_H4yz_F3y_Px_S+ABZ*I_ERI_G4y_F3y_Px_S;
  abcd[176] = I_ERI_H3y2z_F3y_Px_S+ABZ*I_ERI_G3yz_F3y_Px_S;
  abcd[177] = I_ERI_H2y3z_F3y_Px_S+ABZ*I_ERI_G2y2z_F3y_Px_S;
  abcd[178] = I_ERI_Hy4z_F3y_Px_S+ABZ*I_ERI_Gy3z_F3y_Px_S;
  abcd[179] = I_ERI_H5z_F3y_Px_S+ABZ*I_ERI_G4z_F3y_Px_S;
  abcd[180] = I_ERI_H4xz_F2yz_Px_S+ABZ*I_ERI_G4x_F2yz_Px_S;
  abcd[181] = I_ERI_H3xyz_F2yz_Px_S+ABZ*I_ERI_G3xy_F2yz_Px_S;
  abcd[182] = I_ERI_H3x2z_F2yz_Px_S+ABZ*I_ERI_G3xz_F2yz_Px_S;
  abcd[183] = I_ERI_H2x2yz_F2yz_Px_S+ABZ*I_ERI_G2x2y_F2yz_Px_S;
  abcd[184] = I_ERI_H2xy2z_F2yz_Px_S+ABZ*I_ERI_G2xyz_F2yz_Px_S;
  abcd[185] = I_ERI_H2x3z_F2yz_Px_S+ABZ*I_ERI_G2x2z_F2yz_Px_S;
  abcd[186] = I_ERI_Hx3yz_F2yz_Px_S+ABZ*I_ERI_Gx3y_F2yz_Px_S;
  abcd[187] = I_ERI_Hx2y2z_F2yz_Px_S+ABZ*I_ERI_Gx2yz_F2yz_Px_S;
  abcd[188] = I_ERI_Hxy3z_F2yz_Px_S+ABZ*I_ERI_Gxy2z_F2yz_Px_S;
  abcd[189] = I_ERI_Hx4z_F2yz_Px_S+ABZ*I_ERI_Gx3z_F2yz_Px_S;
  abcd[190] = I_ERI_H4yz_F2yz_Px_S+ABZ*I_ERI_G4y_F2yz_Px_S;
  abcd[191] = I_ERI_H3y2z_F2yz_Px_S+ABZ*I_ERI_G3yz_F2yz_Px_S;
  abcd[192] = I_ERI_H2y3z_F2yz_Px_S+ABZ*I_ERI_G2y2z_F2yz_Px_S;
  abcd[193] = I_ERI_Hy4z_F2yz_Px_S+ABZ*I_ERI_Gy3z_F2yz_Px_S;
  abcd[194] = I_ERI_H5z_F2yz_Px_S+ABZ*I_ERI_G4z_F2yz_Px_S;
  abcd[195] = I_ERI_H4xy_F3z_Px_S+ABY*I_ERI_G4x_F3z_Px_S;
  abcd[196] = I_ERI_H3x2y_F3z_Px_S+ABY*I_ERI_G3xy_F3z_Px_S;
  abcd[197] = I_ERI_H3xyz_F3z_Px_S+ABY*I_ERI_G3xz_F3z_Px_S;
  abcd[198] = I_ERI_H2x3y_F3z_Px_S+ABY*I_ERI_G2x2y_F3z_Px_S;
  abcd[199] = I_ERI_H2x2yz_F3z_Px_S+ABY*I_ERI_G2xyz_F3z_Px_S;
  abcd[200] = I_ERI_H2xy2z_F3z_Px_S+ABY*I_ERI_G2x2z_F3z_Px_S;
  abcd[201] = I_ERI_Hx4y_F3z_Px_S+ABY*I_ERI_Gx3y_F3z_Px_S;
  abcd[202] = I_ERI_Hx3yz_F3z_Px_S+ABY*I_ERI_Gx2yz_F3z_Px_S;
  abcd[203] = I_ERI_Hx2y2z_F3z_Px_S+ABY*I_ERI_Gxy2z_F3z_Px_S;
  abcd[204] = I_ERI_Hxy3z_F3z_Px_S+ABY*I_ERI_Gx3z_F3z_Px_S;
  abcd[205] = I_ERI_H5y_F3z_Px_S+ABY*I_ERI_G4y_F3z_Px_S;
  abcd[206] = I_ERI_H4yz_F3z_Px_S+ABY*I_ERI_G3yz_F3z_Px_S;
  abcd[207] = I_ERI_H3y2z_F3z_Px_S+ABY*I_ERI_G2y2z_F3z_Px_S;
  abcd[208] = I_ERI_H2y3z_F3z_Px_S+ABY*I_ERI_Gy3z_F3z_Px_S;
  abcd[209] = I_ERI_Hy4z_F3z_Px_S+ABY*I_ERI_G4z_F3z_Px_S;
  abcd[210] = I_ERI_H4xz_F3z_Px_S+ABZ*I_ERI_G4x_F3z_Px_S;
  abcd[211] = I_ERI_H3xyz_F3z_Px_S+ABZ*I_ERI_G3xy_F3z_Px_S;
  abcd[212] = I_ERI_H3x2z_F3z_Px_S+ABZ*I_ERI_G3xz_F3z_Px_S;
  abcd[213] = I_ERI_H2x2yz_F3z_Px_S+ABZ*I_ERI_G2x2y_F3z_Px_S;
  abcd[214] = I_ERI_H2xy2z_F3z_Px_S+ABZ*I_ERI_G2xyz_F3z_Px_S;
  abcd[215] = I_ERI_H2x3z_F3z_Px_S+ABZ*I_ERI_G2x2z_F3z_Px_S;
  abcd[216] = I_ERI_Hx3yz_F3z_Px_S+ABZ*I_ERI_Gx3y_F3z_Px_S;
  abcd[217] = I_ERI_Hx2y2z_F3z_Px_S+ABZ*I_ERI_Gx2yz_F3z_Px_S;
  abcd[218] = I_ERI_Hxy3z_F3z_Px_S+ABZ*I_ERI_Gxy2z_F3z_Px_S;
  abcd[219] = I_ERI_Hx4z_F3z_Px_S+ABZ*I_ERI_Gx3z_F3z_Px_S;
  abcd[220] = I_ERI_H4yz_F3z_Px_S+ABZ*I_ERI_G4y_F3z_Px_S;
  abcd[221] = I_ERI_H3y2z_F3z_Px_S+ABZ*I_ERI_G3yz_F3z_Px_S;
  abcd[222] = I_ERI_H2y3z_F3z_Px_S+ABZ*I_ERI_G2y2z_F3z_Px_S;
  abcd[223] = I_ERI_Hy4z_F3z_Px_S+ABZ*I_ERI_Gy3z_F3z_Px_S;
  abcd[224] = I_ERI_H5z_F3z_Px_S+ABZ*I_ERI_G4z_F3z_Px_S;
  abcd[225] = I_ERI_H5x_F3x_Py_S+ABX*I_ERI_G4x_F3x_Py_S;
  abcd[226] = I_ERI_H4xy_F3x_Py_S+ABX*I_ERI_G3xy_F3x_Py_S;
  abcd[227] = I_ERI_H4xz_F3x_Py_S+ABX*I_ERI_G3xz_F3x_Py_S;
  abcd[228] = I_ERI_H3x2y_F3x_Py_S+ABX*I_ERI_G2x2y_F3x_Py_S;
  abcd[229] = I_ERI_H3xyz_F3x_Py_S+ABX*I_ERI_G2xyz_F3x_Py_S;
  abcd[230] = I_ERI_H3x2z_F3x_Py_S+ABX*I_ERI_G2x2z_F3x_Py_S;
  abcd[231] = I_ERI_H2x3y_F3x_Py_S+ABX*I_ERI_Gx3y_F3x_Py_S;
  abcd[232] = I_ERI_H2x2yz_F3x_Py_S+ABX*I_ERI_Gx2yz_F3x_Py_S;
  abcd[233] = I_ERI_H2xy2z_F3x_Py_S+ABX*I_ERI_Gxy2z_F3x_Py_S;
  abcd[234] = I_ERI_H2x3z_F3x_Py_S+ABX*I_ERI_Gx3z_F3x_Py_S;
  abcd[235] = I_ERI_Hx4y_F3x_Py_S+ABX*I_ERI_G4y_F3x_Py_S;
  abcd[236] = I_ERI_Hx3yz_F3x_Py_S+ABX*I_ERI_G3yz_F3x_Py_S;
  abcd[237] = I_ERI_Hx2y2z_F3x_Py_S+ABX*I_ERI_G2y2z_F3x_Py_S;
  abcd[238] = I_ERI_Hxy3z_F3x_Py_S+ABX*I_ERI_Gy3z_F3x_Py_S;
  abcd[239] = I_ERI_Hx4z_F3x_Py_S+ABX*I_ERI_G4z_F3x_Py_S;
  abcd[240] = I_ERI_H4xy_F3x_Py_S+ABY*I_ERI_G4x_F3x_Py_S;
  abcd[241] = I_ERI_H3x2y_F3x_Py_S+ABY*I_ERI_G3xy_F3x_Py_S;
  abcd[242] = I_ERI_H3xyz_F3x_Py_S+ABY*I_ERI_G3xz_F3x_Py_S;
  abcd[243] = I_ERI_H2x3y_F3x_Py_S+ABY*I_ERI_G2x2y_F3x_Py_S;
  abcd[244] = I_ERI_H2x2yz_F3x_Py_S+ABY*I_ERI_G2xyz_F3x_Py_S;
  abcd[245] = I_ERI_H2xy2z_F3x_Py_S+ABY*I_ERI_G2x2z_F3x_Py_S;
  abcd[246] = I_ERI_Hx4y_F3x_Py_S+ABY*I_ERI_Gx3y_F3x_Py_S;
  abcd[247] = I_ERI_Hx3yz_F3x_Py_S+ABY*I_ERI_Gx2yz_F3x_Py_S;
  abcd[248] = I_ERI_Hx2y2z_F3x_Py_S+ABY*I_ERI_Gxy2z_F3x_Py_S;
  abcd[249] = I_ERI_Hxy3z_F3x_Py_S+ABY*I_ERI_Gx3z_F3x_Py_S;
  abcd[250] = I_ERI_H5y_F3x_Py_S+ABY*I_ERI_G4y_F3x_Py_S;
  abcd[251] = I_ERI_H4yz_F3x_Py_S+ABY*I_ERI_G3yz_F3x_Py_S;
  abcd[252] = I_ERI_H3y2z_F3x_Py_S+ABY*I_ERI_G2y2z_F3x_Py_S;
  abcd[253] = I_ERI_H2y3z_F3x_Py_S+ABY*I_ERI_Gy3z_F3x_Py_S;
  abcd[254] = I_ERI_Hy4z_F3x_Py_S+ABY*I_ERI_G4z_F3x_Py_S;
  abcd[255] = I_ERI_H4xz_F3x_Py_S+ABZ*I_ERI_G4x_F3x_Py_S;
  abcd[256] = I_ERI_H3xyz_F3x_Py_S+ABZ*I_ERI_G3xy_F3x_Py_S;
  abcd[257] = I_ERI_H3x2z_F3x_Py_S+ABZ*I_ERI_G3xz_F3x_Py_S;
  abcd[258] = I_ERI_H2x2yz_F3x_Py_S+ABZ*I_ERI_G2x2y_F3x_Py_S;
  abcd[259] = I_ERI_H2xy2z_F3x_Py_S+ABZ*I_ERI_G2xyz_F3x_Py_S;
  abcd[260] = I_ERI_H2x3z_F3x_Py_S+ABZ*I_ERI_G2x2z_F3x_Py_S;
  abcd[261] = I_ERI_Hx3yz_F3x_Py_S+ABZ*I_ERI_Gx3y_F3x_Py_S;
  abcd[262] = I_ERI_Hx2y2z_F3x_Py_S+ABZ*I_ERI_Gx2yz_F3x_Py_S;
  abcd[263] = I_ERI_Hxy3z_F3x_Py_S+ABZ*I_ERI_Gxy2z_F3x_Py_S;
  abcd[264] = I_ERI_Hx4z_F3x_Py_S+ABZ*I_ERI_Gx3z_F3x_Py_S;
  abcd[265] = I_ERI_H4yz_F3x_Py_S+ABZ*I_ERI_G4y_F3x_Py_S;
  abcd[266] = I_ERI_H3y2z_F3x_Py_S+ABZ*I_ERI_G3yz_F3x_Py_S;
  abcd[267] = I_ERI_H2y3z_F3x_Py_S+ABZ*I_ERI_G2y2z_F3x_Py_S;
  abcd[268] = I_ERI_Hy4z_F3x_Py_S+ABZ*I_ERI_Gy3z_F3x_Py_S;
  abcd[269] = I_ERI_H5z_F3x_Py_S+ABZ*I_ERI_G4z_F3x_Py_S;
  abcd[270] = I_ERI_H4xy_F2xy_Py_S+ABY*I_ERI_G4x_F2xy_Py_S;
  abcd[271] = I_ERI_H3x2y_F2xy_Py_S+ABY*I_ERI_G3xy_F2xy_Py_S;
  abcd[272] = I_ERI_H3xyz_F2xy_Py_S+ABY*I_ERI_G3xz_F2xy_Py_S;
  abcd[273] = I_ERI_H2x3y_F2xy_Py_S+ABY*I_ERI_G2x2y_F2xy_Py_S;
  abcd[274] = I_ERI_H2x2yz_F2xy_Py_S+ABY*I_ERI_G2xyz_F2xy_Py_S;
  abcd[275] = I_ERI_H2xy2z_F2xy_Py_S+ABY*I_ERI_G2x2z_F2xy_Py_S;
  abcd[276] = I_ERI_Hx4y_F2xy_Py_S+ABY*I_ERI_Gx3y_F2xy_Py_S;
  abcd[277] = I_ERI_Hx3yz_F2xy_Py_S+ABY*I_ERI_Gx2yz_F2xy_Py_S;
  abcd[278] = I_ERI_Hx2y2z_F2xy_Py_S+ABY*I_ERI_Gxy2z_F2xy_Py_S;
  abcd[279] = I_ERI_Hxy3z_F2xy_Py_S+ABY*I_ERI_Gx3z_F2xy_Py_S;
  abcd[280] = I_ERI_H5y_F2xy_Py_S+ABY*I_ERI_G4y_F2xy_Py_S;
  abcd[281] = I_ERI_H4yz_F2xy_Py_S+ABY*I_ERI_G3yz_F2xy_Py_S;
  abcd[282] = I_ERI_H3y2z_F2xy_Py_S+ABY*I_ERI_G2y2z_F2xy_Py_S;
  abcd[283] = I_ERI_H2y3z_F2xy_Py_S+ABY*I_ERI_Gy3z_F2xy_Py_S;
  abcd[284] = I_ERI_Hy4z_F2xy_Py_S+ABY*I_ERI_G4z_F2xy_Py_S;
  abcd[285] = I_ERI_H4xz_F2xy_Py_S+ABZ*I_ERI_G4x_F2xy_Py_S;
  abcd[286] = I_ERI_H3xyz_F2xy_Py_S+ABZ*I_ERI_G3xy_F2xy_Py_S;
  abcd[287] = I_ERI_H3x2z_F2xy_Py_S+ABZ*I_ERI_G3xz_F2xy_Py_S;
  abcd[288] = I_ERI_H2x2yz_F2xy_Py_S+ABZ*I_ERI_G2x2y_F2xy_Py_S;
  abcd[289] = I_ERI_H2xy2z_F2xy_Py_S+ABZ*I_ERI_G2xyz_F2xy_Py_S;
  abcd[290] = I_ERI_H2x3z_F2xy_Py_S+ABZ*I_ERI_G2x2z_F2xy_Py_S;
  abcd[291] = I_ERI_Hx3yz_F2xy_Py_S+ABZ*I_ERI_Gx3y_F2xy_Py_S;
  abcd[292] = I_ERI_Hx2y2z_F2xy_Py_S+ABZ*I_ERI_Gx2yz_F2xy_Py_S;
  abcd[293] = I_ERI_Hxy3z_F2xy_Py_S+ABZ*I_ERI_Gxy2z_F2xy_Py_S;
  abcd[294] = I_ERI_Hx4z_F2xy_Py_S+ABZ*I_ERI_Gx3z_F2xy_Py_S;
  abcd[295] = I_ERI_H4yz_F2xy_Py_S+ABZ*I_ERI_G4y_F2xy_Py_S;
  abcd[296] = I_ERI_H3y2z_F2xy_Py_S+ABZ*I_ERI_G3yz_F2xy_Py_S;
  abcd[297] = I_ERI_H2y3z_F2xy_Py_S+ABZ*I_ERI_G2y2z_F2xy_Py_S;
  abcd[298] = I_ERI_Hy4z_F2xy_Py_S+ABZ*I_ERI_Gy3z_F2xy_Py_S;
  abcd[299] = I_ERI_H5z_F2xy_Py_S+ABZ*I_ERI_G4z_F2xy_Py_S;
  abcd[300] = I_ERI_H4xz_F2xz_Py_S+ABZ*I_ERI_G4x_F2xz_Py_S;
  abcd[301] = I_ERI_H3xyz_F2xz_Py_S+ABZ*I_ERI_G3xy_F2xz_Py_S;
  abcd[302] = I_ERI_H3x2z_F2xz_Py_S+ABZ*I_ERI_G3xz_F2xz_Py_S;
  abcd[303] = I_ERI_H2x2yz_F2xz_Py_S+ABZ*I_ERI_G2x2y_F2xz_Py_S;
  abcd[304] = I_ERI_H2xy2z_F2xz_Py_S+ABZ*I_ERI_G2xyz_F2xz_Py_S;
  abcd[305] = I_ERI_H2x3z_F2xz_Py_S+ABZ*I_ERI_G2x2z_F2xz_Py_S;
  abcd[306] = I_ERI_Hx3yz_F2xz_Py_S+ABZ*I_ERI_Gx3y_F2xz_Py_S;
  abcd[307] = I_ERI_Hx2y2z_F2xz_Py_S+ABZ*I_ERI_Gx2yz_F2xz_Py_S;
  abcd[308] = I_ERI_Hxy3z_F2xz_Py_S+ABZ*I_ERI_Gxy2z_F2xz_Py_S;
  abcd[309] = I_ERI_Hx4z_F2xz_Py_S+ABZ*I_ERI_Gx3z_F2xz_Py_S;
  abcd[310] = I_ERI_H4yz_F2xz_Py_S+ABZ*I_ERI_G4y_F2xz_Py_S;
  abcd[311] = I_ERI_H3y2z_F2xz_Py_S+ABZ*I_ERI_G3yz_F2xz_Py_S;
  abcd[312] = I_ERI_H2y3z_F2xz_Py_S+ABZ*I_ERI_G2y2z_F2xz_Py_S;
  abcd[313] = I_ERI_Hy4z_F2xz_Py_S+ABZ*I_ERI_Gy3z_F2xz_Py_S;
  abcd[314] = I_ERI_H5z_F2xz_Py_S+ABZ*I_ERI_G4z_F2xz_Py_S;
  abcd[315] = I_ERI_H5x_F3y_Py_S+ABX*I_ERI_G4x_F3y_Py_S;
  abcd[316] = I_ERI_H4xy_F3y_Py_S+ABX*I_ERI_G3xy_F3y_Py_S;
  abcd[317] = I_ERI_H4xz_F3y_Py_S+ABX*I_ERI_G3xz_F3y_Py_S;
  abcd[318] = I_ERI_H3x2y_F3y_Py_S+ABX*I_ERI_G2x2y_F3y_Py_S;
  abcd[319] = I_ERI_H3xyz_F3y_Py_S+ABX*I_ERI_G2xyz_F3y_Py_S;
  abcd[320] = I_ERI_H3x2z_F3y_Py_S+ABX*I_ERI_G2x2z_F3y_Py_S;
  abcd[321] = I_ERI_H2x3y_F3y_Py_S+ABX*I_ERI_Gx3y_F3y_Py_S;
  abcd[322] = I_ERI_H2x2yz_F3y_Py_S+ABX*I_ERI_Gx2yz_F3y_Py_S;
  abcd[323] = I_ERI_H2xy2z_F3y_Py_S+ABX*I_ERI_Gxy2z_F3y_Py_S;
  abcd[324] = I_ERI_H2x3z_F3y_Py_S+ABX*I_ERI_Gx3z_F3y_Py_S;
  abcd[325] = I_ERI_Hx4y_F3y_Py_S+ABX*I_ERI_G4y_F3y_Py_S;
  abcd[326] = I_ERI_Hx3yz_F3y_Py_S+ABX*I_ERI_G3yz_F3y_Py_S;
  abcd[327] = I_ERI_Hx2y2z_F3y_Py_S+ABX*I_ERI_G2y2z_F3y_Py_S;
  abcd[328] = I_ERI_Hxy3z_F3y_Py_S+ABX*I_ERI_Gy3z_F3y_Py_S;
  abcd[329] = I_ERI_Hx4z_F3y_Py_S+ABX*I_ERI_G4z_F3y_Py_S;
  abcd[330] = I_ERI_H4xz_Fx2y_Py_S+ABZ*I_ERI_G4x_Fx2y_Py_S;
  abcd[331] = I_ERI_H3xyz_Fx2y_Py_S+ABZ*I_ERI_G3xy_Fx2y_Py_S;
  abcd[332] = I_ERI_H3x2z_Fx2y_Py_S+ABZ*I_ERI_G3xz_Fx2y_Py_S;
  abcd[333] = I_ERI_H2x2yz_Fx2y_Py_S+ABZ*I_ERI_G2x2y_Fx2y_Py_S;
  abcd[334] = I_ERI_H2xy2z_Fx2y_Py_S+ABZ*I_ERI_G2xyz_Fx2y_Py_S;
  abcd[335] = I_ERI_H2x3z_Fx2y_Py_S+ABZ*I_ERI_G2x2z_Fx2y_Py_S;
  abcd[336] = I_ERI_Hx3yz_Fx2y_Py_S+ABZ*I_ERI_Gx3y_Fx2y_Py_S;
  abcd[337] = I_ERI_Hx2y2z_Fx2y_Py_S+ABZ*I_ERI_Gx2yz_Fx2y_Py_S;
  abcd[338] = I_ERI_Hxy3z_Fx2y_Py_S+ABZ*I_ERI_Gxy2z_Fx2y_Py_S;
  abcd[339] = I_ERI_Hx4z_Fx2y_Py_S+ABZ*I_ERI_Gx3z_Fx2y_Py_S;
  abcd[340] = I_ERI_H4yz_Fx2y_Py_S+ABZ*I_ERI_G4y_Fx2y_Py_S;
  abcd[341] = I_ERI_H3y2z_Fx2y_Py_S+ABZ*I_ERI_G3yz_Fx2y_Py_S;
  abcd[342] = I_ERI_H2y3z_Fx2y_Py_S+ABZ*I_ERI_G2y2z_Fx2y_Py_S;
  abcd[343] = I_ERI_Hy4z_Fx2y_Py_S+ABZ*I_ERI_Gy3z_Fx2y_Py_S;
  abcd[344] = I_ERI_H5z_Fx2y_Py_S+ABZ*I_ERI_G4z_Fx2y_Py_S;
  abcd[345] = I_ERI_H4xy_Fx2z_Py_S+ABY*I_ERI_G4x_Fx2z_Py_S;
  abcd[346] = I_ERI_H3x2y_Fx2z_Py_S+ABY*I_ERI_G3xy_Fx2z_Py_S;
  abcd[347] = I_ERI_H3xyz_Fx2z_Py_S+ABY*I_ERI_G3xz_Fx2z_Py_S;
  abcd[348] = I_ERI_H2x3y_Fx2z_Py_S+ABY*I_ERI_G2x2y_Fx2z_Py_S;
  abcd[349] = I_ERI_H2x2yz_Fx2z_Py_S+ABY*I_ERI_G2xyz_Fx2z_Py_S;
  abcd[350] = I_ERI_H2xy2z_Fx2z_Py_S+ABY*I_ERI_G2x2z_Fx2z_Py_S;
  abcd[351] = I_ERI_Hx4y_Fx2z_Py_S+ABY*I_ERI_Gx3y_Fx2z_Py_S;
  abcd[352] = I_ERI_Hx3yz_Fx2z_Py_S+ABY*I_ERI_Gx2yz_Fx2z_Py_S;
  abcd[353] = I_ERI_Hx2y2z_Fx2z_Py_S+ABY*I_ERI_Gxy2z_Fx2z_Py_S;
  abcd[354] = I_ERI_Hxy3z_Fx2z_Py_S+ABY*I_ERI_Gx3z_Fx2z_Py_S;
  abcd[355] = I_ERI_H5y_Fx2z_Py_S+ABY*I_ERI_G4y_Fx2z_Py_S;
  abcd[356] = I_ERI_H4yz_Fx2z_Py_S+ABY*I_ERI_G3yz_Fx2z_Py_S;
  abcd[357] = I_ERI_H3y2z_Fx2z_Py_S+ABY*I_ERI_G2y2z_Fx2z_Py_S;
  abcd[358] = I_ERI_H2y3z_Fx2z_Py_S+ABY*I_ERI_Gy3z_Fx2z_Py_S;
  abcd[359] = I_ERI_Hy4z_Fx2z_Py_S+ABY*I_ERI_G4z_Fx2z_Py_S;
  abcd[360] = I_ERI_H5x_F3z_Py_S+ABX*I_ERI_G4x_F3z_Py_S;
  abcd[361] = I_ERI_H4xy_F3z_Py_S+ABX*I_ERI_G3xy_F3z_Py_S;
  abcd[362] = I_ERI_H4xz_F3z_Py_S+ABX*I_ERI_G3xz_F3z_Py_S;
  abcd[363] = I_ERI_H3x2y_F3z_Py_S+ABX*I_ERI_G2x2y_F3z_Py_S;
  abcd[364] = I_ERI_H3xyz_F3z_Py_S+ABX*I_ERI_G2xyz_F3z_Py_S;
  abcd[365] = I_ERI_H3x2z_F3z_Py_S+ABX*I_ERI_G2x2z_F3z_Py_S;
  abcd[366] = I_ERI_H2x3y_F3z_Py_S+ABX*I_ERI_Gx3y_F3z_Py_S;
  abcd[367] = I_ERI_H2x2yz_F3z_Py_S+ABX*I_ERI_Gx2yz_F3z_Py_S;
  abcd[368] = I_ERI_H2xy2z_F3z_Py_S+ABX*I_ERI_Gxy2z_F3z_Py_S;
  abcd[369] = I_ERI_H2x3z_F3z_Py_S+ABX*I_ERI_Gx3z_F3z_Py_S;
  abcd[370] = I_ERI_Hx4y_F3z_Py_S+ABX*I_ERI_G4y_F3z_Py_S;
  abcd[371] = I_ERI_Hx3yz_F3z_Py_S+ABX*I_ERI_G3yz_F3z_Py_S;
  abcd[372] = I_ERI_Hx2y2z_F3z_Py_S+ABX*I_ERI_G2y2z_F3z_Py_S;
  abcd[373] = I_ERI_Hxy3z_F3z_Py_S+ABX*I_ERI_Gy3z_F3z_Py_S;
  abcd[374] = I_ERI_Hx4z_F3z_Py_S+ABX*I_ERI_G4z_F3z_Py_S;
  abcd[375] = I_ERI_H4xy_F3y_Py_S+ABY*I_ERI_G4x_F3y_Py_S;
  abcd[376] = I_ERI_H3x2y_F3y_Py_S+ABY*I_ERI_G3xy_F3y_Py_S;
  abcd[377] = I_ERI_H3xyz_F3y_Py_S+ABY*I_ERI_G3xz_F3y_Py_S;
  abcd[378] = I_ERI_H2x3y_F3y_Py_S+ABY*I_ERI_G2x2y_F3y_Py_S;
  abcd[379] = I_ERI_H2x2yz_F3y_Py_S+ABY*I_ERI_G2xyz_F3y_Py_S;
  abcd[380] = I_ERI_H2xy2z_F3y_Py_S+ABY*I_ERI_G2x2z_F3y_Py_S;
  abcd[381] = I_ERI_Hx4y_F3y_Py_S+ABY*I_ERI_Gx3y_F3y_Py_S;
  abcd[382] = I_ERI_Hx3yz_F3y_Py_S+ABY*I_ERI_Gx2yz_F3y_Py_S;
  abcd[383] = I_ERI_Hx2y2z_F3y_Py_S+ABY*I_ERI_Gxy2z_F3y_Py_S;
  abcd[384] = I_ERI_Hxy3z_F3y_Py_S+ABY*I_ERI_Gx3z_F3y_Py_S;
  abcd[385] = I_ERI_H5y_F3y_Py_S+ABY*I_ERI_G4y_F3y_Py_S;
  abcd[386] = I_ERI_H4yz_F3y_Py_S+ABY*I_ERI_G3yz_F3y_Py_S;
  abcd[387] = I_ERI_H3y2z_F3y_Py_S+ABY*I_ERI_G2y2z_F3y_Py_S;
  abcd[388] = I_ERI_H2y3z_F3y_Py_S+ABY*I_ERI_Gy3z_F3y_Py_S;
  abcd[389] = I_ERI_Hy4z_F3y_Py_S+ABY*I_ERI_G4z_F3y_Py_S;
  abcd[390] = I_ERI_H4xz_F3y_Py_S+ABZ*I_ERI_G4x_F3y_Py_S;
  abcd[391] = I_ERI_H3xyz_F3y_Py_S+ABZ*I_ERI_G3xy_F3y_Py_S;
  abcd[392] = I_ERI_H3x2z_F3y_Py_S+ABZ*I_ERI_G3xz_F3y_Py_S;
  abcd[393] = I_ERI_H2x2yz_F3y_Py_S+ABZ*I_ERI_G2x2y_F3y_Py_S;
  abcd[394] = I_ERI_H2xy2z_F3y_Py_S+ABZ*I_ERI_G2xyz_F3y_Py_S;
  abcd[395] = I_ERI_H2x3z_F3y_Py_S+ABZ*I_ERI_G2x2z_F3y_Py_S;
  abcd[396] = I_ERI_Hx3yz_F3y_Py_S+ABZ*I_ERI_Gx3y_F3y_Py_S;
  abcd[397] = I_ERI_Hx2y2z_F3y_Py_S+ABZ*I_ERI_Gx2yz_F3y_Py_S;
  abcd[398] = I_ERI_Hxy3z_F3y_Py_S+ABZ*I_ERI_Gxy2z_F3y_Py_S;
  abcd[399] = I_ERI_Hx4z_F3y_Py_S+ABZ*I_ERI_Gx3z_F3y_Py_S;
  abcd[400] = I_ERI_H4yz_F3y_Py_S+ABZ*I_ERI_G4y_F3y_Py_S;
  abcd[401] = I_ERI_H3y2z_F3y_Py_S+ABZ*I_ERI_G3yz_F3y_Py_S;
  abcd[402] = I_ERI_H2y3z_F3y_Py_S+ABZ*I_ERI_G2y2z_F3y_Py_S;
  abcd[403] = I_ERI_Hy4z_F3y_Py_S+ABZ*I_ERI_Gy3z_F3y_Py_S;
  abcd[404] = I_ERI_H5z_F3y_Py_S+ABZ*I_ERI_G4z_F3y_Py_S;
  abcd[405] = I_ERI_H4xz_F2yz_Py_S+ABZ*I_ERI_G4x_F2yz_Py_S;
  abcd[406] = I_ERI_H3xyz_F2yz_Py_S+ABZ*I_ERI_G3xy_F2yz_Py_S;
  abcd[407] = I_ERI_H3x2z_F2yz_Py_S+ABZ*I_ERI_G3xz_F2yz_Py_S;
  abcd[408] = I_ERI_H2x2yz_F2yz_Py_S+ABZ*I_ERI_G2x2y_F2yz_Py_S;
  abcd[409] = I_ERI_H2xy2z_F2yz_Py_S+ABZ*I_ERI_G2xyz_F2yz_Py_S;
  abcd[410] = I_ERI_H2x3z_F2yz_Py_S+ABZ*I_ERI_G2x2z_F2yz_Py_S;
  abcd[411] = I_ERI_Hx3yz_F2yz_Py_S+ABZ*I_ERI_Gx3y_F2yz_Py_S;
  abcd[412] = I_ERI_Hx2y2z_F2yz_Py_S+ABZ*I_ERI_Gx2yz_F2yz_Py_S;
  abcd[413] = I_ERI_Hxy3z_F2yz_Py_S+ABZ*I_ERI_Gxy2z_F2yz_Py_S;
  abcd[414] = I_ERI_Hx4z_F2yz_Py_S+ABZ*I_ERI_Gx3z_F2yz_Py_S;
  abcd[415] = I_ERI_H4yz_F2yz_Py_S+ABZ*I_ERI_G4y_F2yz_Py_S;
  abcd[416] = I_ERI_H3y2z_F2yz_Py_S+ABZ*I_ERI_G3yz_F2yz_Py_S;
  abcd[417] = I_ERI_H2y3z_F2yz_Py_S+ABZ*I_ERI_G2y2z_F2yz_Py_S;
  abcd[418] = I_ERI_Hy4z_F2yz_Py_S+ABZ*I_ERI_Gy3z_F2yz_Py_S;
  abcd[419] = I_ERI_H5z_F2yz_Py_S+ABZ*I_ERI_G4z_F2yz_Py_S;
  abcd[420] = I_ERI_H4xy_F3z_Py_S+ABY*I_ERI_G4x_F3z_Py_S;
  abcd[421] = I_ERI_H3x2y_F3z_Py_S+ABY*I_ERI_G3xy_F3z_Py_S;
  abcd[422] = I_ERI_H3xyz_F3z_Py_S+ABY*I_ERI_G3xz_F3z_Py_S;
  abcd[423] = I_ERI_H2x3y_F3z_Py_S+ABY*I_ERI_G2x2y_F3z_Py_S;
  abcd[424] = I_ERI_H2x2yz_F3z_Py_S+ABY*I_ERI_G2xyz_F3z_Py_S;
  abcd[425] = I_ERI_H2xy2z_F3z_Py_S+ABY*I_ERI_G2x2z_F3z_Py_S;
  abcd[426] = I_ERI_Hx4y_F3z_Py_S+ABY*I_ERI_Gx3y_F3z_Py_S;
  abcd[427] = I_ERI_Hx3yz_F3z_Py_S+ABY*I_ERI_Gx2yz_F3z_Py_S;
  abcd[428] = I_ERI_Hx2y2z_F3z_Py_S+ABY*I_ERI_Gxy2z_F3z_Py_S;
  abcd[429] = I_ERI_Hxy3z_F3z_Py_S+ABY*I_ERI_Gx3z_F3z_Py_S;
  abcd[430] = I_ERI_H5y_F3z_Py_S+ABY*I_ERI_G4y_F3z_Py_S;
  abcd[431] = I_ERI_H4yz_F3z_Py_S+ABY*I_ERI_G3yz_F3z_Py_S;
  abcd[432] = I_ERI_H3y2z_F3z_Py_S+ABY*I_ERI_G2y2z_F3z_Py_S;
  abcd[433] = I_ERI_H2y3z_F3z_Py_S+ABY*I_ERI_Gy3z_F3z_Py_S;
  abcd[434] = I_ERI_Hy4z_F3z_Py_S+ABY*I_ERI_G4z_F3z_Py_S;
  abcd[435] = I_ERI_H4xz_F3z_Py_S+ABZ*I_ERI_G4x_F3z_Py_S;
  abcd[436] = I_ERI_H3xyz_F3z_Py_S+ABZ*I_ERI_G3xy_F3z_Py_S;
  abcd[437] = I_ERI_H3x2z_F3z_Py_S+ABZ*I_ERI_G3xz_F3z_Py_S;
  abcd[438] = I_ERI_H2x2yz_F3z_Py_S+ABZ*I_ERI_G2x2y_F3z_Py_S;
  abcd[439] = I_ERI_H2xy2z_F3z_Py_S+ABZ*I_ERI_G2xyz_F3z_Py_S;
  abcd[440] = I_ERI_H2x3z_F3z_Py_S+ABZ*I_ERI_G2x2z_F3z_Py_S;
  abcd[441] = I_ERI_Hx3yz_F3z_Py_S+ABZ*I_ERI_Gx3y_F3z_Py_S;
  abcd[442] = I_ERI_Hx2y2z_F3z_Py_S+ABZ*I_ERI_Gx2yz_F3z_Py_S;
  abcd[443] = I_ERI_Hxy3z_F3z_Py_S+ABZ*I_ERI_Gxy2z_F3z_Py_S;
  abcd[444] = I_ERI_Hx4z_F3z_Py_S+ABZ*I_ERI_Gx3z_F3z_Py_S;
  abcd[445] = I_ERI_H4yz_F3z_Py_S+ABZ*I_ERI_G4y_F3z_Py_S;
  abcd[446] = I_ERI_H3y2z_F3z_Py_S+ABZ*I_ERI_G3yz_F3z_Py_S;
  abcd[447] = I_ERI_H2y3z_F3z_Py_S+ABZ*I_ERI_G2y2z_F3z_Py_S;
  abcd[448] = I_ERI_Hy4z_F3z_Py_S+ABZ*I_ERI_Gy3z_F3z_Py_S;
  abcd[449] = I_ERI_H5z_F3z_Py_S+ABZ*I_ERI_G4z_F3z_Py_S;
  abcd[450] = I_ERI_H5x_F3x_Pz_S+ABX*I_ERI_G4x_F3x_Pz_S;
  abcd[451] = I_ERI_H4xy_F3x_Pz_S+ABX*I_ERI_G3xy_F3x_Pz_S;
  abcd[452] = I_ERI_H4xz_F3x_Pz_S+ABX*I_ERI_G3xz_F3x_Pz_S;
  abcd[453] = I_ERI_H3x2y_F3x_Pz_S+ABX*I_ERI_G2x2y_F3x_Pz_S;
  abcd[454] = I_ERI_H3xyz_F3x_Pz_S+ABX*I_ERI_G2xyz_F3x_Pz_S;
  abcd[455] = I_ERI_H3x2z_F3x_Pz_S+ABX*I_ERI_G2x2z_F3x_Pz_S;
  abcd[456] = I_ERI_H2x3y_F3x_Pz_S+ABX*I_ERI_Gx3y_F3x_Pz_S;
  abcd[457] = I_ERI_H2x2yz_F3x_Pz_S+ABX*I_ERI_Gx2yz_F3x_Pz_S;
  abcd[458] = I_ERI_H2xy2z_F3x_Pz_S+ABX*I_ERI_Gxy2z_F3x_Pz_S;
  abcd[459] = I_ERI_H2x3z_F3x_Pz_S+ABX*I_ERI_Gx3z_F3x_Pz_S;
  abcd[460] = I_ERI_Hx4y_F3x_Pz_S+ABX*I_ERI_G4y_F3x_Pz_S;
  abcd[461] = I_ERI_Hx3yz_F3x_Pz_S+ABX*I_ERI_G3yz_F3x_Pz_S;
  abcd[462] = I_ERI_Hx2y2z_F3x_Pz_S+ABX*I_ERI_G2y2z_F3x_Pz_S;
  abcd[463] = I_ERI_Hxy3z_F3x_Pz_S+ABX*I_ERI_Gy3z_F3x_Pz_S;
  abcd[464] = I_ERI_Hx4z_F3x_Pz_S+ABX*I_ERI_G4z_F3x_Pz_S;
  abcd[465] = I_ERI_H4xy_F3x_Pz_S+ABY*I_ERI_G4x_F3x_Pz_S;
  abcd[466] = I_ERI_H3x2y_F3x_Pz_S+ABY*I_ERI_G3xy_F3x_Pz_S;
  abcd[467] = I_ERI_H3xyz_F3x_Pz_S+ABY*I_ERI_G3xz_F3x_Pz_S;
  abcd[468] = I_ERI_H2x3y_F3x_Pz_S+ABY*I_ERI_G2x2y_F3x_Pz_S;
  abcd[469] = I_ERI_H2x2yz_F3x_Pz_S+ABY*I_ERI_G2xyz_F3x_Pz_S;
  abcd[470] = I_ERI_H2xy2z_F3x_Pz_S+ABY*I_ERI_G2x2z_F3x_Pz_S;
  abcd[471] = I_ERI_Hx4y_F3x_Pz_S+ABY*I_ERI_Gx3y_F3x_Pz_S;
  abcd[472] = I_ERI_Hx3yz_F3x_Pz_S+ABY*I_ERI_Gx2yz_F3x_Pz_S;
  abcd[473] = I_ERI_Hx2y2z_F3x_Pz_S+ABY*I_ERI_Gxy2z_F3x_Pz_S;
  abcd[474] = I_ERI_Hxy3z_F3x_Pz_S+ABY*I_ERI_Gx3z_F3x_Pz_S;
  abcd[475] = I_ERI_H5y_F3x_Pz_S+ABY*I_ERI_G4y_F3x_Pz_S;
  abcd[476] = I_ERI_H4yz_F3x_Pz_S+ABY*I_ERI_G3yz_F3x_Pz_S;
  abcd[477] = I_ERI_H3y2z_F3x_Pz_S+ABY*I_ERI_G2y2z_F3x_Pz_S;
  abcd[478] = I_ERI_H2y3z_F3x_Pz_S+ABY*I_ERI_Gy3z_F3x_Pz_S;
  abcd[479] = I_ERI_Hy4z_F3x_Pz_S+ABY*I_ERI_G4z_F3x_Pz_S;
  abcd[480] = I_ERI_H4xz_F3x_Pz_S+ABZ*I_ERI_G4x_F3x_Pz_S;
  abcd[481] = I_ERI_H3xyz_F3x_Pz_S+ABZ*I_ERI_G3xy_F3x_Pz_S;
  abcd[482] = I_ERI_H3x2z_F3x_Pz_S+ABZ*I_ERI_G3xz_F3x_Pz_S;
  abcd[483] = I_ERI_H2x2yz_F3x_Pz_S+ABZ*I_ERI_G2x2y_F3x_Pz_S;
  abcd[484] = I_ERI_H2xy2z_F3x_Pz_S+ABZ*I_ERI_G2xyz_F3x_Pz_S;
  abcd[485] = I_ERI_H2x3z_F3x_Pz_S+ABZ*I_ERI_G2x2z_F3x_Pz_S;
  abcd[486] = I_ERI_Hx3yz_F3x_Pz_S+ABZ*I_ERI_Gx3y_F3x_Pz_S;
  abcd[487] = I_ERI_Hx2y2z_F3x_Pz_S+ABZ*I_ERI_Gx2yz_F3x_Pz_S;
  abcd[488] = I_ERI_Hxy3z_F3x_Pz_S+ABZ*I_ERI_Gxy2z_F3x_Pz_S;
  abcd[489] = I_ERI_Hx4z_F3x_Pz_S+ABZ*I_ERI_Gx3z_F3x_Pz_S;
  abcd[490] = I_ERI_H4yz_F3x_Pz_S+ABZ*I_ERI_G4y_F3x_Pz_S;
  abcd[491] = I_ERI_H3y2z_F3x_Pz_S+ABZ*I_ERI_G3yz_F3x_Pz_S;
  abcd[492] = I_ERI_H2y3z_F3x_Pz_S+ABZ*I_ERI_G2y2z_F3x_Pz_S;
  abcd[493] = I_ERI_Hy4z_F3x_Pz_S+ABZ*I_ERI_Gy3z_F3x_Pz_S;
  abcd[494] = I_ERI_H5z_F3x_Pz_S+ABZ*I_ERI_G4z_F3x_Pz_S;
  abcd[495] = I_ERI_H4xy_F2xy_Pz_S+ABY*I_ERI_G4x_F2xy_Pz_S;
  abcd[496] = I_ERI_H3x2y_F2xy_Pz_S+ABY*I_ERI_G3xy_F2xy_Pz_S;
  abcd[497] = I_ERI_H3xyz_F2xy_Pz_S+ABY*I_ERI_G3xz_F2xy_Pz_S;
  abcd[498] = I_ERI_H2x3y_F2xy_Pz_S+ABY*I_ERI_G2x2y_F2xy_Pz_S;
  abcd[499] = I_ERI_H2x2yz_F2xy_Pz_S+ABY*I_ERI_G2xyz_F2xy_Pz_S;
  abcd[500] = I_ERI_H2xy2z_F2xy_Pz_S+ABY*I_ERI_G2x2z_F2xy_Pz_S;
  abcd[501] = I_ERI_Hx4y_F2xy_Pz_S+ABY*I_ERI_Gx3y_F2xy_Pz_S;
  abcd[502] = I_ERI_Hx3yz_F2xy_Pz_S+ABY*I_ERI_Gx2yz_F2xy_Pz_S;
  abcd[503] = I_ERI_Hx2y2z_F2xy_Pz_S+ABY*I_ERI_Gxy2z_F2xy_Pz_S;
  abcd[504] = I_ERI_Hxy3z_F2xy_Pz_S+ABY*I_ERI_Gx3z_F2xy_Pz_S;
  abcd[505] = I_ERI_H5y_F2xy_Pz_S+ABY*I_ERI_G4y_F2xy_Pz_S;
  abcd[506] = I_ERI_H4yz_F2xy_Pz_S+ABY*I_ERI_G3yz_F2xy_Pz_S;
  abcd[507] = I_ERI_H3y2z_F2xy_Pz_S+ABY*I_ERI_G2y2z_F2xy_Pz_S;
  abcd[508] = I_ERI_H2y3z_F2xy_Pz_S+ABY*I_ERI_Gy3z_F2xy_Pz_S;
  abcd[509] = I_ERI_Hy4z_F2xy_Pz_S+ABY*I_ERI_G4z_F2xy_Pz_S;
  abcd[510] = I_ERI_H4xz_F2xy_Pz_S+ABZ*I_ERI_G4x_F2xy_Pz_S;
  abcd[511] = I_ERI_H3xyz_F2xy_Pz_S+ABZ*I_ERI_G3xy_F2xy_Pz_S;
  abcd[512] = I_ERI_H3x2z_F2xy_Pz_S+ABZ*I_ERI_G3xz_F2xy_Pz_S;
  abcd[513] = I_ERI_H2x2yz_F2xy_Pz_S+ABZ*I_ERI_G2x2y_F2xy_Pz_S;
  abcd[514] = I_ERI_H2xy2z_F2xy_Pz_S+ABZ*I_ERI_G2xyz_F2xy_Pz_S;
  abcd[515] = I_ERI_H2x3z_F2xy_Pz_S+ABZ*I_ERI_G2x2z_F2xy_Pz_S;
  abcd[516] = I_ERI_Hx3yz_F2xy_Pz_S+ABZ*I_ERI_Gx3y_F2xy_Pz_S;
  abcd[517] = I_ERI_Hx2y2z_F2xy_Pz_S+ABZ*I_ERI_Gx2yz_F2xy_Pz_S;
  abcd[518] = I_ERI_Hxy3z_F2xy_Pz_S+ABZ*I_ERI_Gxy2z_F2xy_Pz_S;
  abcd[519] = I_ERI_Hx4z_F2xy_Pz_S+ABZ*I_ERI_Gx3z_F2xy_Pz_S;
  abcd[520] = I_ERI_H4yz_F2xy_Pz_S+ABZ*I_ERI_G4y_F2xy_Pz_S;
  abcd[521] = I_ERI_H3y2z_F2xy_Pz_S+ABZ*I_ERI_G3yz_F2xy_Pz_S;
  abcd[522] = I_ERI_H2y3z_F2xy_Pz_S+ABZ*I_ERI_G2y2z_F2xy_Pz_S;
  abcd[523] = I_ERI_Hy4z_F2xy_Pz_S+ABZ*I_ERI_Gy3z_F2xy_Pz_S;
  abcd[524] = I_ERI_H5z_F2xy_Pz_S+ABZ*I_ERI_G4z_F2xy_Pz_S;
  abcd[525] = I_ERI_H4xz_F2xz_Pz_S+ABZ*I_ERI_G4x_F2xz_Pz_S;
  abcd[526] = I_ERI_H3xyz_F2xz_Pz_S+ABZ*I_ERI_G3xy_F2xz_Pz_S;
  abcd[527] = I_ERI_H3x2z_F2xz_Pz_S+ABZ*I_ERI_G3xz_F2xz_Pz_S;
  abcd[528] = I_ERI_H2x2yz_F2xz_Pz_S+ABZ*I_ERI_G2x2y_F2xz_Pz_S;
  abcd[529] = I_ERI_H2xy2z_F2xz_Pz_S+ABZ*I_ERI_G2xyz_F2xz_Pz_S;
  abcd[530] = I_ERI_H2x3z_F2xz_Pz_S+ABZ*I_ERI_G2x2z_F2xz_Pz_S;
  abcd[531] = I_ERI_Hx3yz_F2xz_Pz_S+ABZ*I_ERI_Gx3y_F2xz_Pz_S;
  abcd[532] = I_ERI_Hx2y2z_F2xz_Pz_S+ABZ*I_ERI_Gx2yz_F2xz_Pz_S;
  abcd[533] = I_ERI_Hxy3z_F2xz_Pz_S+ABZ*I_ERI_Gxy2z_F2xz_Pz_S;
  abcd[534] = I_ERI_Hx4z_F2xz_Pz_S+ABZ*I_ERI_Gx3z_F2xz_Pz_S;
  abcd[535] = I_ERI_H4yz_F2xz_Pz_S+ABZ*I_ERI_G4y_F2xz_Pz_S;
  abcd[536] = I_ERI_H3y2z_F2xz_Pz_S+ABZ*I_ERI_G3yz_F2xz_Pz_S;
  abcd[537] = I_ERI_H2y3z_F2xz_Pz_S+ABZ*I_ERI_G2y2z_F2xz_Pz_S;
  abcd[538] = I_ERI_Hy4z_F2xz_Pz_S+ABZ*I_ERI_Gy3z_F2xz_Pz_S;
  abcd[539] = I_ERI_H5z_F2xz_Pz_S+ABZ*I_ERI_G4z_F2xz_Pz_S;
  abcd[540] = I_ERI_H5x_F3y_Pz_S+ABX*I_ERI_G4x_F3y_Pz_S;
  abcd[541] = I_ERI_H4xy_F3y_Pz_S+ABX*I_ERI_G3xy_F3y_Pz_S;
  abcd[542] = I_ERI_H4xz_F3y_Pz_S+ABX*I_ERI_G3xz_F3y_Pz_S;
  abcd[543] = I_ERI_H3x2y_F3y_Pz_S+ABX*I_ERI_G2x2y_F3y_Pz_S;
  abcd[544] = I_ERI_H3xyz_F3y_Pz_S+ABX*I_ERI_G2xyz_F3y_Pz_S;
  abcd[545] = I_ERI_H3x2z_F3y_Pz_S+ABX*I_ERI_G2x2z_F3y_Pz_S;
  abcd[546] = I_ERI_H2x3y_F3y_Pz_S+ABX*I_ERI_Gx3y_F3y_Pz_S;
  abcd[547] = I_ERI_H2x2yz_F3y_Pz_S+ABX*I_ERI_Gx2yz_F3y_Pz_S;
  abcd[548] = I_ERI_H2xy2z_F3y_Pz_S+ABX*I_ERI_Gxy2z_F3y_Pz_S;
  abcd[549] = I_ERI_H2x3z_F3y_Pz_S+ABX*I_ERI_Gx3z_F3y_Pz_S;
  abcd[550] = I_ERI_Hx4y_F3y_Pz_S+ABX*I_ERI_G4y_F3y_Pz_S;
  abcd[551] = I_ERI_Hx3yz_F3y_Pz_S+ABX*I_ERI_G3yz_F3y_Pz_S;
  abcd[552] = I_ERI_Hx2y2z_F3y_Pz_S+ABX*I_ERI_G2y2z_F3y_Pz_S;
  abcd[553] = I_ERI_Hxy3z_F3y_Pz_S+ABX*I_ERI_Gy3z_F3y_Pz_S;
  abcd[554] = I_ERI_Hx4z_F3y_Pz_S+ABX*I_ERI_G4z_F3y_Pz_S;
  abcd[555] = I_ERI_H4xz_Fx2y_Pz_S+ABZ*I_ERI_G4x_Fx2y_Pz_S;
  abcd[556] = I_ERI_H3xyz_Fx2y_Pz_S+ABZ*I_ERI_G3xy_Fx2y_Pz_S;
  abcd[557] = I_ERI_H3x2z_Fx2y_Pz_S+ABZ*I_ERI_G3xz_Fx2y_Pz_S;
  abcd[558] = I_ERI_H2x2yz_Fx2y_Pz_S+ABZ*I_ERI_G2x2y_Fx2y_Pz_S;
  abcd[559] = I_ERI_H2xy2z_Fx2y_Pz_S+ABZ*I_ERI_G2xyz_Fx2y_Pz_S;
  abcd[560] = I_ERI_H2x3z_Fx2y_Pz_S+ABZ*I_ERI_G2x2z_Fx2y_Pz_S;
  abcd[561] = I_ERI_Hx3yz_Fx2y_Pz_S+ABZ*I_ERI_Gx3y_Fx2y_Pz_S;
  abcd[562] = I_ERI_Hx2y2z_Fx2y_Pz_S+ABZ*I_ERI_Gx2yz_Fx2y_Pz_S;
  abcd[563] = I_ERI_Hxy3z_Fx2y_Pz_S+ABZ*I_ERI_Gxy2z_Fx2y_Pz_S;
  abcd[564] = I_ERI_Hx4z_Fx2y_Pz_S+ABZ*I_ERI_Gx3z_Fx2y_Pz_S;
  abcd[565] = I_ERI_H4yz_Fx2y_Pz_S+ABZ*I_ERI_G4y_Fx2y_Pz_S;
  abcd[566] = I_ERI_H3y2z_Fx2y_Pz_S+ABZ*I_ERI_G3yz_Fx2y_Pz_S;
  abcd[567] = I_ERI_H2y3z_Fx2y_Pz_S+ABZ*I_ERI_G2y2z_Fx2y_Pz_S;
  abcd[568] = I_ERI_Hy4z_Fx2y_Pz_S+ABZ*I_ERI_Gy3z_Fx2y_Pz_S;
  abcd[569] = I_ERI_H5z_Fx2y_Pz_S+ABZ*I_ERI_G4z_Fx2y_Pz_S;
  abcd[570] = I_ERI_H4xy_Fx2z_Pz_S+ABY*I_ERI_G4x_Fx2z_Pz_S;
  abcd[571] = I_ERI_H3x2y_Fx2z_Pz_S+ABY*I_ERI_G3xy_Fx2z_Pz_S;
  abcd[572] = I_ERI_H3xyz_Fx2z_Pz_S+ABY*I_ERI_G3xz_Fx2z_Pz_S;
  abcd[573] = I_ERI_H2x3y_Fx2z_Pz_S+ABY*I_ERI_G2x2y_Fx2z_Pz_S;
  abcd[574] = I_ERI_H2x2yz_Fx2z_Pz_S+ABY*I_ERI_G2xyz_Fx2z_Pz_S;
  abcd[575] = I_ERI_H2xy2z_Fx2z_Pz_S+ABY*I_ERI_G2x2z_Fx2z_Pz_S;
  abcd[576] = I_ERI_Hx4y_Fx2z_Pz_S+ABY*I_ERI_Gx3y_Fx2z_Pz_S;
  abcd[577] = I_ERI_Hx3yz_Fx2z_Pz_S+ABY*I_ERI_Gx2yz_Fx2z_Pz_S;
  abcd[578] = I_ERI_Hx2y2z_Fx2z_Pz_S+ABY*I_ERI_Gxy2z_Fx2z_Pz_S;
  abcd[579] = I_ERI_Hxy3z_Fx2z_Pz_S+ABY*I_ERI_Gx3z_Fx2z_Pz_S;
  abcd[580] = I_ERI_H5y_Fx2z_Pz_S+ABY*I_ERI_G4y_Fx2z_Pz_S;
  abcd[581] = I_ERI_H4yz_Fx2z_Pz_S+ABY*I_ERI_G3yz_Fx2z_Pz_S;
  abcd[582] = I_ERI_H3y2z_Fx2z_Pz_S+ABY*I_ERI_G2y2z_Fx2z_Pz_S;
  abcd[583] = I_ERI_H2y3z_Fx2z_Pz_S+ABY*I_ERI_Gy3z_Fx2z_Pz_S;
  abcd[584] = I_ERI_Hy4z_Fx2z_Pz_S+ABY*I_ERI_G4z_Fx2z_Pz_S;
  abcd[585] = I_ERI_H5x_F3z_Pz_S+ABX*I_ERI_G4x_F3z_Pz_S;
  abcd[586] = I_ERI_H4xy_F3z_Pz_S+ABX*I_ERI_G3xy_F3z_Pz_S;
  abcd[587] = I_ERI_H4xz_F3z_Pz_S+ABX*I_ERI_G3xz_F3z_Pz_S;
  abcd[588] = I_ERI_H3x2y_F3z_Pz_S+ABX*I_ERI_G2x2y_F3z_Pz_S;
  abcd[589] = I_ERI_H3xyz_F3z_Pz_S+ABX*I_ERI_G2xyz_F3z_Pz_S;
  abcd[590] = I_ERI_H3x2z_F3z_Pz_S+ABX*I_ERI_G2x2z_F3z_Pz_S;
  abcd[591] = I_ERI_H2x3y_F3z_Pz_S+ABX*I_ERI_Gx3y_F3z_Pz_S;
  abcd[592] = I_ERI_H2x2yz_F3z_Pz_S+ABX*I_ERI_Gx2yz_F3z_Pz_S;
  abcd[593] = I_ERI_H2xy2z_F3z_Pz_S+ABX*I_ERI_Gxy2z_F3z_Pz_S;
  abcd[594] = I_ERI_H2x3z_F3z_Pz_S+ABX*I_ERI_Gx3z_F3z_Pz_S;
  abcd[595] = I_ERI_Hx4y_F3z_Pz_S+ABX*I_ERI_G4y_F3z_Pz_S;
  abcd[596] = I_ERI_Hx3yz_F3z_Pz_S+ABX*I_ERI_G3yz_F3z_Pz_S;
  abcd[597] = I_ERI_Hx2y2z_F3z_Pz_S+ABX*I_ERI_G2y2z_F3z_Pz_S;
  abcd[598] = I_ERI_Hxy3z_F3z_Pz_S+ABX*I_ERI_Gy3z_F3z_Pz_S;
  abcd[599] = I_ERI_Hx4z_F3z_Pz_S+ABX*I_ERI_G4z_F3z_Pz_S;
  abcd[600] = I_ERI_H4xy_F3y_Pz_S+ABY*I_ERI_G4x_F3y_Pz_S;
  abcd[601] = I_ERI_H3x2y_F3y_Pz_S+ABY*I_ERI_G3xy_F3y_Pz_S;
  abcd[602] = I_ERI_H3xyz_F3y_Pz_S+ABY*I_ERI_G3xz_F3y_Pz_S;
  abcd[603] = I_ERI_H2x3y_F3y_Pz_S+ABY*I_ERI_G2x2y_F3y_Pz_S;
  abcd[604] = I_ERI_H2x2yz_F3y_Pz_S+ABY*I_ERI_G2xyz_F3y_Pz_S;
  abcd[605] = I_ERI_H2xy2z_F3y_Pz_S+ABY*I_ERI_G2x2z_F3y_Pz_S;
  abcd[606] = I_ERI_Hx4y_F3y_Pz_S+ABY*I_ERI_Gx3y_F3y_Pz_S;
  abcd[607] = I_ERI_Hx3yz_F3y_Pz_S+ABY*I_ERI_Gx2yz_F3y_Pz_S;
  abcd[608] = I_ERI_Hx2y2z_F3y_Pz_S+ABY*I_ERI_Gxy2z_F3y_Pz_S;
  abcd[609] = I_ERI_Hxy3z_F3y_Pz_S+ABY*I_ERI_Gx3z_F3y_Pz_S;
  abcd[610] = I_ERI_H5y_F3y_Pz_S+ABY*I_ERI_G4y_F3y_Pz_S;
  abcd[611] = I_ERI_H4yz_F3y_Pz_S+ABY*I_ERI_G3yz_F3y_Pz_S;
  abcd[612] = I_ERI_H3y2z_F3y_Pz_S+ABY*I_ERI_G2y2z_F3y_Pz_S;
  abcd[613] = I_ERI_H2y3z_F3y_Pz_S+ABY*I_ERI_Gy3z_F3y_Pz_S;
  abcd[614] = I_ERI_Hy4z_F3y_Pz_S+ABY*I_ERI_G4z_F3y_Pz_S;
  abcd[615] = I_ERI_H4xz_F3y_Pz_S+ABZ*I_ERI_G4x_F3y_Pz_S;
  abcd[616] = I_ERI_H3xyz_F3y_Pz_S+ABZ*I_ERI_G3xy_F3y_Pz_S;
  abcd[617] = I_ERI_H3x2z_F3y_Pz_S+ABZ*I_ERI_G3xz_F3y_Pz_S;
  abcd[618] = I_ERI_H2x2yz_F3y_Pz_S+ABZ*I_ERI_G2x2y_F3y_Pz_S;
  abcd[619] = I_ERI_H2xy2z_F3y_Pz_S+ABZ*I_ERI_G2xyz_F3y_Pz_S;
  abcd[620] = I_ERI_H2x3z_F3y_Pz_S+ABZ*I_ERI_G2x2z_F3y_Pz_S;
  abcd[621] = I_ERI_Hx3yz_F3y_Pz_S+ABZ*I_ERI_Gx3y_F3y_Pz_S;
  abcd[622] = I_ERI_Hx2y2z_F3y_Pz_S+ABZ*I_ERI_Gx2yz_F3y_Pz_S;
  abcd[623] = I_ERI_Hxy3z_F3y_Pz_S+ABZ*I_ERI_Gxy2z_F3y_Pz_S;
  abcd[624] = I_ERI_Hx4z_F3y_Pz_S+ABZ*I_ERI_Gx3z_F3y_Pz_S;
  abcd[625] = I_ERI_H4yz_F3y_Pz_S+ABZ*I_ERI_G4y_F3y_Pz_S;
  abcd[626] = I_ERI_H3y2z_F3y_Pz_S+ABZ*I_ERI_G3yz_F3y_Pz_S;
  abcd[627] = I_ERI_H2y3z_F3y_Pz_S+ABZ*I_ERI_G2y2z_F3y_Pz_S;
  abcd[628] = I_ERI_Hy4z_F3y_Pz_S+ABZ*I_ERI_Gy3z_F3y_Pz_S;
  abcd[629] = I_ERI_H5z_F3y_Pz_S+ABZ*I_ERI_G4z_F3y_Pz_S;
  abcd[630] = I_ERI_H4xz_F2yz_Pz_S+ABZ*I_ERI_G4x_F2yz_Pz_S;
  abcd[631] = I_ERI_H3xyz_F2yz_Pz_S+ABZ*I_ERI_G3xy_F2yz_Pz_S;
  abcd[632] = I_ERI_H3x2z_F2yz_Pz_S+ABZ*I_ERI_G3xz_F2yz_Pz_S;
  abcd[633] = I_ERI_H2x2yz_F2yz_Pz_S+ABZ*I_ERI_G2x2y_F2yz_Pz_S;
  abcd[634] = I_ERI_H2xy2z_F2yz_Pz_S+ABZ*I_ERI_G2xyz_F2yz_Pz_S;
  abcd[635] = I_ERI_H2x3z_F2yz_Pz_S+ABZ*I_ERI_G2x2z_F2yz_Pz_S;
  abcd[636] = I_ERI_Hx3yz_F2yz_Pz_S+ABZ*I_ERI_Gx3y_F2yz_Pz_S;
  abcd[637] = I_ERI_Hx2y2z_F2yz_Pz_S+ABZ*I_ERI_Gx2yz_F2yz_Pz_S;
  abcd[638] = I_ERI_Hxy3z_F2yz_Pz_S+ABZ*I_ERI_Gxy2z_F2yz_Pz_S;
  abcd[639] = I_ERI_Hx4z_F2yz_Pz_S+ABZ*I_ERI_Gx3z_F2yz_Pz_S;
  abcd[640] = I_ERI_H4yz_F2yz_Pz_S+ABZ*I_ERI_G4y_F2yz_Pz_S;
  abcd[641] = I_ERI_H3y2z_F2yz_Pz_S+ABZ*I_ERI_G3yz_F2yz_Pz_S;
  abcd[642] = I_ERI_H2y3z_F2yz_Pz_S+ABZ*I_ERI_G2y2z_F2yz_Pz_S;
  abcd[643] = I_ERI_Hy4z_F2yz_Pz_S+ABZ*I_ERI_Gy3z_F2yz_Pz_S;
  abcd[644] = I_ERI_H5z_F2yz_Pz_S+ABZ*I_ERI_G4z_F2yz_Pz_S;
  abcd[645] = I_ERI_H4xy_F3z_Pz_S+ABY*I_ERI_G4x_F3z_Pz_S;
  abcd[646] = I_ERI_H3x2y_F3z_Pz_S+ABY*I_ERI_G3xy_F3z_Pz_S;
  abcd[647] = I_ERI_H3xyz_F3z_Pz_S+ABY*I_ERI_G3xz_F3z_Pz_S;
  abcd[648] = I_ERI_H2x3y_F3z_Pz_S+ABY*I_ERI_G2x2y_F3z_Pz_S;
  abcd[649] = I_ERI_H2x2yz_F3z_Pz_S+ABY*I_ERI_G2xyz_F3z_Pz_S;
  abcd[650] = I_ERI_H2xy2z_F3z_Pz_S+ABY*I_ERI_G2x2z_F3z_Pz_S;
  abcd[651] = I_ERI_Hx4y_F3z_Pz_S+ABY*I_ERI_Gx3y_F3z_Pz_S;
  abcd[652] = I_ERI_Hx3yz_F3z_Pz_S+ABY*I_ERI_Gx2yz_F3z_Pz_S;
  abcd[653] = I_ERI_Hx2y2z_F3z_Pz_S+ABY*I_ERI_Gxy2z_F3z_Pz_S;
  abcd[654] = I_ERI_Hxy3z_F3z_Pz_S+ABY*I_ERI_Gx3z_F3z_Pz_S;
  abcd[655] = I_ERI_H5y_F3z_Pz_S+ABY*I_ERI_G4y_F3z_Pz_S;
  abcd[656] = I_ERI_H4yz_F3z_Pz_S+ABY*I_ERI_G3yz_F3z_Pz_S;
  abcd[657] = I_ERI_H3y2z_F3z_Pz_S+ABY*I_ERI_G2y2z_F3z_Pz_S;
  abcd[658] = I_ERI_H2y3z_F3z_Pz_S+ABY*I_ERI_Gy3z_F3z_Pz_S;
  abcd[659] = I_ERI_Hy4z_F3z_Pz_S+ABY*I_ERI_G4z_F3z_Pz_S;
  abcd[660] = I_ERI_H4xz_F3z_Pz_S+ABZ*I_ERI_G4x_F3z_Pz_S;
  abcd[661] = I_ERI_H3xyz_F3z_Pz_S+ABZ*I_ERI_G3xy_F3z_Pz_S;
  abcd[662] = I_ERI_H3x2z_F3z_Pz_S+ABZ*I_ERI_G3xz_F3z_Pz_S;
  abcd[663] = I_ERI_H2x2yz_F3z_Pz_S+ABZ*I_ERI_G2x2y_F3z_Pz_S;
  abcd[664] = I_ERI_H2xy2z_F3z_Pz_S+ABZ*I_ERI_G2xyz_F3z_Pz_S;
  abcd[665] = I_ERI_H2x3z_F3z_Pz_S+ABZ*I_ERI_G2x2z_F3z_Pz_S;
  abcd[666] = I_ERI_Hx3yz_F3z_Pz_S+ABZ*I_ERI_Gx3y_F3z_Pz_S;
  abcd[667] = I_ERI_Hx2y2z_F3z_Pz_S+ABZ*I_ERI_Gx2yz_F3z_Pz_S;
  abcd[668] = I_ERI_Hxy3z_F3z_Pz_S+ABZ*I_ERI_Gxy2z_F3z_Pz_S;
  abcd[669] = I_ERI_Hx4z_F3z_Pz_S+ABZ*I_ERI_Gx3z_F3z_Pz_S;
  abcd[670] = I_ERI_H4yz_F3z_Pz_S+ABZ*I_ERI_G4y_F3z_Pz_S;
  abcd[671] = I_ERI_H3y2z_F3z_Pz_S+ABZ*I_ERI_G3yz_F3z_Pz_S;
  abcd[672] = I_ERI_H2y3z_F3z_Pz_S+ABZ*I_ERI_G2y2z_F3z_Pz_S;
  abcd[673] = I_ERI_Hy4z_F3z_Pz_S+ABZ*I_ERI_Gy3z_F3z_Pz_S;
  abcd[674] = I_ERI_H5z_F3z_Pz_S+ABZ*I_ERI_G4z_F3z_Pz_S;
}
