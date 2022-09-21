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

void hgp_os_eri_f_d_p_p(const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd)
{
  // check that whether we use erf(r12)/r12 form operator 
  // such setting also applied for NAI operator etc. 
  bool withErfR12 = false;
  if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;

  //
  // declare the variables as result of VRR process
  //
  Double I_ERI_H5x_S_D2x_S = 0.0E0;
  Double I_ERI_H4xy_S_D2x_S = 0.0E0;
  Double I_ERI_H4xz_S_D2x_S = 0.0E0;
  Double I_ERI_H3x2y_S_D2x_S = 0.0E0;
  Double I_ERI_H3xyz_S_D2x_S = 0.0E0;
  Double I_ERI_H3x2z_S_D2x_S = 0.0E0;
  Double I_ERI_H2x3y_S_D2x_S = 0.0E0;
  Double I_ERI_H2x2yz_S_D2x_S = 0.0E0;
  Double I_ERI_H2xy2z_S_D2x_S = 0.0E0;
  Double I_ERI_H2x3z_S_D2x_S = 0.0E0;
  Double I_ERI_Hx4y_S_D2x_S = 0.0E0;
  Double I_ERI_Hx3yz_S_D2x_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2x_S = 0.0E0;
  Double I_ERI_Hxy3z_S_D2x_S = 0.0E0;
  Double I_ERI_Hx4z_S_D2x_S = 0.0E0;
  Double I_ERI_H5y_S_D2x_S = 0.0E0;
  Double I_ERI_H4yz_S_D2x_S = 0.0E0;
  Double I_ERI_H3y2z_S_D2x_S = 0.0E0;
  Double I_ERI_H2y3z_S_D2x_S = 0.0E0;
  Double I_ERI_Hy4z_S_D2x_S = 0.0E0;
  Double I_ERI_H5z_S_D2x_S = 0.0E0;
  Double I_ERI_H5x_S_Dxy_S = 0.0E0;
  Double I_ERI_H4xy_S_Dxy_S = 0.0E0;
  Double I_ERI_H4xz_S_Dxy_S = 0.0E0;
  Double I_ERI_H3x2y_S_Dxy_S = 0.0E0;
  Double I_ERI_H3xyz_S_Dxy_S = 0.0E0;
  Double I_ERI_H3x2z_S_Dxy_S = 0.0E0;
  Double I_ERI_H2x3y_S_Dxy_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxy_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxy_S = 0.0E0;
  Double I_ERI_H2x3z_S_Dxy_S = 0.0E0;
  Double I_ERI_Hx4y_S_Dxy_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxy_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxy_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxy_S = 0.0E0;
  Double I_ERI_Hx4z_S_Dxy_S = 0.0E0;
  Double I_ERI_H5y_S_Dxy_S = 0.0E0;
  Double I_ERI_H4yz_S_Dxy_S = 0.0E0;
  Double I_ERI_H3y2z_S_Dxy_S = 0.0E0;
  Double I_ERI_H2y3z_S_Dxy_S = 0.0E0;
  Double I_ERI_Hy4z_S_Dxy_S = 0.0E0;
  Double I_ERI_H5z_S_Dxy_S = 0.0E0;
  Double I_ERI_H5x_S_Dxz_S = 0.0E0;
  Double I_ERI_H4xy_S_Dxz_S = 0.0E0;
  Double I_ERI_H4xz_S_Dxz_S = 0.0E0;
  Double I_ERI_H3x2y_S_Dxz_S = 0.0E0;
  Double I_ERI_H3xyz_S_Dxz_S = 0.0E0;
  Double I_ERI_H3x2z_S_Dxz_S = 0.0E0;
  Double I_ERI_H2x3y_S_Dxz_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Dxz_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Dxz_S = 0.0E0;
  Double I_ERI_H2x3z_S_Dxz_S = 0.0E0;
  Double I_ERI_Hx4y_S_Dxz_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Dxz_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dxz_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Dxz_S = 0.0E0;
  Double I_ERI_Hx4z_S_Dxz_S = 0.0E0;
  Double I_ERI_H5y_S_Dxz_S = 0.0E0;
  Double I_ERI_H4yz_S_Dxz_S = 0.0E0;
  Double I_ERI_H3y2z_S_Dxz_S = 0.0E0;
  Double I_ERI_H2y3z_S_Dxz_S = 0.0E0;
  Double I_ERI_Hy4z_S_Dxz_S = 0.0E0;
  Double I_ERI_H5z_S_Dxz_S = 0.0E0;
  Double I_ERI_H5x_S_D2y_S = 0.0E0;
  Double I_ERI_H4xy_S_D2y_S = 0.0E0;
  Double I_ERI_H4xz_S_D2y_S = 0.0E0;
  Double I_ERI_H3x2y_S_D2y_S = 0.0E0;
  Double I_ERI_H3xyz_S_D2y_S = 0.0E0;
  Double I_ERI_H3x2z_S_D2y_S = 0.0E0;
  Double I_ERI_H2x3y_S_D2y_S = 0.0E0;
  Double I_ERI_H2x2yz_S_D2y_S = 0.0E0;
  Double I_ERI_H2xy2z_S_D2y_S = 0.0E0;
  Double I_ERI_H2x3z_S_D2y_S = 0.0E0;
  Double I_ERI_Hx4y_S_D2y_S = 0.0E0;
  Double I_ERI_Hx3yz_S_D2y_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2y_S = 0.0E0;
  Double I_ERI_Hxy3z_S_D2y_S = 0.0E0;
  Double I_ERI_Hx4z_S_D2y_S = 0.0E0;
  Double I_ERI_H5y_S_D2y_S = 0.0E0;
  Double I_ERI_H4yz_S_D2y_S = 0.0E0;
  Double I_ERI_H3y2z_S_D2y_S = 0.0E0;
  Double I_ERI_H2y3z_S_D2y_S = 0.0E0;
  Double I_ERI_Hy4z_S_D2y_S = 0.0E0;
  Double I_ERI_H5z_S_D2y_S = 0.0E0;
  Double I_ERI_H5x_S_Dyz_S = 0.0E0;
  Double I_ERI_H4xy_S_Dyz_S = 0.0E0;
  Double I_ERI_H4xz_S_Dyz_S = 0.0E0;
  Double I_ERI_H3x2y_S_Dyz_S = 0.0E0;
  Double I_ERI_H3xyz_S_Dyz_S = 0.0E0;
  Double I_ERI_H3x2z_S_Dyz_S = 0.0E0;
  Double I_ERI_H2x3y_S_Dyz_S = 0.0E0;
  Double I_ERI_H2x2yz_S_Dyz_S = 0.0E0;
  Double I_ERI_H2xy2z_S_Dyz_S = 0.0E0;
  Double I_ERI_H2x3z_S_Dyz_S = 0.0E0;
  Double I_ERI_Hx4y_S_Dyz_S = 0.0E0;
  Double I_ERI_Hx3yz_S_Dyz_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_Dyz_S = 0.0E0;
  Double I_ERI_Hxy3z_S_Dyz_S = 0.0E0;
  Double I_ERI_Hx4z_S_Dyz_S = 0.0E0;
  Double I_ERI_H5y_S_Dyz_S = 0.0E0;
  Double I_ERI_H4yz_S_Dyz_S = 0.0E0;
  Double I_ERI_H3y2z_S_Dyz_S = 0.0E0;
  Double I_ERI_H2y3z_S_Dyz_S = 0.0E0;
  Double I_ERI_Hy4z_S_Dyz_S = 0.0E0;
  Double I_ERI_H5z_S_Dyz_S = 0.0E0;
  Double I_ERI_H5x_S_D2z_S = 0.0E0;
  Double I_ERI_H4xy_S_D2z_S = 0.0E0;
  Double I_ERI_H4xz_S_D2z_S = 0.0E0;
  Double I_ERI_H3x2y_S_D2z_S = 0.0E0;
  Double I_ERI_H3xyz_S_D2z_S = 0.0E0;
  Double I_ERI_H3x2z_S_D2z_S = 0.0E0;
  Double I_ERI_H2x3y_S_D2z_S = 0.0E0;
  Double I_ERI_H2x2yz_S_D2z_S = 0.0E0;
  Double I_ERI_H2xy2z_S_D2z_S = 0.0E0;
  Double I_ERI_H2x3z_S_D2z_S = 0.0E0;
  Double I_ERI_Hx4y_S_D2z_S = 0.0E0;
  Double I_ERI_Hx3yz_S_D2z_S = 0.0E0;
  Double I_ERI_Hx2y2z_S_D2z_S = 0.0E0;
  Double I_ERI_Hxy3z_S_D2z_S = 0.0E0;
  Double I_ERI_Hx4z_S_D2z_S = 0.0E0;
  Double I_ERI_H5y_S_D2z_S = 0.0E0;
  Double I_ERI_H4yz_S_D2z_S = 0.0E0;
  Double I_ERI_H3y2z_S_D2z_S = 0.0E0;
  Double I_ERI_H2y3z_S_D2z_S = 0.0E0;
  Double I_ERI_Hy4z_S_D2z_S = 0.0E0;
  Double I_ERI_H5z_S_D2z_S = 0.0E0;
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
  Double I_ERI_G4x_S_D2x_S = 0.0E0;
  Double I_ERI_G3xy_S_D2x_S = 0.0E0;
  Double I_ERI_G3xz_S_D2x_S = 0.0E0;
  Double I_ERI_G2x2y_S_D2x_S = 0.0E0;
  Double I_ERI_G2xyz_S_D2x_S = 0.0E0;
  Double I_ERI_G2x2z_S_D2x_S = 0.0E0;
  Double I_ERI_Gx3y_S_D2x_S = 0.0E0;
  Double I_ERI_Gx2yz_S_D2x_S = 0.0E0;
  Double I_ERI_Gxy2z_S_D2x_S = 0.0E0;
  Double I_ERI_Gx3z_S_D2x_S = 0.0E0;
  Double I_ERI_G4y_S_D2x_S = 0.0E0;
  Double I_ERI_G3yz_S_D2x_S = 0.0E0;
  Double I_ERI_G2y2z_S_D2x_S = 0.0E0;
  Double I_ERI_Gy3z_S_D2x_S = 0.0E0;
  Double I_ERI_G4z_S_D2x_S = 0.0E0;
  Double I_ERI_G4x_S_Dxy_S = 0.0E0;
  Double I_ERI_G3xy_S_Dxy_S = 0.0E0;
  Double I_ERI_G3xz_S_Dxy_S = 0.0E0;
  Double I_ERI_G2x2y_S_Dxy_S = 0.0E0;
  Double I_ERI_G2xyz_S_Dxy_S = 0.0E0;
  Double I_ERI_G2x2z_S_Dxy_S = 0.0E0;
  Double I_ERI_Gx3y_S_Dxy_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxy_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxy_S = 0.0E0;
  Double I_ERI_Gx3z_S_Dxy_S = 0.0E0;
  Double I_ERI_G4y_S_Dxy_S = 0.0E0;
  Double I_ERI_G3yz_S_Dxy_S = 0.0E0;
  Double I_ERI_G2y2z_S_Dxy_S = 0.0E0;
  Double I_ERI_Gy3z_S_Dxy_S = 0.0E0;
  Double I_ERI_G4z_S_Dxy_S = 0.0E0;
  Double I_ERI_G4x_S_Dxz_S = 0.0E0;
  Double I_ERI_G3xy_S_Dxz_S = 0.0E0;
  Double I_ERI_G3xz_S_Dxz_S = 0.0E0;
  Double I_ERI_G2x2y_S_Dxz_S = 0.0E0;
  Double I_ERI_G2xyz_S_Dxz_S = 0.0E0;
  Double I_ERI_G2x2z_S_Dxz_S = 0.0E0;
  Double I_ERI_Gx3y_S_Dxz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Dxz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Dxz_S = 0.0E0;
  Double I_ERI_Gx3z_S_Dxz_S = 0.0E0;
  Double I_ERI_G4y_S_Dxz_S = 0.0E0;
  Double I_ERI_G3yz_S_Dxz_S = 0.0E0;
  Double I_ERI_G2y2z_S_Dxz_S = 0.0E0;
  Double I_ERI_Gy3z_S_Dxz_S = 0.0E0;
  Double I_ERI_G4z_S_Dxz_S = 0.0E0;
  Double I_ERI_G4x_S_D2y_S = 0.0E0;
  Double I_ERI_G3xy_S_D2y_S = 0.0E0;
  Double I_ERI_G3xz_S_D2y_S = 0.0E0;
  Double I_ERI_G2x2y_S_D2y_S = 0.0E0;
  Double I_ERI_G2xyz_S_D2y_S = 0.0E0;
  Double I_ERI_G2x2z_S_D2y_S = 0.0E0;
  Double I_ERI_Gx3y_S_D2y_S = 0.0E0;
  Double I_ERI_Gx2yz_S_D2y_S = 0.0E0;
  Double I_ERI_Gxy2z_S_D2y_S = 0.0E0;
  Double I_ERI_Gx3z_S_D2y_S = 0.0E0;
  Double I_ERI_G4y_S_D2y_S = 0.0E0;
  Double I_ERI_G3yz_S_D2y_S = 0.0E0;
  Double I_ERI_G2y2z_S_D2y_S = 0.0E0;
  Double I_ERI_Gy3z_S_D2y_S = 0.0E0;
  Double I_ERI_G4z_S_D2y_S = 0.0E0;
  Double I_ERI_G4x_S_Dyz_S = 0.0E0;
  Double I_ERI_G3xy_S_Dyz_S = 0.0E0;
  Double I_ERI_G3xz_S_Dyz_S = 0.0E0;
  Double I_ERI_G2x2y_S_Dyz_S = 0.0E0;
  Double I_ERI_G2xyz_S_Dyz_S = 0.0E0;
  Double I_ERI_G2x2z_S_Dyz_S = 0.0E0;
  Double I_ERI_Gx3y_S_Dyz_S = 0.0E0;
  Double I_ERI_Gx2yz_S_Dyz_S = 0.0E0;
  Double I_ERI_Gxy2z_S_Dyz_S = 0.0E0;
  Double I_ERI_Gx3z_S_Dyz_S = 0.0E0;
  Double I_ERI_G4y_S_Dyz_S = 0.0E0;
  Double I_ERI_G3yz_S_Dyz_S = 0.0E0;
  Double I_ERI_G2y2z_S_Dyz_S = 0.0E0;
  Double I_ERI_Gy3z_S_Dyz_S = 0.0E0;
  Double I_ERI_G4z_S_Dyz_S = 0.0E0;
  Double I_ERI_G4x_S_D2z_S = 0.0E0;
  Double I_ERI_G3xy_S_D2z_S = 0.0E0;
  Double I_ERI_G3xz_S_D2z_S = 0.0E0;
  Double I_ERI_G2x2y_S_D2z_S = 0.0E0;
  Double I_ERI_G2xyz_S_D2z_S = 0.0E0;
  Double I_ERI_G2x2z_S_D2z_S = 0.0E0;
  Double I_ERI_Gx3y_S_D2z_S = 0.0E0;
  Double I_ERI_Gx2yz_S_D2z_S = 0.0E0;
  Double I_ERI_Gxy2z_S_D2z_S = 0.0E0;
  Double I_ERI_Gx3z_S_D2z_S = 0.0E0;
  Double I_ERI_G4y_S_D2z_S = 0.0E0;
  Double I_ERI_G3yz_S_D2z_S = 0.0E0;
  Double I_ERI_G2y2z_S_D2z_S = 0.0E0;
  Double I_ERI_Gy3z_S_D2z_S = 0.0E0;
  Double I_ERI_G4z_S_D2z_S = 0.0E0;
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
  Double I_ERI_F3x_S_D2x_S = 0.0E0;
  Double I_ERI_F2xy_S_D2x_S = 0.0E0;
  Double I_ERI_F2xz_S_D2x_S = 0.0E0;
  Double I_ERI_Fx2y_S_D2x_S = 0.0E0;
  Double I_ERI_Fxyz_S_D2x_S = 0.0E0;
  Double I_ERI_Fx2z_S_D2x_S = 0.0E0;
  Double I_ERI_F3y_S_D2x_S = 0.0E0;
  Double I_ERI_F2yz_S_D2x_S = 0.0E0;
  Double I_ERI_Fy2z_S_D2x_S = 0.0E0;
  Double I_ERI_F3z_S_D2x_S = 0.0E0;
  Double I_ERI_F3x_S_Dxy_S = 0.0E0;
  Double I_ERI_F2xy_S_Dxy_S = 0.0E0;
  Double I_ERI_F2xz_S_Dxy_S = 0.0E0;
  Double I_ERI_Fx2y_S_Dxy_S = 0.0E0;
  Double I_ERI_Fxyz_S_Dxy_S = 0.0E0;
  Double I_ERI_Fx2z_S_Dxy_S = 0.0E0;
  Double I_ERI_F3y_S_Dxy_S = 0.0E0;
  Double I_ERI_F2yz_S_Dxy_S = 0.0E0;
  Double I_ERI_Fy2z_S_Dxy_S = 0.0E0;
  Double I_ERI_F3z_S_Dxy_S = 0.0E0;
  Double I_ERI_F3x_S_Dxz_S = 0.0E0;
  Double I_ERI_F2xy_S_Dxz_S = 0.0E0;
  Double I_ERI_F2xz_S_Dxz_S = 0.0E0;
  Double I_ERI_Fx2y_S_Dxz_S = 0.0E0;
  Double I_ERI_Fxyz_S_Dxz_S = 0.0E0;
  Double I_ERI_Fx2z_S_Dxz_S = 0.0E0;
  Double I_ERI_F3y_S_Dxz_S = 0.0E0;
  Double I_ERI_F2yz_S_Dxz_S = 0.0E0;
  Double I_ERI_Fy2z_S_Dxz_S = 0.0E0;
  Double I_ERI_F3z_S_Dxz_S = 0.0E0;
  Double I_ERI_F3x_S_D2y_S = 0.0E0;
  Double I_ERI_F2xy_S_D2y_S = 0.0E0;
  Double I_ERI_F2xz_S_D2y_S = 0.0E0;
  Double I_ERI_Fx2y_S_D2y_S = 0.0E0;
  Double I_ERI_Fxyz_S_D2y_S = 0.0E0;
  Double I_ERI_Fx2z_S_D2y_S = 0.0E0;
  Double I_ERI_F3y_S_D2y_S = 0.0E0;
  Double I_ERI_F2yz_S_D2y_S = 0.0E0;
  Double I_ERI_Fy2z_S_D2y_S = 0.0E0;
  Double I_ERI_F3z_S_D2y_S = 0.0E0;
  Double I_ERI_F3x_S_Dyz_S = 0.0E0;
  Double I_ERI_F2xy_S_Dyz_S = 0.0E0;
  Double I_ERI_F2xz_S_Dyz_S = 0.0E0;
  Double I_ERI_Fx2y_S_Dyz_S = 0.0E0;
  Double I_ERI_Fxyz_S_Dyz_S = 0.0E0;
  Double I_ERI_Fx2z_S_Dyz_S = 0.0E0;
  Double I_ERI_F3y_S_Dyz_S = 0.0E0;
  Double I_ERI_F2yz_S_Dyz_S = 0.0E0;
  Double I_ERI_Fy2z_S_Dyz_S = 0.0E0;
  Double I_ERI_F3z_S_Dyz_S = 0.0E0;
  Double I_ERI_F3x_S_D2z_S = 0.0E0;
  Double I_ERI_F2xy_S_D2z_S = 0.0E0;
  Double I_ERI_F2xz_S_D2z_S = 0.0E0;
  Double I_ERI_Fx2y_S_D2z_S = 0.0E0;
  Double I_ERI_Fxyz_S_D2z_S = 0.0E0;
  Double I_ERI_Fx2z_S_D2z_S = 0.0E0;
  Double I_ERI_F3y_S_D2z_S = 0.0E0;
  Double I_ERI_F2yz_S_D2z_S = 0.0E0;
  Double I_ERI_Fy2z_S_D2z_S = 0.0E0;
  Double I_ERI_F3z_S_D2z_S = 0.0E0;
  Double I_ERI_F3x_S_Px_S = 0.0E0;
  Double I_ERI_F2xy_S_Px_S = 0.0E0;
  Double I_ERI_F2xz_S_Px_S = 0.0E0;
  Double I_ERI_Fx2y_S_Px_S = 0.0E0;
  Double I_ERI_Fxyz_S_Px_S = 0.0E0;
  Double I_ERI_Fx2z_S_Px_S = 0.0E0;
  Double I_ERI_F3y_S_Px_S = 0.0E0;
  Double I_ERI_F2yz_S_Px_S = 0.0E0;
  Double I_ERI_Fy2z_S_Px_S = 0.0E0;
  Double I_ERI_F3z_S_Px_S = 0.0E0;
  Double I_ERI_F3x_S_Py_S = 0.0E0;
  Double I_ERI_F2xy_S_Py_S = 0.0E0;
  Double I_ERI_F2xz_S_Py_S = 0.0E0;
  Double I_ERI_Fx2y_S_Py_S = 0.0E0;
  Double I_ERI_Fxyz_S_Py_S = 0.0E0;
  Double I_ERI_Fx2z_S_Py_S = 0.0E0;
  Double I_ERI_F3y_S_Py_S = 0.0E0;
  Double I_ERI_F2yz_S_Py_S = 0.0E0;
  Double I_ERI_Fy2z_S_Py_S = 0.0E0;
  Double I_ERI_F3z_S_Py_S = 0.0E0;
  Double I_ERI_F3x_S_Pz_S = 0.0E0;
  Double I_ERI_F2xy_S_Pz_S = 0.0E0;
  Double I_ERI_F2xz_S_Pz_S = 0.0E0;
  Double I_ERI_Fx2y_S_Pz_S = 0.0E0;
  Double I_ERI_Fxyz_S_Pz_S = 0.0E0;
  Double I_ERI_Fx2z_S_Pz_S = 0.0E0;
  Double I_ERI_F3y_S_Pz_S = 0.0E0;
  Double I_ERI_F2yz_S_Pz_S = 0.0E0;
  Double I_ERI_Fy2z_S_Pz_S = 0.0E0;
  Double I_ERI_F3z_S_Pz_S = 0.0E0;

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
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M5
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M5
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M4_vrr = PAX*I_ERI_D2x_S_S_S_M4_vrr+WPX*I_ERI_D2x_S_S_S_M5_vrr+2*oned2z*I_ERI_Px_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M5_vrr;
      Double I_ERI_F2xy_S_S_S_M4_vrr = PAY*I_ERI_D2x_S_S_S_M4_vrr+WPY*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_F2xz_S_S_S_M4_vrr = PAZ*I_ERI_D2x_S_S_S_M4_vrr+WPZ*I_ERI_D2x_S_S_S_M5_vrr;
      Double I_ERI_F3y_S_S_S_M4_vrr = PAY*I_ERI_D2y_S_S_S_M4_vrr+WPY*I_ERI_D2y_S_S_S_M5_vrr+2*oned2z*I_ERI_Py_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M5_vrr;
      Double I_ERI_F2yz_S_S_S_M4_vrr = PAZ*I_ERI_D2y_S_S_S_M4_vrr+WPZ*I_ERI_D2y_S_S_S_M5_vrr;
      Double I_ERI_F3z_S_S_S_M4_vrr = PAZ*I_ERI_D2z_S_S_S_M4_vrr+WPZ*I_ERI_D2z_S_S_S_M5_vrr+2*oned2z*I_ERI_Pz_S_S_S_M4_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M5_vrr;

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
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_S_S_S_S_M4
       ************************************************************/
      Double I_ERI_D2x_S_S_S_M3_vrr = PAX*I_ERI_Px_S_S_S_M3_vrr+WPX*I_ERI_Px_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_Dxy_S_S_S_M3_vrr = PAY*I_ERI_Px_S_S_S_M3_vrr+WPY*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_D2y_S_S_S_M3_vrr = PAY*I_ERI_Py_S_S_S_M3_vrr+WPY*I_ERI_Py_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;
      Double I_ERI_D2z_S_S_S_M3_vrr = PAZ*I_ERI_Pz_S_S_S_M3_vrr+WPZ*I_ERI_Pz_S_S_S_M4_vrr+oned2z*I_ERI_S_S_S_S_M3_vrr-rhod2zsq*I_ERI_S_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_P_S_S_S_M4
       ************************************************************/
      Double I_ERI_F3x_S_S_S_M3_vrr = PAX*I_ERI_D2x_S_S_S_M3_vrr+WPX*I_ERI_D2x_S_S_S_M4_vrr+2*oned2z*I_ERI_Px_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Px_S_S_S_M4_vrr;
      Double I_ERI_F2xy_S_S_S_M3_vrr = PAY*I_ERI_D2x_S_S_S_M3_vrr+WPY*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_F2xz_S_S_S_M3_vrr = PAZ*I_ERI_D2x_S_S_S_M3_vrr+WPZ*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Fx2y_S_S_S_M3_vrr = PAX*I_ERI_D2y_S_S_S_M3_vrr+WPX*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Fx2z_S_S_S_M3_vrr = PAX*I_ERI_D2z_S_S_S_M3_vrr+WPX*I_ERI_D2z_S_S_S_M4_vrr;
      Double I_ERI_F3y_S_S_S_M3_vrr = PAY*I_ERI_D2y_S_S_S_M3_vrr+WPY*I_ERI_D2y_S_S_S_M4_vrr+2*oned2z*I_ERI_Py_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Py_S_S_S_M4_vrr;
      Double I_ERI_F2yz_S_S_S_M3_vrr = PAZ*I_ERI_D2y_S_S_S_M3_vrr+WPZ*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_F3z_S_S_S_M3_vrr = PAZ*I_ERI_D2z_S_S_S_M3_vrr+WPZ*I_ERI_D2z_S_S_S_M4_vrr+2*oned2z*I_ERI_Pz_S_S_S_M3_vrr-2*rhod2zsq*I_ERI_Pz_S_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_G_S_S_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M4
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_D_S_S_S_M4
       ************************************************************/
      Double I_ERI_G4x_S_S_S_M3_vrr = PAX*I_ERI_F3x_S_S_S_M3_vrr+WPX*I_ERI_F3x_S_S_S_M4_vrr+3*oned2z*I_ERI_D2x_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G3xy_S_S_S_M3_vrr = PAY*I_ERI_F3x_S_S_S_M3_vrr+WPY*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G3xz_S_S_S_M3_vrr = PAZ*I_ERI_F3x_S_S_S_M3_vrr+WPZ*I_ERI_F3x_S_S_S_M4_vrr;
      Double I_ERI_G2x2y_S_S_S_M3_vrr = PAY*I_ERI_F2xy_S_S_S_M3_vrr+WPY*I_ERI_F2xy_S_S_S_M4_vrr+oned2z*I_ERI_D2x_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_G2x2z_S_S_S_M3_vrr = PAZ*I_ERI_F2xz_S_S_S_M3_vrr+WPZ*I_ERI_F2xz_S_S_S_M4_vrr+oned2z*I_ERI_D2x_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2x_S_S_S_M4_vrr;
      Double I_ERI_Gx3y_S_S_S_M3_vrr = PAX*I_ERI_F3y_S_S_S_M3_vrr+WPX*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_Gx3z_S_S_S_M3_vrr = PAX*I_ERI_F3z_S_S_S_M3_vrr+WPX*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4y_S_S_S_M3_vrr = PAY*I_ERI_F3y_S_S_S_M3_vrr+WPY*I_ERI_F3y_S_S_S_M4_vrr+3*oned2z*I_ERI_D2y_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_G3yz_S_S_S_M3_vrr = PAZ*I_ERI_F3y_S_S_S_M3_vrr+WPZ*I_ERI_F3y_S_S_S_M4_vrr;
      Double I_ERI_G2y2z_S_S_S_M3_vrr = PAZ*I_ERI_F2yz_S_S_S_M3_vrr+WPZ*I_ERI_F2yz_S_S_S_M4_vrr+oned2z*I_ERI_D2y_S_S_S_M3_vrr-rhod2zsq*I_ERI_D2y_S_S_S_M4_vrr;
      Double I_ERI_Gy3z_S_S_S_M3_vrr = PAY*I_ERI_F3z_S_S_S_M3_vrr+WPY*I_ERI_F3z_S_S_S_M4_vrr;
      Double I_ERI_G4z_S_S_S_M3_vrr = PAZ*I_ERI_F3z_S_S_S_M3_vrr+WPZ*I_ERI_F3z_S_S_S_M4_vrr+3*oned2z*I_ERI_D2z_S_S_S_M3_vrr-3*rhod2zsq*I_ERI_D2z_S_S_S_M4_vrr;

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
       * shell quartet name: SQ_ERI_H_S_S_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M3
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M3
       ************************************************************/
      Double I_ERI_H5x_S_S_S_M2_vrr = PAX*I_ERI_G4x_S_S_S_M2_vrr+WPX*I_ERI_G4x_S_S_S_M3_vrr+4*oned2z*I_ERI_F3x_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H4xy_S_S_S_M2_vrr = PAY*I_ERI_G4x_S_S_S_M2_vrr+WPY*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_H4xz_S_S_S_M2_vrr = PAZ*I_ERI_G4x_S_S_S_M2_vrr+WPZ*I_ERI_G4x_S_S_S_M3_vrr;
      Double I_ERI_H3x2y_S_S_S_M2_vrr = PAY*I_ERI_G3xy_S_S_S_M2_vrr+WPY*I_ERI_G3xy_S_S_S_M3_vrr+oned2z*I_ERI_F3x_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H3xyz_S_S_S_M2_vrr = PAZ*I_ERI_G3xy_S_S_S_M2_vrr+WPZ*I_ERI_G3xy_S_S_S_M3_vrr;
      Double I_ERI_H3x2z_S_S_S_M2_vrr = PAZ*I_ERI_G3xz_S_S_S_M2_vrr+WPZ*I_ERI_G3xz_S_S_S_M3_vrr+oned2z*I_ERI_F3x_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3x_S_S_S_M3_vrr;
      Double I_ERI_H2x3y_S_S_S_M2_vrr = PAX*I_ERI_Gx3y_S_S_S_M2_vrr+WPX*I_ERI_Gx3y_S_S_S_M3_vrr+oned2z*I_ERI_F3y_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H2x2yz_S_S_S_M2_vrr = PAZ*I_ERI_G2x2y_S_S_S_M2_vrr+WPZ*I_ERI_G2x2y_S_S_S_M3_vrr;
      Double I_ERI_H2xy2z_S_S_S_M2_vrr = PAY*I_ERI_G2x2z_S_S_S_M2_vrr+WPY*I_ERI_G2x2z_S_S_S_M3_vrr;
      Double I_ERI_H2x3z_S_S_S_M2_vrr = PAX*I_ERI_Gx3z_S_S_S_M2_vrr+WPX*I_ERI_Gx3z_S_S_S_M3_vrr+oned2z*I_ERI_F3z_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_Hx4y_S_S_S_M2_vrr = PAX*I_ERI_G4y_S_S_S_M2_vrr+WPX*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_Hx3yz_S_S_S_M2_vrr = PAZ*I_ERI_Gx3y_S_S_S_M2_vrr+WPZ*I_ERI_Gx3y_S_S_S_M3_vrr;
      Double I_ERI_Hx2y2z_S_S_S_M2_vrr = PAX*I_ERI_G2y2z_S_S_S_M2_vrr+WPX*I_ERI_G2y2z_S_S_S_M3_vrr;
      Double I_ERI_Hxy3z_S_S_S_M2_vrr = PAY*I_ERI_Gx3z_S_S_S_M2_vrr+WPY*I_ERI_Gx3z_S_S_S_M3_vrr;
      Double I_ERI_Hx4z_S_S_S_M2_vrr = PAX*I_ERI_G4z_S_S_S_M2_vrr+WPX*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_H5y_S_S_S_M2_vrr = PAY*I_ERI_G4y_S_S_S_M2_vrr+WPY*I_ERI_G4y_S_S_S_M3_vrr+4*oned2z*I_ERI_F3y_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H4yz_S_S_S_M2_vrr = PAZ*I_ERI_G4y_S_S_S_M2_vrr+WPZ*I_ERI_G4y_S_S_S_M3_vrr;
      Double I_ERI_H3y2z_S_S_S_M2_vrr = PAZ*I_ERI_G3yz_S_S_S_M2_vrr+WPZ*I_ERI_G3yz_S_S_S_M3_vrr+oned2z*I_ERI_F3y_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3y_S_S_S_M3_vrr;
      Double I_ERI_H2y3z_S_S_S_M2_vrr = PAY*I_ERI_Gy3z_S_S_S_M2_vrr+WPY*I_ERI_Gy3z_S_S_S_M3_vrr+oned2z*I_ERI_F3z_S_S_S_M2_vrr-rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;
      Double I_ERI_Hy4z_S_S_S_M2_vrr = PAY*I_ERI_G4z_S_S_S_M2_vrr+WPY*I_ERI_G4z_S_S_S_M3_vrr;
      Double I_ERI_H5z_S_S_S_M2_vrr = PAZ*I_ERI_G4z_S_S_S_M2_vrr+WPZ*I_ERI_G4z_S_S_S_M3_vrr+4*oned2z*I_ERI_F3z_S_S_S_M2_vrr-4*rhod2zsq*I_ERI_F3z_S_S_S_M3_vrr;

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
       * shell quartet name: SQ_ERI_G_S_P_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_F_S_S_S_M2
       ************************************************************/
      Double I_ERI_G4x_S_Px_S_M1_vrr = QCX*I_ERI_G4x_S_S_S_M1_vrr+WQX*I_ERI_G4x_S_S_S_M2_vrr+4*oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_Px_S_M1_vrr = QCX*I_ERI_G3xy_S_S_S_M1_vrr+WQX*I_ERI_G3xy_S_S_S_M2_vrr+3*oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_Px_S_M1_vrr = QCX*I_ERI_G3xz_S_S_S_M1_vrr+WQX*I_ERI_G3xz_S_S_S_M2_vrr+3*oned2k*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_Px_S_M1_vrr = QCX*I_ERI_G2x2y_S_S_S_M1_vrr+WQX*I_ERI_G2x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_Px_S_M1_vrr = QCX*I_ERI_G2xyz_S_S_S_M1_vrr+WQX*I_ERI_G2xyz_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_Px_S_M1_vrr = QCX*I_ERI_G2x2z_S_S_S_M1_vrr+WQX*I_ERI_G2x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_Px_S_M1_vrr = QCX*I_ERI_Gx3y_S_S_S_M1_vrr+WQX*I_ERI_Gx3y_S_S_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Px_S_M1_vrr = QCX*I_ERI_Gx2yz_S_S_S_M1_vrr+WQX*I_ERI_Gx2yz_S_S_S_M2_vrr+oned2k*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Px_S_M1_vrr = QCX*I_ERI_Gxy2z_S_S_S_M1_vrr+WQX*I_ERI_Gxy2z_S_S_S_M2_vrr+oned2k*I_ERI_Fy2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_Px_S_M1_vrr = QCX*I_ERI_Gx3z_S_S_S_M1_vrr+WQX*I_ERI_Gx3z_S_S_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_Px_S_M1_vrr = QCX*I_ERI_G4y_S_S_S_M1_vrr+WQX*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_Px_S_M1_vrr = QCX*I_ERI_G3yz_S_S_S_M1_vrr+WQX*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_Px_S_M1_vrr = QCX*I_ERI_G2y2z_S_S_S_M1_vrr+WQX*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_Px_S_M1_vrr = QCX*I_ERI_Gy3z_S_S_S_M1_vrr+WQX*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_Px_S_M1_vrr = QCX*I_ERI_G4z_S_S_S_M1_vrr+WQX*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_G4x_S_Py_S_M1_vrr = QCY*I_ERI_G4x_S_S_S_M1_vrr+WQY*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_Py_S_M1_vrr = QCY*I_ERI_G3xy_S_S_S_M1_vrr+WQY*I_ERI_G3xy_S_S_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_Py_S_M1_vrr = QCY*I_ERI_G3xz_S_S_S_M1_vrr+WQY*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_Py_S_M1_vrr = QCY*I_ERI_G2x2y_S_S_S_M1_vrr+WQY*I_ERI_G2x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_Py_S_M1_vrr = QCY*I_ERI_G2xyz_S_S_S_M1_vrr+WQY*I_ERI_G2xyz_S_S_S_M2_vrr+oned2k*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_Py_S_M1_vrr = QCY*I_ERI_G2x2z_S_S_S_M1_vrr+WQY*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_Py_S_M1_vrr = QCY*I_ERI_Gx3y_S_S_S_M1_vrr+WQY*I_ERI_Gx3y_S_S_S_M2_vrr+3*oned2k*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Py_S_M1_vrr = QCY*I_ERI_Gx2yz_S_S_S_M1_vrr+WQY*I_ERI_Gx2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Py_S_M1_vrr = QCY*I_ERI_Gxy2z_S_S_S_M1_vrr+WQY*I_ERI_Gxy2z_S_S_S_M2_vrr+oned2k*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_Py_S_M1_vrr = QCY*I_ERI_Gx3z_S_S_S_M1_vrr+WQY*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_Py_S_M1_vrr = QCY*I_ERI_G4y_S_S_S_M1_vrr+WQY*I_ERI_G4y_S_S_S_M2_vrr+4*oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_Py_S_M1_vrr = QCY*I_ERI_G3yz_S_S_S_M1_vrr+WQY*I_ERI_G3yz_S_S_S_M2_vrr+3*oned2k*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_Py_S_M1_vrr = QCY*I_ERI_G2y2z_S_S_S_M1_vrr+WQY*I_ERI_G2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fy2z_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_Py_S_M1_vrr = QCY*I_ERI_Gy3z_S_S_S_M1_vrr+WQY*I_ERI_Gy3z_S_S_S_M2_vrr+oned2k*I_ERI_F3z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_Py_S_M1_vrr = QCY*I_ERI_G4z_S_S_S_M1_vrr+WQY*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_G4x_S_Pz_S_M1_vrr = QCZ*I_ERI_G4x_S_S_S_M1_vrr+WQZ*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_G3xy_S_Pz_S_M1_vrr = QCZ*I_ERI_G3xy_S_S_S_M1_vrr+WQZ*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_G3xz_S_Pz_S_M1_vrr = QCZ*I_ERI_G3xz_S_S_S_M1_vrr+WQZ*I_ERI_G3xz_S_S_S_M2_vrr+oned2k*I_ERI_F3x_S_S_S_M2_vrr;
      Double I_ERI_G2x2y_S_Pz_S_M1_vrr = QCZ*I_ERI_G2x2y_S_S_S_M1_vrr+WQZ*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_G2xyz_S_Pz_S_M1_vrr = QCZ*I_ERI_G2xyz_S_S_S_M1_vrr+WQZ*I_ERI_G2xyz_S_S_S_M2_vrr+oned2k*I_ERI_F2xy_S_S_S_M2_vrr;
      Double I_ERI_G2x2z_S_Pz_S_M1_vrr = QCZ*I_ERI_G2x2z_S_S_S_M1_vrr+WQZ*I_ERI_G2x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_F2xz_S_S_S_M2_vrr;
      Double I_ERI_Gx3y_S_Pz_S_M1_vrr = QCZ*I_ERI_Gx3y_S_S_S_M1_vrr+WQZ*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_Gx2yz_S_Pz_S_M1_vrr = QCZ*I_ERI_Gx2yz_S_S_S_M1_vrr+WQZ*I_ERI_Gx2yz_S_S_S_M2_vrr+oned2k*I_ERI_Fx2y_S_S_S_M2_vrr;
      Double I_ERI_Gxy2z_S_Pz_S_M1_vrr = QCZ*I_ERI_Gxy2z_S_S_S_M1_vrr+WQZ*I_ERI_Gxy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Fxyz_S_S_S_M2_vrr;
      Double I_ERI_Gx3z_S_Pz_S_M1_vrr = QCZ*I_ERI_Gx3z_S_S_S_M1_vrr+WQZ*I_ERI_Gx3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Fx2z_S_S_S_M2_vrr;
      Double I_ERI_G4y_S_Pz_S_M1_vrr = QCZ*I_ERI_G4y_S_S_S_M1_vrr+WQZ*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_G3yz_S_Pz_S_M1_vrr = QCZ*I_ERI_G3yz_S_S_S_M1_vrr+WQZ*I_ERI_G3yz_S_S_S_M2_vrr+oned2k*I_ERI_F3y_S_S_S_M2_vrr;
      Double I_ERI_G2y2z_S_Pz_S_M1_vrr = QCZ*I_ERI_G2y2z_S_S_S_M1_vrr+WQZ*I_ERI_G2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_F2yz_S_S_S_M2_vrr;
      Double I_ERI_Gy3z_S_Pz_S_M1_vrr = QCZ*I_ERI_Gy3z_S_S_S_M1_vrr+WQZ*I_ERI_Gy3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Fy2z_S_S_S_M2_vrr;
      Double I_ERI_G4z_S_Pz_S_M1_vrr = QCZ*I_ERI_G4z_S_S_S_M1_vrr+WQZ*I_ERI_G4z_S_S_S_M2_vrr+4*oned2k*I_ERI_F3z_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_H_S_P_S_M1
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M2
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M2
       ************************************************************/
      Double I_ERI_H5x_S_Px_S_M1_vrr = QCX*I_ERI_H5x_S_S_S_M1_vrr+WQX*I_ERI_H5x_S_S_S_M2_vrr+5*oned2k*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H4xy_S_Px_S_M1_vrr = QCX*I_ERI_H4xy_S_S_S_M1_vrr+WQX*I_ERI_H4xy_S_S_S_M2_vrr+4*oned2k*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_H4xz_S_Px_S_M1_vrr = QCX*I_ERI_H4xz_S_S_S_M1_vrr+WQX*I_ERI_H4xz_S_S_S_M2_vrr+4*oned2k*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_H3x2y_S_Px_S_M1_vrr = QCX*I_ERI_H3x2y_S_S_S_M1_vrr+WQX*I_ERI_H3x2y_S_S_S_M2_vrr+3*oned2k*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_H3xyz_S_Px_S_M1_vrr = QCX*I_ERI_H3xyz_S_S_S_M1_vrr+WQX*I_ERI_H3xyz_S_S_S_M2_vrr+3*oned2k*I_ERI_G2xyz_S_S_S_M2_vrr;
      Double I_ERI_H3x2z_S_Px_S_M1_vrr = QCX*I_ERI_H3x2z_S_S_S_M1_vrr+WQX*I_ERI_H3x2z_S_S_S_M2_vrr+3*oned2k*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_H2x3y_S_Px_S_M1_vrr = QCX*I_ERI_H2x3y_S_S_S_M1_vrr+WQX*I_ERI_H2x3y_S_S_S_M2_vrr+2*oned2k*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Px_S_M1_vrr = QCX*I_ERI_H2x2yz_S_S_S_M1_vrr+WQX*I_ERI_H2x2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_Gx2yz_S_S_S_M2_vrr;
      Double I_ERI_H2xy2z_S_Px_S_M1_vrr = QCX*I_ERI_H2xy2z_S_S_S_M1_vrr+WQX*I_ERI_H2xy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gxy2z_S_S_S_M2_vrr;
      Double I_ERI_H2x3z_S_Px_S_M1_vrr = QCX*I_ERI_H2x3z_S_S_S_M1_vrr+WQX*I_ERI_H2x3z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4y_S_Px_S_M1_vrr = QCX*I_ERI_Hx4y_S_S_S_M1_vrr+WQX*I_ERI_Hx4y_S_S_S_M2_vrr+oned2k*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_Hx3yz_S_Px_S_M1_vrr = QCX*I_ERI_Hx3yz_S_S_S_M1_vrr+WQX*I_ERI_Hx3yz_S_S_S_M2_vrr+oned2k*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_Hx2y2z_S_Px_S_M1_vrr = QCX*I_ERI_Hx2y2z_S_S_S_M1_vrr+WQX*I_ERI_Hx2y2z_S_S_S_M2_vrr+oned2k*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_Hxy3z_S_Px_S_M1_vrr = QCX*I_ERI_Hxy3z_S_S_S_M1_vrr+WQX*I_ERI_Hxy3z_S_S_S_M2_vrr+oned2k*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4z_S_Px_S_M1_vrr = QCX*I_ERI_Hx4z_S_S_S_M1_vrr+WQX*I_ERI_Hx4z_S_S_S_M2_vrr+oned2k*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_H5y_S_Px_S_M1_vrr = QCX*I_ERI_H5y_S_S_S_M1_vrr+WQX*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_H4yz_S_Px_S_M1_vrr = QCX*I_ERI_H4yz_S_S_S_M1_vrr+WQX*I_ERI_H4yz_S_S_S_M2_vrr;
      Double I_ERI_H3y2z_S_Px_S_M1_vrr = QCX*I_ERI_H3y2z_S_S_S_M1_vrr+WQX*I_ERI_H3y2z_S_S_S_M2_vrr;
      Double I_ERI_H2y3z_S_Px_S_M1_vrr = QCX*I_ERI_H2y3z_S_S_S_M1_vrr+WQX*I_ERI_H2y3z_S_S_S_M2_vrr;
      Double I_ERI_Hy4z_S_Px_S_M1_vrr = QCX*I_ERI_Hy4z_S_S_S_M1_vrr+WQX*I_ERI_Hy4z_S_S_S_M2_vrr;
      Double I_ERI_H5z_S_Px_S_M1_vrr = QCX*I_ERI_H5z_S_S_S_M1_vrr+WQX*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_H5x_S_Py_S_M1_vrr = QCY*I_ERI_H5x_S_S_S_M1_vrr+WQY*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_H4xy_S_Py_S_M1_vrr = QCY*I_ERI_H4xy_S_S_S_M1_vrr+WQY*I_ERI_H4xy_S_S_S_M2_vrr+oned2k*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H4xz_S_Py_S_M1_vrr = QCY*I_ERI_H4xz_S_S_S_M1_vrr+WQY*I_ERI_H4xz_S_S_S_M2_vrr;
      Double I_ERI_H3x2y_S_Py_S_M1_vrr = QCY*I_ERI_H3x2y_S_S_S_M1_vrr+WQY*I_ERI_H3x2y_S_S_S_M2_vrr+2*oned2k*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_H3xyz_S_Py_S_M1_vrr = QCY*I_ERI_H3xyz_S_S_S_M1_vrr+WQY*I_ERI_H3xyz_S_S_S_M2_vrr+oned2k*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_H3x2z_S_Py_S_M1_vrr = QCY*I_ERI_H3x2z_S_S_S_M1_vrr+WQY*I_ERI_H3x2z_S_S_S_M2_vrr;
      Double I_ERI_H2x3y_S_Py_S_M1_vrr = QCY*I_ERI_H2x3y_S_S_S_M1_vrr+WQY*I_ERI_H2x3y_S_S_S_M2_vrr+3*oned2k*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Py_S_M1_vrr = QCY*I_ERI_H2x2yz_S_S_S_M1_vrr+WQY*I_ERI_H2x2yz_S_S_S_M2_vrr+2*oned2k*I_ERI_G2xyz_S_S_S_M2_vrr;
      Double I_ERI_H2xy2z_S_Py_S_M1_vrr = QCY*I_ERI_H2xy2z_S_S_S_M1_vrr+WQY*I_ERI_H2xy2z_S_S_S_M2_vrr+oned2k*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_H2x3z_S_Py_S_M1_vrr = QCY*I_ERI_H2x3z_S_S_S_M1_vrr+WQY*I_ERI_H2x3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4y_S_Py_S_M1_vrr = QCY*I_ERI_Hx4y_S_S_S_M1_vrr+WQY*I_ERI_Hx4y_S_S_S_M2_vrr+4*oned2k*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_Hx3yz_S_Py_S_M1_vrr = QCY*I_ERI_Hx3yz_S_S_S_M1_vrr+WQY*I_ERI_Hx3yz_S_S_S_M2_vrr+3*oned2k*I_ERI_Gx2yz_S_S_S_M2_vrr;
      Double I_ERI_Hx2y2z_S_Py_S_M1_vrr = QCY*I_ERI_Hx2y2z_S_S_S_M1_vrr+WQY*I_ERI_Hx2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gxy2z_S_S_S_M2_vrr;
      Double I_ERI_Hxy3z_S_Py_S_M1_vrr = QCY*I_ERI_Hxy3z_S_S_S_M1_vrr+WQY*I_ERI_Hxy3z_S_S_S_M2_vrr+oned2k*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_Hx4z_S_Py_S_M1_vrr = QCY*I_ERI_Hx4z_S_S_S_M1_vrr+WQY*I_ERI_Hx4z_S_S_S_M2_vrr;
      Double I_ERI_H5y_S_Py_S_M1_vrr = QCY*I_ERI_H5y_S_S_S_M1_vrr+WQY*I_ERI_H5y_S_S_S_M2_vrr+5*oned2k*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_H4yz_S_Py_S_M1_vrr = QCY*I_ERI_H4yz_S_S_S_M1_vrr+WQY*I_ERI_H4yz_S_S_S_M2_vrr+4*oned2k*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_H3y2z_S_Py_S_M1_vrr = QCY*I_ERI_H3y2z_S_S_S_M1_vrr+WQY*I_ERI_H3y2z_S_S_S_M2_vrr+3*oned2k*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_H2y3z_S_Py_S_M1_vrr = QCY*I_ERI_H2y3z_S_S_S_M1_vrr+WQY*I_ERI_H2y3z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_Hy4z_S_Py_S_M1_vrr = QCY*I_ERI_Hy4z_S_S_S_M1_vrr+WQY*I_ERI_Hy4z_S_S_S_M2_vrr+oned2k*I_ERI_G4z_S_S_S_M2_vrr;
      Double I_ERI_H5z_S_Py_S_M1_vrr = QCY*I_ERI_H5z_S_S_S_M1_vrr+WQY*I_ERI_H5z_S_S_S_M2_vrr;
      Double I_ERI_H5x_S_Pz_S_M1_vrr = QCZ*I_ERI_H5x_S_S_S_M1_vrr+WQZ*I_ERI_H5x_S_S_S_M2_vrr;
      Double I_ERI_H4xy_S_Pz_S_M1_vrr = QCZ*I_ERI_H4xy_S_S_S_M1_vrr+WQZ*I_ERI_H4xy_S_S_S_M2_vrr;
      Double I_ERI_H4xz_S_Pz_S_M1_vrr = QCZ*I_ERI_H4xz_S_S_S_M1_vrr+WQZ*I_ERI_H4xz_S_S_S_M2_vrr+oned2k*I_ERI_G4x_S_S_S_M2_vrr;
      Double I_ERI_H3x2y_S_Pz_S_M1_vrr = QCZ*I_ERI_H3x2y_S_S_S_M1_vrr+WQZ*I_ERI_H3x2y_S_S_S_M2_vrr;
      Double I_ERI_H3xyz_S_Pz_S_M1_vrr = QCZ*I_ERI_H3xyz_S_S_S_M1_vrr+WQZ*I_ERI_H3xyz_S_S_S_M2_vrr+oned2k*I_ERI_G3xy_S_S_S_M2_vrr;
      Double I_ERI_H3x2z_S_Pz_S_M1_vrr = QCZ*I_ERI_H3x2z_S_S_S_M1_vrr+WQZ*I_ERI_H3x2z_S_S_S_M2_vrr+2*oned2k*I_ERI_G3xz_S_S_S_M2_vrr;
      Double I_ERI_H2x3y_S_Pz_S_M1_vrr = QCZ*I_ERI_H2x3y_S_S_S_M1_vrr+WQZ*I_ERI_H2x3y_S_S_S_M2_vrr;
      Double I_ERI_H2x2yz_S_Pz_S_M1_vrr = QCZ*I_ERI_H2x2yz_S_S_S_M1_vrr+WQZ*I_ERI_H2x2yz_S_S_S_M2_vrr+oned2k*I_ERI_G2x2y_S_S_S_M2_vrr;
      Double I_ERI_H2xy2z_S_Pz_S_M1_vrr = QCZ*I_ERI_H2xy2z_S_S_S_M1_vrr+WQZ*I_ERI_H2xy2z_S_S_S_M2_vrr+2*oned2k*I_ERI_G2xyz_S_S_S_M2_vrr;
      Double I_ERI_H2x3z_S_Pz_S_M1_vrr = QCZ*I_ERI_H2x3z_S_S_S_M1_vrr+WQZ*I_ERI_H2x3z_S_S_S_M2_vrr+3*oned2k*I_ERI_G2x2z_S_S_S_M2_vrr;
      Double I_ERI_Hx4y_S_Pz_S_M1_vrr = QCZ*I_ERI_Hx4y_S_S_S_M1_vrr+WQZ*I_ERI_Hx4y_S_S_S_M2_vrr;
      Double I_ERI_Hx3yz_S_Pz_S_M1_vrr = QCZ*I_ERI_Hx3yz_S_S_S_M1_vrr+WQZ*I_ERI_Hx3yz_S_S_S_M2_vrr+oned2k*I_ERI_Gx3y_S_S_S_M2_vrr;
      Double I_ERI_Hx2y2z_S_Pz_S_M1_vrr = QCZ*I_ERI_Hx2y2z_S_S_S_M1_vrr+WQZ*I_ERI_Hx2y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_Gx2yz_S_S_S_M2_vrr;
      Double I_ERI_Hxy3z_S_Pz_S_M1_vrr = QCZ*I_ERI_Hxy3z_S_S_S_M1_vrr+WQZ*I_ERI_Hxy3z_S_S_S_M2_vrr+3*oned2k*I_ERI_Gxy2z_S_S_S_M2_vrr;
      Double I_ERI_Hx4z_S_Pz_S_M1_vrr = QCZ*I_ERI_Hx4z_S_S_S_M1_vrr+WQZ*I_ERI_Hx4z_S_S_S_M2_vrr+4*oned2k*I_ERI_Gx3z_S_S_S_M2_vrr;
      Double I_ERI_H5y_S_Pz_S_M1_vrr = QCZ*I_ERI_H5y_S_S_S_M1_vrr+WQZ*I_ERI_H5y_S_S_S_M2_vrr;
      Double I_ERI_H4yz_S_Pz_S_M1_vrr = QCZ*I_ERI_H4yz_S_S_S_M1_vrr+WQZ*I_ERI_H4yz_S_S_S_M2_vrr+oned2k*I_ERI_G4y_S_S_S_M2_vrr;
      Double I_ERI_H3y2z_S_Pz_S_M1_vrr = QCZ*I_ERI_H3y2z_S_S_S_M1_vrr+WQZ*I_ERI_H3y2z_S_S_S_M2_vrr+2*oned2k*I_ERI_G3yz_S_S_S_M2_vrr;
      Double I_ERI_H2y3z_S_Pz_S_M1_vrr = QCZ*I_ERI_H2y3z_S_S_S_M1_vrr+WQZ*I_ERI_H2y3z_S_S_S_M2_vrr+3*oned2k*I_ERI_G2y2z_S_S_S_M2_vrr;
      Double I_ERI_Hy4z_S_Pz_S_M1_vrr = QCZ*I_ERI_Hy4z_S_S_S_M1_vrr+WQZ*I_ERI_Hy4z_S_S_S_M2_vrr+4*oned2k*I_ERI_Gy3z_S_S_S_M2_vrr;
      Double I_ERI_H5z_S_Pz_S_M1_vrr = QCZ*I_ERI_H5z_S_S_S_M1_vrr+WQZ*I_ERI_H5z_S_S_S_M2_vrr+5*oned2k*I_ERI_G4z_S_S_S_M2_vrr;

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
       * shell quartet name: SQ_ERI_H_S_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_G_S_P_S
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_F_S_P_S
       * RHS shell quartet name: SQ_ERI_F_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_S_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_Px_S_vrr = PAX*I_ERI_G4x_S_Px_S_vrr+WPX*I_ERI_G4x_S_Px_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Px_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Px_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_Px_S_vrr = PAY*I_ERI_G4x_S_Px_S_vrr+WPY*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_Px_S_vrr = PAZ*I_ERI_G4x_S_Px_S_vrr+WPZ*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_Px_S_vrr = PAY*I_ERI_G3xy_S_Px_S_vrr+WPY*I_ERI_G3xy_S_Px_S_M1_vrr+oned2z*I_ERI_F3x_S_Px_S_vrr-rhod2zsq*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_Px_S_vrr = PAZ*I_ERI_G3xy_S_Px_S_vrr+WPZ*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_Px_S_vrr = PAZ*I_ERI_G3xz_S_Px_S_vrr+WPZ*I_ERI_G3xz_S_Px_S_M1_vrr+oned2z*I_ERI_F3x_S_Px_S_vrr-rhod2zsq*I_ERI_F3x_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_Px_S_vrr = PAX*I_ERI_Gx3y_S_Px_S_vrr+WPX*I_ERI_Gx3y_S_Px_S_M1_vrr+oned2z*I_ERI_F3y_S_Px_S_vrr-rhod2zsq*I_ERI_F3y_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Px_S_vrr = PAZ*I_ERI_G2x2y_S_Px_S_vrr+WPZ*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Px_S_vrr = PAY*I_ERI_G2x2z_S_Px_S_vrr+WPY*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_Px_S_vrr = PAX*I_ERI_Gx3z_S_Px_S_vrr+WPX*I_ERI_Gx3z_S_Px_S_M1_vrr+oned2z*I_ERI_F3z_S_Px_S_vrr-rhod2zsq*I_ERI_F3z_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_Px_S_vrr = PAX*I_ERI_G4y_S_Px_S_vrr+WPX*I_ERI_G4y_S_Px_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Px_S_vrr = PAZ*I_ERI_Gx3y_S_Px_S_vrr+WPZ*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Px_S_vrr = PAX*I_ERI_G2y2z_S_Px_S_vrr+WPX*I_ERI_G2y2z_S_Px_S_M1_vrr+oned2k*I_ERI_G2y2z_S_S_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Px_S_vrr = PAY*I_ERI_Gx3z_S_Px_S_vrr+WPY*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_Px_S_vrr = PAX*I_ERI_G4z_S_Px_S_vrr+WPX*I_ERI_G4z_S_Px_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_Px_S_vrr = PAY*I_ERI_G4y_S_Px_S_vrr+WPY*I_ERI_G4y_S_Px_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Px_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_H4yz_S_Px_S_vrr = PAZ*I_ERI_G4y_S_Px_S_vrr+WPZ*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_H3y2z_S_Px_S_vrr = PAZ*I_ERI_G3yz_S_Px_S_vrr+WPZ*I_ERI_G3yz_S_Px_S_M1_vrr+oned2z*I_ERI_F3y_S_Px_S_vrr-rhod2zsq*I_ERI_F3y_S_Px_S_M1_vrr;
      Double I_ERI_H2y3z_S_Px_S_vrr = PAY*I_ERI_Gy3z_S_Px_S_vrr+WPY*I_ERI_Gy3z_S_Px_S_M1_vrr+oned2z*I_ERI_F3z_S_Px_S_vrr-rhod2zsq*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_Hy4z_S_Px_S_vrr = PAY*I_ERI_G4z_S_Px_S_vrr+WPY*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5z_S_Px_S_vrr = PAZ*I_ERI_G4z_S_Px_S_vrr+WPZ*I_ERI_G4z_S_Px_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Px_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Px_S_M1_vrr;
      Double I_ERI_H5x_S_Py_S_vrr = PAX*I_ERI_G4x_S_Py_S_vrr+WPX*I_ERI_G4x_S_Py_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Py_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_H4xy_S_Py_S_vrr = PAY*I_ERI_G4x_S_Py_S_vrr+WPY*I_ERI_G4x_S_Py_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_Py_S_vrr = PAZ*I_ERI_G4x_S_Py_S_vrr+WPZ*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H3x2y_S_Py_S_vrr = PAY*I_ERI_G3xy_S_Py_S_vrr+WPY*I_ERI_G3xy_S_Py_S_M1_vrr+oned2z*I_ERI_F3x_S_Py_S_vrr-rhod2zsq*I_ERI_F3x_S_Py_S_M1_vrr+oned2k*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_Py_S_vrr = PAZ*I_ERI_G3xy_S_Py_S_vrr+WPZ*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_H3x2z_S_Py_S_vrr = PAZ*I_ERI_G3xz_S_Py_S_vrr+WPZ*I_ERI_G3xz_S_Py_S_M1_vrr+oned2z*I_ERI_F3x_S_Py_S_vrr-rhod2zsq*I_ERI_F3x_S_Py_S_M1_vrr;
      Double I_ERI_H2x3y_S_Py_S_vrr = PAX*I_ERI_Gx3y_S_Py_S_vrr+WPX*I_ERI_Gx3y_S_Py_S_M1_vrr+oned2z*I_ERI_F3y_S_Py_S_vrr-rhod2zsq*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Py_S_vrr = PAZ*I_ERI_G2x2y_S_Py_S_vrr+WPZ*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Py_S_vrr = PAY*I_ERI_G2x2z_S_Py_S_vrr+WPY*I_ERI_G2x2z_S_Py_S_M1_vrr+oned2k*I_ERI_G2x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3z_S_Py_S_vrr = PAX*I_ERI_Gx3z_S_Py_S_vrr+WPX*I_ERI_Gx3z_S_Py_S_M1_vrr+oned2z*I_ERI_F3z_S_Py_S_vrr-rhod2zsq*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4y_S_Py_S_vrr = PAX*I_ERI_G4y_S_Py_S_vrr+WPX*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Py_S_vrr = PAZ*I_ERI_Gx3y_S_Py_S_vrr+WPZ*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Py_S_vrr = PAX*I_ERI_G2y2z_S_Py_S_vrr+WPX*I_ERI_G2y2z_S_Py_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Py_S_vrr = PAY*I_ERI_Gx3z_S_Py_S_vrr+WPY*I_ERI_Gx3z_S_Py_S_M1_vrr+oned2k*I_ERI_Gx3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4z_S_Py_S_vrr = PAX*I_ERI_G4z_S_Py_S_vrr+WPX*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5y_S_Py_S_vrr = PAY*I_ERI_G4y_S_Py_S_vrr+WPY*I_ERI_G4y_S_Py_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Py_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Py_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_Py_S_vrr = PAZ*I_ERI_G4y_S_Py_S_vrr+WPZ*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_H3y2z_S_Py_S_vrr = PAZ*I_ERI_G3yz_S_Py_S_vrr+WPZ*I_ERI_G3yz_S_Py_S_M1_vrr+oned2z*I_ERI_F3y_S_Py_S_vrr-rhod2zsq*I_ERI_F3y_S_Py_S_M1_vrr;
      Double I_ERI_H2y3z_S_Py_S_vrr = PAY*I_ERI_Gy3z_S_Py_S_vrr+WPY*I_ERI_Gy3z_S_Py_S_M1_vrr+oned2z*I_ERI_F3z_S_Py_S_vrr-rhod2zsq*I_ERI_F3z_S_Py_S_M1_vrr+oned2k*I_ERI_Gy3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_Py_S_vrr = PAY*I_ERI_G4z_S_Py_S_vrr+WPY*I_ERI_G4z_S_Py_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_Py_S_vrr = PAZ*I_ERI_G4z_S_Py_S_vrr+WPZ*I_ERI_G4z_S_Py_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Py_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Py_S_M1_vrr;
      Double I_ERI_H5x_S_Pz_S_vrr = PAX*I_ERI_G4x_S_Pz_S_vrr+WPX*I_ERI_G4x_S_Pz_S_M1_vrr+4*oned2z*I_ERI_F3x_S_Pz_S_vrr-4*rhod2zsq*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_H4xy_S_Pz_S_vrr = PAY*I_ERI_G4x_S_Pz_S_vrr+WPY*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_H4xz_S_Pz_S_vrr = PAZ*I_ERI_G4x_S_Pz_S_vrr+WPZ*I_ERI_G4x_S_Pz_S_M1_vrr+oned2k*I_ERI_G4x_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_Pz_S_vrr = PAY*I_ERI_G3xy_S_Pz_S_vrr+WPY*I_ERI_G3xy_S_Pz_S_M1_vrr+oned2z*I_ERI_F3x_S_Pz_S_vrr-rhod2zsq*I_ERI_F3x_S_Pz_S_M1_vrr;
      Double I_ERI_H3xyz_S_Pz_S_vrr = PAZ*I_ERI_G3xy_S_Pz_S_vrr+WPZ*I_ERI_G3xy_S_Pz_S_M1_vrr+oned2k*I_ERI_G3xy_S_S_S_M1_vrr;
      Double I_ERI_H3x2z_S_Pz_S_vrr = PAZ*I_ERI_G3xz_S_Pz_S_vrr+WPZ*I_ERI_G3xz_S_Pz_S_M1_vrr+oned2z*I_ERI_F3x_S_Pz_S_vrr-rhod2zsq*I_ERI_F3x_S_Pz_S_M1_vrr+oned2k*I_ERI_G3xz_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_Pz_S_vrr = PAX*I_ERI_Gx3y_S_Pz_S_vrr+WPX*I_ERI_Gx3y_S_Pz_S_M1_vrr+oned2z*I_ERI_F3y_S_Pz_S_vrr-rhod2zsq*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Pz_S_vrr = PAZ*I_ERI_G2x2y_S_Pz_S_vrr+WPZ*I_ERI_G2x2y_S_Pz_S_M1_vrr+oned2k*I_ERI_G2x2y_S_S_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Pz_S_vrr = PAY*I_ERI_G2x2z_S_Pz_S_vrr+WPY*I_ERI_G2x2z_S_Pz_S_M1_vrr;
      Double I_ERI_H2x3z_S_Pz_S_vrr = PAX*I_ERI_Gx3z_S_Pz_S_vrr+WPX*I_ERI_Gx3z_S_Pz_S_M1_vrr+oned2z*I_ERI_F3z_S_Pz_S_vrr-rhod2zsq*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4y_S_Pz_S_vrr = PAX*I_ERI_G4y_S_Pz_S_vrr+WPX*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Pz_S_vrr = PAZ*I_ERI_Gx3y_S_Pz_S_vrr+WPZ*I_ERI_Gx3y_S_Pz_S_M1_vrr+oned2k*I_ERI_Gx3y_S_S_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Pz_S_vrr = PAX*I_ERI_G2y2z_S_Pz_S_vrr+WPX*I_ERI_G2y2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Pz_S_vrr = PAY*I_ERI_Gx3z_S_Pz_S_vrr+WPY*I_ERI_Gx3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4z_S_Pz_S_vrr = PAX*I_ERI_G4z_S_Pz_S_vrr+WPX*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_H5y_S_Pz_S_vrr = PAY*I_ERI_G4y_S_Pz_S_vrr+WPY*I_ERI_G4y_S_Pz_S_M1_vrr+4*oned2z*I_ERI_F3y_S_Pz_S_vrr-4*rhod2zsq*I_ERI_F3y_S_Pz_S_M1_vrr;
      Double I_ERI_H4yz_S_Pz_S_vrr = PAZ*I_ERI_G4y_S_Pz_S_vrr+WPZ*I_ERI_G4y_S_Pz_S_M1_vrr+oned2k*I_ERI_G4y_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_Pz_S_vrr = PAZ*I_ERI_G3yz_S_Pz_S_vrr+WPZ*I_ERI_G3yz_S_Pz_S_M1_vrr+oned2z*I_ERI_F3y_S_Pz_S_vrr-rhod2zsq*I_ERI_F3y_S_Pz_S_M1_vrr+oned2k*I_ERI_G3yz_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_Pz_S_vrr = PAY*I_ERI_Gy3z_S_Pz_S_vrr+WPY*I_ERI_Gy3z_S_Pz_S_M1_vrr+oned2z*I_ERI_F3z_S_Pz_S_vrr-rhod2zsq*I_ERI_F3z_S_Pz_S_M1_vrr;
      Double I_ERI_Hy4z_S_Pz_S_vrr = PAY*I_ERI_G4z_S_Pz_S_vrr+WPY*I_ERI_G4z_S_Pz_S_M1_vrr;
      Double I_ERI_H5z_S_Pz_S_vrr = PAZ*I_ERI_G4z_S_Pz_S_vrr+WPZ*I_ERI_G4z_S_Pz_S_M1_vrr+4*oned2z*I_ERI_F3z_S_Pz_S_vrr-4*rhod2zsq*I_ERI_F3z_S_Pz_S_M1_vrr+oned2k*I_ERI_G4z_S_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_D_S
       * expanding position: KET1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_ERI_H_S_P_S
       * RHS shell quartet name: SQ_ERI_H_S_P_S_M1
       * RHS shell quartet name: SQ_ERI_H_S_S_S
       * RHS shell quartet name: SQ_ERI_H_S_S_S_M1
       * RHS shell quartet name: SQ_ERI_G_S_P_S_M1
       ************************************************************/
      Double I_ERI_H5x_S_D2x_S_vrr = QCX*I_ERI_H5x_S_Px_S_vrr+WQX*I_ERI_H5x_S_Px_S_M1_vrr+oned2e*I_ERI_H5x_S_S_S_vrr-rhod2esq*I_ERI_H5x_S_S_S_M1_vrr+5*oned2k*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H4xy_S_D2x_S_vrr = QCX*I_ERI_H4xy_S_Px_S_vrr+WQX*I_ERI_H4xy_S_Px_S_M1_vrr+oned2e*I_ERI_H4xy_S_S_S_vrr-rhod2esq*I_ERI_H4xy_S_S_S_M1_vrr+4*oned2k*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_D2x_S_vrr = QCX*I_ERI_H4xz_S_Px_S_vrr+WQX*I_ERI_H4xz_S_Px_S_M1_vrr+oned2e*I_ERI_H4xz_S_S_S_vrr-rhod2esq*I_ERI_H4xz_S_S_S_M1_vrr+4*oned2k*I_ERI_G3xz_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_D2x_S_vrr = QCX*I_ERI_H3x2y_S_Px_S_vrr+WQX*I_ERI_H3x2y_S_Px_S_M1_vrr+oned2e*I_ERI_H3x2y_S_S_S_vrr-rhod2esq*I_ERI_H3x2y_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_D2x_S_vrr = QCX*I_ERI_H3xyz_S_Px_S_vrr+WQX*I_ERI_H3xyz_S_Px_S_M1_vrr+oned2e*I_ERI_H3xyz_S_S_S_vrr-rhod2esq*I_ERI_H3xyz_S_S_S_M1_vrr+3*oned2k*I_ERI_G2xyz_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_D2x_S_vrr = QCX*I_ERI_H3x2z_S_Px_S_vrr+WQX*I_ERI_H3x2z_S_Px_S_M1_vrr+oned2e*I_ERI_H3x2z_S_S_S_vrr-rhod2esq*I_ERI_H3x2z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_D2x_S_vrr = QCX*I_ERI_H2x3y_S_Px_S_vrr+WQX*I_ERI_H2x3y_S_Px_S_M1_vrr+oned2e*I_ERI_H2x3y_S_S_S_vrr-rhod2esq*I_ERI_H2x3y_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_H2x2yz_S_D2x_S_vrr = QCX*I_ERI_H2x2yz_S_Px_S_vrr+WQX*I_ERI_H2x2yz_S_Px_S_M1_vrr+oned2e*I_ERI_H2x2yz_S_S_S_vrr-rhod2esq*I_ERI_H2x2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_D2x_S_vrr = QCX*I_ERI_H2xy2z_S_Px_S_vrr+WQX*I_ERI_H2xy2z_S_Px_S_M1_vrr+oned2e*I_ERI_H2xy2z_S_S_S_vrr-rhod2esq*I_ERI_H2xy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gxy2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_D2x_S_vrr = QCX*I_ERI_H2x3z_S_Px_S_vrr+WQX*I_ERI_H2x3z_S_Px_S_M1_vrr+oned2e*I_ERI_H2x3z_S_S_S_vrr-rhod2esq*I_ERI_H2x3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4y_S_D2x_S_vrr = QCX*I_ERI_Hx4y_S_Px_S_vrr+WQX*I_ERI_Hx4y_S_Px_S_M1_vrr+oned2e*I_ERI_Hx4y_S_S_S_vrr-rhod2esq*I_ERI_Hx4y_S_S_S_M1_vrr+oned2k*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_Hx3yz_S_D2x_S_vrr = QCX*I_ERI_Hx3yz_S_Px_S_vrr+WQX*I_ERI_Hx3yz_S_Px_S_M1_vrr+oned2e*I_ERI_Hx3yz_S_S_S_vrr-rhod2esq*I_ERI_Hx3yz_S_S_S_M1_vrr+oned2k*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_D2x_S_vrr = QCX*I_ERI_Hx2y2z_S_Px_S_vrr+WQX*I_ERI_Hx2y2z_S_Px_S_M1_vrr+oned2e*I_ERI_Hx2y2z_S_S_S_vrr-rhod2esq*I_ERI_Hx2y2z_S_S_S_M1_vrr+oned2k*I_ERI_G2y2z_S_Px_S_M1_vrr;
      Double I_ERI_Hxy3z_S_D2x_S_vrr = QCX*I_ERI_Hxy3z_S_Px_S_vrr+WQX*I_ERI_Hxy3z_S_Px_S_M1_vrr+oned2e*I_ERI_Hxy3z_S_S_S_vrr-rhod2esq*I_ERI_Hxy3z_S_S_S_M1_vrr+oned2k*I_ERI_Gy3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_D2x_S_vrr = QCX*I_ERI_Hx4z_S_Px_S_vrr+WQX*I_ERI_Hx4z_S_Px_S_M1_vrr+oned2e*I_ERI_Hx4z_S_S_S_vrr-rhod2esq*I_ERI_Hx4z_S_S_S_M1_vrr+oned2k*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5y_S_D2x_S_vrr = QCX*I_ERI_H5y_S_Px_S_vrr+WQX*I_ERI_H5y_S_Px_S_M1_vrr+oned2e*I_ERI_H5y_S_S_S_vrr-rhod2esq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_D2x_S_vrr = QCX*I_ERI_H4yz_S_Px_S_vrr+WQX*I_ERI_H4yz_S_Px_S_M1_vrr+oned2e*I_ERI_H4yz_S_S_S_vrr-rhod2esq*I_ERI_H4yz_S_S_S_M1_vrr;
      Double I_ERI_H3y2z_S_D2x_S_vrr = QCX*I_ERI_H3y2z_S_Px_S_vrr+WQX*I_ERI_H3y2z_S_Px_S_M1_vrr+oned2e*I_ERI_H3y2z_S_S_S_vrr-rhod2esq*I_ERI_H3y2z_S_S_S_M1_vrr;
      Double I_ERI_H2y3z_S_D2x_S_vrr = QCX*I_ERI_H2y3z_S_Px_S_vrr+WQX*I_ERI_H2y3z_S_Px_S_M1_vrr+oned2e*I_ERI_H2y3z_S_S_S_vrr-rhod2esq*I_ERI_H2y3z_S_S_S_M1_vrr;
      Double I_ERI_Hy4z_S_D2x_S_vrr = QCX*I_ERI_Hy4z_S_Px_S_vrr+WQX*I_ERI_Hy4z_S_Px_S_M1_vrr+oned2e*I_ERI_Hy4z_S_S_S_vrr-rhod2esq*I_ERI_Hy4z_S_S_S_M1_vrr;
      Double I_ERI_H5z_S_D2x_S_vrr = QCX*I_ERI_H5z_S_Px_S_vrr+WQX*I_ERI_H5z_S_Px_S_M1_vrr+oned2e*I_ERI_H5z_S_S_S_vrr-rhod2esq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_H5x_S_Dxy_S_vrr = QCY*I_ERI_H5x_S_Px_S_vrr+WQY*I_ERI_H5x_S_Px_S_M1_vrr;
      Double I_ERI_H4xy_S_Dxy_S_vrr = QCY*I_ERI_H4xy_S_Px_S_vrr+WQY*I_ERI_H4xy_S_Px_S_M1_vrr+oned2k*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_Dxy_S_vrr = QCY*I_ERI_H4xz_S_Px_S_vrr+WQY*I_ERI_H4xz_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_Dxy_S_vrr = QCY*I_ERI_H3x2y_S_Px_S_vrr+WQY*I_ERI_H3x2y_S_Px_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_Dxy_S_vrr = QCY*I_ERI_H3xyz_S_Px_S_vrr+WQY*I_ERI_H3xyz_S_Px_S_M1_vrr+oned2k*I_ERI_G3xz_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_Dxy_S_vrr = QCY*I_ERI_H3x2z_S_Px_S_vrr+WQY*I_ERI_H3x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_Dxy_S_vrr = QCY*I_ERI_H2x3y_S_Px_S_vrr+WQY*I_ERI_H2x3y_S_Px_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Dxy_S_vrr = QCY*I_ERI_H2x2yz_S_Px_S_vrr+WQY*I_ERI_H2x2yz_S_Px_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Dxy_S_vrr = QCY*I_ERI_H2xy2z_S_Px_S_vrr+WQY*I_ERI_H2xy2z_S_Px_S_M1_vrr+oned2k*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_Dxy_S_vrr = QCY*I_ERI_H2x3z_S_Px_S_vrr+WQY*I_ERI_H2x3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4y_S_Dxy_S_vrr = QCY*I_ERI_Hx4y_S_Px_S_vrr+WQY*I_ERI_Hx4y_S_Px_S_M1_vrr+4*oned2k*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Dxy_S_vrr = QCY*I_ERI_Hx3yz_S_Px_S_vrr+WQY*I_ERI_Hx3yz_S_Px_S_M1_vrr+3*oned2k*I_ERI_Gx2yz_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Dxy_S_vrr = QCY*I_ERI_Hx2y2z_S_Px_S_vrr+WQY*I_ERI_Hx2y2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Gxy2z_S_Px_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Dxy_S_vrr = QCY*I_ERI_Hxy3z_S_Px_S_vrr+WQY*I_ERI_Hxy3z_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_Dxy_S_vrr = QCY*I_ERI_Hx4z_S_Px_S_vrr+WQY*I_ERI_Hx4z_S_Px_S_M1_vrr;
      Double I_ERI_H5y_S_Dxy_S_vrr = QCY*I_ERI_H5y_S_Px_S_vrr+WQY*I_ERI_H5y_S_Px_S_M1_vrr+5*oned2k*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_H4yz_S_Dxy_S_vrr = QCY*I_ERI_H4yz_S_Px_S_vrr+WQY*I_ERI_H4yz_S_Px_S_M1_vrr+4*oned2k*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_H3y2z_S_Dxy_S_vrr = QCY*I_ERI_H3y2z_S_Px_S_vrr+WQY*I_ERI_H3y2z_S_Px_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Px_S_M1_vrr;
      Double I_ERI_H2y3z_S_Dxy_S_vrr = QCY*I_ERI_H2y3z_S_Px_S_vrr+WQY*I_ERI_H2y3z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Gy3z_S_Px_S_M1_vrr;
      Double I_ERI_Hy4z_S_Dxy_S_vrr = QCY*I_ERI_Hy4z_S_Px_S_vrr+WQY*I_ERI_Hy4z_S_Px_S_M1_vrr+oned2k*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5z_S_Dxy_S_vrr = QCY*I_ERI_H5z_S_Px_S_vrr+WQY*I_ERI_H5z_S_Px_S_M1_vrr;
      Double I_ERI_H5x_S_Dxz_S_vrr = QCZ*I_ERI_H5x_S_Px_S_vrr+WQZ*I_ERI_H5x_S_Px_S_M1_vrr;
      Double I_ERI_H4xy_S_Dxz_S_vrr = QCZ*I_ERI_H4xy_S_Px_S_vrr+WQZ*I_ERI_H4xy_S_Px_S_M1_vrr;
      Double I_ERI_H4xz_S_Dxz_S_vrr = QCZ*I_ERI_H4xz_S_Px_S_vrr+WQZ*I_ERI_H4xz_S_Px_S_M1_vrr+oned2k*I_ERI_G4x_S_Px_S_M1_vrr;
      Double I_ERI_H3x2y_S_Dxz_S_vrr = QCZ*I_ERI_H3x2y_S_Px_S_vrr+WQZ*I_ERI_H3x2y_S_Px_S_M1_vrr;
      Double I_ERI_H3xyz_S_Dxz_S_vrr = QCZ*I_ERI_H3xyz_S_Px_S_vrr+WQZ*I_ERI_H3xyz_S_Px_S_M1_vrr+oned2k*I_ERI_G3xy_S_Px_S_M1_vrr;
      Double I_ERI_H3x2z_S_Dxz_S_vrr = QCZ*I_ERI_H3x2z_S_Px_S_vrr+WQZ*I_ERI_H3x2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_Px_S_M1_vrr;
      Double I_ERI_H2x3y_S_Dxz_S_vrr = QCZ*I_ERI_H2x3y_S_Px_S_vrr+WQZ*I_ERI_H2x3y_S_Px_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Dxz_S_vrr = QCZ*I_ERI_H2x2yz_S_Px_S_vrr+WQZ*I_ERI_H2x2yz_S_Px_S_M1_vrr+oned2k*I_ERI_G2x2y_S_Px_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Dxz_S_vrr = QCZ*I_ERI_H2xy2z_S_Px_S_vrr+WQZ*I_ERI_H2xy2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Px_S_M1_vrr;
      Double I_ERI_H2x3z_S_Dxz_S_vrr = QCZ*I_ERI_H2x3z_S_Px_S_vrr+WQZ*I_ERI_H2x3z_S_Px_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4y_S_Dxz_S_vrr = QCZ*I_ERI_Hx4y_S_Px_S_vrr+WQZ*I_ERI_Hx4y_S_Px_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Dxz_S_vrr = QCZ*I_ERI_Hx3yz_S_Px_S_vrr+WQZ*I_ERI_Hx3yz_S_Px_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Px_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Dxz_S_vrr = QCZ*I_ERI_Hx2y2z_S_Px_S_vrr+WQZ*I_ERI_Hx2y2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_Px_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Dxz_S_vrr = QCZ*I_ERI_Hxy3z_S_Px_S_vrr+WQZ*I_ERI_Hxy3z_S_Px_S_M1_vrr+3*oned2k*I_ERI_Gxy2z_S_Px_S_M1_vrr;
      Double I_ERI_Hx4z_S_Dxz_S_vrr = QCZ*I_ERI_Hx4z_S_Px_S_vrr+WQZ*I_ERI_Hx4z_S_Px_S_M1_vrr+4*oned2k*I_ERI_Gx3z_S_Px_S_M1_vrr;
      Double I_ERI_H5y_S_Dxz_S_vrr = QCZ*I_ERI_H5y_S_Px_S_vrr+WQZ*I_ERI_H5y_S_Px_S_M1_vrr;
      Double I_ERI_H4yz_S_Dxz_S_vrr = QCZ*I_ERI_H4yz_S_Px_S_vrr+WQZ*I_ERI_H4yz_S_Px_S_M1_vrr+oned2k*I_ERI_G4y_S_Px_S_M1_vrr;
      Double I_ERI_H3y2z_S_Dxz_S_vrr = QCZ*I_ERI_H3y2z_S_Px_S_vrr+WQZ*I_ERI_H3y2z_S_Px_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_Px_S_M1_vrr;
      Double I_ERI_H2y3z_S_Dxz_S_vrr = QCZ*I_ERI_H2y3z_S_Px_S_vrr+WQZ*I_ERI_H2y3z_S_Px_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Px_S_M1_vrr;
      Double I_ERI_Hy4z_S_Dxz_S_vrr = QCZ*I_ERI_Hy4z_S_Px_S_vrr+WQZ*I_ERI_Hy4z_S_Px_S_M1_vrr+4*oned2k*I_ERI_Gy3z_S_Px_S_M1_vrr;
      Double I_ERI_H5z_S_Dxz_S_vrr = QCZ*I_ERI_H5z_S_Px_S_vrr+WQZ*I_ERI_H5z_S_Px_S_M1_vrr+5*oned2k*I_ERI_G4z_S_Px_S_M1_vrr;
      Double I_ERI_H5x_S_D2y_S_vrr = QCY*I_ERI_H5x_S_Py_S_vrr+WQY*I_ERI_H5x_S_Py_S_M1_vrr+oned2e*I_ERI_H5x_S_S_S_vrr-rhod2esq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_D2y_S_vrr = QCY*I_ERI_H4xy_S_Py_S_vrr+WQY*I_ERI_H4xy_S_Py_S_M1_vrr+oned2e*I_ERI_H4xy_S_S_S_vrr-rhod2esq*I_ERI_H4xy_S_S_S_M1_vrr+oned2k*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H4xz_S_D2y_S_vrr = QCY*I_ERI_H4xz_S_Py_S_vrr+WQY*I_ERI_H4xz_S_Py_S_M1_vrr+oned2e*I_ERI_H4xz_S_S_S_vrr-rhod2esq*I_ERI_H4xz_S_S_S_M1_vrr;
      Double I_ERI_H3x2y_S_D2y_S_vrr = QCY*I_ERI_H3x2y_S_Py_S_vrr+WQY*I_ERI_H3x2y_S_Py_S_M1_vrr+oned2e*I_ERI_H3x2y_S_S_S_vrr-rhod2esq*I_ERI_H3x2y_S_S_S_M1_vrr+2*oned2k*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_H3xyz_S_D2y_S_vrr = QCY*I_ERI_H3xyz_S_Py_S_vrr+WQY*I_ERI_H3xyz_S_Py_S_M1_vrr+oned2e*I_ERI_H3xyz_S_S_S_vrr-rhod2esq*I_ERI_H3xyz_S_S_S_M1_vrr+oned2k*I_ERI_G3xz_S_Py_S_M1_vrr;
      Double I_ERI_H3x2z_S_D2y_S_vrr = QCY*I_ERI_H3x2z_S_Py_S_vrr+WQY*I_ERI_H3x2z_S_Py_S_M1_vrr+oned2e*I_ERI_H3x2z_S_S_S_vrr-rhod2esq*I_ERI_H3x2z_S_S_S_M1_vrr;
      Double I_ERI_H2x3y_S_D2y_S_vrr = QCY*I_ERI_H2x3y_S_Py_S_vrr+WQY*I_ERI_H2x3y_S_Py_S_M1_vrr+oned2e*I_ERI_H2x3y_S_S_S_vrr-rhod2esq*I_ERI_H2x3y_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_H2x2yz_S_D2y_S_vrr = QCY*I_ERI_H2x2yz_S_Py_S_vrr+WQY*I_ERI_H2x2yz_S_Py_S_M1_vrr+oned2e*I_ERI_H2x2yz_S_S_S_vrr-rhod2esq*I_ERI_H2x2yz_S_S_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Py_S_M1_vrr;
      Double I_ERI_H2xy2z_S_D2y_S_vrr = QCY*I_ERI_H2xy2z_S_Py_S_vrr+WQY*I_ERI_H2xy2z_S_Py_S_M1_vrr+oned2e*I_ERI_H2xy2z_S_S_S_vrr-rhod2esq*I_ERI_H2xy2z_S_S_S_M1_vrr+oned2k*I_ERI_G2x2z_S_Py_S_M1_vrr;
      Double I_ERI_H2x3z_S_D2y_S_vrr = QCY*I_ERI_H2x3z_S_Py_S_vrr+WQY*I_ERI_H2x3z_S_Py_S_M1_vrr+oned2e*I_ERI_H2x3z_S_S_S_vrr-rhod2esq*I_ERI_H2x3z_S_S_S_M1_vrr;
      Double I_ERI_Hx4y_S_D2y_S_vrr = QCY*I_ERI_Hx4y_S_Py_S_vrr+WQY*I_ERI_Hx4y_S_Py_S_M1_vrr+oned2e*I_ERI_Hx4y_S_S_S_vrr-rhod2esq*I_ERI_Hx4y_S_S_S_M1_vrr+4*oned2k*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_Hx3yz_S_D2y_S_vrr = QCY*I_ERI_Hx3yz_S_Py_S_vrr+WQY*I_ERI_Hx3yz_S_Py_S_M1_vrr+oned2e*I_ERI_Hx3yz_S_S_S_vrr-rhod2esq*I_ERI_Hx3yz_S_S_S_M1_vrr+3*oned2k*I_ERI_Gx2yz_S_Py_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_D2y_S_vrr = QCY*I_ERI_Hx2y2z_S_Py_S_vrr+WQY*I_ERI_Hx2y2z_S_Py_S_M1_vrr+oned2e*I_ERI_Hx2y2z_S_S_S_vrr-rhod2esq*I_ERI_Hx2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gxy2z_S_Py_S_M1_vrr;
      Double I_ERI_Hxy3z_S_D2y_S_vrr = QCY*I_ERI_Hxy3z_S_Py_S_vrr+WQY*I_ERI_Hxy3z_S_Py_S_M1_vrr+oned2e*I_ERI_Hxy3z_S_S_S_vrr-rhod2esq*I_ERI_Hxy3z_S_S_S_M1_vrr+oned2k*I_ERI_Gx3z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4z_S_D2y_S_vrr = QCY*I_ERI_Hx4z_S_Py_S_vrr+WQY*I_ERI_Hx4z_S_Py_S_M1_vrr+oned2e*I_ERI_Hx4z_S_S_S_vrr-rhod2esq*I_ERI_Hx4z_S_S_S_M1_vrr;
      Double I_ERI_H5y_S_D2y_S_vrr = QCY*I_ERI_H5y_S_Py_S_vrr+WQY*I_ERI_H5y_S_Py_S_M1_vrr+oned2e*I_ERI_H5y_S_S_S_vrr-rhod2esq*I_ERI_H5y_S_S_S_M1_vrr+5*oned2k*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_H4yz_S_D2y_S_vrr = QCY*I_ERI_H4yz_S_Py_S_vrr+WQY*I_ERI_H4yz_S_Py_S_M1_vrr+oned2e*I_ERI_H4yz_S_S_S_vrr-rhod2esq*I_ERI_H4yz_S_S_S_M1_vrr+4*oned2k*I_ERI_G3yz_S_Py_S_M1_vrr;
      Double I_ERI_H3y2z_S_D2y_S_vrr = QCY*I_ERI_H3y2z_S_Py_S_vrr+WQY*I_ERI_H3y2z_S_Py_S_M1_vrr+oned2e*I_ERI_H3y2z_S_S_S_vrr-rhod2esq*I_ERI_H3y2z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Py_S_M1_vrr;
      Double I_ERI_H2y3z_S_D2y_S_vrr = QCY*I_ERI_H2y3z_S_Py_S_vrr+WQY*I_ERI_H2y3z_S_Py_S_M1_vrr+oned2e*I_ERI_H2y3z_S_S_S_vrr-rhod2esq*I_ERI_H2y3z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gy3z_S_Py_S_M1_vrr;
      Double I_ERI_Hy4z_S_D2y_S_vrr = QCY*I_ERI_Hy4z_S_Py_S_vrr+WQY*I_ERI_Hy4z_S_Py_S_M1_vrr+oned2e*I_ERI_Hy4z_S_S_S_vrr-rhod2esq*I_ERI_Hy4z_S_S_S_M1_vrr+oned2k*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5z_S_D2y_S_vrr = QCY*I_ERI_H5z_S_Py_S_vrr+WQY*I_ERI_H5z_S_Py_S_M1_vrr+oned2e*I_ERI_H5z_S_S_S_vrr-rhod2esq*I_ERI_H5z_S_S_S_M1_vrr;
      Double I_ERI_H5x_S_Dyz_S_vrr = QCZ*I_ERI_H5x_S_Py_S_vrr+WQZ*I_ERI_H5x_S_Py_S_M1_vrr;
      Double I_ERI_H4xy_S_Dyz_S_vrr = QCZ*I_ERI_H4xy_S_Py_S_vrr+WQZ*I_ERI_H4xy_S_Py_S_M1_vrr;
      Double I_ERI_H4xz_S_Dyz_S_vrr = QCZ*I_ERI_H4xz_S_Py_S_vrr+WQZ*I_ERI_H4xz_S_Py_S_M1_vrr+oned2k*I_ERI_G4x_S_Py_S_M1_vrr;
      Double I_ERI_H3x2y_S_Dyz_S_vrr = QCZ*I_ERI_H3x2y_S_Py_S_vrr+WQZ*I_ERI_H3x2y_S_Py_S_M1_vrr;
      Double I_ERI_H3xyz_S_Dyz_S_vrr = QCZ*I_ERI_H3xyz_S_Py_S_vrr+WQZ*I_ERI_H3xyz_S_Py_S_M1_vrr+oned2k*I_ERI_G3xy_S_Py_S_M1_vrr;
      Double I_ERI_H3x2z_S_Dyz_S_vrr = QCZ*I_ERI_H3x2z_S_Py_S_vrr+WQZ*I_ERI_H3x2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_Py_S_M1_vrr;
      Double I_ERI_H2x3y_S_Dyz_S_vrr = QCZ*I_ERI_H2x3y_S_Py_S_vrr+WQZ*I_ERI_H2x3y_S_Py_S_M1_vrr;
      Double I_ERI_H2x2yz_S_Dyz_S_vrr = QCZ*I_ERI_H2x2yz_S_Py_S_vrr+WQZ*I_ERI_H2x2yz_S_Py_S_M1_vrr+oned2k*I_ERI_G2x2y_S_Py_S_M1_vrr;
      Double I_ERI_H2xy2z_S_Dyz_S_vrr = QCZ*I_ERI_H2xy2z_S_Py_S_vrr+WQZ*I_ERI_H2xy2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Py_S_M1_vrr;
      Double I_ERI_H2x3z_S_Dyz_S_vrr = QCZ*I_ERI_H2x3z_S_Py_S_vrr+WQZ*I_ERI_H2x3z_S_Py_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4y_S_Dyz_S_vrr = QCZ*I_ERI_Hx4y_S_Py_S_vrr+WQZ*I_ERI_Hx4y_S_Py_S_M1_vrr;
      Double I_ERI_Hx3yz_S_Dyz_S_vrr = QCZ*I_ERI_Hx3yz_S_Py_S_vrr+WQZ*I_ERI_Hx3yz_S_Py_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Py_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_Dyz_S_vrr = QCZ*I_ERI_Hx2y2z_S_Py_S_vrr+WQZ*I_ERI_Hx2y2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_Py_S_M1_vrr;
      Double I_ERI_Hxy3z_S_Dyz_S_vrr = QCZ*I_ERI_Hxy3z_S_Py_S_vrr+WQZ*I_ERI_Hxy3z_S_Py_S_M1_vrr+3*oned2k*I_ERI_Gxy2z_S_Py_S_M1_vrr;
      Double I_ERI_Hx4z_S_Dyz_S_vrr = QCZ*I_ERI_Hx4z_S_Py_S_vrr+WQZ*I_ERI_Hx4z_S_Py_S_M1_vrr+4*oned2k*I_ERI_Gx3z_S_Py_S_M1_vrr;
      Double I_ERI_H5y_S_Dyz_S_vrr = QCZ*I_ERI_H5y_S_Py_S_vrr+WQZ*I_ERI_H5y_S_Py_S_M1_vrr;
      Double I_ERI_H4yz_S_Dyz_S_vrr = QCZ*I_ERI_H4yz_S_Py_S_vrr+WQZ*I_ERI_H4yz_S_Py_S_M1_vrr+oned2k*I_ERI_G4y_S_Py_S_M1_vrr;
      Double I_ERI_H3y2z_S_Dyz_S_vrr = QCZ*I_ERI_H3y2z_S_Py_S_vrr+WQZ*I_ERI_H3y2z_S_Py_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_Py_S_M1_vrr;
      Double I_ERI_H2y3z_S_Dyz_S_vrr = QCZ*I_ERI_H2y3z_S_Py_S_vrr+WQZ*I_ERI_H2y3z_S_Py_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Py_S_M1_vrr;
      Double I_ERI_Hy4z_S_Dyz_S_vrr = QCZ*I_ERI_Hy4z_S_Py_S_vrr+WQZ*I_ERI_Hy4z_S_Py_S_M1_vrr+4*oned2k*I_ERI_Gy3z_S_Py_S_M1_vrr;
      Double I_ERI_H5z_S_Dyz_S_vrr = QCZ*I_ERI_H5z_S_Py_S_vrr+WQZ*I_ERI_H5z_S_Py_S_M1_vrr+5*oned2k*I_ERI_G4z_S_Py_S_M1_vrr;
      Double I_ERI_H5x_S_D2z_S_vrr = QCZ*I_ERI_H5x_S_Pz_S_vrr+WQZ*I_ERI_H5x_S_Pz_S_M1_vrr+oned2e*I_ERI_H5x_S_S_S_vrr-rhod2esq*I_ERI_H5x_S_S_S_M1_vrr;
      Double I_ERI_H4xy_S_D2z_S_vrr = QCZ*I_ERI_H4xy_S_Pz_S_vrr+WQZ*I_ERI_H4xy_S_Pz_S_M1_vrr+oned2e*I_ERI_H4xy_S_S_S_vrr-rhod2esq*I_ERI_H4xy_S_S_S_M1_vrr;
      Double I_ERI_H4xz_S_D2z_S_vrr = QCZ*I_ERI_H4xz_S_Pz_S_vrr+WQZ*I_ERI_H4xz_S_Pz_S_M1_vrr+oned2e*I_ERI_H4xz_S_S_S_vrr-rhod2esq*I_ERI_H4xz_S_S_S_M1_vrr+oned2k*I_ERI_G4x_S_Pz_S_M1_vrr;
      Double I_ERI_H3x2y_S_D2z_S_vrr = QCZ*I_ERI_H3x2y_S_Pz_S_vrr+WQZ*I_ERI_H3x2y_S_Pz_S_M1_vrr+oned2e*I_ERI_H3x2y_S_S_S_vrr-rhod2esq*I_ERI_H3x2y_S_S_S_M1_vrr;
      Double I_ERI_H3xyz_S_D2z_S_vrr = QCZ*I_ERI_H3xyz_S_Pz_S_vrr+WQZ*I_ERI_H3xyz_S_Pz_S_M1_vrr+oned2e*I_ERI_H3xyz_S_S_S_vrr-rhod2esq*I_ERI_H3xyz_S_S_S_M1_vrr+oned2k*I_ERI_G3xy_S_Pz_S_M1_vrr;
      Double I_ERI_H3x2z_S_D2z_S_vrr = QCZ*I_ERI_H3x2z_S_Pz_S_vrr+WQZ*I_ERI_H3x2z_S_Pz_S_M1_vrr+oned2e*I_ERI_H3x2z_S_S_S_vrr-rhod2esq*I_ERI_H3x2z_S_S_S_M1_vrr+2*oned2k*I_ERI_G3xz_S_Pz_S_M1_vrr;
      Double I_ERI_H2x3y_S_D2z_S_vrr = QCZ*I_ERI_H2x3y_S_Pz_S_vrr+WQZ*I_ERI_H2x3y_S_Pz_S_M1_vrr+oned2e*I_ERI_H2x3y_S_S_S_vrr-rhod2esq*I_ERI_H2x3y_S_S_S_M1_vrr;
      Double I_ERI_H2x2yz_S_D2z_S_vrr = QCZ*I_ERI_H2x2yz_S_Pz_S_vrr+WQZ*I_ERI_H2x2yz_S_Pz_S_M1_vrr+oned2e*I_ERI_H2x2yz_S_S_S_vrr-rhod2esq*I_ERI_H2x2yz_S_S_S_M1_vrr+oned2k*I_ERI_G2x2y_S_Pz_S_M1_vrr;
      Double I_ERI_H2xy2z_S_D2z_S_vrr = QCZ*I_ERI_H2xy2z_S_Pz_S_vrr+WQZ*I_ERI_H2xy2z_S_Pz_S_M1_vrr+oned2e*I_ERI_H2xy2z_S_S_S_vrr-rhod2esq*I_ERI_H2xy2z_S_S_S_M1_vrr+2*oned2k*I_ERI_G2xyz_S_Pz_S_M1_vrr;
      Double I_ERI_H2x3z_S_D2z_S_vrr = QCZ*I_ERI_H2x3z_S_Pz_S_vrr+WQZ*I_ERI_H2x3z_S_Pz_S_M1_vrr+oned2e*I_ERI_H2x3z_S_S_S_vrr-rhod2esq*I_ERI_H2x3z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2x2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4y_S_D2z_S_vrr = QCZ*I_ERI_Hx4y_S_Pz_S_vrr+WQZ*I_ERI_Hx4y_S_Pz_S_M1_vrr+oned2e*I_ERI_Hx4y_S_S_S_vrr-rhod2esq*I_ERI_Hx4y_S_S_S_M1_vrr;
      Double I_ERI_Hx3yz_S_D2z_S_vrr = QCZ*I_ERI_Hx3yz_S_Pz_S_vrr+WQZ*I_ERI_Hx3yz_S_Pz_S_M1_vrr+oned2e*I_ERI_Hx3yz_S_S_S_vrr-rhod2esq*I_ERI_Hx3yz_S_S_S_M1_vrr+oned2k*I_ERI_Gx3y_S_Pz_S_M1_vrr;
      Double I_ERI_Hx2y2z_S_D2z_S_vrr = QCZ*I_ERI_Hx2y2z_S_Pz_S_vrr+WQZ*I_ERI_Hx2y2z_S_Pz_S_M1_vrr+oned2e*I_ERI_Hx2y2z_S_S_S_vrr-rhod2esq*I_ERI_Hx2y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_Gx2yz_S_Pz_S_M1_vrr;
      Double I_ERI_Hxy3z_S_D2z_S_vrr = QCZ*I_ERI_Hxy3z_S_Pz_S_vrr+WQZ*I_ERI_Hxy3z_S_Pz_S_M1_vrr+oned2e*I_ERI_Hxy3z_S_S_S_vrr-rhod2esq*I_ERI_Hxy3z_S_S_S_M1_vrr+3*oned2k*I_ERI_Gxy2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hx4z_S_D2z_S_vrr = QCZ*I_ERI_Hx4z_S_Pz_S_vrr+WQZ*I_ERI_Hx4z_S_Pz_S_M1_vrr+oned2e*I_ERI_Hx4z_S_S_S_vrr-rhod2esq*I_ERI_Hx4z_S_S_S_M1_vrr+4*oned2k*I_ERI_Gx3z_S_Pz_S_M1_vrr;
      Double I_ERI_H5y_S_D2z_S_vrr = QCZ*I_ERI_H5y_S_Pz_S_vrr+WQZ*I_ERI_H5y_S_Pz_S_M1_vrr+oned2e*I_ERI_H5y_S_S_S_vrr-rhod2esq*I_ERI_H5y_S_S_S_M1_vrr;
      Double I_ERI_H4yz_S_D2z_S_vrr = QCZ*I_ERI_H4yz_S_Pz_S_vrr+WQZ*I_ERI_H4yz_S_Pz_S_M1_vrr+oned2e*I_ERI_H4yz_S_S_S_vrr-rhod2esq*I_ERI_H4yz_S_S_S_M1_vrr+oned2k*I_ERI_G4y_S_Pz_S_M1_vrr;
      Double I_ERI_H3y2z_S_D2z_S_vrr = QCZ*I_ERI_H3y2z_S_Pz_S_vrr+WQZ*I_ERI_H3y2z_S_Pz_S_M1_vrr+oned2e*I_ERI_H3y2z_S_S_S_vrr-rhod2esq*I_ERI_H3y2z_S_S_S_M1_vrr+2*oned2k*I_ERI_G3yz_S_Pz_S_M1_vrr;
      Double I_ERI_H2y3z_S_D2z_S_vrr = QCZ*I_ERI_H2y3z_S_Pz_S_vrr+WQZ*I_ERI_H2y3z_S_Pz_S_M1_vrr+oned2e*I_ERI_H2y3z_S_S_S_vrr-rhod2esq*I_ERI_H2y3z_S_S_S_M1_vrr+3*oned2k*I_ERI_G2y2z_S_Pz_S_M1_vrr;
      Double I_ERI_Hy4z_S_D2z_S_vrr = QCZ*I_ERI_Hy4z_S_Pz_S_vrr+WQZ*I_ERI_Hy4z_S_Pz_S_M1_vrr+oned2e*I_ERI_Hy4z_S_S_S_vrr-rhod2esq*I_ERI_Hy4z_S_S_S_M1_vrr+4*oned2k*I_ERI_Gy3z_S_Pz_S_M1_vrr;
      Double I_ERI_H5z_S_D2z_S_vrr = QCZ*I_ERI_H5z_S_Pz_S_vrr+WQZ*I_ERI_H5z_S_Pz_S_M1_vrr+oned2e*I_ERI_H5z_S_S_S_vrr-rhod2esq*I_ERI_H5z_S_S_S_M1_vrr+5*oned2k*I_ERI_G4z_S_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_H_S_D_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_H5x_S_D2x_S += I_ERI_H5x_S_D2x_S_vrr;
      I_ERI_H4xy_S_D2x_S += I_ERI_H4xy_S_D2x_S_vrr;
      I_ERI_H4xz_S_D2x_S += I_ERI_H4xz_S_D2x_S_vrr;
      I_ERI_H3x2y_S_D2x_S += I_ERI_H3x2y_S_D2x_S_vrr;
      I_ERI_H3xyz_S_D2x_S += I_ERI_H3xyz_S_D2x_S_vrr;
      I_ERI_H3x2z_S_D2x_S += I_ERI_H3x2z_S_D2x_S_vrr;
      I_ERI_H2x3y_S_D2x_S += I_ERI_H2x3y_S_D2x_S_vrr;
      I_ERI_H2x2yz_S_D2x_S += I_ERI_H2x2yz_S_D2x_S_vrr;
      I_ERI_H2xy2z_S_D2x_S += I_ERI_H2xy2z_S_D2x_S_vrr;
      I_ERI_H2x3z_S_D2x_S += I_ERI_H2x3z_S_D2x_S_vrr;
      I_ERI_Hx4y_S_D2x_S += I_ERI_Hx4y_S_D2x_S_vrr;
      I_ERI_Hx3yz_S_D2x_S += I_ERI_Hx3yz_S_D2x_S_vrr;
      I_ERI_Hx2y2z_S_D2x_S += I_ERI_Hx2y2z_S_D2x_S_vrr;
      I_ERI_Hxy3z_S_D2x_S += I_ERI_Hxy3z_S_D2x_S_vrr;
      I_ERI_Hx4z_S_D2x_S += I_ERI_Hx4z_S_D2x_S_vrr;
      I_ERI_H5y_S_D2x_S += I_ERI_H5y_S_D2x_S_vrr;
      I_ERI_H4yz_S_D2x_S += I_ERI_H4yz_S_D2x_S_vrr;
      I_ERI_H3y2z_S_D2x_S += I_ERI_H3y2z_S_D2x_S_vrr;
      I_ERI_H2y3z_S_D2x_S += I_ERI_H2y3z_S_D2x_S_vrr;
      I_ERI_Hy4z_S_D2x_S += I_ERI_Hy4z_S_D2x_S_vrr;
      I_ERI_H5z_S_D2x_S += I_ERI_H5z_S_D2x_S_vrr;
      I_ERI_H5x_S_Dxy_S += I_ERI_H5x_S_Dxy_S_vrr;
      I_ERI_H4xy_S_Dxy_S += I_ERI_H4xy_S_Dxy_S_vrr;
      I_ERI_H4xz_S_Dxy_S += I_ERI_H4xz_S_Dxy_S_vrr;
      I_ERI_H3x2y_S_Dxy_S += I_ERI_H3x2y_S_Dxy_S_vrr;
      I_ERI_H3xyz_S_Dxy_S += I_ERI_H3xyz_S_Dxy_S_vrr;
      I_ERI_H3x2z_S_Dxy_S += I_ERI_H3x2z_S_Dxy_S_vrr;
      I_ERI_H2x3y_S_Dxy_S += I_ERI_H2x3y_S_Dxy_S_vrr;
      I_ERI_H2x2yz_S_Dxy_S += I_ERI_H2x2yz_S_Dxy_S_vrr;
      I_ERI_H2xy2z_S_Dxy_S += I_ERI_H2xy2z_S_Dxy_S_vrr;
      I_ERI_H2x3z_S_Dxy_S += I_ERI_H2x3z_S_Dxy_S_vrr;
      I_ERI_Hx4y_S_Dxy_S += I_ERI_Hx4y_S_Dxy_S_vrr;
      I_ERI_Hx3yz_S_Dxy_S += I_ERI_Hx3yz_S_Dxy_S_vrr;
      I_ERI_Hx2y2z_S_Dxy_S += I_ERI_Hx2y2z_S_Dxy_S_vrr;
      I_ERI_Hxy3z_S_Dxy_S += I_ERI_Hxy3z_S_Dxy_S_vrr;
      I_ERI_Hx4z_S_Dxy_S += I_ERI_Hx4z_S_Dxy_S_vrr;
      I_ERI_H5y_S_Dxy_S += I_ERI_H5y_S_Dxy_S_vrr;
      I_ERI_H4yz_S_Dxy_S += I_ERI_H4yz_S_Dxy_S_vrr;
      I_ERI_H3y2z_S_Dxy_S += I_ERI_H3y2z_S_Dxy_S_vrr;
      I_ERI_H2y3z_S_Dxy_S += I_ERI_H2y3z_S_Dxy_S_vrr;
      I_ERI_Hy4z_S_Dxy_S += I_ERI_Hy4z_S_Dxy_S_vrr;
      I_ERI_H5z_S_Dxy_S += I_ERI_H5z_S_Dxy_S_vrr;
      I_ERI_H5x_S_Dxz_S += I_ERI_H5x_S_Dxz_S_vrr;
      I_ERI_H4xy_S_Dxz_S += I_ERI_H4xy_S_Dxz_S_vrr;
      I_ERI_H4xz_S_Dxz_S += I_ERI_H4xz_S_Dxz_S_vrr;
      I_ERI_H3x2y_S_Dxz_S += I_ERI_H3x2y_S_Dxz_S_vrr;
      I_ERI_H3xyz_S_Dxz_S += I_ERI_H3xyz_S_Dxz_S_vrr;
      I_ERI_H3x2z_S_Dxz_S += I_ERI_H3x2z_S_Dxz_S_vrr;
      I_ERI_H2x3y_S_Dxz_S += I_ERI_H2x3y_S_Dxz_S_vrr;
      I_ERI_H2x2yz_S_Dxz_S += I_ERI_H2x2yz_S_Dxz_S_vrr;
      I_ERI_H2xy2z_S_Dxz_S += I_ERI_H2xy2z_S_Dxz_S_vrr;
      I_ERI_H2x3z_S_Dxz_S += I_ERI_H2x3z_S_Dxz_S_vrr;
      I_ERI_Hx4y_S_Dxz_S += I_ERI_Hx4y_S_Dxz_S_vrr;
      I_ERI_Hx3yz_S_Dxz_S += I_ERI_Hx3yz_S_Dxz_S_vrr;
      I_ERI_Hx2y2z_S_Dxz_S += I_ERI_Hx2y2z_S_Dxz_S_vrr;
      I_ERI_Hxy3z_S_Dxz_S += I_ERI_Hxy3z_S_Dxz_S_vrr;
      I_ERI_Hx4z_S_Dxz_S += I_ERI_Hx4z_S_Dxz_S_vrr;
      I_ERI_H5y_S_Dxz_S += I_ERI_H5y_S_Dxz_S_vrr;
      I_ERI_H4yz_S_Dxz_S += I_ERI_H4yz_S_Dxz_S_vrr;
      I_ERI_H3y2z_S_Dxz_S += I_ERI_H3y2z_S_Dxz_S_vrr;
      I_ERI_H2y3z_S_Dxz_S += I_ERI_H2y3z_S_Dxz_S_vrr;
      I_ERI_Hy4z_S_Dxz_S += I_ERI_Hy4z_S_Dxz_S_vrr;
      I_ERI_H5z_S_Dxz_S += I_ERI_H5z_S_Dxz_S_vrr;
      I_ERI_H5x_S_D2y_S += I_ERI_H5x_S_D2y_S_vrr;
      I_ERI_H4xy_S_D2y_S += I_ERI_H4xy_S_D2y_S_vrr;
      I_ERI_H4xz_S_D2y_S += I_ERI_H4xz_S_D2y_S_vrr;
      I_ERI_H3x2y_S_D2y_S += I_ERI_H3x2y_S_D2y_S_vrr;
      I_ERI_H3xyz_S_D2y_S += I_ERI_H3xyz_S_D2y_S_vrr;
      I_ERI_H3x2z_S_D2y_S += I_ERI_H3x2z_S_D2y_S_vrr;
      I_ERI_H2x3y_S_D2y_S += I_ERI_H2x3y_S_D2y_S_vrr;
      I_ERI_H2x2yz_S_D2y_S += I_ERI_H2x2yz_S_D2y_S_vrr;
      I_ERI_H2xy2z_S_D2y_S += I_ERI_H2xy2z_S_D2y_S_vrr;
      I_ERI_H2x3z_S_D2y_S += I_ERI_H2x3z_S_D2y_S_vrr;
      I_ERI_Hx4y_S_D2y_S += I_ERI_Hx4y_S_D2y_S_vrr;
      I_ERI_Hx3yz_S_D2y_S += I_ERI_Hx3yz_S_D2y_S_vrr;
      I_ERI_Hx2y2z_S_D2y_S += I_ERI_Hx2y2z_S_D2y_S_vrr;
      I_ERI_Hxy3z_S_D2y_S += I_ERI_Hxy3z_S_D2y_S_vrr;
      I_ERI_Hx4z_S_D2y_S += I_ERI_Hx4z_S_D2y_S_vrr;
      I_ERI_H5y_S_D2y_S += I_ERI_H5y_S_D2y_S_vrr;
      I_ERI_H4yz_S_D2y_S += I_ERI_H4yz_S_D2y_S_vrr;
      I_ERI_H3y2z_S_D2y_S += I_ERI_H3y2z_S_D2y_S_vrr;
      I_ERI_H2y3z_S_D2y_S += I_ERI_H2y3z_S_D2y_S_vrr;
      I_ERI_Hy4z_S_D2y_S += I_ERI_Hy4z_S_D2y_S_vrr;
      I_ERI_H5z_S_D2y_S += I_ERI_H5z_S_D2y_S_vrr;
      I_ERI_H5x_S_Dyz_S += I_ERI_H5x_S_Dyz_S_vrr;
      I_ERI_H4xy_S_Dyz_S += I_ERI_H4xy_S_Dyz_S_vrr;
      I_ERI_H4xz_S_Dyz_S += I_ERI_H4xz_S_Dyz_S_vrr;
      I_ERI_H3x2y_S_Dyz_S += I_ERI_H3x2y_S_Dyz_S_vrr;
      I_ERI_H3xyz_S_Dyz_S += I_ERI_H3xyz_S_Dyz_S_vrr;
      I_ERI_H3x2z_S_Dyz_S += I_ERI_H3x2z_S_Dyz_S_vrr;
      I_ERI_H2x3y_S_Dyz_S += I_ERI_H2x3y_S_Dyz_S_vrr;
      I_ERI_H2x2yz_S_Dyz_S += I_ERI_H2x2yz_S_Dyz_S_vrr;
      I_ERI_H2xy2z_S_Dyz_S += I_ERI_H2xy2z_S_Dyz_S_vrr;
      I_ERI_H2x3z_S_Dyz_S += I_ERI_H2x3z_S_Dyz_S_vrr;
      I_ERI_Hx4y_S_Dyz_S += I_ERI_Hx4y_S_Dyz_S_vrr;
      I_ERI_Hx3yz_S_Dyz_S += I_ERI_Hx3yz_S_Dyz_S_vrr;
      I_ERI_Hx2y2z_S_Dyz_S += I_ERI_Hx2y2z_S_Dyz_S_vrr;
      I_ERI_Hxy3z_S_Dyz_S += I_ERI_Hxy3z_S_Dyz_S_vrr;
      I_ERI_Hx4z_S_Dyz_S += I_ERI_Hx4z_S_Dyz_S_vrr;
      I_ERI_H5y_S_Dyz_S += I_ERI_H5y_S_Dyz_S_vrr;
      I_ERI_H4yz_S_Dyz_S += I_ERI_H4yz_S_Dyz_S_vrr;
      I_ERI_H3y2z_S_Dyz_S += I_ERI_H3y2z_S_Dyz_S_vrr;
      I_ERI_H2y3z_S_Dyz_S += I_ERI_H2y3z_S_Dyz_S_vrr;
      I_ERI_Hy4z_S_Dyz_S += I_ERI_Hy4z_S_Dyz_S_vrr;
      I_ERI_H5z_S_Dyz_S += I_ERI_H5z_S_Dyz_S_vrr;
      I_ERI_H5x_S_D2z_S += I_ERI_H5x_S_D2z_S_vrr;
      I_ERI_H4xy_S_D2z_S += I_ERI_H4xy_S_D2z_S_vrr;
      I_ERI_H4xz_S_D2z_S += I_ERI_H4xz_S_D2z_S_vrr;
      I_ERI_H3x2y_S_D2z_S += I_ERI_H3x2y_S_D2z_S_vrr;
      I_ERI_H3xyz_S_D2z_S += I_ERI_H3xyz_S_D2z_S_vrr;
      I_ERI_H3x2z_S_D2z_S += I_ERI_H3x2z_S_D2z_S_vrr;
      I_ERI_H2x3y_S_D2z_S += I_ERI_H2x3y_S_D2z_S_vrr;
      I_ERI_H2x2yz_S_D2z_S += I_ERI_H2x2yz_S_D2z_S_vrr;
      I_ERI_H2xy2z_S_D2z_S += I_ERI_H2xy2z_S_D2z_S_vrr;
      I_ERI_H2x3z_S_D2z_S += I_ERI_H2x3z_S_D2z_S_vrr;
      I_ERI_Hx4y_S_D2z_S += I_ERI_Hx4y_S_D2z_S_vrr;
      I_ERI_Hx3yz_S_D2z_S += I_ERI_Hx3yz_S_D2z_S_vrr;
      I_ERI_Hx2y2z_S_D2z_S += I_ERI_Hx2y2z_S_D2z_S_vrr;
      I_ERI_Hxy3z_S_D2z_S += I_ERI_Hxy3z_S_D2z_S_vrr;
      I_ERI_Hx4z_S_D2z_S += I_ERI_Hx4z_S_D2z_S_vrr;
      I_ERI_H5y_S_D2z_S += I_ERI_H5y_S_D2z_S_vrr;
      I_ERI_H4yz_S_D2z_S += I_ERI_H4yz_S_D2z_S_vrr;
      I_ERI_H3y2z_S_D2z_S += I_ERI_H3y2z_S_D2z_S_vrr;
      I_ERI_H2y3z_S_D2z_S += I_ERI_H2y3z_S_D2z_S_vrr;
      I_ERI_Hy4z_S_D2z_S += I_ERI_Hy4z_S_D2z_S_vrr;
      I_ERI_H5z_S_D2z_S += I_ERI_H5z_S_D2z_S_vrr;

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
       * shell quartet name: SQ_ERI_G_S_D_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_G4x_S_D2x_S += I_ERI_G4x_S_D2x_S_vrr;
      I_ERI_G3xy_S_D2x_S += I_ERI_G3xy_S_D2x_S_vrr;
      I_ERI_G3xz_S_D2x_S += I_ERI_G3xz_S_D2x_S_vrr;
      I_ERI_G2x2y_S_D2x_S += I_ERI_G2x2y_S_D2x_S_vrr;
      I_ERI_G2xyz_S_D2x_S += I_ERI_G2xyz_S_D2x_S_vrr;
      I_ERI_G2x2z_S_D2x_S += I_ERI_G2x2z_S_D2x_S_vrr;
      I_ERI_Gx3y_S_D2x_S += I_ERI_Gx3y_S_D2x_S_vrr;
      I_ERI_Gx2yz_S_D2x_S += I_ERI_Gx2yz_S_D2x_S_vrr;
      I_ERI_Gxy2z_S_D2x_S += I_ERI_Gxy2z_S_D2x_S_vrr;
      I_ERI_Gx3z_S_D2x_S += I_ERI_Gx3z_S_D2x_S_vrr;
      I_ERI_G4y_S_D2x_S += I_ERI_G4y_S_D2x_S_vrr;
      I_ERI_G3yz_S_D2x_S += I_ERI_G3yz_S_D2x_S_vrr;
      I_ERI_G2y2z_S_D2x_S += I_ERI_G2y2z_S_D2x_S_vrr;
      I_ERI_Gy3z_S_D2x_S += I_ERI_Gy3z_S_D2x_S_vrr;
      I_ERI_G4z_S_D2x_S += I_ERI_G4z_S_D2x_S_vrr;
      I_ERI_G4x_S_Dxy_S += I_ERI_G4x_S_Dxy_S_vrr;
      I_ERI_G3xy_S_Dxy_S += I_ERI_G3xy_S_Dxy_S_vrr;
      I_ERI_G3xz_S_Dxy_S += I_ERI_G3xz_S_Dxy_S_vrr;
      I_ERI_G2x2y_S_Dxy_S += I_ERI_G2x2y_S_Dxy_S_vrr;
      I_ERI_G2xyz_S_Dxy_S += I_ERI_G2xyz_S_Dxy_S_vrr;
      I_ERI_G2x2z_S_Dxy_S += I_ERI_G2x2z_S_Dxy_S_vrr;
      I_ERI_Gx3y_S_Dxy_S += I_ERI_Gx3y_S_Dxy_S_vrr;
      I_ERI_Gx2yz_S_Dxy_S += I_ERI_Gx2yz_S_Dxy_S_vrr;
      I_ERI_Gxy2z_S_Dxy_S += I_ERI_Gxy2z_S_Dxy_S_vrr;
      I_ERI_Gx3z_S_Dxy_S += I_ERI_Gx3z_S_Dxy_S_vrr;
      I_ERI_G4y_S_Dxy_S += I_ERI_G4y_S_Dxy_S_vrr;
      I_ERI_G3yz_S_Dxy_S += I_ERI_G3yz_S_Dxy_S_vrr;
      I_ERI_G2y2z_S_Dxy_S += I_ERI_G2y2z_S_Dxy_S_vrr;
      I_ERI_Gy3z_S_Dxy_S += I_ERI_Gy3z_S_Dxy_S_vrr;
      I_ERI_G4z_S_Dxy_S += I_ERI_G4z_S_Dxy_S_vrr;
      I_ERI_G4x_S_Dxz_S += I_ERI_G4x_S_Dxz_S_vrr;
      I_ERI_G3xy_S_Dxz_S += I_ERI_G3xy_S_Dxz_S_vrr;
      I_ERI_G3xz_S_Dxz_S += I_ERI_G3xz_S_Dxz_S_vrr;
      I_ERI_G2x2y_S_Dxz_S += I_ERI_G2x2y_S_Dxz_S_vrr;
      I_ERI_G2xyz_S_Dxz_S += I_ERI_G2xyz_S_Dxz_S_vrr;
      I_ERI_G2x2z_S_Dxz_S += I_ERI_G2x2z_S_Dxz_S_vrr;
      I_ERI_Gx3y_S_Dxz_S += I_ERI_Gx3y_S_Dxz_S_vrr;
      I_ERI_Gx2yz_S_Dxz_S += I_ERI_Gx2yz_S_Dxz_S_vrr;
      I_ERI_Gxy2z_S_Dxz_S += I_ERI_Gxy2z_S_Dxz_S_vrr;
      I_ERI_Gx3z_S_Dxz_S += I_ERI_Gx3z_S_Dxz_S_vrr;
      I_ERI_G4y_S_Dxz_S += I_ERI_G4y_S_Dxz_S_vrr;
      I_ERI_G3yz_S_Dxz_S += I_ERI_G3yz_S_Dxz_S_vrr;
      I_ERI_G2y2z_S_Dxz_S += I_ERI_G2y2z_S_Dxz_S_vrr;
      I_ERI_Gy3z_S_Dxz_S += I_ERI_Gy3z_S_Dxz_S_vrr;
      I_ERI_G4z_S_Dxz_S += I_ERI_G4z_S_Dxz_S_vrr;
      I_ERI_G4x_S_D2y_S += I_ERI_G4x_S_D2y_S_vrr;
      I_ERI_G3xy_S_D2y_S += I_ERI_G3xy_S_D2y_S_vrr;
      I_ERI_G3xz_S_D2y_S += I_ERI_G3xz_S_D2y_S_vrr;
      I_ERI_G2x2y_S_D2y_S += I_ERI_G2x2y_S_D2y_S_vrr;
      I_ERI_G2xyz_S_D2y_S += I_ERI_G2xyz_S_D2y_S_vrr;
      I_ERI_G2x2z_S_D2y_S += I_ERI_G2x2z_S_D2y_S_vrr;
      I_ERI_Gx3y_S_D2y_S += I_ERI_Gx3y_S_D2y_S_vrr;
      I_ERI_Gx2yz_S_D2y_S += I_ERI_Gx2yz_S_D2y_S_vrr;
      I_ERI_Gxy2z_S_D2y_S += I_ERI_Gxy2z_S_D2y_S_vrr;
      I_ERI_Gx3z_S_D2y_S += I_ERI_Gx3z_S_D2y_S_vrr;
      I_ERI_G4y_S_D2y_S += I_ERI_G4y_S_D2y_S_vrr;
      I_ERI_G3yz_S_D2y_S += I_ERI_G3yz_S_D2y_S_vrr;
      I_ERI_G2y2z_S_D2y_S += I_ERI_G2y2z_S_D2y_S_vrr;
      I_ERI_Gy3z_S_D2y_S += I_ERI_Gy3z_S_D2y_S_vrr;
      I_ERI_G4z_S_D2y_S += I_ERI_G4z_S_D2y_S_vrr;
      I_ERI_G4x_S_Dyz_S += I_ERI_G4x_S_Dyz_S_vrr;
      I_ERI_G3xy_S_Dyz_S += I_ERI_G3xy_S_Dyz_S_vrr;
      I_ERI_G3xz_S_Dyz_S += I_ERI_G3xz_S_Dyz_S_vrr;
      I_ERI_G2x2y_S_Dyz_S += I_ERI_G2x2y_S_Dyz_S_vrr;
      I_ERI_G2xyz_S_Dyz_S += I_ERI_G2xyz_S_Dyz_S_vrr;
      I_ERI_G2x2z_S_Dyz_S += I_ERI_G2x2z_S_Dyz_S_vrr;
      I_ERI_Gx3y_S_Dyz_S += I_ERI_Gx3y_S_Dyz_S_vrr;
      I_ERI_Gx2yz_S_Dyz_S += I_ERI_Gx2yz_S_Dyz_S_vrr;
      I_ERI_Gxy2z_S_Dyz_S += I_ERI_Gxy2z_S_Dyz_S_vrr;
      I_ERI_Gx3z_S_Dyz_S += I_ERI_Gx3z_S_Dyz_S_vrr;
      I_ERI_G4y_S_Dyz_S += I_ERI_G4y_S_Dyz_S_vrr;
      I_ERI_G3yz_S_Dyz_S += I_ERI_G3yz_S_Dyz_S_vrr;
      I_ERI_G2y2z_S_Dyz_S += I_ERI_G2y2z_S_Dyz_S_vrr;
      I_ERI_Gy3z_S_Dyz_S += I_ERI_Gy3z_S_Dyz_S_vrr;
      I_ERI_G4z_S_Dyz_S += I_ERI_G4z_S_Dyz_S_vrr;
      I_ERI_G4x_S_D2z_S += I_ERI_G4x_S_D2z_S_vrr;
      I_ERI_G3xy_S_D2z_S += I_ERI_G3xy_S_D2z_S_vrr;
      I_ERI_G3xz_S_D2z_S += I_ERI_G3xz_S_D2z_S_vrr;
      I_ERI_G2x2y_S_D2z_S += I_ERI_G2x2y_S_D2z_S_vrr;
      I_ERI_G2xyz_S_D2z_S += I_ERI_G2xyz_S_D2z_S_vrr;
      I_ERI_G2x2z_S_D2z_S += I_ERI_G2x2z_S_D2z_S_vrr;
      I_ERI_Gx3y_S_D2z_S += I_ERI_Gx3y_S_D2z_S_vrr;
      I_ERI_Gx2yz_S_D2z_S += I_ERI_Gx2yz_S_D2z_S_vrr;
      I_ERI_Gxy2z_S_D2z_S += I_ERI_Gxy2z_S_D2z_S_vrr;
      I_ERI_Gx3z_S_D2z_S += I_ERI_Gx3z_S_D2z_S_vrr;
      I_ERI_G4y_S_D2z_S += I_ERI_G4y_S_D2z_S_vrr;
      I_ERI_G3yz_S_D2z_S += I_ERI_G3yz_S_D2z_S_vrr;
      I_ERI_G2y2z_S_D2z_S += I_ERI_G2y2z_S_D2z_S_vrr;
      I_ERI_Gy3z_S_D2z_S += I_ERI_Gy3z_S_D2z_S_vrr;
      I_ERI_G4z_S_D2z_S += I_ERI_G4z_S_D2z_S_vrr;

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

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_D_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_F3x_S_D2x_S += I_ERI_F3x_S_D2x_S_vrr;
      I_ERI_F2xy_S_D2x_S += I_ERI_F2xy_S_D2x_S_vrr;
      I_ERI_F2xz_S_D2x_S += I_ERI_F2xz_S_D2x_S_vrr;
      I_ERI_Fx2y_S_D2x_S += I_ERI_Fx2y_S_D2x_S_vrr;
      I_ERI_Fxyz_S_D2x_S += I_ERI_Fxyz_S_D2x_S_vrr;
      I_ERI_Fx2z_S_D2x_S += I_ERI_Fx2z_S_D2x_S_vrr;
      I_ERI_F3y_S_D2x_S += I_ERI_F3y_S_D2x_S_vrr;
      I_ERI_F2yz_S_D2x_S += I_ERI_F2yz_S_D2x_S_vrr;
      I_ERI_Fy2z_S_D2x_S += I_ERI_Fy2z_S_D2x_S_vrr;
      I_ERI_F3z_S_D2x_S += I_ERI_F3z_S_D2x_S_vrr;
      I_ERI_F3x_S_Dxy_S += I_ERI_F3x_S_Dxy_S_vrr;
      I_ERI_F2xy_S_Dxy_S += I_ERI_F2xy_S_Dxy_S_vrr;
      I_ERI_F2xz_S_Dxy_S += I_ERI_F2xz_S_Dxy_S_vrr;
      I_ERI_Fx2y_S_Dxy_S += I_ERI_Fx2y_S_Dxy_S_vrr;
      I_ERI_Fxyz_S_Dxy_S += I_ERI_Fxyz_S_Dxy_S_vrr;
      I_ERI_Fx2z_S_Dxy_S += I_ERI_Fx2z_S_Dxy_S_vrr;
      I_ERI_F3y_S_Dxy_S += I_ERI_F3y_S_Dxy_S_vrr;
      I_ERI_F2yz_S_Dxy_S += I_ERI_F2yz_S_Dxy_S_vrr;
      I_ERI_Fy2z_S_Dxy_S += I_ERI_Fy2z_S_Dxy_S_vrr;
      I_ERI_F3z_S_Dxy_S += I_ERI_F3z_S_Dxy_S_vrr;
      I_ERI_F3x_S_Dxz_S += I_ERI_F3x_S_Dxz_S_vrr;
      I_ERI_F2xy_S_Dxz_S += I_ERI_F2xy_S_Dxz_S_vrr;
      I_ERI_F2xz_S_Dxz_S += I_ERI_F2xz_S_Dxz_S_vrr;
      I_ERI_Fx2y_S_Dxz_S += I_ERI_Fx2y_S_Dxz_S_vrr;
      I_ERI_Fxyz_S_Dxz_S += I_ERI_Fxyz_S_Dxz_S_vrr;
      I_ERI_Fx2z_S_Dxz_S += I_ERI_Fx2z_S_Dxz_S_vrr;
      I_ERI_F3y_S_Dxz_S += I_ERI_F3y_S_Dxz_S_vrr;
      I_ERI_F2yz_S_Dxz_S += I_ERI_F2yz_S_Dxz_S_vrr;
      I_ERI_Fy2z_S_Dxz_S += I_ERI_Fy2z_S_Dxz_S_vrr;
      I_ERI_F3z_S_Dxz_S += I_ERI_F3z_S_Dxz_S_vrr;
      I_ERI_F3x_S_D2y_S += I_ERI_F3x_S_D2y_S_vrr;
      I_ERI_F2xy_S_D2y_S += I_ERI_F2xy_S_D2y_S_vrr;
      I_ERI_F2xz_S_D2y_S += I_ERI_F2xz_S_D2y_S_vrr;
      I_ERI_Fx2y_S_D2y_S += I_ERI_Fx2y_S_D2y_S_vrr;
      I_ERI_Fxyz_S_D2y_S += I_ERI_Fxyz_S_D2y_S_vrr;
      I_ERI_Fx2z_S_D2y_S += I_ERI_Fx2z_S_D2y_S_vrr;
      I_ERI_F3y_S_D2y_S += I_ERI_F3y_S_D2y_S_vrr;
      I_ERI_F2yz_S_D2y_S += I_ERI_F2yz_S_D2y_S_vrr;
      I_ERI_Fy2z_S_D2y_S += I_ERI_Fy2z_S_D2y_S_vrr;
      I_ERI_F3z_S_D2y_S += I_ERI_F3z_S_D2y_S_vrr;
      I_ERI_F3x_S_Dyz_S += I_ERI_F3x_S_Dyz_S_vrr;
      I_ERI_F2xy_S_Dyz_S += I_ERI_F2xy_S_Dyz_S_vrr;
      I_ERI_F2xz_S_Dyz_S += I_ERI_F2xz_S_Dyz_S_vrr;
      I_ERI_Fx2y_S_Dyz_S += I_ERI_Fx2y_S_Dyz_S_vrr;
      I_ERI_Fxyz_S_Dyz_S += I_ERI_Fxyz_S_Dyz_S_vrr;
      I_ERI_Fx2z_S_Dyz_S += I_ERI_Fx2z_S_Dyz_S_vrr;
      I_ERI_F3y_S_Dyz_S += I_ERI_F3y_S_Dyz_S_vrr;
      I_ERI_F2yz_S_Dyz_S += I_ERI_F2yz_S_Dyz_S_vrr;
      I_ERI_Fy2z_S_Dyz_S += I_ERI_Fy2z_S_Dyz_S_vrr;
      I_ERI_F3z_S_Dyz_S += I_ERI_F3z_S_Dyz_S_vrr;
      I_ERI_F3x_S_D2z_S += I_ERI_F3x_S_D2z_S_vrr;
      I_ERI_F2xy_S_D2z_S += I_ERI_F2xy_S_D2z_S_vrr;
      I_ERI_F2xz_S_D2z_S += I_ERI_F2xz_S_D2z_S_vrr;
      I_ERI_Fx2y_S_D2z_S += I_ERI_Fx2y_S_D2z_S_vrr;
      I_ERI_Fxyz_S_D2z_S += I_ERI_Fxyz_S_D2z_S_vrr;
      I_ERI_Fx2z_S_D2z_S += I_ERI_Fx2z_S_D2z_S_vrr;
      I_ERI_F3y_S_D2z_S += I_ERI_F3y_S_D2z_S_vrr;
      I_ERI_F2yz_S_D2z_S += I_ERI_F2yz_S_D2z_S_vrr;
      I_ERI_Fy2z_S_D2z_S += I_ERI_Fy2z_S_D2z_S_vrr;
      I_ERI_F3z_S_D2z_S += I_ERI_F3z_S_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_ERI_F_S_P_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_ERI_F3x_S_Px_S += I_ERI_F3x_S_Px_S_vrr;
      I_ERI_F2xy_S_Px_S += I_ERI_F2xy_S_Px_S_vrr;
      I_ERI_F2xz_S_Px_S += I_ERI_F2xz_S_Px_S_vrr;
      I_ERI_Fx2y_S_Px_S += I_ERI_Fx2y_S_Px_S_vrr;
      I_ERI_Fxyz_S_Px_S += I_ERI_Fxyz_S_Px_S_vrr;
      I_ERI_Fx2z_S_Px_S += I_ERI_Fx2z_S_Px_S_vrr;
      I_ERI_F3y_S_Px_S += I_ERI_F3y_S_Px_S_vrr;
      I_ERI_F2yz_S_Px_S += I_ERI_F2yz_S_Px_S_vrr;
      I_ERI_Fy2z_S_Px_S += I_ERI_Fy2z_S_Px_S_vrr;
      I_ERI_F3z_S_Px_S += I_ERI_F3z_S_Px_S_vrr;
      I_ERI_F3x_S_Py_S += I_ERI_F3x_S_Py_S_vrr;
      I_ERI_F2xy_S_Py_S += I_ERI_F2xy_S_Py_S_vrr;
      I_ERI_F2xz_S_Py_S += I_ERI_F2xz_S_Py_S_vrr;
      I_ERI_Fx2y_S_Py_S += I_ERI_Fx2y_S_Py_S_vrr;
      I_ERI_Fxyz_S_Py_S += I_ERI_Fxyz_S_Py_S_vrr;
      I_ERI_Fx2z_S_Py_S += I_ERI_Fx2z_S_Py_S_vrr;
      I_ERI_F3y_S_Py_S += I_ERI_F3y_S_Py_S_vrr;
      I_ERI_F2yz_S_Py_S += I_ERI_F2yz_S_Py_S_vrr;
      I_ERI_Fy2z_S_Py_S += I_ERI_Fy2z_S_Py_S_vrr;
      I_ERI_F3z_S_Py_S += I_ERI_F3z_S_Py_S_vrr;
      I_ERI_F3x_S_Pz_S += I_ERI_F3x_S_Pz_S_vrr;
      I_ERI_F2xy_S_Pz_S += I_ERI_F2xy_S_Pz_S_vrr;
      I_ERI_F2xz_S_Pz_S += I_ERI_F2xz_S_Pz_S_vrr;
      I_ERI_Fx2y_S_Pz_S += I_ERI_Fx2y_S_Pz_S_vrr;
      I_ERI_Fxyz_S_Pz_S += I_ERI_Fxyz_S_Pz_S_vrr;
      I_ERI_Fx2z_S_Pz_S += I_ERI_Fx2z_S_Pz_S_vrr;
      I_ERI_F3y_S_Pz_S += I_ERI_F3y_S_Pz_S_vrr;
      I_ERI_F2yz_S_Pz_S += I_ERI_F2yz_S_Pz_S_vrr;
      I_ERI_Fy2z_S_Pz_S += I_ERI_Fy2z_S_Pz_S_vrr;
      I_ERI_F3z_S_Pz_S += I_ERI_F3z_S_Pz_S_vrr;
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
   * shell quartet name: SQ_ERI_F_S_P_P
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_F_S_D_S
   * RHS shell quartet name: SQ_ERI_F_S_P_S
   ************************************************************/
  Double I_ERI_F3x_S_Px_Px = I_ERI_F3x_S_D2x_S+CDX*I_ERI_F3x_S_Px_S;
  Double I_ERI_F2xy_S_Px_Px = I_ERI_F2xy_S_D2x_S+CDX*I_ERI_F2xy_S_Px_S;
  Double I_ERI_F2xz_S_Px_Px = I_ERI_F2xz_S_D2x_S+CDX*I_ERI_F2xz_S_Px_S;
  Double I_ERI_Fx2y_S_Px_Px = I_ERI_Fx2y_S_D2x_S+CDX*I_ERI_Fx2y_S_Px_S;
  Double I_ERI_Fxyz_S_Px_Px = I_ERI_Fxyz_S_D2x_S+CDX*I_ERI_Fxyz_S_Px_S;
  Double I_ERI_Fx2z_S_Px_Px = I_ERI_Fx2z_S_D2x_S+CDX*I_ERI_Fx2z_S_Px_S;
  Double I_ERI_F3y_S_Px_Px = I_ERI_F3y_S_D2x_S+CDX*I_ERI_F3y_S_Px_S;
  Double I_ERI_F2yz_S_Px_Px = I_ERI_F2yz_S_D2x_S+CDX*I_ERI_F2yz_S_Px_S;
  Double I_ERI_Fy2z_S_Px_Px = I_ERI_Fy2z_S_D2x_S+CDX*I_ERI_Fy2z_S_Px_S;
  Double I_ERI_F3z_S_Px_Px = I_ERI_F3z_S_D2x_S+CDX*I_ERI_F3z_S_Px_S;
  Double I_ERI_F3x_S_Py_Px = I_ERI_F3x_S_Dxy_S+CDX*I_ERI_F3x_S_Py_S;
  Double I_ERI_F2xy_S_Py_Px = I_ERI_F2xy_S_Dxy_S+CDX*I_ERI_F2xy_S_Py_S;
  Double I_ERI_F2xz_S_Py_Px = I_ERI_F2xz_S_Dxy_S+CDX*I_ERI_F2xz_S_Py_S;
  Double I_ERI_Fx2y_S_Py_Px = I_ERI_Fx2y_S_Dxy_S+CDX*I_ERI_Fx2y_S_Py_S;
  Double I_ERI_Fxyz_S_Py_Px = I_ERI_Fxyz_S_Dxy_S+CDX*I_ERI_Fxyz_S_Py_S;
  Double I_ERI_Fx2z_S_Py_Px = I_ERI_Fx2z_S_Dxy_S+CDX*I_ERI_Fx2z_S_Py_S;
  Double I_ERI_F3y_S_Py_Px = I_ERI_F3y_S_Dxy_S+CDX*I_ERI_F3y_S_Py_S;
  Double I_ERI_F2yz_S_Py_Px = I_ERI_F2yz_S_Dxy_S+CDX*I_ERI_F2yz_S_Py_S;
  Double I_ERI_Fy2z_S_Py_Px = I_ERI_Fy2z_S_Dxy_S+CDX*I_ERI_Fy2z_S_Py_S;
  Double I_ERI_F3z_S_Py_Px = I_ERI_F3z_S_Dxy_S+CDX*I_ERI_F3z_S_Py_S;
  Double I_ERI_F3x_S_Pz_Px = I_ERI_F3x_S_Dxz_S+CDX*I_ERI_F3x_S_Pz_S;
  Double I_ERI_F2xy_S_Pz_Px = I_ERI_F2xy_S_Dxz_S+CDX*I_ERI_F2xy_S_Pz_S;
  Double I_ERI_F2xz_S_Pz_Px = I_ERI_F2xz_S_Dxz_S+CDX*I_ERI_F2xz_S_Pz_S;
  Double I_ERI_Fx2y_S_Pz_Px = I_ERI_Fx2y_S_Dxz_S+CDX*I_ERI_Fx2y_S_Pz_S;
  Double I_ERI_Fxyz_S_Pz_Px = I_ERI_Fxyz_S_Dxz_S+CDX*I_ERI_Fxyz_S_Pz_S;
  Double I_ERI_Fx2z_S_Pz_Px = I_ERI_Fx2z_S_Dxz_S+CDX*I_ERI_Fx2z_S_Pz_S;
  Double I_ERI_F3y_S_Pz_Px = I_ERI_F3y_S_Dxz_S+CDX*I_ERI_F3y_S_Pz_S;
  Double I_ERI_F2yz_S_Pz_Px = I_ERI_F2yz_S_Dxz_S+CDX*I_ERI_F2yz_S_Pz_S;
  Double I_ERI_Fy2z_S_Pz_Px = I_ERI_Fy2z_S_Dxz_S+CDX*I_ERI_Fy2z_S_Pz_S;
  Double I_ERI_F3z_S_Pz_Px = I_ERI_F3z_S_Dxz_S+CDX*I_ERI_F3z_S_Pz_S;
  Double I_ERI_F3x_S_Px_Py = I_ERI_F3x_S_Dxy_S+CDY*I_ERI_F3x_S_Px_S;
  Double I_ERI_F2xy_S_Px_Py = I_ERI_F2xy_S_Dxy_S+CDY*I_ERI_F2xy_S_Px_S;
  Double I_ERI_F2xz_S_Px_Py = I_ERI_F2xz_S_Dxy_S+CDY*I_ERI_F2xz_S_Px_S;
  Double I_ERI_Fx2y_S_Px_Py = I_ERI_Fx2y_S_Dxy_S+CDY*I_ERI_Fx2y_S_Px_S;
  Double I_ERI_Fxyz_S_Px_Py = I_ERI_Fxyz_S_Dxy_S+CDY*I_ERI_Fxyz_S_Px_S;
  Double I_ERI_Fx2z_S_Px_Py = I_ERI_Fx2z_S_Dxy_S+CDY*I_ERI_Fx2z_S_Px_S;
  Double I_ERI_F3y_S_Px_Py = I_ERI_F3y_S_Dxy_S+CDY*I_ERI_F3y_S_Px_S;
  Double I_ERI_F2yz_S_Px_Py = I_ERI_F2yz_S_Dxy_S+CDY*I_ERI_F2yz_S_Px_S;
  Double I_ERI_Fy2z_S_Px_Py = I_ERI_Fy2z_S_Dxy_S+CDY*I_ERI_Fy2z_S_Px_S;
  Double I_ERI_F3z_S_Px_Py = I_ERI_F3z_S_Dxy_S+CDY*I_ERI_F3z_S_Px_S;
  Double I_ERI_F3x_S_Py_Py = I_ERI_F3x_S_D2y_S+CDY*I_ERI_F3x_S_Py_S;
  Double I_ERI_F2xy_S_Py_Py = I_ERI_F2xy_S_D2y_S+CDY*I_ERI_F2xy_S_Py_S;
  Double I_ERI_F2xz_S_Py_Py = I_ERI_F2xz_S_D2y_S+CDY*I_ERI_F2xz_S_Py_S;
  Double I_ERI_Fx2y_S_Py_Py = I_ERI_Fx2y_S_D2y_S+CDY*I_ERI_Fx2y_S_Py_S;
  Double I_ERI_Fxyz_S_Py_Py = I_ERI_Fxyz_S_D2y_S+CDY*I_ERI_Fxyz_S_Py_S;
  Double I_ERI_Fx2z_S_Py_Py = I_ERI_Fx2z_S_D2y_S+CDY*I_ERI_Fx2z_S_Py_S;
  Double I_ERI_F3y_S_Py_Py = I_ERI_F3y_S_D2y_S+CDY*I_ERI_F3y_S_Py_S;
  Double I_ERI_F2yz_S_Py_Py = I_ERI_F2yz_S_D2y_S+CDY*I_ERI_F2yz_S_Py_S;
  Double I_ERI_Fy2z_S_Py_Py = I_ERI_Fy2z_S_D2y_S+CDY*I_ERI_Fy2z_S_Py_S;
  Double I_ERI_F3z_S_Py_Py = I_ERI_F3z_S_D2y_S+CDY*I_ERI_F3z_S_Py_S;
  Double I_ERI_F3x_S_Pz_Py = I_ERI_F3x_S_Dyz_S+CDY*I_ERI_F3x_S_Pz_S;
  Double I_ERI_F2xy_S_Pz_Py = I_ERI_F2xy_S_Dyz_S+CDY*I_ERI_F2xy_S_Pz_S;
  Double I_ERI_F2xz_S_Pz_Py = I_ERI_F2xz_S_Dyz_S+CDY*I_ERI_F2xz_S_Pz_S;
  Double I_ERI_Fx2y_S_Pz_Py = I_ERI_Fx2y_S_Dyz_S+CDY*I_ERI_Fx2y_S_Pz_S;
  Double I_ERI_Fxyz_S_Pz_Py = I_ERI_Fxyz_S_Dyz_S+CDY*I_ERI_Fxyz_S_Pz_S;
  Double I_ERI_Fx2z_S_Pz_Py = I_ERI_Fx2z_S_Dyz_S+CDY*I_ERI_Fx2z_S_Pz_S;
  Double I_ERI_F3y_S_Pz_Py = I_ERI_F3y_S_Dyz_S+CDY*I_ERI_F3y_S_Pz_S;
  Double I_ERI_F2yz_S_Pz_Py = I_ERI_F2yz_S_Dyz_S+CDY*I_ERI_F2yz_S_Pz_S;
  Double I_ERI_Fy2z_S_Pz_Py = I_ERI_Fy2z_S_Dyz_S+CDY*I_ERI_Fy2z_S_Pz_S;
  Double I_ERI_F3z_S_Pz_Py = I_ERI_F3z_S_Dyz_S+CDY*I_ERI_F3z_S_Pz_S;
  Double I_ERI_F3x_S_Px_Pz = I_ERI_F3x_S_Dxz_S+CDZ*I_ERI_F3x_S_Px_S;
  Double I_ERI_F2xy_S_Px_Pz = I_ERI_F2xy_S_Dxz_S+CDZ*I_ERI_F2xy_S_Px_S;
  Double I_ERI_F2xz_S_Px_Pz = I_ERI_F2xz_S_Dxz_S+CDZ*I_ERI_F2xz_S_Px_S;
  Double I_ERI_Fx2y_S_Px_Pz = I_ERI_Fx2y_S_Dxz_S+CDZ*I_ERI_Fx2y_S_Px_S;
  Double I_ERI_Fxyz_S_Px_Pz = I_ERI_Fxyz_S_Dxz_S+CDZ*I_ERI_Fxyz_S_Px_S;
  Double I_ERI_Fx2z_S_Px_Pz = I_ERI_Fx2z_S_Dxz_S+CDZ*I_ERI_Fx2z_S_Px_S;
  Double I_ERI_F3y_S_Px_Pz = I_ERI_F3y_S_Dxz_S+CDZ*I_ERI_F3y_S_Px_S;
  Double I_ERI_F2yz_S_Px_Pz = I_ERI_F2yz_S_Dxz_S+CDZ*I_ERI_F2yz_S_Px_S;
  Double I_ERI_Fy2z_S_Px_Pz = I_ERI_Fy2z_S_Dxz_S+CDZ*I_ERI_Fy2z_S_Px_S;
  Double I_ERI_F3z_S_Px_Pz = I_ERI_F3z_S_Dxz_S+CDZ*I_ERI_F3z_S_Px_S;
  Double I_ERI_F3x_S_Py_Pz = I_ERI_F3x_S_Dyz_S+CDZ*I_ERI_F3x_S_Py_S;
  Double I_ERI_F2xy_S_Py_Pz = I_ERI_F2xy_S_Dyz_S+CDZ*I_ERI_F2xy_S_Py_S;
  Double I_ERI_F2xz_S_Py_Pz = I_ERI_F2xz_S_Dyz_S+CDZ*I_ERI_F2xz_S_Py_S;
  Double I_ERI_Fx2y_S_Py_Pz = I_ERI_Fx2y_S_Dyz_S+CDZ*I_ERI_Fx2y_S_Py_S;
  Double I_ERI_Fxyz_S_Py_Pz = I_ERI_Fxyz_S_Dyz_S+CDZ*I_ERI_Fxyz_S_Py_S;
  Double I_ERI_Fx2z_S_Py_Pz = I_ERI_Fx2z_S_Dyz_S+CDZ*I_ERI_Fx2z_S_Py_S;
  Double I_ERI_F3y_S_Py_Pz = I_ERI_F3y_S_Dyz_S+CDZ*I_ERI_F3y_S_Py_S;
  Double I_ERI_F2yz_S_Py_Pz = I_ERI_F2yz_S_Dyz_S+CDZ*I_ERI_F2yz_S_Py_S;
  Double I_ERI_Fy2z_S_Py_Pz = I_ERI_Fy2z_S_Dyz_S+CDZ*I_ERI_Fy2z_S_Py_S;
  Double I_ERI_F3z_S_Py_Pz = I_ERI_F3z_S_Dyz_S+CDZ*I_ERI_F3z_S_Py_S;
  Double I_ERI_F3x_S_Pz_Pz = I_ERI_F3x_S_D2z_S+CDZ*I_ERI_F3x_S_Pz_S;
  Double I_ERI_F2xy_S_Pz_Pz = I_ERI_F2xy_S_D2z_S+CDZ*I_ERI_F2xy_S_Pz_S;
  Double I_ERI_F2xz_S_Pz_Pz = I_ERI_F2xz_S_D2z_S+CDZ*I_ERI_F2xz_S_Pz_S;
  Double I_ERI_Fx2y_S_Pz_Pz = I_ERI_Fx2y_S_D2z_S+CDZ*I_ERI_Fx2y_S_Pz_S;
  Double I_ERI_Fxyz_S_Pz_Pz = I_ERI_Fxyz_S_D2z_S+CDZ*I_ERI_Fxyz_S_Pz_S;
  Double I_ERI_Fx2z_S_Pz_Pz = I_ERI_Fx2z_S_D2z_S+CDZ*I_ERI_Fx2z_S_Pz_S;
  Double I_ERI_F3y_S_Pz_Pz = I_ERI_F3y_S_D2z_S+CDZ*I_ERI_F3y_S_Pz_S;
  Double I_ERI_F2yz_S_Pz_Pz = I_ERI_F2yz_S_D2z_S+CDZ*I_ERI_F2yz_S_Pz_S;
  Double I_ERI_Fy2z_S_Pz_Pz = I_ERI_Fy2z_S_D2z_S+CDZ*I_ERI_Fy2z_S_Pz_S;
  Double I_ERI_F3z_S_Pz_Pz = I_ERI_F3z_S_D2z_S+CDZ*I_ERI_F3z_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_G_S_P_P
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_D_S
   * RHS shell quartet name: SQ_ERI_G_S_P_S
   ************************************************************/
  Double I_ERI_G4x_S_Px_Px = I_ERI_G4x_S_D2x_S+CDX*I_ERI_G4x_S_Px_S;
  Double I_ERI_G3xy_S_Px_Px = I_ERI_G3xy_S_D2x_S+CDX*I_ERI_G3xy_S_Px_S;
  Double I_ERI_G3xz_S_Px_Px = I_ERI_G3xz_S_D2x_S+CDX*I_ERI_G3xz_S_Px_S;
  Double I_ERI_G2x2y_S_Px_Px = I_ERI_G2x2y_S_D2x_S+CDX*I_ERI_G2x2y_S_Px_S;
  Double I_ERI_G2xyz_S_Px_Px = I_ERI_G2xyz_S_D2x_S+CDX*I_ERI_G2xyz_S_Px_S;
  Double I_ERI_G2x2z_S_Px_Px = I_ERI_G2x2z_S_D2x_S+CDX*I_ERI_G2x2z_S_Px_S;
  Double I_ERI_Gx3y_S_Px_Px = I_ERI_Gx3y_S_D2x_S+CDX*I_ERI_Gx3y_S_Px_S;
  Double I_ERI_Gx2yz_S_Px_Px = I_ERI_Gx2yz_S_D2x_S+CDX*I_ERI_Gx2yz_S_Px_S;
  Double I_ERI_Gxy2z_S_Px_Px = I_ERI_Gxy2z_S_D2x_S+CDX*I_ERI_Gxy2z_S_Px_S;
  Double I_ERI_Gx3z_S_Px_Px = I_ERI_Gx3z_S_D2x_S+CDX*I_ERI_Gx3z_S_Px_S;
  Double I_ERI_G4y_S_Px_Px = I_ERI_G4y_S_D2x_S+CDX*I_ERI_G4y_S_Px_S;
  Double I_ERI_G3yz_S_Px_Px = I_ERI_G3yz_S_D2x_S+CDX*I_ERI_G3yz_S_Px_S;
  Double I_ERI_G2y2z_S_Px_Px = I_ERI_G2y2z_S_D2x_S+CDX*I_ERI_G2y2z_S_Px_S;
  Double I_ERI_Gy3z_S_Px_Px = I_ERI_Gy3z_S_D2x_S+CDX*I_ERI_Gy3z_S_Px_S;
  Double I_ERI_G4z_S_Px_Px = I_ERI_G4z_S_D2x_S+CDX*I_ERI_G4z_S_Px_S;
  Double I_ERI_G4x_S_Py_Px = I_ERI_G4x_S_Dxy_S+CDX*I_ERI_G4x_S_Py_S;
  Double I_ERI_G3xy_S_Py_Px = I_ERI_G3xy_S_Dxy_S+CDX*I_ERI_G3xy_S_Py_S;
  Double I_ERI_G3xz_S_Py_Px = I_ERI_G3xz_S_Dxy_S+CDX*I_ERI_G3xz_S_Py_S;
  Double I_ERI_G2x2y_S_Py_Px = I_ERI_G2x2y_S_Dxy_S+CDX*I_ERI_G2x2y_S_Py_S;
  Double I_ERI_G2xyz_S_Py_Px = I_ERI_G2xyz_S_Dxy_S+CDX*I_ERI_G2xyz_S_Py_S;
  Double I_ERI_G2x2z_S_Py_Px = I_ERI_G2x2z_S_Dxy_S+CDX*I_ERI_G2x2z_S_Py_S;
  Double I_ERI_Gx3y_S_Py_Px = I_ERI_Gx3y_S_Dxy_S+CDX*I_ERI_Gx3y_S_Py_S;
  Double I_ERI_Gx2yz_S_Py_Px = I_ERI_Gx2yz_S_Dxy_S+CDX*I_ERI_Gx2yz_S_Py_S;
  Double I_ERI_Gxy2z_S_Py_Px = I_ERI_Gxy2z_S_Dxy_S+CDX*I_ERI_Gxy2z_S_Py_S;
  Double I_ERI_Gx3z_S_Py_Px = I_ERI_Gx3z_S_Dxy_S+CDX*I_ERI_Gx3z_S_Py_S;
  Double I_ERI_G4y_S_Py_Px = I_ERI_G4y_S_Dxy_S+CDX*I_ERI_G4y_S_Py_S;
  Double I_ERI_G3yz_S_Py_Px = I_ERI_G3yz_S_Dxy_S+CDX*I_ERI_G3yz_S_Py_S;
  Double I_ERI_G2y2z_S_Py_Px = I_ERI_G2y2z_S_Dxy_S+CDX*I_ERI_G2y2z_S_Py_S;
  Double I_ERI_Gy3z_S_Py_Px = I_ERI_Gy3z_S_Dxy_S+CDX*I_ERI_Gy3z_S_Py_S;
  Double I_ERI_G4z_S_Py_Px = I_ERI_G4z_S_Dxy_S+CDX*I_ERI_G4z_S_Py_S;
  Double I_ERI_G4x_S_Pz_Px = I_ERI_G4x_S_Dxz_S+CDX*I_ERI_G4x_S_Pz_S;
  Double I_ERI_G3xy_S_Pz_Px = I_ERI_G3xy_S_Dxz_S+CDX*I_ERI_G3xy_S_Pz_S;
  Double I_ERI_G3xz_S_Pz_Px = I_ERI_G3xz_S_Dxz_S+CDX*I_ERI_G3xz_S_Pz_S;
  Double I_ERI_G2x2y_S_Pz_Px = I_ERI_G2x2y_S_Dxz_S+CDX*I_ERI_G2x2y_S_Pz_S;
  Double I_ERI_G2xyz_S_Pz_Px = I_ERI_G2xyz_S_Dxz_S+CDX*I_ERI_G2xyz_S_Pz_S;
  Double I_ERI_G2x2z_S_Pz_Px = I_ERI_G2x2z_S_Dxz_S+CDX*I_ERI_G2x2z_S_Pz_S;
  Double I_ERI_Gx3y_S_Pz_Px = I_ERI_Gx3y_S_Dxz_S+CDX*I_ERI_Gx3y_S_Pz_S;
  Double I_ERI_Gx2yz_S_Pz_Px = I_ERI_Gx2yz_S_Dxz_S+CDX*I_ERI_Gx2yz_S_Pz_S;
  Double I_ERI_Gxy2z_S_Pz_Px = I_ERI_Gxy2z_S_Dxz_S+CDX*I_ERI_Gxy2z_S_Pz_S;
  Double I_ERI_Gx3z_S_Pz_Px = I_ERI_Gx3z_S_Dxz_S+CDX*I_ERI_Gx3z_S_Pz_S;
  Double I_ERI_G4y_S_Pz_Px = I_ERI_G4y_S_Dxz_S+CDX*I_ERI_G4y_S_Pz_S;
  Double I_ERI_G3yz_S_Pz_Px = I_ERI_G3yz_S_Dxz_S+CDX*I_ERI_G3yz_S_Pz_S;
  Double I_ERI_G2y2z_S_Pz_Px = I_ERI_G2y2z_S_Dxz_S+CDX*I_ERI_G2y2z_S_Pz_S;
  Double I_ERI_Gy3z_S_Pz_Px = I_ERI_Gy3z_S_Dxz_S+CDX*I_ERI_Gy3z_S_Pz_S;
  Double I_ERI_G4z_S_Pz_Px = I_ERI_G4z_S_Dxz_S+CDX*I_ERI_G4z_S_Pz_S;
  Double I_ERI_G4x_S_Px_Py = I_ERI_G4x_S_Dxy_S+CDY*I_ERI_G4x_S_Px_S;
  Double I_ERI_G3xy_S_Px_Py = I_ERI_G3xy_S_Dxy_S+CDY*I_ERI_G3xy_S_Px_S;
  Double I_ERI_G3xz_S_Px_Py = I_ERI_G3xz_S_Dxy_S+CDY*I_ERI_G3xz_S_Px_S;
  Double I_ERI_G2x2y_S_Px_Py = I_ERI_G2x2y_S_Dxy_S+CDY*I_ERI_G2x2y_S_Px_S;
  Double I_ERI_G2xyz_S_Px_Py = I_ERI_G2xyz_S_Dxy_S+CDY*I_ERI_G2xyz_S_Px_S;
  Double I_ERI_G2x2z_S_Px_Py = I_ERI_G2x2z_S_Dxy_S+CDY*I_ERI_G2x2z_S_Px_S;
  Double I_ERI_Gx3y_S_Px_Py = I_ERI_Gx3y_S_Dxy_S+CDY*I_ERI_Gx3y_S_Px_S;
  Double I_ERI_Gx2yz_S_Px_Py = I_ERI_Gx2yz_S_Dxy_S+CDY*I_ERI_Gx2yz_S_Px_S;
  Double I_ERI_Gxy2z_S_Px_Py = I_ERI_Gxy2z_S_Dxy_S+CDY*I_ERI_Gxy2z_S_Px_S;
  Double I_ERI_Gx3z_S_Px_Py = I_ERI_Gx3z_S_Dxy_S+CDY*I_ERI_Gx3z_S_Px_S;
  Double I_ERI_G4y_S_Px_Py = I_ERI_G4y_S_Dxy_S+CDY*I_ERI_G4y_S_Px_S;
  Double I_ERI_G3yz_S_Px_Py = I_ERI_G3yz_S_Dxy_S+CDY*I_ERI_G3yz_S_Px_S;
  Double I_ERI_G2y2z_S_Px_Py = I_ERI_G2y2z_S_Dxy_S+CDY*I_ERI_G2y2z_S_Px_S;
  Double I_ERI_Gy3z_S_Px_Py = I_ERI_Gy3z_S_Dxy_S+CDY*I_ERI_Gy3z_S_Px_S;
  Double I_ERI_G4z_S_Px_Py = I_ERI_G4z_S_Dxy_S+CDY*I_ERI_G4z_S_Px_S;
  Double I_ERI_G4x_S_Py_Py = I_ERI_G4x_S_D2y_S+CDY*I_ERI_G4x_S_Py_S;
  Double I_ERI_G3xy_S_Py_Py = I_ERI_G3xy_S_D2y_S+CDY*I_ERI_G3xy_S_Py_S;
  Double I_ERI_G3xz_S_Py_Py = I_ERI_G3xz_S_D2y_S+CDY*I_ERI_G3xz_S_Py_S;
  Double I_ERI_G2x2y_S_Py_Py = I_ERI_G2x2y_S_D2y_S+CDY*I_ERI_G2x2y_S_Py_S;
  Double I_ERI_G2xyz_S_Py_Py = I_ERI_G2xyz_S_D2y_S+CDY*I_ERI_G2xyz_S_Py_S;
  Double I_ERI_G2x2z_S_Py_Py = I_ERI_G2x2z_S_D2y_S+CDY*I_ERI_G2x2z_S_Py_S;
  Double I_ERI_Gx3y_S_Py_Py = I_ERI_Gx3y_S_D2y_S+CDY*I_ERI_Gx3y_S_Py_S;
  Double I_ERI_Gx2yz_S_Py_Py = I_ERI_Gx2yz_S_D2y_S+CDY*I_ERI_Gx2yz_S_Py_S;
  Double I_ERI_Gxy2z_S_Py_Py = I_ERI_Gxy2z_S_D2y_S+CDY*I_ERI_Gxy2z_S_Py_S;
  Double I_ERI_Gx3z_S_Py_Py = I_ERI_Gx3z_S_D2y_S+CDY*I_ERI_Gx3z_S_Py_S;
  Double I_ERI_G4y_S_Py_Py = I_ERI_G4y_S_D2y_S+CDY*I_ERI_G4y_S_Py_S;
  Double I_ERI_G3yz_S_Py_Py = I_ERI_G3yz_S_D2y_S+CDY*I_ERI_G3yz_S_Py_S;
  Double I_ERI_G2y2z_S_Py_Py = I_ERI_G2y2z_S_D2y_S+CDY*I_ERI_G2y2z_S_Py_S;
  Double I_ERI_Gy3z_S_Py_Py = I_ERI_Gy3z_S_D2y_S+CDY*I_ERI_Gy3z_S_Py_S;
  Double I_ERI_G4z_S_Py_Py = I_ERI_G4z_S_D2y_S+CDY*I_ERI_G4z_S_Py_S;
  Double I_ERI_G4x_S_Pz_Py = I_ERI_G4x_S_Dyz_S+CDY*I_ERI_G4x_S_Pz_S;
  Double I_ERI_G3xy_S_Pz_Py = I_ERI_G3xy_S_Dyz_S+CDY*I_ERI_G3xy_S_Pz_S;
  Double I_ERI_G3xz_S_Pz_Py = I_ERI_G3xz_S_Dyz_S+CDY*I_ERI_G3xz_S_Pz_S;
  Double I_ERI_G2x2y_S_Pz_Py = I_ERI_G2x2y_S_Dyz_S+CDY*I_ERI_G2x2y_S_Pz_S;
  Double I_ERI_G2xyz_S_Pz_Py = I_ERI_G2xyz_S_Dyz_S+CDY*I_ERI_G2xyz_S_Pz_S;
  Double I_ERI_G2x2z_S_Pz_Py = I_ERI_G2x2z_S_Dyz_S+CDY*I_ERI_G2x2z_S_Pz_S;
  Double I_ERI_Gx3y_S_Pz_Py = I_ERI_Gx3y_S_Dyz_S+CDY*I_ERI_Gx3y_S_Pz_S;
  Double I_ERI_Gx2yz_S_Pz_Py = I_ERI_Gx2yz_S_Dyz_S+CDY*I_ERI_Gx2yz_S_Pz_S;
  Double I_ERI_Gxy2z_S_Pz_Py = I_ERI_Gxy2z_S_Dyz_S+CDY*I_ERI_Gxy2z_S_Pz_S;
  Double I_ERI_Gx3z_S_Pz_Py = I_ERI_Gx3z_S_Dyz_S+CDY*I_ERI_Gx3z_S_Pz_S;
  Double I_ERI_G4y_S_Pz_Py = I_ERI_G4y_S_Dyz_S+CDY*I_ERI_G4y_S_Pz_S;
  Double I_ERI_G3yz_S_Pz_Py = I_ERI_G3yz_S_Dyz_S+CDY*I_ERI_G3yz_S_Pz_S;
  Double I_ERI_G2y2z_S_Pz_Py = I_ERI_G2y2z_S_Dyz_S+CDY*I_ERI_G2y2z_S_Pz_S;
  Double I_ERI_Gy3z_S_Pz_Py = I_ERI_Gy3z_S_Dyz_S+CDY*I_ERI_Gy3z_S_Pz_S;
  Double I_ERI_G4z_S_Pz_Py = I_ERI_G4z_S_Dyz_S+CDY*I_ERI_G4z_S_Pz_S;
  Double I_ERI_G4x_S_Px_Pz = I_ERI_G4x_S_Dxz_S+CDZ*I_ERI_G4x_S_Px_S;
  Double I_ERI_G3xy_S_Px_Pz = I_ERI_G3xy_S_Dxz_S+CDZ*I_ERI_G3xy_S_Px_S;
  Double I_ERI_G3xz_S_Px_Pz = I_ERI_G3xz_S_Dxz_S+CDZ*I_ERI_G3xz_S_Px_S;
  Double I_ERI_G2x2y_S_Px_Pz = I_ERI_G2x2y_S_Dxz_S+CDZ*I_ERI_G2x2y_S_Px_S;
  Double I_ERI_G2xyz_S_Px_Pz = I_ERI_G2xyz_S_Dxz_S+CDZ*I_ERI_G2xyz_S_Px_S;
  Double I_ERI_G2x2z_S_Px_Pz = I_ERI_G2x2z_S_Dxz_S+CDZ*I_ERI_G2x2z_S_Px_S;
  Double I_ERI_Gx3y_S_Px_Pz = I_ERI_Gx3y_S_Dxz_S+CDZ*I_ERI_Gx3y_S_Px_S;
  Double I_ERI_Gx2yz_S_Px_Pz = I_ERI_Gx2yz_S_Dxz_S+CDZ*I_ERI_Gx2yz_S_Px_S;
  Double I_ERI_Gxy2z_S_Px_Pz = I_ERI_Gxy2z_S_Dxz_S+CDZ*I_ERI_Gxy2z_S_Px_S;
  Double I_ERI_Gx3z_S_Px_Pz = I_ERI_Gx3z_S_Dxz_S+CDZ*I_ERI_Gx3z_S_Px_S;
  Double I_ERI_G4y_S_Px_Pz = I_ERI_G4y_S_Dxz_S+CDZ*I_ERI_G4y_S_Px_S;
  Double I_ERI_G3yz_S_Px_Pz = I_ERI_G3yz_S_Dxz_S+CDZ*I_ERI_G3yz_S_Px_S;
  Double I_ERI_G2y2z_S_Px_Pz = I_ERI_G2y2z_S_Dxz_S+CDZ*I_ERI_G2y2z_S_Px_S;
  Double I_ERI_Gy3z_S_Px_Pz = I_ERI_Gy3z_S_Dxz_S+CDZ*I_ERI_Gy3z_S_Px_S;
  Double I_ERI_G4z_S_Px_Pz = I_ERI_G4z_S_Dxz_S+CDZ*I_ERI_G4z_S_Px_S;
  Double I_ERI_G4x_S_Py_Pz = I_ERI_G4x_S_Dyz_S+CDZ*I_ERI_G4x_S_Py_S;
  Double I_ERI_G3xy_S_Py_Pz = I_ERI_G3xy_S_Dyz_S+CDZ*I_ERI_G3xy_S_Py_S;
  Double I_ERI_G3xz_S_Py_Pz = I_ERI_G3xz_S_Dyz_S+CDZ*I_ERI_G3xz_S_Py_S;
  Double I_ERI_G2x2y_S_Py_Pz = I_ERI_G2x2y_S_Dyz_S+CDZ*I_ERI_G2x2y_S_Py_S;
  Double I_ERI_G2xyz_S_Py_Pz = I_ERI_G2xyz_S_Dyz_S+CDZ*I_ERI_G2xyz_S_Py_S;
  Double I_ERI_G2x2z_S_Py_Pz = I_ERI_G2x2z_S_Dyz_S+CDZ*I_ERI_G2x2z_S_Py_S;
  Double I_ERI_Gx3y_S_Py_Pz = I_ERI_Gx3y_S_Dyz_S+CDZ*I_ERI_Gx3y_S_Py_S;
  Double I_ERI_Gx2yz_S_Py_Pz = I_ERI_Gx2yz_S_Dyz_S+CDZ*I_ERI_Gx2yz_S_Py_S;
  Double I_ERI_Gxy2z_S_Py_Pz = I_ERI_Gxy2z_S_Dyz_S+CDZ*I_ERI_Gxy2z_S_Py_S;
  Double I_ERI_Gx3z_S_Py_Pz = I_ERI_Gx3z_S_Dyz_S+CDZ*I_ERI_Gx3z_S_Py_S;
  Double I_ERI_G4y_S_Py_Pz = I_ERI_G4y_S_Dyz_S+CDZ*I_ERI_G4y_S_Py_S;
  Double I_ERI_G3yz_S_Py_Pz = I_ERI_G3yz_S_Dyz_S+CDZ*I_ERI_G3yz_S_Py_S;
  Double I_ERI_G2y2z_S_Py_Pz = I_ERI_G2y2z_S_Dyz_S+CDZ*I_ERI_G2y2z_S_Py_S;
  Double I_ERI_Gy3z_S_Py_Pz = I_ERI_Gy3z_S_Dyz_S+CDZ*I_ERI_Gy3z_S_Py_S;
  Double I_ERI_G4z_S_Py_Pz = I_ERI_G4z_S_Dyz_S+CDZ*I_ERI_G4z_S_Py_S;
  Double I_ERI_G4x_S_Pz_Pz = I_ERI_G4x_S_D2z_S+CDZ*I_ERI_G4x_S_Pz_S;
  Double I_ERI_G3xy_S_Pz_Pz = I_ERI_G3xy_S_D2z_S+CDZ*I_ERI_G3xy_S_Pz_S;
  Double I_ERI_G3xz_S_Pz_Pz = I_ERI_G3xz_S_D2z_S+CDZ*I_ERI_G3xz_S_Pz_S;
  Double I_ERI_G2x2y_S_Pz_Pz = I_ERI_G2x2y_S_D2z_S+CDZ*I_ERI_G2x2y_S_Pz_S;
  Double I_ERI_G2xyz_S_Pz_Pz = I_ERI_G2xyz_S_D2z_S+CDZ*I_ERI_G2xyz_S_Pz_S;
  Double I_ERI_G2x2z_S_Pz_Pz = I_ERI_G2x2z_S_D2z_S+CDZ*I_ERI_G2x2z_S_Pz_S;
  Double I_ERI_Gx3y_S_Pz_Pz = I_ERI_Gx3y_S_D2z_S+CDZ*I_ERI_Gx3y_S_Pz_S;
  Double I_ERI_Gx2yz_S_Pz_Pz = I_ERI_Gx2yz_S_D2z_S+CDZ*I_ERI_Gx2yz_S_Pz_S;
  Double I_ERI_Gxy2z_S_Pz_Pz = I_ERI_Gxy2z_S_D2z_S+CDZ*I_ERI_Gxy2z_S_Pz_S;
  Double I_ERI_Gx3z_S_Pz_Pz = I_ERI_Gx3z_S_D2z_S+CDZ*I_ERI_Gx3z_S_Pz_S;
  Double I_ERI_G4y_S_Pz_Pz = I_ERI_G4y_S_D2z_S+CDZ*I_ERI_G4y_S_Pz_S;
  Double I_ERI_G3yz_S_Pz_Pz = I_ERI_G3yz_S_D2z_S+CDZ*I_ERI_G3yz_S_Pz_S;
  Double I_ERI_G2y2z_S_Pz_Pz = I_ERI_G2y2z_S_D2z_S+CDZ*I_ERI_G2y2z_S_Pz_S;
  Double I_ERI_Gy3z_S_Pz_Pz = I_ERI_Gy3z_S_D2z_S+CDZ*I_ERI_Gy3z_S_Pz_S;
  Double I_ERI_G4z_S_Pz_Pz = I_ERI_G4z_S_D2z_S+CDZ*I_ERI_G4z_S_Pz_S;

  /************************************************************
   * shell quartet name: SQ_ERI_H_S_P_P
   * expanding position: KET2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_D_S
   * RHS shell quartet name: SQ_ERI_H_S_P_S
   ************************************************************/
  Double I_ERI_H5x_S_Px_Px = I_ERI_H5x_S_D2x_S+CDX*I_ERI_H5x_S_Px_S;
  Double I_ERI_H4xy_S_Px_Px = I_ERI_H4xy_S_D2x_S+CDX*I_ERI_H4xy_S_Px_S;
  Double I_ERI_H4xz_S_Px_Px = I_ERI_H4xz_S_D2x_S+CDX*I_ERI_H4xz_S_Px_S;
  Double I_ERI_H3x2y_S_Px_Px = I_ERI_H3x2y_S_D2x_S+CDX*I_ERI_H3x2y_S_Px_S;
  Double I_ERI_H3xyz_S_Px_Px = I_ERI_H3xyz_S_D2x_S+CDX*I_ERI_H3xyz_S_Px_S;
  Double I_ERI_H3x2z_S_Px_Px = I_ERI_H3x2z_S_D2x_S+CDX*I_ERI_H3x2z_S_Px_S;
  Double I_ERI_H2x3y_S_Px_Px = I_ERI_H2x3y_S_D2x_S+CDX*I_ERI_H2x3y_S_Px_S;
  Double I_ERI_H2x2yz_S_Px_Px = I_ERI_H2x2yz_S_D2x_S+CDX*I_ERI_H2x2yz_S_Px_S;
  Double I_ERI_H2xy2z_S_Px_Px = I_ERI_H2xy2z_S_D2x_S+CDX*I_ERI_H2xy2z_S_Px_S;
  Double I_ERI_H2x3z_S_Px_Px = I_ERI_H2x3z_S_D2x_S+CDX*I_ERI_H2x3z_S_Px_S;
  Double I_ERI_Hx4y_S_Px_Px = I_ERI_Hx4y_S_D2x_S+CDX*I_ERI_Hx4y_S_Px_S;
  Double I_ERI_Hx3yz_S_Px_Px = I_ERI_Hx3yz_S_D2x_S+CDX*I_ERI_Hx3yz_S_Px_S;
  Double I_ERI_Hx2y2z_S_Px_Px = I_ERI_Hx2y2z_S_D2x_S+CDX*I_ERI_Hx2y2z_S_Px_S;
  Double I_ERI_Hxy3z_S_Px_Px = I_ERI_Hxy3z_S_D2x_S+CDX*I_ERI_Hxy3z_S_Px_S;
  Double I_ERI_Hx4z_S_Px_Px = I_ERI_Hx4z_S_D2x_S+CDX*I_ERI_Hx4z_S_Px_S;
  Double I_ERI_H5y_S_Px_Px = I_ERI_H5y_S_D2x_S+CDX*I_ERI_H5y_S_Px_S;
  Double I_ERI_H4yz_S_Px_Px = I_ERI_H4yz_S_D2x_S+CDX*I_ERI_H4yz_S_Px_S;
  Double I_ERI_H3y2z_S_Px_Px = I_ERI_H3y2z_S_D2x_S+CDX*I_ERI_H3y2z_S_Px_S;
  Double I_ERI_H2y3z_S_Px_Px = I_ERI_H2y3z_S_D2x_S+CDX*I_ERI_H2y3z_S_Px_S;
  Double I_ERI_Hy4z_S_Px_Px = I_ERI_Hy4z_S_D2x_S+CDX*I_ERI_Hy4z_S_Px_S;
  Double I_ERI_H5z_S_Px_Px = I_ERI_H5z_S_D2x_S+CDX*I_ERI_H5z_S_Px_S;
  Double I_ERI_H5x_S_Py_Px = I_ERI_H5x_S_Dxy_S+CDX*I_ERI_H5x_S_Py_S;
  Double I_ERI_H4xy_S_Py_Px = I_ERI_H4xy_S_Dxy_S+CDX*I_ERI_H4xy_S_Py_S;
  Double I_ERI_H4xz_S_Py_Px = I_ERI_H4xz_S_Dxy_S+CDX*I_ERI_H4xz_S_Py_S;
  Double I_ERI_H3x2y_S_Py_Px = I_ERI_H3x2y_S_Dxy_S+CDX*I_ERI_H3x2y_S_Py_S;
  Double I_ERI_H3xyz_S_Py_Px = I_ERI_H3xyz_S_Dxy_S+CDX*I_ERI_H3xyz_S_Py_S;
  Double I_ERI_H3x2z_S_Py_Px = I_ERI_H3x2z_S_Dxy_S+CDX*I_ERI_H3x2z_S_Py_S;
  Double I_ERI_H2x3y_S_Py_Px = I_ERI_H2x3y_S_Dxy_S+CDX*I_ERI_H2x3y_S_Py_S;
  Double I_ERI_H2x2yz_S_Py_Px = I_ERI_H2x2yz_S_Dxy_S+CDX*I_ERI_H2x2yz_S_Py_S;
  Double I_ERI_H2xy2z_S_Py_Px = I_ERI_H2xy2z_S_Dxy_S+CDX*I_ERI_H2xy2z_S_Py_S;
  Double I_ERI_H2x3z_S_Py_Px = I_ERI_H2x3z_S_Dxy_S+CDX*I_ERI_H2x3z_S_Py_S;
  Double I_ERI_Hx4y_S_Py_Px = I_ERI_Hx4y_S_Dxy_S+CDX*I_ERI_Hx4y_S_Py_S;
  Double I_ERI_Hx3yz_S_Py_Px = I_ERI_Hx3yz_S_Dxy_S+CDX*I_ERI_Hx3yz_S_Py_S;
  Double I_ERI_Hx2y2z_S_Py_Px = I_ERI_Hx2y2z_S_Dxy_S+CDX*I_ERI_Hx2y2z_S_Py_S;
  Double I_ERI_Hxy3z_S_Py_Px = I_ERI_Hxy3z_S_Dxy_S+CDX*I_ERI_Hxy3z_S_Py_S;
  Double I_ERI_Hx4z_S_Py_Px = I_ERI_Hx4z_S_Dxy_S+CDX*I_ERI_Hx4z_S_Py_S;
  Double I_ERI_H5y_S_Py_Px = I_ERI_H5y_S_Dxy_S+CDX*I_ERI_H5y_S_Py_S;
  Double I_ERI_H4yz_S_Py_Px = I_ERI_H4yz_S_Dxy_S+CDX*I_ERI_H4yz_S_Py_S;
  Double I_ERI_H3y2z_S_Py_Px = I_ERI_H3y2z_S_Dxy_S+CDX*I_ERI_H3y2z_S_Py_S;
  Double I_ERI_H2y3z_S_Py_Px = I_ERI_H2y3z_S_Dxy_S+CDX*I_ERI_H2y3z_S_Py_S;
  Double I_ERI_Hy4z_S_Py_Px = I_ERI_Hy4z_S_Dxy_S+CDX*I_ERI_Hy4z_S_Py_S;
  Double I_ERI_H5z_S_Py_Px = I_ERI_H5z_S_Dxy_S+CDX*I_ERI_H5z_S_Py_S;
  Double I_ERI_H5x_S_Pz_Px = I_ERI_H5x_S_Dxz_S+CDX*I_ERI_H5x_S_Pz_S;
  Double I_ERI_H4xy_S_Pz_Px = I_ERI_H4xy_S_Dxz_S+CDX*I_ERI_H4xy_S_Pz_S;
  Double I_ERI_H4xz_S_Pz_Px = I_ERI_H4xz_S_Dxz_S+CDX*I_ERI_H4xz_S_Pz_S;
  Double I_ERI_H3x2y_S_Pz_Px = I_ERI_H3x2y_S_Dxz_S+CDX*I_ERI_H3x2y_S_Pz_S;
  Double I_ERI_H3xyz_S_Pz_Px = I_ERI_H3xyz_S_Dxz_S+CDX*I_ERI_H3xyz_S_Pz_S;
  Double I_ERI_H3x2z_S_Pz_Px = I_ERI_H3x2z_S_Dxz_S+CDX*I_ERI_H3x2z_S_Pz_S;
  Double I_ERI_H2x3y_S_Pz_Px = I_ERI_H2x3y_S_Dxz_S+CDX*I_ERI_H2x3y_S_Pz_S;
  Double I_ERI_H2x2yz_S_Pz_Px = I_ERI_H2x2yz_S_Dxz_S+CDX*I_ERI_H2x2yz_S_Pz_S;
  Double I_ERI_H2xy2z_S_Pz_Px = I_ERI_H2xy2z_S_Dxz_S+CDX*I_ERI_H2xy2z_S_Pz_S;
  Double I_ERI_H2x3z_S_Pz_Px = I_ERI_H2x3z_S_Dxz_S+CDX*I_ERI_H2x3z_S_Pz_S;
  Double I_ERI_Hx4y_S_Pz_Px = I_ERI_Hx4y_S_Dxz_S+CDX*I_ERI_Hx4y_S_Pz_S;
  Double I_ERI_Hx3yz_S_Pz_Px = I_ERI_Hx3yz_S_Dxz_S+CDX*I_ERI_Hx3yz_S_Pz_S;
  Double I_ERI_Hx2y2z_S_Pz_Px = I_ERI_Hx2y2z_S_Dxz_S+CDX*I_ERI_Hx2y2z_S_Pz_S;
  Double I_ERI_Hxy3z_S_Pz_Px = I_ERI_Hxy3z_S_Dxz_S+CDX*I_ERI_Hxy3z_S_Pz_S;
  Double I_ERI_Hx4z_S_Pz_Px = I_ERI_Hx4z_S_Dxz_S+CDX*I_ERI_Hx4z_S_Pz_S;
  Double I_ERI_H5y_S_Pz_Px = I_ERI_H5y_S_Dxz_S+CDX*I_ERI_H5y_S_Pz_S;
  Double I_ERI_H4yz_S_Pz_Px = I_ERI_H4yz_S_Dxz_S+CDX*I_ERI_H4yz_S_Pz_S;
  Double I_ERI_H3y2z_S_Pz_Px = I_ERI_H3y2z_S_Dxz_S+CDX*I_ERI_H3y2z_S_Pz_S;
  Double I_ERI_H2y3z_S_Pz_Px = I_ERI_H2y3z_S_Dxz_S+CDX*I_ERI_H2y3z_S_Pz_S;
  Double I_ERI_Hy4z_S_Pz_Px = I_ERI_Hy4z_S_Dxz_S+CDX*I_ERI_Hy4z_S_Pz_S;
  Double I_ERI_H5z_S_Pz_Px = I_ERI_H5z_S_Dxz_S+CDX*I_ERI_H5z_S_Pz_S;
  Double I_ERI_H5x_S_Px_Py = I_ERI_H5x_S_Dxy_S+CDY*I_ERI_H5x_S_Px_S;
  Double I_ERI_H4xy_S_Px_Py = I_ERI_H4xy_S_Dxy_S+CDY*I_ERI_H4xy_S_Px_S;
  Double I_ERI_H4xz_S_Px_Py = I_ERI_H4xz_S_Dxy_S+CDY*I_ERI_H4xz_S_Px_S;
  Double I_ERI_H3x2y_S_Px_Py = I_ERI_H3x2y_S_Dxy_S+CDY*I_ERI_H3x2y_S_Px_S;
  Double I_ERI_H3xyz_S_Px_Py = I_ERI_H3xyz_S_Dxy_S+CDY*I_ERI_H3xyz_S_Px_S;
  Double I_ERI_H3x2z_S_Px_Py = I_ERI_H3x2z_S_Dxy_S+CDY*I_ERI_H3x2z_S_Px_S;
  Double I_ERI_H2x3y_S_Px_Py = I_ERI_H2x3y_S_Dxy_S+CDY*I_ERI_H2x3y_S_Px_S;
  Double I_ERI_H2x2yz_S_Px_Py = I_ERI_H2x2yz_S_Dxy_S+CDY*I_ERI_H2x2yz_S_Px_S;
  Double I_ERI_H2xy2z_S_Px_Py = I_ERI_H2xy2z_S_Dxy_S+CDY*I_ERI_H2xy2z_S_Px_S;
  Double I_ERI_H2x3z_S_Px_Py = I_ERI_H2x3z_S_Dxy_S+CDY*I_ERI_H2x3z_S_Px_S;
  Double I_ERI_Hx4y_S_Px_Py = I_ERI_Hx4y_S_Dxy_S+CDY*I_ERI_Hx4y_S_Px_S;
  Double I_ERI_Hx3yz_S_Px_Py = I_ERI_Hx3yz_S_Dxy_S+CDY*I_ERI_Hx3yz_S_Px_S;
  Double I_ERI_Hx2y2z_S_Px_Py = I_ERI_Hx2y2z_S_Dxy_S+CDY*I_ERI_Hx2y2z_S_Px_S;
  Double I_ERI_Hxy3z_S_Px_Py = I_ERI_Hxy3z_S_Dxy_S+CDY*I_ERI_Hxy3z_S_Px_S;
  Double I_ERI_Hx4z_S_Px_Py = I_ERI_Hx4z_S_Dxy_S+CDY*I_ERI_Hx4z_S_Px_S;
  Double I_ERI_H5y_S_Px_Py = I_ERI_H5y_S_Dxy_S+CDY*I_ERI_H5y_S_Px_S;
  Double I_ERI_H4yz_S_Px_Py = I_ERI_H4yz_S_Dxy_S+CDY*I_ERI_H4yz_S_Px_S;
  Double I_ERI_H3y2z_S_Px_Py = I_ERI_H3y2z_S_Dxy_S+CDY*I_ERI_H3y2z_S_Px_S;
  Double I_ERI_H2y3z_S_Px_Py = I_ERI_H2y3z_S_Dxy_S+CDY*I_ERI_H2y3z_S_Px_S;
  Double I_ERI_Hy4z_S_Px_Py = I_ERI_Hy4z_S_Dxy_S+CDY*I_ERI_Hy4z_S_Px_S;
  Double I_ERI_H5z_S_Px_Py = I_ERI_H5z_S_Dxy_S+CDY*I_ERI_H5z_S_Px_S;
  Double I_ERI_H5x_S_Py_Py = I_ERI_H5x_S_D2y_S+CDY*I_ERI_H5x_S_Py_S;
  Double I_ERI_H4xy_S_Py_Py = I_ERI_H4xy_S_D2y_S+CDY*I_ERI_H4xy_S_Py_S;
  Double I_ERI_H4xz_S_Py_Py = I_ERI_H4xz_S_D2y_S+CDY*I_ERI_H4xz_S_Py_S;
  Double I_ERI_H3x2y_S_Py_Py = I_ERI_H3x2y_S_D2y_S+CDY*I_ERI_H3x2y_S_Py_S;
  Double I_ERI_H3xyz_S_Py_Py = I_ERI_H3xyz_S_D2y_S+CDY*I_ERI_H3xyz_S_Py_S;
  Double I_ERI_H3x2z_S_Py_Py = I_ERI_H3x2z_S_D2y_S+CDY*I_ERI_H3x2z_S_Py_S;
  Double I_ERI_H2x3y_S_Py_Py = I_ERI_H2x3y_S_D2y_S+CDY*I_ERI_H2x3y_S_Py_S;
  Double I_ERI_H2x2yz_S_Py_Py = I_ERI_H2x2yz_S_D2y_S+CDY*I_ERI_H2x2yz_S_Py_S;
  Double I_ERI_H2xy2z_S_Py_Py = I_ERI_H2xy2z_S_D2y_S+CDY*I_ERI_H2xy2z_S_Py_S;
  Double I_ERI_H2x3z_S_Py_Py = I_ERI_H2x3z_S_D2y_S+CDY*I_ERI_H2x3z_S_Py_S;
  Double I_ERI_Hx4y_S_Py_Py = I_ERI_Hx4y_S_D2y_S+CDY*I_ERI_Hx4y_S_Py_S;
  Double I_ERI_Hx3yz_S_Py_Py = I_ERI_Hx3yz_S_D2y_S+CDY*I_ERI_Hx3yz_S_Py_S;
  Double I_ERI_Hx2y2z_S_Py_Py = I_ERI_Hx2y2z_S_D2y_S+CDY*I_ERI_Hx2y2z_S_Py_S;
  Double I_ERI_Hxy3z_S_Py_Py = I_ERI_Hxy3z_S_D2y_S+CDY*I_ERI_Hxy3z_S_Py_S;
  Double I_ERI_Hx4z_S_Py_Py = I_ERI_Hx4z_S_D2y_S+CDY*I_ERI_Hx4z_S_Py_S;
  Double I_ERI_H5y_S_Py_Py = I_ERI_H5y_S_D2y_S+CDY*I_ERI_H5y_S_Py_S;
  Double I_ERI_H4yz_S_Py_Py = I_ERI_H4yz_S_D2y_S+CDY*I_ERI_H4yz_S_Py_S;
  Double I_ERI_H3y2z_S_Py_Py = I_ERI_H3y2z_S_D2y_S+CDY*I_ERI_H3y2z_S_Py_S;
  Double I_ERI_H2y3z_S_Py_Py = I_ERI_H2y3z_S_D2y_S+CDY*I_ERI_H2y3z_S_Py_S;
  Double I_ERI_Hy4z_S_Py_Py = I_ERI_Hy4z_S_D2y_S+CDY*I_ERI_Hy4z_S_Py_S;
  Double I_ERI_H5z_S_Py_Py = I_ERI_H5z_S_D2y_S+CDY*I_ERI_H5z_S_Py_S;
  Double I_ERI_H5x_S_Pz_Py = I_ERI_H5x_S_Dyz_S+CDY*I_ERI_H5x_S_Pz_S;
  Double I_ERI_H4xy_S_Pz_Py = I_ERI_H4xy_S_Dyz_S+CDY*I_ERI_H4xy_S_Pz_S;
  Double I_ERI_H4xz_S_Pz_Py = I_ERI_H4xz_S_Dyz_S+CDY*I_ERI_H4xz_S_Pz_S;
  Double I_ERI_H3x2y_S_Pz_Py = I_ERI_H3x2y_S_Dyz_S+CDY*I_ERI_H3x2y_S_Pz_S;
  Double I_ERI_H3xyz_S_Pz_Py = I_ERI_H3xyz_S_Dyz_S+CDY*I_ERI_H3xyz_S_Pz_S;
  Double I_ERI_H3x2z_S_Pz_Py = I_ERI_H3x2z_S_Dyz_S+CDY*I_ERI_H3x2z_S_Pz_S;
  Double I_ERI_H2x3y_S_Pz_Py = I_ERI_H2x3y_S_Dyz_S+CDY*I_ERI_H2x3y_S_Pz_S;
  Double I_ERI_H2x2yz_S_Pz_Py = I_ERI_H2x2yz_S_Dyz_S+CDY*I_ERI_H2x2yz_S_Pz_S;
  Double I_ERI_H2xy2z_S_Pz_Py = I_ERI_H2xy2z_S_Dyz_S+CDY*I_ERI_H2xy2z_S_Pz_S;
  Double I_ERI_H2x3z_S_Pz_Py = I_ERI_H2x3z_S_Dyz_S+CDY*I_ERI_H2x3z_S_Pz_S;
  Double I_ERI_Hx4y_S_Pz_Py = I_ERI_Hx4y_S_Dyz_S+CDY*I_ERI_Hx4y_S_Pz_S;
  Double I_ERI_Hx3yz_S_Pz_Py = I_ERI_Hx3yz_S_Dyz_S+CDY*I_ERI_Hx3yz_S_Pz_S;
  Double I_ERI_Hx2y2z_S_Pz_Py = I_ERI_Hx2y2z_S_Dyz_S+CDY*I_ERI_Hx2y2z_S_Pz_S;
  Double I_ERI_Hxy3z_S_Pz_Py = I_ERI_Hxy3z_S_Dyz_S+CDY*I_ERI_Hxy3z_S_Pz_S;
  Double I_ERI_Hx4z_S_Pz_Py = I_ERI_Hx4z_S_Dyz_S+CDY*I_ERI_Hx4z_S_Pz_S;
  Double I_ERI_H5y_S_Pz_Py = I_ERI_H5y_S_Dyz_S+CDY*I_ERI_H5y_S_Pz_S;
  Double I_ERI_H4yz_S_Pz_Py = I_ERI_H4yz_S_Dyz_S+CDY*I_ERI_H4yz_S_Pz_S;
  Double I_ERI_H3y2z_S_Pz_Py = I_ERI_H3y2z_S_Dyz_S+CDY*I_ERI_H3y2z_S_Pz_S;
  Double I_ERI_H2y3z_S_Pz_Py = I_ERI_H2y3z_S_Dyz_S+CDY*I_ERI_H2y3z_S_Pz_S;
  Double I_ERI_Hy4z_S_Pz_Py = I_ERI_Hy4z_S_Dyz_S+CDY*I_ERI_Hy4z_S_Pz_S;
  Double I_ERI_H5z_S_Pz_Py = I_ERI_H5z_S_Dyz_S+CDY*I_ERI_H5z_S_Pz_S;
  Double I_ERI_H5x_S_Px_Pz = I_ERI_H5x_S_Dxz_S+CDZ*I_ERI_H5x_S_Px_S;
  Double I_ERI_H4xy_S_Px_Pz = I_ERI_H4xy_S_Dxz_S+CDZ*I_ERI_H4xy_S_Px_S;
  Double I_ERI_H4xz_S_Px_Pz = I_ERI_H4xz_S_Dxz_S+CDZ*I_ERI_H4xz_S_Px_S;
  Double I_ERI_H3x2y_S_Px_Pz = I_ERI_H3x2y_S_Dxz_S+CDZ*I_ERI_H3x2y_S_Px_S;
  Double I_ERI_H3xyz_S_Px_Pz = I_ERI_H3xyz_S_Dxz_S+CDZ*I_ERI_H3xyz_S_Px_S;
  Double I_ERI_H3x2z_S_Px_Pz = I_ERI_H3x2z_S_Dxz_S+CDZ*I_ERI_H3x2z_S_Px_S;
  Double I_ERI_H2x3y_S_Px_Pz = I_ERI_H2x3y_S_Dxz_S+CDZ*I_ERI_H2x3y_S_Px_S;
  Double I_ERI_H2x2yz_S_Px_Pz = I_ERI_H2x2yz_S_Dxz_S+CDZ*I_ERI_H2x2yz_S_Px_S;
  Double I_ERI_H2xy2z_S_Px_Pz = I_ERI_H2xy2z_S_Dxz_S+CDZ*I_ERI_H2xy2z_S_Px_S;
  Double I_ERI_H2x3z_S_Px_Pz = I_ERI_H2x3z_S_Dxz_S+CDZ*I_ERI_H2x3z_S_Px_S;
  Double I_ERI_Hx4y_S_Px_Pz = I_ERI_Hx4y_S_Dxz_S+CDZ*I_ERI_Hx4y_S_Px_S;
  Double I_ERI_Hx3yz_S_Px_Pz = I_ERI_Hx3yz_S_Dxz_S+CDZ*I_ERI_Hx3yz_S_Px_S;
  Double I_ERI_Hx2y2z_S_Px_Pz = I_ERI_Hx2y2z_S_Dxz_S+CDZ*I_ERI_Hx2y2z_S_Px_S;
  Double I_ERI_Hxy3z_S_Px_Pz = I_ERI_Hxy3z_S_Dxz_S+CDZ*I_ERI_Hxy3z_S_Px_S;
  Double I_ERI_Hx4z_S_Px_Pz = I_ERI_Hx4z_S_Dxz_S+CDZ*I_ERI_Hx4z_S_Px_S;
  Double I_ERI_H5y_S_Px_Pz = I_ERI_H5y_S_Dxz_S+CDZ*I_ERI_H5y_S_Px_S;
  Double I_ERI_H4yz_S_Px_Pz = I_ERI_H4yz_S_Dxz_S+CDZ*I_ERI_H4yz_S_Px_S;
  Double I_ERI_H3y2z_S_Px_Pz = I_ERI_H3y2z_S_Dxz_S+CDZ*I_ERI_H3y2z_S_Px_S;
  Double I_ERI_H2y3z_S_Px_Pz = I_ERI_H2y3z_S_Dxz_S+CDZ*I_ERI_H2y3z_S_Px_S;
  Double I_ERI_Hy4z_S_Px_Pz = I_ERI_Hy4z_S_Dxz_S+CDZ*I_ERI_Hy4z_S_Px_S;
  Double I_ERI_H5z_S_Px_Pz = I_ERI_H5z_S_Dxz_S+CDZ*I_ERI_H5z_S_Px_S;
  Double I_ERI_H5x_S_Py_Pz = I_ERI_H5x_S_Dyz_S+CDZ*I_ERI_H5x_S_Py_S;
  Double I_ERI_H4xy_S_Py_Pz = I_ERI_H4xy_S_Dyz_S+CDZ*I_ERI_H4xy_S_Py_S;
  Double I_ERI_H4xz_S_Py_Pz = I_ERI_H4xz_S_Dyz_S+CDZ*I_ERI_H4xz_S_Py_S;
  Double I_ERI_H3x2y_S_Py_Pz = I_ERI_H3x2y_S_Dyz_S+CDZ*I_ERI_H3x2y_S_Py_S;
  Double I_ERI_H3xyz_S_Py_Pz = I_ERI_H3xyz_S_Dyz_S+CDZ*I_ERI_H3xyz_S_Py_S;
  Double I_ERI_H3x2z_S_Py_Pz = I_ERI_H3x2z_S_Dyz_S+CDZ*I_ERI_H3x2z_S_Py_S;
  Double I_ERI_H2x3y_S_Py_Pz = I_ERI_H2x3y_S_Dyz_S+CDZ*I_ERI_H2x3y_S_Py_S;
  Double I_ERI_H2x2yz_S_Py_Pz = I_ERI_H2x2yz_S_Dyz_S+CDZ*I_ERI_H2x2yz_S_Py_S;
  Double I_ERI_H2xy2z_S_Py_Pz = I_ERI_H2xy2z_S_Dyz_S+CDZ*I_ERI_H2xy2z_S_Py_S;
  Double I_ERI_H2x3z_S_Py_Pz = I_ERI_H2x3z_S_Dyz_S+CDZ*I_ERI_H2x3z_S_Py_S;
  Double I_ERI_Hx4y_S_Py_Pz = I_ERI_Hx4y_S_Dyz_S+CDZ*I_ERI_Hx4y_S_Py_S;
  Double I_ERI_Hx3yz_S_Py_Pz = I_ERI_Hx3yz_S_Dyz_S+CDZ*I_ERI_Hx3yz_S_Py_S;
  Double I_ERI_Hx2y2z_S_Py_Pz = I_ERI_Hx2y2z_S_Dyz_S+CDZ*I_ERI_Hx2y2z_S_Py_S;
  Double I_ERI_Hxy3z_S_Py_Pz = I_ERI_Hxy3z_S_Dyz_S+CDZ*I_ERI_Hxy3z_S_Py_S;
  Double I_ERI_Hx4z_S_Py_Pz = I_ERI_Hx4z_S_Dyz_S+CDZ*I_ERI_Hx4z_S_Py_S;
  Double I_ERI_H5y_S_Py_Pz = I_ERI_H5y_S_Dyz_S+CDZ*I_ERI_H5y_S_Py_S;
  Double I_ERI_H4yz_S_Py_Pz = I_ERI_H4yz_S_Dyz_S+CDZ*I_ERI_H4yz_S_Py_S;
  Double I_ERI_H3y2z_S_Py_Pz = I_ERI_H3y2z_S_Dyz_S+CDZ*I_ERI_H3y2z_S_Py_S;
  Double I_ERI_H2y3z_S_Py_Pz = I_ERI_H2y3z_S_Dyz_S+CDZ*I_ERI_H2y3z_S_Py_S;
  Double I_ERI_Hy4z_S_Py_Pz = I_ERI_Hy4z_S_Dyz_S+CDZ*I_ERI_Hy4z_S_Py_S;
  Double I_ERI_H5z_S_Py_Pz = I_ERI_H5z_S_Dyz_S+CDZ*I_ERI_H5z_S_Py_S;
  Double I_ERI_H5x_S_Pz_Pz = I_ERI_H5x_S_D2z_S+CDZ*I_ERI_H5x_S_Pz_S;
  Double I_ERI_H4xy_S_Pz_Pz = I_ERI_H4xy_S_D2z_S+CDZ*I_ERI_H4xy_S_Pz_S;
  Double I_ERI_H4xz_S_Pz_Pz = I_ERI_H4xz_S_D2z_S+CDZ*I_ERI_H4xz_S_Pz_S;
  Double I_ERI_H3x2y_S_Pz_Pz = I_ERI_H3x2y_S_D2z_S+CDZ*I_ERI_H3x2y_S_Pz_S;
  Double I_ERI_H3xyz_S_Pz_Pz = I_ERI_H3xyz_S_D2z_S+CDZ*I_ERI_H3xyz_S_Pz_S;
  Double I_ERI_H3x2z_S_Pz_Pz = I_ERI_H3x2z_S_D2z_S+CDZ*I_ERI_H3x2z_S_Pz_S;
  Double I_ERI_H2x3y_S_Pz_Pz = I_ERI_H2x3y_S_D2z_S+CDZ*I_ERI_H2x3y_S_Pz_S;
  Double I_ERI_H2x2yz_S_Pz_Pz = I_ERI_H2x2yz_S_D2z_S+CDZ*I_ERI_H2x2yz_S_Pz_S;
  Double I_ERI_H2xy2z_S_Pz_Pz = I_ERI_H2xy2z_S_D2z_S+CDZ*I_ERI_H2xy2z_S_Pz_S;
  Double I_ERI_H2x3z_S_Pz_Pz = I_ERI_H2x3z_S_D2z_S+CDZ*I_ERI_H2x3z_S_Pz_S;
  Double I_ERI_Hx4y_S_Pz_Pz = I_ERI_Hx4y_S_D2z_S+CDZ*I_ERI_Hx4y_S_Pz_S;
  Double I_ERI_Hx3yz_S_Pz_Pz = I_ERI_Hx3yz_S_D2z_S+CDZ*I_ERI_Hx3yz_S_Pz_S;
  Double I_ERI_Hx2y2z_S_Pz_Pz = I_ERI_Hx2y2z_S_D2z_S+CDZ*I_ERI_Hx2y2z_S_Pz_S;
  Double I_ERI_Hxy3z_S_Pz_Pz = I_ERI_Hxy3z_S_D2z_S+CDZ*I_ERI_Hxy3z_S_Pz_S;
  Double I_ERI_Hx4z_S_Pz_Pz = I_ERI_Hx4z_S_D2z_S+CDZ*I_ERI_Hx4z_S_Pz_S;
  Double I_ERI_H5y_S_Pz_Pz = I_ERI_H5y_S_D2z_S+CDZ*I_ERI_H5y_S_Pz_S;
  Double I_ERI_H4yz_S_Pz_Pz = I_ERI_H4yz_S_D2z_S+CDZ*I_ERI_H4yz_S_Pz_S;
  Double I_ERI_H3y2z_S_Pz_Pz = I_ERI_H3y2z_S_D2z_S+CDZ*I_ERI_H3y2z_S_Pz_S;
  Double I_ERI_H2y3z_S_Pz_Pz = I_ERI_H2y3z_S_D2z_S+CDZ*I_ERI_H2y3z_S_Pz_S;
  Double I_ERI_Hy4z_S_Pz_Pz = I_ERI_Hy4z_S_D2z_S+CDZ*I_ERI_Hy4z_S_Pz_S;
  Double I_ERI_H5z_S_Pz_Pz = I_ERI_H5z_S_D2z_S+CDZ*I_ERI_H5z_S_Pz_S;

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
   * shell quartet name: SQ_ERI_F_P_P_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_S_P_P
   * RHS shell quartet name: SQ_ERI_F_S_P_P
   ************************************************************/
  Double I_ERI_F3x_Px_Px_Px = I_ERI_G4x_S_Px_Px+ABX*I_ERI_F3x_S_Px_Px;
  Double I_ERI_F2xy_Px_Px_Px = I_ERI_G3xy_S_Px_Px+ABX*I_ERI_F2xy_S_Px_Px;
  Double I_ERI_F2xz_Px_Px_Px = I_ERI_G3xz_S_Px_Px+ABX*I_ERI_F2xz_S_Px_Px;
  Double I_ERI_Fx2y_Px_Px_Px = I_ERI_G2x2y_S_Px_Px+ABX*I_ERI_Fx2y_S_Px_Px;
  Double I_ERI_Fxyz_Px_Px_Px = I_ERI_G2xyz_S_Px_Px+ABX*I_ERI_Fxyz_S_Px_Px;
  Double I_ERI_Fx2z_Px_Px_Px = I_ERI_G2x2z_S_Px_Px+ABX*I_ERI_Fx2z_S_Px_Px;
  Double I_ERI_F3y_Px_Px_Px = I_ERI_Gx3y_S_Px_Px+ABX*I_ERI_F3y_S_Px_Px;
  Double I_ERI_F2yz_Px_Px_Px = I_ERI_Gx2yz_S_Px_Px+ABX*I_ERI_F2yz_S_Px_Px;
  Double I_ERI_Fy2z_Px_Px_Px = I_ERI_Gxy2z_S_Px_Px+ABX*I_ERI_Fy2z_S_Px_Px;
  Double I_ERI_F3z_Px_Px_Px = I_ERI_Gx3z_S_Px_Px+ABX*I_ERI_F3z_S_Px_Px;
  Double I_ERI_F3x_Py_Px_Px = I_ERI_G3xy_S_Px_Px+ABY*I_ERI_F3x_S_Px_Px;
  Double I_ERI_F2xy_Py_Px_Px = I_ERI_G2x2y_S_Px_Px+ABY*I_ERI_F2xy_S_Px_Px;
  Double I_ERI_F2xz_Py_Px_Px = I_ERI_G2xyz_S_Px_Px+ABY*I_ERI_F2xz_S_Px_Px;
  Double I_ERI_Fx2y_Py_Px_Px = I_ERI_Gx3y_S_Px_Px+ABY*I_ERI_Fx2y_S_Px_Px;
  Double I_ERI_Fxyz_Py_Px_Px = I_ERI_Gx2yz_S_Px_Px+ABY*I_ERI_Fxyz_S_Px_Px;
  Double I_ERI_Fx2z_Py_Px_Px = I_ERI_Gxy2z_S_Px_Px+ABY*I_ERI_Fx2z_S_Px_Px;
  Double I_ERI_F3y_Py_Px_Px = I_ERI_G4y_S_Px_Px+ABY*I_ERI_F3y_S_Px_Px;
  Double I_ERI_F2yz_Py_Px_Px = I_ERI_G3yz_S_Px_Px+ABY*I_ERI_F2yz_S_Px_Px;
  Double I_ERI_Fy2z_Py_Px_Px = I_ERI_G2y2z_S_Px_Px+ABY*I_ERI_Fy2z_S_Px_Px;
  Double I_ERI_F3z_Py_Px_Px = I_ERI_Gy3z_S_Px_Px+ABY*I_ERI_F3z_S_Px_Px;
  Double I_ERI_F3x_Pz_Px_Px = I_ERI_G3xz_S_Px_Px+ABZ*I_ERI_F3x_S_Px_Px;
  Double I_ERI_F2xy_Pz_Px_Px = I_ERI_G2xyz_S_Px_Px+ABZ*I_ERI_F2xy_S_Px_Px;
  Double I_ERI_F2xz_Pz_Px_Px = I_ERI_G2x2z_S_Px_Px+ABZ*I_ERI_F2xz_S_Px_Px;
  Double I_ERI_Fx2y_Pz_Px_Px = I_ERI_Gx2yz_S_Px_Px+ABZ*I_ERI_Fx2y_S_Px_Px;
  Double I_ERI_Fxyz_Pz_Px_Px = I_ERI_Gxy2z_S_Px_Px+ABZ*I_ERI_Fxyz_S_Px_Px;
  Double I_ERI_Fx2z_Pz_Px_Px = I_ERI_Gx3z_S_Px_Px+ABZ*I_ERI_Fx2z_S_Px_Px;
  Double I_ERI_F3y_Pz_Px_Px = I_ERI_G3yz_S_Px_Px+ABZ*I_ERI_F3y_S_Px_Px;
  Double I_ERI_F2yz_Pz_Px_Px = I_ERI_G2y2z_S_Px_Px+ABZ*I_ERI_F2yz_S_Px_Px;
  Double I_ERI_Fy2z_Pz_Px_Px = I_ERI_Gy3z_S_Px_Px+ABZ*I_ERI_Fy2z_S_Px_Px;
  Double I_ERI_F3z_Pz_Px_Px = I_ERI_G4z_S_Px_Px+ABZ*I_ERI_F3z_S_Px_Px;
  Double I_ERI_F3x_Px_Py_Px = I_ERI_G4x_S_Py_Px+ABX*I_ERI_F3x_S_Py_Px;
  Double I_ERI_F2xy_Px_Py_Px = I_ERI_G3xy_S_Py_Px+ABX*I_ERI_F2xy_S_Py_Px;
  Double I_ERI_F2xz_Px_Py_Px = I_ERI_G3xz_S_Py_Px+ABX*I_ERI_F2xz_S_Py_Px;
  Double I_ERI_Fx2y_Px_Py_Px = I_ERI_G2x2y_S_Py_Px+ABX*I_ERI_Fx2y_S_Py_Px;
  Double I_ERI_Fxyz_Px_Py_Px = I_ERI_G2xyz_S_Py_Px+ABX*I_ERI_Fxyz_S_Py_Px;
  Double I_ERI_Fx2z_Px_Py_Px = I_ERI_G2x2z_S_Py_Px+ABX*I_ERI_Fx2z_S_Py_Px;
  Double I_ERI_F3y_Px_Py_Px = I_ERI_Gx3y_S_Py_Px+ABX*I_ERI_F3y_S_Py_Px;
  Double I_ERI_F2yz_Px_Py_Px = I_ERI_Gx2yz_S_Py_Px+ABX*I_ERI_F2yz_S_Py_Px;
  Double I_ERI_Fy2z_Px_Py_Px = I_ERI_Gxy2z_S_Py_Px+ABX*I_ERI_Fy2z_S_Py_Px;
  Double I_ERI_F3z_Px_Py_Px = I_ERI_Gx3z_S_Py_Px+ABX*I_ERI_F3z_S_Py_Px;
  Double I_ERI_F3x_Py_Py_Px = I_ERI_G3xy_S_Py_Px+ABY*I_ERI_F3x_S_Py_Px;
  Double I_ERI_F2xy_Py_Py_Px = I_ERI_G2x2y_S_Py_Px+ABY*I_ERI_F2xy_S_Py_Px;
  Double I_ERI_F2xz_Py_Py_Px = I_ERI_G2xyz_S_Py_Px+ABY*I_ERI_F2xz_S_Py_Px;
  Double I_ERI_Fx2y_Py_Py_Px = I_ERI_Gx3y_S_Py_Px+ABY*I_ERI_Fx2y_S_Py_Px;
  Double I_ERI_Fxyz_Py_Py_Px = I_ERI_Gx2yz_S_Py_Px+ABY*I_ERI_Fxyz_S_Py_Px;
  Double I_ERI_Fx2z_Py_Py_Px = I_ERI_Gxy2z_S_Py_Px+ABY*I_ERI_Fx2z_S_Py_Px;
  Double I_ERI_F3y_Py_Py_Px = I_ERI_G4y_S_Py_Px+ABY*I_ERI_F3y_S_Py_Px;
  Double I_ERI_F2yz_Py_Py_Px = I_ERI_G3yz_S_Py_Px+ABY*I_ERI_F2yz_S_Py_Px;
  Double I_ERI_Fy2z_Py_Py_Px = I_ERI_G2y2z_S_Py_Px+ABY*I_ERI_Fy2z_S_Py_Px;
  Double I_ERI_F3z_Py_Py_Px = I_ERI_Gy3z_S_Py_Px+ABY*I_ERI_F3z_S_Py_Px;
  Double I_ERI_F3x_Pz_Py_Px = I_ERI_G3xz_S_Py_Px+ABZ*I_ERI_F3x_S_Py_Px;
  Double I_ERI_F2xy_Pz_Py_Px = I_ERI_G2xyz_S_Py_Px+ABZ*I_ERI_F2xy_S_Py_Px;
  Double I_ERI_F2xz_Pz_Py_Px = I_ERI_G2x2z_S_Py_Px+ABZ*I_ERI_F2xz_S_Py_Px;
  Double I_ERI_Fx2y_Pz_Py_Px = I_ERI_Gx2yz_S_Py_Px+ABZ*I_ERI_Fx2y_S_Py_Px;
  Double I_ERI_Fxyz_Pz_Py_Px = I_ERI_Gxy2z_S_Py_Px+ABZ*I_ERI_Fxyz_S_Py_Px;
  Double I_ERI_Fx2z_Pz_Py_Px = I_ERI_Gx3z_S_Py_Px+ABZ*I_ERI_Fx2z_S_Py_Px;
  Double I_ERI_F3y_Pz_Py_Px = I_ERI_G3yz_S_Py_Px+ABZ*I_ERI_F3y_S_Py_Px;
  Double I_ERI_F2yz_Pz_Py_Px = I_ERI_G2y2z_S_Py_Px+ABZ*I_ERI_F2yz_S_Py_Px;
  Double I_ERI_Fy2z_Pz_Py_Px = I_ERI_Gy3z_S_Py_Px+ABZ*I_ERI_Fy2z_S_Py_Px;
  Double I_ERI_F3z_Pz_Py_Px = I_ERI_G4z_S_Py_Px+ABZ*I_ERI_F3z_S_Py_Px;
  Double I_ERI_F3x_Px_Pz_Px = I_ERI_G4x_S_Pz_Px+ABX*I_ERI_F3x_S_Pz_Px;
  Double I_ERI_F2xy_Px_Pz_Px = I_ERI_G3xy_S_Pz_Px+ABX*I_ERI_F2xy_S_Pz_Px;
  Double I_ERI_F2xz_Px_Pz_Px = I_ERI_G3xz_S_Pz_Px+ABX*I_ERI_F2xz_S_Pz_Px;
  Double I_ERI_Fx2y_Px_Pz_Px = I_ERI_G2x2y_S_Pz_Px+ABX*I_ERI_Fx2y_S_Pz_Px;
  Double I_ERI_Fxyz_Px_Pz_Px = I_ERI_G2xyz_S_Pz_Px+ABX*I_ERI_Fxyz_S_Pz_Px;
  Double I_ERI_Fx2z_Px_Pz_Px = I_ERI_G2x2z_S_Pz_Px+ABX*I_ERI_Fx2z_S_Pz_Px;
  Double I_ERI_F3y_Px_Pz_Px = I_ERI_Gx3y_S_Pz_Px+ABX*I_ERI_F3y_S_Pz_Px;
  Double I_ERI_F2yz_Px_Pz_Px = I_ERI_Gx2yz_S_Pz_Px+ABX*I_ERI_F2yz_S_Pz_Px;
  Double I_ERI_Fy2z_Px_Pz_Px = I_ERI_Gxy2z_S_Pz_Px+ABX*I_ERI_Fy2z_S_Pz_Px;
  Double I_ERI_F3z_Px_Pz_Px = I_ERI_Gx3z_S_Pz_Px+ABX*I_ERI_F3z_S_Pz_Px;
  Double I_ERI_F3x_Py_Pz_Px = I_ERI_G3xy_S_Pz_Px+ABY*I_ERI_F3x_S_Pz_Px;
  Double I_ERI_F2xy_Py_Pz_Px = I_ERI_G2x2y_S_Pz_Px+ABY*I_ERI_F2xy_S_Pz_Px;
  Double I_ERI_F2xz_Py_Pz_Px = I_ERI_G2xyz_S_Pz_Px+ABY*I_ERI_F2xz_S_Pz_Px;
  Double I_ERI_Fx2y_Py_Pz_Px = I_ERI_Gx3y_S_Pz_Px+ABY*I_ERI_Fx2y_S_Pz_Px;
  Double I_ERI_Fxyz_Py_Pz_Px = I_ERI_Gx2yz_S_Pz_Px+ABY*I_ERI_Fxyz_S_Pz_Px;
  Double I_ERI_Fx2z_Py_Pz_Px = I_ERI_Gxy2z_S_Pz_Px+ABY*I_ERI_Fx2z_S_Pz_Px;
  Double I_ERI_F3y_Py_Pz_Px = I_ERI_G4y_S_Pz_Px+ABY*I_ERI_F3y_S_Pz_Px;
  Double I_ERI_F2yz_Py_Pz_Px = I_ERI_G3yz_S_Pz_Px+ABY*I_ERI_F2yz_S_Pz_Px;
  Double I_ERI_Fy2z_Py_Pz_Px = I_ERI_G2y2z_S_Pz_Px+ABY*I_ERI_Fy2z_S_Pz_Px;
  Double I_ERI_F3z_Py_Pz_Px = I_ERI_Gy3z_S_Pz_Px+ABY*I_ERI_F3z_S_Pz_Px;
  Double I_ERI_F3x_Pz_Pz_Px = I_ERI_G3xz_S_Pz_Px+ABZ*I_ERI_F3x_S_Pz_Px;
  Double I_ERI_F2xy_Pz_Pz_Px = I_ERI_G2xyz_S_Pz_Px+ABZ*I_ERI_F2xy_S_Pz_Px;
  Double I_ERI_F2xz_Pz_Pz_Px = I_ERI_G2x2z_S_Pz_Px+ABZ*I_ERI_F2xz_S_Pz_Px;
  Double I_ERI_Fx2y_Pz_Pz_Px = I_ERI_Gx2yz_S_Pz_Px+ABZ*I_ERI_Fx2y_S_Pz_Px;
  Double I_ERI_Fxyz_Pz_Pz_Px = I_ERI_Gxy2z_S_Pz_Px+ABZ*I_ERI_Fxyz_S_Pz_Px;
  Double I_ERI_Fx2z_Pz_Pz_Px = I_ERI_Gx3z_S_Pz_Px+ABZ*I_ERI_Fx2z_S_Pz_Px;
  Double I_ERI_F3y_Pz_Pz_Px = I_ERI_G3yz_S_Pz_Px+ABZ*I_ERI_F3y_S_Pz_Px;
  Double I_ERI_F2yz_Pz_Pz_Px = I_ERI_G2y2z_S_Pz_Px+ABZ*I_ERI_F2yz_S_Pz_Px;
  Double I_ERI_Fy2z_Pz_Pz_Px = I_ERI_Gy3z_S_Pz_Px+ABZ*I_ERI_Fy2z_S_Pz_Px;
  Double I_ERI_F3z_Pz_Pz_Px = I_ERI_G4z_S_Pz_Px+ABZ*I_ERI_F3z_S_Pz_Px;
  Double I_ERI_F3x_Px_Px_Py = I_ERI_G4x_S_Px_Py+ABX*I_ERI_F3x_S_Px_Py;
  Double I_ERI_F2xy_Px_Px_Py = I_ERI_G3xy_S_Px_Py+ABX*I_ERI_F2xy_S_Px_Py;
  Double I_ERI_F2xz_Px_Px_Py = I_ERI_G3xz_S_Px_Py+ABX*I_ERI_F2xz_S_Px_Py;
  Double I_ERI_Fx2y_Px_Px_Py = I_ERI_G2x2y_S_Px_Py+ABX*I_ERI_Fx2y_S_Px_Py;
  Double I_ERI_Fxyz_Px_Px_Py = I_ERI_G2xyz_S_Px_Py+ABX*I_ERI_Fxyz_S_Px_Py;
  Double I_ERI_Fx2z_Px_Px_Py = I_ERI_G2x2z_S_Px_Py+ABX*I_ERI_Fx2z_S_Px_Py;
  Double I_ERI_F3y_Px_Px_Py = I_ERI_Gx3y_S_Px_Py+ABX*I_ERI_F3y_S_Px_Py;
  Double I_ERI_F2yz_Px_Px_Py = I_ERI_Gx2yz_S_Px_Py+ABX*I_ERI_F2yz_S_Px_Py;
  Double I_ERI_Fy2z_Px_Px_Py = I_ERI_Gxy2z_S_Px_Py+ABX*I_ERI_Fy2z_S_Px_Py;
  Double I_ERI_F3z_Px_Px_Py = I_ERI_Gx3z_S_Px_Py+ABX*I_ERI_F3z_S_Px_Py;
  Double I_ERI_F3x_Py_Px_Py = I_ERI_G3xy_S_Px_Py+ABY*I_ERI_F3x_S_Px_Py;
  Double I_ERI_F2xy_Py_Px_Py = I_ERI_G2x2y_S_Px_Py+ABY*I_ERI_F2xy_S_Px_Py;
  Double I_ERI_F2xz_Py_Px_Py = I_ERI_G2xyz_S_Px_Py+ABY*I_ERI_F2xz_S_Px_Py;
  Double I_ERI_Fx2y_Py_Px_Py = I_ERI_Gx3y_S_Px_Py+ABY*I_ERI_Fx2y_S_Px_Py;
  Double I_ERI_Fxyz_Py_Px_Py = I_ERI_Gx2yz_S_Px_Py+ABY*I_ERI_Fxyz_S_Px_Py;
  Double I_ERI_Fx2z_Py_Px_Py = I_ERI_Gxy2z_S_Px_Py+ABY*I_ERI_Fx2z_S_Px_Py;
  Double I_ERI_F3y_Py_Px_Py = I_ERI_G4y_S_Px_Py+ABY*I_ERI_F3y_S_Px_Py;
  Double I_ERI_F2yz_Py_Px_Py = I_ERI_G3yz_S_Px_Py+ABY*I_ERI_F2yz_S_Px_Py;
  Double I_ERI_Fy2z_Py_Px_Py = I_ERI_G2y2z_S_Px_Py+ABY*I_ERI_Fy2z_S_Px_Py;
  Double I_ERI_F3z_Py_Px_Py = I_ERI_Gy3z_S_Px_Py+ABY*I_ERI_F3z_S_Px_Py;
  Double I_ERI_F3x_Pz_Px_Py = I_ERI_G3xz_S_Px_Py+ABZ*I_ERI_F3x_S_Px_Py;
  Double I_ERI_F2xy_Pz_Px_Py = I_ERI_G2xyz_S_Px_Py+ABZ*I_ERI_F2xy_S_Px_Py;
  Double I_ERI_F2xz_Pz_Px_Py = I_ERI_G2x2z_S_Px_Py+ABZ*I_ERI_F2xz_S_Px_Py;
  Double I_ERI_Fx2y_Pz_Px_Py = I_ERI_Gx2yz_S_Px_Py+ABZ*I_ERI_Fx2y_S_Px_Py;
  Double I_ERI_Fxyz_Pz_Px_Py = I_ERI_Gxy2z_S_Px_Py+ABZ*I_ERI_Fxyz_S_Px_Py;
  Double I_ERI_Fx2z_Pz_Px_Py = I_ERI_Gx3z_S_Px_Py+ABZ*I_ERI_Fx2z_S_Px_Py;
  Double I_ERI_F3y_Pz_Px_Py = I_ERI_G3yz_S_Px_Py+ABZ*I_ERI_F3y_S_Px_Py;
  Double I_ERI_F2yz_Pz_Px_Py = I_ERI_G2y2z_S_Px_Py+ABZ*I_ERI_F2yz_S_Px_Py;
  Double I_ERI_Fy2z_Pz_Px_Py = I_ERI_Gy3z_S_Px_Py+ABZ*I_ERI_Fy2z_S_Px_Py;
  Double I_ERI_F3z_Pz_Px_Py = I_ERI_G4z_S_Px_Py+ABZ*I_ERI_F3z_S_Px_Py;
  Double I_ERI_F3x_Px_Py_Py = I_ERI_G4x_S_Py_Py+ABX*I_ERI_F3x_S_Py_Py;
  Double I_ERI_F2xy_Px_Py_Py = I_ERI_G3xy_S_Py_Py+ABX*I_ERI_F2xy_S_Py_Py;
  Double I_ERI_F2xz_Px_Py_Py = I_ERI_G3xz_S_Py_Py+ABX*I_ERI_F2xz_S_Py_Py;
  Double I_ERI_Fx2y_Px_Py_Py = I_ERI_G2x2y_S_Py_Py+ABX*I_ERI_Fx2y_S_Py_Py;
  Double I_ERI_Fxyz_Px_Py_Py = I_ERI_G2xyz_S_Py_Py+ABX*I_ERI_Fxyz_S_Py_Py;
  Double I_ERI_Fx2z_Px_Py_Py = I_ERI_G2x2z_S_Py_Py+ABX*I_ERI_Fx2z_S_Py_Py;
  Double I_ERI_F3y_Px_Py_Py = I_ERI_Gx3y_S_Py_Py+ABX*I_ERI_F3y_S_Py_Py;
  Double I_ERI_F2yz_Px_Py_Py = I_ERI_Gx2yz_S_Py_Py+ABX*I_ERI_F2yz_S_Py_Py;
  Double I_ERI_Fy2z_Px_Py_Py = I_ERI_Gxy2z_S_Py_Py+ABX*I_ERI_Fy2z_S_Py_Py;
  Double I_ERI_F3z_Px_Py_Py = I_ERI_Gx3z_S_Py_Py+ABX*I_ERI_F3z_S_Py_Py;
  Double I_ERI_F3x_Py_Py_Py = I_ERI_G3xy_S_Py_Py+ABY*I_ERI_F3x_S_Py_Py;
  Double I_ERI_F2xy_Py_Py_Py = I_ERI_G2x2y_S_Py_Py+ABY*I_ERI_F2xy_S_Py_Py;
  Double I_ERI_F2xz_Py_Py_Py = I_ERI_G2xyz_S_Py_Py+ABY*I_ERI_F2xz_S_Py_Py;
  Double I_ERI_Fx2y_Py_Py_Py = I_ERI_Gx3y_S_Py_Py+ABY*I_ERI_Fx2y_S_Py_Py;
  Double I_ERI_Fxyz_Py_Py_Py = I_ERI_Gx2yz_S_Py_Py+ABY*I_ERI_Fxyz_S_Py_Py;
  Double I_ERI_Fx2z_Py_Py_Py = I_ERI_Gxy2z_S_Py_Py+ABY*I_ERI_Fx2z_S_Py_Py;
  Double I_ERI_F3y_Py_Py_Py = I_ERI_G4y_S_Py_Py+ABY*I_ERI_F3y_S_Py_Py;
  Double I_ERI_F2yz_Py_Py_Py = I_ERI_G3yz_S_Py_Py+ABY*I_ERI_F2yz_S_Py_Py;
  Double I_ERI_Fy2z_Py_Py_Py = I_ERI_G2y2z_S_Py_Py+ABY*I_ERI_Fy2z_S_Py_Py;
  Double I_ERI_F3z_Py_Py_Py = I_ERI_Gy3z_S_Py_Py+ABY*I_ERI_F3z_S_Py_Py;
  Double I_ERI_F3x_Pz_Py_Py = I_ERI_G3xz_S_Py_Py+ABZ*I_ERI_F3x_S_Py_Py;
  Double I_ERI_F2xy_Pz_Py_Py = I_ERI_G2xyz_S_Py_Py+ABZ*I_ERI_F2xy_S_Py_Py;
  Double I_ERI_F2xz_Pz_Py_Py = I_ERI_G2x2z_S_Py_Py+ABZ*I_ERI_F2xz_S_Py_Py;
  Double I_ERI_Fx2y_Pz_Py_Py = I_ERI_Gx2yz_S_Py_Py+ABZ*I_ERI_Fx2y_S_Py_Py;
  Double I_ERI_Fxyz_Pz_Py_Py = I_ERI_Gxy2z_S_Py_Py+ABZ*I_ERI_Fxyz_S_Py_Py;
  Double I_ERI_Fx2z_Pz_Py_Py = I_ERI_Gx3z_S_Py_Py+ABZ*I_ERI_Fx2z_S_Py_Py;
  Double I_ERI_F3y_Pz_Py_Py = I_ERI_G3yz_S_Py_Py+ABZ*I_ERI_F3y_S_Py_Py;
  Double I_ERI_F2yz_Pz_Py_Py = I_ERI_G2y2z_S_Py_Py+ABZ*I_ERI_F2yz_S_Py_Py;
  Double I_ERI_Fy2z_Pz_Py_Py = I_ERI_Gy3z_S_Py_Py+ABZ*I_ERI_Fy2z_S_Py_Py;
  Double I_ERI_F3z_Pz_Py_Py = I_ERI_G4z_S_Py_Py+ABZ*I_ERI_F3z_S_Py_Py;
  Double I_ERI_F3x_Px_Pz_Py = I_ERI_G4x_S_Pz_Py+ABX*I_ERI_F3x_S_Pz_Py;
  Double I_ERI_F2xy_Px_Pz_Py = I_ERI_G3xy_S_Pz_Py+ABX*I_ERI_F2xy_S_Pz_Py;
  Double I_ERI_F2xz_Px_Pz_Py = I_ERI_G3xz_S_Pz_Py+ABX*I_ERI_F2xz_S_Pz_Py;
  Double I_ERI_Fx2y_Px_Pz_Py = I_ERI_G2x2y_S_Pz_Py+ABX*I_ERI_Fx2y_S_Pz_Py;
  Double I_ERI_Fxyz_Px_Pz_Py = I_ERI_G2xyz_S_Pz_Py+ABX*I_ERI_Fxyz_S_Pz_Py;
  Double I_ERI_Fx2z_Px_Pz_Py = I_ERI_G2x2z_S_Pz_Py+ABX*I_ERI_Fx2z_S_Pz_Py;
  Double I_ERI_F3y_Px_Pz_Py = I_ERI_Gx3y_S_Pz_Py+ABX*I_ERI_F3y_S_Pz_Py;
  Double I_ERI_F2yz_Px_Pz_Py = I_ERI_Gx2yz_S_Pz_Py+ABX*I_ERI_F2yz_S_Pz_Py;
  Double I_ERI_Fy2z_Px_Pz_Py = I_ERI_Gxy2z_S_Pz_Py+ABX*I_ERI_Fy2z_S_Pz_Py;
  Double I_ERI_F3z_Px_Pz_Py = I_ERI_Gx3z_S_Pz_Py+ABX*I_ERI_F3z_S_Pz_Py;
  Double I_ERI_F3x_Py_Pz_Py = I_ERI_G3xy_S_Pz_Py+ABY*I_ERI_F3x_S_Pz_Py;
  Double I_ERI_F2xy_Py_Pz_Py = I_ERI_G2x2y_S_Pz_Py+ABY*I_ERI_F2xy_S_Pz_Py;
  Double I_ERI_F2xz_Py_Pz_Py = I_ERI_G2xyz_S_Pz_Py+ABY*I_ERI_F2xz_S_Pz_Py;
  Double I_ERI_Fx2y_Py_Pz_Py = I_ERI_Gx3y_S_Pz_Py+ABY*I_ERI_Fx2y_S_Pz_Py;
  Double I_ERI_Fxyz_Py_Pz_Py = I_ERI_Gx2yz_S_Pz_Py+ABY*I_ERI_Fxyz_S_Pz_Py;
  Double I_ERI_Fx2z_Py_Pz_Py = I_ERI_Gxy2z_S_Pz_Py+ABY*I_ERI_Fx2z_S_Pz_Py;
  Double I_ERI_F3y_Py_Pz_Py = I_ERI_G4y_S_Pz_Py+ABY*I_ERI_F3y_S_Pz_Py;
  Double I_ERI_F2yz_Py_Pz_Py = I_ERI_G3yz_S_Pz_Py+ABY*I_ERI_F2yz_S_Pz_Py;
  Double I_ERI_Fy2z_Py_Pz_Py = I_ERI_G2y2z_S_Pz_Py+ABY*I_ERI_Fy2z_S_Pz_Py;
  Double I_ERI_F3z_Py_Pz_Py = I_ERI_Gy3z_S_Pz_Py+ABY*I_ERI_F3z_S_Pz_Py;
  Double I_ERI_F3x_Pz_Pz_Py = I_ERI_G3xz_S_Pz_Py+ABZ*I_ERI_F3x_S_Pz_Py;
  Double I_ERI_F2xy_Pz_Pz_Py = I_ERI_G2xyz_S_Pz_Py+ABZ*I_ERI_F2xy_S_Pz_Py;
  Double I_ERI_F2xz_Pz_Pz_Py = I_ERI_G2x2z_S_Pz_Py+ABZ*I_ERI_F2xz_S_Pz_Py;
  Double I_ERI_Fx2y_Pz_Pz_Py = I_ERI_Gx2yz_S_Pz_Py+ABZ*I_ERI_Fx2y_S_Pz_Py;
  Double I_ERI_Fxyz_Pz_Pz_Py = I_ERI_Gxy2z_S_Pz_Py+ABZ*I_ERI_Fxyz_S_Pz_Py;
  Double I_ERI_Fx2z_Pz_Pz_Py = I_ERI_Gx3z_S_Pz_Py+ABZ*I_ERI_Fx2z_S_Pz_Py;
  Double I_ERI_F3y_Pz_Pz_Py = I_ERI_G3yz_S_Pz_Py+ABZ*I_ERI_F3y_S_Pz_Py;
  Double I_ERI_F2yz_Pz_Pz_Py = I_ERI_G2y2z_S_Pz_Py+ABZ*I_ERI_F2yz_S_Pz_Py;
  Double I_ERI_Fy2z_Pz_Pz_Py = I_ERI_Gy3z_S_Pz_Py+ABZ*I_ERI_Fy2z_S_Pz_Py;
  Double I_ERI_F3z_Pz_Pz_Py = I_ERI_G4z_S_Pz_Py+ABZ*I_ERI_F3z_S_Pz_Py;
  Double I_ERI_F3x_Px_Px_Pz = I_ERI_G4x_S_Px_Pz+ABX*I_ERI_F3x_S_Px_Pz;
  Double I_ERI_F2xy_Px_Px_Pz = I_ERI_G3xy_S_Px_Pz+ABX*I_ERI_F2xy_S_Px_Pz;
  Double I_ERI_F2xz_Px_Px_Pz = I_ERI_G3xz_S_Px_Pz+ABX*I_ERI_F2xz_S_Px_Pz;
  Double I_ERI_Fx2y_Px_Px_Pz = I_ERI_G2x2y_S_Px_Pz+ABX*I_ERI_Fx2y_S_Px_Pz;
  Double I_ERI_Fxyz_Px_Px_Pz = I_ERI_G2xyz_S_Px_Pz+ABX*I_ERI_Fxyz_S_Px_Pz;
  Double I_ERI_Fx2z_Px_Px_Pz = I_ERI_G2x2z_S_Px_Pz+ABX*I_ERI_Fx2z_S_Px_Pz;
  Double I_ERI_F3y_Px_Px_Pz = I_ERI_Gx3y_S_Px_Pz+ABX*I_ERI_F3y_S_Px_Pz;
  Double I_ERI_F2yz_Px_Px_Pz = I_ERI_Gx2yz_S_Px_Pz+ABX*I_ERI_F2yz_S_Px_Pz;
  Double I_ERI_Fy2z_Px_Px_Pz = I_ERI_Gxy2z_S_Px_Pz+ABX*I_ERI_Fy2z_S_Px_Pz;
  Double I_ERI_F3z_Px_Px_Pz = I_ERI_Gx3z_S_Px_Pz+ABX*I_ERI_F3z_S_Px_Pz;
  Double I_ERI_F3x_Py_Px_Pz = I_ERI_G3xy_S_Px_Pz+ABY*I_ERI_F3x_S_Px_Pz;
  Double I_ERI_F2xy_Py_Px_Pz = I_ERI_G2x2y_S_Px_Pz+ABY*I_ERI_F2xy_S_Px_Pz;
  Double I_ERI_F2xz_Py_Px_Pz = I_ERI_G2xyz_S_Px_Pz+ABY*I_ERI_F2xz_S_Px_Pz;
  Double I_ERI_Fx2y_Py_Px_Pz = I_ERI_Gx3y_S_Px_Pz+ABY*I_ERI_Fx2y_S_Px_Pz;
  Double I_ERI_Fxyz_Py_Px_Pz = I_ERI_Gx2yz_S_Px_Pz+ABY*I_ERI_Fxyz_S_Px_Pz;
  Double I_ERI_Fx2z_Py_Px_Pz = I_ERI_Gxy2z_S_Px_Pz+ABY*I_ERI_Fx2z_S_Px_Pz;
  Double I_ERI_F3y_Py_Px_Pz = I_ERI_G4y_S_Px_Pz+ABY*I_ERI_F3y_S_Px_Pz;
  Double I_ERI_F2yz_Py_Px_Pz = I_ERI_G3yz_S_Px_Pz+ABY*I_ERI_F2yz_S_Px_Pz;
  Double I_ERI_Fy2z_Py_Px_Pz = I_ERI_G2y2z_S_Px_Pz+ABY*I_ERI_Fy2z_S_Px_Pz;
  Double I_ERI_F3z_Py_Px_Pz = I_ERI_Gy3z_S_Px_Pz+ABY*I_ERI_F3z_S_Px_Pz;
  Double I_ERI_F3x_Pz_Px_Pz = I_ERI_G3xz_S_Px_Pz+ABZ*I_ERI_F3x_S_Px_Pz;
  Double I_ERI_F2xy_Pz_Px_Pz = I_ERI_G2xyz_S_Px_Pz+ABZ*I_ERI_F2xy_S_Px_Pz;
  Double I_ERI_F2xz_Pz_Px_Pz = I_ERI_G2x2z_S_Px_Pz+ABZ*I_ERI_F2xz_S_Px_Pz;
  Double I_ERI_Fx2y_Pz_Px_Pz = I_ERI_Gx2yz_S_Px_Pz+ABZ*I_ERI_Fx2y_S_Px_Pz;
  Double I_ERI_Fxyz_Pz_Px_Pz = I_ERI_Gxy2z_S_Px_Pz+ABZ*I_ERI_Fxyz_S_Px_Pz;
  Double I_ERI_Fx2z_Pz_Px_Pz = I_ERI_Gx3z_S_Px_Pz+ABZ*I_ERI_Fx2z_S_Px_Pz;
  Double I_ERI_F3y_Pz_Px_Pz = I_ERI_G3yz_S_Px_Pz+ABZ*I_ERI_F3y_S_Px_Pz;
  Double I_ERI_F2yz_Pz_Px_Pz = I_ERI_G2y2z_S_Px_Pz+ABZ*I_ERI_F2yz_S_Px_Pz;
  Double I_ERI_Fy2z_Pz_Px_Pz = I_ERI_Gy3z_S_Px_Pz+ABZ*I_ERI_Fy2z_S_Px_Pz;
  Double I_ERI_F3z_Pz_Px_Pz = I_ERI_G4z_S_Px_Pz+ABZ*I_ERI_F3z_S_Px_Pz;
  Double I_ERI_F3x_Px_Py_Pz = I_ERI_G4x_S_Py_Pz+ABX*I_ERI_F3x_S_Py_Pz;
  Double I_ERI_F2xy_Px_Py_Pz = I_ERI_G3xy_S_Py_Pz+ABX*I_ERI_F2xy_S_Py_Pz;
  Double I_ERI_F2xz_Px_Py_Pz = I_ERI_G3xz_S_Py_Pz+ABX*I_ERI_F2xz_S_Py_Pz;
  Double I_ERI_Fx2y_Px_Py_Pz = I_ERI_G2x2y_S_Py_Pz+ABX*I_ERI_Fx2y_S_Py_Pz;
  Double I_ERI_Fxyz_Px_Py_Pz = I_ERI_G2xyz_S_Py_Pz+ABX*I_ERI_Fxyz_S_Py_Pz;
  Double I_ERI_Fx2z_Px_Py_Pz = I_ERI_G2x2z_S_Py_Pz+ABX*I_ERI_Fx2z_S_Py_Pz;
  Double I_ERI_F3y_Px_Py_Pz = I_ERI_Gx3y_S_Py_Pz+ABX*I_ERI_F3y_S_Py_Pz;
  Double I_ERI_F2yz_Px_Py_Pz = I_ERI_Gx2yz_S_Py_Pz+ABX*I_ERI_F2yz_S_Py_Pz;
  Double I_ERI_Fy2z_Px_Py_Pz = I_ERI_Gxy2z_S_Py_Pz+ABX*I_ERI_Fy2z_S_Py_Pz;
  Double I_ERI_F3z_Px_Py_Pz = I_ERI_Gx3z_S_Py_Pz+ABX*I_ERI_F3z_S_Py_Pz;
  Double I_ERI_F3x_Py_Py_Pz = I_ERI_G3xy_S_Py_Pz+ABY*I_ERI_F3x_S_Py_Pz;
  Double I_ERI_F2xy_Py_Py_Pz = I_ERI_G2x2y_S_Py_Pz+ABY*I_ERI_F2xy_S_Py_Pz;
  Double I_ERI_F2xz_Py_Py_Pz = I_ERI_G2xyz_S_Py_Pz+ABY*I_ERI_F2xz_S_Py_Pz;
  Double I_ERI_Fx2y_Py_Py_Pz = I_ERI_Gx3y_S_Py_Pz+ABY*I_ERI_Fx2y_S_Py_Pz;
  Double I_ERI_Fxyz_Py_Py_Pz = I_ERI_Gx2yz_S_Py_Pz+ABY*I_ERI_Fxyz_S_Py_Pz;
  Double I_ERI_Fx2z_Py_Py_Pz = I_ERI_Gxy2z_S_Py_Pz+ABY*I_ERI_Fx2z_S_Py_Pz;
  Double I_ERI_F3y_Py_Py_Pz = I_ERI_G4y_S_Py_Pz+ABY*I_ERI_F3y_S_Py_Pz;
  Double I_ERI_F2yz_Py_Py_Pz = I_ERI_G3yz_S_Py_Pz+ABY*I_ERI_F2yz_S_Py_Pz;
  Double I_ERI_Fy2z_Py_Py_Pz = I_ERI_G2y2z_S_Py_Pz+ABY*I_ERI_Fy2z_S_Py_Pz;
  Double I_ERI_F3z_Py_Py_Pz = I_ERI_Gy3z_S_Py_Pz+ABY*I_ERI_F3z_S_Py_Pz;
  Double I_ERI_F3x_Pz_Py_Pz = I_ERI_G3xz_S_Py_Pz+ABZ*I_ERI_F3x_S_Py_Pz;
  Double I_ERI_F2xy_Pz_Py_Pz = I_ERI_G2xyz_S_Py_Pz+ABZ*I_ERI_F2xy_S_Py_Pz;
  Double I_ERI_F2xz_Pz_Py_Pz = I_ERI_G2x2z_S_Py_Pz+ABZ*I_ERI_F2xz_S_Py_Pz;
  Double I_ERI_Fx2y_Pz_Py_Pz = I_ERI_Gx2yz_S_Py_Pz+ABZ*I_ERI_Fx2y_S_Py_Pz;
  Double I_ERI_Fxyz_Pz_Py_Pz = I_ERI_Gxy2z_S_Py_Pz+ABZ*I_ERI_Fxyz_S_Py_Pz;
  Double I_ERI_Fx2z_Pz_Py_Pz = I_ERI_Gx3z_S_Py_Pz+ABZ*I_ERI_Fx2z_S_Py_Pz;
  Double I_ERI_F3y_Pz_Py_Pz = I_ERI_G3yz_S_Py_Pz+ABZ*I_ERI_F3y_S_Py_Pz;
  Double I_ERI_F2yz_Pz_Py_Pz = I_ERI_G2y2z_S_Py_Pz+ABZ*I_ERI_F2yz_S_Py_Pz;
  Double I_ERI_Fy2z_Pz_Py_Pz = I_ERI_Gy3z_S_Py_Pz+ABZ*I_ERI_Fy2z_S_Py_Pz;
  Double I_ERI_F3z_Pz_Py_Pz = I_ERI_G4z_S_Py_Pz+ABZ*I_ERI_F3z_S_Py_Pz;
  Double I_ERI_F3x_Px_Pz_Pz = I_ERI_G4x_S_Pz_Pz+ABX*I_ERI_F3x_S_Pz_Pz;
  Double I_ERI_F2xy_Px_Pz_Pz = I_ERI_G3xy_S_Pz_Pz+ABX*I_ERI_F2xy_S_Pz_Pz;
  Double I_ERI_F2xz_Px_Pz_Pz = I_ERI_G3xz_S_Pz_Pz+ABX*I_ERI_F2xz_S_Pz_Pz;
  Double I_ERI_Fx2y_Px_Pz_Pz = I_ERI_G2x2y_S_Pz_Pz+ABX*I_ERI_Fx2y_S_Pz_Pz;
  Double I_ERI_Fxyz_Px_Pz_Pz = I_ERI_G2xyz_S_Pz_Pz+ABX*I_ERI_Fxyz_S_Pz_Pz;
  Double I_ERI_Fx2z_Px_Pz_Pz = I_ERI_G2x2z_S_Pz_Pz+ABX*I_ERI_Fx2z_S_Pz_Pz;
  Double I_ERI_F3y_Px_Pz_Pz = I_ERI_Gx3y_S_Pz_Pz+ABX*I_ERI_F3y_S_Pz_Pz;
  Double I_ERI_F2yz_Px_Pz_Pz = I_ERI_Gx2yz_S_Pz_Pz+ABX*I_ERI_F2yz_S_Pz_Pz;
  Double I_ERI_Fy2z_Px_Pz_Pz = I_ERI_Gxy2z_S_Pz_Pz+ABX*I_ERI_Fy2z_S_Pz_Pz;
  Double I_ERI_F3z_Px_Pz_Pz = I_ERI_Gx3z_S_Pz_Pz+ABX*I_ERI_F3z_S_Pz_Pz;
  Double I_ERI_F3x_Py_Pz_Pz = I_ERI_G3xy_S_Pz_Pz+ABY*I_ERI_F3x_S_Pz_Pz;
  Double I_ERI_F2xy_Py_Pz_Pz = I_ERI_G2x2y_S_Pz_Pz+ABY*I_ERI_F2xy_S_Pz_Pz;
  Double I_ERI_F2xz_Py_Pz_Pz = I_ERI_G2xyz_S_Pz_Pz+ABY*I_ERI_F2xz_S_Pz_Pz;
  Double I_ERI_Fx2y_Py_Pz_Pz = I_ERI_Gx3y_S_Pz_Pz+ABY*I_ERI_Fx2y_S_Pz_Pz;
  Double I_ERI_Fxyz_Py_Pz_Pz = I_ERI_Gx2yz_S_Pz_Pz+ABY*I_ERI_Fxyz_S_Pz_Pz;
  Double I_ERI_Fx2z_Py_Pz_Pz = I_ERI_Gxy2z_S_Pz_Pz+ABY*I_ERI_Fx2z_S_Pz_Pz;
  Double I_ERI_F3y_Py_Pz_Pz = I_ERI_G4y_S_Pz_Pz+ABY*I_ERI_F3y_S_Pz_Pz;
  Double I_ERI_F2yz_Py_Pz_Pz = I_ERI_G3yz_S_Pz_Pz+ABY*I_ERI_F2yz_S_Pz_Pz;
  Double I_ERI_Fy2z_Py_Pz_Pz = I_ERI_G2y2z_S_Pz_Pz+ABY*I_ERI_Fy2z_S_Pz_Pz;
  Double I_ERI_F3z_Py_Pz_Pz = I_ERI_Gy3z_S_Pz_Pz+ABY*I_ERI_F3z_S_Pz_Pz;
  Double I_ERI_F3x_Pz_Pz_Pz = I_ERI_G3xz_S_Pz_Pz+ABZ*I_ERI_F3x_S_Pz_Pz;
  Double I_ERI_F2xy_Pz_Pz_Pz = I_ERI_G2xyz_S_Pz_Pz+ABZ*I_ERI_F2xy_S_Pz_Pz;
  Double I_ERI_F2xz_Pz_Pz_Pz = I_ERI_G2x2z_S_Pz_Pz+ABZ*I_ERI_F2xz_S_Pz_Pz;
  Double I_ERI_Fx2y_Pz_Pz_Pz = I_ERI_Gx2yz_S_Pz_Pz+ABZ*I_ERI_Fx2y_S_Pz_Pz;
  Double I_ERI_Fxyz_Pz_Pz_Pz = I_ERI_Gxy2z_S_Pz_Pz+ABZ*I_ERI_Fxyz_S_Pz_Pz;
  Double I_ERI_Fx2z_Pz_Pz_Pz = I_ERI_Gx3z_S_Pz_Pz+ABZ*I_ERI_Fx2z_S_Pz_Pz;
  Double I_ERI_F3y_Pz_Pz_Pz = I_ERI_G3yz_S_Pz_Pz+ABZ*I_ERI_F3y_S_Pz_Pz;
  Double I_ERI_F2yz_Pz_Pz_Pz = I_ERI_G2y2z_S_Pz_Pz+ABZ*I_ERI_F2yz_S_Pz_Pz;
  Double I_ERI_Fy2z_Pz_Pz_Pz = I_ERI_Gy3z_S_Pz_Pz+ABZ*I_ERI_Fy2z_S_Pz_Pz;
  Double I_ERI_F3z_Pz_Pz_Pz = I_ERI_G4z_S_Pz_Pz+ABZ*I_ERI_F3z_S_Pz_Pz;

  /************************************************************
   * shell quartet name: SQ_ERI_G_P_P_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 54 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_H_S_P_P
   * RHS shell quartet name: SQ_ERI_G_S_P_P
   ************************************************************/
  Double I_ERI_G4x_Px_Px_Px = I_ERI_H5x_S_Px_Px+ABX*I_ERI_G4x_S_Px_Px;
  Double I_ERI_G3xy_Px_Px_Px = I_ERI_H4xy_S_Px_Px+ABX*I_ERI_G3xy_S_Px_Px;
  Double I_ERI_G3xz_Px_Px_Px = I_ERI_H4xz_S_Px_Px+ABX*I_ERI_G3xz_S_Px_Px;
  Double I_ERI_G2x2y_Px_Px_Px = I_ERI_H3x2y_S_Px_Px+ABX*I_ERI_G2x2y_S_Px_Px;
  Double I_ERI_G2xyz_Px_Px_Px = I_ERI_H3xyz_S_Px_Px+ABX*I_ERI_G2xyz_S_Px_Px;
  Double I_ERI_G2x2z_Px_Px_Px = I_ERI_H3x2z_S_Px_Px+ABX*I_ERI_G2x2z_S_Px_Px;
  Double I_ERI_Gx3y_Px_Px_Px = I_ERI_H2x3y_S_Px_Px+ABX*I_ERI_Gx3y_S_Px_Px;
  Double I_ERI_Gx2yz_Px_Px_Px = I_ERI_H2x2yz_S_Px_Px+ABX*I_ERI_Gx2yz_S_Px_Px;
  Double I_ERI_Gxy2z_Px_Px_Px = I_ERI_H2xy2z_S_Px_Px+ABX*I_ERI_Gxy2z_S_Px_Px;
  Double I_ERI_Gx3z_Px_Px_Px = I_ERI_H2x3z_S_Px_Px+ABX*I_ERI_Gx3z_S_Px_Px;
  Double I_ERI_G4y_Px_Px_Px = I_ERI_Hx4y_S_Px_Px+ABX*I_ERI_G4y_S_Px_Px;
  Double I_ERI_G3yz_Px_Px_Px = I_ERI_Hx3yz_S_Px_Px+ABX*I_ERI_G3yz_S_Px_Px;
  Double I_ERI_G2y2z_Px_Px_Px = I_ERI_Hx2y2z_S_Px_Px+ABX*I_ERI_G2y2z_S_Px_Px;
  Double I_ERI_Gy3z_Px_Px_Px = I_ERI_Hxy3z_S_Px_Px+ABX*I_ERI_Gy3z_S_Px_Px;
  Double I_ERI_G4z_Px_Px_Px = I_ERI_Hx4z_S_Px_Px+ABX*I_ERI_G4z_S_Px_Px;
  Double I_ERI_G3xy_Py_Px_Px = I_ERI_H3x2y_S_Px_Px+ABY*I_ERI_G3xy_S_Px_Px;
  Double I_ERI_G3xz_Py_Px_Px = I_ERI_H3xyz_S_Px_Px+ABY*I_ERI_G3xz_S_Px_Px;
  Double I_ERI_G2x2y_Py_Px_Px = I_ERI_H2x3y_S_Px_Px+ABY*I_ERI_G2x2y_S_Px_Px;
  Double I_ERI_G2xyz_Py_Px_Px = I_ERI_H2x2yz_S_Px_Px+ABY*I_ERI_G2xyz_S_Px_Px;
  Double I_ERI_G2x2z_Py_Px_Px = I_ERI_H2xy2z_S_Px_Px+ABY*I_ERI_G2x2z_S_Px_Px;
  Double I_ERI_Gx3y_Py_Px_Px = I_ERI_Hx4y_S_Px_Px+ABY*I_ERI_Gx3y_S_Px_Px;
  Double I_ERI_Gx2yz_Py_Px_Px = I_ERI_Hx3yz_S_Px_Px+ABY*I_ERI_Gx2yz_S_Px_Px;
  Double I_ERI_Gxy2z_Py_Px_Px = I_ERI_Hx2y2z_S_Px_Px+ABY*I_ERI_Gxy2z_S_Px_Px;
  Double I_ERI_Gx3z_Py_Px_Px = I_ERI_Hxy3z_S_Px_Px+ABY*I_ERI_Gx3z_S_Px_Px;
  Double I_ERI_G4y_Py_Px_Px = I_ERI_H5y_S_Px_Px+ABY*I_ERI_G4y_S_Px_Px;
  Double I_ERI_G3yz_Py_Px_Px = I_ERI_H4yz_S_Px_Px+ABY*I_ERI_G3yz_S_Px_Px;
  Double I_ERI_G2y2z_Py_Px_Px = I_ERI_H3y2z_S_Px_Px+ABY*I_ERI_G2y2z_S_Px_Px;
  Double I_ERI_Gy3z_Py_Px_Px = I_ERI_H2y3z_S_Px_Px+ABY*I_ERI_Gy3z_S_Px_Px;
  Double I_ERI_G4z_Py_Px_Px = I_ERI_Hy4z_S_Px_Px+ABY*I_ERI_G4z_S_Px_Px;
  Double I_ERI_G3xz_Pz_Px_Px = I_ERI_H3x2z_S_Px_Px+ABZ*I_ERI_G3xz_S_Px_Px;
  Double I_ERI_G2xyz_Pz_Px_Px = I_ERI_H2xy2z_S_Px_Px+ABZ*I_ERI_G2xyz_S_Px_Px;
  Double I_ERI_G2x2z_Pz_Px_Px = I_ERI_H2x3z_S_Px_Px+ABZ*I_ERI_G2x2z_S_Px_Px;
  Double I_ERI_Gx2yz_Pz_Px_Px = I_ERI_Hx2y2z_S_Px_Px+ABZ*I_ERI_Gx2yz_S_Px_Px;
  Double I_ERI_Gxy2z_Pz_Px_Px = I_ERI_Hxy3z_S_Px_Px+ABZ*I_ERI_Gxy2z_S_Px_Px;
  Double I_ERI_Gx3z_Pz_Px_Px = I_ERI_Hx4z_S_Px_Px+ABZ*I_ERI_Gx3z_S_Px_Px;
  Double I_ERI_G3yz_Pz_Px_Px = I_ERI_H3y2z_S_Px_Px+ABZ*I_ERI_G3yz_S_Px_Px;
  Double I_ERI_G2y2z_Pz_Px_Px = I_ERI_H2y3z_S_Px_Px+ABZ*I_ERI_G2y2z_S_Px_Px;
  Double I_ERI_Gy3z_Pz_Px_Px = I_ERI_Hy4z_S_Px_Px+ABZ*I_ERI_Gy3z_S_Px_Px;
  Double I_ERI_G4z_Pz_Px_Px = I_ERI_H5z_S_Px_Px+ABZ*I_ERI_G4z_S_Px_Px;
  Double I_ERI_G4x_Px_Py_Px = I_ERI_H5x_S_Py_Px+ABX*I_ERI_G4x_S_Py_Px;
  Double I_ERI_G3xy_Px_Py_Px = I_ERI_H4xy_S_Py_Px+ABX*I_ERI_G3xy_S_Py_Px;
  Double I_ERI_G3xz_Px_Py_Px = I_ERI_H4xz_S_Py_Px+ABX*I_ERI_G3xz_S_Py_Px;
  Double I_ERI_G2x2y_Px_Py_Px = I_ERI_H3x2y_S_Py_Px+ABX*I_ERI_G2x2y_S_Py_Px;
  Double I_ERI_G2xyz_Px_Py_Px = I_ERI_H3xyz_S_Py_Px+ABX*I_ERI_G2xyz_S_Py_Px;
  Double I_ERI_G2x2z_Px_Py_Px = I_ERI_H3x2z_S_Py_Px+ABX*I_ERI_G2x2z_S_Py_Px;
  Double I_ERI_Gx3y_Px_Py_Px = I_ERI_H2x3y_S_Py_Px+ABX*I_ERI_Gx3y_S_Py_Px;
  Double I_ERI_Gx2yz_Px_Py_Px = I_ERI_H2x2yz_S_Py_Px+ABX*I_ERI_Gx2yz_S_Py_Px;
  Double I_ERI_Gxy2z_Px_Py_Px = I_ERI_H2xy2z_S_Py_Px+ABX*I_ERI_Gxy2z_S_Py_Px;
  Double I_ERI_Gx3z_Px_Py_Px = I_ERI_H2x3z_S_Py_Px+ABX*I_ERI_Gx3z_S_Py_Px;
  Double I_ERI_G4y_Px_Py_Px = I_ERI_Hx4y_S_Py_Px+ABX*I_ERI_G4y_S_Py_Px;
  Double I_ERI_G3yz_Px_Py_Px = I_ERI_Hx3yz_S_Py_Px+ABX*I_ERI_G3yz_S_Py_Px;
  Double I_ERI_G2y2z_Px_Py_Px = I_ERI_Hx2y2z_S_Py_Px+ABX*I_ERI_G2y2z_S_Py_Px;
  Double I_ERI_Gy3z_Px_Py_Px = I_ERI_Hxy3z_S_Py_Px+ABX*I_ERI_Gy3z_S_Py_Px;
  Double I_ERI_G4z_Px_Py_Px = I_ERI_Hx4z_S_Py_Px+ABX*I_ERI_G4z_S_Py_Px;
  Double I_ERI_G3xy_Py_Py_Px = I_ERI_H3x2y_S_Py_Px+ABY*I_ERI_G3xy_S_Py_Px;
  Double I_ERI_G3xz_Py_Py_Px = I_ERI_H3xyz_S_Py_Px+ABY*I_ERI_G3xz_S_Py_Px;
  Double I_ERI_G2x2y_Py_Py_Px = I_ERI_H2x3y_S_Py_Px+ABY*I_ERI_G2x2y_S_Py_Px;
  Double I_ERI_G2xyz_Py_Py_Px = I_ERI_H2x2yz_S_Py_Px+ABY*I_ERI_G2xyz_S_Py_Px;
  Double I_ERI_G2x2z_Py_Py_Px = I_ERI_H2xy2z_S_Py_Px+ABY*I_ERI_G2x2z_S_Py_Px;
  Double I_ERI_Gx3y_Py_Py_Px = I_ERI_Hx4y_S_Py_Px+ABY*I_ERI_Gx3y_S_Py_Px;
  Double I_ERI_Gx2yz_Py_Py_Px = I_ERI_Hx3yz_S_Py_Px+ABY*I_ERI_Gx2yz_S_Py_Px;
  Double I_ERI_Gxy2z_Py_Py_Px = I_ERI_Hx2y2z_S_Py_Px+ABY*I_ERI_Gxy2z_S_Py_Px;
  Double I_ERI_Gx3z_Py_Py_Px = I_ERI_Hxy3z_S_Py_Px+ABY*I_ERI_Gx3z_S_Py_Px;
  Double I_ERI_G4y_Py_Py_Px = I_ERI_H5y_S_Py_Px+ABY*I_ERI_G4y_S_Py_Px;
  Double I_ERI_G3yz_Py_Py_Px = I_ERI_H4yz_S_Py_Px+ABY*I_ERI_G3yz_S_Py_Px;
  Double I_ERI_G2y2z_Py_Py_Px = I_ERI_H3y2z_S_Py_Px+ABY*I_ERI_G2y2z_S_Py_Px;
  Double I_ERI_Gy3z_Py_Py_Px = I_ERI_H2y3z_S_Py_Px+ABY*I_ERI_Gy3z_S_Py_Px;
  Double I_ERI_G4z_Py_Py_Px = I_ERI_Hy4z_S_Py_Px+ABY*I_ERI_G4z_S_Py_Px;
  Double I_ERI_G3xz_Pz_Py_Px = I_ERI_H3x2z_S_Py_Px+ABZ*I_ERI_G3xz_S_Py_Px;
  Double I_ERI_G2xyz_Pz_Py_Px = I_ERI_H2xy2z_S_Py_Px+ABZ*I_ERI_G2xyz_S_Py_Px;
  Double I_ERI_G2x2z_Pz_Py_Px = I_ERI_H2x3z_S_Py_Px+ABZ*I_ERI_G2x2z_S_Py_Px;
  Double I_ERI_Gx2yz_Pz_Py_Px = I_ERI_Hx2y2z_S_Py_Px+ABZ*I_ERI_Gx2yz_S_Py_Px;
  Double I_ERI_Gxy2z_Pz_Py_Px = I_ERI_Hxy3z_S_Py_Px+ABZ*I_ERI_Gxy2z_S_Py_Px;
  Double I_ERI_Gx3z_Pz_Py_Px = I_ERI_Hx4z_S_Py_Px+ABZ*I_ERI_Gx3z_S_Py_Px;
  Double I_ERI_G3yz_Pz_Py_Px = I_ERI_H3y2z_S_Py_Px+ABZ*I_ERI_G3yz_S_Py_Px;
  Double I_ERI_G2y2z_Pz_Py_Px = I_ERI_H2y3z_S_Py_Px+ABZ*I_ERI_G2y2z_S_Py_Px;
  Double I_ERI_Gy3z_Pz_Py_Px = I_ERI_Hy4z_S_Py_Px+ABZ*I_ERI_Gy3z_S_Py_Px;
  Double I_ERI_G4z_Pz_Py_Px = I_ERI_H5z_S_Py_Px+ABZ*I_ERI_G4z_S_Py_Px;
  Double I_ERI_G4x_Px_Pz_Px = I_ERI_H5x_S_Pz_Px+ABX*I_ERI_G4x_S_Pz_Px;
  Double I_ERI_G3xy_Px_Pz_Px = I_ERI_H4xy_S_Pz_Px+ABX*I_ERI_G3xy_S_Pz_Px;
  Double I_ERI_G3xz_Px_Pz_Px = I_ERI_H4xz_S_Pz_Px+ABX*I_ERI_G3xz_S_Pz_Px;
  Double I_ERI_G2x2y_Px_Pz_Px = I_ERI_H3x2y_S_Pz_Px+ABX*I_ERI_G2x2y_S_Pz_Px;
  Double I_ERI_G2xyz_Px_Pz_Px = I_ERI_H3xyz_S_Pz_Px+ABX*I_ERI_G2xyz_S_Pz_Px;
  Double I_ERI_G2x2z_Px_Pz_Px = I_ERI_H3x2z_S_Pz_Px+ABX*I_ERI_G2x2z_S_Pz_Px;
  Double I_ERI_Gx3y_Px_Pz_Px = I_ERI_H2x3y_S_Pz_Px+ABX*I_ERI_Gx3y_S_Pz_Px;
  Double I_ERI_Gx2yz_Px_Pz_Px = I_ERI_H2x2yz_S_Pz_Px+ABX*I_ERI_Gx2yz_S_Pz_Px;
  Double I_ERI_Gxy2z_Px_Pz_Px = I_ERI_H2xy2z_S_Pz_Px+ABX*I_ERI_Gxy2z_S_Pz_Px;
  Double I_ERI_Gx3z_Px_Pz_Px = I_ERI_H2x3z_S_Pz_Px+ABX*I_ERI_Gx3z_S_Pz_Px;
  Double I_ERI_G4y_Px_Pz_Px = I_ERI_Hx4y_S_Pz_Px+ABX*I_ERI_G4y_S_Pz_Px;
  Double I_ERI_G3yz_Px_Pz_Px = I_ERI_Hx3yz_S_Pz_Px+ABX*I_ERI_G3yz_S_Pz_Px;
  Double I_ERI_G2y2z_Px_Pz_Px = I_ERI_Hx2y2z_S_Pz_Px+ABX*I_ERI_G2y2z_S_Pz_Px;
  Double I_ERI_Gy3z_Px_Pz_Px = I_ERI_Hxy3z_S_Pz_Px+ABX*I_ERI_Gy3z_S_Pz_Px;
  Double I_ERI_G4z_Px_Pz_Px = I_ERI_Hx4z_S_Pz_Px+ABX*I_ERI_G4z_S_Pz_Px;
  Double I_ERI_G3xy_Py_Pz_Px = I_ERI_H3x2y_S_Pz_Px+ABY*I_ERI_G3xy_S_Pz_Px;
  Double I_ERI_G3xz_Py_Pz_Px = I_ERI_H3xyz_S_Pz_Px+ABY*I_ERI_G3xz_S_Pz_Px;
  Double I_ERI_G2x2y_Py_Pz_Px = I_ERI_H2x3y_S_Pz_Px+ABY*I_ERI_G2x2y_S_Pz_Px;
  Double I_ERI_G2xyz_Py_Pz_Px = I_ERI_H2x2yz_S_Pz_Px+ABY*I_ERI_G2xyz_S_Pz_Px;
  Double I_ERI_G2x2z_Py_Pz_Px = I_ERI_H2xy2z_S_Pz_Px+ABY*I_ERI_G2x2z_S_Pz_Px;
  Double I_ERI_Gx3y_Py_Pz_Px = I_ERI_Hx4y_S_Pz_Px+ABY*I_ERI_Gx3y_S_Pz_Px;
  Double I_ERI_Gx2yz_Py_Pz_Px = I_ERI_Hx3yz_S_Pz_Px+ABY*I_ERI_Gx2yz_S_Pz_Px;
  Double I_ERI_Gxy2z_Py_Pz_Px = I_ERI_Hx2y2z_S_Pz_Px+ABY*I_ERI_Gxy2z_S_Pz_Px;
  Double I_ERI_Gx3z_Py_Pz_Px = I_ERI_Hxy3z_S_Pz_Px+ABY*I_ERI_Gx3z_S_Pz_Px;
  Double I_ERI_G4y_Py_Pz_Px = I_ERI_H5y_S_Pz_Px+ABY*I_ERI_G4y_S_Pz_Px;
  Double I_ERI_G3yz_Py_Pz_Px = I_ERI_H4yz_S_Pz_Px+ABY*I_ERI_G3yz_S_Pz_Px;
  Double I_ERI_G2y2z_Py_Pz_Px = I_ERI_H3y2z_S_Pz_Px+ABY*I_ERI_G2y2z_S_Pz_Px;
  Double I_ERI_Gy3z_Py_Pz_Px = I_ERI_H2y3z_S_Pz_Px+ABY*I_ERI_Gy3z_S_Pz_Px;
  Double I_ERI_G4z_Py_Pz_Px = I_ERI_Hy4z_S_Pz_Px+ABY*I_ERI_G4z_S_Pz_Px;
  Double I_ERI_G3xz_Pz_Pz_Px = I_ERI_H3x2z_S_Pz_Px+ABZ*I_ERI_G3xz_S_Pz_Px;
  Double I_ERI_G2xyz_Pz_Pz_Px = I_ERI_H2xy2z_S_Pz_Px+ABZ*I_ERI_G2xyz_S_Pz_Px;
  Double I_ERI_G2x2z_Pz_Pz_Px = I_ERI_H2x3z_S_Pz_Px+ABZ*I_ERI_G2x2z_S_Pz_Px;
  Double I_ERI_Gx2yz_Pz_Pz_Px = I_ERI_Hx2y2z_S_Pz_Px+ABZ*I_ERI_Gx2yz_S_Pz_Px;
  Double I_ERI_Gxy2z_Pz_Pz_Px = I_ERI_Hxy3z_S_Pz_Px+ABZ*I_ERI_Gxy2z_S_Pz_Px;
  Double I_ERI_Gx3z_Pz_Pz_Px = I_ERI_Hx4z_S_Pz_Px+ABZ*I_ERI_Gx3z_S_Pz_Px;
  Double I_ERI_G3yz_Pz_Pz_Px = I_ERI_H3y2z_S_Pz_Px+ABZ*I_ERI_G3yz_S_Pz_Px;
  Double I_ERI_G2y2z_Pz_Pz_Px = I_ERI_H2y3z_S_Pz_Px+ABZ*I_ERI_G2y2z_S_Pz_Px;
  Double I_ERI_Gy3z_Pz_Pz_Px = I_ERI_Hy4z_S_Pz_Px+ABZ*I_ERI_Gy3z_S_Pz_Px;
  Double I_ERI_G4z_Pz_Pz_Px = I_ERI_H5z_S_Pz_Px+ABZ*I_ERI_G4z_S_Pz_Px;
  Double I_ERI_G4x_Px_Px_Py = I_ERI_H5x_S_Px_Py+ABX*I_ERI_G4x_S_Px_Py;
  Double I_ERI_G3xy_Px_Px_Py = I_ERI_H4xy_S_Px_Py+ABX*I_ERI_G3xy_S_Px_Py;
  Double I_ERI_G3xz_Px_Px_Py = I_ERI_H4xz_S_Px_Py+ABX*I_ERI_G3xz_S_Px_Py;
  Double I_ERI_G2x2y_Px_Px_Py = I_ERI_H3x2y_S_Px_Py+ABX*I_ERI_G2x2y_S_Px_Py;
  Double I_ERI_G2xyz_Px_Px_Py = I_ERI_H3xyz_S_Px_Py+ABX*I_ERI_G2xyz_S_Px_Py;
  Double I_ERI_G2x2z_Px_Px_Py = I_ERI_H3x2z_S_Px_Py+ABX*I_ERI_G2x2z_S_Px_Py;
  Double I_ERI_Gx3y_Px_Px_Py = I_ERI_H2x3y_S_Px_Py+ABX*I_ERI_Gx3y_S_Px_Py;
  Double I_ERI_Gx2yz_Px_Px_Py = I_ERI_H2x2yz_S_Px_Py+ABX*I_ERI_Gx2yz_S_Px_Py;
  Double I_ERI_Gxy2z_Px_Px_Py = I_ERI_H2xy2z_S_Px_Py+ABX*I_ERI_Gxy2z_S_Px_Py;
  Double I_ERI_Gx3z_Px_Px_Py = I_ERI_H2x3z_S_Px_Py+ABX*I_ERI_Gx3z_S_Px_Py;
  Double I_ERI_G4y_Px_Px_Py = I_ERI_Hx4y_S_Px_Py+ABX*I_ERI_G4y_S_Px_Py;
  Double I_ERI_G3yz_Px_Px_Py = I_ERI_Hx3yz_S_Px_Py+ABX*I_ERI_G3yz_S_Px_Py;
  Double I_ERI_G2y2z_Px_Px_Py = I_ERI_Hx2y2z_S_Px_Py+ABX*I_ERI_G2y2z_S_Px_Py;
  Double I_ERI_Gy3z_Px_Px_Py = I_ERI_Hxy3z_S_Px_Py+ABX*I_ERI_Gy3z_S_Px_Py;
  Double I_ERI_G4z_Px_Px_Py = I_ERI_Hx4z_S_Px_Py+ABX*I_ERI_G4z_S_Px_Py;
  Double I_ERI_G3xy_Py_Px_Py = I_ERI_H3x2y_S_Px_Py+ABY*I_ERI_G3xy_S_Px_Py;
  Double I_ERI_G3xz_Py_Px_Py = I_ERI_H3xyz_S_Px_Py+ABY*I_ERI_G3xz_S_Px_Py;
  Double I_ERI_G2x2y_Py_Px_Py = I_ERI_H2x3y_S_Px_Py+ABY*I_ERI_G2x2y_S_Px_Py;
  Double I_ERI_G2xyz_Py_Px_Py = I_ERI_H2x2yz_S_Px_Py+ABY*I_ERI_G2xyz_S_Px_Py;
  Double I_ERI_G2x2z_Py_Px_Py = I_ERI_H2xy2z_S_Px_Py+ABY*I_ERI_G2x2z_S_Px_Py;
  Double I_ERI_Gx3y_Py_Px_Py = I_ERI_Hx4y_S_Px_Py+ABY*I_ERI_Gx3y_S_Px_Py;
  Double I_ERI_Gx2yz_Py_Px_Py = I_ERI_Hx3yz_S_Px_Py+ABY*I_ERI_Gx2yz_S_Px_Py;
  Double I_ERI_Gxy2z_Py_Px_Py = I_ERI_Hx2y2z_S_Px_Py+ABY*I_ERI_Gxy2z_S_Px_Py;
  Double I_ERI_Gx3z_Py_Px_Py = I_ERI_Hxy3z_S_Px_Py+ABY*I_ERI_Gx3z_S_Px_Py;
  Double I_ERI_G4y_Py_Px_Py = I_ERI_H5y_S_Px_Py+ABY*I_ERI_G4y_S_Px_Py;
  Double I_ERI_G3yz_Py_Px_Py = I_ERI_H4yz_S_Px_Py+ABY*I_ERI_G3yz_S_Px_Py;
  Double I_ERI_G2y2z_Py_Px_Py = I_ERI_H3y2z_S_Px_Py+ABY*I_ERI_G2y2z_S_Px_Py;
  Double I_ERI_Gy3z_Py_Px_Py = I_ERI_H2y3z_S_Px_Py+ABY*I_ERI_Gy3z_S_Px_Py;
  Double I_ERI_G4z_Py_Px_Py = I_ERI_Hy4z_S_Px_Py+ABY*I_ERI_G4z_S_Px_Py;
  Double I_ERI_G3xz_Pz_Px_Py = I_ERI_H3x2z_S_Px_Py+ABZ*I_ERI_G3xz_S_Px_Py;
  Double I_ERI_G2xyz_Pz_Px_Py = I_ERI_H2xy2z_S_Px_Py+ABZ*I_ERI_G2xyz_S_Px_Py;
  Double I_ERI_G2x2z_Pz_Px_Py = I_ERI_H2x3z_S_Px_Py+ABZ*I_ERI_G2x2z_S_Px_Py;
  Double I_ERI_Gx2yz_Pz_Px_Py = I_ERI_Hx2y2z_S_Px_Py+ABZ*I_ERI_Gx2yz_S_Px_Py;
  Double I_ERI_Gxy2z_Pz_Px_Py = I_ERI_Hxy3z_S_Px_Py+ABZ*I_ERI_Gxy2z_S_Px_Py;
  Double I_ERI_Gx3z_Pz_Px_Py = I_ERI_Hx4z_S_Px_Py+ABZ*I_ERI_Gx3z_S_Px_Py;
  Double I_ERI_G3yz_Pz_Px_Py = I_ERI_H3y2z_S_Px_Py+ABZ*I_ERI_G3yz_S_Px_Py;
  Double I_ERI_G2y2z_Pz_Px_Py = I_ERI_H2y3z_S_Px_Py+ABZ*I_ERI_G2y2z_S_Px_Py;
  Double I_ERI_Gy3z_Pz_Px_Py = I_ERI_Hy4z_S_Px_Py+ABZ*I_ERI_Gy3z_S_Px_Py;
  Double I_ERI_G4z_Pz_Px_Py = I_ERI_H5z_S_Px_Py+ABZ*I_ERI_G4z_S_Px_Py;
  Double I_ERI_G4x_Px_Py_Py = I_ERI_H5x_S_Py_Py+ABX*I_ERI_G4x_S_Py_Py;
  Double I_ERI_G3xy_Px_Py_Py = I_ERI_H4xy_S_Py_Py+ABX*I_ERI_G3xy_S_Py_Py;
  Double I_ERI_G3xz_Px_Py_Py = I_ERI_H4xz_S_Py_Py+ABX*I_ERI_G3xz_S_Py_Py;
  Double I_ERI_G2x2y_Px_Py_Py = I_ERI_H3x2y_S_Py_Py+ABX*I_ERI_G2x2y_S_Py_Py;
  Double I_ERI_G2xyz_Px_Py_Py = I_ERI_H3xyz_S_Py_Py+ABX*I_ERI_G2xyz_S_Py_Py;
  Double I_ERI_G2x2z_Px_Py_Py = I_ERI_H3x2z_S_Py_Py+ABX*I_ERI_G2x2z_S_Py_Py;
  Double I_ERI_Gx3y_Px_Py_Py = I_ERI_H2x3y_S_Py_Py+ABX*I_ERI_Gx3y_S_Py_Py;
  Double I_ERI_Gx2yz_Px_Py_Py = I_ERI_H2x2yz_S_Py_Py+ABX*I_ERI_Gx2yz_S_Py_Py;
  Double I_ERI_Gxy2z_Px_Py_Py = I_ERI_H2xy2z_S_Py_Py+ABX*I_ERI_Gxy2z_S_Py_Py;
  Double I_ERI_Gx3z_Px_Py_Py = I_ERI_H2x3z_S_Py_Py+ABX*I_ERI_Gx3z_S_Py_Py;
  Double I_ERI_G4y_Px_Py_Py = I_ERI_Hx4y_S_Py_Py+ABX*I_ERI_G4y_S_Py_Py;
  Double I_ERI_G3yz_Px_Py_Py = I_ERI_Hx3yz_S_Py_Py+ABX*I_ERI_G3yz_S_Py_Py;
  Double I_ERI_G2y2z_Px_Py_Py = I_ERI_Hx2y2z_S_Py_Py+ABX*I_ERI_G2y2z_S_Py_Py;
  Double I_ERI_Gy3z_Px_Py_Py = I_ERI_Hxy3z_S_Py_Py+ABX*I_ERI_Gy3z_S_Py_Py;
  Double I_ERI_G4z_Px_Py_Py = I_ERI_Hx4z_S_Py_Py+ABX*I_ERI_G4z_S_Py_Py;
  Double I_ERI_G3xy_Py_Py_Py = I_ERI_H3x2y_S_Py_Py+ABY*I_ERI_G3xy_S_Py_Py;
  Double I_ERI_G3xz_Py_Py_Py = I_ERI_H3xyz_S_Py_Py+ABY*I_ERI_G3xz_S_Py_Py;
  Double I_ERI_G2x2y_Py_Py_Py = I_ERI_H2x3y_S_Py_Py+ABY*I_ERI_G2x2y_S_Py_Py;
  Double I_ERI_G2xyz_Py_Py_Py = I_ERI_H2x2yz_S_Py_Py+ABY*I_ERI_G2xyz_S_Py_Py;
  Double I_ERI_G2x2z_Py_Py_Py = I_ERI_H2xy2z_S_Py_Py+ABY*I_ERI_G2x2z_S_Py_Py;
  Double I_ERI_Gx3y_Py_Py_Py = I_ERI_Hx4y_S_Py_Py+ABY*I_ERI_Gx3y_S_Py_Py;
  Double I_ERI_Gx2yz_Py_Py_Py = I_ERI_Hx3yz_S_Py_Py+ABY*I_ERI_Gx2yz_S_Py_Py;
  Double I_ERI_Gxy2z_Py_Py_Py = I_ERI_Hx2y2z_S_Py_Py+ABY*I_ERI_Gxy2z_S_Py_Py;
  Double I_ERI_Gx3z_Py_Py_Py = I_ERI_Hxy3z_S_Py_Py+ABY*I_ERI_Gx3z_S_Py_Py;
  Double I_ERI_G4y_Py_Py_Py = I_ERI_H5y_S_Py_Py+ABY*I_ERI_G4y_S_Py_Py;
  Double I_ERI_G3yz_Py_Py_Py = I_ERI_H4yz_S_Py_Py+ABY*I_ERI_G3yz_S_Py_Py;
  Double I_ERI_G2y2z_Py_Py_Py = I_ERI_H3y2z_S_Py_Py+ABY*I_ERI_G2y2z_S_Py_Py;
  Double I_ERI_Gy3z_Py_Py_Py = I_ERI_H2y3z_S_Py_Py+ABY*I_ERI_Gy3z_S_Py_Py;
  Double I_ERI_G4z_Py_Py_Py = I_ERI_Hy4z_S_Py_Py+ABY*I_ERI_G4z_S_Py_Py;
  Double I_ERI_G3xz_Pz_Py_Py = I_ERI_H3x2z_S_Py_Py+ABZ*I_ERI_G3xz_S_Py_Py;
  Double I_ERI_G2xyz_Pz_Py_Py = I_ERI_H2xy2z_S_Py_Py+ABZ*I_ERI_G2xyz_S_Py_Py;
  Double I_ERI_G2x2z_Pz_Py_Py = I_ERI_H2x3z_S_Py_Py+ABZ*I_ERI_G2x2z_S_Py_Py;
  Double I_ERI_Gx2yz_Pz_Py_Py = I_ERI_Hx2y2z_S_Py_Py+ABZ*I_ERI_Gx2yz_S_Py_Py;
  Double I_ERI_Gxy2z_Pz_Py_Py = I_ERI_Hxy3z_S_Py_Py+ABZ*I_ERI_Gxy2z_S_Py_Py;
  Double I_ERI_Gx3z_Pz_Py_Py = I_ERI_Hx4z_S_Py_Py+ABZ*I_ERI_Gx3z_S_Py_Py;
  Double I_ERI_G3yz_Pz_Py_Py = I_ERI_H3y2z_S_Py_Py+ABZ*I_ERI_G3yz_S_Py_Py;
  Double I_ERI_G2y2z_Pz_Py_Py = I_ERI_H2y3z_S_Py_Py+ABZ*I_ERI_G2y2z_S_Py_Py;
  Double I_ERI_Gy3z_Pz_Py_Py = I_ERI_Hy4z_S_Py_Py+ABZ*I_ERI_Gy3z_S_Py_Py;
  Double I_ERI_G4z_Pz_Py_Py = I_ERI_H5z_S_Py_Py+ABZ*I_ERI_G4z_S_Py_Py;
  Double I_ERI_G4x_Px_Pz_Py = I_ERI_H5x_S_Pz_Py+ABX*I_ERI_G4x_S_Pz_Py;
  Double I_ERI_G3xy_Px_Pz_Py = I_ERI_H4xy_S_Pz_Py+ABX*I_ERI_G3xy_S_Pz_Py;
  Double I_ERI_G3xz_Px_Pz_Py = I_ERI_H4xz_S_Pz_Py+ABX*I_ERI_G3xz_S_Pz_Py;
  Double I_ERI_G2x2y_Px_Pz_Py = I_ERI_H3x2y_S_Pz_Py+ABX*I_ERI_G2x2y_S_Pz_Py;
  Double I_ERI_G2xyz_Px_Pz_Py = I_ERI_H3xyz_S_Pz_Py+ABX*I_ERI_G2xyz_S_Pz_Py;
  Double I_ERI_G2x2z_Px_Pz_Py = I_ERI_H3x2z_S_Pz_Py+ABX*I_ERI_G2x2z_S_Pz_Py;
  Double I_ERI_Gx3y_Px_Pz_Py = I_ERI_H2x3y_S_Pz_Py+ABX*I_ERI_Gx3y_S_Pz_Py;
  Double I_ERI_Gx2yz_Px_Pz_Py = I_ERI_H2x2yz_S_Pz_Py+ABX*I_ERI_Gx2yz_S_Pz_Py;
  Double I_ERI_Gxy2z_Px_Pz_Py = I_ERI_H2xy2z_S_Pz_Py+ABX*I_ERI_Gxy2z_S_Pz_Py;
  Double I_ERI_Gx3z_Px_Pz_Py = I_ERI_H2x3z_S_Pz_Py+ABX*I_ERI_Gx3z_S_Pz_Py;
  Double I_ERI_G4y_Px_Pz_Py = I_ERI_Hx4y_S_Pz_Py+ABX*I_ERI_G4y_S_Pz_Py;
  Double I_ERI_G3yz_Px_Pz_Py = I_ERI_Hx3yz_S_Pz_Py+ABX*I_ERI_G3yz_S_Pz_Py;
  Double I_ERI_G2y2z_Px_Pz_Py = I_ERI_Hx2y2z_S_Pz_Py+ABX*I_ERI_G2y2z_S_Pz_Py;
  Double I_ERI_Gy3z_Px_Pz_Py = I_ERI_Hxy3z_S_Pz_Py+ABX*I_ERI_Gy3z_S_Pz_Py;
  Double I_ERI_G4z_Px_Pz_Py = I_ERI_Hx4z_S_Pz_Py+ABX*I_ERI_G4z_S_Pz_Py;
  Double I_ERI_G3xy_Py_Pz_Py = I_ERI_H3x2y_S_Pz_Py+ABY*I_ERI_G3xy_S_Pz_Py;
  Double I_ERI_G3xz_Py_Pz_Py = I_ERI_H3xyz_S_Pz_Py+ABY*I_ERI_G3xz_S_Pz_Py;
  Double I_ERI_G2x2y_Py_Pz_Py = I_ERI_H2x3y_S_Pz_Py+ABY*I_ERI_G2x2y_S_Pz_Py;
  Double I_ERI_G2xyz_Py_Pz_Py = I_ERI_H2x2yz_S_Pz_Py+ABY*I_ERI_G2xyz_S_Pz_Py;
  Double I_ERI_G2x2z_Py_Pz_Py = I_ERI_H2xy2z_S_Pz_Py+ABY*I_ERI_G2x2z_S_Pz_Py;
  Double I_ERI_Gx3y_Py_Pz_Py = I_ERI_Hx4y_S_Pz_Py+ABY*I_ERI_Gx3y_S_Pz_Py;
  Double I_ERI_Gx2yz_Py_Pz_Py = I_ERI_Hx3yz_S_Pz_Py+ABY*I_ERI_Gx2yz_S_Pz_Py;
  Double I_ERI_Gxy2z_Py_Pz_Py = I_ERI_Hx2y2z_S_Pz_Py+ABY*I_ERI_Gxy2z_S_Pz_Py;
  Double I_ERI_Gx3z_Py_Pz_Py = I_ERI_Hxy3z_S_Pz_Py+ABY*I_ERI_Gx3z_S_Pz_Py;
  Double I_ERI_G4y_Py_Pz_Py = I_ERI_H5y_S_Pz_Py+ABY*I_ERI_G4y_S_Pz_Py;
  Double I_ERI_G3yz_Py_Pz_Py = I_ERI_H4yz_S_Pz_Py+ABY*I_ERI_G3yz_S_Pz_Py;
  Double I_ERI_G2y2z_Py_Pz_Py = I_ERI_H3y2z_S_Pz_Py+ABY*I_ERI_G2y2z_S_Pz_Py;
  Double I_ERI_Gy3z_Py_Pz_Py = I_ERI_H2y3z_S_Pz_Py+ABY*I_ERI_Gy3z_S_Pz_Py;
  Double I_ERI_G4z_Py_Pz_Py = I_ERI_Hy4z_S_Pz_Py+ABY*I_ERI_G4z_S_Pz_Py;
  Double I_ERI_G3xz_Pz_Pz_Py = I_ERI_H3x2z_S_Pz_Py+ABZ*I_ERI_G3xz_S_Pz_Py;
  Double I_ERI_G2xyz_Pz_Pz_Py = I_ERI_H2xy2z_S_Pz_Py+ABZ*I_ERI_G2xyz_S_Pz_Py;
  Double I_ERI_G2x2z_Pz_Pz_Py = I_ERI_H2x3z_S_Pz_Py+ABZ*I_ERI_G2x2z_S_Pz_Py;
  Double I_ERI_Gx2yz_Pz_Pz_Py = I_ERI_Hx2y2z_S_Pz_Py+ABZ*I_ERI_Gx2yz_S_Pz_Py;
  Double I_ERI_Gxy2z_Pz_Pz_Py = I_ERI_Hxy3z_S_Pz_Py+ABZ*I_ERI_Gxy2z_S_Pz_Py;
  Double I_ERI_Gx3z_Pz_Pz_Py = I_ERI_Hx4z_S_Pz_Py+ABZ*I_ERI_Gx3z_S_Pz_Py;
  Double I_ERI_G3yz_Pz_Pz_Py = I_ERI_H3y2z_S_Pz_Py+ABZ*I_ERI_G3yz_S_Pz_Py;
  Double I_ERI_G2y2z_Pz_Pz_Py = I_ERI_H2y3z_S_Pz_Py+ABZ*I_ERI_G2y2z_S_Pz_Py;
  Double I_ERI_Gy3z_Pz_Pz_Py = I_ERI_Hy4z_S_Pz_Py+ABZ*I_ERI_Gy3z_S_Pz_Py;
  Double I_ERI_G4z_Pz_Pz_Py = I_ERI_H5z_S_Pz_Py+ABZ*I_ERI_G4z_S_Pz_Py;
  Double I_ERI_G4x_Px_Px_Pz = I_ERI_H5x_S_Px_Pz+ABX*I_ERI_G4x_S_Px_Pz;
  Double I_ERI_G3xy_Px_Px_Pz = I_ERI_H4xy_S_Px_Pz+ABX*I_ERI_G3xy_S_Px_Pz;
  Double I_ERI_G3xz_Px_Px_Pz = I_ERI_H4xz_S_Px_Pz+ABX*I_ERI_G3xz_S_Px_Pz;
  Double I_ERI_G2x2y_Px_Px_Pz = I_ERI_H3x2y_S_Px_Pz+ABX*I_ERI_G2x2y_S_Px_Pz;
  Double I_ERI_G2xyz_Px_Px_Pz = I_ERI_H3xyz_S_Px_Pz+ABX*I_ERI_G2xyz_S_Px_Pz;
  Double I_ERI_G2x2z_Px_Px_Pz = I_ERI_H3x2z_S_Px_Pz+ABX*I_ERI_G2x2z_S_Px_Pz;
  Double I_ERI_Gx3y_Px_Px_Pz = I_ERI_H2x3y_S_Px_Pz+ABX*I_ERI_Gx3y_S_Px_Pz;
  Double I_ERI_Gx2yz_Px_Px_Pz = I_ERI_H2x2yz_S_Px_Pz+ABX*I_ERI_Gx2yz_S_Px_Pz;
  Double I_ERI_Gxy2z_Px_Px_Pz = I_ERI_H2xy2z_S_Px_Pz+ABX*I_ERI_Gxy2z_S_Px_Pz;
  Double I_ERI_Gx3z_Px_Px_Pz = I_ERI_H2x3z_S_Px_Pz+ABX*I_ERI_Gx3z_S_Px_Pz;
  Double I_ERI_G4y_Px_Px_Pz = I_ERI_Hx4y_S_Px_Pz+ABX*I_ERI_G4y_S_Px_Pz;
  Double I_ERI_G3yz_Px_Px_Pz = I_ERI_Hx3yz_S_Px_Pz+ABX*I_ERI_G3yz_S_Px_Pz;
  Double I_ERI_G2y2z_Px_Px_Pz = I_ERI_Hx2y2z_S_Px_Pz+ABX*I_ERI_G2y2z_S_Px_Pz;
  Double I_ERI_Gy3z_Px_Px_Pz = I_ERI_Hxy3z_S_Px_Pz+ABX*I_ERI_Gy3z_S_Px_Pz;
  Double I_ERI_G4z_Px_Px_Pz = I_ERI_Hx4z_S_Px_Pz+ABX*I_ERI_G4z_S_Px_Pz;
  Double I_ERI_G3xy_Py_Px_Pz = I_ERI_H3x2y_S_Px_Pz+ABY*I_ERI_G3xy_S_Px_Pz;
  Double I_ERI_G3xz_Py_Px_Pz = I_ERI_H3xyz_S_Px_Pz+ABY*I_ERI_G3xz_S_Px_Pz;
  Double I_ERI_G2x2y_Py_Px_Pz = I_ERI_H2x3y_S_Px_Pz+ABY*I_ERI_G2x2y_S_Px_Pz;
  Double I_ERI_G2xyz_Py_Px_Pz = I_ERI_H2x2yz_S_Px_Pz+ABY*I_ERI_G2xyz_S_Px_Pz;
  Double I_ERI_G2x2z_Py_Px_Pz = I_ERI_H2xy2z_S_Px_Pz+ABY*I_ERI_G2x2z_S_Px_Pz;
  Double I_ERI_Gx3y_Py_Px_Pz = I_ERI_Hx4y_S_Px_Pz+ABY*I_ERI_Gx3y_S_Px_Pz;
  Double I_ERI_Gx2yz_Py_Px_Pz = I_ERI_Hx3yz_S_Px_Pz+ABY*I_ERI_Gx2yz_S_Px_Pz;
  Double I_ERI_Gxy2z_Py_Px_Pz = I_ERI_Hx2y2z_S_Px_Pz+ABY*I_ERI_Gxy2z_S_Px_Pz;
  Double I_ERI_Gx3z_Py_Px_Pz = I_ERI_Hxy3z_S_Px_Pz+ABY*I_ERI_Gx3z_S_Px_Pz;
  Double I_ERI_G4y_Py_Px_Pz = I_ERI_H5y_S_Px_Pz+ABY*I_ERI_G4y_S_Px_Pz;
  Double I_ERI_G3yz_Py_Px_Pz = I_ERI_H4yz_S_Px_Pz+ABY*I_ERI_G3yz_S_Px_Pz;
  Double I_ERI_G2y2z_Py_Px_Pz = I_ERI_H3y2z_S_Px_Pz+ABY*I_ERI_G2y2z_S_Px_Pz;
  Double I_ERI_Gy3z_Py_Px_Pz = I_ERI_H2y3z_S_Px_Pz+ABY*I_ERI_Gy3z_S_Px_Pz;
  Double I_ERI_G4z_Py_Px_Pz = I_ERI_Hy4z_S_Px_Pz+ABY*I_ERI_G4z_S_Px_Pz;
  Double I_ERI_G3xz_Pz_Px_Pz = I_ERI_H3x2z_S_Px_Pz+ABZ*I_ERI_G3xz_S_Px_Pz;
  Double I_ERI_G2xyz_Pz_Px_Pz = I_ERI_H2xy2z_S_Px_Pz+ABZ*I_ERI_G2xyz_S_Px_Pz;
  Double I_ERI_G2x2z_Pz_Px_Pz = I_ERI_H2x3z_S_Px_Pz+ABZ*I_ERI_G2x2z_S_Px_Pz;
  Double I_ERI_Gx2yz_Pz_Px_Pz = I_ERI_Hx2y2z_S_Px_Pz+ABZ*I_ERI_Gx2yz_S_Px_Pz;
  Double I_ERI_Gxy2z_Pz_Px_Pz = I_ERI_Hxy3z_S_Px_Pz+ABZ*I_ERI_Gxy2z_S_Px_Pz;
  Double I_ERI_Gx3z_Pz_Px_Pz = I_ERI_Hx4z_S_Px_Pz+ABZ*I_ERI_Gx3z_S_Px_Pz;
  Double I_ERI_G3yz_Pz_Px_Pz = I_ERI_H3y2z_S_Px_Pz+ABZ*I_ERI_G3yz_S_Px_Pz;
  Double I_ERI_G2y2z_Pz_Px_Pz = I_ERI_H2y3z_S_Px_Pz+ABZ*I_ERI_G2y2z_S_Px_Pz;
  Double I_ERI_Gy3z_Pz_Px_Pz = I_ERI_Hy4z_S_Px_Pz+ABZ*I_ERI_Gy3z_S_Px_Pz;
  Double I_ERI_G4z_Pz_Px_Pz = I_ERI_H5z_S_Px_Pz+ABZ*I_ERI_G4z_S_Px_Pz;
  Double I_ERI_G4x_Px_Py_Pz = I_ERI_H5x_S_Py_Pz+ABX*I_ERI_G4x_S_Py_Pz;
  Double I_ERI_G3xy_Px_Py_Pz = I_ERI_H4xy_S_Py_Pz+ABX*I_ERI_G3xy_S_Py_Pz;
  Double I_ERI_G3xz_Px_Py_Pz = I_ERI_H4xz_S_Py_Pz+ABX*I_ERI_G3xz_S_Py_Pz;
  Double I_ERI_G2x2y_Px_Py_Pz = I_ERI_H3x2y_S_Py_Pz+ABX*I_ERI_G2x2y_S_Py_Pz;
  Double I_ERI_G2xyz_Px_Py_Pz = I_ERI_H3xyz_S_Py_Pz+ABX*I_ERI_G2xyz_S_Py_Pz;
  Double I_ERI_G2x2z_Px_Py_Pz = I_ERI_H3x2z_S_Py_Pz+ABX*I_ERI_G2x2z_S_Py_Pz;
  Double I_ERI_Gx3y_Px_Py_Pz = I_ERI_H2x3y_S_Py_Pz+ABX*I_ERI_Gx3y_S_Py_Pz;
  Double I_ERI_Gx2yz_Px_Py_Pz = I_ERI_H2x2yz_S_Py_Pz+ABX*I_ERI_Gx2yz_S_Py_Pz;
  Double I_ERI_Gxy2z_Px_Py_Pz = I_ERI_H2xy2z_S_Py_Pz+ABX*I_ERI_Gxy2z_S_Py_Pz;
  Double I_ERI_Gx3z_Px_Py_Pz = I_ERI_H2x3z_S_Py_Pz+ABX*I_ERI_Gx3z_S_Py_Pz;
  Double I_ERI_G4y_Px_Py_Pz = I_ERI_Hx4y_S_Py_Pz+ABX*I_ERI_G4y_S_Py_Pz;
  Double I_ERI_G3yz_Px_Py_Pz = I_ERI_Hx3yz_S_Py_Pz+ABX*I_ERI_G3yz_S_Py_Pz;
  Double I_ERI_G2y2z_Px_Py_Pz = I_ERI_Hx2y2z_S_Py_Pz+ABX*I_ERI_G2y2z_S_Py_Pz;
  Double I_ERI_Gy3z_Px_Py_Pz = I_ERI_Hxy3z_S_Py_Pz+ABX*I_ERI_Gy3z_S_Py_Pz;
  Double I_ERI_G4z_Px_Py_Pz = I_ERI_Hx4z_S_Py_Pz+ABX*I_ERI_G4z_S_Py_Pz;
  Double I_ERI_G3xy_Py_Py_Pz = I_ERI_H3x2y_S_Py_Pz+ABY*I_ERI_G3xy_S_Py_Pz;
  Double I_ERI_G3xz_Py_Py_Pz = I_ERI_H3xyz_S_Py_Pz+ABY*I_ERI_G3xz_S_Py_Pz;
  Double I_ERI_G2x2y_Py_Py_Pz = I_ERI_H2x3y_S_Py_Pz+ABY*I_ERI_G2x2y_S_Py_Pz;
  Double I_ERI_G2xyz_Py_Py_Pz = I_ERI_H2x2yz_S_Py_Pz+ABY*I_ERI_G2xyz_S_Py_Pz;
  Double I_ERI_G2x2z_Py_Py_Pz = I_ERI_H2xy2z_S_Py_Pz+ABY*I_ERI_G2x2z_S_Py_Pz;
  Double I_ERI_Gx3y_Py_Py_Pz = I_ERI_Hx4y_S_Py_Pz+ABY*I_ERI_Gx3y_S_Py_Pz;
  Double I_ERI_Gx2yz_Py_Py_Pz = I_ERI_Hx3yz_S_Py_Pz+ABY*I_ERI_Gx2yz_S_Py_Pz;
  Double I_ERI_Gxy2z_Py_Py_Pz = I_ERI_Hx2y2z_S_Py_Pz+ABY*I_ERI_Gxy2z_S_Py_Pz;
  Double I_ERI_Gx3z_Py_Py_Pz = I_ERI_Hxy3z_S_Py_Pz+ABY*I_ERI_Gx3z_S_Py_Pz;
  Double I_ERI_G4y_Py_Py_Pz = I_ERI_H5y_S_Py_Pz+ABY*I_ERI_G4y_S_Py_Pz;
  Double I_ERI_G3yz_Py_Py_Pz = I_ERI_H4yz_S_Py_Pz+ABY*I_ERI_G3yz_S_Py_Pz;
  Double I_ERI_G2y2z_Py_Py_Pz = I_ERI_H3y2z_S_Py_Pz+ABY*I_ERI_G2y2z_S_Py_Pz;
  Double I_ERI_Gy3z_Py_Py_Pz = I_ERI_H2y3z_S_Py_Pz+ABY*I_ERI_Gy3z_S_Py_Pz;
  Double I_ERI_G4z_Py_Py_Pz = I_ERI_Hy4z_S_Py_Pz+ABY*I_ERI_G4z_S_Py_Pz;
  Double I_ERI_G3xz_Pz_Py_Pz = I_ERI_H3x2z_S_Py_Pz+ABZ*I_ERI_G3xz_S_Py_Pz;
  Double I_ERI_G2xyz_Pz_Py_Pz = I_ERI_H2xy2z_S_Py_Pz+ABZ*I_ERI_G2xyz_S_Py_Pz;
  Double I_ERI_G2x2z_Pz_Py_Pz = I_ERI_H2x3z_S_Py_Pz+ABZ*I_ERI_G2x2z_S_Py_Pz;
  Double I_ERI_Gx2yz_Pz_Py_Pz = I_ERI_Hx2y2z_S_Py_Pz+ABZ*I_ERI_Gx2yz_S_Py_Pz;
  Double I_ERI_Gxy2z_Pz_Py_Pz = I_ERI_Hxy3z_S_Py_Pz+ABZ*I_ERI_Gxy2z_S_Py_Pz;
  Double I_ERI_Gx3z_Pz_Py_Pz = I_ERI_Hx4z_S_Py_Pz+ABZ*I_ERI_Gx3z_S_Py_Pz;
  Double I_ERI_G3yz_Pz_Py_Pz = I_ERI_H3y2z_S_Py_Pz+ABZ*I_ERI_G3yz_S_Py_Pz;
  Double I_ERI_G2y2z_Pz_Py_Pz = I_ERI_H2y3z_S_Py_Pz+ABZ*I_ERI_G2y2z_S_Py_Pz;
  Double I_ERI_Gy3z_Pz_Py_Pz = I_ERI_Hy4z_S_Py_Pz+ABZ*I_ERI_Gy3z_S_Py_Pz;
  Double I_ERI_G4z_Pz_Py_Pz = I_ERI_H5z_S_Py_Pz+ABZ*I_ERI_G4z_S_Py_Pz;
  Double I_ERI_G4x_Px_Pz_Pz = I_ERI_H5x_S_Pz_Pz+ABX*I_ERI_G4x_S_Pz_Pz;
  Double I_ERI_G3xy_Px_Pz_Pz = I_ERI_H4xy_S_Pz_Pz+ABX*I_ERI_G3xy_S_Pz_Pz;
  Double I_ERI_G3xz_Px_Pz_Pz = I_ERI_H4xz_S_Pz_Pz+ABX*I_ERI_G3xz_S_Pz_Pz;
  Double I_ERI_G2x2y_Px_Pz_Pz = I_ERI_H3x2y_S_Pz_Pz+ABX*I_ERI_G2x2y_S_Pz_Pz;
  Double I_ERI_G2xyz_Px_Pz_Pz = I_ERI_H3xyz_S_Pz_Pz+ABX*I_ERI_G2xyz_S_Pz_Pz;
  Double I_ERI_G2x2z_Px_Pz_Pz = I_ERI_H3x2z_S_Pz_Pz+ABX*I_ERI_G2x2z_S_Pz_Pz;
  Double I_ERI_Gx3y_Px_Pz_Pz = I_ERI_H2x3y_S_Pz_Pz+ABX*I_ERI_Gx3y_S_Pz_Pz;
  Double I_ERI_Gx2yz_Px_Pz_Pz = I_ERI_H2x2yz_S_Pz_Pz+ABX*I_ERI_Gx2yz_S_Pz_Pz;
  Double I_ERI_Gxy2z_Px_Pz_Pz = I_ERI_H2xy2z_S_Pz_Pz+ABX*I_ERI_Gxy2z_S_Pz_Pz;
  Double I_ERI_Gx3z_Px_Pz_Pz = I_ERI_H2x3z_S_Pz_Pz+ABX*I_ERI_Gx3z_S_Pz_Pz;
  Double I_ERI_G4y_Px_Pz_Pz = I_ERI_Hx4y_S_Pz_Pz+ABX*I_ERI_G4y_S_Pz_Pz;
  Double I_ERI_G3yz_Px_Pz_Pz = I_ERI_Hx3yz_S_Pz_Pz+ABX*I_ERI_G3yz_S_Pz_Pz;
  Double I_ERI_G2y2z_Px_Pz_Pz = I_ERI_Hx2y2z_S_Pz_Pz+ABX*I_ERI_G2y2z_S_Pz_Pz;
  Double I_ERI_Gy3z_Px_Pz_Pz = I_ERI_Hxy3z_S_Pz_Pz+ABX*I_ERI_Gy3z_S_Pz_Pz;
  Double I_ERI_G4z_Px_Pz_Pz = I_ERI_Hx4z_S_Pz_Pz+ABX*I_ERI_G4z_S_Pz_Pz;
  Double I_ERI_G3xy_Py_Pz_Pz = I_ERI_H3x2y_S_Pz_Pz+ABY*I_ERI_G3xy_S_Pz_Pz;
  Double I_ERI_G3xz_Py_Pz_Pz = I_ERI_H3xyz_S_Pz_Pz+ABY*I_ERI_G3xz_S_Pz_Pz;
  Double I_ERI_G2x2y_Py_Pz_Pz = I_ERI_H2x3y_S_Pz_Pz+ABY*I_ERI_G2x2y_S_Pz_Pz;
  Double I_ERI_G2xyz_Py_Pz_Pz = I_ERI_H2x2yz_S_Pz_Pz+ABY*I_ERI_G2xyz_S_Pz_Pz;
  Double I_ERI_G2x2z_Py_Pz_Pz = I_ERI_H2xy2z_S_Pz_Pz+ABY*I_ERI_G2x2z_S_Pz_Pz;
  Double I_ERI_Gx3y_Py_Pz_Pz = I_ERI_Hx4y_S_Pz_Pz+ABY*I_ERI_Gx3y_S_Pz_Pz;
  Double I_ERI_Gx2yz_Py_Pz_Pz = I_ERI_Hx3yz_S_Pz_Pz+ABY*I_ERI_Gx2yz_S_Pz_Pz;
  Double I_ERI_Gxy2z_Py_Pz_Pz = I_ERI_Hx2y2z_S_Pz_Pz+ABY*I_ERI_Gxy2z_S_Pz_Pz;
  Double I_ERI_Gx3z_Py_Pz_Pz = I_ERI_Hxy3z_S_Pz_Pz+ABY*I_ERI_Gx3z_S_Pz_Pz;
  Double I_ERI_G4y_Py_Pz_Pz = I_ERI_H5y_S_Pz_Pz+ABY*I_ERI_G4y_S_Pz_Pz;
  Double I_ERI_G3yz_Py_Pz_Pz = I_ERI_H4yz_S_Pz_Pz+ABY*I_ERI_G3yz_S_Pz_Pz;
  Double I_ERI_G2y2z_Py_Pz_Pz = I_ERI_H3y2z_S_Pz_Pz+ABY*I_ERI_G2y2z_S_Pz_Pz;
  Double I_ERI_Gy3z_Py_Pz_Pz = I_ERI_H2y3z_S_Pz_Pz+ABY*I_ERI_Gy3z_S_Pz_Pz;
  Double I_ERI_G4z_Py_Pz_Pz = I_ERI_Hy4z_S_Pz_Pz+ABY*I_ERI_G4z_S_Pz_Pz;
  Double I_ERI_G3xz_Pz_Pz_Pz = I_ERI_H3x2z_S_Pz_Pz+ABZ*I_ERI_G3xz_S_Pz_Pz;
  Double I_ERI_G2xyz_Pz_Pz_Pz = I_ERI_H2xy2z_S_Pz_Pz+ABZ*I_ERI_G2xyz_S_Pz_Pz;
  Double I_ERI_G2x2z_Pz_Pz_Pz = I_ERI_H2x3z_S_Pz_Pz+ABZ*I_ERI_G2x2z_S_Pz_Pz;
  Double I_ERI_Gx2yz_Pz_Pz_Pz = I_ERI_Hx2y2z_S_Pz_Pz+ABZ*I_ERI_Gx2yz_S_Pz_Pz;
  Double I_ERI_Gxy2z_Pz_Pz_Pz = I_ERI_Hxy3z_S_Pz_Pz+ABZ*I_ERI_Gxy2z_S_Pz_Pz;
  Double I_ERI_Gx3z_Pz_Pz_Pz = I_ERI_Hx4z_S_Pz_Pz+ABZ*I_ERI_Gx3z_S_Pz_Pz;
  Double I_ERI_G3yz_Pz_Pz_Pz = I_ERI_H3y2z_S_Pz_Pz+ABZ*I_ERI_G3yz_S_Pz_Pz;
  Double I_ERI_G2y2z_Pz_Pz_Pz = I_ERI_H2y3z_S_Pz_Pz+ABZ*I_ERI_G2y2z_S_Pz_Pz;
  Double I_ERI_Gy3z_Pz_Pz_Pz = I_ERI_Hy4z_S_Pz_Pz+ABZ*I_ERI_Gy3z_S_Pz_Pz;
  Double I_ERI_G4z_Pz_Pz_Pz = I_ERI_H5z_S_Pz_Pz+ABZ*I_ERI_G4z_S_Pz_Pz;

  /************************************************************
   * shell quartet name: SQ_ERI_F_D_P_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_ERI_G_P_P_P
   * RHS shell quartet name: SQ_ERI_F_P_P_P
   ************************************************************/
  abcd[0] = I_ERI_G4x_Px_Px_Px+ABX*I_ERI_F3x_Px_Px_Px;
  abcd[1] = I_ERI_G3xy_Px_Px_Px+ABX*I_ERI_F2xy_Px_Px_Px;
  abcd[2] = I_ERI_G3xz_Px_Px_Px+ABX*I_ERI_F2xz_Px_Px_Px;
  abcd[3] = I_ERI_G2x2y_Px_Px_Px+ABX*I_ERI_Fx2y_Px_Px_Px;
  abcd[4] = I_ERI_G2xyz_Px_Px_Px+ABX*I_ERI_Fxyz_Px_Px_Px;
  abcd[5] = I_ERI_G2x2z_Px_Px_Px+ABX*I_ERI_Fx2z_Px_Px_Px;
  abcd[6] = I_ERI_Gx3y_Px_Px_Px+ABX*I_ERI_F3y_Px_Px_Px;
  abcd[7] = I_ERI_Gx2yz_Px_Px_Px+ABX*I_ERI_F2yz_Px_Px_Px;
  abcd[8] = I_ERI_Gxy2z_Px_Px_Px+ABX*I_ERI_Fy2z_Px_Px_Px;
  abcd[9] = I_ERI_Gx3z_Px_Px_Px+ABX*I_ERI_F3z_Px_Px_Px;
  abcd[10] = I_ERI_G3xy_Px_Px_Px+ABY*I_ERI_F3x_Px_Px_Px;
  abcd[11] = I_ERI_G2x2y_Px_Px_Px+ABY*I_ERI_F2xy_Px_Px_Px;
  abcd[12] = I_ERI_G2xyz_Px_Px_Px+ABY*I_ERI_F2xz_Px_Px_Px;
  abcd[13] = I_ERI_Gx3y_Px_Px_Px+ABY*I_ERI_Fx2y_Px_Px_Px;
  abcd[14] = I_ERI_Gx2yz_Px_Px_Px+ABY*I_ERI_Fxyz_Px_Px_Px;
  abcd[15] = I_ERI_Gxy2z_Px_Px_Px+ABY*I_ERI_Fx2z_Px_Px_Px;
  abcd[16] = I_ERI_G4y_Px_Px_Px+ABY*I_ERI_F3y_Px_Px_Px;
  abcd[17] = I_ERI_G3yz_Px_Px_Px+ABY*I_ERI_F2yz_Px_Px_Px;
  abcd[18] = I_ERI_G2y2z_Px_Px_Px+ABY*I_ERI_Fy2z_Px_Px_Px;
  abcd[19] = I_ERI_Gy3z_Px_Px_Px+ABY*I_ERI_F3z_Px_Px_Px;
  abcd[20] = I_ERI_G3xz_Px_Px_Px+ABZ*I_ERI_F3x_Px_Px_Px;
  abcd[21] = I_ERI_G2xyz_Px_Px_Px+ABZ*I_ERI_F2xy_Px_Px_Px;
  abcd[22] = I_ERI_G2x2z_Px_Px_Px+ABZ*I_ERI_F2xz_Px_Px_Px;
  abcd[23] = I_ERI_Gx2yz_Px_Px_Px+ABZ*I_ERI_Fx2y_Px_Px_Px;
  abcd[24] = I_ERI_Gxy2z_Px_Px_Px+ABZ*I_ERI_Fxyz_Px_Px_Px;
  abcd[25] = I_ERI_Gx3z_Px_Px_Px+ABZ*I_ERI_Fx2z_Px_Px_Px;
  abcd[26] = I_ERI_G3yz_Px_Px_Px+ABZ*I_ERI_F3y_Px_Px_Px;
  abcd[27] = I_ERI_G2y2z_Px_Px_Px+ABZ*I_ERI_F2yz_Px_Px_Px;
  abcd[28] = I_ERI_Gy3z_Px_Px_Px+ABZ*I_ERI_Fy2z_Px_Px_Px;
  abcd[29] = I_ERI_G4z_Px_Px_Px+ABZ*I_ERI_F3z_Px_Px_Px;
  abcd[30] = I_ERI_G3xy_Py_Px_Px+ABY*I_ERI_F3x_Py_Px_Px;
  abcd[31] = I_ERI_G2x2y_Py_Px_Px+ABY*I_ERI_F2xy_Py_Px_Px;
  abcd[32] = I_ERI_G2xyz_Py_Px_Px+ABY*I_ERI_F2xz_Py_Px_Px;
  abcd[33] = I_ERI_Gx3y_Py_Px_Px+ABY*I_ERI_Fx2y_Py_Px_Px;
  abcd[34] = I_ERI_Gx2yz_Py_Px_Px+ABY*I_ERI_Fxyz_Py_Px_Px;
  abcd[35] = I_ERI_Gxy2z_Py_Px_Px+ABY*I_ERI_Fx2z_Py_Px_Px;
  abcd[36] = I_ERI_G4y_Py_Px_Px+ABY*I_ERI_F3y_Py_Px_Px;
  abcd[37] = I_ERI_G3yz_Py_Px_Px+ABY*I_ERI_F2yz_Py_Px_Px;
  abcd[38] = I_ERI_G2y2z_Py_Px_Px+ABY*I_ERI_Fy2z_Py_Px_Px;
  abcd[39] = I_ERI_Gy3z_Py_Px_Px+ABY*I_ERI_F3z_Py_Px_Px;
  abcd[40] = I_ERI_G3xz_Py_Px_Px+ABZ*I_ERI_F3x_Py_Px_Px;
  abcd[41] = I_ERI_G2xyz_Py_Px_Px+ABZ*I_ERI_F2xy_Py_Px_Px;
  abcd[42] = I_ERI_G2x2z_Py_Px_Px+ABZ*I_ERI_F2xz_Py_Px_Px;
  abcd[43] = I_ERI_Gx2yz_Py_Px_Px+ABZ*I_ERI_Fx2y_Py_Px_Px;
  abcd[44] = I_ERI_Gxy2z_Py_Px_Px+ABZ*I_ERI_Fxyz_Py_Px_Px;
  abcd[45] = I_ERI_Gx3z_Py_Px_Px+ABZ*I_ERI_Fx2z_Py_Px_Px;
  abcd[46] = I_ERI_G3yz_Py_Px_Px+ABZ*I_ERI_F3y_Py_Px_Px;
  abcd[47] = I_ERI_G2y2z_Py_Px_Px+ABZ*I_ERI_F2yz_Py_Px_Px;
  abcd[48] = I_ERI_Gy3z_Py_Px_Px+ABZ*I_ERI_Fy2z_Py_Px_Px;
  abcd[49] = I_ERI_G4z_Py_Px_Px+ABZ*I_ERI_F3z_Py_Px_Px;
  abcd[50] = I_ERI_G3xz_Pz_Px_Px+ABZ*I_ERI_F3x_Pz_Px_Px;
  abcd[51] = I_ERI_G2xyz_Pz_Px_Px+ABZ*I_ERI_F2xy_Pz_Px_Px;
  abcd[52] = I_ERI_G2x2z_Pz_Px_Px+ABZ*I_ERI_F2xz_Pz_Px_Px;
  abcd[53] = I_ERI_Gx2yz_Pz_Px_Px+ABZ*I_ERI_Fx2y_Pz_Px_Px;
  abcd[54] = I_ERI_Gxy2z_Pz_Px_Px+ABZ*I_ERI_Fxyz_Pz_Px_Px;
  abcd[55] = I_ERI_Gx3z_Pz_Px_Px+ABZ*I_ERI_Fx2z_Pz_Px_Px;
  abcd[56] = I_ERI_G3yz_Pz_Px_Px+ABZ*I_ERI_F3y_Pz_Px_Px;
  abcd[57] = I_ERI_G2y2z_Pz_Px_Px+ABZ*I_ERI_F2yz_Pz_Px_Px;
  abcd[58] = I_ERI_Gy3z_Pz_Px_Px+ABZ*I_ERI_Fy2z_Pz_Px_Px;
  abcd[59] = I_ERI_G4z_Pz_Px_Px+ABZ*I_ERI_F3z_Pz_Px_Px;
  abcd[60] = I_ERI_G4x_Px_Py_Px+ABX*I_ERI_F3x_Px_Py_Px;
  abcd[61] = I_ERI_G3xy_Px_Py_Px+ABX*I_ERI_F2xy_Px_Py_Px;
  abcd[62] = I_ERI_G3xz_Px_Py_Px+ABX*I_ERI_F2xz_Px_Py_Px;
  abcd[63] = I_ERI_G2x2y_Px_Py_Px+ABX*I_ERI_Fx2y_Px_Py_Px;
  abcd[64] = I_ERI_G2xyz_Px_Py_Px+ABX*I_ERI_Fxyz_Px_Py_Px;
  abcd[65] = I_ERI_G2x2z_Px_Py_Px+ABX*I_ERI_Fx2z_Px_Py_Px;
  abcd[66] = I_ERI_Gx3y_Px_Py_Px+ABX*I_ERI_F3y_Px_Py_Px;
  abcd[67] = I_ERI_Gx2yz_Px_Py_Px+ABX*I_ERI_F2yz_Px_Py_Px;
  abcd[68] = I_ERI_Gxy2z_Px_Py_Px+ABX*I_ERI_Fy2z_Px_Py_Px;
  abcd[69] = I_ERI_Gx3z_Px_Py_Px+ABX*I_ERI_F3z_Px_Py_Px;
  abcd[70] = I_ERI_G3xy_Px_Py_Px+ABY*I_ERI_F3x_Px_Py_Px;
  abcd[71] = I_ERI_G2x2y_Px_Py_Px+ABY*I_ERI_F2xy_Px_Py_Px;
  abcd[72] = I_ERI_G2xyz_Px_Py_Px+ABY*I_ERI_F2xz_Px_Py_Px;
  abcd[73] = I_ERI_Gx3y_Px_Py_Px+ABY*I_ERI_Fx2y_Px_Py_Px;
  abcd[74] = I_ERI_Gx2yz_Px_Py_Px+ABY*I_ERI_Fxyz_Px_Py_Px;
  abcd[75] = I_ERI_Gxy2z_Px_Py_Px+ABY*I_ERI_Fx2z_Px_Py_Px;
  abcd[76] = I_ERI_G4y_Px_Py_Px+ABY*I_ERI_F3y_Px_Py_Px;
  abcd[77] = I_ERI_G3yz_Px_Py_Px+ABY*I_ERI_F2yz_Px_Py_Px;
  abcd[78] = I_ERI_G2y2z_Px_Py_Px+ABY*I_ERI_Fy2z_Px_Py_Px;
  abcd[79] = I_ERI_Gy3z_Px_Py_Px+ABY*I_ERI_F3z_Px_Py_Px;
  abcd[80] = I_ERI_G3xz_Px_Py_Px+ABZ*I_ERI_F3x_Px_Py_Px;
  abcd[81] = I_ERI_G2xyz_Px_Py_Px+ABZ*I_ERI_F2xy_Px_Py_Px;
  abcd[82] = I_ERI_G2x2z_Px_Py_Px+ABZ*I_ERI_F2xz_Px_Py_Px;
  abcd[83] = I_ERI_Gx2yz_Px_Py_Px+ABZ*I_ERI_Fx2y_Px_Py_Px;
  abcd[84] = I_ERI_Gxy2z_Px_Py_Px+ABZ*I_ERI_Fxyz_Px_Py_Px;
  abcd[85] = I_ERI_Gx3z_Px_Py_Px+ABZ*I_ERI_Fx2z_Px_Py_Px;
  abcd[86] = I_ERI_G3yz_Px_Py_Px+ABZ*I_ERI_F3y_Px_Py_Px;
  abcd[87] = I_ERI_G2y2z_Px_Py_Px+ABZ*I_ERI_F2yz_Px_Py_Px;
  abcd[88] = I_ERI_Gy3z_Px_Py_Px+ABZ*I_ERI_Fy2z_Px_Py_Px;
  abcd[89] = I_ERI_G4z_Px_Py_Px+ABZ*I_ERI_F3z_Px_Py_Px;
  abcd[90] = I_ERI_G3xy_Py_Py_Px+ABY*I_ERI_F3x_Py_Py_Px;
  abcd[91] = I_ERI_G2x2y_Py_Py_Px+ABY*I_ERI_F2xy_Py_Py_Px;
  abcd[92] = I_ERI_G2xyz_Py_Py_Px+ABY*I_ERI_F2xz_Py_Py_Px;
  abcd[93] = I_ERI_Gx3y_Py_Py_Px+ABY*I_ERI_Fx2y_Py_Py_Px;
  abcd[94] = I_ERI_Gx2yz_Py_Py_Px+ABY*I_ERI_Fxyz_Py_Py_Px;
  abcd[95] = I_ERI_Gxy2z_Py_Py_Px+ABY*I_ERI_Fx2z_Py_Py_Px;
  abcd[96] = I_ERI_G4y_Py_Py_Px+ABY*I_ERI_F3y_Py_Py_Px;
  abcd[97] = I_ERI_G3yz_Py_Py_Px+ABY*I_ERI_F2yz_Py_Py_Px;
  abcd[98] = I_ERI_G2y2z_Py_Py_Px+ABY*I_ERI_Fy2z_Py_Py_Px;
  abcd[99] = I_ERI_Gy3z_Py_Py_Px+ABY*I_ERI_F3z_Py_Py_Px;
  abcd[100] = I_ERI_G3xz_Py_Py_Px+ABZ*I_ERI_F3x_Py_Py_Px;
  abcd[101] = I_ERI_G2xyz_Py_Py_Px+ABZ*I_ERI_F2xy_Py_Py_Px;
  abcd[102] = I_ERI_G2x2z_Py_Py_Px+ABZ*I_ERI_F2xz_Py_Py_Px;
  abcd[103] = I_ERI_Gx2yz_Py_Py_Px+ABZ*I_ERI_Fx2y_Py_Py_Px;
  abcd[104] = I_ERI_Gxy2z_Py_Py_Px+ABZ*I_ERI_Fxyz_Py_Py_Px;
  abcd[105] = I_ERI_Gx3z_Py_Py_Px+ABZ*I_ERI_Fx2z_Py_Py_Px;
  abcd[106] = I_ERI_G3yz_Py_Py_Px+ABZ*I_ERI_F3y_Py_Py_Px;
  abcd[107] = I_ERI_G2y2z_Py_Py_Px+ABZ*I_ERI_F2yz_Py_Py_Px;
  abcd[108] = I_ERI_Gy3z_Py_Py_Px+ABZ*I_ERI_Fy2z_Py_Py_Px;
  abcd[109] = I_ERI_G4z_Py_Py_Px+ABZ*I_ERI_F3z_Py_Py_Px;
  abcd[110] = I_ERI_G3xz_Pz_Py_Px+ABZ*I_ERI_F3x_Pz_Py_Px;
  abcd[111] = I_ERI_G2xyz_Pz_Py_Px+ABZ*I_ERI_F2xy_Pz_Py_Px;
  abcd[112] = I_ERI_G2x2z_Pz_Py_Px+ABZ*I_ERI_F2xz_Pz_Py_Px;
  abcd[113] = I_ERI_Gx2yz_Pz_Py_Px+ABZ*I_ERI_Fx2y_Pz_Py_Px;
  abcd[114] = I_ERI_Gxy2z_Pz_Py_Px+ABZ*I_ERI_Fxyz_Pz_Py_Px;
  abcd[115] = I_ERI_Gx3z_Pz_Py_Px+ABZ*I_ERI_Fx2z_Pz_Py_Px;
  abcd[116] = I_ERI_G3yz_Pz_Py_Px+ABZ*I_ERI_F3y_Pz_Py_Px;
  abcd[117] = I_ERI_G2y2z_Pz_Py_Px+ABZ*I_ERI_F2yz_Pz_Py_Px;
  abcd[118] = I_ERI_Gy3z_Pz_Py_Px+ABZ*I_ERI_Fy2z_Pz_Py_Px;
  abcd[119] = I_ERI_G4z_Pz_Py_Px+ABZ*I_ERI_F3z_Pz_Py_Px;
  abcd[120] = I_ERI_G4x_Px_Pz_Px+ABX*I_ERI_F3x_Px_Pz_Px;
  abcd[121] = I_ERI_G3xy_Px_Pz_Px+ABX*I_ERI_F2xy_Px_Pz_Px;
  abcd[122] = I_ERI_G3xz_Px_Pz_Px+ABX*I_ERI_F2xz_Px_Pz_Px;
  abcd[123] = I_ERI_G2x2y_Px_Pz_Px+ABX*I_ERI_Fx2y_Px_Pz_Px;
  abcd[124] = I_ERI_G2xyz_Px_Pz_Px+ABX*I_ERI_Fxyz_Px_Pz_Px;
  abcd[125] = I_ERI_G2x2z_Px_Pz_Px+ABX*I_ERI_Fx2z_Px_Pz_Px;
  abcd[126] = I_ERI_Gx3y_Px_Pz_Px+ABX*I_ERI_F3y_Px_Pz_Px;
  abcd[127] = I_ERI_Gx2yz_Px_Pz_Px+ABX*I_ERI_F2yz_Px_Pz_Px;
  abcd[128] = I_ERI_Gxy2z_Px_Pz_Px+ABX*I_ERI_Fy2z_Px_Pz_Px;
  abcd[129] = I_ERI_Gx3z_Px_Pz_Px+ABX*I_ERI_F3z_Px_Pz_Px;
  abcd[130] = I_ERI_G3xy_Px_Pz_Px+ABY*I_ERI_F3x_Px_Pz_Px;
  abcd[131] = I_ERI_G2x2y_Px_Pz_Px+ABY*I_ERI_F2xy_Px_Pz_Px;
  abcd[132] = I_ERI_G2xyz_Px_Pz_Px+ABY*I_ERI_F2xz_Px_Pz_Px;
  abcd[133] = I_ERI_Gx3y_Px_Pz_Px+ABY*I_ERI_Fx2y_Px_Pz_Px;
  abcd[134] = I_ERI_Gx2yz_Px_Pz_Px+ABY*I_ERI_Fxyz_Px_Pz_Px;
  abcd[135] = I_ERI_Gxy2z_Px_Pz_Px+ABY*I_ERI_Fx2z_Px_Pz_Px;
  abcd[136] = I_ERI_G4y_Px_Pz_Px+ABY*I_ERI_F3y_Px_Pz_Px;
  abcd[137] = I_ERI_G3yz_Px_Pz_Px+ABY*I_ERI_F2yz_Px_Pz_Px;
  abcd[138] = I_ERI_G2y2z_Px_Pz_Px+ABY*I_ERI_Fy2z_Px_Pz_Px;
  abcd[139] = I_ERI_Gy3z_Px_Pz_Px+ABY*I_ERI_F3z_Px_Pz_Px;
  abcd[140] = I_ERI_G3xz_Px_Pz_Px+ABZ*I_ERI_F3x_Px_Pz_Px;
  abcd[141] = I_ERI_G2xyz_Px_Pz_Px+ABZ*I_ERI_F2xy_Px_Pz_Px;
  abcd[142] = I_ERI_G2x2z_Px_Pz_Px+ABZ*I_ERI_F2xz_Px_Pz_Px;
  abcd[143] = I_ERI_Gx2yz_Px_Pz_Px+ABZ*I_ERI_Fx2y_Px_Pz_Px;
  abcd[144] = I_ERI_Gxy2z_Px_Pz_Px+ABZ*I_ERI_Fxyz_Px_Pz_Px;
  abcd[145] = I_ERI_Gx3z_Px_Pz_Px+ABZ*I_ERI_Fx2z_Px_Pz_Px;
  abcd[146] = I_ERI_G3yz_Px_Pz_Px+ABZ*I_ERI_F3y_Px_Pz_Px;
  abcd[147] = I_ERI_G2y2z_Px_Pz_Px+ABZ*I_ERI_F2yz_Px_Pz_Px;
  abcd[148] = I_ERI_Gy3z_Px_Pz_Px+ABZ*I_ERI_Fy2z_Px_Pz_Px;
  abcd[149] = I_ERI_G4z_Px_Pz_Px+ABZ*I_ERI_F3z_Px_Pz_Px;
  abcd[150] = I_ERI_G3xy_Py_Pz_Px+ABY*I_ERI_F3x_Py_Pz_Px;
  abcd[151] = I_ERI_G2x2y_Py_Pz_Px+ABY*I_ERI_F2xy_Py_Pz_Px;
  abcd[152] = I_ERI_G2xyz_Py_Pz_Px+ABY*I_ERI_F2xz_Py_Pz_Px;
  abcd[153] = I_ERI_Gx3y_Py_Pz_Px+ABY*I_ERI_Fx2y_Py_Pz_Px;
  abcd[154] = I_ERI_Gx2yz_Py_Pz_Px+ABY*I_ERI_Fxyz_Py_Pz_Px;
  abcd[155] = I_ERI_Gxy2z_Py_Pz_Px+ABY*I_ERI_Fx2z_Py_Pz_Px;
  abcd[156] = I_ERI_G4y_Py_Pz_Px+ABY*I_ERI_F3y_Py_Pz_Px;
  abcd[157] = I_ERI_G3yz_Py_Pz_Px+ABY*I_ERI_F2yz_Py_Pz_Px;
  abcd[158] = I_ERI_G2y2z_Py_Pz_Px+ABY*I_ERI_Fy2z_Py_Pz_Px;
  abcd[159] = I_ERI_Gy3z_Py_Pz_Px+ABY*I_ERI_F3z_Py_Pz_Px;
  abcd[160] = I_ERI_G3xz_Py_Pz_Px+ABZ*I_ERI_F3x_Py_Pz_Px;
  abcd[161] = I_ERI_G2xyz_Py_Pz_Px+ABZ*I_ERI_F2xy_Py_Pz_Px;
  abcd[162] = I_ERI_G2x2z_Py_Pz_Px+ABZ*I_ERI_F2xz_Py_Pz_Px;
  abcd[163] = I_ERI_Gx2yz_Py_Pz_Px+ABZ*I_ERI_Fx2y_Py_Pz_Px;
  abcd[164] = I_ERI_Gxy2z_Py_Pz_Px+ABZ*I_ERI_Fxyz_Py_Pz_Px;
  abcd[165] = I_ERI_Gx3z_Py_Pz_Px+ABZ*I_ERI_Fx2z_Py_Pz_Px;
  abcd[166] = I_ERI_G3yz_Py_Pz_Px+ABZ*I_ERI_F3y_Py_Pz_Px;
  abcd[167] = I_ERI_G2y2z_Py_Pz_Px+ABZ*I_ERI_F2yz_Py_Pz_Px;
  abcd[168] = I_ERI_Gy3z_Py_Pz_Px+ABZ*I_ERI_Fy2z_Py_Pz_Px;
  abcd[169] = I_ERI_G4z_Py_Pz_Px+ABZ*I_ERI_F3z_Py_Pz_Px;
  abcd[170] = I_ERI_G3xz_Pz_Pz_Px+ABZ*I_ERI_F3x_Pz_Pz_Px;
  abcd[171] = I_ERI_G2xyz_Pz_Pz_Px+ABZ*I_ERI_F2xy_Pz_Pz_Px;
  abcd[172] = I_ERI_G2x2z_Pz_Pz_Px+ABZ*I_ERI_F2xz_Pz_Pz_Px;
  abcd[173] = I_ERI_Gx2yz_Pz_Pz_Px+ABZ*I_ERI_Fx2y_Pz_Pz_Px;
  abcd[174] = I_ERI_Gxy2z_Pz_Pz_Px+ABZ*I_ERI_Fxyz_Pz_Pz_Px;
  abcd[175] = I_ERI_Gx3z_Pz_Pz_Px+ABZ*I_ERI_Fx2z_Pz_Pz_Px;
  abcd[176] = I_ERI_G3yz_Pz_Pz_Px+ABZ*I_ERI_F3y_Pz_Pz_Px;
  abcd[177] = I_ERI_G2y2z_Pz_Pz_Px+ABZ*I_ERI_F2yz_Pz_Pz_Px;
  abcd[178] = I_ERI_Gy3z_Pz_Pz_Px+ABZ*I_ERI_Fy2z_Pz_Pz_Px;
  abcd[179] = I_ERI_G4z_Pz_Pz_Px+ABZ*I_ERI_F3z_Pz_Pz_Px;
  abcd[180] = I_ERI_G4x_Px_Px_Py+ABX*I_ERI_F3x_Px_Px_Py;
  abcd[181] = I_ERI_G3xy_Px_Px_Py+ABX*I_ERI_F2xy_Px_Px_Py;
  abcd[182] = I_ERI_G3xz_Px_Px_Py+ABX*I_ERI_F2xz_Px_Px_Py;
  abcd[183] = I_ERI_G2x2y_Px_Px_Py+ABX*I_ERI_Fx2y_Px_Px_Py;
  abcd[184] = I_ERI_G2xyz_Px_Px_Py+ABX*I_ERI_Fxyz_Px_Px_Py;
  abcd[185] = I_ERI_G2x2z_Px_Px_Py+ABX*I_ERI_Fx2z_Px_Px_Py;
  abcd[186] = I_ERI_Gx3y_Px_Px_Py+ABX*I_ERI_F3y_Px_Px_Py;
  abcd[187] = I_ERI_Gx2yz_Px_Px_Py+ABX*I_ERI_F2yz_Px_Px_Py;
  abcd[188] = I_ERI_Gxy2z_Px_Px_Py+ABX*I_ERI_Fy2z_Px_Px_Py;
  abcd[189] = I_ERI_Gx3z_Px_Px_Py+ABX*I_ERI_F3z_Px_Px_Py;
  abcd[190] = I_ERI_G3xy_Px_Px_Py+ABY*I_ERI_F3x_Px_Px_Py;
  abcd[191] = I_ERI_G2x2y_Px_Px_Py+ABY*I_ERI_F2xy_Px_Px_Py;
  abcd[192] = I_ERI_G2xyz_Px_Px_Py+ABY*I_ERI_F2xz_Px_Px_Py;
  abcd[193] = I_ERI_Gx3y_Px_Px_Py+ABY*I_ERI_Fx2y_Px_Px_Py;
  abcd[194] = I_ERI_Gx2yz_Px_Px_Py+ABY*I_ERI_Fxyz_Px_Px_Py;
  abcd[195] = I_ERI_Gxy2z_Px_Px_Py+ABY*I_ERI_Fx2z_Px_Px_Py;
  abcd[196] = I_ERI_G4y_Px_Px_Py+ABY*I_ERI_F3y_Px_Px_Py;
  abcd[197] = I_ERI_G3yz_Px_Px_Py+ABY*I_ERI_F2yz_Px_Px_Py;
  abcd[198] = I_ERI_G2y2z_Px_Px_Py+ABY*I_ERI_Fy2z_Px_Px_Py;
  abcd[199] = I_ERI_Gy3z_Px_Px_Py+ABY*I_ERI_F3z_Px_Px_Py;
  abcd[200] = I_ERI_G3xz_Px_Px_Py+ABZ*I_ERI_F3x_Px_Px_Py;
  abcd[201] = I_ERI_G2xyz_Px_Px_Py+ABZ*I_ERI_F2xy_Px_Px_Py;
  abcd[202] = I_ERI_G2x2z_Px_Px_Py+ABZ*I_ERI_F2xz_Px_Px_Py;
  abcd[203] = I_ERI_Gx2yz_Px_Px_Py+ABZ*I_ERI_Fx2y_Px_Px_Py;
  abcd[204] = I_ERI_Gxy2z_Px_Px_Py+ABZ*I_ERI_Fxyz_Px_Px_Py;
  abcd[205] = I_ERI_Gx3z_Px_Px_Py+ABZ*I_ERI_Fx2z_Px_Px_Py;
  abcd[206] = I_ERI_G3yz_Px_Px_Py+ABZ*I_ERI_F3y_Px_Px_Py;
  abcd[207] = I_ERI_G2y2z_Px_Px_Py+ABZ*I_ERI_F2yz_Px_Px_Py;
  abcd[208] = I_ERI_Gy3z_Px_Px_Py+ABZ*I_ERI_Fy2z_Px_Px_Py;
  abcd[209] = I_ERI_G4z_Px_Px_Py+ABZ*I_ERI_F3z_Px_Px_Py;
  abcd[210] = I_ERI_G3xy_Py_Px_Py+ABY*I_ERI_F3x_Py_Px_Py;
  abcd[211] = I_ERI_G2x2y_Py_Px_Py+ABY*I_ERI_F2xy_Py_Px_Py;
  abcd[212] = I_ERI_G2xyz_Py_Px_Py+ABY*I_ERI_F2xz_Py_Px_Py;
  abcd[213] = I_ERI_Gx3y_Py_Px_Py+ABY*I_ERI_Fx2y_Py_Px_Py;
  abcd[214] = I_ERI_Gx2yz_Py_Px_Py+ABY*I_ERI_Fxyz_Py_Px_Py;
  abcd[215] = I_ERI_Gxy2z_Py_Px_Py+ABY*I_ERI_Fx2z_Py_Px_Py;
  abcd[216] = I_ERI_G4y_Py_Px_Py+ABY*I_ERI_F3y_Py_Px_Py;
  abcd[217] = I_ERI_G3yz_Py_Px_Py+ABY*I_ERI_F2yz_Py_Px_Py;
  abcd[218] = I_ERI_G2y2z_Py_Px_Py+ABY*I_ERI_Fy2z_Py_Px_Py;
  abcd[219] = I_ERI_Gy3z_Py_Px_Py+ABY*I_ERI_F3z_Py_Px_Py;
  abcd[220] = I_ERI_G3xz_Py_Px_Py+ABZ*I_ERI_F3x_Py_Px_Py;
  abcd[221] = I_ERI_G2xyz_Py_Px_Py+ABZ*I_ERI_F2xy_Py_Px_Py;
  abcd[222] = I_ERI_G2x2z_Py_Px_Py+ABZ*I_ERI_F2xz_Py_Px_Py;
  abcd[223] = I_ERI_Gx2yz_Py_Px_Py+ABZ*I_ERI_Fx2y_Py_Px_Py;
  abcd[224] = I_ERI_Gxy2z_Py_Px_Py+ABZ*I_ERI_Fxyz_Py_Px_Py;
  abcd[225] = I_ERI_Gx3z_Py_Px_Py+ABZ*I_ERI_Fx2z_Py_Px_Py;
  abcd[226] = I_ERI_G3yz_Py_Px_Py+ABZ*I_ERI_F3y_Py_Px_Py;
  abcd[227] = I_ERI_G2y2z_Py_Px_Py+ABZ*I_ERI_F2yz_Py_Px_Py;
  abcd[228] = I_ERI_Gy3z_Py_Px_Py+ABZ*I_ERI_Fy2z_Py_Px_Py;
  abcd[229] = I_ERI_G4z_Py_Px_Py+ABZ*I_ERI_F3z_Py_Px_Py;
  abcd[230] = I_ERI_G3xz_Pz_Px_Py+ABZ*I_ERI_F3x_Pz_Px_Py;
  abcd[231] = I_ERI_G2xyz_Pz_Px_Py+ABZ*I_ERI_F2xy_Pz_Px_Py;
  abcd[232] = I_ERI_G2x2z_Pz_Px_Py+ABZ*I_ERI_F2xz_Pz_Px_Py;
  abcd[233] = I_ERI_Gx2yz_Pz_Px_Py+ABZ*I_ERI_Fx2y_Pz_Px_Py;
  abcd[234] = I_ERI_Gxy2z_Pz_Px_Py+ABZ*I_ERI_Fxyz_Pz_Px_Py;
  abcd[235] = I_ERI_Gx3z_Pz_Px_Py+ABZ*I_ERI_Fx2z_Pz_Px_Py;
  abcd[236] = I_ERI_G3yz_Pz_Px_Py+ABZ*I_ERI_F3y_Pz_Px_Py;
  abcd[237] = I_ERI_G2y2z_Pz_Px_Py+ABZ*I_ERI_F2yz_Pz_Px_Py;
  abcd[238] = I_ERI_Gy3z_Pz_Px_Py+ABZ*I_ERI_Fy2z_Pz_Px_Py;
  abcd[239] = I_ERI_G4z_Pz_Px_Py+ABZ*I_ERI_F3z_Pz_Px_Py;
  abcd[240] = I_ERI_G4x_Px_Py_Py+ABX*I_ERI_F3x_Px_Py_Py;
  abcd[241] = I_ERI_G3xy_Px_Py_Py+ABX*I_ERI_F2xy_Px_Py_Py;
  abcd[242] = I_ERI_G3xz_Px_Py_Py+ABX*I_ERI_F2xz_Px_Py_Py;
  abcd[243] = I_ERI_G2x2y_Px_Py_Py+ABX*I_ERI_Fx2y_Px_Py_Py;
  abcd[244] = I_ERI_G2xyz_Px_Py_Py+ABX*I_ERI_Fxyz_Px_Py_Py;
  abcd[245] = I_ERI_G2x2z_Px_Py_Py+ABX*I_ERI_Fx2z_Px_Py_Py;
  abcd[246] = I_ERI_Gx3y_Px_Py_Py+ABX*I_ERI_F3y_Px_Py_Py;
  abcd[247] = I_ERI_Gx2yz_Px_Py_Py+ABX*I_ERI_F2yz_Px_Py_Py;
  abcd[248] = I_ERI_Gxy2z_Px_Py_Py+ABX*I_ERI_Fy2z_Px_Py_Py;
  abcd[249] = I_ERI_Gx3z_Px_Py_Py+ABX*I_ERI_F3z_Px_Py_Py;
  abcd[250] = I_ERI_G3xy_Px_Py_Py+ABY*I_ERI_F3x_Px_Py_Py;
  abcd[251] = I_ERI_G2x2y_Px_Py_Py+ABY*I_ERI_F2xy_Px_Py_Py;
  abcd[252] = I_ERI_G2xyz_Px_Py_Py+ABY*I_ERI_F2xz_Px_Py_Py;
  abcd[253] = I_ERI_Gx3y_Px_Py_Py+ABY*I_ERI_Fx2y_Px_Py_Py;
  abcd[254] = I_ERI_Gx2yz_Px_Py_Py+ABY*I_ERI_Fxyz_Px_Py_Py;
  abcd[255] = I_ERI_Gxy2z_Px_Py_Py+ABY*I_ERI_Fx2z_Px_Py_Py;
  abcd[256] = I_ERI_G4y_Px_Py_Py+ABY*I_ERI_F3y_Px_Py_Py;
  abcd[257] = I_ERI_G3yz_Px_Py_Py+ABY*I_ERI_F2yz_Px_Py_Py;
  abcd[258] = I_ERI_G2y2z_Px_Py_Py+ABY*I_ERI_Fy2z_Px_Py_Py;
  abcd[259] = I_ERI_Gy3z_Px_Py_Py+ABY*I_ERI_F3z_Px_Py_Py;
  abcd[260] = I_ERI_G3xz_Px_Py_Py+ABZ*I_ERI_F3x_Px_Py_Py;
  abcd[261] = I_ERI_G2xyz_Px_Py_Py+ABZ*I_ERI_F2xy_Px_Py_Py;
  abcd[262] = I_ERI_G2x2z_Px_Py_Py+ABZ*I_ERI_F2xz_Px_Py_Py;
  abcd[263] = I_ERI_Gx2yz_Px_Py_Py+ABZ*I_ERI_Fx2y_Px_Py_Py;
  abcd[264] = I_ERI_Gxy2z_Px_Py_Py+ABZ*I_ERI_Fxyz_Px_Py_Py;
  abcd[265] = I_ERI_Gx3z_Px_Py_Py+ABZ*I_ERI_Fx2z_Px_Py_Py;
  abcd[266] = I_ERI_G3yz_Px_Py_Py+ABZ*I_ERI_F3y_Px_Py_Py;
  abcd[267] = I_ERI_G2y2z_Px_Py_Py+ABZ*I_ERI_F2yz_Px_Py_Py;
  abcd[268] = I_ERI_Gy3z_Px_Py_Py+ABZ*I_ERI_Fy2z_Px_Py_Py;
  abcd[269] = I_ERI_G4z_Px_Py_Py+ABZ*I_ERI_F3z_Px_Py_Py;
  abcd[270] = I_ERI_G3xy_Py_Py_Py+ABY*I_ERI_F3x_Py_Py_Py;
  abcd[271] = I_ERI_G2x2y_Py_Py_Py+ABY*I_ERI_F2xy_Py_Py_Py;
  abcd[272] = I_ERI_G2xyz_Py_Py_Py+ABY*I_ERI_F2xz_Py_Py_Py;
  abcd[273] = I_ERI_Gx3y_Py_Py_Py+ABY*I_ERI_Fx2y_Py_Py_Py;
  abcd[274] = I_ERI_Gx2yz_Py_Py_Py+ABY*I_ERI_Fxyz_Py_Py_Py;
  abcd[275] = I_ERI_Gxy2z_Py_Py_Py+ABY*I_ERI_Fx2z_Py_Py_Py;
  abcd[276] = I_ERI_G4y_Py_Py_Py+ABY*I_ERI_F3y_Py_Py_Py;
  abcd[277] = I_ERI_G3yz_Py_Py_Py+ABY*I_ERI_F2yz_Py_Py_Py;
  abcd[278] = I_ERI_G2y2z_Py_Py_Py+ABY*I_ERI_Fy2z_Py_Py_Py;
  abcd[279] = I_ERI_Gy3z_Py_Py_Py+ABY*I_ERI_F3z_Py_Py_Py;
  abcd[280] = I_ERI_G3xz_Py_Py_Py+ABZ*I_ERI_F3x_Py_Py_Py;
  abcd[281] = I_ERI_G2xyz_Py_Py_Py+ABZ*I_ERI_F2xy_Py_Py_Py;
  abcd[282] = I_ERI_G2x2z_Py_Py_Py+ABZ*I_ERI_F2xz_Py_Py_Py;
  abcd[283] = I_ERI_Gx2yz_Py_Py_Py+ABZ*I_ERI_Fx2y_Py_Py_Py;
  abcd[284] = I_ERI_Gxy2z_Py_Py_Py+ABZ*I_ERI_Fxyz_Py_Py_Py;
  abcd[285] = I_ERI_Gx3z_Py_Py_Py+ABZ*I_ERI_Fx2z_Py_Py_Py;
  abcd[286] = I_ERI_G3yz_Py_Py_Py+ABZ*I_ERI_F3y_Py_Py_Py;
  abcd[287] = I_ERI_G2y2z_Py_Py_Py+ABZ*I_ERI_F2yz_Py_Py_Py;
  abcd[288] = I_ERI_Gy3z_Py_Py_Py+ABZ*I_ERI_Fy2z_Py_Py_Py;
  abcd[289] = I_ERI_G4z_Py_Py_Py+ABZ*I_ERI_F3z_Py_Py_Py;
  abcd[290] = I_ERI_G3xz_Pz_Py_Py+ABZ*I_ERI_F3x_Pz_Py_Py;
  abcd[291] = I_ERI_G2xyz_Pz_Py_Py+ABZ*I_ERI_F2xy_Pz_Py_Py;
  abcd[292] = I_ERI_G2x2z_Pz_Py_Py+ABZ*I_ERI_F2xz_Pz_Py_Py;
  abcd[293] = I_ERI_Gx2yz_Pz_Py_Py+ABZ*I_ERI_Fx2y_Pz_Py_Py;
  abcd[294] = I_ERI_Gxy2z_Pz_Py_Py+ABZ*I_ERI_Fxyz_Pz_Py_Py;
  abcd[295] = I_ERI_Gx3z_Pz_Py_Py+ABZ*I_ERI_Fx2z_Pz_Py_Py;
  abcd[296] = I_ERI_G3yz_Pz_Py_Py+ABZ*I_ERI_F3y_Pz_Py_Py;
  abcd[297] = I_ERI_G2y2z_Pz_Py_Py+ABZ*I_ERI_F2yz_Pz_Py_Py;
  abcd[298] = I_ERI_Gy3z_Pz_Py_Py+ABZ*I_ERI_Fy2z_Pz_Py_Py;
  abcd[299] = I_ERI_G4z_Pz_Py_Py+ABZ*I_ERI_F3z_Pz_Py_Py;
  abcd[300] = I_ERI_G4x_Px_Pz_Py+ABX*I_ERI_F3x_Px_Pz_Py;
  abcd[301] = I_ERI_G3xy_Px_Pz_Py+ABX*I_ERI_F2xy_Px_Pz_Py;
  abcd[302] = I_ERI_G3xz_Px_Pz_Py+ABX*I_ERI_F2xz_Px_Pz_Py;
  abcd[303] = I_ERI_G2x2y_Px_Pz_Py+ABX*I_ERI_Fx2y_Px_Pz_Py;
  abcd[304] = I_ERI_G2xyz_Px_Pz_Py+ABX*I_ERI_Fxyz_Px_Pz_Py;
  abcd[305] = I_ERI_G2x2z_Px_Pz_Py+ABX*I_ERI_Fx2z_Px_Pz_Py;
  abcd[306] = I_ERI_Gx3y_Px_Pz_Py+ABX*I_ERI_F3y_Px_Pz_Py;
  abcd[307] = I_ERI_Gx2yz_Px_Pz_Py+ABX*I_ERI_F2yz_Px_Pz_Py;
  abcd[308] = I_ERI_Gxy2z_Px_Pz_Py+ABX*I_ERI_Fy2z_Px_Pz_Py;
  abcd[309] = I_ERI_Gx3z_Px_Pz_Py+ABX*I_ERI_F3z_Px_Pz_Py;
  abcd[310] = I_ERI_G3xy_Px_Pz_Py+ABY*I_ERI_F3x_Px_Pz_Py;
  abcd[311] = I_ERI_G2x2y_Px_Pz_Py+ABY*I_ERI_F2xy_Px_Pz_Py;
  abcd[312] = I_ERI_G2xyz_Px_Pz_Py+ABY*I_ERI_F2xz_Px_Pz_Py;
  abcd[313] = I_ERI_Gx3y_Px_Pz_Py+ABY*I_ERI_Fx2y_Px_Pz_Py;
  abcd[314] = I_ERI_Gx2yz_Px_Pz_Py+ABY*I_ERI_Fxyz_Px_Pz_Py;
  abcd[315] = I_ERI_Gxy2z_Px_Pz_Py+ABY*I_ERI_Fx2z_Px_Pz_Py;
  abcd[316] = I_ERI_G4y_Px_Pz_Py+ABY*I_ERI_F3y_Px_Pz_Py;
  abcd[317] = I_ERI_G3yz_Px_Pz_Py+ABY*I_ERI_F2yz_Px_Pz_Py;
  abcd[318] = I_ERI_G2y2z_Px_Pz_Py+ABY*I_ERI_Fy2z_Px_Pz_Py;
  abcd[319] = I_ERI_Gy3z_Px_Pz_Py+ABY*I_ERI_F3z_Px_Pz_Py;
  abcd[320] = I_ERI_G3xz_Px_Pz_Py+ABZ*I_ERI_F3x_Px_Pz_Py;
  abcd[321] = I_ERI_G2xyz_Px_Pz_Py+ABZ*I_ERI_F2xy_Px_Pz_Py;
  abcd[322] = I_ERI_G2x2z_Px_Pz_Py+ABZ*I_ERI_F2xz_Px_Pz_Py;
  abcd[323] = I_ERI_Gx2yz_Px_Pz_Py+ABZ*I_ERI_Fx2y_Px_Pz_Py;
  abcd[324] = I_ERI_Gxy2z_Px_Pz_Py+ABZ*I_ERI_Fxyz_Px_Pz_Py;
  abcd[325] = I_ERI_Gx3z_Px_Pz_Py+ABZ*I_ERI_Fx2z_Px_Pz_Py;
  abcd[326] = I_ERI_G3yz_Px_Pz_Py+ABZ*I_ERI_F3y_Px_Pz_Py;
  abcd[327] = I_ERI_G2y2z_Px_Pz_Py+ABZ*I_ERI_F2yz_Px_Pz_Py;
  abcd[328] = I_ERI_Gy3z_Px_Pz_Py+ABZ*I_ERI_Fy2z_Px_Pz_Py;
  abcd[329] = I_ERI_G4z_Px_Pz_Py+ABZ*I_ERI_F3z_Px_Pz_Py;
  abcd[330] = I_ERI_G3xy_Py_Pz_Py+ABY*I_ERI_F3x_Py_Pz_Py;
  abcd[331] = I_ERI_G2x2y_Py_Pz_Py+ABY*I_ERI_F2xy_Py_Pz_Py;
  abcd[332] = I_ERI_G2xyz_Py_Pz_Py+ABY*I_ERI_F2xz_Py_Pz_Py;
  abcd[333] = I_ERI_Gx3y_Py_Pz_Py+ABY*I_ERI_Fx2y_Py_Pz_Py;
  abcd[334] = I_ERI_Gx2yz_Py_Pz_Py+ABY*I_ERI_Fxyz_Py_Pz_Py;
  abcd[335] = I_ERI_Gxy2z_Py_Pz_Py+ABY*I_ERI_Fx2z_Py_Pz_Py;
  abcd[336] = I_ERI_G4y_Py_Pz_Py+ABY*I_ERI_F3y_Py_Pz_Py;
  abcd[337] = I_ERI_G3yz_Py_Pz_Py+ABY*I_ERI_F2yz_Py_Pz_Py;
  abcd[338] = I_ERI_G2y2z_Py_Pz_Py+ABY*I_ERI_Fy2z_Py_Pz_Py;
  abcd[339] = I_ERI_Gy3z_Py_Pz_Py+ABY*I_ERI_F3z_Py_Pz_Py;
  abcd[340] = I_ERI_G3xz_Py_Pz_Py+ABZ*I_ERI_F3x_Py_Pz_Py;
  abcd[341] = I_ERI_G2xyz_Py_Pz_Py+ABZ*I_ERI_F2xy_Py_Pz_Py;
  abcd[342] = I_ERI_G2x2z_Py_Pz_Py+ABZ*I_ERI_F2xz_Py_Pz_Py;
  abcd[343] = I_ERI_Gx2yz_Py_Pz_Py+ABZ*I_ERI_Fx2y_Py_Pz_Py;
  abcd[344] = I_ERI_Gxy2z_Py_Pz_Py+ABZ*I_ERI_Fxyz_Py_Pz_Py;
  abcd[345] = I_ERI_Gx3z_Py_Pz_Py+ABZ*I_ERI_Fx2z_Py_Pz_Py;
  abcd[346] = I_ERI_G3yz_Py_Pz_Py+ABZ*I_ERI_F3y_Py_Pz_Py;
  abcd[347] = I_ERI_G2y2z_Py_Pz_Py+ABZ*I_ERI_F2yz_Py_Pz_Py;
  abcd[348] = I_ERI_Gy3z_Py_Pz_Py+ABZ*I_ERI_Fy2z_Py_Pz_Py;
  abcd[349] = I_ERI_G4z_Py_Pz_Py+ABZ*I_ERI_F3z_Py_Pz_Py;
  abcd[350] = I_ERI_G3xz_Pz_Pz_Py+ABZ*I_ERI_F3x_Pz_Pz_Py;
  abcd[351] = I_ERI_G2xyz_Pz_Pz_Py+ABZ*I_ERI_F2xy_Pz_Pz_Py;
  abcd[352] = I_ERI_G2x2z_Pz_Pz_Py+ABZ*I_ERI_F2xz_Pz_Pz_Py;
  abcd[353] = I_ERI_Gx2yz_Pz_Pz_Py+ABZ*I_ERI_Fx2y_Pz_Pz_Py;
  abcd[354] = I_ERI_Gxy2z_Pz_Pz_Py+ABZ*I_ERI_Fxyz_Pz_Pz_Py;
  abcd[355] = I_ERI_Gx3z_Pz_Pz_Py+ABZ*I_ERI_Fx2z_Pz_Pz_Py;
  abcd[356] = I_ERI_G3yz_Pz_Pz_Py+ABZ*I_ERI_F3y_Pz_Pz_Py;
  abcd[357] = I_ERI_G2y2z_Pz_Pz_Py+ABZ*I_ERI_F2yz_Pz_Pz_Py;
  abcd[358] = I_ERI_Gy3z_Pz_Pz_Py+ABZ*I_ERI_Fy2z_Pz_Pz_Py;
  abcd[359] = I_ERI_G4z_Pz_Pz_Py+ABZ*I_ERI_F3z_Pz_Pz_Py;
  abcd[360] = I_ERI_G4x_Px_Px_Pz+ABX*I_ERI_F3x_Px_Px_Pz;
  abcd[361] = I_ERI_G3xy_Px_Px_Pz+ABX*I_ERI_F2xy_Px_Px_Pz;
  abcd[362] = I_ERI_G3xz_Px_Px_Pz+ABX*I_ERI_F2xz_Px_Px_Pz;
  abcd[363] = I_ERI_G2x2y_Px_Px_Pz+ABX*I_ERI_Fx2y_Px_Px_Pz;
  abcd[364] = I_ERI_G2xyz_Px_Px_Pz+ABX*I_ERI_Fxyz_Px_Px_Pz;
  abcd[365] = I_ERI_G2x2z_Px_Px_Pz+ABX*I_ERI_Fx2z_Px_Px_Pz;
  abcd[366] = I_ERI_Gx3y_Px_Px_Pz+ABX*I_ERI_F3y_Px_Px_Pz;
  abcd[367] = I_ERI_Gx2yz_Px_Px_Pz+ABX*I_ERI_F2yz_Px_Px_Pz;
  abcd[368] = I_ERI_Gxy2z_Px_Px_Pz+ABX*I_ERI_Fy2z_Px_Px_Pz;
  abcd[369] = I_ERI_Gx3z_Px_Px_Pz+ABX*I_ERI_F3z_Px_Px_Pz;
  abcd[370] = I_ERI_G3xy_Px_Px_Pz+ABY*I_ERI_F3x_Px_Px_Pz;
  abcd[371] = I_ERI_G2x2y_Px_Px_Pz+ABY*I_ERI_F2xy_Px_Px_Pz;
  abcd[372] = I_ERI_G2xyz_Px_Px_Pz+ABY*I_ERI_F2xz_Px_Px_Pz;
  abcd[373] = I_ERI_Gx3y_Px_Px_Pz+ABY*I_ERI_Fx2y_Px_Px_Pz;
  abcd[374] = I_ERI_Gx2yz_Px_Px_Pz+ABY*I_ERI_Fxyz_Px_Px_Pz;
  abcd[375] = I_ERI_Gxy2z_Px_Px_Pz+ABY*I_ERI_Fx2z_Px_Px_Pz;
  abcd[376] = I_ERI_G4y_Px_Px_Pz+ABY*I_ERI_F3y_Px_Px_Pz;
  abcd[377] = I_ERI_G3yz_Px_Px_Pz+ABY*I_ERI_F2yz_Px_Px_Pz;
  abcd[378] = I_ERI_G2y2z_Px_Px_Pz+ABY*I_ERI_Fy2z_Px_Px_Pz;
  abcd[379] = I_ERI_Gy3z_Px_Px_Pz+ABY*I_ERI_F3z_Px_Px_Pz;
  abcd[380] = I_ERI_G3xz_Px_Px_Pz+ABZ*I_ERI_F3x_Px_Px_Pz;
  abcd[381] = I_ERI_G2xyz_Px_Px_Pz+ABZ*I_ERI_F2xy_Px_Px_Pz;
  abcd[382] = I_ERI_G2x2z_Px_Px_Pz+ABZ*I_ERI_F2xz_Px_Px_Pz;
  abcd[383] = I_ERI_Gx2yz_Px_Px_Pz+ABZ*I_ERI_Fx2y_Px_Px_Pz;
  abcd[384] = I_ERI_Gxy2z_Px_Px_Pz+ABZ*I_ERI_Fxyz_Px_Px_Pz;
  abcd[385] = I_ERI_Gx3z_Px_Px_Pz+ABZ*I_ERI_Fx2z_Px_Px_Pz;
  abcd[386] = I_ERI_G3yz_Px_Px_Pz+ABZ*I_ERI_F3y_Px_Px_Pz;
  abcd[387] = I_ERI_G2y2z_Px_Px_Pz+ABZ*I_ERI_F2yz_Px_Px_Pz;
  abcd[388] = I_ERI_Gy3z_Px_Px_Pz+ABZ*I_ERI_Fy2z_Px_Px_Pz;
  abcd[389] = I_ERI_G4z_Px_Px_Pz+ABZ*I_ERI_F3z_Px_Px_Pz;
  abcd[390] = I_ERI_G3xy_Py_Px_Pz+ABY*I_ERI_F3x_Py_Px_Pz;
  abcd[391] = I_ERI_G2x2y_Py_Px_Pz+ABY*I_ERI_F2xy_Py_Px_Pz;
  abcd[392] = I_ERI_G2xyz_Py_Px_Pz+ABY*I_ERI_F2xz_Py_Px_Pz;
  abcd[393] = I_ERI_Gx3y_Py_Px_Pz+ABY*I_ERI_Fx2y_Py_Px_Pz;
  abcd[394] = I_ERI_Gx2yz_Py_Px_Pz+ABY*I_ERI_Fxyz_Py_Px_Pz;
  abcd[395] = I_ERI_Gxy2z_Py_Px_Pz+ABY*I_ERI_Fx2z_Py_Px_Pz;
  abcd[396] = I_ERI_G4y_Py_Px_Pz+ABY*I_ERI_F3y_Py_Px_Pz;
  abcd[397] = I_ERI_G3yz_Py_Px_Pz+ABY*I_ERI_F2yz_Py_Px_Pz;
  abcd[398] = I_ERI_G2y2z_Py_Px_Pz+ABY*I_ERI_Fy2z_Py_Px_Pz;
  abcd[399] = I_ERI_Gy3z_Py_Px_Pz+ABY*I_ERI_F3z_Py_Px_Pz;
  abcd[400] = I_ERI_G3xz_Py_Px_Pz+ABZ*I_ERI_F3x_Py_Px_Pz;
  abcd[401] = I_ERI_G2xyz_Py_Px_Pz+ABZ*I_ERI_F2xy_Py_Px_Pz;
  abcd[402] = I_ERI_G2x2z_Py_Px_Pz+ABZ*I_ERI_F2xz_Py_Px_Pz;
  abcd[403] = I_ERI_Gx2yz_Py_Px_Pz+ABZ*I_ERI_Fx2y_Py_Px_Pz;
  abcd[404] = I_ERI_Gxy2z_Py_Px_Pz+ABZ*I_ERI_Fxyz_Py_Px_Pz;
  abcd[405] = I_ERI_Gx3z_Py_Px_Pz+ABZ*I_ERI_Fx2z_Py_Px_Pz;
  abcd[406] = I_ERI_G3yz_Py_Px_Pz+ABZ*I_ERI_F3y_Py_Px_Pz;
  abcd[407] = I_ERI_G2y2z_Py_Px_Pz+ABZ*I_ERI_F2yz_Py_Px_Pz;
  abcd[408] = I_ERI_Gy3z_Py_Px_Pz+ABZ*I_ERI_Fy2z_Py_Px_Pz;
  abcd[409] = I_ERI_G4z_Py_Px_Pz+ABZ*I_ERI_F3z_Py_Px_Pz;
  abcd[410] = I_ERI_G3xz_Pz_Px_Pz+ABZ*I_ERI_F3x_Pz_Px_Pz;
  abcd[411] = I_ERI_G2xyz_Pz_Px_Pz+ABZ*I_ERI_F2xy_Pz_Px_Pz;
  abcd[412] = I_ERI_G2x2z_Pz_Px_Pz+ABZ*I_ERI_F2xz_Pz_Px_Pz;
  abcd[413] = I_ERI_Gx2yz_Pz_Px_Pz+ABZ*I_ERI_Fx2y_Pz_Px_Pz;
  abcd[414] = I_ERI_Gxy2z_Pz_Px_Pz+ABZ*I_ERI_Fxyz_Pz_Px_Pz;
  abcd[415] = I_ERI_Gx3z_Pz_Px_Pz+ABZ*I_ERI_Fx2z_Pz_Px_Pz;
  abcd[416] = I_ERI_G3yz_Pz_Px_Pz+ABZ*I_ERI_F3y_Pz_Px_Pz;
  abcd[417] = I_ERI_G2y2z_Pz_Px_Pz+ABZ*I_ERI_F2yz_Pz_Px_Pz;
  abcd[418] = I_ERI_Gy3z_Pz_Px_Pz+ABZ*I_ERI_Fy2z_Pz_Px_Pz;
  abcd[419] = I_ERI_G4z_Pz_Px_Pz+ABZ*I_ERI_F3z_Pz_Px_Pz;
  abcd[420] = I_ERI_G4x_Px_Py_Pz+ABX*I_ERI_F3x_Px_Py_Pz;
  abcd[421] = I_ERI_G3xy_Px_Py_Pz+ABX*I_ERI_F2xy_Px_Py_Pz;
  abcd[422] = I_ERI_G3xz_Px_Py_Pz+ABX*I_ERI_F2xz_Px_Py_Pz;
  abcd[423] = I_ERI_G2x2y_Px_Py_Pz+ABX*I_ERI_Fx2y_Px_Py_Pz;
  abcd[424] = I_ERI_G2xyz_Px_Py_Pz+ABX*I_ERI_Fxyz_Px_Py_Pz;
  abcd[425] = I_ERI_G2x2z_Px_Py_Pz+ABX*I_ERI_Fx2z_Px_Py_Pz;
  abcd[426] = I_ERI_Gx3y_Px_Py_Pz+ABX*I_ERI_F3y_Px_Py_Pz;
  abcd[427] = I_ERI_Gx2yz_Px_Py_Pz+ABX*I_ERI_F2yz_Px_Py_Pz;
  abcd[428] = I_ERI_Gxy2z_Px_Py_Pz+ABX*I_ERI_Fy2z_Px_Py_Pz;
  abcd[429] = I_ERI_Gx3z_Px_Py_Pz+ABX*I_ERI_F3z_Px_Py_Pz;
  abcd[430] = I_ERI_G3xy_Px_Py_Pz+ABY*I_ERI_F3x_Px_Py_Pz;
  abcd[431] = I_ERI_G2x2y_Px_Py_Pz+ABY*I_ERI_F2xy_Px_Py_Pz;
  abcd[432] = I_ERI_G2xyz_Px_Py_Pz+ABY*I_ERI_F2xz_Px_Py_Pz;
  abcd[433] = I_ERI_Gx3y_Px_Py_Pz+ABY*I_ERI_Fx2y_Px_Py_Pz;
  abcd[434] = I_ERI_Gx2yz_Px_Py_Pz+ABY*I_ERI_Fxyz_Px_Py_Pz;
  abcd[435] = I_ERI_Gxy2z_Px_Py_Pz+ABY*I_ERI_Fx2z_Px_Py_Pz;
  abcd[436] = I_ERI_G4y_Px_Py_Pz+ABY*I_ERI_F3y_Px_Py_Pz;
  abcd[437] = I_ERI_G3yz_Px_Py_Pz+ABY*I_ERI_F2yz_Px_Py_Pz;
  abcd[438] = I_ERI_G2y2z_Px_Py_Pz+ABY*I_ERI_Fy2z_Px_Py_Pz;
  abcd[439] = I_ERI_Gy3z_Px_Py_Pz+ABY*I_ERI_F3z_Px_Py_Pz;
  abcd[440] = I_ERI_G3xz_Px_Py_Pz+ABZ*I_ERI_F3x_Px_Py_Pz;
  abcd[441] = I_ERI_G2xyz_Px_Py_Pz+ABZ*I_ERI_F2xy_Px_Py_Pz;
  abcd[442] = I_ERI_G2x2z_Px_Py_Pz+ABZ*I_ERI_F2xz_Px_Py_Pz;
  abcd[443] = I_ERI_Gx2yz_Px_Py_Pz+ABZ*I_ERI_Fx2y_Px_Py_Pz;
  abcd[444] = I_ERI_Gxy2z_Px_Py_Pz+ABZ*I_ERI_Fxyz_Px_Py_Pz;
  abcd[445] = I_ERI_Gx3z_Px_Py_Pz+ABZ*I_ERI_Fx2z_Px_Py_Pz;
  abcd[446] = I_ERI_G3yz_Px_Py_Pz+ABZ*I_ERI_F3y_Px_Py_Pz;
  abcd[447] = I_ERI_G2y2z_Px_Py_Pz+ABZ*I_ERI_F2yz_Px_Py_Pz;
  abcd[448] = I_ERI_Gy3z_Px_Py_Pz+ABZ*I_ERI_Fy2z_Px_Py_Pz;
  abcd[449] = I_ERI_G4z_Px_Py_Pz+ABZ*I_ERI_F3z_Px_Py_Pz;
  abcd[450] = I_ERI_G3xy_Py_Py_Pz+ABY*I_ERI_F3x_Py_Py_Pz;
  abcd[451] = I_ERI_G2x2y_Py_Py_Pz+ABY*I_ERI_F2xy_Py_Py_Pz;
  abcd[452] = I_ERI_G2xyz_Py_Py_Pz+ABY*I_ERI_F2xz_Py_Py_Pz;
  abcd[453] = I_ERI_Gx3y_Py_Py_Pz+ABY*I_ERI_Fx2y_Py_Py_Pz;
  abcd[454] = I_ERI_Gx2yz_Py_Py_Pz+ABY*I_ERI_Fxyz_Py_Py_Pz;
  abcd[455] = I_ERI_Gxy2z_Py_Py_Pz+ABY*I_ERI_Fx2z_Py_Py_Pz;
  abcd[456] = I_ERI_G4y_Py_Py_Pz+ABY*I_ERI_F3y_Py_Py_Pz;
  abcd[457] = I_ERI_G3yz_Py_Py_Pz+ABY*I_ERI_F2yz_Py_Py_Pz;
  abcd[458] = I_ERI_G2y2z_Py_Py_Pz+ABY*I_ERI_Fy2z_Py_Py_Pz;
  abcd[459] = I_ERI_Gy3z_Py_Py_Pz+ABY*I_ERI_F3z_Py_Py_Pz;
  abcd[460] = I_ERI_G3xz_Py_Py_Pz+ABZ*I_ERI_F3x_Py_Py_Pz;
  abcd[461] = I_ERI_G2xyz_Py_Py_Pz+ABZ*I_ERI_F2xy_Py_Py_Pz;
  abcd[462] = I_ERI_G2x2z_Py_Py_Pz+ABZ*I_ERI_F2xz_Py_Py_Pz;
  abcd[463] = I_ERI_Gx2yz_Py_Py_Pz+ABZ*I_ERI_Fx2y_Py_Py_Pz;
  abcd[464] = I_ERI_Gxy2z_Py_Py_Pz+ABZ*I_ERI_Fxyz_Py_Py_Pz;
  abcd[465] = I_ERI_Gx3z_Py_Py_Pz+ABZ*I_ERI_Fx2z_Py_Py_Pz;
  abcd[466] = I_ERI_G3yz_Py_Py_Pz+ABZ*I_ERI_F3y_Py_Py_Pz;
  abcd[467] = I_ERI_G2y2z_Py_Py_Pz+ABZ*I_ERI_F2yz_Py_Py_Pz;
  abcd[468] = I_ERI_Gy3z_Py_Py_Pz+ABZ*I_ERI_Fy2z_Py_Py_Pz;
  abcd[469] = I_ERI_G4z_Py_Py_Pz+ABZ*I_ERI_F3z_Py_Py_Pz;
  abcd[470] = I_ERI_G3xz_Pz_Py_Pz+ABZ*I_ERI_F3x_Pz_Py_Pz;
  abcd[471] = I_ERI_G2xyz_Pz_Py_Pz+ABZ*I_ERI_F2xy_Pz_Py_Pz;
  abcd[472] = I_ERI_G2x2z_Pz_Py_Pz+ABZ*I_ERI_F2xz_Pz_Py_Pz;
  abcd[473] = I_ERI_Gx2yz_Pz_Py_Pz+ABZ*I_ERI_Fx2y_Pz_Py_Pz;
  abcd[474] = I_ERI_Gxy2z_Pz_Py_Pz+ABZ*I_ERI_Fxyz_Pz_Py_Pz;
  abcd[475] = I_ERI_Gx3z_Pz_Py_Pz+ABZ*I_ERI_Fx2z_Pz_Py_Pz;
  abcd[476] = I_ERI_G3yz_Pz_Py_Pz+ABZ*I_ERI_F3y_Pz_Py_Pz;
  abcd[477] = I_ERI_G2y2z_Pz_Py_Pz+ABZ*I_ERI_F2yz_Pz_Py_Pz;
  abcd[478] = I_ERI_Gy3z_Pz_Py_Pz+ABZ*I_ERI_Fy2z_Pz_Py_Pz;
  abcd[479] = I_ERI_G4z_Pz_Py_Pz+ABZ*I_ERI_F3z_Pz_Py_Pz;
  abcd[480] = I_ERI_G4x_Px_Pz_Pz+ABX*I_ERI_F3x_Px_Pz_Pz;
  abcd[481] = I_ERI_G3xy_Px_Pz_Pz+ABX*I_ERI_F2xy_Px_Pz_Pz;
  abcd[482] = I_ERI_G3xz_Px_Pz_Pz+ABX*I_ERI_F2xz_Px_Pz_Pz;
  abcd[483] = I_ERI_G2x2y_Px_Pz_Pz+ABX*I_ERI_Fx2y_Px_Pz_Pz;
  abcd[484] = I_ERI_G2xyz_Px_Pz_Pz+ABX*I_ERI_Fxyz_Px_Pz_Pz;
  abcd[485] = I_ERI_G2x2z_Px_Pz_Pz+ABX*I_ERI_Fx2z_Px_Pz_Pz;
  abcd[486] = I_ERI_Gx3y_Px_Pz_Pz+ABX*I_ERI_F3y_Px_Pz_Pz;
  abcd[487] = I_ERI_Gx2yz_Px_Pz_Pz+ABX*I_ERI_F2yz_Px_Pz_Pz;
  abcd[488] = I_ERI_Gxy2z_Px_Pz_Pz+ABX*I_ERI_Fy2z_Px_Pz_Pz;
  abcd[489] = I_ERI_Gx3z_Px_Pz_Pz+ABX*I_ERI_F3z_Px_Pz_Pz;
  abcd[490] = I_ERI_G3xy_Px_Pz_Pz+ABY*I_ERI_F3x_Px_Pz_Pz;
  abcd[491] = I_ERI_G2x2y_Px_Pz_Pz+ABY*I_ERI_F2xy_Px_Pz_Pz;
  abcd[492] = I_ERI_G2xyz_Px_Pz_Pz+ABY*I_ERI_F2xz_Px_Pz_Pz;
  abcd[493] = I_ERI_Gx3y_Px_Pz_Pz+ABY*I_ERI_Fx2y_Px_Pz_Pz;
  abcd[494] = I_ERI_Gx2yz_Px_Pz_Pz+ABY*I_ERI_Fxyz_Px_Pz_Pz;
  abcd[495] = I_ERI_Gxy2z_Px_Pz_Pz+ABY*I_ERI_Fx2z_Px_Pz_Pz;
  abcd[496] = I_ERI_G4y_Px_Pz_Pz+ABY*I_ERI_F3y_Px_Pz_Pz;
  abcd[497] = I_ERI_G3yz_Px_Pz_Pz+ABY*I_ERI_F2yz_Px_Pz_Pz;
  abcd[498] = I_ERI_G2y2z_Px_Pz_Pz+ABY*I_ERI_Fy2z_Px_Pz_Pz;
  abcd[499] = I_ERI_Gy3z_Px_Pz_Pz+ABY*I_ERI_F3z_Px_Pz_Pz;
  abcd[500] = I_ERI_G3xz_Px_Pz_Pz+ABZ*I_ERI_F3x_Px_Pz_Pz;
  abcd[501] = I_ERI_G2xyz_Px_Pz_Pz+ABZ*I_ERI_F2xy_Px_Pz_Pz;
  abcd[502] = I_ERI_G2x2z_Px_Pz_Pz+ABZ*I_ERI_F2xz_Px_Pz_Pz;
  abcd[503] = I_ERI_Gx2yz_Px_Pz_Pz+ABZ*I_ERI_Fx2y_Px_Pz_Pz;
  abcd[504] = I_ERI_Gxy2z_Px_Pz_Pz+ABZ*I_ERI_Fxyz_Px_Pz_Pz;
  abcd[505] = I_ERI_Gx3z_Px_Pz_Pz+ABZ*I_ERI_Fx2z_Px_Pz_Pz;
  abcd[506] = I_ERI_G3yz_Px_Pz_Pz+ABZ*I_ERI_F3y_Px_Pz_Pz;
  abcd[507] = I_ERI_G2y2z_Px_Pz_Pz+ABZ*I_ERI_F2yz_Px_Pz_Pz;
  abcd[508] = I_ERI_Gy3z_Px_Pz_Pz+ABZ*I_ERI_Fy2z_Px_Pz_Pz;
  abcd[509] = I_ERI_G4z_Px_Pz_Pz+ABZ*I_ERI_F3z_Px_Pz_Pz;
  abcd[510] = I_ERI_G3xy_Py_Pz_Pz+ABY*I_ERI_F3x_Py_Pz_Pz;
  abcd[511] = I_ERI_G2x2y_Py_Pz_Pz+ABY*I_ERI_F2xy_Py_Pz_Pz;
  abcd[512] = I_ERI_G2xyz_Py_Pz_Pz+ABY*I_ERI_F2xz_Py_Pz_Pz;
  abcd[513] = I_ERI_Gx3y_Py_Pz_Pz+ABY*I_ERI_Fx2y_Py_Pz_Pz;
  abcd[514] = I_ERI_Gx2yz_Py_Pz_Pz+ABY*I_ERI_Fxyz_Py_Pz_Pz;
  abcd[515] = I_ERI_Gxy2z_Py_Pz_Pz+ABY*I_ERI_Fx2z_Py_Pz_Pz;
  abcd[516] = I_ERI_G4y_Py_Pz_Pz+ABY*I_ERI_F3y_Py_Pz_Pz;
  abcd[517] = I_ERI_G3yz_Py_Pz_Pz+ABY*I_ERI_F2yz_Py_Pz_Pz;
  abcd[518] = I_ERI_G2y2z_Py_Pz_Pz+ABY*I_ERI_Fy2z_Py_Pz_Pz;
  abcd[519] = I_ERI_Gy3z_Py_Pz_Pz+ABY*I_ERI_F3z_Py_Pz_Pz;
  abcd[520] = I_ERI_G3xz_Py_Pz_Pz+ABZ*I_ERI_F3x_Py_Pz_Pz;
  abcd[521] = I_ERI_G2xyz_Py_Pz_Pz+ABZ*I_ERI_F2xy_Py_Pz_Pz;
  abcd[522] = I_ERI_G2x2z_Py_Pz_Pz+ABZ*I_ERI_F2xz_Py_Pz_Pz;
  abcd[523] = I_ERI_Gx2yz_Py_Pz_Pz+ABZ*I_ERI_Fx2y_Py_Pz_Pz;
  abcd[524] = I_ERI_Gxy2z_Py_Pz_Pz+ABZ*I_ERI_Fxyz_Py_Pz_Pz;
  abcd[525] = I_ERI_Gx3z_Py_Pz_Pz+ABZ*I_ERI_Fx2z_Py_Pz_Pz;
  abcd[526] = I_ERI_G3yz_Py_Pz_Pz+ABZ*I_ERI_F3y_Py_Pz_Pz;
  abcd[527] = I_ERI_G2y2z_Py_Pz_Pz+ABZ*I_ERI_F2yz_Py_Pz_Pz;
  abcd[528] = I_ERI_Gy3z_Py_Pz_Pz+ABZ*I_ERI_Fy2z_Py_Pz_Pz;
  abcd[529] = I_ERI_G4z_Py_Pz_Pz+ABZ*I_ERI_F3z_Py_Pz_Pz;
  abcd[530] = I_ERI_G3xz_Pz_Pz_Pz+ABZ*I_ERI_F3x_Pz_Pz_Pz;
  abcd[531] = I_ERI_G2xyz_Pz_Pz_Pz+ABZ*I_ERI_F2xy_Pz_Pz_Pz;
  abcd[532] = I_ERI_G2x2z_Pz_Pz_Pz+ABZ*I_ERI_F2xz_Pz_Pz_Pz;
  abcd[533] = I_ERI_Gx2yz_Pz_Pz_Pz+ABZ*I_ERI_Fx2y_Pz_Pz_Pz;
  abcd[534] = I_ERI_Gxy2z_Pz_Pz_Pz+ABZ*I_ERI_Fxyz_Pz_Pz_Pz;
  abcd[535] = I_ERI_Gx3z_Pz_Pz_Pz+ABZ*I_ERI_Fx2z_Pz_Pz_Pz;
  abcd[536] = I_ERI_G3yz_Pz_Pz_Pz+ABZ*I_ERI_F3y_Pz_Pz_Pz;
  abcd[537] = I_ERI_G2y2z_Pz_Pz_Pz+ABZ*I_ERI_F2yz_Pz_Pz_Pz;
  abcd[538] = I_ERI_Gy3z_Pz_Pz_Pz+ABZ*I_ERI_Fy2z_Pz_Pz_Pz;
  abcd[539] = I_ERI_G4z_Pz_Pz_Pz+ABZ*I_ERI_F3z_Pz_Pz_Pz;
}
