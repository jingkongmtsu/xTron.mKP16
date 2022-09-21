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
// BRA2
// X
// Y
// Z
// ####

void hgp_os_nai_f_f_d1(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_NAI_H5x_S = 0.0E0;
  Double I_NAI_H4xy_S = 0.0E0;
  Double I_NAI_H4xz_S = 0.0E0;
  Double I_NAI_H3x2y_S = 0.0E0;
  Double I_NAI_H3xyz_S = 0.0E0;
  Double I_NAI_H3x2z_S = 0.0E0;
  Double I_NAI_H2x3y_S = 0.0E0;
  Double I_NAI_H2x2yz_S = 0.0E0;
  Double I_NAI_H2xy2z_S = 0.0E0;
  Double I_NAI_H2x3z_S = 0.0E0;
  Double I_NAI_Hx4y_S = 0.0E0;
  Double I_NAI_Hx3yz_S = 0.0E0;
  Double I_NAI_Hx2y2z_S = 0.0E0;
  Double I_NAI_Hxy3z_S = 0.0E0;
  Double I_NAI_Hx4z_S = 0.0E0;
  Double I_NAI_H5y_S = 0.0E0;
  Double I_NAI_H4yz_S = 0.0E0;
  Double I_NAI_H3y2z_S = 0.0E0;
  Double I_NAI_H2y3z_S = 0.0E0;
  Double I_NAI_Hy4z_S = 0.0E0;
  Double I_NAI_H5z_S = 0.0E0;
  Double I_NAI_G4x_S = 0.0E0;
  Double I_NAI_G3xy_S = 0.0E0;
  Double I_NAI_G3xz_S = 0.0E0;
  Double I_NAI_G2x2y_S = 0.0E0;
  Double I_NAI_G2xyz_S = 0.0E0;
  Double I_NAI_G2x2z_S = 0.0E0;
  Double I_NAI_Gx3y_S = 0.0E0;
  Double I_NAI_Gx2yz_S = 0.0E0;
  Double I_NAI_Gxy2z_S = 0.0E0;
  Double I_NAI_Gx3z_S = 0.0E0;
  Double I_NAI_G4y_S = 0.0E0;
  Double I_NAI_G3yz_S = 0.0E0;
  Double I_NAI_G2y2z_S = 0.0E0;
  Double I_NAI_Gy3z_S = 0.0E0;
  Double I_NAI_G4z_S = 0.0E0;
  Double I_NAI_F3x_S = 0.0E0;
  Double I_NAI_F2xy_S = 0.0E0;
  Double I_NAI_F2xz_S = 0.0E0;
  Double I_NAI_Fx2y_S = 0.0E0;
  Double I_NAI_Fxyz_S = 0.0E0;
  Double I_NAI_Fx2z_S = 0.0E0;
  Double I_NAI_F3y_S = 0.0E0;
  Double I_NAI_F2yz_S = 0.0E0;
  Double I_NAI_Fy2z_S = 0.0E0;
  Double I_NAI_F3z_S = 0.0E0;
  Double I_NAI_K7x_S_a = 0.0E0;
  Double I_NAI_K6xy_S_a = 0.0E0;
  Double I_NAI_K6xz_S_a = 0.0E0;
  Double I_NAI_K5x2y_S_a = 0.0E0;
  Double I_NAI_K5xyz_S_a = 0.0E0;
  Double I_NAI_K5x2z_S_a = 0.0E0;
  Double I_NAI_K4x3y_S_a = 0.0E0;
  Double I_NAI_K4x2yz_S_a = 0.0E0;
  Double I_NAI_K4xy2z_S_a = 0.0E0;
  Double I_NAI_K4x3z_S_a = 0.0E0;
  Double I_NAI_K3x4y_S_a = 0.0E0;
  Double I_NAI_K3x3yz_S_a = 0.0E0;
  Double I_NAI_K3x2y2z_S_a = 0.0E0;
  Double I_NAI_K3xy3z_S_a = 0.0E0;
  Double I_NAI_K3x4z_S_a = 0.0E0;
  Double I_NAI_K2x5y_S_a = 0.0E0;
  Double I_NAI_K2x4yz_S_a = 0.0E0;
  Double I_NAI_K2x3y2z_S_a = 0.0E0;
  Double I_NAI_K2x2y3z_S_a = 0.0E0;
  Double I_NAI_K2xy4z_S_a = 0.0E0;
  Double I_NAI_K2x5z_S_a = 0.0E0;
  Double I_NAI_Kx6y_S_a = 0.0E0;
  Double I_NAI_Kx5yz_S_a = 0.0E0;
  Double I_NAI_Kx4y2z_S_a = 0.0E0;
  Double I_NAI_Kx3y3z_S_a = 0.0E0;
  Double I_NAI_Kx2y4z_S_a = 0.0E0;
  Double I_NAI_Kxy5z_S_a = 0.0E0;
  Double I_NAI_Kx6z_S_a = 0.0E0;
  Double I_NAI_K7y_S_a = 0.0E0;
  Double I_NAI_K6yz_S_a = 0.0E0;
  Double I_NAI_K5y2z_S_a = 0.0E0;
  Double I_NAI_K4y3z_S_a = 0.0E0;
  Double I_NAI_K3y4z_S_a = 0.0E0;
  Double I_NAI_K2y5z_S_a = 0.0E0;
  Double I_NAI_Ky6z_S_a = 0.0E0;
  Double I_NAI_K7z_S_a = 0.0E0;
  Double I_NAI_I6x_S_a = 0.0E0;
  Double I_NAI_I5xy_S_a = 0.0E0;
  Double I_NAI_I5xz_S_a = 0.0E0;
  Double I_NAI_I4x2y_S_a = 0.0E0;
  Double I_NAI_I4xyz_S_a = 0.0E0;
  Double I_NAI_I4x2z_S_a = 0.0E0;
  Double I_NAI_I3x3y_S_a = 0.0E0;
  Double I_NAI_I3x2yz_S_a = 0.0E0;
  Double I_NAI_I3xy2z_S_a = 0.0E0;
  Double I_NAI_I3x3z_S_a = 0.0E0;
  Double I_NAI_I2x4y_S_a = 0.0E0;
  Double I_NAI_I2x3yz_S_a = 0.0E0;
  Double I_NAI_I2x2y2z_S_a = 0.0E0;
  Double I_NAI_I2xy3z_S_a = 0.0E0;
  Double I_NAI_I2x4z_S_a = 0.0E0;
  Double I_NAI_Ix5y_S_a = 0.0E0;
  Double I_NAI_Ix4yz_S_a = 0.0E0;
  Double I_NAI_Ix3y2z_S_a = 0.0E0;
  Double I_NAI_Ix2y3z_S_a = 0.0E0;
  Double I_NAI_Ixy4z_S_a = 0.0E0;
  Double I_NAI_Ix5z_S_a = 0.0E0;
  Double I_NAI_I6y_S_a = 0.0E0;
  Double I_NAI_I5yz_S_a = 0.0E0;
  Double I_NAI_I4y2z_S_a = 0.0E0;
  Double I_NAI_I3y3z_S_a = 0.0E0;
  Double I_NAI_I2y4z_S_a = 0.0E0;
  Double I_NAI_Iy5z_S_a = 0.0E0;
  Double I_NAI_I6z_S_a = 0.0E0;
  Double I_NAI_H5x_S_a = 0.0E0;
  Double I_NAI_H4xy_S_a = 0.0E0;
  Double I_NAI_H4xz_S_a = 0.0E0;
  Double I_NAI_H3x2y_S_a = 0.0E0;
  Double I_NAI_H3xyz_S_a = 0.0E0;
  Double I_NAI_H3x2z_S_a = 0.0E0;
  Double I_NAI_H2x3y_S_a = 0.0E0;
  Double I_NAI_H2x2yz_S_a = 0.0E0;
  Double I_NAI_H2xy2z_S_a = 0.0E0;
  Double I_NAI_H2x3z_S_a = 0.0E0;
  Double I_NAI_Hx4y_S_a = 0.0E0;
  Double I_NAI_Hx3yz_S_a = 0.0E0;
  Double I_NAI_Hx2y2z_S_a = 0.0E0;
  Double I_NAI_Hxy3z_S_a = 0.0E0;
  Double I_NAI_Hx4z_S_a = 0.0E0;
  Double I_NAI_H5y_S_a = 0.0E0;
  Double I_NAI_H4yz_S_a = 0.0E0;
  Double I_NAI_H3y2z_S_a = 0.0E0;
  Double I_NAI_H2y3z_S_a = 0.0E0;
  Double I_NAI_Hy4z_S_a = 0.0E0;
  Double I_NAI_H5z_S_a = 0.0E0;
  Double I_NAI_G4x_S_a = 0.0E0;
  Double I_NAI_G3xy_S_a = 0.0E0;
  Double I_NAI_G3xz_S_a = 0.0E0;
  Double I_NAI_G2x2y_S_a = 0.0E0;
  Double I_NAI_G2xyz_S_a = 0.0E0;
  Double I_NAI_G2x2z_S_a = 0.0E0;
  Double I_NAI_Gx3y_S_a = 0.0E0;
  Double I_NAI_Gx2yz_S_a = 0.0E0;
  Double I_NAI_Gxy2z_S_a = 0.0E0;
  Double I_NAI_Gx3z_S_a = 0.0E0;
  Double I_NAI_G4y_S_a = 0.0E0;
  Double I_NAI_G3yz_S_a = 0.0E0;
  Double I_NAI_G2y2z_S_a = 0.0E0;
  Double I_NAI_Gy3z_S_a = 0.0E0;
  Double I_NAI_G4z_S_a = 0.0E0;
  Double I_NAI_D2x_S = 0.0E0;
  Double I_NAI_Dxy_S = 0.0E0;
  Double I_NAI_Dxz_S = 0.0E0;
  Double I_NAI_D2y_S = 0.0E0;
  Double I_NAI_Dyz_S = 0.0E0;
  Double I_NAI_D2z_S = 0.0E0;
  Double I_NAI_K7x_S_b = 0.0E0;
  Double I_NAI_K6xy_S_b = 0.0E0;
  Double I_NAI_K6xz_S_b = 0.0E0;
  Double I_NAI_K5x2y_S_b = 0.0E0;
  Double I_NAI_K5xyz_S_b = 0.0E0;
  Double I_NAI_K5x2z_S_b = 0.0E0;
  Double I_NAI_K4x3y_S_b = 0.0E0;
  Double I_NAI_K4x2yz_S_b = 0.0E0;
  Double I_NAI_K4xy2z_S_b = 0.0E0;
  Double I_NAI_K4x3z_S_b = 0.0E0;
  Double I_NAI_K3x4y_S_b = 0.0E0;
  Double I_NAI_K3x3yz_S_b = 0.0E0;
  Double I_NAI_K3x2y2z_S_b = 0.0E0;
  Double I_NAI_K3xy3z_S_b = 0.0E0;
  Double I_NAI_K3x4z_S_b = 0.0E0;
  Double I_NAI_K2x5y_S_b = 0.0E0;
  Double I_NAI_K2x4yz_S_b = 0.0E0;
  Double I_NAI_K2x3y2z_S_b = 0.0E0;
  Double I_NAI_K2x2y3z_S_b = 0.0E0;
  Double I_NAI_K2xy4z_S_b = 0.0E0;
  Double I_NAI_K2x5z_S_b = 0.0E0;
  Double I_NAI_Kx6y_S_b = 0.0E0;
  Double I_NAI_Kx5yz_S_b = 0.0E0;
  Double I_NAI_Kx4y2z_S_b = 0.0E0;
  Double I_NAI_Kx3y3z_S_b = 0.0E0;
  Double I_NAI_Kx2y4z_S_b = 0.0E0;
  Double I_NAI_Kxy5z_S_b = 0.0E0;
  Double I_NAI_Kx6z_S_b = 0.0E0;
  Double I_NAI_K7y_S_b = 0.0E0;
  Double I_NAI_K6yz_S_b = 0.0E0;
  Double I_NAI_K5y2z_S_b = 0.0E0;
  Double I_NAI_K4y3z_S_b = 0.0E0;
  Double I_NAI_K3y4z_S_b = 0.0E0;
  Double I_NAI_K2y5z_S_b = 0.0E0;
  Double I_NAI_Ky6z_S_b = 0.0E0;
  Double I_NAI_K7z_S_b = 0.0E0;
  Double I_NAI_I6x_S_b = 0.0E0;
  Double I_NAI_I5xy_S_b = 0.0E0;
  Double I_NAI_I5xz_S_b = 0.0E0;
  Double I_NAI_I4x2y_S_b = 0.0E0;
  Double I_NAI_I4xyz_S_b = 0.0E0;
  Double I_NAI_I4x2z_S_b = 0.0E0;
  Double I_NAI_I3x3y_S_b = 0.0E0;
  Double I_NAI_I3x2yz_S_b = 0.0E0;
  Double I_NAI_I3xy2z_S_b = 0.0E0;
  Double I_NAI_I3x3z_S_b = 0.0E0;
  Double I_NAI_I2x4y_S_b = 0.0E0;
  Double I_NAI_I2x3yz_S_b = 0.0E0;
  Double I_NAI_I2x2y2z_S_b = 0.0E0;
  Double I_NAI_I2xy3z_S_b = 0.0E0;
  Double I_NAI_I2x4z_S_b = 0.0E0;
  Double I_NAI_Ix5y_S_b = 0.0E0;
  Double I_NAI_Ix4yz_S_b = 0.0E0;
  Double I_NAI_Ix3y2z_S_b = 0.0E0;
  Double I_NAI_Ix2y3z_S_b = 0.0E0;
  Double I_NAI_Ixy4z_S_b = 0.0E0;
  Double I_NAI_Ix5z_S_b = 0.0E0;
  Double I_NAI_I6y_S_b = 0.0E0;
  Double I_NAI_I5yz_S_b = 0.0E0;
  Double I_NAI_I4y2z_S_b = 0.0E0;
  Double I_NAI_I3y3z_S_b = 0.0E0;
  Double I_NAI_I2y4z_S_b = 0.0E0;
  Double I_NAI_Iy5z_S_b = 0.0E0;
  Double I_NAI_I6z_S_b = 0.0E0;
  Double I_NAI_H5x_S_b = 0.0E0;
  Double I_NAI_H4xy_S_b = 0.0E0;
  Double I_NAI_H4xz_S_b = 0.0E0;
  Double I_NAI_H3x2y_S_b = 0.0E0;
  Double I_NAI_H3xyz_S_b = 0.0E0;
  Double I_NAI_H3x2z_S_b = 0.0E0;
  Double I_NAI_H2x3y_S_b = 0.0E0;
  Double I_NAI_H2x2yz_S_b = 0.0E0;
  Double I_NAI_H2xy2z_S_b = 0.0E0;
  Double I_NAI_H2x3z_S_b = 0.0E0;
  Double I_NAI_Hx4y_S_b = 0.0E0;
  Double I_NAI_Hx3yz_S_b = 0.0E0;
  Double I_NAI_Hx2y2z_S_b = 0.0E0;
  Double I_NAI_Hxy3z_S_b = 0.0E0;
  Double I_NAI_Hx4z_S_b = 0.0E0;
  Double I_NAI_H5y_S_b = 0.0E0;
  Double I_NAI_H4yz_S_b = 0.0E0;
  Double I_NAI_H3y2z_S_b = 0.0E0;
  Double I_NAI_H2y3z_S_b = 0.0E0;
  Double I_NAI_Hy4z_S_b = 0.0E0;
  Double I_NAI_H5z_S_b = 0.0E0;
  Double I_NAI_G4x_S_b = 0.0E0;
  Double I_NAI_G3xy_S_b = 0.0E0;
  Double I_NAI_G3xz_S_b = 0.0E0;
  Double I_NAI_G2x2y_S_b = 0.0E0;
  Double I_NAI_G2xyz_S_b = 0.0E0;
  Double I_NAI_G2x2z_S_b = 0.0E0;
  Double I_NAI_Gx3y_S_b = 0.0E0;
  Double I_NAI_Gx2yz_S_b = 0.0E0;
  Double I_NAI_Gxy2z_S_b = 0.0E0;
  Double I_NAI_Gx3z_S_b = 0.0E0;
  Double I_NAI_G4y_S_b = 0.0E0;
  Double I_NAI_G3yz_S_b = 0.0E0;
  Double I_NAI_G2y2z_S_b = 0.0E0;
  Double I_NAI_Gy3z_S_b = 0.0E0;
  Double I_NAI_G4z_S_b = 0.0E0;
  Double I_NAI_F3x_S_b = 0.0E0;
  Double I_NAI_F2xy_S_b = 0.0E0;
  Double I_NAI_F2xz_S_b = 0.0E0;
  Double I_NAI_Fx2y_S_b = 0.0E0;
  Double I_NAI_Fxyz_S_b = 0.0E0;
  Double I_NAI_Fx2z_S_b = 0.0E0;
  Double I_NAI_F3y_S_b = 0.0E0;
  Double I_NAI_F2yz_S_b = 0.0E0;
  Double I_NAI_Fy2z_S_b = 0.0E0;
  Double I_NAI_F3z_S_b = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double onedz = iexp[ip2];
    Double rho   = 1.0E0/onedz;
    Double zeta  = rho;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double sqrho = sqrt(rho);
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
      Double PNX   = PX - N[iAtom*3  ];
      Double PNY   = PY - N[iAtom*3+1];
      Double PNZ   = PZ - N[iAtom*3+2];
      Double PN2   = PNX*PNX+PNY*PNY+PNZ*PNZ;
      Double charge= Z[iAtom];
      Double u     = rho*PN2;
      Double squ   = sqrt(u);
      Double prefactor = -ic2*charge*fbra;

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

      Double I_NAI_S_S_vrr  = 0.0E0;
      Double I_NAI_S_S_M1_vrr  = 0.0E0;
      Double I_NAI_S_S_M2_vrr  = 0.0E0;
      Double I_NAI_S_S_M3_vrr  = 0.0E0;
      Double I_NAI_S_S_M4_vrr  = 0.0E0;
      Double I_NAI_S_S_M5_vrr  = 0.0E0;
      Double I_NAI_S_S_M6_vrr  = 0.0E0;
      Double I_NAI_S_S_M7_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_NAI_S_S_vrr_d  = 0.0E0;
      double I_NAI_S_S_M1_vrr_d  = 0.0E0;
      double I_NAI_S_S_M2_vrr_d  = 0.0E0;
      double I_NAI_S_S_M3_vrr_d  = 0.0E0;
      double I_NAI_S_S_M4_vrr_d  = 0.0E0;
      double I_NAI_S_S_M5_vrr_d  = 0.0E0;
      double I_NAI_S_S_M6_vrr_d  = 0.0E0;
      double I_NAI_S_S_M7_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER49;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER47*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER45*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER43*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER41*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER39*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER37*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER35*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER33*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER31*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER29*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER27*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER25*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER23*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER21*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER19*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = 1.0E0+u2*ONEOVER17*I_NAI_S_S_M7_vrr;
        I_NAI_S_S_M7_vrr = ONEOVER15*I_NAI_S_S_M7_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M7_vrr  = f*I_NAI_S_S_M7_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_NAI_S_S_M6_vrr  = ONEOVER13*(u2*I_NAI_S_S_M7_vrr+f);
        I_NAI_S_S_M5_vrr  = ONEOVER11*(u2*I_NAI_S_S_M6_vrr+f);
        I_NAI_S_S_M4_vrr  = ONEOVER9*(u2*I_NAI_S_S_M5_vrr+f);
        I_NAI_S_S_M3_vrr  = ONEOVER7*(u2*I_NAI_S_S_M4_vrr+f);
        I_NAI_S_S_M2_vrr  = ONEOVER5*(u2*I_NAI_S_S_M3_vrr+f);
        I_NAI_S_S_M1_vrr  = ONEOVER3*(u2*I_NAI_S_S_M2_vrr+f);
        I_NAI_S_S_vrr  = ONEOVER1*(u2*I_NAI_S_S_M1_vrr+f);

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
          I_NAI_S_S_vrr_d = fac_d*sqrho_d*TWOOVERSQRTPI;
        }else{
          I_NAI_S_S_vrr_d = (fac_d*sqrho_d/squ_d)*erf(squ_d);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        double oneO2u  = 0.5E0/u_d;
        double eu      = exp(-u_d);
        double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;
        I_NAI_S_S_M1_vrr_d = oneO2u*(1.0E0*I_NAI_S_S_vrr_d-f);
        I_NAI_S_S_M2_vrr_d = oneO2u*(3.0E0*I_NAI_S_S_M1_vrr_d-f);
        I_NAI_S_S_M3_vrr_d = oneO2u*(5.0E0*I_NAI_S_S_M2_vrr_d-f);
        I_NAI_S_S_M4_vrr_d = oneO2u*(7.0E0*I_NAI_S_S_M3_vrr_d-f);
        I_NAI_S_S_M5_vrr_d = oneO2u*(9.0E0*I_NAI_S_S_M4_vrr_d-f);
        I_NAI_S_S_M6_vrr_d = oneO2u*(11.0E0*I_NAI_S_S_M5_vrr_d-f);
        I_NAI_S_S_M7_vrr_d = oneO2u*(13.0E0*I_NAI_S_S_M6_vrr_d-f);

        // write the double result back to the float var
        I_NAI_S_S_vrr = static_cast<Double>(I_NAI_S_S_vrr_d);
        I_NAI_S_S_M1_vrr = static_cast<Double>(I_NAI_S_S_M1_vrr_d);
        I_NAI_S_S_M2_vrr = static_cast<Double>(I_NAI_S_S_M2_vrr_d);
        I_NAI_S_S_M3_vrr = static_cast<Double>(I_NAI_S_S_M3_vrr_d);
        I_NAI_S_S_M4_vrr = static_cast<Double>(I_NAI_S_S_M4_vrr_d);
        I_NAI_S_S_M5_vrr = static_cast<Double>(I_NAI_S_S_M5_vrr_d);
        I_NAI_S_S_M6_vrr = static_cast<Double>(I_NAI_S_S_M6_vrr_d);
        I_NAI_S_S_M7_vrr = static_cast<Double>(I_NAI_S_S_M7_vrr_d);

#else

        // use erf function to get (SS|SS)^{0}
        if (fabs(u)<THRESHOLD_MATH) {
          I_NAI_S_S_vrr = prefactor*sqrho*TWOOVERSQRTPI;
        }else{
          I_NAI_S_S_vrr = (prefactor*sqrho/squ)*erf(squ);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        Double oneO2u = 0.5E0/u;
        Double eu     = exp(-u);
        Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_NAI_S_S_M1_vrr = oneO2u*(1.0E0*I_NAI_S_S_vrr-f);
        I_NAI_S_S_M2_vrr = oneO2u*(3.0E0*I_NAI_S_S_M1_vrr-f);
        I_NAI_S_S_M3_vrr = oneO2u*(5.0E0*I_NAI_S_S_M2_vrr-f);
        I_NAI_S_S_M4_vrr = oneO2u*(7.0E0*I_NAI_S_S_M3_vrr-f);
        I_NAI_S_S_M5_vrr = oneO2u*(9.0E0*I_NAI_S_S_M4_vrr-f);
        I_NAI_S_S_M6_vrr = oneO2u*(11.0E0*I_NAI_S_S_M5_vrr-f);
        I_NAI_S_S_M7_vrr = oneO2u*(13.0E0*I_NAI_S_S_M6_vrr-f);

#endif

      }


      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M6
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M6
       * RHS shell quartet name: SQ_NAI_S_S_M7
       ************************************************************/
      Double I_NAI_Px_S_M6_vrr = PAX*I_NAI_S_S_M6_vrr-PNX*I_NAI_S_S_M7_vrr;
      Double I_NAI_Py_S_M6_vrr = PAY*I_NAI_S_S_M6_vrr-PNY*I_NAI_S_S_M7_vrr;
      Double I_NAI_Pz_S_M6_vrr = PAZ*I_NAI_S_S_M6_vrr-PNZ*I_NAI_S_S_M7_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M5
       * RHS shell quartet name: SQ_NAI_S_S_M6
       ************************************************************/
      Double I_NAI_Px_S_M5_vrr = PAX*I_NAI_S_S_M5_vrr-PNX*I_NAI_S_S_M6_vrr;
      Double I_NAI_Py_S_M5_vrr = PAY*I_NAI_S_S_M5_vrr-PNY*I_NAI_S_S_M6_vrr;
      Double I_NAI_Pz_S_M5_vrr = PAZ*I_NAI_S_S_M5_vrr-PNZ*I_NAI_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M5
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M5
       * RHS shell quartet name: SQ_NAI_P_S_M6
       * RHS shell quartet name: SQ_NAI_S_S_M5
       * RHS shell quartet name: SQ_NAI_S_S_M6
       ************************************************************/
      Double I_NAI_D2x_S_M5_vrr = PAX*I_NAI_Px_S_M5_vrr-PNX*I_NAI_Px_S_M6_vrr+oned2z*I_NAI_S_S_M5_vrr-oned2z*I_NAI_S_S_M6_vrr;
      Double I_NAI_D2y_S_M5_vrr = PAY*I_NAI_Py_S_M5_vrr-PNY*I_NAI_Py_S_M6_vrr+oned2z*I_NAI_S_S_M5_vrr-oned2z*I_NAI_S_S_M6_vrr;
      Double I_NAI_D2z_S_M5_vrr = PAZ*I_NAI_Pz_S_M5_vrr-PNZ*I_NAI_Pz_S_M6_vrr+oned2z*I_NAI_S_S_M5_vrr-oned2z*I_NAI_S_S_M6_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M4
       * RHS shell quartet name: SQ_NAI_S_S_M5
       ************************************************************/
      Double I_NAI_Px_S_M4_vrr = PAX*I_NAI_S_S_M4_vrr-PNX*I_NAI_S_S_M5_vrr;
      Double I_NAI_Py_S_M4_vrr = PAY*I_NAI_S_S_M4_vrr-PNY*I_NAI_S_S_M5_vrr;
      Double I_NAI_Pz_S_M4_vrr = PAZ*I_NAI_S_S_M4_vrr-PNZ*I_NAI_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M4
       * RHS shell quartet name: SQ_NAI_P_S_M5
       * RHS shell quartet name: SQ_NAI_S_S_M4
       * RHS shell quartet name: SQ_NAI_S_S_M5
       ************************************************************/
      Double I_NAI_D2x_S_M4_vrr = PAX*I_NAI_Px_S_M4_vrr-PNX*I_NAI_Px_S_M5_vrr+oned2z*I_NAI_S_S_M4_vrr-oned2z*I_NAI_S_S_M5_vrr;
      Double I_NAI_D2y_S_M4_vrr = PAY*I_NAI_Py_S_M4_vrr-PNY*I_NAI_Py_S_M5_vrr+oned2z*I_NAI_S_S_M4_vrr-oned2z*I_NAI_S_S_M5_vrr;
      Double I_NAI_D2z_S_M4_vrr = PAZ*I_NAI_Pz_S_M4_vrr-PNZ*I_NAI_Pz_S_M5_vrr+oned2z*I_NAI_S_S_M4_vrr-oned2z*I_NAI_S_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M4
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M4
       * RHS shell quartet name: SQ_NAI_D_S_M5
       * RHS shell quartet name: SQ_NAI_P_S_M4
       * RHS shell quartet name: SQ_NAI_P_S_M5
       ************************************************************/
      Double I_NAI_F3x_S_M4_vrr = PAX*I_NAI_D2x_S_M4_vrr-PNX*I_NAI_D2x_S_M5_vrr+2*oned2z*I_NAI_Px_S_M4_vrr-2*oned2z*I_NAI_Px_S_M5_vrr;
      Double I_NAI_F3y_S_M4_vrr = PAY*I_NAI_D2y_S_M4_vrr-PNY*I_NAI_D2y_S_M5_vrr+2*oned2z*I_NAI_Py_S_M4_vrr-2*oned2z*I_NAI_Py_S_M5_vrr;
      Double I_NAI_F3z_S_M4_vrr = PAZ*I_NAI_D2z_S_M4_vrr-PNZ*I_NAI_D2z_S_M5_vrr+2*oned2z*I_NAI_Pz_S_M4_vrr-2*oned2z*I_NAI_Pz_S_M5_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M3
       * RHS shell quartet name: SQ_NAI_S_S_M4
       ************************************************************/
      Double I_NAI_Px_S_M3_vrr = PAX*I_NAI_S_S_M3_vrr-PNX*I_NAI_S_S_M4_vrr;
      Double I_NAI_Py_S_M3_vrr = PAY*I_NAI_S_S_M3_vrr-PNY*I_NAI_S_S_M4_vrr;
      Double I_NAI_Pz_S_M3_vrr = PAZ*I_NAI_S_S_M3_vrr-PNZ*I_NAI_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M3
       * RHS shell quartet name: SQ_NAI_P_S_M4
       * RHS shell quartet name: SQ_NAI_S_S_M3
       * RHS shell quartet name: SQ_NAI_S_S_M4
       ************************************************************/
      Double I_NAI_D2x_S_M3_vrr = PAX*I_NAI_Px_S_M3_vrr-PNX*I_NAI_Px_S_M4_vrr+oned2z*I_NAI_S_S_M3_vrr-oned2z*I_NAI_S_S_M4_vrr;
      Double I_NAI_D2y_S_M3_vrr = PAY*I_NAI_Py_S_M3_vrr-PNY*I_NAI_Py_S_M4_vrr+oned2z*I_NAI_S_S_M3_vrr-oned2z*I_NAI_S_S_M4_vrr;
      Double I_NAI_D2z_S_M3_vrr = PAZ*I_NAI_Pz_S_M3_vrr-PNZ*I_NAI_Pz_S_M4_vrr+oned2z*I_NAI_S_S_M3_vrr-oned2z*I_NAI_S_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 6 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M3
       * RHS shell quartet name: SQ_NAI_D_S_M4
       * RHS shell quartet name: SQ_NAI_P_S_M3
       * RHS shell quartet name: SQ_NAI_P_S_M4
       ************************************************************/
      Double I_NAI_F3x_S_M3_vrr = PAX*I_NAI_D2x_S_M3_vrr-PNX*I_NAI_D2x_S_M4_vrr+2*oned2z*I_NAI_Px_S_M3_vrr-2*oned2z*I_NAI_Px_S_M4_vrr;
      Double I_NAI_F2xy_S_M3_vrr = PAY*I_NAI_D2x_S_M3_vrr-PNY*I_NAI_D2x_S_M4_vrr;
      Double I_NAI_F3y_S_M3_vrr = PAY*I_NAI_D2y_S_M3_vrr-PNY*I_NAI_D2y_S_M4_vrr+2*oned2z*I_NAI_Py_S_M3_vrr-2*oned2z*I_NAI_Py_S_M4_vrr;
      Double I_NAI_F3z_S_M3_vrr = PAZ*I_NAI_D2z_S_M3_vrr-PNZ*I_NAI_D2z_S_M4_vrr+2*oned2z*I_NAI_Pz_S_M3_vrr-2*oned2z*I_NAI_Pz_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_M3
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M3
       * RHS shell quartet name: SQ_NAI_F_S_M4
       * RHS shell quartet name: SQ_NAI_D_S_M3
       * RHS shell quartet name: SQ_NAI_D_S_M4
       ************************************************************/
      Double I_NAI_G4x_S_M3_vrr = PAX*I_NAI_F3x_S_M3_vrr-PNX*I_NAI_F3x_S_M4_vrr+3*oned2z*I_NAI_D2x_S_M3_vrr-3*oned2z*I_NAI_D2x_S_M4_vrr;
      Double I_NAI_G3xy_S_M3_vrr = PAY*I_NAI_F3x_S_M3_vrr-PNY*I_NAI_F3x_S_M4_vrr;
      Double I_NAI_G3xz_S_M3_vrr = PAZ*I_NAI_F3x_S_M3_vrr-PNZ*I_NAI_F3x_S_M4_vrr;
      Double I_NAI_Gx3y_S_M3_vrr = PAX*I_NAI_F3y_S_M3_vrr-PNX*I_NAI_F3y_S_M4_vrr;
      Double I_NAI_Gx3z_S_M3_vrr = PAX*I_NAI_F3z_S_M3_vrr-PNX*I_NAI_F3z_S_M4_vrr;
      Double I_NAI_G4y_S_M3_vrr = PAY*I_NAI_F3y_S_M3_vrr-PNY*I_NAI_F3y_S_M4_vrr+3*oned2z*I_NAI_D2y_S_M3_vrr-3*oned2z*I_NAI_D2y_S_M4_vrr;
      Double I_NAI_G3yz_S_M3_vrr = PAZ*I_NAI_F3y_S_M3_vrr-PNZ*I_NAI_F3y_S_M4_vrr;
      Double I_NAI_G4z_S_M3_vrr = PAZ*I_NAI_F3z_S_M3_vrr-PNZ*I_NAI_F3z_S_M4_vrr+3*oned2z*I_NAI_D2z_S_M3_vrr-3*oned2z*I_NAI_D2z_S_M4_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M3
       ************************************************************/
      Double I_NAI_Px_S_M2_vrr = PAX*I_NAI_S_S_M2_vrr-PNX*I_NAI_S_S_M3_vrr;
      Double I_NAI_Py_S_M2_vrr = PAY*I_NAI_S_S_M2_vrr-PNY*I_NAI_S_S_M3_vrr;
      Double I_NAI_Pz_S_M2_vrr = PAZ*I_NAI_S_S_M2_vrr-PNZ*I_NAI_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M3
       * RHS shell quartet name: SQ_NAI_S_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M3
       ************************************************************/
      Double I_NAI_D2x_S_M2_vrr = PAX*I_NAI_Px_S_M2_vrr-PNX*I_NAI_Px_S_M3_vrr+oned2z*I_NAI_S_S_M2_vrr-oned2z*I_NAI_S_S_M3_vrr;
      Double I_NAI_D2y_S_M2_vrr = PAY*I_NAI_Py_S_M2_vrr-PNY*I_NAI_Py_S_M3_vrr+oned2z*I_NAI_S_S_M2_vrr-oned2z*I_NAI_S_S_M3_vrr;
      Double I_NAI_D2z_S_M2_vrr = PAZ*I_NAI_Pz_S_M2_vrr-PNZ*I_NAI_Pz_S_M3_vrr+oned2z*I_NAI_S_S_M2_vrr-oned2z*I_NAI_S_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 4 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M2
       * RHS shell quartet name: SQ_NAI_D_S_M3
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M3
       ************************************************************/
      Double I_NAI_F3x_S_M2_vrr = PAX*I_NAI_D2x_S_M2_vrr-PNX*I_NAI_D2x_S_M3_vrr+2*oned2z*I_NAI_Px_S_M2_vrr-2*oned2z*I_NAI_Px_S_M3_vrr;
      Double I_NAI_F2xy_S_M2_vrr = PAY*I_NAI_D2x_S_M2_vrr-PNY*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_F2xz_S_M2_vrr = PAZ*I_NAI_D2x_S_M2_vrr-PNZ*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_F3y_S_M2_vrr = PAY*I_NAI_D2y_S_M2_vrr-PNY*I_NAI_D2y_S_M3_vrr+2*oned2z*I_NAI_Py_S_M2_vrr-2*oned2z*I_NAI_Py_S_M3_vrr;
      Double I_NAI_F2yz_S_M2_vrr = PAZ*I_NAI_D2y_S_M2_vrr-PNZ*I_NAI_D2y_S_M3_vrr;
      Double I_NAI_F3z_S_M2_vrr = PAZ*I_NAI_D2z_S_M2_vrr-PNZ*I_NAI_D2z_S_M3_vrr+2*oned2z*I_NAI_Pz_S_M2_vrr-2*oned2z*I_NAI_Pz_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 5 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M2
       * RHS shell quartet name: SQ_NAI_F_S_M3
       * RHS shell quartet name: SQ_NAI_D_S_M2
       * RHS shell quartet name: SQ_NAI_D_S_M3
       ************************************************************/
      Double I_NAI_G4x_S_M2_vrr = PAX*I_NAI_F3x_S_M2_vrr-PNX*I_NAI_F3x_S_M3_vrr+3*oned2z*I_NAI_D2x_S_M2_vrr-3*oned2z*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_G3xy_S_M2_vrr = PAY*I_NAI_F3x_S_M2_vrr-PNY*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_G3xz_S_M2_vrr = PAZ*I_NAI_F3x_S_M2_vrr-PNZ*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_G2x2y_S_M2_vrr = PAY*I_NAI_F2xy_S_M2_vrr-PNY*I_NAI_F2xy_S_M3_vrr+oned2z*I_NAI_D2x_S_M2_vrr-oned2z*I_NAI_D2x_S_M3_vrr;
      Double I_NAI_Gx3y_S_M2_vrr = PAX*I_NAI_F3y_S_M2_vrr-PNX*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_Gx3z_S_M2_vrr = PAX*I_NAI_F3z_S_M2_vrr-PNX*I_NAI_F3z_S_M3_vrr;
      Double I_NAI_G4y_S_M2_vrr = PAY*I_NAI_F3y_S_M2_vrr-PNY*I_NAI_F3y_S_M3_vrr+3*oned2z*I_NAI_D2y_S_M2_vrr-3*oned2z*I_NAI_D2y_S_M3_vrr;
      Double I_NAI_G3yz_S_M2_vrr = PAZ*I_NAI_F3y_S_M2_vrr-PNZ*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_Gy3z_S_M2_vrr = PAY*I_NAI_F3z_S_M2_vrr-PNY*I_NAI_F3z_S_M3_vrr;
      Double I_NAI_G4z_S_M2_vrr = PAZ*I_NAI_F3z_S_M2_vrr-PNZ*I_NAI_F3z_S_M3_vrr+3*oned2z*I_NAI_D2z_S_M2_vrr-3*oned2z*I_NAI_D2z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_M2
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S_M2
       * RHS shell quartet name: SQ_NAI_G_S_M3
       * RHS shell quartet name: SQ_NAI_F_S_M2
       * RHS shell quartet name: SQ_NAI_F_S_M3
       ************************************************************/
      Double I_NAI_H5x_S_M2_vrr = PAX*I_NAI_G4x_S_M2_vrr-PNX*I_NAI_G4x_S_M3_vrr+4*oned2z*I_NAI_F3x_S_M2_vrr-4*oned2z*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_H4xy_S_M2_vrr = PAY*I_NAI_G4x_S_M2_vrr-PNY*I_NAI_G4x_S_M3_vrr;
      Double I_NAI_H4xz_S_M2_vrr = PAZ*I_NAI_G4x_S_M2_vrr-PNZ*I_NAI_G4x_S_M3_vrr;
      Double I_NAI_H3x2y_S_M2_vrr = PAY*I_NAI_G3xy_S_M2_vrr-PNY*I_NAI_G3xy_S_M3_vrr+oned2z*I_NAI_F3x_S_M2_vrr-oned2z*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_H3x2z_S_M2_vrr = PAZ*I_NAI_G3xz_S_M2_vrr-PNZ*I_NAI_G3xz_S_M3_vrr+oned2z*I_NAI_F3x_S_M2_vrr-oned2z*I_NAI_F3x_S_M3_vrr;
      Double I_NAI_H2x3y_S_M2_vrr = PAX*I_NAI_Gx3y_S_M2_vrr-PNX*I_NAI_Gx3y_S_M3_vrr+oned2z*I_NAI_F3y_S_M2_vrr-oned2z*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_H2x3z_S_M2_vrr = PAX*I_NAI_Gx3z_S_M2_vrr-PNX*I_NAI_Gx3z_S_M3_vrr+oned2z*I_NAI_F3z_S_M2_vrr-oned2z*I_NAI_F3z_S_M3_vrr;
      Double I_NAI_Hx4y_S_M2_vrr = PAX*I_NAI_G4y_S_M2_vrr-PNX*I_NAI_G4y_S_M3_vrr;
      Double I_NAI_Hx4z_S_M2_vrr = PAX*I_NAI_G4z_S_M2_vrr-PNX*I_NAI_G4z_S_M3_vrr;
      Double I_NAI_H5y_S_M2_vrr = PAY*I_NAI_G4y_S_M2_vrr-PNY*I_NAI_G4y_S_M3_vrr+4*oned2z*I_NAI_F3y_S_M2_vrr-4*oned2z*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_H4yz_S_M2_vrr = PAZ*I_NAI_G4y_S_M2_vrr-PNZ*I_NAI_G4y_S_M3_vrr;
      Double I_NAI_H3y2z_S_M2_vrr = PAZ*I_NAI_G3yz_S_M2_vrr-PNZ*I_NAI_G3yz_S_M3_vrr+oned2z*I_NAI_F3y_S_M2_vrr-oned2z*I_NAI_F3y_S_M3_vrr;
      Double I_NAI_Hy4z_S_M2_vrr = PAY*I_NAI_G4z_S_M2_vrr-PNY*I_NAI_G4z_S_M3_vrr;
      Double I_NAI_H5z_S_M2_vrr = PAZ*I_NAI_G4z_S_M2_vrr-PNZ*I_NAI_G4z_S_M3_vrr+4*oned2z*I_NAI_F3z_S_M2_vrr-4*oned2z*I_NAI_F3z_S_M3_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_Px_S_M1_vrr = PAX*I_NAI_S_S_M1_vrr-PNX*I_NAI_S_S_M2_vrr;
      Double I_NAI_Py_S_M1_vrr = PAY*I_NAI_S_S_M1_vrr-PNY*I_NAI_S_S_M2_vrr;
      Double I_NAI_Pz_S_M1_vrr = PAZ*I_NAI_S_S_M1_vrr-PNZ*I_NAI_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_P_S_M2
       * RHS shell quartet name: SQ_NAI_S_S_M1
       * RHS shell quartet name: SQ_NAI_S_S_M2
       ************************************************************/
      Double I_NAI_D2x_S_M1_vrr = PAX*I_NAI_Px_S_M1_vrr-PNX*I_NAI_Px_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_Dxy_S_M1_vrr = PAY*I_NAI_Px_S_M1_vrr-PNY*I_NAI_Px_S_M2_vrr;
      Double I_NAI_D2y_S_M1_vrr = PAY*I_NAI_Py_S_M1_vrr-PNY*I_NAI_Py_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;
      Double I_NAI_D2z_S_M1_vrr = PAZ*I_NAI_Pz_S_M1_vrr-PNZ*I_NAI_Pz_S_M2_vrr+oned2z*I_NAI_S_S_M1_vrr-oned2z*I_NAI_S_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 2 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_D_S_M2
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_P_S_M2
       ************************************************************/
      Double I_NAI_F3x_S_M1_vrr = PAX*I_NAI_D2x_S_M1_vrr-PNX*I_NAI_D2x_S_M2_vrr+2*oned2z*I_NAI_Px_S_M1_vrr-2*oned2z*I_NAI_Px_S_M2_vrr;
      Double I_NAI_F2xy_S_M1_vrr = PAY*I_NAI_D2x_S_M1_vrr-PNY*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_F2xz_S_M1_vrr = PAZ*I_NAI_D2x_S_M1_vrr-PNZ*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_Fx2y_S_M1_vrr = PAX*I_NAI_D2y_S_M1_vrr-PNX*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_Fx2z_S_M1_vrr = PAX*I_NAI_D2z_S_M1_vrr-PNX*I_NAI_D2z_S_M2_vrr;
      Double I_NAI_F3y_S_M1_vrr = PAY*I_NAI_D2y_S_M1_vrr-PNY*I_NAI_D2y_S_M2_vrr+2*oned2z*I_NAI_Py_S_M1_vrr-2*oned2z*I_NAI_Py_S_M2_vrr;
      Double I_NAI_F2yz_S_M1_vrr = PAZ*I_NAI_D2y_S_M1_vrr-PNZ*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_F3z_S_M1_vrr = PAZ*I_NAI_D2z_S_M1_vrr-PNZ*I_NAI_D2z_S_M2_vrr+2*oned2z*I_NAI_Pz_S_M1_vrr-2*oned2z*I_NAI_Pz_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 3 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S_M1
       * RHS shell quartet name: SQ_NAI_F_S_M2
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_D_S_M2
       ************************************************************/
      Double I_NAI_G4x_S_M1_vrr = PAX*I_NAI_F3x_S_M1_vrr-PNX*I_NAI_F3x_S_M2_vrr+3*oned2z*I_NAI_D2x_S_M1_vrr-3*oned2z*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_G3xy_S_M1_vrr = PAY*I_NAI_F3x_S_M1_vrr-PNY*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_G3xz_S_M1_vrr = PAZ*I_NAI_F3x_S_M1_vrr-PNZ*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_G2x2y_S_M1_vrr = PAY*I_NAI_F2xy_S_M1_vrr-PNY*I_NAI_F2xy_S_M2_vrr+oned2z*I_NAI_D2x_S_M1_vrr-oned2z*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_G2x2z_S_M1_vrr = PAZ*I_NAI_F2xz_S_M1_vrr-PNZ*I_NAI_F2xz_S_M2_vrr+oned2z*I_NAI_D2x_S_M1_vrr-oned2z*I_NAI_D2x_S_M2_vrr;
      Double I_NAI_Gx3y_S_M1_vrr = PAX*I_NAI_F3y_S_M1_vrr-PNX*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_Gx3z_S_M1_vrr = PAX*I_NAI_F3z_S_M1_vrr-PNX*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_G4y_S_M1_vrr = PAY*I_NAI_F3y_S_M1_vrr-PNY*I_NAI_F3y_S_M2_vrr+3*oned2z*I_NAI_D2y_S_M1_vrr-3*oned2z*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_G3yz_S_M1_vrr = PAZ*I_NAI_F3y_S_M1_vrr-PNZ*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_G2y2z_S_M1_vrr = PAZ*I_NAI_F2yz_S_M1_vrr-PNZ*I_NAI_F2yz_S_M2_vrr+oned2z*I_NAI_D2y_S_M1_vrr-oned2z*I_NAI_D2y_S_M2_vrr;
      Double I_NAI_Gy3z_S_M1_vrr = PAY*I_NAI_F3z_S_M1_vrr-PNY*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_G4z_S_M1_vrr = PAZ*I_NAI_F3z_S_M1_vrr-PNZ*I_NAI_F3z_S_M2_vrr+3*oned2z*I_NAI_D2z_S_M1_vrr-3*oned2z*I_NAI_D2z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 5 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S_M1
       * RHS shell quartet name: SQ_NAI_G_S_M2
       * RHS shell quartet name: SQ_NAI_F_S_M1
       * RHS shell quartet name: SQ_NAI_F_S_M2
       ************************************************************/
      Double I_NAI_H5x_S_M1_vrr = PAX*I_NAI_G4x_S_M1_vrr-PNX*I_NAI_G4x_S_M2_vrr+4*oned2z*I_NAI_F3x_S_M1_vrr-4*oned2z*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_H4xy_S_M1_vrr = PAY*I_NAI_G4x_S_M1_vrr-PNY*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_H4xz_S_M1_vrr = PAZ*I_NAI_G4x_S_M1_vrr-PNZ*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_H3x2y_S_M1_vrr = PAY*I_NAI_G3xy_S_M1_vrr-PNY*I_NAI_G3xy_S_M2_vrr+oned2z*I_NAI_F3x_S_M1_vrr-oned2z*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_H3x2z_S_M1_vrr = PAZ*I_NAI_G3xz_S_M1_vrr-PNZ*I_NAI_G3xz_S_M2_vrr+oned2z*I_NAI_F3x_S_M1_vrr-oned2z*I_NAI_F3x_S_M2_vrr;
      Double I_NAI_H2x3y_S_M1_vrr = PAX*I_NAI_Gx3y_S_M1_vrr-PNX*I_NAI_Gx3y_S_M2_vrr+oned2z*I_NAI_F3y_S_M1_vrr-oned2z*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_H2x2yz_S_M1_vrr = PAZ*I_NAI_G2x2y_S_M1_vrr-PNZ*I_NAI_G2x2y_S_M2_vrr;
      Double I_NAI_H2x3z_S_M1_vrr = PAX*I_NAI_Gx3z_S_M1_vrr-PNX*I_NAI_Gx3z_S_M2_vrr+oned2z*I_NAI_F3z_S_M1_vrr-oned2z*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_Hx4y_S_M1_vrr = PAX*I_NAI_G4y_S_M1_vrr-PNX*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_Hx4z_S_M1_vrr = PAX*I_NAI_G4z_S_M1_vrr-PNX*I_NAI_G4z_S_M2_vrr;
      Double I_NAI_H5y_S_M1_vrr = PAY*I_NAI_G4y_S_M1_vrr-PNY*I_NAI_G4y_S_M2_vrr+4*oned2z*I_NAI_F3y_S_M1_vrr-4*oned2z*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_H4yz_S_M1_vrr = PAZ*I_NAI_G4y_S_M1_vrr-PNZ*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_H3y2z_S_M1_vrr = PAZ*I_NAI_G3yz_S_M1_vrr-PNZ*I_NAI_G3yz_S_M2_vrr+oned2z*I_NAI_F3y_S_M1_vrr-oned2z*I_NAI_F3y_S_M2_vrr;
      Double I_NAI_H2y3z_S_M1_vrr = PAY*I_NAI_Gy3z_S_M1_vrr-PNY*I_NAI_Gy3z_S_M2_vrr+oned2z*I_NAI_F3z_S_M1_vrr-oned2z*I_NAI_F3z_S_M2_vrr;
      Double I_NAI_Hy4z_S_M1_vrr = PAY*I_NAI_G4z_S_M1_vrr-PNY*I_NAI_G4z_S_M2_vrr;
      Double I_NAI_H5z_S_M1_vrr = PAZ*I_NAI_G4z_S_M1_vrr-PNZ*I_NAI_G4z_S_M2_vrr+4*oned2z*I_NAI_F3z_S_M1_vrr-4*oned2z*I_NAI_F3z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_M1
       * expanding position: BRA1
       * code section is: VRR
       * totally 7 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_H_S_M1
       * RHS shell quartet name: SQ_NAI_H_S_M2
       * RHS shell quartet name: SQ_NAI_G_S_M1
       * RHS shell quartet name: SQ_NAI_G_S_M2
       ************************************************************/
      Double I_NAI_I6x_S_M1_vrr = PAX*I_NAI_H5x_S_M1_vrr-PNX*I_NAI_H5x_S_M2_vrr+5*oned2z*I_NAI_G4x_S_M1_vrr-5*oned2z*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_I5xy_S_M1_vrr = PAY*I_NAI_H5x_S_M1_vrr-PNY*I_NAI_H5x_S_M2_vrr;
      Double I_NAI_I5xz_S_M1_vrr = PAZ*I_NAI_H5x_S_M1_vrr-PNZ*I_NAI_H5x_S_M2_vrr;
      Double I_NAI_I4x2y_S_M1_vrr = PAY*I_NAI_H4xy_S_M1_vrr-PNY*I_NAI_H4xy_S_M2_vrr+oned2z*I_NAI_G4x_S_M1_vrr-oned2z*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_I4x2z_S_M1_vrr = PAZ*I_NAI_H4xz_S_M1_vrr-PNZ*I_NAI_H4xz_S_M2_vrr+oned2z*I_NAI_G4x_S_M1_vrr-oned2z*I_NAI_G4x_S_M2_vrr;
      Double I_NAI_I3x3y_S_M1_vrr = PAY*I_NAI_H3x2y_S_M1_vrr-PNY*I_NAI_H3x2y_S_M2_vrr+2*oned2z*I_NAI_G3xy_S_M1_vrr-2*oned2z*I_NAI_G3xy_S_M2_vrr;
      Double I_NAI_I3x2yz_S_M1_vrr = PAZ*I_NAI_H3x2y_S_M1_vrr-PNZ*I_NAI_H3x2y_S_M2_vrr;
      Double I_NAI_I3x3z_S_M1_vrr = PAZ*I_NAI_H3x2z_S_M1_vrr-PNZ*I_NAI_H3x2z_S_M2_vrr+2*oned2z*I_NAI_G3xz_S_M1_vrr-2*oned2z*I_NAI_G3xz_S_M2_vrr;
      Double I_NAI_I2x4y_S_M1_vrr = PAX*I_NAI_Hx4y_S_M1_vrr-PNX*I_NAI_Hx4y_S_M2_vrr+oned2z*I_NAI_G4y_S_M1_vrr-oned2z*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_I2x3yz_S_M1_vrr = PAZ*I_NAI_H2x3y_S_M1_vrr-PNZ*I_NAI_H2x3y_S_M2_vrr;
      Double I_NAI_I2xy3z_S_M1_vrr = PAY*I_NAI_H2x3z_S_M1_vrr-PNY*I_NAI_H2x3z_S_M2_vrr;
      Double I_NAI_I2x4z_S_M1_vrr = PAX*I_NAI_Hx4z_S_M1_vrr-PNX*I_NAI_Hx4z_S_M2_vrr+oned2z*I_NAI_G4z_S_M1_vrr-oned2z*I_NAI_G4z_S_M2_vrr;
      Double I_NAI_Ix5y_S_M1_vrr = PAX*I_NAI_H5y_S_M1_vrr-PNX*I_NAI_H5y_S_M2_vrr;
      Double I_NAI_Ix5z_S_M1_vrr = PAX*I_NAI_H5z_S_M1_vrr-PNX*I_NAI_H5z_S_M2_vrr;
      Double I_NAI_I6y_S_M1_vrr = PAY*I_NAI_H5y_S_M1_vrr-PNY*I_NAI_H5y_S_M2_vrr+5*oned2z*I_NAI_G4y_S_M1_vrr-5*oned2z*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_I5yz_S_M1_vrr = PAZ*I_NAI_H5y_S_M1_vrr-PNZ*I_NAI_H5y_S_M2_vrr;
      Double I_NAI_I4y2z_S_M1_vrr = PAZ*I_NAI_H4yz_S_M1_vrr-PNZ*I_NAI_H4yz_S_M2_vrr+oned2z*I_NAI_G4y_S_M1_vrr-oned2z*I_NAI_G4y_S_M2_vrr;
      Double I_NAI_I3y3z_S_M1_vrr = PAZ*I_NAI_H3y2z_S_M1_vrr-PNZ*I_NAI_H3y2z_S_M2_vrr+2*oned2z*I_NAI_G3yz_S_M1_vrr-2*oned2z*I_NAI_G3yz_S_M2_vrr;
      Double I_NAI_I2y4z_S_M1_vrr = PAY*I_NAI_Hy4z_S_M1_vrr-PNY*I_NAI_Hy4z_S_M2_vrr+oned2z*I_NAI_G4z_S_M1_vrr-oned2z*I_NAI_G4z_S_M2_vrr;
      Double I_NAI_Iy5z_S_M1_vrr = PAY*I_NAI_H5z_S_M1_vrr-PNY*I_NAI_H5z_S_M2_vrr;
      Double I_NAI_I6z_S_M1_vrr = PAZ*I_NAI_H5z_S_M1_vrr-PNZ*I_NAI_H5z_S_M2_vrr+5*oned2z*I_NAI_G4z_S_M1_vrr-5*oned2z*I_NAI_G4z_S_M2_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_P_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_Px_S_vrr = PAX*I_NAI_S_S_vrr-PNX*I_NAI_S_S_M1_vrr;
      Double I_NAI_Py_S_vrr = PAY*I_NAI_S_S_vrr-PNY*I_NAI_S_S_M1_vrr;
      Double I_NAI_Pz_S_vrr = PAZ*I_NAI_S_S_vrr-PNZ*I_NAI_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       * RHS shell quartet name: SQ_NAI_S_S
       * RHS shell quartet name: SQ_NAI_S_S_M1
       ************************************************************/
      Double I_NAI_D2x_S_vrr = PAX*I_NAI_Px_S_vrr-PNX*I_NAI_Px_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_Dxy_S_vrr = PAY*I_NAI_Px_S_vrr-PNY*I_NAI_Px_S_M1_vrr;
      Double I_NAI_Dxz_S_vrr = PAZ*I_NAI_Px_S_vrr-PNZ*I_NAI_Px_S_M1_vrr;
      Double I_NAI_D2y_S_vrr = PAY*I_NAI_Py_S_vrr-PNY*I_NAI_Py_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;
      Double I_NAI_Dyz_S_vrr = PAZ*I_NAI_Py_S_vrr-PNZ*I_NAI_Py_S_M1_vrr;
      Double I_NAI_D2z_S_vrr = PAZ*I_NAI_Pz_S_vrr-PNZ*I_NAI_Pz_S_M1_vrr+oned2z*I_NAI_S_S_vrr-oned2z*I_NAI_S_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_D_S
       * RHS shell quartet name: SQ_NAI_D_S_M1
       * RHS shell quartet name: SQ_NAI_P_S
       * RHS shell quartet name: SQ_NAI_P_S_M1
       ************************************************************/
      Double I_NAI_F3x_S_vrr = PAX*I_NAI_D2x_S_vrr-PNX*I_NAI_D2x_S_M1_vrr+2*oned2z*I_NAI_Px_S_vrr-2*oned2z*I_NAI_Px_S_M1_vrr;
      Double I_NAI_F2xy_S_vrr = PAY*I_NAI_D2x_S_vrr-PNY*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_F2xz_S_vrr = PAZ*I_NAI_D2x_S_vrr-PNZ*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_Fx2y_S_vrr = PAX*I_NAI_D2y_S_vrr-PNX*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Fxyz_S_vrr = PAZ*I_NAI_Dxy_S_vrr-PNZ*I_NAI_Dxy_S_M1_vrr;
      Double I_NAI_Fx2z_S_vrr = PAX*I_NAI_D2z_S_vrr-PNX*I_NAI_D2z_S_M1_vrr;
      Double I_NAI_F3y_S_vrr = PAY*I_NAI_D2y_S_vrr-PNY*I_NAI_D2y_S_M1_vrr+2*oned2z*I_NAI_Py_S_vrr-2*oned2z*I_NAI_Py_S_M1_vrr;
      Double I_NAI_F2yz_S_vrr = PAZ*I_NAI_D2y_S_vrr-PNZ*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Fy2z_S_vrr = PAY*I_NAI_D2z_S_vrr-PNY*I_NAI_D2z_S_M1_vrr;
      Double I_NAI_F3z_S_vrr = PAZ*I_NAI_D2z_S_vrr-PNZ*I_NAI_D2z_S_M1_vrr+2*oned2z*I_NAI_Pz_S_vrr-2*oned2z*I_NAI_Pz_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_F_S
       * RHS shell quartet name: SQ_NAI_F_S_M1
       * RHS shell quartet name: SQ_NAI_D_S
       * RHS shell quartet name: SQ_NAI_D_S_M1
       ************************************************************/
      Double I_NAI_G4x_S_vrr = PAX*I_NAI_F3x_S_vrr-PNX*I_NAI_F3x_S_M1_vrr+3*oned2z*I_NAI_D2x_S_vrr-3*oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_G3xy_S_vrr = PAY*I_NAI_F3x_S_vrr-PNY*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_G3xz_S_vrr = PAZ*I_NAI_F3x_S_vrr-PNZ*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_G2x2y_S_vrr = PAY*I_NAI_F2xy_S_vrr-PNY*I_NAI_F2xy_S_M1_vrr+oned2z*I_NAI_D2x_S_vrr-oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_G2xyz_S_vrr = PAZ*I_NAI_F2xy_S_vrr-PNZ*I_NAI_F2xy_S_M1_vrr;
      Double I_NAI_G2x2z_S_vrr = PAZ*I_NAI_F2xz_S_vrr-PNZ*I_NAI_F2xz_S_M1_vrr+oned2z*I_NAI_D2x_S_vrr-oned2z*I_NAI_D2x_S_M1_vrr;
      Double I_NAI_Gx3y_S_vrr = PAX*I_NAI_F3y_S_vrr-PNX*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_Gx2yz_S_vrr = PAZ*I_NAI_Fx2y_S_vrr-PNZ*I_NAI_Fx2y_S_M1_vrr;
      Double I_NAI_Gxy2z_S_vrr = PAY*I_NAI_Fx2z_S_vrr-PNY*I_NAI_Fx2z_S_M1_vrr;
      Double I_NAI_Gx3z_S_vrr = PAX*I_NAI_F3z_S_vrr-PNX*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_G4y_S_vrr = PAY*I_NAI_F3y_S_vrr-PNY*I_NAI_F3y_S_M1_vrr+3*oned2z*I_NAI_D2y_S_vrr-3*oned2z*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_G3yz_S_vrr = PAZ*I_NAI_F3y_S_vrr-PNZ*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_G2y2z_S_vrr = PAZ*I_NAI_F2yz_S_vrr-PNZ*I_NAI_F2yz_S_M1_vrr+oned2z*I_NAI_D2y_S_vrr-oned2z*I_NAI_D2y_S_M1_vrr;
      Double I_NAI_Gy3z_S_vrr = PAY*I_NAI_F3z_S_vrr-PNY*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_G4z_S_vrr = PAZ*I_NAI_F3z_S_vrr-PNZ*I_NAI_F3z_S_M1_vrr+3*oned2z*I_NAI_D2z_S_vrr-3*oned2z*I_NAI_D2z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_G_S
       * RHS shell quartet name: SQ_NAI_G_S_M1
       * RHS shell quartet name: SQ_NAI_F_S
       * RHS shell quartet name: SQ_NAI_F_S_M1
       ************************************************************/
      Double I_NAI_H5x_S_vrr = PAX*I_NAI_G4x_S_vrr-PNX*I_NAI_G4x_S_M1_vrr+4*oned2z*I_NAI_F3x_S_vrr-4*oned2z*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_H4xy_S_vrr = PAY*I_NAI_G4x_S_vrr-PNY*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_H4xz_S_vrr = PAZ*I_NAI_G4x_S_vrr-PNZ*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_H3x2y_S_vrr = PAY*I_NAI_G3xy_S_vrr-PNY*I_NAI_G3xy_S_M1_vrr+oned2z*I_NAI_F3x_S_vrr-oned2z*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_H3xyz_S_vrr = PAZ*I_NAI_G3xy_S_vrr-PNZ*I_NAI_G3xy_S_M1_vrr;
      Double I_NAI_H3x2z_S_vrr = PAZ*I_NAI_G3xz_S_vrr-PNZ*I_NAI_G3xz_S_M1_vrr+oned2z*I_NAI_F3x_S_vrr-oned2z*I_NAI_F3x_S_M1_vrr;
      Double I_NAI_H2x3y_S_vrr = PAX*I_NAI_Gx3y_S_vrr-PNX*I_NAI_Gx3y_S_M1_vrr+oned2z*I_NAI_F3y_S_vrr-oned2z*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_H2x2yz_S_vrr = PAZ*I_NAI_G2x2y_S_vrr-PNZ*I_NAI_G2x2y_S_M1_vrr;
      Double I_NAI_H2xy2z_S_vrr = PAY*I_NAI_G2x2z_S_vrr-PNY*I_NAI_G2x2z_S_M1_vrr;
      Double I_NAI_H2x3z_S_vrr = PAX*I_NAI_Gx3z_S_vrr-PNX*I_NAI_Gx3z_S_M1_vrr+oned2z*I_NAI_F3z_S_vrr-oned2z*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_Hx4y_S_vrr = PAX*I_NAI_G4y_S_vrr-PNX*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_Hx3yz_S_vrr = PAZ*I_NAI_Gx3y_S_vrr-PNZ*I_NAI_Gx3y_S_M1_vrr;
      Double I_NAI_Hx2y2z_S_vrr = PAX*I_NAI_G2y2z_S_vrr-PNX*I_NAI_G2y2z_S_M1_vrr;
      Double I_NAI_Hxy3z_S_vrr = PAY*I_NAI_Gx3z_S_vrr-PNY*I_NAI_Gx3z_S_M1_vrr;
      Double I_NAI_Hx4z_S_vrr = PAX*I_NAI_G4z_S_vrr-PNX*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_H5y_S_vrr = PAY*I_NAI_G4y_S_vrr-PNY*I_NAI_G4y_S_M1_vrr+4*oned2z*I_NAI_F3y_S_vrr-4*oned2z*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_H4yz_S_vrr = PAZ*I_NAI_G4y_S_vrr-PNZ*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_H3y2z_S_vrr = PAZ*I_NAI_G3yz_S_vrr-PNZ*I_NAI_G3yz_S_M1_vrr+oned2z*I_NAI_F3y_S_vrr-oned2z*I_NAI_F3y_S_M1_vrr;
      Double I_NAI_H2y3z_S_vrr = PAY*I_NAI_Gy3z_S_vrr-PNY*I_NAI_Gy3z_S_M1_vrr+oned2z*I_NAI_F3z_S_vrr-oned2z*I_NAI_F3z_S_M1_vrr;
      Double I_NAI_Hy4z_S_vrr = PAY*I_NAI_G4z_S_vrr-PNY*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_H5z_S_vrr = PAZ*I_NAI_G4z_S_vrr-PNZ*I_NAI_G4z_S_M1_vrr+4*oned2z*I_NAI_F3z_S_vrr-4*oned2z*I_NAI_F3z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_H_S
       * RHS shell quartet name: SQ_NAI_H_S_M1
       * RHS shell quartet name: SQ_NAI_G_S
       * RHS shell quartet name: SQ_NAI_G_S_M1
       ************************************************************/
      Double I_NAI_I6x_S_vrr = PAX*I_NAI_H5x_S_vrr-PNX*I_NAI_H5x_S_M1_vrr+5*oned2z*I_NAI_G4x_S_vrr-5*oned2z*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_I5xy_S_vrr = PAY*I_NAI_H5x_S_vrr-PNY*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_I5xz_S_vrr = PAZ*I_NAI_H5x_S_vrr-PNZ*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_I4x2y_S_vrr = PAY*I_NAI_H4xy_S_vrr-PNY*I_NAI_H4xy_S_M1_vrr+oned2z*I_NAI_G4x_S_vrr-oned2z*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_I4xyz_S_vrr = PAZ*I_NAI_H4xy_S_vrr-PNZ*I_NAI_H4xy_S_M1_vrr;
      Double I_NAI_I4x2z_S_vrr = PAZ*I_NAI_H4xz_S_vrr-PNZ*I_NAI_H4xz_S_M1_vrr+oned2z*I_NAI_G4x_S_vrr-oned2z*I_NAI_G4x_S_M1_vrr;
      Double I_NAI_I3x3y_S_vrr = PAY*I_NAI_H3x2y_S_vrr-PNY*I_NAI_H3x2y_S_M1_vrr+2*oned2z*I_NAI_G3xy_S_vrr-2*oned2z*I_NAI_G3xy_S_M1_vrr;
      Double I_NAI_I3x2yz_S_vrr = PAZ*I_NAI_H3x2y_S_vrr-PNZ*I_NAI_H3x2y_S_M1_vrr;
      Double I_NAI_I3xy2z_S_vrr = PAY*I_NAI_H3x2z_S_vrr-PNY*I_NAI_H3x2z_S_M1_vrr;
      Double I_NAI_I3x3z_S_vrr = PAZ*I_NAI_H3x2z_S_vrr-PNZ*I_NAI_H3x2z_S_M1_vrr+2*oned2z*I_NAI_G3xz_S_vrr-2*oned2z*I_NAI_G3xz_S_M1_vrr;
      Double I_NAI_I2x4y_S_vrr = PAX*I_NAI_Hx4y_S_vrr-PNX*I_NAI_Hx4y_S_M1_vrr+oned2z*I_NAI_G4y_S_vrr-oned2z*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_I2x3yz_S_vrr = PAZ*I_NAI_H2x3y_S_vrr-PNZ*I_NAI_H2x3y_S_M1_vrr;
      Double I_NAI_I2x2y2z_S_vrr = PAZ*I_NAI_H2x2yz_S_vrr-PNZ*I_NAI_H2x2yz_S_M1_vrr+oned2z*I_NAI_G2x2y_S_vrr-oned2z*I_NAI_G2x2y_S_M1_vrr;
      Double I_NAI_I2xy3z_S_vrr = PAY*I_NAI_H2x3z_S_vrr-PNY*I_NAI_H2x3z_S_M1_vrr;
      Double I_NAI_I2x4z_S_vrr = PAX*I_NAI_Hx4z_S_vrr-PNX*I_NAI_Hx4z_S_M1_vrr+oned2z*I_NAI_G4z_S_vrr-oned2z*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_Ix5y_S_vrr = PAX*I_NAI_H5y_S_vrr-PNX*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_Ix4yz_S_vrr = PAZ*I_NAI_Hx4y_S_vrr-PNZ*I_NAI_Hx4y_S_M1_vrr;
      Double I_NAI_Ix3y2z_S_vrr = PAX*I_NAI_H3y2z_S_vrr-PNX*I_NAI_H3y2z_S_M1_vrr;
      Double I_NAI_Ix2y3z_S_vrr = PAX*I_NAI_H2y3z_S_vrr-PNX*I_NAI_H2y3z_S_M1_vrr;
      Double I_NAI_Ixy4z_S_vrr = PAY*I_NAI_Hx4z_S_vrr-PNY*I_NAI_Hx4z_S_M1_vrr;
      Double I_NAI_Ix5z_S_vrr = PAX*I_NAI_H5z_S_vrr-PNX*I_NAI_H5z_S_M1_vrr;
      Double I_NAI_I6y_S_vrr = PAY*I_NAI_H5y_S_vrr-PNY*I_NAI_H5y_S_M1_vrr+5*oned2z*I_NAI_G4y_S_vrr-5*oned2z*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_I5yz_S_vrr = PAZ*I_NAI_H5y_S_vrr-PNZ*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_I4y2z_S_vrr = PAZ*I_NAI_H4yz_S_vrr-PNZ*I_NAI_H4yz_S_M1_vrr+oned2z*I_NAI_G4y_S_vrr-oned2z*I_NAI_G4y_S_M1_vrr;
      Double I_NAI_I3y3z_S_vrr = PAZ*I_NAI_H3y2z_S_vrr-PNZ*I_NAI_H3y2z_S_M1_vrr+2*oned2z*I_NAI_G3yz_S_vrr-2*oned2z*I_NAI_G3yz_S_M1_vrr;
      Double I_NAI_I2y4z_S_vrr = PAY*I_NAI_Hy4z_S_vrr-PNY*I_NAI_Hy4z_S_M1_vrr+oned2z*I_NAI_G4z_S_vrr-oned2z*I_NAI_G4z_S_M1_vrr;
      Double I_NAI_Iy5z_S_vrr = PAY*I_NAI_H5z_S_vrr-PNY*I_NAI_H5z_S_M1_vrr;
      Double I_NAI_I6z_S_vrr = PAZ*I_NAI_H5z_S_vrr-PNZ*I_NAI_H5z_S_M1_vrr+5*oned2z*I_NAI_G4z_S_vrr-5*oned2z*I_NAI_G4z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S
       * expanding position: BRA1
       * code section is: VRR
       * totally 0 integrals are omitted 
       * RHS shell quartet name: SQ_NAI_I_S
       * RHS shell quartet name: SQ_NAI_I_S_M1
       * RHS shell quartet name: SQ_NAI_H_S
       * RHS shell quartet name: SQ_NAI_H_S_M1
       ************************************************************/
      Double I_NAI_K7x_S_vrr = PAX*I_NAI_I6x_S_vrr-PNX*I_NAI_I6x_S_M1_vrr+6*oned2z*I_NAI_H5x_S_vrr-6*oned2z*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_K6xy_S_vrr = PAY*I_NAI_I6x_S_vrr-PNY*I_NAI_I6x_S_M1_vrr;
      Double I_NAI_K6xz_S_vrr = PAZ*I_NAI_I6x_S_vrr-PNZ*I_NAI_I6x_S_M1_vrr;
      Double I_NAI_K5x2y_S_vrr = PAY*I_NAI_I5xy_S_vrr-PNY*I_NAI_I5xy_S_M1_vrr+oned2z*I_NAI_H5x_S_vrr-oned2z*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_K5xyz_S_vrr = PAZ*I_NAI_I5xy_S_vrr-PNZ*I_NAI_I5xy_S_M1_vrr;
      Double I_NAI_K5x2z_S_vrr = PAZ*I_NAI_I5xz_S_vrr-PNZ*I_NAI_I5xz_S_M1_vrr+oned2z*I_NAI_H5x_S_vrr-oned2z*I_NAI_H5x_S_M1_vrr;
      Double I_NAI_K4x3y_S_vrr = PAY*I_NAI_I4x2y_S_vrr-PNY*I_NAI_I4x2y_S_M1_vrr+2*oned2z*I_NAI_H4xy_S_vrr-2*oned2z*I_NAI_H4xy_S_M1_vrr;
      Double I_NAI_K4x2yz_S_vrr = PAZ*I_NAI_I4x2y_S_vrr-PNZ*I_NAI_I4x2y_S_M1_vrr;
      Double I_NAI_K4xy2z_S_vrr = PAY*I_NAI_I4x2z_S_vrr-PNY*I_NAI_I4x2z_S_M1_vrr;
      Double I_NAI_K4x3z_S_vrr = PAZ*I_NAI_I4x2z_S_vrr-PNZ*I_NAI_I4x2z_S_M1_vrr+2*oned2z*I_NAI_H4xz_S_vrr-2*oned2z*I_NAI_H4xz_S_M1_vrr;
      Double I_NAI_K3x4y_S_vrr = PAX*I_NAI_I2x4y_S_vrr-PNX*I_NAI_I2x4y_S_M1_vrr+2*oned2z*I_NAI_Hx4y_S_vrr-2*oned2z*I_NAI_Hx4y_S_M1_vrr;
      Double I_NAI_K3x3yz_S_vrr = PAZ*I_NAI_I3x3y_S_vrr-PNZ*I_NAI_I3x3y_S_M1_vrr;
      Double I_NAI_K3x2y2z_S_vrr = PAZ*I_NAI_I3x2yz_S_vrr-PNZ*I_NAI_I3x2yz_S_M1_vrr+oned2z*I_NAI_H3x2y_S_vrr-oned2z*I_NAI_H3x2y_S_M1_vrr;
      Double I_NAI_K3xy3z_S_vrr = PAY*I_NAI_I3x3z_S_vrr-PNY*I_NAI_I3x3z_S_M1_vrr;
      Double I_NAI_K3x4z_S_vrr = PAX*I_NAI_I2x4z_S_vrr-PNX*I_NAI_I2x4z_S_M1_vrr+2*oned2z*I_NAI_Hx4z_S_vrr-2*oned2z*I_NAI_Hx4z_S_M1_vrr;
      Double I_NAI_K2x5y_S_vrr = PAX*I_NAI_Ix5y_S_vrr-PNX*I_NAI_Ix5y_S_M1_vrr+oned2z*I_NAI_H5y_S_vrr-oned2z*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_K2x4yz_S_vrr = PAZ*I_NAI_I2x4y_S_vrr-PNZ*I_NAI_I2x4y_S_M1_vrr;
      Double I_NAI_K2x3y2z_S_vrr = PAZ*I_NAI_I2x3yz_S_vrr-PNZ*I_NAI_I2x3yz_S_M1_vrr+oned2z*I_NAI_H2x3y_S_vrr-oned2z*I_NAI_H2x3y_S_M1_vrr;
      Double I_NAI_K2x2y3z_S_vrr = PAY*I_NAI_I2xy3z_S_vrr-PNY*I_NAI_I2xy3z_S_M1_vrr+oned2z*I_NAI_H2x3z_S_vrr-oned2z*I_NAI_H2x3z_S_M1_vrr;
      Double I_NAI_K2xy4z_S_vrr = PAY*I_NAI_I2x4z_S_vrr-PNY*I_NAI_I2x4z_S_M1_vrr;
      Double I_NAI_K2x5z_S_vrr = PAX*I_NAI_Ix5z_S_vrr-PNX*I_NAI_Ix5z_S_M1_vrr+oned2z*I_NAI_H5z_S_vrr-oned2z*I_NAI_H5z_S_M1_vrr;
      Double I_NAI_Kx6y_S_vrr = PAX*I_NAI_I6y_S_vrr-PNX*I_NAI_I6y_S_M1_vrr;
      Double I_NAI_Kx5yz_S_vrr = PAZ*I_NAI_Ix5y_S_vrr-PNZ*I_NAI_Ix5y_S_M1_vrr;
      Double I_NAI_Kx4y2z_S_vrr = PAX*I_NAI_I4y2z_S_vrr-PNX*I_NAI_I4y2z_S_M1_vrr;
      Double I_NAI_Kx3y3z_S_vrr = PAX*I_NAI_I3y3z_S_vrr-PNX*I_NAI_I3y3z_S_M1_vrr;
      Double I_NAI_Kx2y4z_S_vrr = PAX*I_NAI_I2y4z_S_vrr-PNX*I_NAI_I2y4z_S_M1_vrr;
      Double I_NAI_Kxy5z_S_vrr = PAY*I_NAI_Ix5z_S_vrr-PNY*I_NAI_Ix5z_S_M1_vrr;
      Double I_NAI_Kx6z_S_vrr = PAX*I_NAI_I6z_S_vrr-PNX*I_NAI_I6z_S_M1_vrr;
      Double I_NAI_K7y_S_vrr = PAY*I_NAI_I6y_S_vrr-PNY*I_NAI_I6y_S_M1_vrr+6*oned2z*I_NAI_H5y_S_vrr-6*oned2z*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_K6yz_S_vrr = PAZ*I_NAI_I6y_S_vrr-PNZ*I_NAI_I6y_S_M1_vrr;
      Double I_NAI_K5y2z_S_vrr = PAZ*I_NAI_I5yz_S_vrr-PNZ*I_NAI_I5yz_S_M1_vrr+oned2z*I_NAI_H5y_S_vrr-oned2z*I_NAI_H5y_S_M1_vrr;
      Double I_NAI_K4y3z_S_vrr = PAZ*I_NAI_I4y2z_S_vrr-PNZ*I_NAI_I4y2z_S_M1_vrr+2*oned2z*I_NAI_H4yz_S_vrr-2*oned2z*I_NAI_H4yz_S_M1_vrr;
      Double I_NAI_K3y4z_S_vrr = PAY*I_NAI_I2y4z_S_vrr-PNY*I_NAI_I2y4z_S_M1_vrr+2*oned2z*I_NAI_Hy4z_S_vrr-2*oned2z*I_NAI_Hy4z_S_M1_vrr;
      Double I_NAI_K2y5z_S_vrr = PAY*I_NAI_Iy5z_S_vrr-PNY*I_NAI_Iy5z_S_M1_vrr+oned2z*I_NAI_H5z_S_vrr-oned2z*I_NAI_H5z_S_M1_vrr;
      Double I_NAI_Ky6z_S_vrr = PAY*I_NAI_I6z_S_vrr-PNY*I_NAI_I6z_S_M1_vrr;
      Double I_NAI_K7z_S_vrr = PAZ*I_NAI_I6z_S_vrr-PNZ*I_NAI_I6z_S_M1_vrr+6*oned2z*I_NAI_H5z_S_vrr-6*oned2z*I_NAI_H5z_S_M1_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_H5x_S += I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S += I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S += I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S += I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S += I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S += I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S += I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S += I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S += I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S += I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S += I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S += I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S += I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S += I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S += I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S += I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S += I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S += I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S += I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S += I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S += I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_G4x_S += I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S += I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S += I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S += I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S += I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S += I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S += I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S += I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S += I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S += I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S += I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S += I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S += I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S += I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S += I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_F3x_S += I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S += I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S += I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S += I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S += I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S += I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S += I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S += I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S += I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S += I_NAI_F3z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_K_S_a_coefs = alpha;
      I_NAI_K7x_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S_a += SQ_NAI_K_S_a_coefs*I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_a_coefs = alpha;
      I_NAI_I6x_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_a += SQ_NAI_I_S_a_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_a_coefs = alpha;
      I_NAI_H5x_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_a += SQ_NAI_H_S_a_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_a
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_a_coefs = alpha;
      I_NAI_G4x_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_a += SQ_NAI_G_S_a_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_D_S
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      I_NAI_D2x_S += I_NAI_D2x_S_vrr;
      I_NAI_Dxy_S += I_NAI_Dxy_S_vrr;
      I_NAI_Dxz_S += I_NAI_Dxz_S_vrr;
      I_NAI_D2y_S += I_NAI_D2y_S_vrr;
      I_NAI_Dyz_S += I_NAI_Dyz_S_vrr;
      I_NAI_D2z_S += I_NAI_D2z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_K_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_K_S_b_coefs = beta;
      I_NAI_K7x_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K7x_S_vrr;
      I_NAI_K6xy_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K6xy_S_vrr;
      I_NAI_K6xz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K6xz_S_vrr;
      I_NAI_K5x2y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K5x2y_S_vrr;
      I_NAI_K5xyz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K5xyz_S_vrr;
      I_NAI_K5x2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K5x2z_S_vrr;
      I_NAI_K4x3y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4x3y_S_vrr;
      I_NAI_K4x2yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4x2yz_S_vrr;
      I_NAI_K4xy2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4xy2z_S_vrr;
      I_NAI_K4x3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4x3z_S_vrr;
      I_NAI_K3x4y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3x4y_S_vrr;
      I_NAI_K3x3yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3x3yz_S_vrr;
      I_NAI_K3x2y2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3x2y2z_S_vrr;
      I_NAI_K3xy3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3xy3z_S_vrr;
      I_NAI_K3x4z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3x4z_S_vrr;
      I_NAI_K2x5y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x5y_S_vrr;
      I_NAI_K2x4yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x4yz_S_vrr;
      I_NAI_K2x3y2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x3y2z_S_vrr;
      I_NAI_K2x2y3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x2y3z_S_vrr;
      I_NAI_K2xy4z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2xy4z_S_vrr;
      I_NAI_K2x5z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2x5z_S_vrr;
      I_NAI_Kx6y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx6y_S_vrr;
      I_NAI_Kx5yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx5yz_S_vrr;
      I_NAI_Kx4y2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx4y2z_S_vrr;
      I_NAI_Kx3y3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx3y3z_S_vrr;
      I_NAI_Kx2y4z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx2y4z_S_vrr;
      I_NAI_Kxy5z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kxy5z_S_vrr;
      I_NAI_Kx6z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Kx6z_S_vrr;
      I_NAI_K7y_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K7y_S_vrr;
      I_NAI_K6yz_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K6yz_S_vrr;
      I_NAI_K5y2z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K5y2z_S_vrr;
      I_NAI_K4y3z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K4y3z_S_vrr;
      I_NAI_K3y4z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K3y4z_S_vrr;
      I_NAI_K2y5z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K2y5z_S_vrr;
      I_NAI_Ky6z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_Ky6z_S_vrr;
      I_NAI_K7z_S_b += SQ_NAI_K_S_b_coefs*I_NAI_K7z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_I_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_I_S_b_coefs = beta;
      I_NAI_I6x_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I6x_S_vrr;
      I_NAI_I5xy_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I5xy_S_vrr;
      I_NAI_I5xz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I5xz_S_vrr;
      I_NAI_I4x2y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I4x2y_S_vrr;
      I_NAI_I4xyz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I4xyz_S_vrr;
      I_NAI_I4x2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I4x2z_S_vrr;
      I_NAI_I3x3y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3x3y_S_vrr;
      I_NAI_I3x2yz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3x2yz_S_vrr;
      I_NAI_I3xy2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3xy2z_S_vrr;
      I_NAI_I3x3z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3x3z_S_vrr;
      I_NAI_I2x4y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2x4y_S_vrr;
      I_NAI_I2x3yz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2x3yz_S_vrr;
      I_NAI_I2x2y2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2x2y2z_S_vrr;
      I_NAI_I2xy3z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2xy3z_S_vrr;
      I_NAI_I2x4z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2x4z_S_vrr;
      I_NAI_Ix5y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix5y_S_vrr;
      I_NAI_Ix4yz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix4yz_S_vrr;
      I_NAI_Ix3y2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix3y2z_S_vrr;
      I_NAI_Ix2y3z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix2y3z_S_vrr;
      I_NAI_Ixy4z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ixy4z_S_vrr;
      I_NAI_Ix5z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Ix5z_S_vrr;
      I_NAI_I6y_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I6y_S_vrr;
      I_NAI_I5yz_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I5yz_S_vrr;
      I_NAI_I4y2z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I4y2z_S_vrr;
      I_NAI_I3y3z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I3y3z_S_vrr;
      I_NAI_I2y4z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I2y4z_S_vrr;
      I_NAI_Iy5z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_Iy5z_S_vrr;
      I_NAI_I6z_S_b += SQ_NAI_I_S_b_coefs*I_NAI_I6z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_H_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_H_S_b_coefs = beta;
      I_NAI_H5x_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H5x_S_vrr;
      I_NAI_H4xy_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H4xy_S_vrr;
      I_NAI_H4xz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H4xz_S_vrr;
      I_NAI_H3x2y_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H3x2y_S_vrr;
      I_NAI_H3xyz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H3xyz_S_vrr;
      I_NAI_H3x2z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H3x2z_S_vrr;
      I_NAI_H2x3y_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2x3y_S_vrr;
      I_NAI_H2x2yz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2x2yz_S_vrr;
      I_NAI_H2xy2z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2xy2z_S_vrr;
      I_NAI_H2x3z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2x3z_S_vrr;
      I_NAI_Hx4y_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hx4y_S_vrr;
      I_NAI_Hx3yz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hx3yz_S_vrr;
      I_NAI_Hx2y2z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hx2y2z_S_vrr;
      I_NAI_Hxy3z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hxy3z_S_vrr;
      I_NAI_Hx4z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hx4z_S_vrr;
      I_NAI_H5y_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H5y_S_vrr;
      I_NAI_H4yz_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H4yz_S_vrr;
      I_NAI_H3y2z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H3y2z_S_vrr;
      I_NAI_H2y3z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H2y3z_S_vrr;
      I_NAI_Hy4z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_Hy4z_S_vrr;
      I_NAI_H5z_S_b += SQ_NAI_H_S_b_coefs*I_NAI_H5z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_G_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_G_S_b_coefs = beta;
      I_NAI_G4x_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G4x_S_vrr;
      I_NAI_G3xy_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G3xy_S_vrr;
      I_NAI_G3xz_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G3xz_S_vrr;
      I_NAI_G2x2y_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G2x2y_S_vrr;
      I_NAI_G2xyz_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G2xyz_S_vrr;
      I_NAI_G2x2z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G2x2z_S_vrr;
      I_NAI_Gx3y_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gx3y_S_vrr;
      I_NAI_Gx2yz_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gx2yz_S_vrr;
      I_NAI_Gxy2z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gxy2z_S_vrr;
      I_NAI_Gx3z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gx3z_S_vrr;
      I_NAI_G4y_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G4y_S_vrr;
      I_NAI_G3yz_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G3yz_S_vrr;
      I_NAI_G2y2z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G2y2z_S_vrr;
      I_NAI_Gy3z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_Gy3z_S_vrr;
      I_NAI_G4z_S_b += SQ_NAI_G_S_b_coefs*I_NAI_G4z_S_vrr;

      /************************************************************
       * shell quartet name: SQ_NAI_F_S_b
       * doing contraction work for VRR part 
       * totally 0 integrals are omitted 
       ************************************************************/
      Double SQ_NAI_F_S_b_coefs = beta;
      I_NAI_F3x_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F3x_S_vrr;
      I_NAI_F2xy_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F2xy_S_vrr;
      I_NAI_F2xz_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F2xz_S_vrr;
      I_NAI_Fx2y_S_b += SQ_NAI_F_S_b_coefs*I_NAI_Fx2y_S_vrr;
      I_NAI_Fxyz_S_b += SQ_NAI_F_S_b_coefs*I_NAI_Fxyz_S_vrr;
      I_NAI_Fx2z_S_b += SQ_NAI_F_S_b_coefs*I_NAI_Fx2z_S_vrr;
      I_NAI_F3y_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F3y_S_vrr;
      I_NAI_F2yz_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F2yz_S_vrr;
      I_NAI_Fy2z_S_b += SQ_NAI_F_S_b_coefs*I_NAI_Fy2z_S_vrr;
      I_NAI_F3z_S_b += SQ_NAI_F_S_b_coefs*I_NAI_F3z_S_vrr;
    }
  }

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
   * shell quartet name: SQ_NAI_D_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_S
   * RHS shell quartet name: SQ_NAI_D_S
   ************************************************************/
  Double I_NAI_D2x_Px = I_NAI_F3x_S+ABX*I_NAI_D2x_S;
  Double I_NAI_Dxy_Px = I_NAI_F2xy_S+ABX*I_NAI_Dxy_S;
  Double I_NAI_Dxz_Px = I_NAI_F2xz_S+ABX*I_NAI_Dxz_S;
  Double I_NAI_D2y_Px = I_NAI_Fx2y_S+ABX*I_NAI_D2y_S;
  Double I_NAI_Dyz_Px = I_NAI_Fxyz_S+ABX*I_NAI_Dyz_S;
  Double I_NAI_D2z_Px = I_NAI_Fx2z_S+ABX*I_NAI_D2z_S;
  Double I_NAI_D2x_Py = I_NAI_F2xy_S+ABY*I_NAI_D2x_S;
  Double I_NAI_Dxy_Py = I_NAI_Fx2y_S+ABY*I_NAI_Dxy_S;
  Double I_NAI_Dxz_Py = I_NAI_Fxyz_S+ABY*I_NAI_Dxz_S;
  Double I_NAI_D2y_Py = I_NAI_F3y_S+ABY*I_NAI_D2y_S;
  Double I_NAI_Dyz_Py = I_NAI_F2yz_S+ABY*I_NAI_Dyz_S;
  Double I_NAI_D2z_Py = I_NAI_Fy2z_S+ABY*I_NAI_D2z_S;
  Double I_NAI_D2x_Pz = I_NAI_F2xz_S+ABZ*I_NAI_D2x_S;
  Double I_NAI_Dxy_Pz = I_NAI_Fxyz_S+ABZ*I_NAI_Dxy_S;
  Double I_NAI_Dxz_Pz = I_NAI_Fx2z_S+ABZ*I_NAI_Dxz_S;
  Double I_NAI_D2y_Pz = I_NAI_F2yz_S+ABZ*I_NAI_D2y_S;
  Double I_NAI_Dyz_Pz = I_NAI_Fy2z_S+ABZ*I_NAI_Dyz_S;
  Double I_NAI_D2z_Pz = I_NAI_F3z_S+ABZ*I_NAI_D2z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S
   * RHS shell quartet name: SQ_NAI_F_S
   ************************************************************/
  Double I_NAI_F3x_Px = I_NAI_G4x_S+ABX*I_NAI_F3x_S;
  Double I_NAI_F2xy_Px = I_NAI_G3xy_S+ABX*I_NAI_F2xy_S;
  Double I_NAI_F2xz_Px = I_NAI_G3xz_S+ABX*I_NAI_F2xz_S;
  Double I_NAI_Fx2y_Px = I_NAI_G2x2y_S+ABX*I_NAI_Fx2y_S;
  Double I_NAI_Fxyz_Px = I_NAI_G2xyz_S+ABX*I_NAI_Fxyz_S;
  Double I_NAI_Fx2z_Px = I_NAI_G2x2z_S+ABX*I_NAI_Fx2z_S;
  Double I_NAI_F3y_Px = I_NAI_Gx3y_S+ABX*I_NAI_F3y_S;
  Double I_NAI_F2yz_Px = I_NAI_Gx2yz_S+ABX*I_NAI_F2yz_S;
  Double I_NAI_Fy2z_Px = I_NAI_Gxy2z_S+ABX*I_NAI_Fy2z_S;
  Double I_NAI_F3z_Px = I_NAI_Gx3z_S+ABX*I_NAI_F3z_S;
  Double I_NAI_F3x_Py = I_NAI_G3xy_S+ABY*I_NAI_F3x_S;
  Double I_NAI_F2xy_Py = I_NAI_G2x2y_S+ABY*I_NAI_F2xy_S;
  Double I_NAI_F2xz_Py = I_NAI_G2xyz_S+ABY*I_NAI_F2xz_S;
  Double I_NAI_Fx2y_Py = I_NAI_Gx3y_S+ABY*I_NAI_Fx2y_S;
  Double I_NAI_Fxyz_Py = I_NAI_Gx2yz_S+ABY*I_NAI_Fxyz_S;
  Double I_NAI_Fx2z_Py = I_NAI_Gxy2z_S+ABY*I_NAI_Fx2z_S;
  Double I_NAI_F3y_Py = I_NAI_G4y_S+ABY*I_NAI_F3y_S;
  Double I_NAI_F2yz_Py = I_NAI_G3yz_S+ABY*I_NAI_F2yz_S;
  Double I_NAI_Fy2z_Py = I_NAI_G2y2z_S+ABY*I_NAI_Fy2z_S;
  Double I_NAI_F3z_Py = I_NAI_Gy3z_S+ABY*I_NAI_F3z_S;
  Double I_NAI_F3x_Pz = I_NAI_G3xz_S+ABZ*I_NAI_F3x_S;
  Double I_NAI_F2xy_Pz = I_NAI_G2xyz_S+ABZ*I_NAI_F2xy_S;
  Double I_NAI_F2xz_Pz = I_NAI_G2x2z_S+ABZ*I_NAI_F2xz_S;
  Double I_NAI_Fx2y_Pz = I_NAI_Gx2yz_S+ABZ*I_NAI_Fx2y_S;
  Double I_NAI_Fxyz_Pz = I_NAI_Gxy2z_S+ABZ*I_NAI_Fxyz_S;
  Double I_NAI_Fx2z_Pz = I_NAI_Gx3z_S+ABZ*I_NAI_Fx2z_S;
  Double I_NAI_F3y_Pz = I_NAI_G3yz_S+ABZ*I_NAI_F3y_S;
  Double I_NAI_F2yz_Pz = I_NAI_G2y2z_S+ABZ*I_NAI_F2yz_S;
  Double I_NAI_Fy2z_Pz = I_NAI_Gy3z_S+ABZ*I_NAI_Fy2z_S;
  Double I_NAI_F3z_Pz = I_NAI_G4z_S+ABZ*I_NAI_F3z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_D_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_P
   * RHS shell quartet name: SQ_NAI_D_P
   ************************************************************/
  Double I_NAI_D2x_D2x = I_NAI_F3x_Px+ABX*I_NAI_D2x_Px;
  Double I_NAI_Dxy_D2x = I_NAI_F2xy_Px+ABX*I_NAI_Dxy_Px;
  Double I_NAI_Dxz_D2x = I_NAI_F2xz_Px+ABX*I_NAI_Dxz_Px;
  Double I_NAI_D2y_D2x = I_NAI_Fx2y_Px+ABX*I_NAI_D2y_Px;
  Double I_NAI_Dyz_D2x = I_NAI_Fxyz_Px+ABX*I_NAI_Dyz_Px;
  Double I_NAI_D2z_D2x = I_NAI_Fx2z_Px+ABX*I_NAI_D2z_Px;
  Double I_NAI_D2x_Dxy = I_NAI_F2xy_Px+ABY*I_NAI_D2x_Px;
  Double I_NAI_Dxy_Dxy = I_NAI_Fx2y_Px+ABY*I_NAI_Dxy_Px;
  Double I_NAI_Dxz_Dxy = I_NAI_Fxyz_Px+ABY*I_NAI_Dxz_Px;
  Double I_NAI_D2y_Dxy = I_NAI_F3y_Px+ABY*I_NAI_D2y_Px;
  Double I_NAI_Dyz_Dxy = I_NAI_F2yz_Px+ABY*I_NAI_Dyz_Px;
  Double I_NAI_D2z_Dxy = I_NAI_Fy2z_Px+ABY*I_NAI_D2z_Px;
  Double I_NAI_D2x_D2y = I_NAI_F2xy_Py+ABY*I_NAI_D2x_Py;
  Double I_NAI_Dxy_D2y = I_NAI_Fx2y_Py+ABY*I_NAI_Dxy_Py;
  Double I_NAI_Dxz_D2y = I_NAI_Fxyz_Py+ABY*I_NAI_Dxz_Py;
  Double I_NAI_D2y_D2y = I_NAI_F3y_Py+ABY*I_NAI_D2y_Py;
  Double I_NAI_Dyz_D2y = I_NAI_F2yz_Py+ABY*I_NAI_Dyz_Py;
  Double I_NAI_D2z_D2y = I_NAI_Fy2z_Py+ABY*I_NAI_D2z_Py;
  Double I_NAI_D2x_D2z = I_NAI_F2xz_Pz+ABZ*I_NAI_D2x_Pz;
  Double I_NAI_Dxy_D2z = I_NAI_Fxyz_Pz+ABZ*I_NAI_Dxy_Pz;
  Double I_NAI_Dxz_D2z = I_NAI_Fx2z_Pz+ABZ*I_NAI_Dxz_Pz;
  Double I_NAI_D2y_D2z = I_NAI_F2yz_Pz+ABZ*I_NAI_D2y_Pz;
  Double I_NAI_Dyz_D2z = I_NAI_Fy2z_Pz+ABZ*I_NAI_Dyz_Pz;
  Double I_NAI_D2z_D2z = I_NAI_F3z_Pz+ABZ*I_NAI_D2z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S
   * RHS shell quartet name: SQ_NAI_G_S
   ************************************************************/
  Double I_NAI_G4x_Px = I_NAI_H5x_S+ABX*I_NAI_G4x_S;
  Double I_NAI_G3xy_Px = I_NAI_H4xy_S+ABX*I_NAI_G3xy_S;
  Double I_NAI_G3xz_Px = I_NAI_H4xz_S+ABX*I_NAI_G3xz_S;
  Double I_NAI_G2x2y_Px = I_NAI_H3x2y_S+ABX*I_NAI_G2x2y_S;
  Double I_NAI_G2xyz_Px = I_NAI_H3xyz_S+ABX*I_NAI_G2xyz_S;
  Double I_NAI_G2x2z_Px = I_NAI_H3x2z_S+ABX*I_NAI_G2x2z_S;
  Double I_NAI_Gx3y_Px = I_NAI_H2x3y_S+ABX*I_NAI_Gx3y_S;
  Double I_NAI_Gx2yz_Px = I_NAI_H2x2yz_S+ABX*I_NAI_Gx2yz_S;
  Double I_NAI_Gxy2z_Px = I_NAI_H2xy2z_S+ABX*I_NAI_Gxy2z_S;
  Double I_NAI_Gx3z_Px = I_NAI_H2x3z_S+ABX*I_NAI_Gx3z_S;
  Double I_NAI_G4y_Px = I_NAI_Hx4y_S+ABX*I_NAI_G4y_S;
  Double I_NAI_G3yz_Px = I_NAI_Hx3yz_S+ABX*I_NAI_G3yz_S;
  Double I_NAI_G2y2z_Px = I_NAI_Hx2y2z_S+ABX*I_NAI_G2y2z_S;
  Double I_NAI_Gy3z_Px = I_NAI_Hxy3z_S+ABX*I_NAI_Gy3z_S;
  Double I_NAI_G4z_Px = I_NAI_Hx4z_S+ABX*I_NAI_G4z_S;
  Double I_NAI_G3xy_Py = I_NAI_H3x2y_S+ABY*I_NAI_G3xy_S;
  Double I_NAI_G3xz_Py = I_NAI_H3xyz_S+ABY*I_NAI_G3xz_S;
  Double I_NAI_G2x2y_Py = I_NAI_H2x3y_S+ABY*I_NAI_G2x2y_S;
  Double I_NAI_G2xyz_Py = I_NAI_H2x2yz_S+ABY*I_NAI_G2xyz_S;
  Double I_NAI_G2x2z_Py = I_NAI_H2xy2z_S+ABY*I_NAI_G2x2z_S;
  Double I_NAI_Gx3y_Py = I_NAI_Hx4y_S+ABY*I_NAI_Gx3y_S;
  Double I_NAI_Gx2yz_Py = I_NAI_Hx3yz_S+ABY*I_NAI_Gx2yz_S;
  Double I_NAI_Gxy2z_Py = I_NAI_Hx2y2z_S+ABY*I_NAI_Gxy2z_S;
  Double I_NAI_Gx3z_Py = I_NAI_Hxy3z_S+ABY*I_NAI_Gx3z_S;
  Double I_NAI_G4y_Py = I_NAI_H5y_S+ABY*I_NAI_G4y_S;
  Double I_NAI_G3yz_Py = I_NAI_H4yz_S+ABY*I_NAI_G3yz_S;
  Double I_NAI_G2y2z_Py = I_NAI_H3y2z_S+ABY*I_NAI_G2y2z_S;
  Double I_NAI_Gy3z_Py = I_NAI_H2y3z_S+ABY*I_NAI_Gy3z_S;
  Double I_NAI_G4z_Py = I_NAI_Hy4z_S+ABY*I_NAI_G4z_S;
  Double I_NAI_G3xz_Pz = I_NAI_H3x2z_S+ABZ*I_NAI_G3xz_S;
  Double I_NAI_G2xyz_Pz = I_NAI_H2xy2z_S+ABZ*I_NAI_G2xyz_S;
  Double I_NAI_G2x2z_Pz = I_NAI_H2x3z_S+ABZ*I_NAI_G2x2z_S;
  Double I_NAI_Gx2yz_Pz = I_NAI_Hx2y2z_S+ABZ*I_NAI_Gx2yz_S;
  Double I_NAI_Gxy2z_Pz = I_NAI_Hxy3z_S+ABZ*I_NAI_Gxy2z_S;
  Double I_NAI_Gx3z_Pz = I_NAI_Hx4z_S+ABZ*I_NAI_Gx3z_S;
  Double I_NAI_G3yz_Pz = I_NAI_H3y2z_S+ABZ*I_NAI_G3yz_S;
  Double I_NAI_G2y2z_Pz = I_NAI_H2y3z_S+ABZ*I_NAI_G2y2z_S;
  Double I_NAI_Gy3z_Pz = I_NAI_Hy4z_S+ABZ*I_NAI_Gy3z_S;
  Double I_NAI_G4z_Pz = I_NAI_H5z_S+ABZ*I_NAI_G4z_S;

  /************************************************************
   * shell quartet name: SQ_NAI_F_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P
   * RHS shell quartet name: SQ_NAI_F_P
   ************************************************************/
  Double I_NAI_F3x_D2x = I_NAI_G4x_Px+ABX*I_NAI_F3x_Px;
  Double I_NAI_F2xy_D2x = I_NAI_G3xy_Px+ABX*I_NAI_F2xy_Px;
  Double I_NAI_F2xz_D2x = I_NAI_G3xz_Px+ABX*I_NAI_F2xz_Px;
  Double I_NAI_Fx2y_D2x = I_NAI_G2x2y_Px+ABX*I_NAI_Fx2y_Px;
  Double I_NAI_Fxyz_D2x = I_NAI_G2xyz_Px+ABX*I_NAI_Fxyz_Px;
  Double I_NAI_Fx2z_D2x = I_NAI_G2x2z_Px+ABX*I_NAI_Fx2z_Px;
  Double I_NAI_F3y_D2x = I_NAI_Gx3y_Px+ABX*I_NAI_F3y_Px;
  Double I_NAI_F2yz_D2x = I_NAI_Gx2yz_Px+ABX*I_NAI_F2yz_Px;
  Double I_NAI_Fy2z_D2x = I_NAI_Gxy2z_Px+ABX*I_NAI_Fy2z_Px;
  Double I_NAI_F3z_D2x = I_NAI_Gx3z_Px+ABX*I_NAI_F3z_Px;
  Double I_NAI_F3x_Dxy = I_NAI_G3xy_Px+ABY*I_NAI_F3x_Px;
  Double I_NAI_F2xy_Dxy = I_NAI_G2x2y_Px+ABY*I_NAI_F2xy_Px;
  Double I_NAI_F2xz_Dxy = I_NAI_G2xyz_Px+ABY*I_NAI_F2xz_Px;
  Double I_NAI_Fx2y_Dxy = I_NAI_Gx3y_Px+ABY*I_NAI_Fx2y_Px;
  Double I_NAI_Fxyz_Dxy = I_NAI_Gx2yz_Px+ABY*I_NAI_Fxyz_Px;
  Double I_NAI_Fx2z_Dxy = I_NAI_Gxy2z_Px+ABY*I_NAI_Fx2z_Px;
  Double I_NAI_F3y_Dxy = I_NAI_G4y_Px+ABY*I_NAI_F3y_Px;
  Double I_NAI_F2yz_Dxy = I_NAI_G3yz_Px+ABY*I_NAI_F2yz_Px;
  Double I_NAI_Fy2z_Dxy = I_NAI_G2y2z_Px+ABY*I_NAI_Fy2z_Px;
  Double I_NAI_F3z_Dxy = I_NAI_Gy3z_Px+ABY*I_NAI_F3z_Px;
  Double I_NAI_F3x_Dxz = I_NAI_G3xz_Px+ABZ*I_NAI_F3x_Px;
  Double I_NAI_F2xy_Dxz = I_NAI_G2xyz_Px+ABZ*I_NAI_F2xy_Px;
  Double I_NAI_F2xz_Dxz = I_NAI_G2x2z_Px+ABZ*I_NAI_F2xz_Px;
  Double I_NAI_Fx2y_Dxz = I_NAI_Gx2yz_Px+ABZ*I_NAI_Fx2y_Px;
  Double I_NAI_Fxyz_Dxz = I_NAI_Gxy2z_Px+ABZ*I_NAI_Fxyz_Px;
  Double I_NAI_Fx2z_Dxz = I_NAI_Gx3z_Px+ABZ*I_NAI_Fx2z_Px;
  Double I_NAI_F3y_Dxz = I_NAI_G3yz_Px+ABZ*I_NAI_F3y_Px;
  Double I_NAI_F2yz_Dxz = I_NAI_G2y2z_Px+ABZ*I_NAI_F2yz_Px;
  Double I_NAI_Fy2z_Dxz = I_NAI_Gy3z_Px+ABZ*I_NAI_Fy2z_Px;
  Double I_NAI_F3z_Dxz = I_NAI_G4z_Px+ABZ*I_NAI_F3z_Px;
  Double I_NAI_F3x_D2y = I_NAI_G3xy_Py+ABY*I_NAI_F3x_Py;
  Double I_NAI_F2xy_D2y = I_NAI_G2x2y_Py+ABY*I_NAI_F2xy_Py;
  Double I_NAI_F2xz_D2y = I_NAI_G2xyz_Py+ABY*I_NAI_F2xz_Py;
  Double I_NAI_Fx2y_D2y = I_NAI_Gx3y_Py+ABY*I_NAI_Fx2y_Py;
  Double I_NAI_Fxyz_D2y = I_NAI_Gx2yz_Py+ABY*I_NAI_Fxyz_Py;
  Double I_NAI_Fx2z_D2y = I_NAI_Gxy2z_Py+ABY*I_NAI_Fx2z_Py;
  Double I_NAI_F3y_D2y = I_NAI_G4y_Py+ABY*I_NAI_F3y_Py;
  Double I_NAI_F2yz_D2y = I_NAI_G3yz_Py+ABY*I_NAI_F2yz_Py;
  Double I_NAI_Fy2z_D2y = I_NAI_G2y2z_Py+ABY*I_NAI_Fy2z_Py;
  Double I_NAI_F3z_D2y = I_NAI_Gy3z_Py+ABY*I_NAI_F3z_Py;
  Double I_NAI_F3x_Dyz = I_NAI_G3xz_Py+ABZ*I_NAI_F3x_Py;
  Double I_NAI_F2xy_Dyz = I_NAI_G2xyz_Py+ABZ*I_NAI_F2xy_Py;
  Double I_NAI_F2xz_Dyz = I_NAI_G2x2z_Py+ABZ*I_NAI_F2xz_Py;
  Double I_NAI_Fx2y_Dyz = I_NAI_Gx2yz_Py+ABZ*I_NAI_Fx2y_Py;
  Double I_NAI_Fxyz_Dyz = I_NAI_Gxy2z_Py+ABZ*I_NAI_Fxyz_Py;
  Double I_NAI_Fx2z_Dyz = I_NAI_Gx3z_Py+ABZ*I_NAI_Fx2z_Py;
  Double I_NAI_F3y_Dyz = I_NAI_G3yz_Py+ABZ*I_NAI_F3y_Py;
  Double I_NAI_F2yz_Dyz = I_NAI_G2y2z_Py+ABZ*I_NAI_F2yz_Py;
  Double I_NAI_Fy2z_Dyz = I_NAI_Gy3z_Py+ABZ*I_NAI_Fy2z_Py;
  Double I_NAI_F3z_Dyz = I_NAI_G4z_Py+ABZ*I_NAI_F3z_Py;
  Double I_NAI_F3x_D2z = I_NAI_G3xz_Pz+ABZ*I_NAI_F3x_Pz;
  Double I_NAI_F2xy_D2z = I_NAI_G2xyz_Pz+ABZ*I_NAI_F2xy_Pz;
  Double I_NAI_F2xz_D2z = I_NAI_G2x2z_Pz+ABZ*I_NAI_F2xz_Pz;
  Double I_NAI_Fx2y_D2z = I_NAI_Gx2yz_Pz+ABZ*I_NAI_Fx2y_Pz;
  Double I_NAI_Fxyz_D2z = I_NAI_Gxy2z_Pz+ABZ*I_NAI_Fxyz_Pz;
  Double I_NAI_Fx2z_D2z = I_NAI_Gx3z_Pz+ABZ*I_NAI_Fx2z_Pz;
  Double I_NAI_F3y_D2z = I_NAI_G3yz_Pz+ABZ*I_NAI_F3y_Pz;
  Double I_NAI_F2yz_D2z = I_NAI_G2y2z_Pz+ABZ*I_NAI_F2yz_Pz;
  Double I_NAI_Fy2z_D2z = I_NAI_Gy3z_Pz+ABZ*I_NAI_Fy2z_Pz;
  Double I_NAI_F3z_D2z = I_NAI_G4z_Pz+ABZ*I_NAI_F3z_Pz;

  /************************************************************
   * shell quartet name: SQ_NAI_D_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_D
   * RHS shell quartet name: SQ_NAI_D_D
   ************************************************************/
  Double I_NAI_D2x_F3x = I_NAI_F3x_D2x+ABX*I_NAI_D2x_D2x;
  Double I_NAI_Dxy_F3x = I_NAI_F2xy_D2x+ABX*I_NAI_Dxy_D2x;
  Double I_NAI_Dxz_F3x = I_NAI_F2xz_D2x+ABX*I_NAI_Dxz_D2x;
  Double I_NAI_D2y_F3x = I_NAI_Fx2y_D2x+ABX*I_NAI_D2y_D2x;
  Double I_NAI_Dyz_F3x = I_NAI_Fxyz_D2x+ABX*I_NAI_Dyz_D2x;
  Double I_NAI_D2z_F3x = I_NAI_Fx2z_D2x+ABX*I_NAI_D2z_D2x;
  Double I_NAI_D2x_F2xy = I_NAI_F2xy_D2x+ABY*I_NAI_D2x_D2x;
  Double I_NAI_Dxy_F2xy = I_NAI_Fx2y_D2x+ABY*I_NAI_Dxy_D2x;
  Double I_NAI_Dxz_F2xy = I_NAI_Fxyz_D2x+ABY*I_NAI_Dxz_D2x;
  Double I_NAI_D2y_F2xy = I_NAI_F3y_D2x+ABY*I_NAI_D2y_D2x;
  Double I_NAI_Dyz_F2xy = I_NAI_F2yz_D2x+ABY*I_NAI_Dyz_D2x;
  Double I_NAI_D2z_F2xy = I_NAI_Fy2z_D2x+ABY*I_NAI_D2z_D2x;
  Double I_NAI_D2x_F2xz = I_NAI_F2xz_D2x+ABZ*I_NAI_D2x_D2x;
  Double I_NAI_Dxy_F2xz = I_NAI_Fxyz_D2x+ABZ*I_NAI_Dxy_D2x;
  Double I_NAI_Dxz_F2xz = I_NAI_Fx2z_D2x+ABZ*I_NAI_Dxz_D2x;
  Double I_NAI_D2y_F2xz = I_NAI_F2yz_D2x+ABZ*I_NAI_D2y_D2x;
  Double I_NAI_Dyz_F2xz = I_NAI_Fy2z_D2x+ABZ*I_NAI_Dyz_D2x;
  Double I_NAI_D2z_F2xz = I_NAI_F3z_D2x+ABZ*I_NAI_D2z_D2x;
  Double I_NAI_D2x_Fx2y = I_NAI_F3x_D2y+ABX*I_NAI_D2x_D2y;
  Double I_NAI_Dxy_Fx2y = I_NAI_F2xy_D2y+ABX*I_NAI_Dxy_D2y;
  Double I_NAI_Dxz_Fx2y = I_NAI_F2xz_D2y+ABX*I_NAI_Dxz_D2y;
  Double I_NAI_D2y_Fx2y = I_NAI_Fx2y_D2y+ABX*I_NAI_D2y_D2y;
  Double I_NAI_Dyz_Fx2y = I_NAI_Fxyz_D2y+ABX*I_NAI_Dyz_D2y;
  Double I_NAI_D2z_Fx2y = I_NAI_Fx2z_D2y+ABX*I_NAI_D2z_D2y;
  Double I_NAI_D2x_Fxyz = I_NAI_F2xz_Dxy+ABZ*I_NAI_D2x_Dxy;
  Double I_NAI_Dxy_Fxyz = I_NAI_Fxyz_Dxy+ABZ*I_NAI_Dxy_Dxy;
  Double I_NAI_Dxz_Fxyz = I_NAI_Fx2z_Dxy+ABZ*I_NAI_Dxz_Dxy;
  Double I_NAI_D2y_Fxyz = I_NAI_F2yz_Dxy+ABZ*I_NAI_D2y_Dxy;
  Double I_NAI_Dyz_Fxyz = I_NAI_Fy2z_Dxy+ABZ*I_NAI_Dyz_Dxy;
  Double I_NAI_D2z_Fxyz = I_NAI_F3z_Dxy+ABZ*I_NAI_D2z_Dxy;
  Double I_NAI_D2x_Fx2z = I_NAI_F3x_D2z+ABX*I_NAI_D2x_D2z;
  Double I_NAI_Dxy_Fx2z = I_NAI_F2xy_D2z+ABX*I_NAI_Dxy_D2z;
  Double I_NAI_Dxz_Fx2z = I_NAI_F2xz_D2z+ABX*I_NAI_Dxz_D2z;
  Double I_NAI_D2y_Fx2z = I_NAI_Fx2y_D2z+ABX*I_NAI_D2y_D2z;
  Double I_NAI_Dyz_Fx2z = I_NAI_Fxyz_D2z+ABX*I_NAI_Dyz_D2z;
  Double I_NAI_D2z_Fx2z = I_NAI_Fx2z_D2z+ABX*I_NAI_D2z_D2z;
  Double I_NAI_D2x_F3y = I_NAI_F2xy_D2y+ABY*I_NAI_D2x_D2y;
  Double I_NAI_Dxy_F3y = I_NAI_Fx2y_D2y+ABY*I_NAI_Dxy_D2y;
  Double I_NAI_Dxz_F3y = I_NAI_Fxyz_D2y+ABY*I_NAI_Dxz_D2y;
  Double I_NAI_D2y_F3y = I_NAI_F3y_D2y+ABY*I_NAI_D2y_D2y;
  Double I_NAI_Dyz_F3y = I_NAI_F2yz_D2y+ABY*I_NAI_Dyz_D2y;
  Double I_NAI_D2z_F3y = I_NAI_Fy2z_D2y+ABY*I_NAI_D2z_D2y;
  Double I_NAI_D2x_F2yz = I_NAI_F2xz_D2y+ABZ*I_NAI_D2x_D2y;
  Double I_NAI_Dxy_F2yz = I_NAI_Fxyz_D2y+ABZ*I_NAI_Dxy_D2y;
  Double I_NAI_Dxz_F2yz = I_NAI_Fx2z_D2y+ABZ*I_NAI_Dxz_D2y;
  Double I_NAI_D2y_F2yz = I_NAI_F2yz_D2y+ABZ*I_NAI_D2y_D2y;
  Double I_NAI_Dyz_F2yz = I_NAI_Fy2z_D2y+ABZ*I_NAI_Dyz_D2y;
  Double I_NAI_D2z_F2yz = I_NAI_F3z_D2y+ABZ*I_NAI_D2z_D2y;
  Double I_NAI_D2x_Fy2z = I_NAI_F2xy_D2z+ABY*I_NAI_D2x_D2z;
  Double I_NAI_Dxy_Fy2z = I_NAI_Fx2y_D2z+ABY*I_NAI_Dxy_D2z;
  Double I_NAI_Dxz_Fy2z = I_NAI_Fxyz_D2z+ABY*I_NAI_Dxz_D2z;
  Double I_NAI_D2y_Fy2z = I_NAI_F3y_D2z+ABY*I_NAI_D2y_D2z;
  Double I_NAI_Dyz_Fy2z = I_NAI_F2yz_D2z+ABY*I_NAI_Dyz_D2z;
  Double I_NAI_D2z_Fy2z = I_NAI_Fy2z_D2z+ABY*I_NAI_D2z_D2z;
  Double I_NAI_D2x_F3z = I_NAI_F2xz_D2z+ABZ*I_NAI_D2x_D2z;
  Double I_NAI_Dxy_F3z = I_NAI_Fxyz_D2z+ABZ*I_NAI_Dxy_D2z;
  Double I_NAI_Dxz_F3z = I_NAI_Fx2z_D2z+ABZ*I_NAI_Dxz_D2z;
  Double I_NAI_D2y_F3z = I_NAI_F2yz_D2z+ABZ*I_NAI_D2y_D2z;
  Double I_NAI_Dyz_F3z = I_NAI_Fy2z_D2z+ABZ*I_NAI_Dyz_D2z;
  Double I_NAI_D2z_F3z = I_NAI_F3z_D2z+ABZ*I_NAI_D2z_D2z;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_a
   * RHS shell quartet name: SQ_NAI_G_S_a
   ************************************************************/
  Double I_NAI_G4x_Px_a = I_NAI_H5x_S_a+ABX*I_NAI_G4x_S_a;
  Double I_NAI_G3xy_Px_a = I_NAI_H4xy_S_a+ABX*I_NAI_G3xy_S_a;
  Double I_NAI_G3xz_Px_a = I_NAI_H4xz_S_a+ABX*I_NAI_G3xz_S_a;
  Double I_NAI_G2x2y_Px_a = I_NAI_H3x2y_S_a+ABX*I_NAI_G2x2y_S_a;
  Double I_NAI_G2xyz_Px_a = I_NAI_H3xyz_S_a+ABX*I_NAI_G2xyz_S_a;
  Double I_NAI_G2x2z_Px_a = I_NAI_H3x2z_S_a+ABX*I_NAI_G2x2z_S_a;
  Double I_NAI_Gx3y_Px_a = I_NAI_H2x3y_S_a+ABX*I_NAI_Gx3y_S_a;
  Double I_NAI_Gx2yz_Px_a = I_NAI_H2x2yz_S_a+ABX*I_NAI_Gx2yz_S_a;
  Double I_NAI_Gxy2z_Px_a = I_NAI_H2xy2z_S_a+ABX*I_NAI_Gxy2z_S_a;
  Double I_NAI_Gx3z_Px_a = I_NAI_H2x3z_S_a+ABX*I_NAI_Gx3z_S_a;
  Double I_NAI_G4y_Px_a = I_NAI_Hx4y_S_a+ABX*I_NAI_G4y_S_a;
  Double I_NAI_G3yz_Px_a = I_NAI_Hx3yz_S_a+ABX*I_NAI_G3yz_S_a;
  Double I_NAI_G2y2z_Px_a = I_NAI_Hx2y2z_S_a+ABX*I_NAI_G2y2z_S_a;
  Double I_NAI_Gy3z_Px_a = I_NAI_Hxy3z_S_a+ABX*I_NAI_Gy3z_S_a;
  Double I_NAI_G4z_Px_a = I_NAI_Hx4z_S_a+ABX*I_NAI_G4z_S_a;
  Double I_NAI_G4x_Py_a = I_NAI_H4xy_S_a+ABY*I_NAI_G4x_S_a;
  Double I_NAI_G3xy_Py_a = I_NAI_H3x2y_S_a+ABY*I_NAI_G3xy_S_a;
  Double I_NAI_G3xz_Py_a = I_NAI_H3xyz_S_a+ABY*I_NAI_G3xz_S_a;
  Double I_NAI_G2x2y_Py_a = I_NAI_H2x3y_S_a+ABY*I_NAI_G2x2y_S_a;
  Double I_NAI_G2xyz_Py_a = I_NAI_H2x2yz_S_a+ABY*I_NAI_G2xyz_S_a;
  Double I_NAI_G2x2z_Py_a = I_NAI_H2xy2z_S_a+ABY*I_NAI_G2x2z_S_a;
  Double I_NAI_Gx3y_Py_a = I_NAI_Hx4y_S_a+ABY*I_NAI_Gx3y_S_a;
  Double I_NAI_Gx2yz_Py_a = I_NAI_Hx3yz_S_a+ABY*I_NAI_Gx2yz_S_a;
  Double I_NAI_Gxy2z_Py_a = I_NAI_Hx2y2z_S_a+ABY*I_NAI_Gxy2z_S_a;
  Double I_NAI_Gx3z_Py_a = I_NAI_Hxy3z_S_a+ABY*I_NAI_Gx3z_S_a;
  Double I_NAI_G4y_Py_a = I_NAI_H5y_S_a+ABY*I_NAI_G4y_S_a;
  Double I_NAI_G3yz_Py_a = I_NAI_H4yz_S_a+ABY*I_NAI_G3yz_S_a;
  Double I_NAI_G2y2z_Py_a = I_NAI_H3y2z_S_a+ABY*I_NAI_G2y2z_S_a;
  Double I_NAI_Gy3z_Py_a = I_NAI_H2y3z_S_a+ABY*I_NAI_Gy3z_S_a;
  Double I_NAI_G4z_Py_a = I_NAI_Hy4z_S_a+ABY*I_NAI_G4z_S_a;
  Double I_NAI_G4x_Pz_a = I_NAI_H4xz_S_a+ABZ*I_NAI_G4x_S_a;
  Double I_NAI_G3xy_Pz_a = I_NAI_H3xyz_S_a+ABZ*I_NAI_G3xy_S_a;
  Double I_NAI_G3xz_Pz_a = I_NAI_H3x2z_S_a+ABZ*I_NAI_G3xz_S_a;
  Double I_NAI_G2x2y_Pz_a = I_NAI_H2x2yz_S_a+ABZ*I_NAI_G2x2y_S_a;
  Double I_NAI_G2xyz_Pz_a = I_NAI_H2xy2z_S_a+ABZ*I_NAI_G2xyz_S_a;
  Double I_NAI_G2x2z_Pz_a = I_NAI_H2x3z_S_a+ABZ*I_NAI_G2x2z_S_a;
  Double I_NAI_Gx3y_Pz_a = I_NAI_Hx3yz_S_a+ABZ*I_NAI_Gx3y_S_a;
  Double I_NAI_Gx2yz_Pz_a = I_NAI_Hx2y2z_S_a+ABZ*I_NAI_Gx2yz_S_a;
  Double I_NAI_Gxy2z_Pz_a = I_NAI_Hxy3z_S_a+ABZ*I_NAI_Gxy2z_S_a;
  Double I_NAI_Gx3z_Pz_a = I_NAI_Hx4z_S_a+ABZ*I_NAI_Gx3z_S_a;
  Double I_NAI_G4y_Pz_a = I_NAI_H4yz_S_a+ABZ*I_NAI_G4y_S_a;
  Double I_NAI_G3yz_Pz_a = I_NAI_H3y2z_S_a+ABZ*I_NAI_G3yz_S_a;
  Double I_NAI_G2y2z_Pz_a = I_NAI_H2y3z_S_a+ABZ*I_NAI_G2y2z_S_a;
  Double I_NAI_Gy3z_Pz_a = I_NAI_Hy4z_S_a+ABZ*I_NAI_Gy3z_S_a;
  Double I_NAI_G4z_Pz_a = I_NAI_H5z_S_a+ABZ*I_NAI_G4z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_a
   * RHS shell quartet name: SQ_NAI_H_S_a
   ************************************************************/
  Double I_NAI_H5x_Px_a = I_NAI_I6x_S_a+ABX*I_NAI_H5x_S_a;
  Double I_NAI_H4xy_Px_a = I_NAI_I5xy_S_a+ABX*I_NAI_H4xy_S_a;
  Double I_NAI_H4xz_Px_a = I_NAI_I5xz_S_a+ABX*I_NAI_H4xz_S_a;
  Double I_NAI_H3x2y_Px_a = I_NAI_I4x2y_S_a+ABX*I_NAI_H3x2y_S_a;
  Double I_NAI_H3xyz_Px_a = I_NAI_I4xyz_S_a+ABX*I_NAI_H3xyz_S_a;
  Double I_NAI_H3x2z_Px_a = I_NAI_I4x2z_S_a+ABX*I_NAI_H3x2z_S_a;
  Double I_NAI_H2x3y_Px_a = I_NAI_I3x3y_S_a+ABX*I_NAI_H2x3y_S_a;
  Double I_NAI_H2x2yz_Px_a = I_NAI_I3x2yz_S_a+ABX*I_NAI_H2x2yz_S_a;
  Double I_NAI_H2xy2z_Px_a = I_NAI_I3xy2z_S_a+ABX*I_NAI_H2xy2z_S_a;
  Double I_NAI_H2x3z_Px_a = I_NAI_I3x3z_S_a+ABX*I_NAI_H2x3z_S_a;
  Double I_NAI_Hx4y_Px_a = I_NAI_I2x4y_S_a+ABX*I_NAI_Hx4y_S_a;
  Double I_NAI_Hx3yz_Px_a = I_NAI_I2x3yz_S_a+ABX*I_NAI_Hx3yz_S_a;
  Double I_NAI_Hx2y2z_Px_a = I_NAI_I2x2y2z_S_a+ABX*I_NAI_Hx2y2z_S_a;
  Double I_NAI_Hxy3z_Px_a = I_NAI_I2xy3z_S_a+ABX*I_NAI_Hxy3z_S_a;
  Double I_NAI_Hx4z_Px_a = I_NAI_I2x4z_S_a+ABX*I_NAI_Hx4z_S_a;
  Double I_NAI_H5y_Px_a = I_NAI_Ix5y_S_a+ABX*I_NAI_H5y_S_a;
  Double I_NAI_H4yz_Px_a = I_NAI_Ix4yz_S_a+ABX*I_NAI_H4yz_S_a;
  Double I_NAI_H3y2z_Px_a = I_NAI_Ix3y2z_S_a+ABX*I_NAI_H3y2z_S_a;
  Double I_NAI_H2y3z_Px_a = I_NAI_Ix2y3z_S_a+ABX*I_NAI_H2y3z_S_a;
  Double I_NAI_Hy4z_Px_a = I_NAI_Ixy4z_S_a+ABX*I_NAI_Hy4z_S_a;
  Double I_NAI_H5z_Px_a = I_NAI_Ix5z_S_a+ABX*I_NAI_H5z_S_a;
  Double I_NAI_H5x_Py_a = I_NAI_I5xy_S_a+ABY*I_NAI_H5x_S_a;
  Double I_NAI_H4xy_Py_a = I_NAI_I4x2y_S_a+ABY*I_NAI_H4xy_S_a;
  Double I_NAI_H4xz_Py_a = I_NAI_I4xyz_S_a+ABY*I_NAI_H4xz_S_a;
  Double I_NAI_H3x2y_Py_a = I_NAI_I3x3y_S_a+ABY*I_NAI_H3x2y_S_a;
  Double I_NAI_H3xyz_Py_a = I_NAI_I3x2yz_S_a+ABY*I_NAI_H3xyz_S_a;
  Double I_NAI_H3x2z_Py_a = I_NAI_I3xy2z_S_a+ABY*I_NAI_H3x2z_S_a;
  Double I_NAI_H2x3y_Py_a = I_NAI_I2x4y_S_a+ABY*I_NAI_H2x3y_S_a;
  Double I_NAI_H2x2yz_Py_a = I_NAI_I2x3yz_S_a+ABY*I_NAI_H2x2yz_S_a;
  Double I_NAI_H2xy2z_Py_a = I_NAI_I2x2y2z_S_a+ABY*I_NAI_H2xy2z_S_a;
  Double I_NAI_H2x3z_Py_a = I_NAI_I2xy3z_S_a+ABY*I_NAI_H2x3z_S_a;
  Double I_NAI_Hx4y_Py_a = I_NAI_Ix5y_S_a+ABY*I_NAI_Hx4y_S_a;
  Double I_NAI_Hx3yz_Py_a = I_NAI_Ix4yz_S_a+ABY*I_NAI_Hx3yz_S_a;
  Double I_NAI_Hx2y2z_Py_a = I_NAI_Ix3y2z_S_a+ABY*I_NAI_Hx2y2z_S_a;
  Double I_NAI_Hxy3z_Py_a = I_NAI_Ix2y3z_S_a+ABY*I_NAI_Hxy3z_S_a;
  Double I_NAI_Hx4z_Py_a = I_NAI_Ixy4z_S_a+ABY*I_NAI_Hx4z_S_a;
  Double I_NAI_H5y_Py_a = I_NAI_I6y_S_a+ABY*I_NAI_H5y_S_a;
  Double I_NAI_H4yz_Py_a = I_NAI_I5yz_S_a+ABY*I_NAI_H4yz_S_a;
  Double I_NAI_H3y2z_Py_a = I_NAI_I4y2z_S_a+ABY*I_NAI_H3y2z_S_a;
  Double I_NAI_H2y3z_Py_a = I_NAI_I3y3z_S_a+ABY*I_NAI_H2y3z_S_a;
  Double I_NAI_Hy4z_Py_a = I_NAI_I2y4z_S_a+ABY*I_NAI_Hy4z_S_a;
  Double I_NAI_H5z_Py_a = I_NAI_Iy5z_S_a+ABY*I_NAI_H5z_S_a;
  Double I_NAI_H5x_Pz_a = I_NAI_I5xz_S_a+ABZ*I_NAI_H5x_S_a;
  Double I_NAI_H4xy_Pz_a = I_NAI_I4xyz_S_a+ABZ*I_NAI_H4xy_S_a;
  Double I_NAI_H4xz_Pz_a = I_NAI_I4x2z_S_a+ABZ*I_NAI_H4xz_S_a;
  Double I_NAI_H3x2y_Pz_a = I_NAI_I3x2yz_S_a+ABZ*I_NAI_H3x2y_S_a;
  Double I_NAI_H3xyz_Pz_a = I_NAI_I3xy2z_S_a+ABZ*I_NAI_H3xyz_S_a;
  Double I_NAI_H3x2z_Pz_a = I_NAI_I3x3z_S_a+ABZ*I_NAI_H3x2z_S_a;
  Double I_NAI_H2x3y_Pz_a = I_NAI_I2x3yz_S_a+ABZ*I_NAI_H2x3y_S_a;
  Double I_NAI_H2x2yz_Pz_a = I_NAI_I2x2y2z_S_a+ABZ*I_NAI_H2x2yz_S_a;
  Double I_NAI_H2xy2z_Pz_a = I_NAI_I2xy3z_S_a+ABZ*I_NAI_H2xy2z_S_a;
  Double I_NAI_H2x3z_Pz_a = I_NAI_I2x4z_S_a+ABZ*I_NAI_H2x3z_S_a;
  Double I_NAI_Hx4y_Pz_a = I_NAI_Ix4yz_S_a+ABZ*I_NAI_Hx4y_S_a;
  Double I_NAI_Hx3yz_Pz_a = I_NAI_Ix3y2z_S_a+ABZ*I_NAI_Hx3yz_S_a;
  Double I_NAI_Hx2y2z_Pz_a = I_NAI_Ix2y3z_S_a+ABZ*I_NAI_Hx2y2z_S_a;
  Double I_NAI_Hxy3z_Pz_a = I_NAI_Ixy4z_S_a+ABZ*I_NAI_Hxy3z_S_a;
  Double I_NAI_Hx4z_Pz_a = I_NAI_Ix5z_S_a+ABZ*I_NAI_Hx4z_S_a;
  Double I_NAI_H5y_Pz_a = I_NAI_I5yz_S_a+ABZ*I_NAI_H5y_S_a;
  Double I_NAI_H4yz_Pz_a = I_NAI_I4y2z_S_a+ABZ*I_NAI_H4yz_S_a;
  Double I_NAI_H3y2z_Pz_a = I_NAI_I3y3z_S_a+ABZ*I_NAI_H3y2z_S_a;
  Double I_NAI_H2y3z_Pz_a = I_NAI_I2y4z_S_a+ABZ*I_NAI_H2y3z_S_a;
  Double I_NAI_Hy4z_Pz_a = I_NAI_Iy5z_S_a+ABZ*I_NAI_Hy4z_S_a;
  Double I_NAI_H5z_Pz_a = I_NAI_I6z_S_a+ABZ*I_NAI_H5z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_a
   * RHS shell quartet name: SQ_NAI_G_P_a
   ************************************************************/
  Double I_NAI_G4x_D2x_a = I_NAI_H5x_Px_a+ABX*I_NAI_G4x_Px_a;
  Double I_NAI_G3xy_D2x_a = I_NAI_H4xy_Px_a+ABX*I_NAI_G3xy_Px_a;
  Double I_NAI_G3xz_D2x_a = I_NAI_H4xz_Px_a+ABX*I_NAI_G3xz_Px_a;
  Double I_NAI_G2x2y_D2x_a = I_NAI_H3x2y_Px_a+ABX*I_NAI_G2x2y_Px_a;
  Double I_NAI_G2xyz_D2x_a = I_NAI_H3xyz_Px_a+ABX*I_NAI_G2xyz_Px_a;
  Double I_NAI_G2x2z_D2x_a = I_NAI_H3x2z_Px_a+ABX*I_NAI_G2x2z_Px_a;
  Double I_NAI_Gx3y_D2x_a = I_NAI_H2x3y_Px_a+ABX*I_NAI_Gx3y_Px_a;
  Double I_NAI_Gx2yz_D2x_a = I_NAI_H2x2yz_Px_a+ABX*I_NAI_Gx2yz_Px_a;
  Double I_NAI_Gxy2z_D2x_a = I_NAI_H2xy2z_Px_a+ABX*I_NAI_Gxy2z_Px_a;
  Double I_NAI_Gx3z_D2x_a = I_NAI_H2x3z_Px_a+ABX*I_NAI_Gx3z_Px_a;
  Double I_NAI_G4y_D2x_a = I_NAI_Hx4y_Px_a+ABX*I_NAI_G4y_Px_a;
  Double I_NAI_G3yz_D2x_a = I_NAI_Hx3yz_Px_a+ABX*I_NAI_G3yz_Px_a;
  Double I_NAI_G2y2z_D2x_a = I_NAI_Hx2y2z_Px_a+ABX*I_NAI_G2y2z_Px_a;
  Double I_NAI_Gy3z_D2x_a = I_NAI_Hxy3z_Px_a+ABX*I_NAI_Gy3z_Px_a;
  Double I_NAI_G4z_D2x_a = I_NAI_Hx4z_Px_a+ABX*I_NAI_G4z_Px_a;
  Double I_NAI_G4x_Dxy_a = I_NAI_H4xy_Px_a+ABY*I_NAI_G4x_Px_a;
  Double I_NAI_G3xy_Dxy_a = I_NAI_H3x2y_Px_a+ABY*I_NAI_G3xy_Px_a;
  Double I_NAI_G3xz_Dxy_a = I_NAI_H3xyz_Px_a+ABY*I_NAI_G3xz_Px_a;
  Double I_NAI_G2x2y_Dxy_a = I_NAI_H2x3y_Px_a+ABY*I_NAI_G2x2y_Px_a;
  Double I_NAI_G2xyz_Dxy_a = I_NAI_H2x2yz_Px_a+ABY*I_NAI_G2xyz_Px_a;
  Double I_NAI_G2x2z_Dxy_a = I_NAI_H2xy2z_Px_a+ABY*I_NAI_G2x2z_Px_a;
  Double I_NAI_Gx3y_Dxy_a = I_NAI_Hx4y_Px_a+ABY*I_NAI_Gx3y_Px_a;
  Double I_NAI_Gx2yz_Dxy_a = I_NAI_Hx3yz_Px_a+ABY*I_NAI_Gx2yz_Px_a;
  Double I_NAI_Gxy2z_Dxy_a = I_NAI_Hx2y2z_Px_a+ABY*I_NAI_Gxy2z_Px_a;
  Double I_NAI_Gx3z_Dxy_a = I_NAI_Hxy3z_Px_a+ABY*I_NAI_Gx3z_Px_a;
  Double I_NAI_G4y_Dxy_a = I_NAI_H5y_Px_a+ABY*I_NAI_G4y_Px_a;
  Double I_NAI_G3yz_Dxy_a = I_NAI_H4yz_Px_a+ABY*I_NAI_G3yz_Px_a;
  Double I_NAI_G2y2z_Dxy_a = I_NAI_H3y2z_Px_a+ABY*I_NAI_G2y2z_Px_a;
  Double I_NAI_Gy3z_Dxy_a = I_NAI_H2y3z_Px_a+ABY*I_NAI_Gy3z_Px_a;
  Double I_NAI_G4z_Dxy_a = I_NAI_Hy4z_Px_a+ABY*I_NAI_G4z_Px_a;
  Double I_NAI_G4x_D2y_a = I_NAI_H4xy_Py_a+ABY*I_NAI_G4x_Py_a;
  Double I_NAI_G3xy_D2y_a = I_NAI_H3x2y_Py_a+ABY*I_NAI_G3xy_Py_a;
  Double I_NAI_G3xz_D2y_a = I_NAI_H3xyz_Py_a+ABY*I_NAI_G3xz_Py_a;
  Double I_NAI_G2x2y_D2y_a = I_NAI_H2x3y_Py_a+ABY*I_NAI_G2x2y_Py_a;
  Double I_NAI_G2xyz_D2y_a = I_NAI_H2x2yz_Py_a+ABY*I_NAI_G2xyz_Py_a;
  Double I_NAI_G2x2z_D2y_a = I_NAI_H2xy2z_Py_a+ABY*I_NAI_G2x2z_Py_a;
  Double I_NAI_Gx3y_D2y_a = I_NAI_Hx4y_Py_a+ABY*I_NAI_Gx3y_Py_a;
  Double I_NAI_Gx2yz_D2y_a = I_NAI_Hx3yz_Py_a+ABY*I_NAI_Gx2yz_Py_a;
  Double I_NAI_Gxy2z_D2y_a = I_NAI_Hx2y2z_Py_a+ABY*I_NAI_Gxy2z_Py_a;
  Double I_NAI_Gx3z_D2y_a = I_NAI_Hxy3z_Py_a+ABY*I_NAI_Gx3z_Py_a;
  Double I_NAI_G4y_D2y_a = I_NAI_H5y_Py_a+ABY*I_NAI_G4y_Py_a;
  Double I_NAI_G3yz_D2y_a = I_NAI_H4yz_Py_a+ABY*I_NAI_G3yz_Py_a;
  Double I_NAI_G2y2z_D2y_a = I_NAI_H3y2z_Py_a+ABY*I_NAI_G2y2z_Py_a;
  Double I_NAI_Gy3z_D2y_a = I_NAI_H2y3z_Py_a+ABY*I_NAI_Gy3z_Py_a;
  Double I_NAI_G4z_D2y_a = I_NAI_Hy4z_Py_a+ABY*I_NAI_G4z_Py_a;
  Double I_NAI_G4x_D2z_a = I_NAI_H4xz_Pz_a+ABZ*I_NAI_G4x_Pz_a;
  Double I_NAI_G3xy_D2z_a = I_NAI_H3xyz_Pz_a+ABZ*I_NAI_G3xy_Pz_a;
  Double I_NAI_G3xz_D2z_a = I_NAI_H3x2z_Pz_a+ABZ*I_NAI_G3xz_Pz_a;
  Double I_NAI_G2x2y_D2z_a = I_NAI_H2x2yz_Pz_a+ABZ*I_NAI_G2x2y_Pz_a;
  Double I_NAI_G2xyz_D2z_a = I_NAI_H2xy2z_Pz_a+ABZ*I_NAI_G2xyz_Pz_a;
  Double I_NAI_G2x2z_D2z_a = I_NAI_H2x3z_Pz_a+ABZ*I_NAI_G2x2z_Pz_a;
  Double I_NAI_Gx3y_D2z_a = I_NAI_Hx3yz_Pz_a+ABZ*I_NAI_Gx3y_Pz_a;
  Double I_NAI_Gx2yz_D2z_a = I_NAI_Hx2y2z_Pz_a+ABZ*I_NAI_Gx2yz_Pz_a;
  Double I_NAI_Gxy2z_D2z_a = I_NAI_Hxy3z_Pz_a+ABZ*I_NAI_Gxy2z_Pz_a;
  Double I_NAI_Gx3z_D2z_a = I_NAI_Hx4z_Pz_a+ABZ*I_NAI_Gx3z_Pz_a;
  Double I_NAI_G4y_D2z_a = I_NAI_H4yz_Pz_a+ABZ*I_NAI_G4y_Pz_a;
  Double I_NAI_G3yz_D2z_a = I_NAI_H3y2z_Pz_a+ABZ*I_NAI_G3yz_Pz_a;
  Double I_NAI_G2y2z_D2z_a = I_NAI_H2y3z_Pz_a+ABZ*I_NAI_G2y2z_Pz_a;
  Double I_NAI_Gy3z_D2z_a = I_NAI_Hy4z_Pz_a+ABZ*I_NAI_Gy3z_Pz_a;
  Double I_NAI_G4z_D2z_a = I_NAI_H5z_Pz_a+ABZ*I_NAI_G4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 16 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S_a
   * RHS shell quartet name: SQ_NAI_I_S_a
   ************************************************************/
  Double I_NAI_I6x_Px_a = I_NAI_K7x_S_a+ABX*I_NAI_I6x_S_a;
  Double I_NAI_I5xy_Px_a = I_NAI_K6xy_S_a+ABX*I_NAI_I5xy_S_a;
  Double I_NAI_I5xz_Px_a = I_NAI_K6xz_S_a+ABX*I_NAI_I5xz_S_a;
  Double I_NAI_I4x2y_Px_a = I_NAI_K5x2y_S_a+ABX*I_NAI_I4x2y_S_a;
  Double I_NAI_I4xyz_Px_a = I_NAI_K5xyz_S_a+ABX*I_NAI_I4xyz_S_a;
  Double I_NAI_I4x2z_Px_a = I_NAI_K5x2z_S_a+ABX*I_NAI_I4x2z_S_a;
  Double I_NAI_I3x3y_Px_a = I_NAI_K4x3y_S_a+ABX*I_NAI_I3x3y_S_a;
  Double I_NAI_I3x2yz_Px_a = I_NAI_K4x2yz_S_a+ABX*I_NAI_I3x2yz_S_a;
  Double I_NAI_I3xy2z_Px_a = I_NAI_K4xy2z_S_a+ABX*I_NAI_I3xy2z_S_a;
  Double I_NAI_I3x3z_Px_a = I_NAI_K4x3z_S_a+ABX*I_NAI_I3x3z_S_a;
  Double I_NAI_I2x4y_Px_a = I_NAI_K3x4y_S_a+ABX*I_NAI_I2x4y_S_a;
  Double I_NAI_I2x3yz_Px_a = I_NAI_K3x3yz_S_a+ABX*I_NAI_I2x3yz_S_a;
  Double I_NAI_I2x2y2z_Px_a = I_NAI_K3x2y2z_S_a+ABX*I_NAI_I2x2y2z_S_a;
  Double I_NAI_I2xy3z_Px_a = I_NAI_K3xy3z_S_a+ABX*I_NAI_I2xy3z_S_a;
  Double I_NAI_I2x4z_Px_a = I_NAI_K3x4z_S_a+ABX*I_NAI_I2x4z_S_a;
  Double I_NAI_Ix5y_Px_a = I_NAI_K2x5y_S_a+ABX*I_NAI_Ix5y_S_a;
  Double I_NAI_Ix4yz_Px_a = I_NAI_K2x4yz_S_a+ABX*I_NAI_Ix4yz_S_a;
  Double I_NAI_Ix3y2z_Px_a = I_NAI_K2x3y2z_S_a+ABX*I_NAI_Ix3y2z_S_a;
  Double I_NAI_Ix2y3z_Px_a = I_NAI_K2x2y3z_S_a+ABX*I_NAI_Ix2y3z_S_a;
  Double I_NAI_Ixy4z_Px_a = I_NAI_K2xy4z_S_a+ABX*I_NAI_Ixy4z_S_a;
  Double I_NAI_Ix5z_Px_a = I_NAI_K2x5z_S_a+ABX*I_NAI_Ix5z_S_a;
  Double I_NAI_I5yz_Px_a = I_NAI_Kx5yz_S_a+ABX*I_NAI_I5yz_S_a;
  Double I_NAI_I4y2z_Px_a = I_NAI_Kx4y2z_S_a+ABX*I_NAI_I4y2z_S_a;
  Double I_NAI_I3y3z_Px_a = I_NAI_Kx3y3z_S_a+ABX*I_NAI_I3y3z_S_a;
  Double I_NAI_I2y4z_Px_a = I_NAI_Kx2y4z_S_a+ABX*I_NAI_I2y4z_S_a;
  Double I_NAI_Iy5z_Px_a = I_NAI_Kxy5z_S_a+ABX*I_NAI_Iy5z_S_a;
  Double I_NAI_I5xy_Py_a = I_NAI_K5x2y_S_a+ABY*I_NAI_I5xy_S_a;
  Double I_NAI_I4x2y_Py_a = I_NAI_K4x3y_S_a+ABY*I_NAI_I4x2y_S_a;
  Double I_NAI_I4xyz_Py_a = I_NAI_K4x2yz_S_a+ABY*I_NAI_I4xyz_S_a;
  Double I_NAI_I3x3y_Py_a = I_NAI_K3x4y_S_a+ABY*I_NAI_I3x3y_S_a;
  Double I_NAI_I3x2yz_Py_a = I_NAI_K3x3yz_S_a+ABY*I_NAI_I3x2yz_S_a;
  Double I_NAI_I3xy2z_Py_a = I_NAI_K3x2y2z_S_a+ABY*I_NAI_I3xy2z_S_a;
  Double I_NAI_I2x4y_Py_a = I_NAI_K2x5y_S_a+ABY*I_NAI_I2x4y_S_a;
  Double I_NAI_I2x3yz_Py_a = I_NAI_K2x4yz_S_a+ABY*I_NAI_I2x3yz_S_a;
  Double I_NAI_I2x2y2z_Py_a = I_NAI_K2x3y2z_S_a+ABY*I_NAI_I2x2y2z_S_a;
  Double I_NAI_I2xy3z_Py_a = I_NAI_K2x2y3z_S_a+ABY*I_NAI_I2xy3z_S_a;
  Double I_NAI_Ix5y_Py_a = I_NAI_Kx6y_S_a+ABY*I_NAI_Ix5y_S_a;
  Double I_NAI_Ix4yz_Py_a = I_NAI_Kx5yz_S_a+ABY*I_NAI_Ix4yz_S_a;
  Double I_NAI_Ix3y2z_Py_a = I_NAI_Kx4y2z_S_a+ABY*I_NAI_Ix3y2z_S_a;
  Double I_NAI_Ix2y3z_Py_a = I_NAI_Kx3y3z_S_a+ABY*I_NAI_Ix2y3z_S_a;
  Double I_NAI_Ixy4z_Py_a = I_NAI_Kx2y4z_S_a+ABY*I_NAI_Ixy4z_S_a;
  Double I_NAI_I6y_Py_a = I_NAI_K7y_S_a+ABY*I_NAI_I6y_S_a;
  Double I_NAI_I5yz_Py_a = I_NAI_K6yz_S_a+ABY*I_NAI_I5yz_S_a;
  Double I_NAI_I4y2z_Py_a = I_NAI_K5y2z_S_a+ABY*I_NAI_I4y2z_S_a;
  Double I_NAI_I3y3z_Py_a = I_NAI_K4y3z_S_a+ABY*I_NAI_I3y3z_S_a;
  Double I_NAI_I2y4z_Py_a = I_NAI_K3y4z_S_a+ABY*I_NAI_I2y4z_S_a;
  Double I_NAI_Iy5z_Py_a = I_NAI_K2y5z_S_a+ABY*I_NAI_Iy5z_S_a;
  Double I_NAI_I5xz_Pz_a = I_NAI_K5x2z_S_a+ABZ*I_NAI_I5xz_S_a;
  Double I_NAI_I4xyz_Pz_a = I_NAI_K4xy2z_S_a+ABZ*I_NAI_I4xyz_S_a;
  Double I_NAI_I4x2z_Pz_a = I_NAI_K4x3z_S_a+ABZ*I_NAI_I4x2z_S_a;
  Double I_NAI_I3x2yz_Pz_a = I_NAI_K3x2y2z_S_a+ABZ*I_NAI_I3x2yz_S_a;
  Double I_NAI_I3xy2z_Pz_a = I_NAI_K3xy3z_S_a+ABZ*I_NAI_I3xy2z_S_a;
  Double I_NAI_I3x3z_Pz_a = I_NAI_K3x4z_S_a+ABZ*I_NAI_I3x3z_S_a;
  Double I_NAI_I2x3yz_Pz_a = I_NAI_K2x3y2z_S_a+ABZ*I_NAI_I2x3yz_S_a;
  Double I_NAI_I2x2y2z_Pz_a = I_NAI_K2x2y3z_S_a+ABZ*I_NAI_I2x2y2z_S_a;
  Double I_NAI_I2xy3z_Pz_a = I_NAI_K2xy4z_S_a+ABZ*I_NAI_I2xy3z_S_a;
  Double I_NAI_I2x4z_Pz_a = I_NAI_K2x5z_S_a+ABZ*I_NAI_I2x4z_S_a;
  Double I_NAI_Ix4yz_Pz_a = I_NAI_Kx4y2z_S_a+ABZ*I_NAI_Ix4yz_S_a;
  Double I_NAI_Ix3y2z_Pz_a = I_NAI_Kx3y3z_S_a+ABZ*I_NAI_Ix3y2z_S_a;
  Double I_NAI_Ix2y3z_Pz_a = I_NAI_Kx2y4z_S_a+ABZ*I_NAI_Ix2y3z_S_a;
  Double I_NAI_Ixy4z_Pz_a = I_NAI_Kxy5z_S_a+ABZ*I_NAI_Ixy4z_S_a;
  Double I_NAI_Ix5z_Pz_a = I_NAI_Kx6z_S_a+ABZ*I_NAI_Ix5z_S_a;
  Double I_NAI_I5yz_Pz_a = I_NAI_K5y2z_S_a+ABZ*I_NAI_I5yz_S_a;
  Double I_NAI_I4y2z_Pz_a = I_NAI_K4y3z_S_a+ABZ*I_NAI_I4y2z_S_a;
  Double I_NAI_I3y3z_Pz_a = I_NAI_K3y4z_S_a+ABZ*I_NAI_I3y3z_S_a;
  Double I_NAI_I2y4z_Pz_a = I_NAI_K2y5z_S_a+ABZ*I_NAI_I2y4z_S_a;
  Double I_NAI_Iy5z_Pz_a = I_NAI_Ky6z_S_a+ABZ*I_NAI_Iy5z_S_a;
  Double I_NAI_I6z_Pz_a = I_NAI_K7z_S_a+ABZ*I_NAI_I6z_S_a;

  /************************************************************
   * shell quartet name: SQ_NAI_H_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 48 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_a
   * RHS shell quartet name: SQ_NAI_H_P_a
   ************************************************************/
  Double I_NAI_H5x_D2x_a = I_NAI_I6x_Px_a+ABX*I_NAI_H5x_Px_a;
  Double I_NAI_H4xy_D2x_a = I_NAI_I5xy_Px_a+ABX*I_NAI_H4xy_Px_a;
  Double I_NAI_H4xz_D2x_a = I_NAI_I5xz_Px_a+ABX*I_NAI_H4xz_Px_a;
  Double I_NAI_H3x2y_D2x_a = I_NAI_I4x2y_Px_a+ABX*I_NAI_H3x2y_Px_a;
  Double I_NAI_H3xyz_D2x_a = I_NAI_I4xyz_Px_a+ABX*I_NAI_H3xyz_Px_a;
  Double I_NAI_H3x2z_D2x_a = I_NAI_I4x2z_Px_a+ABX*I_NAI_H3x2z_Px_a;
  Double I_NAI_H2x3y_D2x_a = I_NAI_I3x3y_Px_a+ABX*I_NAI_H2x3y_Px_a;
  Double I_NAI_H2x2yz_D2x_a = I_NAI_I3x2yz_Px_a+ABX*I_NAI_H2x2yz_Px_a;
  Double I_NAI_H2xy2z_D2x_a = I_NAI_I3xy2z_Px_a+ABX*I_NAI_H2xy2z_Px_a;
  Double I_NAI_H2x3z_D2x_a = I_NAI_I3x3z_Px_a+ABX*I_NAI_H2x3z_Px_a;
  Double I_NAI_Hx4y_D2x_a = I_NAI_I2x4y_Px_a+ABX*I_NAI_Hx4y_Px_a;
  Double I_NAI_Hx3yz_D2x_a = I_NAI_I2x3yz_Px_a+ABX*I_NAI_Hx3yz_Px_a;
  Double I_NAI_Hx2y2z_D2x_a = I_NAI_I2x2y2z_Px_a+ABX*I_NAI_Hx2y2z_Px_a;
  Double I_NAI_Hxy3z_D2x_a = I_NAI_I2xy3z_Px_a+ABX*I_NAI_Hxy3z_Px_a;
  Double I_NAI_Hx4z_D2x_a = I_NAI_I2x4z_Px_a+ABX*I_NAI_Hx4z_Px_a;
  Double I_NAI_H5y_D2x_a = I_NAI_Ix5y_Px_a+ABX*I_NAI_H5y_Px_a;
  Double I_NAI_H4yz_D2x_a = I_NAI_Ix4yz_Px_a+ABX*I_NAI_H4yz_Px_a;
  Double I_NAI_H3y2z_D2x_a = I_NAI_Ix3y2z_Px_a+ABX*I_NAI_H3y2z_Px_a;
  Double I_NAI_H2y3z_D2x_a = I_NAI_Ix2y3z_Px_a+ABX*I_NAI_H2y3z_Px_a;
  Double I_NAI_Hy4z_D2x_a = I_NAI_Ixy4z_Px_a+ABX*I_NAI_Hy4z_Px_a;
  Double I_NAI_H5z_D2x_a = I_NAI_Ix5z_Px_a+ABX*I_NAI_H5z_Px_a;
  Double I_NAI_H4xz_Dxy_a = I_NAI_I4xyz_Px_a+ABY*I_NAI_H4xz_Px_a;
  Double I_NAI_H3xyz_Dxy_a = I_NAI_I3x2yz_Px_a+ABY*I_NAI_H3xyz_Px_a;
  Double I_NAI_H3x2z_Dxy_a = I_NAI_I3xy2z_Px_a+ABY*I_NAI_H3x2z_Px_a;
  Double I_NAI_H2x2yz_Dxy_a = I_NAI_I2x3yz_Px_a+ABY*I_NAI_H2x2yz_Px_a;
  Double I_NAI_H2xy2z_Dxy_a = I_NAI_I2x2y2z_Px_a+ABY*I_NAI_H2xy2z_Px_a;
  Double I_NAI_H2x3z_Dxy_a = I_NAI_I2xy3z_Px_a+ABY*I_NAI_H2x3z_Px_a;
  Double I_NAI_Hx3yz_Dxy_a = I_NAI_Ix4yz_Px_a+ABY*I_NAI_Hx3yz_Px_a;
  Double I_NAI_Hx2y2z_Dxy_a = I_NAI_Ix3y2z_Px_a+ABY*I_NAI_Hx2y2z_Px_a;
  Double I_NAI_Hxy3z_Dxy_a = I_NAI_Ix2y3z_Px_a+ABY*I_NAI_Hxy3z_Px_a;
  Double I_NAI_Hx4z_Dxy_a = I_NAI_Ixy4z_Px_a+ABY*I_NAI_Hx4z_Px_a;
  Double I_NAI_H4yz_Dxy_a = I_NAI_I5yz_Px_a+ABY*I_NAI_H4yz_Px_a;
  Double I_NAI_H3y2z_Dxy_a = I_NAI_I4y2z_Px_a+ABY*I_NAI_H3y2z_Px_a;
  Double I_NAI_H2y3z_Dxy_a = I_NAI_I3y3z_Px_a+ABY*I_NAI_H2y3z_Px_a;
  Double I_NAI_Hy4z_Dxy_a = I_NAI_I2y4z_Px_a+ABY*I_NAI_Hy4z_Px_a;
  Double I_NAI_H5z_Dxy_a = I_NAI_Iy5z_Px_a+ABY*I_NAI_H5z_Px_a;
  Double I_NAI_H5x_D2y_a = I_NAI_I5xy_Py_a+ABY*I_NAI_H5x_Py_a;
  Double I_NAI_H4xy_D2y_a = I_NAI_I4x2y_Py_a+ABY*I_NAI_H4xy_Py_a;
  Double I_NAI_H4xz_D2y_a = I_NAI_I4xyz_Py_a+ABY*I_NAI_H4xz_Py_a;
  Double I_NAI_H3x2y_D2y_a = I_NAI_I3x3y_Py_a+ABY*I_NAI_H3x2y_Py_a;
  Double I_NAI_H3xyz_D2y_a = I_NAI_I3x2yz_Py_a+ABY*I_NAI_H3xyz_Py_a;
  Double I_NAI_H3x2z_D2y_a = I_NAI_I3xy2z_Py_a+ABY*I_NAI_H3x2z_Py_a;
  Double I_NAI_H2x3y_D2y_a = I_NAI_I2x4y_Py_a+ABY*I_NAI_H2x3y_Py_a;
  Double I_NAI_H2x2yz_D2y_a = I_NAI_I2x3yz_Py_a+ABY*I_NAI_H2x2yz_Py_a;
  Double I_NAI_H2xy2z_D2y_a = I_NAI_I2x2y2z_Py_a+ABY*I_NAI_H2xy2z_Py_a;
  Double I_NAI_H2x3z_D2y_a = I_NAI_I2xy3z_Py_a+ABY*I_NAI_H2x3z_Py_a;
  Double I_NAI_Hx4y_D2y_a = I_NAI_Ix5y_Py_a+ABY*I_NAI_Hx4y_Py_a;
  Double I_NAI_Hx3yz_D2y_a = I_NAI_Ix4yz_Py_a+ABY*I_NAI_Hx3yz_Py_a;
  Double I_NAI_Hx2y2z_D2y_a = I_NAI_Ix3y2z_Py_a+ABY*I_NAI_Hx2y2z_Py_a;
  Double I_NAI_Hxy3z_D2y_a = I_NAI_Ix2y3z_Py_a+ABY*I_NAI_Hxy3z_Py_a;
  Double I_NAI_Hx4z_D2y_a = I_NAI_Ixy4z_Py_a+ABY*I_NAI_Hx4z_Py_a;
  Double I_NAI_H5y_D2y_a = I_NAI_I6y_Py_a+ABY*I_NAI_H5y_Py_a;
  Double I_NAI_H4yz_D2y_a = I_NAI_I5yz_Py_a+ABY*I_NAI_H4yz_Py_a;
  Double I_NAI_H3y2z_D2y_a = I_NAI_I4y2z_Py_a+ABY*I_NAI_H3y2z_Py_a;
  Double I_NAI_H2y3z_D2y_a = I_NAI_I3y3z_Py_a+ABY*I_NAI_H2y3z_Py_a;
  Double I_NAI_Hy4z_D2y_a = I_NAI_I2y4z_Py_a+ABY*I_NAI_Hy4z_Py_a;
  Double I_NAI_H5z_D2y_a = I_NAI_Iy5z_Py_a+ABY*I_NAI_H5z_Py_a;
  Double I_NAI_H5x_D2z_a = I_NAI_I5xz_Pz_a+ABZ*I_NAI_H5x_Pz_a;
  Double I_NAI_H4xy_D2z_a = I_NAI_I4xyz_Pz_a+ABZ*I_NAI_H4xy_Pz_a;
  Double I_NAI_H4xz_D2z_a = I_NAI_I4x2z_Pz_a+ABZ*I_NAI_H4xz_Pz_a;
  Double I_NAI_H3x2y_D2z_a = I_NAI_I3x2yz_Pz_a+ABZ*I_NAI_H3x2y_Pz_a;
  Double I_NAI_H3xyz_D2z_a = I_NAI_I3xy2z_Pz_a+ABZ*I_NAI_H3xyz_Pz_a;
  Double I_NAI_H3x2z_D2z_a = I_NAI_I3x3z_Pz_a+ABZ*I_NAI_H3x2z_Pz_a;
  Double I_NAI_H2x3y_D2z_a = I_NAI_I2x3yz_Pz_a+ABZ*I_NAI_H2x3y_Pz_a;
  Double I_NAI_H2x2yz_D2z_a = I_NAI_I2x2y2z_Pz_a+ABZ*I_NAI_H2x2yz_Pz_a;
  Double I_NAI_H2xy2z_D2z_a = I_NAI_I2xy3z_Pz_a+ABZ*I_NAI_H2xy2z_Pz_a;
  Double I_NAI_H2x3z_D2z_a = I_NAI_I2x4z_Pz_a+ABZ*I_NAI_H2x3z_Pz_a;
  Double I_NAI_Hx4y_D2z_a = I_NAI_Ix4yz_Pz_a+ABZ*I_NAI_Hx4y_Pz_a;
  Double I_NAI_Hx3yz_D2z_a = I_NAI_Ix3y2z_Pz_a+ABZ*I_NAI_Hx3yz_Pz_a;
  Double I_NAI_Hx2y2z_D2z_a = I_NAI_Ix2y3z_Pz_a+ABZ*I_NAI_Hx2y2z_Pz_a;
  Double I_NAI_Hxy3z_D2z_a = I_NAI_Ixy4z_Pz_a+ABZ*I_NAI_Hxy3z_Pz_a;
  Double I_NAI_Hx4z_D2z_a = I_NAI_Ix5z_Pz_a+ABZ*I_NAI_Hx4z_Pz_a;
  Double I_NAI_H5y_D2z_a = I_NAI_I5yz_Pz_a+ABZ*I_NAI_H5y_Pz_a;
  Double I_NAI_H4yz_D2z_a = I_NAI_I4y2z_Pz_a+ABZ*I_NAI_H4yz_Pz_a;
  Double I_NAI_H3y2z_D2z_a = I_NAI_I3y3z_Pz_a+ABZ*I_NAI_H3y2z_Pz_a;
  Double I_NAI_H2y3z_D2z_a = I_NAI_I2y4z_Pz_a+ABZ*I_NAI_H2y3z_Pz_a;
  Double I_NAI_Hy4z_D2z_a = I_NAI_Iy5z_Pz_a+ABZ*I_NAI_Hy4z_Pz_a;
  Double I_NAI_H5z_D2z_a = I_NAI_I6z_Pz_a+ABZ*I_NAI_H5z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_D_a
   * RHS shell quartet name: SQ_NAI_G_D_a
   ************************************************************/
  Double I_NAI_G4x_F3x_a = I_NAI_H5x_D2x_a+ABX*I_NAI_G4x_D2x_a;
  Double I_NAI_G3xy_F3x_a = I_NAI_H4xy_D2x_a+ABX*I_NAI_G3xy_D2x_a;
  Double I_NAI_G3xz_F3x_a = I_NAI_H4xz_D2x_a+ABX*I_NAI_G3xz_D2x_a;
  Double I_NAI_G2x2y_F3x_a = I_NAI_H3x2y_D2x_a+ABX*I_NAI_G2x2y_D2x_a;
  Double I_NAI_G2xyz_F3x_a = I_NAI_H3xyz_D2x_a+ABX*I_NAI_G2xyz_D2x_a;
  Double I_NAI_G2x2z_F3x_a = I_NAI_H3x2z_D2x_a+ABX*I_NAI_G2x2z_D2x_a;
  Double I_NAI_Gx3y_F3x_a = I_NAI_H2x3y_D2x_a+ABX*I_NAI_Gx3y_D2x_a;
  Double I_NAI_Gx2yz_F3x_a = I_NAI_H2x2yz_D2x_a+ABX*I_NAI_Gx2yz_D2x_a;
  Double I_NAI_Gxy2z_F3x_a = I_NAI_H2xy2z_D2x_a+ABX*I_NAI_Gxy2z_D2x_a;
  Double I_NAI_Gx3z_F3x_a = I_NAI_H2x3z_D2x_a+ABX*I_NAI_Gx3z_D2x_a;
  Double I_NAI_G4y_F3x_a = I_NAI_Hx4y_D2x_a+ABX*I_NAI_G4y_D2x_a;
  Double I_NAI_G3yz_F3x_a = I_NAI_Hx3yz_D2x_a+ABX*I_NAI_G3yz_D2x_a;
  Double I_NAI_G2y2z_F3x_a = I_NAI_Hx2y2z_D2x_a+ABX*I_NAI_G2y2z_D2x_a;
  Double I_NAI_Gy3z_F3x_a = I_NAI_Hxy3z_D2x_a+ABX*I_NAI_Gy3z_D2x_a;
  Double I_NAI_G4z_F3x_a = I_NAI_Hx4z_D2x_a+ABX*I_NAI_G4z_D2x_a;
  Double I_NAI_G4x_F2xy_a = I_NAI_H4xy_D2x_a+ABY*I_NAI_G4x_D2x_a;
  Double I_NAI_G3xy_F2xy_a = I_NAI_H3x2y_D2x_a+ABY*I_NAI_G3xy_D2x_a;
  Double I_NAI_G3xz_F2xy_a = I_NAI_H3xyz_D2x_a+ABY*I_NAI_G3xz_D2x_a;
  Double I_NAI_G2x2y_F2xy_a = I_NAI_H2x3y_D2x_a+ABY*I_NAI_G2x2y_D2x_a;
  Double I_NAI_G2xyz_F2xy_a = I_NAI_H2x2yz_D2x_a+ABY*I_NAI_G2xyz_D2x_a;
  Double I_NAI_G2x2z_F2xy_a = I_NAI_H2xy2z_D2x_a+ABY*I_NAI_G2x2z_D2x_a;
  Double I_NAI_Gx3y_F2xy_a = I_NAI_Hx4y_D2x_a+ABY*I_NAI_Gx3y_D2x_a;
  Double I_NAI_Gx2yz_F2xy_a = I_NAI_Hx3yz_D2x_a+ABY*I_NAI_Gx2yz_D2x_a;
  Double I_NAI_Gxy2z_F2xy_a = I_NAI_Hx2y2z_D2x_a+ABY*I_NAI_Gxy2z_D2x_a;
  Double I_NAI_Gx3z_F2xy_a = I_NAI_Hxy3z_D2x_a+ABY*I_NAI_Gx3z_D2x_a;
  Double I_NAI_G4y_F2xy_a = I_NAI_H5y_D2x_a+ABY*I_NAI_G4y_D2x_a;
  Double I_NAI_G3yz_F2xy_a = I_NAI_H4yz_D2x_a+ABY*I_NAI_G3yz_D2x_a;
  Double I_NAI_G2y2z_F2xy_a = I_NAI_H3y2z_D2x_a+ABY*I_NAI_G2y2z_D2x_a;
  Double I_NAI_Gy3z_F2xy_a = I_NAI_H2y3z_D2x_a+ABY*I_NAI_Gy3z_D2x_a;
  Double I_NAI_G4z_F2xy_a = I_NAI_Hy4z_D2x_a+ABY*I_NAI_G4z_D2x_a;
  Double I_NAI_G4x_F2xz_a = I_NAI_H4xz_D2x_a+ABZ*I_NAI_G4x_D2x_a;
  Double I_NAI_G3xy_F2xz_a = I_NAI_H3xyz_D2x_a+ABZ*I_NAI_G3xy_D2x_a;
  Double I_NAI_G3xz_F2xz_a = I_NAI_H3x2z_D2x_a+ABZ*I_NAI_G3xz_D2x_a;
  Double I_NAI_G2x2y_F2xz_a = I_NAI_H2x2yz_D2x_a+ABZ*I_NAI_G2x2y_D2x_a;
  Double I_NAI_G2xyz_F2xz_a = I_NAI_H2xy2z_D2x_a+ABZ*I_NAI_G2xyz_D2x_a;
  Double I_NAI_G2x2z_F2xz_a = I_NAI_H2x3z_D2x_a+ABZ*I_NAI_G2x2z_D2x_a;
  Double I_NAI_Gx3y_F2xz_a = I_NAI_Hx3yz_D2x_a+ABZ*I_NAI_Gx3y_D2x_a;
  Double I_NAI_Gx2yz_F2xz_a = I_NAI_Hx2y2z_D2x_a+ABZ*I_NAI_Gx2yz_D2x_a;
  Double I_NAI_Gxy2z_F2xz_a = I_NAI_Hxy3z_D2x_a+ABZ*I_NAI_Gxy2z_D2x_a;
  Double I_NAI_Gx3z_F2xz_a = I_NAI_Hx4z_D2x_a+ABZ*I_NAI_Gx3z_D2x_a;
  Double I_NAI_G4y_F2xz_a = I_NAI_H4yz_D2x_a+ABZ*I_NAI_G4y_D2x_a;
  Double I_NAI_G3yz_F2xz_a = I_NAI_H3y2z_D2x_a+ABZ*I_NAI_G3yz_D2x_a;
  Double I_NAI_G2y2z_F2xz_a = I_NAI_H2y3z_D2x_a+ABZ*I_NAI_G2y2z_D2x_a;
  Double I_NAI_Gy3z_F2xz_a = I_NAI_Hy4z_D2x_a+ABZ*I_NAI_Gy3z_D2x_a;
  Double I_NAI_G4z_F2xz_a = I_NAI_H5z_D2x_a+ABZ*I_NAI_G4z_D2x_a;
  Double I_NAI_G4x_Fx2y_a = I_NAI_H5x_D2y_a+ABX*I_NAI_G4x_D2y_a;
  Double I_NAI_G3xy_Fx2y_a = I_NAI_H4xy_D2y_a+ABX*I_NAI_G3xy_D2y_a;
  Double I_NAI_G3xz_Fx2y_a = I_NAI_H4xz_D2y_a+ABX*I_NAI_G3xz_D2y_a;
  Double I_NAI_G2x2y_Fx2y_a = I_NAI_H3x2y_D2y_a+ABX*I_NAI_G2x2y_D2y_a;
  Double I_NAI_G2xyz_Fx2y_a = I_NAI_H3xyz_D2y_a+ABX*I_NAI_G2xyz_D2y_a;
  Double I_NAI_G2x2z_Fx2y_a = I_NAI_H3x2z_D2y_a+ABX*I_NAI_G2x2z_D2y_a;
  Double I_NAI_Gx3y_Fx2y_a = I_NAI_H2x3y_D2y_a+ABX*I_NAI_Gx3y_D2y_a;
  Double I_NAI_Gx2yz_Fx2y_a = I_NAI_H2x2yz_D2y_a+ABX*I_NAI_Gx2yz_D2y_a;
  Double I_NAI_Gxy2z_Fx2y_a = I_NAI_H2xy2z_D2y_a+ABX*I_NAI_Gxy2z_D2y_a;
  Double I_NAI_Gx3z_Fx2y_a = I_NAI_H2x3z_D2y_a+ABX*I_NAI_Gx3z_D2y_a;
  Double I_NAI_G4y_Fx2y_a = I_NAI_Hx4y_D2y_a+ABX*I_NAI_G4y_D2y_a;
  Double I_NAI_G3yz_Fx2y_a = I_NAI_Hx3yz_D2y_a+ABX*I_NAI_G3yz_D2y_a;
  Double I_NAI_G2y2z_Fx2y_a = I_NAI_Hx2y2z_D2y_a+ABX*I_NAI_G2y2z_D2y_a;
  Double I_NAI_Gy3z_Fx2y_a = I_NAI_Hxy3z_D2y_a+ABX*I_NAI_Gy3z_D2y_a;
  Double I_NAI_G4z_Fx2y_a = I_NAI_Hx4z_D2y_a+ABX*I_NAI_G4z_D2y_a;
  Double I_NAI_G4x_Fxyz_a = I_NAI_H4xz_Dxy_a+ABZ*I_NAI_G4x_Dxy_a;
  Double I_NAI_G3xy_Fxyz_a = I_NAI_H3xyz_Dxy_a+ABZ*I_NAI_G3xy_Dxy_a;
  Double I_NAI_G3xz_Fxyz_a = I_NAI_H3x2z_Dxy_a+ABZ*I_NAI_G3xz_Dxy_a;
  Double I_NAI_G2x2y_Fxyz_a = I_NAI_H2x2yz_Dxy_a+ABZ*I_NAI_G2x2y_Dxy_a;
  Double I_NAI_G2xyz_Fxyz_a = I_NAI_H2xy2z_Dxy_a+ABZ*I_NAI_G2xyz_Dxy_a;
  Double I_NAI_G2x2z_Fxyz_a = I_NAI_H2x3z_Dxy_a+ABZ*I_NAI_G2x2z_Dxy_a;
  Double I_NAI_Gx3y_Fxyz_a = I_NAI_Hx3yz_Dxy_a+ABZ*I_NAI_Gx3y_Dxy_a;
  Double I_NAI_Gx2yz_Fxyz_a = I_NAI_Hx2y2z_Dxy_a+ABZ*I_NAI_Gx2yz_Dxy_a;
  Double I_NAI_Gxy2z_Fxyz_a = I_NAI_Hxy3z_Dxy_a+ABZ*I_NAI_Gxy2z_Dxy_a;
  Double I_NAI_Gx3z_Fxyz_a = I_NAI_Hx4z_Dxy_a+ABZ*I_NAI_Gx3z_Dxy_a;
  Double I_NAI_G4y_Fxyz_a = I_NAI_H4yz_Dxy_a+ABZ*I_NAI_G4y_Dxy_a;
  Double I_NAI_G3yz_Fxyz_a = I_NAI_H3y2z_Dxy_a+ABZ*I_NAI_G3yz_Dxy_a;
  Double I_NAI_G2y2z_Fxyz_a = I_NAI_H2y3z_Dxy_a+ABZ*I_NAI_G2y2z_Dxy_a;
  Double I_NAI_Gy3z_Fxyz_a = I_NAI_Hy4z_Dxy_a+ABZ*I_NAI_Gy3z_Dxy_a;
  Double I_NAI_G4z_Fxyz_a = I_NAI_H5z_Dxy_a+ABZ*I_NAI_G4z_Dxy_a;
  Double I_NAI_G4x_Fx2z_a = I_NAI_H5x_D2z_a+ABX*I_NAI_G4x_D2z_a;
  Double I_NAI_G3xy_Fx2z_a = I_NAI_H4xy_D2z_a+ABX*I_NAI_G3xy_D2z_a;
  Double I_NAI_G3xz_Fx2z_a = I_NAI_H4xz_D2z_a+ABX*I_NAI_G3xz_D2z_a;
  Double I_NAI_G2x2y_Fx2z_a = I_NAI_H3x2y_D2z_a+ABX*I_NAI_G2x2y_D2z_a;
  Double I_NAI_G2xyz_Fx2z_a = I_NAI_H3xyz_D2z_a+ABX*I_NAI_G2xyz_D2z_a;
  Double I_NAI_G2x2z_Fx2z_a = I_NAI_H3x2z_D2z_a+ABX*I_NAI_G2x2z_D2z_a;
  Double I_NAI_Gx3y_Fx2z_a = I_NAI_H2x3y_D2z_a+ABX*I_NAI_Gx3y_D2z_a;
  Double I_NAI_Gx2yz_Fx2z_a = I_NAI_H2x2yz_D2z_a+ABX*I_NAI_Gx2yz_D2z_a;
  Double I_NAI_Gxy2z_Fx2z_a = I_NAI_H2xy2z_D2z_a+ABX*I_NAI_Gxy2z_D2z_a;
  Double I_NAI_Gx3z_Fx2z_a = I_NAI_H2x3z_D2z_a+ABX*I_NAI_Gx3z_D2z_a;
  Double I_NAI_G4y_Fx2z_a = I_NAI_Hx4y_D2z_a+ABX*I_NAI_G4y_D2z_a;
  Double I_NAI_G3yz_Fx2z_a = I_NAI_Hx3yz_D2z_a+ABX*I_NAI_G3yz_D2z_a;
  Double I_NAI_G2y2z_Fx2z_a = I_NAI_Hx2y2z_D2z_a+ABX*I_NAI_G2y2z_D2z_a;
  Double I_NAI_Gy3z_Fx2z_a = I_NAI_Hxy3z_D2z_a+ABX*I_NAI_Gy3z_D2z_a;
  Double I_NAI_G4z_Fx2z_a = I_NAI_Hx4z_D2z_a+ABX*I_NAI_G4z_D2z_a;
  Double I_NAI_G4x_F3y_a = I_NAI_H4xy_D2y_a+ABY*I_NAI_G4x_D2y_a;
  Double I_NAI_G3xy_F3y_a = I_NAI_H3x2y_D2y_a+ABY*I_NAI_G3xy_D2y_a;
  Double I_NAI_G3xz_F3y_a = I_NAI_H3xyz_D2y_a+ABY*I_NAI_G3xz_D2y_a;
  Double I_NAI_G2x2y_F3y_a = I_NAI_H2x3y_D2y_a+ABY*I_NAI_G2x2y_D2y_a;
  Double I_NAI_G2xyz_F3y_a = I_NAI_H2x2yz_D2y_a+ABY*I_NAI_G2xyz_D2y_a;
  Double I_NAI_G2x2z_F3y_a = I_NAI_H2xy2z_D2y_a+ABY*I_NAI_G2x2z_D2y_a;
  Double I_NAI_Gx3y_F3y_a = I_NAI_Hx4y_D2y_a+ABY*I_NAI_Gx3y_D2y_a;
  Double I_NAI_Gx2yz_F3y_a = I_NAI_Hx3yz_D2y_a+ABY*I_NAI_Gx2yz_D2y_a;
  Double I_NAI_Gxy2z_F3y_a = I_NAI_Hx2y2z_D2y_a+ABY*I_NAI_Gxy2z_D2y_a;
  Double I_NAI_Gx3z_F3y_a = I_NAI_Hxy3z_D2y_a+ABY*I_NAI_Gx3z_D2y_a;
  Double I_NAI_G4y_F3y_a = I_NAI_H5y_D2y_a+ABY*I_NAI_G4y_D2y_a;
  Double I_NAI_G3yz_F3y_a = I_NAI_H4yz_D2y_a+ABY*I_NAI_G3yz_D2y_a;
  Double I_NAI_G2y2z_F3y_a = I_NAI_H3y2z_D2y_a+ABY*I_NAI_G2y2z_D2y_a;
  Double I_NAI_Gy3z_F3y_a = I_NAI_H2y3z_D2y_a+ABY*I_NAI_Gy3z_D2y_a;
  Double I_NAI_G4z_F3y_a = I_NAI_Hy4z_D2y_a+ABY*I_NAI_G4z_D2y_a;
  Double I_NAI_G4x_F2yz_a = I_NAI_H4xz_D2y_a+ABZ*I_NAI_G4x_D2y_a;
  Double I_NAI_G3xy_F2yz_a = I_NAI_H3xyz_D2y_a+ABZ*I_NAI_G3xy_D2y_a;
  Double I_NAI_G3xz_F2yz_a = I_NAI_H3x2z_D2y_a+ABZ*I_NAI_G3xz_D2y_a;
  Double I_NAI_G2x2y_F2yz_a = I_NAI_H2x2yz_D2y_a+ABZ*I_NAI_G2x2y_D2y_a;
  Double I_NAI_G2xyz_F2yz_a = I_NAI_H2xy2z_D2y_a+ABZ*I_NAI_G2xyz_D2y_a;
  Double I_NAI_G2x2z_F2yz_a = I_NAI_H2x3z_D2y_a+ABZ*I_NAI_G2x2z_D2y_a;
  Double I_NAI_Gx3y_F2yz_a = I_NAI_Hx3yz_D2y_a+ABZ*I_NAI_Gx3y_D2y_a;
  Double I_NAI_Gx2yz_F2yz_a = I_NAI_Hx2y2z_D2y_a+ABZ*I_NAI_Gx2yz_D2y_a;
  Double I_NAI_Gxy2z_F2yz_a = I_NAI_Hxy3z_D2y_a+ABZ*I_NAI_Gxy2z_D2y_a;
  Double I_NAI_Gx3z_F2yz_a = I_NAI_Hx4z_D2y_a+ABZ*I_NAI_Gx3z_D2y_a;
  Double I_NAI_G4y_F2yz_a = I_NAI_H4yz_D2y_a+ABZ*I_NAI_G4y_D2y_a;
  Double I_NAI_G3yz_F2yz_a = I_NAI_H3y2z_D2y_a+ABZ*I_NAI_G3yz_D2y_a;
  Double I_NAI_G2y2z_F2yz_a = I_NAI_H2y3z_D2y_a+ABZ*I_NAI_G2y2z_D2y_a;
  Double I_NAI_Gy3z_F2yz_a = I_NAI_Hy4z_D2y_a+ABZ*I_NAI_Gy3z_D2y_a;
  Double I_NAI_G4z_F2yz_a = I_NAI_H5z_D2y_a+ABZ*I_NAI_G4z_D2y_a;
  Double I_NAI_G4x_Fy2z_a = I_NAI_H4xy_D2z_a+ABY*I_NAI_G4x_D2z_a;
  Double I_NAI_G3xy_Fy2z_a = I_NAI_H3x2y_D2z_a+ABY*I_NAI_G3xy_D2z_a;
  Double I_NAI_G3xz_Fy2z_a = I_NAI_H3xyz_D2z_a+ABY*I_NAI_G3xz_D2z_a;
  Double I_NAI_G2x2y_Fy2z_a = I_NAI_H2x3y_D2z_a+ABY*I_NAI_G2x2y_D2z_a;
  Double I_NAI_G2xyz_Fy2z_a = I_NAI_H2x2yz_D2z_a+ABY*I_NAI_G2xyz_D2z_a;
  Double I_NAI_G2x2z_Fy2z_a = I_NAI_H2xy2z_D2z_a+ABY*I_NAI_G2x2z_D2z_a;
  Double I_NAI_Gx3y_Fy2z_a = I_NAI_Hx4y_D2z_a+ABY*I_NAI_Gx3y_D2z_a;
  Double I_NAI_Gx2yz_Fy2z_a = I_NAI_Hx3yz_D2z_a+ABY*I_NAI_Gx2yz_D2z_a;
  Double I_NAI_Gxy2z_Fy2z_a = I_NAI_Hx2y2z_D2z_a+ABY*I_NAI_Gxy2z_D2z_a;
  Double I_NAI_Gx3z_Fy2z_a = I_NAI_Hxy3z_D2z_a+ABY*I_NAI_Gx3z_D2z_a;
  Double I_NAI_G4y_Fy2z_a = I_NAI_H5y_D2z_a+ABY*I_NAI_G4y_D2z_a;
  Double I_NAI_G3yz_Fy2z_a = I_NAI_H4yz_D2z_a+ABY*I_NAI_G3yz_D2z_a;
  Double I_NAI_G2y2z_Fy2z_a = I_NAI_H3y2z_D2z_a+ABY*I_NAI_G2y2z_D2z_a;
  Double I_NAI_Gy3z_Fy2z_a = I_NAI_H2y3z_D2z_a+ABY*I_NAI_Gy3z_D2z_a;
  Double I_NAI_G4z_Fy2z_a = I_NAI_Hy4z_D2z_a+ABY*I_NAI_G4z_D2z_a;
  Double I_NAI_G4x_F3z_a = I_NAI_H4xz_D2z_a+ABZ*I_NAI_G4x_D2z_a;
  Double I_NAI_G3xy_F3z_a = I_NAI_H3xyz_D2z_a+ABZ*I_NAI_G3xy_D2z_a;
  Double I_NAI_G3xz_F3z_a = I_NAI_H3x2z_D2z_a+ABZ*I_NAI_G3xz_D2z_a;
  Double I_NAI_G2x2y_F3z_a = I_NAI_H2x2yz_D2z_a+ABZ*I_NAI_G2x2y_D2z_a;
  Double I_NAI_G2xyz_F3z_a = I_NAI_H2xy2z_D2z_a+ABZ*I_NAI_G2xyz_D2z_a;
  Double I_NAI_G2x2z_F3z_a = I_NAI_H2x3z_D2z_a+ABZ*I_NAI_G2x2z_D2z_a;
  Double I_NAI_Gx3y_F3z_a = I_NAI_Hx3yz_D2z_a+ABZ*I_NAI_Gx3y_D2z_a;
  Double I_NAI_Gx2yz_F3z_a = I_NAI_Hx2y2z_D2z_a+ABZ*I_NAI_Gx2yz_D2z_a;
  Double I_NAI_Gxy2z_F3z_a = I_NAI_Hxy3z_D2z_a+ABZ*I_NAI_Gxy2z_D2z_a;
  Double I_NAI_Gx3z_F3z_a = I_NAI_Hx4z_D2z_a+ABZ*I_NAI_Gx3z_D2z_a;
  Double I_NAI_G4y_F3z_a = I_NAI_H4yz_D2z_a+ABZ*I_NAI_G4y_D2z_a;
  Double I_NAI_G3yz_F3z_a = I_NAI_H3y2z_D2z_a+ABZ*I_NAI_G3yz_D2z_a;
  Double I_NAI_G2y2z_F3z_a = I_NAI_H2y3z_D2z_a+ABZ*I_NAI_G2y2z_D2z_a;
  Double I_NAI_Gy3z_F3z_a = I_NAI_Hy4z_D2z_a+ABZ*I_NAI_Gy3z_D2z_a;
  Double I_NAI_G4z_F3z_a = I_NAI_H5z_D2z_a+ABZ*I_NAI_G4z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_S_b
   * RHS shell quartet name: SQ_NAI_F_S_b
   ************************************************************/
  Double I_NAI_F3x_Px_b = I_NAI_G4x_S_b+ABX*I_NAI_F3x_S_b;
  Double I_NAI_F2xy_Px_b = I_NAI_G3xy_S_b+ABX*I_NAI_F2xy_S_b;
  Double I_NAI_F2xz_Px_b = I_NAI_G3xz_S_b+ABX*I_NAI_F2xz_S_b;
  Double I_NAI_Fx2y_Px_b = I_NAI_G2x2y_S_b+ABX*I_NAI_Fx2y_S_b;
  Double I_NAI_Fxyz_Px_b = I_NAI_G2xyz_S_b+ABX*I_NAI_Fxyz_S_b;
  Double I_NAI_Fx2z_Px_b = I_NAI_G2x2z_S_b+ABX*I_NAI_Fx2z_S_b;
  Double I_NAI_F3y_Px_b = I_NAI_Gx3y_S_b+ABX*I_NAI_F3y_S_b;
  Double I_NAI_F2yz_Px_b = I_NAI_Gx2yz_S_b+ABX*I_NAI_F2yz_S_b;
  Double I_NAI_Fy2z_Px_b = I_NAI_Gxy2z_S_b+ABX*I_NAI_Fy2z_S_b;
  Double I_NAI_F3z_Px_b = I_NAI_Gx3z_S_b+ABX*I_NAI_F3z_S_b;
  Double I_NAI_F3x_Py_b = I_NAI_G3xy_S_b+ABY*I_NAI_F3x_S_b;
  Double I_NAI_F2xy_Py_b = I_NAI_G2x2y_S_b+ABY*I_NAI_F2xy_S_b;
  Double I_NAI_F2xz_Py_b = I_NAI_G2xyz_S_b+ABY*I_NAI_F2xz_S_b;
  Double I_NAI_Fx2y_Py_b = I_NAI_Gx3y_S_b+ABY*I_NAI_Fx2y_S_b;
  Double I_NAI_Fxyz_Py_b = I_NAI_Gx2yz_S_b+ABY*I_NAI_Fxyz_S_b;
  Double I_NAI_Fx2z_Py_b = I_NAI_Gxy2z_S_b+ABY*I_NAI_Fx2z_S_b;
  Double I_NAI_F3y_Py_b = I_NAI_G4y_S_b+ABY*I_NAI_F3y_S_b;
  Double I_NAI_F2yz_Py_b = I_NAI_G3yz_S_b+ABY*I_NAI_F2yz_S_b;
  Double I_NAI_Fy2z_Py_b = I_NAI_G2y2z_S_b+ABY*I_NAI_Fy2z_S_b;
  Double I_NAI_F3z_Py_b = I_NAI_Gy3z_S_b+ABY*I_NAI_F3z_S_b;
  Double I_NAI_F3x_Pz_b = I_NAI_G3xz_S_b+ABZ*I_NAI_F3x_S_b;
  Double I_NAI_F2xy_Pz_b = I_NAI_G2xyz_S_b+ABZ*I_NAI_F2xy_S_b;
  Double I_NAI_F2xz_Pz_b = I_NAI_G2x2z_S_b+ABZ*I_NAI_F2xz_S_b;
  Double I_NAI_Fx2y_Pz_b = I_NAI_Gx2yz_S_b+ABZ*I_NAI_Fx2y_S_b;
  Double I_NAI_Fxyz_Pz_b = I_NAI_Gxy2z_S_b+ABZ*I_NAI_Fxyz_S_b;
  Double I_NAI_Fx2z_Pz_b = I_NAI_Gx3z_S_b+ABZ*I_NAI_Fx2z_S_b;
  Double I_NAI_F3y_Pz_b = I_NAI_G3yz_S_b+ABZ*I_NAI_F3y_S_b;
  Double I_NAI_F2yz_Pz_b = I_NAI_G2y2z_S_b+ABZ*I_NAI_F2yz_S_b;
  Double I_NAI_Fy2z_Pz_b = I_NAI_Gy3z_S_b+ABZ*I_NAI_Fy2z_S_b;
  Double I_NAI_F3z_Pz_b = I_NAI_G4z_S_b+ABZ*I_NAI_F3z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_S_b
   * RHS shell quartet name: SQ_NAI_G_S_b
   ************************************************************/
  Double I_NAI_G4x_Px_b = I_NAI_H5x_S_b+ABX*I_NAI_G4x_S_b;
  Double I_NAI_G3xy_Px_b = I_NAI_H4xy_S_b+ABX*I_NAI_G3xy_S_b;
  Double I_NAI_G3xz_Px_b = I_NAI_H4xz_S_b+ABX*I_NAI_G3xz_S_b;
  Double I_NAI_G2x2y_Px_b = I_NAI_H3x2y_S_b+ABX*I_NAI_G2x2y_S_b;
  Double I_NAI_G2xyz_Px_b = I_NAI_H3xyz_S_b+ABX*I_NAI_G2xyz_S_b;
  Double I_NAI_G2x2z_Px_b = I_NAI_H3x2z_S_b+ABX*I_NAI_G2x2z_S_b;
  Double I_NAI_Gx3y_Px_b = I_NAI_H2x3y_S_b+ABX*I_NAI_Gx3y_S_b;
  Double I_NAI_Gx2yz_Px_b = I_NAI_H2x2yz_S_b+ABX*I_NAI_Gx2yz_S_b;
  Double I_NAI_Gxy2z_Px_b = I_NAI_H2xy2z_S_b+ABX*I_NAI_Gxy2z_S_b;
  Double I_NAI_Gx3z_Px_b = I_NAI_H2x3z_S_b+ABX*I_NAI_Gx3z_S_b;
  Double I_NAI_G4y_Px_b = I_NAI_Hx4y_S_b+ABX*I_NAI_G4y_S_b;
  Double I_NAI_G3yz_Px_b = I_NAI_Hx3yz_S_b+ABX*I_NAI_G3yz_S_b;
  Double I_NAI_G2y2z_Px_b = I_NAI_Hx2y2z_S_b+ABX*I_NAI_G2y2z_S_b;
  Double I_NAI_Gy3z_Px_b = I_NAI_Hxy3z_S_b+ABX*I_NAI_Gy3z_S_b;
  Double I_NAI_G4z_Px_b = I_NAI_Hx4z_S_b+ABX*I_NAI_G4z_S_b;
  Double I_NAI_G4x_Py_b = I_NAI_H4xy_S_b+ABY*I_NAI_G4x_S_b;
  Double I_NAI_G3xy_Py_b = I_NAI_H3x2y_S_b+ABY*I_NAI_G3xy_S_b;
  Double I_NAI_G3xz_Py_b = I_NAI_H3xyz_S_b+ABY*I_NAI_G3xz_S_b;
  Double I_NAI_G2x2y_Py_b = I_NAI_H2x3y_S_b+ABY*I_NAI_G2x2y_S_b;
  Double I_NAI_G2xyz_Py_b = I_NAI_H2x2yz_S_b+ABY*I_NAI_G2xyz_S_b;
  Double I_NAI_G2x2z_Py_b = I_NAI_H2xy2z_S_b+ABY*I_NAI_G2x2z_S_b;
  Double I_NAI_Gx3y_Py_b = I_NAI_Hx4y_S_b+ABY*I_NAI_Gx3y_S_b;
  Double I_NAI_Gx2yz_Py_b = I_NAI_Hx3yz_S_b+ABY*I_NAI_Gx2yz_S_b;
  Double I_NAI_Gxy2z_Py_b = I_NAI_Hx2y2z_S_b+ABY*I_NAI_Gxy2z_S_b;
  Double I_NAI_Gx3z_Py_b = I_NAI_Hxy3z_S_b+ABY*I_NAI_Gx3z_S_b;
  Double I_NAI_G4y_Py_b = I_NAI_H5y_S_b+ABY*I_NAI_G4y_S_b;
  Double I_NAI_G3yz_Py_b = I_NAI_H4yz_S_b+ABY*I_NAI_G3yz_S_b;
  Double I_NAI_G2y2z_Py_b = I_NAI_H3y2z_S_b+ABY*I_NAI_G2y2z_S_b;
  Double I_NAI_Gy3z_Py_b = I_NAI_H2y3z_S_b+ABY*I_NAI_Gy3z_S_b;
  Double I_NAI_G4z_Py_b = I_NAI_Hy4z_S_b+ABY*I_NAI_G4z_S_b;
  Double I_NAI_G4x_Pz_b = I_NAI_H4xz_S_b+ABZ*I_NAI_G4x_S_b;
  Double I_NAI_G3xy_Pz_b = I_NAI_H3xyz_S_b+ABZ*I_NAI_G3xy_S_b;
  Double I_NAI_G3xz_Pz_b = I_NAI_H3x2z_S_b+ABZ*I_NAI_G3xz_S_b;
  Double I_NAI_G2x2y_Pz_b = I_NAI_H2x2yz_S_b+ABZ*I_NAI_G2x2y_S_b;
  Double I_NAI_G2xyz_Pz_b = I_NAI_H2xy2z_S_b+ABZ*I_NAI_G2xyz_S_b;
  Double I_NAI_G2x2z_Pz_b = I_NAI_H2x3z_S_b+ABZ*I_NAI_G2x2z_S_b;
  Double I_NAI_Gx3y_Pz_b = I_NAI_Hx3yz_S_b+ABZ*I_NAI_Gx3y_S_b;
  Double I_NAI_Gx2yz_Pz_b = I_NAI_Hx2y2z_S_b+ABZ*I_NAI_Gx2yz_S_b;
  Double I_NAI_Gxy2z_Pz_b = I_NAI_Hxy3z_S_b+ABZ*I_NAI_Gxy2z_S_b;
  Double I_NAI_Gx3z_Pz_b = I_NAI_Hx4z_S_b+ABZ*I_NAI_Gx3z_S_b;
  Double I_NAI_G4y_Pz_b = I_NAI_H4yz_S_b+ABZ*I_NAI_G4y_S_b;
  Double I_NAI_G3yz_Pz_b = I_NAI_H3y2z_S_b+ABZ*I_NAI_G3yz_S_b;
  Double I_NAI_G2y2z_Pz_b = I_NAI_H2y3z_S_b+ABZ*I_NAI_G2y2z_S_b;
  Double I_NAI_Gy3z_Pz_b = I_NAI_Hy4z_S_b+ABZ*I_NAI_Gy3z_S_b;
  Double I_NAI_G4z_Pz_b = I_NAI_H5z_S_b+ABZ*I_NAI_G4z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_P_b
   * RHS shell quartet name: SQ_NAI_F_P_b
   ************************************************************/
  Double I_NAI_F3x_D2x_b = I_NAI_G4x_Px_b+ABX*I_NAI_F3x_Px_b;
  Double I_NAI_F2xy_D2x_b = I_NAI_G3xy_Px_b+ABX*I_NAI_F2xy_Px_b;
  Double I_NAI_F2xz_D2x_b = I_NAI_G3xz_Px_b+ABX*I_NAI_F2xz_Px_b;
  Double I_NAI_Fx2y_D2x_b = I_NAI_G2x2y_Px_b+ABX*I_NAI_Fx2y_Px_b;
  Double I_NAI_Fxyz_D2x_b = I_NAI_G2xyz_Px_b+ABX*I_NAI_Fxyz_Px_b;
  Double I_NAI_Fx2z_D2x_b = I_NAI_G2x2z_Px_b+ABX*I_NAI_Fx2z_Px_b;
  Double I_NAI_F3y_D2x_b = I_NAI_Gx3y_Px_b+ABX*I_NAI_F3y_Px_b;
  Double I_NAI_F2yz_D2x_b = I_NAI_Gx2yz_Px_b+ABX*I_NAI_F2yz_Px_b;
  Double I_NAI_Fy2z_D2x_b = I_NAI_Gxy2z_Px_b+ABX*I_NAI_Fy2z_Px_b;
  Double I_NAI_F3z_D2x_b = I_NAI_Gx3z_Px_b+ABX*I_NAI_F3z_Px_b;
  Double I_NAI_F3x_D2y_b = I_NAI_G3xy_Py_b+ABY*I_NAI_F3x_Py_b;
  Double I_NAI_F2xy_D2y_b = I_NAI_G2x2y_Py_b+ABY*I_NAI_F2xy_Py_b;
  Double I_NAI_F2xz_D2y_b = I_NAI_G2xyz_Py_b+ABY*I_NAI_F2xz_Py_b;
  Double I_NAI_Fx2y_D2y_b = I_NAI_Gx3y_Py_b+ABY*I_NAI_Fx2y_Py_b;
  Double I_NAI_Fxyz_D2y_b = I_NAI_Gx2yz_Py_b+ABY*I_NAI_Fxyz_Py_b;
  Double I_NAI_Fx2z_D2y_b = I_NAI_Gxy2z_Py_b+ABY*I_NAI_Fx2z_Py_b;
  Double I_NAI_F3y_D2y_b = I_NAI_G4y_Py_b+ABY*I_NAI_F3y_Py_b;
  Double I_NAI_F2yz_D2y_b = I_NAI_G3yz_Py_b+ABY*I_NAI_F2yz_Py_b;
  Double I_NAI_Fy2z_D2y_b = I_NAI_G2y2z_Py_b+ABY*I_NAI_Fy2z_Py_b;
  Double I_NAI_F3z_D2y_b = I_NAI_Gy3z_Py_b+ABY*I_NAI_F3z_Py_b;
  Double I_NAI_F3x_D2z_b = I_NAI_G3xz_Pz_b+ABZ*I_NAI_F3x_Pz_b;
  Double I_NAI_F2xy_D2z_b = I_NAI_G2xyz_Pz_b+ABZ*I_NAI_F2xy_Pz_b;
  Double I_NAI_F2xz_D2z_b = I_NAI_G2x2z_Pz_b+ABZ*I_NAI_F2xz_Pz_b;
  Double I_NAI_Fx2y_D2z_b = I_NAI_Gx2yz_Pz_b+ABZ*I_NAI_Fx2y_Pz_b;
  Double I_NAI_Fxyz_D2z_b = I_NAI_Gxy2z_Pz_b+ABZ*I_NAI_Fxyz_Pz_b;
  Double I_NAI_Fx2z_D2z_b = I_NAI_Gx3z_Pz_b+ABZ*I_NAI_Fx2z_Pz_b;
  Double I_NAI_F3y_D2z_b = I_NAI_G3yz_Pz_b+ABZ*I_NAI_F3y_Pz_b;
  Double I_NAI_F2yz_D2z_b = I_NAI_G2y2z_Pz_b+ABZ*I_NAI_F2yz_Pz_b;
  Double I_NAI_Fy2z_D2z_b = I_NAI_Gy3z_Pz_b+ABZ*I_NAI_Fy2z_Pz_b;
  Double I_NAI_F3z_D2z_b = I_NAI_G4z_Pz_b+ABZ*I_NAI_F3z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_S_b
   * RHS shell quartet name: SQ_NAI_H_S_b
   ************************************************************/
  Double I_NAI_H5x_Px_b = I_NAI_I6x_S_b+ABX*I_NAI_H5x_S_b;
  Double I_NAI_H4xy_Px_b = I_NAI_I5xy_S_b+ABX*I_NAI_H4xy_S_b;
  Double I_NAI_H4xz_Px_b = I_NAI_I5xz_S_b+ABX*I_NAI_H4xz_S_b;
  Double I_NAI_H3x2y_Px_b = I_NAI_I4x2y_S_b+ABX*I_NAI_H3x2y_S_b;
  Double I_NAI_H3xyz_Px_b = I_NAI_I4xyz_S_b+ABX*I_NAI_H3xyz_S_b;
  Double I_NAI_H3x2z_Px_b = I_NAI_I4x2z_S_b+ABX*I_NAI_H3x2z_S_b;
  Double I_NAI_H2x3y_Px_b = I_NAI_I3x3y_S_b+ABX*I_NAI_H2x3y_S_b;
  Double I_NAI_H2x2yz_Px_b = I_NAI_I3x2yz_S_b+ABX*I_NAI_H2x2yz_S_b;
  Double I_NAI_H2xy2z_Px_b = I_NAI_I3xy2z_S_b+ABX*I_NAI_H2xy2z_S_b;
  Double I_NAI_H2x3z_Px_b = I_NAI_I3x3z_S_b+ABX*I_NAI_H2x3z_S_b;
  Double I_NAI_Hx4y_Px_b = I_NAI_I2x4y_S_b+ABX*I_NAI_Hx4y_S_b;
  Double I_NAI_Hx3yz_Px_b = I_NAI_I2x3yz_S_b+ABX*I_NAI_Hx3yz_S_b;
  Double I_NAI_Hx2y2z_Px_b = I_NAI_I2x2y2z_S_b+ABX*I_NAI_Hx2y2z_S_b;
  Double I_NAI_Hxy3z_Px_b = I_NAI_I2xy3z_S_b+ABX*I_NAI_Hxy3z_S_b;
  Double I_NAI_Hx4z_Px_b = I_NAI_I2x4z_S_b+ABX*I_NAI_Hx4z_S_b;
  Double I_NAI_H5y_Px_b = I_NAI_Ix5y_S_b+ABX*I_NAI_H5y_S_b;
  Double I_NAI_H4yz_Px_b = I_NAI_Ix4yz_S_b+ABX*I_NAI_H4yz_S_b;
  Double I_NAI_H3y2z_Px_b = I_NAI_Ix3y2z_S_b+ABX*I_NAI_H3y2z_S_b;
  Double I_NAI_H2y3z_Px_b = I_NAI_Ix2y3z_S_b+ABX*I_NAI_H2y3z_S_b;
  Double I_NAI_Hy4z_Px_b = I_NAI_Ixy4z_S_b+ABX*I_NAI_Hy4z_S_b;
  Double I_NAI_H5z_Px_b = I_NAI_Ix5z_S_b+ABX*I_NAI_H5z_S_b;
  Double I_NAI_H4xy_Py_b = I_NAI_I4x2y_S_b+ABY*I_NAI_H4xy_S_b;
  Double I_NAI_H4xz_Py_b = I_NAI_I4xyz_S_b+ABY*I_NAI_H4xz_S_b;
  Double I_NAI_H3x2y_Py_b = I_NAI_I3x3y_S_b+ABY*I_NAI_H3x2y_S_b;
  Double I_NAI_H3xyz_Py_b = I_NAI_I3x2yz_S_b+ABY*I_NAI_H3xyz_S_b;
  Double I_NAI_H3x2z_Py_b = I_NAI_I3xy2z_S_b+ABY*I_NAI_H3x2z_S_b;
  Double I_NAI_H2x3y_Py_b = I_NAI_I2x4y_S_b+ABY*I_NAI_H2x3y_S_b;
  Double I_NAI_H2x2yz_Py_b = I_NAI_I2x3yz_S_b+ABY*I_NAI_H2x2yz_S_b;
  Double I_NAI_H2xy2z_Py_b = I_NAI_I2x2y2z_S_b+ABY*I_NAI_H2xy2z_S_b;
  Double I_NAI_H2x3z_Py_b = I_NAI_I2xy3z_S_b+ABY*I_NAI_H2x3z_S_b;
  Double I_NAI_Hx4y_Py_b = I_NAI_Ix5y_S_b+ABY*I_NAI_Hx4y_S_b;
  Double I_NAI_Hx3yz_Py_b = I_NAI_Ix4yz_S_b+ABY*I_NAI_Hx3yz_S_b;
  Double I_NAI_Hx2y2z_Py_b = I_NAI_Ix3y2z_S_b+ABY*I_NAI_Hx2y2z_S_b;
  Double I_NAI_Hxy3z_Py_b = I_NAI_Ix2y3z_S_b+ABY*I_NAI_Hxy3z_S_b;
  Double I_NAI_Hx4z_Py_b = I_NAI_Ixy4z_S_b+ABY*I_NAI_Hx4z_S_b;
  Double I_NAI_H5y_Py_b = I_NAI_I6y_S_b+ABY*I_NAI_H5y_S_b;
  Double I_NAI_H4yz_Py_b = I_NAI_I5yz_S_b+ABY*I_NAI_H4yz_S_b;
  Double I_NAI_H3y2z_Py_b = I_NAI_I4y2z_S_b+ABY*I_NAI_H3y2z_S_b;
  Double I_NAI_H2y3z_Py_b = I_NAI_I3y3z_S_b+ABY*I_NAI_H2y3z_S_b;
  Double I_NAI_Hy4z_Py_b = I_NAI_I2y4z_S_b+ABY*I_NAI_Hy4z_S_b;
  Double I_NAI_H5z_Py_b = I_NAI_Iy5z_S_b+ABY*I_NAI_H5z_S_b;
  Double I_NAI_H4xy_Pz_b = I_NAI_I4xyz_S_b+ABZ*I_NAI_H4xy_S_b;
  Double I_NAI_H4xz_Pz_b = I_NAI_I4x2z_S_b+ABZ*I_NAI_H4xz_S_b;
  Double I_NAI_H3x2y_Pz_b = I_NAI_I3x2yz_S_b+ABZ*I_NAI_H3x2y_S_b;
  Double I_NAI_H3xyz_Pz_b = I_NAI_I3xy2z_S_b+ABZ*I_NAI_H3xyz_S_b;
  Double I_NAI_H3x2z_Pz_b = I_NAI_I3x3z_S_b+ABZ*I_NAI_H3x2z_S_b;
  Double I_NAI_H2x3y_Pz_b = I_NAI_I2x3yz_S_b+ABZ*I_NAI_H2x3y_S_b;
  Double I_NAI_H2x2yz_Pz_b = I_NAI_I2x2y2z_S_b+ABZ*I_NAI_H2x2yz_S_b;
  Double I_NAI_H2xy2z_Pz_b = I_NAI_I2xy3z_S_b+ABZ*I_NAI_H2xy2z_S_b;
  Double I_NAI_H2x3z_Pz_b = I_NAI_I2x4z_S_b+ABZ*I_NAI_H2x3z_S_b;
  Double I_NAI_Hx4y_Pz_b = I_NAI_Ix4yz_S_b+ABZ*I_NAI_Hx4y_S_b;
  Double I_NAI_Hx3yz_Pz_b = I_NAI_Ix3y2z_S_b+ABZ*I_NAI_Hx3yz_S_b;
  Double I_NAI_Hx2y2z_Pz_b = I_NAI_Ix2y3z_S_b+ABZ*I_NAI_Hx2y2z_S_b;
  Double I_NAI_Hxy3z_Pz_b = I_NAI_Ixy4z_S_b+ABZ*I_NAI_Hxy3z_S_b;
  Double I_NAI_Hx4z_Pz_b = I_NAI_Ix5z_S_b+ABZ*I_NAI_Hx4z_S_b;
  Double I_NAI_H4yz_Pz_b = I_NAI_I4y2z_S_b+ABZ*I_NAI_H4yz_S_b;
  Double I_NAI_H3y2z_Pz_b = I_NAI_I3y3z_S_b+ABZ*I_NAI_H3y2z_S_b;
  Double I_NAI_H2y3z_Pz_b = I_NAI_I2y4z_S_b+ABZ*I_NAI_H2y3z_S_b;
  Double I_NAI_Hy4z_Pz_b = I_NAI_Iy5z_S_b+ABZ*I_NAI_Hy4z_S_b;
  Double I_NAI_H5z_Pz_b = I_NAI_I6z_S_b+ABZ*I_NAI_H5z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 45 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_P_b
   * RHS shell quartet name: SQ_NAI_G_P_b
   ************************************************************/
  Double I_NAI_G4x_D2x_b = I_NAI_H5x_Px_b+ABX*I_NAI_G4x_Px_b;
  Double I_NAI_G3xy_D2x_b = I_NAI_H4xy_Px_b+ABX*I_NAI_G3xy_Px_b;
  Double I_NAI_G3xz_D2x_b = I_NAI_H4xz_Px_b+ABX*I_NAI_G3xz_Px_b;
  Double I_NAI_G2x2y_D2x_b = I_NAI_H3x2y_Px_b+ABX*I_NAI_G2x2y_Px_b;
  Double I_NAI_G2xyz_D2x_b = I_NAI_H3xyz_Px_b+ABX*I_NAI_G2xyz_Px_b;
  Double I_NAI_G2x2z_D2x_b = I_NAI_H3x2z_Px_b+ABX*I_NAI_G2x2z_Px_b;
  Double I_NAI_Gx3y_D2x_b = I_NAI_H2x3y_Px_b+ABX*I_NAI_Gx3y_Px_b;
  Double I_NAI_Gx2yz_D2x_b = I_NAI_H2x2yz_Px_b+ABX*I_NAI_Gx2yz_Px_b;
  Double I_NAI_Gxy2z_D2x_b = I_NAI_H2xy2z_Px_b+ABX*I_NAI_Gxy2z_Px_b;
  Double I_NAI_Gx3z_D2x_b = I_NAI_H2x3z_Px_b+ABX*I_NAI_Gx3z_Px_b;
  Double I_NAI_G4y_D2x_b = I_NAI_Hx4y_Px_b+ABX*I_NAI_G4y_Px_b;
  Double I_NAI_G3yz_D2x_b = I_NAI_Hx3yz_Px_b+ABX*I_NAI_G3yz_Px_b;
  Double I_NAI_G2y2z_D2x_b = I_NAI_Hx2y2z_Px_b+ABX*I_NAI_G2y2z_Px_b;
  Double I_NAI_Gy3z_D2x_b = I_NAI_Hxy3z_Px_b+ABX*I_NAI_Gy3z_Px_b;
  Double I_NAI_G4z_D2x_b = I_NAI_Hx4z_Px_b+ABX*I_NAI_G4z_Px_b;
  Double I_NAI_G4x_D2y_b = I_NAI_H4xy_Py_b+ABY*I_NAI_G4x_Py_b;
  Double I_NAI_G3xy_D2y_b = I_NAI_H3x2y_Py_b+ABY*I_NAI_G3xy_Py_b;
  Double I_NAI_G3xz_D2y_b = I_NAI_H3xyz_Py_b+ABY*I_NAI_G3xz_Py_b;
  Double I_NAI_G2x2y_D2y_b = I_NAI_H2x3y_Py_b+ABY*I_NAI_G2x2y_Py_b;
  Double I_NAI_G2xyz_D2y_b = I_NAI_H2x2yz_Py_b+ABY*I_NAI_G2xyz_Py_b;
  Double I_NAI_G2x2z_D2y_b = I_NAI_H2xy2z_Py_b+ABY*I_NAI_G2x2z_Py_b;
  Double I_NAI_Gx3y_D2y_b = I_NAI_Hx4y_Py_b+ABY*I_NAI_Gx3y_Py_b;
  Double I_NAI_Gx2yz_D2y_b = I_NAI_Hx3yz_Py_b+ABY*I_NAI_Gx2yz_Py_b;
  Double I_NAI_Gxy2z_D2y_b = I_NAI_Hx2y2z_Py_b+ABY*I_NAI_Gxy2z_Py_b;
  Double I_NAI_Gx3z_D2y_b = I_NAI_Hxy3z_Py_b+ABY*I_NAI_Gx3z_Py_b;
  Double I_NAI_G4y_D2y_b = I_NAI_H5y_Py_b+ABY*I_NAI_G4y_Py_b;
  Double I_NAI_G3yz_D2y_b = I_NAI_H4yz_Py_b+ABY*I_NAI_G3yz_Py_b;
  Double I_NAI_G2y2z_D2y_b = I_NAI_H3y2z_Py_b+ABY*I_NAI_G2y2z_Py_b;
  Double I_NAI_Gy3z_D2y_b = I_NAI_H2y3z_Py_b+ABY*I_NAI_Gy3z_Py_b;
  Double I_NAI_G4z_D2y_b = I_NAI_Hy4z_Py_b+ABY*I_NAI_G4z_Py_b;
  Double I_NAI_G4x_D2z_b = I_NAI_H4xz_Pz_b+ABZ*I_NAI_G4x_Pz_b;
  Double I_NAI_G3xy_D2z_b = I_NAI_H3xyz_Pz_b+ABZ*I_NAI_G3xy_Pz_b;
  Double I_NAI_G3xz_D2z_b = I_NAI_H3x2z_Pz_b+ABZ*I_NAI_G3xz_Pz_b;
  Double I_NAI_G2x2y_D2z_b = I_NAI_H2x2yz_Pz_b+ABZ*I_NAI_G2x2y_Pz_b;
  Double I_NAI_G2xyz_D2z_b = I_NAI_H2xy2z_Pz_b+ABZ*I_NAI_G2xyz_Pz_b;
  Double I_NAI_G2x2z_D2z_b = I_NAI_H2x3z_Pz_b+ABZ*I_NAI_G2x2z_Pz_b;
  Double I_NAI_Gx3y_D2z_b = I_NAI_Hx3yz_Pz_b+ABZ*I_NAI_Gx3y_Pz_b;
  Double I_NAI_Gx2yz_D2z_b = I_NAI_Hx2y2z_Pz_b+ABZ*I_NAI_Gx2yz_Pz_b;
  Double I_NAI_Gxy2z_D2z_b = I_NAI_Hxy3z_Pz_b+ABZ*I_NAI_Gxy2z_Pz_b;
  Double I_NAI_Gx3z_D2z_b = I_NAI_Hx4z_Pz_b+ABZ*I_NAI_Gx3z_Pz_b;
  Double I_NAI_G4y_D2z_b = I_NAI_H4yz_Pz_b+ABZ*I_NAI_G4y_Pz_b;
  Double I_NAI_G3yz_D2z_b = I_NAI_H3y2z_Pz_b+ABZ*I_NAI_G3yz_Pz_b;
  Double I_NAI_G2y2z_D2z_b = I_NAI_H2y3z_Pz_b+ABZ*I_NAI_G2y2z_Pz_b;
  Double I_NAI_Gy3z_D2z_b = I_NAI_Hy4z_Pz_b+ABZ*I_NAI_Gy3z_Pz_b;
  Double I_NAI_G4z_D2z_b = I_NAI_H5z_Pz_b+ABZ*I_NAI_G4z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_D_b
   * RHS shell quartet name: SQ_NAI_F_D_b
   ************************************************************/
  Double I_NAI_F3x_F3x_b = I_NAI_G4x_D2x_b+ABX*I_NAI_F3x_D2x_b;
  Double I_NAI_F2xy_F3x_b = I_NAI_G3xy_D2x_b+ABX*I_NAI_F2xy_D2x_b;
  Double I_NAI_F2xz_F3x_b = I_NAI_G3xz_D2x_b+ABX*I_NAI_F2xz_D2x_b;
  Double I_NAI_Fx2y_F3x_b = I_NAI_G2x2y_D2x_b+ABX*I_NAI_Fx2y_D2x_b;
  Double I_NAI_Fxyz_F3x_b = I_NAI_G2xyz_D2x_b+ABX*I_NAI_Fxyz_D2x_b;
  Double I_NAI_Fx2z_F3x_b = I_NAI_G2x2z_D2x_b+ABX*I_NAI_Fx2z_D2x_b;
  Double I_NAI_F3y_F3x_b = I_NAI_Gx3y_D2x_b+ABX*I_NAI_F3y_D2x_b;
  Double I_NAI_F2yz_F3x_b = I_NAI_Gx2yz_D2x_b+ABX*I_NAI_F2yz_D2x_b;
  Double I_NAI_Fy2z_F3x_b = I_NAI_Gxy2z_D2x_b+ABX*I_NAI_Fy2z_D2x_b;
  Double I_NAI_F3z_F3x_b = I_NAI_Gx3z_D2x_b+ABX*I_NAI_F3z_D2x_b;
  Double I_NAI_F3x_F2xy_b = I_NAI_G3xy_D2x_b+ABY*I_NAI_F3x_D2x_b;
  Double I_NAI_F2xy_F2xy_b = I_NAI_G2x2y_D2x_b+ABY*I_NAI_F2xy_D2x_b;
  Double I_NAI_F2xz_F2xy_b = I_NAI_G2xyz_D2x_b+ABY*I_NAI_F2xz_D2x_b;
  Double I_NAI_Fx2y_F2xy_b = I_NAI_Gx3y_D2x_b+ABY*I_NAI_Fx2y_D2x_b;
  Double I_NAI_Fxyz_F2xy_b = I_NAI_Gx2yz_D2x_b+ABY*I_NAI_Fxyz_D2x_b;
  Double I_NAI_Fx2z_F2xy_b = I_NAI_Gxy2z_D2x_b+ABY*I_NAI_Fx2z_D2x_b;
  Double I_NAI_F3y_F2xy_b = I_NAI_G4y_D2x_b+ABY*I_NAI_F3y_D2x_b;
  Double I_NAI_F2yz_F2xy_b = I_NAI_G3yz_D2x_b+ABY*I_NAI_F2yz_D2x_b;
  Double I_NAI_Fy2z_F2xy_b = I_NAI_G2y2z_D2x_b+ABY*I_NAI_Fy2z_D2x_b;
  Double I_NAI_F3z_F2xy_b = I_NAI_Gy3z_D2x_b+ABY*I_NAI_F3z_D2x_b;
  Double I_NAI_F3x_F2xz_b = I_NAI_G3xz_D2x_b+ABZ*I_NAI_F3x_D2x_b;
  Double I_NAI_F2xy_F2xz_b = I_NAI_G2xyz_D2x_b+ABZ*I_NAI_F2xy_D2x_b;
  Double I_NAI_F2xz_F2xz_b = I_NAI_G2x2z_D2x_b+ABZ*I_NAI_F2xz_D2x_b;
  Double I_NAI_Fx2y_F2xz_b = I_NAI_Gx2yz_D2x_b+ABZ*I_NAI_Fx2y_D2x_b;
  Double I_NAI_Fxyz_F2xz_b = I_NAI_Gxy2z_D2x_b+ABZ*I_NAI_Fxyz_D2x_b;
  Double I_NAI_Fx2z_F2xz_b = I_NAI_Gx3z_D2x_b+ABZ*I_NAI_Fx2z_D2x_b;
  Double I_NAI_F3y_F2xz_b = I_NAI_G3yz_D2x_b+ABZ*I_NAI_F3y_D2x_b;
  Double I_NAI_F2yz_F2xz_b = I_NAI_G2y2z_D2x_b+ABZ*I_NAI_F2yz_D2x_b;
  Double I_NAI_Fy2z_F2xz_b = I_NAI_Gy3z_D2x_b+ABZ*I_NAI_Fy2z_D2x_b;
  Double I_NAI_F3z_F2xz_b = I_NAI_G4z_D2x_b+ABZ*I_NAI_F3z_D2x_b;
  Double I_NAI_F3x_Fx2y_b = I_NAI_G4x_D2y_b+ABX*I_NAI_F3x_D2y_b;
  Double I_NAI_F2xy_Fx2y_b = I_NAI_G3xy_D2y_b+ABX*I_NAI_F2xy_D2y_b;
  Double I_NAI_F2xz_Fx2y_b = I_NAI_G3xz_D2y_b+ABX*I_NAI_F2xz_D2y_b;
  Double I_NAI_Fx2y_Fx2y_b = I_NAI_G2x2y_D2y_b+ABX*I_NAI_Fx2y_D2y_b;
  Double I_NAI_Fxyz_Fx2y_b = I_NAI_G2xyz_D2y_b+ABX*I_NAI_Fxyz_D2y_b;
  Double I_NAI_Fx2z_Fx2y_b = I_NAI_G2x2z_D2y_b+ABX*I_NAI_Fx2z_D2y_b;
  Double I_NAI_F3y_Fx2y_b = I_NAI_Gx3y_D2y_b+ABX*I_NAI_F3y_D2y_b;
  Double I_NAI_F2yz_Fx2y_b = I_NAI_Gx2yz_D2y_b+ABX*I_NAI_F2yz_D2y_b;
  Double I_NAI_Fy2z_Fx2y_b = I_NAI_Gxy2z_D2y_b+ABX*I_NAI_Fy2z_D2y_b;
  Double I_NAI_F3z_Fx2y_b = I_NAI_Gx3z_D2y_b+ABX*I_NAI_F3z_D2y_b;
  Double I_NAI_F3x_Fx2z_b = I_NAI_G4x_D2z_b+ABX*I_NAI_F3x_D2z_b;
  Double I_NAI_F2xy_Fx2z_b = I_NAI_G3xy_D2z_b+ABX*I_NAI_F2xy_D2z_b;
  Double I_NAI_F2xz_Fx2z_b = I_NAI_G3xz_D2z_b+ABX*I_NAI_F2xz_D2z_b;
  Double I_NAI_Fx2y_Fx2z_b = I_NAI_G2x2y_D2z_b+ABX*I_NAI_Fx2y_D2z_b;
  Double I_NAI_Fxyz_Fx2z_b = I_NAI_G2xyz_D2z_b+ABX*I_NAI_Fxyz_D2z_b;
  Double I_NAI_Fx2z_Fx2z_b = I_NAI_G2x2z_D2z_b+ABX*I_NAI_Fx2z_D2z_b;
  Double I_NAI_F3y_Fx2z_b = I_NAI_Gx3y_D2z_b+ABX*I_NAI_F3y_D2z_b;
  Double I_NAI_F2yz_Fx2z_b = I_NAI_Gx2yz_D2z_b+ABX*I_NAI_F2yz_D2z_b;
  Double I_NAI_Fy2z_Fx2z_b = I_NAI_Gxy2z_D2z_b+ABX*I_NAI_Fy2z_D2z_b;
  Double I_NAI_F3z_Fx2z_b = I_NAI_Gx3z_D2z_b+ABX*I_NAI_F3z_D2z_b;
  Double I_NAI_F3x_F3y_b = I_NAI_G3xy_D2y_b+ABY*I_NAI_F3x_D2y_b;
  Double I_NAI_F2xy_F3y_b = I_NAI_G2x2y_D2y_b+ABY*I_NAI_F2xy_D2y_b;
  Double I_NAI_F2xz_F3y_b = I_NAI_G2xyz_D2y_b+ABY*I_NAI_F2xz_D2y_b;
  Double I_NAI_Fx2y_F3y_b = I_NAI_Gx3y_D2y_b+ABY*I_NAI_Fx2y_D2y_b;
  Double I_NAI_Fxyz_F3y_b = I_NAI_Gx2yz_D2y_b+ABY*I_NAI_Fxyz_D2y_b;
  Double I_NAI_Fx2z_F3y_b = I_NAI_Gxy2z_D2y_b+ABY*I_NAI_Fx2z_D2y_b;
  Double I_NAI_F3y_F3y_b = I_NAI_G4y_D2y_b+ABY*I_NAI_F3y_D2y_b;
  Double I_NAI_F2yz_F3y_b = I_NAI_G3yz_D2y_b+ABY*I_NAI_F2yz_D2y_b;
  Double I_NAI_Fy2z_F3y_b = I_NAI_G2y2z_D2y_b+ABY*I_NAI_Fy2z_D2y_b;
  Double I_NAI_F3z_F3y_b = I_NAI_Gy3z_D2y_b+ABY*I_NAI_F3z_D2y_b;
  Double I_NAI_F3x_F2yz_b = I_NAI_G3xz_D2y_b+ABZ*I_NAI_F3x_D2y_b;
  Double I_NAI_F2xy_F2yz_b = I_NAI_G2xyz_D2y_b+ABZ*I_NAI_F2xy_D2y_b;
  Double I_NAI_F2xz_F2yz_b = I_NAI_G2x2z_D2y_b+ABZ*I_NAI_F2xz_D2y_b;
  Double I_NAI_Fx2y_F2yz_b = I_NAI_Gx2yz_D2y_b+ABZ*I_NAI_Fx2y_D2y_b;
  Double I_NAI_Fxyz_F2yz_b = I_NAI_Gxy2z_D2y_b+ABZ*I_NAI_Fxyz_D2y_b;
  Double I_NAI_Fx2z_F2yz_b = I_NAI_Gx3z_D2y_b+ABZ*I_NAI_Fx2z_D2y_b;
  Double I_NAI_F3y_F2yz_b = I_NAI_G3yz_D2y_b+ABZ*I_NAI_F3y_D2y_b;
  Double I_NAI_F2yz_F2yz_b = I_NAI_G2y2z_D2y_b+ABZ*I_NAI_F2yz_D2y_b;
  Double I_NAI_Fy2z_F2yz_b = I_NAI_Gy3z_D2y_b+ABZ*I_NAI_Fy2z_D2y_b;
  Double I_NAI_F3z_F2yz_b = I_NAI_G4z_D2y_b+ABZ*I_NAI_F3z_D2y_b;
  Double I_NAI_F3x_F3z_b = I_NAI_G3xz_D2z_b+ABZ*I_NAI_F3x_D2z_b;
  Double I_NAI_F2xy_F3z_b = I_NAI_G2xyz_D2z_b+ABZ*I_NAI_F2xy_D2z_b;
  Double I_NAI_F2xz_F3z_b = I_NAI_G2x2z_D2z_b+ABZ*I_NAI_F2xz_D2z_b;
  Double I_NAI_Fx2y_F3z_b = I_NAI_Gx2yz_D2z_b+ABZ*I_NAI_Fx2y_D2z_b;
  Double I_NAI_Fxyz_F3z_b = I_NAI_Gxy2z_D2z_b+ABZ*I_NAI_Fxyz_D2z_b;
  Double I_NAI_Fx2z_F3z_b = I_NAI_Gx3z_D2z_b+ABZ*I_NAI_Fx2z_D2z_b;
  Double I_NAI_F3y_F3z_b = I_NAI_G3yz_D2z_b+ABZ*I_NAI_F3y_D2z_b;
  Double I_NAI_F2yz_F3z_b = I_NAI_G2y2z_D2z_b+ABZ*I_NAI_F2yz_D2z_b;
  Double I_NAI_Fy2z_F3z_b = I_NAI_Gy3z_D2z_b+ABZ*I_NAI_Fy2z_D2z_b;
  Double I_NAI_F3z_F3z_b = I_NAI_G4z_D2z_b+ABZ*I_NAI_F3z_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_I_P_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 24 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_K_S_b
   * RHS shell quartet name: SQ_NAI_I_S_b
   ************************************************************/
  Double I_NAI_I6x_Px_b = I_NAI_K7x_S_b+ABX*I_NAI_I6x_S_b;
  Double I_NAI_I5xy_Px_b = I_NAI_K6xy_S_b+ABX*I_NAI_I5xy_S_b;
  Double I_NAI_I5xz_Px_b = I_NAI_K6xz_S_b+ABX*I_NAI_I5xz_S_b;
  Double I_NAI_I4x2y_Px_b = I_NAI_K5x2y_S_b+ABX*I_NAI_I4x2y_S_b;
  Double I_NAI_I4xyz_Px_b = I_NAI_K5xyz_S_b+ABX*I_NAI_I4xyz_S_b;
  Double I_NAI_I4x2z_Px_b = I_NAI_K5x2z_S_b+ABX*I_NAI_I4x2z_S_b;
  Double I_NAI_I3x3y_Px_b = I_NAI_K4x3y_S_b+ABX*I_NAI_I3x3y_S_b;
  Double I_NAI_I3x2yz_Px_b = I_NAI_K4x2yz_S_b+ABX*I_NAI_I3x2yz_S_b;
  Double I_NAI_I3xy2z_Px_b = I_NAI_K4xy2z_S_b+ABX*I_NAI_I3xy2z_S_b;
  Double I_NAI_I3x3z_Px_b = I_NAI_K4x3z_S_b+ABX*I_NAI_I3x3z_S_b;
  Double I_NAI_I2x4y_Px_b = I_NAI_K3x4y_S_b+ABX*I_NAI_I2x4y_S_b;
  Double I_NAI_I2x3yz_Px_b = I_NAI_K3x3yz_S_b+ABX*I_NAI_I2x3yz_S_b;
  Double I_NAI_I2x2y2z_Px_b = I_NAI_K3x2y2z_S_b+ABX*I_NAI_I2x2y2z_S_b;
  Double I_NAI_I2xy3z_Px_b = I_NAI_K3xy3z_S_b+ABX*I_NAI_I2xy3z_S_b;
  Double I_NAI_I2x4z_Px_b = I_NAI_K3x4z_S_b+ABX*I_NAI_I2x4z_S_b;
  Double I_NAI_Ix5y_Px_b = I_NAI_K2x5y_S_b+ABX*I_NAI_Ix5y_S_b;
  Double I_NAI_Ix4yz_Px_b = I_NAI_K2x4yz_S_b+ABX*I_NAI_Ix4yz_S_b;
  Double I_NAI_Ix3y2z_Px_b = I_NAI_K2x3y2z_S_b+ABX*I_NAI_Ix3y2z_S_b;
  Double I_NAI_Ix2y3z_Px_b = I_NAI_K2x2y3z_S_b+ABX*I_NAI_Ix2y3z_S_b;
  Double I_NAI_Ixy4z_Px_b = I_NAI_K2xy4z_S_b+ABX*I_NAI_Ixy4z_S_b;
  Double I_NAI_Ix5z_Px_b = I_NAI_K2x5z_S_b+ABX*I_NAI_Ix5z_S_b;
  Double I_NAI_I4x2y_Py_b = I_NAI_K4x3y_S_b+ABY*I_NAI_I4x2y_S_b;
  Double I_NAI_I4xyz_Py_b = I_NAI_K4x2yz_S_b+ABY*I_NAI_I4xyz_S_b;
  Double I_NAI_I3x3y_Py_b = I_NAI_K3x4y_S_b+ABY*I_NAI_I3x3y_S_b;
  Double I_NAI_I3x2yz_Py_b = I_NAI_K3x3yz_S_b+ABY*I_NAI_I3x2yz_S_b;
  Double I_NAI_I3xy2z_Py_b = I_NAI_K3x2y2z_S_b+ABY*I_NAI_I3xy2z_S_b;
  Double I_NAI_I2x4y_Py_b = I_NAI_K2x5y_S_b+ABY*I_NAI_I2x4y_S_b;
  Double I_NAI_I2x3yz_Py_b = I_NAI_K2x4yz_S_b+ABY*I_NAI_I2x3yz_S_b;
  Double I_NAI_I2x2y2z_Py_b = I_NAI_K2x3y2z_S_b+ABY*I_NAI_I2x2y2z_S_b;
  Double I_NAI_I2xy3z_Py_b = I_NAI_K2x2y3z_S_b+ABY*I_NAI_I2xy3z_S_b;
  Double I_NAI_Ix5y_Py_b = I_NAI_Kx6y_S_b+ABY*I_NAI_Ix5y_S_b;
  Double I_NAI_Ix4yz_Py_b = I_NAI_Kx5yz_S_b+ABY*I_NAI_Ix4yz_S_b;
  Double I_NAI_Ix3y2z_Py_b = I_NAI_Kx4y2z_S_b+ABY*I_NAI_Ix3y2z_S_b;
  Double I_NAI_Ix2y3z_Py_b = I_NAI_Kx3y3z_S_b+ABY*I_NAI_Ix2y3z_S_b;
  Double I_NAI_Ixy4z_Py_b = I_NAI_Kx2y4z_S_b+ABY*I_NAI_Ixy4z_S_b;
  Double I_NAI_I6y_Py_b = I_NAI_K7y_S_b+ABY*I_NAI_I6y_S_b;
  Double I_NAI_I5yz_Py_b = I_NAI_K6yz_S_b+ABY*I_NAI_I5yz_S_b;
  Double I_NAI_I4y2z_Py_b = I_NAI_K5y2z_S_b+ABY*I_NAI_I4y2z_S_b;
  Double I_NAI_I3y3z_Py_b = I_NAI_K4y3z_S_b+ABY*I_NAI_I3y3z_S_b;
  Double I_NAI_I2y4z_Py_b = I_NAI_K3y4z_S_b+ABY*I_NAI_I2y4z_S_b;
  Double I_NAI_Iy5z_Py_b = I_NAI_K2y5z_S_b+ABY*I_NAI_Iy5z_S_b;
  Double I_NAI_I4xyz_Pz_b = I_NAI_K4xy2z_S_b+ABZ*I_NAI_I4xyz_S_b;
  Double I_NAI_I4x2z_Pz_b = I_NAI_K4x3z_S_b+ABZ*I_NAI_I4x2z_S_b;
  Double I_NAI_I3x2yz_Pz_b = I_NAI_K3x2y2z_S_b+ABZ*I_NAI_I3x2yz_S_b;
  Double I_NAI_I3xy2z_Pz_b = I_NAI_K3xy3z_S_b+ABZ*I_NAI_I3xy2z_S_b;
  Double I_NAI_I3x3z_Pz_b = I_NAI_K3x4z_S_b+ABZ*I_NAI_I3x3z_S_b;
  Double I_NAI_I2x3yz_Pz_b = I_NAI_K2x3y2z_S_b+ABZ*I_NAI_I2x3yz_S_b;
  Double I_NAI_I2x2y2z_Pz_b = I_NAI_K2x2y3z_S_b+ABZ*I_NAI_I2x2y2z_S_b;
  Double I_NAI_I2xy3z_Pz_b = I_NAI_K2xy4z_S_b+ABZ*I_NAI_I2xy3z_S_b;
  Double I_NAI_I2x4z_Pz_b = I_NAI_K2x5z_S_b+ABZ*I_NAI_I2x4z_S_b;
  Double I_NAI_Ix4yz_Pz_b = I_NAI_Kx4y2z_S_b+ABZ*I_NAI_Ix4yz_S_b;
  Double I_NAI_Ix3y2z_Pz_b = I_NAI_Kx3y3z_S_b+ABZ*I_NAI_Ix3y2z_S_b;
  Double I_NAI_Ix2y3z_Pz_b = I_NAI_Kx2y4z_S_b+ABZ*I_NAI_Ix2y3z_S_b;
  Double I_NAI_Ixy4z_Pz_b = I_NAI_Kxy5z_S_b+ABZ*I_NAI_Ixy4z_S_b;
  Double I_NAI_Ix5z_Pz_b = I_NAI_Kx6z_S_b+ABZ*I_NAI_Ix5z_S_b;
  Double I_NAI_I4y2z_Pz_b = I_NAI_K4y3z_S_b+ABZ*I_NAI_I4y2z_S_b;
  Double I_NAI_I3y3z_Pz_b = I_NAI_K3y4z_S_b+ABZ*I_NAI_I3y3z_S_b;
  Double I_NAI_I2y4z_Pz_b = I_NAI_K2y5z_S_b+ABZ*I_NAI_I2y4z_S_b;
  Double I_NAI_Iy5z_Pz_b = I_NAI_Ky6z_S_b+ABZ*I_NAI_Iy5z_S_b;
  Double I_NAI_I6z_Pz_b = I_NAI_K7z_S_b+ABZ*I_NAI_I6z_S_b;

  /************************************************************
   * shell quartet name: SQ_NAI_H_D_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 66 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_I_P_b
   * RHS shell quartet name: SQ_NAI_H_P_b
   ************************************************************/
  Double I_NAI_H5x_D2x_b = I_NAI_I6x_Px_b+ABX*I_NAI_H5x_Px_b;
  Double I_NAI_H4xy_D2x_b = I_NAI_I5xy_Px_b+ABX*I_NAI_H4xy_Px_b;
  Double I_NAI_H4xz_D2x_b = I_NAI_I5xz_Px_b+ABX*I_NAI_H4xz_Px_b;
  Double I_NAI_H3x2y_D2x_b = I_NAI_I4x2y_Px_b+ABX*I_NAI_H3x2y_Px_b;
  Double I_NAI_H3xyz_D2x_b = I_NAI_I4xyz_Px_b+ABX*I_NAI_H3xyz_Px_b;
  Double I_NAI_H3x2z_D2x_b = I_NAI_I4x2z_Px_b+ABX*I_NAI_H3x2z_Px_b;
  Double I_NAI_H2x3y_D2x_b = I_NAI_I3x3y_Px_b+ABX*I_NAI_H2x3y_Px_b;
  Double I_NAI_H2x2yz_D2x_b = I_NAI_I3x2yz_Px_b+ABX*I_NAI_H2x2yz_Px_b;
  Double I_NAI_H2xy2z_D2x_b = I_NAI_I3xy2z_Px_b+ABX*I_NAI_H2xy2z_Px_b;
  Double I_NAI_H2x3z_D2x_b = I_NAI_I3x3z_Px_b+ABX*I_NAI_H2x3z_Px_b;
  Double I_NAI_Hx4y_D2x_b = I_NAI_I2x4y_Px_b+ABX*I_NAI_Hx4y_Px_b;
  Double I_NAI_Hx3yz_D2x_b = I_NAI_I2x3yz_Px_b+ABX*I_NAI_Hx3yz_Px_b;
  Double I_NAI_Hx2y2z_D2x_b = I_NAI_I2x2y2z_Px_b+ABX*I_NAI_Hx2y2z_Px_b;
  Double I_NAI_Hxy3z_D2x_b = I_NAI_I2xy3z_Px_b+ABX*I_NAI_Hxy3z_Px_b;
  Double I_NAI_Hx4z_D2x_b = I_NAI_I2x4z_Px_b+ABX*I_NAI_Hx4z_Px_b;
  Double I_NAI_H5y_D2x_b = I_NAI_Ix5y_Px_b+ABX*I_NAI_H5y_Px_b;
  Double I_NAI_H4yz_D2x_b = I_NAI_Ix4yz_Px_b+ABX*I_NAI_H4yz_Px_b;
  Double I_NAI_H3y2z_D2x_b = I_NAI_Ix3y2z_Px_b+ABX*I_NAI_H3y2z_Px_b;
  Double I_NAI_H2y3z_D2x_b = I_NAI_Ix2y3z_Px_b+ABX*I_NAI_H2y3z_Px_b;
  Double I_NAI_Hy4z_D2x_b = I_NAI_Ixy4z_Px_b+ABX*I_NAI_Hy4z_Px_b;
  Double I_NAI_H5z_D2x_b = I_NAI_Ix5z_Px_b+ABX*I_NAI_H5z_Px_b;
  Double I_NAI_H4xy_D2y_b = I_NAI_I4x2y_Py_b+ABY*I_NAI_H4xy_Py_b;
  Double I_NAI_H4xz_D2y_b = I_NAI_I4xyz_Py_b+ABY*I_NAI_H4xz_Py_b;
  Double I_NAI_H3x2y_D2y_b = I_NAI_I3x3y_Py_b+ABY*I_NAI_H3x2y_Py_b;
  Double I_NAI_H3xyz_D2y_b = I_NAI_I3x2yz_Py_b+ABY*I_NAI_H3xyz_Py_b;
  Double I_NAI_H3x2z_D2y_b = I_NAI_I3xy2z_Py_b+ABY*I_NAI_H3x2z_Py_b;
  Double I_NAI_H2x3y_D2y_b = I_NAI_I2x4y_Py_b+ABY*I_NAI_H2x3y_Py_b;
  Double I_NAI_H2x2yz_D2y_b = I_NAI_I2x3yz_Py_b+ABY*I_NAI_H2x2yz_Py_b;
  Double I_NAI_H2xy2z_D2y_b = I_NAI_I2x2y2z_Py_b+ABY*I_NAI_H2xy2z_Py_b;
  Double I_NAI_H2x3z_D2y_b = I_NAI_I2xy3z_Py_b+ABY*I_NAI_H2x3z_Py_b;
  Double I_NAI_Hx4y_D2y_b = I_NAI_Ix5y_Py_b+ABY*I_NAI_Hx4y_Py_b;
  Double I_NAI_Hx3yz_D2y_b = I_NAI_Ix4yz_Py_b+ABY*I_NAI_Hx3yz_Py_b;
  Double I_NAI_Hx2y2z_D2y_b = I_NAI_Ix3y2z_Py_b+ABY*I_NAI_Hx2y2z_Py_b;
  Double I_NAI_Hxy3z_D2y_b = I_NAI_Ix2y3z_Py_b+ABY*I_NAI_Hxy3z_Py_b;
  Double I_NAI_Hx4z_D2y_b = I_NAI_Ixy4z_Py_b+ABY*I_NAI_Hx4z_Py_b;
  Double I_NAI_H5y_D2y_b = I_NAI_I6y_Py_b+ABY*I_NAI_H5y_Py_b;
  Double I_NAI_H4yz_D2y_b = I_NAI_I5yz_Py_b+ABY*I_NAI_H4yz_Py_b;
  Double I_NAI_H3y2z_D2y_b = I_NAI_I4y2z_Py_b+ABY*I_NAI_H3y2z_Py_b;
  Double I_NAI_H2y3z_D2y_b = I_NAI_I3y3z_Py_b+ABY*I_NAI_H2y3z_Py_b;
  Double I_NAI_Hy4z_D2y_b = I_NAI_I2y4z_Py_b+ABY*I_NAI_Hy4z_Py_b;
  Double I_NAI_H5z_D2y_b = I_NAI_Iy5z_Py_b+ABY*I_NAI_H5z_Py_b;
  Double I_NAI_H4xy_D2z_b = I_NAI_I4xyz_Pz_b+ABZ*I_NAI_H4xy_Pz_b;
  Double I_NAI_H4xz_D2z_b = I_NAI_I4x2z_Pz_b+ABZ*I_NAI_H4xz_Pz_b;
  Double I_NAI_H3x2y_D2z_b = I_NAI_I3x2yz_Pz_b+ABZ*I_NAI_H3x2y_Pz_b;
  Double I_NAI_H3xyz_D2z_b = I_NAI_I3xy2z_Pz_b+ABZ*I_NAI_H3xyz_Pz_b;
  Double I_NAI_H3x2z_D2z_b = I_NAI_I3x3z_Pz_b+ABZ*I_NAI_H3x2z_Pz_b;
  Double I_NAI_H2x3y_D2z_b = I_NAI_I2x3yz_Pz_b+ABZ*I_NAI_H2x3y_Pz_b;
  Double I_NAI_H2x2yz_D2z_b = I_NAI_I2x2y2z_Pz_b+ABZ*I_NAI_H2x2yz_Pz_b;
  Double I_NAI_H2xy2z_D2z_b = I_NAI_I2xy3z_Pz_b+ABZ*I_NAI_H2xy2z_Pz_b;
  Double I_NAI_H2x3z_D2z_b = I_NAI_I2x4z_Pz_b+ABZ*I_NAI_H2x3z_Pz_b;
  Double I_NAI_Hx4y_D2z_b = I_NAI_Ix4yz_Pz_b+ABZ*I_NAI_Hx4y_Pz_b;
  Double I_NAI_Hx3yz_D2z_b = I_NAI_Ix3y2z_Pz_b+ABZ*I_NAI_Hx3yz_Pz_b;
  Double I_NAI_Hx2y2z_D2z_b = I_NAI_Ix2y3z_Pz_b+ABZ*I_NAI_Hx2y2z_Pz_b;
  Double I_NAI_Hxy3z_D2z_b = I_NAI_Ixy4z_Pz_b+ABZ*I_NAI_Hxy3z_Pz_b;
  Double I_NAI_Hx4z_D2z_b = I_NAI_Ix5z_Pz_b+ABZ*I_NAI_Hx4z_Pz_b;
  Double I_NAI_H4yz_D2z_b = I_NAI_I4y2z_Pz_b+ABZ*I_NAI_H4yz_Pz_b;
  Double I_NAI_H3y2z_D2z_b = I_NAI_I3y3z_Pz_b+ABZ*I_NAI_H3y2z_Pz_b;
  Double I_NAI_H2y3z_D2z_b = I_NAI_I2y4z_Pz_b+ABZ*I_NAI_H2y3z_Pz_b;
  Double I_NAI_Hy4z_D2z_b = I_NAI_Iy5z_Pz_b+ABZ*I_NAI_Hy4z_Pz_b;
  Double I_NAI_H5z_D2z_b = I_NAI_I6z_Pz_b+ABZ*I_NAI_H5z_Pz_b;

  /************************************************************
   * shell quartet name: SQ_NAI_G_F_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 51 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_H_D_b
   * RHS shell quartet name: SQ_NAI_G_D_b
   ************************************************************/
  Double I_NAI_G4x_F3x_b = I_NAI_H5x_D2x_b+ABX*I_NAI_G4x_D2x_b;
  Double I_NAI_G3xy_F3x_b = I_NAI_H4xy_D2x_b+ABX*I_NAI_G3xy_D2x_b;
  Double I_NAI_G3xz_F3x_b = I_NAI_H4xz_D2x_b+ABX*I_NAI_G3xz_D2x_b;
  Double I_NAI_G2x2y_F3x_b = I_NAI_H3x2y_D2x_b+ABX*I_NAI_G2x2y_D2x_b;
  Double I_NAI_G2xyz_F3x_b = I_NAI_H3xyz_D2x_b+ABX*I_NAI_G2xyz_D2x_b;
  Double I_NAI_G2x2z_F3x_b = I_NAI_H3x2z_D2x_b+ABX*I_NAI_G2x2z_D2x_b;
  Double I_NAI_Gx3y_F3x_b = I_NAI_H2x3y_D2x_b+ABX*I_NAI_Gx3y_D2x_b;
  Double I_NAI_Gx2yz_F3x_b = I_NAI_H2x2yz_D2x_b+ABX*I_NAI_Gx2yz_D2x_b;
  Double I_NAI_Gxy2z_F3x_b = I_NAI_H2xy2z_D2x_b+ABX*I_NAI_Gxy2z_D2x_b;
  Double I_NAI_Gx3z_F3x_b = I_NAI_H2x3z_D2x_b+ABX*I_NAI_Gx3z_D2x_b;
  Double I_NAI_G4y_F3x_b = I_NAI_Hx4y_D2x_b+ABX*I_NAI_G4y_D2x_b;
  Double I_NAI_G3yz_F3x_b = I_NAI_Hx3yz_D2x_b+ABX*I_NAI_G3yz_D2x_b;
  Double I_NAI_G2y2z_F3x_b = I_NAI_Hx2y2z_D2x_b+ABX*I_NAI_G2y2z_D2x_b;
  Double I_NAI_Gy3z_F3x_b = I_NAI_Hxy3z_D2x_b+ABX*I_NAI_Gy3z_D2x_b;
  Double I_NAI_G4z_F3x_b = I_NAI_Hx4z_D2x_b+ABX*I_NAI_G4z_D2x_b;
  Double I_NAI_G3xy_F2xy_b = I_NAI_H3x2y_D2x_b+ABY*I_NAI_G3xy_D2x_b;
  Double I_NAI_G3xz_F2xy_b = I_NAI_H3xyz_D2x_b+ABY*I_NAI_G3xz_D2x_b;
  Double I_NAI_G2x2y_F2xy_b = I_NAI_H2x3y_D2x_b+ABY*I_NAI_G2x2y_D2x_b;
  Double I_NAI_G2xyz_F2xy_b = I_NAI_H2x2yz_D2x_b+ABY*I_NAI_G2xyz_D2x_b;
  Double I_NAI_G2x2z_F2xy_b = I_NAI_H2xy2z_D2x_b+ABY*I_NAI_G2x2z_D2x_b;
  Double I_NAI_Gx3y_F2xy_b = I_NAI_Hx4y_D2x_b+ABY*I_NAI_Gx3y_D2x_b;
  Double I_NAI_Gx2yz_F2xy_b = I_NAI_Hx3yz_D2x_b+ABY*I_NAI_Gx2yz_D2x_b;
  Double I_NAI_Gxy2z_F2xy_b = I_NAI_Hx2y2z_D2x_b+ABY*I_NAI_Gxy2z_D2x_b;
  Double I_NAI_Gx3z_F2xy_b = I_NAI_Hxy3z_D2x_b+ABY*I_NAI_Gx3z_D2x_b;
  Double I_NAI_G4y_F2xy_b = I_NAI_H5y_D2x_b+ABY*I_NAI_G4y_D2x_b;
  Double I_NAI_G3yz_F2xy_b = I_NAI_H4yz_D2x_b+ABY*I_NAI_G3yz_D2x_b;
  Double I_NAI_G2y2z_F2xy_b = I_NAI_H3y2z_D2x_b+ABY*I_NAI_G2y2z_D2x_b;
  Double I_NAI_Gy3z_F2xy_b = I_NAI_H2y3z_D2x_b+ABY*I_NAI_Gy3z_D2x_b;
  Double I_NAI_G4z_F2xy_b = I_NAI_Hy4z_D2x_b+ABY*I_NAI_G4z_D2x_b;
  Double I_NAI_G3xz_F2xz_b = I_NAI_H3x2z_D2x_b+ABZ*I_NAI_G3xz_D2x_b;
  Double I_NAI_G2xyz_F2xz_b = I_NAI_H2xy2z_D2x_b+ABZ*I_NAI_G2xyz_D2x_b;
  Double I_NAI_G2x2z_F2xz_b = I_NAI_H2x3z_D2x_b+ABZ*I_NAI_G2x2z_D2x_b;
  Double I_NAI_Gx2yz_F2xz_b = I_NAI_Hx2y2z_D2x_b+ABZ*I_NAI_Gx2yz_D2x_b;
  Double I_NAI_Gxy2z_F2xz_b = I_NAI_Hxy3z_D2x_b+ABZ*I_NAI_Gxy2z_D2x_b;
  Double I_NAI_Gx3z_F2xz_b = I_NAI_Hx4z_D2x_b+ABZ*I_NAI_Gx3z_D2x_b;
  Double I_NAI_G3yz_F2xz_b = I_NAI_H3y2z_D2x_b+ABZ*I_NAI_G3yz_D2x_b;
  Double I_NAI_G2y2z_F2xz_b = I_NAI_H2y3z_D2x_b+ABZ*I_NAI_G2y2z_D2x_b;
  Double I_NAI_Gy3z_F2xz_b = I_NAI_Hy4z_D2x_b+ABZ*I_NAI_Gy3z_D2x_b;
  Double I_NAI_G4z_F2xz_b = I_NAI_H5z_D2x_b+ABZ*I_NAI_G4z_D2x_b;
  Double I_NAI_G3xz_Fx2y_b = I_NAI_H4xz_D2y_b+ABX*I_NAI_G3xz_D2y_b;
  Double I_NAI_G2xyz_Fx2y_b = I_NAI_H3xyz_D2y_b+ABX*I_NAI_G2xyz_D2y_b;
  Double I_NAI_G2x2z_Fx2y_b = I_NAI_H3x2z_D2y_b+ABX*I_NAI_G2x2z_D2y_b;
  Double I_NAI_Gx2yz_Fx2y_b = I_NAI_H2x2yz_D2y_b+ABX*I_NAI_Gx2yz_D2y_b;
  Double I_NAI_Gxy2z_Fx2y_b = I_NAI_H2xy2z_D2y_b+ABX*I_NAI_Gxy2z_D2y_b;
  Double I_NAI_Gx3z_Fx2y_b = I_NAI_H2x3z_D2y_b+ABX*I_NAI_Gx3z_D2y_b;
  Double I_NAI_G3yz_Fx2y_b = I_NAI_Hx3yz_D2y_b+ABX*I_NAI_G3yz_D2y_b;
  Double I_NAI_G2y2z_Fx2y_b = I_NAI_Hx2y2z_D2y_b+ABX*I_NAI_G2y2z_D2y_b;
  Double I_NAI_Gy3z_Fx2y_b = I_NAI_Hxy3z_D2y_b+ABX*I_NAI_Gy3z_D2y_b;
  Double I_NAI_G4z_Fx2y_b = I_NAI_Hx4z_D2y_b+ABX*I_NAI_G4z_D2y_b;
  Double I_NAI_G3xy_Fx2z_b = I_NAI_H4xy_D2z_b+ABX*I_NAI_G3xy_D2z_b;
  Double I_NAI_G2x2y_Fx2z_b = I_NAI_H3x2y_D2z_b+ABX*I_NAI_G2x2y_D2z_b;
  Double I_NAI_G2xyz_Fx2z_b = I_NAI_H3xyz_D2z_b+ABX*I_NAI_G2xyz_D2z_b;
  Double I_NAI_Gx3y_Fx2z_b = I_NAI_H2x3y_D2z_b+ABX*I_NAI_Gx3y_D2z_b;
  Double I_NAI_Gx2yz_Fx2z_b = I_NAI_H2x2yz_D2z_b+ABX*I_NAI_Gx2yz_D2z_b;
  Double I_NAI_Gxy2z_Fx2z_b = I_NAI_H2xy2z_D2z_b+ABX*I_NAI_Gxy2z_D2z_b;
  Double I_NAI_G4y_Fx2z_b = I_NAI_Hx4y_D2z_b+ABX*I_NAI_G4y_D2z_b;
  Double I_NAI_G3yz_Fx2z_b = I_NAI_Hx3yz_D2z_b+ABX*I_NAI_G3yz_D2z_b;
  Double I_NAI_G2y2z_Fx2z_b = I_NAI_Hx2y2z_D2z_b+ABX*I_NAI_G2y2z_D2z_b;
  Double I_NAI_Gy3z_Fx2z_b = I_NAI_Hxy3z_D2z_b+ABX*I_NAI_Gy3z_D2z_b;
  Double I_NAI_G4x_F3y_b = I_NAI_H4xy_D2y_b+ABY*I_NAI_G4x_D2y_b;
  Double I_NAI_G3xy_F3y_b = I_NAI_H3x2y_D2y_b+ABY*I_NAI_G3xy_D2y_b;
  Double I_NAI_G3xz_F3y_b = I_NAI_H3xyz_D2y_b+ABY*I_NAI_G3xz_D2y_b;
  Double I_NAI_G2x2y_F3y_b = I_NAI_H2x3y_D2y_b+ABY*I_NAI_G2x2y_D2y_b;
  Double I_NAI_G2xyz_F3y_b = I_NAI_H2x2yz_D2y_b+ABY*I_NAI_G2xyz_D2y_b;
  Double I_NAI_G2x2z_F3y_b = I_NAI_H2xy2z_D2y_b+ABY*I_NAI_G2x2z_D2y_b;
  Double I_NAI_Gx3y_F3y_b = I_NAI_Hx4y_D2y_b+ABY*I_NAI_Gx3y_D2y_b;
  Double I_NAI_Gx2yz_F3y_b = I_NAI_Hx3yz_D2y_b+ABY*I_NAI_Gx2yz_D2y_b;
  Double I_NAI_Gxy2z_F3y_b = I_NAI_Hx2y2z_D2y_b+ABY*I_NAI_Gxy2z_D2y_b;
  Double I_NAI_Gx3z_F3y_b = I_NAI_Hxy3z_D2y_b+ABY*I_NAI_Gx3z_D2y_b;
  Double I_NAI_G4y_F3y_b = I_NAI_H5y_D2y_b+ABY*I_NAI_G4y_D2y_b;
  Double I_NAI_G3yz_F3y_b = I_NAI_H4yz_D2y_b+ABY*I_NAI_G3yz_D2y_b;
  Double I_NAI_G2y2z_F3y_b = I_NAI_H3y2z_D2y_b+ABY*I_NAI_G2y2z_D2y_b;
  Double I_NAI_Gy3z_F3y_b = I_NAI_H2y3z_D2y_b+ABY*I_NAI_Gy3z_D2y_b;
  Double I_NAI_G4z_F3y_b = I_NAI_Hy4z_D2y_b+ABY*I_NAI_G4z_D2y_b;
  Double I_NAI_G3xz_F2yz_b = I_NAI_H3x2z_D2y_b+ABZ*I_NAI_G3xz_D2y_b;
  Double I_NAI_G2xyz_F2yz_b = I_NAI_H2xy2z_D2y_b+ABZ*I_NAI_G2xyz_D2y_b;
  Double I_NAI_G2x2z_F2yz_b = I_NAI_H2x3z_D2y_b+ABZ*I_NAI_G2x2z_D2y_b;
  Double I_NAI_Gx2yz_F2yz_b = I_NAI_Hx2y2z_D2y_b+ABZ*I_NAI_Gx2yz_D2y_b;
  Double I_NAI_Gxy2z_F2yz_b = I_NAI_Hxy3z_D2y_b+ABZ*I_NAI_Gxy2z_D2y_b;
  Double I_NAI_Gx3z_F2yz_b = I_NAI_Hx4z_D2y_b+ABZ*I_NAI_Gx3z_D2y_b;
  Double I_NAI_G3yz_F2yz_b = I_NAI_H3y2z_D2y_b+ABZ*I_NAI_G3yz_D2y_b;
  Double I_NAI_G2y2z_F2yz_b = I_NAI_H2y3z_D2y_b+ABZ*I_NAI_G2y2z_D2y_b;
  Double I_NAI_Gy3z_F2yz_b = I_NAI_Hy4z_D2y_b+ABZ*I_NAI_Gy3z_D2y_b;
  Double I_NAI_G4z_F2yz_b = I_NAI_H5z_D2y_b+ABZ*I_NAI_G4z_D2y_b;
  Double I_NAI_G4x_F3z_b = I_NAI_H4xz_D2z_b+ABZ*I_NAI_G4x_D2z_b;
  Double I_NAI_G3xy_F3z_b = I_NAI_H3xyz_D2z_b+ABZ*I_NAI_G3xy_D2z_b;
  Double I_NAI_G3xz_F3z_b = I_NAI_H3x2z_D2z_b+ABZ*I_NAI_G3xz_D2z_b;
  Double I_NAI_G2x2y_F3z_b = I_NAI_H2x2yz_D2z_b+ABZ*I_NAI_G2x2y_D2z_b;
  Double I_NAI_G2xyz_F3z_b = I_NAI_H2xy2z_D2z_b+ABZ*I_NAI_G2xyz_D2z_b;
  Double I_NAI_G2x2z_F3z_b = I_NAI_H2x3z_D2z_b+ABZ*I_NAI_G2x2z_D2z_b;
  Double I_NAI_Gx3y_F3z_b = I_NAI_Hx3yz_D2z_b+ABZ*I_NAI_Gx3y_D2z_b;
  Double I_NAI_Gx2yz_F3z_b = I_NAI_Hx2y2z_D2z_b+ABZ*I_NAI_Gx2yz_D2z_b;
  Double I_NAI_Gxy2z_F3z_b = I_NAI_Hxy3z_D2z_b+ABZ*I_NAI_Gxy2z_D2z_b;
  Double I_NAI_Gx3z_F3z_b = I_NAI_Hx4z_D2z_b+ABZ*I_NAI_Gx3z_D2z_b;
  Double I_NAI_G4y_F3z_b = I_NAI_H4yz_D2z_b+ABZ*I_NAI_G4y_D2z_b;
  Double I_NAI_G3yz_F3z_b = I_NAI_H3y2z_D2z_b+ABZ*I_NAI_G3yz_D2z_b;
  Double I_NAI_G2y2z_F3z_b = I_NAI_H2y3z_D2z_b+ABZ*I_NAI_G2y2z_D2z_b;
  Double I_NAI_Gy3z_F3z_b = I_NAI_Hy4z_D2z_b+ABZ*I_NAI_Gy3z_D2z_b;
  Double I_NAI_G4z_F3z_b = I_NAI_H5z_D2z_b+ABZ*I_NAI_G4z_D2z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_G_b
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_F_b
   * RHS shell quartet name: SQ_NAI_F_F_b
   ************************************************************/
  Double I_NAI_F3x_G4x_b = I_NAI_G4x_F3x_b+ABX*I_NAI_F3x_F3x_b;
  Double I_NAI_F2xy_G4x_b = I_NAI_G3xy_F3x_b+ABX*I_NAI_F2xy_F3x_b;
  Double I_NAI_F2xz_G4x_b = I_NAI_G3xz_F3x_b+ABX*I_NAI_F2xz_F3x_b;
  Double I_NAI_Fx2y_G4x_b = I_NAI_G2x2y_F3x_b+ABX*I_NAI_Fx2y_F3x_b;
  Double I_NAI_Fxyz_G4x_b = I_NAI_G2xyz_F3x_b+ABX*I_NAI_Fxyz_F3x_b;
  Double I_NAI_Fx2z_G4x_b = I_NAI_G2x2z_F3x_b+ABX*I_NAI_Fx2z_F3x_b;
  Double I_NAI_F3y_G4x_b = I_NAI_Gx3y_F3x_b+ABX*I_NAI_F3y_F3x_b;
  Double I_NAI_F2yz_G4x_b = I_NAI_Gx2yz_F3x_b+ABX*I_NAI_F2yz_F3x_b;
  Double I_NAI_Fy2z_G4x_b = I_NAI_Gxy2z_F3x_b+ABX*I_NAI_Fy2z_F3x_b;
  Double I_NAI_F3z_G4x_b = I_NAI_Gx3z_F3x_b+ABX*I_NAI_F3z_F3x_b;
  Double I_NAI_F3x_G3xy_b = I_NAI_G3xy_F3x_b+ABY*I_NAI_F3x_F3x_b;
  Double I_NAI_F2xy_G3xy_b = I_NAI_G2x2y_F3x_b+ABY*I_NAI_F2xy_F3x_b;
  Double I_NAI_F2xz_G3xy_b = I_NAI_G2xyz_F3x_b+ABY*I_NAI_F2xz_F3x_b;
  Double I_NAI_Fx2y_G3xy_b = I_NAI_Gx3y_F3x_b+ABY*I_NAI_Fx2y_F3x_b;
  Double I_NAI_Fxyz_G3xy_b = I_NAI_Gx2yz_F3x_b+ABY*I_NAI_Fxyz_F3x_b;
  Double I_NAI_Fx2z_G3xy_b = I_NAI_Gxy2z_F3x_b+ABY*I_NAI_Fx2z_F3x_b;
  Double I_NAI_F3y_G3xy_b = I_NAI_G4y_F3x_b+ABY*I_NAI_F3y_F3x_b;
  Double I_NAI_F2yz_G3xy_b = I_NAI_G3yz_F3x_b+ABY*I_NAI_F2yz_F3x_b;
  Double I_NAI_Fy2z_G3xy_b = I_NAI_G2y2z_F3x_b+ABY*I_NAI_Fy2z_F3x_b;
  Double I_NAI_F3z_G3xy_b = I_NAI_Gy3z_F3x_b+ABY*I_NAI_F3z_F3x_b;
  Double I_NAI_F3x_G3xz_b = I_NAI_G3xz_F3x_b+ABZ*I_NAI_F3x_F3x_b;
  Double I_NAI_F2xy_G3xz_b = I_NAI_G2xyz_F3x_b+ABZ*I_NAI_F2xy_F3x_b;
  Double I_NAI_F2xz_G3xz_b = I_NAI_G2x2z_F3x_b+ABZ*I_NAI_F2xz_F3x_b;
  Double I_NAI_Fx2y_G3xz_b = I_NAI_Gx2yz_F3x_b+ABZ*I_NAI_Fx2y_F3x_b;
  Double I_NAI_Fxyz_G3xz_b = I_NAI_Gxy2z_F3x_b+ABZ*I_NAI_Fxyz_F3x_b;
  Double I_NAI_Fx2z_G3xz_b = I_NAI_Gx3z_F3x_b+ABZ*I_NAI_Fx2z_F3x_b;
  Double I_NAI_F3y_G3xz_b = I_NAI_G3yz_F3x_b+ABZ*I_NAI_F3y_F3x_b;
  Double I_NAI_F2yz_G3xz_b = I_NAI_G2y2z_F3x_b+ABZ*I_NAI_F2yz_F3x_b;
  Double I_NAI_Fy2z_G3xz_b = I_NAI_Gy3z_F3x_b+ABZ*I_NAI_Fy2z_F3x_b;
  Double I_NAI_F3z_G3xz_b = I_NAI_G4z_F3x_b+ABZ*I_NAI_F3z_F3x_b;
  Double I_NAI_F3x_G2x2y_b = I_NAI_G3xy_F2xy_b+ABY*I_NAI_F3x_F2xy_b;
  Double I_NAI_F2xy_G2x2y_b = I_NAI_G2x2y_F2xy_b+ABY*I_NAI_F2xy_F2xy_b;
  Double I_NAI_F2xz_G2x2y_b = I_NAI_G2xyz_F2xy_b+ABY*I_NAI_F2xz_F2xy_b;
  Double I_NAI_Fx2y_G2x2y_b = I_NAI_Gx3y_F2xy_b+ABY*I_NAI_Fx2y_F2xy_b;
  Double I_NAI_Fxyz_G2x2y_b = I_NAI_Gx2yz_F2xy_b+ABY*I_NAI_Fxyz_F2xy_b;
  Double I_NAI_Fx2z_G2x2y_b = I_NAI_Gxy2z_F2xy_b+ABY*I_NAI_Fx2z_F2xy_b;
  Double I_NAI_F3y_G2x2y_b = I_NAI_G4y_F2xy_b+ABY*I_NAI_F3y_F2xy_b;
  Double I_NAI_F2yz_G2x2y_b = I_NAI_G3yz_F2xy_b+ABY*I_NAI_F2yz_F2xy_b;
  Double I_NAI_Fy2z_G2x2y_b = I_NAI_G2y2z_F2xy_b+ABY*I_NAI_Fy2z_F2xy_b;
  Double I_NAI_F3z_G2x2y_b = I_NAI_Gy3z_F2xy_b+ABY*I_NAI_F3z_F2xy_b;
  Double I_NAI_F3x_G2xyz_b = I_NAI_G3xz_F2xy_b+ABZ*I_NAI_F3x_F2xy_b;
  Double I_NAI_F2xy_G2xyz_b = I_NAI_G2xyz_F2xy_b+ABZ*I_NAI_F2xy_F2xy_b;
  Double I_NAI_F2xz_G2xyz_b = I_NAI_G2x2z_F2xy_b+ABZ*I_NAI_F2xz_F2xy_b;
  Double I_NAI_Fx2y_G2xyz_b = I_NAI_Gx2yz_F2xy_b+ABZ*I_NAI_Fx2y_F2xy_b;
  Double I_NAI_Fxyz_G2xyz_b = I_NAI_Gxy2z_F2xy_b+ABZ*I_NAI_Fxyz_F2xy_b;
  Double I_NAI_Fx2z_G2xyz_b = I_NAI_Gx3z_F2xy_b+ABZ*I_NAI_Fx2z_F2xy_b;
  Double I_NAI_F3y_G2xyz_b = I_NAI_G3yz_F2xy_b+ABZ*I_NAI_F3y_F2xy_b;
  Double I_NAI_F2yz_G2xyz_b = I_NAI_G2y2z_F2xy_b+ABZ*I_NAI_F2yz_F2xy_b;
  Double I_NAI_Fy2z_G2xyz_b = I_NAI_Gy3z_F2xy_b+ABZ*I_NAI_Fy2z_F2xy_b;
  Double I_NAI_F3z_G2xyz_b = I_NAI_G4z_F2xy_b+ABZ*I_NAI_F3z_F2xy_b;
  Double I_NAI_F3x_G2x2z_b = I_NAI_G3xz_F2xz_b+ABZ*I_NAI_F3x_F2xz_b;
  Double I_NAI_F2xy_G2x2z_b = I_NAI_G2xyz_F2xz_b+ABZ*I_NAI_F2xy_F2xz_b;
  Double I_NAI_F2xz_G2x2z_b = I_NAI_G2x2z_F2xz_b+ABZ*I_NAI_F2xz_F2xz_b;
  Double I_NAI_Fx2y_G2x2z_b = I_NAI_Gx2yz_F2xz_b+ABZ*I_NAI_Fx2y_F2xz_b;
  Double I_NAI_Fxyz_G2x2z_b = I_NAI_Gxy2z_F2xz_b+ABZ*I_NAI_Fxyz_F2xz_b;
  Double I_NAI_Fx2z_G2x2z_b = I_NAI_Gx3z_F2xz_b+ABZ*I_NAI_Fx2z_F2xz_b;
  Double I_NAI_F3y_G2x2z_b = I_NAI_G3yz_F2xz_b+ABZ*I_NAI_F3y_F2xz_b;
  Double I_NAI_F2yz_G2x2z_b = I_NAI_G2y2z_F2xz_b+ABZ*I_NAI_F2yz_F2xz_b;
  Double I_NAI_Fy2z_G2x2z_b = I_NAI_Gy3z_F2xz_b+ABZ*I_NAI_Fy2z_F2xz_b;
  Double I_NAI_F3z_G2x2z_b = I_NAI_G4z_F2xz_b+ABZ*I_NAI_F3z_F2xz_b;
  Double I_NAI_F3x_Gx3y_b = I_NAI_G4x_F3y_b+ABX*I_NAI_F3x_F3y_b;
  Double I_NAI_F2xy_Gx3y_b = I_NAI_G3xy_F3y_b+ABX*I_NAI_F2xy_F3y_b;
  Double I_NAI_F2xz_Gx3y_b = I_NAI_G3xz_F3y_b+ABX*I_NAI_F2xz_F3y_b;
  Double I_NAI_Fx2y_Gx3y_b = I_NAI_G2x2y_F3y_b+ABX*I_NAI_Fx2y_F3y_b;
  Double I_NAI_Fxyz_Gx3y_b = I_NAI_G2xyz_F3y_b+ABX*I_NAI_Fxyz_F3y_b;
  Double I_NAI_Fx2z_Gx3y_b = I_NAI_G2x2z_F3y_b+ABX*I_NAI_Fx2z_F3y_b;
  Double I_NAI_F3y_Gx3y_b = I_NAI_Gx3y_F3y_b+ABX*I_NAI_F3y_F3y_b;
  Double I_NAI_F2yz_Gx3y_b = I_NAI_Gx2yz_F3y_b+ABX*I_NAI_F2yz_F3y_b;
  Double I_NAI_Fy2z_Gx3y_b = I_NAI_Gxy2z_F3y_b+ABX*I_NAI_Fy2z_F3y_b;
  Double I_NAI_F3z_Gx3y_b = I_NAI_Gx3z_F3y_b+ABX*I_NAI_F3z_F3y_b;
  Double I_NAI_F3x_Gx2yz_b = I_NAI_G3xz_Fx2y_b+ABZ*I_NAI_F3x_Fx2y_b;
  Double I_NAI_F2xy_Gx2yz_b = I_NAI_G2xyz_Fx2y_b+ABZ*I_NAI_F2xy_Fx2y_b;
  Double I_NAI_F2xz_Gx2yz_b = I_NAI_G2x2z_Fx2y_b+ABZ*I_NAI_F2xz_Fx2y_b;
  Double I_NAI_Fx2y_Gx2yz_b = I_NAI_Gx2yz_Fx2y_b+ABZ*I_NAI_Fx2y_Fx2y_b;
  Double I_NAI_Fxyz_Gx2yz_b = I_NAI_Gxy2z_Fx2y_b+ABZ*I_NAI_Fxyz_Fx2y_b;
  Double I_NAI_Fx2z_Gx2yz_b = I_NAI_Gx3z_Fx2y_b+ABZ*I_NAI_Fx2z_Fx2y_b;
  Double I_NAI_F3y_Gx2yz_b = I_NAI_G3yz_Fx2y_b+ABZ*I_NAI_F3y_Fx2y_b;
  Double I_NAI_F2yz_Gx2yz_b = I_NAI_G2y2z_Fx2y_b+ABZ*I_NAI_F2yz_Fx2y_b;
  Double I_NAI_Fy2z_Gx2yz_b = I_NAI_Gy3z_Fx2y_b+ABZ*I_NAI_Fy2z_Fx2y_b;
  Double I_NAI_F3z_Gx2yz_b = I_NAI_G4z_Fx2y_b+ABZ*I_NAI_F3z_Fx2y_b;
  Double I_NAI_F3x_Gxy2z_b = I_NAI_G3xy_Fx2z_b+ABY*I_NAI_F3x_Fx2z_b;
  Double I_NAI_F2xy_Gxy2z_b = I_NAI_G2x2y_Fx2z_b+ABY*I_NAI_F2xy_Fx2z_b;
  Double I_NAI_F2xz_Gxy2z_b = I_NAI_G2xyz_Fx2z_b+ABY*I_NAI_F2xz_Fx2z_b;
  Double I_NAI_Fx2y_Gxy2z_b = I_NAI_Gx3y_Fx2z_b+ABY*I_NAI_Fx2y_Fx2z_b;
  Double I_NAI_Fxyz_Gxy2z_b = I_NAI_Gx2yz_Fx2z_b+ABY*I_NAI_Fxyz_Fx2z_b;
  Double I_NAI_Fx2z_Gxy2z_b = I_NAI_Gxy2z_Fx2z_b+ABY*I_NAI_Fx2z_Fx2z_b;
  Double I_NAI_F3y_Gxy2z_b = I_NAI_G4y_Fx2z_b+ABY*I_NAI_F3y_Fx2z_b;
  Double I_NAI_F2yz_Gxy2z_b = I_NAI_G3yz_Fx2z_b+ABY*I_NAI_F2yz_Fx2z_b;
  Double I_NAI_Fy2z_Gxy2z_b = I_NAI_G2y2z_Fx2z_b+ABY*I_NAI_Fy2z_Fx2z_b;
  Double I_NAI_F3z_Gxy2z_b = I_NAI_Gy3z_Fx2z_b+ABY*I_NAI_F3z_Fx2z_b;
  Double I_NAI_F3x_Gx3z_b = I_NAI_G4x_F3z_b+ABX*I_NAI_F3x_F3z_b;
  Double I_NAI_F2xy_Gx3z_b = I_NAI_G3xy_F3z_b+ABX*I_NAI_F2xy_F3z_b;
  Double I_NAI_F2xz_Gx3z_b = I_NAI_G3xz_F3z_b+ABX*I_NAI_F2xz_F3z_b;
  Double I_NAI_Fx2y_Gx3z_b = I_NAI_G2x2y_F3z_b+ABX*I_NAI_Fx2y_F3z_b;
  Double I_NAI_Fxyz_Gx3z_b = I_NAI_G2xyz_F3z_b+ABX*I_NAI_Fxyz_F3z_b;
  Double I_NAI_Fx2z_Gx3z_b = I_NAI_G2x2z_F3z_b+ABX*I_NAI_Fx2z_F3z_b;
  Double I_NAI_F3y_Gx3z_b = I_NAI_Gx3y_F3z_b+ABX*I_NAI_F3y_F3z_b;
  Double I_NAI_F2yz_Gx3z_b = I_NAI_Gx2yz_F3z_b+ABX*I_NAI_F2yz_F3z_b;
  Double I_NAI_Fy2z_Gx3z_b = I_NAI_Gxy2z_F3z_b+ABX*I_NAI_Fy2z_F3z_b;
  Double I_NAI_F3z_Gx3z_b = I_NAI_Gx3z_F3z_b+ABX*I_NAI_F3z_F3z_b;
  Double I_NAI_F3x_G4y_b = I_NAI_G3xy_F3y_b+ABY*I_NAI_F3x_F3y_b;
  Double I_NAI_F2xy_G4y_b = I_NAI_G2x2y_F3y_b+ABY*I_NAI_F2xy_F3y_b;
  Double I_NAI_F2xz_G4y_b = I_NAI_G2xyz_F3y_b+ABY*I_NAI_F2xz_F3y_b;
  Double I_NAI_Fx2y_G4y_b = I_NAI_Gx3y_F3y_b+ABY*I_NAI_Fx2y_F3y_b;
  Double I_NAI_Fxyz_G4y_b = I_NAI_Gx2yz_F3y_b+ABY*I_NAI_Fxyz_F3y_b;
  Double I_NAI_Fx2z_G4y_b = I_NAI_Gxy2z_F3y_b+ABY*I_NAI_Fx2z_F3y_b;
  Double I_NAI_F3y_G4y_b = I_NAI_G4y_F3y_b+ABY*I_NAI_F3y_F3y_b;
  Double I_NAI_F2yz_G4y_b = I_NAI_G3yz_F3y_b+ABY*I_NAI_F2yz_F3y_b;
  Double I_NAI_Fy2z_G4y_b = I_NAI_G2y2z_F3y_b+ABY*I_NAI_Fy2z_F3y_b;
  Double I_NAI_F3z_G4y_b = I_NAI_Gy3z_F3y_b+ABY*I_NAI_F3z_F3y_b;
  Double I_NAI_F3x_G3yz_b = I_NAI_G3xz_F3y_b+ABZ*I_NAI_F3x_F3y_b;
  Double I_NAI_F2xy_G3yz_b = I_NAI_G2xyz_F3y_b+ABZ*I_NAI_F2xy_F3y_b;
  Double I_NAI_F2xz_G3yz_b = I_NAI_G2x2z_F3y_b+ABZ*I_NAI_F2xz_F3y_b;
  Double I_NAI_Fx2y_G3yz_b = I_NAI_Gx2yz_F3y_b+ABZ*I_NAI_Fx2y_F3y_b;
  Double I_NAI_Fxyz_G3yz_b = I_NAI_Gxy2z_F3y_b+ABZ*I_NAI_Fxyz_F3y_b;
  Double I_NAI_Fx2z_G3yz_b = I_NAI_Gx3z_F3y_b+ABZ*I_NAI_Fx2z_F3y_b;
  Double I_NAI_F3y_G3yz_b = I_NAI_G3yz_F3y_b+ABZ*I_NAI_F3y_F3y_b;
  Double I_NAI_F2yz_G3yz_b = I_NAI_G2y2z_F3y_b+ABZ*I_NAI_F2yz_F3y_b;
  Double I_NAI_Fy2z_G3yz_b = I_NAI_Gy3z_F3y_b+ABZ*I_NAI_Fy2z_F3y_b;
  Double I_NAI_F3z_G3yz_b = I_NAI_G4z_F3y_b+ABZ*I_NAI_F3z_F3y_b;
  Double I_NAI_F3x_G2y2z_b = I_NAI_G3xz_F2yz_b+ABZ*I_NAI_F3x_F2yz_b;
  Double I_NAI_F2xy_G2y2z_b = I_NAI_G2xyz_F2yz_b+ABZ*I_NAI_F2xy_F2yz_b;
  Double I_NAI_F2xz_G2y2z_b = I_NAI_G2x2z_F2yz_b+ABZ*I_NAI_F2xz_F2yz_b;
  Double I_NAI_Fx2y_G2y2z_b = I_NAI_Gx2yz_F2yz_b+ABZ*I_NAI_Fx2y_F2yz_b;
  Double I_NAI_Fxyz_G2y2z_b = I_NAI_Gxy2z_F2yz_b+ABZ*I_NAI_Fxyz_F2yz_b;
  Double I_NAI_Fx2z_G2y2z_b = I_NAI_Gx3z_F2yz_b+ABZ*I_NAI_Fx2z_F2yz_b;
  Double I_NAI_F3y_G2y2z_b = I_NAI_G3yz_F2yz_b+ABZ*I_NAI_F3y_F2yz_b;
  Double I_NAI_F2yz_G2y2z_b = I_NAI_G2y2z_F2yz_b+ABZ*I_NAI_F2yz_F2yz_b;
  Double I_NAI_Fy2z_G2y2z_b = I_NAI_Gy3z_F2yz_b+ABZ*I_NAI_Fy2z_F2yz_b;
  Double I_NAI_F3z_G2y2z_b = I_NAI_G4z_F2yz_b+ABZ*I_NAI_F3z_F2yz_b;
  Double I_NAI_F3x_Gy3z_b = I_NAI_G3xy_F3z_b+ABY*I_NAI_F3x_F3z_b;
  Double I_NAI_F2xy_Gy3z_b = I_NAI_G2x2y_F3z_b+ABY*I_NAI_F2xy_F3z_b;
  Double I_NAI_F2xz_Gy3z_b = I_NAI_G2xyz_F3z_b+ABY*I_NAI_F2xz_F3z_b;
  Double I_NAI_Fx2y_Gy3z_b = I_NAI_Gx3y_F3z_b+ABY*I_NAI_Fx2y_F3z_b;
  Double I_NAI_Fxyz_Gy3z_b = I_NAI_Gx2yz_F3z_b+ABY*I_NAI_Fxyz_F3z_b;
  Double I_NAI_Fx2z_Gy3z_b = I_NAI_Gxy2z_F3z_b+ABY*I_NAI_Fx2z_F3z_b;
  Double I_NAI_F3y_Gy3z_b = I_NAI_G4y_F3z_b+ABY*I_NAI_F3y_F3z_b;
  Double I_NAI_F2yz_Gy3z_b = I_NAI_G3yz_F3z_b+ABY*I_NAI_F2yz_F3z_b;
  Double I_NAI_Fy2z_Gy3z_b = I_NAI_G2y2z_F3z_b+ABY*I_NAI_Fy2z_F3z_b;
  Double I_NAI_F3z_Gy3z_b = I_NAI_Gy3z_F3z_b+ABY*I_NAI_F3z_F3z_b;
  Double I_NAI_F3x_G4z_b = I_NAI_G3xz_F3z_b+ABZ*I_NAI_F3x_F3z_b;
  Double I_NAI_F2xy_G4z_b = I_NAI_G2xyz_F3z_b+ABZ*I_NAI_F2xy_F3z_b;
  Double I_NAI_F2xz_G4z_b = I_NAI_G2x2z_F3z_b+ABZ*I_NAI_F2xz_F3z_b;
  Double I_NAI_Fx2y_G4z_b = I_NAI_Gx2yz_F3z_b+ABZ*I_NAI_Fx2y_F3z_b;
  Double I_NAI_Fxyz_G4z_b = I_NAI_Gxy2z_F3z_b+ABZ*I_NAI_Fxyz_F3z_b;
  Double I_NAI_Fx2z_G4z_b = I_NAI_Gx3z_F3z_b+ABZ*I_NAI_Fx2z_F3z_b;
  Double I_NAI_F3y_G4z_b = I_NAI_G3yz_F3z_b+ABZ*I_NAI_F3y_F3z_b;
  Double I_NAI_F2yz_G4z_b = I_NAI_G2y2z_F3z_b+ABZ*I_NAI_F2yz_F3z_b;
  Double I_NAI_Fy2z_G4z_b = I_NAI_Gy3z_F3z_b+ABZ*I_NAI_Fy2z_F3z_b;
  Double I_NAI_F3z_G4z_b = I_NAI_G4z_F3z_b+ABZ*I_NAI_F3z_F3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_D_F
   ************************************************************/
  abcd[0] = 2.0E0*I_NAI_G4x_F3x_a-3*I_NAI_D2x_F3x;
  abcd[1] = 2.0E0*I_NAI_G3xy_F3x_a-2*I_NAI_Dxy_F3x;
  abcd[2] = 2.0E0*I_NAI_G3xz_F3x_a-2*I_NAI_Dxz_F3x;
  abcd[3] = 2.0E0*I_NAI_G2x2y_F3x_a-1*I_NAI_D2y_F3x;
  abcd[4] = 2.0E0*I_NAI_G2xyz_F3x_a-1*I_NAI_Dyz_F3x;
  abcd[5] = 2.0E0*I_NAI_G2x2z_F3x_a-1*I_NAI_D2z_F3x;
  abcd[6] = 2.0E0*I_NAI_Gx3y_F3x_a;
  abcd[7] = 2.0E0*I_NAI_Gx2yz_F3x_a;
  abcd[8] = 2.0E0*I_NAI_Gxy2z_F3x_a;
  abcd[9] = 2.0E0*I_NAI_Gx3z_F3x_a;
  abcd[10] = 2.0E0*I_NAI_G4x_F2xy_a-3*I_NAI_D2x_F2xy;
  abcd[11] = 2.0E0*I_NAI_G3xy_F2xy_a-2*I_NAI_Dxy_F2xy;
  abcd[12] = 2.0E0*I_NAI_G3xz_F2xy_a-2*I_NAI_Dxz_F2xy;
  abcd[13] = 2.0E0*I_NAI_G2x2y_F2xy_a-1*I_NAI_D2y_F2xy;
  abcd[14] = 2.0E0*I_NAI_G2xyz_F2xy_a-1*I_NAI_Dyz_F2xy;
  abcd[15] = 2.0E0*I_NAI_G2x2z_F2xy_a-1*I_NAI_D2z_F2xy;
  abcd[16] = 2.0E0*I_NAI_Gx3y_F2xy_a;
  abcd[17] = 2.0E0*I_NAI_Gx2yz_F2xy_a;
  abcd[18] = 2.0E0*I_NAI_Gxy2z_F2xy_a;
  abcd[19] = 2.0E0*I_NAI_Gx3z_F2xy_a;
  abcd[20] = 2.0E0*I_NAI_G4x_F2xz_a-3*I_NAI_D2x_F2xz;
  abcd[21] = 2.0E0*I_NAI_G3xy_F2xz_a-2*I_NAI_Dxy_F2xz;
  abcd[22] = 2.0E0*I_NAI_G3xz_F2xz_a-2*I_NAI_Dxz_F2xz;
  abcd[23] = 2.0E0*I_NAI_G2x2y_F2xz_a-1*I_NAI_D2y_F2xz;
  abcd[24] = 2.0E0*I_NAI_G2xyz_F2xz_a-1*I_NAI_Dyz_F2xz;
  abcd[25] = 2.0E0*I_NAI_G2x2z_F2xz_a-1*I_NAI_D2z_F2xz;
  abcd[26] = 2.0E0*I_NAI_Gx3y_F2xz_a;
  abcd[27] = 2.0E0*I_NAI_Gx2yz_F2xz_a;
  abcd[28] = 2.0E0*I_NAI_Gxy2z_F2xz_a;
  abcd[29] = 2.0E0*I_NAI_Gx3z_F2xz_a;
  abcd[30] = 2.0E0*I_NAI_G4x_Fx2y_a-3*I_NAI_D2x_Fx2y;
  abcd[31] = 2.0E0*I_NAI_G3xy_Fx2y_a-2*I_NAI_Dxy_Fx2y;
  abcd[32] = 2.0E0*I_NAI_G3xz_Fx2y_a-2*I_NAI_Dxz_Fx2y;
  abcd[33] = 2.0E0*I_NAI_G2x2y_Fx2y_a-1*I_NAI_D2y_Fx2y;
  abcd[34] = 2.0E0*I_NAI_G2xyz_Fx2y_a-1*I_NAI_Dyz_Fx2y;
  abcd[35] = 2.0E0*I_NAI_G2x2z_Fx2y_a-1*I_NAI_D2z_Fx2y;
  abcd[36] = 2.0E0*I_NAI_Gx3y_Fx2y_a;
  abcd[37] = 2.0E0*I_NAI_Gx2yz_Fx2y_a;
  abcd[38] = 2.0E0*I_NAI_Gxy2z_Fx2y_a;
  abcd[39] = 2.0E0*I_NAI_Gx3z_Fx2y_a;
  abcd[40] = 2.0E0*I_NAI_G4x_Fxyz_a-3*I_NAI_D2x_Fxyz;
  abcd[41] = 2.0E0*I_NAI_G3xy_Fxyz_a-2*I_NAI_Dxy_Fxyz;
  abcd[42] = 2.0E0*I_NAI_G3xz_Fxyz_a-2*I_NAI_Dxz_Fxyz;
  abcd[43] = 2.0E0*I_NAI_G2x2y_Fxyz_a-1*I_NAI_D2y_Fxyz;
  abcd[44] = 2.0E0*I_NAI_G2xyz_Fxyz_a-1*I_NAI_Dyz_Fxyz;
  abcd[45] = 2.0E0*I_NAI_G2x2z_Fxyz_a-1*I_NAI_D2z_Fxyz;
  abcd[46] = 2.0E0*I_NAI_Gx3y_Fxyz_a;
  abcd[47] = 2.0E0*I_NAI_Gx2yz_Fxyz_a;
  abcd[48] = 2.0E0*I_NAI_Gxy2z_Fxyz_a;
  abcd[49] = 2.0E0*I_NAI_Gx3z_Fxyz_a;
  abcd[50] = 2.0E0*I_NAI_G4x_Fx2z_a-3*I_NAI_D2x_Fx2z;
  abcd[51] = 2.0E0*I_NAI_G3xy_Fx2z_a-2*I_NAI_Dxy_Fx2z;
  abcd[52] = 2.0E0*I_NAI_G3xz_Fx2z_a-2*I_NAI_Dxz_Fx2z;
  abcd[53] = 2.0E0*I_NAI_G2x2y_Fx2z_a-1*I_NAI_D2y_Fx2z;
  abcd[54] = 2.0E0*I_NAI_G2xyz_Fx2z_a-1*I_NAI_Dyz_Fx2z;
  abcd[55] = 2.0E0*I_NAI_G2x2z_Fx2z_a-1*I_NAI_D2z_Fx2z;
  abcd[56] = 2.0E0*I_NAI_Gx3y_Fx2z_a;
  abcd[57] = 2.0E0*I_NAI_Gx2yz_Fx2z_a;
  abcd[58] = 2.0E0*I_NAI_Gxy2z_Fx2z_a;
  abcd[59] = 2.0E0*I_NAI_Gx3z_Fx2z_a;
  abcd[60] = 2.0E0*I_NAI_G4x_F3y_a-3*I_NAI_D2x_F3y;
  abcd[61] = 2.0E0*I_NAI_G3xy_F3y_a-2*I_NAI_Dxy_F3y;
  abcd[62] = 2.0E0*I_NAI_G3xz_F3y_a-2*I_NAI_Dxz_F3y;
  abcd[63] = 2.0E0*I_NAI_G2x2y_F3y_a-1*I_NAI_D2y_F3y;
  abcd[64] = 2.0E0*I_NAI_G2xyz_F3y_a-1*I_NAI_Dyz_F3y;
  abcd[65] = 2.0E0*I_NAI_G2x2z_F3y_a-1*I_NAI_D2z_F3y;
  abcd[66] = 2.0E0*I_NAI_Gx3y_F3y_a;
  abcd[67] = 2.0E0*I_NAI_Gx2yz_F3y_a;
  abcd[68] = 2.0E0*I_NAI_Gxy2z_F3y_a;
  abcd[69] = 2.0E0*I_NAI_Gx3z_F3y_a;
  abcd[70] = 2.0E0*I_NAI_G4x_F2yz_a-3*I_NAI_D2x_F2yz;
  abcd[71] = 2.0E0*I_NAI_G3xy_F2yz_a-2*I_NAI_Dxy_F2yz;
  abcd[72] = 2.0E0*I_NAI_G3xz_F2yz_a-2*I_NAI_Dxz_F2yz;
  abcd[73] = 2.0E0*I_NAI_G2x2y_F2yz_a-1*I_NAI_D2y_F2yz;
  abcd[74] = 2.0E0*I_NAI_G2xyz_F2yz_a-1*I_NAI_Dyz_F2yz;
  abcd[75] = 2.0E0*I_NAI_G2x2z_F2yz_a-1*I_NAI_D2z_F2yz;
  abcd[76] = 2.0E0*I_NAI_Gx3y_F2yz_a;
  abcd[77] = 2.0E0*I_NAI_Gx2yz_F2yz_a;
  abcd[78] = 2.0E0*I_NAI_Gxy2z_F2yz_a;
  abcd[79] = 2.0E0*I_NAI_Gx3z_F2yz_a;
  abcd[80] = 2.0E0*I_NAI_G4x_Fy2z_a-3*I_NAI_D2x_Fy2z;
  abcd[81] = 2.0E0*I_NAI_G3xy_Fy2z_a-2*I_NAI_Dxy_Fy2z;
  abcd[82] = 2.0E0*I_NAI_G3xz_Fy2z_a-2*I_NAI_Dxz_Fy2z;
  abcd[83] = 2.0E0*I_NAI_G2x2y_Fy2z_a-1*I_NAI_D2y_Fy2z;
  abcd[84] = 2.0E0*I_NAI_G2xyz_Fy2z_a-1*I_NAI_Dyz_Fy2z;
  abcd[85] = 2.0E0*I_NAI_G2x2z_Fy2z_a-1*I_NAI_D2z_Fy2z;
  abcd[86] = 2.0E0*I_NAI_Gx3y_Fy2z_a;
  abcd[87] = 2.0E0*I_NAI_Gx2yz_Fy2z_a;
  abcd[88] = 2.0E0*I_NAI_Gxy2z_Fy2z_a;
  abcd[89] = 2.0E0*I_NAI_Gx3z_Fy2z_a;
  abcd[90] = 2.0E0*I_NAI_G4x_F3z_a-3*I_NAI_D2x_F3z;
  abcd[91] = 2.0E0*I_NAI_G3xy_F3z_a-2*I_NAI_Dxy_F3z;
  abcd[92] = 2.0E0*I_NAI_G3xz_F3z_a-2*I_NAI_Dxz_F3z;
  abcd[93] = 2.0E0*I_NAI_G2x2y_F3z_a-1*I_NAI_D2y_F3z;
  abcd[94] = 2.0E0*I_NAI_G2xyz_F3z_a-1*I_NAI_Dyz_F3z;
  abcd[95] = 2.0E0*I_NAI_G2x2z_F3z_a-1*I_NAI_D2z_F3z;
  abcd[96] = 2.0E0*I_NAI_Gx3y_F3z_a;
  abcd[97] = 2.0E0*I_NAI_Gx2yz_F3z_a;
  abcd[98] = 2.0E0*I_NAI_Gxy2z_F3z_a;
  abcd[99] = 2.0E0*I_NAI_Gx3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_D_F
   ************************************************************/
  abcd[100] = 2.0E0*I_NAI_G3xy_F3x_a;
  abcd[101] = 2.0E0*I_NAI_G2x2y_F3x_a-1*I_NAI_D2x_F3x;
  abcd[102] = 2.0E0*I_NAI_G2xyz_F3x_a;
  abcd[103] = 2.0E0*I_NAI_Gx3y_F3x_a-2*I_NAI_Dxy_F3x;
  abcd[104] = 2.0E0*I_NAI_Gx2yz_F3x_a-1*I_NAI_Dxz_F3x;
  abcd[105] = 2.0E0*I_NAI_Gxy2z_F3x_a;
  abcd[106] = 2.0E0*I_NAI_G4y_F3x_a-3*I_NAI_D2y_F3x;
  abcd[107] = 2.0E0*I_NAI_G3yz_F3x_a-2*I_NAI_Dyz_F3x;
  abcd[108] = 2.0E0*I_NAI_G2y2z_F3x_a-1*I_NAI_D2z_F3x;
  abcd[109] = 2.0E0*I_NAI_Gy3z_F3x_a;
  abcd[110] = 2.0E0*I_NAI_G3xy_F2xy_a;
  abcd[111] = 2.0E0*I_NAI_G2x2y_F2xy_a-1*I_NAI_D2x_F2xy;
  abcd[112] = 2.0E0*I_NAI_G2xyz_F2xy_a;
  abcd[113] = 2.0E0*I_NAI_Gx3y_F2xy_a-2*I_NAI_Dxy_F2xy;
  abcd[114] = 2.0E0*I_NAI_Gx2yz_F2xy_a-1*I_NAI_Dxz_F2xy;
  abcd[115] = 2.0E0*I_NAI_Gxy2z_F2xy_a;
  abcd[116] = 2.0E0*I_NAI_G4y_F2xy_a-3*I_NAI_D2y_F2xy;
  abcd[117] = 2.0E0*I_NAI_G3yz_F2xy_a-2*I_NAI_Dyz_F2xy;
  abcd[118] = 2.0E0*I_NAI_G2y2z_F2xy_a-1*I_NAI_D2z_F2xy;
  abcd[119] = 2.0E0*I_NAI_Gy3z_F2xy_a;
  abcd[120] = 2.0E0*I_NAI_G3xy_F2xz_a;
  abcd[121] = 2.0E0*I_NAI_G2x2y_F2xz_a-1*I_NAI_D2x_F2xz;
  abcd[122] = 2.0E0*I_NAI_G2xyz_F2xz_a;
  abcd[123] = 2.0E0*I_NAI_Gx3y_F2xz_a-2*I_NAI_Dxy_F2xz;
  abcd[124] = 2.0E0*I_NAI_Gx2yz_F2xz_a-1*I_NAI_Dxz_F2xz;
  abcd[125] = 2.0E0*I_NAI_Gxy2z_F2xz_a;
  abcd[126] = 2.0E0*I_NAI_G4y_F2xz_a-3*I_NAI_D2y_F2xz;
  abcd[127] = 2.0E0*I_NAI_G3yz_F2xz_a-2*I_NAI_Dyz_F2xz;
  abcd[128] = 2.0E0*I_NAI_G2y2z_F2xz_a-1*I_NAI_D2z_F2xz;
  abcd[129] = 2.0E0*I_NAI_Gy3z_F2xz_a;
  abcd[130] = 2.0E0*I_NAI_G3xy_Fx2y_a;
  abcd[131] = 2.0E0*I_NAI_G2x2y_Fx2y_a-1*I_NAI_D2x_Fx2y;
  abcd[132] = 2.0E0*I_NAI_G2xyz_Fx2y_a;
  abcd[133] = 2.0E0*I_NAI_Gx3y_Fx2y_a-2*I_NAI_Dxy_Fx2y;
  abcd[134] = 2.0E0*I_NAI_Gx2yz_Fx2y_a-1*I_NAI_Dxz_Fx2y;
  abcd[135] = 2.0E0*I_NAI_Gxy2z_Fx2y_a;
  abcd[136] = 2.0E0*I_NAI_G4y_Fx2y_a-3*I_NAI_D2y_Fx2y;
  abcd[137] = 2.0E0*I_NAI_G3yz_Fx2y_a-2*I_NAI_Dyz_Fx2y;
  abcd[138] = 2.0E0*I_NAI_G2y2z_Fx2y_a-1*I_NAI_D2z_Fx2y;
  abcd[139] = 2.0E0*I_NAI_Gy3z_Fx2y_a;
  abcd[140] = 2.0E0*I_NAI_G3xy_Fxyz_a;
  abcd[141] = 2.0E0*I_NAI_G2x2y_Fxyz_a-1*I_NAI_D2x_Fxyz;
  abcd[142] = 2.0E0*I_NAI_G2xyz_Fxyz_a;
  abcd[143] = 2.0E0*I_NAI_Gx3y_Fxyz_a-2*I_NAI_Dxy_Fxyz;
  abcd[144] = 2.0E0*I_NAI_Gx2yz_Fxyz_a-1*I_NAI_Dxz_Fxyz;
  abcd[145] = 2.0E0*I_NAI_Gxy2z_Fxyz_a;
  abcd[146] = 2.0E0*I_NAI_G4y_Fxyz_a-3*I_NAI_D2y_Fxyz;
  abcd[147] = 2.0E0*I_NAI_G3yz_Fxyz_a-2*I_NAI_Dyz_Fxyz;
  abcd[148] = 2.0E0*I_NAI_G2y2z_Fxyz_a-1*I_NAI_D2z_Fxyz;
  abcd[149] = 2.0E0*I_NAI_Gy3z_Fxyz_a;
  abcd[150] = 2.0E0*I_NAI_G3xy_Fx2z_a;
  abcd[151] = 2.0E0*I_NAI_G2x2y_Fx2z_a-1*I_NAI_D2x_Fx2z;
  abcd[152] = 2.0E0*I_NAI_G2xyz_Fx2z_a;
  abcd[153] = 2.0E0*I_NAI_Gx3y_Fx2z_a-2*I_NAI_Dxy_Fx2z;
  abcd[154] = 2.0E0*I_NAI_Gx2yz_Fx2z_a-1*I_NAI_Dxz_Fx2z;
  abcd[155] = 2.0E0*I_NAI_Gxy2z_Fx2z_a;
  abcd[156] = 2.0E0*I_NAI_G4y_Fx2z_a-3*I_NAI_D2y_Fx2z;
  abcd[157] = 2.0E0*I_NAI_G3yz_Fx2z_a-2*I_NAI_Dyz_Fx2z;
  abcd[158] = 2.0E0*I_NAI_G2y2z_Fx2z_a-1*I_NAI_D2z_Fx2z;
  abcd[159] = 2.0E0*I_NAI_Gy3z_Fx2z_a;
  abcd[160] = 2.0E0*I_NAI_G3xy_F3y_a;
  abcd[161] = 2.0E0*I_NAI_G2x2y_F3y_a-1*I_NAI_D2x_F3y;
  abcd[162] = 2.0E0*I_NAI_G2xyz_F3y_a;
  abcd[163] = 2.0E0*I_NAI_Gx3y_F3y_a-2*I_NAI_Dxy_F3y;
  abcd[164] = 2.0E0*I_NAI_Gx2yz_F3y_a-1*I_NAI_Dxz_F3y;
  abcd[165] = 2.0E0*I_NAI_Gxy2z_F3y_a;
  abcd[166] = 2.0E0*I_NAI_G4y_F3y_a-3*I_NAI_D2y_F3y;
  abcd[167] = 2.0E0*I_NAI_G3yz_F3y_a-2*I_NAI_Dyz_F3y;
  abcd[168] = 2.0E0*I_NAI_G2y2z_F3y_a-1*I_NAI_D2z_F3y;
  abcd[169] = 2.0E0*I_NAI_Gy3z_F3y_a;
  abcd[170] = 2.0E0*I_NAI_G3xy_F2yz_a;
  abcd[171] = 2.0E0*I_NAI_G2x2y_F2yz_a-1*I_NAI_D2x_F2yz;
  abcd[172] = 2.0E0*I_NAI_G2xyz_F2yz_a;
  abcd[173] = 2.0E0*I_NAI_Gx3y_F2yz_a-2*I_NAI_Dxy_F2yz;
  abcd[174] = 2.0E0*I_NAI_Gx2yz_F2yz_a-1*I_NAI_Dxz_F2yz;
  abcd[175] = 2.0E0*I_NAI_Gxy2z_F2yz_a;
  abcd[176] = 2.0E0*I_NAI_G4y_F2yz_a-3*I_NAI_D2y_F2yz;
  abcd[177] = 2.0E0*I_NAI_G3yz_F2yz_a-2*I_NAI_Dyz_F2yz;
  abcd[178] = 2.0E0*I_NAI_G2y2z_F2yz_a-1*I_NAI_D2z_F2yz;
  abcd[179] = 2.0E0*I_NAI_Gy3z_F2yz_a;
  abcd[180] = 2.0E0*I_NAI_G3xy_Fy2z_a;
  abcd[181] = 2.0E0*I_NAI_G2x2y_Fy2z_a-1*I_NAI_D2x_Fy2z;
  abcd[182] = 2.0E0*I_NAI_G2xyz_Fy2z_a;
  abcd[183] = 2.0E0*I_NAI_Gx3y_Fy2z_a-2*I_NAI_Dxy_Fy2z;
  abcd[184] = 2.0E0*I_NAI_Gx2yz_Fy2z_a-1*I_NAI_Dxz_Fy2z;
  abcd[185] = 2.0E0*I_NAI_Gxy2z_Fy2z_a;
  abcd[186] = 2.0E0*I_NAI_G4y_Fy2z_a-3*I_NAI_D2y_Fy2z;
  abcd[187] = 2.0E0*I_NAI_G3yz_Fy2z_a-2*I_NAI_Dyz_Fy2z;
  abcd[188] = 2.0E0*I_NAI_G2y2z_Fy2z_a-1*I_NAI_D2z_Fy2z;
  abcd[189] = 2.0E0*I_NAI_Gy3z_Fy2z_a;
  abcd[190] = 2.0E0*I_NAI_G3xy_F3z_a;
  abcd[191] = 2.0E0*I_NAI_G2x2y_F3z_a-1*I_NAI_D2x_F3z;
  abcd[192] = 2.0E0*I_NAI_G2xyz_F3z_a;
  abcd[193] = 2.0E0*I_NAI_Gx3y_F3z_a-2*I_NAI_Dxy_F3z;
  abcd[194] = 2.0E0*I_NAI_Gx2yz_F3z_a-1*I_NAI_Dxz_F3z;
  abcd[195] = 2.0E0*I_NAI_Gxy2z_F3z_a;
  abcd[196] = 2.0E0*I_NAI_G4y_F3z_a-3*I_NAI_D2y_F3z;
  abcd[197] = 2.0E0*I_NAI_G3yz_F3z_a-2*I_NAI_Dyz_F3z;
  abcd[198] = 2.0E0*I_NAI_G2y2z_F3z_a-1*I_NAI_D2z_F3z;
  abcd[199] = 2.0E0*I_NAI_Gy3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_G_F_a
   * RHS shell quartet name: SQ_NAI_D_F
   ************************************************************/
  abcd[200] = 2.0E0*I_NAI_G3xz_F3x_a;
  abcd[201] = 2.0E0*I_NAI_G2xyz_F3x_a;
  abcd[202] = 2.0E0*I_NAI_G2x2z_F3x_a-1*I_NAI_D2x_F3x;
  abcd[203] = 2.0E0*I_NAI_Gx2yz_F3x_a;
  abcd[204] = 2.0E0*I_NAI_Gxy2z_F3x_a-1*I_NAI_Dxy_F3x;
  abcd[205] = 2.0E0*I_NAI_Gx3z_F3x_a-2*I_NAI_Dxz_F3x;
  abcd[206] = 2.0E0*I_NAI_G3yz_F3x_a;
  abcd[207] = 2.0E0*I_NAI_G2y2z_F3x_a-1*I_NAI_D2y_F3x;
  abcd[208] = 2.0E0*I_NAI_Gy3z_F3x_a-2*I_NAI_Dyz_F3x;
  abcd[209] = 2.0E0*I_NAI_G4z_F3x_a-3*I_NAI_D2z_F3x;
  abcd[210] = 2.0E0*I_NAI_G3xz_F2xy_a;
  abcd[211] = 2.0E0*I_NAI_G2xyz_F2xy_a;
  abcd[212] = 2.0E0*I_NAI_G2x2z_F2xy_a-1*I_NAI_D2x_F2xy;
  abcd[213] = 2.0E0*I_NAI_Gx2yz_F2xy_a;
  abcd[214] = 2.0E0*I_NAI_Gxy2z_F2xy_a-1*I_NAI_Dxy_F2xy;
  abcd[215] = 2.0E0*I_NAI_Gx3z_F2xy_a-2*I_NAI_Dxz_F2xy;
  abcd[216] = 2.0E0*I_NAI_G3yz_F2xy_a;
  abcd[217] = 2.0E0*I_NAI_G2y2z_F2xy_a-1*I_NAI_D2y_F2xy;
  abcd[218] = 2.0E0*I_NAI_Gy3z_F2xy_a-2*I_NAI_Dyz_F2xy;
  abcd[219] = 2.0E0*I_NAI_G4z_F2xy_a-3*I_NAI_D2z_F2xy;
  abcd[220] = 2.0E0*I_NAI_G3xz_F2xz_a;
  abcd[221] = 2.0E0*I_NAI_G2xyz_F2xz_a;
  abcd[222] = 2.0E0*I_NAI_G2x2z_F2xz_a-1*I_NAI_D2x_F2xz;
  abcd[223] = 2.0E0*I_NAI_Gx2yz_F2xz_a;
  abcd[224] = 2.0E0*I_NAI_Gxy2z_F2xz_a-1*I_NAI_Dxy_F2xz;
  abcd[225] = 2.0E0*I_NAI_Gx3z_F2xz_a-2*I_NAI_Dxz_F2xz;
  abcd[226] = 2.0E0*I_NAI_G3yz_F2xz_a;
  abcd[227] = 2.0E0*I_NAI_G2y2z_F2xz_a-1*I_NAI_D2y_F2xz;
  abcd[228] = 2.0E0*I_NAI_Gy3z_F2xz_a-2*I_NAI_Dyz_F2xz;
  abcd[229] = 2.0E0*I_NAI_G4z_F2xz_a-3*I_NAI_D2z_F2xz;
  abcd[230] = 2.0E0*I_NAI_G3xz_Fx2y_a;
  abcd[231] = 2.0E0*I_NAI_G2xyz_Fx2y_a;
  abcd[232] = 2.0E0*I_NAI_G2x2z_Fx2y_a-1*I_NAI_D2x_Fx2y;
  abcd[233] = 2.0E0*I_NAI_Gx2yz_Fx2y_a;
  abcd[234] = 2.0E0*I_NAI_Gxy2z_Fx2y_a-1*I_NAI_Dxy_Fx2y;
  abcd[235] = 2.0E0*I_NAI_Gx3z_Fx2y_a-2*I_NAI_Dxz_Fx2y;
  abcd[236] = 2.0E0*I_NAI_G3yz_Fx2y_a;
  abcd[237] = 2.0E0*I_NAI_G2y2z_Fx2y_a-1*I_NAI_D2y_Fx2y;
  abcd[238] = 2.0E0*I_NAI_Gy3z_Fx2y_a-2*I_NAI_Dyz_Fx2y;
  abcd[239] = 2.0E0*I_NAI_G4z_Fx2y_a-3*I_NAI_D2z_Fx2y;
  abcd[240] = 2.0E0*I_NAI_G3xz_Fxyz_a;
  abcd[241] = 2.0E0*I_NAI_G2xyz_Fxyz_a;
  abcd[242] = 2.0E0*I_NAI_G2x2z_Fxyz_a-1*I_NAI_D2x_Fxyz;
  abcd[243] = 2.0E0*I_NAI_Gx2yz_Fxyz_a;
  abcd[244] = 2.0E0*I_NAI_Gxy2z_Fxyz_a-1*I_NAI_Dxy_Fxyz;
  abcd[245] = 2.0E0*I_NAI_Gx3z_Fxyz_a-2*I_NAI_Dxz_Fxyz;
  abcd[246] = 2.0E0*I_NAI_G3yz_Fxyz_a;
  abcd[247] = 2.0E0*I_NAI_G2y2z_Fxyz_a-1*I_NAI_D2y_Fxyz;
  abcd[248] = 2.0E0*I_NAI_Gy3z_Fxyz_a-2*I_NAI_Dyz_Fxyz;
  abcd[249] = 2.0E0*I_NAI_G4z_Fxyz_a-3*I_NAI_D2z_Fxyz;
  abcd[250] = 2.0E0*I_NAI_G3xz_Fx2z_a;
  abcd[251] = 2.0E0*I_NAI_G2xyz_Fx2z_a;
  abcd[252] = 2.0E0*I_NAI_G2x2z_Fx2z_a-1*I_NAI_D2x_Fx2z;
  abcd[253] = 2.0E0*I_NAI_Gx2yz_Fx2z_a;
  abcd[254] = 2.0E0*I_NAI_Gxy2z_Fx2z_a-1*I_NAI_Dxy_Fx2z;
  abcd[255] = 2.0E0*I_NAI_Gx3z_Fx2z_a-2*I_NAI_Dxz_Fx2z;
  abcd[256] = 2.0E0*I_NAI_G3yz_Fx2z_a;
  abcd[257] = 2.0E0*I_NAI_G2y2z_Fx2z_a-1*I_NAI_D2y_Fx2z;
  abcd[258] = 2.0E0*I_NAI_Gy3z_Fx2z_a-2*I_NAI_Dyz_Fx2z;
  abcd[259] = 2.0E0*I_NAI_G4z_Fx2z_a-3*I_NAI_D2z_Fx2z;
  abcd[260] = 2.0E0*I_NAI_G3xz_F3y_a;
  abcd[261] = 2.0E0*I_NAI_G2xyz_F3y_a;
  abcd[262] = 2.0E0*I_NAI_G2x2z_F3y_a-1*I_NAI_D2x_F3y;
  abcd[263] = 2.0E0*I_NAI_Gx2yz_F3y_a;
  abcd[264] = 2.0E0*I_NAI_Gxy2z_F3y_a-1*I_NAI_Dxy_F3y;
  abcd[265] = 2.0E0*I_NAI_Gx3z_F3y_a-2*I_NAI_Dxz_F3y;
  abcd[266] = 2.0E0*I_NAI_G3yz_F3y_a;
  abcd[267] = 2.0E0*I_NAI_G2y2z_F3y_a-1*I_NAI_D2y_F3y;
  abcd[268] = 2.0E0*I_NAI_Gy3z_F3y_a-2*I_NAI_Dyz_F3y;
  abcd[269] = 2.0E0*I_NAI_G4z_F3y_a-3*I_NAI_D2z_F3y;
  abcd[270] = 2.0E0*I_NAI_G3xz_F2yz_a;
  abcd[271] = 2.0E0*I_NAI_G2xyz_F2yz_a;
  abcd[272] = 2.0E0*I_NAI_G2x2z_F2yz_a-1*I_NAI_D2x_F2yz;
  abcd[273] = 2.0E0*I_NAI_Gx2yz_F2yz_a;
  abcd[274] = 2.0E0*I_NAI_Gxy2z_F2yz_a-1*I_NAI_Dxy_F2yz;
  abcd[275] = 2.0E0*I_NAI_Gx3z_F2yz_a-2*I_NAI_Dxz_F2yz;
  abcd[276] = 2.0E0*I_NAI_G3yz_F2yz_a;
  abcd[277] = 2.0E0*I_NAI_G2y2z_F2yz_a-1*I_NAI_D2y_F2yz;
  abcd[278] = 2.0E0*I_NAI_Gy3z_F2yz_a-2*I_NAI_Dyz_F2yz;
  abcd[279] = 2.0E0*I_NAI_G4z_F2yz_a-3*I_NAI_D2z_F2yz;
  abcd[280] = 2.0E0*I_NAI_G3xz_Fy2z_a;
  abcd[281] = 2.0E0*I_NAI_G2xyz_Fy2z_a;
  abcd[282] = 2.0E0*I_NAI_G2x2z_Fy2z_a-1*I_NAI_D2x_Fy2z;
  abcd[283] = 2.0E0*I_NAI_Gx2yz_Fy2z_a;
  abcd[284] = 2.0E0*I_NAI_Gxy2z_Fy2z_a-1*I_NAI_Dxy_Fy2z;
  abcd[285] = 2.0E0*I_NAI_Gx3z_Fy2z_a-2*I_NAI_Dxz_Fy2z;
  abcd[286] = 2.0E0*I_NAI_G3yz_Fy2z_a;
  abcd[287] = 2.0E0*I_NAI_G2y2z_Fy2z_a-1*I_NAI_D2y_Fy2z;
  abcd[288] = 2.0E0*I_NAI_Gy3z_Fy2z_a-2*I_NAI_Dyz_Fy2z;
  abcd[289] = 2.0E0*I_NAI_G4z_Fy2z_a-3*I_NAI_D2z_Fy2z;
  abcd[290] = 2.0E0*I_NAI_G3xz_F3z_a;
  abcd[291] = 2.0E0*I_NAI_G2xyz_F3z_a;
  abcd[292] = 2.0E0*I_NAI_G2x2z_F3z_a-1*I_NAI_D2x_F3z;
  abcd[293] = 2.0E0*I_NAI_Gx2yz_F3z_a;
  abcd[294] = 2.0E0*I_NAI_Gxy2z_F3z_a-1*I_NAI_Dxy_F3z;
  abcd[295] = 2.0E0*I_NAI_Gx3z_F3z_a-2*I_NAI_Dxz_F3z;
  abcd[296] = 2.0E0*I_NAI_G3yz_F3z_a;
  abcd[297] = 2.0E0*I_NAI_G2y2z_F3z_a-1*I_NAI_D2y_F3z;
  abcd[298] = 2.0E0*I_NAI_Gy3z_F3z_a-2*I_NAI_Dyz_F3z;
  abcd[299] = 2.0E0*I_NAI_G4z_F3z_a-3*I_NAI_D2z_F3z;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[300] = 2.0E0*I_NAI_F3x_G4x_b-3*I_NAI_F3x_D2x;
  abcd[301] = 2.0E0*I_NAI_F2xy_G4x_b-3*I_NAI_F2xy_D2x;
  abcd[302] = 2.0E0*I_NAI_F2xz_G4x_b-3*I_NAI_F2xz_D2x;
  abcd[303] = 2.0E0*I_NAI_Fx2y_G4x_b-3*I_NAI_Fx2y_D2x;
  abcd[304] = 2.0E0*I_NAI_Fxyz_G4x_b-3*I_NAI_Fxyz_D2x;
  abcd[305] = 2.0E0*I_NAI_Fx2z_G4x_b-3*I_NAI_Fx2z_D2x;
  abcd[306] = 2.0E0*I_NAI_F3y_G4x_b-3*I_NAI_F3y_D2x;
  abcd[307] = 2.0E0*I_NAI_F2yz_G4x_b-3*I_NAI_F2yz_D2x;
  abcd[308] = 2.0E0*I_NAI_Fy2z_G4x_b-3*I_NAI_Fy2z_D2x;
  abcd[309] = 2.0E0*I_NAI_F3z_G4x_b-3*I_NAI_F3z_D2x;
  abcd[310] = 2.0E0*I_NAI_F3x_G3xy_b-2*I_NAI_F3x_Dxy;
  abcd[311] = 2.0E0*I_NAI_F2xy_G3xy_b-2*I_NAI_F2xy_Dxy;
  abcd[312] = 2.0E0*I_NAI_F2xz_G3xy_b-2*I_NAI_F2xz_Dxy;
  abcd[313] = 2.0E0*I_NAI_Fx2y_G3xy_b-2*I_NAI_Fx2y_Dxy;
  abcd[314] = 2.0E0*I_NAI_Fxyz_G3xy_b-2*I_NAI_Fxyz_Dxy;
  abcd[315] = 2.0E0*I_NAI_Fx2z_G3xy_b-2*I_NAI_Fx2z_Dxy;
  abcd[316] = 2.0E0*I_NAI_F3y_G3xy_b-2*I_NAI_F3y_Dxy;
  abcd[317] = 2.0E0*I_NAI_F2yz_G3xy_b-2*I_NAI_F2yz_Dxy;
  abcd[318] = 2.0E0*I_NAI_Fy2z_G3xy_b-2*I_NAI_Fy2z_Dxy;
  abcd[319] = 2.0E0*I_NAI_F3z_G3xy_b-2*I_NAI_F3z_Dxy;
  abcd[320] = 2.0E0*I_NAI_F3x_G3xz_b-2*I_NAI_F3x_Dxz;
  abcd[321] = 2.0E0*I_NAI_F2xy_G3xz_b-2*I_NAI_F2xy_Dxz;
  abcd[322] = 2.0E0*I_NAI_F2xz_G3xz_b-2*I_NAI_F2xz_Dxz;
  abcd[323] = 2.0E0*I_NAI_Fx2y_G3xz_b-2*I_NAI_Fx2y_Dxz;
  abcd[324] = 2.0E0*I_NAI_Fxyz_G3xz_b-2*I_NAI_Fxyz_Dxz;
  abcd[325] = 2.0E0*I_NAI_Fx2z_G3xz_b-2*I_NAI_Fx2z_Dxz;
  abcd[326] = 2.0E0*I_NAI_F3y_G3xz_b-2*I_NAI_F3y_Dxz;
  abcd[327] = 2.0E0*I_NAI_F2yz_G3xz_b-2*I_NAI_F2yz_Dxz;
  abcd[328] = 2.0E0*I_NAI_Fy2z_G3xz_b-2*I_NAI_Fy2z_Dxz;
  abcd[329] = 2.0E0*I_NAI_F3z_G3xz_b-2*I_NAI_F3z_Dxz;
  abcd[330] = 2.0E0*I_NAI_F3x_G2x2y_b-1*I_NAI_F3x_D2y;
  abcd[331] = 2.0E0*I_NAI_F2xy_G2x2y_b-1*I_NAI_F2xy_D2y;
  abcd[332] = 2.0E0*I_NAI_F2xz_G2x2y_b-1*I_NAI_F2xz_D2y;
  abcd[333] = 2.0E0*I_NAI_Fx2y_G2x2y_b-1*I_NAI_Fx2y_D2y;
  abcd[334] = 2.0E0*I_NAI_Fxyz_G2x2y_b-1*I_NAI_Fxyz_D2y;
  abcd[335] = 2.0E0*I_NAI_Fx2z_G2x2y_b-1*I_NAI_Fx2z_D2y;
  abcd[336] = 2.0E0*I_NAI_F3y_G2x2y_b-1*I_NAI_F3y_D2y;
  abcd[337] = 2.0E0*I_NAI_F2yz_G2x2y_b-1*I_NAI_F2yz_D2y;
  abcd[338] = 2.0E0*I_NAI_Fy2z_G2x2y_b-1*I_NAI_Fy2z_D2y;
  abcd[339] = 2.0E0*I_NAI_F3z_G2x2y_b-1*I_NAI_F3z_D2y;
  abcd[340] = 2.0E0*I_NAI_F3x_G2xyz_b-1*I_NAI_F3x_Dyz;
  abcd[341] = 2.0E0*I_NAI_F2xy_G2xyz_b-1*I_NAI_F2xy_Dyz;
  abcd[342] = 2.0E0*I_NAI_F2xz_G2xyz_b-1*I_NAI_F2xz_Dyz;
  abcd[343] = 2.0E0*I_NAI_Fx2y_G2xyz_b-1*I_NAI_Fx2y_Dyz;
  abcd[344] = 2.0E0*I_NAI_Fxyz_G2xyz_b-1*I_NAI_Fxyz_Dyz;
  abcd[345] = 2.0E0*I_NAI_Fx2z_G2xyz_b-1*I_NAI_Fx2z_Dyz;
  abcd[346] = 2.0E0*I_NAI_F3y_G2xyz_b-1*I_NAI_F3y_Dyz;
  abcd[347] = 2.0E0*I_NAI_F2yz_G2xyz_b-1*I_NAI_F2yz_Dyz;
  abcd[348] = 2.0E0*I_NAI_Fy2z_G2xyz_b-1*I_NAI_Fy2z_Dyz;
  abcd[349] = 2.0E0*I_NAI_F3z_G2xyz_b-1*I_NAI_F3z_Dyz;
  abcd[350] = 2.0E0*I_NAI_F3x_G2x2z_b-1*I_NAI_F3x_D2z;
  abcd[351] = 2.0E0*I_NAI_F2xy_G2x2z_b-1*I_NAI_F2xy_D2z;
  abcd[352] = 2.0E0*I_NAI_F2xz_G2x2z_b-1*I_NAI_F2xz_D2z;
  abcd[353] = 2.0E0*I_NAI_Fx2y_G2x2z_b-1*I_NAI_Fx2y_D2z;
  abcd[354] = 2.0E0*I_NAI_Fxyz_G2x2z_b-1*I_NAI_Fxyz_D2z;
  abcd[355] = 2.0E0*I_NAI_Fx2z_G2x2z_b-1*I_NAI_Fx2z_D2z;
  abcd[356] = 2.0E0*I_NAI_F3y_G2x2z_b-1*I_NAI_F3y_D2z;
  abcd[357] = 2.0E0*I_NAI_F2yz_G2x2z_b-1*I_NAI_F2yz_D2z;
  abcd[358] = 2.0E0*I_NAI_Fy2z_G2x2z_b-1*I_NAI_Fy2z_D2z;
  abcd[359] = 2.0E0*I_NAI_F3z_G2x2z_b-1*I_NAI_F3z_D2z;
  abcd[360] = 2.0E0*I_NAI_F3x_Gx3y_b;
  abcd[361] = 2.0E0*I_NAI_F2xy_Gx3y_b;
  abcd[362] = 2.0E0*I_NAI_F2xz_Gx3y_b;
  abcd[363] = 2.0E0*I_NAI_Fx2y_Gx3y_b;
  abcd[364] = 2.0E0*I_NAI_Fxyz_Gx3y_b;
  abcd[365] = 2.0E0*I_NAI_Fx2z_Gx3y_b;
  abcd[366] = 2.0E0*I_NAI_F3y_Gx3y_b;
  abcd[367] = 2.0E0*I_NAI_F2yz_Gx3y_b;
  abcd[368] = 2.0E0*I_NAI_Fy2z_Gx3y_b;
  abcd[369] = 2.0E0*I_NAI_F3z_Gx3y_b;
  abcd[370] = 2.0E0*I_NAI_F3x_Gx2yz_b;
  abcd[371] = 2.0E0*I_NAI_F2xy_Gx2yz_b;
  abcd[372] = 2.0E0*I_NAI_F2xz_Gx2yz_b;
  abcd[373] = 2.0E0*I_NAI_Fx2y_Gx2yz_b;
  abcd[374] = 2.0E0*I_NAI_Fxyz_Gx2yz_b;
  abcd[375] = 2.0E0*I_NAI_Fx2z_Gx2yz_b;
  abcd[376] = 2.0E0*I_NAI_F3y_Gx2yz_b;
  abcd[377] = 2.0E0*I_NAI_F2yz_Gx2yz_b;
  abcd[378] = 2.0E0*I_NAI_Fy2z_Gx2yz_b;
  abcd[379] = 2.0E0*I_NAI_F3z_Gx2yz_b;
  abcd[380] = 2.0E0*I_NAI_F3x_Gxy2z_b;
  abcd[381] = 2.0E0*I_NAI_F2xy_Gxy2z_b;
  abcd[382] = 2.0E0*I_NAI_F2xz_Gxy2z_b;
  abcd[383] = 2.0E0*I_NAI_Fx2y_Gxy2z_b;
  abcd[384] = 2.0E0*I_NAI_Fxyz_Gxy2z_b;
  abcd[385] = 2.0E0*I_NAI_Fx2z_Gxy2z_b;
  abcd[386] = 2.0E0*I_NAI_F3y_Gxy2z_b;
  abcd[387] = 2.0E0*I_NAI_F2yz_Gxy2z_b;
  abcd[388] = 2.0E0*I_NAI_Fy2z_Gxy2z_b;
  abcd[389] = 2.0E0*I_NAI_F3z_Gxy2z_b;
  abcd[390] = 2.0E0*I_NAI_F3x_Gx3z_b;
  abcd[391] = 2.0E0*I_NAI_F2xy_Gx3z_b;
  abcd[392] = 2.0E0*I_NAI_F2xz_Gx3z_b;
  abcd[393] = 2.0E0*I_NAI_Fx2y_Gx3z_b;
  abcd[394] = 2.0E0*I_NAI_Fxyz_Gx3z_b;
  abcd[395] = 2.0E0*I_NAI_Fx2z_Gx3z_b;
  abcd[396] = 2.0E0*I_NAI_F3y_Gx3z_b;
  abcd[397] = 2.0E0*I_NAI_F2yz_Gx3z_b;
  abcd[398] = 2.0E0*I_NAI_Fy2z_Gx3z_b;
  abcd[399] = 2.0E0*I_NAI_F3z_Gx3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[400] = 2.0E0*I_NAI_F3x_G3xy_b;
  abcd[401] = 2.0E0*I_NAI_F2xy_G3xy_b;
  abcd[402] = 2.0E0*I_NAI_F2xz_G3xy_b;
  abcd[403] = 2.0E0*I_NAI_Fx2y_G3xy_b;
  abcd[404] = 2.0E0*I_NAI_Fxyz_G3xy_b;
  abcd[405] = 2.0E0*I_NAI_Fx2z_G3xy_b;
  abcd[406] = 2.0E0*I_NAI_F3y_G3xy_b;
  abcd[407] = 2.0E0*I_NAI_F2yz_G3xy_b;
  abcd[408] = 2.0E0*I_NAI_Fy2z_G3xy_b;
  abcd[409] = 2.0E0*I_NAI_F3z_G3xy_b;
  abcd[410] = 2.0E0*I_NAI_F3x_G2x2y_b-1*I_NAI_F3x_D2x;
  abcd[411] = 2.0E0*I_NAI_F2xy_G2x2y_b-1*I_NAI_F2xy_D2x;
  abcd[412] = 2.0E0*I_NAI_F2xz_G2x2y_b-1*I_NAI_F2xz_D2x;
  abcd[413] = 2.0E0*I_NAI_Fx2y_G2x2y_b-1*I_NAI_Fx2y_D2x;
  abcd[414] = 2.0E0*I_NAI_Fxyz_G2x2y_b-1*I_NAI_Fxyz_D2x;
  abcd[415] = 2.0E0*I_NAI_Fx2z_G2x2y_b-1*I_NAI_Fx2z_D2x;
  abcd[416] = 2.0E0*I_NAI_F3y_G2x2y_b-1*I_NAI_F3y_D2x;
  abcd[417] = 2.0E0*I_NAI_F2yz_G2x2y_b-1*I_NAI_F2yz_D2x;
  abcd[418] = 2.0E0*I_NAI_Fy2z_G2x2y_b-1*I_NAI_Fy2z_D2x;
  abcd[419] = 2.0E0*I_NAI_F3z_G2x2y_b-1*I_NAI_F3z_D2x;
  abcd[420] = 2.0E0*I_NAI_F3x_G2xyz_b;
  abcd[421] = 2.0E0*I_NAI_F2xy_G2xyz_b;
  abcd[422] = 2.0E0*I_NAI_F2xz_G2xyz_b;
  abcd[423] = 2.0E0*I_NAI_Fx2y_G2xyz_b;
  abcd[424] = 2.0E0*I_NAI_Fxyz_G2xyz_b;
  abcd[425] = 2.0E0*I_NAI_Fx2z_G2xyz_b;
  abcd[426] = 2.0E0*I_NAI_F3y_G2xyz_b;
  abcd[427] = 2.0E0*I_NAI_F2yz_G2xyz_b;
  abcd[428] = 2.0E0*I_NAI_Fy2z_G2xyz_b;
  abcd[429] = 2.0E0*I_NAI_F3z_G2xyz_b;
  abcd[430] = 2.0E0*I_NAI_F3x_Gx3y_b-2*I_NAI_F3x_Dxy;
  abcd[431] = 2.0E0*I_NAI_F2xy_Gx3y_b-2*I_NAI_F2xy_Dxy;
  abcd[432] = 2.0E0*I_NAI_F2xz_Gx3y_b-2*I_NAI_F2xz_Dxy;
  abcd[433] = 2.0E0*I_NAI_Fx2y_Gx3y_b-2*I_NAI_Fx2y_Dxy;
  abcd[434] = 2.0E0*I_NAI_Fxyz_Gx3y_b-2*I_NAI_Fxyz_Dxy;
  abcd[435] = 2.0E0*I_NAI_Fx2z_Gx3y_b-2*I_NAI_Fx2z_Dxy;
  abcd[436] = 2.0E0*I_NAI_F3y_Gx3y_b-2*I_NAI_F3y_Dxy;
  abcd[437] = 2.0E0*I_NAI_F2yz_Gx3y_b-2*I_NAI_F2yz_Dxy;
  abcd[438] = 2.0E0*I_NAI_Fy2z_Gx3y_b-2*I_NAI_Fy2z_Dxy;
  abcd[439] = 2.0E0*I_NAI_F3z_Gx3y_b-2*I_NAI_F3z_Dxy;
  abcd[440] = 2.0E0*I_NAI_F3x_Gx2yz_b-1*I_NAI_F3x_Dxz;
  abcd[441] = 2.0E0*I_NAI_F2xy_Gx2yz_b-1*I_NAI_F2xy_Dxz;
  abcd[442] = 2.0E0*I_NAI_F2xz_Gx2yz_b-1*I_NAI_F2xz_Dxz;
  abcd[443] = 2.0E0*I_NAI_Fx2y_Gx2yz_b-1*I_NAI_Fx2y_Dxz;
  abcd[444] = 2.0E0*I_NAI_Fxyz_Gx2yz_b-1*I_NAI_Fxyz_Dxz;
  abcd[445] = 2.0E0*I_NAI_Fx2z_Gx2yz_b-1*I_NAI_Fx2z_Dxz;
  abcd[446] = 2.0E0*I_NAI_F3y_Gx2yz_b-1*I_NAI_F3y_Dxz;
  abcd[447] = 2.0E0*I_NAI_F2yz_Gx2yz_b-1*I_NAI_F2yz_Dxz;
  abcd[448] = 2.0E0*I_NAI_Fy2z_Gx2yz_b-1*I_NAI_Fy2z_Dxz;
  abcd[449] = 2.0E0*I_NAI_F3z_Gx2yz_b-1*I_NAI_F3z_Dxz;
  abcd[450] = 2.0E0*I_NAI_F3x_Gxy2z_b;
  abcd[451] = 2.0E0*I_NAI_F2xy_Gxy2z_b;
  abcd[452] = 2.0E0*I_NAI_F2xz_Gxy2z_b;
  abcd[453] = 2.0E0*I_NAI_Fx2y_Gxy2z_b;
  abcd[454] = 2.0E0*I_NAI_Fxyz_Gxy2z_b;
  abcd[455] = 2.0E0*I_NAI_Fx2z_Gxy2z_b;
  abcd[456] = 2.0E0*I_NAI_F3y_Gxy2z_b;
  abcd[457] = 2.0E0*I_NAI_F2yz_Gxy2z_b;
  abcd[458] = 2.0E0*I_NAI_Fy2z_Gxy2z_b;
  abcd[459] = 2.0E0*I_NAI_F3z_Gxy2z_b;
  abcd[460] = 2.0E0*I_NAI_F3x_G4y_b-3*I_NAI_F3x_D2y;
  abcd[461] = 2.0E0*I_NAI_F2xy_G4y_b-3*I_NAI_F2xy_D2y;
  abcd[462] = 2.0E0*I_NAI_F2xz_G4y_b-3*I_NAI_F2xz_D2y;
  abcd[463] = 2.0E0*I_NAI_Fx2y_G4y_b-3*I_NAI_Fx2y_D2y;
  abcd[464] = 2.0E0*I_NAI_Fxyz_G4y_b-3*I_NAI_Fxyz_D2y;
  abcd[465] = 2.0E0*I_NAI_Fx2z_G4y_b-3*I_NAI_Fx2z_D2y;
  abcd[466] = 2.0E0*I_NAI_F3y_G4y_b-3*I_NAI_F3y_D2y;
  abcd[467] = 2.0E0*I_NAI_F2yz_G4y_b-3*I_NAI_F2yz_D2y;
  abcd[468] = 2.0E0*I_NAI_Fy2z_G4y_b-3*I_NAI_Fy2z_D2y;
  abcd[469] = 2.0E0*I_NAI_F3z_G4y_b-3*I_NAI_F3z_D2y;
  abcd[470] = 2.0E0*I_NAI_F3x_G3yz_b-2*I_NAI_F3x_Dyz;
  abcd[471] = 2.0E0*I_NAI_F2xy_G3yz_b-2*I_NAI_F2xy_Dyz;
  abcd[472] = 2.0E0*I_NAI_F2xz_G3yz_b-2*I_NAI_F2xz_Dyz;
  abcd[473] = 2.0E0*I_NAI_Fx2y_G3yz_b-2*I_NAI_Fx2y_Dyz;
  abcd[474] = 2.0E0*I_NAI_Fxyz_G3yz_b-2*I_NAI_Fxyz_Dyz;
  abcd[475] = 2.0E0*I_NAI_Fx2z_G3yz_b-2*I_NAI_Fx2z_Dyz;
  abcd[476] = 2.0E0*I_NAI_F3y_G3yz_b-2*I_NAI_F3y_Dyz;
  abcd[477] = 2.0E0*I_NAI_F2yz_G3yz_b-2*I_NAI_F2yz_Dyz;
  abcd[478] = 2.0E0*I_NAI_Fy2z_G3yz_b-2*I_NAI_Fy2z_Dyz;
  abcd[479] = 2.0E0*I_NAI_F3z_G3yz_b-2*I_NAI_F3z_Dyz;
  abcd[480] = 2.0E0*I_NAI_F3x_G2y2z_b-1*I_NAI_F3x_D2z;
  abcd[481] = 2.0E0*I_NAI_F2xy_G2y2z_b-1*I_NAI_F2xy_D2z;
  abcd[482] = 2.0E0*I_NAI_F2xz_G2y2z_b-1*I_NAI_F2xz_D2z;
  abcd[483] = 2.0E0*I_NAI_Fx2y_G2y2z_b-1*I_NAI_Fx2y_D2z;
  abcd[484] = 2.0E0*I_NAI_Fxyz_G2y2z_b-1*I_NAI_Fxyz_D2z;
  abcd[485] = 2.0E0*I_NAI_Fx2z_G2y2z_b-1*I_NAI_Fx2z_D2z;
  abcd[486] = 2.0E0*I_NAI_F3y_G2y2z_b-1*I_NAI_F3y_D2z;
  abcd[487] = 2.0E0*I_NAI_F2yz_G2y2z_b-1*I_NAI_F2yz_D2z;
  abcd[488] = 2.0E0*I_NAI_Fy2z_G2y2z_b-1*I_NAI_Fy2z_D2z;
  abcd[489] = 2.0E0*I_NAI_F3z_G2y2z_b-1*I_NAI_F3z_D2z;
  abcd[490] = 2.0E0*I_NAI_F3x_Gy3z_b;
  abcd[491] = 2.0E0*I_NAI_F2xy_Gy3z_b;
  abcd[492] = 2.0E0*I_NAI_F2xz_Gy3z_b;
  abcd[493] = 2.0E0*I_NAI_Fx2y_Gy3z_b;
  abcd[494] = 2.0E0*I_NAI_Fxyz_Gy3z_b;
  abcd[495] = 2.0E0*I_NAI_Fx2z_Gy3z_b;
  abcd[496] = 2.0E0*I_NAI_F3y_Gy3z_b;
  abcd[497] = 2.0E0*I_NAI_F2yz_Gy3z_b;
  abcd[498] = 2.0E0*I_NAI_Fy2z_Gy3z_b;
  abcd[499] = 2.0E0*I_NAI_F3z_Gy3z_b;

  /************************************************************
   * shell quartet name: SQ_NAI_F_F_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_NAI_F_G_b
   * RHS shell quartet name: SQ_NAI_F_D
   ************************************************************/
  abcd[500] = 2.0E0*I_NAI_F3x_G3xz_b;
  abcd[501] = 2.0E0*I_NAI_F2xy_G3xz_b;
  abcd[502] = 2.0E0*I_NAI_F2xz_G3xz_b;
  abcd[503] = 2.0E0*I_NAI_Fx2y_G3xz_b;
  abcd[504] = 2.0E0*I_NAI_Fxyz_G3xz_b;
  abcd[505] = 2.0E0*I_NAI_Fx2z_G3xz_b;
  abcd[506] = 2.0E0*I_NAI_F3y_G3xz_b;
  abcd[507] = 2.0E0*I_NAI_F2yz_G3xz_b;
  abcd[508] = 2.0E0*I_NAI_Fy2z_G3xz_b;
  abcd[509] = 2.0E0*I_NAI_F3z_G3xz_b;
  abcd[510] = 2.0E0*I_NAI_F3x_G2xyz_b;
  abcd[511] = 2.0E0*I_NAI_F2xy_G2xyz_b;
  abcd[512] = 2.0E0*I_NAI_F2xz_G2xyz_b;
  abcd[513] = 2.0E0*I_NAI_Fx2y_G2xyz_b;
  abcd[514] = 2.0E0*I_NAI_Fxyz_G2xyz_b;
  abcd[515] = 2.0E0*I_NAI_Fx2z_G2xyz_b;
  abcd[516] = 2.0E0*I_NAI_F3y_G2xyz_b;
  abcd[517] = 2.0E0*I_NAI_F2yz_G2xyz_b;
  abcd[518] = 2.0E0*I_NAI_Fy2z_G2xyz_b;
  abcd[519] = 2.0E0*I_NAI_F3z_G2xyz_b;
  abcd[520] = 2.0E0*I_NAI_F3x_G2x2z_b-1*I_NAI_F3x_D2x;
  abcd[521] = 2.0E0*I_NAI_F2xy_G2x2z_b-1*I_NAI_F2xy_D2x;
  abcd[522] = 2.0E0*I_NAI_F2xz_G2x2z_b-1*I_NAI_F2xz_D2x;
  abcd[523] = 2.0E0*I_NAI_Fx2y_G2x2z_b-1*I_NAI_Fx2y_D2x;
  abcd[524] = 2.0E0*I_NAI_Fxyz_G2x2z_b-1*I_NAI_Fxyz_D2x;
  abcd[525] = 2.0E0*I_NAI_Fx2z_G2x2z_b-1*I_NAI_Fx2z_D2x;
  abcd[526] = 2.0E0*I_NAI_F3y_G2x2z_b-1*I_NAI_F3y_D2x;
  abcd[527] = 2.0E0*I_NAI_F2yz_G2x2z_b-1*I_NAI_F2yz_D2x;
  abcd[528] = 2.0E0*I_NAI_Fy2z_G2x2z_b-1*I_NAI_Fy2z_D2x;
  abcd[529] = 2.0E0*I_NAI_F3z_G2x2z_b-1*I_NAI_F3z_D2x;
  abcd[530] = 2.0E0*I_NAI_F3x_Gx2yz_b;
  abcd[531] = 2.0E0*I_NAI_F2xy_Gx2yz_b;
  abcd[532] = 2.0E0*I_NAI_F2xz_Gx2yz_b;
  abcd[533] = 2.0E0*I_NAI_Fx2y_Gx2yz_b;
  abcd[534] = 2.0E0*I_NAI_Fxyz_Gx2yz_b;
  abcd[535] = 2.0E0*I_NAI_Fx2z_Gx2yz_b;
  abcd[536] = 2.0E0*I_NAI_F3y_Gx2yz_b;
  abcd[537] = 2.0E0*I_NAI_F2yz_Gx2yz_b;
  abcd[538] = 2.0E0*I_NAI_Fy2z_Gx2yz_b;
  abcd[539] = 2.0E0*I_NAI_F3z_Gx2yz_b;
  abcd[540] = 2.0E0*I_NAI_F3x_Gxy2z_b-1*I_NAI_F3x_Dxy;
  abcd[541] = 2.0E0*I_NAI_F2xy_Gxy2z_b-1*I_NAI_F2xy_Dxy;
  abcd[542] = 2.0E0*I_NAI_F2xz_Gxy2z_b-1*I_NAI_F2xz_Dxy;
  abcd[543] = 2.0E0*I_NAI_Fx2y_Gxy2z_b-1*I_NAI_Fx2y_Dxy;
  abcd[544] = 2.0E0*I_NAI_Fxyz_Gxy2z_b-1*I_NAI_Fxyz_Dxy;
  abcd[545] = 2.0E0*I_NAI_Fx2z_Gxy2z_b-1*I_NAI_Fx2z_Dxy;
  abcd[546] = 2.0E0*I_NAI_F3y_Gxy2z_b-1*I_NAI_F3y_Dxy;
  abcd[547] = 2.0E0*I_NAI_F2yz_Gxy2z_b-1*I_NAI_F2yz_Dxy;
  abcd[548] = 2.0E0*I_NAI_Fy2z_Gxy2z_b-1*I_NAI_Fy2z_Dxy;
  abcd[549] = 2.0E0*I_NAI_F3z_Gxy2z_b-1*I_NAI_F3z_Dxy;
  abcd[550] = 2.0E0*I_NAI_F3x_Gx3z_b-2*I_NAI_F3x_Dxz;
  abcd[551] = 2.0E0*I_NAI_F2xy_Gx3z_b-2*I_NAI_F2xy_Dxz;
  abcd[552] = 2.0E0*I_NAI_F2xz_Gx3z_b-2*I_NAI_F2xz_Dxz;
  abcd[553] = 2.0E0*I_NAI_Fx2y_Gx3z_b-2*I_NAI_Fx2y_Dxz;
  abcd[554] = 2.0E0*I_NAI_Fxyz_Gx3z_b-2*I_NAI_Fxyz_Dxz;
  abcd[555] = 2.0E0*I_NAI_Fx2z_Gx3z_b-2*I_NAI_Fx2z_Dxz;
  abcd[556] = 2.0E0*I_NAI_F3y_Gx3z_b-2*I_NAI_F3y_Dxz;
  abcd[557] = 2.0E0*I_NAI_F2yz_Gx3z_b-2*I_NAI_F2yz_Dxz;
  abcd[558] = 2.0E0*I_NAI_Fy2z_Gx3z_b-2*I_NAI_Fy2z_Dxz;
  abcd[559] = 2.0E0*I_NAI_F3z_Gx3z_b-2*I_NAI_F3z_Dxz;
  abcd[560] = 2.0E0*I_NAI_F3x_G3yz_b;
  abcd[561] = 2.0E0*I_NAI_F2xy_G3yz_b;
  abcd[562] = 2.0E0*I_NAI_F2xz_G3yz_b;
  abcd[563] = 2.0E0*I_NAI_Fx2y_G3yz_b;
  abcd[564] = 2.0E0*I_NAI_Fxyz_G3yz_b;
  abcd[565] = 2.0E0*I_NAI_Fx2z_G3yz_b;
  abcd[566] = 2.0E0*I_NAI_F3y_G3yz_b;
  abcd[567] = 2.0E0*I_NAI_F2yz_G3yz_b;
  abcd[568] = 2.0E0*I_NAI_Fy2z_G3yz_b;
  abcd[569] = 2.0E0*I_NAI_F3z_G3yz_b;
  abcd[570] = 2.0E0*I_NAI_F3x_G2y2z_b-1*I_NAI_F3x_D2y;
  abcd[571] = 2.0E0*I_NAI_F2xy_G2y2z_b-1*I_NAI_F2xy_D2y;
  abcd[572] = 2.0E0*I_NAI_F2xz_G2y2z_b-1*I_NAI_F2xz_D2y;
  abcd[573] = 2.0E0*I_NAI_Fx2y_G2y2z_b-1*I_NAI_Fx2y_D2y;
  abcd[574] = 2.0E0*I_NAI_Fxyz_G2y2z_b-1*I_NAI_Fxyz_D2y;
  abcd[575] = 2.0E0*I_NAI_Fx2z_G2y2z_b-1*I_NAI_Fx2z_D2y;
  abcd[576] = 2.0E0*I_NAI_F3y_G2y2z_b-1*I_NAI_F3y_D2y;
  abcd[577] = 2.0E0*I_NAI_F2yz_G2y2z_b-1*I_NAI_F2yz_D2y;
  abcd[578] = 2.0E0*I_NAI_Fy2z_G2y2z_b-1*I_NAI_Fy2z_D2y;
  abcd[579] = 2.0E0*I_NAI_F3z_G2y2z_b-1*I_NAI_F3z_D2y;
  abcd[580] = 2.0E0*I_NAI_F3x_Gy3z_b-2*I_NAI_F3x_Dyz;
  abcd[581] = 2.0E0*I_NAI_F2xy_Gy3z_b-2*I_NAI_F2xy_Dyz;
  abcd[582] = 2.0E0*I_NAI_F2xz_Gy3z_b-2*I_NAI_F2xz_Dyz;
  abcd[583] = 2.0E0*I_NAI_Fx2y_Gy3z_b-2*I_NAI_Fx2y_Dyz;
  abcd[584] = 2.0E0*I_NAI_Fxyz_Gy3z_b-2*I_NAI_Fxyz_Dyz;
  abcd[585] = 2.0E0*I_NAI_Fx2z_Gy3z_b-2*I_NAI_Fx2z_Dyz;
  abcd[586] = 2.0E0*I_NAI_F3y_Gy3z_b-2*I_NAI_F3y_Dyz;
  abcd[587] = 2.0E0*I_NAI_F2yz_Gy3z_b-2*I_NAI_F2yz_Dyz;
  abcd[588] = 2.0E0*I_NAI_Fy2z_Gy3z_b-2*I_NAI_Fy2z_Dyz;
  abcd[589] = 2.0E0*I_NAI_F3z_Gy3z_b-2*I_NAI_F3z_Dyz;
  abcd[590] = 2.0E0*I_NAI_F3x_G4z_b-3*I_NAI_F3x_D2z;
  abcd[591] = 2.0E0*I_NAI_F2xy_G4z_b-3*I_NAI_F2xy_D2z;
  abcd[592] = 2.0E0*I_NAI_F2xz_G4z_b-3*I_NAI_F2xz_D2z;
  abcd[593] = 2.0E0*I_NAI_Fx2y_G4z_b-3*I_NAI_Fx2y_D2z;
  abcd[594] = 2.0E0*I_NAI_Fxyz_G4z_b-3*I_NAI_Fxyz_D2z;
  abcd[595] = 2.0E0*I_NAI_Fx2z_G4z_b-3*I_NAI_Fx2z_D2z;
  abcd[596] = 2.0E0*I_NAI_F3y_G4z_b-3*I_NAI_F3y_D2z;
  abcd[597] = 2.0E0*I_NAI_F2yz_G4z_b-3*I_NAI_F2yz_D2z;
  abcd[598] = 2.0E0*I_NAI_Fy2z_G4z_b-3*I_NAI_Fy2z_D2z;
  abcd[599] = 2.0E0*I_NAI_F3z_G4z_b-3*I_NAI_F3z_D2z;
}
