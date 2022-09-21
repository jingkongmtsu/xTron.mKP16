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
// BRA1 as redundant position, total RHS integrals evaluated as: 6499
// BRA2 as redundant position, total RHS integrals evaluated as: 5580
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: BRA2
//

//
// @@@@ derivative position-direction information
// BRA1  BRA1
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_esp_h_d_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_M9x_S_aa = 0.0E0;
    Double I_ESP_M8xy_S_aa = 0.0E0;
    Double I_ESP_M8xz_S_aa = 0.0E0;
    Double I_ESP_M7x2y_S_aa = 0.0E0;
    Double I_ESP_M7xyz_S_aa = 0.0E0;
    Double I_ESP_M7x2z_S_aa = 0.0E0;
    Double I_ESP_M6x3y_S_aa = 0.0E0;
    Double I_ESP_M6x2yz_S_aa = 0.0E0;
    Double I_ESP_M6xy2z_S_aa = 0.0E0;
    Double I_ESP_M6x3z_S_aa = 0.0E0;
    Double I_ESP_M5x4y_S_aa = 0.0E0;
    Double I_ESP_M5x3yz_S_aa = 0.0E0;
    Double I_ESP_M5x2y2z_S_aa = 0.0E0;
    Double I_ESP_M5xy3z_S_aa = 0.0E0;
    Double I_ESP_M5x4z_S_aa = 0.0E0;
    Double I_ESP_M4x5y_S_aa = 0.0E0;
    Double I_ESP_M4x4yz_S_aa = 0.0E0;
    Double I_ESP_M4x3y2z_S_aa = 0.0E0;
    Double I_ESP_M4x2y3z_S_aa = 0.0E0;
    Double I_ESP_M4xy4z_S_aa = 0.0E0;
    Double I_ESP_M4x5z_S_aa = 0.0E0;
    Double I_ESP_M3x6y_S_aa = 0.0E0;
    Double I_ESP_M3x5yz_S_aa = 0.0E0;
    Double I_ESP_M3x4y2z_S_aa = 0.0E0;
    Double I_ESP_M3x3y3z_S_aa = 0.0E0;
    Double I_ESP_M3x2y4z_S_aa = 0.0E0;
    Double I_ESP_M3xy5z_S_aa = 0.0E0;
    Double I_ESP_M3x6z_S_aa = 0.0E0;
    Double I_ESP_M2x7y_S_aa = 0.0E0;
    Double I_ESP_M2x6yz_S_aa = 0.0E0;
    Double I_ESP_M2x5y2z_S_aa = 0.0E0;
    Double I_ESP_M2x4y3z_S_aa = 0.0E0;
    Double I_ESP_M2x3y4z_S_aa = 0.0E0;
    Double I_ESP_M2x2y5z_S_aa = 0.0E0;
    Double I_ESP_M2xy6z_S_aa = 0.0E0;
    Double I_ESP_M2x7z_S_aa = 0.0E0;
    Double I_ESP_Mx8y_S_aa = 0.0E0;
    Double I_ESP_Mx7yz_S_aa = 0.0E0;
    Double I_ESP_Mx6y2z_S_aa = 0.0E0;
    Double I_ESP_Mx5y3z_S_aa = 0.0E0;
    Double I_ESP_Mx4y4z_S_aa = 0.0E0;
    Double I_ESP_Mx3y5z_S_aa = 0.0E0;
    Double I_ESP_Mx2y6z_S_aa = 0.0E0;
    Double I_ESP_Mxy7z_S_aa = 0.0E0;
    Double I_ESP_Mx8z_S_aa = 0.0E0;
    Double I_ESP_M9y_S_aa = 0.0E0;
    Double I_ESP_M8yz_S_aa = 0.0E0;
    Double I_ESP_M7y2z_S_aa = 0.0E0;
    Double I_ESP_M6y3z_S_aa = 0.0E0;
    Double I_ESP_M5y4z_S_aa = 0.0E0;
    Double I_ESP_M4y5z_S_aa = 0.0E0;
    Double I_ESP_M3y6z_S_aa = 0.0E0;
    Double I_ESP_M2y7z_S_aa = 0.0E0;
    Double I_ESP_My8z_S_aa = 0.0E0;
    Double I_ESP_M9z_S_aa = 0.0E0;
    Double I_ESP_L8x_S_aa = 0.0E0;
    Double I_ESP_L7xy_S_aa = 0.0E0;
    Double I_ESP_L7xz_S_aa = 0.0E0;
    Double I_ESP_L6x2y_S_aa = 0.0E0;
    Double I_ESP_L6xyz_S_aa = 0.0E0;
    Double I_ESP_L6x2z_S_aa = 0.0E0;
    Double I_ESP_L5x3y_S_aa = 0.0E0;
    Double I_ESP_L5x2yz_S_aa = 0.0E0;
    Double I_ESP_L5xy2z_S_aa = 0.0E0;
    Double I_ESP_L5x3z_S_aa = 0.0E0;
    Double I_ESP_L4x4y_S_aa = 0.0E0;
    Double I_ESP_L4x3yz_S_aa = 0.0E0;
    Double I_ESP_L4x2y2z_S_aa = 0.0E0;
    Double I_ESP_L4xy3z_S_aa = 0.0E0;
    Double I_ESP_L4x4z_S_aa = 0.0E0;
    Double I_ESP_L3x5y_S_aa = 0.0E0;
    Double I_ESP_L3x4yz_S_aa = 0.0E0;
    Double I_ESP_L3x3y2z_S_aa = 0.0E0;
    Double I_ESP_L3x2y3z_S_aa = 0.0E0;
    Double I_ESP_L3xy4z_S_aa = 0.0E0;
    Double I_ESP_L3x5z_S_aa = 0.0E0;
    Double I_ESP_L2x6y_S_aa = 0.0E0;
    Double I_ESP_L2x5yz_S_aa = 0.0E0;
    Double I_ESP_L2x4y2z_S_aa = 0.0E0;
    Double I_ESP_L2x3y3z_S_aa = 0.0E0;
    Double I_ESP_L2x2y4z_S_aa = 0.0E0;
    Double I_ESP_L2xy5z_S_aa = 0.0E0;
    Double I_ESP_L2x6z_S_aa = 0.0E0;
    Double I_ESP_Lx7y_S_aa = 0.0E0;
    Double I_ESP_Lx6yz_S_aa = 0.0E0;
    Double I_ESP_Lx5y2z_S_aa = 0.0E0;
    Double I_ESP_Lx4y3z_S_aa = 0.0E0;
    Double I_ESP_Lx3y4z_S_aa = 0.0E0;
    Double I_ESP_Lx2y5z_S_aa = 0.0E0;
    Double I_ESP_Lxy6z_S_aa = 0.0E0;
    Double I_ESP_Lx7z_S_aa = 0.0E0;
    Double I_ESP_L8y_S_aa = 0.0E0;
    Double I_ESP_L7yz_S_aa = 0.0E0;
    Double I_ESP_L6y2z_S_aa = 0.0E0;
    Double I_ESP_L5y3z_S_aa = 0.0E0;
    Double I_ESP_L4y4z_S_aa = 0.0E0;
    Double I_ESP_L3y5z_S_aa = 0.0E0;
    Double I_ESP_L2y6z_S_aa = 0.0E0;
    Double I_ESP_Ly7z_S_aa = 0.0E0;
    Double I_ESP_L8z_S_aa = 0.0E0;
    Double I_ESP_K7x_S_aa = 0.0E0;
    Double I_ESP_K6xy_S_aa = 0.0E0;
    Double I_ESP_K6xz_S_aa = 0.0E0;
    Double I_ESP_K5x2y_S_aa = 0.0E0;
    Double I_ESP_K5xyz_S_aa = 0.0E0;
    Double I_ESP_K5x2z_S_aa = 0.0E0;
    Double I_ESP_K4x3y_S_aa = 0.0E0;
    Double I_ESP_K4x2yz_S_aa = 0.0E0;
    Double I_ESP_K4xy2z_S_aa = 0.0E0;
    Double I_ESP_K4x3z_S_aa = 0.0E0;
    Double I_ESP_K3x4y_S_aa = 0.0E0;
    Double I_ESP_K3x3yz_S_aa = 0.0E0;
    Double I_ESP_K3x2y2z_S_aa = 0.0E0;
    Double I_ESP_K3xy3z_S_aa = 0.0E0;
    Double I_ESP_K3x4z_S_aa = 0.0E0;
    Double I_ESP_K2x5y_S_aa = 0.0E0;
    Double I_ESP_K2x4yz_S_aa = 0.0E0;
    Double I_ESP_K2x3y2z_S_aa = 0.0E0;
    Double I_ESP_K2x2y3z_S_aa = 0.0E0;
    Double I_ESP_K2xy4z_S_aa = 0.0E0;
    Double I_ESP_K2x5z_S_aa = 0.0E0;
    Double I_ESP_Kx6y_S_aa = 0.0E0;
    Double I_ESP_Kx5yz_S_aa = 0.0E0;
    Double I_ESP_Kx4y2z_S_aa = 0.0E0;
    Double I_ESP_Kx3y3z_S_aa = 0.0E0;
    Double I_ESP_Kx2y4z_S_aa = 0.0E0;
    Double I_ESP_Kxy5z_S_aa = 0.0E0;
    Double I_ESP_Kx6z_S_aa = 0.0E0;
    Double I_ESP_K7y_S_aa = 0.0E0;
    Double I_ESP_K6yz_S_aa = 0.0E0;
    Double I_ESP_K5y2z_S_aa = 0.0E0;
    Double I_ESP_K4y3z_S_aa = 0.0E0;
    Double I_ESP_K3y4z_S_aa = 0.0E0;
    Double I_ESP_K2y5z_S_aa = 0.0E0;
    Double I_ESP_Ky6z_S_aa = 0.0E0;
    Double I_ESP_K7z_S_aa = 0.0E0;
    Double I_ESP_K7x_S_a = 0.0E0;
    Double I_ESP_K6xy_S_a = 0.0E0;
    Double I_ESP_K6xz_S_a = 0.0E0;
    Double I_ESP_K5x2y_S_a = 0.0E0;
    Double I_ESP_K5xyz_S_a = 0.0E0;
    Double I_ESP_K5x2z_S_a = 0.0E0;
    Double I_ESP_K4x3y_S_a = 0.0E0;
    Double I_ESP_K4x2yz_S_a = 0.0E0;
    Double I_ESP_K4xy2z_S_a = 0.0E0;
    Double I_ESP_K4x3z_S_a = 0.0E0;
    Double I_ESP_K3x4y_S_a = 0.0E0;
    Double I_ESP_K3x3yz_S_a = 0.0E0;
    Double I_ESP_K3x2y2z_S_a = 0.0E0;
    Double I_ESP_K3xy3z_S_a = 0.0E0;
    Double I_ESP_K3x4z_S_a = 0.0E0;
    Double I_ESP_K2x5y_S_a = 0.0E0;
    Double I_ESP_K2x4yz_S_a = 0.0E0;
    Double I_ESP_K2x3y2z_S_a = 0.0E0;
    Double I_ESP_K2x2y3z_S_a = 0.0E0;
    Double I_ESP_K2xy4z_S_a = 0.0E0;
    Double I_ESP_K2x5z_S_a = 0.0E0;
    Double I_ESP_Kx6y_S_a = 0.0E0;
    Double I_ESP_Kx5yz_S_a = 0.0E0;
    Double I_ESP_Kx4y2z_S_a = 0.0E0;
    Double I_ESP_Kx3y3z_S_a = 0.0E0;
    Double I_ESP_Kx2y4z_S_a = 0.0E0;
    Double I_ESP_Kxy5z_S_a = 0.0E0;
    Double I_ESP_Kx6z_S_a = 0.0E0;
    Double I_ESP_K7y_S_a = 0.0E0;
    Double I_ESP_K6yz_S_a = 0.0E0;
    Double I_ESP_K5y2z_S_a = 0.0E0;
    Double I_ESP_K4y3z_S_a = 0.0E0;
    Double I_ESP_K3y4z_S_a = 0.0E0;
    Double I_ESP_K2y5z_S_a = 0.0E0;
    Double I_ESP_Ky6z_S_a = 0.0E0;
    Double I_ESP_K7z_S_a = 0.0E0;
    Double I_ESP_I6x_S_a = 0.0E0;
    Double I_ESP_I5xy_S_a = 0.0E0;
    Double I_ESP_I5xz_S_a = 0.0E0;
    Double I_ESP_I4x2y_S_a = 0.0E0;
    Double I_ESP_I4xyz_S_a = 0.0E0;
    Double I_ESP_I4x2z_S_a = 0.0E0;
    Double I_ESP_I3x3y_S_a = 0.0E0;
    Double I_ESP_I3x2yz_S_a = 0.0E0;
    Double I_ESP_I3xy2z_S_a = 0.0E0;
    Double I_ESP_I3x3z_S_a = 0.0E0;
    Double I_ESP_I2x4y_S_a = 0.0E0;
    Double I_ESP_I2x3yz_S_a = 0.0E0;
    Double I_ESP_I2x2y2z_S_a = 0.0E0;
    Double I_ESP_I2xy3z_S_a = 0.0E0;
    Double I_ESP_I2x4z_S_a = 0.0E0;
    Double I_ESP_Ix5y_S_a = 0.0E0;
    Double I_ESP_Ix4yz_S_a = 0.0E0;
    Double I_ESP_Ix3y2z_S_a = 0.0E0;
    Double I_ESP_Ix2y3z_S_a = 0.0E0;
    Double I_ESP_Ixy4z_S_a = 0.0E0;
    Double I_ESP_Ix5z_S_a = 0.0E0;
    Double I_ESP_I6y_S_a = 0.0E0;
    Double I_ESP_I5yz_S_a = 0.0E0;
    Double I_ESP_I4y2z_S_a = 0.0E0;
    Double I_ESP_I3y3z_S_a = 0.0E0;
    Double I_ESP_I2y4z_S_a = 0.0E0;
    Double I_ESP_Iy5z_S_a = 0.0E0;
    Double I_ESP_I6z_S_a = 0.0E0;
    Double I_ESP_H5x_S_a = 0.0E0;
    Double I_ESP_H4xy_S_a = 0.0E0;
    Double I_ESP_H4xz_S_a = 0.0E0;
    Double I_ESP_H3x2y_S_a = 0.0E0;
    Double I_ESP_H3xyz_S_a = 0.0E0;
    Double I_ESP_H3x2z_S_a = 0.0E0;
    Double I_ESP_H2x3y_S_a = 0.0E0;
    Double I_ESP_H2x2yz_S_a = 0.0E0;
    Double I_ESP_H2xy2z_S_a = 0.0E0;
    Double I_ESP_H2x3z_S_a = 0.0E0;
    Double I_ESP_Hx4y_S_a = 0.0E0;
    Double I_ESP_Hx3yz_S_a = 0.0E0;
    Double I_ESP_Hx2y2z_S_a = 0.0E0;
    Double I_ESP_Hxy3z_S_a = 0.0E0;
    Double I_ESP_Hx4z_S_a = 0.0E0;
    Double I_ESP_H5y_S_a = 0.0E0;
    Double I_ESP_H4yz_S_a = 0.0E0;
    Double I_ESP_H3y2z_S_a = 0.0E0;
    Double I_ESP_H2y3z_S_a = 0.0E0;
    Double I_ESP_Hy4z_S_a = 0.0E0;
    Double I_ESP_H5z_S_a = 0.0E0;
    Double I_ESP_H5x_S = 0.0E0;
    Double I_ESP_H4xy_S = 0.0E0;
    Double I_ESP_H4xz_S = 0.0E0;
    Double I_ESP_H3x2y_S = 0.0E0;
    Double I_ESP_H3xyz_S = 0.0E0;
    Double I_ESP_H3x2z_S = 0.0E0;
    Double I_ESP_H2x3y_S = 0.0E0;
    Double I_ESP_H2x2yz_S = 0.0E0;
    Double I_ESP_H2xy2z_S = 0.0E0;
    Double I_ESP_H2x3z_S = 0.0E0;
    Double I_ESP_Hx4y_S = 0.0E0;
    Double I_ESP_Hx3yz_S = 0.0E0;
    Double I_ESP_Hx2y2z_S = 0.0E0;
    Double I_ESP_Hxy3z_S = 0.0E0;
    Double I_ESP_Hx4z_S = 0.0E0;
    Double I_ESP_H5y_S = 0.0E0;
    Double I_ESP_H4yz_S = 0.0E0;
    Double I_ESP_H3y2z_S = 0.0E0;
    Double I_ESP_H2y3z_S = 0.0E0;
    Double I_ESP_Hy4z_S = 0.0E0;
    Double I_ESP_H5z_S = 0.0E0;
    Double I_ESP_G4x_S = 0.0E0;
    Double I_ESP_G3xy_S = 0.0E0;
    Double I_ESP_G3xz_S = 0.0E0;
    Double I_ESP_G2x2y_S = 0.0E0;
    Double I_ESP_G2xyz_S = 0.0E0;
    Double I_ESP_G2x2z_S = 0.0E0;
    Double I_ESP_Gx3y_S = 0.0E0;
    Double I_ESP_Gx2yz_S = 0.0E0;
    Double I_ESP_Gxy2z_S = 0.0E0;
    Double I_ESP_Gx3z_S = 0.0E0;
    Double I_ESP_G4y_S = 0.0E0;
    Double I_ESP_G3yz_S = 0.0E0;
    Double I_ESP_G2y2z_S = 0.0E0;
    Double I_ESP_Gy3z_S = 0.0E0;
    Double I_ESP_G4z_S = 0.0E0;
    Double I_ESP_F3x_S = 0.0E0;
    Double I_ESP_F2xy_S = 0.0E0;
    Double I_ESP_F2xz_S = 0.0E0;
    Double I_ESP_Fx2y_S = 0.0E0;
    Double I_ESP_Fxyz_S = 0.0E0;
    Double I_ESP_Fx2z_S = 0.0E0;
    Double I_ESP_F3y_S = 0.0E0;
    Double I_ESP_F2yz_S = 0.0E0;
    Double I_ESP_Fy2z_S = 0.0E0;
    Double I_ESP_F3z_S = 0.0E0;

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
      Double PRX   = PX - R[iGrid*3  ];
      Double PRY   = PY - R[iGrid*3+1];
      Double PRZ   = PZ - R[iGrid*3+2];
      Double PR2   = PRX*PRX+PRY*PRY+PRZ*PRZ;
      Double u     = rho*PR2;
      Double squ   = sqrt(u);
      Double prefactor = ic2*fbra;

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

      Double I_ESP_S_S_vrr  = 0.0E0;
      Double I_ESP_S_S_M1_vrr  = 0.0E0;
      Double I_ESP_S_S_M2_vrr  = 0.0E0;
      Double I_ESP_S_S_M3_vrr  = 0.0E0;
      Double I_ESP_S_S_M4_vrr  = 0.0E0;
      Double I_ESP_S_S_M5_vrr  = 0.0E0;
      Double I_ESP_S_S_M6_vrr  = 0.0E0;
      Double I_ESP_S_S_M7_vrr  = 0.0E0;
      Double I_ESP_S_S_M8_vrr  = 0.0E0;
      Double I_ESP_S_S_M9_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ESP_S_S_vrr_d  = 0.0E0;
      double I_ESP_S_S_M1_vrr_d  = 0.0E0;
      double I_ESP_S_S_M2_vrr_d  = 0.0E0;
      double I_ESP_S_S_M3_vrr_d  = 0.0E0;
      double I_ESP_S_S_M4_vrr_d  = 0.0E0;
      double I_ESP_S_S_M5_vrr_d  = 0.0E0;
      double I_ESP_S_S_M6_vrr_d  = 0.0E0;
      double I_ESP_S_S_M7_vrr_d  = 0.0E0;
      double I_ESP_S_S_M8_vrr_d  = 0.0E0;
      double I_ESP_S_S_M9_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER53;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER51*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER49*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER47*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER45*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER43*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER41*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER39*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER37*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER35*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER33*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER31*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER29*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER27*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER25*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER23*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = 1.0E0+u2*ONEOVER21*I_ESP_S_S_M9_vrr;
        I_ESP_S_S_M9_vrr = ONEOVER19*I_ESP_S_S_M9_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M9_vrr  = f*I_ESP_S_S_M9_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ESP_S_S_M8_vrr  = ONEOVER17*(u2*I_ESP_S_S_M9_vrr+f);
        I_ESP_S_S_M7_vrr  = ONEOVER15*(u2*I_ESP_S_S_M8_vrr+f);
        I_ESP_S_S_M6_vrr  = ONEOVER13*(u2*I_ESP_S_S_M7_vrr+f);
        I_ESP_S_S_M5_vrr  = ONEOVER11*(u2*I_ESP_S_S_M6_vrr+f);
        I_ESP_S_S_M4_vrr  = ONEOVER9*(u2*I_ESP_S_S_M5_vrr+f);
        I_ESP_S_S_M3_vrr  = ONEOVER7*(u2*I_ESP_S_S_M4_vrr+f);
        I_ESP_S_S_M2_vrr  = ONEOVER5*(u2*I_ESP_S_S_M3_vrr+f);
        I_ESP_S_S_M1_vrr  = ONEOVER3*(u2*I_ESP_S_S_M2_vrr+f);
        I_ESP_S_S_vrr  = ONEOVER1*(u2*I_ESP_S_S_M1_vrr+f);

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
          I_ESP_S_S_vrr_d = fac_d*sqrho_d*TWOOVERSQRTPI;
        }else{
          I_ESP_S_S_vrr_d = (fac_d*sqrho_d/squ_d)*erf(squ_d);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        double oneO2u  = 0.5E0/u_d;
        double eu      = exp(-u_d);
        double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;
        I_ESP_S_S_M1_vrr_d = oneO2u*(1.0E0*I_ESP_S_S_vrr_d-f);
        I_ESP_S_S_M2_vrr_d = oneO2u*(3.0E0*I_ESP_S_S_M1_vrr_d-f);
        I_ESP_S_S_M3_vrr_d = oneO2u*(5.0E0*I_ESP_S_S_M2_vrr_d-f);
        I_ESP_S_S_M4_vrr_d = oneO2u*(7.0E0*I_ESP_S_S_M3_vrr_d-f);
        I_ESP_S_S_M5_vrr_d = oneO2u*(9.0E0*I_ESP_S_S_M4_vrr_d-f);
        I_ESP_S_S_M6_vrr_d = oneO2u*(11.0E0*I_ESP_S_S_M5_vrr_d-f);
        I_ESP_S_S_M7_vrr_d = oneO2u*(13.0E0*I_ESP_S_S_M6_vrr_d-f);
        I_ESP_S_S_M8_vrr_d = oneO2u*(15.0E0*I_ESP_S_S_M7_vrr_d-f);
        I_ESP_S_S_M9_vrr_d = oneO2u*(17.0E0*I_ESP_S_S_M8_vrr_d-f);

        // write the double result back to the float var
        I_ESP_S_S_vrr = static_cast<Double>(I_ESP_S_S_vrr_d);
        I_ESP_S_S_M1_vrr = static_cast<Double>(I_ESP_S_S_M1_vrr_d);
        I_ESP_S_S_M2_vrr = static_cast<Double>(I_ESP_S_S_M2_vrr_d);
        I_ESP_S_S_M3_vrr = static_cast<Double>(I_ESP_S_S_M3_vrr_d);
        I_ESP_S_S_M4_vrr = static_cast<Double>(I_ESP_S_S_M4_vrr_d);
        I_ESP_S_S_M5_vrr = static_cast<Double>(I_ESP_S_S_M5_vrr_d);
        I_ESP_S_S_M6_vrr = static_cast<Double>(I_ESP_S_S_M6_vrr_d);
        I_ESP_S_S_M7_vrr = static_cast<Double>(I_ESP_S_S_M7_vrr_d);
        I_ESP_S_S_M8_vrr = static_cast<Double>(I_ESP_S_S_M8_vrr_d);
        I_ESP_S_S_M9_vrr = static_cast<Double>(I_ESP_S_S_M9_vrr_d);

#else

        // use erf function to get (SS|SS)^{0}
        if (fabs(u)<THRESHOLD_MATH) {
          I_ESP_S_S_vrr = prefactor*sqrho*TWOOVERSQRTPI;
        }else{
          I_ESP_S_S_vrr = (prefactor*sqrho/squ)*erf(squ);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        Double oneO2u = 0.5E0/u;
        Double eu     = exp(-u);
        Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M1_vrr = oneO2u*(1.0E0*I_ESP_S_S_vrr-f);
        I_ESP_S_S_M2_vrr = oneO2u*(3.0E0*I_ESP_S_S_M1_vrr-f);
        I_ESP_S_S_M3_vrr = oneO2u*(5.0E0*I_ESP_S_S_M2_vrr-f);
        I_ESP_S_S_M4_vrr = oneO2u*(7.0E0*I_ESP_S_S_M3_vrr-f);
        I_ESP_S_S_M5_vrr = oneO2u*(9.0E0*I_ESP_S_S_M4_vrr-f);
        I_ESP_S_S_M6_vrr = oneO2u*(11.0E0*I_ESP_S_S_M5_vrr-f);
        I_ESP_S_S_M7_vrr = oneO2u*(13.0E0*I_ESP_S_S_M6_vrr-f);
        I_ESP_S_S_M8_vrr = oneO2u*(15.0E0*I_ESP_S_S_M7_vrr-f);
        I_ESP_S_S_M9_vrr = oneO2u*(17.0E0*I_ESP_S_S_M8_vrr-f);

#endif

      }


        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M8
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M8
         * RHS shell quartet name: SQ_ESP_S_S_M9
         ************************************************************/
        Double I_ESP_Px_S_M8_vrr = PAX*I_ESP_S_S_M8_vrr-PRX*I_ESP_S_S_M9_vrr;
        Double I_ESP_Py_S_M8_vrr = PAY*I_ESP_S_S_M8_vrr-PRY*I_ESP_S_S_M9_vrr;
        Double I_ESP_Pz_S_M8_vrr = PAZ*I_ESP_S_S_M8_vrr-PRZ*I_ESP_S_S_M9_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M7
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M7
         * RHS shell quartet name: SQ_ESP_S_S_M8
         ************************************************************/
        Double I_ESP_Px_S_M7_vrr = PAX*I_ESP_S_S_M7_vrr-PRX*I_ESP_S_S_M8_vrr;
        Double I_ESP_Py_S_M7_vrr = PAY*I_ESP_S_S_M7_vrr-PRY*I_ESP_S_S_M8_vrr;
        Double I_ESP_Pz_S_M7_vrr = PAZ*I_ESP_S_S_M7_vrr-PRZ*I_ESP_S_S_M8_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M7
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M7
         * RHS shell quartet name: SQ_ESP_P_S_M8
         * RHS shell quartet name: SQ_ESP_S_S_M7
         * RHS shell quartet name: SQ_ESP_S_S_M8
         ************************************************************/
        Double I_ESP_D2x_S_M7_vrr = PAX*I_ESP_Px_S_M7_vrr-PRX*I_ESP_Px_S_M8_vrr+oned2z*I_ESP_S_S_M7_vrr-oned2z*I_ESP_S_S_M8_vrr;
        Double I_ESP_D2y_S_M7_vrr = PAY*I_ESP_Py_S_M7_vrr-PRY*I_ESP_Py_S_M8_vrr+oned2z*I_ESP_S_S_M7_vrr-oned2z*I_ESP_S_S_M8_vrr;
        Double I_ESP_D2z_S_M7_vrr = PAZ*I_ESP_Pz_S_M7_vrr-PRZ*I_ESP_Pz_S_M8_vrr+oned2z*I_ESP_S_S_M7_vrr-oned2z*I_ESP_S_S_M8_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M6
         * RHS shell quartet name: SQ_ESP_S_S_M7
         ************************************************************/
        Double I_ESP_Px_S_M6_vrr = PAX*I_ESP_S_S_M6_vrr-PRX*I_ESP_S_S_M7_vrr;
        Double I_ESP_Py_S_M6_vrr = PAY*I_ESP_S_S_M6_vrr-PRY*I_ESP_S_S_M7_vrr;
        Double I_ESP_Pz_S_M6_vrr = PAZ*I_ESP_S_S_M6_vrr-PRZ*I_ESP_S_S_M7_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M6
         * RHS shell quartet name: SQ_ESP_P_S_M7
         * RHS shell quartet name: SQ_ESP_S_S_M6
         * RHS shell quartet name: SQ_ESP_S_S_M7
         ************************************************************/
        Double I_ESP_D2x_S_M6_vrr = PAX*I_ESP_Px_S_M6_vrr-PRX*I_ESP_Px_S_M7_vrr+oned2z*I_ESP_S_S_M6_vrr-oned2z*I_ESP_S_S_M7_vrr;
        Double I_ESP_D2y_S_M6_vrr = PAY*I_ESP_Py_S_M6_vrr-PRY*I_ESP_Py_S_M7_vrr+oned2z*I_ESP_S_S_M6_vrr-oned2z*I_ESP_S_S_M7_vrr;
        Double I_ESP_D2z_S_M6_vrr = PAZ*I_ESP_Pz_S_M6_vrr-PRZ*I_ESP_Pz_S_M7_vrr+oned2z*I_ESP_S_S_M6_vrr-oned2z*I_ESP_S_S_M7_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M6
         * RHS shell quartet name: SQ_ESP_D_S_M7
         * RHS shell quartet name: SQ_ESP_P_S_M6
         * RHS shell quartet name: SQ_ESP_P_S_M7
         ************************************************************/
        Double I_ESP_F3x_S_M6_vrr = PAX*I_ESP_D2x_S_M6_vrr-PRX*I_ESP_D2x_S_M7_vrr+2*oned2z*I_ESP_Px_S_M6_vrr-2*oned2z*I_ESP_Px_S_M7_vrr;
        Double I_ESP_F3y_S_M6_vrr = PAY*I_ESP_D2y_S_M6_vrr-PRY*I_ESP_D2y_S_M7_vrr+2*oned2z*I_ESP_Py_S_M6_vrr-2*oned2z*I_ESP_Py_S_M7_vrr;
        Double I_ESP_F3z_S_M6_vrr = PAZ*I_ESP_D2z_S_M6_vrr-PRZ*I_ESP_D2z_S_M7_vrr+2*oned2z*I_ESP_Pz_S_M6_vrr-2*oned2z*I_ESP_Pz_S_M7_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M5
         * RHS shell quartet name: SQ_ESP_S_S_M6
         ************************************************************/
        Double I_ESP_Px_S_M5_vrr = PAX*I_ESP_S_S_M5_vrr-PRX*I_ESP_S_S_M6_vrr;
        Double I_ESP_Py_S_M5_vrr = PAY*I_ESP_S_S_M5_vrr-PRY*I_ESP_S_S_M6_vrr;
        Double I_ESP_Pz_S_M5_vrr = PAZ*I_ESP_S_S_M5_vrr-PRZ*I_ESP_S_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M5
         * RHS shell quartet name: SQ_ESP_P_S_M6
         * RHS shell quartet name: SQ_ESP_S_S_M5
         * RHS shell quartet name: SQ_ESP_S_S_M6
         ************************************************************/
        Double I_ESP_D2x_S_M5_vrr = PAX*I_ESP_Px_S_M5_vrr-PRX*I_ESP_Px_S_M6_vrr+oned2z*I_ESP_S_S_M5_vrr-oned2z*I_ESP_S_S_M6_vrr;
        Double I_ESP_D2y_S_M5_vrr = PAY*I_ESP_Py_S_M5_vrr-PRY*I_ESP_Py_S_M6_vrr+oned2z*I_ESP_S_S_M5_vrr-oned2z*I_ESP_S_S_M6_vrr;
        Double I_ESP_D2z_S_M5_vrr = PAZ*I_ESP_Pz_S_M5_vrr-PRZ*I_ESP_Pz_S_M6_vrr+oned2z*I_ESP_S_S_M5_vrr-oned2z*I_ESP_S_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M5
         * RHS shell quartet name: SQ_ESP_D_S_M6
         * RHS shell quartet name: SQ_ESP_P_S_M5
         * RHS shell quartet name: SQ_ESP_P_S_M6
         ************************************************************/
        Double I_ESP_F3x_S_M5_vrr = PAX*I_ESP_D2x_S_M5_vrr-PRX*I_ESP_D2x_S_M6_vrr+2*oned2z*I_ESP_Px_S_M5_vrr-2*oned2z*I_ESP_Px_S_M6_vrr;
        Double I_ESP_F3y_S_M5_vrr = PAY*I_ESP_D2y_S_M5_vrr-PRY*I_ESP_D2y_S_M6_vrr+2*oned2z*I_ESP_Py_S_M5_vrr-2*oned2z*I_ESP_Py_S_M6_vrr;
        Double I_ESP_F3z_S_M5_vrr = PAZ*I_ESP_D2z_S_M5_vrr-PRZ*I_ESP_D2z_S_M6_vrr+2*oned2z*I_ESP_Pz_S_M5_vrr-2*oned2z*I_ESP_Pz_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 11 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M5
         * RHS shell quartet name: SQ_ESP_F_S_M6
         * RHS shell quartet name: SQ_ESP_D_S_M5
         * RHS shell quartet name: SQ_ESP_D_S_M6
         ************************************************************/
        Double I_ESP_G4x_S_M5_vrr = PAX*I_ESP_F3x_S_M5_vrr-PRX*I_ESP_F3x_S_M6_vrr+3*oned2z*I_ESP_D2x_S_M5_vrr-3*oned2z*I_ESP_D2x_S_M6_vrr;
        Double I_ESP_G3xy_S_M5_vrr = PAY*I_ESP_F3x_S_M5_vrr-PRY*I_ESP_F3x_S_M6_vrr;
        Double I_ESP_G4y_S_M5_vrr = PAY*I_ESP_F3y_S_M5_vrr-PRY*I_ESP_F3y_S_M6_vrr+3*oned2z*I_ESP_D2y_S_M5_vrr-3*oned2z*I_ESP_D2y_S_M6_vrr;
        Double I_ESP_G4z_S_M5_vrr = PAZ*I_ESP_F3z_S_M5_vrr-PRZ*I_ESP_F3z_S_M6_vrr+3*oned2z*I_ESP_D2z_S_M5_vrr-3*oned2z*I_ESP_D2z_S_M6_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M5
         ************************************************************/
        Double I_ESP_Px_S_M4_vrr = PAX*I_ESP_S_S_M4_vrr-PRX*I_ESP_S_S_M5_vrr;
        Double I_ESP_Py_S_M4_vrr = PAY*I_ESP_S_S_M4_vrr-PRY*I_ESP_S_S_M5_vrr;
        Double I_ESP_Pz_S_M4_vrr = PAZ*I_ESP_S_S_M4_vrr-PRZ*I_ESP_S_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_P_S_M5
         * RHS shell quartet name: SQ_ESP_S_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M5
         ************************************************************/
        Double I_ESP_D2x_S_M4_vrr = PAX*I_ESP_Px_S_M4_vrr-PRX*I_ESP_Px_S_M5_vrr+oned2z*I_ESP_S_S_M4_vrr-oned2z*I_ESP_S_S_M5_vrr;
        Double I_ESP_D2y_S_M4_vrr = PAY*I_ESP_Py_S_M4_vrr-PRY*I_ESP_Py_S_M5_vrr+oned2z*I_ESP_S_S_M4_vrr-oned2z*I_ESP_S_S_M5_vrr;
        Double I_ESP_D2z_S_M4_vrr = PAZ*I_ESP_Pz_S_M4_vrr-PRZ*I_ESP_Pz_S_M5_vrr+oned2z*I_ESP_S_S_M4_vrr-oned2z*I_ESP_S_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M4
         * RHS shell quartet name: SQ_ESP_D_S_M5
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_P_S_M5
         ************************************************************/
        Double I_ESP_F3x_S_M4_vrr = PAX*I_ESP_D2x_S_M4_vrr-PRX*I_ESP_D2x_S_M5_vrr+2*oned2z*I_ESP_Px_S_M4_vrr-2*oned2z*I_ESP_Px_S_M5_vrr;
        Double I_ESP_F3y_S_M4_vrr = PAY*I_ESP_D2y_S_M4_vrr-PRY*I_ESP_D2y_S_M5_vrr+2*oned2z*I_ESP_Py_S_M4_vrr-2*oned2z*I_ESP_Py_S_M5_vrr;
        Double I_ESP_F3z_S_M4_vrr = PAZ*I_ESP_D2z_S_M4_vrr-PRZ*I_ESP_D2z_S_M5_vrr+2*oned2z*I_ESP_Pz_S_M4_vrr-2*oned2z*I_ESP_Pz_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 9 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M4
         * RHS shell quartet name: SQ_ESP_F_S_M5
         * RHS shell quartet name: SQ_ESP_D_S_M4
         * RHS shell quartet name: SQ_ESP_D_S_M5
         ************************************************************/
        Double I_ESP_G4x_S_M4_vrr = PAX*I_ESP_F3x_S_M4_vrr-PRX*I_ESP_F3x_S_M5_vrr+3*oned2z*I_ESP_D2x_S_M4_vrr-3*oned2z*I_ESP_D2x_S_M5_vrr;
        Double I_ESP_G3xy_S_M4_vrr = PAY*I_ESP_F3x_S_M4_vrr-PRY*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_G3xz_S_M4_vrr = PAZ*I_ESP_F3x_S_M4_vrr-PRZ*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_G4y_S_M4_vrr = PAY*I_ESP_F3y_S_M4_vrr-PRY*I_ESP_F3y_S_M5_vrr+3*oned2z*I_ESP_D2y_S_M4_vrr-3*oned2z*I_ESP_D2y_S_M5_vrr;
        Double I_ESP_G3yz_S_M4_vrr = PAZ*I_ESP_F3y_S_M4_vrr-PRZ*I_ESP_F3y_S_M5_vrr;
        Double I_ESP_G4z_S_M4_vrr = PAZ*I_ESP_F3z_S_M4_vrr-PRZ*I_ESP_F3z_S_M5_vrr+3*oned2z*I_ESP_D2z_S_M4_vrr-3*oned2z*I_ESP_D2z_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 11 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M4
         * RHS shell quartet name: SQ_ESP_G_S_M5
         * RHS shell quartet name: SQ_ESP_F_S_M4
         * RHS shell quartet name: SQ_ESP_F_S_M5
         ************************************************************/
        Double I_ESP_H5x_S_M4_vrr = PAX*I_ESP_G4x_S_M4_vrr-PRX*I_ESP_G4x_S_M5_vrr+4*oned2z*I_ESP_F3x_S_M4_vrr-4*oned2z*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_H4xy_S_M4_vrr = PAY*I_ESP_G4x_S_M4_vrr-PRY*I_ESP_G4x_S_M5_vrr;
        Double I_ESP_H4xz_S_M4_vrr = PAZ*I_ESP_G4x_S_M4_vrr-PRZ*I_ESP_G4x_S_M5_vrr;
        Double I_ESP_H3x2y_S_M4_vrr = PAY*I_ESP_G3xy_S_M4_vrr-PRY*I_ESP_G3xy_S_M5_vrr+oned2z*I_ESP_F3x_S_M4_vrr-oned2z*I_ESP_F3x_S_M5_vrr;
        Double I_ESP_Hx4y_S_M4_vrr = PAX*I_ESP_G4y_S_M4_vrr-PRX*I_ESP_G4y_S_M5_vrr;
        Double I_ESP_Hx4z_S_M4_vrr = PAX*I_ESP_G4z_S_M4_vrr-PRX*I_ESP_G4z_S_M5_vrr;
        Double I_ESP_H5y_S_M4_vrr = PAY*I_ESP_G4y_S_M4_vrr-PRY*I_ESP_G4y_S_M5_vrr+4*oned2z*I_ESP_F3y_S_M4_vrr-4*oned2z*I_ESP_F3y_S_M5_vrr;
        Double I_ESP_H4yz_S_M4_vrr = PAZ*I_ESP_G4y_S_M4_vrr-PRZ*I_ESP_G4y_S_M5_vrr;
        Double I_ESP_Hy4z_S_M4_vrr = PAY*I_ESP_G4z_S_M4_vrr-PRY*I_ESP_G4z_S_M5_vrr;
        Double I_ESP_H5z_S_M4_vrr = PAZ*I_ESP_G4z_S_M4_vrr-PRZ*I_ESP_G4z_S_M5_vrr+4*oned2z*I_ESP_F3z_S_M4_vrr-4*oned2z*I_ESP_F3z_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M4
         ************************************************************/
        Double I_ESP_Px_S_M3_vrr = PAX*I_ESP_S_S_M3_vrr-PRX*I_ESP_S_S_M4_vrr;
        Double I_ESP_Py_S_M3_vrr = PAY*I_ESP_S_S_M3_vrr-PRY*I_ESP_S_S_M4_vrr;
        Double I_ESP_Pz_S_M3_vrr = PAZ*I_ESP_S_S_M3_vrr-PRZ*I_ESP_S_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M4
         ************************************************************/
        Double I_ESP_D2x_S_M3_vrr = PAX*I_ESP_Px_S_M3_vrr-PRX*I_ESP_Px_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;
        Double I_ESP_D2y_S_M3_vrr = PAY*I_ESP_Py_S_M3_vrr-PRY*I_ESP_Py_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;
        Double I_ESP_D2z_S_M3_vrr = PAZ*I_ESP_Pz_S_M3_vrr-PRZ*I_ESP_Pz_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 6 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_D_S_M4
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M4
         ************************************************************/
        Double I_ESP_F3x_S_M3_vrr = PAX*I_ESP_D2x_S_M3_vrr-PRX*I_ESP_D2x_S_M4_vrr+2*oned2z*I_ESP_Px_S_M3_vrr-2*oned2z*I_ESP_Px_S_M4_vrr;
        Double I_ESP_F2xy_S_M3_vrr = PAY*I_ESP_D2x_S_M3_vrr-PRY*I_ESP_D2x_S_M4_vrr;
        Double I_ESP_F3y_S_M3_vrr = PAY*I_ESP_D2y_S_M3_vrr-PRY*I_ESP_D2y_S_M4_vrr+2*oned2z*I_ESP_Py_S_M3_vrr-2*oned2z*I_ESP_Py_S_M4_vrr;
        Double I_ESP_F3z_S_M3_vrr = PAZ*I_ESP_D2z_S_M3_vrr-PRZ*I_ESP_D2z_S_M4_vrr+2*oned2z*I_ESP_Pz_S_M3_vrr-2*oned2z*I_ESP_Pz_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M3
         * RHS shell quartet name: SQ_ESP_F_S_M4
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_D_S_M4
         ************************************************************/
        Double I_ESP_G4x_S_M3_vrr = PAX*I_ESP_F3x_S_M3_vrr-PRX*I_ESP_F3x_S_M4_vrr+3*oned2z*I_ESP_D2x_S_M3_vrr-3*oned2z*I_ESP_D2x_S_M4_vrr;
        Double I_ESP_G3xy_S_M3_vrr = PAY*I_ESP_F3x_S_M3_vrr-PRY*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_G3xz_S_M3_vrr = PAZ*I_ESP_F3x_S_M3_vrr-PRZ*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_Gx3y_S_M3_vrr = PAX*I_ESP_F3y_S_M3_vrr-PRX*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_Gx3z_S_M3_vrr = PAX*I_ESP_F3z_S_M3_vrr-PRX*I_ESP_F3z_S_M4_vrr;
        Double I_ESP_G4y_S_M3_vrr = PAY*I_ESP_F3y_S_M3_vrr-PRY*I_ESP_F3y_S_M4_vrr+3*oned2z*I_ESP_D2y_S_M3_vrr-3*oned2z*I_ESP_D2y_S_M4_vrr;
        Double I_ESP_G3yz_S_M3_vrr = PAZ*I_ESP_F3y_S_M3_vrr-PRZ*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_G4z_S_M3_vrr = PAZ*I_ESP_F3z_S_M3_vrr-PRZ*I_ESP_F3z_S_M4_vrr+3*oned2z*I_ESP_D2z_S_M3_vrr-3*oned2z*I_ESP_D2z_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 9 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M3
         * RHS shell quartet name: SQ_ESP_G_S_M4
         * RHS shell quartet name: SQ_ESP_F_S_M3
         * RHS shell quartet name: SQ_ESP_F_S_M4
         ************************************************************/
        Double I_ESP_H5x_S_M3_vrr = PAX*I_ESP_G4x_S_M3_vrr-PRX*I_ESP_G4x_S_M4_vrr+4*oned2z*I_ESP_F3x_S_M3_vrr-4*oned2z*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_H4xy_S_M3_vrr = PAY*I_ESP_G4x_S_M3_vrr-PRY*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_H4xz_S_M3_vrr = PAZ*I_ESP_G4x_S_M3_vrr-PRZ*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_H3x2y_S_M3_vrr = PAY*I_ESP_G3xy_S_M3_vrr-PRY*I_ESP_G3xy_S_M4_vrr+oned2z*I_ESP_F3x_S_M3_vrr-oned2z*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_H3x2z_S_M3_vrr = PAZ*I_ESP_G3xz_S_M3_vrr-PRZ*I_ESP_G3xz_S_M4_vrr+oned2z*I_ESP_F3x_S_M3_vrr-oned2z*I_ESP_F3x_S_M4_vrr;
        Double I_ESP_Hx4y_S_M3_vrr = PAX*I_ESP_G4y_S_M3_vrr-PRX*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_Hx4z_S_M3_vrr = PAX*I_ESP_G4z_S_M3_vrr-PRX*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_H5y_S_M3_vrr = PAY*I_ESP_G4y_S_M3_vrr-PRY*I_ESP_G4y_S_M4_vrr+4*oned2z*I_ESP_F3y_S_M3_vrr-4*oned2z*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_H4yz_S_M3_vrr = PAZ*I_ESP_G4y_S_M3_vrr-PRZ*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_H3y2z_S_M3_vrr = PAZ*I_ESP_G3yz_S_M3_vrr-PRZ*I_ESP_G3yz_S_M4_vrr+oned2z*I_ESP_F3y_S_M3_vrr-oned2z*I_ESP_F3y_S_M4_vrr;
        Double I_ESP_Hy4z_S_M3_vrr = PAY*I_ESP_G4z_S_M3_vrr-PRY*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_H5z_S_M3_vrr = PAZ*I_ESP_G4z_S_M3_vrr-PRZ*I_ESP_G4z_S_M4_vrr+4*oned2z*I_ESP_F3z_S_M3_vrr-4*oned2z*I_ESP_F3z_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 12 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M3
         * RHS shell quartet name: SQ_ESP_H_S_M4
         * RHS shell quartet name: SQ_ESP_G_S_M3
         * RHS shell quartet name: SQ_ESP_G_S_M4
         ************************************************************/
        Double I_ESP_I6x_S_M3_vrr = PAX*I_ESP_H5x_S_M3_vrr-PRX*I_ESP_H5x_S_M4_vrr+5*oned2z*I_ESP_G4x_S_M3_vrr-5*oned2z*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_I5xy_S_M3_vrr = PAY*I_ESP_H5x_S_M3_vrr-PRY*I_ESP_H5x_S_M4_vrr;
        Double I_ESP_I5xz_S_M3_vrr = PAZ*I_ESP_H5x_S_M3_vrr-PRZ*I_ESP_H5x_S_M4_vrr;
        Double I_ESP_I4x2y_S_M3_vrr = PAY*I_ESP_H4xy_S_M3_vrr-PRY*I_ESP_H4xy_S_M4_vrr+oned2z*I_ESP_G4x_S_M3_vrr-oned2z*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_I4x2z_S_M3_vrr = PAZ*I_ESP_H4xz_S_M3_vrr-PRZ*I_ESP_H4xz_S_M4_vrr+oned2z*I_ESP_G4x_S_M3_vrr-oned2z*I_ESP_G4x_S_M4_vrr;
        Double I_ESP_I3x3y_S_M3_vrr = PAY*I_ESP_H3x2y_S_M3_vrr-PRY*I_ESP_H3x2y_S_M4_vrr+2*oned2z*I_ESP_G3xy_S_M3_vrr-2*oned2z*I_ESP_G3xy_S_M4_vrr;
        Double I_ESP_I2x4y_S_M3_vrr = PAX*I_ESP_Hx4y_S_M3_vrr-PRX*I_ESP_Hx4y_S_M4_vrr+oned2z*I_ESP_G4y_S_M3_vrr-oned2z*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_I2x4z_S_M3_vrr = PAX*I_ESP_Hx4z_S_M3_vrr-PRX*I_ESP_Hx4z_S_M4_vrr+oned2z*I_ESP_G4z_S_M3_vrr-oned2z*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_Ix5y_S_M3_vrr = PAX*I_ESP_H5y_S_M3_vrr-PRX*I_ESP_H5y_S_M4_vrr;
        Double I_ESP_Ix5z_S_M3_vrr = PAX*I_ESP_H5z_S_M3_vrr-PRX*I_ESP_H5z_S_M4_vrr;
        Double I_ESP_I6y_S_M3_vrr = PAY*I_ESP_H5y_S_M3_vrr-PRY*I_ESP_H5y_S_M4_vrr+5*oned2z*I_ESP_G4y_S_M3_vrr-5*oned2z*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_I5yz_S_M3_vrr = PAZ*I_ESP_H5y_S_M3_vrr-PRZ*I_ESP_H5y_S_M4_vrr;
        Double I_ESP_I4y2z_S_M3_vrr = PAZ*I_ESP_H4yz_S_M3_vrr-PRZ*I_ESP_H4yz_S_M4_vrr+oned2z*I_ESP_G4y_S_M3_vrr-oned2z*I_ESP_G4y_S_M4_vrr;
        Double I_ESP_I2y4z_S_M3_vrr = PAY*I_ESP_Hy4z_S_M3_vrr-PRY*I_ESP_Hy4z_S_M4_vrr+oned2z*I_ESP_G4z_S_M3_vrr-oned2z*I_ESP_G4z_S_M4_vrr;
        Double I_ESP_Iy5z_S_M3_vrr = PAY*I_ESP_H5z_S_M3_vrr-PRY*I_ESP_H5z_S_M4_vrr;
        Double I_ESP_I6z_S_M3_vrr = PAZ*I_ESP_H5z_S_M3_vrr-PRZ*I_ESP_H5z_S_M4_vrr+5*oned2z*I_ESP_G4z_S_M3_vrr-5*oned2z*I_ESP_G4z_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M3
         ************************************************************/
        Double I_ESP_Px_S_M2_vrr = PAX*I_ESP_S_S_M2_vrr-PRX*I_ESP_S_S_M3_vrr;
        Double I_ESP_Py_S_M2_vrr = PAY*I_ESP_S_S_M2_vrr-PRY*I_ESP_S_S_M3_vrr;
        Double I_ESP_Pz_S_M2_vrr = PAZ*I_ESP_S_S_M2_vrr-PRZ*I_ESP_S_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M3
         ************************************************************/
        Double I_ESP_D2x_S_M2_vrr = PAX*I_ESP_Px_S_M2_vrr-PRX*I_ESP_Px_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;
        Double I_ESP_D2y_S_M2_vrr = PAY*I_ESP_Py_S_M2_vrr-PRY*I_ESP_Py_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;
        Double I_ESP_D2z_S_M2_vrr = PAZ*I_ESP_Pz_S_M2_vrr-PRZ*I_ESP_Pz_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 4 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M3
         ************************************************************/
        Double I_ESP_F3x_S_M2_vrr = PAX*I_ESP_D2x_S_M2_vrr-PRX*I_ESP_D2x_S_M3_vrr+2*oned2z*I_ESP_Px_S_M2_vrr-2*oned2z*I_ESP_Px_S_M3_vrr;
        Double I_ESP_F2xy_S_M2_vrr = PAY*I_ESP_D2x_S_M2_vrr-PRY*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_F2xz_S_M2_vrr = PAZ*I_ESP_D2x_S_M2_vrr-PRZ*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_F3y_S_M2_vrr = PAY*I_ESP_D2y_S_M2_vrr-PRY*I_ESP_D2y_S_M3_vrr+2*oned2z*I_ESP_Py_S_M2_vrr-2*oned2z*I_ESP_Py_S_M3_vrr;
        Double I_ESP_F2yz_S_M2_vrr = PAZ*I_ESP_D2y_S_M2_vrr-PRZ*I_ESP_D2y_S_M3_vrr;
        Double I_ESP_F3z_S_M2_vrr = PAZ*I_ESP_D2z_S_M2_vrr-PRZ*I_ESP_D2z_S_M3_vrr+2*oned2z*I_ESP_Pz_S_M2_vrr-2*oned2z*I_ESP_Pz_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 5 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_F_S_M3
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M3
         ************************************************************/
        Double I_ESP_G4x_S_M2_vrr = PAX*I_ESP_F3x_S_M2_vrr-PRX*I_ESP_F3x_S_M3_vrr+3*oned2z*I_ESP_D2x_S_M2_vrr-3*oned2z*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_G3xy_S_M2_vrr = PAY*I_ESP_F3x_S_M2_vrr-PRY*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_G3xz_S_M2_vrr = PAZ*I_ESP_F3x_S_M2_vrr-PRZ*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_G2x2y_S_M2_vrr = PAY*I_ESP_F2xy_S_M2_vrr-PRY*I_ESP_F2xy_S_M3_vrr+oned2z*I_ESP_D2x_S_M2_vrr-oned2z*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_Gx3y_S_M2_vrr = PAX*I_ESP_F3y_S_M2_vrr-PRX*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_Gx3z_S_M2_vrr = PAX*I_ESP_F3z_S_M2_vrr-PRX*I_ESP_F3z_S_M3_vrr;
        Double I_ESP_G4y_S_M2_vrr = PAY*I_ESP_F3y_S_M2_vrr-PRY*I_ESP_F3y_S_M3_vrr+3*oned2z*I_ESP_D2y_S_M2_vrr-3*oned2z*I_ESP_D2y_S_M3_vrr;
        Double I_ESP_G3yz_S_M2_vrr = PAZ*I_ESP_F3y_S_M2_vrr-PRZ*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_Gy3z_S_M2_vrr = PAY*I_ESP_F3z_S_M2_vrr-PRY*I_ESP_F3z_S_M3_vrr;
        Double I_ESP_G4z_S_M2_vrr = PAZ*I_ESP_F3z_S_M2_vrr-PRZ*I_ESP_F3z_S_M3_vrr+3*oned2z*I_ESP_D2z_S_M2_vrr-3*oned2z*I_ESP_D2z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M2
         * RHS shell quartet name: SQ_ESP_G_S_M3
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_F_S_M3
         ************************************************************/
        Double I_ESP_H5x_S_M2_vrr = PAX*I_ESP_G4x_S_M2_vrr-PRX*I_ESP_G4x_S_M3_vrr+4*oned2z*I_ESP_F3x_S_M2_vrr-4*oned2z*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_H4xy_S_M2_vrr = PAY*I_ESP_G4x_S_M2_vrr-PRY*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_H4xz_S_M2_vrr = PAZ*I_ESP_G4x_S_M2_vrr-PRZ*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_H3x2y_S_M2_vrr = PAY*I_ESP_G3xy_S_M2_vrr-PRY*I_ESP_G3xy_S_M3_vrr+oned2z*I_ESP_F3x_S_M2_vrr-oned2z*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_H3x2z_S_M2_vrr = PAZ*I_ESP_G3xz_S_M2_vrr-PRZ*I_ESP_G3xz_S_M3_vrr+oned2z*I_ESP_F3x_S_M2_vrr-oned2z*I_ESP_F3x_S_M3_vrr;
        Double I_ESP_H2x3y_S_M2_vrr = PAX*I_ESP_Gx3y_S_M2_vrr-PRX*I_ESP_Gx3y_S_M3_vrr+oned2z*I_ESP_F3y_S_M2_vrr-oned2z*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_H2x3z_S_M2_vrr = PAX*I_ESP_Gx3z_S_M2_vrr-PRX*I_ESP_Gx3z_S_M3_vrr+oned2z*I_ESP_F3z_S_M2_vrr-oned2z*I_ESP_F3z_S_M3_vrr;
        Double I_ESP_Hx4y_S_M2_vrr = PAX*I_ESP_G4y_S_M2_vrr-PRX*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_Hx4z_S_M2_vrr = PAX*I_ESP_G4z_S_M2_vrr-PRX*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_H5y_S_M2_vrr = PAY*I_ESP_G4y_S_M2_vrr-PRY*I_ESP_G4y_S_M3_vrr+4*oned2z*I_ESP_F3y_S_M2_vrr-4*oned2z*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_H4yz_S_M2_vrr = PAZ*I_ESP_G4y_S_M2_vrr-PRZ*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_H3y2z_S_M2_vrr = PAZ*I_ESP_G3yz_S_M2_vrr-PRZ*I_ESP_G3yz_S_M3_vrr+oned2z*I_ESP_F3y_S_M2_vrr-oned2z*I_ESP_F3y_S_M3_vrr;
        Double I_ESP_Hy4z_S_M2_vrr = PAY*I_ESP_G4z_S_M2_vrr-PRY*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_H5z_S_M2_vrr = PAZ*I_ESP_G4z_S_M2_vrr-PRZ*I_ESP_G4z_S_M3_vrr+4*oned2z*I_ESP_F3z_S_M2_vrr-4*oned2z*I_ESP_F3z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 10 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M2
         * RHS shell quartet name: SQ_ESP_H_S_M3
         * RHS shell quartet name: SQ_ESP_G_S_M2
         * RHS shell quartet name: SQ_ESP_G_S_M3
         ************************************************************/
        Double I_ESP_I6x_S_M2_vrr = PAX*I_ESP_H5x_S_M2_vrr-PRX*I_ESP_H5x_S_M3_vrr+5*oned2z*I_ESP_G4x_S_M2_vrr-5*oned2z*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_I5xy_S_M2_vrr = PAY*I_ESP_H5x_S_M2_vrr-PRY*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_I5xz_S_M2_vrr = PAZ*I_ESP_H5x_S_M2_vrr-PRZ*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_I4x2y_S_M2_vrr = PAY*I_ESP_H4xy_S_M2_vrr-PRY*I_ESP_H4xy_S_M3_vrr+oned2z*I_ESP_G4x_S_M2_vrr-oned2z*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_I4x2z_S_M2_vrr = PAZ*I_ESP_H4xz_S_M2_vrr-PRZ*I_ESP_H4xz_S_M3_vrr+oned2z*I_ESP_G4x_S_M2_vrr-oned2z*I_ESP_G4x_S_M3_vrr;
        Double I_ESP_I3x3y_S_M2_vrr = PAY*I_ESP_H3x2y_S_M2_vrr-PRY*I_ESP_H3x2y_S_M3_vrr+2*oned2z*I_ESP_G3xy_S_M2_vrr-2*oned2z*I_ESP_G3xy_S_M3_vrr;
        Double I_ESP_I3x3z_S_M2_vrr = PAZ*I_ESP_H3x2z_S_M2_vrr-PRZ*I_ESP_H3x2z_S_M3_vrr+2*oned2z*I_ESP_G3xz_S_M2_vrr-2*oned2z*I_ESP_G3xz_S_M3_vrr;
        Double I_ESP_I2x4y_S_M2_vrr = PAX*I_ESP_Hx4y_S_M2_vrr-PRX*I_ESP_Hx4y_S_M3_vrr+oned2z*I_ESP_G4y_S_M2_vrr-oned2z*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_I2x4z_S_M2_vrr = PAX*I_ESP_Hx4z_S_M2_vrr-PRX*I_ESP_Hx4z_S_M3_vrr+oned2z*I_ESP_G4z_S_M2_vrr-oned2z*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_Ix5y_S_M2_vrr = PAX*I_ESP_H5y_S_M2_vrr-PRX*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_Ix5z_S_M2_vrr = PAX*I_ESP_H5z_S_M2_vrr-PRX*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_I6y_S_M2_vrr = PAY*I_ESP_H5y_S_M2_vrr-PRY*I_ESP_H5y_S_M3_vrr+5*oned2z*I_ESP_G4y_S_M2_vrr-5*oned2z*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_I5yz_S_M2_vrr = PAZ*I_ESP_H5y_S_M2_vrr-PRZ*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_I4y2z_S_M2_vrr = PAZ*I_ESP_H4yz_S_M2_vrr-PRZ*I_ESP_H4yz_S_M3_vrr+oned2z*I_ESP_G4y_S_M2_vrr-oned2z*I_ESP_G4y_S_M3_vrr;
        Double I_ESP_I3y3z_S_M2_vrr = PAZ*I_ESP_H3y2z_S_M2_vrr-PRZ*I_ESP_H3y2z_S_M3_vrr+2*oned2z*I_ESP_G3yz_S_M2_vrr-2*oned2z*I_ESP_G3yz_S_M3_vrr;
        Double I_ESP_I2y4z_S_M2_vrr = PAY*I_ESP_Hy4z_S_M2_vrr-PRY*I_ESP_Hy4z_S_M3_vrr+oned2z*I_ESP_G4z_S_M2_vrr-oned2z*I_ESP_G4z_S_M3_vrr;
        Double I_ESP_Iy5z_S_M2_vrr = PAY*I_ESP_H5z_S_M2_vrr-PRY*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_I6z_S_M2_vrr = PAZ*I_ESP_H5z_S_M2_vrr-PRZ*I_ESP_H5z_S_M3_vrr+5*oned2z*I_ESP_G4z_S_M2_vrr-5*oned2z*I_ESP_G4z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 14 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S_M2
         * RHS shell quartet name: SQ_ESP_I_S_M3
         * RHS shell quartet name: SQ_ESP_H_S_M2
         * RHS shell quartet name: SQ_ESP_H_S_M3
         ************************************************************/
        Double I_ESP_K7x_S_M2_vrr = PAX*I_ESP_I6x_S_M2_vrr-PRX*I_ESP_I6x_S_M3_vrr+6*oned2z*I_ESP_H5x_S_M2_vrr-6*oned2z*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_K6xy_S_M2_vrr = PAY*I_ESP_I6x_S_M2_vrr-PRY*I_ESP_I6x_S_M3_vrr;
        Double I_ESP_K6xz_S_M2_vrr = PAZ*I_ESP_I6x_S_M2_vrr-PRZ*I_ESP_I6x_S_M3_vrr;
        Double I_ESP_K5x2y_S_M2_vrr = PAY*I_ESP_I5xy_S_M2_vrr-PRY*I_ESP_I5xy_S_M3_vrr+oned2z*I_ESP_H5x_S_M2_vrr-oned2z*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_K5x2z_S_M2_vrr = PAZ*I_ESP_I5xz_S_M2_vrr-PRZ*I_ESP_I5xz_S_M3_vrr+oned2z*I_ESP_H5x_S_M2_vrr-oned2z*I_ESP_H5x_S_M3_vrr;
        Double I_ESP_K4x3y_S_M2_vrr = PAY*I_ESP_I4x2y_S_M2_vrr-PRY*I_ESP_I4x2y_S_M3_vrr+2*oned2z*I_ESP_H4xy_S_M2_vrr-2*oned2z*I_ESP_H4xy_S_M3_vrr;
        Double I_ESP_K4x3z_S_M2_vrr = PAZ*I_ESP_I4x2z_S_M2_vrr-PRZ*I_ESP_I4x2z_S_M3_vrr+2*oned2z*I_ESP_H4xz_S_M2_vrr-2*oned2z*I_ESP_H4xz_S_M3_vrr;
        Double I_ESP_K3x4y_S_M2_vrr = PAX*I_ESP_I2x4y_S_M2_vrr-PRX*I_ESP_I2x4y_S_M3_vrr+2*oned2z*I_ESP_Hx4y_S_M2_vrr-2*oned2z*I_ESP_Hx4y_S_M3_vrr;
        Double I_ESP_K3x3yz_S_M2_vrr = PAZ*I_ESP_I3x3y_S_M2_vrr-PRZ*I_ESP_I3x3y_S_M3_vrr;
        Double I_ESP_K3x4z_S_M2_vrr = PAX*I_ESP_I2x4z_S_M2_vrr-PRX*I_ESP_I2x4z_S_M3_vrr+2*oned2z*I_ESP_Hx4z_S_M2_vrr-2*oned2z*I_ESP_Hx4z_S_M3_vrr;
        Double I_ESP_K2x5y_S_M2_vrr = PAX*I_ESP_Ix5y_S_M2_vrr-PRX*I_ESP_Ix5y_S_M3_vrr+oned2z*I_ESP_H5y_S_M2_vrr-oned2z*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_K2x5z_S_M2_vrr = PAX*I_ESP_Ix5z_S_M2_vrr-PRX*I_ESP_Ix5z_S_M3_vrr+oned2z*I_ESP_H5z_S_M2_vrr-oned2z*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_Kx6y_S_M2_vrr = PAX*I_ESP_I6y_S_M2_vrr-PRX*I_ESP_I6y_S_M3_vrr;
        Double I_ESP_Kx6z_S_M2_vrr = PAX*I_ESP_I6z_S_M2_vrr-PRX*I_ESP_I6z_S_M3_vrr;
        Double I_ESP_K7y_S_M2_vrr = PAY*I_ESP_I6y_S_M2_vrr-PRY*I_ESP_I6y_S_M3_vrr+6*oned2z*I_ESP_H5y_S_M2_vrr-6*oned2z*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_K6yz_S_M2_vrr = PAZ*I_ESP_I6y_S_M2_vrr-PRZ*I_ESP_I6y_S_M3_vrr;
        Double I_ESP_K5y2z_S_M2_vrr = PAZ*I_ESP_I5yz_S_M2_vrr-PRZ*I_ESP_I5yz_S_M3_vrr+oned2z*I_ESP_H5y_S_M2_vrr-oned2z*I_ESP_H5y_S_M3_vrr;
        Double I_ESP_K4y3z_S_M2_vrr = PAZ*I_ESP_I4y2z_S_M2_vrr-PRZ*I_ESP_I4y2z_S_M3_vrr+2*oned2z*I_ESP_H4yz_S_M2_vrr-2*oned2z*I_ESP_H4yz_S_M3_vrr;
        Double I_ESP_K3y4z_S_M2_vrr = PAY*I_ESP_I2y4z_S_M2_vrr-PRY*I_ESP_I2y4z_S_M3_vrr+2*oned2z*I_ESP_Hy4z_S_M2_vrr-2*oned2z*I_ESP_Hy4z_S_M3_vrr;
        Double I_ESP_K2y5z_S_M2_vrr = PAY*I_ESP_Iy5z_S_M2_vrr-PRY*I_ESP_Iy5z_S_M3_vrr+oned2z*I_ESP_H5z_S_M2_vrr-oned2z*I_ESP_H5z_S_M3_vrr;
        Double I_ESP_Ky6z_S_M2_vrr = PAY*I_ESP_I6z_S_M2_vrr-PRY*I_ESP_I6z_S_M3_vrr;
        Double I_ESP_K7z_S_M2_vrr = PAZ*I_ESP_I6z_S_M2_vrr-PRZ*I_ESP_I6z_S_M3_vrr+6*oned2z*I_ESP_H5z_S_M2_vrr-6*oned2z*I_ESP_H5z_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_Px_S_M1_vrr = PAX*I_ESP_S_S_M1_vrr-PRX*I_ESP_S_S_M2_vrr;
        Double I_ESP_Py_S_M1_vrr = PAY*I_ESP_S_S_M1_vrr-PRY*I_ESP_S_S_M2_vrr;
        Double I_ESP_Pz_S_M1_vrr = PAZ*I_ESP_S_S_M1_vrr-PRZ*I_ESP_S_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_D2x_S_M1_vrr = PAX*I_ESP_Px_S_M1_vrr-PRX*I_ESP_Px_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_Dxy_S_M1_vrr = PAY*I_ESP_Px_S_M1_vrr-PRY*I_ESP_Px_S_M2_vrr;
        Double I_ESP_D2y_S_M1_vrr = PAY*I_ESP_Py_S_M1_vrr-PRY*I_ESP_Py_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_D2z_S_M1_vrr = PAZ*I_ESP_Pz_S_M1_vrr-PRZ*I_ESP_Pz_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         ************************************************************/
        Double I_ESP_F3x_S_M1_vrr = PAX*I_ESP_D2x_S_M1_vrr-PRX*I_ESP_D2x_S_M2_vrr+2*oned2z*I_ESP_Px_S_M1_vrr-2*oned2z*I_ESP_Px_S_M2_vrr;
        Double I_ESP_F2xy_S_M1_vrr = PAY*I_ESP_D2x_S_M1_vrr-PRY*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_F2xz_S_M1_vrr = PAZ*I_ESP_D2x_S_M1_vrr-PRZ*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_Fx2y_S_M1_vrr = PAX*I_ESP_D2y_S_M1_vrr-PRX*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_Fx2z_S_M1_vrr = PAX*I_ESP_D2z_S_M1_vrr-PRX*I_ESP_D2z_S_M2_vrr;
        Double I_ESP_F3y_S_M1_vrr = PAY*I_ESP_D2y_S_M1_vrr-PRY*I_ESP_D2y_S_M2_vrr+2*oned2z*I_ESP_Py_S_M1_vrr-2*oned2z*I_ESP_Py_S_M2_vrr;
        Double I_ESP_F2yz_S_M1_vrr = PAZ*I_ESP_D2y_S_M1_vrr-PRZ*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_F3z_S_M1_vrr = PAZ*I_ESP_D2z_S_M1_vrr-PRZ*I_ESP_D2z_S_M2_vrr+2*oned2z*I_ESP_Pz_S_M1_vrr-2*oned2z*I_ESP_Pz_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         ************************************************************/
        Double I_ESP_G4x_S_M1_vrr = PAX*I_ESP_F3x_S_M1_vrr-PRX*I_ESP_F3x_S_M2_vrr+3*oned2z*I_ESP_D2x_S_M1_vrr-3*oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_G3xy_S_M1_vrr = PAY*I_ESP_F3x_S_M1_vrr-PRY*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_G3xz_S_M1_vrr = PAZ*I_ESP_F3x_S_M1_vrr-PRZ*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_G2x2y_S_M1_vrr = PAY*I_ESP_F2xy_S_M1_vrr-PRY*I_ESP_F2xy_S_M2_vrr+oned2z*I_ESP_D2x_S_M1_vrr-oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_G2x2z_S_M1_vrr = PAZ*I_ESP_F2xz_S_M1_vrr-PRZ*I_ESP_F2xz_S_M2_vrr+oned2z*I_ESP_D2x_S_M1_vrr-oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_Gx3y_S_M1_vrr = PAX*I_ESP_F3y_S_M1_vrr-PRX*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_Gx3z_S_M1_vrr = PAX*I_ESP_F3z_S_M1_vrr-PRX*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_G4y_S_M1_vrr = PAY*I_ESP_F3y_S_M1_vrr-PRY*I_ESP_F3y_S_M2_vrr+3*oned2z*I_ESP_D2y_S_M1_vrr-3*oned2z*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_G3yz_S_M1_vrr = PAZ*I_ESP_F3y_S_M1_vrr-PRZ*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_G2y2z_S_M1_vrr = PAZ*I_ESP_F2yz_S_M1_vrr-PRZ*I_ESP_F2yz_S_M2_vrr+oned2z*I_ESP_D2y_S_M1_vrr-oned2z*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_Gy3z_S_M1_vrr = PAY*I_ESP_F3z_S_M1_vrr-PRY*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_G4z_S_M1_vrr = PAZ*I_ESP_F3z_S_M1_vrr-PRZ*I_ESP_F3z_S_M2_vrr+3*oned2z*I_ESP_D2z_S_M1_vrr-3*oned2z*I_ESP_D2z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 5 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_G_S_M2
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_F_S_M2
         ************************************************************/
        Double I_ESP_H5x_S_M1_vrr = PAX*I_ESP_G4x_S_M1_vrr-PRX*I_ESP_G4x_S_M2_vrr+4*oned2z*I_ESP_F3x_S_M1_vrr-4*oned2z*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_H4xy_S_M1_vrr = PAY*I_ESP_G4x_S_M1_vrr-PRY*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_H4xz_S_M1_vrr = PAZ*I_ESP_G4x_S_M1_vrr-PRZ*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_H3x2y_S_M1_vrr = PAY*I_ESP_G3xy_S_M1_vrr-PRY*I_ESP_G3xy_S_M2_vrr+oned2z*I_ESP_F3x_S_M1_vrr-oned2z*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_H3x2z_S_M1_vrr = PAZ*I_ESP_G3xz_S_M1_vrr-PRZ*I_ESP_G3xz_S_M2_vrr+oned2z*I_ESP_F3x_S_M1_vrr-oned2z*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_H2x3y_S_M1_vrr = PAX*I_ESP_Gx3y_S_M1_vrr-PRX*I_ESP_Gx3y_S_M2_vrr+oned2z*I_ESP_F3y_S_M1_vrr-oned2z*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_H2x2yz_S_M1_vrr = PAZ*I_ESP_G2x2y_S_M1_vrr-PRZ*I_ESP_G2x2y_S_M2_vrr;
        Double I_ESP_H2x3z_S_M1_vrr = PAX*I_ESP_Gx3z_S_M1_vrr-PRX*I_ESP_Gx3z_S_M2_vrr+oned2z*I_ESP_F3z_S_M1_vrr-oned2z*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_Hx4y_S_M1_vrr = PAX*I_ESP_G4y_S_M1_vrr-PRX*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_Hx4z_S_M1_vrr = PAX*I_ESP_G4z_S_M1_vrr-PRX*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_H5y_S_M1_vrr = PAY*I_ESP_G4y_S_M1_vrr-PRY*I_ESP_G4y_S_M2_vrr+4*oned2z*I_ESP_F3y_S_M1_vrr-4*oned2z*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_H4yz_S_M1_vrr = PAZ*I_ESP_G4y_S_M1_vrr-PRZ*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_H3y2z_S_M1_vrr = PAZ*I_ESP_G3yz_S_M1_vrr-PRZ*I_ESP_G3yz_S_M2_vrr+oned2z*I_ESP_F3y_S_M1_vrr-oned2z*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_H2y3z_S_M1_vrr = PAY*I_ESP_Gy3z_S_M1_vrr-PRY*I_ESP_Gy3z_S_M2_vrr+oned2z*I_ESP_F3z_S_M1_vrr-oned2z*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_Hy4z_S_M1_vrr = PAY*I_ESP_G4z_S_M1_vrr-PRY*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_H5z_S_M1_vrr = PAZ*I_ESP_G4z_S_M1_vrr-PRZ*I_ESP_G4z_S_M2_vrr+4*oned2z*I_ESP_F3z_S_M1_vrr-4*oned2z*I_ESP_F3z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M1
         * RHS shell quartet name: SQ_ESP_H_S_M2
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_G_S_M2
         ************************************************************/
        Double I_ESP_I6x_S_M1_vrr = PAX*I_ESP_H5x_S_M1_vrr-PRX*I_ESP_H5x_S_M2_vrr+5*oned2z*I_ESP_G4x_S_M1_vrr-5*oned2z*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_I5xy_S_M1_vrr = PAY*I_ESP_H5x_S_M1_vrr-PRY*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_I5xz_S_M1_vrr = PAZ*I_ESP_H5x_S_M1_vrr-PRZ*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_I4x2y_S_M1_vrr = PAY*I_ESP_H4xy_S_M1_vrr-PRY*I_ESP_H4xy_S_M2_vrr+oned2z*I_ESP_G4x_S_M1_vrr-oned2z*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_I4x2z_S_M1_vrr = PAZ*I_ESP_H4xz_S_M1_vrr-PRZ*I_ESP_H4xz_S_M2_vrr+oned2z*I_ESP_G4x_S_M1_vrr-oned2z*I_ESP_G4x_S_M2_vrr;
        Double I_ESP_I3x3y_S_M1_vrr = PAY*I_ESP_H3x2y_S_M1_vrr-PRY*I_ESP_H3x2y_S_M2_vrr+2*oned2z*I_ESP_G3xy_S_M1_vrr-2*oned2z*I_ESP_G3xy_S_M2_vrr;
        Double I_ESP_I3x2yz_S_M1_vrr = PAZ*I_ESP_H3x2y_S_M1_vrr-PRZ*I_ESP_H3x2y_S_M2_vrr;
        Double I_ESP_I3x3z_S_M1_vrr = PAZ*I_ESP_H3x2z_S_M1_vrr-PRZ*I_ESP_H3x2z_S_M2_vrr+2*oned2z*I_ESP_G3xz_S_M1_vrr-2*oned2z*I_ESP_G3xz_S_M2_vrr;
        Double I_ESP_I2x4y_S_M1_vrr = PAX*I_ESP_Hx4y_S_M1_vrr-PRX*I_ESP_Hx4y_S_M2_vrr+oned2z*I_ESP_G4y_S_M1_vrr-oned2z*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_I2x3yz_S_M1_vrr = PAZ*I_ESP_H2x3y_S_M1_vrr-PRZ*I_ESP_H2x3y_S_M2_vrr;
        Double I_ESP_I2xy3z_S_M1_vrr = PAY*I_ESP_H2x3z_S_M1_vrr-PRY*I_ESP_H2x3z_S_M2_vrr;
        Double I_ESP_I2x4z_S_M1_vrr = PAX*I_ESP_Hx4z_S_M1_vrr-PRX*I_ESP_Hx4z_S_M2_vrr+oned2z*I_ESP_G4z_S_M1_vrr-oned2z*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_Ix5y_S_M1_vrr = PAX*I_ESP_H5y_S_M1_vrr-PRX*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_Ix5z_S_M1_vrr = PAX*I_ESP_H5z_S_M1_vrr-PRX*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_I6y_S_M1_vrr = PAY*I_ESP_H5y_S_M1_vrr-PRY*I_ESP_H5y_S_M2_vrr+5*oned2z*I_ESP_G4y_S_M1_vrr-5*oned2z*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_I5yz_S_M1_vrr = PAZ*I_ESP_H5y_S_M1_vrr-PRZ*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_I4y2z_S_M1_vrr = PAZ*I_ESP_H4yz_S_M1_vrr-PRZ*I_ESP_H4yz_S_M2_vrr+oned2z*I_ESP_G4y_S_M1_vrr-oned2z*I_ESP_G4y_S_M2_vrr;
        Double I_ESP_I3y3z_S_M1_vrr = PAZ*I_ESP_H3y2z_S_M1_vrr-PRZ*I_ESP_H3y2z_S_M2_vrr+2*oned2z*I_ESP_G3yz_S_M1_vrr-2*oned2z*I_ESP_G3yz_S_M2_vrr;
        Double I_ESP_I2y4z_S_M1_vrr = PAY*I_ESP_Hy4z_S_M1_vrr-PRY*I_ESP_Hy4z_S_M2_vrr+oned2z*I_ESP_G4z_S_M1_vrr-oned2z*I_ESP_G4z_S_M2_vrr;
        Double I_ESP_Iy5z_S_M1_vrr = PAY*I_ESP_H5z_S_M1_vrr-PRY*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_I6z_S_M1_vrr = PAZ*I_ESP_H5z_S_M1_vrr-PRZ*I_ESP_H5z_S_M2_vrr+5*oned2z*I_ESP_G4z_S_M1_vrr-5*oned2z*I_ESP_G4z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 9 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S_M1
         * RHS shell quartet name: SQ_ESP_I_S_M2
         * RHS shell quartet name: SQ_ESP_H_S_M1
         * RHS shell quartet name: SQ_ESP_H_S_M2
         ************************************************************/
        Double I_ESP_K7x_S_M1_vrr = PAX*I_ESP_I6x_S_M1_vrr-PRX*I_ESP_I6x_S_M2_vrr+6*oned2z*I_ESP_H5x_S_M1_vrr-6*oned2z*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_K6xy_S_M1_vrr = PAY*I_ESP_I6x_S_M1_vrr-PRY*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_K6xz_S_M1_vrr = PAZ*I_ESP_I6x_S_M1_vrr-PRZ*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_K5x2y_S_M1_vrr = PAY*I_ESP_I5xy_S_M1_vrr-PRY*I_ESP_I5xy_S_M2_vrr+oned2z*I_ESP_H5x_S_M1_vrr-oned2z*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_K5x2z_S_M1_vrr = PAZ*I_ESP_I5xz_S_M1_vrr-PRZ*I_ESP_I5xz_S_M2_vrr+oned2z*I_ESP_H5x_S_M1_vrr-oned2z*I_ESP_H5x_S_M2_vrr;
        Double I_ESP_K4x3y_S_M1_vrr = PAY*I_ESP_I4x2y_S_M1_vrr-PRY*I_ESP_I4x2y_S_M2_vrr+2*oned2z*I_ESP_H4xy_S_M1_vrr-2*oned2z*I_ESP_H4xy_S_M2_vrr;
        Double I_ESP_K4x2yz_S_M1_vrr = PAZ*I_ESP_I4x2y_S_M1_vrr-PRZ*I_ESP_I4x2y_S_M2_vrr;
        Double I_ESP_K4x3z_S_M1_vrr = PAZ*I_ESP_I4x2z_S_M1_vrr-PRZ*I_ESP_I4x2z_S_M2_vrr+2*oned2z*I_ESP_H4xz_S_M1_vrr-2*oned2z*I_ESP_H4xz_S_M2_vrr;
        Double I_ESP_K3x4y_S_M1_vrr = PAX*I_ESP_I2x4y_S_M1_vrr-PRX*I_ESP_I2x4y_S_M2_vrr+2*oned2z*I_ESP_Hx4y_S_M1_vrr-2*oned2z*I_ESP_Hx4y_S_M2_vrr;
        Double I_ESP_K3x3yz_S_M1_vrr = PAZ*I_ESP_I3x3y_S_M1_vrr-PRZ*I_ESP_I3x3y_S_M2_vrr;
        Double I_ESP_K3xy3z_S_M1_vrr = PAY*I_ESP_I3x3z_S_M1_vrr-PRY*I_ESP_I3x3z_S_M2_vrr;
        Double I_ESP_K3x4z_S_M1_vrr = PAX*I_ESP_I2x4z_S_M1_vrr-PRX*I_ESP_I2x4z_S_M2_vrr+2*oned2z*I_ESP_Hx4z_S_M1_vrr-2*oned2z*I_ESP_Hx4z_S_M2_vrr;
        Double I_ESP_K2x5y_S_M1_vrr = PAX*I_ESP_Ix5y_S_M1_vrr-PRX*I_ESP_Ix5y_S_M2_vrr+oned2z*I_ESP_H5y_S_M1_vrr-oned2z*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_K2x4yz_S_M1_vrr = PAZ*I_ESP_I2x4y_S_M1_vrr-PRZ*I_ESP_I2x4y_S_M2_vrr;
        Double I_ESP_K2xy4z_S_M1_vrr = PAY*I_ESP_I2x4z_S_M1_vrr-PRY*I_ESP_I2x4z_S_M2_vrr;
        Double I_ESP_K2x5z_S_M1_vrr = PAX*I_ESP_Ix5z_S_M1_vrr-PRX*I_ESP_Ix5z_S_M2_vrr+oned2z*I_ESP_H5z_S_M1_vrr-oned2z*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_Kx6y_S_M1_vrr = PAX*I_ESP_I6y_S_M1_vrr-PRX*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_Kx3y3z_S_M1_vrr = PAX*I_ESP_I3y3z_S_M1_vrr-PRX*I_ESP_I3y3z_S_M2_vrr;
        Double I_ESP_Kx6z_S_M1_vrr = PAX*I_ESP_I6z_S_M1_vrr-PRX*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_K7y_S_M1_vrr = PAY*I_ESP_I6y_S_M1_vrr-PRY*I_ESP_I6y_S_M2_vrr+6*oned2z*I_ESP_H5y_S_M1_vrr-6*oned2z*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_K6yz_S_M1_vrr = PAZ*I_ESP_I6y_S_M1_vrr-PRZ*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_K5y2z_S_M1_vrr = PAZ*I_ESP_I5yz_S_M1_vrr-PRZ*I_ESP_I5yz_S_M2_vrr+oned2z*I_ESP_H5y_S_M1_vrr-oned2z*I_ESP_H5y_S_M2_vrr;
        Double I_ESP_K4y3z_S_M1_vrr = PAZ*I_ESP_I4y2z_S_M1_vrr-PRZ*I_ESP_I4y2z_S_M2_vrr+2*oned2z*I_ESP_H4yz_S_M1_vrr-2*oned2z*I_ESP_H4yz_S_M2_vrr;
        Double I_ESP_K3y4z_S_M1_vrr = PAY*I_ESP_I2y4z_S_M1_vrr-PRY*I_ESP_I2y4z_S_M2_vrr+2*oned2z*I_ESP_Hy4z_S_M1_vrr-2*oned2z*I_ESP_Hy4z_S_M2_vrr;
        Double I_ESP_K2y5z_S_M1_vrr = PAY*I_ESP_Iy5z_S_M1_vrr-PRY*I_ESP_Iy5z_S_M2_vrr+oned2z*I_ESP_H5z_S_M1_vrr-oned2z*I_ESP_H5z_S_M2_vrr;
        Double I_ESP_Ky6z_S_M1_vrr = PAY*I_ESP_I6z_S_M1_vrr-PRY*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_K7z_S_M1_vrr = PAZ*I_ESP_I6z_S_M1_vrr-PRZ*I_ESP_I6z_S_M2_vrr+6*oned2z*I_ESP_H5z_S_M1_vrr-6*oned2z*I_ESP_H5z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_L_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 11 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_K_S_M1
         * RHS shell quartet name: SQ_ESP_K_S_M2
         * RHS shell quartet name: SQ_ESP_I_S_M1
         * RHS shell quartet name: SQ_ESP_I_S_M2
         ************************************************************/
        Double I_ESP_L8x_S_M1_vrr = PAX*I_ESP_K7x_S_M1_vrr-PRX*I_ESP_K7x_S_M2_vrr+7*oned2z*I_ESP_I6x_S_M1_vrr-7*oned2z*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_L7xy_S_M1_vrr = PAY*I_ESP_K7x_S_M1_vrr-PRY*I_ESP_K7x_S_M2_vrr;
        Double I_ESP_L7xz_S_M1_vrr = PAZ*I_ESP_K7x_S_M1_vrr-PRZ*I_ESP_K7x_S_M2_vrr;
        Double I_ESP_L6x2y_S_M1_vrr = PAY*I_ESP_K6xy_S_M1_vrr-PRY*I_ESP_K6xy_S_M2_vrr+oned2z*I_ESP_I6x_S_M1_vrr-oned2z*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_L6x2z_S_M1_vrr = PAZ*I_ESP_K6xz_S_M1_vrr-PRZ*I_ESP_K6xz_S_M2_vrr+oned2z*I_ESP_I6x_S_M1_vrr-oned2z*I_ESP_I6x_S_M2_vrr;
        Double I_ESP_L5x3y_S_M1_vrr = PAY*I_ESP_K5x2y_S_M1_vrr-PRY*I_ESP_K5x2y_S_M2_vrr+2*oned2z*I_ESP_I5xy_S_M1_vrr-2*oned2z*I_ESP_I5xy_S_M2_vrr;
        Double I_ESP_L5x2yz_S_M1_vrr = PAZ*I_ESP_K5x2y_S_M1_vrr-PRZ*I_ESP_K5x2y_S_M2_vrr;
        Double I_ESP_L5x3z_S_M1_vrr = PAZ*I_ESP_K5x2z_S_M1_vrr-PRZ*I_ESP_K5x2z_S_M2_vrr+2*oned2z*I_ESP_I5xz_S_M1_vrr-2*oned2z*I_ESP_I5xz_S_M2_vrr;
        Double I_ESP_L4x4y_S_M1_vrr = PAY*I_ESP_K4x3y_S_M1_vrr-PRY*I_ESP_K4x3y_S_M2_vrr+3*oned2z*I_ESP_I4x2y_S_M1_vrr-3*oned2z*I_ESP_I4x2y_S_M2_vrr;
        Double I_ESP_L4x3yz_S_M1_vrr = PAZ*I_ESP_K4x3y_S_M1_vrr-PRZ*I_ESP_K4x3y_S_M2_vrr;
        Double I_ESP_L4xy3z_S_M1_vrr = PAY*I_ESP_K4x3z_S_M1_vrr-PRY*I_ESP_K4x3z_S_M2_vrr;
        Double I_ESP_L4x4z_S_M1_vrr = PAZ*I_ESP_K4x3z_S_M1_vrr-PRZ*I_ESP_K4x3z_S_M2_vrr+3*oned2z*I_ESP_I4x2z_S_M1_vrr-3*oned2z*I_ESP_I4x2z_S_M2_vrr;
        Double I_ESP_L3x5y_S_M1_vrr = PAX*I_ESP_K2x5y_S_M1_vrr-PRX*I_ESP_K2x5y_S_M2_vrr+2*oned2z*I_ESP_Ix5y_S_M1_vrr-2*oned2z*I_ESP_Ix5y_S_M2_vrr;
        Double I_ESP_L3x4yz_S_M1_vrr = PAZ*I_ESP_K3x4y_S_M1_vrr-PRZ*I_ESP_K3x4y_S_M2_vrr;
        Double I_ESP_L3x3y2z_S_M1_vrr = PAZ*I_ESP_K3x3yz_S_M1_vrr-PRZ*I_ESP_K3x3yz_S_M2_vrr+oned2z*I_ESP_I3x3y_S_M1_vrr-oned2z*I_ESP_I3x3y_S_M2_vrr;
        Double I_ESP_L3xy4z_S_M1_vrr = PAY*I_ESP_K3x4z_S_M1_vrr-PRY*I_ESP_K3x4z_S_M2_vrr;
        Double I_ESP_L3x5z_S_M1_vrr = PAX*I_ESP_K2x5z_S_M1_vrr-PRX*I_ESP_K2x5z_S_M2_vrr+2*oned2z*I_ESP_Ix5z_S_M1_vrr-2*oned2z*I_ESP_Ix5z_S_M2_vrr;
        Double I_ESP_L2x6y_S_M1_vrr = PAX*I_ESP_Kx6y_S_M1_vrr-PRX*I_ESP_Kx6y_S_M2_vrr+oned2z*I_ESP_I6y_S_M1_vrr-oned2z*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_L2x5yz_S_M1_vrr = PAZ*I_ESP_K2x5y_S_M1_vrr-PRZ*I_ESP_K2x5y_S_M2_vrr;
        Double I_ESP_L2xy5z_S_M1_vrr = PAY*I_ESP_K2x5z_S_M1_vrr-PRY*I_ESP_K2x5z_S_M2_vrr;
        Double I_ESP_L2x6z_S_M1_vrr = PAX*I_ESP_Kx6z_S_M1_vrr-PRX*I_ESP_Kx6z_S_M2_vrr+oned2z*I_ESP_I6z_S_M1_vrr-oned2z*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_Lx7y_S_M1_vrr = PAX*I_ESP_K7y_S_M1_vrr-PRX*I_ESP_K7y_S_M2_vrr;
        Double I_ESP_Lx4y3z_S_M1_vrr = PAX*I_ESP_K4y3z_S_M1_vrr-PRX*I_ESP_K4y3z_S_M2_vrr;
        Double I_ESP_Lx3y4z_S_M1_vrr = PAX*I_ESP_K3y4z_S_M1_vrr-PRX*I_ESP_K3y4z_S_M2_vrr;
        Double I_ESP_Lx7z_S_M1_vrr = PAX*I_ESP_K7z_S_M1_vrr-PRX*I_ESP_K7z_S_M2_vrr;
        Double I_ESP_L8y_S_M1_vrr = PAY*I_ESP_K7y_S_M1_vrr-PRY*I_ESP_K7y_S_M2_vrr+7*oned2z*I_ESP_I6y_S_M1_vrr-7*oned2z*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_L7yz_S_M1_vrr = PAZ*I_ESP_K7y_S_M1_vrr-PRZ*I_ESP_K7y_S_M2_vrr;
        Double I_ESP_L6y2z_S_M1_vrr = PAZ*I_ESP_K6yz_S_M1_vrr-PRZ*I_ESP_K6yz_S_M2_vrr+oned2z*I_ESP_I6y_S_M1_vrr-oned2z*I_ESP_I6y_S_M2_vrr;
        Double I_ESP_L5y3z_S_M1_vrr = PAZ*I_ESP_K5y2z_S_M1_vrr-PRZ*I_ESP_K5y2z_S_M2_vrr+2*oned2z*I_ESP_I5yz_S_M1_vrr-2*oned2z*I_ESP_I5yz_S_M2_vrr;
        Double I_ESP_L4y4z_S_M1_vrr = PAZ*I_ESP_K4y3z_S_M1_vrr-PRZ*I_ESP_K4y3z_S_M2_vrr+3*oned2z*I_ESP_I4y2z_S_M1_vrr-3*oned2z*I_ESP_I4y2z_S_M2_vrr;
        Double I_ESP_L3y5z_S_M1_vrr = PAY*I_ESP_K2y5z_S_M1_vrr-PRY*I_ESP_K2y5z_S_M2_vrr+2*oned2z*I_ESP_Iy5z_S_M1_vrr-2*oned2z*I_ESP_Iy5z_S_M2_vrr;
        Double I_ESP_L2y6z_S_M1_vrr = PAY*I_ESP_Ky6z_S_M1_vrr-PRY*I_ESP_Ky6z_S_M2_vrr+oned2z*I_ESP_I6z_S_M1_vrr-oned2z*I_ESP_I6z_S_M2_vrr;
        Double I_ESP_Ly7z_S_M1_vrr = PAY*I_ESP_K7z_S_M1_vrr-PRY*I_ESP_K7z_S_M2_vrr;
        Double I_ESP_L8z_S_M1_vrr = PAZ*I_ESP_K7z_S_M1_vrr-PRZ*I_ESP_K7z_S_M2_vrr+7*oned2z*I_ESP_I6z_S_M1_vrr-7*oned2z*I_ESP_I6z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_Px_S_vrr = PAX*I_ESP_S_S_vrr-PRX*I_ESP_S_S_M1_vrr;
        Double I_ESP_Py_S_vrr = PAY*I_ESP_S_S_vrr-PRY*I_ESP_S_S_M1_vrr;
        Double I_ESP_Pz_S_vrr = PAZ*I_ESP_S_S_vrr-PRZ*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_D2x_S_vrr = PAX*I_ESP_Px_S_vrr-PRX*I_ESP_Px_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_Dxy_S_vrr = PAY*I_ESP_Px_S_vrr-PRY*I_ESP_Px_S_M1_vrr;
        Double I_ESP_D2y_S_vrr = PAY*I_ESP_Py_S_vrr-PRY*I_ESP_Py_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_D2z_S_vrr = PAZ*I_ESP_Pz_S_vrr-PRZ*I_ESP_Pz_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         ************************************************************/
        Double I_ESP_F3x_S_vrr = PAX*I_ESP_D2x_S_vrr-PRX*I_ESP_D2x_S_M1_vrr+2*oned2z*I_ESP_Px_S_vrr-2*oned2z*I_ESP_Px_S_M1_vrr;
        Double I_ESP_F2xy_S_vrr = PAY*I_ESP_D2x_S_vrr-PRY*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F2xz_S_vrr = PAZ*I_ESP_D2x_S_vrr-PRZ*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Fx2y_S_vrr = PAX*I_ESP_D2y_S_vrr-PRX*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fxyz_S_vrr = PAZ*I_ESP_Dxy_S_vrr-PRZ*I_ESP_Dxy_S_M1_vrr;
        Double I_ESP_Fx2z_S_vrr = PAX*I_ESP_D2z_S_vrr-PRX*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3y_S_vrr = PAY*I_ESP_D2y_S_vrr-PRY*I_ESP_D2y_S_M1_vrr+2*oned2z*I_ESP_Py_S_vrr-2*oned2z*I_ESP_Py_S_M1_vrr;
        Double I_ESP_F2yz_S_vrr = PAZ*I_ESP_D2y_S_vrr-PRZ*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fy2z_S_vrr = PAY*I_ESP_D2z_S_vrr-PRY*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3z_S_vrr = PAZ*I_ESP_D2z_S_vrr-PRZ*I_ESP_D2z_S_M1_vrr+2*oned2z*I_ESP_Pz_S_vrr-2*oned2z*I_ESP_Pz_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         ************************************************************/
        Double I_ESP_G4x_S_vrr = PAX*I_ESP_F3x_S_vrr-PRX*I_ESP_F3x_S_M1_vrr+3*oned2z*I_ESP_D2x_S_vrr-3*oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G3xy_S_vrr = PAY*I_ESP_F3x_S_vrr-PRY*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G3xz_S_vrr = PAZ*I_ESP_F3x_S_vrr-PRZ*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G2x2y_S_vrr = PAY*I_ESP_F2xy_S_vrr-PRY*I_ESP_F2xy_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G2xyz_S_vrr = PAZ*I_ESP_F2xy_S_vrr-PRZ*I_ESP_F2xy_S_M1_vrr;
        Double I_ESP_G2x2z_S_vrr = PAZ*I_ESP_F2xz_S_vrr-PRZ*I_ESP_F2xz_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Gx3y_S_vrr = PAX*I_ESP_F3y_S_vrr-PRX*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_Gx2yz_S_vrr = PAZ*I_ESP_Fx2y_S_vrr-PRZ*I_ESP_Fx2y_S_M1_vrr;
        Double I_ESP_Gxy2z_S_vrr = PAY*I_ESP_Fx2z_S_vrr-PRY*I_ESP_Fx2z_S_M1_vrr;
        Double I_ESP_Gx3z_S_vrr = PAX*I_ESP_F3z_S_vrr-PRX*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_G4y_S_vrr = PAY*I_ESP_F3y_S_vrr-PRY*I_ESP_F3y_S_M1_vrr+3*oned2z*I_ESP_D2y_S_vrr-3*oned2z*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_G3yz_S_vrr = PAZ*I_ESP_F3y_S_vrr-PRZ*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_G2y2z_S_vrr = PAZ*I_ESP_F2yz_S_vrr-PRZ*I_ESP_F2yz_S_M1_vrr+oned2z*I_ESP_D2y_S_vrr-oned2z*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Gy3z_S_vrr = PAY*I_ESP_F3z_S_vrr-PRY*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_G4z_S_vrr = PAZ*I_ESP_F3z_S_vrr-PRZ*I_ESP_F3z_S_M1_vrr+3*oned2z*I_ESP_D2z_S_vrr-3*oned2z*I_ESP_D2z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         ************************************************************/
        Double I_ESP_H5x_S_vrr = PAX*I_ESP_G4x_S_vrr-PRX*I_ESP_G4x_S_M1_vrr+4*oned2z*I_ESP_F3x_S_vrr-4*oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H4xy_S_vrr = PAY*I_ESP_G4x_S_vrr-PRY*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_H4xz_S_vrr = PAZ*I_ESP_G4x_S_vrr-PRZ*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_H3x2y_S_vrr = PAY*I_ESP_G3xy_S_vrr-PRY*I_ESP_G3xy_S_M1_vrr+oned2z*I_ESP_F3x_S_vrr-oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H3xyz_S_vrr = PAZ*I_ESP_G3xy_S_vrr-PRZ*I_ESP_G3xy_S_M1_vrr;
        Double I_ESP_H3x2z_S_vrr = PAZ*I_ESP_G3xz_S_vrr-PRZ*I_ESP_G3xz_S_M1_vrr+oned2z*I_ESP_F3x_S_vrr-oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H2x3y_S_vrr = PAX*I_ESP_Gx3y_S_vrr-PRX*I_ESP_Gx3y_S_M1_vrr+oned2z*I_ESP_F3y_S_vrr-oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H2x2yz_S_vrr = PAZ*I_ESP_G2x2y_S_vrr-PRZ*I_ESP_G2x2y_S_M1_vrr;
        Double I_ESP_H2xy2z_S_vrr = PAY*I_ESP_G2x2z_S_vrr-PRY*I_ESP_G2x2z_S_M1_vrr;
        Double I_ESP_H2x3z_S_vrr = PAX*I_ESP_Gx3z_S_vrr-PRX*I_ESP_Gx3z_S_M1_vrr+oned2z*I_ESP_F3z_S_vrr-oned2z*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_Hx4y_S_vrr = PAX*I_ESP_G4y_S_vrr-PRX*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_Hx3yz_S_vrr = PAZ*I_ESP_Gx3y_S_vrr-PRZ*I_ESP_Gx3y_S_M1_vrr;
        Double I_ESP_Hx2y2z_S_vrr = PAX*I_ESP_G2y2z_S_vrr-PRX*I_ESP_G2y2z_S_M1_vrr;
        Double I_ESP_Hxy3z_S_vrr = PAY*I_ESP_Gx3z_S_vrr-PRY*I_ESP_Gx3z_S_M1_vrr;
        Double I_ESP_Hx4z_S_vrr = PAX*I_ESP_G4z_S_vrr-PRX*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_H5y_S_vrr = PAY*I_ESP_G4y_S_vrr-PRY*I_ESP_G4y_S_M1_vrr+4*oned2z*I_ESP_F3y_S_vrr-4*oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H4yz_S_vrr = PAZ*I_ESP_G4y_S_vrr-PRZ*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_H3y2z_S_vrr = PAZ*I_ESP_G3yz_S_vrr-PRZ*I_ESP_G3yz_S_M1_vrr+oned2z*I_ESP_F3y_S_vrr-oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H2y3z_S_vrr = PAY*I_ESP_Gy3z_S_vrr-PRY*I_ESP_Gy3z_S_M1_vrr+oned2z*I_ESP_F3z_S_vrr-oned2z*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_Hy4z_S_vrr = PAY*I_ESP_G4z_S_vrr-PRY*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_H5z_S_vrr = PAZ*I_ESP_G4z_S_vrr-PRZ*I_ESP_G4z_S_M1_vrr+4*oned2z*I_ESP_F3z_S_vrr-4*oned2z*I_ESP_F3z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S
         * RHS shell quartet name: SQ_ESP_H_S_M1
         * RHS shell quartet name: SQ_ESP_G_S
         * RHS shell quartet name: SQ_ESP_G_S_M1
         ************************************************************/
        Double I_ESP_I6x_S_vrr = PAX*I_ESP_H5x_S_vrr-PRX*I_ESP_H5x_S_M1_vrr+5*oned2z*I_ESP_G4x_S_vrr-5*oned2z*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_I5xy_S_vrr = PAY*I_ESP_H5x_S_vrr-PRY*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_I5xz_S_vrr = PAZ*I_ESP_H5x_S_vrr-PRZ*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_I4x2y_S_vrr = PAY*I_ESP_H4xy_S_vrr-PRY*I_ESP_H4xy_S_M1_vrr+oned2z*I_ESP_G4x_S_vrr-oned2z*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_I4xyz_S_vrr = PAZ*I_ESP_H4xy_S_vrr-PRZ*I_ESP_H4xy_S_M1_vrr;
        Double I_ESP_I4x2z_S_vrr = PAZ*I_ESP_H4xz_S_vrr-PRZ*I_ESP_H4xz_S_M1_vrr+oned2z*I_ESP_G4x_S_vrr-oned2z*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_I3x3y_S_vrr = PAY*I_ESP_H3x2y_S_vrr-PRY*I_ESP_H3x2y_S_M1_vrr+2*oned2z*I_ESP_G3xy_S_vrr-2*oned2z*I_ESP_G3xy_S_M1_vrr;
        Double I_ESP_I3x2yz_S_vrr = PAZ*I_ESP_H3x2y_S_vrr-PRZ*I_ESP_H3x2y_S_M1_vrr;
        Double I_ESP_I3xy2z_S_vrr = PAY*I_ESP_H3x2z_S_vrr-PRY*I_ESP_H3x2z_S_M1_vrr;
        Double I_ESP_I3x3z_S_vrr = PAZ*I_ESP_H3x2z_S_vrr-PRZ*I_ESP_H3x2z_S_M1_vrr+2*oned2z*I_ESP_G3xz_S_vrr-2*oned2z*I_ESP_G3xz_S_M1_vrr;
        Double I_ESP_I2x4y_S_vrr = PAX*I_ESP_Hx4y_S_vrr-PRX*I_ESP_Hx4y_S_M1_vrr+oned2z*I_ESP_G4y_S_vrr-oned2z*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_I2x3yz_S_vrr = PAZ*I_ESP_H2x3y_S_vrr-PRZ*I_ESP_H2x3y_S_M1_vrr;
        Double I_ESP_I2x2y2z_S_vrr = PAZ*I_ESP_H2x2yz_S_vrr-PRZ*I_ESP_H2x2yz_S_M1_vrr+oned2z*I_ESP_G2x2y_S_vrr-oned2z*I_ESP_G2x2y_S_M1_vrr;
        Double I_ESP_I2xy3z_S_vrr = PAY*I_ESP_H2x3z_S_vrr-PRY*I_ESP_H2x3z_S_M1_vrr;
        Double I_ESP_I2x4z_S_vrr = PAX*I_ESP_Hx4z_S_vrr-PRX*I_ESP_Hx4z_S_M1_vrr+oned2z*I_ESP_G4z_S_vrr-oned2z*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_Ix5y_S_vrr = PAX*I_ESP_H5y_S_vrr-PRX*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_Ix4yz_S_vrr = PAZ*I_ESP_Hx4y_S_vrr-PRZ*I_ESP_Hx4y_S_M1_vrr;
        Double I_ESP_Ix3y2z_S_vrr = PAX*I_ESP_H3y2z_S_vrr-PRX*I_ESP_H3y2z_S_M1_vrr;
        Double I_ESP_Ix2y3z_S_vrr = PAX*I_ESP_H2y3z_S_vrr-PRX*I_ESP_H2y3z_S_M1_vrr;
        Double I_ESP_Ixy4z_S_vrr = PAY*I_ESP_Hx4z_S_vrr-PRY*I_ESP_Hx4z_S_M1_vrr;
        Double I_ESP_Ix5z_S_vrr = PAX*I_ESP_H5z_S_vrr-PRX*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_I6y_S_vrr = PAY*I_ESP_H5y_S_vrr-PRY*I_ESP_H5y_S_M1_vrr+5*oned2z*I_ESP_G4y_S_vrr-5*oned2z*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_I5yz_S_vrr = PAZ*I_ESP_H5y_S_vrr-PRZ*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_I4y2z_S_vrr = PAZ*I_ESP_H4yz_S_vrr-PRZ*I_ESP_H4yz_S_M1_vrr+oned2z*I_ESP_G4y_S_vrr-oned2z*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_I3y3z_S_vrr = PAZ*I_ESP_H3y2z_S_vrr-PRZ*I_ESP_H3y2z_S_M1_vrr+2*oned2z*I_ESP_G3yz_S_vrr-2*oned2z*I_ESP_G3yz_S_M1_vrr;
        Double I_ESP_I2y4z_S_vrr = PAY*I_ESP_Hy4z_S_vrr-PRY*I_ESP_Hy4z_S_M1_vrr+oned2z*I_ESP_G4z_S_vrr-oned2z*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_Iy5z_S_vrr = PAY*I_ESP_H5z_S_vrr-PRY*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_I6z_S_vrr = PAZ*I_ESP_H5z_S_vrr-PRZ*I_ESP_H5z_S_M1_vrr+5*oned2z*I_ESP_G4z_S_vrr-5*oned2z*I_ESP_G4z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S
         * RHS shell quartet name: SQ_ESP_I_S_M1
         * RHS shell quartet name: SQ_ESP_H_S
         * RHS shell quartet name: SQ_ESP_H_S_M1
         ************************************************************/
        Double I_ESP_K7x_S_vrr = PAX*I_ESP_I6x_S_vrr-PRX*I_ESP_I6x_S_M1_vrr+6*oned2z*I_ESP_H5x_S_vrr-6*oned2z*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_K6xy_S_vrr = PAY*I_ESP_I6x_S_vrr-PRY*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_K6xz_S_vrr = PAZ*I_ESP_I6x_S_vrr-PRZ*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_K5x2y_S_vrr = PAY*I_ESP_I5xy_S_vrr-PRY*I_ESP_I5xy_S_M1_vrr+oned2z*I_ESP_H5x_S_vrr-oned2z*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_K5xyz_S_vrr = PAZ*I_ESP_I5xy_S_vrr-PRZ*I_ESP_I5xy_S_M1_vrr;
        Double I_ESP_K5x2z_S_vrr = PAZ*I_ESP_I5xz_S_vrr-PRZ*I_ESP_I5xz_S_M1_vrr+oned2z*I_ESP_H5x_S_vrr-oned2z*I_ESP_H5x_S_M1_vrr;
        Double I_ESP_K4x3y_S_vrr = PAY*I_ESP_I4x2y_S_vrr-PRY*I_ESP_I4x2y_S_M1_vrr+2*oned2z*I_ESP_H4xy_S_vrr-2*oned2z*I_ESP_H4xy_S_M1_vrr;
        Double I_ESP_K4x2yz_S_vrr = PAZ*I_ESP_I4x2y_S_vrr-PRZ*I_ESP_I4x2y_S_M1_vrr;
        Double I_ESP_K4xy2z_S_vrr = PAY*I_ESP_I4x2z_S_vrr-PRY*I_ESP_I4x2z_S_M1_vrr;
        Double I_ESP_K4x3z_S_vrr = PAZ*I_ESP_I4x2z_S_vrr-PRZ*I_ESP_I4x2z_S_M1_vrr+2*oned2z*I_ESP_H4xz_S_vrr-2*oned2z*I_ESP_H4xz_S_M1_vrr;
        Double I_ESP_K3x4y_S_vrr = PAX*I_ESP_I2x4y_S_vrr-PRX*I_ESP_I2x4y_S_M1_vrr+2*oned2z*I_ESP_Hx4y_S_vrr-2*oned2z*I_ESP_Hx4y_S_M1_vrr;
        Double I_ESP_K3x3yz_S_vrr = PAZ*I_ESP_I3x3y_S_vrr-PRZ*I_ESP_I3x3y_S_M1_vrr;
        Double I_ESP_K3x2y2z_S_vrr = PAZ*I_ESP_I3x2yz_S_vrr-PRZ*I_ESP_I3x2yz_S_M1_vrr+oned2z*I_ESP_H3x2y_S_vrr-oned2z*I_ESP_H3x2y_S_M1_vrr;
        Double I_ESP_K3xy3z_S_vrr = PAY*I_ESP_I3x3z_S_vrr-PRY*I_ESP_I3x3z_S_M1_vrr;
        Double I_ESP_K3x4z_S_vrr = PAX*I_ESP_I2x4z_S_vrr-PRX*I_ESP_I2x4z_S_M1_vrr+2*oned2z*I_ESP_Hx4z_S_vrr-2*oned2z*I_ESP_Hx4z_S_M1_vrr;
        Double I_ESP_K2x5y_S_vrr = PAX*I_ESP_Ix5y_S_vrr-PRX*I_ESP_Ix5y_S_M1_vrr+oned2z*I_ESP_H5y_S_vrr-oned2z*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_K2x4yz_S_vrr = PAZ*I_ESP_I2x4y_S_vrr-PRZ*I_ESP_I2x4y_S_M1_vrr;
        Double I_ESP_K2x3y2z_S_vrr = PAZ*I_ESP_I2x3yz_S_vrr-PRZ*I_ESP_I2x3yz_S_M1_vrr+oned2z*I_ESP_H2x3y_S_vrr-oned2z*I_ESP_H2x3y_S_M1_vrr;
        Double I_ESP_K2x2y3z_S_vrr = PAY*I_ESP_I2xy3z_S_vrr-PRY*I_ESP_I2xy3z_S_M1_vrr+oned2z*I_ESP_H2x3z_S_vrr-oned2z*I_ESP_H2x3z_S_M1_vrr;
        Double I_ESP_K2xy4z_S_vrr = PAY*I_ESP_I2x4z_S_vrr-PRY*I_ESP_I2x4z_S_M1_vrr;
        Double I_ESP_K2x5z_S_vrr = PAX*I_ESP_Ix5z_S_vrr-PRX*I_ESP_Ix5z_S_M1_vrr+oned2z*I_ESP_H5z_S_vrr-oned2z*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_Kx6y_S_vrr = PAX*I_ESP_I6y_S_vrr-PRX*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_Kx5yz_S_vrr = PAZ*I_ESP_Ix5y_S_vrr-PRZ*I_ESP_Ix5y_S_M1_vrr;
        Double I_ESP_Kx4y2z_S_vrr = PAX*I_ESP_I4y2z_S_vrr-PRX*I_ESP_I4y2z_S_M1_vrr;
        Double I_ESP_Kx3y3z_S_vrr = PAX*I_ESP_I3y3z_S_vrr-PRX*I_ESP_I3y3z_S_M1_vrr;
        Double I_ESP_Kx2y4z_S_vrr = PAX*I_ESP_I2y4z_S_vrr-PRX*I_ESP_I2y4z_S_M1_vrr;
        Double I_ESP_Kxy5z_S_vrr = PAY*I_ESP_Ix5z_S_vrr-PRY*I_ESP_Ix5z_S_M1_vrr;
        Double I_ESP_Kx6z_S_vrr = PAX*I_ESP_I6z_S_vrr-PRX*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_K7y_S_vrr = PAY*I_ESP_I6y_S_vrr-PRY*I_ESP_I6y_S_M1_vrr+6*oned2z*I_ESP_H5y_S_vrr-6*oned2z*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_K6yz_S_vrr = PAZ*I_ESP_I6y_S_vrr-PRZ*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_K5y2z_S_vrr = PAZ*I_ESP_I5yz_S_vrr-PRZ*I_ESP_I5yz_S_M1_vrr+oned2z*I_ESP_H5y_S_vrr-oned2z*I_ESP_H5y_S_M1_vrr;
        Double I_ESP_K4y3z_S_vrr = PAZ*I_ESP_I4y2z_S_vrr-PRZ*I_ESP_I4y2z_S_M1_vrr+2*oned2z*I_ESP_H4yz_S_vrr-2*oned2z*I_ESP_H4yz_S_M1_vrr;
        Double I_ESP_K3y4z_S_vrr = PAY*I_ESP_I2y4z_S_vrr-PRY*I_ESP_I2y4z_S_M1_vrr+2*oned2z*I_ESP_Hy4z_S_vrr-2*oned2z*I_ESP_Hy4z_S_M1_vrr;
        Double I_ESP_K2y5z_S_vrr = PAY*I_ESP_Iy5z_S_vrr-PRY*I_ESP_Iy5z_S_M1_vrr+oned2z*I_ESP_H5z_S_vrr-oned2z*I_ESP_H5z_S_M1_vrr;
        Double I_ESP_Ky6z_S_vrr = PAY*I_ESP_I6z_S_vrr-PRY*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_K7z_S_vrr = PAZ*I_ESP_I6z_S_vrr-PRZ*I_ESP_I6z_S_M1_vrr+6*oned2z*I_ESP_H5z_S_vrr-6*oned2z*I_ESP_H5z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_L_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_K_S
         * RHS shell quartet name: SQ_ESP_K_S_M1
         * RHS shell quartet name: SQ_ESP_I_S
         * RHS shell quartet name: SQ_ESP_I_S_M1
         ************************************************************/
        Double I_ESP_L8x_S_vrr = PAX*I_ESP_K7x_S_vrr-PRX*I_ESP_K7x_S_M1_vrr+7*oned2z*I_ESP_I6x_S_vrr-7*oned2z*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_L7xy_S_vrr = PAY*I_ESP_K7x_S_vrr-PRY*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_L7xz_S_vrr = PAZ*I_ESP_K7x_S_vrr-PRZ*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_L6x2y_S_vrr = PAY*I_ESP_K6xy_S_vrr-PRY*I_ESP_K6xy_S_M1_vrr+oned2z*I_ESP_I6x_S_vrr-oned2z*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_L6xyz_S_vrr = PAZ*I_ESP_K6xy_S_vrr-PRZ*I_ESP_K6xy_S_M1_vrr;
        Double I_ESP_L6x2z_S_vrr = PAZ*I_ESP_K6xz_S_vrr-PRZ*I_ESP_K6xz_S_M1_vrr+oned2z*I_ESP_I6x_S_vrr-oned2z*I_ESP_I6x_S_M1_vrr;
        Double I_ESP_L5x3y_S_vrr = PAY*I_ESP_K5x2y_S_vrr-PRY*I_ESP_K5x2y_S_M1_vrr+2*oned2z*I_ESP_I5xy_S_vrr-2*oned2z*I_ESP_I5xy_S_M1_vrr;
        Double I_ESP_L5x2yz_S_vrr = PAZ*I_ESP_K5x2y_S_vrr-PRZ*I_ESP_K5x2y_S_M1_vrr;
        Double I_ESP_L5xy2z_S_vrr = PAY*I_ESP_K5x2z_S_vrr-PRY*I_ESP_K5x2z_S_M1_vrr;
        Double I_ESP_L5x3z_S_vrr = PAZ*I_ESP_K5x2z_S_vrr-PRZ*I_ESP_K5x2z_S_M1_vrr+2*oned2z*I_ESP_I5xz_S_vrr-2*oned2z*I_ESP_I5xz_S_M1_vrr;
        Double I_ESP_L4x4y_S_vrr = PAY*I_ESP_K4x3y_S_vrr-PRY*I_ESP_K4x3y_S_M1_vrr+3*oned2z*I_ESP_I4x2y_S_vrr-3*oned2z*I_ESP_I4x2y_S_M1_vrr;
        Double I_ESP_L4x3yz_S_vrr = PAZ*I_ESP_K4x3y_S_vrr-PRZ*I_ESP_K4x3y_S_M1_vrr;
        Double I_ESP_L4x2y2z_S_vrr = PAZ*I_ESP_K4x2yz_S_vrr-PRZ*I_ESP_K4x2yz_S_M1_vrr+oned2z*I_ESP_I4x2y_S_vrr-oned2z*I_ESP_I4x2y_S_M1_vrr;
        Double I_ESP_L4xy3z_S_vrr = PAY*I_ESP_K4x3z_S_vrr-PRY*I_ESP_K4x3z_S_M1_vrr;
        Double I_ESP_L4x4z_S_vrr = PAZ*I_ESP_K4x3z_S_vrr-PRZ*I_ESP_K4x3z_S_M1_vrr+3*oned2z*I_ESP_I4x2z_S_vrr-3*oned2z*I_ESP_I4x2z_S_M1_vrr;
        Double I_ESP_L3x5y_S_vrr = PAX*I_ESP_K2x5y_S_vrr-PRX*I_ESP_K2x5y_S_M1_vrr+2*oned2z*I_ESP_Ix5y_S_vrr-2*oned2z*I_ESP_Ix5y_S_M1_vrr;
        Double I_ESP_L3x4yz_S_vrr = PAZ*I_ESP_K3x4y_S_vrr-PRZ*I_ESP_K3x4y_S_M1_vrr;
        Double I_ESP_L3x3y2z_S_vrr = PAZ*I_ESP_K3x3yz_S_vrr-PRZ*I_ESP_K3x3yz_S_M1_vrr+oned2z*I_ESP_I3x3y_S_vrr-oned2z*I_ESP_I3x3y_S_M1_vrr;
        Double I_ESP_L3x2y3z_S_vrr = PAY*I_ESP_K3xy3z_S_vrr-PRY*I_ESP_K3xy3z_S_M1_vrr+oned2z*I_ESP_I3x3z_S_vrr-oned2z*I_ESP_I3x3z_S_M1_vrr;
        Double I_ESP_L3xy4z_S_vrr = PAY*I_ESP_K3x4z_S_vrr-PRY*I_ESP_K3x4z_S_M1_vrr;
        Double I_ESP_L3x5z_S_vrr = PAX*I_ESP_K2x5z_S_vrr-PRX*I_ESP_K2x5z_S_M1_vrr+2*oned2z*I_ESP_Ix5z_S_vrr-2*oned2z*I_ESP_Ix5z_S_M1_vrr;
        Double I_ESP_L2x6y_S_vrr = PAX*I_ESP_Kx6y_S_vrr-PRX*I_ESP_Kx6y_S_M1_vrr+oned2z*I_ESP_I6y_S_vrr-oned2z*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_L2x5yz_S_vrr = PAZ*I_ESP_K2x5y_S_vrr-PRZ*I_ESP_K2x5y_S_M1_vrr;
        Double I_ESP_L2x4y2z_S_vrr = PAZ*I_ESP_K2x4yz_S_vrr-PRZ*I_ESP_K2x4yz_S_M1_vrr+oned2z*I_ESP_I2x4y_S_vrr-oned2z*I_ESP_I2x4y_S_M1_vrr;
        Double I_ESP_L2x3y3z_S_vrr = PAX*I_ESP_Kx3y3z_S_vrr-PRX*I_ESP_Kx3y3z_S_M1_vrr+oned2z*I_ESP_I3y3z_S_vrr-oned2z*I_ESP_I3y3z_S_M1_vrr;
        Double I_ESP_L2x2y4z_S_vrr = PAY*I_ESP_K2xy4z_S_vrr-PRY*I_ESP_K2xy4z_S_M1_vrr+oned2z*I_ESP_I2x4z_S_vrr-oned2z*I_ESP_I2x4z_S_M1_vrr;
        Double I_ESP_L2xy5z_S_vrr = PAY*I_ESP_K2x5z_S_vrr-PRY*I_ESP_K2x5z_S_M1_vrr;
        Double I_ESP_L2x6z_S_vrr = PAX*I_ESP_Kx6z_S_vrr-PRX*I_ESP_Kx6z_S_M1_vrr+oned2z*I_ESP_I6z_S_vrr-oned2z*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_Lx7y_S_vrr = PAX*I_ESP_K7y_S_vrr-PRX*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_Lx6yz_S_vrr = PAZ*I_ESP_Kx6y_S_vrr-PRZ*I_ESP_Kx6y_S_M1_vrr;
        Double I_ESP_Lx5y2z_S_vrr = PAX*I_ESP_K5y2z_S_vrr-PRX*I_ESP_K5y2z_S_M1_vrr;
        Double I_ESP_Lx4y3z_S_vrr = PAX*I_ESP_K4y3z_S_vrr-PRX*I_ESP_K4y3z_S_M1_vrr;
        Double I_ESP_Lx3y4z_S_vrr = PAX*I_ESP_K3y4z_S_vrr-PRX*I_ESP_K3y4z_S_M1_vrr;
        Double I_ESP_Lx2y5z_S_vrr = PAX*I_ESP_K2y5z_S_vrr-PRX*I_ESP_K2y5z_S_M1_vrr;
        Double I_ESP_Lxy6z_S_vrr = PAY*I_ESP_Kx6z_S_vrr-PRY*I_ESP_Kx6z_S_M1_vrr;
        Double I_ESP_Lx7z_S_vrr = PAX*I_ESP_K7z_S_vrr-PRX*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_L8y_S_vrr = PAY*I_ESP_K7y_S_vrr-PRY*I_ESP_K7y_S_M1_vrr+7*oned2z*I_ESP_I6y_S_vrr-7*oned2z*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_L7yz_S_vrr = PAZ*I_ESP_K7y_S_vrr-PRZ*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_L6y2z_S_vrr = PAZ*I_ESP_K6yz_S_vrr-PRZ*I_ESP_K6yz_S_M1_vrr+oned2z*I_ESP_I6y_S_vrr-oned2z*I_ESP_I6y_S_M1_vrr;
        Double I_ESP_L5y3z_S_vrr = PAZ*I_ESP_K5y2z_S_vrr-PRZ*I_ESP_K5y2z_S_M1_vrr+2*oned2z*I_ESP_I5yz_S_vrr-2*oned2z*I_ESP_I5yz_S_M1_vrr;
        Double I_ESP_L4y4z_S_vrr = PAZ*I_ESP_K4y3z_S_vrr-PRZ*I_ESP_K4y3z_S_M1_vrr+3*oned2z*I_ESP_I4y2z_S_vrr-3*oned2z*I_ESP_I4y2z_S_M1_vrr;
        Double I_ESP_L3y5z_S_vrr = PAY*I_ESP_K2y5z_S_vrr-PRY*I_ESP_K2y5z_S_M1_vrr+2*oned2z*I_ESP_Iy5z_S_vrr-2*oned2z*I_ESP_Iy5z_S_M1_vrr;
        Double I_ESP_L2y6z_S_vrr = PAY*I_ESP_Ky6z_S_vrr-PRY*I_ESP_Ky6z_S_M1_vrr+oned2z*I_ESP_I6z_S_vrr-oned2z*I_ESP_I6z_S_M1_vrr;
        Double I_ESP_Ly7z_S_vrr = PAY*I_ESP_K7z_S_vrr-PRY*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_L8z_S_vrr = PAZ*I_ESP_K7z_S_vrr-PRZ*I_ESP_K7z_S_M1_vrr+7*oned2z*I_ESP_I6z_S_vrr-7*oned2z*I_ESP_I6z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_M_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_L_S
         * RHS shell quartet name: SQ_ESP_L_S_M1
         * RHS shell quartet name: SQ_ESP_K_S
         * RHS shell quartet name: SQ_ESP_K_S_M1
         ************************************************************/
        Double I_ESP_M9x_S_vrr = PAX*I_ESP_L8x_S_vrr-PRX*I_ESP_L8x_S_M1_vrr+8*oned2z*I_ESP_K7x_S_vrr-8*oned2z*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_M8xy_S_vrr = PAY*I_ESP_L8x_S_vrr-PRY*I_ESP_L8x_S_M1_vrr;
        Double I_ESP_M8xz_S_vrr = PAZ*I_ESP_L8x_S_vrr-PRZ*I_ESP_L8x_S_M1_vrr;
        Double I_ESP_M7x2y_S_vrr = PAY*I_ESP_L7xy_S_vrr-PRY*I_ESP_L7xy_S_M1_vrr+oned2z*I_ESP_K7x_S_vrr-oned2z*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_M7xyz_S_vrr = PAZ*I_ESP_L7xy_S_vrr-PRZ*I_ESP_L7xy_S_M1_vrr;
        Double I_ESP_M7x2z_S_vrr = PAZ*I_ESP_L7xz_S_vrr-PRZ*I_ESP_L7xz_S_M1_vrr+oned2z*I_ESP_K7x_S_vrr-oned2z*I_ESP_K7x_S_M1_vrr;
        Double I_ESP_M6x3y_S_vrr = PAY*I_ESP_L6x2y_S_vrr-PRY*I_ESP_L6x2y_S_M1_vrr+2*oned2z*I_ESP_K6xy_S_vrr-2*oned2z*I_ESP_K6xy_S_M1_vrr;
        Double I_ESP_M6x2yz_S_vrr = PAZ*I_ESP_L6x2y_S_vrr-PRZ*I_ESP_L6x2y_S_M1_vrr;
        Double I_ESP_M6xy2z_S_vrr = PAY*I_ESP_L6x2z_S_vrr-PRY*I_ESP_L6x2z_S_M1_vrr;
        Double I_ESP_M6x3z_S_vrr = PAZ*I_ESP_L6x2z_S_vrr-PRZ*I_ESP_L6x2z_S_M1_vrr+2*oned2z*I_ESP_K6xz_S_vrr-2*oned2z*I_ESP_K6xz_S_M1_vrr;
        Double I_ESP_M5x4y_S_vrr = PAY*I_ESP_L5x3y_S_vrr-PRY*I_ESP_L5x3y_S_M1_vrr+3*oned2z*I_ESP_K5x2y_S_vrr-3*oned2z*I_ESP_K5x2y_S_M1_vrr;
        Double I_ESP_M5x3yz_S_vrr = PAZ*I_ESP_L5x3y_S_vrr-PRZ*I_ESP_L5x3y_S_M1_vrr;
        Double I_ESP_M5x2y2z_S_vrr = PAZ*I_ESP_L5x2yz_S_vrr-PRZ*I_ESP_L5x2yz_S_M1_vrr+oned2z*I_ESP_K5x2y_S_vrr-oned2z*I_ESP_K5x2y_S_M1_vrr;
        Double I_ESP_M5xy3z_S_vrr = PAY*I_ESP_L5x3z_S_vrr-PRY*I_ESP_L5x3z_S_M1_vrr;
        Double I_ESP_M5x4z_S_vrr = PAZ*I_ESP_L5x3z_S_vrr-PRZ*I_ESP_L5x3z_S_M1_vrr+3*oned2z*I_ESP_K5x2z_S_vrr-3*oned2z*I_ESP_K5x2z_S_M1_vrr;
        Double I_ESP_M4x5y_S_vrr = PAX*I_ESP_L3x5y_S_vrr-PRX*I_ESP_L3x5y_S_M1_vrr+3*oned2z*I_ESP_K2x5y_S_vrr-3*oned2z*I_ESP_K2x5y_S_M1_vrr;
        Double I_ESP_M4x4yz_S_vrr = PAZ*I_ESP_L4x4y_S_vrr-PRZ*I_ESP_L4x4y_S_M1_vrr;
        Double I_ESP_M4x3y2z_S_vrr = PAZ*I_ESP_L4x3yz_S_vrr-PRZ*I_ESP_L4x3yz_S_M1_vrr+oned2z*I_ESP_K4x3y_S_vrr-oned2z*I_ESP_K4x3y_S_M1_vrr;
        Double I_ESP_M4x2y3z_S_vrr = PAY*I_ESP_L4xy3z_S_vrr-PRY*I_ESP_L4xy3z_S_M1_vrr+oned2z*I_ESP_K4x3z_S_vrr-oned2z*I_ESP_K4x3z_S_M1_vrr;
        Double I_ESP_M4xy4z_S_vrr = PAY*I_ESP_L4x4z_S_vrr-PRY*I_ESP_L4x4z_S_M1_vrr;
        Double I_ESP_M4x5z_S_vrr = PAX*I_ESP_L3x5z_S_vrr-PRX*I_ESP_L3x5z_S_M1_vrr+3*oned2z*I_ESP_K2x5z_S_vrr-3*oned2z*I_ESP_K2x5z_S_M1_vrr;
        Double I_ESP_M3x6y_S_vrr = PAX*I_ESP_L2x6y_S_vrr-PRX*I_ESP_L2x6y_S_M1_vrr+2*oned2z*I_ESP_Kx6y_S_vrr-2*oned2z*I_ESP_Kx6y_S_M1_vrr;
        Double I_ESP_M3x5yz_S_vrr = PAZ*I_ESP_L3x5y_S_vrr-PRZ*I_ESP_L3x5y_S_M1_vrr;
        Double I_ESP_M3x4y2z_S_vrr = PAZ*I_ESP_L3x4yz_S_vrr-PRZ*I_ESP_L3x4yz_S_M1_vrr+oned2z*I_ESP_K3x4y_S_vrr-oned2z*I_ESP_K3x4y_S_M1_vrr;
        Double I_ESP_M3x3y3z_S_vrr = PAZ*I_ESP_L3x3y2z_S_vrr-PRZ*I_ESP_L3x3y2z_S_M1_vrr+2*oned2z*I_ESP_K3x3yz_S_vrr-2*oned2z*I_ESP_K3x3yz_S_M1_vrr;
        Double I_ESP_M3x2y4z_S_vrr = PAY*I_ESP_L3xy4z_S_vrr-PRY*I_ESP_L3xy4z_S_M1_vrr+oned2z*I_ESP_K3x4z_S_vrr-oned2z*I_ESP_K3x4z_S_M1_vrr;
        Double I_ESP_M3xy5z_S_vrr = PAY*I_ESP_L3x5z_S_vrr-PRY*I_ESP_L3x5z_S_M1_vrr;
        Double I_ESP_M3x6z_S_vrr = PAX*I_ESP_L2x6z_S_vrr-PRX*I_ESP_L2x6z_S_M1_vrr+2*oned2z*I_ESP_Kx6z_S_vrr-2*oned2z*I_ESP_Kx6z_S_M1_vrr;
        Double I_ESP_M2x7y_S_vrr = PAX*I_ESP_Lx7y_S_vrr-PRX*I_ESP_Lx7y_S_M1_vrr+oned2z*I_ESP_K7y_S_vrr-oned2z*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_M2x6yz_S_vrr = PAZ*I_ESP_L2x6y_S_vrr-PRZ*I_ESP_L2x6y_S_M1_vrr;
        Double I_ESP_M2x5y2z_S_vrr = PAZ*I_ESP_L2x5yz_S_vrr-PRZ*I_ESP_L2x5yz_S_M1_vrr+oned2z*I_ESP_K2x5y_S_vrr-oned2z*I_ESP_K2x5y_S_M1_vrr;
        Double I_ESP_M2x4y3z_S_vrr = PAX*I_ESP_Lx4y3z_S_vrr-PRX*I_ESP_Lx4y3z_S_M1_vrr+oned2z*I_ESP_K4y3z_S_vrr-oned2z*I_ESP_K4y3z_S_M1_vrr;
        Double I_ESP_M2x3y4z_S_vrr = PAX*I_ESP_Lx3y4z_S_vrr-PRX*I_ESP_Lx3y4z_S_M1_vrr+oned2z*I_ESP_K3y4z_S_vrr-oned2z*I_ESP_K3y4z_S_M1_vrr;
        Double I_ESP_M2x2y5z_S_vrr = PAY*I_ESP_L2xy5z_S_vrr-PRY*I_ESP_L2xy5z_S_M1_vrr+oned2z*I_ESP_K2x5z_S_vrr-oned2z*I_ESP_K2x5z_S_M1_vrr;
        Double I_ESP_M2xy6z_S_vrr = PAY*I_ESP_L2x6z_S_vrr-PRY*I_ESP_L2x6z_S_M1_vrr;
        Double I_ESP_M2x7z_S_vrr = PAX*I_ESP_Lx7z_S_vrr-PRX*I_ESP_Lx7z_S_M1_vrr+oned2z*I_ESP_K7z_S_vrr-oned2z*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_Mx8y_S_vrr = PAX*I_ESP_L8y_S_vrr-PRX*I_ESP_L8y_S_M1_vrr;
        Double I_ESP_Mx7yz_S_vrr = PAZ*I_ESP_Lx7y_S_vrr-PRZ*I_ESP_Lx7y_S_M1_vrr;
        Double I_ESP_Mx6y2z_S_vrr = PAX*I_ESP_L6y2z_S_vrr-PRX*I_ESP_L6y2z_S_M1_vrr;
        Double I_ESP_Mx5y3z_S_vrr = PAX*I_ESP_L5y3z_S_vrr-PRX*I_ESP_L5y3z_S_M1_vrr;
        Double I_ESP_Mx4y4z_S_vrr = PAX*I_ESP_L4y4z_S_vrr-PRX*I_ESP_L4y4z_S_M1_vrr;
        Double I_ESP_Mx3y5z_S_vrr = PAX*I_ESP_L3y5z_S_vrr-PRX*I_ESP_L3y5z_S_M1_vrr;
        Double I_ESP_Mx2y6z_S_vrr = PAX*I_ESP_L2y6z_S_vrr-PRX*I_ESP_L2y6z_S_M1_vrr;
        Double I_ESP_Mxy7z_S_vrr = PAY*I_ESP_Lx7z_S_vrr-PRY*I_ESP_Lx7z_S_M1_vrr;
        Double I_ESP_Mx8z_S_vrr = PAX*I_ESP_L8z_S_vrr-PRX*I_ESP_L8z_S_M1_vrr;
        Double I_ESP_M9y_S_vrr = PAY*I_ESP_L8y_S_vrr-PRY*I_ESP_L8y_S_M1_vrr+8*oned2z*I_ESP_K7y_S_vrr-8*oned2z*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_M8yz_S_vrr = PAZ*I_ESP_L8y_S_vrr-PRZ*I_ESP_L8y_S_M1_vrr;
        Double I_ESP_M7y2z_S_vrr = PAZ*I_ESP_L7yz_S_vrr-PRZ*I_ESP_L7yz_S_M1_vrr+oned2z*I_ESP_K7y_S_vrr-oned2z*I_ESP_K7y_S_M1_vrr;
        Double I_ESP_M6y3z_S_vrr = PAZ*I_ESP_L6y2z_S_vrr-PRZ*I_ESP_L6y2z_S_M1_vrr+2*oned2z*I_ESP_K6yz_S_vrr-2*oned2z*I_ESP_K6yz_S_M1_vrr;
        Double I_ESP_M5y4z_S_vrr = PAZ*I_ESP_L5y3z_S_vrr-PRZ*I_ESP_L5y3z_S_M1_vrr+3*oned2z*I_ESP_K5y2z_S_vrr-3*oned2z*I_ESP_K5y2z_S_M1_vrr;
        Double I_ESP_M4y5z_S_vrr = PAY*I_ESP_L3y5z_S_vrr-PRY*I_ESP_L3y5z_S_M1_vrr+3*oned2z*I_ESP_K2y5z_S_vrr-3*oned2z*I_ESP_K2y5z_S_M1_vrr;
        Double I_ESP_M3y6z_S_vrr = PAY*I_ESP_L2y6z_S_vrr-PRY*I_ESP_L2y6z_S_M1_vrr+2*oned2z*I_ESP_Ky6z_S_vrr-2*oned2z*I_ESP_Ky6z_S_M1_vrr;
        Double I_ESP_M2y7z_S_vrr = PAY*I_ESP_Ly7z_S_vrr-PRY*I_ESP_Ly7z_S_M1_vrr+oned2z*I_ESP_K7z_S_vrr-oned2z*I_ESP_K7z_S_M1_vrr;
        Double I_ESP_My8z_S_vrr = PAY*I_ESP_L8z_S_vrr-PRY*I_ESP_L8z_S_M1_vrr;
        Double I_ESP_M9z_S_vrr = PAZ*I_ESP_L8z_S_vrr-PRZ*I_ESP_L8z_S_M1_vrr+8*oned2z*I_ESP_K7z_S_vrr-8*oned2z*I_ESP_K7z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_M_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_M_S_aa_coefs = alpha*alpha;
        I_ESP_M9x_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M9x_S_vrr;
        I_ESP_M8xy_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M8xy_S_vrr;
        I_ESP_M8xz_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M8xz_S_vrr;
        I_ESP_M7x2y_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M7x2y_S_vrr;
        I_ESP_M7xyz_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M7xyz_S_vrr;
        I_ESP_M7x2z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M7x2z_S_vrr;
        I_ESP_M6x3y_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M6x3y_S_vrr;
        I_ESP_M6x2yz_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M6x2yz_S_vrr;
        I_ESP_M6xy2z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M6xy2z_S_vrr;
        I_ESP_M6x3z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M6x3z_S_vrr;
        I_ESP_M5x4y_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M5x4y_S_vrr;
        I_ESP_M5x3yz_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M5x3yz_S_vrr;
        I_ESP_M5x2y2z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M5x2y2z_S_vrr;
        I_ESP_M5xy3z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M5xy3z_S_vrr;
        I_ESP_M5x4z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M5x4z_S_vrr;
        I_ESP_M4x5y_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M4x5y_S_vrr;
        I_ESP_M4x4yz_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M4x4yz_S_vrr;
        I_ESP_M4x3y2z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M4x3y2z_S_vrr;
        I_ESP_M4x2y3z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M4x2y3z_S_vrr;
        I_ESP_M4xy4z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M4xy4z_S_vrr;
        I_ESP_M4x5z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M4x5z_S_vrr;
        I_ESP_M3x6y_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M3x6y_S_vrr;
        I_ESP_M3x5yz_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M3x5yz_S_vrr;
        I_ESP_M3x4y2z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M3x4y2z_S_vrr;
        I_ESP_M3x3y3z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M3x3y3z_S_vrr;
        I_ESP_M3x2y4z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M3x2y4z_S_vrr;
        I_ESP_M3xy5z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M3xy5z_S_vrr;
        I_ESP_M3x6z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M3x6z_S_vrr;
        I_ESP_M2x7y_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M2x7y_S_vrr;
        I_ESP_M2x6yz_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M2x6yz_S_vrr;
        I_ESP_M2x5y2z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M2x5y2z_S_vrr;
        I_ESP_M2x4y3z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M2x4y3z_S_vrr;
        I_ESP_M2x3y4z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M2x3y4z_S_vrr;
        I_ESP_M2x2y5z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M2x2y5z_S_vrr;
        I_ESP_M2xy6z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M2xy6z_S_vrr;
        I_ESP_M2x7z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M2x7z_S_vrr;
        I_ESP_Mx8y_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_Mx8y_S_vrr;
        I_ESP_Mx7yz_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_Mx7yz_S_vrr;
        I_ESP_Mx6y2z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_Mx6y2z_S_vrr;
        I_ESP_Mx5y3z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_Mx5y3z_S_vrr;
        I_ESP_Mx4y4z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_Mx4y4z_S_vrr;
        I_ESP_Mx3y5z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_Mx3y5z_S_vrr;
        I_ESP_Mx2y6z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_Mx2y6z_S_vrr;
        I_ESP_Mxy7z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_Mxy7z_S_vrr;
        I_ESP_Mx8z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_Mx8z_S_vrr;
        I_ESP_M9y_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M9y_S_vrr;
        I_ESP_M8yz_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M8yz_S_vrr;
        I_ESP_M7y2z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M7y2z_S_vrr;
        I_ESP_M6y3z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M6y3z_S_vrr;
        I_ESP_M5y4z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M5y4z_S_vrr;
        I_ESP_M4y5z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M4y5z_S_vrr;
        I_ESP_M3y6z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M3y6z_S_vrr;
        I_ESP_M2y7z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M2y7z_S_vrr;
        I_ESP_My8z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_My8z_S_vrr;
        I_ESP_M9z_S_aa += SQ_ESP_M_S_aa_coefs*I_ESP_M9z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_L_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_L_S_aa_coefs = alpha*alpha;
        I_ESP_L8x_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L8x_S_vrr;
        I_ESP_L7xy_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L7xy_S_vrr;
        I_ESP_L7xz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L7xz_S_vrr;
        I_ESP_L6x2y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L6x2y_S_vrr;
        I_ESP_L6xyz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L6xyz_S_vrr;
        I_ESP_L6x2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L6x2z_S_vrr;
        I_ESP_L5x3y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5x3y_S_vrr;
        I_ESP_L5x2yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5x2yz_S_vrr;
        I_ESP_L5xy2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5xy2z_S_vrr;
        I_ESP_L5x3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5x3z_S_vrr;
        I_ESP_L4x4y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4x4y_S_vrr;
        I_ESP_L4x3yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4x3yz_S_vrr;
        I_ESP_L4x2y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4x2y2z_S_vrr;
        I_ESP_L4xy3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4xy3z_S_vrr;
        I_ESP_L4x4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4x4z_S_vrr;
        I_ESP_L3x5y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x5y_S_vrr;
        I_ESP_L3x4yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x4yz_S_vrr;
        I_ESP_L3x3y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x3y2z_S_vrr;
        I_ESP_L3x2y3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x2y3z_S_vrr;
        I_ESP_L3xy4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3xy4z_S_vrr;
        I_ESP_L3x5z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3x5z_S_vrr;
        I_ESP_L2x6y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x6y_S_vrr;
        I_ESP_L2x5yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x5yz_S_vrr;
        I_ESP_L2x4y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x4y2z_S_vrr;
        I_ESP_L2x3y3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x3y3z_S_vrr;
        I_ESP_L2x2y4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x2y4z_S_vrr;
        I_ESP_L2xy5z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2xy5z_S_vrr;
        I_ESP_L2x6z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2x6z_S_vrr;
        I_ESP_Lx7y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx7y_S_vrr;
        I_ESP_Lx6yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx6yz_S_vrr;
        I_ESP_Lx5y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx5y2z_S_vrr;
        I_ESP_Lx4y3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx4y3z_S_vrr;
        I_ESP_Lx3y4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx3y4z_S_vrr;
        I_ESP_Lx2y5z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx2y5z_S_vrr;
        I_ESP_Lxy6z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lxy6z_S_vrr;
        I_ESP_Lx7z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Lx7z_S_vrr;
        I_ESP_L8y_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L8y_S_vrr;
        I_ESP_L7yz_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L7yz_S_vrr;
        I_ESP_L6y2z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L6y2z_S_vrr;
        I_ESP_L5y3z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L5y3z_S_vrr;
        I_ESP_L4y4z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L4y4z_S_vrr;
        I_ESP_L3y5z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L3y5z_S_vrr;
        I_ESP_L2y6z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L2y6z_S_vrr;
        I_ESP_Ly7z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_Ly7z_S_vrr;
        I_ESP_L8z_S_aa += SQ_ESP_L_S_aa_coefs*I_ESP_L8z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_K_S_aa_coefs = alpha*alpha;
        I_ESP_K7x_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K7x_S_vrr;
        I_ESP_K6xy_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K6xy_S_vrr;
        I_ESP_K6xz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K6xz_S_vrr;
        I_ESP_K5x2y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K5x2y_S_vrr;
        I_ESP_K5xyz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K5xyz_S_vrr;
        I_ESP_K5x2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K5x2z_S_vrr;
        I_ESP_K4x3y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4x3y_S_vrr;
        I_ESP_K4x2yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4x2yz_S_vrr;
        I_ESP_K4xy2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4xy2z_S_vrr;
        I_ESP_K4x3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4x3z_S_vrr;
        I_ESP_K3x4y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3x4y_S_vrr;
        I_ESP_K3x3yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3x3yz_S_vrr;
        I_ESP_K3x2y2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3x2y2z_S_vrr;
        I_ESP_K3xy3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3xy3z_S_vrr;
        I_ESP_K3x4z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3x4z_S_vrr;
        I_ESP_K2x5y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x5y_S_vrr;
        I_ESP_K2x4yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x4yz_S_vrr;
        I_ESP_K2x3y2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x3y2z_S_vrr;
        I_ESP_K2x2y3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x2y3z_S_vrr;
        I_ESP_K2xy4z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2xy4z_S_vrr;
        I_ESP_K2x5z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2x5z_S_vrr;
        I_ESP_Kx6y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx6y_S_vrr;
        I_ESP_Kx5yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx5yz_S_vrr;
        I_ESP_Kx4y2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx4y2z_S_vrr;
        I_ESP_Kx3y3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx3y3z_S_vrr;
        I_ESP_Kx2y4z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx2y4z_S_vrr;
        I_ESP_Kxy5z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kxy5z_S_vrr;
        I_ESP_Kx6z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Kx6z_S_vrr;
        I_ESP_K7y_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K7y_S_vrr;
        I_ESP_K6yz_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K6yz_S_vrr;
        I_ESP_K5y2z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K5y2z_S_vrr;
        I_ESP_K4y3z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K4y3z_S_vrr;
        I_ESP_K3y4z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K3y4z_S_vrr;
        I_ESP_K2y5z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K2y5z_S_vrr;
        I_ESP_Ky6z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_Ky6z_S_vrr;
        I_ESP_K7z_S_aa += SQ_ESP_K_S_aa_coefs*I_ESP_K7z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_K_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_K_S_a_coefs = alpha;
        I_ESP_K7x_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K7x_S_vrr;
        I_ESP_K6xy_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K6xy_S_vrr;
        I_ESP_K6xz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K6xz_S_vrr;
        I_ESP_K5x2y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K5x2y_S_vrr;
        I_ESP_K5xyz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K5xyz_S_vrr;
        I_ESP_K5x2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K5x2z_S_vrr;
        I_ESP_K4x3y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4x3y_S_vrr;
        I_ESP_K4x2yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4x2yz_S_vrr;
        I_ESP_K4xy2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4xy2z_S_vrr;
        I_ESP_K4x3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4x3z_S_vrr;
        I_ESP_K3x4y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3x4y_S_vrr;
        I_ESP_K3x3yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3x3yz_S_vrr;
        I_ESP_K3x2y2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3x2y2z_S_vrr;
        I_ESP_K3xy3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3xy3z_S_vrr;
        I_ESP_K3x4z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3x4z_S_vrr;
        I_ESP_K2x5y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x5y_S_vrr;
        I_ESP_K2x4yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x4yz_S_vrr;
        I_ESP_K2x3y2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x3y2z_S_vrr;
        I_ESP_K2x2y3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x2y3z_S_vrr;
        I_ESP_K2xy4z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2xy4z_S_vrr;
        I_ESP_K2x5z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2x5z_S_vrr;
        I_ESP_Kx6y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx6y_S_vrr;
        I_ESP_Kx5yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx5yz_S_vrr;
        I_ESP_Kx4y2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx4y2z_S_vrr;
        I_ESP_Kx3y3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx3y3z_S_vrr;
        I_ESP_Kx2y4z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx2y4z_S_vrr;
        I_ESP_Kxy5z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kxy5z_S_vrr;
        I_ESP_Kx6z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Kx6z_S_vrr;
        I_ESP_K7y_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K7y_S_vrr;
        I_ESP_K6yz_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K6yz_S_vrr;
        I_ESP_K5y2z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K5y2z_S_vrr;
        I_ESP_K4y3z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K4y3z_S_vrr;
        I_ESP_K3y4z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K3y4z_S_vrr;
        I_ESP_K2y5z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K2y5z_S_vrr;
        I_ESP_Ky6z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_Ky6z_S_vrr;
        I_ESP_K7z_S_a += SQ_ESP_K_S_a_coefs*I_ESP_K7z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_I_S_a_coefs = alpha;
        I_ESP_I6x_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S_a += SQ_ESP_I_S_a_coefs*I_ESP_I6z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_H_S_a_coefs = alpha;
        I_ESP_H5x_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H5x_S_vrr;
        I_ESP_H4xy_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H4xy_S_vrr;
        I_ESP_H4xz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H4xz_S_vrr;
        I_ESP_H3x2y_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H3x2y_S_vrr;
        I_ESP_H3xyz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H3xyz_S_vrr;
        I_ESP_H3x2z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H3x2z_S_vrr;
        I_ESP_H2x3y_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2x3y_S_vrr;
        I_ESP_H2x2yz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2x2yz_S_vrr;
        I_ESP_H2xy2z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2xy2z_S_vrr;
        I_ESP_H2x3z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2x3z_S_vrr;
        I_ESP_Hx4y_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hx4y_S_vrr;
        I_ESP_Hx3yz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hx3yz_S_vrr;
        I_ESP_Hx2y2z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hx2y2z_S_vrr;
        I_ESP_Hxy3z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hxy3z_S_vrr;
        I_ESP_Hx4z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hx4z_S_vrr;
        I_ESP_H5y_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H5y_S_vrr;
        I_ESP_H4yz_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H4yz_S_vrr;
        I_ESP_H3y2z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H3y2z_S_vrr;
        I_ESP_H2y3z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H2y3z_S_vrr;
        I_ESP_Hy4z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_Hy4z_S_vrr;
        I_ESP_H5z_S_a += SQ_ESP_H_S_a_coefs*I_ESP_H5z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_H5x_S += I_ESP_H5x_S_vrr;
        I_ESP_H4xy_S += I_ESP_H4xy_S_vrr;
        I_ESP_H4xz_S += I_ESP_H4xz_S_vrr;
        I_ESP_H3x2y_S += I_ESP_H3x2y_S_vrr;
        I_ESP_H3xyz_S += I_ESP_H3xyz_S_vrr;
        I_ESP_H3x2z_S += I_ESP_H3x2z_S_vrr;
        I_ESP_H2x3y_S += I_ESP_H2x3y_S_vrr;
        I_ESP_H2x2yz_S += I_ESP_H2x2yz_S_vrr;
        I_ESP_H2xy2z_S += I_ESP_H2xy2z_S_vrr;
        I_ESP_H2x3z_S += I_ESP_H2x3z_S_vrr;
        I_ESP_Hx4y_S += I_ESP_Hx4y_S_vrr;
        I_ESP_Hx3yz_S += I_ESP_Hx3yz_S_vrr;
        I_ESP_Hx2y2z_S += I_ESP_Hx2y2z_S_vrr;
        I_ESP_Hxy3z_S += I_ESP_Hxy3z_S_vrr;
        I_ESP_Hx4z_S += I_ESP_Hx4z_S_vrr;
        I_ESP_H5y_S += I_ESP_H5y_S_vrr;
        I_ESP_H4yz_S += I_ESP_H4yz_S_vrr;
        I_ESP_H3y2z_S += I_ESP_H3y2z_S_vrr;
        I_ESP_H2y3z_S += I_ESP_H2y3z_S_vrr;
        I_ESP_Hy4z_S += I_ESP_Hy4z_S_vrr;
        I_ESP_H5z_S += I_ESP_H5z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_G4x_S += I_ESP_G4x_S_vrr;
        I_ESP_G3xy_S += I_ESP_G3xy_S_vrr;
        I_ESP_G3xz_S += I_ESP_G3xz_S_vrr;
        I_ESP_G2x2y_S += I_ESP_G2x2y_S_vrr;
        I_ESP_G2xyz_S += I_ESP_G2xyz_S_vrr;
        I_ESP_G2x2z_S += I_ESP_G2x2z_S_vrr;
        I_ESP_Gx3y_S += I_ESP_Gx3y_S_vrr;
        I_ESP_Gx2yz_S += I_ESP_Gx2yz_S_vrr;
        I_ESP_Gxy2z_S += I_ESP_Gxy2z_S_vrr;
        I_ESP_Gx3z_S += I_ESP_Gx3z_S_vrr;
        I_ESP_G4y_S += I_ESP_G4y_S_vrr;
        I_ESP_G3yz_S += I_ESP_G3yz_S_vrr;
        I_ESP_G2y2z_S += I_ESP_G2y2z_S_vrr;
        I_ESP_Gy3z_S += I_ESP_Gy3z_S_vrr;
        I_ESP_G4z_S += I_ESP_G4z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_F3x_S += I_ESP_F3x_S_vrr;
        I_ESP_F2xy_S += I_ESP_F2xy_S_vrr;
        I_ESP_F2xz_S += I_ESP_F2xz_S_vrr;
        I_ESP_Fx2y_S += I_ESP_Fx2y_S_vrr;
        I_ESP_Fxyz_S += I_ESP_Fxyz_S_vrr;
        I_ESP_Fx2z_S += I_ESP_Fx2z_S_vrr;
        I_ESP_F3y_S += I_ESP_F3y_S_vrr;
        I_ESP_F2yz_S += I_ESP_F2yz_S_vrr;
        I_ESP_Fy2z_S += I_ESP_Fy2z_S_vrr;
        I_ESP_F3z_S += I_ESP_F3z_S_vrr;
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
     * shell quartet name: SQ_ESP_F_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_S
     * RHS shell quartet name: SQ_ESP_F_S
     ************************************************************/
    Double I_ESP_F3x_Px = I_ESP_G4x_S+ABX*I_ESP_F3x_S;
    Double I_ESP_F2xy_Px = I_ESP_G3xy_S+ABX*I_ESP_F2xy_S;
    Double I_ESP_F2xz_Px = I_ESP_G3xz_S+ABX*I_ESP_F2xz_S;
    Double I_ESP_Fx2y_Px = I_ESP_G2x2y_S+ABX*I_ESP_Fx2y_S;
    Double I_ESP_Fxyz_Px = I_ESP_G2xyz_S+ABX*I_ESP_Fxyz_S;
    Double I_ESP_Fx2z_Px = I_ESP_G2x2z_S+ABX*I_ESP_Fx2z_S;
    Double I_ESP_F3y_Px = I_ESP_Gx3y_S+ABX*I_ESP_F3y_S;
    Double I_ESP_F2yz_Px = I_ESP_Gx2yz_S+ABX*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Px = I_ESP_Gxy2z_S+ABX*I_ESP_Fy2z_S;
    Double I_ESP_F3z_Px = I_ESP_Gx3z_S+ABX*I_ESP_F3z_S;
    Double I_ESP_F3x_Py = I_ESP_G3xy_S+ABY*I_ESP_F3x_S;
    Double I_ESP_F2xy_Py = I_ESP_G2x2y_S+ABY*I_ESP_F2xy_S;
    Double I_ESP_F2xz_Py = I_ESP_G2xyz_S+ABY*I_ESP_F2xz_S;
    Double I_ESP_Fx2y_Py = I_ESP_Gx3y_S+ABY*I_ESP_Fx2y_S;
    Double I_ESP_Fxyz_Py = I_ESP_Gx2yz_S+ABY*I_ESP_Fxyz_S;
    Double I_ESP_Fx2z_Py = I_ESP_Gxy2z_S+ABY*I_ESP_Fx2z_S;
    Double I_ESP_F3y_Py = I_ESP_G4y_S+ABY*I_ESP_F3y_S;
    Double I_ESP_F2yz_Py = I_ESP_G3yz_S+ABY*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Py = I_ESP_G2y2z_S+ABY*I_ESP_Fy2z_S;
    Double I_ESP_F3z_Py = I_ESP_Gy3z_S+ABY*I_ESP_F3z_S;
    Double I_ESP_F3x_Pz = I_ESP_G3xz_S+ABZ*I_ESP_F3x_S;
    Double I_ESP_F2xy_Pz = I_ESP_G2xyz_S+ABZ*I_ESP_F2xy_S;
    Double I_ESP_F2xz_Pz = I_ESP_G2x2z_S+ABZ*I_ESP_F2xz_S;
    Double I_ESP_Fx2y_Pz = I_ESP_Gx2yz_S+ABZ*I_ESP_Fx2y_S;
    Double I_ESP_Fxyz_Pz = I_ESP_Gxy2z_S+ABZ*I_ESP_Fxyz_S;
    Double I_ESP_Fx2z_Pz = I_ESP_Gx3z_S+ABZ*I_ESP_Fx2z_S;
    Double I_ESP_F3y_Pz = I_ESP_G3yz_S+ABZ*I_ESP_F3y_S;
    Double I_ESP_F2yz_Pz = I_ESP_G2y2z_S+ABZ*I_ESP_F2yz_S;
    Double I_ESP_Fy2z_Pz = I_ESP_Gy3z_S+ABZ*I_ESP_Fy2z_S;
    Double I_ESP_F3z_Pz = I_ESP_G4z_S+ABZ*I_ESP_F3z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S
     * RHS shell quartet name: SQ_ESP_G_S
     ************************************************************/
    Double I_ESP_G4x_Px = I_ESP_H5x_S+ABX*I_ESP_G4x_S;
    Double I_ESP_G3xy_Px = I_ESP_H4xy_S+ABX*I_ESP_G3xy_S;
    Double I_ESP_G3xz_Px = I_ESP_H4xz_S+ABX*I_ESP_G3xz_S;
    Double I_ESP_G2x2y_Px = I_ESP_H3x2y_S+ABX*I_ESP_G2x2y_S;
    Double I_ESP_G2xyz_Px = I_ESP_H3xyz_S+ABX*I_ESP_G2xyz_S;
    Double I_ESP_G2x2z_Px = I_ESP_H3x2z_S+ABX*I_ESP_G2x2z_S;
    Double I_ESP_Gx3y_Px = I_ESP_H2x3y_S+ABX*I_ESP_Gx3y_S;
    Double I_ESP_Gx2yz_Px = I_ESP_H2x2yz_S+ABX*I_ESP_Gx2yz_S;
    Double I_ESP_Gxy2z_Px = I_ESP_H2xy2z_S+ABX*I_ESP_Gxy2z_S;
    Double I_ESP_Gx3z_Px = I_ESP_H2x3z_S+ABX*I_ESP_Gx3z_S;
    Double I_ESP_G4y_Px = I_ESP_Hx4y_S+ABX*I_ESP_G4y_S;
    Double I_ESP_G3yz_Px = I_ESP_Hx3yz_S+ABX*I_ESP_G3yz_S;
    Double I_ESP_G2y2z_Px = I_ESP_Hx2y2z_S+ABX*I_ESP_G2y2z_S;
    Double I_ESP_Gy3z_Px = I_ESP_Hxy3z_S+ABX*I_ESP_Gy3z_S;
    Double I_ESP_G4z_Px = I_ESP_Hx4z_S+ABX*I_ESP_G4z_S;
    Double I_ESP_G3xy_Py = I_ESP_H3x2y_S+ABY*I_ESP_G3xy_S;
    Double I_ESP_G3xz_Py = I_ESP_H3xyz_S+ABY*I_ESP_G3xz_S;
    Double I_ESP_G2x2y_Py = I_ESP_H2x3y_S+ABY*I_ESP_G2x2y_S;
    Double I_ESP_G2xyz_Py = I_ESP_H2x2yz_S+ABY*I_ESP_G2xyz_S;
    Double I_ESP_G2x2z_Py = I_ESP_H2xy2z_S+ABY*I_ESP_G2x2z_S;
    Double I_ESP_Gx3y_Py = I_ESP_Hx4y_S+ABY*I_ESP_Gx3y_S;
    Double I_ESP_Gx2yz_Py = I_ESP_Hx3yz_S+ABY*I_ESP_Gx2yz_S;
    Double I_ESP_Gxy2z_Py = I_ESP_Hx2y2z_S+ABY*I_ESP_Gxy2z_S;
    Double I_ESP_Gx3z_Py = I_ESP_Hxy3z_S+ABY*I_ESP_Gx3z_S;
    Double I_ESP_G4y_Py = I_ESP_H5y_S+ABY*I_ESP_G4y_S;
    Double I_ESP_G3yz_Py = I_ESP_H4yz_S+ABY*I_ESP_G3yz_S;
    Double I_ESP_G2y2z_Py = I_ESP_H3y2z_S+ABY*I_ESP_G2y2z_S;
    Double I_ESP_Gy3z_Py = I_ESP_H2y3z_S+ABY*I_ESP_Gy3z_S;
    Double I_ESP_G4z_Py = I_ESP_Hy4z_S+ABY*I_ESP_G4z_S;
    Double I_ESP_G3xz_Pz = I_ESP_H3x2z_S+ABZ*I_ESP_G3xz_S;
    Double I_ESP_G2xyz_Pz = I_ESP_H2xy2z_S+ABZ*I_ESP_G2xyz_S;
    Double I_ESP_G2x2z_Pz = I_ESP_H2x3z_S+ABZ*I_ESP_G2x2z_S;
    Double I_ESP_Gx2yz_Pz = I_ESP_Hx2y2z_S+ABZ*I_ESP_Gx2yz_S;
    Double I_ESP_Gxy2z_Pz = I_ESP_Hxy3z_S+ABZ*I_ESP_Gxy2z_S;
    Double I_ESP_Gx3z_Pz = I_ESP_Hx4z_S+ABZ*I_ESP_Gx3z_S;
    Double I_ESP_G3yz_Pz = I_ESP_H3y2z_S+ABZ*I_ESP_G3yz_S;
    Double I_ESP_G2y2z_Pz = I_ESP_H2y3z_S+ABZ*I_ESP_G2y2z_S;
    Double I_ESP_Gy3z_Pz = I_ESP_Hy4z_S+ABZ*I_ESP_Gy3z_S;
    Double I_ESP_G4z_Pz = I_ESP_H5z_S+ABZ*I_ESP_G4z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_F_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_P
     * RHS shell quartet name: SQ_ESP_F_P
     ************************************************************/
    Double I_ESP_F3x_D2x = I_ESP_G4x_Px+ABX*I_ESP_F3x_Px;
    Double I_ESP_F2xy_D2x = I_ESP_G3xy_Px+ABX*I_ESP_F2xy_Px;
    Double I_ESP_F2xz_D2x = I_ESP_G3xz_Px+ABX*I_ESP_F2xz_Px;
    Double I_ESP_Fx2y_D2x = I_ESP_G2x2y_Px+ABX*I_ESP_Fx2y_Px;
    Double I_ESP_Fxyz_D2x = I_ESP_G2xyz_Px+ABX*I_ESP_Fxyz_Px;
    Double I_ESP_Fx2z_D2x = I_ESP_G2x2z_Px+ABX*I_ESP_Fx2z_Px;
    Double I_ESP_F3y_D2x = I_ESP_Gx3y_Px+ABX*I_ESP_F3y_Px;
    Double I_ESP_F2yz_D2x = I_ESP_Gx2yz_Px+ABX*I_ESP_F2yz_Px;
    Double I_ESP_Fy2z_D2x = I_ESP_Gxy2z_Px+ABX*I_ESP_Fy2z_Px;
    Double I_ESP_F3z_D2x = I_ESP_Gx3z_Px+ABX*I_ESP_F3z_Px;
    Double I_ESP_F3x_Dxy = I_ESP_G3xy_Px+ABY*I_ESP_F3x_Px;
    Double I_ESP_F2xy_Dxy = I_ESP_G2x2y_Px+ABY*I_ESP_F2xy_Px;
    Double I_ESP_F2xz_Dxy = I_ESP_G2xyz_Px+ABY*I_ESP_F2xz_Px;
    Double I_ESP_Fx2y_Dxy = I_ESP_Gx3y_Px+ABY*I_ESP_Fx2y_Px;
    Double I_ESP_Fxyz_Dxy = I_ESP_Gx2yz_Px+ABY*I_ESP_Fxyz_Px;
    Double I_ESP_Fx2z_Dxy = I_ESP_Gxy2z_Px+ABY*I_ESP_Fx2z_Px;
    Double I_ESP_F3y_Dxy = I_ESP_G4y_Px+ABY*I_ESP_F3y_Px;
    Double I_ESP_F2yz_Dxy = I_ESP_G3yz_Px+ABY*I_ESP_F2yz_Px;
    Double I_ESP_Fy2z_Dxy = I_ESP_G2y2z_Px+ABY*I_ESP_Fy2z_Px;
    Double I_ESP_F3z_Dxy = I_ESP_Gy3z_Px+ABY*I_ESP_F3z_Px;
    Double I_ESP_F3x_Dxz = I_ESP_G3xz_Px+ABZ*I_ESP_F3x_Px;
    Double I_ESP_F2xy_Dxz = I_ESP_G2xyz_Px+ABZ*I_ESP_F2xy_Px;
    Double I_ESP_F2xz_Dxz = I_ESP_G2x2z_Px+ABZ*I_ESP_F2xz_Px;
    Double I_ESP_Fx2y_Dxz = I_ESP_Gx2yz_Px+ABZ*I_ESP_Fx2y_Px;
    Double I_ESP_Fxyz_Dxz = I_ESP_Gxy2z_Px+ABZ*I_ESP_Fxyz_Px;
    Double I_ESP_Fx2z_Dxz = I_ESP_Gx3z_Px+ABZ*I_ESP_Fx2z_Px;
    Double I_ESP_F3y_Dxz = I_ESP_G3yz_Px+ABZ*I_ESP_F3y_Px;
    Double I_ESP_F2yz_Dxz = I_ESP_G2y2z_Px+ABZ*I_ESP_F2yz_Px;
    Double I_ESP_Fy2z_Dxz = I_ESP_Gy3z_Px+ABZ*I_ESP_Fy2z_Px;
    Double I_ESP_F3z_Dxz = I_ESP_G4z_Px+ABZ*I_ESP_F3z_Px;
    Double I_ESP_F3x_D2y = I_ESP_G3xy_Py+ABY*I_ESP_F3x_Py;
    Double I_ESP_F2xy_D2y = I_ESP_G2x2y_Py+ABY*I_ESP_F2xy_Py;
    Double I_ESP_F2xz_D2y = I_ESP_G2xyz_Py+ABY*I_ESP_F2xz_Py;
    Double I_ESP_Fx2y_D2y = I_ESP_Gx3y_Py+ABY*I_ESP_Fx2y_Py;
    Double I_ESP_Fxyz_D2y = I_ESP_Gx2yz_Py+ABY*I_ESP_Fxyz_Py;
    Double I_ESP_Fx2z_D2y = I_ESP_Gxy2z_Py+ABY*I_ESP_Fx2z_Py;
    Double I_ESP_F3y_D2y = I_ESP_G4y_Py+ABY*I_ESP_F3y_Py;
    Double I_ESP_F2yz_D2y = I_ESP_G3yz_Py+ABY*I_ESP_F2yz_Py;
    Double I_ESP_Fy2z_D2y = I_ESP_G2y2z_Py+ABY*I_ESP_Fy2z_Py;
    Double I_ESP_F3z_D2y = I_ESP_Gy3z_Py+ABY*I_ESP_F3z_Py;
    Double I_ESP_F3x_Dyz = I_ESP_G3xz_Py+ABZ*I_ESP_F3x_Py;
    Double I_ESP_F2xy_Dyz = I_ESP_G2xyz_Py+ABZ*I_ESP_F2xy_Py;
    Double I_ESP_F2xz_Dyz = I_ESP_G2x2z_Py+ABZ*I_ESP_F2xz_Py;
    Double I_ESP_Fx2y_Dyz = I_ESP_Gx2yz_Py+ABZ*I_ESP_Fx2y_Py;
    Double I_ESP_Fxyz_Dyz = I_ESP_Gxy2z_Py+ABZ*I_ESP_Fxyz_Py;
    Double I_ESP_Fx2z_Dyz = I_ESP_Gx3z_Py+ABZ*I_ESP_Fx2z_Py;
    Double I_ESP_F3y_Dyz = I_ESP_G3yz_Py+ABZ*I_ESP_F3y_Py;
    Double I_ESP_F2yz_Dyz = I_ESP_G2y2z_Py+ABZ*I_ESP_F2yz_Py;
    Double I_ESP_Fy2z_Dyz = I_ESP_Gy3z_Py+ABZ*I_ESP_Fy2z_Py;
    Double I_ESP_F3z_Dyz = I_ESP_G4z_Py+ABZ*I_ESP_F3z_Py;
    Double I_ESP_F3x_D2z = I_ESP_G3xz_Pz+ABZ*I_ESP_F3x_Pz;
    Double I_ESP_F2xy_D2z = I_ESP_G2xyz_Pz+ABZ*I_ESP_F2xy_Pz;
    Double I_ESP_F2xz_D2z = I_ESP_G2x2z_Pz+ABZ*I_ESP_F2xz_Pz;
    Double I_ESP_Fx2y_D2z = I_ESP_Gx2yz_Pz+ABZ*I_ESP_Fx2y_Pz;
    Double I_ESP_Fxyz_D2z = I_ESP_Gxy2z_Pz+ABZ*I_ESP_Fxyz_Pz;
    Double I_ESP_Fx2z_D2z = I_ESP_Gx3z_Pz+ABZ*I_ESP_Fx2z_Pz;
    Double I_ESP_F3y_D2z = I_ESP_G3yz_Pz+ABZ*I_ESP_F3y_Pz;
    Double I_ESP_F2yz_D2z = I_ESP_G2y2z_Pz+ABZ*I_ESP_F2yz_Pz;
    Double I_ESP_Fy2z_D2z = I_ESP_Gy3z_Pz+ABZ*I_ESP_Fy2z_Pz;
    Double I_ESP_F3z_D2z = I_ESP_G4z_Pz+ABZ*I_ESP_F3z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S_a
     * RHS shell quartet name: SQ_ESP_H_S_a
     ************************************************************/
    Double I_ESP_H5x_Px_a = I_ESP_I6x_S_a+ABX*I_ESP_H5x_S_a;
    Double I_ESP_H4xy_Px_a = I_ESP_I5xy_S_a+ABX*I_ESP_H4xy_S_a;
    Double I_ESP_H4xz_Px_a = I_ESP_I5xz_S_a+ABX*I_ESP_H4xz_S_a;
    Double I_ESP_H3x2y_Px_a = I_ESP_I4x2y_S_a+ABX*I_ESP_H3x2y_S_a;
    Double I_ESP_H3xyz_Px_a = I_ESP_I4xyz_S_a+ABX*I_ESP_H3xyz_S_a;
    Double I_ESP_H3x2z_Px_a = I_ESP_I4x2z_S_a+ABX*I_ESP_H3x2z_S_a;
    Double I_ESP_H2x3y_Px_a = I_ESP_I3x3y_S_a+ABX*I_ESP_H2x3y_S_a;
    Double I_ESP_H2x2yz_Px_a = I_ESP_I3x2yz_S_a+ABX*I_ESP_H2x2yz_S_a;
    Double I_ESP_H2xy2z_Px_a = I_ESP_I3xy2z_S_a+ABX*I_ESP_H2xy2z_S_a;
    Double I_ESP_H2x3z_Px_a = I_ESP_I3x3z_S_a+ABX*I_ESP_H2x3z_S_a;
    Double I_ESP_Hx4y_Px_a = I_ESP_I2x4y_S_a+ABX*I_ESP_Hx4y_S_a;
    Double I_ESP_Hx3yz_Px_a = I_ESP_I2x3yz_S_a+ABX*I_ESP_Hx3yz_S_a;
    Double I_ESP_Hx2y2z_Px_a = I_ESP_I2x2y2z_S_a+ABX*I_ESP_Hx2y2z_S_a;
    Double I_ESP_Hxy3z_Px_a = I_ESP_I2xy3z_S_a+ABX*I_ESP_Hxy3z_S_a;
    Double I_ESP_Hx4z_Px_a = I_ESP_I2x4z_S_a+ABX*I_ESP_Hx4z_S_a;
    Double I_ESP_H5y_Px_a = I_ESP_Ix5y_S_a+ABX*I_ESP_H5y_S_a;
    Double I_ESP_H4yz_Px_a = I_ESP_Ix4yz_S_a+ABX*I_ESP_H4yz_S_a;
    Double I_ESP_H3y2z_Px_a = I_ESP_Ix3y2z_S_a+ABX*I_ESP_H3y2z_S_a;
    Double I_ESP_H2y3z_Px_a = I_ESP_Ix2y3z_S_a+ABX*I_ESP_H2y3z_S_a;
    Double I_ESP_Hy4z_Px_a = I_ESP_Ixy4z_S_a+ABX*I_ESP_Hy4z_S_a;
    Double I_ESP_H5z_Px_a = I_ESP_Ix5z_S_a+ABX*I_ESP_H5z_S_a;
    Double I_ESP_H5x_Py_a = I_ESP_I5xy_S_a+ABY*I_ESP_H5x_S_a;
    Double I_ESP_H4xy_Py_a = I_ESP_I4x2y_S_a+ABY*I_ESP_H4xy_S_a;
    Double I_ESP_H4xz_Py_a = I_ESP_I4xyz_S_a+ABY*I_ESP_H4xz_S_a;
    Double I_ESP_H3x2y_Py_a = I_ESP_I3x3y_S_a+ABY*I_ESP_H3x2y_S_a;
    Double I_ESP_H3xyz_Py_a = I_ESP_I3x2yz_S_a+ABY*I_ESP_H3xyz_S_a;
    Double I_ESP_H3x2z_Py_a = I_ESP_I3xy2z_S_a+ABY*I_ESP_H3x2z_S_a;
    Double I_ESP_H2x3y_Py_a = I_ESP_I2x4y_S_a+ABY*I_ESP_H2x3y_S_a;
    Double I_ESP_H2x2yz_Py_a = I_ESP_I2x3yz_S_a+ABY*I_ESP_H2x2yz_S_a;
    Double I_ESP_H2xy2z_Py_a = I_ESP_I2x2y2z_S_a+ABY*I_ESP_H2xy2z_S_a;
    Double I_ESP_H2x3z_Py_a = I_ESP_I2xy3z_S_a+ABY*I_ESP_H2x3z_S_a;
    Double I_ESP_Hx4y_Py_a = I_ESP_Ix5y_S_a+ABY*I_ESP_Hx4y_S_a;
    Double I_ESP_Hx3yz_Py_a = I_ESP_Ix4yz_S_a+ABY*I_ESP_Hx3yz_S_a;
    Double I_ESP_Hx2y2z_Py_a = I_ESP_Ix3y2z_S_a+ABY*I_ESP_Hx2y2z_S_a;
    Double I_ESP_Hxy3z_Py_a = I_ESP_Ix2y3z_S_a+ABY*I_ESP_Hxy3z_S_a;
    Double I_ESP_Hx4z_Py_a = I_ESP_Ixy4z_S_a+ABY*I_ESP_Hx4z_S_a;
    Double I_ESP_H5y_Py_a = I_ESP_I6y_S_a+ABY*I_ESP_H5y_S_a;
    Double I_ESP_H4yz_Py_a = I_ESP_I5yz_S_a+ABY*I_ESP_H4yz_S_a;
    Double I_ESP_H3y2z_Py_a = I_ESP_I4y2z_S_a+ABY*I_ESP_H3y2z_S_a;
    Double I_ESP_H2y3z_Py_a = I_ESP_I3y3z_S_a+ABY*I_ESP_H2y3z_S_a;
    Double I_ESP_Hy4z_Py_a = I_ESP_I2y4z_S_a+ABY*I_ESP_Hy4z_S_a;
    Double I_ESP_H5z_Py_a = I_ESP_Iy5z_S_a+ABY*I_ESP_H5z_S_a;
    Double I_ESP_H5x_Pz_a = I_ESP_I5xz_S_a+ABZ*I_ESP_H5x_S_a;
    Double I_ESP_H4xy_Pz_a = I_ESP_I4xyz_S_a+ABZ*I_ESP_H4xy_S_a;
    Double I_ESP_H4xz_Pz_a = I_ESP_I4x2z_S_a+ABZ*I_ESP_H4xz_S_a;
    Double I_ESP_H3x2y_Pz_a = I_ESP_I3x2yz_S_a+ABZ*I_ESP_H3x2y_S_a;
    Double I_ESP_H3xyz_Pz_a = I_ESP_I3xy2z_S_a+ABZ*I_ESP_H3xyz_S_a;
    Double I_ESP_H3x2z_Pz_a = I_ESP_I3x3z_S_a+ABZ*I_ESP_H3x2z_S_a;
    Double I_ESP_H2x3y_Pz_a = I_ESP_I2x3yz_S_a+ABZ*I_ESP_H2x3y_S_a;
    Double I_ESP_H2x2yz_Pz_a = I_ESP_I2x2y2z_S_a+ABZ*I_ESP_H2x2yz_S_a;
    Double I_ESP_H2xy2z_Pz_a = I_ESP_I2xy3z_S_a+ABZ*I_ESP_H2xy2z_S_a;
    Double I_ESP_H2x3z_Pz_a = I_ESP_I2x4z_S_a+ABZ*I_ESP_H2x3z_S_a;
    Double I_ESP_Hx4y_Pz_a = I_ESP_Ix4yz_S_a+ABZ*I_ESP_Hx4y_S_a;
    Double I_ESP_Hx3yz_Pz_a = I_ESP_Ix3y2z_S_a+ABZ*I_ESP_Hx3yz_S_a;
    Double I_ESP_Hx2y2z_Pz_a = I_ESP_Ix2y3z_S_a+ABZ*I_ESP_Hx2y2z_S_a;
    Double I_ESP_Hxy3z_Pz_a = I_ESP_Ixy4z_S_a+ABZ*I_ESP_Hxy3z_S_a;
    Double I_ESP_Hx4z_Pz_a = I_ESP_Ix5z_S_a+ABZ*I_ESP_Hx4z_S_a;
    Double I_ESP_H5y_Pz_a = I_ESP_I5yz_S_a+ABZ*I_ESP_H5y_S_a;
    Double I_ESP_H4yz_Pz_a = I_ESP_I4y2z_S_a+ABZ*I_ESP_H4yz_S_a;
    Double I_ESP_H3y2z_Pz_a = I_ESP_I3y3z_S_a+ABZ*I_ESP_H3y2z_S_a;
    Double I_ESP_H2y3z_Pz_a = I_ESP_I2y4z_S_a+ABZ*I_ESP_H2y3z_S_a;
    Double I_ESP_Hy4z_Pz_a = I_ESP_Iy5z_S_a+ABZ*I_ESP_Hy4z_S_a;
    Double I_ESP_H5z_Pz_a = I_ESP_I6z_S_a+ABZ*I_ESP_H5z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 8 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_S_a
     * RHS shell quartet name: SQ_ESP_I_S_a
     ************************************************************/
    Double I_ESP_I6x_Px_a = I_ESP_K7x_S_a+ABX*I_ESP_I6x_S_a;
    Double I_ESP_I5xy_Px_a = I_ESP_K6xy_S_a+ABX*I_ESP_I5xy_S_a;
    Double I_ESP_I5xz_Px_a = I_ESP_K6xz_S_a+ABX*I_ESP_I5xz_S_a;
    Double I_ESP_I4x2y_Px_a = I_ESP_K5x2y_S_a+ABX*I_ESP_I4x2y_S_a;
    Double I_ESP_I4xyz_Px_a = I_ESP_K5xyz_S_a+ABX*I_ESP_I4xyz_S_a;
    Double I_ESP_I4x2z_Px_a = I_ESP_K5x2z_S_a+ABX*I_ESP_I4x2z_S_a;
    Double I_ESP_I3x3y_Px_a = I_ESP_K4x3y_S_a+ABX*I_ESP_I3x3y_S_a;
    Double I_ESP_I3x2yz_Px_a = I_ESP_K4x2yz_S_a+ABX*I_ESP_I3x2yz_S_a;
    Double I_ESP_I3xy2z_Px_a = I_ESP_K4xy2z_S_a+ABX*I_ESP_I3xy2z_S_a;
    Double I_ESP_I3x3z_Px_a = I_ESP_K4x3z_S_a+ABX*I_ESP_I3x3z_S_a;
    Double I_ESP_I2x4y_Px_a = I_ESP_K3x4y_S_a+ABX*I_ESP_I2x4y_S_a;
    Double I_ESP_I2x3yz_Px_a = I_ESP_K3x3yz_S_a+ABX*I_ESP_I2x3yz_S_a;
    Double I_ESP_I2x2y2z_Px_a = I_ESP_K3x2y2z_S_a+ABX*I_ESP_I2x2y2z_S_a;
    Double I_ESP_I2xy3z_Px_a = I_ESP_K3xy3z_S_a+ABX*I_ESP_I2xy3z_S_a;
    Double I_ESP_I2x4z_Px_a = I_ESP_K3x4z_S_a+ABX*I_ESP_I2x4z_S_a;
    Double I_ESP_Ix5y_Px_a = I_ESP_K2x5y_S_a+ABX*I_ESP_Ix5y_S_a;
    Double I_ESP_Ix4yz_Px_a = I_ESP_K2x4yz_S_a+ABX*I_ESP_Ix4yz_S_a;
    Double I_ESP_Ix3y2z_Px_a = I_ESP_K2x3y2z_S_a+ABX*I_ESP_Ix3y2z_S_a;
    Double I_ESP_Ix2y3z_Px_a = I_ESP_K2x2y3z_S_a+ABX*I_ESP_Ix2y3z_S_a;
    Double I_ESP_Ixy4z_Px_a = I_ESP_K2xy4z_S_a+ABX*I_ESP_Ixy4z_S_a;
    Double I_ESP_Ix5z_Px_a = I_ESP_K2x5z_S_a+ABX*I_ESP_Ix5z_S_a;
    Double I_ESP_I6y_Px_a = I_ESP_Kx6y_S_a+ABX*I_ESP_I6y_S_a;
    Double I_ESP_I5yz_Px_a = I_ESP_Kx5yz_S_a+ABX*I_ESP_I5yz_S_a;
    Double I_ESP_I4y2z_Px_a = I_ESP_Kx4y2z_S_a+ABX*I_ESP_I4y2z_S_a;
    Double I_ESP_I3y3z_Px_a = I_ESP_Kx3y3z_S_a+ABX*I_ESP_I3y3z_S_a;
    Double I_ESP_I2y4z_Px_a = I_ESP_Kx2y4z_S_a+ABX*I_ESP_I2y4z_S_a;
    Double I_ESP_Iy5z_Px_a = I_ESP_Kxy5z_S_a+ABX*I_ESP_Iy5z_S_a;
    Double I_ESP_I6z_Px_a = I_ESP_Kx6z_S_a+ABX*I_ESP_I6z_S_a;
    Double I_ESP_I5xy_Py_a = I_ESP_K5x2y_S_a+ABY*I_ESP_I5xy_S_a;
    Double I_ESP_I5xz_Py_a = I_ESP_K5xyz_S_a+ABY*I_ESP_I5xz_S_a;
    Double I_ESP_I4x2y_Py_a = I_ESP_K4x3y_S_a+ABY*I_ESP_I4x2y_S_a;
    Double I_ESP_I4xyz_Py_a = I_ESP_K4x2yz_S_a+ABY*I_ESP_I4xyz_S_a;
    Double I_ESP_I4x2z_Py_a = I_ESP_K4xy2z_S_a+ABY*I_ESP_I4x2z_S_a;
    Double I_ESP_I3x3y_Py_a = I_ESP_K3x4y_S_a+ABY*I_ESP_I3x3y_S_a;
    Double I_ESP_I3x2yz_Py_a = I_ESP_K3x3yz_S_a+ABY*I_ESP_I3x2yz_S_a;
    Double I_ESP_I3xy2z_Py_a = I_ESP_K3x2y2z_S_a+ABY*I_ESP_I3xy2z_S_a;
    Double I_ESP_I3x3z_Py_a = I_ESP_K3xy3z_S_a+ABY*I_ESP_I3x3z_S_a;
    Double I_ESP_I2x4y_Py_a = I_ESP_K2x5y_S_a+ABY*I_ESP_I2x4y_S_a;
    Double I_ESP_I2x3yz_Py_a = I_ESP_K2x4yz_S_a+ABY*I_ESP_I2x3yz_S_a;
    Double I_ESP_I2x2y2z_Py_a = I_ESP_K2x3y2z_S_a+ABY*I_ESP_I2x2y2z_S_a;
    Double I_ESP_I2xy3z_Py_a = I_ESP_K2x2y3z_S_a+ABY*I_ESP_I2xy3z_S_a;
    Double I_ESP_I2x4z_Py_a = I_ESP_K2xy4z_S_a+ABY*I_ESP_I2x4z_S_a;
    Double I_ESP_Ix5y_Py_a = I_ESP_Kx6y_S_a+ABY*I_ESP_Ix5y_S_a;
    Double I_ESP_Ix4yz_Py_a = I_ESP_Kx5yz_S_a+ABY*I_ESP_Ix4yz_S_a;
    Double I_ESP_Ix3y2z_Py_a = I_ESP_Kx4y2z_S_a+ABY*I_ESP_Ix3y2z_S_a;
    Double I_ESP_Ix2y3z_Py_a = I_ESP_Kx3y3z_S_a+ABY*I_ESP_Ix2y3z_S_a;
    Double I_ESP_Ixy4z_Py_a = I_ESP_Kx2y4z_S_a+ABY*I_ESP_Ixy4z_S_a;
    Double I_ESP_Ix5z_Py_a = I_ESP_Kxy5z_S_a+ABY*I_ESP_Ix5z_S_a;
    Double I_ESP_I6y_Py_a = I_ESP_K7y_S_a+ABY*I_ESP_I6y_S_a;
    Double I_ESP_I5yz_Py_a = I_ESP_K6yz_S_a+ABY*I_ESP_I5yz_S_a;
    Double I_ESP_I4y2z_Py_a = I_ESP_K5y2z_S_a+ABY*I_ESP_I4y2z_S_a;
    Double I_ESP_I3y3z_Py_a = I_ESP_K4y3z_S_a+ABY*I_ESP_I3y3z_S_a;
    Double I_ESP_I2y4z_Py_a = I_ESP_K3y4z_S_a+ABY*I_ESP_I2y4z_S_a;
    Double I_ESP_Iy5z_Py_a = I_ESP_K2y5z_S_a+ABY*I_ESP_Iy5z_S_a;
    Double I_ESP_I6z_Py_a = I_ESP_Ky6z_S_a+ABY*I_ESP_I6z_S_a;
    Double I_ESP_I5xz_Pz_a = I_ESP_K5x2z_S_a+ABZ*I_ESP_I5xz_S_a;
    Double I_ESP_I4xyz_Pz_a = I_ESP_K4xy2z_S_a+ABZ*I_ESP_I4xyz_S_a;
    Double I_ESP_I4x2z_Pz_a = I_ESP_K4x3z_S_a+ABZ*I_ESP_I4x2z_S_a;
    Double I_ESP_I3x2yz_Pz_a = I_ESP_K3x2y2z_S_a+ABZ*I_ESP_I3x2yz_S_a;
    Double I_ESP_I3xy2z_Pz_a = I_ESP_K3xy3z_S_a+ABZ*I_ESP_I3xy2z_S_a;
    Double I_ESP_I3x3z_Pz_a = I_ESP_K3x4z_S_a+ABZ*I_ESP_I3x3z_S_a;
    Double I_ESP_I2x3yz_Pz_a = I_ESP_K2x3y2z_S_a+ABZ*I_ESP_I2x3yz_S_a;
    Double I_ESP_I2x2y2z_Pz_a = I_ESP_K2x2y3z_S_a+ABZ*I_ESP_I2x2y2z_S_a;
    Double I_ESP_I2xy3z_Pz_a = I_ESP_K2xy4z_S_a+ABZ*I_ESP_I2xy3z_S_a;
    Double I_ESP_I2x4z_Pz_a = I_ESP_K2x5z_S_a+ABZ*I_ESP_I2x4z_S_a;
    Double I_ESP_Ix4yz_Pz_a = I_ESP_Kx4y2z_S_a+ABZ*I_ESP_Ix4yz_S_a;
    Double I_ESP_Ix3y2z_Pz_a = I_ESP_Kx3y3z_S_a+ABZ*I_ESP_Ix3y2z_S_a;
    Double I_ESP_Ix2y3z_Pz_a = I_ESP_Kx2y4z_S_a+ABZ*I_ESP_Ix2y3z_S_a;
    Double I_ESP_Ixy4z_Pz_a = I_ESP_Kxy5z_S_a+ABZ*I_ESP_Ixy4z_S_a;
    Double I_ESP_Ix5z_Pz_a = I_ESP_Kx6z_S_a+ABZ*I_ESP_Ix5z_S_a;
    Double I_ESP_I5yz_Pz_a = I_ESP_K5y2z_S_a+ABZ*I_ESP_I5yz_S_a;
    Double I_ESP_I4y2z_Pz_a = I_ESP_K4y3z_S_a+ABZ*I_ESP_I4y2z_S_a;
    Double I_ESP_I3y3z_Pz_a = I_ESP_K3y4z_S_a+ABZ*I_ESP_I3y3z_S_a;
    Double I_ESP_I2y4z_Pz_a = I_ESP_K2y5z_S_a+ABZ*I_ESP_I2y4z_S_a;
    Double I_ESP_Iy5z_Pz_a = I_ESP_Ky6z_S_a+ABZ*I_ESP_Iy5z_S_a;
    Double I_ESP_I6z_Pz_a = I_ESP_K7z_S_a+ABZ*I_ESP_I6z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P_a
     * RHS shell quartet name: SQ_ESP_H_P_a
     ************************************************************/
    Double I_ESP_H5x_D2x_a = I_ESP_I6x_Px_a+ABX*I_ESP_H5x_Px_a;
    Double I_ESP_H4xy_D2x_a = I_ESP_I5xy_Px_a+ABX*I_ESP_H4xy_Px_a;
    Double I_ESP_H4xz_D2x_a = I_ESP_I5xz_Px_a+ABX*I_ESP_H4xz_Px_a;
    Double I_ESP_H3x2y_D2x_a = I_ESP_I4x2y_Px_a+ABX*I_ESP_H3x2y_Px_a;
    Double I_ESP_H3xyz_D2x_a = I_ESP_I4xyz_Px_a+ABX*I_ESP_H3xyz_Px_a;
    Double I_ESP_H3x2z_D2x_a = I_ESP_I4x2z_Px_a+ABX*I_ESP_H3x2z_Px_a;
    Double I_ESP_H2x3y_D2x_a = I_ESP_I3x3y_Px_a+ABX*I_ESP_H2x3y_Px_a;
    Double I_ESP_H2x2yz_D2x_a = I_ESP_I3x2yz_Px_a+ABX*I_ESP_H2x2yz_Px_a;
    Double I_ESP_H2xy2z_D2x_a = I_ESP_I3xy2z_Px_a+ABX*I_ESP_H2xy2z_Px_a;
    Double I_ESP_H2x3z_D2x_a = I_ESP_I3x3z_Px_a+ABX*I_ESP_H2x3z_Px_a;
    Double I_ESP_Hx4y_D2x_a = I_ESP_I2x4y_Px_a+ABX*I_ESP_Hx4y_Px_a;
    Double I_ESP_Hx3yz_D2x_a = I_ESP_I2x3yz_Px_a+ABX*I_ESP_Hx3yz_Px_a;
    Double I_ESP_Hx2y2z_D2x_a = I_ESP_I2x2y2z_Px_a+ABX*I_ESP_Hx2y2z_Px_a;
    Double I_ESP_Hxy3z_D2x_a = I_ESP_I2xy3z_Px_a+ABX*I_ESP_Hxy3z_Px_a;
    Double I_ESP_Hx4z_D2x_a = I_ESP_I2x4z_Px_a+ABX*I_ESP_Hx4z_Px_a;
    Double I_ESP_H5y_D2x_a = I_ESP_Ix5y_Px_a+ABX*I_ESP_H5y_Px_a;
    Double I_ESP_H4yz_D2x_a = I_ESP_Ix4yz_Px_a+ABX*I_ESP_H4yz_Px_a;
    Double I_ESP_H3y2z_D2x_a = I_ESP_Ix3y2z_Px_a+ABX*I_ESP_H3y2z_Px_a;
    Double I_ESP_H2y3z_D2x_a = I_ESP_Ix2y3z_Px_a+ABX*I_ESP_H2y3z_Px_a;
    Double I_ESP_Hy4z_D2x_a = I_ESP_Ixy4z_Px_a+ABX*I_ESP_Hy4z_Px_a;
    Double I_ESP_H5z_D2x_a = I_ESP_Ix5z_Px_a+ABX*I_ESP_H5z_Px_a;
    Double I_ESP_H5x_Dxy_a = I_ESP_I5xy_Px_a+ABY*I_ESP_H5x_Px_a;
    Double I_ESP_H4xy_Dxy_a = I_ESP_I4x2y_Px_a+ABY*I_ESP_H4xy_Px_a;
    Double I_ESP_H4xz_Dxy_a = I_ESP_I4xyz_Px_a+ABY*I_ESP_H4xz_Px_a;
    Double I_ESP_H3x2y_Dxy_a = I_ESP_I3x3y_Px_a+ABY*I_ESP_H3x2y_Px_a;
    Double I_ESP_H3xyz_Dxy_a = I_ESP_I3x2yz_Px_a+ABY*I_ESP_H3xyz_Px_a;
    Double I_ESP_H3x2z_Dxy_a = I_ESP_I3xy2z_Px_a+ABY*I_ESP_H3x2z_Px_a;
    Double I_ESP_H2x3y_Dxy_a = I_ESP_I2x4y_Px_a+ABY*I_ESP_H2x3y_Px_a;
    Double I_ESP_H2x2yz_Dxy_a = I_ESP_I2x3yz_Px_a+ABY*I_ESP_H2x2yz_Px_a;
    Double I_ESP_H2xy2z_Dxy_a = I_ESP_I2x2y2z_Px_a+ABY*I_ESP_H2xy2z_Px_a;
    Double I_ESP_H2x3z_Dxy_a = I_ESP_I2xy3z_Px_a+ABY*I_ESP_H2x3z_Px_a;
    Double I_ESP_Hx4y_Dxy_a = I_ESP_Ix5y_Px_a+ABY*I_ESP_Hx4y_Px_a;
    Double I_ESP_Hx3yz_Dxy_a = I_ESP_Ix4yz_Px_a+ABY*I_ESP_Hx3yz_Px_a;
    Double I_ESP_Hx2y2z_Dxy_a = I_ESP_Ix3y2z_Px_a+ABY*I_ESP_Hx2y2z_Px_a;
    Double I_ESP_Hxy3z_Dxy_a = I_ESP_Ix2y3z_Px_a+ABY*I_ESP_Hxy3z_Px_a;
    Double I_ESP_Hx4z_Dxy_a = I_ESP_Ixy4z_Px_a+ABY*I_ESP_Hx4z_Px_a;
    Double I_ESP_H5y_Dxy_a = I_ESP_I6y_Px_a+ABY*I_ESP_H5y_Px_a;
    Double I_ESP_H4yz_Dxy_a = I_ESP_I5yz_Px_a+ABY*I_ESP_H4yz_Px_a;
    Double I_ESP_H3y2z_Dxy_a = I_ESP_I4y2z_Px_a+ABY*I_ESP_H3y2z_Px_a;
    Double I_ESP_H2y3z_Dxy_a = I_ESP_I3y3z_Px_a+ABY*I_ESP_H2y3z_Px_a;
    Double I_ESP_Hy4z_Dxy_a = I_ESP_I2y4z_Px_a+ABY*I_ESP_Hy4z_Px_a;
    Double I_ESP_H5z_Dxy_a = I_ESP_Iy5z_Px_a+ABY*I_ESP_H5z_Px_a;
    Double I_ESP_H5x_Dxz_a = I_ESP_I5xz_Px_a+ABZ*I_ESP_H5x_Px_a;
    Double I_ESP_H4xy_Dxz_a = I_ESP_I4xyz_Px_a+ABZ*I_ESP_H4xy_Px_a;
    Double I_ESP_H4xz_Dxz_a = I_ESP_I4x2z_Px_a+ABZ*I_ESP_H4xz_Px_a;
    Double I_ESP_H3x2y_Dxz_a = I_ESP_I3x2yz_Px_a+ABZ*I_ESP_H3x2y_Px_a;
    Double I_ESP_H3xyz_Dxz_a = I_ESP_I3xy2z_Px_a+ABZ*I_ESP_H3xyz_Px_a;
    Double I_ESP_H3x2z_Dxz_a = I_ESP_I3x3z_Px_a+ABZ*I_ESP_H3x2z_Px_a;
    Double I_ESP_H2x3y_Dxz_a = I_ESP_I2x3yz_Px_a+ABZ*I_ESP_H2x3y_Px_a;
    Double I_ESP_H2x2yz_Dxz_a = I_ESP_I2x2y2z_Px_a+ABZ*I_ESP_H2x2yz_Px_a;
    Double I_ESP_H2xy2z_Dxz_a = I_ESP_I2xy3z_Px_a+ABZ*I_ESP_H2xy2z_Px_a;
    Double I_ESP_H2x3z_Dxz_a = I_ESP_I2x4z_Px_a+ABZ*I_ESP_H2x3z_Px_a;
    Double I_ESP_Hx4y_Dxz_a = I_ESP_Ix4yz_Px_a+ABZ*I_ESP_Hx4y_Px_a;
    Double I_ESP_Hx3yz_Dxz_a = I_ESP_Ix3y2z_Px_a+ABZ*I_ESP_Hx3yz_Px_a;
    Double I_ESP_Hx2y2z_Dxz_a = I_ESP_Ix2y3z_Px_a+ABZ*I_ESP_Hx2y2z_Px_a;
    Double I_ESP_Hxy3z_Dxz_a = I_ESP_Ixy4z_Px_a+ABZ*I_ESP_Hxy3z_Px_a;
    Double I_ESP_Hx4z_Dxz_a = I_ESP_Ix5z_Px_a+ABZ*I_ESP_Hx4z_Px_a;
    Double I_ESP_H5y_Dxz_a = I_ESP_I5yz_Px_a+ABZ*I_ESP_H5y_Px_a;
    Double I_ESP_H4yz_Dxz_a = I_ESP_I4y2z_Px_a+ABZ*I_ESP_H4yz_Px_a;
    Double I_ESP_H3y2z_Dxz_a = I_ESP_I3y3z_Px_a+ABZ*I_ESP_H3y2z_Px_a;
    Double I_ESP_H2y3z_Dxz_a = I_ESP_I2y4z_Px_a+ABZ*I_ESP_H2y3z_Px_a;
    Double I_ESP_Hy4z_Dxz_a = I_ESP_Iy5z_Px_a+ABZ*I_ESP_Hy4z_Px_a;
    Double I_ESP_H5z_Dxz_a = I_ESP_I6z_Px_a+ABZ*I_ESP_H5z_Px_a;
    Double I_ESP_H5x_D2y_a = I_ESP_I5xy_Py_a+ABY*I_ESP_H5x_Py_a;
    Double I_ESP_H4xy_D2y_a = I_ESP_I4x2y_Py_a+ABY*I_ESP_H4xy_Py_a;
    Double I_ESP_H4xz_D2y_a = I_ESP_I4xyz_Py_a+ABY*I_ESP_H4xz_Py_a;
    Double I_ESP_H3x2y_D2y_a = I_ESP_I3x3y_Py_a+ABY*I_ESP_H3x2y_Py_a;
    Double I_ESP_H3xyz_D2y_a = I_ESP_I3x2yz_Py_a+ABY*I_ESP_H3xyz_Py_a;
    Double I_ESP_H3x2z_D2y_a = I_ESP_I3xy2z_Py_a+ABY*I_ESP_H3x2z_Py_a;
    Double I_ESP_H2x3y_D2y_a = I_ESP_I2x4y_Py_a+ABY*I_ESP_H2x3y_Py_a;
    Double I_ESP_H2x2yz_D2y_a = I_ESP_I2x3yz_Py_a+ABY*I_ESP_H2x2yz_Py_a;
    Double I_ESP_H2xy2z_D2y_a = I_ESP_I2x2y2z_Py_a+ABY*I_ESP_H2xy2z_Py_a;
    Double I_ESP_H2x3z_D2y_a = I_ESP_I2xy3z_Py_a+ABY*I_ESP_H2x3z_Py_a;
    Double I_ESP_Hx4y_D2y_a = I_ESP_Ix5y_Py_a+ABY*I_ESP_Hx4y_Py_a;
    Double I_ESP_Hx3yz_D2y_a = I_ESP_Ix4yz_Py_a+ABY*I_ESP_Hx3yz_Py_a;
    Double I_ESP_Hx2y2z_D2y_a = I_ESP_Ix3y2z_Py_a+ABY*I_ESP_Hx2y2z_Py_a;
    Double I_ESP_Hxy3z_D2y_a = I_ESP_Ix2y3z_Py_a+ABY*I_ESP_Hxy3z_Py_a;
    Double I_ESP_Hx4z_D2y_a = I_ESP_Ixy4z_Py_a+ABY*I_ESP_Hx4z_Py_a;
    Double I_ESP_H5y_D2y_a = I_ESP_I6y_Py_a+ABY*I_ESP_H5y_Py_a;
    Double I_ESP_H4yz_D2y_a = I_ESP_I5yz_Py_a+ABY*I_ESP_H4yz_Py_a;
    Double I_ESP_H3y2z_D2y_a = I_ESP_I4y2z_Py_a+ABY*I_ESP_H3y2z_Py_a;
    Double I_ESP_H2y3z_D2y_a = I_ESP_I3y3z_Py_a+ABY*I_ESP_H2y3z_Py_a;
    Double I_ESP_Hy4z_D2y_a = I_ESP_I2y4z_Py_a+ABY*I_ESP_Hy4z_Py_a;
    Double I_ESP_H5z_D2y_a = I_ESP_Iy5z_Py_a+ABY*I_ESP_H5z_Py_a;
    Double I_ESP_H5x_Dyz_a = I_ESP_I5xz_Py_a+ABZ*I_ESP_H5x_Py_a;
    Double I_ESP_H4xy_Dyz_a = I_ESP_I4xyz_Py_a+ABZ*I_ESP_H4xy_Py_a;
    Double I_ESP_H4xz_Dyz_a = I_ESP_I4x2z_Py_a+ABZ*I_ESP_H4xz_Py_a;
    Double I_ESP_H3x2y_Dyz_a = I_ESP_I3x2yz_Py_a+ABZ*I_ESP_H3x2y_Py_a;
    Double I_ESP_H3xyz_Dyz_a = I_ESP_I3xy2z_Py_a+ABZ*I_ESP_H3xyz_Py_a;
    Double I_ESP_H3x2z_Dyz_a = I_ESP_I3x3z_Py_a+ABZ*I_ESP_H3x2z_Py_a;
    Double I_ESP_H2x3y_Dyz_a = I_ESP_I2x3yz_Py_a+ABZ*I_ESP_H2x3y_Py_a;
    Double I_ESP_H2x2yz_Dyz_a = I_ESP_I2x2y2z_Py_a+ABZ*I_ESP_H2x2yz_Py_a;
    Double I_ESP_H2xy2z_Dyz_a = I_ESP_I2xy3z_Py_a+ABZ*I_ESP_H2xy2z_Py_a;
    Double I_ESP_H2x3z_Dyz_a = I_ESP_I2x4z_Py_a+ABZ*I_ESP_H2x3z_Py_a;
    Double I_ESP_Hx4y_Dyz_a = I_ESP_Ix4yz_Py_a+ABZ*I_ESP_Hx4y_Py_a;
    Double I_ESP_Hx3yz_Dyz_a = I_ESP_Ix3y2z_Py_a+ABZ*I_ESP_Hx3yz_Py_a;
    Double I_ESP_Hx2y2z_Dyz_a = I_ESP_Ix2y3z_Py_a+ABZ*I_ESP_Hx2y2z_Py_a;
    Double I_ESP_Hxy3z_Dyz_a = I_ESP_Ixy4z_Py_a+ABZ*I_ESP_Hxy3z_Py_a;
    Double I_ESP_Hx4z_Dyz_a = I_ESP_Ix5z_Py_a+ABZ*I_ESP_Hx4z_Py_a;
    Double I_ESP_H5y_Dyz_a = I_ESP_I5yz_Py_a+ABZ*I_ESP_H5y_Py_a;
    Double I_ESP_H4yz_Dyz_a = I_ESP_I4y2z_Py_a+ABZ*I_ESP_H4yz_Py_a;
    Double I_ESP_H3y2z_Dyz_a = I_ESP_I3y3z_Py_a+ABZ*I_ESP_H3y2z_Py_a;
    Double I_ESP_H2y3z_Dyz_a = I_ESP_I2y4z_Py_a+ABZ*I_ESP_H2y3z_Py_a;
    Double I_ESP_Hy4z_Dyz_a = I_ESP_Iy5z_Py_a+ABZ*I_ESP_Hy4z_Py_a;
    Double I_ESP_H5z_Dyz_a = I_ESP_I6z_Py_a+ABZ*I_ESP_H5z_Py_a;
    Double I_ESP_H5x_D2z_a = I_ESP_I5xz_Pz_a+ABZ*I_ESP_H5x_Pz_a;
    Double I_ESP_H4xy_D2z_a = I_ESP_I4xyz_Pz_a+ABZ*I_ESP_H4xy_Pz_a;
    Double I_ESP_H4xz_D2z_a = I_ESP_I4x2z_Pz_a+ABZ*I_ESP_H4xz_Pz_a;
    Double I_ESP_H3x2y_D2z_a = I_ESP_I3x2yz_Pz_a+ABZ*I_ESP_H3x2y_Pz_a;
    Double I_ESP_H3xyz_D2z_a = I_ESP_I3xy2z_Pz_a+ABZ*I_ESP_H3xyz_Pz_a;
    Double I_ESP_H3x2z_D2z_a = I_ESP_I3x3z_Pz_a+ABZ*I_ESP_H3x2z_Pz_a;
    Double I_ESP_H2x3y_D2z_a = I_ESP_I2x3yz_Pz_a+ABZ*I_ESP_H2x3y_Pz_a;
    Double I_ESP_H2x2yz_D2z_a = I_ESP_I2x2y2z_Pz_a+ABZ*I_ESP_H2x2yz_Pz_a;
    Double I_ESP_H2xy2z_D2z_a = I_ESP_I2xy3z_Pz_a+ABZ*I_ESP_H2xy2z_Pz_a;
    Double I_ESP_H2x3z_D2z_a = I_ESP_I2x4z_Pz_a+ABZ*I_ESP_H2x3z_Pz_a;
    Double I_ESP_Hx4y_D2z_a = I_ESP_Ix4yz_Pz_a+ABZ*I_ESP_Hx4y_Pz_a;
    Double I_ESP_Hx3yz_D2z_a = I_ESP_Ix3y2z_Pz_a+ABZ*I_ESP_Hx3yz_Pz_a;
    Double I_ESP_Hx2y2z_D2z_a = I_ESP_Ix2y3z_Pz_a+ABZ*I_ESP_Hx2y2z_Pz_a;
    Double I_ESP_Hxy3z_D2z_a = I_ESP_Ixy4z_Pz_a+ABZ*I_ESP_Hxy3z_Pz_a;
    Double I_ESP_Hx4z_D2z_a = I_ESP_Ix5z_Pz_a+ABZ*I_ESP_Hx4z_Pz_a;
    Double I_ESP_H5y_D2z_a = I_ESP_I5yz_Pz_a+ABZ*I_ESP_H5y_Pz_a;
    Double I_ESP_H4yz_D2z_a = I_ESP_I4y2z_Pz_a+ABZ*I_ESP_H4yz_Pz_a;
    Double I_ESP_H3y2z_D2z_a = I_ESP_I3y3z_Pz_a+ABZ*I_ESP_H3y2z_Pz_a;
    Double I_ESP_H2y3z_D2z_a = I_ESP_I2y4z_Pz_a+ABZ*I_ESP_H2y3z_Pz_a;
    Double I_ESP_Hy4z_D2z_a = I_ESP_Iy5z_Pz_a+ABZ*I_ESP_Hy4z_Pz_a;
    Double I_ESP_H5z_D2z_a = I_ESP_I6z_Pz_a+ABZ*I_ESP_H5z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_K_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_S_aa
     * RHS shell quartet name: SQ_ESP_K_S_aa
     ************************************************************/
    Double I_ESP_K7x_Px_aa = I_ESP_L8x_S_aa+ABX*I_ESP_K7x_S_aa;
    Double I_ESP_K6xy_Px_aa = I_ESP_L7xy_S_aa+ABX*I_ESP_K6xy_S_aa;
    Double I_ESP_K6xz_Px_aa = I_ESP_L7xz_S_aa+ABX*I_ESP_K6xz_S_aa;
    Double I_ESP_K5x2y_Px_aa = I_ESP_L6x2y_S_aa+ABX*I_ESP_K5x2y_S_aa;
    Double I_ESP_K5xyz_Px_aa = I_ESP_L6xyz_S_aa+ABX*I_ESP_K5xyz_S_aa;
    Double I_ESP_K5x2z_Px_aa = I_ESP_L6x2z_S_aa+ABX*I_ESP_K5x2z_S_aa;
    Double I_ESP_K4x3y_Px_aa = I_ESP_L5x3y_S_aa+ABX*I_ESP_K4x3y_S_aa;
    Double I_ESP_K4x2yz_Px_aa = I_ESP_L5x2yz_S_aa+ABX*I_ESP_K4x2yz_S_aa;
    Double I_ESP_K4xy2z_Px_aa = I_ESP_L5xy2z_S_aa+ABX*I_ESP_K4xy2z_S_aa;
    Double I_ESP_K4x3z_Px_aa = I_ESP_L5x3z_S_aa+ABX*I_ESP_K4x3z_S_aa;
    Double I_ESP_K3x4y_Px_aa = I_ESP_L4x4y_S_aa+ABX*I_ESP_K3x4y_S_aa;
    Double I_ESP_K3x3yz_Px_aa = I_ESP_L4x3yz_S_aa+ABX*I_ESP_K3x3yz_S_aa;
    Double I_ESP_K3x2y2z_Px_aa = I_ESP_L4x2y2z_S_aa+ABX*I_ESP_K3x2y2z_S_aa;
    Double I_ESP_K3xy3z_Px_aa = I_ESP_L4xy3z_S_aa+ABX*I_ESP_K3xy3z_S_aa;
    Double I_ESP_K3x4z_Px_aa = I_ESP_L4x4z_S_aa+ABX*I_ESP_K3x4z_S_aa;
    Double I_ESP_K2x5y_Px_aa = I_ESP_L3x5y_S_aa+ABX*I_ESP_K2x5y_S_aa;
    Double I_ESP_K2x4yz_Px_aa = I_ESP_L3x4yz_S_aa+ABX*I_ESP_K2x4yz_S_aa;
    Double I_ESP_K2x3y2z_Px_aa = I_ESP_L3x3y2z_S_aa+ABX*I_ESP_K2x3y2z_S_aa;
    Double I_ESP_K2x2y3z_Px_aa = I_ESP_L3x2y3z_S_aa+ABX*I_ESP_K2x2y3z_S_aa;
    Double I_ESP_K2xy4z_Px_aa = I_ESP_L3xy4z_S_aa+ABX*I_ESP_K2xy4z_S_aa;
    Double I_ESP_K2x5z_Px_aa = I_ESP_L3x5z_S_aa+ABX*I_ESP_K2x5z_S_aa;
    Double I_ESP_Kx6y_Px_aa = I_ESP_L2x6y_S_aa+ABX*I_ESP_Kx6y_S_aa;
    Double I_ESP_Kx5yz_Px_aa = I_ESP_L2x5yz_S_aa+ABX*I_ESP_Kx5yz_S_aa;
    Double I_ESP_Kx4y2z_Px_aa = I_ESP_L2x4y2z_S_aa+ABX*I_ESP_Kx4y2z_S_aa;
    Double I_ESP_Kx3y3z_Px_aa = I_ESP_L2x3y3z_S_aa+ABX*I_ESP_Kx3y3z_S_aa;
    Double I_ESP_Kx2y4z_Px_aa = I_ESP_L2x2y4z_S_aa+ABX*I_ESP_Kx2y4z_S_aa;
    Double I_ESP_Kxy5z_Px_aa = I_ESP_L2xy5z_S_aa+ABX*I_ESP_Kxy5z_S_aa;
    Double I_ESP_Kx6z_Px_aa = I_ESP_L2x6z_S_aa+ABX*I_ESP_Kx6z_S_aa;
    Double I_ESP_K7y_Px_aa = I_ESP_Lx7y_S_aa+ABX*I_ESP_K7y_S_aa;
    Double I_ESP_K6yz_Px_aa = I_ESP_Lx6yz_S_aa+ABX*I_ESP_K6yz_S_aa;
    Double I_ESP_K5y2z_Px_aa = I_ESP_Lx5y2z_S_aa+ABX*I_ESP_K5y2z_S_aa;
    Double I_ESP_K4y3z_Px_aa = I_ESP_Lx4y3z_S_aa+ABX*I_ESP_K4y3z_S_aa;
    Double I_ESP_K3y4z_Px_aa = I_ESP_Lx3y4z_S_aa+ABX*I_ESP_K3y4z_S_aa;
    Double I_ESP_K2y5z_Px_aa = I_ESP_Lx2y5z_S_aa+ABX*I_ESP_K2y5z_S_aa;
    Double I_ESP_Ky6z_Px_aa = I_ESP_Lxy6z_S_aa+ABX*I_ESP_Ky6z_S_aa;
    Double I_ESP_K7z_Px_aa = I_ESP_Lx7z_S_aa+ABX*I_ESP_K7z_S_aa;
    Double I_ESP_K7x_Py_aa = I_ESP_L7xy_S_aa+ABY*I_ESP_K7x_S_aa;
    Double I_ESP_K6xy_Py_aa = I_ESP_L6x2y_S_aa+ABY*I_ESP_K6xy_S_aa;
    Double I_ESP_K6xz_Py_aa = I_ESP_L6xyz_S_aa+ABY*I_ESP_K6xz_S_aa;
    Double I_ESP_K5x2y_Py_aa = I_ESP_L5x3y_S_aa+ABY*I_ESP_K5x2y_S_aa;
    Double I_ESP_K5xyz_Py_aa = I_ESP_L5x2yz_S_aa+ABY*I_ESP_K5xyz_S_aa;
    Double I_ESP_K5x2z_Py_aa = I_ESP_L5xy2z_S_aa+ABY*I_ESP_K5x2z_S_aa;
    Double I_ESP_K4x3y_Py_aa = I_ESP_L4x4y_S_aa+ABY*I_ESP_K4x3y_S_aa;
    Double I_ESP_K4x2yz_Py_aa = I_ESP_L4x3yz_S_aa+ABY*I_ESP_K4x2yz_S_aa;
    Double I_ESP_K4xy2z_Py_aa = I_ESP_L4x2y2z_S_aa+ABY*I_ESP_K4xy2z_S_aa;
    Double I_ESP_K4x3z_Py_aa = I_ESP_L4xy3z_S_aa+ABY*I_ESP_K4x3z_S_aa;
    Double I_ESP_K3x4y_Py_aa = I_ESP_L3x5y_S_aa+ABY*I_ESP_K3x4y_S_aa;
    Double I_ESP_K3x3yz_Py_aa = I_ESP_L3x4yz_S_aa+ABY*I_ESP_K3x3yz_S_aa;
    Double I_ESP_K3x2y2z_Py_aa = I_ESP_L3x3y2z_S_aa+ABY*I_ESP_K3x2y2z_S_aa;
    Double I_ESP_K3xy3z_Py_aa = I_ESP_L3x2y3z_S_aa+ABY*I_ESP_K3xy3z_S_aa;
    Double I_ESP_K3x4z_Py_aa = I_ESP_L3xy4z_S_aa+ABY*I_ESP_K3x4z_S_aa;
    Double I_ESP_K2x5y_Py_aa = I_ESP_L2x6y_S_aa+ABY*I_ESP_K2x5y_S_aa;
    Double I_ESP_K2x4yz_Py_aa = I_ESP_L2x5yz_S_aa+ABY*I_ESP_K2x4yz_S_aa;
    Double I_ESP_K2x3y2z_Py_aa = I_ESP_L2x4y2z_S_aa+ABY*I_ESP_K2x3y2z_S_aa;
    Double I_ESP_K2x2y3z_Py_aa = I_ESP_L2x3y3z_S_aa+ABY*I_ESP_K2x2y3z_S_aa;
    Double I_ESP_K2xy4z_Py_aa = I_ESP_L2x2y4z_S_aa+ABY*I_ESP_K2xy4z_S_aa;
    Double I_ESP_K2x5z_Py_aa = I_ESP_L2xy5z_S_aa+ABY*I_ESP_K2x5z_S_aa;
    Double I_ESP_Kx6y_Py_aa = I_ESP_Lx7y_S_aa+ABY*I_ESP_Kx6y_S_aa;
    Double I_ESP_Kx5yz_Py_aa = I_ESP_Lx6yz_S_aa+ABY*I_ESP_Kx5yz_S_aa;
    Double I_ESP_Kx4y2z_Py_aa = I_ESP_Lx5y2z_S_aa+ABY*I_ESP_Kx4y2z_S_aa;
    Double I_ESP_Kx3y3z_Py_aa = I_ESP_Lx4y3z_S_aa+ABY*I_ESP_Kx3y3z_S_aa;
    Double I_ESP_Kx2y4z_Py_aa = I_ESP_Lx3y4z_S_aa+ABY*I_ESP_Kx2y4z_S_aa;
    Double I_ESP_Kxy5z_Py_aa = I_ESP_Lx2y5z_S_aa+ABY*I_ESP_Kxy5z_S_aa;
    Double I_ESP_Kx6z_Py_aa = I_ESP_Lxy6z_S_aa+ABY*I_ESP_Kx6z_S_aa;
    Double I_ESP_K7y_Py_aa = I_ESP_L8y_S_aa+ABY*I_ESP_K7y_S_aa;
    Double I_ESP_K6yz_Py_aa = I_ESP_L7yz_S_aa+ABY*I_ESP_K6yz_S_aa;
    Double I_ESP_K5y2z_Py_aa = I_ESP_L6y2z_S_aa+ABY*I_ESP_K5y2z_S_aa;
    Double I_ESP_K4y3z_Py_aa = I_ESP_L5y3z_S_aa+ABY*I_ESP_K4y3z_S_aa;
    Double I_ESP_K3y4z_Py_aa = I_ESP_L4y4z_S_aa+ABY*I_ESP_K3y4z_S_aa;
    Double I_ESP_K2y5z_Py_aa = I_ESP_L3y5z_S_aa+ABY*I_ESP_K2y5z_S_aa;
    Double I_ESP_Ky6z_Py_aa = I_ESP_L2y6z_S_aa+ABY*I_ESP_Ky6z_S_aa;
    Double I_ESP_K7z_Py_aa = I_ESP_Ly7z_S_aa+ABY*I_ESP_K7z_S_aa;
    Double I_ESP_K7x_Pz_aa = I_ESP_L7xz_S_aa+ABZ*I_ESP_K7x_S_aa;
    Double I_ESP_K6xy_Pz_aa = I_ESP_L6xyz_S_aa+ABZ*I_ESP_K6xy_S_aa;
    Double I_ESP_K6xz_Pz_aa = I_ESP_L6x2z_S_aa+ABZ*I_ESP_K6xz_S_aa;
    Double I_ESP_K5x2y_Pz_aa = I_ESP_L5x2yz_S_aa+ABZ*I_ESP_K5x2y_S_aa;
    Double I_ESP_K5xyz_Pz_aa = I_ESP_L5xy2z_S_aa+ABZ*I_ESP_K5xyz_S_aa;
    Double I_ESP_K5x2z_Pz_aa = I_ESP_L5x3z_S_aa+ABZ*I_ESP_K5x2z_S_aa;
    Double I_ESP_K4x3y_Pz_aa = I_ESP_L4x3yz_S_aa+ABZ*I_ESP_K4x3y_S_aa;
    Double I_ESP_K4x2yz_Pz_aa = I_ESP_L4x2y2z_S_aa+ABZ*I_ESP_K4x2yz_S_aa;
    Double I_ESP_K4xy2z_Pz_aa = I_ESP_L4xy3z_S_aa+ABZ*I_ESP_K4xy2z_S_aa;
    Double I_ESP_K4x3z_Pz_aa = I_ESP_L4x4z_S_aa+ABZ*I_ESP_K4x3z_S_aa;
    Double I_ESP_K3x4y_Pz_aa = I_ESP_L3x4yz_S_aa+ABZ*I_ESP_K3x4y_S_aa;
    Double I_ESP_K3x3yz_Pz_aa = I_ESP_L3x3y2z_S_aa+ABZ*I_ESP_K3x3yz_S_aa;
    Double I_ESP_K3x2y2z_Pz_aa = I_ESP_L3x2y3z_S_aa+ABZ*I_ESP_K3x2y2z_S_aa;
    Double I_ESP_K3xy3z_Pz_aa = I_ESP_L3xy4z_S_aa+ABZ*I_ESP_K3xy3z_S_aa;
    Double I_ESP_K3x4z_Pz_aa = I_ESP_L3x5z_S_aa+ABZ*I_ESP_K3x4z_S_aa;
    Double I_ESP_K2x5y_Pz_aa = I_ESP_L2x5yz_S_aa+ABZ*I_ESP_K2x5y_S_aa;
    Double I_ESP_K2x4yz_Pz_aa = I_ESP_L2x4y2z_S_aa+ABZ*I_ESP_K2x4yz_S_aa;
    Double I_ESP_K2x3y2z_Pz_aa = I_ESP_L2x3y3z_S_aa+ABZ*I_ESP_K2x3y2z_S_aa;
    Double I_ESP_K2x2y3z_Pz_aa = I_ESP_L2x2y4z_S_aa+ABZ*I_ESP_K2x2y3z_S_aa;
    Double I_ESP_K2xy4z_Pz_aa = I_ESP_L2xy5z_S_aa+ABZ*I_ESP_K2xy4z_S_aa;
    Double I_ESP_K2x5z_Pz_aa = I_ESP_L2x6z_S_aa+ABZ*I_ESP_K2x5z_S_aa;
    Double I_ESP_Kx6y_Pz_aa = I_ESP_Lx6yz_S_aa+ABZ*I_ESP_Kx6y_S_aa;
    Double I_ESP_Kx5yz_Pz_aa = I_ESP_Lx5y2z_S_aa+ABZ*I_ESP_Kx5yz_S_aa;
    Double I_ESP_Kx4y2z_Pz_aa = I_ESP_Lx4y3z_S_aa+ABZ*I_ESP_Kx4y2z_S_aa;
    Double I_ESP_Kx3y3z_Pz_aa = I_ESP_Lx3y4z_S_aa+ABZ*I_ESP_Kx3y3z_S_aa;
    Double I_ESP_Kx2y4z_Pz_aa = I_ESP_Lx2y5z_S_aa+ABZ*I_ESP_Kx2y4z_S_aa;
    Double I_ESP_Kxy5z_Pz_aa = I_ESP_Lxy6z_S_aa+ABZ*I_ESP_Kxy5z_S_aa;
    Double I_ESP_Kx6z_Pz_aa = I_ESP_Lx7z_S_aa+ABZ*I_ESP_Kx6z_S_aa;
    Double I_ESP_K7y_Pz_aa = I_ESP_L7yz_S_aa+ABZ*I_ESP_K7y_S_aa;
    Double I_ESP_K6yz_Pz_aa = I_ESP_L6y2z_S_aa+ABZ*I_ESP_K6yz_S_aa;
    Double I_ESP_K5y2z_Pz_aa = I_ESP_L5y3z_S_aa+ABZ*I_ESP_K5y2z_S_aa;
    Double I_ESP_K4y3z_Pz_aa = I_ESP_L4y4z_S_aa+ABZ*I_ESP_K4y3z_S_aa;
    Double I_ESP_K3y4z_Pz_aa = I_ESP_L3y5z_S_aa+ABZ*I_ESP_K3y4z_S_aa;
    Double I_ESP_K2y5z_Pz_aa = I_ESP_L2y6z_S_aa+ABZ*I_ESP_K2y5z_S_aa;
    Double I_ESP_Ky6z_Pz_aa = I_ESP_Ly7z_S_aa+ABZ*I_ESP_Ky6z_S_aa;
    Double I_ESP_K7z_Pz_aa = I_ESP_L8z_S_aa+ABZ*I_ESP_K7z_S_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_L_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 10 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_M_S_aa
     * RHS shell quartet name: SQ_ESP_L_S_aa
     ************************************************************/
    Double I_ESP_L8x_Px_aa = I_ESP_M9x_S_aa+ABX*I_ESP_L8x_S_aa;
    Double I_ESP_L7xy_Px_aa = I_ESP_M8xy_S_aa+ABX*I_ESP_L7xy_S_aa;
    Double I_ESP_L7xz_Px_aa = I_ESP_M8xz_S_aa+ABX*I_ESP_L7xz_S_aa;
    Double I_ESP_L6x2y_Px_aa = I_ESP_M7x2y_S_aa+ABX*I_ESP_L6x2y_S_aa;
    Double I_ESP_L6xyz_Px_aa = I_ESP_M7xyz_S_aa+ABX*I_ESP_L6xyz_S_aa;
    Double I_ESP_L6x2z_Px_aa = I_ESP_M7x2z_S_aa+ABX*I_ESP_L6x2z_S_aa;
    Double I_ESP_L5x3y_Px_aa = I_ESP_M6x3y_S_aa+ABX*I_ESP_L5x3y_S_aa;
    Double I_ESP_L5x2yz_Px_aa = I_ESP_M6x2yz_S_aa+ABX*I_ESP_L5x2yz_S_aa;
    Double I_ESP_L5xy2z_Px_aa = I_ESP_M6xy2z_S_aa+ABX*I_ESP_L5xy2z_S_aa;
    Double I_ESP_L5x3z_Px_aa = I_ESP_M6x3z_S_aa+ABX*I_ESP_L5x3z_S_aa;
    Double I_ESP_L4x4y_Px_aa = I_ESP_M5x4y_S_aa+ABX*I_ESP_L4x4y_S_aa;
    Double I_ESP_L4x3yz_Px_aa = I_ESP_M5x3yz_S_aa+ABX*I_ESP_L4x3yz_S_aa;
    Double I_ESP_L4x2y2z_Px_aa = I_ESP_M5x2y2z_S_aa+ABX*I_ESP_L4x2y2z_S_aa;
    Double I_ESP_L4xy3z_Px_aa = I_ESP_M5xy3z_S_aa+ABX*I_ESP_L4xy3z_S_aa;
    Double I_ESP_L4x4z_Px_aa = I_ESP_M5x4z_S_aa+ABX*I_ESP_L4x4z_S_aa;
    Double I_ESP_L3x5y_Px_aa = I_ESP_M4x5y_S_aa+ABX*I_ESP_L3x5y_S_aa;
    Double I_ESP_L3x4yz_Px_aa = I_ESP_M4x4yz_S_aa+ABX*I_ESP_L3x4yz_S_aa;
    Double I_ESP_L3x3y2z_Px_aa = I_ESP_M4x3y2z_S_aa+ABX*I_ESP_L3x3y2z_S_aa;
    Double I_ESP_L3x2y3z_Px_aa = I_ESP_M4x2y3z_S_aa+ABX*I_ESP_L3x2y3z_S_aa;
    Double I_ESP_L3xy4z_Px_aa = I_ESP_M4xy4z_S_aa+ABX*I_ESP_L3xy4z_S_aa;
    Double I_ESP_L3x5z_Px_aa = I_ESP_M4x5z_S_aa+ABX*I_ESP_L3x5z_S_aa;
    Double I_ESP_L2x6y_Px_aa = I_ESP_M3x6y_S_aa+ABX*I_ESP_L2x6y_S_aa;
    Double I_ESP_L2x5yz_Px_aa = I_ESP_M3x5yz_S_aa+ABX*I_ESP_L2x5yz_S_aa;
    Double I_ESP_L2x4y2z_Px_aa = I_ESP_M3x4y2z_S_aa+ABX*I_ESP_L2x4y2z_S_aa;
    Double I_ESP_L2x3y3z_Px_aa = I_ESP_M3x3y3z_S_aa+ABX*I_ESP_L2x3y3z_S_aa;
    Double I_ESP_L2x2y4z_Px_aa = I_ESP_M3x2y4z_S_aa+ABX*I_ESP_L2x2y4z_S_aa;
    Double I_ESP_L2xy5z_Px_aa = I_ESP_M3xy5z_S_aa+ABX*I_ESP_L2xy5z_S_aa;
    Double I_ESP_L2x6z_Px_aa = I_ESP_M3x6z_S_aa+ABX*I_ESP_L2x6z_S_aa;
    Double I_ESP_Lx7y_Px_aa = I_ESP_M2x7y_S_aa+ABX*I_ESP_Lx7y_S_aa;
    Double I_ESP_Lx6yz_Px_aa = I_ESP_M2x6yz_S_aa+ABX*I_ESP_Lx6yz_S_aa;
    Double I_ESP_Lx5y2z_Px_aa = I_ESP_M2x5y2z_S_aa+ABX*I_ESP_Lx5y2z_S_aa;
    Double I_ESP_Lx4y3z_Px_aa = I_ESP_M2x4y3z_S_aa+ABX*I_ESP_Lx4y3z_S_aa;
    Double I_ESP_Lx3y4z_Px_aa = I_ESP_M2x3y4z_S_aa+ABX*I_ESP_Lx3y4z_S_aa;
    Double I_ESP_Lx2y5z_Px_aa = I_ESP_M2x2y5z_S_aa+ABX*I_ESP_Lx2y5z_S_aa;
    Double I_ESP_Lxy6z_Px_aa = I_ESP_M2xy6z_S_aa+ABX*I_ESP_Lxy6z_S_aa;
    Double I_ESP_Lx7z_Px_aa = I_ESP_M2x7z_S_aa+ABX*I_ESP_Lx7z_S_aa;
    Double I_ESP_L8y_Px_aa = I_ESP_Mx8y_S_aa+ABX*I_ESP_L8y_S_aa;
    Double I_ESP_L7yz_Px_aa = I_ESP_Mx7yz_S_aa+ABX*I_ESP_L7yz_S_aa;
    Double I_ESP_L6y2z_Px_aa = I_ESP_Mx6y2z_S_aa+ABX*I_ESP_L6y2z_S_aa;
    Double I_ESP_L5y3z_Px_aa = I_ESP_Mx5y3z_S_aa+ABX*I_ESP_L5y3z_S_aa;
    Double I_ESP_L4y4z_Px_aa = I_ESP_Mx4y4z_S_aa+ABX*I_ESP_L4y4z_S_aa;
    Double I_ESP_L3y5z_Px_aa = I_ESP_Mx3y5z_S_aa+ABX*I_ESP_L3y5z_S_aa;
    Double I_ESP_L2y6z_Px_aa = I_ESP_Mx2y6z_S_aa+ABX*I_ESP_L2y6z_S_aa;
    Double I_ESP_Ly7z_Px_aa = I_ESP_Mxy7z_S_aa+ABX*I_ESP_Ly7z_S_aa;
    Double I_ESP_L8z_Px_aa = I_ESP_Mx8z_S_aa+ABX*I_ESP_L8z_S_aa;
    Double I_ESP_L7xy_Py_aa = I_ESP_M7x2y_S_aa+ABY*I_ESP_L7xy_S_aa;
    Double I_ESP_L7xz_Py_aa = I_ESP_M7xyz_S_aa+ABY*I_ESP_L7xz_S_aa;
    Double I_ESP_L6x2y_Py_aa = I_ESP_M6x3y_S_aa+ABY*I_ESP_L6x2y_S_aa;
    Double I_ESP_L6xyz_Py_aa = I_ESP_M6x2yz_S_aa+ABY*I_ESP_L6xyz_S_aa;
    Double I_ESP_L6x2z_Py_aa = I_ESP_M6xy2z_S_aa+ABY*I_ESP_L6x2z_S_aa;
    Double I_ESP_L5x3y_Py_aa = I_ESP_M5x4y_S_aa+ABY*I_ESP_L5x3y_S_aa;
    Double I_ESP_L5x2yz_Py_aa = I_ESP_M5x3yz_S_aa+ABY*I_ESP_L5x2yz_S_aa;
    Double I_ESP_L5xy2z_Py_aa = I_ESP_M5x2y2z_S_aa+ABY*I_ESP_L5xy2z_S_aa;
    Double I_ESP_L5x3z_Py_aa = I_ESP_M5xy3z_S_aa+ABY*I_ESP_L5x3z_S_aa;
    Double I_ESP_L4x4y_Py_aa = I_ESP_M4x5y_S_aa+ABY*I_ESP_L4x4y_S_aa;
    Double I_ESP_L4x3yz_Py_aa = I_ESP_M4x4yz_S_aa+ABY*I_ESP_L4x3yz_S_aa;
    Double I_ESP_L4x2y2z_Py_aa = I_ESP_M4x3y2z_S_aa+ABY*I_ESP_L4x2y2z_S_aa;
    Double I_ESP_L4xy3z_Py_aa = I_ESP_M4x2y3z_S_aa+ABY*I_ESP_L4xy3z_S_aa;
    Double I_ESP_L4x4z_Py_aa = I_ESP_M4xy4z_S_aa+ABY*I_ESP_L4x4z_S_aa;
    Double I_ESP_L3x5y_Py_aa = I_ESP_M3x6y_S_aa+ABY*I_ESP_L3x5y_S_aa;
    Double I_ESP_L3x4yz_Py_aa = I_ESP_M3x5yz_S_aa+ABY*I_ESP_L3x4yz_S_aa;
    Double I_ESP_L3x3y2z_Py_aa = I_ESP_M3x4y2z_S_aa+ABY*I_ESP_L3x3y2z_S_aa;
    Double I_ESP_L3x2y3z_Py_aa = I_ESP_M3x3y3z_S_aa+ABY*I_ESP_L3x2y3z_S_aa;
    Double I_ESP_L3xy4z_Py_aa = I_ESP_M3x2y4z_S_aa+ABY*I_ESP_L3xy4z_S_aa;
    Double I_ESP_L3x5z_Py_aa = I_ESP_M3xy5z_S_aa+ABY*I_ESP_L3x5z_S_aa;
    Double I_ESP_L2x6y_Py_aa = I_ESP_M2x7y_S_aa+ABY*I_ESP_L2x6y_S_aa;
    Double I_ESP_L2x5yz_Py_aa = I_ESP_M2x6yz_S_aa+ABY*I_ESP_L2x5yz_S_aa;
    Double I_ESP_L2x4y2z_Py_aa = I_ESP_M2x5y2z_S_aa+ABY*I_ESP_L2x4y2z_S_aa;
    Double I_ESP_L2x3y3z_Py_aa = I_ESP_M2x4y3z_S_aa+ABY*I_ESP_L2x3y3z_S_aa;
    Double I_ESP_L2x2y4z_Py_aa = I_ESP_M2x3y4z_S_aa+ABY*I_ESP_L2x2y4z_S_aa;
    Double I_ESP_L2xy5z_Py_aa = I_ESP_M2x2y5z_S_aa+ABY*I_ESP_L2xy5z_S_aa;
    Double I_ESP_L2x6z_Py_aa = I_ESP_M2xy6z_S_aa+ABY*I_ESP_L2x6z_S_aa;
    Double I_ESP_Lx7y_Py_aa = I_ESP_Mx8y_S_aa+ABY*I_ESP_Lx7y_S_aa;
    Double I_ESP_Lx6yz_Py_aa = I_ESP_Mx7yz_S_aa+ABY*I_ESP_Lx6yz_S_aa;
    Double I_ESP_Lx5y2z_Py_aa = I_ESP_Mx6y2z_S_aa+ABY*I_ESP_Lx5y2z_S_aa;
    Double I_ESP_Lx4y3z_Py_aa = I_ESP_Mx5y3z_S_aa+ABY*I_ESP_Lx4y3z_S_aa;
    Double I_ESP_Lx3y4z_Py_aa = I_ESP_Mx4y4z_S_aa+ABY*I_ESP_Lx3y4z_S_aa;
    Double I_ESP_Lx2y5z_Py_aa = I_ESP_Mx3y5z_S_aa+ABY*I_ESP_Lx2y5z_S_aa;
    Double I_ESP_Lxy6z_Py_aa = I_ESP_Mx2y6z_S_aa+ABY*I_ESP_Lxy6z_S_aa;
    Double I_ESP_Lx7z_Py_aa = I_ESP_Mxy7z_S_aa+ABY*I_ESP_Lx7z_S_aa;
    Double I_ESP_L8y_Py_aa = I_ESP_M9y_S_aa+ABY*I_ESP_L8y_S_aa;
    Double I_ESP_L7yz_Py_aa = I_ESP_M8yz_S_aa+ABY*I_ESP_L7yz_S_aa;
    Double I_ESP_L6y2z_Py_aa = I_ESP_M7y2z_S_aa+ABY*I_ESP_L6y2z_S_aa;
    Double I_ESP_L5y3z_Py_aa = I_ESP_M6y3z_S_aa+ABY*I_ESP_L5y3z_S_aa;
    Double I_ESP_L4y4z_Py_aa = I_ESP_M5y4z_S_aa+ABY*I_ESP_L4y4z_S_aa;
    Double I_ESP_L3y5z_Py_aa = I_ESP_M4y5z_S_aa+ABY*I_ESP_L3y5z_S_aa;
    Double I_ESP_L2y6z_Py_aa = I_ESP_M3y6z_S_aa+ABY*I_ESP_L2y6z_S_aa;
    Double I_ESP_Ly7z_Py_aa = I_ESP_M2y7z_S_aa+ABY*I_ESP_Ly7z_S_aa;
    Double I_ESP_L8z_Py_aa = I_ESP_My8z_S_aa+ABY*I_ESP_L8z_S_aa;
    Double I_ESP_L7xz_Pz_aa = I_ESP_M7x2z_S_aa+ABZ*I_ESP_L7xz_S_aa;
    Double I_ESP_L6xyz_Pz_aa = I_ESP_M6xy2z_S_aa+ABZ*I_ESP_L6xyz_S_aa;
    Double I_ESP_L6x2z_Pz_aa = I_ESP_M6x3z_S_aa+ABZ*I_ESP_L6x2z_S_aa;
    Double I_ESP_L5x2yz_Pz_aa = I_ESP_M5x2y2z_S_aa+ABZ*I_ESP_L5x2yz_S_aa;
    Double I_ESP_L5xy2z_Pz_aa = I_ESP_M5xy3z_S_aa+ABZ*I_ESP_L5xy2z_S_aa;
    Double I_ESP_L5x3z_Pz_aa = I_ESP_M5x4z_S_aa+ABZ*I_ESP_L5x3z_S_aa;
    Double I_ESP_L4x3yz_Pz_aa = I_ESP_M4x3y2z_S_aa+ABZ*I_ESP_L4x3yz_S_aa;
    Double I_ESP_L4x2y2z_Pz_aa = I_ESP_M4x2y3z_S_aa+ABZ*I_ESP_L4x2y2z_S_aa;
    Double I_ESP_L4xy3z_Pz_aa = I_ESP_M4xy4z_S_aa+ABZ*I_ESP_L4xy3z_S_aa;
    Double I_ESP_L4x4z_Pz_aa = I_ESP_M4x5z_S_aa+ABZ*I_ESP_L4x4z_S_aa;
    Double I_ESP_L3x4yz_Pz_aa = I_ESP_M3x4y2z_S_aa+ABZ*I_ESP_L3x4yz_S_aa;
    Double I_ESP_L3x3y2z_Pz_aa = I_ESP_M3x3y3z_S_aa+ABZ*I_ESP_L3x3y2z_S_aa;
    Double I_ESP_L3x2y3z_Pz_aa = I_ESP_M3x2y4z_S_aa+ABZ*I_ESP_L3x2y3z_S_aa;
    Double I_ESP_L3xy4z_Pz_aa = I_ESP_M3xy5z_S_aa+ABZ*I_ESP_L3xy4z_S_aa;
    Double I_ESP_L3x5z_Pz_aa = I_ESP_M3x6z_S_aa+ABZ*I_ESP_L3x5z_S_aa;
    Double I_ESP_L2x5yz_Pz_aa = I_ESP_M2x5y2z_S_aa+ABZ*I_ESP_L2x5yz_S_aa;
    Double I_ESP_L2x4y2z_Pz_aa = I_ESP_M2x4y3z_S_aa+ABZ*I_ESP_L2x4y2z_S_aa;
    Double I_ESP_L2x3y3z_Pz_aa = I_ESP_M2x3y4z_S_aa+ABZ*I_ESP_L2x3y3z_S_aa;
    Double I_ESP_L2x2y4z_Pz_aa = I_ESP_M2x2y5z_S_aa+ABZ*I_ESP_L2x2y4z_S_aa;
    Double I_ESP_L2xy5z_Pz_aa = I_ESP_M2xy6z_S_aa+ABZ*I_ESP_L2xy5z_S_aa;
    Double I_ESP_L2x6z_Pz_aa = I_ESP_M2x7z_S_aa+ABZ*I_ESP_L2x6z_S_aa;
    Double I_ESP_Lx6yz_Pz_aa = I_ESP_Mx6y2z_S_aa+ABZ*I_ESP_Lx6yz_S_aa;
    Double I_ESP_Lx5y2z_Pz_aa = I_ESP_Mx5y3z_S_aa+ABZ*I_ESP_Lx5y2z_S_aa;
    Double I_ESP_Lx4y3z_Pz_aa = I_ESP_Mx4y4z_S_aa+ABZ*I_ESP_Lx4y3z_S_aa;
    Double I_ESP_Lx3y4z_Pz_aa = I_ESP_Mx3y5z_S_aa+ABZ*I_ESP_Lx3y4z_S_aa;
    Double I_ESP_Lx2y5z_Pz_aa = I_ESP_Mx2y6z_S_aa+ABZ*I_ESP_Lx2y5z_S_aa;
    Double I_ESP_Lxy6z_Pz_aa = I_ESP_Mxy7z_S_aa+ABZ*I_ESP_Lxy6z_S_aa;
    Double I_ESP_Lx7z_Pz_aa = I_ESP_Mx8z_S_aa+ABZ*I_ESP_Lx7z_S_aa;
    Double I_ESP_L7yz_Pz_aa = I_ESP_M7y2z_S_aa+ABZ*I_ESP_L7yz_S_aa;
    Double I_ESP_L6y2z_Pz_aa = I_ESP_M6y3z_S_aa+ABZ*I_ESP_L6y2z_S_aa;
    Double I_ESP_L5y3z_Pz_aa = I_ESP_M5y4z_S_aa+ABZ*I_ESP_L5y3z_S_aa;
    Double I_ESP_L4y4z_Pz_aa = I_ESP_M4y5z_S_aa+ABZ*I_ESP_L4y4z_S_aa;
    Double I_ESP_L3y5z_Pz_aa = I_ESP_M3y6z_S_aa+ABZ*I_ESP_L3y5z_S_aa;
    Double I_ESP_L2y6z_Pz_aa = I_ESP_M2y7z_S_aa+ABZ*I_ESP_L2y6z_S_aa;
    Double I_ESP_Ly7z_Pz_aa = I_ESP_My8z_S_aa+ABZ*I_ESP_Ly7z_S_aa;
    Double I_ESP_L8z_Pz_aa = I_ESP_M9z_S_aa+ABZ*I_ESP_L8z_S_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_K_D_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_P_aa
     * RHS shell quartet name: SQ_ESP_K_P_aa
     ************************************************************/
    Double I_ESP_K7x_D2x_aa = I_ESP_L8x_Px_aa+ABX*I_ESP_K7x_Px_aa;
    Double I_ESP_K6xy_D2x_aa = I_ESP_L7xy_Px_aa+ABX*I_ESP_K6xy_Px_aa;
    Double I_ESP_K6xz_D2x_aa = I_ESP_L7xz_Px_aa+ABX*I_ESP_K6xz_Px_aa;
    Double I_ESP_K5x2y_D2x_aa = I_ESP_L6x2y_Px_aa+ABX*I_ESP_K5x2y_Px_aa;
    Double I_ESP_K5xyz_D2x_aa = I_ESP_L6xyz_Px_aa+ABX*I_ESP_K5xyz_Px_aa;
    Double I_ESP_K5x2z_D2x_aa = I_ESP_L6x2z_Px_aa+ABX*I_ESP_K5x2z_Px_aa;
    Double I_ESP_K4x3y_D2x_aa = I_ESP_L5x3y_Px_aa+ABX*I_ESP_K4x3y_Px_aa;
    Double I_ESP_K4x2yz_D2x_aa = I_ESP_L5x2yz_Px_aa+ABX*I_ESP_K4x2yz_Px_aa;
    Double I_ESP_K4xy2z_D2x_aa = I_ESP_L5xy2z_Px_aa+ABX*I_ESP_K4xy2z_Px_aa;
    Double I_ESP_K4x3z_D2x_aa = I_ESP_L5x3z_Px_aa+ABX*I_ESP_K4x3z_Px_aa;
    Double I_ESP_K3x4y_D2x_aa = I_ESP_L4x4y_Px_aa+ABX*I_ESP_K3x4y_Px_aa;
    Double I_ESP_K3x3yz_D2x_aa = I_ESP_L4x3yz_Px_aa+ABX*I_ESP_K3x3yz_Px_aa;
    Double I_ESP_K3x2y2z_D2x_aa = I_ESP_L4x2y2z_Px_aa+ABX*I_ESP_K3x2y2z_Px_aa;
    Double I_ESP_K3xy3z_D2x_aa = I_ESP_L4xy3z_Px_aa+ABX*I_ESP_K3xy3z_Px_aa;
    Double I_ESP_K3x4z_D2x_aa = I_ESP_L4x4z_Px_aa+ABX*I_ESP_K3x4z_Px_aa;
    Double I_ESP_K2x5y_D2x_aa = I_ESP_L3x5y_Px_aa+ABX*I_ESP_K2x5y_Px_aa;
    Double I_ESP_K2x4yz_D2x_aa = I_ESP_L3x4yz_Px_aa+ABX*I_ESP_K2x4yz_Px_aa;
    Double I_ESP_K2x3y2z_D2x_aa = I_ESP_L3x3y2z_Px_aa+ABX*I_ESP_K2x3y2z_Px_aa;
    Double I_ESP_K2x2y3z_D2x_aa = I_ESP_L3x2y3z_Px_aa+ABX*I_ESP_K2x2y3z_Px_aa;
    Double I_ESP_K2xy4z_D2x_aa = I_ESP_L3xy4z_Px_aa+ABX*I_ESP_K2xy4z_Px_aa;
    Double I_ESP_K2x5z_D2x_aa = I_ESP_L3x5z_Px_aa+ABX*I_ESP_K2x5z_Px_aa;
    Double I_ESP_Kx6y_D2x_aa = I_ESP_L2x6y_Px_aa+ABX*I_ESP_Kx6y_Px_aa;
    Double I_ESP_Kx5yz_D2x_aa = I_ESP_L2x5yz_Px_aa+ABX*I_ESP_Kx5yz_Px_aa;
    Double I_ESP_Kx4y2z_D2x_aa = I_ESP_L2x4y2z_Px_aa+ABX*I_ESP_Kx4y2z_Px_aa;
    Double I_ESP_Kx3y3z_D2x_aa = I_ESP_L2x3y3z_Px_aa+ABX*I_ESP_Kx3y3z_Px_aa;
    Double I_ESP_Kx2y4z_D2x_aa = I_ESP_L2x2y4z_Px_aa+ABX*I_ESP_Kx2y4z_Px_aa;
    Double I_ESP_Kxy5z_D2x_aa = I_ESP_L2xy5z_Px_aa+ABX*I_ESP_Kxy5z_Px_aa;
    Double I_ESP_Kx6z_D2x_aa = I_ESP_L2x6z_Px_aa+ABX*I_ESP_Kx6z_Px_aa;
    Double I_ESP_K7y_D2x_aa = I_ESP_Lx7y_Px_aa+ABX*I_ESP_K7y_Px_aa;
    Double I_ESP_K6yz_D2x_aa = I_ESP_Lx6yz_Px_aa+ABX*I_ESP_K6yz_Px_aa;
    Double I_ESP_K5y2z_D2x_aa = I_ESP_Lx5y2z_Px_aa+ABX*I_ESP_K5y2z_Px_aa;
    Double I_ESP_K4y3z_D2x_aa = I_ESP_Lx4y3z_Px_aa+ABX*I_ESP_K4y3z_Px_aa;
    Double I_ESP_K3y4z_D2x_aa = I_ESP_Lx3y4z_Px_aa+ABX*I_ESP_K3y4z_Px_aa;
    Double I_ESP_K2y5z_D2x_aa = I_ESP_Lx2y5z_Px_aa+ABX*I_ESP_K2y5z_Px_aa;
    Double I_ESP_Ky6z_D2x_aa = I_ESP_Lxy6z_Px_aa+ABX*I_ESP_Ky6z_Px_aa;
    Double I_ESP_K7z_D2x_aa = I_ESP_Lx7z_Px_aa+ABX*I_ESP_K7z_Px_aa;
    Double I_ESP_K7x_Dxy_aa = I_ESP_L7xy_Px_aa+ABY*I_ESP_K7x_Px_aa;
    Double I_ESP_K6xy_Dxy_aa = I_ESP_L6x2y_Px_aa+ABY*I_ESP_K6xy_Px_aa;
    Double I_ESP_K6xz_Dxy_aa = I_ESP_L6xyz_Px_aa+ABY*I_ESP_K6xz_Px_aa;
    Double I_ESP_K5x2y_Dxy_aa = I_ESP_L5x3y_Px_aa+ABY*I_ESP_K5x2y_Px_aa;
    Double I_ESP_K5xyz_Dxy_aa = I_ESP_L5x2yz_Px_aa+ABY*I_ESP_K5xyz_Px_aa;
    Double I_ESP_K5x2z_Dxy_aa = I_ESP_L5xy2z_Px_aa+ABY*I_ESP_K5x2z_Px_aa;
    Double I_ESP_K4x3y_Dxy_aa = I_ESP_L4x4y_Px_aa+ABY*I_ESP_K4x3y_Px_aa;
    Double I_ESP_K4x2yz_Dxy_aa = I_ESP_L4x3yz_Px_aa+ABY*I_ESP_K4x2yz_Px_aa;
    Double I_ESP_K4xy2z_Dxy_aa = I_ESP_L4x2y2z_Px_aa+ABY*I_ESP_K4xy2z_Px_aa;
    Double I_ESP_K4x3z_Dxy_aa = I_ESP_L4xy3z_Px_aa+ABY*I_ESP_K4x3z_Px_aa;
    Double I_ESP_K3x4y_Dxy_aa = I_ESP_L3x5y_Px_aa+ABY*I_ESP_K3x4y_Px_aa;
    Double I_ESP_K3x3yz_Dxy_aa = I_ESP_L3x4yz_Px_aa+ABY*I_ESP_K3x3yz_Px_aa;
    Double I_ESP_K3x2y2z_Dxy_aa = I_ESP_L3x3y2z_Px_aa+ABY*I_ESP_K3x2y2z_Px_aa;
    Double I_ESP_K3xy3z_Dxy_aa = I_ESP_L3x2y3z_Px_aa+ABY*I_ESP_K3xy3z_Px_aa;
    Double I_ESP_K3x4z_Dxy_aa = I_ESP_L3xy4z_Px_aa+ABY*I_ESP_K3x4z_Px_aa;
    Double I_ESP_K2x5y_Dxy_aa = I_ESP_L2x6y_Px_aa+ABY*I_ESP_K2x5y_Px_aa;
    Double I_ESP_K2x4yz_Dxy_aa = I_ESP_L2x5yz_Px_aa+ABY*I_ESP_K2x4yz_Px_aa;
    Double I_ESP_K2x3y2z_Dxy_aa = I_ESP_L2x4y2z_Px_aa+ABY*I_ESP_K2x3y2z_Px_aa;
    Double I_ESP_K2x2y3z_Dxy_aa = I_ESP_L2x3y3z_Px_aa+ABY*I_ESP_K2x2y3z_Px_aa;
    Double I_ESP_K2xy4z_Dxy_aa = I_ESP_L2x2y4z_Px_aa+ABY*I_ESP_K2xy4z_Px_aa;
    Double I_ESP_K2x5z_Dxy_aa = I_ESP_L2xy5z_Px_aa+ABY*I_ESP_K2x5z_Px_aa;
    Double I_ESP_Kx6y_Dxy_aa = I_ESP_Lx7y_Px_aa+ABY*I_ESP_Kx6y_Px_aa;
    Double I_ESP_Kx5yz_Dxy_aa = I_ESP_Lx6yz_Px_aa+ABY*I_ESP_Kx5yz_Px_aa;
    Double I_ESP_Kx4y2z_Dxy_aa = I_ESP_Lx5y2z_Px_aa+ABY*I_ESP_Kx4y2z_Px_aa;
    Double I_ESP_Kx3y3z_Dxy_aa = I_ESP_Lx4y3z_Px_aa+ABY*I_ESP_Kx3y3z_Px_aa;
    Double I_ESP_Kx2y4z_Dxy_aa = I_ESP_Lx3y4z_Px_aa+ABY*I_ESP_Kx2y4z_Px_aa;
    Double I_ESP_Kxy5z_Dxy_aa = I_ESP_Lx2y5z_Px_aa+ABY*I_ESP_Kxy5z_Px_aa;
    Double I_ESP_Kx6z_Dxy_aa = I_ESP_Lxy6z_Px_aa+ABY*I_ESP_Kx6z_Px_aa;
    Double I_ESP_K7y_Dxy_aa = I_ESP_L8y_Px_aa+ABY*I_ESP_K7y_Px_aa;
    Double I_ESP_K6yz_Dxy_aa = I_ESP_L7yz_Px_aa+ABY*I_ESP_K6yz_Px_aa;
    Double I_ESP_K5y2z_Dxy_aa = I_ESP_L6y2z_Px_aa+ABY*I_ESP_K5y2z_Px_aa;
    Double I_ESP_K4y3z_Dxy_aa = I_ESP_L5y3z_Px_aa+ABY*I_ESP_K4y3z_Px_aa;
    Double I_ESP_K3y4z_Dxy_aa = I_ESP_L4y4z_Px_aa+ABY*I_ESP_K3y4z_Px_aa;
    Double I_ESP_K2y5z_Dxy_aa = I_ESP_L3y5z_Px_aa+ABY*I_ESP_K2y5z_Px_aa;
    Double I_ESP_Ky6z_Dxy_aa = I_ESP_L2y6z_Px_aa+ABY*I_ESP_Ky6z_Px_aa;
    Double I_ESP_K7z_Dxy_aa = I_ESP_Ly7z_Px_aa+ABY*I_ESP_K7z_Px_aa;
    Double I_ESP_K7x_Dxz_aa = I_ESP_L7xz_Px_aa+ABZ*I_ESP_K7x_Px_aa;
    Double I_ESP_K6xy_Dxz_aa = I_ESP_L6xyz_Px_aa+ABZ*I_ESP_K6xy_Px_aa;
    Double I_ESP_K6xz_Dxz_aa = I_ESP_L6x2z_Px_aa+ABZ*I_ESP_K6xz_Px_aa;
    Double I_ESP_K5x2y_Dxz_aa = I_ESP_L5x2yz_Px_aa+ABZ*I_ESP_K5x2y_Px_aa;
    Double I_ESP_K5xyz_Dxz_aa = I_ESP_L5xy2z_Px_aa+ABZ*I_ESP_K5xyz_Px_aa;
    Double I_ESP_K5x2z_Dxz_aa = I_ESP_L5x3z_Px_aa+ABZ*I_ESP_K5x2z_Px_aa;
    Double I_ESP_K4x3y_Dxz_aa = I_ESP_L4x3yz_Px_aa+ABZ*I_ESP_K4x3y_Px_aa;
    Double I_ESP_K4x2yz_Dxz_aa = I_ESP_L4x2y2z_Px_aa+ABZ*I_ESP_K4x2yz_Px_aa;
    Double I_ESP_K4xy2z_Dxz_aa = I_ESP_L4xy3z_Px_aa+ABZ*I_ESP_K4xy2z_Px_aa;
    Double I_ESP_K4x3z_Dxz_aa = I_ESP_L4x4z_Px_aa+ABZ*I_ESP_K4x3z_Px_aa;
    Double I_ESP_K3x4y_Dxz_aa = I_ESP_L3x4yz_Px_aa+ABZ*I_ESP_K3x4y_Px_aa;
    Double I_ESP_K3x3yz_Dxz_aa = I_ESP_L3x3y2z_Px_aa+ABZ*I_ESP_K3x3yz_Px_aa;
    Double I_ESP_K3x2y2z_Dxz_aa = I_ESP_L3x2y3z_Px_aa+ABZ*I_ESP_K3x2y2z_Px_aa;
    Double I_ESP_K3xy3z_Dxz_aa = I_ESP_L3xy4z_Px_aa+ABZ*I_ESP_K3xy3z_Px_aa;
    Double I_ESP_K3x4z_Dxz_aa = I_ESP_L3x5z_Px_aa+ABZ*I_ESP_K3x4z_Px_aa;
    Double I_ESP_K2x5y_Dxz_aa = I_ESP_L2x5yz_Px_aa+ABZ*I_ESP_K2x5y_Px_aa;
    Double I_ESP_K2x4yz_Dxz_aa = I_ESP_L2x4y2z_Px_aa+ABZ*I_ESP_K2x4yz_Px_aa;
    Double I_ESP_K2x3y2z_Dxz_aa = I_ESP_L2x3y3z_Px_aa+ABZ*I_ESP_K2x3y2z_Px_aa;
    Double I_ESP_K2x2y3z_Dxz_aa = I_ESP_L2x2y4z_Px_aa+ABZ*I_ESP_K2x2y3z_Px_aa;
    Double I_ESP_K2xy4z_Dxz_aa = I_ESP_L2xy5z_Px_aa+ABZ*I_ESP_K2xy4z_Px_aa;
    Double I_ESP_K2x5z_Dxz_aa = I_ESP_L2x6z_Px_aa+ABZ*I_ESP_K2x5z_Px_aa;
    Double I_ESP_Kx6y_Dxz_aa = I_ESP_Lx6yz_Px_aa+ABZ*I_ESP_Kx6y_Px_aa;
    Double I_ESP_Kx5yz_Dxz_aa = I_ESP_Lx5y2z_Px_aa+ABZ*I_ESP_Kx5yz_Px_aa;
    Double I_ESP_Kx4y2z_Dxz_aa = I_ESP_Lx4y3z_Px_aa+ABZ*I_ESP_Kx4y2z_Px_aa;
    Double I_ESP_Kx3y3z_Dxz_aa = I_ESP_Lx3y4z_Px_aa+ABZ*I_ESP_Kx3y3z_Px_aa;
    Double I_ESP_Kx2y4z_Dxz_aa = I_ESP_Lx2y5z_Px_aa+ABZ*I_ESP_Kx2y4z_Px_aa;
    Double I_ESP_Kxy5z_Dxz_aa = I_ESP_Lxy6z_Px_aa+ABZ*I_ESP_Kxy5z_Px_aa;
    Double I_ESP_Kx6z_Dxz_aa = I_ESP_Lx7z_Px_aa+ABZ*I_ESP_Kx6z_Px_aa;
    Double I_ESP_K7y_Dxz_aa = I_ESP_L7yz_Px_aa+ABZ*I_ESP_K7y_Px_aa;
    Double I_ESP_K6yz_Dxz_aa = I_ESP_L6y2z_Px_aa+ABZ*I_ESP_K6yz_Px_aa;
    Double I_ESP_K5y2z_Dxz_aa = I_ESP_L5y3z_Px_aa+ABZ*I_ESP_K5y2z_Px_aa;
    Double I_ESP_K4y3z_Dxz_aa = I_ESP_L4y4z_Px_aa+ABZ*I_ESP_K4y3z_Px_aa;
    Double I_ESP_K3y4z_Dxz_aa = I_ESP_L3y5z_Px_aa+ABZ*I_ESP_K3y4z_Px_aa;
    Double I_ESP_K2y5z_Dxz_aa = I_ESP_L2y6z_Px_aa+ABZ*I_ESP_K2y5z_Px_aa;
    Double I_ESP_Ky6z_Dxz_aa = I_ESP_Ly7z_Px_aa+ABZ*I_ESP_Ky6z_Px_aa;
    Double I_ESP_K7z_Dxz_aa = I_ESP_L8z_Px_aa+ABZ*I_ESP_K7z_Px_aa;
    Double I_ESP_K7x_D2y_aa = I_ESP_L7xy_Py_aa+ABY*I_ESP_K7x_Py_aa;
    Double I_ESP_K6xy_D2y_aa = I_ESP_L6x2y_Py_aa+ABY*I_ESP_K6xy_Py_aa;
    Double I_ESP_K6xz_D2y_aa = I_ESP_L6xyz_Py_aa+ABY*I_ESP_K6xz_Py_aa;
    Double I_ESP_K5x2y_D2y_aa = I_ESP_L5x3y_Py_aa+ABY*I_ESP_K5x2y_Py_aa;
    Double I_ESP_K5xyz_D2y_aa = I_ESP_L5x2yz_Py_aa+ABY*I_ESP_K5xyz_Py_aa;
    Double I_ESP_K5x2z_D2y_aa = I_ESP_L5xy2z_Py_aa+ABY*I_ESP_K5x2z_Py_aa;
    Double I_ESP_K4x3y_D2y_aa = I_ESP_L4x4y_Py_aa+ABY*I_ESP_K4x3y_Py_aa;
    Double I_ESP_K4x2yz_D2y_aa = I_ESP_L4x3yz_Py_aa+ABY*I_ESP_K4x2yz_Py_aa;
    Double I_ESP_K4xy2z_D2y_aa = I_ESP_L4x2y2z_Py_aa+ABY*I_ESP_K4xy2z_Py_aa;
    Double I_ESP_K4x3z_D2y_aa = I_ESP_L4xy3z_Py_aa+ABY*I_ESP_K4x3z_Py_aa;
    Double I_ESP_K3x4y_D2y_aa = I_ESP_L3x5y_Py_aa+ABY*I_ESP_K3x4y_Py_aa;
    Double I_ESP_K3x3yz_D2y_aa = I_ESP_L3x4yz_Py_aa+ABY*I_ESP_K3x3yz_Py_aa;
    Double I_ESP_K3x2y2z_D2y_aa = I_ESP_L3x3y2z_Py_aa+ABY*I_ESP_K3x2y2z_Py_aa;
    Double I_ESP_K3xy3z_D2y_aa = I_ESP_L3x2y3z_Py_aa+ABY*I_ESP_K3xy3z_Py_aa;
    Double I_ESP_K3x4z_D2y_aa = I_ESP_L3xy4z_Py_aa+ABY*I_ESP_K3x4z_Py_aa;
    Double I_ESP_K2x5y_D2y_aa = I_ESP_L2x6y_Py_aa+ABY*I_ESP_K2x5y_Py_aa;
    Double I_ESP_K2x4yz_D2y_aa = I_ESP_L2x5yz_Py_aa+ABY*I_ESP_K2x4yz_Py_aa;
    Double I_ESP_K2x3y2z_D2y_aa = I_ESP_L2x4y2z_Py_aa+ABY*I_ESP_K2x3y2z_Py_aa;
    Double I_ESP_K2x2y3z_D2y_aa = I_ESP_L2x3y3z_Py_aa+ABY*I_ESP_K2x2y3z_Py_aa;
    Double I_ESP_K2xy4z_D2y_aa = I_ESP_L2x2y4z_Py_aa+ABY*I_ESP_K2xy4z_Py_aa;
    Double I_ESP_K2x5z_D2y_aa = I_ESP_L2xy5z_Py_aa+ABY*I_ESP_K2x5z_Py_aa;
    Double I_ESP_Kx6y_D2y_aa = I_ESP_Lx7y_Py_aa+ABY*I_ESP_Kx6y_Py_aa;
    Double I_ESP_Kx5yz_D2y_aa = I_ESP_Lx6yz_Py_aa+ABY*I_ESP_Kx5yz_Py_aa;
    Double I_ESP_Kx4y2z_D2y_aa = I_ESP_Lx5y2z_Py_aa+ABY*I_ESP_Kx4y2z_Py_aa;
    Double I_ESP_Kx3y3z_D2y_aa = I_ESP_Lx4y3z_Py_aa+ABY*I_ESP_Kx3y3z_Py_aa;
    Double I_ESP_Kx2y4z_D2y_aa = I_ESP_Lx3y4z_Py_aa+ABY*I_ESP_Kx2y4z_Py_aa;
    Double I_ESP_Kxy5z_D2y_aa = I_ESP_Lx2y5z_Py_aa+ABY*I_ESP_Kxy5z_Py_aa;
    Double I_ESP_Kx6z_D2y_aa = I_ESP_Lxy6z_Py_aa+ABY*I_ESP_Kx6z_Py_aa;
    Double I_ESP_K7y_D2y_aa = I_ESP_L8y_Py_aa+ABY*I_ESP_K7y_Py_aa;
    Double I_ESP_K6yz_D2y_aa = I_ESP_L7yz_Py_aa+ABY*I_ESP_K6yz_Py_aa;
    Double I_ESP_K5y2z_D2y_aa = I_ESP_L6y2z_Py_aa+ABY*I_ESP_K5y2z_Py_aa;
    Double I_ESP_K4y3z_D2y_aa = I_ESP_L5y3z_Py_aa+ABY*I_ESP_K4y3z_Py_aa;
    Double I_ESP_K3y4z_D2y_aa = I_ESP_L4y4z_Py_aa+ABY*I_ESP_K3y4z_Py_aa;
    Double I_ESP_K2y5z_D2y_aa = I_ESP_L3y5z_Py_aa+ABY*I_ESP_K2y5z_Py_aa;
    Double I_ESP_Ky6z_D2y_aa = I_ESP_L2y6z_Py_aa+ABY*I_ESP_Ky6z_Py_aa;
    Double I_ESP_K7z_D2y_aa = I_ESP_Ly7z_Py_aa+ABY*I_ESP_K7z_Py_aa;
    Double I_ESP_K7x_Dyz_aa = I_ESP_L7xz_Py_aa+ABZ*I_ESP_K7x_Py_aa;
    Double I_ESP_K6xy_Dyz_aa = I_ESP_L6xyz_Py_aa+ABZ*I_ESP_K6xy_Py_aa;
    Double I_ESP_K6xz_Dyz_aa = I_ESP_L6x2z_Py_aa+ABZ*I_ESP_K6xz_Py_aa;
    Double I_ESP_K5x2y_Dyz_aa = I_ESP_L5x2yz_Py_aa+ABZ*I_ESP_K5x2y_Py_aa;
    Double I_ESP_K5xyz_Dyz_aa = I_ESP_L5xy2z_Py_aa+ABZ*I_ESP_K5xyz_Py_aa;
    Double I_ESP_K5x2z_Dyz_aa = I_ESP_L5x3z_Py_aa+ABZ*I_ESP_K5x2z_Py_aa;
    Double I_ESP_K4x3y_Dyz_aa = I_ESP_L4x3yz_Py_aa+ABZ*I_ESP_K4x3y_Py_aa;
    Double I_ESP_K4x2yz_Dyz_aa = I_ESP_L4x2y2z_Py_aa+ABZ*I_ESP_K4x2yz_Py_aa;
    Double I_ESP_K4xy2z_Dyz_aa = I_ESP_L4xy3z_Py_aa+ABZ*I_ESP_K4xy2z_Py_aa;
    Double I_ESP_K4x3z_Dyz_aa = I_ESP_L4x4z_Py_aa+ABZ*I_ESP_K4x3z_Py_aa;
    Double I_ESP_K3x4y_Dyz_aa = I_ESP_L3x4yz_Py_aa+ABZ*I_ESP_K3x4y_Py_aa;
    Double I_ESP_K3x3yz_Dyz_aa = I_ESP_L3x3y2z_Py_aa+ABZ*I_ESP_K3x3yz_Py_aa;
    Double I_ESP_K3x2y2z_Dyz_aa = I_ESP_L3x2y3z_Py_aa+ABZ*I_ESP_K3x2y2z_Py_aa;
    Double I_ESP_K3xy3z_Dyz_aa = I_ESP_L3xy4z_Py_aa+ABZ*I_ESP_K3xy3z_Py_aa;
    Double I_ESP_K3x4z_Dyz_aa = I_ESP_L3x5z_Py_aa+ABZ*I_ESP_K3x4z_Py_aa;
    Double I_ESP_K2x5y_Dyz_aa = I_ESP_L2x5yz_Py_aa+ABZ*I_ESP_K2x5y_Py_aa;
    Double I_ESP_K2x4yz_Dyz_aa = I_ESP_L2x4y2z_Py_aa+ABZ*I_ESP_K2x4yz_Py_aa;
    Double I_ESP_K2x3y2z_Dyz_aa = I_ESP_L2x3y3z_Py_aa+ABZ*I_ESP_K2x3y2z_Py_aa;
    Double I_ESP_K2x2y3z_Dyz_aa = I_ESP_L2x2y4z_Py_aa+ABZ*I_ESP_K2x2y3z_Py_aa;
    Double I_ESP_K2xy4z_Dyz_aa = I_ESP_L2xy5z_Py_aa+ABZ*I_ESP_K2xy4z_Py_aa;
    Double I_ESP_K2x5z_Dyz_aa = I_ESP_L2x6z_Py_aa+ABZ*I_ESP_K2x5z_Py_aa;
    Double I_ESP_Kx6y_Dyz_aa = I_ESP_Lx6yz_Py_aa+ABZ*I_ESP_Kx6y_Py_aa;
    Double I_ESP_Kx5yz_Dyz_aa = I_ESP_Lx5y2z_Py_aa+ABZ*I_ESP_Kx5yz_Py_aa;
    Double I_ESP_Kx4y2z_Dyz_aa = I_ESP_Lx4y3z_Py_aa+ABZ*I_ESP_Kx4y2z_Py_aa;
    Double I_ESP_Kx3y3z_Dyz_aa = I_ESP_Lx3y4z_Py_aa+ABZ*I_ESP_Kx3y3z_Py_aa;
    Double I_ESP_Kx2y4z_Dyz_aa = I_ESP_Lx2y5z_Py_aa+ABZ*I_ESP_Kx2y4z_Py_aa;
    Double I_ESP_Kxy5z_Dyz_aa = I_ESP_Lxy6z_Py_aa+ABZ*I_ESP_Kxy5z_Py_aa;
    Double I_ESP_Kx6z_Dyz_aa = I_ESP_Lx7z_Py_aa+ABZ*I_ESP_Kx6z_Py_aa;
    Double I_ESP_K7y_Dyz_aa = I_ESP_L7yz_Py_aa+ABZ*I_ESP_K7y_Py_aa;
    Double I_ESP_K6yz_Dyz_aa = I_ESP_L6y2z_Py_aa+ABZ*I_ESP_K6yz_Py_aa;
    Double I_ESP_K5y2z_Dyz_aa = I_ESP_L5y3z_Py_aa+ABZ*I_ESP_K5y2z_Py_aa;
    Double I_ESP_K4y3z_Dyz_aa = I_ESP_L4y4z_Py_aa+ABZ*I_ESP_K4y3z_Py_aa;
    Double I_ESP_K3y4z_Dyz_aa = I_ESP_L3y5z_Py_aa+ABZ*I_ESP_K3y4z_Py_aa;
    Double I_ESP_K2y5z_Dyz_aa = I_ESP_L2y6z_Py_aa+ABZ*I_ESP_K2y5z_Py_aa;
    Double I_ESP_Ky6z_Dyz_aa = I_ESP_Ly7z_Py_aa+ABZ*I_ESP_Ky6z_Py_aa;
    Double I_ESP_K7z_Dyz_aa = I_ESP_L8z_Py_aa+ABZ*I_ESP_K7z_Py_aa;
    Double I_ESP_K7x_D2z_aa = I_ESP_L7xz_Pz_aa+ABZ*I_ESP_K7x_Pz_aa;
    Double I_ESP_K6xy_D2z_aa = I_ESP_L6xyz_Pz_aa+ABZ*I_ESP_K6xy_Pz_aa;
    Double I_ESP_K6xz_D2z_aa = I_ESP_L6x2z_Pz_aa+ABZ*I_ESP_K6xz_Pz_aa;
    Double I_ESP_K5x2y_D2z_aa = I_ESP_L5x2yz_Pz_aa+ABZ*I_ESP_K5x2y_Pz_aa;
    Double I_ESP_K5xyz_D2z_aa = I_ESP_L5xy2z_Pz_aa+ABZ*I_ESP_K5xyz_Pz_aa;
    Double I_ESP_K5x2z_D2z_aa = I_ESP_L5x3z_Pz_aa+ABZ*I_ESP_K5x2z_Pz_aa;
    Double I_ESP_K4x3y_D2z_aa = I_ESP_L4x3yz_Pz_aa+ABZ*I_ESP_K4x3y_Pz_aa;
    Double I_ESP_K4x2yz_D2z_aa = I_ESP_L4x2y2z_Pz_aa+ABZ*I_ESP_K4x2yz_Pz_aa;
    Double I_ESP_K4xy2z_D2z_aa = I_ESP_L4xy3z_Pz_aa+ABZ*I_ESP_K4xy2z_Pz_aa;
    Double I_ESP_K4x3z_D2z_aa = I_ESP_L4x4z_Pz_aa+ABZ*I_ESP_K4x3z_Pz_aa;
    Double I_ESP_K3x4y_D2z_aa = I_ESP_L3x4yz_Pz_aa+ABZ*I_ESP_K3x4y_Pz_aa;
    Double I_ESP_K3x3yz_D2z_aa = I_ESP_L3x3y2z_Pz_aa+ABZ*I_ESP_K3x3yz_Pz_aa;
    Double I_ESP_K3x2y2z_D2z_aa = I_ESP_L3x2y3z_Pz_aa+ABZ*I_ESP_K3x2y2z_Pz_aa;
    Double I_ESP_K3xy3z_D2z_aa = I_ESP_L3xy4z_Pz_aa+ABZ*I_ESP_K3xy3z_Pz_aa;
    Double I_ESP_K3x4z_D2z_aa = I_ESP_L3x5z_Pz_aa+ABZ*I_ESP_K3x4z_Pz_aa;
    Double I_ESP_K2x5y_D2z_aa = I_ESP_L2x5yz_Pz_aa+ABZ*I_ESP_K2x5y_Pz_aa;
    Double I_ESP_K2x4yz_D2z_aa = I_ESP_L2x4y2z_Pz_aa+ABZ*I_ESP_K2x4yz_Pz_aa;
    Double I_ESP_K2x3y2z_D2z_aa = I_ESP_L2x3y3z_Pz_aa+ABZ*I_ESP_K2x3y2z_Pz_aa;
    Double I_ESP_K2x2y3z_D2z_aa = I_ESP_L2x2y4z_Pz_aa+ABZ*I_ESP_K2x2y3z_Pz_aa;
    Double I_ESP_K2xy4z_D2z_aa = I_ESP_L2xy5z_Pz_aa+ABZ*I_ESP_K2xy4z_Pz_aa;
    Double I_ESP_K2x5z_D2z_aa = I_ESP_L2x6z_Pz_aa+ABZ*I_ESP_K2x5z_Pz_aa;
    Double I_ESP_Kx6y_D2z_aa = I_ESP_Lx6yz_Pz_aa+ABZ*I_ESP_Kx6y_Pz_aa;
    Double I_ESP_Kx5yz_D2z_aa = I_ESP_Lx5y2z_Pz_aa+ABZ*I_ESP_Kx5yz_Pz_aa;
    Double I_ESP_Kx4y2z_D2z_aa = I_ESP_Lx4y3z_Pz_aa+ABZ*I_ESP_Kx4y2z_Pz_aa;
    Double I_ESP_Kx3y3z_D2z_aa = I_ESP_Lx3y4z_Pz_aa+ABZ*I_ESP_Kx3y3z_Pz_aa;
    Double I_ESP_Kx2y4z_D2z_aa = I_ESP_Lx2y5z_Pz_aa+ABZ*I_ESP_Kx2y4z_Pz_aa;
    Double I_ESP_Kxy5z_D2z_aa = I_ESP_Lxy6z_Pz_aa+ABZ*I_ESP_Kxy5z_Pz_aa;
    Double I_ESP_Kx6z_D2z_aa = I_ESP_Lx7z_Pz_aa+ABZ*I_ESP_Kx6z_Pz_aa;
    Double I_ESP_K7y_D2z_aa = I_ESP_L7yz_Pz_aa+ABZ*I_ESP_K7y_Pz_aa;
    Double I_ESP_K6yz_D2z_aa = I_ESP_L6y2z_Pz_aa+ABZ*I_ESP_K6yz_Pz_aa;
    Double I_ESP_K5y2z_D2z_aa = I_ESP_L5y3z_Pz_aa+ABZ*I_ESP_K5y2z_Pz_aa;
    Double I_ESP_K4y3z_D2z_aa = I_ESP_L4y4z_Pz_aa+ABZ*I_ESP_K4y3z_Pz_aa;
    Double I_ESP_K3y4z_D2z_aa = I_ESP_L3y5z_Pz_aa+ABZ*I_ESP_K3y4z_Pz_aa;
    Double I_ESP_K2y5z_D2z_aa = I_ESP_L2y6z_Pz_aa+ABZ*I_ESP_K2y5z_Pz_aa;
    Double I_ESP_Ky6z_D2z_aa = I_ESP_Ly7z_Pz_aa+ABZ*I_ESP_Ky6z_Pz_aa;
    Double I_ESP_K7z_D2z_aa = I_ESP_L8z_Pz_aa+ABZ*I_ESP_K7z_Pz_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D_aa
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_F_D
     ************************************************************/
    abcd[iGrid*756+0] = 4.0E0*I_ESP_K7x_D2x_aa-2.0E0*5*I_ESP_H5x_D2x_a-2.0E0*6*I_ESP_H5x_D2x_a+5*4*I_ESP_F3x_D2x;
    abcd[iGrid*756+1] = 4.0E0*I_ESP_K6xy_D2x_aa-2.0E0*4*I_ESP_H4xy_D2x_a-2.0E0*5*I_ESP_H4xy_D2x_a+4*3*I_ESP_F2xy_D2x;
    abcd[iGrid*756+2] = 4.0E0*I_ESP_K6xz_D2x_aa-2.0E0*4*I_ESP_H4xz_D2x_a-2.0E0*5*I_ESP_H4xz_D2x_a+4*3*I_ESP_F2xz_D2x;
    abcd[iGrid*756+3] = 4.0E0*I_ESP_K5x2y_D2x_aa-2.0E0*3*I_ESP_H3x2y_D2x_a-2.0E0*4*I_ESP_H3x2y_D2x_a+3*2*I_ESP_Fx2y_D2x;
    abcd[iGrid*756+4] = 4.0E0*I_ESP_K5xyz_D2x_aa-2.0E0*3*I_ESP_H3xyz_D2x_a-2.0E0*4*I_ESP_H3xyz_D2x_a+3*2*I_ESP_Fxyz_D2x;
    abcd[iGrid*756+5] = 4.0E0*I_ESP_K5x2z_D2x_aa-2.0E0*3*I_ESP_H3x2z_D2x_a-2.0E0*4*I_ESP_H3x2z_D2x_a+3*2*I_ESP_Fx2z_D2x;
    abcd[iGrid*756+6] = 4.0E0*I_ESP_K4x3y_D2x_aa-2.0E0*2*I_ESP_H2x3y_D2x_a-2.0E0*3*I_ESP_H2x3y_D2x_a+2*1*I_ESP_F3y_D2x;
    abcd[iGrid*756+7] = 4.0E0*I_ESP_K4x2yz_D2x_aa-2.0E0*2*I_ESP_H2x2yz_D2x_a-2.0E0*3*I_ESP_H2x2yz_D2x_a+2*1*I_ESP_F2yz_D2x;
    abcd[iGrid*756+8] = 4.0E0*I_ESP_K4xy2z_D2x_aa-2.0E0*2*I_ESP_H2xy2z_D2x_a-2.0E0*3*I_ESP_H2xy2z_D2x_a+2*1*I_ESP_Fy2z_D2x;
    abcd[iGrid*756+9] = 4.0E0*I_ESP_K4x3z_D2x_aa-2.0E0*2*I_ESP_H2x3z_D2x_a-2.0E0*3*I_ESP_H2x3z_D2x_a+2*1*I_ESP_F3z_D2x;
    abcd[iGrid*756+10] = 4.0E0*I_ESP_K3x4y_D2x_aa-2.0E0*1*I_ESP_Hx4y_D2x_a-2.0E0*2*I_ESP_Hx4y_D2x_a;
    abcd[iGrid*756+11] = 4.0E0*I_ESP_K3x3yz_D2x_aa-2.0E0*1*I_ESP_Hx3yz_D2x_a-2.0E0*2*I_ESP_Hx3yz_D2x_a;
    abcd[iGrid*756+12] = 4.0E0*I_ESP_K3x2y2z_D2x_aa-2.0E0*1*I_ESP_Hx2y2z_D2x_a-2.0E0*2*I_ESP_Hx2y2z_D2x_a;
    abcd[iGrid*756+13] = 4.0E0*I_ESP_K3xy3z_D2x_aa-2.0E0*1*I_ESP_Hxy3z_D2x_a-2.0E0*2*I_ESP_Hxy3z_D2x_a;
    abcd[iGrid*756+14] = 4.0E0*I_ESP_K3x4z_D2x_aa-2.0E0*1*I_ESP_Hx4z_D2x_a-2.0E0*2*I_ESP_Hx4z_D2x_a;
    abcd[iGrid*756+15] = 4.0E0*I_ESP_K2x5y_D2x_aa-2.0E0*1*I_ESP_H5y_D2x_a;
    abcd[iGrid*756+16] = 4.0E0*I_ESP_K2x4yz_D2x_aa-2.0E0*1*I_ESP_H4yz_D2x_a;
    abcd[iGrid*756+17] = 4.0E0*I_ESP_K2x3y2z_D2x_aa-2.0E0*1*I_ESP_H3y2z_D2x_a;
    abcd[iGrid*756+18] = 4.0E0*I_ESP_K2x2y3z_D2x_aa-2.0E0*1*I_ESP_H2y3z_D2x_a;
    abcd[iGrid*756+19] = 4.0E0*I_ESP_K2xy4z_D2x_aa-2.0E0*1*I_ESP_Hy4z_D2x_a;
    abcd[iGrid*756+20] = 4.0E0*I_ESP_K2x5z_D2x_aa-2.0E0*1*I_ESP_H5z_D2x_a;
    abcd[iGrid*756+21] = 4.0E0*I_ESP_K7x_Dxy_aa-2.0E0*5*I_ESP_H5x_Dxy_a-2.0E0*6*I_ESP_H5x_Dxy_a+5*4*I_ESP_F3x_Dxy;
    abcd[iGrid*756+22] = 4.0E0*I_ESP_K6xy_Dxy_aa-2.0E0*4*I_ESP_H4xy_Dxy_a-2.0E0*5*I_ESP_H4xy_Dxy_a+4*3*I_ESP_F2xy_Dxy;
    abcd[iGrid*756+23] = 4.0E0*I_ESP_K6xz_Dxy_aa-2.0E0*4*I_ESP_H4xz_Dxy_a-2.0E0*5*I_ESP_H4xz_Dxy_a+4*3*I_ESP_F2xz_Dxy;
    abcd[iGrid*756+24] = 4.0E0*I_ESP_K5x2y_Dxy_aa-2.0E0*3*I_ESP_H3x2y_Dxy_a-2.0E0*4*I_ESP_H3x2y_Dxy_a+3*2*I_ESP_Fx2y_Dxy;
    abcd[iGrid*756+25] = 4.0E0*I_ESP_K5xyz_Dxy_aa-2.0E0*3*I_ESP_H3xyz_Dxy_a-2.0E0*4*I_ESP_H3xyz_Dxy_a+3*2*I_ESP_Fxyz_Dxy;
    abcd[iGrid*756+26] = 4.0E0*I_ESP_K5x2z_Dxy_aa-2.0E0*3*I_ESP_H3x2z_Dxy_a-2.0E0*4*I_ESP_H3x2z_Dxy_a+3*2*I_ESP_Fx2z_Dxy;
    abcd[iGrid*756+27] = 4.0E0*I_ESP_K4x3y_Dxy_aa-2.0E0*2*I_ESP_H2x3y_Dxy_a-2.0E0*3*I_ESP_H2x3y_Dxy_a+2*1*I_ESP_F3y_Dxy;
    abcd[iGrid*756+28] = 4.0E0*I_ESP_K4x2yz_Dxy_aa-2.0E0*2*I_ESP_H2x2yz_Dxy_a-2.0E0*3*I_ESP_H2x2yz_Dxy_a+2*1*I_ESP_F2yz_Dxy;
    abcd[iGrid*756+29] = 4.0E0*I_ESP_K4xy2z_Dxy_aa-2.0E0*2*I_ESP_H2xy2z_Dxy_a-2.0E0*3*I_ESP_H2xy2z_Dxy_a+2*1*I_ESP_Fy2z_Dxy;
    abcd[iGrid*756+30] = 4.0E0*I_ESP_K4x3z_Dxy_aa-2.0E0*2*I_ESP_H2x3z_Dxy_a-2.0E0*3*I_ESP_H2x3z_Dxy_a+2*1*I_ESP_F3z_Dxy;
    abcd[iGrid*756+31] = 4.0E0*I_ESP_K3x4y_Dxy_aa-2.0E0*1*I_ESP_Hx4y_Dxy_a-2.0E0*2*I_ESP_Hx4y_Dxy_a;
    abcd[iGrid*756+32] = 4.0E0*I_ESP_K3x3yz_Dxy_aa-2.0E0*1*I_ESP_Hx3yz_Dxy_a-2.0E0*2*I_ESP_Hx3yz_Dxy_a;
    abcd[iGrid*756+33] = 4.0E0*I_ESP_K3x2y2z_Dxy_aa-2.0E0*1*I_ESP_Hx2y2z_Dxy_a-2.0E0*2*I_ESP_Hx2y2z_Dxy_a;
    abcd[iGrid*756+34] = 4.0E0*I_ESP_K3xy3z_Dxy_aa-2.0E0*1*I_ESP_Hxy3z_Dxy_a-2.0E0*2*I_ESP_Hxy3z_Dxy_a;
    abcd[iGrid*756+35] = 4.0E0*I_ESP_K3x4z_Dxy_aa-2.0E0*1*I_ESP_Hx4z_Dxy_a-2.0E0*2*I_ESP_Hx4z_Dxy_a;
    abcd[iGrid*756+36] = 4.0E0*I_ESP_K2x5y_Dxy_aa-2.0E0*1*I_ESP_H5y_Dxy_a;
    abcd[iGrid*756+37] = 4.0E0*I_ESP_K2x4yz_Dxy_aa-2.0E0*1*I_ESP_H4yz_Dxy_a;
    abcd[iGrid*756+38] = 4.0E0*I_ESP_K2x3y2z_Dxy_aa-2.0E0*1*I_ESP_H3y2z_Dxy_a;
    abcd[iGrid*756+39] = 4.0E0*I_ESP_K2x2y3z_Dxy_aa-2.0E0*1*I_ESP_H2y3z_Dxy_a;
    abcd[iGrid*756+40] = 4.0E0*I_ESP_K2xy4z_Dxy_aa-2.0E0*1*I_ESP_Hy4z_Dxy_a;
    abcd[iGrid*756+41] = 4.0E0*I_ESP_K2x5z_Dxy_aa-2.0E0*1*I_ESP_H5z_Dxy_a;
    abcd[iGrid*756+42] = 4.0E0*I_ESP_K7x_Dxz_aa-2.0E0*5*I_ESP_H5x_Dxz_a-2.0E0*6*I_ESP_H5x_Dxz_a+5*4*I_ESP_F3x_Dxz;
    abcd[iGrid*756+43] = 4.0E0*I_ESP_K6xy_Dxz_aa-2.0E0*4*I_ESP_H4xy_Dxz_a-2.0E0*5*I_ESP_H4xy_Dxz_a+4*3*I_ESP_F2xy_Dxz;
    abcd[iGrid*756+44] = 4.0E0*I_ESP_K6xz_Dxz_aa-2.0E0*4*I_ESP_H4xz_Dxz_a-2.0E0*5*I_ESP_H4xz_Dxz_a+4*3*I_ESP_F2xz_Dxz;
    abcd[iGrid*756+45] = 4.0E0*I_ESP_K5x2y_Dxz_aa-2.0E0*3*I_ESP_H3x2y_Dxz_a-2.0E0*4*I_ESP_H3x2y_Dxz_a+3*2*I_ESP_Fx2y_Dxz;
    abcd[iGrid*756+46] = 4.0E0*I_ESP_K5xyz_Dxz_aa-2.0E0*3*I_ESP_H3xyz_Dxz_a-2.0E0*4*I_ESP_H3xyz_Dxz_a+3*2*I_ESP_Fxyz_Dxz;
    abcd[iGrid*756+47] = 4.0E0*I_ESP_K5x2z_Dxz_aa-2.0E0*3*I_ESP_H3x2z_Dxz_a-2.0E0*4*I_ESP_H3x2z_Dxz_a+3*2*I_ESP_Fx2z_Dxz;
    abcd[iGrid*756+48] = 4.0E0*I_ESP_K4x3y_Dxz_aa-2.0E0*2*I_ESP_H2x3y_Dxz_a-2.0E0*3*I_ESP_H2x3y_Dxz_a+2*1*I_ESP_F3y_Dxz;
    abcd[iGrid*756+49] = 4.0E0*I_ESP_K4x2yz_Dxz_aa-2.0E0*2*I_ESP_H2x2yz_Dxz_a-2.0E0*3*I_ESP_H2x2yz_Dxz_a+2*1*I_ESP_F2yz_Dxz;
    abcd[iGrid*756+50] = 4.0E0*I_ESP_K4xy2z_Dxz_aa-2.0E0*2*I_ESP_H2xy2z_Dxz_a-2.0E0*3*I_ESP_H2xy2z_Dxz_a+2*1*I_ESP_Fy2z_Dxz;
    abcd[iGrid*756+51] = 4.0E0*I_ESP_K4x3z_Dxz_aa-2.0E0*2*I_ESP_H2x3z_Dxz_a-2.0E0*3*I_ESP_H2x3z_Dxz_a+2*1*I_ESP_F3z_Dxz;
    abcd[iGrid*756+52] = 4.0E0*I_ESP_K3x4y_Dxz_aa-2.0E0*1*I_ESP_Hx4y_Dxz_a-2.0E0*2*I_ESP_Hx4y_Dxz_a;
    abcd[iGrid*756+53] = 4.0E0*I_ESP_K3x3yz_Dxz_aa-2.0E0*1*I_ESP_Hx3yz_Dxz_a-2.0E0*2*I_ESP_Hx3yz_Dxz_a;
    abcd[iGrid*756+54] = 4.0E0*I_ESP_K3x2y2z_Dxz_aa-2.0E0*1*I_ESP_Hx2y2z_Dxz_a-2.0E0*2*I_ESP_Hx2y2z_Dxz_a;
    abcd[iGrid*756+55] = 4.0E0*I_ESP_K3xy3z_Dxz_aa-2.0E0*1*I_ESP_Hxy3z_Dxz_a-2.0E0*2*I_ESP_Hxy3z_Dxz_a;
    abcd[iGrid*756+56] = 4.0E0*I_ESP_K3x4z_Dxz_aa-2.0E0*1*I_ESP_Hx4z_Dxz_a-2.0E0*2*I_ESP_Hx4z_Dxz_a;
    abcd[iGrid*756+57] = 4.0E0*I_ESP_K2x5y_Dxz_aa-2.0E0*1*I_ESP_H5y_Dxz_a;
    abcd[iGrid*756+58] = 4.0E0*I_ESP_K2x4yz_Dxz_aa-2.0E0*1*I_ESP_H4yz_Dxz_a;
    abcd[iGrid*756+59] = 4.0E0*I_ESP_K2x3y2z_Dxz_aa-2.0E0*1*I_ESP_H3y2z_Dxz_a;
    abcd[iGrid*756+60] = 4.0E0*I_ESP_K2x2y3z_Dxz_aa-2.0E0*1*I_ESP_H2y3z_Dxz_a;
    abcd[iGrid*756+61] = 4.0E0*I_ESP_K2xy4z_Dxz_aa-2.0E0*1*I_ESP_Hy4z_Dxz_a;
    abcd[iGrid*756+62] = 4.0E0*I_ESP_K2x5z_Dxz_aa-2.0E0*1*I_ESP_H5z_Dxz_a;
    abcd[iGrid*756+63] = 4.0E0*I_ESP_K7x_D2y_aa-2.0E0*5*I_ESP_H5x_D2y_a-2.0E0*6*I_ESP_H5x_D2y_a+5*4*I_ESP_F3x_D2y;
    abcd[iGrid*756+64] = 4.0E0*I_ESP_K6xy_D2y_aa-2.0E0*4*I_ESP_H4xy_D2y_a-2.0E0*5*I_ESP_H4xy_D2y_a+4*3*I_ESP_F2xy_D2y;
    abcd[iGrid*756+65] = 4.0E0*I_ESP_K6xz_D2y_aa-2.0E0*4*I_ESP_H4xz_D2y_a-2.0E0*5*I_ESP_H4xz_D2y_a+4*3*I_ESP_F2xz_D2y;
    abcd[iGrid*756+66] = 4.0E0*I_ESP_K5x2y_D2y_aa-2.0E0*3*I_ESP_H3x2y_D2y_a-2.0E0*4*I_ESP_H3x2y_D2y_a+3*2*I_ESP_Fx2y_D2y;
    abcd[iGrid*756+67] = 4.0E0*I_ESP_K5xyz_D2y_aa-2.0E0*3*I_ESP_H3xyz_D2y_a-2.0E0*4*I_ESP_H3xyz_D2y_a+3*2*I_ESP_Fxyz_D2y;
    abcd[iGrid*756+68] = 4.0E0*I_ESP_K5x2z_D2y_aa-2.0E0*3*I_ESP_H3x2z_D2y_a-2.0E0*4*I_ESP_H3x2z_D2y_a+3*2*I_ESP_Fx2z_D2y;
    abcd[iGrid*756+69] = 4.0E0*I_ESP_K4x3y_D2y_aa-2.0E0*2*I_ESP_H2x3y_D2y_a-2.0E0*3*I_ESP_H2x3y_D2y_a+2*1*I_ESP_F3y_D2y;
    abcd[iGrid*756+70] = 4.0E0*I_ESP_K4x2yz_D2y_aa-2.0E0*2*I_ESP_H2x2yz_D2y_a-2.0E0*3*I_ESP_H2x2yz_D2y_a+2*1*I_ESP_F2yz_D2y;
    abcd[iGrid*756+71] = 4.0E0*I_ESP_K4xy2z_D2y_aa-2.0E0*2*I_ESP_H2xy2z_D2y_a-2.0E0*3*I_ESP_H2xy2z_D2y_a+2*1*I_ESP_Fy2z_D2y;
    abcd[iGrid*756+72] = 4.0E0*I_ESP_K4x3z_D2y_aa-2.0E0*2*I_ESP_H2x3z_D2y_a-2.0E0*3*I_ESP_H2x3z_D2y_a+2*1*I_ESP_F3z_D2y;
    abcd[iGrid*756+73] = 4.0E0*I_ESP_K3x4y_D2y_aa-2.0E0*1*I_ESP_Hx4y_D2y_a-2.0E0*2*I_ESP_Hx4y_D2y_a;
    abcd[iGrid*756+74] = 4.0E0*I_ESP_K3x3yz_D2y_aa-2.0E0*1*I_ESP_Hx3yz_D2y_a-2.0E0*2*I_ESP_Hx3yz_D2y_a;
    abcd[iGrid*756+75] = 4.0E0*I_ESP_K3x2y2z_D2y_aa-2.0E0*1*I_ESP_Hx2y2z_D2y_a-2.0E0*2*I_ESP_Hx2y2z_D2y_a;
    abcd[iGrid*756+76] = 4.0E0*I_ESP_K3xy3z_D2y_aa-2.0E0*1*I_ESP_Hxy3z_D2y_a-2.0E0*2*I_ESP_Hxy3z_D2y_a;
    abcd[iGrid*756+77] = 4.0E0*I_ESP_K3x4z_D2y_aa-2.0E0*1*I_ESP_Hx4z_D2y_a-2.0E0*2*I_ESP_Hx4z_D2y_a;
    abcd[iGrid*756+78] = 4.0E0*I_ESP_K2x5y_D2y_aa-2.0E0*1*I_ESP_H5y_D2y_a;
    abcd[iGrid*756+79] = 4.0E0*I_ESP_K2x4yz_D2y_aa-2.0E0*1*I_ESP_H4yz_D2y_a;
    abcd[iGrid*756+80] = 4.0E0*I_ESP_K2x3y2z_D2y_aa-2.0E0*1*I_ESP_H3y2z_D2y_a;
    abcd[iGrid*756+81] = 4.0E0*I_ESP_K2x2y3z_D2y_aa-2.0E0*1*I_ESP_H2y3z_D2y_a;
    abcd[iGrid*756+82] = 4.0E0*I_ESP_K2xy4z_D2y_aa-2.0E0*1*I_ESP_Hy4z_D2y_a;
    abcd[iGrid*756+83] = 4.0E0*I_ESP_K2x5z_D2y_aa-2.0E0*1*I_ESP_H5z_D2y_a;
    abcd[iGrid*756+84] = 4.0E0*I_ESP_K7x_Dyz_aa-2.0E0*5*I_ESP_H5x_Dyz_a-2.0E0*6*I_ESP_H5x_Dyz_a+5*4*I_ESP_F3x_Dyz;
    abcd[iGrid*756+85] = 4.0E0*I_ESP_K6xy_Dyz_aa-2.0E0*4*I_ESP_H4xy_Dyz_a-2.0E0*5*I_ESP_H4xy_Dyz_a+4*3*I_ESP_F2xy_Dyz;
    abcd[iGrid*756+86] = 4.0E0*I_ESP_K6xz_Dyz_aa-2.0E0*4*I_ESP_H4xz_Dyz_a-2.0E0*5*I_ESP_H4xz_Dyz_a+4*3*I_ESP_F2xz_Dyz;
    abcd[iGrid*756+87] = 4.0E0*I_ESP_K5x2y_Dyz_aa-2.0E0*3*I_ESP_H3x2y_Dyz_a-2.0E0*4*I_ESP_H3x2y_Dyz_a+3*2*I_ESP_Fx2y_Dyz;
    abcd[iGrid*756+88] = 4.0E0*I_ESP_K5xyz_Dyz_aa-2.0E0*3*I_ESP_H3xyz_Dyz_a-2.0E0*4*I_ESP_H3xyz_Dyz_a+3*2*I_ESP_Fxyz_Dyz;
    abcd[iGrid*756+89] = 4.0E0*I_ESP_K5x2z_Dyz_aa-2.0E0*3*I_ESP_H3x2z_Dyz_a-2.0E0*4*I_ESP_H3x2z_Dyz_a+3*2*I_ESP_Fx2z_Dyz;
    abcd[iGrid*756+90] = 4.0E0*I_ESP_K4x3y_Dyz_aa-2.0E0*2*I_ESP_H2x3y_Dyz_a-2.0E0*3*I_ESP_H2x3y_Dyz_a+2*1*I_ESP_F3y_Dyz;
    abcd[iGrid*756+91] = 4.0E0*I_ESP_K4x2yz_Dyz_aa-2.0E0*2*I_ESP_H2x2yz_Dyz_a-2.0E0*3*I_ESP_H2x2yz_Dyz_a+2*1*I_ESP_F2yz_Dyz;
    abcd[iGrid*756+92] = 4.0E0*I_ESP_K4xy2z_Dyz_aa-2.0E0*2*I_ESP_H2xy2z_Dyz_a-2.0E0*3*I_ESP_H2xy2z_Dyz_a+2*1*I_ESP_Fy2z_Dyz;
    abcd[iGrid*756+93] = 4.0E0*I_ESP_K4x3z_Dyz_aa-2.0E0*2*I_ESP_H2x3z_Dyz_a-2.0E0*3*I_ESP_H2x3z_Dyz_a+2*1*I_ESP_F3z_Dyz;
    abcd[iGrid*756+94] = 4.0E0*I_ESP_K3x4y_Dyz_aa-2.0E0*1*I_ESP_Hx4y_Dyz_a-2.0E0*2*I_ESP_Hx4y_Dyz_a;
    abcd[iGrid*756+95] = 4.0E0*I_ESP_K3x3yz_Dyz_aa-2.0E0*1*I_ESP_Hx3yz_Dyz_a-2.0E0*2*I_ESP_Hx3yz_Dyz_a;
    abcd[iGrid*756+96] = 4.0E0*I_ESP_K3x2y2z_Dyz_aa-2.0E0*1*I_ESP_Hx2y2z_Dyz_a-2.0E0*2*I_ESP_Hx2y2z_Dyz_a;
    abcd[iGrid*756+97] = 4.0E0*I_ESP_K3xy3z_Dyz_aa-2.0E0*1*I_ESP_Hxy3z_Dyz_a-2.0E0*2*I_ESP_Hxy3z_Dyz_a;
    abcd[iGrid*756+98] = 4.0E0*I_ESP_K3x4z_Dyz_aa-2.0E0*1*I_ESP_Hx4z_Dyz_a-2.0E0*2*I_ESP_Hx4z_Dyz_a;
    abcd[iGrid*756+99] = 4.0E0*I_ESP_K2x5y_Dyz_aa-2.0E0*1*I_ESP_H5y_Dyz_a;
    abcd[iGrid*756+100] = 4.0E0*I_ESP_K2x4yz_Dyz_aa-2.0E0*1*I_ESP_H4yz_Dyz_a;
    abcd[iGrid*756+101] = 4.0E0*I_ESP_K2x3y2z_Dyz_aa-2.0E0*1*I_ESP_H3y2z_Dyz_a;
    abcd[iGrid*756+102] = 4.0E0*I_ESP_K2x2y3z_Dyz_aa-2.0E0*1*I_ESP_H2y3z_Dyz_a;
    abcd[iGrid*756+103] = 4.0E0*I_ESP_K2xy4z_Dyz_aa-2.0E0*1*I_ESP_Hy4z_Dyz_a;
    abcd[iGrid*756+104] = 4.0E0*I_ESP_K2x5z_Dyz_aa-2.0E0*1*I_ESP_H5z_Dyz_a;
    abcd[iGrid*756+105] = 4.0E0*I_ESP_K7x_D2z_aa-2.0E0*5*I_ESP_H5x_D2z_a-2.0E0*6*I_ESP_H5x_D2z_a+5*4*I_ESP_F3x_D2z;
    abcd[iGrid*756+106] = 4.0E0*I_ESP_K6xy_D2z_aa-2.0E0*4*I_ESP_H4xy_D2z_a-2.0E0*5*I_ESP_H4xy_D2z_a+4*3*I_ESP_F2xy_D2z;
    abcd[iGrid*756+107] = 4.0E0*I_ESP_K6xz_D2z_aa-2.0E0*4*I_ESP_H4xz_D2z_a-2.0E0*5*I_ESP_H4xz_D2z_a+4*3*I_ESP_F2xz_D2z;
    abcd[iGrid*756+108] = 4.0E0*I_ESP_K5x2y_D2z_aa-2.0E0*3*I_ESP_H3x2y_D2z_a-2.0E0*4*I_ESP_H3x2y_D2z_a+3*2*I_ESP_Fx2y_D2z;
    abcd[iGrid*756+109] = 4.0E0*I_ESP_K5xyz_D2z_aa-2.0E0*3*I_ESP_H3xyz_D2z_a-2.0E0*4*I_ESP_H3xyz_D2z_a+3*2*I_ESP_Fxyz_D2z;
    abcd[iGrid*756+110] = 4.0E0*I_ESP_K5x2z_D2z_aa-2.0E0*3*I_ESP_H3x2z_D2z_a-2.0E0*4*I_ESP_H3x2z_D2z_a+3*2*I_ESP_Fx2z_D2z;
    abcd[iGrid*756+111] = 4.0E0*I_ESP_K4x3y_D2z_aa-2.0E0*2*I_ESP_H2x3y_D2z_a-2.0E0*3*I_ESP_H2x3y_D2z_a+2*1*I_ESP_F3y_D2z;
    abcd[iGrid*756+112] = 4.0E0*I_ESP_K4x2yz_D2z_aa-2.0E0*2*I_ESP_H2x2yz_D2z_a-2.0E0*3*I_ESP_H2x2yz_D2z_a+2*1*I_ESP_F2yz_D2z;
    abcd[iGrid*756+113] = 4.0E0*I_ESP_K4xy2z_D2z_aa-2.0E0*2*I_ESP_H2xy2z_D2z_a-2.0E0*3*I_ESP_H2xy2z_D2z_a+2*1*I_ESP_Fy2z_D2z;
    abcd[iGrid*756+114] = 4.0E0*I_ESP_K4x3z_D2z_aa-2.0E0*2*I_ESP_H2x3z_D2z_a-2.0E0*3*I_ESP_H2x3z_D2z_a+2*1*I_ESP_F3z_D2z;
    abcd[iGrid*756+115] = 4.0E0*I_ESP_K3x4y_D2z_aa-2.0E0*1*I_ESP_Hx4y_D2z_a-2.0E0*2*I_ESP_Hx4y_D2z_a;
    abcd[iGrid*756+116] = 4.0E0*I_ESP_K3x3yz_D2z_aa-2.0E0*1*I_ESP_Hx3yz_D2z_a-2.0E0*2*I_ESP_Hx3yz_D2z_a;
    abcd[iGrid*756+117] = 4.0E0*I_ESP_K3x2y2z_D2z_aa-2.0E0*1*I_ESP_Hx2y2z_D2z_a-2.0E0*2*I_ESP_Hx2y2z_D2z_a;
    abcd[iGrid*756+118] = 4.0E0*I_ESP_K3xy3z_D2z_aa-2.0E0*1*I_ESP_Hxy3z_D2z_a-2.0E0*2*I_ESP_Hxy3z_D2z_a;
    abcd[iGrid*756+119] = 4.0E0*I_ESP_K3x4z_D2z_aa-2.0E0*1*I_ESP_Hx4z_D2z_a-2.0E0*2*I_ESP_Hx4z_D2z_a;
    abcd[iGrid*756+120] = 4.0E0*I_ESP_K2x5y_D2z_aa-2.0E0*1*I_ESP_H5y_D2z_a;
    abcd[iGrid*756+121] = 4.0E0*I_ESP_K2x4yz_D2z_aa-2.0E0*1*I_ESP_H4yz_D2z_a;
    abcd[iGrid*756+122] = 4.0E0*I_ESP_K2x3y2z_D2z_aa-2.0E0*1*I_ESP_H3y2z_D2z_a;
    abcd[iGrid*756+123] = 4.0E0*I_ESP_K2x2y3z_D2z_aa-2.0E0*1*I_ESP_H2y3z_D2z_a;
    abcd[iGrid*756+124] = 4.0E0*I_ESP_K2xy4z_D2z_aa-2.0E0*1*I_ESP_Hy4z_D2z_a;
    abcd[iGrid*756+125] = 4.0E0*I_ESP_K2x5z_D2z_aa-2.0E0*1*I_ESP_H5z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D_aa
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_F_D
     ************************************************************/
    abcd[iGrid*756+126] = 4.0E0*I_ESP_K6xy_D2x_aa-2.0E0*5*I_ESP_H4xy_D2x_a;
    abcd[iGrid*756+127] = 4.0E0*I_ESP_K5x2y_D2x_aa-2.0E0*1*I_ESP_H5x_D2x_a-2.0E0*4*I_ESP_H3x2y_D2x_a+4*1*I_ESP_F3x_D2x;
    abcd[iGrid*756+128] = 4.0E0*I_ESP_K5xyz_D2x_aa-2.0E0*4*I_ESP_H3xyz_D2x_a;
    abcd[iGrid*756+129] = 4.0E0*I_ESP_K4x3y_D2x_aa-2.0E0*2*I_ESP_H4xy_D2x_a-2.0E0*3*I_ESP_H2x3y_D2x_a+3*2*I_ESP_F2xy_D2x;
    abcd[iGrid*756+130] = 4.0E0*I_ESP_K4x2yz_D2x_aa-2.0E0*1*I_ESP_H4xz_D2x_a-2.0E0*3*I_ESP_H2x2yz_D2x_a+3*1*I_ESP_F2xz_D2x;
    abcd[iGrid*756+131] = 4.0E0*I_ESP_K4xy2z_D2x_aa-2.0E0*3*I_ESP_H2xy2z_D2x_a;
    abcd[iGrid*756+132] = 4.0E0*I_ESP_K3x4y_D2x_aa-2.0E0*3*I_ESP_H3x2y_D2x_a-2.0E0*2*I_ESP_Hx4y_D2x_a+2*3*I_ESP_Fx2y_D2x;
    abcd[iGrid*756+133] = 4.0E0*I_ESP_K3x3yz_D2x_aa-2.0E0*2*I_ESP_H3xyz_D2x_a-2.0E0*2*I_ESP_Hx3yz_D2x_a+2*2*I_ESP_Fxyz_D2x;
    abcd[iGrid*756+134] = 4.0E0*I_ESP_K3x2y2z_D2x_aa-2.0E0*1*I_ESP_H3x2z_D2x_a-2.0E0*2*I_ESP_Hx2y2z_D2x_a+2*1*I_ESP_Fx2z_D2x;
    abcd[iGrid*756+135] = 4.0E0*I_ESP_K3xy3z_D2x_aa-2.0E0*2*I_ESP_Hxy3z_D2x_a;
    abcd[iGrid*756+136] = 4.0E0*I_ESP_K2x5y_D2x_aa-2.0E0*4*I_ESP_H2x3y_D2x_a-2.0E0*1*I_ESP_H5y_D2x_a+4*I_ESP_F3y_D2x;
    abcd[iGrid*756+137] = 4.0E0*I_ESP_K2x4yz_D2x_aa-2.0E0*3*I_ESP_H2x2yz_D2x_a-2.0E0*1*I_ESP_H4yz_D2x_a+3*I_ESP_F2yz_D2x;
    abcd[iGrid*756+138] = 4.0E0*I_ESP_K2x3y2z_D2x_aa-2.0E0*2*I_ESP_H2xy2z_D2x_a-2.0E0*1*I_ESP_H3y2z_D2x_a+2*I_ESP_Fy2z_D2x;
    abcd[iGrid*756+139] = 4.0E0*I_ESP_K2x2y3z_D2x_aa-2.0E0*1*I_ESP_H2x3z_D2x_a-2.0E0*1*I_ESP_H2y3z_D2x_a+1*I_ESP_F3z_D2x;
    abcd[iGrid*756+140] = 4.0E0*I_ESP_K2xy4z_D2x_aa-2.0E0*1*I_ESP_Hy4z_D2x_a;
    abcd[iGrid*756+141] = 4.0E0*I_ESP_Kx6y_D2x_aa-2.0E0*5*I_ESP_Hx4y_D2x_a;
    abcd[iGrid*756+142] = 4.0E0*I_ESP_Kx5yz_D2x_aa-2.0E0*4*I_ESP_Hx3yz_D2x_a;
    abcd[iGrid*756+143] = 4.0E0*I_ESP_Kx4y2z_D2x_aa-2.0E0*3*I_ESP_Hx2y2z_D2x_a;
    abcd[iGrid*756+144] = 4.0E0*I_ESP_Kx3y3z_D2x_aa-2.0E0*2*I_ESP_Hxy3z_D2x_a;
    abcd[iGrid*756+145] = 4.0E0*I_ESP_Kx2y4z_D2x_aa-2.0E0*1*I_ESP_Hx4z_D2x_a;
    abcd[iGrid*756+146] = 4.0E0*I_ESP_Kxy5z_D2x_aa;
    abcd[iGrid*756+147] = 4.0E0*I_ESP_K6xy_Dxy_aa-2.0E0*5*I_ESP_H4xy_Dxy_a;
    abcd[iGrid*756+148] = 4.0E0*I_ESP_K5x2y_Dxy_aa-2.0E0*1*I_ESP_H5x_Dxy_a-2.0E0*4*I_ESP_H3x2y_Dxy_a+4*1*I_ESP_F3x_Dxy;
    abcd[iGrid*756+149] = 4.0E0*I_ESP_K5xyz_Dxy_aa-2.0E0*4*I_ESP_H3xyz_Dxy_a;
    abcd[iGrid*756+150] = 4.0E0*I_ESP_K4x3y_Dxy_aa-2.0E0*2*I_ESP_H4xy_Dxy_a-2.0E0*3*I_ESP_H2x3y_Dxy_a+3*2*I_ESP_F2xy_Dxy;
    abcd[iGrid*756+151] = 4.0E0*I_ESP_K4x2yz_Dxy_aa-2.0E0*1*I_ESP_H4xz_Dxy_a-2.0E0*3*I_ESP_H2x2yz_Dxy_a+3*1*I_ESP_F2xz_Dxy;
    abcd[iGrid*756+152] = 4.0E0*I_ESP_K4xy2z_Dxy_aa-2.0E0*3*I_ESP_H2xy2z_Dxy_a;
    abcd[iGrid*756+153] = 4.0E0*I_ESP_K3x4y_Dxy_aa-2.0E0*3*I_ESP_H3x2y_Dxy_a-2.0E0*2*I_ESP_Hx4y_Dxy_a+2*3*I_ESP_Fx2y_Dxy;
    abcd[iGrid*756+154] = 4.0E0*I_ESP_K3x3yz_Dxy_aa-2.0E0*2*I_ESP_H3xyz_Dxy_a-2.0E0*2*I_ESP_Hx3yz_Dxy_a+2*2*I_ESP_Fxyz_Dxy;
    abcd[iGrid*756+155] = 4.0E0*I_ESP_K3x2y2z_Dxy_aa-2.0E0*1*I_ESP_H3x2z_Dxy_a-2.0E0*2*I_ESP_Hx2y2z_Dxy_a+2*1*I_ESP_Fx2z_Dxy;
    abcd[iGrid*756+156] = 4.0E0*I_ESP_K3xy3z_Dxy_aa-2.0E0*2*I_ESP_Hxy3z_Dxy_a;
    abcd[iGrid*756+157] = 4.0E0*I_ESP_K2x5y_Dxy_aa-2.0E0*4*I_ESP_H2x3y_Dxy_a-2.0E0*1*I_ESP_H5y_Dxy_a+4*I_ESP_F3y_Dxy;
    abcd[iGrid*756+158] = 4.0E0*I_ESP_K2x4yz_Dxy_aa-2.0E0*3*I_ESP_H2x2yz_Dxy_a-2.0E0*1*I_ESP_H4yz_Dxy_a+3*I_ESP_F2yz_Dxy;
    abcd[iGrid*756+159] = 4.0E0*I_ESP_K2x3y2z_Dxy_aa-2.0E0*2*I_ESP_H2xy2z_Dxy_a-2.0E0*1*I_ESP_H3y2z_Dxy_a+2*I_ESP_Fy2z_Dxy;
    abcd[iGrid*756+160] = 4.0E0*I_ESP_K2x2y3z_Dxy_aa-2.0E0*1*I_ESP_H2x3z_Dxy_a-2.0E0*1*I_ESP_H2y3z_Dxy_a+1*I_ESP_F3z_Dxy;
    abcd[iGrid*756+161] = 4.0E0*I_ESP_K2xy4z_Dxy_aa-2.0E0*1*I_ESP_Hy4z_Dxy_a;
    abcd[iGrid*756+162] = 4.0E0*I_ESP_Kx6y_Dxy_aa-2.0E0*5*I_ESP_Hx4y_Dxy_a;
    abcd[iGrid*756+163] = 4.0E0*I_ESP_Kx5yz_Dxy_aa-2.0E0*4*I_ESP_Hx3yz_Dxy_a;
    abcd[iGrid*756+164] = 4.0E0*I_ESP_Kx4y2z_Dxy_aa-2.0E0*3*I_ESP_Hx2y2z_Dxy_a;
    abcd[iGrid*756+165] = 4.0E0*I_ESP_Kx3y3z_Dxy_aa-2.0E0*2*I_ESP_Hxy3z_Dxy_a;
    abcd[iGrid*756+166] = 4.0E0*I_ESP_Kx2y4z_Dxy_aa-2.0E0*1*I_ESP_Hx4z_Dxy_a;
    abcd[iGrid*756+167] = 4.0E0*I_ESP_Kxy5z_Dxy_aa;
    abcd[iGrid*756+168] = 4.0E0*I_ESP_K6xy_Dxz_aa-2.0E0*5*I_ESP_H4xy_Dxz_a;
    abcd[iGrid*756+169] = 4.0E0*I_ESP_K5x2y_Dxz_aa-2.0E0*1*I_ESP_H5x_Dxz_a-2.0E0*4*I_ESP_H3x2y_Dxz_a+4*1*I_ESP_F3x_Dxz;
    abcd[iGrid*756+170] = 4.0E0*I_ESP_K5xyz_Dxz_aa-2.0E0*4*I_ESP_H3xyz_Dxz_a;
    abcd[iGrid*756+171] = 4.0E0*I_ESP_K4x3y_Dxz_aa-2.0E0*2*I_ESP_H4xy_Dxz_a-2.0E0*3*I_ESP_H2x3y_Dxz_a+3*2*I_ESP_F2xy_Dxz;
    abcd[iGrid*756+172] = 4.0E0*I_ESP_K4x2yz_Dxz_aa-2.0E0*1*I_ESP_H4xz_Dxz_a-2.0E0*3*I_ESP_H2x2yz_Dxz_a+3*1*I_ESP_F2xz_Dxz;
    abcd[iGrid*756+173] = 4.0E0*I_ESP_K4xy2z_Dxz_aa-2.0E0*3*I_ESP_H2xy2z_Dxz_a;
    abcd[iGrid*756+174] = 4.0E0*I_ESP_K3x4y_Dxz_aa-2.0E0*3*I_ESP_H3x2y_Dxz_a-2.0E0*2*I_ESP_Hx4y_Dxz_a+2*3*I_ESP_Fx2y_Dxz;
    abcd[iGrid*756+175] = 4.0E0*I_ESP_K3x3yz_Dxz_aa-2.0E0*2*I_ESP_H3xyz_Dxz_a-2.0E0*2*I_ESP_Hx3yz_Dxz_a+2*2*I_ESP_Fxyz_Dxz;
    abcd[iGrid*756+176] = 4.0E0*I_ESP_K3x2y2z_Dxz_aa-2.0E0*1*I_ESP_H3x2z_Dxz_a-2.0E0*2*I_ESP_Hx2y2z_Dxz_a+2*1*I_ESP_Fx2z_Dxz;
    abcd[iGrid*756+177] = 4.0E0*I_ESP_K3xy3z_Dxz_aa-2.0E0*2*I_ESP_Hxy3z_Dxz_a;
    abcd[iGrid*756+178] = 4.0E0*I_ESP_K2x5y_Dxz_aa-2.0E0*4*I_ESP_H2x3y_Dxz_a-2.0E0*1*I_ESP_H5y_Dxz_a+4*I_ESP_F3y_Dxz;
    abcd[iGrid*756+179] = 4.0E0*I_ESP_K2x4yz_Dxz_aa-2.0E0*3*I_ESP_H2x2yz_Dxz_a-2.0E0*1*I_ESP_H4yz_Dxz_a+3*I_ESP_F2yz_Dxz;
    abcd[iGrid*756+180] = 4.0E0*I_ESP_K2x3y2z_Dxz_aa-2.0E0*2*I_ESP_H2xy2z_Dxz_a-2.0E0*1*I_ESP_H3y2z_Dxz_a+2*I_ESP_Fy2z_Dxz;
    abcd[iGrid*756+181] = 4.0E0*I_ESP_K2x2y3z_Dxz_aa-2.0E0*1*I_ESP_H2x3z_Dxz_a-2.0E0*1*I_ESP_H2y3z_Dxz_a+1*I_ESP_F3z_Dxz;
    abcd[iGrid*756+182] = 4.0E0*I_ESP_K2xy4z_Dxz_aa-2.0E0*1*I_ESP_Hy4z_Dxz_a;
    abcd[iGrid*756+183] = 4.0E0*I_ESP_Kx6y_Dxz_aa-2.0E0*5*I_ESP_Hx4y_Dxz_a;
    abcd[iGrid*756+184] = 4.0E0*I_ESP_Kx5yz_Dxz_aa-2.0E0*4*I_ESP_Hx3yz_Dxz_a;
    abcd[iGrid*756+185] = 4.0E0*I_ESP_Kx4y2z_Dxz_aa-2.0E0*3*I_ESP_Hx2y2z_Dxz_a;
    abcd[iGrid*756+186] = 4.0E0*I_ESP_Kx3y3z_Dxz_aa-2.0E0*2*I_ESP_Hxy3z_Dxz_a;
    abcd[iGrid*756+187] = 4.0E0*I_ESP_Kx2y4z_Dxz_aa-2.0E0*1*I_ESP_Hx4z_Dxz_a;
    abcd[iGrid*756+188] = 4.0E0*I_ESP_Kxy5z_Dxz_aa;
    abcd[iGrid*756+189] = 4.0E0*I_ESP_K6xy_D2y_aa-2.0E0*5*I_ESP_H4xy_D2y_a;
    abcd[iGrid*756+190] = 4.0E0*I_ESP_K5x2y_D2y_aa-2.0E0*1*I_ESP_H5x_D2y_a-2.0E0*4*I_ESP_H3x2y_D2y_a+4*1*I_ESP_F3x_D2y;
    abcd[iGrid*756+191] = 4.0E0*I_ESP_K5xyz_D2y_aa-2.0E0*4*I_ESP_H3xyz_D2y_a;
    abcd[iGrid*756+192] = 4.0E0*I_ESP_K4x3y_D2y_aa-2.0E0*2*I_ESP_H4xy_D2y_a-2.0E0*3*I_ESP_H2x3y_D2y_a+3*2*I_ESP_F2xy_D2y;
    abcd[iGrid*756+193] = 4.0E0*I_ESP_K4x2yz_D2y_aa-2.0E0*1*I_ESP_H4xz_D2y_a-2.0E0*3*I_ESP_H2x2yz_D2y_a+3*1*I_ESP_F2xz_D2y;
    abcd[iGrid*756+194] = 4.0E0*I_ESP_K4xy2z_D2y_aa-2.0E0*3*I_ESP_H2xy2z_D2y_a;
    abcd[iGrid*756+195] = 4.0E0*I_ESP_K3x4y_D2y_aa-2.0E0*3*I_ESP_H3x2y_D2y_a-2.0E0*2*I_ESP_Hx4y_D2y_a+2*3*I_ESP_Fx2y_D2y;
    abcd[iGrid*756+196] = 4.0E0*I_ESP_K3x3yz_D2y_aa-2.0E0*2*I_ESP_H3xyz_D2y_a-2.0E0*2*I_ESP_Hx3yz_D2y_a+2*2*I_ESP_Fxyz_D2y;
    abcd[iGrid*756+197] = 4.0E0*I_ESP_K3x2y2z_D2y_aa-2.0E0*1*I_ESP_H3x2z_D2y_a-2.0E0*2*I_ESP_Hx2y2z_D2y_a+2*1*I_ESP_Fx2z_D2y;
    abcd[iGrid*756+198] = 4.0E0*I_ESP_K3xy3z_D2y_aa-2.0E0*2*I_ESP_Hxy3z_D2y_a;
    abcd[iGrid*756+199] = 4.0E0*I_ESP_K2x5y_D2y_aa-2.0E0*4*I_ESP_H2x3y_D2y_a-2.0E0*1*I_ESP_H5y_D2y_a+4*I_ESP_F3y_D2y;
    abcd[iGrid*756+200] = 4.0E0*I_ESP_K2x4yz_D2y_aa-2.0E0*3*I_ESP_H2x2yz_D2y_a-2.0E0*1*I_ESP_H4yz_D2y_a+3*I_ESP_F2yz_D2y;
    abcd[iGrid*756+201] = 4.0E0*I_ESP_K2x3y2z_D2y_aa-2.0E0*2*I_ESP_H2xy2z_D2y_a-2.0E0*1*I_ESP_H3y2z_D2y_a+2*I_ESP_Fy2z_D2y;
    abcd[iGrid*756+202] = 4.0E0*I_ESP_K2x2y3z_D2y_aa-2.0E0*1*I_ESP_H2x3z_D2y_a-2.0E0*1*I_ESP_H2y3z_D2y_a+1*I_ESP_F3z_D2y;
    abcd[iGrid*756+203] = 4.0E0*I_ESP_K2xy4z_D2y_aa-2.0E0*1*I_ESP_Hy4z_D2y_a;
    abcd[iGrid*756+204] = 4.0E0*I_ESP_Kx6y_D2y_aa-2.0E0*5*I_ESP_Hx4y_D2y_a;
    abcd[iGrid*756+205] = 4.0E0*I_ESP_Kx5yz_D2y_aa-2.0E0*4*I_ESP_Hx3yz_D2y_a;
    abcd[iGrid*756+206] = 4.0E0*I_ESP_Kx4y2z_D2y_aa-2.0E0*3*I_ESP_Hx2y2z_D2y_a;
    abcd[iGrid*756+207] = 4.0E0*I_ESP_Kx3y3z_D2y_aa-2.0E0*2*I_ESP_Hxy3z_D2y_a;
    abcd[iGrid*756+208] = 4.0E0*I_ESP_Kx2y4z_D2y_aa-2.0E0*1*I_ESP_Hx4z_D2y_a;
    abcd[iGrid*756+209] = 4.0E0*I_ESP_Kxy5z_D2y_aa;
    abcd[iGrid*756+210] = 4.0E0*I_ESP_K6xy_Dyz_aa-2.0E0*5*I_ESP_H4xy_Dyz_a;
    abcd[iGrid*756+211] = 4.0E0*I_ESP_K5x2y_Dyz_aa-2.0E0*1*I_ESP_H5x_Dyz_a-2.0E0*4*I_ESP_H3x2y_Dyz_a+4*1*I_ESP_F3x_Dyz;
    abcd[iGrid*756+212] = 4.0E0*I_ESP_K5xyz_Dyz_aa-2.0E0*4*I_ESP_H3xyz_Dyz_a;
    abcd[iGrid*756+213] = 4.0E0*I_ESP_K4x3y_Dyz_aa-2.0E0*2*I_ESP_H4xy_Dyz_a-2.0E0*3*I_ESP_H2x3y_Dyz_a+3*2*I_ESP_F2xy_Dyz;
    abcd[iGrid*756+214] = 4.0E0*I_ESP_K4x2yz_Dyz_aa-2.0E0*1*I_ESP_H4xz_Dyz_a-2.0E0*3*I_ESP_H2x2yz_Dyz_a+3*1*I_ESP_F2xz_Dyz;
    abcd[iGrid*756+215] = 4.0E0*I_ESP_K4xy2z_Dyz_aa-2.0E0*3*I_ESP_H2xy2z_Dyz_a;
    abcd[iGrid*756+216] = 4.0E0*I_ESP_K3x4y_Dyz_aa-2.0E0*3*I_ESP_H3x2y_Dyz_a-2.0E0*2*I_ESP_Hx4y_Dyz_a+2*3*I_ESP_Fx2y_Dyz;
    abcd[iGrid*756+217] = 4.0E0*I_ESP_K3x3yz_Dyz_aa-2.0E0*2*I_ESP_H3xyz_Dyz_a-2.0E0*2*I_ESP_Hx3yz_Dyz_a+2*2*I_ESP_Fxyz_Dyz;
    abcd[iGrid*756+218] = 4.0E0*I_ESP_K3x2y2z_Dyz_aa-2.0E0*1*I_ESP_H3x2z_Dyz_a-2.0E0*2*I_ESP_Hx2y2z_Dyz_a+2*1*I_ESP_Fx2z_Dyz;
    abcd[iGrid*756+219] = 4.0E0*I_ESP_K3xy3z_Dyz_aa-2.0E0*2*I_ESP_Hxy3z_Dyz_a;
    abcd[iGrid*756+220] = 4.0E0*I_ESP_K2x5y_Dyz_aa-2.0E0*4*I_ESP_H2x3y_Dyz_a-2.0E0*1*I_ESP_H5y_Dyz_a+4*I_ESP_F3y_Dyz;
    abcd[iGrid*756+221] = 4.0E0*I_ESP_K2x4yz_Dyz_aa-2.0E0*3*I_ESP_H2x2yz_Dyz_a-2.0E0*1*I_ESP_H4yz_Dyz_a+3*I_ESP_F2yz_Dyz;
    abcd[iGrid*756+222] = 4.0E0*I_ESP_K2x3y2z_Dyz_aa-2.0E0*2*I_ESP_H2xy2z_Dyz_a-2.0E0*1*I_ESP_H3y2z_Dyz_a+2*I_ESP_Fy2z_Dyz;
    abcd[iGrid*756+223] = 4.0E0*I_ESP_K2x2y3z_Dyz_aa-2.0E0*1*I_ESP_H2x3z_Dyz_a-2.0E0*1*I_ESP_H2y3z_Dyz_a+1*I_ESP_F3z_Dyz;
    abcd[iGrid*756+224] = 4.0E0*I_ESP_K2xy4z_Dyz_aa-2.0E0*1*I_ESP_Hy4z_Dyz_a;
    abcd[iGrid*756+225] = 4.0E0*I_ESP_Kx6y_Dyz_aa-2.0E0*5*I_ESP_Hx4y_Dyz_a;
    abcd[iGrid*756+226] = 4.0E0*I_ESP_Kx5yz_Dyz_aa-2.0E0*4*I_ESP_Hx3yz_Dyz_a;
    abcd[iGrid*756+227] = 4.0E0*I_ESP_Kx4y2z_Dyz_aa-2.0E0*3*I_ESP_Hx2y2z_Dyz_a;
    abcd[iGrid*756+228] = 4.0E0*I_ESP_Kx3y3z_Dyz_aa-2.0E0*2*I_ESP_Hxy3z_Dyz_a;
    abcd[iGrid*756+229] = 4.0E0*I_ESP_Kx2y4z_Dyz_aa-2.0E0*1*I_ESP_Hx4z_Dyz_a;
    abcd[iGrid*756+230] = 4.0E0*I_ESP_Kxy5z_Dyz_aa;
    abcd[iGrid*756+231] = 4.0E0*I_ESP_K6xy_D2z_aa-2.0E0*5*I_ESP_H4xy_D2z_a;
    abcd[iGrid*756+232] = 4.0E0*I_ESP_K5x2y_D2z_aa-2.0E0*1*I_ESP_H5x_D2z_a-2.0E0*4*I_ESP_H3x2y_D2z_a+4*1*I_ESP_F3x_D2z;
    abcd[iGrid*756+233] = 4.0E0*I_ESP_K5xyz_D2z_aa-2.0E0*4*I_ESP_H3xyz_D2z_a;
    abcd[iGrid*756+234] = 4.0E0*I_ESP_K4x3y_D2z_aa-2.0E0*2*I_ESP_H4xy_D2z_a-2.0E0*3*I_ESP_H2x3y_D2z_a+3*2*I_ESP_F2xy_D2z;
    abcd[iGrid*756+235] = 4.0E0*I_ESP_K4x2yz_D2z_aa-2.0E0*1*I_ESP_H4xz_D2z_a-2.0E0*3*I_ESP_H2x2yz_D2z_a+3*1*I_ESP_F2xz_D2z;
    abcd[iGrid*756+236] = 4.0E0*I_ESP_K4xy2z_D2z_aa-2.0E0*3*I_ESP_H2xy2z_D2z_a;
    abcd[iGrid*756+237] = 4.0E0*I_ESP_K3x4y_D2z_aa-2.0E0*3*I_ESP_H3x2y_D2z_a-2.0E0*2*I_ESP_Hx4y_D2z_a+2*3*I_ESP_Fx2y_D2z;
    abcd[iGrid*756+238] = 4.0E0*I_ESP_K3x3yz_D2z_aa-2.0E0*2*I_ESP_H3xyz_D2z_a-2.0E0*2*I_ESP_Hx3yz_D2z_a+2*2*I_ESP_Fxyz_D2z;
    abcd[iGrid*756+239] = 4.0E0*I_ESP_K3x2y2z_D2z_aa-2.0E0*1*I_ESP_H3x2z_D2z_a-2.0E0*2*I_ESP_Hx2y2z_D2z_a+2*1*I_ESP_Fx2z_D2z;
    abcd[iGrid*756+240] = 4.0E0*I_ESP_K3xy3z_D2z_aa-2.0E0*2*I_ESP_Hxy3z_D2z_a;
    abcd[iGrid*756+241] = 4.0E0*I_ESP_K2x5y_D2z_aa-2.0E0*4*I_ESP_H2x3y_D2z_a-2.0E0*1*I_ESP_H5y_D2z_a+4*I_ESP_F3y_D2z;
    abcd[iGrid*756+242] = 4.0E0*I_ESP_K2x4yz_D2z_aa-2.0E0*3*I_ESP_H2x2yz_D2z_a-2.0E0*1*I_ESP_H4yz_D2z_a+3*I_ESP_F2yz_D2z;
    abcd[iGrid*756+243] = 4.0E0*I_ESP_K2x3y2z_D2z_aa-2.0E0*2*I_ESP_H2xy2z_D2z_a-2.0E0*1*I_ESP_H3y2z_D2z_a+2*I_ESP_Fy2z_D2z;
    abcd[iGrid*756+244] = 4.0E0*I_ESP_K2x2y3z_D2z_aa-2.0E0*1*I_ESP_H2x3z_D2z_a-2.0E0*1*I_ESP_H2y3z_D2z_a+1*I_ESP_F3z_D2z;
    abcd[iGrid*756+245] = 4.0E0*I_ESP_K2xy4z_D2z_aa-2.0E0*1*I_ESP_Hy4z_D2z_a;
    abcd[iGrid*756+246] = 4.0E0*I_ESP_Kx6y_D2z_aa-2.0E0*5*I_ESP_Hx4y_D2z_a;
    abcd[iGrid*756+247] = 4.0E0*I_ESP_Kx5yz_D2z_aa-2.0E0*4*I_ESP_Hx3yz_D2z_a;
    abcd[iGrid*756+248] = 4.0E0*I_ESP_Kx4y2z_D2z_aa-2.0E0*3*I_ESP_Hx2y2z_D2z_a;
    abcd[iGrid*756+249] = 4.0E0*I_ESP_Kx3y3z_D2z_aa-2.0E0*2*I_ESP_Hxy3z_D2z_a;
    abcd[iGrid*756+250] = 4.0E0*I_ESP_Kx2y4z_D2z_aa-2.0E0*1*I_ESP_Hx4z_D2z_a;
    abcd[iGrid*756+251] = 4.0E0*I_ESP_Kxy5z_D2z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D_aa
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_F_D
     ************************************************************/
    abcd[iGrid*756+252] = 4.0E0*I_ESP_K6xz_D2x_aa-2.0E0*5*I_ESP_H4xz_D2x_a;
    abcd[iGrid*756+253] = 4.0E0*I_ESP_K5xyz_D2x_aa-2.0E0*4*I_ESP_H3xyz_D2x_a;
    abcd[iGrid*756+254] = 4.0E0*I_ESP_K5x2z_D2x_aa-2.0E0*1*I_ESP_H5x_D2x_a-2.0E0*4*I_ESP_H3x2z_D2x_a+4*1*I_ESP_F3x_D2x;
    abcd[iGrid*756+255] = 4.0E0*I_ESP_K4x2yz_D2x_aa-2.0E0*3*I_ESP_H2x2yz_D2x_a;
    abcd[iGrid*756+256] = 4.0E0*I_ESP_K4xy2z_D2x_aa-2.0E0*1*I_ESP_H4xy_D2x_a-2.0E0*3*I_ESP_H2xy2z_D2x_a+3*1*I_ESP_F2xy_D2x;
    abcd[iGrid*756+257] = 4.0E0*I_ESP_K4x3z_D2x_aa-2.0E0*2*I_ESP_H4xz_D2x_a-2.0E0*3*I_ESP_H2x3z_D2x_a+3*2*I_ESP_F2xz_D2x;
    abcd[iGrid*756+258] = 4.0E0*I_ESP_K3x3yz_D2x_aa-2.0E0*2*I_ESP_Hx3yz_D2x_a;
    abcd[iGrid*756+259] = 4.0E0*I_ESP_K3x2y2z_D2x_aa-2.0E0*1*I_ESP_H3x2y_D2x_a-2.0E0*2*I_ESP_Hx2y2z_D2x_a+2*1*I_ESP_Fx2y_D2x;
    abcd[iGrid*756+260] = 4.0E0*I_ESP_K3xy3z_D2x_aa-2.0E0*2*I_ESP_H3xyz_D2x_a-2.0E0*2*I_ESP_Hxy3z_D2x_a+2*2*I_ESP_Fxyz_D2x;
    abcd[iGrid*756+261] = 4.0E0*I_ESP_K3x4z_D2x_aa-2.0E0*3*I_ESP_H3x2z_D2x_a-2.0E0*2*I_ESP_Hx4z_D2x_a+2*3*I_ESP_Fx2z_D2x;
    abcd[iGrid*756+262] = 4.0E0*I_ESP_K2x4yz_D2x_aa-2.0E0*1*I_ESP_H4yz_D2x_a;
    abcd[iGrid*756+263] = 4.0E0*I_ESP_K2x3y2z_D2x_aa-2.0E0*1*I_ESP_H2x3y_D2x_a-2.0E0*1*I_ESP_H3y2z_D2x_a+1*I_ESP_F3y_D2x;
    abcd[iGrid*756+264] = 4.0E0*I_ESP_K2x2y3z_D2x_aa-2.0E0*2*I_ESP_H2x2yz_D2x_a-2.0E0*1*I_ESP_H2y3z_D2x_a+2*I_ESP_F2yz_D2x;
    abcd[iGrid*756+265] = 4.0E0*I_ESP_K2xy4z_D2x_aa-2.0E0*3*I_ESP_H2xy2z_D2x_a-2.0E0*1*I_ESP_Hy4z_D2x_a+3*I_ESP_Fy2z_D2x;
    abcd[iGrid*756+266] = 4.0E0*I_ESP_K2x5z_D2x_aa-2.0E0*4*I_ESP_H2x3z_D2x_a-2.0E0*1*I_ESP_H5z_D2x_a+4*I_ESP_F3z_D2x;
    abcd[iGrid*756+267] = 4.0E0*I_ESP_Kx5yz_D2x_aa;
    abcd[iGrid*756+268] = 4.0E0*I_ESP_Kx4y2z_D2x_aa-2.0E0*1*I_ESP_Hx4y_D2x_a;
    abcd[iGrid*756+269] = 4.0E0*I_ESP_Kx3y3z_D2x_aa-2.0E0*2*I_ESP_Hx3yz_D2x_a;
    abcd[iGrid*756+270] = 4.0E0*I_ESP_Kx2y4z_D2x_aa-2.0E0*3*I_ESP_Hx2y2z_D2x_a;
    abcd[iGrid*756+271] = 4.0E0*I_ESP_Kxy5z_D2x_aa-2.0E0*4*I_ESP_Hxy3z_D2x_a;
    abcd[iGrid*756+272] = 4.0E0*I_ESP_Kx6z_D2x_aa-2.0E0*5*I_ESP_Hx4z_D2x_a;
    abcd[iGrid*756+273] = 4.0E0*I_ESP_K6xz_Dxy_aa-2.0E0*5*I_ESP_H4xz_Dxy_a;
    abcd[iGrid*756+274] = 4.0E0*I_ESP_K5xyz_Dxy_aa-2.0E0*4*I_ESP_H3xyz_Dxy_a;
    abcd[iGrid*756+275] = 4.0E0*I_ESP_K5x2z_Dxy_aa-2.0E0*1*I_ESP_H5x_Dxy_a-2.0E0*4*I_ESP_H3x2z_Dxy_a+4*1*I_ESP_F3x_Dxy;
    abcd[iGrid*756+276] = 4.0E0*I_ESP_K4x2yz_Dxy_aa-2.0E0*3*I_ESP_H2x2yz_Dxy_a;
    abcd[iGrid*756+277] = 4.0E0*I_ESP_K4xy2z_Dxy_aa-2.0E0*1*I_ESP_H4xy_Dxy_a-2.0E0*3*I_ESP_H2xy2z_Dxy_a+3*1*I_ESP_F2xy_Dxy;
    abcd[iGrid*756+278] = 4.0E0*I_ESP_K4x3z_Dxy_aa-2.0E0*2*I_ESP_H4xz_Dxy_a-2.0E0*3*I_ESP_H2x3z_Dxy_a+3*2*I_ESP_F2xz_Dxy;
    abcd[iGrid*756+279] = 4.0E0*I_ESP_K3x3yz_Dxy_aa-2.0E0*2*I_ESP_Hx3yz_Dxy_a;
    abcd[iGrid*756+280] = 4.0E0*I_ESP_K3x2y2z_Dxy_aa-2.0E0*1*I_ESP_H3x2y_Dxy_a-2.0E0*2*I_ESP_Hx2y2z_Dxy_a+2*1*I_ESP_Fx2y_Dxy;
    abcd[iGrid*756+281] = 4.0E0*I_ESP_K3xy3z_Dxy_aa-2.0E0*2*I_ESP_H3xyz_Dxy_a-2.0E0*2*I_ESP_Hxy3z_Dxy_a+2*2*I_ESP_Fxyz_Dxy;
    abcd[iGrid*756+282] = 4.0E0*I_ESP_K3x4z_Dxy_aa-2.0E0*3*I_ESP_H3x2z_Dxy_a-2.0E0*2*I_ESP_Hx4z_Dxy_a+2*3*I_ESP_Fx2z_Dxy;
    abcd[iGrid*756+283] = 4.0E0*I_ESP_K2x4yz_Dxy_aa-2.0E0*1*I_ESP_H4yz_Dxy_a;
    abcd[iGrid*756+284] = 4.0E0*I_ESP_K2x3y2z_Dxy_aa-2.0E0*1*I_ESP_H2x3y_Dxy_a-2.0E0*1*I_ESP_H3y2z_Dxy_a+1*I_ESP_F3y_Dxy;
    abcd[iGrid*756+285] = 4.0E0*I_ESP_K2x2y3z_Dxy_aa-2.0E0*2*I_ESP_H2x2yz_Dxy_a-2.0E0*1*I_ESP_H2y3z_Dxy_a+2*I_ESP_F2yz_Dxy;
    abcd[iGrid*756+286] = 4.0E0*I_ESP_K2xy4z_Dxy_aa-2.0E0*3*I_ESP_H2xy2z_Dxy_a-2.0E0*1*I_ESP_Hy4z_Dxy_a+3*I_ESP_Fy2z_Dxy;
    abcd[iGrid*756+287] = 4.0E0*I_ESP_K2x5z_Dxy_aa-2.0E0*4*I_ESP_H2x3z_Dxy_a-2.0E0*1*I_ESP_H5z_Dxy_a+4*I_ESP_F3z_Dxy;
    abcd[iGrid*756+288] = 4.0E0*I_ESP_Kx5yz_Dxy_aa;
    abcd[iGrid*756+289] = 4.0E0*I_ESP_Kx4y2z_Dxy_aa-2.0E0*1*I_ESP_Hx4y_Dxy_a;
    abcd[iGrid*756+290] = 4.0E0*I_ESP_Kx3y3z_Dxy_aa-2.0E0*2*I_ESP_Hx3yz_Dxy_a;
    abcd[iGrid*756+291] = 4.0E0*I_ESP_Kx2y4z_Dxy_aa-2.0E0*3*I_ESP_Hx2y2z_Dxy_a;
    abcd[iGrid*756+292] = 4.0E0*I_ESP_Kxy5z_Dxy_aa-2.0E0*4*I_ESP_Hxy3z_Dxy_a;
    abcd[iGrid*756+293] = 4.0E0*I_ESP_Kx6z_Dxy_aa-2.0E0*5*I_ESP_Hx4z_Dxy_a;
    abcd[iGrid*756+294] = 4.0E0*I_ESP_K6xz_Dxz_aa-2.0E0*5*I_ESP_H4xz_Dxz_a;
    abcd[iGrid*756+295] = 4.0E0*I_ESP_K5xyz_Dxz_aa-2.0E0*4*I_ESP_H3xyz_Dxz_a;
    abcd[iGrid*756+296] = 4.0E0*I_ESP_K5x2z_Dxz_aa-2.0E0*1*I_ESP_H5x_Dxz_a-2.0E0*4*I_ESP_H3x2z_Dxz_a+4*1*I_ESP_F3x_Dxz;
    abcd[iGrid*756+297] = 4.0E0*I_ESP_K4x2yz_Dxz_aa-2.0E0*3*I_ESP_H2x2yz_Dxz_a;
    abcd[iGrid*756+298] = 4.0E0*I_ESP_K4xy2z_Dxz_aa-2.0E0*1*I_ESP_H4xy_Dxz_a-2.0E0*3*I_ESP_H2xy2z_Dxz_a+3*1*I_ESP_F2xy_Dxz;
    abcd[iGrid*756+299] = 4.0E0*I_ESP_K4x3z_Dxz_aa-2.0E0*2*I_ESP_H4xz_Dxz_a-2.0E0*3*I_ESP_H2x3z_Dxz_a+3*2*I_ESP_F2xz_Dxz;
    abcd[iGrid*756+300] = 4.0E0*I_ESP_K3x3yz_Dxz_aa-2.0E0*2*I_ESP_Hx3yz_Dxz_a;
    abcd[iGrid*756+301] = 4.0E0*I_ESP_K3x2y2z_Dxz_aa-2.0E0*1*I_ESP_H3x2y_Dxz_a-2.0E0*2*I_ESP_Hx2y2z_Dxz_a+2*1*I_ESP_Fx2y_Dxz;
    abcd[iGrid*756+302] = 4.0E0*I_ESP_K3xy3z_Dxz_aa-2.0E0*2*I_ESP_H3xyz_Dxz_a-2.0E0*2*I_ESP_Hxy3z_Dxz_a+2*2*I_ESP_Fxyz_Dxz;
    abcd[iGrid*756+303] = 4.0E0*I_ESP_K3x4z_Dxz_aa-2.0E0*3*I_ESP_H3x2z_Dxz_a-2.0E0*2*I_ESP_Hx4z_Dxz_a+2*3*I_ESP_Fx2z_Dxz;
    abcd[iGrid*756+304] = 4.0E0*I_ESP_K2x4yz_Dxz_aa-2.0E0*1*I_ESP_H4yz_Dxz_a;
    abcd[iGrid*756+305] = 4.0E0*I_ESP_K2x3y2z_Dxz_aa-2.0E0*1*I_ESP_H2x3y_Dxz_a-2.0E0*1*I_ESP_H3y2z_Dxz_a+1*I_ESP_F3y_Dxz;
    abcd[iGrid*756+306] = 4.0E0*I_ESP_K2x2y3z_Dxz_aa-2.0E0*2*I_ESP_H2x2yz_Dxz_a-2.0E0*1*I_ESP_H2y3z_Dxz_a+2*I_ESP_F2yz_Dxz;
    abcd[iGrid*756+307] = 4.0E0*I_ESP_K2xy4z_Dxz_aa-2.0E0*3*I_ESP_H2xy2z_Dxz_a-2.0E0*1*I_ESP_Hy4z_Dxz_a+3*I_ESP_Fy2z_Dxz;
    abcd[iGrid*756+308] = 4.0E0*I_ESP_K2x5z_Dxz_aa-2.0E0*4*I_ESP_H2x3z_Dxz_a-2.0E0*1*I_ESP_H5z_Dxz_a+4*I_ESP_F3z_Dxz;
    abcd[iGrid*756+309] = 4.0E0*I_ESP_Kx5yz_Dxz_aa;
    abcd[iGrid*756+310] = 4.0E0*I_ESP_Kx4y2z_Dxz_aa-2.0E0*1*I_ESP_Hx4y_Dxz_a;
    abcd[iGrid*756+311] = 4.0E0*I_ESP_Kx3y3z_Dxz_aa-2.0E0*2*I_ESP_Hx3yz_Dxz_a;
    abcd[iGrid*756+312] = 4.0E0*I_ESP_Kx2y4z_Dxz_aa-2.0E0*3*I_ESP_Hx2y2z_Dxz_a;
    abcd[iGrid*756+313] = 4.0E0*I_ESP_Kxy5z_Dxz_aa-2.0E0*4*I_ESP_Hxy3z_Dxz_a;
    abcd[iGrid*756+314] = 4.0E0*I_ESP_Kx6z_Dxz_aa-2.0E0*5*I_ESP_Hx4z_Dxz_a;
    abcd[iGrid*756+315] = 4.0E0*I_ESP_K6xz_D2y_aa-2.0E0*5*I_ESP_H4xz_D2y_a;
    abcd[iGrid*756+316] = 4.0E0*I_ESP_K5xyz_D2y_aa-2.0E0*4*I_ESP_H3xyz_D2y_a;
    abcd[iGrid*756+317] = 4.0E0*I_ESP_K5x2z_D2y_aa-2.0E0*1*I_ESP_H5x_D2y_a-2.0E0*4*I_ESP_H3x2z_D2y_a+4*1*I_ESP_F3x_D2y;
    abcd[iGrid*756+318] = 4.0E0*I_ESP_K4x2yz_D2y_aa-2.0E0*3*I_ESP_H2x2yz_D2y_a;
    abcd[iGrid*756+319] = 4.0E0*I_ESP_K4xy2z_D2y_aa-2.0E0*1*I_ESP_H4xy_D2y_a-2.0E0*3*I_ESP_H2xy2z_D2y_a+3*1*I_ESP_F2xy_D2y;
    abcd[iGrid*756+320] = 4.0E0*I_ESP_K4x3z_D2y_aa-2.0E0*2*I_ESP_H4xz_D2y_a-2.0E0*3*I_ESP_H2x3z_D2y_a+3*2*I_ESP_F2xz_D2y;
    abcd[iGrid*756+321] = 4.0E0*I_ESP_K3x3yz_D2y_aa-2.0E0*2*I_ESP_Hx3yz_D2y_a;
    abcd[iGrid*756+322] = 4.0E0*I_ESP_K3x2y2z_D2y_aa-2.0E0*1*I_ESP_H3x2y_D2y_a-2.0E0*2*I_ESP_Hx2y2z_D2y_a+2*1*I_ESP_Fx2y_D2y;
    abcd[iGrid*756+323] = 4.0E0*I_ESP_K3xy3z_D2y_aa-2.0E0*2*I_ESP_H3xyz_D2y_a-2.0E0*2*I_ESP_Hxy3z_D2y_a+2*2*I_ESP_Fxyz_D2y;
    abcd[iGrid*756+324] = 4.0E0*I_ESP_K3x4z_D2y_aa-2.0E0*3*I_ESP_H3x2z_D2y_a-2.0E0*2*I_ESP_Hx4z_D2y_a+2*3*I_ESP_Fx2z_D2y;
    abcd[iGrid*756+325] = 4.0E0*I_ESP_K2x4yz_D2y_aa-2.0E0*1*I_ESP_H4yz_D2y_a;
    abcd[iGrid*756+326] = 4.0E0*I_ESP_K2x3y2z_D2y_aa-2.0E0*1*I_ESP_H2x3y_D2y_a-2.0E0*1*I_ESP_H3y2z_D2y_a+1*I_ESP_F3y_D2y;
    abcd[iGrid*756+327] = 4.0E0*I_ESP_K2x2y3z_D2y_aa-2.0E0*2*I_ESP_H2x2yz_D2y_a-2.0E0*1*I_ESP_H2y3z_D2y_a+2*I_ESP_F2yz_D2y;
    abcd[iGrid*756+328] = 4.0E0*I_ESP_K2xy4z_D2y_aa-2.0E0*3*I_ESP_H2xy2z_D2y_a-2.0E0*1*I_ESP_Hy4z_D2y_a+3*I_ESP_Fy2z_D2y;
    abcd[iGrid*756+329] = 4.0E0*I_ESP_K2x5z_D2y_aa-2.0E0*4*I_ESP_H2x3z_D2y_a-2.0E0*1*I_ESP_H5z_D2y_a+4*I_ESP_F3z_D2y;
    abcd[iGrid*756+330] = 4.0E0*I_ESP_Kx5yz_D2y_aa;
    abcd[iGrid*756+331] = 4.0E0*I_ESP_Kx4y2z_D2y_aa-2.0E0*1*I_ESP_Hx4y_D2y_a;
    abcd[iGrid*756+332] = 4.0E0*I_ESP_Kx3y3z_D2y_aa-2.0E0*2*I_ESP_Hx3yz_D2y_a;
    abcd[iGrid*756+333] = 4.0E0*I_ESP_Kx2y4z_D2y_aa-2.0E0*3*I_ESP_Hx2y2z_D2y_a;
    abcd[iGrid*756+334] = 4.0E0*I_ESP_Kxy5z_D2y_aa-2.0E0*4*I_ESP_Hxy3z_D2y_a;
    abcd[iGrid*756+335] = 4.0E0*I_ESP_Kx6z_D2y_aa-2.0E0*5*I_ESP_Hx4z_D2y_a;
    abcd[iGrid*756+336] = 4.0E0*I_ESP_K6xz_Dyz_aa-2.0E0*5*I_ESP_H4xz_Dyz_a;
    abcd[iGrid*756+337] = 4.0E0*I_ESP_K5xyz_Dyz_aa-2.0E0*4*I_ESP_H3xyz_Dyz_a;
    abcd[iGrid*756+338] = 4.0E0*I_ESP_K5x2z_Dyz_aa-2.0E0*1*I_ESP_H5x_Dyz_a-2.0E0*4*I_ESP_H3x2z_Dyz_a+4*1*I_ESP_F3x_Dyz;
    abcd[iGrid*756+339] = 4.0E0*I_ESP_K4x2yz_Dyz_aa-2.0E0*3*I_ESP_H2x2yz_Dyz_a;
    abcd[iGrid*756+340] = 4.0E0*I_ESP_K4xy2z_Dyz_aa-2.0E0*1*I_ESP_H4xy_Dyz_a-2.0E0*3*I_ESP_H2xy2z_Dyz_a+3*1*I_ESP_F2xy_Dyz;
    abcd[iGrid*756+341] = 4.0E0*I_ESP_K4x3z_Dyz_aa-2.0E0*2*I_ESP_H4xz_Dyz_a-2.0E0*3*I_ESP_H2x3z_Dyz_a+3*2*I_ESP_F2xz_Dyz;
    abcd[iGrid*756+342] = 4.0E0*I_ESP_K3x3yz_Dyz_aa-2.0E0*2*I_ESP_Hx3yz_Dyz_a;
    abcd[iGrid*756+343] = 4.0E0*I_ESP_K3x2y2z_Dyz_aa-2.0E0*1*I_ESP_H3x2y_Dyz_a-2.0E0*2*I_ESP_Hx2y2z_Dyz_a+2*1*I_ESP_Fx2y_Dyz;
    abcd[iGrid*756+344] = 4.0E0*I_ESP_K3xy3z_Dyz_aa-2.0E0*2*I_ESP_H3xyz_Dyz_a-2.0E0*2*I_ESP_Hxy3z_Dyz_a+2*2*I_ESP_Fxyz_Dyz;
    abcd[iGrid*756+345] = 4.0E0*I_ESP_K3x4z_Dyz_aa-2.0E0*3*I_ESP_H3x2z_Dyz_a-2.0E0*2*I_ESP_Hx4z_Dyz_a+2*3*I_ESP_Fx2z_Dyz;
    abcd[iGrid*756+346] = 4.0E0*I_ESP_K2x4yz_Dyz_aa-2.0E0*1*I_ESP_H4yz_Dyz_a;
    abcd[iGrid*756+347] = 4.0E0*I_ESP_K2x3y2z_Dyz_aa-2.0E0*1*I_ESP_H2x3y_Dyz_a-2.0E0*1*I_ESP_H3y2z_Dyz_a+1*I_ESP_F3y_Dyz;
    abcd[iGrid*756+348] = 4.0E0*I_ESP_K2x2y3z_Dyz_aa-2.0E0*2*I_ESP_H2x2yz_Dyz_a-2.0E0*1*I_ESP_H2y3z_Dyz_a+2*I_ESP_F2yz_Dyz;
    abcd[iGrid*756+349] = 4.0E0*I_ESP_K2xy4z_Dyz_aa-2.0E0*3*I_ESP_H2xy2z_Dyz_a-2.0E0*1*I_ESP_Hy4z_Dyz_a+3*I_ESP_Fy2z_Dyz;
    abcd[iGrid*756+350] = 4.0E0*I_ESP_K2x5z_Dyz_aa-2.0E0*4*I_ESP_H2x3z_Dyz_a-2.0E0*1*I_ESP_H5z_Dyz_a+4*I_ESP_F3z_Dyz;
    abcd[iGrid*756+351] = 4.0E0*I_ESP_Kx5yz_Dyz_aa;
    abcd[iGrid*756+352] = 4.0E0*I_ESP_Kx4y2z_Dyz_aa-2.0E0*1*I_ESP_Hx4y_Dyz_a;
    abcd[iGrid*756+353] = 4.0E0*I_ESP_Kx3y3z_Dyz_aa-2.0E0*2*I_ESP_Hx3yz_Dyz_a;
    abcd[iGrid*756+354] = 4.0E0*I_ESP_Kx2y4z_Dyz_aa-2.0E0*3*I_ESP_Hx2y2z_Dyz_a;
    abcd[iGrid*756+355] = 4.0E0*I_ESP_Kxy5z_Dyz_aa-2.0E0*4*I_ESP_Hxy3z_Dyz_a;
    abcd[iGrid*756+356] = 4.0E0*I_ESP_Kx6z_Dyz_aa-2.0E0*5*I_ESP_Hx4z_Dyz_a;
    abcd[iGrid*756+357] = 4.0E0*I_ESP_K6xz_D2z_aa-2.0E0*5*I_ESP_H4xz_D2z_a;
    abcd[iGrid*756+358] = 4.0E0*I_ESP_K5xyz_D2z_aa-2.0E0*4*I_ESP_H3xyz_D2z_a;
    abcd[iGrid*756+359] = 4.0E0*I_ESP_K5x2z_D2z_aa-2.0E0*1*I_ESP_H5x_D2z_a-2.0E0*4*I_ESP_H3x2z_D2z_a+4*1*I_ESP_F3x_D2z;
    abcd[iGrid*756+360] = 4.0E0*I_ESP_K4x2yz_D2z_aa-2.0E0*3*I_ESP_H2x2yz_D2z_a;
    abcd[iGrid*756+361] = 4.0E0*I_ESP_K4xy2z_D2z_aa-2.0E0*1*I_ESP_H4xy_D2z_a-2.0E0*3*I_ESP_H2xy2z_D2z_a+3*1*I_ESP_F2xy_D2z;
    abcd[iGrid*756+362] = 4.0E0*I_ESP_K4x3z_D2z_aa-2.0E0*2*I_ESP_H4xz_D2z_a-2.0E0*3*I_ESP_H2x3z_D2z_a+3*2*I_ESP_F2xz_D2z;
    abcd[iGrid*756+363] = 4.0E0*I_ESP_K3x3yz_D2z_aa-2.0E0*2*I_ESP_Hx3yz_D2z_a;
    abcd[iGrid*756+364] = 4.0E0*I_ESP_K3x2y2z_D2z_aa-2.0E0*1*I_ESP_H3x2y_D2z_a-2.0E0*2*I_ESP_Hx2y2z_D2z_a+2*1*I_ESP_Fx2y_D2z;
    abcd[iGrid*756+365] = 4.0E0*I_ESP_K3xy3z_D2z_aa-2.0E0*2*I_ESP_H3xyz_D2z_a-2.0E0*2*I_ESP_Hxy3z_D2z_a+2*2*I_ESP_Fxyz_D2z;
    abcd[iGrid*756+366] = 4.0E0*I_ESP_K3x4z_D2z_aa-2.0E0*3*I_ESP_H3x2z_D2z_a-2.0E0*2*I_ESP_Hx4z_D2z_a+2*3*I_ESP_Fx2z_D2z;
    abcd[iGrid*756+367] = 4.0E0*I_ESP_K2x4yz_D2z_aa-2.0E0*1*I_ESP_H4yz_D2z_a;
    abcd[iGrid*756+368] = 4.0E0*I_ESP_K2x3y2z_D2z_aa-2.0E0*1*I_ESP_H2x3y_D2z_a-2.0E0*1*I_ESP_H3y2z_D2z_a+1*I_ESP_F3y_D2z;
    abcd[iGrid*756+369] = 4.0E0*I_ESP_K2x2y3z_D2z_aa-2.0E0*2*I_ESP_H2x2yz_D2z_a-2.0E0*1*I_ESP_H2y3z_D2z_a+2*I_ESP_F2yz_D2z;
    abcd[iGrid*756+370] = 4.0E0*I_ESP_K2xy4z_D2z_aa-2.0E0*3*I_ESP_H2xy2z_D2z_a-2.0E0*1*I_ESP_Hy4z_D2z_a+3*I_ESP_Fy2z_D2z;
    abcd[iGrid*756+371] = 4.0E0*I_ESP_K2x5z_D2z_aa-2.0E0*4*I_ESP_H2x3z_D2z_a-2.0E0*1*I_ESP_H5z_D2z_a+4*I_ESP_F3z_D2z;
    abcd[iGrid*756+372] = 4.0E0*I_ESP_Kx5yz_D2z_aa;
    abcd[iGrid*756+373] = 4.0E0*I_ESP_Kx4y2z_D2z_aa-2.0E0*1*I_ESP_Hx4y_D2z_a;
    abcd[iGrid*756+374] = 4.0E0*I_ESP_Kx3y3z_D2z_aa-2.0E0*2*I_ESP_Hx3yz_D2z_a;
    abcd[iGrid*756+375] = 4.0E0*I_ESP_Kx2y4z_D2z_aa-2.0E0*3*I_ESP_Hx2y2z_D2z_a;
    abcd[iGrid*756+376] = 4.0E0*I_ESP_Kxy5z_D2z_aa-2.0E0*4*I_ESP_Hxy3z_D2z_a;
    abcd[iGrid*756+377] = 4.0E0*I_ESP_Kx6z_D2z_aa-2.0E0*5*I_ESP_Hx4z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D_aa
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_F_D
     ************************************************************/
    abcd[iGrid*756+378] = 4.0E0*I_ESP_K5x2y_D2x_aa-2.0E0*1*I_ESP_H5x_D2x_a;
    abcd[iGrid*756+379] = 4.0E0*I_ESP_K4x3y_D2x_aa-2.0E0*1*I_ESP_H4xy_D2x_a-2.0E0*2*I_ESP_H4xy_D2x_a;
    abcd[iGrid*756+380] = 4.0E0*I_ESP_K4x2yz_D2x_aa-2.0E0*1*I_ESP_H4xz_D2x_a;
    abcd[iGrid*756+381] = 4.0E0*I_ESP_K3x4y_D2x_aa-2.0E0*2*I_ESP_H3x2y_D2x_a-2.0E0*3*I_ESP_H3x2y_D2x_a+2*1*I_ESP_F3x_D2x;
    abcd[iGrid*756+382] = 4.0E0*I_ESP_K3x3yz_D2x_aa-2.0E0*1*I_ESP_H3xyz_D2x_a-2.0E0*2*I_ESP_H3xyz_D2x_a;
    abcd[iGrid*756+383] = 4.0E0*I_ESP_K3x2y2z_D2x_aa-2.0E0*1*I_ESP_H3x2z_D2x_a;
    abcd[iGrid*756+384] = 4.0E0*I_ESP_K2x5y_D2x_aa-2.0E0*3*I_ESP_H2x3y_D2x_a-2.0E0*4*I_ESP_H2x3y_D2x_a+3*2*I_ESP_F2xy_D2x;
    abcd[iGrid*756+385] = 4.0E0*I_ESP_K2x4yz_D2x_aa-2.0E0*2*I_ESP_H2x2yz_D2x_a-2.0E0*3*I_ESP_H2x2yz_D2x_a+2*1*I_ESP_F2xz_D2x;
    abcd[iGrid*756+386] = 4.0E0*I_ESP_K2x3y2z_D2x_aa-2.0E0*1*I_ESP_H2xy2z_D2x_a-2.0E0*2*I_ESP_H2xy2z_D2x_a;
    abcd[iGrid*756+387] = 4.0E0*I_ESP_K2x2y3z_D2x_aa-2.0E0*1*I_ESP_H2x3z_D2x_a;
    abcd[iGrid*756+388] = 4.0E0*I_ESP_Kx6y_D2x_aa-2.0E0*4*I_ESP_Hx4y_D2x_a-2.0E0*5*I_ESP_Hx4y_D2x_a+4*3*I_ESP_Fx2y_D2x;
    abcd[iGrid*756+389] = 4.0E0*I_ESP_Kx5yz_D2x_aa-2.0E0*3*I_ESP_Hx3yz_D2x_a-2.0E0*4*I_ESP_Hx3yz_D2x_a+3*2*I_ESP_Fxyz_D2x;
    abcd[iGrid*756+390] = 4.0E0*I_ESP_Kx4y2z_D2x_aa-2.0E0*2*I_ESP_Hx2y2z_D2x_a-2.0E0*3*I_ESP_Hx2y2z_D2x_a+2*1*I_ESP_Fx2z_D2x;
    abcd[iGrid*756+391] = 4.0E0*I_ESP_Kx3y3z_D2x_aa-2.0E0*1*I_ESP_Hxy3z_D2x_a-2.0E0*2*I_ESP_Hxy3z_D2x_a;
    abcd[iGrid*756+392] = 4.0E0*I_ESP_Kx2y4z_D2x_aa-2.0E0*1*I_ESP_Hx4z_D2x_a;
    abcd[iGrid*756+393] = 4.0E0*I_ESP_K7y_D2x_aa-2.0E0*5*I_ESP_H5y_D2x_a-2.0E0*6*I_ESP_H5y_D2x_a+5*4*I_ESP_F3y_D2x;
    abcd[iGrid*756+394] = 4.0E0*I_ESP_K6yz_D2x_aa-2.0E0*4*I_ESP_H4yz_D2x_a-2.0E0*5*I_ESP_H4yz_D2x_a+4*3*I_ESP_F2yz_D2x;
    abcd[iGrid*756+395] = 4.0E0*I_ESP_K5y2z_D2x_aa-2.0E0*3*I_ESP_H3y2z_D2x_a-2.0E0*4*I_ESP_H3y2z_D2x_a+3*2*I_ESP_Fy2z_D2x;
    abcd[iGrid*756+396] = 4.0E0*I_ESP_K4y3z_D2x_aa-2.0E0*2*I_ESP_H2y3z_D2x_a-2.0E0*3*I_ESP_H2y3z_D2x_a+2*1*I_ESP_F3z_D2x;
    abcd[iGrid*756+397] = 4.0E0*I_ESP_K3y4z_D2x_aa-2.0E0*1*I_ESP_Hy4z_D2x_a-2.0E0*2*I_ESP_Hy4z_D2x_a;
    abcd[iGrid*756+398] = 4.0E0*I_ESP_K2y5z_D2x_aa-2.0E0*1*I_ESP_H5z_D2x_a;
    abcd[iGrid*756+399] = 4.0E0*I_ESP_K5x2y_Dxy_aa-2.0E0*1*I_ESP_H5x_Dxy_a;
    abcd[iGrid*756+400] = 4.0E0*I_ESP_K4x3y_Dxy_aa-2.0E0*1*I_ESP_H4xy_Dxy_a-2.0E0*2*I_ESP_H4xy_Dxy_a;
    abcd[iGrid*756+401] = 4.0E0*I_ESP_K4x2yz_Dxy_aa-2.0E0*1*I_ESP_H4xz_Dxy_a;
    abcd[iGrid*756+402] = 4.0E0*I_ESP_K3x4y_Dxy_aa-2.0E0*2*I_ESP_H3x2y_Dxy_a-2.0E0*3*I_ESP_H3x2y_Dxy_a+2*1*I_ESP_F3x_Dxy;
    abcd[iGrid*756+403] = 4.0E0*I_ESP_K3x3yz_Dxy_aa-2.0E0*1*I_ESP_H3xyz_Dxy_a-2.0E0*2*I_ESP_H3xyz_Dxy_a;
    abcd[iGrid*756+404] = 4.0E0*I_ESP_K3x2y2z_Dxy_aa-2.0E0*1*I_ESP_H3x2z_Dxy_a;
    abcd[iGrid*756+405] = 4.0E0*I_ESP_K2x5y_Dxy_aa-2.0E0*3*I_ESP_H2x3y_Dxy_a-2.0E0*4*I_ESP_H2x3y_Dxy_a+3*2*I_ESP_F2xy_Dxy;
    abcd[iGrid*756+406] = 4.0E0*I_ESP_K2x4yz_Dxy_aa-2.0E0*2*I_ESP_H2x2yz_Dxy_a-2.0E0*3*I_ESP_H2x2yz_Dxy_a+2*1*I_ESP_F2xz_Dxy;
    abcd[iGrid*756+407] = 4.0E0*I_ESP_K2x3y2z_Dxy_aa-2.0E0*1*I_ESP_H2xy2z_Dxy_a-2.0E0*2*I_ESP_H2xy2z_Dxy_a;
    abcd[iGrid*756+408] = 4.0E0*I_ESP_K2x2y3z_Dxy_aa-2.0E0*1*I_ESP_H2x3z_Dxy_a;
    abcd[iGrid*756+409] = 4.0E0*I_ESP_Kx6y_Dxy_aa-2.0E0*4*I_ESP_Hx4y_Dxy_a-2.0E0*5*I_ESP_Hx4y_Dxy_a+4*3*I_ESP_Fx2y_Dxy;
    abcd[iGrid*756+410] = 4.0E0*I_ESP_Kx5yz_Dxy_aa-2.0E0*3*I_ESP_Hx3yz_Dxy_a-2.0E0*4*I_ESP_Hx3yz_Dxy_a+3*2*I_ESP_Fxyz_Dxy;
    abcd[iGrid*756+411] = 4.0E0*I_ESP_Kx4y2z_Dxy_aa-2.0E0*2*I_ESP_Hx2y2z_Dxy_a-2.0E0*3*I_ESP_Hx2y2z_Dxy_a+2*1*I_ESP_Fx2z_Dxy;
    abcd[iGrid*756+412] = 4.0E0*I_ESP_Kx3y3z_Dxy_aa-2.0E0*1*I_ESP_Hxy3z_Dxy_a-2.0E0*2*I_ESP_Hxy3z_Dxy_a;
    abcd[iGrid*756+413] = 4.0E0*I_ESP_Kx2y4z_Dxy_aa-2.0E0*1*I_ESP_Hx4z_Dxy_a;
    abcd[iGrid*756+414] = 4.0E0*I_ESP_K7y_Dxy_aa-2.0E0*5*I_ESP_H5y_Dxy_a-2.0E0*6*I_ESP_H5y_Dxy_a+5*4*I_ESP_F3y_Dxy;
    abcd[iGrid*756+415] = 4.0E0*I_ESP_K6yz_Dxy_aa-2.0E0*4*I_ESP_H4yz_Dxy_a-2.0E0*5*I_ESP_H4yz_Dxy_a+4*3*I_ESP_F2yz_Dxy;
    abcd[iGrid*756+416] = 4.0E0*I_ESP_K5y2z_Dxy_aa-2.0E0*3*I_ESP_H3y2z_Dxy_a-2.0E0*4*I_ESP_H3y2z_Dxy_a+3*2*I_ESP_Fy2z_Dxy;
    abcd[iGrid*756+417] = 4.0E0*I_ESP_K4y3z_Dxy_aa-2.0E0*2*I_ESP_H2y3z_Dxy_a-2.0E0*3*I_ESP_H2y3z_Dxy_a+2*1*I_ESP_F3z_Dxy;
    abcd[iGrid*756+418] = 4.0E0*I_ESP_K3y4z_Dxy_aa-2.0E0*1*I_ESP_Hy4z_Dxy_a-2.0E0*2*I_ESP_Hy4z_Dxy_a;
    abcd[iGrid*756+419] = 4.0E0*I_ESP_K2y5z_Dxy_aa-2.0E0*1*I_ESP_H5z_Dxy_a;
    abcd[iGrid*756+420] = 4.0E0*I_ESP_K5x2y_Dxz_aa-2.0E0*1*I_ESP_H5x_Dxz_a;
    abcd[iGrid*756+421] = 4.0E0*I_ESP_K4x3y_Dxz_aa-2.0E0*1*I_ESP_H4xy_Dxz_a-2.0E0*2*I_ESP_H4xy_Dxz_a;
    abcd[iGrid*756+422] = 4.0E0*I_ESP_K4x2yz_Dxz_aa-2.0E0*1*I_ESP_H4xz_Dxz_a;
    abcd[iGrid*756+423] = 4.0E0*I_ESP_K3x4y_Dxz_aa-2.0E0*2*I_ESP_H3x2y_Dxz_a-2.0E0*3*I_ESP_H3x2y_Dxz_a+2*1*I_ESP_F3x_Dxz;
    abcd[iGrid*756+424] = 4.0E0*I_ESP_K3x3yz_Dxz_aa-2.0E0*1*I_ESP_H3xyz_Dxz_a-2.0E0*2*I_ESP_H3xyz_Dxz_a;
    abcd[iGrid*756+425] = 4.0E0*I_ESP_K3x2y2z_Dxz_aa-2.0E0*1*I_ESP_H3x2z_Dxz_a;
    abcd[iGrid*756+426] = 4.0E0*I_ESP_K2x5y_Dxz_aa-2.0E0*3*I_ESP_H2x3y_Dxz_a-2.0E0*4*I_ESP_H2x3y_Dxz_a+3*2*I_ESP_F2xy_Dxz;
    abcd[iGrid*756+427] = 4.0E0*I_ESP_K2x4yz_Dxz_aa-2.0E0*2*I_ESP_H2x2yz_Dxz_a-2.0E0*3*I_ESP_H2x2yz_Dxz_a+2*1*I_ESP_F2xz_Dxz;
    abcd[iGrid*756+428] = 4.0E0*I_ESP_K2x3y2z_Dxz_aa-2.0E0*1*I_ESP_H2xy2z_Dxz_a-2.0E0*2*I_ESP_H2xy2z_Dxz_a;
    abcd[iGrid*756+429] = 4.0E0*I_ESP_K2x2y3z_Dxz_aa-2.0E0*1*I_ESP_H2x3z_Dxz_a;
    abcd[iGrid*756+430] = 4.0E0*I_ESP_Kx6y_Dxz_aa-2.0E0*4*I_ESP_Hx4y_Dxz_a-2.0E0*5*I_ESP_Hx4y_Dxz_a+4*3*I_ESP_Fx2y_Dxz;
    abcd[iGrid*756+431] = 4.0E0*I_ESP_Kx5yz_Dxz_aa-2.0E0*3*I_ESP_Hx3yz_Dxz_a-2.0E0*4*I_ESP_Hx3yz_Dxz_a+3*2*I_ESP_Fxyz_Dxz;
    abcd[iGrid*756+432] = 4.0E0*I_ESP_Kx4y2z_Dxz_aa-2.0E0*2*I_ESP_Hx2y2z_Dxz_a-2.0E0*3*I_ESP_Hx2y2z_Dxz_a+2*1*I_ESP_Fx2z_Dxz;
    abcd[iGrid*756+433] = 4.0E0*I_ESP_Kx3y3z_Dxz_aa-2.0E0*1*I_ESP_Hxy3z_Dxz_a-2.0E0*2*I_ESP_Hxy3z_Dxz_a;
    abcd[iGrid*756+434] = 4.0E0*I_ESP_Kx2y4z_Dxz_aa-2.0E0*1*I_ESP_Hx4z_Dxz_a;
    abcd[iGrid*756+435] = 4.0E0*I_ESP_K7y_Dxz_aa-2.0E0*5*I_ESP_H5y_Dxz_a-2.0E0*6*I_ESP_H5y_Dxz_a+5*4*I_ESP_F3y_Dxz;
    abcd[iGrid*756+436] = 4.0E0*I_ESP_K6yz_Dxz_aa-2.0E0*4*I_ESP_H4yz_Dxz_a-2.0E0*5*I_ESP_H4yz_Dxz_a+4*3*I_ESP_F2yz_Dxz;
    abcd[iGrid*756+437] = 4.0E0*I_ESP_K5y2z_Dxz_aa-2.0E0*3*I_ESP_H3y2z_Dxz_a-2.0E0*4*I_ESP_H3y2z_Dxz_a+3*2*I_ESP_Fy2z_Dxz;
    abcd[iGrid*756+438] = 4.0E0*I_ESP_K4y3z_Dxz_aa-2.0E0*2*I_ESP_H2y3z_Dxz_a-2.0E0*3*I_ESP_H2y3z_Dxz_a+2*1*I_ESP_F3z_Dxz;
    abcd[iGrid*756+439] = 4.0E0*I_ESP_K3y4z_Dxz_aa-2.0E0*1*I_ESP_Hy4z_Dxz_a-2.0E0*2*I_ESP_Hy4z_Dxz_a;
    abcd[iGrid*756+440] = 4.0E0*I_ESP_K2y5z_Dxz_aa-2.0E0*1*I_ESP_H5z_Dxz_a;
    abcd[iGrid*756+441] = 4.0E0*I_ESP_K5x2y_D2y_aa-2.0E0*1*I_ESP_H5x_D2y_a;
    abcd[iGrid*756+442] = 4.0E0*I_ESP_K4x3y_D2y_aa-2.0E0*1*I_ESP_H4xy_D2y_a-2.0E0*2*I_ESP_H4xy_D2y_a;
    abcd[iGrid*756+443] = 4.0E0*I_ESP_K4x2yz_D2y_aa-2.0E0*1*I_ESP_H4xz_D2y_a;
    abcd[iGrid*756+444] = 4.0E0*I_ESP_K3x4y_D2y_aa-2.0E0*2*I_ESP_H3x2y_D2y_a-2.0E0*3*I_ESP_H3x2y_D2y_a+2*1*I_ESP_F3x_D2y;
    abcd[iGrid*756+445] = 4.0E0*I_ESP_K3x3yz_D2y_aa-2.0E0*1*I_ESP_H3xyz_D2y_a-2.0E0*2*I_ESP_H3xyz_D2y_a;
    abcd[iGrid*756+446] = 4.0E0*I_ESP_K3x2y2z_D2y_aa-2.0E0*1*I_ESP_H3x2z_D2y_a;
    abcd[iGrid*756+447] = 4.0E0*I_ESP_K2x5y_D2y_aa-2.0E0*3*I_ESP_H2x3y_D2y_a-2.0E0*4*I_ESP_H2x3y_D2y_a+3*2*I_ESP_F2xy_D2y;
    abcd[iGrid*756+448] = 4.0E0*I_ESP_K2x4yz_D2y_aa-2.0E0*2*I_ESP_H2x2yz_D2y_a-2.0E0*3*I_ESP_H2x2yz_D2y_a+2*1*I_ESP_F2xz_D2y;
    abcd[iGrid*756+449] = 4.0E0*I_ESP_K2x3y2z_D2y_aa-2.0E0*1*I_ESP_H2xy2z_D2y_a-2.0E0*2*I_ESP_H2xy2z_D2y_a;
    abcd[iGrid*756+450] = 4.0E0*I_ESP_K2x2y3z_D2y_aa-2.0E0*1*I_ESP_H2x3z_D2y_a;
    abcd[iGrid*756+451] = 4.0E0*I_ESP_Kx6y_D2y_aa-2.0E0*4*I_ESP_Hx4y_D2y_a-2.0E0*5*I_ESP_Hx4y_D2y_a+4*3*I_ESP_Fx2y_D2y;
    abcd[iGrid*756+452] = 4.0E0*I_ESP_Kx5yz_D2y_aa-2.0E0*3*I_ESP_Hx3yz_D2y_a-2.0E0*4*I_ESP_Hx3yz_D2y_a+3*2*I_ESP_Fxyz_D2y;
    abcd[iGrid*756+453] = 4.0E0*I_ESP_Kx4y2z_D2y_aa-2.0E0*2*I_ESP_Hx2y2z_D2y_a-2.0E0*3*I_ESP_Hx2y2z_D2y_a+2*1*I_ESP_Fx2z_D2y;
    abcd[iGrid*756+454] = 4.0E0*I_ESP_Kx3y3z_D2y_aa-2.0E0*1*I_ESP_Hxy3z_D2y_a-2.0E0*2*I_ESP_Hxy3z_D2y_a;
    abcd[iGrid*756+455] = 4.0E0*I_ESP_Kx2y4z_D2y_aa-2.0E0*1*I_ESP_Hx4z_D2y_a;
    abcd[iGrid*756+456] = 4.0E0*I_ESP_K7y_D2y_aa-2.0E0*5*I_ESP_H5y_D2y_a-2.0E0*6*I_ESP_H5y_D2y_a+5*4*I_ESP_F3y_D2y;
    abcd[iGrid*756+457] = 4.0E0*I_ESP_K6yz_D2y_aa-2.0E0*4*I_ESP_H4yz_D2y_a-2.0E0*5*I_ESP_H4yz_D2y_a+4*3*I_ESP_F2yz_D2y;
    abcd[iGrid*756+458] = 4.0E0*I_ESP_K5y2z_D2y_aa-2.0E0*3*I_ESP_H3y2z_D2y_a-2.0E0*4*I_ESP_H3y2z_D2y_a+3*2*I_ESP_Fy2z_D2y;
    abcd[iGrid*756+459] = 4.0E0*I_ESP_K4y3z_D2y_aa-2.0E0*2*I_ESP_H2y3z_D2y_a-2.0E0*3*I_ESP_H2y3z_D2y_a+2*1*I_ESP_F3z_D2y;
    abcd[iGrid*756+460] = 4.0E0*I_ESP_K3y4z_D2y_aa-2.0E0*1*I_ESP_Hy4z_D2y_a-2.0E0*2*I_ESP_Hy4z_D2y_a;
    abcd[iGrid*756+461] = 4.0E0*I_ESP_K2y5z_D2y_aa-2.0E0*1*I_ESP_H5z_D2y_a;
    abcd[iGrid*756+462] = 4.0E0*I_ESP_K5x2y_Dyz_aa-2.0E0*1*I_ESP_H5x_Dyz_a;
    abcd[iGrid*756+463] = 4.0E0*I_ESP_K4x3y_Dyz_aa-2.0E0*1*I_ESP_H4xy_Dyz_a-2.0E0*2*I_ESP_H4xy_Dyz_a;
    abcd[iGrid*756+464] = 4.0E0*I_ESP_K4x2yz_Dyz_aa-2.0E0*1*I_ESP_H4xz_Dyz_a;
    abcd[iGrid*756+465] = 4.0E0*I_ESP_K3x4y_Dyz_aa-2.0E0*2*I_ESP_H3x2y_Dyz_a-2.0E0*3*I_ESP_H3x2y_Dyz_a+2*1*I_ESP_F3x_Dyz;
    abcd[iGrid*756+466] = 4.0E0*I_ESP_K3x3yz_Dyz_aa-2.0E0*1*I_ESP_H3xyz_Dyz_a-2.0E0*2*I_ESP_H3xyz_Dyz_a;
    abcd[iGrid*756+467] = 4.0E0*I_ESP_K3x2y2z_Dyz_aa-2.0E0*1*I_ESP_H3x2z_Dyz_a;
    abcd[iGrid*756+468] = 4.0E0*I_ESP_K2x5y_Dyz_aa-2.0E0*3*I_ESP_H2x3y_Dyz_a-2.0E0*4*I_ESP_H2x3y_Dyz_a+3*2*I_ESP_F2xy_Dyz;
    abcd[iGrid*756+469] = 4.0E0*I_ESP_K2x4yz_Dyz_aa-2.0E0*2*I_ESP_H2x2yz_Dyz_a-2.0E0*3*I_ESP_H2x2yz_Dyz_a+2*1*I_ESP_F2xz_Dyz;
    abcd[iGrid*756+470] = 4.0E0*I_ESP_K2x3y2z_Dyz_aa-2.0E0*1*I_ESP_H2xy2z_Dyz_a-2.0E0*2*I_ESP_H2xy2z_Dyz_a;
    abcd[iGrid*756+471] = 4.0E0*I_ESP_K2x2y3z_Dyz_aa-2.0E0*1*I_ESP_H2x3z_Dyz_a;
    abcd[iGrid*756+472] = 4.0E0*I_ESP_Kx6y_Dyz_aa-2.0E0*4*I_ESP_Hx4y_Dyz_a-2.0E0*5*I_ESP_Hx4y_Dyz_a+4*3*I_ESP_Fx2y_Dyz;
    abcd[iGrid*756+473] = 4.0E0*I_ESP_Kx5yz_Dyz_aa-2.0E0*3*I_ESP_Hx3yz_Dyz_a-2.0E0*4*I_ESP_Hx3yz_Dyz_a+3*2*I_ESP_Fxyz_Dyz;
    abcd[iGrid*756+474] = 4.0E0*I_ESP_Kx4y2z_Dyz_aa-2.0E0*2*I_ESP_Hx2y2z_Dyz_a-2.0E0*3*I_ESP_Hx2y2z_Dyz_a+2*1*I_ESP_Fx2z_Dyz;
    abcd[iGrid*756+475] = 4.0E0*I_ESP_Kx3y3z_Dyz_aa-2.0E0*1*I_ESP_Hxy3z_Dyz_a-2.0E0*2*I_ESP_Hxy3z_Dyz_a;
    abcd[iGrid*756+476] = 4.0E0*I_ESP_Kx2y4z_Dyz_aa-2.0E0*1*I_ESP_Hx4z_Dyz_a;
    abcd[iGrid*756+477] = 4.0E0*I_ESP_K7y_Dyz_aa-2.0E0*5*I_ESP_H5y_Dyz_a-2.0E0*6*I_ESP_H5y_Dyz_a+5*4*I_ESP_F3y_Dyz;
    abcd[iGrid*756+478] = 4.0E0*I_ESP_K6yz_Dyz_aa-2.0E0*4*I_ESP_H4yz_Dyz_a-2.0E0*5*I_ESP_H4yz_Dyz_a+4*3*I_ESP_F2yz_Dyz;
    abcd[iGrid*756+479] = 4.0E0*I_ESP_K5y2z_Dyz_aa-2.0E0*3*I_ESP_H3y2z_Dyz_a-2.0E0*4*I_ESP_H3y2z_Dyz_a+3*2*I_ESP_Fy2z_Dyz;
    abcd[iGrid*756+480] = 4.0E0*I_ESP_K4y3z_Dyz_aa-2.0E0*2*I_ESP_H2y3z_Dyz_a-2.0E0*3*I_ESP_H2y3z_Dyz_a+2*1*I_ESP_F3z_Dyz;
    abcd[iGrid*756+481] = 4.0E0*I_ESP_K3y4z_Dyz_aa-2.0E0*1*I_ESP_Hy4z_Dyz_a-2.0E0*2*I_ESP_Hy4z_Dyz_a;
    abcd[iGrid*756+482] = 4.0E0*I_ESP_K2y5z_Dyz_aa-2.0E0*1*I_ESP_H5z_Dyz_a;
    abcd[iGrid*756+483] = 4.0E0*I_ESP_K5x2y_D2z_aa-2.0E0*1*I_ESP_H5x_D2z_a;
    abcd[iGrid*756+484] = 4.0E0*I_ESP_K4x3y_D2z_aa-2.0E0*1*I_ESP_H4xy_D2z_a-2.0E0*2*I_ESP_H4xy_D2z_a;
    abcd[iGrid*756+485] = 4.0E0*I_ESP_K4x2yz_D2z_aa-2.0E0*1*I_ESP_H4xz_D2z_a;
    abcd[iGrid*756+486] = 4.0E0*I_ESP_K3x4y_D2z_aa-2.0E0*2*I_ESP_H3x2y_D2z_a-2.0E0*3*I_ESP_H3x2y_D2z_a+2*1*I_ESP_F3x_D2z;
    abcd[iGrid*756+487] = 4.0E0*I_ESP_K3x3yz_D2z_aa-2.0E0*1*I_ESP_H3xyz_D2z_a-2.0E0*2*I_ESP_H3xyz_D2z_a;
    abcd[iGrid*756+488] = 4.0E0*I_ESP_K3x2y2z_D2z_aa-2.0E0*1*I_ESP_H3x2z_D2z_a;
    abcd[iGrid*756+489] = 4.0E0*I_ESP_K2x5y_D2z_aa-2.0E0*3*I_ESP_H2x3y_D2z_a-2.0E0*4*I_ESP_H2x3y_D2z_a+3*2*I_ESP_F2xy_D2z;
    abcd[iGrid*756+490] = 4.0E0*I_ESP_K2x4yz_D2z_aa-2.0E0*2*I_ESP_H2x2yz_D2z_a-2.0E0*3*I_ESP_H2x2yz_D2z_a+2*1*I_ESP_F2xz_D2z;
    abcd[iGrid*756+491] = 4.0E0*I_ESP_K2x3y2z_D2z_aa-2.0E0*1*I_ESP_H2xy2z_D2z_a-2.0E0*2*I_ESP_H2xy2z_D2z_a;
    abcd[iGrid*756+492] = 4.0E0*I_ESP_K2x2y3z_D2z_aa-2.0E0*1*I_ESP_H2x3z_D2z_a;
    abcd[iGrid*756+493] = 4.0E0*I_ESP_Kx6y_D2z_aa-2.0E0*4*I_ESP_Hx4y_D2z_a-2.0E0*5*I_ESP_Hx4y_D2z_a+4*3*I_ESP_Fx2y_D2z;
    abcd[iGrid*756+494] = 4.0E0*I_ESP_Kx5yz_D2z_aa-2.0E0*3*I_ESP_Hx3yz_D2z_a-2.0E0*4*I_ESP_Hx3yz_D2z_a+3*2*I_ESP_Fxyz_D2z;
    abcd[iGrid*756+495] = 4.0E0*I_ESP_Kx4y2z_D2z_aa-2.0E0*2*I_ESP_Hx2y2z_D2z_a-2.0E0*3*I_ESP_Hx2y2z_D2z_a+2*1*I_ESP_Fx2z_D2z;
    abcd[iGrid*756+496] = 4.0E0*I_ESP_Kx3y3z_D2z_aa-2.0E0*1*I_ESP_Hxy3z_D2z_a-2.0E0*2*I_ESP_Hxy3z_D2z_a;
    abcd[iGrid*756+497] = 4.0E0*I_ESP_Kx2y4z_D2z_aa-2.0E0*1*I_ESP_Hx4z_D2z_a;
    abcd[iGrid*756+498] = 4.0E0*I_ESP_K7y_D2z_aa-2.0E0*5*I_ESP_H5y_D2z_a-2.0E0*6*I_ESP_H5y_D2z_a+5*4*I_ESP_F3y_D2z;
    abcd[iGrid*756+499] = 4.0E0*I_ESP_K6yz_D2z_aa-2.0E0*4*I_ESP_H4yz_D2z_a-2.0E0*5*I_ESP_H4yz_D2z_a+4*3*I_ESP_F2yz_D2z;
    abcd[iGrid*756+500] = 4.0E0*I_ESP_K5y2z_D2z_aa-2.0E0*3*I_ESP_H3y2z_D2z_a-2.0E0*4*I_ESP_H3y2z_D2z_a+3*2*I_ESP_Fy2z_D2z;
    abcd[iGrid*756+501] = 4.0E0*I_ESP_K4y3z_D2z_aa-2.0E0*2*I_ESP_H2y3z_D2z_a-2.0E0*3*I_ESP_H2y3z_D2z_a+2*1*I_ESP_F3z_D2z;
    abcd[iGrid*756+502] = 4.0E0*I_ESP_K3y4z_D2z_aa-2.0E0*1*I_ESP_Hy4z_D2z_a-2.0E0*2*I_ESP_Hy4z_D2z_a;
    abcd[iGrid*756+503] = 4.0E0*I_ESP_K2y5z_D2z_aa-2.0E0*1*I_ESP_H5z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D_aa
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_F_D
     ************************************************************/
    abcd[iGrid*756+504] = 4.0E0*I_ESP_K5xyz_D2x_aa;
    abcd[iGrid*756+505] = 4.0E0*I_ESP_K4x2yz_D2x_aa-2.0E0*1*I_ESP_H4xz_D2x_a;
    abcd[iGrid*756+506] = 4.0E0*I_ESP_K4xy2z_D2x_aa-2.0E0*1*I_ESP_H4xy_D2x_a;
    abcd[iGrid*756+507] = 4.0E0*I_ESP_K3x3yz_D2x_aa-2.0E0*2*I_ESP_H3xyz_D2x_a;
    abcd[iGrid*756+508] = 4.0E0*I_ESP_K3x2y2z_D2x_aa-2.0E0*1*I_ESP_H3x2y_D2x_a-2.0E0*1*I_ESP_H3x2z_D2x_a+1*I_ESP_F3x_D2x;
    abcd[iGrid*756+509] = 4.0E0*I_ESP_K3xy3z_D2x_aa-2.0E0*2*I_ESP_H3xyz_D2x_a;
    abcd[iGrid*756+510] = 4.0E0*I_ESP_K2x4yz_D2x_aa-2.0E0*3*I_ESP_H2x2yz_D2x_a;
    abcd[iGrid*756+511] = 4.0E0*I_ESP_K2x3y2z_D2x_aa-2.0E0*1*I_ESP_H2x3y_D2x_a-2.0E0*2*I_ESP_H2xy2z_D2x_a+2*1*I_ESP_F2xy_D2x;
    abcd[iGrid*756+512] = 4.0E0*I_ESP_K2x2y3z_D2x_aa-2.0E0*2*I_ESP_H2x2yz_D2x_a-2.0E0*1*I_ESP_H2x3z_D2x_a+2*I_ESP_F2xz_D2x;
    abcd[iGrid*756+513] = 4.0E0*I_ESP_K2xy4z_D2x_aa-2.0E0*3*I_ESP_H2xy2z_D2x_a;
    abcd[iGrid*756+514] = 4.0E0*I_ESP_Kx5yz_D2x_aa-2.0E0*4*I_ESP_Hx3yz_D2x_a;
    abcd[iGrid*756+515] = 4.0E0*I_ESP_Kx4y2z_D2x_aa-2.0E0*1*I_ESP_Hx4y_D2x_a-2.0E0*3*I_ESP_Hx2y2z_D2x_a+3*1*I_ESP_Fx2y_D2x;
    abcd[iGrid*756+516] = 4.0E0*I_ESP_Kx3y3z_D2x_aa-2.0E0*2*I_ESP_Hx3yz_D2x_a-2.0E0*2*I_ESP_Hxy3z_D2x_a+2*2*I_ESP_Fxyz_D2x;
    abcd[iGrid*756+517] = 4.0E0*I_ESP_Kx2y4z_D2x_aa-2.0E0*3*I_ESP_Hx2y2z_D2x_a-2.0E0*1*I_ESP_Hx4z_D2x_a+3*I_ESP_Fx2z_D2x;
    abcd[iGrid*756+518] = 4.0E0*I_ESP_Kxy5z_D2x_aa-2.0E0*4*I_ESP_Hxy3z_D2x_a;
    abcd[iGrid*756+519] = 4.0E0*I_ESP_K6yz_D2x_aa-2.0E0*5*I_ESP_H4yz_D2x_a;
    abcd[iGrid*756+520] = 4.0E0*I_ESP_K5y2z_D2x_aa-2.0E0*1*I_ESP_H5y_D2x_a-2.0E0*4*I_ESP_H3y2z_D2x_a+4*1*I_ESP_F3y_D2x;
    abcd[iGrid*756+521] = 4.0E0*I_ESP_K4y3z_D2x_aa-2.0E0*2*I_ESP_H4yz_D2x_a-2.0E0*3*I_ESP_H2y3z_D2x_a+3*2*I_ESP_F2yz_D2x;
    abcd[iGrid*756+522] = 4.0E0*I_ESP_K3y4z_D2x_aa-2.0E0*3*I_ESP_H3y2z_D2x_a-2.0E0*2*I_ESP_Hy4z_D2x_a+2*3*I_ESP_Fy2z_D2x;
    abcd[iGrid*756+523] = 4.0E0*I_ESP_K2y5z_D2x_aa-2.0E0*4*I_ESP_H2y3z_D2x_a-2.0E0*1*I_ESP_H5z_D2x_a+4*I_ESP_F3z_D2x;
    abcd[iGrid*756+524] = 4.0E0*I_ESP_Ky6z_D2x_aa-2.0E0*5*I_ESP_Hy4z_D2x_a;
    abcd[iGrid*756+525] = 4.0E0*I_ESP_K5xyz_Dxy_aa;
    abcd[iGrid*756+526] = 4.0E0*I_ESP_K4x2yz_Dxy_aa-2.0E0*1*I_ESP_H4xz_Dxy_a;
    abcd[iGrid*756+527] = 4.0E0*I_ESP_K4xy2z_Dxy_aa-2.0E0*1*I_ESP_H4xy_Dxy_a;
    abcd[iGrid*756+528] = 4.0E0*I_ESP_K3x3yz_Dxy_aa-2.0E0*2*I_ESP_H3xyz_Dxy_a;
    abcd[iGrid*756+529] = 4.0E0*I_ESP_K3x2y2z_Dxy_aa-2.0E0*1*I_ESP_H3x2y_Dxy_a-2.0E0*1*I_ESP_H3x2z_Dxy_a+1*I_ESP_F3x_Dxy;
    abcd[iGrid*756+530] = 4.0E0*I_ESP_K3xy3z_Dxy_aa-2.0E0*2*I_ESP_H3xyz_Dxy_a;
    abcd[iGrid*756+531] = 4.0E0*I_ESP_K2x4yz_Dxy_aa-2.0E0*3*I_ESP_H2x2yz_Dxy_a;
    abcd[iGrid*756+532] = 4.0E0*I_ESP_K2x3y2z_Dxy_aa-2.0E0*1*I_ESP_H2x3y_Dxy_a-2.0E0*2*I_ESP_H2xy2z_Dxy_a+2*1*I_ESP_F2xy_Dxy;
    abcd[iGrid*756+533] = 4.0E0*I_ESP_K2x2y3z_Dxy_aa-2.0E0*2*I_ESP_H2x2yz_Dxy_a-2.0E0*1*I_ESP_H2x3z_Dxy_a+2*I_ESP_F2xz_Dxy;
    abcd[iGrid*756+534] = 4.0E0*I_ESP_K2xy4z_Dxy_aa-2.0E0*3*I_ESP_H2xy2z_Dxy_a;
    abcd[iGrid*756+535] = 4.0E0*I_ESP_Kx5yz_Dxy_aa-2.0E0*4*I_ESP_Hx3yz_Dxy_a;
    abcd[iGrid*756+536] = 4.0E0*I_ESP_Kx4y2z_Dxy_aa-2.0E0*1*I_ESP_Hx4y_Dxy_a-2.0E0*3*I_ESP_Hx2y2z_Dxy_a+3*1*I_ESP_Fx2y_Dxy;
    abcd[iGrid*756+537] = 4.0E0*I_ESP_Kx3y3z_Dxy_aa-2.0E0*2*I_ESP_Hx3yz_Dxy_a-2.0E0*2*I_ESP_Hxy3z_Dxy_a+2*2*I_ESP_Fxyz_Dxy;
    abcd[iGrid*756+538] = 4.0E0*I_ESP_Kx2y4z_Dxy_aa-2.0E0*3*I_ESP_Hx2y2z_Dxy_a-2.0E0*1*I_ESP_Hx4z_Dxy_a+3*I_ESP_Fx2z_Dxy;
    abcd[iGrid*756+539] = 4.0E0*I_ESP_Kxy5z_Dxy_aa-2.0E0*4*I_ESP_Hxy3z_Dxy_a;
    abcd[iGrid*756+540] = 4.0E0*I_ESP_K6yz_Dxy_aa-2.0E0*5*I_ESP_H4yz_Dxy_a;
    abcd[iGrid*756+541] = 4.0E0*I_ESP_K5y2z_Dxy_aa-2.0E0*1*I_ESP_H5y_Dxy_a-2.0E0*4*I_ESP_H3y2z_Dxy_a+4*1*I_ESP_F3y_Dxy;
    abcd[iGrid*756+542] = 4.0E0*I_ESP_K4y3z_Dxy_aa-2.0E0*2*I_ESP_H4yz_Dxy_a-2.0E0*3*I_ESP_H2y3z_Dxy_a+3*2*I_ESP_F2yz_Dxy;
    abcd[iGrid*756+543] = 4.0E0*I_ESP_K3y4z_Dxy_aa-2.0E0*3*I_ESP_H3y2z_Dxy_a-2.0E0*2*I_ESP_Hy4z_Dxy_a+2*3*I_ESP_Fy2z_Dxy;
    abcd[iGrid*756+544] = 4.0E0*I_ESP_K2y5z_Dxy_aa-2.0E0*4*I_ESP_H2y3z_Dxy_a-2.0E0*1*I_ESP_H5z_Dxy_a+4*I_ESP_F3z_Dxy;
    abcd[iGrid*756+545] = 4.0E0*I_ESP_Ky6z_Dxy_aa-2.0E0*5*I_ESP_Hy4z_Dxy_a;
    abcd[iGrid*756+546] = 4.0E0*I_ESP_K5xyz_Dxz_aa;
    abcd[iGrid*756+547] = 4.0E0*I_ESP_K4x2yz_Dxz_aa-2.0E0*1*I_ESP_H4xz_Dxz_a;
    abcd[iGrid*756+548] = 4.0E0*I_ESP_K4xy2z_Dxz_aa-2.0E0*1*I_ESP_H4xy_Dxz_a;
    abcd[iGrid*756+549] = 4.0E0*I_ESP_K3x3yz_Dxz_aa-2.0E0*2*I_ESP_H3xyz_Dxz_a;
    abcd[iGrid*756+550] = 4.0E0*I_ESP_K3x2y2z_Dxz_aa-2.0E0*1*I_ESP_H3x2y_Dxz_a-2.0E0*1*I_ESP_H3x2z_Dxz_a+1*I_ESP_F3x_Dxz;
    abcd[iGrid*756+551] = 4.0E0*I_ESP_K3xy3z_Dxz_aa-2.0E0*2*I_ESP_H3xyz_Dxz_a;
    abcd[iGrid*756+552] = 4.0E0*I_ESP_K2x4yz_Dxz_aa-2.0E0*3*I_ESP_H2x2yz_Dxz_a;
    abcd[iGrid*756+553] = 4.0E0*I_ESP_K2x3y2z_Dxz_aa-2.0E0*1*I_ESP_H2x3y_Dxz_a-2.0E0*2*I_ESP_H2xy2z_Dxz_a+2*1*I_ESP_F2xy_Dxz;
    abcd[iGrid*756+554] = 4.0E0*I_ESP_K2x2y3z_Dxz_aa-2.0E0*2*I_ESP_H2x2yz_Dxz_a-2.0E0*1*I_ESP_H2x3z_Dxz_a+2*I_ESP_F2xz_Dxz;
    abcd[iGrid*756+555] = 4.0E0*I_ESP_K2xy4z_Dxz_aa-2.0E0*3*I_ESP_H2xy2z_Dxz_a;
    abcd[iGrid*756+556] = 4.0E0*I_ESP_Kx5yz_Dxz_aa-2.0E0*4*I_ESP_Hx3yz_Dxz_a;
    abcd[iGrid*756+557] = 4.0E0*I_ESP_Kx4y2z_Dxz_aa-2.0E0*1*I_ESP_Hx4y_Dxz_a-2.0E0*3*I_ESP_Hx2y2z_Dxz_a+3*1*I_ESP_Fx2y_Dxz;
    abcd[iGrid*756+558] = 4.0E0*I_ESP_Kx3y3z_Dxz_aa-2.0E0*2*I_ESP_Hx3yz_Dxz_a-2.0E0*2*I_ESP_Hxy3z_Dxz_a+2*2*I_ESP_Fxyz_Dxz;
    abcd[iGrid*756+559] = 4.0E0*I_ESP_Kx2y4z_Dxz_aa-2.0E0*3*I_ESP_Hx2y2z_Dxz_a-2.0E0*1*I_ESP_Hx4z_Dxz_a+3*I_ESP_Fx2z_Dxz;
    abcd[iGrid*756+560] = 4.0E0*I_ESP_Kxy5z_Dxz_aa-2.0E0*4*I_ESP_Hxy3z_Dxz_a;
    abcd[iGrid*756+561] = 4.0E0*I_ESP_K6yz_Dxz_aa-2.0E0*5*I_ESP_H4yz_Dxz_a;
    abcd[iGrid*756+562] = 4.0E0*I_ESP_K5y2z_Dxz_aa-2.0E0*1*I_ESP_H5y_Dxz_a-2.0E0*4*I_ESP_H3y2z_Dxz_a+4*1*I_ESP_F3y_Dxz;
    abcd[iGrid*756+563] = 4.0E0*I_ESP_K4y3z_Dxz_aa-2.0E0*2*I_ESP_H4yz_Dxz_a-2.0E0*3*I_ESP_H2y3z_Dxz_a+3*2*I_ESP_F2yz_Dxz;
    abcd[iGrid*756+564] = 4.0E0*I_ESP_K3y4z_Dxz_aa-2.0E0*3*I_ESP_H3y2z_Dxz_a-2.0E0*2*I_ESP_Hy4z_Dxz_a+2*3*I_ESP_Fy2z_Dxz;
    abcd[iGrid*756+565] = 4.0E0*I_ESP_K2y5z_Dxz_aa-2.0E0*4*I_ESP_H2y3z_Dxz_a-2.0E0*1*I_ESP_H5z_Dxz_a+4*I_ESP_F3z_Dxz;
    abcd[iGrid*756+566] = 4.0E0*I_ESP_Ky6z_Dxz_aa-2.0E0*5*I_ESP_Hy4z_Dxz_a;
    abcd[iGrid*756+567] = 4.0E0*I_ESP_K5xyz_D2y_aa;
    abcd[iGrid*756+568] = 4.0E0*I_ESP_K4x2yz_D2y_aa-2.0E0*1*I_ESP_H4xz_D2y_a;
    abcd[iGrid*756+569] = 4.0E0*I_ESP_K4xy2z_D2y_aa-2.0E0*1*I_ESP_H4xy_D2y_a;
    abcd[iGrid*756+570] = 4.0E0*I_ESP_K3x3yz_D2y_aa-2.0E0*2*I_ESP_H3xyz_D2y_a;
    abcd[iGrid*756+571] = 4.0E0*I_ESP_K3x2y2z_D2y_aa-2.0E0*1*I_ESP_H3x2y_D2y_a-2.0E0*1*I_ESP_H3x2z_D2y_a+1*I_ESP_F3x_D2y;
    abcd[iGrid*756+572] = 4.0E0*I_ESP_K3xy3z_D2y_aa-2.0E0*2*I_ESP_H3xyz_D2y_a;
    abcd[iGrid*756+573] = 4.0E0*I_ESP_K2x4yz_D2y_aa-2.0E0*3*I_ESP_H2x2yz_D2y_a;
    abcd[iGrid*756+574] = 4.0E0*I_ESP_K2x3y2z_D2y_aa-2.0E0*1*I_ESP_H2x3y_D2y_a-2.0E0*2*I_ESP_H2xy2z_D2y_a+2*1*I_ESP_F2xy_D2y;
    abcd[iGrid*756+575] = 4.0E0*I_ESP_K2x2y3z_D2y_aa-2.0E0*2*I_ESP_H2x2yz_D2y_a-2.0E0*1*I_ESP_H2x3z_D2y_a+2*I_ESP_F2xz_D2y;
    abcd[iGrid*756+576] = 4.0E0*I_ESP_K2xy4z_D2y_aa-2.0E0*3*I_ESP_H2xy2z_D2y_a;
    abcd[iGrid*756+577] = 4.0E0*I_ESP_Kx5yz_D2y_aa-2.0E0*4*I_ESP_Hx3yz_D2y_a;
    abcd[iGrid*756+578] = 4.0E0*I_ESP_Kx4y2z_D2y_aa-2.0E0*1*I_ESP_Hx4y_D2y_a-2.0E0*3*I_ESP_Hx2y2z_D2y_a+3*1*I_ESP_Fx2y_D2y;
    abcd[iGrid*756+579] = 4.0E0*I_ESP_Kx3y3z_D2y_aa-2.0E0*2*I_ESP_Hx3yz_D2y_a-2.0E0*2*I_ESP_Hxy3z_D2y_a+2*2*I_ESP_Fxyz_D2y;
    abcd[iGrid*756+580] = 4.0E0*I_ESP_Kx2y4z_D2y_aa-2.0E0*3*I_ESP_Hx2y2z_D2y_a-2.0E0*1*I_ESP_Hx4z_D2y_a+3*I_ESP_Fx2z_D2y;
    abcd[iGrid*756+581] = 4.0E0*I_ESP_Kxy5z_D2y_aa-2.0E0*4*I_ESP_Hxy3z_D2y_a;
    abcd[iGrid*756+582] = 4.0E0*I_ESP_K6yz_D2y_aa-2.0E0*5*I_ESP_H4yz_D2y_a;
    abcd[iGrid*756+583] = 4.0E0*I_ESP_K5y2z_D2y_aa-2.0E0*1*I_ESP_H5y_D2y_a-2.0E0*4*I_ESP_H3y2z_D2y_a+4*1*I_ESP_F3y_D2y;
    abcd[iGrid*756+584] = 4.0E0*I_ESP_K4y3z_D2y_aa-2.0E0*2*I_ESP_H4yz_D2y_a-2.0E0*3*I_ESP_H2y3z_D2y_a+3*2*I_ESP_F2yz_D2y;
    abcd[iGrid*756+585] = 4.0E0*I_ESP_K3y4z_D2y_aa-2.0E0*3*I_ESP_H3y2z_D2y_a-2.0E0*2*I_ESP_Hy4z_D2y_a+2*3*I_ESP_Fy2z_D2y;
    abcd[iGrid*756+586] = 4.0E0*I_ESP_K2y5z_D2y_aa-2.0E0*4*I_ESP_H2y3z_D2y_a-2.0E0*1*I_ESP_H5z_D2y_a+4*I_ESP_F3z_D2y;
    abcd[iGrid*756+587] = 4.0E0*I_ESP_Ky6z_D2y_aa-2.0E0*5*I_ESP_Hy4z_D2y_a;
    abcd[iGrid*756+588] = 4.0E0*I_ESP_K5xyz_Dyz_aa;
    abcd[iGrid*756+589] = 4.0E0*I_ESP_K4x2yz_Dyz_aa-2.0E0*1*I_ESP_H4xz_Dyz_a;
    abcd[iGrid*756+590] = 4.0E0*I_ESP_K4xy2z_Dyz_aa-2.0E0*1*I_ESP_H4xy_Dyz_a;
    abcd[iGrid*756+591] = 4.0E0*I_ESP_K3x3yz_Dyz_aa-2.0E0*2*I_ESP_H3xyz_Dyz_a;
    abcd[iGrid*756+592] = 4.0E0*I_ESP_K3x2y2z_Dyz_aa-2.0E0*1*I_ESP_H3x2y_Dyz_a-2.0E0*1*I_ESP_H3x2z_Dyz_a+1*I_ESP_F3x_Dyz;
    abcd[iGrid*756+593] = 4.0E0*I_ESP_K3xy3z_Dyz_aa-2.0E0*2*I_ESP_H3xyz_Dyz_a;
    abcd[iGrid*756+594] = 4.0E0*I_ESP_K2x4yz_Dyz_aa-2.0E0*3*I_ESP_H2x2yz_Dyz_a;
    abcd[iGrid*756+595] = 4.0E0*I_ESP_K2x3y2z_Dyz_aa-2.0E0*1*I_ESP_H2x3y_Dyz_a-2.0E0*2*I_ESP_H2xy2z_Dyz_a+2*1*I_ESP_F2xy_Dyz;
    abcd[iGrid*756+596] = 4.0E0*I_ESP_K2x2y3z_Dyz_aa-2.0E0*2*I_ESP_H2x2yz_Dyz_a-2.0E0*1*I_ESP_H2x3z_Dyz_a+2*I_ESP_F2xz_Dyz;
    abcd[iGrid*756+597] = 4.0E0*I_ESP_K2xy4z_Dyz_aa-2.0E0*3*I_ESP_H2xy2z_Dyz_a;
    abcd[iGrid*756+598] = 4.0E0*I_ESP_Kx5yz_Dyz_aa-2.0E0*4*I_ESP_Hx3yz_Dyz_a;
    abcd[iGrid*756+599] = 4.0E0*I_ESP_Kx4y2z_Dyz_aa-2.0E0*1*I_ESP_Hx4y_Dyz_a-2.0E0*3*I_ESP_Hx2y2z_Dyz_a+3*1*I_ESP_Fx2y_Dyz;
    abcd[iGrid*756+600] = 4.0E0*I_ESP_Kx3y3z_Dyz_aa-2.0E0*2*I_ESP_Hx3yz_Dyz_a-2.0E0*2*I_ESP_Hxy3z_Dyz_a+2*2*I_ESP_Fxyz_Dyz;
    abcd[iGrid*756+601] = 4.0E0*I_ESP_Kx2y4z_Dyz_aa-2.0E0*3*I_ESP_Hx2y2z_Dyz_a-2.0E0*1*I_ESP_Hx4z_Dyz_a+3*I_ESP_Fx2z_Dyz;
    abcd[iGrid*756+602] = 4.0E0*I_ESP_Kxy5z_Dyz_aa-2.0E0*4*I_ESP_Hxy3z_Dyz_a;
    abcd[iGrid*756+603] = 4.0E0*I_ESP_K6yz_Dyz_aa-2.0E0*5*I_ESP_H4yz_Dyz_a;
    abcd[iGrid*756+604] = 4.0E0*I_ESP_K5y2z_Dyz_aa-2.0E0*1*I_ESP_H5y_Dyz_a-2.0E0*4*I_ESP_H3y2z_Dyz_a+4*1*I_ESP_F3y_Dyz;
    abcd[iGrid*756+605] = 4.0E0*I_ESP_K4y3z_Dyz_aa-2.0E0*2*I_ESP_H4yz_Dyz_a-2.0E0*3*I_ESP_H2y3z_Dyz_a+3*2*I_ESP_F2yz_Dyz;
    abcd[iGrid*756+606] = 4.0E0*I_ESP_K3y4z_Dyz_aa-2.0E0*3*I_ESP_H3y2z_Dyz_a-2.0E0*2*I_ESP_Hy4z_Dyz_a+2*3*I_ESP_Fy2z_Dyz;
    abcd[iGrid*756+607] = 4.0E0*I_ESP_K2y5z_Dyz_aa-2.0E0*4*I_ESP_H2y3z_Dyz_a-2.0E0*1*I_ESP_H5z_Dyz_a+4*I_ESP_F3z_Dyz;
    abcd[iGrid*756+608] = 4.0E0*I_ESP_Ky6z_Dyz_aa-2.0E0*5*I_ESP_Hy4z_Dyz_a;
    abcd[iGrid*756+609] = 4.0E0*I_ESP_K5xyz_D2z_aa;
    abcd[iGrid*756+610] = 4.0E0*I_ESP_K4x2yz_D2z_aa-2.0E0*1*I_ESP_H4xz_D2z_a;
    abcd[iGrid*756+611] = 4.0E0*I_ESP_K4xy2z_D2z_aa-2.0E0*1*I_ESP_H4xy_D2z_a;
    abcd[iGrid*756+612] = 4.0E0*I_ESP_K3x3yz_D2z_aa-2.0E0*2*I_ESP_H3xyz_D2z_a;
    abcd[iGrid*756+613] = 4.0E0*I_ESP_K3x2y2z_D2z_aa-2.0E0*1*I_ESP_H3x2y_D2z_a-2.0E0*1*I_ESP_H3x2z_D2z_a+1*I_ESP_F3x_D2z;
    abcd[iGrid*756+614] = 4.0E0*I_ESP_K3xy3z_D2z_aa-2.0E0*2*I_ESP_H3xyz_D2z_a;
    abcd[iGrid*756+615] = 4.0E0*I_ESP_K2x4yz_D2z_aa-2.0E0*3*I_ESP_H2x2yz_D2z_a;
    abcd[iGrid*756+616] = 4.0E0*I_ESP_K2x3y2z_D2z_aa-2.0E0*1*I_ESP_H2x3y_D2z_a-2.0E0*2*I_ESP_H2xy2z_D2z_a+2*1*I_ESP_F2xy_D2z;
    abcd[iGrid*756+617] = 4.0E0*I_ESP_K2x2y3z_D2z_aa-2.0E0*2*I_ESP_H2x2yz_D2z_a-2.0E0*1*I_ESP_H2x3z_D2z_a+2*I_ESP_F2xz_D2z;
    abcd[iGrid*756+618] = 4.0E0*I_ESP_K2xy4z_D2z_aa-2.0E0*3*I_ESP_H2xy2z_D2z_a;
    abcd[iGrid*756+619] = 4.0E0*I_ESP_Kx5yz_D2z_aa-2.0E0*4*I_ESP_Hx3yz_D2z_a;
    abcd[iGrid*756+620] = 4.0E0*I_ESP_Kx4y2z_D2z_aa-2.0E0*1*I_ESP_Hx4y_D2z_a-2.0E0*3*I_ESP_Hx2y2z_D2z_a+3*1*I_ESP_Fx2y_D2z;
    abcd[iGrid*756+621] = 4.0E0*I_ESP_Kx3y3z_D2z_aa-2.0E0*2*I_ESP_Hx3yz_D2z_a-2.0E0*2*I_ESP_Hxy3z_D2z_a+2*2*I_ESP_Fxyz_D2z;
    abcd[iGrid*756+622] = 4.0E0*I_ESP_Kx2y4z_D2z_aa-2.0E0*3*I_ESP_Hx2y2z_D2z_a-2.0E0*1*I_ESP_Hx4z_D2z_a+3*I_ESP_Fx2z_D2z;
    abcd[iGrid*756+623] = 4.0E0*I_ESP_Kxy5z_D2z_aa-2.0E0*4*I_ESP_Hxy3z_D2z_a;
    abcd[iGrid*756+624] = 4.0E0*I_ESP_K6yz_D2z_aa-2.0E0*5*I_ESP_H4yz_D2z_a;
    abcd[iGrid*756+625] = 4.0E0*I_ESP_K5y2z_D2z_aa-2.0E0*1*I_ESP_H5y_D2z_a-2.0E0*4*I_ESP_H3y2z_D2z_a+4*1*I_ESP_F3y_D2z;
    abcd[iGrid*756+626] = 4.0E0*I_ESP_K4y3z_D2z_aa-2.0E0*2*I_ESP_H4yz_D2z_a-2.0E0*3*I_ESP_H2y3z_D2z_a+3*2*I_ESP_F2yz_D2z;
    abcd[iGrid*756+627] = 4.0E0*I_ESP_K3y4z_D2z_aa-2.0E0*3*I_ESP_H3y2z_D2z_a-2.0E0*2*I_ESP_Hy4z_D2z_a+2*3*I_ESP_Fy2z_D2z;
    abcd[iGrid*756+628] = 4.0E0*I_ESP_K2y5z_D2z_aa-2.0E0*4*I_ESP_H2y3z_D2z_a-2.0E0*1*I_ESP_H5z_D2z_a+4*I_ESP_F3z_D2z;
    abcd[iGrid*756+629] = 4.0E0*I_ESP_Ky6z_D2z_aa-2.0E0*5*I_ESP_Hy4z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D_aa
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_F_D
     ************************************************************/
    abcd[iGrid*756+630] = 4.0E0*I_ESP_K5x2z_D2x_aa-2.0E0*1*I_ESP_H5x_D2x_a;
    abcd[iGrid*756+631] = 4.0E0*I_ESP_K4xy2z_D2x_aa-2.0E0*1*I_ESP_H4xy_D2x_a;
    abcd[iGrid*756+632] = 4.0E0*I_ESP_K4x3z_D2x_aa-2.0E0*1*I_ESP_H4xz_D2x_a-2.0E0*2*I_ESP_H4xz_D2x_a;
    abcd[iGrid*756+633] = 4.0E0*I_ESP_K3x2y2z_D2x_aa-2.0E0*1*I_ESP_H3x2y_D2x_a;
    abcd[iGrid*756+634] = 4.0E0*I_ESP_K3xy3z_D2x_aa-2.0E0*1*I_ESP_H3xyz_D2x_a-2.0E0*2*I_ESP_H3xyz_D2x_a;
    abcd[iGrid*756+635] = 4.0E0*I_ESP_K3x4z_D2x_aa-2.0E0*2*I_ESP_H3x2z_D2x_a-2.0E0*3*I_ESP_H3x2z_D2x_a+2*1*I_ESP_F3x_D2x;
    abcd[iGrid*756+636] = 4.0E0*I_ESP_K2x3y2z_D2x_aa-2.0E0*1*I_ESP_H2x3y_D2x_a;
    abcd[iGrid*756+637] = 4.0E0*I_ESP_K2x2y3z_D2x_aa-2.0E0*1*I_ESP_H2x2yz_D2x_a-2.0E0*2*I_ESP_H2x2yz_D2x_a;
    abcd[iGrid*756+638] = 4.0E0*I_ESP_K2xy4z_D2x_aa-2.0E0*2*I_ESP_H2xy2z_D2x_a-2.0E0*3*I_ESP_H2xy2z_D2x_a+2*1*I_ESP_F2xy_D2x;
    abcd[iGrid*756+639] = 4.0E0*I_ESP_K2x5z_D2x_aa-2.0E0*3*I_ESP_H2x3z_D2x_a-2.0E0*4*I_ESP_H2x3z_D2x_a+3*2*I_ESP_F2xz_D2x;
    abcd[iGrid*756+640] = 4.0E0*I_ESP_Kx4y2z_D2x_aa-2.0E0*1*I_ESP_Hx4y_D2x_a;
    abcd[iGrid*756+641] = 4.0E0*I_ESP_Kx3y3z_D2x_aa-2.0E0*1*I_ESP_Hx3yz_D2x_a-2.0E0*2*I_ESP_Hx3yz_D2x_a;
    abcd[iGrid*756+642] = 4.0E0*I_ESP_Kx2y4z_D2x_aa-2.0E0*2*I_ESP_Hx2y2z_D2x_a-2.0E0*3*I_ESP_Hx2y2z_D2x_a+2*1*I_ESP_Fx2y_D2x;
    abcd[iGrid*756+643] = 4.0E0*I_ESP_Kxy5z_D2x_aa-2.0E0*3*I_ESP_Hxy3z_D2x_a-2.0E0*4*I_ESP_Hxy3z_D2x_a+3*2*I_ESP_Fxyz_D2x;
    abcd[iGrid*756+644] = 4.0E0*I_ESP_Kx6z_D2x_aa-2.0E0*4*I_ESP_Hx4z_D2x_a-2.0E0*5*I_ESP_Hx4z_D2x_a+4*3*I_ESP_Fx2z_D2x;
    abcd[iGrid*756+645] = 4.0E0*I_ESP_K5y2z_D2x_aa-2.0E0*1*I_ESP_H5y_D2x_a;
    abcd[iGrid*756+646] = 4.0E0*I_ESP_K4y3z_D2x_aa-2.0E0*1*I_ESP_H4yz_D2x_a-2.0E0*2*I_ESP_H4yz_D2x_a;
    abcd[iGrid*756+647] = 4.0E0*I_ESP_K3y4z_D2x_aa-2.0E0*2*I_ESP_H3y2z_D2x_a-2.0E0*3*I_ESP_H3y2z_D2x_a+2*1*I_ESP_F3y_D2x;
    abcd[iGrid*756+648] = 4.0E0*I_ESP_K2y5z_D2x_aa-2.0E0*3*I_ESP_H2y3z_D2x_a-2.0E0*4*I_ESP_H2y3z_D2x_a+3*2*I_ESP_F2yz_D2x;
    abcd[iGrid*756+649] = 4.0E0*I_ESP_Ky6z_D2x_aa-2.0E0*4*I_ESP_Hy4z_D2x_a-2.0E0*5*I_ESP_Hy4z_D2x_a+4*3*I_ESP_Fy2z_D2x;
    abcd[iGrid*756+650] = 4.0E0*I_ESP_K7z_D2x_aa-2.0E0*5*I_ESP_H5z_D2x_a-2.0E0*6*I_ESP_H5z_D2x_a+5*4*I_ESP_F3z_D2x;
    abcd[iGrid*756+651] = 4.0E0*I_ESP_K5x2z_Dxy_aa-2.0E0*1*I_ESP_H5x_Dxy_a;
    abcd[iGrid*756+652] = 4.0E0*I_ESP_K4xy2z_Dxy_aa-2.0E0*1*I_ESP_H4xy_Dxy_a;
    abcd[iGrid*756+653] = 4.0E0*I_ESP_K4x3z_Dxy_aa-2.0E0*1*I_ESP_H4xz_Dxy_a-2.0E0*2*I_ESP_H4xz_Dxy_a;
    abcd[iGrid*756+654] = 4.0E0*I_ESP_K3x2y2z_Dxy_aa-2.0E0*1*I_ESP_H3x2y_Dxy_a;
    abcd[iGrid*756+655] = 4.0E0*I_ESP_K3xy3z_Dxy_aa-2.0E0*1*I_ESP_H3xyz_Dxy_a-2.0E0*2*I_ESP_H3xyz_Dxy_a;
    abcd[iGrid*756+656] = 4.0E0*I_ESP_K3x4z_Dxy_aa-2.0E0*2*I_ESP_H3x2z_Dxy_a-2.0E0*3*I_ESP_H3x2z_Dxy_a+2*1*I_ESP_F3x_Dxy;
    abcd[iGrid*756+657] = 4.0E0*I_ESP_K2x3y2z_Dxy_aa-2.0E0*1*I_ESP_H2x3y_Dxy_a;
    abcd[iGrid*756+658] = 4.0E0*I_ESP_K2x2y3z_Dxy_aa-2.0E0*1*I_ESP_H2x2yz_Dxy_a-2.0E0*2*I_ESP_H2x2yz_Dxy_a;
    abcd[iGrid*756+659] = 4.0E0*I_ESP_K2xy4z_Dxy_aa-2.0E0*2*I_ESP_H2xy2z_Dxy_a-2.0E0*3*I_ESP_H2xy2z_Dxy_a+2*1*I_ESP_F2xy_Dxy;
    abcd[iGrid*756+660] = 4.0E0*I_ESP_K2x5z_Dxy_aa-2.0E0*3*I_ESP_H2x3z_Dxy_a-2.0E0*4*I_ESP_H2x3z_Dxy_a+3*2*I_ESP_F2xz_Dxy;
    abcd[iGrid*756+661] = 4.0E0*I_ESP_Kx4y2z_Dxy_aa-2.0E0*1*I_ESP_Hx4y_Dxy_a;
    abcd[iGrid*756+662] = 4.0E0*I_ESP_Kx3y3z_Dxy_aa-2.0E0*1*I_ESP_Hx3yz_Dxy_a-2.0E0*2*I_ESP_Hx3yz_Dxy_a;
    abcd[iGrid*756+663] = 4.0E0*I_ESP_Kx2y4z_Dxy_aa-2.0E0*2*I_ESP_Hx2y2z_Dxy_a-2.0E0*3*I_ESP_Hx2y2z_Dxy_a+2*1*I_ESP_Fx2y_Dxy;
    abcd[iGrid*756+664] = 4.0E0*I_ESP_Kxy5z_Dxy_aa-2.0E0*3*I_ESP_Hxy3z_Dxy_a-2.0E0*4*I_ESP_Hxy3z_Dxy_a+3*2*I_ESP_Fxyz_Dxy;
    abcd[iGrid*756+665] = 4.0E0*I_ESP_Kx6z_Dxy_aa-2.0E0*4*I_ESP_Hx4z_Dxy_a-2.0E0*5*I_ESP_Hx4z_Dxy_a+4*3*I_ESP_Fx2z_Dxy;
    abcd[iGrid*756+666] = 4.0E0*I_ESP_K5y2z_Dxy_aa-2.0E0*1*I_ESP_H5y_Dxy_a;
    abcd[iGrid*756+667] = 4.0E0*I_ESP_K4y3z_Dxy_aa-2.0E0*1*I_ESP_H4yz_Dxy_a-2.0E0*2*I_ESP_H4yz_Dxy_a;
    abcd[iGrid*756+668] = 4.0E0*I_ESP_K3y4z_Dxy_aa-2.0E0*2*I_ESP_H3y2z_Dxy_a-2.0E0*3*I_ESP_H3y2z_Dxy_a+2*1*I_ESP_F3y_Dxy;
    abcd[iGrid*756+669] = 4.0E0*I_ESP_K2y5z_Dxy_aa-2.0E0*3*I_ESP_H2y3z_Dxy_a-2.0E0*4*I_ESP_H2y3z_Dxy_a+3*2*I_ESP_F2yz_Dxy;
    abcd[iGrid*756+670] = 4.0E0*I_ESP_Ky6z_Dxy_aa-2.0E0*4*I_ESP_Hy4z_Dxy_a-2.0E0*5*I_ESP_Hy4z_Dxy_a+4*3*I_ESP_Fy2z_Dxy;
    abcd[iGrid*756+671] = 4.0E0*I_ESP_K7z_Dxy_aa-2.0E0*5*I_ESP_H5z_Dxy_a-2.0E0*6*I_ESP_H5z_Dxy_a+5*4*I_ESP_F3z_Dxy;
    abcd[iGrid*756+672] = 4.0E0*I_ESP_K5x2z_Dxz_aa-2.0E0*1*I_ESP_H5x_Dxz_a;
    abcd[iGrid*756+673] = 4.0E0*I_ESP_K4xy2z_Dxz_aa-2.0E0*1*I_ESP_H4xy_Dxz_a;
    abcd[iGrid*756+674] = 4.0E0*I_ESP_K4x3z_Dxz_aa-2.0E0*1*I_ESP_H4xz_Dxz_a-2.0E0*2*I_ESP_H4xz_Dxz_a;
    abcd[iGrid*756+675] = 4.0E0*I_ESP_K3x2y2z_Dxz_aa-2.0E0*1*I_ESP_H3x2y_Dxz_a;
    abcd[iGrid*756+676] = 4.0E0*I_ESP_K3xy3z_Dxz_aa-2.0E0*1*I_ESP_H3xyz_Dxz_a-2.0E0*2*I_ESP_H3xyz_Dxz_a;
    abcd[iGrid*756+677] = 4.0E0*I_ESP_K3x4z_Dxz_aa-2.0E0*2*I_ESP_H3x2z_Dxz_a-2.0E0*3*I_ESP_H3x2z_Dxz_a+2*1*I_ESP_F3x_Dxz;
    abcd[iGrid*756+678] = 4.0E0*I_ESP_K2x3y2z_Dxz_aa-2.0E0*1*I_ESP_H2x3y_Dxz_a;
    abcd[iGrid*756+679] = 4.0E0*I_ESP_K2x2y3z_Dxz_aa-2.0E0*1*I_ESP_H2x2yz_Dxz_a-2.0E0*2*I_ESP_H2x2yz_Dxz_a;
    abcd[iGrid*756+680] = 4.0E0*I_ESP_K2xy4z_Dxz_aa-2.0E0*2*I_ESP_H2xy2z_Dxz_a-2.0E0*3*I_ESP_H2xy2z_Dxz_a+2*1*I_ESP_F2xy_Dxz;
    abcd[iGrid*756+681] = 4.0E0*I_ESP_K2x5z_Dxz_aa-2.0E0*3*I_ESP_H2x3z_Dxz_a-2.0E0*4*I_ESP_H2x3z_Dxz_a+3*2*I_ESP_F2xz_Dxz;
    abcd[iGrid*756+682] = 4.0E0*I_ESP_Kx4y2z_Dxz_aa-2.0E0*1*I_ESP_Hx4y_Dxz_a;
    abcd[iGrid*756+683] = 4.0E0*I_ESP_Kx3y3z_Dxz_aa-2.0E0*1*I_ESP_Hx3yz_Dxz_a-2.0E0*2*I_ESP_Hx3yz_Dxz_a;
    abcd[iGrid*756+684] = 4.0E0*I_ESP_Kx2y4z_Dxz_aa-2.0E0*2*I_ESP_Hx2y2z_Dxz_a-2.0E0*3*I_ESP_Hx2y2z_Dxz_a+2*1*I_ESP_Fx2y_Dxz;
    abcd[iGrid*756+685] = 4.0E0*I_ESP_Kxy5z_Dxz_aa-2.0E0*3*I_ESP_Hxy3z_Dxz_a-2.0E0*4*I_ESP_Hxy3z_Dxz_a+3*2*I_ESP_Fxyz_Dxz;
    abcd[iGrid*756+686] = 4.0E0*I_ESP_Kx6z_Dxz_aa-2.0E0*4*I_ESP_Hx4z_Dxz_a-2.0E0*5*I_ESP_Hx4z_Dxz_a+4*3*I_ESP_Fx2z_Dxz;
    abcd[iGrid*756+687] = 4.0E0*I_ESP_K5y2z_Dxz_aa-2.0E0*1*I_ESP_H5y_Dxz_a;
    abcd[iGrid*756+688] = 4.0E0*I_ESP_K4y3z_Dxz_aa-2.0E0*1*I_ESP_H4yz_Dxz_a-2.0E0*2*I_ESP_H4yz_Dxz_a;
    abcd[iGrid*756+689] = 4.0E0*I_ESP_K3y4z_Dxz_aa-2.0E0*2*I_ESP_H3y2z_Dxz_a-2.0E0*3*I_ESP_H3y2z_Dxz_a+2*1*I_ESP_F3y_Dxz;
    abcd[iGrid*756+690] = 4.0E0*I_ESP_K2y5z_Dxz_aa-2.0E0*3*I_ESP_H2y3z_Dxz_a-2.0E0*4*I_ESP_H2y3z_Dxz_a+3*2*I_ESP_F2yz_Dxz;
    abcd[iGrid*756+691] = 4.0E0*I_ESP_Ky6z_Dxz_aa-2.0E0*4*I_ESP_Hy4z_Dxz_a-2.0E0*5*I_ESP_Hy4z_Dxz_a+4*3*I_ESP_Fy2z_Dxz;
    abcd[iGrid*756+692] = 4.0E0*I_ESP_K7z_Dxz_aa-2.0E0*5*I_ESP_H5z_Dxz_a-2.0E0*6*I_ESP_H5z_Dxz_a+5*4*I_ESP_F3z_Dxz;
    abcd[iGrid*756+693] = 4.0E0*I_ESP_K5x2z_D2y_aa-2.0E0*1*I_ESP_H5x_D2y_a;
    abcd[iGrid*756+694] = 4.0E0*I_ESP_K4xy2z_D2y_aa-2.0E0*1*I_ESP_H4xy_D2y_a;
    abcd[iGrid*756+695] = 4.0E0*I_ESP_K4x3z_D2y_aa-2.0E0*1*I_ESP_H4xz_D2y_a-2.0E0*2*I_ESP_H4xz_D2y_a;
    abcd[iGrid*756+696] = 4.0E0*I_ESP_K3x2y2z_D2y_aa-2.0E0*1*I_ESP_H3x2y_D2y_a;
    abcd[iGrid*756+697] = 4.0E0*I_ESP_K3xy3z_D2y_aa-2.0E0*1*I_ESP_H3xyz_D2y_a-2.0E0*2*I_ESP_H3xyz_D2y_a;
    abcd[iGrid*756+698] = 4.0E0*I_ESP_K3x4z_D2y_aa-2.0E0*2*I_ESP_H3x2z_D2y_a-2.0E0*3*I_ESP_H3x2z_D2y_a+2*1*I_ESP_F3x_D2y;
    abcd[iGrid*756+699] = 4.0E0*I_ESP_K2x3y2z_D2y_aa-2.0E0*1*I_ESP_H2x3y_D2y_a;
    abcd[iGrid*756+700] = 4.0E0*I_ESP_K2x2y3z_D2y_aa-2.0E0*1*I_ESP_H2x2yz_D2y_a-2.0E0*2*I_ESP_H2x2yz_D2y_a;
    abcd[iGrid*756+701] = 4.0E0*I_ESP_K2xy4z_D2y_aa-2.0E0*2*I_ESP_H2xy2z_D2y_a-2.0E0*3*I_ESP_H2xy2z_D2y_a+2*1*I_ESP_F2xy_D2y;
    abcd[iGrid*756+702] = 4.0E0*I_ESP_K2x5z_D2y_aa-2.0E0*3*I_ESP_H2x3z_D2y_a-2.0E0*4*I_ESP_H2x3z_D2y_a+3*2*I_ESP_F2xz_D2y;
    abcd[iGrid*756+703] = 4.0E0*I_ESP_Kx4y2z_D2y_aa-2.0E0*1*I_ESP_Hx4y_D2y_a;
    abcd[iGrid*756+704] = 4.0E0*I_ESP_Kx3y3z_D2y_aa-2.0E0*1*I_ESP_Hx3yz_D2y_a-2.0E0*2*I_ESP_Hx3yz_D2y_a;
    abcd[iGrid*756+705] = 4.0E0*I_ESP_Kx2y4z_D2y_aa-2.0E0*2*I_ESP_Hx2y2z_D2y_a-2.0E0*3*I_ESP_Hx2y2z_D2y_a+2*1*I_ESP_Fx2y_D2y;
    abcd[iGrid*756+706] = 4.0E0*I_ESP_Kxy5z_D2y_aa-2.0E0*3*I_ESP_Hxy3z_D2y_a-2.0E0*4*I_ESP_Hxy3z_D2y_a+3*2*I_ESP_Fxyz_D2y;
    abcd[iGrid*756+707] = 4.0E0*I_ESP_Kx6z_D2y_aa-2.0E0*4*I_ESP_Hx4z_D2y_a-2.0E0*5*I_ESP_Hx4z_D2y_a+4*3*I_ESP_Fx2z_D2y;
    abcd[iGrid*756+708] = 4.0E0*I_ESP_K5y2z_D2y_aa-2.0E0*1*I_ESP_H5y_D2y_a;
    abcd[iGrid*756+709] = 4.0E0*I_ESP_K4y3z_D2y_aa-2.0E0*1*I_ESP_H4yz_D2y_a-2.0E0*2*I_ESP_H4yz_D2y_a;
    abcd[iGrid*756+710] = 4.0E0*I_ESP_K3y4z_D2y_aa-2.0E0*2*I_ESP_H3y2z_D2y_a-2.0E0*3*I_ESP_H3y2z_D2y_a+2*1*I_ESP_F3y_D2y;
    abcd[iGrid*756+711] = 4.0E0*I_ESP_K2y5z_D2y_aa-2.0E0*3*I_ESP_H2y3z_D2y_a-2.0E0*4*I_ESP_H2y3z_D2y_a+3*2*I_ESP_F2yz_D2y;
    abcd[iGrid*756+712] = 4.0E0*I_ESP_Ky6z_D2y_aa-2.0E0*4*I_ESP_Hy4z_D2y_a-2.0E0*5*I_ESP_Hy4z_D2y_a+4*3*I_ESP_Fy2z_D2y;
    abcd[iGrid*756+713] = 4.0E0*I_ESP_K7z_D2y_aa-2.0E0*5*I_ESP_H5z_D2y_a-2.0E0*6*I_ESP_H5z_D2y_a+5*4*I_ESP_F3z_D2y;
    abcd[iGrid*756+714] = 4.0E0*I_ESP_K5x2z_Dyz_aa-2.0E0*1*I_ESP_H5x_Dyz_a;
    abcd[iGrid*756+715] = 4.0E0*I_ESP_K4xy2z_Dyz_aa-2.0E0*1*I_ESP_H4xy_Dyz_a;
    abcd[iGrid*756+716] = 4.0E0*I_ESP_K4x3z_Dyz_aa-2.0E0*1*I_ESP_H4xz_Dyz_a-2.0E0*2*I_ESP_H4xz_Dyz_a;
    abcd[iGrid*756+717] = 4.0E0*I_ESP_K3x2y2z_Dyz_aa-2.0E0*1*I_ESP_H3x2y_Dyz_a;
    abcd[iGrid*756+718] = 4.0E0*I_ESP_K3xy3z_Dyz_aa-2.0E0*1*I_ESP_H3xyz_Dyz_a-2.0E0*2*I_ESP_H3xyz_Dyz_a;
    abcd[iGrid*756+719] = 4.0E0*I_ESP_K3x4z_Dyz_aa-2.0E0*2*I_ESP_H3x2z_Dyz_a-2.0E0*3*I_ESP_H3x2z_Dyz_a+2*1*I_ESP_F3x_Dyz;
    abcd[iGrid*756+720] = 4.0E0*I_ESP_K2x3y2z_Dyz_aa-2.0E0*1*I_ESP_H2x3y_Dyz_a;
    abcd[iGrid*756+721] = 4.0E0*I_ESP_K2x2y3z_Dyz_aa-2.0E0*1*I_ESP_H2x2yz_Dyz_a-2.0E0*2*I_ESP_H2x2yz_Dyz_a;
    abcd[iGrid*756+722] = 4.0E0*I_ESP_K2xy4z_Dyz_aa-2.0E0*2*I_ESP_H2xy2z_Dyz_a-2.0E0*3*I_ESP_H2xy2z_Dyz_a+2*1*I_ESP_F2xy_Dyz;
    abcd[iGrid*756+723] = 4.0E0*I_ESP_K2x5z_Dyz_aa-2.0E0*3*I_ESP_H2x3z_Dyz_a-2.0E0*4*I_ESP_H2x3z_Dyz_a+3*2*I_ESP_F2xz_Dyz;
    abcd[iGrid*756+724] = 4.0E0*I_ESP_Kx4y2z_Dyz_aa-2.0E0*1*I_ESP_Hx4y_Dyz_a;
    abcd[iGrid*756+725] = 4.0E0*I_ESP_Kx3y3z_Dyz_aa-2.0E0*1*I_ESP_Hx3yz_Dyz_a-2.0E0*2*I_ESP_Hx3yz_Dyz_a;
    abcd[iGrid*756+726] = 4.0E0*I_ESP_Kx2y4z_Dyz_aa-2.0E0*2*I_ESP_Hx2y2z_Dyz_a-2.0E0*3*I_ESP_Hx2y2z_Dyz_a+2*1*I_ESP_Fx2y_Dyz;
    abcd[iGrid*756+727] = 4.0E0*I_ESP_Kxy5z_Dyz_aa-2.0E0*3*I_ESP_Hxy3z_Dyz_a-2.0E0*4*I_ESP_Hxy3z_Dyz_a+3*2*I_ESP_Fxyz_Dyz;
    abcd[iGrid*756+728] = 4.0E0*I_ESP_Kx6z_Dyz_aa-2.0E0*4*I_ESP_Hx4z_Dyz_a-2.0E0*5*I_ESP_Hx4z_Dyz_a+4*3*I_ESP_Fx2z_Dyz;
    abcd[iGrid*756+729] = 4.0E0*I_ESP_K5y2z_Dyz_aa-2.0E0*1*I_ESP_H5y_Dyz_a;
    abcd[iGrid*756+730] = 4.0E0*I_ESP_K4y3z_Dyz_aa-2.0E0*1*I_ESP_H4yz_Dyz_a-2.0E0*2*I_ESP_H4yz_Dyz_a;
    abcd[iGrid*756+731] = 4.0E0*I_ESP_K3y4z_Dyz_aa-2.0E0*2*I_ESP_H3y2z_Dyz_a-2.0E0*3*I_ESP_H3y2z_Dyz_a+2*1*I_ESP_F3y_Dyz;
    abcd[iGrid*756+732] = 4.0E0*I_ESP_K2y5z_Dyz_aa-2.0E0*3*I_ESP_H2y3z_Dyz_a-2.0E0*4*I_ESP_H2y3z_Dyz_a+3*2*I_ESP_F2yz_Dyz;
    abcd[iGrid*756+733] = 4.0E0*I_ESP_Ky6z_Dyz_aa-2.0E0*4*I_ESP_Hy4z_Dyz_a-2.0E0*5*I_ESP_Hy4z_Dyz_a+4*3*I_ESP_Fy2z_Dyz;
    abcd[iGrid*756+734] = 4.0E0*I_ESP_K7z_Dyz_aa-2.0E0*5*I_ESP_H5z_Dyz_a-2.0E0*6*I_ESP_H5z_Dyz_a+5*4*I_ESP_F3z_Dyz;
    abcd[iGrid*756+735] = 4.0E0*I_ESP_K5x2z_D2z_aa-2.0E0*1*I_ESP_H5x_D2z_a;
    abcd[iGrid*756+736] = 4.0E0*I_ESP_K4xy2z_D2z_aa-2.0E0*1*I_ESP_H4xy_D2z_a;
    abcd[iGrid*756+737] = 4.0E0*I_ESP_K4x3z_D2z_aa-2.0E0*1*I_ESP_H4xz_D2z_a-2.0E0*2*I_ESP_H4xz_D2z_a;
    abcd[iGrid*756+738] = 4.0E0*I_ESP_K3x2y2z_D2z_aa-2.0E0*1*I_ESP_H3x2y_D2z_a;
    abcd[iGrid*756+739] = 4.0E0*I_ESP_K3xy3z_D2z_aa-2.0E0*1*I_ESP_H3xyz_D2z_a-2.0E0*2*I_ESP_H3xyz_D2z_a;
    abcd[iGrid*756+740] = 4.0E0*I_ESP_K3x4z_D2z_aa-2.0E0*2*I_ESP_H3x2z_D2z_a-2.0E0*3*I_ESP_H3x2z_D2z_a+2*1*I_ESP_F3x_D2z;
    abcd[iGrid*756+741] = 4.0E0*I_ESP_K2x3y2z_D2z_aa-2.0E0*1*I_ESP_H2x3y_D2z_a;
    abcd[iGrid*756+742] = 4.0E0*I_ESP_K2x2y3z_D2z_aa-2.0E0*1*I_ESP_H2x2yz_D2z_a-2.0E0*2*I_ESP_H2x2yz_D2z_a;
    abcd[iGrid*756+743] = 4.0E0*I_ESP_K2xy4z_D2z_aa-2.0E0*2*I_ESP_H2xy2z_D2z_a-2.0E0*3*I_ESP_H2xy2z_D2z_a+2*1*I_ESP_F2xy_D2z;
    abcd[iGrid*756+744] = 4.0E0*I_ESP_K2x5z_D2z_aa-2.0E0*3*I_ESP_H2x3z_D2z_a-2.0E0*4*I_ESP_H2x3z_D2z_a+3*2*I_ESP_F2xz_D2z;
    abcd[iGrid*756+745] = 4.0E0*I_ESP_Kx4y2z_D2z_aa-2.0E0*1*I_ESP_Hx4y_D2z_a;
    abcd[iGrid*756+746] = 4.0E0*I_ESP_Kx3y3z_D2z_aa-2.0E0*1*I_ESP_Hx3yz_D2z_a-2.0E0*2*I_ESP_Hx3yz_D2z_a;
    abcd[iGrid*756+747] = 4.0E0*I_ESP_Kx2y4z_D2z_aa-2.0E0*2*I_ESP_Hx2y2z_D2z_a-2.0E0*3*I_ESP_Hx2y2z_D2z_a+2*1*I_ESP_Fx2y_D2z;
    abcd[iGrid*756+748] = 4.0E0*I_ESP_Kxy5z_D2z_aa-2.0E0*3*I_ESP_Hxy3z_D2z_a-2.0E0*4*I_ESP_Hxy3z_D2z_a+3*2*I_ESP_Fxyz_D2z;
    abcd[iGrid*756+749] = 4.0E0*I_ESP_Kx6z_D2z_aa-2.0E0*4*I_ESP_Hx4z_D2z_a-2.0E0*5*I_ESP_Hx4z_D2z_a+4*3*I_ESP_Fx2z_D2z;
    abcd[iGrid*756+750] = 4.0E0*I_ESP_K5y2z_D2z_aa-2.0E0*1*I_ESP_H5y_D2z_a;
    abcd[iGrid*756+751] = 4.0E0*I_ESP_K4y3z_D2z_aa-2.0E0*1*I_ESP_H4yz_D2z_a-2.0E0*2*I_ESP_H4yz_D2z_a;
    abcd[iGrid*756+752] = 4.0E0*I_ESP_K3y4z_D2z_aa-2.0E0*2*I_ESP_H3y2z_D2z_a-2.0E0*3*I_ESP_H3y2z_D2z_a+2*1*I_ESP_F3y_D2z;
    abcd[iGrid*756+753] = 4.0E0*I_ESP_K2y5z_D2z_aa-2.0E0*3*I_ESP_H2y3z_D2z_a-2.0E0*4*I_ESP_H2y3z_D2z_a+3*2*I_ESP_F2yz_D2z;
    abcd[iGrid*756+754] = 4.0E0*I_ESP_Ky6z_D2z_aa-2.0E0*4*I_ESP_Hy4z_D2z_a-2.0E0*5*I_ESP_Hy4z_D2z_a+4*3*I_ESP_Fy2z_D2z;
    abcd[iGrid*756+755] = 4.0E0*I_ESP_K7z_D2z_aa-2.0E0*5*I_ESP_H5z_D2z_a-2.0E0*6*I_ESP_H5z_D2z_a+5*4*I_ESP_F3z_D2z;
  }
}
