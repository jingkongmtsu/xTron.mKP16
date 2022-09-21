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
// BRA1 as redundant position, total RHS integrals evaluated as: 8394
// BRA2 as redundant position, total RHS integrals evaluated as: 7204
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

void hgp_os_esp_g_f_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
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
    Double I_ESP_I6x_S_aa = 0.0E0;
    Double I_ESP_I5xy_S_aa = 0.0E0;
    Double I_ESP_I5xz_S_aa = 0.0E0;
    Double I_ESP_I4x2y_S_aa = 0.0E0;
    Double I_ESP_I4xyz_S_aa = 0.0E0;
    Double I_ESP_I4x2z_S_aa = 0.0E0;
    Double I_ESP_I3x3y_S_aa = 0.0E0;
    Double I_ESP_I3x2yz_S_aa = 0.0E0;
    Double I_ESP_I3xy2z_S_aa = 0.0E0;
    Double I_ESP_I3x3z_S_aa = 0.0E0;
    Double I_ESP_I2x4y_S_aa = 0.0E0;
    Double I_ESP_I2x3yz_S_aa = 0.0E0;
    Double I_ESP_I2x2y2z_S_aa = 0.0E0;
    Double I_ESP_I2xy3z_S_aa = 0.0E0;
    Double I_ESP_I2x4z_S_aa = 0.0E0;
    Double I_ESP_Ix5y_S_aa = 0.0E0;
    Double I_ESP_Ix4yz_S_aa = 0.0E0;
    Double I_ESP_Ix3y2z_S_aa = 0.0E0;
    Double I_ESP_Ix2y3z_S_aa = 0.0E0;
    Double I_ESP_Ixy4z_S_aa = 0.0E0;
    Double I_ESP_Ix5z_S_aa = 0.0E0;
    Double I_ESP_I6y_S_aa = 0.0E0;
    Double I_ESP_I5yz_S_aa = 0.0E0;
    Double I_ESP_I4y2z_S_aa = 0.0E0;
    Double I_ESP_I3y3z_S_aa = 0.0E0;
    Double I_ESP_I2y4z_S_aa = 0.0E0;
    Double I_ESP_Iy5z_S_aa = 0.0E0;
    Double I_ESP_I6z_S_aa = 0.0E0;
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
    Double I_ESP_G4x_S_a = 0.0E0;
    Double I_ESP_G3xy_S_a = 0.0E0;
    Double I_ESP_G3xz_S_a = 0.0E0;
    Double I_ESP_G2x2y_S_a = 0.0E0;
    Double I_ESP_G2xyz_S_a = 0.0E0;
    Double I_ESP_G2x2z_S_a = 0.0E0;
    Double I_ESP_Gx3y_S_a = 0.0E0;
    Double I_ESP_Gx2yz_S_a = 0.0E0;
    Double I_ESP_Gxy2z_S_a = 0.0E0;
    Double I_ESP_Gx3z_S_a = 0.0E0;
    Double I_ESP_G4y_S_a = 0.0E0;
    Double I_ESP_G3yz_S_a = 0.0E0;
    Double I_ESP_G2y2z_S_a = 0.0E0;
    Double I_ESP_Gy3z_S_a = 0.0E0;
    Double I_ESP_G4z_S_a = 0.0E0;
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
    Double I_ESP_D2x_S = 0.0E0;
    Double I_ESP_Dxy_S = 0.0E0;
    Double I_ESP_Dxz_S = 0.0E0;
    Double I_ESP_D2y_S = 0.0E0;
    Double I_ESP_Dyz_S = 0.0E0;
    Double I_ESP_D2z_S = 0.0E0;

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
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_D2x_S_vrr = PAX*I_ESP_Px_S_vrr-PRX*I_ESP_Px_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_Dxy_S_vrr = PAY*I_ESP_Px_S_vrr-PRY*I_ESP_Px_S_M1_vrr;
        Double I_ESP_Dxz_S_vrr = PAZ*I_ESP_Px_S_vrr-PRZ*I_ESP_Px_S_M1_vrr;
        Double I_ESP_D2y_S_vrr = PAY*I_ESP_Py_S_vrr-PRY*I_ESP_Py_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_Dyz_S_vrr = PAZ*I_ESP_Py_S_vrr-PRZ*I_ESP_Py_S_M1_vrr;
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
         * shell quartet name: SQ_ESP_I_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_I_S_aa_coefs = alpha*alpha;
        I_ESP_I6x_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S_aa += SQ_ESP_I_S_aa_coefs*I_ESP_I6z_S_vrr;

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
         * shell quartet name: SQ_ESP_G_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_G_S_a_coefs = alpha;
        I_ESP_G4x_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G4x_S_vrr;
        I_ESP_G3xy_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G3xy_S_vrr;
        I_ESP_G3xz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G3xz_S_vrr;
        I_ESP_G2x2y_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2x2y_S_vrr;
        I_ESP_G2xyz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2xyz_S_vrr;
        I_ESP_G2x2z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2x2z_S_vrr;
        I_ESP_Gx3y_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gx3y_S_vrr;
        I_ESP_Gx2yz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gx2yz_S_vrr;
        I_ESP_Gxy2z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gxy2z_S_vrr;
        I_ESP_Gx3z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gx3z_S_vrr;
        I_ESP_G4y_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G4y_S_vrr;
        I_ESP_G3yz_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G3yz_S_vrr;
        I_ESP_G2y2z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G2y2z_S_vrr;
        I_ESP_Gy3z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_Gy3z_S_vrr;
        I_ESP_G4z_S_a += SQ_ESP_G_S_a_coefs*I_ESP_G4z_S_vrr;

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

        /************************************************************
         * shell quartet name: SQ_ESP_D_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_D2x_S += I_ESP_D2x_S_vrr;
        I_ESP_Dxy_S += I_ESP_Dxy_S_vrr;
        I_ESP_Dxz_S += I_ESP_Dxz_S_vrr;
        I_ESP_D2y_S += I_ESP_D2y_S_vrr;
        I_ESP_Dyz_S += I_ESP_Dyz_S_vrr;
        I_ESP_D2z_S += I_ESP_D2z_S_vrr;
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
     * shell quartet name: SQ_ESP_D_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_F_S
     * RHS shell quartet name: SQ_ESP_D_S
     ************************************************************/
    Double I_ESP_D2x_Px = I_ESP_F3x_S+ABX*I_ESP_D2x_S;
    Double I_ESP_Dxy_Px = I_ESP_F2xy_S+ABX*I_ESP_Dxy_S;
    Double I_ESP_Dxz_Px = I_ESP_F2xz_S+ABX*I_ESP_Dxz_S;
    Double I_ESP_D2y_Px = I_ESP_Fx2y_S+ABX*I_ESP_D2y_S;
    Double I_ESP_Dyz_Px = I_ESP_Fxyz_S+ABX*I_ESP_Dyz_S;
    Double I_ESP_D2z_Px = I_ESP_Fx2z_S+ABX*I_ESP_D2z_S;
    Double I_ESP_D2x_Py = I_ESP_F2xy_S+ABY*I_ESP_D2x_S;
    Double I_ESP_Dxy_Py = I_ESP_Fx2y_S+ABY*I_ESP_Dxy_S;
    Double I_ESP_Dxz_Py = I_ESP_Fxyz_S+ABY*I_ESP_Dxz_S;
    Double I_ESP_D2y_Py = I_ESP_F3y_S+ABY*I_ESP_D2y_S;
    Double I_ESP_Dyz_Py = I_ESP_F2yz_S+ABY*I_ESP_Dyz_S;
    Double I_ESP_D2z_Py = I_ESP_Fy2z_S+ABY*I_ESP_D2z_S;
    Double I_ESP_D2x_Pz = I_ESP_F2xz_S+ABZ*I_ESP_D2x_S;
    Double I_ESP_Dxy_Pz = I_ESP_Fxyz_S+ABZ*I_ESP_Dxy_S;
    Double I_ESP_Dxz_Pz = I_ESP_Fx2z_S+ABZ*I_ESP_Dxz_S;
    Double I_ESP_D2y_Pz = I_ESP_F2yz_S+ABZ*I_ESP_D2y_S;
    Double I_ESP_Dyz_Pz = I_ESP_Fy2z_S+ABZ*I_ESP_Dyz_S;
    Double I_ESP_D2z_Pz = I_ESP_F3z_S+ABZ*I_ESP_D2z_S;

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
     * shell quartet name: SQ_ESP_D_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 12 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_F_P
     * RHS shell quartet name: SQ_ESP_D_P
     ************************************************************/
    Double I_ESP_D2x_D2x = I_ESP_F3x_Px+ABX*I_ESP_D2x_Px;
    Double I_ESP_Dxy_D2x = I_ESP_F2xy_Px+ABX*I_ESP_Dxy_Px;
    Double I_ESP_Dxz_D2x = I_ESP_F2xz_Px+ABX*I_ESP_Dxz_Px;
    Double I_ESP_D2y_D2x = I_ESP_Fx2y_Px+ABX*I_ESP_D2y_Px;
    Double I_ESP_Dyz_D2x = I_ESP_Fxyz_Px+ABX*I_ESP_Dyz_Px;
    Double I_ESP_D2z_D2x = I_ESP_Fx2z_Px+ABX*I_ESP_D2z_Px;
    Double I_ESP_D2x_Dxy = I_ESP_F2xy_Px+ABY*I_ESP_D2x_Px;
    Double I_ESP_Dxy_Dxy = I_ESP_Fx2y_Px+ABY*I_ESP_Dxy_Px;
    Double I_ESP_Dxz_Dxy = I_ESP_Fxyz_Px+ABY*I_ESP_Dxz_Px;
    Double I_ESP_D2y_Dxy = I_ESP_F3y_Px+ABY*I_ESP_D2y_Px;
    Double I_ESP_Dyz_Dxy = I_ESP_F2yz_Px+ABY*I_ESP_Dyz_Px;
    Double I_ESP_D2z_Dxy = I_ESP_Fy2z_Px+ABY*I_ESP_D2z_Px;
    Double I_ESP_D2x_D2y = I_ESP_F2xy_Py+ABY*I_ESP_D2x_Py;
    Double I_ESP_Dxy_D2y = I_ESP_Fx2y_Py+ABY*I_ESP_Dxy_Py;
    Double I_ESP_Dxz_D2y = I_ESP_Fxyz_Py+ABY*I_ESP_Dxz_Py;
    Double I_ESP_D2y_D2y = I_ESP_F3y_Py+ABY*I_ESP_D2y_Py;
    Double I_ESP_Dyz_D2y = I_ESP_F2yz_Py+ABY*I_ESP_Dyz_Py;
    Double I_ESP_D2z_D2y = I_ESP_Fy2z_Py+ABY*I_ESP_D2z_Py;
    Double I_ESP_D2x_D2z = I_ESP_F2xz_Pz+ABZ*I_ESP_D2x_Pz;
    Double I_ESP_Dxy_D2z = I_ESP_Fxyz_Pz+ABZ*I_ESP_Dxy_Pz;
    Double I_ESP_Dxz_D2z = I_ESP_Fx2z_Pz+ABZ*I_ESP_Dxz_Pz;
    Double I_ESP_D2y_D2z = I_ESP_F2yz_Pz+ABZ*I_ESP_D2y_Pz;
    Double I_ESP_Dyz_D2z = I_ESP_Fy2z_Pz+ABZ*I_ESP_Dyz_Pz;
    Double I_ESP_D2z_D2z = I_ESP_F3z_Pz+ABZ*I_ESP_D2z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 12 integrals are omitted 
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
    Double I_ESP_G3yz_Px = I_ESP_Hx3yz_S+ABX*I_ESP_G3yz_S;
    Double I_ESP_G2y2z_Px = I_ESP_Hx2y2z_S+ABX*I_ESP_G2y2z_S;
    Double I_ESP_Gy3z_Px = I_ESP_Hxy3z_S+ABX*I_ESP_Gy3z_S;
    Double I_ESP_G3xy_Py = I_ESP_H3x2y_S+ABY*I_ESP_G3xy_S;
    Double I_ESP_G2x2y_Py = I_ESP_H2x3y_S+ABY*I_ESP_G2x2y_S;
    Double I_ESP_G2xyz_Py = I_ESP_H2x2yz_S+ABY*I_ESP_G2xyz_S;
    Double I_ESP_Gx3y_Py = I_ESP_Hx4y_S+ABY*I_ESP_Gx3y_S;
    Double I_ESP_Gx2yz_Py = I_ESP_Hx3yz_S+ABY*I_ESP_Gx2yz_S;
    Double I_ESP_Gxy2z_Py = I_ESP_Hx2y2z_S+ABY*I_ESP_Gxy2z_S;
    Double I_ESP_G4y_Py = I_ESP_H5y_S+ABY*I_ESP_G4y_S;
    Double I_ESP_G3yz_Py = I_ESP_H4yz_S+ABY*I_ESP_G3yz_S;
    Double I_ESP_G2y2z_Py = I_ESP_H3y2z_S+ABY*I_ESP_G2y2z_S;
    Double I_ESP_Gy3z_Py = I_ESP_H2y3z_S+ABY*I_ESP_Gy3z_S;
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
     * totally 24 integrals are omitted 
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
    Double I_ESP_F2xz_Dxy = I_ESP_G2xyz_Px+ABY*I_ESP_F2xz_Px;
    Double I_ESP_Fxyz_Dxy = I_ESP_Gx2yz_Px+ABY*I_ESP_Fxyz_Px;
    Double I_ESP_Fx2z_Dxy = I_ESP_Gxy2z_Px+ABY*I_ESP_Fx2z_Px;
    Double I_ESP_F2yz_Dxy = I_ESP_G3yz_Px+ABY*I_ESP_F2yz_Px;
    Double I_ESP_Fy2z_Dxy = I_ESP_G2y2z_Px+ABY*I_ESP_Fy2z_Px;
    Double I_ESP_F3z_Dxy = I_ESP_Gy3z_Px+ABY*I_ESP_F3z_Px;
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
     * shell quartet name: SQ_ESP_D_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_F_D
     * RHS shell quartet name: SQ_ESP_D_D
     ************************************************************/
    Double I_ESP_D2x_F3x = I_ESP_F3x_D2x+ABX*I_ESP_D2x_D2x;
    Double I_ESP_Dxy_F3x = I_ESP_F2xy_D2x+ABX*I_ESP_Dxy_D2x;
    Double I_ESP_Dxz_F3x = I_ESP_F2xz_D2x+ABX*I_ESP_Dxz_D2x;
    Double I_ESP_D2y_F3x = I_ESP_Fx2y_D2x+ABX*I_ESP_D2y_D2x;
    Double I_ESP_Dyz_F3x = I_ESP_Fxyz_D2x+ABX*I_ESP_Dyz_D2x;
    Double I_ESP_D2z_F3x = I_ESP_Fx2z_D2x+ABX*I_ESP_D2z_D2x;
    Double I_ESP_D2x_F2xy = I_ESP_F2xy_D2x+ABY*I_ESP_D2x_D2x;
    Double I_ESP_Dxy_F2xy = I_ESP_Fx2y_D2x+ABY*I_ESP_Dxy_D2x;
    Double I_ESP_Dxz_F2xy = I_ESP_Fxyz_D2x+ABY*I_ESP_Dxz_D2x;
    Double I_ESP_D2y_F2xy = I_ESP_F3y_D2x+ABY*I_ESP_D2y_D2x;
    Double I_ESP_Dyz_F2xy = I_ESP_F2yz_D2x+ABY*I_ESP_Dyz_D2x;
    Double I_ESP_D2z_F2xy = I_ESP_Fy2z_D2x+ABY*I_ESP_D2z_D2x;
    Double I_ESP_D2x_F2xz = I_ESP_F2xz_D2x+ABZ*I_ESP_D2x_D2x;
    Double I_ESP_Dxy_F2xz = I_ESP_Fxyz_D2x+ABZ*I_ESP_Dxy_D2x;
    Double I_ESP_Dxz_F2xz = I_ESP_Fx2z_D2x+ABZ*I_ESP_Dxz_D2x;
    Double I_ESP_D2y_F2xz = I_ESP_F2yz_D2x+ABZ*I_ESP_D2y_D2x;
    Double I_ESP_Dyz_F2xz = I_ESP_Fy2z_D2x+ABZ*I_ESP_Dyz_D2x;
    Double I_ESP_D2z_F2xz = I_ESP_F3z_D2x+ABZ*I_ESP_D2z_D2x;
    Double I_ESP_D2x_Fx2y = I_ESP_F3x_D2y+ABX*I_ESP_D2x_D2y;
    Double I_ESP_Dxy_Fx2y = I_ESP_F2xy_D2y+ABX*I_ESP_Dxy_D2y;
    Double I_ESP_Dxz_Fx2y = I_ESP_F2xz_D2y+ABX*I_ESP_Dxz_D2y;
    Double I_ESP_D2y_Fx2y = I_ESP_Fx2y_D2y+ABX*I_ESP_D2y_D2y;
    Double I_ESP_Dyz_Fx2y = I_ESP_Fxyz_D2y+ABX*I_ESP_Dyz_D2y;
    Double I_ESP_D2z_Fx2y = I_ESP_Fx2z_D2y+ABX*I_ESP_D2z_D2y;
    Double I_ESP_D2x_Fxyz = I_ESP_F2xz_Dxy+ABZ*I_ESP_D2x_Dxy;
    Double I_ESP_Dxy_Fxyz = I_ESP_Fxyz_Dxy+ABZ*I_ESP_Dxy_Dxy;
    Double I_ESP_Dxz_Fxyz = I_ESP_Fx2z_Dxy+ABZ*I_ESP_Dxz_Dxy;
    Double I_ESP_D2y_Fxyz = I_ESP_F2yz_Dxy+ABZ*I_ESP_D2y_Dxy;
    Double I_ESP_Dyz_Fxyz = I_ESP_Fy2z_Dxy+ABZ*I_ESP_Dyz_Dxy;
    Double I_ESP_D2z_Fxyz = I_ESP_F3z_Dxy+ABZ*I_ESP_D2z_Dxy;
    Double I_ESP_D2x_Fx2z = I_ESP_F3x_D2z+ABX*I_ESP_D2x_D2z;
    Double I_ESP_Dxy_Fx2z = I_ESP_F2xy_D2z+ABX*I_ESP_Dxy_D2z;
    Double I_ESP_Dxz_Fx2z = I_ESP_F2xz_D2z+ABX*I_ESP_Dxz_D2z;
    Double I_ESP_D2y_Fx2z = I_ESP_Fx2y_D2z+ABX*I_ESP_D2y_D2z;
    Double I_ESP_Dyz_Fx2z = I_ESP_Fxyz_D2z+ABX*I_ESP_Dyz_D2z;
    Double I_ESP_D2z_Fx2z = I_ESP_Fx2z_D2z+ABX*I_ESP_D2z_D2z;
    Double I_ESP_D2x_F3y = I_ESP_F2xy_D2y+ABY*I_ESP_D2x_D2y;
    Double I_ESP_Dxy_F3y = I_ESP_Fx2y_D2y+ABY*I_ESP_Dxy_D2y;
    Double I_ESP_Dxz_F3y = I_ESP_Fxyz_D2y+ABY*I_ESP_Dxz_D2y;
    Double I_ESP_D2y_F3y = I_ESP_F3y_D2y+ABY*I_ESP_D2y_D2y;
    Double I_ESP_Dyz_F3y = I_ESP_F2yz_D2y+ABY*I_ESP_Dyz_D2y;
    Double I_ESP_D2z_F3y = I_ESP_Fy2z_D2y+ABY*I_ESP_D2z_D2y;
    Double I_ESP_D2x_F2yz = I_ESP_F2xz_D2y+ABZ*I_ESP_D2x_D2y;
    Double I_ESP_Dxy_F2yz = I_ESP_Fxyz_D2y+ABZ*I_ESP_Dxy_D2y;
    Double I_ESP_Dxz_F2yz = I_ESP_Fx2z_D2y+ABZ*I_ESP_Dxz_D2y;
    Double I_ESP_D2y_F2yz = I_ESP_F2yz_D2y+ABZ*I_ESP_D2y_D2y;
    Double I_ESP_Dyz_F2yz = I_ESP_Fy2z_D2y+ABZ*I_ESP_Dyz_D2y;
    Double I_ESP_D2z_F2yz = I_ESP_F3z_D2y+ABZ*I_ESP_D2z_D2y;
    Double I_ESP_D2x_Fy2z = I_ESP_F2xy_D2z+ABY*I_ESP_D2x_D2z;
    Double I_ESP_Dxy_Fy2z = I_ESP_Fx2y_D2z+ABY*I_ESP_Dxy_D2z;
    Double I_ESP_Dxz_Fy2z = I_ESP_Fxyz_D2z+ABY*I_ESP_Dxz_D2z;
    Double I_ESP_D2y_Fy2z = I_ESP_F3y_D2z+ABY*I_ESP_D2y_D2z;
    Double I_ESP_Dyz_Fy2z = I_ESP_F2yz_D2z+ABY*I_ESP_Dyz_D2z;
    Double I_ESP_D2z_Fy2z = I_ESP_Fy2z_D2z+ABY*I_ESP_D2z_D2z;
    Double I_ESP_D2x_F3z = I_ESP_F2xz_D2z+ABZ*I_ESP_D2x_D2z;
    Double I_ESP_Dxy_F3z = I_ESP_Fxyz_D2z+ABZ*I_ESP_Dxy_D2z;
    Double I_ESP_Dxz_F3z = I_ESP_Fx2z_D2z+ABZ*I_ESP_Dxz_D2z;
    Double I_ESP_D2y_F3z = I_ESP_F2yz_D2z+ABZ*I_ESP_D2y_D2z;
    Double I_ESP_Dyz_F3z = I_ESP_Fy2z_D2z+ABZ*I_ESP_Dyz_D2z;
    Double I_ESP_D2z_F3z = I_ESP_F3z_D2z+ABZ*I_ESP_D2z_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_G_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S_a
     * RHS shell quartet name: SQ_ESP_G_S_a
     ************************************************************/
    Double I_ESP_G4x_Px_a = I_ESP_H5x_S_a+ABX*I_ESP_G4x_S_a;
    Double I_ESP_G3xy_Px_a = I_ESP_H4xy_S_a+ABX*I_ESP_G3xy_S_a;
    Double I_ESP_G3xz_Px_a = I_ESP_H4xz_S_a+ABX*I_ESP_G3xz_S_a;
    Double I_ESP_G2x2y_Px_a = I_ESP_H3x2y_S_a+ABX*I_ESP_G2x2y_S_a;
    Double I_ESP_G2xyz_Px_a = I_ESP_H3xyz_S_a+ABX*I_ESP_G2xyz_S_a;
    Double I_ESP_G2x2z_Px_a = I_ESP_H3x2z_S_a+ABX*I_ESP_G2x2z_S_a;
    Double I_ESP_Gx3y_Px_a = I_ESP_H2x3y_S_a+ABX*I_ESP_Gx3y_S_a;
    Double I_ESP_Gx2yz_Px_a = I_ESP_H2x2yz_S_a+ABX*I_ESP_Gx2yz_S_a;
    Double I_ESP_Gxy2z_Px_a = I_ESP_H2xy2z_S_a+ABX*I_ESP_Gxy2z_S_a;
    Double I_ESP_Gx3z_Px_a = I_ESP_H2x3z_S_a+ABX*I_ESP_Gx3z_S_a;
    Double I_ESP_G4y_Px_a = I_ESP_Hx4y_S_a+ABX*I_ESP_G4y_S_a;
    Double I_ESP_G3yz_Px_a = I_ESP_Hx3yz_S_a+ABX*I_ESP_G3yz_S_a;
    Double I_ESP_G2y2z_Px_a = I_ESP_Hx2y2z_S_a+ABX*I_ESP_G2y2z_S_a;
    Double I_ESP_Gy3z_Px_a = I_ESP_Hxy3z_S_a+ABX*I_ESP_Gy3z_S_a;
    Double I_ESP_G4z_Px_a = I_ESP_Hx4z_S_a+ABX*I_ESP_G4z_S_a;
    Double I_ESP_G4x_Py_a = I_ESP_H4xy_S_a+ABY*I_ESP_G4x_S_a;
    Double I_ESP_G3xy_Py_a = I_ESP_H3x2y_S_a+ABY*I_ESP_G3xy_S_a;
    Double I_ESP_G3xz_Py_a = I_ESP_H3xyz_S_a+ABY*I_ESP_G3xz_S_a;
    Double I_ESP_G2x2y_Py_a = I_ESP_H2x3y_S_a+ABY*I_ESP_G2x2y_S_a;
    Double I_ESP_G2xyz_Py_a = I_ESP_H2x2yz_S_a+ABY*I_ESP_G2xyz_S_a;
    Double I_ESP_G2x2z_Py_a = I_ESP_H2xy2z_S_a+ABY*I_ESP_G2x2z_S_a;
    Double I_ESP_Gx3y_Py_a = I_ESP_Hx4y_S_a+ABY*I_ESP_Gx3y_S_a;
    Double I_ESP_Gx2yz_Py_a = I_ESP_Hx3yz_S_a+ABY*I_ESP_Gx2yz_S_a;
    Double I_ESP_Gxy2z_Py_a = I_ESP_Hx2y2z_S_a+ABY*I_ESP_Gxy2z_S_a;
    Double I_ESP_Gx3z_Py_a = I_ESP_Hxy3z_S_a+ABY*I_ESP_Gx3z_S_a;
    Double I_ESP_G4y_Py_a = I_ESP_H5y_S_a+ABY*I_ESP_G4y_S_a;
    Double I_ESP_G3yz_Py_a = I_ESP_H4yz_S_a+ABY*I_ESP_G3yz_S_a;
    Double I_ESP_G2y2z_Py_a = I_ESP_H3y2z_S_a+ABY*I_ESP_G2y2z_S_a;
    Double I_ESP_Gy3z_Py_a = I_ESP_H2y3z_S_a+ABY*I_ESP_Gy3z_S_a;
    Double I_ESP_G4z_Py_a = I_ESP_Hy4z_S_a+ABY*I_ESP_G4z_S_a;
    Double I_ESP_G4x_Pz_a = I_ESP_H4xz_S_a+ABZ*I_ESP_G4x_S_a;
    Double I_ESP_G3xy_Pz_a = I_ESP_H3xyz_S_a+ABZ*I_ESP_G3xy_S_a;
    Double I_ESP_G3xz_Pz_a = I_ESP_H3x2z_S_a+ABZ*I_ESP_G3xz_S_a;
    Double I_ESP_G2x2y_Pz_a = I_ESP_H2x2yz_S_a+ABZ*I_ESP_G2x2y_S_a;
    Double I_ESP_G2xyz_Pz_a = I_ESP_H2xy2z_S_a+ABZ*I_ESP_G2xyz_S_a;
    Double I_ESP_G2x2z_Pz_a = I_ESP_H2x3z_S_a+ABZ*I_ESP_G2x2z_S_a;
    Double I_ESP_Gx3y_Pz_a = I_ESP_Hx3yz_S_a+ABZ*I_ESP_Gx3y_S_a;
    Double I_ESP_Gx2yz_Pz_a = I_ESP_Hx2y2z_S_a+ABZ*I_ESP_Gx2yz_S_a;
    Double I_ESP_Gxy2z_Pz_a = I_ESP_Hxy3z_S_a+ABZ*I_ESP_Gxy2z_S_a;
    Double I_ESP_Gx3z_Pz_a = I_ESP_Hx4z_S_a+ABZ*I_ESP_Gx3z_S_a;
    Double I_ESP_G4y_Pz_a = I_ESP_H4yz_S_a+ABZ*I_ESP_G4y_S_a;
    Double I_ESP_G3yz_Pz_a = I_ESP_H3y2z_S_a+ABZ*I_ESP_G3yz_S_a;
    Double I_ESP_G2y2z_Pz_a = I_ESP_H2y3z_S_a+ABZ*I_ESP_G2y2z_S_a;
    Double I_ESP_Gy3z_Pz_a = I_ESP_Hy4z_S_a+ABZ*I_ESP_Gy3z_S_a;
    Double I_ESP_G4z_Pz_a = I_ESP_H5z_S_a+ABZ*I_ESP_G4z_S_a;

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
     * shell quartet name: SQ_ESP_G_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P_a
     * RHS shell quartet name: SQ_ESP_G_P_a
     ************************************************************/
    Double I_ESP_G4x_D2x_a = I_ESP_H5x_Px_a+ABX*I_ESP_G4x_Px_a;
    Double I_ESP_G3xy_D2x_a = I_ESP_H4xy_Px_a+ABX*I_ESP_G3xy_Px_a;
    Double I_ESP_G3xz_D2x_a = I_ESP_H4xz_Px_a+ABX*I_ESP_G3xz_Px_a;
    Double I_ESP_G2x2y_D2x_a = I_ESP_H3x2y_Px_a+ABX*I_ESP_G2x2y_Px_a;
    Double I_ESP_G2xyz_D2x_a = I_ESP_H3xyz_Px_a+ABX*I_ESP_G2xyz_Px_a;
    Double I_ESP_G2x2z_D2x_a = I_ESP_H3x2z_Px_a+ABX*I_ESP_G2x2z_Px_a;
    Double I_ESP_Gx3y_D2x_a = I_ESP_H2x3y_Px_a+ABX*I_ESP_Gx3y_Px_a;
    Double I_ESP_Gx2yz_D2x_a = I_ESP_H2x2yz_Px_a+ABX*I_ESP_Gx2yz_Px_a;
    Double I_ESP_Gxy2z_D2x_a = I_ESP_H2xy2z_Px_a+ABX*I_ESP_Gxy2z_Px_a;
    Double I_ESP_Gx3z_D2x_a = I_ESP_H2x3z_Px_a+ABX*I_ESP_Gx3z_Px_a;
    Double I_ESP_G4y_D2x_a = I_ESP_Hx4y_Px_a+ABX*I_ESP_G4y_Px_a;
    Double I_ESP_G3yz_D2x_a = I_ESP_Hx3yz_Px_a+ABX*I_ESP_G3yz_Px_a;
    Double I_ESP_G2y2z_D2x_a = I_ESP_Hx2y2z_Px_a+ABX*I_ESP_G2y2z_Px_a;
    Double I_ESP_Gy3z_D2x_a = I_ESP_Hxy3z_Px_a+ABX*I_ESP_Gy3z_Px_a;
    Double I_ESP_G4z_D2x_a = I_ESP_Hx4z_Px_a+ABX*I_ESP_G4z_Px_a;
    Double I_ESP_G4x_Dxy_a = I_ESP_H4xy_Px_a+ABY*I_ESP_G4x_Px_a;
    Double I_ESP_G3xy_Dxy_a = I_ESP_H3x2y_Px_a+ABY*I_ESP_G3xy_Px_a;
    Double I_ESP_G3xz_Dxy_a = I_ESP_H3xyz_Px_a+ABY*I_ESP_G3xz_Px_a;
    Double I_ESP_G2x2y_Dxy_a = I_ESP_H2x3y_Px_a+ABY*I_ESP_G2x2y_Px_a;
    Double I_ESP_G2xyz_Dxy_a = I_ESP_H2x2yz_Px_a+ABY*I_ESP_G2xyz_Px_a;
    Double I_ESP_G2x2z_Dxy_a = I_ESP_H2xy2z_Px_a+ABY*I_ESP_G2x2z_Px_a;
    Double I_ESP_Gx3y_Dxy_a = I_ESP_Hx4y_Px_a+ABY*I_ESP_Gx3y_Px_a;
    Double I_ESP_Gx2yz_Dxy_a = I_ESP_Hx3yz_Px_a+ABY*I_ESP_Gx2yz_Px_a;
    Double I_ESP_Gxy2z_Dxy_a = I_ESP_Hx2y2z_Px_a+ABY*I_ESP_Gxy2z_Px_a;
    Double I_ESP_Gx3z_Dxy_a = I_ESP_Hxy3z_Px_a+ABY*I_ESP_Gx3z_Px_a;
    Double I_ESP_G4y_Dxy_a = I_ESP_H5y_Px_a+ABY*I_ESP_G4y_Px_a;
    Double I_ESP_G3yz_Dxy_a = I_ESP_H4yz_Px_a+ABY*I_ESP_G3yz_Px_a;
    Double I_ESP_G2y2z_Dxy_a = I_ESP_H3y2z_Px_a+ABY*I_ESP_G2y2z_Px_a;
    Double I_ESP_Gy3z_Dxy_a = I_ESP_H2y3z_Px_a+ABY*I_ESP_Gy3z_Px_a;
    Double I_ESP_G4z_Dxy_a = I_ESP_Hy4z_Px_a+ABY*I_ESP_G4z_Px_a;
    Double I_ESP_G4x_D2y_a = I_ESP_H4xy_Py_a+ABY*I_ESP_G4x_Py_a;
    Double I_ESP_G3xy_D2y_a = I_ESP_H3x2y_Py_a+ABY*I_ESP_G3xy_Py_a;
    Double I_ESP_G3xz_D2y_a = I_ESP_H3xyz_Py_a+ABY*I_ESP_G3xz_Py_a;
    Double I_ESP_G2x2y_D2y_a = I_ESP_H2x3y_Py_a+ABY*I_ESP_G2x2y_Py_a;
    Double I_ESP_G2xyz_D2y_a = I_ESP_H2x2yz_Py_a+ABY*I_ESP_G2xyz_Py_a;
    Double I_ESP_G2x2z_D2y_a = I_ESP_H2xy2z_Py_a+ABY*I_ESP_G2x2z_Py_a;
    Double I_ESP_Gx3y_D2y_a = I_ESP_Hx4y_Py_a+ABY*I_ESP_Gx3y_Py_a;
    Double I_ESP_Gx2yz_D2y_a = I_ESP_Hx3yz_Py_a+ABY*I_ESP_Gx2yz_Py_a;
    Double I_ESP_Gxy2z_D2y_a = I_ESP_Hx2y2z_Py_a+ABY*I_ESP_Gxy2z_Py_a;
    Double I_ESP_Gx3z_D2y_a = I_ESP_Hxy3z_Py_a+ABY*I_ESP_Gx3z_Py_a;
    Double I_ESP_G4y_D2y_a = I_ESP_H5y_Py_a+ABY*I_ESP_G4y_Py_a;
    Double I_ESP_G3yz_D2y_a = I_ESP_H4yz_Py_a+ABY*I_ESP_G3yz_Py_a;
    Double I_ESP_G2y2z_D2y_a = I_ESP_H3y2z_Py_a+ABY*I_ESP_G2y2z_Py_a;
    Double I_ESP_Gy3z_D2y_a = I_ESP_H2y3z_Py_a+ABY*I_ESP_Gy3z_Py_a;
    Double I_ESP_G4z_D2y_a = I_ESP_Hy4z_Py_a+ABY*I_ESP_G4z_Py_a;
    Double I_ESP_G4x_D2z_a = I_ESP_H4xz_Pz_a+ABZ*I_ESP_G4x_Pz_a;
    Double I_ESP_G3xy_D2z_a = I_ESP_H3xyz_Pz_a+ABZ*I_ESP_G3xy_Pz_a;
    Double I_ESP_G3xz_D2z_a = I_ESP_H3x2z_Pz_a+ABZ*I_ESP_G3xz_Pz_a;
    Double I_ESP_G2x2y_D2z_a = I_ESP_H2x2yz_Pz_a+ABZ*I_ESP_G2x2y_Pz_a;
    Double I_ESP_G2xyz_D2z_a = I_ESP_H2xy2z_Pz_a+ABZ*I_ESP_G2xyz_Pz_a;
    Double I_ESP_G2x2z_D2z_a = I_ESP_H2x3z_Pz_a+ABZ*I_ESP_G2x2z_Pz_a;
    Double I_ESP_Gx3y_D2z_a = I_ESP_Hx3yz_Pz_a+ABZ*I_ESP_Gx3y_Pz_a;
    Double I_ESP_Gx2yz_D2z_a = I_ESP_Hx2y2z_Pz_a+ABZ*I_ESP_Gx2yz_Pz_a;
    Double I_ESP_Gxy2z_D2z_a = I_ESP_Hxy3z_Pz_a+ABZ*I_ESP_Gxy2z_Pz_a;
    Double I_ESP_Gx3z_D2z_a = I_ESP_Hx4z_Pz_a+ABZ*I_ESP_Gx3z_Pz_a;
    Double I_ESP_G4y_D2z_a = I_ESP_H4yz_Pz_a+ABZ*I_ESP_G4y_Pz_a;
    Double I_ESP_G3yz_D2z_a = I_ESP_H3y2z_Pz_a+ABZ*I_ESP_G3yz_Pz_a;
    Double I_ESP_G2y2z_D2z_a = I_ESP_H2y3z_Pz_a+ABZ*I_ESP_G2y2z_Pz_a;
    Double I_ESP_Gy3z_D2z_a = I_ESP_Hy4z_Pz_a+ABZ*I_ESP_Gy3z_Pz_a;
    Double I_ESP_G4z_D2z_a = I_ESP_H5z_Pz_a+ABZ*I_ESP_G4z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 16 integrals are omitted 
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
    Double I_ESP_I5yz_Px_a = I_ESP_Kx5yz_S_a+ABX*I_ESP_I5yz_S_a;
    Double I_ESP_I4y2z_Px_a = I_ESP_Kx4y2z_S_a+ABX*I_ESP_I4y2z_S_a;
    Double I_ESP_I3y3z_Px_a = I_ESP_Kx3y3z_S_a+ABX*I_ESP_I3y3z_S_a;
    Double I_ESP_I2y4z_Px_a = I_ESP_Kx2y4z_S_a+ABX*I_ESP_I2y4z_S_a;
    Double I_ESP_Iy5z_Px_a = I_ESP_Kxy5z_S_a+ABX*I_ESP_Iy5z_S_a;
    Double I_ESP_I5xy_Py_a = I_ESP_K5x2y_S_a+ABY*I_ESP_I5xy_S_a;
    Double I_ESP_I4x2y_Py_a = I_ESP_K4x3y_S_a+ABY*I_ESP_I4x2y_S_a;
    Double I_ESP_I4xyz_Py_a = I_ESP_K4x2yz_S_a+ABY*I_ESP_I4xyz_S_a;
    Double I_ESP_I3x3y_Py_a = I_ESP_K3x4y_S_a+ABY*I_ESP_I3x3y_S_a;
    Double I_ESP_I3x2yz_Py_a = I_ESP_K3x3yz_S_a+ABY*I_ESP_I3x2yz_S_a;
    Double I_ESP_I3xy2z_Py_a = I_ESP_K3x2y2z_S_a+ABY*I_ESP_I3xy2z_S_a;
    Double I_ESP_I2x4y_Py_a = I_ESP_K2x5y_S_a+ABY*I_ESP_I2x4y_S_a;
    Double I_ESP_I2x3yz_Py_a = I_ESP_K2x4yz_S_a+ABY*I_ESP_I2x3yz_S_a;
    Double I_ESP_I2x2y2z_Py_a = I_ESP_K2x3y2z_S_a+ABY*I_ESP_I2x2y2z_S_a;
    Double I_ESP_I2xy3z_Py_a = I_ESP_K2x2y3z_S_a+ABY*I_ESP_I2xy3z_S_a;
    Double I_ESP_Ix5y_Py_a = I_ESP_Kx6y_S_a+ABY*I_ESP_Ix5y_S_a;
    Double I_ESP_Ix4yz_Py_a = I_ESP_Kx5yz_S_a+ABY*I_ESP_Ix4yz_S_a;
    Double I_ESP_Ix3y2z_Py_a = I_ESP_Kx4y2z_S_a+ABY*I_ESP_Ix3y2z_S_a;
    Double I_ESP_Ix2y3z_Py_a = I_ESP_Kx3y3z_S_a+ABY*I_ESP_Ix2y3z_S_a;
    Double I_ESP_Ixy4z_Py_a = I_ESP_Kx2y4z_S_a+ABY*I_ESP_Ixy4z_S_a;
    Double I_ESP_I6y_Py_a = I_ESP_K7y_S_a+ABY*I_ESP_I6y_S_a;
    Double I_ESP_I5yz_Py_a = I_ESP_K6yz_S_a+ABY*I_ESP_I5yz_S_a;
    Double I_ESP_I4y2z_Py_a = I_ESP_K5y2z_S_a+ABY*I_ESP_I4y2z_S_a;
    Double I_ESP_I3y3z_Py_a = I_ESP_K4y3z_S_a+ABY*I_ESP_I3y3z_S_a;
    Double I_ESP_I2y4z_Py_a = I_ESP_K3y4z_S_a+ABY*I_ESP_I2y4z_S_a;
    Double I_ESP_Iy5z_Py_a = I_ESP_K2y5z_S_a+ABY*I_ESP_Iy5z_S_a;
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
     * totally 48 integrals are omitted 
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
    Double I_ESP_H4xz_Dxy_a = I_ESP_I4xyz_Px_a+ABY*I_ESP_H4xz_Px_a;
    Double I_ESP_H3xyz_Dxy_a = I_ESP_I3x2yz_Px_a+ABY*I_ESP_H3xyz_Px_a;
    Double I_ESP_H3x2z_Dxy_a = I_ESP_I3xy2z_Px_a+ABY*I_ESP_H3x2z_Px_a;
    Double I_ESP_H2x2yz_Dxy_a = I_ESP_I2x3yz_Px_a+ABY*I_ESP_H2x2yz_Px_a;
    Double I_ESP_H2xy2z_Dxy_a = I_ESP_I2x2y2z_Px_a+ABY*I_ESP_H2xy2z_Px_a;
    Double I_ESP_H2x3z_Dxy_a = I_ESP_I2xy3z_Px_a+ABY*I_ESP_H2x3z_Px_a;
    Double I_ESP_Hx3yz_Dxy_a = I_ESP_Ix4yz_Px_a+ABY*I_ESP_Hx3yz_Px_a;
    Double I_ESP_Hx2y2z_Dxy_a = I_ESP_Ix3y2z_Px_a+ABY*I_ESP_Hx2y2z_Px_a;
    Double I_ESP_Hxy3z_Dxy_a = I_ESP_Ix2y3z_Px_a+ABY*I_ESP_Hxy3z_Px_a;
    Double I_ESP_Hx4z_Dxy_a = I_ESP_Ixy4z_Px_a+ABY*I_ESP_Hx4z_Px_a;
    Double I_ESP_H4yz_Dxy_a = I_ESP_I5yz_Px_a+ABY*I_ESP_H4yz_Px_a;
    Double I_ESP_H3y2z_Dxy_a = I_ESP_I4y2z_Px_a+ABY*I_ESP_H3y2z_Px_a;
    Double I_ESP_H2y3z_Dxy_a = I_ESP_I3y3z_Px_a+ABY*I_ESP_H2y3z_Px_a;
    Double I_ESP_Hy4z_Dxy_a = I_ESP_I2y4z_Px_a+ABY*I_ESP_Hy4z_Px_a;
    Double I_ESP_H5z_Dxy_a = I_ESP_Iy5z_Px_a+ABY*I_ESP_H5z_Px_a;
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
     * shell quartet name: SQ_ESP_G_F_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_D_a
     * RHS shell quartet name: SQ_ESP_G_D_a
     ************************************************************/
    Double I_ESP_G4x_F3x_a = I_ESP_H5x_D2x_a+ABX*I_ESP_G4x_D2x_a;
    Double I_ESP_G3xy_F3x_a = I_ESP_H4xy_D2x_a+ABX*I_ESP_G3xy_D2x_a;
    Double I_ESP_G3xz_F3x_a = I_ESP_H4xz_D2x_a+ABX*I_ESP_G3xz_D2x_a;
    Double I_ESP_G2x2y_F3x_a = I_ESP_H3x2y_D2x_a+ABX*I_ESP_G2x2y_D2x_a;
    Double I_ESP_G2xyz_F3x_a = I_ESP_H3xyz_D2x_a+ABX*I_ESP_G2xyz_D2x_a;
    Double I_ESP_G2x2z_F3x_a = I_ESP_H3x2z_D2x_a+ABX*I_ESP_G2x2z_D2x_a;
    Double I_ESP_Gx3y_F3x_a = I_ESP_H2x3y_D2x_a+ABX*I_ESP_Gx3y_D2x_a;
    Double I_ESP_Gx2yz_F3x_a = I_ESP_H2x2yz_D2x_a+ABX*I_ESP_Gx2yz_D2x_a;
    Double I_ESP_Gxy2z_F3x_a = I_ESP_H2xy2z_D2x_a+ABX*I_ESP_Gxy2z_D2x_a;
    Double I_ESP_Gx3z_F3x_a = I_ESP_H2x3z_D2x_a+ABX*I_ESP_Gx3z_D2x_a;
    Double I_ESP_G4y_F3x_a = I_ESP_Hx4y_D2x_a+ABX*I_ESP_G4y_D2x_a;
    Double I_ESP_G3yz_F3x_a = I_ESP_Hx3yz_D2x_a+ABX*I_ESP_G3yz_D2x_a;
    Double I_ESP_G2y2z_F3x_a = I_ESP_Hx2y2z_D2x_a+ABX*I_ESP_G2y2z_D2x_a;
    Double I_ESP_Gy3z_F3x_a = I_ESP_Hxy3z_D2x_a+ABX*I_ESP_Gy3z_D2x_a;
    Double I_ESP_G4z_F3x_a = I_ESP_Hx4z_D2x_a+ABX*I_ESP_G4z_D2x_a;
    Double I_ESP_G4x_F2xy_a = I_ESP_H4xy_D2x_a+ABY*I_ESP_G4x_D2x_a;
    Double I_ESP_G3xy_F2xy_a = I_ESP_H3x2y_D2x_a+ABY*I_ESP_G3xy_D2x_a;
    Double I_ESP_G3xz_F2xy_a = I_ESP_H3xyz_D2x_a+ABY*I_ESP_G3xz_D2x_a;
    Double I_ESP_G2x2y_F2xy_a = I_ESP_H2x3y_D2x_a+ABY*I_ESP_G2x2y_D2x_a;
    Double I_ESP_G2xyz_F2xy_a = I_ESP_H2x2yz_D2x_a+ABY*I_ESP_G2xyz_D2x_a;
    Double I_ESP_G2x2z_F2xy_a = I_ESP_H2xy2z_D2x_a+ABY*I_ESP_G2x2z_D2x_a;
    Double I_ESP_Gx3y_F2xy_a = I_ESP_Hx4y_D2x_a+ABY*I_ESP_Gx3y_D2x_a;
    Double I_ESP_Gx2yz_F2xy_a = I_ESP_Hx3yz_D2x_a+ABY*I_ESP_Gx2yz_D2x_a;
    Double I_ESP_Gxy2z_F2xy_a = I_ESP_Hx2y2z_D2x_a+ABY*I_ESP_Gxy2z_D2x_a;
    Double I_ESP_Gx3z_F2xy_a = I_ESP_Hxy3z_D2x_a+ABY*I_ESP_Gx3z_D2x_a;
    Double I_ESP_G4y_F2xy_a = I_ESP_H5y_D2x_a+ABY*I_ESP_G4y_D2x_a;
    Double I_ESP_G3yz_F2xy_a = I_ESP_H4yz_D2x_a+ABY*I_ESP_G3yz_D2x_a;
    Double I_ESP_G2y2z_F2xy_a = I_ESP_H3y2z_D2x_a+ABY*I_ESP_G2y2z_D2x_a;
    Double I_ESP_Gy3z_F2xy_a = I_ESP_H2y3z_D2x_a+ABY*I_ESP_Gy3z_D2x_a;
    Double I_ESP_G4z_F2xy_a = I_ESP_Hy4z_D2x_a+ABY*I_ESP_G4z_D2x_a;
    Double I_ESP_G4x_F2xz_a = I_ESP_H4xz_D2x_a+ABZ*I_ESP_G4x_D2x_a;
    Double I_ESP_G3xy_F2xz_a = I_ESP_H3xyz_D2x_a+ABZ*I_ESP_G3xy_D2x_a;
    Double I_ESP_G3xz_F2xz_a = I_ESP_H3x2z_D2x_a+ABZ*I_ESP_G3xz_D2x_a;
    Double I_ESP_G2x2y_F2xz_a = I_ESP_H2x2yz_D2x_a+ABZ*I_ESP_G2x2y_D2x_a;
    Double I_ESP_G2xyz_F2xz_a = I_ESP_H2xy2z_D2x_a+ABZ*I_ESP_G2xyz_D2x_a;
    Double I_ESP_G2x2z_F2xz_a = I_ESP_H2x3z_D2x_a+ABZ*I_ESP_G2x2z_D2x_a;
    Double I_ESP_Gx3y_F2xz_a = I_ESP_Hx3yz_D2x_a+ABZ*I_ESP_Gx3y_D2x_a;
    Double I_ESP_Gx2yz_F2xz_a = I_ESP_Hx2y2z_D2x_a+ABZ*I_ESP_Gx2yz_D2x_a;
    Double I_ESP_Gxy2z_F2xz_a = I_ESP_Hxy3z_D2x_a+ABZ*I_ESP_Gxy2z_D2x_a;
    Double I_ESP_Gx3z_F2xz_a = I_ESP_Hx4z_D2x_a+ABZ*I_ESP_Gx3z_D2x_a;
    Double I_ESP_G4y_F2xz_a = I_ESP_H4yz_D2x_a+ABZ*I_ESP_G4y_D2x_a;
    Double I_ESP_G3yz_F2xz_a = I_ESP_H3y2z_D2x_a+ABZ*I_ESP_G3yz_D2x_a;
    Double I_ESP_G2y2z_F2xz_a = I_ESP_H2y3z_D2x_a+ABZ*I_ESP_G2y2z_D2x_a;
    Double I_ESP_Gy3z_F2xz_a = I_ESP_Hy4z_D2x_a+ABZ*I_ESP_Gy3z_D2x_a;
    Double I_ESP_G4z_F2xz_a = I_ESP_H5z_D2x_a+ABZ*I_ESP_G4z_D2x_a;
    Double I_ESP_G4x_Fx2y_a = I_ESP_H5x_D2y_a+ABX*I_ESP_G4x_D2y_a;
    Double I_ESP_G3xy_Fx2y_a = I_ESP_H4xy_D2y_a+ABX*I_ESP_G3xy_D2y_a;
    Double I_ESP_G3xz_Fx2y_a = I_ESP_H4xz_D2y_a+ABX*I_ESP_G3xz_D2y_a;
    Double I_ESP_G2x2y_Fx2y_a = I_ESP_H3x2y_D2y_a+ABX*I_ESP_G2x2y_D2y_a;
    Double I_ESP_G2xyz_Fx2y_a = I_ESP_H3xyz_D2y_a+ABX*I_ESP_G2xyz_D2y_a;
    Double I_ESP_G2x2z_Fx2y_a = I_ESP_H3x2z_D2y_a+ABX*I_ESP_G2x2z_D2y_a;
    Double I_ESP_Gx3y_Fx2y_a = I_ESP_H2x3y_D2y_a+ABX*I_ESP_Gx3y_D2y_a;
    Double I_ESP_Gx2yz_Fx2y_a = I_ESP_H2x2yz_D2y_a+ABX*I_ESP_Gx2yz_D2y_a;
    Double I_ESP_Gxy2z_Fx2y_a = I_ESP_H2xy2z_D2y_a+ABX*I_ESP_Gxy2z_D2y_a;
    Double I_ESP_Gx3z_Fx2y_a = I_ESP_H2x3z_D2y_a+ABX*I_ESP_Gx3z_D2y_a;
    Double I_ESP_G4y_Fx2y_a = I_ESP_Hx4y_D2y_a+ABX*I_ESP_G4y_D2y_a;
    Double I_ESP_G3yz_Fx2y_a = I_ESP_Hx3yz_D2y_a+ABX*I_ESP_G3yz_D2y_a;
    Double I_ESP_G2y2z_Fx2y_a = I_ESP_Hx2y2z_D2y_a+ABX*I_ESP_G2y2z_D2y_a;
    Double I_ESP_Gy3z_Fx2y_a = I_ESP_Hxy3z_D2y_a+ABX*I_ESP_Gy3z_D2y_a;
    Double I_ESP_G4z_Fx2y_a = I_ESP_Hx4z_D2y_a+ABX*I_ESP_G4z_D2y_a;
    Double I_ESP_G4x_Fxyz_a = I_ESP_H4xz_Dxy_a+ABZ*I_ESP_G4x_Dxy_a;
    Double I_ESP_G3xy_Fxyz_a = I_ESP_H3xyz_Dxy_a+ABZ*I_ESP_G3xy_Dxy_a;
    Double I_ESP_G3xz_Fxyz_a = I_ESP_H3x2z_Dxy_a+ABZ*I_ESP_G3xz_Dxy_a;
    Double I_ESP_G2x2y_Fxyz_a = I_ESP_H2x2yz_Dxy_a+ABZ*I_ESP_G2x2y_Dxy_a;
    Double I_ESP_G2xyz_Fxyz_a = I_ESP_H2xy2z_Dxy_a+ABZ*I_ESP_G2xyz_Dxy_a;
    Double I_ESP_G2x2z_Fxyz_a = I_ESP_H2x3z_Dxy_a+ABZ*I_ESP_G2x2z_Dxy_a;
    Double I_ESP_Gx3y_Fxyz_a = I_ESP_Hx3yz_Dxy_a+ABZ*I_ESP_Gx3y_Dxy_a;
    Double I_ESP_Gx2yz_Fxyz_a = I_ESP_Hx2y2z_Dxy_a+ABZ*I_ESP_Gx2yz_Dxy_a;
    Double I_ESP_Gxy2z_Fxyz_a = I_ESP_Hxy3z_Dxy_a+ABZ*I_ESP_Gxy2z_Dxy_a;
    Double I_ESP_Gx3z_Fxyz_a = I_ESP_Hx4z_Dxy_a+ABZ*I_ESP_Gx3z_Dxy_a;
    Double I_ESP_G4y_Fxyz_a = I_ESP_H4yz_Dxy_a+ABZ*I_ESP_G4y_Dxy_a;
    Double I_ESP_G3yz_Fxyz_a = I_ESP_H3y2z_Dxy_a+ABZ*I_ESP_G3yz_Dxy_a;
    Double I_ESP_G2y2z_Fxyz_a = I_ESP_H2y3z_Dxy_a+ABZ*I_ESP_G2y2z_Dxy_a;
    Double I_ESP_Gy3z_Fxyz_a = I_ESP_Hy4z_Dxy_a+ABZ*I_ESP_Gy3z_Dxy_a;
    Double I_ESP_G4z_Fxyz_a = I_ESP_H5z_Dxy_a+ABZ*I_ESP_G4z_Dxy_a;
    Double I_ESP_G4x_Fx2z_a = I_ESP_H5x_D2z_a+ABX*I_ESP_G4x_D2z_a;
    Double I_ESP_G3xy_Fx2z_a = I_ESP_H4xy_D2z_a+ABX*I_ESP_G3xy_D2z_a;
    Double I_ESP_G3xz_Fx2z_a = I_ESP_H4xz_D2z_a+ABX*I_ESP_G3xz_D2z_a;
    Double I_ESP_G2x2y_Fx2z_a = I_ESP_H3x2y_D2z_a+ABX*I_ESP_G2x2y_D2z_a;
    Double I_ESP_G2xyz_Fx2z_a = I_ESP_H3xyz_D2z_a+ABX*I_ESP_G2xyz_D2z_a;
    Double I_ESP_G2x2z_Fx2z_a = I_ESP_H3x2z_D2z_a+ABX*I_ESP_G2x2z_D2z_a;
    Double I_ESP_Gx3y_Fx2z_a = I_ESP_H2x3y_D2z_a+ABX*I_ESP_Gx3y_D2z_a;
    Double I_ESP_Gx2yz_Fx2z_a = I_ESP_H2x2yz_D2z_a+ABX*I_ESP_Gx2yz_D2z_a;
    Double I_ESP_Gxy2z_Fx2z_a = I_ESP_H2xy2z_D2z_a+ABX*I_ESP_Gxy2z_D2z_a;
    Double I_ESP_Gx3z_Fx2z_a = I_ESP_H2x3z_D2z_a+ABX*I_ESP_Gx3z_D2z_a;
    Double I_ESP_G4y_Fx2z_a = I_ESP_Hx4y_D2z_a+ABX*I_ESP_G4y_D2z_a;
    Double I_ESP_G3yz_Fx2z_a = I_ESP_Hx3yz_D2z_a+ABX*I_ESP_G3yz_D2z_a;
    Double I_ESP_G2y2z_Fx2z_a = I_ESP_Hx2y2z_D2z_a+ABX*I_ESP_G2y2z_D2z_a;
    Double I_ESP_Gy3z_Fx2z_a = I_ESP_Hxy3z_D2z_a+ABX*I_ESP_Gy3z_D2z_a;
    Double I_ESP_G4z_Fx2z_a = I_ESP_Hx4z_D2z_a+ABX*I_ESP_G4z_D2z_a;
    Double I_ESP_G4x_F3y_a = I_ESP_H4xy_D2y_a+ABY*I_ESP_G4x_D2y_a;
    Double I_ESP_G3xy_F3y_a = I_ESP_H3x2y_D2y_a+ABY*I_ESP_G3xy_D2y_a;
    Double I_ESP_G3xz_F3y_a = I_ESP_H3xyz_D2y_a+ABY*I_ESP_G3xz_D2y_a;
    Double I_ESP_G2x2y_F3y_a = I_ESP_H2x3y_D2y_a+ABY*I_ESP_G2x2y_D2y_a;
    Double I_ESP_G2xyz_F3y_a = I_ESP_H2x2yz_D2y_a+ABY*I_ESP_G2xyz_D2y_a;
    Double I_ESP_G2x2z_F3y_a = I_ESP_H2xy2z_D2y_a+ABY*I_ESP_G2x2z_D2y_a;
    Double I_ESP_Gx3y_F3y_a = I_ESP_Hx4y_D2y_a+ABY*I_ESP_Gx3y_D2y_a;
    Double I_ESP_Gx2yz_F3y_a = I_ESP_Hx3yz_D2y_a+ABY*I_ESP_Gx2yz_D2y_a;
    Double I_ESP_Gxy2z_F3y_a = I_ESP_Hx2y2z_D2y_a+ABY*I_ESP_Gxy2z_D2y_a;
    Double I_ESP_Gx3z_F3y_a = I_ESP_Hxy3z_D2y_a+ABY*I_ESP_Gx3z_D2y_a;
    Double I_ESP_G4y_F3y_a = I_ESP_H5y_D2y_a+ABY*I_ESP_G4y_D2y_a;
    Double I_ESP_G3yz_F3y_a = I_ESP_H4yz_D2y_a+ABY*I_ESP_G3yz_D2y_a;
    Double I_ESP_G2y2z_F3y_a = I_ESP_H3y2z_D2y_a+ABY*I_ESP_G2y2z_D2y_a;
    Double I_ESP_Gy3z_F3y_a = I_ESP_H2y3z_D2y_a+ABY*I_ESP_Gy3z_D2y_a;
    Double I_ESP_G4z_F3y_a = I_ESP_Hy4z_D2y_a+ABY*I_ESP_G4z_D2y_a;
    Double I_ESP_G4x_F2yz_a = I_ESP_H4xz_D2y_a+ABZ*I_ESP_G4x_D2y_a;
    Double I_ESP_G3xy_F2yz_a = I_ESP_H3xyz_D2y_a+ABZ*I_ESP_G3xy_D2y_a;
    Double I_ESP_G3xz_F2yz_a = I_ESP_H3x2z_D2y_a+ABZ*I_ESP_G3xz_D2y_a;
    Double I_ESP_G2x2y_F2yz_a = I_ESP_H2x2yz_D2y_a+ABZ*I_ESP_G2x2y_D2y_a;
    Double I_ESP_G2xyz_F2yz_a = I_ESP_H2xy2z_D2y_a+ABZ*I_ESP_G2xyz_D2y_a;
    Double I_ESP_G2x2z_F2yz_a = I_ESP_H2x3z_D2y_a+ABZ*I_ESP_G2x2z_D2y_a;
    Double I_ESP_Gx3y_F2yz_a = I_ESP_Hx3yz_D2y_a+ABZ*I_ESP_Gx3y_D2y_a;
    Double I_ESP_Gx2yz_F2yz_a = I_ESP_Hx2y2z_D2y_a+ABZ*I_ESP_Gx2yz_D2y_a;
    Double I_ESP_Gxy2z_F2yz_a = I_ESP_Hxy3z_D2y_a+ABZ*I_ESP_Gxy2z_D2y_a;
    Double I_ESP_Gx3z_F2yz_a = I_ESP_Hx4z_D2y_a+ABZ*I_ESP_Gx3z_D2y_a;
    Double I_ESP_G4y_F2yz_a = I_ESP_H4yz_D2y_a+ABZ*I_ESP_G4y_D2y_a;
    Double I_ESP_G3yz_F2yz_a = I_ESP_H3y2z_D2y_a+ABZ*I_ESP_G3yz_D2y_a;
    Double I_ESP_G2y2z_F2yz_a = I_ESP_H2y3z_D2y_a+ABZ*I_ESP_G2y2z_D2y_a;
    Double I_ESP_Gy3z_F2yz_a = I_ESP_Hy4z_D2y_a+ABZ*I_ESP_Gy3z_D2y_a;
    Double I_ESP_G4z_F2yz_a = I_ESP_H5z_D2y_a+ABZ*I_ESP_G4z_D2y_a;
    Double I_ESP_G4x_Fy2z_a = I_ESP_H4xy_D2z_a+ABY*I_ESP_G4x_D2z_a;
    Double I_ESP_G3xy_Fy2z_a = I_ESP_H3x2y_D2z_a+ABY*I_ESP_G3xy_D2z_a;
    Double I_ESP_G3xz_Fy2z_a = I_ESP_H3xyz_D2z_a+ABY*I_ESP_G3xz_D2z_a;
    Double I_ESP_G2x2y_Fy2z_a = I_ESP_H2x3y_D2z_a+ABY*I_ESP_G2x2y_D2z_a;
    Double I_ESP_G2xyz_Fy2z_a = I_ESP_H2x2yz_D2z_a+ABY*I_ESP_G2xyz_D2z_a;
    Double I_ESP_G2x2z_Fy2z_a = I_ESP_H2xy2z_D2z_a+ABY*I_ESP_G2x2z_D2z_a;
    Double I_ESP_Gx3y_Fy2z_a = I_ESP_Hx4y_D2z_a+ABY*I_ESP_Gx3y_D2z_a;
    Double I_ESP_Gx2yz_Fy2z_a = I_ESP_Hx3yz_D2z_a+ABY*I_ESP_Gx2yz_D2z_a;
    Double I_ESP_Gxy2z_Fy2z_a = I_ESP_Hx2y2z_D2z_a+ABY*I_ESP_Gxy2z_D2z_a;
    Double I_ESP_Gx3z_Fy2z_a = I_ESP_Hxy3z_D2z_a+ABY*I_ESP_Gx3z_D2z_a;
    Double I_ESP_G4y_Fy2z_a = I_ESP_H5y_D2z_a+ABY*I_ESP_G4y_D2z_a;
    Double I_ESP_G3yz_Fy2z_a = I_ESP_H4yz_D2z_a+ABY*I_ESP_G3yz_D2z_a;
    Double I_ESP_G2y2z_Fy2z_a = I_ESP_H3y2z_D2z_a+ABY*I_ESP_G2y2z_D2z_a;
    Double I_ESP_Gy3z_Fy2z_a = I_ESP_H2y3z_D2z_a+ABY*I_ESP_Gy3z_D2z_a;
    Double I_ESP_G4z_Fy2z_a = I_ESP_Hy4z_D2z_a+ABY*I_ESP_G4z_D2z_a;
    Double I_ESP_G4x_F3z_a = I_ESP_H4xz_D2z_a+ABZ*I_ESP_G4x_D2z_a;
    Double I_ESP_G3xy_F3z_a = I_ESP_H3xyz_D2z_a+ABZ*I_ESP_G3xy_D2z_a;
    Double I_ESP_G3xz_F3z_a = I_ESP_H3x2z_D2z_a+ABZ*I_ESP_G3xz_D2z_a;
    Double I_ESP_G2x2y_F3z_a = I_ESP_H2x2yz_D2z_a+ABZ*I_ESP_G2x2y_D2z_a;
    Double I_ESP_G2xyz_F3z_a = I_ESP_H2xy2z_D2z_a+ABZ*I_ESP_G2xyz_D2z_a;
    Double I_ESP_G2x2z_F3z_a = I_ESP_H2x3z_D2z_a+ABZ*I_ESP_G2x2z_D2z_a;
    Double I_ESP_Gx3y_F3z_a = I_ESP_Hx3yz_D2z_a+ABZ*I_ESP_Gx3y_D2z_a;
    Double I_ESP_Gx2yz_F3z_a = I_ESP_Hx2y2z_D2z_a+ABZ*I_ESP_Gx2yz_D2z_a;
    Double I_ESP_Gxy2z_F3z_a = I_ESP_Hxy3z_D2z_a+ABZ*I_ESP_Gxy2z_D2z_a;
    Double I_ESP_Gx3z_F3z_a = I_ESP_Hx4z_D2z_a+ABZ*I_ESP_Gx3z_D2z_a;
    Double I_ESP_G4y_F3z_a = I_ESP_H4yz_D2z_a+ABZ*I_ESP_G4y_D2z_a;
    Double I_ESP_G3yz_F3z_a = I_ESP_H3y2z_D2z_a+ABZ*I_ESP_G3yz_D2z_a;
    Double I_ESP_G2y2z_F3z_a = I_ESP_H2y3z_D2z_a+ABZ*I_ESP_G2y2z_D2z_a;
    Double I_ESP_Gy3z_F3z_a = I_ESP_Hy4z_D2z_a+ABZ*I_ESP_Gy3z_D2z_a;
    Double I_ESP_G4z_F3z_a = I_ESP_H5z_D2z_a+ABZ*I_ESP_G4z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_S_aa
     * RHS shell quartet name: SQ_ESP_I_S_aa
     ************************************************************/
    Double I_ESP_I6x_Px_aa = I_ESP_K7x_S_aa+ABX*I_ESP_I6x_S_aa;
    Double I_ESP_I5xy_Px_aa = I_ESP_K6xy_S_aa+ABX*I_ESP_I5xy_S_aa;
    Double I_ESP_I5xz_Px_aa = I_ESP_K6xz_S_aa+ABX*I_ESP_I5xz_S_aa;
    Double I_ESP_I4x2y_Px_aa = I_ESP_K5x2y_S_aa+ABX*I_ESP_I4x2y_S_aa;
    Double I_ESP_I4xyz_Px_aa = I_ESP_K5xyz_S_aa+ABX*I_ESP_I4xyz_S_aa;
    Double I_ESP_I4x2z_Px_aa = I_ESP_K5x2z_S_aa+ABX*I_ESP_I4x2z_S_aa;
    Double I_ESP_I3x3y_Px_aa = I_ESP_K4x3y_S_aa+ABX*I_ESP_I3x3y_S_aa;
    Double I_ESP_I3x2yz_Px_aa = I_ESP_K4x2yz_S_aa+ABX*I_ESP_I3x2yz_S_aa;
    Double I_ESP_I3xy2z_Px_aa = I_ESP_K4xy2z_S_aa+ABX*I_ESP_I3xy2z_S_aa;
    Double I_ESP_I3x3z_Px_aa = I_ESP_K4x3z_S_aa+ABX*I_ESP_I3x3z_S_aa;
    Double I_ESP_I2x4y_Px_aa = I_ESP_K3x4y_S_aa+ABX*I_ESP_I2x4y_S_aa;
    Double I_ESP_I2x3yz_Px_aa = I_ESP_K3x3yz_S_aa+ABX*I_ESP_I2x3yz_S_aa;
    Double I_ESP_I2x2y2z_Px_aa = I_ESP_K3x2y2z_S_aa+ABX*I_ESP_I2x2y2z_S_aa;
    Double I_ESP_I2xy3z_Px_aa = I_ESP_K3xy3z_S_aa+ABX*I_ESP_I2xy3z_S_aa;
    Double I_ESP_I2x4z_Px_aa = I_ESP_K3x4z_S_aa+ABX*I_ESP_I2x4z_S_aa;
    Double I_ESP_Ix5y_Px_aa = I_ESP_K2x5y_S_aa+ABX*I_ESP_Ix5y_S_aa;
    Double I_ESP_Ix4yz_Px_aa = I_ESP_K2x4yz_S_aa+ABX*I_ESP_Ix4yz_S_aa;
    Double I_ESP_Ix3y2z_Px_aa = I_ESP_K2x3y2z_S_aa+ABX*I_ESP_Ix3y2z_S_aa;
    Double I_ESP_Ix2y3z_Px_aa = I_ESP_K2x2y3z_S_aa+ABX*I_ESP_Ix2y3z_S_aa;
    Double I_ESP_Ixy4z_Px_aa = I_ESP_K2xy4z_S_aa+ABX*I_ESP_Ixy4z_S_aa;
    Double I_ESP_Ix5z_Px_aa = I_ESP_K2x5z_S_aa+ABX*I_ESP_Ix5z_S_aa;
    Double I_ESP_I6y_Px_aa = I_ESP_Kx6y_S_aa+ABX*I_ESP_I6y_S_aa;
    Double I_ESP_I5yz_Px_aa = I_ESP_Kx5yz_S_aa+ABX*I_ESP_I5yz_S_aa;
    Double I_ESP_I4y2z_Px_aa = I_ESP_Kx4y2z_S_aa+ABX*I_ESP_I4y2z_S_aa;
    Double I_ESP_I3y3z_Px_aa = I_ESP_Kx3y3z_S_aa+ABX*I_ESP_I3y3z_S_aa;
    Double I_ESP_I2y4z_Px_aa = I_ESP_Kx2y4z_S_aa+ABX*I_ESP_I2y4z_S_aa;
    Double I_ESP_Iy5z_Px_aa = I_ESP_Kxy5z_S_aa+ABX*I_ESP_Iy5z_S_aa;
    Double I_ESP_I6z_Px_aa = I_ESP_Kx6z_S_aa+ABX*I_ESP_I6z_S_aa;
    Double I_ESP_I6x_Py_aa = I_ESP_K6xy_S_aa+ABY*I_ESP_I6x_S_aa;
    Double I_ESP_I5xy_Py_aa = I_ESP_K5x2y_S_aa+ABY*I_ESP_I5xy_S_aa;
    Double I_ESP_I5xz_Py_aa = I_ESP_K5xyz_S_aa+ABY*I_ESP_I5xz_S_aa;
    Double I_ESP_I4x2y_Py_aa = I_ESP_K4x3y_S_aa+ABY*I_ESP_I4x2y_S_aa;
    Double I_ESP_I4xyz_Py_aa = I_ESP_K4x2yz_S_aa+ABY*I_ESP_I4xyz_S_aa;
    Double I_ESP_I4x2z_Py_aa = I_ESP_K4xy2z_S_aa+ABY*I_ESP_I4x2z_S_aa;
    Double I_ESP_I3x3y_Py_aa = I_ESP_K3x4y_S_aa+ABY*I_ESP_I3x3y_S_aa;
    Double I_ESP_I3x2yz_Py_aa = I_ESP_K3x3yz_S_aa+ABY*I_ESP_I3x2yz_S_aa;
    Double I_ESP_I3xy2z_Py_aa = I_ESP_K3x2y2z_S_aa+ABY*I_ESP_I3xy2z_S_aa;
    Double I_ESP_I3x3z_Py_aa = I_ESP_K3xy3z_S_aa+ABY*I_ESP_I3x3z_S_aa;
    Double I_ESP_I2x4y_Py_aa = I_ESP_K2x5y_S_aa+ABY*I_ESP_I2x4y_S_aa;
    Double I_ESP_I2x3yz_Py_aa = I_ESP_K2x4yz_S_aa+ABY*I_ESP_I2x3yz_S_aa;
    Double I_ESP_I2x2y2z_Py_aa = I_ESP_K2x3y2z_S_aa+ABY*I_ESP_I2x2y2z_S_aa;
    Double I_ESP_I2xy3z_Py_aa = I_ESP_K2x2y3z_S_aa+ABY*I_ESP_I2xy3z_S_aa;
    Double I_ESP_I2x4z_Py_aa = I_ESP_K2xy4z_S_aa+ABY*I_ESP_I2x4z_S_aa;
    Double I_ESP_Ix5y_Py_aa = I_ESP_Kx6y_S_aa+ABY*I_ESP_Ix5y_S_aa;
    Double I_ESP_Ix4yz_Py_aa = I_ESP_Kx5yz_S_aa+ABY*I_ESP_Ix4yz_S_aa;
    Double I_ESP_Ix3y2z_Py_aa = I_ESP_Kx4y2z_S_aa+ABY*I_ESP_Ix3y2z_S_aa;
    Double I_ESP_Ix2y3z_Py_aa = I_ESP_Kx3y3z_S_aa+ABY*I_ESP_Ix2y3z_S_aa;
    Double I_ESP_Ixy4z_Py_aa = I_ESP_Kx2y4z_S_aa+ABY*I_ESP_Ixy4z_S_aa;
    Double I_ESP_Ix5z_Py_aa = I_ESP_Kxy5z_S_aa+ABY*I_ESP_Ix5z_S_aa;
    Double I_ESP_I6y_Py_aa = I_ESP_K7y_S_aa+ABY*I_ESP_I6y_S_aa;
    Double I_ESP_I5yz_Py_aa = I_ESP_K6yz_S_aa+ABY*I_ESP_I5yz_S_aa;
    Double I_ESP_I4y2z_Py_aa = I_ESP_K5y2z_S_aa+ABY*I_ESP_I4y2z_S_aa;
    Double I_ESP_I3y3z_Py_aa = I_ESP_K4y3z_S_aa+ABY*I_ESP_I3y3z_S_aa;
    Double I_ESP_I2y4z_Py_aa = I_ESP_K3y4z_S_aa+ABY*I_ESP_I2y4z_S_aa;
    Double I_ESP_Iy5z_Py_aa = I_ESP_K2y5z_S_aa+ABY*I_ESP_Iy5z_S_aa;
    Double I_ESP_I6z_Py_aa = I_ESP_Ky6z_S_aa+ABY*I_ESP_I6z_S_aa;
    Double I_ESP_I6x_Pz_aa = I_ESP_K6xz_S_aa+ABZ*I_ESP_I6x_S_aa;
    Double I_ESP_I5xy_Pz_aa = I_ESP_K5xyz_S_aa+ABZ*I_ESP_I5xy_S_aa;
    Double I_ESP_I5xz_Pz_aa = I_ESP_K5x2z_S_aa+ABZ*I_ESP_I5xz_S_aa;
    Double I_ESP_I4x2y_Pz_aa = I_ESP_K4x2yz_S_aa+ABZ*I_ESP_I4x2y_S_aa;
    Double I_ESP_I4xyz_Pz_aa = I_ESP_K4xy2z_S_aa+ABZ*I_ESP_I4xyz_S_aa;
    Double I_ESP_I4x2z_Pz_aa = I_ESP_K4x3z_S_aa+ABZ*I_ESP_I4x2z_S_aa;
    Double I_ESP_I3x3y_Pz_aa = I_ESP_K3x3yz_S_aa+ABZ*I_ESP_I3x3y_S_aa;
    Double I_ESP_I3x2yz_Pz_aa = I_ESP_K3x2y2z_S_aa+ABZ*I_ESP_I3x2yz_S_aa;
    Double I_ESP_I3xy2z_Pz_aa = I_ESP_K3xy3z_S_aa+ABZ*I_ESP_I3xy2z_S_aa;
    Double I_ESP_I3x3z_Pz_aa = I_ESP_K3x4z_S_aa+ABZ*I_ESP_I3x3z_S_aa;
    Double I_ESP_I2x4y_Pz_aa = I_ESP_K2x4yz_S_aa+ABZ*I_ESP_I2x4y_S_aa;
    Double I_ESP_I2x3yz_Pz_aa = I_ESP_K2x3y2z_S_aa+ABZ*I_ESP_I2x3yz_S_aa;
    Double I_ESP_I2x2y2z_Pz_aa = I_ESP_K2x2y3z_S_aa+ABZ*I_ESP_I2x2y2z_S_aa;
    Double I_ESP_I2xy3z_Pz_aa = I_ESP_K2xy4z_S_aa+ABZ*I_ESP_I2xy3z_S_aa;
    Double I_ESP_I2x4z_Pz_aa = I_ESP_K2x5z_S_aa+ABZ*I_ESP_I2x4z_S_aa;
    Double I_ESP_Ix5y_Pz_aa = I_ESP_Kx5yz_S_aa+ABZ*I_ESP_Ix5y_S_aa;
    Double I_ESP_Ix4yz_Pz_aa = I_ESP_Kx4y2z_S_aa+ABZ*I_ESP_Ix4yz_S_aa;
    Double I_ESP_Ix3y2z_Pz_aa = I_ESP_Kx3y3z_S_aa+ABZ*I_ESP_Ix3y2z_S_aa;
    Double I_ESP_Ix2y3z_Pz_aa = I_ESP_Kx2y4z_S_aa+ABZ*I_ESP_Ix2y3z_S_aa;
    Double I_ESP_Ixy4z_Pz_aa = I_ESP_Kxy5z_S_aa+ABZ*I_ESP_Ixy4z_S_aa;
    Double I_ESP_Ix5z_Pz_aa = I_ESP_Kx6z_S_aa+ABZ*I_ESP_Ix5z_S_aa;
    Double I_ESP_I6y_Pz_aa = I_ESP_K6yz_S_aa+ABZ*I_ESP_I6y_S_aa;
    Double I_ESP_I5yz_Pz_aa = I_ESP_K5y2z_S_aa+ABZ*I_ESP_I5yz_S_aa;
    Double I_ESP_I4y2z_Pz_aa = I_ESP_K4y3z_S_aa+ABZ*I_ESP_I4y2z_S_aa;
    Double I_ESP_I3y3z_Pz_aa = I_ESP_K3y4z_S_aa+ABZ*I_ESP_I3y3z_S_aa;
    Double I_ESP_I2y4z_Pz_aa = I_ESP_K2y5z_S_aa+ABZ*I_ESP_I2y4z_S_aa;
    Double I_ESP_Iy5z_Pz_aa = I_ESP_Ky6z_S_aa+ABZ*I_ESP_Iy5z_S_aa;
    Double I_ESP_I6z_Pz_aa = I_ESP_K7z_S_aa+ABZ*I_ESP_I6z_S_aa;

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
     * shell quartet name: SQ_ESP_I_D_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 56 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_aa
     * RHS shell quartet name: SQ_ESP_I_P_aa
     ************************************************************/
    Double I_ESP_I6x_D2x_aa = I_ESP_K7x_Px_aa+ABX*I_ESP_I6x_Px_aa;
    Double I_ESP_I5xy_D2x_aa = I_ESP_K6xy_Px_aa+ABX*I_ESP_I5xy_Px_aa;
    Double I_ESP_I5xz_D2x_aa = I_ESP_K6xz_Px_aa+ABX*I_ESP_I5xz_Px_aa;
    Double I_ESP_I4x2y_D2x_aa = I_ESP_K5x2y_Px_aa+ABX*I_ESP_I4x2y_Px_aa;
    Double I_ESP_I4xyz_D2x_aa = I_ESP_K5xyz_Px_aa+ABX*I_ESP_I4xyz_Px_aa;
    Double I_ESP_I4x2z_D2x_aa = I_ESP_K5x2z_Px_aa+ABX*I_ESP_I4x2z_Px_aa;
    Double I_ESP_I3x3y_D2x_aa = I_ESP_K4x3y_Px_aa+ABX*I_ESP_I3x3y_Px_aa;
    Double I_ESP_I3x2yz_D2x_aa = I_ESP_K4x2yz_Px_aa+ABX*I_ESP_I3x2yz_Px_aa;
    Double I_ESP_I3xy2z_D2x_aa = I_ESP_K4xy2z_Px_aa+ABX*I_ESP_I3xy2z_Px_aa;
    Double I_ESP_I3x3z_D2x_aa = I_ESP_K4x3z_Px_aa+ABX*I_ESP_I3x3z_Px_aa;
    Double I_ESP_I2x4y_D2x_aa = I_ESP_K3x4y_Px_aa+ABX*I_ESP_I2x4y_Px_aa;
    Double I_ESP_I2x3yz_D2x_aa = I_ESP_K3x3yz_Px_aa+ABX*I_ESP_I2x3yz_Px_aa;
    Double I_ESP_I2x2y2z_D2x_aa = I_ESP_K3x2y2z_Px_aa+ABX*I_ESP_I2x2y2z_Px_aa;
    Double I_ESP_I2xy3z_D2x_aa = I_ESP_K3xy3z_Px_aa+ABX*I_ESP_I2xy3z_Px_aa;
    Double I_ESP_I2x4z_D2x_aa = I_ESP_K3x4z_Px_aa+ABX*I_ESP_I2x4z_Px_aa;
    Double I_ESP_Ix5y_D2x_aa = I_ESP_K2x5y_Px_aa+ABX*I_ESP_Ix5y_Px_aa;
    Double I_ESP_Ix4yz_D2x_aa = I_ESP_K2x4yz_Px_aa+ABX*I_ESP_Ix4yz_Px_aa;
    Double I_ESP_Ix3y2z_D2x_aa = I_ESP_K2x3y2z_Px_aa+ABX*I_ESP_Ix3y2z_Px_aa;
    Double I_ESP_Ix2y3z_D2x_aa = I_ESP_K2x2y3z_Px_aa+ABX*I_ESP_Ix2y3z_Px_aa;
    Double I_ESP_Ixy4z_D2x_aa = I_ESP_K2xy4z_Px_aa+ABX*I_ESP_Ixy4z_Px_aa;
    Double I_ESP_Ix5z_D2x_aa = I_ESP_K2x5z_Px_aa+ABX*I_ESP_Ix5z_Px_aa;
    Double I_ESP_I6y_D2x_aa = I_ESP_Kx6y_Px_aa+ABX*I_ESP_I6y_Px_aa;
    Double I_ESP_I5yz_D2x_aa = I_ESP_Kx5yz_Px_aa+ABX*I_ESP_I5yz_Px_aa;
    Double I_ESP_I4y2z_D2x_aa = I_ESP_Kx4y2z_Px_aa+ABX*I_ESP_I4y2z_Px_aa;
    Double I_ESP_I3y3z_D2x_aa = I_ESP_Kx3y3z_Px_aa+ABX*I_ESP_I3y3z_Px_aa;
    Double I_ESP_I2y4z_D2x_aa = I_ESP_Kx2y4z_Px_aa+ABX*I_ESP_I2y4z_Px_aa;
    Double I_ESP_Iy5z_D2x_aa = I_ESP_Kxy5z_Px_aa+ABX*I_ESP_Iy5z_Px_aa;
    Double I_ESP_I6z_D2x_aa = I_ESP_Kx6z_Px_aa+ABX*I_ESP_I6z_Px_aa;
    Double I_ESP_I6x_Dxy_aa = I_ESP_K6xy_Px_aa+ABY*I_ESP_I6x_Px_aa;
    Double I_ESP_I5xy_Dxy_aa = I_ESP_K5x2y_Px_aa+ABY*I_ESP_I5xy_Px_aa;
    Double I_ESP_I5xz_Dxy_aa = I_ESP_K5xyz_Px_aa+ABY*I_ESP_I5xz_Px_aa;
    Double I_ESP_I4x2y_Dxy_aa = I_ESP_K4x3y_Px_aa+ABY*I_ESP_I4x2y_Px_aa;
    Double I_ESP_I4xyz_Dxy_aa = I_ESP_K4x2yz_Px_aa+ABY*I_ESP_I4xyz_Px_aa;
    Double I_ESP_I4x2z_Dxy_aa = I_ESP_K4xy2z_Px_aa+ABY*I_ESP_I4x2z_Px_aa;
    Double I_ESP_I3x3y_Dxy_aa = I_ESP_K3x4y_Px_aa+ABY*I_ESP_I3x3y_Px_aa;
    Double I_ESP_I3x2yz_Dxy_aa = I_ESP_K3x3yz_Px_aa+ABY*I_ESP_I3x2yz_Px_aa;
    Double I_ESP_I3xy2z_Dxy_aa = I_ESP_K3x2y2z_Px_aa+ABY*I_ESP_I3xy2z_Px_aa;
    Double I_ESP_I3x3z_Dxy_aa = I_ESP_K3xy3z_Px_aa+ABY*I_ESP_I3x3z_Px_aa;
    Double I_ESP_I2x4y_Dxy_aa = I_ESP_K2x5y_Px_aa+ABY*I_ESP_I2x4y_Px_aa;
    Double I_ESP_I2x3yz_Dxy_aa = I_ESP_K2x4yz_Px_aa+ABY*I_ESP_I2x3yz_Px_aa;
    Double I_ESP_I2x2y2z_Dxy_aa = I_ESP_K2x3y2z_Px_aa+ABY*I_ESP_I2x2y2z_Px_aa;
    Double I_ESP_I2xy3z_Dxy_aa = I_ESP_K2x2y3z_Px_aa+ABY*I_ESP_I2xy3z_Px_aa;
    Double I_ESP_I2x4z_Dxy_aa = I_ESP_K2xy4z_Px_aa+ABY*I_ESP_I2x4z_Px_aa;
    Double I_ESP_Ix5y_Dxy_aa = I_ESP_Kx6y_Px_aa+ABY*I_ESP_Ix5y_Px_aa;
    Double I_ESP_Ix4yz_Dxy_aa = I_ESP_Kx5yz_Px_aa+ABY*I_ESP_Ix4yz_Px_aa;
    Double I_ESP_Ix3y2z_Dxy_aa = I_ESP_Kx4y2z_Px_aa+ABY*I_ESP_Ix3y2z_Px_aa;
    Double I_ESP_Ix2y3z_Dxy_aa = I_ESP_Kx3y3z_Px_aa+ABY*I_ESP_Ix2y3z_Px_aa;
    Double I_ESP_Ixy4z_Dxy_aa = I_ESP_Kx2y4z_Px_aa+ABY*I_ESP_Ixy4z_Px_aa;
    Double I_ESP_Ix5z_Dxy_aa = I_ESP_Kxy5z_Px_aa+ABY*I_ESP_Ix5z_Px_aa;
    Double I_ESP_I6y_Dxy_aa = I_ESP_K7y_Px_aa+ABY*I_ESP_I6y_Px_aa;
    Double I_ESP_I5yz_Dxy_aa = I_ESP_K6yz_Px_aa+ABY*I_ESP_I5yz_Px_aa;
    Double I_ESP_I4y2z_Dxy_aa = I_ESP_K5y2z_Px_aa+ABY*I_ESP_I4y2z_Px_aa;
    Double I_ESP_I3y3z_Dxy_aa = I_ESP_K4y3z_Px_aa+ABY*I_ESP_I3y3z_Px_aa;
    Double I_ESP_I2y4z_Dxy_aa = I_ESP_K3y4z_Px_aa+ABY*I_ESP_I2y4z_Px_aa;
    Double I_ESP_Iy5z_Dxy_aa = I_ESP_K2y5z_Px_aa+ABY*I_ESP_Iy5z_Px_aa;
    Double I_ESP_I6z_Dxy_aa = I_ESP_Ky6z_Px_aa+ABY*I_ESP_I6z_Px_aa;
    Double I_ESP_I6x_D2y_aa = I_ESP_K6xy_Py_aa+ABY*I_ESP_I6x_Py_aa;
    Double I_ESP_I5xy_D2y_aa = I_ESP_K5x2y_Py_aa+ABY*I_ESP_I5xy_Py_aa;
    Double I_ESP_I5xz_D2y_aa = I_ESP_K5xyz_Py_aa+ABY*I_ESP_I5xz_Py_aa;
    Double I_ESP_I4x2y_D2y_aa = I_ESP_K4x3y_Py_aa+ABY*I_ESP_I4x2y_Py_aa;
    Double I_ESP_I4xyz_D2y_aa = I_ESP_K4x2yz_Py_aa+ABY*I_ESP_I4xyz_Py_aa;
    Double I_ESP_I4x2z_D2y_aa = I_ESP_K4xy2z_Py_aa+ABY*I_ESP_I4x2z_Py_aa;
    Double I_ESP_I3x3y_D2y_aa = I_ESP_K3x4y_Py_aa+ABY*I_ESP_I3x3y_Py_aa;
    Double I_ESP_I3x2yz_D2y_aa = I_ESP_K3x3yz_Py_aa+ABY*I_ESP_I3x2yz_Py_aa;
    Double I_ESP_I3xy2z_D2y_aa = I_ESP_K3x2y2z_Py_aa+ABY*I_ESP_I3xy2z_Py_aa;
    Double I_ESP_I3x3z_D2y_aa = I_ESP_K3xy3z_Py_aa+ABY*I_ESP_I3x3z_Py_aa;
    Double I_ESP_I2x4y_D2y_aa = I_ESP_K2x5y_Py_aa+ABY*I_ESP_I2x4y_Py_aa;
    Double I_ESP_I2x3yz_D2y_aa = I_ESP_K2x4yz_Py_aa+ABY*I_ESP_I2x3yz_Py_aa;
    Double I_ESP_I2x2y2z_D2y_aa = I_ESP_K2x3y2z_Py_aa+ABY*I_ESP_I2x2y2z_Py_aa;
    Double I_ESP_I2xy3z_D2y_aa = I_ESP_K2x2y3z_Py_aa+ABY*I_ESP_I2xy3z_Py_aa;
    Double I_ESP_I2x4z_D2y_aa = I_ESP_K2xy4z_Py_aa+ABY*I_ESP_I2x4z_Py_aa;
    Double I_ESP_Ix5y_D2y_aa = I_ESP_Kx6y_Py_aa+ABY*I_ESP_Ix5y_Py_aa;
    Double I_ESP_Ix4yz_D2y_aa = I_ESP_Kx5yz_Py_aa+ABY*I_ESP_Ix4yz_Py_aa;
    Double I_ESP_Ix3y2z_D2y_aa = I_ESP_Kx4y2z_Py_aa+ABY*I_ESP_Ix3y2z_Py_aa;
    Double I_ESP_Ix2y3z_D2y_aa = I_ESP_Kx3y3z_Py_aa+ABY*I_ESP_Ix2y3z_Py_aa;
    Double I_ESP_Ixy4z_D2y_aa = I_ESP_Kx2y4z_Py_aa+ABY*I_ESP_Ixy4z_Py_aa;
    Double I_ESP_Ix5z_D2y_aa = I_ESP_Kxy5z_Py_aa+ABY*I_ESP_Ix5z_Py_aa;
    Double I_ESP_I6y_D2y_aa = I_ESP_K7y_Py_aa+ABY*I_ESP_I6y_Py_aa;
    Double I_ESP_I5yz_D2y_aa = I_ESP_K6yz_Py_aa+ABY*I_ESP_I5yz_Py_aa;
    Double I_ESP_I4y2z_D2y_aa = I_ESP_K5y2z_Py_aa+ABY*I_ESP_I4y2z_Py_aa;
    Double I_ESP_I3y3z_D2y_aa = I_ESP_K4y3z_Py_aa+ABY*I_ESP_I3y3z_Py_aa;
    Double I_ESP_I2y4z_D2y_aa = I_ESP_K3y4z_Py_aa+ABY*I_ESP_I2y4z_Py_aa;
    Double I_ESP_Iy5z_D2y_aa = I_ESP_K2y5z_Py_aa+ABY*I_ESP_Iy5z_Py_aa;
    Double I_ESP_I6z_D2y_aa = I_ESP_Ky6z_Py_aa+ABY*I_ESP_I6z_Py_aa;
    Double I_ESP_I6x_D2z_aa = I_ESP_K6xz_Pz_aa+ABZ*I_ESP_I6x_Pz_aa;
    Double I_ESP_I5xy_D2z_aa = I_ESP_K5xyz_Pz_aa+ABZ*I_ESP_I5xy_Pz_aa;
    Double I_ESP_I5xz_D2z_aa = I_ESP_K5x2z_Pz_aa+ABZ*I_ESP_I5xz_Pz_aa;
    Double I_ESP_I4x2y_D2z_aa = I_ESP_K4x2yz_Pz_aa+ABZ*I_ESP_I4x2y_Pz_aa;
    Double I_ESP_I4xyz_D2z_aa = I_ESP_K4xy2z_Pz_aa+ABZ*I_ESP_I4xyz_Pz_aa;
    Double I_ESP_I4x2z_D2z_aa = I_ESP_K4x3z_Pz_aa+ABZ*I_ESP_I4x2z_Pz_aa;
    Double I_ESP_I3x3y_D2z_aa = I_ESP_K3x3yz_Pz_aa+ABZ*I_ESP_I3x3y_Pz_aa;
    Double I_ESP_I3x2yz_D2z_aa = I_ESP_K3x2y2z_Pz_aa+ABZ*I_ESP_I3x2yz_Pz_aa;
    Double I_ESP_I3xy2z_D2z_aa = I_ESP_K3xy3z_Pz_aa+ABZ*I_ESP_I3xy2z_Pz_aa;
    Double I_ESP_I3x3z_D2z_aa = I_ESP_K3x4z_Pz_aa+ABZ*I_ESP_I3x3z_Pz_aa;
    Double I_ESP_I2x4y_D2z_aa = I_ESP_K2x4yz_Pz_aa+ABZ*I_ESP_I2x4y_Pz_aa;
    Double I_ESP_I2x3yz_D2z_aa = I_ESP_K2x3y2z_Pz_aa+ABZ*I_ESP_I2x3yz_Pz_aa;
    Double I_ESP_I2x2y2z_D2z_aa = I_ESP_K2x2y3z_Pz_aa+ABZ*I_ESP_I2x2y2z_Pz_aa;
    Double I_ESP_I2xy3z_D2z_aa = I_ESP_K2xy4z_Pz_aa+ABZ*I_ESP_I2xy3z_Pz_aa;
    Double I_ESP_I2x4z_D2z_aa = I_ESP_K2x5z_Pz_aa+ABZ*I_ESP_I2x4z_Pz_aa;
    Double I_ESP_Ix5y_D2z_aa = I_ESP_Kx5yz_Pz_aa+ABZ*I_ESP_Ix5y_Pz_aa;
    Double I_ESP_Ix4yz_D2z_aa = I_ESP_Kx4y2z_Pz_aa+ABZ*I_ESP_Ix4yz_Pz_aa;
    Double I_ESP_Ix3y2z_D2z_aa = I_ESP_Kx3y3z_Pz_aa+ABZ*I_ESP_Ix3y2z_Pz_aa;
    Double I_ESP_Ix2y3z_D2z_aa = I_ESP_Kx2y4z_Pz_aa+ABZ*I_ESP_Ix2y3z_Pz_aa;
    Double I_ESP_Ixy4z_D2z_aa = I_ESP_Kxy5z_Pz_aa+ABZ*I_ESP_Ixy4z_Pz_aa;
    Double I_ESP_Ix5z_D2z_aa = I_ESP_Kx6z_Pz_aa+ABZ*I_ESP_Ix5z_Pz_aa;
    Double I_ESP_I6y_D2z_aa = I_ESP_K6yz_Pz_aa+ABZ*I_ESP_I6y_Pz_aa;
    Double I_ESP_I5yz_D2z_aa = I_ESP_K5y2z_Pz_aa+ABZ*I_ESP_I5yz_Pz_aa;
    Double I_ESP_I4y2z_D2z_aa = I_ESP_K4y3z_Pz_aa+ABZ*I_ESP_I4y2z_Pz_aa;
    Double I_ESP_I3y3z_D2z_aa = I_ESP_K3y4z_Pz_aa+ABZ*I_ESP_I3y3z_Pz_aa;
    Double I_ESP_I2y4z_D2z_aa = I_ESP_K2y5z_Pz_aa+ABZ*I_ESP_I2y4z_Pz_aa;
    Double I_ESP_Iy5z_D2z_aa = I_ESP_Ky6z_Pz_aa+ABZ*I_ESP_Iy5z_Pz_aa;
    Double I_ESP_I6z_D2z_aa = I_ESP_K7z_Pz_aa+ABZ*I_ESP_I6z_Pz_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_L_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 20 integrals are omitted 
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
    Double I_ESP_L7yz_Px_aa = I_ESP_Mx7yz_S_aa+ABX*I_ESP_L7yz_S_aa;
    Double I_ESP_L6y2z_Px_aa = I_ESP_Mx6y2z_S_aa+ABX*I_ESP_L6y2z_S_aa;
    Double I_ESP_L5y3z_Px_aa = I_ESP_Mx5y3z_S_aa+ABX*I_ESP_L5y3z_S_aa;
    Double I_ESP_L4y4z_Px_aa = I_ESP_Mx4y4z_S_aa+ABX*I_ESP_L4y4z_S_aa;
    Double I_ESP_L3y5z_Px_aa = I_ESP_Mx3y5z_S_aa+ABX*I_ESP_L3y5z_S_aa;
    Double I_ESP_L2y6z_Px_aa = I_ESP_Mx2y6z_S_aa+ABX*I_ESP_L2y6z_S_aa;
    Double I_ESP_Ly7z_Px_aa = I_ESP_Mxy7z_S_aa+ABX*I_ESP_Ly7z_S_aa;
    Double I_ESP_L7xy_Py_aa = I_ESP_M7x2y_S_aa+ABY*I_ESP_L7xy_S_aa;
    Double I_ESP_L6x2y_Py_aa = I_ESP_M6x3y_S_aa+ABY*I_ESP_L6x2y_S_aa;
    Double I_ESP_L6xyz_Py_aa = I_ESP_M6x2yz_S_aa+ABY*I_ESP_L6xyz_S_aa;
    Double I_ESP_L5x3y_Py_aa = I_ESP_M5x4y_S_aa+ABY*I_ESP_L5x3y_S_aa;
    Double I_ESP_L5x2yz_Py_aa = I_ESP_M5x3yz_S_aa+ABY*I_ESP_L5x2yz_S_aa;
    Double I_ESP_L5xy2z_Py_aa = I_ESP_M5x2y2z_S_aa+ABY*I_ESP_L5xy2z_S_aa;
    Double I_ESP_L4x4y_Py_aa = I_ESP_M4x5y_S_aa+ABY*I_ESP_L4x4y_S_aa;
    Double I_ESP_L4x3yz_Py_aa = I_ESP_M4x4yz_S_aa+ABY*I_ESP_L4x3yz_S_aa;
    Double I_ESP_L4x2y2z_Py_aa = I_ESP_M4x3y2z_S_aa+ABY*I_ESP_L4x2y2z_S_aa;
    Double I_ESP_L4xy3z_Py_aa = I_ESP_M4x2y3z_S_aa+ABY*I_ESP_L4xy3z_S_aa;
    Double I_ESP_L3x5y_Py_aa = I_ESP_M3x6y_S_aa+ABY*I_ESP_L3x5y_S_aa;
    Double I_ESP_L3x4yz_Py_aa = I_ESP_M3x5yz_S_aa+ABY*I_ESP_L3x4yz_S_aa;
    Double I_ESP_L3x3y2z_Py_aa = I_ESP_M3x4y2z_S_aa+ABY*I_ESP_L3x3y2z_S_aa;
    Double I_ESP_L3x2y3z_Py_aa = I_ESP_M3x3y3z_S_aa+ABY*I_ESP_L3x2y3z_S_aa;
    Double I_ESP_L3xy4z_Py_aa = I_ESP_M3x2y4z_S_aa+ABY*I_ESP_L3xy4z_S_aa;
    Double I_ESP_L2x6y_Py_aa = I_ESP_M2x7y_S_aa+ABY*I_ESP_L2x6y_S_aa;
    Double I_ESP_L2x5yz_Py_aa = I_ESP_M2x6yz_S_aa+ABY*I_ESP_L2x5yz_S_aa;
    Double I_ESP_L2x4y2z_Py_aa = I_ESP_M2x5y2z_S_aa+ABY*I_ESP_L2x4y2z_S_aa;
    Double I_ESP_L2x3y3z_Py_aa = I_ESP_M2x4y3z_S_aa+ABY*I_ESP_L2x3y3z_S_aa;
    Double I_ESP_L2x2y4z_Py_aa = I_ESP_M2x3y4z_S_aa+ABY*I_ESP_L2x2y4z_S_aa;
    Double I_ESP_L2xy5z_Py_aa = I_ESP_M2x2y5z_S_aa+ABY*I_ESP_L2xy5z_S_aa;
    Double I_ESP_Lx7y_Py_aa = I_ESP_Mx8y_S_aa+ABY*I_ESP_Lx7y_S_aa;
    Double I_ESP_Lx6yz_Py_aa = I_ESP_Mx7yz_S_aa+ABY*I_ESP_Lx6yz_S_aa;
    Double I_ESP_Lx5y2z_Py_aa = I_ESP_Mx6y2z_S_aa+ABY*I_ESP_Lx5y2z_S_aa;
    Double I_ESP_Lx4y3z_Py_aa = I_ESP_Mx5y3z_S_aa+ABY*I_ESP_Lx4y3z_S_aa;
    Double I_ESP_Lx3y4z_Py_aa = I_ESP_Mx4y4z_S_aa+ABY*I_ESP_Lx3y4z_S_aa;
    Double I_ESP_Lx2y5z_Py_aa = I_ESP_Mx3y5z_S_aa+ABY*I_ESP_Lx2y5z_S_aa;
    Double I_ESP_Lxy6z_Py_aa = I_ESP_Mx2y6z_S_aa+ABY*I_ESP_Lxy6z_S_aa;
    Double I_ESP_L8y_Py_aa = I_ESP_M9y_S_aa+ABY*I_ESP_L8y_S_aa;
    Double I_ESP_L7yz_Py_aa = I_ESP_M8yz_S_aa+ABY*I_ESP_L7yz_S_aa;
    Double I_ESP_L6y2z_Py_aa = I_ESP_M7y2z_S_aa+ABY*I_ESP_L6y2z_S_aa;
    Double I_ESP_L5y3z_Py_aa = I_ESP_M6y3z_S_aa+ABY*I_ESP_L5y3z_S_aa;
    Double I_ESP_L4y4z_Py_aa = I_ESP_M5y4z_S_aa+ABY*I_ESP_L4y4z_S_aa;
    Double I_ESP_L3y5z_Py_aa = I_ESP_M4y5z_S_aa+ABY*I_ESP_L3y5z_S_aa;
    Double I_ESP_L2y6z_Py_aa = I_ESP_M3y6z_S_aa+ABY*I_ESP_L2y6z_S_aa;
    Double I_ESP_Ly7z_Py_aa = I_ESP_M2y7z_S_aa+ABY*I_ESP_Ly7z_S_aa;
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
     * totally 80 integrals are omitted 
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
    Double I_ESP_K6xz_Dxy_aa = I_ESP_L6xyz_Px_aa+ABY*I_ESP_K6xz_Px_aa;
    Double I_ESP_K5xyz_Dxy_aa = I_ESP_L5x2yz_Px_aa+ABY*I_ESP_K5xyz_Px_aa;
    Double I_ESP_K5x2z_Dxy_aa = I_ESP_L5xy2z_Px_aa+ABY*I_ESP_K5x2z_Px_aa;
    Double I_ESP_K4x2yz_Dxy_aa = I_ESP_L4x3yz_Px_aa+ABY*I_ESP_K4x2yz_Px_aa;
    Double I_ESP_K4xy2z_Dxy_aa = I_ESP_L4x2y2z_Px_aa+ABY*I_ESP_K4xy2z_Px_aa;
    Double I_ESP_K4x3z_Dxy_aa = I_ESP_L4xy3z_Px_aa+ABY*I_ESP_K4x3z_Px_aa;
    Double I_ESP_K3x3yz_Dxy_aa = I_ESP_L3x4yz_Px_aa+ABY*I_ESP_K3x3yz_Px_aa;
    Double I_ESP_K3x2y2z_Dxy_aa = I_ESP_L3x3y2z_Px_aa+ABY*I_ESP_K3x2y2z_Px_aa;
    Double I_ESP_K3xy3z_Dxy_aa = I_ESP_L3x2y3z_Px_aa+ABY*I_ESP_K3xy3z_Px_aa;
    Double I_ESP_K3x4z_Dxy_aa = I_ESP_L3xy4z_Px_aa+ABY*I_ESP_K3x4z_Px_aa;
    Double I_ESP_K2x4yz_Dxy_aa = I_ESP_L2x5yz_Px_aa+ABY*I_ESP_K2x4yz_Px_aa;
    Double I_ESP_K2x3y2z_Dxy_aa = I_ESP_L2x4y2z_Px_aa+ABY*I_ESP_K2x3y2z_Px_aa;
    Double I_ESP_K2x2y3z_Dxy_aa = I_ESP_L2x3y3z_Px_aa+ABY*I_ESP_K2x2y3z_Px_aa;
    Double I_ESP_K2xy4z_Dxy_aa = I_ESP_L2x2y4z_Px_aa+ABY*I_ESP_K2xy4z_Px_aa;
    Double I_ESP_K2x5z_Dxy_aa = I_ESP_L2xy5z_Px_aa+ABY*I_ESP_K2x5z_Px_aa;
    Double I_ESP_Kx5yz_Dxy_aa = I_ESP_Lx6yz_Px_aa+ABY*I_ESP_Kx5yz_Px_aa;
    Double I_ESP_Kx4y2z_Dxy_aa = I_ESP_Lx5y2z_Px_aa+ABY*I_ESP_Kx4y2z_Px_aa;
    Double I_ESP_Kx3y3z_Dxy_aa = I_ESP_Lx4y3z_Px_aa+ABY*I_ESP_Kx3y3z_Px_aa;
    Double I_ESP_Kx2y4z_Dxy_aa = I_ESP_Lx3y4z_Px_aa+ABY*I_ESP_Kx2y4z_Px_aa;
    Double I_ESP_Kxy5z_Dxy_aa = I_ESP_Lx2y5z_Px_aa+ABY*I_ESP_Kxy5z_Px_aa;
    Double I_ESP_Kx6z_Dxy_aa = I_ESP_Lxy6z_Px_aa+ABY*I_ESP_Kx6z_Px_aa;
    Double I_ESP_K6yz_Dxy_aa = I_ESP_L7yz_Px_aa+ABY*I_ESP_K6yz_Px_aa;
    Double I_ESP_K5y2z_Dxy_aa = I_ESP_L6y2z_Px_aa+ABY*I_ESP_K5y2z_Px_aa;
    Double I_ESP_K4y3z_Dxy_aa = I_ESP_L5y3z_Px_aa+ABY*I_ESP_K4y3z_Px_aa;
    Double I_ESP_K3y4z_Dxy_aa = I_ESP_L4y4z_Px_aa+ABY*I_ESP_K3y4z_Px_aa;
    Double I_ESP_K2y5z_Dxy_aa = I_ESP_L3y5z_Px_aa+ABY*I_ESP_K2y5z_Px_aa;
    Double I_ESP_Ky6z_Dxy_aa = I_ESP_L2y6z_Px_aa+ABY*I_ESP_Ky6z_Px_aa;
    Double I_ESP_K7z_Dxy_aa = I_ESP_Ly7z_Px_aa+ABY*I_ESP_K7z_Px_aa;
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
     * shell quartet name: SQ_ESP_I_F_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D_aa
     * RHS shell quartet name: SQ_ESP_I_D_aa
     ************************************************************/
    Double I_ESP_I6x_F3x_aa = I_ESP_K7x_D2x_aa+ABX*I_ESP_I6x_D2x_aa;
    Double I_ESP_I5xy_F3x_aa = I_ESP_K6xy_D2x_aa+ABX*I_ESP_I5xy_D2x_aa;
    Double I_ESP_I5xz_F3x_aa = I_ESP_K6xz_D2x_aa+ABX*I_ESP_I5xz_D2x_aa;
    Double I_ESP_I4x2y_F3x_aa = I_ESP_K5x2y_D2x_aa+ABX*I_ESP_I4x2y_D2x_aa;
    Double I_ESP_I4xyz_F3x_aa = I_ESP_K5xyz_D2x_aa+ABX*I_ESP_I4xyz_D2x_aa;
    Double I_ESP_I4x2z_F3x_aa = I_ESP_K5x2z_D2x_aa+ABX*I_ESP_I4x2z_D2x_aa;
    Double I_ESP_I3x3y_F3x_aa = I_ESP_K4x3y_D2x_aa+ABX*I_ESP_I3x3y_D2x_aa;
    Double I_ESP_I3x2yz_F3x_aa = I_ESP_K4x2yz_D2x_aa+ABX*I_ESP_I3x2yz_D2x_aa;
    Double I_ESP_I3xy2z_F3x_aa = I_ESP_K4xy2z_D2x_aa+ABX*I_ESP_I3xy2z_D2x_aa;
    Double I_ESP_I3x3z_F3x_aa = I_ESP_K4x3z_D2x_aa+ABX*I_ESP_I3x3z_D2x_aa;
    Double I_ESP_I2x4y_F3x_aa = I_ESP_K3x4y_D2x_aa+ABX*I_ESP_I2x4y_D2x_aa;
    Double I_ESP_I2x3yz_F3x_aa = I_ESP_K3x3yz_D2x_aa+ABX*I_ESP_I2x3yz_D2x_aa;
    Double I_ESP_I2x2y2z_F3x_aa = I_ESP_K3x2y2z_D2x_aa+ABX*I_ESP_I2x2y2z_D2x_aa;
    Double I_ESP_I2xy3z_F3x_aa = I_ESP_K3xy3z_D2x_aa+ABX*I_ESP_I2xy3z_D2x_aa;
    Double I_ESP_I2x4z_F3x_aa = I_ESP_K3x4z_D2x_aa+ABX*I_ESP_I2x4z_D2x_aa;
    Double I_ESP_Ix5y_F3x_aa = I_ESP_K2x5y_D2x_aa+ABX*I_ESP_Ix5y_D2x_aa;
    Double I_ESP_Ix4yz_F3x_aa = I_ESP_K2x4yz_D2x_aa+ABX*I_ESP_Ix4yz_D2x_aa;
    Double I_ESP_Ix3y2z_F3x_aa = I_ESP_K2x3y2z_D2x_aa+ABX*I_ESP_Ix3y2z_D2x_aa;
    Double I_ESP_Ix2y3z_F3x_aa = I_ESP_K2x2y3z_D2x_aa+ABX*I_ESP_Ix2y3z_D2x_aa;
    Double I_ESP_Ixy4z_F3x_aa = I_ESP_K2xy4z_D2x_aa+ABX*I_ESP_Ixy4z_D2x_aa;
    Double I_ESP_Ix5z_F3x_aa = I_ESP_K2x5z_D2x_aa+ABX*I_ESP_Ix5z_D2x_aa;
    Double I_ESP_I6y_F3x_aa = I_ESP_Kx6y_D2x_aa+ABX*I_ESP_I6y_D2x_aa;
    Double I_ESP_I5yz_F3x_aa = I_ESP_Kx5yz_D2x_aa+ABX*I_ESP_I5yz_D2x_aa;
    Double I_ESP_I4y2z_F3x_aa = I_ESP_Kx4y2z_D2x_aa+ABX*I_ESP_I4y2z_D2x_aa;
    Double I_ESP_I3y3z_F3x_aa = I_ESP_Kx3y3z_D2x_aa+ABX*I_ESP_I3y3z_D2x_aa;
    Double I_ESP_I2y4z_F3x_aa = I_ESP_Kx2y4z_D2x_aa+ABX*I_ESP_I2y4z_D2x_aa;
    Double I_ESP_Iy5z_F3x_aa = I_ESP_Kxy5z_D2x_aa+ABX*I_ESP_Iy5z_D2x_aa;
    Double I_ESP_I6z_F3x_aa = I_ESP_Kx6z_D2x_aa+ABX*I_ESP_I6z_D2x_aa;
    Double I_ESP_I6x_F2xy_aa = I_ESP_K6xy_D2x_aa+ABY*I_ESP_I6x_D2x_aa;
    Double I_ESP_I5xy_F2xy_aa = I_ESP_K5x2y_D2x_aa+ABY*I_ESP_I5xy_D2x_aa;
    Double I_ESP_I5xz_F2xy_aa = I_ESP_K5xyz_D2x_aa+ABY*I_ESP_I5xz_D2x_aa;
    Double I_ESP_I4x2y_F2xy_aa = I_ESP_K4x3y_D2x_aa+ABY*I_ESP_I4x2y_D2x_aa;
    Double I_ESP_I4xyz_F2xy_aa = I_ESP_K4x2yz_D2x_aa+ABY*I_ESP_I4xyz_D2x_aa;
    Double I_ESP_I4x2z_F2xy_aa = I_ESP_K4xy2z_D2x_aa+ABY*I_ESP_I4x2z_D2x_aa;
    Double I_ESP_I3x3y_F2xy_aa = I_ESP_K3x4y_D2x_aa+ABY*I_ESP_I3x3y_D2x_aa;
    Double I_ESP_I3x2yz_F2xy_aa = I_ESP_K3x3yz_D2x_aa+ABY*I_ESP_I3x2yz_D2x_aa;
    Double I_ESP_I3xy2z_F2xy_aa = I_ESP_K3x2y2z_D2x_aa+ABY*I_ESP_I3xy2z_D2x_aa;
    Double I_ESP_I3x3z_F2xy_aa = I_ESP_K3xy3z_D2x_aa+ABY*I_ESP_I3x3z_D2x_aa;
    Double I_ESP_I2x4y_F2xy_aa = I_ESP_K2x5y_D2x_aa+ABY*I_ESP_I2x4y_D2x_aa;
    Double I_ESP_I2x3yz_F2xy_aa = I_ESP_K2x4yz_D2x_aa+ABY*I_ESP_I2x3yz_D2x_aa;
    Double I_ESP_I2x2y2z_F2xy_aa = I_ESP_K2x3y2z_D2x_aa+ABY*I_ESP_I2x2y2z_D2x_aa;
    Double I_ESP_I2xy3z_F2xy_aa = I_ESP_K2x2y3z_D2x_aa+ABY*I_ESP_I2xy3z_D2x_aa;
    Double I_ESP_I2x4z_F2xy_aa = I_ESP_K2xy4z_D2x_aa+ABY*I_ESP_I2x4z_D2x_aa;
    Double I_ESP_Ix5y_F2xy_aa = I_ESP_Kx6y_D2x_aa+ABY*I_ESP_Ix5y_D2x_aa;
    Double I_ESP_Ix4yz_F2xy_aa = I_ESP_Kx5yz_D2x_aa+ABY*I_ESP_Ix4yz_D2x_aa;
    Double I_ESP_Ix3y2z_F2xy_aa = I_ESP_Kx4y2z_D2x_aa+ABY*I_ESP_Ix3y2z_D2x_aa;
    Double I_ESP_Ix2y3z_F2xy_aa = I_ESP_Kx3y3z_D2x_aa+ABY*I_ESP_Ix2y3z_D2x_aa;
    Double I_ESP_Ixy4z_F2xy_aa = I_ESP_Kx2y4z_D2x_aa+ABY*I_ESP_Ixy4z_D2x_aa;
    Double I_ESP_Ix5z_F2xy_aa = I_ESP_Kxy5z_D2x_aa+ABY*I_ESP_Ix5z_D2x_aa;
    Double I_ESP_I6y_F2xy_aa = I_ESP_K7y_D2x_aa+ABY*I_ESP_I6y_D2x_aa;
    Double I_ESP_I5yz_F2xy_aa = I_ESP_K6yz_D2x_aa+ABY*I_ESP_I5yz_D2x_aa;
    Double I_ESP_I4y2z_F2xy_aa = I_ESP_K5y2z_D2x_aa+ABY*I_ESP_I4y2z_D2x_aa;
    Double I_ESP_I3y3z_F2xy_aa = I_ESP_K4y3z_D2x_aa+ABY*I_ESP_I3y3z_D2x_aa;
    Double I_ESP_I2y4z_F2xy_aa = I_ESP_K3y4z_D2x_aa+ABY*I_ESP_I2y4z_D2x_aa;
    Double I_ESP_Iy5z_F2xy_aa = I_ESP_K2y5z_D2x_aa+ABY*I_ESP_Iy5z_D2x_aa;
    Double I_ESP_I6z_F2xy_aa = I_ESP_Ky6z_D2x_aa+ABY*I_ESP_I6z_D2x_aa;
    Double I_ESP_I6x_F2xz_aa = I_ESP_K6xz_D2x_aa+ABZ*I_ESP_I6x_D2x_aa;
    Double I_ESP_I5xy_F2xz_aa = I_ESP_K5xyz_D2x_aa+ABZ*I_ESP_I5xy_D2x_aa;
    Double I_ESP_I5xz_F2xz_aa = I_ESP_K5x2z_D2x_aa+ABZ*I_ESP_I5xz_D2x_aa;
    Double I_ESP_I4x2y_F2xz_aa = I_ESP_K4x2yz_D2x_aa+ABZ*I_ESP_I4x2y_D2x_aa;
    Double I_ESP_I4xyz_F2xz_aa = I_ESP_K4xy2z_D2x_aa+ABZ*I_ESP_I4xyz_D2x_aa;
    Double I_ESP_I4x2z_F2xz_aa = I_ESP_K4x3z_D2x_aa+ABZ*I_ESP_I4x2z_D2x_aa;
    Double I_ESP_I3x3y_F2xz_aa = I_ESP_K3x3yz_D2x_aa+ABZ*I_ESP_I3x3y_D2x_aa;
    Double I_ESP_I3x2yz_F2xz_aa = I_ESP_K3x2y2z_D2x_aa+ABZ*I_ESP_I3x2yz_D2x_aa;
    Double I_ESP_I3xy2z_F2xz_aa = I_ESP_K3xy3z_D2x_aa+ABZ*I_ESP_I3xy2z_D2x_aa;
    Double I_ESP_I3x3z_F2xz_aa = I_ESP_K3x4z_D2x_aa+ABZ*I_ESP_I3x3z_D2x_aa;
    Double I_ESP_I2x4y_F2xz_aa = I_ESP_K2x4yz_D2x_aa+ABZ*I_ESP_I2x4y_D2x_aa;
    Double I_ESP_I2x3yz_F2xz_aa = I_ESP_K2x3y2z_D2x_aa+ABZ*I_ESP_I2x3yz_D2x_aa;
    Double I_ESP_I2x2y2z_F2xz_aa = I_ESP_K2x2y3z_D2x_aa+ABZ*I_ESP_I2x2y2z_D2x_aa;
    Double I_ESP_I2xy3z_F2xz_aa = I_ESP_K2xy4z_D2x_aa+ABZ*I_ESP_I2xy3z_D2x_aa;
    Double I_ESP_I2x4z_F2xz_aa = I_ESP_K2x5z_D2x_aa+ABZ*I_ESP_I2x4z_D2x_aa;
    Double I_ESP_Ix5y_F2xz_aa = I_ESP_Kx5yz_D2x_aa+ABZ*I_ESP_Ix5y_D2x_aa;
    Double I_ESP_Ix4yz_F2xz_aa = I_ESP_Kx4y2z_D2x_aa+ABZ*I_ESP_Ix4yz_D2x_aa;
    Double I_ESP_Ix3y2z_F2xz_aa = I_ESP_Kx3y3z_D2x_aa+ABZ*I_ESP_Ix3y2z_D2x_aa;
    Double I_ESP_Ix2y3z_F2xz_aa = I_ESP_Kx2y4z_D2x_aa+ABZ*I_ESP_Ix2y3z_D2x_aa;
    Double I_ESP_Ixy4z_F2xz_aa = I_ESP_Kxy5z_D2x_aa+ABZ*I_ESP_Ixy4z_D2x_aa;
    Double I_ESP_Ix5z_F2xz_aa = I_ESP_Kx6z_D2x_aa+ABZ*I_ESP_Ix5z_D2x_aa;
    Double I_ESP_I6y_F2xz_aa = I_ESP_K6yz_D2x_aa+ABZ*I_ESP_I6y_D2x_aa;
    Double I_ESP_I5yz_F2xz_aa = I_ESP_K5y2z_D2x_aa+ABZ*I_ESP_I5yz_D2x_aa;
    Double I_ESP_I4y2z_F2xz_aa = I_ESP_K4y3z_D2x_aa+ABZ*I_ESP_I4y2z_D2x_aa;
    Double I_ESP_I3y3z_F2xz_aa = I_ESP_K3y4z_D2x_aa+ABZ*I_ESP_I3y3z_D2x_aa;
    Double I_ESP_I2y4z_F2xz_aa = I_ESP_K2y5z_D2x_aa+ABZ*I_ESP_I2y4z_D2x_aa;
    Double I_ESP_Iy5z_F2xz_aa = I_ESP_Ky6z_D2x_aa+ABZ*I_ESP_Iy5z_D2x_aa;
    Double I_ESP_I6z_F2xz_aa = I_ESP_K7z_D2x_aa+ABZ*I_ESP_I6z_D2x_aa;
    Double I_ESP_I6x_Fx2y_aa = I_ESP_K7x_D2y_aa+ABX*I_ESP_I6x_D2y_aa;
    Double I_ESP_I5xy_Fx2y_aa = I_ESP_K6xy_D2y_aa+ABX*I_ESP_I5xy_D2y_aa;
    Double I_ESP_I5xz_Fx2y_aa = I_ESP_K6xz_D2y_aa+ABX*I_ESP_I5xz_D2y_aa;
    Double I_ESP_I4x2y_Fx2y_aa = I_ESP_K5x2y_D2y_aa+ABX*I_ESP_I4x2y_D2y_aa;
    Double I_ESP_I4xyz_Fx2y_aa = I_ESP_K5xyz_D2y_aa+ABX*I_ESP_I4xyz_D2y_aa;
    Double I_ESP_I4x2z_Fx2y_aa = I_ESP_K5x2z_D2y_aa+ABX*I_ESP_I4x2z_D2y_aa;
    Double I_ESP_I3x3y_Fx2y_aa = I_ESP_K4x3y_D2y_aa+ABX*I_ESP_I3x3y_D2y_aa;
    Double I_ESP_I3x2yz_Fx2y_aa = I_ESP_K4x2yz_D2y_aa+ABX*I_ESP_I3x2yz_D2y_aa;
    Double I_ESP_I3xy2z_Fx2y_aa = I_ESP_K4xy2z_D2y_aa+ABX*I_ESP_I3xy2z_D2y_aa;
    Double I_ESP_I3x3z_Fx2y_aa = I_ESP_K4x3z_D2y_aa+ABX*I_ESP_I3x3z_D2y_aa;
    Double I_ESP_I2x4y_Fx2y_aa = I_ESP_K3x4y_D2y_aa+ABX*I_ESP_I2x4y_D2y_aa;
    Double I_ESP_I2x3yz_Fx2y_aa = I_ESP_K3x3yz_D2y_aa+ABX*I_ESP_I2x3yz_D2y_aa;
    Double I_ESP_I2x2y2z_Fx2y_aa = I_ESP_K3x2y2z_D2y_aa+ABX*I_ESP_I2x2y2z_D2y_aa;
    Double I_ESP_I2xy3z_Fx2y_aa = I_ESP_K3xy3z_D2y_aa+ABX*I_ESP_I2xy3z_D2y_aa;
    Double I_ESP_I2x4z_Fx2y_aa = I_ESP_K3x4z_D2y_aa+ABX*I_ESP_I2x4z_D2y_aa;
    Double I_ESP_Ix5y_Fx2y_aa = I_ESP_K2x5y_D2y_aa+ABX*I_ESP_Ix5y_D2y_aa;
    Double I_ESP_Ix4yz_Fx2y_aa = I_ESP_K2x4yz_D2y_aa+ABX*I_ESP_Ix4yz_D2y_aa;
    Double I_ESP_Ix3y2z_Fx2y_aa = I_ESP_K2x3y2z_D2y_aa+ABX*I_ESP_Ix3y2z_D2y_aa;
    Double I_ESP_Ix2y3z_Fx2y_aa = I_ESP_K2x2y3z_D2y_aa+ABX*I_ESP_Ix2y3z_D2y_aa;
    Double I_ESP_Ixy4z_Fx2y_aa = I_ESP_K2xy4z_D2y_aa+ABX*I_ESP_Ixy4z_D2y_aa;
    Double I_ESP_Ix5z_Fx2y_aa = I_ESP_K2x5z_D2y_aa+ABX*I_ESP_Ix5z_D2y_aa;
    Double I_ESP_I6y_Fx2y_aa = I_ESP_Kx6y_D2y_aa+ABX*I_ESP_I6y_D2y_aa;
    Double I_ESP_I5yz_Fx2y_aa = I_ESP_Kx5yz_D2y_aa+ABX*I_ESP_I5yz_D2y_aa;
    Double I_ESP_I4y2z_Fx2y_aa = I_ESP_Kx4y2z_D2y_aa+ABX*I_ESP_I4y2z_D2y_aa;
    Double I_ESP_I3y3z_Fx2y_aa = I_ESP_Kx3y3z_D2y_aa+ABX*I_ESP_I3y3z_D2y_aa;
    Double I_ESP_I2y4z_Fx2y_aa = I_ESP_Kx2y4z_D2y_aa+ABX*I_ESP_I2y4z_D2y_aa;
    Double I_ESP_Iy5z_Fx2y_aa = I_ESP_Kxy5z_D2y_aa+ABX*I_ESP_Iy5z_D2y_aa;
    Double I_ESP_I6z_Fx2y_aa = I_ESP_Kx6z_D2y_aa+ABX*I_ESP_I6z_D2y_aa;
    Double I_ESP_I6x_Fxyz_aa = I_ESP_K6xz_Dxy_aa+ABZ*I_ESP_I6x_Dxy_aa;
    Double I_ESP_I5xy_Fxyz_aa = I_ESP_K5xyz_Dxy_aa+ABZ*I_ESP_I5xy_Dxy_aa;
    Double I_ESP_I5xz_Fxyz_aa = I_ESP_K5x2z_Dxy_aa+ABZ*I_ESP_I5xz_Dxy_aa;
    Double I_ESP_I4x2y_Fxyz_aa = I_ESP_K4x2yz_Dxy_aa+ABZ*I_ESP_I4x2y_Dxy_aa;
    Double I_ESP_I4xyz_Fxyz_aa = I_ESP_K4xy2z_Dxy_aa+ABZ*I_ESP_I4xyz_Dxy_aa;
    Double I_ESP_I4x2z_Fxyz_aa = I_ESP_K4x3z_Dxy_aa+ABZ*I_ESP_I4x2z_Dxy_aa;
    Double I_ESP_I3x3y_Fxyz_aa = I_ESP_K3x3yz_Dxy_aa+ABZ*I_ESP_I3x3y_Dxy_aa;
    Double I_ESP_I3x2yz_Fxyz_aa = I_ESP_K3x2y2z_Dxy_aa+ABZ*I_ESP_I3x2yz_Dxy_aa;
    Double I_ESP_I3xy2z_Fxyz_aa = I_ESP_K3xy3z_Dxy_aa+ABZ*I_ESP_I3xy2z_Dxy_aa;
    Double I_ESP_I3x3z_Fxyz_aa = I_ESP_K3x4z_Dxy_aa+ABZ*I_ESP_I3x3z_Dxy_aa;
    Double I_ESP_I2x4y_Fxyz_aa = I_ESP_K2x4yz_Dxy_aa+ABZ*I_ESP_I2x4y_Dxy_aa;
    Double I_ESP_I2x3yz_Fxyz_aa = I_ESP_K2x3y2z_Dxy_aa+ABZ*I_ESP_I2x3yz_Dxy_aa;
    Double I_ESP_I2x2y2z_Fxyz_aa = I_ESP_K2x2y3z_Dxy_aa+ABZ*I_ESP_I2x2y2z_Dxy_aa;
    Double I_ESP_I2xy3z_Fxyz_aa = I_ESP_K2xy4z_Dxy_aa+ABZ*I_ESP_I2xy3z_Dxy_aa;
    Double I_ESP_I2x4z_Fxyz_aa = I_ESP_K2x5z_Dxy_aa+ABZ*I_ESP_I2x4z_Dxy_aa;
    Double I_ESP_Ix5y_Fxyz_aa = I_ESP_Kx5yz_Dxy_aa+ABZ*I_ESP_Ix5y_Dxy_aa;
    Double I_ESP_Ix4yz_Fxyz_aa = I_ESP_Kx4y2z_Dxy_aa+ABZ*I_ESP_Ix4yz_Dxy_aa;
    Double I_ESP_Ix3y2z_Fxyz_aa = I_ESP_Kx3y3z_Dxy_aa+ABZ*I_ESP_Ix3y2z_Dxy_aa;
    Double I_ESP_Ix2y3z_Fxyz_aa = I_ESP_Kx2y4z_Dxy_aa+ABZ*I_ESP_Ix2y3z_Dxy_aa;
    Double I_ESP_Ixy4z_Fxyz_aa = I_ESP_Kxy5z_Dxy_aa+ABZ*I_ESP_Ixy4z_Dxy_aa;
    Double I_ESP_Ix5z_Fxyz_aa = I_ESP_Kx6z_Dxy_aa+ABZ*I_ESP_Ix5z_Dxy_aa;
    Double I_ESP_I6y_Fxyz_aa = I_ESP_K6yz_Dxy_aa+ABZ*I_ESP_I6y_Dxy_aa;
    Double I_ESP_I5yz_Fxyz_aa = I_ESP_K5y2z_Dxy_aa+ABZ*I_ESP_I5yz_Dxy_aa;
    Double I_ESP_I4y2z_Fxyz_aa = I_ESP_K4y3z_Dxy_aa+ABZ*I_ESP_I4y2z_Dxy_aa;
    Double I_ESP_I3y3z_Fxyz_aa = I_ESP_K3y4z_Dxy_aa+ABZ*I_ESP_I3y3z_Dxy_aa;
    Double I_ESP_I2y4z_Fxyz_aa = I_ESP_K2y5z_Dxy_aa+ABZ*I_ESP_I2y4z_Dxy_aa;
    Double I_ESP_Iy5z_Fxyz_aa = I_ESP_Ky6z_Dxy_aa+ABZ*I_ESP_Iy5z_Dxy_aa;
    Double I_ESP_I6z_Fxyz_aa = I_ESP_K7z_Dxy_aa+ABZ*I_ESP_I6z_Dxy_aa;
    Double I_ESP_I6x_Fx2z_aa = I_ESP_K7x_D2z_aa+ABX*I_ESP_I6x_D2z_aa;
    Double I_ESP_I5xy_Fx2z_aa = I_ESP_K6xy_D2z_aa+ABX*I_ESP_I5xy_D2z_aa;
    Double I_ESP_I5xz_Fx2z_aa = I_ESP_K6xz_D2z_aa+ABX*I_ESP_I5xz_D2z_aa;
    Double I_ESP_I4x2y_Fx2z_aa = I_ESP_K5x2y_D2z_aa+ABX*I_ESP_I4x2y_D2z_aa;
    Double I_ESP_I4xyz_Fx2z_aa = I_ESP_K5xyz_D2z_aa+ABX*I_ESP_I4xyz_D2z_aa;
    Double I_ESP_I4x2z_Fx2z_aa = I_ESP_K5x2z_D2z_aa+ABX*I_ESP_I4x2z_D2z_aa;
    Double I_ESP_I3x3y_Fx2z_aa = I_ESP_K4x3y_D2z_aa+ABX*I_ESP_I3x3y_D2z_aa;
    Double I_ESP_I3x2yz_Fx2z_aa = I_ESP_K4x2yz_D2z_aa+ABX*I_ESP_I3x2yz_D2z_aa;
    Double I_ESP_I3xy2z_Fx2z_aa = I_ESP_K4xy2z_D2z_aa+ABX*I_ESP_I3xy2z_D2z_aa;
    Double I_ESP_I3x3z_Fx2z_aa = I_ESP_K4x3z_D2z_aa+ABX*I_ESP_I3x3z_D2z_aa;
    Double I_ESP_I2x4y_Fx2z_aa = I_ESP_K3x4y_D2z_aa+ABX*I_ESP_I2x4y_D2z_aa;
    Double I_ESP_I2x3yz_Fx2z_aa = I_ESP_K3x3yz_D2z_aa+ABX*I_ESP_I2x3yz_D2z_aa;
    Double I_ESP_I2x2y2z_Fx2z_aa = I_ESP_K3x2y2z_D2z_aa+ABX*I_ESP_I2x2y2z_D2z_aa;
    Double I_ESP_I2xy3z_Fx2z_aa = I_ESP_K3xy3z_D2z_aa+ABX*I_ESP_I2xy3z_D2z_aa;
    Double I_ESP_I2x4z_Fx2z_aa = I_ESP_K3x4z_D2z_aa+ABX*I_ESP_I2x4z_D2z_aa;
    Double I_ESP_Ix5y_Fx2z_aa = I_ESP_K2x5y_D2z_aa+ABX*I_ESP_Ix5y_D2z_aa;
    Double I_ESP_Ix4yz_Fx2z_aa = I_ESP_K2x4yz_D2z_aa+ABX*I_ESP_Ix4yz_D2z_aa;
    Double I_ESP_Ix3y2z_Fx2z_aa = I_ESP_K2x3y2z_D2z_aa+ABX*I_ESP_Ix3y2z_D2z_aa;
    Double I_ESP_Ix2y3z_Fx2z_aa = I_ESP_K2x2y3z_D2z_aa+ABX*I_ESP_Ix2y3z_D2z_aa;
    Double I_ESP_Ixy4z_Fx2z_aa = I_ESP_K2xy4z_D2z_aa+ABX*I_ESP_Ixy4z_D2z_aa;
    Double I_ESP_Ix5z_Fx2z_aa = I_ESP_K2x5z_D2z_aa+ABX*I_ESP_Ix5z_D2z_aa;
    Double I_ESP_I6y_Fx2z_aa = I_ESP_Kx6y_D2z_aa+ABX*I_ESP_I6y_D2z_aa;
    Double I_ESP_I5yz_Fx2z_aa = I_ESP_Kx5yz_D2z_aa+ABX*I_ESP_I5yz_D2z_aa;
    Double I_ESP_I4y2z_Fx2z_aa = I_ESP_Kx4y2z_D2z_aa+ABX*I_ESP_I4y2z_D2z_aa;
    Double I_ESP_I3y3z_Fx2z_aa = I_ESP_Kx3y3z_D2z_aa+ABX*I_ESP_I3y3z_D2z_aa;
    Double I_ESP_I2y4z_Fx2z_aa = I_ESP_Kx2y4z_D2z_aa+ABX*I_ESP_I2y4z_D2z_aa;
    Double I_ESP_Iy5z_Fx2z_aa = I_ESP_Kxy5z_D2z_aa+ABX*I_ESP_Iy5z_D2z_aa;
    Double I_ESP_I6z_Fx2z_aa = I_ESP_Kx6z_D2z_aa+ABX*I_ESP_I6z_D2z_aa;
    Double I_ESP_I6x_F3y_aa = I_ESP_K6xy_D2y_aa+ABY*I_ESP_I6x_D2y_aa;
    Double I_ESP_I5xy_F3y_aa = I_ESP_K5x2y_D2y_aa+ABY*I_ESP_I5xy_D2y_aa;
    Double I_ESP_I5xz_F3y_aa = I_ESP_K5xyz_D2y_aa+ABY*I_ESP_I5xz_D2y_aa;
    Double I_ESP_I4x2y_F3y_aa = I_ESP_K4x3y_D2y_aa+ABY*I_ESP_I4x2y_D2y_aa;
    Double I_ESP_I4xyz_F3y_aa = I_ESP_K4x2yz_D2y_aa+ABY*I_ESP_I4xyz_D2y_aa;
    Double I_ESP_I4x2z_F3y_aa = I_ESP_K4xy2z_D2y_aa+ABY*I_ESP_I4x2z_D2y_aa;
    Double I_ESP_I3x3y_F3y_aa = I_ESP_K3x4y_D2y_aa+ABY*I_ESP_I3x3y_D2y_aa;
    Double I_ESP_I3x2yz_F3y_aa = I_ESP_K3x3yz_D2y_aa+ABY*I_ESP_I3x2yz_D2y_aa;
    Double I_ESP_I3xy2z_F3y_aa = I_ESP_K3x2y2z_D2y_aa+ABY*I_ESP_I3xy2z_D2y_aa;
    Double I_ESP_I3x3z_F3y_aa = I_ESP_K3xy3z_D2y_aa+ABY*I_ESP_I3x3z_D2y_aa;
    Double I_ESP_I2x4y_F3y_aa = I_ESP_K2x5y_D2y_aa+ABY*I_ESP_I2x4y_D2y_aa;
    Double I_ESP_I2x3yz_F3y_aa = I_ESP_K2x4yz_D2y_aa+ABY*I_ESP_I2x3yz_D2y_aa;
    Double I_ESP_I2x2y2z_F3y_aa = I_ESP_K2x3y2z_D2y_aa+ABY*I_ESP_I2x2y2z_D2y_aa;
    Double I_ESP_I2xy3z_F3y_aa = I_ESP_K2x2y3z_D2y_aa+ABY*I_ESP_I2xy3z_D2y_aa;
    Double I_ESP_I2x4z_F3y_aa = I_ESP_K2xy4z_D2y_aa+ABY*I_ESP_I2x4z_D2y_aa;
    Double I_ESP_Ix5y_F3y_aa = I_ESP_Kx6y_D2y_aa+ABY*I_ESP_Ix5y_D2y_aa;
    Double I_ESP_Ix4yz_F3y_aa = I_ESP_Kx5yz_D2y_aa+ABY*I_ESP_Ix4yz_D2y_aa;
    Double I_ESP_Ix3y2z_F3y_aa = I_ESP_Kx4y2z_D2y_aa+ABY*I_ESP_Ix3y2z_D2y_aa;
    Double I_ESP_Ix2y3z_F3y_aa = I_ESP_Kx3y3z_D2y_aa+ABY*I_ESP_Ix2y3z_D2y_aa;
    Double I_ESP_Ixy4z_F3y_aa = I_ESP_Kx2y4z_D2y_aa+ABY*I_ESP_Ixy4z_D2y_aa;
    Double I_ESP_Ix5z_F3y_aa = I_ESP_Kxy5z_D2y_aa+ABY*I_ESP_Ix5z_D2y_aa;
    Double I_ESP_I6y_F3y_aa = I_ESP_K7y_D2y_aa+ABY*I_ESP_I6y_D2y_aa;
    Double I_ESP_I5yz_F3y_aa = I_ESP_K6yz_D2y_aa+ABY*I_ESP_I5yz_D2y_aa;
    Double I_ESP_I4y2z_F3y_aa = I_ESP_K5y2z_D2y_aa+ABY*I_ESP_I4y2z_D2y_aa;
    Double I_ESP_I3y3z_F3y_aa = I_ESP_K4y3z_D2y_aa+ABY*I_ESP_I3y3z_D2y_aa;
    Double I_ESP_I2y4z_F3y_aa = I_ESP_K3y4z_D2y_aa+ABY*I_ESP_I2y4z_D2y_aa;
    Double I_ESP_Iy5z_F3y_aa = I_ESP_K2y5z_D2y_aa+ABY*I_ESP_Iy5z_D2y_aa;
    Double I_ESP_I6z_F3y_aa = I_ESP_Ky6z_D2y_aa+ABY*I_ESP_I6z_D2y_aa;
    Double I_ESP_I6x_F2yz_aa = I_ESP_K6xz_D2y_aa+ABZ*I_ESP_I6x_D2y_aa;
    Double I_ESP_I5xy_F2yz_aa = I_ESP_K5xyz_D2y_aa+ABZ*I_ESP_I5xy_D2y_aa;
    Double I_ESP_I5xz_F2yz_aa = I_ESP_K5x2z_D2y_aa+ABZ*I_ESP_I5xz_D2y_aa;
    Double I_ESP_I4x2y_F2yz_aa = I_ESP_K4x2yz_D2y_aa+ABZ*I_ESP_I4x2y_D2y_aa;
    Double I_ESP_I4xyz_F2yz_aa = I_ESP_K4xy2z_D2y_aa+ABZ*I_ESP_I4xyz_D2y_aa;
    Double I_ESP_I4x2z_F2yz_aa = I_ESP_K4x3z_D2y_aa+ABZ*I_ESP_I4x2z_D2y_aa;
    Double I_ESP_I3x3y_F2yz_aa = I_ESP_K3x3yz_D2y_aa+ABZ*I_ESP_I3x3y_D2y_aa;
    Double I_ESP_I3x2yz_F2yz_aa = I_ESP_K3x2y2z_D2y_aa+ABZ*I_ESP_I3x2yz_D2y_aa;
    Double I_ESP_I3xy2z_F2yz_aa = I_ESP_K3xy3z_D2y_aa+ABZ*I_ESP_I3xy2z_D2y_aa;
    Double I_ESP_I3x3z_F2yz_aa = I_ESP_K3x4z_D2y_aa+ABZ*I_ESP_I3x3z_D2y_aa;
    Double I_ESP_I2x4y_F2yz_aa = I_ESP_K2x4yz_D2y_aa+ABZ*I_ESP_I2x4y_D2y_aa;
    Double I_ESP_I2x3yz_F2yz_aa = I_ESP_K2x3y2z_D2y_aa+ABZ*I_ESP_I2x3yz_D2y_aa;
    Double I_ESP_I2x2y2z_F2yz_aa = I_ESP_K2x2y3z_D2y_aa+ABZ*I_ESP_I2x2y2z_D2y_aa;
    Double I_ESP_I2xy3z_F2yz_aa = I_ESP_K2xy4z_D2y_aa+ABZ*I_ESP_I2xy3z_D2y_aa;
    Double I_ESP_I2x4z_F2yz_aa = I_ESP_K2x5z_D2y_aa+ABZ*I_ESP_I2x4z_D2y_aa;
    Double I_ESP_Ix5y_F2yz_aa = I_ESP_Kx5yz_D2y_aa+ABZ*I_ESP_Ix5y_D2y_aa;
    Double I_ESP_Ix4yz_F2yz_aa = I_ESP_Kx4y2z_D2y_aa+ABZ*I_ESP_Ix4yz_D2y_aa;
    Double I_ESP_Ix3y2z_F2yz_aa = I_ESP_Kx3y3z_D2y_aa+ABZ*I_ESP_Ix3y2z_D2y_aa;
    Double I_ESP_Ix2y3z_F2yz_aa = I_ESP_Kx2y4z_D2y_aa+ABZ*I_ESP_Ix2y3z_D2y_aa;
    Double I_ESP_Ixy4z_F2yz_aa = I_ESP_Kxy5z_D2y_aa+ABZ*I_ESP_Ixy4z_D2y_aa;
    Double I_ESP_Ix5z_F2yz_aa = I_ESP_Kx6z_D2y_aa+ABZ*I_ESP_Ix5z_D2y_aa;
    Double I_ESP_I6y_F2yz_aa = I_ESP_K6yz_D2y_aa+ABZ*I_ESP_I6y_D2y_aa;
    Double I_ESP_I5yz_F2yz_aa = I_ESP_K5y2z_D2y_aa+ABZ*I_ESP_I5yz_D2y_aa;
    Double I_ESP_I4y2z_F2yz_aa = I_ESP_K4y3z_D2y_aa+ABZ*I_ESP_I4y2z_D2y_aa;
    Double I_ESP_I3y3z_F2yz_aa = I_ESP_K3y4z_D2y_aa+ABZ*I_ESP_I3y3z_D2y_aa;
    Double I_ESP_I2y4z_F2yz_aa = I_ESP_K2y5z_D2y_aa+ABZ*I_ESP_I2y4z_D2y_aa;
    Double I_ESP_Iy5z_F2yz_aa = I_ESP_Ky6z_D2y_aa+ABZ*I_ESP_Iy5z_D2y_aa;
    Double I_ESP_I6z_F2yz_aa = I_ESP_K7z_D2y_aa+ABZ*I_ESP_I6z_D2y_aa;
    Double I_ESP_I6x_Fy2z_aa = I_ESP_K6xy_D2z_aa+ABY*I_ESP_I6x_D2z_aa;
    Double I_ESP_I5xy_Fy2z_aa = I_ESP_K5x2y_D2z_aa+ABY*I_ESP_I5xy_D2z_aa;
    Double I_ESP_I5xz_Fy2z_aa = I_ESP_K5xyz_D2z_aa+ABY*I_ESP_I5xz_D2z_aa;
    Double I_ESP_I4x2y_Fy2z_aa = I_ESP_K4x3y_D2z_aa+ABY*I_ESP_I4x2y_D2z_aa;
    Double I_ESP_I4xyz_Fy2z_aa = I_ESP_K4x2yz_D2z_aa+ABY*I_ESP_I4xyz_D2z_aa;
    Double I_ESP_I4x2z_Fy2z_aa = I_ESP_K4xy2z_D2z_aa+ABY*I_ESP_I4x2z_D2z_aa;
    Double I_ESP_I3x3y_Fy2z_aa = I_ESP_K3x4y_D2z_aa+ABY*I_ESP_I3x3y_D2z_aa;
    Double I_ESP_I3x2yz_Fy2z_aa = I_ESP_K3x3yz_D2z_aa+ABY*I_ESP_I3x2yz_D2z_aa;
    Double I_ESP_I3xy2z_Fy2z_aa = I_ESP_K3x2y2z_D2z_aa+ABY*I_ESP_I3xy2z_D2z_aa;
    Double I_ESP_I3x3z_Fy2z_aa = I_ESP_K3xy3z_D2z_aa+ABY*I_ESP_I3x3z_D2z_aa;
    Double I_ESP_I2x4y_Fy2z_aa = I_ESP_K2x5y_D2z_aa+ABY*I_ESP_I2x4y_D2z_aa;
    Double I_ESP_I2x3yz_Fy2z_aa = I_ESP_K2x4yz_D2z_aa+ABY*I_ESP_I2x3yz_D2z_aa;
    Double I_ESP_I2x2y2z_Fy2z_aa = I_ESP_K2x3y2z_D2z_aa+ABY*I_ESP_I2x2y2z_D2z_aa;
    Double I_ESP_I2xy3z_Fy2z_aa = I_ESP_K2x2y3z_D2z_aa+ABY*I_ESP_I2xy3z_D2z_aa;
    Double I_ESP_I2x4z_Fy2z_aa = I_ESP_K2xy4z_D2z_aa+ABY*I_ESP_I2x4z_D2z_aa;
    Double I_ESP_Ix5y_Fy2z_aa = I_ESP_Kx6y_D2z_aa+ABY*I_ESP_Ix5y_D2z_aa;
    Double I_ESP_Ix4yz_Fy2z_aa = I_ESP_Kx5yz_D2z_aa+ABY*I_ESP_Ix4yz_D2z_aa;
    Double I_ESP_Ix3y2z_Fy2z_aa = I_ESP_Kx4y2z_D2z_aa+ABY*I_ESP_Ix3y2z_D2z_aa;
    Double I_ESP_Ix2y3z_Fy2z_aa = I_ESP_Kx3y3z_D2z_aa+ABY*I_ESP_Ix2y3z_D2z_aa;
    Double I_ESP_Ixy4z_Fy2z_aa = I_ESP_Kx2y4z_D2z_aa+ABY*I_ESP_Ixy4z_D2z_aa;
    Double I_ESP_Ix5z_Fy2z_aa = I_ESP_Kxy5z_D2z_aa+ABY*I_ESP_Ix5z_D2z_aa;
    Double I_ESP_I6y_Fy2z_aa = I_ESP_K7y_D2z_aa+ABY*I_ESP_I6y_D2z_aa;
    Double I_ESP_I5yz_Fy2z_aa = I_ESP_K6yz_D2z_aa+ABY*I_ESP_I5yz_D2z_aa;
    Double I_ESP_I4y2z_Fy2z_aa = I_ESP_K5y2z_D2z_aa+ABY*I_ESP_I4y2z_D2z_aa;
    Double I_ESP_I3y3z_Fy2z_aa = I_ESP_K4y3z_D2z_aa+ABY*I_ESP_I3y3z_D2z_aa;
    Double I_ESP_I2y4z_Fy2z_aa = I_ESP_K3y4z_D2z_aa+ABY*I_ESP_I2y4z_D2z_aa;
    Double I_ESP_Iy5z_Fy2z_aa = I_ESP_K2y5z_D2z_aa+ABY*I_ESP_Iy5z_D2z_aa;
    Double I_ESP_I6z_Fy2z_aa = I_ESP_Ky6z_D2z_aa+ABY*I_ESP_I6z_D2z_aa;
    Double I_ESP_I6x_F3z_aa = I_ESP_K6xz_D2z_aa+ABZ*I_ESP_I6x_D2z_aa;
    Double I_ESP_I5xy_F3z_aa = I_ESP_K5xyz_D2z_aa+ABZ*I_ESP_I5xy_D2z_aa;
    Double I_ESP_I5xz_F3z_aa = I_ESP_K5x2z_D2z_aa+ABZ*I_ESP_I5xz_D2z_aa;
    Double I_ESP_I4x2y_F3z_aa = I_ESP_K4x2yz_D2z_aa+ABZ*I_ESP_I4x2y_D2z_aa;
    Double I_ESP_I4xyz_F3z_aa = I_ESP_K4xy2z_D2z_aa+ABZ*I_ESP_I4xyz_D2z_aa;
    Double I_ESP_I4x2z_F3z_aa = I_ESP_K4x3z_D2z_aa+ABZ*I_ESP_I4x2z_D2z_aa;
    Double I_ESP_I3x3y_F3z_aa = I_ESP_K3x3yz_D2z_aa+ABZ*I_ESP_I3x3y_D2z_aa;
    Double I_ESP_I3x2yz_F3z_aa = I_ESP_K3x2y2z_D2z_aa+ABZ*I_ESP_I3x2yz_D2z_aa;
    Double I_ESP_I3xy2z_F3z_aa = I_ESP_K3xy3z_D2z_aa+ABZ*I_ESP_I3xy2z_D2z_aa;
    Double I_ESP_I3x3z_F3z_aa = I_ESP_K3x4z_D2z_aa+ABZ*I_ESP_I3x3z_D2z_aa;
    Double I_ESP_I2x4y_F3z_aa = I_ESP_K2x4yz_D2z_aa+ABZ*I_ESP_I2x4y_D2z_aa;
    Double I_ESP_I2x3yz_F3z_aa = I_ESP_K2x3y2z_D2z_aa+ABZ*I_ESP_I2x3yz_D2z_aa;
    Double I_ESP_I2x2y2z_F3z_aa = I_ESP_K2x2y3z_D2z_aa+ABZ*I_ESP_I2x2y2z_D2z_aa;
    Double I_ESP_I2xy3z_F3z_aa = I_ESP_K2xy4z_D2z_aa+ABZ*I_ESP_I2xy3z_D2z_aa;
    Double I_ESP_I2x4z_F3z_aa = I_ESP_K2x5z_D2z_aa+ABZ*I_ESP_I2x4z_D2z_aa;
    Double I_ESP_Ix5y_F3z_aa = I_ESP_Kx5yz_D2z_aa+ABZ*I_ESP_Ix5y_D2z_aa;
    Double I_ESP_Ix4yz_F3z_aa = I_ESP_Kx4y2z_D2z_aa+ABZ*I_ESP_Ix4yz_D2z_aa;
    Double I_ESP_Ix3y2z_F3z_aa = I_ESP_Kx3y3z_D2z_aa+ABZ*I_ESP_Ix3y2z_D2z_aa;
    Double I_ESP_Ix2y3z_F3z_aa = I_ESP_Kx2y4z_D2z_aa+ABZ*I_ESP_Ix2y3z_D2z_aa;
    Double I_ESP_Ixy4z_F3z_aa = I_ESP_Kxy5z_D2z_aa+ABZ*I_ESP_Ixy4z_D2z_aa;
    Double I_ESP_Ix5z_F3z_aa = I_ESP_Kx6z_D2z_aa+ABZ*I_ESP_Ix5z_D2z_aa;
    Double I_ESP_I6y_F3z_aa = I_ESP_K6yz_D2z_aa+ABZ*I_ESP_I6y_D2z_aa;
    Double I_ESP_I5yz_F3z_aa = I_ESP_K5y2z_D2z_aa+ABZ*I_ESP_I5yz_D2z_aa;
    Double I_ESP_I4y2z_F3z_aa = I_ESP_K4y3z_D2z_aa+ABZ*I_ESP_I4y2z_D2z_aa;
    Double I_ESP_I3y3z_F3z_aa = I_ESP_K3y4z_D2z_aa+ABZ*I_ESP_I3y3z_D2z_aa;
    Double I_ESP_I2y4z_F3z_aa = I_ESP_K2y5z_D2z_aa+ABZ*I_ESP_I2y4z_D2z_aa;
    Double I_ESP_Iy5z_F3z_aa = I_ESP_Ky6z_D2z_aa+ABZ*I_ESP_Iy5z_D2z_aa;
    Double I_ESP_I6z_F3z_aa = I_ESP_K7z_D2z_aa+ABZ*I_ESP_I6z_D2z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_aa
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_D_F
     ************************************************************/
    abcd[iGrid*900+0] = 4.0E0*I_ESP_I6x_F3x_aa-2.0E0*4*I_ESP_G4x_F3x_a-2.0E0*5*I_ESP_G4x_F3x_a+4*3*I_ESP_D2x_F3x;
    abcd[iGrid*900+1] = 4.0E0*I_ESP_I5xy_F3x_aa-2.0E0*3*I_ESP_G3xy_F3x_a-2.0E0*4*I_ESP_G3xy_F3x_a+3*2*I_ESP_Dxy_F3x;
    abcd[iGrid*900+2] = 4.0E0*I_ESP_I5xz_F3x_aa-2.0E0*3*I_ESP_G3xz_F3x_a-2.0E0*4*I_ESP_G3xz_F3x_a+3*2*I_ESP_Dxz_F3x;
    abcd[iGrid*900+3] = 4.0E0*I_ESP_I4x2y_F3x_aa-2.0E0*2*I_ESP_G2x2y_F3x_a-2.0E0*3*I_ESP_G2x2y_F3x_a+2*1*I_ESP_D2y_F3x;
    abcd[iGrid*900+4] = 4.0E0*I_ESP_I4xyz_F3x_aa-2.0E0*2*I_ESP_G2xyz_F3x_a-2.0E0*3*I_ESP_G2xyz_F3x_a+2*1*I_ESP_Dyz_F3x;
    abcd[iGrid*900+5] = 4.0E0*I_ESP_I4x2z_F3x_aa-2.0E0*2*I_ESP_G2x2z_F3x_a-2.0E0*3*I_ESP_G2x2z_F3x_a+2*1*I_ESP_D2z_F3x;
    abcd[iGrid*900+6] = 4.0E0*I_ESP_I3x3y_F3x_aa-2.0E0*1*I_ESP_Gx3y_F3x_a-2.0E0*2*I_ESP_Gx3y_F3x_a;
    abcd[iGrid*900+7] = 4.0E0*I_ESP_I3x2yz_F3x_aa-2.0E0*1*I_ESP_Gx2yz_F3x_a-2.0E0*2*I_ESP_Gx2yz_F3x_a;
    abcd[iGrid*900+8] = 4.0E0*I_ESP_I3xy2z_F3x_aa-2.0E0*1*I_ESP_Gxy2z_F3x_a-2.0E0*2*I_ESP_Gxy2z_F3x_a;
    abcd[iGrid*900+9] = 4.0E0*I_ESP_I3x3z_F3x_aa-2.0E0*1*I_ESP_Gx3z_F3x_a-2.0E0*2*I_ESP_Gx3z_F3x_a;
    abcd[iGrid*900+10] = 4.0E0*I_ESP_I2x4y_F3x_aa-2.0E0*1*I_ESP_G4y_F3x_a;
    abcd[iGrid*900+11] = 4.0E0*I_ESP_I2x3yz_F3x_aa-2.0E0*1*I_ESP_G3yz_F3x_a;
    abcd[iGrid*900+12] = 4.0E0*I_ESP_I2x2y2z_F3x_aa-2.0E0*1*I_ESP_G2y2z_F3x_a;
    abcd[iGrid*900+13] = 4.0E0*I_ESP_I2xy3z_F3x_aa-2.0E0*1*I_ESP_Gy3z_F3x_a;
    abcd[iGrid*900+14] = 4.0E0*I_ESP_I2x4z_F3x_aa-2.0E0*1*I_ESP_G4z_F3x_a;
    abcd[iGrid*900+15] = 4.0E0*I_ESP_I6x_F2xy_aa-2.0E0*4*I_ESP_G4x_F2xy_a-2.0E0*5*I_ESP_G4x_F2xy_a+4*3*I_ESP_D2x_F2xy;
    abcd[iGrid*900+16] = 4.0E0*I_ESP_I5xy_F2xy_aa-2.0E0*3*I_ESP_G3xy_F2xy_a-2.0E0*4*I_ESP_G3xy_F2xy_a+3*2*I_ESP_Dxy_F2xy;
    abcd[iGrid*900+17] = 4.0E0*I_ESP_I5xz_F2xy_aa-2.0E0*3*I_ESP_G3xz_F2xy_a-2.0E0*4*I_ESP_G3xz_F2xy_a+3*2*I_ESP_Dxz_F2xy;
    abcd[iGrid*900+18] = 4.0E0*I_ESP_I4x2y_F2xy_aa-2.0E0*2*I_ESP_G2x2y_F2xy_a-2.0E0*3*I_ESP_G2x2y_F2xy_a+2*1*I_ESP_D2y_F2xy;
    abcd[iGrid*900+19] = 4.0E0*I_ESP_I4xyz_F2xy_aa-2.0E0*2*I_ESP_G2xyz_F2xy_a-2.0E0*3*I_ESP_G2xyz_F2xy_a+2*1*I_ESP_Dyz_F2xy;
    abcd[iGrid*900+20] = 4.0E0*I_ESP_I4x2z_F2xy_aa-2.0E0*2*I_ESP_G2x2z_F2xy_a-2.0E0*3*I_ESP_G2x2z_F2xy_a+2*1*I_ESP_D2z_F2xy;
    abcd[iGrid*900+21] = 4.0E0*I_ESP_I3x3y_F2xy_aa-2.0E0*1*I_ESP_Gx3y_F2xy_a-2.0E0*2*I_ESP_Gx3y_F2xy_a;
    abcd[iGrid*900+22] = 4.0E0*I_ESP_I3x2yz_F2xy_aa-2.0E0*1*I_ESP_Gx2yz_F2xy_a-2.0E0*2*I_ESP_Gx2yz_F2xy_a;
    abcd[iGrid*900+23] = 4.0E0*I_ESP_I3xy2z_F2xy_aa-2.0E0*1*I_ESP_Gxy2z_F2xy_a-2.0E0*2*I_ESP_Gxy2z_F2xy_a;
    abcd[iGrid*900+24] = 4.0E0*I_ESP_I3x3z_F2xy_aa-2.0E0*1*I_ESP_Gx3z_F2xy_a-2.0E0*2*I_ESP_Gx3z_F2xy_a;
    abcd[iGrid*900+25] = 4.0E0*I_ESP_I2x4y_F2xy_aa-2.0E0*1*I_ESP_G4y_F2xy_a;
    abcd[iGrid*900+26] = 4.0E0*I_ESP_I2x3yz_F2xy_aa-2.0E0*1*I_ESP_G3yz_F2xy_a;
    abcd[iGrid*900+27] = 4.0E0*I_ESP_I2x2y2z_F2xy_aa-2.0E0*1*I_ESP_G2y2z_F2xy_a;
    abcd[iGrid*900+28] = 4.0E0*I_ESP_I2xy3z_F2xy_aa-2.0E0*1*I_ESP_Gy3z_F2xy_a;
    abcd[iGrid*900+29] = 4.0E0*I_ESP_I2x4z_F2xy_aa-2.0E0*1*I_ESP_G4z_F2xy_a;
    abcd[iGrid*900+30] = 4.0E0*I_ESP_I6x_F2xz_aa-2.0E0*4*I_ESP_G4x_F2xz_a-2.0E0*5*I_ESP_G4x_F2xz_a+4*3*I_ESP_D2x_F2xz;
    abcd[iGrid*900+31] = 4.0E0*I_ESP_I5xy_F2xz_aa-2.0E0*3*I_ESP_G3xy_F2xz_a-2.0E0*4*I_ESP_G3xy_F2xz_a+3*2*I_ESP_Dxy_F2xz;
    abcd[iGrid*900+32] = 4.0E0*I_ESP_I5xz_F2xz_aa-2.0E0*3*I_ESP_G3xz_F2xz_a-2.0E0*4*I_ESP_G3xz_F2xz_a+3*2*I_ESP_Dxz_F2xz;
    abcd[iGrid*900+33] = 4.0E0*I_ESP_I4x2y_F2xz_aa-2.0E0*2*I_ESP_G2x2y_F2xz_a-2.0E0*3*I_ESP_G2x2y_F2xz_a+2*1*I_ESP_D2y_F2xz;
    abcd[iGrid*900+34] = 4.0E0*I_ESP_I4xyz_F2xz_aa-2.0E0*2*I_ESP_G2xyz_F2xz_a-2.0E0*3*I_ESP_G2xyz_F2xz_a+2*1*I_ESP_Dyz_F2xz;
    abcd[iGrid*900+35] = 4.0E0*I_ESP_I4x2z_F2xz_aa-2.0E0*2*I_ESP_G2x2z_F2xz_a-2.0E0*3*I_ESP_G2x2z_F2xz_a+2*1*I_ESP_D2z_F2xz;
    abcd[iGrid*900+36] = 4.0E0*I_ESP_I3x3y_F2xz_aa-2.0E0*1*I_ESP_Gx3y_F2xz_a-2.0E0*2*I_ESP_Gx3y_F2xz_a;
    abcd[iGrid*900+37] = 4.0E0*I_ESP_I3x2yz_F2xz_aa-2.0E0*1*I_ESP_Gx2yz_F2xz_a-2.0E0*2*I_ESP_Gx2yz_F2xz_a;
    abcd[iGrid*900+38] = 4.0E0*I_ESP_I3xy2z_F2xz_aa-2.0E0*1*I_ESP_Gxy2z_F2xz_a-2.0E0*2*I_ESP_Gxy2z_F2xz_a;
    abcd[iGrid*900+39] = 4.0E0*I_ESP_I3x3z_F2xz_aa-2.0E0*1*I_ESP_Gx3z_F2xz_a-2.0E0*2*I_ESP_Gx3z_F2xz_a;
    abcd[iGrid*900+40] = 4.0E0*I_ESP_I2x4y_F2xz_aa-2.0E0*1*I_ESP_G4y_F2xz_a;
    abcd[iGrid*900+41] = 4.0E0*I_ESP_I2x3yz_F2xz_aa-2.0E0*1*I_ESP_G3yz_F2xz_a;
    abcd[iGrid*900+42] = 4.0E0*I_ESP_I2x2y2z_F2xz_aa-2.0E0*1*I_ESP_G2y2z_F2xz_a;
    abcd[iGrid*900+43] = 4.0E0*I_ESP_I2xy3z_F2xz_aa-2.0E0*1*I_ESP_Gy3z_F2xz_a;
    abcd[iGrid*900+44] = 4.0E0*I_ESP_I2x4z_F2xz_aa-2.0E0*1*I_ESP_G4z_F2xz_a;
    abcd[iGrid*900+45] = 4.0E0*I_ESP_I6x_Fx2y_aa-2.0E0*4*I_ESP_G4x_Fx2y_a-2.0E0*5*I_ESP_G4x_Fx2y_a+4*3*I_ESP_D2x_Fx2y;
    abcd[iGrid*900+46] = 4.0E0*I_ESP_I5xy_Fx2y_aa-2.0E0*3*I_ESP_G3xy_Fx2y_a-2.0E0*4*I_ESP_G3xy_Fx2y_a+3*2*I_ESP_Dxy_Fx2y;
    abcd[iGrid*900+47] = 4.0E0*I_ESP_I5xz_Fx2y_aa-2.0E0*3*I_ESP_G3xz_Fx2y_a-2.0E0*4*I_ESP_G3xz_Fx2y_a+3*2*I_ESP_Dxz_Fx2y;
    abcd[iGrid*900+48] = 4.0E0*I_ESP_I4x2y_Fx2y_aa-2.0E0*2*I_ESP_G2x2y_Fx2y_a-2.0E0*3*I_ESP_G2x2y_Fx2y_a+2*1*I_ESP_D2y_Fx2y;
    abcd[iGrid*900+49] = 4.0E0*I_ESP_I4xyz_Fx2y_aa-2.0E0*2*I_ESP_G2xyz_Fx2y_a-2.0E0*3*I_ESP_G2xyz_Fx2y_a+2*1*I_ESP_Dyz_Fx2y;
    abcd[iGrid*900+50] = 4.0E0*I_ESP_I4x2z_Fx2y_aa-2.0E0*2*I_ESP_G2x2z_Fx2y_a-2.0E0*3*I_ESP_G2x2z_Fx2y_a+2*1*I_ESP_D2z_Fx2y;
    abcd[iGrid*900+51] = 4.0E0*I_ESP_I3x3y_Fx2y_aa-2.0E0*1*I_ESP_Gx3y_Fx2y_a-2.0E0*2*I_ESP_Gx3y_Fx2y_a;
    abcd[iGrid*900+52] = 4.0E0*I_ESP_I3x2yz_Fx2y_aa-2.0E0*1*I_ESP_Gx2yz_Fx2y_a-2.0E0*2*I_ESP_Gx2yz_Fx2y_a;
    abcd[iGrid*900+53] = 4.0E0*I_ESP_I3xy2z_Fx2y_aa-2.0E0*1*I_ESP_Gxy2z_Fx2y_a-2.0E0*2*I_ESP_Gxy2z_Fx2y_a;
    abcd[iGrid*900+54] = 4.0E0*I_ESP_I3x3z_Fx2y_aa-2.0E0*1*I_ESP_Gx3z_Fx2y_a-2.0E0*2*I_ESP_Gx3z_Fx2y_a;
    abcd[iGrid*900+55] = 4.0E0*I_ESP_I2x4y_Fx2y_aa-2.0E0*1*I_ESP_G4y_Fx2y_a;
    abcd[iGrid*900+56] = 4.0E0*I_ESP_I2x3yz_Fx2y_aa-2.0E0*1*I_ESP_G3yz_Fx2y_a;
    abcd[iGrid*900+57] = 4.0E0*I_ESP_I2x2y2z_Fx2y_aa-2.0E0*1*I_ESP_G2y2z_Fx2y_a;
    abcd[iGrid*900+58] = 4.0E0*I_ESP_I2xy3z_Fx2y_aa-2.0E0*1*I_ESP_Gy3z_Fx2y_a;
    abcd[iGrid*900+59] = 4.0E0*I_ESP_I2x4z_Fx2y_aa-2.0E0*1*I_ESP_G4z_Fx2y_a;
    abcd[iGrid*900+60] = 4.0E0*I_ESP_I6x_Fxyz_aa-2.0E0*4*I_ESP_G4x_Fxyz_a-2.0E0*5*I_ESP_G4x_Fxyz_a+4*3*I_ESP_D2x_Fxyz;
    abcd[iGrid*900+61] = 4.0E0*I_ESP_I5xy_Fxyz_aa-2.0E0*3*I_ESP_G3xy_Fxyz_a-2.0E0*4*I_ESP_G3xy_Fxyz_a+3*2*I_ESP_Dxy_Fxyz;
    abcd[iGrid*900+62] = 4.0E0*I_ESP_I5xz_Fxyz_aa-2.0E0*3*I_ESP_G3xz_Fxyz_a-2.0E0*4*I_ESP_G3xz_Fxyz_a+3*2*I_ESP_Dxz_Fxyz;
    abcd[iGrid*900+63] = 4.0E0*I_ESP_I4x2y_Fxyz_aa-2.0E0*2*I_ESP_G2x2y_Fxyz_a-2.0E0*3*I_ESP_G2x2y_Fxyz_a+2*1*I_ESP_D2y_Fxyz;
    abcd[iGrid*900+64] = 4.0E0*I_ESP_I4xyz_Fxyz_aa-2.0E0*2*I_ESP_G2xyz_Fxyz_a-2.0E0*3*I_ESP_G2xyz_Fxyz_a+2*1*I_ESP_Dyz_Fxyz;
    abcd[iGrid*900+65] = 4.0E0*I_ESP_I4x2z_Fxyz_aa-2.0E0*2*I_ESP_G2x2z_Fxyz_a-2.0E0*3*I_ESP_G2x2z_Fxyz_a+2*1*I_ESP_D2z_Fxyz;
    abcd[iGrid*900+66] = 4.0E0*I_ESP_I3x3y_Fxyz_aa-2.0E0*1*I_ESP_Gx3y_Fxyz_a-2.0E0*2*I_ESP_Gx3y_Fxyz_a;
    abcd[iGrid*900+67] = 4.0E0*I_ESP_I3x2yz_Fxyz_aa-2.0E0*1*I_ESP_Gx2yz_Fxyz_a-2.0E0*2*I_ESP_Gx2yz_Fxyz_a;
    abcd[iGrid*900+68] = 4.0E0*I_ESP_I3xy2z_Fxyz_aa-2.0E0*1*I_ESP_Gxy2z_Fxyz_a-2.0E0*2*I_ESP_Gxy2z_Fxyz_a;
    abcd[iGrid*900+69] = 4.0E0*I_ESP_I3x3z_Fxyz_aa-2.0E0*1*I_ESP_Gx3z_Fxyz_a-2.0E0*2*I_ESP_Gx3z_Fxyz_a;
    abcd[iGrid*900+70] = 4.0E0*I_ESP_I2x4y_Fxyz_aa-2.0E0*1*I_ESP_G4y_Fxyz_a;
    abcd[iGrid*900+71] = 4.0E0*I_ESP_I2x3yz_Fxyz_aa-2.0E0*1*I_ESP_G3yz_Fxyz_a;
    abcd[iGrid*900+72] = 4.0E0*I_ESP_I2x2y2z_Fxyz_aa-2.0E0*1*I_ESP_G2y2z_Fxyz_a;
    abcd[iGrid*900+73] = 4.0E0*I_ESP_I2xy3z_Fxyz_aa-2.0E0*1*I_ESP_Gy3z_Fxyz_a;
    abcd[iGrid*900+74] = 4.0E0*I_ESP_I2x4z_Fxyz_aa-2.0E0*1*I_ESP_G4z_Fxyz_a;
    abcd[iGrid*900+75] = 4.0E0*I_ESP_I6x_Fx2z_aa-2.0E0*4*I_ESP_G4x_Fx2z_a-2.0E0*5*I_ESP_G4x_Fx2z_a+4*3*I_ESP_D2x_Fx2z;
    abcd[iGrid*900+76] = 4.0E0*I_ESP_I5xy_Fx2z_aa-2.0E0*3*I_ESP_G3xy_Fx2z_a-2.0E0*4*I_ESP_G3xy_Fx2z_a+3*2*I_ESP_Dxy_Fx2z;
    abcd[iGrid*900+77] = 4.0E0*I_ESP_I5xz_Fx2z_aa-2.0E0*3*I_ESP_G3xz_Fx2z_a-2.0E0*4*I_ESP_G3xz_Fx2z_a+3*2*I_ESP_Dxz_Fx2z;
    abcd[iGrid*900+78] = 4.0E0*I_ESP_I4x2y_Fx2z_aa-2.0E0*2*I_ESP_G2x2y_Fx2z_a-2.0E0*3*I_ESP_G2x2y_Fx2z_a+2*1*I_ESP_D2y_Fx2z;
    abcd[iGrid*900+79] = 4.0E0*I_ESP_I4xyz_Fx2z_aa-2.0E0*2*I_ESP_G2xyz_Fx2z_a-2.0E0*3*I_ESP_G2xyz_Fx2z_a+2*1*I_ESP_Dyz_Fx2z;
    abcd[iGrid*900+80] = 4.0E0*I_ESP_I4x2z_Fx2z_aa-2.0E0*2*I_ESP_G2x2z_Fx2z_a-2.0E0*3*I_ESP_G2x2z_Fx2z_a+2*1*I_ESP_D2z_Fx2z;
    abcd[iGrid*900+81] = 4.0E0*I_ESP_I3x3y_Fx2z_aa-2.0E0*1*I_ESP_Gx3y_Fx2z_a-2.0E0*2*I_ESP_Gx3y_Fx2z_a;
    abcd[iGrid*900+82] = 4.0E0*I_ESP_I3x2yz_Fx2z_aa-2.0E0*1*I_ESP_Gx2yz_Fx2z_a-2.0E0*2*I_ESP_Gx2yz_Fx2z_a;
    abcd[iGrid*900+83] = 4.0E0*I_ESP_I3xy2z_Fx2z_aa-2.0E0*1*I_ESP_Gxy2z_Fx2z_a-2.0E0*2*I_ESP_Gxy2z_Fx2z_a;
    abcd[iGrid*900+84] = 4.0E0*I_ESP_I3x3z_Fx2z_aa-2.0E0*1*I_ESP_Gx3z_Fx2z_a-2.0E0*2*I_ESP_Gx3z_Fx2z_a;
    abcd[iGrid*900+85] = 4.0E0*I_ESP_I2x4y_Fx2z_aa-2.0E0*1*I_ESP_G4y_Fx2z_a;
    abcd[iGrid*900+86] = 4.0E0*I_ESP_I2x3yz_Fx2z_aa-2.0E0*1*I_ESP_G3yz_Fx2z_a;
    abcd[iGrid*900+87] = 4.0E0*I_ESP_I2x2y2z_Fx2z_aa-2.0E0*1*I_ESP_G2y2z_Fx2z_a;
    abcd[iGrid*900+88] = 4.0E0*I_ESP_I2xy3z_Fx2z_aa-2.0E0*1*I_ESP_Gy3z_Fx2z_a;
    abcd[iGrid*900+89] = 4.0E0*I_ESP_I2x4z_Fx2z_aa-2.0E0*1*I_ESP_G4z_Fx2z_a;
    abcd[iGrid*900+90] = 4.0E0*I_ESP_I6x_F3y_aa-2.0E0*4*I_ESP_G4x_F3y_a-2.0E0*5*I_ESP_G4x_F3y_a+4*3*I_ESP_D2x_F3y;
    abcd[iGrid*900+91] = 4.0E0*I_ESP_I5xy_F3y_aa-2.0E0*3*I_ESP_G3xy_F3y_a-2.0E0*4*I_ESP_G3xy_F3y_a+3*2*I_ESP_Dxy_F3y;
    abcd[iGrid*900+92] = 4.0E0*I_ESP_I5xz_F3y_aa-2.0E0*3*I_ESP_G3xz_F3y_a-2.0E0*4*I_ESP_G3xz_F3y_a+3*2*I_ESP_Dxz_F3y;
    abcd[iGrid*900+93] = 4.0E0*I_ESP_I4x2y_F3y_aa-2.0E0*2*I_ESP_G2x2y_F3y_a-2.0E0*3*I_ESP_G2x2y_F3y_a+2*1*I_ESP_D2y_F3y;
    abcd[iGrid*900+94] = 4.0E0*I_ESP_I4xyz_F3y_aa-2.0E0*2*I_ESP_G2xyz_F3y_a-2.0E0*3*I_ESP_G2xyz_F3y_a+2*1*I_ESP_Dyz_F3y;
    abcd[iGrid*900+95] = 4.0E0*I_ESP_I4x2z_F3y_aa-2.0E0*2*I_ESP_G2x2z_F3y_a-2.0E0*3*I_ESP_G2x2z_F3y_a+2*1*I_ESP_D2z_F3y;
    abcd[iGrid*900+96] = 4.0E0*I_ESP_I3x3y_F3y_aa-2.0E0*1*I_ESP_Gx3y_F3y_a-2.0E0*2*I_ESP_Gx3y_F3y_a;
    abcd[iGrid*900+97] = 4.0E0*I_ESP_I3x2yz_F3y_aa-2.0E0*1*I_ESP_Gx2yz_F3y_a-2.0E0*2*I_ESP_Gx2yz_F3y_a;
    abcd[iGrid*900+98] = 4.0E0*I_ESP_I3xy2z_F3y_aa-2.0E0*1*I_ESP_Gxy2z_F3y_a-2.0E0*2*I_ESP_Gxy2z_F3y_a;
    abcd[iGrid*900+99] = 4.0E0*I_ESP_I3x3z_F3y_aa-2.0E0*1*I_ESP_Gx3z_F3y_a-2.0E0*2*I_ESP_Gx3z_F3y_a;
    abcd[iGrid*900+100] = 4.0E0*I_ESP_I2x4y_F3y_aa-2.0E0*1*I_ESP_G4y_F3y_a;
    abcd[iGrid*900+101] = 4.0E0*I_ESP_I2x3yz_F3y_aa-2.0E0*1*I_ESP_G3yz_F3y_a;
    abcd[iGrid*900+102] = 4.0E0*I_ESP_I2x2y2z_F3y_aa-2.0E0*1*I_ESP_G2y2z_F3y_a;
    abcd[iGrid*900+103] = 4.0E0*I_ESP_I2xy3z_F3y_aa-2.0E0*1*I_ESP_Gy3z_F3y_a;
    abcd[iGrid*900+104] = 4.0E0*I_ESP_I2x4z_F3y_aa-2.0E0*1*I_ESP_G4z_F3y_a;
    abcd[iGrid*900+105] = 4.0E0*I_ESP_I6x_F2yz_aa-2.0E0*4*I_ESP_G4x_F2yz_a-2.0E0*5*I_ESP_G4x_F2yz_a+4*3*I_ESP_D2x_F2yz;
    abcd[iGrid*900+106] = 4.0E0*I_ESP_I5xy_F2yz_aa-2.0E0*3*I_ESP_G3xy_F2yz_a-2.0E0*4*I_ESP_G3xy_F2yz_a+3*2*I_ESP_Dxy_F2yz;
    abcd[iGrid*900+107] = 4.0E0*I_ESP_I5xz_F2yz_aa-2.0E0*3*I_ESP_G3xz_F2yz_a-2.0E0*4*I_ESP_G3xz_F2yz_a+3*2*I_ESP_Dxz_F2yz;
    abcd[iGrid*900+108] = 4.0E0*I_ESP_I4x2y_F2yz_aa-2.0E0*2*I_ESP_G2x2y_F2yz_a-2.0E0*3*I_ESP_G2x2y_F2yz_a+2*1*I_ESP_D2y_F2yz;
    abcd[iGrid*900+109] = 4.0E0*I_ESP_I4xyz_F2yz_aa-2.0E0*2*I_ESP_G2xyz_F2yz_a-2.0E0*3*I_ESP_G2xyz_F2yz_a+2*1*I_ESP_Dyz_F2yz;
    abcd[iGrid*900+110] = 4.0E0*I_ESP_I4x2z_F2yz_aa-2.0E0*2*I_ESP_G2x2z_F2yz_a-2.0E0*3*I_ESP_G2x2z_F2yz_a+2*1*I_ESP_D2z_F2yz;
    abcd[iGrid*900+111] = 4.0E0*I_ESP_I3x3y_F2yz_aa-2.0E0*1*I_ESP_Gx3y_F2yz_a-2.0E0*2*I_ESP_Gx3y_F2yz_a;
    abcd[iGrid*900+112] = 4.0E0*I_ESP_I3x2yz_F2yz_aa-2.0E0*1*I_ESP_Gx2yz_F2yz_a-2.0E0*2*I_ESP_Gx2yz_F2yz_a;
    abcd[iGrid*900+113] = 4.0E0*I_ESP_I3xy2z_F2yz_aa-2.0E0*1*I_ESP_Gxy2z_F2yz_a-2.0E0*2*I_ESP_Gxy2z_F2yz_a;
    abcd[iGrid*900+114] = 4.0E0*I_ESP_I3x3z_F2yz_aa-2.0E0*1*I_ESP_Gx3z_F2yz_a-2.0E0*2*I_ESP_Gx3z_F2yz_a;
    abcd[iGrid*900+115] = 4.0E0*I_ESP_I2x4y_F2yz_aa-2.0E0*1*I_ESP_G4y_F2yz_a;
    abcd[iGrid*900+116] = 4.0E0*I_ESP_I2x3yz_F2yz_aa-2.0E0*1*I_ESP_G3yz_F2yz_a;
    abcd[iGrid*900+117] = 4.0E0*I_ESP_I2x2y2z_F2yz_aa-2.0E0*1*I_ESP_G2y2z_F2yz_a;
    abcd[iGrid*900+118] = 4.0E0*I_ESP_I2xy3z_F2yz_aa-2.0E0*1*I_ESP_Gy3z_F2yz_a;
    abcd[iGrid*900+119] = 4.0E0*I_ESP_I2x4z_F2yz_aa-2.0E0*1*I_ESP_G4z_F2yz_a;
    abcd[iGrid*900+120] = 4.0E0*I_ESP_I6x_Fy2z_aa-2.0E0*4*I_ESP_G4x_Fy2z_a-2.0E0*5*I_ESP_G4x_Fy2z_a+4*3*I_ESP_D2x_Fy2z;
    abcd[iGrid*900+121] = 4.0E0*I_ESP_I5xy_Fy2z_aa-2.0E0*3*I_ESP_G3xy_Fy2z_a-2.0E0*4*I_ESP_G3xy_Fy2z_a+3*2*I_ESP_Dxy_Fy2z;
    abcd[iGrid*900+122] = 4.0E0*I_ESP_I5xz_Fy2z_aa-2.0E0*3*I_ESP_G3xz_Fy2z_a-2.0E0*4*I_ESP_G3xz_Fy2z_a+3*2*I_ESP_Dxz_Fy2z;
    abcd[iGrid*900+123] = 4.0E0*I_ESP_I4x2y_Fy2z_aa-2.0E0*2*I_ESP_G2x2y_Fy2z_a-2.0E0*3*I_ESP_G2x2y_Fy2z_a+2*1*I_ESP_D2y_Fy2z;
    abcd[iGrid*900+124] = 4.0E0*I_ESP_I4xyz_Fy2z_aa-2.0E0*2*I_ESP_G2xyz_Fy2z_a-2.0E0*3*I_ESP_G2xyz_Fy2z_a+2*1*I_ESP_Dyz_Fy2z;
    abcd[iGrid*900+125] = 4.0E0*I_ESP_I4x2z_Fy2z_aa-2.0E0*2*I_ESP_G2x2z_Fy2z_a-2.0E0*3*I_ESP_G2x2z_Fy2z_a+2*1*I_ESP_D2z_Fy2z;
    abcd[iGrid*900+126] = 4.0E0*I_ESP_I3x3y_Fy2z_aa-2.0E0*1*I_ESP_Gx3y_Fy2z_a-2.0E0*2*I_ESP_Gx3y_Fy2z_a;
    abcd[iGrid*900+127] = 4.0E0*I_ESP_I3x2yz_Fy2z_aa-2.0E0*1*I_ESP_Gx2yz_Fy2z_a-2.0E0*2*I_ESP_Gx2yz_Fy2z_a;
    abcd[iGrid*900+128] = 4.0E0*I_ESP_I3xy2z_Fy2z_aa-2.0E0*1*I_ESP_Gxy2z_Fy2z_a-2.0E0*2*I_ESP_Gxy2z_Fy2z_a;
    abcd[iGrid*900+129] = 4.0E0*I_ESP_I3x3z_Fy2z_aa-2.0E0*1*I_ESP_Gx3z_Fy2z_a-2.0E0*2*I_ESP_Gx3z_Fy2z_a;
    abcd[iGrid*900+130] = 4.0E0*I_ESP_I2x4y_Fy2z_aa-2.0E0*1*I_ESP_G4y_Fy2z_a;
    abcd[iGrid*900+131] = 4.0E0*I_ESP_I2x3yz_Fy2z_aa-2.0E0*1*I_ESP_G3yz_Fy2z_a;
    abcd[iGrid*900+132] = 4.0E0*I_ESP_I2x2y2z_Fy2z_aa-2.0E0*1*I_ESP_G2y2z_Fy2z_a;
    abcd[iGrid*900+133] = 4.0E0*I_ESP_I2xy3z_Fy2z_aa-2.0E0*1*I_ESP_Gy3z_Fy2z_a;
    abcd[iGrid*900+134] = 4.0E0*I_ESP_I2x4z_Fy2z_aa-2.0E0*1*I_ESP_G4z_Fy2z_a;
    abcd[iGrid*900+135] = 4.0E0*I_ESP_I6x_F3z_aa-2.0E0*4*I_ESP_G4x_F3z_a-2.0E0*5*I_ESP_G4x_F3z_a+4*3*I_ESP_D2x_F3z;
    abcd[iGrid*900+136] = 4.0E0*I_ESP_I5xy_F3z_aa-2.0E0*3*I_ESP_G3xy_F3z_a-2.0E0*4*I_ESP_G3xy_F3z_a+3*2*I_ESP_Dxy_F3z;
    abcd[iGrid*900+137] = 4.0E0*I_ESP_I5xz_F3z_aa-2.0E0*3*I_ESP_G3xz_F3z_a-2.0E0*4*I_ESP_G3xz_F3z_a+3*2*I_ESP_Dxz_F3z;
    abcd[iGrid*900+138] = 4.0E0*I_ESP_I4x2y_F3z_aa-2.0E0*2*I_ESP_G2x2y_F3z_a-2.0E0*3*I_ESP_G2x2y_F3z_a+2*1*I_ESP_D2y_F3z;
    abcd[iGrid*900+139] = 4.0E0*I_ESP_I4xyz_F3z_aa-2.0E0*2*I_ESP_G2xyz_F3z_a-2.0E0*3*I_ESP_G2xyz_F3z_a+2*1*I_ESP_Dyz_F3z;
    abcd[iGrid*900+140] = 4.0E0*I_ESP_I4x2z_F3z_aa-2.0E0*2*I_ESP_G2x2z_F3z_a-2.0E0*3*I_ESP_G2x2z_F3z_a+2*1*I_ESP_D2z_F3z;
    abcd[iGrid*900+141] = 4.0E0*I_ESP_I3x3y_F3z_aa-2.0E0*1*I_ESP_Gx3y_F3z_a-2.0E0*2*I_ESP_Gx3y_F3z_a;
    abcd[iGrid*900+142] = 4.0E0*I_ESP_I3x2yz_F3z_aa-2.0E0*1*I_ESP_Gx2yz_F3z_a-2.0E0*2*I_ESP_Gx2yz_F3z_a;
    abcd[iGrid*900+143] = 4.0E0*I_ESP_I3xy2z_F3z_aa-2.0E0*1*I_ESP_Gxy2z_F3z_a-2.0E0*2*I_ESP_Gxy2z_F3z_a;
    abcd[iGrid*900+144] = 4.0E0*I_ESP_I3x3z_F3z_aa-2.0E0*1*I_ESP_Gx3z_F3z_a-2.0E0*2*I_ESP_Gx3z_F3z_a;
    abcd[iGrid*900+145] = 4.0E0*I_ESP_I2x4y_F3z_aa-2.0E0*1*I_ESP_G4y_F3z_a;
    abcd[iGrid*900+146] = 4.0E0*I_ESP_I2x3yz_F3z_aa-2.0E0*1*I_ESP_G3yz_F3z_a;
    abcd[iGrid*900+147] = 4.0E0*I_ESP_I2x2y2z_F3z_aa-2.0E0*1*I_ESP_G2y2z_F3z_a;
    abcd[iGrid*900+148] = 4.0E0*I_ESP_I2xy3z_F3z_aa-2.0E0*1*I_ESP_Gy3z_F3z_a;
    abcd[iGrid*900+149] = 4.0E0*I_ESP_I2x4z_F3z_aa-2.0E0*1*I_ESP_G4z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_aa
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_D_F
     ************************************************************/
    abcd[iGrid*900+150] = 4.0E0*I_ESP_I5xy_F3x_aa-2.0E0*4*I_ESP_G3xy_F3x_a;
    abcd[iGrid*900+151] = 4.0E0*I_ESP_I4x2y_F3x_aa-2.0E0*1*I_ESP_G4x_F3x_a-2.0E0*3*I_ESP_G2x2y_F3x_a+3*1*I_ESP_D2x_F3x;
    abcd[iGrid*900+152] = 4.0E0*I_ESP_I4xyz_F3x_aa-2.0E0*3*I_ESP_G2xyz_F3x_a;
    abcd[iGrid*900+153] = 4.0E0*I_ESP_I3x3y_F3x_aa-2.0E0*2*I_ESP_G3xy_F3x_a-2.0E0*2*I_ESP_Gx3y_F3x_a+2*2*I_ESP_Dxy_F3x;
    abcd[iGrid*900+154] = 4.0E0*I_ESP_I3x2yz_F3x_aa-2.0E0*1*I_ESP_G3xz_F3x_a-2.0E0*2*I_ESP_Gx2yz_F3x_a+2*1*I_ESP_Dxz_F3x;
    abcd[iGrid*900+155] = 4.0E0*I_ESP_I3xy2z_F3x_aa-2.0E0*2*I_ESP_Gxy2z_F3x_a;
    abcd[iGrid*900+156] = 4.0E0*I_ESP_I2x4y_F3x_aa-2.0E0*3*I_ESP_G2x2y_F3x_a-2.0E0*1*I_ESP_G4y_F3x_a+3*I_ESP_D2y_F3x;
    abcd[iGrid*900+157] = 4.0E0*I_ESP_I2x3yz_F3x_aa-2.0E0*2*I_ESP_G2xyz_F3x_a-2.0E0*1*I_ESP_G3yz_F3x_a+2*I_ESP_Dyz_F3x;
    abcd[iGrid*900+158] = 4.0E0*I_ESP_I2x2y2z_F3x_aa-2.0E0*1*I_ESP_G2x2z_F3x_a-2.0E0*1*I_ESP_G2y2z_F3x_a+1*I_ESP_D2z_F3x;
    abcd[iGrid*900+159] = 4.0E0*I_ESP_I2xy3z_F3x_aa-2.0E0*1*I_ESP_Gy3z_F3x_a;
    abcd[iGrid*900+160] = 4.0E0*I_ESP_Ix5y_F3x_aa-2.0E0*4*I_ESP_Gx3y_F3x_a;
    abcd[iGrid*900+161] = 4.0E0*I_ESP_Ix4yz_F3x_aa-2.0E0*3*I_ESP_Gx2yz_F3x_a;
    abcd[iGrid*900+162] = 4.0E0*I_ESP_Ix3y2z_F3x_aa-2.0E0*2*I_ESP_Gxy2z_F3x_a;
    abcd[iGrid*900+163] = 4.0E0*I_ESP_Ix2y3z_F3x_aa-2.0E0*1*I_ESP_Gx3z_F3x_a;
    abcd[iGrid*900+164] = 4.0E0*I_ESP_Ixy4z_F3x_aa;
    abcd[iGrid*900+165] = 4.0E0*I_ESP_I5xy_F2xy_aa-2.0E0*4*I_ESP_G3xy_F2xy_a;
    abcd[iGrid*900+166] = 4.0E0*I_ESP_I4x2y_F2xy_aa-2.0E0*1*I_ESP_G4x_F2xy_a-2.0E0*3*I_ESP_G2x2y_F2xy_a+3*1*I_ESP_D2x_F2xy;
    abcd[iGrid*900+167] = 4.0E0*I_ESP_I4xyz_F2xy_aa-2.0E0*3*I_ESP_G2xyz_F2xy_a;
    abcd[iGrid*900+168] = 4.0E0*I_ESP_I3x3y_F2xy_aa-2.0E0*2*I_ESP_G3xy_F2xy_a-2.0E0*2*I_ESP_Gx3y_F2xy_a+2*2*I_ESP_Dxy_F2xy;
    abcd[iGrid*900+169] = 4.0E0*I_ESP_I3x2yz_F2xy_aa-2.0E0*1*I_ESP_G3xz_F2xy_a-2.0E0*2*I_ESP_Gx2yz_F2xy_a+2*1*I_ESP_Dxz_F2xy;
    abcd[iGrid*900+170] = 4.0E0*I_ESP_I3xy2z_F2xy_aa-2.0E0*2*I_ESP_Gxy2z_F2xy_a;
    abcd[iGrid*900+171] = 4.0E0*I_ESP_I2x4y_F2xy_aa-2.0E0*3*I_ESP_G2x2y_F2xy_a-2.0E0*1*I_ESP_G4y_F2xy_a+3*I_ESP_D2y_F2xy;
    abcd[iGrid*900+172] = 4.0E0*I_ESP_I2x3yz_F2xy_aa-2.0E0*2*I_ESP_G2xyz_F2xy_a-2.0E0*1*I_ESP_G3yz_F2xy_a+2*I_ESP_Dyz_F2xy;
    abcd[iGrid*900+173] = 4.0E0*I_ESP_I2x2y2z_F2xy_aa-2.0E0*1*I_ESP_G2x2z_F2xy_a-2.0E0*1*I_ESP_G2y2z_F2xy_a+1*I_ESP_D2z_F2xy;
    abcd[iGrid*900+174] = 4.0E0*I_ESP_I2xy3z_F2xy_aa-2.0E0*1*I_ESP_Gy3z_F2xy_a;
    abcd[iGrid*900+175] = 4.0E0*I_ESP_Ix5y_F2xy_aa-2.0E0*4*I_ESP_Gx3y_F2xy_a;
    abcd[iGrid*900+176] = 4.0E0*I_ESP_Ix4yz_F2xy_aa-2.0E0*3*I_ESP_Gx2yz_F2xy_a;
    abcd[iGrid*900+177] = 4.0E0*I_ESP_Ix3y2z_F2xy_aa-2.0E0*2*I_ESP_Gxy2z_F2xy_a;
    abcd[iGrid*900+178] = 4.0E0*I_ESP_Ix2y3z_F2xy_aa-2.0E0*1*I_ESP_Gx3z_F2xy_a;
    abcd[iGrid*900+179] = 4.0E0*I_ESP_Ixy4z_F2xy_aa;
    abcd[iGrid*900+180] = 4.0E0*I_ESP_I5xy_F2xz_aa-2.0E0*4*I_ESP_G3xy_F2xz_a;
    abcd[iGrid*900+181] = 4.0E0*I_ESP_I4x2y_F2xz_aa-2.0E0*1*I_ESP_G4x_F2xz_a-2.0E0*3*I_ESP_G2x2y_F2xz_a+3*1*I_ESP_D2x_F2xz;
    abcd[iGrid*900+182] = 4.0E0*I_ESP_I4xyz_F2xz_aa-2.0E0*3*I_ESP_G2xyz_F2xz_a;
    abcd[iGrid*900+183] = 4.0E0*I_ESP_I3x3y_F2xz_aa-2.0E0*2*I_ESP_G3xy_F2xz_a-2.0E0*2*I_ESP_Gx3y_F2xz_a+2*2*I_ESP_Dxy_F2xz;
    abcd[iGrid*900+184] = 4.0E0*I_ESP_I3x2yz_F2xz_aa-2.0E0*1*I_ESP_G3xz_F2xz_a-2.0E0*2*I_ESP_Gx2yz_F2xz_a+2*1*I_ESP_Dxz_F2xz;
    abcd[iGrid*900+185] = 4.0E0*I_ESP_I3xy2z_F2xz_aa-2.0E0*2*I_ESP_Gxy2z_F2xz_a;
    abcd[iGrid*900+186] = 4.0E0*I_ESP_I2x4y_F2xz_aa-2.0E0*3*I_ESP_G2x2y_F2xz_a-2.0E0*1*I_ESP_G4y_F2xz_a+3*I_ESP_D2y_F2xz;
    abcd[iGrid*900+187] = 4.0E0*I_ESP_I2x3yz_F2xz_aa-2.0E0*2*I_ESP_G2xyz_F2xz_a-2.0E0*1*I_ESP_G3yz_F2xz_a+2*I_ESP_Dyz_F2xz;
    abcd[iGrid*900+188] = 4.0E0*I_ESP_I2x2y2z_F2xz_aa-2.0E0*1*I_ESP_G2x2z_F2xz_a-2.0E0*1*I_ESP_G2y2z_F2xz_a+1*I_ESP_D2z_F2xz;
    abcd[iGrid*900+189] = 4.0E0*I_ESP_I2xy3z_F2xz_aa-2.0E0*1*I_ESP_Gy3z_F2xz_a;
    abcd[iGrid*900+190] = 4.0E0*I_ESP_Ix5y_F2xz_aa-2.0E0*4*I_ESP_Gx3y_F2xz_a;
    abcd[iGrid*900+191] = 4.0E0*I_ESP_Ix4yz_F2xz_aa-2.0E0*3*I_ESP_Gx2yz_F2xz_a;
    abcd[iGrid*900+192] = 4.0E0*I_ESP_Ix3y2z_F2xz_aa-2.0E0*2*I_ESP_Gxy2z_F2xz_a;
    abcd[iGrid*900+193] = 4.0E0*I_ESP_Ix2y3z_F2xz_aa-2.0E0*1*I_ESP_Gx3z_F2xz_a;
    abcd[iGrid*900+194] = 4.0E0*I_ESP_Ixy4z_F2xz_aa;
    abcd[iGrid*900+195] = 4.0E0*I_ESP_I5xy_Fx2y_aa-2.0E0*4*I_ESP_G3xy_Fx2y_a;
    abcd[iGrid*900+196] = 4.0E0*I_ESP_I4x2y_Fx2y_aa-2.0E0*1*I_ESP_G4x_Fx2y_a-2.0E0*3*I_ESP_G2x2y_Fx2y_a+3*1*I_ESP_D2x_Fx2y;
    abcd[iGrid*900+197] = 4.0E0*I_ESP_I4xyz_Fx2y_aa-2.0E0*3*I_ESP_G2xyz_Fx2y_a;
    abcd[iGrid*900+198] = 4.0E0*I_ESP_I3x3y_Fx2y_aa-2.0E0*2*I_ESP_G3xy_Fx2y_a-2.0E0*2*I_ESP_Gx3y_Fx2y_a+2*2*I_ESP_Dxy_Fx2y;
    abcd[iGrid*900+199] = 4.0E0*I_ESP_I3x2yz_Fx2y_aa-2.0E0*1*I_ESP_G3xz_Fx2y_a-2.0E0*2*I_ESP_Gx2yz_Fx2y_a+2*1*I_ESP_Dxz_Fx2y;
    abcd[iGrid*900+200] = 4.0E0*I_ESP_I3xy2z_Fx2y_aa-2.0E0*2*I_ESP_Gxy2z_Fx2y_a;
    abcd[iGrid*900+201] = 4.0E0*I_ESP_I2x4y_Fx2y_aa-2.0E0*3*I_ESP_G2x2y_Fx2y_a-2.0E0*1*I_ESP_G4y_Fx2y_a+3*I_ESP_D2y_Fx2y;
    abcd[iGrid*900+202] = 4.0E0*I_ESP_I2x3yz_Fx2y_aa-2.0E0*2*I_ESP_G2xyz_Fx2y_a-2.0E0*1*I_ESP_G3yz_Fx2y_a+2*I_ESP_Dyz_Fx2y;
    abcd[iGrid*900+203] = 4.0E0*I_ESP_I2x2y2z_Fx2y_aa-2.0E0*1*I_ESP_G2x2z_Fx2y_a-2.0E0*1*I_ESP_G2y2z_Fx2y_a+1*I_ESP_D2z_Fx2y;
    abcd[iGrid*900+204] = 4.0E0*I_ESP_I2xy3z_Fx2y_aa-2.0E0*1*I_ESP_Gy3z_Fx2y_a;
    abcd[iGrid*900+205] = 4.0E0*I_ESP_Ix5y_Fx2y_aa-2.0E0*4*I_ESP_Gx3y_Fx2y_a;
    abcd[iGrid*900+206] = 4.0E0*I_ESP_Ix4yz_Fx2y_aa-2.0E0*3*I_ESP_Gx2yz_Fx2y_a;
    abcd[iGrid*900+207] = 4.0E0*I_ESP_Ix3y2z_Fx2y_aa-2.0E0*2*I_ESP_Gxy2z_Fx2y_a;
    abcd[iGrid*900+208] = 4.0E0*I_ESP_Ix2y3z_Fx2y_aa-2.0E0*1*I_ESP_Gx3z_Fx2y_a;
    abcd[iGrid*900+209] = 4.0E0*I_ESP_Ixy4z_Fx2y_aa;
    abcd[iGrid*900+210] = 4.0E0*I_ESP_I5xy_Fxyz_aa-2.0E0*4*I_ESP_G3xy_Fxyz_a;
    abcd[iGrid*900+211] = 4.0E0*I_ESP_I4x2y_Fxyz_aa-2.0E0*1*I_ESP_G4x_Fxyz_a-2.0E0*3*I_ESP_G2x2y_Fxyz_a+3*1*I_ESP_D2x_Fxyz;
    abcd[iGrid*900+212] = 4.0E0*I_ESP_I4xyz_Fxyz_aa-2.0E0*3*I_ESP_G2xyz_Fxyz_a;
    abcd[iGrid*900+213] = 4.0E0*I_ESP_I3x3y_Fxyz_aa-2.0E0*2*I_ESP_G3xy_Fxyz_a-2.0E0*2*I_ESP_Gx3y_Fxyz_a+2*2*I_ESP_Dxy_Fxyz;
    abcd[iGrid*900+214] = 4.0E0*I_ESP_I3x2yz_Fxyz_aa-2.0E0*1*I_ESP_G3xz_Fxyz_a-2.0E0*2*I_ESP_Gx2yz_Fxyz_a+2*1*I_ESP_Dxz_Fxyz;
    abcd[iGrid*900+215] = 4.0E0*I_ESP_I3xy2z_Fxyz_aa-2.0E0*2*I_ESP_Gxy2z_Fxyz_a;
    abcd[iGrid*900+216] = 4.0E0*I_ESP_I2x4y_Fxyz_aa-2.0E0*3*I_ESP_G2x2y_Fxyz_a-2.0E0*1*I_ESP_G4y_Fxyz_a+3*I_ESP_D2y_Fxyz;
    abcd[iGrid*900+217] = 4.0E0*I_ESP_I2x3yz_Fxyz_aa-2.0E0*2*I_ESP_G2xyz_Fxyz_a-2.0E0*1*I_ESP_G3yz_Fxyz_a+2*I_ESP_Dyz_Fxyz;
    abcd[iGrid*900+218] = 4.0E0*I_ESP_I2x2y2z_Fxyz_aa-2.0E0*1*I_ESP_G2x2z_Fxyz_a-2.0E0*1*I_ESP_G2y2z_Fxyz_a+1*I_ESP_D2z_Fxyz;
    abcd[iGrid*900+219] = 4.0E0*I_ESP_I2xy3z_Fxyz_aa-2.0E0*1*I_ESP_Gy3z_Fxyz_a;
    abcd[iGrid*900+220] = 4.0E0*I_ESP_Ix5y_Fxyz_aa-2.0E0*4*I_ESP_Gx3y_Fxyz_a;
    abcd[iGrid*900+221] = 4.0E0*I_ESP_Ix4yz_Fxyz_aa-2.0E0*3*I_ESP_Gx2yz_Fxyz_a;
    abcd[iGrid*900+222] = 4.0E0*I_ESP_Ix3y2z_Fxyz_aa-2.0E0*2*I_ESP_Gxy2z_Fxyz_a;
    abcd[iGrid*900+223] = 4.0E0*I_ESP_Ix2y3z_Fxyz_aa-2.0E0*1*I_ESP_Gx3z_Fxyz_a;
    abcd[iGrid*900+224] = 4.0E0*I_ESP_Ixy4z_Fxyz_aa;
    abcd[iGrid*900+225] = 4.0E0*I_ESP_I5xy_Fx2z_aa-2.0E0*4*I_ESP_G3xy_Fx2z_a;
    abcd[iGrid*900+226] = 4.0E0*I_ESP_I4x2y_Fx2z_aa-2.0E0*1*I_ESP_G4x_Fx2z_a-2.0E0*3*I_ESP_G2x2y_Fx2z_a+3*1*I_ESP_D2x_Fx2z;
    abcd[iGrid*900+227] = 4.0E0*I_ESP_I4xyz_Fx2z_aa-2.0E0*3*I_ESP_G2xyz_Fx2z_a;
    abcd[iGrid*900+228] = 4.0E0*I_ESP_I3x3y_Fx2z_aa-2.0E0*2*I_ESP_G3xy_Fx2z_a-2.0E0*2*I_ESP_Gx3y_Fx2z_a+2*2*I_ESP_Dxy_Fx2z;
    abcd[iGrid*900+229] = 4.0E0*I_ESP_I3x2yz_Fx2z_aa-2.0E0*1*I_ESP_G3xz_Fx2z_a-2.0E0*2*I_ESP_Gx2yz_Fx2z_a+2*1*I_ESP_Dxz_Fx2z;
    abcd[iGrid*900+230] = 4.0E0*I_ESP_I3xy2z_Fx2z_aa-2.0E0*2*I_ESP_Gxy2z_Fx2z_a;
    abcd[iGrid*900+231] = 4.0E0*I_ESP_I2x4y_Fx2z_aa-2.0E0*3*I_ESP_G2x2y_Fx2z_a-2.0E0*1*I_ESP_G4y_Fx2z_a+3*I_ESP_D2y_Fx2z;
    abcd[iGrid*900+232] = 4.0E0*I_ESP_I2x3yz_Fx2z_aa-2.0E0*2*I_ESP_G2xyz_Fx2z_a-2.0E0*1*I_ESP_G3yz_Fx2z_a+2*I_ESP_Dyz_Fx2z;
    abcd[iGrid*900+233] = 4.0E0*I_ESP_I2x2y2z_Fx2z_aa-2.0E0*1*I_ESP_G2x2z_Fx2z_a-2.0E0*1*I_ESP_G2y2z_Fx2z_a+1*I_ESP_D2z_Fx2z;
    abcd[iGrid*900+234] = 4.0E0*I_ESP_I2xy3z_Fx2z_aa-2.0E0*1*I_ESP_Gy3z_Fx2z_a;
    abcd[iGrid*900+235] = 4.0E0*I_ESP_Ix5y_Fx2z_aa-2.0E0*4*I_ESP_Gx3y_Fx2z_a;
    abcd[iGrid*900+236] = 4.0E0*I_ESP_Ix4yz_Fx2z_aa-2.0E0*3*I_ESP_Gx2yz_Fx2z_a;
    abcd[iGrid*900+237] = 4.0E0*I_ESP_Ix3y2z_Fx2z_aa-2.0E0*2*I_ESP_Gxy2z_Fx2z_a;
    abcd[iGrid*900+238] = 4.0E0*I_ESP_Ix2y3z_Fx2z_aa-2.0E0*1*I_ESP_Gx3z_Fx2z_a;
    abcd[iGrid*900+239] = 4.0E0*I_ESP_Ixy4z_Fx2z_aa;
    abcd[iGrid*900+240] = 4.0E0*I_ESP_I5xy_F3y_aa-2.0E0*4*I_ESP_G3xy_F3y_a;
    abcd[iGrid*900+241] = 4.0E0*I_ESP_I4x2y_F3y_aa-2.0E0*1*I_ESP_G4x_F3y_a-2.0E0*3*I_ESP_G2x2y_F3y_a+3*1*I_ESP_D2x_F3y;
    abcd[iGrid*900+242] = 4.0E0*I_ESP_I4xyz_F3y_aa-2.0E0*3*I_ESP_G2xyz_F3y_a;
    abcd[iGrid*900+243] = 4.0E0*I_ESP_I3x3y_F3y_aa-2.0E0*2*I_ESP_G3xy_F3y_a-2.0E0*2*I_ESP_Gx3y_F3y_a+2*2*I_ESP_Dxy_F3y;
    abcd[iGrid*900+244] = 4.0E0*I_ESP_I3x2yz_F3y_aa-2.0E0*1*I_ESP_G3xz_F3y_a-2.0E0*2*I_ESP_Gx2yz_F3y_a+2*1*I_ESP_Dxz_F3y;
    abcd[iGrid*900+245] = 4.0E0*I_ESP_I3xy2z_F3y_aa-2.0E0*2*I_ESP_Gxy2z_F3y_a;
    abcd[iGrid*900+246] = 4.0E0*I_ESP_I2x4y_F3y_aa-2.0E0*3*I_ESP_G2x2y_F3y_a-2.0E0*1*I_ESP_G4y_F3y_a+3*I_ESP_D2y_F3y;
    abcd[iGrid*900+247] = 4.0E0*I_ESP_I2x3yz_F3y_aa-2.0E0*2*I_ESP_G2xyz_F3y_a-2.0E0*1*I_ESP_G3yz_F3y_a+2*I_ESP_Dyz_F3y;
    abcd[iGrid*900+248] = 4.0E0*I_ESP_I2x2y2z_F3y_aa-2.0E0*1*I_ESP_G2x2z_F3y_a-2.0E0*1*I_ESP_G2y2z_F3y_a+1*I_ESP_D2z_F3y;
    abcd[iGrid*900+249] = 4.0E0*I_ESP_I2xy3z_F3y_aa-2.0E0*1*I_ESP_Gy3z_F3y_a;
    abcd[iGrid*900+250] = 4.0E0*I_ESP_Ix5y_F3y_aa-2.0E0*4*I_ESP_Gx3y_F3y_a;
    abcd[iGrid*900+251] = 4.0E0*I_ESP_Ix4yz_F3y_aa-2.0E0*3*I_ESP_Gx2yz_F3y_a;
    abcd[iGrid*900+252] = 4.0E0*I_ESP_Ix3y2z_F3y_aa-2.0E0*2*I_ESP_Gxy2z_F3y_a;
    abcd[iGrid*900+253] = 4.0E0*I_ESP_Ix2y3z_F3y_aa-2.0E0*1*I_ESP_Gx3z_F3y_a;
    abcd[iGrid*900+254] = 4.0E0*I_ESP_Ixy4z_F3y_aa;
    abcd[iGrid*900+255] = 4.0E0*I_ESP_I5xy_F2yz_aa-2.0E0*4*I_ESP_G3xy_F2yz_a;
    abcd[iGrid*900+256] = 4.0E0*I_ESP_I4x2y_F2yz_aa-2.0E0*1*I_ESP_G4x_F2yz_a-2.0E0*3*I_ESP_G2x2y_F2yz_a+3*1*I_ESP_D2x_F2yz;
    abcd[iGrid*900+257] = 4.0E0*I_ESP_I4xyz_F2yz_aa-2.0E0*3*I_ESP_G2xyz_F2yz_a;
    abcd[iGrid*900+258] = 4.0E0*I_ESP_I3x3y_F2yz_aa-2.0E0*2*I_ESP_G3xy_F2yz_a-2.0E0*2*I_ESP_Gx3y_F2yz_a+2*2*I_ESP_Dxy_F2yz;
    abcd[iGrid*900+259] = 4.0E0*I_ESP_I3x2yz_F2yz_aa-2.0E0*1*I_ESP_G3xz_F2yz_a-2.0E0*2*I_ESP_Gx2yz_F2yz_a+2*1*I_ESP_Dxz_F2yz;
    abcd[iGrid*900+260] = 4.0E0*I_ESP_I3xy2z_F2yz_aa-2.0E0*2*I_ESP_Gxy2z_F2yz_a;
    abcd[iGrid*900+261] = 4.0E0*I_ESP_I2x4y_F2yz_aa-2.0E0*3*I_ESP_G2x2y_F2yz_a-2.0E0*1*I_ESP_G4y_F2yz_a+3*I_ESP_D2y_F2yz;
    abcd[iGrid*900+262] = 4.0E0*I_ESP_I2x3yz_F2yz_aa-2.0E0*2*I_ESP_G2xyz_F2yz_a-2.0E0*1*I_ESP_G3yz_F2yz_a+2*I_ESP_Dyz_F2yz;
    abcd[iGrid*900+263] = 4.0E0*I_ESP_I2x2y2z_F2yz_aa-2.0E0*1*I_ESP_G2x2z_F2yz_a-2.0E0*1*I_ESP_G2y2z_F2yz_a+1*I_ESP_D2z_F2yz;
    abcd[iGrid*900+264] = 4.0E0*I_ESP_I2xy3z_F2yz_aa-2.0E0*1*I_ESP_Gy3z_F2yz_a;
    abcd[iGrid*900+265] = 4.0E0*I_ESP_Ix5y_F2yz_aa-2.0E0*4*I_ESP_Gx3y_F2yz_a;
    abcd[iGrid*900+266] = 4.0E0*I_ESP_Ix4yz_F2yz_aa-2.0E0*3*I_ESP_Gx2yz_F2yz_a;
    abcd[iGrid*900+267] = 4.0E0*I_ESP_Ix3y2z_F2yz_aa-2.0E0*2*I_ESP_Gxy2z_F2yz_a;
    abcd[iGrid*900+268] = 4.0E0*I_ESP_Ix2y3z_F2yz_aa-2.0E0*1*I_ESP_Gx3z_F2yz_a;
    abcd[iGrid*900+269] = 4.0E0*I_ESP_Ixy4z_F2yz_aa;
    abcd[iGrid*900+270] = 4.0E0*I_ESP_I5xy_Fy2z_aa-2.0E0*4*I_ESP_G3xy_Fy2z_a;
    abcd[iGrid*900+271] = 4.0E0*I_ESP_I4x2y_Fy2z_aa-2.0E0*1*I_ESP_G4x_Fy2z_a-2.0E0*3*I_ESP_G2x2y_Fy2z_a+3*1*I_ESP_D2x_Fy2z;
    abcd[iGrid*900+272] = 4.0E0*I_ESP_I4xyz_Fy2z_aa-2.0E0*3*I_ESP_G2xyz_Fy2z_a;
    abcd[iGrid*900+273] = 4.0E0*I_ESP_I3x3y_Fy2z_aa-2.0E0*2*I_ESP_G3xy_Fy2z_a-2.0E0*2*I_ESP_Gx3y_Fy2z_a+2*2*I_ESP_Dxy_Fy2z;
    abcd[iGrid*900+274] = 4.0E0*I_ESP_I3x2yz_Fy2z_aa-2.0E0*1*I_ESP_G3xz_Fy2z_a-2.0E0*2*I_ESP_Gx2yz_Fy2z_a+2*1*I_ESP_Dxz_Fy2z;
    abcd[iGrid*900+275] = 4.0E0*I_ESP_I3xy2z_Fy2z_aa-2.0E0*2*I_ESP_Gxy2z_Fy2z_a;
    abcd[iGrid*900+276] = 4.0E0*I_ESP_I2x4y_Fy2z_aa-2.0E0*3*I_ESP_G2x2y_Fy2z_a-2.0E0*1*I_ESP_G4y_Fy2z_a+3*I_ESP_D2y_Fy2z;
    abcd[iGrid*900+277] = 4.0E0*I_ESP_I2x3yz_Fy2z_aa-2.0E0*2*I_ESP_G2xyz_Fy2z_a-2.0E0*1*I_ESP_G3yz_Fy2z_a+2*I_ESP_Dyz_Fy2z;
    abcd[iGrid*900+278] = 4.0E0*I_ESP_I2x2y2z_Fy2z_aa-2.0E0*1*I_ESP_G2x2z_Fy2z_a-2.0E0*1*I_ESP_G2y2z_Fy2z_a+1*I_ESP_D2z_Fy2z;
    abcd[iGrid*900+279] = 4.0E0*I_ESP_I2xy3z_Fy2z_aa-2.0E0*1*I_ESP_Gy3z_Fy2z_a;
    abcd[iGrid*900+280] = 4.0E0*I_ESP_Ix5y_Fy2z_aa-2.0E0*4*I_ESP_Gx3y_Fy2z_a;
    abcd[iGrid*900+281] = 4.0E0*I_ESP_Ix4yz_Fy2z_aa-2.0E0*3*I_ESP_Gx2yz_Fy2z_a;
    abcd[iGrid*900+282] = 4.0E0*I_ESP_Ix3y2z_Fy2z_aa-2.0E0*2*I_ESP_Gxy2z_Fy2z_a;
    abcd[iGrid*900+283] = 4.0E0*I_ESP_Ix2y3z_Fy2z_aa-2.0E0*1*I_ESP_Gx3z_Fy2z_a;
    abcd[iGrid*900+284] = 4.0E0*I_ESP_Ixy4z_Fy2z_aa;
    abcd[iGrid*900+285] = 4.0E0*I_ESP_I5xy_F3z_aa-2.0E0*4*I_ESP_G3xy_F3z_a;
    abcd[iGrid*900+286] = 4.0E0*I_ESP_I4x2y_F3z_aa-2.0E0*1*I_ESP_G4x_F3z_a-2.0E0*3*I_ESP_G2x2y_F3z_a+3*1*I_ESP_D2x_F3z;
    abcd[iGrid*900+287] = 4.0E0*I_ESP_I4xyz_F3z_aa-2.0E0*3*I_ESP_G2xyz_F3z_a;
    abcd[iGrid*900+288] = 4.0E0*I_ESP_I3x3y_F3z_aa-2.0E0*2*I_ESP_G3xy_F3z_a-2.0E0*2*I_ESP_Gx3y_F3z_a+2*2*I_ESP_Dxy_F3z;
    abcd[iGrid*900+289] = 4.0E0*I_ESP_I3x2yz_F3z_aa-2.0E0*1*I_ESP_G3xz_F3z_a-2.0E0*2*I_ESP_Gx2yz_F3z_a+2*1*I_ESP_Dxz_F3z;
    abcd[iGrid*900+290] = 4.0E0*I_ESP_I3xy2z_F3z_aa-2.0E0*2*I_ESP_Gxy2z_F3z_a;
    abcd[iGrid*900+291] = 4.0E0*I_ESP_I2x4y_F3z_aa-2.0E0*3*I_ESP_G2x2y_F3z_a-2.0E0*1*I_ESP_G4y_F3z_a+3*I_ESP_D2y_F3z;
    abcd[iGrid*900+292] = 4.0E0*I_ESP_I2x3yz_F3z_aa-2.0E0*2*I_ESP_G2xyz_F3z_a-2.0E0*1*I_ESP_G3yz_F3z_a+2*I_ESP_Dyz_F3z;
    abcd[iGrid*900+293] = 4.0E0*I_ESP_I2x2y2z_F3z_aa-2.0E0*1*I_ESP_G2x2z_F3z_a-2.0E0*1*I_ESP_G2y2z_F3z_a+1*I_ESP_D2z_F3z;
    abcd[iGrid*900+294] = 4.0E0*I_ESP_I2xy3z_F3z_aa-2.0E0*1*I_ESP_Gy3z_F3z_a;
    abcd[iGrid*900+295] = 4.0E0*I_ESP_Ix5y_F3z_aa-2.0E0*4*I_ESP_Gx3y_F3z_a;
    abcd[iGrid*900+296] = 4.0E0*I_ESP_Ix4yz_F3z_aa-2.0E0*3*I_ESP_Gx2yz_F3z_a;
    abcd[iGrid*900+297] = 4.0E0*I_ESP_Ix3y2z_F3z_aa-2.0E0*2*I_ESP_Gxy2z_F3z_a;
    abcd[iGrid*900+298] = 4.0E0*I_ESP_Ix2y3z_F3z_aa-2.0E0*1*I_ESP_Gx3z_F3z_a;
    abcd[iGrid*900+299] = 4.0E0*I_ESP_Ixy4z_F3z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_aa
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_D_F
     ************************************************************/
    abcd[iGrid*900+300] = 4.0E0*I_ESP_I5xz_F3x_aa-2.0E0*4*I_ESP_G3xz_F3x_a;
    abcd[iGrid*900+301] = 4.0E0*I_ESP_I4xyz_F3x_aa-2.0E0*3*I_ESP_G2xyz_F3x_a;
    abcd[iGrid*900+302] = 4.0E0*I_ESP_I4x2z_F3x_aa-2.0E0*1*I_ESP_G4x_F3x_a-2.0E0*3*I_ESP_G2x2z_F3x_a+3*1*I_ESP_D2x_F3x;
    abcd[iGrid*900+303] = 4.0E0*I_ESP_I3x2yz_F3x_aa-2.0E0*2*I_ESP_Gx2yz_F3x_a;
    abcd[iGrid*900+304] = 4.0E0*I_ESP_I3xy2z_F3x_aa-2.0E0*1*I_ESP_G3xy_F3x_a-2.0E0*2*I_ESP_Gxy2z_F3x_a+2*1*I_ESP_Dxy_F3x;
    abcd[iGrid*900+305] = 4.0E0*I_ESP_I3x3z_F3x_aa-2.0E0*2*I_ESP_G3xz_F3x_a-2.0E0*2*I_ESP_Gx3z_F3x_a+2*2*I_ESP_Dxz_F3x;
    abcd[iGrid*900+306] = 4.0E0*I_ESP_I2x3yz_F3x_aa-2.0E0*1*I_ESP_G3yz_F3x_a;
    abcd[iGrid*900+307] = 4.0E0*I_ESP_I2x2y2z_F3x_aa-2.0E0*1*I_ESP_G2x2y_F3x_a-2.0E0*1*I_ESP_G2y2z_F3x_a+1*I_ESP_D2y_F3x;
    abcd[iGrid*900+308] = 4.0E0*I_ESP_I2xy3z_F3x_aa-2.0E0*2*I_ESP_G2xyz_F3x_a-2.0E0*1*I_ESP_Gy3z_F3x_a+2*I_ESP_Dyz_F3x;
    abcd[iGrid*900+309] = 4.0E0*I_ESP_I2x4z_F3x_aa-2.0E0*3*I_ESP_G2x2z_F3x_a-2.0E0*1*I_ESP_G4z_F3x_a+3*I_ESP_D2z_F3x;
    abcd[iGrid*900+310] = 4.0E0*I_ESP_Ix4yz_F3x_aa;
    abcd[iGrid*900+311] = 4.0E0*I_ESP_Ix3y2z_F3x_aa-2.0E0*1*I_ESP_Gx3y_F3x_a;
    abcd[iGrid*900+312] = 4.0E0*I_ESP_Ix2y3z_F3x_aa-2.0E0*2*I_ESP_Gx2yz_F3x_a;
    abcd[iGrid*900+313] = 4.0E0*I_ESP_Ixy4z_F3x_aa-2.0E0*3*I_ESP_Gxy2z_F3x_a;
    abcd[iGrid*900+314] = 4.0E0*I_ESP_Ix5z_F3x_aa-2.0E0*4*I_ESP_Gx3z_F3x_a;
    abcd[iGrid*900+315] = 4.0E0*I_ESP_I5xz_F2xy_aa-2.0E0*4*I_ESP_G3xz_F2xy_a;
    abcd[iGrid*900+316] = 4.0E0*I_ESP_I4xyz_F2xy_aa-2.0E0*3*I_ESP_G2xyz_F2xy_a;
    abcd[iGrid*900+317] = 4.0E0*I_ESP_I4x2z_F2xy_aa-2.0E0*1*I_ESP_G4x_F2xy_a-2.0E0*3*I_ESP_G2x2z_F2xy_a+3*1*I_ESP_D2x_F2xy;
    abcd[iGrid*900+318] = 4.0E0*I_ESP_I3x2yz_F2xy_aa-2.0E0*2*I_ESP_Gx2yz_F2xy_a;
    abcd[iGrid*900+319] = 4.0E0*I_ESP_I3xy2z_F2xy_aa-2.0E0*1*I_ESP_G3xy_F2xy_a-2.0E0*2*I_ESP_Gxy2z_F2xy_a+2*1*I_ESP_Dxy_F2xy;
    abcd[iGrid*900+320] = 4.0E0*I_ESP_I3x3z_F2xy_aa-2.0E0*2*I_ESP_G3xz_F2xy_a-2.0E0*2*I_ESP_Gx3z_F2xy_a+2*2*I_ESP_Dxz_F2xy;
    abcd[iGrid*900+321] = 4.0E0*I_ESP_I2x3yz_F2xy_aa-2.0E0*1*I_ESP_G3yz_F2xy_a;
    abcd[iGrid*900+322] = 4.0E0*I_ESP_I2x2y2z_F2xy_aa-2.0E0*1*I_ESP_G2x2y_F2xy_a-2.0E0*1*I_ESP_G2y2z_F2xy_a+1*I_ESP_D2y_F2xy;
    abcd[iGrid*900+323] = 4.0E0*I_ESP_I2xy3z_F2xy_aa-2.0E0*2*I_ESP_G2xyz_F2xy_a-2.0E0*1*I_ESP_Gy3z_F2xy_a+2*I_ESP_Dyz_F2xy;
    abcd[iGrid*900+324] = 4.0E0*I_ESP_I2x4z_F2xy_aa-2.0E0*3*I_ESP_G2x2z_F2xy_a-2.0E0*1*I_ESP_G4z_F2xy_a+3*I_ESP_D2z_F2xy;
    abcd[iGrid*900+325] = 4.0E0*I_ESP_Ix4yz_F2xy_aa;
    abcd[iGrid*900+326] = 4.0E0*I_ESP_Ix3y2z_F2xy_aa-2.0E0*1*I_ESP_Gx3y_F2xy_a;
    abcd[iGrid*900+327] = 4.0E0*I_ESP_Ix2y3z_F2xy_aa-2.0E0*2*I_ESP_Gx2yz_F2xy_a;
    abcd[iGrid*900+328] = 4.0E0*I_ESP_Ixy4z_F2xy_aa-2.0E0*3*I_ESP_Gxy2z_F2xy_a;
    abcd[iGrid*900+329] = 4.0E0*I_ESP_Ix5z_F2xy_aa-2.0E0*4*I_ESP_Gx3z_F2xy_a;
    abcd[iGrid*900+330] = 4.0E0*I_ESP_I5xz_F2xz_aa-2.0E0*4*I_ESP_G3xz_F2xz_a;
    abcd[iGrid*900+331] = 4.0E0*I_ESP_I4xyz_F2xz_aa-2.0E0*3*I_ESP_G2xyz_F2xz_a;
    abcd[iGrid*900+332] = 4.0E0*I_ESP_I4x2z_F2xz_aa-2.0E0*1*I_ESP_G4x_F2xz_a-2.0E0*3*I_ESP_G2x2z_F2xz_a+3*1*I_ESP_D2x_F2xz;
    abcd[iGrid*900+333] = 4.0E0*I_ESP_I3x2yz_F2xz_aa-2.0E0*2*I_ESP_Gx2yz_F2xz_a;
    abcd[iGrid*900+334] = 4.0E0*I_ESP_I3xy2z_F2xz_aa-2.0E0*1*I_ESP_G3xy_F2xz_a-2.0E0*2*I_ESP_Gxy2z_F2xz_a+2*1*I_ESP_Dxy_F2xz;
    abcd[iGrid*900+335] = 4.0E0*I_ESP_I3x3z_F2xz_aa-2.0E0*2*I_ESP_G3xz_F2xz_a-2.0E0*2*I_ESP_Gx3z_F2xz_a+2*2*I_ESP_Dxz_F2xz;
    abcd[iGrid*900+336] = 4.0E0*I_ESP_I2x3yz_F2xz_aa-2.0E0*1*I_ESP_G3yz_F2xz_a;
    abcd[iGrid*900+337] = 4.0E0*I_ESP_I2x2y2z_F2xz_aa-2.0E0*1*I_ESP_G2x2y_F2xz_a-2.0E0*1*I_ESP_G2y2z_F2xz_a+1*I_ESP_D2y_F2xz;
    abcd[iGrid*900+338] = 4.0E0*I_ESP_I2xy3z_F2xz_aa-2.0E0*2*I_ESP_G2xyz_F2xz_a-2.0E0*1*I_ESP_Gy3z_F2xz_a+2*I_ESP_Dyz_F2xz;
    abcd[iGrid*900+339] = 4.0E0*I_ESP_I2x4z_F2xz_aa-2.0E0*3*I_ESP_G2x2z_F2xz_a-2.0E0*1*I_ESP_G4z_F2xz_a+3*I_ESP_D2z_F2xz;
    abcd[iGrid*900+340] = 4.0E0*I_ESP_Ix4yz_F2xz_aa;
    abcd[iGrid*900+341] = 4.0E0*I_ESP_Ix3y2z_F2xz_aa-2.0E0*1*I_ESP_Gx3y_F2xz_a;
    abcd[iGrid*900+342] = 4.0E0*I_ESP_Ix2y3z_F2xz_aa-2.0E0*2*I_ESP_Gx2yz_F2xz_a;
    abcd[iGrid*900+343] = 4.0E0*I_ESP_Ixy4z_F2xz_aa-2.0E0*3*I_ESP_Gxy2z_F2xz_a;
    abcd[iGrid*900+344] = 4.0E0*I_ESP_Ix5z_F2xz_aa-2.0E0*4*I_ESP_Gx3z_F2xz_a;
    abcd[iGrid*900+345] = 4.0E0*I_ESP_I5xz_Fx2y_aa-2.0E0*4*I_ESP_G3xz_Fx2y_a;
    abcd[iGrid*900+346] = 4.0E0*I_ESP_I4xyz_Fx2y_aa-2.0E0*3*I_ESP_G2xyz_Fx2y_a;
    abcd[iGrid*900+347] = 4.0E0*I_ESP_I4x2z_Fx2y_aa-2.0E0*1*I_ESP_G4x_Fx2y_a-2.0E0*3*I_ESP_G2x2z_Fx2y_a+3*1*I_ESP_D2x_Fx2y;
    abcd[iGrid*900+348] = 4.0E0*I_ESP_I3x2yz_Fx2y_aa-2.0E0*2*I_ESP_Gx2yz_Fx2y_a;
    abcd[iGrid*900+349] = 4.0E0*I_ESP_I3xy2z_Fx2y_aa-2.0E0*1*I_ESP_G3xy_Fx2y_a-2.0E0*2*I_ESP_Gxy2z_Fx2y_a+2*1*I_ESP_Dxy_Fx2y;
    abcd[iGrid*900+350] = 4.0E0*I_ESP_I3x3z_Fx2y_aa-2.0E0*2*I_ESP_G3xz_Fx2y_a-2.0E0*2*I_ESP_Gx3z_Fx2y_a+2*2*I_ESP_Dxz_Fx2y;
    abcd[iGrid*900+351] = 4.0E0*I_ESP_I2x3yz_Fx2y_aa-2.0E0*1*I_ESP_G3yz_Fx2y_a;
    abcd[iGrid*900+352] = 4.0E0*I_ESP_I2x2y2z_Fx2y_aa-2.0E0*1*I_ESP_G2x2y_Fx2y_a-2.0E0*1*I_ESP_G2y2z_Fx2y_a+1*I_ESP_D2y_Fx2y;
    abcd[iGrid*900+353] = 4.0E0*I_ESP_I2xy3z_Fx2y_aa-2.0E0*2*I_ESP_G2xyz_Fx2y_a-2.0E0*1*I_ESP_Gy3z_Fx2y_a+2*I_ESP_Dyz_Fx2y;
    abcd[iGrid*900+354] = 4.0E0*I_ESP_I2x4z_Fx2y_aa-2.0E0*3*I_ESP_G2x2z_Fx2y_a-2.0E0*1*I_ESP_G4z_Fx2y_a+3*I_ESP_D2z_Fx2y;
    abcd[iGrid*900+355] = 4.0E0*I_ESP_Ix4yz_Fx2y_aa;
    abcd[iGrid*900+356] = 4.0E0*I_ESP_Ix3y2z_Fx2y_aa-2.0E0*1*I_ESP_Gx3y_Fx2y_a;
    abcd[iGrid*900+357] = 4.0E0*I_ESP_Ix2y3z_Fx2y_aa-2.0E0*2*I_ESP_Gx2yz_Fx2y_a;
    abcd[iGrid*900+358] = 4.0E0*I_ESP_Ixy4z_Fx2y_aa-2.0E0*3*I_ESP_Gxy2z_Fx2y_a;
    abcd[iGrid*900+359] = 4.0E0*I_ESP_Ix5z_Fx2y_aa-2.0E0*4*I_ESP_Gx3z_Fx2y_a;
    abcd[iGrid*900+360] = 4.0E0*I_ESP_I5xz_Fxyz_aa-2.0E0*4*I_ESP_G3xz_Fxyz_a;
    abcd[iGrid*900+361] = 4.0E0*I_ESP_I4xyz_Fxyz_aa-2.0E0*3*I_ESP_G2xyz_Fxyz_a;
    abcd[iGrid*900+362] = 4.0E0*I_ESP_I4x2z_Fxyz_aa-2.0E0*1*I_ESP_G4x_Fxyz_a-2.0E0*3*I_ESP_G2x2z_Fxyz_a+3*1*I_ESP_D2x_Fxyz;
    abcd[iGrid*900+363] = 4.0E0*I_ESP_I3x2yz_Fxyz_aa-2.0E0*2*I_ESP_Gx2yz_Fxyz_a;
    abcd[iGrid*900+364] = 4.0E0*I_ESP_I3xy2z_Fxyz_aa-2.0E0*1*I_ESP_G3xy_Fxyz_a-2.0E0*2*I_ESP_Gxy2z_Fxyz_a+2*1*I_ESP_Dxy_Fxyz;
    abcd[iGrid*900+365] = 4.0E0*I_ESP_I3x3z_Fxyz_aa-2.0E0*2*I_ESP_G3xz_Fxyz_a-2.0E0*2*I_ESP_Gx3z_Fxyz_a+2*2*I_ESP_Dxz_Fxyz;
    abcd[iGrid*900+366] = 4.0E0*I_ESP_I2x3yz_Fxyz_aa-2.0E0*1*I_ESP_G3yz_Fxyz_a;
    abcd[iGrid*900+367] = 4.0E0*I_ESP_I2x2y2z_Fxyz_aa-2.0E0*1*I_ESP_G2x2y_Fxyz_a-2.0E0*1*I_ESP_G2y2z_Fxyz_a+1*I_ESP_D2y_Fxyz;
    abcd[iGrid*900+368] = 4.0E0*I_ESP_I2xy3z_Fxyz_aa-2.0E0*2*I_ESP_G2xyz_Fxyz_a-2.0E0*1*I_ESP_Gy3z_Fxyz_a+2*I_ESP_Dyz_Fxyz;
    abcd[iGrid*900+369] = 4.0E0*I_ESP_I2x4z_Fxyz_aa-2.0E0*3*I_ESP_G2x2z_Fxyz_a-2.0E0*1*I_ESP_G4z_Fxyz_a+3*I_ESP_D2z_Fxyz;
    abcd[iGrid*900+370] = 4.0E0*I_ESP_Ix4yz_Fxyz_aa;
    abcd[iGrid*900+371] = 4.0E0*I_ESP_Ix3y2z_Fxyz_aa-2.0E0*1*I_ESP_Gx3y_Fxyz_a;
    abcd[iGrid*900+372] = 4.0E0*I_ESP_Ix2y3z_Fxyz_aa-2.0E0*2*I_ESP_Gx2yz_Fxyz_a;
    abcd[iGrid*900+373] = 4.0E0*I_ESP_Ixy4z_Fxyz_aa-2.0E0*3*I_ESP_Gxy2z_Fxyz_a;
    abcd[iGrid*900+374] = 4.0E0*I_ESP_Ix5z_Fxyz_aa-2.0E0*4*I_ESP_Gx3z_Fxyz_a;
    abcd[iGrid*900+375] = 4.0E0*I_ESP_I5xz_Fx2z_aa-2.0E0*4*I_ESP_G3xz_Fx2z_a;
    abcd[iGrid*900+376] = 4.0E0*I_ESP_I4xyz_Fx2z_aa-2.0E0*3*I_ESP_G2xyz_Fx2z_a;
    abcd[iGrid*900+377] = 4.0E0*I_ESP_I4x2z_Fx2z_aa-2.0E0*1*I_ESP_G4x_Fx2z_a-2.0E0*3*I_ESP_G2x2z_Fx2z_a+3*1*I_ESP_D2x_Fx2z;
    abcd[iGrid*900+378] = 4.0E0*I_ESP_I3x2yz_Fx2z_aa-2.0E0*2*I_ESP_Gx2yz_Fx2z_a;
    abcd[iGrid*900+379] = 4.0E0*I_ESP_I3xy2z_Fx2z_aa-2.0E0*1*I_ESP_G3xy_Fx2z_a-2.0E0*2*I_ESP_Gxy2z_Fx2z_a+2*1*I_ESP_Dxy_Fx2z;
    abcd[iGrid*900+380] = 4.0E0*I_ESP_I3x3z_Fx2z_aa-2.0E0*2*I_ESP_G3xz_Fx2z_a-2.0E0*2*I_ESP_Gx3z_Fx2z_a+2*2*I_ESP_Dxz_Fx2z;
    abcd[iGrid*900+381] = 4.0E0*I_ESP_I2x3yz_Fx2z_aa-2.0E0*1*I_ESP_G3yz_Fx2z_a;
    abcd[iGrid*900+382] = 4.0E0*I_ESP_I2x2y2z_Fx2z_aa-2.0E0*1*I_ESP_G2x2y_Fx2z_a-2.0E0*1*I_ESP_G2y2z_Fx2z_a+1*I_ESP_D2y_Fx2z;
    abcd[iGrid*900+383] = 4.0E0*I_ESP_I2xy3z_Fx2z_aa-2.0E0*2*I_ESP_G2xyz_Fx2z_a-2.0E0*1*I_ESP_Gy3z_Fx2z_a+2*I_ESP_Dyz_Fx2z;
    abcd[iGrid*900+384] = 4.0E0*I_ESP_I2x4z_Fx2z_aa-2.0E0*3*I_ESP_G2x2z_Fx2z_a-2.0E0*1*I_ESP_G4z_Fx2z_a+3*I_ESP_D2z_Fx2z;
    abcd[iGrid*900+385] = 4.0E0*I_ESP_Ix4yz_Fx2z_aa;
    abcd[iGrid*900+386] = 4.0E0*I_ESP_Ix3y2z_Fx2z_aa-2.0E0*1*I_ESP_Gx3y_Fx2z_a;
    abcd[iGrid*900+387] = 4.0E0*I_ESP_Ix2y3z_Fx2z_aa-2.0E0*2*I_ESP_Gx2yz_Fx2z_a;
    abcd[iGrid*900+388] = 4.0E0*I_ESP_Ixy4z_Fx2z_aa-2.0E0*3*I_ESP_Gxy2z_Fx2z_a;
    abcd[iGrid*900+389] = 4.0E0*I_ESP_Ix5z_Fx2z_aa-2.0E0*4*I_ESP_Gx3z_Fx2z_a;
    abcd[iGrid*900+390] = 4.0E0*I_ESP_I5xz_F3y_aa-2.0E0*4*I_ESP_G3xz_F3y_a;
    abcd[iGrid*900+391] = 4.0E0*I_ESP_I4xyz_F3y_aa-2.0E0*3*I_ESP_G2xyz_F3y_a;
    abcd[iGrid*900+392] = 4.0E0*I_ESP_I4x2z_F3y_aa-2.0E0*1*I_ESP_G4x_F3y_a-2.0E0*3*I_ESP_G2x2z_F3y_a+3*1*I_ESP_D2x_F3y;
    abcd[iGrid*900+393] = 4.0E0*I_ESP_I3x2yz_F3y_aa-2.0E0*2*I_ESP_Gx2yz_F3y_a;
    abcd[iGrid*900+394] = 4.0E0*I_ESP_I3xy2z_F3y_aa-2.0E0*1*I_ESP_G3xy_F3y_a-2.0E0*2*I_ESP_Gxy2z_F3y_a+2*1*I_ESP_Dxy_F3y;
    abcd[iGrid*900+395] = 4.0E0*I_ESP_I3x3z_F3y_aa-2.0E0*2*I_ESP_G3xz_F3y_a-2.0E0*2*I_ESP_Gx3z_F3y_a+2*2*I_ESP_Dxz_F3y;
    abcd[iGrid*900+396] = 4.0E0*I_ESP_I2x3yz_F3y_aa-2.0E0*1*I_ESP_G3yz_F3y_a;
    abcd[iGrid*900+397] = 4.0E0*I_ESP_I2x2y2z_F3y_aa-2.0E0*1*I_ESP_G2x2y_F3y_a-2.0E0*1*I_ESP_G2y2z_F3y_a+1*I_ESP_D2y_F3y;
    abcd[iGrid*900+398] = 4.0E0*I_ESP_I2xy3z_F3y_aa-2.0E0*2*I_ESP_G2xyz_F3y_a-2.0E0*1*I_ESP_Gy3z_F3y_a+2*I_ESP_Dyz_F3y;
    abcd[iGrid*900+399] = 4.0E0*I_ESP_I2x4z_F3y_aa-2.0E0*3*I_ESP_G2x2z_F3y_a-2.0E0*1*I_ESP_G4z_F3y_a+3*I_ESP_D2z_F3y;
    abcd[iGrid*900+400] = 4.0E0*I_ESP_Ix4yz_F3y_aa;
    abcd[iGrid*900+401] = 4.0E0*I_ESP_Ix3y2z_F3y_aa-2.0E0*1*I_ESP_Gx3y_F3y_a;
    abcd[iGrid*900+402] = 4.0E0*I_ESP_Ix2y3z_F3y_aa-2.0E0*2*I_ESP_Gx2yz_F3y_a;
    abcd[iGrid*900+403] = 4.0E0*I_ESP_Ixy4z_F3y_aa-2.0E0*3*I_ESP_Gxy2z_F3y_a;
    abcd[iGrid*900+404] = 4.0E0*I_ESP_Ix5z_F3y_aa-2.0E0*4*I_ESP_Gx3z_F3y_a;
    abcd[iGrid*900+405] = 4.0E0*I_ESP_I5xz_F2yz_aa-2.0E0*4*I_ESP_G3xz_F2yz_a;
    abcd[iGrid*900+406] = 4.0E0*I_ESP_I4xyz_F2yz_aa-2.0E0*3*I_ESP_G2xyz_F2yz_a;
    abcd[iGrid*900+407] = 4.0E0*I_ESP_I4x2z_F2yz_aa-2.0E0*1*I_ESP_G4x_F2yz_a-2.0E0*3*I_ESP_G2x2z_F2yz_a+3*1*I_ESP_D2x_F2yz;
    abcd[iGrid*900+408] = 4.0E0*I_ESP_I3x2yz_F2yz_aa-2.0E0*2*I_ESP_Gx2yz_F2yz_a;
    abcd[iGrid*900+409] = 4.0E0*I_ESP_I3xy2z_F2yz_aa-2.0E0*1*I_ESP_G3xy_F2yz_a-2.0E0*2*I_ESP_Gxy2z_F2yz_a+2*1*I_ESP_Dxy_F2yz;
    abcd[iGrid*900+410] = 4.0E0*I_ESP_I3x3z_F2yz_aa-2.0E0*2*I_ESP_G3xz_F2yz_a-2.0E0*2*I_ESP_Gx3z_F2yz_a+2*2*I_ESP_Dxz_F2yz;
    abcd[iGrid*900+411] = 4.0E0*I_ESP_I2x3yz_F2yz_aa-2.0E0*1*I_ESP_G3yz_F2yz_a;
    abcd[iGrid*900+412] = 4.0E0*I_ESP_I2x2y2z_F2yz_aa-2.0E0*1*I_ESP_G2x2y_F2yz_a-2.0E0*1*I_ESP_G2y2z_F2yz_a+1*I_ESP_D2y_F2yz;
    abcd[iGrid*900+413] = 4.0E0*I_ESP_I2xy3z_F2yz_aa-2.0E0*2*I_ESP_G2xyz_F2yz_a-2.0E0*1*I_ESP_Gy3z_F2yz_a+2*I_ESP_Dyz_F2yz;
    abcd[iGrid*900+414] = 4.0E0*I_ESP_I2x4z_F2yz_aa-2.0E0*3*I_ESP_G2x2z_F2yz_a-2.0E0*1*I_ESP_G4z_F2yz_a+3*I_ESP_D2z_F2yz;
    abcd[iGrid*900+415] = 4.0E0*I_ESP_Ix4yz_F2yz_aa;
    abcd[iGrid*900+416] = 4.0E0*I_ESP_Ix3y2z_F2yz_aa-2.0E0*1*I_ESP_Gx3y_F2yz_a;
    abcd[iGrid*900+417] = 4.0E0*I_ESP_Ix2y3z_F2yz_aa-2.0E0*2*I_ESP_Gx2yz_F2yz_a;
    abcd[iGrid*900+418] = 4.0E0*I_ESP_Ixy4z_F2yz_aa-2.0E0*3*I_ESP_Gxy2z_F2yz_a;
    abcd[iGrid*900+419] = 4.0E0*I_ESP_Ix5z_F2yz_aa-2.0E0*4*I_ESP_Gx3z_F2yz_a;
    abcd[iGrid*900+420] = 4.0E0*I_ESP_I5xz_Fy2z_aa-2.0E0*4*I_ESP_G3xz_Fy2z_a;
    abcd[iGrid*900+421] = 4.0E0*I_ESP_I4xyz_Fy2z_aa-2.0E0*3*I_ESP_G2xyz_Fy2z_a;
    abcd[iGrid*900+422] = 4.0E0*I_ESP_I4x2z_Fy2z_aa-2.0E0*1*I_ESP_G4x_Fy2z_a-2.0E0*3*I_ESP_G2x2z_Fy2z_a+3*1*I_ESP_D2x_Fy2z;
    abcd[iGrid*900+423] = 4.0E0*I_ESP_I3x2yz_Fy2z_aa-2.0E0*2*I_ESP_Gx2yz_Fy2z_a;
    abcd[iGrid*900+424] = 4.0E0*I_ESP_I3xy2z_Fy2z_aa-2.0E0*1*I_ESP_G3xy_Fy2z_a-2.0E0*2*I_ESP_Gxy2z_Fy2z_a+2*1*I_ESP_Dxy_Fy2z;
    abcd[iGrid*900+425] = 4.0E0*I_ESP_I3x3z_Fy2z_aa-2.0E0*2*I_ESP_G3xz_Fy2z_a-2.0E0*2*I_ESP_Gx3z_Fy2z_a+2*2*I_ESP_Dxz_Fy2z;
    abcd[iGrid*900+426] = 4.0E0*I_ESP_I2x3yz_Fy2z_aa-2.0E0*1*I_ESP_G3yz_Fy2z_a;
    abcd[iGrid*900+427] = 4.0E0*I_ESP_I2x2y2z_Fy2z_aa-2.0E0*1*I_ESP_G2x2y_Fy2z_a-2.0E0*1*I_ESP_G2y2z_Fy2z_a+1*I_ESP_D2y_Fy2z;
    abcd[iGrid*900+428] = 4.0E0*I_ESP_I2xy3z_Fy2z_aa-2.0E0*2*I_ESP_G2xyz_Fy2z_a-2.0E0*1*I_ESP_Gy3z_Fy2z_a+2*I_ESP_Dyz_Fy2z;
    abcd[iGrid*900+429] = 4.0E0*I_ESP_I2x4z_Fy2z_aa-2.0E0*3*I_ESP_G2x2z_Fy2z_a-2.0E0*1*I_ESP_G4z_Fy2z_a+3*I_ESP_D2z_Fy2z;
    abcd[iGrid*900+430] = 4.0E0*I_ESP_Ix4yz_Fy2z_aa;
    abcd[iGrid*900+431] = 4.0E0*I_ESP_Ix3y2z_Fy2z_aa-2.0E0*1*I_ESP_Gx3y_Fy2z_a;
    abcd[iGrid*900+432] = 4.0E0*I_ESP_Ix2y3z_Fy2z_aa-2.0E0*2*I_ESP_Gx2yz_Fy2z_a;
    abcd[iGrid*900+433] = 4.0E0*I_ESP_Ixy4z_Fy2z_aa-2.0E0*3*I_ESP_Gxy2z_Fy2z_a;
    abcd[iGrid*900+434] = 4.0E0*I_ESP_Ix5z_Fy2z_aa-2.0E0*4*I_ESP_Gx3z_Fy2z_a;
    abcd[iGrid*900+435] = 4.0E0*I_ESP_I5xz_F3z_aa-2.0E0*4*I_ESP_G3xz_F3z_a;
    abcd[iGrid*900+436] = 4.0E0*I_ESP_I4xyz_F3z_aa-2.0E0*3*I_ESP_G2xyz_F3z_a;
    abcd[iGrid*900+437] = 4.0E0*I_ESP_I4x2z_F3z_aa-2.0E0*1*I_ESP_G4x_F3z_a-2.0E0*3*I_ESP_G2x2z_F3z_a+3*1*I_ESP_D2x_F3z;
    abcd[iGrid*900+438] = 4.0E0*I_ESP_I3x2yz_F3z_aa-2.0E0*2*I_ESP_Gx2yz_F3z_a;
    abcd[iGrid*900+439] = 4.0E0*I_ESP_I3xy2z_F3z_aa-2.0E0*1*I_ESP_G3xy_F3z_a-2.0E0*2*I_ESP_Gxy2z_F3z_a+2*1*I_ESP_Dxy_F3z;
    abcd[iGrid*900+440] = 4.0E0*I_ESP_I3x3z_F3z_aa-2.0E0*2*I_ESP_G3xz_F3z_a-2.0E0*2*I_ESP_Gx3z_F3z_a+2*2*I_ESP_Dxz_F3z;
    abcd[iGrid*900+441] = 4.0E0*I_ESP_I2x3yz_F3z_aa-2.0E0*1*I_ESP_G3yz_F3z_a;
    abcd[iGrid*900+442] = 4.0E0*I_ESP_I2x2y2z_F3z_aa-2.0E0*1*I_ESP_G2x2y_F3z_a-2.0E0*1*I_ESP_G2y2z_F3z_a+1*I_ESP_D2y_F3z;
    abcd[iGrid*900+443] = 4.0E0*I_ESP_I2xy3z_F3z_aa-2.0E0*2*I_ESP_G2xyz_F3z_a-2.0E0*1*I_ESP_Gy3z_F3z_a+2*I_ESP_Dyz_F3z;
    abcd[iGrid*900+444] = 4.0E0*I_ESP_I2x4z_F3z_aa-2.0E0*3*I_ESP_G2x2z_F3z_a-2.0E0*1*I_ESP_G4z_F3z_a+3*I_ESP_D2z_F3z;
    abcd[iGrid*900+445] = 4.0E0*I_ESP_Ix4yz_F3z_aa;
    abcd[iGrid*900+446] = 4.0E0*I_ESP_Ix3y2z_F3z_aa-2.0E0*1*I_ESP_Gx3y_F3z_a;
    abcd[iGrid*900+447] = 4.0E0*I_ESP_Ix2y3z_F3z_aa-2.0E0*2*I_ESP_Gx2yz_F3z_a;
    abcd[iGrid*900+448] = 4.0E0*I_ESP_Ixy4z_F3z_aa-2.0E0*3*I_ESP_Gxy2z_F3z_a;
    abcd[iGrid*900+449] = 4.0E0*I_ESP_Ix5z_F3z_aa-2.0E0*4*I_ESP_Gx3z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_aa
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_D_F
     ************************************************************/
    abcd[iGrid*900+450] = 4.0E0*I_ESP_I4x2y_F3x_aa-2.0E0*1*I_ESP_G4x_F3x_a;
    abcd[iGrid*900+451] = 4.0E0*I_ESP_I3x3y_F3x_aa-2.0E0*1*I_ESP_G3xy_F3x_a-2.0E0*2*I_ESP_G3xy_F3x_a;
    abcd[iGrid*900+452] = 4.0E0*I_ESP_I3x2yz_F3x_aa-2.0E0*1*I_ESP_G3xz_F3x_a;
    abcd[iGrid*900+453] = 4.0E0*I_ESP_I2x4y_F3x_aa-2.0E0*2*I_ESP_G2x2y_F3x_a-2.0E0*3*I_ESP_G2x2y_F3x_a+2*1*I_ESP_D2x_F3x;
    abcd[iGrid*900+454] = 4.0E0*I_ESP_I2x3yz_F3x_aa-2.0E0*1*I_ESP_G2xyz_F3x_a-2.0E0*2*I_ESP_G2xyz_F3x_a;
    abcd[iGrid*900+455] = 4.0E0*I_ESP_I2x2y2z_F3x_aa-2.0E0*1*I_ESP_G2x2z_F3x_a;
    abcd[iGrid*900+456] = 4.0E0*I_ESP_Ix5y_F3x_aa-2.0E0*3*I_ESP_Gx3y_F3x_a-2.0E0*4*I_ESP_Gx3y_F3x_a+3*2*I_ESP_Dxy_F3x;
    abcd[iGrid*900+457] = 4.0E0*I_ESP_Ix4yz_F3x_aa-2.0E0*2*I_ESP_Gx2yz_F3x_a-2.0E0*3*I_ESP_Gx2yz_F3x_a+2*1*I_ESP_Dxz_F3x;
    abcd[iGrid*900+458] = 4.0E0*I_ESP_Ix3y2z_F3x_aa-2.0E0*1*I_ESP_Gxy2z_F3x_a-2.0E0*2*I_ESP_Gxy2z_F3x_a;
    abcd[iGrid*900+459] = 4.0E0*I_ESP_Ix2y3z_F3x_aa-2.0E0*1*I_ESP_Gx3z_F3x_a;
    abcd[iGrid*900+460] = 4.0E0*I_ESP_I6y_F3x_aa-2.0E0*4*I_ESP_G4y_F3x_a-2.0E0*5*I_ESP_G4y_F3x_a+4*3*I_ESP_D2y_F3x;
    abcd[iGrid*900+461] = 4.0E0*I_ESP_I5yz_F3x_aa-2.0E0*3*I_ESP_G3yz_F3x_a-2.0E0*4*I_ESP_G3yz_F3x_a+3*2*I_ESP_Dyz_F3x;
    abcd[iGrid*900+462] = 4.0E0*I_ESP_I4y2z_F3x_aa-2.0E0*2*I_ESP_G2y2z_F3x_a-2.0E0*3*I_ESP_G2y2z_F3x_a+2*1*I_ESP_D2z_F3x;
    abcd[iGrid*900+463] = 4.0E0*I_ESP_I3y3z_F3x_aa-2.0E0*1*I_ESP_Gy3z_F3x_a-2.0E0*2*I_ESP_Gy3z_F3x_a;
    abcd[iGrid*900+464] = 4.0E0*I_ESP_I2y4z_F3x_aa-2.0E0*1*I_ESP_G4z_F3x_a;
    abcd[iGrid*900+465] = 4.0E0*I_ESP_I4x2y_F2xy_aa-2.0E0*1*I_ESP_G4x_F2xy_a;
    abcd[iGrid*900+466] = 4.0E0*I_ESP_I3x3y_F2xy_aa-2.0E0*1*I_ESP_G3xy_F2xy_a-2.0E0*2*I_ESP_G3xy_F2xy_a;
    abcd[iGrid*900+467] = 4.0E0*I_ESP_I3x2yz_F2xy_aa-2.0E0*1*I_ESP_G3xz_F2xy_a;
    abcd[iGrid*900+468] = 4.0E0*I_ESP_I2x4y_F2xy_aa-2.0E0*2*I_ESP_G2x2y_F2xy_a-2.0E0*3*I_ESP_G2x2y_F2xy_a+2*1*I_ESP_D2x_F2xy;
    abcd[iGrid*900+469] = 4.0E0*I_ESP_I2x3yz_F2xy_aa-2.0E0*1*I_ESP_G2xyz_F2xy_a-2.0E0*2*I_ESP_G2xyz_F2xy_a;
    abcd[iGrid*900+470] = 4.0E0*I_ESP_I2x2y2z_F2xy_aa-2.0E0*1*I_ESP_G2x2z_F2xy_a;
    abcd[iGrid*900+471] = 4.0E0*I_ESP_Ix5y_F2xy_aa-2.0E0*3*I_ESP_Gx3y_F2xy_a-2.0E0*4*I_ESP_Gx3y_F2xy_a+3*2*I_ESP_Dxy_F2xy;
    abcd[iGrid*900+472] = 4.0E0*I_ESP_Ix4yz_F2xy_aa-2.0E0*2*I_ESP_Gx2yz_F2xy_a-2.0E0*3*I_ESP_Gx2yz_F2xy_a+2*1*I_ESP_Dxz_F2xy;
    abcd[iGrid*900+473] = 4.0E0*I_ESP_Ix3y2z_F2xy_aa-2.0E0*1*I_ESP_Gxy2z_F2xy_a-2.0E0*2*I_ESP_Gxy2z_F2xy_a;
    abcd[iGrid*900+474] = 4.0E0*I_ESP_Ix2y3z_F2xy_aa-2.0E0*1*I_ESP_Gx3z_F2xy_a;
    abcd[iGrid*900+475] = 4.0E0*I_ESP_I6y_F2xy_aa-2.0E0*4*I_ESP_G4y_F2xy_a-2.0E0*5*I_ESP_G4y_F2xy_a+4*3*I_ESP_D2y_F2xy;
    abcd[iGrid*900+476] = 4.0E0*I_ESP_I5yz_F2xy_aa-2.0E0*3*I_ESP_G3yz_F2xy_a-2.0E0*4*I_ESP_G3yz_F2xy_a+3*2*I_ESP_Dyz_F2xy;
    abcd[iGrid*900+477] = 4.0E0*I_ESP_I4y2z_F2xy_aa-2.0E0*2*I_ESP_G2y2z_F2xy_a-2.0E0*3*I_ESP_G2y2z_F2xy_a+2*1*I_ESP_D2z_F2xy;
    abcd[iGrid*900+478] = 4.0E0*I_ESP_I3y3z_F2xy_aa-2.0E0*1*I_ESP_Gy3z_F2xy_a-2.0E0*2*I_ESP_Gy3z_F2xy_a;
    abcd[iGrid*900+479] = 4.0E0*I_ESP_I2y4z_F2xy_aa-2.0E0*1*I_ESP_G4z_F2xy_a;
    abcd[iGrid*900+480] = 4.0E0*I_ESP_I4x2y_F2xz_aa-2.0E0*1*I_ESP_G4x_F2xz_a;
    abcd[iGrid*900+481] = 4.0E0*I_ESP_I3x3y_F2xz_aa-2.0E0*1*I_ESP_G3xy_F2xz_a-2.0E0*2*I_ESP_G3xy_F2xz_a;
    abcd[iGrid*900+482] = 4.0E0*I_ESP_I3x2yz_F2xz_aa-2.0E0*1*I_ESP_G3xz_F2xz_a;
    abcd[iGrid*900+483] = 4.0E0*I_ESP_I2x4y_F2xz_aa-2.0E0*2*I_ESP_G2x2y_F2xz_a-2.0E0*3*I_ESP_G2x2y_F2xz_a+2*1*I_ESP_D2x_F2xz;
    abcd[iGrid*900+484] = 4.0E0*I_ESP_I2x3yz_F2xz_aa-2.0E0*1*I_ESP_G2xyz_F2xz_a-2.0E0*2*I_ESP_G2xyz_F2xz_a;
    abcd[iGrid*900+485] = 4.0E0*I_ESP_I2x2y2z_F2xz_aa-2.0E0*1*I_ESP_G2x2z_F2xz_a;
    abcd[iGrid*900+486] = 4.0E0*I_ESP_Ix5y_F2xz_aa-2.0E0*3*I_ESP_Gx3y_F2xz_a-2.0E0*4*I_ESP_Gx3y_F2xz_a+3*2*I_ESP_Dxy_F2xz;
    abcd[iGrid*900+487] = 4.0E0*I_ESP_Ix4yz_F2xz_aa-2.0E0*2*I_ESP_Gx2yz_F2xz_a-2.0E0*3*I_ESP_Gx2yz_F2xz_a+2*1*I_ESP_Dxz_F2xz;
    abcd[iGrid*900+488] = 4.0E0*I_ESP_Ix3y2z_F2xz_aa-2.0E0*1*I_ESP_Gxy2z_F2xz_a-2.0E0*2*I_ESP_Gxy2z_F2xz_a;
    abcd[iGrid*900+489] = 4.0E0*I_ESP_Ix2y3z_F2xz_aa-2.0E0*1*I_ESP_Gx3z_F2xz_a;
    abcd[iGrid*900+490] = 4.0E0*I_ESP_I6y_F2xz_aa-2.0E0*4*I_ESP_G4y_F2xz_a-2.0E0*5*I_ESP_G4y_F2xz_a+4*3*I_ESP_D2y_F2xz;
    abcd[iGrid*900+491] = 4.0E0*I_ESP_I5yz_F2xz_aa-2.0E0*3*I_ESP_G3yz_F2xz_a-2.0E0*4*I_ESP_G3yz_F2xz_a+3*2*I_ESP_Dyz_F2xz;
    abcd[iGrid*900+492] = 4.0E0*I_ESP_I4y2z_F2xz_aa-2.0E0*2*I_ESP_G2y2z_F2xz_a-2.0E0*3*I_ESP_G2y2z_F2xz_a+2*1*I_ESP_D2z_F2xz;
    abcd[iGrid*900+493] = 4.0E0*I_ESP_I3y3z_F2xz_aa-2.0E0*1*I_ESP_Gy3z_F2xz_a-2.0E0*2*I_ESP_Gy3z_F2xz_a;
    abcd[iGrid*900+494] = 4.0E0*I_ESP_I2y4z_F2xz_aa-2.0E0*1*I_ESP_G4z_F2xz_a;
    abcd[iGrid*900+495] = 4.0E0*I_ESP_I4x2y_Fx2y_aa-2.0E0*1*I_ESP_G4x_Fx2y_a;
    abcd[iGrid*900+496] = 4.0E0*I_ESP_I3x3y_Fx2y_aa-2.0E0*1*I_ESP_G3xy_Fx2y_a-2.0E0*2*I_ESP_G3xy_Fx2y_a;
    abcd[iGrid*900+497] = 4.0E0*I_ESP_I3x2yz_Fx2y_aa-2.0E0*1*I_ESP_G3xz_Fx2y_a;
    abcd[iGrid*900+498] = 4.0E0*I_ESP_I2x4y_Fx2y_aa-2.0E0*2*I_ESP_G2x2y_Fx2y_a-2.0E0*3*I_ESP_G2x2y_Fx2y_a+2*1*I_ESP_D2x_Fx2y;
    abcd[iGrid*900+499] = 4.0E0*I_ESP_I2x3yz_Fx2y_aa-2.0E0*1*I_ESP_G2xyz_Fx2y_a-2.0E0*2*I_ESP_G2xyz_Fx2y_a;
    abcd[iGrid*900+500] = 4.0E0*I_ESP_I2x2y2z_Fx2y_aa-2.0E0*1*I_ESP_G2x2z_Fx2y_a;
    abcd[iGrid*900+501] = 4.0E0*I_ESP_Ix5y_Fx2y_aa-2.0E0*3*I_ESP_Gx3y_Fx2y_a-2.0E0*4*I_ESP_Gx3y_Fx2y_a+3*2*I_ESP_Dxy_Fx2y;
    abcd[iGrid*900+502] = 4.0E0*I_ESP_Ix4yz_Fx2y_aa-2.0E0*2*I_ESP_Gx2yz_Fx2y_a-2.0E0*3*I_ESP_Gx2yz_Fx2y_a+2*1*I_ESP_Dxz_Fx2y;
    abcd[iGrid*900+503] = 4.0E0*I_ESP_Ix3y2z_Fx2y_aa-2.0E0*1*I_ESP_Gxy2z_Fx2y_a-2.0E0*2*I_ESP_Gxy2z_Fx2y_a;
    abcd[iGrid*900+504] = 4.0E0*I_ESP_Ix2y3z_Fx2y_aa-2.0E0*1*I_ESP_Gx3z_Fx2y_a;
    abcd[iGrid*900+505] = 4.0E0*I_ESP_I6y_Fx2y_aa-2.0E0*4*I_ESP_G4y_Fx2y_a-2.0E0*5*I_ESP_G4y_Fx2y_a+4*3*I_ESP_D2y_Fx2y;
    abcd[iGrid*900+506] = 4.0E0*I_ESP_I5yz_Fx2y_aa-2.0E0*3*I_ESP_G3yz_Fx2y_a-2.0E0*4*I_ESP_G3yz_Fx2y_a+3*2*I_ESP_Dyz_Fx2y;
    abcd[iGrid*900+507] = 4.0E0*I_ESP_I4y2z_Fx2y_aa-2.0E0*2*I_ESP_G2y2z_Fx2y_a-2.0E0*3*I_ESP_G2y2z_Fx2y_a+2*1*I_ESP_D2z_Fx2y;
    abcd[iGrid*900+508] = 4.0E0*I_ESP_I3y3z_Fx2y_aa-2.0E0*1*I_ESP_Gy3z_Fx2y_a-2.0E0*2*I_ESP_Gy3z_Fx2y_a;
    abcd[iGrid*900+509] = 4.0E0*I_ESP_I2y4z_Fx2y_aa-2.0E0*1*I_ESP_G4z_Fx2y_a;
    abcd[iGrid*900+510] = 4.0E0*I_ESP_I4x2y_Fxyz_aa-2.0E0*1*I_ESP_G4x_Fxyz_a;
    abcd[iGrid*900+511] = 4.0E0*I_ESP_I3x3y_Fxyz_aa-2.0E0*1*I_ESP_G3xy_Fxyz_a-2.0E0*2*I_ESP_G3xy_Fxyz_a;
    abcd[iGrid*900+512] = 4.0E0*I_ESP_I3x2yz_Fxyz_aa-2.0E0*1*I_ESP_G3xz_Fxyz_a;
    abcd[iGrid*900+513] = 4.0E0*I_ESP_I2x4y_Fxyz_aa-2.0E0*2*I_ESP_G2x2y_Fxyz_a-2.0E0*3*I_ESP_G2x2y_Fxyz_a+2*1*I_ESP_D2x_Fxyz;
    abcd[iGrid*900+514] = 4.0E0*I_ESP_I2x3yz_Fxyz_aa-2.0E0*1*I_ESP_G2xyz_Fxyz_a-2.0E0*2*I_ESP_G2xyz_Fxyz_a;
    abcd[iGrid*900+515] = 4.0E0*I_ESP_I2x2y2z_Fxyz_aa-2.0E0*1*I_ESP_G2x2z_Fxyz_a;
    abcd[iGrid*900+516] = 4.0E0*I_ESP_Ix5y_Fxyz_aa-2.0E0*3*I_ESP_Gx3y_Fxyz_a-2.0E0*4*I_ESP_Gx3y_Fxyz_a+3*2*I_ESP_Dxy_Fxyz;
    abcd[iGrid*900+517] = 4.0E0*I_ESP_Ix4yz_Fxyz_aa-2.0E0*2*I_ESP_Gx2yz_Fxyz_a-2.0E0*3*I_ESP_Gx2yz_Fxyz_a+2*1*I_ESP_Dxz_Fxyz;
    abcd[iGrid*900+518] = 4.0E0*I_ESP_Ix3y2z_Fxyz_aa-2.0E0*1*I_ESP_Gxy2z_Fxyz_a-2.0E0*2*I_ESP_Gxy2z_Fxyz_a;
    abcd[iGrid*900+519] = 4.0E0*I_ESP_Ix2y3z_Fxyz_aa-2.0E0*1*I_ESP_Gx3z_Fxyz_a;
    abcd[iGrid*900+520] = 4.0E0*I_ESP_I6y_Fxyz_aa-2.0E0*4*I_ESP_G4y_Fxyz_a-2.0E0*5*I_ESP_G4y_Fxyz_a+4*3*I_ESP_D2y_Fxyz;
    abcd[iGrid*900+521] = 4.0E0*I_ESP_I5yz_Fxyz_aa-2.0E0*3*I_ESP_G3yz_Fxyz_a-2.0E0*4*I_ESP_G3yz_Fxyz_a+3*2*I_ESP_Dyz_Fxyz;
    abcd[iGrid*900+522] = 4.0E0*I_ESP_I4y2z_Fxyz_aa-2.0E0*2*I_ESP_G2y2z_Fxyz_a-2.0E0*3*I_ESP_G2y2z_Fxyz_a+2*1*I_ESP_D2z_Fxyz;
    abcd[iGrid*900+523] = 4.0E0*I_ESP_I3y3z_Fxyz_aa-2.0E0*1*I_ESP_Gy3z_Fxyz_a-2.0E0*2*I_ESP_Gy3z_Fxyz_a;
    abcd[iGrid*900+524] = 4.0E0*I_ESP_I2y4z_Fxyz_aa-2.0E0*1*I_ESP_G4z_Fxyz_a;
    abcd[iGrid*900+525] = 4.0E0*I_ESP_I4x2y_Fx2z_aa-2.0E0*1*I_ESP_G4x_Fx2z_a;
    abcd[iGrid*900+526] = 4.0E0*I_ESP_I3x3y_Fx2z_aa-2.0E0*1*I_ESP_G3xy_Fx2z_a-2.0E0*2*I_ESP_G3xy_Fx2z_a;
    abcd[iGrid*900+527] = 4.0E0*I_ESP_I3x2yz_Fx2z_aa-2.0E0*1*I_ESP_G3xz_Fx2z_a;
    abcd[iGrid*900+528] = 4.0E0*I_ESP_I2x4y_Fx2z_aa-2.0E0*2*I_ESP_G2x2y_Fx2z_a-2.0E0*3*I_ESP_G2x2y_Fx2z_a+2*1*I_ESP_D2x_Fx2z;
    abcd[iGrid*900+529] = 4.0E0*I_ESP_I2x3yz_Fx2z_aa-2.0E0*1*I_ESP_G2xyz_Fx2z_a-2.0E0*2*I_ESP_G2xyz_Fx2z_a;
    abcd[iGrid*900+530] = 4.0E0*I_ESP_I2x2y2z_Fx2z_aa-2.0E0*1*I_ESP_G2x2z_Fx2z_a;
    abcd[iGrid*900+531] = 4.0E0*I_ESP_Ix5y_Fx2z_aa-2.0E0*3*I_ESP_Gx3y_Fx2z_a-2.0E0*4*I_ESP_Gx3y_Fx2z_a+3*2*I_ESP_Dxy_Fx2z;
    abcd[iGrid*900+532] = 4.0E0*I_ESP_Ix4yz_Fx2z_aa-2.0E0*2*I_ESP_Gx2yz_Fx2z_a-2.0E0*3*I_ESP_Gx2yz_Fx2z_a+2*1*I_ESP_Dxz_Fx2z;
    abcd[iGrid*900+533] = 4.0E0*I_ESP_Ix3y2z_Fx2z_aa-2.0E0*1*I_ESP_Gxy2z_Fx2z_a-2.0E0*2*I_ESP_Gxy2z_Fx2z_a;
    abcd[iGrid*900+534] = 4.0E0*I_ESP_Ix2y3z_Fx2z_aa-2.0E0*1*I_ESP_Gx3z_Fx2z_a;
    abcd[iGrid*900+535] = 4.0E0*I_ESP_I6y_Fx2z_aa-2.0E0*4*I_ESP_G4y_Fx2z_a-2.0E0*5*I_ESP_G4y_Fx2z_a+4*3*I_ESP_D2y_Fx2z;
    abcd[iGrid*900+536] = 4.0E0*I_ESP_I5yz_Fx2z_aa-2.0E0*3*I_ESP_G3yz_Fx2z_a-2.0E0*4*I_ESP_G3yz_Fx2z_a+3*2*I_ESP_Dyz_Fx2z;
    abcd[iGrid*900+537] = 4.0E0*I_ESP_I4y2z_Fx2z_aa-2.0E0*2*I_ESP_G2y2z_Fx2z_a-2.0E0*3*I_ESP_G2y2z_Fx2z_a+2*1*I_ESP_D2z_Fx2z;
    abcd[iGrid*900+538] = 4.0E0*I_ESP_I3y3z_Fx2z_aa-2.0E0*1*I_ESP_Gy3z_Fx2z_a-2.0E0*2*I_ESP_Gy3z_Fx2z_a;
    abcd[iGrid*900+539] = 4.0E0*I_ESP_I2y4z_Fx2z_aa-2.0E0*1*I_ESP_G4z_Fx2z_a;
    abcd[iGrid*900+540] = 4.0E0*I_ESP_I4x2y_F3y_aa-2.0E0*1*I_ESP_G4x_F3y_a;
    abcd[iGrid*900+541] = 4.0E0*I_ESP_I3x3y_F3y_aa-2.0E0*1*I_ESP_G3xy_F3y_a-2.0E0*2*I_ESP_G3xy_F3y_a;
    abcd[iGrid*900+542] = 4.0E0*I_ESP_I3x2yz_F3y_aa-2.0E0*1*I_ESP_G3xz_F3y_a;
    abcd[iGrid*900+543] = 4.0E0*I_ESP_I2x4y_F3y_aa-2.0E0*2*I_ESP_G2x2y_F3y_a-2.0E0*3*I_ESP_G2x2y_F3y_a+2*1*I_ESP_D2x_F3y;
    abcd[iGrid*900+544] = 4.0E0*I_ESP_I2x3yz_F3y_aa-2.0E0*1*I_ESP_G2xyz_F3y_a-2.0E0*2*I_ESP_G2xyz_F3y_a;
    abcd[iGrid*900+545] = 4.0E0*I_ESP_I2x2y2z_F3y_aa-2.0E0*1*I_ESP_G2x2z_F3y_a;
    abcd[iGrid*900+546] = 4.0E0*I_ESP_Ix5y_F3y_aa-2.0E0*3*I_ESP_Gx3y_F3y_a-2.0E0*4*I_ESP_Gx3y_F3y_a+3*2*I_ESP_Dxy_F3y;
    abcd[iGrid*900+547] = 4.0E0*I_ESP_Ix4yz_F3y_aa-2.0E0*2*I_ESP_Gx2yz_F3y_a-2.0E0*3*I_ESP_Gx2yz_F3y_a+2*1*I_ESP_Dxz_F3y;
    abcd[iGrid*900+548] = 4.0E0*I_ESP_Ix3y2z_F3y_aa-2.0E0*1*I_ESP_Gxy2z_F3y_a-2.0E0*2*I_ESP_Gxy2z_F3y_a;
    abcd[iGrid*900+549] = 4.0E0*I_ESP_Ix2y3z_F3y_aa-2.0E0*1*I_ESP_Gx3z_F3y_a;
    abcd[iGrid*900+550] = 4.0E0*I_ESP_I6y_F3y_aa-2.0E0*4*I_ESP_G4y_F3y_a-2.0E0*5*I_ESP_G4y_F3y_a+4*3*I_ESP_D2y_F3y;
    abcd[iGrid*900+551] = 4.0E0*I_ESP_I5yz_F3y_aa-2.0E0*3*I_ESP_G3yz_F3y_a-2.0E0*4*I_ESP_G3yz_F3y_a+3*2*I_ESP_Dyz_F3y;
    abcd[iGrid*900+552] = 4.0E0*I_ESP_I4y2z_F3y_aa-2.0E0*2*I_ESP_G2y2z_F3y_a-2.0E0*3*I_ESP_G2y2z_F3y_a+2*1*I_ESP_D2z_F3y;
    abcd[iGrid*900+553] = 4.0E0*I_ESP_I3y3z_F3y_aa-2.0E0*1*I_ESP_Gy3z_F3y_a-2.0E0*2*I_ESP_Gy3z_F3y_a;
    abcd[iGrid*900+554] = 4.0E0*I_ESP_I2y4z_F3y_aa-2.0E0*1*I_ESP_G4z_F3y_a;
    abcd[iGrid*900+555] = 4.0E0*I_ESP_I4x2y_F2yz_aa-2.0E0*1*I_ESP_G4x_F2yz_a;
    abcd[iGrid*900+556] = 4.0E0*I_ESP_I3x3y_F2yz_aa-2.0E0*1*I_ESP_G3xy_F2yz_a-2.0E0*2*I_ESP_G3xy_F2yz_a;
    abcd[iGrid*900+557] = 4.0E0*I_ESP_I3x2yz_F2yz_aa-2.0E0*1*I_ESP_G3xz_F2yz_a;
    abcd[iGrid*900+558] = 4.0E0*I_ESP_I2x4y_F2yz_aa-2.0E0*2*I_ESP_G2x2y_F2yz_a-2.0E0*3*I_ESP_G2x2y_F2yz_a+2*1*I_ESP_D2x_F2yz;
    abcd[iGrid*900+559] = 4.0E0*I_ESP_I2x3yz_F2yz_aa-2.0E0*1*I_ESP_G2xyz_F2yz_a-2.0E0*2*I_ESP_G2xyz_F2yz_a;
    abcd[iGrid*900+560] = 4.0E0*I_ESP_I2x2y2z_F2yz_aa-2.0E0*1*I_ESP_G2x2z_F2yz_a;
    abcd[iGrid*900+561] = 4.0E0*I_ESP_Ix5y_F2yz_aa-2.0E0*3*I_ESP_Gx3y_F2yz_a-2.0E0*4*I_ESP_Gx3y_F2yz_a+3*2*I_ESP_Dxy_F2yz;
    abcd[iGrid*900+562] = 4.0E0*I_ESP_Ix4yz_F2yz_aa-2.0E0*2*I_ESP_Gx2yz_F2yz_a-2.0E0*3*I_ESP_Gx2yz_F2yz_a+2*1*I_ESP_Dxz_F2yz;
    abcd[iGrid*900+563] = 4.0E0*I_ESP_Ix3y2z_F2yz_aa-2.0E0*1*I_ESP_Gxy2z_F2yz_a-2.0E0*2*I_ESP_Gxy2z_F2yz_a;
    abcd[iGrid*900+564] = 4.0E0*I_ESP_Ix2y3z_F2yz_aa-2.0E0*1*I_ESP_Gx3z_F2yz_a;
    abcd[iGrid*900+565] = 4.0E0*I_ESP_I6y_F2yz_aa-2.0E0*4*I_ESP_G4y_F2yz_a-2.0E0*5*I_ESP_G4y_F2yz_a+4*3*I_ESP_D2y_F2yz;
    abcd[iGrid*900+566] = 4.0E0*I_ESP_I5yz_F2yz_aa-2.0E0*3*I_ESP_G3yz_F2yz_a-2.0E0*4*I_ESP_G3yz_F2yz_a+3*2*I_ESP_Dyz_F2yz;
    abcd[iGrid*900+567] = 4.0E0*I_ESP_I4y2z_F2yz_aa-2.0E0*2*I_ESP_G2y2z_F2yz_a-2.0E0*3*I_ESP_G2y2z_F2yz_a+2*1*I_ESP_D2z_F2yz;
    abcd[iGrid*900+568] = 4.0E0*I_ESP_I3y3z_F2yz_aa-2.0E0*1*I_ESP_Gy3z_F2yz_a-2.0E0*2*I_ESP_Gy3z_F2yz_a;
    abcd[iGrid*900+569] = 4.0E0*I_ESP_I2y4z_F2yz_aa-2.0E0*1*I_ESP_G4z_F2yz_a;
    abcd[iGrid*900+570] = 4.0E0*I_ESP_I4x2y_Fy2z_aa-2.0E0*1*I_ESP_G4x_Fy2z_a;
    abcd[iGrid*900+571] = 4.0E0*I_ESP_I3x3y_Fy2z_aa-2.0E0*1*I_ESP_G3xy_Fy2z_a-2.0E0*2*I_ESP_G3xy_Fy2z_a;
    abcd[iGrid*900+572] = 4.0E0*I_ESP_I3x2yz_Fy2z_aa-2.0E0*1*I_ESP_G3xz_Fy2z_a;
    abcd[iGrid*900+573] = 4.0E0*I_ESP_I2x4y_Fy2z_aa-2.0E0*2*I_ESP_G2x2y_Fy2z_a-2.0E0*3*I_ESP_G2x2y_Fy2z_a+2*1*I_ESP_D2x_Fy2z;
    abcd[iGrid*900+574] = 4.0E0*I_ESP_I2x3yz_Fy2z_aa-2.0E0*1*I_ESP_G2xyz_Fy2z_a-2.0E0*2*I_ESP_G2xyz_Fy2z_a;
    abcd[iGrid*900+575] = 4.0E0*I_ESP_I2x2y2z_Fy2z_aa-2.0E0*1*I_ESP_G2x2z_Fy2z_a;
    abcd[iGrid*900+576] = 4.0E0*I_ESP_Ix5y_Fy2z_aa-2.0E0*3*I_ESP_Gx3y_Fy2z_a-2.0E0*4*I_ESP_Gx3y_Fy2z_a+3*2*I_ESP_Dxy_Fy2z;
    abcd[iGrid*900+577] = 4.0E0*I_ESP_Ix4yz_Fy2z_aa-2.0E0*2*I_ESP_Gx2yz_Fy2z_a-2.0E0*3*I_ESP_Gx2yz_Fy2z_a+2*1*I_ESP_Dxz_Fy2z;
    abcd[iGrid*900+578] = 4.0E0*I_ESP_Ix3y2z_Fy2z_aa-2.0E0*1*I_ESP_Gxy2z_Fy2z_a-2.0E0*2*I_ESP_Gxy2z_Fy2z_a;
    abcd[iGrid*900+579] = 4.0E0*I_ESP_Ix2y3z_Fy2z_aa-2.0E0*1*I_ESP_Gx3z_Fy2z_a;
    abcd[iGrid*900+580] = 4.0E0*I_ESP_I6y_Fy2z_aa-2.0E0*4*I_ESP_G4y_Fy2z_a-2.0E0*5*I_ESP_G4y_Fy2z_a+4*3*I_ESP_D2y_Fy2z;
    abcd[iGrid*900+581] = 4.0E0*I_ESP_I5yz_Fy2z_aa-2.0E0*3*I_ESP_G3yz_Fy2z_a-2.0E0*4*I_ESP_G3yz_Fy2z_a+3*2*I_ESP_Dyz_Fy2z;
    abcd[iGrid*900+582] = 4.0E0*I_ESP_I4y2z_Fy2z_aa-2.0E0*2*I_ESP_G2y2z_Fy2z_a-2.0E0*3*I_ESP_G2y2z_Fy2z_a+2*1*I_ESP_D2z_Fy2z;
    abcd[iGrid*900+583] = 4.0E0*I_ESP_I3y3z_Fy2z_aa-2.0E0*1*I_ESP_Gy3z_Fy2z_a-2.0E0*2*I_ESP_Gy3z_Fy2z_a;
    abcd[iGrid*900+584] = 4.0E0*I_ESP_I2y4z_Fy2z_aa-2.0E0*1*I_ESP_G4z_Fy2z_a;
    abcd[iGrid*900+585] = 4.0E0*I_ESP_I4x2y_F3z_aa-2.0E0*1*I_ESP_G4x_F3z_a;
    abcd[iGrid*900+586] = 4.0E0*I_ESP_I3x3y_F3z_aa-2.0E0*1*I_ESP_G3xy_F3z_a-2.0E0*2*I_ESP_G3xy_F3z_a;
    abcd[iGrid*900+587] = 4.0E0*I_ESP_I3x2yz_F3z_aa-2.0E0*1*I_ESP_G3xz_F3z_a;
    abcd[iGrid*900+588] = 4.0E0*I_ESP_I2x4y_F3z_aa-2.0E0*2*I_ESP_G2x2y_F3z_a-2.0E0*3*I_ESP_G2x2y_F3z_a+2*1*I_ESP_D2x_F3z;
    abcd[iGrid*900+589] = 4.0E0*I_ESP_I2x3yz_F3z_aa-2.0E0*1*I_ESP_G2xyz_F3z_a-2.0E0*2*I_ESP_G2xyz_F3z_a;
    abcd[iGrid*900+590] = 4.0E0*I_ESP_I2x2y2z_F3z_aa-2.0E0*1*I_ESP_G2x2z_F3z_a;
    abcd[iGrid*900+591] = 4.0E0*I_ESP_Ix5y_F3z_aa-2.0E0*3*I_ESP_Gx3y_F3z_a-2.0E0*4*I_ESP_Gx3y_F3z_a+3*2*I_ESP_Dxy_F3z;
    abcd[iGrid*900+592] = 4.0E0*I_ESP_Ix4yz_F3z_aa-2.0E0*2*I_ESP_Gx2yz_F3z_a-2.0E0*3*I_ESP_Gx2yz_F3z_a+2*1*I_ESP_Dxz_F3z;
    abcd[iGrid*900+593] = 4.0E0*I_ESP_Ix3y2z_F3z_aa-2.0E0*1*I_ESP_Gxy2z_F3z_a-2.0E0*2*I_ESP_Gxy2z_F3z_a;
    abcd[iGrid*900+594] = 4.0E0*I_ESP_Ix2y3z_F3z_aa-2.0E0*1*I_ESP_Gx3z_F3z_a;
    abcd[iGrid*900+595] = 4.0E0*I_ESP_I6y_F3z_aa-2.0E0*4*I_ESP_G4y_F3z_a-2.0E0*5*I_ESP_G4y_F3z_a+4*3*I_ESP_D2y_F3z;
    abcd[iGrid*900+596] = 4.0E0*I_ESP_I5yz_F3z_aa-2.0E0*3*I_ESP_G3yz_F3z_a-2.0E0*4*I_ESP_G3yz_F3z_a+3*2*I_ESP_Dyz_F3z;
    abcd[iGrid*900+597] = 4.0E0*I_ESP_I4y2z_F3z_aa-2.0E0*2*I_ESP_G2y2z_F3z_a-2.0E0*3*I_ESP_G2y2z_F3z_a+2*1*I_ESP_D2z_F3z;
    abcd[iGrid*900+598] = 4.0E0*I_ESP_I3y3z_F3z_aa-2.0E0*1*I_ESP_Gy3z_F3z_a-2.0E0*2*I_ESP_Gy3z_F3z_a;
    abcd[iGrid*900+599] = 4.0E0*I_ESP_I2y4z_F3z_aa-2.0E0*1*I_ESP_G4z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_aa
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_D_F
     ************************************************************/
    abcd[iGrid*900+600] = 4.0E0*I_ESP_I4xyz_F3x_aa;
    abcd[iGrid*900+601] = 4.0E0*I_ESP_I3x2yz_F3x_aa-2.0E0*1*I_ESP_G3xz_F3x_a;
    abcd[iGrid*900+602] = 4.0E0*I_ESP_I3xy2z_F3x_aa-2.0E0*1*I_ESP_G3xy_F3x_a;
    abcd[iGrid*900+603] = 4.0E0*I_ESP_I2x3yz_F3x_aa-2.0E0*2*I_ESP_G2xyz_F3x_a;
    abcd[iGrid*900+604] = 4.0E0*I_ESP_I2x2y2z_F3x_aa-2.0E0*1*I_ESP_G2x2y_F3x_a-2.0E0*1*I_ESP_G2x2z_F3x_a+1*I_ESP_D2x_F3x;
    abcd[iGrid*900+605] = 4.0E0*I_ESP_I2xy3z_F3x_aa-2.0E0*2*I_ESP_G2xyz_F3x_a;
    abcd[iGrid*900+606] = 4.0E0*I_ESP_Ix4yz_F3x_aa-2.0E0*3*I_ESP_Gx2yz_F3x_a;
    abcd[iGrid*900+607] = 4.0E0*I_ESP_Ix3y2z_F3x_aa-2.0E0*1*I_ESP_Gx3y_F3x_a-2.0E0*2*I_ESP_Gxy2z_F3x_a+2*1*I_ESP_Dxy_F3x;
    abcd[iGrid*900+608] = 4.0E0*I_ESP_Ix2y3z_F3x_aa-2.0E0*2*I_ESP_Gx2yz_F3x_a-2.0E0*1*I_ESP_Gx3z_F3x_a+2*I_ESP_Dxz_F3x;
    abcd[iGrid*900+609] = 4.0E0*I_ESP_Ixy4z_F3x_aa-2.0E0*3*I_ESP_Gxy2z_F3x_a;
    abcd[iGrid*900+610] = 4.0E0*I_ESP_I5yz_F3x_aa-2.0E0*4*I_ESP_G3yz_F3x_a;
    abcd[iGrid*900+611] = 4.0E0*I_ESP_I4y2z_F3x_aa-2.0E0*1*I_ESP_G4y_F3x_a-2.0E0*3*I_ESP_G2y2z_F3x_a+3*1*I_ESP_D2y_F3x;
    abcd[iGrid*900+612] = 4.0E0*I_ESP_I3y3z_F3x_aa-2.0E0*2*I_ESP_G3yz_F3x_a-2.0E0*2*I_ESP_Gy3z_F3x_a+2*2*I_ESP_Dyz_F3x;
    abcd[iGrid*900+613] = 4.0E0*I_ESP_I2y4z_F3x_aa-2.0E0*3*I_ESP_G2y2z_F3x_a-2.0E0*1*I_ESP_G4z_F3x_a+3*I_ESP_D2z_F3x;
    abcd[iGrid*900+614] = 4.0E0*I_ESP_Iy5z_F3x_aa-2.0E0*4*I_ESP_Gy3z_F3x_a;
    abcd[iGrid*900+615] = 4.0E0*I_ESP_I4xyz_F2xy_aa;
    abcd[iGrid*900+616] = 4.0E0*I_ESP_I3x2yz_F2xy_aa-2.0E0*1*I_ESP_G3xz_F2xy_a;
    abcd[iGrid*900+617] = 4.0E0*I_ESP_I3xy2z_F2xy_aa-2.0E0*1*I_ESP_G3xy_F2xy_a;
    abcd[iGrid*900+618] = 4.0E0*I_ESP_I2x3yz_F2xy_aa-2.0E0*2*I_ESP_G2xyz_F2xy_a;
    abcd[iGrid*900+619] = 4.0E0*I_ESP_I2x2y2z_F2xy_aa-2.0E0*1*I_ESP_G2x2y_F2xy_a-2.0E0*1*I_ESP_G2x2z_F2xy_a+1*I_ESP_D2x_F2xy;
    abcd[iGrid*900+620] = 4.0E0*I_ESP_I2xy3z_F2xy_aa-2.0E0*2*I_ESP_G2xyz_F2xy_a;
    abcd[iGrid*900+621] = 4.0E0*I_ESP_Ix4yz_F2xy_aa-2.0E0*3*I_ESP_Gx2yz_F2xy_a;
    abcd[iGrid*900+622] = 4.0E0*I_ESP_Ix3y2z_F2xy_aa-2.0E0*1*I_ESP_Gx3y_F2xy_a-2.0E0*2*I_ESP_Gxy2z_F2xy_a+2*1*I_ESP_Dxy_F2xy;
    abcd[iGrid*900+623] = 4.0E0*I_ESP_Ix2y3z_F2xy_aa-2.0E0*2*I_ESP_Gx2yz_F2xy_a-2.0E0*1*I_ESP_Gx3z_F2xy_a+2*I_ESP_Dxz_F2xy;
    abcd[iGrid*900+624] = 4.0E0*I_ESP_Ixy4z_F2xy_aa-2.0E0*3*I_ESP_Gxy2z_F2xy_a;
    abcd[iGrid*900+625] = 4.0E0*I_ESP_I5yz_F2xy_aa-2.0E0*4*I_ESP_G3yz_F2xy_a;
    abcd[iGrid*900+626] = 4.0E0*I_ESP_I4y2z_F2xy_aa-2.0E0*1*I_ESP_G4y_F2xy_a-2.0E0*3*I_ESP_G2y2z_F2xy_a+3*1*I_ESP_D2y_F2xy;
    abcd[iGrid*900+627] = 4.0E0*I_ESP_I3y3z_F2xy_aa-2.0E0*2*I_ESP_G3yz_F2xy_a-2.0E0*2*I_ESP_Gy3z_F2xy_a+2*2*I_ESP_Dyz_F2xy;
    abcd[iGrid*900+628] = 4.0E0*I_ESP_I2y4z_F2xy_aa-2.0E0*3*I_ESP_G2y2z_F2xy_a-2.0E0*1*I_ESP_G4z_F2xy_a+3*I_ESP_D2z_F2xy;
    abcd[iGrid*900+629] = 4.0E0*I_ESP_Iy5z_F2xy_aa-2.0E0*4*I_ESP_Gy3z_F2xy_a;
    abcd[iGrid*900+630] = 4.0E0*I_ESP_I4xyz_F2xz_aa;
    abcd[iGrid*900+631] = 4.0E0*I_ESP_I3x2yz_F2xz_aa-2.0E0*1*I_ESP_G3xz_F2xz_a;
    abcd[iGrid*900+632] = 4.0E0*I_ESP_I3xy2z_F2xz_aa-2.0E0*1*I_ESP_G3xy_F2xz_a;
    abcd[iGrid*900+633] = 4.0E0*I_ESP_I2x3yz_F2xz_aa-2.0E0*2*I_ESP_G2xyz_F2xz_a;
    abcd[iGrid*900+634] = 4.0E0*I_ESP_I2x2y2z_F2xz_aa-2.0E0*1*I_ESP_G2x2y_F2xz_a-2.0E0*1*I_ESP_G2x2z_F2xz_a+1*I_ESP_D2x_F2xz;
    abcd[iGrid*900+635] = 4.0E0*I_ESP_I2xy3z_F2xz_aa-2.0E0*2*I_ESP_G2xyz_F2xz_a;
    abcd[iGrid*900+636] = 4.0E0*I_ESP_Ix4yz_F2xz_aa-2.0E0*3*I_ESP_Gx2yz_F2xz_a;
    abcd[iGrid*900+637] = 4.0E0*I_ESP_Ix3y2z_F2xz_aa-2.0E0*1*I_ESP_Gx3y_F2xz_a-2.0E0*2*I_ESP_Gxy2z_F2xz_a+2*1*I_ESP_Dxy_F2xz;
    abcd[iGrid*900+638] = 4.0E0*I_ESP_Ix2y3z_F2xz_aa-2.0E0*2*I_ESP_Gx2yz_F2xz_a-2.0E0*1*I_ESP_Gx3z_F2xz_a+2*I_ESP_Dxz_F2xz;
    abcd[iGrid*900+639] = 4.0E0*I_ESP_Ixy4z_F2xz_aa-2.0E0*3*I_ESP_Gxy2z_F2xz_a;
    abcd[iGrid*900+640] = 4.0E0*I_ESP_I5yz_F2xz_aa-2.0E0*4*I_ESP_G3yz_F2xz_a;
    abcd[iGrid*900+641] = 4.0E0*I_ESP_I4y2z_F2xz_aa-2.0E0*1*I_ESP_G4y_F2xz_a-2.0E0*3*I_ESP_G2y2z_F2xz_a+3*1*I_ESP_D2y_F2xz;
    abcd[iGrid*900+642] = 4.0E0*I_ESP_I3y3z_F2xz_aa-2.0E0*2*I_ESP_G3yz_F2xz_a-2.0E0*2*I_ESP_Gy3z_F2xz_a+2*2*I_ESP_Dyz_F2xz;
    abcd[iGrid*900+643] = 4.0E0*I_ESP_I2y4z_F2xz_aa-2.0E0*3*I_ESP_G2y2z_F2xz_a-2.0E0*1*I_ESP_G4z_F2xz_a+3*I_ESP_D2z_F2xz;
    abcd[iGrid*900+644] = 4.0E0*I_ESP_Iy5z_F2xz_aa-2.0E0*4*I_ESP_Gy3z_F2xz_a;
    abcd[iGrid*900+645] = 4.0E0*I_ESP_I4xyz_Fx2y_aa;
    abcd[iGrid*900+646] = 4.0E0*I_ESP_I3x2yz_Fx2y_aa-2.0E0*1*I_ESP_G3xz_Fx2y_a;
    abcd[iGrid*900+647] = 4.0E0*I_ESP_I3xy2z_Fx2y_aa-2.0E0*1*I_ESP_G3xy_Fx2y_a;
    abcd[iGrid*900+648] = 4.0E0*I_ESP_I2x3yz_Fx2y_aa-2.0E0*2*I_ESP_G2xyz_Fx2y_a;
    abcd[iGrid*900+649] = 4.0E0*I_ESP_I2x2y2z_Fx2y_aa-2.0E0*1*I_ESP_G2x2y_Fx2y_a-2.0E0*1*I_ESP_G2x2z_Fx2y_a+1*I_ESP_D2x_Fx2y;
    abcd[iGrid*900+650] = 4.0E0*I_ESP_I2xy3z_Fx2y_aa-2.0E0*2*I_ESP_G2xyz_Fx2y_a;
    abcd[iGrid*900+651] = 4.0E0*I_ESP_Ix4yz_Fx2y_aa-2.0E0*3*I_ESP_Gx2yz_Fx2y_a;
    abcd[iGrid*900+652] = 4.0E0*I_ESP_Ix3y2z_Fx2y_aa-2.0E0*1*I_ESP_Gx3y_Fx2y_a-2.0E0*2*I_ESP_Gxy2z_Fx2y_a+2*1*I_ESP_Dxy_Fx2y;
    abcd[iGrid*900+653] = 4.0E0*I_ESP_Ix2y3z_Fx2y_aa-2.0E0*2*I_ESP_Gx2yz_Fx2y_a-2.0E0*1*I_ESP_Gx3z_Fx2y_a+2*I_ESP_Dxz_Fx2y;
    abcd[iGrid*900+654] = 4.0E0*I_ESP_Ixy4z_Fx2y_aa-2.0E0*3*I_ESP_Gxy2z_Fx2y_a;
    abcd[iGrid*900+655] = 4.0E0*I_ESP_I5yz_Fx2y_aa-2.0E0*4*I_ESP_G3yz_Fx2y_a;
    abcd[iGrid*900+656] = 4.0E0*I_ESP_I4y2z_Fx2y_aa-2.0E0*1*I_ESP_G4y_Fx2y_a-2.0E0*3*I_ESP_G2y2z_Fx2y_a+3*1*I_ESP_D2y_Fx2y;
    abcd[iGrid*900+657] = 4.0E0*I_ESP_I3y3z_Fx2y_aa-2.0E0*2*I_ESP_G3yz_Fx2y_a-2.0E0*2*I_ESP_Gy3z_Fx2y_a+2*2*I_ESP_Dyz_Fx2y;
    abcd[iGrid*900+658] = 4.0E0*I_ESP_I2y4z_Fx2y_aa-2.0E0*3*I_ESP_G2y2z_Fx2y_a-2.0E0*1*I_ESP_G4z_Fx2y_a+3*I_ESP_D2z_Fx2y;
    abcd[iGrid*900+659] = 4.0E0*I_ESP_Iy5z_Fx2y_aa-2.0E0*4*I_ESP_Gy3z_Fx2y_a;
    abcd[iGrid*900+660] = 4.0E0*I_ESP_I4xyz_Fxyz_aa;
    abcd[iGrid*900+661] = 4.0E0*I_ESP_I3x2yz_Fxyz_aa-2.0E0*1*I_ESP_G3xz_Fxyz_a;
    abcd[iGrid*900+662] = 4.0E0*I_ESP_I3xy2z_Fxyz_aa-2.0E0*1*I_ESP_G3xy_Fxyz_a;
    abcd[iGrid*900+663] = 4.0E0*I_ESP_I2x3yz_Fxyz_aa-2.0E0*2*I_ESP_G2xyz_Fxyz_a;
    abcd[iGrid*900+664] = 4.0E0*I_ESP_I2x2y2z_Fxyz_aa-2.0E0*1*I_ESP_G2x2y_Fxyz_a-2.0E0*1*I_ESP_G2x2z_Fxyz_a+1*I_ESP_D2x_Fxyz;
    abcd[iGrid*900+665] = 4.0E0*I_ESP_I2xy3z_Fxyz_aa-2.0E0*2*I_ESP_G2xyz_Fxyz_a;
    abcd[iGrid*900+666] = 4.0E0*I_ESP_Ix4yz_Fxyz_aa-2.0E0*3*I_ESP_Gx2yz_Fxyz_a;
    abcd[iGrid*900+667] = 4.0E0*I_ESP_Ix3y2z_Fxyz_aa-2.0E0*1*I_ESP_Gx3y_Fxyz_a-2.0E0*2*I_ESP_Gxy2z_Fxyz_a+2*1*I_ESP_Dxy_Fxyz;
    abcd[iGrid*900+668] = 4.0E0*I_ESP_Ix2y3z_Fxyz_aa-2.0E0*2*I_ESP_Gx2yz_Fxyz_a-2.0E0*1*I_ESP_Gx3z_Fxyz_a+2*I_ESP_Dxz_Fxyz;
    abcd[iGrid*900+669] = 4.0E0*I_ESP_Ixy4z_Fxyz_aa-2.0E0*3*I_ESP_Gxy2z_Fxyz_a;
    abcd[iGrid*900+670] = 4.0E0*I_ESP_I5yz_Fxyz_aa-2.0E0*4*I_ESP_G3yz_Fxyz_a;
    abcd[iGrid*900+671] = 4.0E0*I_ESP_I4y2z_Fxyz_aa-2.0E0*1*I_ESP_G4y_Fxyz_a-2.0E0*3*I_ESP_G2y2z_Fxyz_a+3*1*I_ESP_D2y_Fxyz;
    abcd[iGrid*900+672] = 4.0E0*I_ESP_I3y3z_Fxyz_aa-2.0E0*2*I_ESP_G3yz_Fxyz_a-2.0E0*2*I_ESP_Gy3z_Fxyz_a+2*2*I_ESP_Dyz_Fxyz;
    abcd[iGrid*900+673] = 4.0E0*I_ESP_I2y4z_Fxyz_aa-2.0E0*3*I_ESP_G2y2z_Fxyz_a-2.0E0*1*I_ESP_G4z_Fxyz_a+3*I_ESP_D2z_Fxyz;
    abcd[iGrid*900+674] = 4.0E0*I_ESP_Iy5z_Fxyz_aa-2.0E0*4*I_ESP_Gy3z_Fxyz_a;
    abcd[iGrid*900+675] = 4.0E0*I_ESP_I4xyz_Fx2z_aa;
    abcd[iGrid*900+676] = 4.0E0*I_ESP_I3x2yz_Fx2z_aa-2.0E0*1*I_ESP_G3xz_Fx2z_a;
    abcd[iGrid*900+677] = 4.0E0*I_ESP_I3xy2z_Fx2z_aa-2.0E0*1*I_ESP_G3xy_Fx2z_a;
    abcd[iGrid*900+678] = 4.0E0*I_ESP_I2x3yz_Fx2z_aa-2.0E0*2*I_ESP_G2xyz_Fx2z_a;
    abcd[iGrid*900+679] = 4.0E0*I_ESP_I2x2y2z_Fx2z_aa-2.0E0*1*I_ESP_G2x2y_Fx2z_a-2.0E0*1*I_ESP_G2x2z_Fx2z_a+1*I_ESP_D2x_Fx2z;
    abcd[iGrid*900+680] = 4.0E0*I_ESP_I2xy3z_Fx2z_aa-2.0E0*2*I_ESP_G2xyz_Fx2z_a;
    abcd[iGrid*900+681] = 4.0E0*I_ESP_Ix4yz_Fx2z_aa-2.0E0*3*I_ESP_Gx2yz_Fx2z_a;
    abcd[iGrid*900+682] = 4.0E0*I_ESP_Ix3y2z_Fx2z_aa-2.0E0*1*I_ESP_Gx3y_Fx2z_a-2.0E0*2*I_ESP_Gxy2z_Fx2z_a+2*1*I_ESP_Dxy_Fx2z;
    abcd[iGrid*900+683] = 4.0E0*I_ESP_Ix2y3z_Fx2z_aa-2.0E0*2*I_ESP_Gx2yz_Fx2z_a-2.0E0*1*I_ESP_Gx3z_Fx2z_a+2*I_ESP_Dxz_Fx2z;
    abcd[iGrid*900+684] = 4.0E0*I_ESP_Ixy4z_Fx2z_aa-2.0E0*3*I_ESP_Gxy2z_Fx2z_a;
    abcd[iGrid*900+685] = 4.0E0*I_ESP_I5yz_Fx2z_aa-2.0E0*4*I_ESP_G3yz_Fx2z_a;
    abcd[iGrid*900+686] = 4.0E0*I_ESP_I4y2z_Fx2z_aa-2.0E0*1*I_ESP_G4y_Fx2z_a-2.0E0*3*I_ESP_G2y2z_Fx2z_a+3*1*I_ESP_D2y_Fx2z;
    abcd[iGrid*900+687] = 4.0E0*I_ESP_I3y3z_Fx2z_aa-2.0E0*2*I_ESP_G3yz_Fx2z_a-2.0E0*2*I_ESP_Gy3z_Fx2z_a+2*2*I_ESP_Dyz_Fx2z;
    abcd[iGrid*900+688] = 4.0E0*I_ESP_I2y4z_Fx2z_aa-2.0E0*3*I_ESP_G2y2z_Fx2z_a-2.0E0*1*I_ESP_G4z_Fx2z_a+3*I_ESP_D2z_Fx2z;
    abcd[iGrid*900+689] = 4.0E0*I_ESP_Iy5z_Fx2z_aa-2.0E0*4*I_ESP_Gy3z_Fx2z_a;
    abcd[iGrid*900+690] = 4.0E0*I_ESP_I4xyz_F3y_aa;
    abcd[iGrid*900+691] = 4.0E0*I_ESP_I3x2yz_F3y_aa-2.0E0*1*I_ESP_G3xz_F3y_a;
    abcd[iGrid*900+692] = 4.0E0*I_ESP_I3xy2z_F3y_aa-2.0E0*1*I_ESP_G3xy_F3y_a;
    abcd[iGrid*900+693] = 4.0E0*I_ESP_I2x3yz_F3y_aa-2.0E0*2*I_ESP_G2xyz_F3y_a;
    abcd[iGrid*900+694] = 4.0E0*I_ESP_I2x2y2z_F3y_aa-2.0E0*1*I_ESP_G2x2y_F3y_a-2.0E0*1*I_ESP_G2x2z_F3y_a+1*I_ESP_D2x_F3y;
    abcd[iGrid*900+695] = 4.0E0*I_ESP_I2xy3z_F3y_aa-2.0E0*2*I_ESP_G2xyz_F3y_a;
    abcd[iGrid*900+696] = 4.0E0*I_ESP_Ix4yz_F3y_aa-2.0E0*3*I_ESP_Gx2yz_F3y_a;
    abcd[iGrid*900+697] = 4.0E0*I_ESP_Ix3y2z_F3y_aa-2.0E0*1*I_ESP_Gx3y_F3y_a-2.0E0*2*I_ESP_Gxy2z_F3y_a+2*1*I_ESP_Dxy_F3y;
    abcd[iGrid*900+698] = 4.0E0*I_ESP_Ix2y3z_F3y_aa-2.0E0*2*I_ESP_Gx2yz_F3y_a-2.0E0*1*I_ESP_Gx3z_F3y_a+2*I_ESP_Dxz_F3y;
    abcd[iGrid*900+699] = 4.0E0*I_ESP_Ixy4z_F3y_aa-2.0E0*3*I_ESP_Gxy2z_F3y_a;
    abcd[iGrid*900+700] = 4.0E0*I_ESP_I5yz_F3y_aa-2.0E0*4*I_ESP_G3yz_F3y_a;
    abcd[iGrid*900+701] = 4.0E0*I_ESP_I4y2z_F3y_aa-2.0E0*1*I_ESP_G4y_F3y_a-2.0E0*3*I_ESP_G2y2z_F3y_a+3*1*I_ESP_D2y_F3y;
    abcd[iGrid*900+702] = 4.0E0*I_ESP_I3y3z_F3y_aa-2.0E0*2*I_ESP_G3yz_F3y_a-2.0E0*2*I_ESP_Gy3z_F3y_a+2*2*I_ESP_Dyz_F3y;
    abcd[iGrid*900+703] = 4.0E0*I_ESP_I2y4z_F3y_aa-2.0E0*3*I_ESP_G2y2z_F3y_a-2.0E0*1*I_ESP_G4z_F3y_a+3*I_ESP_D2z_F3y;
    abcd[iGrid*900+704] = 4.0E0*I_ESP_Iy5z_F3y_aa-2.0E0*4*I_ESP_Gy3z_F3y_a;
    abcd[iGrid*900+705] = 4.0E0*I_ESP_I4xyz_F2yz_aa;
    abcd[iGrid*900+706] = 4.0E0*I_ESP_I3x2yz_F2yz_aa-2.0E0*1*I_ESP_G3xz_F2yz_a;
    abcd[iGrid*900+707] = 4.0E0*I_ESP_I3xy2z_F2yz_aa-2.0E0*1*I_ESP_G3xy_F2yz_a;
    abcd[iGrid*900+708] = 4.0E0*I_ESP_I2x3yz_F2yz_aa-2.0E0*2*I_ESP_G2xyz_F2yz_a;
    abcd[iGrid*900+709] = 4.0E0*I_ESP_I2x2y2z_F2yz_aa-2.0E0*1*I_ESP_G2x2y_F2yz_a-2.0E0*1*I_ESP_G2x2z_F2yz_a+1*I_ESP_D2x_F2yz;
    abcd[iGrid*900+710] = 4.0E0*I_ESP_I2xy3z_F2yz_aa-2.0E0*2*I_ESP_G2xyz_F2yz_a;
    abcd[iGrid*900+711] = 4.0E0*I_ESP_Ix4yz_F2yz_aa-2.0E0*3*I_ESP_Gx2yz_F2yz_a;
    abcd[iGrid*900+712] = 4.0E0*I_ESP_Ix3y2z_F2yz_aa-2.0E0*1*I_ESP_Gx3y_F2yz_a-2.0E0*2*I_ESP_Gxy2z_F2yz_a+2*1*I_ESP_Dxy_F2yz;
    abcd[iGrid*900+713] = 4.0E0*I_ESP_Ix2y3z_F2yz_aa-2.0E0*2*I_ESP_Gx2yz_F2yz_a-2.0E0*1*I_ESP_Gx3z_F2yz_a+2*I_ESP_Dxz_F2yz;
    abcd[iGrid*900+714] = 4.0E0*I_ESP_Ixy4z_F2yz_aa-2.0E0*3*I_ESP_Gxy2z_F2yz_a;
    abcd[iGrid*900+715] = 4.0E0*I_ESP_I5yz_F2yz_aa-2.0E0*4*I_ESP_G3yz_F2yz_a;
    abcd[iGrid*900+716] = 4.0E0*I_ESP_I4y2z_F2yz_aa-2.0E0*1*I_ESP_G4y_F2yz_a-2.0E0*3*I_ESP_G2y2z_F2yz_a+3*1*I_ESP_D2y_F2yz;
    abcd[iGrid*900+717] = 4.0E0*I_ESP_I3y3z_F2yz_aa-2.0E0*2*I_ESP_G3yz_F2yz_a-2.0E0*2*I_ESP_Gy3z_F2yz_a+2*2*I_ESP_Dyz_F2yz;
    abcd[iGrid*900+718] = 4.0E0*I_ESP_I2y4z_F2yz_aa-2.0E0*3*I_ESP_G2y2z_F2yz_a-2.0E0*1*I_ESP_G4z_F2yz_a+3*I_ESP_D2z_F2yz;
    abcd[iGrid*900+719] = 4.0E0*I_ESP_Iy5z_F2yz_aa-2.0E0*4*I_ESP_Gy3z_F2yz_a;
    abcd[iGrid*900+720] = 4.0E0*I_ESP_I4xyz_Fy2z_aa;
    abcd[iGrid*900+721] = 4.0E0*I_ESP_I3x2yz_Fy2z_aa-2.0E0*1*I_ESP_G3xz_Fy2z_a;
    abcd[iGrid*900+722] = 4.0E0*I_ESP_I3xy2z_Fy2z_aa-2.0E0*1*I_ESP_G3xy_Fy2z_a;
    abcd[iGrid*900+723] = 4.0E0*I_ESP_I2x3yz_Fy2z_aa-2.0E0*2*I_ESP_G2xyz_Fy2z_a;
    abcd[iGrid*900+724] = 4.0E0*I_ESP_I2x2y2z_Fy2z_aa-2.0E0*1*I_ESP_G2x2y_Fy2z_a-2.0E0*1*I_ESP_G2x2z_Fy2z_a+1*I_ESP_D2x_Fy2z;
    abcd[iGrid*900+725] = 4.0E0*I_ESP_I2xy3z_Fy2z_aa-2.0E0*2*I_ESP_G2xyz_Fy2z_a;
    abcd[iGrid*900+726] = 4.0E0*I_ESP_Ix4yz_Fy2z_aa-2.0E0*3*I_ESP_Gx2yz_Fy2z_a;
    abcd[iGrid*900+727] = 4.0E0*I_ESP_Ix3y2z_Fy2z_aa-2.0E0*1*I_ESP_Gx3y_Fy2z_a-2.0E0*2*I_ESP_Gxy2z_Fy2z_a+2*1*I_ESP_Dxy_Fy2z;
    abcd[iGrid*900+728] = 4.0E0*I_ESP_Ix2y3z_Fy2z_aa-2.0E0*2*I_ESP_Gx2yz_Fy2z_a-2.0E0*1*I_ESP_Gx3z_Fy2z_a+2*I_ESP_Dxz_Fy2z;
    abcd[iGrid*900+729] = 4.0E0*I_ESP_Ixy4z_Fy2z_aa-2.0E0*3*I_ESP_Gxy2z_Fy2z_a;
    abcd[iGrid*900+730] = 4.0E0*I_ESP_I5yz_Fy2z_aa-2.0E0*4*I_ESP_G3yz_Fy2z_a;
    abcd[iGrid*900+731] = 4.0E0*I_ESP_I4y2z_Fy2z_aa-2.0E0*1*I_ESP_G4y_Fy2z_a-2.0E0*3*I_ESP_G2y2z_Fy2z_a+3*1*I_ESP_D2y_Fy2z;
    abcd[iGrid*900+732] = 4.0E0*I_ESP_I3y3z_Fy2z_aa-2.0E0*2*I_ESP_G3yz_Fy2z_a-2.0E0*2*I_ESP_Gy3z_Fy2z_a+2*2*I_ESP_Dyz_Fy2z;
    abcd[iGrid*900+733] = 4.0E0*I_ESP_I2y4z_Fy2z_aa-2.0E0*3*I_ESP_G2y2z_Fy2z_a-2.0E0*1*I_ESP_G4z_Fy2z_a+3*I_ESP_D2z_Fy2z;
    abcd[iGrid*900+734] = 4.0E0*I_ESP_Iy5z_Fy2z_aa-2.0E0*4*I_ESP_Gy3z_Fy2z_a;
    abcd[iGrid*900+735] = 4.0E0*I_ESP_I4xyz_F3z_aa;
    abcd[iGrid*900+736] = 4.0E0*I_ESP_I3x2yz_F3z_aa-2.0E0*1*I_ESP_G3xz_F3z_a;
    abcd[iGrid*900+737] = 4.0E0*I_ESP_I3xy2z_F3z_aa-2.0E0*1*I_ESP_G3xy_F3z_a;
    abcd[iGrid*900+738] = 4.0E0*I_ESP_I2x3yz_F3z_aa-2.0E0*2*I_ESP_G2xyz_F3z_a;
    abcd[iGrid*900+739] = 4.0E0*I_ESP_I2x2y2z_F3z_aa-2.0E0*1*I_ESP_G2x2y_F3z_a-2.0E0*1*I_ESP_G2x2z_F3z_a+1*I_ESP_D2x_F3z;
    abcd[iGrid*900+740] = 4.0E0*I_ESP_I2xy3z_F3z_aa-2.0E0*2*I_ESP_G2xyz_F3z_a;
    abcd[iGrid*900+741] = 4.0E0*I_ESP_Ix4yz_F3z_aa-2.0E0*3*I_ESP_Gx2yz_F3z_a;
    abcd[iGrid*900+742] = 4.0E0*I_ESP_Ix3y2z_F3z_aa-2.0E0*1*I_ESP_Gx3y_F3z_a-2.0E0*2*I_ESP_Gxy2z_F3z_a+2*1*I_ESP_Dxy_F3z;
    abcd[iGrid*900+743] = 4.0E0*I_ESP_Ix2y3z_F3z_aa-2.0E0*2*I_ESP_Gx2yz_F3z_a-2.0E0*1*I_ESP_Gx3z_F3z_a+2*I_ESP_Dxz_F3z;
    abcd[iGrid*900+744] = 4.0E0*I_ESP_Ixy4z_F3z_aa-2.0E0*3*I_ESP_Gxy2z_F3z_a;
    abcd[iGrid*900+745] = 4.0E0*I_ESP_I5yz_F3z_aa-2.0E0*4*I_ESP_G3yz_F3z_a;
    abcd[iGrid*900+746] = 4.0E0*I_ESP_I4y2z_F3z_aa-2.0E0*1*I_ESP_G4y_F3z_a-2.0E0*3*I_ESP_G2y2z_F3z_a+3*1*I_ESP_D2y_F3z;
    abcd[iGrid*900+747] = 4.0E0*I_ESP_I3y3z_F3z_aa-2.0E0*2*I_ESP_G3yz_F3z_a-2.0E0*2*I_ESP_Gy3z_F3z_a+2*2*I_ESP_Dyz_F3z;
    abcd[iGrid*900+748] = 4.0E0*I_ESP_I2y4z_F3z_aa-2.0E0*3*I_ESP_G2y2z_F3z_a-2.0E0*1*I_ESP_G4z_F3z_a+3*I_ESP_D2z_F3z;
    abcd[iGrid*900+749] = 4.0E0*I_ESP_Iy5z_F3z_aa-2.0E0*4*I_ESP_Gy3z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_aa
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_G_F_a
     * RHS shell quartet name: SQ_ESP_D_F
     ************************************************************/
    abcd[iGrid*900+750] = 4.0E0*I_ESP_I4x2z_F3x_aa-2.0E0*1*I_ESP_G4x_F3x_a;
    abcd[iGrid*900+751] = 4.0E0*I_ESP_I3xy2z_F3x_aa-2.0E0*1*I_ESP_G3xy_F3x_a;
    abcd[iGrid*900+752] = 4.0E0*I_ESP_I3x3z_F3x_aa-2.0E0*1*I_ESP_G3xz_F3x_a-2.0E0*2*I_ESP_G3xz_F3x_a;
    abcd[iGrid*900+753] = 4.0E0*I_ESP_I2x2y2z_F3x_aa-2.0E0*1*I_ESP_G2x2y_F3x_a;
    abcd[iGrid*900+754] = 4.0E0*I_ESP_I2xy3z_F3x_aa-2.0E0*1*I_ESP_G2xyz_F3x_a-2.0E0*2*I_ESP_G2xyz_F3x_a;
    abcd[iGrid*900+755] = 4.0E0*I_ESP_I2x4z_F3x_aa-2.0E0*2*I_ESP_G2x2z_F3x_a-2.0E0*3*I_ESP_G2x2z_F3x_a+2*1*I_ESP_D2x_F3x;
    abcd[iGrid*900+756] = 4.0E0*I_ESP_Ix3y2z_F3x_aa-2.0E0*1*I_ESP_Gx3y_F3x_a;
    abcd[iGrid*900+757] = 4.0E0*I_ESP_Ix2y3z_F3x_aa-2.0E0*1*I_ESP_Gx2yz_F3x_a-2.0E0*2*I_ESP_Gx2yz_F3x_a;
    abcd[iGrid*900+758] = 4.0E0*I_ESP_Ixy4z_F3x_aa-2.0E0*2*I_ESP_Gxy2z_F3x_a-2.0E0*3*I_ESP_Gxy2z_F3x_a+2*1*I_ESP_Dxy_F3x;
    abcd[iGrid*900+759] = 4.0E0*I_ESP_Ix5z_F3x_aa-2.0E0*3*I_ESP_Gx3z_F3x_a-2.0E0*4*I_ESP_Gx3z_F3x_a+3*2*I_ESP_Dxz_F3x;
    abcd[iGrid*900+760] = 4.0E0*I_ESP_I4y2z_F3x_aa-2.0E0*1*I_ESP_G4y_F3x_a;
    abcd[iGrid*900+761] = 4.0E0*I_ESP_I3y3z_F3x_aa-2.0E0*1*I_ESP_G3yz_F3x_a-2.0E0*2*I_ESP_G3yz_F3x_a;
    abcd[iGrid*900+762] = 4.0E0*I_ESP_I2y4z_F3x_aa-2.0E0*2*I_ESP_G2y2z_F3x_a-2.0E0*3*I_ESP_G2y2z_F3x_a+2*1*I_ESP_D2y_F3x;
    abcd[iGrid*900+763] = 4.0E0*I_ESP_Iy5z_F3x_aa-2.0E0*3*I_ESP_Gy3z_F3x_a-2.0E0*4*I_ESP_Gy3z_F3x_a+3*2*I_ESP_Dyz_F3x;
    abcd[iGrid*900+764] = 4.0E0*I_ESP_I6z_F3x_aa-2.0E0*4*I_ESP_G4z_F3x_a-2.0E0*5*I_ESP_G4z_F3x_a+4*3*I_ESP_D2z_F3x;
    abcd[iGrid*900+765] = 4.0E0*I_ESP_I4x2z_F2xy_aa-2.0E0*1*I_ESP_G4x_F2xy_a;
    abcd[iGrid*900+766] = 4.0E0*I_ESP_I3xy2z_F2xy_aa-2.0E0*1*I_ESP_G3xy_F2xy_a;
    abcd[iGrid*900+767] = 4.0E0*I_ESP_I3x3z_F2xy_aa-2.0E0*1*I_ESP_G3xz_F2xy_a-2.0E0*2*I_ESP_G3xz_F2xy_a;
    abcd[iGrid*900+768] = 4.0E0*I_ESP_I2x2y2z_F2xy_aa-2.0E0*1*I_ESP_G2x2y_F2xy_a;
    abcd[iGrid*900+769] = 4.0E0*I_ESP_I2xy3z_F2xy_aa-2.0E0*1*I_ESP_G2xyz_F2xy_a-2.0E0*2*I_ESP_G2xyz_F2xy_a;
    abcd[iGrid*900+770] = 4.0E0*I_ESP_I2x4z_F2xy_aa-2.0E0*2*I_ESP_G2x2z_F2xy_a-2.0E0*3*I_ESP_G2x2z_F2xy_a+2*1*I_ESP_D2x_F2xy;
    abcd[iGrid*900+771] = 4.0E0*I_ESP_Ix3y2z_F2xy_aa-2.0E0*1*I_ESP_Gx3y_F2xy_a;
    abcd[iGrid*900+772] = 4.0E0*I_ESP_Ix2y3z_F2xy_aa-2.0E0*1*I_ESP_Gx2yz_F2xy_a-2.0E0*2*I_ESP_Gx2yz_F2xy_a;
    abcd[iGrid*900+773] = 4.0E0*I_ESP_Ixy4z_F2xy_aa-2.0E0*2*I_ESP_Gxy2z_F2xy_a-2.0E0*3*I_ESP_Gxy2z_F2xy_a+2*1*I_ESP_Dxy_F2xy;
    abcd[iGrid*900+774] = 4.0E0*I_ESP_Ix5z_F2xy_aa-2.0E0*3*I_ESP_Gx3z_F2xy_a-2.0E0*4*I_ESP_Gx3z_F2xy_a+3*2*I_ESP_Dxz_F2xy;
    abcd[iGrid*900+775] = 4.0E0*I_ESP_I4y2z_F2xy_aa-2.0E0*1*I_ESP_G4y_F2xy_a;
    abcd[iGrid*900+776] = 4.0E0*I_ESP_I3y3z_F2xy_aa-2.0E0*1*I_ESP_G3yz_F2xy_a-2.0E0*2*I_ESP_G3yz_F2xy_a;
    abcd[iGrid*900+777] = 4.0E0*I_ESP_I2y4z_F2xy_aa-2.0E0*2*I_ESP_G2y2z_F2xy_a-2.0E0*3*I_ESP_G2y2z_F2xy_a+2*1*I_ESP_D2y_F2xy;
    abcd[iGrid*900+778] = 4.0E0*I_ESP_Iy5z_F2xy_aa-2.0E0*3*I_ESP_Gy3z_F2xy_a-2.0E0*4*I_ESP_Gy3z_F2xy_a+3*2*I_ESP_Dyz_F2xy;
    abcd[iGrid*900+779] = 4.0E0*I_ESP_I6z_F2xy_aa-2.0E0*4*I_ESP_G4z_F2xy_a-2.0E0*5*I_ESP_G4z_F2xy_a+4*3*I_ESP_D2z_F2xy;
    abcd[iGrid*900+780] = 4.0E0*I_ESP_I4x2z_F2xz_aa-2.0E0*1*I_ESP_G4x_F2xz_a;
    abcd[iGrid*900+781] = 4.0E0*I_ESP_I3xy2z_F2xz_aa-2.0E0*1*I_ESP_G3xy_F2xz_a;
    abcd[iGrid*900+782] = 4.0E0*I_ESP_I3x3z_F2xz_aa-2.0E0*1*I_ESP_G3xz_F2xz_a-2.0E0*2*I_ESP_G3xz_F2xz_a;
    abcd[iGrid*900+783] = 4.0E0*I_ESP_I2x2y2z_F2xz_aa-2.0E0*1*I_ESP_G2x2y_F2xz_a;
    abcd[iGrid*900+784] = 4.0E0*I_ESP_I2xy3z_F2xz_aa-2.0E0*1*I_ESP_G2xyz_F2xz_a-2.0E0*2*I_ESP_G2xyz_F2xz_a;
    abcd[iGrid*900+785] = 4.0E0*I_ESP_I2x4z_F2xz_aa-2.0E0*2*I_ESP_G2x2z_F2xz_a-2.0E0*3*I_ESP_G2x2z_F2xz_a+2*1*I_ESP_D2x_F2xz;
    abcd[iGrid*900+786] = 4.0E0*I_ESP_Ix3y2z_F2xz_aa-2.0E0*1*I_ESP_Gx3y_F2xz_a;
    abcd[iGrid*900+787] = 4.0E0*I_ESP_Ix2y3z_F2xz_aa-2.0E0*1*I_ESP_Gx2yz_F2xz_a-2.0E0*2*I_ESP_Gx2yz_F2xz_a;
    abcd[iGrid*900+788] = 4.0E0*I_ESP_Ixy4z_F2xz_aa-2.0E0*2*I_ESP_Gxy2z_F2xz_a-2.0E0*3*I_ESP_Gxy2z_F2xz_a+2*1*I_ESP_Dxy_F2xz;
    abcd[iGrid*900+789] = 4.0E0*I_ESP_Ix5z_F2xz_aa-2.0E0*3*I_ESP_Gx3z_F2xz_a-2.0E0*4*I_ESP_Gx3z_F2xz_a+3*2*I_ESP_Dxz_F2xz;
    abcd[iGrid*900+790] = 4.0E0*I_ESP_I4y2z_F2xz_aa-2.0E0*1*I_ESP_G4y_F2xz_a;
    abcd[iGrid*900+791] = 4.0E0*I_ESP_I3y3z_F2xz_aa-2.0E0*1*I_ESP_G3yz_F2xz_a-2.0E0*2*I_ESP_G3yz_F2xz_a;
    abcd[iGrid*900+792] = 4.0E0*I_ESP_I2y4z_F2xz_aa-2.0E0*2*I_ESP_G2y2z_F2xz_a-2.0E0*3*I_ESP_G2y2z_F2xz_a+2*1*I_ESP_D2y_F2xz;
    abcd[iGrid*900+793] = 4.0E0*I_ESP_Iy5z_F2xz_aa-2.0E0*3*I_ESP_Gy3z_F2xz_a-2.0E0*4*I_ESP_Gy3z_F2xz_a+3*2*I_ESP_Dyz_F2xz;
    abcd[iGrid*900+794] = 4.0E0*I_ESP_I6z_F2xz_aa-2.0E0*4*I_ESP_G4z_F2xz_a-2.0E0*5*I_ESP_G4z_F2xz_a+4*3*I_ESP_D2z_F2xz;
    abcd[iGrid*900+795] = 4.0E0*I_ESP_I4x2z_Fx2y_aa-2.0E0*1*I_ESP_G4x_Fx2y_a;
    abcd[iGrid*900+796] = 4.0E0*I_ESP_I3xy2z_Fx2y_aa-2.0E0*1*I_ESP_G3xy_Fx2y_a;
    abcd[iGrid*900+797] = 4.0E0*I_ESP_I3x3z_Fx2y_aa-2.0E0*1*I_ESP_G3xz_Fx2y_a-2.0E0*2*I_ESP_G3xz_Fx2y_a;
    abcd[iGrid*900+798] = 4.0E0*I_ESP_I2x2y2z_Fx2y_aa-2.0E0*1*I_ESP_G2x2y_Fx2y_a;
    abcd[iGrid*900+799] = 4.0E0*I_ESP_I2xy3z_Fx2y_aa-2.0E0*1*I_ESP_G2xyz_Fx2y_a-2.0E0*2*I_ESP_G2xyz_Fx2y_a;
    abcd[iGrid*900+800] = 4.0E0*I_ESP_I2x4z_Fx2y_aa-2.0E0*2*I_ESP_G2x2z_Fx2y_a-2.0E0*3*I_ESP_G2x2z_Fx2y_a+2*1*I_ESP_D2x_Fx2y;
    abcd[iGrid*900+801] = 4.0E0*I_ESP_Ix3y2z_Fx2y_aa-2.0E0*1*I_ESP_Gx3y_Fx2y_a;
    abcd[iGrid*900+802] = 4.0E0*I_ESP_Ix2y3z_Fx2y_aa-2.0E0*1*I_ESP_Gx2yz_Fx2y_a-2.0E0*2*I_ESP_Gx2yz_Fx2y_a;
    abcd[iGrid*900+803] = 4.0E0*I_ESP_Ixy4z_Fx2y_aa-2.0E0*2*I_ESP_Gxy2z_Fx2y_a-2.0E0*3*I_ESP_Gxy2z_Fx2y_a+2*1*I_ESP_Dxy_Fx2y;
    abcd[iGrid*900+804] = 4.0E0*I_ESP_Ix5z_Fx2y_aa-2.0E0*3*I_ESP_Gx3z_Fx2y_a-2.0E0*4*I_ESP_Gx3z_Fx2y_a+3*2*I_ESP_Dxz_Fx2y;
    abcd[iGrid*900+805] = 4.0E0*I_ESP_I4y2z_Fx2y_aa-2.0E0*1*I_ESP_G4y_Fx2y_a;
    abcd[iGrid*900+806] = 4.0E0*I_ESP_I3y3z_Fx2y_aa-2.0E0*1*I_ESP_G3yz_Fx2y_a-2.0E0*2*I_ESP_G3yz_Fx2y_a;
    abcd[iGrid*900+807] = 4.0E0*I_ESP_I2y4z_Fx2y_aa-2.0E0*2*I_ESP_G2y2z_Fx2y_a-2.0E0*3*I_ESP_G2y2z_Fx2y_a+2*1*I_ESP_D2y_Fx2y;
    abcd[iGrid*900+808] = 4.0E0*I_ESP_Iy5z_Fx2y_aa-2.0E0*3*I_ESP_Gy3z_Fx2y_a-2.0E0*4*I_ESP_Gy3z_Fx2y_a+3*2*I_ESP_Dyz_Fx2y;
    abcd[iGrid*900+809] = 4.0E0*I_ESP_I6z_Fx2y_aa-2.0E0*4*I_ESP_G4z_Fx2y_a-2.0E0*5*I_ESP_G4z_Fx2y_a+4*3*I_ESP_D2z_Fx2y;
    abcd[iGrid*900+810] = 4.0E0*I_ESP_I4x2z_Fxyz_aa-2.0E0*1*I_ESP_G4x_Fxyz_a;
    abcd[iGrid*900+811] = 4.0E0*I_ESP_I3xy2z_Fxyz_aa-2.0E0*1*I_ESP_G3xy_Fxyz_a;
    abcd[iGrid*900+812] = 4.0E0*I_ESP_I3x3z_Fxyz_aa-2.0E0*1*I_ESP_G3xz_Fxyz_a-2.0E0*2*I_ESP_G3xz_Fxyz_a;
    abcd[iGrid*900+813] = 4.0E0*I_ESP_I2x2y2z_Fxyz_aa-2.0E0*1*I_ESP_G2x2y_Fxyz_a;
    abcd[iGrid*900+814] = 4.0E0*I_ESP_I2xy3z_Fxyz_aa-2.0E0*1*I_ESP_G2xyz_Fxyz_a-2.0E0*2*I_ESP_G2xyz_Fxyz_a;
    abcd[iGrid*900+815] = 4.0E0*I_ESP_I2x4z_Fxyz_aa-2.0E0*2*I_ESP_G2x2z_Fxyz_a-2.0E0*3*I_ESP_G2x2z_Fxyz_a+2*1*I_ESP_D2x_Fxyz;
    abcd[iGrid*900+816] = 4.0E0*I_ESP_Ix3y2z_Fxyz_aa-2.0E0*1*I_ESP_Gx3y_Fxyz_a;
    abcd[iGrid*900+817] = 4.0E0*I_ESP_Ix2y3z_Fxyz_aa-2.0E0*1*I_ESP_Gx2yz_Fxyz_a-2.0E0*2*I_ESP_Gx2yz_Fxyz_a;
    abcd[iGrid*900+818] = 4.0E0*I_ESP_Ixy4z_Fxyz_aa-2.0E0*2*I_ESP_Gxy2z_Fxyz_a-2.0E0*3*I_ESP_Gxy2z_Fxyz_a+2*1*I_ESP_Dxy_Fxyz;
    abcd[iGrid*900+819] = 4.0E0*I_ESP_Ix5z_Fxyz_aa-2.0E0*3*I_ESP_Gx3z_Fxyz_a-2.0E0*4*I_ESP_Gx3z_Fxyz_a+3*2*I_ESP_Dxz_Fxyz;
    abcd[iGrid*900+820] = 4.0E0*I_ESP_I4y2z_Fxyz_aa-2.0E0*1*I_ESP_G4y_Fxyz_a;
    abcd[iGrid*900+821] = 4.0E0*I_ESP_I3y3z_Fxyz_aa-2.0E0*1*I_ESP_G3yz_Fxyz_a-2.0E0*2*I_ESP_G3yz_Fxyz_a;
    abcd[iGrid*900+822] = 4.0E0*I_ESP_I2y4z_Fxyz_aa-2.0E0*2*I_ESP_G2y2z_Fxyz_a-2.0E0*3*I_ESP_G2y2z_Fxyz_a+2*1*I_ESP_D2y_Fxyz;
    abcd[iGrid*900+823] = 4.0E0*I_ESP_Iy5z_Fxyz_aa-2.0E0*3*I_ESP_Gy3z_Fxyz_a-2.0E0*4*I_ESP_Gy3z_Fxyz_a+3*2*I_ESP_Dyz_Fxyz;
    abcd[iGrid*900+824] = 4.0E0*I_ESP_I6z_Fxyz_aa-2.0E0*4*I_ESP_G4z_Fxyz_a-2.0E0*5*I_ESP_G4z_Fxyz_a+4*3*I_ESP_D2z_Fxyz;
    abcd[iGrid*900+825] = 4.0E0*I_ESP_I4x2z_Fx2z_aa-2.0E0*1*I_ESP_G4x_Fx2z_a;
    abcd[iGrid*900+826] = 4.0E0*I_ESP_I3xy2z_Fx2z_aa-2.0E0*1*I_ESP_G3xy_Fx2z_a;
    abcd[iGrid*900+827] = 4.0E0*I_ESP_I3x3z_Fx2z_aa-2.0E0*1*I_ESP_G3xz_Fx2z_a-2.0E0*2*I_ESP_G3xz_Fx2z_a;
    abcd[iGrid*900+828] = 4.0E0*I_ESP_I2x2y2z_Fx2z_aa-2.0E0*1*I_ESP_G2x2y_Fx2z_a;
    abcd[iGrid*900+829] = 4.0E0*I_ESP_I2xy3z_Fx2z_aa-2.0E0*1*I_ESP_G2xyz_Fx2z_a-2.0E0*2*I_ESP_G2xyz_Fx2z_a;
    abcd[iGrid*900+830] = 4.0E0*I_ESP_I2x4z_Fx2z_aa-2.0E0*2*I_ESP_G2x2z_Fx2z_a-2.0E0*3*I_ESP_G2x2z_Fx2z_a+2*1*I_ESP_D2x_Fx2z;
    abcd[iGrid*900+831] = 4.0E0*I_ESP_Ix3y2z_Fx2z_aa-2.0E0*1*I_ESP_Gx3y_Fx2z_a;
    abcd[iGrid*900+832] = 4.0E0*I_ESP_Ix2y3z_Fx2z_aa-2.0E0*1*I_ESP_Gx2yz_Fx2z_a-2.0E0*2*I_ESP_Gx2yz_Fx2z_a;
    abcd[iGrid*900+833] = 4.0E0*I_ESP_Ixy4z_Fx2z_aa-2.0E0*2*I_ESP_Gxy2z_Fx2z_a-2.0E0*3*I_ESP_Gxy2z_Fx2z_a+2*1*I_ESP_Dxy_Fx2z;
    abcd[iGrid*900+834] = 4.0E0*I_ESP_Ix5z_Fx2z_aa-2.0E0*3*I_ESP_Gx3z_Fx2z_a-2.0E0*4*I_ESP_Gx3z_Fx2z_a+3*2*I_ESP_Dxz_Fx2z;
    abcd[iGrid*900+835] = 4.0E0*I_ESP_I4y2z_Fx2z_aa-2.0E0*1*I_ESP_G4y_Fx2z_a;
    abcd[iGrid*900+836] = 4.0E0*I_ESP_I3y3z_Fx2z_aa-2.0E0*1*I_ESP_G3yz_Fx2z_a-2.0E0*2*I_ESP_G3yz_Fx2z_a;
    abcd[iGrid*900+837] = 4.0E0*I_ESP_I2y4z_Fx2z_aa-2.0E0*2*I_ESP_G2y2z_Fx2z_a-2.0E0*3*I_ESP_G2y2z_Fx2z_a+2*1*I_ESP_D2y_Fx2z;
    abcd[iGrid*900+838] = 4.0E0*I_ESP_Iy5z_Fx2z_aa-2.0E0*3*I_ESP_Gy3z_Fx2z_a-2.0E0*4*I_ESP_Gy3z_Fx2z_a+3*2*I_ESP_Dyz_Fx2z;
    abcd[iGrid*900+839] = 4.0E0*I_ESP_I6z_Fx2z_aa-2.0E0*4*I_ESP_G4z_Fx2z_a-2.0E0*5*I_ESP_G4z_Fx2z_a+4*3*I_ESP_D2z_Fx2z;
    abcd[iGrid*900+840] = 4.0E0*I_ESP_I4x2z_F3y_aa-2.0E0*1*I_ESP_G4x_F3y_a;
    abcd[iGrid*900+841] = 4.0E0*I_ESP_I3xy2z_F3y_aa-2.0E0*1*I_ESP_G3xy_F3y_a;
    abcd[iGrid*900+842] = 4.0E0*I_ESP_I3x3z_F3y_aa-2.0E0*1*I_ESP_G3xz_F3y_a-2.0E0*2*I_ESP_G3xz_F3y_a;
    abcd[iGrid*900+843] = 4.0E0*I_ESP_I2x2y2z_F3y_aa-2.0E0*1*I_ESP_G2x2y_F3y_a;
    abcd[iGrid*900+844] = 4.0E0*I_ESP_I2xy3z_F3y_aa-2.0E0*1*I_ESP_G2xyz_F3y_a-2.0E0*2*I_ESP_G2xyz_F3y_a;
    abcd[iGrid*900+845] = 4.0E0*I_ESP_I2x4z_F3y_aa-2.0E0*2*I_ESP_G2x2z_F3y_a-2.0E0*3*I_ESP_G2x2z_F3y_a+2*1*I_ESP_D2x_F3y;
    abcd[iGrid*900+846] = 4.0E0*I_ESP_Ix3y2z_F3y_aa-2.0E0*1*I_ESP_Gx3y_F3y_a;
    abcd[iGrid*900+847] = 4.0E0*I_ESP_Ix2y3z_F3y_aa-2.0E0*1*I_ESP_Gx2yz_F3y_a-2.0E0*2*I_ESP_Gx2yz_F3y_a;
    abcd[iGrid*900+848] = 4.0E0*I_ESP_Ixy4z_F3y_aa-2.0E0*2*I_ESP_Gxy2z_F3y_a-2.0E0*3*I_ESP_Gxy2z_F3y_a+2*1*I_ESP_Dxy_F3y;
    abcd[iGrid*900+849] = 4.0E0*I_ESP_Ix5z_F3y_aa-2.0E0*3*I_ESP_Gx3z_F3y_a-2.0E0*4*I_ESP_Gx3z_F3y_a+3*2*I_ESP_Dxz_F3y;
    abcd[iGrid*900+850] = 4.0E0*I_ESP_I4y2z_F3y_aa-2.0E0*1*I_ESP_G4y_F3y_a;
    abcd[iGrid*900+851] = 4.0E0*I_ESP_I3y3z_F3y_aa-2.0E0*1*I_ESP_G3yz_F3y_a-2.0E0*2*I_ESP_G3yz_F3y_a;
    abcd[iGrid*900+852] = 4.0E0*I_ESP_I2y4z_F3y_aa-2.0E0*2*I_ESP_G2y2z_F3y_a-2.0E0*3*I_ESP_G2y2z_F3y_a+2*1*I_ESP_D2y_F3y;
    abcd[iGrid*900+853] = 4.0E0*I_ESP_Iy5z_F3y_aa-2.0E0*3*I_ESP_Gy3z_F3y_a-2.0E0*4*I_ESP_Gy3z_F3y_a+3*2*I_ESP_Dyz_F3y;
    abcd[iGrid*900+854] = 4.0E0*I_ESP_I6z_F3y_aa-2.0E0*4*I_ESP_G4z_F3y_a-2.0E0*5*I_ESP_G4z_F3y_a+4*3*I_ESP_D2z_F3y;
    abcd[iGrid*900+855] = 4.0E0*I_ESP_I4x2z_F2yz_aa-2.0E0*1*I_ESP_G4x_F2yz_a;
    abcd[iGrid*900+856] = 4.0E0*I_ESP_I3xy2z_F2yz_aa-2.0E0*1*I_ESP_G3xy_F2yz_a;
    abcd[iGrid*900+857] = 4.0E0*I_ESP_I3x3z_F2yz_aa-2.0E0*1*I_ESP_G3xz_F2yz_a-2.0E0*2*I_ESP_G3xz_F2yz_a;
    abcd[iGrid*900+858] = 4.0E0*I_ESP_I2x2y2z_F2yz_aa-2.0E0*1*I_ESP_G2x2y_F2yz_a;
    abcd[iGrid*900+859] = 4.0E0*I_ESP_I2xy3z_F2yz_aa-2.0E0*1*I_ESP_G2xyz_F2yz_a-2.0E0*2*I_ESP_G2xyz_F2yz_a;
    abcd[iGrid*900+860] = 4.0E0*I_ESP_I2x4z_F2yz_aa-2.0E0*2*I_ESP_G2x2z_F2yz_a-2.0E0*3*I_ESP_G2x2z_F2yz_a+2*1*I_ESP_D2x_F2yz;
    abcd[iGrid*900+861] = 4.0E0*I_ESP_Ix3y2z_F2yz_aa-2.0E0*1*I_ESP_Gx3y_F2yz_a;
    abcd[iGrid*900+862] = 4.0E0*I_ESP_Ix2y3z_F2yz_aa-2.0E0*1*I_ESP_Gx2yz_F2yz_a-2.0E0*2*I_ESP_Gx2yz_F2yz_a;
    abcd[iGrid*900+863] = 4.0E0*I_ESP_Ixy4z_F2yz_aa-2.0E0*2*I_ESP_Gxy2z_F2yz_a-2.0E0*3*I_ESP_Gxy2z_F2yz_a+2*1*I_ESP_Dxy_F2yz;
    abcd[iGrid*900+864] = 4.0E0*I_ESP_Ix5z_F2yz_aa-2.0E0*3*I_ESP_Gx3z_F2yz_a-2.0E0*4*I_ESP_Gx3z_F2yz_a+3*2*I_ESP_Dxz_F2yz;
    abcd[iGrid*900+865] = 4.0E0*I_ESP_I4y2z_F2yz_aa-2.0E0*1*I_ESP_G4y_F2yz_a;
    abcd[iGrid*900+866] = 4.0E0*I_ESP_I3y3z_F2yz_aa-2.0E0*1*I_ESP_G3yz_F2yz_a-2.0E0*2*I_ESP_G3yz_F2yz_a;
    abcd[iGrid*900+867] = 4.0E0*I_ESP_I2y4z_F2yz_aa-2.0E0*2*I_ESP_G2y2z_F2yz_a-2.0E0*3*I_ESP_G2y2z_F2yz_a+2*1*I_ESP_D2y_F2yz;
    abcd[iGrid*900+868] = 4.0E0*I_ESP_Iy5z_F2yz_aa-2.0E0*3*I_ESP_Gy3z_F2yz_a-2.0E0*4*I_ESP_Gy3z_F2yz_a+3*2*I_ESP_Dyz_F2yz;
    abcd[iGrid*900+869] = 4.0E0*I_ESP_I6z_F2yz_aa-2.0E0*4*I_ESP_G4z_F2yz_a-2.0E0*5*I_ESP_G4z_F2yz_a+4*3*I_ESP_D2z_F2yz;
    abcd[iGrid*900+870] = 4.0E0*I_ESP_I4x2z_Fy2z_aa-2.0E0*1*I_ESP_G4x_Fy2z_a;
    abcd[iGrid*900+871] = 4.0E0*I_ESP_I3xy2z_Fy2z_aa-2.0E0*1*I_ESP_G3xy_Fy2z_a;
    abcd[iGrid*900+872] = 4.0E0*I_ESP_I3x3z_Fy2z_aa-2.0E0*1*I_ESP_G3xz_Fy2z_a-2.0E0*2*I_ESP_G3xz_Fy2z_a;
    abcd[iGrid*900+873] = 4.0E0*I_ESP_I2x2y2z_Fy2z_aa-2.0E0*1*I_ESP_G2x2y_Fy2z_a;
    abcd[iGrid*900+874] = 4.0E0*I_ESP_I2xy3z_Fy2z_aa-2.0E0*1*I_ESP_G2xyz_Fy2z_a-2.0E0*2*I_ESP_G2xyz_Fy2z_a;
    abcd[iGrid*900+875] = 4.0E0*I_ESP_I2x4z_Fy2z_aa-2.0E0*2*I_ESP_G2x2z_Fy2z_a-2.0E0*3*I_ESP_G2x2z_Fy2z_a+2*1*I_ESP_D2x_Fy2z;
    abcd[iGrid*900+876] = 4.0E0*I_ESP_Ix3y2z_Fy2z_aa-2.0E0*1*I_ESP_Gx3y_Fy2z_a;
    abcd[iGrid*900+877] = 4.0E0*I_ESP_Ix2y3z_Fy2z_aa-2.0E0*1*I_ESP_Gx2yz_Fy2z_a-2.0E0*2*I_ESP_Gx2yz_Fy2z_a;
    abcd[iGrid*900+878] = 4.0E0*I_ESP_Ixy4z_Fy2z_aa-2.0E0*2*I_ESP_Gxy2z_Fy2z_a-2.0E0*3*I_ESP_Gxy2z_Fy2z_a+2*1*I_ESP_Dxy_Fy2z;
    abcd[iGrid*900+879] = 4.0E0*I_ESP_Ix5z_Fy2z_aa-2.0E0*3*I_ESP_Gx3z_Fy2z_a-2.0E0*4*I_ESP_Gx3z_Fy2z_a+3*2*I_ESP_Dxz_Fy2z;
    abcd[iGrid*900+880] = 4.0E0*I_ESP_I4y2z_Fy2z_aa-2.0E0*1*I_ESP_G4y_Fy2z_a;
    abcd[iGrid*900+881] = 4.0E0*I_ESP_I3y3z_Fy2z_aa-2.0E0*1*I_ESP_G3yz_Fy2z_a-2.0E0*2*I_ESP_G3yz_Fy2z_a;
    abcd[iGrid*900+882] = 4.0E0*I_ESP_I2y4z_Fy2z_aa-2.0E0*2*I_ESP_G2y2z_Fy2z_a-2.0E0*3*I_ESP_G2y2z_Fy2z_a+2*1*I_ESP_D2y_Fy2z;
    abcd[iGrid*900+883] = 4.0E0*I_ESP_Iy5z_Fy2z_aa-2.0E0*3*I_ESP_Gy3z_Fy2z_a-2.0E0*4*I_ESP_Gy3z_Fy2z_a+3*2*I_ESP_Dyz_Fy2z;
    abcd[iGrid*900+884] = 4.0E0*I_ESP_I6z_Fy2z_aa-2.0E0*4*I_ESP_G4z_Fy2z_a-2.0E0*5*I_ESP_G4z_Fy2z_a+4*3*I_ESP_D2z_Fy2z;
    abcd[iGrid*900+885] = 4.0E0*I_ESP_I4x2z_F3z_aa-2.0E0*1*I_ESP_G4x_F3z_a;
    abcd[iGrid*900+886] = 4.0E0*I_ESP_I3xy2z_F3z_aa-2.0E0*1*I_ESP_G3xy_F3z_a;
    abcd[iGrid*900+887] = 4.0E0*I_ESP_I3x3z_F3z_aa-2.0E0*1*I_ESP_G3xz_F3z_a-2.0E0*2*I_ESP_G3xz_F3z_a;
    abcd[iGrid*900+888] = 4.0E0*I_ESP_I2x2y2z_F3z_aa-2.0E0*1*I_ESP_G2x2y_F3z_a;
    abcd[iGrid*900+889] = 4.0E0*I_ESP_I2xy3z_F3z_aa-2.0E0*1*I_ESP_G2xyz_F3z_a-2.0E0*2*I_ESP_G2xyz_F3z_a;
    abcd[iGrid*900+890] = 4.0E0*I_ESP_I2x4z_F3z_aa-2.0E0*2*I_ESP_G2x2z_F3z_a-2.0E0*3*I_ESP_G2x2z_F3z_a+2*1*I_ESP_D2x_F3z;
    abcd[iGrid*900+891] = 4.0E0*I_ESP_Ix3y2z_F3z_aa-2.0E0*1*I_ESP_Gx3y_F3z_a;
    abcd[iGrid*900+892] = 4.0E0*I_ESP_Ix2y3z_F3z_aa-2.0E0*1*I_ESP_Gx2yz_F3z_a-2.0E0*2*I_ESP_Gx2yz_F3z_a;
    abcd[iGrid*900+893] = 4.0E0*I_ESP_Ixy4z_F3z_aa-2.0E0*2*I_ESP_Gxy2z_F3z_a-2.0E0*3*I_ESP_Gxy2z_F3z_a+2*1*I_ESP_Dxy_F3z;
    abcd[iGrid*900+894] = 4.0E0*I_ESP_Ix5z_F3z_aa-2.0E0*3*I_ESP_Gx3z_F3z_a-2.0E0*4*I_ESP_Gx3z_F3z_a+3*2*I_ESP_Dxz_F3z;
    abcd[iGrid*900+895] = 4.0E0*I_ESP_I4y2z_F3z_aa-2.0E0*1*I_ESP_G4y_F3z_a;
    abcd[iGrid*900+896] = 4.0E0*I_ESP_I3y3z_F3z_aa-2.0E0*1*I_ESP_G3yz_F3z_a-2.0E0*2*I_ESP_G3yz_F3z_a;
    abcd[iGrid*900+897] = 4.0E0*I_ESP_I2y4z_F3z_aa-2.0E0*2*I_ESP_G2y2z_F3z_a-2.0E0*3*I_ESP_G2y2z_F3z_a+2*1*I_ESP_D2y_F3z;
    abcd[iGrid*900+898] = 4.0E0*I_ESP_Iy5z_F3z_aa-2.0E0*3*I_ESP_Gy3z_F3z_a-2.0E0*4*I_ESP_Gy3z_F3z_a+3*2*I_ESP_Dyz_F3z;
    abcd[iGrid*900+899] = 4.0E0*I_ESP_I6z_F3z_aa-2.0E0*4*I_ESP_G4z_F3z_a-2.0E0*5*I_ESP_G4z_F3z_a+4*3*I_ESP_D2z_F3z;
  }
}
