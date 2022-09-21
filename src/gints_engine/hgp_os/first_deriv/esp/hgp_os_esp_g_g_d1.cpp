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

void hgp_os_esp_g_g_d1(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_M9x_S_a = 0.0E0;
    Double I_ESP_M8xy_S_a = 0.0E0;
    Double I_ESP_M8xz_S_a = 0.0E0;
    Double I_ESP_M7x2y_S_a = 0.0E0;
    Double I_ESP_M7xyz_S_a = 0.0E0;
    Double I_ESP_M7x2z_S_a = 0.0E0;
    Double I_ESP_M6x3y_S_a = 0.0E0;
    Double I_ESP_M6x2yz_S_a = 0.0E0;
    Double I_ESP_M6xy2z_S_a = 0.0E0;
    Double I_ESP_M6x3z_S_a = 0.0E0;
    Double I_ESP_M5x4y_S_a = 0.0E0;
    Double I_ESP_M5x3yz_S_a = 0.0E0;
    Double I_ESP_M5x2y2z_S_a = 0.0E0;
    Double I_ESP_M5xy3z_S_a = 0.0E0;
    Double I_ESP_M5x4z_S_a = 0.0E0;
    Double I_ESP_M4x5y_S_a = 0.0E0;
    Double I_ESP_M4x4yz_S_a = 0.0E0;
    Double I_ESP_M4x3y2z_S_a = 0.0E0;
    Double I_ESP_M4x2y3z_S_a = 0.0E0;
    Double I_ESP_M4xy4z_S_a = 0.0E0;
    Double I_ESP_M4x5z_S_a = 0.0E0;
    Double I_ESP_M3x6y_S_a = 0.0E0;
    Double I_ESP_M3x5yz_S_a = 0.0E0;
    Double I_ESP_M3x4y2z_S_a = 0.0E0;
    Double I_ESP_M3x3y3z_S_a = 0.0E0;
    Double I_ESP_M3x2y4z_S_a = 0.0E0;
    Double I_ESP_M3xy5z_S_a = 0.0E0;
    Double I_ESP_M3x6z_S_a = 0.0E0;
    Double I_ESP_M2x7y_S_a = 0.0E0;
    Double I_ESP_M2x6yz_S_a = 0.0E0;
    Double I_ESP_M2x5y2z_S_a = 0.0E0;
    Double I_ESP_M2x4y3z_S_a = 0.0E0;
    Double I_ESP_M2x3y4z_S_a = 0.0E0;
    Double I_ESP_M2x2y5z_S_a = 0.0E0;
    Double I_ESP_M2xy6z_S_a = 0.0E0;
    Double I_ESP_M2x7z_S_a = 0.0E0;
    Double I_ESP_Mx8y_S_a = 0.0E0;
    Double I_ESP_Mx7yz_S_a = 0.0E0;
    Double I_ESP_Mx6y2z_S_a = 0.0E0;
    Double I_ESP_Mx5y3z_S_a = 0.0E0;
    Double I_ESP_Mx4y4z_S_a = 0.0E0;
    Double I_ESP_Mx3y5z_S_a = 0.0E0;
    Double I_ESP_Mx2y6z_S_a = 0.0E0;
    Double I_ESP_Mxy7z_S_a = 0.0E0;
    Double I_ESP_Mx8z_S_a = 0.0E0;
    Double I_ESP_M9y_S_a = 0.0E0;
    Double I_ESP_M8yz_S_a = 0.0E0;
    Double I_ESP_M7y2z_S_a = 0.0E0;
    Double I_ESP_M6y3z_S_a = 0.0E0;
    Double I_ESP_M5y4z_S_a = 0.0E0;
    Double I_ESP_M4y5z_S_a = 0.0E0;
    Double I_ESP_M3y6z_S_a = 0.0E0;
    Double I_ESP_M2y7z_S_a = 0.0E0;
    Double I_ESP_My8z_S_a = 0.0E0;
    Double I_ESP_M9z_S_a = 0.0E0;
    Double I_ESP_L8x_S_a = 0.0E0;
    Double I_ESP_L7xy_S_a = 0.0E0;
    Double I_ESP_L7xz_S_a = 0.0E0;
    Double I_ESP_L6x2y_S_a = 0.0E0;
    Double I_ESP_L6xyz_S_a = 0.0E0;
    Double I_ESP_L6x2z_S_a = 0.0E0;
    Double I_ESP_L5x3y_S_a = 0.0E0;
    Double I_ESP_L5x2yz_S_a = 0.0E0;
    Double I_ESP_L5xy2z_S_a = 0.0E0;
    Double I_ESP_L5x3z_S_a = 0.0E0;
    Double I_ESP_L4x4y_S_a = 0.0E0;
    Double I_ESP_L4x3yz_S_a = 0.0E0;
    Double I_ESP_L4x2y2z_S_a = 0.0E0;
    Double I_ESP_L4xy3z_S_a = 0.0E0;
    Double I_ESP_L4x4z_S_a = 0.0E0;
    Double I_ESP_L3x5y_S_a = 0.0E0;
    Double I_ESP_L3x4yz_S_a = 0.0E0;
    Double I_ESP_L3x3y2z_S_a = 0.0E0;
    Double I_ESP_L3x2y3z_S_a = 0.0E0;
    Double I_ESP_L3xy4z_S_a = 0.0E0;
    Double I_ESP_L3x5z_S_a = 0.0E0;
    Double I_ESP_L2x6y_S_a = 0.0E0;
    Double I_ESP_L2x5yz_S_a = 0.0E0;
    Double I_ESP_L2x4y2z_S_a = 0.0E0;
    Double I_ESP_L2x3y3z_S_a = 0.0E0;
    Double I_ESP_L2x2y4z_S_a = 0.0E0;
    Double I_ESP_L2xy5z_S_a = 0.0E0;
    Double I_ESP_L2x6z_S_a = 0.0E0;
    Double I_ESP_Lx7y_S_a = 0.0E0;
    Double I_ESP_Lx6yz_S_a = 0.0E0;
    Double I_ESP_Lx5y2z_S_a = 0.0E0;
    Double I_ESP_Lx4y3z_S_a = 0.0E0;
    Double I_ESP_Lx3y4z_S_a = 0.0E0;
    Double I_ESP_Lx2y5z_S_a = 0.0E0;
    Double I_ESP_Lxy6z_S_a = 0.0E0;
    Double I_ESP_Lx7z_S_a = 0.0E0;
    Double I_ESP_L8y_S_a = 0.0E0;
    Double I_ESP_L7yz_S_a = 0.0E0;
    Double I_ESP_L6y2z_S_a = 0.0E0;
    Double I_ESP_L5y3z_S_a = 0.0E0;
    Double I_ESP_L4y4z_S_a = 0.0E0;
    Double I_ESP_L3y5z_S_a = 0.0E0;
    Double I_ESP_L2y6z_S_a = 0.0E0;
    Double I_ESP_Ly7z_S_a = 0.0E0;
    Double I_ESP_L8z_S_a = 0.0E0;
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
    Double I_ESP_K7x_S = 0.0E0;
    Double I_ESP_K6xy_S = 0.0E0;
    Double I_ESP_K6xz_S = 0.0E0;
    Double I_ESP_K5x2y_S = 0.0E0;
    Double I_ESP_K5xyz_S = 0.0E0;
    Double I_ESP_K5x2z_S = 0.0E0;
    Double I_ESP_K4x3y_S = 0.0E0;
    Double I_ESP_K4x2yz_S = 0.0E0;
    Double I_ESP_K4xy2z_S = 0.0E0;
    Double I_ESP_K4x3z_S = 0.0E0;
    Double I_ESP_K3x4y_S = 0.0E0;
    Double I_ESP_K3x3yz_S = 0.0E0;
    Double I_ESP_K3x2y2z_S = 0.0E0;
    Double I_ESP_K3xy3z_S = 0.0E0;
    Double I_ESP_K3x4z_S = 0.0E0;
    Double I_ESP_K2x5y_S = 0.0E0;
    Double I_ESP_K2x4yz_S = 0.0E0;
    Double I_ESP_K2x3y2z_S = 0.0E0;
    Double I_ESP_K2x2y3z_S = 0.0E0;
    Double I_ESP_K2xy4z_S = 0.0E0;
    Double I_ESP_K2x5z_S = 0.0E0;
    Double I_ESP_Kx6y_S = 0.0E0;
    Double I_ESP_Kx5yz_S = 0.0E0;
    Double I_ESP_Kx4y2z_S = 0.0E0;
    Double I_ESP_Kx3y3z_S = 0.0E0;
    Double I_ESP_Kx2y4z_S = 0.0E0;
    Double I_ESP_Kxy5z_S = 0.0E0;
    Double I_ESP_Kx6z_S = 0.0E0;
    Double I_ESP_K7y_S = 0.0E0;
    Double I_ESP_K6yz_S = 0.0E0;
    Double I_ESP_K5y2z_S = 0.0E0;
    Double I_ESP_K4y3z_S = 0.0E0;
    Double I_ESP_K3y4z_S = 0.0E0;
    Double I_ESP_K2y5z_S = 0.0E0;
    Double I_ESP_Ky6z_S = 0.0E0;
    Double I_ESP_K7z_S = 0.0E0;
    Double I_ESP_I6x_S = 0.0E0;
    Double I_ESP_I5xy_S = 0.0E0;
    Double I_ESP_I5xz_S = 0.0E0;
    Double I_ESP_I4x2y_S = 0.0E0;
    Double I_ESP_I4xyz_S = 0.0E0;
    Double I_ESP_I4x2z_S = 0.0E0;
    Double I_ESP_I3x3y_S = 0.0E0;
    Double I_ESP_I3x2yz_S = 0.0E0;
    Double I_ESP_I3xy2z_S = 0.0E0;
    Double I_ESP_I3x3z_S = 0.0E0;
    Double I_ESP_I2x4y_S = 0.0E0;
    Double I_ESP_I2x3yz_S = 0.0E0;
    Double I_ESP_I2x2y2z_S = 0.0E0;
    Double I_ESP_I2xy3z_S = 0.0E0;
    Double I_ESP_I2x4z_S = 0.0E0;
    Double I_ESP_Ix5y_S = 0.0E0;
    Double I_ESP_Ix4yz_S = 0.0E0;
    Double I_ESP_Ix3y2z_S = 0.0E0;
    Double I_ESP_Ix2y3z_S = 0.0E0;
    Double I_ESP_Ixy4z_S = 0.0E0;
    Double I_ESP_Ix5z_S = 0.0E0;
    Double I_ESP_I6y_S = 0.0E0;
    Double I_ESP_I5yz_S = 0.0E0;
    Double I_ESP_I4y2z_S = 0.0E0;
    Double I_ESP_I3y3z_S = 0.0E0;
    Double I_ESP_I2y4z_S = 0.0E0;
    Double I_ESP_Iy5z_S = 0.0E0;
    Double I_ESP_I6z_S = 0.0E0;
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
         * shell quartet name: SQ_ESP_M_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_M_S_a_coefs = alpha;
        I_ESP_M9x_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M9x_S_vrr;
        I_ESP_M8xy_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M8xy_S_vrr;
        I_ESP_M8xz_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M8xz_S_vrr;
        I_ESP_M7x2y_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M7x2y_S_vrr;
        I_ESP_M7xyz_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M7xyz_S_vrr;
        I_ESP_M7x2z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M7x2z_S_vrr;
        I_ESP_M6x3y_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M6x3y_S_vrr;
        I_ESP_M6x2yz_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M6x2yz_S_vrr;
        I_ESP_M6xy2z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M6xy2z_S_vrr;
        I_ESP_M6x3z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M6x3z_S_vrr;
        I_ESP_M5x4y_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M5x4y_S_vrr;
        I_ESP_M5x3yz_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M5x3yz_S_vrr;
        I_ESP_M5x2y2z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M5x2y2z_S_vrr;
        I_ESP_M5xy3z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M5xy3z_S_vrr;
        I_ESP_M5x4z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M5x4z_S_vrr;
        I_ESP_M4x5y_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M4x5y_S_vrr;
        I_ESP_M4x4yz_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M4x4yz_S_vrr;
        I_ESP_M4x3y2z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M4x3y2z_S_vrr;
        I_ESP_M4x2y3z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M4x2y3z_S_vrr;
        I_ESP_M4xy4z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M4xy4z_S_vrr;
        I_ESP_M4x5z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M4x5z_S_vrr;
        I_ESP_M3x6y_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M3x6y_S_vrr;
        I_ESP_M3x5yz_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M3x5yz_S_vrr;
        I_ESP_M3x4y2z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M3x4y2z_S_vrr;
        I_ESP_M3x3y3z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M3x3y3z_S_vrr;
        I_ESP_M3x2y4z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M3x2y4z_S_vrr;
        I_ESP_M3xy5z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M3xy5z_S_vrr;
        I_ESP_M3x6z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M3x6z_S_vrr;
        I_ESP_M2x7y_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M2x7y_S_vrr;
        I_ESP_M2x6yz_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M2x6yz_S_vrr;
        I_ESP_M2x5y2z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M2x5y2z_S_vrr;
        I_ESP_M2x4y3z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M2x4y3z_S_vrr;
        I_ESP_M2x3y4z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M2x3y4z_S_vrr;
        I_ESP_M2x2y5z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M2x2y5z_S_vrr;
        I_ESP_M2xy6z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M2xy6z_S_vrr;
        I_ESP_M2x7z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M2x7z_S_vrr;
        I_ESP_Mx8y_S_a += SQ_ESP_M_S_a_coefs*I_ESP_Mx8y_S_vrr;
        I_ESP_Mx7yz_S_a += SQ_ESP_M_S_a_coefs*I_ESP_Mx7yz_S_vrr;
        I_ESP_Mx6y2z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_Mx6y2z_S_vrr;
        I_ESP_Mx5y3z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_Mx5y3z_S_vrr;
        I_ESP_Mx4y4z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_Mx4y4z_S_vrr;
        I_ESP_Mx3y5z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_Mx3y5z_S_vrr;
        I_ESP_Mx2y6z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_Mx2y6z_S_vrr;
        I_ESP_Mxy7z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_Mxy7z_S_vrr;
        I_ESP_Mx8z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_Mx8z_S_vrr;
        I_ESP_M9y_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M9y_S_vrr;
        I_ESP_M8yz_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M8yz_S_vrr;
        I_ESP_M7y2z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M7y2z_S_vrr;
        I_ESP_M6y3z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M6y3z_S_vrr;
        I_ESP_M5y4z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M5y4z_S_vrr;
        I_ESP_M4y5z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M4y5z_S_vrr;
        I_ESP_M3y6z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M3y6z_S_vrr;
        I_ESP_M2y7z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M2y7z_S_vrr;
        I_ESP_My8z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_My8z_S_vrr;
        I_ESP_M9z_S_a += SQ_ESP_M_S_a_coefs*I_ESP_M9z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_L_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_L_S_a_coefs = alpha;
        I_ESP_L8x_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L8x_S_vrr;
        I_ESP_L7xy_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L7xy_S_vrr;
        I_ESP_L7xz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L7xz_S_vrr;
        I_ESP_L6x2y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L6x2y_S_vrr;
        I_ESP_L6xyz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L6xyz_S_vrr;
        I_ESP_L6x2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L6x2z_S_vrr;
        I_ESP_L5x3y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5x3y_S_vrr;
        I_ESP_L5x2yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5x2yz_S_vrr;
        I_ESP_L5xy2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5xy2z_S_vrr;
        I_ESP_L5x3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5x3z_S_vrr;
        I_ESP_L4x4y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4x4y_S_vrr;
        I_ESP_L4x3yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4x3yz_S_vrr;
        I_ESP_L4x2y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4x2y2z_S_vrr;
        I_ESP_L4xy3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4xy3z_S_vrr;
        I_ESP_L4x4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4x4z_S_vrr;
        I_ESP_L3x5y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x5y_S_vrr;
        I_ESP_L3x4yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x4yz_S_vrr;
        I_ESP_L3x3y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x3y2z_S_vrr;
        I_ESP_L3x2y3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x2y3z_S_vrr;
        I_ESP_L3xy4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3xy4z_S_vrr;
        I_ESP_L3x5z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3x5z_S_vrr;
        I_ESP_L2x6y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x6y_S_vrr;
        I_ESP_L2x5yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x5yz_S_vrr;
        I_ESP_L2x4y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x4y2z_S_vrr;
        I_ESP_L2x3y3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x3y3z_S_vrr;
        I_ESP_L2x2y4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x2y4z_S_vrr;
        I_ESP_L2xy5z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2xy5z_S_vrr;
        I_ESP_L2x6z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2x6z_S_vrr;
        I_ESP_Lx7y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx7y_S_vrr;
        I_ESP_Lx6yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx6yz_S_vrr;
        I_ESP_Lx5y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx5y2z_S_vrr;
        I_ESP_Lx4y3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx4y3z_S_vrr;
        I_ESP_Lx3y4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx3y4z_S_vrr;
        I_ESP_Lx2y5z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx2y5z_S_vrr;
        I_ESP_Lxy6z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lxy6z_S_vrr;
        I_ESP_Lx7z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Lx7z_S_vrr;
        I_ESP_L8y_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L8y_S_vrr;
        I_ESP_L7yz_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L7yz_S_vrr;
        I_ESP_L6y2z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L6y2z_S_vrr;
        I_ESP_L5y3z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L5y3z_S_vrr;
        I_ESP_L4y4z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L4y4z_S_vrr;
        I_ESP_L3y5z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L3y5z_S_vrr;
        I_ESP_L2y6z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L2y6z_S_vrr;
        I_ESP_Ly7z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_Ly7z_S_vrr;
        I_ESP_L8z_S_a += SQ_ESP_L_S_a_coefs*I_ESP_L8z_S_vrr;

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
         * shell quartet name: SQ_ESP_K_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_K7x_S += I_ESP_K7x_S_vrr;
        I_ESP_K6xy_S += I_ESP_K6xy_S_vrr;
        I_ESP_K6xz_S += I_ESP_K6xz_S_vrr;
        I_ESP_K5x2y_S += I_ESP_K5x2y_S_vrr;
        I_ESP_K5xyz_S += I_ESP_K5xyz_S_vrr;
        I_ESP_K5x2z_S += I_ESP_K5x2z_S_vrr;
        I_ESP_K4x3y_S += I_ESP_K4x3y_S_vrr;
        I_ESP_K4x2yz_S += I_ESP_K4x2yz_S_vrr;
        I_ESP_K4xy2z_S += I_ESP_K4xy2z_S_vrr;
        I_ESP_K4x3z_S += I_ESP_K4x3z_S_vrr;
        I_ESP_K3x4y_S += I_ESP_K3x4y_S_vrr;
        I_ESP_K3x3yz_S += I_ESP_K3x3yz_S_vrr;
        I_ESP_K3x2y2z_S += I_ESP_K3x2y2z_S_vrr;
        I_ESP_K3xy3z_S += I_ESP_K3xy3z_S_vrr;
        I_ESP_K3x4z_S += I_ESP_K3x4z_S_vrr;
        I_ESP_K2x5y_S += I_ESP_K2x5y_S_vrr;
        I_ESP_K2x4yz_S += I_ESP_K2x4yz_S_vrr;
        I_ESP_K2x3y2z_S += I_ESP_K2x3y2z_S_vrr;
        I_ESP_K2x2y3z_S += I_ESP_K2x2y3z_S_vrr;
        I_ESP_K2xy4z_S += I_ESP_K2xy4z_S_vrr;
        I_ESP_K2x5z_S += I_ESP_K2x5z_S_vrr;
        I_ESP_Kx6y_S += I_ESP_Kx6y_S_vrr;
        I_ESP_Kx5yz_S += I_ESP_Kx5yz_S_vrr;
        I_ESP_Kx4y2z_S += I_ESP_Kx4y2z_S_vrr;
        I_ESP_Kx3y3z_S += I_ESP_Kx3y3z_S_vrr;
        I_ESP_Kx2y4z_S += I_ESP_Kx2y4z_S_vrr;
        I_ESP_Kxy5z_S += I_ESP_Kxy5z_S_vrr;
        I_ESP_Kx6z_S += I_ESP_Kx6z_S_vrr;
        I_ESP_K7y_S += I_ESP_K7y_S_vrr;
        I_ESP_K6yz_S += I_ESP_K6yz_S_vrr;
        I_ESP_K5y2z_S += I_ESP_K5y2z_S_vrr;
        I_ESP_K4y3z_S += I_ESP_K4y3z_S_vrr;
        I_ESP_K3y4z_S += I_ESP_K3y4z_S_vrr;
        I_ESP_K2y5z_S += I_ESP_K2y5z_S_vrr;
        I_ESP_Ky6z_S += I_ESP_Ky6z_S_vrr;
        I_ESP_K7z_S += I_ESP_K7z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_I_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_I6x_S += I_ESP_I6x_S_vrr;
        I_ESP_I5xy_S += I_ESP_I5xy_S_vrr;
        I_ESP_I5xz_S += I_ESP_I5xz_S_vrr;
        I_ESP_I4x2y_S += I_ESP_I4x2y_S_vrr;
        I_ESP_I4xyz_S += I_ESP_I4xyz_S_vrr;
        I_ESP_I4x2z_S += I_ESP_I4x2z_S_vrr;
        I_ESP_I3x3y_S += I_ESP_I3x3y_S_vrr;
        I_ESP_I3x2yz_S += I_ESP_I3x2yz_S_vrr;
        I_ESP_I3xy2z_S += I_ESP_I3xy2z_S_vrr;
        I_ESP_I3x3z_S += I_ESP_I3x3z_S_vrr;
        I_ESP_I2x4y_S += I_ESP_I2x4y_S_vrr;
        I_ESP_I2x3yz_S += I_ESP_I2x3yz_S_vrr;
        I_ESP_I2x2y2z_S += I_ESP_I2x2y2z_S_vrr;
        I_ESP_I2xy3z_S += I_ESP_I2xy3z_S_vrr;
        I_ESP_I2x4z_S += I_ESP_I2x4z_S_vrr;
        I_ESP_Ix5y_S += I_ESP_Ix5y_S_vrr;
        I_ESP_Ix4yz_S += I_ESP_Ix4yz_S_vrr;
        I_ESP_Ix3y2z_S += I_ESP_Ix3y2z_S_vrr;
        I_ESP_Ix2y3z_S += I_ESP_Ix2y3z_S_vrr;
        I_ESP_Ixy4z_S += I_ESP_Ixy4z_S_vrr;
        I_ESP_Ix5z_S += I_ESP_Ix5z_S_vrr;
        I_ESP_I6y_S += I_ESP_I6y_S_vrr;
        I_ESP_I5yz_S += I_ESP_I5yz_S_vrr;
        I_ESP_I4y2z_S += I_ESP_I4y2z_S_vrr;
        I_ESP_I3y3z_S += I_ESP_I3y3z_S_vrr;
        I_ESP_I2y4z_S += I_ESP_I2y4z_S_vrr;
        I_ESP_Iy5z_S += I_ESP_Iy5z_S_vrr;
        I_ESP_I6z_S += I_ESP_I6z_S_vrr;

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
     * totally 0 integrals are omitted 
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
    Double I_ESP_G4x_Py = I_ESP_H4xy_S+ABY*I_ESP_G4x_S;
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
    Double I_ESP_G4x_Pz = I_ESP_H4xz_S+ABZ*I_ESP_G4x_S;
    Double I_ESP_G3xy_Pz = I_ESP_H3xyz_S+ABZ*I_ESP_G3xy_S;
    Double I_ESP_G3xz_Pz = I_ESP_H3x2z_S+ABZ*I_ESP_G3xz_S;
    Double I_ESP_G2x2y_Pz = I_ESP_H2x2yz_S+ABZ*I_ESP_G2x2y_S;
    Double I_ESP_G2xyz_Pz = I_ESP_H2xy2z_S+ABZ*I_ESP_G2xyz_S;
    Double I_ESP_G2x2z_Pz = I_ESP_H2x3z_S+ABZ*I_ESP_G2x2z_S;
    Double I_ESP_Gx3y_Pz = I_ESP_Hx3yz_S+ABZ*I_ESP_Gx3y_S;
    Double I_ESP_Gx2yz_Pz = I_ESP_Hx2y2z_S+ABZ*I_ESP_Gx2yz_S;
    Double I_ESP_Gxy2z_Pz = I_ESP_Hxy3z_S+ABZ*I_ESP_Gxy2z_S;
    Double I_ESP_Gx3z_Pz = I_ESP_Hx4z_S+ABZ*I_ESP_Gx3z_S;
    Double I_ESP_G4y_Pz = I_ESP_H4yz_S+ABZ*I_ESP_G4y_S;
    Double I_ESP_G3yz_Pz = I_ESP_H3y2z_S+ABZ*I_ESP_G3yz_S;
    Double I_ESP_G2y2z_Pz = I_ESP_H2y3z_S+ABZ*I_ESP_G2y2z_S;
    Double I_ESP_Gy3z_Pz = I_ESP_Hy4z_S+ABZ*I_ESP_Gy3z_S;
    Double I_ESP_G4z_Pz = I_ESP_H5z_S+ABZ*I_ESP_G4z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_F_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 30 integrals are omitted 
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
     * shell quartet name: SQ_ESP_H_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_S
     * RHS shell quartet name: SQ_ESP_H_S
     ************************************************************/
    Double I_ESP_H5x_Px = I_ESP_I6x_S+ABX*I_ESP_H5x_S;
    Double I_ESP_H4xy_Px = I_ESP_I5xy_S+ABX*I_ESP_H4xy_S;
    Double I_ESP_H4xz_Px = I_ESP_I5xz_S+ABX*I_ESP_H4xz_S;
    Double I_ESP_H3x2y_Px = I_ESP_I4x2y_S+ABX*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Px = I_ESP_I4xyz_S+ABX*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Px = I_ESP_I4x2z_S+ABX*I_ESP_H3x2z_S;
    Double I_ESP_H2x3y_Px = I_ESP_I3x3y_S+ABX*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Px = I_ESP_I3x2yz_S+ABX*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Px = I_ESP_I3xy2z_S+ABX*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Px = I_ESP_I3x3z_S+ABX*I_ESP_H2x3z_S;
    Double I_ESP_Hx4y_Px = I_ESP_I2x4y_S+ABX*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Px = I_ESP_I2x3yz_S+ABX*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Px = I_ESP_I2x2y2z_S+ABX*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Px = I_ESP_I2xy3z_S+ABX*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Px = I_ESP_I2x4z_S+ABX*I_ESP_Hx4z_S;
    Double I_ESP_H5y_Px = I_ESP_Ix5y_S+ABX*I_ESP_H5y_S;
    Double I_ESP_H4yz_Px = I_ESP_Ix4yz_S+ABX*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Px = I_ESP_Ix3y2z_S+ABX*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Px = I_ESP_Ix2y3z_S+ABX*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Px = I_ESP_Ixy4z_S+ABX*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Px = I_ESP_Ix5z_S+ABX*I_ESP_H5z_S;
    Double I_ESP_H4xy_Py = I_ESP_I4x2y_S+ABY*I_ESP_H4xy_S;
    Double I_ESP_H4xz_Py = I_ESP_I4xyz_S+ABY*I_ESP_H4xz_S;
    Double I_ESP_H3x2y_Py = I_ESP_I3x3y_S+ABY*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Py = I_ESP_I3x2yz_S+ABY*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Py = I_ESP_I3xy2z_S+ABY*I_ESP_H3x2z_S;
    Double I_ESP_H2x3y_Py = I_ESP_I2x4y_S+ABY*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Py = I_ESP_I2x3yz_S+ABY*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Py = I_ESP_I2x2y2z_S+ABY*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Py = I_ESP_I2xy3z_S+ABY*I_ESP_H2x3z_S;
    Double I_ESP_Hx4y_Py = I_ESP_Ix5y_S+ABY*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Py = I_ESP_Ix4yz_S+ABY*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Py = I_ESP_Ix3y2z_S+ABY*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Py = I_ESP_Ix2y3z_S+ABY*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Py = I_ESP_Ixy4z_S+ABY*I_ESP_Hx4z_S;
    Double I_ESP_H5y_Py = I_ESP_I6y_S+ABY*I_ESP_H5y_S;
    Double I_ESP_H4yz_Py = I_ESP_I5yz_S+ABY*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Py = I_ESP_I4y2z_S+ABY*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Py = I_ESP_I3y3z_S+ABY*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Py = I_ESP_I2y4z_S+ABY*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Py = I_ESP_Iy5z_S+ABY*I_ESP_H5z_S;
    Double I_ESP_H4xy_Pz = I_ESP_I4xyz_S+ABZ*I_ESP_H4xy_S;
    Double I_ESP_H4xz_Pz = I_ESP_I4x2z_S+ABZ*I_ESP_H4xz_S;
    Double I_ESP_H3x2y_Pz = I_ESP_I3x2yz_S+ABZ*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Pz = I_ESP_I3xy2z_S+ABZ*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Pz = I_ESP_I3x3z_S+ABZ*I_ESP_H3x2z_S;
    Double I_ESP_H2x3y_Pz = I_ESP_I2x3yz_S+ABZ*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Pz = I_ESP_I2x2y2z_S+ABZ*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Pz = I_ESP_I2xy3z_S+ABZ*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Pz = I_ESP_I2x4z_S+ABZ*I_ESP_H2x3z_S;
    Double I_ESP_Hx4y_Pz = I_ESP_Ix4yz_S+ABZ*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Pz = I_ESP_Ix3y2z_S+ABZ*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Pz = I_ESP_Ix2y3z_S+ABZ*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Pz = I_ESP_Ixy4z_S+ABZ*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Pz = I_ESP_Ix5z_S+ABZ*I_ESP_Hx4z_S;
    Double I_ESP_H4yz_Pz = I_ESP_I4y2z_S+ABZ*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Pz = I_ESP_I3y3z_S+ABZ*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Pz = I_ESP_I2y4z_S+ABZ*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Pz = I_ESP_Iy5z_S+ABZ*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Pz = I_ESP_I6z_S+ABZ*I_ESP_H5z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 45 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_P
     * RHS shell quartet name: SQ_ESP_G_P
     ************************************************************/
    Double I_ESP_G4x_D2x = I_ESP_H5x_Px+ABX*I_ESP_G4x_Px;
    Double I_ESP_G3xy_D2x = I_ESP_H4xy_Px+ABX*I_ESP_G3xy_Px;
    Double I_ESP_G3xz_D2x = I_ESP_H4xz_Px+ABX*I_ESP_G3xz_Px;
    Double I_ESP_G2x2y_D2x = I_ESP_H3x2y_Px+ABX*I_ESP_G2x2y_Px;
    Double I_ESP_G2xyz_D2x = I_ESP_H3xyz_Px+ABX*I_ESP_G2xyz_Px;
    Double I_ESP_G2x2z_D2x = I_ESP_H3x2z_Px+ABX*I_ESP_G2x2z_Px;
    Double I_ESP_Gx3y_D2x = I_ESP_H2x3y_Px+ABX*I_ESP_Gx3y_Px;
    Double I_ESP_Gx2yz_D2x = I_ESP_H2x2yz_Px+ABX*I_ESP_Gx2yz_Px;
    Double I_ESP_Gxy2z_D2x = I_ESP_H2xy2z_Px+ABX*I_ESP_Gxy2z_Px;
    Double I_ESP_Gx3z_D2x = I_ESP_H2x3z_Px+ABX*I_ESP_Gx3z_Px;
    Double I_ESP_G4y_D2x = I_ESP_Hx4y_Px+ABX*I_ESP_G4y_Px;
    Double I_ESP_G3yz_D2x = I_ESP_Hx3yz_Px+ABX*I_ESP_G3yz_Px;
    Double I_ESP_G2y2z_D2x = I_ESP_Hx2y2z_Px+ABX*I_ESP_G2y2z_Px;
    Double I_ESP_Gy3z_D2x = I_ESP_Hxy3z_Px+ABX*I_ESP_Gy3z_Px;
    Double I_ESP_G4z_D2x = I_ESP_Hx4z_Px+ABX*I_ESP_G4z_Px;
    Double I_ESP_G4x_D2y = I_ESP_H4xy_Py+ABY*I_ESP_G4x_Py;
    Double I_ESP_G3xy_D2y = I_ESP_H3x2y_Py+ABY*I_ESP_G3xy_Py;
    Double I_ESP_G3xz_D2y = I_ESP_H3xyz_Py+ABY*I_ESP_G3xz_Py;
    Double I_ESP_G2x2y_D2y = I_ESP_H2x3y_Py+ABY*I_ESP_G2x2y_Py;
    Double I_ESP_G2xyz_D2y = I_ESP_H2x2yz_Py+ABY*I_ESP_G2xyz_Py;
    Double I_ESP_G2x2z_D2y = I_ESP_H2xy2z_Py+ABY*I_ESP_G2x2z_Py;
    Double I_ESP_Gx3y_D2y = I_ESP_Hx4y_Py+ABY*I_ESP_Gx3y_Py;
    Double I_ESP_Gx2yz_D2y = I_ESP_Hx3yz_Py+ABY*I_ESP_Gx2yz_Py;
    Double I_ESP_Gxy2z_D2y = I_ESP_Hx2y2z_Py+ABY*I_ESP_Gxy2z_Py;
    Double I_ESP_Gx3z_D2y = I_ESP_Hxy3z_Py+ABY*I_ESP_Gx3z_Py;
    Double I_ESP_G4y_D2y = I_ESP_H5y_Py+ABY*I_ESP_G4y_Py;
    Double I_ESP_G3yz_D2y = I_ESP_H4yz_Py+ABY*I_ESP_G3yz_Py;
    Double I_ESP_G2y2z_D2y = I_ESP_H3y2z_Py+ABY*I_ESP_G2y2z_Py;
    Double I_ESP_Gy3z_D2y = I_ESP_H2y3z_Py+ABY*I_ESP_Gy3z_Py;
    Double I_ESP_G4z_D2y = I_ESP_Hy4z_Py+ABY*I_ESP_G4z_Py;
    Double I_ESP_G4x_D2z = I_ESP_H4xz_Pz+ABZ*I_ESP_G4x_Pz;
    Double I_ESP_G3xy_D2z = I_ESP_H3xyz_Pz+ABZ*I_ESP_G3xy_Pz;
    Double I_ESP_G3xz_D2z = I_ESP_H3x2z_Pz+ABZ*I_ESP_G3xz_Pz;
    Double I_ESP_G2x2y_D2z = I_ESP_H2x2yz_Pz+ABZ*I_ESP_G2x2y_Pz;
    Double I_ESP_G2xyz_D2z = I_ESP_H2xy2z_Pz+ABZ*I_ESP_G2xyz_Pz;
    Double I_ESP_G2x2z_D2z = I_ESP_H2x3z_Pz+ABZ*I_ESP_G2x2z_Pz;
    Double I_ESP_Gx3y_D2z = I_ESP_Hx3yz_Pz+ABZ*I_ESP_Gx3y_Pz;
    Double I_ESP_Gx2yz_D2z = I_ESP_Hx2y2z_Pz+ABZ*I_ESP_Gx2yz_Pz;
    Double I_ESP_Gxy2z_D2z = I_ESP_Hxy3z_Pz+ABZ*I_ESP_Gxy2z_Pz;
    Double I_ESP_Gx3z_D2z = I_ESP_Hx4z_Pz+ABZ*I_ESP_Gx3z_Pz;
    Double I_ESP_G4y_D2z = I_ESP_H4yz_Pz+ABZ*I_ESP_G4y_Pz;
    Double I_ESP_G3yz_D2z = I_ESP_H3y2z_Pz+ABZ*I_ESP_G3yz_Pz;
    Double I_ESP_G2y2z_D2z = I_ESP_H2y3z_Pz+ABZ*I_ESP_G2y2z_Pz;
    Double I_ESP_Gy3z_D2z = I_ESP_Hy4z_Pz+ABZ*I_ESP_Gy3z_Pz;
    Double I_ESP_G4z_D2z = I_ESP_H5z_Pz+ABZ*I_ESP_G4z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_F_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 20 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_D
     * RHS shell quartet name: SQ_ESP_F_D
     ************************************************************/
    Double I_ESP_F3x_F3x = I_ESP_G4x_D2x+ABX*I_ESP_F3x_D2x;
    Double I_ESP_F2xy_F3x = I_ESP_G3xy_D2x+ABX*I_ESP_F2xy_D2x;
    Double I_ESP_F2xz_F3x = I_ESP_G3xz_D2x+ABX*I_ESP_F2xz_D2x;
    Double I_ESP_Fx2y_F3x = I_ESP_G2x2y_D2x+ABX*I_ESP_Fx2y_D2x;
    Double I_ESP_Fxyz_F3x = I_ESP_G2xyz_D2x+ABX*I_ESP_Fxyz_D2x;
    Double I_ESP_Fx2z_F3x = I_ESP_G2x2z_D2x+ABX*I_ESP_Fx2z_D2x;
    Double I_ESP_F3y_F3x = I_ESP_Gx3y_D2x+ABX*I_ESP_F3y_D2x;
    Double I_ESP_F2yz_F3x = I_ESP_Gx2yz_D2x+ABX*I_ESP_F2yz_D2x;
    Double I_ESP_Fy2z_F3x = I_ESP_Gxy2z_D2x+ABX*I_ESP_Fy2z_D2x;
    Double I_ESP_F3z_F3x = I_ESP_Gx3z_D2x+ABX*I_ESP_F3z_D2x;
    Double I_ESP_F3x_F2xy = I_ESP_G3xy_D2x+ABY*I_ESP_F3x_D2x;
    Double I_ESP_F2xy_F2xy = I_ESP_G2x2y_D2x+ABY*I_ESP_F2xy_D2x;
    Double I_ESP_F2xz_F2xy = I_ESP_G2xyz_D2x+ABY*I_ESP_F2xz_D2x;
    Double I_ESP_Fx2y_F2xy = I_ESP_Gx3y_D2x+ABY*I_ESP_Fx2y_D2x;
    Double I_ESP_Fxyz_F2xy = I_ESP_Gx2yz_D2x+ABY*I_ESP_Fxyz_D2x;
    Double I_ESP_Fx2z_F2xy = I_ESP_Gxy2z_D2x+ABY*I_ESP_Fx2z_D2x;
    Double I_ESP_F3y_F2xy = I_ESP_G4y_D2x+ABY*I_ESP_F3y_D2x;
    Double I_ESP_F2yz_F2xy = I_ESP_G3yz_D2x+ABY*I_ESP_F2yz_D2x;
    Double I_ESP_Fy2z_F2xy = I_ESP_G2y2z_D2x+ABY*I_ESP_Fy2z_D2x;
    Double I_ESP_F3z_F2xy = I_ESP_Gy3z_D2x+ABY*I_ESP_F3z_D2x;
    Double I_ESP_F3x_F2xz = I_ESP_G3xz_D2x+ABZ*I_ESP_F3x_D2x;
    Double I_ESP_F2xy_F2xz = I_ESP_G2xyz_D2x+ABZ*I_ESP_F2xy_D2x;
    Double I_ESP_F2xz_F2xz = I_ESP_G2x2z_D2x+ABZ*I_ESP_F2xz_D2x;
    Double I_ESP_Fx2y_F2xz = I_ESP_Gx2yz_D2x+ABZ*I_ESP_Fx2y_D2x;
    Double I_ESP_Fxyz_F2xz = I_ESP_Gxy2z_D2x+ABZ*I_ESP_Fxyz_D2x;
    Double I_ESP_Fx2z_F2xz = I_ESP_Gx3z_D2x+ABZ*I_ESP_Fx2z_D2x;
    Double I_ESP_F3y_F2xz = I_ESP_G3yz_D2x+ABZ*I_ESP_F3y_D2x;
    Double I_ESP_F2yz_F2xz = I_ESP_G2y2z_D2x+ABZ*I_ESP_F2yz_D2x;
    Double I_ESP_Fy2z_F2xz = I_ESP_Gy3z_D2x+ABZ*I_ESP_Fy2z_D2x;
    Double I_ESP_F3z_F2xz = I_ESP_G4z_D2x+ABZ*I_ESP_F3z_D2x;
    Double I_ESP_F3x_Fx2y = I_ESP_G4x_D2y+ABX*I_ESP_F3x_D2y;
    Double I_ESP_F2xy_Fx2y = I_ESP_G3xy_D2y+ABX*I_ESP_F2xy_D2y;
    Double I_ESP_F2xz_Fx2y = I_ESP_G3xz_D2y+ABX*I_ESP_F2xz_D2y;
    Double I_ESP_Fx2y_Fx2y = I_ESP_G2x2y_D2y+ABX*I_ESP_Fx2y_D2y;
    Double I_ESP_Fxyz_Fx2y = I_ESP_G2xyz_D2y+ABX*I_ESP_Fxyz_D2y;
    Double I_ESP_Fx2z_Fx2y = I_ESP_G2x2z_D2y+ABX*I_ESP_Fx2z_D2y;
    Double I_ESP_F3y_Fx2y = I_ESP_Gx3y_D2y+ABX*I_ESP_F3y_D2y;
    Double I_ESP_F2yz_Fx2y = I_ESP_Gx2yz_D2y+ABX*I_ESP_F2yz_D2y;
    Double I_ESP_Fy2z_Fx2y = I_ESP_Gxy2z_D2y+ABX*I_ESP_Fy2z_D2y;
    Double I_ESP_F3z_Fx2y = I_ESP_Gx3z_D2y+ABX*I_ESP_F3z_D2y;
    Double I_ESP_F3x_Fx2z = I_ESP_G4x_D2z+ABX*I_ESP_F3x_D2z;
    Double I_ESP_F2xy_Fx2z = I_ESP_G3xy_D2z+ABX*I_ESP_F2xy_D2z;
    Double I_ESP_F2xz_Fx2z = I_ESP_G3xz_D2z+ABX*I_ESP_F2xz_D2z;
    Double I_ESP_Fx2y_Fx2z = I_ESP_G2x2y_D2z+ABX*I_ESP_Fx2y_D2z;
    Double I_ESP_Fxyz_Fx2z = I_ESP_G2xyz_D2z+ABX*I_ESP_Fxyz_D2z;
    Double I_ESP_Fx2z_Fx2z = I_ESP_G2x2z_D2z+ABX*I_ESP_Fx2z_D2z;
    Double I_ESP_F3y_Fx2z = I_ESP_Gx3y_D2z+ABX*I_ESP_F3y_D2z;
    Double I_ESP_F2yz_Fx2z = I_ESP_Gx2yz_D2z+ABX*I_ESP_F2yz_D2z;
    Double I_ESP_Fy2z_Fx2z = I_ESP_Gxy2z_D2z+ABX*I_ESP_Fy2z_D2z;
    Double I_ESP_F3z_Fx2z = I_ESP_Gx3z_D2z+ABX*I_ESP_F3z_D2z;
    Double I_ESP_F3x_F3y = I_ESP_G3xy_D2y+ABY*I_ESP_F3x_D2y;
    Double I_ESP_F2xy_F3y = I_ESP_G2x2y_D2y+ABY*I_ESP_F2xy_D2y;
    Double I_ESP_F2xz_F3y = I_ESP_G2xyz_D2y+ABY*I_ESP_F2xz_D2y;
    Double I_ESP_Fx2y_F3y = I_ESP_Gx3y_D2y+ABY*I_ESP_Fx2y_D2y;
    Double I_ESP_Fxyz_F3y = I_ESP_Gx2yz_D2y+ABY*I_ESP_Fxyz_D2y;
    Double I_ESP_Fx2z_F3y = I_ESP_Gxy2z_D2y+ABY*I_ESP_Fx2z_D2y;
    Double I_ESP_F3y_F3y = I_ESP_G4y_D2y+ABY*I_ESP_F3y_D2y;
    Double I_ESP_F2yz_F3y = I_ESP_G3yz_D2y+ABY*I_ESP_F2yz_D2y;
    Double I_ESP_Fy2z_F3y = I_ESP_G2y2z_D2y+ABY*I_ESP_Fy2z_D2y;
    Double I_ESP_F3z_F3y = I_ESP_Gy3z_D2y+ABY*I_ESP_F3z_D2y;
    Double I_ESP_F3x_F2yz = I_ESP_G3xz_D2y+ABZ*I_ESP_F3x_D2y;
    Double I_ESP_F2xy_F2yz = I_ESP_G2xyz_D2y+ABZ*I_ESP_F2xy_D2y;
    Double I_ESP_F2xz_F2yz = I_ESP_G2x2z_D2y+ABZ*I_ESP_F2xz_D2y;
    Double I_ESP_Fx2y_F2yz = I_ESP_Gx2yz_D2y+ABZ*I_ESP_Fx2y_D2y;
    Double I_ESP_Fxyz_F2yz = I_ESP_Gxy2z_D2y+ABZ*I_ESP_Fxyz_D2y;
    Double I_ESP_Fx2z_F2yz = I_ESP_Gx3z_D2y+ABZ*I_ESP_Fx2z_D2y;
    Double I_ESP_F3y_F2yz = I_ESP_G3yz_D2y+ABZ*I_ESP_F3y_D2y;
    Double I_ESP_F2yz_F2yz = I_ESP_G2y2z_D2y+ABZ*I_ESP_F2yz_D2y;
    Double I_ESP_Fy2z_F2yz = I_ESP_Gy3z_D2y+ABZ*I_ESP_Fy2z_D2y;
    Double I_ESP_F3z_F2yz = I_ESP_G4z_D2y+ABZ*I_ESP_F3z_D2y;
    Double I_ESP_F3x_F3z = I_ESP_G3xz_D2z+ABZ*I_ESP_F3x_D2z;
    Double I_ESP_F2xy_F3z = I_ESP_G2xyz_D2z+ABZ*I_ESP_F2xy_D2z;
    Double I_ESP_F2xz_F3z = I_ESP_G2x2z_D2z+ABZ*I_ESP_F2xz_D2z;
    Double I_ESP_Fx2y_F3z = I_ESP_Gx2yz_D2z+ABZ*I_ESP_Fx2y_D2z;
    Double I_ESP_Fxyz_F3z = I_ESP_Gxy2z_D2z+ABZ*I_ESP_Fxyz_D2z;
    Double I_ESP_Fx2z_F3z = I_ESP_Gx3z_D2z+ABZ*I_ESP_Fx2z_D2z;
    Double I_ESP_F3y_F3z = I_ESP_G3yz_D2z+ABZ*I_ESP_F3y_D2z;
    Double I_ESP_F2yz_F3z = I_ESP_G2y2z_D2z+ABZ*I_ESP_F2yz_D2z;
    Double I_ESP_Fy2z_F3z = I_ESP_Gy3z_D2z+ABZ*I_ESP_Fy2z_D2z;
    Double I_ESP_F3z_F3z = I_ESP_G4z_D2z+ABZ*I_ESP_F3z_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_I_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 24 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_S
     * RHS shell quartet name: SQ_ESP_I_S
     ************************************************************/
    Double I_ESP_I6x_Px = I_ESP_K7x_S+ABX*I_ESP_I6x_S;
    Double I_ESP_I5xy_Px = I_ESP_K6xy_S+ABX*I_ESP_I5xy_S;
    Double I_ESP_I5xz_Px = I_ESP_K6xz_S+ABX*I_ESP_I5xz_S;
    Double I_ESP_I4x2y_Px = I_ESP_K5x2y_S+ABX*I_ESP_I4x2y_S;
    Double I_ESP_I4xyz_Px = I_ESP_K5xyz_S+ABX*I_ESP_I4xyz_S;
    Double I_ESP_I4x2z_Px = I_ESP_K5x2z_S+ABX*I_ESP_I4x2z_S;
    Double I_ESP_I3x3y_Px = I_ESP_K4x3y_S+ABX*I_ESP_I3x3y_S;
    Double I_ESP_I3x2yz_Px = I_ESP_K4x2yz_S+ABX*I_ESP_I3x2yz_S;
    Double I_ESP_I3xy2z_Px = I_ESP_K4xy2z_S+ABX*I_ESP_I3xy2z_S;
    Double I_ESP_I3x3z_Px = I_ESP_K4x3z_S+ABX*I_ESP_I3x3z_S;
    Double I_ESP_I2x4y_Px = I_ESP_K3x4y_S+ABX*I_ESP_I2x4y_S;
    Double I_ESP_I2x3yz_Px = I_ESP_K3x3yz_S+ABX*I_ESP_I2x3yz_S;
    Double I_ESP_I2x2y2z_Px = I_ESP_K3x2y2z_S+ABX*I_ESP_I2x2y2z_S;
    Double I_ESP_I2xy3z_Px = I_ESP_K3xy3z_S+ABX*I_ESP_I2xy3z_S;
    Double I_ESP_I2x4z_Px = I_ESP_K3x4z_S+ABX*I_ESP_I2x4z_S;
    Double I_ESP_Ix5y_Px = I_ESP_K2x5y_S+ABX*I_ESP_Ix5y_S;
    Double I_ESP_Ix4yz_Px = I_ESP_K2x4yz_S+ABX*I_ESP_Ix4yz_S;
    Double I_ESP_Ix3y2z_Px = I_ESP_K2x3y2z_S+ABX*I_ESP_Ix3y2z_S;
    Double I_ESP_Ix2y3z_Px = I_ESP_K2x2y3z_S+ABX*I_ESP_Ix2y3z_S;
    Double I_ESP_Ixy4z_Px = I_ESP_K2xy4z_S+ABX*I_ESP_Ixy4z_S;
    Double I_ESP_Ix5z_Px = I_ESP_K2x5z_S+ABX*I_ESP_Ix5z_S;
    Double I_ESP_I4x2y_Py = I_ESP_K4x3y_S+ABY*I_ESP_I4x2y_S;
    Double I_ESP_I4xyz_Py = I_ESP_K4x2yz_S+ABY*I_ESP_I4xyz_S;
    Double I_ESP_I3x3y_Py = I_ESP_K3x4y_S+ABY*I_ESP_I3x3y_S;
    Double I_ESP_I3x2yz_Py = I_ESP_K3x3yz_S+ABY*I_ESP_I3x2yz_S;
    Double I_ESP_I3xy2z_Py = I_ESP_K3x2y2z_S+ABY*I_ESP_I3xy2z_S;
    Double I_ESP_I2x4y_Py = I_ESP_K2x5y_S+ABY*I_ESP_I2x4y_S;
    Double I_ESP_I2x3yz_Py = I_ESP_K2x4yz_S+ABY*I_ESP_I2x3yz_S;
    Double I_ESP_I2x2y2z_Py = I_ESP_K2x3y2z_S+ABY*I_ESP_I2x2y2z_S;
    Double I_ESP_I2xy3z_Py = I_ESP_K2x2y3z_S+ABY*I_ESP_I2xy3z_S;
    Double I_ESP_Ix5y_Py = I_ESP_Kx6y_S+ABY*I_ESP_Ix5y_S;
    Double I_ESP_Ix4yz_Py = I_ESP_Kx5yz_S+ABY*I_ESP_Ix4yz_S;
    Double I_ESP_Ix3y2z_Py = I_ESP_Kx4y2z_S+ABY*I_ESP_Ix3y2z_S;
    Double I_ESP_Ix2y3z_Py = I_ESP_Kx3y3z_S+ABY*I_ESP_Ix2y3z_S;
    Double I_ESP_Ixy4z_Py = I_ESP_Kx2y4z_S+ABY*I_ESP_Ixy4z_S;
    Double I_ESP_I6y_Py = I_ESP_K7y_S+ABY*I_ESP_I6y_S;
    Double I_ESP_I5yz_Py = I_ESP_K6yz_S+ABY*I_ESP_I5yz_S;
    Double I_ESP_I4y2z_Py = I_ESP_K5y2z_S+ABY*I_ESP_I4y2z_S;
    Double I_ESP_I3y3z_Py = I_ESP_K4y3z_S+ABY*I_ESP_I3y3z_S;
    Double I_ESP_I2y4z_Py = I_ESP_K3y4z_S+ABY*I_ESP_I2y4z_S;
    Double I_ESP_Iy5z_Py = I_ESP_K2y5z_S+ABY*I_ESP_Iy5z_S;
    Double I_ESP_I4xyz_Pz = I_ESP_K4xy2z_S+ABZ*I_ESP_I4xyz_S;
    Double I_ESP_I4x2z_Pz = I_ESP_K4x3z_S+ABZ*I_ESP_I4x2z_S;
    Double I_ESP_I3x2yz_Pz = I_ESP_K3x2y2z_S+ABZ*I_ESP_I3x2yz_S;
    Double I_ESP_I3xy2z_Pz = I_ESP_K3xy3z_S+ABZ*I_ESP_I3xy2z_S;
    Double I_ESP_I3x3z_Pz = I_ESP_K3x4z_S+ABZ*I_ESP_I3x3z_S;
    Double I_ESP_I2x3yz_Pz = I_ESP_K2x3y2z_S+ABZ*I_ESP_I2x3yz_S;
    Double I_ESP_I2x2y2z_Pz = I_ESP_K2x2y3z_S+ABZ*I_ESP_I2x2y2z_S;
    Double I_ESP_I2xy3z_Pz = I_ESP_K2xy4z_S+ABZ*I_ESP_I2xy3z_S;
    Double I_ESP_I2x4z_Pz = I_ESP_K2x5z_S+ABZ*I_ESP_I2x4z_S;
    Double I_ESP_Ix4yz_Pz = I_ESP_Kx4y2z_S+ABZ*I_ESP_Ix4yz_S;
    Double I_ESP_Ix3y2z_Pz = I_ESP_Kx3y3z_S+ABZ*I_ESP_Ix3y2z_S;
    Double I_ESP_Ix2y3z_Pz = I_ESP_Kx2y4z_S+ABZ*I_ESP_Ix2y3z_S;
    Double I_ESP_Ixy4z_Pz = I_ESP_Kxy5z_S+ABZ*I_ESP_Ixy4z_S;
    Double I_ESP_Ix5z_Pz = I_ESP_Kx6z_S+ABZ*I_ESP_Ix5z_S;
    Double I_ESP_I4y2z_Pz = I_ESP_K4y3z_S+ABZ*I_ESP_I4y2z_S;
    Double I_ESP_I3y3z_Pz = I_ESP_K3y4z_S+ABZ*I_ESP_I3y3z_S;
    Double I_ESP_I2y4z_Pz = I_ESP_K2y5z_S+ABZ*I_ESP_I2y4z_S;
    Double I_ESP_Iy5z_Pz = I_ESP_Ky6z_S+ABZ*I_ESP_Iy5z_S;
    Double I_ESP_I6z_Pz = I_ESP_K7z_S+ABZ*I_ESP_I6z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 66 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_P
     * RHS shell quartet name: SQ_ESP_H_P
     ************************************************************/
    Double I_ESP_H5x_D2x = I_ESP_I6x_Px+ABX*I_ESP_H5x_Px;
    Double I_ESP_H4xy_D2x = I_ESP_I5xy_Px+ABX*I_ESP_H4xy_Px;
    Double I_ESP_H4xz_D2x = I_ESP_I5xz_Px+ABX*I_ESP_H4xz_Px;
    Double I_ESP_H3x2y_D2x = I_ESP_I4x2y_Px+ABX*I_ESP_H3x2y_Px;
    Double I_ESP_H3xyz_D2x = I_ESP_I4xyz_Px+ABX*I_ESP_H3xyz_Px;
    Double I_ESP_H3x2z_D2x = I_ESP_I4x2z_Px+ABX*I_ESP_H3x2z_Px;
    Double I_ESP_H2x3y_D2x = I_ESP_I3x3y_Px+ABX*I_ESP_H2x3y_Px;
    Double I_ESP_H2x2yz_D2x = I_ESP_I3x2yz_Px+ABX*I_ESP_H2x2yz_Px;
    Double I_ESP_H2xy2z_D2x = I_ESP_I3xy2z_Px+ABX*I_ESP_H2xy2z_Px;
    Double I_ESP_H2x3z_D2x = I_ESP_I3x3z_Px+ABX*I_ESP_H2x3z_Px;
    Double I_ESP_Hx4y_D2x = I_ESP_I2x4y_Px+ABX*I_ESP_Hx4y_Px;
    Double I_ESP_Hx3yz_D2x = I_ESP_I2x3yz_Px+ABX*I_ESP_Hx3yz_Px;
    Double I_ESP_Hx2y2z_D2x = I_ESP_I2x2y2z_Px+ABX*I_ESP_Hx2y2z_Px;
    Double I_ESP_Hxy3z_D2x = I_ESP_I2xy3z_Px+ABX*I_ESP_Hxy3z_Px;
    Double I_ESP_Hx4z_D2x = I_ESP_I2x4z_Px+ABX*I_ESP_Hx4z_Px;
    Double I_ESP_H5y_D2x = I_ESP_Ix5y_Px+ABX*I_ESP_H5y_Px;
    Double I_ESP_H4yz_D2x = I_ESP_Ix4yz_Px+ABX*I_ESP_H4yz_Px;
    Double I_ESP_H3y2z_D2x = I_ESP_Ix3y2z_Px+ABX*I_ESP_H3y2z_Px;
    Double I_ESP_H2y3z_D2x = I_ESP_Ix2y3z_Px+ABX*I_ESP_H2y3z_Px;
    Double I_ESP_Hy4z_D2x = I_ESP_Ixy4z_Px+ABX*I_ESP_Hy4z_Px;
    Double I_ESP_H5z_D2x = I_ESP_Ix5z_Px+ABX*I_ESP_H5z_Px;
    Double I_ESP_H4xy_D2y = I_ESP_I4x2y_Py+ABY*I_ESP_H4xy_Py;
    Double I_ESP_H4xz_D2y = I_ESP_I4xyz_Py+ABY*I_ESP_H4xz_Py;
    Double I_ESP_H3x2y_D2y = I_ESP_I3x3y_Py+ABY*I_ESP_H3x2y_Py;
    Double I_ESP_H3xyz_D2y = I_ESP_I3x2yz_Py+ABY*I_ESP_H3xyz_Py;
    Double I_ESP_H3x2z_D2y = I_ESP_I3xy2z_Py+ABY*I_ESP_H3x2z_Py;
    Double I_ESP_H2x3y_D2y = I_ESP_I2x4y_Py+ABY*I_ESP_H2x3y_Py;
    Double I_ESP_H2x2yz_D2y = I_ESP_I2x3yz_Py+ABY*I_ESP_H2x2yz_Py;
    Double I_ESP_H2xy2z_D2y = I_ESP_I2x2y2z_Py+ABY*I_ESP_H2xy2z_Py;
    Double I_ESP_H2x3z_D2y = I_ESP_I2xy3z_Py+ABY*I_ESP_H2x3z_Py;
    Double I_ESP_Hx4y_D2y = I_ESP_Ix5y_Py+ABY*I_ESP_Hx4y_Py;
    Double I_ESP_Hx3yz_D2y = I_ESP_Ix4yz_Py+ABY*I_ESP_Hx3yz_Py;
    Double I_ESP_Hx2y2z_D2y = I_ESP_Ix3y2z_Py+ABY*I_ESP_Hx2y2z_Py;
    Double I_ESP_Hxy3z_D2y = I_ESP_Ix2y3z_Py+ABY*I_ESP_Hxy3z_Py;
    Double I_ESP_Hx4z_D2y = I_ESP_Ixy4z_Py+ABY*I_ESP_Hx4z_Py;
    Double I_ESP_H5y_D2y = I_ESP_I6y_Py+ABY*I_ESP_H5y_Py;
    Double I_ESP_H4yz_D2y = I_ESP_I5yz_Py+ABY*I_ESP_H4yz_Py;
    Double I_ESP_H3y2z_D2y = I_ESP_I4y2z_Py+ABY*I_ESP_H3y2z_Py;
    Double I_ESP_H2y3z_D2y = I_ESP_I3y3z_Py+ABY*I_ESP_H2y3z_Py;
    Double I_ESP_Hy4z_D2y = I_ESP_I2y4z_Py+ABY*I_ESP_Hy4z_Py;
    Double I_ESP_H5z_D2y = I_ESP_Iy5z_Py+ABY*I_ESP_H5z_Py;
    Double I_ESP_H4xy_D2z = I_ESP_I4xyz_Pz+ABZ*I_ESP_H4xy_Pz;
    Double I_ESP_H4xz_D2z = I_ESP_I4x2z_Pz+ABZ*I_ESP_H4xz_Pz;
    Double I_ESP_H3x2y_D2z = I_ESP_I3x2yz_Pz+ABZ*I_ESP_H3x2y_Pz;
    Double I_ESP_H3xyz_D2z = I_ESP_I3xy2z_Pz+ABZ*I_ESP_H3xyz_Pz;
    Double I_ESP_H3x2z_D2z = I_ESP_I3x3z_Pz+ABZ*I_ESP_H3x2z_Pz;
    Double I_ESP_H2x3y_D2z = I_ESP_I2x3yz_Pz+ABZ*I_ESP_H2x3y_Pz;
    Double I_ESP_H2x2yz_D2z = I_ESP_I2x2y2z_Pz+ABZ*I_ESP_H2x2yz_Pz;
    Double I_ESP_H2xy2z_D2z = I_ESP_I2xy3z_Pz+ABZ*I_ESP_H2xy2z_Pz;
    Double I_ESP_H2x3z_D2z = I_ESP_I2x4z_Pz+ABZ*I_ESP_H2x3z_Pz;
    Double I_ESP_Hx4y_D2z = I_ESP_Ix4yz_Pz+ABZ*I_ESP_Hx4y_Pz;
    Double I_ESP_Hx3yz_D2z = I_ESP_Ix3y2z_Pz+ABZ*I_ESP_Hx3yz_Pz;
    Double I_ESP_Hx2y2z_D2z = I_ESP_Ix2y3z_Pz+ABZ*I_ESP_Hx2y2z_Pz;
    Double I_ESP_Hxy3z_D2z = I_ESP_Ixy4z_Pz+ABZ*I_ESP_Hxy3z_Pz;
    Double I_ESP_Hx4z_D2z = I_ESP_Ix5z_Pz+ABZ*I_ESP_Hx4z_Pz;
    Double I_ESP_H4yz_D2z = I_ESP_I4y2z_Pz+ABZ*I_ESP_H4yz_Pz;
    Double I_ESP_H3y2z_D2z = I_ESP_I3y3z_Pz+ABZ*I_ESP_H3y2z_Pz;
    Double I_ESP_H2y3z_D2z = I_ESP_I2y4z_Pz+ABZ*I_ESP_H2y3z_Pz;
    Double I_ESP_Hy4z_D2z = I_ESP_Iy5z_Pz+ABZ*I_ESP_Hy4z_Pz;
    Double I_ESP_H5z_D2z = I_ESP_I6z_Pz+ABZ*I_ESP_H5z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 51 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_D
     * RHS shell quartet name: SQ_ESP_G_D
     ************************************************************/
    Double I_ESP_G4x_F3x = I_ESP_H5x_D2x+ABX*I_ESP_G4x_D2x;
    Double I_ESP_G3xy_F3x = I_ESP_H4xy_D2x+ABX*I_ESP_G3xy_D2x;
    Double I_ESP_G3xz_F3x = I_ESP_H4xz_D2x+ABX*I_ESP_G3xz_D2x;
    Double I_ESP_G2x2y_F3x = I_ESP_H3x2y_D2x+ABX*I_ESP_G2x2y_D2x;
    Double I_ESP_G2xyz_F3x = I_ESP_H3xyz_D2x+ABX*I_ESP_G2xyz_D2x;
    Double I_ESP_G2x2z_F3x = I_ESP_H3x2z_D2x+ABX*I_ESP_G2x2z_D2x;
    Double I_ESP_Gx3y_F3x = I_ESP_H2x3y_D2x+ABX*I_ESP_Gx3y_D2x;
    Double I_ESP_Gx2yz_F3x = I_ESP_H2x2yz_D2x+ABX*I_ESP_Gx2yz_D2x;
    Double I_ESP_Gxy2z_F3x = I_ESP_H2xy2z_D2x+ABX*I_ESP_Gxy2z_D2x;
    Double I_ESP_Gx3z_F3x = I_ESP_H2x3z_D2x+ABX*I_ESP_Gx3z_D2x;
    Double I_ESP_G4y_F3x = I_ESP_Hx4y_D2x+ABX*I_ESP_G4y_D2x;
    Double I_ESP_G3yz_F3x = I_ESP_Hx3yz_D2x+ABX*I_ESP_G3yz_D2x;
    Double I_ESP_G2y2z_F3x = I_ESP_Hx2y2z_D2x+ABX*I_ESP_G2y2z_D2x;
    Double I_ESP_Gy3z_F3x = I_ESP_Hxy3z_D2x+ABX*I_ESP_Gy3z_D2x;
    Double I_ESP_G4z_F3x = I_ESP_Hx4z_D2x+ABX*I_ESP_G4z_D2x;
    Double I_ESP_G3xy_F2xy = I_ESP_H3x2y_D2x+ABY*I_ESP_G3xy_D2x;
    Double I_ESP_G3xz_F2xy = I_ESP_H3xyz_D2x+ABY*I_ESP_G3xz_D2x;
    Double I_ESP_G2x2y_F2xy = I_ESP_H2x3y_D2x+ABY*I_ESP_G2x2y_D2x;
    Double I_ESP_G2xyz_F2xy = I_ESP_H2x2yz_D2x+ABY*I_ESP_G2xyz_D2x;
    Double I_ESP_G2x2z_F2xy = I_ESP_H2xy2z_D2x+ABY*I_ESP_G2x2z_D2x;
    Double I_ESP_Gx3y_F2xy = I_ESP_Hx4y_D2x+ABY*I_ESP_Gx3y_D2x;
    Double I_ESP_Gx2yz_F2xy = I_ESP_Hx3yz_D2x+ABY*I_ESP_Gx2yz_D2x;
    Double I_ESP_Gxy2z_F2xy = I_ESP_Hx2y2z_D2x+ABY*I_ESP_Gxy2z_D2x;
    Double I_ESP_Gx3z_F2xy = I_ESP_Hxy3z_D2x+ABY*I_ESP_Gx3z_D2x;
    Double I_ESP_G4y_F2xy = I_ESP_H5y_D2x+ABY*I_ESP_G4y_D2x;
    Double I_ESP_G3yz_F2xy = I_ESP_H4yz_D2x+ABY*I_ESP_G3yz_D2x;
    Double I_ESP_G2y2z_F2xy = I_ESP_H3y2z_D2x+ABY*I_ESP_G2y2z_D2x;
    Double I_ESP_Gy3z_F2xy = I_ESP_H2y3z_D2x+ABY*I_ESP_Gy3z_D2x;
    Double I_ESP_G4z_F2xy = I_ESP_Hy4z_D2x+ABY*I_ESP_G4z_D2x;
    Double I_ESP_G3xz_F2xz = I_ESP_H3x2z_D2x+ABZ*I_ESP_G3xz_D2x;
    Double I_ESP_G2xyz_F2xz = I_ESP_H2xy2z_D2x+ABZ*I_ESP_G2xyz_D2x;
    Double I_ESP_G2x2z_F2xz = I_ESP_H2x3z_D2x+ABZ*I_ESP_G2x2z_D2x;
    Double I_ESP_Gx2yz_F2xz = I_ESP_Hx2y2z_D2x+ABZ*I_ESP_Gx2yz_D2x;
    Double I_ESP_Gxy2z_F2xz = I_ESP_Hxy3z_D2x+ABZ*I_ESP_Gxy2z_D2x;
    Double I_ESP_Gx3z_F2xz = I_ESP_Hx4z_D2x+ABZ*I_ESP_Gx3z_D2x;
    Double I_ESP_G3yz_F2xz = I_ESP_H3y2z_D2x+ABZ*I_ESP_G3yz_D2x;
    Double I_ESP_G2y2z_F2xz = I_ESP_H2y3z_D2x+ABZ*I_ESP_G2y2z_D2x;
    Double I_ESP_Gy3z_F2xz = I_ESP_Hy4z_D2x+ABZ*I_ESP_Gy3z_D2x;
    Double I_ESP_G4z_F2xz = I_ESP_H5z_D2x+ABZ*I_ESP_G4z_D2x;
    Double I_ESP_G3xz_Fx2y = I_ESP_H4xz_D2y+ABX*I_ESP_G3xz_D2y;
    Double I_ESP_G2xyz_Fx2y = I_ESP_H3xyz_D2y+ABX*I_ESP_G2xyz_D2y;
    Double I_ESP_G2x2z_Fx2y = I_ESP_H3x2z_D2y+ABX*I_ESP_G2x2z_D2y;
    Double I_ESP_Gx2yz_Fx2y = I_ESP_H2x2yz_D2y+ABX*I_ESP_Gx2yz_D2y;
    Double I_ESP_Gxy2z_Fx2y = I_ESP_H2xy2z_D2y+ABX*I_ESP_Gxy2z_D2y;
    Double I_ESP_Gx3z_Fx2y = I_ESP_H2x3z_D2y+ABX*I_ESP_Gx3z_D2y;
    Double I_ESP_G3yz_Fx2y = I_ESP_Hx3yz_D2y+ABX*I_ESP_G3yz_D2y;
    Double I_ESP_G2y2z_Fx2y = I_ESP_Hx2y2z_D2y+ABX*I_ESP_G2y2z_D2y;
    Double I_ESP_Gy3z_Fx2y = I_ESP_Hxy3z_D2y+ABX*I_ESP_Gy3z_D2y;
    Double I_ESP_G4z_Fx2y = I_ESP_Hx4z_D2y+ABX*I_ESP_G4z_D2y;
    Double I_ESP_G3xy_Fx2z = I_ESP_H4xy_D2z+ABX*I_ESP_G3xy_D2z;
    Double I_ESP_G2x2y_Fx2z = I_ESP_H3x2y_D2z+ABX*I_ESP_G2x2y_D2z;
    Double I_ESP_G2xyz_Fx2z = I_ESP_H3xyz_D2z+ABX*I_ESP_G2xyz_D2z;
    Double I_ESP_Gx3y_Fx2z = I_ESP_H2x3y_D2z+ABX*I_ESP_Gx3y_D2z;
    Double I_ESP_Gx2yz_Fx2z = I_ESP_H2x2yz_D2z+ABX*I_ESP_Gx2yz_D2z;
    Double I_ESP_Gxy2z_Fx2z = I_ESP_H2xy2z_D2z+ABX*I_ESP_Gxy2z_D2z;
    Double I_ESP_G4y_Fx2z = I_ESP_Hx4y_D2z+ABX*I_ESP_G4y_D2z;
    Double I_ESP_G3yz_Fx2z = I_ESP_Hx3yz_D2z+ABX*I_ESP_G3yz_D2z;
    Double I_ESP_G2y2z_Fx2z = I_ESP_Hx2y2z_D2z+ABX*I_ESP_G2y2z_D2z;
    Double I_ESP_Gy3z_Fx2z = I_ESP_Hxy3z_D2z+ABX*I_ESP_Gy3z_D2z;
    Double I_ESP_G4x_F3y = I_ESP_H4xy_D2y+ABY*I_ESP_G4x_D2y;
    Double I_ESP_G3xy_F3y = I_ESP_H3x2y_D2y+ABY*I_ESP_G3xy_D2y;
    Double I_ESP_G3xz_F3y = I_ESP_H3xyz_D2y+ABY*I_ESP_G3xz_D2y;
    Double I_ESP_G2x2y_F3y = I_ESP_H2x3y_D2y+ABY*I_ESP_G2x2y_D2y;
    Double I_ESP_G2xyz_F3y = I_ESP_H2x2yz_D2y+ABY*I_ESP_G2xyz_D2y;
    Double I_ESP_G2x2z_F3y = I_ESP_H2xy2z_D2y+ABY*I_ESP_G2x2z_D2y;
    Double I_ESP_Gx3y_F3y = I_ESP_Hx4y_D2y+ABY*I_ESP_Gx3y_D2y;
    Double I_ESP_Gx2yz_F3y = I_ESP_Hx3yz_D2y+ABY*I_ESP_Gx2yz_D2y;
    Double I_ESP_Gxy2z_F3y = I_ESP_Hx2y2z_D2y+ABY*I_ESP_Gxy2z_D2y;
    Double I_ESP_Gx3z_F3y = I_ESP_Hxy3z_D2y+ABY*I_ESP_Gx3z_D2y;
    Double I_ESP_G4y_F3y = I_ESP_H5y_D2y+ABY*I_ESP_G4y_D2y;
    Double I_ESP_G3yz_F3y = I_ESP_H4yz_D2y+ABY*I_ESP_G3yz_D2y;
    Double I_ESP_G2y2z_F3y = I_ESP_H3y2z_D2y+ABY*I_ESP_G2y2z_D2y;
    Double I_ESP_Gy3z_F3y = I_ESP_H2y3z_D2y+ABY*I_ESP_Gy3z_D2y;
    Double I_ESP_G4z_F3y = I_ESP_Hy4z_D2y+ABY*I_ESP_G4z_D2y;
    Double I_ESP_G3xz_F2yz = I_ESP_H3x2z_D2y+ABZ*I_ESP_G3xz_D2y;
    Double I_ESP_G2xyz_F2yz = I_ESP_H2xy2z_D2y+ABZ*I_ESP_G2xyz_D2y;
    Double I_ESP_G2x2z_F2yz = I_ESP_H2x3z_D2y+ABZ*I_ESP_G2x2z_D2y;
    Double I_ESP_Gx2yz_F2yz = I_ESP_Hx2y2z_D2y+ABZ*I_ESP_Gx2yz_D2y;
    Double I_ESP_Gxy2z_F2yz = I_ESP_Hxy3z_D2y+ABZ*I_ESP_Gxy2z_D2y;
    Double I_ESP_Gx3z_F2yz = I_ESP_Hx4z_D2y+ABZ*I_ESP_Gx3z_D2y;
    Double I_ESP_G3yz_F2yz = I_ESP_H3y2z_D2y+ABZ*I_ESP_G3yz_D2y;
    Double I_ESP_G2y2z_F2yz = I_ESP_H2y3z_D2y+ABZ*I_ESP_G2y2z_D2y;
    Double I_ESP_Gy3z_F2yz = I_ESP_Hy4z_D2y+ABZ*I_ESP_Gy3z_D2y;
    Double I_ESP_G4z_F2yz = I_ESP_H5z_D2y+ABZ*I_ESP_G4z_D2y;
    Double I_ESP_G4x_F3z = I_ESP_H4xz_D2z+ABZ*I_ESP_G4x_D2z;
    Double I_ESP_G3xy_F3z = I_ESP_H3xyz_D2z+ABZ*I_ESP_G3xy_D2z;
    Double I_ESP_G3xz_F3z = I_ESP_H3x2z_D2z+ABZ*I_ESP_G3xz_D2z;
    Double I_ESP_G2x2y_F3z = I_ESP_H2x2yz_D2z+ABZ*I_ESP_G2x2y_D2z;
    Double I_ESP_G2xyz_F3z = I_ESP_H2xy2z_D2z+ABZ*I_ESP_G2xyz_D2z;
    Double I_ESP_G2x2z_F3z = I_ESP_H2x3z_D2z+ABZ*I_ESP_G2x2z_D2z;
    Double I_ESP_Gx3y_F3z = I_ESP_Hx3yz_D2z+ABZ*I_ESP_Gx3y_D2z;
    Double I_ESP_Gx2yz_F3z = I_ESP_Hx2y2z_D2z+ABZ*I_ESP_Gx2yz_D2z;
    Double I_ESP_Gxy2z_F3z = I_ESP_Hxy3z_D2z+ABZ*I_ESP_Gxy2z_D2z;
    Double I_ESP_Gx3z_F3z = I_ESP_Hx4z_D2z+ABZ*I_ESP_Gx3z_D2z;
    Double I_ESP_G4y_F3z = I_ESP_H4yz_D2z+ABZ*I_ESP_G4y_D2z;
    Double I_ESP_G3yz_F3z = I_ESP_H3y2z_D2z+ABZ*I_ESP_G3yz_D2z;
    Double I_ESP_G2y2z_F3z = I_ESP_H2y3z_D2z+ABZ*I_ESP_G2y2z_D2z;
    Double I_ESP_Gy3z_F3z = I_ESP_Hy4z_D2z+ABZ*I_ESP_Gy3z_D2z;
    Double I_ESP_G4z_F3z = I_ESP_H5z_D2z+ABZ*I_ESP_G4z_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_F_G
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_G_F
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    Double I_ESP_F3x_G4x = I_ESP_G4x_F3x+ABX*I_ESP_F3x_F3x;
    Double I_ESP_F2xy_G4x = I_ESP_G3xy_F3x+ABX*I_ESP_F2xy_F3x;
    Double I_ESP_F2xz_G4x = I_ESP_G3xz_F3x+ABX*I_ESP_F2xz_F3x;
    Double I_ESP_Fx2y_G4x = I_ESP_G2x2y_F3x+ABX*I_ESP_Fx2y_F3x;
    Double I_ESP_Fxyz_G4x = I_ESP_G2xyz_F3x+ABX*I_ESP_Fxyz_F3x;
    Double I_ESP_Fx2z_G4x = I_ESP_G2x2z_F3x+ABX*I_ESP_Fx2z_F3x;
    Double I_ESP_F3y_G4x = I_ESP_Gx3y_F3x+ABX*I_ESP_F3y_F3x;
    Double I_ESP_F2yz_G4x = I_ESP_Gx2yz_F3x+ABX*I_ESP_F2yz_F3x;
    Double I_ESP_Fy2z_G4x = I_ESP_Gxy2z_F3x+ABX*I_ESP_Fy2z_F3x;
    Double I_ESP_F3z_G4x = I_ESP_Gx3z_F3x+ABX*I_ESP_F3z_F3x;
    Double I_ESP_F3x_G3xy = I_ESP_G3xy_F3x+ABY*I_ESP_F3x_F3x;
    Double I_ESP_F2xy_G3xy = I_ESP_G2x2y_F3x+ABY*I_ESP_F2xy_F3x;
    Double I_ESP_F2xz_G3xy = I_ESP_G2xyz_F3x+ABY*I_ESP_F2xz_F3x;
    Double I_ESP_Fx2y_G3xy = I_ESP_Gx3y_F3x+ABY*I_ESP_Fx2y_F3x;
    Double I_ESP_Fxyz_G3xy = I_ESP_Gx2yz_F3x+ABY*I_ESP_Fxyz_F3x;
    Double I_ESP_Fx2z_G3xy = I_ESP_Gxy2z_F3x+ABY*I_ESP_Fx2z_F3x;
    Double I_ESP_F3y_G3xy = I_ESP_G4y_F3x+ABY*I_ESP_F3y_F3x;
    Double I_ESP_F2yz_G3xy = I_ESP_G3yz_F3x+ABY*I_ESP_F2yz_F3x;
    Double I_ESP_Fy2z_G3xy = I_ESP_G2y2z_F3x+ABY*I_ESP_Fy2z_F3x;
    Double I_ESP_F3z_G3xy = I_ESP_Gy3z_F3x+ABY*I_ESP_F3z_F3x;
    Double I_ESP_F3x_G3xz = I_ESP_G3xz_F3x+ABZ*I_ESP_F3x_F3x;
    Double I_ESP_F2xy_G3xz = I_ESP_G2xyz_F3x+ABZ*I_ESP_F2xy_F3x;
    Double I_ESP_F2xz_G3xz = I_ESP_G2x2z_F3x+ABZ*I_ESP_F2xz_F3x;
    Double I_ESP_Fx2y_G3xz = I_ESP_Gx2yz_F3x+ABZ*I_ESP_Fx2y_F3x;
    Double I_ESP_Fxyz_G3xz = I_ESP_Gxy2z_F3x+ABZ*I_ESP_Fxyz_F3x;
    Double I_ESP_Fx2z_G3xz = I_ESP_Gx3z_F3x+ABZ*I_ESP_Fx2z_F3x;
    Double I_ESP_F3y_G3xz = I_ESP_G3yz_F3x+ABZ*I_ESP_F3y_F3x;
    Double I_ESP_F2yz_G3xz = I_ESP_G2y2z_F3x+ABZ*I_ESP_F2yz_F3x;
    Double I_ESP_Fy2z_G3xz = I_ESP_Gy3z_F3x+ABZ*I_ESP_Fy2z_F3x;
    Double I_ESP_F3z_G3xz = I_ESP_G4z_F3x+ABZ*I_ESP_F3z_F3x;
    Double I_ESP_F3x_G2x2y = I_ESP_G3xy_F2xy+ABY*I_ESP_F3x_F2xy;
    Double I_ESP_F2xy_G2x2y = I_ESP_G2x2y_F2xy+ABY*I_ESP_F2xy_F2xy;
    Double I_ESP_F2xz_G2x2y = I_ESP_G2xyz_F2xy+ABY*I_ESP_F2xz_F2xy;
    Double I_ESP_Fx2y_G2x2y = I_ESP_Gx3y_F2xy+ABY*I_ESP_Fx2y_F2xy;
    Double I_ESP_Fxyz_G2x2y = I_ESP_Gx2yz_F2xy+ABY*I_ESP_Fxyz_F2xy;
    Double I_ESP_Fx2z_G2x2y = I_ESP_Gxy2z_F2xy+ABY*I_ESP_Fx2z_F2xy;
    Double I_ESP_F3y_G2x2y = I_ESP_G4y_F2xy+ABY*I_ESP_F3y_F2xy;
    Double I_ESP_F2yz_G2x2y = I_ESP_G3yz_F2xy+ABY*I_ESP_F2yz_F2xy;
    Double I_ESP_Fy2z_G2x2y = I_ESP_G2y2z_F2xy+ABY*I_ESP_Fy2z_F2xy;
    Double I_ESP_F3z_G2x2y = I_ESP_Gy3z_F2xy+ABY*I_ESP_F3z_F2xy;
    Double I_ESP_F3x_G2xyz = I_ESP_G3xz_F2xy+ABZ*I_ESP_F3x_F2xy;
    Double I_ESP_F2xy_G2xyz = I_ESP_G2xyz_F2xy+ABZ*I_ESP_F2xy_F2xy;
    Double I_ESP_F2xz_G2xyz = I_ESP_G2x2z_F2xy+ABZ*I_ESP_F2xz_F2xy;
    Double I_ESP_Fx2y_G2xyz = I_ESP_Gx2yz_F2xy+ABZ*I_ESP_Fx2y_F2xy;
    Double I_ESP_Fxyz_G2xyz = I_ESP_Gxy2z_F2xy+ABZ*I_ESP_Fxyz_F2xy;
    Double I_ESP_Fx2z_G2xyz = I_ESP_Gx3z_F2xy+ABZ*I_ESP_Fx2z_F2xy;
    Double I_ESP_F3y_G2xyz = I_ESP_G3yz_F2xy+ABZ*I_ESP_F3y_F2xy;
    Double I_ESP_F2yz_G2xyz = I_ESP_G2y2z_F2xy+ABZ*I_ESP_F2yz_F2xy;
    Double I_ESP_Fy2z_G2xyz = I_ESP_Gy3z_F2xy+ABZ*I_ESP_Fy2z_F2xy;
    Double I_ESP_F3z_G2xyz = I_ESP_G4z_F2xy+ABZ*I_ESP_F3z_F2xy;
    Double I_ESP_F3x_G2x2z = I_ESP_G3xz_F2xz+ABZ*I_ESP_F3x_F2xz;
    Double I_ESP_F2xy_G2x2z = I_ESP_G2xyz_F2xz+ABZ*I_ESP_F2xy_F2xz;
    Double I_ESP_F2xz_G2x2z = I_ESP_G2x2z_F2xz+ABZ*I_ESP_F2xz_F2xz;
    Double I_ESP_Fx2y_G2x2z = I_ESP_Gx2yz_F2xz+ABZ*I_ESP_Fx2y_F2xz;
    Double I_ESP_Fxyz_G2x2z = I_ESP_Gxy2z_F2xz+ABZ*I_ESP_Fxyz_F2xz;
    Double I_ESP_Fx2z_G2x2z = I_ESP_Gx3z_F2xz+ABZ*I_ESP_Fx2z_F2xz;
    Double I_ESP_F3y_G2x2z = I_ESP_G3yz_F2xz+ABZ*I_ESP_F3y_F2xz;
    Double I_ESP_F2yz_G2x2z = I_ESP_G2y2z_F2xz+ABZ*I_ESP_F2yz_F2xz;
    Double I_ESP_Fy2z_G2x2z = I_ESP_Gy3z_F2xz+ABZ*I_ESP_Fy2z_F2xz;
    Double I_ESP_F3z_G2x2z = I_ESP_G4z_F2xz+ABZ*I_ESP_F3z_F2xz;
    Double I_ESP_F3x_Gx3y = I_ESP_G4x_F3y+ABX*I_ESP_F3x_F3y;
    Double I_ESP_F2xy_Gx3y = I_ESP_G3xy_F3y+ABX*I_ESP_F2xy_F3y;
    Double I_ESP_F2xz_Gx3y = I_ESP_G3xz_F3y+ABX*I_ESP_F2xz_F3y;
    Double I_ESP_Fx2y_Gx3y = I_ESP_G2x2y_F3y+ABX*I_ESP_Fx2y_F3y;
    Double I_ESP_Fxyz_Gx3y = I_ESP_G2xyz_F3y+ABX*I_ESP_Fxyz_F3y;
    Double I_ESP_Fx2z_Gx3y = I_ESP_G2x2z_F3y+ABX*I_ESP_Fx2z_F3y;
    Double I_ESP_F3y_Gx3y = I_ESP_Gx3y_F3y+ABX*I_ESP_F3y_F3y;
    Double I_ESP_F2yz_Gx3y = I_ESP_Gx2yz_F3y+ABX*I_ESP_F2yz_F3y;
    Double I_ESP_Fy2z_Gx3y = I_ESP_Gxy2z_F3y+ABX*I_ESP_Fy2z_F3y;
    Double I_ESP_F3z_Gx3y = I_ESP_Gx3z_F3y+ABX*I_ESP_F3z_F3y;
    Double I_ESP_F3x_Gx2yz = I_ESP_G3xz_Fx2y+ABZ*I_ESP_F3x_Fx2y;
    Double I_ESP_F2xy_Gx2yz = I_ESP_G2xyz_Fx2y+ABZ*I_ESP_F2xy_Fx2y;
    Double I_ESP_F2xz_Gx2yz = I_ESP_G2x2z_Fx2y+ABZ*I_ESP_F2xz_Fx2y;
    Double I_ESP_Fx2y_Gx2yz = I_ESP_Gx2yz_Fx2y+ABZ*I_ESP_Fx2y_Fx2y;
    Double I_ESP_Fxyz_Gx2yz = I_ESP_Gxy2z_Fx2y+ABZ*I_ESP_Fxyz_Fx2y;
    Double I_ESP_Fx2z_Gx2yz = I_ESP_Gx3z_Fx2y+ABZ*I_ESP_Fx2z_Fx2y;
    Double I_ESP_F3y_Gx2yz = I_ESP_G3yz_Fx2y+ABZ*I_ESP_F3y_Fx2y;
    Double I_ESP_F2yz_Gx2yz = I_ESP_G2y2z_Fx2y+ABZ*I_ESP_F2yz_Fx2y;
    Double I_ESP_Fy2z_Gx2yz = I_ESP_Gy3z_Fx2y+ABZ*I_ESP_Fy2z_Fx2y;
    Double I_ESP_F3z_Gx2yz = I_ESP_G4z_Fx2y+ABZ*I_ESP_F3z_Fx2y;
    Double I_ESP_F3x_Gxy2z = I_ESP_G3xy_Fx2z+ABY*I_ESP_F3x_Fx2z;
    Double I_ESP_F2xy_Gxy2z = I_ESP_G2x2y_Fx2z+ABY*I_ESP_F2xy_Fx2z;
    Double I_ESP_F2xz_Gxy2z = I_ESP_G2xyz_Fx2z+ABY*I_ESP_F2xz_Fx2z;
    Double I_ESP_Fx2y_Gxy2z = I_ESP_Gx3y_Fx2z+ABY*I_ESP_Fx2y_Fx2z;
    Double I_ESP_Fxyz_Gxy2z = I_ESP_Gx2yz_Fx2z+ABY*I_ESP_Fxyz_Fx2z;
    Double I_ESP_Fx2z_Gxy2z = I_ESP_Gxy2z_Fx2z+ABY*I_ESP_Fx2z_Fx2z;
    Double I_ESP_F3y_Gxy2z = I_ESP_G4y_Fx2z+ABY*I_ESP_F3y_Fx2z;
    Double I_ESP_F2yz_Gxy2z = I_ESP_G3yz_Fx2z+ABY*I_ESP_F2yz_Fx2z;
    Double I_ESP_Fy2z_Gxy2z = I_ESP_G2y2z_Fx2z+ABY*I_ESP_Fy2z_Fx2z;
    Double I_ESP_F3z_Gxy2z = I_ESP_Gy3z_Fx2z+ABY*I_ESP_F3z_Fx2z;
    Double I_ESP_F3x_Gx3z = I_ESP_G4x_F3z+ABX*I_ESP_F3x_F3z;
    Double I_ESP_F2xy_Gx3z = I_ESP_G3xy_F3z+ABX*I_ESP_F2xy_F3z;
    Double I_ESP_F2xz_Gx3z = I_ESP_G3xz_F3z+ABX*I_ESP_F2xz_F3z;
    Double I_ESP_Fx2y_Gx3z = I_ESP_G2x2y_F3z+ABX*I_ESP_Fx2y_F3z;
    Double I_ESP_Fxyz_Gx3z = I_ESP_G2xyz_F3z+ABX*I_ESP_Fxyz_F3z;
    Double I_ESP_Fx2z_Gx3z = I_ESP_G2x2z_F3z+ABX*I_ESP_Fx2z_F3z;
    Double I_ESP_F3y_Gx3z = I_ESP_Gx3y_F3z+ABX*I_ESP_F3y_F3z;
    Double I_ESP_F2yz_Gx3z = I_ESP_Gx2yz_F3z+ABX*I_ESP_F2yz_F3z;
    Double I_ESP_Fy2z_Gx3z = I_ESP_Gxy2z_F3z+ABX*I_ESP_Fy2z_F3z;
    Double I_ESP_F3z_Gx3z = I_ESP_Gx3z_F3z+ABX*I_ESP_F3z_F3z;
    Double I_ESP_F3x_G4y = I_ESP_G3xy_F3y+ABY*I_ESP_F3x_F3y;
    Double I_ESP_F2xy_G4y = I_ESP_G2x2y_F3y+ABY*I_ESP_F2xy_F3y;
    Double I_ESP_F2xz_G4y = I_ESP_G2xyz_F3y+ABY*I_ESP_F2xz_F3y;
    Double I_ESP_Fx2y_G4y = I_ESP_Gx3y_F3y+ABY*I_ESP_Fx2y_F3y;
    Double I_ESP_Fxyz_G4y = I_ESP_Gx2yz_F3y+ABY*I_ESP_Fxyz_F3y;
    Double I_ESP_Fx2z_G4y = I_ESP_Gxy2z_F3y+ABY*I_ESP_Fx2z_F3y;
    Double I_ESP_F3y_G4y = I_ESP_G4y_F3y+ABY*I_ESP_F3y_F3y;
    Double I_ESP_F2yz_G4y = I_ESP_G3yz_F3y+ABY*I_ESP_F2yz_F3y;
    Double I_ESP_Fy2z_G4y = I_ESP_G2y2z_F3y+ABY*I_ESP_Fy2z_F3y;
    Double I_ESP_F3z_G4y = I_ESP_Gy3z_F3y+ABY*I_ESP_F3z_F3y;
    Double I_ESP_F3x_G3yz = I_ESP_G3xz_F3y+ABZ*I_ESP_F3x_F3y;
    Double I_ESP_F2xy_G3yz = I_ESP_G2xyz_F3y+ABZ*I_ESP_F2xy_F3y;
    Double I_ESP_F2xz_G3yz = I_ESP_G2x2z_F3y+ABZ*I_ESP_F2xz_F3y;
    Double I_ESP_Fx2y_G3yz = I_ESP_Gx2yz_F3y+ABZ*I_ESP_Fx2y_F3y;
    Double I_ESP_Fxyz_G3yz = I_ESP_Gxy2z_F3y+ABZ*I_ESP_Fxyz_F3y;
    Double I_ESP_Fx2z_G3yz = I_ESP_Gx3z_F3y+ABZ*I_ESP_Fx2z_F3y;
    Double I_ESP_F3y_G3yz = I_ESP_G3yz_F3y+ABZ*I_ESP_F3y_F3y;
    Double I_ESP_F2yz_G3yz = I_ESP_G2y2z_F3y+ABZ*I_ESP_F2yz_F3y;
    Double I_ESP_Fy2z_G3yz = I_ESP_Gy3z_F3y+ABZ*I_ESP_Fy2z_F3y;
    Double I_ESP_F3z_G3yz = I_ESP_G4z_F3y+ABZ*I_ESP_F3z_F3y;
    Double I_ESP_F3x_G2y2z = I_ESP_G3xz_F2yz+ABZ*I_ESP_F3x_F2yz;
    Double I_ESP_F2xy_G2y2z = I_ESP_G2xyz_F2yz+ABZ*I_ESP_F2xy_F2yz;
    Double I_ESP_F2xz_G2y2z = I_ESP_G2x2z_F2yz+ABZ*I_ESP_F2xz_F2yz;
    Double I_ESP_Fx2y_G2y2z = I_ESP_Gx2yz_F2yz+ABZ*I_ESP_Fx2y_F2yz;
    Double I_ESP_Fxyz_G2y2z = I_ESP_Gxy2z_F2yz+ABZ*I_ESP_Fxyz_F2yz;
    Double I_ESP_Fx2z_G2y2z = I_ESP_Gx3z_F2yz+ABZ*I_ESP_Fx2z_F2yz;
    Double I_ESP_F3y_G2y2z = I_ESP_G3yz_F2yz+ABZ*I_ESP_F3y_F2yz;
    Double I_ESP_F2yz_G2y2z = I_ESP_G2y2z_F2yz+ABZ*I_ESP_F2yz_F2yz;
    Double I_ESP_Fy2z_G2y2z = I_ESP_Gy3z_F2yz+ABZ*I_ESP_Fy2z_F2yz;
    Double I_ESP_F3z_G2y2z = I_ESP_G4z_F2yz+ABZ*I_ESP_F3z_F2yz;
    Double I_ESP_F3x_Gy3z = I_ESP_G3xy_F3z+ABY*I_ESP_F3x_F3z;
    Double I_ESP_F2xy_Gy3z = I_ESP_G2x2y_F3z+ABY*I_ESP_F2xy_F3z;
    Double I_ESP_F2xz_Gy3z = I_ESP_G2xyz_F3z+ABY*I_ESP_F2xz_F3z;
    Double I_ESP_Fx2y_Gy3z = I_ESP_Gx3y_F3z+ABY*I_ESP_Fx2y_F3z;
    Double I_ESP_Fxyz_Gy3z = I_ESP_Gx2yz_F3z+ABY*I_ESP_Fxyz_F3z;
    Double I_ESP_Fx2z_Gy3z = I_ESP_Gxy2z_F3z+ABY*I_ESP_Fx2z_F3z;
    Double I_ESP_F3y_Gy3z = I_ESP_G4y_F3z+ABY*I_ESP_F3y_F3z;
    Double I_ESP_F2yz_Gy3z = I_ESP_G3yz_F3z+ABY*I_ESP_F2yz_F3z;
    Double I_ESP_Fy2z_Gy3z = I_ESP_G2y2z_F3z+ABY*I_ESP_Fy2z_F3z;
    Double I_ESP_F3z_Gy3z = I_ESP_Gy3z_F3z+ABY*I_ESP_F3z_F3z;
    Double I_ESP_F3x_G4z = I_ESP_G3xz_F3z+ABZ*I_ESP_F3x_F3z;
    Double I_ESP_F2xy_G4z = I_ESP_G2xyz_F3z+ABZ*I_ESP_F2xy_F3z;
    Double I_ESP_F2xz_G4z = I_ESP_G2x2z_F3z+ABZ*I_ESP_F2xz_F3z;
    Double I_ESP_Fx2y_G4z = I_ESP_Gx2yz_F3z+ABZ*I_ESP_Fx2y_F3z;
    Double I_ESP_Fxyz_G4z = I_ESP_Gxy2z_F3z+ABZ*I_ESP_Fxyz_F3z;
    Double I_ESP_Fx2z_G4z = I_ESP_Gx3z_F3z+ABZ*I_ESP_Fx2z_F3z;
    Double I_ESP_F3y_G4z = I_ESP_G3yz_F3z+ABZ*I_ESP_F3y_F3z;
    Double I_ESP_F2yz_G4z = I_ESP_G2y2z_F3z+ABZ*I_ESP_F2yz_F3z;
    Double I_ESP_Fy2z_G4z = I_ESP_Gy3z_F3z+ABZ*I_ESP_Fy2z_F3z;
    Double I_ESP_F3z_G4z = I_ESP_G4z_F3z+ABZ*I_ESP_F3z_F3z;

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
     * totally 0 integrals are omitted 
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
    Double I_ESP_I6x_Py_a = I_ESP_K6xy_S_a+ABY*I_ESP_I6x_S_a;
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
    Double I_ESP_I6x_Pz_a = I_ESP_K6xz_S_a+ABZ*I_ESP_I6x_S_a;
    Double I_ESP_I5xy_Pz_a = I_ESP_K5xyz_S_a+ABZ*I_ESP_I5xy_S_a;
    Double I_ESP_I5xz_Pz_a = I_ESP_K5x2z_S_a+ABZ*I_ESP_I5xz_S_a;
    Double I_ESP_I4x2y_Pz_a = I_ESP_K4x2yz_S_a+ABZ*I_ESP_I4x2y_S_a;
    Double I_ESP_I4xyz_Pz_a = I_ESP_K4xy2z_S_a+ABZ*I_ESP_I4xyz_S_a;
    Double I_ESP_I4x2z_Pz_a = I_ESP_K4x3z_S_a+ABZ*I_ESP_I4x2z_S_a;
    Double I_ESP_I3x3y_Pz_a = I_ESP_K3x3yz_S_a+ABZ*I_ESP_I3x3y_S_a;
    Double I_ESP_I3x2yz_Pz_a = I_ESP_K3x2y2z_S_a+ABZ*I_ESP_I3x2yz_S_a;
    Double I_ESP_I3xy2z_Pz_a = I_ESP_K3xy3z_S_a+ABZ*I_ESP_I3xy2z_S_a;
    Double I_ESP_I3x3z_Pz_a = I_ESP_K3x4z_S_a+ABZ*I_ESP_I3x3z_S_a;
    Double I_ESP_I2x4y_Pz_a = I_ESP_K2x4yz_S_a+ABZ*I_ESP_I2x4y_S_a;
    Double I_ESP_I2x3yz_Pz_a = I_ESP_K2x3y2z_S_a+ABZ*I_ESP_I2x3yz_S_a;
    Double I_ESP_I2x2y2z_Pz_a = I_ESP_K2x2y3z_S_a+ABZ*I_ESP_I2x2y2z_S_a;
    Double I_ESP_I2xy3z_Pz_a = I_ESP_K2xy4z_S_a+ABZ*I_ESP_I2xy3z_S_a;
    Double I_ESP_I2x4z_Pz_a = I_ESP_K2x5z_S_a+ABZ*I_ESP_I2x4z_S_a;
    Double I_ESP_Ix5y_Pz_a = I_ESP_Kx5yz_S_a+ABZ*I_ESP_Ix5y_S_a;
    Double I_ESP_Ix4yz_Pz_a = I_ESP_Kx4y2z_S_a+ABZ*I_ESP_Ix4yz_S_a;
    Double I_ESP_Ix3y2z_Pz_a = I_ESP_Kx3y3z_S_a+ABZ*I_ESP_Ix3y2z_S_a;
    Double I_ESP_Ix2y3z_Pz_a = I_ESP_Kx2y4z_S_a+ABZ*I_ESP_Ix2y3z_S_a;
    Double I_ESP_Ixy4z_Pz_a = I_ESP_Kxy5z_S_a+ABZ*I_ESP_Ixy4z_S_a;
    Double I_ESP_Ix5z_Pz_a = I_ESP_Kx6z_S_a+ABZ*I_ESP_Ix5z_S_a;
    Double I_ESP_I6y_Pz_a = I_ESP_K6yz_S_a+ABZ*I_ESP_I6y_S_a;
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
     * totally 63 integrals are omitted 
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
     * shell quartet name: SQ_ESP_K_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_S_a
     * RHS shell quartet name: SQ_ESP_K_S_a
     ************************************************************/
    Double I_ESP_K7x_Px_a = I_ESP_L8x_S_a+ABX*I_ESP_K7x_S_a;
    Double I_ESP_K6xy_Px_a = I_ESP_L7xy_S_a+ABX*I_ESP_K6xy_S_a;
    Double I_ESP_K6xz_Px_a = I_ESP_L7xz_S_a+ABX*I_ESP_K6xz_S_a;
    Double I_ESP_K5x2y_Px_a = I_ESP_L6x2y_S_a+ABX*I_ESP_K5x2y_S_a;
    Double I_ESP_K5xyz_Px_a = I_ESP_L6xyz_S_a+ABX*I_ESP_K5xyz_S_a;
    Double I_ESP_K5x2z_Px_a = I_ESP_L6x2z_S_a+ABX*I_ESP_K5x2z_S_a;
    Double I_ESP_K4x3y_Px_a = I_ESP_L5x3y_S_a+ABX*I_ESP_K4x3y_S_a;
    Double I_ESP_K4x2yz_Px_a = I_ESP_L5x2yz_S_a+ABX*I_ESP_K4x2yz_S_a;
    Double I_ESP_K4xy2z_Px_a = I_ESP_L5xy2z_S_a+ABX*I_ESP_K4xy2z_S_a;
    Double I_ESP_K4x3z_Px_a = I_ESP_L5x3z_S_a+ABX*I_ESP_K4x3z_S_a;
    Double I_ESP_K3x4y_Px_a = I_ESP_L4x4y_S_a+ABX*I_ESP_K3x4y_S_a;
    Double I_ESP_K3x3yz_Px_a = I_ESP_L4x3yz_S_a+ABX*I_ESP_K3x3yz_S_a;
    Double I_ESP_K3x2y2z_Px_a = I_ESP_L4x2y2z_S_a+ABX*I_ESP_K3x2y2z_S_a;
    Double I_ESP_K3xy3z_Px_a = I_ESP_L4xy3z_S_a+ABX*I_ESP_K3xy3z_S_a;
    Double I_ESP_K3x4z_Px_a = I_ESP_L4x4z_S_a+ABX*I_ESP_K3x4z_S_a;
    Double I_ESP_K2x5y_Px_a = I_ESP_L3x5y_S_a+ABX*I_ESP_K2x5y_S_a;
    Double I_ESP_K2x4yz_Px_a = I_ESP_L3x4yz_S_a+ABX*I_ESP_K2x4yz_S_a;
    Double I_ESP_K2x3y2z_Px_a = I_ESP_L3x3y2z_S_a+ABX*I_ESP_K2x3y2z_S_a;
    Double I_ESP_K2x2y3z_Px_a = I_ESP_L3x2y3z_S_a+ABX*I_ESP_K2x2y3z_S_a;
    Double I_ESP_K2xy4z_Px_a = I_ESP_L3xy4z_S_a+ABX*I_ESP_K2xy4z_S_a;
    Double I_ESP_K2x5z_Px_a = I_ESP_L3x5z_S_a+ABX*I_ESP_K2x5z_S_a;
    Double I_ESP_Kx6y_Px_a = I_ESP_L2x6y_S_a+ABX*I_ESP_Kx6y_S_a;
    Double I_ESP_Kx5yz_Px_a = I_ESP_L2x5yz_S_a+ABX*I_ESP_Kx5yz_S_a;
    Double I_ESP_Kx4y2z_Px_a = I_ESP_L2x4y2z_S_a+ABX*I_ESP_Kx4y2z_S_a;
    Double I_ESP_Kx3y3z_Px_a = I_ESP_L2x3y3z_S_a+ABX*I_ESP_Kx3y3z_S_a;
    Double I_ESP_Kx2y4z_Px_a = I_ESP_L2x2y4z_S_a+ABX*I_ESP_Kx2y4z_S_a;
    Double I_ESP_Kxy5z_Px_a = I_ESP_L2xy5z_S_a+ABX*I_ESP_Kxy5z_S_a;
    Double I_ESP_Kx6z_Px_a = I_ESP_L2x6z_S_a+ABX*I_ESP_Kx6z_S_a;
    Double I_ESP_K7y_Px_a = I_ESP_Lx7y_S_a+ABX*I_ESP_K7y_S_a;
    Double I_ESP_K6yz_Px_a = I_ESP_Lx6yz_S_a+ABX*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Px_a = I_ESP_Lx5y2z_S_a+ABX*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Px_a = I_ESP_Lx4y3z_S_a+ABX*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Px_a = I_ESP_Lx3y4z_S_a+ABX*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Px_a = I_ESP_Lx2y5z_S_a+ABX*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Px_a = I_ESP_Lxy6z_S_a+ABX*I_ESP_Ky6z_S_a;
    Double I_ESP_K7z_Px_a = I_ESP_Lx7z_S_a+ABX*I_ESP_K7z_S_a;
    Double I_ESP_K6xy_Py_a = I_ESP_L6x2y_S_a+ABY*I_ESP_K6xy_S_a;
    Double I_ESP_K6xz_Py_a = I_ESP_L6xyz_S_a+ABY*I_ESP_K6xz_S_a;
    Double I_ESP_K5x2y_Py_a = I_ESP_L5x3y_S_a+ABY*I_ESP_K5x2y_S_a;
    Double I_ESP_K5xyz_Py_a = I_ESP_L5x2yz_S_a+ABY*I_ESP_K5xyz_S_a;
    Double I_ESP_K5x2z_Py_a = I_ESP_L5xy2z_S_a+ABY*I_ESP_K5x2z_S_a;
    Double I_ESP_K4x3y_Py_a = I_ESP_L4x4y_S_a+ABY*I_ESP_K4x3y_S_a;
    Double I_ESP_K4x2yz_Py_a = I_ESP_L4x3yz_S_a+ABY*I_ESP_K4x2yz_S_a;
    Double I_ESP_K4xy2z_Py_a = I_ESP_L4x2y2z_S_a+ABY*I_ESP_K4xy2z_S_a;
    Double I_ESP_K4x3z_Py_a = I_ESP_L4xy3z_S_a+ABY*I_ESP_K4x3z_S_a;
    Double I_ESP_K3x4y_Py_a = I_ESP_L3x5y_S_a+ABY*I_ESP_K3x4y_S_a;
    Double I_ESP_K3x3yz_Py_a = I_ESP_L3x4yz_S_a+ABY*I_ESP_K3x3yz_S_a;
    Double I_ESP_K3x2y2z_Py_a = I_ESP_L3x3y2z_S_a+ABY*I_ESP_K3x2y2z_S_a;
    Double I_ESP_K3xy3z_Py_a = I_ESP_L3x2y3z_S_a+ABY*I_ESP_K3xy3z_S_a;
    Double I_ESP_K3x4z_Py_a = I_ESP_L3xy4z_S_a+ABY*I_ESP_K3x4z_S_a;
    Double I_ESP_K2x5y_Py_a = I_ESP_L2x6y_S_a+ABY*I_ESP_K2x5y_S_a;
    Double I_ESP_K2x4yz_Py_a = I_ESP_L2x5yz_S_a+ABY*I_ESP_K2x4yz_S_a;
    Double I_ESP_K2x3y2z_Py_a = I_ESP_L2x4y2z_S_a+ABY*I_ESP_K2x3y2z_S_a;
    Double I_ESP_K2x2y3z_Py_a = I_ESP_L2x3y3z_S_a+ABY*I_ESP_K2x2y3z_S_a;
    Double I_ESP_K2xy4z_Py_a = I_ESP_L2x2y4z_S_a+ABY*I_ESP_K2xy4z_S_a;
    Double I_ESP_K2x5z_Py_a = I_ESP_L2xy5z_S_a+ABY*I_ESP_K2x5z_S_a;
    Double I_ESP_Kx6y_Py_a = I_ESP_Lx7y_S_a+ABY*I_ESP_Kx6y_S_a;
    Double I_ESP_Kx5yz_Py_a = I_ESP_Lx6yz_S_a+ABY*I_ESP_Kx5yz_S_a;
    Double I_ESP_Kx4y2z_Py_a = I_ESP_Lx5y2z_S_a+ABY*I_ESP_Kx4y2z_S_a;
    Double I_ESP_Kx3y3z_Py_a = I_ESP_Lx4y3z_S_a+ABY*I_ESP_Kx3y3z_S_a;
    Double I_ESP_Kx2y4z_Py_a = I_ESP_Lx3y4z_S_a+ABY*I_ESP_Kx2y4z_S_a;
    Double I_ESP_Kxy5z_Py_a = I_ESP_Lx2y5z_S_a+ABY*I_ESP_Kxy5z_S_a;
    Double I_ESP_Kx6z_Py_a = I_ESP_Lxy6z_S_a+ABY*I_ESP_Kx6z_S_a;
    Double I_ESP_K7y_Py_a = I_ESP_L8y_S_a+ABY*I_ESP_K7y_S_a;
    Double I_ESP_K6yz_Py_a = I_ESP_L7yz_S_a+ABY*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Py_a = I_ESP_L6y2z_S_a+ABY*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Py_a = I_ESP_L5y3z_S_a+ABY*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Py_a = I_ESP_L4y4z_S_a+ABY*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Py_a = I_ESP_L3y5z_S_a+ABY*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Py_a = I_ESP_L2y6z_S_a+ABY*I_ESP_Ky6z_S_a;
    Double I_ESP_K7z_Py_a = I_ESP_Ly7z_S_a+ABY*I_ESP_K7z_S_a;
    Double I_ESP_K6xy_Pz_a = I_ESP_L6xyz_S_a+ABZ*I_ESP_K6xy_S_a;
    Double I_ESP_K6xz_Pz_a = I_ESP_L6x2z_S_a+ABZ*I_ESP_K6xz_S_a;
    Double I_ESP_K5x2y_Pz_a = I_ESP_L5x2yz_S_a+ABZ*I_ESP_K5x2y_S_a;
    Double I_ESP_K5xyz_Pz_a = I_ESP_L5xy2z_S_a+ABZ*I_ESP_K5xyz_S_a;
    Double I_ESP_K5x2z_Pz_a = I_ESP_L5x3z_S_a+ABZ*I_ESP_K5x2z_S_a;
    Double I_ESP_K4x3y_Pz_a = I_ESP_L4x3yz_S_a+ABZ*I_ESP_K4x3y_S_a;
    Double I_ESP_K4x2yz_Pz_a = I_ESP_L4x2y2z_S_a+ABZ*I_ESP_K4x2yz_S_a;
    Double I_ESP_K4xy2z_Pz_a = I_ESP_L4xy3z_S_a+ABZ*I_ESP_K4xy2z_S_a;
    Double I_ESP_K4x3z_Pz_a = I_ESP_L4x4z_S_a+ABZ*I_ESP_K4x3z_S_a;
    Double I_ESP_K3x4y_Pz_a = I_ESP_L3x4yz_S_a+ABZ*I_ESP_K3x4y_S_a;
    Double I_ESP_K3x3yz_Pz_a = I_ESP_L3x3y2z_S_a+ABZ*I_ESP_K3x3yz_S_a;
    Double I_ESP_K3x2y2z_Pz_a = I_ESP_L3x2y3z_S_a+ABZ*I_ESP_K3x2y2z_S_a;
    Double I_ESP_K3xy3z_Pz_a = I_ESP_L3xy4z_S_a+ABZ*I_ESP_K3xy3z_S_a;
    Double I_ESP_K3x4z_Pz_a = I_ESP_L3x5z_S_a+ABZ*I_ESP_K3x4z_S_a;
    Double I_ESP_K2x5y_Pz_a = I_ESP_L2x5yz_S_a+ABZ*I_ESP_K2x5y_S_a;
    Double I_ESP_K2x4yz_Pz_a = I_ESP_L2x4y2z_S_a+ABZ*I_ESP_K2x4yz_S_a;
    Double I_ESP_K2x3y2z_Pz_a = I_ESP_L2x3y3z_S_a+ABZ*I_ESP_K2x3y2z_S_a;
    Double I_ESP_K2x2y3z_Pz_a = I_ESP_L2x2y4z_S_a+ABZ*I_ESP_K2x2y3z_S_a;
    Double I_ESP_K2xy4z_Pz_a = I_ESP_L2xy5z_S_a+ABZ*I_ESP_K2xy4z_S_a;
    Double I_ESP_K2x5z_Pz_a = I_ESP_L2x6z_S_a+ABZ*I_ESP_K2x5z_S_a;
    Double I_ESP_Kx6y_Pz_a = I_ESP_Lx6yz_S_a+ABZ*I_ESP_Kx6y_S_a;
    Double I_ESP_Kx5yz_Pz_a = I_ESP_Lx5y2z_S_a+ABZ*I_ESP_Kx5yz_S_a;
    Double I_ESP_Kx4y2z_Pz_a = I_ESP_Lx4y3z_S_a+ABZ*I_ESP_Kx4y2z_S_a;
    Double I_ESP_Kx3y3z_Pz_a = I_ESP_Lx3y4z_S_a+ABZ*I_ESP_Kx3y3z_S_a;
    Double I_ESP_Kx2y4z_Pz_a = I_ESP_Lx2y5z_S_a+ABZ*I_ESP_Kx2y4z_S_a;
    Double I_ESP_Kxy5z_Pz_a = I_ESP_Lxy6z_S_a+ABZ*I_ESP_Kxy5z_S_a;
    Double I_ESP_Kx6z_Pz_a = I_ESP_Lx7z_S_a+ABZ*I_ESP_Kx6z_S_a;
    Double I_ESP_K6yz_Pz_a = I_ESP_L6y2z_S_a+ABZ*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Pz_a = I_ESP_L5y3z_S_a+ABZ*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Pz_a = I_ESP_L4y4z_S_a+ABZ*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Pz_a = I_ESP_L3y5z_S_a+ABZ*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Pz_a = I_ESP_L2y6z_S_a+ABZ*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Pz_a = I_ESP_Ly7z_S_a+ABZ*I_ESP_Ky6z_S_a;
    Double I_ESP_K7z_Pz_a = I_ESP_L8z_S_a+ABZ*I_ESP_K7z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 84 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P_a
     * RHS shell quartet name: SQ_ESP_I_P_a
     ************************************************************/
    Double I_ESP_I6x_D2x_a = I_ESP_K7x_Px_a+ABX*I_ESP_I6x_Px_a;
    Double I_ESP_I5xy_D2x_a = I_ESP_K6xy_Px_a+ABX*I_ESP_I5xy_Px_a;
    Double I_ESP_I5xz_D2x_a = I_ESP_K6xz_Px_a+ABX*I_ESP_I5xz_Px_a;
    Double I_ESP_I4x2y_D2x_a = I_ESP_K5x2y_Px_a+ABX*I_ESP_I4x2y_Px_a;
    Double I_ESP_I4xyz_D2x_a = I_ESP_K5xyz_Px_a+ABX*I_ESP_I4xyz_Px_a;
    Double I_ESP_I4x2z_D2x_a = I_ESP_K5x2z_Px_a+ABX*I_ESP_I4x2z_Px_a;
    Double I_ESP_I3x3y_D2x_a = I_ESP_K4x3y_Px_a+ABX*I_ESP_I3x3y_Px_a;
    Double I_ESP_I3x2yz_D2x_a = I_ESP_K4x2yz_Px_a+ABX*I_ESP_I3x2yz_Px_a;
    Double I_ESP_I3xy2z_D2x_a = I_ESP_K4xy2z_Px_a+ABX*I_ESP_I3xy2z_Px_a;
    Double I_ESP_I3x3z_D2x_a = I_ESP_K4x3z_Px_a+ABX*I_ESP_I3x3z_Px_a;
    Double I_ESP_I2x4y_D2x_a = I_ESP_K3x4y_Px_a+ABX*I_ESP_I2x4y_Px_a;
    Double I_ESP_I2x3yz_D2x_a = I_ESP_K3x3yz_Px_a+ABX*I_ESP_I2x3yz_Px_a;
    Double I_ESP_I2x2y2z_D2x_a = I_ESP_K3x2y2z_Px_a+ABX*I_ESP_I2x2y2z_Px_a;
    Double I_ESP_I2xy3z_D2x_a = I_ESP_K3xy3z_Px_a+ABX*I_ESP_I2xy3z_Px_a;
    Double I_ESP_I2x4z_D2x_a = I_ESP_K3x4z_Px_a+ABX*I_ESP_I2x4z_Px_a;
    Double I_ESP_Ix5y_D2x_a = I_ESP_K2x5y_Px_a+ABX*I_ESP_Ix5y_Px_a;
    Double I_ESP_Ix4yz_D2x_a = I_ESP_K2x4yz_Px_a+ABX*I_ESP_Ix4yz_Px_a;
    Double I_ESP_Ix3y2z_D2x_a = I_ESP_K2x3y2z_Px_a+ABX*I_ESP_Ix3y2z_Px_a;
    Double I_ESP_Ix2y3z_D2x_a = I_ESP_K2x2y3z_Px_a+ABX*I_ESP_Ix2y3z_Px_a;
    Double I_ESP_Ixy4z_D2x_a = I_ESP_K2xy4z_Px_a+ABX*I_ESP_Ixy4z_Px_a;
    Double I_ESP_Ix5z_D2x_a = I_ESP_K2x5z_Px_a+ABX*I_ESP_Ix5z_Px_a;
    Double I_ESP_I6y_D2x_a = I_ESP_Kx6y_Px_a+ABX*I_ESP_I6y_Px_a;
    Double I_ESP_I5yz_D2x_a = I_ESP_Kx5yz_Px_a+ABX*I_ESP_I5yz_Px_a;
    Double I_ESP_I4y2z_D2x_a = I_ESP_Kx4y2z_Px_a+ABX*I_ESP_I4y2z_Px_a;
    Double I_ESP_I3y3z_D2x_a = I_ESP_Kx3y3z_Px_a+ABX*I_ESP_I3y3z_Px_a;
    Double I_ESP_I2y4z_D2x_a = I_ESP_Kx2y4z_Px_a+ABX*I_ESP_I2y4z_Px_a;
    Double I_ESP_Iy5z_D2x_a = I_ESP_Kxy5z_Px_a+ABX*I_ESP_Iy5z_Px_a;
    Double I_ESP_I6z_D2x_a = I_ESP_Kx6z_Px_a+ABX*I_ESP_I6z_Px_a;
    Double I_ESP_I6x_D2y_a = I_ESP_K6xy_Py_a+ABY*I_ESP_I6x_Py_a;
    Double I_ESP_I5xy_D2y_a = I_ESP_K5x2y_Py_a+ABY*I_ESP_I5xy_Py_a;
    Double I_ESP_I5xz_D2y_a = I_ESP_K5xyz_Py_a+ABY*I_ESP_I5xz_Py_a;
    Double I_ESP_I4x2y_D2y_a = I_ESP_K4x3y_Py_a+ABY*I_ESP_I4x2y_Py_a;
    Double I_ESP_I4xyz_D2y_a = I_ESP_K4x2yz_Py_a+ABY*I_ESP_I4xyz_Py_a;
    Double I_ESP_I4x2z_D2y_a = I_ESP_K4xy2z_Py_a+ABY*I_ESP_I4x2z_Py_a;
    Double I_ESP_I3x3y_D2y_a = I_ESP_K3x4y_Py_a+ABY*I_ESP_I3x3y_Py_a;
    Double I_ESP_I3x2yz_D2y_a = I_ESP_K3x3yz_Py_a+ABY*I_ESP_I3x2yz_Py_a;
    Double I_ESP_I3xy2z_D2y_a = I_ESP_K3x2y2z_Py_a+ABY*I_ESP_I3xy2z_Py_a;
    Double I_ESP_I3x3z_D2y_a = I_ESP_K3xy3z_Py_a+ABY*I_ESP_I3x3z_Py_a;
    Double I_ESP_I2x4y_D2y_a = I_ESP_K2x5y_Py_a+ABY*I_ESP_I2x4y_Py_a;
    Double I_ESP_I2x3yz_D2y_a = I_ESP_K2x4yz_Py_a+ABY*I_ESP_I2x3yz_Py_a;
    Double I_ESP_I2x2y2z_D2y_a = I_ESP_K2x3y2z_Py_a+ABY*I_ESP_I2x2y2z_Py_a;
    Double I_ESP_I2xy3z_D2y_a = I_ESP_K2x2y3z_Py_a+ABY*I_ESP_I2xy3z_Py_a;
    Double I_ESP_I2x4z_D2y_a = I_ESP_K2xy4z_Py_a+ABY*I_ESP_I2x4z_Py_a;
    Double I_ESP_Ix5y_D2y_a = I_ESP_Kx6y_Py_a+ABY*I_ESP_Ix5y_Py_a;
    Double I_ESP_Ix4yz_D2y_a = I_ESP_Kx5yz_Py_a+ABY*I_ESP_Ix4yz_Py_a;
    Double I_ESP_Ix3y2z_D2y_a = I_ESP_Kx4y2z_Py_a+ABY*I_ESP_Ix3y2z_Py_a;
    Double I_ESP_Ix2y3z_D2y_a = I_ESP_Kx3y3z_Py_a+ABY*I_ESP_Ix2y3z_Py_a;
    Double I_ESP_Ixy4z_D2y_a = I_ESP_Kx2y4z_Py_a+ABY*I_ESP_Ixy4z_Py_a;
    Double I_ESP_Ix5z_D2y_a = I_ESP_Kxy5z_Py_a+ABY*I_ESP_Ix5z_Py_a;
    Double I_ESP_I6y_D2y_a = I_ESP_K7y_Py_a+ABY*I_ESP_I6y_Py_a;
    Double I_ESP_I5yz_D2y_a = I_ESP_K6yz_Py_a+ABY*I_ESP_I5yz_Py_a;
    Double I_ESP_I4y2z_D2y_a = I_ESP_K5y2z_Py_a+ABY*I_ESP_I4y2z_Py_a;
    Double I_ESP_I3y3z_D2y_a = I_ESP_K4y3z_Py_a+ABY*I_ESP_I3y3z_Py_a;
    Double I_ESP_I2y4z_D2y_a = I_ESP_K3y4z_Py_a+ABY*I_ESP_I2y4z_Py_a;
    Double I_ESP_Iy5z_D2y_a = I_ESP_K2y5z_Py_a+ABY*I_ESP_Iy5z_Py_a;
    Double I_ESP_I6z_D2y_a = I_ESP_Ky6z_Py_a+ABY*I_ESP_I6z_Py_a;
    Double I_ESP_I6x_D2z_a = I_ESP_K6xz_Pz_a+ABZ*I_ESP_I6x_Pz_a;
    Double I_ESP_I5xy_D2z_a = I_ESP_K5xyz_Pz_a+ABZ*I_ESP_I5xy_Pz_a;
    Double I_ESP_I5xz_D2z_a = I_ESP_K5x2z_Pz_a+ABZ*I_ESP_I5xz_Pz_a;
    Double I_ESP_I4x2y_D2z_a = I_ESP_K4x2yz_Pz_a+ABZ*I_ESP_I4x2y_Pz_a;
    Double I_ESP_I4xyz_D2z_a = I_ESP_K4xy2z_Pz_a+ABZ*I_ESP_I4xyz_Pz_a;
    Double I_ESP_I4x2z_D2z_a = I_ESP_K4x3z_Pz_a+ABZ*I_ESP_I4x2z_Pz_a;
    Double I_ESP_I3x3y_D2z_a = I_ESP_K3x3yz_Pz_a+ABZ*I_ESP_I3x3y_Pz_a;
    Double I_ESP_I3x2yz_D2z_a = I_ESP_K3x2y2z_Pz_a+ABZ*I_ESP_I3x2yz_Pz_a;
    Double I_ESP_I3xy2z_D2z_a = I_ESP_K3xy3z_Pz_a+ABZ*I_ESP_I3xy2z_Pz_a;
    Double I_ESP_I3x3z_D2z_a = I_ESP_K3x4z_Pz_a+ABZ*I_ESP_I3x3z_Pz_a;
    Double I_ESP_I2x4y_D2z_a = I_ESP_K2x4yz_Pz_a+ABZ*I_ESP_I2x4y_Pz_a;
    Double I_ESP_I2x3yz_D2z_a = I_ESP_K2x3y2z_Pz_a+ABZ*I_ESP_I2x3yz_Pz_a;
    Double I_ESP_I2x2y2z_D2z_a = I_ESP_K2x2y3z_Pz_a+ABZ*I_ESP_I2x2y2z_Pz_a;
    Double I_ESP_I2xy3z_D2z_a = I_ESP_K2xy4z_Pz_a+ABZ*I_ESP_I2xy3z_Pz_a;
    Double I_ESP_I2x4z_D2z_a = I_ESP_K2x5z_Pz_a+ABZ*I_ESP_I2x4z_Pz_a;
    Double I_ESP_Ix5y_D2z_a = I_ESP_Kx5yz_Pz_a+ABZ*I_ESP_Ix5y_Pz_a;
    Double I_ESP_Ix4yz_D2z_a = I_ESP_Kx4y2z_Pz_a+ABZ*I_ESP_Ix4yz_Pz_a;
    Double I_ESP_Ix3y2z_D2z_a = I_ESP_Kx3y3z_Pz_a+ABZ*I_ESP_Ix3y2z_Pz_a;
    Double I_ESP_Ix2y3z_D2z_a = I_ESP_Kx2y4z_Pz_a+ABZ*I_ESP_Ix2y3z_Pz_a;
    Double I_ESP_Ixy4z_D2z_a = I_ESP_Kxy5z_Pz_a+ABZ*I_ESP_Ixy4z_Pz_a;
    Double I_ESP_Ix5z_D2z_a = I_ESP_Kx6z_Pz_a+ABZ*I_ESP_Ix5z_Pz_a;
    Double I_ESP_I6y_D2z_a = I_ESP_K6yz_Pz_a+ABZ*I_ESP_I6y_Pz_a;
    Double I_ESP_I5yz_D2z_a = I_ESP_K5y2z_Pz_a+ABZ*I_ESP_I5yz_Pz_a;
    Double I_ESP_I4y2z_D2z_a = I_ESP_K4y3z_Pz_a+ABZ*I_ESP_I4y2z_Pz_a;
    Double I_ESP_I3y3z_D2z_a = I_ESP_K3y4z_Pz_a+ABZ*I_ESP_I3y3z_Pz_a;
    Double I_ESP_I2y4z_D2z_a = I_ESP_K2y5z_Pz_a+ABZ*I_ESP_I2y4z_Pz_a;
    Double I_ESP_Iy5z_D2z_a = I_ESP_Ky6z_Pz_a+ABZ*I_ESP_Iy5z_Pz_a;
    Double I_ESP_I6z_D2z_a = I_ESP_K7z_Pz_a+ABZ*I_ESP_I6z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 42 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D_a
     * RHS shell quartet name: SQ_ESP_H_D_a
     ************************************************************/
    Double I_ESP_H5x_F3x_a = I_ESP_I6x_D2x_a+ABX*I_ESP_H5x_D2x_a;
    Double I_ESP_H4xy_F3x_a = I_ESP_I5xy_D2x_a+ABX*I_ESP_H4xy_D2x_a;
    Double I_ESP_H4xz_F3x_a = I_ESP_I5xz_D2x_a+ABX*I_ESP_H4xz_D2x_a;
    Double I_ESP_H3x2y_F3x_a = I_ESP_I4x2y_D2x_a+ABX*I_ESP_H3x2y_D2x_a;
    Double I_ESP_H3xyz_F3x_a = I_ESP_I4xyz_D2x_a+ABX*I_ESP_H3xyz_D2x_a;
    Double I_ESP_H3x2z_F3x_a = I_ESP_I4x2z_D2x_a+ABX*I_ESP_H3x2z_D2x_a;
    Double I_ESP_H2x3y_F3x_a = I_ESP_I3x3y_D2x_a+ABX*I_ESP_H2x3y_D2x_a;
    Double I_ESP_H2x2yz_F3x_a = I_ESP_I3x2yz_D2x_a+ABX*I_ESP_H2x2yz_D2x_a;
    Double I_ESP_H2xy2z_F3x_a = I_ESP_I3xy2z_D2x_a+ABX*I_ESP_H2xy2z_D2x_a;
    Double I_ESP_H2x3z_F3x_a = I_ESP_I3x3z_D2x_a+ABX*I_ESP_H2x3z_D2x_a;
    Double I_ESP_Hx4y_F3x_a = I_ESP_I2x4y_D2x_a+ABX*I_ESP_Hx4y_D2x_a;
    Double I_ESP_Hx3yz_F3x_a = I_ESP_I2x3yz_D2x_a+ABX*I_ESP_Hx3yz_D2x_a;
    Double I_ESP_Hx2y2z_F3x_a = I_ESP_I2x2y2z_D2x_a+ABX*I_ESP_Hx2y2z_D2x_a;
    Double I_ESP_Hxy3z_F3x_a = I_ESP_I2xy3z_D2x_a+ABX*I_ESP_Hxy3z_D2x_a;
    Double I_ESP_Hx4z_F3x_a = I_ESP_I2x4z_D2x_a+ABX*I_ESP_Hx4z_D2x_a;
    Double I_ESP_H5y_F3x_a = I_ESP_Ix5y_D2x_a+ABX*I_ESP_H5y_D2x_a;
    Double I_ESP_H4yz_F3x_a = I_ESP_Ix4yz_D2x_a+ABX*I_ESP_H4yz_D2x_a;
    Double I_ESP_H3y2z_F3x_a = I_ESP_Ix3y2z_D2x_a+ABX*I_ESP_H3y2z_D2x_a;
    Double I_ESP_H2y3z_F3x_a = I_ESP_Ix2y3z_D2x_a+ABX*I_ESP_H2y3z_D2x_a;
    Double I_ESP_Hy4z_F3x_a = I_ESP_Ixy4z_D2x_a+ABX*I_ESP_Hy4z_D2x_a;
    Double I_ESP_H5z_F3x_a = I_ESP_Ix5z_D2x_a+ABX*I_ESP_H5z_D2x_a;
    Double I_ESP_H5x_F2xy_a = I_ESP_I5xy_D2x_a+ABY*I_ESP_H5x_D2x_a;
    Double I_ESP_H4xy_F2xy_a = I_ESP_I4x2y_D2x_a+ABY*I_ESP_H4xy_D2x_a;
    Double I_ESP_H4xz_F2xy_a = I_ESP_I4xyz_D2x_a+ABY*I_ESP_H4xz_D2x_a;
    Double I_ESP_H3x2y_F2xy_a = I_ESP_I3x3y_D2x_a+ABY*I_ESP_H3x2y_D2x_a;
    Double I_ESP_H3xyz_F2xy_a = I_ESP_I3x2yz_D2x_a+ABY*I_ESP_H3xyz_D2x_a;
    Double I_ESP_H3x2z_F2xy_a = I_ESP_I3xy2z_D2x_a+ABY*I_ESP_H3x2z_D2x_a;
    Double I_ESP_H2x3y_F2xy_a = I_ESP_I2x4y_D2x_a+ABY*I_ESP_H2x3y_D2x_a;
    Double I_ESP_H2x2yz_F2xy_a = I_ESP_I2x3yz_D2x_a+ABY*I_ESP_H2x2yz_D2x_a;
    Double I_ESP_H2xy2z_F2xy_a = I_ESP_I2x2y2z_D2x_a+ABY*I_ESP_H2xy2z_D2x_a;
    Double I_ESP_H2x3z_F2xy_a = I_ESP_I2xy3z_D2x_a+ABY*I_ESP_H2x3z_D2x_a;
    Double I_ESP_Hx4y_F2xy_a = I_ESP_Ix5y_D2x_a+ABY*I_ESP_Hx4y_D2x_a;
    Double I_ESP_Hx3yz_F2xy_a = I_ESP_Ix4yz_D2x_a+ABY*I_ESP_Hx3yz_D2x_a;
    Double I_ESP_Hx2y2z_F2xy_a = I_ESP_Ix3y2z_D2x_a+ABY*I_ESP_Hx2y2z_D2x_a;
    Double I_ESP_Hxy3z_F2xy_a = I_ESP_Ix2y3z_D2x_a+ABY*I_ESP_Hxy3z_D2x_a;
    Double I_ESP_Hx4z_F2xy_a = I_ESP_Ixy4z_D2x_a+ABY*I_ESP_Hx4z_D2x_a;
    Double I_ESP_H5y_F2xy_a = I_ESP_I6y_D2x_a+ABY*I_ESP_H5y_D2x_a;
    Double I_ESP_H4yz_F2xy_a = I_ESP_I5yz_D2x_a+ABY*I_ESP_H4yz_D2x_a;
    Double I_ESP_H3y2z_F2xy_a = I_ESP_I4y2z_D2x_a+ABY*I_ESP_H3y2z_D2x_a;
    Double I_ESP_H2y3z_F2xy_a = I_ESP_I3y3z_D2x_a+ABY*I_ESP_H2y3z_D2x_a;
    Double I_ESP_Hy4z_F2xy_a = I_ESP_I2y4z_D2x_a+ABY*I_ESP_Hy4z_D2x_a;
    Double I_ESP_H5z_F2xy_a = I_ESP_Iy5z_D2x_a+ABY*I_ESP_H5z_D2x_a;
    Double I_ESP_H5x_F2xz_a = I_ESP_I5xz_D2x_a+ABZ*I_ESP_H5x_D2x_a;
    Double I_ESP_H4xy_F2xz_a = I_ESP_I4xyz_D2x_a+ABZ*I_ESP_H4xy_D2x_a;
    Double I_ESP_H4xz_F2xz_a = I_ESP_I4x2z_D2x_a+ABZ*I_ESP_H4xz_D2x_a;
    Double I_ESP_H3x2y_F2xz_a = I_ESP_I3x2yz_D2x_a+ABZ*I_ESP_H3x2y_D2x_a;
    Double I_ESP_H3xyz_F2xz_a = I_ESP_I3xy2z_D2x_a+ABZ*I_ESP_H3xyz_D2x_a;
    Double I_ESP_H3x2z_F2xz_a = I_ESP_I3x3z_D2x_a+ABZ*I_ESP_H3x2z_D2x_a;
    Double I_ESP_H2x3y_F2xz_a = I_ESP_I2x3yz_D2x_a+ABZ*I_ESP_H2x3y_D2x_a;
    Double I_ESP_H2x2yz_F2xz_a = I_ESP_I2x2y2z_D2x_a+ABZ*I_ESP_H2x2yz_D2x_a;
    Double I_ESP_H2xy2z_F2xz_a = I_ESP_I2xy3z_D2x_a+ABZ*I_ESP_H2xy2z_D2x_a;
    Double I_ESP_H2x3z_F2xz_a = I_ESP_I2x4z_D2x_a+ABZ*I_ESP_H2x3z_D2x_a;
    Double I_ESP_Hx4y_F2xz_a = I_ESP_Ix4yz_D2x_a+ABZ*I_ESP_Hx4y_D2x_a;
    Double I_ESP_Hx3yz_F2xz_a = I_ESP_Ix3y2z_D2x_a+ABZ*I_ESP_Hx3yz_D2x_a;
    Double I_ESP_Hx2y2z_F2xz_a = I_ESP_Ix2y3z_D2x_a+ABZ*I_ESP_Hx2y2z_D2x_a;
    Double I_ESP_Hxy3z_F2xz_a = I_ESP_Ixy4z_D2x_a+ABZ*I_ESP_Hxy3z_D2x_a;
    Double I_ESP_Hx4z_F2xz_a = I_ESP_Ix5z_D2x_a+ABZ*I_ESP_Hx4z_D2x_a;
    Double I_ESP_H5y_F2xz_a = I_ESP_I5yz_D2x_a+ABZ*I_ESP_H5y_D2x_a;
    Double I_ESP_H4yz_F2xz_a = I_ESP_I4y2z_D2x_a+ABZ*I_ESP_H4yz_D2x_a;
    Double I_ESP_H3y2z_F2xz_a = I_ESP_I3y3z_D2x_a+ABZ*I_ESP_H3y2z_D2x_a;
    Double I_ESP_H2y3z_F2xz_a = I_ESP_I2y4z_D2x_a+ABZ*I_ESP_H2y3z_D2x_a;
    Double I_ESP_Hy4z_F2xz_a = I_ESP_Iy5z_D2x_a+ABZ*I_ESP_Hy4z_D2x_a;
    Double I_ESP_H5z_F2xz_a = I_ESP_I6z_D2x_a+ABZ*I_ESP_H5z_D2x_a;
    Double I_ESP_H5x_Fx2y_a = I_ESP_I6x_D2y_a+ABX*I_ESP_H5x_D2y_a;
    Double I_ESP_H4xy_Fx2y_a = I_ESP_I5xy_D2y_a+ABX*I_ESP_H4xy_D2y_a;
    Double I_ESP_H4xz_Fx2y_a = I_ESP_I5xz_D2y_a+ABX*I_ESP_H4xz_D2y_a;
    Double I_ESP_H3x2y_Fx2y_a = I_ESP_I4x2y_D2y_a+ABX*I_ESP_H3x2y_D2y_a;
    Double I_ESP_H3xyz_Fx2y_a = I_ESP_I4xyz_D2y_a+ABX*I_ESP_H3xyz_D2y_a;
    Double I_ESP_H3x2z_Fx2y_a = I_ESP_I4x2z_D2y_a+ABX*I_ESP_H3x2z_D2y_a;
    Double I_ESP_H2x3y_Fx2y_a = I_ESP_I3x3y_D2y_a+ABX*I_ESP_H2x3y_D2y_a;
    Double I_ESP_H2x2yz_Fx2y_a = I_ESP_I3x2yz_D2y_a+ABX*I_ESP_H2x2yz_D2y_a;
    Double I_ESP_H2xy2z_Fx2y_a = I_ESP_I3xy2z_D2y_a+ABX*I_ESP_H2xy2z_D2y_a;
    Double I_ESP_H2x3z_Fx2y_a = I_ESP_I3x3z_D2y_a+ABX*I_ESP_H2x3z_D2y_a;
    Double I_ESP_Hx4y_Fx2y_a = I_ESP_I2x4y_D2y_a+ABX*I_ESP_Hx4y_D2y_a;
    Double I_ESP_Hx3yz_Fx2y_a = I_ESP_I2x3yz_D2y_a+ABX*I_ESP_Hx3yz_D2y_a;
    Double I_ESP_Hx2y2z_Fx2y_a = I_ESP_I2x2y2z_D2y_a+ABX*I_ESP_Hx2y2z_D2y_a;
    Double I_ESP_Hxy3z_Fx2y_a = I_ESP_I2xy3z_D2y_a+ABX*I_ESP_Hxy3z_D2y_a;
    Double I_ESP_Hx4z_Fx2y_a = I_ESP_I2x4z_D2y_a+ABX*I_ESP_Hx4z_D2y_a;
    Double I_ESP_H5y_Fx2y_a = I_ESP_Ix5y_D2y_a+ABX*I_ESP_H5y_D2y_a;
    Double I_ESP_H4yz_Fx2y_a = I_ESP_Ix4yz_D2y_a+ABX*I_ESP_H4yz_D2y_a;
    Double I_ESP_H3y2z_Fx2y_a = I_ESP_Ix3y2z_D2y_a+ABX*I_ESP_H3y2z_D2y_a;
    Double I_ESP_H2y3z_Fx2y_a = I_ESP_Ix2y3z_D2y_a+ABX*I_ESP_H2y3z_D2y_a;
    Double I_ESP_Hy4z_Fx2y_a = I_ESP_Ixy4z_D2y_a+ABX*I_ESP_Hy4z_D2y_a;
    Double I_ESP_H5z_Fx2y_a = I_ESP_Ix5z_D2y_a+ABX*I_ESP_H5z_D2y_a;
    Double I_ESP_H5x_Fx2z_a = I_ESP_I6x_D2z_a+ABX*I_ESP_H5x_D2z_a;
    Double I_ESP_H4xy_Fx2z_a = I_ESP_I5xy_D2z_a+ABX*I_ESP_H4xy_D2z_a;
    Double I_ESP_H4xz_Fx2z_a = I_ESP_I5xz_D2z_a+ABX*I_ESP_H4xz_D2z_a;
    Double I_ESP_H3x2y_Fx2z_a = I_ESP_I4x2y_D2z_a+ABX*I_ESP_H3x2y_D2z_a;
    Double I_ESP_H3xyz_Fx2z_a = I_ESP_I4xyz_D2z_a+ABX*I_ESP_H3xyz_D2z_a;
    Double I_ESP_H3x2z_Fx2z_a = I_ESP_I4x2z_D2z_a+ABX*I_ESP_H3x2z_D2z_a;
    Double I_ESP_H2x3y_Fx2z_a = I_ESP_I3x3y_D2z_a+ABX*I_ESP_H2x3y_D2z_a;
    Double I_ESP_H2x2yz_Fx2z_a = I_ESP_I3x2yz_D2z_a+ABX*I_ESP_H2x2yz_D2z_a;
    Double I_ESP_H2xy2z_Fx2z_a = I_ESP_I3xy2z_D2z_a+ABX*I_ESP_H2xy2z_D2z_a;
    Double I_ESP_H2x3z_Fx2z_a = I_ESP_I3x3z_D2z_a+ABX*I_ESP_H2x3z_D2z_a;
    Double I_ESP_Hx4y_Fx2z_a = I_ESP_I2x4y_D2z_a+ABX*I_ESP_Hx4y_D2z_a;
    Double I_ESP_Hx3yz_Fx2z_a = I_ESP_I2x3yz_D2z_a+ABX*I_ESP_Hx3yz_D2z_a;
    Double I_ESP_Hx2y2z_Fx2z_a = I_ESP_I2x2y2z_D2z_a+ABX*I_ESP_Hx2y2z_D2z_a;
    Double I_ESP_Hxy3z_Fx2z_a = I_ESP_I2xy3z_D2z_a+ABX*I_ESP_Hxy3z_D2z_a;
    Double I_ESP_Hx4z_Fx2z_a = I_ESP_I2x4z_D2z_a+ABX*I_ESP_Hx4z_D2z_a;
    Double I_ESP_H5y_Fx2z_a = I_ESP_Ix5y_D2z_a+ABX*I_ESP_H5y_D2z_a;
    Double I_ESP_H4yz_Fx2z_a = I_ESP_Ix4yz_D2z_a+ABX*I_ESP_H4yz_D2z_a;
    Double I_ESP_H3y2z_Fx2z_a = I_ESP_Ix3y2z_D2z_a+ABX*I_ESP_H3y2z_D2z_a;
    Double I_ESP_H2y3z_Fx2z_a = I_ESP_Ix2y3z_D2z_a+ABX*I_ESP_H2y3z_D2z_a;
    Double I_ESP_Hy4z_Fx2z_a = I_ESP_Ixy4z_D2z_a+ABX*I_ESP_Hy4z_D2z_a;
    Double I_ESP_H5z_Fx2z_a = I_ESP_Ix5z_D2z_a+ABX*I_ESP_H5z_D2z_a;
    Double I_ESP_H5x_F3y_a = I_ESP_I5xy_D2y_a+ABY*I_ESP_H5x_D2y_a;
    Double I_ESP_H4xy_F3y_a = I_ESP_I4x2y_D2y_a+ABY*I_ESP_H4xy_D2y_a;
    Double I_ESP_H4xz_F3y_a = I_ESP_I4xyz_D2y_a+ABY*I_ESP_H4xz_D2y_a;
    Double I_ESP_H3x2y_F3y_a = I_ESP_I3x3y_D2y_a+ABY*I_ESP_H3x2y_D2y_a;
    Double I_ESP_H3xyz_F3y_a = I_ESP_I3x2yz_D2y_a+ABY*I_ESP_H3xyz_D2y_a;
    Double I_ESP_H3x2z_F3y_a = I_ESP_I3xy2z_D2y_a+ABY*I_ESP_H3x2z_D2y_a;
    Double I_ESP_H2x3y_F3y_a = I_ESP_I2x4y_D2y_a+ABY*I_ESP_H2x3y_D2y_a;
    Double I_ESP_H2x2yz_F3y_a = I_ESP_I2x3yz_D2y_a+ABY*I_ESP_H2x2yz_D2y_a;
    Double I_ESP_H2xy2z_F3y_a = I_ESP_I2x2y2z_D2y_a+ABY*I_ESP_H2xy2z_D2y_a;
    Double I_ESP_H2x3z_F3y_a = I_ESP_I2xy3z_D2y_a+ABY*I_ESP_H2x3z_D2y_a;
    Double I_ESP_Hx4y_F3y_a = I_ESP_Ix5y_D2y_a+ABY*I_ESP_Hx4y_D2y_a;
    Double I_ESP_Hx3yz_F3y_a = I_ESP_Ix4yz_D2y_a+ABY*I_ESP_Hx3yz_D2y_a;
    Double I_ESP_Hx2y2z_F3y_a = I_ESP_Ix3y2z_D2y_a+ABY*I_ESP_Hx2y2z_D2y_a;
    Double I_ESP_Hxy3z_F3y_a = I_ESP_Ix2y3z_D2y_a+ABY*I_ESP_Hxy3z_D2y_a;
    Double I_ESP_Hx4z_F3y_a = I_ESP_Ixy4z_D2y_a+ABY*I_ESP_Hx4z_D2y_a;
    Double I_ESP_H5y_F3y_a = I_ESP_I6y_D2y_a+ABY*I_ESP_H5y_D2y_a;
    Double I_ESP_H4yz_F3y_a = I_ESP_I5yz_D2y_a+ABY*I_ESP_H4yz_D2y_a;
    Double I_ESP_H3y2z_F3y_a = I_ESP_I4y2z_D2y_a+ABY*I_ESP_H3y2z_D2y_a;
    Double I_ESP_H2y3z_F3y_a = I_ESP_I3y3z_D2y_a+ABY*I_ESP_H2y3z_D2y_a;
    Double I_ESP_Hy4z_F3y_a = I_ESP_I2y4z_D2y_a+ABY*I_ESP_Hy4z_D2y_a;
    Double I_ESP_H5z_F3y_a = I_ESP_Iy5z_D2y_a+ABY*I_ESP_H5z_D2y_a;
    Double I_ESP_H5x_F2yz_a = I_ESP_I5xz_D2y_a+ABZ*I_ESP_H5x_D2y_a;
    Double I_ESP_H4xy_F2yz_a = I_ESP_I4xyz_D2y_a+ABZ*I_ESP_H4xy_D2y_a;
    Double I_ESP_H4xz_F2yz_a = I_ESP_I4x2z_D2y_a+ABZ*I_ESP_H4xz_D2y_a;
    Double I_ESP_H3x2y_F2yz_a = I_ESP_I3x2yz_D2y_a+ABZ*I_ESP_H3x2y_D2y_a;
    Double I_ESP_H3xyz_F2yz_a = I_ESP_I3xy2z_D2y_a+ABZ*I_ESP_H3xyz_D2y_a;
    Double I_ESP_H3x2z_F2yz_a = I_ESP_I3x3z_D2y_a+ABZ*I_ESP_H3x2z_D2y_a;
    Double I_ESP_H2x3y_F2yz_a = I_ESP_I2x3yz_D2y_a+ABZ*I_ESP_H2x3y_D2y_a;
    Double I_ESP_H2x2yz_F2yz_a = I_ESP_I2x2y2z_D2y_a+ABZ*I_ESP_H2x2yz_D2y_a;
    Double I_ESP_H2xy2z_F2yz_a = I_ESP_I2xy3z_D2y_a+ABZ*I_ESP_H2xy2z_D2y_a;
    Double I_ESP_H2x3z_F2yz_a = I_ESP_I2x4z_D2y_a+ABZ*I_ESP_H2x3z_D2y_a;
    Double I_ESP_Hx4y_F2yz_a = I_ESP_Ix4yz_D2y_a+ABZ*I_ESP_Hx4y_D2y_a;
    Double I_ESP_Hx3yz_F2yz_a = I_ESP_Ix3y2z_D2y_a+ABZ*I_ESP_Hx3yz_D2y_a;
    Double I_ESP_Hx2y2z_F2yz_a = I_ESP_Ix2y3z_D2y_a+ABZ*I_ESP_Hx2y2z_D2y_a;
    Double I_ESP_Hxy3z_F2yz_a = I_ESP_Ixy4z_D2y_a+ABZ*I_ESP_Hxy3z_D2y_a;
    Double I_ESP_Hx4z_F2yz_a = I_ESP_Ix5z_D2y_a+ABZ*I_ESP_Hx4z_D2y_a;
    Double I_ESP_H5y_F2yz_a = I_ESP_I5yz_D2y_a+ABZ*I_ESP_H5y_D2y_a;
    Double I_ESP_H4yz_F2yz_a = I_ESP_I4y2z_D2y_a+ABZ*I_ESP_H4yz_D2y_a;
    Double I_ESP_H3y2z_F2yz_a = I_ESP_I3y3z_D2y_a+ABZ*I_ESP_H3y2z_D2y_a;
    Double I_ESP_H2y3z_F2yz_a = I_ESP_I2y4z_D2y_a+ABZ*I_ESP_H2y3z_D2y_a;
    Double I_ESP_Hy4z_F2yz_a = I_ESP_Iy5z_D2y_a+ABZ*I_ESP_Hy4z_D2y_a;
    Double I_ESP_H5z_F2yz_a = I_ESP_I6z_D2y_a+ABZ*I_ESP_H5z_D2y_a;
    Double I_ESP_H5x_F3z_a = I_ESP_I5xz_D2z_a+ABZ*I_ESP_H5x_D2z_a;
    Double I_ESP_H4xy_F3z_a = I_ESP_I4xyz_D2z_a+ABZ*I_ESP_H4xy_D2z_a;
    Double I_ESP_H4xz_F3z_a = I_ESP_I4x2z_D2z_a+ABZ*I_ESP_H4xz_D2z_a;
    Double I_ESP_H3x2y_F3z_a = I_ESP_I3x2yz_D2z_a+ABZ*I_ESP_H3x2y_D2z_a;
    Double I_ESP_H3xyz_F3z_a = I_ESP_I3xy2z_D2z_a+ABZ*I_ESP_H3xyz_D2z_a;
    Double I_ESP_H3x2z_F3z_a = I_ESP_I3x3z_D2z_a+ABZ*I_ESP_H3x2z_D2z_a;
    Double I_ESP_H2x3y_F3z_a = I_ESP_I2x3yz_D2z_a+ABZ*I_ESP_H2x3y_D2z_a;
    Double I_ESP_H2x2yz_F3z_a = I_ESP_I2x2y2z_D2z_a+ABZ*I_ESP_H2x2yz_D2z_a;
    Double I_ESP_H2xy2z_F3z_a = I_ESP_I2xy3z_D2z_a+ABZ*I_ESP_H2xy2z_D2z_a;
    Double I_ESP_H2x3z_F3z_a = I_ESP_I2x4z_D2z_a+ABZ*I_ESP_H2x3z_D2z_a;
    Double I_ESP_Hx4y_F3z_a = I_ESP_Ix4yz_D2z_a+ABZ*I_ESP_Hx4y_D2z_a;
    Double I_ESP_Hx3yz_F3z_a = I_ESP_Ix3y2z_D2z_a+ABZ*I_ESP_Hx3yz_D2z_a;
    Double I_ESP_Hx2y2z_F3z_a = I_ESP_Ix2y3z_D2z_a+ABZ*I_ESP_Hx2y2z_D2z_a;
    Double I_ESP_Hxy3z_F3z_a = I_ESP_Ixy4z_D2z_a+ABZ*I_ESP_Hxy3z_D2z_a;
    Double I_ESP_Hx4z_F3z_a = I_ESP_Ix5z_D2z_a+ABZ*I_ESP_Hx4z_D2z_a;
    Double I_ESP_H5y_F3z_a = I_ESP_I5yz_D2z_a+ABZ*I_ESP_H5y_D2z_a;
    Double I_ESP_H4yz_F3z_a = I_ESP_I4y2z_D2z_a+ABZ*I_ESP_H4yz_D2z_a;
    Double I_ESP_H3y2z_F3z_a = I_ESP_I3y3z_D2z_a+ABZ*I_ESP_H3y2z_D2z_a;
    Double I_ESP_H2y3z_F3z_a = I_ESP_I2y4z_D2z_a+ABZ*I_ESP_H2y3z_D2z_a;
    Double I_ESP_Hy4z_F3z_a = I_ESP_Iy5z_D2z_a+ABZ*I_ESP_Hy4z_D2z_a;
    Double I_ESP_H5z_F3z_a = I_ESP_I6z_D2z_a+ABZ*I_ESP_H5z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_L_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_M_S_a
     * RHS shell quartet name: SQ_ESP_L_S_a
     ************************************************************/
    Double I_ESP_L8x_Px_a = I_ESP_M9x_S_a+ABX*I_ESP_L8x_S_a;
    Double I_ESP_L7xy_Px_a = I_ESP_M8xy_S_a+ABX*I_ESP_L7xy_S_a;
    Double I_ESP_L7xz_Px_a = I_ESP_M8xz_S_a+ABX*I_ESP_L7xz_S_a;
    Double I_ESP_L6x2y_Px_a = I_ESP_M7x2y_S_a+ABX*I_ESP_L6x2y_S_a;
    Double I_ESP_L6xyz_Px_a = I_ESP_M7xyz_S_a+ABX*I_ESP_L6xyz_S_a;
    Double I_ESP_L6x2z_Px_a = I_ESP_M7x2z_S_a+ABX*I_ESP_L6x2z_S_a;
    Double I_ESP_L5x3y_Px_a = I_ESP_M6x3y_S_a+ABX*I_ESP_L5x3y_S_a;
    Double I_ESP_L5x2yz_Px_a = I_ESP_M6x2yz_S_a+ABX*I_ESP_L5x2yz_S_a;
    Double I_ESP_L5xy2z_Px_a = I_ESP_M6xy2z_S_a+ABX*I_ESP_L5xy2z_S_a;
    Double I_ESP_L5x3z_Px_a = I_ESP_M6x3z_S_a+ABX*I_ESP_L5x3z_S_a;
    Double I_ESP_L4x4y_Px_a = I_ESP_M5x4y_S_a+ABX*I_ESP_L4x4y_S_a;
    Double I_ESP_L4x3yz_Px_a = I_ESP_M5x3yz_S_a+ABX*I_ESP_L4x3yz_S_a;
    Double I_ESP_L4x2y2z_Px_a = I_ESP_M5x2y2z_S_a+ABX*I_ESP_L4x2y2z_S_a;
    Double I_ESP_L4xy3z_Px_a = I_ESP_M5xy3z_S_a+ABX*I_ESP_L4xy3z_S_a;
    Double I_ESP_L4x4z_Px_a = I_ESP_M5x4z_S_a+ABX*I_ESP_L4x4z_S_a;
    Double I_ESP_L3x5y_Px_a = I_ESP_M4x5y_S_a+ABX*I_ESP_L3x5y_S_a;
    Double I_ESP_L3x4yz_Px_a = I_ESP_M4x4yz_S_a+ABX*I_ESP_L3x4yz_S_a;
    Double I_ESP_L3x3y2z_Px_a = I_ESP_M4x3y2z_S_a+ABX*I_ESP_L3x3y2z_S_a;
    Double I_ESP_L3x2y3z_Px_a = I_ESP_M4x2y3z_S_a+ABX*I_ESP_L3x2y3z_S_a;
    Double I_ESP_L3xy4z_Px_a = I_ESP_M4xy4z_S_a+ABX*I_ESP_L3xy4z_S_a;
    Double I_ESP_L3x5z_Px_a = I_ESP_M4x5z_S_a+ABX*I_ESP_L3x5z_S_a;
    Double I_ESP_L2x6y_Px_a = I_ESP_M3x6y_S_a+ABX*I_ESP_L2x6y_S_a;
    Double I_ESP_L2x5yz_Px_a = I_ESP_M3x5yz_S_a+ABX*I_ESP_L2x5yz_S_a;
    Double I_ESP_L2x4y2z_Px_a = I_ESP_M3x4y2z_S_a+ABX*I_ESP_L2x4y2z_S_a;
    Double I_ESP_L2x3y3z_Px_a = I_ESP_M3x3y3z_S_a+ABX*I_ESP_L2x3y3z_S_a;
    Double I_ESP_L2x2y4z_Px_a = I_ESP_M3x2y4z_S_a+ABX*I_ESP_L2x2y4z_S_a;
    Double I_ESP_L2xy5z_Px_a = I_ESP_M3xy5z_S_a+ABX*I_ESP_L2xy5z_S_a;
    Double I_ESP_L2x6z_Px_a = I_ESP_M3x6z_S_a+ABX*I_ESP_L2x6z_S_a;
    Double I_ESP_Lx7y_Px_a = I_ESP_M2x7y_S_a+ABX*I_ESP_Lx7y_S_a;
    Double I_ESP_Lx6yz_Px_a = I_ESP_M2x6yz_S_a+ABX*I_ESP_Lx6yz_S_a;
    Double I_ESP_Lx5y2z_Px_a = I_ESP_M2x5y2z_S_a+ABX*I_ESP_Lx5y2z_S_a;
    Double I_ESP_Lx4y3z_Px_a = I_ESP_M2x4y3z_S_a+ABX*I_ESP_Lx4y3z_S_a;
    Double I_ESP_Lx3y4z_Px_a = I_ESP_M2x3y4z_S_a+ABX*I_ESP_Lx3y4z_S_a;
    Double I_ESP_Lx2y5z_Px_a = I_ESP_M2x2y5z_S_a+ABX*I_ESP_Lx2y5z_S_a;
    Double I_ESP_Lxy6z_Px_a = I_ESP_M2xy6z_S_a+ABX*I_ESP_Lxy6z_S_a;
    Double I_ESP_Lx7z_Px_a = I_ESP_M2x7z_S_a+ABX*I_ESP_Lx7z_S_a;
    Double I_ESP_L6x2y_Py_a = I_ESP_M6x3y_S_a+ABY*I_ESP_L6x2y_S_a;
    Double I_ESP_L6xyz_Py_a = I_ESP_M6x2yz_S_a+ABY*I_ESP_L6xyz_S_a;
    Double I_ESP_L5x3y_Py_a = I_ESP_M5x4y_S_a+ABY*I_ESP_L5x3y_S_a;
    Double I_ESP_L5x2yz_Py_a = I_ESP_M5x3yz_S_a+ABY*I_ESP_L5x2yz_S_a;
    Double I_ESP_L5xy2z_Py_a = I_ESP_M5x2y2z_S_a+ABY*I_ESP_L5xy2z_S_a;
    Double I_ESP_L4x4y_Py_a = I_ESP_M4x5y_S_a+ABY*I_ESP_L4x4y_S_a;
    Double I_ESP_L4x3yz_Py_a = I_ESP_M4x4yz_S_a+ABY*I_ESP_L4x3yz_S_a;
    Double I_ESP_L4x2y2z_Py_a = I_ESP_M4x3y2z_S_a+ABY*I_ESP_L4x2y2z_S_a;
    Double I_ESP_L4xy3z_Py_a = I_ESP_M4x2y3z_S_a+ABY*I_ESP_L4xy3z_S_a;
    Double I_ESP_L3x5y_Py_a = I_ESP_M3x6y_S_a+ABY*I_ESP_L3x5y_S_a;
    Double I_ESP_L3x4yz_Py_a = I_ESP_M3x5yz_S_a+ABY*I_ESP_L3x4yz_S_a;
    Double I_ESP_L3x3y2z_Py_a = I_ESP_M3x4y2z_S_a+ABY*I_ESP_L3x3y2z_S_a;
    Double I_ESP_L3x2y3z_Py_a = I_ESP_M3x3y3z_S_a+ABY*I_ESP_L3x2y3z_S_a;
    Double I_ESP_L3xy4z_Py_a = I_ESP_M3x2y4z_S_a+ABY*I_ESP_L3xy4z_S_a;
    Double I_ESP_L2x6y_Py_a = I_ESP_M2x7y_S_a+ABY*I_ESP_L2x6y_S_a;
    Double I_ESP_L2x5yz_Py_a = I_ESP_M2x6yz_S_a+ABY*I_ESP_L2x5yz_S_a;
    Double I_ESP_L2x4y2z_Py_a = I_ESP_M2x5y2z_S_a+ABY*I_ESP_L2x4y2z_S_a;
    Double I_ESP_L2x3y3z_Py_a = I_ESP_M2x4y3z_S_a+ABY*I_ESP_L2x3y3z_S_a;
    Double I_ESP_L2x2y4z_Py_a = I_ESP_M2x3y4z_S_a+ABY*I_ESP_L2x2y4z_S_a;
    Double I_ESP_L2xy5z_Py_a = I_ESP_M2x2y5z_S_a+ABY*I_ESP_L2xy5z_S_a;
    Double I_ESP_Lx7y_Py_a = I_ESP_Mx8y_S_a+ABY*I_ESP_Lx7y_S_a;
    Double I_ESP_Lx6yz_Py_a = I_ESP_Mx7yz_S_a+ABY*I_ESP_Lx6yz_S_a;
    Double I_ESP_Lx5y2z_Py_a = I_ESP_Mx6y2z_S_a+ABY*I_ESP_Lx5y2z_S_a;
    Double I_ESP_Lx4y3z_Py_a = I_ESP_Mx5y3z_S_a+ABY*I_ESP_Lx4y3z_S_a;
    Double I_ESP_Lx3y4z_Py_a = I_ESP_Mx4y4z_S_a+ABY*I_ESP_Lx3y4z_S_a;
    Double I_ESP_Lx2y5z_Py_a = I_ESP_Mx3y5z_S_a+ABY*I_ESP_Lx2y5z_S_a;
    Double I_ESP_Lxy6z_Py_a = I_ESP_Mx2y6z_S_a+ABY*I_ESP_Lxy6z_S_a;
    Double I_ESP_L8y_Py_a = I_ESP_M9y_S_a+ABY*I_ESP_L8y_S_a;
    Double I_ESP_L7yz_Py_a = I_ESP_M8yz_S_a+ABY*I_ESP_L7yz_S_a;
    Double I_ESP_L6y2z_Py_a = I_ESP_M7y2z_S_a+ABY*I_ESP_L6y2z_S_a;
    Double I_ESP_L5y3z_Py_a = I_ESP_M6y3z_S_a+ABY*I_ESP_L5y3z_S_a;
    Double I_ESP_L4y4z_Py_a = I_ESP_M5y4z_S_a+ABY*I_ESP_L4y4z_S_a;
    Double I_ESP_L3y5z_Py_a = I_ESP_M4y5z_S_a+ABY*I_ESP_L3y5z_S_a;
    Double I_ESP_L2y6z_Py_a = I_ESP_M3y6z_S_a+ABY*I_ESP_L2y6z_S_a;
    Double I_ESP_Ly7z_Py_a = I_ESP_M2y7z_S_a+ABY*I_ESP_Ly7z_S_a;
    Double I_ESP_L6xyz_Pz_a = I_ESP_M6xy2z_S_a+ABZ*I_ESP_L6xyz_S_a;
    Double I_ESP_L6x2z_Pz_a = I_ESP_M6x3z_S_a+ABZ*I_ESP_L6x2z_S_a;
    Double I_ESP_L5x2yz_Pz_a = I_ESP_M5x2y2z_S_a+ABZ*I_ESP_L5x2yz_S_a;
    Double I_ESP_L5xy2z_Pz_a = I_ESP_M5xy3z_S_a+ABZ*I_ESP_L5xy2z_S_a;
    Double I_ESP_L5x3z_Pz_a = I_ESP_M5x4z_S_a+ABZ*I_ESP_L5x3z_S_a;
    Double I_ESP_L4x3yz_Pz_a = I_ESP_M4x3y2z_S_a+ABZ*I_ESP_L4x3yz_S_a;
    Double I_ESP_L4x2y2z_Pz_a = I_ESP_M4x2y3z_S_a+ABZ*I_ESP_L4x2y2z_S_a;
    Double I_ESP_L4xy3z_Pz_a = I_ESP_M4xy4z_S_a+ABZ*I_ESP_L4xy3z_S_a;
    Double I_ESP_L4x4z_Pz_a = I_ESP_M4x5z_S_a+ABZ*I_ESP_L4x4z_S_a;
    Double I_ESP_L3x4yz_Pz_a = I_ESP_M3x4y2z_S_a+ABZ*I_ESP_L3x4yz_S_a;
    Double I_ESP_L3x3y2z_Pz_a = I_ESP_M3x3y3z_S_a+ABZ*I_ESP_L3x3y2z_S_a;
    Double I_ESP_L3x2y3z_Pz_a = I_ESP_M3x2y4z_S_a+ABZ*I_ESP_L3x2y3z_S_a;
    Double I_ESP_L3xy4z_Pz_a = I_ESP_M3xy5z_S_a+ABZ*I_ESP_L3xy4z_S_a;
    Double I_ESP_L3x5z_Pz_a = I_ESP_M3x6z_S_a+ABZ*I_ESP_L3x5z_S_a;
    Double I_ESP_L2x5yz_Pz_a = I_ESP_M2x5y2z_S_a+ABZ*I_ESP_L2x5yz_S_a;
    Double I_ESP_L2x4y2z_Pz_a = I_ESP_M2x4y3z_S_a+ABZ*I_ESP_L2x4y2z_S_a;
    Double I_ESP_L2x3y3z_Pz_a = I_ESP_M2x3y4z_S_a+ABZ*I_ESP_L2x3y3z_S_a;
    Double I_ESP_L2x2y4z_Pz_a = I_ESP_M2x2y5z_S_a+ABZ*I_ESP_L2x2y4z_S_a;
    Double I_ESP_L2xy5z_Pz_a = I_ESP_M2xy6z_S_a+ABZ*I_ESP_L2xy5z_S_a;
    Double I_ESP_L2x6z_Pz_a = I_ESP_M2x7z_S_a+ABZ*I_ESP_L2x6z_S_a;
    Double I_ESP_Lx6yz_Pz_a = I_ESP_Mx6y2z_S_a+ABZ*I_ESP_Lx6yz_S_a;
    Double I_ESP_Lx5y2z_Pz_a = I_ESP_Mx5y3z_S_a+ABZ*I_ESP_Lx5y2z_S_a;
    Double I_ESP_Lx4y3z_Pz_a = I_ESP_Mx4y4z_S_a+ABZ*I_ESP_Lx4y3z_S_a;
    Double I_ESP_Lx3y4z_Pz_a = I_ESP_Mx3y5z_S_a+ABZ*I_ESP_Lx3y4z_S_a;
    Double I_ESP_Lx2y5z_Pz_a = I_ESP_Mx2y6z_S_a+ABZ*I_ESP_Lx2y5z_S_a;
    Double I_ESP_Lxy6z_Pz_a = I_ESP_Mxy7z_S_a+ABZ*I_ESP_Lxy6z_S_a;
    Double I_ESP_Lx7z_Pz_a = I_ESP_Mx8z_S_a+ABZ*I_ESP_Lx7z_S_a;
    Double I_ESP_L6y2z_Pz_a = I_ESP_M6y3z_S_a+ABZ*I_ESP_L6y2z_S_a;
    Double I_ESP_L5y3z_Pz_a = I_ESP_M5y4z_S_a+ABZ*I_ESP_L5y3z_S_a;
    Double I_ESP_L4y4z_Pz_a = I_ESP_M4y5z_S_a+ABZ*I_ESP_L4y4z_S_a;
    Double I_ESP_L3y5z_Pz_a = I_ESP_M3y6z_S_a+ABZ*I_ESP_L3y5z_S_a;
    Double I_ESP_L2y6z_Pz_a = I_ESP_M2y7z_S_a+ABZ*I_ESP_L2y6z_S_a;
    Double I_ESP_Ly7z_Pz_a = I_ESP_My8z_S_a+ABZ*I_ESP_Ly7z_S_a;
    Double I_ESP_L8z_Pz_a = I_ESP_M9z_S_a+ABZ*I_ESP_L8z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_K_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 111 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_P_a
     * RHS shell quartet name: SQ_ESP_K_P_a
     ************************************************************/
    Double I_ESP_K7x_D2x_a = I_ESP_L8x_Px_a+ABX*I_ESP_K7x_Px_a;
    Double I_ESP_K6xy_D2x_a = I_ESP_L7xy_Px_a+ABX*I_ESP_K6xy_Px_a;
    Double I_ESP_K6xz_D2x_a = I_ESP_L7xz_Px_a+ABX*I_ESP_K6xz_Px_a;
    Double I_ESP_K5x2y_D2x_a = I_ESP_L6x2y_Px_a+ABX*I_ESP_K5x2y_Px_a;
    Double I_ESP_K5xyz_D2x_a = I_ESP_L6xyz_Px_a+ABX*I_ESP_K5xyz_Px_a;
    Double I_ESP_K5x2z_D2x_a = I_ESP_L6x2z_Px_a+ABX*I_ESP_K5x2z_Px_a;
    Double I_ESP_K4x3y_D2x_a = I_ESP_L5x3y_Px_a+ABX*I_ESP_K4x3y_Px_a;
    Double I_ESP_K4x2yz_D2x_a = I_ESP_L5x2yz_Px_a+ABX*I_ESP_K4x2yz_Px_a;
    Double I_ESP_K4xy2z_D2x_a = I_ESP_L5xy2z_Px_a+ABX*I_ESP_K4xy2z_Px_a;
    Double I_ESP_K4x3z_D2x_a = I_ESP_L5x3z_Px_a+ABX*I_ESP_K4x3z_Px_a;
    Double I_ESP_K3x4y_D2x_a = I_ESP_L4x4y_Px_a+ABX*I_ESP_K3x4y_Px_a;
    Double I_ESP_K3x3yz_D2x_a = I_ESP_L4x3yz_Px_a+ABX*I_ESP_K3x3yz_Px_a;
    Double I_ESP_K3x2y2z_D2x_a = I_ESP_L4x2y2z_Px_a+ABX*I_ESP_K3x2y2z_Px_a;
    Double I_ESP_K3xy3z_D2x_a = I_ESP_L4xy3z_Px_a+ABX*I_ESP_K3xy3z_Px_a;
    Double I_ESP_K3x4z_D2x_a = I_ESP_L4x4z_Px_a+ABX*I_ESP_K3x4z_Px_a;
    Double I_ESP_K2x5y_D2x_a = I_ESP_L3x5y_Px_a+ABX*I_ESP_K2x5y_Px_a;
    Double I_ESP_K2x4yz_D2x_a = I_ESP_L3x4yz_Px_a+ABX*I_ESP_K2x4yz_Px_a;
    Double I_ESP_K2x3y2z_D2x_a = I_ESP_L3x3y2z_Px_a+ABX*I_ESP_K2x3y2z_Px_a;
    Double I_ESP_K2x2y3z_D2x_a = I_ESP_L3x2y3z_Px_a+ABX*I_ESP_K2x2y3z_Px_a;
    Double I_ESP_K2xy4z_D2x_a = I_ESP_L3xy4z_Px_a+ABX*I_ESP_K2xy4z_Px_a;
    Double I_ESP_K2x5z_D2x_a = I_ESP_L3x5z_Px_a+ABX*I_ESP_K2x5z_Px_a;
    Double I_ESP_Kx6y_D2x_a = I_ESP_L2x6y_Px_a+ABX*I_ESP_Kx6y_Px_a;
    Double I_ESP_Kx5yz_D2x_a = I_ESP_L2x5yz_Px_a+ABX*I_ESP_Kx5yz_Px_a;
    Double I_ESP_Kx4y2z_D2x_a = I_ESP_L2x4y2z_Px_a+ABX*I_ESP_Kx4y2z_Px_a;
    Double I_ESP_Kx3y3z_D2x_a = I_ESP_L2x3y3z_Px_a+ABX*I_ESP_Kx3y3z_Px_a;
    Double I_ESP_Kx2y4z_D2x_a = I_ESP_L2x2y4z_Px_a+ABX*I_ESP_Kx2y4z_Px_a;
    Double I_ESP_Kxy5z_D2x_a = I_ESP_L2xy5z_Px_a+ABX*I_ESP_Kxy5z_Px_a;
    Double I_ESP_Kx6z_D2x_a = I_ESP_L2x6z_Px_a+ABX*I_ESP_Kx6z_Px_a;
    Double I_ESP_K7y_D2x_a = I_ESP_Lx7y_Px_a+ABX*I_ESP_K7y_Px_a;
    Double I_ESP_K6yz_D2x_a = I_ESP_Lx6yz_Px_a+ABX*I_ESP_K6yz_Px_a;
    Double I_ESP_K5y2z_D2x_a = I_ESP_Lx5y2z_Px_a+ABX*I_ESP_K5y2z_Px_a;
    Double I_ESP_K4y3z_D2x_a = I_ESP_Lx4y3z_Px_a+ABX*I_ESP_K4y3z_Px_a;
    Double I_ESP_K3y4z_D2x_a = I_ESP_Lx3y4z_Px_a+ABX*I_ESP_K3y4z_Px_a;
    Double I_ESP_K2y5z_D2x_a = I_ESP_Lx2y5z_Px_a+ABX*I_ESP_K2y5z_Px_a;
    Double I_ESP_Ky6z_D2x_a = I_ESP_Lxy6z_Px_a+ABX*I_ESP_Ky6z_Px_a;
    Double I_ESP_K7z_D2x_a = I_ESP_Lx7z_Px_a+ABX*I_ESP_K7z_Px_a;
    Double I_ESP_K6xy_D2y_a = I_ESP_L6x2y_Py_a+ABY*I_ESP_K6xy_Py_a;
    Double I_ESP_K6xz_D2y_a = I_ESP_L6xyz_Py_a+ABY*I_ESP_K6xz_Py_a;
    Double I_ESP_K5x2y_D2y_a = I_ESP_L5x3y_Py_a+ABY*I_ESP_K5x2y_Py_a;
    Double I_ESP_K5xyz_D2y_a = I_ESP_L5x2yz_Py_a+ABY*I_ESP_K5xyz_Py_a;
    Double I_ESP_K5x2z_D2y_a = I_ESP_L5xy2z_Py_a+ABY*I_ESP_K5x2z_Py_a;
    Double I_ESP_K4x3y_D2y_a = I_ESP_L4x4y_Py_a+ABY*I_ESP_K4x3y_Py_a;
    Double I_ESP_K4x2yz_D2y_a = I_ESP_L4x3yz_Py_a+ABY*I_ESP_K4x2yz_Py_a;
    Double I_ESP_K4xy2z_D2y_a = I_ESP_L4x2y2z_Py_a+ABY*I_ESP_K4xy2z_Py_a;
    Double I_ESP_K4x3z_D2y_a = I_ESP_L4xy3z_Py_a+ABY*I_ESP_K4x3z_Py_a;
    Double I_ESP_K3x4y_D2y_a = I_ESP_L3x5y_Py_a+ABY*I_ESP_K3x4y_Py_a;
    Double I_ESP_K3x3yz_D2y_a = I_ESP_L3x4yz_Py_a+ABY*I_ESP_K3x3yz_Py_a;
    Double I_ESP_K3x2y2z_D2y_a = I_ESP_L3x3y2z_Py_a+ABY*I_ESP_K3x2y2z_Py_a;
    Double I_ESP_K3xy3z_D2y_a = I_ESP_L3x2y3z_Py_a+ABY*I_ESP_K3xy3z_Py_a;
    Double I_ESP_K3x4z_D2y_a = I_ESP_L3xy4z_Py_a+ABY*I_ESP_K3x4z_Py_a;
    Double I_ESP_K2x5y_D2y_a = I_ESP_L2x6y_Py_a+ABY*I_ESP_K2x5y_Py_a;
    Double I_ESP_K2x4yz_D2y_a = I_ESP_L2x5yz_Py_a+ABY*I_ESP_K2x4yz_Py_a;
    Double I_ESP_K2x3y2z_D2y_a = I_ESP_L2x4y2z_Py_a+ABY*I_ESP_K2x3y2z_Py_a;
    Double I_ESP_K2x2y3z_D2y_a = I_ESP_L2x3y3z_Py_a+ABY*I_ESP_K2x2y3z_Py_a;
    Double I_ESP_K2xy4z_D2y_a = I_ESP_L2x2y4z_Py_a+ABY*I_ESP_K2xy4z_Py_a;
    Double I_ESP_K2x5z_D2y_a = I_ESP_L2xy5z_Py_a+ABY*I_ESP_K2x5z_Py_a;
    Double I_ESP_Kx6y_D2y_a = I_ESP_Lx7y_Py_a+ABY*I_ESP_Kx6y_Py_a;
    Double I_ESP_Kx5yz_D2y_a = I_ESP_Lx6yz_Py_a+ABY*I_ESP_Kx5yz_Py_a;
    Double I_ESP_Kx4y2z_D2y_a = I_ESP_Lx5y2z_Py_a+ABY*I_ESP_Kx4y2z_Py_a;
    Double I_ESP_Kx3y3z_D2y_a = I_ESP_Lx4y3z_Py_a+ABY*I_ESP_Kx3y3z_Py_a;
    Double I_ESP_Kx2y4z_D2y_a = I_ESP_Lx3y4z_Py_a+ABY*I_ESP_Kx2y4z_Py_a;
    Double I_ESP_Kxy5z_D2y_a = I_ESP_Lx2y5z_Py_a+ABY*I_ESP_Kxy5z_Py_a;
    Double I_ESP_Kx6z_D2y_a = I_ESP_Lxy6z_Py_a+ABY*I_ESP_Kx6z_Py_a;
    Double I_ESP_K7y_D2y_a = I_ESP_L8y_Py_a+ABY*I_ESP_K7y_Py_a;
    Double I_ESP_K6yz_D2y_a = I_ESP_L7yz_Py_a+ABY*I_ESP_K6yz_Py_a;
    Double I_ESP_K5y2z_D2y_a = I_ESP_L6y2z_Py_a+ABY*I_ESP_K5y2z_Py_a;
    Double I_ESP_K4y3z_D2y_a = I_ESP_L5y3z_Py_a+ABY*I_ESP_K4y3z_Py_a;
    Double I_ESP_K3y4z_D2y_a = I_ESP_L4y4z_Py_a+ABY*I_ESP_K3y4z_Py_a;
    Double I_ESP_K2y5z_D2y_a = I_ESP_L3y5z_Py_a+ABY*I_ESP_K2y5z_Py_a;
    Double I_ESP_Ky6z_D2y_a = I_ESP_L2y6z_Py_a+ABY*I_ESP_Ky6z_Py_a;
    Double I_ESP_K7z_D2y_a = I_ESP_Ly7z_Py_a+ABY*I_ESP_K7z_Py_a;
    Double I_ESP_K6xy_D2z_a = I_ESP_L6xyz_Pz_a+ABZ*I_ESP_K6xy_Pz_a;
    Double I_ESP_K6xz_D2z_a = I_ESP_L6x2z_Pz_a+ABZ*I_ESP_K6xz_Pz_a;
    Double I_ESP_K5x2y_D2z_a = I_ESP_L5x2yz_Pz_a+ABZ*I_ESP_K5x2y_Pz_a;
    Double I_ESP_K5xyz_D2z_a = I_ESP_L5xy2z_Pz_a+ABZ*I_ESP_K5xyz_Pz_a;
    Double I_ESP_K5x2z_D2z_a = I_ESP_L5x3z_Pz_a+ABZ*I_ESP_K5x2z_Pz_a;
    Double I_ESP_K4x3y_D2z_a = I_ESP_L4x3yz_Pz_a+ABZ*I_ESP_K4x3y_Pz_a;
    Double I_ESP_K4x2yz_D2z_a = I_ESP_L4x2y2z_Pz_a+ABZ*I_ESP_K4x2yz_Pz_a;
    Double I_ESP_K4xy2z_D2z_a = I_ESP_L4xy3z_Pz_a+ABZ*I_ESP_K4xy2z_Pz_a;
    Double I_ESP_K4x3z_D2z_a = I_ESP_L4x4z_Pz_a+ABZ*I_ESP_K4x3z_Pz_a;
    Double I_ESP_K3x4y_D2z_a = I_ESP_L3x4yz_Pz_a+ABZ*I_ESP_K3x4y_Pz_a;
    Double I_ESP_K3x3yz_D2z_a = I_ESP_L3x3y2z_Pz_a+ABZ*I_ESP_K3x3yz_Pz_a;
    Double I_ESP_K3x2y2z_D2z_a = I_ESP_L3x2y3z_Pz_a+ABZ*I_ESP_K3x2y2z_Pz_a;
    Double I_ESP_K3xy3z_D2z_a = I_ESP_L3xy4z_Pz_a+ABZ*I_ESP_K3xy3z_Pz_a;
    Double I_ESP_K3x4z_D2z_a = I_ESP_L3x5z_Pz_a+ABZ*I_ESP_K3x4z_Pz_a;
    Double I_ESP_K2x5y_D2z_a = I_ESP_L2x5yz_Pz_a+ABZ*I_ESP_K2x5y_Pz_a;
    Double I_ESP_K2x4yz_D2z_a = I_ESP_L2x4y2z_Pz_a+ABZ*I_ESP_K2x4yz_Pz_a;
    Double I_ESP_K2x3y2z_D2z_a = I_ESP_L2x3y3z_Pz_a+ABZ*I_ESP_K2x3y2z_Pz_a;
    Double I_ESP_K2x2y3z_D2z_a = I_ESP_L2x2y4z_Pz_a+ABZ*I_ESP_K2x2y3z_Pz_a;
    Double I_ESP_K2xy4z_D2z_a = I_ESP_L2xy5z_Pz_a+ABZ*I_ESP_K2xy4z_Pz_a;
    Double I_ESP_K2x5z_D2z_a = I_ESP_L2x6z_Pz_a+ABZ*I_ESP_K2x5z_Pz_a;
    Double I_ESP_Kx6y_D2z_a = I_ESP_Lx6yz_Pz_a+ABZ*I_ESP_Kx6y_Pz_a;
    Double I_ESP_Kx5yz_D2z_a = I_ESP_Lx5y2z_Pz_a+ABZ*I_ESP_Kx5yz_Pz_a;
    Double I_ESP_Kx4y2z_D2z_a = I_ESP_Lx4y3z_Pz_a+ABZ*I_ESP_Kx4y2z_Pz_a;
    Double I_ESP_Kx3y3z_D2z_a = I_ESP_Lx3y4z_Pz_a+ABZ*I_ESP_Kx3y3z_Pz_a;
    Double I_ESP_Kx2y4z_D2z_a = I_ESP_Lx2y5z_Pz_a+ABZ*I_ESP_Kx2y4z_Pz_a;
    Double I_ESP_Kxy5z_D2z_a = I_ESP_Lxy6z_Pz_a+ABZ*I_ESP_Kxy5z_Pz_a;
    Double I_ESP_Kx6z_D2z_a = I_ESP_Lx7z_Pz_a+ABZ*I_ESP_Kx6z_Pz_a;
    Double I_ESP_K6yz_D2z_a = I_ESP_L6y2z_Pz_a+ABZ*I_ESP_K6yz_Pz_a;
    Double I_ESP_K5y2z_D2z_a = I_ESP_L5y3z_Pz_a+ABZ*I_ESP_K5y2z_Pz_a;
    Double I_ESP_K4y3z_D2z_a = I_ESP_L4y4z_Pz_a+ABZ*I_ESP_K4y3z_Pz_a;
    Double I_ESP_K3y4z_D2z_a = I_ESP_L3y5z_Pz_a+ABZ*I_ESP_K3y4z_Pz_a;
    Double I_ESP_K2y5z_D2z_a = I_ESP_L2y6z_Pz_a+ABZ*I_ESP_K2y5z_Pz_a;
    Double I_ESP_Ky6z_D2z_a = I_ESP_Ly7z_Pz_a+ABZ*I_ESP_Ky6z_Pz_a;
    Double I_ESP_K7z_D2z_a = I_ESP_L8z_Pz_a+ABZ*I_ESP_K7z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_F_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 85 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D_a
     * RHS shell quartet name: SQ_ESP_I_D_a
     ************************************************************/
    Double I_ESP_I6x_F3x_a = I_ESP_K7x_D2x_a+ABX*I_ESP_I6x_D2x_a;
    Double I_ESP_I5xy_F3x_a = I_ESP_K6xy_D2x_a+ABX*I_ESP_I5xy_D2x_a;
    Double I_ESP_I5xz_F3x_a = I_ESP_K6xz_D2x_a+ABX*I_ESP_I5xz_D2x_a;
    Double I_ESP_I4x2y_F3x_a = I_ESP_K5x2y_D2x_a+ABX*I_ESP_I4x2y_D2x_a;
    Double I_ESP_I4xyz_F3x_a = I_ESP_K5xyz_D2x_a+ABX*I_ESP_I4xyz_D2x_a;
    Double I_ESP_I4x2z_F3x_a = I_ESP_K5x2z_D2x_a+ABX*I_ESP_I4x2z_D2x_a;
    Double I_ESP_I3x3y_F3x_a = I_ESP_K4x3y_D2x_a+ABX*I_ESP_I3x3y_D2x_a;
    Double I_ESP_I3x2yz_F3x_a = I_ESP_K4x2yz_D2x_a+ABX*I_ESP_I3x2yz_D2x_a;
    Double I_ESP_I3xy2z_F3x_a = I_ESP_K4xy2z_D2x_a+ABX*I_ESP_I3xy2z_D2x_a;
    Double I_ESP_I3x3z_F3x_a = I_ESP_K4x3z_D2x_a+ABX*I_ESP_I3x3z_D2x_a;
    Double I_ESP_I2x4y_F3x_a = I_ESP_K3x4y_D2x_a+ABX*I_ESP_I2x4y_D2x_a;
    Double I_ESP_I2x3yz_F3x_a = I_ESP_K3x3yz_D2x_a+ABX*I_ESP_I2x3yz_D2x_a;
    Double I_ESP_I2x2y2z_F3x_a = I_ESP_K3x2y2z_D2x_a+ABX*I_ESP_I2x2y2z_D2x_a;
    Double I_ESP_I2xy3z_F3x_a = I_ESP_K3xy3z_D2x_a+ABX*I_ESP_I2xy3z_D2x_a;
    Double I_ESP_I2x4z_F3x_a = I_ESP_K3x4z_D2x_a+ABX*I_ESP_I2x4z_D2x_a;
    Double I_ESP_Ix5y_F3x_a = I_ESP_K2x5y_D2x_a+ABX*I_ESP_Ix5y_D2x_a;
    Double I_ESP_Ix4yz_F3x_a = I_ESP_K2x4yz_D2x_a+ABX*I_ESP_Ix4yz_D2x_a;
    Double I_ESP_Ix3y2z_F3x_a = I_ESP_K2x3y2z_D2x_a+ABX*I_ESP_Ix3y2z_D2x_a;
    Double I_ESP_Ix2y3z_F3x_a = I_ESP_K2x2y3z_D2x_a+ABX*I_ESP_Ix2y3z_D2x_a;
    Double I_ESP_Ixy4z_F3x_a = I_ESP_K2xy4z_D2x_a+ABX*I_ESP_Ixy4z_D2x_a;
    Double I_ESP_Ix5z_F3x_a = I_ESP_K2x5z_D2x_a+ABX*I_ESP_Ix5z_D2x_a;
    Double I_ESP_I6y_F3x_a = I_ESP_Kx6y_D2x_a+ABX*I_ESP_I6y_D2x_a;
    Double I_ESP_I5yz_F3x_a = I_ESP_Kx5yz_D2x_a+ABX*I_ESP_I5yz_D2x_a;
    Double I_ESP_I4y2z_F3x_a = I_ESP_Kx4y2z_D2x_a+ABX*I_ESP_I4y2z_D2x_a;
    Double I_ESP_I3y3z_F3x_a = I_ESP_Kx3y3z_D2x_a+ABX*I_ESP_I3y3z_D2x_a;
    Double I_ESP_I2y4z_F3x_a = I_ESP_Kx2y4z_D2x_a+ABX*I_ESP_I2y4z_D2x_a;
    Double I_ESP_Iy5z_F3x_a = I_ESP_Kxy5z_D2x_a+ABX*I_ESP_Iy5z_D2x_a;
    Double I_ESP_I6z_F3x_a = I_ESP_Kx6z_D2x_a+ABX*I_ESP_I6z_D2x_a;
    Double I_ESP_I5xy_F2xy_a = I_ESP_K5x2y_D2x_a+ABY*I_ESP_I5xy_D2x_a;
    Double I_ESP_I5xz_F2xy_a = I_ESP_K5xyz_D2x_a+ABY*I_ESP_I5xz_D2x_a;
    Double I_ESP_I4x2y_F2xy_a = I_ESP_K4x3y_D2x_a+ABY*I_ESP_I4x2y_D2x_a;
    Double I_ESP_I4xyz_F2xy_a = I_ESP_K4x2yz_D2x_a+ABY*I_ESP_I4xyz_D2x_a;
    Double I_ESP_I4x2z_F2xy_a = I_ESP_K4xy2z_D2x_a+ABY*I_ESP_I4x2z_D2x_a;
    Double I_ESP_I3x3y_F2xy_a = I_ESP_K3x4y_D2x_a+ABY*I_ESP_I3x3y_D2x_a;
    Double I_ESP_I3x2yz_F2xy_a = I_ESP_K3x3yz_D2x_a+ABY*I_ESP_I3x2yz_D2x_a;
    Double I_ESP_I3xy2z_F2xy_a = I_ESP_K3x2y2z_D2x_a+ABY*I_ESP_I3xy2z_D2x_a;
    Double I_ESP_I3x3z_F2xy_a = I_ESP_K3xy3z_D2x_a+ABY*I_ESP_I3x3z_D2x_a;
    Double I_ESP_I2x4y_F2xy_a = I_ESP_K2x5y_D2x_a+ABY*I_ESP_I2x4y_D2x_a;
    Double I_ESP_I2x3yz_F2xy_a = I_ESP_K2x4yz_D2x_a+ABY*I_ESP_I2x3yz_D2x_a;
    Double I_ESP_I2x2y2z_F2xy_a = I_ESP_K2x3y2z_D2x_a+ABY*I_ESP_I2x2y2z_D2x_a;
    Double I_ESP_I2xy3z_F2xy_a = I_ESP_K2x2y3z_D2x_a+ABY*I_ESP_I2xy3z_D2x_a;
    Double I_ESP_I2x4z_F2xy_a = I_ESP_K2xy4z_D2x_a+ABY*I_ESP_I2x4z_D2x_a;
    Double I_ESP_Ix5y_F2xy_a = I_ESP_Kx6y_D2x_a+ABY*I_ESP_Ix5y_D2x_a;
    Double I_ESP_Ix4yz_F2xy_a = I_ESP_Kx5yz_D2x_a+ABY*I_ESP_Ix4yz_D2x_a;
    Double I_ESP_Ix3y2z_F2xy_a = I_ESP_Kx4y2z_D2x_a+ABY*I_ESP_Ix3y2z_D2x_a;
    Double I_ESP_Ix2y3z_F2xy_a = I_ESP_Kx3y3z_D2x_a+ABY*I_ESP_Ix2y3z_D2x_a;
    Double I_ESP_Ixy4z_F2xy_a = I_ESP_Kx2y4z_D2x_a+ABY*I_ESP_Ixy4z_D2x_a;
    Double I_ESP_Ix5z_F2xy_a = I_ESP_Kxy5z_D2x_a+ABY*I_ESP_Ix5z_D2x_a;
    Double I_ESP_I6y_F2xy_a = I_ESP_K7y_D2x_a+ABY*I_ESP_I6y_D2x_a;
    Double I_ESP_I5yz_F2xy_a = I_ESP_K6yz_D2x_a+ABY*I_ESP_I5yz_D2x_a;
    Double I_ESP_I4y2z_F2xy_a = I_ESP_K5y2z_D2x_a+ABY*I_ESP_I4y2z_D2x_a;
    Double I_ESP_I3y3z_F2xy_a = I_ESP_K4y3z_D2x_a+ABY*I_ESP_I3y3z_D2x_a;
    Double I_ESP_I2y4z_F2xy_a = I_ESP_K3y4z_D2x_a+ABY*I_ESP_I2y4z_D2x_a;
    Double I_ESP_Iy5z_F2xy_a = I_ESP_K2y5z_D2x_a+ABY*I_ESP_Iy5z_D2x_a;
    Double I_ESP_I6z_F2xy_a = I_ESP_Ky6z_D2x_a+ABY*I_ESP_I6z_D2x_a;
    Double I_ESP_I5xz_F2xz_a = I_ESP_K5x2z_D2x_a+ABZ*I_ESP_I5xz_D2x_a;
    Double I_ESP_I4xyz_F2xz_a = I_ESP_K4xy2z_D2x_a+ABZ*I_ESP_I4xyz_D2x_a;
    Double I_ESP_I4x2z_F2xz_a = I_ESP_K4x3z_D2x_a+ABZ*I_ESP_I4x2z_D2x_a;
    Double I_ESP_I3x2yz_F2xz_a = I_ESP_K3x2y2z_D2x_a+ABZ*I_ESP_I3x2yz_D2x_a;
    Double I_ESP_I3xy2z_F2xz_a = I_ESP_K3xy3z_D2x_a+ABZ*I_ESP_I3xy2z_D2x_a;
    Double I_ESP_I3x3z_F2xz_a = I_ESP_K3x4z_D2x_a+ABZ*I_ESP_I3x3z_D2x_a;
    Double I_ESP_I2x3yz_F2xz_a = I_ESP_K2x3y2z_D2x_a+ABZ*I_ESP_I2x3yz_D2x_a;
    Double I_ESP_I2x2y2z_F2xz_a = I_ESP_K2x2y3z_D2x_a+ABZ*I_ESP_I2x2y2z_D2x_a;
    Double I_ESP_I2xy3z_F2xz_a = I_ESP_K2xy4z_D2x_a+ABZ*I_ESP_I2xy3z_D2x_a;
    Double I_ESP_I2x4z_F2xz_a = I_ESP_K2x5z_D2x_a+ABZ*I_ESP_I2x4z_D2x_a;
    Double I_ESP_Ix4yz_F2xz_a = I_ESP_Kx4y2z_D2x_a+ABZ*I_ESP_Ix4yz_D2x_a;
    Double I_ESP_Ix3y2z_F2xz_a = I_ESP_Kx3y3z_D2x_a+ABZ*I_ESP_Ix3y2z_D2x_a;
    Double I_ESP_Ix2y3z_F2xz_a = I_ESP_Kx2y4z_D2x_a+ABZ*I_ESP_Ix2y3z_D2x_a;
    Double I_ESP_Ixy4z_F2xz_a = I_ESP_Kxy5z_D2x_a+ABZ*I_ESP_Ixy4z_D2x_a;
    Double I_ESP_Ix5z_F2xz_a = I_ESP_Kx6z_D2x_a+ABZ*I_ESP_Ix5z_D2x_a;
    Double I_ESP_I5yz_F2xz_a = I_ESP_K5y2z_D2x_a+ABZ*I_ESP_I5yz_D2x_a;
    Double I_ESP_I4y2z_F2xz_a = I_ESP_K4y3z_D2x_a+ABZ*I_ESP_I4y2z_D2x_a;
    Double I_ESP_I3y3z_F2xz_a = I_ESP_K3y4z_D2x_a+ABZ*I_ESP_I3y3z_D2x_a;
    Double I_ESP_I2y4z_F2xz_a = I_ESP_K2y5z_D2x_a+ABZ*I_ESP_I2y4z_D2x_a;
    Double I_ESP_Iy5z_F2xz_a = I_ESP_Ky6z_D2x_a+ABZ*I_ESP_Iy5z_D2x_a;
    Double I_ESP_I6z_F2xz_a = I_ESP_K7z_D2x_a+ABZ*I_ESP_I6z_D2x_a;
    Double I_ESP_I5xz_Fx2y_a = I_ESP_K6xz_D2y_a+ABX*I_ESP_I5xz_D2y_a;
    Double I_ESP_I4xyz_Fx2y_a = I_ESP_K5xyz_D2y_a+ABX*I_ESP_I4xyz_D2y_a;
    Double I_ESP_I4x2z_Fx2y_a = I_ESP_K5x2z_D2y_a+ABX*I_ESP_I4x2z_D2y_a;
    Double I_ESP_I3x2yz_Fx2y_a = I_ESP_K4x2yz_D2y_a+ABX*I_ESP_I3x2yz_D2y_a;
    Double I_ESP_I3xy2z_Fx2y_a = I_ESP_K4xy2z_D2y_a+ABX*I_ESP_I3xy2z_D2y_a;
    Double I_ESP_I3x3z_Fx2y_a = I_ESP_K4x3z_D2y_a+ABX*I_ESP_I3x3z_D2y_a;
    Double I_ESP_I2x3yz_Fx2y_a = I_ESP_K3x3yz_D2y_a+ABX*I_ESP_I2x3yz_D2y_a;
    Double I_ESP_I2x2y2z_Fx2y_a = I_ESP_K3x2y2z_D2y_a+ABX*I_ESP_I2x2y2z_D2y_a;
    Double I_ESP_I2xy3z_Fx2y_a = I_ESP_K3xy3z_D2y_a+ABX*I_ESP_I2xy3z_D2y_a;
    Double I_ESP_I2x4z_Fx2y_a = I_ESP_K3x4z_D2y_a+ABX*I_ESP_I2x4z_D2y_a;
    Double I_ESP_Ix4yz_Fx2y_a = I_ESP_K2x4yz_D2y_a+ABX*I_ESP_Ix4yz_D2y_a;
    Double I_ESP_Ix3y2z_Fx2y_a = I_ESP_K2x3y2z_D2y_a+ABX*I_ESP_Ix3y2z_D2y_a;
    Double I_ESP_Ix2y3z_Fx2y_a = I_ESP_K2x2y3z_D2y_a+ABX*I_ESP_Ix2y3z_D2y_a;
    Double I_ESP_Ixy4z_Fx2y_a = I_ESP_K2xy4z_D2y_a+ABX*I_ESP_Ixy4z_D2y_a;
    Double I_ESP_Ix5z_Fx2y_a = I_ESP_K2x5z_D2y_a+ABX*I_ESP_Ix5z_D2y_a;
    Double I_ESP_I5yz_Fx2y_a = I_ESP_Kx5yz_D2y_a+ABX*I_ESP_I5yz_D2y_a;
    Double I_ESP_I4y2z_Fx2y_a = I_ESP_Kx4y2z_D2y_a+ABX*I_ESP_I4y2z_D2y_a;
    Double I_ESP_I3y3z_Fx2y_a = I_ESP_Kx3y3z_D2y_a+ABX*I_ESP_I3y3z_D2y_a;
    Double I_ESP_I2y4z_Fx2y_a = I_ESP_Kx2y4z_D2y_a+ABX*I_ESP_I2y4z_D2y_a;
    Double I_ESP_Iy5z_Fx2y_a = I_ESP_Kxy5z_D2y_a+ABX*I_ESP_Iy5z_D2y_a;
    Double I_ESP_I6z_Fx2y_a = I_ESP_Kx6z_D2y_a+ABX*I_ESP_I6z_D2y_a;
    Double I_ESP_I5xy_Fx2z_a = I_ESP_K6xy_D2z_a+ABX*I_ESP_I5xy_D2z_a;
    Double I_ESP_I4x2y_Fx2z_a = I_ESP_K5x2y_D2z_a+ABX*I_ESP_I4x2y_D2z_a;
    Double I_ESP_I4xyz_Fx2z_a = I_ESP_K5xyz_D2z_a+ABX*I_ESP_I4xyz_D2z_a;
    Double I_ESP_I3x3y_Fx2z_a = I_ESP_K4x3y_D2z_a+ABX*I_ESP_I3x3y_D2z_a;
    Double I_ESP_I3x2yz_Fx2z_a = I_ESP_K4x2yz_D2z_a+ABX*I_ESP_I3x2yz_D2z_a;
    Double I_ESP_I3xy2z_Fx2z_a = I_ESP_K4xy2z_D2z_a+ABX*I_ESP_I3xy2z_D2z_a;
    Double I_ESP_I2x4y_Fx2z_a = I_ESP_K3x4y_D2z_a+ABX*I_ESP_I2x4y_D2z_a;
    Double I_ESP_I2x3yz_Fx2z_a = I_ESP_K3x3yz_D2z_a+ABX*I_ESP_I2x3yz_D2z_a;
    Double I_ESP_I2x2y2z_Fx2z_a = I_ESP_K3x2y2z_D2z_a+ABX*I_ESP_I2x2y2z_D2z_a;
    Double I_ESP_I2xy3z_Fx2z_a = I_ESP_K3xy3z_D2z_a+ABX*I_ESP_I2xy3z_D2z_a;
    Double I_ESP_Ix5y_Fx2z_a = I_ESP_K2x5y_D2z_a+ABX*I_ESP_Ix5y_D2z_a;
    Double I_ESP_Ix4yz_Fx2z_a = I_ESP_K2x4yz_D2z_a+ABX*I_ESP_Ix4yz_D2z_a;
    Double I_ESP_Ix3y2z_Fx2z_a = I_ESP_K2x3y2z_D2z_a+ABX*I_ESP_Ix3y2z_D2z_a;
    Double I_ESP_Ix2y3z_Fx2z_a = I_ESP_K2x2y3z_D2z_a+ABX*I_ESP_Ix2y3z_D2z_a;
    Double I_ESP_Ixy4z_Fx2z_a = I_ESP_K2xy4z_D2z_a+ABX*I_ESP_Ixy4z_D2z_a;
    Double I_ESP_I6y_Fx2z_a = I_ESP_Kx6y_D2z_a+ABX*I_ESP_I6y_D2z_a;
    Double I_ESP_I5yz_Fx2z_a = I_ESP_Kx5yz_D2z_a+ABX*I_ESP_I5yz_D2z_a;
    Double I_ESP_I4y2z_Fx2z_a = I_ESP_Kx4y2z_D2z_a+ABX*I_ESP_I4y2z_D2z_a;
    Double I_ESP_I3y3z_Fx2z_a = I_ESP_Kx3y3z_D2z_a+ABX*I_ESP_I3y3z_D2z_a;
    Double I_ESP_I2y4z_Fx2z_a = I_ESP_Kx2y4z_D2z_a+ABX*I_ESP_I2y4z_D2z_a;
    Double I_ESP_Iy5z_Fx2z_a = I_ESP_Kxy5z_D2z_a+ABX*I_ESP_Iy5z_D2z_a;
    Double I_ESP_I6x_F3y_a = I_ESP_K6xy_D2y_a+ABY*I_ESP_I6x_D2y_a;
    Double I_ESP_I5xy_F3y_a = I_ESP_K5x2y_D2y_a+ABY*I_ESP_I5xy_D2y_a;
    Double I_ESP_I5xz_F3y_a = I_ESP_K5xyz_D2y_a+ABY*I_ESP_I5xz_D2y_a;
    Double I_ESP_I4x2y_F3y_a = I_ESP_K4x3y_D2y_a+ABY*I_ESP_I4x2y_D2y_a;
    Double I_ESP_I4xyz_F3y_a = I_ESP_K4x2yz_D2y_a+ABY*I_ESP_I4xyz_D2y_a;
    Double I_ESP_I4x2z_F3y_a = I_ESP_K4xy2z_D2y_a+ABY*I_ESP_I4x2z_D2y_a;
    Double I_ESP_I3x3y_F3y_a = I_ESP_K3x4y_D2y_a+ABY*I_ESP_I3x3y_D2y_a;
    Double I_ESP_I3x2yz_F3y_a = I_ESP_K3x3yz_D2y_a+ABY*I_ESP_I3x2yz_D2y_a;
    Double I_ESP_I3xy2z_F3y_a = I_ESP_K3x2y2z_D2y_a+ABY*I_ESP_I3xy2z_D2y_a;
    Double I_ESP_I3x3z_F3y_a = I_ESP_K3xy3z_D2y_a+ABY*I_ESP_I3x3z_D2y_a;
    Double I_ESP_I2x4y_F3y_a = I_ESP_K2x5y_D2y_a+ABY*I_ESP_I2x4y_D2y_a;
    Double I_ESP_I2x3yz_F3y_a = I_ESP_K2x4yz_D2y_a+ABY*I_ESP_I2x3yz_D2y_a;
    Double I_ESP_I2x2y2z_F3y_a = I_ESP_K2x3y2z_D2y_a+ABY*I_ESP_I2x2y2z_D2y_a;
    Double I_ESP_I2xy3z_F3y_a = I_ESP_K2x2y3z_D2y_a+ABY*I_ESP_I2xy3z_D2y_a;
    Double I_ESP_I2x4z_F3y_a = I_ESP_K2xy4z_D2y_a+ABY*I_ESP_I2x4z_D2y_a;
    Double I_ESP_Ix5y_F3y_a = I_ESP_Kx6y_D2y_a+ABY*I_ESP_Ix5y_D2y_a;
    Double I_ESP_Ix4yz_F3y_a = I_ESP_Kx5yz_D2y_a+ABY*I_ESP_Ix4yz_D2y_a;
    Double I_ESP_Ix3y2z_F3y_a = I_ESP_Kx4y2z_D2y_a+ABY*I_ESP_Ix3y2z_D2y_a;
    Double I_ESP_Ix2y3z_F3y_a = I_ESP_Kx3y3z_D2y_a+ABY*I_ESP_Ix2y3z_D2y_a;
    Double I_ESP_Ixy4z_F3y_a = I_ESP_Kx2y4z_D2y_a+ABY*I_ESP_Ixy4z_D2y_a;
    Double I_ESP_Ix5z_F3y_a = I_ESP_Kxy5z_D2y_a+ABY*I_ESP_Ix5z_D2y_a;
    Double I_ESP_I6y_F3y_a = I_ESP_K7y_D2y_a+ABY*I_ESP_I6y_D2y_a;
    Double I_ESP_I5yz_F3y_a = I_ESP_K6yz_D2y_a+ABY*I_ESP_I5yz_D2y_a;
    Double I_ESP_I4y2z_F3y_a = I_ESP_K5y2z_D2y_a+ABY*I_ESP_I4y2z_D2y_a;
    Double I_ESP_I3y3z_F3y_a = I_ESP_K4y3z_D2y_a+ABY*I_ESP_I3y3z_D2y_a;
    Double I_ESP_I2y4z_F3y_a = I_ESP_K3y4z_D2y_a+ABY*I_ESP_I2y4z_D2y_a;
    Double I_ESP_Iy5z_F3y_a = I_ESP_K2y5z_D2y_a+ABY*I_ESP_Iy5z_D2y_a;
    Double I_ESP_I6z_F3y_a = I_ESP_Ky6z_D2y_a+ABY*I_ESP_I6z_D2y_a;
    Double I_ESP_I5xz_F2yz_a = I_ESP_K5x2z_D2y_a+ABZ*I_ESP_I5xz_D2y_a;
    Double I_ESP_I4xyz_F2yz_a = I_ESP_K4xy2z_D2y_a+ABZ*I_ESP_I4xyz_D2y_a;
    Double I_ESP_I4x2z_F2yz_a = I_ESP_K4x3z_D2y_a+ABZ*I_ESP_I4x2z_D2y_a;
    Double I_ESP_I3x2yz_F2yz_a = I_ESP_K3x2y2z_D2y_a+ABZ*I_ESP_I3x2yz_D2y_a;
    Double I_ESP_I3xy2z_F2yz_a = I_ESP_K3xy3z_D2y_a+ABZ*I_ESP_I3xy2z_D2y_a;
    Double I_ESP_I3x3z_F2yz_a = I_ESP_K3x4z_D2y_a+ABZ*I_ESP_I3x3z_D2y_a;
    Double I_ESP_I2x3yz_F2yz_a = I_ESP_K2x3y2z_D2y_a+ABZ*I_ESP_I2x3yz_D2y_a;
    Double I_ESP_I2x2y2z_F2yz_a = I_ESP_K2x2y3z_D2y_a+ABZ*I_ESP_I2x2y2z_D2y_a;
    Double I_ESP_I2xy3z_F2yz_a = I_ESP_K2xy4z_D2y_a+ABZ*I_ESP_I2xy3z_D2y_a;
    Double I_ESP_I2x4z_F2yz_a = I_ESP_K2x5z_D2y_a+ABZ*I_ESP_I2x4z_D2y_a;
    Double I_ESP_Ix4yz_F2yz_a = I_ESP_Kx4y2z_D2y_a+ABZ*I_ESP_Ix4yz_D2y_a;
    Double I_ESP_Ix3y2z_F2yz_a = I_ESP_Kx3y3z_D2y_a+ABZ*I_ESP_Ix3y2z_D2y_a;
    Double I_ESP_Ix2y3z_F2yz_a = I_ESP_Kx2y4z_D2y_a+ABZ*I_ESP_Ix2y3z_D2y_a;
    Double I_ESP_Ixy4z_F2yz_a = I_ESP_Kxy5z_D2y_a+ABZ*I_ESP_Ixy4z_D2y_a;
    Double I_ESP_Ix5z_F2yz_a = I_ESP_Kx6z_D2y_a+ABZ*I_ESP_Ix5z_D2y_a;
    Double I_ESP_I5yz_F2yz_a = I_ESP_K5y2z_D2y_a+ABZ*I_ESP_I5yz_D2y_a;
    Double I_ESP_I4y2z_F2yz_a = I_ESP_K4y3z_D2y_a+ABZ*I_ESP_I4y2z_D2y_a;
    Double I_ESP_I3y3z_F2yz_a = I_ESP_K3y4z_D2y_a+ABZ*I_ESP_I3y3z_D2y_a;
    Double I_ESP_I2y4z_F2yz_a = I_ESP_K2y5z_D2y_a+ABZ*I_ESP_I2y4z_D2y_a;
    Double I_ESP_Iy5z_F2yz_a = I_ESP_Ky6z_D2y_a+ABZ*I_ESP_Iy5z_D2y_a;
    Double I_ESP_I6z_F2yz_a = I_ESP_K7z_D2y_a+ABZ*I_ESP_I6z_D2y_a;
    Double I_ESP_I6x_F3z_a = I_ESP_K6xz_D2z_a+ABZ*I_ESP_I6x_D2z_a;
    Double I_ESP_I5xy_F3z_a = I_ESP_K5xyz_D2z_a+ABZ*I_ESP_I5xy_D2z_a;
    Double I_ESP_I5xz_F3z_a = I_ESP_K5x2z_D2z_a+ABZ*I_ESP_I5xz_D2z_a;
    Double I_ESP_I4x2y_F3z_a = I_ESP_K4x2yz_D2z_a+ABZ*I_ESP_I4x2y_D2z_a;
    Double I_ESP_I4xyz_F3z_a = I_ESP_K4xy2z_D2z_a+ABZ*I_ESP_I4xyz_D2z_a;
    Double I_ESP_I4x2z_F3z_a = I_ESP_K4x3z_D2z_a+ABZ*I_ESP_I4x2z_D2z_a;
    Double I_ESP_I3x3y_F3z_a = I_ESP_K3x3yz_D2z_a+ABZ*I_ESP_I3x3y_D2z_a;
    Double I_ESP_I3x2yz_F3z_a = I_ESP_K3x2y2z_D2z_a+ABZ*I_ESP_I3x2yz_D2z_a;
    Double I_ESP_I3xy2z_F3z_a = I_ESP_K3xy3z_D2z_a+ABZ*I_ESP_I3xy2z_D2z_a;
    Double I_ESP_I3x3z_F3z_a = I_ESP_K3x4z_D2z_a+ABZ*I_ESP_I3x3z_D2z_a;
    Double I_ESP_I2x4y_F3z_a = I_ESP_K2x4yz_D2z_a+ABZ*I_ESP_I2x4y_D2z_a;
    Double I_ESP_I2x3yz_F3z_a = I_ESP_K2x3y2z_D2z_a+ABZ*I_ESP_I2x3yz_D2z_a;
    Double I_ESP_I2x2y2z_F3z_a = I_ESP_K2x2y3z_D2z_a+ABZ*I_ESP_I2x2y2z_D2z_a;
    Double I_ESP_I2xy3z_F3z_a = I_ESP_K2xy4z_D2z_a+ABZ*I_ESP_I2xy3z_D2z_a;
    Double I_ESP_I2x4z_F3z_a = I_ESP_K2x5z_D2z_a+ABZ*I_ESP_I2x4z_D2z_a;
    Double I_ESP_Ix5y_F3z_a = I_ESP_Kx5yz_D2z_a+ABZ*I_ESP_Ix5y_D2z_a;
    Double I_ESP_Ix4yz_F3z_a = I_ESP_Kx4y2z_D2z_a+ABZ*I_ESP_Ix4yz_D2z_a;
    Double I_ESP_Ix3y2z_F3z_a = I_ESP_Kx3y3z_D2z_a+ABZ*I_ESP_Ix3y2z_D2z_a;
    Double I_ESP_Ix2y3z_F3z_a = I_ESP_Kx2y4z_D2z_a+ABZ*I_ESP_Ix2y3z_D2z_a;
    Double I_ESP_Ixy4z_F3z_a = I_ESP_Kxy5z_D2z_a+ABZ*I_ESP_Ixy4z_D2z_a;
    Double I_ESP_Ix5z_F3z_a = I_ESP_Kx6z_D2z_a+ABZ*I_ESP_Ix5z_D2z_a;
    Double I_ESP_I6y_F3z_a = I_ESP_K6yz_D2z_a+ABZ*I_ESP_I6y_D2z_a;
    Double I_ESP_I5yz_F3z_a = I_ESP_K5y2z_D2z_a+ABZ*I_ESP_I5yz_D2z_a;
    Double I_ESP_I4y2z_F3z_a = I_ESP_K4y3z_D2z_a+ABZ*I_ESP_I4y2z_D2z_a;
    Double I_ESP_I3y3z_F3z_a = I_ESP_K3y4z_D2z_a+ABZ*I_ESP_I3y3z_D2z_a;
    Double I_ESP_I2y4z_F3z_a = I_ESP_K2y5z_D2z_a+ABZ*I_ESP_I2y4z_D2z_a;
    Double I_ESP_Iy5z_F3z_a = I_ESP_Ky6z_D2z_a+ABZ*I_ESP_Iy5z_D2z_a;
    Double I_ESP_I6z_F3z_a = I_ESP_K7z_D2z_a+ABZ*I_ESP_I6z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_a
     * RHS shell quartet name: SQ_ESP_H_F_a
     ************************************************************/
    Double I_ESP_H5x_G4x_a = I_ESP_I6x_F3x_a+ABX*I_ESP_H5x_F3x_a;
    Double I_ESP_H4xy_G4x_a = I_ESP_I5xy_F3x_a+ABX*I_ESP_H4xy_F3x_a;
    Double I_ESP_H4xz_G4x_a = I_ESP_I5xz_F3x_a+ABX*I_ESP_H4xz_F3x_a;
    Double I_ESP_H3x2y_G4x_a = I_ESP_I4x2y_F3x_a+ABX*I_ESP_H3x2y_F3x_a;
    Double I_ESP_H3xyz_G4x_a = I_ESP_I4xyz_F3x_a+ABX*I_ESP_H3xyz_F3x_a;
    Double I_ESP_H3x2z_G4x_a = I_ESP_I4x2z_F3x_a+ABX*I_ESP_H3x2z_F3x_a;
    Double I_ESP_H2x3y_G4x_a = I_ESP_I3x3y_F3x_a+ABX*I_ESP_H2x3y_F3x_a;
    Double I_ESP_H2x2yz_G4x_a = I_ESP_I3x2yz_F3x_a+ABX*I_ESP_H2x2yz_F3x_a;
    Double I_ESP_H2xy2z_G4x_a = I_ESP_I3xy2z_F3x_a+ABX*I_ESP_H2xy2z_F3x_a;
    Double I_ESP_H2x3z_G4x_a = I_ESP_I3x3z_F3x_a+ABX*I_ESP_H2x3z_F3x_a;
    Double I_ESP_Hx4y_G4x_a = I_ESP_I2x4y_F3x_a+ABX*I_ESP_Hx4y_F3x_a;
    Double I_ESP_Hx3yz_G4x_a = I_ESP_I2x3yz_F3x_a+ABX*I_ESP_Hx3yz_F3x_a;
    Double I_ESP_Hx2y2z_G4x_a = I_ESP_I2x2y2z_F3x_a+ABX*I_ESP_Hx2y2z_F3x_a;
    Double I_ESP_Hxy3z_G4x_a = I_ESP_I2xy3z_F3x_a+ABX*I_ESP_Hxy3z_F3x_a;
    Double I_ESP_Hx4z_G4x_a = I_ESP_I2x4z_F3x_a+ABX*I_ESP_Hx4z_F3x_a;
    Double I_ESP_H5y_G4x_a = I_ESP_Ix5y_F3x_a+ABX*I_ESP_H5y_F3x_a;
    Double I_ESP_H4yz_G4x_a = I_ESP_Ix4yz_F3x_a+ABX*I_ESP_H4yz_F3x_a;
    Double I_ESP_H3y2z_G4x_a = I_ESP_Ix3y2z_F3x_a+ABX*I_ESP_H3y2z_F3x_a;
    Double I_ESP_H2y3z_G4x_a = I_ESP_Ix2y3z_F3x_a+ABX*I_ESP_H2y3z_F3x_a;
    Double I_ESP_Hy4z_G4x_a = I_ESP_Ixy4z_F3x_a+ABX*I_ESP_Hy4z_F3x_a;
    Double I_ESP_H5z_G4x_a = I_ESP_Ix5z_F3x_a+ABX*I_ESP_H5z_F3x_a;
    Double I_ESP_H5x_G3xy_a = I_ESP_I5xy_F3x_a+ABY*I_ESP_H5x_F3x_a;
    Double I_ESP_H4xy_G3xy_a = I_ESP_I4x2y_F3x_a+ABY*I_ESP_H4xy_F3x_a;
    Double I_ESP_H4xz_G3xy_a = I_ESP_I4xyz_F3x_a+ABY*I_ESP_H4xz_F3x_a;
    Double I_ESP_H3x2y_G3xy_a = I_ESP_I3x3y_F3x_a+ABY*I_ESP_H3x2y_F3x_a;
    Double I_ESP_H3xyz_G3xy_a = I_ESP_I3x2yz_F3x_a+ABY*I_ESP_H3xyz_F3x_a;
    Double I_ESP_H3x2z_G3xy_a = I_ESP_I3xy2z_F3x_a+ABY*I_ESP_H3x2z_F3x_a;
    Double I_ESP_H2x3y_G3xy_a = I_ESP_I2x4y_F3x_a+ABY*I_ESP_H2x3y_F3x_a;
    Double I_ESP_H2x2yz_G3xy_a = I_ESP_I2x3yz_F3x_a+ABY*I_ESP_H2x2yz_F3x_a;
    Double I_ESP_H2xy2z_G3xy_a = I_ESP_I2x2y2z_F3x_a+ABY*I_ESP_H2xy2z_F3x_a;
    Double I_ESP_H2x3z_G3xy_a = I_ESP_I2xy3z_F3x_a+ABY*I_ESP_H2x3z_F3x_a;
    Double I_ESP_Hx4y_G3xy_a = I_ESP_Ix5y_F3x_a+ABY*I_ESP_Hx4y_F3x_a;
    Double I_ESP_Hx3yz_G3xy_a = I_ESP_Ix4yz_F3x_a+ABY*I_ESP_Hx3yz_F3x_a;
    Double I_ESP_Hx2y2z_G3xy_a = I_ESP_Ix3y2z_F3x_a+ABY*I_ESP_Hx2y2z_F3x_a;
    Double I_ESP_Hxy3z_G3xy_a = I_ESP_Ix2y3z_F3x_a+ABY*I_ESP_Hxy3z_F3x_a;
    Double I_ESP_Hx4z_G3xy_a = I_ESP_Ixy4z_F3x_a+ABY*I_ESP_Hx4z_F3x_a;
    Double I_ESP_H5y_G3xy_a = I_ESP_I6y_F3x_a+ABY*I_ESP_H5y_F3x_a;
    Double I_ESP_H4yz_G3xy_a = I_ESP_I5yz_F3x_a+ABY*I_ESP_H4yz_F3x_a;
    Double I_ESP_H3y2z_G3xy_a = I_ESP_I4y2z_F3x_a+ABY*I_ESP_H3y2z_F3x_a;
    Double I_ESP_H2y3z_G3xy_a = I_ESP_I3y3z_F3x_a+ABY*I_ESP_H2y3z_F3x_a;
    Double I_ESP_Hy4z_G3xy_a = I_ESP_I2y4z_F3x_a+ABY*I_ESP_Hy4z_F3x_a;
    Double I_ESP_H5z_G3xy_a = I_ESP_Iy5z_F3x_a+ABY*I_ESP_H5z_F3x_a;
    Double I_ESP_H5x_G3xz_a = I_ESP_I5xz_F3x_a+ABZ*I_ESP_H5x_F3x_a;
    Double I_ESP_H4xy_G3xz_a = I_ESP_I4xyz_F3x_a+ABZ*I_ESP_H4xy_F3x_a;
    Double I_ESP_H4xz_G3xz_a = I_ESP_I4x2z_F3x_a+ABZ*I_ESP_H4xz_F3x_a;
    Double I_ESP_H3x2y_G3xz_a = I_ESP_I3x2yz_F3x_a+ABZ*I_ESP_H3x2y_F3x_a;
    Double I_ESP_H3xyz_G3xz_a = I_ESP_I3xy2z_F3x_a+ABZ*I_ESP_H3xyz_F3x_a;
    Double I_ESP_H3x2z_G3xz_a = I_ESP_I3x3z_F3x_a+ABZ*I_ESP_H3x2z_F3x_a;
    Double I_ESP_H2x3y_G3xz_a = I_ESP_I2x3yz_F3x_a+ABZ*I_ESP_H2x3y_F3x_a;
    Double I_ESP_H2x2yz_G3xz_a = I_ESP_I2x2y2z_F3x_a+ABZ*I_ESP_H2x2yz_F3x_a;
    Double I_ESP_H2xy2z_G3xz_a = I_ESP_I2xy3z_F3x_a+ABZ*I_ESP_H2xy2z_F3x_a;
    Double I_ESP_H2x3z_G3xz_a = I_ESP_I2x4z_F3x_a+ABZ*I_ESP_H2x3z_F3x_a;
    Double I_ESP_Hx4y_G3xz_a = I_ESP_Ix4yz_F3x_a+ABZ*I_ESP_Hx4y_F3x_a;
    Double I_ESP_Hx3yz_G3xz_a = I_ESP_Ix3y2z_F3x_a+ABZ*I_ESP_Hx3yz_F3x_a;
    Double I_ESP_Hx2y2z_G3xz_a = I_ESP_Ix2y3z_F3x_a+ABZ*I_ESP_Hx2y2z_F3x_a;
    Double I_ESP_Hxy3z_G3xz_a = I_ESP_Ixy4z_F3x_a+ABZ*I_ESP_Hxy3z_F3x_a;
    Double I_ESP_Hx4z_G3xz_a = I_ESP_Ix5z_F3x_a+ABZ*I_ESP_Hx4z_F3x_a;
    Double I_ESP_H5y_G3xz_a = I_ESP_I5yz_F3x_a+ABZ*I_ESP_H5y_F3x_a;
    Double I_ESP_H4yz_G3xz_a = I_ESP_I4y2z_F3x_a+ABZ*I_ESP_H4yz_F3x_a;
    Double I_ESP_H3y2z_G3xz_a = I_ESP_I3y3z_F3x_a+ABZ*I_ESP_H3y2z_F3x_a;
    Double I_ESP_H2y3z_G3xz_a = I_ESP_I2y4z_F3x_a+ABZ*I_ESP_H2y3z_F3x_a;
    Double I_ESP_Hy4z_G3xz_a = I_ESP_Iy5z_F3x_a+ABZ*I_ESP_Hy4z_F3x_a;
    Double I_ESP_H5z_G3xz_a = I_ESP_I6z_F3x_a+ABZ*I_ESP_H5z_F3x_a;
    Double I_ESP_H5x_G2x2y_a = I_ESP_I5xy_F2xy_a+ABY*I_ESP_H5x_F2xy_a;
    Double I_ESP_H4xy_G2x2y_a = I_ESP_I4x2y_F2xy_a+ABY*I_ESP_H4xy_F2xy_a;
    Double I_ESP_H4xz_G2x2y_a = I_ESP_I4xyz_F2xy_a+ABY*I_ESP_H4xz_F2xy_a;
    Double I_ESP_H3x2y_G2x2y_a = I_ESP_I3x3y_F2xy_a+ABY*I_ESP_H3x2y_F2xy_a;
    Double I_ESP_H3xyz_G2x2y_a = I_ESP_I3x2yz_F2xy_a+ABY*I_ESP_H3xyz_F2xy_a;
    Double I_ESP_H3x2z_G2x2y_a = I_ESP_I3xy2z_F2xy_a+ABY*I_ESP_H3x2z_F2xy_a;
    Double I_ESP_H2x3y_G2x2y_a = I_ESP_I2x4y_F2xy_a+ABY*I_ESP_H2x3y_F2xy_a;
    Double I_ESP_H2x2yz_G2x2y_a = I_ESP_I2x3yz_F2xy_a+ABY*I_ESP_H2x2yz_F2xy_a;
    Double I_ESP_H2xy2z_G2x2y_a = I_ESP_I2x2y2z_F2xy_a+ABY*I_ESP_H2xy2z_F2xy_a;
    Double I_ESP_H2x3z_G2x2y_a = I_ESP_I2xy3z_F2xy_a+ABY*I_ESP_H2x3z_F2xy_a;
    Double I_ESP_Hx4y_G2x2y_a = I_ESP_Ix5y_F2xy_a+ABY*I_ESP_Hx4y_F2xy_a;
    Double I_ESP_Hx3yz_G2x2y_a = I_ESP_Ix4yz_F2xy_a+ABY*I_ESP_Hx3yz_F2xy_a;
    Double I_ESP_Hx2y2z_G2x2y_a = I_ESP_Ix3y2z_F2xy_a+ABY*I_ESP_Hx2y2z_F2xy_a;
    Double I_ESP_Hxy3z_G2x2y_a = I_ESP_Ix2y3z_F2xy_a+ABY*I_ESP_Hxy3z_F2xy_a;
    Double I_ESP_Hx4z_G2x2y_a = I_ESP_Ixy4z_F2xy_a+ABY*I_ESP_Hx4z_F2xy_a;
    Double I_ESP_H5y_G2x2y_a = I_ESP_I6y_F2xy_a+ABY*I_ESP_H5y_F2xy_a;
    Double I_ESP_H4yz_G2x2y_a = I_ESP_I5yz_F2xy_a+ABY*I_ESP_H4yz_F2xy_a;
    Double I_ESP_H3y2z_G2x2y_a = I_ESP_I4y2z_F2xy_a+ABY*I_ESP_H3y2z_F2xy_a;
    Double I_ESP_H2y3z_G2x2y_a = I_ESP_I3y3z_F2xy_a+ABY*I_ESP_H2y3z_F2xy_a;
    Double I_ESP_Hy4z_G2x2y_a = I_ESP_I2y4z_F2xy_a+ABY*I_ESP_Hy4z_F2xy_a;
    Double I_ESP_H5z_G2x2y_a = I_ESP_Iy5z_F2xy_a+ABY*I_ESP_H5z_F2xy_a;
    Double I_ESP_H5x_G2xyz_a = I_ESP_I5xz_F2xy_a+ABZ*I_ESP_H5x_F2xy_a;
    Double I_ESP_H4xy_G2xyz_a = I_ESP_I4xyz_F2xy_a+ABZ*I_ESP_H4xy_F2xy_a;
    Double I_ESP_H4xz_G2xyz_a = I_ESP_I4x2z_F2xy_a+ABZ*I_ESP_H4xz_F2xy_a;
    Double I_ESP_H3x2y_G2xyz_a = I_ESP_I3x2yz_F2xy_a+ABZ*I_ESP_H3x2y_F2xy_a;
    Double I_ESP_H3xyz_G2xyz_a = I_ESP_I3xy2z_F2xy_a+ABZ*I_ESP_H3xyz_F2xy_a;
    Double I_ESP_H3x2z_G2xyz_a = I_ESP_I3x3z_F2xy_a+ABZ*I_ESP_H3x2z_F2xy_a;
    Double I_ESP_H2x3y_G2xyz_a = I_ESP_I2x3yz_F2xy_a+ABZ*I_ESP_H2x3y_F2xy_a;
    Double I_ESP_H2x2yz_G2xyz_a = I_ESP_I2x2y2z_F2xy_a+ABZ*I_ESP_H2x2yz_F2xy_a;
    Double I_ESP_H2xy2z_G2xyz_a = I_ESP_I2xy3z_F2xy_a+ABZ*I_ESP_H2xy2z_F2xy_a;
    Double I_ESP_H2x3z_G2xyz_a = I_ESP_I2x4z_F2xy_a+ABZ*I_ESP_H2x3z_F2xy_a;
    Double I_ESP_Hx4y_G2xyz_a = I_ESP_Ix4yz_F2xy_a+ABZ*I_ESP_Hx4y_F2xy_a;
    Double I_ESP_Hx3yz_G2xyz_a = I_ESP_Ix3y2z_F2xy_a+ABZ*I_ESP_Hx3yz_F2xy_a;
    Double I_ESP_Hx2y2z_G2xyz_a = I_ESP_Ix2y3z_F2xy_a+ABZ*I_ESP_Hx2y2z_F2xy_a;
    Double I_ESP_Hxy3z_G2xyz_a = I_ESP_Ixy4z_F2xy_a+ABZ*I_ESP_Hxy3z_F2xy_a;
    Double I_ESP_Hx4z_G2xyz_a = I_ESP_Ix5z_F2xy_a+ABZ*I_ESP_Hx4z_F2xy_a;
    Double I_ESP_H5y_G2xyz_a = I_ESP_I5yz_F2xy_a+ABZ*I_ESP_H5y_F2xy_a;
    Double I_ESP_H4yz_G2xyz_a = I_ESP_I4y2z_F2xy_a+ABZ*I_ESP_H4yz_F2xy_a;
    Double I_ESP_H3y2z_G2xyz_a = I_ESP_I3y3z_F2xy_a+ABZ*I_ESP_H3y2z_F2xy_a;
    Double I_ESP_H2y3z_G2xyz_a = I_ESP_I2y4z_F2xy_a+ABZ*I_ESP_H2y3z_F2xy_a;
    Double I_ESP_Hy4z_G2xyz_a = I_ESP_Iy5z_F2xy_a+ABZ*I_ESP_Hy4z_F2xy_a;
    Double I_ESP_H5z_G2xyz_a = I_ESP_I6z_F2xy_a+ABZ*I_ESP_H5z_F2xy_a;
    Double I_ESP_H5x_G2x2z_a = I_ESP_I5xz_F2xz_a+ABZ*I_ESP_H5x_F2xz_a;
    Double I_ESP_H4xy_G2x2z_a = I_ESP_I4xyz_F2xz_a+ABZ*I_ESP_H4xy_F2xz_a;
    Double I_ESP_H4xz_G2x2z_a = I_ESP_I4x2z_F2xz_a+ABZ*I_ESP_H4xz_F2xz_a;
    Double I_ESP_H3x2y_G2x2z_a = I_ESP_I3x2yz_F2xz_a+ABZ*I_ESP_H3x2y_F2xz_a;
    Double I_ESP_H3xyz_G2x2z_a = I_ESP_I3xy2z_F2xz_a+ABZ*I_ESP_H3xyz_F2xz_a;
    Double I_ESP_H3x2z_G2x2z_a = I_ESP_I3x3z_F2xz_a+ABZ*I_ESP_H3x2z_F2xz_a;
    Double I_ESP_H2x3y_G2x2z_a = I_ESP_I2x3yz_F2xz_a+ABZ*I_ESP_H2x3y_F2xz_a;
    Double I_ESP_H2x2yz_G2x2z_a = I_ESP_I2x2y2z_F2xz_a+ABZ*I_ESP_H2x2yz_F2xz_a;
    Double I_ESP_H2xy2z_G2x2z_a = I_ESP_I2xy3z_F2xz_a+ABZ*I_ESP_H2xy2z_F2xz_a;
    Double I_ESP_H2x3z_G2x2z_a = I_ESP_I2x4z_F2xz_a+ABZ*I_ESP_H2x3z_F2xz_a;
    Double I_ESP_Hx4y_G2x2z_a = I_ESP_Ix4yz_F2xz_a+ABZ*I_ESP_Hx4y_F2xz_a;
    Double I_ESP_Hx3yz_G2x2z_a = I_ESP_Ix3y2z_F2xz_a+ABZ*I_ESP_Hx3yz_F2xz_a;
    Double I_ESP_Hx2y2z_G2x2z_a = I_ESP_Ix2y3z_F2xz_a+ABZ*I_ESP_Hx2y2z_F2xz_a;
    Double I_ESP_Hxy3z_G2x2z_a = I_ESP_Ixy4z_F2xz_a+ABZ*I_ESP_Hxy3z_F2xz_a;
    Double I_ESP_Hx4z_G2x2z_a = I_ESP_Ix5z_F2xz_a+ABZ*I_ESP_Hx4z_F2xz_a;
    Double I_ESP_H5y_G2x2z_a = I_ESP_I5yz_F2xz_a+ABZ*I_ESP_H5y_F2xz_a;
    Double I_ESP_H4yz_G2x2z_a = I_ESP_I4y2z_F2xz_a+ABZ*I_ESP_H4yz_F2xz_a;
    Double I_ESP_H3y2z_G2x2z_a = I_ESP_I3y3z_F2xz_a+ABZ*I_ESP_H3y2z_F2xz_a;
    Double I_ESP_H2y3z_G2x2z_a = I_ESP_I2y4z_F2xz_a+ABZ*I_ESP_H2y3z_F2xz_a;
    Double I_ESP_Hy4z_G2x2z_a = I_ESP_Iy5z_F2xz_a+ABZ*I_ESP_Hy4z_F2xz_a;
    Double I_ESP_H5z_G2x2z_a = I_ESP_I6z_F2xz_a+ABZ*I_ESP_H5z_F2xz_a;
    Double I_ESP_H5x_Gx3y_a = I_ESP_I6x_F3y_a+ABX*I_ESP_H5x_F3y_a;
    Double I_ESP_H4xy_Gx3y_a = I_ESP_I5xy_F3y_a+ABX*I_ESP_H4xy_F3y_a;
    Double I_ESP_H4xz_Gx3y_a = I_ESP_I5xz_F3y_a+ABX*I_ESP_H4xz_F3y_a;
    Double I_ESP_H3x2y_Gx3y_a = I_ESP_I4x2y_F3y_a+ABX*I_ESP_H3x2y_F3y_a;
    Double I_ESP_H3xyz_Gx3y_a = I_ESP_I4xyz_F3y_a+ABX*I_ESP_H3xyz_F3y_a;
    Double I_ESP_H3x2z_Gx3y_a = I_ESP_I4x2z_F3y_a+ABX*I_ESP_H3x2z_F3y_a;
    Double I_ESP_H2x3y_Gx3y_a = I_ESP_I3x3y_F3y_a+ABX*I_ESP_H2x3y_F3y_a;
    Double I_ESP_H2x2yz_Gx3y_a = I_ESP_I3x2yz_F3y_a+ABX*I_ESP_H2x2yz_F3y_a;
    Double I_ESP_H2xy2z_Gx3y_a = I_ESP_I3xy2z_F3y_a+ABX*I_ESP_H2xy2z_F3y_a;
    Double I_ESP_H2x3z_Gx3y_a = I_ESP_I3x3z_F3y_a+ABX*I_ESP_H2x3z_F3y_a;
    Double I_ESP_Hx4y_Gx3y_a = I_ESP_I2x4y_F3y_a+ABX*I_ESP_Hx4y_F3y_a;
    Double I_ESP_Hx3yz_Gx3y_a = I_ESP_I2x3yz_F3y_a+ABX*I_ESP_Hx3yz_F3y_a;
    Double I_ESP_Hx2y2z_Gx3y_a = I_ESP_I2x2y2z_F3y_a+ABX*I_ESP_Hx2y2z_F3y_a;
    Double I_ESP_Hxy3z_Gx3y_a = I_ESP_I2xy3z_F3y_a+ABX*I_ESP_Hxy3z_F3y_a;
    Double I_ESP_Hx4z_Gx3y_a = I_ESP_I2x4z_F3y_a+ABX*I_ESP_Hx4z_F3y_a;
    Double I_ESP_H5y_Gx3y_a = I_ESP_Ix5y_F3y_a+ABX*I_ESP_H5y_F3y_a;
    Double I_ESP_H4yz_Gx3y_a = I_ESP_Ix4yz_F3y_a+ABX*I_ESP_H4yz_F3y_a;
    Double I_ESP_H3y2z_Gx3y_a = I_ESP_Ix3y2z_F3y_a+ABX*I_ESP_H3y2z_F3y_a;
    Double I_ESP_H2y3z_Gx3y_a = I_ESP_Ix2y3z_F3y_a+ABX*I_ESP_H2y3z_F3y_a;
    Double I_ESP_Hy4z_Gx3y_a = I_ESP_Ixy4z_F3y_a+ABX*I_ESP_Hy4z_F3y_a;
    Double I_ESP_H5z_Gx3y_a = I_ESP_Ix5z_F3y_a+ABX*I_ESP_H5z_F3y_a;
    Double I_ESP_H5x_Gx2yz_a = I_ESP_I5xz_Fx2y_a+ABZ*I_ESP_H5x_Fx2y_a;
    Double I_ESP_H4xy_Gx2yz_a = I_ESP_I4xyz_Fx2y_a+ABZ*I_ESP_H4xy_Fx2y_a;
    Double I_ESP_H4xz_Gx2yz_a = I_ESP_I4x2z_Fx2y_a+ABZ*I_ESP_H4xz_Fx2y_a;
    Double I_ESP_H3x2y_Gx2yz_a = I_ESP_I3x2yz_Fx2y_a+ABZ*I_ESP_H3x2y_Fx2y_a;
    Double I_ESP_H3xyz_Gx2yz_a = I_ESP_I3xy2z_Fx2y_a+ABZ*I_ESP_H3xyz_Fx2y_a;
    Double I_ESP_H3x2z_Gx2yz_a = I_ESP_I3x3z_Fx2y_a+ABZ*I_ESP_H3x2z_Fx2y_a;
    Double I_ESP_H2x3y_Gx2yz_a = I_ESP_I2x3yz_Fx2y_a+ABZ*I_ESP_H2x3y_Fx2y_a;
    Double I_ESP_H2x2yz_Gx2yz_a = I_ESP_I2x2y2z_Fx2y_a+ABZ*I_ESP_H2x2yz_Fx2y_a;
    Double I_ESP_H2xy2z_Gx2yz_a = I_ESP_I2xy3z_Fx2y_a+ABZ*I_ESP_H2xy2z_Fx2y_a;
    Double I_ESP_H2x3z_Gx2yz_a = I_ESP_I2x4z_Fx2y_a+ABZ*I_ESP_H2x3z_Fx2y_a;
    Double I_ESP_Hx4y_Gx2yz_a = I_ESP_Ix4yz_Fx2y_a+ABZ*I_ESP_Hx4y_Fx2y_a;
    Double I_ESP_Hx3yz_Gx2yz_a = I_ESP_Ix3y2z_Fx2y_a+ABZ*I_ESP_Hx3yz_Fx2y_a;
    Double I_ESP_Hx2y2z_Gx2yz_a = I_ESP_Ix2y3z_Fx2y_a+ABZ*I_ESP_Hx2y2z_Fx2y_a;
    Double I_ESP_Hxy3z_Gx2yz_a = I_ESP_Ixy4z_Fx2y_a+ABZ*I_ESP_Hxy3z_Fx2y_a;
    Double I_ESP_Hx4z_Gx2yz_a = I_ESP_Ix5z_Fx2y_a+ABZ*I_ESP_Hx4z_Fx2y_a;
    Double I_ESP_H5y_Gx2yz_a = I_ESP_I5yz_Fx2y_a+ABZ*I_ESP_H5y_Fx2y_a;
    Double I_ESP_H4yz_Gx2yz_a = I_ESP_I4y2z_Fx2y_a+ABZ*I_ESP_H4yz_Fx2y_a;
    Double I_ESP_H3y2z_Gx2yz_a = I_ESP_I3y3z_Fx2y_a+ABZ*I_ESP_H3y2z_Fx2y_a;
    Double I_ESP_H2y3z_Gx2yz_a = I_ESP_I2y4z_Fx2y_a+ABZ*I_ESP_H2y3z_Fx2y_a;
    Double I_ESP_Hy4z_Gx2yz_a = I_ESP_Iy5z_Fx2y_a+ABZ*I_ESP_Hy4z_Fx2y_a;
    Double I_ESP_H5z_Gx2yz_a = I_ESP_I6z_Fx2y_a+ABZ*I_ESP_H5z_Fx2y_a;
    Double I_ESP_H5x_Gxy2z_a = I_ESP_I5xy_Fx2z_a+ABY*I_ESP_H5x_Fx2z_a;
    Double I_ESP_H4xy_Gxy2z_a = I_ESP_I4x2y_Fx2z_a+ABY*I_ESP_H4xy_Fx2z_a;
    Double I_ESP_H4xz_Gxy2z_a = I_ESP_I4xyz_Fx2z_a+ABY*I_ESP_H4xz_Fx2z_a;
    Double I_ESP_H3x2y_Gxy2z_a = I_ESP_I3x3y_Fx2z_a+ABY*I_ESP_H3x2y_Fx2z_a;
    Double I_ESP_H3xyz_Gxy2z_a = I_ESP_I3x2yz_Fx2z_a+ABY*I_ESP_H3xyz_Fx2z_a;
    Double I_ESP_H3x2z_Gxy2z_a = I_ESP_I3xy2z_Fx2z_a+ABY*I_ESP_H3x2z_Fx2z_a;
    Double I_ESP_H2x3y_Gxy2z_a = I_ESP_I2x4y_Fx2z_a+ABY*I_ESP_H2x3y_Fx2z_a;
    Double I_ESP_H2x2yz_Gxy2z_a = I_ESP_I2x3yz_Fx2z_a+ABY*I_ESP_H2x2yz_Fx2z_a;
    Double I_ESP_H2xy2z_Gxy2z_a = I_ESP_I2x2y2z_Fx2z_a+ABY*I_ESP_H2xy2z_Fx2z_a;
    Double I_ESP_H2x3z_Gxy2z_a = I_ESP_I2xy3z_Fx2z_a+ABY*I_ESP_H2x3z_Fx2z_a;
    Double I_ESP_Hx4y_Gxy2z_a = I_ESP_Ix5y_Fx2z_a+ABY*I_ESP_Hx4y_Fx2z_a;
    Double I_ESP_Hx3yz_Gxy2z_a = I_ESP_Ix4yz_Fx2z_a+ABY*I_ESP_Hx3yz_Fx2z_a;
    Double I_ESP_Hx2y2z_Gxy2z_a = I_ESP_Ix3y2z_Fx2z_a+ABY*I_ESP_Hx2y2z_Fx2z_a;
    Double I_ESP_Hxy3z_Gxy2z_a = I_ESP_Ix2y3z_Fx2z_a+ABY*I_ESP_Hxy3z_Fx2z_a;
    Double I_ESP_Hx4z_Gxy2z_a = I_ESP_Ixy4z_Fx2z_a+ABY*I_ESP_Hx4z_Fx2z_a;
    Double I_ESP_H5y_Gxy2z_a = I_ESP_I6y_Fx2z_a+ABY*I_ESP_H5y_Fx2z_a;
    Double I_ESP_H4yz_Gxy2z_a = I_ESP_I5yz_Fx2z_a+ABY*I_ESP_H4yz_Fx2z_a;
    Double I_ESP_H3y2z_Gxy2z_a = I_ESP_I4y2z_Fx2z_a+ABY*I_ESP_H3y2z_Fx2z_a;
    Double I_ESP_H2y3z_Gxy2z_a = I_ESP_I3y3z_Fx2z_a+ABY*I_ESP_H2y3z_Fx2z_a;
    Double I_ESP_Hy4z_Gxy2z_a = I_ESP_I2y4z_Fx2z_a+ABY*I_ESP_Hy4z_Fx2z_a;
    Double I_ESP_H5z_Gxy2z_a = I_ESP_Iy5z_Fx2z_a+ABY*I_ESP_H5z_Fx2z_a;
    Double I_ESP_H5x_Gx3z_a = I_ESP_I6x_F3z_a+ABX*I_ESP_H5x_F3z_a;
    Double I_ESP_H4xy_Gx3z_a = I_ESP_I5xy_F3z_a+ABX*I_ESP_H4xy_F3z_a;
    Double I_ESP_H4xz_Gx3z_a = I_ESP_I5xz_F3z_a+ABX*I_ESP_H4xz_F3z_a;
    Double I_ESP_H3x2y_Gx3z_a = I_ESP_I4x2y_F3z_a+ABX*I_ESP_H3x2y_F3z_a;
    Double I_ESP_H3xyz_Gx3z_a = I_ESP_I4xyz_F3z_a+ABX*I_ESP_H3xyz_F3z_a;
    Double I_ESP_H3x2z_Gx3z_a = I_ESP_I4x2z_F3z_a+ABX*I_ESP_H3x2z_F3z_a;
    Double I_ESP_H2x3y_Gx3z_a = I_ESP_I3x3y_F3z_a+ABX*I_ESP_H2x3y_F3z_a;
    Double I_ESP_H2x2yz_Gx3z_a = I_ESP_I3x2yz_F3z_a+ABX*I_ESP_H2x2yz_F3z_a;
    Double I_ESP_H2xy2z_Gx3z_a = I_ESP_I3xy2z_F3z_a+ABX*I_ESP_H2xy2z_F3z_a;
    Double I_ESP_H2x3z_Gx3z_a = I_ESP_I3x3z_F3z_a+ABX*I_ESP_H2x3z_F3z_a;
    Double I_ESP_Hx4y_Gx3z_a = I_ESP_I2x4y_F3z_a+ABX*I_ESP_Hx4y_F3z_a;
    Double I_ESP_Hx3yz_Gx3z_a = I_ESP_I2x3yz_F3z_a+ABX*I_ESP_Hx3yz_F3z_a;
    Double I_ESP_Hx2y2z_Gx3z_a = I_ESP_I2x2y2z_F3z_a+ABX*I_ESP_Hx2y2z_F3z_a;
    Double I_ESP_Hxy3z_Gx3z_a = I_ESP_I2xy3z_F3z_a+ABX*I_ESP_Hxy3z_F3z_a;
    Double I_ESP_Hx4z_Gx3z_a = I_ESP_I2x4z_F3z_a+ABX*I_ESP_Hx4z_F3z_a;
    Double I_ESP_H5y_Gx3z_a = I_ESP_Ix5y_F3z_a+ABX*I_ESP_H5y_F3z_a;
    Double I_ESP_H4yz_Gx3z_a = I_ESP_Ix4yz_F3z_a+ABX*I_ESP_H4yz_F3z_a;
    Double I_ESP_H3y2z_Gx3z_a = I_ESP_Ix3y2z_F3z_a+ABX*I_ESP_H3y2z_F3z_a;
    Double I_ESP_H2y3z_Gx3z_a = I_ESP_Ix2y3z_F3z_a+ABX*I_ESP_H2y3z_F3z_a;
    Double I_ESP_Hy4z_Gx3z_a = I_ESP_Ixy4z_F3z_a+ABX*I_ESP_Hy4z_F3z_a;
    Double I_ESP_H5z_Gx3z_a = I_ESP_Ix5z_F3z_a+ABX*I_ESP_H5z_F3z_a;
    Double I_ESP_H5x_G4y_a = I_ESP_I5xy_F3y_a+ABY*I_ESP_H5x_F3y_a;
    Double I_ESP_H4xy_G4y_a = I_ESP_I4x2y_F3y_a+ABY*I_ESP_H4xy_F3y_a;
    Double I_ESP_H4xz_G4y_a = I_ESP_I4xyz_F3y_a+ABY*I_ESP_H4xz_F3y_a;
    Double I_ESP_H3x2y_G4y_a = I_ESP_I3x3y_F3y_a+ABY*I_ESP_H3x2y_F3y_a;
    Double I_ESP_H3xyz_G4y_a = I_ESP_I3x2yz_F3y_a+ABY*I_ESP_H3xyz_F3y_a;
    Double I_ESP_H3x2z_G4y_a = I_ESP_I3xy2z_F3y_a+ABY*I_ESP_H3x2z_F3y_a;
    Double I_ESP_H2x3y_G4y_a = I_ESP_I2x4y_F3y_a+ABY*I_ESP_H2x3y_F3y_a;
    Double I_ESP_H2x2yz_G4y_a = I_ESP_I2x3yz_F3y_a+ABY*I_ESP_H2x2yz_F3y_a;
    Double I_ESP_H2xy2z_G4y_a = I_ESP_I2x2y2z_F3y_a+ABY*I_ESP_H2xy2z_F3y_a;
    Double I_ESP_H2x3z_G4y_a = I_ESP_I2xy3z_F3y_a+ABY*I_ESP_H2x3z_F3y_a;
    Double I_ESP_Hx4y_G4y_a = I_ESP_Ix5y_F3y_a+ABY*I_ESP_Hx4y_F3y_a;
    Double I_ESP_Hx3yz_G4y_a = I_ESP_Ix4yz_F3y_a+ABY*I_ESP_Hx3yz_F3y_a;
    Double I_ESP_Hx2y2z_G4y_a = I_ESP_Ix3y2z_F3y_a+ABY*I_ESP_Hx2y2z_F3y_a;
    Double I_ESP_Hxy3z_G4y_a = I_ESP_Ix2y3z_F3y_a+ABY*I_ESP_Hxy3z_F3y_a;
    Double I_ESP_Hx4z_G4y_a = I_ESP_Ixy4z_F3y_a+ABY*I_ESP_Hx4z_F3y_a;
    Double I_ESP_H5y_G4y_a = I_ESP_I6y_F3y_a+ABY*I_ESP_H5y_F3y_a;
    Double I_ESP_H4yz_G4y_a = I_ESP_I5yz_F3y_a+ABY*I_ESP_H4yz_F3y_a;
    Double I_ESP_H3y2z_G4y_a = I_ESP_I4y2z_F3y_a+ABY*I_ESP_H3y2z_F3y_a;
    Double I_ESP_H2y3z_G4y_a = I_ESP_I3y3z_F3y_a+ABY*I_ESP_H2y3z_F3y_a;
    Double I_ESP_Hy4z_G4y_a = I_ESP_I2y4z_F3y_a+ABY*I_ESP_Hy4z_F3y_a;
    Double I_ESP_H5z_G4y_a = I_ESP_Iy5z_F3y_a+ABY*I_ESP_H5z_F3y_a;
    Double I_ESP_H5x_G3yz_a = I_ESP_I5xz_F3y_a+ABZ*I_ESP_H5x_F3y_a;
    Double I_ESP_H4xy_G3yz_a = I_ESP_I4xyz_F3y_a+ABZ*I_ESP_H4xy_F3y_a;
    Double I_ESP_H4xz_G3yz_a = I_ESP_I4x2z_F3y_a+ABZ*I_ESP_H4xz_F3y_a;
    Double I_ESP_H3x2y_G3yz_a = I_ESP_I3x2yz_F3y_a+ABZ*I_ESP_H3x2y_F3y_a;
    Double I_ESP_H3xyz_G3yz_a = I_ESP_I3xy2z_F3y_a+ABZ*I_ESP_H3xyz_F3y_a;
    Double I_ESP_H3x2z_G3yz_a = I_ESP_I3x3z_F3y_a+ABZ*I_ESP_H3x2z_F3y_a;
    Double I_ESP_H2x3y_G3yz_a = I_ESP_I2x3yz_F3y_a+ABZ*I_ESP_H2x3y_F3y_a;
    Double I_ESP_H2x2yz_G3yz_a = I_ESP_I2x2y2z_F3y_a+ABZ*I_ESP_H2x2yz_F3y_a;
    Double I_ESP_H2xy2z_G3yz_a = I_ESP_I2xy3z_F3y_a+ABZ*I_ESP_H2xy2z_F3y_a;
    Double I_ESP_H2x3z_G3yz_a = I_ESP_I2x4z_F3y_a+ABZ*I_ESP_H2x3z_F3y_a;
    Double I_ESP_Hx4y_G3yz_a = I_ESP_Ix4yz_F3y_a+ABZ*I_ESP_Hx4y_F3y_a;
    Double I_ESP_Hx3yz_G3yz_a = I_ESP_Ix3y2z_F3y_a+ABZ*I_ESP_Hx3yz_F3y_a;
    Double I_ESP_Hx2y2z_G3yz_a = I_ESP_Ix2y3z_F3y_a+ABZ*I_ESP_Hx2y2z_F3y_a;
    Double I_ESP_Hxy3z_G3yz_a = I_ESP_Ixy4z_F3y_a+ABZ*I_ESP_Hxy3z_F3y_a;
    Double I_ESP_Hx4z_G3yz_a = I_ESP_Ix5z_F3y_a+ABZ*I_ESP_Hx4z_F3y_a;
    Double I_ESP_H5y_G3yz_a = I_ESP_I5yz_F3y_a+ABZ*I_ESP_H5y_F3y_a;
    Double I_ESP_H4yz_G3yz_a = I_ESP_I4y2z_F3y_a+ABZ*I_ESP_H4yz_F3y_a;
    Double I_ESP_H3y2z_G3yz_a = I_ESP_I3y3z_F3y_a+ABZ*I_ESP_H3y2z_F3y_a;
    Double I_ESP_H2y3z_G3yz_a = I_ESP_I2y4z_F3y_a+ABZ*I_ESP_H2y3z_F3y_a;
    Double I_ESP_Hy4z_G3yz_a = I_ESP_Iy5z_F3y_a+ABZ*I_ESP_Hy4z_F3y_a;
    Double I_ESP_H5z_G3yz_a = I_ESP_I6z_F3y_a+ABZ*I_ESP_H5z_F3y_a;
    Double I_ESP_H5x_G2y2z_a = I_ESP_I5xz_F2yz_a+ABZ*I_ESP_H5x_F2yz_a;
    Double I_ESP_H4xy_G2y2z_a = I_ESP_I4xyz_F2yz_a+ABZ*I_ESP_H4xy_F2yz_a;
    Double I_ESP_H4xz_G2y2z_a = I_ESP_I4x2z_F2yz_a+ABZ*I_ESP_H4xz_F2yz_a;
    Double I_ESP_H3x2y_G2y2z_a = I_ESP_I3x2yz_F2yz_a+ABZ*I_ESP_H3x2y_F2yz_a;
    Double I_ESP_H3xyz_G2y2z_a = I_ESP_I3xy2z_F2yz_a+ABZ*I_ESP_H3xyz_F2yz_a;
    Double I_ESP_H3x2z_G2y2z_a = I_ESP_I3x3z_F2yz_a+ABZ*I_ESP_H3x2z_F2yz_a;
    Double I_ESP_H2x3y_G2y2z_a = I_ESP_I2x3yz_F2yz_a+ABZ*I_ESP_H2x3y_F2yz_a;
    Double I_ESP_H2x2yz_G2y2z_a = I_ESP_I2x2y2z_F2yz_a+ABZ*I_ESP_H2x2yz_F2yz_a;
    Double I_ESP_H2xy2z_G2y2z_a = I_ESP_I2xy3z_F2yz_a+ABZ*I_ESP_H2xy2z_F2yz_a;
    Double I_ESP_H2x3z_G2y2z_a = I_ESP_I2x4z_F2yz_a+ABZ*I_ESP_H2x3z_F2yz_a;
    Double I_ESP_Hx4y_G2y2z_a = I_ESP_Ix4yz_F2yz_a+ABZ*I_ESP_Hx4y_F2yz_a;
    Double I_ESP_Hx3yz_G2y2z_a = I_ESP_Ix3y2z_F2yz_a+ABZ*I_ESP_Hx3yz_F2yz_a;
    Double I_ESP_Hx2y2z_G2y2z_a = I_ESP_Ix2y3z_F2yz_a+ABZ*I_ESP_Hx2y2z_F2yz_a;
    Double I_ESP_Hxy3z_G2y2z_a = I_ESP_Ixy4z_F2yz_a+ABZ*I_ESP_Hxy3z_F2yz_a;
    Double I_ESP_Hx4z_G2y2z_a = I_ESP_Ix5z_F2yz_a+ABZ*I_ESP_Hx4z_F2yz_a;
    Double I_ESP_H5y_G2y2z_a = I_ESP_I5yz_F2yz_a+ABZ*I_ESP_H5y_F2yz_a;
    Double I_ESP_H4yz_G2y2z_a = I_ESP_I4y2z_F2yz_a+ABZ*I_ESP_H4yz_F2yz_a;
    Double I_ESP_H3y2z_G2y2z_a = I_ESP_I3y3z_F2yz_a+ABZ*I_ESP_H3y2z_F2yz_a;
    Double I_ESP_H2y3z_G2y2z_a = I_ESP_I2y4z_F2yz_a+ABZ*I_ESP_H2y3z_F2yz_a;
    Double I_ESP_Hy4z_G2y2z_a = I_ESP_Iy5z_F2yz_a+ABZ*I_ESP_Hy4z_F2yz_a;
    Double I_ESP_H5z_G2y2z_a = I_ESP_I6z_F2yz_a+ABZ*I_ESP_H5z_F2yz_a;
    Double I_ESP_H5x_Gy3z_a = I_ESP_I5xy_F3z_a+ABY*I_ESP_H5x_F3z_a;
    Double I_ESP_H4xy_Gy3z_a = I_ESP_I4x2y_F3z_a+ABY*I_ESP_H4xy_F3z_a;
    Double I_ESP_H4xz_Gy3z_a = I_ESP_I4xyz_F3z_a+ABY*I_ESP_H4xz_F3z_a;
    Double I_ESP_H3x2y_Gy3z_a = I_ESP_I3x3y_F3z_a+ABY*I_ESP_H3x2y_F3z_a;
    Double I_ESP_H3xyz_Gy3z_a = I_ESP_I3x2yz_F3z_a+ABY*I_ESP_H3xyz_F3z_a;
    Double I_ESP_H3x2z_Gy3z_a = I_ESP_I3xy2z_F3z_a+ABY*I_ESP_H3x2z_F3z_a;
    Double I_ESP_H2x3y_Gy3z_a = I_ESP_I2x4y_F3z_a+ABY*I_ESP_H2x3y_F3z_a;
    Double I_ESP_H2x2yz_Gy3z_a = I_ESP_I2x3yz_F3z_a+ABY*I_ESP_H2x2yz_F3z_a;
    Double I_ESP_H2xy2z_Gy3z_a = I_ESP_I2x2y2z_F3z_a+ABY*I_ESP_H2xy2z_F3z_a;
    Double I_ESP_H2x3z_Gy3z_a = I_ESP_I2xy3z_F3z_a+ABY*I_ESP_H2x3z_F3z_a;
    Double I_ESP_Hx4y_Gy3z_a = I_ESP_Ix5y_F3z_a+ABY*I_ESP_Hx4y_F3z_a;
    Double I_ESP_Hx3yz_Gy3z_a = I_ESP_Ix4yz_F3z_a+ABY*I_ESP_Hx3yz_F3z_a;
    Double I_ESP_Hx2y2z_Gy3z_a = I_ESP_Ix3y2z_F3z_a+ABY*I_ESP_Hx2y2z_F3z_a;
    Double I_ESP_Hxy3z_Gy3z_a = I_ESP_Ix2y3z_F3z_a+ABY*I_ESP_Hxy3z_F3z_a;
    Double I_ESP_Hx4z_Gy3z_a = I_ESP_Ixy4z_F3z_a+ABY*I_ESP_Hx4z_F3z_a;
    Double I_ESP_H5y_Gy3z_a = I_ESP_I6y_F3z_a+ABY*I_ESP_H5y_F3z_a;
    Double I_ESP_H4yz_Gy3z_a = I_ESP_I5yz_F3z_a+ABY*I_ESP_H4yz_F3z_a;
    Double I_ESP_H3y2z_Gy3z_a = I_ESP_I4y2z_F3z_a+ABY*I_ESP_H3y2z_F3z_a;
    Double I_ESP_H2y3z_Gy3z_a = I_ESP_I3y3z_F3z_a+ABY*I_ESP_H2y3z_F3z_a;
    Double I_ESP_Hy4z_Gy3z_a = I_ESP_I2y4z_F3z_a+ABY*I_ESP_Hy4z_F3z_a;
    Double I_ESP_H5z_Gy3z_a = I_ESP_Iy5z_F3z_a+ABY*I_ESP_H5z_F3z_a;
    Double I_ESP_H5x_G4z_a = I_ESP_I5xz_F3z_a+ABZ*I_ESP_H5x_F3z_a;
    Double I_ESP_H4xy_G4z_a = I_ESP_I4xyz_F3z_a+ABZ*I_ESP_H4xy_F3z_a;
    Double I_ESP_H4xz_G4z_a = I_ESP_I4x2z_F3z_a+ABZ*I_ESP_H4xz_F3z_a;
    Double I_ESP_H3x2y_G4z_a = I_ESP_I3x2yz_F3z_a+ABZ*I_ESP_H3x2y_F3z_a;
    Double I_ESP_H3xyz_G4z_a = I_ESP_I3xy2z_F3z_a+ABZ*I_ESP_H3xyz_F3z_a;
    Double I_ESP_H3x2z_G4z_a = I_ESP_I3x3z_F3z_a+ABZ*I_ESP_H3x2z_F3z_a;
    Double I_ESP_H2x3y_G4z_a = I_ESP_I2x3yz_F3z_a+ABZ*I_ESP_H2x3y_F3z_a;
    Double I_ESP_H2x2yz_G4z_a = I_ESP_I2x2y2z_F3z_a+ABZ*I_ESP_H2x2yz_F3z_a;
    Double I_ESP_H2xy2z_G4z_a = I_ESP_I2xy3z_F3z_a+ABZ*I_ESP_H2xy2z_F3z_a;
    Double I_ESP_H2x3z_G4z_a = I_ESP_I2x4z_F3z_a+ABZ*I_ESP_H2x3z_F3z_a;
    Double I_ESP_Hx4y_G4z_a = I_ESP_Ix4yz_F3z_a+ABZ*I_ESP_Hx4y_F3z_a;
    Double I_ESP_Hx3yz_G4z_a = I_ESP_Ix3y2z_F3z_a+ABZ*I_ESP_Hx3yz_F3z_a;
    Double I_ESP_Hx2y2z_G4z_a = I_ESP_Ix2y3z_F3z_a+ABZ*I_ESP_Hx2y2z_F3z_a;
    Double I_ESP_Hxy3z_G4z_a = I_ESP_Ixy4z_F3z_a+ABZ*I_ESP_Hxy3z_F3z_a;
    Double I_ESP_Hx4z_G4z_a = I_ESP_Ix5z_F3z_a+ABZ*I_ESP_Hx4z_F3z_a;
    Double I_ESP_H5y_G4z_a = I_ESP_I5yz_F3z_a+ABZ*I_ESP_H5y_F3z_a;
    Double I_ESP_H4yz_G4z_a = I_ESP_I4y2z_F3z_a+ABZ*I_ESP_H4yz_F3z_a;
    Double I_ESP_H3y2z_G4z_a = I_ESP_I3y3z_F3z_a+ABZ*I_ESP_H3y2z_F3z_a;
    Double I_ESP_H2y3z_G4z_a = I_ESP_I2y4z_F3z_a+ABZ*I_ESP_H2y3z_F3z_a;
    Double I_ESP_Hy4z_G4z_a = I_ESP_Iy5z_F3z_a+ABZ*I_ESP_Hy4z_F3z_a;
    Double I_ESP_H5z_G4z_a = I_ESP_I6z_F3z_a+ABZ*I_ESP_H5z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_F_G
     ************************************************************/
    abcd[iGrid*675+0] = 2.0E0*I_ESP_H5x_G4x_a-4*I_ESP_F3x_G4x;
    abcd[iGrid*675+1] = 2.0E0*I_ESP_H4xy_G4x_a-3*I_ESP_F2xy_G4x;
    abcd[iGrid*675+2] = 2.0E0*I_ESP_H4xz_G4x_a-3*I_ESP_F2xz_G4x;
    abcd[iGrid*675+3] = 2.0E0*I_ESP_H3x2y_G4x_a-2*I_ESP_Fx2y_G4x;
    abcd[iGrid*675+4] = 2.0E0*I_ESP_H3xyz_G4x_a-2*I_ESP_Fxyz_G4x;
    abcd[iGrid*675+5] = 2.0E0*I_ESP_H3x2z_G4x_a-2*I_ESP_Fx2z_G4x;
    abcd[iGrid*675+6] = 2.0E0*I_ESP_H2x3y_G4x_a-1*I_ESP_F3y_G4x;
    abcd[iGrid*675+7] = 2.0E0*I_ESP_H2x2yz_G4x_a-1*I_ESP_F2yz_G4x;
    abcd[iGrid*675+8] = 2.0E0*I_ESP_H2xy2z_G4x_a-1*I_ESP_Fy2z_G4x;
    abcd[iGrid*675+9] = 2.0E0*I_ESP_H2x3z_G4x_a-1*I_ESP_F3z_G4x;
    abcd[iGrid*675+10] = 2.0E0*I_ESP_Hx4y_G4x_a;
    abcd[iGrid*675+11] = 2.0E0*I_ESP_Hx3yz_G4x_a;
    abcd[iGrid*675+12] = 2.0E0*I_ESP_Hx2y2z_G4x_a;
    abcd[iGrid*675+13] = 2.0E0*I_ESP_Hxy3z_G4x_a;
    abcd[iGrid*675+14] = 2.0E0*I_ESP_Hx4z_G4x_a;
    abcd[iGrid*675+15] = 2.0E0*I_ESP_H5x_G3xy_a-4*I_ESP_F3x_G3xy;
    abcd[iGrid*675+16] = 2.0E0*I_ESP_H4xy_G3xy_a-3*I_ESP_F2xy_G3xy;
    abcd[iGrid*675+17] = 2.0E0*I_ESP_H4xz_G3xy_a-3*I_ESP_F2xz_G3xy;
    abcd[iGrid*675+18] = 2.0E0*I_ESP_H3x2y_G3xy_a-2*I_ESP_Fx2y_G3xy;
    abcd[iGrid*675+19] = 2.0E0*I_ESP_H3xyz_G3xy_a-2*I_ESP_Fxyz_G3xy;
    abcd[iGrid*675+20] = 2.0E0*I_ESP_H3x2z_G3xy_a-2*I_ESP_Fx2z_G3xy;
    abcd[iGrid*675+21] = 2.0E0*I_ESP_H2x3y_G3xy_a-1*I_ESP_F3y_G3xy;
    abcd[iGrid*675+22] = 2.0E0*I_ESP_H2x2yz_G3xy_a-1*I_ESP_F2yz_G3xy;
    abcd[iGrid*675+23] = 2.0E0*I_ESP_H2xy2z_G3xy_a-1*I_ESP_Fy2z_G3xy;
    abcd[iGrid*675+24] = 2.0E0*I_ESP_H2x3z_G3xy_a-1*I_ESP_F3z_G3xy;
    abcd[iGrid*675+25] = 2.0E0*I_ESP_Hx4y_G3xy_a;
    abcd[iGrid*675+26] = 2.0E0*I_ESP_Hx3yz_G3xy_a;
    abcd[iGrid*675+27] = 2.0E0*I_ESP_Hx2y2z_G3xy_a;
    abcd[iGrid*675+28] = 2.0E0*I_ESP_Hxy3z_G3xy_a;
    abcd[iGrid*675+29] = 2.0E0*I_ESP_Hx4z_G3xy_a;
    abcd[iGrid*675+30] = 2.0E0*I_ESP_H5x_G3xz_a-4*I_ESP_F3x_G3xz;
    abcd[iGrid*675+31] = 2.0E0*I_ESP_H4xy_G3xz_a-3*I_ESP_F2xy_G3xz;
    abcd[iGrid*675+32] = 2.0E0*I_ESP_H4xz_G3xz_a-3*I_ESP_F2xz_G3xz;
    abcd[iGrid*675+33] = 2.0E0*I_ESP_H3x2y_G3xz_a-2*I_ESP_Fx2y_G3xz;
    abcd[iGrid*675+34] = 2.0E0*I_ESP_H3xyz_G3xz_a-2*I_ESP_Fxyz_G3xz;
    abcd[iGrid*675+35] = 2.0E0*I_ESP_H3x2z_G3xz_a-2*I_ESP_Fx2z_G3xz;
    abcd[iGrid*675+36] = 2.0E0*I_ESP_H2x3y_G3xz_a-1*I_ESP_F3y_G3xz;
    abcd[iGrid*675+37] = 2.0E0*I_ESP_H2x2yz_G3xz_a-1*I_ESP_F2yz_G3xz;
    abcd[iGrid*675+38] = 2.0E0*I_ESP_H2xy2z_G3xz_a-1*I_ESP_Fy2z_G3xz;
    abcd[iGrid*675+39] = 2.0E0*I_ESP_H2x3z_G3xz_a-1*I_ESP_F3z_G3xz;
    abcd[iGrid*675+40] = 2.0E0*I_ESP_Hx4y_G3xz_a;
    abcd[iGrid*675+41] = 2.0E0*I_ESP_Hx3yz_G3xz_a;
    abcd[iGrid*675+42] = 2.0E0*I_ESP_Hx2y2z_G3xz_a;
    abcd[iGrid*675+43] = 2.0E0*I_ESP_Hxy3z_G3xz_a;
    abcd[iGrid*675+44] = 2.0E0*I_ESP_Hx4z_G3xz_a;
    abcd[iGrid*675+45] = 2.0E0*I_ESP_H5x_G2x2y_a-4*I_ESP_F3x_G2x2y;
    abcd[iGrid*675+46] = 2.0E0*I_ESP_H4xy_G2x2y_a-3*I_ESP_F2xy_G2x2y;
    abcd[iGrid*675+47] = 2.0E0*I_ESP_H4xz_G2x2y_a-3*I_ESP_F2xz_G2x2y;
    abcd[iGrid*675+48] = 2.0E0*I_ESP_H3x2y_G2x2y_a-2*I_ESP_Fx2y_G2x2y;
    abcd[iGrid*675+49] = 2.0E0*I_ESP_H3xyz_G2x2y_a-2*I_ESP_Fxyz_G2x2y;
    abcd[iGrid*675+50] = 2.0E0*I_ESP_H3x2z_G2x2y_a-2*I_ESP_Fx2z_G2x2y;
    abcd[iGrid*675+51] = 2.0E0*I_ESP_H2x3y_G2x2y_a-1*I_ESP_F3y_G2x2y;
    abcd[iGrid*675+52] = 2.0E0*I_ESP_H2x2yz_G2x2y_a-1*I_ESP_F2yz_G2x2y;
    abcd[iGrid*675+53] = 2.0E0*I_ESP_H2xy2z_G2x2y_a-1*I_ESP_Fy2z_G2x2y;
    abcd[iGrid*675+54] = 2.0E0*I_ESP_H2x3z_G2x2y_a-1*I_ESP_F3z_G2x2y;
    abcd[iGrid*675+55] = 2.0E0*I_ESP_Hx4y_G2x2y_a;
    abcd[iGrid*675+56] = 2.0E0*I_ESP_Hx3yz_G2x2y_a;
    abcd[iGrid*675+57] = 2.0E0*I_ESP_Hx2y2z_G2x2y_a;
    abcd[iGrid*675+58] = 2.0E0*I_ESP_Hxy3z_G2x2y_a;
    abcd[iGrid*675+59] = 2.0E0*I_ESP_Hx4z_G2x2y_a;
    abcd[iGrid*675+60] = 2.0E0*I_ESP_H5x_G2xyz_a-4*I_ESP_F3x_G2xyz;
    abcd[iGrid*675+61] = 2.0E0*I_ESP_H4xy_G2xyz_a-3*I_ESP_F2xy_G2xyz;
    abcd[iGrid*675+62] = 2.0E0*I_ESP_H4xz_G2xyz_a-3*I_ESP_F2xz_G2xyz;
    abcd[iGrid*675+63] = 2.0E0*I_ESP_H3x2y_G2xyz_a-2*I_ESP_Fx2y_G2xyz;
    abcd[iGrid*675+64] = 2.0E0*I_ESP_H3xyz_G2xyz_a-2*I_ESP_Fxyz_G2xyz;
    abcd[iGrid*675+65] = 2.0E0*I_ESP_H3x2z_G2xyz_a-2*I_ESP_Fx2z_G2xyz;
    abcd[iGrid*675+66] = 2.0E0*I_ESP_H2x3y_G2xyz_a-1*I_ESP_F3y_G2xyz;
    abcd[iGrid*675+67] = 2.0E0*I_ESP_H2x2yz_G2xyz_a-1*I_ESP_F2yz_G2xyz;
    abcd[iGrid*675+68] = 2.0E0*I_ESP_H2xy2z_G2xyz_a-1*I_ESP_Fy2z_G2xyz;
    abcd[iGrid*675+69] = 2.0E0*I_ESP_H2x3z_G2xyz_a-1*I_ESP_F3z_G2xyz;
    abcd[iGrid*675+70] = 2.0E0*I_ESP_Hx4y_G2xyz_a;
    abcd[iGrid*675+71] = 2.0E0*I_ESP_Hx3yz_G2xyz_a;
    abcd[iGrid*675+72] = 2.0E0*I_ESP_Hx2y2z_G2xyz_a;
    abcd[iGrid*675+73] = 2.0E0*I_ESP_Hxy3z_G2xyz_a;
    abcd[iGrid*675+74] = 2.0E0*I_ESP_Hx4z_G2xyz_a;
    abcd[iGrid*675+75] = 2.0E0*I_ESP_H5x_G2x2z_a-4*I_ESP_F3x_G2x2z;
    abcd[iGrid*675+76] = 2.0E0*I_ESP_H4xy_G2x2z_a-3*I_ESP_F2xy_G2x2z;
    abcd[iGrid*675+77] = 2.0E0*I_ESP_H4xz_G2x2z_a-3*I_ESP_F2xz_G2x2z;
    abcd[iGrid*675+78] = 2.0E0*I_ESP_H3x2y_G2x2z_a-2*I_ESP_Fx2y_G2x2z;
    abcd[iGrid*675+79] = 2.0E0*I_ESP_H3xyz_G2x2z_a-2*I_ESP_Fxyz_G2x2z;
    abcd[iGrid*675+80] = 2.0E0*I_ESP_H3x2z_G2x2z_a-2*I_ESP_Fx2z_G2x2z;
    abcd[iGrid*675+81] = 2.0E0*I_ESP_H2x3y_G2x2z_a-1*I_ESP_F3y_G2x2z;
    abcd[iGrid*675+82] = 2.0E0*I_ESP_H2x2yz_G2x2z_a-1*I_ESP_F2yz_G2x2z;
    abcd[iGrid*675+83] = 2.0E0*I_ESP_H2xy2z_G2x2z_a-1*I_ESP_Fy2z_G2x2z;
    abcd[iGrid*675+84] = 2.0E0*I_ESP_H2x3z_G2x2z_a-1*I_ESP_F3z_G2x2z;
    abcd[iGrid*675+85] = 2.0E0*I_ESP_Hx4y_G2x2z_a;
    abcd[iGrid*675+86] = 2.0E0*I_ESP_Hx3yz_G2x2z_a;
    abcd[iGrid*675+87] = 2.0E0*I_ESP_Hx2y2z_G2x2z_a;
    abcd[iGrid*675+88] = 2.0E0*I_ESP_Hxy3z_G2x2z_a;
    abcd[iGrid*675+89] = 2.0E0*I_ESP_Hx4z_G2x2z_a;
    abcd[iGrid*675+90] = 2.0E0*I_ESP_H5x_Gx3y_a-4*I_ESP_F3x_Gx3y;
    abcd[iGrid*675+91] = 2.0E0*I_ESP_H4xy_Gx3y_a-3*I_ESP_F2xy_Gx3y;
    abcd[iGrid*675+92] = 2.0E0*I_ESP_H4xz_Gx3y_a-3*I_ESP_F2xz_Gx3y;
    abcd[iGrid*675+93] = 2.0E0*I_ESP_H3x2y_Gx3y_a-2*I_ESP_Fx2y_Gx3y;
    abcd[iGrid*675+94] = 2.0E0*I_ESP_H3xyz_Gx3y_a-2*I_ESP_Fxyz_Gx3y;
    abcd[iGrid*675+95] = 2.0E0*I_ESP_H3x2z_Gx3y_a-2*I_ESP_Fx2z_Gx3y;
    abcd[iGrid*675+96] = 2.0E0*I_ESP_H2x3y_Gx3y_a-1*I_ESP_F3y_Gx3y;
    abcd[iGrid*675+97] = 2.0E0*I_ESP_H2x2yz_Gx3y_a-1*I_ESP_F2yz_Gx3y;
    abcd[iGrid*675+98] = 2.0E0*I_ESP_H2xy2z_Gx3y_a-1*I_ESP_Fy2z_Gx3y;
    abcd[iGrid*675+99] = 2.0E0*I_ESP_H2x3z_Gx3y_a-1*I_ESP_F3z_Gx3y;
    abcd[iGrid*675+100] = 2.0E0*I_ESP_Hx4y_Gx3y_a;
    abcd[iGrid*675+101] = 2.0E0*I_ESP_Hx3yz_Gx3y_a;
    abcd[iGrid*675+102] = 2.0E0*I_ESP_Hx2y2z_Gx3y_a;
    abcd[iGrid*675+103] = 2.0E0*I_ESP_Hxy3z_Gx3y_a;
    abcd[iGrid*675+104] = 2.0E0*I_ESP_Hx4z_Gx3y_a;
    abcd[iGrid*675+105] = 2.0E0*I_ESP_H5x_Gx2yz_a-4*I_ESP_F3x_Gx2yz;
    abcd[iGrid*675+106] = 2.0E0*I_ESP_H4xy_Gx2yz_a-3*I_ESP_F2xy_Gx2yz;
    abcd[iGrid*675+107] = 2.0E0*I_ESP_H4xz_Gx2yz_a-3*I_ESP_F2xz_Gx2yz;
    abcd[iGrid*675+108] = 2.0E0*I_ESP_H3x2y_Gx2yz_a-2*I_ESP_Fx2y_Gx2yz;
    abcd[iGrid*675+109] = 2.0E0*I_ESP_H3xyz_Gx2yz_a-2*I_ESP_Fxyz_Gx2yz;
    abcd[iGrid*675+110] = 2.0E0*I_ESP_H3x2z_Gx2yz_a-2*I_ESP_Fx2z_Gx2yz;
    abcd[iGrid*675+111] = 2.0E0*I_ESP_H2x3y_Gx2yz_a-1*I_ESP_F3y_Gx2yz;
    abcd[iGrid*675+112] = 2.0E0*I_ESP_H2x2yz_Gx2yz_a-1*I_ESP_F2yz_Gx2yz;
    abcd[iGrid*675+113] = 2.0E0*I_ESP_H2xy2z_Gx2yz_a-1*I_ESP_Fy2z_Gx2yz;
    abcd[iGrid*675+114] = 2.0E0*I_ESP_H2x3z_Gx2yz_a-1*I_ESP_F3z_Gx2yz;
    abcd[iGrid*675+115] = 2.0E0*I_ESP_Hx4y_Gx2yz_a;
    abcd[iGrid*675+116] = 2.0E0*I_ESP_Hx3yz_Gx2yz_a;
    abcd[iGrid*675+117] = 2.0E0*I_ESP_Hx2y2z_Gx2yz_a;
    abcd[iGrid*675+118] = 2.0E0*I_ESP_Hxy3z_Gx2yz_a;
    abcd[iGrid*675+119] = 2.0E0*I_ESP_Hx4z_Gx2yz_a;
    abcd[iGrid*675+120] = 2.0E0*I_ESP_H5x_Gxy2z_a-4*I_ESP_F3x_Gxy2z;
    abcd[iGrid*675+121] = 2.0E0*I_ESP_H4xy_Gxy2z_a-3*I_ESP_F2xy_Gxy2z;
    abcd[iGrid*675+122] = 2.0E0*I_ESP_H4xz_Gxy2z_a-3*I_ESP_F2xz_Gxy2z;
    abcd[iGrid*675+123] = 2.0E0*I_ESP_H3x2y_Gxy2z_a-2*I_ESP_Fx2y_Gxy2z;
    abcd[iGrid*675+124] = 2.0E0*I_ESP_H3xyz_Gxy2z_a-2*I_ESP_Fxyz_Gxy2z;
    abcd[iGrid*675+125] = 2.0E0*I_ESP_H3x2z_Gxy2z_a-2*I_ESP_Fx2z_Gxy2z;
    abcd[iGrid*675+126] = 2.0E0*I_ESP_H2x3y_Gxy2z_a-1*I_ESP_F3y_Gxy2z;
    abcd[iGrid*675+127] = 2.0E0*I_ESP_H2x2yz_Gxy2z_a-1*I_ESP_F2yz_Gxy2z;
    abcd[iGrid*675+128] = 2.0E0*I_ESP_H2xy2z_Gxy2z_a-1*I_ESP_Fy2z_Gxy2z;
    abcd[iGrid*675+129] = 2.0E0*I_ESP_H2x3z_Gxy2z_a-1*I_ESP_F3z_Gxy2z;
    abcd[iGrid*675+130] = 2.0E0*I_ESP_Hx4y_Gxy2z_a;
    abcd[iGrid*675+131] = 2.0E0*I_ESP_Hx3yz_Gxy2z_a;
    abcd[iGrid*675+132] = 2.0E0*I_ESP_Hx2y2z_Gxy2z_a;
    abcd[iGrid*675+133] = 2.0E0*I_ESP_Hxy3z_Gxy2z_a;
    abcd[iGrid*675+134] = 2.0E0*I_ESP_Hx4z_Gxy2z_a;
    abcd[iGrid*675+135] = 2.0E0*I_ESP_H5x_Gx3z_a-4*I_ESP_F3x_Gx3z;
    abcd[iGrid*675+136] = 2.0E0*I_ESP_H4xy_Gx3z_a-3*I_ESP_F2xy_Gx3z;
    abcd[iGrid*675+137] = 2.0E0*I_ESP_H4xz_Gx3z_a-3*I_ESP_F2xz_Gx3z;
    abcd[iGrid*675+138] = 2.0E0*I_ESP_H3x2y_Gx3z_a-2*I_ESP_Fx2y_Gx3z;
    abcd[iGrid*675+139] = 2.0E0*I_ESP_H3xyz_Gx3z_a-2*I_ESP_Fxyz_Gx3z;
    abcd[iGrid*675+140] = 2.0E0*I_ESP_H3x2z_Gx3z_a-2*I_ESP_Fx2z_Gx3z;
    abcd[iGrid*675+141] = 2.0E0*I_ESP_H2x3y_Gx3z_a-1*I_ESP_F3y_Gx3z;
    abcd[iGrid*675+142] = 2.0E0*I_ESP_H2x2yz_Gx3z_a-1*I_ESP_F2yz_Gx3z;
    abcd[iGrid*675+143] = 2.0E0*I_ESP_H2xy2z_Gx3z_a-1*I_ESP_Fy2z_Gx3z;
    abcd[iGrid*675+144] = 2.0E0*I_ESP_H2x3z_Gx3z_a-1*I_ESP_F3z_Gx3z;
    abcd[iGrid*675+145] = 2.0E0*I_ESP_Hx4y_Gx3z_a;
    abcd[iGrid*675+146] = 2.0E0*I_ESP_Hx3yz_Gx3z_a;
    abcd[iGrid*675+147] = 2.0E0*I_ESP_Hx2y2z_Gx3z_a;
    abcd[iGrid*675+148] = 2.0E0*I_ESP_Hxy3z_Gx3z_a;
    abcd[iGrid*675+149] = 2.0E0*I_ESP_Hx4z_Gx3z_a;
    abcd[iGrid*675+150] = 2.0E0*I_ESP_H5x_G4y_a-4*I_ESP_F3x_G4y;
    abcd[iGrid*675+151] = 2.0E0*I_ESP_H4xy_G4y_a-3*I_ESP_F2xy_G4y;
    abcd[iGrid*675+152] = 2.0E0*I_ESP_H4xz_G4y_a-3*I_ESP_F2xz_G4y;
    abcd[iGrid*675+153] = 2.0E0*I_ESP_H3x2y_G4y_a-2*I_ESP_Fx2y_G4y;
    abcd[iGrid*675+154] = 2.0E0*I_ESP_H3xyz_G4y_a-2*I_ESP_Fxyz_G4y;
    abcd[iGrid*675+155] = 2.0E0*I_ESP_H3x2z_G4y_a-2*I_ESP_Fx2z_G4y;
    abcd[iGrid*675+156] = 2.0E0*I_ESP_H2x3y_G4y_a-1*I_ESP_F3y_G4y;
    abcd[iGrid*675+157] = 2.0E0*I_ESP_H2x2yz_G4y_a-1*I_ESP_F2yz_G4y;
    abcd[iGrid*675+158] = 2.0E0*I_ESP_H2xy2z_G4y_a-1*I_ESP_Fy2z_G4y;
    abcd[iGrid*675+159] = 2.0E0*I_ESP_H2x3z_G4y_a-1*I_ESP_F3z_G4y;
    abcd[iGrid*675+160] = 2.0E0*I_ESP_Hx4y_G4y_a;
    abcd[iGrid*675+161] = 2.0E0*I_ESP_Hx3yz_G4y_a;
    abcd[iGrid*675+162] = 2.0E0*I_ESP_Hx2y2z_G4y_a;
    abcd[iGrid*675+163] = 2.0E0*I_ESP_Hxy3z_G4y_a;
    abcd[iGrid*675+164] = 2.0E0*I_ESP_Hx4z_G4y_a;
    abcd[iGrid*675+165] = 2.0E0*I_ESP_H5x_G3yz_a-4*I_ESP_F3x_G3yz;
    abcd[iGrid*675+166] = 2.0E0*I_ESP_H4xy_G3yz_a-3*I_ESP_F2xy_G3yz;
    abcd[iGrid*675+167] = 2.0E0*I_ESP_H4xz_G3yz_a-3*I_ESP_F2xz_G3yz;
    abcd[iGrid*675+168] = 2.0E0*I_ESP_H3x2y_G3yz_a-2*I_ESP_Fx2y_G3yz;
    abcd[iGrid*675+169] = 2.0E0*I_ESP_H3xyz_G3yz_a-2*I_ESP_Fxyz_G3yz;
    abcd[iGrid*675+170] = 2.0E0*I_ESP_H3x2z_G3yz_a-2*I_ESP_Fx2z_G3yz;
    abcd[iGrid*675+171] = 2.0E0*I_ESP_H2x3y_G3yz_a-1*I_ESP_F3y_G3yz;
    abcd[iGrid*675+172] = 2.0E0*I_ESP_H2x2yz_G3yz_a-1*I_ESP_F2yz_G3yz;
    abcd[iGrid*675+173] = 2.0E0*I_ESP_H2xy2z_G3yz_a-1*I_ESP_Fy2z_G3yz;
    abcd[iGrid*675+174] = 2.0E0*I_ESP_H2x3z_G3yz_a-1*I_ESP_F3z_G3yz;
    abcd[iGrid*675+175] = 2.0E0*I_ESP_Hx4y_G3yz_a;
    abcd[iGrid*675+176] = 2.0E0*I_ESP_Hx3yz_G3yz_a;
    abcd[iGrid*675+177] = 2.0E0*I_ESP_Hx2y2z_G3yz_a;
    abcd[iGrid*675+178] = 2.0E0*I_ESP_Hxy3z_G3yz_a;
    abcd[iGrid*675+179] = 2.0E0*I_ESP_Hx4z_G3yz_a;
    abcd[iGrid*675+180] = 2.0E0*I_ESP_H5x_G2y2z_a-4*I_ESP_F3x_G2y2z;
    abcd[iGrid*675+181] = 2.0E0*I_ESP_H4xy_G2y2z_a-3*I_ESP_F2xy_G2y2z;
    abcd[iGrid*675+182] = 2.0E0*I_ESP_H4xz_G2y2z_a-3*I_ESP_F2xz_G2y2z;
    abcd[iGrid*675+183] = 2.0E0*I_ESP_H3x2y_G2y2z_a-2*I_ESP_Fx2y_G2y2z;
    abcd[iGrid*675+184] = 2.0E0*I_ESP_H3xyz_G2y2z_a-2*I_ESP_Fxyz_G2y2z;
    abcd[iGrid*675+185] = 2.0E0*I_ESP_H3x2z_G2y2z_a-2*I_ESP_Fx2z_G2y2z;
    abcd[iGrid*675+186] = 2.0E0*I_ESP_H2x3y_G2y2z_a-1*I_ESP_F3y_G2y2z;
    abcd[iGrid*675+187] = 2.0E0*I_ESP_H2x2yz_G2y2z_a-1*I_ESP_F2yz_G2y2z;
    abcd[iGrid*675+188] = 2.0E0*I_ESP_H2xy2z_G2y2z_a-1*I_ESP_Fy2z_G2y2z;
    abcd[iGrid*675+189] = 2.0E0*I_ESP_H2x3z_G2y2z_a-1*I_ESP_F3z_G2y2z;
    abcd[iGrid*675+190] = 2.0E0*I_ESP_Hx4y_G2y2z_a;
    abcd[iGrid*675+191] = 2.0E0*I_ESP_Hx3yz_G2y2z_a;
    abcd[iGrid*675+192] = 2.0E0*I_ESP_Hx2y2z_G2y2z_a;
    abcd[iGrid*675+193] = 2.0E0*I_ESP_Hxy3z_G2y2z_a;
    abcd[iGrid*675+194] = 2.0E0*I_ESP_Hx4z_G2y2z_a;
    abcd[iGrid*675+195] = 2.0E0*I_ESP_H5x_Gy3z_a-4*I_ESP_F3x_Gy3z;
    abcd[iGrid*675+196] = 2.0E0*I_ESP_H4xy_Gy3z_a-3*I_ESP_F2xy_Gy3z;
    abcd[iGrid*675+197] = 2.0E0*I_ESP_H4xz_Gy3z_a-3*I_ESP_F2xz_Gy3z;
    abcd[iGrid*675+198] = 2.0E0*I_ESP_H3x2y_Gy3z_a-2*I_ESP_Fx2y_Gy3z;
    abcd[iGrid*675+199] = 2.0E0*I_ESP_H3xyz_Gy3z_a-2*I_ESP_Fxyz_Gy3z;
    abcd[iGrid*675+200] = 2.0E0*I_ESP_H3x2z_Gy3z_a-2*I_ESP_Fx2z_Gy3z;
    abcd[iGrid*675+201] = 2.0E0*I_ESP_H2x3y_Gy3z_a-1*I_ESP_F3y_Gy3z;
    abcd[iGrid*675+202] = 2.0E0*I_ESP_H2x2yz_Gy3z_a-1*I_ESP_F2yz_Gy3z;
    abcd[iGrid*675+203] = 2.0E0*I_ESP_H2xy2z_Gy3z_a-1*I_ESP_Fy2z_Gy3z;
    abcd[iGrid*675+204] = 2.0E0*I_ESP_H2x3z_Gy3z_a-1*I_ESP_F3z_Gy3z;
    abcd[iGrid*675+205] = 2.0E0*I_ESP_Hx4y_Gy3z_a;
    abcd[iGrid*675+206] = 2.0E0*I_ESP_Hx3yz_Gy3z_a;
    abcd[iGrid*675+207] = 2.0E0*I_ESP_Hx2y2z_Gy3z_a;
    abcd[iGrid*675+208] = 2.0E0*I_ESP_Hxy3z_Gy3z_a;
    abcd[iGrid*675+209] = 2.0E0*I_ESP_Hx4z_Gy3z_a;
    abcd[iGrid*675+210] = 2.0E0*I_ESP_H5x_G4z_a-4*I_ESP_F3x_G4z;
    abcd[iGrid*675+211] = 2.0E0*I_ESP_H4xy_G4z_a-3*I_ESP_F2xy_G4z;
    abcd[iGrid*675+212] = 2.0E0*I_ESP_H4xz_G4z_a-3*I_ESP_F2xz_G4z;
    abcd[iGrid*675+213] = 2.0E0*I_ESP_H3x2y_G4z_a-2*I_ESP_Fx2y_G4z;
    abcd[iGrid*675+214] = 2.0E0*I_ESP_H3xyz_G4z_a-2*I_ESP_Fxyz_G4z;
    abcd[iGrid*675+215] = 2.0E0*I_ESP_H3x2z_G4z_a-2*I_ESP_Fx2z_G4z;
    abcd[iGrid*675+216] = 2.0E0*I_ESP_H2x3y_G4z_a-1*I_ESP_F3y_G4z;
    abcd[iGrid*675+217] = 2.0E0*I_ESP_H2x2yz_G4z_a-1*I_ESP_F2yz_G4z;
    abcd[iGrid*675+218] = 2.0E0*I_ESP_H2xy2z_G4z_a-1*I_ESP_Fy2z_G4z;
    abcd[iGrid*675+219] = 2.0E0*I_ESP_H2x3z_G4z_a-1*I_ESP_F3z_G4z;
    abcd[iGrid*675+220] = 2.0E0*I_ESP_Hx4y_G4z_a;
    abcd[iGrid*675+221] = 2.0E0*I_ESP_Hx3yz_G4z_a;
    abcd[iGrid*675+222] = 2.0E0*I_ESP_Hx2y2z_G4z_a;
    abcd[iGrid*675+223] = 2.0E0*I_ESP_Hxy3z_G4z_a;
    abcd[iGrid*675+224] = 2.0E0*I_ESP_Hx4z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_F_G
     ************************************************************/
    abcd[iGrid*675+225] = 2.0E0*I_ESP_H4xy_G4x_a;
    abcd[iGrid*675+226] = 2.0E0*I_ESP_H3x2y_G4x_a-1*I_ESP_F3x_G4x;
    abcd[iGrid*675+227] = 2.0E0*I_ESP_H3xyz_G4x_a;
    abcd[iGrid*675+228] = 2.0E0*I_ESP_H2x3y_G4x_a-2*I_ESP_F2xy_G4x;
    abcd[iGrid*675+229] = 2.0E0*I_ESP_H2x2yz_G4x_a-1*I_ESP_F2xz_G4x;
    abcd[iGrid*675+230] = 2.0E0*I_ESP_H2xy2z_G4x_a;
    abcd[iGrid*675+231] = 2.0E0*I_ESP_Hx4y_G4x_a-3*I_ESP_Fx2y_G4x;
    abcd[iGrid*675+232] = 2.0E0*I_ESP_Hx3yz_G4x_a-2*I_ESP_Fxyz_G4x;
    abcd[iGrid*675+233] = 2.0E0*I_ESP_Hx2y2z_G4x_a-1*I_ESP_Fx2z_G4x;
    abcd[iGrid*675+234] = 2.0E0*I_ESP_Hxy3z_G4x_a;
    abcd[iGrid*675+235] = 2.0E0*I_ESP_H5y_G4x_a-4*I_ESP_F3y_G4x;
    abcd[iGrid*675+236] = 2.0E0*I_ESP_H4yz_G4x_a-3*I_ESP_F2yz_G4x;
    abcd[iGrid*675+237] = 2.0E0*I_ESP_H3y2z_G4x_a-2*I_ESP_Fy2z_G4x;
    abcd[iGrid*675+238] = 2.0E0*I_ESP_H2y3z_G4x_a-1*I_ESP_F3z_G4x;
    abcd[iGrid*675+239] = 2.0E0*I_ESP_Hy4z_G4x_a;
    abcd[iGrid*675+240] = 2.0E0*I_ESP_H4xy_G3xy_a;
    abcd[iGrid*675+241] = 2.0E0*I_ESP_H3x2y_G3xy_a-1*I_ESP_F3x_G3xy;
    abcd[iGrid*675+242] = 2.0E0*I_ESP_H3xyz_G3xy_a;
    abcd[iGrid*675+243] = 2.0E0*I_ESP_H2x3y_G3xy_a-2*I_ESP_F2xy_G3xy;
    abcd[iGrid*675+244] = 2.0E0*I_ESP_H2x2yz_G3xy_a-1*I_ESP_F2xz_G3xy;
    abcd[iGrid*675+245] = 2.0E0*I_ESP_H2xy2z_G3xy_a;
    abcd[iGrid*675+246] = 2.0E0*I_ESP_Hx4y_G3xy_a-3*I_ESP_Fx2y_G3xy;
    abcd[iGrid*675+247] = 2.0E0*I_ESP_Hx3yz_G3xy_a-2*I_ESP_Fxyz_G3xy;
    abcd[iGrid*675+248] = 2.0E0*I_ESP_Hx2y2z_G3xy_a-1*I_ESP_Fx2z_G3xy;
    abcd[iGrid*675+249] = 2.0E0*I_ESP_Hxy3z_G3xy_a;
    abcd[iGrid*675+250] = 2.0E0*I_ESP_H5y_G3xy_a-4*I_ESP_F3y_G3xy;
    abcd[iGrid*675+251] = 2.0E0*I_ESP_H4yz_G3xy_a-3*I_ESP_F2yz_G3xy;
    abcd[iGrid*675+252] = 2.0E0*I_ESP_H3y2z_G3xy_a-2*I_ESP_Fy2z_G3xy;
    abcd[iGrid*675+253] = 2.0E0*I_ESP_H2y3z_G3xy_a-1*I_ESP_F3z_G3xy;
    abcd[iGrid*675+254] = 2.0E0*I_ESP_Hy4z_G3xy_a;
    abcd[iGrid*675+255] = 2.0E0*I_ESP_H4xy_G3xz_a;
    abcd[iGrid*675+256] = 2.0E0*I_ESP_H3x2y_G3xz_a-1*I_ESP_F3x_G3xz;
    abcd[iGrid*675+257] = 2.0E0*I_ESP_H3xyz_G3xz_a;
    abcd[iGrid*675+258] = 2.0E0*I_ESP_H2x3y_G3xz_a-2*I_ESP_F2xy_G3xz;
    abcd[iGrid*675+259] = 2.0E0*I_ESP_H2x2yz_G3xz_a-1*I_ESP_F2xz_G3xz;
    abcd[iGrid*675+260] = 2.0E0*I_ESP_H2xy2z_G3xz_a;
    abcd[iGrid*675+261] = 2.0E0*I_ESP_Hx4y_G3xz_a-3*I_ESP_Fx2y_G3xz;
    abcd[iGrid*675+262] = 2.0E0*I_ESP_Hx3yz_G3xz_a-2*I_ESP_Fxyz_G3xz;
    abcd[iGrid*675+263] = 2.0E0*I_ESP_Hx2y2z_G3xz_a-1*I_ESP_Fx2z_G3xz;
    abcd[iGrid*675+264] = 2.0E0*I_ESP_Hxy3z_G3xz_a;
    abcd[iGrid*675+265] = 2.0E0*I_ESP_H5y_G3xz_a-4*I_ESP_F3y_G3xz;
    abcd[iGrid*675+266] = 2.0E0*I_ESP_H4yz_G3xz_a-3*I_ESP_F2yz_G3xz;
    abcd[iGrid*675+267] = 2.0E0*I_ESP_H3y2z_G3xz_a-2*I_ESP_Fy2z_G3xz;
    abcd[iGrid*675+268] = 2.0E0*I_ESP_H2y3z_G3xz_a-1*I_ESP_F3z_G3xz;
    abcd[iGrid*675+269] = 2.0E0*I_ESP_Hy4z_G3xz_a;
    abcd[iGrid*675+270] = 2.0E0*I_ESP_H4xy_G2x2y_a;
    abcd[iGrid*675+271] = 2.0E0*I_ESP_H3x2y_G2x2y_a-1*I_ESP_F3x_G2x2y;
    abcd[iGrid*675+272] = 2.0E0*I_ESP_H3xyz_G2x2y_a;
    abcd[iGrid*675+273] = 2.0E0*I_ESP_H2x3y_G2x2y_a-2*I_ESP_F2xy_G2x2y;
    abcd[iGrid*675+274] = 2.0E0*I_ESP_H2x2yz_G2x2y_a-1*I_ESP_F2xz_G2x2y;
    abcd[iGrid*675+275] = 2.0E0*I_ESP_H2xy2z_G2x2y_a;
    abcd[iGrid*675+276] = 2.0E0*I_ESP_Hx4y_G2x2y_a-3*I_ESP_Fx2y_G2x2y;
    abcd[iGrid*675+277] = 2.0E0*I_ESP_Hx3yz_G2x2y_a-2*I_ESP_Fxyz_G2x2y;
    abcd[iGrid*675+278] = 2.0E0*I_ESP_Hx2y2z_G2x2y_a-1*I_ESP_Fx2z_G2x2y;
    abcd[iGrid*675+279] = 2.0E0*I_ESP_Hxy3z_G2x2y_a;
    abcd[iGrid*675+280] = 2.0E0*I_ESP_H5y_G2x2y_a-4*I_ESP_F3y_G2x2y;
    abcd[iGrid*675+281] = 2.0E0*I_ESP_H4yz_G2x2y_a-3*I_ESP_F2yz_G2x2y;
    abcd[iGrid*675+282] = 2.0E0*I_ESP_H3y2z_G2x2y_a-2*I_ESP_Fy2z_G2x2y;
    abcd[iGrid*675+283] = 2.0E0*I_ESP_H2y3z_G2x2y_a-1*I_ESP_F3z_G2x2y;
    abcd[iGrid*675+284] = 2.0E0*I_ESP_Hy4z_G2x2y_a;
    abcd[iGrid*675+285] = 2.0E0*I_ESP_H4xy_G2xyz_a;
    abcd[iGrid*675+286] = 2.0E0*I_ESP_H3x2y_G2xyz_a-1*I_ESP_F3x_G2xyz;
    abcd[iGrid*675+287] = 2.0E0*I_ESP_H3xyz_G2xyz_a;
    abcd[iGrid*675+288] = 2.0E0*I_ESP_H2x3y_G2xyz_a-2*I_ESP_F2xy_G2xyz;
    abcd[iGrid*675+289] = 2.0E0*I_ESP_H2x2yz_G2xyz_a-1*I_ESP_F2xz_G2xyz;
    abcd[iGrid*675+290] = 2.0E0*I_ESP_H2xy2z_G2xyz_a;
    abcd[iGrid*675+291] = 2.0E0*I_ESP_Hx4y_G2xyz_a-3*I_ESP_Fx2y_G2xyz;
    abcd[iGrid*675+292] = 2.0E0*I_ESP_Hx3yz_G2xyz_a-2*I_ESP_Fxyz_G2xyz;
    abcd[iGrid*675+293] = 2.0E0*I_ESP_Hx2y2z_G2xyz_a-1*I_ESP_Fx2z_G2xyz;
    abcd[iGrid*675+294] = 2.0E0*I_ESP_Hxy3z_G2xyz_a;
    abcd[iGrid*675+295] = 2.0E0*I_ESP_H5y_G2xyz_a-4*I_ESP_F3y_G2xyz;
    abcd[iGrid*675+296] = 2.0E0*I_ESP_H4yz_G2xyz_a-3*I_ESP_F2yz_G2xyz;
    abcd[iGrid*675+297] = 2.0E0*I_ESP_H3y2z_G2xyz_a-2*I_ESP_Fy2z_G2xyz;
    abcd[iGrid*675+298] = 2.0E0*I_ESP_H2y3z_G2xyz_a-1*I_ESP_F3z_G2xyz;
    abcd[iGrid*675+299] = 2.0E0*I_ESP_Hy4z_G2xyz_a;
    abcd[iGrid*675+300] = 2.0E0*I_ESP_H4xy_G2x2z_a;
    abcd[iGrid*675+301] = 2.0E0*I_ESP_H3x2y_G2x2z_a-1*I_ESP_F3x_G2x2z;
    abcd[iGrid*675+302] = 2.0E0*I_ESP_H3xyz_G2x2z_a;
    abcd[iGrid*675+303] = 2.0E0*I_ESP_H2x3y_G2x2z_a-2*I_ESP_F2xy_G2x2z;
    abcd[iGrid*675+304] = 2.0E0*I_ESP_H2x2yz_G2x2z_a-1*I_ESP_F2xz_G2x2z;
    abcd[iGrid*675+305] = 2.0E0*I_ESP_H2xy2z_G2x2z_a;
    abcd[iGrid*675+306] = 2.0E0*I_ESP_Hx4y_G2x2z_a-3*I_ESP_Fx2y_G2x2z;
    abcd[iGrid*675+307] = 2.0E0*I_ESP_Hx3yz_G2x2z_a-2*I_ESP_Fxyz_G2x2z;
    abcd[iGrid*675+308] = 2.0E0*I_ESP_Hx2y2z_G2x2z_a-1*I_ESP_Fx2z_G2x2z;
    abcd[iGrid*675+309] = 2.0E0*I_ESP_Hxy3z_G2x2z_a;
    abcd[iGrid*675+310] = 2.0E0*I_ESP_H5y_G2x2z_a-4*I_ESP_F3y_G2x2z;
    abcd[iGrid*675+311] = 2.0E0*I_ESP_H4yz_G2x2z_a-3*I_ESP_F2yz_G2x2z;
    abcd[iGrid*675+312] = 2.0E0*I_ESP_H3y2z_G2x2z_a-2*I_ESP_Fy2z_G2x2z;
    abcd[iGrid*675+313] = 2.0E0*I_ESP_H2y3z_G2x2z_a-1*I_ESP_F3z_G2x2z;
    abcd[iGrid*675+314] = 2.0E0*I_ESP_Hy4z_G2x2z_a;
    abcd[iGrid*675+315] = 2.0E0*I_ESP_H4xy_Gx3y_a;
    abcd[iGrid*675+316] = 2.0E0*I_ESP_H3x2y_Gx3y_a-1*I_ESP_F3x_Gx3y;
    abcd[iGrid*675+317] = 2.0E0*I_ESP_H3xyz_Gx3y_a;
    abcd[iGrid*675+318] = 2.0E0*I_ESP_H2x3y_Gx3y_a-2*I_ESP_F2xy_Gx3y;
    abcd[iGrid*675+319] = 2.0E0*I_ESP_H2x2yz_Gx3y_a-1*I_ESP_F2xz_Gx3y;
    abcd[iGrid*675+320] = 2.0E0*I_ESP_H2xy2z_Gx3y_a;
    abcd[iGrid*675+321] = 2.0E0*I_ESP_Hx4y_Gx3y_a-3*I_ESP_Fx2y_Gx3y;
    abcd[iGrid*675+322] = 2.0E0*I_ESP_Hx3yz_Gx3y_a-2*I_ESP_Fxyz_Gx3y;
    abcd[iGrid*675+323] = 2.0E0*I_ESP_Hx2y2z_Gx3y_a-1*I_ESP_Fx2z_Gx3y;
    abcd[iGrid*675+324] = 2.0E0*I_ESP_Hxy3z_Gx3y_a;
    abcd[iGrid*675+325] = 2.0E0*I_ESP_H5y_Gx3y_a-4*I_ESP_F3y_Gx3y;
    abcd[iGrid*675+326] = 2.0E0*I_ESP_H4yz_Gx3y_a-3*I_ESP_F2yz_Gx3y;
    abcd[iGrid*675+327] = 2.0E0*I_ESP_H3y2z_Gx3y_a-2*I_ESP_Fy2z_Gx3y;
    abcd[iGrid*675+328] = 2.0E0*I_ESP_H2y3z_Gx3y_a-1*I_ESP_F3z_Gx3y;
    abcd[iGrid*675+329] = 2.0E0*I_ESP_Hy4z_Gx3y_a;
    abcd[iGrid*675+330] = 2.0E0*I_ESP_H4xy_Gx2yz_a;
    abcd[iGrid*675+331] = 2.0E0*I_ESP_H3x2y_Gx2yz_a-1*I_ESP_F3x_Gx2yz;
    abcd[iGrid*675+332] = 2.0E0*I_ESP_H3xyz_Gx2yz_a;
    abcd[iGrid*675+333] = 2.0E0*I_ESP_H2x3y_Gx2yz_a-2*I_ESP_F2xy_Gx2yz;
    abcd[iGrid*675+334] = 2.0E0*I_ESP_H2x2yz_Gx2yz_a-1*I_ESP_F2xz_Gx2yz;
    abcd[iGrid*675+335] = 2.0E0*I_ESP_H2xy2z_Gx2yz_a;
    abcd[iGrid*675+336] = 2.0E0*I_ESP_Hx4y_Gx2yz_a-3*I_ESP_Fx2y_Gx2yz;
    abcd[iGrid*675+337] = 2.0E0*I_ESP_Hx3yz_Gx2yz_a-2*I_ESP_Fxyz_Gx2yz;
    abcd[iGrid*675+338] = 2.0E0*I_ESP_Hx2y2z_Gx2yz_a-1*I_ESP_Fx2z_Gx2yz;
    abcd[iGrid*675+339] = 2.0E0*I_ESP_Hxy3z_Gx2yz_a;
    abcd[iGrid*675+340] = 2.0E0*I_ESP_H5y_Gx2yz_a-4*I_ESP_F3y_Gx2yz;
    abcd[iGrid*675+341] = 2.0E0*I_ESP_H4yz_Gx2yz_a-3*I_ESP_F2yz_Gx2yz;
    abcd[iGrid*675+342] = 2.0E0*I_ESP_H3y2z_Gx2yz_a-2*I_ESP_Fy2z_Gx2yz;
    abcd[iGrid*675+343] = 2.0E0*I_ESP_H2y3z_Gx2yz_a-1*I_ESP_F3z_Gx2yz;
    abcd[iGrid*675+344] = 2.0E0*I_ESP_Hy4z_Gx2yz_a;
    abcd[iGrid*675+345] = 2.0E0*I_ESP_H4xy_Gxy2z_a;
    abcd[iGrid*675+346] = 2.0E0*I_ESP_H3x2y_Gxy2z_a-1*I_ESP_F3x_Gxy2z;
    abcd[iGrid*675+347] = 2.0E0*I_ESP_H3xyz_Gxy2z_a;
    abcd[iGrid*675+348] = 2.0E0*I_ESP_H2x3y_Gxy2z_a-2*I_ESP_F2xy_Gxy2z;
    abcd[iGrid*675+349] = 2.0E0*I_ESP_H2x2yz_Gxy2z_a-1*I_ESP_F2xz_Gxy2z;
    abcd[iGrid*675+350] = 2.0E0*I_ESP_H2xy2z_Gxy2z_a;
    abcd[iGrid*675+351] = 2.0E0*I_ESP_Hx4y_Gxy2z_a-3*I_ESP_Fx2y_Gxy2z;
    abcd[iGrid*675+352] = 2.0E0*I_ESP_Hx3yz_Gxy2z_a-2*I_ESP_Fxyz_Gxy2z;
    abcd[iGrid*675+353] = 2.0E0*I_ESP_Hx2y2z_Gxy2z_a-1*I_ESP_Fx2z_Gxy2z;
    abcd[iGrid*675+354] = 2.0E0*I_ESP_Hxy3z_Gxy2z_a;
    abcd[iGrid*675+355] = 2.0E0*I_ESP_H5y_Gxy2z_a-4*I_ESP_F3y_Gxy2z;
    abcd[iGrid*675+356] = 2.0E0*I_ESP_H4yz_Gxy2z_a-3*I_ESP_F2yz_Gxy2z;
    abcd[iGrid*675+357] = 2.0E0*I_ESP_H3y2z_Gxy2z_a-2*I_ESP_Fy2z_Gxy2z;
    abcd[iGrid*675+358] = 2.0E0*I_ESP_H2y3z_Gxy2z_a-1*I_ESP_F3z_Gxy2z;
    abcd[iGrid*675+359] = 2.0E0*I_ESP_Hy4z_Gxy2z_a;
    abcd[iGrid*675+360] = 2.0E0*I_ESP_H4xy_Gx3z_a;
    abcd[iGrid*675+361] = 2.0E0*I_ESP_H3x2y_Gx3z_a-1*I_ESP_F3x_Gx3z;
    abcd[iGrid*675+362] = 2.0E0*I_ESP_H3xyz_Gx3z_a;
    abcd[iGrid*675+363] = 2.0E0*I_ESP_H2x3y_Gx3z_a-2*I_ESP_F2xy_Gx3z;
    abcd[iGrid*675+364] = 2.0E0*I_ESP_H2x2yz_Gx3z_a-1*I_ESP_F2xz_Gx3z;
    abcd[iGrid*675+365] = 2.0E0*I_ESP_H2xy2z_Gx3z_a;
    abcd[iGrid*675+366] = 2.0E0*I_ESP_Hx4y_Gx3z_a-3*I_ESP_Fx2y_Gx3z;
    abcd[iGrid*675+367] = 2.0E0*I_ESP_Hx3yz_Gx3z_a-2*I_ESP_Fxyz_Gx3z;
    abcd[iGrid*675+368] = 2.0E0*I_ESP_Hx2y2z_Gx3z_a-1*I_ESP_Fx2z_Gx3z;
    abcd[iGrid*675+369] = 2.0E0*I_ESP_Hxy3z_Gx3z_a;
    abcd[iGrid*675+370] = 2.0E0*I_ESP_H5y_Gx3z_a-4*I_ESP_F3y_Gx3z;
    abcd[iGrid*675+371] = 2.0E0*I_ESP_H4yz_Gx3z_a-3*I_ESP_F2yz_Gx3z;
    abcd[iGrid*675+372] = 2.0E0*I_ESP_H3y2z_Gx3z_a-2*I_ESP_Fy2z_Gx3z;
    abcd[iGrid*675+373] = 2.0E0*I_ESP_H2y3z_Gx3z_a-1*I_ESP_F3z_Gx3z;
    abcd[iGrid*675+374] = 2.0E0*I_ESP_Hy4z_Gx3z_a;
    abcd[iGrid*675+375] = 2.0E0*I_ESP_H4xy_G4y_a;
    abcd[iGrid*675+376] = 2.0E0*I_ESP_H3x2y_G4y_a-1*I_ESP_F3x_G4y;
    abcd[iGrid*675+377] = 2.0E0*I_ESP_H3xyz_G4y_a;
    abcd[iGrid*675+378] = 2.0E0*I_ESP_H2x3y_G4y_a-2*I_ESP_F2xy_G4y;
    abcd[iGrid*675+379] = 2.0E0*I_ESP_H2x2yz_G4y_a-1*I_ESP_F2xz_G4y;
    abcd[iGrid*675+380] = 2.0E0*I_ESP_H2xy2z_G4y_a;
    abcd[iGrid*675+381] = 2.0E0*I_ESP_Hx4y_G4y_a-3*I_ESP_Fx2y_G4y;
    abcd[iGrid*675+382] = 2.0E0*I_ESP_Hx3yz_G4y_a-2*I_ESP_Fxyz_G4y;
    abcd[iGrid*675+383] = 2.0E0*I_ESP_Hx2y2z_G4y_a-1*I_ESP_Fx2z_G4y;
    abcd[iGrid*675+384] = 2.0E0*I_ESP_Hxy3z_G4y_a;
    abcd[iGrid*675+385] = 2.0E0*I_ESP_H5y_G4y_a-4*I_ESP_F3y_G4y;
    abcd[iGrid*675+386] = 2.0E0*I_ESP_H4yz_G4y_a-3*I_ESP_F2yz_G4y;
    abcd[iGrid*675+387] = 2.0E0*I_ESP_H3y2z_G4y_a-2*I_ESP_Fy2z_G4y;
    abcd[iGrid*675+388] = 2.0E0*I_ESP_H2y3z_G4y_a-1*I_ESP_F3z_G4y;
    abcd[iGrid*675+389] = 2.0E0*I_ESP_Hy4z_G4y_a;
    abcd[iGrid*675+390] = 2.0E0*I_ESP_H4xy_G3yz_a;
    abcd[iGrid*675+391] = 2.0E0*I_ESP_H3x2y_G3yz_a-1*I_ESP_F3x_G3yz;
    abcd[iGrid*675+392] = 2.0E0*I_ESP_H3xyz_G3yz_a;
    abcd[iGrid*675+393] = 2.0E0*I_ESP_H2x3y_G3yz_a-2*I_ESP_F2xy_G3yz;
    abcd[iGrid*675+394] = 2.0E0*I_ESP_H2x2yz_G3yz_a-1*I_ESP_F2xz_G3yz;
    abcd[iGrid*675+395] = 2.0E0*I_ESP_H2xy2z_G3yz_a;
    abcd[iGrid*675+396] = 2.0E0*I_ESP_Hx4y_G3yz_a-3*I_ESP_Fx2y_G3yz;
    abcd[iGrid*675+397] = 2.0E0*I_ESP_Hx3yz_G3yz_a-2*I_ESP_Fxyz_G3yz;
    abcd[iGrid*675+398] = 2.0E0*I_ESP_Hx2y2z_G3yz_a-1*I_ESP_Fx2z_G3yz;
    abcd[iGrid*675+399] = 2.0E0*I_ESP_Hxy3z_G3yz_a;
    abcd[iGrid*675+400] = 2.0E0*I_ESP_H5y_G3yz_a-4*I_ESP_F3y_G3yz;
    abcd[iGrid*675+401] = 2.0E0*I_ESP_H4yz_G3yz_a-3*I_ESP_F2yz_G3yz;
    abcd[iGrid*675+402] = 2.0E0*I_ESP_H3y2z_G3yz_a-2*I_ESP_Fy2z_G3yz;
    abcd[iGrid*675+403] = 2.0E0*I_ESP_H2y3z_G3yz_a-1*I_ESP_F3z_G3yz;
    abcd[iGrid*675+404] = 2.0E0*I_ESP_Hy4z_G3yz_a;
    abcd[iGrid*675+405] = 2.0E0*I_ESP_H4xy_G2y2z_a;
    abcd[iGrid*675+406] = 2.0E0*I_ESP_H3x2y_G2y2z_a-1*I_ESP_F3x_G2y2z;
    abcd[iGrid*675+407] = 2.0E0*I_ESP_H3xyz_G2y2z_a;
    abcd[iGrid*675+408] = 2.0E0*I_ESP_H2x3y_G2y2z_a-2*I_ESP_F2xy_G2y2z;
    abcd[iGrid*675+409] = 2.0E0*I_ESP_H2x2yz_G2y2z_a-1*I_ESP_F2xz_G2y2z;
    abcd[iGrid*675+410] = 2.0E0*I_ESP_H2xy2z_G2y2z_a;
    abcd[iGrid*675+411] = 2.0E0*I_ESP_Hx4y_G2y2z_a-3*I_ESP_Fx2y_G2y2z;
    abcd[iGrid*675+412] = 2.0E0*I_ESP_Hx3yz_G2y2z_a-2*I_ESP_Fxyz_G2y2z;
    abcd[iGrid*675+413] = 2.0E0*I_ESP_Hx2y2z_G2y2z_a-1*I_ESP_Fx2z_G2y2z;
    abcd[iGrid*675+414] = 2.0E0*I_ESP_Hxy3z_G2y2z_a;
    abcd[iGrid*675+415] = 2.0E0*I_ESP_H5y_G2y2z_a-4*I_ESP_F3y_G2y2z;
    abcd[iGrid*675+416] = 2.0E0*I_ESP_H4yz_G2y2z_a-3*I_ESP_F2yz_G2y2z;
    abcd[iGrid*675+417] = 2.0E0*I_ESP_H3y2z_G2y2z_a-2*I_ESP_Fy2z_G2y2z;
    abcd[iGrid*675+418] = 2.0E0*I_ESP_H2y3z_G2y2z_a-1*I_ESP_F3z_G2y2z;
    abcd[iGrid*675+419] = 2.0E0*I_ESP_Hy4z_G2y2z_a;
    abcd[iGrid*675+420] = 2.0E0*I_ESP_H4xy_Gy3z_a;
    abcd[iGrid*675+421] = 2.0E0*I_ESP_H3x2y_Gy3z_a-1*I_ESP_F3x_Gy3z;
    abcd[iGrid*675+422] = 2.0E0*I_ESP_H3xyz_Gy3z_a;
    abcd[iGrid*675+423] = 2.0E0*I_ESP_H2x3y_Gy3z_a-2*I_ESP_F2xy_Gy3z;
    abcd[iGrid*675+424] = 2.0E0*I_ESP_H2x2yz_Gy3z_a-1*I_ESP_F2xz_Gy3z;
    abcd[iGrid*675+425] = 2.0E0*I_ESP_H2xy2z_Gy3z_a;
    abcd[iGrid*675+426] = 2.0E0*I_ESP_Hx4y_Gy3z_a-3*I_ESP_Fx2y_Gy3z;
    abcd[iGrid*675+427] = 2.0E0*I_ESP_Hx3yz_Gy3z_a-2*I_ESP_Fxyz_Gy3z;
    abcd[iGrid*675+428] = 2.0E0*I_ESP_Hx2y2z_Gy3z_a-1*I_ESP_Fx2z_Gy3z;
    abcd[iGrid*675+429] = 2.0E0*I_ESP_Hxy3z_Gy3z_a;
    abcd[iGrid*675+430] = 2.0E0*I_ESP_H5y_Gy3z_a-4*I_ESP_F3y_Gy3z;
    abcd[iGrid*675+431] = 2.0E0*I_ESP_H4yz_Gy3z_a-3*I_ESP_F2yz_Gy3z;
    abcd[iGrid*675+432] = 2.0E0*I_ESP_H3y2z_Gy3z_a-2*I_ESP_Fy2z_Gy3z;
    abcd[iGrid*675+433] = 2.0E0*I_ESP_H2y3z_Gy3z_a-1*I_ESP_F3z_Gy3z;
    abcd[iGrid*675+434] = 2.0E0*I_ESP_Hy4z_Gy3z_a;
    abcd[iGrid*675+435] = 2.0E0*I_ESP_H4xy_G4z_a;
    abcd[iGrid*675+436] = 2.0E0*I_ESP_H3x2y_G4z_a-1*I_ESP_F3x_G4z;
    abcd[iGrid*675+437] = 2.0E0*I_ESP_H3xyz_G4z_a;
    abcd[iGrid*675+438] = 2.0E0*I_ESP_H2x3y_G4z_a-2*I_ESP_F2xy_G4z;
    abcd[iGrid*675+439] = 2.0E0*I_ESP_H2x2yz_G4z_a-1*I_ESP_F2xz_G4z;
    abcd[iGrid*675+440] = 2.0E0*I_ESP_H2xy2z_G4z_a;
    abcd[iGrid*675+441] = 2.0E0*I_ESP_Hx4y_G4z_a-3*I_ESP_Fx2y_G4z;
    abcd[iGrid*675+442] = 2.0E0*I_ESP_Hx3yz_G4z_a-2*I_ESP_Fxyz_G4z;
    abcd[iGrid*675+443] = 2.0E0*I_ESP_Hx2y2z_G4z_a-1*I_ESP_Fx2z_G4z;
    abcd[iGrid*675+444] = 2.0E0*I_ESP_Hxy3z_G4z_a;
    abcd[iGrid*675+445] = 2.0E0*I_ESP_H5y_G4z_a-4*I_ESP_F3y_G4z;
    abcd[iGrid*675+446] = 2.0E0*I_ESP_H4yz_G4z_a-3*I_ESP_F2yz_G4z;
    abcd[iGrid*675+447] = 2.0E0*I_ESP_H3y2z_G4z_a-2*I_ESP_Fy2z_G4z;
    abcd[iGrid*675+448] = 2.0E0*I_ESP_H2y3z_G4z_a-1*I_ESP_F3z_G4z;
    abcd[iGrid*675+449] = 2.0E0*I_ESP_Hy4z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_F_G
     ************************************************************/
    abcd[iGrid*675+450] = 2.0E0*I_ESP_H4xz_G4x_a;
    abcd[iGrid*675+451] = 2.0E0*I_ESP_H3xyz_G4x_a;
    abcd[iGrid*675+452] = 2.0E0*I_ESP_H3x2z_G4x_a-1*I_ESP_F3x_G4x;
    abcd[iGrid*675+453] = 2.0E0*I_ESP_H2x2yz_G4x_a;
    abcd[iGrid*675+454] = 2.0E0*I_ESP_H2xy2z_G4x_a-1*I_ESP_F2xy_G4x;
    abcd[iGrid*675+455] = 2.0E0*I_ESP_H2x3z_G4x_a-2*I_ESP_F2xz_G4x;
    abcd[iGrid*675+456] = 2.0E0*I_ESP_Hx3yz_G4x_a;
    abcd[iGrid*675+457] = 2.0E0*I_ESP_Hx2y2z_G4x_a-1*I_ESP_Fx2y_G4x;
    abcd[iGrid*675+458] = 2.0E0*I_ESP_Hxy3z_G4x_a-2*I_ESP_Fxyz_G4x;
    abcd[iGrid*675+459] = 2.0E0*I_ESP_Hx4z_G4x_a-3*I_ESP_Fx2z_G4x;
    abcd[iGrid*675+460] = 2.0E0*I_ESP_H4yz_G4x_a;
    abcd[iGrid*675+461] = 2.0E0*I_ESP_H3y2z_G4x_a-1*I_ESP_F3y_G4x;
    abcd[iGrid*675+462] = 2.0E0*I_ESP_H2y3z_G4x_a-2*I_ESP_F2yz_G4x;
    abcd[iGrid*675+463] = 2.0E0*I_ESP_Hy4z_G4x_a-3*I_ESP_Fy2z_G4x;
    abcd[iGrid*675+464] = 2.0E0*I_ESP_H5z_G4x_a-4*I_ESP_F3z_G4x;
    abcd[iGrid*675+465] = 2.0E0*I_ESP_H4xz_G3xy_a;
    abcd[iGrid*675+466] = 2.0E0*I_ESP_H3xyz_G3xy_a;
    abcd[iGrid*675+467] = 2.0E0*I_ESP_H3x2z_G3xy_a-1*I_ESP_F3x_G3xy;
    abcd[iGrid*675+468] = 2.0E0*I_ESP_H2x2yz_G3xy_a;
    abcd[iGrid*675+469] = 2.0E0*I_ESP_H2xy2z_G3xy_a-1*I_ESP_F2xy_G3xy;
    abcd[iGrid*675+470] = 2.0E0*I_ESP_H2x3z_G3xy_a-2*I_ESP_F2xz_G3xy;
    abcd[iGrid*675+471] = 2.0E0*I_ESP_Hx3yz_G3xy_a;
    abcd[iGrid*675+472] = 2.0E0*I_ESP_Hx2y2z_G3xy_a-1*I_ESP_Fx2y_G3xy;
    abcd[iGrid*675+473] = 2.0E0*I_ESP_Hxy3z_G3xy_a-2*I_ESP_Fxyz_G3xy;
    abcd[iGrid*675+474] = 2.0E0*I_ESP_Hx4z_G3xy_a-3*I_ESP_Fx2z_G3xy;
    abcd[iGrid*675+475] = 2.0E0*I_ESP_H4yz_G3xy_a;
    abcd[iGrid*675+476] = 2.0E0*I_ESP_H3y2z_G3xy_a-1*I_ESP_F3y_G3xy;
    abcd[iGrid*675+477] = 2.0E0*I_ESP_H2y3z_G3xy_a-2*I_ESP_F2yz_G3xy;
    abcd[iGrid*675+478] = 2.0E0*I_ESP_Hy4z_G3xy_a-3*I_ESP_Fy2z_G3xy;
    abcd[iGrid*675+479] = 2.0E0*I_ESP_H5z_G3xy_a-4*I_ESP_F3z_G3xy;
    abcd[iGrid*675+480] = 2.0E0*I_ESP_H4xz_G3xz_a;
    abcd[iGrid*675+481] = 2.0E0*I_ESP_H3xyz_G3xz_a;
    abcd[iGrid*675+482] = 2.0E0*I_ESP_H3x2z_G3xz_a-1*I_ESP_F3x_G3xz;
    abcd[iGrid*675+483] = 2.0E0*I_ESP_H2x2yz_G3xz_a;
    abcd[iGrid*675+484] = 2.0E0*I_ESP_H2xy2z_G3xz_a-1*I_ESP_F2xy_G3xz;
    abcd[iGrid*675+485] = 2.0E0*I_ESP_H2x3z_G3xz_a-2*I_ESP_F2xz_G3xz;
    abcd[iGrid*675+486] = 2.0E0*I_ESP_Hx3yz_G3xz_a;
    abcd[iGrid*675+487] = 2.0E0*I_ESP_Hx2y2z_G3xz_a-1*I_ESP_Fx2y_G3xz;
    abcd[iGrid*675+488] = 2.0E0*I_ESP_Hxy3z_G3xz_a-2*I_ESP_Fxyz_G3xz;
    abcd[iGrid*675+489] = 2.0E0*I_ESP_Hx4z_G3xz_a-3*I_ESP_Fx2z_G3xz;
    abcd[iGrid*675+490] = 2.0E0*I_ESP_H4yz_G3xz_a;
    abcd[iGrid*675+491] = 2.0E0*I_ESP_H3y2z_G3xz_a-1*I_ESP_F3y_G3xz;
    abcd[iGrid*675+492] = 2.0E0*I_ESP_H2y3z_G3xz_a-2*I_ESP_F2yz_G3xz;
    abcd[iGrid*675+493] = 2.0E0*I_ESP_Hy4z_G3xz_a-3*I_ESP_Fy2z_G3xz;
    abcd[iGrid*675+494] = 2.0E0*I_ESP_H5z_G3xz_a-4*I_ESP_F3z_G3xz;
    abcd[iGrid*675+495] = 2.0E0*I_ESP_H4xz_G2x2y_a;
    abcd[iGrid*675+496] = 2.0E0*I_ESP_H3xyz_G2x2y_a;
    abcd[iGrid*675+497] = 2.0E0*I_ESP_H3x2z_G2x2y_a-1*I_ESP_F3x_G2x2y;
    abcd[iGrid*675+498] = 2.0E0*I_ESP_H2x2yz_G2x2y_a;
    abcd[iGrid*675+499] = 2.0E0*I_ESP_H2xy2z_G2x2y_a-1*I_ESP_F2xy_G2x2y;
    abcd[iGrid*675+500] = 2.0E0*I_ESP_H2x3z_G2x2y_a-2*I_ESP_F2xz_G2x2y;
    abcd[iGrid*675+501] = 2.0E0*I_ESP_Hx3yz_G2x2y_a;
    abcd[iGrid*675+502] = 2.0E0*I_ESP_Hx2y2z_G2x2y_a-1*I_ESP_Fx2y_G2x2y;
    abcd[iGrid*675+503] = 2.0E0*I_ESP_Hxy3z_G2x2y_a-2*I_ESP_Fxyz_G2x2y;
    abcd[iGrid*675+504] = 2.0E0*I_ESP_Hx4z_G2x2y_a-3*I_ESP_Fx2z_G2x2y;
    abcd[iGrid*675+505] = 2.0E0*I_ESP_H4yz_G2x2y_a;
    abcd[iGrid*675+506] = 2.0E0*I_ESP_H3y2z_G2x2y_a-1*I_ESP_F3y_G2x2y;
    abcd[iGrid*675+507] = 2.0E0*I_ESP_H2y3z_G2x2y_a-2*I_ESP_F2yz_G2x2y;
    abcd[iGrid*675+508] = 2.0E0*I_ESP_Hy4z_G2x2y_a-3*I_ESP_Fy2z_G2x2y;
    abcd[iGrid*675+509] = 2.0E0*I_ESP_H5z_G2x2y_a-4*I_ESP_F3z_G2x2y;
    abcd[iGrid*675+510] = 2.0E0*I_ESP_H4xz_G2xyz_a;
    abcd[iGrid*675+511] = 2.0E0*I_ESP_H3xyz_G2xyz_a;
    abcd[iGrid*675+512] = 2.0E0*I_ESP_H3x2z_G2xyz_a-1*I_ESP_F3x_G2xyz;
    abcd[iGrid*675+513] = 2.0E0*I_ESP_H2x2yz_G2xyz_a;
    abcd[iGrid*675+514] = 2.0E0*I_ESP_H2xy2z_G2xyz_a-1*I_ESP_F2xy_G2xyz;
    abcd[iGrid*675+515] = 2.0E0*I_ESP_H2x3z_G2xyz_a-2*I_ESP_F2xz_G2xyz;
    abcd[iGrid*675+516] = 2.0E0*I_ESP_Hx3yz_G2xyz_a;
    abcd[iGrid*675+517] = 2.0E0*I_ESP_Hx2y2z_G2xyz_a-1*I_ESP_Fx2y_G2xyz;
    abcd[iGrid*675+518] = 2.0E0*I_ESP_Hxy3z_G2xyz_a-2*I_ESP_Fxyz_G2xyz;
    abcd[iGrid*675+519] = 2.0E0*I_ESP_Hx4z_G2xyz_a-3*I_ESP_Fx2z_G2xyz;
    abcd[iGrid*675+520] = 2.0E0*I_ESP_H4yz_G2xyz_a;
    abcd[iGrid*675+521] = 2.0E0*I_ESP_H3y2z_G2xyz_a-1*I_ESP_F3y_G2xyz;
    abcd[iGrid*675+522] = 2.0E0*I_ESP_H2y3z_G2xyz_a-2*I_ESP_F2yz_G2xyz;
    abcd[iGrid*675+523] = 2.0E0*I_ESP_Hy4z_G2xyz_a-3*I_ESP_Fy2z_G2xyz;
    abcd[iGrid*675+524] = 2.0E0*I_ESP_H5z_G2xyz_a-4*I_ESP_F3z_G2xyz;
    abcd[iGrid*675+525] = 2.0E0*I_ESP_H4xz_G2x2z_a;
    abcd[iGrid*675+526] = 2.0E0*I_ESP_H3xyz_G2x2z_a;
    abcd[iGrid*675+527] = 2.0E0*I_ESP_H3x2z_G2x2z_a-1*I_ESP_F3x_G2x2z;
    abcd[iGrid*675+528] = 2.0E0*I_ESP_H2x2yz_G2x2z_a;
    abcd[iGrid*675+529] = 2.0E0*I_ESP_H2xy2z_G2x2z_a-1*I_ESP_F2xy_G2x2z;
    abcd[iGrid*675+530] = 2.0E0*I_ESP_H2x3z_G2x2z_a-2*I_ESP_F2xz_G2x2z;
    abcd[iGrid*675+531] = 2.0E0*I_ESP_Hx3yz_G2x2z_a;
    abcd[iGrid*675+532] = 2.0E0*I_ESP_Hx2y2z_G2x2z_a-1*I_ESP_Fx2y_G2x2z;
    abcd[iGrid*675+533] = 2.0E0*I_ESP_Hxy3z_G2x2z_a-2*I_ESP_Fxyz_G2x2z;
    abcd[iGrid*675+534] = 2.0E0*I_ESP_Hx4z_G2x2z_a-3*I_ESP_Fx2z_G2x2z;
    abcd[iGrid*675+535] = 2.0E0*I_ESP_H4yz_G2x2z_a;
    abcd[iGrid*675+536] = 2.0E0*I_ESP_H3y2z_G2x2z_a-1*I_ESP_F3y_G2x2z;
    abcd[iGrid*675+537] = 2.0E0*I_ESP_H2y3z_G2x2z_a-2*I_ESP_F2yz_G2x2z;
    abcd[iGrid*675+538] = 2.0E0*I_ESP_Hy4z_G2x2z_a-3*I_ESP_Fy2z_G2x2z;
    abcd[iGrid*675+539] = 2.0E0*I_ESP_H5z_G2x2z_a-4*I_ESP_F3z_G2x2z;
    abcd[iGrid*675+540] = 2.0E0*I_ESP_H4xz_Gx3y_a;
    abcd[iGrid*675+541] = 2.0E0*I_ESP_H3xyz_Gx3y_a;
    abcd[iGrid*675+542] = 2.0E0*I_ESP_H3x2z_Gx3y_a-1*I_ESP_F3x_Gx3y;
    abcd[iGrid*675+543] = 2.0E0*I_ESP_H2x2yz_Gx3y_a;
    abcd[iGrid*675+544] = 2.0E0*I_ESP_H2xy2z_Gx3y_a-1*I_ESP_F2xy_Gx3y;
    abcd[iGrid*675+545] = 2.0E0*I_ESP_H2x3z_Gx3y_a-2*I_ESP_F2xz_Gx3y;
    abcd[iGrid*675+546] = 2.0E0*I_ESP_Hx3yz_Gx3y_a;
    abcd[iGrid*675+547] = 2.0E0*I_ESP_Hx2y2z_Gx3y_a-1*I_ESP_Fx2y_Gx3y;
    abcd[iGrid*675+548] = 2.0E0*I_ESP_Hxy3z_Gx3y_a-2*I_ESP_Fxyz_Gx3y;
    abcd[iGrid*675+549] = 2.0E0*I_ESP_Hx4z_Gx3y_a-3*I_ESP_Fx2z_Gx3y;
    abcd[iGrid*675+550] = 2.0E0*I_ESP_H4yz_Gx3y_a;
    abcd[iGrid*675+551] = 2.0E0*I_ESP_H3y2z_Gx3y_a-1*I_ESP_F3y_Gx3y;
    abcd[iGrid*675+552] = 2.0E0*I_ESP_H2y3z_Gx3y_a-2*I_ESP_F2yz_Gx3y;
    abcd[iGrid*675+553] = 2.0E0*I_ESP_Hy4z_Gx3y_a-3*I_ESP_Fy2z_Gx3y;
    abcd[iGrid*675+554] = 2.0E0*I_ESP_H5z_Gx3y_a-4*I_ESP_F3z_Gx3y;
    abcd[iGrid*675+555] = 2.0E0*I_ESP_H4xz_Gx2yz_a;
    abcd[iGrid*675+556] = 2.0E0*I_ESP_H3xyz_Gx2yz_a;
    abcd[iGrid*675+557] = 2.0E0*I_ESP_H3x2z_Gx2yz_a-1*I_ESP_F3x_Gx2yz;
    abcd[iGrid*675+558] = 2.0E0*I_ESP_H2x2yz_Gx2yz_a;
    abcd[iGrid*675+559] = 2.0E0*I_ESP_H2xy2z_Gx2yz_a-1*I_ESP_F2xy_Gx2yz;
    abcd[iGrid*675+560] = 2.0E0*I_ESP_H2x3z_Gx2yz_a-2*I_ESP_F2xz_Gx2yz;
    abcd[iGrid*675+561] = 2.0E0*I_ESP_Hx3yz_Gx2yz_a;
    abcd[iGrid*675+562] = 2.0E0*I_ESP_Hx2y2z_Gx2yz_a-1*I_ESP_Fx2y_Gx2yz;
    abcd[iGrid*675+563] = 2.0E0*I_ESP_Hxy3z_Gx2yz_a-2*I_ESP_Fxyz_Gx2yz;
    abcd[iGrid*675+564] = 2.0E0*I_ESP_Hx4z_Gx2yz_a-3*I_ESP_Fx2z_Gx2yz;
    abcd[iGrid*675+565] = 2.0E0*I_ESP_H4yz_Gx2yz_a;
    abcd[iGrid*675+566] = 2.0E0*I_ESP_H3y2z_Gx2yz_a-1*I_ESP_F3y_Gx2yz;
    abcd[iGrid*675+567] = 2.0E0*I_ESP_H2y3z_Gx2yz_a-2*I_ESP_F2yz_Gx2yz;
    abcd[iGrid*675+568] = 2.0E0*I_ESP_Hy4z_Gx2yz_a-3*I_ESP_Fy2z_Gx2yz;
    abcd[iGrid*675+569] = 2.0E0*I_ESP_H5z_Gx2yz_a-4*I_ESP_F3z_Gx2yz;
    abcd[iGrid*675+570] = 2.0E0*I_ESP_H4xz_Gxy2z_a;
    abcd[iGrid*675+571] = 2.0E0*I_ESP_H3xyz_Gxy2z_a;
    abcd[iGrid*675+572] = 2.0E0*I_ESP_H3x2z_Gxy2z_a-1*I_ESP_F3x_Gxy2z;
    abcd[iGrid*675+573] = 2.0E0*I_ESP_H2x2yz_Gxy2z_a;
    abcd[iGrid*675+574] = 2.0E0*I_ESP_H2xy2z_Gxy2z_a-1*I_ESP_F2xy_Gxy2z;
    abcd[iGrid*675+575] = 2.0E0*I_ESP_H2x3z_Gxy2z_a-2*I_ESP_F2xz_Gxy2z;
    abcd[iGrid*675+576] = 2.0E0*I_ESP_Hx3yz_Gxy2z_a;
    abcd[iGrid*675+577] = 2.0E0*I_ESP_Hx2y2z_Gxy2z_a-1*I_ESP_Fx2y_Gxy2z;
    abcd[iGrid*675+578] = 2.0E0*I_ESP_Hxy3z_Gxy2z_a-2*I_ESP_Fxyz_Gxy2z;
    abcd[iGrid*675+579] = 2.0E0*I_ESP_Hx4z_Gxy2z_a-3*I_ESP_Fx2z_Gxy2z;
    abcd[iGrid*675+580] = 2.0E0*I_ESP_H4yz_Gxy2z_a;
    abcd[iGrid*675+581] = 2.0E0*I_ESP_H3y2z_Gxy2z_a-1*I_ESP_F3y_Gxy2z;
    abcd[iGrid*675+582] = 2.0E0*I_ESP_H2y3z_Gxy2z_a-2*I_ESP_F2yz_Gxy2z;
    abcd[iGrid*675+583] = 2.0E0*I_ESP_Hy4z_Gxy2z_a-3*I_ESP_Fy2z_Gxy2z;
    abcd[iGrid*675+584] = 2.0E0*I_ESP_H5z_Gxy2z_a-4*I_ESP_F3z_Gxy2z;
    abcd[iGrid*675+585] = 2.0E0*I_ESP_H4xz_Gx3z_a;
    abcd[iGrid*675+586] = 2.0E0*I_ESP_H3xyz_Gx3z_a;
    abcd[iGrid*675+587] = 2.0E0*I_ESP_H3x2z_Gx3z_a-1*I_ESP_F3x_Gx3z;
    abcd[iGrid*675+588] = 2.0E0*I_ESP_H2x2yz_Gx3z_a;
    abcd[iGrid*675+589] = 2.0E0*I_ESP_H2xy2z_Gx3z_a-1*I_ESP_F2xy_Gx3z;
    abcd[iGrid*675+590] = 2.0E0*I_ESP_H2x3z_Gx3z_a-2*I_ESP_F2xz_Gx3z;
    abcd[iGrid*675+591] = 2.0E0*I_ESP_Hx3yz_Gx3z_a;
    abcd[iGrid*675+592] = 2.0E0*I_ESP_Hx2y2z_Gx3z_a-1*I_ESP_Fx2y_Gx3z;
    abcd[iGrid*675+593] = 2.0E0*I_ESP_Hxy3z_Gx3z_a-2*I_ESP_Fxyz_Gx3z;
    abcd[iGrid*675+594] = 2.0E0*I_ESP_Hx4z_Gx3z_a-3*I_ESP_Fx2z_Gx3z;
    abcd[iGrid*675+595] = 2.0E0*I_ESP_H4yz_Gx3z_a;
    abcd[iGrid*675+596] = 2.0E0*I_ESP_H3y2z_Gx3z_a-1*I_ESP_F3y_Gx3z;
    abcd[iGrid*675+597] = 2.0E0*I_ESP_H2y3z_Gx3z_a-2*I_ESP_F2yz_Gx3z;
    abcd[iGrid*675+598] = 2.0E0*I_ESP_Hy4z_Gx3z_a-3*I_ESP_Fy2z_Gx3z;
    abcd[iGrid*675+599] = 2.0E0*I_ESP_H5z_Gx3z_a-4*I_ESP_F3z_Gx3z;
    abcd[iGrid*675+600] = 2.0E0*I_ESP_H4xz_G4y_a;
    abcd[iGrid*675+601] = 2.0E0*I_ESP_H3xyz_G4y_a;
    abcd[iGrid*675+602] = 2.0E0*I_ESP_H3x2z_G4y_a-1*I_ESP_F3x_G4y;
    abcd[iGrid*675+603] = 2.0E0*I_ESP_H2x2yz_G4y_a;
    abcd[iGrid*675+604] = 2.0E0*I_ESP_H2xy2z_G4y_a-1*I_ESP_F2xy_G4y;
    abcd[iGrid*675+605] = 2.0E0*I_ESP_H2x3z_G4y_a-2*I_ESP_F2xz_G4y;
    abcd[iGrid*675+606] = 2.0E0*I_ESP_Hx3yz_G4y_a;
    abcd[iGrid*675+607] = 2.0E0*I_ESP_Hx2y2z_G4y_a-1*I_ESP_Fx2y_G4y;
    abcd[iGrid*675+608] = 2.0E0*I_ESP_Hxy3z_G4y_a-2*I_ESP_Fxyz_G4y;
    abcd[iGrid*675+609] = 2.0E0*I_ESP_Hx4z_G4y_a-3*I_ESP_Fx2z_G4y;
    abcd[iGrid*675+610] = 2.0E0*I_ESP_H4yz_G4y_a;
    abcd[iGrid*675+611] = 2.0E0*I_ESP_H3y2z_G4y_a-1*I_ESP_F3y_G4y;
    abcd[iGrid*675+612] = 2.0E0*I_ESP_H2y3z_G4y_a-2*I_ESP_F2yz_G4y;
    abcd[iGrid*675+613] = 2.0E0*I_ESP_Hy4z_G4y_a-3*I_ESP_Fy2z_G4y;
    abcd[iGrid*675+614] = 2.0E0*I_ESP_H5z_G4y_a-4*I_ESP_F3z_G4y;
    abcd[iGrid*675+615] = 2.0E0*I_ESP_H4xz_G3yz_a;
    abcd[iGrid*675+616] = 2.0E0*I_ESP_H3xyz_G3yz_a;
    abcd[iGrid*675+617] = 2.0E0*I_ESP_H3x2z_G3yz_a-1*I_ESP_F3x_G3yz;
    abcd[iGrid*675+618] = 2.0E0*I_ESP_H2x2yz_G3yz_a;
    abcd[iGrid*675+619] = 2.0E0*I_ESP_H2xy2z_G3yz_a-1*I_ESP_F2xy_G3yz;
    abcd[iGrid*675+620] = 2.0E0*I_ESP_H2x3z_G3yz_a-2*I_ESP_F2xz_G3yz;
    abcd[iGrid*675+621] = 2.0E0*I_ESP_Hx3yz_G3yz_a;
    abcd[iGrid*675+622] = 2.0E0*I_ESP_Hx2y2z_G3yz_a-1*I_ESP_Fx2y_G3yz;
    abcd[iGrid*675+623] = 2.0E0*I_ESP_Hxy3z_G3yz_a-2*I_ESP_Fxyz_G3yz;
    abcd[iGrid*675+624] = 2.0E0*I_ESP_Hx4z_G3yz_a-3*I_ESP_Fx2z_G3yz;
    abcd[iGrid*675+625] = 2.0E0*I_ESP_H4yz_G3yz_a;
    abcd[iGrid*675+626] = 2.0E0*I_ESP_H3y2z_G3yz_a-1*I_ESP_F3y_G3yz;
    abcd[iGrid*675+627] = 2.0E0*I_ESP_H2y3z_G3yz_a-2*I_ESP_F2yz_G3yz;
    abcd[iGrid*675+628] = 2.0E0*I_ESP_Hy4z_G3yz_a-3*I_ESP_Fy2z_G3yz;
    abcd[iGrid*675+629] = 2.0E0*I_ESP_H5z_G3yz_a-4*I_ESP_F3z_G3yz;
    abcd[iGrid*675+630] = 2.0E0*I_ESP_H4xz_G2y2z_a;
    abcd[iGrid*675+631] = 2.0E0*I_ESP_H3xyz_G2y2z_a;
    abcd[iGrid*675+632] = 2.0E0*I_ESP_H3x2z_G2y2z_a-1*I_ESP_F3x_G2y2z;
    abcd[iGrid*675+633] = 2.0E0*I_ESP_H2x2yz_G2y2z_a;
    abcd[iGrid*675+634] = 2.0E0*I_ESP_H2xy2z_G2y2z_a-1*I_ESP_F2xy_G2y2z;
    abcd[iGrid*675+635] = 2.0E0*I_ESP_H2x3z_G2y2z_a-2*I_ESP_F2xz_G2y2z;
    abcd[iGrid*675+636] = 2.0E0*I_ESP_Hx3yz_G2y2z_a;
    abcd[iGrid*675+637] = 2.0E0*I_ESP_Hx2y2z_G2y2z_a-1*I_ESP_Fx2y_G2y2z;
    abcd[iGrid*675+638] = 2.0E0*I_ESP_Hxy3z_G2y2z_a-2*I_ESP_Fxyz_G2y2z;
    abcd[iGrid*675+639] = 2.0E0*I_ESP_Hx4z_G2y2z_a-3*I_ESP_Fx2z_G2y2z;
    abcd[iGrid*675+640] = 2.0E0*I_ESP_H4yz_G2y2z_a;
    abcd[iGrid*675+641] = 2.0E0*I_ESP_H3y2z_G2y2z_a-1*I_ESP_F3y_G2y2z;
    abcd[iGrid*675+642] = 2.0E0*I_ESP_H2y3z_G2y2z_a-2*I_ESP_F2yz_G2y2z;
    abcd[iGrid*675+643] = 2.0E0*I_ESP_Hy4z_G2y2z_a-3*I_ESP_Fy2z_G2y2z;
    abcd[iGrid*675+644] = 2.0E0*I_ESP_H5z_G2y2z_a-4*I_ESP_F3z_G2y2z;
    abcd[iGrid*675+645] = 2.0E0*I_ESP_H4xz_Gy3z_a;
    abcd[iGrid*675+646] = 2.0E0*I_ESP_H3xyz_Gy3z_a;
    abcd[iGrid*675+647] = 2.0E0*I_ESP_H3x2z_Gy3z_a-1*I_ESP_F3x_Gy3z;
    abcd[iGrid*675+648] = 2.0E0*I_ESP_H2x2yz_Gy3z_a;
    abcd[iGrid*675+649] = 2.0E0*I_ESP_H2xy2z_Gy3z_a-1*I_ESP_F2xy_Gy3z;
    abcd[iGrid*675+650] = 2.0E0*I_ESP_H2x3z_Gy3z_a-2*I_ESP_F2xz_Gy3z;
    abcd[iGrid*675+651] = 2.0E0*I_ESP_Hx3yz_Gy3z_a;
    abcd[iGrid*675+652] = 2.0E0*I_ESP_Hx2y2z_Gy3z_a-1*I_ESP_Fx2y_Gy3z;
    abcd[iGrid*675+653] = 2.0E0*I_ESP_Hxy3z_Gy3z_a-2*I_ESP_Fxyz_Gy3z;
    abcd[iGrid*675+654] = 2.0E0*I_ESP_Hx4z_Gy3z_a-3*I_ESP_Fx2z_Gy3z;
    abcd[iGrid*675+655] = 2.0E0*I_ESP_H4yz_Gy3z_a;
    abcd[iGrid*675+656] = 2.0E0*I_ESP_H3y2z_Gy3z_a-1*I_ESP_F3y_Gy3z;
    abcd[iGrid*675+657] = 2.0E0*I_ESP_H2y3z_Gy3z_a-2*I_ESP_F2yz_Gy3z;
    abcd[iGrid*675+658] = 2.0E0*I_ESP_Hy4z_Gy3z_a-3*I_ESP_Fy2z_Gy3z;
    abcd[iGrid*675+659] = 2.0E0*I_ESP_H5z_Gy3z_a-4*I_ESP_F3z_Gy3z;
    abcd[iGrid*675+660] = 2.0E0*I_ESP_H4xz_G4z_a;
    abcd[iGrid*675+661] = 2.0E0*I_ESP_H3xyz_G4z_a;
    abcd[iGrid*675+662] = 2.0E0*I_ESP_H3x2z_G4z_a-1*I_ESP_F3x_G4z;
    abcd[iGrid*675+663] = 2.0E0*I_ESP_H2x2yz_G4z_a;
    abcd[iGrid*675+664] = 2.0E0*I_ESP_H2xy2z_G4z_a-1*I_ESP_F2xy_G4z;
    abcd[iGrid*675+665] = 2.0E0*I_ESP_H2x3z_G4z_a-2*I_ESP_F2xz_G4z;
    abcd[iGrid*675+666] = 2.0E0*I_ESP_Hx3yz_G4z_a;
    abcd[iGrid*675+667] = 2.0E0*I_ESP_Hx2y2z_G4z_a-1*I_ESP_Fx2y_G4z;
    abcd[iGrid*675+668] = 2.0E0*I_ESP_Hxy3z_G4z_a-2*I_ESP_Fxyz_G4z;
    abcd[iGrid*675+669] = 2.0E0*I_ESP_Hx4z_G4z_a-3*I_ESP_Fx2z_G4z;
    abcd[iGrid*675+670] = 2.0E0*I_ESP_H4yz_G4z_a;
    abcd[iGrid*675+671] = 2.0E0*I_ESP_H3y2z_G4z_a-1*I_ESP_F3y_G4z;
    abcd[iGrid*675+672] = 2.0E0*I_ESP_H2y3z_G4z_a-2*I_ESP_F2yz_G4z;
    abcd[iGrid*675+673] = 2.0E0*I_ESP_Hy4z_G4z_a-3*I_ESP_Fy2z_G4z;
    abcd[iGrid*675+674] = 2.0E0*I_ESP_H5z_G4z_a-4*I_ESP_F3z_G4z;
  }
}
