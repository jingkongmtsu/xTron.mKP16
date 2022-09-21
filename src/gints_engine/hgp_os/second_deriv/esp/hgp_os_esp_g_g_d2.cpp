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
// BRA1  BRA1
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_esp_g_g_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_N10x_S_aa = 0.0E0;
    Double I_ESP_N9xy_S_aa = 0.0E0;
    Double I_ESP_N9xz_S_aa = 0.0E0;
    Double I_ESP_N8x2y_S_aa = 0.0E0;
    Double I_ESP_N8xyz_S_aa = 0.0E0;
    Double I_ESP_N8x2z_S_aa = 0.0E0;
    Double I_ESP_N7x3y_S_aa = 0.0E0;
    Double I_ESP_N7x2yz_S_aa = 0.0E0;
    Double I_ESP_N7xy2z_S_aa = 0.0E0;
    Double I_ESP_N7x3z_S_aa = 0.0E0;
    Double I_ESP_N6x4y_S_aa = 0.0E0;
    Double I_ESP_N6x3yz_S_aa = 0.0E0;
    Double I_ESP_N6x2y2z_S_aa = 0.0E0;
    Double I_ESP_N6xy3z_S_aa = 0.0E0;
    Double I_ESP_N6x4z_S_aa = 0.0E0;
    Double I_ESP_N5x5y_S_aa = 0.0E0;
    Double I_ESP_N5x4yz_S_aa = 0.0E0;
    Double I_ESP_N5x3y2z_S_aa = 0.0E0;
    Double I_ESP_N5x2y3z_S_aa = 0.0E0;
    Double I_ESP_N5xy4z_S_aa = 0.0E0;
    Double I_ESP_N5x5z_S_aa = 0.0E0;
    Double I_ESP_N4x6y_S_aa = 0.0E0;
    Double I_ESP_N4x5yz_S_aa = 0.0E0;
    Double I_ESP_N4x4y2z_S_aa = 0.0E0;
    Double I_ESP_N4x3y3z_S_aa = 0.0E0;
    Double I_ESP_N4x2y4z_S_aa = 0.0E0;
    Double I_ESP_N4xy5z_S_aa = 0.0E0;
    Double I_ESP_N4x6z_S_aa = 0.0E0;
    Double I_ESP_N3x7y_S_aa = 0.0E0;
    Double I_ESP_N3x6yz_S_aa = 0.0E0;
    Double I_ESP_N3x5y2z_S_aa = 0.0E0;
    Double I_ESP_N3x4y3z_S_aa = 0.0E0;
    Double I_ESP_N3x3y4z_S_aa = 0.0E0;
    Double I_ESP_N3x2y5z_S_aa = 0.0E0;
    Double I_ESP_N3xy6z_S_aa = 0.0E0;
    Double I_ESP_N3x7z_S_aa = 0.0E0;
    Double I_ESP_N2x8y_S_aa = 0.0E0;
    Double I_ESP_N2x7yz_S_aa = 0.0E0;
    Double I_ESP_N2x6y2z_S_aa = 0.0E0;
    Double I_ESP_N2x5y3z_S_aa = 0.0E0;
    Double I_ESP_N2x4y4z_S_aa = 0.0E0;
    Double I_ESP_N2x3y5z_S_aa = 0.0E0;
    Double I_ESP_N2x2y6z_S_aa = 0.0E0;
    Double I_ESP_N2xy7z_S_aa = 0.0E0;
    Double I_ESP_N2x8z_S_aa = 0.0E0;
    Double I_ESP_Nx9y_S_aa = 0.0E0;
    Double I_ESP_Nx8yz_S_aa = 0.0E0;
    Double I_ESP_Nx7y2z_S_aa = 0.0E0;
    Double I_ESP_Nx6y3z_S_aa = 0.0E0;
    Double I_ESP_Nx5y4z_S_aa = 0.0E0;
    Double I_ESP_Nx4y5z_S_aa = 0.0E0;
    Double I_ESP_Nx3y6z_S_aa = 0.0E0;
    Double I_ESP_Nx2y7z_S_aa = 0.0E0;
    Double I_ESP_Nxy8z_S_aa = 0.0E0;
    Double I_ESP_Nx9z_S_aa = 0.0E0;
    Double I_ESP_N10y_S_aa = 0.0E0;
    Double I_ESP_N9yz_S_aa = 0.0E0;
    Double I_ESP_N8y2z_S_aa = 0.0E0;
    Double I_ESP_N7y3z_S_aa = 0.0E0;
    Double I_ESP_N6y4z_S_aa = 0.0E0;
    Double I_ESP_N5y5z_S_aa = 0.0E0;
    Double I_ESP_N4y6z_S_aa = 0.0E0;
    Double I_ESP_N3y7z_S_aa = 0.0E0;
    Double I_ESP_N2y8z_S_aa = 0.0E0;
    Double I_ESP_Ny9z_S_aa = 0.0E0;
    Double I_ESP_N10z_S_aa = 0.0E0;
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
      Double I_ESP_S_S_M10_vrr  = 0.0E0;

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
      double I_ESP_S_S_M10_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER55;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER53*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER51*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER49*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER47*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER45*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER43*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER41*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER39*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER37*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER35*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER33*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER31*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER29*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER27*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER25*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = 1.0E0+u2*ONEOVER23*I_ESP_S_S_M10_vrr;
        I_ESP_S_S_M10_vrr = ONEOVER21*I_ESP_S_S_M10_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M10_vrr  = f*I_ESP_S_S_M10_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ESP_S_S_M9_vrr  = ONEOVER19*(u2*I_ESP_S_S_M10_vrr+f);
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
        I_ESP_S_S_M10_vrr_d = oneO2u*(19.0E0*I_ESP_S_S_M9_vrr_d-f);

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
        I_ESP_S_S_M10_vrr = static_cast<Double>(I_ESP_S_S_M10_vrr_d);

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
        I_ESP_S_S_M10_vrr = oneO2u*(19.0E0*I_ESP_S_S_M9_vrr-f);

#endif

      }


        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M9
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M9
         * RHS shell quartet name: SQ_ESP_S_S_M10
         ************************************************************/
        Double I_ESP_Px_S_M9_vrr = PAX*I_ESP_S_S_M9_vrr-PRX*I_ESP_S_S_M10_vrr;
        Double I_ESP_Py_S_M9_vrr = PAY*I_ESP_S_S_M9_vrr-PRY*I_ESP_S_S_M10_vrr;
        Double I_ESP_Pz_S_M9_vrr = PAZ*I_ESP_S_S_M9_vrr-PRZ*I_ESP_S_S_M10_vrr;

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
         * shell quartet name: SQ_ESP_D_S_M8
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M8
         * RHS shell quartet name: SQ_ESP_P_S_M9
         * RHS shell quartet name: SQ_ESP_S_S_M8
         * RHS shell quartet name: SQ_ESP_S_S_M9
         ************************************************************/
        Double I_ESP_D2x_S_M8_vrr = PAX*I_ESP_Px_S_M8_vrr-PRX*I_ESP_Px_S_M9_vrr+oned2z*I_ESP_S_S_M8_vrr-oned2z*I_ESP_S_S_M9_vrr;
        Double I_ESP_D2y_S_M8_vrr = PAY*I_ESP_Py_S_M8_vrr-PRY*I_ESP_Py_S_M9_vrr+oned2z*I_ESP_S_S_M8_vrr-oned2z*I_ESP_S_S_M9_vrr;
        Double I_ESP_D2z_S_M8_vrr = PAZ*I_ESP_Pz_S_M8_vrr-PRZ*I_ESP_Pz_S_M9_vrr+oned2z*I_ESP_S_S_M8_vrr-oned2z*I_ESP_S_S_M9_vrr;

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
         * shell quartet name: SQ_ESP_F_S_M7
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M7
         * RHS shell quartet name: SQ_ESP_D_S_M8
         * RHS shell quartet name: SQ_ESP_P_S_M7
         * RHS shell quartet name: SQ_ESP_P_S_M8
         ************************************************************/
        Double I_ESP_F3x_S_M7_vrr = PAX*I_ESP_D2x_S_M7_vrr-PRX*I_ESP_D2x_S_M8_vrr+2*oned2z*I_ESP_Px_S_M7_vrr-2*oned2z*I_ESP_Px_S_M8_vrr;
        Double I_ESP_F3y_S_M7_vrr = PAY*I_ESP_D2y_S_M7_vrr-PRY*I_ESP_D2y_S_M8_vrr+2*oned2z*I_ESP_Py_S_M7_vrr-2*oned2z*I_ESP_Py_S_M8_vrr;
        Double I_ESP_F3z_S_M7_vrr = PAZ*I_ESP_D2z_S_M7_vrr-PRZ*I_ESP_D2z_S_M8_vrr+2*oned2z*I_ESP_Pz_S_M7_vrr-2*oned2z*I_ESP_Pz_S_M8_vrr;

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
         * shell quartet name: SQ_ESP_G_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 12 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M6
         * RHS shell quartet name: SQ_ESP_F_S_M7
         * RHS shell quartet name: SQ_ESP_D_S_M6
         * RHS shell quartet name: SQ_ESP_D_S_M7
         ************************************************************/
        Double I_ESP_G4x_S_M6_vrr = PAX*I_ESP_F3x_S_M6_vrr-PRX*I_ESP_F3x_S_M7_vrr+3*oned2z*I_ESP_D2x_S_M6_vrr-3*oned2z*I_ESP_D2x_S_M7_vrr;
        Double I_ESP_G4y_S_M6_vrr = PAY*I_ESP_F3y_S_M6_vrr-PRY*I_ESP_F3y_S_M7_vrr+3*oned2z*I_ESP_D2y_S_M6_vrr-3*oned2z*I_ESP_D2y_S_M7_vrr;
        Double I_ESP_G4z_S_M6_vrr = PAZ*I_ESP_F3z_S_M6_vrr-PRZ*I_ESP_F3z_S_M7_vrr+3*oned2z*I_ESP_D2z_S_M6_vrr-3*oned2z*I_ESP_D2z_S_M7_vrr;

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
         * shell quartet name: SQ_ESP_H_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 13 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M5
         * RHS shell quartet name: SQ_ESP_G_S_M6
         * RHS shell quartet name: SQ_ESP_F_S_M5
         * RHS shell quartet name: SQ_ESP_F_S_M6
         ************************************************************/
        Double I_ESP_H5x_S_M5_vrr = PAX*I_ESP_G4x_S_M5_vrr-PRX*I_ESP_G4x_S_M6_vrr+4*oned2z*I_ESP_F3x_S_M5_vrr-4*oned2z*I_ESP_F3x_S_M6_vrr;
        Double I_ESP_H4xy_S_M5_vrr = PAY*I_ESP_G4x_S_M5_vrr-PRY*I_ESP_G4x_S_M6_vrr;
        Double I_ESP_H4xz_S_M5_vrr = PAZ*I_ESP_G4x_S_M5_vrr-PRZ*I_ESP_G4x_S_M6_vrr;
        Double I_ESP_Hx4y_S_M5_vrr = PAX*I_ESP_G4y_S_M5_vrr-PRX*I_ESP_G4y_S_M6_vrr;
        Double I_ESP_Hx4z_S_M5_vrr = PAX*I_ESP_G4z_S_M5_vrr-PRX*I_ESP_G4z_S_M6_vrr;
        Double I_ESP_H5y_S_M5_vrr = PAY*I_ESP_G4y_S_M5_vrr-PRY*I_ESP_G4y_S_M6_vrr+4*oned2z*I_ESP_F3y_S_M5_vrr-4*oned2z*I_ESP_F3y_S_M6_vrr;
        Double I_ESP_H4yz_S_M5_vrr = PAZ*I_ESP_G4y_S_M5_vrr-PRZ*I_ESP_G4y_S_M6_vrr;
        Double I_ESP_H5z_S_M5_vrr = PAZ*I_ESP_G4z_S_M5_vrr-PRZ*I_ESP_G4z_S_M6_vrr+4*oned2z*I_ESP_F3z_S_M5_vrr-4*oned2z*I_ESP_F3z_S_M6_vrr;

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
         * shell quartet name: SQ_ESP_I_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 14 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M4
         * RHS shell quartet name: SQ_ESP_H_S_M5
         * RHS shell quartet name: SQ_ESP_G_S_M4
         * RHS shell quartet name: SQ_ESP_G_S_M5
         ************************************************************/
        Double I_ESP_I6x_S_M4_vrr = PAX*I_ESP_H5x_S_M4_vrr-PRX*I_ESP_H5x_S_M5_vrr+5*oned2z*I_ESP_G4x_S_M4_vrr-5*oned2z*I_ESP_G4x_S_M5_vrr;
        Double I_ESP_I5xy_S_M4_vrr = PAY*I_ESP_H5x_S_M4_vrr-PRY*I_ESP_H5x_S_M5_vrr;
        Double I_ESP_I5xz_S_M4_vrr = PAZ*I_ESP_H5x_S_M4_vrr-PRZ*I_ESP_H5x_S_M5_vrr;
        Double I_ESP_I4x2y_S_M4_vrr = PAY*I_ESP_H4xy_S_M4_vrr-PRY*I_ESP_H4xy_S_M5_vrr+oned2z*I_ESP_G4x_S_M4_vrr-oned2z*I_ESP_G4x_S_M5_vrr;
        Double I_ESP_I4x2z_S_M4_vrr = PAZ*I_ESP_H4xz_S_M4_vrr-PRZ*I_ESP_H4xz_S_M5_vrr+oned2z*I_ESP_G4x_S_M4_vrr-oned2z*I_ESP_G4x_S_M5_vrr;
        Double I_ESP_I2x4y_S_M4_vrr = PAX*I_ESP_Hx4y_S_M4_vrr-PRX*I_ESP_Hx4y_S_M5_vrr+oned2z*I_ESP_G4y_S_M4_vrr-oned2z*I_ESP_G4y_S_M5_vrr;
        Double I_ESP_I2x4z_S_M4_vrr = PAX*I_ESP_Hx4z_S_M4_vrr-PRX*I_ESP_Hx4z_S_M5_vrr+oned2z*I_ESP_G4z_S_M4_vrr-oned2z*I_ESP_G4z_S_M5_vrr;
        Double I_ESP_Ix5y_S_M4_vrr = PAX*I_ESP_H5y_S_M4_vrr-PRX*I_ESP_H5y_S_M5_vrr;
        Double I_ESP_Ix5z_S_M4_vrr = PAX*I_ESP_H5z_S_M4_vrr-PRX*I_ESP_H5z_S_M5_vrr;
        Double I_ESP_I6y_S_M4_vrr = PAY*I_ESP_H5y_S_M4_vrr-PRY*I_ESP_H5y_S_M5_vrr+5*oned2z*I_ESP_G4y_S_M4_vrr-5*oned2z*I_ESP_G4y_S_M5_vrr;
        Double I_ESP_I5yz_S_M4_vrr = PAZ*I_ESP_H5y_S_M4_vrr-PRZ*I_ESP_H5y_S_M5_vrr;
        Double I_ESP_I4y2z_S_M4_vrr = PAZ*I_ESP_H4yz_S_M4_vrr-PRZ*I_ESP_H4yz_S_M5_vrr+oned2z*I_ESP_G4y_S_M4_vrr-oned2z*I_ESP_G4y_S_M5_vrr;
        Double I_ESP_Iy5z_S_M4_vrr = PAY*I_ESP_H5z_S_M4_vrr-PRY*I_ESP_H5z_S_M5_vrr;
        Double I_ESP_I6z_S_M4_vrr = PAZ*I_ESP_H5z_S_M4_vrr-PRZ*I_ESP_H5z_S_M5_vrr+5*oned2z*I_ESP_G4z_S_M4_vrr-5*oned2z*I_ESP_G4z_S_M5_vrr;

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
         * shell quartet name: SQ_ESP_K_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 16 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S_M3
         * RHS shell quartet name: SQ_ESP_I_S_M4
         * RHS shell quartet name: SQ_ESP_H_S_M3
         * RHS shell quartet name: SQ_ESP_H_S_M4
         ************************************************************/
        Double I_ESP_K7x_S_M3_vrr = PAX*I_ESP_I6x_S_M3_vrr-PRX*I_ESP_I6x_S_M4_vrr+6*oned2z*I_ESP_H5x_S_M3_vrr-6*oned2z*I_ESP_H5x_S_M4_vrr;
        Double I_ESP_K6xy_S_M3_vrr = PAY*I_ESP_I6x_S_M3_vrr-PRY*I_ESP_I6x_S_M4_vrr;
        Double I_ESP_K6xz_S_M3_vrr = PAZ*I_ESP_I6x_S_M3_vrr-PRZ*I_ESP_I6x_S_M4_vrr;
        Double I_ESP_K5x2y_S_M3_vrr = PAY*I_ESP_I5xy_S_M3_vrr-PRY*I_ESP_I5xy_S_M4_vrr+oned2z*I_ESP_H5x_S_M3_vrr-oned2z*I_ESP_H5x_S_M4_vrr;
        Double I_ESP_K5x2z_S_M3_vrr = PAZ*I_ESP_I5xz_S_M3_vrr-PRZ*I_ESP_I5xz_S_M4_vrr+oned2z*I_ESP_H5x_S_M3_vrr-oned2z*I_ESP_H5x_S_M4_vrr;
        Double I_ESP_K4x3y_S_M3_vrr = PAY*I_ESP_I4x2y_S_M3_vrr-PRY*I_ESP_I4x2y_S_M4_vrr+2*oned2z*I_ESP_H4xy_S_M3_vrr-2*oned2z*I_ESP_H4xy_S_M4_vrr;
        Double I_ESP_K4x3z_S_M3_vrr = PAZ*I_ESP_I4x2z_S_M3_vrr-PRZ*I_ESP_I4x2z_S_M4_vrr+2*oned2z*I_ESP_H4xz_S_M3_vrr-2*oned2z*I_ESP_H4xz_S_M4_vrr;
        Double I_ESP_K3x4y_S_M3_vrr = PAX*I_ESP_I2x4y_S_M3_vrr-PRX*I_ESP_I2x4y_S_M4_vrr+2*oned2z*I_ESP_Hx4y_S_M3_vrr-2*oned2z*I_ESP_Hx4y_S_M4_vrr;
        Double I_ESP_K3x4z_S_M3_vrr = PAX*I_ESP_I2x4z_S_M3_vrr-PRX*I_ESP_I2x4z_S_M4_vrr+2*oned2z*I_ESP_Hx4z_S_M3_vrr-2*oned2z*I_ESP_Hx4z_S_M4_vrr;
        Double I_ESP_K2x5y_S_M3_vrr = PAX*I_ESP_Ix5y_S_M3_vrr-PRX*I_ESP_Ix5y_S_M4_vrr+oned2z*I_ESP_H5y_S_M3_vrr-oned2z*I_ESP_H5y_S_M4_vrr;
        Double I_ESP_K2x5z_S_M3_vrr = PAX*I_ESP_Ix5z_S_M3_vrr-PRX*I_ESP_Ix5z_S_M4_vrr+oned2z*I_ESP_H5z_S_M3_vrr-oned2z*I_ESP_H5z_S_M4_vrr;
        Double I_ESP_Kx6y_S_M3_vrr = PAX*I_ESP_I6y_S_M3_vrr-PRX*I_ESP_I6y_S_M4_vrr;
        Double I_ESP_Kx6z_S_M3_vrr = PAX*I_ESP_I6z_S_M3_vrr-PRX*I_ESP_I6z_S_M4_vrr;
        Double I_ESP_K7y_S_M3_vrr = PAY*I_ESP_I6y_S_M3_vrr-PRY*I_ESP_I6y_S_M4_vrr+6*oned2z*I_ESP_H5y_S_M3_vrr-6*oned2z*I_ESP_H5y_S_M4_vrr;
        Double I_ESP_K6yz_S_M3_vrr = PAZ*I_ESP_I6y_S_M3_vrr-PRZ*I_ESP_I6y_S_M4_vrr;
        Double I_ESP_K5y2z_S_M3_vrr = PAZ*I_ESP_I5yz_S_M3_vrr-PRZ*I_ESP_I5yz_S_M4_vrr+oned2z*I_ESP_H5y_S_M3_vrr-oned2z*I_ESP_H5y_S_M4_vrr;
        Double I_ESP_K4y3z_S_M3_vrr = PAZ*I_ESP_I4y2z_S_M3_vrr-PRZ*I_ESP_I4y2z_S_M4_vrr+2*oned2z*I_ESP_H4yz_S_M3_vrr-2*oned2z*I_ESP_H4yz_S_M4_vrr;
        Double I_ESP_K2y5z_S_M3_vrr = PAY*I_ESP_Iy5z_S_M3_vrr-PRY*I_ESP_Iy5z_S_M4_vrr+oned2z*I_ESP_H5z_S_M3_vrr-oned2z*I_ESP_H5z_S_M4_vrr;
        Double I_ESP_Ky6z_S_M3_vrr = PAY*I_ESP_I6z_S_M3_vrr-PRY*I_ESP_I6z_S_M4_vrr;
        Double I_ESP_K7z_S_M3_vrr = PAZ*I_ESP_I6z_S_M3_vrr-PRZ*I_ESP_I6z_S_M4_vrr+6*oned2z*I_ESP_H5z_S_M3_vrr-6*oned2z*I_ESP_H5z_S_M4_vrr;

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
         * shell quartet name: SQ_ESP_L_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 18 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_K_S_M2
         * RHS shell quartet name: SQ_ESP_K_S_M3
         * RHS shell quartet name: SQ_ESP_I_S_M2
         * RHS shell quartet name: SQ_ESP_I_S_M3
         ************************************************************/
        Double I_ESP_L8x_S_M2_vrr = PAX*I_ESP_K7x_S_M2_vrr-PRX*I_ESP_K7x_S_M3_vrr+7*oned2z*I_ESP_I6x_S_M2_vrr-7*oned2z*I_ESP_I6x_S_M3_vrr;
        Double I_ESP_L7xy_S_M2_vrr = PAY*I_ESP_K7x_S_M2_vrr-PRY*I_ESP_K7x_S_M3_vrr;
        Double I_ESP_L7xz_S_M2_vrr = PAZ*I_ESP_K7x_S_M2_vrr-PRZ*I_ESP_K7x_S_M3_vrr;
        Double I_ESP_L6x2y_S_M2_vrr = PAY*I_ESP_K6xy_S_M2_vrr-PRY*I_ESP_K6xy_S_M3_vrr+oned2z*I_ESP_I6x_S_M2_vrr-oned2z*I_ESP_I6x_S_M3_vrr;
        Double I_ESP_L6x2z_S_M2_vrr = PAZ*I_ESP_K6xz_S_M2_vrr-PRZ*I_ESP_K6xz_S_M3_vrr+oned2z*I_ESP_I6x_S_M2_vrr-oned2z*I_ESP_I6x_S_M3_vrr;
        Double I_ESP_L5x3y_S_M2_vrr = PAY*I_ESP_K5x2y_S_M2_vrr-PRY*I_ESP_K5x2y_S_M3_vrr+2*oned2z*I_ESP_I5xy_S_M2_vrr-2*oned2z*I_ESP_I5xy_S_M3_vrr;
        Double I_ESP_L5x3z_S_M2_vrr = PAZ*I_ESP_K5x2z_S_M2_vrr-PRZ*I_ESP_K5x2z_S_M3_vrr+2*oned2z*I_ESP_I5xz_S_M2_vrr-2*oned2z*I_ESP_I5xz_S_M3_vrr;
        Double I_ESP_L4x4y_S_M2_vrr = PAY*I_ESP_K4x3y_S_M2_vrr-PRY*I_ESP_K4x3y_S_M3_vrr+3*oned2z*I_ESP_I4x2y_S_M2_vrr-3*oned2z*I_ESP_I4x2y_S_M3_vrr;
        Double I_ESP_L4x3yz_S_M2_vrr = PAZ*I_ESP_K4x3y_S_M2_vrr-PRZ*I_ESP_K4x3y_S_M3_vrr;
        Double I_ESP_L4x4z_S_M2_vrr = PAZ*I_ESP_K4x3z_S_M2_vrr-PRZ*I_ESP_K4x3z_S_M3_vrr+3*oned2z*I_ESP_I4x2z_S_M2_vrr-3*oned2z*I_ESP_I4x2z_S_M3_vrr;
        Double I_ESP_L3x5y_S_M2_vrr = PAX*I_ESP_K2x5y_S_M2_vrr-PRX*I_ESP_K2x5y_S_M3_vrr+2*oned2z*I_ESP_Ix5y_S_M2_vrr-2*oned2z*I_ESP_Ix5y_S_M3_vrr;
        Double I_ESP_L3x4yz_S_M2_vrr = PAZ*I_ESP_K3x4y_S_M2_vrr-PRZ*I_ESP_K3x4y_S_M3_vrr;
        Double I_ESP_L3xy4z_S_M2_vrr = PAY*I_ESP_K3x4z_S_M2_vrr-PRY*I_ESP_K3x4z_S_M3_vrr;
        Double I_ESP_L3x5z_S_M2_vrr = PAX*I_ESP_K2x5z_S_M2_vrr-PRX*I_ESP_K2x5z_S_M3_vrr+2*oned2z*I_ESP_Ix5z_S_M2_vrr-2*oned2z*I_ESP_Ix5z_S_M3_vrr;
        Double I_ESP_L2x6y_S_M2_vrr = PAX*I_ESP_Kx6y_S_M2_vrr-PRX*I_ESP_Kx6y_S_M3_vrr+oned2z*I_ESP_I6y_S_M2_vrr-oned2z*I_ESP_I6y_S_M3_vrr;
        Double I_ESP_L2x6z_S_M2_vrr = PAX*I_ESP_Kx6z_S_M2_vrr-PRX*I_ESP_Kx6z_S_M3_vrr+oned2z*I_ESP_I6z_S_M2_vrr-oned2z*I_ESP_I6z_S_M3_vrr;
        Double I_ESP_Lx7y_S_M2_vrr = PAX*I_ESP_K7y_S_M2_vrr-PRX*I_ESP_K7y_S_M3_vrr;
        Double I_ESP_Lx7z_S_M2_vrr = PAX*I_ESP_K7z_S_M2_vrr-PRX*I_ESP_K7z_S_M3_vrr;
        Double I_ESP_L8y_S_M2_vrr = PAY*I_ESP_K7y_S_M2_vrr-PRY*I_ESP_K7y_S_M3_vrr+7*oned2z*I_ESP_I6y_S_M2_vrr-7*oned2z*I_ESP_I6y_S_M3_vrr;
        Double I_ESP_L7yz_S_M2_vrr = PAZ*I_ESP_K7y_S_M2_vrr-PRZ*I_ESP_K7y_S_M3_vrr;
        Double I_ESP_L6y2z_S_M2_vrr = PAZ*I_ESP_K6yz_S_M2_vrr-PRZ*I_ESP_K6yz_S_M3_vrr+oned2z*I_ESP_I6y_S_M2_vrr-oned2z*I_ESP_I6y_S_M3_vrr;
        Double I_ESP_L5y3z_S_M2_vrr = PAZ*I_ESP_K5y2z_S_M2_vrr-PRZ*I_ESP_K5y2z_S_M3_vrr+2*oned2z*I_ESP_I5yz_S_M2_vrr-2*oned2z*I_ESP_I5yz_S_M3_vrr;
        Double I_ESP_L4y4z_S_M2_vrr = PAZ*I_ESP_K4y3z_S_M2_vrr-PRZ*I_ESP_K4y3z_S_M3_vrr+3*oned2z*I_ESP_I4y2z_S_M2_vrr-3*oned2z*I_ESP_I4y2z_S_M3_vrr;
        Double I_ESP_L3y5z_S_M2_vrr = PAY*I_ESP_K2y5z_S_M2_vrr-PRY*I_ESP_K2y5z_S_M3_vrr+2*oned2z*I_ESP_Iy5z_S_M2_vrr-2*oned2z*I_ESP_Iy5z_S_M3_vrr;
        Double I_ESP_L2y6z_S_M2_vrr = PAY*I_ESP_Ky6z_S_M2_vrr-PRY*I_ESP_Ky6z_S_M3_vrr+oned2z*I_ESP_I6z_S_M2_vrr-oned2z*I_ESP_I6z_S_M3_vrr;
        Double I_ESP_Ly7z_S_M2_vrr = PAY*I_ESP_K7z_S_M2_vrr-PRY*I_ESP_K7z_S_M3_vrr;
        Double I_ESP_L8z_S_M2_vrr = PAZ*I_ESP_K7z_S_M2_vrr-PRZ*I_ESP_K7z_S_M3_vrr+7*oned2z*I_ESP_I6z_S_M2_vrr-7*oned2z*I_ESP_I6z_S_M3_vrr;

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
         * shell quartet name: SQ_ESP_M_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 13 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_L_S_M1
         * RHS shell quartet name: SQ_ESP_L_S_M2
         * RHS shell quartet name: SQ_ESP_K_S_M1
         * RHS shell quartet name: SQ_ESP_K_S_M2
         ************************************************************/
        Double I_ESP_M9x_S_M1_vrr = PAX*I_ESP_L8x_S_M1_vrr-PRX*I_ESP_L8x_S_M2_vrr+8*oned2z*I_ESP_K7x_S_M1_vrr-8*oned2z*I_ESP_K7x_S_M2_vrr;
        Double I_ESP_M8xy_S_M1_vrr = PAY*I_ESP_L8x_S_M1_vrr-PRY*I_ESP_L8x_S_M2_vrr;
        Double I_ESP_M8xz_S_M1_vrr = PAZ*I_ESP_L8x_S_M1_vrr-PRZ*I_ESP_L8x_S_M2_vrr;
        Double I_ESP_M7x2y_S_M1_vrr = PAY*I_ESP_L7xy_S_M1_vrr-PRY*I_ESP_L7xy_S_M2_vrr+oned2z*I_ESP_K7x_S_M1_vrr-oned2z*I_ESP_K7x_S_M2_vrr;
        Double I_ESP_M7x2z_S_M1_vrr = PAZ*I_ESP_L7xz_S_M1_vrr-PRZ*I_ESP_L7xz_S_M2_vrr+oned2z*I_ESP_K7x_S_M1_vrr-oned2z*I_ESP_K7x_S_M2_vrr;
        Double I_ESP_M6x3y_S_M1_vrr = PAY*I_ESP_L6x2y_S_M1_vrr-PRY*I_ESP_L6x2y_S_M2_vrr+2*oned2z*I_ESP_K6xy_S_M1_vrr-2*oned2z*I_ESP_K6xy_S_M2_vrr;
        Double I_ESP_M6x2yz_S_M1_vrr = PAZ*I_ESP_L6x2y_S_M1_vrr-PRZ*I_ESP_L6x2y_S_M2_vrr;
        Double I_ESP_M6x3z_S_M1_vrr = PAZ*I_ESP_L6x2z_S_M1_vrr-PRZ*I_ESP_L6x2z_S_M2_vrr+2*oned2z*I_ESP_K6xz_S_M1_vrr-2*oned2z*I_ESP_K6xz_S_M2_vrr;
        Double I_ESP_M5x4y_S_M1_vrr = PAY*I_ESP_L5x3y_S_M1_vrr-PRY*I_ESP_L5x3y_S_M2_vrr+3*oned2z*I_ESP_K5x2y_S_M1_vrr-3*oned2z*I_ESP_K5x2y_S_M2_vrr;
        Double I_ESP_M5x3yz_S_M1_vrr = PAZ*I_ESP_L5x3y_S_M1_vrr-PRZ*I_ESP_L5x3y_S_M2_vrr;
        Double I_ESP_M5xy3z_S_M1_vrr = PAY*I_ESP_L5x3z_S_M1_vrr-PRY*I_ESP_L5x3z_S_M2_vrr;
        Double I_ESP_M5x4z_S_M1_vrr = PAZ*I_ESP_L5x3z_S_M1_vrr-PRZ*I_ESP_L5x3z_S_M2_vrr+3*oned2z*I_ESP_K5x2z_S_M1_vrr-3*oned2z*I_ESP_K5x2z_S_M2_vrr;
        Double I_ESP_M4x5y_S_M1_vrr = PAX*I_ESP_L3x5y_S_M1_vrr-PRX*I_ESP_L3x5y_S_M2_vrr+3*oned2z*I_ESP_K2x5y_S_M1_vrr-3*oned2z*I_ESP_K2x5y_S_M2_vrr;
        Double I_ESP_M4x4yz_S_M1_vrr = PAZ*I_ESP_L4x4y_S_M1_vrr-PRZ*I_ESP_L4x4y_S_M2_vrr;
        Double I_ESP_M4x3y2z_S_M1_vrr = PAZ*I_ESP_L4x3yz_S_M1_vrr-PRZ*I_ESP_L4x3yz_S_M2_vrr+oned2z*I_ESP_K4x3y_S_M1_vrr-oned2z*I_ESP_K4x3y_S_M2_vrr;
        Double I_ESP_M4xy4z_S_M1_vrr = PAY*I_ESP_L4x4z_S_M1_vrr-PRY*I_ESP_L4x4z_S_M2_vrr;
        Double I_ESP_M4x5z_S_M1_vrr = PAX*I_ESP_L3x5z_S_M1_vrr-PRX*I_ESP_L3x5z_S_M2_vrr+3*oned2z*I_ESP_K2x5z_S_M1_vrr-3*oned2z*I_ESP_K2x5z_S_M2_vrr;
        Double I_ESP_M3x6y_S_M1_vrr = PAX*I_ESP_L2x6y_S_M1_vrr-PRX*I_ESP_L2x6y_S_M2_vrr+2*oned2z*I_ESP_Kx6y_S_M1_vrr-2*oned2z*I_ESP_Kx6y_S_M2_vrr;
        Double I_ESP_M3x5yz_S_M1_vrr = PAZ*I_ESP_L3x5y_S_M1_vrr-PRZ*I_ESP_L3x5y_S_M2_vrr;
        Double I_ESP_M3x4y2z_S_M1_vrr = PAZ*I_ESP_L3x4yz_S_M1_vrr-PRZ*I_ESP_L3x4yz_S_M2_vrr+oned2z*I_ESP_K3x4y_S_M1_vrr-oned2z*I_ESP_K3x4y_S_M2_vrr;
        Double I_ESP_M3x2y4z_S_M1_vrr = PAY*I_ESP_L3xy4z_S_M1_vrr-PRY*I_ESP_L3xy4z_S_M2_vrr+oned2z*I_ESP_K3x4z_S_M1_vrr-oned2z*I_ESP_K3x4z_S_M2_vrr;
        Double I_ESP_M3xy5z_S_M1_vrr = PAY*I_ESP_L3x5z_S_M1_vrr-PRY*I_ESP_L3x5z_S_M2_vrr;
        Double I_ESP_M3x6z_S_M1_vrr = PAX*I_ESP_L2x6z_S_M1_vrr-PRX*I_ESP_L2x6z_S_M2_vrr+2*oned2z*I_ESP_Kx6z_S_M1_vrr-2*oned2z*I_ESP_Kx6z_S_M2_vrr;
        Double I_ESP_M2x7y_S_M1_vrr = PAX*I_ESP_Lx7y_S_M1_vrr-PRX*I_ESP_Lx7y_S_M2_vrr+oned2z*I_ESP_K7y_S_M1_vrr-oned2z*I_ESP_K7y_S_M2_vrr;
        Double I_ESP_M2x6yz_S_M1_vrr = PAZ*I_ESP_L2x6y_S_M1_vrr-PRZ*I_ESP_L2x6y_S_M2_vrr;
        Double I_ESP_M2xy6z_S_M1_vrr = PAY*I_ESP_L2x6z_S_M1_vrr-PRY*I_ESP_L2x6z_S_M2_vrr;
        Double I_ESP_M2x7z_S_M1_vrr = PAX*I_ESP_Lx7z_S_M1_vrr-PRX*I_ESP_Lx7z_S_M2_vrr+oned2z*I_ESP_K7z_S_M1_vrr-oned2z*I_ESP_K7z_S_M2_vrr;
        Double I_ESP_Mx8y_S_M1_vrr = PAX*I_ESP_L8y_S_M1_vrr-PRX*I_ESP_L8y_S_M2_vrr;
        Double I_ESP_Mx5y3z_S_M1_vrr = PAX*I_ESP_L5y3z_S_M1_vrr-PRX*I_ESP_L5y3z_S_M2_vrr;
        Double I_ESP_Mx4y4z_S_M1_vrr = PAX*I_ESP_L4y4z_S_M1_vrr-PRX*I_ESP_L4y4z_S_M2_vrr;
        Double I_ESP_Mx3y5z_S_M1_vrr = PAX*I_ESP_L3y5z_S_M1_vrr-PRX*I_ESP_L3y5z_S_M2_vrr;
        Double I_ESP_Mx8z_S_M1_vrr = PAX*I_ESP_L8z_S_M1_vrr-PRX*I_ESP_L8z_S_M2_vrr;
        Double I_ESP_M9y_S_M1_vrr = PAY*I_ESP_L8y_S_M1_vrr-PRY*I_ESP_L8y_S_M2_vrr+8*oned2z*I_ESP_K7y_S_M1_vrr-8*oned2z*I_ESP_K7y_S_M2_vrr;
        Double I_ESP_M8yz_S_M1_vrr = PAZ*I_ESP_L8y_S_M1_vrr-PRZ*I_ESP_L8y_S_M2_vrr;
        Double I_ESP_M7y2z_S_M1_vrr = PAZ*I_ESP_L7yz_S_M1_vrr-PRZ*I_ESP_L7yz_S_M2_vrr+oned2z*I_ESP_K7y_S_M1_vrr-oned2z*I_ESP_K7y_S_M2_vrr;
        Double I_ESP_M6y3z_S_M1_vrr = PAZ*I_ESP_L6y2z_S_M1_vrr-PRZ*I_ESP_L6y2z_S_M2_vrr+2*oned2z*I_ESP_K6yz_S_M1_vrr-2*oned2z*I_ESP_K6yz_S_M2_vrr;
        Double I_ESP_M5y4z_S_M1_vrr = PAZ*I_ESP_L5y3z_S_M1_vrr-PRZ*I_ESP_L5y3z_S_M2_vrr+3*oned2z*I_ESP_K5y2z_S_M1_vrr-3*oned2z*I_ESP_K5y2z_S_M2_vrr;
        Double I_ESP_M4y5z_S_M1_vrr = PAY*I_ESP_L3y5z_S_M1_vrr-PRY*I_ESP_L3y5z_S_M2_vrr+3*oned2z*I_ESP_K2y5z_S_M1_vrr-3*oned2z*I_ESP_K2y5z_S_M2_vrr;
        Double I_ESP_M3y6z_S_M1_vrr = PAY*I_ESP_L2y6z_S_M1_vrr-PRY*I_ESP_L2y6z_S_M2_vrr+2*oned2z*I_ESP_Ky6z_S_M1_vrr-2*oned2z*I_ESP_Ky6z_S_M2_vrr;
        Double I_ESP_M2y7z_S_M1_vrr = PAY*I_ESP_Ly7z_S_M1_vrr-PRY*I_ESP_Ly7z_S_M2_vrr+oned2z*I_ESP_K7z_S_M1_vrr-oned2z*I_ESP_K7z_S_M2_vrr;
        Double I_ESP_My8z_S_M1_vrr = PAY*I_ESP_L8z_S_M1_vrr-PRY*I_ESP_L8z_S_M2_vrr;
        Double I_ESP_M9z_S_M1_vrr = PAZ*I_ESP_L8z_S_M1_vrr-PRZ*I_ESP_L8z_S_M2_vrr+8*oned2z*I_ESP_K7z_S_M1_vrr-8*oned2z*I_ESP_K7z_S_M2_vrr;

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
         * shell quartet name: SQ_ESP_N_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_M_S
         * RHS shell quartet name: SQ_ESP_M_S_M1
         * RHS shell quartet name: SQ_ESP_L_S
         * RHS shell quartet name: SQ_ESP_L_S_M1
         ************************************************************/
        Double I_ESP_N10x_S_vrr = PAX*I_ESP_M9x_S_vrr-PRX*I_ESP_M9x_S_M1_vrr+9*oned2z*I_ESP_L8x_S_vrr-9*oned2z*I_ESP_L8x_S_M1_vrr;
        Double I_ESP_N9xy_S_vrr = PAY*I_ESP_M9x_S_vrr-PRY*I_ESP_M9x_S_M1_vrr;
        Double I_ESP_N9xz_S_vrr = PAZ*I_ESP_M9x_S_vrr-PRZ*I_ESP_M9x_S_M1_vrr;
        Double I_ESP_N8x2y_S_vrr = PAY*I_ESP_M8xy_S_vrr-PRY*I_ESP_M8xy_S_M1_vrr+oned2z*I_ESP_L8x_S_vrr-oned2z*I_ESP_L8x_S_M1_vrr;
        Double I_ESP_N8xyz_S_vrr = PAZ*I_ESP_M8xy_S_vrr-PRZ*I_ESP_M8xy_S_M1_vrr;
        Double I_ESP_N8x2z_S_vrr = PAZ*I_ESP_M8xz_S_vrr-PRZ*I_ESP_M8xz_S_M1_vrr+oned2z*I_ESP_L8x_S_vrr-oned2z*I_ESP_L8x_S_M1_vrr;
        Double I_ESP_N7x3y_S_vrr = PAY*I_ESP_M7x2y_S_vrr-PRY*I_ESP_M7x2y_S_M1_vrr+2*oned2z*I_ESP_L7xy_S_vrr-2*oned2z*I_ESP_L7xy_S_M1_vrr;
        Double I_ESP_N7x2yz_S_vrr = PAZ*I_ESP_M7x2y_S_vrr-PRZ*I_ESP_M7x2y_S_M1_vrr;
        Double I_ESP_N7xy2z_S_vrr = PAY*I_ESP_M7x2z_S_vrr-PRY*I_ESP_M7x2z_S_M1_vrr;
        Double I_ESP_N7x3z_S_vrr = PAZ*I_ESP_M7x2z_S_vrr-PRZ*I_ESP_M7x2z_S_M1_vrr+2*oned2z*I_ESP_L7xz_S_vrr-2*oned2z*I_ESP_L7xz_S_M1_vrr;
        Double I_ESP_N6x4y_S_vrr = PAY*I_ESP_M6x3y_S_vrr-PRY*I_ESP_M6x3y_S_M1_vrr+3*oned2z*I_ESP_L6x2y_S_vrr-3*oned2z*I_ESP_L6x2y_S_M1_vrr;
        Double I_ESP_N6x3yz_S_vrr = PAZ*I_ESP_M6x3y_S_vrr-PRZ*I_ESP_M6x3y_S_M1_vrr;
        Double I_ESP_N6x2y2z_S_vrr = PAZ*I_ESP_M6x2yz_S_vrr-PRZ*I_ESP_M6x2yz_S_M1_vrr+oned2z*I_ESP_L6x2y_S_vrr-oned2z*I_ESP_L6x2y_S_M1_vrr;
        Double I_ESP_N6xy3z_S_vrr = PAY*I_ESP_M6x3z_S_vrr-PRY*I_ESP_M6x3z_S_M1_vrr;
        Double I_ESP_N6x4z_S_vrr = PAZ*I_ESP_M6x3z_S_vrr-PRZ*I_ESP_M6x3z_S_M1_vrr+3*oned2z*I_ESP_L6x2z_S_vrr-3*oned2z*I_ESP_L6x2z_S_M1_vrr;
        Double I_ESP_N5x5y_S_vrr = PAY*I_ESP_M5x4y_S_vrr-PRY*I_ESP_M5x4y_S_M1_vrr+4*oned2z*I_ESP_L5x3y_S_vrr-4*oned2z*I_ESP_L5x3y_S_M1_vrr;
        Double I_ESP_N5x4yz_S_vrr = PAZ*I_ESP_M5x4y_S_vrr-PRZ*I_ESP_M5x4y_S_M1_vrr;
        Double I_ESP_N5x3y2z_S_vrr = PAZ*I_ESP_M5x3yz_S_vrr-PRZ*I_ESP_M5x3yz_S_M1_vrr+oned2z*I_ESP_L5x3y_S_vrr-oned2z*I_ESP_L5x3y_S_M1_vrr;
        Double I_ESP_N5x2y3z_S_vrr = PAY*I_ESP_M5xy3z_S_vrr-PRY*I_ESP_M5xy3z_S_M1_vrr+oned2z*I_ESP_L5x3z_S_vrr-oned2z*I_ESP_L5x3z_S_M1_vrr;
        Double I_ESP_N5xy4z_S_vrr = PAY*I_ESP_M5x4z_S_vrr-PRY*I_ESP_M5x4z_S_M1_vrr;
        Double I_ESP_N5x5z_S_vrr = PAZ*I_ESP_M5x4z_S_vrr-PRZ*I_ESP_M5x4z_S_M1_vrr+4*oned2z*I_ESP_L5x3z_S_vrr-4*oned2z*I_ESP_L5x3z_S_M1_vrr;
        Double I_ESP_N4x6y_S_vrr = PAX*I_ESP_M3x6y_S_vrr-PRX*I_ESP_M3x6y_S_M1_vrr+3*oned2z*I_ESP_L2x6y_S_vrr-3*oned2z*I_ESP_L2x6y_S_M1_vrr;
        Double I_ESP_N4x5yz_S_vrr = PAZ*I_ESP_M4x5y_S_vrr-PRZ*I_ESP_M4x5y_S_M1_vrr;
        Double I_ESP_N4x4y2z_S_vrr = PAZ*I_ESP_M4x4yz_S_vrr-PRZ*I_ESP_M4x4yz_S_M1_vrr+oned2z*I_ESP_L4x4y_S_vrr-oned2z*I_ESP_L4x4y_S_M1_vrr;
        Double I_ESP_N4x3y3z_S_vrr = PAZ*I_ESP_M4x3y2z_S_vrr-PRZ*I_ESP_M4x3y2z_S_M1_vrr+2*oned2z*I_ESP_L4x3yz_S_vrr-2*oned2z*I_ESP_L4x3yz_S_M1_vrr;
        Double I_ESP_N4x2y4z_S_vrr = PAY*I_ESP_M4xy4z_S_vrr-PRY*I_ESP_M4xy4z_S_M1_vrr+oned2z*I_ESP_L4x4z_S_vrr-oned2z*I_ESP_L4x4z_S_M1_vrr;
        Double I_ESP_N4xy5z_S_vrr = PAY*I_ESP_M4x5z_S_vrr-PRY*I_ESP_M4x5z_S_M1_vrr;
        Double I_ESP_N4x6z_S_vrr = PAX*I_ESP_M3x6z_S_vrr-PRX*I_ESP_M3x6z_S_M1_vrr+3*oned2z*I_ESP_L2x6z_S_vrr-3*oned2z*I_ESP_L2x6z_S_M1_vrr;
        Double I_ESP_N3x7y_S_vrr = PAX*I_ESP_M2x7y_S_vrr-PRX*I_ESP_M2x7y_S_M1_vrr+2*oned2z*I_ESP_Lx7y_S_vrr-2*oned2z*I_ESP_Lx7y_S_M1_vrr;
        Double I_ESP_N3x6yz_S_vrr = PAZ*I_ESP_M3x6y_S_vrr-PRZ*I_ESP_M3x6y_S_M1_vrr;
        Double I_ESP_N3x5y2z_S_vrr = PAZ*I_ESP_M3x5yz_S_vrr-PRZ*I_ESP_M3x5yz_S_M1_vrr+oned2z*I_ESP_L3x5y_S_vrr-oned2z*I_ESP_L3x5y_S_M1_vrr;
        Double I_ESP_N3x4y3z_S_vrr = PAZ*I_ESP_M3x4y2z_S_vrr-PRZ*I_ESP_M3x4y2z_S_M1_vrr+2*oned2z*I_ESP_L3x4yz_S_vrr-2*oned2z*I_ESP_L3x4yz_S_M1_vrr;
        Double I_ESP_N3x3y4z_S_vrr = PAY*I_ESP_M3x2y4z_S_vrr-PRY*I_ESP_M3x2y4z_S_M1_vrr+2*oned2z*I_ESP_L3xy4z_S_vrr-2*oned2z*I_ESP_L3xy4z_S_M1_vrr;
        Double I_ESP_N3x2y5z_S_vrr = PAY*I_ESP_M3xy5z_S_vrr-PRY*I_ESP_M3xy5z_S_M1_vrr+oned2z*I_ESP_L3x5z_S_vrr-oned2z*I_ESP_L3x5z_S_M1_vrr;
        Double I_ESP_N3xy6z_S_vrr = PAY*I_ESP_M3x6z_S_vrr-PRY*I_ESP_M3x6z_S_M1_vrr;
        Double I_ESP_N3x7z_S_vrr = PAX*I_ESP_M2x7z_S_vrr-PRX*I_ESP_M2x7z_S_M1_vrr+2*oned2z*I_ESP_Lx7z_S_vrr-2*oned2z*I_ESP_Lx7z_S_M1_vrr;
        Double I_ESP_N2x8y_S_vrr = PAX*I_ESP_Mx8y_S_vrr-PRX*I_ESP_Mx8y_S_M1_vrr+oned2z*I_ESP_L8y_S_vrr-oned2z*I_ESP_L8y_S_M1_vrr;
        Double I_ESP_N2x7yz_S_vrr = PAZ*I_ESP_M2x7y_S_vrr-PRZ*I_ESP_M2x7y_S_M1_vrr;
        Double I_ESP_N2x6y2z_S_vrr = PAZ*I_ESP_M2x6yz_S_vrr-PRZ*I_ESP_M2x6yz_S_M1_vrr+oned2z*I_ESP_L2x6y_S_vrr-oned2z*I_ESP_L2x6y_S_M1_vrr;
        Double I_ESP_N2x5y3z_S_vrr = PAX*I_ESP_Mx5y3z_S_vrr-PRX*I_ESP_Mx5y3z_S_M1_vrr+oned2z*I_ESP_L5y3z_S_vrr-oned2z*I_ESP_L5y3z_S_M1_vrr;
        Double I_ESP_N2x4y4z_S_vrr = PAX*I_ESP_Mx4y4z_S_vrr-PRX*I_ESP_Mx4y4z_S_M1_vrr+oned2z*I_ESP_L4y4z_S_vrr-oned2z*I_ESP_L4y4z_S_M1_vrr;
        Double I_ESP_N2x3y5z_S_vrr = PAX*I_ESP_Mx3y5z_S_vrr-PRX*I_ESP_Mx3y5z_S_M1_vrr+oned2z*I_ESP_L3y5z_S_vrr-oned2z*I_ESP_L3y5z_S_M1_vrr;
        Double I_ESP_N2x2y6z_S_vrr = PAY*I_ESP_M2xy6z_S_vrr-PRY*I_ESP_M2xy6z_S_M1_vrr+oned2z*I_ESP_L2x6z_S_vrr-oned2z*I_ESP_L2x6z_S_M1_vrr;
        Double I_ESP_N2xy7z_S_vrr = PAY*I_ESP_M2x7z_S_vrr-PRY*I_ESP_M2x7z_S_M1_vrr;
        Double I_ESP_N2x8z_S_vrr = PAX*I_ESP_Mx8z_S_vrr-PRX*I_ESP_Mx8z_S_M1_vrr+oned2z*I_ESP_L8z_S_vrr-oned2z*I_ESP_L8z_S_M1_vrr;
        Double I_ESP_Nx9y_S_vrr = PAX*I_ESP_M9y_S_vrr-PRX*I_ESP_M9y_S_M1_vrr;
        Double I_ESP_Nx8yz_S_vrr = PAZ*I_ESP_Mx8y_S_vrr-PRZ*I_ESP_Mx8y_S_M1_vrr;
        Double I_ESP_Nx7y2z_S_vrr = PAX*I_ESP_M7y2z_S_vrr-PRX*I_ESP_M7y2z_S_M1_vrr;
        Double I_ESP_Nx6y3z_S_vrr = PAX*I_ESP_M6y3z_S_vrr-PRX*I_ESP_M6y3z_S_M1_vrr;
        Double I_ESP_Nx5y4z_S_vrr = PAX*I_ESP_M5y4z_S_vrr-PRX*I_ESP_M5y4z_S_M1_vrr;
        Double I_ESP_Nx4y5z_S_vrr = PAX*I_ESP_M4y5z_S_vrr-PRX*I_ESP_M4y5z_S_M1_vrr;
        Double I_ESP_Nx3y6z_S_vrr = PAX*I_ESP_M3y6z_S_vrr-PRX*I_ESP_M3y6z_S_M1_vrr;
        Double I_ESP_Nx2y7z_S_vrr = PAX*I_ESP_M2y7z_S_vrr-PRX*I_ESP_M2y7z_S_M1_vrr;
        Double I_ESP_Nxy8z_S_vrr = PAY*I_ESP_Mx8z_S_vrr-PRY*I_ESP_Mx8z_S_M1_vrr;
        Double I_ESP_Nx9z_S_vrr = PAX*I_ESP_M9z_S_vrr-PRX*I_ESP_M9z_S_M1_vrr;
        Double I_ESP_N10y_S_vrr = PAY*I_ESP_M9y_S_vrr-PRY*I_ESP_M9y_S_M1_vrr+9*oned2z*I_ESP_L8y_S_vrr-9*oned2z*I_ESP_L8y_S_M1_vrr;
        Double I_ESP_N9yz_S_vrr = PAZ*I_ESP_M9y_S_vrr-PRZ*I_ESP_M9y_S_M1_vrr;
        Double I_ESP_N8y2z_S_vrr = PAZ*I_ESP_M8yz_S_vrr-PRZ*I_ESP_M8yz_S_M1_vrr+oned2z*I_ESP_L8y_S_vrr-oned2z*I_ESP_L8y_S_M1_vrr;
        Double I_ESP_N7y3z_S_vrr = PAZ*I_ESP_M7y2z_S_vrr-PRZ*I_ESP_M7y2z_S_M1_vrr+2*oned2z*I_ESP_L7yz_S_vrr-2*oned2z*I_ESP_L7yz_S_M1_vrr;
        Double I_ESP_N6y4z_S_vrr = PAZ*I_ESP_M6y3z_S_vrr-PRZ*I_ESP_M6y3z_S_M1_vrr+3*oned2z*I_ESP_L6y2z_S_vrr-3*oned2z*I_ESP_L6y2z_S_M1_vrr;
        Double I_ESP_N5y5z_S_vrr = PAZ*I_ESP_M5y4z_S_vrr-PRZ*I_ESP_M5y4z_S_M1_vrr+4*oned2z*I_ESP_L5y3z_S_vrr-4*oned2z*I_ESP_L5y3z_S_M1_vrr;
        Double I_ESP_N4y6z_S_vrr = PAY*I_ESP_M3y6z_S_vrr-PRY*I_ESP_M3y6z_S_M1_vrr+3*oned2z*I_ESP_L2y6z_S_vrr-3*oned2z*I_ESP_L2y6z_S_M1_vrr;
        Double I_ESP_N3y7z_S_vrr = PAY*I_ESP_M2y7z_S_vrr-PRY*I_ESP_M2y7z_S_M1_vrr+2*oned2z*I_ESP_Ly7z_S_vrr-2*oned2z*I_ESP_Ly7z_S_M1_vrr;
        Double I_ESP_N2y8z_S_vrr = PAY*I_ESP_My8z_S_vrr-PRY*I_ESP_My8z_S_M1_vrr+oned2z*I_ESP_L8z_S_vrr-oned2z*I_ESP_L8z_S_M1_vrr;
        Double I_ESP_Ny9z_S_vrr = PAY*I_ESP_M9z_S_vrr-PRY*I_ESP_M9z_S_M1_vrr;
        Double I_ESP_N10z_S_vrr = PAZ*I_ESP_M9z_S_vrr-PRZ*I_ESP_M9z_S_M1_vrr+9*oned2z*I_ESP_L8z_S_vrr-9*oned2z*I_ESP_L8z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_N_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_N_S_aa_coefs = alpha*alpha;
        I_ESP_N10x_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N10x_S_vrr;
        I_ESP_N9xy_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N9xy_S_vrr;
        I_ESP_N9xz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N9xz_S_vrr;
        I_ESP_N8x2y_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N8x2y_S_vrr;
        I_ESP_N8xyz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N8xyz_S_vrr;
        I_ESP_N8x2z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N8x2z_S_vrr;
        I_ESP_N7x3y_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N7x3y_S_vrr;
        I_ESP_N7x2yz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N7x2yz_S_vrr;
        I_ESP_N7xy2z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N7xy2z_S_vrr;
        I_ESP_N7x3z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N7x3z_S_vrr;
        I_ESP_N6x4y_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N6x4y_S_vrr;
        I_ESP_N6x3yz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N6x3yz_S_vrr;
        I_ESP_N6x2y2z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N6x2y2z_S_vrr;
        I_ESP_N6xy3z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N6xy3z_S_vrr;
        I_ESP_N6x4z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N6x4z_S_vrr;
        I_ESP_N5x5y_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N5x5y_S_vrr;
        I_ESP_N5x4yz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N5x4yz_S_vrr;
        I_ESP_N5x3y2z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N5x3y2z_S_vrr;
        I_ESP_N5x2y3z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N5x2y3z_S_vrr;
        I_ESP_N5xy4z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N5xy4z_S_vrr;
        I_ESP_N5x5z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N5x5z_S_vrr;
        I_ESP_N4x6y_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N4x6y_S_vrr;
        I_ESP_N4x5yz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N4x5yz_S_vrr;
        I_ESP_N4x4y2z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N4x4y2z_S_vrr;
        I_ESP_N4x3y3z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N4x3y3z_S_vrr;
        I_ESP_N4x2y4z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N4x2y4z_S_vrr;
        I_ESP_N4xy5z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N4xy5z_S_vrr;
        I_ESP_N4x6z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N4x6z_S_vrr;
        I_ESP_N3x7y_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N3x7y_S_vrr;
        I_ESP_N3x6yz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N3x6yz_S_vrr;
        I_ESP_N3x5y2z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N3x5y2z_S_vrr;
        I_ESP_N3x4y3z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N3x4y3z_S_vrr;
        I_ESP_N3x3y4z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N3x3y4z_S_vrr;
        I_ESP_N3x2y5z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N3x2y5z_S_vrr;
        I_ESP_N3xy6z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N3xy6z_S_vrr;
        I_ESP_N3x7z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N3x7z_S_vrr;
        I_ESP_N2x8y_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2x8y_S_vrr;
        I_ESP_N2x7yz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2x7yz_S_vrr;
        I_ESP_N2x6y2z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2x6y2z_S_vrr;
        I_ESP_N2x5y3z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2x5y3z_S_vrr;
        I_ESP_N2x4y4z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2x4y4z_S_vrr;
        I_ESP_N2x3y5z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2x3y5z_S_vrr;
        I_ESP_N2x2y6z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2x2y6z_S_vrr;
        I_ESP_N2xy7z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2xy7z_S_vrr;
        I_ESP_N2x8z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2x8z_S_vrr;
        I_ESP_Nx9y_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nx9y_S_vrr;
        I_ESP_Nx8yz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nx8yz_S_vrr;
        I_ESP_Nx7y2z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nx7y2z_S_vrr;
        I_ESP_Nx6y3z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nx6y3z_S_vrr;
        I_ESP_Nx5y4z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nx5y4z_S_vrr;
        I_ESP_Nx4y5z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nx4y5z_S_vrr;
        I_ESP_Nx3y6z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nx3y6z_S_vrr;
        I_ESP_Nx2y7z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nx2y7z_S_vrr;
        I_ESP_Nxy8z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nxy8z_S_vrr;
        I_ESP_Nx9z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Nx9z_S_vrr;
        I_ESP_N10y_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N10y_S_vrr;
        I_ESP_N9yz_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N9yz_S_vrr;
        I_ESP_N8y2z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N8y2z_S_vrr;
        I_ESP_N7y3z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N7y3z_S_vrr;
        I_ESP_N6y4z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N6y4z_S_vrr;
        I_ESP_N5y5z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N5y5z_S_vrr;
        I_ESP_N4y6z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N4y6z_S_vrr;
        I_ESP_N3y7z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N3y7z_S_vrr;
        I_ESP_N2y8z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N2y8z_S_vrr;
        I_ESP_Ny9z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_Ny9z_S_vrr;
        I_ESP_N10z_S_aa += SQ_ESP_N_S_aa_coefs*I_ESP_N10z_S_vrr;

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
     * totally 18 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_F_P
     * RHS shell quartet name: SQ_ESP_D_P
     ************************************************************/
    Double I_ESP_D2x_D2x = I_ESP_F3x_Px+ABX*I_ESP_D2x_Px;
    Double I_ESP_Dxy_D2x = I_ESP_F2xy_Px+ABX*I_ESP_Dxy_Px;
    Double I_ESP_Dxz_D2x = I_ESP_F2xz_Px+ABX*I_ESP_Dxz_Px;
    Double I_ESP_D2y_D2x = I_ESP_Fx2y_Px+ABX*I_ESP_D2y_Px;
    Double I_ESP_Dyz_D2x = I_ESP_Fxyz_Px+ABX*I_ESP_Dyz_Px;
    Double I_ESP_D2z_D2x = I_ESP_Fx2z_Px+ABX*I_ESP_D2z_Px;
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
     * totally 3 integrals are omitted 
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
    Double I_ESP_G3xy_Pz = I_ESP_H3xyz_S+ABZ*I_ESP_G3xy_S;
    Double I_ESP_G3xz_Pz = I_ESP_H3x2z_S+ABZ*I_ESP_G3xz_S;
    Double I_ESP_G2x2y_Pz = I_ESP_H2x2yz_S+ABZ*I_ESP_G2x2y_S;
    Double I_ESP_G2xyz_Pz = I_ESP_H2xy2z_S+ABZ*I_ESP_G2xyz_S;
    Double I_ESP_G2x2z_Pz = I_ESP_H2x3z_S+ABZ*I_ESP_G2x2z_S;
    Double I_ESP_Gx3y_Pz = I_ESP_Hx3yz_S+ABZ*I_ESP_Gx3y_S;
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
     * shell quartet name: SQ_ESP_D_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 12 integrals are omitted 
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
    Double I_ESP_D2x_F3z = I_ESP_F2xz_D2z+ABZ*I_ESP_D2x_D2z;
    Double I_ESP_Dxy_F3z = I_ESP_Fxyz_D2z+ABZ*I_ESP_Dxy_D2z;
    Double I_ESP_Dxz_F3z = I_ESP_Fx2z_D2z+ABZ*I_ESP_Dxz_D2z;
    Double I_ESP_D2y_F3z = I_ESP_F2yz_D2z+ABZ*I_ESP_D2y_D2z;
    Double I_ESP_Dyz_F3z = I_ESP_Fy2z_D2z+ABZ*I_ESP_Dyz_D2z;
    Double I_ESP_D2z_F3z = I_ESP_F3z_D2z+ABZ*I_ESP_D2z_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_H_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 21 integrals are omitted 
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
    Double I_ESP_H3x2y_Py = I_ESP_I3x3y_S+ABY*I_ESP_H3x2y_S;
    Double I_ESP_H3xyz_Py = I_ESP_I3x2yz_S+ABY*I_ESP_H3xyz_S;
    Double I_ESP_H2x3y_Py = I_ESP_I2x4y_S+ABY*I_ESP_H2x3y_S;
    Double I_ESP_H2x2yz_Py = I_ESP_I2x3yz_S+ABY*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Py = I_ESP_I2x2y2z_S+ABY*I_ESP_H2xy2z_S;
    Double I_ESP_Hx4y_Py = I_ESP_Ix5y_S+ABY*I_ESP_Hx4y_S;
    Double I_ESP_Hx3yz_Py = I_ESP_Ix4yz_S+ABY*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Py = I_ESP_Ix3y2z_S+ABY*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Py = I_ESP_Ix2y3z_S+ABY*I_ESP_Hxy3z_S;
    Double I_ESP_H5y_Py = I_ESP_I6y_S+ABY*I_ESP_H5y_S;
    Double I_ESP_H4yz_Py = I_ESP_I5yz_S+ABY*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Py = I_ESP_I4y2z_S+ABY*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Py = I_ESP_I3y3z_S+ABY*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Py = I_ESP_I2y4z_S+ABY*I_ESP_Hy4z_S;
    Double I_ESP_H3xyz_Pz = I_ESP_I3xy2z_S+ABZ*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Pz = I_ESP_I3x3z_S+ABZ*I_ESP_H3x2z_S;
    Double I_ESP_H2x2yz_Pz = I_ESP_I2x2y2z_S+ABZ*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Pz = I_ESP_I2xy3z_S+ABZ*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Pz = I_ESP_I2x4z_S+ABZ*I_ESP_H2x3z_S;
    Double I_ESP_Hx3yz_Pz = I_ESP_Ix3y2z_S+ABZ*I_ESP_Hx3yz_S;
    Double I_ESP_Hx2y2z_Pz = I_ESP_Ix2y3z_S+ABZ*I_ESP_Hx2y2z_S;
    Double I_ESP_Hxy3z_Pz = I_ESP_Ixy4z_S+ABZ*I_ESP_Hxy3z_S;
    Double I_ESP_Hx4z_Pz = I_ESP_Ix5z_S+ABZ*I_ESP_Hx4z_S;
    Double I_ESP_H3y2z_Pz = I_ESP_I3y3z_S+ABZ*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Pz = I_ESP_I2y4z_S+ABZ*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Pz = I_ESP_Iy5z_S+ABZ*I_ESP_Hy4z_S;
    Double I_ESP_H5z_Pz = I_ESP_I6z_S+ABZ*I_ESP_H5z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_G_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 48 integrals are omitted 
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
    Double I_ESP_G3xy_D2z = I_ESP_H3xyz_Pz+ABZ*I_ESP_G3xy_Pz;
    Double I_ESP_G3xz_D2z = I_ESP_H3x2z_Pz+ABZ*I_ESP_G3xz_Pz;
    Double I_ESP_G2x2y_D2z = I_ESP_H2x2yz_Pz+ABZ*I_ESP_G2x2y_Pz;
    Double I_ESP_G2xyz_D2z = I_ESP_H2xy2z_Pz+ABZ*I_ESP_G2xyz_Pz;
    Double I_ESP_G2x2z_D2z = I_ESP_H2x3z_Pz+ABZ*I_ESP_G2x2z_Pz;
    Double I_ESP_Gx3y_D2z = I_ESP_Hx3yz_Pz+ABZ*I_ESP_Gx3y_Pz;
    Double I_ESP_Gx2yz_D2z = I_ESP_Hx2y2z_Pz+ABZ*I_ESP_Gx2yz_Pz;
    Double I_ESP_Gxy2z_D2z = I_ESP_Hxy3z_Pz+ABZ*I_ESP_Gxy2z_Pz;
    Double I_ESP_Gx3z_D2z = I_ESP_Hx4z_Pz+ABZ*I_ESP_Gx3z_Pz;
    Double I_ESP_G3yz_D2z = I_ESP_H3y2z_Pz+ABZ*I_ESP_G3yz_Pz;
    Double I_ESP_G2y2z_D2z = I_ESP_H2y3z_Pz+ABZ*I_ESP_G2y2z_Pz;
    Double I_ESP_Gy3z_D2z = I_ESP_Hy4z_Pz+ABZ*I_ESP_Gy3z_Pz;
    Double I_ESP_G4z_D2z = I_ESP_H5z_Pz+ABZ*I_ESP_G4z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_F_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 37 integrals are omitted 
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
    Double I_ESP_F2xy_F2xy = I_ESP_G2x2y_D2x+ABY*I_ESP_F2xy_D2x;
    Double I_ESP_F2xz_F2xy = I_ESP_G2xyz_D2x+ABY*I_ESP_F2xz_D2x;
    Double I_ESP_Fx2y_F2xy = I_ESP_Gx3y_D2x+ABY*I_ESP_Fx2y_D2x;
    Double I_ESP_Fxyz_F2xy = I_ESP_Gx2yz_D2x+ABY*I_ESP_Fxyz_D2x;
    Double I_ESP_Fx2z_F2xy = I_ESP_Gxy2z_D2x+ABY*I_ESP_Fx2z_D2x;
    Double I_ESP_F3y_F2xy = I_ESP_G4y_D2x+ABY*I_ESP_F3y_D2x;
    Double I_ESP_F2yz_F2xy = I_ESP_G3yz_D2x+ABY*I_ESP_F2yz_D2x;
    Double I_ESP_Fy2z_F2xy = I_ESP_G2y2z_D2x+ABY*I_ESP_Fy2z_D2x;
    Double I_ESP_F3z_F2xy = I_ESP_Gy3z_D2x+ABY*I_ESP_F3z_D2x;
    Double I_ESP_F2xz_F2xz = I_ESP_G2x2z_D2x+ABZ*I_ESP_F2xz_D2x;
    Double I_ESP_Fxyz_F2xz = I_ESP_Gxy2z_D2x+ABZ*I_ESP_Fxyz_D2x;
    Double I_ESP_Fx2z_F2xz = I_ESP_Gx3z_D2x+ABZ*I_ESP_Fx2z_D2x;
    Double I_ESP_F2yz_F2xz = I_ESP_G2y2z_D2x+ABZ*I_ESP_F2yz_D2x;
    Double I_ESP_Fy2z_F2xz = I_ESP_Gy3z_D2x+ABZ*I_ESP_Fy2z_D2x;
    Double I_ESP_F3z_F2xz = I_ESP_G4z_D2x+ABZ*I_ESP_F3z_D2x;
    Double I_ESP_F2xz_Fx2y = I_ESP_G3xz_D2y+ABX*I_ESP_F2xz_D2y;
    Double I_ESP_Fxyz_Fx2y = I_ESP_G2xyz_D2y+ABX*I_ESP_Fxyz_D2y;
    Double I_ESP_Fx2z_Fx2y = I_ESP_G2x2z_D2y+ABX*I_ESP_Fx2z_D2y;
    Double I_ESP_F2yz_Fx2y = I_ESP_Gx2yz_D2y+ABX*I_ESP_F2yz_D2y;
    Double I_ESP_Fy2z_Fx2y = I_ESP_Gxy2z_D2y+ABX*I_ESP_Fy2z_D2y;
    Double I_ESP_F3z_Fx2y = I_ESP_Gx3z_D2y+ABX*I_ESP_F3z_D2y;
    Double I_ESP_F2xy_Fx2z = I_ESP_G3xy_D2z+ABX*I_ESP_F2xy_D2z;
    Double I_ESP_Fx2y_Fx2z = I_ESP_G2x2y_D2z+ABX*I_ESP_Fx2y_D2z;
    Double I_ESP_Fxyz_Fx2z = I_ESP_G2xyz_D2z+ABX*I_ESP_Fxyz_D2z;
    Double I_ESP_F3y_Fx2z = I_ESP_Gx3y_D2z+ABX*I_ESP_F3y_D2z;
    Double I_ESP_F2yz_Fx2z = I_ESP_Gx2yz_D2z+ABX*I_ESP_F2yz_D2z;
    Double I_ESP_Fy2z_Fx2z = I_ESP_Gxy2z_D2z+ABX*I_ESP_Fy2z_D2z;
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
    Double I_ESP_F2xz_F2yz = I_ESP_G2x2z_D2y+ABZ*I_ESP_F2xz_D2y;
    Double I_ESP_Fxyz_F2yz = I_ESP_Gxy2z_D2y+ABZ*I_ESP_Fxyz_D2y;
    Double I_ESP_Fx2z_F2yz = I_ESP_Gx3z_D2y+ABZ*I_ESP_Fx2z_D2y;
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
     * shell quartet name: SQ_ESP_D_G
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_F_F
     * RHS shell quartet name: SQ_ESP_D_F
     ************************************************************/
    Double I_ESP_D2x_G4x = I_ESP_F3x_F3x+ABX*I_ESP_D2x_F3x;
    Double I_ESP_Dxy_G4x = I_ESP_F2xy_F3x+ABX*I_ESP_Dxy_F3x;
    Double I_ESP_Dxz_G4x = I_ESP_F2xz_F3x+ABX*I_ESP_Dxz_F3x;
    Double I_ESP_D2y_G4x = I_ESP_Fx2y_F3x+ABX*I_ESP_D2y_F3x;
    Double I_ESP_Dyz_G4x = I_ESP_Fxyz_F3x+ABX*I_ESP_Dyz_F3x;
    Double I_ESP_D2z_G4x = I_ESP_Fx2z_F3x+ABX*I_ESP_D2z_F3x;
    Double I_ESP_D2x_G3xy = I_ESP_F2xy_F3x+ABY*I_ESP_D2x_F3x;
    Double I_ESP_Dxy_G3xy = I_ESP_Fx2y_F3x+ABY*I_ESP_Dxy_F3x;
    Double I_ESP_Dxz_G3xy = I_ESP_Fxyz_F3x+ABY*I_ESP_Dxz_F3x;
    Double I_ESP_D2y_G3xy = I_ESP_F3y_F3x+ABY*I_ESP_D2y_F3x;
    Double I_ESP_Dyz_G3xy = I_ESP_F2yz_F3x+ABY*I_ESP_Dyz_F3x;
    Double I_ESP_D2z_G3xy = I_ESP_Fy2z_F3x+ABY*I_ESP_D2z_F3x;
    Double I_ESP_D2x_G3xz = I_ESP_F2xz_F3x+ABZ*I_ESP_D2x_F3x;
    Double I_ESP_Dxy_G3xz = I_ESP_Fxyz_F3x+ABZ*I_ESP_Dxy_F3x;
    Double I_ESP_Dxz_G3xz = I_ESP_Fx2z_F3x+ABZ*I_ESP_Dxz_F3x;
    Double I_ESP_D2y_G3xz = I_ESP_F2yz_F3x+ABZ*I_ESP_D2y_F3x;
    Double I_ESP_Dyz_G3xz = I_ESP_Fy2z_F3x+ABZ*I_ESP_Dyz_F3x;
    Double I_ESP_D2z_G3xz = I_ESP_F3z_F3x+ABZ*I_ESP_D2z_F3x;
    Double I_ESP_D2x_G2x2y = I_ESP_F2xy_F2xy+ABY*I_ESP_D2x_F2xy;
    Double I_ESP_Dxy_G2x2y = I_ESP_Fx2y_F2xy+ABY*I_ESP_Dxy_F2xy;
    Double I_ESP_Dxz_G2x2y = I_ESP_Fxyz_F2xy+ABY*I_ESP_Dxz_F2xy;
    Double I_ESP_D2y_G2x2y = I_ESP_F3y_F2xy+ABY*I_ESP_D2y_F2xy;
    Double I_ESP_Dyz_G2x2y = I_ESP_F2yz_F2xy+ABY*I_ESP_Dyz_F2xy;
    Double I_ESP_D2z_G2x2y = I_ESP_Fy2z_F2xy+ABY*I_ESP_D2z_F2xy;
    Double I_ESP_D2x_G2xyz = I_ESP_F2xz_F2xy+ABZ*I_ESP_D2x_F2xy;
    Double I_ESP_Dxy_G2xyz = I_ESP_Fxyz_F2xy+ABZ*I_ESP_Dxy_F2xy;
    Double I_ESP_Dxz_G2xyz = I_ESP_Fx2z_F2xy+ABZ*I_ESP_Dxz_F2xy;
    Double I_ESP_D2y_G2xyz = I_ESP_F2yz_F2xy+ABZ*I_ESP_D2y_F2xy;
    Double I_ESP_Dyz_G2xyz = I_ESP_Fy2z_F2xy+ABZ*I_ESP_Dyz_F2xy;
    Double I_ESP_D2z_G2xyz = I_ESP_F3z_F2xy+ABZ*I_ESP_D2z_F2xy;
    Double I_ESP_D2x_G2x2z = I_ESP_F2xz_F2xz+ABZ*I_ESP_D2x_F2xz;
    Double I_ESP_Dxy_G2x2z = I_ESP_Fxyz_F2xz+ABZ*I_ESP_Dxy_F2xz;
    Double I_ESP_Dxz_G2x2z = I_ESP_Fx2z_F2xz+ABZ*I_ESP_Dxz_F2xz;
    Double I_ESP_D2y_G2x2z = I_ESP_F2yz_F2xz+ABZ*I_ESP_D2y_F2xz;
    Double I_ESP_Dyz_G2x2z = I_ESP_Fy2z_F2xz+ABZ*I_ESP_Dyz_F2xz;
    Double I_ESP_D2z_G2x2z = I_ESP_F3z_F2xz+ABZ*I_ESP_D2z_F2xz;
    Double I_ESP_D2x_Gx3y = I_ESP_F3x_F3y+ABX*I_ESP_D2x_F3y;
    Double I_ESP_Dxy_Gx3y = I_ESP_F2xy_F3y+ABX*I_ESP_Dxy_F3y;
    Double I_ESP_Dxz_Gx3y = I_ESP_F2xz_F3y+ABX*I_ESP_Dxz_F3y;
    Double I_ESP_D2y_Gx3y = I_ESP_Fx2y_F3y+ABX*I_ESP_D2y_F3y;
    Double I_ESP_Dyz_Gx3y = I_ESP_Fxyz_F3y+ABX*I_ESP_Dyz_F3y;
    Double I_ESP_D2z_Gx3y = I_ESP_Fx2z_F3y+ABX*I_ESP_D2z_F3y;
    Double I_ESP_D2x_Gx2yz = I_ESP_F2xz_Fx2y+ABZ*I_ESP_D2x_Fx2y;
    Double I_ESP_Dxy_Gx2yz = I_ESP_Fxyz_Fx2y+ABZ*I_ESP_Dxy_Fx2y;
    Double I_ESP_Dxz_Gx2yz = I_ESP_Fx2z_Fx2y+ABZ*I_ESP_Dxz_Fx2y;
    Double I_ESP_D2y_Gx2yz = I_ESP_F2yz_Fx2y+ABZ*I_ESP_D2y_Fx2y;
    Double I_ESP_Dyz_Gx2yz = I_ESP_Fy2z_Fx2y+ABZ*I_ESP_Dyz_Fx2y;
    Double I_ESP_D2z_Gx2yz = I_ESP_F3z_Fx2y+ABZ*I_ESP_D2z_Fx2y;
    Double I_ESP_D2x_Gxy2z = I_ESP_F2xy_Fx2z+ABY*I_ESP_D2x_Fx2z;
    Double I_ESP_Dxy_Gxy2z = I_ESP_Fx2y_Fx2z+ABY*I_ESP_Dxy_Fx2z;
    Double I_ESP_Dxz_Gxy2z = I_ESP_Fxyz_Fx2z+ABY*I_ESP_Dxz_Fx2z;
    Double I_ESP_D2y_Gxy2z = I_ESP_F3y_Fx2z+ABY*I_ESP_D2y_Fx2z;
    Double I_ESP_Dyz_Gxy2z = I_ESP_F2yz_Fx2z+ABY*I_ESP_Dyz_Fx2z;
    Double I_ESP_D2z_Gxy2z = I_ESP_Fy2z_Fx2z+ABY*I_ESP_D2z_Fx2z;
    Double I_ESP_D2x_Gx3z = I_ESP_F3x_F3z+ABX*I_ESP_D2x_F3z;
    Double I_ESP_Dxy_Gx3z = I_ESP_F2xy_F3z+ABX*I_ESP_Dxy_F3z;
    Double I_ESP_Dxz_Gx3z = I_ESP_F2xz_F3z+ABX*I_ESP_Dxz_F3z;
    Double I_ESP_D2y_Gx3z = I_ESP_Fx2y_F3z+ABX*I_ESP_D2y_F3z;
    Double I_ESP_Dyz_Gx3z = I_ESP_Fxyz_F3z+ABX*I_ESP_Dyz_F3z;
    Double I_ESP_D2z_Gx3z = I_ESP_Fx2z_F3z+ABX*I_ESP_D2z_F3z;
    Double I_ESP_D2x_G4y = I_ESP_F2xy_F3y+ABY*I_ESP_D2x_F3y;
    Double I_ESP_Dxy_G4y = I_ESP_Fx2y_F3y+ABY*I_ESP_Dxy_F3y;
    Double I_ESP_Dxz_G4y = I_ESP_Fxyz_F3y+ABY*I_ESP_Dxz_F3y;
    Double I_ESP_D2y_G4y = I_ESP_F3y_F3y+ABY*I_ESP_D2y_F3y;
    Double I_ESP_Dyz_G4y = I_ESP_F2yz_F3y+ABY*I_ESP_Dyz_F3y;
    Double I_ESP_D2z_G4y = I_ESP_Fy2z_F3y+ABY*I_ESP_D2z_F3y;
    Double I_ESP_D2x_G3yz = I_ESP_F2xz_F3y+ABZ*I_ESP_D2x_F3y;
    Double I_ESP_Dxy_G3yz = I_ESP_Fxyz_F3y+ABZ*I_ESP_Dxy_F3y;
    Double I_ESP_Dxz_G3yz = I_ESP_Fx2z_F3y+ABZ*I_ESP_Dxz_F3y;
    Double I_ESP_D2y_G3yz = I_ESP_F2yz_F3y+ABZ*I_ESP_D2y_F3y;
    Double I_ESP_Dyz_G3yz = I_ESP_Fy2z_F3y+ABZ*I_ESP_Dyz_F3y;
    Double I_ESP_D2z_G3yz = I_ESP_F3z_F3y+ABZ*I_ESP_D2z_F3y;
    Double I_ESP_D2x_G2y2z = I_ESP_F2xz_F2yz+ABZ*I_ESP_D2x_F2yz;
    Double I_ESP_Dxy_G2y2z = I_ESP_Fxyz_F2yz+ABZ*I_ESP_Dxy_F2yz;
    Double I_ESP_Dxz_G2y2z = I_ESP_Fx2z_F2yz+ABZ*I_ESP_Dxz_F2yz;
    Double I_ESP_D2y_G2y2z = I_ESP_F2yz_F2yz+ABZ*I_ESP_D2y_F2yz;
    Double I_ESP_Dyz_G2y2z = I_ESP_Fy2z_F2yz+ABZ*I_ESP_Dyz_F2yz;
    Double I_ESP_D2z_G2y2z = I_ESP_F3z_F2yz+ABZ*I_ESP_D2z_F2yz;
    Double I_ESP_D2x_Gy3z = I_ESP_F2xy_F3z+ABY*I_ESP_D2x_F3z;
    Double I_ESP_Dxy_Gy3z = I_ESP_Fx2y_F3z+ABY*I_ESP_Dxy_F3z;
    Double I_ESP_Dxz_Gy3z = I_ESP_Fxyz_F3z+ABY*I_ESP_Dxz_F3z;
    Double I_ESP_D2y_Gy3z = I_ESP_F3y_F3z+ABY*I_ESP_D2y_F3z;
    Double I_ESP_Dyz_Gy3z = I_ESP_F2yz_F3z+ABY*I_ESP_Dyz_F3z;
    Double I_ESP_D2z_Gy3z = I_ESP_Fy2z_F3z+ABY*I_ESP_D2z_F3z;
    Double I_ESP_D2x_G4z = I_ESP_F2xz_F3z+ABZ*I_ESP_D2x_F3z;
    Double I_ESP_Dxy_G4z = I_ESP_Fxyz_F3z+ABZ*I_ESP_Dxy_F3z;
    Double I_ESP_Dxz_G4z = I_ESP_Fx2z_F3z+ABZ*I_ESP_Dxz_F3z;
    Double I_ESP_D2y_G4z = I_ESP_F2yz_F3z+ABZ*I_ESP_D2y_F3z;
    Double I_ESP_Dyz_G4z = I_ESP_Fy2z_F3z+ABZ*I_ESP_Dyz_F3z;
    Double I_ESP_D2z_G4z = I_ESP_F3z_F3z+ABZ*I_ESP_D2z_F3z;

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
     * totally 45 integrals are omitted 
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
     * totally 3 integrals are omitted 
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
     * shell quartet name: SQ_ESP_G_F_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 30 integrals are omitted 
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
     * shell quartet name: SQ_ESP_K_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 27 integrals are omitted 
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
    Double I_ESP_K5x2y_Py_a = I_ESP_L5x3y_S_a+ABY*I_ESP_K5x2y_S_a;
    Double I_ESP_K5xyz_Py_a = I_ESP_L5x2yz_S_a+ABY*I_ESP_K5xyz_S_a;
    Double I_ESP_K4x3y_Py_a = I_ESP_L4x4y_S_a+ABY*I_ESP_K4x3y_S_a;
    Double I_ESP_K4x2yz_Py_a = I_ESP_L4x3yz_S_a+ABY*I_ESP_K4x2yz_S_a;
    Double I_ESP_K4xy2z_Py_a = I_ESP_L4x2y2z_S_a+ABY*I_ESP_K4xy2z_S_a;
    Double I_ESP_K3x4y_Py_a = I_ESP_L3x5y_S_a+ABY*I_ESP_K3x4y_S_a;
    Double I_ESP_K3x3yz_Py_a = I_ESP_L3x4yz_S_a+ABY*I_ESP_K3x3yz_S_a;
    Double I_ESP_K3x2y2z_Py_a = I_ESP_L3x3y2z_S_a+ABY*I_ESP_K3x2y2z_S_a;
    Double I_ESP_K3xy3z_Py_a = I_ESP_L3x2y3z_S_a+ABY*I_ESP_K3xy3z_S_a;
    Double I_ESP_K2x5y_Py_a = I_ESP_L2x6y_S_a+ABY*I_ESP_K2x5y_S_a;
    Double I_ESP_K2x4yz_Py_a = I_ESP_L2x5yz_S_a+ABY*I_ESP_K2x4yz_S_a;
    Double I_ESP_K2x3y2z_Py_a = I_ESP_L2x4y2z_S_a+ABY*I_ESP_K2x3y2z_S_a;
    Double I_ESP_K2x2y3z_Py_a = I_ESP_L2x3y3z_S_a+ABY*I_ESP_K2x2y3z_S_a;
    Double I_ESP_K2xy4z_Py_a = I_ESP_L2x2y4z_S_a+ABY*I_ESP_K2xy4z_S_a;
    Double I_ESP_Kx6y_Py_a = I_ESP_Lx7y_S_a+ABY*I_ESP_Kx6y_S_a;
    Double I_ESP_Kx5yz_Py_a = I_ESP_Lx6yz_S_a+ABY*I_ESP_Kx5yz_S_a;
    Double I_ESP_Kx4y2z_Py_a = I_ESP_Lx5y2z_S_a+ABY*I_ESP_Kx4y2z_S_a;
    Double I_ESP_Kx3y3z_Py_a = I_ESP_Lx4y3z_S_a+ABY*I_ESP_Kx3y3z_S_a;
    Double I_ESP_Kx2y4z_Py_a = I_ESP_Lx3y4z_S_a+ABY*I_ESP_Kx2y4z_S_a;
    Double I_ESP_Kxy5z_Py_a = I_ESP_Lx2y5z_S_a+ABY*I_ESP_Kxy5z_S_a;
    Double I_ESP_K7y_Py_a = I_ESP_L8y_S_a+ABY*I_ESP_K7y_S_a;
    Double I_ESP_K6yz_Py_a = I_ESP_L7yz_S_a+ABY*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Py_a = I_ESP_L6y2z_S_a+ABY*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Py_a = I_ESP_L5y3z_S_a+ABY*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Py_a = I_ESP_L4y4z_S_a+ABY*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Py_a = I_ESP_L3y5z_S_a+ABY*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Py_a = I_ESP_L2y6z_S_a+ABY*I_ESP_Ky6z_S_a;
    Double I_ESP_K5xyz_Pz_a = I_ESP_L5xy2z_S_a+ABZ*I_ESP_K5xyz_S_a;
    Double I_ESP_K5x2z_Pz_a = I_ESP_L5x3z_S_a+ABZ*I_ESP_K5x2z_S_a;
    Double I_ESP_K4x2yz_Pz_a = I_ESP_L4x2y2z_S_a+ABZ*I_ESP_K4x2yz_S_a;
    Double I_ESP_K4xy2z_Pz_a = I_ESP_L4xy3z_S_a+ABZ*I_ESP_K4xy2z_S_a;
    Double I_ESP_K4x3z_Pz_a = I_ESP_L4x4z_S_a+ABZ*I_ESP_K4x3z_S_a;
    Double I_ESP_K3x3yz_Pz_a = I_ESP_L3x3y2z_S_a+ABZ*I_ESP_K3x3yz_S_a;
    Double I_ESP_K3x2y2z_Pz_a = I_ESP_L3x2y3z_S_a+ABZ*I_ESP_K3x2y2z_S_a;
    Double I_ESP_K3xy3z_Pz_a = I_ESP_L3xy4z_S_a+ABZ*I_ESP_K3xy3z_S_a;
    Double I_ESP_K3x4z_Pz_a = I_ESP_L3x5z_S_a+ABZ*I_ESP_K3x4z_S_a;
    Double I_ESP_K2x4yz_Pz_a = I_ESP_L2x4y2z_S_a+ABZ*I_ESP_K2x4yz_S_a;
    Double I_ESP_K2x3y2z_Pz_a = I_ESP_L2x3y3z_S_a+ABZ*I_ESP_K2x3y2z_S_a;
    Double I_ESP_K2x2y3z_Pz_a = I_ESP_L2x2y4z_S_a+ABZ*I_ESP_K2x2y3z_S_a;
    Double I_ESP_K2xy4z_Pz_a = I_ESP_L2xy5z_S_a+ABZ*I_ESP_K2xy4z_S_a;
    Double I_ESP_K2x5z_Pz_a = I_ESP_L2x6z_S_a+ABZ*I_ESP_K2x5z_S_a;
    Double I_ESP_Kx5yz_Pz_a = I_ESP_Lx5y2z_S_a+ABZ*I_ESP_Kx5yz_S_a;
    Double I_ESP_Kx4y2z_Pz_a = I_ESP_Lx4y3z_S_a+ABZ*I_ESP_Kx4y2z_S_a;
    Double I_ESP_Kx3y3z_Pz_a = I_ESP_Lx3y4z_S_a+ABZ*I_ESP_Kx3y3z_S_a;
    Double I_ESP_Kx2y4z_Pz_a = I_ESP_Lx2y5z_S_a+ABZ*I_ESP_Kx2y4z_S_a;
    Double I_ESP_Kxy5z_Pz_a = I_ESP_Lxy6z_S_a+ABZ*I_ESP_Kxy5z_S_a;
    Double I_ESP_Kx6z_Pz_a = I_ESP_Lx7z_S_a+ABZ*I_ESP_Kx6z_S_a;
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
     * totally 87 integrals are omitted 
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
     * totally 67 integrals are omitted 
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
    Double I_ESP_H4xz_F2xz_a = I_ESP_I4x2z_D2x_a+ABZ*I_ESP_H4xz_D2x_a;
    Double I_ESP_H3xyz_F2xz_a = I_ESP_I3xy2z_D2x_a+ABZ*I_ESP_H3xyz_D2x_a;
    Double I_ESP_H3x2z_F2xz_a = I_ESP_I3x3z_D2x_a+ABZ*I_ESP_H3x2z_D2x_a;
    Double I_ESP_H2x2yz_F2xz_a = I_ESP_I2x2y2z_D2x_a+ABZ*I_ESP_H2x2yz_D2x_a;
    Double I_ESP_H2xy2z_F2xz_a = I_ESP_I2xy3z_D2x_a+ABZ*I_ESP_H2xy2z_D2x_a;
    Double I_ESP_H2x3z_F2xz_a = I_ESP_I2x4z_D2x_a+ABZ*I_ESP_H2x3z_D2x_a;
    Double I_ESP_Hx3yz_F2xz_a = I_ESP_Ix3y2z_D2x_a+ABZ*I_ESP_Hx3yz_D2x_a;
    Double I_ESP_Hx2y2z_F2xz_a = I_ESP_Ix2y3z_D2x_a+ABZ*I_ESP_Hx2y2z_D2x_a;
    Double I_ESP_Hxy3z_F2xz_a = I_ESP_Ixy4z_D2x_a+ABZ*I_ESP_Hxy3z_D2x_a;
    Double I_ESP_Hx4z_F2xz_a = I_ESP_Ix5z_D2x_a+ABZ*I_ESP_Hx4z_D2x_a;
    Double I_ESP_H4yz_F2xz_a = I_ESP_I4y2z_D2x_a+ABZ*I_ESP_H4yz_D2x_a;
    Double I_ESP_H3y2z_F2xz_a = I_ESP_I3y3z_D2x_a+ABZ*I_ESP_H3y2z_D2x_a;
    Double I_ESP_H2y3z_F2xz_a = I_ESP_I2y4z_D2x_a+ABZ*I_ESP_H2y3z_D2x_a;
    Double I_ESP_Hy4z_F2xz_a = I_ESP_Iy5z_D2x_a+ABZ*I_ESP_Hy4z_D2x_a;
    Double I_ESP_H5z_F2xz_a = I_ESP_I6z_D2x_a+ABZ*I_ESP_H5z_D2x_a;
    Double I_ESP_H4xz_Fx2y_a = I_ESP_I5xz_D2y_a+ABX*I_ESP_H4xz_D2y_a;
    Double I_ESP_H3xyz_Fx2y_a = I_ESP_I4xyz_D2y_a+ABX*I_ESP_H3xyz_D2y_a;
    Double I_ESP_H3x2z_Fx2y_a = I_ESP_I4x2z_D2y_a+ABX*I_ESP_H3x2z_D2y_a;
    Double I_ESP_H2x2yz_Fx2y_a = I_ESP_I3x2yz_D2y_a+ABX*I_ESP_H2x2yz_D2y_a;
    Double I_ESP_H2xy2z_Fx2y_a = I_ESP_I3xy2z_D2y_a+ABX*I_ESP_H2xy2z_D2y_a;
    Double I_ESP_H2x3z_Fx2y_a = I_ESP_I3x3z_D2y_a+ABX*I_ESP_H2x3z_D2y_a;
    Double I_ESP_Hx3yz_Fx2y_a = I_ESP_I2x3yz_D2y_a+ABX*I_ESP_Hx3yz_D2y_a;
    Double I_ESP_Hx2y2z_Fx2y_a = I_ESP_I2x2y2z_D2y_a+ABX*I_ESP_Hx2y2z_D2y_a;
    Double I_ESP_Hxy3z_Fx2y_a = I_ESP_I2xy3z_D2y_a+ABX*I_ESP_Hxy3z_D2y_a;
    Double I_ESP_Hx4z_Fx2y_a = I_ESP_I2x4z_D2y_a+ABX*I_ESP_Hx4z_D2y_a;
    Double I_ESP_H4yz_Fx2y_a = I_ESP_Ix4yz_D2y_a+ABX*I_ESP_H4yz_D2y_a;
    Double I_ESP_H3y2z_Fx2y_a = I_ESP_Ix3y2z_D2y_a+ABX*I_ESP_H3y2z_D2y_a;
    Double I_ESP_H2y3z_Fx2y_a = I_ESP_Ix2y3z_D2y_a+ABX*I_ESP_H2y3z_D2y_a;
    Double I_ESP_Hy4z_Fx2y_a = I_ESP_Ixy4z_D2y_a+ABX*I_ESP_Hy4z_D2y_a;
    Double I_ESP_H5z_Fx2y_a = I_ESP_Ix5z_D2y_a+ABX*I_ESP_H5z_D2y_a;
    Double I_ESP_H4xy_Fx2z_a = I_ESP_I5xy_D2z_a+ABX*I_ESP_H4xy_D2z_a;
    Double I_ESP_H3x2y_Fx2z_a = I_ESP_I4x2y_D2z_a+ABX*I_ESP_H3x2y_D2z_a;
    Double I_ESP_H3xyz_Fx2z_a = I_ESP_I4xyz_D2z_a+ABX*I_ESP_H3xyz_D2z_a;
    Double I_ESP_H2x3y_Fx2z_a = I_ESP_I3x3y_D2z_a+ABX*I_ESP_H2x3y_D2z_a;
    Double I_ESP_H2x2yz_Fx2z_a = I_ESP_I3x2yz_D2z_a+ABX*I_ESP_H2x2yz_D2z_a;
    Double I_ESP_H2xy2z_Fx2z_a = I_ESP_I3xy2z_D2z_a+ABX*I_ESP_H2xy2z_D2z_a;
    Double I_ESP_Hx4y_Fx2z_a = I_ESP_I2x4y_D2z_a+ABX*I_ESP_Hx4y_D2z_a;
    Double I_ESP_Hx3yz_Fx2z_a = I_ESP_I2x3yz_D2z_a+ABX*I_ESP_Hx3yz_D2z_a;
    Double I_ESP_Hx2y2z_Fx2z_a = I_ESP_I2x2y2z_D2z_a+ABX*I_ESP_Hx2y2z_D2z_a;
    Double I_ESP_Hxy3z_Fx2z_a = I_ESP_I2xy3z_D2z_a+ABX*I_ESP_Hxy3z_D2z_a;
    Double I_ESP_H5y_Fx2z_a = I_ESP_Ix5y_D2z_a+ABX*I_ESP_H5y_D2z_a;
    Double I_ESP_H4yz_Fx2z_a = I_ESP_Ix4yz_D2z_a+ABX*I_ESP_H4yz_D2z_a;
    Double I_ESP_H3y2z_Fx2z_a = I_ESP_Ix3y2z_D2z_a+ABX*I_ESP_H3y2z_D2z_a;
    Double I_ESP_H2y3z_Fx2z_a = I_ESP_Ix2y3z_D2z_a+ABX*I_ESP_H2y3z_D2z_a;
    Double I_ESP_Hy4z_Fx2z_a = I_ESP_Ixy4z_D2z_a+ABX*I_ESP_Hy4z_D2z_a;
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
    Double I_ESP_H4xz_F2yz_a = I_ESP_I4x2z_D2y_a+ABZ*I_ESP_H4xz_D2y_a;
    Double I_ESP_H3xyz_F2yz_a = I_ESP_I3xy2z_D2y_a+ABZ*I_ESP_H3xyz_D2y_a;
    Double I_ESP_H3x2z_F2yz_a = I_ESP_I3x3z_D2y_a+ABZ*I_ESP_H3x2z_D2y_a;
    Double I_ESP_H2x2yz_F2yz_a = I_ESP_I2x2y2z_D2y_a+ABZ*I_ESP_H2x2yz_D2y_a;
    Double I_ESP_H2xy2z_F2yz_a = I_ESP_I2xy3z_D2y_a+ABZ*I_ESP_H2xy2z_D2y_a;
    Double I_ESP_H2x3z_F2yz_a = I_ESP_I2x4z_D2y_a+ABZ*I_ESP_H2x3z_D2y_a;
    Double I_ESP_Hx3yz_F2yz_a = I_ESP_Ix3y2z_D2y_a+ABZ*I_ESP_Hx3yz_D2y_a;
    Double I_ESP_Hx2y2z_F2yz_a = I_ESP_Ix2y3z_D2y_a+ABZ*I_ESP_Hx2y2z_D2y_a;
    Double I_ESP_Hxy3z_F2yz_a = I_ESP_Ixy4z_D2y_a+ABZ*I_ESP_Hxy3z_D2y_a;
    Double I_ESP_Hx4z_F2yz_a = I_ESP_Ix5z_D2y_a+ABZ*I_ESP_Hx4z_D2y_a;
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
     * shell quartet name: SQ_ESP_G_G_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_G_F_a
     ************************************************************/
    Double I_ESP_G4x_G4x_a = I_ESP_H5x_F3x_a+ABX*I_ESP_G4x_F3x_a;
    Double I_ESP_G3xy_G4x_a = I_ESP_H4xy_F3x_a+ABX*I_ESP_G3xy_F3x_a;
    Double I_ESP_G3xz_G4x_a = I_ESP_H4xz_F3x_a+ABX*I_ESP_G3xz_F3x_a;
    Double I_ESP_G2x2y_G4x_a = I_ESP_H3x2y_F3x_a+ABX*I_ESP_G2x2y_F3x_a;
    Double I_ESP_G2xyz_G4x_a = I_ESP_H3xyz_F3x_a+ABX*I_ESP_G2xyz_F3x_a;
    Double I_ESP_G2x2z_G4x_a = I_ESP_H3x2z_F3x_a+ABX*I_ESP_G2x2z_F3x_a;
    Double I_ESP_Gx3y_G4x_a = I_ESP_H2x3y_F3x_a+ABX*I_ESP_Gx3y_F3x_a;
    Double I_ESP_Gx2yz_G4x_a = I_ESP_H2x2yz_F3x_a+ABX*I_ESP_Gx2yz_F3x_a;
    Double I_ESP_Gxy2z_G4x_a = I_ESP_H2xy2z_F3x_a+ABX*I_ESP_Gxy2z_F3x_a;
    Double I_ESP_Gx3z_G4x_a = I_ESP_H2x3z_F3x_a+ABX*I_ESP_Gx3z_F3x_a;
    Double I_ESP_G4y_G4x_a = I_ESP_Hx4y_F3x_a+ABX*I_ESP_G4y_F3x_a;
    Double I_ESP_G3yz_G4x_a = I_ESP_Hx3yz_F3x_a+ABX*I_ESP_G3yz_F3x_a;
    Double I_ESP_G2y2z_G4x_a = I_ESP_Hx2y2z_F3x_a+ABX*I_ESP_G2y2z_F3x_a;
    Double I_ESP_Gy3z_G4x_a = I_ESP_Hxy3z_F3x_a+ABX*I_ESP_Gy3z_F3x_a;
    Double I_ESP_G4z_G4x_a = I_ESP_Hx4z_F3x_a+ABX*I_ESP_G4z_F3x_a;
    Double I_ESP_G4x_G3xy_a = I_ESP_H4xy_F3x_a+ABY*I_ESP_G4x_F3x_a;
    Double I_ESP_G3xy_G3xy_a = I_ESP_H3x2y_F3x_a+ABY*I_ESP_G3xy_F3x_a;
    Double I_ESP_G3xz_G3xy_a = I_ESP_H3xyz_F3x_a+ABY*I_ESP_G3xz_F3x_a;
    Double I_ESP_G2x2y_G3xy_a = I_ESP_H2x3y_F3x_a+ABY*I_ESP_G2x2y_F3x_a;
    Double I_ESP_G2xyz_G3xy_a = I_ESP_H2x2yz_F3x_a+ABY*I_ESP_G2xyz_F3x_a;
    Double I_ESP_G2x2z_G3xy_a = I_ESP_H2xy2z_F3x_a+ABY*I_ESP_G2x2z_F3x_a;
    Double I_ESP_Gx3y_G3xy_a = I_ESP_Hx4y_F3x_a+ABY*I_ESP_Gx3y_F3x_a;
    Double I_ESP_Gx2yz_G3xy_a = I_ESP_Hx3yz_F3x_a+ABY*I_ESP_Gx2yz_F3x_a;
    Double I_ESP_Gxy2z_G3xy_a = I_ESP_Hx2y2z_F3x_a+ABY*I_ESP_Gxy2z_F3x_a;
    Double I_ESP_Gx3z_G3xy_a = I_ESP_Hxy3z_F3x_a+ABY*I_ESP_Gx3z_F3x_a;
    Double I_ESP_G4y_G3xy_a = I_ESP_H5y_F3x_a+ABY*I_ESP_G4y_F3x_a;
    Double I_ESP_G3yz_G3xy_a = I_ESP_H4yz_F3x_a+ABY*I_ESP_G3yz_F3x_a;
    Double I_ESP_G2y2z_G3xy_a = I_ESP_H3y2z_F3x_a+ABY*I_ESP_G2y2z_F3x_a;
    Double I_ESP_Gy3z_G3xy_a = I_ESP_H2y3z_F3x_a+ABY*I_ESP_Gy3z_F3x_a;
    Double I_ESP_G4z_G3xy_a = I_ESP_Hy4z_F3x_a+ABY*I_ESP_G4z_F3x_a;
    Double I_ESP_G4x_G3xz_a = I_ESP_H4xz_F3x_a+ABZ*I_ESP_G4x_F3x_a;
    Double I_ESP_G3xy_G3xz_a = I_ESP_H3xyz_F3x_a+ABZ*I_ESP_G3xy_F3x_a;
    Double I_ESP_G3xz_G3xz_a = I_ESP_H3x2z_F3x_a+ABZ*I_ESP_G3xz_F3x_a;
    Double I_ESP_G2x2y_G3xz_a = I_ESP_H2x2yz_F3x_a+ABZ*I_ESP_G2x2y_F3x_a;
    Double I_ESP_G2xyz_G3xz_a = I_ESP_H2xy2z_F3x_a+ABZ*I_ESP_G2xyz_F3x_a;
    Double I_ESP_G2x2z_G3xz_a = I_ESP_H2x3z_F3x_a+ABZ*I_ESP_G2x2z_F3x_a;
    Double I_ESP_Gx3y_G3xz_a = I_ESP_Hx3yz_F3x_a+ABZ*I_ESP_Gx3y_F3x_a;
    Double I_ESP_Gx2yz_G3xz_a = I_ESP_Hx2y2z_F3x_a+ABZ*I_ESP_Gx2yz_F3x_a;
    Double I_ESP_Gxy2z_G3xz_a = I_ESP_Hxy3z_F3x_a+ABZ*I_ESP_Gxy2z_F3x_a;
    Double I_ESP_Gx3z_G3xz_a = I_ESP_Hx4z_F3x_a+ABZ*I_ESP_Gx3z_F3x_a;
    Double I_ESP_G4y_G3xz_a = I_ESP_H4yz_F3x_a+ABZ*I_ESP_G4y_F3x_a;
    Double I_ESP_G3yz_G3xz_a = I_ESP_H3y2z_F3x_a+ABZ*I_ESP_G3yz_F3x_a;
    Double I_ESP_G2y2z_G3xz_a = I_ESP_H2y3z_F3x_a+ABZ*I_ESP_G2y2z_F3x_a;
    Double I_ESP_Gy3z_G3xz_a = I_ESP_Hy4z_F3x_a+ABZ*I_ESP_Gy3z_F3x_a;
    Double I_ESP_G4z_G3xz_a = I_ESP_H5z_F3x_a+ABZ*I_ESP_G4z_F3x_a;
    Double I_ESP_G4x_G2x2y_a = I_ESP_H4xy_F2xy_a+ABY*I_ESP_G4x_F2xy_a;
    Double I_ESP_G3xy_G2x2y_a = I_ESP_H3x2y_F2xy_a+ABY*I_ESP_G3xy_F2xy_a;
    Double I_ESP_G3xz_G2x2y_a = I_ESP_H3xyz_F2xy_a+ABY*I_ESP_G3xz_F2xy_a;
    Double I_ESP_G2x2y_G2x2y_a = I_ESP_H2x3y_F2xy_a+ABY*I_ESP_G2x2y_F2xy_a;
    Double I_ESP_G2xyz_G2x2y_a = I_ESP_H2x2yz_F2xy_a+ABY*I_ESP_G2xyz_F2xy_a;
    Double I_ESP_G2x2z_G2x2y_a = I_ESP_H2xy2z_F2xy_a+ABY*I_ESP_G2x2z_F2xy_a;
    Double I_ESP_Gx3y_G2x2y_a = I_ESP_Hx4y_F2xy_a+ABY*I_ESP_Gx3y_F2xy_a;
    Double I_ESP_Gx2yz_G2x2y_a = I_ESP_Hx3yz_F2xy_a+ABY*I_ESP_Gx2yz_F2xy_a;
    Double I_ESP_Gxy2z_G2x2y_a = I_ESP_Hx2y2z_F2xy_a+ABY*I_ESP_Gxy2z_F2xy_a;
    Double I_ESP_Gx3z_G2x2y_a = I_ESP_Hxy3z_F2xy_a+ABY*I_ESP_Gx3z_F2xy_a;
    Double I_ESP_G4y_G2x2y_a = I_ESP_H5y_F2xy_a+ABY*I_ESP_G4y_F2xy_a;
    Double I_ESP_G3yz_G2x2y_a = I_ESP_H4yz_F2xy_a+ABY*I_ESP_G3yz_F2xy_a;
    Double I_ESP_G2y2z_G2x2y_a = I_ESP_H3y2z_F2xy_a+ABY*I_ESP_G2y2z_F2xy_a;
    Double I_ESP_Gy3z_G2x2y_a = I_ESP_H2y3z_F2xy_a+ABY*I_ESP_Gy3z_F2xy_a;
    Double I_ESP_G4z_G2x2y_a = I_ESP_Hy4z_F2xy_a+ABY*I_ESP_G4z_F2xy_a;
    Double I_ESP_G4x_G2xyz_a = I_ESP_H4xz_F2xy_a+ABZ*I_ESP_G4x_F2xy_a;
    Double I_ESP_G3xy_G2xyz_a = I_ESP_H3xyz_F2xy_a+ABZ*I_ESP_G3xy_F2xy_a;
    Double I_ESP_G3xz_G2xyz_a = I_ESP_H3x2z_F2xy_a+ABZ*I_ESP_G3xz_F2xy_a;
    Double I_ESP_G2x2y_G2xyz_a = I_ESP_H2x2yz_F2xy_a+ABZ*I_ESP_G2x2y_F2xy_a;
    Double I_ESP_G2xyz_G2xyz_a = I_ESP_H2xy2z_F2xy_a+ABZ*I_ESP_G2xyz_F2xy_a;
    Double I_ESP_G2x2z_G2xyz_a = I_ESP_H2x3z_F2xy_a+ABZ*I_ESP_G2x2z_F2xy_a;
    Double I_ESP_Gx3y_G2xyz_a = I_ESP_Hx3yz_F2xy_a+ABZ*I_ESP_Gx3y_F2xy_a;
    Double I_ESP_Gx2yz_G2xyz_a = I_ESP_Hx2y2z_F2xy_a+ABZ*I_ESP_Gx2yz_F2xy_a;
    Double I_ESP_Gxy2z_G2xyz_a = I_ESP_Hxy3z_F2xy_a+ABZ*I_ESP_Gxy2z_F2xy_a;
    Double I_ESP_Gx3z_G2xyz_a = I_ESP_Hx4z_F2xy_a+ABZ*I_ESP_Gx3z_F2xy_a;
    Double I_ESP_G4y_G2xyz_a = I_ESP_H4yz_F2xy_a+ABZ*I_ESP_G4y_F2xy_a;
    Double I_ESP_G3yz_G2xyz_a = I_ESP_H3y2z_F2xy_a+ABZ*I_ESP_G3yz_F2xy_a;
    Double I_ESP_G2y2z_G2xyz_a = I_ESP_H2y3z_F2xy_a+ABZ*I_ESP_G2y2z_F2xy_a;
    Double I_ESP_Gy3z_G2xyz_a = I_ESP_Hy4z_F2xy_a+ABZ*I_ESP_Gy3z_F2xy_a;
    Double I_ESP_G4z_G2xyz_a = I_ESP_H5z_F2xy_a+ABZ*I_ESP_G4z_F2xy_a;
    Double I_ESP_G4x_G2x2z_a = I_ESP_H4xz_F2xz_a+ABZ*I_ESP_G4x_F2xz_a;
    Double I_ESP_G3xy_G2x2z_a = I_ESP_H3xyz_F2xz_a+ABZ*I_ESP_G3xy_F2xz_a;
    Double I_ESP_G3xz_G2x2z_a = I_ESP_H3x2z_F2xz_a+ABZ*I_ESP_G3xz_F2xz_a;
    Double I_ESP_G2x2y_G2x2z_a = I_ESP_H2x2yz_F2xz_a+ABZ*I_ESP_G2x2y_F2xz_a;
    Double I_ESP_G2xyz_G2x2z_a = I_ESP_H2xy2z_F2xz_a+ABZ*I_ESP_G2xyz_F2xz_a;
    Double I_ESP_G2x2z_G2x2z_a = I_ESP_H2x3z_F2xz_a+ABZ*I_ESP_G2x2z_F2xz_a;
    Double I_ESP_Gx3y_G2x2z_a = I_ESP_Hx3yz_F2xz_a+ABZ*I_ESP_Gx3y_F2xz_a;
    Double I_ESP_Gx2yz_G2x2z_a = I_ESP_Hx2y2z_F2xz_a+ABZ*I_ESP_Gx2yz_F2xz_a;
    Double I_ESP_Gxy2z_G2x2z_a = I_ESP_Hxy3z_F2xz_a+ABZ*I_ESP_Gxy2z_F2xz_a;
    Double I_ESP_Gx3z_G2x2z_a = I_ESP_Hx4z_F2xz_a+ABZ*I_ESP_Gx3z_F2xz_a;
    Double I_ESP_G4y_G2x2z_a = I_ESP_H4yz_F2xz_a+ABZ*I_ESP_G4y_F2xz_a;
    Double I_ESP_G3yz_G2x2z_a = I_ESP_H3y2z_F2xz_a+ABZ*I_ESP_G3yz_F2xz_a;
    Double I_ESP_G2y2z_G2x2z_a = I_ESP_H2y3z_F2xz_a+ABZ*I_ESP_G2y2z_F2xz_a;
    Double I_ESP_Gy3z_G2x2z_a = I_ESP_Hy4z_F2xz_a+ABZ*I_ESP_Gy3z_F2xz_a;
    Double I_ESP_G4z_G2x2z_a = I_ESP_H5z_F2xz_a+ABZ*I_ESP_G4z_F2xz_a;
    Double I_ESP_G4x_Gx3y_a = I_ESP_H5x_F3y_a+ABX*I_ESP_G4x_F3y_a;
    Double I_ESP_G3xy_Gx3y_a = I_ESP_H4xy_F3y_a+ABX*I_ESP_G3xy_F3y_a;
    Double I_ESP_G3xz_Gx3y_a = I_ESP_H4xz_F3y_a+ABX*I_ESP_G3xz_F3y_a;
    Double I_ESP_G2x2y_Gx3y_a = I_ESP_H3x2y_F3y_a+ABX*I_ESP_G2x2y_F3y_a;
    Double I_ESP_G2xyz_Gx3y_a = I_ESP_H3xyz_F3y_a+ABX*I_ESP_G2xyz_F3y_a;
    Double I_ESP_G2x2z_Gx3y_a = I_ESP_H3x2z_F3y_a+ABX*I_ESP_G2x2z_F3y_a;
    Double I_ESP_Gx3y_Gx3y_a = I_ESP_H2x3y_F3y_a+ABX*I_ESP_Gx3y_F3y_a;
    Double I_ESP_Gx2yz_Gx3y_a = I_ESP_H2x2yz_F3y_a+ABX*I_ESP_Gx2yz_F3y_a;
    Double I_ESP_Gxy2z_Gx3y_a = I_ESP_H2xy2z_F3y_a+ABX*I_ESP_Gxy2z_F3y_a;
    Double I_ESP_Gx3z_Gx3y_a = I_ESP_H2x3z_F3y_a+ABX*I_ESP_Gx3z_F3y_a;
    Double I_ESP_G4y_Gx3y_a = I_ESP_Hx4y_F3y_a+ABX*I_ESP_G4y_F3y_a;
    Double I_ESP_G3yz_Gx3y_a = I_ESP_Hx3yz_F3y_a+ABX*I_ESP_G3yz_F3y_a;
    Double I_ESP_G2y2z_Gx3y_a = I_ESP_Hx2y2z_F3y_a+ABX*I_ESP_G2y2z_F3y_a;
    Double I_ESP_Gy3z_Gx3y_a = I_ESP_Hxy3z_F3y_a+ABX*I_ESP_Gy3z_F3y_a;
    Double I_ESP_G4z_Gx3y_a = I_ESP_Hx4z_F3y_a+ABX*I_ESP_G4z_F3y_a;
    Double I_ESP_G4x_Gx2yz_a = I_ESP_H4xz_Fx2y_a+ABZ*I_ESP_G4x_Fx2y_a;
    Double I_ESP_G3xy_Gx2yz_a = I_ESP_H3xyz_Fx2y_a+ABZ*I_ESP_G3xy_Fx2y_a;
    Double I_ESP_G3xz_Gx2yz_a = I_ESP_H3x2z_Fx2y_a+ABZ*I_ESP_G3xz_Fx2y_a;
    Double I_ESP_G2x2y_Gx2yz_a = I_ESP_H2x2yz_Fx2y_a+ABZ*I_ESP_G2x2y_Fx2y_a;
    Double I_ESP_G2xyz_Gx2yz_a = I_ESP_H2xy2z_Fx2y_a+ABZ*I_ESP_G2xyz_Fx2y_a;
    Double I_ESP_G2x2z_Gx2yz_a = I_ESP_H2x3z_Fx2y_a+ABZ*I_ESP_G2x2z_Fx2y_a;
    Double I_ESP_Gx3y_Gx2yz_a = I_ESP_Hx3yz_Fx2y_a+ABZ*I_ESP_Gx3y_Fx2y_a;
    Double I_ESP_Gx2yz_Gx2yz_a = I_ESP_Hx2y2z_Fx2y_a+ABZ*I_ESP_Gx2yz_Fx2y_a;
    Double I_ESP_Gxy2z_Gx2yz_a = I_ESP_Hxy3z_Fx2y_a+ABZ*I_ESP_Gxy2z_Fx2y_a;
    Double I_ESP_Gx3z_Gx2yz_a = I_ESP_Hx4z_Fx2y_a+ABZ*I_ESP_Gx3z_Fx2y_a;
    Double I_ESP_G4y_Gx2yz_a = I_ESP_H4yz_Fx2y_a+ABZ*I_ESP_G4y_Fx2y_a;
    Double I_ESP_G3yz_Gx2yz_a = I_ESP_H3y2z_Fx2y_a+ABZ*I_ESP_G3yz_Fx2y_a;
    Double I_ESP_G2y2z_Gx2yz_a = I_ESP_H2y3z_Fx2y_a+ABZ*I_ESP_G2y2z_Fx2y_a;
    Double I_ESP_Gy3z_Gx2yz_a = I_ESP_Hy4z_Fx2y_a+ABZ*I_ESP_Gy3z_Fx2y_a;
    Double I_ESP_G4z_Gx2yz_a = I_ESP_H5z_Fx2y_a+ABZ*I_ESP_G4z_Fx2y_a;
    Double I_ESP_G4x_Gxy2z_a = I_ESP_H4xy_Fx2z_a+ABY*I_ESP_G4x_Fx2z_a;
    Double I_ESP_G3xy_Gxy2z_a = I_ESP_H3x2y_Fx2z_a+ABY*I_ESP_G3xy_Fx2z_a;
    Double I_ESP_G3xz_Gxy2z_a = I_ESP_H3xyz_Fx2z_a+ABY*I_ESP_G3xz_Fx2z_a;
    Double I_ESP_G2x2y_Gxy2z_a = I_ESP_H2x3y_Fx2z_a+ABY*I_ESP_G2x2y_Fx2z_a;
    Double I_ESP_G2xyz_Gxy2z_a = I_ESP_H2x2yz_Fx2z_a+ABY*I_ESP_G2xyz_Fx2z_a;
    Double I_ESP_G2x2z_Gxy2z_a = I_ESP_H2xy2z_Fx2z_a+ABY*I_ESP_G2x2z_Fx2z_a;
    Double I_ESP_Gx3y_Gxy2z_a = I_ESP_Hx4y_Fx2z_a+ABY*I_ESP_Gx3y_Fx2z_a;
    Double I_ESP_Gx2yz_Gxy2z_a = I_ESP_Hx3yz_Fx2z_a+ABY*I_ESP_Gx2yz_Fx2z_a;
    Double I_ESP_Gxy2z_Gxy2z_a = I_ESP_Hx2y2z_Fx2z_a+ABY*I_ESP_Gxy2z_Fx2z_a;
    Double I_ESP_Gx3z_Gxy2z_a = I_ESP_Hxy3z_Fx2z_a+ABY*I_ESP_Gx3z_Fx2z_a;
    Double I_ESP_G4y_Gxy2z_a = I_ESP_H5y_Fx2z_a+ABY*I_ESP_G4y_Fx2z_a;
    Double I_ESP_G3yz_Gxy2z_a = I_ESP_H4yz_Fx2z_a+ABY*I_ESP_G3yz_Fx2z_a;
    Double I_ESP_G2y2z_Gxy2z_a = I_ESP_H3y2z_Fx2z_a+ABY*I_ESP_G2y2z_Fx2z_a;
    Double I_ESP_Gy3z_Gxy2z_a = I_ESP_H2y3z_Fx2z_a+ABY*I_ESP_Gy3z_Fx2z_a;
    Double I_ESP_G4z_Gxy2z_a = I_ESP_Hy4z_Fx2z_a+ABY*I_ESP_G4z_Fx2z_a;
    Double I_ESP_G4x_Gx3z_a = I_ESP_H5x_F3z_a+ABX*I_ESP_G4x_F3z_a;
    Double I_ESP_G3xy_Gx3z_a = I_ESP_H4xy_F3z_a+ABX*I_ESP_G3xy_F3z_a;
    Double I_ESP_G3xz_Gx3z_a = I_ESP_H4xz_F3z_a+ABX*I_ESP_G3xz_F3z_a;
    Double I_ESP_G2x2y_Gx3z_a = I_ESP_H3x2y_F3z_a+ABX*I_ESP_G2x2y_F3z_a;
    Double I_ESP_G2xyz_Gx3z_a = I_ESP_H3xyz_F3z_a+ABX*I_ESP_G2xyz_F3z_a;
    Double I_ESP_G2x2z_Gx3z_a = I_ESP_H3x2z_F3z_a+ABX*I_ESP_G2x2z_F3z_a;
    Double I_ESP_Gx3y_Gx3z_a = I_ESP_H2x3y_F3z_a+ABX*I_ESP_Gx3y_F3z_a;
    Double I_ESP_Gx2yz_Gx3z_a = I_ESP_H2x2yz_F3z_a+ABX*I_ESP_Gx2yz_F3z_a;
    Double I_ESP_Gxy2z_Gx3z_a = I_ESP_H2xy2z_F3z_a+ABX*I_ESP_Gxy2z_F3z_a;
    Double I_ESP_Gx3z_Gx3z_a = I_ESP_H2x3z_F3z_a+ABX*I_ESP_Gx3z_F3z_a;
    Double I_ESP_G4y_Gx3z_a = I_ESP_Hx4y_F3z_a+ABX*I_ESP_G4y_F3z_a;
    Double I_ESP_G3yz_Gx3z_a = I_ESP_Hx3yz_F3z_a+ABX*I_ESP_G3yz_F3z_a;
    Double I_ESP_G2y2z_Gx3z_a = I_ESP_Hx2y2z_F3z_a+ABX*I_ESP_G2y2z_F3z_a;
    Double I_ESP_Gy3z_Gx3z_a = I_ESP_Hxy3z_F3z_a+ABX*I_ESP_Gy3z_F3z_a;
    Double I_ESP_G4z_Gx3z_a = I_ESP_Hx4z_F3z_a+ABX*I_ESP_G4z_F3z_a;
    Double I_ESP_G4x_G4y_a = I_ESP_H4xy_F3y_a+ABY*I_ESP_G4x_F3y_a;
    Double I_ESP_G3xy_G4y_a = I_ESP_H3x2y_F3y_a+ABY*I_ESP_G3xy_F3y_a;
    Double I_ESP_G3xz_G4y_a = I_ESP_H3xyz_F3y_a+ABY*I_ESP_G3xz_F3y_a;
    Double I_ESP_G2x2y_G4y_a = I_ESP_H2x3y_F3y_a+ABY*I_ESP_G2x2y_F3y_a;
    Double I_ESP_G2xyz_G4y_a = I_ESP_H2x2yz_F3y_a+ABY*I_ESP_G2xyz_F3y_a;
    Double I_ESP_G2x2z_G4y_a = I_ESP_H2xy2z_F3y_a+ABY*I_ESP_G2x2z_F3y_a;
    Double I_ESP_Gx3y_G4y_a = I_ESP_Hx4y_F3y_a+ABY*I_ESP_Gx3y_F3y_a;
    Double I_ESP_Gx2yz_G4y_a = I_ESP_Hx3yz_F3y_a+ABY*I_ESP_Gx2yz_F3y_a;
    Double I_ESP_Gxy2z_G4y_a = I_ESP_Hx2y2z_F3y_a+ABY*I_ESP_Gxy2z_F3y_a;
    Double I_ESP_Gx3z_G4y_a = I_ESP_Hxy3z_F3y_a+ABY*I_ESP_Gx3z_F3y_a;
    Double I_ESP_G4y_G4y_a = I_ESP_H5y_F3y_a+ABY*I_ESP_G4y_F3y_a;
    Double I_ESP_G3yz_G4y_a = I_ESP_H4yz_F3y_a+ABY*I_ESP_G3yz_F3y_a;
    Double I_ESP_G2y2z_G4y_a = I_ESP_H3y2z_F3y_a+ABY*I_ESP_G2y2z_F3y_a;
    Double I_ESP_Gy3z_G4y_a = I_ESP_H2y3z_F3y_a+ABY*I_ESP_Gy3z_F3y_a;
    Double I_ESP_G4z_G4y_a = I_ESP_Hy4z_F3y_a+ABY*I_ESP_G4z_F3y_a;
    Double I_ESP_G4x_G3yz_a = I_ESP_H4xz_F3y_a+ABZ*I_ESP_G4x_F3y_a;
    Double I_ESP_G3xy_G3yz_a = I_ESP_H3xyz_F3y_a+ABZ*I_ESP_G3xy_F3y_a;
    Double I_ESP_G3xz_G3yz_a = I_ESP_H3x2z_F3y_a+ABZ*I_ESP_G3xz_F3y_a;
    Double I_ESP_G2x2y_G3yz_a = I_ESP_H2x2yz_F3y_a+ABZ*I_ESP_G2x2y_F3y_a;
    Double I_ESP_G2xyz_G3yz_a = I_ESP_H2xy2z_F3y_a+ABZ*I_ESP_G2xyz_F3y_a;
    Double I_ESP_G2x2z_G3yz_a = I_ESP_H2x3z_F3y_a+ABZ*I_ESP_G2x2z_F3y_a;
    Double I_ESP_Gx3y_G3yz_a = I_ESP_Hx3yz_F3y_a+ABZ*I_ESP_Gx3y_F3y_a;
    Double I_ESP_Gx2yz_G3yz_a = I_ESP_Hx2y2z_F3y_a+ABZ*I_ESP_Gx2yz_F3y_a;
    Double I_ESP_Gxy2z_G3yz_a = I_ESP_Hxy3z_F3y_a+ABZ*I_ESP_Gxy2z_F3y_a;
    Double I_ESP_Gx3z_G3yz_a = I_ESP_Hx4z_F3y_a+ABZ*I_ESP_Gx3z_F3y_a;
    Double I_ESP_G4y_G3yz_a = I_ESP_H4yz_F3y_a+ABZ*I_ESP_G4y_F3y_a;
    Double I_ESP_G3yz_G3yz_a = I_ESP_H3y2z_F3y_a+ABZ*I_ESP_G3yz_F3y_a;
    Double I_ESP_G2y2z_G3yz_a = I_ESP_H2y3z_F3y_a+ABZ*I_ESP_G2y2z_F3y_a;
    Double I_ESP_Gy3z_G3yz_a = I_ESP_Hy4z_F3y_a+ABZ*I_ESP_Gy3z_F3y_a;
    Double I_ESP_G4z_G3yz_a = I_ESP_H5z_F3y_a+ABZ*I_ESP_G4z_F3y_a;
    Double I_ESP_G4x_G2y2z_a = I_ESP_H4xz_F2yz_a+ABZ*I_ESP_G4x_F2yz_a;
    Double I_ESP_G3xy_G2y2z_a = I_ESP_H3xyz_F2yz_a+ABZ*I_ESP_G3xy_F2yz_a;
    Double I_ESP_G3xz_G2y2z_a = I_ESP_H3x2z_F2yz_a+ABZ*I_ESP_G3xz_F2yz_a;
    Double I_ESP_G2x2y_G2y2z_a = I_ESP_H2x2yz_F2yz_a+ABZ*I_ESP_G2x2y_F2yz_a;
    Double I_ESP_G2xyz_G2y2z_a = I_ESP_H2xy2z_F2yz_a+ABZ*I_ESP_G2xyz_F2yz_a;
    Double I_ESP_G2x2z_G2y2z_a = I_ESP_H2x3z_F2yz_a+ABZ*I_ESP_G2x2z_F2yz_a;
    Double I_ESP_Gx3y_G2y2z_a = I_ESP_Hx3yz_F2yz_a+ABZ*I_ESP_Gx3y_F2yz_a;
    Double I_ESP_Gx2yz_G2y2z_a = I_ESP_Hx2y2z_F2yz_a+ABZ*I_ESP_Gx2yz_F2yz_a;
    Double I_ESP_Gxy2z_G2y2z_a = I_ESP_Hxy3z_F2yz_a+ABZ*I_ESP_Gxy2z_F2yz_a;
    Double I_ESP_Gx3z_G2y2z_a = I_ESP_Hx4z_F2yz_a+ABZ*I_ESP_Gx3z_F2yz_a;
    Double I_ESP_G4y_G2y2z_a = I_ESP_H4yz_F2yz_a+ABZ*I_ESP_G4y_F2yz_a;
    Double I_ESP_G3yz_G2y2z_a = I_ESP_H3y2z_F2yz_a+ABZ*I_ESP_G3yz_F2yz_a;
    Double I_ESP_G2y2z_G2y2z_a = I_ESP_H2y3z_F2yz_a+ABZ*I_ESP_G2y2z_F2yz_a;
    Double I_ESP_Gy3z_G2y2z_a = I_ESP_Hy4z_F2yz_a+ABZ*I_ESP_Gy3z_F2yz_a;
    Double I_ESP_G4z_G2y2z_a = I_ESP_H5z_F2yz_a+ABZ*I_ESP_G4z_F2yz_a;
    Double I_ESP_G4x_Gy3z_a = I_ESP_H4xy_F3z_a+ABY*I_ESP_G4x_F3z_a;
    Double I_ESP_G3xy_Gy3z_a = I_ESP_H3x2y_F3z_a+ABY*I_ESP_G3xy_F3z_a;
    Double I_ESP_G3xz_Gy3z_a = I_ESP_H3xyz_F3z_a+ABY*I_ESP_G3xz_F3z_a;
    Double I_ESP_G2x2y_Gy3z_a = I_ESP_H2x3y_F3z_a+ABY*I_ESP_G2x2y_F3z_a;
    Double I_ESP_G2xyz_Gy3z_a = I_ESP_H2x2yz_F3z_a+ABY*I_ESP_G2xyz_F3z_a;
    Double I_ESP_G2x2z_Gy3z_a = I_ESP_H2xy2z_F3z_a+ABY*I_ESP_G2x2z_F3z_a;
    Double I_ESP_Gx3y_Gy3z_a = I_ESP_Hx4y_F3z_a+ABY*I_ESP_Gx3y_F3z_a;
    Double I_ESP_Gx2yz_Gy3z_a = I_ESP_Hx3yz_F3z_a+ABY*I_ESP_Gx2yz_F3z_a;
    Double I_ESP_Gxy2z_Gy3z_a = I_ESP_Hx2y2z_F3z_a+ABY*I_ESP_Gxy2z_F3z_a;
    Double I_ESP_Gx3z_Gy3z_a = I_ESP_Hxy3z_F3z_a+ABY*I_ESP_Gx3z_F3z_a;
    Double I_ESP_G4y_Gy3z_a = I_ESP_H5y_F3z_a+ABY*I_ESP_G4y_F3z_a;
    Double I_ESP_G3yz_Gy3z_a = I_ESP_H4yz_F3z_a+ABY*I_ESP_G3yz_F3z_a;
    Double I_ESP_G2y2z_Gy3z_a = I_ESP_H3y2z_F3z_a+ABY*I_ESP_G2y2z_F3z_a;
    Double I_ESP_Gy3z_Gy3z_a = I_ESP_H2y3z_F3z_a+ABY*I_ESP_Gy3z_F3z_a;
    Double I_ESP_G4z_Gy3z_a = I_ESP_Hy4z_F3z_a+ABY*I_ESP_G4z_F3z_a;
    Double I_ESP_G4x_G4z_a = I_ESP_H4xz_F3z_a+ABZ*I_ESP_G4x_F3z_a;
    Double I_ESP_G3xy_G4z_a = I_ESP_H3xyz_F3z_a+ABZ*I_ESP_G3xy_F3z_a;
    Double I_ESP_G3xz_G4z_a = I_ESP_H3x2z_F3z_a+ABZ*I_ESP_G3xz_F3z_a;
    Double I_ESP_G2x2y_G4z_a = I_ESP_H2x2yz_F3z_a+ABZ*I_ESP_G2x2y_F3z_a;
    Double I_ESP_G2xyz_G4z_a = I_ESP_H2xy2z_F3z_a+ABZ*I_ESP_G2xyz_F3z_a;
    Double I_ESP_G2x2z_G4z_a = I_ESP_H2x3z_F3z_a+ABZ*I_ESP_G2x2z_F3z_a;
    Double I_ESP_Gx3y_G4z_a = I_ESP_Hx3yz_F3z_a+ABZ*I_ESP_Gx3y_F3z_a;
    Double I_ESP_Gx2yz_G4z_a = I_ESP_Hx2y2z_F3z_a+ABZ*I_ESP_Gx2yz_F3z_a;
    Double I_ESP_Gxy2z_G4z_a = I_ESP_Hxy3z_F3z_a+ABZ*I_ESP_Gxy2z_F3z_a;
    Double I_ESP_Gx3z_G4z_a = I_ESP_Hx4z_F3z_a+ABZ*I_ESP_Gx3z_F3z_a;
    Double I_ESP_G4y_G4z_a = I_ESP_H4yz_F3z_a+ABZ*I_ESP_G4y_F3z_a;
    Double I_ESP_G3yz_G4z_a = I_ESP_H3y2z_F3z_a+ABZ*I_ESP_G3yz_F3z_a;
    Double I_ESP_G2y2z_G4z_a = I_ESP_H2y3z_F3z_a+ABZ*I_ESP_G2y2z_F3z_a;
    Double I_ESP_Gy3z_G4z_a = I_ESP_Hy4z_F3z_a+ABZ*I_ESP_Gy3z_F3z_a;
    Double I_ESP_G4z_G4z_a = I_ESP_H5z_F3z_a+ABZ*I_ESP_G4z_F3z_a;

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
     * totally 84 integrals are omitted 
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
     * totally 3 integrals are omitted 
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
    Double I_ESP_L7xy_Pz_aa = I_ESP_M7xyz_S_aa+ABZ*I_ESP_L7xy_S_aa;
    Double I_ESP_L7xz_Pz_aa = I_ESP_M7x2z_S_aa+ABZ*I_ESP_L7xz_S_aa;
    Double I_ESP_L6x2y_Pz_aa = I_ESP_M6x2yz_S_aa+ABZ*I_ESP_L6x2y_S_aa;
    Double I_ESP_L6xyz_Pz_aa = I_ESP_M6xy2z_S_aa+ABZ*I_ESP_L6xyz_S_aa;
    Double I_ESP_L6x2z_Pz_aa = I_ESP_M6x3z_S_aa+ABZ*I_ESP_L6x2z_S_aa;
    Double I_ESP_L5x3y_Pz_aa = I_ESP_M5x3yz_S_aa+ABZ*I_ESP_L5x3y_S_aa;
    Double I_ESP_L5x2yz_Pz_aa = I_ESP_M5x2y2z_S_aa+ABZ*I_ESP_L5x2yz_S_aa;
    Double I_ESP_L5xy2z_Pz_aa = I_ESP_M5xy3z_S_aa+ABZ*I_ESP_L5xy2z_S_aa;
    Double I_ESP_L5x3z_Pz_aa = I_ESP_M5x4z_S_aa+ABZ*I_ESP_L5x3z_S_aa;
    Double I_ESP_L4x4y_Pz_aa = I_ESP_M4x4yz_S_aa+ABZ*I_ESP_L4x4y_S_aa;
    Double I_ESP_L4x3yz_Pz_aa = I_ESP_M4x3y2z_S_aa+ABZ*I_ESP_L4x3yz_S_aa;
    Double I_ESP_L4x2y2z_Pz_aa = I_ESP_M4x2y3z_S_aa+ABZ*I_ESP_L4x2y2z_S_aa;
    Double I_ESP_L4xy3z_Pz_aa = I_ESP_M4xy4z_S_aa+ABZ*I_ESP_L4xy3z_S_aa;
    Double I_ESP_L4x4z_Pz_aa = I_ESP_M4x5z_S_aa+ABZ*I_ESP_L4x4z_S_aa;
    Double I_ESP_L3x5y_Pz_aa = I_ESP_M3x5yz_S_aa+ABZ*I_ESP_L3x5y_S_aa;
    Double I_ESP_L3x4yz_Pz_aa = I_ESP_M3x4y2z_S_aa+ABZ*I_ESP_L3x4yz_S_aa;
    Double I_ESP_L3x3y2z_Pz_aa = I_ESP_M3x3y3z_S_aa+ABZ*I_ESP_L3x3y2z_S_aa;
    Double I_ESP_L3x2y3z_Pz_aa = I_ESP_M3x2y4z_S_aa+ABZ*I_ESP_L3x2y3z_S_aa;
    Double I_ESP_L3xy4z_Pz_aa = I_ESP_M3xy5z_S_aa+ABZ*I_ESP_L3xy4z_S_aa;
    Double I_ESP_L3x5z_Pz_aa = I_ESP_M3x6z_S_aa+ABZ*I_ESP_L3x5z_S_aa;
    Double I_ESP_L2x6y_Pz_aa = I_ESP_M2x6yz_S_aa+ABZ*I_ESP_L2x6y_S_aa;
    Double I_ESP_L2x5yz_Pz_aa = I_ESP_M2x5y2z_S_aa+ABZ*I_ESP_L2x5yz_S_aa;
    Double I_ESP_L2x4y2z_Pz_aa = I_ESP_M2x4y3z_S_aa+ABZ*I_ESP_L2x4y2z_S_aa;
    Double I_ESP_L2x3y3z_Pz_aa = I_ESP_M2x3y4z_S_aa+ABZ*I_ESP_L2x3y3z_S_aa;
    Double I_ESP_L2x2y4z_Pz_aa = I_ESP_M2x2y5z_S_aa+ABZ*I_ESP_L2x2y4z_S_aa;
    Double I_ESP_L2xy5z_Pz_aa = I_ESP_M2xy6z_S_aa+ABZ*I_ESP_L2xy5z_S_aa;
    Double I_ESP_L2x6z_Pz_aa = I_ESP_M2x7z_S_aa+ABZ*I_ESP_L2x6z_S_aa;
    Double I_ESP_Lx7y_Pz_aa = I_ESP_Mx7yz_S_aa+ABZ*I_ESP_Lx7y_S_aa;
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
     * totally 108 integrals are omitted 
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
     * totally 56 integrals are omitted 
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
     * shell quartet name: SQ_ESP_M_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 33 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_N_S_aa
     * RHS shell quartet name: SQ_ESP_M_S_aa
     ************************************************************/
    Double I_ESP_M9x_Px_aa = I_ESP_N10x_S_aa+ABX*I_ESP_M9x_S_aa;
    Double I_ESP_M8xy_Px_aa = I_ESP_N9xy_S_aa+ABX*I_ESP_M8xy_S_aa;
    Double I_ESP_M8xz_Px_aa = I_ESP_N9xz_S_aa+ABX*I_ESP_M8xz_S_aa;
    Double I_ESP_M7x2y_Px_aa = I_ESP_N8x2y_S_aa+ABX*I_ESP_M7x2y_S_aa;
    Double I_ESP_M7xyz_Px_aa = I_ESP_N8xyz_S_aa+ABX*I_ESP_M7xyz_S_aa;
    Double I_ESP_M7x2z_Px_aa = I_ESP_N8x2z_S_aa+ABX*I_ESP_M7x2z_S_aa;
    Double I_ESP_M6x3y_Px_aa = I_ESP_N7x3y_S_aa+ABX*I_ESP_M6x3y_S_aa;
    Double I_ESP_M6x2yz_Px_aa = I_ESP_N7x2yz_S_aa+ABX*I_ESP_M6x2yz_S_aa;
    Double I_ESP_M6xy2z_Px_aa = I_ESP_N7xy2z_S_aa+ABX*I_ESP_M6xy2z_S_aa;
    Double I_ESP_M6x3z_Px_aa = I_ESP_N7x3z_S_aa+ABX*I_ESP_M6x3z_S_aa;
    Double I_ESP_M5x4y_Px_aa = I_ESP_N6x4y_S_aa+ABX*I_ESP_M5x4y_S_aa;
    Double I_ESP_M5x3yz_Px_aa = I_ESP_N6x3yz_S_aa+ABX*I_ESP_M5x3yz_S_aa;
    Double I_ESP_M5x2y2z_Px_aa = I_ESP_N6x2y2z_S_aa+ABX*I_ESP_M5x2y2z_S_aa;
    Double I_ESP_M5xy3z_Px_aa = I_ESP_N6xy3z_S_aa+ABX*I_ESP_M5xy3z_S_aa;
    Double I_ESP_M5x4z_Px_aa = I_ESP_N6x4z_S_aa+ABX*I_ESP_M5x4z_S_aa;
    Double I_ESP_M4x5y_Px_aa = I_ESP_N5x5y_S_aa+ABX*I_ESP_M4x5y_S_aa;
    Double I_ESP_M4x4yz_Px_aa = I_ESP_N5x4yz_S_aa+ABX*I_ESP_M4x4yz_S_aa;
    Double I_ESP_M4x3y2z_Px_aa = I_ESP_N5x3y2z_S_aa+ABX*I_ESP_M4x3y2z_S_aa;
    Double I_ESP_M4x2y3z_Px_aa = I_ESP_N5x2y3z_S_aa+ABX*I_ESP_M4x2y3z_S_aa;
    Double I_ESP_M4xy4z_Px_aa = I_ESP_N5xy4z_S_aa+ABX*I_ESP_M4xy4z_S_aa;
    Double I_ESP_M4x5z_Px_aa = I_ESP_N5x5z_S_aa+ABX*I_ESP_M4x5z_S_aa;
    Double I_ESP_M3x6y_Px_aa = I_ESP_N4x6y_S_aa+ABX*I_ESP_M3x6y_S_aa;
    Double I_ESP_M3x5yz_Px_aa = I_ESP_N4x5yz_S_aa+ABX*I_ESP_M3x5yz_S_aa;
    Double I_ESP_M3x4y2z_Px_aa = I_ESP_N4x4y2z_S_aa+ABX*I_ESP_M3x4y2z_S_aa;
    Double I_ESP_M3x3y3z_Px_aa = I_ESP_N4x3y3z_S_aa+ABX*I_ESP_M3x3y3z_S_aa;
    Double I_ESP_M3x2y4z_Px_aa = I_ESP_N4x2y4z_S_aa+ABX*I_ESP_M3x2y4z_S_aa;
    Double I_ESP_M3xy5z_Px_aa = I_ESP_N4xy5z_S_aa+ABX*I_ESP_M3xy5z_S_aa;
    Double I_ESP_M3x6z_Px_aa = I_ESP_N4x6z_S_aa+ABX*I_ESP_M3x6z_S_aa;
    Double I_ESP_M2x7y_Px_aa = I_ESP_N3x7y_S_aa+ABX*I_ESP_M2x7y_S_aa;
    Double I_ESP_M2x6yz_Px_aa = I_ESP_N3x6yz_S_aa+ABX*I_ESP_M2x6yz_S_aa;
    Double I_ESP_M2x5y2z_Px_aa = I_ESP_N3x5y2z_S_aa+ABX*I_ESP_M2x5y2z_S_aa;
    Double I_ESP_M2x4y3z_Px_aa = I_ESP_N3x4y3z_S_aa+ABX*I_ESP_M2x4y3z_S_aa;
    Double I_ESP_M2x3y4z_Px_aa = I_ESP_N3x3y4z_S_aa+ABX*I_ESP_M2x3y4z_S_aa;
    Double I_ESP_M2x2y5z_Px_aa = I_ESP_N3x2y5z_S_aa+ABX*I_ESP_M2x2y5z_S_aa;
    Double I_ESP_M2xy6z_Px_aa = I_ESP_N3xy6z_S_aa+ABX*I_ESP_M2xy6z_S_aa;
    Double I_ESP_M2x7z_Px_aa = I_ESP_N3x7z_S_aa+ABX*I_ESP_M2x7z_S_aa;
    Double I_ESP_Mx8y_Px_aa = I_ESP_N2x8y_S_aa+ABX*I_ESP_Mx8y_S_aa;
    Double I_ESP_Mx7yz_Px_aa = I_ESP_N2x7yz_S_aa+ABX*I_ESP_Mx7yz_S_aa;
    Double I_ESP_Mx6y2z_Px_aa = I_ESP_N2x6y2z_S_aa+ABX*I_ESP_Mx6y2z_S_aa;
    Double I_ESP_Mx5y3z_Px_aa = I_ESP_N2x5y3z_S_aa+ABX*I_ESP_Mx5y3z_S_aa;
    Double I_ESP_Mx4y4z_Px_aa = I_ESP_N2x4y4z_S_aa+ABX*I_ESP_Mx4y4z_S_aa;
    Double I_ESP_Mx3y5z_Px_aa = I_ESP_N2x3y5z_S_aa+ABX*I_ESP_Mx3y5z_S_aa;
    Double I_ESP_Mx2y6z_Px_aa = I_ESP_N2x2y6z_S_aa+ABX*I_ESP_Mx2y6z_S_aa;
    Double I_ESP_Mxy7z_Px_aa = I_ESP_N2xy7z_S_aa+ABX*I_ESP_Mxy7z_S_aa;
    Double I_ESP_Mx8z_Px_aa = I_ESP_N2x8z_S_aa+ABX*I_ESP_Mx8z_S_aa;
    Double I_ESP_M7x2y_Py_aa = I_ESP_N7x3y_S_aa+ABY*I_ESP_M7x2y_S_aa;
    Double I_ESP_M7xyz_Py_aa = I_ESP_N7x2yz_S_aa+ABY*I_ESP_M7xyz_S_aa;
    Double I_ESP_M6x3y_Py_aa = I_ESP_N6x4y_S_aa+ABY*I_ESP_M6x3y_S_aa;
    Double I_ESP_M6x2yz_Py_aa = I_ESP_N6x3yz_S_aa+ABY*I_ESP_M6x2yz_S_aa;
    Double I_ESP_M6xy2z_Py_aa = I_ESP_N6x2y2z_S_aa+ABY*I_ESP_M6xy2z_S_aa;
    Double I_ESP_M5x4y_Py_aa = I_ESP_N5x5y_S_aa+ABY*I_ESP_M5x4y_S_aa;
    Double I_ESP_M5x3yz_Py_aa = I_ESP_N5x4yz_S_aa+ABY*I_ESP_M5x3yz_S_aa;
    Double I_ESP_M5x2y2z_Py_aa = I_ESP_N5x3y2z_S_aa+ABY*I_ESP_M5x2y2z_S_aa;
    Double I_ESP_M5xy3z_Py_aa = I_ESP_N5x2y3z_S_aa+ABY*I_ESP_M5xy3z_S_aa;
    Double I_ESP_M4x5y_Py_aa = I_ESP_N4x6y_S_aa+ABY*I_ESP_M4x5y_S_aa;
    Double I_ESP_M4x4yz_Py_aa = I_ESP_N4x5yz_S_aa+ABY*I_ESP_M4x4yz_S_aa;
    Double I_ESP_M4x3y2z_Py_aa = I_ESP_N4x4y2z_S_aa+ABY*I_ESP_M4x3y2z_S_aa;
    Double I_ESP_M4x2y3z_Py_aa = I_ESP_N4x3y3z_S_aa+ABY*I_ESP_M4x2y3z_S_aa;
    Double I_ESP_M4xy4z_Py_aa = I_ESP_N4x2y4z_S_aa+ABY*I_ESP_M4xy4z_S_aa;
    Double I_ESP_M3x6y_Py_aa = I_ESP_N3x7y_S_aa+ABY*I_ESP_M3x6y_S_aa;
    Double I_ESP_M3x5yz_Py_aa = I_ESP_N3x6yz_S_aa+ABY*I_ESP_M3x5yz_S_aa;
    Double I_ESP_M3x4y2z_Py_aa = I_ESP_N3x5y2z_S_aa+ABY*I_ESP_M3x4y2z_S_aa;
    Double I_ESP_M3x3y3z_Py_aa = I_ESP_N3x4y3z_S_aa+ABY*I_ESP_M3x3y3z_S_aa;
    Double I_ESP_M3x2y4z_Py_aa = I_ESP_N3x3y4z_S_aa+ABY*I_ESP_M3x2y4z_S_aa;
    Double I_ESP_M3xy5z_Py_aa = I_ESP_N3x2y5z_S_aa+ABY*I_ESP_M3xy5z_S_aa;
    Double I_ESP_M2x7y_Py_aa = I_ESP_N2x8y_S_aa+ABY*I_ESP_M2x7y_S_aa;
    Double I_ESP_M2x6yz_Py_aa = I_ESP_N2x7yz_S_aa+ABY*I_ESP_M2x6yz_S_aa;
    Double I_ESP_M2x5y2z_Py_aa = I_ESP_N2x6y2z_S_aa+ABY*I_ESP_M2x5y2z_S_aa;
    Double I_ESP_M2x4y3z_Py_aa = I_ESP_N2x5y3z_S_aa+ABY*I_ESP_M2x4y3z_S_aa;
    Double I_ESP_M2x3y4z_Py_aa = I_ESP_N2x4y4z_S_aa+ABY*I_ESP_M2x3y4z_S_aa;
    Double I_ESP_M2x2y5z_Py_aa = I_ESP_N2x3y5z_S_aa+ABY*I_ESP_M2x2y5z_S_aa;
    Double I_ESP_M2xy6z_Py_aa = I_ESP_N2x2y6z_S_aa+ABY*I_ESP_M2xy6z_S_aa;
    Double I_ESP_Mx8y_Py_aa = I_ESP_Nx9y_S_aa+ABY*I_ESP_Mx8y_S_aa;
    Double I_ESP_Mx7yz_Py_aa = I_ESP_Nx8yz_S_aa+ABY*I_ESP_Mx7yz_S_aa;
    Double I_ESP_Mx6y2z_Py_aa = I_ESP_Nx7y2z_S_aa+ABY*I_ESP_Mx6y2z_S_aa;
    Double I_ESP_Mx5y3z_Py_aa = I_ESP_Nx6y3z_S_aa+ABY*I_ESP_Mx5y3z_S_aa;
    Double I_ESP_Mx4y4z_Py_aa = I_ESP_Nx5y4z_S_aa+ABY*I_ESP_Mx4y4z_S_aa;
    Double I_ESP_Mx3y5z_Py_aa = I_ESP_Nx4y5z_S_aa+ABY*I_ESP_Mx3y5z_S_aa;
    Double I_ESP_Mx2y6z_Py_aa = I_ESP_Nx3y6z_S_aa+ABY*I_ESP_Mx2y6z_S_aa;
    Double I_ESP_Mxy7z_Py_aa = I_ESP_Nx2y7z_S_aa+ABY*I_ESP_Mxy7z_S_aa;
    Double I_ESP_M9y_Py_aa = I_ESP_N10y_S_aa+ABY*I_ESP_M9y_S_aa;
    Double I_ESP_M8yz_Py_aa = I_ESP_N9yz_S_aa+ABY*I_ESP_M8yz_S_aa;
    Double I_ESP_M7y2z_Py_aa = I_ESP_N8y2z_S_aa+ABY*I_ESP_M7y2z_S_aa;
    Double I_ESP_M6y3z_Py_aa = I_ESP_N7y3z_S_aa+ABY*I_ESP_M6y3z_S_aa;
    Double I_ESP_M5y4z_Py_aa = I_ESP_N6y4z_S_aa+ABY*I_ESP_M5y4z_S_aa;
    Double I_ESP_M4y5z_Py_aa = I_ESP_N5y5z_S_aa+ABY*I_ESP_M4y5z_S_aa;
    Double I_ESP_M3y6z_Py_aa = I_ESP_N4y6z_S_aa+ABY*I_ESP_M3y6z_S_aa;
    Double I_ESP_M2y7z_Py_aa = I_ESP_N3y7z_S_aa+ABY*I_ESP_M2y7z_S_aa;
    Double I_ESP_My8z_Py_aa = I_ESP_N2y8z_S_aa+ABY*I_ESP_My8z_S_aa;
    Double I_ESP_M7xyz_Pz_aa = I_ESP_N7xy2z_S_aa+ABZ*I_ESP_M7xyz_S_aa;
    Double I_ESP_M7x2z_Pz_aa = I_ESP_N7x3z_S_aa+ABZ*I_ESP_M7x2z_S_aa;
    Double I_ESP_M6x2yz_Pz_aa = I_ESP_N6x2y2z_S_aa+ABZ*I_ESP_M6x2yz_S_aa;
    Double I_ESP_M6xy2z_Pz_aa = I_ESP_N6xy3z_S_aa+ABZ*I_ESP_M6xy2z_S_aa;
    Double I_ESP_M6x3z_Pz_aa = I_ESP_N6x4z_S_aa+ABZ*I_ESP_M6x3z_S_aa;
    Double I_ESP_M5x3yz_Pz_aa = I_ESP_N5x3y2z_S_aa+ABZ*I_ESP_M5x3yz_S_aa;
    Double I_ESP_M5x2y2z_Pz_aa = I_ESP_N5x2y3z_S_aa+ABZ*I_ESP_M5x2y2z_S_aa;
    Double I_ESP_M5xy3z_Pz_aa = I_ESP_N5xy4z_S_aa+ABZ*I_ESP_M5xy3z_S_aa;
    Double I_ESP_M5x4z_Pz_aa = I_ESP_N5x5z_S_aa+ABZ*I_ESP_M5x4z_S_aa;
    Double I_ESP_M4x4yz_Pz_aa = I_ESP_N4x4y2z_S_aa+ABZ*I_ESP_M4x4yz_S_aa;
    Double I_ESP_M4x3y2z_Pz_aa = I_ESP_N4x3y3z_S_aa+ABZ*I_ESP_M4x3y2z_S_aa;
    Double I_ESP_M4x2y3z_Pz_aa = I_ESP_N4x2y4z_S_aa+ABZ*I_ESP_M4x2y3z_S_aa;
    Double I_ESP_M4xy4z_Pz_aa = I_ESP_N4xy5z_S_aa+ABZ*I_ESP_M4xy4z_S_aa;
    Double I_ESP_M4x5z_Pz_aa = I_ESP_N4x6z_S_aa+ABZ*I_ESP_M4x5z_S_aa;
    Double I_ESP_M3x5yz_Pz_aa = I_ESP_N3x5y2z_S_aa+ABZ*I_ESP_M3x5yz_S_aa;
    Double I_ESP_M3x4y2z_Pz_aa = I_ESP_N3x4y3z_S_aa+ABZ*I_ESP_M3x4y2z_S_aa;
    Double I_ESP_M3x3y3z_Pz_aa = I_ESP_N3x3y4z_S_aa+ABZ*I_ESP_M3x3y3z_S_aa;
    Double I_ESP_M3x2y4z_Pz_aa = I_ESP_N3x2y5z_S_aa+ABZ*I_ESP_M3x2y4z_S_aa;
    Double I_ESP_M3xy5z_Pz_aa = I_ESP_N3xy6z_S_aa+ABZ*I_ESP_M3xy5z_S_aa;
    Double I_ESP_M3x6z_Pz_aa = I_ESP_N3x7z_S_aa+ABZ*I_ESP_M3x6z_S_aa;
    Double I_ESP_M2x6yz_Pz_aa = I_ESP_N2x6y2z_S_aa+ABZ*I_ESP_M2x6yz_S_aa;
    Double I_ESP_M2x5y2z_Pz_aa = I_ESP_N2x5y3z_S_aa+ABZ*I_ESP_M2x5y2z_S_aa;
    Double I_ESP_M2x4y3z_Pz_aa = I_ESP_N2x4y4z_S_aa+ABZ*I_ESP_M2x4y3z_S_aa;
    Double I_ESP_M2x3y4z_Pz_aa = I_ESP_N2x3y5z_S_aa+ABZ*I_ESP_M2x3y4z_S_aa;
    Double I_ESP_M2x2y5z_Pz_aa = I_ESP_N2x2y6z_S_aa+ABZ*I_ESP_M2x2y5z_S_aa;
    Double I_ESP_M2xy6z_Pz_aa = I_ESP_N2xy7z_S_aa+ABZ*I_ESP_M2xy6z_S_aa;
    Double I_ESP_M2x7z_Pz_aa = I_ESP_N2x8z_S_aa+ABZ*I_ESP_M2x7z_S_aa;
    Double I_ESP_Mx7yz_Pz_aa = I_ESP_Nx7y2z_S_aa+ABZ*I_ESP_Mx7yz_S_aa;
    Double I_ESP_Mx6y2z_Pz_aa = I_ESP_Nx6y3z_S_aa+ABZ*I_ESP_Mx6y2z_S_aa;
    Double I_ESP_Mx5y3z_Pz_aa = I_ESP_Nx5y4z_S_aa+ABZ*I_ESP_Mx5y3z_S_aa;
    Double I_ESP_Mx4y4z_Pz_aa = I_ESP_Nx4y5z_S_aa+ABZ*I_ESP_Mx4y4z_S_aa;
    Double I_ESP_Mx3y5z_Pz_aa = I_ESP_Nx3y6z_S_aa+ABZ*I_ESP_Mx3y5z_S_aa;
    Double I_ESP_Mx2y6z_Pz_aa = I_ESP_Nx2y7z_S_aa+ABZ*I_ESP_Mx2y6z_S_aa;
    Double I_ESP_Mxy7z_Pz_aa = I_ESP_Nxy8z_S_aa+ABZ*I_ESP_Mxy7z_S_aa;
    Double I_ESP_Mx8z_Pz_aa = I_ESP_Nx9z_S_aa+ABZ*I_ESP_Mx8z_S_aa;
    Double I_ESP_M7y2z_Pz_aa = I_ESP_N7y3z_S_aa+ABZ*I_ESP_M7y2z_S_aa;
    Double I_ESP_M6y3z_Pz_aa = I_ESP_N6y4z_S_aa+ABZ*I_ESP_M6y3z_S_aa;
    Double I_ESP_M5y4z_Pz_aa = I_ESP_N5y5z_S_aa+ABZ*I_ESP_M5y4z_S_aa;
    Double I_ESP_M4y5z_Pz_aa = I_ESP_N4y6z_S_aa+ABZ*I_ESP_M4y5z_S_aa;
    Double I_ESP_M3y6z_Pz_aa = I_ESP_N3y7z_S_aa+ABZ*I_ESP_M3y6z_S_aa;
    Double I_ESP_M2y7z_Pz_aa = I_ESP_N2y8z_S_aa+ABZ*I_ESP_M2y7z_S_aa;
    Double I_ESP_My8z_Pz_aa = I_ESP_Ny9z_S_aa+ABZ*I_ESP_My8z_S_aa;
    Double I_ESP_M9z_Pz_aa = I_ESP_N10z_S_aa+ABZ*I_ESP_M9z_S_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_L_D_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 138 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_M_P_aa
     * RHS shell quartet name: SQ_ESP_L_P_aa
     ************************************************************/
    Double I_ESP_L8x_D2x_aa = I_ESP_M9x_Px_aa+ABX*I_ESP_L8x_Px_aa;
    Double I_ESP_L7xy_D2x_aa = I_ESP_M8xy_Px_aa+ABX*I_ESP_L7xy_Px_aa;
    Double I_ESP_L7xz_D2x_aa = I_ESP_M8xz_Px_aa+ABX*I_ESP_L7xz_Px_aa;
    Double I_ESP_L6x2y_D2x_aa = I_ESP_M7x2y_Px_aa+ABX*I_ESP_L6x2y_Px_aa;
    Double I_ESP_L6xyz_D2x_aa = I_ESP_M7xyz_Px_aa+ABX*I_ESP_L6xyz_Px_aa;
    Double I_ESP_L6x2z_D2x_aa = I_ESP_M7x2z_Px_aa+ABX*I_ESP_L6x2z_Px_aa;
    Double I_ESP_L5x3y_D2x_aa = I_ESP_M6x3y_Px_aa+ABX*I_ESP_L5x3y_Px_aa;
    Double I_ESP_L5x2yz_D2x_aa = I_ESP_M6x2yz_Px_aa+ABX*I_ESP_L5x2yz_Px_aa;
    Double I_ESP_L5xy2z_D2x_aa = I_ESP_M6xy2z_Px_aa+ABX*I_ESP_L5xy2z_Px_aa;
    Double I_ESP_L5x3z_D2x_aa = I_ESP_M6x3z_Px_aa+ABX*I_ESP_L5x3z_Px_aa;
    Double I_ESP_L4x4y_D2x_aa = I_ESP_M5x4y_Px_aa+ABX*I_ESP_L4x4y_Px_aa;
    Double I_ESP_L4x3yz_D2x_aa = I_ESP_M5x3yz_Px_aa+ABX*I_ESP_L4x3yz_Px_aa;
    Double I_ESP_L4x2y2z_D2x_aa = I_ESP_M5x2y2z_Px_aa+ABX*I_ESP_L4x2y2z_Px_aa;
    Double I_ESP_L4xy3z_D2x_aa = I_ESP_M5xy3z_Px_aa+ABX*I_ESP_L4xy3z_Px_aa;
    Double I_ESP_L4x4z_D2x_aa = I_ESP_M5x4z_Px_aa+ABX*I_ESP_L4x4z_Px_aa;
    Double I_ESP_L3x5y_D2x_aa = I_ESP_M4x5y_Px_aa+ABX*I_ESP_L3x5y_Px_aa;
    Double I_ESP_L3x4yz_D2x_aa = I_ESP_M4x4yz_Px_aa+ABX*I_ESP_L3x4yz_Px_aa;
    Double I_ESP_L3x3y2z_D2x_aa = I_ESP_M4x3y2z_Px_aa+ABX*I_ESP_L3x3y2z_Px_aa;
    Double I_ESP_L3x2y3z_D2x_aa = I_ESP_M4x2y3z_Px_aa+ABX*I_ESP_L3x2y3z_Px_aa;
    Double I_ESP_L3xy4z_D2x_aa = I_ESP_M4xy4z_Px_aa+ABX*I_ESP_L3xy4z_Px_aa;
    Double I_ESP_L3x5z_D2x_aa = I_ESP_M4x5z_Px_aa+ABX*I_ESP_L3x5z_Px_aa;
    Double I_ESP_L2x6y_D2x_aa = I_ESP_M3x6y_Px_aa+ABX*I_ESP_L2x6y_Px_aa;
    Double I_ESP_L2x5yz_D2x_aa = I_ESP_M3x5yz_Px_aa+ABX*I_ESP_L2x5yz_Px_aa;
    Double I_ESP_L2x4y2z_D2x_aa = I_ESP_M3x4y2z_Px_aa+ABX*I_ESP_L2x4y2z_Px_aa;
    Double I_ESP_L2x3y3z_D2x_aa = I_ESP_M3x3y3z_Px_aa+ABX*I_ESP_L2x3y3z_Px_aa;
    Double I_ESP_L2x2y4z_D2x_aa = I_ESP_M3x2y4z_Px_aa+ABX*I_ESP_L2x2y4z_Px_aa;
    Double I_ESP_L2xy5z_D2x_aa = I_ESP_M3xy5z_Px_aa+ABX*I_ESP_L2xy5z_Px_aa;
    Double I_ESP_L2x6z_D2x_aa = I_ESP_M3x6z_Px_aa+ABX*I_ESP_L2x6z_Px_aa;
    Double I_ESP_Lx7y_D2x_aa = I_ESP_M2x7y_Px_aa+ABX*I_ESP_Lx7y_Px_aa;
    Double I_ESP_Lx6yz_D2x_aa = I_ESP_M2x6yz_Px_aa+ABX*I_ESP_Lx6yz_Px_aa;
    Double I_ESP_Lx5y2z_D2x_aa = I_ESP_M2x5y2z_Px_aa+ABX*I_ESP_Lx5y2z_Px_aa;
    Double I_ESP_Lx4y3z_D2x_aa = I_ESP_M2x4y3z_Px_aa+ABX*I_ESP_Lx4y3z_Px_aa;
    Double I_ESP_Lx3y4z_D2x_aa = I_ESP_M2x3y4z_Px_aa+ABX*I_ESP_Lx3y4z_Px_aa;
    Double I_ESP_Lx2y5z_D2x_aa = I_ESP_M2x2y5z_Px_aa+ABX*I_ESP_Lx2y5z_Px_aa;
    Double I_ESP_Lxy6z_D2x_aa = I_ESP_M2xy6z_Px_aa+ABX*I_ESP_Lxy6z_Px_aa;
    Double I_ESP_Lx7z_D2x_aa = I_ESP_M2x7z_Px_aa+ABX*I_ESP_Lx7z_Px_aa;
    Double I_ESP_L8y_D2x_aa = I_ESP_Mx8y_Px_aa+ABX*I_ESP_L8y_Px_aa;
    Double I_ESP_L7yz_D2x_aa = I_ESP_Mx7yz_Px_aa+ABX*I_ESP_L7yz_Px_aa;
    Double I_ESP_L6y2z_D2x_aa = I_ESP_Mx6y2z_Px_aa+ABX*I_ESP_L6y2z_Px_aa;
    Double I_ESP_L5y3z_D2x_aa = I_ESP_Mx5y3z_Px_aa+ABX*I_ESP_L5y3z_Px_aa;
    Double I_ESP_L4y4z_D2x_aa = I_ESP_Mx4y4z_Px_aa+ABX*I_ESP_L4y4z_Px_aa;
    Double I_ESP_L3y5z_D2x_aa = I_ESP_Mx3y5z_Px_aa+ABX*I_ESP_L3y5z_Px_aa;
    Double I_ESP_L2y6z_D2x_aa = I_ESP_Mx2y6z_Px_aa+ABX*I_ESP_L2y6z_Px_aa;
    Double I_ESP_Ly7z_D2x_aa = I_ESP_Mxy7z_Px_aa+ABX*I_ESP_Ly7z_Px_aa;
    Double I_ESP_L8z_D2x_aa = I_ESP_Mx8z_Px_aa+ABX*I_ESP_L8z_Px_aa;
    Double I_ESP_L7xy_D2y_aa = I_ESP_M7x2y_Py_aa+ABY*I_ESP_L7xy_Py_aa;
    Double I_ESP_L7xz_D2y_aa = I_ESP_M7xyz_Py_aa+ABY*I_ESP_L7xz_Py_aa;
    Double I_ESP_L6x2y_D2y_aa = I_ESP_M6x3y_Py_aa+ABY*I_ESP_L6x2y_Py_aa;
    Double I_ESP_L6xyz_D2y_aa = I_ESP_M6x2yz_Py_aa+ABY*I_ESP_L6xyz_Py_aa;
    Double I_ESP_L6x2z_D2y_aa = I_ESP_M6xy2z_Py_aa+ABY*I_ESP_L6x2z_Py_aa;
    Double I_ESP_L5x3y_D2y_aa = I_ESP_M5x4y_Py_aa+ABY*I_ESP_L5x3y_Py_aa;
    Double I_ESP_L5x2yz_D2y_aa = I_ESP_M5x3yz_Py_aa+ABY*I_ESP_L5x2yz_Py_aa;
    Double I_ESP_L5xy2z_D2y_aa = I_ESP_M5x2y2z_Py_aa+ABY*I_ESP_L5xy2z_Py_aa;
    Double I_ESP_L5x3z_D2y_aa = I_ESP_M5xy3z_Py_aa+ABY*I_ESP_L5x3z_Py_aa;
    Double I_ESP_L4x4y_D2y_aa = I_ESP_M4x5y_Py_aa+ABY*I_ESP_L4x4y_Py_aa;
    Double I_ESP_L4x3yz_D2y_aa = I_ESP_M4x4yz_Py_aa+ABY*I_ESP_L4x3yz_Py_aa;
    Double I_ESP_L4x2y2z_D2y_aa = I_ESP_M4x3y2z_Py_aa+ABY*I_ESP_L4x2y2z_Py_aa;
    Double I_ESP_L4xy3z_D2y_aa = I_ESP_M4x2y3z_Py_aa+ABY*I_ESP_L4xy3z_Py_aa;
    Double I_ESP_L4x4z_D2y_aa = I_ESP_M4xy4z_Py_aa+ABY*I_ESP_L4x4z_Py_aa;
    Double I_ESP_L3x5y_D2y_aa = I_ESP_M3x6y_Py_aa+ABY*I_ESP_L3x5y_Py_aa;
    Double I_ESP_L3x4yz_D2y_aa = I_ESP_M3x5yz_Py_aa+ABY*I_ESP_L3x4yz_Py_aa;
    Double I_ESP_L3x3y2z_D2y_aa = I_ESP_M3x4y2z_Py_aa+ABY*I_ESP_L3x3y2z_Py_aa;
    Double I_ESP_L3x2y3z_D2y_aa = I_ESP_M3x3y3z_Py_aa+ABY*I_ESP_L3x2y3z_Py_aa;
    Double I_ESP_L3xy4z_D2y_aa = I_ESP_M3x2y4z_Py_aa+ABY*I_ESP_L3xy4z_Py_aa;
    Double I_ESP_L3x5z_D2y_aa = I_ESP_M3xy5z_Py_aa+ABY*I_ESP_L3x5z_Py_aa;
    Double I_ESP_L2x6y_D2y_aa = I_ESP_M2x7y_Py_aa+ABY*I_ESP_L2x6y_Py_aa;
    Double I_ESP_L2x5yz_D2y_aa = I_ESP_M2x6yz_Py_aa+ABY*I_ESP_L2x5yz_Py_aa;
    Double I_ESP_L2x4y2z_D2y_aa = I_ESP_M2x5y2z_Py_aa+ABY*I_ESP_L2x4y2z_Py_aa;
    Double I_ESP_L2x3y3z_D2y_aa = I_ESP_M2x4y3z_Py_aa+ABY*I_ESP_L2x3y3z_Py_aa;
    Double I_ESP_L2x2y4z_D2y_aa = I_ESP_M2x3y4z_Py_aa+ABY*I_ESP_L2x2y4z_Py_aa;
    Double I_ESP_L2xy5z_D2y_aa = I_ESP_M2x2y5z_Py_aa+ABY*I_ESP_L2xy5z_Py_aa;
    Double I_ESP_L2x6z_D2y_aa = I_ESP_M2xy6z_Py_aa+ABY*I_ESP_L2x6z_Py_aa;
    Double I_ESP_Lx7y_D2y_aa = I_ESP_Mx8y_Py_aa+ABY*I_ESP_Lx7y_Py_aa;
    Double I_ESP_Lx6yz_D2y_aa = I_ESP_Mx7yz_Py_aa+ABY*I_ESP_Lx6yz_Py_aa;
    Double I_ESP_Lx5y2z_D2y_aa = I_ESP_Mx6y2z_Py_aa+ABY*I_ESP_Lx5y2z_Py_aa;
    Double I_ESP_Lx4y3z_D2y_aa = I_ESP_Mx5y3z_Py_aa+ABY*I_ESP_Lx4y3z_Py_aa;
    Double I_ESP_Lx3y4z_D2y_aa = I_ESP_Mx4y4z_Py_aa+ABY*I_ESP_Lx3y4z_Py_aa;
    Double I_ESP_Lx2y5z_D2y_aa = I_ESP_Mx3y5z_Py_aa+ABY*I_ESP_Lx2y5z_Py_aa;
    Double I_ESP_Lxy6z_D2y_aa = I_ESP_Mx2y6z_Py_aa+ABY*I_ESP_Lxy6z_Py_aa;
    Double I_ESP_Lx7z_D2y_aa = I_ESP_Mxy7z_Py_aa+ABY*I_ESP_Lx7z_Py_aa;
    Double I_ESP_L8y_D2y_aa = I_ESP_M9y_Py_aa+ABY*I_ESP_L8y_Py_aa;
    Double I_ESP_L7yz_D2y_aa = I_ESP_M8yz_Py_aa+ABY*I_ESP_L7yz_Py_aa;
    Double I_ESP_L6y2z_D2y_aa = I_ESP_M7y2z_Py_aa+ABY*I_ESP_L6y2z_Py_aa;
    Double I_ESP_L5y3z_D2y_aa = I_ESP_M6y3z_Py_aa+ABY*I_ESP_L5y3z_Py_aa;
    Double I_ESP_L4y4z_D2y_aa = I_ESP_M5y4z_Py_aa+ABY*I_ESP_L4y4z_Py_aa;
    Double I_ESP_L3y5z_D2y_aa = I_ESP_M4y5z_Py_aa+ABY*I_ESP_L3y5z_Py_aa;
    Double I_ESP_L2y6z_D2y_aa = I_ESP_M3y6z_Py_aa+ABY*I_ESP_L2y6z_Py_aa;
    Double I_ESP_Ly7z_D2y_aa = I_ESP_M2y7z_Py_aa+ABY*I_ESP_Ly7z_Py_aa;
    Double I_ESP_L8z_D2y_aa = I_ESP_My8z_Py_aa+ABY*I_ESP_L8z_Py_aa;
    Double I_ESP_L7xy_D2z_aa = I_ESP_M7xyz_Pz_aa+ABZ*I_ESP_L7xy_Pz_aa;
    Double I_ESP_L7xz_D2z_aa = I_ESP_M7x2z_Pz_aa+ABZ*I_ESP_L7xz_Pz_aa;
    Double I_ESP_L6x2y_D2z_aa = I_ESP_M6x2yz_Pz_aa+ABZ*I_ESP_L6x2y_Pz_aa;
    Double I_ESP_L6xyz_D2z_aa = I_ESP_M6xy2z_Pz_aa+ABZ*I_ESP_L6xyz_Pz_aa;
    Double I_ESP_L6x2z_D2z_aa = I_ESP_M6x3z_Pz_aa+ABZ*I_ESP_L6x2z_Pz_aa;
    Double I_ESP_L5x3y_D2z_aa = I_ESP_M5x3yz_Pz_aa+ABZ*I_ESP_L5x3y_Pz_aa;
    Double I_ESP_L5x2yz_D2z_aa = I_ESP_M5x2y2z_Pz_aa+ABZ*I_ESP_L5x2yz_Pz_aa;
    Double I_ESP_L5xy2z_D2z_aa = I_ESP_M5xy3z_Pz_aa+ABZ*I_ESP_L5xy2z_Pz_aa;
    Double I_ESP_L5x3z_D2z_aa = I_ESP_M5x4z_Pz_aa+ABZ*I_ESP_L5x3z_Pz_aa;
    Double I_ESP_L4x4y_D2z_aa = I_ESP_M4x4yz_Pz_aa+ABZ*I_ESP_L4x4y_Pz_aa;
    Double I_ESP_L4x3yz_D2z_aa = I_ESP_M4x3y2z_Pz_aa+ABZ*I_ESP_L4x3yz_Pz_aa;
    Double I_ESP_L4x2y2z_D2z_aa = I_ESP_M4x2y3z_Pz_aa+ABZ*I_ESP_L4x2y2z_Pz_aa;
    Double I_ESP_L4xy3z_D2z_aa = I_ESP_M4xy4z_Pz_aa+ABZ*I_ESP_L4xy3z_Pz_aa;
    Double I_ESP_L4x4z_D2z_aa = I_ESP_M4x5z_Pz_aa+ABZ*I_ESP_L4x4z_Pz_aa;
    Double I_ESP_L3x5y_D2z_aa = I_ESP_M3x5yz_Pz_aa+ABZ*I_ESP_L3x5y_Pz_aa;
    Double I_ESP_L3x4yz_D2z_aa = I_ESP_M3x4y2z_Pz_aa+ABZ*I_ESP_L3x4yz_Pz_aa;
    Double I_ESP_L3x3y2z_D2z_aa = I_ESP_M3x3y3z_Pz_aa+ABZ*I_ESP_L3x3y2z_Pz_aa;
    Double I_ESP_L3x2y3z_D2z_aa = I_ESP_M3x2y4z_Pz_aa+ABZ*I_ESP_L3x2y3z_Pz_aa;
    Double I_ESP_L3xy4z_D2z_aa = I_ESP_M3xy5z_Pz_aa+ABZ*I_ESP_L3xy4z_Pz_aa;
    Double I_ESP_L3x5z_D2z_aa = I_ESP_M3x6z_Pz_aa+ABZ*I_ESP_L3x5z_Pz_aa;
    Double I_ESP_L2x6y_D2z_aa = I_ESP_M2x6yz_Pz_aa+ABZ*I_ESP_L2x6y_Pz_aa;
    Double I_ESP_L2x5yz_D2z_aa = I_ESP_M2x5y2z_Pz_aa+ABZ*I_ESP_L2x5yz_Pz_aa;
    Double I_ESP_L2x4y2z_D2z_aa = I_ESP_M2x4y3z_Pz_aa+ABZ*I_ESP_L2x4y2z_Pz_aa;
    Double I_ESP_L2x3y3z_D2z_aa = I_ESP_M2x3y4z_Pz_aa+ABZ*I_ESP_L2x3y3z_Pz_aa;
    Double I_ESP_L2x2y4z_D2z_aa = I_ESP_M2x2y5z_Pz_aa+ABZ*I_ESP_L2x2y4z_Pz_aa;
    Double I_ESP_L2xy5z_D2z_aa = I_ESP_M2xy6z_Pz_aa+ABZ*I_ESP_L2xy5z_Pz_aa;
    Double I_ESP_L2x6z_D2z_aa = I_ESP_M2x7z_Pz_aa+ABZ*I_ESP_L2x6z_Pz_aa;
    Double I_ESP_Lx7y_D2z_aa = I_ESP_Mx7yz_Pz_aa+ABZ*I_ESP_Lx7y_Pz_aa;
    Double I_ESP_Lx6yz_D2z_aa = I_ESP_Mx6y2z_Pz_aa+ABZ*I_ESP_Lx6yz_Pz_aa;
    Double I_ESP_Lx5y2z_D2z_aa = I_ESP_Mx5y3z_Pz_aa+ABZ*I_ESP_Lx5y2z_Pz_aa;
    Double I_ESP_Lx4y3z_D2z_aa = I_ESP_Mx4y4z_Pz_aa+ABZ*I_ESP_Lx4y3z_Pz_aa;
    Double I_ESP_Lx3y4z_D2z_aa = I_ESP_Mx3y5z_Pz_aa+ABZ*I_ESP_Lx3y4z_Pz_aa;
    Double I_ESP_Lx2y5z_D2z_aa = I_ESP_Mx2y6z_Pz_aa+ABZ*I_ESP_Lx2y5z_Pz_aa;
    Double I_ESP_Lxy6z_D2z_aa = I_ESP_Mxy7z_Pz_aa+ABZ*I_ESP_Lxy6z_Pz_aa;
    Double I_ESP_Lx7z_D2z_aa = I_ESP_Mx8z_Pz_aa+ABZ*I_ESP_Lx7z_Pz_aa;
    Double I_ESP_L7yz_D2z_aa = I_ESP_M7y2z_Pz_aa+ABZ*I_ESP_L7yz_Pz_aa;
    Double I_ESP_L6y2z_D2z_aa = I_ESP_M6y3z_Pz_aa+ABZ*I_ESP_L6y2z_Pz_aa;
    Double I_ESP_L5y3z_D2z_aa = I_ESP_M5y4z_Pz_aa+ABZ*I_ESP_L5y3z_Pz_aa;
    Double I_ESP_L4y4z_D2z_aa = I_ESP_M4y5z_Pz_aa+ABZ*I_ESP_L4y4z_Pz_aa;
    Double I_ESP_L3y5z_D2z_aa = I_ESP_M3y6z_Pz_aa+ABZ*I_ESP_L3y5z_Pz_aa;
    Double I_ESP_L2y6z_D2z_aa = I_ESP_M2y7z_Pz_aa+ABZ*I_ESP_L2y6z_Pz_aa;
    Double I_ESP_Ly7z_D2z_aa = I_ESP_My8z_Pz_aa+ABZ*I_ESP_Ly7z_Pz_aa;
    Double I_ESP_L8z_D2z_aa = I_ESP_M9z_Pz_aa+ABZ*I_ESP_L8z_Pz_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_K_F_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 105 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_D_aa
     * RHS shell quartet name: SQ_ESP_K_D_aa
     ************************************************************/
    Double I_ESP_K7x_F3x_aa = I_ESP_L8x_D2x_aa+ABX*I_ESP_K7x_D2x_aa;
    Double I_ESP_K6xy_F3x_aa = I_ESP_L7xy_D2x_aa+ABX*I_ESP_K6xy_D2x_aa;
    Double I_ESP_K6xz_F3x_aa = I_ESP_L7xz_D2x_aa+ABX*I_ESP_K6xz_D2x_aa;
    Double I_ESP_K5x2y_F3x_aa = I_ESP_L6x2y_D2x_aa+ABX*I_ESP_K5x2y_D2x_aa;
    Double I_ESP_K5xyz_F3x_aa = I_ESP_L6xyz_D2x_aa+ABX*I_ESP_K5xyz_D2x_aa;
    Double I_ESP_K5x2z_F3x_aa = I_ESP_L6x2z_D2x_aa+ABX*I_ESP_K5x2z_D2x_aa;
    Double I_ESP_K4x3y_F3x_aa = I_ESP_L5x3y_D2x_aa+ABX*I_ESP_K4x3y_D2x_aa;
    Double I_ESP_K4x2yz_F3x_aa = I_ESP_L5x2yz_D2x_aa+ABX*I_ESP_K4x2yz_D2x_aa;
    Double I_ESP_K4xy2z_F3x_aa = I_ESP_L5xy2z_D2x_aa+ABX*I_ESP_K4xy2z_D2x_aa;
    Double I_ESP_K4x3z_F3x_aa = I_ESP_L5x3z_D2x_aa+ABX*I_ESP_K4x3z_D2x_aa;
    Double I_ESP_K3x4y_F3x_aa = I_ESP_L4x4y_D2x_aa+ABX*I_ESP_K3x4y_D2x_aa;
    Double I_ESP_K3x3yz_F3x_aa = I_ESP_L4x3yz_D2x_aa+ABX*I_ESP_K3x3yz_D2x_aa;
    Double I_ESP_K3x2y2z_F3x_aa = I_ESP_L4x2y2z_D2x_aa+ABX*I_ESP_K3x2y2z_D2x_aa;
    Double I_ESP_K3xy3z_F3x_aa = I_ESP_L4xy3z_D2x_aa+ABX*I_ESP_K3xy3z_D2x_aa;
    Double I_ESP_K3x4z_F3x_aa = I_ESP_L4x4z_D2x_aa+ABX*I_ESP_K3x4z_D2x_aa;
    Double I_ESP_K2x5y_F3x_aa = I_ESP_L3x5y_D2x_aa+ABX*I_ESP_K2x5y_D2x_aa;
    Double I_ESP_K2x4yz_F3x_aa = I_ESP_L3x4yz_D2x_aa+ABX*I_ESP_K2x4yz_D2x_aa;
    Double I_ESP_K2x3y2z_F3x_aa = I_ESP_L3x3y2z_D2x_aa+ABX*I_ESP_K2x3y2z_D2x_aa;
    Double I_ESP_K2x2y3z_F3x_aa = I_ESP_L3x2y3z_D2x_aa+ABX*I_ESP_K2x2y3z_D2x_aa;
    Double I_ESP_K2xy4z_F3x_aa = I_ESP_L3xy4z_D2x_aa+ABX*I_ESP_K2xy4z_D2x_aa;
    Double I_ESP_K2x5z_F3x_aa = I_ESP_L3x5z_D2x_aa+ABX*I_ESP_K2x5z_D2x_aa;
    Double I_ESP_Kx6y_F3x_aa = I_ESP_L2x6y_D2x_aa+ABX*I_ESP_Kx6y_D2x_aa;
    Double I_ESP_Kx5yz_F3x_aa = I_ESP_L2x5yz_D2x_aa+ABX*I_ESP_Kx5yz_D2x_aa;
    Double I_ESP_Kx4y2z_F3x_aa = I_ESP_L2x4y2z_D2x_aa+ABX*I_ESP_Kx4y2z_D2x_aa;
    Double I_ESP_Kx3y3z_F3x_aa = I_ESP_L2x3y3z_D2x_aa+ABX*I_ESP_Kx3y3z_D2x_aa;
    Double I_ESP_Kx2y4z_F3x_aa = I_ESP_L2x2y4z_D2x_aa+ABX*I_ESP_Kx2y4z_D2x_aa;
    Double I_ESP_Kxy5z_F3x_aa = I_ESP_L2xy5z_D2x_aa+ABX*I_ESP_Kxy5z_D2x_aa;
    Double I_ESP_Kx6z_F3x_aa = I_ESP_L2x6z_D2x_aa+ABX*I_ESP_Kx6z_D2x_aa;
    Double I_ESP_K7y_F3x_aa = I_ESP_Lx7y_D2x_aa+ABX*I_ESP_K7y_D2x_aa;
    Double I_ESP_K6yz_F3x_aa = I_ESP_Lx6yz_D2x_aa+ABX*I_ESP_K6yz_D2x_aa;
    Double I_ESP_K5y2z_F3x_aa = I_ESP_Lx5y2z_D2x_aa+ABX*I_ESP_K5y2z_D2x_aa;
    Double I_ESP_K4y3z_F3x_aa = I_ESP_Lx4y3z_D2x_aa+ABX*I_ESP_K4y3z_D2x_aa;
    Double I_ESP_K3y4z_F3x_aa = I_ESP_Lx3y4z_D2x_aa+ABX*I_ESP_K3y4z_D2x_aa;
    Double I_ESP_K2y5z_F3x_aa = I_ESP_Lx2y5z_D2x_aa+ABX*I_ESP_K2y5z_D2x_aa;
    Double I_ESP_Ky6z_F3x_aa = I_ESP_Lxy6z_D2x_aa+ABX*I_ESP_Ky6z_D2x_aa;
    Double I_ESP_K7z_F3x_aa = I_ESP_Lx7z_D2x_aa+ABX*I_ESP_K7z_D2x_aa;
    Double I_ESP_K6xy_F2xy_aa = I_ESP_L6x2y_D2x_aa+ABY*I_ESP_K6xy_D2x_aa;
    Double I_ESP_K6xz_F2xy_aa = I_ESP_L6xyz_D2x_aa+ABY*I_ESP_K6xz_D2x_aa;
    Double I_ESP_K5x2y_F2xy_aa = I_ESP_L5x3y_D2x_aa+ABY*I_ESP_K5x2y_D2x_aa;
    Double I_ESP_K5xyz_F2xy_aa = I_ESP_L5x2yz_D2x_aa+ABY*I_ESP_K5xyz_D2x_aa;
    Double I_ESP_K5x2z_F2xy_aa = I_ESP_L5xy2z_D2x_aa+ABY*I_ESP_K5x2z_D2x_aa;
    Double I_ESP_K4x3y_F2xy_aa = I_ESP_L4x4y_D2x_aa+ABY*I_ESP_K4x3y_D2x_aa;
    Double I_ESP_K4x2yz_F2xy_aa = I_ESP_L4x3yz_D2x_aa+ABY*I_ESP_K4x2yz_D2x_aa;
    Double I_ESP_K4xy2z_F2xy_aa = I_ESP_L4x2y2z_D2x_aa+ABY*I_ESP_K4xy2z_D2x_aa;
    Double I_ESP_K4x3z_F2xy_aa = I_ESP_L4xy3z_D2x_aa+ABY*I_ESP_K4x3z_D2x_aa;
    Double I_ESP_K3x4y_F2xy_aa = I_ESP_L3x5y_D2x_aa+ABY*I_ESP_K3x4y_D2x_aa;
    Double I_ESP_K3x3yz_F2xy_aa = I_ESP_L3x4yz_D2x_aa+ABY*I_ESP_K3x3yz_D2x_aa;
    Double I_ESP_K3x2y2z_F2xy_aa = I_ESP_L3x3y2z_D2x_aa+ABY*I_ESP_K3x2y2z_D2x_aa;
    Double I_ESP_K3xy3z_F2xy_aa = I_ESP_L3x2y3z_D2x_aa+ABY*I_ESP_K3xy3z_D2x_aa;
    Double I_ESP_K3x4z_F2xy_aa = I_ESP_L3xy4z_D2x_aa+ABY*I_ESP_K3x4z_D2x_aa;
    Double I_ESP_K2x5y_F2xy_aa = I_ESP_L2x6y_D2x_aa+ABY*I_ESP_K2x5y_D2x_aa;
    Double I_ESP_K2x4yz_F2xy_aa = I_ESP_L2x5yz_D2x_aa+ABY*I_ESP_K2x4yz_D2x_aa;
    Double I_ESP_K2x3y2z_F2xy_aa = I_ESP_L2x4y2z_D2x_aa+ABY*I_ESP_K2x3y2z_D2x_aa;
    Double I_ESP_K2x2y3z_F2xy_aa = I_ESP_L2x3y3z_D2x_aa+ABY*I_ESP_K2x2y3z_D2x_aa;
    Double I_ESP_K2xy4z_F2xy_aa = I_ESP_L2x2y4z_D2x_aa+ABY*I_ESP_K2xy4z_D2x_aa;
    Double I_ESP_K2x5z_F2xy_aa = I_ESP_L2xy5z_D2x_aa+ABY*I_ESP_K2x5z_D2x_aa;
    Double I_ESP_Kx6y_F2xy_aa = I_ESP_Lx7y_D2x_aa+ABY*I_ESP_Kx6y_D2x_aa;
    Double I_ESP_Kx5yz_F2xy_aa = I_ESP_Lx6yz_D2x_aa+ABY*I_ESP_Kx5yz_D2x_aa;
    Double I_ESP_Kx4y2z_F2xy_aa = I_ESP_Lx5y2z_D2x_aa+ABY*I_ESP_Kx4y2z_D2x_aa;
    Double I_ESP_Kx3y3z_F2xy_aa = I_ESP_Lx4y3z_D2x_aa+ABY*I_ESP_Kx3y3z_D2x_aa;
    Double I_ESP_Kx2y4z_F2xy_aa = I_ESP_Lx3y4z_D2x_aa+ABY*I_ESP_Kx2y4z_D2x_aa;
    Double I_ESP_Kxy5z_F2xy_aa = I_ESP_Lx2y5z_D2x_aa+ABY*I_ESP_Kxy5z_D2x_aa;
    Double I_ESP_Kx6z_F2xy_aa = I_ESP_Lxy6z_D2x_aa+ABY*I_ESP_Kx6z_D2x_aa;
    Double I_ESP_K7y_F2xy_aa = I_ESP_L8y_D2x_aa+ABY*I_ESP_K7y_D2x_aa;
    Double I_ESP_K6yz_F2xy_aa = I_ESP_L7yz_D2x_aa+ABY*I_ESP_K6yz_D2x_aa;
    Double I_ESP_K5y2z_F2xy_aa = I_ESP_L6y2z_D2x_aa+ABY*I_ESP_K5y2z_D2x_aa;
    Double I_ESP_K4y3z_F2xy_aa = I_ESP_L5y3z_D2x_aa+ABY*I_ESP_K4y3z_D2x_aa;
    Double I_ESP_K3y4z_F2xy_aa = I_ESP_L4y4z_D2x_aa+ABY*I_ESP_K3y4z_D2x_aa;
    Double I_ESP_K2y5z_F2xy_aa = I_ESP_L3y5z_D2x_aa+ABY*I_ESP_K2y5z_D2x_aa;
    Double I_ESP_Ky6z_F2xy_aa = I_ESP_L2y6z_D2x_aa+ABY*I_ESP_Ky6z_D2x_aa;
    Double I_ESP_K7z_F2xy_aa = I_ESP_Ly7z_D2x_aa+ABY*I_ESP_K7z_D2x_aa;
    Double I_ESP_K6xz_F2xz_aa = I_ESP_L6x2z_D2x_aa+ABZ*I_ESP_K6xz_D2x_aa;
    Double I_ESP_K5xyz_F2xz_aa = I_ESP_L5xy2z_D2x_aa+ABZ*I_ESP_K5xyz_D2x_aa;
    Double I_ESP_K5x2z_F2xz_aa = I_ESP_L5x3z_D2x_aa+ABZ*I_ESP_K5x2z_D2x_aa;
    Double I_ESP_K4x2yz_F2xz_aa = I_ESP_L4x2y2z_D2x_aa+ABZ*I_ESP_K4x2yz_D2x_aa;
    Double I_ESP_K4xy2z_F2xz_aa = I_ESP_L4xy3z_D2x_aa+ABZ*I_ESP_K4xy2z_D2x_aa;
    Double I_ESP_K4x3z_F2xz_aa = I_ESP_L4x4z_D2x_aa+ABZ*I_ESP_K4x3z_D2x_aa;
    Double I_ESP_K3x3yz_F2xz_aa = I_ESP_L3x3y2z_D2x_aa+ABZ*I_ESP_K3x3yz_D2x_aa;
    Double I_ESP_K3x2y2z_F2xz_aa = I_ESP_L3x2y3z_D2x_aa+ABZ*I_ESP_K3x2y2z_D2x_aa;
    Double I_ESP_K3xy3z_F2xz_aa = I_ESP_L3xy4z_D2x_aa+ABZ*I_ESP_K3xy3z_D2x_aa;
    Double I_ESP_K3x4z_F2xz_aa = I_ESP_L3x5z_D2x_aa+ABZ*I_ESP_K3x4z_D2x_aa;
    Double I_ESP_K2x4yz_F2xz_aa = I_ESP_L2x4y2z_D2x_aa+ABZ*I_ESP_K2x4yz_D2x_aa;
    Double I_ESP_K2x3y2z_F2xz_aa = I_ESP_L2x3y3z_D2x_aa+ABZ*I_ESP_K2x3y2z_D2x_aa;
    Double I_ESP_K2x2y3z_F2xz_aa = I_ESP_L2x2y4z_D2x_aa+ABZ*I_ESP_K2x2y3z_D2x_aa;
    Double I_ESP_K2xy4z_F2xz_aa = I_ESP_L2xy5z_D2x_aa+ABZ*I_ESP_K2xy4z_D2x_aa;
    Double I_ESP_K2x5z_F2xz_aa = I_ESP_L2x6z_D2x_aa+ABZ*I_ESP_K2x5z_D2x_aa;
    Double I_ESP_Kx5yz_F2xz_aa = I_ESP_Lx5y2z_D2x_aa+ABZ*I_ESP_Kx5yz_D2x_aa;
    Double I_ESP_Kx4y2z_F2xz_aa = I_ESP_Lx4y3z_D2x_aa+ABZ*I_ESP_Kx4y2z_D2x_aa;
    Double I_ESP_Kx3y3z_F2xz_aa = I_ESP_Lx3y4z_D2x_aa+ABZ*I_ESP_Kx3y3z_D2x_aa;
    Double I_ESP_Kx2y4z_F2xz_aa = I_ESP_Lx2y5z_D2x_aa+ABZ*I_ESP_Kx2y4z_D2x_aa;
    Double I_ESP_Kxy5z_F2xz_aa = I_ESP_Lxy6z_D2x_aa+ABZ*I_ESP_Kxy5z_D2x_aa;
    Double I_ESP_Kx6z_F2xz_aa = I_ESP_Lx7z_D2x_aa+ABZ*I_ESP_Kx6z_D2x_aa;
    Double I_ESP_K6yz_F2xz_aa = I_ESP_L6y2z_D2x_aa+ABZ*I_ESP_K6yz_D2x_aa;
    Double I_ESP_K5y2z_F2xz_aa = I_ESP_L5y3z_D2x_aa+ABZ*I_ESP_K5y2z_D2x_aa;
    Double I_ESP_K4y3z_F2xz_aa = I_ESP_L4y4z_D2x_aa+ABZ*I_ESP_K4y3z_D2x_aa;
    Double I_ESP_K3y4z_F2xz_aa = I_ESP_L3y5z_D2x_aa+ABZ*I_ESP_K3y4z_D2x_aa;
    Double I_ESP_K2y5z_F2xz_aa = I_ESP_L2y6z_D2x_aa+ABZ*I_ESP_K2y5z_D2x_aa;
    Double I_ESP_Ky6z_F2xz_aa = I_ESP_Ly7z_D2x_aa+ABZ*I_ESP_Ky6z_D2x_aa;
    Double I_ESP_K7z_F2xz_aa = I_ESP_L8z_D2x_aa+ABZ*I_ESP_K7z_D2x_aa;
    Double I_ESP_K6xz_Fx2y_aa = I_ESP_L7xz_D2y_aa+ABX*I_ESP_K6xz_D2y_aa;
    Double I_ESP_K5xyz_Fx2y_aa = I_ESP_L6xyz_D2y_aa+ABX*I_ESP_K5xyz_D2y_aa;
    Double I_ESP_K5x2z_Fx2y_aa = I_ESP_L6x2z_D2y_aa+ABX*I_ESP_K5x2z_D2y_aa;
    Double I_ESP_K4x2yz_Fx2y_aa = I_ESP_L5x2yz_D2y_aa+ABX*I_ESP_K4x2yz_D2y_aa;
    Double I_ESP_K4xy2z_Fx2y_aa = I_ESP_L5xy2z_D2y_aa+ABX*I_ESP_K4xy2z_D2y_aa;
    Double I_ESP_K4x3z_Fx2y_aa = I_ESP_L5x3z_D2y_aa+ABX*I_ESP_K4x3z_D2y_aa;
    Double I_ESP_K3x3yz_Fx2y_aa = I_ESP_L4x3yz_D2y_aa+ABX*I_ESP_K3x3yz_D2y_aa;
    Double I_ESP_K3x2y2z_Fx2y_aa = I_ESP_L4x2y2z_D2y_aa+ABX*I_ESP_K3x2y2z_D2y_aa;
    Double I_ESP_K3xy3z_Fx2y_aa = I_ESP_L4xy3z_D2y_aa+ABX*I_ESP_K3xy3z_D2y_aa;
    Double I_ESP_K3x4z_Fx2y_aa = I_ESP_L4x4z_D2y_aa+ABX*I_ESP_K3x4z_D2y_aa;
    Double I_ESP_K2x4yz_Fx2y_aa = I_ESP_L3x4yz_D2y_aa+ABX*I_ESP_K2x4yz_D2y_aa;
    Double I_ESP_K2x3y2z_Fx2y_aa = I_ESP_L3x3y2z_D2y_aa+ABX*I_ESP_K2x3y2z_D2y_aa;
    Double I_ESP_K2x2y3z_Fx2y_aa = I_ESP_L3x2y3z_D2y_aa+ABX*I_ESP_K2x2y3z_D2y_aa;
    Double I_ESP_K2xy4z_Fx2y_aa = I_ESP_L3xy4z_D2y_aa+ABX*I_ESP_K2xy4z_D2y_aa;
    Double I_ESP_K2x5z_Fx2y_aa = I_ESP_L3x5z_D2y_aa+ABX*I_ESP_K2x5z_D2y_aa;
    Double I_ESP_Kx5yz_Fx2y_aa = I_ESP_L2x5yz_D2y_aa+ABX*I_ESP_Kx5yz_D2y_aa;
    Double I_ESP_Kx4y2z_Fx2y_aa = I_ESP_L2x4y2z_D2y_aa+ABX*I_ESP_Kx4y2z_D2y_aa;
    Double I_ESP_Kx3y3z_Fx2y_aa = I_ESP_L2x3y3z_D2y_aa+ABX*I_ESP_Kx3y3z_D2y_aa;
    Double I_ESP_Kx2y4z_Fx2y_aa = I_ESP_L2x2y4z_D2y_aa+ABX*I_ESP_Kx2y4z_D2y_aa;
    Double I_ESP_Kxy5z_Fx2y_aa = I_ESP_L2xy5z_D2y_aa+ABX*I_ESP_Kxy5z_D2y_aa;
    Double I_ESP_Kx6z_Fx2y_aa = I_ESP_L2x6z_D2y_aa+ABX*I_ESP_Kx6z_D2y_aa;
    Double I_ESP_K6yz_Fx2y_aa = I_ESP_Lx6yz_D2y_aa+ABX*I_ESP_K6yz_D2y_aa;
    Double I_ESP_K5y2z_Fx2y_aa = I_ESP_Lx5y2z_D2y_aa+ABX*I_ESP_K5y2z_D2y_aa;
    Double I_ESP_K4y3z_Fx2y_aa = I_ESP_Lx4y3z_D2y_aa+ABX*I_ESP_K4y3z_D2y_aa;
    Double I_ESP_K3y4z_Fx2y_aa = I_ESP_Lx3y4z_D2y_aa+ABX*I_ESP_K3y4z_D2y_aa;
    Double I_ESP_K2y5z_Fx2y_aa = I_ESP_Lx2y5z_D2y_aa+ABX*I_ESP_K2y5z_D2y_aa;
    Double I_ESP_Ky6z_Fx2y_aa = I_ESP_Lxy6z_D2y_aa+ABX*I_ESP_Ky6z_D2y_aa;
    Double I_ESP_K7z_Fx2y_aa = I_ESP_Lx7z_D2y_aa+ABX*I_ESP_K7z_D2y_aa;
    Double I_ESP_K6xy_Fx2z_aa = I_ESP_L7xy_D2z_aa+ABX*I_ESP_K6xy_D2z_aa;
    Double I_ESP_K5x2y_Fx2z_aa = I_ESP_L6x2y_D2z_aa+ABX*I_ESP_K5x2y_D2z_aa;
    Double I_ESP_K5xyz_Fx2z_aa = I_ESP_L6xyz_D2z_aa+ABX*I_ESP_K5xyz_D2z_aa;
    Double I_ESP_K4x3y_Fx2z_aa = I_ESP_L5x3y_D2z_aa+ABX*I_ESP_K4x3y_D2z_aa;
    Double I_ESP_K4x2yz_Fx2z_aa = I_ESP_L5x2yz_D2z_aa+ABX*I_ESP_K4x2yz_D2z_aa;
    Double I_ESP_K4xy2z_Fx2z_aa = I_ESP_L5xy2z_D2z_aa+ABX*I_ESP_K4xy2z_D2z_aa;
    Double I_ESP_K3x4y_Fx2z_aa = I_ESP_L4x4y_D2z_aa+ABX*I_ESP_K3x4y_D2z_aa;
    Double I_ESP_K3x3yz_Fx2z_aa = I_ESP_L4x3yz_D2z_aa+ABX*I_ESP_K3x3yz_D2z_aa;
    Double I_ESP_K3x2y2z_Fx2z_aa = I_ESP_L4x2y2z_D2z_aa+ABX*I_ESP_K3x2y2z_D2z_aa;
    Double I_ESP_K3xy3z_Fx2z_aa = I_ESP_L4xy3z_D2z_aa+ABX*I_ESP_K3xy3z_D2z_aa;
    Double I_ESP_K2x5y_Fx2z_aa = I_ESP_L3x5y_D2z_aa+ABX*I_ESP_K2x5y_D2z_aa;
    Double I_ESP_K2x4yz_Fx2z_aa = I_ESP_L3x4yz_D2z_aa+ABX*I_ESP_K2x4yz_D2z_aa;
    Double I_ESP_K2x3y2z_Fx2z_aa = I_ESP_L3x3y2z_D2z_aa+ABX*I_ESP_K2x3y2z_D2z_aa;
    Double I_ESP_K2x2y3z_Fx2z_aa = I_ESP_L3x2y3z_D2z_aa+ABX*I_ESP_K2x2y3z_D2z_aa;
    Double I_ESP_K2xy4z_Fx2z_aa = I_ESP_L3xy4z_D2z_aa+ABX*I_ESP_K2xy4z_D2z_aa;
    Double I_ESP_Kx6y_Fx2z_aa = I_ESP_L2x6y_D2z_aa+ABX*I_ESP_Kx6y_D2z_aa;
    Double I_ESP_Kx5yz_Fx2z_aa = I_ESP_L2x5yz_D2z_aa+ABX*I_ESP_Kx5yz_D2z_aa;
    Double I_ESP_Kx4y2z_Fx2z_aa = I_ESP_L2x4y2z_D2z_aa+ABX*I_ESP_Kx4y2z_D2z_aa;
    Double I_ESP_Kx3y3z_Fx2z_aa = I_ESP_L2x3y3z_D2z_aa+ABX*I_ESP_Kx3y3z_D2z_aa;
    Double I_ESP_Kx2y4z_Fx2z_aa = I_ESP_L2x2y4z_D2z_aa+ABX*I_ESP_Kx2y4z_D2z_aa;
    Double I_ESP_Kxy5z_Fx2z_aa = I_ESP_L2xy5z_D2z_aa+ABX*I_ESP_Kxy5z_D2z_aa;
    Double I_ESP_K7y_Fx2z_aa = I_ESP_Lx7y_D2z_aa+ABX*I_ESP_K7y_D2z_aa;
    Double I_ESP_K6yz_Fx2z_aa = I_ESP_Lx6yz_D2z_aa+ABX*I_ESP_K6yz_D2z_aa;
    Double I_ESP_K5y2z_Fx2z_aa = I_ESP_Lx5y2z_D2z_aa+ABX*I_ESP_K5y2z_D2z_aa;
    Double I_ESP_K4y3z_Fx2z_aa = I_ESP_Lx4y3z_D2z_aa+ABX*I_ESP_K4y3z_D2z_aa;
    Double I_ESP_K3y4z_Fx2z_aa = I_ESP_Lx3y4z_D2z_aa+ABX*I_ESP_K3y4z_D2z_aa;
    Double I_ESP_K2y5z_Fx2z_aa = I_ESP_Lx2y5z_D2z_aa+ABX*I_ESP_K2y5z_D2z_aa;
    Double I_ESP_Ky6z_Fx2z_aa = I_ESP_Lxy6z_D2z_aa+ABX*I_ESP_Ky6z_D2z_aa;
    Double I_ESP_K7x_F3y_aa = I_ESP_L7xy_D2y_aa+ABY*I_ESP_K7x_D2y_aa;
    Double I_ESP_K6xy_F3y_aa = I_ESP_L6x2y_D2y_aa+ABY*I_ESP_K6xy_D2y_aa;
    Double I_ESP_K6xz_F3y_aa = I_ESP_L6xyz_D2y_aa+ABY*I_ESP_K6xz_D2y_aa;
    Double I_ESP_K5x2y_F3y_aa = I_ESP_L5x3y_D2y_aa+ABY*I_ESP_K5x2y_D2y_aa;
    Double I_ESP_K5xyz_F3y_aa = I_ESP_L5x2yz_D2y_aa+ABY*I_ESP_K5xyz_D2y_aa;
    Double I_ESP_K5x2z_F3y_aa = I_ESP_L5xy2z_D2y_aa+ABY*I_ESP_K5x2z_D2y_aa;
    Double I_ESP_K4x3y_F3y_aa = I_ESP_L4x4y_D2y_aa+ABY*I_ESP_K4x3y_D2y_aa;
    Double I_ESP_K4x2yz_F3y_aa = I_ESP_L4x3yz_D2y_aa+ABY*I_ESP_K4x2yz_D2y_aa;
    Double I_ESP_K4xy2z_F3y_aa = I_ESP_L4x2y2z_D2y_aa+ABY*I_ESP_K4xy2z_D2y_aa;
    Double I_ESP_K4x3z_F3y_aa = I_ESP_L4xy3z_D2y_aa+ABY*I_ESP_K4x3z_D2y_aa;
    Double I_ESP_K3x4y_F3y_aa = I_ESP_L3x5y_D2y_aa+ABY*I_ESP_K3x4y_D2y_aa;
    Double I_ESP_K3x3yz_F3y_aa = I_ESP_L3x4yz_D2y_aa+ABY*I_ESP_K3x3yz_D2y_aa;
    Double I_ESP_K3x2y2z_F3y_aa = I_ESP_L3x3y2z_D2y_aa+ABY*I_ESP_K3x2y2z_D2y_aa;
    Double I_ESP_K3xy3z_F3y_aa = I_ESP_L3x2y3z_D2y_aa+ABY*I_ESP_K3xy3z_D2y_aa;
    Double I_ESP_K3x4z_F3y_aa = I_ESP_L3xy4z_D2y_aa+ABY*I_ESP_K3x4z_D2y_aa;
    Double I_ESP_K2x5y_F3y_aa = I_ESP_L2x6y_D2y_aa+ABY*I_ESP_K2x5y_D2y_aa;
    Double I_ESP_K2x4yz_F3y_aa = I_ESP_L2x5yz_D2y_aa+ABY*I_ESP_K2x4yz_D2y_aa;
    Double I_ESP_K2x3y2z_F3y_aa = I_ESP_L2x4y2z_D2y_aa+ABY*I_ESP_K2x3y2z_D2y_aa;
    Double I_ESP_K2x2y3z_F3y_aa = I_ESP_L2x3y3z_D2y_aa+ABY*I_ESP_K2x2y3z_D2y_aa;
    Double I_ESP_K2xy4z_F3y_aa = I_ESP_L2x2y4z_D2y_aa+ABY*I_ESP_K2xy4z_D2y_aa;
    Double I_ESP_K2x5z_F3y_aa = I_ESP_L2xy5z_D2y_aa+ABY*I_ESP_K2x5z_D2y_aa;
    Double I_ESP_Kx6y_F3y_aa = I_ESP_Lx7y_D2y_aa+ABY*I_ESP_Kx6y_D2y_aa;
    Double I_ESP_Kx5yz_F3y_aa = I_ESP_Lx6yz_D2y_aa+ABY*I_ESP_Kx5yz_D2y_aa;
    Double I_ESP_Kx4y2z_F3y_aa = I_ESP_Lx5y2z_D2y_aa+ABY*I_ESP_Kx4y2z_D2y_aa;
    Double I_ESP_Kx3y3z_F3y_aa = I_ESP_Lx4y3z_D2y_aa+ABY*I_ESP_Kx3y3z_D2y_aa;
    Double I_ESP_Kx2y4z_F3y_aa = I_ESP_Lx3y4z_D2y_aa+ABY*I_ESP_Kx2y4z_D2y_aa;
    Double I_ESP_Kxy5z_F3y_aa = I_ESP_Lx2y5z_D2y_aa+ABY*I_ESP_Kxy5z_D2y_aa;
    Double I_ESP_Kx6z_F3y_aa = I_ESP_Lxy6z_D2y_aa+ABY*I_ESP_Kx6z_D2y_aa;
    Double I_ESP_K7y_F3y_aa = I_ESP_L8y_D2y_aa+ABY*I_ESP_K7y_D2y_aa;
    Double I_ESP_K6yz_F3y_aa = I_ESP_L7yz_D2y_aa+ABY*I_ESP_K6yz_D2y_aa;
    Double I_ESP_K5y2z_F3y_aa = I_ESP_L6y2z_D2y_aa+ABY*I_ESP_K5y2z_D2y_aa;
    Double I_ESP_K4y3z_F3y_aa = I_ESP_L5y3z_D2y_aa+ABY*I_ESP_K4y3z_D2y_aa;
    Double I_ESP_K3y4z_F3y_aa = I_ESP_L4y4z_D2y_aa+ABY*I_ESP_K3y4z_D2y_aa;
    Double I_ESP_K2y5z_F3y_aa = I_ESP_L3y5z_D2y_aa+ABY*I_ESP_K2y5z_D2y_aa;
    Double I_ESP_Ky6z_F3y_aa = I_ESP_L2y6z_D2y_aa+ABY*I_ESP_Ky6z_D2y_aa;
    Double I_ESP_K7z_F3y_aa = I_ESP_Ly7z_D2y_aa+ABY*I_ESP_K7z_D2y_aa;
    Double I_ESP_K6xz_F2yz_aa = I_ESP_L6x2z_D2y_aa+ABZ*I_ESP_K6xz_D2y_aa;
    Double I_ESP_K5xyz_F2yz_aa = I_ESP_L5xy2z_D2y_aa+ABZ*I_ESP_K5xyz_D2y_aa;
    Double I_ESP_K5x2z_F2yz_aa = I_ESP_L5x3z_D2y_aa+ABZ*I_ESP_K5x2z_D2y_aa;
    Double I_ESP_K4x2yz_F2yz_aa = I_ESP_L4x2y2z_D2y_aa+ABZ*I_ESP_K4x2yz_D2y_aa;
    Double I_ESP_K4xy2z_F2yz_aa = I_ESP_L4xy3z_D2y_aa+ABZ*I_ESP_K4xy2z_D2y_aa;
    Double I_ESP_K4x3z_F2yz_aa = I_ESP_L4x4z_D2y_aa+ABZ*I_ESP_K4x3z_D2y_aa;
    Double I_ESP_K3x3yz_F2yz_aa = I_ESP_L3x3y2z_D2y_aa+ABZ*I_ESP_K3x3yz_D2y_aa;
    Double I_ESP_K3x2y2z_F2yz_aa = I_ESP_L3x2y3z_D2y_aa+ABZ*I_ESP_K3x2y2z_D2y_aa;
    Double I_ESP_K3xy3z_F2yz_aa = I_ESP_L3xy4z_D2y_aa+ABZ*I_ESP_K3xy3z_D2y_aa;
    Double I_ESP_K3x4z_F2yz_aa = I_ESP_L3x5z_D2y_aa+ABZ*I_ESP_K3x4z_D2y_aa;
    Double I_ESP_K2x4yz_F2yz_aa = I_ESP_L2x4y2z_D2y_aa+ABZ*I_ESP_K2x4yz_D2y_aa;
    Double I_ESP_K2x3y2z_F2yz_aa = I_ESP_L2x3y3z_D2y_aa+ABZ*I_ESP_K2x3y2z_D2y_aa;
    Double I_ESP_K2x2y3z_F2yz_aa = I_ESP_L2x2y4z_D2y_aa+ABZ*I_ESP_K2x2y3z_D2y_aa;
    Double I_ESP_K2xy4z_F2yz_aa = I_ESP_L2xy5z_D2y_aa+ABZ*I_ESP_K2xy4z_D2y_aa;
    Double I_ESP_K2x5z_F2yz_aa = I_ESP_L2x6z_D2y_aa+ABZ*I_ESP_K2x5z_D2y_aa;
    Double I_ESP_Kx5yz_F2yz_aa = I_ESP_Lx5y2z_D2y_aa+ABZ*I_ESP_Kx5yz_D2y_aa;
    Double I_ESP_Kx4y2z_F2yz_aa = I_ESP_Lx4y3z_D2y_aa+ABZ*I_ESP_Kx4y2z_D2y_aa;
    Double I_ESP_Kx3y3z_F2yz_aa = I_ESP_Lx3y4z_D2y_aa+ABZ*I_ESP_Kx3y3z_D2y_aa;
    Double I_ESP_Kx2y4z_F2yz_aa = I_ESP_Lx2y5z_D2y_aa+ABZ*I_ESP_Kx2y4z_D2y_aa;
    Double I_ESP_Kxy5z_F2yz_aa = I_ESP_Lxy6z_D2y_aa+ABZ*I_ESP_Kxy5z_D2y_aa;
    Double I_ESP_Kx6z_F2yz_aa = I_ESP_Lx7z_D2y_aa+ABZ*I_ESP_Kx6z_D2y_aa;
    Double I_ESP_K6yz_F2yz_aa = I_ESP_L6y2z_D2y_aa+ABZ*I_ESP_K6yz_D2y_aa;
    Double I_ESP_K5y2z_F2yz_aa = I_ESP_L5y3z_D2y_aa+ABZ*I_ESP_K5y2z_D2y_aa;
    Double I_ESP_K4y3z_F2yz_aa = I_ESP_L4y4z_D2y_aa+ABZ*I_ESP_K4y3z_D2y_aa;
    Double I_ESP_K3y4z_F2yz_aa = I_ESP_L3y5z_D2y_aa+ABZ*I_ESP_K3y4z_D2y_aa;
    Double I_ESP_K2y5z_F2yz_aa = I_ESP_L2y6z_D2y_aa+ABZ*I_ESP_K2y5z_D2y_aa;
    Double I_ESP_Ky6z_F2yz_aa = I_ESP_Ly7z_D2y_aa+ABZ*I_ESP_Ky6z_D2y_aa;
    Double I_ESP_K7z_F2yz_aa = I_ESP_L8z_D2y_aa+ABZ*I_ESP_K7z_D2y_aa;
    Double I_ESP_K7x_F3z_aa = I_ESP_L7xz_D2z_aa+ABZ*I_ESP_K7x_D2z_aa;
    Double I_ESP_K6xy_F3z_aa = I_ESP_L6xyz_D2z_aa+ABZ*I_ESP_K6xy_D2z_aa;
    Double I_ESP_K6xz_F3z_aa = I_ESP_L6x2z_D2z_aa+ABZ*I_ESP_K6xz_D2z_aa;
    Double I_ESP_K5x2y_F3z_aa = I_ESP_L5x2yz_D2z_aa+ABZ*I_ESP_K5x2y_D2z_aa;
    Double I_ESP_K5xyz_F3z_aa = I_ESP_L5xy2z_D2z_aa+ABZ*I_ESP_K5xyz_D2z_aa;
    Double I_ESP_K5x2z_F3z_aa = I_ESP_L5x3z_D2z_aa+ABZ*I_ESP_K5x2z_D2z_aa;
    Double I_ESP_K4x3y_F3z_aa = I_ESP_L4x3yz_D2z_aa+ABZ*I_ESP_K4x3y_D2z_aa;
    Double I_ESP_K4x2yz_F3z_aa = I_ESP_L4x2y2z_D2z_aa+ABZ*I_ESP_K4x2yz_D2z_aa;
    Double I_ESP_K4xy2z_F3z_aa = I_ESP_L4xy3z_D2z_aa+ABZ*I_ESP_K4xy2z_D2z_aa;
    Double I_ESP_K4x3z_F3z_aa = I_ESP_L4x4z_D2z_aa+ABZ*I_ESP_K4x3z_D2z_aa;
    Double I_ESP_K3x4y_F3z_aa = I_ESP_L3x4yz_D2z_aa+ABZ*I_ESP_K3x4y_D2z_aa;
    Double I_ESP_K3x3yz_F3z_aa = I_ESP_L3x3y2z_D2z_aa+ABZ*I_ESP_K3x3yz_D2z_aa;
    Double I_ESP_K3x2y2z_F3z_aa = I_ESP_L3x2y3z_D2z_aa+ABZ*I_ESP_K3x2y2z_D2z_aa;
    Double I_ESP_K3xy3z_F3z_aa = I_ESP_L3xy4z_D2z_aa+ABZ*I_ESP_K3xy3z_D2z_aa;
    Double I_ESP_K3x4z_F3z_aa = I_ESP_L3x5z_D2z_aa+ABZ*I_ESP_K3x4z_D2z_aa;
    Double I_ESP_K2x5y_F3z_aa = I_ESP_L2x5yz_D2z_aa+ABZ*I_ESP_K2x5y_D2z_aa;
    Double I_ESP_K2x4yz_F3z_aa = I_ESP_L2x4y2z_D2z_aa+ABZ*I_ESP_K2x4yz_D2z_aa;
    Double I_ESP_K2x3y2z_F3z_aa = I_ESP_L2x3y3z_D2z_aa+ABZ*I_ESP_K2x3y2z_D2z_aa;
    Double I_ESP_K2x2y3z_F3z_aa = I_ESP_L2x2y4z_D2z_aa+ABZ*I_ESP_K2x2y3z_D2z_aa;
    Double I_ESP_K2xy4z_F3z_aa = I_ESP_L2xy5z_D2z_aa+ABZ*I_ESP_K2xy4z_D2z_aa;
    Double I_ESP_K2x5z_F3z_aa = I_ESP_L2x6z_D2z_aa+ABZ*I_ESP_K2x5z_D2z_aa;
    Double I_ESP_Kx6y_F3z_aa = I_ESP_Lx6yz_D2z_aa+ABZ*I_ESP_Kx6y_D2z_aa;
    Double I_ESP_Kx5yz_F3z_aa = I_ESP_Lx5y2z_D2z_aa+ABZ*I_ESP_Kx5yz_D2z_aa;
    Double I_ESP_Kx4y2z_F3z_aa = I_ESP_Lx4y3z_D2z_aa+ABZ*I_ESP_Kx4y2z_D2z_aa;
    Double I_ESP_Kx3y3z_F3z_aa = I_ESP_Lx3y4z_D2z_aa+ABZ*I_ESP_Kx3y3z_D2z_aa;
    Double I_ESP_Kx2y4z_F3z_aa = I_ESP_Lx2y5z_D2z_aa+ABZ*I_ESP_Kx2y4z_D2z_aa;
    Double I_ESP_Kxy5z_F3z_aa = I_ESP_Lxy6z_D2z_aa+ABZ*I_ESP_Kxy5z_D2z_aa;
    Double I_ESP_Kx6z_F3z_aa = I_ESP_Lx7z_D2z_aa+ABZ*I_ESP_Kx6z_D2z_aa;
    Double I_ESP_K7y_F3z_aa = I_ESP_L7yz_D2z_aa+ABZ*I_ESP_K7y_D2z_aa;
    Double I_ESP_K6yz_F3z_aa = I_ESP_L6y2z_D2z_aa+ABZ*I_ESP_K6yz_D2z_aa;
    Double I_ESP_K5y2z_F3z_aa = I_ESP_L5y3z_D2z_aa+ABZ*I_ESP_K5y2z_D2z_aa;
    Double I_ESP_K4y3z_F3z_aa = I_ESP_L4y4z_D2z_aa+ABZ*I_ESP_K4y3z_D2z_aa;
    Double I_ESP_K3y4z_F3z_aa = I_ESP_L3y5z_D2z_aa+ABZ*I_ESP_K3y4z_D2z_aa;
    Double I_ESP_K2y5z_F3z_aa = I_ESP_L2y6z_D2z_aa+ABZ*I_ESP_K2y5z_D2z_aa;
    Double I_ESP_Ky6z_F3z_aa = I_ESP_Ly7z_D2z_aa+ABZ*I_ESP_Ky6z_D2z_aa;
    Double I_ESP_K7z_F3z_aa = I_ESP_L8z_D2z_aa+ABZ*I_ESP_K7z_D2z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_I_G_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_F_aa
     * RHS shell quartet name: SQ_ESP_I_F_aa
     ************************************************************/
    Double I_ESP_I6x_G4x_aa = I_ESP_K7x_F3x_aa+ABX*I_ESP_I6x_F3x_aa;
    Double I_ESP_I5xy_G4x_aa = I_ESP_K6xy_F3x_aa+ABX*I_ESP_I5xy_F3x_aa;
    Double I_ESP_I5xz_G4x_aa = I_ESP_K6xz_F3x_aa+ABX*I_ESP_I5xz_F3x_aa;
    Double I_ESP_I4x2y_G4x_aa = I_ESP_K5x2y_F3x_aa+ABX*I_ESP_I4x2y_F3x_aa;
    Double I_ESP_I4xyz_G4x_aa = I_ESP_K5xyz_F3x_aa+ABX*I_ESP_I4xyz_F3x_aa;
    Double I_ESP_I4x2z_G4x_aa = I_ESP_K5x2z_F3x_aa+ABX*I_ESP_I4x2z_F3x_aa;
    Double I_ESP_I3x3y_G4x_aa = I_ESP_K4x3y_F3x_aa+ABX*I_ESP_I3x3y_F3x_aa;
    Double I_ESP_I3x2yz_G4x_aa = I_ESP_K4x2yz_F3x_aa+ABX*I_ESP_I3x2yz_F3x_aa;
    Double I_ESP_I3xy2z_G4x_aa = I_ESP_K4xy2z_F3x_aa+ABX*I_ESP_I3xy2z_F3x_aa;
    Double I_ESP_I3x3z_G4x_aa = I_ESP_K4x3z_F3x_aa+ABX*I_ESP_I3x3z_F3x_aa;
    Double I_ESP_I2x4y_G4x_aa = I_ESP_K3x4y_F3x_aa+ABX*I_ESP_I2x4y_F3x_aa;
    Double I_ESP_I2x3yz_G4x_aa = I_ESP_K3x3yz_F3x_aa+ABX*I_ESP_I2x3yz_F3x_aa;
    Double I_ESP_I2x2y2z_G4x_aa = I_ESP_K3x2y2z_F3x_aa+ABX*I_ESP_I2x2y2z_F3x_aa;
    Double I_ESP_I2xy3z_G4x_aa = I_ESP_K3xy3z_F3x_aa+ABX*I_ESP_I2xy3z_F3x_aa;
    Double I_ESP_I2x4z_G4x_aa = I_ESP_K3x4z_F3x_aa+ABX*I_ESP_I2x4z_F3x_aa;
    Double I_ESP_Ix5y_G4x_aa = I_ESP_K2x5y_F3x_aa+ABX*I_ESP_Ix5y_F3x_aa;
    Double I_ESP_Ix4yz_G4x_aa = I_ESP_K2x4yz_F3x_aa+ABX*I_ESP_Ix4yz_F3x_aa;
    Double I_ESP_Ix3y2z_G4x_aa = I_ESP_K2x3y2z_F3x_aa+ABX*I_ESP_Ix3y2z_F3x_aa;
    Double I_ESP_Ix2y3z_G4x_aa = I_ESP_K2x2y3z_F3x_aa+ABX*I_ESP_Ix2y3z_F3x_aa;
    Double I_ESP_Ixy4z_G4x_aa = I_ESP_K2xy4z_F3x_aa+ABX*I_ESP_Ixy4z_F3x_aa;
    Double I_ESP_Ix5z_G4x_aa = I_ESP_K2x5z_F3x_aa+ABX*I_ESP_Ix5z_F3x_aa;
    Double I_ESP_I6y_G4x_aa = I_ESP_Kx6y_F3x_aa+ABX*I_ESP_I6y_F3x_aa;
    Double I_ESP_I5yz_G4x_aa = I_ESP_Kx5yz_F3x_aa+ABX*I_ESP_I5yz_F3x_aa;
    Double I_ESP_I4y2z_G4x_aa = I_ESP_Kx4y2z_F3x_aa+ABX*I_ESP_I4y2z_F3x_aa;
    Double I_ESP_I3y3z_G4x_aa = I_ESP_Kx3y3z_F3x_aa+ABX*I_ESP_I3y3z_F3x_aa;
    Double I_ESP_I2y4z_G4x_aa = I_ESP_Kx2y4z_F3x_aa+ABX*I_ESP_I2y4z_F3x_aa;
    Double I_ESP_Iy5z_G4x_aa = I_ESP_Kxy5z_F3x_aa+ABX*I_ESP_Iy5z_F3x_aa;
    Double I_ESP_I6z_G4x_aa = I_ESP_Kx6z_F3x_aa+ABX*I_ESP_I6z_F3x_aa;
    Double I_ESP_I6x_G3xy_aa = I_ESP_K6xy_F3x_aa+ABY*I_ESP_I6x_F3x_aa;
    Double I_ESP_I5xy_G3xy_aa = I_ESP_K5x2y_F3x_aa+ABY*I_ESP_I5xy_F3x_aa;
    Double I_ESP_I5xz_G3xy_aa = I_ESP_K5xyz_F3x_aa+ABY*I_ESP_I5xz_F3x_aa;
    Double I_ESP_I4x2y_G3xy_aa = I_ESP_K4x3y_F3x_aa+ABY*I_ESP_I4x2y_F3x_aa;
    Double I_ESP_I4xyz_G3xy_aa = I_ESP_K4x2yz_F3x_aa+ABY*I_ESP_I4xyz_F3x_aa;
    Double I_ESP_I4x2z_G3xy_aa = I_ESP_K4xy2z_F3x_aa+ABY*I_ESP_I4x2z_F3x_aa;
    Double I_ESP_I3x3y_G3xy_aa = I_ESP_K3x4y_F3x_aa+ABY*I_ESP_I3x3y_F3x_aa;
    Double I_ESP_I3x2yz_G3xy_aa = I_ESP_K3x3yz_F3x_aa+ABY*I_ESP_I3x2yz_F3x_aa;
    Double I_ESP_I3xy2z_G3xy_aa = I_ESP_K3x2y2z_F3x_aa+ABY*I_ESP_I3xy2z_F3x_aa;
    Double I_ESP_I3x3z_G3xy_aa = I_ESP_K3xy3z_F3x_aa+ABY*I_ESP_I3x3z_F3x_aa;
    Double I_ESP_I2x4y_G3xy_aa = I_ESP_K2x5y_F3x_aa+ABY*I_ESP_I2x4y_F3x_aa;
    Double I_ESP_I2x3yz_G3xy_aa = I_ESP_K2x4yz_F3x_aa+ABY*I_ESP_I2x3yz_F3x_aa;
    Double I_ESP_I2x2y2z_G3xy_aa = I_ESP_K2x3y2z_F3x_aa+ABY*I_ESP_I2x2y2z_F3x_aa;
    Double I_ESP_I2xy3z_G3xy_aa = I_ESP_K2x2y3z_F3x_aa+ABY*I_ESP_I2xy3z_F3x_aa;
    Double I_ESP_I2x4z_G3xy_aa = I_ESP_K2xy4z_F3x_aa+ABY*I_ESP_I2x4z_F3x_aa;
    Double I_ESP_Ix5y_G3xy_aa = I_ESP_Kx6y_F3x_aa+ABY*I_ESP_Ix5y_F3x_aa;
    Double I_ESP_Ix4yz_G3xy_aa = I_ESP_Kx5yz_F3x_aa+ABY*I_ESP_Ix4yz_F3x_aa;
    Double I_ESP_Ix3y2z_G3xy_aa = I_ESP_Kx4y2z_F3x_aa+ABY*I_ESP_Ix3y2z_F3x_aa;
    Double I_ESP_Ix2y3z_G3xy_aa = I_ESP_Kx3y3z_F3x_aa+ABY*I_ESP_Ix2y3z_F3x_aa;
    Double I_ESP_Ixy4z_G3xy_aa = I_ESP_Kx2y4z_F3x_aa+ABY*I_ESP_Ixy4z_F3x_aa;
    Double I_ESP_Ix5z_G3xy_aa = I_ESP_Kxy5z_F3x_aa+ABY*I_ESP_Ix5z_F3x_aa;
    Double I_ESP_I6y_G3xy_aa = I_ESP_K7y_F3x_aa+ABY*I_ESP_I6y_F3x_aa;
    Double I_ESP_I5yz_G3xy_aa = I_ESP_K6yz_F3x_aa+ABY*I_ESP_I5yz_F3x_aa;
    Double I_ESP_I4y2z_G3xy_aa = I_ESP_K5y2z_F3x_aa+ABY*I_ESP_I4y2z_F3x_aa;
    Double I_ESP_I3y3z_G3xy_aa = I_ESP_K4y3z_F3x_aa+ABY*I_ESP_I3y3z_F3x_aa;
    Double I_ESP_I2y4z_G3xy_aa = I_ESP_K3y4z_F3x_aa+ABY*I_ESP_I2y4z_F3x_aa;
    Double I_ESP_Iy5z_G3xy_aa = I_ESP_K2y5z_F3x_aa+ABY*I_ESP_Iy5z_F3x_aa;
    Double I_ESP_I6z_G3xy_aa = I_ESP_Ky6z_F3x_aa+ABY*I_ESP_I6z_F3x_aa;
    Double I_ESP_I6x_G3xz_aa = I_ESP_K6xz_F3x_aa+ABZ*I_ESP_I6x_F3x_aa;
    Double I_ESP_I5xy_G3xz_aa = I_ESP_K5xyz_F3x_aa+ABZ*I_ESP_I5xy_F3x_aa;
    Double I_ESP_I5xz_G3xz_aa = I_ESP_K5x2z_F3x_aa+ABZ*I_ESP_I5xz_F3x_aa;
    Double I_ESP_I4x2y_G3xz_aa = I_ESP_K4x2yz_F3x_aa+ABZ*I_ESP_I4x2y_F3x_aa;
    Double I_ESP_I4xyz_G3xz_aa = I_ESP_K4xy2z_F3x_aa+ABZ*I_ESP_I4xyz_F3x_aa;
    Double I_ESP_I4x2z_G3xz_aa = I_ESP_K4x3z_F3x_aa+ABZ*I_ESP_I4x2z_F3x_aa;
    Double I_ESP_I3x3y_G3xz_aa = I_ESP_K3x3yz_F3x_aa+ABZ*I_ESP_I3x3y_F3x_aa;
    Double I_ESP_I3x2yz_G3xz_aa = I_ESP_K3x2y2z_F3x_aa+ABZ*I_ESP_I3x2yz_F3x_aa;
    Double I_ESP_I3xy2z_G3xz_aa = I_ESP_K3xy3z_F3x_aa+ABZ*I_ESP_I3xy2z_F3x_aa;
    Double I_ESP_I3x3z_G3xz_aa = I_ESP_K3x4z_F3x_aa+ABZ*I_ESP_I3x3z_F3x_aa;
    Double I_ESP_I2x4y_G3xz_aa = I_ESP_K2x4yz_F3x_aa+ABZ*I_ESP_I2x4y_F3x_aa;
    Double I_ESP_I2x3yz_G3xz_aa = I_ESP_K2x3y2z_F3x_aa+ABZ*I_ESP_I2x3yz_F3x_aa;
    Double I_ESP_I2x2y2z_G3xz_aa = I_ESP_K2x2y3z_F3x_aa+ABZ*I_ESP_I2x2y2z_F3x_aa;
    Double I_ESP_I2xy3z_G3xz_aa = I_ESP_K2xy4z_F3x_aa+ABZ*I_ESP_I2xy3z_F3x_aa;
    Double I_ESP_I2x4z_G3xz_aa = I_ESP_K2x5z_F3x_aa+ABZ*I_ESP_I2x4z_F3x_aa;
    Double I_ESP_Ix5y_G3xz_aa = I_ESP_Kx5yz_F3x_aa+ABZ*I_ESP_Ix5y_F3x_aa;
    Double I_ESP_Ix4yz_G3xz_aa = I_ESP_Kx4y2z_F3x_aa+ABZ*I_ESP_Ix4yz_F3x_aa;
    Double I_ESP_Ix3y2z_G3xz_aa = I_ESP_Kx3y3z_F3x_aa+ABZ*I_ESP_Ix3y2z_F3x_aa;
    Double I_ESP_Ix2y3z_G3xz_aa = I_ESP_Kx2y4z_F3x_aa+ABZ*I_ESP_Ix2y3z_F3x_aa;
    Double I_ESP_Ixy4z_G3xz_aa = I_ESP_Kxy5z_F3x_aa+ABZ*I_ESP_Ixy4z_F3x_aa;
    Double I_ESP_Ix5z_G3xz_aa = I_ESP_Kx6z_F3x_aa+ABZ*I_ESP_Ix5z_F3x_aa;
    Double I_ESP_I6y_G3xz_aa = I_ESP_K6yz_F3x_aa+ABZ*I_ESP_I6y_F3x_aa;
    Double I_ESP_I5yz_G3xz_aa = I_ESP_K5y2z_F3x_aa+ABZ*I_ESP_I5yz_F3x_aa;
    Double I_ESP_I4y2z_G3xz_aa = I_ESP_K4y3z_F3x_aa+ABZ*I_ESP_I4y2z_F3x_aa;
    Double I_ESP_I3y3z_G3xz_aa = I_ESP_K3y4z_F3x_aa+ABZ*I_ESP_I3y3z_F3x_aa;
    Double I_ESP_I2y4z_G3xz_aa = I_ESP_K2y5z_F3x_aa+ABZ*I_ESP_I2y4z_F3x_aa;
    Double I_ESP_Iy5z_G3xz_aa = I_ESP_Ky6z_F3x_aa+ABZ*I_ESP_Iy5z_F3x_aa;
    Double I_ESP_I6z_G3xz_aa = I_ESP_K7z_F3x_aa+ABZ*I_ESP_I6z_F3x_aa;
    Double I_ESP_I6x_G2x2y_aa = I_ESP_K6xy_F2xy_aa+ABY*I_ESP_I6x_F2xy_aa;
    Double I_ESP_I5xy_G2x2y_aa = I_ESP_K5x2y_F2xy_aa+ABY*I_ESP_I5xy_F2xy_aa;
    Double I_ESP_I5xz_G2x2y_aa = I_ESP_K5xyz_F2xy_aa+ABY*I_ESP_I5xz_F2xy_aa;
    Double I_ESP_I4x2y_G2x2y_aa = I_ESP_K4x3y_F2xy_aa+ABY*I_ESP_I4x2y_F2xy_aa;
    Double I_ESP_I4xyz_G2x2y_aa = I_ESP_K4x2yz_F2xy_aa+ABY*I_ESP_I4xyz_F2xy_aa;
    Double I_ESP_I4x2z_G2x2y_aa = I_ESP_K4xy2z_F2xy_aa+ABY*I_ESP_I4x2z_F2xy_aa;
    Double I_ESP_I3x3y_G2x2y_aa = I_ESP_K3x4y_F2xy_aa+ABY*I_ESP_I3x3y_F2xy_aa;
    Double I_ESP_I3x2yz_G2x2y_aa = I_ESP_K3x3yz_F2xy_aa+ABY*I_ESP_I3x2yz_F2xy_aa;
    Double I_ESP_I3xy2z_G2x2y_aa = I_ESP_K3x2y2z_F2xy_aa+ABY*I_ESP_I3xy2z_F2xy_aa;
    Double I_ESP_I3x3z_G2x2y_aa = I_ESP_K3xy3z_F2xy_aa+ABY*I_ESP_I3x3z_F2xy_aa;
    Double I_ESP_I2x4y_G2x2y_aa = I_ESP_K2x5y_F2xy_aa+ABY*I_ESP_I2x4y_F2xy_aa;
    Double I_ESP_I2x3yz_G2x2y_aa = I_ESP_K2x4yz_F2xy_aa+ABY*I_ESP_I2x3yz_F2xy_aa;
    Double I_ESP_I2x2y2z_G2x2y_aa = I_ESP_K2x3y2z_F2xy_aa+ABY*I_ESP_I2x2y2z_F2xy_aa;
    Double I_ESP_I2xy3z_G2x2y_aa = I_ESP_K2x2y3z_F2xy_aa+ABY*I_ESP_I2xy3z_F2xy_aa;
    Double I_ESP_I2x4z_G2x2y_aa = I_ESP_K2xy4z_F2xy_aa+ABY*I_ESP_I2x4z_F2xy_aa;
    Double I_ESP_Ix5y_G2x2y_aa = I_ESP_Kx6y_F2xy_aa+ABY*I_ESP_Ix5y_F2xy_aa;
    Double I_ESP_Ix4yz_G2x2y_aa = I_ESP_Kx5yz_F2xy_aa+ABY*I_ESP_Ix4yz_F2xy_aa;
    Double I_ESP_Ix3y2z_G2x2y_aa = I_ESP_Kx4y2z_F2xy_aa+ABY*I_ESP_Ix3y2z_F2xy_aa;
    Double I_ESP_Ix2y3z_G2x2y_aa = I_ESP_Kx3y3z_F2xy_aa+ABY*I_ESP_Ix2y3z_F2xy_aa;
    Double I_ESP_Ixy4z_G2x2y_aa = I_ESP_Kx2y4z_F2xy_aa+ABY*I_ESP_Ixy4z_F2xy_aa;
    Double I_ESP_Ix5z_G2x2y_aa = I_ESP_Kxy5z_F2xy_aa+ABY*I_ESP_Ix5z_F2xy_aa;
    Double I_ESP_I6y_G2x2y_aa = I_ESP_K7y_F2xy_aa+ABY*I_ESP_I6y_F2xy_aa;
    Double I_ESP_I5yz_G2x2y_aa = I_ESP_K6yz_F2xy_aa+ABY*I_ESP_I5yz_F2xy_aa;
    Double I_ESP_I4y2z_G2x2y_aa = I_ESP_K5y2z_F2xy_aa+ABY*I_ESP_I4y2z_F2xy_aa;
    Double I_ESP_I3y3z_G2x2y_aa = I_ESP_K4y3z_F2xy_aa+ABY*I_ESP_I3y3z_F2xy_aa;
    Double I_ESP_I2y4z_G2x2y_aa = I_ESP_K3y4z_F2xy_aa+ABY*I_ESP_I2y4z_F2xy_aa;
    Double I_ESP_Iy5z_G2x2y_aa = I_ESP_K2y5z_F2xy_aa+ABY*I_ESP_Iy5z_F2xy_aa;
    Double I_ESP_I6z_G2x2y_aa = I_ESP_Ky6z_F2xy_aa+ABY*I_ESP_I6z_F2xy_aa;
    Double I_ESP_I6x_G2xyz_aa = I_ESP_K6xz_F2xy_aa+ABZ*I_ESP_I6x_F2xy_aa;
    Double I_ESP_I5xy_G2xyz_aa = I_ESP_K5xyz_F2xy_aa+ABZ*I_ESP_I5xy_F2xy_aa;
    Double I_ESP_I5xz_G2xyz_aa = I_ESP_K5x2z_F2xy_aa+ABZ*I_ESP_I5xz_F2xy_aa;
    Double I_ESP_I4x2y_G2xyz_aa = I_ESP_K4x2yz_F2xy_aa+ABZ*I_ESP_I4x2y_F2xy_aa;
    Double I_ESP_I4xyz_G2xyz_aa = I_ESP_K4xy2z_F2xy_aa+ABZ*I_ESP_I4xyz_F2xy_aa;
    Double I_ESP_I4x2z_G2xyz_aa = I_ESP_K4x3z_F2xy_aa+ABZ*I_ESP_I4x2z_F2xy_aa;
    Double I_ESP_I3x3y_G2xyz_aa = I_ESP_K3x3yz_F2xy_aa+ABZ*I_ESP_I3x3y_F2xy_aa;
    Double I_ESP_I3x2yz_G2xyz_aa = I_ESP_K3x2y2z_F2xy_aa+ABZ*I_ESP_I3x2yz_F2xy_aa;
    Double I_ESP_I3xy2z_G2xyz_aa = I_ESP_K3xy3z_F2xy_aa+ABZ*I_ESP_I3xy2z_F2xy_aa;
    Double I_ESP_I3x3z_G2xyz_aa = I_ESP_K3x4z_F2xy_aa+ABZ*I_ESP_I3x3z_F2xy_aa;
    Double I_ESP_I2x4y_G2xyz_aa = I_ESP_K2x4yz_F2xy_aa+ABZ*I_ESP_I2x4y_F2xy_aa;
    Double I_ESP_I2x3yz_G2xyz_aa = I_ESP_K2x3y2z_F2xy_aa+ABZ*I_ESP_I2x3yz_F2xy_aa;
    Double I_ESP_I2x2y2z_G2xyz_aa = I_ESP_K2x2y3z_F2xy_aa+ABZ*I_ESP_I2x2y2z_F2xy_aa;
    Double I_ESP_I2xy3z_G2xyz_aa = I_ESP_K2xy4z_F2xy_aa+ABZ*I_ESP_I2xy3z_F2xy_aa;
    Double I_ESP_I2x4z_G2xyz_aa = I_ESP_K2x5z_F2xy_aa+ABZ*I_ESP_I2x4z_F2xy_aa;
    Double I_ESP_Ix5y_G2xyz_aa = I_ESP_Kx5yz_F2xy_aa+ABZ*I_ESP_Ix5y_F2xy_aa;
    Double I_ESP_Ix4yz_G2xyz_aa = I_ESP_Kx4y2z_F2xy_aa+ABZ*I_ESP_Ix4yz_F2xy_aa;
    Double I_ESP_Ix3y2z_G2xyz_aa = I_ESP_Kx3y3z_F2xy_aa+ABZ*I_ESP_Ix3y2z_F2xy_aa;
    Double I_ESP_Ix2y3z_G2xyz_aa = I_ESP_Kx2y4z_F2xy_aa+ABZ*I_ESP_Ix2y3z_F2xy_aa;
    Double I_ESP_Ixy4z_G2xyz_aa = I_ESP_Kxy5z_F2xy_aa+ABZ*I_ESP_Ixy4z_F2xy_aa;
    Double I_ESP_Ix5z_G2xyz_aa = I_ESP_Kx6z_F2xy_aa+ABZ*I_ESP_Ix5z_F2xy_aa;
    Double I_ESP_I6y_G2xyz_aa = I_ESP_K6yz_F2xy_aa+ABZ*I_ESP_I6y_F2xy_aa;
    Double I_ESP_I5yz_G2xyz_aa = I_ESP_K5y2z_F2xy_aa+ABZ*I_ESP_I5yz_F2xy_aa;
    Double I_ESP_I4y2z_G2xyz_aa = I_ESP_K4y3z_F2xy_aa+ABZ*I_ESP_I4y2z_F2xy_aa;
    Double I_ESP_I3y3z_G2xyz_aa = I_ESP_K3y4z_F2xy_aa+ABZ*I_ESP_I3y3z_F2xy_aa;
    Double I_ESP_I2y4z_G2xyz_aa = I_ESP_K2y5z_F2xy_aa+ABZ*I_ESP_I2y4z_F2xy_aa;
    Double I_ESP_Iy5z_G2xyz_aa = I_ESP_Ky6z_F2xy_aa+ABZ*I_ESP_Iy5z_F2xy_aa;
    Double I_ESP_I6z_G2xyz_aa = I_ESP_K7z_F2xy_aa+ABZ*I_ESP_I6z_F2xy_aa;
    Double I_ESP_I6x_G2x2z_aa = I_ESP_K6xz_F2xz_aa+ABZ*I_ESP_I6x_F2xz_aa;
    Double I_ESP_I5xy_G2x2z_aa = I_ESP_K5xyz_F2xz_aa+ABZ*I_ESP_I5xy_F2xz_aa;
    Double I_ESP_I5xz_G2x2z_aa = I_ESP_K5x2z_F2xz_aa+ABZ*I_ESP_I5xz_F2xz_aa;
    Double I_ESP_I4x2y_G2x2z_aa = I_ESP_K4x2yz_F2xz_aa+ABZ*I_ESP_I4x2y_F2xz_aa;
    Double I_ESP_I4xyz_G2x2z_aa = I_ESP_K4xy2z_F2xz_aa+ABZ*I_ESP_I4xyz_F2xz_aa;
    Double I_ESP_I4x2z_G2x2z_aa = I_ESP_K4x3z_F2xz_aa+ABZ*I_ESP_I4x2z_F2xz_aa;
    Double I_ESP_I3x3y_G2x2z_aa = I_ESP_K3x3yz_F2xz_aa+ABZ*I_ESP_I3x3y_F2xz_aa;
    Double I_ESP_I3x2yz_G2x2z_aa = I_ESP_K3x2y2z_F2xz_aa+ABZ*I_ESP_I3x2yz_F2xz_aa;
    Double I_ESP_I3xy2z_G2x2z_aa = I_ESP_K3xy3z_F2xz_aa+ABZ*I_ESP_I3xy2z_F2xz_aa;
    Double I_ESP_I3x3z_G2x2z_aa = I_ESP_K3x4z_F2xz_aa+ABZ*I_ESP_I3x3z_F2xz_aa;
    Double I_ESP_I2x4y_G2x2z_aa = I_ESP_K2x4yz_F2xz_aa+ABZ*I_ESP_I2x4y_F2xz_aa;
    Double I_ESP_I2x3yz_G2x2z_aa = I_ESP_K2x3y2z_F2xz_aa+ABZ*I_ESP_I2x3yz_F2xz_aa;
    Double I_ESP_I2x2y2z_G2x2z_aa = I_ESP_K2x2y3z_F2xz_aa+ABZ*I_ESP_I2x2y2z_F2xz_aa;
    Double I_ESP_I2xy3z_G2x2z_aa = I_ESP_K2xy4z_F2xz_aa+ABZ*I_ESP_I2xy3z_F2xz_aa;
    Double I_ESP_I2x4z_G2x2z_aa = I_ESP_K2x5z_F2xz_aa+ABZ*I_ESP_I2x4z_F2xz_aa;
    Double I_ESP_Ix5y_G2x2z_aa = I_ESP_Kx5yz_F2xz_aa+ABZ*I_ESP_Ix5y_F2xz_aa;
    Double I_ESP_Ix4yz_G2x2z_aa = I_ESP_Kx4y2z_F2xz_aa+ABZ*I_ESP_Ix4yz_F2xz_aa;
    Double I_ESP_Ix3y2z_G2x2z_aa = I_ESP_Kx3y3z_F2xz_aa+ABZ*I_ESP_Ix3y2z_F2xz_aa;
    Double I_ESP_Ix2y3z_G2x2z_aa = I_ESP_Kx2y4z_F2xz_aa+ABZ*I_ESP_Ix2y3z_F2xz_aa;
    Double I_ESP_Ixy4z_G2x2z_aa = I_ESP_Kxy5z_F2xz_aa+ABZ*I_ESP_Ixy4z_F2xz_aa;
    Double I_ESP_Ix5z_G2x2z_aa = I_ESP_Kx6z_F2xz_aa+ABZ*I_ESP_Ix5z_F2xz_aa;
    Double I_ESP_I6y_G2x2z_aa = I_ESP_K6yz_F2xz_aa+ABZ*I_ESP_I6y_F2xz_aa;
    Double I_ESP_I5yz_G2x2z_aa = I_ESP_K5y2z_F2xz_aa+ABZ*I_ESP_I5yz_F2xz_aa;
    Double I_ESP_I4y2z_G2x2z_aa = I_ESP_K4y3z_F2xz_aa+ABZ*I_ESP_I4y2z_F2xz_aa;
    Double I_ESP_I3y3z_G2x2z_aa = I_ESP_K3y4z_F2xz_aa+ABZ*I_ESP_I3y3z_F2xz_aa;
    Double I_ESP_I2y4z_G2x2z_aa = I_ESP_K2y5z_F2xz_aa+ABZ*I_ESP_I2y4z_F2xz_aa;
    Double I_ESP_Iy5z_G2x2z_aa = I_ESP_Ky6z_F2xz_aa+ABZ*I_ESP_Iy5z_F2xz_aa;
    Double I_ESP_I6z_G2x2z_aa = I_ESP_K7z_F2xz_aa+ABZ*I_ESP_I6z_F2xz_aa;
    Double I_ESP_I6x_Gx3y_aa = I_ESP_K7x_F3y_aa+ABX*I_ESP_I6x_F3y_aa;
    Double I_ESP_I5xy_Gx3y_aa = I_ESP_K6xy_F3y_aa+ABX*I_ESP_I5xy_F3y_aa;
    Double I_ESP_I5xz_Gx3y_aa = I_ESP_K6xz_F3y_aa+ABX*I_ESP_I5xz_F3y_aa;
    Double I_ESP_I4x2y_Gx3y_aa = I_ESP_K5x2y_F3y_aa+ABX*I_ESP_I4x2y_F3y_aa;
    Double I_ESP_I4xyz_Gx3y_aa = I_ESP_K5xyz_F3y_aa+ABX*I_ESP_I4xyz_F3y_aa;
    Double I_ESP_I4x2z_Gx3y_aa = I_ESP_K5x2z_F3y_aa+ABX*I_ESP_I4x2z_F3y_aa;
    Double I_ESP_I3x3y_Gx3y_aa = I_ESP_K4x3y_F3y_aa+ABX*I_ESP_I3x3y_F3y_aa;
    Double I_ESP_I3x2yz_Gx3y_aa = I_ESP_K4x2yz_F3y_aa+ABX*I_ESP_I3x2yz_F3y_aa;
    Double I_ESP_I3xy2z_Gx3y_aa = I_ESP_K4xy2z_F3y_aa+ABX*I_ESP_I3xy2z_F3y_aa;
    Double I_ESP_I3x3z_Gx3y_aa = I_ESP_K4x3z_F3y_aa+ABX*I_ESP_I3x3z_F3y_aa;
    Double I_ESP_I2x4y_Gx3y_aa = I_ESP_K3x4y_F3y_aa+ABX*I_ESP_I2x4y_F3y_aa;
    Double I_ESP_I2x3yz_Gx3y_aa = I_ESP_K3x3yz_F3y_aa+ABX*I_ESP_I2x3yz_F3y_aa;
    Double I_ESP_I2x2y2z_Gx3y_aa = I_ESP_K3x2y2z_F3y_aa+ABX*I_ESP_I2x2y2z_F3y_aa;
    Double I_ESP_I2xy3z_Gx3y_aa = I_ESP_K3xy3z_F3y_aa+ABX*I_ESP_I2xy3z_F3y_aa;
    Double I_ESP_I2x4z_Gx3y_aa = I_ESP_K3x4z_F3y_aa+ABX*I_ESP_I2x4z_F3y_aa;
    Double I_ESP_Ix5y_Gx3y_aa = I_ESP_K2x5y_F3y_aa+ABX*I_ESP_Ix5y_F3y_aa;
    Double I_ESP_Ix4yz_Gx3y_aa = I_ESP_K2x4yz_F3y_aa+ABX*I_ESP_Ix4yz_F3y_aa;
    Double I_ESP_Ix3y2z_Gx3y_aa = I_ESP_K2x3y2z_F3y_aa+ABX*I_ESP_Ix3y2z_F3y_aa;
    Double I_ESP_Ix2y3z_Gx3y_aa = I_ESP_K2x2y3z_F3y_aa+ABX*I_ESP_Ix2y3z_F3y_aa;
    Double I_ESP_Ixy4z_Gx3y_aa = I_ESP_K2xy4z_F3y_aa+ABX*I_ESP_Ixy4z_F3y_aa;
    Double I_ESP_Ix5z_Gx3y_aa = I_ESP_K2x5z_F3y_aa+ABX*I_ESP_Ix5z_F3y_aa;
    Double I_ESP_I6y_Gx3y_aa = I_ESP_Kx6y_F3y_aa+ABX*I_ESP_I6y_F3y_aa;
    Double I_ESP_I5yz_Gx3y_aa = I_ESP_Kx5yz_F3y_aa+ABX*I_ESP_I5yz_F3y_aa;
    Double I_ESP_I4y2z_Gx3y_aa = I_ESP_Kx4y2z_F3y_aa+ABX*I_ESP_I4y2z_F3y_aa;
    Double I_ESP_I3y3z_Gx3y_aa = I_ESP_Kx3y3z_F3y_aa+ABX*I_ESP_I3y3z_F3y_aa;
    Double I_ESP_I2y4z_Gx3y_aa = I_ESP_Kx2y4z_F3y_aa+ABX*I_ESP_I2y4z_F3y_aa;
    Double I_ESP_Iy5z_Gx3y_aa = I_ESP_Kxy5z_F3y_aa+ABX*I_ESP_Iy5z_F3y_aa;
    Double I_ESP_I6z_Gx3y_aa = I_ESP_Kx6z_F3y_aa+ABX*I_ESP_I6z_F3y_aa;
    Double I_ESP_I6x_Gx2yz_aa = I_ESP_K6xz_Fx2y_aa+ABZ*I_ESP_I6x_Fx2y_aa;
    Double I_ESP_I5xy_Gx2yz_aa = I_ESP_K5xyz_Fx2y_aa+ABZ*I_ESP_I5xy_Fx2y_aa;
    Double I_ESP_I5xz_Gx2yz_aa = I_ESP_K5x2z_Fx2y_aa+ABZ*I_ESP_I5xz_Fx2y_aa;
    Double I_ESP_I4x2y_Gx2yz_aa = I_ESP_K4x2yz_Fx2y_aa+ABZ*I_ESP_I4x2y_Fx2y_aa;
    Double I_ESP_I4xyz_Gx2yz_aa = I_ESP_K4xy2z_Fx2y_aa+ABZ*I_ESP_I4xyz_Fx2y_aa;
    Double I_ESP_I4x2z_Gx2yz_aa = I_ESP_K4x3z_Fx2y_aa+ABZ*I_ESP_I4x2z_Fx2y_aa;
    Double I_ESP_I3x3y_Gx2yz_aa = I_ESP_K3x3yz_Fx2y_aa+ABZ*I_ESP_I3x3y_Fx2y_aa;
    Double I_ESP_I3x2yz_Gx2yz_aa = I_ESP_K3x2y2z_Fx2y_aa+ABZ*I_ESP_I3x2yz_Fx2y_aa;
    Double I_ESP_I3xy2z_Gx2yz_aa = I_ESP_K3xy3z_Fx2y_aa+ABZ*I_ESP_I3xy2z_Fx2y_aa;
    Double I_ESP_I3x3z_Gx2yz_aa = I_ESP_K3x4z_Fx2y_aa+ABZ*I_ESP_I3x3z_Fx2y_aa;
    Double I_ESP_I2x4y_Gx2yz_aa = I_ESP_K2x4yz_Fx2y_aa+ABZ*I_ESP_I2x4y_Fx2y_aa;
    Double I_ESP_I2x3yz_Gx2yz_aa = I_ESP_K2x3y2z_Fx2y_aa+ABZ*I_ESP_I2x3yz_Fx2y_aa;
    Double I_ESP_I2x2y2z_Gx2yz_aa = I_ESP_K2x2y3z_Fx2y_aa+ABZ*I_ESP_I2x2y2z_Fx2y_aa;
    Double I_ESP_I2xy3z_Gx2yz_aa = I_ESP_K2xy4z_Fx2y_aa+ABZ*I_ESP_I2xy3z_Fx2y_aa;
    Double I_ESP_I2x4z_Gx2yz_aa = I_ESP_K2x5z_Fx2y_aa+ABZ*I_ESP_I2x4z_Fx2y_aa;
    Double I_ESP_Ix5y_Gx2yz_aa = I_ESP_Kx5yz_Fx2y_aa+ABZ*I_ESP_Ix5y_Fx2y_aa;
    Double I_ESP_Ix4yz_Gx2yz_aa = I_ESP_Kx4y2z_Fx2y_aa+ABZ*I_ESP_Ix4yz_Fx2y_aa;
    Double I_ESP_Ix3y2z_Gx2yz_aa = I_ESP_Kx3y3z_Fx2y_aa+ABZ*I_ESP_Ix3y2z_Fx2y_aa;
    Double I_ESP_Ix2y3z_Gx2yz_aa = I_ESP_Kx2y4z_Fx2y_aa+ABZ*I_ESP_Ix2y3z_Fx2y_aa;
    Double I_ESP_Ixy4z_Gx2yz_aa = I_ESP_Kxy5z_Fx2y_aa+ABZ*I_ESP_Ixy4z_Fx2y_aa;
    Double I_ESP_Ix5z_Gx2yz_aa = I_ESP_Kx6z_Fx2y_aa+ABZ*I_ESP_Ix5z_Fx2y_aa;
    Double I_ESP_I6y_Gx2yz_aa = I_ESP_K6yz_Fx2y_aa+ABZ*I_ESP_I6y_Fx2y_aa;
    Double I_ESP_I5yz_Gx2yz_aa = I_ESP_K5y2z_Fx2y_aa+ABZ*I_ESP_I5yz_Fx2y_aa;
    Double I_ESP_I4y2z_Gx2yz_aa = I_ESP_K4y3z_Fx2y_aa+ABZ*I_ESP_I4y2z_Fx2y_aa;
    Double I_ESP_I3y3z_Gx2yz_aa = I_ESP_K3y4z_Fx2y_aa+ABZ*I_ESP_I3y3z_Fx2y_aa;
    Double I_ESP_I2y4z_Gx2yz_aa = I_ESP_K2y5z_Fx2y_aa+ABZ*I_ESP_I2y4z_Fx2y_aa;
    Double I_ESP_Iy5z_Gx2yz_aa = I_ESP_Ky6z_Fx2y_aa+ABZ*I_ESP_Iy5z_Fx2y_aa;
    Double I_ESP_I6z_Gx2yz_aa = I_ESP_K7z_Fx2y_aa+ABZ*I_ESP_I6z_Fx2y_aa;
    Double I_ESP_I6x_Gxy2z_aa = I_ESP_K6xy_Fx2z_aa+ABY*I_ESP_I6x_Fx2z_aa;
    Double I_ESP_I5xy_Gxy2z_aa = I_ESP_K5x2y_Fx2z_aa+ABY*I_ESP_I5xy_Fx2z_aa;
    Double I_ESP_I5xz_Gxy2z_aa = I_ESP_K5xyz_Fx2z_aa+ABY*I_ESP_I5xz_Fx2z_aa;
    Double I_ESP_I4x2y_Gxy2z_aa = I_ESP_K4x3y_Fx2z_aa+ABY*I_ESP_I4x2y_Fx2z_aa;
    Double I_ESP_I4xyz_Gxy2z_aa = I_ESP_K4x2yz_Fx2z_aa+ABY*I_ESP_I4xyz_Fx2z_aa;
    Double I_ESP_I4x2z_Gxy2z_aa = I_ESP_K4xy2z_Fx2z_aa+ABY*I_ESP_I4x2z_Fx2z_aa;
    Double I_ESP_I3x3y_Gxy2z_aa = I_ESP_K3x4y_Fx2z_aa+ABY*I_ESP_I3x3y_Fx2z_aa;
    Double I_ESP_I3x2yz_Gxy2z_aa = I_ESP_K3x3yz_Fx2z_aa+ABY*I_ESP_I3x2yz_Fx2z_aa;
    Double I_ESP_I3xy2z_Gxy2z_aa = I_ESP_K3x2y2z_Fx2z_aa+ABY*I_ESP_I3xy2z_Fx2z_aa;
    Double I_ESP_I3x3z_Gxy2z_aa = I_ESP_K3xy3z_Fx2z_aa+ABY*I_ESP_I3x3z_Fx2z_aa;
    Double I_ESP_I2x4y_Gxy2z_aa = I_ESP_K2x5y_Fx2z_aa+ABY*I_ESP_I2x4y_Fx2z_aa;
    Double I_ESP_I2x3yz_Gxy2z_aa = I_ESP_K2x4yz_Fx2z_aa+ABY*I_ESP_I2x3yz_Fx2z_aa;
    Double I_ESP_I2x2y2z_Gxy2z_aa = I_ESP_K2x3y2z_Fx2z_aa+ABY*I_ESP_I2x2y2z_Fx2z_aa;
    Double I_ESP_I2xy3z_Gxy2z_aa = I_ESP_K2x2y3z_Fx2z_aa+ABY*I_ESP_I2xy3z_Fx2z_aa;
    Double I_ESP_I2x4z_Gxy2z_aa = I_ESP_K2xy4z_Fx2z_aa+ABY*I_ESP_I2x4z_Fx2z_aa;
    Double I_ESP_Ix5y_Gxy2z_aa = I_ESP_Kx6y_Fx2z_aa+ABY*I_ESP_Ix5y_Fx2z_aa;
    Double I_ESP_Ix4yz_Gxy2z_aa = I_ESP_Kx5yz_Fx2z_aa+ABY*I_ESP_Ix4yz_Fx2z_aa;
    Double I_ESP_Ix3y2z_Gxy2z_aa = I_ESP_Kx4y2z_Fx2z_aa+ABY*I_ESP_Ix3y2z_Fx2z_aa;
    Double I_ESP_Ix2y3z_Gxy2z_aa = I_ESP_Kx3y3z_Fx2z_aa+ABY*I_ESP_Ix2y3z_Fx2z_aa;
    Double I_ESP_Ixy4z_Gxy2z_aa = I_ESP_Kx2y4z_Fx2z_aa+ABY*I_ESP_Ixy4z_Fx2z_aa;
    Double I_ESP_Ix5z_Gxy2z_aa = I_ESP_Kxy5z_Fx2z_aa+ABY*I_ESP_Ix5z_Fx2z_aa;
    Double I_ESP_I6y_Gxy2z_aa = I_ESP_K7y_Fx2z_aa+ABY*I_ESP_I6y_Fx2z_aa;
    Double I_ESP_I5yz_Gxy2z_aa = I_ESP_K6yz_Fx2z_aa+ABY*I_ESP_I5yz_Fx2z_aa;
    Double I_ESP_I4y2z_Gxy2z_aa = I_ESP_K5y2z_Fx2z_aa+ABY*I_ESP_I4y2z_Fx2z_aa;
    Double I_ESP_I3y3z_Gxy2z_aa = I_ESP_K4y3z_Fx2z_aa+ABY*I_ESP_I3y3z_Fx2z_aa;
    Double I_ESP_I2y4z_Gxy2z_aa = I_ESP_K3y4z_Fx2z_aa+ABY*I_ESP_I2y4z_Fx2z_aa;
    Double I_ESP_Iy5z_Gxy2z_aa = I_ESP_K2y5z_Fx2z_aa+ABY*I_ESP_Iy5z_Fx2z_aa;
    Double I_ESP_I6z_Gxy2z_aa = I_ESP_Ky6z_Fx2z_aa+ABY*I_ESP_I6z_Fx2z_aa;
    Double I_ESP_I6x_Gx3z_aa = I_ESP_K7x_F3z_aa+ABX*I_ESP_I6x_F3z_aa;
    Double I_ESP_I5xy_Gx3z_aa = I_ESP_K6xy_F3z_aa+ABX*I_ESP_I5xy_F3z_aa;
    Double I_ESP_I5xz_Gx3z_aa = I_ESP_K6xz_F3z_aa+ABX*I_ESP_I5xz_F3z_aa;
    Double I_ESP_I4x2y_Gx3z_aa = I_ESP_K5x2y_F3z_aa+ABX*I_ESP_I4x2y_F3z_aa;
    Double I_ESP_I4xyz_Gx3z_aa = I_ESP_K5xyz_F3z_aa+ABX*I_ESP_I4xyz_F3z_aa;
    Double I_ESP_I4x2z_Gx3z_aa = I_ESP_K5x2z_F3z_aa+ABX*I_ESP_I4x2z_F3z_aa;
    Double I_ESP_I3x3y_Gx3z_aa = I_ESP_K4x3y_F3z_aa+ABX*I_ESP_I3x3y_F3z_aa;
    Double I_ESP_I3x2yz_Gx3z_aa = I_ESP_K4x2yz_F3z_aa+ABX*I_ESP_I3x2yz_F3z_aa;
    Double I_ESP_I3xy2z_Gx3z_aa = I_ESP_K4xy2z_F3z_aa+ABX*I_ESP_I3xy2z_F3z_aa;
    Double I_ESP_I3x3z_Gx3z_aa = I_ESP_K4x3z_F3z_aa+ABX*I_ESP_I3x3z_F3z_aa;
    Double I_ESP_I2x4y_Gx3z_aa = I_ESP_K3x4y_F3z_aa+ABX*I_ESP_I2x4y_F3z_aa;
    Double I_ESP_I2x3yz_Gx3z_aa = I_ESP_K3x3yz_F3z_aa+ABX*I_ESP_I2x3yz_F3z_aa;
    Double I_ESP_I2x2y2z_Gx3z_aa = I_ESP_K3x2y2z_F3z_aa+ABX*I_ESP_I2x2y2z_F3z_aa;
    Double I_ESP_I2xy3z_Gx3z_aa = I_ESP_K3xy3z_F3z_aa+ABX*I_ESP_I2xy3z_F3z_aa;
    Double I_ESP_I2x4z_Gx3z_aa = I_ESP_K3x4z_F3z_aa+ABX*I_ESP_I2x4z_F3z_aa;
    Double I_ESP_Ix5y_Gx3z_aa = I_ESP_K2x5y_F3z_aa+ABX*I_ESP_Ix5y_F3z_aa;
    Double I_ESP_Ix4yz_Gx3z_aa = I_ESP_K2x4yz_F3z_aa+ABX*I_ESP_Ix4yz_F3z_aa;
    Double I_ESP_Ix3y2z_Gx3z_aa = I_ESP_K2x3y2z_F3z_aa+ABX*I_ESP_Ix3y2z_F3z_aa;
    Double I_ESP_Ix2y3z_Gx3z_aa = I_ESP_K2x2y3z_F3z_aa+ABX*I_ESP_Ix2y3z_F3z_aa;
    Double I_ESP_Ixy4z_Gx3z_aa = I_ESP_K2xy4z_F3z_aa+ABX*I_ESP_Ixy4z_F3z_aa;
    Double I_ESP_Ix5z_Gx3z_aa = I_ESP_K2x5z_F3z_aa+ABX*I_ESP_Ix5z_F3z_aa;
    Double I_ESP_I6y_Gx3z_aa = I_ESP_Kx6y_F3z_aa+ABX*I_ESP_I6y_F3z_aa;
    Double I_ESP_I5yz_Gx3z_aa = I_ESP_Kx5yz_F3z_aa+ABX*I_ESP_I5yz_F3z_aa;
    Double I_ESP_I4y2z_Gx3z_aa = I_ESP_Kx4y2z_F3z_aa+ABX*I_ESP_I4y2z_F3z_aa;
    Double I_ESP_I3y3z_Gx3z_aa = I_ESP_Kx3y3z_F3z_aa+ABX*I_ESP_I3y3z_F3z_aa;
    Double I_ESP_I2y4z_Gx3z_aa = I_ESP_Kx2y4z_F3z_aa+ABX*I_ESP_I2y4z_F3z_aa;
    Double I_ESP_Iy5z_Gx3z_aa = I_ESP_Kxy5z_F3z_aa+ABX*I_ESP_Iy5z_F3z_aa;
    Double I_ESP_I6z_Gx3z_aa = I_ESP_Kx6z_F3z_aa+ABX*I_ESP_I6z_F3z_aa;
    Double I_ESP_I6x_G4y_aa = I_ESP_K6xy_F3y_aa+ABY*I_ESP_I6x_F3y_aa;
    Double I_ESP_I5xy_G4y_aa = I_ESP_K5x2y_F3y_aa+ABY*I_ESP_I5xy_F3y_aa;
    Double I_ESP_I5xz_G4y_aa = I_ESP_K5xyz_F3y_aa+ABY*I_ESP_I5xz_F3y_aa;
    Double I_ESP_I4x2y_G4y_aa = I_ESP_K4x3y_F3y_aa+ABY*I_ESP_I4x2y_F3y_aa;
    Double I_ESP_I4xyz_G4y_aa = I_ESP_K4x2yz_F3y_aa+ABY*I_ESP_I4xyz_F3y_aa;
    Double I_ESP_I4x2z_G4y_aa = I_ESP_K4xy2z_F3y_aa+ABY*I_ESP_I4x2z_F3y_aa;
    Double I_ESP_I3x3y_G4y_aa = I_ESP_K3x4y_F3y_aa+ABY*I_ESP_I3x3y_F3y_aa;
    Double I_ESP_I3x2yz_G4y_aa = I_ESP_K3x3yz_F3y_aa+ABY*I_ESP_I3x2yz_F3y_aa;
    Double I_ESP_I3xy2z_G4y_aa = I_ESP_K3x2y2z_F3y_aa+ABY*I_ESP_I3xy2z_F3y_aa;
    Double I_ESP_I3x3z_G4y_aa = I_ESP_K3xy3z_F3y_aa+ABY*I_ESP_I3x3z_F3y_aa;
    Double I_ESP_I2x4y_G4y_aa = I_ESP_K2x5y_F3y_aa+ABY*I_ESP_I2x4y_F3y_aa;
    Double I_ESP_I2x3yz_G4y_aa = I_ESP_K2x4yz_F3y_aa+ABY*I_ESP_I2x3yz_F3y_aa;
    Double I_ESP_I2x2y2z_G4y_aa = I_ESP_K2x3y2z_F3y_aa+ABY*I_ESP_I2x2y2z_F3y_aa;
    Double I_ESP_I2xy3z_G4y_aa = I_ESP_K2x2y3z_F3y_aa+ABY*I_ESP_I2xy3z_F3y_aa;
    Double I_ESP_I2x4z_G4y_aa = I_ESP_K2xy4z_F3y_aa+ABY*I_ESP_I2x4z_F3y_aa;
    Double I_ESP_Ix5y_G4y_aa = I_ESP_Kx6y_F3y_aa+ABY*I_ESP_Ix5y_F3y_aa;
    Double I_ESP_Ix4yz_G4y_aa = I_ESP_Kx5yz_F3y_aa+ABY*I_ESP_Ix4yz_F3y_aa;
    Double I_ESP_Ix3y2z_G4y_aa = I_ESP_Kx4y2z_F3y_aa+ABY*I_ESP_Ix3y2z_F3y_aa;
    Double I_ESP_Ix2y3z_G4y_aa = I_ESP_Kx3y3z_F3y_aa+ABY*I_ESP_Ix2y3z_F3y_aa;
    Double I_ESP_Ixy4z_G4y_aa = I_ESP_Kx2y4z_F3y_aa+ABY*I_ESP_Ixy4z_F3y_aa;
    Double I_ESP_Ix5z_G4y_aa = I_ESP_Kxy5z_F3y_aa+ABY*I_ESP_Ix5z_F3y_aa;
    Double I_ESP_I6y_G4y_aa = I_ESP_K7y_F3y_aa+ABY*I_ESP_I6y_F3y_aa;
    Double I_ESP_I5yz_G4y_aa = I_ESP_K6yz_F3y_aa+ABY*I_ESP_I5yz_F3y_aa;
    Double I_ESP_I4y2z_G4y_aa = I_ESP_K5y2z_F3y_aa+ABY*I_ESP_I4y2z_F3y_aa;
    Double I_ESP_I3y3z_G4y_aa = I_ESP_K4y3z_F3y_aa+ABY*I_ESP_I3y3z_F3y_aa;
    Double I_ESP_I2y4z_G4y_aa = I_ESP_K3y4z_F3y_aa+ABY*I_ESP_I2y4z_F3y_aa;
    Double I_ESP_Iy5z_G4y_aa = I_ESP_K2y5z_F3y_aa+ABY*I_ESP_Iy5z_F3y_aa;
    Double I_ESP_I6z_G4y_aa = I_ESP_Ky6z_F3y_aa+ABY*I_ESP_I6z_F3y_aa;
    Double I_ESP_I6x_G3yz_aa = I_ESP_K6xz_F3y_aa+ABZ*I_ESP_I6x_F3y_aa;
    Double I_ESP_I5xy_G3yz_aa = I_ESP_K5xyz_F3y_aa+ABZ*I_ESP_I5xy_F3y_aa;
    Double I_ESP_I5xz_G3yz_aa = I_ESP_K5x2z_F3y_aa+ABZ*I_ESP_I5xz_F3y_aa;
    Double I_ESP_I4x2y_G3yz_aa = I_ESP_K4x2yz_F3y_aa+ABZ*I_ESP_I4x2y_F3y_aa;
    Double I_ESP_I4xyz_G3yz_aa = I_ESP_K4xy2z_F3y_aa+ABZ*I_ESP_I4xyz_F3y_aa;
    Double I_ESP_I4x2z_G3yz_aa = I_ESP_K4x3z_F3y_aa+ABZ*I_ESP_I4x2z_F3y_aa;
    Double I_ESP_I3x3y_G3yz_aa = I_ESP_K3x3yz_F3y_aa+ABZ*I_ESP_I3x3y_F3y_aa;
    Double I_ESP_I3x2yz_G3yz_aa = I_ESP_K3x2y2z_F3y_aa+ABZ*I_ESP_I3x2yz_F3y_aa;
    Double I_ESP_I3xy2z_G3yz_aa = I_ESP_K3xy3z_F3y_aa+ABZ*I_ESP_I3xy2z_F3y_aa;
    Double I_ESP_I3x3z_G3yz_aa = I_ESP_K3x4z_F3y_aa+ABZ*I_ESP_I3x3z_F3y_aa;
    Double I_ESP_I2x4y_G3yz_aa = I_ESP_K2x4yz_F3y_aa+ABZ*I_ESP_I2x4y_F3y_aa;
    Double I_ESP_I2x3yz_G3yz_aa = I_ESP_K2x3y2z_F3y_aa+ABZ*I_ESP_I2x3yz_F3y_aa;
    Double I_ESP_I2x2y2z_G3yz_aa = I_ESP_K2x2y3z_F3y_aa+ABZ*I_ESP_I2x2y2z_F3y_aa;
    Double I_ESP_I2xy3z_G3yz_aa = I_ESP_K2xy4z_F3y_aa+ABZ*I_ESP_I2xy3z_F3y_aa;
    Double I_ESP_I2x4z_G3yz_aa = I_ESP_K2x5z_F3y_aa+ABZ*I_ESP_I2x4z_F3y_aa;
    Double I_ESP_Ix5y_G3yz_aa = I_ESP_Kx5yz_F3y_aa+ABZ*I_ESP_Ix5y_F3y_aa;
    Double I_ESP_Ix4yz_G3yz_aa = I_ESP_Kx4y2z_F3y_aa+ABZ*I_ESP_Ix4yz_F3y_aa;
    Double I_ESP_Ix3y2z_G3yz_aa = I_ESP_Kx3y3z_F3y_aa+ABZ*I_ESP_Ix3y2z_F3y_aa;
    Double I_ESP_Ix2y3z_G3yz_aa = I_ESP_Kx2y4z_F3y_aa+ABZ*I_ESP_Ix2y3z_F3y_aa;
    Double I_ESP_Ixy4z_G3yz_aa = I_ESP_Kxy5z_F3y_aa+ABZ*I_ESP_Ixy4z_F3y_aa;
    Double I_ESP_Ix5z_G3yz_aa = I_ESP_Kx6z_F3y_aa+ABZ*I_ESP_Ix5z_F3y_aa;
    Double I_ESP_I6y_G3yz_aa = I_ESP_K6yz_F3y_aa+ABZ*I_ESP_I6y_F3y_aa;
    Double I_ESP_I5yz_G3yz_aa = I_ESP_K5y2z_F3y_aa+ABZ*I_ESP_I5yz_F3y_aa;
    Double I_ESP_I4y2z_G3yz_aa = I_ESP_K4y3z_F3y_aa+ABZ*I_ESP_I4y2z_F3y_aa;
    Double I_ESP_I3y3z_G3yz_aa = I_ESP_K3y4z_F3y_aa+ABZ*I_ESP_I3y3z_F3y_aa;
    Double I_ESP_I2y4z_G3yz_aa = I_ESP_K2y5z_F3y_aa+ABZ*I_ESP_I2y4z_F3y_aa;
    Double I_ESP_Iy5z_G3yz_aa = I_ESP_Ky6z_F3y_aa+ABZ*I_ESP_Iy5z_F3y_aa;
    Double I_ESP_I6z_G3yz_aa = I_ESP_K7z_F3y_aa+ABZ*I_ESP_I6z_F3y_aa;
    Double I_ESP_I6x_G2y2z_aa = I_ESP_K6xz_F2yz_aa+ABZ*I_ESP_I6x_F2yz_aa;
    Double I_ESP_I5xy_G2y2z_aa = I_ESP_K5xyz_F2yz_aa+ABZ*I_ESP_I5xy_F2yz_aa;
    Double I_ESP_I5xz_G2y2z_aa = I_ESP_K5x2z_F2yz_aa+ABZ*I_ESP_I5xz_F2yz_aa;
    Double I_ESP_I4x2y_G2y2z_aa = I_ESP_K4x2yz_F2yz_aa+ABZ*I_ESP_I4x2y_F2yz_aa;
    Double I_ESP_I4xyz_G2y2z_aa = I_ESP_K4xy2z_F2yz_aa+ABZ*I_ESP_I4xyz_F2yz_aa;
    Double I_ESP_I4x2z_G2y2z_aa = I_ESP_K4x3z_F2yz_aa+ABZ*I_ESP_I4x2z_F2yz_aa;
    Double I_ESP_I3x3y_G2y2z_aa = I_ESP_K3x3yz_F2yz_aa+ABZ*I_ESP_I3x3y_F2yz_aa;
    Double I_ESP_I3x2yz_G2y2z_aa = I_ESP_K3x2y2z_F2yz_aa+ABZ*I_ESP_I3x2yz_F2yz_aa;
    Double I_ESP_I3xy2z_G2y2z_aa = I_ESP_K3xy3z_F2yz_aa+ABZ*I_ESP_I3xy2z_F2yz_aa;
    Double I_ESP_I3x3z_G2y2z_aa = I_ESP_K3x4z_F2yz_aa+ABZ*I_ESP_I3x3z_F2yz_aa;
    Double I_ESP_I2x4y_G2y2z_aa = I_ESP_K2x4yz_F2yz_aa+ABZ*I_ESP_I2x4y_F2yz_aa;
    Double I_ESP_I2x3yz_G2y2z_aa = I_ESP_K2x3y2z_F2yz_aa+ABZ*I_ESP_I2x3yz_F2yz_aa;
    Double I_ESP_I2x2y2z_G2y2z_aa = I_ESP_K2x2y3z_F2yz_aa+ABZ*I_ESP_I2x2y2z_F2yz_aa;
    Double I_ESP_I2xy3z_G2y2z_aa = I_ESP_K2xy4z_F2yz_aa+ABZ*I_ESP_I2xy3z_F2yz_aa;
    Double I_ESP_I2x4z_G2y2z_aa = I_ESP_K2x5z_F2yz_aa+ABZ*I_ESP_I2x4z_F2yz_aa;
    Double I_ESP_Ix5y_G2y2z_aa = I_ESP_Kx5yz_F2yz_aa+ABZ*I_ESP_Ix5y_F2yz_aa;
    Double I_ESP_Ix4yz_G2y2z_aa = I_ESP_Kx4y2z_F2yz_aa+ABZ*I_ESP_Ix4yz_F2yz_aa;
    Double I_ESP_Ix3y2z_G2y2z_aa = I_ESP_Kx3y3z_F2yz_aa+ABZ*I_ESP_Ix3y2z_F2yz_aa;
    Double I_ESP_Ix2y3z_G2y2z_aa = I_ESP_Kx2y4z_F2yz_aa+ABZ*I_ESP_Ix2y3z_F2yz_aa;
    Double I_ESP_Ixy4z_G2y2z_aa = I_ESP_Kxy5z_F2yz_aa+ABZ*I_ESP_Ixy4z_F2yz_aa;
    Double I_ESP_Ix5z_G2y2z_aa = I_ESP_Kx6z_F2yz_aa+ABZ*I_ESP_Ix5z_F2yz_aa;
    Double I_ESP_I6y_G2y2z_aa = I_ESP_K6yz_F2yz_aa+ABZ*I_ESP_I6y_F2yz_aa;
    Double I_ESP_I5yz_G2y2z_aa = I_ESP_K5y2z_F2yz_aa+ABZ*I_ESP_I5yz_F2yz_aa;
    Double I_ESP_I4y2z_G2y2z_aa = I_ESP_K4y3z_F2yz_aa+ABZ*I_ESP_I4y2z_F2yz_aa;
    Double I_ESP_I3y3z_G2y2z_aa = I_ESP_K3y4z_F2yz_aa+ABZ*I_ESP_I3y3z_F2yz_aa;
    Double I_ESP_I2y4z_G2y2z_aa = I_ESP_K2y5z_F2yz_aa+ABZ*I_ESP_I2y4z_F2yz_aa;
    Double I_ESP_Iy5z_G2y2z_aa = I_ESP_Ky6z_F2yz_aa+ABZ*I_ESP_Iy5z_F2yz_aa;
    Double I_ESP_I6z_G2y2z_aa = I_ESP_K7z_F2yz_aa+ABZ*I_ESP_I6z_F2yz_aa;
    Double I_ESP_I6x_Gy3z_aa = I_ESP_K6xy_F3z_aa+ABY*I_ESP_I6x_F3z_aa;
    Double I_ESP_I5xy_Gy3z_aa = I_ESP_K5x2y_F3z_aa+ABY*I_ESP_I5xy_F3z_aa;
    Double I_ESP_I5xz_Gy3z_aa = I_ESP_K5xyz_F3z_aa+ABY*I_ESP_I5xz_F3z_aa;
    Double I_ESP_I4x2y_Gy3z_aa = I_ESP_K4x3y_F3z_aa+ABY*I_ESP_I4x2y_F3z_aa;
    Double I_ESP_I4xyz_Gy3z_aa = I_ESP_K4x2yz_F3z_aa+ABY*I_ESP_I4xyz_F3z_aa;
    Double I_ESP_I4x2z_Gy3z_aa = I_ESP_K4xy2z_F3z_aa+ABY*I_ESP_I4x2z_F3z_aa;
    Double I_ESP_I3x3y_Gy3z_aa = I_ESP_K3x4y_F3z_aa+ABY*I_ESP_I3x3y_F3z_aa;
    Double I_ESP_I3x2yz_Gy3z_aa = I_ESP_K3x3yz_F3z_aa+ABY*I_ESP_I3x2yz_F3z_aa;
    Double I_ESP_I3xy2z_Gy3z_aa = I_ESP_K3x2y2z_F3z_aa+ABY*I_ESP_I3xy2z_F3z_aa;
    Double I_ESP_I3x3z_Gy3z_aa = I_ESP_K3xy3z_F3z_aa+ABY*I_ESP_I3x3z_F3z_aa;
    Double I_ESP_I2x4y_Gy3z_aa = I_ESP_K2x5y_F3z_aa+ABY*I_ESP_I2x4y_F3z_aa;
    Double I_ESP_I2x3yz_Gy3z_aa = I_ESP_K2x4yz_F3z_aa+ABY*I_ESP_I2x3yz_F3z_aa;
    Double I_ESP_I2x2y2z_Gy3z_aa = I_ESP_K2x3y2z_F3z_aa+ABY*I_ESP_I2x2y2z_F3z_aa;
    Double I_ESP_I2xy3z_Gy3z_aa = I_ESP_K2x2y3z_F3z_aa+ABY*I_ESP_I2xy3z_F3z_aa;
    Double I_ESP_I2x4z_Gy3z_aa = I_ESP_K2xy4z_F3z_aa+ABY*I_ESP_I2x4z_F3z_aa;
    Double I_ESP_Ix5y_Gy3z_aa = I_ESP_Kx6y_F3z_aa+ABY*I_ESP_Ix5y_F3z_aa;
    Double I_ESP_Ix4yz_Gy3z_aa = I_ESP_Kx5yz_F3z_aa+ABY*I_ESP_Ix4yz_F3z_aa;
    Double I_ESP_Ix3y2z_Gy3z_aa = I_ESP_Kx4y2z_F3z_aa+ABY*I_ESP_Ix3y2z_F3z_aa;
    Double I_ESP_Ix2y3z_Gy3z_aa = I_ESP_Kx3y3z_F3z_aa+ABY*I_ESP_Ix2y3z_F3z_aa;
    Double I_ESP_Ixy4z_Gy3z_aa = I_ESP_Kx2y4z_F3z_aa+ABY*I_ESP_Ixy4z_F3z_aa;
    Double I_ESP_Ix5z_Gy3z_aa = I_ESP_Kxy5z_F3z_aa+ABY*I_ESP_Ix5z_F3z_aa;
    Double I_ESP_I6y_Gy3z_aa = I_ESP_K7y_F3z_aa+ABY*I_ESP_I6y_F3z_aa;
    Double I_ESP_I5yz_Gy3z_aa = I_ESP_K6yz_F3z_aa+ABY*I_ESP_I5yz_F3z_aa;
    Double I_ESP_I4y2z_Gy3z_aa = I_ESP_K5y2z_F3z_aa+ABY*I_ESP_I4y2z_F3z_aa;
    Double I_ESP_I3y3z_Gy3z_aa = I_ESP_K4y3z_F3z_aa+ABY*I_ESP_I3y3z_F3z_aa;
    Double I_ESP_I2y4z_Gy3z_aa = I_ESP_K3y4z_F3z_aa+ABY*I_ESP_I2y4z_F3z_aa;
    Double I_ESP_Iy5z_Gy3z_aa = I_ESP_K2y5z_F3z_aa+ABY*I_ESP_Iy5z_F3z_aa;
    Double I_ESP_I6z_Gy3z_aa = I_ESP_Ky6z_F3z_aa+ABY*I_ESP_I6z_F3z_aa;
    Double I_ESP_I6x_G4z_aa = I_ESP_K6xz_F3z_aa+ABZ*I_ESP_I6x_F3z_aa;
    Double I_ESP_I5xy_G4z_aa = I_ESP_K5xyz_F3z_aa+ABZ*I_ESP_I5xy_F3z_aa;
    Double I_ESP_I5xz_G4z_aa = I_ESP_K5x2z_F3z_aa+ABZ*I_ESP_I5xz_F3z_aa;
    Double I_ESP_I4x2y_G4z_aa = I_ESP_K4x2yz_F3z_aa+ABZ*I_ESP_I4x2y_F3z_aa;
    Double I_ESP_I4xyz_G4z_aa = I_ESP_K4xy2z_F3z_aa+ABZ*I_ESP_I4xyz_F3z_aa;
    Double I_ESP_I4x2z_G4z_aa = I_ESP_K4x3z_F3z_aa+ABZ*I_ESP_I4x2z_F3z_aa;
    Double I_ESP_I3x3y_G4z_aa = I_ESP_K3x3yz_F3z_aa+ABZ*I_ESP_I3x3y_F3z_aa;
    Double I_ESP_I3x2yz_G4z_aa = I_ESP_K3x2y2z_F3z_aa+ABZ*I_ESP_I3x2yz_F3z_aa;
    Double I_ESP_I3xy2z_G4z_aa = I_ESP_K3xy3z_F3z_aa+ABZ*I_ESP_I3xy2z_F3z_aa;
    Double I_ESP_I3x3z_G4z_aa = I_ESP_K3x4z_F3z_aa+ABZ*I_ESP_I3x3z_F3z_aa;
    Double I_ESP_I2x4y_G4z_aa = I_ESP_K2x4yz_F3z_aa+ABZ*I_ESP_I2x4y_F3z_aa;
    Double I_ESP_I2x3yz_G4z_aa = I_ESP_K2x3y2z_F3z_aa+ABZ*I_ESP_I2x3yz_F3z_aa;
    Double I_ESP_I2x2y2z_G4z_aa = I_ESP_K2x2y3z_F3z_aa+ABZ*I_ESP_I2x2y2z_F3z_aa;
    Double I_ESP_I2xy3z_G4z_aa = I_ESP_K2xy4z_F3z_aa+ABZ*I_ESP_I2xy3z_F3z_aa;
    Double I_ESP_I2x4z_G4z_aa = I_ESP_K2x5z_F3z_aa+ABZ*I_ESP_I2x4z_F3z_aa;
    Double I_ESP_Ix5y_G4z_aa = I_ESP_Kx5yz_F3z_aa+ABZ*I_ESP_Ix5y_F3z_aa;
    Double I_ESP_Ix4yz_G4z_aa = I_ESP_Kx4y2z_F3z_aa+ABZ*I_ESP_Ix4yz_F3z_aa;
    Double I_ESP_Ix3y2z_G4z_aa = I_ESP_Kx3y3z_F3z_aa+ABZ*I_ESP_Ix3y2z_F3z_aa;
    Double I_ESP_Ix2y3z_G4z_aa = I_ESP_Kx2y4z_F3z_aa+ABZ*I_ESP_Ix2y3z_F3z_aa;
    Double I_ESP_Ixy4z_G4z_aa = I_ESP_Kxy5z_F3z_aa+ABZ*I_ESP_Ixy4z_F3z_aa;
    Double I_ESP_Ix5z_G4z_aa = I_ESP_Kx6z_F3z_aa+ABZ*I_ESP_Ix5z_F3z_aa;
    Double I_ESP_I6y_G4z_aa = I_ESP_K6yz_F3z_aa+ABZ*I_ESP_I6y_F3z_aa;
    Double I_ESP_I5yz_G4z_aa = I_ESP_K5y2z_F3z_aa+ABZ*I_ESP_I5yz_F3z_aa;
    Double I_ESP_I4y2z_G4z_aa = I_ESP_K4y3z_F3z_aa+ABZ*I_ESP_I4y2z_F3z_aa;
    Double I_ESP_I3y3z_G4z_aa = I_ESP_K3y4z_F3z_aa+ABZ*I_ESP_I3y3z_F3z_aa;
    Double I_ESP_I2y4z_G4z_aa = I_ESP_K2y5z_F3z_aa+ABZ*I_ESP_I2y4z_F3z_aa;
    Double I_ESP_Iy5z_G4z_aa = I_ESP_Ky6z_F3z_aa+ABZ*I_ESP_Iy5z_F3z_aa;
    Double I_ESP_I6z_G4z_aa = I_ESP_K7z_F3z_aa+ABZ*I_ESP_I6z_F3z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G_aa
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_D_G
     ************************************************************/
    abcd[iGrid*1350+0] = 4.0E0*I_ESP_I6x_G4x_aa-2.0E0*4*I_ESP_G4x_G4x_a-2.0E0*5*I_ESP_G4x_G4x_a+4*3*I_ESP_D2x_G4x;
    abcd[iGrid*1350+1] = 4.0E0*I_ESP_I5xy_G4x_aa-2.0E0*3*I_ESP_G3xy_G4x_a-2.0E0*4*I_ESP_G3xy_G4x_a+3*2*I_ESP_Dxy_G4x;
    abcd[iGrid*1350+2] = 4.0E0*I_ESP_I5xz_G4x_aa-2.0E0*3*I_ESP_G3xz_G4x_a-2.0E0*4*I_ESP_G3xz_G4x_a+3*2*I_ESP_Dxz_G4x;
    abcd[iGrid*1350+3] = 4.0E0*I_ESP_I4x2y_G4x_aa-2.0E0*2*I_ESP_G2x2y_G4x_a-2.0E0*3*I_ESP_G2x2y_G4x_a+2*1*I_ESP_D2y_G4x;
    abcd[iGrid*1350+4] = 4.0E0*I_ESP_I4xyz_G4x_aa-2.0E0*2*I_ESP_G2xyz_G4x_a-2.0E0*3*I_ESP_G2xyz_G4x_a+2*1*I_ESP_Dyz_G4x;
    abcd[iGrid*1350+5] = 4.0E0*I_ESP_I4x2z_G4x_aa-2.0E0*2*I_ESP_G2x2z_G4x_a-2.0E0*3*I_ESP_G2x2z_G4x_a+2*1*I_ESP_D2z_G4x;
    abcd[iGrid*1350+6] = 4.0E0*I_ESP_I3x3y_G4x_aa-2.0E0*1*I_ESP_Gx3y_G4x_a-2.0E0*2*I_ESP_Gx3y_G4x_a;
    abcd[iGrid*1350+7] = 4.0E0*I_ESP_I3x2yz_G4x_aa-2.0E0*1*I_ESP_Gx2yz_G4x_a-2.0E0*2*I_ESP_Gx2yz_G4x_a;
    abcd[iGrid*1350+8] = 4.0E0*I_ESP_I3xy2z_G4x_aa-2.0E0*1*I_ESP_Gxy2z_G4x_a-2.0E0*2*I_ESP_Gxy2z_G4x_a;
    abcd[iGrid*1350+9] = 4.0E0*I_ESP_I3x3z_G4x_aa-2.0E0*1*I_ESP_Gx3z_G4x_a-2.0E0*2*I_ESP_Gx3z_G4x_a;
    abcd[iGrid*1350+10] = 4.0E0*I_ESP_I2x4y_G4x_aa-2.0E0*1*I_ESP_G4y_G4x_a;
    abcd[iGrid*1350+11] = 4.0E0*I_ESP_I2x3yz_G4x_aa-2.0E0*1*I_ESP_G3yz_G4x_a;
    abcd[iGrid*1350+12] = 4.0E0*I_ESP_I2x2y2z_G4x_aa-2.0E0*1*I_ESP_G2y2z_G4x_a;
    abcd[iGrid*1350+13] = 4.0E0*I_ESP_I2xy3z_G4x_aa-2.0E0*1*I_ESP_Gy3z_G4x_a;
    abcd[iGrid*1350+14] = 4.0E0*I_ESP_I2x4z_G4x_aa-2.0E0*1*I_ESP_G4z_G4x_a;
    abcd[iGrid*1350+15] = 4.0E0*I_ESP_I6x_G3xy_aa-2.0E0*4*I_ESP_G4x_G3xy_a-2.0E0*5*I_ESP_G4x_G3xy_a+4*3*I_ESP_D2x_G3xy;
    abcd[iGrid*1350+16] = 4.0E0*I_ESP_I5xy_G3xy_aa-2.0E0*3*I_ESP_G3xy_G3xy_a-2.0E0*4*I_ESP_G3xy_G3xy_a+3*2*I_ESP_Dxy_G3xy;
    abcd[iGrid*1350+17] = 4.0E0*I_ESP_I5xz_G3xy_aa-2.0E0*3*I_ESP_G3xz_G3xy_a-2.0E0*4*I_ESP_G3xz_G3xy_a+3*2*I_ESP_Dxz_G3xy;
    abcd[iGrid*1350+18] = 4.0E0*I_ESP_I4x2y_G3xy_aa-2.0E0*2*I_ESP_G2x2y_G3xy_a-2.0E0*3*I_ESP_G2x2y_G3xy_a+2*1*I_ESP_D2y_G3xy;
    abcd[iGrid*1350+19] = 4.0E0*I_ESP_I4xyz_G3xy_aa-2.0E0*2*I_ESP_G2xyz_G3xy_a-2.0E0*3*I_ESP_G2xyz_G3xy_a+2*1*I_ESP_Dyz_G3xy;
    abcd[iGrid*1350+20] = 4.0E0*I_ESP_I4x2z_G3xy_aa-2.0E0*2*I_ESP_G2x2z_G3xy_a-2.0E0*3*I_ESP_G2x2z_G3xy_a+2*1*I_ESP_D2z_G3xy;
    abcd[iGrid*1350+21] = 4.0E0*I_ESP_I3x3y_G3xy_aa-2.0E0*1*I_ESP_Gx3y_G3xy_a-2.0E0*2*I_ESP_Gx3y_G3xy_a;
    abcd[iGrid*1350+22] = 4.0E0*I_ESP_I3x2yz_G3xy_aa-2.0E0*1*I_ESP_Gx2yz_G3xy_a-2.0E0*2*I_ESP_Gx2yz_G3xy_a;
    abcd[iGrid*1350+23] = 4.0E0*I_ESP_I3xy2z_G3xy_aa-2.0E0*1*I_ESP_Gxy2z_G3xy_a-2.0E0*2*I_ESP_Gxy2z_G3xy_a;
    abcd[iGrid*1350+24] = 4.0E0*I_ESP_I3x3z_G3xy_aa-2.0E0*1*I_ESP_Gx3z_G3xy_a-2.0E0*2*I_ESP_Gx3z_G3xy_a;
    abcd[iGrid*1350+25] = 4.0E0*I_ESP_I2x4y_G3xy_aa-2.0E0*1*I_ESP_G4y_G3xy_a;
    abcd[iGrid*1350+26] = 4.0E0*I_ESP_I2x3yz_G3xy_aa-2.0E0*1*I_ESP_G3yz_G3xy_a;
    abcd[iGrid*1350+27] = 4.0E0*I_ESP_I2x2y2z_G3xy_aa-2.0E0*1*I_ESP_G2y2z_G3xy_a;
    abcd[iGrid*1350+28] = 4.0E0*I_ESP_I2xy3z_G3xy_aa-2.0E0*1*I_ESP_Gy3z_G3xy_a;
    abcd[iGrid*1350+29] = 4.0E0*I_ESP_I2x4z_G3xy_aa-2.0E0*1*I_ESP_G4z_G3xy_a;
    abcd[iGrid*1350+30] = 4.0E0*I_ESP_I6x_G3xz_aa-2.0E0*4*I_ESP_G4x_G3xz_a-2.0E0*5*I_ESP_G4x_G3xz_a+4*3*I_ESP_D2x_G3xz;
    abcd[iGrid*1350+31] = 4.0E0*I_ESP_I5xy_G3xz_aa-2.0E0*3*I_ESP_G3xy_G3xz_a-2.0E0*4*I_ESP_G3xy_G3xz_a+3*2*I_ESP_Dxy_G3xz;
    abcd[iGrid*1350+32] = 4.0E0*I_ESP_I5xz_G3xz_aa-2.0E0*3*I_ESP_G3xz_G3xz_a-2.0E0*4*I_ESP_G3xz_G3xz_a+3*2*I_ESP_Dxz_G3xz;
    abcd[iGrid*1350+33] = 4.0E0*I_ESP_I4x2y_G3xz_aa-2.0E0*2*I_ESP_G2x2y_G3xz_a-2.0E0*3*I_ESP_G2x2y_G3xz_a+2*1*I_ESP_D2y_G3xz;
    abcd[iGrid*1350+34] = 4.0E0*I_ESP_I4xyz_G3xz_aa-2.0E0*2*I_ESP_G2xyz_G3xz_a-2.0E0*3*I_ESP_G2xyz_G3xz_a+2*1*I_ESP_Dyz_G3xz;
    abcd[iGrid*1350+35] = 4.0E0*I_ESP_I4x2z_G3xz_aa-2.0E0*2*I_ESP_G2x2z_G3xz_a-2.0E0*3*I_ESP_G2x2z_G3xz_a+2*1*I_ESP_D2z_G3xz;
    abcd[iGrid*1350+36] = 4.0E0*I_ESP_I3x3y_G3xz_aa-2.0E0*1*I_ESP_Gx3y_G3xz_a-2.0E0*2*I_ESP_Gx3y_G3xz_a;
    abcd[iGrid*1350+37] = 4.0E0*I_ESP_I3x2yz_G3xz_aa-2.0E0*1*I_ESP_Gx2yz_G3xz_a-2.0E0*2*I_ESP_Gx2yz_G3xz_a;
    abcd[iGrid*1350+38] = 4.0E0*I_ESP_I3xy2z_G3xz_aa-2.0E0*1*I_ESP_Gxy2z_G3xz_a-2.0E0*2*I_ESP_Gxy2z_G3xz_a;
    abcd[iGrid*1350+39] = 4.0E0*I_ESP_I3x3z_G3xz_aa-2.0E0*1*I_ESP_Gx3z_G3xz_a-2.0E0*2*I_ESP_Gx3z_G3xz_a;
    abcd[iGrid*1350+40] = 4.0E0*I_ESP_I2x4y_G3xz_aa-2.0E0*1*I_ESP_G4y_G3xz_a;
    abcd[iGrid*1350+41] = 4.0E0*I_ESP_I2x3yz_G3xz_aa-2.0E0*1*I_ESP_G3yz_G3xz_a;
    abcd[iGrid*1350+42] = 4.0E0*I_ESP_I2x2y2z_G3xz_aa-2.0E0*1*I_ESP_G2y2z_G3xz_a;
    abcd[iGrid*1350+43] = 4.0E0*I_ESP_I2xy3z_G3xz_aa-2.0E0*1*I_ESP_Gy3z_G3xz_a;
    abcd[iGrid*1350+44] = 4.0E0*I_ESP_I2x4z_G3xz_aa-2.0E0*1*I_ESP_G4z_G3xz_a;
    abcd[iGrid*1350+45] = 4.0E0*I_ESP_I6x_G2x2y_aa-2.0E0*4*I_ESP_G4x_G2x2y_a-2.0E0*5*I_ESP_G4x_G2x2y_a+4*3*I_ESP_D2x_G2x2y;
    abcd[iGrid*1350+46] = 4.0E0*I_ESP_I5xy_G2x2y_aa-2.0E0*3*I_ESP_G3xy_G2x2y_a-2.0E0*4*I_ESP_G3xy_G2x2y_a+3*2*I_ESP_Dxy_G2x2y;
    abcd[iGrid*1350+47] = 4.0E0*I_ESP_I5xz_G2x2y_aa-2.0E0*3*I_ESP_G3xz_G2x2y_a-2.0E0*4*I_ESP_G3xz_G2x2y_a+3*2*I_ESP_Dxz_G2x2y;
    abcd[iGrid*1350+48] = 4.0E0*I_ESP_I4x2y_G2x2y_aa-2.0E0*2*I_ESP_G2x2y_G2x2y_a-2.0E0*3*I_ESP_G2x2y_G2x2y_a+2*1*I_ESP_D2y_G2x2y;
    abcd[iGrid*1350+49] = 4.0E0*I_ESP_I4xyz_G2x2y_aa-2.0E0*2*I_ESP_G2xyz_G2x2y_a-2.0E0*3*I_ESP_G2xyz_G2x2y_a+2*1*I_ESP_Dyz_G2x2y;
    abcd[iGrid*1350+50] = 4.0E0*I_ESP_I4x2z_G2x2y_aa-2.0E0*2*I_ESP_G2x2z_G2x2y_a-2.0E0*3*I_ESP_G2x2z_G2x2y_a+2*1*I_ESP_D2z_G2x2y;
    abcd[iGrid*1350+51] = 4.0E0*I_ESP_I3x3y_G2x2y_aa-2.0E0*1*I_ESP_Gx3y_G2x2y_a-2.0E0*2*I_ESP_Gx3y_G2x2y_a;
    abcd[iGrid*1350+52] = 4.0E0*I_ESP_I3x2yz_G2x2y_aa-2.0E0*1*I_ESP_Gx2yz_G2x2y_a-2.0E0*2*I_ESP_Gx2yz_G2x2y_a;
    abcd[iGrid*1350+53] = 4.0E0*I_ESP_I3xy2z_G2x2y_aa-2.0E0*1*I_ESP_Gxy2z_G2x2y_a-2.0E0*2*I_ESP_Gxy2z_G2x2y_a;
    abcd[iGrid*1350+54] = 4.0E0*I_ESP_I3x3z_G2x2y_aa-2.0E0*1*I_ESP_Gx3z_G2x2y_a-2.0E0*2*I_ESP_Gx3z_G2x2y_a;
    abcd[iGrid*1350+55] = 4.0E0*I_ESP_I2x4y_G2x2y_aa-2.0E0*1*I_ESP_G4y_G2x2y_a;
    abcd[iGrid*1350+56] = 4.0E0*I_ESP_I2x3yz_G2x2y_aa-2.0E0*1*I_ESP_G3yz_G2x2y_a;
    abcd[iGrid*1350+57] = 4.0E0*I_ESP_I2x2y2z_G2x2y_aa-2.0E0*1*I_ESP_G2y2z_G2x2y_a;
    abcd[iGrid*1350+58] = 4.0E0*I_ESP_I2xy3z_G2x2y_aa-2.0E0*1*I_ESP_Gy3z_G2x2y_a;
    abcd[iGrid*1350+59] = 4.0E0*I_ESP_I2x4z_G2x2y_aa-2.0E0*1*I_ESP_G4z_G2x2y_a;
    abcd[iGrid*1350+60] = 4.0E0*I_ESP_I6x_G2xyz_aa-2.0E0*4*I_ESP_G4x_G2xyz_a-2.0E0*5*I_ESP_G4x_G2xyz_a+4*3*I_ESP_D2x_G2xyz;
    abcd[iGrid*1350+61] = 4.0E0*I_ESP_I5xy_G2xyz_aa-2.0E0*3*I_ESP_G3xy_G2xyz_a-2.0E0*4*I_ESP_G3xy_G2xyz_a+3*2*I_ESP_Dxy_G2xyz;
    abcd[iGrid*1350+62] = 4.0E0*I_ESP_I5xz_G2xyz_aa-2.0E0*3*I_ESP_G3xz_G2xyz_a-2.0E0*4*I_ESP_G3xz_G2xyz_a+3*2*I_ESP_Dxz_G2xyz;
    abcd[iGrid*1350+63] = 4.0E0*I_ESP_I4x2y_G2xyz_aa-2.0E0*2*I_ESP_G2x2y_G2xyz_a-2.0E0*3*I_ESP_G2x2y_G2xyz_a+2*1*I_ESP_D2y_G2xyz;
    abcd[iGrid*1350+64] = 4.0E0*I_ESP_I4xyz_G2xyz_aa-2.0E0*2*I_ESP_G2xyz_G2xyz_a-2.0E0*3*I_ESP_G2xyz_G2xyz_a+2*1*I_ESP_Dyz_G2xyz;
    abcd[iGrid*1350+65] = 4.0E0*I_ESP_I4x2z_G2xyz_aa-2.0E0*2*I_ESP_G2x2z_G2xyz_a-2.0E0*3*I_ESP_G2x2z_G2xyz_a+2*1*I_ESP_D2z_G2xyz;
    abcd[iGrid*1350+66] = 4.0E0*I_ESP_I3x3y_G2xyz_aa-2.0E0*1*I_ESP_Gx3y_G2xyz_a-2.0E0*2*I_ESP_Gx3y_G2xyz_a;
    abcd[iGrid*1350+67] = 4.0E0*I_ESP_I3x2yz_G2xyz_aa-2.0E0*1*I_ESP_Gx2yz_G2xyz_a-2.0E0*2*I_ESP_Gx2yz_G2xyz_a;
    abcd[iGrid*1350+68] = 4.0E0*I_ESP_I3xy2z_G2xyz_aa-2.0E0*1*I_ESP_Gxy2z_G2xyz_a-2.0E0*2*I_ESP_Gxy2z_G2xyz_a;
    abcd[iGrid*1350+69] = 4.0E0*I_ESP_I3x3z_G2xyz_aa-2.0E0*1*I_ESP_Gx3z_G2xyz_a-2.0E0*2*I_ESP_Gx3z_G2xyz_a;
    abcd[iGrid*1350+70] = 4.0E0*I_ESP_I2x4y_G2xyz_aa-2.0E0*1*I_ESP_G4y_G2xyz_a;
    abcd[iGrid*1350+71] = 4.0E0*I_ESP_I2x3yz_G2xyz_aa-2.0E0*1*I_ESP_G3yz_G2xyz_a;
    abcd[iGrid*1350+72] = 4.0E0*I_ESP_I2x2y2z_G2xyz_aa-2.0E0*1*I_ESP_G2y2z_G2xyz_a;
    abcd[iGrid*1350+73] = 4.0E0*I_ESP_I2xy3z_G2xyz_aa-2.0E0*1*I_ESP_Gy3z_G2xyz_a;
    abcd[iGrid*1350+74] = 4.0E0*I_ESP_I2x4z_G2xyz_aa-2.0E0*1*I_ESP_G4z_G2xyz_a;
    abcd[iGrid*1350+75] = 4.0E0*I_ESP_I6x_G2x2z_aa-2.0E0*4*I_ESP_G4x_G2x2z_a-2.0E0*5*I_ESP_G4x_G2x2z_a+4*3*I_ESP_D2x_G2x2z;
    abcd[iGrid*1350+76] = 4.0E0*I_ESP_I5xy_G2x2z_aa-2.0E0*3*I_ESP_G3xy_G2x2z_a-2.0E0*4*I_ESP_G3xy_G2x2z_a+3*2*I_ESP_Dxy_G2x2z;
    abcd[iGrid*1350+77] = 4.0E0*I_ESP_I5xz_G2x2z_aa-2.0E0*3*I_ESP_G3xz_G2x2z_a-2.0E0*4*I_ESP_G3xz_G2x2z_a+3*2*I_ESP_Dxz_G2x2z;
    abcd[iGrid*1350+78] = 4.0E0*I_ESP_I4x2y_G2x2z_aa-2.0E0*2*I_ESP_G2x2y_G2x2z_a-2.0E0*3*I_ESP_G2x2y_G2x2z_a+2*1*I_ESP_D2y_G2x2z;
    abcd[iGrid*1350+79] = 4.0E0*I_ESP_I4xyz_G2x2z_aa-2.0E0*2*I_ESP_G2xyz_G2x2z_a-2.0E0*3*I_ESP_G2xyz_G2x2z_a+2*1*I_ESP_Dyz_G2x2z;
    abcd[iGrid*1350+80] = 4.0E0*I_ESP_I4x2z_G2x2z_aa-2.0E0*2*I_ESP_G2x2z_G2x2z_a-2.0E0*3*I_ESP_G2x2z_G2x2z_a+2*1*I_ESP_D2z_G2x2z;
    abcd[iGrid*1350+81] = 4.0E0*I_ESP_I3x3y_G2x2z_aa-2.0E0*1*I_ESP_Gx3y_G2x2z_a-2.0E0*2*I_ESP_Gx3y_G2x2z_a;
    abcd[iGrid*1350+82] = 4.0E0*I_ESP_I3x2yz_G2x2z_aa-2.0E0*1*I_ESP_Gx2yz_G2x2z_a-2.0E0*2*I_ESP_Gx2yz_G2x2z_a;
    abcd[iGrid*1350+83] = 4.0E0*I_ESP_I3xy2z_G2x2z_aa-2.0E0*1*I_ESP_Gxy2z_G2x2z_a-2.0E0*2*I_ESP_Gxy2z_G2x2z_a;
    abcd[iGrid*1350+84] = 4.0E0*I_ESP_I3x3z_G2x2z_aa-2.0E0*1*I_ESP_Gx3z_G2x2z_a-2.0E0*2*I_ESP_Gx3z_G2x2z_a;
    abcd[iGrid*1350+85] = 4.0E0*I_ESP_I2x4y_G2x2z_aa-2.0E0*1*I_ESP_G4y_G2x2z_a;
    abcd[iGrid*1350+86] = 4.0E0*I_ESP_I2x3yz_G2x2z_aa-2.0E0*1*I_ESP_G3yz_G2x2z_a;
    abcd[iGrid*1350+87] = 4.0E0*I_ESP_I2x2y2z_G2x2z_aa-2.0E0*1*I_ESP_G2y2z_G2x2z_a;
    abcd[iGrid*1350+88] = 4.0E0*I_ESP_I2xy3z_G2x2z_aa-2.0E0*1*I_ESP_Gy3z_G2x2z_a;
    abcd[iGrid*1350+89] = 4.0E0*I_ESP_I2x4z_G2x2z_aa-2.0E0*1*I_ESP_G4z_G2x2z_a;
    abcd[iGrid*1350+90] = 4.0E0*I_ESP_I6x_Gx3y_aa-2.0E0*4*I_ESP_G4x_Gx3y_a-2.0E0*5*I_ESP_G4x_Gx3y_a+4*3*I_ESP_D2x_Gx3y;
    abcd[iGrid*1350+91] = 4.0E0*I_ESP_I5xy_Gx3y_aa-2.0E0*3*I_ESP_G3xy_Gx3y_a-2.0E0*4*I_ESP_G3xy_Gx3y_a+3*2*I_ESP_Dxy_Gx3y;
    abcd[iGrid*1350+92] = 4.0E0*I_ESP_I5xz_Gx3y_aa-2.0E0*3*I_ESP_G3xz_Gx3y_a-2.0E0*4*I_ESP_G3xz_Gx3y_a+3*2*I_ESP_Dxz_Gx3y;
    abcd[iGrid*1350+93] = 4.0E0*I_ESP_I4x2y_Gx3y_aa-2.0E0*2*I_ESP_G2x2y_Gx3y_a-2.0E0*3*I_ESP_G2x2y_Gx3y_a+2*1*I_ESP_D2y_Gx3y;
    abcd[iGrid*1350+94] = 4.0E0*I_ESP_I4xyz_Gx3y_aa-2.0E0*2*I_ESP_G2xyz_Gx3y_a-2.0E0*3*I_ESP_G2xyz_Gx3y_a+2*1*I_ESP_Dyz_Gx3y;
    abcd[iGrid*1350+95] = 4.0E0*I_ESP_I4x2z_Gx3y_aa-2.0E0*2*I_ESP_G2x2z_Gx3y_a-2.0E0*3*I_ESP_G2x2z_Gx3y_a+2*1*I_ESP_D2z_Gx3y;
    abcd[iGrid*1350+96] = 4.0E0*I_ESP_I3x3y_Gx3y_aa-2.0E0*1*I_ESP_Gx3y_Gx3y_a-2.0E0*2*I_ESP_Gx3y_Gx3y_a;
    abcd[iGrid*1350+97] = 4.0E0*I_ESP_I3x2yz_Gx3y_aa-2.0E0*1*I_ESP_Gx2yz_Gx3y_a-2.0E0*2*I_ESP_Gx2yz_Gx3y_a;
    abcd[iGrid*1350+98] = 4.0E0*I_ESP_I3xy2z_Gx3y_aa-2.0E0*1*I_ESP_Gxy2z_Gx3y_a-2.0E0*2*I_ESP_Gxy2z_Gx3y_a;
    abcd[iGrid*1350+99] = 4.0E0*I_ESP_I3x3z_Gx3y_aa-2.0E0*1*I_ESP_Gx3z_Gx3y_a-2.0E0*2*I_ESP_Gx3z_Gx3y_a;
    abcd[iGrid*1350+100] = 4.0E0*I_ESP_I2x4y_Gx3y_aa-2.0E0*1*I_ESP_G4y_Gx3y_a;
    abcd[iGrid*1350+101] = 4.0E0*I_ESP_I2x3yz_Gx3y_aa-2.0E0*1*I_ESP_G3yz_Gx3y_a;
    abcd[iGrid*1350+102] = 4.0E0*I_ESP_I2x2y2z_Gx3y_aa-2.0E0*1*I_ESP_G2y2z_Gx3y_a;
    abcd[iGrid*1350+103] = 4.0E0*I_ESP_I2xy3z_Gx3y_aa-2.0E0*1*I_ESP_Gy3z_Gx3y_a;
    abcd[iGrid*1350+104] = 4.0E0*I_ESP_I2x4z_Gx3y_aa-2.0E0*1*I_ESP_G4z_Gx3y_a;
    abcd[iGrid*1350+105] = 4.0E0*I_ESP_I6x_Gx2yz_aa-2.0E0*4*I_ESP_G4x_Gx2yz_a-2.0E0*5*I_ESP_G4x_Gx2yz_a+4*3*I_ESP_D2x_Gx2yz;
    abcd[iGrid*1350+106] = 4.0E0*I_ESP_I5xy_Gx2yz_aa-2.0E0*3*I_ESP_G3xy_Gx2yz_a-2.0E0*4*I_ESP_G3xy_Gx2yz_a+3*2*I_ESP_Dxy_Gx2yz;
    abcd[iGrid*1350+107] = 4.0E0*I_ESP_I5xz_Gx2yz_aa-2.0E0*3*I_ESP_G3xz_Gx2yz_a-2.0E0*4*I_ESP_G3xz_Gx2yz_a+3*2*I_ESP_Dxz_Gx2yz;
    abcd[iGrid*1350+108] = 4.0E0*I_ESP_I4x2y_Gx2yz_aa-2.0E0*2*I_ESP_G2x2y_Gx2yz_a-2.0E0*3*I_ESP_G2x2y_Gx2yz_a+2*1*I_ESP_D2y_Gx2yz;
    abcd[iGrid*1350+109] = 4.0E0*I_ESP_I4xyz_Gx2yz_aa-2.0E0*2*I_ESP_G2xyz_Gx2yz_a-2.0E0*3*I_ESP_G2xyz_Gx2yz_a+2*1*I_ESP_Dyz_Gx2yz;
    abcd[iGrid*1350+110] = 4.0E0*I_ESP_I4x2z_Gx2yz_aa-2.0E0*2*I_ESP_G2x2z_Gx2yz_a-2.0E0*3*I_ESP_G2x2z_Gx2yz_a+2*1*I_ESP_D2z_Gx2yz;
    abcd[iGrid*1350+111] = 4.0E0*I_ESP_I3x3y_Gx2yz_aa-2.0E0*1*I_ESP_Gx3y_Gx2yz_a-2.0E0*2*I_ESP_Gx3y_Gx2yz_a;
    abcd[iGrid*1350+112] = 4.0E0*I_ESP_I3x2yz_Gx2yz_aa-2.0E0*1*I_ESP_Gx2yz_Gx2yz_a-2.0E0*2*I_ESP_Gx2yz_Gx2yz_a;
    abcd[iGrid*1350+113] = 4.0E0*I_ESP_I3xy2z_Gx2yz_aa-2.0E0*1*I_ESP_Gxy2z_Gx2yz_a-2.0E0*2*I_ESP_Gxy2z_Gx2yz_a;
    abcd[iGrid*1350+114] = 4.0E0*I_ESP_I3x3z_Gx2yz_aa-2.0E0*1*I_ESP_Gx3z_Gx2yz_a-2.0E0*2*I_ESP_Gx3z_Gx2yz_a;
    abcd[iGrid*1350+115] = 4.0E0*I_ESP_I2x4y_Gx2yz_aa-2.0E0*1*I_ESP_G4y_Gx2yz_a;
    abcd[iGrid*1350+116] = 4.0E0*I_ESP_I2x3yz_Gx2yz_aa-2.0E0*1*I_ESP_G3yz_Gx2yz_a;
    abcd[iGrid*1350+117] = 4.0E0*I_ESP_I2x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_G2y2z_Gx2yz_a;
    abcd[iGrid*1350+118] = 4.0E0*I_ESP_I2xy3z_Gx2yz_aa-2.0E0*1*I_ESP_Gy3z_Gx2yz_a;
    abcd[iGrid*1350+119] = 4.0E0*I_ESP_I2x4z_Gx2yz_aa-2.0E0*1*I_ESP_G4z_Gx2yz_a;
    abcd[iGrid*1350+120] = 4.0E0*I_ESP_I6x_Gxy2z_aa-2.0E0*4*I_ESP_G4x_Gxy2z_a-2.0E0*5*I_ESP_G4x_Gxy2z_a+4*3*I_ESP_D2x_Gxy2z;
    abcd[iGrid*1350+121] = 4.0E0*I_ESP_I5xy_Gxy2z_aa-2.0E0*3*I_ESP_G3xy_Gxy2z_a-2.0E0*4*I_ESP_G3xy_Gxy2z_a+3*2*I_ESP_Dxy_Gxy2z;
    abcd[iGrid*1350+122] = 4.0E0*I_ESP_I5xz_Gxy2z_aa-2.0E0*3*I_ESP_G3xz_Gxy2z_a-2.0E0*4*I_ESP_G3xz_Gxy2z_a+3*2*I_ESP_Dxz_Gxy2z;
    abcd[iGrid*1350+123] = 4.0E0*I_ESP_I4x2y_Gxy2z_aa-2.0E0*2*I_ESP_G2x2y_Gxy2z_a-2.0E0*3*I_ESP_G2x2y_Gxy2z_a+2*1*I_ESP_D2y_Gxy2z;
    abcd[iGrid*1350+124] = 4.0E0*I_ESP_I4xyz_Gxy2z_aa-2.0E0*2*I_ESP_G2xyz_Gxy2z_a-2.0E0*3*I_ESP_G2xyz_Gxy2z_a+2*1*I_ESP_Dyz_Gxy2z;
    abcd[iGrid*1350+125] = 4.0E0*I_ESP_I4x2z_Gxy2z_aa-2.0E0*2*I_ESP_G2x2z_Gxy2z_a-2.0E0*3*I_ESP_G2x2z_Gxy2z_a+2*1*I_ESP_D2z_Gxy2z;
    abcd[iGrid*1350+126] = 4.0E0*I_ESP_I3x3y_Gxy2z_aa-2.0E0*1*I_ESP_Gx3y_Gxy2z_a-2.0E0*2*I_ESP_Gx3y_Gxy2z_a;
    abcd[iGrid*1350+127] = 4.0E0*I_ESP_I3x2yz_Gxy2z_aa-2.0E0*1*I_ESP_Gx2yz_Gxy2z_a-2.0E0*2*I_ESP_Gx2yz_Gxy2z_a;
    abcd[iGrid*1350+128] = 4.0E0*I_ESP_I3xy2z_Gxy2z_aa-2.0E0*1*I_ESP_Gxy2z_Gxy2z_a-2.0E0*2*I_ESP_Gxy2z_Gxy2z_a;
    abcd[iGrid*1350+129] = 4.0E0*I_ESP_I3x3z_Gxy2z_aa-2.0E0*1*I_ESP_Gx3z_Gxy2z_a-2.0E0*2*I_ESP_Gx3z_Gxy2z_a;
    abcd[iGrid*1350+130] = 4.0E0*I_ESP_I2x4y_Gxy2z_aa-2.0E0*1*I_ESP_G4y_Gxy2z_a;
    abcd[iGrid*1350+131] = 4.0E0*I_ESP_I2x3yz_Gxy2z_aa-2.0E0*1*I_ESP_G3yz_Gxy2z_a;
    abcd[iGrid*1350+132] = 4.0E0*I_ESP_I2x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_G2y2z_Gxy2z_a;
    abcd[iGrid*1350+133] = 4.0E0*I_ESP_I2xy3z_Gxy2z_aa-2.0E0*1*I_ESP_Gy3z_Gxy2z_a;
    abcd[iGrid*1350+134] = 4.0E0*I_ESP_I2x4z_Gxy2z_aa-2.0E0*1*I_ESP_G4z_Gxy2z_a;
    abcd[iGrid*1350+135] = 4.0E0*I_ESP_I6x_Gx3z_aa-2.0E0*4*I_ESP_G4x_Gx3z_a-2.0E0*5*I_ESP_G4x_Gx3z_a+4*3*I_ESP_D2x_Gx3z;
    abcd[iGrid*1350+136] = 4.0E0*I_ESP_I5xy_Gx3z_aa-2.0E0*3*I_ESP_G3xy_Gx3z_a-2.0E0*4*I_ESP_G3xy_Gx3z_a+3*2*I_ESP_Dxy_Gx3z;
    abcd[iGrid*1350+137] = 4.0E0*I_ESP_I5xz_Gx3z_aa-2.0E0*3*I_ESP_G3xz_Gx3z_a-2.0E0*4*I_ESP_G3xz_Gx3z_a+3*2*I_ESP_Dxz_Gx3z;
    abcd[iGrid*1350+138] = 4.0E0*I_ESP_I4x2y_Gx3z_aa-2.0E0*2*I_ESP_G2x2y_Gx3z_a-2.0E0*3*I_ESP_G2x2y_Gx3z_a+2*1*I_ESP_D2y_Gx3z;
    abcd[iGrid*1350+139] = 4.0E0*I_ESP_I4xyz_Gx3z_aa-2.0E0*2*I_ESP_G2xyz_Gx3z_a-2.0E0*3*I_ESP_G2xyz_Gx3z_a+2*1*I_ESP_Dyz_Gx3z;
    abcd[iGrid*1350+140] = 4.0E0*I_ESP_I4x2z_Gx3z_aa-2.0E0*2*I_ESP_G2x2z_Gx3z_a-2.0E0*3*I_ESP_G2x2z_Gx3z_a+2*1*I_ESP_D2z_Gx3z;
    abcd[iGrid*1350+141] = 4.0E0*I_ESP_I3x3y_Gx3z_aa-2.0E0*1*I_ESP_Gx3y_Gx3z_a-2.0E0*2*I_ESP_Gx3y_Gx3z_a;
    abcd[iGrid*1350+142] = 4.0E0*I_ESP_I3x2yz_Gx3z_aa-2.0E0*1*I_ESP_Gx2yz_Gx3z_a-2.0E0*2*I_ESP_Gx2yz_Gx3z_a;
    abcd[iGrid*1350+143] = 4.0E0*I_ESP_I3xy2z_Gx3z_aa-2.0E0*1*I_ESP_Gxy2z_Gx3z_a-2.0E0*2*I_ESP_Gxy2z_Gx3z_a;
    abcd[iGrid*1350+144] = 4.0E0*I_ESP_I3x3z_Gx3z_aa-2.0E0*1*I_ESP_Gx3z_Gx3z_a-2.0E0*2*I_ESP_Gx3z_Gx3z_a;
    abcd[iGrid*1350+145] = 4.0E0*I_ESP_I2x4y_Gx3z_aa-2.0E0*1*I_ESP_G4y_Gx3z_a;
    abcd[iGrid*1350+146] = 4.0E0*I_ESP_I2x3yz_Gx3z_aa-2.0E0*1*I_ESP_G3yz_Gx3z_a;
    abcd[iGrid*1350+147] = 4.0E0*I_ESP_I2x2y2z_Gx3z_aa-2.0E0*1*I_ESP_G2y2z_Gx3z_a;
    abcd[iGrid*1350+148] = 4.0E0*I_ESP_I2xy3z_Gx3z_aa-2.0E0*1*I_ESP_Gy3z_Gx3z_a;
    abcd[iGrid*1350+149] = 4.0E0*I_ESP_I2x4z_Gx3z_aa-2.0E0*1*I_ESP_G4z_Gx3z_a;
    abcd[iGrid*1350+150] = 4.0E0*I_ESP_I6x_G4y_aa-2.0E0*4*I_ESP_G4x_G4y_a-2.0E0*5*I_ESP_G4x_G4y_a+4*3*I_ESP_D2x_G4y;
    abcd[iGrid*1350+151] = 4.0E0*I_ESP_I5xy_G4y_aa-2.0E0*3*I_ESP_G3xy_G4y_a-2.0E0*4*I_ESP_G3xy_G4y_a+3*2*I_ESP_Dxy_G4y;
    abcd[iGrid*1350+152] = 4.0E0*I_ESP_I5xz_G4y_aa-2.0E0*3*I_ESP_G3xz_G4y_a-2.0E0*4*I_ESP_G3xz_G4y_a+3*2*I_ESP_Dxz_G4y;
    abcd[iGrid*1350+153] = 4.0E0*I_ESP_I4x2y_G4y_aa-2.0E0*2*I_ESP_G2x2y_G4y_a-2.0E0*3*I_ESP_G2x2y_G4y_a+2*1*I_ESP_D2y_G4y;
    abcd[iGrid*1350+154] = 4.0E0*I_ESP_I4xyz_G4y_aa-2.0E0*2*I_ESP_G2xyz_G4y_a-2.0E0*3*I_ESP_G2xyz_G4y_a+2*1*I_ESP_Dyz_G4y;
    abcd[iGrid*1350+155] = 4.0E0*I_ESP_I4x2z_G4y_aa-2.0E0*2*I_ESP_G2x2z_G4y_a-2.0E0*3*I_ESP_G2x2z_G4y_a+2*1*I_ESP_D2z_G4y;
    abcd[iGrid*1350+156] = 4.0E0*I_ESP_I3x3y_G4y_aa-2.0E0*1*I_ESP_Gx3y_G4y_a-2.0E0*2*I_ESP_Gx3y_G4y_a;
    abcd[iGrid*1350+157] = 4.0E0*I_ESP_I3x2yz_G4y_aa-2.0E0*1*I_ESP_Gx2yz_G4y_a-2.0E0*2*I_ESP_Gx2yz_G4y_a;
    abcd[iGrid*1350+158] = 4.0E0*I_ESP_I3xy2z_G4y_aa-2.0E0*1*I_ESP_Gxy2z_G4y_a-2.0E0*2*I_ESP_Gxy2z_G4y_a;
    abcd[iGrid*1350+159] = 4.0E0*I_ESP_I3x3z_G4y_aa-2.0E0*1*I_ESP_Gx3z_G4y_a-2.0E0*2*I_ESP_Gx3z_G4y_a;
    abcd[iGrid*1350+160] = 4.0E0*I_ESP_I2x4y_G4y_aa-2.0E0*1*I_ESP_G4y_G4y_a;
    abcd[iGrid*1350+161] = 4.0E0*I_ESP_I2x3yz_G4y_aa-2.0E0*1*I_ESP_G3yz_G4y_a;
    abcd[iGrid*1350+162] = 4.0E0*I_ESP_I2x2y2z_G4y_aa-2.0E0*1*I_ESP_G2y2z_G4y_a;
    abcd[iGrid*1350+163] = 4.0E0*I_ESP_I2xy3z_G4y_aa-2.0E0*1*I_ESP_Gy3z_G4y_a;
    abcd[iGrid*1350+164] = 4.0E0*I_ESP_I2x4z_G4y_aa-2.0E0*1*I_ESP_G4z_G4y_a;
    abcd[iGrid*1350+165] = 4.0E0*I_ESP_I6x_G3yz_aa-2.0E0*4*I_ESP_G4x_G3yz_a-2.0E0*5*I_ESP_G4x_G3yz_a+4*3*I_ESP_D2x_G3yz;
    abcd[iGrid*1350+166] = 4.0E0*I_ESP_I5xy_G3yz_aa-2.0E0*3*I_ESP_G3xy_G3yz_a-2.0E0*4*I_ESP_G3xy_G3yz_a+3*2*I_ESP_Dxy_G3yz;
    abcd[iGrid*1350+167] = 4.0E0*I_ESP_I5xz_G3yz_aa-2.0E0*3*I_ESP_G3xz_G3yz_a-2.0E0*4*I_ESP_G3xz_G3yz_a+3*2*I_ESP_Dxz_G3yz;
    abcd[iGrid*1350+168] = 4.0E0*I_ESP_I4x2y_G3yz_aa-2.0E0*2*I_ESP_G2x2y_G3yz_a-2.0E0*3*I_ESP_G2x2y_G3yz_a+2*1*I_ESP_D2y_G3yz;
    abcd[iGrid*1350+169] = 4.0E0*I_ESP_I4xyz_G3yz_aa-2.0E0*2*I_ESP_G2xyz_G3yz_a-2.0E0*3*I_ESP_G2xyz_G3yz_a+2*1*I_ESP_Dyz_G3yz;
    abcd[iGrid*1350+170] = 4.0E0*I_ESP_I4x2z_G3yz_aa-2.0E0*2*I_ESP_G2x2z_G3yz_a-2.0E0*3*I_ESP_G2x2z_G3yz_a+2*1*I_ESP_D2z_G3yz;
    abcd[iGrid*1350+171] = 4.0E0*I_ESP_I3x3y_G3yz_aa-2.0E0*1*I_ESP_Gx3y_G3yz_a-2.0E0*2*I_ESP_Gx3y_G3yz_a;
    abcd[iGrid*1350+172] = 4.0E0*I_ESP_I3x2yz_G3yz_aa-2.0E0*1*I_ESP_Gx2yz_G3yz_a-2.0E0*2*I_ESP_Gx2yz_G3yz_a;
    abcd[iGrid*1350+173] = 4.0E0*I_ESP_I3xy2z_G3yz_aa-2.0E0*1*I_ESP_Gxy2z_G3yz_a-2.0E0*2*I_ESP_Gxy2z_G3yz_a;
    abcd[iGrid*1350+174] = 4.0E0*I_ESP_I3x3z_G3yz_aa-2.0E0*1*I_ESP_Gx3z_G3yz_a-2.0E0*2*I_ESP_Gx3z_G3yz_a;
    abcd[iGrid*1350+175] = 4.0E0*I_ESP_I2x4y_G3yz_aa-2.0E0*1*I_ESP_G4y_G3yz_a;
    abcd[iGrid*1350+176] = 4.0E0*I_ESP_I2x3yz_G3yz_aa-2.0E0*1*I_ESP_G3yz_G3yz_a;
    abcd[iGrid*1350+177] = 4.0E0*I_ESP_I2x2y2z_G3yz_aa-2.0E0*1*I_ESP_G2y2z_G3yz_a;
    abcd[iGrid*1350+178] = 4.0E0*I_ESP_I2xy3z_G3yz_aa-2.0E0*1*I_ESP_Gy3z_G3yz_a;
    abcd[iGrid*1350+179] = 4.0E0*I_ESP_I2x4z_G3yz_aa-2.0E0*1*I_ESP_G4z_G3yz_a;
    abcd[iGrid*1350+180] = 4.0E0*I_ESP_I6x_G2y2z_aa-2.0E0*4*I_ESP_G4x_G2y2z_a-2.0E0*5*I_ESP_G4x_G2y2z_a+4*3*I_ESP_D2x_G2y2z;
    abcd[iGrid*1350+181] = 4.0E0*I_ESP_I5xy_G2y2z_aa-2.0E0*3*I_ESP_G3xy_G2y2z_a-2.0E0*4*I_ESP_G3xy_G2y2z_a+3*2*I_ESP_Dxy_G2y2z;
    abcd[iGrid*1350+182] = 4.0E0*I_ESP_I5xz_G2y2z_aa-2.0E0*3*I_ESP_G3xz_G2y2z_a-2.0E0*4*I_ESP_G3xz_G2y2z_a+3*2*I_ESP_Dxz_G2y2z;
    abcd[iGrid*1350+183] = 4.0E0*I_ESP_I4x2y_G2y2z_aa-2.0E0*2*I_ESP_G2x2y_G2y2z_a-2.0E0*3*I_ESP_G2x2y_G2y2z_a+2*1*I_ESP_D2y_G2y2z;
    abcd[iGrid*1350+184] = 4.0E0*I_ESP_I4xyz_G2y2z_aa-2.0E0*2*I_ESP_G2xyz_G2y2z_a-2.0E0*3*I_ESP_G2xyz_G2y2z_a+2*1*I_ESP_Dyz_G2y2z;
    abcd[iGrid*1350+185] = 4.0E0*I_ESP_I4x2z_G2y2z_aa-2.0E0*2*I_ESP_G2x2z_G2y2z_a-2.0E0*3*I_ESP_G2x2z_G2y2z_a+2*1*I_ESP_D2z_G2y2z;
    abcd[iGrid*1350+186] = 4.0E0*I_ESP_I3x3y_G2y2z_aa-2.0E0*1*I_ESP_Gx3y_G2y2z_a-2.0E0*2*I_ESP_Gx3y_G2y2z_a;
    abcd[iGrid*1350+187] = 4.0E0*I_ESP_I3x2yz_G2y2z_aa-2.0E0*1*I_ESP_Gx2yz_G2y2z_a-2.0E0*2*I_ESP_Gx2yz_G2y2z_a;
    abcd[iGrid*1350+188] = 4.0E0*I_ESP_I3xy2z_G2y2z_aa-2.0E0*1*I_ESP_Gxy2z_G2y2z_a-2.0E0*2*I_ESP_Gxy2z_G2y2z_a;
    abcd[iGrid*1350+189] = 4.0E0*I_ESP_I3x3z_G2y2z_aa-2.0E0*1*I_ESP_Gx3z_G2y2z_a-2.0E0*2*I_ESP_Gx3z_G2y2z_a;
    abcd[iGrid*1350+190] = 4.0E0*I_ESP_I2x4y_G2y2z_aa-2.0E0*1*I_ESP_G4y_G2y2z_a;
    abcd[iGrid*1350+191] = 4.0E0*I_ESP_I2x3yz_G2y2z_aa-2.0E0*1*I_ESP_G3yz_G2y2z_a;
    abcd[iGrid*1350+192] = 4.0E0*I_ESP_I2x2y2z_G2y2z_aa-2.0E0*1*I_ESP_G2y2z_G2y2z_a;
    abcd[iGrid*1350+193] = 4.0E0*I_ESP_I2xy3z_G2y2z_aa-2.0E0*1*I_ESP_Gy3z_G2y2z_a;
    abcd[iGrid*1350+194] = 4.0E0*I_ESP_I2x4z_G2y2z_aa-2.0E0*1*I_ESP_G4z_G2y2z_a;
    abcd[iGrid*1350+195] = 4.0E0*I_ESP_I6x_Gy3z_aa-2.0E0*4*I_ESP_G4x_Gy3z_a-2.0E0*5*I_ESP_G4x_Gy3z_a+4*3*I_ESP_D2x_Gy3z;
    abcd[iGrid*1350+196] = 4.0E0*I_ESP_I5xy_Gy3z_aa-2.0E0*3*I_ESP_G3xy_Gy3z_a-2.0E0*4*I_ESP_G3xy_Gy3z_a+3*2*I_ESP_Dxy_Gy3z;
    abcd[iGrid*1350+197] = 4.0E0*I_ESP_I5xz_Gy3z_aa-2.0E0*3*I_ESP_G3xz_Gy3z_a-2.0E0*4*I_ESP_G3xz_Gy3z_a+3*2*I_ESP_Dxz_Gy3z;
    abcd[iGrid*1350+198] = 4.0E0*I_ESP_I4x2y_Gy3z_aa-2.0E0*2*I_ESP_G2x2y_Gy3z_a-2.0E0*3*I_ESP_G2x2y_Gy3z_a+2*1*I_ESP_D2y_Gy3z;
    abcd[iGrid*1350+199] = 4.0E0*I_ESP_I4xyz_Gy3z_aa-2.0E0*2*I_ESP_G2xyz_Gy3z_a-2.0E0*3*I_ESP_G2xyz_Gy3z_a+2*1*I_ESP_Dyz_Gy3z;
    abcd[iGrid*1350+200] = 4.0E0*I_ESP_I4x2z_Gy3z_aa-2.0E0*2*I_ESP_G2x2z_Gy3z_a-2.0E0*3*I_ESP_G2x2z_Gy3z_a+2*1*I_ESP_D2z_Gy3z;
    abcd[iGrid*1350+201] = 4.0E0*I_ESP_I3x3y_Gy3z_aa-2.0E0*1*I_ESP_Gx3y_Gy3z_a-2.0E0*2*I_ESP_Gx3y_Gy3z_a;
    abcd[iGrid*1350+202] = 4.0E0*I_ESP_I3x2yz_Gy3z_aa-2.0E0*1*I_ESP_Gx2yz_Gy3z_a-2.0E0*2*I_ESP_Gx2yz_Gy3z_a;
    abcd[iGrid*1350+203] = 4.0E0*I_ESP_I3xy2z_Gy3z_aa-2.0E0*1*I_ESP_Gxy2z_Gy3z_a-2.0E0*2*I_ESP_Gxy2z_Gy3z_a;
    abcd[iGrid*1350+204] = 4.0E0*I_ESP_I3x3z_Gy3z_aa-2.0E0*1*I_ESP_Gx3z_Gy3z_a-2.0E0*2*I_ESP_Gx3z_Gy3z_a;
    abcd[iGrid*1350+205] = 4.0E0*I_ESP_I2x4y_Gy3z_aa-2.0E0*1*I_ESP_G4y_Gy3z_a;
    abcd[iGrid*1350+206] = 4.0E0*I_ESP_I2x3yz_Gy3z_aa-2.0E0*1*I_ESP_G3yz_Gy3z_a;
    abcd[iGrid*1350+207] = 4.0E0*I_ESP_I2x2y2z_Gy3z_aa-2.0E0*1*I_ESP_G2y2z_Gy3z_a;
    abcd[iGrid*1350+208] = 4.0E0*I_ESP_I2xy3z_Gy3z_aa-2.0E0*1*I_ESP_Gy3z_Gy3z_a;
    abcd[iGrid*1350+209] = 4.0E0*I_ESP_I2x4z_Gy3z_aa-2.0E0*1*I_ESP_G4z_Gy3z_a;
    abcd[iGrid*1350+210] = 4.0E0*I_ESP_I6x_G4z_aa-2.0E0*4*I_ESP_G4x_G4z_a-2.0E0*5*I_ESP_G4x_G4z_a+4*3*I_ESP_D2x_G4z;
    abcd[iGrid*1350+211] = 4.0E0*I_ESP_I5xy_G4z_aa-2.0E0*3*I_ESP_G3xy_G4z_a-2.0E0*4*I_ESP_G3xy_G4z_a+3*2*I_ESP_Dxy_G4z;
    abcd[iGrid*1350+212] = 4.0E0*I_ESP_I5xz_G4z_aa-2.0E0*3*I_ESP_G3xz_G4z_a-2.0E0*4*I_ESP_G3xz_G4z_a+3*2*I_ESP_Dxz_G4z;
    abcd[iGrid*1350+213] = 4.0E0*I_ESP_I4x2y_G4z_aa-2.0E0*2*I_ESP_G2x2y_G4z_a-2.0E0*3*I_ESP_G2x2y_G4z_a+2*1*I_ESP_D2y_G4z;
    abcd[iGrid*1350+214] = 4.0E0*I_ESP_I4xyz_G4z_aa-2.0E0*2*I_ESP_G2xyz_G4z_a-2.0E0*3*I_ESP_G2xyz_G4z_a+2*1*I_ESP_Dyz_G4z;
    abcd[iGrid*1350+215] = 4.0E0*I_ESP_I4x2z_G4z_aa-2.0E0*2*I_ESP_G2x2z_G4z_a-2.0E0*3*I_ESP_G2x2z_G4z_a+2*1*I_ESP_D2z_G4z;
    abcd[iGrid*1350+216] = 4.0E0*I_ESP_I3x3y_G4z_aa-2.0E0*1*I_ESP_Gx3y_G4z_a-2.0E0*2*I_ESP_Gx3y_G4z_a;
    abcd[iGrid*1350+217] = 4.0E0*I_ESP_I3x2yz_G4z_aa-2.0E0*1*I_ESP_Gx2yz_G4z_a-2.0E0*2*I_ESP_Gx2yz_G4z_a;
    abcd[iGrid*1350+218] = 4.0E0*I_ESP_I3xy2z_G4z_aa-2.0E0*1*I_ESP_Gxy2z_G4z_a-2.0E0*2*I_ESP_Gxy2z_G4z_a;
    abcd[iGrid*1350+219] = 4.0E0*I_ESP_I3x3z_G4z_aa-2.0E0*1*I_ESP_Gx3z_G4z_a-2.0E0*2*I_ESP_Gx3z_G4z_a;
    abcd[iGrid*1350+220] = 4.0E0*I_ESP_I2x4y_G4z_aa-2.0E0*1*I_ESP_G4y_G4z_a;
    abcd[iGrid*1350+221] = 4.0E0*I_ESP_I2x3yz_G4z_aa-2.0E0*1*I_ESP_G3yz_G4z_a;
    abcd[iGrid*1350+222] = 4.0E0*I_ESP_I2x2y2z_G4z_aa-2.0E0*1*I_ESP_G2y2z_G4z_a;
    abcd[iGrid*1350+223] = 4.0E0*I_ESP_I2xy3z_G4z_aa-2.0E0*1*I_ESP_Gy3z_G4z_a;
    abcd[iGrid*1350+224] = 4.0E0*I_ESP_I2x4z_G4z_aa-2.0E0*1*I_ESP_G4z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G_aa
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_D_G
     ************************************************************/
    abcd[iGrid*1350+225] = 4.0E0*I_ESP_I5xy_G4x_aa-2.0E0*4*I_ESP_G3xy_G4x_a;
    abcd[iGrid*1350+226] = 4.0E0*I_ESP_I4x2y_G4x_aa-2.0E0*1*I_ESP_G4x_G4x_a-2.0E0*3*I_ESP_G2x2y_G4x_a+3*1*I_ESP_D2x_G4x;
    abcd[iGrid*1350+227] = 4.0E0*I_ESP_I4xyz_G4x_aa-2.0E0*3*I_ESP_G2xyz_G4x_a;
    abcd[iGrid*1350+228] = 4.0E0*I_ESP_I3x3y_G4x_aa-2.0E0*2*I_ESP_G3xy_G4x_a-2.0E0*2*I_ESP_Gx3y_G4x_a+2*2*I_ESP_Dxy_G4x;
    abcd[iGrid*1350+229] = 4.0E0*I_ESP_I3x2yz_G4x_aa-2.0E0*1*I_ESP_G3xz_G4x_a-2.0E0*2*I_ESP_Gx2yz_G4x_a+2*1*I_ESP_Dxz_G4x;
    abcd[iGrid*1350+230] = 4.0E0*I_ESP_I3xy2z_G4x_aa-2.0E0*2*I_ESP_Gxy2z_G4x_a;
    abcd[iGrid*1350+231] = 4.0E0*I_ESP_I2x4y_G4x_aa-2.0E0*3*I_ESP_G2x2y_G4x_a-2.0E0*1*I_ESP_G4y_G4x_a+3*I_ESP_D2y_G4x;
    abcd[iGrid*1350+232] = 4.0E0*I_ESP_I2x3yz_G4x_aa-2.0E0*2*I_ESP_G2xyz_G4x_a-2.0E0*1*I_ESP_G3yz_G4x_a+2*I_ESP_Dyz_G4x;
    abcd[iGrid*1350+233] = 4.0E0*I_ESP_I2x2y2z_G4x_aa-2.0E0*1*I_ESP_G2x2z_G4x_a-2.0E0*1*I_ESP_G2y2z_G4x_a+1*I_ESP_D2z_G4x;
    abcd[iGrid*1350+234] = 4.0E0*I_ESP_I2xy3z_G4x_aa-2.0E0*1*I_ESP_Gy3z_G4x_a;
    abcd[iGrid*1350+235] = 4.0E0*I_ESP_Ix5y_G4x_aa-2.0E0*4*I_ESP_Gx3y_G4x_a;
    abcd[iGrid*1350+236] = 4.0E0*I_ESP_Ix4yz_G4x_aa-2.0E0*3*I_ESP_Gx2yz_G4x_a;
    abcd[iGrid*1350+237] = 4.0E0*I_ESP_Ix3y2z_G4x_aa-2.0E0*2*I_ESP_Gxy2z_G4x_a;
    abcd[iGrid*1350+238] = 4.0E0*I_ESP_Ix2y3z_G4x_aa-2.0E0*1*I_ESP_Gx3z_G4x_a;
    abcd[iGrid*1350+239] = 4.0E0*I_ESP_Ixy4z_G4x_aa;
    abcd[iGrid*1350+240] = 4.0E0*I_ESP_I5xy_G3xy_aa-2.0E0*4*I_ESP_G3xy_G3xy_a;
    abcd[iGrid*1350+241] = 4.0E0*I_ESP_I4x2y_G3xy_aa-2.0E0*1*I_ESP_G4x_G3xy_a-2.0E0*3*I_ESP_G2x2y_G3xy_a+3*1*I_ESP_D2x_G3xy;
    abcd[iGrid*1350+242] = 4.0E0*I_ESP_I4xyz_G3xy_aa-2.0E0*3*I_ESP_G2xyz_G3xy_a;
    abcd[iGrid*1350+243] = 4.0E0*I_ESP_I3x3y_G3xy_aa-2.0E0*2*I_ESP_G3xy_G3xy_a-2.0E0*2*I_ESP_Gx3y_G3xy_a+2*2*I_ESP_Dxy_G3xy;
    abcd[iGrid*1350+244] = 4.0E0*I_ESP_I3x2yz_G3xy_aa-2.0E0*1*I_ESP_G3xz_G3xy_a-2.0E0*2*I_ESP_Gx2yz_G3xy_a+2*1*I_ESP_Dxz_G3xy;
    abcd[iGrid*1350+245] = 4.0E0*I_ESP_I3xy2z_G3xy_aa-2.0E0*2*I_ESP_Gxy2z_G3xy_a;
    abcd[iGrid*1350+246] = 4.0E0*I_ESP_I2x4y_G3xy_aa-2.0E0*3*I_ESP_G2x2y_G3xy_a-2.0E0*1*I_ESP_G4y_G3xy_a+3*I_ESP_D2y_G3xy;
    abcd[iGrid*1350+247] = 4.0E0*I_ESP_I2x3yz_G3xy_aa-2.0E0*2*I_ESP_G2xyz_G3xy_a-2.0E0*1*I_ESP_G3yz_G3xy_a+2*I_ESP_Dyz_G3xy;
    abcd[iGrid*1350+248] = 4.0E0*I_ESP_I2x2y2z_G3xy_aa-2.0E0*1*I_ESP_G2x2z_G3xy_a-2.0E0*1*I_ESP_G2y2z_G3xy_a+1*I_ESP_D2z_G3xy;
    abcd[iGrid*1350+249] = 4.0E0*I_ESP_I2xy3z_G3xy_aa-2.0E0*1*I_ESP_Gy3z_G3xy_a;
    abcd[iGrid*1350+250] = 4.0E0*I_ESP_Ix5y_G3xy_aa-2.0E0*4*I_ESP_Gx3y_G3xy_a;
    abcd[iGrid*1350+251] = 4.0E0*I_ESP_Ix4yz_G3xy_aa-2.0E0*3*I_ESP_Gx2yz_G3xy_a;
    abcd[iGrid*1350+252] = 4.0E0*I_ESP_Ix3y2z_G3xy_aa-2.0E0*2*I_ESP_Gxy2z_G3xy_a;
    abcd[iGrid*1350+253] = 4.0E0*I_ESP_Ix2y3z_G3xy_aa-2.0E0*1*I_ESP_Gx3z_G3xy_a;
    abcd[iGrid*1350+254] = 4.0E0*I_ESP_Ixy4z_G3xy_aa;
    abcd[iGrid*1350+255] = 4.0E0*I_ESP_I5xy_G3xz_aa-2.0E0*4*I_ESP_G3xy_G3xz_a;
    abcd[iGrid*1350+256] = 4.0E0*I_ESP_I4x2y_G3xz_aa-2.0E0*1*I_ESP_G4x_G3xz_a-2.0E0*3*I_ESP_G2x2y_G3xz_a+3*1*I_ESP_D2x_G3xz;
    abcd[iGrid*1350+257] = 4.0E0*I_ESP_I4xyz_G3xz_aa-2.0E0*3*I_ESP_G2xyz_G3xz_a;
    abcd[iGrid*1350+258] = 4.0E0*I_ESP_I3x3y_G3xz_aa-2.0E0*2*I_ESP_G3xy_G3xz_a-2.0E0*2*I_ESP_Gx3y_G3xz_a+2*2*I_ESP_Dxy_G3xz;
    abcd[iGrid*1350+259] = 4.0E0*I_ESP_I3x2yz_G3xz_aa-2.0E0*1*I_ESP_G3xz_G3xz_a-2.0E0*2*I_ESP_Gx2yz_G3xz_a+2*1*I_ESP_Dxz_G3xz;
    abcd[iGrid*1350+260] = 4.0E0*I_ESP_I3xy2z_G3xz_aa-2.0E0*2*I_ESP_Gxy2z_G3xz_a;
    abcd[iGrid*1350+261] = 4.0E0*I_ESP_I2x4y_G3xz_aa-2.0E0*3*I_ESP_G2x2y_G3xz_a-2.0E0*1*I_ESP_G4y_G3xz_a+3*I_ESP_D2y_G3xz;
    abcd[iGrid*1350+262] = 4.0E0*I_ESP_I2x3yz_G3xz_aa-2.0E0*2*I_ESP_G2xyz_G3xz_a-2.0E0*1*I_ESP_G3yz_G3xz_a+2*I_ESP_Dyz_G3xz;
    abcd[iGrid*1350+263] = 4.0E0*I_ESP_I2x2y2z_G3xz_aa-2.0E0*1*I_ESP_G2x2z_G3xz_a-2.0E0*1*I_ESP_G2y2z_G3xz_a+1*I_ESP_D2z_G3xz;
    abcd[iGrid*1350+264] = 4.0E0*I_ESP_I2xy3z_G3xz_aa-2.0E0*1*I_ESP_Gy3z_G3xz_a;
    abcd[iGrid*1350+265] = 4.0E0*I_ESP_Ix5y_G3xz_aa-2.0E0*4*I_ESP_Gx3y_G3xz_a;
    abcd[iGrid*1350+266] = 4.0E0*I_ESP_Ix4yz_G3xz_aa-2.0E0*3*I_ESP_Gx2yz_G3xz_a;
    abcd[iGrid*1350+267] = 4.0E0*I_ESP_Ix3y2z_G3xz_aa-2.0E0*2*I_ESP_Gxy2z_G3xz_a;
    abcd[iGrid*1350+268] = 4.0E0*I_ESP_Ix2y3z_G3xz_aa-2.0E0*1*I_ESP_Gx3z_G3xz_a;
    abcd[iGrid*1350+269] = 4.0E0*I_ESP_Ixy4z_G3xz_aa;
    abcd[iGrid*1350+270] = 4.0E0*I_ESP_I5xy_G2x2y_aa-2.0E0*4*I_ESP_G3xy_G2x2y_a;
    abcd[iGrid*1350+271] = 4.0E0*I_ESP_I4x2y_G2x2y_aa-2.0E0*1*I_ESP_G4x_G2x2y_a-2.0E0*3*I_ESP_G2x2y_G2x2y_a+3*1*I_ESP_D2x_G2x2y;
    abcd[iGrid*1350+272] = 4.0E0*I_ESP_I4xyz_G2x2y_aa-2.0E0*3*I_ESP_G2xyz_G2x2y_a;
    abcd[iGrid*1350+273] = 4.0E0*I_ESP_I3x3y_G2x2y_aa-2.0E0*2*I_ESP_G3xy_G2x2y_a-2.0E0*2*I_ESP_Gx3y_G2x2y_a+2*2*I_ESP_Dxy_G2x2y;
    abcd[iGrid*1350+274] = 4.0E0*I_ESP_I3x2yz_G2x2y_aa-2.0E0*1*I_ESP_G3xz_G2x2y_a-2.0E0*2*I_ESP_Gx2yz_G2x2y_a+2*1*I_ESP_Dxz_G2x2y;
    abcd[iGrid*1350+275] = 4.0E0*I_ESP_I3xy2z_G2x2y_aa-2.0E0*2*I_ESP_Gxy2z_G2x2y_a;
    abcd[iGrid*1350+276] = 4.0E0*I_ESP_I2x4y_G2x2y_aa-2.0E0*3*I_ESP_G2x2y_G2x2y_a-2.0E0*1*I_ESP_G4y_G2x2y_a+3*I_ESP_D2y_G2x2y;
    abcd[iGrid*1350+277] = 4.0E0*I_ESP_I2x3yz_G2x2y_aa-2.0E0*2*I_ESP_G2xyz_G2x2y_a-2.0E0*1*I_ESP_G3yz_G2x2y_a+2*I_ESP_Dyz_G2x2y;
    abcd[iGrid*1350+278] = 4.0E0*I_ESP_I2x2y2z_G2x2y_aa-2.0E0*1*I_ESP_G2x2z_G2x2y_a-2.0E0*1*I_ESP_G2y2z_G2x2y_a+1*I_ESP_D2z_G2x2y;
    abcd[iGrid*1350+279] = 4.0E0*I_ESP_I2xy3z_G2x2y_aa-2.0E0*1*I_ESP_Gy3z_G2x2y_a;
    abcd[iGrid*1350+280] = 4.0E0*I_ESP_Ix5y_G2x2y_aa-2.0E0*4*I_ESP_Gx3y_G2x2y_a;
    abcd[iGrid*1350+281] = 4.0E0*I_ESP_Ix4yz_G2x2y_aa-2.0E0*3*I_ESP_Gx2yz_G2x2y_a;
    abcd[iGrid*1350+282] = 4.0E0*I_ESP_Ix3y2z_G2x2y_aa-2.0E0*2*I_ESP_Gxy2z_G2x2y_a;
    abcd[iGrid*1350+283] = 4.0E0*I_ESP_Ix2y3z_G2x2y_aa-2.0E0*1*I_ESP_Gx3z_G2x2y_a;
    abcd[iGrid*1350+284] = 4.0E0*I_ESP_Ixy4z_G2x2y_aa;
    abcd[iGrid*1350+285] = 4.0E0*I_ESP_I5xy_G2xyz_aa-2.0E0*4*I_ESP_G3xy_G2xyz_a;
    abcd[iGrid*1350+286] = 4.0E0*I_ESP_I4x2y_G2xyz_aa-2.0E0*1*I_ESP_G4x_G2xyz_a-2.0E0*3*I_ESP_G2x2y_G2xyz_a+3*1*I_ESP_D2x_G2xyz;
    abcd[iGrid*1350+287] = 4.0E0*I_ESP_I4xyz_G2xyz_aa-2.0E0*3*I_ESP_G2xyz_G2xyz_a;
    abcd[iGrid*1350+288] = 4.0E0*I_ESP_I3x3y_G2xyz_aa-2.0E0*2*I_ESP_G3xy_G2xyz_a-2.0E0*2*I_ESP_Gx3y_G2xyz_a+2*2*I_ESP_Dxy_G2xyz;
    abcd[iGrid*1350+289] = 4.0E0*I_ESP_I3x2yz_G2xyz_aa-2.0E0*1*I_ESP_G3xz_G2xyz_a-2.0E0*2*I_ESP_Gx2yz_G2xyz_a+2*1*I_ESP_Dxz_G2xyz;
    abcd[iGrid*1350+290] = 4.0E0*I_ESP_I3xy2z_G2xyz_aa-2.0E0*2*I_ESP_Gxy2z_G2xyz_a;
    abcd[iGrid*1350+291] = 4.0E0*I_ESP_I2x4y_G2xyz_aa-2.0E0*3*I_ESP_G2x2y_G2xyz_a-2.0E0*1*I_ESP_G4y_G2xyz_a+3*I_ESP_D2y_G2xyz;
    abcd[iGrid*1350+292] = 4.0E0*I_ESP_I2x3yz_G2xyz_aa-2.0E0*2*I_ESP_G2xyz_G2xyz_a-2.0E0*1*I_ESP_G3yz_G2xyz_a+2*I_ESP_Dyz_G2xyz;
    abcd[iGrid*1350+293] = 4.0E0*I_ESP_I2x2y2z_G2xyz_aa-2.0E0*1*I_ESP_G2x2z_G2xyz_a-2.0E0*1*I_ESP_G2y2z_G2xyz_a+1*I_ESP_D2z_G2xyz;
    abcd[iGrid*1350+294] = 4.0E0*I_ESP_I2xy3z_G2xyz_aa-2.0E0*1*I_ESP_Gy3z_G2xyz_a;
    abcd[iGrid*1350+295] = 4.0E0*I_ESP_Ix5y_G2xyz_aa-2.0E0*4*I_ESP_Gx3y_G2xyz_a;
    abcd[iGrid*1350+296] = 4.0E0*I_ESP_Ix4yz_G2xyz_aa-2.0E0*3*I_ESP_Gx2yz_G2xyz_a;
    abcd[iGrid*1350+297] = 4.0E0*I_ESP_Ix3y2z_G2xyz_aa-2.0E0*2*I_ESP_Gxy2z_G2xyz_a;
    abcd[iGrid*1350+298] = 4.0E0*I_ESP_Ix2y3z_G2xyz_aa-2.0E0*1*I_ESP_Gx3z_G2xyz_a;
    abcd[iGrid*1350+299] = 4.0E0*I_ESP_Ixy4z_G2xyz_aa;
    abcd[iGrid*1350+300] = 4.0E0*I_ESP_I5xy_G2x2z_aa-2.0E0*4*I_ESP_G3xy_G2x2z_a;
    abcd[iGrid*1350+301] = 4.0E0*I_ESP_I4x2y_G2x2z_aa-2.0E0*1*I_ESP_G4x_G2x2z_a-2.0E0*3*I_ESP_G2x2y_G2x2z_a+3*1*I_ESP_D2x_G2x2z;
    abcd[iGrid*1350+302] = 4.0E0*I_ESP_I4xyz_G2x2z_aa-2.0E0*3*I_ESP_G2xyz_G2x2z_a;
    abcd[iGrid*1350+303] = 4.0E0*I_ESP_I3x3y_G2x2z_aa-2.0E0*2*I_ESP_G3xy_G2x2z_a-2.0E0*2*I_ESP_Gx3y_G2x2z_a+2*2*I_ESP_Dxy_G2x2z;
    abcd[iGrid*1350+304] = 4.0E0*I_ESP_I3x2yz_G2x2z_aa-2.0E0*1*I_ESP_G3xz_G2x2z_a-2.0E0*2*I_ESP_Gx2yz_G2x2z_a+2*1*I_ESP_Dxz_G2x2z;
    abcd[iGrid*1350+305] = 4.0E0*I_ESP_I3xy2z_G2x2z_aa-2.0E0*2*I_ESP_Gxy2z_G2x2z_a;
    abcd[iGrid*1350+306] = 4.0E0*I_ESP_I2x4y_G2x2z_aa-2.0E0*3*I_ESP_G2x2y_G2x2z_a-2.0E0*1*I_ESP_G4y_G2x2z_a+3*I_ESP_D2y_G2x2z;
    abcd[iGrid*1350+307] = 4.0E0*I_ESP_I2x3yz_G2x2z_aa-2.0E0*2*I_ESP_G2xyz_G2x2z_a-2.0E0*1*I_ESP_G3yz_G2x2z_a+2*I_ESP_Dyz_G2x2z;
    abcd[iGrid*1350+308] = 4.0E0*I_ESP_I2x2y2z_G2x2z_aa-2.0E0*1*I_ESP_G2x2z_G2x2z_a-2.0E0*1*I_ESP_G2y2z_G2x2z_a+1*I_ESP_D2z_G2x2z;
    abcd[iGrid*1350+309] = 4.0E0*I_ESP_I2xy3z_G2x2z_aa-2.0E0*1*I_ESP_Gy3z_G2x2z_a;
    abcd[iGrid*1350+310] = 4.0E0*I_ESP_Ix5y_G2x2z_aa-2.0E0*4*I_ESP_Gx3y_G2x2z_a;
    abcd[iGrid*1350+311] = 4.0E0*I_ESP_Ix4yz_G2x2z_aa-2.0E0*3*I_ESP_Gx2yz_G2x2z_a;
    abcd[iGrid*1350+312] = 4.0E0*I_ESP_Ix3y2z_G2x2z_aa-2.0E0*2*I_ESP_Gxy2z_G2x2z_a;
    abcd[iGrid*1350+313] = 4.0E0*I_ESP_Ix2y3z_G2x2z_aa-2.0E0*1*I_ESP_Gx3z_G2x2z_a;
    abcd[iGrid*1350+314] = 4.0E0*I_ESP_Ixy4z_G2x2z_aa;
    abcd[iGrid*1350+315] = 4.0E0*I_ESP_I5xy_Gx3y_aa-2.0E0*4*I_ESP_G3xy_Gx3y_a;
    abcd[iGrid*1350+316] = 4.0E0*I_ESP_I4x2y_Gx3y_aa-2.0E0*1*I_ESP_G4x_Gx3y_a-2.0E0*3*I_ESP_G2x2y_Gx3y_a+3*1*I_ESP_D2x_Gx3y;
    abcd[iGrid*1350+317] = 4.0E0*I_ESP_I4xyz_Gx3y_aa-2.0E0*3*I_ESP_G2xyz_Gx3y_a;
    abcd[iGrid*1350+318] = 4.0E0*I_ESP_I3x3y_Gx3y_aa-2.0E0*2*I_ESP_G3xy_Gx3y_a-2.0E0*2*I_ESP_Gx3y_Gx3y_a+2*2*I_ESP_Dxy_Gx3y;
    abcd[iGrid*1350+319] = 4.0E0*I_ESP_I3x2yz_Gx3y_aa-2.0E0*1*I_ESP_G3xz_Gx3y_a-2.0E0*2*I_ESP_Gx2yz_Gx3y_a+2*1*I_ESP_Dxz_Gx3y;
    abcd[iGrid*1350+320] = 4.0E0*I_ESP_I3xy2z_Gx3y_aa-2.0E0*2*I_ESP_Gxy2z_Gx3y_a;
    abcd[iGrid*1350+321] = 4.0E0*I_ESP_I2x4y_Gx3y_aa-2.0E0*3*I_ESP_G2x2y_Gx3y_a-2.0E0*1*I_ESP_G4y_Gx3y_a+3*I_ESP_D2y_Gx3y;
    abcd[iGrid*1350+322] = 4.0E0*I_ESP_I2x3yz_Gx3y_aa-2.0E0*2*I_ESP_G2xyz_Gx3y_a-2.0E0*1*I_ESP_G3yz_Gx3y_a+2*I_ESP_Dyz_Gx3y;
    abcd[iGrid*1350+323] = 4.0E0*I_ESP_I2x2y2z_Gx3y_aa-2.0E0*1*I_ESP_G2x2z_Gx3y_a-2.0E0*1*I_ESP_G2y2z_Gx3y_a+1*I_ESP_D2z_Gx3y;
    abcd[iGrid*1350+324] = 4.0E0*I_ESP_I2xy3z_Gx3y_aa-2.0E0*1*I_ESP_Gy3z_Gx3y_a;
    abcd[iGrid*1350+325] = 4.0E0*I_ESP_Ix5y_Gx3y_aa-2.0E0*4*I_ESP_Gx3y_Gx3y_a;
    abcd[iGrid*1350+326] = 4.0E0*I_ESP_Ix4yz_Gx3y_aa-2.0E0*3*I_ESP_Gx2yz_Gx3y_a;
    abcd[iGrid*1350+327] = 4.0E0*I_ESP_Ix3y2z_Gx3y_aa-2.0E0*2*I_ESP_Gxy2z_Gx3y_a;
    abcd[iGrid*1350+328] = 4.0E0*I_ESP_Ix2y3z_Gx3y_aa-2.0E0*1*I_ESP_Gx3z_Gx3y_a;
    abcd[iGrid*1350+329] = 4.0E0*I_ESP_Ixy4z_Gx3y_aa;
    abcd[iGrid*1350+330] = 4.0E0*I_ESP_I5xy_Gx2yz_aa-2.0E0*4*I_ESP_G3xy_Gx2yz_a;
    abcd[iGrid*1350+331] = 4.0E0*I_ESP_I4x2y_Gx2yz_aa-2.0E0*1*I_ESP_G4x_Gx2yz_a-2.0E0*3*I_ESP_G2x2y_Gx2yz_a+3*1*I_ESP_D2x_Gx2yz;
    abcd[iGrid*1350+332] = 4.0E0*I_ESP_I4xyz_Gx2yz_aa-2.0E0*3*I_ESP_G2xyz_Gx2yz_a;
    abcd[iGrid*1350+333] = 4.0E0*I_ESP_I3x3y_Gx2yz_aa-2.0E0*2*I_ESP_G3xy_Gx2yz_a-2.0E0*2*I_ESP_Gx3y_Gx2yz_a+2*2*I_ESP_Dxy_Gx2yz;
    abcd[iGrid*1350+334] = 4.0E0*I_ESP_I3x2yz_Gx2yz_aa-2.0E0*1*I_ESP_G3xz_Gx2yz_a-2.0E0*2*I_ESP_Gx2yz_Gx2yz_a+2*1*I_ESP_Dxz_Gx2yz;
    abcd[iGrid*1350+335] = 4.0E0*I_ESP_I3xy2z_Gx2yz_aa-2.0E0*2*I_ESP_Gxy2z_Gx2yz_a;
    abcd[iGrid*1350+336] = 4.0E0*I_ESP_I2x4y_Gx2yz_aa-2.0E0*3*I_ESP_G2x2y_Gx2yz_a-2.0E0*1*I_ESP_G4y_Gx2yz_a+3*I_ESP_D2y_Gx2yz;
    abcd[iGrid*1350+337] = 4.0E0*I_ESP_I2x3yz_Gx2yz_aa-2.0E0*2*I_ESP_G2xyz_Gx2yz_a-2.0E0*1*I_ESP_G3yz_Gx2yz_a+2*I_ESP_Dyz_Gx2yz;
    abcd[iGrid*1350+338] = 4.0E0*I_ESP_I2x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_G2x2z_Gx2yz_a-2.0E0*1*I_ESP_G2y2z_Gx2yz_a+1*I_ESP_D2z_Gx2yz;
    abcd[iGrid*1350+339] = 4.0E0*I_ESP_I2xy3z_Gx2yz_aa-2.0E0*1*I_ESP_Gy3z_Gx2yz_a;
    abcd[iGrid*1350+340] = 4.0E0*I_ESP_Ix5y_Gx2yz_aa-2.0E0*4*I_ESP_Gx3y_Gx2yz_a;
    abcd[iGrid*1350+341] = 4.0E0*I_ESP_Ix4yz_Gx2yz_aa-2.0E0*3*I_ESP_Gx2yz_Gx2yz_a;
    abcd[iGrid*1350+342] = 4.0E0*I_ESP_Ix3y2z_Gx2yz_aa-2.0E0*2*I_ESP_Gxy2z_Gx2yz_a;
    abcd[iGrid*1350+343] = 4.0E0*I_ESP_Ix2y3z_Gx2yz_aa-2.0E0*1*I_ESP_Gx3z_Gx2yz_a;
    abcd[iGrid*1350+344] = 4.0E0*I_ESP_Ixy4z_Gx2yz_aa;
    abcd[iGrid*1350+345] = 4.0E0*I_ESP_I5xy_Gxy2z_aa-2.0E0*4*I_ESP_G3xy_Gxy2z_a;
    abcd[iGrid*1350+346] = 4.0E0*I_ESP_I4x2y_Gxy2z_aa-2.0E0*1*I_ESP_G4x_Gxy2z_a-2.0E0*3*I_ESP_G2x2y_Gxy2z_a+3*1*I_ESP_D2x_Gxy2z;
    abcd[iGrid*1350+347] = 4.0E0*I_ESP_I4xyz_Gxy2z_aa-2.0E0*3*I_ESP_G2xyz_Gxy2z_a;
    abcd[iGrid*1350+348] = 4.0E0*I_ESP_I3x3y_Gxy2z_aa-2.0E0*2*I_ESP_G3xy_Gxy2z_a-2.0E0*2*I_ESP_Gx3y_Gxy2z_a+2*2*I_ESP_Dxy_Gxy2z;
    abcd[iGrid*1350+349] = 4.0E0*I_ESP_I3x2yz_Gxy2z_aa-2.0E0*1*I_ESP_G3xz_Gxy2z_a-2.0E0*2*I_ESP_Gx2yz_Gxy2z_a+2*1*I_ESP_Dxz_Gxy2z;
    abcd[iGrid*1350+350] = 4.0E0*I_ESP_I3xy2z_Gxy2z_aa-2.0E0*2*I_ESP_Gxy2z_Gxy2z_a;
    abcd[iGrid*1350+351] = 4.0E0*I_ESP_I2x4y_Gxy2z_aa-2.0E0*3*I_ESP_G2x2y_Gxy2z_a-2.0E0*1*I_ESP_G4y_Gxy2z_a+3*I_ESP_D2y_Gxy2z;
    abcd[iGrid*1350+352] = 4.0E0*I_ESP_I2x3yz_Gxy2z_aa-2.0E0*2*I_ESP_G2xyz_Gxy2z_a-2.0E0*1*I_ESP_G3yz_Gxy2z_a+2*I_ESP_Dyz_Gxy2z;
    abcd[iGrid*1350+353] = 4.0E0*I_ESP_I2x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_G2x2z_Gxy2z_a-2.0E0*1*I_ESP_G2y2z_Gxy2z_a+1*I_ESP_D2z_Gxy2z;
    abcd[iGrid*1350+354] = 4.0E0*I_ESP_I2xy3z_Gxy2z_aa-2.0E0*1*I_ESP_Gy3z_Gxy2z_a;
    abcd[iGrid*1350+355] = 4.0E0*I_ESP_Ix5y_Gxy2z_aa-2.0E0*4*I_ESP_Gx3y_Gxy2z_a;
    abcd[iGrid*1350+356] = 4.0E0*I_ESP_Ix4yz_Gxy2z_aa-2.0E0*3*I_ESP_Gx2yz_Gxy2z_a;
    abcd[iGrid*1350+357] = 4.0E0*I_ESP_Ix3y2z_Gxy2z_aa-2.0E0*2*I_ESP_Gxy2z_Gxy2z_a;
    abcd[iGrid*1350+358] = 4.0E0*I_ESP_Ix2y3z_Gxy2z_aa-2.0E0*1*I_ESP_Gx3z_Gxy2z_a;
    abcd[iGrid*1350+359] = 4.0E0*I_ESP_Ixy4z_Gxy2z_aa;
    abcd[iGrid*1350+360] = 4.0E0*I_ESP_I5xy_Gx3z_aa-2.0E0*4*I_ESP_G3xy_Gx3z_a;
    abcd[iGrid*1350+361] = 4.0E0*I_ESP_I4x2y_Gx3z_aa-2.0E0*1*I_ESP_G4x_Gx3z_a-2.0E0*3*I_ESP_G2x2y_Gx3z_a+3*1*I_ESP_D2x_Gx3z;
    abcd[iGrid*1350+362] = 4.0E0*I_ESP_I4xyz_Gx3z_aa-2.0E0*3*I_ESP_G2xyz_Gx3z_a;
    abcd[iGrid*1350+363] = 4.0E0*I_ESP_I3x3y_Gx3z_aa-2.0E0*2*I_ESP_G3xy_Gx3z_a-2.0E0*2*I_ESP_Gx3y_Gx3z_a+2*2*I_ESP_Dxy_Gx3z;
    abcd[iGrid*1350+364] = 4.0E0*I_ESP_I3x2yz_Gx3z_aa-2.0E0*1*I_ESP_G3xz_Gx3z_a-2.0E0*2*I_ESP_Gx2yz_Gx3z_a+2*1*I_ESP_Dxz_Gx3z;
    abcd[iGrid*1350+365] = 4.0E0*I_ESP_I3xy2z_Gx3z_aa-2.0E0*2*I_ESP_Gxy2z_Gx3z_a;
    abcd[iGrid*1350+366] = 4.0E0*I_ESP_I2x4y_Gx3z_aa-2.0E0*3*I_ESP_G2x2y_Gx3z_a-2.0E0*1*I_ESP_G4y_Gx3z_a+3*I_ESP_D2y_Gx3z;
    abcd[iGrid*1350+367] = 4.0E0*I_ESP_I2x3yz_Gx3z_aa-2.0E0*2*I_ESP_G2xyz_Gx3z_a-2.0E0*1*I_ESP_G3yz_Gx3z_a+2*I_ESP_Dyz_Gx3z;
    abcd[iGrid*1350+368] = 4.0E0*I_ESP_I2x2y2z_Gx3z_aa-2.0E0*1*I_ESP_G2x2z_Gx3z_a-2.0E0*1*I_ESP_G2y2z_Gx3z_a+1*I_ESP_D2z_Gx3z;
    abcd[iGrid*1350+369] = 4.0E0*I_ESP_I2xy3z_Gx3z_aa-2.0E0*1*I_ESP_Gy3z_Gx3z_a;
    abcd[iGrid*1350+370] = 4.0E0*I_ESP_Ix5y_Gx3z_aa-2.0E0*4*I_ESP_Gx3y_Gx3z_a;
    abcd[iGrid*1350+371] = 4.0E0*I_ESP_Ix4yz_Gx3z_aa-2.0E0*3*I_ESP_Gx2yz_Gx3z_a;
    abcd[iGrid*1350+372] = 4.0E0*I_ESP_Ix3y2z_Gx3z_aa-2.0E0*2*I_ESP_Gxy2z_Gx3z_a;
    abcd[iGrid*1350+373] = 4.0E0*I_ESP_Ix2y3z_Gx3z_aa-2.0E0*1*I_ESP_Gx3z_Gx3z_a;
    abcd[iGrid*1350+374] = 4.0E0*I_ESP_Ixy4z_Gx3z_aa;
    abcd[iGrid*1350+375] = 4.0E0*I_ESP_I5xy_G4y_aa-2.0E0*4*I_ESP_G3xy_G4y_a;
    abcd[iGrid*1350+376] = 4.0E0*I_ESP_I4x2y_G4y_aa-2.0E0*1*I_ESP_G4x_G4y_a-2.0E0*3*I_ESP_G2x2y_G4y_a+3*1*I_ESP_D2x_G4y;
    abcd[iGrid*1350+377] = 4.0E0*I_ESP_I4xyz_G4y_aa-2.0E0*3*I_ESP_G2xyz_G4y_a;
    abcd[iGrid*1350+378] = 4.0E0*I_ESP_I3x3y_G4y_aa-2.0E0*2*I_ESP_G3xy_G4y_a-2.0E0*2*I_ESP_Gx3y_G4y_a+2*2*I_ESP_Dxy_G4y;
    abcd[iGrid*1350+379] = 4.0E0*I_ESP_I3x2yz_G4y_aa-2.0E0*1*I_ESP_G3xz_G4y_a-2.0E0*2*I_ESP_Gx2yz_G4y_a+2*1*I_ESP_Dxz_G4y;
    abcd[iGrid*1350+380] = 4.0E0*I_ESP_I3xy2z_G4y_aa-2.0E0*2*I_ESP_Gxy2z_G4y_a;
    abcd[iGrid*1350+381] = 4.0E0*I_ESP_I2x4y_G4y_aa-2.0E0*3*I_ESP_G2x2y_G4y_a-2.0E0*1*I_ESP_G4y_G4y_a+3*I_ESP_D2y_G4y;
    abcd[iGrid*1350+382] = 4.0E0*I_ESP_I2x3yz_G4y_aa-2.0E0*2*I_ESP_G2xyz_G4y_a-2.0E0*1*I_ESP_G3yz_G4y_a+2*I_ESP_Dyz_G4y;
    abcd[iGrid*1350+383] = 4.0E0*I_ESP_I2x2y2z_G4y_aa-2.0E0*1*I_ESP_G2x2z_G4y_a-2.0E0*1*I_ESP_G2y2z_G4y_a+1*I_ESP_D2z_G4y;
    abcd[iGrid*1350+384] = 4.0E0*I_ESP_I2xy3z_G4y_aa-2.0E0*1*I_ESP_Gy3z_G4y_a;
    abcd[iGrid*1350+385] = 4.0E0*I_ESP_Ix5y_G4y_aa-2.0E0*4*I_ESP_Gx3y_G4y_a;
    abcd[iGrid*1350+386] = 4.0E0*I_ESP_Ix4yz_G4y_aa-2.0E0*3*I_ESP_Gx2yz_G4y_a;
    abcd[iGrid*1350+387] = 4.0E0*I_ESP_Ix3y2z_G4y_aa-2.0E0*2*I_ESP_Gxy2z_G4y_a;
    abcd[iGrid*1350+388] = 4.0E0*I_ESP_Ix2y3z_G4y_aa-2.0E0*1*I_ESP_Gx3z_G4y_a;
    abcd[iGrid*1350+389] = 4.0E0*I_ESP_Ixy4z_G4y_aa;
    abcd[iGrid*1350+390] = 4.0E0*I_ESP_I5xy_G3yz_aa-2.0E0*4*I_ESP_G3xy_G3yz_a;
    abcd[iGrid*1350+391] = 4.0E0*I_ESP_I4x2y_G3yz_aa-2.0E0*1*I_ESP_G4x_G3yz_a-2.0E0*3*I_ESP_G2x2y_G3yz_a+3*1*I_ESP_D2x_G3yz;
    abcd[iGrid*1350+392] = 4.0E0*I_ESP_I4xyz_G3yz_aa-2.0E0*3*I_ESP_G2xyz_G3yz_a;
    abcd[iGrid*1350+393] = 4.0E0*I_ESP_I3x3y_G3yz_aa-2.0E0*2*I_ESP_G3xy_G3yz_a-2.0E0*2*I_ESP_Gx3y_G3yz_a+2*2*I_ESP_Dxy_G3yz;
    abcd[iGrid*1350+394] = 4.0E0*I_ESP_I3x2yz_G3yz_aa-2.0E0*1*I_ESP_G3xz_G3yz_a-2.0E0*2*I_ESP_Gx2yz_G3yz_a+2*1*I_ESP_Dxz_G3yz;
    abcd[iGrid*1350+395] = 4.0E0*I_ESP_I3xy2z_G3yz_aa-2.0E0*2*I_ESP_Gxy2z_G3yz_a;
    abcd[iGrid*1350+396] = 4.0E0*I_ESP_I2x4y_G3yz_aa-2.0E0*3*I_ESP_G2x2y_G3yz_a-2.0E0*1*I_ESP_G4y_G3yz_a+3*I_ESP_D2y_G3yz;
    abcd[iGrid*1350+397] = 4.0E0*I_ESP_I2x3yz_G3yz_aa-2.0E0*2*I_ESP_G2xyz_G3yz_a-2.0E0*1*I_ESP_G3yz_G3yz_a+2*I_ESP_Dyz_G3yz;
    abcd[iGrid*1350+398] = 4.0E0*I_ESP_I2x2y2z_G3yz_aa-2.0E0*1*I_ESP_G2x2z_G3yz_a-2.0E0*1*I_ESP_G2y2z_G3yz_a+1*I_ESP_D2z_G3yz;
    abcd[iGrid*1350+399] = 4.0E0*I_ESP_I2xy3z_G3yz_aa-2.0E0*1*I_ESP_Gy3z_G3yz_a;
    abcd[iGrid*1350+400] = 4.0E0*I_ESP_Ix5y_G3yz_aa-2.0E0*4*I_ESP_Gx3y_G3yz_a;
    abcd[iGrid*1350+401] = 4.0E0*I_ESP_Ix4yz_G3yz_aa-2.0E0*3*I_ESP_Gx2yz_G3yz_a;
    abcd[iGrid*1350+402] = 4.0E0*I_ESP_Ix3y2z_G3yz_aa-2.0E0*2*I_ESP_Gxy2z_G3yz_a;
    abcd[iGrid*1350+403] = 4.0E0*I_ESP_Ix2y3z_G3yz_aa-2.0E0*1*I_ESP_Gx3z_G3yz_a;
    abcd[iGrid*1350+404] = 4.0E0*I_ESP_Ixy4z_G3yz_aa;
    abcd[iGrid*1350+405] = 4.0E0*I_ESP_I5xy_G2y2z_aa-2.0E0*4*I_ESP_G3xy_G2y2z_a;
    abcd[iGrid*1350+406] = 4.0E0*I_ESP_I4x2y_G2y2z_aa-2.0E0*1*I_ESP_G4x_G2y2z_a-2.0E0*3*I_ESP_G2x2y_G2y2z_a+3*1*I_ESP_D2x_G2y2z;
    abcd[iGrid*1350+407] = 4.0E0*I_ESP_I4xyz_G2y2z_aa-2.0E0*3*I_ESP_G2xyz_G2y2z_a;
    abcd[iGrid*1350+408] = 4.0E0*I_ESP_I3x3y_G2y2z_aa-2.0E0*2*I_ESP_G3xy_G2y2z_a-2.0E0*2*I_ESP_Gx3y_G2y2z_a+2*2*I_ESP_Dxy_G2y2z;
    abcd[iGrid*1350+409] = 4.0E0*I_ESP_I3x2yz_G2y2z_aa-2.0E0*1*I_ESP_G3xz_G2y2z_a-2.0E0*2*I_ESP_Gx2yz_G2y2z_a+2*1*I_ESP_Dxz_G2y2z;
    abcd[iGrid*1350+410] = 4.0E0*I_ESP_I3xy2z_G2y2z_aa-2.0E0*2*I_ESP_Gxy2z_G2y2z_a;
    abcd[iGrid*1350+411] = 4.0E0*I_ESP_I2x4y_G2y2z_aa-2.0E0*3*I_ESP_G2x2y_G2y2z_a-2.0E0*1*I_ESP_G4y_G2y2z_a+3*I_ESP_D2y_G2y2z;
    abcd[iGrid*1350+412] = 4.0E0*I_ESP_I2x3yz_G2y2z_aa-2.0E0*2*I_ESP_G2xyz_G2y2z_a-2.0E0*1*I_ESP_G3yz_G2y2z_a+2*I_ESP_Dyz_G2y2z;
    abcd[iGrid*1350+413] = 4.0E0*I_ESP_I2x2y2z_G2y2z_aa-2.0E0*1*I_ESP_G2x2z_G2y2z_a-2.0E0*1*I_ESP_G2y2z_G2y2z_a+1*I_ESP_D2z_G2y2z;
    abcd[iGrid*1350+414] = 4.0E0*I_ESP_I2xy3z_G2y2z_aa-2.0E0*1*I_ESP_Gy3z_G2y2z_a;
    abcd[iGrid*1350+415] = 4.0E0*I_ESP_Ix5y_G2y2z_aa-2.0E0*4*I_ESP_Gx3y_G2y2z_a;
    abcd[iGrid*1350+416] = 4.0E0*I_ESP_Ix4yz_G2y2z_aa-2.0E0*3*I_ESP_Gx2yz_G2y2z_a;
    abcd[iGrid*1350+417] = 4.0E0*I_ESP_Ix3y2z_G2y2z_aa-2.0E0*2*I_ESP_Gxy2z_G2y2z_a;
    abcd[iGrid*1350+418] = 4.0E0*I_ESP_Ix2y3z_G2y2z_aa-2.0E0*1*I_ESP_Gx3z_G2y2z_a;
    abcd[iGrid*1350+419] = 4.0E0*I_ESP_Ixy4z_G2y2z_aa;
    abcd[iGrid*1350+420] = 4.0E0*I_ESP_I5xy_Gy3z_aa-2.0E0*4*I_ESP_G3xy_Gy3z_a;
    abcd[iGrid*1350+421] = 4.0E0*I_ESP_I4x2y_Gy3z_aa-2.0E0*1*I_ESP_G4x_Gy3z_a-2.0E0*3*I_ESP_G2x2y_Gy3z_a+3*1*I_ESP_D2x_Gy3z;
    abcd[iGrid*1350+422] = 4.0E0*I_ESP_I4xyz_Gy3z_aa-2.0E0*3*I_ESP_G2xyz_Gy3z_a;
    abcd[iGrid*1350+423] = 4.0E0*I_ESP_I3x3y_Gy3z_aa-2.0E0*2*I_ESP_G3xy_Gy3z_a-2.0E0*2*I_ESP_Gx3y_Gy3z_a+2*2*I_ESP_Dxy_Gy3z;
    abcd[iGrid*1350+424] = 4.0E0*I_ESP_I3x2yz_Gy3z_aa-2.0E0*1*I_ESP_G3xz_Gy3z_a-2.0E0*2*I_ESP_Gx2yz_Gy3z_a+2*1*I_ESP_Dxz_Gy3z;
    abcd[iGrid*1350+425] = 4.0E0*I_ESP_I3xy2z_Gy3z_aa-2.0E0*2*I_ESP_Gxy2z_Gy3z_a;
    abcd[iGrid*1350+426] = 4.0E0*I_ESP_I2x4y_Gy3z_aa-2.0E0*3*I_ESP_G2x2y_Gy3z_a-2.0E0*1*I_ESP_G4y_Gy3z_a+3*I_ESP_D2y_Gy3z;
    abcd[iGrid*1350+427] = 4.0E0*I_ESP_I2x3yz_Gy3z_aa-2.0E0*2*I_ESP_G2xyz_Gy3z_a-2.0E0*1*I_ESP_G3yz_Gy3z_a+2*I_ESP_Dyz_Gy3z;
    abcd[iGrid*1350+428] = 4.0E0*I_ESP_I2x2y2z_Gy3z_aa-2.0E0*1*I_ESP_G2x2z_Gy3z_a-2.0E0*1*I_ESP_G2y2z_Gy3z_a+1*I_ESP_D2z_Gy3z;
    abcd[iGrid*1350+429] = 4.0E0*I_ESP_I2xy3z_Gy3z_aa-2.0E0*1*I_ESP_Gy3z_Gy3z_a;
    abcd[iGrid*1350+430] = 4.0E0*I_ESP_Ix5y_Gy3z_aa-2.0E0*4*I_ESP_Gx3y_Gy3z_a;
    abcd[iGrid*1350+431] = 4.0E0*I_ESP_Ix4yz_Gy3z_aa-2.0E0*3*I_ESP_Gx2yz_Gy3z_a;
    abcd[iGrid*1350+432] = 4.0E0*I_ESP_Ix3y2z_Gy3z_aa-2.0E0*2*I_ESP_Gxy2z_Gy3z_a;
    abcd[iGrid*1350+433] = 4.0E0*I_ESP_Ix2y3z_Gy3z_aa-2.0E0*1*I_ESP_Gx3z_Gy3z_a;
    abcd[iGrid*1350+434] = 4.0E0*I_ESP_Ixy4z_Gy3z_aa;
    abcd[iGrid*1350+435] = 4.0E0*I_ESP_I5xy_G4z_aa-2.0E0*4*I_ESP_G3xy_G4z_a;
    abcd[iGrid*1350+436] = 4.0E0*I_ESP_I4x2y_G4z_aa-2.0E0*1*I_ESP_G4x_G4z_a-2.0E0*3*I_ESP_G2x2y_G4z_a+3*1*I_ESP_D2x_G4z;
    abcd[iGrid*1350+437] = 4.0E0*I_ESP_I4xyz_G4z_aa-2.0E0*3*I_ESP_G2xyz_G4z_a;
    abcd[iGrid*1350+438] = 4.0E0*I_ESP_I3x3y_G4z_aa-2.0E0*2*I_ESP_G3xy_G4z_a-2.0E0*2*I_ESP_Gx3y_G4z_a+2*2*I_ESP_Dxy_G4z;
    abcd[iGrid*1350+439] = 4.0E0*I_ESP_I3x2yz_G4z_aa-2.0E0*1*I_ESP_G3xz_G4z_a-2.0E0*2*I_ESP_Gx2yz_G4z_a+2*1*I_ESP_Dxz_G4z;
    abcd[iGrid*1350+440] = 4.0E0*I_ESP_I3xy2z_G4z_aa-2.0E0*2*I_ESP_Gxy2z_G4z_a;
    abcd[iGrid*1350+441] = 4.0E0*I_ESP_I2x4y_G4z_aa-2.0E0*3*I_ESP_G2x2y_G4z_a-2.0E0*1*I_ESP_G4y_G4z_a+3*I_ESP_D2y_G4z;
    abcd[iGrid*1350+442] = 4.0E0*I_ESP_I2x3yz_G4z_aa-2.0E0*2*I_ESP_G2xyz_G4z_a-2.0E0*1*I_ESP_G3yz_G4z_a+2*I_ESP_Dyz_G4z;
    abcd[iGrid*1350+443] = 4.0E0*I_ESP_I2x2y2z_G4z_aa-2.0E0*1*I_ESP_G2x2z_G4z_a-2.0E0*1*I_ESP_G2y2z_G4z_a+1*I_ESP_D2z_G4z;
    abcd[iGrid*1350+444] = 4.0E0*I_ESP_I2xy3z_G4z_aa-2.0E0*1*I_ESP_Gy3z_G4z_a;
    abcd[iGrid*1350+445] = 4.0E0*I_ESP_Ix5y_G4z_aa-2.0E0*4*I_ESP_Gx3y_G4z_a;
    abcd[iGrid*1350+446] = 4.0E0*I_ESP_Ix4yz_G4z_aa-2.0E0*3*I_ESP_Gx2yz_G4z_a;
    abcd[iGrid*1350+447] = 4.0E0*I_ESP_Ix3y2z_G4z_aa-2.0E0*2*I_ESP_Gxy2z_G4z_a;
    abcd[iGrid*1350+448] = 4.0E0*I_ESP_Ix2y3z_G4z_aa-2.0E0*1*I_ESP_Gx3z_G4z_a;
    abcd[iGrid*1350+449] = 4.0E0*I_ESP_Ixy4z_G4z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G_aa
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_D_G
     ************************************************************/
    abcd[iGrid*1350+450] = 4.0E0*I_ESP_I5xz_G4x_aa-2.0E0*4*I_ESP_G3xz_G4x_a;
    abcd[iGrid*1350+451] = 4.0E0*I_ESP_I4xyz_G4x_aa-2.0E0*3*I_ESP_G2xyz_G4x_a;
    abcd[iGrid*1350+452] = 4.0E0*I_ESP_I4x2z_G4x_aa-2.0E0*1*I_ESP_G4x_G4x_a-2.0E0*3*I_ESP_G2x2z_G4x_a+3*1*I_ESP_D2x_G4x;
    abcd[iGrid*1350+453] = 4.0E0*I_ESP_I3x2yz_G4x_aa-2.0E0*2*I_ESP_Gx2yz_G4x_a;
    abcd[iGrid*1350+454] = 4.0E0*I_ESP_I3xy2z_G4x_aa-2.0E0*1*I_ESP_G3xy_G4x_a-2.0E0*2*I_ESP_Gxy2z_G4x_a+2*1*I_ESP_Dxy_G4x;
    abcd[iGrid*1350+455] = 4.0E0*I_ESP_I3x3z_G4x_aa-2.0E0*2*I_ESP_G3xz_G4x_a-2.0E0*2*I_ESP_Gx3z_G4x_a+2*2*I_ESP_Dxz_G4x;
    abcd[iGrid*1350+456] = 4.0E0*I_ESP_I2x3yz_G4x_aa-2.0E0*1*I_ESP_G3yz_G4x_a;
    abcd[iGrid*1350+457] = 4.0E0*I_ESP_I2x2y2z_G4x_aa-2.0E0*1*I_ESP_G2x2y_G4x_a-2.0E0*1*I_ESP_G2y2z_G4x_a+1*I_ESP_D2y_G4x;
    abcd[iGrid*1350+458] = 4.0E0*I_ESP_I2xy3z_G4x_aa-2.0E0*2*I_ESP_G2xyz_G4x_a-2.0E0*1*I_ESP_Gy3z_G4x_a+2*I_ESP_Dyz_G4x;
    abcd[iGrid*1350+459] = 4.0E0*I_ESP_I2x4z_G4x_aa-2.0E0*3*I_ESP_G2x2z_G4x_a-2.0E0*1*I_ESP_G4z_G4x_a+3*I_ESP_D2z_G4x;
    abcd[iGrid*1350+460] = 4.0E0*I_ESP_Ix4yz_G4x_aa;
    abcd[iGrid*1350+461] = 4.0E0*I_ESP_Ix3y2z_G4x_aa-2.0E0*1*I_ESP_Gx3y_G4x_a;
    abcd[iGrid*1350+462] = 4.0E0*I_ESP_Ix2y3z_G4x_aa-2.0E0*2*I_ESP_Gx2yz_G4x_a;
    abcd[iGrid*1350+463] = 4.0E0*I_ESP_Ixy4z_G4x_aa-2.0E0*3*I_ESP_Gxy2z_G4x_a;
    abcd[iGrid*1350+464] = 4.0E0*I_ESP_Ix5z_G4x_aa-2.0E0*4*I_ESP_Gx3z_G4x_a;
    abcd[iGrid*1350+465] = 4.0E0*I_ESP_I5xz_G3xy_aa-2.0E0*4*I_ESP_G3xz_G3xy_a;
    abcd[iGrid*1350+466] = 4.0E0*I_ESP_I4xyz_G3xy_aa-2.0E0*3*I_ESP_G2xyz_G3xy_a;
    abcd[iGrid*1350+467] = 4.0E0*I_ESP_I4x2z_G3xy_aa-2.0E0*1*I_ESP_G4x_G3xy_a-2.0E0*3*I_ESP_G2x2z_G3xy_a+3*1*I_ESP_D2x_G3xy;
    abcd[iGrid*1350+468] = 4.0E0*I_ESP_I3x2yz_G3xy_aa-2.0E0*2*I_ESP_Gx2yz_G3xy_a;
    abcd[iGrid*1350+469] = 4.0E0*I_ESP_I3xy2z_G3xy_aa-2.0E0*1*I_ESP_G3xy_G3xy_a-2.0E0*2*I_ESP_Gxy2z_G3xy_a+2*1*I_ESP_Dxy_G3xy;
    abcd[iGrid*1350+470] = 4.0E0*I_ESP_I3x3z_G3xy_aa-2.0E0*2*I_ESP_G3xz_G3xy_a-2.0E0*2*I_ESP_Gx3z_G3xy_a+2*2*I_ESP_Dxz_G3xy;
    abcd[iGrid*1350+471] = 4.0E0*I_ESP_I2x3yz_G3xy_aa-2.0E0*1*I_ESP_G3yz_G3xy_a;
    abcd[iGrid*1350+472] = 4.0E0*I_ESP_I2x2y2z_G3xy_aa-2.0E0*1*I_ESP_G2x2y_G3xy_a-2.0E0*1*I_ESP_G2y2z_G3xy_a+1*I_ESP_D2y_G3xy;
    abcd[iGrid*1350+473] = 4.0E0*I_ESP_I2xy3z_G3xy_aa-2.0E0*2*I_ESP_G2xyz_G3xy_a-2.0E0*1*I_ESP_Gy3z_G3xy_a+2*I_ESP_Dyz_G3xy;
    abcd[iGrid*1350+474] = 4.0E0*I_ESP_I2x4z_G3xy_aa-2.0E0*3*I_ESP_G2x2z_G3xy_a-2.0E0*1*I_ESP_G4z_G3xy_a+3*I_ESP_D2z_G3xy;
    abcd[iGrid*1350+475] = 4.0E0*I_ESP_Ix4yz_G3xy_aa;
    abcd[iGrid*1350+476] = 4.0E0*I_ESP_Ix3y2z_G3xy_aa-2.0E0*1*I_ESP_Gx3y_G3xy_a;
    abcd[iGrid*1350+477] = 4.0E0*I_ESP_Ix2y3z_G3xy_aa-2.0E0*2*I_ESP_Gx2yz_G3xy_a;
    abcd[iGrid*1350+478] = 4.0E0*I_ESP_Ixy4z_G3xy_aa-2.0E0*3*I_ESP_Gxy2z_G3xy_a;
    abcd[iGrid*1350+479] = 4.0E0*I_ESP_Ix5z_G3xy_aa-2.0E0*4*I_ESP_Gx3z_G3xy_a;
    abcd[iGrid*1350+480] = 4.0E0*I_ESP_I5xz_G3xz_aa-2.0E0*4*I_ESP_G3xz_G3xz_a;
    abcd[iGrid*1350+481] = 4.0E0*I_ESP_I4xyz_G3xz_aa-2.0E0*3*I_ESP_G2xyz_G3xz_a;
    abcd[iGrid*1350+482] = 4.0E0*I_ESP_I4x2z_G3xz_aa-2.0E0*1*I_ESP_G4x_G3xz_a-2.0E0*3*I_ESP_G2x2z_G3xz_a+3*1*I_ESP_D2x_G3xz;
    abcd[iGrid*1350+483] = 4.0E0*I_ESP_I3x2yz_G3xz_aa-2.0E0*2*I_ESP_Gx2yz_G3xz_a;
    abcd[iGrid*1350+484] = 4.0E0*I_ESP_I3xy2z_G3xz_aa-2.0E0*1*I_ESP_G3xy_G3xz_a-2.0E0*2*I_ESP_Gxy2z_G3xz_a+2*1*I_ESP_Dxy_G3xz;
    abcd[iGrid*1350+485] = 4.0E0*I_ESP_I3x3z_G3xz_aa-2.0E0*2*I_ESP_G3xz_G3xz_a-2.0E0*2*I_ESP_Gx3z_G3xz_a+2*2*I_ESP_Dxz_G3xz;
    abcd[iGrid*1350+486] = 4.0E0*I_ESP_I2x3yz_G3xz_aa-2.0E0*1*I_ESP_G3yz_G3xz_a;
    abcd[iGrid*1350+487] = 4.0E0*I_ESP_I2x2y2z_G3xz_aa-2.0E0*1*I_ESP_G2x2y_G3xz_a-2.0E0*1*I_ESP_G2y2z_G3xz_a+1*I_ESP_D2y_G3xz;
    abcd[iGrid*1350+488] = 4.0E0*I_ESP_I2xy3z_G3xz_aa-2.0E0*2*I_ESP_G2xyz_G3xz_a-2.0E0*1*I_ESP_Gy3z_G3xz_a+2*I_ESP_Dyz_G3xz;
    abcd[iGrid*1350+489] = 4.0E0*I_ESP_I2x4z_G3xz_aa-2.0E0*3*I_ESP_G2x2z_G3xz_a-2.0E0*1*I_ESP_G4z_G3xz_a+3*I_ESP_D2z_G3xz;
    abcd[iGrid*1350+490] = 4.0E0*I_ESP_Ix4yz_G3xz_aa;
    abcd[iGrid*1350+491] = 4.0E0*I_ESP_Ix3y2z_G3xz_aa-2.0E0*1*I_ESP_Gx3y_G3xz_a;
    abcd[iGrid*1350+492] = 4.0E0*I_ESP_Ix2y3z_G3xz_aa-2.0E0*2*I_ESP_Gx2yz_G3xz_a;
    abcd[iGrid*1350+493] = 4.0E0*I_ESP_Ixy4z_G3xz_aa-2.0E0*3*I_ESP_Gxy2z_G3xz_a;
    abcd[iGrid*1350+494] = 4.0E0*I_ESP_Ix5z_G3xz_aa-2.0E0*4*I_ESP_Gx3z_G3xz_a;
    abcd[iGrid*1350+495] = 4.0E0*I_ESP_I5xz_G2x2y_aa-2.0E0*4*I_ESP_G3xz_G2x2y_a;
    abcd[iGrid*1350+496] = 4.0E0*I_ESP_I4xyz_G2x2y_aa-2.0E0*3*I_ESP_G2xyz_G2x2y_a;
    abcd[iGrid*1350+497] = 4.0E0*I_ESP_I4x2z_G2x2y_aa-2.0E0*1*I_ESP_G4x_G2x2y_a-2.0E0*3*I_ESP_G2x2z_G2x2y_a+3*1*I_ESP_D2x_G2x2y;
    abcd[iGrid*1350+498] = 4.0E0*I_ESP_I3x2yz_G2x2y_aa-2.0E0*2*I_ESP_Gx2yz_G2x2y_a;
    abcd[iGrid*1350+499] = 4.0E0*I_ESP_I3xy2z_G2x2y_aa-2.0E0*1*I_ESP_G3xy_G2x2y_a-2.0E0*2*I_ESP_Gxy2z_G2x2y_a+2*1*I_ESP_Dxy_G2x2y;
    abcd[iGrid*1350+500] = 4.0E0*I_ESP_I3x3z_G2x2y_aa-2.0E0*2*I_ESP_G3xz_G2x2y_a-2.0E0*2*I_ESP_Gx3z_G2x2y_a+2*2*I_ESP_Dxz_G2x2y;
    abcd[iGrid*1350+501] = 4.0E0*I_ESP_I2x3yz_G2x2y_aa-2.0E0*1*I_ESP_G3yz_G2x2y_a;
    abcd[iGrid*1350+502] = 4.0E0*I_ESP_I2x2y2z_G2x2y_aa-2.0E0*1*I_ESP_G2x2y_G2x2y_a-2.0E0*1*I_ESP_G2y2z_G2x2y_a+1*I_ESP_D2y_G2x2y;
    abcd[iGrid*1350+503] = 4.0E0*I_ESP_I2xy3z_G2x2y_aa-2.0E0*2*I_ESP_G2xyz_G2x2y_a-2.0E0*1*I_ESP_Gy3z_G2x2y_a+2*I_ESP_Dyz_G2x2y;
    abcd[iGrid*1350+504] = 4.0E0*I_ESP_I2x4z_G2x2y_aa-2.0E0*3*I_ESP_G2x2z_G2x2y_a-2.0E0*1*I_ESP_G4z_G2x2y_a+3*I_ESP_D2z_G2x2y;
    abcd[iGrid*1350+505] = 4.0E0*I_ESP_Ix4yz_G2x2y_aa;
    abcd[iGrid*1350+506] = 4.0E0*I_ESP_Ix3y2z_G2x2y_aa-2.0E0*1*I_ESP_Gx3y_G2x2y_a;
    abcd[iGrid*1350+507] = 4.0E0*I_ESP_Ix2y3z_G2x2y_aa-2.0E0*2*I_ESP_Gx2yz_G2x2y_a;
    abcd[iGrid*1350+508] = 4.0E0*I_ESP_Ixy4z_G2x2y_aa-2.0E0*3*I_ESP_Gxy2z_G2x2y_a;
    abcd[iGrid*1350+509] = 4.0E0*I_ESP_Ix5z_G2x2y_aa-2.0E0*4*I_ESP_Gx3z_G2x2y_a;
    abcd[iGrid*1350+510] = 4.0E0*I_ESP_I5xz_G2xyz_aa-2.0E0*4*I_ESP_G3xz_G2xyz_a;
    abcd[iGrid*1350+511] = 4.0E0*I_ESP_I4xyz_G2xyz_aa-2.0E0*3*I_ESP_G2xyz_G2xyz_a;
    abcd[iGrid*1350+512] = 4.0E0*I_ESP_I4x2z_G2xyz_aa-2.0E0*1*I_ESP_G4x_G2xyz_a-2.0E0*3*I_ESP_G2x2z_G2xyz_a+3*1*I_ESP_D2x_G2xyz;
    abcd[iGrid*1350+513] = 4.0E0*I_ESP_I3x2yz_G2xyz_aa-2.0E0*2*I_ESP_Gx2yz_G2xyz_a;
    abcd[iGrid*1350+514] = 4.0E0*I_ESP_I3xy2z_G2xyz_aa-2.0E0*1*I_ESP_G3xy_G2xyz_a-2.0E0*2*I_ESP_Gxy2z_G2xyz_a+2*1*I_ESP_Dxy_G2xyz;
    abcd[iGrid*1350+515] = 4.0E0*I_ESP_I3x3z_G2xyz_aa-2.0E0*2*I_ESP_G3xz_G2xyz_a-2.0E0*2*I_ESP_Gx3z_G2xyz_a+2*2*I_ESP_Dxz_G2xyz;
    abcd[iGrid*1350+516] = 4.0E0*I_ESP_I2x3yz_G2xyz_aa-2.0E0*1*I_ESP_G3yz_G2xyz_a;
    abcd[iGrid*1350+517] = 4.0E0*I_ESP_I2x2y2z_G2xyz_aa-2.0E0*1*I_ESP_G2x2y_G2xyz_a-2.0E0*1*I_ESP_G2y2z_G2xyz_a+1*I_ESP_D2y_G2xyz;
    abcd[iGrid*1350+518] = 4.0E0*I_ESP_I2xy3z_G2xyz_aa-2.0E0*2*I_ESP_G2xyz_G2xyz_a-2.0E0*1*I_ESP_Gy3z_G2xyz_a+2*I_ESP_Dyz_G2xyz;
    abcd[iGrid*1350+519] = 4.0E0*I_ESP_I2x4z_G2xyz_aa-2.0E0*3*I_ESP_G2x2z_G2xyz_a-2.0E0*1*I_ESP_G4z_G2xyz_a+3*I_ESP_D2z_G2xyz;
    abcd[iGrid*1350+520] = 4.0E0*I_ESP_Ix4yz_G2xyz_aa;
    abcd[iGrid*1350+521] = 4.0E0*I_ESP_Ix3y2z_G2xyz_aa-2.0E0*1*I_ESP_Gx3y_G2xyz_a;
    abcd[iGrid*1350+522] = 4.0E0*I_ESP_Ix2y3z_G2xyz_aa-2.0E0*2*I_ESP_Gx2yz_G2xyz_a;
    abcd[iGrid*1350+523] = 4.0E0*I_ESP_Ixy4z_G2xyz_aa-2.0E0*3*I_ESP_Gxy2z_G2xyz_a;
    abcd[iGrid*1350+524] = 4.0E0*I_ESP_Ix5z_G2xyz_aa-2.0E0*4*I_ESP_Gx3z_G2xyz_a;
    abcd[iGrid*1350+525] = 4.0E0*I_ESP_I5xz_G2x2z_aa-2.0E0*4*I_ESP_G3xz_G2x2z_a;
    abcd[iGrid*1350+526] = 4.0E0*I_ESP_I4xyz_G2x2z_aa-2.0E0*3*I_ESP_G2xyz_G2x2z_a;
    abcd[iGrid*1350+527] = 4.0E0*I_ESP_I4x2z_G2x2z_aa-2.0E0*1*I_ESP_G4x_G2x2z_a-2.0E0*3*I_ESP_G2x2z_G2x2z_a+3*1*I_ESP_D2x_G2x2z;
    abcd[iGrid*1350+528] = 4.0E0*I_ESP_I3x2yz_G2x2z_aa-2.0E0*2*I_ESP_Gx2yz_G2x2z_a;
    abcd[iGrid*1350+529] = 4.0E0*I_ESP_I3xy2z_G2x2z_aa-2.0E0*1*I_ESP_G3xy_G2x2z_a-2.0E0*2*I_ESP_Gxy2z_G2x2z_a+2*1*I_ESP_Dxy_G2x2z;
    abcd[iGrid*1350+530] = 4.0E0*I_ESP_I3x3z_G2x2z_aa-2.0E0*2*I_ESP_G3xz_G2x2z_a-2.0E0*2*I_ESP_Gx3z_G2x2z_a+2*2*I_ESP_Dxz_G2x2z;
    abcd[iGrid*1350+531] = 4.0E0*I_ESP_I2x3yz_G2x2z_aa-2.0E0*1*I_ESP_G3yz_G2x2z_a;
    abcd[iGrid*1350+532] = 4.0E0*I_ESP_I2x2y2z_G2x2z_aa-2.0E0*1*I_ESP_G2x2y_G2x2z_a-2.0E0*1*I_ESP_G2y2z_G2x2z_a+1*I_ESP_D2y_G2x2z;
    abcd[iGrid*1350+533] = 4.0E0*I_ESP_I2xy3z_G2x2z_aa-2.0E0*2*I_ESP_G2xyz_G2x2z_a-2.0E0*1*I_ESP_Gy3z_G2x2z_a+2*I_ESP_Dyz_G2x2z;
    abcd[iGrid*1350+534] = 4.0E0*I_ESP_I2x4z_G2x2z_aa-2.0E0*3*I_ESP_G2x2z_G2x2z_a-2.0E0*1*I_ESP_G4z_G2x2z_a+3*I_ESP_D2z_G2x2z;
    abcd[iGrid*1350+535] = 4.0E0*I_ESP_Ix4yz_G2x2z_aa;
    abcd[iGrid*1350+536] = 4.0E0*I_ESP_Ix3y2z_G2x2z_aa-2.0E0*1*I_ESP_Gx3y_G2x2z_a;
    abcd[iGrid*1350+537] = 4.0E0*I_ESP_Ix2y3z_G2x2z_aa-2.0E0*2*I_ESP_Gx2yz_G2x2z_a;
    abcd[iGrid*1350+538] = 4.0E0*I_ESP_Ixy4z_G2x2z_aa-2.0E0*3*I_ESP_Gxy2z_G2x2z_a;
    abcd[iGrid*1350+539] = 4.0E0*I_ESP_Ix5z_G2x2z_aa-2.0E0*4*I_ESP_Gx3z_G2x2z_a;
    abcd[iGrid*1350+540] = 4.0E0*I_ESP_I5xz_Gx3y_aa-2.0E0*4*I_ESP_G3xz_Gx3y_a;
    abcd[iGrid*1350+541] = 4.0E0*I_ESP_I4xyz_Gx3y_aa-2.0E0*3*I_ESP_G2xyz_Gx3y_a;
    abcd[iGrid*1350+542] = 4.0E0*I_ESP_I4x2z_Gx3y_aa-2.0E0*1*I_ESP_G4x_Gx3y_a-2.0E0*3*I_ESP_G2x2z_Gx3y_a+3*1*I_ESP_D2x_Gx3y;
    abcd[iGrid*1350+543] = 4.0E0*I_ESP_I3x2yz_Gx3y_aa-2.0E0*2*I_ESP_Gx2yz_Gx3y_a;
    abcd[iGrid*1350+544] = 4.0E0*I_ESP_I3xy2z_Gx3y_aa-2.0E0*1*I_ESP_G3xy_Gx3y_a-2.0E0*2*I_ESP_Gxy2z_Gx3y_a+2*1*I_ESP_Dxy_Gx3y;
    abcd[iGrid*1350+545] = 4.0E0*I_ESP_I3x3z_Gx3y_aa-2.0E0*2*I_ESP_G3xz_Gx3y_a-2.0E0*2*I_ESP_Gx3z_Gx3y_a+2*2*I_ESP_Dxz_Gx3y;
    abcd[iGrid*1350+546] = 4.0E0*I_ESP_I2x3yz_Gx3y_aa-2.0E0*1*I_ESP_G3yz_Gx3y_a;
    abcd[iGrid*1350+547] = 4.0E0*I_ESP_I2x2y2z_Gx3y_aa-2.0E0*1*I_ESP_G2x2y_Gx3y_a-2.0E0*1*I_ESP_G2y2z_Gx3y_a+1*I_ESP_D2y_Gx3y;
    abcd[iGrid*1350+548] = 4.0E0*I_ESP_I2xy3z_Gx3y_aa-2.0E0*2*I_ESP_G2xyz_Gx3y_a-2.0E0*1*I_ESP_Gy3z_Gx3y_a+2*I_ESP_Dyz_Gx3y;
    abcd[iGrid*1350+549] = 4.0E0*I_ESP_I2x4z_Gx3y_aa-2.0E0*3*I_ESP_G2x2z_Gx3y_a-2.0E0*1*I_ESP_G4z_Gx3y_a+3*I_ESP_D2z_Gx3y;
    abcd[iGrid*1350+550] = 4.0E0*I_ESP_Ix4yz_Gx3y_aa;
    abcd[iGrid*1350+551] = 4.0E0*I_ESP_Ix3y2z_Gx3y_aa-2.0E0*1*I_ESP_Gx3y_Gx3y_a;
    abcd[iGrid*1350+552] = 4.0E0*I_ESP_Ix2y3z_Gx3y_aa-2.0E0*2*I_ESP_Gx2yz_Gx3y_a;
    abcd[iGrid*1350+553] = 4.0E0*I_ESP_Ixy4z_Gx3y_aa-2.0E0*3*I_ESP_Gxy2z_Gx3y_a;
    abcd[iGrid*1350+554] = 4.0E0*I_ESP_Ix5z_Gx3y_aa-2.0E0*4*I_ESP_Gx3z_Gx3y_a;
    abcd[iGrid*1350+555] = 4.0E0*I_ESP_I5xz_Gx2yz_aa-2.0E0*4*I_ESP_G3xz_Gx2yz_a;
    abcd[iGrid*1350+556] = 4.0E0*I_ESP_I4xyz_Gx2yz_aa-2.0E0*3*I_ESP_G2xyz_Gx2yz_a;
    abcd[iGrid*1350+557] = 4.0E0*I_ESP_I4x2z_Gx2yz_aa-2.0E0*1*I_ESP_G4x_Gx2yz_a-2.0E0*3*I_ESP_G2x2z_Gx2yz_a+3*1*I_ESP_D2x_Gx2yz;
    abcd[iGrid*1350+558] = 4.0E0*I_ESP_I3x2yz_Gx2yz_aa-2.0E0*2*I_ESP_Gx2yz_Gx2yz_a;
    abcd[iGrid*1350+559] = 4.0E0*I_ESP_I3xy2z_Gx2yz_aa-2.0E0*1*I_ESP_G3xy_Gx2yz_a-2.0E0*2*I_ESP_Gxy2z_Gx2yz_a+2*1*I_ESP_Dxy_Gx2yz;
    abcd[iGrid*1350+560] = 4.0E0*I_ESP_I3x3z_Gx2yz_aa-2.0E0*2*I_ESP_G3xz_Gx2yz_a-2.0E0*2*I_ESP_Gx3z_Gx2yz_a+2*2*I_ESP_Dxz_Gx2yz;
    abcd[iGrid*1350+561] = 4.0E0*I_ESP_I2x3yz_Gx2yz_aa-2.0E0*1*I_ESP_G3yz_Gx2yz_a;
    abcd[iGrid*1350+562] = 4.0E0*I_ESP_I2x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_G2x2y_Gx2yz_a-2.0E0*1*I_ESP_G2y2z_Gx2yz_a+1*I_ESP_D2y_Gx2yz;
    abcd[iGrid*1350+563] = 4.0E0*I_ESP_I2xy3z_Gx2yz_aa-2.0E0*2*I_ESP_G2xyz_Gx2yz_a-2.0E0*1*I_ESP_Gy3z_Gx2yz_a+2*I_ESP_Dyz_Gx2yz;
    abcd[iGrid*1350+564] = 4.0E0*I_ESP_I2x4z_Gx2yz_aa-2.0E0*3*I_ESP_G2x2z_Gx2yz_a-2.0E0*1*I_ESP_G4z_Gx2yz_a+3*I_ESP_D2z_Gx2yz;
    abcd[iGrid*1350+565] = 4.0E0*I_ESP_Ix4yz_Gx2yz_aa;
    abcd[iGrid*1350+566] = 4.0E0*I_ESP_Ix3y2z_Gx2yz_aa-2.0E0*1*I_ESP_Gx3y_Gx2yz_a;
    abcd[iGrid*1350+567] = 4.0E0*I_ESP_Ix2y3z_Gx2yz_aa-2.0E0*2*I_ESP_Gx2yz_Gx2yz_a;
    abcd[iGrid*1350+568] = 4.0E0*I_ESP_Ixy4z_Gx2yz_aa-2.0E0*3*I_ESP_Gxy2z_Gx2yz_a;
    abcd[iGrid*1350+569] = 4.0E0*I_ESP_Ix5z_Gx2yz_aa-2.0E0*4*I_ESP_Gx3z_Gx2yz_a;
    abcd[iGrid*1350+570] = 4.0E0*I_ESP_I5xz_Gxy2z_aa-2.0E0*4*I_ESP_G3xz_Gxy2z_a;
    abcd[iGrid*1350+571] = 4.0E0*I_ESP_I4xyz_Gxy2z_aa-2.0E0*3*I_ESP_G2xyz_Gxy2z_a;
    abcd[iGrid*1350+572] = 4.0E0*I_ESP_I4x2z_Gxy2z_aa-2.0E0*1*I_ESP_G4x_Gxy2z_a-2.0E0*3*I_ESP_G2x2z_Gxy2z_a+3*1*I_ESP_D2x_Gxy2z;
    abcd[iGrid*1350+573] = 4.0E0*I_ESP_I3x2yz_Gxy2z_aa-2.0E0*2*I_ESP_Gx2yz_Gxy2z_a;
    abcd[iGrid*1350+574] = 4.0E0*I_ESP_I3xy2z_Gxy2z_aa-2.0E0*1*I_ESP_G3xy_Gxy2z_a-2.0E0*2*I_ESP_Gxy2z_Gxy2z_a+2*1*I_ESP_Dxy_Gxy2z;
    abcd[iGrid*1350+575] = 4.0E0*I_ESP_I3x3z_Gxy2z_aa-2.0E0*2*I_ESP_G3xz_Gxy2z_a-2.0E0*2*I_ESP_Gx3z_Gxy2z_a+2*2*I_ESP_Dxz_Gxy2z;
    abcd[iGrid*1350+576] = 4.0E0*I_ESP_I2x3yz_Gxy2z_aa-2.0E0*1*I_ESP_G3yz_Gxy2z_a;
    abcd[iGrid*1350+577] = 4.0E0*I_ESP_I2x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_G2x2y_Gxy2z_a-2.0E0*1*I_ESP_G2y2z_Gxy2z_a+1*I_ESP_D2y_Gxy2z;
    abcd[iGrid*1350+578] = 4.0E0*I_ESP_I2xy3z_Gxy2z_aa-2.0E0*2*I_ESP_G2xyz_Gxy2z_a-2.0E0*1*I_ESP_Gy3z_Gxy2z_a+2*I_ESP_Dyz_Gxy2z;
    abcd[iGrid*1350+579] = 4.0E0*I_ESP_I2x4z_Gxy2z_aa-2.0E0*3*I_ESP_G2x2z_Gxy2z_a-2.0E0*1*I_ESP_G4z_Gxy2z_a+3*I_ESP_D2z_Gxy2z;
    abcd[iGrid*1350+580] = 4.0E0*I_ESP_Ix4yz_Gxy2z_aa;
    abcd[iGrid*1350+581] = 4.0E0*I_ESP_Ix3y2z_Gxy2z_aa-2.0E0*1*I_ESP_Gx3y_Gxy2z_a;
    abcd[iGrid*1350+582] = 4.0E0*I_ESP_Ix2y3z_Gxy2z_aa-2.0E0*2*I_ESP_Gx2yz_Gxy2z_a;
    abcd[iGrid*1350+583] = 4.0E0*I_ESP_Ixy4z_Gxy2z_aa-2.0E0*3*I_ESP_Gxy2z_Gxy2z_a;
    abcd[iGrid*1350+584] = 4.0E0*I_ESP_Ix5z_Gxy2z_aa-2.0E0*4*I_ESP_Gx3z_Gxy2z_a;
    abcd[iGrid*1350+585] = 4.0E0*I_ESP_I5xz_Gx3z_aa-2.0E0*4*I_ESP_G3xz_Gx3z_a;
    abcd[iGrid*1350+586] = 4.0E0*I_ESP_I4xyz_Gx3z_aa-2.0E0*3*I_ESP_G2xyz_Gx3z_a;
    abcd[iGrid*1350+587] = 4.0E0*I_ESP_I4x2z_Gx3z_aa-2.0E0*1*I_ESP_G4x_Gx3z_a-2.0E0*3*I_ESP_G2x2z_Gx3z_a+3*1*I_ESP_D2x_Gx3z;
    abcd[iGrid*1350+588] = 4.0E0*I_ESP_I3x2yz_Gx3z_aa-2.0E0*2*I_ESP_Gx2yz_Gx3z_a;
    abcd[iGrid*1350+589] = 4.0E0*I_ESP_I3xy2z_Gx3z_aa-2.0E0*1*I_ESP_G3xy_Gx3z_a-2.0E0*2*I_ESP_Gxy2z_Gx3z_a+2*1*I_ESP_Dxy_Gx3z;
    abcd[iGrid*1350+590] = 4.0E0*I_ESP_I3x3z_Gx3z_aa-2.0E0*2*I_ESP_G3xz_Gx3z_a-2.0E0*2*I_ESP_Gx3z_Gx3z_a+2*2*I_ESP_Dxz_Gx3z;
    abcd[iGrid*1350+591] = 4.0E0*I_ESP_I2x3yz_Gx3z_aa-2.0E0*1*I_ESP_G3yz_Gx3z_a;
    abcd[iGrid*1350+592] = 4.0E0*I_ESP_I2x2y2z_Gx3z_aa-2.0E0*1*I_ESP_G2x2y_Gx3z_a-2.0E0*1*I_ESP_G2y2z_Gx3z_a+1*I_ESP_D2y_Gx3z;
    abcd[iGrid*1350+593] = 4.0E0*I_ESP_I2xy3z_Gx3z_aa-2.0E0*2*I_ESP_G2xyz_Gx3z_a-2.0E0*1*I_ESP_Gy3z_Gx3z_a+2*I_ESP_Dyz_Gx3z;
    abcd[iGrid*1350+594] = 4.0E0*I_ESP_I2x4z_Gx3z_aa-2.0E0*3*I_ESP_G2x2z_Gx3z_a-2.0E0*1*I_ESP_G4z_Gx3z_a+3*I_ESP_D2z_Gx3z;
    abcd[iGrid*1350+595] = 4.0E0*I_ESP_Ix4yz_Gx3z_aa;
    abcd[iGrid*1350+596] = 4.0E0*I_ESP_Ix3y2z_Gx3z_aa-2.0E0*1*I_ESP_Gx3y_Gx3z_a;
    abcd[iGrid*1350+597] = 4.0E0*I_ESP_Ix2y3z_Gx3z_aa-2.0E0*2*I_ESP_Gx2yz_Gx3z_a;
    abcd[iGrid*1350+598] = 4.0E0*I_ESP_Ixy4z_Gx3z_aa-2.0E0*3*I_ESP_Gxy2z_Gx3z_a;
    abcd[iGrid*1350+599] = 4.0E0*I_ESP_Ix5z_Gx3z_aa-2.0E0*4*I_ESP_Gx3z_Gx3z_a;
    abcd[iGrid*1350+600] = 4.0E0*I_ESP_I5xz_G4y_aa-2.0E0*4*I_ESP_G3xz_G4y_a;
    abcd[iGrid*1350+601] = 4.0E0*I_ESP_I4xyz_G4y_aa-2.0E0*3*I_ESP_G2xyz_G4y_a;
    abcd[iGrid*1350+602] = 4.0E0*I_ESP_I4x2z_G4y_aa-2.0E0*1*I_ESP_G4x_G4y_a-2.0E0*3*I_ESP_G2x2z_G4y_a+3*1*I_ESP_D2x_G4y;
    abcd[iGrid*1350+603] = 4.0E0*I_ESP_I3x2yz_G4y_aa-2.0E0*2*I_ESP_Gx2yz_G4y_a;
    abcd[iGrid*1350+604] = 4.0E0*I_ESP_I3xy2z_G4y_aa-2.0E0*1*I_ESP_G3xy_G4y_a-2.0E0*2*I_ESP_Gxy2z_G4y_a+2*1*I_ESP_Dxy_G4y;
    abcd[iGrid*1350+605] = 4.0E0*I_ESP_I3x3z_G4y_aa-2.0E0*2*I_ESP_G3xz_G4y_a-2.0E0*2*I_ESP_Gx3z_G4y_a+2*2*I_ESP_Dxz_G4y;
    abcd[iGrid*1350+606] = 4.0E0*I_ESP_I2x3yz_G4y_aa-2.0E0*1*I_ESP_G3yz_G4y_a;
    abcd[iGrid*1350+607] = 4.0E0*I_ESP_I2x2y2z_G4y_aa-2.0E0*1*I_ESP_G2x2y_G4y_a-2.0E0*1*I_ESP_G2y2z_G4y_a+1*I_ESP_D2y_G4y;
    abcd[iGrid*1350+608] = 4.0E0*I_ESP_I2xy3z_G4y_aa-2.0E0*2*I_ESP_G2xyz_G4y_a-2.0E0*1*I_ESP_Gy3z_G4y_a+2*I_ESP_Dyz_G4y;
    abcd[iGrid*1350+609] = 4.0E0*I_ESP_I2x4z_G4y_aa-2.0E0*3*I_ESP_G2x2z_G4y_a-2.0E0*1*I_ESP_G4z_G4y_a+3*I_ESP_D2z_G4y;
    abcd[iGrid*1350+610] = 4.0E0*I_ESP_Ix4yz_G4y_aa;
    abcd[iGrid*1350+611] = 4.0E0*I_ESP_Ix3y2z_G4y_aa-2.0E0*1*I_ESP_Gx3y_G4y_a;
    abcd[iGrid*1350+612] = 4.0E0*I_ESP_Ix2y3z_G4y_aa-2.0E0*2*I_ESP_Gx2yz_G4y_a;
    abcd[iGrid*1350+613] = 4.0E0*I_ESP_Ixy4z_G4y_aa-2.0E0*3*I_ESP_Gxy2z_G4y_a;
    abcd[iGrid*1350+614] = 4.0E0*I_ESP_Ix5z_G4y_aa-2.0E0*4*I_ESP_Gx3z_G4y_a;
    abcd[iGrid*1350+615] = 4.0E0*I_ESP_I5xz_G3yz_aa-2.0E0*4*I_ESP_G3xz_G3yz_a;
    abcd[iGrid*1350+616] = 4.0E0*I_ESP_I4xyz_G3yz_aa-2.0E0*3*I_ESP_G2xyz_G3yz_a;
    abcd[iGrid*1350+617] = 4.0E0*I_ESP_I4x2z_G3yz_aa-2.0E0*1*I_ESP_G4x_G3yz_a-2.0E0*3*I_ESP_G2x2z_G3yz_a+3*1*I_ESP_D2x_G3yz;
    abcd[iGrid*1350+618] = 4.0E0*I_ESP_I3x2yz_G3yz_aa-2.0E0*2*I_ESP_Gx2yz_G3yz_a;
    abcd[iGrid*1350+619] = 4.0E0*I_ESP_I3xy2z_G3yz_aa-2.0E0*1*I_ESP_G3xy_G3yz_a-2.0E0*2*I_ESP_Gxy2z_G3yz_a+2*1*I_ESP_Dxy_G3yz;
    abcd[iGrid*1350+620] = 4.0E0*I_ESP_I3x3z_G3yz_aa-2.0E0*2*I_ESP_G3xz_G3yz_a-2.0E0*2*I_ESP_Gx3z_G3yz_a+2*2*I_ESP_Dxz_G3yz;
    abcd[iGrid*1350+621] = 4.0E0*I_ESP_I2x3yz_G3yz_aa-2.0E0*1*I_ESP_G3yz_G3yz_a;
    abcd[iGrid*1350+622] = 4.0E0*I_ESP_I2x2y2z_G3yz_aa-2.0E0*1*I_ESP_G2x2y_G3yz_a-2.0E0*1*I_ESP_G2y2z_G3yz_a+1*I_ESP_D2y_G3yz;
    abcd[iGrid*1350+623] = 4.0E0*I_ESP_I2xy3z_G3yz_aa-2.0E0*2*I_ESP_G2xyz_G3yz_a-2.0E0*1*I_ESP_Gy3z_G3yz_a+2*I_ESP_Dyz_G3yz;
    abcd[iGrid*1350+624] = 4.0E0*I_ESP_I2x4z_G3yz_aa-2.0E0*3*I_ESP_G2x2z_G3yz_a-2.0E0*1*I_ESP_G4z_G3yz_a+3*I_ESP_D2z_G3yz;
    abcd[iGrid*1350+625] = 4.0E0*I_ESP_Ix4yz_G3yz_aa;
    abcd[iGrid*1350+626] = 4.0E0*I_ESP_Ix3y2z_G3yz_aa-2.0E0*1*I_ESP_Gx3y_G3yz_a;
    abcd[iGrid*1350+627] = 4.0E0*I_ESP_Ix2y3z_G3yz_aa-2.0E0*2*I_ESP_Gx2yz_G3yz_a;
    abcd[iGrid*1350+628] = 4.0E0*I_ESP_Ixy4z_G3yz_aa-2.0E0*3*I_ESP_Gxy2z_G3yz_a;
    abcd[iGrid*1350+629] = 4.0E0*I_ESP_Ix5z_G3yz_aa-2.0E0*4*I_ESP_Gx3z_G3yz_a;
    abcd[iGrid*1350+630] = 4.0E0*I_ESP_I5xz_G2y2z_aa-2.0E0*4*I_ESP_G3xz_G2y2z_a;
    abcd[iGrid*1350+631] = 4.0E0*I_ESP_I4xyz_G2y2z_aa-2.0E0*3*I_ESP_G2xyz_G2y2z_a;
    abcd[iGrid*1350+632] = 4.0E0*I_ESP_I4x2z_G2y2z_aa-2.0E0*1*I_ESP_G4x_G2y2z_a-2.0E0*3*I_ESP_G2x2z_G2y2z_a+3*1*I_ESP_D2x_G2y2z;
    abcd[iGrid*1350+633] = 4.0E0*I_ESP_I3x2yz_G2y2z_aa-2.0E0*2*I_ESP_Gx2yz_G2y2z_a;
    abcd[iGrid*1350+634] = 4.0E0*I_ESP_I3xy2z_G2y2z_aa-2.0E0*1*I_ESP_G3xy_G2y2z_a-2.0E0*2*I_ESP_Gxy2z_G2y2z_a+2*1*I_ESP_Dxy_G2y2z;
    abcd[iGrid*1350+635] = 4.0E0*I_ESP_I3x3z_G2y2z_aa-2.0E0*2*I_ESP_G3xz_G2y2z_a-2.0E0*2*I_ESP_Gx3z_G2y2z_a+2*2*I_ESP_Dxz_G2y2z;
    abcd[iGrid*1350+636] = 4.0E0*I_ESP_I2x3yz_G2y2z_aa-2.0E0*1*I_ESP_G3yz_G2y2z_a;
    abcd[iGrid*1350+637] = 4.0E0*I_ESP_I2x2y2z_G2y2z_aa-2.0E0*1*I_ESP_G2x2y_G2y2z_a-2.0E0*1*I_ESP_G2y2z_G2y2z_a+1*I_ESP_D2y_G2y2z;
    abcd[iGrid*1350+638] = 4.0E0*I_ESP_I2xy3z_G2y2z_aa-2.0E0*2*I_ESP_G2xyz_G2y2z_a-2.0E0*1*I_ESP_Gy3z_G2y2z_a+2*I_ESP_Dyz_G2y2z;
    abcd[iGrid*1350+639] = 4.0E0*I_ESP_I2x4z_G2y2z_aa-2.0E0*3*I_ESP_G2x2z_G2y2z_a-2.0E0*1*I_ESP_G4z_G2y2z_a+3*I_ESP_D2z_G2y2z;
    abcd[iGrid*1350+640] = 4.0E0*I_ESP_Ix4yz_G2y2z_aa;
    abcd[iGrid*1350+641] = 4.0E0*I_ESP_Ix3y2z_G2y2z_aa-2.0E0*1*I_ESP_Gx3y_G2y2z_a;
    abcd[iGrid*1350+642] = 4.0E0*I_ESP_Ix2y3z_G2y2z_aa-2.0E0*2*I_ESP_Gx2yz_G2y2z_a;
    abcd[iGrid*1350+643] = 4.0E0*I_ESP_Ixy4z_G2y2z_aa-2.0E0*3*I_ESP_Gxy2z_G2y2z_a;
    abcd[iGrid*1350+644] = 4.0E0*I_ESP_Ix5z_G2y2z_aa-2.0E0*4*I_ESP_Gx3z_G2y2z_a;
    abcd[iGrid*1350+645] = 4.0E0*I_ESP_I5xz_Gy3z_aa-2.0E0*4*I_ESP_G3xz_Gy3z_a;
    abcd[iGrid*1350+646] = 4.0E0*I_ESP_I4xyz_Gy3z_aa-2.0E0*3*I_ESP_G2xyz_Gy3z_a;
    abcd[iGrid*1350+647] = 4.0E0*I_ESP_I4x2z_Gy3z_aa-2.0E0*1*I_ESP_G4x_Gy3z_a-2.0E0*3*I_ESP_G2x2z_Gy3z_a+3*1*I_ESP_D2x_Gy3z;
    abcd[iGrid*1350+648] = 4.0E0*I_ESP_I3x2yz_Gy3z_aa-2.0E0*2*I_ESP_Gx2yz_Gy3z_a;
    abcd[iGrid*1350+649] = 4.0E0*I_ESP_I3xy2z_Gy3z_aa-2.0E0*1*I_ESP_G3xy_Gy3z_a-2.0E0*2*I_ESP_Gxy2z_Gy3z_a+2*1*I_ESP_Dxy_Gy3z;
    abcd[iGrid*1350+650] = 4.0E0*I_ESP_I3x3z_Gy3z_aa-2.0E0*2*I_ESP_G3xz_Gy3z_a-2.0E0*2*I_ESP_Gx3z_Gy3z_a+2*2*I_ESP_Dxz_Gy3z;
    abcd[iGrid*1350+651] = 4.0E0*I_ESP_I2x3yz_Gy3z_aa-2.0E0*1*I_ESP_G3yz_Gy3z_a;
    abcd[iGrid*1350+652] = 4.0E0*I_ESP_I2x2y2z_Gy3z_aa-2.0E0*1*I_ESP_G2x2y_Gy3z_a-2.0E0*1*I_ESP_G2y2z_Gy3z_a+1*I_ESP_D2y_Gy3z;
    abcd[iGrid*1350+653] = 4.0E0*I_ESP_I2xy3z_Gy3z_aa-2.0E0*2*I_ESP_G2xyz_Gy3z_a-2.0E0*1*I_ESP_Gy3z_Gy3z_a+2*I_ESP_Dyz_Gy3z;
    abcd[iGrid*1350+654] = 4.0E0*I_ESP_I2x4z_Gy3z_aa-2.0E0*3*I_ESP_G2x2z_Gy3z_a-2.0E0*1*I_ESP_G4z_Gy3z_a+3*I_ESP_D2z_Gy3z;
    abcd[iGrid*1350+655] = 4.0E0*I_ESP_Ix4yz_Gy3z_aa;
    abcd[iGrid*1350+656] = 4.0E0*I_ESP_Ix3y2z_Gy3z_aa-2.0E0*1*I_ESP_Gx3y_Gy3z_a;
    abcd[iGrid*1350+657] = 4.0E0*I_ESP_Ix2y3z_Gy3z_aa-2.0E0*2*I_ESP_Gx2yz_Gy3z_a;
    abcd[iGrid*1350+658] = 4.0E0*I_ESP_Ixy4z_Gy3z_aa-2.0E0*3*I_ESP_Gxy2z_Gy3z_a;
    abcd[iGrid*1350+659] = 4.0E0*I_ESP_Ix5z_Gy3z_aa-2.0E0*4*I_ESP_Gx3z_Gy3z_a;
    abcd[iGrid*1350+660] = 4.0E0*I_ESP_I5xz_G4z_aa-2.0E0*4*I_ESP_G3xz_G4z_a;
    abcd[iGrid*1350+661] = 4.0E0*I_ESP_I4xyz_G4z_aa-2.0E0*3*I_ESP_G2xyz_G4z_a;
    abcd[iGrid*1350+662] = 4.0E0*I_ESP_I4x2z_G4z_aa-2.0E0*1*I_ESP_G4x_G4z_a-2.0E0*3*I_ESP_G2x2z_G4z_a+3*1*I_ESP_D2x_G4z;
    abcd[iGrid*1350+663] = 4.0E0*I_ESP_I3x2yz_G4z_aa-2.0E0*2*I_ESP_Gx2yz_G4z_a;
    abcd[iGrid*1350+664] = 4.0E0*I_ESP_I3xy2z_G4z_aa-2.0E0*1*I_ESP_G3xy_G4z_a-2.0E0*2*I_ESP_Gxy2z_G4z_a+2*1*I_ESP_Dxy_G4z;
    abcd[iGrid*1350+665] = 4.0E0*I_ESP_I3x3z_G4z_aa-2.0E0*2*I_ESP_G3xz_G4z_a-2.0E0*2*I_ESP_Gx3z_G4z_a+2*2*I_ESP_Dxz_G4z;
    abcd[iGrid*1350+666] = 4.0E0*I_ESP_I2x3yz_G4z_aa-2.0E0*1*I_ESP_G3yz_G4z_a;
    abcd[iGrid*1350+667] = 4.0E0*I_ESP_I2x2y2z_G4z_aa-2.0E0*1*I_ESP_G2x2y_G4z_a-2.0E0*1*I_ESP_G2y2z_G4z_a+1*I_ESP_D2y_G4z;
    abcd[iGrid*1350+668] = 4.0E0*I_ESP_I2xy3z_G4z_aa-2.0E0*2*I_ESP_G2xyz_G4z_a-2.0E0*1*I_ESP_Gy3z_G4z_a+2*I_ESP_Dyz_G4z;
    abcd[iGrid*1350+669] = 4.0E0*I_ESP_I2x4z_G4z_aa-2.0E0*3*I_ESP_G2x2z_G4z_a-2.0E0*1*I_ESP_G4z_G4z_a+3*I_ESP_D2z_G4z;
    abcd[iGrid*1350+670] = 4.0E0*I_ESP_Ix4yz_G4z_aa;
    abcd[iGrid*1350+671] = 4.0E0*I_ESP_Ix3y2z_G4z_aa-2.0E0*1*I_ESP_Gx3y_G4z_a;
    abcd[iGrid*1350+672] = 4.0E0*I_ESP_Ix2y3z_G4z_aa-2.0E0*2*I_ESP_Gx2yz_G4z_a;
    abcd[iGrid*1350+673] = 4.0E0*I_ESP_Ixy4z_G4z_aa-2.0E0*3*I_ESP_Gxy2z_G4z_a;
    abcd[iGrid*1350+674] = 4.0E0*I_ESP_Ix5z_G4z_aa-2.0E0*4*I_ESP_Gx3z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G_aa
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_D_G
     ************************************************************/
    abcd[iGrid*1350+675] = 4.0E0*I_ESP_I4x2y_G4x_aa-2.0E0*1*I_ESP_G4x_G4x_a;
    abcd[iGrid*1350+676] = 4.0E0*I_ESP_I3x3y_G4x_aa-2.0E0*1*I_ESP_G3xy_G4x_a-2.0E0*2*I_ESP_G3xy_G4x_a;
    abcd[iGrid*1350+677] = 4.0E0*I_ESP_I3x2yz_G4x_aa-2.0E0*1*I_ESP_G3xz_G4x_a;
    abcd[iGrid*1350+678] = 4.0E0*I_ESP_I2x4y_G4x_aa-2.0E0*2*I_ESP_G2x2y_G4x_a-2.0E0*3*I_ESP_G2x2y_G4x_a+2*1*I_ESP_D2x_G4x;
    abcd[iGrid*1350+679] = 4.0E0*I_ESP_I2x3yz_G4x_aa-2.0E0*1*I_ESP_G2xyz_G4x_a-2.0E0*2*I_ESP_G2xyz_G4x_a;
    abcd[iGrid*1350+680] = 4.0E0*I_ESP_I2x2y2z_G4x_aa-2.0E0*1*I_ESP_G2x2z_G4x_a;
    abcd[iGrid*1350+681] = 4.0E0*I_ESP_Ix5y_G4x_aa-2.0E0*3*I_ESP_Gx3y_G4x_a-2.0E0*4*I_ESP_Gx3y_G4x_a+3*2*I_ESP_Dxy_G4x;
    abcd[iGrid*1350+682] = 4.0E0*I_ESP_Ix4yz_G4x_aa-2.0E0*2*I_ESP_Gx2yz_G4x_a-2.0E0*3*I_ESP_Gx2yz_G4x_a+2*1*I_ESP_Dxz_G4x;
    abcd[iGrid*1350+683] = 4.0E0*I_ESP_Ix3y2z_G4x_aa-2.0E0*1*I_ESP_Gxy2z_G4x_a-2.0E0*2*I_ESP_Gxy2z_G4x_a;
    abcd[iGrid*1350+684] = 4.0E0*I_ESP_Ix2y3z_G4x_aa-2.0E0*1*I_ESP_Gx3z_G4x_a;
    abcd[iGrid*1350+685] = 4.0E0*I_ESP_I6y_G4x_aa-2.0E0*4*I_ESP_G4y_G4x_a-2.0E0*5*I_ESP_G4y_G4x_a+4*3*I_ESP_D2y_G4x;
    abcd[iGrid*1350+686] = 4.0E0*I_ESP_I5yz_G4x_aa-2.0E0*3*I_ESP_G3yz_G4x_a-2.0E0*4*I_ESP_G3yz_G4x_a+3*2*I_ESP_Dyz_G4x;
    abcd[iGrid*1350+687] = 4.0E0*I_ESP_I4y2z_G4x_aa-2.0E0*2*I_ESP_G2y2z_G4x_a-2.0E0*3*I_ESP_G2y2z_G4x_a+2*1*I_ESP_D2z_G4x;
    abcd[iGrid*1350+688] = 4.0E0*I_ESP_I3y3z_G4x_aa-2.0E0*1*I_ESP_Gy3z_G4x_a-2.0E0*2*I_ESP_Gy3z_G4x_a;
    abcd[iGrid*1350+689] = 4.0E0*I_ESP_I2y4z_G4x_aa-2.0E0*1*I_ESP_G4z_G4x_a;
    abcd[iGrid*1350+690] = 4.0E0*I_ESP_I4x2y_G3xy_aa-2.0E0*1*I_ESP_G4x_G3xy_a;
    abcd[iGrid*1350+691] = 4.0E0*I_ESP_I3x3y_G3xy_aa-2.0E0*1*I_ESP_G3xy_G3xy_a-2.0E0*2*I_ESP_G3xy_G3xy_a;
    abcd[iGrid*1350+692] = 4.0E0*I_ESP_I3x2yz_G3xy_aa-2.0E0*1*I_ESP_G3xz_G3xy_a;
    abcd[iGrid*1350+693] = 4.0E0*I_ESP_I2x4y_G3xy_aa-2.0E0*2*I_ESP_G2x2y_G3xy_a-2.0E0*3*I_ESP_G2x2y_G3xy_a+2*1*I_ESP_D2x_G3xy;
    abcd[iGrid*1350+694] = 4.0E0*I_ESP_I2x3yz_G3xy_aa-2.0E0*1*I_ESP_G2xyz_G3xy_a-2.0E0*2*I_ESP_G2xyz_G3xy_a;
    abcd[iGrid*1350+695] = 4.0E0*I_ESP_I2x2y2z_G3xy_aa-2.0E0*1*I_ESP_G2x2z_G3xy_a;
    abcd[iGrid*1350+696] = 4.0E0*I_ESP_Ix5y_G3xy_aa-2.0E0*3*I_ESP_Gx3y_G3xy_a-2.0E0*4*I_ESP_Gx3y_G3xy_a+3*2*I_ESP_Dxy_G3xy;
    abcd[iGrid*1350+697] = 4.0E0*I_ESP_Ix4yz_G3xy_aa-2.0E0*2*I_ESP_Gx2yz_G3xy_a-2.0E0*3*I_ESP_Gx2yz_G3xy_a+2*1*I_ESP_Dxz_G3xy;
    abcd[iGrid*1350+698] = 4.0E0*I_ESP_Ix3y2z_G3xy_aa-2.0E0*1*I_ESP_Gxy2z_G3xy_a-2.0E0*2*I_ESP_Gxy2z_G3xy_a;
    abcd[iGrid*1350+699] = 4.0E0*I_ESP_Ix2y3z_G3xy_aa-2.0E0*1*I_ESP_Gx3z_G3xy_a;
    abcd[iGrid*1350+700] = 4.0E0*I_ESP_I6y_G3xy_aa-2.0E0*4*I_ESP_G4y_G3xy_a-2.0E0*5*I_ESP_G4y_G3xy_a+4*3*I_ESP_D2y_G3xy;
    abcd[iGrid*1350+701] = 4.0E0*I_ESP_I5yz_G3xy_aa-2.0E0*3*I_ESP_G3yz_G3xy_a-2.0E0*4*I_ESP_G3yz_G3xy_a+3*2*I_ESP_Dyz_G3xy;
    abcd[iGrid*1350+702] = 4.0E0*I_ESP_I4y2z_G3xy_aa-2.0E0*2*I_ESP_G2y2z_G3xy_a-2.0E0*3*I_ESP_G2y2z_G3xy_a+2*1*I_ESP_D2z_G3xy;
    abcd[iGrid*1350+703] = 4.0E0*I_ESP_I3y3z_G3xy_aa-2.0E0*1*I_ESP_Gy3z_G3xy_a-2.0E0*2*I_ESP_Gy3z_G3xy_a;
    abcd[iGrid*1350+704] = 4.0E0*I_ESP_I2y4z_G3xy_aa-2.0E0*1*I_ESP_G4z_G3xy_a;
    abcd[iGrid*1350+705] = 4.0E0*I_ESP_I4x2y_G3xz_aa-2.0E0*1*I_ESP_G4x_G3xz_a;
    abcd[iGrid*1350+706] = 4.0E0*I_ESP_I3x3y_G3xz_aa-2.0E0*1*I_ESP_G3xy_G3xz_a-2.0E0*2*I_ESP_G3xy_G3xz_a;
    abcd[iGrid*1350+707] = 4.0E0*I_ESP_I3x2yz_G3xz_aa-2.0E0*1*I_ESP_G3xz_G3xz_a;
    abcd[iGrid*1350+708] = 4.0E0*I_ESP_I2x4y_G3xz_aa-2.0E0*2*I_ESP_G2x2y_G3xz_a-2.0E0*3*I_ESP_G2x2y_G3xz_a+2*1*I_ESP_D2x_G3xz;
    abcd[iGrid*1350+709] = 4.0E0*I_ESP_I2x3yz_G3xz_aa-2.0E0*1*I_ESP_G2xyz_G3xz_a-2.0E0*2*I_ESP_G2xyz_G3xz_a;
    abcd[iGrid*1350+710] = 4.0E0*I_ESP_I2x2y2z_G3xz_aa-2.0E0*1*I_ESP_G2x2z_G3xz_a;
    abcd[iGrid*1350+711] = 4.0E0*I_ESP_Ix5y_G3xz_aa-2.0E0*3*I_ESP_Gx3y_G3xz_a-2.0E0*4*I_ESP_Gx3y_G3xz_a+3*2*I_ESP_Dxy_G3xz;
    abcd[iGrid*1350+712] = 4.0E0*I_ESP_Ix4yz_G3xz_aa-2.0E0*2*I_ESP_Gx2yz_G3xz_a-2.0E0*3*I_ESP_Gx2yz_G3xz_a+2*1*I_ESP_Dxz_G3xz;
    abcd[iGrid*1350+713] = 4.0E0*I_ESP_Ix3y2z_G3xz_aa-2.0E0*1*I_ESP_Gxy2z_G3xz_a-2.0E0*2*I_ESP_Gxy2z_G3xz_a;
    abcd[iGrid*1350+714] = 4.0E0*I_ESP_Ix2y3z_G3xz_aa-2.0E0*1*I_ESP_Gx3z_G3xz_a;
    abcd[iGrid*1350+715] = 4.0E0*I_ESP_I6y_G3xz_aa-2.0E0*4*I_ESP_G4y_G3xz_a-2.0E0*5*I_ESP_G4y_G3xz_a+4*3*I_ESP_D2y_G3xz;
    abcd[iGrid*1350+716] = 4.0E0*I_ESP_I5yz_G3xz_aa-2.0E0*3*I_ESP_G3yz_G3xz_a-2.0E0*4*I_ESP_G3yz_G3xz_a+3*2*I_ESP_Dyz_G3xz;
    abcd[iGrid*1350+717] = 4.0E0*I_ESP_I4y2z_G3xz_aa-2.0E0*2*I_ESP_G2y2z_G3xz_a-2.0E0*3*I_ESP_G2y2z_G3xz_a+2*1*I_ESP_D2z_G3xz;
    abcd[iGrid*1350+718] = 4.0E0*I_ESP_I3y3z_G3xz_aa-2.0E0*1*I_ESP_Gy3z_G3xz_a-2.0E0*2*I_ESP_Gy3z_G3xz_a;
    abcd[iGrid*1350+719] = 4.0E0*I_ESP_I2y4z_G3xz_aa-2.0E0*1*I_ESP_G4z_G3xz_a;
    abcd[iGrid*1350+720] = 4.0E0*I_ESP_I4x2y_G2x2y_aa-2.0E0*1*I_ESP_G4x_G2x2y_a;
    abcd[iGrid*1350+721] = 4.0E0*I_ESP_I3x3y_G2x2y_aa-2.0E0*1*I_ESP_G3xy_G2x2y_a-2.0E0*2*I_ESP_G3xy_G2x2y_a;
    abcd[iGrid*1350+722] = 4.0E0*I_ESP_I3x2yz_G2x2y_aa-2.0E0*1*I_ESP_G3xz_G2x2y_a;
    abcd[iGrid*1350+723] = 4.0E0*I_ESP_I2x4y_G2x2y_aa-2.0E0*2*I_ESP_G2x2y_G2x2y_a-2.0E0*3*I_ESP_G2x2y_G2x2y_a+2*1*I_ESP_D2x_G2x2y;
    abcd[iGrid*1350+724] = 4.0E0*I_ESP_I2x3yz_G2x2y_aa-2.0E0*1*I_ESP_G2xyz_G2x2y_a-2.0E0*2*I_ESP_G2xyz_G2x2y_a;
    abcd[iGrid*1350+725] = 4.0E0*I_ESP_I2x2y2z_G2x2y_aa-2.0E0*1*I_ESP_G2x2z_G2x2y_a;
    abcd[iGrid*1350+726] = 4.0E0*I_ESP_Ix5y_G2x2y_aa-2.0E0*3*I_ESP_Gx3y_G2x2y_a-2.0E0*4*I_ESP_Gx3y_G2x2y_a+3*2*I_ESP_Dxy_G2x2y;
    abcd[iGrid*1350+727] = 4.0E0*I_ESP_Ix4yz_G2x2y_aa-2.0E0*2*I_ESP_Gx2yz_G2x2y_a-2.0E0*3*I_ESP_Gx2yz_G2x2y_a+2*1*I_ESP_Dxz_G2x2y;
    abcd[iGrid*1350+728] = 4.0E0*I_ESP_Ix3y2z_G2x2y_aa-2.0E0*1*I_ESP_Gxy2z_G2x2y_a-2.0E0*2*I_ESP_Gxy2z_G2x2y_a;
    abcd[iGrid*1350+729] = 4.0E0*I_ESP_Ix2y3z_G2x2y_aa-2.0E0*1*I_ESP_Gx3z_G2x2y_a;
    abcd[iGrid*1350+730] = 4.0E0*I_ESP_I6y_G2x2y_aa-2.0E0*4*I_ESP_G4y_G2x2y_a-2.0E0*5*I_ESP_G4y_G2x2y_a+4*3*I_ESP_D2y_G2x2y;
    abcd[iGrid*1350+731] = 4.0E0*I_ESP_I5yz_G2x2y_aa-2.0E0*3*I_ESP_G3yz_G2x2y_a-2.0E0*4*I_ESP_G3yz_G2x2y_a+3*2*I_ESP_Dyz_G2x2y;
    abcd[iGrid*1350+732] = 4.0E0*I_ESP_I4y2z_G2x2y_aa-2.0E0*2*I_ESP_G2y2z_G2x2y_a-2.0E0*3*I_ESP_G2y2z_G2x2y_a+2*1*I_ESP_D2z_G2x2y;
    abcd[iGrid*1350+733] = 4.0E0*I_ESP_I3y3z_G2x2y_aa-2.0E0*1*I_ESP_Gy3z_G2x2y_a-2.0E0*2*I_ESP_Gy3z_G2x2y_a;
    abcd[iGrid*1350+734] = 4.0E0*I_ESP_I2y4z_G2x2y_aa-2.0E0*1*I_ESP_G4z_G2x2y_a;
    abcd[iGrid*1350+735] = 4.0E0*I_ESP_I4x2y_G2xyz_aa-2.0E0*1*I_ESP_G4x_G2xyz_a;
    abcd[iGrid*1350+736] = 4.0E0*I_ESP_I3x3y_G2xyz_aa-2.0E0*1*I_ESP_G3xy_G2xyz_a-2.0E0*2*I_ESP_G3xy_G2xyz_a;
    abcd[iGrid*1350+737] = 4.0E0*I_ESP_I3x2yz_G2xyz_aa-2.0E0*1*I_ESP_G3xz_G2xyz_a;
    abcd[iGrid*1350+738] = 4.0E0*I_ESP_I2x4y_G2xyz_aa-2.0E0*2*I_ESP_G2x2y_G2xyz_a-2.0E0*3*I_ESP_G2x2y_G2xyz_a+2*1*I_ESP_D2x_G2xyz;
    abcd[iGrid*1350+739] = 4.0E0*I_ESP_I2x3yz_G2xyz_aa-2.0E0*1*I_ESP_G2xyz_G2xyz_a-2.0E0*2*I_ESP_G2xyz_G2xyz_a;
    abcd[iGrid*1350+740] = 4.0E0*I_ESP_I2x2y2z_G2xyz_aa-2.0E0*1*I_ESP_G2x2z_G2xyz_a;
    abcd[iGrid*1350+741] = 4.0E0*I_ESP_Ix5y_G2xyz_aa-2.0E0*3*I_ESP_Gx3y_G2xyz_a-2.0E0*4*I_ESP_Gx3y_G2xyz_a+3*2*I_ESP_Dxy_G2xyz;
    abcd[iGrid*1350+742] = 4.0E0*I_ESP_Ix4yz_G2xyz_aa-2.0E0*2*I_ESP_Gx2yz_G2xyz_a-2.0E0*3*I_ESP_Gx2yz_G2xyz_a+2*1*I_ESP_Dxz_G2xyz;
    abcd[iGrid*1350+743] = 4.0E0*I_ESP_Ix3y2z_G2xyz_aa-2.0E0*1*I_ESP_Gxy2z_G2xyz_a-2.0E0*2*I_ESP_Gxy2z_G2xyz_a;
    abcd[iGrid*1350+744] = 4.0E0*I_ESP_Ix2y3z_G2xyz_aa-2.0E0*1*I_ESP_Gx3z_G2xyz_a;
    abcd[iGrid*1350+745] = 4.0E0*I_ESP_I6y_G2xyz_aa-2.0E0*4*I_ESP_G4y_G2xyz_a-2.0E0*5*I_ESP_G4y_G2xyz_a+4*3*I_ESP_D2y_G2xyz;
    abcd[iGrid*1350+746] = 4.0E0*I_ESP_I5yz_G2xyz_aa-2.0E0*3*I_ESP_G3yz_G2xyz_a-2.0E0*4*I_ESP_G3yz_G2xyz_a+3*2*I_ESP_Dyz_G2xyz;
    abcd[iGrid*1350+747] = 4.0E0*I_ESP_I4y2z_G2xyz_aa-2.0E0*2*I_ESP_G2y2z_G2xyz_a-2.0E0*3*I_ESP_G2y2z_G2xyz_a+2*1*I_ESP_D2z_G2xyz;
    abcd[iGrid*1350+748] = 4.0E0*I_ESP_I3y3z_G2xyz_aa-2.0E0*1*I_ESP_Gy3z_G2xyz_a-2.0E0*2*I_ESP_Gy3z_G2xyz_a;
    abcd[iGrid*1350+749] = 4.0E0*I_ESP_I2y4z_G2xyz_aa-2.0E0*1*I_ESP_G4z_G2xyz_a;
    abcd[iGrid*1350+750] = 4.0E0*I_ESP_I4x2y_G2x2z_aa-2.0E0*1*I_ESP_G4x_G2x2z_a;
    abcd[iGrid*1350+751] = 4.0E0*I_ESP_I3x3y_G2x2z_aa-2.0E0*1*I_ESP_G3xy_G2x2z_a-2.0E0*2*I_ESP_G3xy_G2x2z_a;
    abcd[iGrid*1350+752] = 4.0E0*I_ESP_I3x2yz_G2x2z_aa-2.0E0*1*I_ESP_G3xz_G2x2z_a;
    abcd[iGrid*1350+753] = 4.0E0*I_ESP_I2x4y_G2x2z_aa-2.0E0*2*I_ESP_G2x2y_G2x2z_a-2.0E0*3*I_ESP_G2x2y_G2x2z_a+2*1*I_ESP_D2x_G2x2z;
    abcd[iGrid*1350+754] = 4.0E0*I_ESP_I2x3yz_G2x2z_aa-2.0E0*1*I_ESP_G2xyz_G2x2z_a-2.0E0*2*I_ESP_G2xyz_G2x2z_a;
    abcd[iGrid*1350+755] = 4.0E0*I_ESP_I2x2y2z_G2x2z_aa-2.0E0*1*I_ESP_G2x2z_G2x2z_a;
    abcd[iGrid*1350+756] = 4.0E0*I_ESP_Ix5y_G2x2z_aa-2.0E0*3*I_ESP_Gx3y_G2x2z_a-2.0E0*4*I_ESP_Gx3y_G2x2z_a+3*2*I_ESP_Dxy_G2x2z;
    abcd[iGrid*1350+757] = 4.0E0*I_ESP_Ix4yz_G2x2z_aa-2.0E0*2*I_ESP_Gx2yz_G2x2z_a-2.0E0*3*I_ESP_Gx2yz_G2x2z_a+2*1*I_ESP_Dxz_G2x2z;
    abcd[iGrid*1350+758] = 4.0E0*I_ESP_Ix3y2z_G2x2z_aa-2.0E0*1*I_ESP_Gxy2z_G2x2z_a-2.0E0*2*I_ESP_Gxy2z_G2x2z_a;
    abcd[iGrid*1350+759] = 4.0E0*I_ESP_Ix2y3z_G2x2z_aa-2.0E0*1*I_ESP_Gx3z_G2x2z_a;
    abcd[iGrid*1350+760] = 4.0E0*I_ESP_I6y_G2x2z_aa-2.0E0*4*I_ESP_G4y_G2x2z_a-2.0E0*5*I_ESP_G4y_G2x2z_a+4*3*I_ESP_D2y_G2x2z;
    abcd[iGrid*1350+761] = 4.0E0*I_ESP_I5yz_G2x2z_aa-2.0E0*3*I_ESP_G3yz_G2x2z_a-2.0E0*4*I_ESP_G3yz_G2x2z_a+3*2*I_ESP_Dyz_G2x2z;
    abcd[iGrid*1350+762] = 4.0E0*I_ESP_I4y2z_G2x2z_aa-2.0E0*2*I_ESP_G2y2z_G2x2z_a-2.0E0*3*I_ESP_G2y2z_G2x2z_a+2*1*I_ESP_D2z_G2x2z;
    abcd[iGrid*1350+763] = 4.0E0*I_ESP_I3y3z_G2x2z_aa-2.0E0*1*I_ESP_Gy3z_G2x2z_a-2.0E0*2*I_ESP_Gy3z_G2x2z_a;
    abcd[iGrid*1350+764] = 4.0E0*I_ESP_I2y4z_G2x2z_aa-2.0E0*1*I_ESP_G4z_G2x2z_a;
    abcd[iGrid*1350+765] = 4.0E0*I_ESP_I4x2y_Gx3y_aa-2.0E0*1*I_ESP_G4x_Gx3y_a;
    abcd[iGrid*1350+766] = 4.0E0*I_ESP_I3x3y_Gx3y_aa-2.0E0*1*I_ESP_G3xy_Gx3y_a-2.0E0*2*I_ESP_G3xy_Gx3y_a;
    abcd[iGrid*1350+767] = 4.0E0*I_ESP_I3x2yz_Gx3y_aa-2.0E0*1*I_ESP_G3xz_Gx3y_a;
    abcd[iGrid*1350+768] = 4.0E0*I_ESP_I2x4y_Gx3y_aa-2.0E0*2*I_ESP_G2x2y_Gx3y_a-2.0E0*3*I_ESP_G2x2y_Gx3y_a+2*1*I_ESP_D2x_Gx3y;
    abcd[iGrid*1350+769] = 4.0E0*I_ESP_I2x3yz_Gx3y_aa-2.0E0*1*I_ESP_G2xyz_Gx3y_a-2.0E0*2*I_ESP_G2xyz_Gx3y_a;
    abcd[iGrid*1350+770] = 4.0E0*I_ESP_I2x2y2z_Gx3y_aa-2.0E0*1*I_ESP_G2x2z_Gx3y_a;
    abcd[iGrid*1350+771] = 4.0E0*I_ESP_Ix5y_Gx3y_aa-2.0E0*3*I_ESP_Gx3y_Gx3y_a-2.0E0*4*I_ESP_Gx3y_Gx3y_a+3*2*I_ESP_Dxy_Gx3y;
    abcd[iGrid*1350+772] = 4.0E0*I_ESP_Ix4yz_Gx3y_aa-2.0E0*2*I_ESP_Gx2yz_Gx3y_a-2.0E0*3*I_ESP_Gx2yz_Gx3y_a+2*1*I_ESP_Dxz_Gx3y;
    abcd[iGrid*1350+773] = 4.0E0*I_ESP_Ix3y2z_Gx3y_aa-2.0E0*1*I_ESP_Gxy2z_Gx3y_a-2.0E0*2*I_ESP_Gxy2z_Gx3y_a;
    abcd[iGrid*1350+774] = 4.0E0*I_ESP_Ix2y3z_Gx3y_aa-2.0E0*1*I_ESP_Gx3z_Gx3y_a;
    abcd[iGrid*1350+775] = 4.0E0*I_ESP_I6y_Gx3y_aa-2.0E0*4*I_ESP_G4y_Gx3y_a-2.0E0*5*I_ESP_G4y_Gx3y_a+4*3*I_ESP_D2y_Gx3y;
    abcd[iGrid*1350+776] = 4.0E0*I_ESP_I5yz_Gx3y_aa-2.0E0*3*I_ESP_G3yz_Gx3y_a-2.0E0*4*I_ESP_G3yz_Gx3y_a+3*2*I_ESP_Dyz_Gx3y;
    abcd[iGrid*1350+777] = 4.0E0*I_ESP_I4y2z_Gx3y_aa-2.0E0*2*I_ESP_G2y2z_Gx3y_a-2.0E0*3*I_ESP_G2y2z_Gx3y_a+2*1*I_ESP_D2z_Gx3y;
    abcd[iGrid*1350+778] = 4.0E0*I_ESP_I3y3z_Gx3y_aa-2.0E0*1*I_ESP_Gy3z_Gx3y_a-2.0E0*2*I_ESP_Gy3z_Gx3y_a;
    abcd[iGrid*1350+779] = 4.0E0*I_ESP_I2y4z_Gx3y_aa-2.0E0*1*I_ESP_G4z_Gx3y_a;
    abcd[iGrid*1350+780] = 4.0E0*I_ESP_I4x2y_Gx2yz_aa-2.0E0*1*I_ESP_G4x_Gx2yz_a;
    abcd[iGrid*1350+781] = 4.0E0*I_ESP_I3x3y_Gx2yz_aa-2.0E0*1*I_ESP_G3xy_Gx2yz_a-2.0E0*2*I_ESP_G3xy_Gx2yz_a;
    abcd[iGrid*1350+782] = 4.0E0*I_ESP_I3x2yz_Gx2yz_aa-2.0E0*1*I_ESP_G3xz_Gx2yz_a;
    abcd[iGrid*1350+783] = 4.0E0*I_ESP_I2x4y_Gx2yz_aa-2.0E0*2*I_ESP_G2x2y_Gx2yz_a-2.0E0*3*I_ESP_G2x2y_Gx2yz_a+2*1*I_ESP_D2x_Gx2yz;
    abcd[iGrid*1350+784] = 4.0E0*I_ESP_I2x3yz_Gx2yz_aa-2.0E0*1*I_ESP_G2xyz_Gx2yz_a-2.0E0*2*I_ESP_G2xyz_Gx2yz_a;
    abcd[iGrid*1350+785] = 4.0E0*I_ESP_I2x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_G2x2z_Gx2yz_a;
    abcd[iGrid*1350+786] = 4.0E0*I_ESP_Ix5y_Gx2yz_aa-2.0E0*3*I_ESP_Gx3y_Gx2yz_a-2.0E0*4*I_ESP_Gx3y_Gx2yz_a+3*2*I_ESP_Dxy_Gx2yz;
    abcd[iGrid*1350+787] = 4.0E0*I_ESP_Ix4yz_Gx2yz_aa-2.0E0*2*I_ESP_Gx2yz_Gx2yz_a-2.0E0*3*I_ESP_Gx2yz_Gx2yz_a+2*1*I_ESP_Dxz_Gx2yz;
    abcd[iGrid*1350+788] = 4.0E0*I_ESP_Ix3y2z_Gx2yz_aa-2.0E0*1*I_ESP_Gxy2z_Gx2yz_a-2.0E0*2*I_ESP_Gxy2z_Gx2yz_a;
    abcd[iGrid*1350+789] = 4.0E0*I_ESP_Ix2y3z_Gx2yz_aa-2.0E0*1*I_ESP_Gx3z_Gx2yz_a;
    abcd[iGrid*1350+790] = 4.0E0*I_ESP_I6y_Gx2yz_aa-2.0E0*4*I_ESP_G4y_Gx2yz_a-2.0E0*5*I_ESP_G4y_Gx2yz_a+4*3*I_ESP_D2y_Gx2yz;
    abcd[iGrid*1350+791] = 4.0E0*I_ESP_I5yz_Gx2yz_aa-2.0E0*3*I_ESP_G3yz_Gx2yz_a-2.0E0*4*I_ESP_G3yz_Gx2yz_a+3*2*I_ESP_Dyz_Gx2yz;
    abcd[iGrid*1350+792] = 4.0E0*I_ESP_I4y2z_Gx2yz_aa-2.0E0*2*I_ESP_G2y2z_Gx2yz_a-2.0E0*3*I_ESP_G2y2z_Gx2yz_a+2*1*I_ESP_D2z_Gx2yz;
    abcd[iGrid*1350+793] = 4.0E0*I_ESP_I3y3z_Gx2yz_aa-2.0E0*1*I_ESP_Gy3z_Gx2yz_a-2.0E0*2*I_ESP_Gy3z_Gx2yz_a;
    abcd[iGrid*1350+794] = 4.0E0*I_ESP_I2y4z_Gx2yz_aa-2.0E0*1*I_ESP_G4z_Gx2yz_a;
    abcd[iGrid*1350+795] = 4.0E0*I_ESP_I4x2y_Gxy2z_aa-2.0E0*1*I_ESP_G4x_Gxy2z_a;
    abcd[iGrid*1350+796] = 4.0E0*I_ESP_I3x3y_Gxy2z_aa-2.0E0*1*I_ESP_G3xy_Gxy2z_a-2.0E0*2*I_ESP_G3xy_Gxy2z_a;
    abcd[iGrid*1350+797] = 4.0E0*I_ESP_I3x2yz_Gxy2z_aa-2.0E0*1*I_ESP_G3xz_Gxy2z_a;
    abcd[iGrid*1350+798] = 4.0E0*I_ESP_I2x4y_Gxy2z_aa-2.0E0*2*I_ESP_G2x2y_Gxy2z_a-2.0E0*3*I_ESP_G2x2y_Gxy2z_a+2*1*I_ESP_D2x_Gxy2z;
    abcd[iGrid*1350+799] = 4.0E0*I_ESP_I2x3yz_Gxy2z_aa-2.0E0*1*I_ESP_G2xyz_Gxy2z_a-2.0E0*2*I_ESP_G2xyz_Gxy2z_a;
    abcd[iGrid*1350+800] = 4.0E0*I_ESP_I2x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_G2x2z_Gxy2z_a;
    abcd[iGrid*1350+801] = 4.0E0*I_ESP_Ix5y_Gxy2z_aa-2.0E0*3*I_ESP_Gx3y_Gxy2z_a-2.0E0*4*I_ESP_Gx3y_Gxy2z_a+3*2*I_ESP_Dxy_Gxy2z;
    abcd[iGrid*1350+802] = 4.0E0*I_ESP_Ix4yz_Gxy2z_aa-2.0E0*2*I_ESP_Gx2yz_Gxy2z_a-2.0E0*3*I_ESP_Gx2yz_Gxy2z_a+2*1*I_ESP_Dxz_Gxy2z;
    abcd[iGrid*1350+803] = 4.0E0*I_ESP_Ix3y2z_Gxy2z_aa-2.0E0*1*I_ESP_Gxy2z_Gxy2z_a-2.0E0*2*I_ESP_Gxy2z_Gxy2z_a;
    abcd[iGrid*1350+804] = 4.0E0*I_ESP_Ix2y3z_Gxy2z_aa-2.0E0*1*I_ESP_Gx3z_Gxy2z_a;
    abcd[iGrid*1350+805] = 4.0E0*I_ESP_I6y_Gxy2z_aa-2.0E0*4*I_ESP_G4y_Gxy2z_a-2.0E0*5*I_ESP_G4y_Gxy2z_a+4*3*I_ESP_D2y_Gxy2z;
    abcd[iGrid*1350+806] = 4.0E0*I_ESP_I5yz_Gxy2z_aa-2.0E0*3*I_ESP_G3yz_Gxy2z_a-2.0E0*4*I_ESP_G3yz_Gxy2z_a+3*2*I_ESP_Dyz_Gxy2z;
    abcd[iGrid*1350+807] = 4.0E0*I_ESP_I4y2z_Gxy2z_aa-2.0E0*2*I_ESP_G2y2z_Gxy2z_a-2.0E0*3*I_ESP_G2y2z_Gxy2z_a+2*1*I_ESP_D2z_Gxy2z;
    abcd[iGrid*1350+808] = 4.0E0*I_ESP_I3y3z_Gxy2z_aa-2.0E0*1*I_ESP_Gy3z_Gxy2z_a-2.0E0*2*I_ESP_Gy3z_Gxy2z_a;
    abcd[iGrid*1350+809] = 4.0E0*I_ESP_I2y4z_Gxy2z_aa-2.0E0*1*I_ESP_G4z_Gxy2z_a;
    abcd[iGrid*1350+810] = 4.0E0*I_ESP_I4x2y_Gx3z_aa-2.0E0*1*I_ESP_G4x_Gx3z_a;
    abcd[iGrid*1350+811] = 4.0E0*I_ESP_I3x3y_Gx3z_aa-2.0E0*1*I_ESP_G3xy_Gx3z_a-2.0E0*2*I_ESP_G3xy_Gx3z_a;
    abcd[iGrid*1350+812] = 4.0E0*I_ESP_I3x2yz_Gx3z_aa-2.0E0*1*I_ESP_G3xz_Gx3z_a;
    abcd[iGrid*1350+813] = 4.0E0*I_ESP_I2x4y_Gx3z_aa-2.0E0*2*I_ESP_G2x2y_Gx3z_a-2.0E0*3*I_ESP_G2x2y_Gx3z_a+2*1*I_ESP_D2x_Gx3z;
    abcd[iGrid*1350+814] = 4.0E0*I_ESP_I2x3yz_Gx3z_aa-2.0E0*1*I_ESP_G2xyz_Gx3z_a-2.0E0*2*I_ESP_G2xyz_Gx3z_a;
    abcd[iGrid*1350+815] = 4.0E0*I_ESP_I2x2y2z_Gx3z_aa-2.0E0*1*I_ESP_G2x2z_Gx3z_a;
    abcd[iGrid*1350+816] = 4.0E0*I_ESP_Ix5y_Gx3z_aa-2.0E0*3*I_ESP_Gx3y_Gx3z_a-2.0E0*4*I_ESP_Gx3y_Gx3z_a+3*2*I_ESP_Dxy_Gx3z;
    abcd[iGrid*1350+817] = 4.0E0*I_ESP_Ix4yz_Gx3z_aa-2.0E0*2*I_ESP_Gx2yz_Gx3z_a-2.0E0*3*I_ESP_Gx2yz_Gx3z_a+2*1*I_ESP_Dxz_Gx3z;
    abcd[iGrid*1350+818] = 4.0E0*I_ESP_Ix3y2z_Gx3z_aa-2.0E0*1*I_ESP_Gxy2z_Gx3z_a-2.0E0*2*I_ESP_Gxy2z_Gx3z_a;
    abcd[iGrid*1350+819] = 4.0E0*I_ESP_Ix2y3z_Gx3z_aa-2.0E0*1*I_ESP_Gx3z_Gx3z_a;
    abcd[iGrid*1350+820] = 4.0E0*I_ESP_I6y_Gx3z_aa-2.0E0*4*I_ESP_G4y_Gx3z_a-2.0E0*5*I_ESP_G4y_Gx3z_a+4*3*I_ESP_D2y_Gx3z;
    abcd[iGrid*1350+821] = 4.0E0*I_ESP_I5yz_Gx3z_aa-2.0E0*3*I_ESP_G3yz_Gx3z_a-2.0E0*4*I_ESP_G3yz_Gx3z_a+3*2*I_ESP_Dyz_Gx3z;
    abcd[iGrid*1350+822] = 4.0E0*I_ESP_I4y2z_Gx3z_aa-2.0E0*2*I_ESP_G2y2z_Gx3z_a-2.0E0*3*I_ESP_G2y2z_Gx3z_a+2*1*I_ESP_D2z_Gx3z;
    abcd[iGrid*1350+823] = 4.0E0*I_ESP_I3y3z_Gx3z_aa-2.0E0*1*I_ESP_Gy3z_Gx3z_a-2.0E0*2*I_ESP_Gy3z_Gx3z_a;
    abcd[iGrid*1350+824] = 4.0E0*I_ESP_I2y4z_Gx3z_aa-2.0E0*1*I_ESP_G4z_Gx3z_a;
    abcd[iGrid*1350+825] = 4.0E0*I_ESP_I4x2y_G4y_aa-2.0E0*1*I_ESP_G4x_G4y_a;
    abcd[iGrid*1350+826] = 4.0E0*I_ESP_I3x3y_G4y_aa-2.0E0*1*I_ESP_G3xy_G4y_a-2.0E0*2*I_ESP_G3xy_G4y_a;
    abcd[iGrid*1350+827] = 4.0E0*I_ESP_I3x2yz_G4y_aa-2.0E0*1*I_ESP_G3xz_G4y_a;
    abcd[iGrid*1350+828] = 4.0E0*I_ESP_I2x4y_G4y_aa-2.0E0*2*I_ESP_G2x2y_G4y_a-2.0E0*3*I_ESP_G2x2y_G4y_a+2*1*I_ESP_D2x_G4y;
    abcd[iGrid*1350+829] = 4.0E0*I_ESP_I2x3yz_G4y_aa-2.0E0*1*I_ESP_G2xyz_G4y_a-2.0E0*2*I_ESP_G2xyz_G4y_a;
    abcd[iGrid*1350+830] = 4.0E0*I_ESP_I2x2y2z_G4y_aa-2.0E0*1*I_ESP_G2x2z_G4y_a;
    abcd[iGrid*1350+831] = 4.0E0*I_ESP_Ix5y_G4y_aa-2.0E0*3*I_ESP_Gx3y_G4y_a-2.0E0*4*I_ESP_Gx3y_G4y_a+3*2*I_ESP_Dxy_G4y;
    abcd[iGrid*1350+832] = 4.0E0*I_ESP_Ix4yz_G4y_aa-2.0E0*2*I_ESP_Gx2yz_G4y_a-2.0E0*3*I_ESP_Gx2yz_G4y_a+2*1*I_ESP_Dxz_G4y;
    abcd[iGrid*1350+833] = 4.0E0*I_ESP_Ix3y2z_G4y_aa-2.0E0*1*I_ESP_Gxy2z_G4y_a-2.0E0*2*I_ESP_Gxy2z_G4y_a;
    abcd[iGrid*1350+834] = 4.0E0*I_ESP_Ix2y3z_G4y_aa-2.0E0*1*I_ESP_Gx3z_G4y_a;
    abcd[iGrid*1350+835] = 4.0E0*I_ESP_I6y_G4y_aa-2.0E0*4*I_ESP_G4y_G4y_a-2.0E0*5*I_ESP_G4y_G4y_a+4*3*I_ESP_D2y_G4y;
    abcd[iGrid*1350+836] = 4.0E0*I_ESP_I5yz_G4y_aa-2.0E0*3*I_ESP_G3yz_G4y_a-2.0E0*4*I_ESP_G3yz_G4y_a+3*2*I_ESP_Dyz_G4y;
    abcd[iGrid*1350+837] = 4.0E0*I_ESP_I4y2z_G4y_aa-2.0E0*2*I_ESP_G2y2z_G4y_a-2.0E0*3*I_ESP_G2y2z_G4y_a+2*1*I_ESP_D2z_G4y;
    abcd[iGrid*1350+838] = 4.0E0*I_ESP_I3y3z_G4y_aa-2.0E0*1*I_ESP_Gy3z_G4y_a-2.0E0*2*I_ESP_Gy3z_G4y_a;
    abcd[iGrid*1350+839] = 4.0E0*I_ESP_I2y4z_G4y_aa-2.0E0*1*I_ESP_G4z_G4y_a;
    abcd[iGrid*1350+840] = 4.0E0*I_ESP_I4x2y_G3yz_aa-2.0E0*1*I_ESP_G4x_G3yz_a;
    abcd[iGrid*1350+841] = 4.0E0*I_ESP_I3x3y_G3yz_aa-2.0E0*1*I_ESP_G3xy_G3yz_a-2.0E0*2*I_ESP_G3xy_G3yz_a;
    abcd[iGrid*1350+842] = 4.0E0*I_ESP_I3x2yz_G3yz_aa-2.0E0*1*I_ESP_G3xz_G3yz_a;
    abcd[iGrid*1350+843] = 4.0E0*I_ESP_I2x4y_G3yz_aa-2.0E0*2*I_ESP_G2x2y_G3yz_a-2.0E0*3*I_ESP_G2x2y_G3yz_a+2*1*I_ESP_D2x_G3yz;
    abcd[iGrid*1350+844] = 4.0E0*I_ESP_I2x3yz_G3yz_aa-2.0E0*1*I_ESP_G2xyz_G3yz_a-2.0E0*2*I_ESP_G2xyz_G3yz_a;
    abcd[iGrid*1350+845] = 4.0E0*I_ESP_I2x2y2z_G3yz_aa-2.0E0*1*I_ESP_G2x2z_G3yz_a;
    abcd[iGrid*1350+846] = 4.0E0*I_ESP_Ix5y_G3yz_aa-2.0E0*3*I_ESP_Gx3y_G3yz_a-2.0E0*4*I_ESP_Gx3y_G3yz_a+3*2*I_ESP_Dxy_G3yz;
    abcd[iGrid*1350+847] = 4.0E0*I_ESP_Ix4yz_G3yz_aa-2.0E0*2*I_ESP_Gx2yz_G3yz_a-2.0E0*3*I_ESP_Gx2yz_G3yz_a+2*1*I_ESP_Dxz_G3yz;
    abcd[iGrid*1350+848] = 4.0E0*I_ESP_Ix3y2z_G3yz_aa-2.0E0*1*I_ESP_Gxy2z_G3yz_a-2.0E0*2*I_ESP_Gxy2z_G3yz_a;
    abcd[iGrid*1350+849] = 4.0E0*I_ESP_Ix2y3z_G3yz_aa-2.0E0*1*I_ESP_Gx3z_G3yz_a;
    abcd[iGrid*1350+850] = 4.0E0*I_ESP_I6y_G3yz_aa-2.0E0*4*I_ESP_G4y_G3yz_a-2.0E0*5*I_ESP_G4y_G3yz_a+4*3*I_ESP_D2y_G3yz;
    abcd[iGrid*1350+851] = 4.0E0*I_ESP_I5yz_G3yz_aa-2.0E0*3*I_ESP_G3yz_G3yz_a-2.0E0*4*I_ESP_G3yz_G3yz_a+3*2*I_ESP_Dyz_G3yz;
    abcd[iGrid*1350+852] = 4.0E0*I_ESP_I4y2z_G3yz_aa-2.0E0*2*I_ESP_G2y2z_G3yz_a-2.0E0*3*I_ESP_G2y2z_G3yz_a+2*1*I_ESP_D2z_G3yz;
    abcd[iGrid*1350+853] = 4.0E0*I_ESP_I3y3z_G3yz_aa-2.0E0*1*I_ESP_Gy3z_G3yz_a-2.0E0*2*I_ESP_Gy3z_G3yz_a;
    abcd[iGrid*1350+854] = 4.0E0*I_ESP_I2y4z_G3yz_aa-2.0E0*1*I_ESP_G4z_G3yz_a;
    abcd[iGrid*1350+855] = 4.0E0*I_ESP_I4x2y_G2y2z_aa-2.0E0*1*I_ESP_G4x_G2y2z_a;
    abcd[iGrid*1350+856] = 4.0E0*I_ESP_I3x3y_G2y2z_aa-2.0E0*1*I_ESP_G3xy_G2y2z_a-2.0E0*2*I_ESP_G3xy_G2y2z_a;
    abcd[iGrid*1350+857] = 4.0E0*I_ESP_I3x2yz_G2y2z_aa-2.0E0*1*I_ESP_G3xz_G2y2z_a;
    abcd[iGrid*1350+858] = 4.0E0*I_ESP_I2x4y_G2y2z_aa-2.0E0*2*I_ESP_G2x2y_G2y2z_a-2.0E0*3*I_ESP_G2x2y_G2y2z_a+2*1*I_ESP_D2x_G2y2z;
    abcd[iGrid*1350+859] = 4.0E0*I_ESP_I2x3yz_G2y2z_aa-2.0E0*1*I_ESP_G2xyz_G2y2z_a-2.0E0*2*I_ESP_G2xyz_G2y2z_a;
    abcd[iGrid*1350+860] = 4.0E0*I_ESP_I2x2y2z_G2y2z_aa-2.0E0*1*I_ESP_G2x2z_G2y2z_a;
    abcd[iGrid*1350+861] = 4.0E0*I_ESP_Ix5y_G2y2z_aa-2.0E0*3*I_ESP_Gx3y_G2y2z_a-2.0E0*4*I_ESP_Gx3y_G2y2z_a+3*2*I_ESP_Dxy_G2y2z;
    abcd[iGrid*1350+862] = 4.0E0*I_ESP_Ix4yz_G2y2z_aa-2.0E0*2*I_ESP_Gx2yz_G2y2z_a-2.0E0*3*I_ESP_Gx2yz_G2y2z_a+2*1*I_ESP_Dxz_G2y2z;
    abcd[iGrid*1350+863] = 4.0E0*I_ESP_Ix3y2z_G2y2z_aa-2.0E0*1*I_ESP_Gxy2z_G2y2z_a-2.0E0*2*I_ESP_Gxy2z_G2y2z_a;
    abcd[iGrid*1350+864] = 4.0E0*I_ESP_Ix2y3z_G2y2z_aa-2.0E0*1*I_ESP_Gx3z_G2y2z_a;
    abcd[iGrid*1350+865] = 4.0E0*I_ESP_I6y_G2y2z_aa-2.0E0*4*I_ESP_G4y_G2y2z_a-2.0E0*5*I_ESP_G4y_G2y2z_a+4*3*I_ESP_D2y_G2y2z;
    abcd[iGrid*1350+866] = 4.0E0*I_ESP_I5yz_G2y2z_aa-2.0E0*3*I_ESP_G3yz_G2y2z_a-2.0E0*4*I_ESP_G3yz_G2y2z_a+3*2*I_ESP_Dyz_G2y2z;
    abcd[iGrid*1350+867] = 4.0E0*I_ESP_I4y2z_G2y2z_aa-2.0E0*2*I_ESP_G2y2z_G2y2z_a-2.0E0*3*I_ESP_G2y2z_G2y2z_a+2*1*I_ESP_D2z_G2y2z;
    abcd[iGrid*1350+868] = 4.0E0*I_ESP_I3y3z_G2y2z_aa-2.0E0*1*I_ESP_Gy3z_G2y2z_a-2.0E0*2*I_ESP_Gy3z_G2y2z_a;
    abcd[iGrid*1350+869] = 4.0E0*I_ESP_I2y4z_G2y2z_aa-2.0E0*1*I_ESP_G4z_G2y2z_a;
    abcd[iGrid*1350+870] = 4.0E0*I_ESP_I4x2y_Gy3z_aa-2.0E0*1*I_ESP_G4x_Gy3z_a;
    abcd[iGrid*1350+871] = 4.0E0*I_ESP_I3x3y_Gy3z_aa-2.0E0*1*I_ESP_G3xy_Gy3z_a-2.0E0*2*I_ESP_G3xy_Gy3z_a;
    abcd[iGrid*1350+872] = 4.0E0*I_ESP_I3x2yz_Gy3z_aa-2.0E0*1*I_ESP_G3xz_Gy3z_a;
    abcd[iGrid*1350+873] = 4.0E0*I_ESP_I2x4y_Gy3z_aa-2.0E0*2*I_ESP_G2x2y_Gy3z_a-2.0E0*3*I_ESP_G2x2y_Gy3z_a+2*1*I_ESP_D2x_Gy3z;
    abcd[iGrid*1350+874] = 4.0E0*I_ESP_I2x3yz_Gy3z_aa-2.0E0*1*I_ESP_G2xyz_Gy3z_a-2.0E0*2*I_ESP_G2xyz_Gy3z_a;
    abcd[iGrid*1350+875] = 4.0E0*I_ESP_I2x2y2z_Gy3z_aa-2.0E0*1*I_ESP_G2x2z_Gy3z_a;
    abcd[iGrid*1350+876] = 4.0E0*I_ESP_Ix5y_Gy3z_aa-2.0E0*3*I_ESP_Gx3y_Gy3z_a-2.0E0*4*I_ESP_Gx3y_Gy3z_a+3*2*I_ESP_Dxy_Gy3z;
    abcd[iGrid*1350+877] = 4.0E0*I_ESP_Ix4yz_Gy3z_aa-2.0E0*2*I_ESP_Gx2yz_Gy3z_a-2.0E0*3*I_ESP_Gx2yz_Gy3z_a+2*1*I_ESP_Dxz_Gy3z;
    abcd[iGrid*1350+878] = 4.0E0*I_ESP_Ix3y2z_Gy3z_aa-2.0E0*1*I_ESP_Gxy2z_Gy3z_a-2.0E0*2*I_ESP_Gxy2z_Gy3z_a;
    abcd[iGrid*1350+879] = 4.0E0*I_ESP_Ix2y3z_Gy3z_aa-2.0E0*1*I_ESP_Gx3z_Gy3z_a;
    abcd[iGrid*1350+880] = 4.0E0*I_ESP_I6y_Gy3z_aa-2.0E0*4*I_ESP_G4y_Gy3z_a-2.0E0*5*I_ESP_G4y_Gy3z_a+4*3*I_ESP_D2y_Gy3z;
    abcd[iGrid*1350+881] = 4.0E0*I_ESP_I5yz_Gy3z_aa-2.0E0*3*I_ESP_G3yz_Gy3z_a-2.0E0*4*I_ESP_G3yz_Gy3z_a+3*2*I_ESP_Dyz_Gy3z;
    abcd[iGrid*1350+882] = 4.0E0*I_ESP_I4y2z_Gy3z_aa-2.0E0*2*I_ESP_G2y2z_Gy3z_a-2.0E0*3*I_ESP_G2y2z_Gy3z_a+2*1*I_ESP_D2z_Gy3z;
    abcd[iGrid*1350+883] = 4.0E0*I_ESP_I3y3z_Gy3z_aa-2.0E0*1*I_ESP_Gy3z_Gy3z_a-2.0E0*2*I_ESP_Gy3z_Gy3z_a;
    abcd[iGrid*1350+884] = 4.0E0*I_ESP_I2y4z_Gy3z_aa-2.0E0*1*I_ESP_G4z_Gy3z_a;
    abcd[iGrid*1350+885] = 4.0E0*I_ESP_I4x2y_G4z_aa-2.0E0*1*I_ESP_G4x_G4z_a;
    abcd[iGrid*1350+886] = 4.0E0*I_ESP_I3x3y_G4z_aa-2.0E0*1*I_ESP_G3xy_G4z_a-2.0E0*2*I_ESP_G3xy_G4z_a;
    abcd[iGrid*1350+887] = 4.0E0*I_ESP_I3x2yz_G4z_aa-2.0E0*1*I_ESP_G3xz_G4z_a;
    abcd[iGrid*1350+888] = 4.0E0*I_ESP_I2x4y_G4z_aa-2.0E0*2*I_ESP_G2x2y_G4z_a-2.0E0*3*I_ESP_G2x2y_G4z_a+2*1*I_ESP_D2x_G4z;
    abcd[iGrid*1350+889] = 4.0E0*I_ESP_I2x3yz_G4z_aa-2.0E0*1*I_ESP_G2xyz_G4z_a-2.0E0*2*I_ESP_G2xyz_G4z_a;
    abcd[iGrid*1350+890] = 4.0E0*I_ESP_I2x2y2z_G4z_aa-2.0E0*1*I_ESP_G2x2z_G4z_a;
    abcd[iGrid*1350+891] = 4.0E0*I_ESP_Ix5y_G4z_aa-2.0E0*3*I_ESP_Gx3y_G4z_a-2.0E0*4*I_ESP_Gx3y_G4z_a+3*2*I_ESP_Dxy_G4z;
    abcd[iGrid*1350+892] = 4.0E0*I_ESP_Ix4yz_G4z_aa-2.0E0*2*I_ESP_Gx2yz_G4z_a-2.0E0*3*I_ESP_Gx2yz_G4z_a+2*1*I_ESP_Dxz_G4z;
    abcd[iGrid*1350+893] = 4.0E0*I_ESP_Ix3y2z_G4z_aa-2.0E0*1*I_ESP_Gxy2z_G4z_a-2.0E0*2*I_ESP_Gxy2z_G4z_a;
    abcd[iGrid*1350+894] = 4.0E0*I_ESP_Ix2y3z_G4z_aa-2.0E0*1*I_ESP_Gx3z_G4z_a;
    abcd[iGrid*1350+895] = 4.0E0*I_ESP_I6y_G4z_aa-2.0E0*4*I_ESP_G4y_G4z_a-2.0E0*5*I_ESP_G4y_G4z_a+4*3*I_ESP_D2y_G4z;
    abcd[iGrid*1350+896] = 4.0E0*I_ESP_I5yz_G4z_aa-2.0E0*3*I_ESP_G3yz_G4z_a-2.0E0*4*I_ESP_G3yz_G4z_a+3*2*I_ESP_Dyz_G4z;
    abcd[iGrid*1350+897] = 4.0E0*I_ESP_I4y2z_G4z_aa-2.0E0*2*I_ESP_G2y2z_G4z_a-2.0E0*3*I_ESP_G2y2z_G4z_a+2*1*I_ESP_D2z_G4z;
    abcd[iGrid*1350+898] = 4.0E0*I_ESP_I3y3z_G4z_aa-2.0E0*1*I_ESP_Gy3z_G4z_a-2.0E0*2*I_ESP_Gy3z_G4z_a;
    abcd[iGrid*1350+899] = 4.0E0*I_ESP_I2y4z_G4z_aa-2.0E0*1*I_ESP_G4z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G_aa
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_D_G
     ************************************************************/
    abcd[iGrid*1350+900] = 4.0E0*I_ESP_I4xyz_G4x_aa;
    abcd[iGrid*1350+901] = 4.0E0*I_ESP_I3x2yz_G4x_aa-2.0E0*1*I_ESP_G3xz_G4x_a;
    abcd[iGrid*1350+902] = 4.0E0*I_ESP_I3xy2z_G4x_aa-2.0E0*1*I_ESP_G3xy_G4x_a;
    abcd[iGrid*1350+903] = 4.0E0*I_ESP_I2x3yz_G4x_aa-2.0E0*2*I_ESP_G2xyz_G4x_a;
    abcd[iGrid*1350+904] = 4.0E0*I_ESP_I2x2y2z_G4x_aa-2.0E0*1*I_ESP_G2x2y_G4x_a-2.0E0*1*I_ESP_G2x2z_G4x_a+1*I_ESP_D2x_G4x;
    abcd[iGrid*1350+905] = 4.0E0*I_ESP_I2xy3z_G4x_aa-2.0E0*2*I_ESP_G2xyz_G4x_a;
    abcd[iGrid*1350+906] = 4.0E0*I_ESP_Ix4yz_G4x_aa-2.0E0*3*I_ESP_Gx2yz_G4x_a;
    abcd[iGrid*1350+907] = 4.0E0*I_ESP_Ix3y2z_G4x_aa-2.0E0*1*I_ESP_Gx3y_G4x_a-2.0E0*2*I_ESP_Gxy2z_G4x_a+2*1*I_ESP_Dxy_G4x;
    abcd[iGrid*1350+908] = 4.0E0*I_ESP_Ix2y3z_G4x_aa-2.0E0*2*I_ESP_Gx2yz_G4x_a-2.0E0*1*I_ESP_Gx3z_G4x_a+2*I_ESP_Dxz_G4x;
    abcd[iGrid*1350+909] = 4.0E0*I_ESP_Ixy4z_G4x_aa-2.0E0*3*I_ESP_Gxy2z_G4x_a;
    abcd[iGrid*1350+910] = 4.0E0*I_ESP_I5yz_G4x_aa-2.0E0*4*I_ESP_G3yz_G4x_a;
    abcd[iGrid*1350+911] = 4.0E0*I_ESP_I4y2z_G4x_aa-2.0E0*1*I_ESP_G4y_G4x_a-2.0E0*3*I_ESP_G2y2z_G4x_a+3*1*I_ESP_D2y_G4x;
    abcd[iGrid*1350+912] = 4.0E0*I_ESP_I3y3z_G4x_aa-2.0E0*2*I_ESP_G3yz_G4x_a-2.0E0*2*I_ESP_Gy3z_G4x_a+2*2*I_ESP_Dyz_G4x;
    abcd[iGrid*1350+913] = 4.0E0*I_ESP_I2y4z_G4x_aa-2.0E0*3*I_ESP_G2y2z_G4x_a-2.0E0*1*I_ESP_G4z_G4x_a+3*I_ESP_D2z_G4x;
    abcd[iGrid*1350+914] = 4.0E0*I_ESP_Iy5z_G4x_aa-2.0E0*4*I_ESP_Gy3z_G4x_a;
    abcd[iGrid*1350+915] = 4.0E0*I_ESP_I4xyz_G3xy_aa;
    abcd[iGrid*1350+916] = 4.0E0*I_ESP_I3x2yz_G3xy_aa-2.0E0*1*I_ESP_G3xz_G3xy_a;
    abcd[iGrid*1350+917] = 4.0E0*I_ESP_I3xy2z_G3xy_aa-2.0E0*1*I_ESP_G3xy_G3xy_a;
    abcd[iGrid*1350+918] = 4.0E0*I_ESP_I2x3yz_G3xy_aa-2.0E0*2*I_ESP_G2xyz_G3xy_a;
    abcd[iGrid*1350+919] = 4.0E0*I_ESP_I2x2y2z_G3xy_aa-2.0E0*1*I_ESP_G2x2y_G3xy_a-2.0E0*1*I_ESP_G2x2z_G3xy_a+1*I_ESP_D2x_G3xy;
    abcd[iGrid*1350+920] = 4.0E0*I_ESP_I2xy3z_G3xy_aa-2.0E0*2*I_ESP_G2xyz_G3xy_a;
    abcd[iGrid*1350+921] = 4.0E0*I_ESP_Ix4yz_G3xy_aa-2.0E0*3*I_ESP_Gx2yz_G3xy_a;
    abcd[iGrid*1350+922] = 4.0E0*I_ESP_Ix3y2z_G3xy_aa-2.0E0*1*I_ESP_Gx3y_G3xy_a-2.0E0*2*I_ESP_Gxy2z_G3xy_a+2*1*I_ESP_Dxy_G3xy;
    abcd[iGrid*1350+923] = 4.0E0*I_ESP_Ix2y3z_G3xy_aa-2.0E0*2*I_ESP_Gx2yz_G3xy_a-2.0E0*1*I_ESP_Gx3z_G3xy_a+2*I_ESP_Dxz_G3xy;
    abcd[iGrid*1350+924] = 4.0E0*I_ESP_Ixy4z_G3xy_aa-2.0E0*3*I_ESP_Gxy2z_G3xy_a;
    abcd[iGrid*1350+925] = 4.0E0*I_ESP_I5yz_G3xy_aa-2.0E0*4*I_ESP_G3yz_G3xy_a;
    abcd[iGrid*1350+926] = 4.0E0*I_ESP_I4y2z_G3xy_aa-2.0E0*1*I_ESP_G4y_G3xy_a-2.0E0*3*I_ESP_G2y2z_G3xy_a+3*1*I_ESP_D2y_G3xy;
    abcd[iGrid*1350+927] = 4.0E0*I_ESP_I3y3z_G3xy_aa-2.0E0*2*I_ESP_G3yz_G3xy_a-2.0E0*2*I_ESP_Gy3z_G3xy_a+2*2*I_ESP_Dyz_G3xy;
    abcd[iGrid*1350+928] = 4.0E0*I_ESP_I2y4z_G3xy_aa-2.0E0*3*I_ESP_G2y2z_G3xy_a-2.0E0*1*I_ESP_G4z_G3xy_a+3*I_ESP_D2z_G3xy;
    abcd[iGrid*1350+929] = 4.0E0*I_ESP_Iy5z_G3xy_aa-2.0E0*4*I_ESP_Gy3z_G3xy_a;
    abcd[iGrid*1350+930] = 4.0E0*I_ESP_I4xyz_G3xz_aa;
    abcd[iGrid*1350+931] = 4.0E0*I_ESP_I3x2yz_G3xz_aa-2.0E0*1*I_ESP_G3xz_G3xz_a;
    abcd[iGrid*1350+932] = 4.0E0*I_ESP_I3xy2z_G3xz_aa-2.0E0*1*I_ESP_G3xy_G3xz_a;
    abcd[iGrid*1350+933] = 4.0E0*I_ESP_I2x3yz_G3xz_aa-2.0E0*2*I_ESP_G2xyz_G3xz_a;
    abcd[iGrid*1350+934] = 4.0E0*I_ESP_I2x2y2z_G3xz_aa-2.0E0*1*I_ESP_G2x2y_G3xz_a-2.0E0*1*I_ESP_G2x2z_G3xz_a+1*I_ESP_D2x_G3xz;
    abcd[iGrid*1350+935] = 4.0E0*I_ESP_I2xy3z_G3xz_aa-2.0E0*2*I_ESP_G2xyz_G3xz_a;
    abcd[iGrid*1350+936] = 4.0E0*I_ESP_Ix4yz_G3xz_aa-2.0E0*3*I_ESP_Gx2yz_G3xz_a;
    abcd[iGrid*1350+937] = 4.0E0*I_ESP_Ix3y2z_G3xz_aa-2.0E0*1*I_ESP_Gx3y_G3xz_a-2.0E0*2*I_ESP_Gxy2z_G3xz_a+2*1*I_ESP_Dxy_G3xz;
    abcd[iGrid*1350+938] = 4.0E0*I_ESP_Ix2y3z_G3xz_aa-2.0E0*2*I_ESP_Gx2yz_G3xz_a-2.0E0*1*I_ESP_Gx3z_G3xz_a+2*I_ESP_Dxz_G3xz;
    abcd[iGrid*1350+939] = 4.0E0*I_ESP_Ixy4z_G3xz_aa-2.0E0*3*I_ESP_Gxy2z_G3xz_a;
    abcd[iGrid*1350+940] = 4.0E0*I_ESP_I5yz_G3xz_aa-2.0E0*4*I_ESP_G3yz_G3xz_a;
    abcd[iGrid*1350+941] = 4.0E0*I_ESP_I4y2z_G3xz_aa-2.0E0*1*I_ESP_G4y_G3xz_a-2.0E0*3*I_ESP_G2y2z_G3xz_a+3*1*I_ESP_D2y_G3xz;
    abcd[iGrid*1350+942] = 4.0E0*I_ESP_I3y3z_G3xz_aa-2.0E0*2*I_ESP_G3yz_G3xz_a-2.0E0*2*I_ESP_Gy3z_G3xz_a+2*2*I_ESP_Dyz_G3xz;
    abcd[iGrid*1350+943] = 4.0E0*I_ESP_I2y4z_G3xz_aa-2.0E0*3*I_ESP_G2y2z_G3xz_a-2.0E0*1*I_ESP_G4z_G3xz_a+3*I_ESP_D2z_G3xz;
    abcd[iGrid*1350+944] = 4.0E0*I_ESP_Iy5z_G3xz_aa-2.0E0*4*I_ESP_Gy3z_G3xz_a;
    abcd[iGrid*1350+945] = 4.0E0*I_ESP_I4xyz_G2x2y_aa;
    abcd[iGrid*1350+946] = 4.0E0*I_ESP_I3x2yz_G2x2y_aa-2.0E0*1*I_ESP_G3xz_G2x2y_a;
    abcd[iGrid*1350+947] = 4.0E0*I_ESP_I3xy2z_G2x2y_aa-2.0E0*1*I_ESP_G3xy_G2x2y_a;
    abcd[iGrid*1350+948] = 4.0E0*I_ESP_I2x3yz_G2x2y_aa-2.0E0*2*I_ESP_G2xyz_G2x2y_a;
    abcd[iGrid*1350+949] = 4.0E0*I_ESP_I2x2y2z_G2x2y_aa-2.0E0*1*I_ESP_G2x2y_G2x2y_a-2.0E0*1*I_ESP_G2x2z_G2x2y_a+1*I_ESP_D2x_G2x2y;
    abcd[iGrid*1350+950] = 4.0E0*I_ESP_I2xy3z_G2x2y_aa-2.0E0*2*I_ESP_G2xyz_G2x2y_a;
    abcd[iGrid*1350+951] = 4.0E0*I_ESP_Ix4yz_G2x2y_aa-2.0E0*3*I_ESP_Gx2yz_G2x2y_a;
    abcd[iGrid*1350+952] = 4.0E0*I_ESP_Ix3y2z_G2x2y_aa-2.0E0*1*I_ESP_Gx3y_G2x2y_a-2.0E0*2*I_ESP_Gxy2z_G2x2y_a+2*1*I_ESP_Dxy_G2x2y;
    abcd[iGrid*1350+953] = 4.0E0*I_ESP_Ix2y3z_G2x2y_aa-2.0E0*2*I_ESP_Gx2yz_G2x2y_a-2.0E0*1*I_ESP_Gx3z_G2x2y_a+2*I_ESP_Dxz_G2x2y;
    abcd[iGrid*1350+954] = 4.0E0*I_ESP_Ixy4z_G2x2y_aa-2.0E0*3*I_ESP_Gxy2z_G2x2y_a;
    abcd[iGrid*1350+955] = 4.0E0*I_ESP_I5yz_G2x2y_aa-2.0E0*4*I_ESP_G3yz_G2x2y_a;
    abcd[iGrid*1350+956] = 4.0E0*I_ESP_I4y2z_G2x2y_aa-2.0E0*1*I_ESP_G4y_G2x2y_a-2.0E0*3*I_ESP_G2y2z_G2x2y_a+3*1*I_ESP_D2y_G2x2y;
    abcd[iGrid*1350+957] = 4.0E0*I_ESP_I3y3z_G2x2y_aa-2.0E0*2*I_ESP_G3yz_G2x2y_a-2.0E0*2*I_ESP_Gy3z_G2x2y_a+2*2*I_ESP_Dyz_G2x2y;
    abcd[iGrid*1350+958] = 4.0E0*I_ESP_I2y4z_G2x2y_aa-2.0E0*3*I_ESP_G2y2z_G2x2y_a-2.0E0*1*I_ESP_G4z_G2x2y_a+3*I_ESP_D2z_G2x2y;
    abcd[iGrid*1350+959] = 4.0E0*I_ESP_Iy5z_G2x2y_aa-2.0E0*4*I_ESP_Gy3z_G2x2y_a;
    abcd[iGrid*1350+960] = 4.0E0*I_ESP_I4xyz_G2xyz_aa;
    abcd[iGrid*1350+961] = 4.0E0*I_ESP_I3x2yz_G2xyz_aa-2.0E0*1*I_ESP_G3xz_G2xyz_a;
    abcd[iGrid*1350+962] = 4.0E0*I_ESP_I3xy2z_G2xyz_aa-2.0E0*1*I_ESP_G3xy_G2xyz_a;
    abcd[iGrid*1350+963] = 4.0E0*I_ESP_I2x3yz_G2xyz_aa-2.0E0*2*I_ESP_G2xyz_G2xyz_a;
    abcd[iGrid*1350+964] = 4.0E0*I_ESP_I2x2y2z_G2xyz_aa-2.0E0*1*I_ESP_G2x2y_G2xyz_a-2.0E0*1*I_ESP_G2x2z_G2xyz_a+1*I_ESP_D2x_G2xyz;
    abcd[iGrid*1350+965] = 4.0E0*I_ESP_I2xy3z_G2xyz_aa-2.0E0*2*I_ESP_G2xyz_G2xyz_a;
    abcd[iGrid*1350+966] = 4.0E0*I_ESP_Ix4yz_G2xyz_aa-2.0E0*3*I_ESP_Gx2yz_G2xyz_a;
    abcd[iGrid*1350+967] = 4.0E0*I_ESP_Ix3y2z_G2xyz_aa-2.0E0*1*I_ESP_Gx3y_G2xyz_a-2.0E0*2*I_ESP_Gxy2z_G2xyz_a+2*1*I_ESP_Dxy_G2xyz;
    abcd[iGrid*1350+968] = 4.0E0*I_ESP_Ix2y3z_G2xyz_aa-2.0E0*2*I_ESP_Gx2yz_G2xyz_a-2.0E0*1*I_ESP_Gx3z_G2xyz_a+2*I_ESP_Dxz_G2xyz;
    abcd[iGrid*1350+969] = 4.0E0*I_ESP_Ixy4z_G2xyz_aa-2.0E0*3*I_ESP_Gxy2z_G2xyz_a;
    abcd[iGrid*1350+970] = 4.0E0*I_ESP_I5yz_G2xyz_aa-2.0E0*4*I_ESP_G3yz_G2xyz_a;
    abcd[iGrid*1350+971] = 4.0E0*I_ESP_I4y2z_G2xyz_aa-2.0E0*1*I_ESP_G4y_G2xyz_a-2.0E0*3*I_ESP_G2y2z_G2xyz_a+3*1*I_ESP_D2y_G2xyz;
    abcd[iGrid*1350+972] = 4.0E0*I_ESP_I3y3z_G2xyz_aa-2.0E0*2*I_ESP_G3yz_G2xyz_a-2.0E0*2*I_ESP_Gy3z_G2xyz_a+2*2*I_ESP_Dyz_G2xyz;
    abcd[iGrid*1350+973] = 4.0E0*I_ESP_I2y4z_G2xyz_aa-2.0E0*3*I_ESP_G2y2z_G2xyz_a-2.0E0*1*I_ESP_G4z_G2xyz_a+3*I_ESP_D2z_G2xyz;
    abcd[iGrid*1350+974] = 4.0E0*I_ESP_Iy5z_G2xyz_aa-2.0E0*4*I_ESP_Gy3z_G2xyz_a;
    abcd[iGrid*1350+975] = 4.0E0*I_ESP_I4xyz_G2x2z_aa;
    abcd[iGrid*1350+976] = 4.0E0*I_ESP_I3x2yz_G2x2z_aa-2.0E0*1*I_ESP_G3xz_G2x2z_a;
    abcd[iGrid*1350+977] = 4.0E0*I_ESP_I3xy2z_G2x2z_aa-2.0E0*1*I_ESP_G3xy_G2x2z_a;
    abcd[iGrid*1350+978] = 4.0E0*I_ESP_I2x3yz_G2x2z_aa-2.0E0*2*I_ESP_G2xyz_G2x2z_a;
    abcd[iGrid*1350+979] = 4.0E0*I_ESP_I2x2y2z_G2x2z_aa-2.0E0*1*I_ESP_G2x2y_G2x2z_a-2.0E0*1*I_ESP_G2x2z_G2x2z_a+1*I_ESP_D2x_G2x2z;
    abcd[iGrid*1350+980] = 4.0E0*I_ESP_I2xy3z_G2x2z_aa-2.0E0*2*I_ESP_G2xyz_G2x2z_a;
    abcd[iGrid*1350+981] = 4.0E0*I_ESP_Ix4yz_G2x2z_aa-2.0E0*3*I_ESP_Gx2yz_G2x2z_a;
    abcd[iGrid*1350+982] = 4.0E0*I_ESP_Ix3y2z_G2x2z_aa-2.0E0*1*I_ESP_Gx3y_G2x2z_a-2.0E0*2*I_ESP_Gxy2z_G2x2z_a+2*1*I_ESP_Dxy_G2x2z;
    abcd[iGrid*1350+983] = 4.0E0*I_ESP_Ix2y3z_G2x2z_aa-2.0E0*2*I_ESP_Gx2yz_G2x2z_a-2.0E0*1*I_ESP_Gx3z_G2x2z_a+2*I_ESP_Dxz_G2x2z;
    abcd[iGrid*1350+984] = 4.0E0*I_ESP_Ixy4z_G2x2z_aa-2.0E0*3*I_ESP_Gxy2z_G2x2z_a;
    abcd[iGrid*1350+985] = 4.0E0*I_ESP_I5yz_G2x2z_aa-2.0E0*4*I_ESP_G3yz_G2x2z_a;
    abcd[iGrid*1350+986] = 4.0E0*I_ESP_I4y2z_G2x2z_aa-2.0E0*1*I_ESP_G4y_G2x2z_a-2.0E0*3*I_ESP_G2y2z_G2x2z_a+3*1*I_ESP_D2y_G2x2z;
    abcd[iGrid*1350+987] = 4.0E0*I_ESP_I3y3z_G2x2z_aa-2.0E0*2*I_ESP_G3yz_G2x2z_a-2.0E0*2*I_ESP_Gy3z_G2x2z_a+2*2*I_ESP_Dyz_G2x2z;
    abcd[iGrid*1350+988] = 4.0E0*I_ESP_I2y4z_G2x2z_aa-2.0E0*3*I_ESP_G2y2z_G2x2z_a-2.0E0*1*I_ESP_G4z_G2x2z_a+3*I_ESP_D2z_G2x2z;
    abcd[iGrid*1350+989] = 4.0E0*I_ESP_Iy5z_G2x2z_aa-2.0E0*4*I_ESP_Gy3z_G2x2z_a;
    abcd[iGrid*1350+990] = 4.0E0*I_ESP_I4xyz_Gx3y_aa;
    abcd[iGrid*1350+991] = 4.0E0*I_ESP_I3x2yz_Gx3y_aa-2.0E0*1*I_ESP_G3xz_Gx3y_a;
    abcd[iGrid*1350+992] = 4.0E0*I_ESP_I3xy2z_Gx3y_aa-2.0E0*1*I_ESP_G3xy_Gx3y_a;
    abcd[iGrid*1350+993] = 4.0E0*I_ESP_I2x3yz_Gx3y_aa-2.0E0*2*I_ESP_G2xyz_Gx3y_a;
    abcd[iGrid*1350+994] = 4.0E0*I_ESP_I2x2y2z_Gx3y_aa-2.0E0*1*I_ESP_G2x2y_Gx3y_a-2.0E0*1*I_ESP_G2x2z_Gx3y_a+1*I_ESP_D2x_Gx3y;
    abcd[iGrid*1350+995] = 4.0E0*I_ESP_I2xy3z_Gx3y_aa-2.0E0*2*I_ESP_G2xyz_Gx3y_a;
    abcd[iGrid*1350+996] = 4.0E0*I_ESP_Ix4yz_Gx3y_aa-2.0E0*3*I_ESP_Gx2yz_Gx3y_a;
    abcd[iGrid*1350+997] = 4.0E0*I_ESP_Ix3y2z_Gx3y_aa-2.0E0*1*I_ESP_Gx3y_Gx3y_a-2.0E0*2*I_ESP_Gxy2z_Gx3y_a+2*1*I_ESP_Dxy_Gx3y;
    abcd[iGrid*1350+998] = 4.0E0*I_ESP_Ix2y3z_Gx3y_aa-2.0E0*2*I_ESP_Gx2yz_Gx3y_a-2.0E0*1*I_ESP_Gx3z_Gx3y_a+2*I_ESP_Dxz_Gx3y;
    abcd[iGrid*1350+999] = 4.0E0*I_ESP_Ixy4z_Gx3y_aa-2.0E0*3*I_ESP_Gxy2z_Gx3y_a;
    abcd[iGrid*1350+1000] = 4.0E0*I_ESP_I5yz_Gx3y_aa-2.0E0*4*I_ESP_G3yz_Gx3y_a;
    abcd[iGrid*1350+1001] = 4.0E0*I_ESP_I4y2z_Gx3y_aa-2.0E0*1*I_ESP_G4y_Gx3y_a-2.0E0*3*I_ESP_G2y2z_Gx3y_a+3*1*I_ESP_D2y_Gx3y;
    abcd[iGrid*1350+1002] = 4.0E0*I_ESP_I3y3z_Gx3y_aa-2.0E0*2*I_ESP_G3yz_Gx3y_a-2.0E0*2*I_ESP_Gy3z_Gx3y_a+2*2*I_ESP_Dyz_Gx3y;
    abcd[iGrid*1350+1003] = 4.0E0*I_ESP_I2y4z_Gx3y_aa-2.0E0*3*I_ESP_G2y2z_Gx3y_a-2.0E0*1*I_ESP_G4z_Gx3y_a+3*I_ESP_D2z_Gx3y;
    abcd[iGrid*1350+1004] = 4.0E0*I_ESP_Iy5z_Gx3y_aa-2.0E0*4*I_ESP_Gy3z_Gx3y_a;
    abcd[iGrid*1350+1005] = 4.0E0*I_ESP_I4xyz_Gx2yz_aa;
    abcd[iGrid*1350+1006] = 4.0E0*I_ESP_I3x2yz_Gx2yz_aa-2.0E0*1*I_ESP_G3xz_Gx2yz_a;
    abcd[iGrid*1350+1007] = 4.0E0*I_ESP_I3xy2z_Gx2yz_aa-2.0E0*1*I_ESP_G3xy_Gx2yz_a;
    abcd[iGrid*1350+1008] = 4.0E0*I_ESP_I2x3yz_Gx2yz_aa-2.0E0*2*I_ESP_G2xyz_Gx2yz_a;
    abcd[iGrid*1350+1009] = 4.0E0*I_ESP_I2x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_G2x2y_Gx2yz_a-2.0E0*1*I_ESP_G2x2z_Gx2yz_a+1*I_ESP_D2x_Gx2yz;
    abcd[iGrid*1350+1010] = 4.0E0*I_ESP_I2xy3z_Gx2yz_aa-2.0E0*2*I_ESP_G2xyz_Gx2yz_a;
    abcd[iGrid*1350+1011] = 4.0E0*I_ESP_Ix4yz_Gx2yz_aa-2.0E0*3*I_ESP_Gx2yz_Gx2yz_a;
    abcd[iGrid*1350+1012] = 4.0E0*I_ESP_Ix3y2z_Gx2yz_aa-2.0E0*1*I_ESP_Gx3y_Gx2yz_a-2.0E0*2*I_ESP_Gxy2z_Gx2yz_a+2*1*I_ESP_Dxy_Gx2yz;
    abcd[iGrid*1350+1013] = 4.0E0*I_ESP_Ix2y3z_Gx2yz_aa-2.0E0*2*I_ESP_Gx2yz_Gx2yz_a-2.0E0*1*I_ESP_Gx3z_Gx2yz_a+2*I_ESP_Dxz_Gx2yz;
    abcd[iGrid*1350+1014] = 4.0E0*I_ESP_Ixy4z_Gx2yz_aa-2.0E0*3*I_ESP_Gxy2z_Gx2yz_a;
    abcd[iGrid*1350+1015] = 4.0E0*I_ESP_I5yz_Gx2yz_aa-2.0E0*4*I_ESP_G3yz_Gx2yz_a;
    abcd[iGrid*1350+1016] = 4.0E0*I_ESP_I4y2z_Gx2yz_aa-2.0E0*1*I_ESP_G4y_Gx2yz_a-2.0E0*3*I_ESP_G2y2z_Gx2yz_a+3*1*I_ESP_D2y_Gx2yz;
    abcd[iGrid*1350+1017] = 4.0E0*I_ESP_I3y3z_Gx2yz_aa-2.0E0*2*I_ESP_G3yz_Gx2yz_a-2.0E0*2*I_ESP_Gy3z_Gx2yz_a+2*2*I_ESP_Dyz_Gx2yz;
    abcd[iGrid*1350+1018] = 4.0E0*I_ESP_I2y4z_Gx2yz_aa-2.0E0*3*I_ESP_G2y2z_Gx2yz_a-2.0E0*1*I_ESP_G4z_Gx2yz_a+3*I_ESP_D2z_Gx2yz;
    abcd[iGrid*1350+1019] = 4.0E0*I_ESP_Iy5z_Gx2yz_aa-2.0E0*4*I_ESP_Gy3z_Gx2yz_a;
    abcd[iGrid*1350+1020] = 4.0E0*I_ESP_I4xyz_Gxy2z_aa;
    abcd[iGrid*1350+1021] = 4.0E0*I_ESP_I3x2yz_Gxy2z_aa-2.0E0*1*I_ESP_G3xz_Gxy2z_a;
    abcd[iGrid*1350+1022] = 4.0E0*I_ESP_I3xy2z_Gxy2z_aa-2.0E0*1*I_ESP_G3xy_Gxy2z_a;
    abcd[iGrid*1350+1023] = 4.0E0*I_ESP_I2x3yz_Gxy2z_aa-2.0E0*2*I_ESP_G2xyz_Gxy2z_a;
    abcd[iGrid*1350+1024] = 4.0E0*I_ESP_I2x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_G2x2y_Gxy2z_a-2.0E0*1*I_ESP_G2x2z_Gxy2z_a+1*I_ESP_D2x_Gxy2z;
    abcd[iGrid*1350+1025] = 4.0E0*I_ESP_I2xy3z_Gxy2z_aa-2.0E0*2*I_ESP_G2xyz_Gxy2z_a;
    abcd[iGrid*1350+1026] = 4.0E0*I_ESP_Ix4yz_Gxy2z_aa-2.0E0*3*I_ESP_Gx2yz_Gxy2z_a;
    abcd[iGrid*1350+1027] = 4.0E0*I_ESP_Ix3y2z_Gxy2z_aa-2.0E0*1*I_ESP_Gx3y_Gxy2z_a-2.0E0*2*I_ESP_Gxy2z_Gxy2z_a+2*1*I_ESP_Dxy_Gxy2z;
    abcd[iGrid*1350+1028] = 4.0E0*I_ESP_Ix2y3z_Gxy2z_aa-2.0E0*2*I_ESP_Gx2yz_Gxy2z_a-2.0E0*1*I_ESP_Gx3z_Gxy2z_a+2*I_ESP_Dxz_Gxy2z;
    abcd[iGrid*1350+1029] = 4.0E0*I_ESP_Ixy4z_Gxy2z_aa-2.0E0*3*I_ESP_Gxy2z_Gxy2z_a;
    abcd[iGrid*1350+1030] = 4.0E0*I_ESP_I5yz_Gxy2z_aa-2.0E0*4*I_ESP_G3yz_Gxy2z_a;
    abcd[iGrid*1350+1031] = 4.0E0*I_ESP_I4y2z_Gxy2z_aa-2.0E0*1*I_ESP_G4y_Gxy2z_a-2.0E0*3*I_ESP_G2y2z_Gxy2z_a+3*1*I_ESP_D2y_Gxy2z;
    abcd[iGrid*1350+1032] = 4.0E0*I_ESP_I3y3z_Gxy2z_aa-2.0E0*2*I_ESP_G3yz_Gxy2z_a-2.0E0*2*I_ESP_Gy3z_Gxy2z_a+2*2*I_ESP_Dyz_Gxy2z;
    abcd[iGrid*1350+1033] = 4.0E0*I_ESP_I2y4z_Gxy2z_aa-2.0E0*3*I_ESP_G2y2z_Gxy2z_a-2.0E0*1*I_ESP_G4z_Gxy2z_a+3*I_ESP_D2z_Gxy2z;
    abcd[iGrid*1350+1034] = 4.0E0*I_ESP_Iy5z_Gxy2z_aa-2.0E0*4*I_ESP_Gy3z_Gxy2z_a;
    abcd[iGrid*1350+1035] = 4.0E0*I_ESP_I4xyz_Gx3z_aa;
    abcd[iGrid*1350+1036] = 4.0E0*I_ESP_I3x2yz_Gx3z_aa-2.0E0*1*I_ESP_G3xz_Gx3z_a;
    abcd[iGrid*1350+1037] = 4.0E0*I_ESP_I3xy2z_Gx3z_aa-2.0E0*1*I_ESP_G3xy_Gx3z_a;
    abcd[iGrid*1350+1038] = 4.0E0*I_ESP_I2x3yz_Gx3z_aa-2.0E0*2*I_ESP_G2xyz_Gx3z_a;
    abcd[iGrid*1350+1039] = 4.0E0*I_ESP_I2x2y2z_Gx3z_aa-2.0E0*1*I_ESP_G2x2y_Gx3z_a-2.0E0*1*I_ESP_G2x2z_Gx3z_a+1*I_ESP_D2x_Gx3z;
    abcd[iGrid*1350+1040] = 4.0E0*I_ESP_I2xy3z_Gx3z_aa-2.0E0*2*I_ESP_G2xyz_Gx3z_a;
    abcd[iGrid*1350+1041] = 4.0E0*I_ESP_Ix4yz_Gx3z_aa-2.0E0*3*I_ESP_Gx2yz_Gx3z_a;
    abcd[iGrid*1350+1042] = 4.0E0*I_ESP_Ix3y2z_Gx3z_aa-2.0E0*1*I_ESP_Gx3y_Gx3z_a-2.0E0*2*I_ESP_Gxy2z_Gx3z_a+2*1*I_ESP_Dxy_Gx3z;
    abcd[iGrid*1350+1043] = 4.0E0*I_ESP_Ix2y3z_Gx3z_aa-2.0E0*2*I_ESP_Gx2yz_Gx3z_a-2.0E0*1*I_ESP_Gx3z_Gx3z_a+2*I_ESP_Dxz_Gx3z;
    abcd[iGrid*1350+1044] = 4.0E0*I_ESP_Ixy4z_Gx3z_aa-2.0E0*3*I_ESP_Gxy2z_Gx3z_a;
    abcd[iGrid*1350+1045] = 4.0E0*I_ESP_I5yz_Gx3z_aa-2.0E0*4*I_ESP_G3yz_Gx3z_a;
    abcd[iGrid*1350+1046] = 4.0E0*I_ESP_I4y2z_Gx3z_aa-2.0E0*1*I_ESP_G4y_Gx3z_a-2.0E0*3*I_ESP_G2y2z_Gx3z_a+3*1*I_ESP_D2y_Gx3z;
    abcd[iGrid*1350+1047] = 4.0E0*I_ESP_I3y3z_Gx3z_aa-2.0E0*2*I_ESP_G3yz_Gx3z_a-2.0E0*2*I_ESP_Gy3z_Gx3z_a+2*2*I_ESP_Dyz_Gx3z;
    abcd[iGrid*1350+1048] = 4.0E0*I_ESP_I2y4z_Gx3z_aa-2.0E0*3*I_ESP_G2y2z_Gx3z_a-2.0E0*1*I_ESP_G4z_Gx3z_a+3*I_ESP_D2z_Gx3z;
    abcd[iGrid*1350+1049] = 4.0E0*I_ESP_Iy5z_Gx3z_aa-2.0E0*4*I_ESP_Gy3z_Gx3z_a;
    abcd[iGrid*1350+1050] = 4.0E0*I_ESP_I4xyz_G4y_aa;
    abcd[iGrid*1350+1051] = 4.0E0*I_ESP_I3x2yz_G4y_aa-2.0E0*1*I_ESP_G3xz_G4y_a;
    abcd[iGrid*1350+1052] = 4.0E0*I_ESP_I3xy2z_G4y_aa-2.0E0*1*I_ESP_G3xy_G4y_a;
    abcd[iGrid*1350+1053] = 4.0E0*I_ESP_I2x3yz_G4y_aa-2.0E0*2*I_ESP_G2xyz_G4y_a;
    abcd[iGrid*1350+1054] = 4.0E0*I_ESP_I2x2y2z_G4y_aa-2.0E0*1*I_ESP_G2x2y_G4y_a-2.0E0*1*I_ESP_G2x2z_G4y_a+1*I_ESP_D2x_G4y;
    abcd[iGrid*1350+1055] = 4.0E0*I_ESP_I2xy3z_G4y_aa-2.0E0*2*I_ESP_G2xyz_G4y_a;
    abcd[iGrid*1350+1056] = 4.0E0*I_ESP_Ix4yz_G4y_aa-2.0E0*3*I_ESP_Gx2yz_G4y_a;
    abcd[iGrid*1350+1057] = 4.0E0*I_ESP_Ix3y2z_G4y_aa-2.0E0*1*I_ESP_Gx3y_G4y_a-2.0E0*2*I_ESP_Gxy2z_G4y_a+2*1*I_ESP_Dxy_G4y;
    abcd[iGrid*1350+1058] = 4.0E0*I_ESP_Ix2y3z_G4y_aa-2.0E0*2*I_ESP_Gx2yz_G4y_a-2.0E0*1*I_ESP_Gx3z_G4y_a+2*I_ESP_Dxz_G4y;
    abcd[iGrid*1350+1059] = 4.0E0*I_ESP_Ixy4z_G4y_aa-2.0E0*3*I_ESP_Gxy2z_G4y_a;
    abcd[iGrid*1350+1060] = 4.0E0*I_ESP_I5yz_G4y_aa-2.0E0*4*I_ESP_G3yz_G4y_a;
    abcd[iGrid*1350+1061] = 4.0E0*I_ESP_I4y2z_G4y_aa-2.0E0*1*I_ESP_G4y_G4y_a-2.0E0*3*I_ESP_G2y2z_G4y_a+3*1*I_ESP_D2y_G4y;
    abcd[iGrid*1350+1062] = 4.0E0*I_ESP_I3y3z_G4y_aa-2.0E0*2*I_ESP_G3yz_G4y_a-2.0E0*2*I_ESP_Gy3z_G4y_a+2*2*I_ESP_Dyz_G4y;
    abcd[iGrid*1350+1063] = 4.0E0*I_ESP_I2y4z_G4y_aa-2.0E0*3*I_ESP_G2y2z_G4y_a-2.0E0*1*I_ESP_G4z_G4y_a+3*I_ESP_D2z_G4y;
    abcd[iGrid*1350+1064] = 4.0E0*I_ESP_Iy5z_G4y_aa-2.0E0*4*I_ESP_Gy3z_G4y_a;
    abcd[iGrid*1350+1065] = 4.0E0*I_ESP_I4xyz_G3yz_aa;
    abcd[iGrid*1350+1066] = 4.0E0*I_ESP_I3x2yz_G3yz_aa-2.0E0*1*I_ESP_G3xz_G3yz_a;
    abcd[iGrid*1350+1067] = 4.0E0*I_ESP_I3xy2z_G3yz_aa-2.0E0*1*I_ESP_G3xy_G3yz_a;
    abcd[iGrid*1350+1068] = 4.0E0*I_ESP_I2x3yz_G3yz_aa-2.0E0*2*I_ESP_G2xyz_G3yz_a;
    abcd[iGrid*1350+1069] = 4.0E0*I_ESP_I2x2y2z_G3yz_aa-2.0E0*1*I_ESP_G2x2y_G3yz_a-2.0E0*1*I_ESP_G2x2z_G3yz_a+1*I_ESP_D2x_G3yz;
    abcd[iGrid*1350+1070] = 4.0E0*I_ESP_I2xy3z_G3yz_aa-2.0E0*2*I_ESP_G2xyz_G3yz_a;
    abcd[iGrid*1350+1071] = 4.0E0*I_ESP_Ix4yz_G3yz_aa-2.0E0*3*I_ESP_Gx2yz_G3yz_a;
    abcd[iGrid*1350+1072] = 4.0E0*I_ESP_Ix3y2z_G3yz_aa-2.0E0*1*I_ESP_Gx3y_G3yz_a-2.0E0*2*I_ESP_Gxy2z_G3yz_a+2*1*I_ESP_Dxy_G3yz;
    abcd[iGrid*1350+1073] = 4.0E0*I_ESP_Ix2y3z_G3yz_aa-2.0E0*2*I_ESP_Gx2yz_G3yz_a-2.0E0*1*I_ESP_Gx3z_G3yz_a+2*I_ESP_Dxz_G3yz;
    abcd[iGrid*1350+1074] = 4.0E0*I_ESP_Ixy4z_G3yz_aa-2.0E0*3*I_ESP_Gxy2z_G3yz_a;
    abcd[iGrid*1350+1075] = 4.0E0*I_ESP_I5yz_G3yz_aa-2.0E0*4*I_ESP_G3yz_G3yz_a;
    abcd[iGrid*1350+1076] = 4.0E0*I_ESP_I4y2z_G3yz_aa-2.0E0*1*I_ESP_G4y_G3yz_a-2.0E0*3*I_ESP_G2y2z_G3yz_a+3*1*I_ESP_D2y_G3yz;
    abcd[iGrid*1350+1077] = 4.0E0*I_ESP_I3y3z_G3yz_aa-2.0E0*2*I_ESP_G3yz_G3yz_a-2.0E0*2*I_ESP_Gy3z_G3yz_a+2*2*I_ESP_Dyz_G3yz;
    abcd[iGrid*1350+1078] = 4.0E0*I_ESP_I2y4z_G3yz_aa-2.0E0*3*I_ESP_G2y2z_G3yz_a-2.0E0*1*I_ESP_G4z_G3yz_a+3*I_ESP_D2z_G3yz;
    abcd[iGrid*1350+1079] = 4.0E0*I_ESP_Iy5z_G3yz_aa-2.0E0*4*I_ESP_Gy3z_G3yz_a;
    abcd[iGrid*1350+1080] = 4.0E0*I_ESP_I4xyz_G2y2z_aa;
    abcd[iGrid*1350+1081] = 4.0E0*I_ESP_I3x2yz_G2y2z_aa-2.0E0*1*I_ESP_G3xz_G2y2z_a;
    abcd[iGrid*1350+1082] = 4.0E0*I_ESP_I3xy2z_G2y2z_aa-2.0E0*1*I_ESP_G3xy_G2y2z_a;
    abcd[iGrid*1350+1083] = 4.0E0*I_ESP_I2x3yz_G2y2z_aa-2.0E0*2*I_ESP_G2xyz_G2y2z_a;
    abcd[iGrid*1350+1084] = 4.0E0*I_ESP_I2x2y2z_G2y2z_aa-2.0E0*1*I_ESP_G2x2y_G2y2z_a-2.0E0*1*I_ESP_G2x2z_G2y2z_a+1*I_ESP_D2x_G2y2z;
    abcd[iGrid*1350+1085] = 4.0E0*I_ESP_I2xy3z_G2y2z_aa-2.0E0*2*I_ESP_G2xyz_G2y2z_a;
    abcd[iGrid*1350+1086] = 4.0E0*I_ESP_Ix4yz_G2y2z_aa-2.0E0*3*I_ESP_Gx2yz_G2y2z_a;
    abcd[iGrid*1350+1087] = 4.0E0*I_ESP_Ix3y2z_G2y2z_aa-2.0E0*1*I_ESP_Gx3y_G2y2z_a-2.0E0*2*I_ESP_Gxy2z_G2y2z_a+2*1*I_ESP_Dxy_G2y2z;
    abcd[iGrid*1350+1088] = 4.0E0*I_ESP_Ix2y3z_G2y2z_aa-2.0E0*2*I_ESP_Gx2yz_G2y2z_a-2.0E0*1*I_ESP_Gx3z_G2y2z_a+2*I_ESP_Dxz_G2y2z;
    abcd[iGrid*1350+1089] = 4.0E0*I_ESP_Ixy4z_G2y2z_aa-2.0E0*3*I_ESP_Gxy2z_G2y2z_a;
    abcd[iGrid*1350+1090] = 4.0E0*I_ESP_I5yz_G2y2z_aa-2.0E0*4*I_ESP_G3yz_G2y2z_a;
    abcd[iGrid*1350+1091] = 4.0E0*I_ESP_I4y2z_G2y2z_aa-2.0E0*1*I_ESP_G4y_G2y2z_a-2.0E0*3*I_ESP_G2y2z_G2y2z_a+3*1*I_ESP_D2y_G2y2z;
    abcd[iGrid*1350+1092] = 4.0E0*I_ESP_I3y3z_G2y2z_aa-2.0E0*2*I_ESP_G3yz_G2y2z_a-2.0E0*2*I_ESP_Gy3z_G2y2z_a+2*2*I_ESP_Dyz_G2y2z;
    abcd[iGrid*1350+1093] = 4.0E0*I_ESP_I2y4z_G2y2z_aa-2.0E0*3*I_ESP_G2y2z_G2y2z_a-2.0E0*1*I_ESP_G4z_G2y2z_a+3*I_ESP_D2z_G2y2z;
    abcd[iGrid*1350+1094] = 4.0E0*I_ESP_Iy5z_G2y2z_aa-2.0E0*4*I_ESP_Gy3z_G2y2z_a;
    abcd[iGrid*1350+1095] = 4.0E0*I_ESP_I4xyz_Gy3z_aa;
    abcd[iGrid*1350+1096] = 4.0E0*I_ESP_I3x2yz_Gy3z_aa-2.0E0*1*I_ESP_G3xz_Gy3z_a;
    abcd[iGrid*1350+1097] = 4.0E0*I_ESP_I3xy2z_Gy3z_aa-2.0E0*1*I_ESP_G3xy_Gy3z_a;
    abcd[iGrid*1350+1098] = 4.0E0*I_ESP_I2x3yz_Gy3z_aa-2.0E0*2*I_ESP_G2xyz_Gy3z_a;
    abcd[iGrid*1350+1099] = 4.0E0*I_ESP_I2x2y2z_Gy3z_aa-2.0E0*1*I_ESP_G2x2y_Gy3z_a-2.0E0*1*I_ESP_G2x2z_Gy3z_a+1*I_ESP_D2x_Gy3z;
    abcd[iGrid*1350+1100] = 4.0E0*I_ESP_I2xy3z_Gy3z_aa-2.0E0*2*I_ESP_G2xyz_Gy3z_a;
    abcd[iGrid*1350+1101] = 4.0E0*I_ESP_Ix4yz_Gy3z_aa-2.0E0*3*I_ESP_Gx2yz_Gy3z_a;
    abcd[iGrid*1350+1102] = 4.0E0*I_ESP_Ix3y2z_Gy3z_aa-2.0E0*1*I_ESP_Gx3y_Gy3z_a-2.0E0*2*I_ESP_Gxy2z_Gy3z_a+2*1*I_ESP_Dxy_Gy3z;
    abcd[iGrid*1350+1103] = 4.0E0*I_ESP_Ix2y3z_Gy3z_aa-2.0E0*2*I_ESP_Gx2yz_Gy3z_a-2.0E0*1*I_ESP_Gx3z_Gy3z_a+2*I_ESP_Dxz_Gy3z;
    abcd[iGrid*1350+1104] = 4.0E0*I_ESP_Ixy4z_Gy3z_aa-2.0E0*3*I_ESP_Gxy2z_Gy3z_a;
    abcd[iGrid*1350+1105] = 4.0E0*I_ESP_I5yz_Gy3z_aa-2.0E0*4*I_ESP_G3yz_Gy3z_a;
    abcd[iGrid*1350+1106] = 4.0E0*I_ESP_I4y2z_Gy3z_aa-2.0E0*1*I_ESP_G4y_Gy3z_a-2.0E0*3*I_ESP_G2y2z_Gy3z_a+3*1*I_ESP_D2y_Gy3z;
    abcd[iGrid*1350+1107] = 4.0E0*I_ESP_I3y3z_Gy3z_aa-2.0E0*2*I_ESP_G3yz_Gy3z_a-2.0E0*2*I_ESP_Gy3z_Gy3z_a+2*2*I_ESP_Dyz_Gy3z;
    abcd[iGrid*1350+1108] = 4.0E0*I_ESP_I2y4z_Gy3z_aa-2.0E0*3*I_ESP_G2y2z_Gy3z_a-2.0E0*1*I_ESP_G4z_Gy3z_a+3*I_ESP_D2z_Gy3z;
    abcd[iGrid*1350+1109] = 4.0E0*I_ESP_Iy5z_Gy3z_aa-2.0E0*4*I_ESP_Gy3z_Gy3z_a;
    abcd[iGrid*1350+1110] = 4.0E0*I_ESP_I4xyz_G4z_aa;
    abcd[iGrid*1350+1111] = 4.0E0*I_ESP_I3x2yz_G4z_aa-2.0E0*1*I_ESP_G3xz_G4z_a;
    abcd[iGrid*1350+1112] = 4.0E0*I_ESP_I3xy2z_G4z_aa-2.0E0*1*I_ESP_G3xy_G4z_a;
    abcd[iGrid*1350+1113] = 4.0E0*I_ESP_I2x3yz_G4z_aa-2.0E0*2*I_ESP_G2xyz_G4z_a;
    abcd[iGrid*1350+1114] = 4.0E0*I_ESP_I2x2y2z_G4z_aa-2.0E0*1*I_ESP_G2x2y_G4z_a-2.0E0*1*I_ESP_G2x2z_G4z_a+1*I_ESP_D2x_G4z;
    abcd[iGrid*1350+1115] = 4.0E0*I_ESP_I2xy3z_G4z_aa-2.0E0*2*I_ESP_G2xyz_G4z_a;
    abcd[iGrid*1350+1116] = 4.0E0*I_ESP_Ix4yz_G4z_aa-2.0E0*3*I_ESP_Gx2yz_G4z_a;
    abcd[iGrid*1350+1117] = 4.0E0*I_ESP_Ix3y2z_G4z_aa-2.0E0*1*I_ESP_Gx3y_G4z_a-2.0E0*2*I_ESP_Gxy2z_G4z_a+2*1*I_ESP_Dxy_G4z;
    abcd[iGrid*1350+1118] = 4.0E0*I_ESP_Ix2y3z_G4z_aa-2.0E0*2*I_ESP_Gx2yz_G4z_a-2.0E0*1*I_ESP_Gx3z_G4z_a+2*I_ESP_Dxz_G4z;
    abcd[iGrid*1350+1119] = 4.0E0*I_ESP_Ixy4z_G4z_aa-2.0E0*3*I_ESP_Gxy2z_G4z_a;
    abcd[iGrid*1350+1120] = 4.0E0*I_ESP_I5yz_G4z_aa-2.0E0*4*I_ESP_G3yz_G4z_a;
    abcd[iGrid*1350+1121] = 4.0E0*I_ESP_I4y2z_G4z_aa-2.0E0*1*I_ESP_G4y_G4z_a-2.0E0*3*I_ESP_G2y2z_G4z_a+3*1*I_ESP_D2y_G4z;
    abcd[iGrid*1350+1122] = 4.0E0*I_ESP_I3y3z_G4z_aa-2.0E0*2*I_ESP_G3yz_G4z_a-2.0E0*2*I_ESP_Gy3z_G4z_a+2*2*I_ESP_Dyz_G4z;
    abcd[iGrid*1350+1123] = 4.0E0*I_ESP_I2y4z_G4z_aa-2.0E0*3*I_ESP_G2y2z_G4z_a-2.0E0*1*I_ESP_G4z_G4z_a+3*I_ESP_D2z_G4z;
    abcd[iGrid*1350+1124] = 4.0E0*I_ESP_Iy5z_G4z_aa-2.0E0*4*I_ESP_Gy3z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G_aa
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_G_G_a
     * RHS shell quartet name: SQ_ESP_D_G
     ************************************************************/
    abcd[iGrid*1350+1125] = 4.0E0*I_ESP_I4x2z_G4x_aa-2.0E0*1*I_ESP_G4x_G4x_a;
    abcd[iGrid*1350+1126] = 4.0E0*I_ESP_I3xy2z_G4x_aa-2.0E0*1*I_ESP_G3xy_G4x_a;
    abcd[iGrid*1350+1127] = 4.0E0*I_ESP_I3x3z_G4x_aa-2.0E0*1*I_ESP_G3xz_G4x_a-2.0E0*2*I_ESP_G3xz_G4x_a;
    abcd[iGrid*1350+1128] = 4.0E0*I_ESP_I2x2y2z_G4x_aa-2.0E0*1*I_ESP_G2x2y_G4x_a;
    abcd[iGrid*1350+1129] = 4.0E0*I_ESP_I2xy3z_G4x_aa-2.0E0*1*I_ESP_G2xyz_G4x_a-2.0E0*2*I_ESP_G2xyz_G4x_a;
    abcd[iGrid*1350+1130] = 4.0E0*I_ESP_I2x4z_G4x_aa-2.0E0*2*I_ESP_G2x2z_G4x_a-2.0E0*3*I_ESP_G2x2z_G4x_a+2*1*I_ESP_D2x_G4x;
    abcd[iGrid*1350+1131] = 4.0E0*I_ESP_Ix3y2z_G4x_aa-2.0E0*1*I_ESP_Gx3y_G4x_a;
    abcd[iGrid*1350+1132] = 4.0E0*I_ESP_Ix2y3z_G4x_aa-2.0E0*1*I_ESP_Gx2yz_G4x_a-2.0E0*2*I_ESP_Gx2yz_G4x_a;
    abcd[iGrid*1350+1133] = 4.0E0*I_ESP_Ixy4z_G4x_aa-2.0E0*2*I_ESP_Gxy2z_G4x_a-2.0E0*3*I_ESP_Gxy2z_G4x_a+2*1*I_ESP_Dxy_G4x;
    abcd[iGrid*1350+1134] = 4.0E0*I_ESP_Ix5z_G4x_aa-2.0E0*3*I_ESP_Gx3z_G4x_a-2.0E0*4*I_ESP_Gx3z_G4x_a+3*2*I_ESP_Dxz_G4x;
    abcd[iGrid*1350+1135] = 4.0E0*I_ESP_I4y2z_G4x_aa-2.0E0*1*I_ESP_G4y_G4x_a;
    abcd[iGrid*1350+1136] = 4.0E0*I_ESP_I3y3z_G4x_aa-2.0E0*1*I_ESP_G3yz_G4x_a-2.0E0*2*I_ESP_G3yz_G4x_a;
    abcd[iGrid*1350+1137] = 4.0E0*I_ESP_I2y4z_G4x_aa-2.0E0*2*I_ESP_G2y2z_G4x_a-2.0E0*3*I_ESP_G2y2z_G4x_a+2*1*I_ESP_D2y_G4x;
    abcd[iGrid*1350+1138] = 4.0E0*I_ESP_Iy5z_G4x_aa-2.0E0*3*I_ESP_Gy3z_G4x_a-2.0E0*4*I_ESP_Gy3z_G4x_a+3*2*I_ESP_Dyz_G4x;
    abcd[iGrid*1350+1139] = 4.0E0*I_ESP_I6z_G4x_aa-2.0E0*4*I_ESP_G4z_G4x_a-2.0E0*5*I_ESP_G4z_G4x_a+4*3*I_ESP_D2z_G4x;
    abcd[iGrid*1350+1140] = 4.0E0*I_ESP_I4x2z_G3xy_aa-2.0E0*1*I_ESP_G4x_G3xy_a;
    abcd[iGrid*1350+1141] = 4.0E0*I_ESP_I3xy2z_G3xy_aa-2.0E0*1*I_ESP_G3xy_G3xy_a;
    abcd[iGrid*1350+1142] = 4.0E0*I_ESP_I3x3z_G3xy_aa-2.0E0*1*I_ESP_G3xz_G3xy_a-2.0E0*2*I_ESP_G3xz_G3xy_a;
    abcd[iGrid*1350+1143] = 4.0E0*I_ESP_I2x2y2z_G3xy_aa-2.0E0*1*I_ESP_G2x2y_G3xy_a;
    abcd[iGrid*1350+1144] = 4.0E0*I_ESP_I2xy3z_G3xy_aa-2.0E0*1*I_ESP_G2xyz_G3xy_a-2.0E0*2*I_ESP_G2xyz_G3xy_a;
    abcd[iGrid*1350+1145] = 4.0E0*I_ESP_I2x4z_G3xy_aa-2.0E0*2*I_ESP_G2x2z_G3xy_a-2.0E0*3*I_ESP_G2x2z_G3xy_a+2*1*I_ESP_D2x_G3xy;
    abcd[iGrid*1350+1146] = 4.0E0*I_ESP_Ix3y2z_G3xy_aa-2.0E0*1*I_ESP_Gx3y_G3xy_a;
    abcd[iGrid*1350+1147] = 4.0E0*I_ESP_Ix2y3z_G3xy_aa-2.0E0*1*I_ESP_Gx2yz_G3xy_a-2.0E0*2*I_ESP_Gx2yz_G3xy_a;
    abcd[iGrid*1350+1148] = 4.0E0*I_ESP_Ixy4z_G3xy_aa-2.0E0*2*I_ESP_Gxy2z_G3xy_a-2.0E0*3*I_ESP_Gxy2z_G3xy_a+2*1*I_ESP_Dxy_G3xy;
    abcd[iGrid*1350+1149] = 4.0E0*I_ESP_Ix5z_G3xy_aa-2.0E0*3*I_ESP_Gx3z_G3xy_a-2.0E0*4*I_ESP_Gx3z_G3xy_a+3*2*I_ESP_Dxz_G3xy;
    abcd[iGrid*1350+1150] = 4.0E0*I_ESP_I4y2z_G3xy_aa-2.0E0*1*I_ESP_G4y_G3xy_a;
    abcd[iGrid*1350+1151] = 4.0E0*I_ESP_I3y3z_G3xy_aa-2.0E0*1*I_ESP_G3yz_G3xy_a-2.0E0*2*I_ESP_G3yz_G3xy_a;
    abcd[iGrid*1350+1152] = 4.0E0*I_ESP_I2y4z_G3xy_aa-2.0E0*2*I_ESP_G2y2z_G3xy_a-2.0E0*3*I_ESP_G2y2z_G3xy_a+2*1*I_ESP_D2y_G3xy;
    abcd[iGrid*1350+1153] = 4.0E0*I_ESP_Iy5z_G3xy_aa-2.0E0*3*I_ESP_Gy3z_G3xy_a-2.0E0*4*I_ESP_Gy3z_G3xy_a+3*2*I_ESP_Dyz_G3xy;
    abcd[iGrid*1350+1154] = 4.0E0*I_ESP_I6z_G3xy_aa-2.0E0*4*I_ESP_G4z_G3xy_a-2.0E0*5*I_ESP_G4z_G3xy_a+4*3*I_ESP_D2z_G3xy;
    abcd[iGrid*1350+1155] = 4.0E0*I_ESP_I4x2z_G3xz_aa-2.0E0*1*I_ESP_G4x_G3xz_a;
    abcd[iGrid*1350+1156] = 4.0E0*I_ESP_I3xy2z_G3xz_aa-2.0E0*1*I_ESP_G3xy_G3xz_a;
    abcd[iGrid*1350+1157] = 4.0E0*I_ESP_I3x3z_G3xz_aa-2.0E0*1*I_ESP_G3xz_G3xz_a-2.0E0*2*I_ESP_G3xz_G3xz_a;
    abcd[iGrid*1350+1158] = 4.0E0*I_ESP_I2x2y2z_G3xz_aa-2.0E0*1*I_ESP_G2x2y_G3xz_a;
    abcd[iGrid*1350+1159] = 4.0E0*I_ESP_I2xy3z_G3xz_aa-2.0E0*1*I_ESP_G2xyz_G3xz_a-2.0E0*2*I_ESP_G2xyz_G3xz_a;
    abcd[iGrid*1350+1160] = 4.0E0*I_ESP_I2x4z_G3xz_aa-2.0E0*2*I_ESP_G2x2z_G3xz_a-2.0E0*3*I_ESP_G2x2z_G3xz_a+2*1*I_ESP_D2x_G3xz;
    abcd[iGrid*1350+1161] = 4.0E0*I_ESP_Ix3y2z_G3xz_aa-2.0E0*1*I_ESP_Gx3y_G3xz_a;
    abcd[iGrid*1350+1162] = 4.0E0*I_ESP_Ix2y3z_G3xz_aa-2.0E0*1*I_ESP_Gx2yz_G3xz_a-2.0E0*2*I_ESP_Gx2yz_G3xz_a;
    abcd[iGrid*1350+1163] = 4.0E0*I_ESP_Ixy4z_G3xz_aa-2.0E0*2*I_ESP_Gxy2z_G3xz_a-2.0E0*3*I_ESP_Gxy2z_G3xz_a+2*1*I_ESP_Dxy_G3xz;
    abcd[iGrid*1350+1164] = 4.0E0*I_ESP_Ix5z_G3xz_aa-2.0E0*3*I_ESP_Gx3z_G3xz_a-2.0E0*4*I_ESP_Gx3z_G3xz_a+3*2*I_ESP_Dxz_G3xz;
    abcd[iGrid*1350+1165] = 4.0E0*I_ESP_I4y2z_G3xz_aa-2.0E0*1*I_ESP_G4y_G3xz_a;
    abcd[iGrid*1350+1166] = 4.0E0*I_ESP_I3y3z_G3xz_aa-2.0E0*1*I_ESP_G3yz_G3xz_a-2.0E0*2*I_ESP_G3yz_G3xz_a;
    abcd[iGrid*1350+1167] = 4.0E0*I_ESP_I2y4z_G3xz_aa-2.0E0*2*I_ESP_G2y2z_G3xz_a-2.0E0*3*I_ESP_G2y2z_G3xz_a+2*1*I_ESP_D2y_G3xz;
    abcd[iGrid*1350+1168] = 4.0E0*I_ESP_Iy5z_G3xz_aa-2.0E0*3*I_ESP_Gy3z_G3xz_a-2.0E0*4*I_ESP_Gy3z_G3xz_a+3*2*I_ESP_Dyz_G3xz;
    abcd[iGrid*1350+1169] = 4.0E0*I_ESP_I6z_G3xz_aa-2.0E0*4*I_ESP_G4z_G3xz_a-2.0E0*5*I_ESP_G4z_G3xz_a+4*3*I_ESP_D2z_G3xz;
    abcd[iGrid*1350+1170] = 4.0E0*I_ESP_I4x2z_G2x2y_aa-2.0E0*1*I_ESP_G4x_G2x2y_a;
    abcd[iGrid*1350+1171] = 4.0E0*I_ESP_I3xy2z_G2x2y_aa-2.0E0*1*I_ESP_G3xy_G2x2y_a;
    abcd[iGrid*1350+1172] = 4.0E0*I_ESP_I3x3z_G2x2y_aa-2.0E0*1*I_ESP_G3xz_G2x2y_a-2.0E0*2*I_ESP_G3xz_G2x2y_a;
    abcd[iGrid*1350+1173] = 4.0E0*I_ESP_I2x2y2z_G2x2y_aa-2.0E0*1*I_ESP_G2x2y_G2x2y_a;
    abcd[iGrid*1350+1174] = 4.0E0*I_ESP_I2xy3z_G2x2y_aa-2.0E0*1*I_ESP_G2xyz_G2x2y_a-2.0E0*2*I_ESP_G2xyz_G2x2y_a;
    abcd[iGrid*1350+1175] = 4.0E0*I_ESP_I2x4z_G2x2y_aa-2.0E0*2*I_ESP_G2x2z_G2x2y_a-2.0E0*3*I_ESP_G2x2z_G2x2y_a+2*1*I_ESP_D2x_G2x2y;
    abcd[iGrid*1350+1176] = 4.0E0*I_ESP_Ix3y2z_G2x2y_aa-2.0E0*1*I_ESP_Gx3y_G2x2y_a;
    abcd[iGrid*1350+1177] = 4.0E0*I_ESP_Ix2y3z_G2x2y_aa-2.0E0*1*I_ESP_Gx2yz_G2x2y_a-2.0E0*2*I_ESP_Gx2yz_G2x2y_a;
    abcd[iGrid*1350+1178] = 4.0E0*I_ESP_Ixy4z_G2x2y_aa-2.0E0*2*I_ESP_Gxy2z_G2x2y_a-2.0E0*3*I_ESP_Gxy2z_G2x2y_a+2*1*I_ESP_Dxy_G2x2y;
    abcd[iGrid*1350+1179] = 4.0E0*I_ESP_Ix5z_G2x2y_aa-2.0E0*3*I_ESP_Gx3z_G2x2y_a-2.0E0*4*I_ESP_Gx3z_G2x2y_a+3*2*I_ESP_Dxz_G2x2y;
    abcd[iGrid*1350+1180] = 4.0E0*I_ESP_I4y2z_G2x2y_aa-2.0E0*1*I_ESP_G4y_G2x2y_a;
    abcd[iGrid*1350+1181] = 4.0E0*I_ESP_I3y3z_G2x2y_aa-2.0E0*1*I_ESP_G3yz_G2x2y_a-2.0E0*2*I_ESP_G3yz_G2x2y_a;
    abcd[iGrid*1350+1182] = 4.0E0*I_ESP_I2y4z_G2x2y_aa-2.0E0*2*I_ESP_G2y2z_G2x2y_a-2.0E0*3*I_ESP_G2y2z_G2x2y_a+2*1*I_ESP_D2y_G2x2y;
    abcd[iGrid*1350+1183] = 4.0E0*I_ESP_Iy5z_G2x2y_aa-2.0E0*3*I_ESP_Gy3z_G2x2y_a-2.0E0*4*I_ESP_Gy3z_G2x2y_a+3*2*I_ESP_Dyz_G2x2y;
    abcd[iGrid*1350+1184] = 4.0E0*I_ESP_I6z_G2x2y_aa-2.0E0*4*I_ESP_G4z_G2x2y_a-2.0E0*5*I_ESP_G4z_G2x2y_a+4*3*I_ESP_D2z_G2x2y;
    abcd[iGrid*1350+1185] = 4.0E0*I_ESP_I4x2z_G2xyz_aa-2.0E0*1*I_ESP_G4x_G2xyz_a;
    abcd[iGrid*1350+1186] = 4.0E0*I_ESP_I3xy2z_G2xyz_aa-2.0E0*1*I_ESP_G3xy_G2xyz_a;
    abcd[iGrid*1350+1187] = 4.0E0*I_ESP_I3x3z_G2xyz_aa-2.0E0*1*I_ESP_G3xz_G2xyz_a-2.0E0*2*I_ESP_G3xz_G2xyz_a;
    abcd[iGrid*1350+1188] = 4.0E0*I_ESP_I2x2y2z_G2xyz_aa-2.0E0*1*I_ESP_G2x2y_G2xyz_a;
    abcd[iGrid*1350+1189] = 4.0E0*I_ESP_I2xy3z_G2xyz_aa-2.0E0*1*I_ESP_G2xyz_G2xyz_a-2.0E0*2*I_ESP_G2xyz_G2xyz_a;
    abcd[iGrid*1350+1190] = 4.0E0*I_ESP_I2x4z_G2xyz_aa-2.0E0*2*I_ESP_G2x2z_G2xyz_a-2.0E0*3*I_ESP_G2x2z_G2xyz_a+2*1*I_ESP_D2x_G2xyz;
    abcd[iGrid*1350+1191] = 4.0E0*I_ESP_Ix3y2z_G2xyz_aa-2.0E0*1*I_ESP_Gx3y_G2xyz_a;
    abcd[iGrid*1350+1192] = 4.0E0*I_ESP_Ix2y3z_G2xyz_aa-2.0E0*1*I_ESP_Gx2yz_G2xyz_a-2.0E0*2*I_ESP_Gx2yz_G2xyz_a;
    abcd[iGrid*1350+1193] = 4.0E0*I_ESP_Ixy4z_G2xyz_aa-2.0E0*2*I_ESP_Gxy2z_G2xyz_a-2.0E0*3*I_ESP_Gxy2z_G2xyz_a+2*1*I_ESP_Dxy_G2xyz;
    abcd[iGrid*1350+1194] = 4.0E0*I_ESP_Ix5z_G2xyz_aa-2.0E0*3*I_ESP_Gx3z_G2xyz_a-2.0E0*4*I_ESP_Gx3z_G2xyz_a+3*2*I_ESP_Dxz_G2xyz;
    abcd[iGrid*1350+1195] = 4.0E0*I_ESP_I4y2z_G2xyz_aa-2.0E0*1*I_ESP_G4y_G2xyz_a;
    abcd[iGrid*1350+1196] = 4.0E0*I_ESP_I3y3z_G2xyz_aa-2.0E0*1*I_ESP_G3yz_G2xyz_a-2.0E0*2*I_ESP_G3yz_G2xyz_a;
    abcd[iGrid*1350+1197] = 4.0E0*I_ESP_I2y4z_G2xyz_aa-2.0E0*2*I_ESP_G2y2z_G2xyz_a-2.0E0*3*I_ESP_G2y2z_G2xyz_a+2*1*I_ESP_D2y_G2xyz;
    abcd[iGrid*1350+1198] = 4.0E0*I_ESP_Iy5z_G2xyz_aa-2.0E0*3*I_ESP_Gy3z_G2xyz_a-2.0E0*4*I_ESP_Gy3z_G2xyz_a+3*2*I_ESP_Dyz_G2xyz;
    abcd[iGrid*1350+1199] = 4.0E0*I_ESP_I6z_G2xyz_aa-2.0E0*4*I_ESP_G4z_G2xyz_a-2.0E0*5*I_ESP_G4z_G2xyz_a+4*3*I_ESP_D2z_G2xyz;
    abcd[iGrid*1350+1200] = 4.0E0*I_ESP_I4x2z_G2x2z_aa-2.0E0*1*I_ESP_G4x_G2x2z_a;
    abcd[iGrid*1350+1201] = 4.0E0*I_ESP_I3xy2z_G2x2z_aa-2.0E0*1*I_ESP_G3xy_G2x2z_a;
    abcd[iGrid*1350+1202] = 4.0E0*I_ESP_I3x3z_G2x2z_aa-2.0E0*1*I_ESP_G3xz_G2x2z_a-2.0E0*2*I_ESP_G3xz_G2x2z_a;
    abcd[iGrid*1350+1203] = 4.0E0*I_ESP_I2x2y2z_G2x2z_aa-2.0E0*1*I_ESP_G2x2y_G2x2z_a;
    abcd[iGrid*1350+1204] = 4.0E0*I_ESP_I2xy3z_G2x2z_aa-2.0E0*1*I_ESP_G2xyz_G2x2z_a-2.0E0*2*I_ESP_G2xyz_G2x2z_a;
    abcd[iGrid*1350+1205] = 4.0E0*I_ESP_I2x4z_G2x2z_aa-2.0E0*2*I_ESP_G2x2z_G2x2z_a-2.0E0*3*I_ESP_G2x2z_G2x2z_a+2*1*I_ESP_D2x_G2x2z;
    abcd[iGrid*1350+1206] = 4.0E0*I_ESP_Ix3y2z_G2x2z_aa-2.0E0*1*I_ESP_Gx3y_G2x2z_a;
    abcd[iGrid*1350+1207] = 4.0E0*I_ESP_Ix2y3z_G2x2z_aa-2.0E0*1*I_ESP_Gx2yz_G2x2z_a-2.0E0*2*I_ESP_Gx2yz_G2x2z_a;
    abcd[iGrid*1350+1208] = 4.0E0*I_ESP_Ixy4z_G2x2z_aa-2.0E0*2*I_ESP_Gxy2z_G2x2z_a-2.0E0*3*I_ESP_Gxy2z_G2x2z_a+2*1*I_ESP_Dxy_G2x2z;
    abcd[iGrid*1350+1209] = 4.0E0*I_ESP_Ix5z_G2x2z_aa-2.0E0*3*I_ESP_Gx3z_G2x2z_a-2.0E0*4*I_ESP_Gx3z_G2x2z_a+3*2*I_ESP_Dxz_G2x2z;
    abcd[iGrid*1350+1210] = 4.0E0*I_ESP_I4y2z_G2x2z_aa-2.0E0*1*I_ESP_G4y_G2x2z_a;
    abcd[iGrid*1350+1211] = 4.0E0*I_ESP_I3y3z_G2x2z_aa-2.0E0*1*I_ESP_G3yz_G2x2z_a-2.0E0*2*I_ESP_G3yz_G2x2z_a;
    abcd[iGrid*1350+1212] = 4.0E0*I_ESP_I2y4z_G2x2z_aa-2.0E0*2*I_ESP_G2y2z_G2x2z_a-2.0E0*3*I_ESP_G2y2z_G2x2z_a+2*1*I_ESP_D2y_G2x2z;
    abcd[iGrid*1350+1213] = 4.0E0*I_ESP_Iy5z_G2x2z_aa-2.0E0*3*I_ESP_Gy3z_G2x2z_a-2.0E0*4*I_ESP_Gy3z_G2x2z_a+3*2*I_ESP_Dyz_G2x2z;
    abcd[iGrid*1350+1214] = 4.0E0*I_ESP_I6z_G2x2z_aa-2.0E0*4*I_ESP_G4z_G2x2z_a-2.0E0*5*I_ESP_G4z_G2x2z_a+4*3*I_ESP_D2z_G2x2z;
    abcd[iGrid*1350+1215] = 4.0E0*I_ESP_I4x2z_Gx3y_aa-2.0E0*1*I_ESP_G4x_Gx3y_a;
    abcd[iGrid*1350+1216] = 4.0E0*I_ESP_I3xy2z_Gx3y_aa-2.0E0*1*I_ESP_G3xy_Gx3y_a;
    abcd[iGrid*1350+1217] = 4.0E0*I_ESP_I3x3z_Gx3y_aa-2.0E0*1*I_ESP_G3xz_Gx3y_a-2.0E0*2*I_ESP_G3xz_Gx3y_a;
    abcd[iGrid*1350+1218] = 4.0E0*I_ESP_I2x2y2z_Gx3y_aa-2.0E0*1*I_ESP_G2x2y_Gx3y_a;
    abcd[iGrid*1350+1219] = 4.0E0*I_ESP_I2xy3z_Gx3y_aa-2.0E0*1*I_ESP_G2xyz_Gx3y_a-2.0E0*2*I_ESP_G2xyz_Gx3y_a;
    abcd[iGrid*1350+1220] = 4.0E0*I_ESP_I2x4z_Gx3y_aa-2.0E0*2*I_ESP_G2x2z_Gx3y_a-2.0E0*3*I_ESP_G2x2z_Gx3y_a+2*1*I_ESP_D2x_Gx3y;
    abcd[iGrid*1350+1221] = 4.0E0*I_ESP_Ix3y2z_Gx3y_aa-2.0E0*1*I_ESP_Gx3y_Gx3y_a;
    abcd[iGrid*1350+1222] = 4.0E0*I_ESP_Ix2y3z_Gx3y_aa-2.0E0*1*I_ESP_Gx2yz_Gx3y_a-2.0E0*2*I_ESP_Gx2yz_Gx3y_a;
    abcd[iGrid*1350+1223] = 4.0E0*I_ESP_Ixy4z_Gx3y_aa-2.0E0*2*I_ESP_Gxy2z_Gx3y_a-2.0E0*3*I_ESP_Gxy2z_Gx3y_a+2*1*I_ESP_Dxy_Gx3y;
    abcd[iGrid*1350+1224] = 4.0E0*I_ESP_Ix5z_Gx3y_aa-2.0E0*3*I_ESP_Gx3z_Gx3y_a-2.0E0*4*I_ESP_Gx3z_Gx3y_a+3*2*I_ESP_Dxz_Gx3y;
    abcd[iGrid*1350+1225] = 4.0E0*I_ESP_I4y2z_Gx3y_aa-2.0E0*1*I_ESP_G4y_Gx3y_a;
    abcd[iGrid*1350+1226] = 4.0E0*I_ESP_I3y3z_Gx3y_aa-2.0E0*1*I_ESP_G3yz_Gx3y_a-2.0E0*2*I_ESP_G3yz_Gx3y_a;
    abcd[iGrid*1350+1227] = 4.0E0*I_ESP_I2y4z_Gx3y_aa-2.0E0*2*I_ESP_G2y2z_Gx3y_a-2.0E0*3*I_ESP_G2y2z_Gx3y_a+2*1*I_ESP_D2y_Gx3y;
    abcd[iGrid*1350+1228] = 4.0E0*I_ESP_Iy5z_Gx3y_aa-2.0E0*3*I_ESP_Gy3z_Gx3y_a-2.0E0*4*I_ESP_Gy3z_Gx3y_a+3*2*I_ESP_Dyz_Gx3y;
    abcd[iGrid*1350+1229] = 4.0E0*I_ESP_I6z_Gx3y_aa-2.0E0*4*I_ESP_G4z_Gx3y_a-2.0E0*5*I_ESP_G4z_Gx3y_a+4*3*I_ESP_D2z_Gx3y;
    abcd[iGrid*1350+1230] = 4.0E0*I_ESP_I4x2z_Gx2yz_aa-2.0E0*1*I_ESP_G4x_Gx2yz_a;
    abcd[iGrid*1350+1231] = 4.0E0*I_ESP_I3xy2z_Gx2yz_aa-2.0E0*1*I_ESP_G3xy_Gx2yz_a;
    abcd[iGrid*1350+1232] = 4.0E0*I_ESP_I3x3z_Gx2yz_aa-2.0E0*1*I_ESP_G3xz_Gx2yz_a-2.0E0*2*I_ESP_G3xz_Gx2yz_a;
    abcd[iGrid*1350+1233] = 4.0E0*I_ESP_I2x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_G2x2y_Gx2yz_a;
    abcd[iGrid*1350+1234] = 4.0E0*I_ESP_I2xy3z_Gx2yz_aa-2.0E0*1*I_ESP_G2xyz_Gx2yz_a-2.0E0*2*I_ESP_G2xyz_Gx2yz_a;
    abcd[iGrid*1350+1235] = 4.0E0*I_ESP_I2x4z_Gx2yz_aa-2.0E0*2*I_ESP_G2x2z_Gx2yz_a-2.0E0*3*I_ESP_G2x2z_Gx2yz_a+2*1*I_ESP_D2x_Gx2yz;
    abcd[iGrid*1350+1236] = 4.0E0*I_ESP_Ix3y2z_Gx2yz_aa-2.0E0*1*I_ESP_Gx3y_Gx2yz_a;
    abcd[iGrid*1350+1237] = 4.0E0*I_ESP_Ix2y3z_Gx2yz_aa-2.0E0*1*I_ESP_Gx2yz_Gx2yz_a-2.0E0*2*I_ESP_Gx2yz_Gx2yz_a;
    abcd[iGrid*1350+1238] = 4.0E0*I_ESP_Ixy4z_Gx2yz_aa-2.0E0*2*I_ESP_Gxy2z_Gx2yz_a-2.0E0*3*I_ESP_Gxy2z_Gx2yz_a+2*1*I_ESP_Dxy_Gx2yz;
    abcd[iGrid*1350+1239] = 4.0E0*I_ESP_Ix5z_Gx2yz_aa-2.0E0*3*I_ESP_Gx3z_Gx2yz_a-2.0E0*4*I_ESP_Gx3z_Gx2yz_a+3*2*I_ESP_Dxz_Gx2yz;
    abcd[iGrid*1350+1240] = 4.0E0*I_ESP_I4y2z_Gx2yz_aa-2.0E0*1*I_ESP_G4y_Gx2yz_a;
    abcd[iGrid*1350+1241] = 4.0E0*I_ESP_I3y3z_Gx2yz_aa-2.0E0*1*I_ESP_G3yz_Gx2yz_a-2.0E0*2*I_ESP_G3yz_Gx2yz_a;
    abcd[iGrid*1350+1242] = 4.0E0*I_ESP_I2y4z_Gx2yz_aa-2.0E0*2*I_ESP_G2y2z_Gx2yz_a-2.0E0*3*I_ESP_G2y2z_Gx2yz_a+2*1*I_ESP_D2y_Gx2yz;
    abcd[iGrid*1350+1243] = 4.0E0*I_ESP_Iy5z_Gx2yz_aa-2.0E0*3*I_ESP_Gy3z_Gx2yz_a-2.0E0*4*I_ESP_Gy3z_Gx2yz_a+3*2*I_ESP_Dyz_Gx2yz;
    abcd[iGrid*1350+1244] = 4.0E0*I_ESP_I6z_Gx2yz_aa-2.0E0*4*I_ESP_G4z_Gx2yz_a-2.0E0*5*I_ESP_G4z_Gx2yz_a+4*3*I_ESP_D2z_Gx2yz;
    abcd[iGrid*1350+1245] = 4.0E0*I_ESP_I4x2z_Gxy2z_aa-2.0E0*1*I_ESP_G4x_Gxy2z_a;
    abcd[iGrid*1350+1246] = 4.0E0*I_ESP_I3xy2z_Gxy2z_aa-2.0E0*1*I_ESP_G3xy_Gxy2z_a;
    abcd[iGrid*1350+1247] = 4.0E0*I_ESP_I3x3z_Gxy2z_aa-2.0E0*1*I_ESP_G3xz_Gxy2z_a-2.0E0*2*I_ESP_G3xz_Gxy2z_a;
    abcd[iGrid*1350+1248] = 4.0E0*I_ESP_I2x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_G2x2y_Gxy2z_a;
    abcd[iGrid*1350+1249] = 4.0E0*I_ESP_I2xy3z_Gxy2z_aa-2.0E0*1*I_ESP_G2xyz_Gxy2z_a-2.0E0*2*I_ESP_G2xyz_Gxy2z_a;
    abcd[iGrid*1350+1250] = 4.0E0*I_ESP_I2x4z_Gxy2z_aa-2.0E0*2*I_ESP_G2x2z_Gxy2z_a-2.0E0*3*I_ESP_G2x2z_Gxy2z_a+2*1*I_ESP_D2x_Gxy2z;
    abcd[iGrid*1350+1251] = 4.0E0*I_ESP_Ix3y2z_Gxy2z_aa-2.0E0*1*I_ESP_Gx3y_Gxy2z_a;
    abcd[iGrid*1350+1252] = 4.0E0*I_ESP_Ix2y3z_Gxy2z_aa-2.0E0*1*I_ESP_Gx2yz_Gxy2z_a-2.0E0*2*I_ESP_Gx2yz_Gxy2z_a;
    abcd[iGrid*1350+1253] = 4.0E0*I_ESP_Ixy4z_Gxy2z_aa-2.0E0*2*I_ESP_Gxy2z_Gxy2z_a-2.0E0*3*I_ESP_Gxy2z_Gxy2z_a+2*1*I_ESP_Dxy_Gxy2z;
    abcd[iGrid*1350+1254] = 4.0E0*I_ESP_Ix5z_Gxy2z_aa-2.0E0*3*I_ESP_Gx3z_Gxy2z_a-2.0E0*4*I_ESP_Gx3z_Gxy2z_a+3*2*I_ESP_Dxz_Gxy2z;
    abcd[iGrid*1350+1255] = 4.0E0*I_ESP_I4y2z_Gxy2z_aa-2.0E0*1*I_ESP_G4y_Gxy2z_a;
    abcd[iGrid*1350+1256] = 4.0E0*I_ESP_I3y3z_Gxy2z_aa-2.0E0*1*I_ESP_G3yz_Gxy2z_a-2.0E0*2*I_ESP_G3yz_Gxy2z_a;
    abcd[iGrid*1350+1257] = 4.0E0*I_ESP_I2y4z_Gxy2z_aa-2.0E0*2*I_ESP_G2y2z_Gxy2z_a-2.0E0*3*I_ESP_G2y2z_Gxy2z_a+2*1*I_ESP_D2y_Gxy2z;
    abcd[iGrid*1350+1258] = 4.0E0*I_ESP_Iy5z_Gxy2z_aa-2.0E0*3*I_ESP_Gy3z_Gxy2z_a-2.0E0*4*I_ESP_Gy3z_Gxy2z_a+3*2*I_ESP_Dyz_Gxy2z;
    abcd[iGrid*1350+1259] = 4.0E0*I_ESP_I6z_Gxy2z_aa-2.0E0*4*I_ESP_G4z_Gxy2z_a-2.0E0*5*I_ESP_G4z_Gxy2z_a+4*3*I_ESP_D2z_Gxy2z;
    abcd[iGrid*1350+1260] = 4.0E0*I_ESP_I4x2z_Gx3z_aa-2.0E0*1*I_ESP_G4x_Gx3z_a;
    abcd[iGrid*1350+1261] = 4.0E0*I_ESP_I3xy2z_Gx3z_aa-2.0E0*1*I_ESP_G3xy_Gx3z_a;
    abcd[iGrid*1350+1262] = 4.0E0*I_ESP_I3x3z_Gx3z_aa-2.0E0*1*I_ESP_G3xz_Gx3z_a-2.0E0*2*I_ESP_G3xz_Gx3z_a;
    abcd[iGrid*1350+1263] = 4.0E0*I_ESP_I2x2y2z_Gx3z_aa-2.0E0*1*I_ESP_G2x2y_Gx3z_a;
    abcd[iGrid*1350+1264] = 4.0E0*I_ESP_I2xy3z_Gx3z_aa-2.0E0*1*I_ESP_G2xyz_Gx3z_a-2.0E0*2*I_ESP_G2xyz_Gx3z_a;
    abcd[iGrid*1350+1265] = 4.0E0*I_ESP_I2x4z_Gx3z_aa-2.0E0*2*I_ESP_G2x2z_Gx3z_a-2.0E0*3*I_ESP_G2x2z_Gx3z_a+2*1*I_ESP_D2x_Gx3z;
    abcd[iGrid*1350+1266] = 4.0E0*I_ESP_Ix3y2z_Gx3z_aa-2.0E0*1*I_ESP_Gx3y_Gx3z_a;
    abcd[iGrid*1350+1267] = 4.0E0*I_ESP_Ix2y3z_Gx3z_aa-2.0E0*1*I_ESP_Gx2yz_Gx3z_a-2.0E0*2*I_ESP_Gx2yz_Gx3z_a;
    abcd[iGrid*1350+1268] = 4.0E0*I_ESP_Ixy4z_Gx3z_aa-2.0E0*2*I_ESP_Gxy2z_Gx3z_a-2.0E0*3*I_ESP_Gxy2z_Gx3z_a+2*1*I_ESP_Dxy_Gx3z;
    abcd[iGrid*1350+1269] = 4.0E0*I_ESP_Ix5z_Gx3z_aa-2.0E0*3*I_ESP_Gx3z_Gx3z_a-2.0E0*4*I_ESP_Gx3z_Gx3z_a+3*2*I_ESP_Dxz_Gx3z;
    abcd[iGrid*1350+1270] = 4.0E0*I_ESP_I4y2z_Gx3z_aa-2.0E0*1*I_ESP_G4y_Gx3z_a;
    abcd[iGrid*1350+1271] = 4.0E0*I_ESP_I3y3z_Gx3z_aa-2.0E0*1*I_ESP_G3yz_Gx3z_a-2.0E0*2*I_ESP_G3yz_Gx3z_a;
    abcd[iGrid*1350+1272] = 4.0E0*I_ESP_I2y4z_Gx3z_aa-2.0E0*2*I_ESP_G2y2z_Gx3z_a-2.0E0*3*I_ESP_G2y2z_Gx3z_a+2*1*I_ESP_D2y_Gx3z;
    abcd[iGrid*1350+1273] = 4.0E0*I_ESP_Iy5z_Gx3z_aa-2.0E0*3*I_ESP_Gy3z_Gx3z_a-2.0E0*4*I_ESP_Gy3z_Gx3z_a+3*2*I_ESP_Dyz_Gx3z;
    abcd[iGrid*1350+1274] = 4.0E0*I_ESP_I6z_Gx3z_aa-2.0E0*4*I_ESP_G4z_Gx3z_a-2.0E0*5*I_ESP_G4z_Gx3z_a+4*3*I_ESP_D2z_Gx3z;
    abcd[iGrid*1350+1275] = 4.0E0*I_ESP_I4x2z_G4y_aa-2.0E0*1*I_ESP_G4x_G4y_a;
    abcd[iGrid*1350+1276] = 4.0E0*I_ESP_I3xy2z_G4y_aa-2.0E0*1*I_ESP_G3xy_G4y_a;
    abcd[iGrid*1350+1277] = 4.0E0*I_ESP_I3x3z_G4y_aa-2.0E0*1*I_ESP_G3xz_G4y_a-2.0E0*2*I_ESP_G3xz_G4y_a;
    abcd[iGrid*1350+1278] = 4.0E0*I_ESP_I2x2y2z_G4y_aa-2.0E0*1*I_ESP_G2x2y_G4y_a;
    abcd[iGrid*1350+1279] = 4.0E0*I_ESP_I2xy3z_G4y_aa-2.0E0*1*I_ESP_G2xyz_G4y_a-2.0E0*2*I_ESP_G2xyz_G4y_a;
    abcd[iGrid*1350+1280] = 4.0E0*I_ESP_I2x4z_G4y_aa-2.0E0*2*I_ESP_G2x2z_G4y_a-2.0E0*3*I_ESP_G2x2z_G4y_a+2*1*I_ESP_D2x_G4y;
    abcd[iGrid*1350+1281] = 4.0E0*I_ESP_Ix3y2z_G4y_aa-2.0E0*1*I_ESP_Gx3y_G4y_a;
    abcd[iGrid*1350+1282] = 4.0E0*I_ESP_Ix2y3z_G4y_aa-2.0E0*1*I_ESP_Gx2yz_G4y_a-2.0E0*2*I_ESP_Gx2yz_G4y_a;
    abcd[iGrid*1350+1283] = 4.0E0*I_ESP_Ixy4z_G4y_aa-2.0E0*2*I_ESP_Gxy2z_G4y_a-2.0E0*3*I_ESP_Gxy2z_G4y_a+2*1*I_ESP_Dxy_G4y;
    abcd[iGrid*1350+1284] = 4.0E0*I_ESP_Ix5z_G4y_aa-2.0E0*3*I_ESP_Gx3z_G4y_a-2.0E0*4*I_ESP_Gx3z_G4y_a+3*2*I_ESP_Dxz_G4y;
    abcd[iGrid*1350+1285] = 4.0E0*I_ESP_I4y2z_G4y_aa-2.0E0*1*I_ESP_G4y_G4y_a;
    abcd[iGrid*1350+1286] = 4.0E0*I_ESP_I3y3z_G4y_aa-2.0E0*1*I_ESP_G3yz_G4y_a-2.0E0*2*I_ESP_G3yz_G4y_a;
    abcd[iGrid*1350+1287] = 4.0E0*I_ESP_I2y4z_G4y_aa-2.0E0*2*I_ESP_G2y2z_G4y_a-2.0E0*3*I_ESP_G2y2z_G4y_a+2*1*I_ESP_D2y_G4y;
    abcd[iGrid*1350+1288] = 4.0E0*I_ESP_Iy5z_G4y_aa-2.0E0*3*I_ESP_Gy3z_G4y_a-2.0E0*4*I_ESP_Gy3z_G4y_a+3*2*I_ESP_Dyz_G4y;
    abcd[iGrid*1350+1289] = 4.0E0*I_ESP_I6z_G4y_aa-2.0E0*4*I_ESP_G4z_G4y_a-2.0E0*5*I_ESP_G4z_G4y_a+4*3*I_ESP_D2z_G4y;
    abcd[iGrid*1350+1290] = 4.0E0*I_ESP_I4x2z_G3yz_aa-2.0E0*1*I_ESP_G4x_G3yz_a;
    abcd[iGrid*1350+1291] = 4.0E0*I_ESP_I3xy2z_G3yz_aa-2.0E0*1*I_ESP_G3xy_G3yz_a;
    abcd[iGrid*1350+1292] = 4.0E0*I_ESP_I3x3z_G3yz_aa-2.0E0*1*I_ESP_G3xz_G3yz_a-2.0E0*2*I_ESP_G3xz_G3yz_a;
    abcd[iGrid*1350+1293] = 4.0E0*I_ESP_I2x2y2z_G3yz_aa-2.0E0*1*I_ESP_G2x2y_G3yz_a;
    abcd[iGrid*1350+1294] = 4.0E0*I_ESP_I2xy3z_G3yz_aa-2.0E0*1*I_ESP_G2xyz_G3yz_a-2.0E0*2*I_ESP_G2xyz_G3yz_a;
    abcd[iGrid*1350+1295] = 4.0E0*I_ESP_I2x4z_G3yz_aa-2.0E0*2*I_ESP_G2x2z_G3yz_a-2.0E0*3*I_ESP_G2x2z_G3yz_a+2*1*I_ESP_D2x_G3yz;
    abcd[iGrid*1350+1296] = 4.0E0*I_ESP_Ix3y2z_G3yz_aa-2.0E0*1*I_ESP_Gx3y_G3yz_a;
    abcd[iGrid*1350+1297] = 4.0E0*I_ESP_Ix2y3z_G3yz_aa-2.0E0*1*I_ESP_Gx2yz_G3yz_a-2.0E0*2*I_ESP_Gx2yz_G3yz_a;
    abcd[iGrid*1350+1298] = 4.0E0*I_ESP_Ixy4z_G3yz_aa-2.0E0*2*I_ESP_Gxy2z_G3yz_a-2.0E0*3*I_ESP_Gxy2z_G3yz_a+2*1*I_ESP_Dxy_G3yz;
    abcd[iGrid*1350+1299] = 4.0E0*I_ESP_Ix5z_G3yz_aa-2.0E0*3*I_ESP_Gx3z_G3yz_a-2.0E0*4*I_ESP_Gx3z_G3yz_a+3*2*I_ESP_Dxz_G3yz;
    abcd[iGrid*1350+1300] = 4.0E0*I_ESP_I4y2z_G3yz_aa-2.0E0*1*I_ESP_G4y_G3yz_a;
    abcd[iGrid*1350+1301] = 4.0E0*I_ESP_I3y3z_G3yz_aa-2.0E0*1*I_ESP_G3yz_G3yz_a-2.0E0*2*I_ESP_G3yz_G3yz_a;
    abcd[iGrid*1350+1302] = 4.0E0*I_ESP_I2y4z_G3yz_aa-2.0E0*2*I_ESP_G2y2z_G3yz_a-2.0E0*3*I_ESP_G2y2z_G3yz_a+2*1*I_ESP_D2y_G3yz;
    abcd[iGrid*1350+1303] = 4.0E0*I_ESP_Iy5z_G3yz_aa-2.0E0*3*I_ESP_Gy3z_G3yz_a-2.0E0*4*I_ESP_Gy3z_G3yz_a+3*2*I_ESP_Dyz_G3yz;
    abcd[iGrid*1350+1304] = 4.0E0*I_ESP_I6z_G3yz_aa-2.0E0*4*I_ESP_G4z_G3yz_a-2.0E0*5*I_ESP_G4z_G3yz_a+4*3*I_ESP_D2z_G3yz;
    abcd[iGrid*1350+1305] = 4.0E0*I_ESP_I4x2z_G2y2z_aa-2.0E0*1*I_ESP_G4x_G2y2z_a;
    abcd[iGrid*1350+1306] = 4.0E0*I_ESP_I3xy2z_G2y2z_aa-2.0E0*1*I_ESP_G3xy_G2y2z_a;
    abcd[iGrid*1350+1307] = 4.0E0*I_ESP_I3x3z_G2y2z_aa-2.0E0*1*I_ESP_G3xz_G2y2z_a-2.0E0*2*I_ESP_G3xz_G2y2z_a;
    abcd[iGrid*1350+1308] = 4.0E0*I_ESP_I2x2y2z_G2y2z_aa-2.0E0*1*I_ESP_G2x2y_G2y2z_a;
    abcd[iGrid*1350+1309] = 4.0E0*I_ESP_I2xy3z_G2y2z_aa-2.0E0*1*I_ESP_G2xyz_G2y2z_a-2.0E0*2*I_ESP_G2xyz_G2y2z_a;
    abcd[iGrid*1350+1310] = 4.0E0*I_ESP_I2x4z_G2y2z_aa-2.0E0*2*I_ESP_G2x2z_G2y2z_a-2.0E0*3*I_ESP_G2x2z_G2y2z_a+2*1*I_ESP_D2x_G2y2z;
    abcd[iGrid*1350+1311] = 4.0E0*I_ESP_Ix3y2z_G2y2z_aa-2.0E0*1*I_ESP_Gx3y_G2y2z_a;
    abcd[iGrid*1350+1312] = 4.0E0*I_ESP_Ix2y3z_G2y2z_aa-2.0E0*1*I_ESP_Gx2yz_G2y2z_a-2.0E0*2*I_ESP_Gx2yz_G2y2z_a;
    abcd[iGrid*1350+1313] = 4.0E0*I_ESP_Ixy4z_G2y2z_aa-2.0E0*2*I_ESP_Gxy2z_G2y2z_a-2.0E0*3*I_ESP_Gxy2z_G2y2z_a+2*1*I_ESP_Dxy_G2y2z;
    abcd[iGrid*1350+1314] = 4.0E0*I_ESP_Ix5z_G2y2z_aa-2.0E0*3*I_ESP_Gx3z_G2y2z_a-2.0E0*4*I_ESP_Gx3z_G2y2z_a+3*2*I_ESP_Dxz_G2y2z;
    abcd[iGrid*1350+1315] = 4.0E0*I_ESP_I4y2z_G2y2z_aa-2.0E0*1*I_ESP_G4y_G2y2z_a;
    abcd[iGrid*1350+1316] = 4.0E0*I_ESP_I3y3z_G2y2z_aa-2.0E0*1*I_ESP_G3yz_G2y2z_a-2.0E0*2*I_ESP_G3yz_G2y2z_a;
    abcd[iGrid*1350+1317] = 4.0E0*I_ESP_I2y4z_G2y2z_aa-2.0E0*2*I_ESP_G2y2z_G2y2z_a-2.0E0*3*I_ESP_G2y2z_G2y2z_a+2*1*I_ESP_D2y_G2y2z;
    abcd[iGrid*1350+1318] = 4.0E0*I_ESP_Iy5z_G2y2z_aa-2.0E0*3*I_ESP_Gy3z_G2y2z_a-2.0E0*4*I_ESP_Gy3z_G2y2z_a+3*2*I_ESP_Dyz_G2y2z;
    abcd[iGrid*1350+1319] = 4.0E0*I_ESP_I6z_G2y2z_aa-2.0E0*4*I_ESP_G4z_G2y2z_a-2.0E0*5*I_ESP_G4z_G2y2z_a+4*3*I_ESP_D2z_G2y2z;
    abcd[iGrid*1350+1320] = 4.0E0*I_ESP_I4x2z_Gy3z_aa-2.0E0*1*I_ESP_G4x_Gy3z_a;
    abcd[iGrid*1350+1321] = 4.0E0*I_ESP_I3xy2z_Gy3z_aa-2.0E0*1*I_ESP_G3xy_Gy3z_a;
    abcd[iGrid*1350+1322] = 4.0E0*I_ESP_I3x3z_Gy3z_aa-2.0E0*1*I_ESP_G3xz_Gy3z_a-2.0E0*2*I_ESP_G3xz_Gy3z_a;
    abcd[iGrid*1350+1323] = 4.0E0*I_ESP_I2x2y2z_Gy3z_aa-2.0E0*1*I_ESP_G2x2y_Gy3z_a;
    abcd[iGrid*1350+1324] = 4.0E0*I_ESP_I2xy3z_Gy3z_aa-2.0E0*1*I_ESP_G2xyz_Gy3z_a-2.0E0*2*I_ESP_G2xyz_Gy3z_a;
    abcd[iGrid*1350+1325] = 4.0E0*I_ESP_I2x4z_Gy3z_aa-2.0E0*2*I_ESP_G2x2z_Gy3z_a-2.0E0*3*I_ESP_G2x2z_Gy3z_a+2*1*I_ESP_D2x_Gy3z;
    abcd[iGrid*1350+1326] = 4.0E0*I_ESP_Ix3y2z_Gy3z_aa-2.0E0*1*I_ESP_Gx3y_Gy3z_a;
    abcd[iGrid*1350+1327] = 4.0E0*I_ESP_Ix2y3z_Gy3z_aa-2.0E0*1*I_ESP_Gx2yz_Gy3z_a-2.0E0*2*I_ESP_Gx2yz_Gy3z_a;
    abcd[iGrid*1350+1328] = 4.0E0*I_ESP_Ixy4z_Gy3z_aa-2.0E0*2*I_ESP_Gxy2z_Gy3z_a-2.0E0*3*I_ESP_Gxy2z_Gy3z_a+2*1*I_ESP_Dxy_Gy3z;
    abcd[iGrid*1350+1329] = 4.0E0*I_ESP_Ix5z_Gy3z_aa-2.0E0*3*I_ESP_Gx3z_Gy3z_a-2.0E0*4*I_ESP_Gx3z_Gy3z_a+3*2*I_ESP_Dxz_Gy3z;
    abcd[iGrid*1350+1330] = 4.0E0*I_ESP_I4y2z_Gy3z_aa-2.0E0*1*I_ESP_G4y_Gy3z_a;
    abcd[iGrid*1350+1331] = 4.0E0*I_ESP_I3y3z_Gy3z_aa-2.0E0*1*I_ESP_G3yz_Gy3z_a-2.0E0*2*I_ESP_G3yz_Gy3z_a;
    abcd[iGrid*1350+1332] = 4.0E0*I_ESP_I2y4z_Gy3z_aa-2.0E0*2*I_ESP_G2y2z_Gy3z_a-2.0E0*3*I_ESP_G2y2z_Gy3z_a+2*1*I_ESP_D2y_Gy3z;
    abcd[iGrid*1350+1333] = 4.0E0*I_ESP_Iy5z_Gy3z_aa-2.0E0*3*I_ESP_Gy3z_Gy3z_a-2.0E0*4*I_ESP_Gy3z_Gy3z_a+3*2*I_ESP_Dyz_Gy3z;
    abcd[iGrid*1350+1334] = 4.0E0*I_ESP_I6z_Gy3z_aa-2.0E0*4*I_ESP_G4z_Gy3z_a-2.0E0*5*I_ESP_G4z_Gy3z_a+4*3*I_ESP_D2z_Gy3z;
    abcd[iGrid*1350+1335] = 4.0E0*I_ESP_I4x2z_G4z_aa-2.0E0*1*I_ESP_G4x_G4z_a;
    abcd[iGrid*1350+1336] = 4.0E0*I_ESP_I3xy2z_G4z_aa-2.0E0*1*I_ESP_G3xy_G4z_a;
    abcd[iGrid*1350+1337] = 4.0E0*I_ESP_I3x3z_G4z_aa-2.0E0*1*I_ESP_G3xz_G4z_a-2.0E0*2*I_ESP_G3xz_G4z_a;
    abcd[iGrid*1350+1338] = 4.0E0*I_ESP_I2x2y2z_G4z_aa-2.0E0*1*I_ESP_G2x2y_G4z_a;
    abcd[iGrid*1350+1339] = 4.0E0*I_ESP_I2xy3z_G4z_aa-2.0E0*1*I_ESP_G2xyz_G4z_a-2.0E0*2*I_ESP_G2xyz_G4z_a;
    abcd[iGrid*1350+1340] = 4.0E0*I_ESP_I2x4z_G4z_aa-2.0E0*2*I_ESP_G2x2z_G4z_a-2.0E0*3*I_ESP_G2x2z_G4z_a+2*1*I_ESP_D2x_G4z;
    abcd[iGrid*1350+1341] = 4.0E0*I_ESP_Ix3y2z_G4z_aa-2.0E0*1*I_ESP_Gx3y_G4z_a;
    abcd[iGrid*1350+1342] = 4.0E0*I_ESP_Ix2y3z_G4z_aa-2.0E0*1*I_ESP_Gx2yz_G4z_a-2.0E0*2*I_ESP_Gx2yz_G4z_a;
    abcd[iGrid*1350+1343] = 4.0E0*I_ESP_Ixy4z_G4z_aa-2.0E0*2*I_ESP_Gxy2z_G4z_a-2.0E0*3*I_ESP_Gxy2z_G4z_a+2*1*I_ESP_Dxy_G4z;
    abcd[iGrid*1350+1344] = 4.0E0*I_ESP_Ix5z_G4z_aa-2.0E0*3*I_ESP_Gx3z_G4z_a-2.0E0*4*I_ESP_Gx3z_G4z_a+3*2*I_ESP_Dxz_G4z;
    abcd[iGrid*1350+1345] = 4.0E0*I_ESP_I4y2z_G4z_aa-2.0E0*1*I_ESP_G4y_G4z_a;
    abcd[iGrid*1350+1346] = 4.0E0*I_ESP_I3y3z_G4z_aa-2.0E0*1*I_ESP_G3yz_G4z_a-2.0E0*2*I_ESP_G3yz_G4z_a;
    abcd[iGrid*1350+1347] = 4.0E0*I_ESP_I2y4z_G4z_aa-2.0E0*2*I_ESP_G2y2z_G4z_a-2.0E0*3*I_ESP_G2y2z_G4z_a+2*1*I_ESP_D2y_G4z;
    abcd[iGrid*1350+1348] = 4.0E0*I_ESP_Iy5z_G4z_aa-2.0E0*3*I_ESP_Gy3z_G4z_a-2.0E0*4*I_ESP_Gy3z_G4z_a+3*2*I_ESP_Dyz_G4z;
    abcd[iGrid*1350+1349] = 4.0E0*I_ESP_I6z_G4z_aa-2.0E0*4*I_ESP_G4z_G4z_a-2.0E0*5*I_ESP_G4z_G4z_a+4*3*I_ESP_D2z_G4z;
  }
}
