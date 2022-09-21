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
// BRA1 as redundant position, total RHS integrals evaluated as: 9661
// BRA2 as redundant position, total RHS integrals evaluated as: 9062
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: BRA2
//

//
// @@@@ derivative position-direction information
// BRA1
// X
// Y
// Z
// ####

void hgp_os_esp_h_g_d1(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_N10x_S_a = 0.0E0;
    Double I_ESP_N9xy_S_a = 0.0E0;
    Double I_ESP_N9xz_S_a = 0.0E0;
    Double I_ESP_N8x2y_S_a = 0.0E0;
    Double I_ESP_N8xyz_S_a = 0.0E0;
    Double I_ESP_N8x2z_S_a = 0.0E0;
    Double I_ESP_N7x3y_S_a = 0.0E0;
    Double I_ESP_N7x2yz_S_a = 0.0E0;
    Double I_ESP_N7xy2z_S_a = 0.0E0;
    Double I_ESP_N7x3z_S_a = 0.0E0;
    Double I_ESP_N6x4y_S_a = 0.0E0;
    Double I_ESP_N6x3yz_S_a = 0.0E0;
    Double I_ESP_N6x2y2z_S_a = 0.0E0;
    Double I_ESP_N6xy3z_S_a = 0.0E0;
    Double I_ESP_N6x4z_S_a = 0.0E0;
    Double I_ESP_N5x5y_S_a = 0.0E0;
    Double I_ESP_N5x4yz_S_a = 0.0E0;
    Double I_ESP_N5x3y2z_S_a = 0.0E0;
    Double I_ESP_N5x2y3z_S_a = 0.0E0;
    Double I_ESP_N5xy4z_S_a = 0.0E0;
    Double I_ESP_N5x5z_S_a = 0.0E0;
    Double I_ESP_N4x6y_S_a = 0.0E0;
    Double I_ESP_N4x5yz_S_a = 0.0E0;
    Double I_ESP_N4x4y2z_S_a = 0.0E0;
    Double I_ESP_N4x3y3z_S_a = 0.0E0;
    Double I_ESP_N4x2y4z_S_a = 0.0E0;
    Double I_ESP_N4xy5z_S_a = 0.0E0;
    Double I_ESP_N4x6z_S_a = 0.0E0;
    Double I_ESP_N3x7y_S_a = 0.0E0;
    Double I_ESP_N3x6yz_S_a = 0.0E0;
    Double I_ESP_N3x5y2z_S_a = 0.0E0;
    Double I_ESP_N3x4y3z_S_a = 0.0E0;
    Double I_ESP_N3x3y4z_S_a = 0.0E0;
    Double I_ESP_N3x2y5z_S_a = 0.0E0;
    Double I_ESP_N3xy6z_S_a = 0.0E0;
    Double I_ESP_N3x7z_S_a = 0.0E0;
    Double I_ESP_N2x8y_S_a = 0.0E0;
    Double I_ESP_N2x7yz_S_a = 0.0E0;
    Double I_ESP_N2x6y2z_S_a = 0.0E0;
    Double I_ESP_N2x5y3z_S_a = 0.0E0;
    Double I_ESP_N2x4y4z_S_a = 0.0E0;
    Double I_ESP_N2x3y5z_S_a = 0.0E0;
    Double I_ESP_N2x2y6z_S_a = 0.0E0;
    Double I_ESP_N2xy7z_S_a = 0.0E0;
    Double I_ESP_N2x8z_S_a = 0.0E0;
    Double I_ESP_Nx9y_S_a = 0.0E0;
    Double I_ESP_Nx8yz_S_a = 0.0E0;
    Double I_ESP_Nx7y2z_S_a = 0.0E0;
    Double I_ESP_Nx6y3z_S_a = 0.0E0;
    Double I_ESP_Nx5y4z_S_a = 0.0E0;
    Double I_ESP_Nx4y5z_S_a = 0.0E0;
    Double I_ESP_Nx3y6z_S_a = 0.0E0;
    Double I_ESP_Nx2y7z_S_a = 0.0E0;
    Double I_ESP_Nxy8z_S_a = 0.0E0;
    Double I_ESP_Nx9z_S_a = 0.0E0;
    Double I_ESP_N10y_S_a = 0.0E0;
    Double I_ESP_N9yz_S_a = 0.0E0;
    Double I_ESP_N8y2z_S_a = 0.0E0;
    Double I_ESP_N7y3z_S_a = 0.0E0;
    Double I_ESP_N6y4z_S_a = 0.0E0;
    Double I_ESP_N5y5z_S_a = 0.0E0;
    Double I_ESP_N4y6z_S_a = 0.0E0;
    Double I_ESP_N3y7z_S_a = 0.0E0;
    Double I_ESP_N2y8z_S_a = 0.0E0;
    Double I_ESP_Ny9z_S_a = 0.0E0;
    Double I_ESP_N10z_S_a = 0.0E0;
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
    Double I_ESP_L8x_S = 0.0E0;
    Double I_ESP_L7xy_S = 0.0E0;
    Double I_ESP_L7xz_S = 0.0E0;
    Double I_ESP_L6x2y_S = 0.0E0;
    Double I_ESP_L6xyz_S = 0.0E0;
    Double I_ESP_L6x2z_S = 0.0E0;
    Double I_ESP_L5x3y_S = 0.0E0;
    Double I_ESP_L5x2yz_S = 0.0E0;
    Double I_ESP_L5xy2z_S = 0.0E0;
    Double I_ESP_L5x3z_S = 0.0E0;
    Double I_ESP_L4x4y_S = 0.0E0;
    Double I_ESP_L4x3yz_S = 0.0E0;
    Double I_ESP_L4x2y2z_S = 0.0E0;
    Double I_ESP_L4xy3z_S = 0.0E0;
    Double I_ESP_L4x4z_S = 0.0E0;
    Double I_ESP_L3x5y_S = 0.0E0;
    Double I_ESP_L3x4yz_S = 0.0E0;
    Double I_ESP_L3x3y2z_S = 0.0E0;
    Double I_ESP_L3x2y3z_S = 0.0E0;
    Double I_ESP_L3xy4z_S = 0.0E0;
    Double I_ESP_L3x5z_S = 0.0E0;
    Double I_ESP_L2x6y_S = 0.0E0;
    Double I_ESP_L2x5yz_S = 0.0E0;
    Double I_ESP_L2x4y2z_S = 0.0E0;
    Double I_ESP_L2x3y3z_S = 0.0E0;
    Double I_ESP_L2x2y4z_S = 0.0E0;
    Double I_ESP_L2xy5z_S = 0.0E0;
    Double I_ESP_L2x6z_S = 0.0E0;
    Double I_ESP_Lx7y_S = 0.0E0;
    Double I_ESP_Lx6yz_S = 0.0E0;
    Double I_ESP_Lx5y2z_S = 0.0E0;
    Double I_ESP_Lx4y3z_S = 0.0E0;
    Double I_ESP_Lx3y4z_S = 0.0E0;
    Double I_ESP_Lx2y5z_S = 0.0E0;
    Double I_ESP_Lxy6z_S = 0.0E0;
    Double I_ESP_Lx7z_S = 0.0E0;
    Double I_ESP_L8y_S = 0.0E0;
    Double I_ESP_L7yz_S = 0.0E0;
    Double I_ESP_L6y2z_S = 0.0E0;
    Double I_ESP_L5y3z_S = 0.0E0;
    Double I_ESP_L4y4z_S = 0.0E0;
    Double I_ESP_L3y5z_S = 0.0E0;
    Double I_ESP_L2y6z_S = 0.0E0;
    Double I_ESP_Ly7z_S = 0.0E0;
    Double I_ESP_L8z_S = 0.0E0;
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
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_D2x_S_M1_vrr = PAX*I_ESP_Px_S_M1_vrr-PRX*I_ESP_Px_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
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
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_D2x_S_vrr = PAX*I_ESP_Px_S_vrr-PRX*I_ESP_Px_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_D2y_S_vrr = PAY*I_ESP_Py_S_vrr-PRY*I_ESP_Py_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_D2z_S_vrr = PAZ*I_ESP_Pz_S_vrr-PRZ*I_ESP_Pz_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         ************************************************************/
        Double I_ESP_F3x_S_vrr = PAX*I_ESP_D2x_S_vrr-PRX*I_ESP_D2x_S_M1_vrr+2*oned2z*I_ESP_Px_S_vrr-2*oned2z*I_ESP_Px_S_M1_vrr;
        Double I_ESP_F2xy_S_vrr = PAY*I_ESP_D2x_S_vrr-PRY*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F2xz_S_vrr = PAZ*I_ESP_D2x_S_vrr-PRZ*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Fx2y_S_vrr = PAX*I_ESP_D2y_S_vrr-PRX*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fx2z_S_vrr = PAX*I_ESP_D2z_S_vrr-PRX*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3y_S_vrr = PAY*I_ESP_D2y_S_vrr-PRY*I_ESP_D2y_S_M1_vrr+2*oned2z*I_ESP_Py_S_vrr-2*oned2z*I_ESP_Py_S_M1_vrr;
        Double I_ESP_F2yz_S_vrr = PAZ*I_ESP_D2y_S_vrr-PRZ*I_ESP_D2y_S_M1_vrr;
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
         * shell quartet name: SQ_ESP_N_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_N_S_a_coefs = alpha;
        I_ESP_N10x_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N10x_S_vrr;
        I_ESP_N9xy_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N9xy_S_vrr;
        I_ESP_N9xz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N9xz_S_vrr;
        I_ESP_N8x2y_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N8x2y_S_vrr;
        I_ESP_N8xyz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N8xyz_S_vrr;
        I_ESP_N8x2z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N8x2z_S_vrr;
        I_ESP_N7x3y_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N7x3y_S_vrr;
        I_ESP_N7x2yz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N7x2yz_S_vrr;
        I_ESP_N7xy2z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N7xy2z_S_vrr;
        I_ESP_N7x3z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N7x3z_S_vrr;
        I_ESP_N6x4y_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N6x4y_S_vrr;
        I_ESP_N6x3yz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N6x3yz_S_vrr;
        I_ESP_N6x2y2z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N6x2y2z_S_vrr;
        I_ESP_N6xy3z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N6xy3z_S_vrr;
        I_ESP_N6x4z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N6x4z_S_vrr;
        I_ESP_N5x5y_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N5x5y_S_vrr;
        I_ESP_N5x4yz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N5x4yz_S_vrr;
        I_ESP_N5x3y2z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N5x3y2z_S_vrr;
        I_ESP_N5x2y3z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N5x2y3z_S_vrr;
        I_ESP_N5xy4z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N5xy4z_S_vrr;
        I_ESP_N5x5z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N5x5z_S_vrr;
        I_ESP_N4x6y_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N4x6y_S_vrr;
        I_ESP_N4x5yz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N4x5yz_S_vrr;
        I_ESP_N4x4y2z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N4x4y2z_S_vrr;
        I_ESP_N4x3y3z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N4x3y3z_S_vrr;
        I_ESP_N4x2y4z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N4x2y4z_S_vrr;
        I_ESP_N4xy5z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N4xy5z_S_vrr;
        I_ESP_N4x6z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N4x6z_S_vrr;
        I_ESP_N3x7y_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N3x7y_S_vrr;
        I_ESP_N3x6yz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N3x6yz_S_vrr;
        I_ESP_N3x5y2z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N3x5y2z_S_vrr;
        I_ESP_N3x4y3z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N3x4y3z_S_vrr;
        I_ESP_N3x3y4z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N3x3y4z_S_vrr;
        I_ESP_N3x2y5z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N3x2y5z_S_vrr;
        I_ESP_N3xy6z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N3xy6z_S_vrr;
        I_ESP_N3x7z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N3x7z_S_vrr;
        I_ESP_N2x8y_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2x8y_S_vrr;
        I_ESP_N2x7yz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2x7yz_S_vrr;
        I_ESP_N2x6y2z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2x6y2z_S_vrr;
        I_ESP_N2x5y3z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2x5y3z_S_vrr;
        I_ESP_N2x4y4z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2x4y4z_S_vrr;
        I_ESP_N2x3y5z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2x3y5z_S_vrr;
        I_ESP_N2x2y6z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2x2y6z_S_vrr;
        I_ESP_N2xy7z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2xy7z_S_vrr;
        I_ESP_N2x8z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2x8z_S_vrr;
        I_ESP_Nx9y_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nx9y_S_vrr;
        I_ESP_Nx8yz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nx8yz_S_vrr;
        I_ESP_Nx7y2z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nx7y2z_S_vrr;
        I_ESP_Nx6y3z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nx6y3z_S_vrr;
        I_ESP_Nx5y4z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nx5y4z_S_vrr;
        I_ESP_Nx4y5z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nx4y5z_S_vrr;
        I_ESP_Nx3y6z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nx3y6z_S_vrr;
        I_ESP_Nx2y7z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nx2y7z_S_vrr;
        I_ESP_Nxy8z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nxy8z_S_vrr;
        I_ESP_Nx9z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Nx9z_S_vrr;
        I_ESP_N10y_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N10y_S_vrr;
        I_ESP_N9yz_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N9yz_S_vrr;
        I_ESP_N8y2z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N8y2z_S_vrr;
        I_ESP_N7y3z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N7y3z_S_vrr;
        I_ESP_N6y4z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N6y4z_S_vrr;
        I_ESP_N5y5z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N5y5z_S_vrr;
        I_ESP_N4y6z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N4y6z_S_vrr;
        I_ESP_N3y7z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N3y7z_S_vrr;
        I_ESP_N2y8z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N2y8z_S_vrr;
        I_ESP_Ny9z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_Ny9z_S_vrr;
        I_ESP_N10z_S_a += SQ_ESP_N_S_a_coefs*I_ESP_N10z_S_vrr;

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
         * shell quartet name: SQ_ESP_L_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_L8x_S += I_ESP_L8x_S_vrr;
        I_ESP_L7xy_S += I_ESP_L7xy_S_vrr;
        I_ESP_L7xz_S += I_ESP_L7xz_S_vrr;
        I_ESP_L6x2y_S += I_ESP_L6x2y_S_vrr;
        I_ESP_L6xyz_S += I_ESP_L6xyz_S_vrr;
        I_ESP_L6x2z_S += I_ESP_L6x2z_S_vrr;
        I_ESP_L5x3y_S += I_ESP_L5x3y_S_vrr;
        I_ESP_L5x2yz_S += I_ESP_L5x2yz_S_vrr;
        I_ESP_L5xy2z_S += I_ESP_L5xy2z_S_vrr;
        I_ESP_L5x3z_S += I_ESP_L5x3z_S_vrr;
        I_ESP_L4x4y_S += I_ESP_L4x4y_S_vrr;
        I_ESP_L4x3yz_S += I_ESP_L4x3yz_S_vrr;
        I_ESP_L4x2y2z_S += I_ESP_L4x2y2z_S_vrr;
        I_ESP_L4xy3z_S += I_ESP_L4xy3z_S_vrr;
        I_ESP_L4x4z_S += I_ESP_L4x4z_S_vrr;
        I_ESP_L3x5y_S += I_ESP_L3x5y_S_vrr;
        I_ESP_L3x4yz_S += I_ESP_L3x4yz_S_vrr;
        I_ESP_L3x3y2z_S += I_ESP_L3x3y2z_S_vrr;
        I_ESP_L3x2y3z_S += I_ESP_L3x2y3z_S_vrr;
        I_ESP_L3xy4z_S += I_ESP_L3xy4z_S_vrr;
        I_ESP_L3x5z_S += I_ESP_L3x5z_S_vrr;
        I_ESP_L2x6y_S += I_ESP_L2x6y_S_vrr;
        I_ESP_L2x5yz_S += I_ESP_L2x5yz_S_vrr;
        I_ESP_L2x4y2z_S += I_ESP_L2x4y2z_S_vrr;
        I_ESP_L2x3y3z_S += I_ESP_L2x3y3z_S_vrr;
        I_ESP_L2x2y4z_S += I_ESP_L2x2y4z_S_vrr;
        I_ESP_L2xy5z_S += I_ESP_L2xy5z_S_vrr;
        I_ESP_L2x6z_S += I_ESP_L2x6z_S_vrr;
        I_ESP_Lx7y_S += I_ESP_Lx7y_S_vrr;
        I_ESP_Lx6yz_S += I_ESP_Lx6yz_S_vrr;
        I_ESP_Lx5y2z_S += I_ESP_Lx5y2z_S_vrr;
        I_ESP_Lx4y3z_S += I_ESP_Lx4y3z_S_vrr;
        I_ESP_Lx3y4z_S += I_ESP_Lx3y4z_S_vrr;
        I_ESP_Lx2y5z_S += I_ESP_Lx2y5z_S_vrr;
        I_ESP_Lxy6z_S += I_ESP_Lxy6z_S_vrr;
        I_ESP_Lx7z_S += I_ESP_Lx7z_S_vrr;
        I_ESP_L8y_S += I_ESP_L8y_S_vrr;
        I_ESP_L7yz_S += I_ESP_L7yz_S_vrr;
        I_ESP_L6y2z_S += I_ESP_L6y2z_S_vrr;
        I_ESP_L5y3z_S += I_ESP_L5y3z_S_vrr;
        I_ESP_L4y4z_S += I_ESP_L4y4z_S_vrr;
        I_ESP_L3y5z_S += I_ESP_L3y5z_S_vrr;
        I_ESP_L2y6z_S += I_ESP_L2y6z_S_vrr;
        I_ESP_Ly7z_S += I_ESP_Ly7z_S_vrr;
        I_ESP_L8z_S += I_ESP_L8z_S_vrr;

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
     * shell quartet name: SQ_ESP_H_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
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
    Double I_ESP_H5x_Py = I_ESP_I5xy_S+ABY*I_ESP_H5x_S;
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
    Double I_ESP_H5x_Pz = I_ESP_I5xz_S+ABZ*I_ESP_H5x_S;
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
    Double I_ESP_H5y_Pz = I_ESP_I5yz_S+ABZ*I_ESP_H5y_S;
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
     * shell quartet name: SQ_ESP_I_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 3 integrals are omitted 
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
    Double I_ESP_I6y_Px = I_ESP_Kx6y_S+ABX*I_ESP_I6y_S;
    Double I_ESP_I5yz_Px = I_ESP_Kx5yz_S+ABX*I_ESP_I5yz_S;
    Double I_ESP_I4y2z_Px = I_ESP_Kx4y2z_S+ABX*I_ESP_I4y2z_S;
    Double I_ESP_I3y3z_Px = I_ESP_Kx3y3z_S+ABX*I_ESP_I3y3z_S;
    Double I_ESP_I2y4z_Px = I_ESP_Kx2y4z_S+ABX*I_ESP_I2y4z_S;
    Double I_ESP_Iy5z_Px = I_ESP_Kxy5z_S+ABX*I_ESP_Iy5z_S;
    Double I_ESP_I6z_Px = I_ESP_Kx6z_S+ABX*I_ESP_I6z_S;
    Double I_ESP_I5xy_Py = I_ESP_K5x2y_S+ABY*I_ESP_I5xy_S;
    Double I_ESP_I5xz_Py = I_ESP_K5xyz_S+ABY*I_ESP_I5xz_S;
    Double I_ESP_I4x2y_Py = I_ESP_K4x3y_S+ABY*I_ESP_I4x2y_S;
    Double I_ESP_I4xyz_Py = I_ESP_K4x2yz_S+ABY*I_ESP_I4xyz_S;
    Double I_ESP_I4x2z_Py = I_ESP_K4xy2z_S+ABY*I_ESP_I4x2z_S;
    Double I_ESP_I3x3y_Py = I_ESP_K3x4y_S+ABY*I_ESP_I3x3y_S;
    Double I_ESP_I3x2yz_Py = I_ESP_K3x3yz_S+ABY*I_ESP_I3x2yz_S;
    Double I_ESP_I3xy2z_Py = I_ESP_K3x2y2z_S+ABY*I_ESP_I3xy2z_S;
    Double I_ESP_I3x3z_Py = I_ESP_K3xy3z_S+ABY*I_ESP_I3x3z_S;
    Double I_ESP_I2x4y_Py = I_ESP_K2x5y_S+ABY*I_ESP_I2x4y_S;
    Double I_ESP_I2x3yz_Py = I_ESP_K2x4yz_S+ABY*I_ESP_I2x3yz_S;
    Double I_ESP_I2x2y2z_Py = I_ESP_K2x3y2z_S+ABY*I_ESP_I2x2y2z_S;
    Double I_ESP_I2xy3z_Py = I_ESP_K2x2y3z_S+ABY*I_ESP_I2xy3z_S;
    Double I_ESP_I2x4z_Py = I_ESP_K2xy4z_S+ABY*I_ESP_I2x4z_S;
    Double I_ESP_Ix5y_Py = I_ESP_Kx6y_S+ABY*I_ESP_Ix5y_S;
    Double I_ESP_Ix4yz_Py = I_ESP_Kx5yz_S+ABY*I_ESP_Ix4yz_S;
    Double I_ESP_Ix3y2z_Py = I_ESP_Kx4y2z_S+ABY*I_ESP_Ix3y2z_S;
    Double I_ESP_Ix2y3z_Py = I_ESP_Kx3y3z_S+ABY*I_ESP_Ix2y3z_S;
    Double I_ESP_Ixy4z_Py = I_ESP_Kx2y4z_S+ABY*I_ESP_Ixy4z_S;
    Double I_ESP_Ix5z_Py = I_ESP_Kxy5z_S+ABY*I_ESP_Ix5z_S;
    Double I_ESP_I6y_Py = I_ESP_K7y_S+ABY*I_ESP_I6y_S;
    Double I_ESP_I5yz_Py = I_ESP_K6yz_S+ABY*I_ESP_I5yz_S;
    Double I_ESP_I4y2z_Py = I_ESP_K5y2z_S+ABY*I_ESP_I4y2z_S;
    Double I_ESP_I3y3z_Py = I_ESP_K4y3z_S+ABY*I_ESP_I3y3z_S;
    Double I_ESP_I2y4z_Py = I_ESP_K3y4z_S+ABY*I_ESP_I2y4z_S;
    Double I_ESP_Iy5z_Py = I_ESP_K2y5z_S+ABY*I_ESP_Iy5z_S;
    Double I_ESP_I6z_Py = I_ESP_Ky6z_S+ABY*I_ESP_I6z_S;
    Double I_ESP_I5xy_Pz = I_ESP_K5xyz_S+ABZ*I_ESP_I5xy_S;
    Double I_ESP_I5xz_Pz = I_ESP_K5x2z_S+ABZ*I_ESP_I5xz_S;
    Double I_ESP_I4x2y_Pz = I_ESP_K4x2yz_S+ABZ*I_ESP_I4x2y_S;
    Double I_ESP_I4xyz_Pz = I_ESP_K4xy2z_S+ABZ*I_ESP_I4xyz_S;
    Double I_ESP_I4x2z_Pz = I_ESP_K4x3z_S+ABZ*I_ESP_I4x2z_S;
    Double I_ESP_I3x3y_Pz = I_ESP_K3x3yz_S+ABZ*I_ESP_I3x3y_S;
    Double I_ESP_I3x2yz_Pz = I_ESP_K3x2y2z_S+ABZ*I_ESP_I3x2yz_S;
    Double I_ESP_I3xy2z_Pz = I_ESP_K3xy3z_S+ABZ*I_ESP_I3xy2z_S;
    Double I_ESP_I3x3z_Pz = I_ESP_K3x4z_S+ABZ*I_ESP_I3x3z_S;
    Double I_ESP_I2x4y_Pz = I_ESP_K2x4yz_S+ABZ*I_ESP_I2x4y_S;
    Double I_ESP_I2x3yz_Pz = I_ESP_K2x3y2z_S+ABZ*I_ESP_I2x3yz_S;
    Double I_ESP_I2x2y2z_Pz = I_ESP_K2x2y3z_S+ABZ*I_ESP_I2x2y2z_S;
    Double I_ESP_I2xy3z_Pz = I_ESP_K2xy4z_S+ABZ*I_ESP_I2xy3z_S;
    Double I_ESP_I2x4z_Pz = I_ESP_K2x5z_S+ABZ*I_ESP_I2x4z_S;
    Double I_ESP_Ix5y_Pz = I_ESP_Kx5yz_S+ABZ*I_ESP_Ix5y_S;
    Double I_ESP_Ix4yz_Pz = I_ESP_Kx4y2z_S+ABZ*I_ESP_Ix4yz_S;
    Double I_ESP_Ix3y2z_Pz = I_ESP_Kx3y3z_S+ABZ*I_ESP_Ix3y2z_S;
    Double I_ESP_Ix2y3z_Pz = I_ESP_Kx2y4z_S+ABZ*I_ESP_Ix2y3z_S;
    Double I_ESP_Ixy4z_Pz = I_ESP_Kxy5z_S+ABZ*I_ESP_Ixy4z_S;
    Double I_ESP_Ix5z_Pz = I_ESP_Kx6z_S+ABZ*I_ESP_Ix5z_S;
    Double I_ESP_I5yz_Pz = I_ESP_K5y2z_S+ABZ*I_ESP_I5yz_S;
    Double I_ESP_I4y2z_Pz = I_ESP_K4y3z_S+ABZ*I_ESP_I4y2z_S;
    Double I_ESP_I3y3z_Pz = I_ESP_K3y4z_S+ABZ*I_ESP_I3y3z_S;
    Double I_ESP_I2y4z_Pz = I_ESP_K2y5z_S+ABZ*I_ESP_I2y4z_S;
    Double I_ESP_Iy5z_Pz = I_ESP_Ky6z_S+ABZ*I_ESP_Iy5z_S;
    Double I_ESP_I6z_Pz = I_ESP_K7z_S+ABZ*I_ESP_I6z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_H_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 63 integrals are omitted 
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
    Double I_ESP_H5x_D2y = I_ESP_I5xy_Py+ABY*I_ESP_H5x_Py;
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
    Double I_ESP_H5x_D2z = I_ESP_I5xz_Pz+ABZ*I_ESP_H5x_Pz;
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
    Double I_ESP_H5y_D2z = I_ESP_I5yz_Pz+ABZ*I_ESP_H5y_Pz;
    Double I_ESP_H4yz_D2z = I_ESP_I4y2z_Pz+ABZ*I_ESP_H4yz_Pz;
    Double I_ESP_H3y2z_D2z = I_ESP_I3y3z_Pz+ABZ*I_ESP_H3y2z_Pz;
    Double I_ESP_H2y3z_D2z = I_ESP_I2y4z_Pz+ABZ*I_ESP_H2y3z_Pz;
    Double I_ESP_Hy4z_D2z = I_ESP_Iy5z_Pz+ABZ*I_ESP_Hy4z_Pz;
    Double I_ESP_H5z_D2z = I_ESP_I6z_Pz+ABZ*I_ESP_H5z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_G_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 30 integrals are omitted 
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
    Double I_ESP_G4x_F2xy = I_ESP_H4xy_D2x+ABY*I_ESP_G4x_D2x;
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
    Double I_ESP_G4x_F2xz = I_ESP_H4xz_D2x+ABZ*I_ESP_G4x_D2x;
    Double I_ESP_G3xy_F2xz = I_ESP_H3xyz_D2x+ABZ*I_ESP_G3xy_D2x;
    Double I_ESP_G3xz_F2xz = I_ESP_H3x2z_D2x+ABZ*I_ESP_G3xz_D2x;
    Double I_ESP_G2x2y_F2xz = I_ESP_H2x2yz_D2x+ABZ*I_ESP_G2x2y_D2x;
    Double I_ESP_G2xyz_F2xz = I_ESP_H2xy2z_D2x+ABZ*I_ESP_G2xyz_D2x;
    Double I_ESP_G2x2z_F2xz = I_ESP_H2x3z_D2x+ABZ*I_ESP_G2x2z_D2x;
    Double I_ESP_Gx3y_F2xz = I_ESP_Hx3yz_D2x+ABZ*I_ESP_Gx3y_D2x;
    Double I_ESP_Gx2yz_F2xz = I_ESP_Hx2y2z_D2x+ABZ*I_ESP_Gx2yz_D2x;
    Double I_ESP_Gxy2z_F2xz = I_ESP_Hxy3z_D2x+ABZ*I_ESP_Gxy2z_D2x;
    Double I_ESP_Gx3z_F2xz = I_ESP_Hx4z_D2x+ABZ*I_ESP_Gx3z_D2x;
    Double I_ESP_G4y_F2xz = I_ESP_H4yz_D2x+ABZ*I_ESP_G4y_D2x;
    Double I_ESP_G3yz_F2xz = I_ESP_H3y2z_D2x+ABZ*I_ESP_G3yz_D2x;
    Double I_ESP_G2y2z_F2xz = I_ESP_H2y3z_D2x+ABZ*I_ESP_G2y2z_D2x;
    Double I_ESP_Gy3z_F2xz = I_ESP_Hy4z_D2x+ABZ*I_ESP_Gy3z_D2x;
    Double I_ESP_G4z_F2xz = I_ESP_H5z_D2x+ABZ*I_ESP_G4z_D2x;
    Double I_ESP_G4x_Fx2y = I_ESP_H5x_D2y+ABX*I_ESP_G4x_D2y;
    Double I_ESP_G3xy_Fx2y = I_ESP_H4xy_D2y+ABX*I_ESP_G3xy_D2y;
    Double I_ESP_G3xz_Fx2y = I_ESP_H4xz_D2y+ABX*I_ESP_G3xz_D2y;
    Double I_ESP_G2x2y_Fx2y = I_ESP_H3x2y_D2y+ABX*I_ESP_G2x2y_D2y;
    Double I_ESP_G2xyz_Fx2y = I_ESP_H3xyz_D2y+ABX*I_ESP_G2xyz_D2y;
    Double I_ESP_G2x2z_Fx2y = I_ESP_H3x2z_D2y+ABX*I_ESP_G2x2z_D2y;
    Double I_ESP_Gx3y_Fx2y = I_ESP_H2x3y_D2y+ABX*I_ESP_Gx3y_D2y;
    Double I_ESP_Gx2yz_Fx2y = I_ESP_H2x2yz_D2y+ABX*I_ESP_Gx2yz_D2y;
    Double I_ESP_Gxy2z_Fx2y = I_ESP_H2xy2z_D2y+ABX*I_ESP_Gxy2z_D2y;
    Double I_ESP_Gx3z_Fx2y = I_ESP_H2x3z_D2y+ABX*I_ESP_Gx3z_D2y;
    Double I_ESP_G4y_Fx2y = I_ESP_Hx4y_D2y+ABX*I_ESP_G4y_D2y;
    Double I_ESP_G3yz_Fx2y = I_ESP_Hx3yz_D2y+ABX*I_ESP_G3yz_D2y;
    Double I_ESP_G2y2z_Fx2y = I_ESP_Hx2y2z_D2y+ABX*I_ESP_G2y2z_D2y;
    Double I_ESP_Gy3z_Fx2y = I_ESP_Hxy3z_D2y+ABX*I_ESP_Gy3z_D2y;
    Double I_ESP_G4z_Fx2y = I_ESP_Hx4z_D2y+ABX*I_ESP_G4z_D2y;
    Double I_ESP_G4x_Fx2z = I_ESP_H5x_D2z+ABX*I_ESP_G4x_D2z;
    Double I_ESP_G3xy_Fx2z = I_ESP_H4xy_D2z+ABX*I_ESP_G3xy_D2z;
    Double I_ESP_G3xz_Fx2z = I_ESP_H4xz_D2z+ABX*I_ESP_G3xz_D2z;
    Double I_ESP_G2x2y_Fx2z = I_ESP_H3x2y_D2z+ABX*I_ESP_G2x2y_D2z;
    Double I_ESP_G2xyz_Fx2z = I_ESP_H3xyz_D2z+ABX*I_ESP_G2xyz_D2z;
    Double I_ESP_G2x2z_Fx2z = I_ESP_H3x2z_D2z+ABX*I_ESP_G2x2z_D2z;
    Double I_ESP_Gx3y_Fx2z = I_ESP_H2x3y_D2z+ABX*I_ESP_Gx3y_D2z;
    Double I_ESP_Gx2yz_Fx2z = I_ESP_H2x2yz_D2z+ABX*I_ESP_Gx2yz_D2z;
    Double I_ESP_Gxy2z_Fx2z = I_ESP_H2xy2z_D2z+ABX*I_ESP_Gxy2z_D2z;
    Double I_ESP_Gx3z_Fx2z = I_ESP_H2x3z_D2z+ABX*I_ESP_Gx3z_D2z;
    Double I_ESP_G4y_Fx2z = I_ESP_Hx4y_D2z+ABX*I_ESP_G4y_D2z;
    Double I_ESP_G3yz_Fx2z = I_ESP_Hx3yz_D2z+ABX*I_ESP_G3yz_D2z;
    Double I_ESP_G2y2z_Fx2z = I_ESP_Hx2y2z_D2z+ABX*I_ESP_G2y2z_D2z;
    Double I_ESP_Gy3z_Fx2z = I_ESP_Hxy3z_D2z+ABX*I_ESP_Gy3z_D2z;
    Double I_ESP_G4z_Fx2z = I_ESP_Hx4z_D2z+ABX*I_ESP_G4z_D2z;
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
    Double I_ESP_G4x_F2yz = I_ESP_H4xz_D2y+ABZ*I_ESP_G4x_D2y;
    Double I_ESP_G3xy_F2yz = I_ESP_H3xyz_D2y+ABZ*I_ESP_G3xy_D2y;
    Double I_ESP_G3xz_F2yz = I_ESP_H3x2z_D2y+ABZ*I_ESP_G3xz_D2y;
    Double I_ESP_G2x2y_F2yz = I_ESP_H2x2yz_D2y+ABZ*I_ESP_G2x2y_D2y;
    Double I_ESP_G2xyz_F2yz = I_ESP_H2xy2z_D2y+ABZ*I_ESP_G2xyz_D2y;
    Double I_ESP_G2x2z_F2yz = I_ESP_H2x3z_D2y+ABZ*I_ESP_G2x2z_D2y;
    Double I_ESP_Gx3y_F2yz = I_ESP_Hx3yz_D2y+ABZ*I_ESP_Gx3y_D2y;
    Double I_ESP_Gx2yz_F2yz = I_ESP_Hx2y2z_D2y+ABZ*I_ESP_Gx2yz_D2y;
    Double I_ESP_Gxy2z_F2yz = I_ESP_Hxy3z_D2y+ABZ*I_ESP_Gxy2z_D2y;
    Double I_ESP_Gx3z_F2yz = I_ESP_Hx4z_D2y+ABZ*I_ESP_Gx3z_D2y;
    Double I_ESP_G4y_F2yz = I_ESP_H4yz_D2y+ABZ*I_ESP_G4y_D2y;
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
     * shell quartet name: SQ_ESP_K_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 27 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_S
     * RHS shell quartet name: SQ_ESP_K_S
     ************************************************************/
    Double I_ESP_K7x_Px = I_ESP_L8x_S+ABX*I_ESP_K7x_S;
    Double I_ESP_K6xy_Px = I_ESP_L7xy_S+ABX*I_ESP_K6xy_S;
    Double I_ESP_K6xz_Px = I_ESP_L7xz_S+ABX*I_ESP_K6xz_S;
    Double I_ESP_K5x2y_Px = I_ESP_L6x2y_S+ABX*I_ESP_K5x2y_S;
    Double I_ESP_K5xyz_Px = I_ESP_L6xyz_S+ABX*I_ESP_K5xyz_S;
    Double I_ESP_K5x2z_Px = I_ESP_L6x2z_S+ABX*I_ESP_K5x2z_S;
    Double I_ESP_K4x3y_Px = I_ESP_L5x3y_S+ABX*I_ESP_K4x3y_S;
    Double I_ESP_K4x2yz_Px = I_ESP_L5x2yz_S+ABX*I_ESP_K4x2yz_S;
    Double I_ESP_K4xy2z_Px = I_ESP_L5xy2z_S+ABX*I_ESP_K4xy2z_S;
    Double I_ESP_K4x3z_Px = I_ESP_L5x3z_S+ABX*I_ESP_K4x3z_S;
    Double I_ESP_K3x4y_Px = I_ESP_L4x4y_S+ABX*I_ESP_K3x4y_S;
    Double I_ESP_K3x3yz_Px = I_ESP_L4x3yz_S+ABX*I_ESP_K3x3yz_S;
    Double I_ESP_K3x2y2z_Px = I_ESP_L4x2y2z_S+ABX*I_ESP_K3x2y2z_S;
    Double I_ESP_K3xy3z_Px = I_ESP_L4xy3z_S+ABX*I_ESP_K3xy3z_S;
    Double I_ESP_K3x4z_Px = I_ESP_L4x4z_S+ABX*I_ESP_K3x4z_S;
    Double I_ESP_K2x5y_Px = I_ESP_L3x5y_S+ABX*I_ESP_K2x5y_S;
    Double I_ESP_K2x4yz_Px = I_ESP_L3x4yz_S+ABX*I_ESP_K2x4yz_S;
    Double I_ESP_K2x3y2z_Px = I_ESP_L3x3y2z_S+ABX*I_ESP_K2x3y2z_S;
    Double I_ESP_K2x2y3z_Px = I_ESP_L3x2y3z_S+ABX*I_ESP_K2x2y3z_S;
    Double I_ESP_K2xy4z_Px = I_ESP_L3xy4z_S+ABX*I_ESP_K2xy4z_S;
    Double I_ESP_K2x5z_Px = I_ESP_L3x5z_S+ABX*I_ESP_K2x5z_S;
    Double I_ESP_Kx6y_Px = I_ESP_L2x6y_S+ABX*I_ESP_Kx6y_S;
    Double I_ESP_Kx5yz_Px = I_ESP_L2x5yz_S+ABX*I_ESP_Kx5yz_S;
    Double I_ESP_Kx4y2z_Px = I_ESP_L2x4y2z_S+ABX*I_ESP_Kx4y2z_S;
    Double I_ESP_Kx3y3z_Px = I_ESP_L2x3y3z_S+ABX*I_ESP_Kx3y3z_S;
    Double I_ESP_Kx2y4z_Px = I_ESP_L2x2y4z_S+ABX*I_ESP_Kx2y4z_S;
    Double I_ESP_Kxy5z_Px = I_ESP_L2xy5z_S+ABX*I_ESP_Kxy5z_S;
    Double I_ESP_Kx6z_Px = I_ESP_L2x6z_S+ABX*I_ESP_Kx6z_S;
    Double I_ESP_K5x2y_Py = I_ESP_L5x3y_S+ABY*I_ESP_K5x2y_S;
    Double I_ESP_K5xyz_Py = I_ESP_L5x2yz_S+ABY*I_ESP_K5xyz_S;
    Double I_ESP_K4x3y_Py = I_ESP_L4x4y_S+ABY*I_ESP_K4x3y_S;
    Double I_ESP_K4x2yz_Py = I_ESP_L4x3yz_S+ABY*I_ESP_K4x2yz_S;
    Double I_ESP_K4xy2z_Py = I_ESP_L4x2y2z_S+ABY*I_ESP_K4xy2z_S;
    Double I_ESP_K3x4y_Py = I_ESP_L3x5y_S+ABY*I_ESP_K3x4y_S;
    Double I_ESP_K3x3yz_Py = I_ESP_L3x4yz_S+ABY*I_ESP_K3x3yz_S;
    Double I_ESP_K3x2y2z_Py = I_ESP_L3x3y2z_S+ABY*I_ESP_K3x2y2z_S;
    Double I_ESP_K3xy3z_Py = I_ESP_L3x2y3z_S+ABY*I_ESP_K3xy3z_S;
    Double I_ESP_K2x5y_Py = I_ESP_L2x6y_S+ABY*I_ESP_K2x5y_S;
    Double I_ESP_K2x4yz_Py = I_ESP_L2x5yz_S+ABY*I_ESP_K2x4yz_S;
    Double I_ESP_K2x3y2z_Py = I_ESP_L2x4y2z_S+ABY*I_ESP_K2x3y2z_S;
    Double I_ESP_K2x2y3z_Py = I_ESP_L2x3y3z_S+ABY*I_ESP_K2x2y3z_S;
    Double I_ESP_K2xy4z_Py = I_ESP_L2x2y4z_S+ABY*I_ESP_K2xy4z_S;
    Double I_ESP_Kx6y_Py = I_ESP_Lx7y_S+ABY*I_ESP_Kx6y_S;
    Double I_ESP_Kx5yz_Py = I_ESP_Lx6yz_S+ABY*I_ESP_Kx5yz_S;
    Double I_ESP_Kx4y2z_Py = I_ESP_Lx5y2z_S+ABY*I_ESP_Kx4y2z_S;
    Double I_ESP_Kx3y3z_Py = I_ESP_Lx4y3z_S+ABY*I_ESP_Kx3y3z_S;
    Double I_ESP_Kx2y4z_Py = I_ESP_Lx3y4z_S+ABY*I_ESP_Kx2y4z_S;
    Double I_ESP_Kxy5z_Py = I_ESP_Lx2y5z_S+ABY*I_ESP_Kxy5z_S;
    Double I_ESP_K7y_Py = I_ESP_L8y_S+ABY*I_ESP_K7y_S;
    Double I_ESP_K6yz_Py = I_ESP_L7yz_S+ABY*I_ESP_K6yz_S;
    Double I_ESP_K5y2z_Py = I_ESP_L6y2z_S+ABY*I_ESP_K5y2z_S;
    Double I_ESP_K4y3z_Py = I_ESP_L5y3z_S+ABY*I_ESP_K4y3z_S;
    Double I_ESP_K3y4z_Py = I_ESP_L4y4z_S+ABY*I_ESP_K3y4z_S;
    Double I_ESP_K2y5z_Py = I_ESP_L3y5z_S+ABY*I_ESP_K2y5z_S;
    Double I_ESP_Ky6z_Py = I_ESP_L2y6z_S+ABY*I_ESP_Ky6z_S;
    Double I_ESP_K5xyz_Pz = I_ESP_L5xy2z_S+ABZ*I_ESP_K5xyz_S;
    Double I_ESP_K5x2z_Pz = I_ESP_L5x3z_S+ABZ*I_ESP_K5x2z_S;
    Double I_ESP_K4x2yz_Pz = I_ESP_L4x2y2z_S+ABZ*I_ESP_K4x2yz_S;
    Double I_ESP_K4xy2z_Pz = I_ESP_L4xy3z_S+ABZ*I_ESP_K4xy2z_S;
    Double I_ESP_K4x3z_Pz = I_ESP_L4x4z_S+ABZ*I_ESP_K4x3z_S;
    Double I_ESP_K3x3yz_Pz = I_ESP_L3x3y2z_S+ABZ*I_ESP_K3x3yz_S;
    Double I_ESP_K3x2y2z_Pz = I_ESP_L3x2y3z_S+ABZ*I_ESP_K3x2y2z_S;
    Double I_ESP_K3xy3z_Pz = I_ESP_L3xy4z_S+ABZ*I_ESP_K3xy3z_S;
    Double I_ESP_K3x4z_Pz = I_ESP_L3x5z_S+ABZ*I_ESP_K3x4z_S;
    Double I_ESP_K2x4yz_Pz = I_ESP_L2x4y2z_S+ABZ*I_ESP_K2x4yz_S;
    Double I_ESP_K2x3y2z_Pz = I_ESP_L2x3y3z_S+ABZ*I_ESP_K2x3y2z_S;
    Double I_ESP_K2x2y3z_Pz = I_ESP_L2x2y4z_S+ABZ*I_ESP_K2x2y3z_S;
    Double I_ESP_K2xy4z_Pz = I_ESP_L2xy5z_S+ABZ*I_ESP_K2xy4z_S;
    Double I_ESP_K2x5z_Pz = I_ESP_L2x6z_S+ABZ*I_ESP_K2x5z_S;
    Double I_ESP_Kx5yz_Pz = I_ESP_Lx5y2z_S+ABZ*I_ESP_Kx5yz_S;
    Double I_ESP_Kx4y2z_Pz = I_ESP_Lx4y3z_S+ABZ*I_ESP_Kx4y2z_S;
    Double I_ESP_Kx3y3z_Pz = I_ESP_Lx3y4z_S+ABZ*I_ESP_Kx3y3z_S;
    Double I_ESP_Kx2y4z_Pz = I_ESP_Lx2y5z_S+ABZ*I_ESP_Kx2y4z_S;
    Double I_ESP_Kxy5z_Pz = I_ESP_Lxy6z_S+ABZ*I_ESP_Kxy5z_S;
    Double I_ESP_Kx6z_Pz = I_ESP_Lx7z_S+ABZ*I_ESP_Kx6z_S;
    Double I_ESP_K5y2z_Pz = I_ESP_L5y3z_S+ABZ*I_ESP_K5y2z_S;
    Double I_ESP_K4y3z_Pz = I_ESP_L4y4z_S+ABZ*I_ESP_K4y3z_S;
    Double I_ESP_K3y4z_Pz = I_ESP_L3y5z_S+ABZ*I_ESP_K3y4z_S;
    Double I_ESP_K2y5z_Pz = I_ESP_L2y6z_S+ABZ*I_ESP_K2y5z_S;
    Double I_ESP_Ky6z_Pz = I_ESP_Ly7z_S+ABZ*I_ESP_Ky6z_S;
    Double I_ESP_K7z_Pz = I_ESP_L8z_S+ABZ*I_ESP_K7z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_I_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 87 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_P
     * RHS shell quartet name: SQ_ESP_I_P
     ************************************************************/
    Double I_ESP_I6x_D2x = I_ESP_K7x_Px+ABX*I_ESP_I6x_Px;
    Double I_ESP_I5xy_D2x = I_ESP_K6xy_Px+ABX*I_ESP_I5xy_Px;
    Double I_ESP_I5xz_D2x = I_ESP_K6xz_Px+ABX*I_ESP_I5xz_Px;
    Double I_ESP_I4x2y_D2x = I_ESP_K5x2y_Px+ABX*I_ESP_I4x2y_Px;
    Double I_ESP_I4xyz_D2x = I_ESP_K5xyz_Px+ABX*I_ESP_I4xyz_Px;
    Double I_ESP_I4x2z_D2x = I_ESP_K5x2z_Px+ABX*I_ESP_I4x2z_Px;
    Double I_ESP_I3x3y_D2x = I_ESP_K4x3y_Px+ABX*I_ESP_I3x3y_Px;
    Double I_ESP_I3x2yz_D2x = I_ESP_K4x2yz_Px+ABX*I_ESP_I3x2yz_Px;
    Double I_ESP_I3xy2z_D2x = I_ESP_K4xy2z_Px+ABX*I_ESP_I3xy2z_Px;
    Double I_ESP_I3x3z_D2x = I_ESP_K4x3z_Px+ABX*I_ESP_I3x3z_Px;
    Double I_ESP_I2x4y_D2x = I_ESP_K3x4y_Px+ABX*I_ESP_I2x4y_Px;
    Double I_ESP_I2x3yz_D2x = I_ESP_K3x3yz_Px+ABX*I_ESP_I2x3yz_Px;
    Double I_ESP_I2x2y2z_D2x = I_ESP_K3x2y2z_Px+ABX*I_ESP_I2x2y2z_Px;
    Double I_ESP_I2xy3z_D2x = I_ESP_K3xy3z_Px+ABX*I_ESP_I2xy3z_Px;
    Double I_ESP_I2x4z_D2x = I_ESP_K3x4z_Px+ABX*I_ESP_I2x4z_Px;
    Double I_ESP_Ix5y_D2x = I_ESP_K2x5y_Px+ABX*I_ESP_Ix5y_Px;
    Double I_ESP_Ix4yz_D2x = I_ESP_K2x4yz_Px+ABX*I_ESP_Ix4yz_Px;
    Double I_ESP_Ix3y2z_D2x = I_ESP_K2x3y2z_Px+ABX*I_ESP_Ix3y2z_Px;
    Double I_ESP_Ix2y3z_D2x = I_ESP_K2x2y3z_Px+ABX*I_ESP_Ix2y3z_Px;
    Double I_ESP_Ixy4z_D2x = I_ESP_K2xy4z_Px+ABX*I_ESP_Ixy4z_Px;
    Double I_ESP_Ix5z_D2x = I_ESP_K2x5z_Px+ABX*I_ESP_Ix5z_Px;
    Double I_ESP_I6y_D2x = I_ESP_Kx6y_Px+ABX*I_ESP_I6y_Px;
    Double I_ESP_I5yz_D2x = I_ESP_Kx5yz_Px+ABX*I_ESP_I5yz_Px;
    Double I_ESP_I4y2z_D2x = I_ESP_Kx4y2z_Px+ABX*I_ESP_I4y2z_Px;
    Double I_ESP_I3y3z_D2x = I_ESP_Kx3y3z_Px+ABX*I_ESP_I3y3z_Px;
    Double I_ESP_I2y4z_D2x = I_ESP_Kx2y4z_Px+ABX*I_ESP_I2y4z_Px;
    Double I_ESP_Iy5z_D2x = I_ESP_Kxy5z_Px+ABX*I_ESP_Iy5z_Px;
    Double I_ESP_I6z_D2x = I_ESP_Kx6z_Px+ABX*I_ESP_I6z_Px;
    Double I_ESP_I5xy_D2y = I_ESP_K5x2y_Py+ABY*I_ESP_I5xy_Py;
    Double I_ESP_I5xz_D2y = I_ESP_K5xyz_Py+ABY*I_ESP_I5xz_Py;
    Double I_ESP_I4x2y_D2y = I_ESP_K4x3y_Py+ABY*I_ESP_I4x2y_Py;
    Double I_ESP_I4xyz_D2y = I_ESP_K4x2yz_Py+ABY*I_ESP_I4xyz_Py;
    Double I_ESP_I4x2z_D2y = I_ESP_K4xy2z_Py+ABY*I_ESP_I4x2z_Py;
    Double I_ESP_I3x3y_D2y = I_ESP_K3x4y_Py+ABY*I_ESP_I3x3y_Py;
    Double I_ESP_I3x2yz_D2y = I_ESP_K3x3yz_Py+ABY*I_ESP_I3x2yz_Py;
    Double I_ESP_I3xy2z_D2y = I_ESP_K3x2y2z_Py+ABY*I_ESP_I3xy2z_Py;
    Double I_ESP_I3x3z_D2y = I_ESP_K3xy3z_Py+ABY*I_ESP_I3x3z_Py;
    Double I_ESP_I2x4y_D2y = I_ESP_K2x5y_Py+ABY*I_ESP_I2x4y_Py;
    Double I_ESP_I2x3yz_D2y = I_ESP_K2x4yz_Py+ABY*I_ESP_I2x3yz_Py;
    Double I_ESP_I2x2y2z_D2y = I_ESP_K2x3y2z_Py+ABY*I_ESP_I2x2y2z_Py;
    Double I_ESP_I2xy3z_D2y = I_ESP_K2x2y3z_Py+ABY*I_ESP_I2xy3z_Py;
    Double I_ESP_I2x4z_D2y = I_ESP_K2xy4z_Py+ABY*I_ESP_I2x4z_Py;
    Double I_ESP_Ix5y_D2y = I_ESP_Kx6y_Py+ABY*I_ESP_Ix5y_Py;
    Double I_ESP_Ix4yz_D2y = I_ESP_Kx5yz_Py+ABY*I_ESP_Ix4yz_Py;
    Double I_ESP_Ix3y2z_D2y = I_ESP_Kx4y2z_Py+ABY*I_ESP_Ix3y2z_Py;
    Double I_ESP_Ix2y3z_D2y = I_ESP_Kx3y3z_Py+ABY*I_ESP_Ix2y3z_Py;
    Double I_ESP_Ixy4z_D2y = I_ESP_Kx2y4z_Py+ABY*I_ESP_Ixy4z_Py;
    Double I_ESP_Ix5z_D2y = I_ESP_Kxy5z_Py+ABY*I_ESP_Ix5z_Py;
    Double I_ESP_I6y_D2y = I_ESP_K7y_Py+ABY*I_ESP_I6y_Py;
    Double I_ESP_I5yz_D2y = I_ESP_K6yz_Py+ABY*I_ESP_I5yz_Py;
    Double I_ESP_I4y2z_D2y = I_ESP_K5y2z_Py+ABY*I_ESP_I4y2z_Py;
    Double I_ESP_I3y3z_D2y = I_ESP_K4y3z_Py+ABY*I_ESP_I3y3z_Py;
    Double I_ESP_I2y4z_D2y = I_ESP_K3y4z_Py+ABY*I_ESP_I2y4z_Py;
    Double I_ESP_Iy5z_D2y = I_ESP_K2y5z_Py+ABY*I_ESP_Iy5z_Py;
    Double I_ESP_I6z_D2y = I_ESP_Ky6z_Py+ABY*I_ESP_I6z_Py;
    Double I_ESP_I5xy_D2z = I_ESP_K5xyz_Pz+ABZ*I_ESP_I5xy_Pz;
    Double I_ESP_I5xz_D2z = I_ESP_K5x2z_Pz+ABZ*I_ESP_I5xz_Pz;
    Double I_ESP_I4x2y_D2z = I_ESP_K4x2yz_Pz+ABZ*I_ESP_I4x2y_Pz;
    Double I_ESP_I4xyz_D2z = I_ESP_K4xy2z_Pz+ABZ*I_ESP_I4xyz_Pz;
    Double I_ESP_I4x2z_D2z = I_ESP_K4x3z_Pz+ABZ*I_ESP_I4x2z_Pz;
    Double I_ESP_I3x3y_D2z = I_ESP_K3x3yz_Pz+ABZ*I_ESP_I3x3y_Pz;
    Double I_ESP_I3x2yz_D2z = I_ESP_K3x2y2z_Pz+ABZ*I_ESP_I3x2yz_Pz;
    Double I_ESP_I3xy2z_D2z = I_ESP_K3xy3z_Pz+ABZ*I_ESP_I3xy2z_Pz;
    Double I_ESP_I3x3z_D2z = I_ESP_K3x4z_Pz+ABZ*I_ESP_I3x3z_Pz;
    Double I_ESP_I2x4y_D2z = I_ESP_K2x4yz_Pz+ABZ*I_ESP_I2x4y_Pz;
    Double I_ESP_I2x3yz_D2z = I_ESP_K2x3y2z_Pz+ABZ*I_ESP_I2x3yz_Pz;
    Double I_ESP_I2x2y2z_D2z = I_ESP_K2x2y3z_Pz+ABZ*I_ESP_I2x2y2z_Pz;
    Double I_ESP_I2xy3z_D2z = I_ESP_K2xy4z_Pz+ABZ*I_ESP_I2xy3z_Pz;
    Double I_ESP_I2x4z_D2z = I_ESP_K2x5z_Pz+ABZ*I_ESP_I2x4z_Pz;
    Double I_ESP_Ix5y_D2z = I_ESP_Kx5yz_Pz+ABZ*I_ESP_Ix5y_Pz;
    Double I_ESP_Ix4yz_D2z = I_ESP_Kx4y2z_Pz+ABZ*I_ESP_Ix4yz_Pz;
    Double I_ESP_Ix3y2z_D2z = I_ESP_Kx3y3z_Pz+ABZ*I_ESP_Ix3y2z_Pz;
    Double I_ESP_Ix2y3z_D2z = I_ESP_Kx2y4z_Pz+ABZ*I_ESP_Ix2y3z_Pz;
    Double I_ESP_Ixy4z_D2z = I_ESP_Kxy5z_Pz+ABZ*I_ESP_Ixy4z_Pz;
    Double I_ESP_Ix5z_D2z = I_ESP_Kx6z_Pz+ABZ*I_ESP_Ix5z_Pz;
    Double I_ESP_I5yz_D2z = I_ESP_K5y2z_Pz+ABZ*I_ESP_I5yz_Pz;
    Double I_ESP_I4y2z_D2z = I_ESP_K4y3z_Pz+ABZ*I_ESP_I4y2z_Pz;
    Double I_ESP_I3y3z_D2z = I_ESP_K3y4z_Pz+ABZ*I_ESP_I3y3z_Pz;
    Double I_ESP_I2y4z_D2z = I_ESP_K2y5z_Pz+ABZ*I_ESP_I2y4z_Pz;
    Double I_ESP_Iy5z_D2z = I_ESP_Ky6z_Pz+ABZ*I_ESP_Iy5z_Pz;
    Double I_ESP_I6z_D2z = I_ESP_K7z_Pz+ABZ*I_ESP_I6z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 67 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_D
     * RHS shell quartet name: SQ_ESP_H_D
     ************************************************************/
    Double I_ESP_H5x_F3x = I_ESP_I6x_D2x+ABX*I_ESP_H5x_D2x;
    Double I_ESP_H4xy_F3x = I_ESP_I5xy_D2x+ABX*I_ESP_H4xy_D2x;
    Double I_ESP_H4xz_F3x = I_ESP_I5xz_D2x+ABX*I_ESP_H4xz_D2x;
    Double I_ESP_H3x2y_F3x = I_ESP_I4x2y_D2x+ABX*I_ESP_H3x2y_D2x;
    Double I_ESP_H3xyz_F3x = I_ESP_I4xyz_D2x+ABX*I_ESP_H3xyz_D2x;
    Double I_ESP_H3x2z_F3x = I_ESP_I4x2z_D2x+ABX*I_ESP_H3x2z_D2x;
    Double I_ESP_H2x3y_F3x = I_ESP_I3x3y_D2x+ABX*I_ESP_H2x3y_D2x;
    Double I_ESP_H2x2yz_F3x = I_ESP_I3x2yz_D2x+ABX*I_ESP_H2x2yz_D2x;
    Double I_ESP_H2xy2z_F3x = I_ESP_I3xy2z_D2x+ABX*I_ESP_H2xy2z_D2x;
    Double I_ESP_H2x3z_F3x = I_ESP_I3x3z_D2x+ABX*I_ESP_H2x3z_D2x;
    Double I_ESP_Hx4y_F3x = I_ESP_I2x4y_D2x+ABX*I_ESP_Hx4y_D2x;
    Double I_ESP_Hx3yz_F3x = I_ESP_I2x3yz_D2x+ABX*I_ESP_Hx3yz_D2x;
    Double I_ESP_Hx2y2z_F3x = I_ESP_I2x2y2z_D2x+ABX*I_ESP_Hx2y2z_D2x;
    Double I_ESP_Hxy3z_F3x = I_ESP_I2xy3z_D2x+ABX*I_ESP_Hxy3z_D2x;
    Double I_ESP_Hx4z_F3x = I_ESP_I2x4z_D2x+ABX*I_ESP_Hx4z_D2x;
    Double I_ESP_H5y_F3x = I_ESP_Ix5y_D2x+ABX*I_ESP_H5y_D2x;
    Double I_ESP_H4yz_F3x = I_ESP_Ix4yz_D2x+ABX*I_ESP_H4yz_D2x;
    Double I_ESP_H3y2z_F3x = I_ESP_Ix3y2z_D2x+ABX*I_ESP_H3y2z_D2x;
    Double I_ESP_H2y3z_F3x = I_ESP_Ix2y3z_D2x+ABX*I_ESP_H2y3z_D2x;
    Double I_ESP_Hy4z_F3x = I_ESP_Ixy4z_D2x+ABX*I_ESP_Hy4z_D2x;
    Double I_ESP_H5z_F3x = I_ESP_Ix5z_D2x+ABX*I_ESP_H5z_D2x;
    Double I_ESP_H4xy_F2xy = I_ESP_I4x2y_D2x+ABY*I_ESP_H4xy_D2x;
    Double I_ESP_H4xz_F2xy = I_ESP_I4xyz_D2x+ABY*I_ESP_H4xz_D2x;
    Double I_ESP_H3x2y_F2xy = I_ESP_I3x3y_D2x+ABY*I_ESP_H3x2y_D2x;
    Double I_ESP_H3xyz_F2xy = I_ESP_I3x2yz_D2x+ABY*I_ESP_H3xyz_D2x;
    Double I_ESP_H3x2z_F2xy = I_ESP_I3xy2z_D2x+ABY*I_ESP_H3x2z_D2x;
    Double I_ESP_H2x3y_F2xy = I_ESP_I2x4y_D2x+ABY*I_ESP_H2x3y_D2x;
    Double I_ESP_H2x2yz_F2xy = I_ESP_I2x3yz_D2x+ABY*I_ESP_H2x2yz_D2x;
    Double I_ESP_H2xy2z_F2xy = I_ESP_I2x2y2z_D2x+ABY*I_ESP_H2xy2z_D2x;
    Double I_ESP_H2x3z_F2xy = I_ESP_I2xy3z_D2x+ABY*I_ESP_H2x3z_D2x;
    Double I_ESP_Hx4y_F2xy = I_ESP_Ix5y_D2x+ABY*I_ESP_Hx4y_D2x;
    Double I_ESP_Hx3yz_F2xy = I_ESP_Ix4yz_D2x+ABY*I_ESP_Hx3yz_D2x;
    Double I_ESP_Hx2y2z_F2xy = I_ESP_Ix3y2z_D2x+ABY*I_ESP_Hx2y2z_D2x;
    Double I_ESP_Hxy3z_F2xy = I_ESP_Ix2y3z_D2x+ABY*I_ESP_Hxy3z_D2x;
    Double I_ESP_Hx4z_F2xy = I_ESP_Ixy4z_D2x+ABY*I_ESP_Hx4z_D2x;
    Double I_ESP_H5y_F2xy = I_ESP_I6y_D2x+ABY*I_ESP_H5y_D2x;
    Double I_ESP_H4yz_F2xy = I_ESP_I5yz_D2x+ABY*I_ESP_H4yz_D2x;
    Double I_ESP_H3y2z_F2xy = I_ESP_I4y2z_D2x+ABY*I_ESP_H3y2z_D2x;
    Double I_ESP_H2y3z_F2xy = I_ESP_I3y3z_D2x+ABY*I_ESP_H2y3z_D2x;
    Double I_ESP_Hy4z_F2xy = I_ESP_I2y4z_D2x+ABY*I_ESP_Hy4z_D2x;
    Double I_ESP_H5z_F2xy = I_ESP_Iy5z_D2x+ABY*I_ESP_H5z_D2x;
    Double I_ESP_H4xz_F2xz = I_ESP_I4x2z_D2x+ABZ*I_ESP_H4xz_D2x;
    Double I_ESP_H3xyz_F2xz = I_ESP_I3xy2z_D2x+ABZ*I_ESP_H3xyz_D2x;
    Double I_ESP_H3x2z_F2xz = I_ESP_I3x3z_D2x+ABZ*I_ESP_H3x2z_D2x;
    Double I_ESP_H2x2yz_F2xz = I_ESP_I2x2y2z_D2x+ABZ*I_ESP_H2x2yz_D2x;
    Double I_ESP_H2xy2z_F2xz = I_ESP_I2xy3z_D2x+ABZ*I_ESP_H2xy2z_D2x;
    Double I_ESP_H2x3z_F2xz = I_ESP_I2x4z_D2x+ABZ*I_ESP_H2x3z_D2x;
    Double I_ESP_Hx3yz_F2xz = I_ESP_Ix3y2z_D2x+ABZ*I_ESP_Hx3yz_D2x;
    Double I_ESP_Hx2y2z_F2xz = I_ESP_Ix2y3z_D2x+ABZ*I_ESP_Hx2y2z_D2x;
    Double I_ESP_Hxy3z_F2xz = I_ESP_Ixy4z_D2x+ABZ*I_ESP_Hxy3z_D2x;
    Double I_ESP_Hx4z_F2xz = I_ESP_Ix5z_D2x+ABZ*I_ESP_Hx4z_D2x;
    Double I_ESP_H4yz_F2xz = I_ESP_I4y2z_D2x+ABZ*I_ESP_H4yz_D2x;
    Double I_ESP_H3y2z_F2xz = I_ESP_I3y3z_D2x+ABZ*I_ESP_H3y2z_D2x;
    Double I_ESP_H2y3z_F2xz = I_ESP_I2y4z_D2x+ABZ*I_ESP_H2y3z_D2x;
    Double I_ESP_Hy4z_F2xz = I_ESP_Iy5z_D2x+ABZ*I_ESP_Hy4z_D2x;
    Double I_ESP_H5z_F2xz = I_ESP_I6z_D2x+ABZ*I_ESP_H5z_D2x;
    Double I_ESP_H4xz_Fx2y = I_ESP_I5xz_D2y+ABX*I_ESP_H4xz_D2y;
    Double I_ESP_H3xyz_Fx2y = I_ESP_I4xyz_D2y+ABX*I_ESP_H3xyz_D2y;
    Double I_ESP_H3x2z_Fx2y = I_ESP_I4x2z_D2y+ABX*I_ESP_H3x2z_D2y;
    Double I_ESP_H2x2yz_Fx2y = I_ESP_I3x2yz_D2y+ABX*I_ESP_H2x2yz_D2y;
    Double I_ESP_H2xy2z_Fx2y = I_ESP_I3xy2z_D2y+ABX*I_ESP_H2xy2z_D2y;
    Double I_ESP_H2x3z_Fx2y = I_ESP_I3x3z_D2y+ABX*I_ESP_H2x3z_D2y;
    Double I_ESP_Hx3yz_Fx2y = I_ESP_I2x3yz_D2y+ABX*I_ESP_Hx3yz_D2y;
    Double I_ESP_Hx2y2z_Fx2y = I_ESP_I2x2y2z_D2y+ABX*I_ESP_Hx2y2z_D2y;
    Double I_ESP_Hxy3z_Fx2y = I_ESP_I2xy3z_D2y+ABX*I_ESP_Hxy3z_D2y;
    Double I_ESP_Hx4z_Fx2y = I_ESP_I2x4z_D2y+ABX*I_ESP_Hx4z_D2y;
    Double I_ESP_H4yz_Fx2y = I_ESP_Ix4yz_D2y+ABX*I_ESP_H4yz_D2y;
    Double I_ESP_H3y2z_Fx2y = I_ESP_Ix3y2z_D2y+ABX*I_ESP_H3y2z_D2y;
    Double I_ESP_H2y3z_Fx2y = I_ESP_Ix2y3z_D2y+ABX*I_ESP_H2y3z_D2y;
    Double I_ESP_Hy4z_Fx2y = I_ESP_Ixy4z_D2y+ABX*I_ESP_Hy4z_D2y;
    Double I_ESP_H5z_Fx2y = I_ESP_Ix5z_D2y+ABX*I_ESP_H5z_D2y;
    Double I_ESP_H4xy_Fx2z = I_ESP_I5xy_D2z+ABX*I_ESP_H4xy_D2z;
    Double I_ESP_H3x2y_Fx2z = I_ESP_I4x2y_D2z+ABX*I_ESP_H3x2y_D2z;
    Double I_ESP_H3xyz_Fx2z = I_ESP_I4xyz_D2z+ABX*I_ESP_H3xyz_D2z;
    Double I_ESP_H2x3y_Fx2z = I_ESP_I3x3y_D2z+ABX*I_ESP_H2x3y_D2z;
    Double I_ESP_H2x2yz_Fx2z = I_ESP_I3x2yz_D2z+ABX*I_ESP_H2x2yz_D2z;
    Double I_ESP_H2xy2z_Fx2z = I_ESP_I3xy2z_D2z+ABX*I_ESP_H2xy2z_D2z;
    Double I_ESP_Hx4y_Fx2z = I_ESP_I2x4y_D2z+ABX*I_ESP_Hx4y_D2z;
    Double I_ESP_Hx3yz_Fx2z = I_ESP_I2x3yz_D2z+ABX*I_ESP_Hx3yz_D2z;
    Double I_ESP_Hx2y2z_Fx2z = I_ESP_I2x2y2z_D2z+ABX*I_ESP_Hx2y2z_D2z;
    Double I_ESP_Hxy3z_Fx2z = I_ESP_I2xy3z_D2z+ABX*I_ESP_Hxy3z_D2z;
    Double I_ESP_H5y_Fx2z = I_ESP_Ix5y_D2z+ABX*I_ESP_H5y_D2z;
    Double I_ESP_H4yz_Fx2z = I_ESP_Ix4yz_D2z+ABX*I_ESP_H4yz_D2z;
    Double I_ESP_H3y2z_Fx2z = I_ESP_Ix3y2z_D2z+ABX*I_ESP_H3y2z_D2z;
    Double I_ESP_H2y3z_Fx2z = I_ESP_Ix2y3z_D2z+ABX*I_ESP_H2y3z_D2z;
    Double I_ESP_Hy4z_Fx2z = I_ESP_Ixy4z_D2z+ABX*I_ESP_Hy4z_D2z;
    Double I_ESP_H5x_F3y = I_ESP_I5xy_D2y+ABY*I_ESP_H5x_D2y;
    Double I_ESP_H4xy_F3y = I_ESP_I4x2y_D2y+ABY*I_ESP_H4xy_D2y;
    Double I_ESP_H4xz_F3y = I_ESP_I4xyz_D2y+ABY*I_ESP_H4xz_D2y;
    Double I_ESP_H3x2y_F3y = I_ESP_I3x3y_D2y+ABY*I_ESP_H3x2y_D2y;
    Double I_ESP_H3xyz_F3y = I_ESP_I3x2yz_D2y+ABY*I_ESP_H3xyz_D2y;
    Double I_ESP_H3x2z_F3y = I_ESP_I3xy2z_D2y+ABY*I_ESP_H3x2z_D2y;
    Double I_ESP_H2x3y_F3y = I_ESP_I2x4y_D2y+ABY*I_ESP_H2x3y_D2y;
    Double I_ESP_H2x2yz_F3y = I_ESP_I2x3yz_D2y+ABY*I_ESP_H2x2yz_D2y;
    Double I_ESP_H2xy2z_F3y = I_ESP_I2x2y2z_D2y+ABY*I_ESP_H2xy2z_D2y;
    Double I_ESP_H2x3z_F3y = I_ESP_I2xy3z_D2y+ABY*I_ESP_H2x3z_D2y;
    Double I_ESP_Hx4y_F3y = I_ESP_Ix5y_D2y+ABY*I_ESP_Hx4y_D2y;
    Double I_ESP_Hx3yz_F3y = I_ESP_Ix4yz_D2y+ABY*I_ESP_Hx3yz_D2y;
    Double I_ESP_Hx2y2z_F3y = I_ESP_Ix3y2z_D2y+ABY*I_ESP_Hx2y2z_D2y;
    Double I_ESP_Hxy3z_F3y = I_ESP_Ix2y3z_D2y+ABY*I_ESP_Hxy3z_D2y;
    Double I_ESP_Hx4z_F3y = I_ESP_Ixy4z_D2y+ABY*I_ESP_Hx4z_D2y;
    Double I_ESP_H5y_F3y = I_ESP_I6y_D2y+ABY*I_ESP_H5y_D2y;
    Double I_ESP_H4yz_F3y = I_ESP_I5yz_D2y+ABY*I_ESP_H4yz_D2y;
    Double I_ESP_H3y2z_F3y = I_ESP_I4y2z_D2y+ABY*I_ESP_H3y2z_D2y;
    Double I_ESP_H2y3z_F3y = I_ESP_I3y3z_D2y+ABY*I_ESP_H2y3z_D2y;
    Double I_ESP_Hy4z_F3y = I_ESP_I2y4z_D2y+ABY*I_ESP_Hy4z_D2y;
    Double I_ESP_H5z_F3y = I_ESP_Iy5z_D2y+ABY*I_ESP_H5z_D2y;
    Double I_ESP_H4xz_F2yz = I_ESP_I4x2z_D2y+ABZ*I_ESP_H4xz_D2y;
    Double I_ESP_H3xyz_F2yz = I_ESP_I3xy2z_D2y+ABZ*I_ESP_H3xyz_D2y;
    Double I_ESP_H3x2z_F2yz = I_ESP_I3x3z_D2y+ABZ*I_ESP_H3x2z_D2y;
    Double I_ESP_H2x2yz_F2yz = I_ESP_I2x2y2z_D2y+ABZ*I_ESP_H2x2yz_D2y;
    Double I_ESP_H2xy2z_F2yz = I_ESP_I2xy3z_D2y+ABZ*I_ESP_H2xy2z_D2y;
    Double I_ESP_H2x3z_F2yz = I_ESP_I2x4z_D2y+ABZ*I_ESP_H2x3z_D2y;
    Double I_ESP_Hx3yz_F2yz = I_ESP_Ix3y2z_D2y+ABZ*I_ESP_Hx3yz_D2y;
    Double I_ESP_Hx2y2z_F2yz = I_ESP_Ix2y3z_D2y+ABZ*I_ESP_Hx2y2z_D2y;
    Double I_ESP_Hxy3z_F2yz = I_ESP_Ixy4z_D2y+ABZ*I_ESP_Hxy3z_D2y;
    Double I_ESP_Hx4z_F2yz = I_ESP_Ix5z_D2y+ABZ*I_ESP_Hx4z_D2y;
    Double I_ESP_H4yz_F2yz = I_ESP_I4y2z_D2y+ABZ*I_ESP_H4yz_D2y;
    Double I_ESP_H3y2z_F2yz = I_ESP_I3y3z_D2y+ABZ*I_ESP_H3y2z_D2y;
    Double I_ESP_H2y3z_F2yz = I_ESP_I2y4z_D2y+ABZ*I_ESP_H2y3z_D2y;
    Double I_ESP_Hy4z_F2yz = I_ESP_Iy5z_D2y+ABZ*I_ESP_Hy4z_D2y;
    Double I_ESP_H5z_F2yz = I_ESP_I6z_D2y+ABZ*I_ESP_H5z_D2y;
    Double I_ESP_H5x_F3z = I_ESP_I5xz_D2z+ABZ*I_ESP_H5x_D2z;
    Double I_ESP_H4xy_F3z = I_ESP_I4xyz_D2z+ABZ*I_ESP_H4xy_D2z;
    Double I_ESP_H4xz_F3z = I_ESP_I4x2z_D2z+ABZ*I_ESP_H4xz_D2z;
    Double I_ESP_H3x2y_F3z = I_ESP_I3x2yz_D2z+ABZ*I_ESP_H3x2y_D2z;
    Double I_ESP_H3xyz_F3z = I_ESP_I3xy2z_D2z+ABZ*I_ESP_H3xyz_D2z;
    Double I_ESP_H3x2z_F3z = I_ESP_I3x3z_D2z+ABZ*I_ESP_H3x2z_D2z;
    Double I_ESP_H2x3y_F3z = I_ESP_I2x3yz_D2z+ABZ*I_ESP_H2x3y_D2z;
    Double I_ESP_H2x2yz_F3z = I_ESP_I2x2y2z_D2z+ABZ*I_ESP_H2x2yz_D2z;
    Double I_ESP_H2xy2z_F3z = I_ESP_I2xy3z_D2z+ABZ*I_ESP_H2xy2z_D2z;
    Double I_ESP_H2x3z_F3z = I_ESP_I2x4z_D2z+ABZ*I_ESP_H2x3z_D2z;
    Double I_ESP_Hx4y_F3z = I_ESP_Ix4yz_D2z+ABZ*I_ESP_Hx4y_D2z;
    Double I_ESP_Hx3yz_F3z = I_ESP_Ix3y2z_D2z+ABZ*I_ESP_Hx3yz_D2z;
    Double I_ESP_Hx2y2z_F3z = I_ESP_Ix2y3z_D2z+ABZ*I_ESP_Hx2y2z_D2z;
    Double I_ESP_Hxy3z_F3z = I_ESP_Ixy4z_D2z+ABZ*I_ESP_Hxy3z_D2z;
    Double I_ESP_Hx4z_F3z = I_ESP_Ix5z_D2z+ABZ*I_ESP_Hx4z_D2z;
    Double I_ESP_H5y_F3z = I_ESP_I5yz_D2z+ABZ*I_ESP_H5y_D2z;
    Double I_ESP_H4yz_F3z = I_ESP_I4y2z_D2z+ABZ*I_ESP_H4yz_D2z;
    Double I_ESP_H3y2z_F3z = I_ESP_I3y3z_D2z+ABZ*I_ESP_H3y2z_D2z;
    Double I_ESP_H2y3z_F3z = I_ESP_I2y4z_D2z+ABZ*I_ESP_H2y3z_D2z;
    Double I_ESP_Hy4z_F3z = I_ESP_Iy5z_D2z+ABZ*I_ESP_Hy4z_D2z;
    Double I_ESP_H5z_F3z = I_ESP_I6z_D2z+ABZ*I_ESP_H5z_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_G_G
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_F
     * RHS shell quartet name: SQ_ESP_G_F
     ************************************************************/
    Double I_ESP_G4x_G4x = I_ESP_H5x_F3x+ABX*I_ESP_G4x_F3x;
    Double I_ESP_G3xy_G4x = I_ESP_H4xy_F3x+ABX*I_ESP_G3xy_F3x;
    Double I_ESP_G3xz_G4x = I_ESP_H4xz_F3x+ABX*I_ESP_G3xz_F3x;
    Double I_ESP_G2x2y_G4x = I_ESP_H3x2y_F3x+ABX*I_ESP_G2x2y_F3x;
    Double I_ESP_G2xyz_G4x = I_ESP_H3xyz_F3x+ABX*I_ESP_G2xyz_F3x;
    Double I_ESP_G2x2z_G4x = I_ESP_H3x2z_F3x+ABX*I_ESP_G2x2z_F3x;
    Double I_ESP_Gx3y_G4x = I_ESP_H2x3y_F3x+ABX*I_ESP_Gx3y_F3x;
    Double I_ESP_Gx2yz_G4x = I_ESP_H2x2yz_F3x+ABX*I_ESP_Gx2yz_F3x;
    Double I_ESP_Gxy2z_G4x = I_ESP_H2xy2z_F3x+ABX*I_ESP_Gxy2z_F3x;
    Double I_ESP_Gx3z_G4x = I_ESP_H2x3z_F3x+ABX*I_ESP_Gx3z_F3x;
    Double I_ESP_G4y_G4x = I_ESP_Hx4y_F3x+ABX*I_ESP_G4y_F3x;
    Double I_ESP_G3yz_G4x = I_ESP_Hx3yz_F3x+ABX*I_ESP_G3yz_F3x;
    Double I_ESP_G2y2z_G4x = I_ESP_Hx2y2z_F3x+ABX*I_ESP_G2y2z_F3x;
    Double I_ESP_Gy3z_G4x = I_ESP_Hxy3z_F3x+ABX*I_ESP_Gy3z_F3x;
    Double I_ESP_G4z_G4x = I_ESP_Hx4z_F3x+ABX*I_ESP_G4z_F3x;
    Double I_ESP_G4x_G3xy = I_ESP_H4xy_F3x+ABY*I_ESP_G4x_F3x;
    Double I_ESP_G3xy_G3xy = I_ESP_H3x2y_F3x+ABY*I_ESP_G3xy_F3x;
    Double I_ESP_G3xz_G3xy = I_ESP_H3xyz_F3x+ABY*I_ESP_G3xz_F3x;
    Double I_ESP_G2x2y_G3xy = I_ESP_H2x3y_F3x+ABY*I_ESP_G2x2y_F3x;
    Double I_ESP_G2xyz_G3xy = I_ESP_H2x2yz_F3x+ABY*I_ESP_G2xyz_F3x;
    Double I_ESP_G2x2z_G3xy = I_ESP_H2xy2z_F3x+ABY*I_ESP_G2x2z_F3x;
    Double I_ESP_Gx3y_G3xy = I_ESP_Hx4y_F3x+ABY*I_ESP_Gx3y_F3x;
    Double I_ESP_Gx2yz_G3xy = I_ESP_Hx3yz_F3x+ABY*I_ESP_Gx2yz_F3x;
    Double I_ESP_Gxy2z_G3xy = I_ESP_Hx2y2z_F3x+ABY*I_ESP_Gxy2z_F3x;
    Double I_ESP_Gx3z_G3xy = I_ESP_Hxy3z_F3x+ABY*I_ESP_Gx3z_F3x;
    Double I_ESP_G4y_G3xy = I_ESP_H5y_F3x+ABY*I_ESP_G4y_F3x;
    Double I_ESP_G3yz_G3xy = I_ESP_H4yz_F3x+ABY*I_ESP_G3yz_F3x;
    Double I_ESP_G2y2z_G3xy = I_ESP_H3y2z_F3x+ABY*I_ESP_G2y2z_F3x;
    Double I_ESP_Gy3z_G3xy = I_ESP_H2y3z_F3x+ABY*I_ESP_Gy3z_F3x;
    Double I_ESP_G4z_G3xy = I_ESP_Hy4z_F3x+ABY*I_ESP_G4z_F3x;
    Double I_ESP_G4x_G3xz = I_ESP_H4xz_F3x+ABZ*I_ESP_G4x_F3x;
    Double I_ESP_G3xy_G3xz = I_ESP_H3xyz_F3x+ABZ*I_ESP_G3xy_F3x;
    Double I_ESP_G3xz_G3xz = I_ESP_H3x2z_F3x+ABZ*I_ESP_G3xz_F3x;
    Double I_ESP_G2x2y_G3xz = I_ESP_H2x2yz_F3x+ABZ*I_ESP_G2x2y_F3x;
    Double I_ESP_G2xyz_G3xz = I_ESP_H2xy2z_F3x+ABZ*I_ESP_G2xyz_F3x;
    Double I_ESP_G2x2z_G3xz = I_ESP_H2x3z_F3x+ABZ*I_ESP_G2x2z_F3x;
    Double I_ESP_Gx3y_G3xz = I_ESP_Hx3yz_F3x+ABZ*I_ESP_Gx3y_F3x;
    Double I_ESP_Gx2yz_G3xz = I_ESP_Hx2y2z_F3x+ABZ*I_ESP_Gx2yz_F3x;
    Double I_ESP_Gxy2z_G3xz = I_ESP_Hxy3z_F3x+ABZ*I_ESP_Gxy2z_F3x;
    Double I_ESP_Gx3z_G3xz = I_ESP_Hx4z_F3x+ABZ*I_ESP_Gx3z_F3x;
    Double I_ESP_G4y_G3xz = I_ESP_H4yz_F3x+ABZ*I_ESP_G4y_F3x;
    Double I_ESP_G3yz_G3xz = I_ESP_H3y2z_F3x+ABZ*I_ESP_G3yz_F3x;
    Double I_ESP_G2y2z_G3xz = I_ESP_H2y3z_F3x+ABZ*I_ESP_G2y2z_F3x;
    Double I_ESP_Gy3z_G3xz = I_ESP_Hy4z_F3x+ABZ*I_ESP_Gy3z_F3x;
    Double I_ESP_G4z_G3xz = I_ESP_H5z_F3x+ABZ*I_ESP_G4z_F3x;
    Double I_ESP_G4x_G2x2y = I_ESP_H4xy_F2xy+ABY*I_ESP_G4x_F2xy;
    Double I_ESP_G3xy_G2x2y = I_ESP_H3x2y_F2xy+ABY*I_ESP_G3xy_F2xy;
    Double I_ESP_G3xz_G2x2y = I_ESP_H3xyz_F2xy+ABY*I_ESP_G3xz_F2xy;
    Double I_ESP_G2x2y_G2x2y = I_ESP_H2x3y_F2xy+ABY*I_ESP_G2x2y_F2xy;
    Double I_ESP_G2xyz_G2x2y = I_ESP_H2x2yz_F2xy+ABY*I_ESP_G2xyz_F2xy;
    Double I_ESP_G2x2z_G2x2y = I_ESP_H2xy2z_F2xy+ABY*I_ESP_G2x2z_F2xy;
    Double I_ESP_Gx3y_G2x2y = I_ESP_Hx4y_F2xy+ABY*I_ESP_Gx3y_F2xy;
    Double I_ESP_Gx2yz_G2x2y = I_ESP_Hx3yz_F2xy+ABY*I_ESP_Gx2yz_F2xy;
    Double I_ESP_Gxy2z_G2x2y = I_ESP_Hx2y2z_F2xy+ABY*I_ESP_Gxy2z_F2xy;
    Double I_ESP_Gx3z_G2x2y = I_ESP_Hxy3z_F2xy+ABY*I_ESP_Gx3z_F2xy;
    Double I_ESP_G4y_G2x2y = I_ESP_H5y_F2xy+ABY*I_ESP_G4y_F2xy;
    Double I_ESP_G3yz_G2x2y = I_ESP_H4yz_F2xy+ABY*I_ESP_G3yz_F2xy;
    Double I_ESP_G2y2z_G2x2y = I_ESP_H3y2z_F2xy+ABY*I_ESP_G2y2z_F2xy;
    Double I_ESP_Gy3z_G2x2y = I_ESP_H2y3z_F2xy+ABY*I_ESP_Gy3z_F2xy;
    Double I_ESP_G4z_G2x2y = I_ESP_Hy4z_F2xy+ABY*I_ESP_G4z_F2xy;
    Double I_ESP_G4x_G2xyz = I_ESP_H4xz_F2xy+ABZ*I_ESP_G4x_F2xy;
    Double I_ESP_G3xy_G2xyz = I_ESP_H3xyz_F2xy+ABZ*I_ESP_G3xy_F2xy;
    Double I_ESP_G3xz_G2xyz = I_ESP_H3x2z_F2xy+ABZ*I_ESP_G3xz_F2xy;
    Double I_ESP_G2x2y_G2xyz = I_ESP_H2x2yz_F2xy+ABZ*I_ESP_G2x2y_F2xy;
    Double I_ESP_G2xyz_G2xyz = I_ESP_H2xy2z_F2xy+ABZ*I_ESP_G2xyz_F2xy;
    Double I_ESP_G2x2z_G2xyz = I_ESP_H2x3z_F2xy+ABZ*I_ESP_G2x2z_F2xy;
    Double I_ESP_Gx3y_G2xyz = I_ESP_Hx3yz_F2xy+ABZ*I_ESP_Gx3y_F2xy;
    Double I_ESP_Gx2yz_G2xyz = I_ESP_Hx2y2z_F2xy+ABZ*I_ESP_Gx2yz_F2xy;
    Double I_ESP_Gxy2z_G2xyz = I_ESP_Hxy3z_F2xy+ABZ*I_ESP_Gxy2z_F2xy;
    Double I_ESP_Gx3z_G2xyz = I_ESP_Hx4z_F2xy+ABZ*I_ESP_Gx3z_F2xy;
    Double I_ESP_G4y_G2xyz = I_ESP_H4yz_F2xy+ABZ*I_ESP_G4y_F2xy;
    Double I_ESP_G3yz_G2xyz = I_ESP_H3y2z_F2xy+ABZ*I_ESP_G3yz_F2xy;
    Double I_ESP_G2y2z_G2xyz = I_ESP_H2y3z_F2xy+ABZ*I_ESP_G2y2z_F2xy;
    Double I_ESP_Gy3z_G2xyz = I_ESP_Hy4z_F2xy+ABZ*I_ESP_Gy3z_F2xy;
    Double I_ESP_G4z_G2xyz = I_ESP_H5z_F2xy+ABZ*I_ESP_G4z_F2xy;
    Double I_ESP_G4x_G2x2z = I_ESP_H4xz_F2xz+ABZ*I_ESP_G4x_F2xz;
    Double I_ESP_G3xy_G2x2z = I_ESP_H3xyz_F2xz+ABZ*I_ESP_G3xy_F2xz;
    Double I_ESP_G3xz_G2x2z = I_ESP_H3x2z_F2xz+ABZ*I_ESP_G3xz_F2xz;
    Double I_ESP_G2x2y_G2x2z = I_ESP_H2x2yz_F2xz+ABZ*I_ESP_G2x2y_F2xz;
    Double I_ESP_G2xyz_G2x2z = I_ESP_H2xy2z_F2xz+ABZ*I_ESP_G2xyz_F2xz;
    Double I_ESP_G2x2z_G2x2z = I_ESP_H2x3z_F2xz+ABZ*I_ESP_G2x2z_F2xz;
    Double I_ESP_Gx3y_G2x2z = I_ESP_Hx3yz_F2xz+ABZ*I_ESP_Gx3y_F2xz;
    Double I_ESP_Gx2yz_G2x2z = I_ESP_Hx2y2z_F2xz+ABZ*I_ESP_Gx2yz_F2xz;
    Double I_ESP_Gxy2z_G2x2z = I_ESP_Hxy3z_F2xz+ABZ*I_ESP_Gxy2z_F2xz;
    Double I_ESP_Gx3z_G2x2z = I_ESP_Hx4z_F2xz+ABZ*I_ESP_Gx3z_F2xz;
    Double I_ESP_G4y_G2x2z = I_ESP_H4yz_F2xz+ABZ*I_ESP_G4y_F2xz;
    Double I_ESP_G3yz_G2x2z = I_ESP_H3y2z_F2xz+ABZ*I_ESP_G3yz_F2xz;
    Double I_ESP_G2y2z_G2x2z = I_ESP_H2y3z_F2xz+ABZ*I_ESP_G2y2z_F2xz;
    Double I_ESP_Gy3z_G2x2z = I_ESP_Hy4z_F2xz+ABZ*I_ESP_Gy3z_F2xz;
    Double I_ESP_G4z_G2x2z = I_ESP_H5z_F2xz+ABZ*I_ESP_G4z_F2xz;
    Double I_ESP_G4x_Gx3y = I_ESP_H5x_F3y+ABX*I_ESP_G4x_F3y;
    Double I_ESP_G3xy_Gx3y = I_ESP_H4xy_F3y+ABX*I_ESP_G3xy_F3y;
    Double I_ESP_G3xz_Gx3y = I_ESP_H4xz_F3y+ABX*I_ESP_G3xz_F3y;
    Double I_ESP_G2x2y_Gx3y = I_ESP_H3x2y_F3y+ABX*I_ESP_G2x2y_F3y;
    Double I_ESP_G2xyz_Gx3y = I_ESP_H3xyz_F3y+ABX*I_ESP_G2xyz_F3y;
    Double I_ESP_G2x2z_Gx3y = I_ESP_H3x2z_F3y+ABX*I_ESP_G2x2z_F3y;
    Double I_ESP_Gx3y_Gx3y = I_ESP_H2x3y_F3y+ABX*I_ESP_Gx3y_F3y;
    Double I_ESP_Gx2yz_Gx3y = I_ESP_H2x2yz_F3y+ABX*I_ESP_Gx2yz_F3y;
    Double I_ESP_Gxy2z_Gx3y = I_ESP_H2xy2z_F3y+ABX*I_ESP_Gxy2z_F3y;
    Double I_ESP_Gx3z_Gx3y = I_ESP_H2x3z_F3y+ABX*I_ESP_Gx3z_F3y;
    Double I_ESP_G4y_Gx3y = I_ESP_Hx4y_F3y+ABX*I_ESP_G4y_F3y;
    Double I_ESP_G3yz_Gx3y = I_ESP_Hx3yz_F3y+ABX*I_ESP_G3yz_F3y;
    Double I_ESP_G2y2z_Gx3y = I_ESP_Hx2y2z_F3y+ABX*I_ESP_G2y2z_F3y;
    Double I_ESP_Gy3z_Gx3y = I_ESP_Hxy3z_F3y+ABX*I_ESP_Gy3z_F3y;
    Double I_ESP_G4z_Gx3y = I_ESP_Hx4z_F3y+ABX*I_ESP_G4z_F3y;
    Double I_ESP_G4x_Gx2yz = I_ESP_H4xz_Fx2y+ABZ*I_ESP_G4x_Fx2y;
    Double I_ESP_G3xy_Gx2yz = I_ESP_H3xyz_Fx2y+ABZ*I_ESP_G3xy_Fx2y;
    Double I_ESP_G3xz_Gx2yz = I_ESP_H3x2z_Fx2y+ABZ*I_ESP_G3xz_Fx2y;
    Double I_ESP_G2x2y_Gx2yz = I_ESP_H2x2yz_Fx2y+ABZ*I_ESP_G2x2y_Fx2y;
    Double I_ESP_G2xyz_Gx2yz = I_ESP_H2xy2z_Fx2y+ABZ*I_ESP_G2xyz_Fx2y;
    Double I_ESP_G2x2z_Gx2yz = I_ESP_H2x3z_Fx2y+ABZ*I_ESP_G2x2z_Fx2y;
    Double I_ESP_Gx3y_Gx2yz = I_ESP_Hx3yz_Fx2y+ABZ*I_ESP_Gx3y_Fx2y;
    Double I_ESP_Gx2yz_Gx2yz = I_ESP_Hx2y2z_Fx2y+ABZ*I_ESP_Gx2yz_Fx2y;
    Double I_ESP_Gxy2z_Gx2yz = I_ESP_Hxy3z_Fx2y+ABZ*I_ESP_Gxy2z_Fx2y;
    Double I_ESP_Gx3z_Gx2yz = I_ESP_Hx4z_Fx2y+ABZ*I_ESP_Gx3z_Fx2y;
    Double I_ESP_G4y_Gx2yz = I_ESP_H4yz_Fx2y+ABZ*I_ESP_G4y_Fx2y;
    Double I_ESP_G3yz_Gx2yz = I_ESP_H3y2z_Fx2y+ABZ*I_ESP_G3yz_Fx2y;
    Double I_ESP_G2y2z_Gx2yz = I_ESP_H2y3z_Fx2y+ABZ*I_ESP_G2y2z_Fx2y;
    Double I_ESP_Gy3z_Gx2yz = I_ESP_Hy4z_Fx2y+ABZ*I_ESP_Gy3z_Fx2y;
    Double I_ESP_G4z_Gx2yz = I_ESP_H5z_Fx2y+ABZ*I_ESP_G4z_Fx2y;
    Double I_ESP_G4x_Gxy2z = I_ESP_H4xy_Fx2z+ABY*I_ESP_G4x_Fx2z;
    Double I_ESP_G3xy_Gxy2z = I_ESP_H3x2y_Fx2z+ABY*I_ESP_G3xy_Fx2z;
    Double I_ESP_G3xz_Gxy2z = I_ESP_H3xyz_Fx2z+ABY*I_ESP_G3xz_Fx2z;
    Double I_ESP_G2x2y_Gxy2z = I_ESP_H2x3y_Fx2z+ABY*I_ESP_G2x2y_Fx2z;
    Double I_ESP_G2xyz_Gxy2z = I_ESP_H2x2yz_Fx2z+ABY*I_ESP_G2xyz_Fx2z;
    Double I_ESP_G2x2z_Gxy2z = I_ESP_H2xy2z_Fx2z+ABY*I_ESP_G2x2z_Fx2z;
    Double I_ESP_Gx3y_Gxy2z = I_ESP_Hx4y_Fx2z+ABY*I_ESP_Gx3y_Fx2z;
    Double I_ESP_Gx2yz_Gxy2z = I_ESP_Hx3yz_Fx2z+ABY*I_ESP_Gx2yz_Fx2z;
    Double I_ESP_Gxy2z_Gxy2z = I_ESP_Hx2y2z_Fx2z+ABY*I_ESP_Gxy2z_Fx2z;
    Double I_ESP_Gx3z_Gxy2z = I_ESP_Hxy3z_Fx2z+ABY*I_ESP_Gx3z_Fx2z;
    Double I_ESP_G4y_Gxy2z = I_ESP_H5y_Fx2z+ABY*I_ESP_G4y_Fx2z;
    Double I_ESP_G3yz_Gxy2z = I_ESP_H4yz_Fx2z+ABY*I_ESP_G3yz_Fx2z;
    Double I_ESP_G2y2z_Gxy2z = I_ESP_H3y2z_Fx2z+ABY*I_ESP_G2y2z_Fx2z;
    Double I_ESP_Gy3z_Gxy2z = I_ESP_H2y3z_Fx2z+ABY*I_ESP_Gy3z_Fx2z;
    Double I_ESP_G4z_Gxy2z = I_ESP_Hy4z_Fx2z+ABY*I_ESP_G4z_Fx2z;
    Double I_ESP_G4x_Gx3z = I_ESP_H5x_F3z+ABX*I_ESP_G4x_F3z;
    Double I_ESP_G3xy_Gx3z = I_ESP_H4xy_F3z+ABX*I_ESP_G3xy_F3z;
    Double I_ESP_G3xz_Gx3z = I_ESP_H4xz_F3z+ABX*I_ESP_G3xz_F3z;
    Double I_ESP_G2x2y_Gx3z = I_ESP_H3x2y_F3z+ABX*I_ESP_G2x2y_F3z;
    Double I_ESP_G2xyz_Gx3z = I_ESP_H3xyz_F3z+ABX*I_ESP_G2xyz_F3z;
    Double I_ESP_G2x2z_Gx3z = I_ESP_H3x2z_F3z+ABX*I_ESP_G2x2z_F3z;
    Double I_ESP_Gx3y_Gx3z = I_ESP_H2x3y_F3z+ABX*I_ESP_Gx3y_F3z;
    Double I_ESP_Gx2yz_Gx3z = I_ESP_H2x2yz_F3z+ABX*I_ESP_Gx2yz_F3z;
    Double I_ESP_Gxy2z_Gx3z = I_ESP_H2xy2z_F3z+ABX*I_ESP_Gxy2z_F3z;
    Double I_ESP_Gx3z_Gx3z = I_ESP_H2x3z_F3z+ABX*I_ESP_Gx3z_F3z;
    Double I_ESP_G4y_Gx3z = I_ESP_Hx4y_F3z+ABX*I_ESP_G4y_F3z;
    Double I_ESP_G3yz_Gx3z = I_ESP_Hx3yz_F3z+ABX*I_ESP_G3yz_F3z;
    Double I_ESP_G2y2z_Gx3z = I_ESP_Hx2y2z_F3z+ABX*I_ESP_G2y2z_F3z;
    Double I_ESP_Gy3z_Gx3z = I_ESP_Hxy3z_F3z+ABX*I_ESP_Gy3z_F3z;
    Double I_ESP_G4z_Gx3z = I_ESP_Hx4z_F3z+ABX*I_ESP_G4z_F3z;
    Double I_ESP_G4x_G4y = I_ESP_H4xy_F3y+ABY*I_ESP_G4x_F3y;
    Double I_ESP_G3xy_G4y = I_ESP_H3x2y_F3y+ABY*I_ESP_G3xy_F3y;
    Double I_ESP_G3xz_G4y = I_ESP_H3xyz_F3y+ABY*I_ESP_G3xz_F3y;
    Double I_ESP_G2x2y_G4y = I_ESP_H2x3y_F3y+ABY*I_ESP_G2x2y_F3y;
    Double I_ESP_G2xyz_G4y = I_ESP_H2x2yz_F3y+ABY*I_ESP_G2xyz_F3y;
    Double I_ESP_G2x2z_G4y = I_ESP_H2xy2z_F3y+ABY*I_ESP_G2x2z_F3y;
    Double I_ESP_Gx3y_G4y = I_ESP_Hx4y_F3y+ABY*I_ESP_Gx3y_F3y;
    Double I_ESP_Gx2yz_G4y = I_ESP_Hx3yz_F3y+ABY*I_ESP_Gx2yz_F3y;
    Double I_ESP_Gxy2z_G4y = I_ESP_Hx2y2z_F3y+ABY*I_ESP_Gxy2z_F3y;
    Double I_ESP_Gx3z_G4y = I_ESP_Hxy3z_F3y+ABY*I_ESP_Gx3z_F3y;
    Double I_ESP_G4y_G4y = I_ESP_H5y_F3y+ABY*I_ESP_G4y_F3y;
    Double I_ESP_G3yz_G4y = I_ESP_H4yz_F3y+ABY*I_ESP_G3yz_F3y;
    Double I_ESP_G2y2z_G4y = I_ESP_H3y2z_F3y+ABY*I_ESP_G2y2z_F3y;
    Double I_ESP_Gy3z_G4y = I_ESP_H2y3z_F3y+ABY*I_ESP_Gy3z_F3y;
    Double I_ESP_G4z_G4y = I_ESP_Hy4z_F3y+ABY*I_ESP_G4z_F3y;
    Double I_ESP_G4x_G3yz = I_ESP_H4xz_F3y+ABZ*I_ESP_G4x_F3y;
    Double I_ESP_G3xy_G3yz = I_ESP_H3xyz_F3y+ABZ*I_ESP_G3xy_F3y;
    Double I_ESP_G3xz_G3yz = I_ESP_H3x2z_F3y+ABZ*I_ESP_G3xz_F3y;
    Double I_ESP_G2x2y_G3yz = I_ESP_H2x2yz_F3y+ABZ*I_ESP_G2x2y_F3y;
    Double I_ESP_G2xyz_G3yz = I_ESP_H2xy2z_F3y+ABZ*I_ESP_G2xyz_F3y;
    Double I_ESP_G2x2z_G3yz = I_ESP_H2x3z_F3y+ABZ*I_ESP_G2x2z_F3y;
    Double I_ESP_Gx3y_G3yz = I_ESP_Hx3yz_F3y+ABZ*I_ESP_Gx3y_F3y;
    Double I_ESP_Gx2yz_G3yz = I_ESP_Hx2y2z_F3y+ABZ*I_ESP_Gx2yz_F3y;
    Double I_ESP_Gxy2z_G3yz = I_ESP_Hxy3z_F3y+ABZ*I_ESP_Gxy2z_F3y;
    Double I_ESP_Gx3z_G3yz = I_ESP_Hx4z_F3y+ABZ*I_ESP_Gx3z_F3y;
    Double I_ESP_G4y_G3yz = I_ESP_H4yz_F3y+ABZ*I_ESP_G4y_F3y;
    Double I_ESP_G3yz_G3yz = I_ESP_H3y2z_F3y+ABZ*I_ESP_G3yz_F3y;
    Double I_ESP_G2y2z_G3yz = I_ESP_H2y3z_F3y+ABZ*I_ESP_G2y2z_F3y;
    Double I_ESP_Gy3z_G3yz = I_ESP_Hy4z_F3y+ABZ*I_ESP_Gy3z_F3y;
    Double I_ESP_G4z_G3yz = I_ESP_H5z_F3y+ABZ*I_ESP_G4z_F3y;
    Double I_ESP_G4x_G2y2z = I_ESP_H4xz_F2yz+ABZ*I_ESP_G4x_F2yz;
    Double I_ESP_G3xy_G2y2z = I_ESP_H3xyz_F2yz+ABZ*I_ESP_G3xy_F2yz;
    Double I_ESP_G3xz_G2y2z = I_ESP_H3x2z_F2yz+ABZ*I_ESP_G3xz_F2yz;
    Double I_ESP_G2x2y_G2y2z = I_ESP_H2x2yz_F2yz+ABZ*I_ESP_G2x2y_F2yz;
    Double I_ESP_G2xyz_G2y2z = I_ESP_H2xy2z_F2yz+ABZ*I_ESP_G2xyz_F2yz;
    Double I_ESP_G2x2z_G2y2z = I_ESP_H2x3z_F2yz+ABZ*I_ESP_G2x2z_F2yz;
    Double I_ESP_Gx3y_G2y2z = I_ESP_Hx3yz_F2yz+ABZ*I_ESP_Gx3y_F2yz;
    Double I_ESP_Gx2yz_G2y2z = I_ESP_Hx2y2z_F2yz+ABZ*I_ESP_Gx2yz_F2yz;
    Double I_ESP_Gxy2z_G2y2z = I_ESP_Hxy3z_F2yz+ABZ*I_ESP_Gxy2z_F2yz;
    Double I_ESP_Gx3z_G2y2z = I_ESP_Hx4z_F2yz+ABZ*I_ESP_Gx3z_F2yz;
    Double I_ESP_G4y_G2y2z = I_ESP_H4yz_F2yz+ABZ*I_ESP_G4y_F2yz;
    Double I_ESP_G3yz_G2y2z = I_ESP_H3y2z_F2yz+ABZ*I_ESP_G3yz_F2yz;
    Double I_ESP_G2y2z_G2y2z = I_ESP_H2y3z_F2yz+ABZ*I_ESP_G2y2z_F2yz;
    Double I_ESP_Gy3z_G2y2z = I_ESP_Hy4z_F2yz+ABZ*I_ESP_Gy3z_F2yz;
    Double I_ESP_G4z_G2y2z = I_ESP_H5z_F2yz+ABZ*I_ESP_G4z_F2yz;
    Double I_ESP_G4x_Gy3z = I_ESP_H4xy_F3z+ABY*I_ESP_G4x_F3z;
    Double I_ESP_G3xy_Gy3z = I_ESP_H3x2y_F3z+ABY*I_ESP_G3xy_F3z;
    Double I_ESP_G3xz_Gy3z = I_ESP_H3xyz_F3z+ABY*I_ESP_G3xz_F3z;
    Double I_ESP_G2x2y_Gy3z = I_ESP_H2x3y_F3z+ABY*I_ESP_G2x2y_F3z;
    Double I_ESP_G2xyz_Gy3z = I_ESP_H2x2yz_F3z+ABY*I_ESP_G2xyz_F3z;
    Double I_ESP_G2x2z_Gy3z = I_ESP_H2xy2z_F3z+ABY*I_ESP_G2x2z_F3z;
    Double I_ESP_Gx3y_Gy3z = I_ESP_Hx4y_F3z+ABY*I_ESP_Gx3y_F3z;
    Double I_ESP_Gx2yz_Gy3z = I_ESP_Hx3yz_F3z+ABY*I_ESP_Gx2yz_F3z;
    Double I_ESP_Gxy2z_Gy3z = I_ESP_Hx2y2z_F3z+ABY*I_ESP_Gxy2z_F3z;
    Double I_ESP_Gx3z_Gy3z = I_ESP_Hxy3z_F3z+ABY*I_ESP_Gx3z_F3z;
    Double I_ESP_G4y_Gy3z = I_ESP_H5y_F3z+ABY*I_ESP_G4y_F3z;
    Double I_ESP_G3yz_Gy3z = I_ESP_H4yz_F3z+ABY*I_ESP_G3yz_F3z;
    Double I_ESP_G2y2z_Gy3z = I_ESP_H3y2z_F3z+ABY*I_ESP_G2y2z_F3z;
    Double I_ESP_Gy3z_Gy3z = I_ESP_H2y3z_F3z+ABY*I_ESP_Gy3z_F3z;
    Double I_ESP_G4z_Gy3z = I_ESP_Hy4z_F3z+ABY*I_ESP_G4z_F3z;
    Double I_ESP_G4x_G4z = I_ESP_H4xz_F3z+ABZ*I_ESP_G4x_F3z;
    Double I_ESP_G3xy_G4z = I_ESP_H3xyz_F3z+ABZ*I_ESP_G3xy_F3z;
    Double I_ESP_G3xz_G4z = I_ESP_H3x2z_F3z+ABZ*I_ESP_G3xz_F3z;
    Double I_ESP_G2x2y_G4z = I_ESP_H2x2yz_F3z+ABZ*I_ESP_G2x2y_F3z;
    Double I_ESP_G2xyz_G4z = I_ESP_H2xy2z_F3z+ABZ*I_ESP_G2xyz_F3z;
    Double I_ESP_G2x2z_G4z = I_ESP_H2x3z_F3z+ABZ*I_ESP_G2x2z_F3z;
    Double I_ESP_Gx3y_G4z = I_ESP_Hx3yz_F3z+ABZ*I_ESP_Gx3y_F3z;
    Double I_ESP_Gx2yz_G4z = I_ESP_Hx2y2z_F3z+ABZ*I_ESP_Gx2yz_F3z;
    Double I_ESP_Gxy2z_G4z = I_ESP_Hxy3z_F3z+ABZ*I_ESP_Gxy2z_F3z;
    Double I_ESP_Gx3z_G4z = I_ESP_Hx4z_F3z+ABZ*I_ESP_Gx3z_F3z;
    Double I_ESP_G4y_G4z = I_ESP_H4yz_F3z+ABZ*I_ESP_G4y_F3z;
    Double I_ESP_G3yz_G4z = I_ESP_H3y2z_F3z+ABZ*I_ESP_G3yz_F3z;
    Double I_ESP_G2y2z_G4z = I_ESP_H2y3z_F3z+ABZ*I_ESP_G2y2z_F3z;
    Double I_ESP_Gy3z_G4z = I_ESP_Hy4z_F3z+ABZ*I_ESP_Gy3z_F3z;
    Double I_ESP_G4z_G4z = I_ESP_H5z_F3z+ABZ*I_ESP_G4z_F3z;

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
     * shell quartet name: SQ_ESP_K_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
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
    Double I_ESP_K7x_Py_a = I_ESP_L7xy_S_a+ABY*I_ESP_K7x_S_a;
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
    Double I_ESP_K7x_Pz_a = I_ESP_L7xz_S_a+ABZ*I_ESP_K7x_S_a;
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
    Double I_ESP_K7y_Pz_a = I_ESP_L7yz_S_a+ABZ*I_ESP_K7y_S_a;
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
     * shell quartet name: SQ_ESP_L_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 3 integrals are omitted 
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
    Double I_ESP_L8y_Px_a = I_ESP_Mx8y_S_a+ABX*I_ESP_L8y_S_a;
    Double I_ESP_L7yz_Px_a = I_ESP_Mx7yz_S_a+ABX*I_ESP_L7yz_S_a;
    Double I_ESP_L6y2z_Px_a = I_ESP_Mx6y2z_S_a+ABX*I_ESP_L6y2z_S_a;
    Double I_ESP_L5y3z_Px_a = I_ESP_Mx5y3z_S_a+ABX*I_ESP_L5y3z_S_a;
    Double I_ESP_L4y4z_Px_a = I_ESP_Mx4y4z_S_a+ABX*I_ESP_L4y4z_S_a;
    Double I_ESP_L3y5z_Px_a = I_ESP_Mx3y5z_S_a+ABX*I_ESP_L3y5z_S_a;
    Double I_ESP_L2y6z_Px_a = I_ESP_Mx2y6z_S_a+ABX*I_ESP_L2y6z_S_a;
    Double I_ESP_Ly7z_Px_a = I_ESP_Mxy7z_S_a+ABX*I_ESP_Ly7z_S_a;
    Double I_ESP_L8z_Px_a = I_ESP_Mx8z_S_a+ABX*I_ESP_L8z_S_a;
    Double I_ESP_L7xy_Py_a = I_ESP_M7x2y_S_a+ABY*I_ESP_L7xy_S_a;
    Double I_ESP_L7xz_Py_a = I_ESP_M7xyz_S_a+ABY*I_ESP_L7xz_S_a;
    Double I_ESP_L6x2y_Py_a = I_ESP_M6x3y_S_a+ABY*I_ESP_L6x2y_S_a;
    Double I_ESP_L6xyz_Py_a = I_ESP_M6x2yz_S_a+ABY*I_ESP_L6xyz_S_a;
    Double I_ESP_L6x2z_Py_a = I_ESP_M6xy2z_S_a+ABY*I_ESP_L6x2z_S_a;
    Double I_ESP_L5x3y_Py_a = I_ESP_M5x4y_S_a+ABY*I_ESP_L5x3y_S_a;
    Double I_ESP_L5x2yz_Py_a = I_ESP_M5x3yz_S_a+ABY*I_ESP_L5x2yz_S_a;
    Double I_ESP_L5xy2z_Py_a = I_ESP_M5x2y2z_S_a+ABY*I_ESP_L5xy2z_S_a;
    Double I_ESP_L5x3z_Py_a = I_ESP_M5xy3z_S_a+ABY*I_ESP_L5x3z_S_a;
    Double I_ESP_L4x4y_Py_a = I_ESP_M4x5y_S_a+ABY*I_ESP_L4x4y_S_a;
    Double I_ESP_L4x3yz_Py_a = I_ESP_M4x4yz_S_a+ABY*I_ESP_L4x3yz_S_a;
    Double I_ESP_L4x2y2z_Py_a = I_ESP_M4x3y2z_S_a+ABY*I_ESP_L4x2y2z_S_a;
    Double I_ESP_L4xy3z_Py_a = I_ESP_M4x2y3z_S_a+ABY*I_ESP_L4xy3z_S_a;
    Double I_ESP_L4x4z_Py_a = I_ESP_M4xy4z_S_a+ABY*I_ESP_L4x4z_S_a;
    Double I_ESP_L3x5y_Py_a = I_ESP_M3x6y_S_a+ABY*I_ESP_L3x5y_S_a;
    Double I_ESP_L3x4yz_Py_a = I_ESP_M3x5yz_S_a+ABY*I_ESP_L3x4yz_S_a;
    Double I_ESP_L3x3y2z_Py_a = I_ESP_M3x4y2z_S_a+ABY*I_ESP_L3x3y2z_S_a;
    Double I_ESP_L3x2y3z_Py_a = I_ESP_M3x3y3z_S_a+ABY*I_ESP_L3x2y3z_S_a;
    Double I_ESP_L3xy4z_Py_a = I_ESP_M3x2y4z_S_a+ABY*I_ESP_L3xy4z_S_a;
    Double I_ESP_L3x5z_Py_a = I_ESP_M3xy5z_S_a+ABY*I_ESP_L3x5z_S_a;
    Double I_ESP_L2x6y_Py_a = I_ESP_M2x7y_S_a+ABY*I_ESP_L2x6y_S_a;
    Double I_ESP_L2x5yz_Py_a = I_ESP_M2x6yz_S_a+ABY*I_ESP_L2x5yz_S_a;
    Double I_ESP_L2x4y2z_Py_a = I_ESP_M2x5y2z_S_a+ABY*I_ESP_L2x4y2z_S_a;
    Double I_ESP_L2x3y3z_Py_a = I_ESP_M2x4y3z_S_a+ABY*I_ESP_L2x3y3z_S_a;
    Double I_ESP_L2x2y4z_Py_a = I_ESP_M2x3y4z_S_a+ABY*I_ESP_L2x2y4z_S_a;
    Double I_ESP_L2xy5z_Py_a = I_ESP_M2x2y5z_S_a+ABY*I_ESP_L2xy5z_S_a;
    Double I_ESP_L2x6z_Py_a = I_ESP_M2xy6z_S_a+ABY*I_ESP_L2x6z_S_a;
    Double I_ESP_Lx7y_Py_a = I_ESP_Mx8y_S_a+ABY*I_ESP_Lx7y_S_a;
    Double I_ESP_Lx6yz_Py_a = I_ESP_Mx7yz_S_a+ABY*I_ESP_Lx6yz_S_a;
    Double I_ESP_Lx5y2z_Py_a = I_ESP_Mx6y2z_S_a+ABY*I_ESP_Lx5y2z_S_a;
    Double I_ESP_Lx4y3z_Py_a = I_ESP_Mx5y3z_S_a+ABY*I_ESP_Lx4y3z_S_a;
    Double I_ESP_Lx3y4z_Py_a = I_ESP_Mx4y4z_S_a+ABY*I_ESP_Lx3y4z_S_a;
    Double I_ESP_Lx2y5z_Py_a = I_ESP_Mx3y5z_S_a+ABY*I_ESP_Lx2y5z_S_a;
    Double I_ESP_Lxy6z_Py_a = I_ESP_Mx2y6z_S_a+ABY*I_ESP_Lxy6z_S_a;
    Double I_ESP_Lx7z_Py_a = I_ESP_Mxy7z_S_a+ABY*I_ESP_Lx7z_S_a;
    Double I_ESP_L8y_Py_a = I_ESP_M9y_S_a+ABY*I_ESP_L8y_S_a;
    Double I_ESP_L7yz_Py_a = I_ESP_M8yz_S_a+ABY*I_ESP_L7yz_S_a;
    Double I_ESP_L6y2z_Py_a = I_ESP_M7y2z_S_a+ABY*I_ESP_L6y2z_S_a;
    Double I_ESP_L5y3z_Py_a = I_ESP_M6y3z_S_a+ABY*I_ESP_L5y3z_S_a;
    Double I_ESP_L4y4z_Py_a = I_ESP_M5y4z_S_a+ABY*I_ESP_L4y4z_S_a;
    Double I_ESP_L3y5z_Py_a = I_ESP_M4y5z_S_a+ABY*I_ESP_L3y5z_S_a;
    Double I_ESP_L2y6z_Py_a = I_ESP_M3y6z_S_a+ABY*I_ESP_L2y6z_S_a;
    Double I_ESP_Ly7z_Py_a = I_ESP_M2y7z_S_a+ABY*I_ESP_Ly7z_S_a;
    Double I_ESP_L8z_Py_a = I_ESP_My8z_S_a+ABY*I_ESP_L8z_S_a;
    Double I_ESP_L7xy_Pz_a = I_ESP_M7xyz_S_a+ABZ*I_ESP_L7xy_S_a;
    Double I_ESP_L7xz_Pz_a = I_ESP_M7x2z_S_a+ABZ*I_ESP_L7xz_S_a;
    Double I_ESP_L6x2y_Pz_a = I_ESP_M6x2yz_S_a+ABZ*I_ESP_L6x2y_S_a;
    Double I_ESP_L6xyz_Pz_a = I_ESP_M6xy2z_S_a+ABZ*I_ESP_L6xyz_S_a;
    Double I_ESP_L6x2z_Pz_a = I_ESP_M6x3z_S_a+ABZ*I_ESP_L6x2z_S_a;
    Double I_ESP_L5x3y_Pz_a = I_ESP_M5x3yz_S_a+ABZ*I_ESP_L5x3y_S_a;
    Double I_ESP_L5x2yz_Pz_a = I_ESP_M5x2y2z_S_a+ABZ*I_ESP_L5x2yz_S_a;
    Double I_ESP_L5xy2z_Pz_a = I_ESP_M5xy3z_S_a+ABZ*I_ESP_L5xy2z_S_a;
    Double I_ESP_L5x3z_Pz_a = I_ESP_M5x4z_S_a+ABZ*I_ESP_L5x3z_S_a;
    Double I_ESP_L4x4y_Pz_a = I_ESP_M4x4yz_S_a+ABZ*I_ESP_L4x4y_S_a;
    Double I_ESP_L4x3yz_Pz_a = I_ESP_M4x3y2z_S_a+ABZ*I_ESP_L4x3yz_S_a;
    Double I_ESP_L4x2y2z_Pz_a = I_ESP_M4x2y3z_S_a+ABZ*I_ESP_L4x2y2z_S_a;
    Double I_ESP_L4xy3z_Pz_a = I_ESP_M4xy4z_S_a+ABZ*I_ESP_L4xy3z_S_a;
    Double I_ESP_L4x4z_Pz_a = I_ESP_M4x5z_S_a+ABZ*I_ESP_L4x4z_S_a;
    Double I_ESP_L3x5y_Pz_a = I_ESP_M3x5yz_S_a+ABZ*I_ESP_L3x5y_S_a;
    Double I_ESP_L3x4yz_Pz_a = I_ESP_M3x4y2z_S_a+ABZ*I_ESP_L3x4yz_S_a;
    Double I_ESP_L3x3y2z_Pz_a = I_ESP_M3x3y3z_S_a+ABZ*I_ESP_L3x3y2z_S_a;
    Double I_ESP_L3x2y3z_Pz_a = I_ESP_M3x2y4z_S_a+ABZ*I_ESP_L3x2y3z_S_a;
    Double I_ESP_L3xy4z_Pz_a = I_ESP_M3xy5z_S_a+ABZ*I_ESP_L3xy4z_S_a;
    Double I_ESP_L3x5z_Pz_a = I_ESP_M3x6z_S_a+ABZ*I_ESP_L3x5z_S_a;
    Double I_ESP_L2x6y_Pz_a = I_ESP_M2x6yz_S_a+ABZ*I_ESP_L2x6y_S_a;
    Double I_ESP_L2x5yz_Pz_a = I_ESP_M2x5y2z_S_a+ABZ*I_ESP_L2x5yz_S_a;
    Double I_ESP_L2x4y2z_Pz_a = I_ESP_M2x4y3z_S_a+ABZ*I_ESP_L2x4y2z_S_a;
    Double I_ESP_L2x3y3z_Pz_a = I_ESP_M2x3y4z_S_a+ABZ*I_ESP_L2x3y3z_S_a;
    Double I_ESP_L2x2y4z_Pz_a = I_ESP_M2x2y5z_S_a+ABZ*I_ESP_L2x2y4z_S_a;
    Double I_ESP_L2xy5z_Pz_a = I_ESP_M2xy6z_S_a+ABZ*I_ESP_L2xy5z_S_a;
    Double I_ESP_L2x6z_Pz_a = I_ESP_M2x7z_S_a+ABZ*I_ESP_L2x6z_S_a;
    Double I_ESP_Lx7y_Pz_a = I_ESP_Mx7yz_S_a+ABZ*I_ESP_Lx7y_S_a;
    Double I_ESP_Lx6yz_Pz_a = I_ESP_Mx6y2z_S_a+ABZ*I_ESP_Lx6yz_S_a;
    Double I_ESP_Lx5y2z_Pz_a = I_ESP_Mx5y3z_S_a+ABZ*I_ESP_Lx5y2z_S_a;
    Double I_ESP_Lx4y3z_Pz_a = I_ESP_Mx4y4z_S_a+ABZ*I_ESP_Lx4y3z_S_a;
    Double I_ESP_Lx3y4z_Pz_a = I_ESP_Mx3y5z_S_a+ABZ*I_ESP_Lx3y4z_S_a;
    Double I_ESP_Lx2y5z_Pz_a = I_ESP_Mx2y6z_S_a+ABZ*I_ESP_Lx2y5z_S_a;
    Double I_ESP_Lxy6z_Pz_a = I_ESP_Mxy7z_S_a+ABZ*I_ESP_Lxy6z_S_a;
    Double I_ESP_Lx7z_Pz_a = I_ESP_Mx8z_S_a+ABZ*I_ESP_Lx7z_S_a;
    Double I_ESP_L7yz_Pz_a = I_ESP_M7y2z_S_a+ABZ*I_ESP_L7yz_S_a;
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
     * totally 108 integrals are omitted 
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
    Double I_ESP_K7x_D2y_a = I_ESP_L7xy_Py_a+ABY*I_ESP_K7x_Py_a;
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
    Double I_ESP_K7x_D2z_a = I_ESP_L7xz_Pz_a+ABZ*I_ESP_K7x_Pz_a;
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
    Double I_ESP_K7y_D2z_a = I_ESP_L7yz_Pz_a+ABZ*I_ESP_K7y_Pz_a;
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
     * totally 56 integrals are omitted 
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
    Double I_ESP_I6x_F2xy_a = I_ESP_K6xy_D2x_a+ABY*I_ESP_I6x_D2x_a;
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
    Double I_ESP_I6x_F2xz_a = I_ESP_K6xz_D2x_a+ABZ*I_ESP_I6x_D2x_a;
    Double I_ESP_I5xy_F2xz_a = I_ESP_K5xyz_D2x_a+ABZ*I_ESP_I5xy_D2x_a;
    Double I_ESP_I5xz_F2xz_a = I_ESP_K5x2z_D2x_a+ABZ*I_ESP_I5xz_D2x_a;
    Double I_ESP_I4x2y_F2xz_a = I_ESP_K4x2yz_D2x_a+ABZ*I_ESP_I4x2y_D2x_a;
    Double I_ESP_I4xyz_F2xz_a = I_ESP_K4xy2z_D2x_a+ABZ*I_ESP_I4xyz_D2x_a;
    Double I_ESP_I4x2z_F2xz_a = I_ESP_K4x3z_D2x_a+ABZ*I_ESP_I4x2z_D2x_a;
    Double I_ESP_I3x3y_F2xz_a = I_ESP_K3x3yz_D2x_a+ABZ*I_ESP_I3x3y_D2x_a;
    Double I_ESP_I3x2yz_F2xz_a = I_ESP_K3x2y2z_D2x_a+ABZ*I_ESP_I3x2yz_D2x_a;
    Double I_ESP_I3xy2z_F2xz_a = I_ESP_K3xy3z_D2x_a+ABZ*I_ESP_I3xy2z_D2x_a;
    Double I_ESP_I3x3z_F2xz_a = I_ESP_K3x4z_D2x_a+ABZ*I_ESP_I3x3z_D2x_a;
    Double I_ESP_I2x4y_F2xz_a = I_ESP_K2x4yz_D2x_a+ABZ*I_ESP_I2x4y_D2x_a;
    Double I_ESP_I2x3yz_F2xz_a = I_ESP_K2x3y2z_D2x_a+ABZ*I_ESP_I2x3yz_D2x_a;
    Double I_ESP_I2x2y2z_F2xz_a = I_ESP_K2x2y3z_D2x_a+ABZ*I_ESP_I2x2y2z_D2x_a;
    Double I_ESP_I2xy3z_F2xz_a = I_ESP_K2xy4z_D2x_a+ABZ*I_ESP_I2xy3z_D2x_a;
    Double I_ESP_I2x4z_F2xz_a = I_ESP_K2x5z_D2x_a+ABZ*I_ESP_I2x4z_D2x_a;
    Double I_ESP_Ix5y_F2xz_a = I_ESP_Kx5yz_D2x_a+ABZ*I_ESP_Ix5y_D2x_a;
    Double I_ESP_Ix4yz_F2xz_a = I_ESP_Kx4y2z_D2x_a+ABZ*I_ESP_Ix4yz_D2x_a;
    Double I_ESP_Ix3y2z_F2xz_a = I_ESP_Kx3y3z_D2x_a+ABZ*I_ESP_Ix3y2z_D2x_a;
    Double I_ESP_Ix2y3z_F2xz_a = I_ESP_Kx2y4z_D2x_a+ABZ*I_ESP_Ix2y3z_D2x_a;
    Double I_ESP_Ixy4z_F2xz_a = I_ESP_Kxy5z_D2x_a+ABZ*I_ESP_Ixy4z_D2x_a;
    Double I_ESP_Ix5z_F2xz_a = I_ESP_Kx6z_D2x_a+ABZ*I_ESP_Ix5z_D2x_a;
    Double I_ESP_I6y_F2xz_a = I_ESP_K6yz_D2x_a+ABZ*I_ESP_I6y_D2x_a;
    Double I_ESP_I5yz_F2xz_a = I_ESP_K5y2z_D2x_a+ABZ*I_ESP_I5yz_D2x_a;
    Double I_ESP_I4y2z_F2xz_a = I_ESP_K4y3z_D2x_a+ABZ*I_ESP_I4y2z_D2x_a;
    Double I_ESP_I3y3z_F2xz_a = I_ESP_K3y4z_D2x_a+ABZ*I_ESP_I3y3z_D2x_a;
    Double I_ESP_I2y4z_F2xz_a = I_ESP_K2y5z_D2x_a+ABZ*I_ESP_I2y4z_D2x_a;
    Double I_ESP_Iy5z_F2xz_a = I_ESP_Ky6z_D2x_a+ABZ*I_ESP_Iy5z_D2x_a;
    Double I_ESP_I6z_F2xz_a = I_ESP_K7z_D2x_a+ABZ*I_ESP_I6z_D2x_a;
    Double I_ESP_I6x_Fx2y_a = I_ESP_K7x_D2y_a+ABX*I_ESP_I6x_D2y_a;
    Double I_ESP_I5xy_Fx2y_a = I_ESP_K6xy_D2y_a+ABX*I_ESP_I5xy_D2y_a;
    Double I_ESP_I5xz_Fx2y_a = I_ESP_K6xz_D2y_a+ABX*I_ESP_I5xz_D2y_a;
    Double I_ESP_I4x2y_Fx2y_a = I_ESP_K5x2y_D2y_a+ABX*I_ESP_I4x2y_D2y_a;
    Double I_ESP_I4xyz_Fx2y_a = I_ESP_K5xyz_D2y_a+ABX*I_ESP_I4xyz_D2y_a;
    Double I_ESP_I4x2z_Fx2y_a = I_ESP_K5x2z_D2y_a+ABX*I_ESP_I4x2z_D2y_a;
    Double I_ESP_I3x3y_Fx2y_a = I_ESP_K4x3y_D2y_a+ABX*I_ESP_I3x3y_D2y_a;
    Double I_ESP_I3x2yz_Fx2y_a = I_ESP_K4x2yz_D2y_a+ABX*I_ESP_I3x2yz_D2y_a;
    Double I_ESP_I3xy2z_Fx2y_a = I_ESP_K4xy2z_D2y_a+ABX*I_ESP_I3xy2z_D2y_a;
    Double I_ESP_I3x3z_Fx2y_a = I_ESP_K4x3z_D2y_a+ABX*I_ESP_I3x3z_D2y_a;
    Double I_ESP_I2x4y_Fx2y_a = I_ESP_K3x4y_D2y_a+ABX*I_ESP_I2x4y_D2y_a;
    Double I_ESP_I2x3yz_Fx2y_a = I_ESP_K3x3yz_D2y_a+ABX*I_ESP_I2x3yz_D2y_a;
    Double I_ESP_I2x2y2z_Fx2y_a = I_ESP_K3x2y2z_D2y_a+ABX*I_ESP_I2x2y2z_D2y_a;
    Double I_ESP_I2xy3z_Fx2y_a = I_ESP_K3xy3z_D2y_a+ABX*I_ESP_I2xy3z_D2y_a;
    Double I_ESP_I2x4z_Fx2y_a = I_ESP_K3x4z_D2y_a+ABX*I_ESP_I2x4z_D2y_a;
    Double I_ESP_Ix5y_Fx2y_a = I_ESP_K2x5y_D2y_a+ABX*I_ESP_Ix5y_D2y_a;
    Double I_ESP_Ix4yz_Fx2y_a = I_ESP_K2x4yz_D2y_a+ABX*I_ESP_Ix4yz_D2y_a;
    Double I_ESP_Ix3y2z_Fx2y_a = I_ESP_K2x3y2z_D2y_a+ABX*I_ESP_Ix3y2z_D2y_a;
    Double I_ESP_Ix2y3z_Fx2y_a = I_ESP_K2x2y3z_D2y_a+ABX*I_ESP_Ix2y3z_D2y_a;
    Double I_ESP_Ixy4z_Fx2y_a = I_ESP_K2xy4z_D2y_a+ABX*I_ESP_Ixy4z_D2y_a;
    Double I_ESP_Ix5z_Fx2y_a = I_ESP_K2x5z_D2y_a+ABX*I_ESP_Ix5z_D2y_a;
    Double I_ESP_I6y_Fx2y_a = I_ESP_Kx6y_D2y_a+ABX*I_ESP_I6y_D2y_a;
    Double I_ESP_I5yz_Fx2y_a = I_ESP_Kx5yz_D2y_a+ABX*I_ESP_I5yz_D2y_a;
    Double I_ESP_I4y2z_Fx2y_a = I_ESP_Kx4y2z_D2y_a+ABX*I_ESP_I4y2z_D2y_a;
    Double I_ESP_I3y3z_Fx2y_a = I_ESP_Kx3y3z_D2y_a+ABX*I_ESP_I3y3z_D2y_a;
    Double I_ESP_I2y4z_Fx2y_a = I_ESP_Kx2y4z_D2y_a+ABX*I_ESP_I2y4z_D2y_a;
    Double I_ESP_Iy5z_Fx2y_a = I_ESP_Kxy5z_D2y_a+ABX*I_ESP_Iy5z_D2y_a;
    Double I_ESP_I6z_Fx2y_a = I_ESP_Kx6z_D2y_a+ABX*I_ESP_I6z_D2y_a;
    Double I_ESP_I6x_Fx2z_a = I_ESP_K7x_D2z_a+ABX*I_ESP_I6x_D2z_a;
    Double I_ESP_I5xy_Fx2z_a = I_ESP_K6xy_D2z_a+ABX*I_ESP_I5xy_D2z_a;
    Double I_ESP_I5xz_Fx2z_a = I_ESP_K6xz_D2z_a+ABX*I_ESP_I5xz_D2z_a;
    Double I_ESP_I4x2y_Fx2z_a = I_ESP_K5x2y_D2z_a+ABX*I_ESP_I4x2y_D2z_a;
    Double I_ESP_I4xyz_Fx2z_a = I_ESP_K5xyz_D2z_a+ABX*I_ESP_I4xyz_D2z_a;
    Double I_ESP_I4x2z_Fx2z_a = I_ESP_K5x2z_D2z_a+ABX*I_ESP_I4x2z_D2z_a;
    Double I_ESP_I3x3y_Fx2z_a = I_ESP_K4x3y_D2z_a+ABX*I_ESP_I3x3y_D2z_a;
    Double I_ESP_I3x2yz_Fx2z_a = I_ESP_K4x2yz_D2z_a+ABX*I_ESP_I3x2yz_D2z_a;
    Double I_ESP_I3xy2z_Fx2z_a = I_ESP_K4xy2z_D2z_a+ABX*I_ESP_I3xy2z_D2z_a;
    Double I_ESP_I3x3z_Fx2z_a = I_ESP_K4x3z_D2z_a+ABX*I_ESP_I3x3z_D2z_a;
    Double I_ESP_I2x4y_Fx2z_a = I_ESP_K3x4y_D2z_a+ABX*I_ESP_I2x4y_D2z_a;
    Double I_ESP_I2x3yz_Fx2z_a = I_ESP_K3x3yz_D2z_a+ABX*I_ESP_I2x3yz_D2z_a;
    Double I_ESP_I2x2y2z_Fx2z_a = I_ESP_K3x2y2z_D2z_a+ABX*I_ESP_I2x2y2z_D2z_a;
    Double I_ESP_I2xy3z_Fx2z_a = I_ESP_K3xy3z_D2z_a+ABX*I_ESP_I2xy3z_D2z_a;
    Double I_ESP_I2x4z_Fx2z_a = I_ESP_K3x4z_D2z_a+ABX*I_ESP_I2x4z_D2z_a;
    Double I_ESP_Ix5y_Fx2z_a = I_ESP_K2x5y_D2z_a+ABX*I_ESP_Ix5y_D2z_a;
    Double I_ESP_Ix4yz_Fx2z_a = I_ESP_K2x4yz_D2z_a+ABX*I_ESP_Ix4yz_D2z_a;
    Double I_ESP_Ix3y2z_Fx2z_a = I_ESP_K2x3y2z_D2z_a+ABX*I_ESP_Ix3y2z_D2z_a;
    Double I_ESP_Ix2y3z_Fx2z_a = I_ESP_K2x2y3z_D2z_a+ABX*I_ESP_Ix2y3z_D2z_a;
    Double I_ESP_Ixy4z_Fx2z_a = I_ESP_K2xy4z_D2z_a+ABX*I_ESP_Ixy4z_D2z_a;
    Double I_ESP_Ix5z_Fx2z_a = I_ESP_K2x5z_D2z_a+ABX*I_ESP_Ix5z_D2z_a;
    Double I_ESP_I6y_Fx2z_a = I_ESP_Kx6y_D2z_a+ABX*I_ESP_I6y_D2z_a;
    Double I_ESP_I5yz_Fx2z_a = I_ESP_Kx5yz_D2z_a+ABX*I_ESP_I5yz_D2z_a;
    Double I_ESP_I4y2z_Fx2z_a = I_ESP_Kx4y2z_D2z_a+ABX*I_ESP_I4y2z_D2z_a;
    Double I_ESP_I3y3z_Fx2z_a = I_ESP_Kx3y3z_D2z_a+ABX*I_ESP_I3y3z_D2z_a;
    Double I_ESP_I2y4z_Fx2z_a = I_ESP_Kx2y4z_D2z_a+ABX*I_ESP_I2y4z_D2z_a;
    Double I_ESP_Iy5z_Fx2z_a = I_ESP_Kxy5z_D2z_a+ABX*I_ESP_Iy5z_D2z_a;
    Double I_ESP_I6z_Fx2z_a = I_ESP_Kx6z_D2z_a+ABX*I_ESP_I6z_D2z_a;
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
    Double I_ESP_I6x_F2yz_a = I_ESP_K6xz_D2y_a+ABZ*I_ESP_I6x_D2y_a;
    Double I_ESP_I5xy_F2yz_a = I_ESP_K5xyz_D2y_a+ABZ*I_ESP_I5xy_D2y_a;
    Double I_ESP_I5xz_F2yz_a = I_ESP_K5x2z_D2y_a+ABZ*I_ESP_I5xz_D2y_a;
    Double I_ESP_I4x2y_F2yz_a = I_ESP_K4x2yz_D2y_a+ABZ*I_ESP_I4x2y_D2y_a;
    Double I_ESP_I4xyz_F2yz_a = I_ESP_K4xy2z_D2y_a+ABZ*I_ESP_I4xyz_D2y_a;
    Double I_ESP_I4x2z_F2yz_a = I_ESP_K4x3z_D2y_a+ABZ*I_ESP_I4x2z_D2y_a;
    Double I_ESP_I3x3y_F2yz_a = I_ESP_K3x3yz_D2y_a+ABZ*I_ESP_I3x3y_D2y_a;
    Double I_ESP_I3x2yz_F2yz_a = I_ESP_K3x2y2z_D2y_a+ABZ*I_ESP_I3x2yz_D2y_a;
    Double I_ESP_I3xy2z_F2yz_a = I_ESP_K3xy3z_D2y_a+ABZ*I_ESP_I3xy2z_D2y_a;
    Double I_ESP_I3x3z_F2yz_a = I_ESP_K3x4z_D2y_a+ABZ*I_ESP_I3x3z_D2y_a;
    Double I_ESP_I2x4y_F2yz_a = I_ESP_K2x4yz_D2y_a+ABZ*I_ESP_I2x4y_D2y_a;
    Double I_ESP_I2x3yz_F2yz_a = I_ESP_K2x3y2z_D2y_a+ABZ*I_ESP_I2x3yz_D2y_a;
    Double I_ESP_I2x2y2z_F2yz_a = I_ESP_K2x2y3z_D2y_a+ABZ*I_ESP_I2x2y2z_D2y_a;
    Double I_ESP_I2xy3z_F2yz_a = I_ESP_K2xy4z_D2y_a+ABZ*I_ESP_I2xy3z_D2y_a;
    Double I_ESP_I2x4z_F2yz_a = I_ESP_K2x5z_D2y_a+ABZ*I_ESP_I2x4z_D2y_a;
    Double I_ESP_Ix5y_F2yz_a = I_ESP_Kx5yz_D2y_a+ABZ*I_ESP_Ix5y_D2y_a;
    Double I_ESP_Ix4yz_F2yz_a = I_ESP_Kx4y2z_D2y_a+ABZ*I_ESP_Ix4yz_D2y_a;
    Double I_ESP_Ix3y2z_F2yz_a = I_ESP_Kx3y3z_D2y_a+ABZ*I_ESP_Ix3y2z_D2y_a;
    Double I_ESP_Ix2y3z_F2yz_a = I_ESP_Kx2y4z_D2y_a+ABZ*I_ESP_Ix2y3z_D2y_a;
    Double I_ESP_Ixy4z_F2yz_a = I_ESP_Kxy5z_D2y_a+ABZ*I_ESP_Ixy4z_D2y_a;
    Double I_ESP_Ix5z_F2yz_a = I_ESP_Kx6z_D2y_a+ABZ*I_ESP_Ix5z_D2y_a;
    Double I_ESP_I6y_F2yz_a = I_ESP_K6yz_D2y_a+ABZ*I_ESP_I6y_D2y_a;
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
     * shell quartet name: SQ_ESP_M_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 33 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_N_S_a
     * RHS shell quartet name: SQ_ESP_M_S_a
     ************************************************************/
    Double I_ESP_M9x_Px_a = I_ESP_N10x_S_a+ABX*I_ESP_M9x_S_a;
    Double I_ESP_M8xy_Px_a = I_ESP_N9xy_S_a+ABX*I_ESP_M8xy_S_a;
    Double I_ESP_M8xz_Px_a = I_ESP_N9xz_S_a+ABX*I_ESP_M8xz_S_a;
    Double I_ESP_M7x2y_Px_a = I_ESP_N8x2y_S_a+ABX*I_ESP_M7x2y_S_a;
    Double I_ESP_M7xyz_Px_a = I_ESP_N8xyz_S_a+ABX*I_ESP_M7xyz_S_a;
    Double I_ESP_M7x2z_Px_a = I_ESP_N8x2z_S_a+ABX*I_ESP_M7x2z_S_a;
    Double I_ESP_M6x3y_Px_a = I_ESP_N7x3y_S_a+ABX*I_ESP_M6x3y_S_a;
    Double I_ESP_M6x2yz_Px_a = I_ESP_N7x2yz_S_a+ABX*I_ESP_M6x2yz_S_a;
    Double I_ESP_M6xy2z_Px_a = I_ESP_N7xy2z_S_a+ABX*I_ESP_M6xy2z_S_a;
    Double I_ESP_M6x3z_Px_a = I_ESP_N7x3z_S_a+ABX*I_ESP_M6x3z_S_a;
    Double I_ESP_M5x4y_Px_a = I_ESP_N6x4y_S_a+ABX*I_ESP_M5x4y_S_a;
    Double I_ESP_M5x3yz_Px_a = I_ESP_N6x3yz_S_a+ABX*I_ESP_M5x3yz_S_a;
    Double I_ESP_M5x2y2z_Px_a = I_ESP_N6x2y2z_S_a+ABX*I_ESP_M5x2y2z_S_a;
    Double I_ESP_M5xy3z_Px_a = I_ESP_N6xy3z_S_a+ABX*I_ESP_M5xy3z_S_a;
    Double I_ESP_M5x4z_Px_a = I_ESP_N6x4z_S_a+ABX*I_ESP_M5x4z_S_a;
    Double I_ESP_M4x5y_Px_a = I_ESP_N5x5y_S_a+ABX*I_ESP_M4x5y_S_a;
    Double I_ESP_M4x4yz_Px_a = I_ESP_N5x4yz_S_a+ABX*I_ESP_M4x4yz_S_a;
    Double I_ESP_M4x3y2z_Px_a = I_ESP_N5x3y2z_S_a+ABX*I_ESP_M4x3y2z_S_a;
    Double I_ESP_M4x2y3z_Px_a = I_ESP_N5x2y3z_S_a+ABX*I_ESP_M4x2y3z_S_a;
    Double I_ESP_M4xy4z_Px_a = I_ESP_N5xy4z_S_a+ABX*I_ESP_M4xy4z_S_a;
    Double I_ESP_M4x5z_Px_a = I_ESP_N5x5z_S_a+ABX*I_ESP_M4x5z_S_a;
    Double I_ESP_M3x6y_Px_a = I_ESP_N4x6y_S_a+ABX*I_ESP_M3x6y_S_a;
    Double I_ESP_M3x5yz_Px_a = I_ESP_N4x5yz_S_a+ABX*I_ESP_M3x5yz_S_a;
    Double I_ESP_M3x4y2z_Px_a = I_ESP_N4x4y2z_S_a+ABX*I_ESP_M3x4y2z_S_a;
    Double I_ESP_M3x3y3z_Px_a = I_ESP_N4x3y3z_S_a+ABX*I_ESP_M3x3y3z_S_a;
    Double I_ESP_M3x2y4z_Px_a = I_ESP_N4x2y4z_S_a+ABX*I_ESP_M3x2y4z_S_a;
    Double I_ESP_M3xy5z_Px_a = I_ESP_N4xy5z_S_a+ABX*I_ESP_M3xy5z_S_a;
    Double I_ESP_M3x6z_Px_a = I_ESP_N4x6z_S_a+ABX*I_ESP_M3x6z_S_a;
    Double I_ESP_M2x7y_Px_a = I_ESP_N3x7y_S_a+ABX*I_ESP_M2x7y_S_a;
    Double I_ESP_M2x6yz_Px_a = I_ESP_N3x6yz_S_a+ABX*I_ESP_M2x6yz_S_a;
    Double I_ESP_M2x5y2z_Px_a = I_ESP_N3x5y2z_S_a+ABX*I_ESP_M2x5y2z_S_a;
    Double I_ESP_M2x4y3z_Px_a = I_ESP_N3x4y3z_S_a+ABX*I_ESP_M2x4y3z_S_a;
    Double I_ESP_M2x3y4z_Px_a = I_ESP_N3x3y4z_S_a+ABX*I_ESP_M2x3y4z_S_a;
    Double I_ESP_M2x2y5z_Px_a = I_ESP_N3x2y5z_S_a+ABX*I_ESP_M2x2y5z_S_a;
    Double I_ESP_M2xy6z_Px_a = I_ESP_N3xy6z_S_a+ABX*I_ESP_M2xy6z_S_a;
    Double I_ESP_M2x7z_Px_a = I_ESP_N3x7z_S_a+ABX*I_ESP_M2x7z_S_a;
    Double I_ESP_Mx8y_Px_a = I_ESP_N2x8y_S_a+ABX*I_ESP_Mx8y_S_a;
    Double I_ESP_Mx7yz_Px_a = I_ESP_N2x7yz_S_a+ABX*I_ESP_Mx7yz_S_a;
    Double I_ESP_Mx6y2z_Px_a = I_ESP_N2x6y2z_S_a+ABX*I_ESP_Mx6y2z_S_a;
    Double I_ESP_Mx5y3z_Px_a = I_ESP_N2x5y3z_S_a+ABX*I_ESP_Mx5y3z_S_a;
    Double I_ESP_Mx4y4z_Px_a = I_ESP_N2x4y4z_S_a+ABX*I_ESP_Mx4y4z_S_a;
    Double I_ESP_Mx3y5z_Px_a = I_ESP_N2x3y5z_S_a+ABX*I_ESP_Mx3y5z_S_a;
    Double I_ESP_Mx2y6z_Px_a = I_ESP_N2x2y6z_S_a+ABX*I_ESP_Mx2y6z_S_a;
    Double I_ESP_Mxy7z_Px_a = I_ESP_N2xy7z_S_a+ABX*I_ESP_Mxy7z_S_a;
    Double I_ESP_Mx8z_Px_a = I_ESP_N2x8z_S_a+ABX*I_ESP_Mx8z_S_a;
    Double I_ESP_M7x2y_Py_a = I_ESP_N7x3y_S_a+ABY*I_ESP_M7x2y_S_a;
    Double I_ESP_M7xyz_Py_a = I_ESP_N7x2yz_S_a+ABY*I_ESP_M7xyz_S_a;
    Double I_ESP_M6x3y_Py_a = I_ESP_N6x4y_S_a+ABY*I_ESP_M6x3y_S_a;
    Double I_ESP_M6x2yz_Py_a = I_ESP_N6x3yz_S_a+ABY*I_ESP_M6x2yz_S_a;
    Double I_ESP_M6xy2z_Py_a = I_ESP_N6x2y2z_S_a+ABY*I_ESP_M6xy2z_S_a;
    Double I_ESP_M5x4y_Py_a = I_ESP_N5x5y_S_a+ABY*I_ESP_M5x4y_S_a;
    Double I_ESP_M5x3yz_Py_a = I_ESP_N5x4yz_S_a+ABY*I_ESP_M5x3yz_S_a;
    Double I_ESP_M5x2y2z_Py_a = I_ESP_N5x3y2z_S_a+ABY*I_ESP_M5x2y2z_S_a;
    Double I_ESP_M5xy3z_Py_a = I_ESP_N5x2y3z_S_a+ABY*I_ESP_M5xy3z_S_a;
    Double I_ESP_M4x5y_Py_a = I_ESP_N4x6y_S_a+ABY*I_ESP_M4x5y_S_a;
    Double I_ESP_M4x4yz_Py_a = I_ESP_N4x5yz_S_a+ABY*I_ESP_M4x4yz_S_a;
    Double I_ESP_M4x3y2z_Py_a = I_ESP_N4x4y2z_S_a+ABY*I_ESP_M4x3y2z_S_a;
    Double I_ESP_M4x2y3z_Py_a = I_ESP_N4x3y3z_S_a+ABY*I_ESP_M4x2y3z_S_a;
    Double I_ESP_M4xy4z_Py_a = I_ESP_N4x2y4z_S_a+ABY*I_ESP_M4xy4z_S_a;
    Double I_ESP_M3x6y_Py_a = I_ESP_N3x7y_S_a+ABY*I_ESP_M3x6y_S_a;
    Double I_ESP_M3x5yz_Py_a = I_ESP_N3x6yz_S_a+ABY*I_ESP_M3x5yz_S_a;
    Double I_ESP_M3x4y2z_Py_a = I_ESP_N3x5y2z_S_a+ABY*I_ESP_M3x4y2z_S_a;
    Double I_ESP_M3x3y3z_Py_a = I_ESP_N3x4y3z_S_a+ABY*I_ESP_M3x3y3z_S_a;
    Double I_ESP_M3x2y4z_Py_a = I_ESP_N3x3y4z_S_a+ABY*I_ESP_M3x2y4z_S_a;
    Double I_ESP_M3xy5z_Py_a = I_ESP_N3x2y5z_S_a+ABY*I_ESP_M3xy5z_S_a;
    Double I_ESP_M2x7y_Py_a = I_ESP_N2x8y_S_a+ABY*I_ESP_M2x7y_S_a;
    Double I_ESP_M2x6yz_Py_a = I_ESP_N2x7yz_S_a+ABY*I_ESP_M2x6yz_S_a;
    Double I_ESP_M2x5y2z_Py_a = I_ESP_N2x6y2z_S_a+ABY*I_ESP_M2x5y2z_S_a;
    Double I_ESP_M2x4y3z_Py_a = I_ESP_N2x5y3z_S_a+ABY*I_ESP_M2x4y3z_S_a;
    Double I_ESP_M2x3y4z_Py_a = I_ESP_N2x4y4z_S_a+ABY*I_ESP_M2x3y4z_S_a;
    Double I_ESP_M2x2y5z_Py_a = I_ESP_N2x3y5z_S_a+ABY*I_ESP_M2x2y5z_S_a;
    Double I_ESP_M2xy6z_Py_a = I_ESP_N2x2y6z_S_a+ABY*I_ESP_M2xy6z_S_a;
    Double I_ESP_Mx8y_Py_a = I_ESP_Nx9y_S_a+ABY*I_ESP_Mx8y_S_a;
    Double I_ESP_Mx7yz_Py_a = I_ESP_Nx8yz_S_a+ABY*I_ESP_Mx7yz_S_a;
    Double I_ESP_Mx6y2z_Py_a = I_ESP_Nx7y2z_S_a+ABY*I_ESP_Mx6y2z_S_a;
    Double I_ESP_Mx5y3z_Py_a = I_ESP_Nx6y3z_S_a+ABY*I_ESP_Mx5y3z_S_a;
    Double I_ESP_Mx4y4z_Py_a = I_ESP_Nx5y4z_S_a+ABY*I_ESP_Mx4y4z_S_a;
    Double I_ESP_Mx3y5z_Py_a = I_ESP_Nx4y5z_S_a+ABY*I_ESP_Mx3y5z_S_a;
    Double I_ESP_Mx2y6z_Py_a = I_ESP_Nx3y6z_S_a+ABY*I_ESP_Mx2y6z_S_a;
    Double I_ESP_Mxy7z_Py_a = I_ESP_Nx2y7z_S_a+ABY*I_ESP_Mxy7z_S_a;
    Double I_ESP_M9y_Py_a = I_ESP_N10y_S_a+ABY*I_ESP_M9y_S_a;
    Double I_ESP_M8yz_Py_a = I_ESP_N9yz_S_a+ABY*I_ESP_M8yz_S_a;
    Double I_ESP_M7y2z_Py_a = I_ESP_N8y2z_S_a+ABY*I_ESP_M7y2z_S_a;
    Double I_ESP_M6y3z_Py_a = I_ESP_N7y3z_S_a+ABY*I_ESP_M6y3z_S_a;
    Double I_ESP_M5y4z_Py_a = I_ESP_N6y4z_S_a+ABY*I_ESP_M5y4z_S_a;
    Double I_ESP_M4y5z_Py_a = I_ESP_N5y5z_S_a+ABY*I_ESP_M4y5z_S_a;
    Double I_ESP_M3y6z_Py_a = I_ESP_N4y6z_S_a+ABY*I_ESP_M3y6z_S_a;
    Double I_ESP_M2y7z_Py_a = I_ESP_N3y7z_S_a+ABY*I_ESP_M2y7z_S_a;
    Double I_ESP_My8z_Py_a = I_ESP_N2y8z_S_a+ABY*I_ESP_My8z_S_a;
    Double I_ESP_M7xyz_Pz_a = I_ESP_N7xy2z_S_a+ABZ*I_ESP_M7xyz_S_a;
    Double I_ESP_M7x2z_Pz_a = I_ESP_N7x3z_S_a+ABZ*I_ESP_M7x2z_S_a;
    Double I_ESP_M6x2yz_Pz_a = I_ESP_N6x2y2z_S_a+ABZ*I_ESP_M6x2yz_S_a;
    Double I_ESP_M6xy2z_Pz_a = I_ESP_N6xy3z_S_a+ABZ*I_ESP_M6xy2z_S_a;
    Double I_ESP_M6x3z_Pz_a = I_ESP_N6x4z_S_a+ABZ*I_ESP_M6x3z_S_a;
    Double I_ESP_M5x3yz_Pz_a = I_ESP_N5x3y2z_S_a+ABZ*I_ESP_M5x3yz_S_a;
    Double I_ESP_M5x2y2z_Pz_a = I_ESP_N5x2y3z_S_a+ABZ*I_ESP_M5x2y2z_S_a;
    Double I_ESP_M5xy3z_Pz_a = I_ESP_N5xy4z_S_a+ABZ*I_ESP_M5xy3z_S_a;
    Double I_ESP_M5x4z_Pz_a = I_ESP_N5x5z_S_a+ABZ*I_ESP_M5x4z_S_a;
    Double I_ESP_M4x4yz_Pz_a = I_ESP_N4x4y2z_S_a+ABZ*I_ESP_M4x4yz_S_a;
    Double I_ESP_M4x3y2z_Pz_a = I_ESP_N4x3y3z_S_a+ABZ*I_ESP_M4x3y2z_S_a;
    Double I_ESP_M4x2y3z_Pz_a = I_ESP_N4x2y4z_S_a+ABZ*I_ESP_M4x2y3z_S_a;
    Double I_ESP_M4xy4z_Pz_a = I_ESP_N4xy5z_S_a+ABZ*I_ESP_M4xy4z_S_a;
    Double I_ESP_M4x5z_Pz_a = I_ESP_N4x6z_S_a+ABZ*I_ESP_M4x5z_S_a;
    Double I_ESP_M3x5yz_Pz_a = I_ESP_N3x5y2z_S_a+ABZ*I_ESP_M3x5yz_S_a;
    Double I_ESP_M3x4y2z_Pz_a = I_ESP_N3x4y3z_S_a+ABZ*I_ESP_M3x4y2z_S_a;
    Double I_ESP_M3x3y3z_Pz_a = I_ESP_N3x3y4z_S_a+ABZ*I_ESP_M3x3y3z_S_a;
    Double I_ESP_M3x2y4z_Pz_a = I_ESP_N3x2y5z_S_a+ABZ*I_ESP_M3x2y4z_S_a;
    Double I_ESP_M3xy5z_Pz_a = I_ESP_N3xy6z_S_a+ABZ*I_ESP_M3xy5z_S_a;
    Double I_ESP_M3x6z_Pz_a = I_ESP_N3x7z_S_a+ABZ*I_ESP_M3x6z_S_a;
    Double I_ESP_M2x6yz_Pz_a = I_ESP_N2x6y2z_S_a+ABZ*I_ESP_M2x6yz_S_a;
    Double I_ESP_M2x5y2z_Pz_a = I_ESP_N2x5y3z_S_a+ABZ*I_ESP_M2x5y2z_S_a;
    Double I_ESP_M2x4y3z_Pz_a = I_ESP_N2x4y4z_S_a+ABZ*I_ESP_M2x4y3z_S_a;
    Double I_ESP_M2x3y4z_Pz_a = I_ESP_N2x3y5z_S_a+ABZ*I_ESP_M2x3y4z_S_a;
    Double I_ESP_M2x2y5z_Pz_a = I_ESP_N2x2y6z_S_a+ABZ*I_ESP_M2x2y5z_S_a;
    Double I_ESP_M2xy6z_Pz_a = I_ESP_N2xy7z_S_a+ABZ*I_ESP_M2xy6z_S_a;
    Double I_ESP_M2x7z_Pz_a = I_ESP_N2x8z_S_a+ABZ*I_ESP_M2x7z_S_a;
    Double I_ESP_Mx7yz_Pz_a = I_ESP_Nx7y2z_S_a+ABZ*I_ESP_Mx7yz_S_a;
    Double I_ESP_Mx6y2z_Pz_a = I_ESP_Nx6y3z_S_a+ABZ*I_ESP_Mx6y2z_S_a;
    Double I_ESP_Mx5y3z_Pz_a = I_ESP_Nx5y4z_S_a+ABZ*I_ESP_Mx5y3z_S_a;
    Double I_ESP_Mx4y4z_Pz_a = I_ESP_Nx4y5z_S_a+ABZ*I_ESP_Mx4y4z_S_a;
    Double I_ESP_Mx3y5z_Pz_a = I_ESP_Nx3y6z_S_a+ABZ*I_ESP_Mx3y5z_S_a;
    Double I_ESP_Mx2y6z_Pz_a = I_ESP_Nx2y7z_S_a+ABZ*I_ESP_Mx2y6z_S_a;
    Double I_ESP_Mxy7z_Pz_a = I_ESP_Nxy8z_S_a+ABZ*I_ESP_Mxy7z_S_a;
    Double I_ESP_Mx8z_Pz_a = I_ESP_Nx9z_S_a+ABZ*I_ESP_Mx8z_S_a;
    Double I_ESP_M7y2z_Pz_a = I_ESP_N7y3z_S_a+ABZ*I_ESP_M7y2z_S_a;
    Double I_ESP_M6y3z_Pz_a = I_ESP_N6y4z_S_a+ABZ*I_ESP_M6y3z_S_a;
    Double I_ESP_M5y4z_Pz_a = I_ESP_N5y5z_S_a+ABZ*I_ESP_M5y4z_S_a;
    Double I_ESP_M4y5z_Pz_a = I_ESP_N4y6z_S_a+ABZ*I_ESP_M4y5z_S_a;
    Double I_ESP_M3y6z_Pz_a = I_ESP_N3y7z_S_a+ABZ*I_ESP_M3y6z_S_a;
    Double I_ESP_M2y7z_Pz_a = I_ESP_N2y8z_S_a+ABZ*I_ESP_M2y7z_S_a;
    Double I_ESP_My8z_Pz_a = I_ESP_Ny9z_S_a+ABZ*I_ESP_My8z_S_a;
    Double I_ESP_M9z_Pz_a = I_ESP_N10z_S_a+ABZ*I_ESP_M9z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_L_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 138 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_M_P_a
     * RHS shell quartet name: SQ_ESP_L_P_a
     ************************************************************/
    Double I_ESP_L8x_D2x_a = I_ESP_M9x_Px_a+ABX*I_ESP_L8x_Px_a;
    Double I_ESP_L7xy_D2x_a = I_ESP_M8xy_Px_a+ABX*I_ESP_L7xy_Px_a;
    Double I_ESP_L7xz_D2x_a = I_ESP_M8xz_Px_a+ABX*I_ESP_L7xz_Px_a;
    Double I_ESP_L6x2y_D2x_a = I_ESP_M7x2y_Px_a+ABX*I_ESP_L6x2y_Px_a;
    Double I_ESP_L6xyz_D2x_a = I_ESP_M7xyz_Px_a+ABX*I_ESP_L6xyz_Px_a;
    Double I_ESP_L6x2z_D2x_a = I_ESP_M7x2z_Px_a+ABX*I_ESP_L6x2z_Px_a;
    Double I_ESP_L5x3y_D2x_a = I_ESP_M6x3y_Px_a+ABX*I_ESP_L5x3y_Px_a;
    Double I_ESP_L5x2yz_D2x_a = I_ESP_M6x2yz_Px_a+ABX*I_ESP_L5x2yz_Px_a;
    Double I_ESP_L5xy2z_D2x_a = I_ESP_M6xy2z_Px_a+ABX*I_ESP_L5xy2z_Px_a;
    Double I_ESP_L5x3z_D2x_a = I_ESP_M6x3z_Px_a+ABX*I_ESP_L5x3z_Px_a;
    Double I_ESP_L4x4y_D2x_a = I_ESP_M5x4y_Px_a+ABX*I_ESP_L4x4y_Px_a;
    Double I_ESP_L4x3yz_D2x_a = I_ESP_M5x3yz_Px_a+ABX*I_ESP_L4x3yz_Px_a;
    Double I_ESP_L4x2y2z_D2x_a = I_ESP_M5x2y2z_Px_a+ABX*I_ESP_L4x2y2z_Px_a;
    Double I_ESP_L4xy3z_D2x_a = I_ESP_M5xy3z_Px_a+ABX*I_ESP_L4xy3z_Px_a;
    Double I_ESP_L4x4z_D2x_a = I_ESP_M5x4z_Px_a+ABX*I_ESP_L4x4z_Px_a;
    Double I_ESP_L3x5y_D2x_a = I_ESP_M4x5y_Px_a+ABX*I_ESP_L3x5y_Px_a;
    Double I_ESP_L3x4yz_D2x_a = I_ESP_M4x4yz_Px_a+ABX*I_ESP_L3x4yz_Px_a;
    Double I_ESP_L3x3y2z_D2x_a = I_ESP_M4x3y2z_Px_a+ABX*I_ESP_L3x3y2z_Px_a;
    Double I_ESP_L3x2y3z_D2x_a = I_ESP_M4x2y3z_Px_a+ABX*I_ESP_L3x2y3z_Px_a;
    Double I_ESP_L3xy4z_D2x_a = I_ESP_M4xy4z_Px_a+ABX*I_ESP_L3xy4z_Px_a;
    Double I_ESP_L3x5z_D2x_a = I_ESP_M4x5z_Px_a+ABX*I_ESP_L3x5z_Px_a;
    Double I_ESP_L2x6y_D2x_a = I_ESP_M3x6y_Px_a+ABX*I_ESP_L2x6y_Px_a;
    Double I_ESP_L2x5yz_D2x_a = I_ESP_M3x5yz_Px_a+ABX*I_ESP_L2x5yz_Px_a;
    Double I_ESP_L2x4y2z_D2x_a = I_ESP_M3x4y2z_Px_a+ABX*I_ESP_L2x4y2z_Px_a;
    Double I_ESP_L2x3y3z_D2x_a = I_ESP_M3x3y3z_Px_a+ABX*I_ESP_L2x3y3z_Px_a;
    Double I_ESP_L2x2y4z_D2x_a = I_ESP_M3x2y4z_Px_a+ABX*I_ESP_L2x2y4z_Px_a;
    Double I_ESP_L2xy5z_D2x_a = I_ESP_M3xy5z_Px_a+ABX*I_ESP_L2xy5z_Px_a;
    Double I_ESP_L2x6z_D2x_a = I_ESP_M3x6z_Px_a+ABX*I_ESP_L2x6z_Px_a;
    Double I_ESP_Lx7y_D2x_a = I_ESP_M2x7y_Px_a+ABX*I_ESP_Lx7y_Px_a;
    Double I_ESP_Lx6yz_D2x_a = I_ESP_M2x6yz_Px_a+ABX*I_ESP_Lx6yz_Px_a;
    Double I_ESP_Lx5y2z_D2x_a = I_ESP_M2x5y2z_Px_a+ABX*I_ESP_Lx5y2z_Px_a;
    Double I_ESP_Lx4y3z_D2x_a = I_ESP_M2x4y3z_Px_a+ABX*I_ESP_Lx4y3z_Px_a;
    Double I_ESP_Lx3y4z_D2x_a = I_ESP_M2x3y4z_Px_a+ABX*I_ESP_Lx3y4z_Px_a;
    Double I_ESP_Lx2y5z_D2x_a = I_ESP_M2x2y5z_Px_a+ABX*I_ESP_Lx2y5z_Px_a;
    Double I_ESP_Lxy6z_D2x_a = I_ESP_M2xy6z_Px_a+ABX*I_ESP_Lxy6z_Px_a;
    Double I_ESP_Lx7z_D2x_a = I_ESP_M2x7z_Px_a+ABX*I_ESP_Lx7z_Px_a;
    Double I_ESP_L8y_D2x_a = I_ESP_Mx8y_Px_a+ABX*I_ESP_L8y_Px_a;
    Double I_ESP_L7yz_D2x_a = I_ESP_Mx7yz_Px_a+ABX*I_ESP_L7yz_Px_a;
    Double I_ESP_L6y2z_D2x_a = I_ESP_Mx6y2z_Px_a+ABX*I_ESP_L6y2z_Px_a;
    Double I_ESP_L5y3z_D2x_a = I_ESP_Mx5y3z_Px_a+ABX*I_ESP_L5y3z_Px_a;
    Double I_ESP_L4y4z_D2x_a = I_ESP_Mx4y4z_Px_a+ABX*I_ESP_L4y4z_Px_a;
    Double I_ESP_L3y5z_D2x_a = I_ESP_Mx3y5z_Px_a+ABX*I_ESP_L3y5z_Px_a;
    Double I_ESP_L2y6z_D2x_a = I_ESP_Mx2y6z_Px_a+ABX*I_ESP_L2y6z_Px_a;
    Double I_ESP_Ly7z_D2x_a = I_ESP_Mxy7z_Px_a+ABX*I_ESP_Ly7z_Px_a;
    Double I_ESP_L8z_D2x_a = I_ESP_Mx8z_Px_a+ABX*I_ESP_L8z_Px_a;
    Double I_ESP_L7xy_D2y_a = I_ESP_M7x2y_Py_a+ABY*I_ESP_L7xy_Py_a;
    Double I_ESP_L7xz_D2y_a = I_ESP_M7xyz_Py_a+ABY*I_ESP_L7xz_Py_a;
    Double I_ESP_L6x2y_D2y_a = I_ESP_M6x3y_Py_a+ABY*I_ESP_L6x2y_Py_a;
    Double I_ESP_L6xyz_D2y_a = I_ESP_M6x2yz_Py_a+ABY*I_ESP_L6xyz_Py_a;
    Double I_ESP_L6x2z_D2y_a = I_ESP_M6xy2z_Py_a+ABY*I_ESP_L6x2z_Py_a;
    Double I_ESP_L5x3y_D2y_a = I_ESP_M5x4y_Py_a+ABY*I_ESP_L5x3y_Py_a;
    Double I_ESP_L5x2yz_D2y_a = I_ESP_M5x3yz_Py_a+ABY*I_ESP_L5x2yz_Py_a;
    Double I_ESP_L5xy2z_D2y_a = I_ESP_M5x2y2z_Py_a+ABY*I_ESP_L5xy2z_Py_a;
    Double I_ESP_L5x3z_D2y_a = I_ESP_M5xy3z_Py_a+ABY*I_ESP_L5x3z_Py_a;
    Double I_ESP_L4x4y_D2y_a = I_ESP_M4x5y_Py_a+ABY*I_ESP_L4x4y_Py_a;
    Double I_ESP_L4x3yz_D2y_a = I_ESP_M4x4yz_Py_a+ABY*I_ESP_L4x3yz_Py_a;
    Double I_ESP_L4x2y2z_D2y_a = I_ESP_M4x3y2z_Py_a+ABY*I_ESP_L4x2y2z_Py_a;
    Double I_ESP_L4xy3z_D2y_a = I_ESP_M4x2y3z_Py_a+ABY*I_ESP_L4xy3z_Py_a;
    Double I_ESP_L4x4z_D2y_a = I_ESP_M4xy4z_Py_a+ABY*I_ESP_L4x4z_Py_a;
    Double I_ESP_L3x5y_D2y_a = I_ESP_M3x6y_Py_a+ABY*I_ESP_L3x5y_Py_a;
    Double I_ESP_L3x4yz_D2y_a = I_ESP_M3x5yz_Py_a+ABY*I_ESP_L3x4yz_Py_a;
    Double I_ESP_L3x3y2z_D2y_a = I_ESP_M3x4y2z_Py_a+ABY*I_ESP_L3x3y2z_Py_a;
    Double I_ESP_L3x2y3z_D2y_a = I_ESP_M3x3y3z_Py_a+ABY*I_ESP_L3x2y3z_Py_a;
    Double I_ESP_L3xy4z_D2y_a = I_ESP_M3x2y4z_Py_a+ABY*I_ESP_L3xy4z_Py_a;
    Double I_ESP_L3x5z_D2y_a = I_ESP_M3xy5z_Py_a+ABY*I_ESP_L3x5z_Py_a;
    Double I_ESP_L2x6y_D2y_a = I_ESP_M2x7y_Py_a+ABY*I_ESP_L2x6y_Py_a;
    Double I_ESP_L2x5yz_D2y_a = I_ESP_M2x6yz_Py_a+ABY*I_ESP_L2x5yz_Py_a;
    Double I_ESP_L2x4y2z_D2y_a = I_ESP_M2x5y2z_Py_a+ABY*I_ESP_L2x4y2z_Py_a;
    Double I_ESP_L2x3y3z_D2y_a = I_ESP_M2x4y3z_Py_a+ABY*I_ESP_L2x3y3z_Py_a;
    Double I_ESP_L2x2y4z_D2y_a = I_ESP_M2x3y4z_Py_a+ABY*I_ESP_L2x2y4z_Py_a;
    Double I_ESP_L2xy5z_D2y_a = I_ESP_M2x2y5z_Py_a+ABY*I_ESP_L2xy5z_Py_a;
    Double I_ESP_L2x6z_D2y_a = I_ESP_M2xy6z_Py_a+ABY*I_ESP_L2x6z_Py_a;
    Double I_ESP_Lx7y_D2y_a = I_ESP_Mx8y_Py_a+ABY*I_ESP_Lx7y_Py_a;
    Double I_ESP_Lx6yz_D2y_a = I_ESP_Mx7yz_Py_a+ABY*I_ESP_Lx6yz_Py_a;
    Double I_ESP_Lx5y2z_D2y_a = I_ESP_Mx6y2z_Py_a+ABY*I_ESP_Lx5y2z_Py_a;
    Double I_ESP_Lx4y3z_D2y_a = I_ESP_Mx5y3z_Py_a+ABY*I_ESP_Lx4y3z_Py_a;
    Double I_ESP_Lx3y4z_D2y_a = I_ESP_Mx4y4z_Py_a+ABY*I_ESP_Lx3y4z_Py_a;
    Double I_ESP_Lx2y5z_D2y_a = I_ESP_Mx3y5z_Py_a+ABY*I_ESP_Lx2y5z_Py_a;
    Double I_ESP_Lxy6z_D2y_a = I_ESP_Mx2y6z_Py_a+ABY*I_ESP_Lxy6z_Py_a;
    Double I_ESP_Lx7z_D2y_a = I_ESP_Mxy7z_Py_a+ABY*I_ESP_Lx7z_Py_a;
    Double I_ESP_L8y_D2y_a = I_ESP_M9y_Py_a+ABY*I_ESP_L8y_Py_a;
    Double I_ESP_L7yz_D2y_a = I_ESP_M8yz_Py_a+ABY*I_ESP_L7yz_Py_a;
    Double I_ESP_L6y2z_D2y_a = I_ESP_M7y2z_Py_a+ABY*I_ESP_L6y2z_Py_a;
    Double I_ESP_L5y3z_D2y_a = I_ESP_M6y3z_Py_a+ABY*I_ESP_L5y3z_Py_a;
    Double I_ESP_L4y4z_D2y_a = I_ESP_M5y4z_Py_a+ABY*I_ESP_L4y4z_Py_a;
    Double I_ESP_L3y5z_D2y_a = I_ESP_M4y5z_Py_a+ABY*I_ESP_L3y5z_Py_a;
    Double I_ESP_L2y6z_D2y_a = I_ESP_M3y6z_Py_a+ABY*I_ESP_L2y6z_Py_a;
    Double I_ESP_Ly7z_D2y_a = I_ESP_M2y7z_Py_a+ABY*I_ESP_Ly7z_Py_a;
    Double I_ESP_L8z_D2y_a = I_ESP_My8z_Py_a+ABY*I_ESP_L8z_Py_a;
    Double I_ESP_L7xy_D2z_a = I_ESP_M7xyz_Pz_a+ABZ*I_ESP_L7xy_Pz_a;
    Double I_ESP_L7xz_D2z_a = I_ESP_M7x2z_Pz_a+ABZ*I_ESP_L7xz_Pz_a;
    Double I_ESP_L6x2y_D2z_a = I_ESP_M6x2yz_Pz_a+ABZ*I_ESP_L6x2y_Pz_a;
    Double I_ESP_L6xyz_D2z_a = I_ESP_M6xy2z_Pz_a+ABZ*I_ESP_L6xyz_Pz_a;
    Double I_ESP_L6x2z_D2z_a = I_ESP_M6x3z_Pz_a+ABZ*I_ESP_L6x2z_Pz_a;
    Double I_ESP_L5x3y_D2z_a = I_ESP_M5x3yz_Pz_a+ABZ*I_ESP_L5x3y_Pz_a;
    Double I_ESP_L5x2yz_D2z_a = I_ESP_M5x2y2z_Pz_a+ABZ*I_ESP_L5x2yz_Pz_a;
    Double I_ESP_L5xy2z_D2z_a = I_ESP_M5xy3z_Pz_a+ABZ*I_ESP_L5xy2z_Pz_a;
    Double I_ESP_L5x3z_D2z_a = I_ESP_M5x4z_Pz_a+ABZ*I_ESP_L5x3z_Pz_a;
    Double I_ESP_L4x4y_D2z_a = I_ESP_M4x4yz_Pz_a+ABZ*I_ESP_L4x4y_Pz_a;
    Double I_ESP_L4x3yz_D2z_a = I_ESP_M4x3y2z_Pz_a+ABZ*I_ESP_L4x3yz_Pz_a;
    Double I_ESP_L4x2y2z_D2z_a = I_ESP_M4x2y3z_Pz_a+ABZ*I_ESP_L4x2y2z_Pz_a;
    Double I_ESP_L4xy3z_D2z_a = I_ESP_M4xy4z_Pz_a+ABZ*I_ESP_L4xy3z_Pz_a;
    Double I_ESP_L4x4z_D2z_a = I_ESP_M4x5z_Pz_a+ABZ*I_ESP_L4x4z_Pz_a;
    Double I_ESP_L3x5y_D2z_a = I_ESP_M3x5yz_Pz_a+ABZ*I_ESP_L3x5y_Pz_a;
    Double I_ESP_L3x4yz_D2z_a = I_ESP_M3x4y2z_Pz_a+ABZ*I_ESP_L3x4yz_Pz_a;
    Double I_ESP_L3x3y2z_D2z_a = I_ESP_M3x3y3z_Pz_a+ABZ*I_ESP_L3x3y2z_Pz_a;
    Double I_ESP_L3x2y3z_D2z_a = I_ESP_M3x2y4z_Pz_a+ABZ*I_ESP_L3x2y3z_Pz_a;
    Double I_ESP_L3xy4z_D2z_a = I_ESP_M3xy5z_Pz_a+ABZ*I_ESP_L3xy4z_Pz_a;
    Double I_ESP_L3x5z_D2z_a = I_ESP_M3x6z_Pz_a+ABZ*I_ESP_L3x5z_Pz_a;
    Double I_ESP_L2x6y_D2z_a = I_ESP_M2x6yz_Pz_a+ABZ*I_ESP_L2x6y_Pz_a;
    Double I_ESP_L2x5yz_D2z_a = I_ESP_M2x5y2z_Pz_a+ABZ*I_ESP_L2x5yz_Pz_a;
    Double I_ESP_L2x4y2z_D2z_a = I_ESP_M2x4y3z_Pz_a+ABZ*I_ESP_L2x4y2z_Pz_a;
    Double I_ESP_L2x3y3z_D2z_a = I_ESP_M2x3y4z_Pz_a+ABZ*I_ESP_L2x3y3z_Pz_a;
    Double I_ESP_L2x2y4z_D2z_a = I_ESP_M2x2y5z_Pz_a+ABZ*I_ESP_L2x2y4z_Pz_a;
    Double I_ESP_L2xy5z_D2z_a = I_ESP_M2xy6z_Pz_a+ABZ*I_ESP_L2xy5z_Pz_a;
    Double I_ESP_L2x6z_D2z_a = I_ESP_M2x7z_Pz_a+ABZ*I_ESP_L2x6z_Pz_a;
    Double I_ESP_Lx7y_D2z_a = I_ESP_Mx7yz_Pz_a+ABZ*I_ESP_Lx7y_Pz_a;
    Double I_ESP_Lx6yz_D2z_a = I_ESP_Mx6y2z_Pz_a+ABZ*I_ESP_Lx6yz_Pz_a;
    Double I_ESP_Lx5y2z_D2z_a = I_ESP_Mx5y3z_Pz_a+ABZ*I_ESP_Lx5y2z_Pz_a;
    Double I_ESP_Lx4y3z_D2z_a = I_ESP_Mx4y4z_Pz_a+ABZ*I_ESP_Lx4y3z_Pz_a;
    Double I_ESP_Lx3y4z_D2z_a = I_ESP_Mx3y5z_Pz_a+ABZ*I_ESP_Lx3y4z_Pz_a;
    Double I_ESP_Lx2y5z_D2z_a = I_ESP_Mx2y6z_Pz_a+ABZ*I_ESP_Lx2y5z_Pz_a;
    Double I_ESP_Lxy6z_D2z_a = I_ESP_Mxy7z_Pz_a+ABZ*I_ESP_Lxy6z_Pz_a;
    Double I_ESP_Lx7z_D2z_a = I_ESP_Mx8z_Pz_a+ABZ*I_ESP_Lx7z_Pz_a;
    Double I_ESP_L7yz_D2z_a = I_ESP_M7y2z_Pz_a+ABZ*I_ESP_L7yz_Pz_a;
    Double I_ESP_L6y2z_D2z_a = I_ESP_M6y3z_Pz_a+ABZ*I_ESP_L6y2z_Pz_a;
    Double I_ESP_L5y3z_D2z_a = I_ESP_M5y4z_Pz_a+ABZ*I_ESP_L5y3z_Pz_a;
    Double I_ESP_L4y4z_D2z_a = I_ESP_M4y5z_Pz_a+ABZ*I_ESP_L4y4z_Pz_a;
    Double I_ESP_L3y5z_D2z_a = I_ESP_M3y6z_Pz_a+ABZ*I_ESP_L3y5z_Pz_a;
    Double I_ESP_L2y6z_D2z_a = I_ESP_M2y7z_Pz_a+ABZ*I_ESP_L2y6z_Pz_a;
    Double I_ESP_Ly7z_D2z_a = I_ESP_My8z_Pz_a+ABZ*I_ESP_Ly7z_Pz_a;
    Double I_ESP_L8z_D2z_a = I_ESP_M9z_Pz_a+ABZ*I_ESP_L8z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_K_F_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 105 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_D_a
     * RHS shell quartet name: SQ_ESP_K_D_a
     ************************************************************/
    Double I_ESP_K7x_F3x_a = I_ESP_L8x_D2x_a+ABX*I_ESP_K7x_D2x_a;
    Double I_ESP_K6xy_F3x_a = I_ESP_L7xy_D2x_a+ABX*I_ESP_K6xy_D2x_a;
    Double I_ESP_K6xz_F3x_a = I_ESP_L7xz_D2x_a+ABX*I_ESP_K6xz_D2x_a;
    Double I_ESP_K5x2y_F3x_a = I_ESP_L6x2y_D2x_a+ABX*I_ESP_K5x2y_D2x_a;
    Double I_ESP_K5xyz_F3x_a = I_ESP_L6xyz_D2x_a+ABX*I_ESP_K5xyz_D2x_a;
    Double I_ESP_K5x2z_F3x_a = I_ESP_L6x2z_D2x_a+ABX*I_ESP_K5x2z_D2x_a;
    Double I_ESP_K4x3y_F3x_a = I_ESP_L5x3y_D2x_a+ABX*I_ESP_K4x3y_D2x_a;
    Double I_ESP_K4x2yz_F3x_a = I_ESP_L5x2yz_D2x_a+ABX*I_ESP_K4x2yz_D2x_a;
    Double I_ESP_K4xy2z_F3x_a = I_ESP_L5xy2z_D2x_a+ABX*I_ESP_K4xy2z_D2x_a;
    Double I_ESP_K4x3z_F3x_a = I_ESP_L5x3z_D2x_a+ABX*I_ESP_K4x3z_D2x_a;
    Double I_ESP_K3x4y_F3x_a = I_ESP_L4x4y_D2x_a+ABX*I_ESP_K3x4y_D2x_a;
    Double I_ESP_K3x3yz_F3x_a = I_ESP_L4x3yz_D2x_a+ABX*I_ESP_K3x3yz_D2x_a;
    Double I_ESP_K3x2y2z_F3x_a = I_ESP_L4x2y2z_D2x_a+ABX*I_ESP_K3x2y2z_D2x_a;
    Double I_ESP_K3xy3z_F3x_a = I_ESP_L4xy3z_D2x_a+ABX*I_ESP_K3xy3z_D2x_a;
    Double I_ESP_K3x4z_F3x_a = I_ESP_L4x4z_D2x_a+ABX*I_ESP_K3x4z_D2x_a;
    Double I_ESP_K2x5y_F3x_a = I_ESP_L3x5y_D2x_a+ABX*I_ESP_K2x5y_D2x_a;
    Double I_ESP_K2x4yz_F3x_a = I_ESP_L3x4yz_D2x_a+ABX*I_ESP_K2x4yz_D2x_a;
    Double I_ESP_K2x3y2z_F3x_a = I_ESP_L3x3y2z_D2x_a+ABX*I_ESP_K2x3y2z_D2x_a;
    Double I_ESP_K2x2y3z_F3x_a = I_ESP_L3x2y3z_D2x_a+ABX*I_ESP_K2x2y3z_D2x_a;
    Double I_ESP_K2xy4z_F3x_a = I_ESP_L3xy4z_D2x_a+ABX*I_ESP_K2xy4z_D2x_a;
    Double I_ESP_K2x5z_F3x_a = I_ESP_L3x5z_D2x_a+ABX*I_ESP_K2x5z_D2x_a;
    Double I_ESP_Kx6y_F3x_a = I_ESP_L2x6y_D2x_a+ABX*I_ESP_Kx6y_D2x_a;
    Double I_ESP_Kx5yz_F3x_a = I_ESP_L2x5yz_D2x_a+ABX*I_ESP_Kx5yz_D2x_a;
    Double I_ESP_Kx4y2z_F3x_a = I_ESP_L2x4y2z_D2x_a+ABX*I_ESP_Kx4y2z_D2x_a;
    Double I_ESP_Kx3y3z_F3x_a = I_ESP_L2x3y3z_D2x_a+ABX*I_ESP_Kx3y3z_D2x_a;
    Double I_ESP_Kx2y4z_F3x_a = I_ESP_L2x2y4z_D2x_a+ABX*I_ESP_Kx2y4z_D2x_a;
    Double I_ESP_Kxy5z_F3x_a = I_ESP_L2xy5z_D2x_a+ABX*I_ESP_Kxy5z_D2x_a;
    Double I_ESP_Kx6z_F3x_a = I_ESP_L2x6z_D2x_a+ABX*I_ESP_Kx6z_D2x_a;
    Double I_ESP_K7y_F3x_a = I_ESP_Lx7y_D2x_a+ABX*I_ESP_K7y_D2x_a;
    Double I_ESP_K6yz_F3x_a = I_ESP_Lx6yz_D2x_a+ABX*I_ESP_K6yz_D2x_a;
    Double I_ESP_K5y2z_F3x_a = I_ESP_Lx5y2z_D2x_a+ABX*I_ESP_K5y2z_D2x_a;
    Double I_ESP_K4y3z_F3x_a = I_ESP_Lx4y3z_D2x_a+ABX*I_ESP_K4y3z_D2x_a;
    Double I_ESP_K3y4z_F3x_a = I_ESP_Lx3y4z_D2x_a+ABX*I_ESP_K3y4z_D2x_a;
    Double I_ESP_K2y5z_F3x_a = I_ESP_Lx2y5z_D2x_a+ABX*I_ESP_K2y5z_D2x_a;
    Double I_ESP_Ky6z_F3x_a = I_ESP_Lxy6z_D2x_a+ABX*I_ESP_Ky6z_D2x_a;
    Double I_ESP_K7z_F3x_a = I_ESP_Lx7z_D2x_a+ABX*I_ESP_K7z_D2x_a;
    Double I_ESP_K6xy_F2xy_a = I_ESP_L6x2y_D2x_a+ABY*I_ESP_K6xy_D2x_a;
    Double I_ESP_K6xz_F2xy_a = I_ESP_L6xyz_D2x_a+ABY*I_ESP_K6xz_D2x_a;
    Double I_ESP_K5x2y_F2xy_a = I_ESP_L5x3y_D2x_a+ABY*I_ESP_K5x2y_D2x_a;
    Double I_ESP_K5xyz_F2xy_a = I_ESP_L5x2yz_D2x_a+ABY*I_ESP_K5xyz_D2x_a;
    Double I_ESP_K5x2z_F2xy_a = I_ESP_L5xy2z_D2x_a+ABY*I_ESP_K5x2z_D2x_a;
    Double I_ESP_K4x3y_F2xy_a = I_ESP_L4x4y_D2x_a+ABY*I_ESP_K4x3y_D2x_a;
    Double I_ESP_K4x2yz_F2xy_a = I_ESP_L4x3yz_D2x_a+ABY*I_ESP_K4x2yz_D2x_a;
    Double I_ESP_K4xy2z_F2xy_a = I_ESP_L4x2y2z_D2x_a+ABY*I_ESP_K4xy2z_D2x_a;
    Double I_ESP_K4x3z_F2xy_a = I_ESP_L4xy3z_D2x_a+ABY*I_ESP_K4x3z_D2x_a;
    Double I_ESP_K3x4y_F2xy_a = I_ESP_L3x5y_D2x_a+ABY*I_ESP_K3x4y_D2x_a;
    Double I_ESP_K3x3yz_F2xy_a = I_ESP_L3x4yz_D2x_a+ABY*I_ESP_K3x3yz_D2x_a;
    Double I_ESP_K3x2y2z_F2xy_a = I_ESP_L3x3y2z_D2x_a+ABY*I_ESP_K3x2y2z_D2x_a;
    Double I_ESP_K3xy3z_F2xy_a = I_ESP_L3x2y3z_D2x_a+ABY*I_ESP_K3xy3z_D2x_a;
    Double I_ESP_K3x4z_F2xy_a = I_ESP_L3xy4z_D2x_a+ABY*I_ESP_K3x4z_D2x_a;
    Double I_ESP_K2x5y_F2xy_a = I_ESP_L2x6y_D2x_a+ABY*I_ESP_K2x5y_D2x_a;
    Double I_ESP_K2x4yz_F2xy_a = I_ESP_L2x5yz_D2x_a+ABY*I_ESP_K2x4yz_D2x_a;
    Double I_ESP_K2x3y2z_F2xy_a = I_ESP_L2x4y2z_D2x_a+ABY*I_ESP_K2x3y2z_D2x_a;
    Double I_ESP_K2x2y3z_F2xy_a = I_ESP_L2x3y3z_D2x_a+ABY*I_ESP_K2x2y3z_D2x_a;
    Double I_ESP_K2xy4z_F2xy_a = I_ESP_L2x2y4z_D2x_a+ABY*I_ESP_K2xy4z_D2x_a;
    Double I_ESP_K2x5z_F2xy_a = I_ESP_L2xy5z_D2x_a+ABY*I_ESP_K2x5z_D2x_a;
    Double I_ESP_Kx6y_F2xy_a = I_ESP_Lx7y_D2x_a+ABY*I_ESP_Kx6y_D2x_a;
    Double I_ESP_Kx5yz_F2xy_a = I_ESP_Lx6yz_D2x_a+ABY*I_ESP_Kx5yz_D2x_a;
    Double I_ESP_Kx4y2z_F2xy_a = I_ESP_Lx5y2z_D2x_a+ABY*I_ESP_Kx4y2z_D2x_a;
    Double I_ESP_Kx3y3z_F2xy_a = I_ESP_Lx4y3z_D2x_a+ABY*I_ESP_Kx3y3z_D2x_a;
    Double I_ESP_Kx2y4z_F2xy_a = I_ESP_Lx3y4z_D2x_a+ABY*I_ESP_Kx2y4z_D2x_a;
    Double I_ESP_Kxy5z_F2xy_a = I_ESP_Lx2y5z_D2x_a+ABY*I_ESP_Kxy5z_D2x_a;
    Double I_ESP_Kx6z_F2xy_a = I_ESP_Lxy6z_D2x_a+ABY*I_ESP_Kx6z_D2x_a;
    Double I_ESP_K7y_F2xy_a = I_ESP_L8y_D2x_a+ABY*I_ESP_K7y_D2x_a;
    Double I_ESP_K6yz_F2xy_a = I_ESP_L7yz_D2x_a+ABY*I_ESP_K6yz_D2x_a;
    Double I_ESP_K5y2z_F2xy_a = I_ESP_L6y2z_D2x_a+ABY*I_ESP_K5y2z_D2x_a;
    Double I_ESP_K4y3z_F2xy_a = I_ESP_L5y3z_D2x_a+ABY*I_ESP_K4y3z_D2x_a;
    Double I_ESP_K3y4z_F2xy_a = I_ESP_L4y4z_D2x_a+ABY*I_ESP_K3y4z_D2x_a;
    Double I_ESP_K2y5z_F2xy_a = I_ESP_L3y5z_D2x_a+ABY*I_ESP_K2y5z_D2x_a;
    Double I_ESP_Ky6z_F2xy_a = I_ESP_L2y6z_D2x_a+ABY*I_ESP_Ky6z_D2x_a;
    Double I_ESP_K7z_F2xy_a = I_ESP_Ly7z_D2x_a+ABY*I_ESP_K7z_D2x_a;
    Double I_ESP_K6xz_F2xz_a = I_ESP_L6x2z_D2x_a+ABZ*I_ESP_K6xz_D2x_a;
    Double I_ESP_K5xyz_F2xz_a = I_ESP_L5xy2z_D2x_a+ABZ*I_ESP_K5xyz_D2x_a;
    Double I_ESP_K5x2z_F2xz_a = I_ESP_L5x3z_D2x_a+ABZ*I_ESP_K5x2z_D2x_a;
    Double I_ESP_K4x2yz_F2xz_a = I_ESP_L4x2y2z_D2x_a+ABZ*I_ESP_K4x2yz_D2x_a;
    Double I_ESP_K4xy2z_F2xz_a = I_ESP_L4xy3z_D2x_a+ABZ*I_ESP_K4xy2z_D2x_a;
    Double I_ESP_K4x3z_F2xz_a = I_ESP_L4x4z_D2x_a+ABZ*I_ESP_K4x3z_D2x_a;
    Double I_ESP_K3x3yz_F2xz_a = I_ESP_L3x3y2z_D2x_a+ABZ*I_ESP_K3x3yz_D2x_a;
    Double I_ESP_K3x2y2z_F2xz_a = I_ESP_L3x2y3z_D2x_a+ABZ*I_ESP_K3x2y2z_D2x_a;
    Double I_ESP_K3xy3z_F2xz_a = I_ESP_L3xy4z_D2x_a+ABZ*I_ESP_K3xy3z_D2x_a;
    Double I_ESP_K3x4z_F2xz_a = I_ESP_L3x5z_D2x_a+ABZ*I_ESP_K3x4z_D2x_a;
    Double I_ESP_K2x4yz_F2xz_a = I_ESP_L2x4y2z_D2x_a+ABZ*I_ESP_K2x4yz_D2x_a;
    Double I_ESP_K2x3y2z_F2xz_a = I_ESP_L2x3y3z_D2x_a+ABZ*I_ESP_K2x3y2z_D2x_a;
    Double I_ESP_K2x2y3z_F2xz_a = I_ESP_L2x2y4z_D2x_a+ABZ*I_ESP_K2x2y3z_D2x_a;
    Double I_ESP_K2xy4z_F2xz_a = I_ESP_L2xy5z_D2x_a+ABZ*I_ESP_K2xy4z_D2x_a;
    Double I_ESP_K2x5z_F2xz_a = I_ESP_L2x6z_D2x_a+ABZ*I_ESP_K2x5z_D2x_a;
    Double I_ESP_Kx5yz_F2xz_a = I_ESP_Lx5y2z_D2x_a+ABZ*I_ESP_Kx5yz_D2x_a;
    Double I_ESP_Kx4y2z_F2xz_a = I_ESP_Lx4y3z_D2x_a+ABZ*I_ESP_Kx4y2z_D2x_a;
    Double I_ESP_Kx3y3z_F2xz_a = I_ESP_Lx3y4z_D2x_a+ABZ*I_ESP_Kx3y3z_D2x_a;
    Double I_ESP_Kx2y4z_F2xz_a = I_ESP_Lx2y5z_D2x_a+ABZ*I_ESP_Kx2y4z_D2x_a;
    Double I_ESP_Kxy5z_F2xz_a = I_ESP_Lxy6z_D2x_a+ABZ*I_ESP_Kxy5z_D2x_a;
    Double I_ESP_Kx6z_F2xz_a = I_ESP_Lx7z_D2x_a+ABZ*I_ESP_Kx6z_D2x_a;
    Double I_ESP_K6yz_F2xz_a = I_ESP_L6y2z_D2x_a+ABZ*I_ESP_K6yz_D2x_a;
    Double I_ESP_K5y2z_F2xz_a = I_ESP_L5y3z_D2x_a+ABZ*I_ESP_K5y2z_D2x_a;
    Double I_ESP_K4y3z_F2xz_a = I_ESP_L4y4z_D2x_a+ABZ*I_ESP_K4y3z_D2x_a;
    Double I_ESP_K3y4z_F2xz_a = I_ESP_L3y5z_D2x_a+ABZ*I_ESP_K3y4z_D2x_a;
    Double I_ESP_K2y5z_F2xz_a = I_ESP_L2y6z_D2x_a+ABZ*I_ESP_K2y5z_D2x_a;
    Double I_ESP_Ky6z_F2xz_a = I_ESP_Ly7z_D2x_a+ABZ*I_ESP_Ky6z_D2x_a;
    Double I_ESP_K7z_F2xz_a = I_ESP_L8z_D2x_a+ABZ*I_ESP_K7z_D2x_a;
    Double I_ESP_K6xz_Fx2y_a = I_ESP_L7xz_D2y_a+ABX*I_ESP_K6xz_D2y_a;
    Double I_ESP_K5xyz_Fx2y_a = I_ESP_L6xyz_D2y_a+ABX*I_ESP_K5xyz_D2y_a;
    Double I_ESP_K5x2z_Fx2y_a = I_ESP_L6x2z_D2y_a+ABX*I_ESP_K5x2z_D2y_a;
    Double I_ESP_K4x2yz_Fx2y_a = I_ESP_L5x2yz_D2y_a+ABX*I_ESP_K4x2yz_D2y_a;
    Double I_ESP_K4xy2z_Fx2y_a = I_ESP_L5xy2z_D2y_a+ABX*I_ESP_K4xy2z_D2y_a;
    Double I_ESP_K4x3z_Fx2y_a = I_ESP_L5x3z_D2y_a+ABX*I_ESP_K4x3z_D2y_a;
    Double I_ESP_K3x3yz_Fx2y_a = I_ESP_L4x3yz_D2y_a+ABX*I_ESP_K3x3yz_D2y_a;
    Double I_ESP_K3x2y2z_Fx2y_a = I_ESP_L4x2y2z_D2y_a+ABX*I_ESP_K3x2y2z_D2y_a;
    Double I_ESP_K3xy3z_Fx2y_a = I_ESP_L4xy3z_D2y_a+ABX*I_ESP_K3xy3z_D2y_a;
    Double I_ESP_K3x4z_Fx2y_a = I_ESP_L4x4z_D2y_a+ABX*I_ESP_K3x4z_D2y_a;
    Double I_ESP_K2x4yz_Fx2y_a = I_ESP_L3x4yz_D2y_a+ABX*I_ESP_K2x4yz_D2y_a;
    Double I_ESP_K2x3y2z_Fx2y_a = I_ESP_L3x3y2z_D2y_a+ABX*I_ESP_K2x3y2z_D2y_a;
    Double I_ESP_K2x2y3z_Fx2y_a = I_ESP_L3x2y3z_D2y_a+ABX*I_ESP_K2x2y3z_D2y_a;
    Double I_ESP_K2xy4z_Fx2y_a = I_ESP_L3xy4z_D2y_a+ABX*I_ESP_K2xy4z_D2y_a;
    Double I_ESP_K2x5z_Fx2y_a = I_ESP_L3x5z_D2y_a+ABX*I_ESP_K2x5z_D2y_a;
    Double I_ESP_Kx5yz_Fx2y_a = I_ESP_L2x5yz_D2y_a+ABX*I_ESP_Kx5yz_D2y_a;
    Double I_ESP_Kx4y2z_Fx2y_a = I_ESP_L2x4y2z_D2y_a+ABX*I_ESP_Kx4y2z_D2y_a;
    Double I_ESP_Kx3y3z_Fx2y_a = I_ESP_L2x3y3z_D2y_a+ABX*I_ESP_Kx3y3z_D2y_a;
    Double I_ESP_Kx2y4z_Fx2y_a = I_ESP_L2x2y4z_D2y_a+ABX*I_ESP_Kx2y4z_D2y_a;
    Double I_ESP_Kxy5z_Fx2y_a = I_ESP_L2xy5z_D2y_a+ABX*I_ESP_Kxy5z_D2y_a;
    Double I_ESP_Kx6z_Fx2y_a = I_ESP_L2x6z_D2y_a+ABX*I_ESP_Kx6z_D2y_a;
    Double I_ESP_K6yz_Fx2y_a = I_ESP_Lx6yz_D2y_a+ABX*I_ESP_K6yz_D2y_a;
    Double I_ESP_K5y2z_Fx2y_a = I_ESP_Lx5y2z_D2y_a+ABX*I_ESP_K5y2z_D2y_a;
    Double I_ESP_K4y3z_Fx2y_a = I_ESP_Lx4y3z_D2y_a+ABX*I_ESP_K4y3z_D2y_a;
    Double I_ESP_K3y4z_Fx2y_a = I_ESP_Lx3y4z_D2y_a+ABX*I_ESP_K3y4z_D2y_a;
    Double I_ESP_K2y5z_Fx2y_a = I_ESP_Lx2y5z_D2y_a+ABX*I_ESP_K2y5z_D2y_a;
    Double I_ESP_Ky6z_Fx2y_a = I_ESP_Lxy6z_D2y_a+ABX*I_ESP_Ky6z_D2y_a;
    Double I_ESP_K7z_Fx2y_a = I_ESP_Lx7z_D2y_a+ABX*I_ESP_K7z_D2y_a;
    Double I_ESP_K6xy_Fx2z_a = I_ESP_L7xy_D2z_a+ABX*I_ESP_K6xy_D2z_a;
    Double I_ESP_K5x2y_Fx2z_a = I_ESP_L6x2y_D2z_a+ABX*I_ESP_K5x2y_D2z_a;
    Double I_ESP_K5xyz_Fx2z_a = I_ESP_L6xyz_D2z_a+ABX*I_ESP_K5xyz_D2z_a;
    Double I_ESP_K4x3y_Fx2z_a = I_ESP_L5x3y_D2z_a+ABX*I_ESP_K4x3y_D2z_a;
    Double I_ESP_K4x2yz_Fx2z_a = I_ESP_L5x2yz_D2z_a+ABX*I_ESP_K4x2yz_D2z_a;
    Double I_ESP_K4xy2z_Fx2z_a = I_ESP_L5xy2z_D2z_a+ABX*I_ESP_K4xy2z_D2z_a;
    Double I_ESP_K3x4y_Fx2z_a = I_ESP_L4x4y_D2z_a+ABX*I_ESP_K3x4y_D2z_a;
    Double I_ESP_K3x3yz_Fx2z_a = I_ESP_L4x3yz_D2z_a+ABX*I_ESP_K3x3yz_D2z_a;
    Double I_ESP_K3x2y2z_Fx2z_a = I_ESP_L4x2y2z_D2z_a+ABX*I_ESP_K3x2y2z_D2z_a;
    Double I_ESP_K3xy3z_Fx2z_a = I_ESP_L4xy3z_D2z_a+ABX*I_ESP_K3xy3z_D2z_a;
    Double I_ESP_K2x5y_Fx2z_a = I_ESP_L3x5y_D2z_a+ABX*I_ESP_K2x5y_D2z_a;
    Double I_ESP_K2x4yz_Fx2z_a = I_ESP_L3x4yz_D2z_a+ABX*I_ESP_K2x4yz_D2z_a;
    Double I_ESP_K2x3y2z_Fx2z_a = I_ESP_L3x3y2z_D2z_a+ABX*I_ESP_K2x3y2z_D2z_a;
    Double I_ESP_K2x2y3z_Fx2z_a = I_ESP_L3x2y3z_D2z_a+ABX*I_ESP_K2x2y3z_D2z_a;
    Double I_ESP_K2xy4z_Fx2z_a = I_ESP_L3xy4z_D2z_a+ABX*I_ESP_K2xy4z_D2z_a;
    Double I_ESP_Kx6y_Fx2z_a = I_ESP_L2x6y_D2z_a+ABX*I_ESP_Kx6y_D2z_a;
    Double I_ESP_Kx5yz_Fx2z_a = I_ESP_L2x5yz_D2z_a+ABX*I_ESP_Kx5yz_D2z_a;
    Double I_ESP_Kx4y2z_Fx2z_a = I_ESP_L2x4y2z_D2z_a+ABX*I_ESP_Kx4y2z_D2z_a;
    Double I_ESP_Kx3y3z_Fx2z_a = I_ESP_L2x3y3z_D2z_a+ABX*I_ESP_Kx3y3z_D2z_a;
    Double I_ESP_Kx2y4z_Fx2z_a = I_ESP_L2x2y4z_D2z_a+ABX*I_ESP_Kx2y4z_D2z_a;
    Double I_ESP_Kxy5z_Fx2z_a = I_ESP_L2xy5z_D2z_a+ABX*I_ESP_Kxy5z_D2z_a;
    Double I_ESP_K7y_Fx2z_a = I_ESP_Lx7y_D2z_a+ABX*I_ESP_K7y_D2z_a;
    Double I_ESP_K6yz_Fx2z_a = I_ESP_Lx6yz_D2z_a+ABX*I_ESP_K6yz_D2z_a;
    Double I_ESP_K5y2z_Fx2z_a = I_ESP_Lx5y2z_D2z_a+ABX*I_ESP_K5y2z_D2z_a;
    Double I_ESP_K4y3z_Fx2z_a = I_ESP_Lx4y3z_D2z_a+ABX*I_ESP_K4y3z_D2z_a;
    Double I_ESP_K3y4z_Fx2z_a = I_ESP_Lx3y4z_D2z_a+ABX*I_ESP_K3y4z_D2z_a;
    Double I_ESP_K2y5z_Fx2z_a = I_ESP_Lx2y5z_D2z_a+ABX*I_ESP_K2y5z_D2z_a;
    Double I_ESP_Ky6z_Fx2z_a = I_ESP_Lxy6z_D2z_a+ABX*I_ESP_Ky6z_D2z_a;
    Double I_ESP_K7x_F3y_a = I_ESP_L7xy_D2y_a+ABY*I_ESP_K7x_D2y_a;
    Double I_ESP_K6xy_F3y_a = I_ESP_L6x2y_D2y_a+ABY*I_ESP_K6xy_D2y_a;
    Double I_ESP_K6xz_F3y_a = I_ESP_L6xyz_D2y_a+ABY*I_ESP_K6xz_D2y_a;
    Double I_ESP_K5x2y_F3y_a = I_ESP_L5x3y_D2y_a+ABY*I_ESP_K5x2y_D2y_a;
    Double I_ESP_K5xyz_F3y_a = I_ESP_L5x2yz_D2y_a+ABY*I_ESP_K5xyz_D2y_a;
    Double I_ESP_K5x2z_F3y_a = I_ESP_L5xy2z_D2y_a+ABY*I_ESP_K5x2z_D2y_a;
    Double I_ESP_K4x3y_F3y_a = I_ESP_L4x4y_D2y_a+ABY*I_ESP_K4x3y_D2y_a;
    Double I_ESP_K4x2yz_F3y_a = I_ESP_L4x3yz_D2y_a+ABY*I_ESP_K4x2yz_D2y_a;
    Double I_ESP_K4xy2z_F3y_a = I_ESP_L4x2y2z_D2y_a+ABY*I_ESP_K4xy2z_D2y_a;
    Double I_ESP_K4x3z_F3y_a = I_ESP_L4xy3z_D2y_a+ABY*I_ESP_K4x3z_D2y_a;
    Double I_ESP_K3x4y_F3y_a = I_ESP_L3x5y_D2y_a+ABY*I_ESP_K3x4y_D2y_a;
    Double I_ESP_K3x3yz_F3y_a = I_ESP_L3x4yz_D2y_a+ABY*I_ESP_K3x3yz_D2y_a;
    Double I_ESP_K3x2y2z_F3y_a = I_ESP_L3x3y2z_D2y_a+ABY*I_ESP_K3x2y2z_D2y_a;
    Double I_ESP_K3xy3z_F3y_a = I_ESP_L3x2y3z_D2y_a+ABY*I_ESP_K3xy3z_D2y_a;
    Double I_ESP_K3x4z_F3y_a = I_ESP_L3xy4z_D2y_a+ABY*I_ESP_K3x4z_D2y_a;
    Double I_ESP_K2x5y_F3y_a = I_ESP_L2x6y_D2y_a+ABY*I_ESP_K2x5y_D2y_a;
    Double I_ESP_K2x4yz_F3y_a = I_ESP_L2x5yz_D2y_a+ABY*I_ESP_K2x4yz_D2y_a;
    Double I_ESP_K2x3y2z_F3y_a = I_ESP_L2x4y2z_D2y_a+ABY*I_ESP_K2x3y2z_D2y_a;
    Double I_ESP_K2x2y3z_F3y_a = I_ESP_L2x3y3z_D2y_a+ABY*I_ESP_K2x2y3z_D2y_a;
    Double I_ESP_K2xy4z_F3y_a = I_ESP_L2x2y4z_D2y_a+ABY*I_ESP_K2xy4z_D2y_a;
    Double I_ESP_K2x5z_F3y_a = I_ESP_L2xy5z_D2y_a+ABY*I_ESP_K2x5z_D2y_a;
    Double I_ESP_Kx6y_F3y_a = I_ESP_Lx7y_D2y_a+ABY*I_ESP_Kx6y_D2y_a;
    Double I_ESP_Kx5yz_F3y_a = I_ESP_Lx6yz_D2y_a+ABY*I_ESP_Kx5yz_D2y_a;
    Double I_ESP_Kx4y2z_F3y_a = I_ESP_Lx5y2z_D2y_a+ABY*I_ESP_Kx4y2z_D2y_a;
    Double I_ESP_Kx3y3z_F3y_a = I_ESP_Lx4y3z_D2y_a+ABY*I_ESP_Kx3y3z_D2y_a;
    Double I_ESP_Kx2y4z_F3y_a = I_ESP_Lx3y4z_D2y_a+ABY*I_ESP_Kx2y4z_D2y_a;
    Double I_ESP_Kxy5z_F3y_a = I_ESP_Lx2y5z_D2y_a+ABY*I_ESP_Kxy5z_D2y_a;
    Double I_ESP_Kx6z_F3y_a = I_ESP_Lxy6z_D2y_a+ABY*I_ESP_Kx6z_D2y_a;
    Double I_ESP_K7y_F3y_a = I_ESP_L8y_D2y_a+ABY*I_ESP_K7y_D2y_a;
    Double I_ESP_K6yz_F3y_a = I_ESP_L7yz_D2y_a+ABY*I_ESP_K6yz_D2y_a;
    Double I_ESP_K5y2z_F3y_a = I_ESP_L6y2z_D2y_a+ABY*I_ESP_K5y2z_D2y_a;
    Double I_ESP_K4y3z_F3y_a = I_ESP_L5y3z_D2y_a+ABY*I_ESP_K4y3z_D2y_a;
    Double I_ESP_K3y4z_F3y_a = I_ESP_L4y4z_D2y_a+ABY*I_ESP_K3y4z_D2y_a;
    Double I_ESP_K2y5z_F3y_a = I_ESP_L3y5z_D2y_a+ABY*I_ESP_K2y5z_D2y_a;
    Double I_ESP_Ky6z_F3y_a = I_ESP_L2y6z_D2y_a+ABY*I_ESP_Ky6z_D2y_a;
    Double I_ESP_K7z_F3y_a = I_ESP_Ly7z_D2y_a+ABY*I_ESP_K7z_D2y_a;
    Double I_ESP_K6xz_F2yz_a = I_ESP_L6x2z_D2y_a+ABZ*I_ESP_K6xz_D2y_a;
    Double I_ESP_K5xyz_F2yz_a = I_ESP_L5xy2z_D2y_a+ABZ*I_ESP_K5xyz_D2y_a;
    Double I_ESP_K5x2z_F2yz_a = I_ESP_L5x3z_D2y_a+ABZ*I_ESP_K5x2z_D2y_a;
    Double I_ESP_K4x2yz_F2yz_a = I_ESP_L4x2y2z_D2y_a+ABZ*I_ESP_K4x2yz_D2y_a;
    Double I_ESP_K4xy2z_F2yz_a = I_ESP_L4xy3z_D2y_a+ABZ*I_ESP_K4xy2z_D2y_a;
    Double I_ESP_K4x3z_F2yz_a = I_ESP_L4x4z_D2y_a+ABZ*I_ESP_K4x3z_D2y_a;
    Double I_ESP_K3x3yz_F2yz_a = I_ESP_L3x3y2z_D2y_a+ABZ*I_ESP_K3x3yz_D2y_a;
    Double I_ESP_K3x2y2z_F2yz_a = I_ESP_L3x2y3z_D2y_a+ABZ*I_ESP_K3x2y2z_D2y_a;
    Double I_ESP_K3xy3z_F2yz_a = I_ESP_L3xy4z_D2y_a+ABZ*I_ESP_K3xy3z_D2y_a;
    Double I_ESP_K3x4z_F2yz_a = I_ESP_L3x5z_D2y_a+ABZ*I_ESP_K3x4z_D2y_a;
    Double I_ESP_K2x4yz_F2yz_a = I_ESP_L2x4y2z_D2y_a+ABZ*I_ESP_K2x4yz_D2y_a;
    Double I_ESP_K2x3y2z_F2yz_a = I_ESP_L2x3y3z_D2y_a+ABZ*I_ESP_K2x3y2z_D2y_a;
    Double I_ESP_K2x2y3z_F2yz_a = I_ESP_L2x2y4z_D2y_a+ABZ*I_ESP_K2x2y3z_D2y_a;
    Double I_ESP_K2xy4z_F2yz_a = I_ESP_L2xy5z_D2y_a+ABZ*I_ESP_K2xy4z_D2y_a;
    Double I_ESP_K2x5z_F2yz_a = I_ESP_L2x6z_D2y_a+ABZ*I_ESP_K2x5z_D2y_a;
    Double I_ESP_Kx5yz_F2yz_a = I_ESP_Lx5y2z_D2y_a+ABZ*I_ESP_Kx5yz_D2y_a;
    Double I_ESP_Kx4y2z_F2yz_a = I_ESP_Lx4y3z_D2y_a+ABZ*I_ESP_Kx4y2z_D2y_a;
    Double I_ESP_Kx3y3z_F2yz_a = I_ESP_Lx3y4z_D2y_a+ABZ*I_ESP_Kx3y3z_D2y_a;
    Double I_ESP_Kx2y4z_F2yz_a = I_ESP_Lx2y5z_D2y_a+ABZ*I_ESP_Kx2y4z_D2y_a;
    Double I_ESP_Kxy5z_F2yz_a = I_ESP_Lxy6z_D2y_a+ABZ*I_ESP_Kxy5z_D2y_a;
    Double I_ESP_Kx6z_F2yz_a = I_ESP_Lx7z_D2y_a+ABZ*I_ESP_Kx6z_D2y_a;
    Double I_ESP_K6yz_F2yz_a = I_ESP_L6y2z_D2y_a+ABZ*I_ESP_K6yz_D2y_a;
    Double I_ESP_K5y2z_F2yz_a = I_ESP_L5y3z_D2y_a+ABZ*I_ESP_K5y2z_D2y_a;
    Double I_ESP_K4y3z_F2yz_a = I_ESP_L4y4z_D2y_a+ABZ*I_ESP_K4y3z_D2y_a;
    Double I_ESP_K3y4z_F2yz_a = I_ESP_L3y5z_D2y_a+ABZ*I_ESP_K3y4z_D2y_a;
    Double I_ESP_K2y5z_F2yz_a = I_ESP_L2y6z_D2y_a+ABZ*I_ESP_K2y5z_D2y_a;
    Double I_ESP_Ky6z_F2yz_a = I_ESP_Ly7z_D2y_a+ABZ*I_ESP_Ky6z_D2y_a;
    Double I_ESP_K7z_F2yz_a = I_ESP_L8z_D2y_a+ABZ*I_ESP_K7z_D2y_a;
    Double I_ESP_K7x_F3z_a = I_ESP_L7xz_D2z_a+ABZ*I_ESP_K7x_D2z_a;
    Double I_ESP_K6xy_F3z_a = I_ESP_L6xyz_D2z_a+ABZ*I_ESP_K6xy_D2z_a;
    Double I_ESP_K6xz_F3z_a = I_ESP_L6x2z_D2z_a+ABZ*I_ESP_K6xz_D2z_a;
    Double I_ESP_K5x2y_F3z_a = I_ESP_L5x2yz_D2z_a+ABZ*I_ESP_K5x2y_D2z_a;
    Double I_ESP_K5xyz_F3z_a = I_ESP_L5xy2z_D2z_a+ABZ*I_ESP_K5xyz_D2z_a;
    Double I_ESP_K5x2z_F3z_a = I_ESP_L5x3z_D2z_a+ABZ*I_ESP_K5x2z_D2z_a;
    Double I_ESP_K4x3y_F3z_a = I_ESP_L4x3yz_D2z_a+ABZ*I_ESP_K4x3y_D2z_a;
    Double I_ESP_K4x2yz_F3z_a = I_ESP_L4x2y2z_D2z_a+ABZ*I_ESP_K4x2yz_D2z_a;
    Double I_ESP_K4xy2z_F3z_a = I_ESP_L4xy3z_D2z_a+ABZ*I_ESP_K4xy2z_D2z_a;
    Double I_ESP_K4x3z_F3z_a = I_ESP_L4x4z_D2z_a+ABZ*I_ESP_K4x3z_D2z_a;
    Double I_ESP_K3x4y_F3z_a = I_ESP_L3x4yz_D2z_a+ABZ*I_ESP_K3x4y_D2z_a;
    Double I_ESP_K3x3yz_F3z_a = I_ESP_L3x3y2z_D2z_a+ABZ*I_ESP_K3x3yz_D2z_a;
    Double I_ESP_K3x2y2z_F3z_a = I_ESP_L3x2y3z_D2z_a+ABZ*I_ESP_K3x2y2z_D2z_a;
    Double I_ESP_K3xy3z_F3z_a = I_ESP_L3xy4z_D2z_a+ABZ*I_ESP_K3xy3z_D2z_a;
    Double I_ESP_K3x4z_F3z_a = I_ESP_L3x5z_D2z_a+ABZ*I_ESP_K3x4z_D2z_a;
    Double I_ESP_K2x5y_F3z_a = I_ESP_L2x5yz_D2z_a+ABZ*I_ESP_K2x5y_D2z_a;
    Double I_ESP_K2x4yz_F3z_a = I_ESP_L2x4y2z_D2z_a+ABZ*I_ESP_K2x4yz_D2z_a;
    Double I_ESP_K2x3y2z_F3z_a = I_ESP_L2x3y3z_D2z_a+ABZ*I_ESP_K2x3y2z_D2z_a;
    Double I_ESP_K2x2y3z_F3z_a = I_ESP_L2x2y4z_D2z_a+ABZ*I_ESP_K2x2y3z_D2z_a;
    Double I_ESP_K2xy4z_F3z_a = I_ESP_L2xy5z_D2z_a+ABZ*I_ESP_K2xy4z_D2z_a;
    Double I_ESP_K2x5z_F3z_a = I_ESP_L2x6z_D2z_a+ABZ*I_ESP_K2x5z_D2z_a;
    Double I_ESP_Kx6y_F3z_a = I_ESP_Lx6yz_D2z_a+ABZ*I_ESP_Kx6y_D2z_a;
    Double I_ESP_Kx5yz_F3z_a = I_ESP_Lx5y2z_D2z_a+ABZ*I_ESP_Kx5yz_D2z_a;
    Double I_ESP_Kx4y2z_F3z_a = I_ESP_Lx4y3z_D2z_a+ABZ*I_ESP_Kx4y2z_D2z_a;
    Double I_ESP_Kx3y3z_F3z_a = I_ESP_Lx3y4z_D2z_a+ABZ*I_ESP_Kx3y3z_D2z_a;
    Double I_ESP_Kx2y4z_F3z_a = I_ESP_Lx2y5z_D2z_a+ABZ*I_ESP_Kx2y4z_D2z_a;
    Double I_ESP_Kxy5z_F3z_a = I_ESP_Lxy6z_D2z_a+ABZ*I_ESP_Kxy5z_D2z_a;
    Double I_ESP_Kx6z_F3z_a = I_ESP_Lx7z_D2z_a+ABZ*I_ESP_Kx6z_D2z_a;
    Double I_ESP_K7y_F3z_a = I_ESP_L7yz_D2z_a+ABZ*I_ESP_K7y_D2z_a;
    Double I_ESP_K6yz_F3z_a = I_ESP_L6y2z_D2z_a+ABZ*I_ESP_K6yz_D2z_a;
    Double I_ESP_K5y2z_F3z_a = I_ESP_L5y3z_D2z_a+ABZ*I_ESP_K5y2z_D2z_a;
    Double I_ESP_K4y3z_F3z_a = I_ESP_L4y4z_D2z_a+ABZ*I_ESP_K4y3z_D2z_a;
    Double I_ESP_K3y4z_F3z_a = I_ESP_L3y5z_D2z_a+ABZ*I_ESP_K3y4z_D2z_a;
    Double I_ESP_K2y5z_F3z_a = I_ESP_L2y6z_D2z_a+ABZ*I_ESP_K2y5z_D2z_a;
    Double I_ESP_Ky6z_F3z_a = I_ESP_Ly7z_D2z_a+ABZ*I_ESP_Ky6z_D2z_a;
    Double I_ESP_K7z_F3z_a = I_ESP_L8z_D2z_a+ABZ*I_ESP_K7z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_G_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_F_a
     * RHS shell quartet name: SQ_ESP_I_F_a
     ************************************************************/
    Double I_ESP_I6x_G4x_a = I_ESP_K7x_F3x_a+ABX*I_ESP_I6x_F3x_a;
    Double I_ESP_I5xy_G4x_a = I_ESP_K6xy_F3x_a+ABX*I_ESP_I5xy_F3x_a;
    Double I_ESP_I5xz_G4x_a = I_ESP_K6xz_F3x_a+ABX*I_ESP_I5xz_F3x_a;
    Double I_ESP_I4x2y_G4x_a = I_ESP_K5x2y_F3x_a+ABX*I_ESP_I4x2y_F3x_a;
    Double I_ESP_I4xyz_G4x_a = I_ESP_K5xyz_F3x_a+ABX*I_ESP_I4xyz_F3x_a;
    Double I_ESP_I4x2z_G4x_a = I_ESP_K5x2z_F3x_a+ABX*I_ESP_I4x2z_F3x_a;
    Double I_ESP_I3x3y_G4x_a = I_ESP_K4x3y_F3x_a+ABX*I_ESP_I3x3y_F3x_a;
    Double I_ESP_I3x2yz_G4x_a = I_ESP_K4x2yz_F3x_a+ABX*I_ESP_I3x2yz_F3x_a;
    Double I_ESP_I3xy2z_G4x_a = I_ESP_K4xy2z_F3x_a+ABX*I_ESP_I3xy2z_F3x_a;
    Double I_ESP_I3x3z_G4x_a = I_ESP_K4x3z_F3x_a+ABX*I_ESP_I3x3z_F3x_a;
    Double I_ESP_I2x4y_G4x_a = I_ESP_K3x4y_F3x_a+ABX*I_ESP_I2x4y_F3x_a;
    Double I_ESP_I2x3yz_G4x_a = I_ESP_K3x3yz_F3x_a+ABX*I_ESP_I2x3yz_F3x_a;
    Double I_ESP_I2x2y2z_G4x_a = I_ESP_K3x2y2z_F3x_a+ABX*I_ESP_I2x2y2z_F3x_a;
    Double I_ESP_I2xy3z_G4x_a = I_ESP_K3xy3z_F3x_a+ABX*I_ESP_I2xy3z_F3x_a;
    Double I_ESP_I2x4z_G4x_a = I_ESP_K3x4z_F3x_a+ABX*I_ESP_I2x4z_F3x_a;
    Double I_ESP_Ix5y_G4x_a = I_ESP_K2x5y_F3x_a+ABX*I_ESP_Ix5y_F3x_a;
    Double I_ESP_Ix4yz_G4x_a = I_ESP_K2x4yz_F3x_a+ABX*I_ESP_Ix4yz_F3x_a;
    Double I_ESP_Ix3y2z_G4x_a = I_ESP_K2x3y2z_F3x_a+ABX*I_ESP_Ix3y2z_F3x_a;
    Double I_ESP_Ix2y3z_G4x_a = I_ESP_K2x2y3z_F3x_a+ABX*I_ESP_Ix2y3z_F3x_a;
    Double I_ESP_Ixy4z_G4x_a = I_ESP_K2xy4z_F3x_a+ABX*I_ESP_Ixy4z_F3x_a;
    Double I_ESP_Ix5z_G4x_a = I_ESP_K2x5z_F3x_a+ABX*I_ESP_Ix5z_F3x_a;
    Double I_ESP_I6y_G4x_a = I_ESP_Kx6y_F3x_a+ABX*I_ESP_I6y_F3x_a;
    Double I_ESP_I5yz_G4x_a = I_ESP_Kx5yz_F3x_a+ABX*I_ESP_I5yz_F3x_a;
    Double I_ESP_I4y2z_G4x_a = I_ESP_Kx4y2z_F3x_a+ABX*I_ESP_I4y2z_F3x_a;
    Double I_ESP_I3y3z_G4x_a = I_ESP_Kx3y3z_F3x_a+ABX*I_ESP_I3y3z_F3x_a;
    Double I_ESP_I2y4z_G4x_a = I_ESP_Kx2y4z_F3x_a+ABX*I_ESP_I2y4z_F3x_a;
    Double I_ESP_Iy5z_G4x_a = I_ESP_Kxy5z_F3x_a+ABX*I_ESP_Iy5z_F3x_a;
    Double I_ESP_I6z_G4x_a = I_ESP_Kx6z_F3x_a+ABX*I_ESP_I6z_F3x_a;
    Double I_ESP_I6x_G3xy_a = I_ESP_K6xy_F3x_a+ABY*I_ESP_I6x_F3x_a;
    Double I_ESP_I5xy_G3xy_a = I_ESP_K5x2y_F3x_a+ABY*I_ESP_I5xy_F3x_a;
    Double I_ESP_I5xz_G3xy_a = I_ESP_K5xyz_F3x_a+ABY*I_ESP_I5xz_F3x_a;
    Double I_ESP_I4x2y_G3xy_a = I_ESP_K4x3y_F3x_a+ABY*I_ESP_I4x2y_F3x_a;
    Double I_ESP_I4xyz_G3xy_a = I_ESP_K4x2yz_F3x_a+ABY*I_ESP_I4xyz_F3x_a;
    Double I_ESP_I4x2z_G3xy_a = I_ESP_K4xy2z_F3x_a+ABY*I_ESP_I4x2z_F3x_a;
    Double I_ESP_I3x3y_G3xy_a = I_ESP_K3x4y_F3x_a+ABY*I_ESP_I3x3y_F3x_a;
    Double I_ESP_I3x2yz_G3xy_a = I_ESP_K3x3yz_F3x_a+ABY*I_ESP_I3x2yz_F3x_a;
    Double I_ESP_I3xy2z_G3xy_a = I_ESP_K3x2y2z_F3x_a+ABY*I_ESP_I3xy2z_F3x_a;
    Double I_ESP_I3x3z_G3xy_a = I_ESP_K3xy3z_F3x_a+ABY*I_ESP_I3x3z_F3x_a;
    Double I_ESP_I2x4y_G3xy_a = I_ESP_K2x5y_F3x_a+ABY*I_ESP_I2x4y_F3x_a;
    Double I_ESP_I2x3yz_G3xy_a = I_ESP_K2x4yz_F3x_a+ABY*I_ESP_I2x3yz_F3x_a;
    Double I_ESP_I2x2y2z_G3xy_a = I_ESP_K2x3y2z_F3x_a+ABY*I_ESP_I2x2y2z_F3x_a;
    Double I_ESP_I2xy3z_G3xy_a = I_ESP_K2x2y3z_F3x_a+ABY*I_ESP_I2xy3z_F3x_a;
    Double I_ESP_I2x4z_G3xy_a = I_ESP_K2xy4z_F3x_a+ABY*I_ESP_I2x4z_F3x_a;
    Double I_ESP_Ix5y_G3xy_a = I_ESP_Kx6y_F3x_a+ABY*I_ESP_Ix5y_F3x_a;
    Double I_ESP_Ix4yz_G3xy_a = I_ESP_Kx5yz_F3x_a+ABY*I_ESP_Ix4yz_F3x_a;
    Double I_ESP_Ix3y2z_G3xy_a = I_ESP_Kx4y2z_F3x_a+ABY*I_ESP_Ix3y2z_F3x_a;
    Double I_ESP_Ix2y3z_G3xy_a = I_ESP_Kx3y3z_F3x_a+ABY*I_ESP_Ix2y3z_F3x_a;
    Double I_ESP_Ixy4z_G3xy_a = I_ESP_Kx2y4z_F3x_a+ABY*I_ESP_Ixy4z_F3x_a;
    Double I_ESP_Ix5z_G3xy_a = I_ESP_Kxy5z_F3x_a+ABY*I_ESP_Ix5z_F3x_a;
    Double I_ESP_I6y_G3xy_a = I_ESP_K7y_F3x_a+ABY*I_ESP_I6y_F3x_a;
    Double I_ESP_I5yz_G3xy_a = I_ESP_K6yz_F3x_a+ABY*I_ESP_I5yz_F3x_a;
    Double I_ESP_I4y2z_G3xy_a = I_ESP_K5y2z_F3x_a+ABY*I_ESP_I4y2z_F3x_a;
    Double I_ESP_I3y3z_G3xy_a = I_ESP_K4y3z_F3x_a+ABY*I_ESP_I3y3z_F3x_a;
    Double I_ESP_I2y4z_G3xy_a = I_ESP_K3y4z_F3x_a+ABY*I_ESP_I2y4z_F3x_a;
    Double I_ESP_Iy5z_G3xy_a = I_ESP_K2y5z_F3x_a+ABY*I_ESP_Iy5z_F3x_a;
    Double I_ESP_I6z_G3xy_a = I_ESP_Ky6z_F3x_a+ABY*I_ESP_I6z_F3x_a;
    Double I_ESP_I6x_G3xz_a = I_ESP_K6xz_F3x_a+ABZ*I_ESP_I6x_F3x_a;
    Double I_ESP_I5xy_G3xz_a = I_ESP_K5xyz_F3x_a+ABZ*I_ESP_I5xy_F3x_a;
    Double I_ESP_I5xz_G3xz_a = I_ESP_K5x2z_F3x_a+ABZ*I_ESP_I5xz_F3x_a;
    Double I_ESP_I4x2y_G3xz_a = I_ESP_K4x2yz_F3x_a+ABZ*I_ESP_I4x2y_F3x_a;
    Double I_ESP_I4xyz_G3xz_a = I_ESP_K4xy2z_F3x_a+ABZ*I_ESP_I4xyz_F3x_a;
    Double I_ESP_I4x2z_G3xz_a = I_ESP_K4x3z_F3x_a+ABZ*I_ESP_I4x2z_F3x_a;
    Double I_ESP_I3x3y_G3xz_a = I_ESP_K3x3yz_F3x_a+ABZ*I_ESP_I3x3y_F3x_a;
    Double I_ESP_I3x2yz_G3xz_a = I_ESP_K3x2y2z_F3x_a+ABZ*I_ESP_I3x2yz_F3x_a;
    Double I_ESP_I3xy2z_G3xz_a = I_ESP_K3xy3z_F3x_a+ABZ*I_ESP_I3xy2z_F3x_a;
    Double I_ESP_I3x3z_G3xz_a = I_ESP_K3x4z_F3x_a+ABZ*I_ESP_I3x3z_F3x_a;
    Double I_ESP_I2x4y_G3xz_a = I_ESP_K2x4yz_F3x_a+ABZ*I_ESP_I2x4y_F3x_a;
    Double I_ESP_I2x3yz_G3xz_a = I_ESP_K2x3y2z_F3x_a+ABZ*I_ESP_I2x3yz_F3x_a;
    Double I_ESP_I2x2y2z_G3xz_a = I_ESP_K2x2y3z_F3x_a+ABZ*I_ESP_I2x2y2z_F3x_a;
    Double I_ESP_I2xy3z_G3xz_a = I_ESP_K2xy4z_F3x_a+ABZ*I_ESP_I2xy3z_F3x_a;
    Double I_ESP_I2x4z_G3xz_a = I_ESP_K2x5z_F3x_a+ABZ*I_ESP_I2x4z_F3x_a;
    Double I_ESP_Ix5y_G3xz_a = I_ESP_Kx5yz_F3x_a+ABZ*I_ESP_Ix5y_F3x_a;
    Double I_ESP_Ix4yz_G3xz_a = I_ESP_Kx4y2z_F3x_a+ABZ*I_ESP_Ix4yz_F3x_a;
    Double I_ESP_Ix3y2z_G3xz_a = I_ESP_Kx3y3z_F3x_a+ABZ*I_ESP_Ix3y2z_F3x_a;
    Double I_ESP_Ix2y3z_G3xz_a = I_ESP_Kx2y4z_F3x_a+ABZ*I_ESP_Ix2y3z_F3x_a;
    Double I_ESP_Ixy4z_G3xz_a = I_ESP_Kxy5z_F3x_a+ABZ*I_ESP_Ixy4z_F3x_a;
    Double I_ESP_Ix5z_G3xz_a = I_ESP_Kx6z_F3x_a+ABZ*I_ESP_Ix5z_F3x_a;
    Double I_ESP_I6y_G3xz_a = I_ESP_K6yz_F3x_a+ABZ*I_ESP_I6y_F3x_a;
    Double I_ESP_I5yz_G3xz_a = I_ESP_K5y2z_F3x_a+ABZ*I_ESP_I5yz_F3x_a;
    Double I_ESP_I4y2z_G3xz_a = I_ESP_K4y3z_F3x_a+ABZ*I_ESP_I4y2z_F3x_a;
    Double I_ESP_I3y3z_G3xz_a = I_ESP_K3y4z_F3x_a+ABZ*I_ESP_I3y3z_F3x_a;
    Double I_ESP_I2y4z_G3xz_a = I_ESP_K2y5z_F3x_a+ABZ*I_ESP_I2y4z_F3x_a;
    Double I_ESP_Iy5z_G3xz_a = I_ESP_Ky6z_F3x_a+ABZ*I_ESP_Iy5z_F3x_a;
    Double I_ESP_I6z_G3xz_a = I_ESP_K7z_F3x_a+ABZ*I_ESP_I6z_F3x_a;
    Double I_ESP_I6x_G2x2y_a = I_ESP_K6xy_F2xy_a+ABY*I_ESP_I6x_F2xy_a;
    Double I_ESP_I5xy_G2x2y_a = I_ESP_K5x2y_F2xy_a+ABY*I_ESP_I5xy_F2xy_a;
    Double I_ESP_I5xz_G2x2y_a = I_ESP_K5xyz_F2xy_a+ABY*I_ESP_I5xz_F2xy_a;
    Double I_ESP_I4x2y_G2x2y_a = I_ESP_K4x3y_F2xy_a+ABY*I_ESP_I4x2y_F2xy_a;
    Double I_ESP_I4xyz_G2x2y_a = I_ESP_K4x2yz_F2xy_a+ABY*I_ESP_I4xyz_F2xy_a;
    Double I_ESP_I4x2z_G2x2y_a = I_ESP_K4xy2z_F2xy_a+ABY*I_ESP_I4x2z_F2xy_a;
    Double I_ESP_I3x3y_G2x2y_a = I_ESP_K3x4y_F2xy_a+ABY*I_ESP_I3x3y_F2xy_a;
    Double I_ESP_I3x2yz_G2x2y_a = I_ESP_K3x3yz_F2xy_a+ABY*I_ESP_I3x2yz_F2xy_a;
    Double I_ESP_I3xy2z_G2x2y_a = I_ESP_K3x2y2z_F2xy_a+ABY*I_ESP_I3xy2z_F2xy_a;
    Double I_ESP_I3x3z_G2x2y_a = I_ESP_K3xy3z_F2xy_a+ABY*I_ESP_I3x3z_F2xy_a;
    Double I_ESP_I2x4y_G2x2y_a = I_ESP_K2x5y_F2xy_a+ABY*I_ESP_I2x4y_F2xy_a;
    Double I_ESP_I2x3yz_G2x2y_a = I_ESP_K2x4yz_F2xy_a+ABY*I_ESP_I2x3yz_F2xy_a;
    Double I_ESP_I2x2y2z_G2x2y_a = I_ESP_K2x3y2z_F2xy_a+ABY*I_ESP_I2x2y2z_F2xy_a;
    Double I_ESP_I2xy3z_G2x2y_a = I_ESP_K2x2y3z_F2xy_a+ABY*I_ESP_I2xy3z_F2xy_a;
    Double I_ESP_I2x4z_G2x2y_a = I_ESP_K2xy4z_F2xy_a+ABY*I_ESP_I2x4z_F2xy_a;
    Double I_ESP_Ix5y_G2x2y_a = I_ESP_Kx6y_F2xy_a+ABY*I_ESP_Ix5y_F2xy_a;
    Double I_ESP_Ix4yz_G2x2y_a = I_ESP_Kx5yz_F2xy_a+ABY*I_ESP_Ix4yz_F2xy_a;
    Double I_ESP_Ix3y2z_G2x2y_a = I_ESP_Kx4y2z_F2xy_a+ABY*I_ESP_Ix3y2z_F2xy_a;
    Double I_ESP_Ix2y3z_G2x2y_a = I_ESP_Kx3y3z_F2xy_a+ABY*I_ESP_Ix2y3z_F2xy_a;
    Double I_ESP_Ixy4z_G2x2y_a = I_ESP_Kx2y4z_F2xy_a+ABY*I_ESP_Ixy4z_F2xy_a;
    Double I_ESP_Ix5z_G2x2y_a = I_ESP_Kxy5z_F2xy_a+ABY*I_ESP_Ix5z_F2xy_a;
    Double I_ESP_I6y_G2x2y_a = I_ESP_K7y_F2xy_a+ABY*I_ESP_I6y_F2xy_a;
    Double I_ESP_I5yz_G2x2y_a = I_ESP_K6yz_F2xy_a+ABY*I_ESP_I5yz_F2xy_a;
    Double I_ESP_I4y2z_G2x2y_a = I_ESP_K5y2z_F2xy_a+ABY*I_ESP_I4y2z_F2xy_a;
    Double I_ESP_I3y3z_G2x2y_a = I_ESP_K4y3z_F2xy_a+ABY*I_ESP_I3y3z_F2xy_a;
    Double I_ESP_I2y4z_G2x2y_a = I_ESP_K3y4z_F2xy_a+ABY*I_ESP_I2y4z_F2xy_a;
    Double I_ESP_Iy5z_G2x2y_a = I_ESP_K2y5z_F2xy_a+ABY*I_ESP_Iy5z_F2xy_a;
    Double I_ESP_I6z_G2x2y_a = I_ESP_Ky6z_F2xy_a+ABY*I_ESP_I6z_F2xy_a;
    Double I_ESP_I6x_G2xyz_a = I_ESP_K6xz_F2xy_a+ABZ*I_ESP_I6x_F2xy_a;
    Double I_ESP_I5xy_G2xyz_a = I_ESP_K5xyz_F2xy_a+ABZ*I_ESP_I5xy_F2xy_a;
    Double I_ESP_I5xz_G2xyz_a = I_ESP_K5x2z_F2xy_a+ABZ*I_ESP_I5xz_F2xy_a;
    Double I_ESP_I4x2y_G2xyz_a = I_ESP_K4x2yz_F2xy_a+ABZ*I_ESP_I4x2y_F2xy_a;
    Double I_ESP_I4xyz_G2xyz_a = I_ESP_K4xy2z_F2xy_a+ABZ*I_ESP_I4xyz_F2xy_a;
    Double I_ESP_I4x2z_G2xyz_a = I_ESP_K4x3z_F2xy_a+ABZ*I_ESP_I4x2z_F2xy_a;
    Double I_ESP_I3x3y_G2xyz_a = I_ESP_K3x3yz_F2xy_a+ABZ*I_ESP_I3x3y_F2xy_a;
    Double I_ESP_I3x2yz_G2xyz_a = I_ESP_K3x2y2z_F2xy_a+ABZ*I_ESP_I3x2yz_F2xy_a;
    Double I_ESP_I3xy2z_G2xyz_a = I_ESP_K3xy3z_F2xy_a+ABZ*I_ESP_I3xy2z_F2xy_a;
    Double I_ESP_I3x3z_G2xyz_a = I_ESP_K3x4z_F2xy_a+ABZ*I_ESP_I3x3z_F2xy_a;
    Double I_ESP_I2x4y_G2xyz_a = I_ESP_K2x4yz_F2xy_a+ABZ*I_ESP_I2x4y_F2xy_a;
    Double I_ESP_I2x3yz_G2xyz_a = I_ESP_K2x3y2z_F2xy_a+ABZ*I_ESP_I2x3yz_F2xy_a;
    Double I_ESP_I2x2y2z_G2xyz_a = I_ESP_K2x2y3z_F2xy_a+ABZ*I_ESP_I2x2y2z_F2xy_a;
    Double I_ESP_I2xy3z_G2xyz_a = I_ESP_K2xy4z_F2xy_a+ABZ*I_ESP_I2xy3z_F2xy_a;
    Double I_ESP_I2x4z_G2xyz_a = I_ESP_K2x5z_F2xy_a+ABZ*I_ESP_I2x4z_F2xy_a;
    Double I_ESP_Ix5y_G2xyz_a = I_ESP_Kx5yz_F2xy_a+ABZ*I_ESP_Ix5y_F2xy_a;
    Double I_ESP_Ix4yz_G2xyz_a = I_ESP_Kx4y2z_F2xy_a+ABZ*I_ESP_Ix4yz_F2xy_a;
    Double I_ESP_Ix3y2z_G2xyz_a = I_ESP_Kx3y3z_F2xy_a+ABZ*I_ESP_Ix3y2z_F2xy_a;
    Double I_ESP_Ix2y3z_G2xyz_a = I_ESP_Kx2y4z_F2xy_a+ABZ*I_ESP_Ix2y3z_F2xy_a;
    Double I_ESP_Ixy4z_G2xyz_a = I_ESP_Kxy5z_F2xy_a+ABZ*I_ESP_Ixy4z_F2xy_a;
    Double I_ESP_Ix5z_G2xyz_a = I_ESP_Kx6z_F2xy_a+ABZ*I_ESP_Ix5z_F2xy_a;
    Double I_ESP_I6y_G2xyz_a = I_ESP_K6yz_F2xy_a+ABZ*I_ESP_I6y_F2xy_a;
    Double I_ESP_I5yz_G2xyz_a = I_ESP_K5y2z_F2xy_a+ABZ*I_ESP_I5yz_F2xy_a;
    Double I_ESP_I4y2z_G2xyz_a = I_ESP_K4y3z_F2xy_a+ABZ*I_ESP_I4y2z_F2xy_a;
    Double I_ESP_I3y3z_G2xyz_a = I_ESP_K3y4z_F2xy_a+ABZ*I_ESP_I3y3z_F2xy_a;
    Double I_ESP_I2y4z_G2xyz_a = I_ESP_K2y5z_F2xy_a+ABZ*I_ESP_I2y4z_F2xy_a;
    Double I_ESP_Iy5z_G2xyz_a = I_ESP_Ky6z_F2xy_a+ABZ*I_ESP_Iy5z_F2xy_a;
    Double I_ESP_I6z_G2xyz_a = I_ESP_K7z_F2xy_a+ABZ*I_ESP_I6z_F2xy_a;
    Double I_ESP_I6x_G2x2z_a = I_ESP_K6xz_F2xz_a+ABZ*I_ESP_I6x_F2xz_a;
    Double I_ESP_I5xy_G2x2z_a = I_ESP_K5xyz_F2xz_a+ABZ*I_ESP_I5xy_F2xz_a;
    Double I_ESP_I5xz_G2x2z_a = I_ESP_K5x2z_F2xz_a+ABZ*I_ESP_I5xz_F2xz_a;
    Double I_ESP_I4x2y_G2x2z_a = I_ESP_K4x2yz_F2xz_a+ABZ*I_ESP_I4x2y_F2xz_a;
    Double I_ESP_I4xyz_G2x2z_a = I_ESP_K4xy2z_F2xz_a+ABZ*I_ESP_I4xyz_F2xz_a;
    Double I_ESP_I4x2z_G2x2z_a = I_ESP_K4x3z_F2xz_a+ABZ*I_ESP_I4x2z_F2xz_a;
    Double I_ESP_I3x3y_G2x2z_a = I_ESP_K3x3yz_F2xz_a+ABZ*I_ESP_I3x3y_F2xz_a;
    Double I_ESP_I3x2yz_G2x2z_a = I_ESP_K3x2y2z_F2xz_a+ABZ*I_ESP_I3x2yz_F2xz_a;
    Double I_ESP_I3xy2z_G2x2z_a = I_ESP_K3xy3z_F2xz_a+ABZ*I_ESP_I3xy2z_F2xz_a;
    Double I_ESP_I3x3z_G2x2z_a = I_ESP_K3x4z_F2xz_a+ABZ*I_ESP_I3x3z_F2xz_a;
    Double I_ESP_I2x4y_G2x2z_a = I_ESP_K2x4yz_F2xz_a+ABZ*I_ESP_I2x4y_F2xz_a;
    Double I_ESP_I2x3yz_G2x2z_a = I_ESP_K2x3y2z_F2xz_a+ABZ*I_ESP_I2x3yz_F2xz_a;
    Double I_ESP_I2x2y2z_G2x2z_a = I_ESP_K2x2y3z_F2xz_a+ABZ*I_ESP_I2x2y2z_F2xz_a;
    Double I_ESP_I2xy3z_G2x2z_a = I_ESP_K2xy4z_F2xz_a+ABZ*I_ESP_I2xy3z_F2xz_a;
    Double I_ESP_I2x4z_G2x2z_a = I_ESP_K2x5z_F2xz_a+ABZ*I_ESP_I2x4z_F2xz_a;
    Double I_ESP_Ix5y_G2x2z_a = I_ESP_Kx5yz_F2xz_a+ABZ*I_ESP_Ix5y_F2xz_a;
    Double I_ESP_Ix4yz_G2x2z_a = I_ESP_Kx4y2z_F2xz_a+ABZ*I_ESP_Ix4yz_F2xz_a;
    Double I_ESP_Ix3y2z_G2x2z_a = I_ESP_Kx3y3z_F2xz_a+ABZ*I_ESP_Ix3y2z_F2xz_a;
    Double I_ESP_Ix2y3z_G2x2z_a = I_ESP_Kx2y4z_F2xz_a+ABZ*I_ESP_Ix2y3z_F2xz_a;
    Double I_ESP_Ixy4z_G2x2z_a = I_ESP_Kxy5z_F2xz_a+ABZ*I_ESP_Ixy4z_F2xz_a;
    Double I_ESP_Ix5z_G2x2z_a = I_ESP_Kx6z_F2xz_a+ABZ*I_ESP_Ix5z_F2xz_a;
    Double I_ESP_I6y_G2x2z_a = I_ESP_K6yz_F2xz_a+ABZ*I_ESP_I6y_F2xz_a;
    Double I_ESP_I5yz_G2x2z_a = I_ESP_K5y2z_F2xz_a+ABZ*I_ESP_I5yz_F2xz_a;
    Double I_ESP_I4y2z_G2x2z_a = I_ESP_K4y3z_F2xz_a+ABZ*I_ESP_I4y2z_F2xz_a;
    Double I_ESP_I3y3z_G2x2z_a = I_ESP_K3y4z_F2xz_a+ABZ*I_ESP_I3y3z_F2xz_a;
    Double I_ESP_I2y4z_G2x2z_a = I_ESP_K2y5z_F2xz_a+ABZ*I_ESP_I2y4z_F2xz_a;
    Double I_ESP_Iy5z_G2x2z_a = I_ESP_Ky6z_F2xz_a+ABZ*I_ESP_Iy5z_F2xz_a;
    Double I_ESP_I6z_G2x2z_a = I_ESP_K7z_F2xz_a+ABZ*I_ESP_I6z_F2xz_a;
    Double I_ESP_I6x_Gx3y_a = I_ESP_K7x_F3y_a+ABX*I_ESP_I6x_F3y_a;
    Double I_ESP_I5xy_Gx3y_a = I_ESP_K6xy_F3y_a+ABX*I_ESP_I5xy_F3y_a;
    Double I_ESP_I5xz_Gx3y_a = I_ESP_K6xz_F3y_a+ABX*I_ESP_I5xz_F3y_a;
    Double I_ESP_I4x2y_Gx3y_a = I_ESP_K5x2y_F3y_a+ABX*I_ESP_I4x2y_F3y_a;
    Double I_ESP_I4xyz_Gx3y_a = I_ESP_K5xyz_F3y_a+ABX*I_ESP_I4xyz_F3y_a;
    Double I_ESP_I4x2z_Gx3y_a = I_ESP_K5x2z_F3y_a+ABX*I_ESP_I4x2z_F3y_a;
    Double I_ESP_I3x3y_Gx3y_a = I_ESP_K4x3y_F3y_a+ABX*I_ESP_I3x3y_F3y_a;
    Double I_ESP_I3x2yz_Gx3y_a = I_ESP_K4x2yz_F3y_a+ABX*I_ESP_I3x2yz_F3y_a;
    Double I_ESP_I3xy2z_Gx3y_a = I_ESP_K4xy2z_F3y_a+ABX*I_ESP_I3xy2z_F3y_a;
    Double I_ESP_I3x3z_Gx3y_a = I_ESP_K4x3z_F3y_a+ABX*I_ESP_I3x3z_F3y_a;
    Double I_ESP_I2x4y_Gx3y_a = I_ESP_K3x4y_F3y_a+ABX*I_ESP_I2x4y_F3y_a;
    Double I_ESP_I2x3yz_Gx3y_a = I_ESP_K3x3yz_F3y_a+ABX*I_ESP_I2x3yz_F3y_a;
    Double I_ESP_I2x2y2z_Gx3y_a = I_ESP_K3x2y2z_F3y_a+ABX*I_ESP_I2x2y2z_F3y_a;
    Double I_ESP_I2xy3z_Gx3y_a = I_ESP_K3xy3z_F3y_a+ABX*I_ESP_I2xy3z_F3y_a;
    Double I_ESP_I2x4z_Gx3y_a = I_ESP_K3x4z_F3y_a+ABX*I_ESP_I2x4z_F3y_a;
    Double I_ESP_Ix5y_Gx3y_a = I_ESP_K2x5y_F3y_a+ABX*I_ESP_Ix5y_F3y_a;
    Double I_ESP_Ix4yz_Gx3y_a = I_ESP_K2x4yz_F3y_a+ABX*I_ESP_Ix4yz_F3y_a;
    Double I_ESP_Ix3y2z_Gx3y_a = I_ESP_K2x3y2z_F3y_a+ABX*I_ESP_Ix3y2z_F3y_a;
    Double I_ESP_Ix2y3z_Gx3y_a = I_ESP_K2x2y3z_F3y_a+ABX*I_ESP_Ix2y3z_F3y_a;
    Double I_ESP_Ixy4z_Gx3y_a = I_ESP_K2xy4z_F3y_a+ABX*I_ESP_Ixy4z_F3y_a;
    Double I_ESP_Ix5z_Gx3y_a = I_ESP_K2x5z_F3y_a+ABX*I_ESP_Ix5z_F3y_a;
    Double I_ESP_I6y_Gx3y_a = I_ESP_Kx6y_F3y_a+ABX*I_ESP_I6y_F3y_a;
    Double I_ESP_I5yz_Gx3y_a = I_ESP_Kx5yz_F3y_a+ABX*I_ESP_I5yz_F3y_a;
    Double I_ESP_I4y2z_Gx3y_a = I_ESP_Kx4y2z_F3y_a+ABX*I_ESP_I4y2z_F3y_a;
    Double I_ESP_I3y3z_Gx3y_a = I_ESP_Kx3y3z_F3y_a+ABX*I_ESP_I3y3z_F3y_a;
    Double I_ESP_I2y4z_Gx3y_a = I_ESP_Kx2y4z_F3y_a+ABX*I_ESP_I2y4z_F3y_a;
    Double I_ESP_Iy5z_Gx3y_a = I_ESP_Kxy5z_F3y_a+ABX*I_ESP_Iy5z_F3y_a;
    Double I_ESP_I6z_Gx3y_a = I_ESP_Kx6z_F3y_a+ABX*I_ESP_I6z_F3y_a;
    Double I_ESP_I6x_Gx2yz_a = I_ESP_K6xz_Fx2y_a+ABZ*I_ESP_I6x_Fx2y_a;
    Double I_ESP_I5xy_Gx2yz_a = I_ESP_K5xyz_Fx2y_a+ABZ*I_ESP_I5xy_Fx2y_a;
    Double I_ESP_I5xz_Gx2yz_a = I_ESP_K5x2z_Fx2y_a+ABZ*I_ESP_I5xz_Fx2y_a;
    Double I_ESP_I4x2y_Gx2yz_a = I_ESP_K4x2yz_Fx2y_a+ABZ*I_ESP_I4x2y_Fx2y_a;
    Double I_ESP_I4xyz_Gx2yz_a = I_ESP_K4xy2z_Fx2y_a+ABZ*I_ESP_I4xyz_Fx2y_a;
    Double I_ESP_I4x2z_Gx2yz_a = I_ESP_K4x3z_Fx2y_a+ABZ*I_ESP_I4x2z_Fx2y_a;
    Double I_ESP_I3x3y_Gx2yz_a = I_ESP_K3x3yz_Fx2y_a+ABZ*I_ESP_I3x3y_Fx2y_a;
    Double I_ESP_I3x2yz_Gx2yz_a = I_ESP_K3x2y2z_Fx2y_a+ABZ*I_ESP_I3x2yz_Fx2y_a;
    Double I_ESP_I3xy2z_Gx2yz_a = I_ESP_K3xy3z_Fx2y_a+ABZ*I_ESP_I3xy2z_Fx2y_a;
    Double I_ESP_I3x3z_Gx2yz_a = I_ESP_K3x4z_Fx2y_a+ABZ*I_ESP_I3x3z_Fx2y_a;
    Double I_ESP_I2x4y_Gx2yz_a = I_ESP_K2x4yz_Fx2y_a+ABZ*I_ESP_I2x4y_Fx2y_a;
    Double I_ESP_I2x3yz_Gx2yz_a = I_ESP_K2x3y2z_Fx2y_a+ABZ*I_ESP_I2x3yz_Fx2y_a;
    Double I_ESP_I2x2y2z_Gx2yz_a = I_ESP_K2x2y3z_Fx2y_a+ABZ*I_ESP_I2x2y2z_Fx2y_a;
    Double I_ESP_I2xy3z_Gx2yz_a = I_ESP_K2xy4z_Fx2y_a+ABZ*I_ESP_I2xy3z_Fx2y_a;
    Double I_ESP_I2x4z_Gx2yz_a = I_ESP_K2x5z_Fx2y_a+ABZ*I_ESP_I2x4z_Fx2y_a;
    Double I_ESP_Ix5y_Gx2yz_a = I_ESP_Kx5yz_Fx2y_a+ABZ*I_ESP_Ix5y_Fx2y_a;
    Double I_ESP_Ix4yz_Gx2yz_a = I_ESP_Kx4y2z_Fx2y_a+ABZ*I_ESP_Ix4yz_Fx2y_a;
    Double I_ESP_Ix3y2z_Gx2yz_a = I_ESP_Kx3y3z_Fx2y_a+ABZ*I_ESP_Ix3y2z_Fx2y_a;
    Double I_ESP_Ix2y3z_Gx2yz_a = I_ESP_Kx2y4z_Fx2y_a+ABZ*I_ESP_Ix2y3z_Fx2y_a;
    Double I_ESP_Ixy4z_Gx2yz_a = I_ESP_Kxy5z_Fx2y_a+ABZ*I_ESP_Ixy4z_Fx2y_a;
    Double I_ESP_Ix5z_Gx2yz_a = I_ESP_Kx6z_Fx2y_a+ABZ*I_ESP_Ix5z_Fx2y_a;
    Double I_ESP_I6y_Gx2yz_a = I_ESP_K6yz_Fx2y_a+ABZ*I_ESP_I6y_Fx2y_a;
    Double I_ESP_I5yz_Gx2yz_a = I_ESP_K5y2z_Fx2y_a+ABZ*I_ESP_I5yz_Fx2y_a;
    Double I_ESP_I4y2z_Gx2yz_a = I_ESP_K4y3z_Fx2y_a+ABZ*I_ESP_I4y2z_Fx2y_a;
    Double I_ESP_I3y3z_Gx2yz_a = I_ESP_K3y4z_Fx2y_a+ABZ*I_ESP_I3y3z_Fx2y_a;
    Double I_ESP_I2y4z_Gx2yz_a = I_ESP_K2y5z_Fx2y_a+ABZ*I_ESP_I2y4z_Fx2y_a;
    Double I_ESP_Iy5z_Gx2yz_a = I_ESP_Ky6z_Fx2y_a+ABZ*I_ESP_Iy5z_Fx2y_a;
    Double I_ESP_I6z_Gx2yz_a = I_ESP_K7z_Fx2y_a+ABZ*I_ESP_I6z_Fx2y_a;
    Double I_ESP_I6x_Gxy2z_a = I_ESP_K6xy_Fx2z_a+ABY*I_ESP_I6x_Fx2z_a;
    Double I_ESP_I5xy_Gxy2z_a = I_ESP_K5x2y_Fx2z_a+ABY*I_ESP_I5xy_Fx2z_a;
    Double I_ESP_I5xz_Gxy2z_a = I_ESP_K5xyz_Fx2z_a+ABY*I_ESP_I5xz_Fx2z_a;
    Double I_ESP_I4x2y_Gxy2z_a = I_ESP_K4x3y_Fx2z_a+ABY*I_ESP_I4x2y_Fx2z_a;
    Double I_ESP_I4xyz_Gxy2z_a = I_ESP_K4x2yz_Fx2z_a+ABY*I_ESP_I4xyz_Fx2z_a;
    Double I_ESP_I4x2z_Gxy2z_a = I_ESP_K4xy2z_Fx2z_a+ABY*I_ESP_I4x2z_Fx2z_a;
    Double I_ESP_I3x3y_Gxy2z_a = I_ESP_K3x4y_Fx2z_a+ABY*I_ESP_I3x3y_Fx2z_a;
    Double I_ESP_I3x2yz_Gxy2z_a = I_ESP_K3x3yz_Fx2z_a+ABY*I_ESP_I3x2yz_Fx2z_a;
    Double I_ESP_I3xy2z_Gxy2z_a = I_ESP_K3x2y2z_Fx2z_a+ABY*I_ESP_I3xy2z_Fx2z_a;
    Double I_ESP_I3x3z_Gxy2z_a = I_ESP_K3xy3z_Fx2z_a+ABY*I_ESP_I3x3z_Fx2z_a;
    Double I_ESP_I2x4y_Gxy2z_a = I_ESP_K2x5y_Fx2z_a+ABY*I_ESP_I2x4y_Fx2z_a;
    Double I_ESP_I2x3yz_Gxy2z_a = I_ESP_K2x4yz_Fx2z_a+ABY*I_ESP_I2x3yz_Fx2z_a;
    Double I_ESP_I2x2y2z_Gxy2z_a = I_ESP_K2x3y2z_Fx2z_a+ABY*I_ESP_I2x2y2z_Fx2z_a;
    Double I_ESP_I2xy3z_Gxy2z_a = I_ESP_K2x2y3z_Fx2z_a+ABY*I_ESP_I2xy3z_Fx2z_a;
    Double I_ESP_I2x4z_Gxy2z_a = I_ESP_K2xy4z_Fx2z_a+ABY*I_ESP_I2x4z_Fx2z_a;
    Double I_ESP_Ix5y_Gxy2z_a = I_ESP_Kx6y_Fx2z_a+ABY*I_ESP_Ix5y_Fx2z_a;
    Double I_ESP_Ix4yz_Gxy2z_a = I_ESP_Kx5yz_Fx2z_a+ABY*I_ESP_Ix4yz_Fx2z_a;
    Double I_ESP_Ix3y2z_Gxy2z_a = I_ESP_Kx4y2z_Fx2z_a+ABY*I_ESP_Ix3y2z_Fx2z_a;
    Double I_ESP_Ix2y3z_Gxy2z_a = I_ESP_Kx3y3z_Fx2z_a+ABY*I_ESP_Ix2y3z_Fx2z_a;
    Double I_ESP_Ixy4z_Gxy2z_a = I_ESP_Kx2y4z_Fx2z_a+ABY*I_ESP_Ixy4z_Fx2z_a;
    Double I_ESP_Ix5z_Gxy2z_a = I_ESP_Kxy5z_Fx2z_a+ABY*I_ESP_Ix5z_Fx2z_a;
    Double I_ESP_I6y_Gxy2z_a = I_ESP_K7y_Fx2z_a+ABY*I_ESP_I6y_Fx2z_a;
    Double I_ESP_I5yz_Gxy2z_a = I_ESP_K6yz_Fx2z_a+ABY*I_ESP_I5yz_Fx2z_a;
    Double I_ESP_I4y2z_Gxy2z_a = I_ESP_K5y2z_Fx2z_a+ABY*I_ESP_I4y2z_Fx2z_a;
    Double I_ESP_I3y3z_Gxy2z_a = I_ESP_K4y3z_Fx2z_a+ABY*I_ESP_I3y3z_Fx2z_a;
    Double I_ESP_I2y4z_Gxy2z_a = I_ESP_K3y4z_Fx2z_a+ABY*I_ESP_I2y4z_Fx2z_a;
    Double I_ESP_Iy5z_Gxy2z_a = I_ESP_K2y5z_Fx2z_a+ABY*I_ESP_Iy5z_Fx2z_a;
    Double I_ESP_I6z_Gxy2z_a = I_ESP_Ky6z_Fx2z_a+ABY*I_ESP_I6z_Fx2z_a;
    Double I_ESP_I6x_Gx3z_a = I_ESP_K7x_F3z_a+ABX*I_ESP_I6x_F3z_a;
    Double I_ESP_I5xy_Gx3z_a = I_ESP_K6xy_F3z_a+ABX*I_ESP_I5xy_F3z_a;
    Double I_ESP_I5xz_Gx3z_a = I_ESP_K6xz_F3z_a+ABX*I_ESP_I5xz_F3z_a;
    Double I_ESP_I4x2y_Gx3z_a = I_ESP_K5x2y_F3z_a+ABX*I_ESP_I4x2y_F3z_a;
    Double I_ESP_I4xyz_Gx3z_a = I_ESP_K5xyz_F3z_a+ABX*I_ESP_I4xyz_F3z_a;
    Double I_ESP_I4x2z_Gx3z_a = I_ESP_K5x2z_F3z_a+ABX*I_ESP_I4x2z_F3z_a;
    Double I_ESP_I3x3y_Gx3z_a = I_ESP_K4x3y_F3z_a+ABX*I_ESP_I3x3y_F3z_a;
    Double I_ESP_I3x2yz_Gx3z_a = I_ESP_K4x2yz_F3z_a+ABX*I_ESP_I3x2yz_F3z_a;
    Double I_ESP_I3xy2z_Gx3z_a = I_ESP_K4xy2z_F3z_a+ABX*I_ESP_I3xy2z_F3z_a;
    Double I_ESP_I3x3z_Gx3z_a = I_ESP_K4x3z_F3z_a+ABX*I_ESP_I3x3z_F3z_a;
    Double I_ESP_I2x4y_Gx3z_a = I_ESP_K3x4y_F3z_a+ABX*I_ESP_I2x4y_F3z_a;
    Double I_ESP_I2x3yz_Gx3z_a = I_ESP_K3x3yz_F3z_a+ABX*I_ESP_I2x3yz_F3z_a;
    Double I_ESP_I2x2y2z_Gx3z_a = I_ESP_K3x2y2z_F3z_a+ABX*I_ESP_I2x2y2z_F3z_a;
    Double I_ESP_I2xy3z_Gx3z_a = I_ESP_K3xy3z_F3z_a+ABX*I_ESP_I2xy3z_F3z_a;
    Double I_ESP_I2x4z_Gx3z_a = I_ESP_K3x4z_F3z_a+ABX*I_ESP_I2x4z_F3z_a;
    Double I_ESP_Ix5y_Gx3z_a = I_ESP_K2x5y_F3z_a+ABX*I_ESP_Ix5y_F3z_a;
    Double I_ESP_Ix4yz_Gx3z_a = I_ESP_K2x4yz_F3z_a+ABX*I_ESP_Ix4yz_F3z_a;
    Double I_ESP_Ix3y2z_Gx3z_a = I_ESP_K2x3y2z_F3z_a+ABX*I_ESP_Ix3y2z_F3z_a;
    Double I_ESP_Ix2y3z_Gx3z_a = I_ESP_K2x2y3z_F3z_a+ABX*I_ESP_Ix2y3z_F3z_a;
    Double I_ESP_Ixy4z_Gx3z_a = I_ESP_K2xy4z_F3z_a+ABX*I_ESP_Ixy4z_F3z_a;
    Double I_ESP_Ix5z_Gx3z_a = I_ESP_K2x5z_F3z_a+ABX*I_ESP_Ix5z_F3z_a;
    Double I_ESP_I6y_Gx3z_a = I_ESP_Kx6y_F3z_a+ABX*I_ESP_I6y_F3z_a;
    Double I_ESP_I5yz_Gx3z_a = I_ESP_Kx5yz_F3z_a+ABX*I_ESP_I5yz_F3z_a;
    Double I_ESP_I4y2z_Gx3z_a = I_ESP_Kx4y2z_F3z_a+ABX*I_ESP_I4y2z_F3z_a;
    Double I_ESP_I3y3z_Gx3z_a = I_ESP_Kx3y3z_F3z_a+ABX*I_ESP_I3y3z_F3z_a;
    Double I_ESP_I2y4z_Gx3z_a = I_ESP_Kx2y4z_F3z_a+ABX*I_ESP_I2y4z_F3z_a;
    Double I_ESP_Iy5z_Gx3z_a = I_ESP_Kxy5z_F3z_a+ABX*I_ESP_Iy5z_F3z_a;
    Double I_ESP_I6z_Gx3z_a = I_ESP_Kx6z_F3z_a+ABX*I_ESP_I6z_F3z_a;
    Double I_ESP_I6x_G4y_a = I_ESP_K6xy_F3y_a+ABY*I_ESP_I6x_F3y_a;
    Double I_ESP_I5xy_G4y_a = I_ESP_K5x2y_F3y_a+ABY*I_ESP_I5xy_F3y_a;
    Double I_ESP_I5xz_G4y_a = I_ESP_K5xyz_F3y_a+ABY*I_ESP_I5xz_F3y_a;
    Double I_ESP_I4x2y_G4y_a = I_ESP_K4x3y_F3y_a+ABY*I_ESP_I4x2y_F3y_a;
    Double I_ESP_I4xyz_G4y_a = I_ESP_K4x2yz_F3y_a+ABY*I_ESP_I4xyz_F3y_a;
    Double I_ESP_I4x2z_G4y_a = I_ESP_K4xy2z_F3y_a+ABY*I_ESP_I4x2z_F3y_a;
    Double I_ESP_I3x3y_G4y_a = I_ESP_K3x4y_F3y_a+ABY*I_ESP_I3x3y_F3y_a;
    Double I_ESP_I3x2yz_G4y_a = I_ESP_K3x3yz_F3y_a+ABY*I_ESP_I3x2yz_F3y_a;
    Double I_ESP_I3xy2z_G4y_a = I_ESP_K3x2y2z_F3y_a+ABY*I_ESP_I3xy2z_F3y_a;
    Double I_ESP_I3x3z_G4y_a = I_ESP_K3xy3z_F3y_a+ABY*I_ESP_I3x3z_F3y_a;
    Double I_ESP_I2x4y_G4y_a = I_ESP_K2x5y_F3y_a+ABY*I_ESP_I2x4y_F3y_a;
    Double I_ESP_I2x3yz_G4y_a = I_ESP_K2x4yz_F3y_a+ABY*I_ESP_I2x3yz_F3y_a;
    Double I_ESP_I2x2y2z_G4y_a = I_ESP_K2x3y2z_F3y_a+ABY*I_ESP_I2x2y2z_F3y_a;
    Double I_ESP_I2xy3z_G4y_a = I_ESP_K2x2y3z_F3y_a+ABY*I_ESP_I2xy3z_F3y_a;
    Double I_ESP_I2x4z_G4y_a = I_ESP_K2xy4z_F3y_a+ABY*I_ESP_I2x4z_F3y_a;
    Double I_ESP_Ix5y_G4y_a = I_ESP_Kx6y_F3y_a+ABY*I_ESP_Ix5y_F3y_a;
    Double I_ESP_Ix4yz_G4y_a = I_ESP_Kx5yz_F3y_a+ABY*I_ESP_Ix4yz_F3y_a;
    Double I_ESP_Ix3y2z_G4y_a = I_ESP_Kx4y2z_F3y_a+ABY*I_ESP_Ix3y2z_F3y_a;
    Double I_ESP_Ix2y3z_G4y_a = I_ESP_Kx3y3z_F3y_a+ABY*I_ESP_Ix2y3z_F3y_a;
    Double I_ESP_Ixy4z_G4y_a = I_ESP_Kx2y4z_F3y_a+ABY*I_ESP_Ixy4z_F3y_a;
    Double I_ESP_Ix5z_G4y_a = I_ESP_Kxy5z_F3y_a+ABY*I_ESP_Ix5z_F3y_a;
    Double I_ESP_I6y_G4y_a = I_ESP_K7y_F3y_a+ABY*I_ESP_I6y_F3y_a;
    Double I_ESP_I5yz_G4y_a = I_ESP_K6yz_F3y_a+ABY*I_ESP_I5yz_F3y_a;
    Double I_ESP_I4y2z_G4y_a = I_ESP_K5y2z_F3y_a+ABY*I_ESP_I4y2z_F3y_a;
    Double I_ESP_I3y3z_G4y_a = I_ESP_K4y3z_F3y_a+ABY*I_ESP_I3y3z_F3y_a;
    Double I_ESP_I2y4z_G4y_a = I_ESP_K3y4z_F3y_a+ABY*I_ESP_I2y4z_F3y_a;
    Double I_ESP_Iy5z_G4y_a = I_ESP_K2y5z_F3y_a+ABY*I_ESP_Iy5z_F3y_a;
    Double I_ESP_I6z_G4y_a = I_ESP_Ky6z_F3y_a+ABY*I_ESP_I6z_F3y_a;
    Double I_ESP_I6x_G3yz_a = I_ESP_K6xz_F3y_a+ABZ*I_ESP_I6x_F3y_a;
    Double I_ESP_I5xy_G3yz_a = I_ESP_K5xyz_F3y_a+ABZ*I_ESP_I5xy_F3y_a;
    Double I_ESP_I5xz_G3yz_a = I_ESP_K5x2z_F3y_a+ABZ*I_ESP_I5xz_F3y_a;
    Double I_ESP_I4x2y_G3yz_a = I_ESP_K4x2yz_F3y_a+ABZ*I_ESP_I4x2y_F3y_a;
    Double I_ESP_I4xyz_G3yz_a = I_ESP_K4xy2z_F3y_a+ABZ*I_ESP_I4xyz_F3y_a;
    Double I_ESP_I4x2z_G3yz_a = I_ESP_K4x3z_F3y_a+ABZ*I_ESP_I4x2z_F3y_a;
    Double I_ESP_I3x3y_G3yz_a = I_ESP_K3x3yz_F3y_a+ABZ*I_ESP_I3x3y_F3y_a;
    Double I_ESP_I3x2yz_G3yz_a = I_ESP_K3x2y2z_F3y_a+ABZ*I_ESP_I3x2yz_F3y_a;
    Double I_ESP_I3xy2z_G3yz_a = I_ESP_K3xy3z_F3y_a+ABZ*I_ESP_I3xy2z_F3y_a;
    Double I_ESP_I3x3z_G3yz_a = I_ESP_K3x4z_F3y_a+ABZ*I_ESP_I3x3z_F3y_a;
    Double I_ESP_I2x4y_G3yz_a = I_ESP_K2x4yz_F3y_a+ABZ*I_ESP_I2x4y_F3y_a;
    Double I_ESP_I2x3yz_G3yz_a = I_ESP_K2x3y2z_F3y_a+ABZ*I_ESP_I2x3yz_F3y_a;
    Double I_ESP_I2x2y2z_G3yz_a = I_ESP_K2x2y3z_F3y_a+ABZ*I_ESP_I2x2y2z_F3y_a;
    Double I_ESP_I2xy3z_G3yz_a = I_ESP_K2xy4z_F3y_a+ABZ*I_ESP_I2xy3z_F3y_a;
    Double I_ESP_I2x4z_G3yz_a = I_ESP_K2x5z_F3y_a+ABZ*I_ESP_I2x4z_F3y_a;
    Double I_ESP_Ix5y_G3yz_a = I_ESP_Kx5yz_F3y_a+ABZ*I_ESP_Ix5y_F3y_a;
    Double I_ESP_Ix4yz_G3yz_a = I_ESP_Kx4y2z_F3y_a+ABZ*I_ESP_Ix4yz_F3y_a;
    Double I_ESP_Ix3y2z_G3yz_a = I_ESP_Kx3y3z_F3y_a+ABZ*I_ESP_Ix3y2z_F3y_a;
    Double I_ESP_Ix2y3z_G3yz_a = I_ESP_Kx2y4z_F3y_a+ABZ*I_ESP_Ix2y3z_F3y_a;
    Double I_ESP_Ixy4z_G3yz_a = I_ESP_Kxy5z_F3y_a+ABZ*I_ESP_Ixy4z_F3y_a;
    Double I_ESP_Ix5z_G3yz_a = I_ESP_Kx6z_F3y_a+ABZ*I_ESP_Ix5z_F3y_a;
    Double I_ESP_I6y_G3yz_a = I_ESP_K6yz_F3y_a+ABZ*I_ESP_I6y_F3y_a;
    Double I_ESP_I5yz_G3yz_a = I_ESP_K5y2z_F3y_a+ABZ*I_ESP_I5yz_F3y_a;
    Double I_ESP_I4y2z_G3yz_a = I_ESP_K4y3z_F3y_a+ABZ*I_ESP_I4y2z_F3y_a;
    Double I_ESP_I3y3z_G3yz_a = I_ESP_K3y4z_F3y_a+ABZ*I_ESP_I3y3z_F3y_a;
    Double I_ESP_I2y4z_G3yz_a = I_ESP_K2y5z_F3y_a+ABZ*I_ESP_I2y4z_F3y_a;
    Double I_ESP_Iy5z_G3yz_a = I_ESP_Ky6z_F3y_a+ABZ*I_ESP_Iy5z_F3y_a;
    Double I_ESP_I6z_G3yz_a = I_ESP_K7z_F3y_a+ABZ*I_ESP_I6z_F3y_a;
    Double I_ESP_I6x_G2y2z_a = I_ESP_K6xz_F2yz_a+ABZ*I_ESP_I6x_F2yz_a;
    Double I_ESP_I5xy_G2y2z_a = I_ESP_K5xyz_F2yz_a+ABZ*I_ESP_I5xy_F2yz_a;
    Double I_ESP_I5xz_G2y2z_a = I_ESP_K5x2z_F2yz_a+ABZ*I_ESP_I5xz_F2yz_a;
    Double I_ESP_I4x2y_G2y2z_a = I_ESP_K4x2yz_F2yz_a+ABZ*I_ESP_I4x2y_F2yz_a;
    Double I_ESP_I4xyz_G2y2z_a = I_ESP_K4xy2z_F2yz_a+ABZ*I_ESP_I4xyz_F2yz_a;
    Double I_ESP_I4x2z_G2y2z_a = I_ESP_K4x3z_F2yz_a+ABZ*I_ESP_I4x2z_F2yz_a;
    Double I_ESP_I3x3y_G2y2z_a = I_ESP_K3x3yz_F2yz_a+ABZ*I_ESP_I3x3y_F2yz_a;
    Double I_ESP_I3x2yz_G2y2z_a = I_ESP_K3x2y2z_F2yz_a+ABZ*I_ESP_I3x2yz_F2yz_a;
    Double I_ESP_I3xy2z_G2y2z_a = I_ESP_K3xy3z_F2yz_a+ABZ*I_ESP_I3xy2z_F2yz_a;
    Double I_ESP_I3x3z_G2y2z_a = I_ESP_K3x4z_F2yz_a+ABZ*I_ESP_I3x3z_F2yz_a;
    Double I_ESP_I2x4y_G2y2z_a = I_ESP_K2x4yz_F2yz_a+ABZ*I_ESP_I2x4y_F2yz_a;
    Double I_ESP_I2x3yz_G2y2z_a = I_ESP_K2x3y2z_F2yz_a+ABZ*I_ESP_I2x3yz_F2yz_a;
    Double I_ESP_I2x2y2z_G2y2z_a = I_ESP_K2x2y3z_F2yz_a+ABZ*I_ESP_I2x2y2z_F2yz_a;
    Double I_ESP_I2xy3z_G2y2z_a = I_ESP_K2xy4z_F2yz_a+ABZ*I_ESP_I2xy3z_F2yz_a;
    Double I_ESP_I2x4z_G2y2z_a = I_ESP_K2x5z_F2yz_a+ABZ*I_ESP_I2x4z_F2yz_a;
    Double I_ESP_Ix5y_G2y2z_a = I_ESP_Kx5yz_F2yz_a+ABZ*I_ESP_Ix5y_F2yz_a;
    Double I_ESP_Ix4yz_G2y2z_a = I_ESP_Kx4y2z_F2yz_a+ABZ*I_ESP_Ix4yz_F2yz_a;
    Double I_ESP_Ix3y2z_G2y2z_a = I_ESP_Kx3y3z_F2yz_a+ABZ*I_ESP_Ix3y2z_F2yz_a;
    Double I_ESP_Ix2y3z_G2y2z_a = I_ESP_Kx2y4z_F2yz_a+ABZ*I_ESP_Ix2y3z_F2yz_a;
    Double I_ESP_Ixy4z_G2y2z_a = I_ESP_Kxy5z_F2yz_a+ABZ*I_ESP_Ixy4z_F2yz_a;
    Double I_ESP_Ix5z_G2y2z_a = I_ESP_Kx6z_F2yz_a+ABZ*I_ESP_Ix5z_F2yz_a;
    Double I_ESP_I6y_G2y2z_a = I_ESP_K6yz_F2yz_a+ABZ*I_ESP_I6y_F2yz_a;
    Double I_ESP_I5yz_G2y2z_a = I_ESP_K5y2z_F2yz_a+ABZ*I_ESP_I5yz_F2yz_a;
    Double I_ESP_I4y2z_G2y2z_a = I_ESP_K4y3z_F2yz_a+ABZ*I_ESP_I4y2z_F2yz_a;
    Double I_ESP_I3y3z_G2y2z_a = I_ESP_K3y4z_F2yz_a+ABZ*I_ESP_I3y3z_F2yz_a;
    Double I_ESP_I2y4z_G2y2z_a = I_ESP_K2y5z_F2yz_a+ABZ*I_ESP_I2y4z_F2yz_a;
    Double I_ESP_Iy5z_G2y2z_a = I_ESP_Ky6z_F2yz_a+ABZ*I_ESP_Iy5z_F2yz_a;
    Double I_ESP_I6z_G2y2z_a = I_ESP_K7z_F2yz_a+ABZ*I_ESP_I6z_F2yz_a;
    Double I_ESP_I6x_Gy3z_a = I_ESP_K6xy_F3z_a+ABY*I_ESP_I6x_F3z_a;
    Double I_ESP_I5xy_Gy3z_a = I_ESP_K5x2y_F3z_a+ABY*I_ESP_I5xy_F3z_a;
    Double I_ESP_I5xz_Gy3z_a = I_ESP_K5xyz_F3z_a+ABY*I_ESP_I5xz_F3z_a;
    Double I_ESP_I4x2y_Gy3z_a = I_ESP_K4x3y_F3z_a+ABY*I_ESP_I4x2y_F3z_a;
    Double I_ESP_I4xyz_Gy3z_a = I_ESP_K4x2yz_F3z_a+ABY*I_ESP_I4xyz_F3z_a;
    Double I_ESP_I4x2z_Gy3z_a = I_ESP_K4xy2z_F3z_a+ABY*I_ESP_I4x2z_F3z_a;
    Double I_ESP_I3x3y_Gy3z_a = I_ESP_K3x4y_F3z_a+ABY*I_ESP_I3x3y_F3z_a;
    Double I_ESP_I3x2yz_Gy3z_a = I_ESP_K3x3yz_F3z_a+ABY*I_ESP_I3x2yz_F3z_a;
    Double I_ESP_I3xy2z_Gy3z_a = I_ESP_K3x2y2z_F3z_a+ABY*I_ESP_I3xy2z_F3z_a;
    Double I_ESP_I3x3z_Gy3z_a = I_ESP_K3xy3z_F3z_a+ABY*I_ESP_I3x3z_F3z_a;
    Double I_ESP_I2x4y_Gy3z_a = I_ESP_K2x5y_F3z_a+ABY*I_ESP_I2x4y_F3z_a;
    Double I_ESP_I2x3yz_Gy3z_a = I_ESP_K2x4yz_F3z_a+ABY*I_ESP_I2x3yz_F3z_a;
    Double I_ESP_I2x2y2z_Gy3z_a = I_ESP_K2x3y2z_F3z_a+ABY*I_ESP_I2x2y2z_F3z_a;
    Double I_ESP_I2xy3z_Gy3z_a = I_ESP_K2x2y3z_F3z_a+ABY*I_ESP_I2xy3z_F3z_a;
    Double I_ESP_I2x4z_Gy3z_a = I_ESP_K2xy4z_F3z_a+ABY*I_ESP_I2x4z_F3z_a;
    Double I_ESP_Ix5y_Gy3z_a = I_ESP_Kx6y_F3z_a+ABY*I_ESP_Ix5y_F3z_a;
    Double I_ESP_Ix4yz_Gy3z_a = I_ESP_Kx5yz_F3z_a+ABY*I_ESP_Ix4yz_F3z_a;
    Double I_ESP_Ix3y2z_Gy3z_a = I_ESP_Kx4y2z_F3z_a+ABY*I_ESP_Ix3y2z_F3z_a;
    Double I_ESP_Ix2y3z_Gy3z_a = I_ESP_Kx3y3z_F3z_a+ABY*I_ESP_Ix2y3z_F3z_a;
    Double I_ESP_Ixy4z_Gy3z_a = I_ESP_Kx2y4z_F3z_a+ABY*I_ESP_Ixy4z_F3z_a;
    Double I_ESP_Ix5z_Gy3z_a = I_ESP_Kxy5z_F3z_a+ABY*I_ESP_Ix5z_F3z_a;
    Double I_ESP_I6y_Gy3z_a = I_ESP_K7y_F3z_a+ABY*I_ESP_I6y_F3z_a;
    Double I_ESP_I5yz_Gy3z_a = I_ESP_K6yz_F3z_a+ABY*I_ESP_I5yz_F3z_a;
    Double I_ESP_I4y2z_Gy3z_a = I_ESP_K5y2z_F3z_a+ABY*I_ESP_I4y2z_F3z_a;
    Double I_ESP_I3y3z_Gy3z_a = I_ESP_K4y3z_F3z_a+ABY*I_ESP_I3y3z_F3z_a;
    Double I_ESP_I2y4z_Gy3z_a = I_ESP_K3y4z_F3z_a+ABY*I_ESP_I2y4z_F3z_a;
    Double I_ESP_Iy5z_Gy3z_a = I_ESP_K2y5z_F3z_a+ABY*I_ESP_Iy5z_F3z_a;
    Double I_ESP_I6z_Gy3z_a = I_ESP_Ky6z_F3z_a+ABY*I_ESP_I6z_F3z_a;
    Double I_ESP_I6x_G4z_a = I_ESP_K6xz_F3z_a+ABZ*I_ESP_I6x_F3z_a;
    Double I_ESP_I5xy_G4z_a = I_ESP_K5xyz_F3z_a+ABZ*I_ESP_I5xy_F3z_a;
    Double I_ESP_I5xz_G4z_a = I_ESP_K5x2z_F3z_a+ABZ*I_ESP_I5xz_F3z_a;
    Double I_ESP_I4x2y_G4z_a = I_ESP_K4x2yz_F3z_a+ABZ*I_ESP_I4x2y_F3z_a;
    Double I_ESP_I4xyz_G4z_a = I_ESP_K4xy2z_F3z_a+ABZ*I_ESP_I4xyz_F3z_a;
    Double I_ESP_I4x2z_G4z_a = I_ESP_K4x3z_F3z_a+ABZ*I_ESP_I4x2z_F3z_a;
    Double I_ESP_I3x3y_G4z_a = I_ESP_K3x3yz_F3z_a+ABZ*I_ESP_I3x3y_F3z_a;
    Double I_ESP_I3x2yz_G4z_a = I_ESP_K3x2y2z_F3z_a+ABZ*I_ESP_I3x2yz_F3z_a;
    Double I_ESP_I3xy2z_G4z_a = I_ESP_K3xy3z_F3z_a+ABZ*I_ESP_I3xy2z_F3z_a;
    Double I_ESP_I3x3z_G4z_a = I_ESP_K3x4z_F3z_a+ABZ*I_ESP_I3x3z_F3z_a;
    Double I_ESP_I2x4y_G4z_a = I_ESP_K2x4yz_F3z_a+ABZ*I_ESP_I2x4y_F3z_a;
    Double I_ESP_I2x3yz_G4z_a = I_ESP_K2x3y2z_F3z_a+ABZ*I_ESP_I2x3yz_F3z_a;
    Double I_ESP_I2x2y2z_G4z_a = I_ESP_K2x2y3z_F3z_a+ABZ*I_ESP_I2x2y2z_F3z_a;
    Double I_ESP_I2xy3z_G4z_a = I_ESP_K2xy4z_F3z_a+ABZ*I_ESP_I2xy3z_F3z_a;
    Double I_ESP_I2x4z_G4z_a = I_ESP_K2x5z_F3z_a+ABZ*I_ESP_I2x4z_F3z_a;
    Double I_ESP_Ix5y_G4z_a = I_ESP_Kx5yz_F3z_a+ABZ*I_ESP_Ix5y_F3z_a;
    Double I_ESP_Ix4yz_G4z_a = I_ESP_Kx4y2z_F3z_a+ABZ*I_ESP_Ix4yz_F3z_a;
    Double I_ESP_Ix3y2z_G4z_a = I_ESP_Kx3y3z_F3z_a+ABZ*I_ESP_Ix3y2z_F3z_a;
    Double I_ESP_Ix2y3z_G4z_a = I_ESP_Kx2y4z_F3z_a+ABZ*I_ESP_Ix2y3z_F3z_a;
    Double I_ESP_Ixy4z_G4z_a = I_ESP_Kxy5z_F3z_a+ABZ*I_ESP_Ixy4z_F3z_a;
    Double I_ESP_Ix5z_G4z_a = I_ESP_Kx6z_F3z_a+ABZ*I_ESP_Ix5z_F3z_a;
    Double I_ESP_I6y_G4z_a = I_ESP_K6yz_F3z_a+ABZ*I_ESP_I6y_F3z_a;
    Double I_ESP_I5yz_G4z_a = I_ESP_K5y2z_F3z_a+ABZ*I_ESP_I5yz_F3z_a;
    Double I_ESP_I4y2z_G4z_a = I_ESP_K4y3z_F3z_a+ABZ*I_ESP_I4y2z_F3z_a;
    Double I_ESP_I3y3z_G4z_a = I_ESP_K3y4z_F3z_a+ABZ*I_ESP_I3y3z_F3z_a;
    Double I_ESP_I2y4z_G4z_a = I_ESP_K2y5z_F3z_a+ABZ*I_ESP_I2y4z_F3z_a;
    Double I_ESP_Iy5z_G4z_a = I_ESP_Ky6z_F3z_a+ABZ*I_ESP_Iy5z_F3z_a;
    Double I_ESP_I6z_G4z_a = I_ESP_K7z_F3z_a+ABZ*I_ESP_I6z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G_a
     * RHS shell quartet name: SQ_ESP_G_G
     ************************************************************/
    abcd[iGrid*945+0] = 2.0E0*I_ESP_I6x_G4x_a-5*I_ESP_G4x_G4x;
    abcd[iGrid*945+1] = 2.0E0*I_ESP_I5xy_G4x_a-4*I_ESP_G3xy_G4x;
    abcd[iGrid*945+2] = 2.0E0*I_ESP_I5xz_G4x_a-4*I_ESP_G3xz_G4x;
    abcd[iGrid*945+3] = 2.0E0*I_ESP_I4x2y_G4x_a-3*I_ESP_G2x2y_G4x;
    abcd[iGrid*945+4] = 2.0E0*I_ESP_I4xyz_G4x_a-3*I_ESP_G2xyz_G4x;
    abcd[iGrid*945+5] = 2.0E0*I_ESP_I4x2z_G4x_a-3*I_ESP_G2x2z_G4x;
    abcd[iGrid*945+6] = 2.0E0*I_ESP_I3x3y_G4x_a-2*I_ESP_Gx3y_G4x;
    abcd[iGrid*945+7] = 2.0E0*I_ESP_I3x2yz_G4x_a-2*I_ESP_Gx2yz_G4x;
    abcd[iGrid*945+8] = 2.0E0*I_ESP_I3xy2z_G4x_a-2*I_ESP_Gxy2z_G4x;
    abcd[iGrid*945+9] = 2.0E0*I_ESP_I3x3z_G4x_a-2*I_ESP_Gx3z_G4x;
    abcd[iGrid*945+10] = 2.0E0*I_ESP_I2x4y_G4x_a-1*I_ESP_G4y_G4x;
    abcd[iGrid*945+11] = 2.0E0*I_ESP_I2x3yz_G4x_a-1*I_ESP_G3yz_G4x;
    abcd[iGrid*945+12] = 2.0E0*I_ESP_I2x2y2z_G4x_a-1*I_ESP_G2y2z_G4x;
    abcd[iGrid*945+13] = 2.0E0*I_ESP_I2xy3z_G4x_a-1*I_ESP_Gy3z_G4x;
    abcd[iGrid*945+14] = 2.0E0*I_ESP_I2x4z_G4x_a-1*I_ESP_G4z_G4x;
    abcd[iGrid*945+15] = 2.0E0*I_ESP_Ix5y_G4x_a;
    abcd[iGrid*945+16] = 2.0E0*I_ESP_Ix4yz_G4x_a;
    abcd[iGrid*945+17] = 2.0E0*I_ESP_Ix3y2z_G4x_a;
    abcd[iGrid*945+18] = 2.0E0*I_ESP_Ix2y3z_G4x_a;
    abcd[iGrid*945+19] = 2.0E0*I_ESP_Ixy4z_G4x_a;
    abcd[iGrid*945+20] = 2.0E0*I_ESP_Ix5z_G4x_a;
    abcd[iGrid*945+21] = 2.0E0*I_ESP_I6x_G3xy_a-5*I_ESP_G4x_G3xy;
    abcd[iGrid*945+22] = 2.0E0*I_ESP_I5xy_G3xy_a-4*I_ESP_G3xy_G3xy;
    abcd[iGrid*945+23] = 2.0E0*I_ESP_I5xz_G3xy_a-4*I_ESP_G3xz_G3xy;
    abcd[iGrid*945+24] = 2.0E0*I_ESP_I4x2y_G3xy_a-3*I_ESP_G2x2y_G3xy;
    abcd[iGrid*945+25] = 2.0E0*I_ESP_I4xyz_G3xy_a-3*I_ESP_G2xyz_G3xy;
    abcd[iGrid*945+26] = 2.0E0*I_ESP_I4x2z_G3xy_a-3*I_ESP_G2x2z_G3xy;
    abcd[iGrid*945+27] = 2.0E0*I_ESP_I3x3y_G3xy_a-2*I_ESP_Gx3y_G3xy;
    abcd[iGrid*945+28] = 2.0E0*I_ESP_I3x2yz_G3xy_a-2*I_ESP_Gx2yz_G3xy;
    abcd[iGrid*945+29] = 2.0E0*I_ESP_I3xy2z_G3xy_a-2*I_ESP_Gxy2z_G3xy;
    abcd[iGrid*945+30] = 2.0E0*I_ESP_I3x3z_G3xy_a-2*I_ESP_Gx3z_G3xy;
    abcd[iGrid*945+31] = 2.0E0*I_ESP_I2x4y_G3xy_a-1*I_ESP_G4y_G3xy;
    abcd[iGrid*945+32] = 2.0E0*I_ESP_I2x3yz_G3xy_a-1*I_ESP_G3yz_G3xy;
    abcd[iGrid*945+33] = 2.0E0*I_ESP_I2x2y2z_G3xy_a-1*I_ESP_G2y2z_G3xy;
    abcd[iGrid*945+34] = 2.0E0*I_ESP_I2xy3z_G3xy_a-1*I_ESP_Gy3z_G3xy;
    abcd[iGrid*945+35] = 2.0E0*I_ESP_I2x4z_G3xy_a-1*I_ESP_G4z_G3xy;
    abcd[iGrid*945+36] = 2.0E0*I_ESP_Ix5y_G3xy_a;
    abcd[iGrid*945+37] = 2.0E0*I_ESP_Ix4yz_G3xy_a;
    abcd[iGrid*945+38] = 2.0E0*I_ESP_Ix3y2z_G3xy_a;
    abcd[iGrid*945+39] = 2.0E0*I_ESP_Ix2y3z_G3xy_a;
    abcd[iGrid*945+40] = 2.0E0*I_ESP_Ixy4z_G3xy_a;
    abcd[iGrid*945+41] = 2.0E0*I_ESP_Ix5z_G3xy_a;
    abcd[iGrid*945+42] = 2.0E0*I_ESP_I6x_G3xz_a-5*I_ESP_G4x_G3xz;
    abcd[iGrid*945+43] = 2.0E0*I_ESP_I5xy_G3xz_a-4*I_ESP_G3xy_G3xz;
    abcd[iGrid*945+44] = 2.0E0*I_ESP_I5xz_G3xz_a-4*I_ESP_G3xz_G3xz;
    abcd[iGrid*945+45] = 2.0E0*I_ESP_I4x2y_G3xz_a-3*I_ESP_G2x2y_G3xz;
    abcd[iGrid*945+46] = 2.0E0*I_ESP_I4xyz_G3xz_a-3*I_ESP_G2xyz_G3xz;
    abcd[iGrid*945+47] = 2.0E0*I_ESP_I4x2z_G3xz_a-3*I_ESP_G2x2z_G3xz;
    abcd[iGrid*945+48] = 2.0E0*I_ESP_I3x3y_G3xz_a-2*I_ESP_Gx3y_G3xz;
    abcd[iGrid*945+49] = 2.0E0*I_ESP_I3x2yz_G3xz_a-2*I_ESP_Gx2yz_G3xz;
    abcd[iGrid*945+50] = 2.0E0*I_ESP_I3xy2z_G3xz_a-2*I_ESP_Gxy2z_G3xz;
    abcd[iGrid*945+51] = 2.0E0*I_ESP_I3x3z_G3xz_a-2*I_ESP_Gx3z_G3xz;
    abcd[iGrid*945+52] = 2.0E0*I_ESP_I2x4y_G3xz_a-1*I_ESP_G4y_G3xz;
    abcd[iGrid*945+53] = 2.0E0*I_ESP_I2x3yz_G3xz_a-1*I_ESP_G3yz_G3xz;
    abcd[iGrid*945+54] = 2.0E0*I_ESP_I2x2y2z_G3xz_a-1*I_ESP_G2y2z_G3xz;
    abcd[iGrid*945+55] = 2.0E0*I_ESP_I2xy3z_G3xz_a-1*I_ESP_Gy3z_G3xz;
    abcd[iGrid*945+56] = 2.0E0*I_ESP_I2x4z_G3xz_a-1*I_ESP_G4z_G3xz;
    abcd[iGrid*945+57] = 2.0E0*I_ESP_Ix5y_G3xz_a;
    abcd[iGrid*945+58] = 2.0E0*I_ESP_Ix4yz_G3xz_a;
    abcd[iGrid*945+59] = 2.0E0*I_ESP_Ix3y2z_G3xz_a;
    abcd[iGrid*945+60] = 2.0E0*I_ESP_Ix2y3z_G3xz_a;
    abcd[iGrid*945+61] = 2.0E0*I_ESP_Ixy4z_G3xz_a;
    abcd[iGrid*945+62] = 2.0E0*I_ESP_Ix5z_G3xz_a;
    abcd[iGrid*945+63] = 2.0E0*I_ESP_I6x_G2x2y_a-5*I_ESP_G4x_G2x2y;
    abcd[iGrid*945+64] = 2.0E0*I_ESP_I5xy_G2x2y_a-4*I_ESP_G3xy_G2x2y;
    abcd[iGrid*945+65] = 2.0E0*I_ESP_I5xz_G2x2y_a-4*I_ESP_G3xz_G2x2y;
    abcd[iGrid*945+66] = 2.0E0*I_ESP_I4x2y_G2x2y_a-3*I_ESP_G2x2y_G2x2y;
    abcd[iGrid*945+67] = 2.0E0*I_ESP_I4xyz_G2x2y_a-3*I_ESP_G2xyz_G2x2y;
    abcd[iGrid*945+68] = 2.0E0*I_ESP_I4x2z_G2x2y_a-3*I_ESP_G2x2z_G2x2y;
    abcd[iGrid*945+69] = 2.0E0*I_ESP_I3x3y_G2x2y_a-2*I_ESP_Gx3y_G2x2y;
    abcd[iGrid*945+70] = 2.0E0*I_ESP_I3x2yz_G2x2y_a-2*I_ESP_Gx2yz_G2x2y;
    abcd[iGrid*945+71] = 2.0E0*I_ESP_I3xy2z_G2x2y_a-2*I_ESP_Gxy2z_G2x2y;
    abcd[iGrid*945+72] = 2.0E0*I_ESP_I3x3z_G2x2y_a-2*I_ESP_Gx3z_G2x2y;
    abcd[iGrid*945+73] = 2.0E0*I_ESP_I2x4y_G2x2y_a-1*I_ESP_G4y_G2x2y;
    abcd[iGrid*945+74] = 2.0E0*I_ESP_I2x3yz_G2x2y_a-1*I_ESP_G3yz_G2x2y;
    abcd[iGrid*945+75] = 2.0E0*I_ESP_I2x2y2z_G2x2y_a-1*I_ESP_G2y2z_G2x2y;
    abcd[iGrid*945+76] = 2.0E0*I_ESP_I2xy3z_G2x2y_a-1*I_ESP_Gy3z_G2x2y;
    abcd[iGrid*945+77] = 2.0E0*I_ESP_I2x4z_G2x2y_a-1*I_ESP_G4z_G2x2y;
    abcd[iGrid*945+78] = 2.0E0*I_ESP_Ix5y_G2x2y_a;
    abcd[iGrid*945+79] = 2.0E0*I_ESP_Ix4yz_G2x2y_a;
    abcd[iGrid*945+80] = 2.0E0*I_ESP_Ix3y2z_G2x2y_a;
    abcd[iGrid*945+81] = 2.0E0*I_ESP_Ix2y3z_G2x2y_a;
    abcd[iGrid*945+82] = 2.0E0*I_ESP_Ixy4z_G2x2y_a;
    abcd[iGrid*945+83] = 2.0E0*I_ESP_Ix5z_G2x2y_a;
    abcd[iGrid*945+84] = 2.0E0*I_ESP_I6x_G2xyz_a-5*I_ESP_G4x_G2xyz;
    abcd[iGrid*945+85] = 2.0E0*I_ESP_I5xy_G2xyz_a-4*I_ESP_G3xy_G2xyz;
    abcd[iGrid*945+86] = 2.0E0*I_ESP_I5xz_G2xyz_a-4*I_ESP_G3xz_G2xyz;
    abcd[iGrid*945+87] = 2.0E0*I_ESP_I4x2y_G2xyz_a-3*I_ESP_G2x2y_G2xyz;
    abcd[iGrid*945+88] = 2.0E0*I_ESP_I4xyz_G2xyz_a-3*I_ESP_G2xyz_G2xyz;
    abcd[iGrid*945+89] = 2.0E0*I_ESP_I4x2z_G2xyz_a-3*I_ESP_G2x2z_G2xyz;
    abcd[iGrid*945+90] = 2.0E0*I_ESP_I3x3y_G2xyz_a-2*I_ESP_Gx3y_G2xyz;
    abcd[iGrid*945+91] = 2.0E0*I_ESP_I3x2yz_G2xyz_a-2*I_ESP_Gx2yz_G2xyz;
    abcd[iGrid*945+92] = 2.0E0*I_ESP_I3xy2z_G2xyz_a-2*I_ESP_Gxy2z_G2xyz;
    abcd[iGrid*945+93] = 2.0E0*I_ESP_I3x3z_G2xyz_a-2*I_ESP_Gx3z_G2xyz;
    abcd[iGrid*945+94] = 2.0E0*I_ESP_I2x4y_G2xyz_a-1*I_ESP_G4y_G2xyz;
    abcd[iGrid*945+95] = 2.0E0*I_ESP_I2x3yz_G2xyz_a-1*I_ESP_G3yz_G2xyz;
    abcd[iGrid*945+96] = 2.0E0*I_ESP_I2x2y2z_G2xyz_a-1*I_ESP_G2y2z_G2xyz;
    abcd[iGrid*945+97] = 2.0E0*I_ESP_I2xy3z_G2xyz_a-1*I_ESP_Gy3z_G2xyz;
    abcd[iGrid*945+98] = 2.0E0*I_ESP_I2x4z_G2xyz_a-1*I_ESP_G4z_G2xyz;
    abcd[iGrid*945+99] = 2.0E0*I_ESP_Ix5y_G2xyz_a;
    abcd[iGrid*945+100] = 2.0E0*I_ESP_Ix4yz_G2xyz_a;
    abcd[iGrid*945+101] = 2.0E0*I_ESP_Ix3y2z_G2xyz_a;
    abcd[iGrid*945+102] = 2.0E0*I_ESP_Ix2y3z_G2xyz_a;
    abcd[iGrid*945+103] = 2.0E0*I_ESP_Ixy4z_G2xyz_a;
    abcd[iGrid*945+104] = 2.0E0*I_ESP_Ix5z_G2xyz_a;
    abcd[iGrid*945+105] = 2.0E0*I_ESP_I6x_G2x2z_a-5*I_ESP_G4x_G2x2z;
    abcd[iGrid*945+106] = 2.0E0*I_ESP_I5xy_G2x2z_a-4*I_ESP_G3xy_G2x2z;
    abcd[iGrid*945+107] = 2.0E0*I_ESP_I5xz_G2x2z_a-4*I_ESP_G3xz_G2x2z;
    abcd[iGrid*945+108] = 2.0E0*I_ESP_I4x2y_G2x2z_a-3*I_ESP_G2x2y_G2x2z;
    abcd[iGrid*945+109] = 2.0E0*I_ESP_I4xyz_G2x2z_a-3*I_ESP_G2xyz_G2x2z;
    abcd[iGrid*945+110] = 2.0E0*I_ESP_I4x2z_G2x2z_a-3*I_ESP_G2x2z_G2x2z;
    abcd[iGrid*945+111] = 2.0E0*I_ESP_I3x3y_G2x2z_a-2*I_ESP_Gx3y_G2x2z;
    abcd[iGrid*945+112] = 2.0E0*I_ESP_I3x2yz_G2x2z_a-2*I_ESP_Gx2yz_G2x2z;
    abcd[iGrid*945+113] = 2.0E0*I_ESP_I3xy2z_G2x2z_a-2*I_ESP_Gxy2z_G2x2z;
    abcd[iGrid*945+114] = 2.0E0*I_ESP_I3x3z_G2x2z_a-2*I_ESP_Gx3z_G2x2z;
    abcd[iGrid*945+115] = 2.0E0*I_ESP_I2x4y_G2x2z_a-1*I_ESP_G4y_G2x2z;
    abcd[iGrid*945+116] = 2.0E0*I_ESP_I2x3yz_G2x2z_a-1*I_ESP_G3yz_G2x2z;
    abcd[iGrid*945+117] = 2.0E0*I_ESP_I2x2y2z_G2x2z_a-1*I_ESP_G2y2z_G2x2z;
    abcd[iGrid*945+118] = 2.0E0*I_ESP_I2xy3z_G2x2z_a-1*I_ESP_Gy3z_G2x2z;
    abcd[iGrid*945+119] = 2.0E0*I_ESP_I2x4z_G2x2z_a-1*I_ESP_G4z_G2x2z;
    abcd[iGrid*945+120] = 2.0E0*I_ESP_Ix5y_G2x2z_a;
    abcd[iGrid*945+121] = 2.0E0*I_ESP_Ix4yz_G2x2z_a;
    abcd[iGrid*945+122] = 2.0E0*I_ESP_Ix3y2z_G2x2z_a;
    abcd[iGrid*945+123] = 2.0E0*I_ESP_Ix2y3z_G2x2z_a;
    abcd[iGrid*945+124] = 2.0E0*I_ESP_Ixy4z_G2x2z_a;
    abcd[iGrid*945+125] = 2.0E0*I_ESP_Ix5z_G2x2z_a;
    abcd[iGrid*945+126] = 2.0E0*I_ESP_I6x_Gx3y_a-5*I_ESP_G4x_Gx3y;
    abcd[iGrid*945+127] = 2.0E0*I_ESP_I5xy_Gx3y_a-4*I_ESP_G3xy_Gx3y;
    abcd[iGrid*945+128] = 2.0E0*I_ESP_I5xz_Gx3y_a-4*I_ESP_G3xz_Gx3y;
    abcd[iGrid*945+129] = 2.0E0*I_ESP_I4x2y_Gx3y_a-3*I_ESP_G2x2y_Gx3y;
    abcd[iGrid*945+130] = 2.0E0*I_ESP_I4xyz_Gx3y_a-3*I_ESP_G2xyz_Gx3y;
    abcd[iGrid*945+131] = 2.0E0*I_ESP_I4x2z_Gx3y_a-3*I_ESP_G2x2z_Gx3y;
    abcd[iGrid*945+132] = 2.0E0*I_ESP_I3x3y_Gx3y_a-2*I_ESP_Gx3y_Gx3y;
    abcd[iGrid*945+133] = 2.0E0*I_ESP_I3x2yz_Gx3y_a-2*I_ESP_Gx2yz_Gx3y;
    abcd[iGrid*945+134] = 2.0E0*I_ESP_I3xy2z_Gx3y_a-2*I_ESP_Gxy2z_Gx3y;
    abcd[iGrid*945+135] = 2.0E0*I_ESP_I3x3z_Gx3y_a-2*I_ESP_Gx3z_Gx3y;
    abcd[iGrid*945+136] = 2.0E0*I_ESP_I2x4y_Gx3y_a-1*I_ESP_G4y_Gx3y;
    abcd[iGrid*945+137] = 2.0E0*I_ESP_I2x3yz_Gx3y_a-1*I_ESP_G3yz_Gx3y;
    abcd[iGrid*945+138] = 2.0E0*I_ESP_I2x2y2z_Gx3y_a-1*I_ESP_G2y2z_Gx3y;
    abcd[iGrid*945+139] = 2.0E0*I_ESP_I2xy3z_Gx3y_a-1*I_ESP_Gy3z_Gx3y;
    abcd[iGrid*945+140] = 2.0E0*I_ESP_I2x4z_Gx3y_a-1*I_ESP_G4z_Gx3y;
    abcd[iGrid*945+141] = 2.0E0*I_ESP_Ix5y_Gx3y_a;
    abcd[iGrid*945+142] = 2.0E0*I_ESP_Ix4yz_Gx3y_a;
    abcd[iGrid*945+143] = 2.0E0*I_ESP_Ix3y2z_Gx3y_a;
    abcd[iGrid*945+144] = 2.0E0*I_ESP_Ix2y3z_Gx3y_a;
    abcd[iGrid*945+145] = 2.0E0*I_ESP_Ixy4z_Gx3y_a;
    abcd[iGrid*945+146] = 2.0E0*I_ESP_Ix5z_Gx3y_a;
    abcd[iGrid*945+147] = 2.0E0*I_ESP_I6x_Gx2yz_a-5*I_ESP_G4x_Gx2yz;
    abcd[iGrid*945+148] = 2.0E0*I_ESP_I5xy_Gx2yz_a-4*I_ESP_G3xy_Gx2yz;
    abcd[iGrid*945+149] = 2.0E0*I_ESP_I5xz_Gx2yz_a-4*I_ESP_G3xz_Gx2yz;
    abcd[iGrid*945+150] = 2.0E0*I_ESP_I4x2y_Gx2yz_a-3*I_ESP_G2x2y_Gx2yz;
    abcd[iGrid*945+151] = 2.0E0*I_ESP_I4xyz_Gx2yz_a-3*I_ESP_G2xyz_Gx2yz;
    abcd[iGrid*945+152] = 2.0E0*I_ESP_I4x2z_Gx2yz_a-3*I_ESP_G2x2z_Gx2yz;
    abcd[iGrid*945+153] = 2.0E0*I_ESP_I3x3y_Gx2yz_a-2*I_ESP_Gx3y_Gx2yz;
    abcd[iGrid*945+154] = 2.0E0*I_ESP_I3x2yz_Gx2yz_a-2*I_ESP_Gx2yz_Gx2yz;
    abcd[iGrid*945+155] = 2.0E0*I_ESP_I3xy2z_Gx2yz_a-2*I_ESP_Gxy2z_Gx2yz;
    abcd[iGrid*945+156] = 2.0E0*I_ESP_I3x3z_Gx2yz_a-2*I_ESP_Gx3z_Gx2yz;
    abcd[iGrid*945+157] = 2.0E0*I_ESP_I2x4y_Gx2yz_a-1*I_ESP_G4y_Gx2yz;
    abcd[iGrid*945+158] = 2.0E0*I_ESP_I2x3yz_Gx2yz_a-1*I_ESP_G3yz_Gx2yz;
    abcd[iGrid*945+159] = 2.0E0*I_ESP_I2x2y2z_Gx2yz_a-1*I_ESP_G2y2z_Gx2yz;
    abcd[iGrid*945+160] = 2.0E0*I_ESP_I2xy3z_Gx2yz_a-1*I_ESP_Gy3z_Gx2yz;
    abcd[iGrid*945+161] = 2.0E0*I_ESP_I2x4z_Gx2yz_a-1*I_ESP_G4z_Gx2yz;
    abcd[iGrid*945+162] = 2.0E0*I_ESP_Ix5y_Gx2yz_a;
    abcd[iGrid*945+163] = 2.0E0*I_ESP_Ix4yz_Gx2yz_a;
    abcd[iGrid*945+164] = 2.0E0*I_ESP_Ix3y2z_Gx2yz_a;
    abcd[iGrid*945+165] = 2.0E0*I_ESP_Ix2y3z_Gx2yz_a;
    abcd[iGrid*945+166] = 2.0E0*I_ESP_Ixy4z_Gx2yz_a;
    abcd[iGrid*945+167] = 2.0E0*I_ESP_Ix5z_Gx2yz_a;
    abcd[iGrid*945+168] = 2.0E0*I_ESP_I6x_Gxy2z_a-5*I_ESP_G4x_Gxy2z;
    abcd[iGrid*945+169] = 2.0E0*I_ESP_I5xy_Gxy2z_a-4*I_ESP_G3xy_Gxy2z;
    abcd[iGrid*945+170] = 2.0E0*I_ESP_I5xz_Gxy2z_a-4*I_ESP_G3xz_Gxy2z;
    abcd[iGrid*945+171] = 2.0E0*I_ESP_I4x2y_Gxy2z_a-3*I_ESP_G2x2y_Gxy2z;
    abcd[iGrid*945+172] = 2.0E0*I_ESP_I4xyz_Gxy2z_a-3*I_ESP_G2xyz_Gxy2z;
    abcd[iGrid*945+173] = 2.0E0*I_ESP_I4x2z_Gxy2z_a-3*I_ESP_G2x2z_Gxy2z;
    abcd[iGrid*945+174] = 2.0E0*I_ESP_I3x3y_Gxy2z_a-2*I_ESP_Gx3y_Gxy2z;
    abcd[iGrid*945+175] = 2.0E0*I_ESP_I3x2yz_Gxy2z_a-2*I_ESP_Gx2yz_Gxy2z;
    abcd[iGrid*945+176] = 2.0E0*I_ESP_I3xy2z_Gxy2z_a-2*I_ESP_Gxy2z_Gxy2z;
    abcd[iGrid*945+177] = 2.0E0*I_ESP_I3x3z_Gxy2z_a-2*I_ESP_Gx3z_Gxy2z;
    abcd[iGrid*945+178] = 2.0E0*I_ESP_I2x4y_Gxy2z_a-1*I_ESP_G4y_Gxy2z;
    abcd[iGrid*945+179] = 2.0E0*I_ESP_I2x3yz_Gxy2z_a-1*I_ESP_G3yz_Gxy2z;
    abcd[iGrid*945+180] = 2.0E0*I_ESP_I2x2y2z_Gxy2z_a-1*I_ESP_G2y2z_Gxy2z;
    abcd[iGrid*945+181] = 2.0E0*I_ESP_I2xy3z_Gxy2z_a-1*I_ESP_Gy3z_Gxy2z;
    abcd[iGrid*945+182] = 2.0E0*I_ESP_I2x4z_Gxy2z_a-1*I_ESP_G4z_Gxy2z;
    abcd[iGrid*945+183] = 2.0E0*I_ESP_Ix5y_Gxy2z_a;
    abcd[iGrid*945+184] = 2.0E0*I_ESP_Ix4yz_Gxy2z_a;
    abcd[iGrid*945+185] = 2.0E0*I_ESP_Ix3y2z_Gxy2z_a;
    abcd[iGrid*945+186] = 2.0E0*I_ESP_Ix2y3z_Gxy2z_a;
    abcd[iGrid*945+187] = 2.0E0*I_ESP_Ixy4z_Gxy2z_a;
    abcd[iGrid*945+188] = 2.0E0*I_ESP_Ix5z_Gxy2z_a;
    abcd[iGrid*945+189] = 2.0E0*I_ESP_I6x_Gx3z_a-5*I_ESP_G4x_Gx3z;
    abcd[iGrid*945+190] = 2.0E0*I_ESP_I5xy_Gx3z_a-4*I_ESP_G3xy_Gx3z;
    abcd[iGrid*945+191] = 2.0E0*I_ESP_I5xz_Gx3z_a-4*I_ESP_G3xz_Gx3z;
    abcd[iGrid*945+192] = 2.0E0*I_ESP_I4x2y_Gx3z_a-3*I_ESP_G2x2y_Gx3z;
    abcd[iGrid*945+193] = 2.0E0*I_ESP_I4xyz_Gx3z_a-3*I_ESP_G2xyz_Gx3z;
    abcd[iGrid*945+194] = 2.0E0*I_ESP_I4x2z_Gx3z_a-3*I_ESP_G2x2z_Gx3z;
    abcd[iGrid*945+195] = 2.0E0*I_ESP_I3x3y_Gx3z_a-2*I_ESP_Gx3y_Gx3z;
    abcd[iGrid*945+196] = 2.0E0*I_ESP_I3x2yz_Gx3z_a-2*I_ESP_Gx2yz_Gx3z;
    abcd[iGrid*945+197] = 2.0E0*I_ESP_I3xy2z_Gx3z_a-2*I_ESP_Gxy2z_Gx3z;
    abcd[iGrid*945+198] = 2.0E0*I_ESP_I3x3z_Gx3z_a-2*I_ESP_Gx3z_Gx3z;
    abcd[iGrid*945+199] = 2.0E0*I_ESP_I2x4y_Gx3z_a-1*I_ESP_G4y_Gx3z;
    abcd[iGrid*945+200] = 2.0E0*I_ESP_I2x3yz_Gx3z_a-1*I_ESP_G3yz_Gx3z;
    abcd[iGrid*945+201] = 2.0E0*I_ESP_I2x2y2z_Gx3z_a-1*I_ESP_G2y2z_Gx3z;
    abcd[iGrid*945+202] = 2.0E0*I_ESP_I2xy3z_Gx3z_a-1*I_ESP_Gy3z_Gx3z;
    abcd[iGrid*945+203] = 2.0E0*I_ESP_I2x4z_Gx3z_a-1*I_ESP_G4z_Gx3z;
    abcd[iGrid*945+204] = 2.0E0*I_ESP_Ix5y_Gx3z_a;
    abcd[iGrid*945+205] = 2.0E0*I_ESP_Ix4yz_Gx3z_a;
    abcd[iGrid*945+206] = 2.0E0*I_ESP_Ix3y2z_Gx3z_a;
    abcd[iGrid*945+207] = 2.0E0*I_ESP_Ix2y3z_Gx3z_a;
    abcd[iGrid*945+208] = 2.0E0*I_ESP_Ixy4z_Gx3z_a;
    abcd[iGrid*945+209] = 2.0E0*I_ESP_Ix5z_Gx3z_a;
    abcd[iGrid*945+210] = 2.0E0*I_ESP_I6x_G4y_a-5*I_ESP_G4x_G4y;
    abcd[iGrid*945+211] = 2.0E0*I_ESP_I5xy_G4y_a-4*I_ESP_G3xy_G4y;
    abcd[iGrid*945+212] = 2.0E0*I_ESP_I5xz_G4y_a-4*I_ESP_G3xz_G4y;
    abcd[iGrid*945+213] = 2.0E0*I_ESP_I4x2y_G4y_a-3*I_ESP_G2x2y_G4y;
    abcd[iGrid*945+214] = 2.0E0*I_ESP_I4xyz_G4y_a-3*I_ESP_G2xyz_G4y;
    abcd[iGrid*945+215] = 2.0E0*I_ESP_I4x2z_G4y_a-3*I_ESP_G2x2z_G4y;
    abcd[iGrid*945+216] = 2.0E0*I_ESP_I3x3y_G4y_a-2*I_ESP_Gx3y_G4y;
    abcd[iGrid*945+217] = 2.0E0*I_ESP_I3x2yz_G4y_a-2*I_ESP_Gx2yz_G4y;
    abcd[iGrid*945+218] = 2.0E0*I_ESP_I3xy2z_G4y_a-2*I_ESP_Gxy2z_G4y;
    abcd[iGrid*945+219] = 2.0E0*I_ESP_I3x3z_G4y_a-2*I_ESP_Gx3z_G4y;
    abcd[iGrid*945+220] = 2.0E0*I_ESP_I2x4y_G4y_a-1*I_ESP_G4y_G4y;
    abcd[iGrid*945+221] = 2.0E0*I_ESP_I2x3yz_G4y_a-1*I_ESP_G3yz_G4y;
    abcd[iGrid*945+222] = 2.0E0*I_ESP_I2x2y2z_G4y_a-1*I_ESP_G2y2z_G4y;
    abcd[iGrid*945+223] = 2.0E0*I_ESP_I2xy3z_G4y_a-1*I_ESP_Gy3z_G4y;
    abcd[iGrid*945+224] = 2.0E0*I_ESP_I2x4z_G4y_a-1*I_ESP_G4z_G4y;
    abcd[iGrid*945+225] = 2.0E0*I_ESP_Ix5y_G4y_a;
    abcd[iGrid*945+226] = 2.0E0*I_ESP_Ix4yz_G4y_a;
    abcd[iGrid*945+227] = 2.0E0*I_ESP_Ix3y2z_G4y_a;
    abcd[iGrid*945+228] = 2.0E0*I_ESP_Ix2y3z_G4y_a;
    abcd[iGrid*945+229] = 2.0E0*I_ESP_Ixy4z_G4y_a;
    abcd[iGrid*945+230] = 2.0E0*I_ESP_Ix5z_G4y_a;
    abcd[iGrid*945+231] = 2.0E0*I_ESP_I6x_G3yz_a-5*I_ESP_G4x_G3yz;
    abcd[iGrid*945+232] = 2.0E0*I_ESP_I5xy_G3yz_a-4*I_ESP_G3xy_G3yz;
    abcd[iGrid*945+233] = 2.0E0*I_ESP_I5xz_G3yz_a-4*I_ESP_G3xz_G3yz;
    abcd[iGrid*945+234] = 2.0E0*I_ESP_I4x2y_G3yz_a-3*I_ESP_G2x2y_G3yz;
    abcd[iGrid*945+235] = 2.0E0*I_ESP_I4xyz_G3yz_a-3*I_ESP_G2xyz_G3yz;
    abcd[iGrid*945+236] = 2.0E0*I_ESP_I4x2z_G3yz_a-3*I_ESP_G2x2z_G3yz;
    abcd[iGrid*945+237] = 2.0E0*I_ESP_I3x3y_G3yz_a-2*I_ESP_Gx3y_G3yz;
    abcd[iGrid*945+238] = 2.0E0*I_ESP_I3x2yz_G3yz_a-2*I_ESP_Gx2yz_G3yz;
    abcd[iGrid*945+239] = 2.0E0*I_ESP_I3xy2z_G3yz_a-2*I_ESP_Gxy2z_G3yz;
    abcd[iGrid*945+240] = 2.0E0*I_ESP_I3x3z_G3yz_a-2*I_ESP_Gx3z_G3yz;
    abcd[iGrid*945+241] = 2.0E0*I_ESP_I2x4y_G3yz_a-1*I_ESP_G4y_G3yz;
    abcd[iGrid*945+242] = 2.0E0*I_ESP_I2x3yz_G3yz_a-1*I_ESP_G3yz_G3yz;
    abcd[iGrid*945+243] = 2.0E0*I_ESP_I2x2y2z_G3yz_a-1*I_ESP_G2y2z_G3yz;
    abcd[iGrid*945+244] = 2.0E0*I_ESP_I2xy3z_G3yz_a-1*I_ESP_Gy3z_G3yz;
    abcd[iGrid*945+245] = 2.0E0*I_ESP_I2x4z_G3yz_a-1*I_ESP_G4z_G3yz;
    abcd[iGrid*945+246] = 2.0E0*I_ESP_Ix5y_G3yz_a;
    abcd[iGrid*945+247] = 2.0E0*I_ESP_Ix4yz_G3yz_a;
    abcd[iGrid*945+248] = 2.0E0*I_ESP_Ix3y2z_G3yz_a;
    abcd[iGrid*945+249] = 2.0E0*I_ESP_Ix2y3z_G3yz_a;
    abcd[iGrid*945+250] = 2.0E0*I_ESP_Ixy4z_G3yz_a;
    abcd[iGrid*945+251] = 2.0E0*I_ESP_Ix5z_G3yz_a;
    abcd[iGrid*945+252] = 2.0E0*I_ESP_I6x_G2y2z_a-5*I_ESP_G4x_G2y2z;
    abcd[iGrid*945+253] = 2.0E0*I_ESP_I5xy_G2y2z_a-4*I_ESP_G3xy_G2y2z;
    abcd[iGrid*945+254] = 2.0E0*I_ESP_I5xz_G2y2z_a-4*I_ESP_G3xz_G2y2z;
    abcd[iGrid*945+255] = 2.0E0*I_ESP_I4x2y_G2y2z_a-3*I_ESP_G2x2y_G2y2z;
    abcd[iGrid*945+256] = 2.0E0*I_ESP_I4xyz_G2y2z_a-3*I_ESP_G2xyz_G2y2z;
    abcd[iGrid*945+257] = 2.0E0*I_ESP_I4x2z_G2y2z_a-3*I_ESP_G2x2z_G2y2z;
    abcd[iGrid*945+258] = 2.0E0*I_ESP_I3x3y_G2y2z_a-2*I_ESP_Gx3y_G2y2z;
    abcd[iGrid*945+259] = 2.0E0*I_ESP_I3x2yz_G2y2z_a-2*I_ESP_Gx2yz_G2y2z;
    abcd[iGrid*945+260] = 2.0E0*I_ESP_I3xy2z_G2y2z_a-2*I_ESP_Gxy2z_G2y2z;
    abcd[iGrid*945+261] = 2.0E0*I_ESP_I3x3z_G2y2z_a-2*I_ESP_Gx3z_G2y2z;
    abcd[iGrid*945+262] = 2.0E0*I_ESP_I2x4y_G2y2z_a-1*I_ESP_G4y_G2y2z;
    abcd[iGrid*945+263] = 2.0E0*I_ESP_I2x3yz_G2y2z_a-1*I_ESP_G3yz_G2y2z;
    abcd[iGrid*945+264] = 2.0E0*I_ESP_I2x2y2z_G2y2z_a-1*I_ESP_G2y2z_G2y2z;
    abcd[iGrid*945+265] = 2.0E0*I_ESP_I2xy3z_G2y2z_a-1*I_ESP_Gy3z_G2y2z;
    abcd[iGrid*945+266] = 2.0E0*I_ESP_I2x4z_G2y2z_a-1*I_ESP_G4z_G2y2z;
    abcd[iGrid*945+267] = 2.0E0*I_ESP_Ix5y_G2y2z_a;
    abcd[iGrid*945+268] = 2.0E0*I_ESP_Ix4yz_G2y2z_a;
    abcd[iGrid*945+269] = 2.0E0*I_ESP_Ix3y2z_G2y2z_a;
    abcd[iGrid*945+270] = 2.0E0*I_ESP_Ix2y3z_G2y2z_a;
    abcd[iGrid*945+271] = 2.0E0*I_ESP_Ixy4z_G2y2z_a;
    abcd[iGrid*945+272] = 2.0E0*I_ESP_Ix5z_G2y2z_a;
    abcd[iGrid*945+273] = 2.0E0*I_ESP_I6x_Gy3z_a-5*I_ESP_G4x_Gy3z;
    abcd[iGrid*945+274] = 2.0E0*I_ESP_I5xy_Gy3z_a-4*I_ESP_G3xy_Gy3z;
    abcd[iGrid*945+275] = 2.0E0*I_ESP_I5xz_Gy3z_a-4*I_ESP_G3xz_Gy3z;
    abcd[iGrid*945+276] = 2.0E0*I_ESP_I4x2y_Gy3z_a-3*I_ESP_G2x2y_Gy3z;
    abcd[iGrid*945+277] = 2.0E0*I_ESP_I4xyz_Gy3z_a-3*I_ESP_G2xyz_Gy3z;
    abcd[iGrid*945+278] = 2.0E0*I_ESP_I4x2z_Gy3z_a-3*I_ESP_G2x2z_Gy3z;
    abcd[iGrid*945+279] = 2.0E0*I_ESP_I3x3y_Gy3z_a-2*I_ESP_Gx3y_Gy3z;
    abcd[iGrid*945+280] = 2.0E0*I_ESP_I3x2yz_Gy3z_a-2*I_ESP_Gx2yz_Gy3z;
    abcd[iGrid*945+281] = 2.0E0*I_ESP_I3xy2z_Gy3z_a-2*I_ESP_Gxy2z_Gy3z;
    abcd[iGrid*945+282] = 2.0E0*I_ESP_I3x3z_Gy3z_a-2*I_ESP_Gx3z_Gy3z;
    abcd[iGrid*945+283] = 2.0E0*I_ESP_I2x4y_Gy3z_a-1*I_ESP_G4y_Gy3z;
    abcd[iGrid*945+284] = 2.0E0*I_ESP_I2x3yz_Gy3z_a-1*I_ESP_G3yz_Gy3z;
    abcd[iGrid*945+285] = 2.0E0*I_ESP_I2x2y2z_Gy3z_a-1*I_ESP_G2y2z_Gy3z;
    abcd[iGrid*945+286] = 2.0E0*I_ESP_I2xy3z_Gy3z_a-1*I_ESP_Gy3z_Gy3z;
    abcd[iGrid*945+287] = 2.0E0*I_ESP_I2x4z_Gy3z_a-1*I_ESP_G4z_Gy3z;
    abcd[iGrid*945+288] = 2.0E0*I_ESP_Ix5y_Gy3z_a;
    abcd[iGrid*945+289] = 2.0E0*I_ESP_Ix4yz_Gy3z_a;
    abcd[iGrid*945+290] = 2.0E0*I_ESP_Ix3y2z_Gy3z_a;
    abcd[iGrid*945+291] = 2.0E0*I_ESP_Ix2y3z_Gy3z_a;
    abcd[iGrid*945+292] = 2.0E0*I_ESP_Ixy4z_Gy3z_a;
    abcd[iGrid*945+293] = 2.0E0*I_ESP_Ix5z_Gy3z_a;
    abcd[iGrid*945+294] = 2.0E0*I_ESP_I6x_G4z_a-5*I_ESP_G4x_G4z;
    abcd[iGrid*945+295] = 2.0E0*I_ESP_I5xy_G4z_a-4*I_ESP_G3xy_G4z;
    abcd[iGrid*945+296] = 2.0E0*I_ESP_I5xz_G4z_a-4*I_ESP_G3xz_G4z;
    abcd[iGrid*945+297] = 2.0E0*I_ESP_I4x2y_G4z_a-3*I_ESP_G2x2y_G4z;
    abcd[iGrid*945+298] = 2.0E0*I_ESP_I4xyz_G4z_a-3*I_ESP_G2xyz_G4z;
    abcd[iGrid*945+299] = 2.0E0*I_ESP_I4x2z_G4z_a-3*I_ESP_G2x2z_G4z;
    abcd[iGrid*945+300] = 2.0E0*I_ESP_I3x3y_G4z_a-2*I_ESP_Gx3y_G4z;
    abcd[iGrid*945+301] = 2.0E0*I_ESP_I3x2yz_G4z_a-2*I_ESP_Gx2yz_G4z;
    abcd[iGrid*945+302] = 2.0E0*I_ESP_I3xy2z_G4z_a-2*I_ESP_Gxy2z_G4z;
    abcd[iGrid*945+303] = 2.0E0*I_ESP_I3x3z_G4z_a-2*I_ESP_Gx3z_G4z;
    abcd[iGrid*945+304] = 2.0E0*I_ESP_I2x4y_G4z_a-1*I_ESP_G4y_G4z;
    abcd[iGrid*945+305] = 2.0E0*I_ESP_I2x3yz_G4z_a-1*I_ESP_G3yz_G4z;
    abcd[iGrid*945+306] = 2.0E0*I_ESP_I2x2y2z_G4z_a-1*I_ESP_G2y2z_G4z;
    abcd[iGrid*945+307] = 2.0E0*I_ESP_I2xy3z_G4z_a-1*I_ESP_Gy3z_G4z;
    abcd[iGrid*945+308] = 2.0E0*I_ESP_I2x4z_G4z_a-1*I_ESP_G4z_G4z;
    abcd[iGrid*945+309] = 2.0E0*I_ESP_Ix5y_G4z_a;
    abcd[iGrid*945+310] = 2.0E0*I_ESP_Ix4yz_G4z_a;
    abcd[iGrid*945+311] = 2.0E0*I_ESP_Ix3y2z_G4z_a;
    abcd[iGrid*945+312] = 2.0E0*I_ESP_Ix2y3z_G4z_a;
    abcd[iGrid*945+313] = 2.0E0*I_ESP_Ixy4z_G4z_a;
    abcd[iGrid*945+314] = 2.0E0*I_ESP_Ix5z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G_a
     * RHS shell quartet name: SQ_ESP_G_G
     ************************************************************/
    abcd[iGrid*945+315] = 2.0E0*I_ESP_I5xy_G4x_a;
    abcd[iGrid*945+316] = 2.0E0*I_ESP_I4x2y_G4x_a-1*I_ESP_G4x_G4x;
    abcd[iGrid*945+317] = 2.0E0*I_ESP_I4xyz_G4x_a;
    abcd[iGrid*945+318] = 2.0E0*I_ESP_I3x3y_G4x_a-2*I_ESP_G3xy_G4x;
    abcd[iGrid*945+319] = 2.0E0*I_ESP_I3x2yz_G4x_a-1*I_ESP_G3xz_G4x;
    abcd[iGrid*945+320] = 2.0E0*I_ESP_I3xy2z_G4x_a;
    abcd[iGrid*945+321] = 2.0E0*I_ESP_I2x4y_G4x_a-3*I_ESP_G2x2y_G4x;
    abcd[iGrid*945+322] = 2.0E0*I_ESP_I2x3yz_G4x_a-2*I_ESP_G2xyz_G4x;
    abcd[iGrid*945+323] = 2.0E0*I_ESP_I2x2y2z_G4x_a-1*I_ESP_G2x2z_G4x;
    abcd[iGrid*945+324] = 2.0E0*I_ESP_I2xy3z_G4x_a;
    abcd[iGrid*945+325] = 2.0E0*I_ESP_Ix5y_G4x_a-4*I_ESP_Gx3y_G4x;
    abcd[iGrid*945+326] = 2.0E0*I_ESP_Ix4yz_G4x_a-3*I_ESP_Gx2yz_G4x;
    abcd[iGrid*945+327] = 2.0E0*I_ESP_Ix3y2z_G4x_a-2*I_ESP_Gxy2z_G4x;
    abcd[iGrid*945+328] = 2.0E0*I_ESP_Ix2y3z_G4x_a-1*I_ESP_Gx3z_G4x;
    abcd[iGrid*945+329] = 2.0E0*I_ESP_Ixy4z_G4x_a;
    abcd[iGrid*945+330] = 2.0E0*I_ESP_I6y_G4x_a-5*I_ESP_G4y_G4x;
    abcd[iGrid*945+331] = 2.0E0*I_ESP_I5yz_G4x_a-4*I_ESP_G3yz_G4x;
    abcd[iGrid*945+332] = 2.0E0*I_ESP_I4y2z_G4x_a-3*I_ESP_G2y2z_G4x;
    abcd[iGrid*945+333] = 2.0E0*I_ESP_I3y3z_G4x_a-2*I_ESP_Gy3z_G4x;
    abcd[iGrid*945+334] = 2.0E0*I_ESP_I2y4z_G4x_a-1*I_ESP_G4z_G4x;
    abcd[iGrid*945+335] = 2.0E0*I_ESP_Iy5z_G4x_a;
    abcd[iGrid*945+336] = 2.0E0*I_ESP_I5xy_G3xy_a;
    abcd[iGrid*945+337] = 2.0E0*I_ESP_I4x2y_G3xy_a-1*I_ESP_G4x_G3xy;
    abcd[iGrid*945+338] = 2.0E0*I_ESP_I4xyz_G3xy_a;
    abcd[iGrid*945+339] = 2.0E0*I_ESP_I3x3y_G3xy_a-2*I_ESP_G3xy_G3xy;
    abcd[iGrid*945+340] = 2.0E0*I_ESP_I3x2yz_G3xy_a-1*I_ESP_G3xz_G3xy;
    abcd[iGrid*945+341] = 2.0E0*I_ESP_I3xy2z_G3xy_a;
    abcd[iGrid*945+342] = 2.0E0*I_ESP_I2x4y_G3xy_a-3*I_ESP_G2x2y_G3xy;
    abcd[iGrid*945+343] = 2.0E0*I_ESP_I2x3yz_G3xy_a-2*I_ESP_G2xyz_G3xy;
    abcd[iGrid*945+344] = 2.0E0*I_ESP_I2x2y2z_G3xy_a-1*I_ESP_G2x2z_G3xy;
    abcd[iGrid*945+345] = 2.0E0*I_ESP_I2xy3z_G3xy_a;
    abcd[iGrid*945+346] = 2.0E0*I_ESP_Ix5y_G3xy_a-4*I_ESP_Gx3y_G3xy;
    abcd[iGrid*945+347] = 2.0E0*I_ESP_Ix4yz_G3xy_a-3*I_ESP_Gx2yz_G3xy;
    abcd[iGrid*945+348] = 2.0E0*I_ESP_Ix3y2z_G3xy_a-2*I_ESP_Gxy2z_G3xy;
    abcd[iGrid*945+349] = 2.0E0*I_ESP_Ix2y3z_G3xy_a-1*I_ESP_Gx3z_G3xy;
    abcd[iGrid*945+350] = 2.0E0*I_ESP_Ixy4z_G3xy_a;
    abcd[iGrid*945+351] = 2.0E0*I_ESP_I6y_G3xy_a-5*I_ESP_G4y_G3xy;
    abcd[iGrid*945+352] = 2.0E0*I_ESP_I5yz_G3xy_a-4*I_ESP_G3yz_G3xy;
    abcd[iGrid*945+353] = 2.0E0*I_ESP_I4y2z_G3xy_a-3*I_ESP_G2y2z_G3xy;
    abcd[iGrid*945+354] = 2.0E0*I_ESP_I3y3z_G3xy_a-2*I_ESP_Gy3z_G3xy;
    abcd[iGrid*945+355] = 2.0E0*I_ESP_I2y4z_G3xy_a-1*I_ESP_G4z_G3xy;
    abcd[iGrid*945+356] = 2.0E0*I_ESP_Iy5z_G3xy_a;
    abcd[iGrid*945+357] = 2.0E0*I_ESP_I5xy_G3xz_a;
    abcd[iGrid*945+358] = 2.0E0*I_ESP_I4x2y_G3xz_a-1*I_ESP_G4x_G3xz;
    abcd[iGrid*945+359] = 2.0E0*I_ESP_I4xyz_G3xz_a;
    abcd[iGrid*945+360] = 2.0E0*I_ESP_I3x3y_G3xz_a-2*I_ESP_G3xy_G3xz;
    abcd[iGrid*945+361] = 2.0E0*I_ESP_I3x2yz_G3xz_a-1*I_ESP_G3xz_G3xz;
    abcd[iGrid*945+362] = 2.0E0*I_ESP_I3xy2z_G3xz_a;
    abcd[iGrid*945+363] = 2.0E0*I_ESP_I2x4y_G3xz_a-3*I_ESP_G2x2y_G3xz;
    abcd[iGrid*945+364] = 2.0E0*I_ESP_I2x3yz_G3xz_a-2*I_ESP_G2xyz_G3xz;
    abcd[iGrid*945+365] = 2.0E0*I_ESP_I2x2y2z_G3xz_a-1*I_ESP_G2x2z_G3xz;
    abcd[iGrid*945+366] = 2.0E0*I_ESP_I2xy3z_G3xz_a;
    abcd[iGrid*945+367] = 2.0E0*I_ESP_Ix5y_G3xz_a-4*I_ESP_Gx3y_G3xz;
    abcd[iGrid*945+368] = 2.0E0*I_ESP_Ix4yz_G3xz_a-3*I_ESP_Gx2yz_G3xz;
    abcd[iGrid*945+369] = 2.0E0*I_ESP_Ix3y2z_G3xz_a-2*I_ESP_Gxy2z_G3xz;
    abcd[iGrid*945+370] = 2.0E0*I_ESP_Ix2y3z_G3xz_a-1*I_ESP_Gx3z_G3xz;
    abcd[iGrid*945+371] = 2.0E0*I_ESP_Ixy4z_G3xz_a;
    abcd[iGrid*945+372] = 2.0E0*I_ESP_I6y_G3xz_a-5*I_ESP_G4y_G3xz;
    abcd[iGrid*945+373] = 2.0E0*I_ESP_I5yz_G3xz_a-4*I_ESP_G3yz_G3xz;
    abcd[iGrid*945+374] = 2.0E0*I_ESP_I4y2z_G3xz_a-3*I_ESP_G2y2z_G3xz;
    abcd[iGrid*945+375] = 2.0E0*I_ESP_I3y3z_G3xz_a-2*I_ESP_Gy3z_G3xz;
    abcd[iGrid*945+376] = 2.0E0*I_ESP_I2y4z_G3xz_a-1*I_ESP_G4z_G3xz;
    abcd[iGrid*945+377] = 2.0E0*I_ESP_Iy5z_G3xz_a;
    abcd[iGrid*945+378] = 2.0E0*I_ESP_I5xy_G2x2y_a;
    abcd[iGrid*945+379] = 2.0E0*I_ESP_I4x2y_G2x2y_a-1*I_ESP_G4x_G2x2y;
    abcd[iGrid*945+380] = 2.0E0*I_ESP_I4xyz_G2x2y_a;
    abcd[iGrid*945+381] = 2.0E0*I_ESP_I3x3y_G2x2y_a-2*I_ESP_G3xy_G2x2y;
    abcd[iGrid*945+382] = 2.0E0*I_ESP_I3x2yz_G2x2y_a-1*I_ESP_G3xz_G2x2y;
    abcd[iGrid*945+383] = 2.0E0*I_ESP_I3xy2z_G2x2y_a;
    abcd[iGrid*945+384] = 2.0E0*I_ESP_I2x4y_G2x2y_a-3*I_ESP_G2x2y_G2x2y;
    abcd[iGrid*945+385] = 2.0E0*I_ESP_I2x3yz_G2x2y_a-2*I_ESP_G2xyz_G2x2y;
    abcd[iGrid*945+386] = 2.0E0*I_ESP_I2x2y2z_G2x2y_a-1*I_ESP_G2x2z_G2x2y;
    abcd[iGrid*945+387] = 2.0E0*I_ESP_I2xy3z_G2x2y_a;
    abcd[iGrid*945+388] = 2.0E0*I_ESP_Ix5y_G2x2y_a-4*I_ESP_Gx3y_G2x2y;
    abcd[iGrid*945+389] = 2.0E0*I_ESP_Ix4yz_G2x2y_a-3*I_ESP_Gx2yz_G2x2y;
    abcd[iGrid*945+390] = 2.0E0*I_ESP_Ix3y2z_G2x2y_a-2*I_ESP_Gxy2z_G2x2y;
    abcd[iGrid*945+391] = 2.0E0*I_ESP_Ix2y3z_G2x2y_a-1*I_ESP_Gx3z_G2x2y;
    abcd[iGrid*945+392] = 2.0E0*I_ESP_Ixy4z_G2x2y_a;
    abcd[iGrid*945+393] = 2.0E0*I_ESP_I6y_G2x2y_a-5*I_ESP_G4y_G2x2y;
    abcd[iGrid*945+394] = 2.0E0*I_ESP_I5yz_G2x2y_a-4*I_ESP_G3yz_G2x2y;
    abcd[iGrid*945+395] = 2.0E0*I_ESP_I4y2z_G2x2y_a-3*I_ESP_G2y2z_G2x2y;
    abcd[iGrid*945+396] = 2.0E0*I_ESP_I3y3z_G2x2y_a-2*I_ESP_Gy3z_G2x2y;
    abcd[iGrid*945+397] = 2.0E0*I_ESP_I2y4z_G2x2y_a-1*I_ESP_G4z_G2x2y;
    abcd[iGrid*945+398] = 2.0E0*I_ESP_Iy5z_G2x2y_a;
    abcd[iGrid*945+399] = 2.0E0*I_ESP_I5xy_G2xyz_a;
    abcd[iGrid*945+400] = 2.0E0*I_ESP_I4x2y_G2xyz_a-1*I_ESP_G4x_G2xyz;
    abcd[iGrid*945+401] = 2.0E0*I_ESP_I4xyz_G2xyz_a;
    abcd[iGrid*945+402] = 2.0E0*I_ESP_I3x3y_G2xyz_a-2*I_ESP_G3xy_G2xyz;
    abcd[iGrid*945+403] = 2.0E0*I_ESP_I3x2yz_G2xyz_a-1*I_ESP_G3xz_G2xyz;
    abcd[iGrid*945+404] = 2.0E0*I_ESP_I3xy2z_G2xyz_a;
    abcd[iGrid*945+405] = 2.0E0*I_ESP_I2x4y_G2xyz_a-3*I_ESP_G2x2y_G2xyz;
    abcd[iGrid*945+406] = 2.0E0*I_ESP_I2x3yz_G2xyz_a-2*I_ESP_G2xyz_G2xyz;
    abcd[iGrid*945+407] = 2.0E0*I_ESP_I2x2y2z_G2xyz_a-1*I_ESP_G2x2z_G2xyz;
    abcd[iGrid*945+408] = 2.0E0*I_ESP_I2xy3z_G2xyz_a;
    abcd[iGrid*945+409] = 2.0E0*I_ESP_Ix5y_G2xyz_a-4*I_ESP_Gx3y_G2xyz;
    abcd[iGrid*945+410] = 2.0E0*I_ESP_Ix4yz_G2xyz_a-3*I_ESP_Gx2yz_G2xyz;
    abcd[iGrid*945+411] = 2.0E0*I_ESP_Ix3y2z_G2xyz_a-2*I_ESP_Gxy2z_G2xyz;
    abcd[iGrid*945+412] = 2.0E0*I_ESP_Ix2y3z_G2xyz_a-1*I_ESP_Gx3z_G2xyz;
    abcd[iGrid*945+413] = 2.0E0*I_ESP_Ixy4z_G2xyz_a;
    abcd[iGrid*945+414] = 2.0E0*I_ESP_I6y_G2xyz_a-5*I_ESP_G4y_G2xyz;
    abcd[iGrid*945+415] = 2.0E0*I_ESP_I5yz_G2xyz_a-4*I_ESP_G3yz_G2xyz;
    abcd[iGrid*945+416] = 2.0E0*I_ESP_I4y2z_G2xyz_a-3*I_ESP_G2y2z_G2xyz;
    abcd[iGrid*945+417] = 2.0E0*I_ESP_I3y3z_G2xyz_a-2*I_ESP_Gy3z_G2xyz;
    abcd[iGrid*945+418] = 2.0E0*I_ESP_I2y4z_G2xyz_a-1*I_ESP_G4z_G2xyz;
    abcd[iGrid*945+419] = 2.0E0*I_ESP_Iy5z_G2xyz_a;
    abcd[iGrid*945+420] = 2.0E0*I_ESP_I5xy_G2x2z_a;
    abcd[iGrid*945+421] = 2.0E0*I_ESP_I4x2y_G2x2z_a-1*I_ESP_G4x_G2x2z;
    abcd[iGrid*945+422] = 2.0E0*I_ESP_I4xyz_G2x2z_a;
    abcd[iGrid*945+423] = 2.0E0*I_ESP_I3x3y_G2x2z_a-2*I_ESP_G3xy_G2x2z;
    abcd[iGrid*945+424] = 2.0E0*I_ESP_I3x2yz_G2x2z_a-1*I_ESP_G3xz_G2x2z;
    abcd[iGrid*945+425] = 2.0E0*I_ESP_I3xy2z_G2x2z_a;
    abcd[iGrid*945+426] = 2.0E0*I_ESP_I2x4y_G2x2z_a-3*I_ESP_G2x2y_G2x2z;
    abcd[iGrid*945+427] = 2.0E0*I_ESP_I2x3yz_G2x2z_a-2*I_ESP_G2xyz_G2x2z;
    abcd[iGrid*945+428] = 2.0E0*I_ESP_I2x2y2z_G2x2z_a-1*I_ESP_G2x2z_G2x2z;
    abcd[iGrid*945+429] = 2.0E0*I_ESP_I2xy3z_G2x2z_a;
    abcd[iGrid*945+430] = 2.0E0*I_ESP_Ix5y_G2x2z_a-4*I_ESP_Gx3y_G2x2z;
    abcd[iGrid*945+431] = 2.0E0*I_ESP_Ix4yz_G2x2z_a-3*I_ESP_Gx2yz_G2x2z;
    abcd[iGrid*945+432] = 2.0E0*I_ESP_Ix3y2z_G2x2z_a-2*I_ESP_Gxy2z_G2x2z;
    abcd[iGrid*945+433] = 2.0E0*I_ESP_Ix2y3z_G2x2z_a-1*I_ESP_Gx3z_G2x2z;
    abcd[iGrid*945+434] = 2.0E0*I_ESP_Ixy4z_G2x2z_a;
    abcd[iGrid*945+435] = 2.0E0*I_ESP_I6y_G2x2z_a-5*I_ESP_G4y_G2x2z;
    abcd[iGrid*945+436] = 2.0E0*I_ESP_I5yz_G2x2z_a-4*I_ESP_G3yz_G2x2z;
    abcd[iGrid*945+437] = 2.0E0*I_ESP_I4y2z_G2x2z_a-3*I_ESP_G2y2z_G2x2z;
    abcd[iGrid*945+438] = 2.0E0*I_ESP_I3y3z_G2x2z_a-2*I_ESP_Gy3z_G2x2z;
    abcd[iGrid*945+439] = 2.0E0*I_ESP_I2y4z_G2x2z_a-1*I_ESP_G4z_G2x2z;
    abcd[iGrid*945+440] = 2.0E0*I_ESP_Iy5z_G2x2z_a;
    abcd[iGrid*945+441] = 2.0E0*I_ESP_I5xy_Gx3y_a;
    abcd[iGrid*945+442] = 2.0E0*I_ESP_I4x2y_Gx3y_a-1*I_ESP_G4x_Gx3y;
    abcd[iGrid*945+443] = 2.0E0*I_ESP_I4xyz_Gx3y_a;
    abcd[iGrid*945+444] = 2.0E0*I_ESP_I3x3y_Gx3y_a-2*I_ESP_G3xy_Gx3y;
    abcd[iGrid*945+445] = 2.0E0*I_ESP_I3x2yz_Gx3y_a-1*I_ESP_G3xz_Gx3y;
    abcd[iGrid*945+446] = 2.0E0*I_ESP_I3xy2z_Gx3y_a;
    abcd[iGrid*945+447] = 2.0E0*I_ESP_I2x4y_Gx3y_a-3*I_ESP_G2x2y_Gx3y;
    abcd[iGrid*945+448] = 2.0E0*I_ESP_I2x3yz_Gx3y_a-2*I_ESP_G2xyz_Gx3y;
    abcd[iGrid*945+449] = 2.0E0*I_ESP_I2x2y2z_Gx3y_a-1*I_ESP_G2x2z_Gx3y;
    abcd[iGrid*945+450] = 2.0E0*I_ESP_I2xy3z_Gx3y_a;
    abcd[iGrid*945+451] = 2.0E0*I_ESP_Ix5y_Gx3y_a-4*I_ESP_Gx3y_Gx3y;
    abcd[iGrid*945+452] = 2.0E0*I_ESP_Ix4yz_Gx3y_a-3*I_ESP_Gx2yz_Gx3y;
    abcd[iGrid*945+453] = 2.0E0*I_ESP_Ix3y2z_Gx3y_a-2*I_ESP_Gxy2z_Gx3y;
    abcd[iGrid*945+454] = 2.0E0*I_ESP_Ix2y3z_Gx3y_a-1*I_ESP_Gx3z_Gx3y;
    abcd[iGrid*945+455] = 2.0E0*I_ESP_Ixy4z_Gx3y_a;
    abcd[iGrid*945+456] = 2.0E0*I_ESP_I6y_Gx3y_a-5*I_ESP_G4y_Gx3y;
    abcd[iGrid*945+457] = 2.0E0*I_ESP_I5yz_Gx3y_a-4*I_ESP_G3yz_Gx3y;
    abcd[iGrid*945+458] = 2.0E0*I_ESP_I4y2z_Gx3y_a-3*I_ESP_G2y2z_Gx3y;
    abcd[iGrid*945+459] = 2.0E0*I_ESP_I3y3z_Gx3y_a-2*I_ESP_Gy3z_Gx3y;
    abcd[iGrid*945+460] = 2.0E0*I_ESP_I2y4z_Gx3y_a-1*I_ESP_G4z_Gx3y;
    abcd[iGrid*945+461] = 2.0E0*I_ESP_Iy5z_Gx3y_a;
    abcd[iGrid*945+462] = 2.0E0*I_ESP_I5xy_Gx2yz_a;
    abcd[iGrid*945+463] = 2.0E0*I_ESP_I4x2y_Gx2yz_a-1*I_ESP_G4x_Gx2yz;
    abcd[iGrid*945+464] = 2.0E0*I_ESP_I4xyz_Gx2yz_a;
    abcd[iGrid*945+465] = 2.0E0*I_ESP_I3x3y_Gx2yz_a-2*I_ESP_G3xy_Gx2yz;
    abcd[iGrid*945+466] = 2.0E0*I_ESP_I3x2yz_Gx2yz_a-1*I_ESP_G3xz_Gx2yz;
    abcd[iGrid*945+467] = 2.0E0*I_ESP_I3xy2z_Gx2yz_a;
    abcd[iGrid*945+468] = 2.0E0*I_ESP_I2x4y_Gx2yz_a-3*I_ESP_G2x2y_Gx2yz;
    abcd[iGrid*945+469] = 2.0E0*I_ESP_I2x3yz_Gx2yz_a-2*I_ESP_G2xyz_Gx2yz;
    abcd[iGrid*945+470] = 2.0E0*I_ESP_I2x2y2z_Gx2yz_a-1*I_ESP_G2x2z_Gx2yz;
    abcd[iGrid*945+471] = 2.0E0*I_ESP_I2xy3z_Gx2yz_a;
    abcd[iGrid*945+472] = 2.0E0*I_ESP_Ix5y_Gx2yz_a-4*I_ESP_Gx3y_Gx2yz;
    abcd[iGrid*945+473] = 2.0E0*I_ESP_Ix4yz_Gx2yz_a-3*I_ESP_Gx2yz_Gx2yz;
    abcd[iGrid*945+474] = 2.0E0*I_ESP_Ix3y2z_Gx2yz_a-2*I_ESP_Gxy2z_Gx2yz;
    abcd[iGrid*945+475] = 2.0E0*I_ESP_Ix2y3z_Gx2yz_a-1*I_ESP_Gx3z_Gx2yz;
    abcd[iGrid*945+476] = 2.0E0*I_ESP_Ixy4z_Gx2yz_a;
    abcd[iGrid*945+477] = 2.0E0*I_ESP_I6y_Gx2yz_a-5*I_ESP_G4y_Gx2yz;
    abcd[iGrid*945+478] = 2.0E0*I_ESP_I5yz_Gx2yz_a-4*I_ESP_G3yz_Gx2yz;
    abcd[iGrid*945+479] = 2.0E0*I_ESP_I4y2z_Gx2yz_a-3*I_ESP_G2y2z_Gx2yz;
    abcd[iGrid*945+480] = 2.0E0*I_ESP_I3y3z_Gx2yz_a-2*I_ESP_Gy3z_Gx2yz;
    abcd[iGrid*945+481] = 2.0E0*I_ESP_I2y4z_Gx2yz_a-1*I_ESP_G4z_Gx2yz;
    abcd[iGrid*945+482] = 2.0E0*I_ESP_Iy5z_Gx2yz_a;
    abcd[iGrid*945+483] = 2.0E0*I_ESP_I5xy_Gxy2z_a;
    abcd[iGrid*945+484] = 2.0E0*I_ESP_I4x2y_Gxy2z_a-1*I_ESP_G4x_Gxy2z;
    abcd[iGrid*945+485] = 2.0E0*I_ESP_I4xyz_Gxy2z_a;
    abcd[iGrid*945+486] = 2.0E0*I_ESP_I3x3y_Gxy2z_a-2*I_ESP_G3xy_Gxy2z;
    abcd[iGrid*945+487] = 2.0E0*I_ESP_I3x2yz_Gxy2z_a-1*I_ESP_G3xz_Gxy2z;
    abcd[iGrid*945+488] = 2.0E0*I_ESP_I3xy2z_Gxy2z_a;
    abcd[iGrid*945+489] = 2.0E0*I_ESP_I2x4y_Gxy2z_a-3*I_ESP_G2x2y_Gxy2z;
    abcd[iGrid*945+490] = 2.0E0*I_ESP_I2x3yz_Gxy2z_a-2*I_ESP_G2xyz_Gxy2z;
    abcd[iGrid*945+491] = 2.0E0*I_ESP_I2x2y2z_Gxy2z_a-1*I_ESP_G2x2z_Gxy2z;
    abcd[iGrid*945+492] = 2.0E0*I_ESP_I2xy3z_Gxy2z_a;
    abcd[iGrid*945+493] = 2.0E0*I_ESP_Ix5y_Gxy2z_a-4*I_ESP_Gx3y_Gxy2z;
    abcd[iGrid*945+494] = 2.0E0*I_ESP_Ix4yz_Gxy2z_a-3*I_ESP_Gx2yz_Gxy2z;
    abcd[iGrid*945+495] = 2.0E0*I_ESP_Ix3y2z_Gxy2z_a-2*I_ESP_Gxy2z_Gxy2z;
    abcd[iGrid*945+496] = 2.0E0*I_ESP_Ix2y3z_Gxy2z_a-1*I_ESP_Gx3z_Gxy2z;
    abcd[iGrid*945+497] = 2.0E0*I_ESP_Ixy4z_Gxy2z_a;
    abcd[iGrid*945+498] = 2.0E0*I_ESP_I6y_Gxy2z_a-5*I_ESP_G4y_Gxy2z;
    abcd[iGrid*945+499] = 2.0E0*I_ESP_I5yz_Gxy2z_a-4*I_ESP_G3yz_Gxy2z;
    abcd[iGrid*945+500] = 2.0E0*I_ESP_I4y2z_Gxy2z_a-3*I_ESP_G2y2z_Gxy2z;
    abcd[iGrid*945+501] = 2.0E0*I_ESP_I3y3z_Gxy2z_a-2*I_ESP_Gy3z_Gxy2z;
    abcd[iGrid*945+502] = 2.0E0*I_ESP_I2y4z_Gxy2z_a-1*I_ESP_G4z_Gxy2z;
    abcd[iGrid*945+503] = 2.0E0*I_ESP_Iy5z_Gxy2z_a;
    abcd[iGrid*945+504] = 2.0E0*I_ESP_I5xy_Gx3z_a;
    abcd[iGrid*945+505] = 2.0E0*I_ESP_I4x2y_Gx3z_a-1*I_ESP_G4x_Gx3z;
    abcd[iGrid*945+506] = 2.0E0*I_ESP_I4xyz_Gx3z_a;
    abcd[iGrid*945+507] = 2.0E0*I_ESP_I3x3y_Gx3z_a-2*I_ESP_G3xy_Gx3z;
    abcd[iGrid*945+508] = 2.0E0*I_ESP_I3x2yz_Gx3z_a-1*I_ESP_G3xz_Gx3z;
    abcd[iGrid*945+509] = 2.0E0*I_ESP_I3xy2z_Gx3z_a;
    abcd[iGrid*945+510] = 2.0E0*I_ESP_I2x4y_Gx3z_a-3*I_ESP_G2x2y_Gx3z;
    abcd[iGrid*945+511] = 2.0E0*I_ESP_I2x3yz_Gx3z_a-2*I_ESP_G2xyz_Gx3z;
    abcd[iGrid*945+512] = 2.0E0*I_ESP_I2x2y2z_Gx3z_a-1*I_ESP_G2x2z_Gx3z;
    abcd[iGrid*945+513] = 2.0E0*I_ESP_I2xy3z_Gx3z_a;
    abcd[iGrid*945+514] = 2.0E0*I_ESP_Ix5y_Gx3z_a-4*I_ESP_Gx3y_Gx3z;
    abcd[iGrid*945+515] = 2.0E0*I_ESP_Ix4yz_Gx3z_a-3*I_ESP_Gx2yz_Gx3z;
    abcd[iGrid*945+516] = 2.0E0*I_ESP_Ix3y2z_Gx3z_a-2*I_ESP_Gxy2z_Gx3z;
    abcd[iGrid*945+517] = 2.0E0*I_ESP_Ix2y3z_Gx3z_a-1*I_ESP_Gx3z_Gx3z;
    abcd[iGrid*945+518] = 2.0E0*I_ESP_Ixy4z_Gx3z_a;
    abcd[iGrid*945+519] = 2.0E0*I_ESP_I6y_Gx3z_a-5*I_ESP_G4y_Gx3z;
    abcd[iGrid*945+520] = 2.0E0*I_ESP_I5yz_Gx3z_a-4*I_ESP_G3yz_Gx3z;
    abcd[iGrid*945+521] = 2.0E0*I_ESP_I4y2z_Gx3z_a-3*I_ESP_G2y2z_Gx3z;
    abcd[iGrid*945+522] = 2.0E0*I_ESP_I3y3z_Gx3z_a-2*I_ESP_Gy3z_Gx3z;
    abcd[iGrid*945+523] = 2.0E0*I_ESP_I2y4z_Gx3z_a-1*I_ESP_G4z_Gx3z;
    abcd[iGrid*945+524] = 2.0E0*I_ESP_Iy5z_Gx3z_a;
    abcd[iGrid*945+525] = 2.0E0*I_ESP_I5xy_G4y_a;
    abcd[iGrid*945+526] = 2.0E0*I_ESP_I4x2y_G4y_a-1*I_ESP_G4x_G4y;
    abcd[iGrid*945+527] = 2.0E0*I_ESP_I4xyz_G4y_a;
    abcd[iGrid*945+528] = 2.0E0*I_ESP_I3x3y_G4y_a-2*I_ESP_G3xy_G4y;
    abcd[iGrid*945+529] = 2.0E0*I_ESP_I3x2yz_G4y_a-1*I_ESP_G3xz_G4y;
    abcd[iGrid*945+530] = 2.0E0*I_ESP_I3xy2z_G4y_a;
    abcd[iGrid*945+531] = 2.0E0*I_ESP_I2x4y_G4y_a-3*I_ESP_G2x2y_G4y;
    abcd[iGrid*945+532] = 2.0E0*I_ESP_I2x3yz_G4y_a-2*I_ESP_G2xyz_G4y;
    abcd[iGrid*945+533] = 2.0E0*I_ESP_I2x2y2z_G4y_a-1*I_ESP_G2x2z_G4y;
    abcd[iGrid*945+534] = 2.0E0*I_ESP_I2xy3z_G4y_a;
    abcd[iGrid*945+535] = 2.0E0*I_ESP_Ix5y_G4y_a-4*I_ESP_Gx3y_G4y;
    abcd[iGrid*945+536] = 2.0E0*I_ESP_Ix4yz_G4y_a-3*I_ESP_Gx2yz_G4y;
    abcd[iGrid*945+537] = 2.0E0*I_ESP_Ix3y2z_G4y_a-2*I_ESP_Gxy2z_G4y;
    abcd[iGrid*945+538] = 2.0E0*I_ESP_Ix2y3z_G4y_a-1*I_ESP_Gx3z_G4y;
    abcd[iGrid*945+539] = 2.0E0*I_ESP_Ixy4z_G4y_a;
    abcd[iGrid*945+540] = 2.0E0*I_ESP_I6y_G4y_a-5*I_ESP_G4y_G4y;
    abcd[iGrid*945+541] = 2.0E0*I_ESP_I5yz_G4y_a-4*I_ESP_G3yz_G4y;
    abcd[iGrid*945+542] = 2.0E0*I_ESP_I4y2z_G4y_a-3*I_ESP_G2y2z_G4y;
    abcd[iGrid*945+543] = 2.0E0*I_ESP_I3y3z_G4y_a-2*I_ESP_Gy3z_G4y;
    abcd[iGrid*945+544] = 2.0E0*I_ESP_I2y4z_G4y_a-1*I_ESP_G4z_G4y;
    abcd[iGrid*945+545] = 2.0E0*I_ESP_Iy5z_G4y_a;
    abcd[iGrid*945+546] = 2.0E0*I_ESP_I5xy_G3yz_a;
    abcd[iGrid*945+547] = 2.0E0*I_ESP_I4x2y_G3yz_a-1*I_ESP_G4x_G3yz;
    abcd[iGrid*945+548] = 2.0E0*I_ESP_I4xyz_G3yz_a;
    abcd[iGrid*945+549] = 2.0E0*I_ESP_I3x3y_G3yz_a-2*I_ESP_G3xy_G3yz;
    abcd[iGrid*945+550] = 2.0E0*I_ESP_I3x2yz_G3yz_a-1*I_ESP_G3xz_G3yz;
    abcd[iGrid*945+551] = 2.0E0*I_ESP_I3xy2z_G3yz_a;
    abcd[iGrid*945+552] = 2.0E0*I_ESP_I2x4y_G3yz_a-3*I_ESP_G2x2y_G3yz;
    abcd[iGrid*945+553] = 2.0E0*I_ESP_I2x3yz_G3yz_a-2*I_ESP_G2xyz_G3yz;
    abcd[iGrid*945+554] = 2.0E0*I_ESP_I2x2y2z_G3yz_a-1*I_ESP_G2x2z_G3yz;
    abcd[iGrid*945+555] = 2.0E0*I_ESP_I2xy3z_G3yz_a;
    abcd[iGrid*945+556] = 2.0E0*I_ESP_Ix5y_G3yz_a-4*I_ESP_Gx3y_G3yz;
    abcd[iGrid*945+557] = 2.0E0*I_ESP_Ix4yz_G3yz_a-3*I_ESP_Gx2yz_G3yz;
    abcd[iGrid*945+558] = 2.0E0*I_ESP_Ix3y2z_G3yz_a-2*I_ESP_Gxy2z_G3yz;
    abcd[iGrid*945+559] = 2.0E0*I_ESP_Ix2y3z_G3yz_a-1*I_ESP_Gx3z_G3yz;
    abcd[iGrid*945+560] = 2.0E0*I_ESP_Ixy4z_G3yz_a;
    abcd[iGrid*945+561] = 2.0E0*I_ESP_I6y_G3yz_a-5*I_ESP_G4y_G3yz;
    abcd[iGrid*945+562] = 2.0E0*I_ESP_I5yz_G3yz_a-4*I_ESP_G3yz_G3yz;
    abcd[iGrid*945+563] = 2.0E0*I_ESP_I4y2z_G3yz_a-3*I_ESP_G2y2z_G3yz;
    abcd[iGrid*945+564] = 2.0E0*I_ESP_I3y3z_G3yz_a-2*I_ESP_Gy3z_G3yz;
    abcd[iGrid*945+565] = 2.0E0*I_ESP_I2y4z_G3yz_a-1*I_ESP_G4z_G3yz;
    abcd[iGrid*945+566] = 2.0E0*I_ESP_Iy5z_G3yz_a;
    abcd[iGrid*945+567] = 2.0E0*I_ESP_I5xy_G2y2z_a;
    abcd[iGrid*945+568] = 2.0E0*I_ESP_I4x2y_G2y2z_a-1*I_ESP_G4x_G2y2z;
    abcd[iGrid*945+569] = 2.0E0*I_ESP_I4xyz_G2y2z_a;
    abcd[iGrid*945+570] = 2.0E0*I_ESP_I3x3y_G2y2z_a-2*I_ESP_G3xy_G2y2z;
    abcd[iGrid*945+571] = 2.0E0*I_ESP_I3x2yz_G2y2z_a-1*I_ESP_G3xz_G2y2z;
    abcd[iGrid*945+572] = 2.0E0*I_ESP_I3xy2z_G2y2z_a;
    abcd[iGrid*945+573] = 2.0E0*I_ESP_I2x4y_G2y2z_a-3*I_ESP_G2x2y_G2y2z;
    abcd[iGrid*945+574] = 2.0E0*I_ESP_I2x3yz_G2y2z_a-2*I_ESP_G2xyz_G2y2z;
    abcd[iGrid*945+575] = 2.0E0*I_ESP_I2x2y2z_G2y2z_a-1*I_ESP_G2x2z_G2y2z;
    abcd[iGrid*945+576] = 2.0E0*I_ESP_I2xy3z_G2y2z_a;
    abcd[iGrid*945+577] = 2.0E0*I_ESP_Ix5y_G2y2z_a-4*I_ESP_Gx3y_G2y2z;
    abcd[iGrid*945+578] = 2.0E0*I_ESP_Ix4yz_G2y2z_a-3*I_ESP_Gx2yz_G2y2z;
    abcd[iGrid*945+579] = 2.0E0*I_ESP_Ix3y2z_G2y2z_a-2*I_ESP_Gxy2z_G2y2z;
    abcd[iGrid*945+580] = 2.0E0*I_ESP_Ix2y3z_G2y2z_a-1*I_ESP_Gx3z_G2y2z;
    abcd[iGrid*945+581] = 2.0E0*I_ESP_Ixy4z_G2y2z_a;
    abcd[iGrid*945+582] = 2.0E0*I_ESP_I6y_G2y2z_a-5*I_ESP_G4y_G2y2z;
    abcd[iGrid*945+583] = 2.0E0*I_ESP_I5yz_G2y2z_a-4*I_ESP_G3yz_G2y2z;
    abcd[iGrid*945+584] = 2.0E0*I_ESP_I4y2z_G2y2z_a-3*I_ESP_G2y2z_G2y2z;
    abcd[iGrid*945+585] = 2.0E0*I_ESP_I3y3z_G2y2z_a-2*I_ESP_Gy3z_G2y2z;
    abcd[iGrid*945+586] = 2.0E0*I_ESP_I2y4z_G2y2z_a-1*I_ESP_G4z_G2y2z;
    abcd[iGrid*945+587] = 2.0E0*I_ESP_Iy5z_G2y2z_a;
    abcd[iGrid*945+588] = 2.0E0*I_ESP_I5xy_Gy3z_a;
    abcd[iGrid*945+589] = 2.0E0*I_ESP_I4x2y_Gy3z_a-1*I_ESP_G4x_Gy3z;
    abcd[iGrid*945+590] = 2.0E0*I_ESP_I4xyz_Gy3z_a;
    abcd[iGrid*945+591] = 2.0E0*I_ESP_I3x3y_Gy3z_a-2*I_ESP_G3xy_Gy3z;
    abcd[iGrid*945+592] = 2.0E0*I_ESP_I3x2yz_Gy3z_a-1*I_ESP_G3xz_Gy3z;
    abcd[iGrid*945+593] = 2.0E0*I_ESP_I3xy2z_Gy3z_a;
    abcd[iGrid*945+594] = 2.0E0*I_ESP_I2x4y_Gy3z_a-3*I_ESP_G2x2y_Gy3z;
    abcd[iGrid*945+595] = 2.0E0*I_ESP_I2x3yz_Gy3z_a-2*I_ESP_G2xyz_Gy3z;
    abcd[iGrid*945+596] = 2.0E0*I_ESP_I2x2y2z_Gy3z_a-1*I_ESP_G2x2z_Gy3z;
    abcd[iGrid*945+597] = 2.0E0*I_ESP_I2xy3z_Gy3z_a;
    abcd[iGrid*945+598] = 2.0E0*I_ESP_Ix5y_Gy3z_a-4*I_ESP_Gx3y_Gy3z;
    abcd[iGrid*945+599] = 2.0E0*I_ESP_Ix4yz_Gy3z_a-3*I_ESP_Gx2yz_Gy3z;
    abcd[iGrid*945+600] = 2.0E0*I_ESP_Ix3y2z_Gy3z_a-2*I_ESP_Gxy2z_Gy3z;
    abcd[iGrid*945+601] = 2.0E0*I_ESP_Ix2y3z_Gy3z_a-1*I_ESP_Gx3z_Gy3z;
    abcd[iGrid*945+602] = 2.0E0*I_ESP_Ixy4z_Gy3z_a;
    abcd[iGrid*945+603] = 2.0E0*I_ESP_I6y_Gy3z_a-5*I_ESP_G4y_Gy3z;
    abcd[iGrid*945+604] = 2.0E0*I_ESP_I5yz_Gy3z_a-4*I_ESP_G3yz_Gy3z;
    abcd[iGrid*945+605] = 2.0E0*I_ESP_I4y2z_Gy3z_a-3*I_ESP_G2y2z_Gy3z;
    abcd[iGrid*945+606] = 2.0E0*I_ESP_I3y3z_Gy3z_a-2*I_ESP_Gy3z_Gy3z;
    abcd[iGrid*945+607] = 2.0E0*I_ESP_I2y4z_Gy3z_a-1*I_ESP_G4z_Gy3z;
    abcd[iGrid*945+608] = 2.0E0*I_ESP_Iy5z_Gy3z_a;
    abcd[iGrid*945+609] = 2.0E0*I_ESP_I5xy_G4z_a;
    abcd[iGrid*945+610] = 2.0E0*I_ESP_I4x2y_G4z_a-1*I_ESP_G4x_G4z;
    abcd[iGrid*945+611] = 2.0E0*I_ESP_I4xyz_G4z_a;
    abcd[iGrid*945+612] = 2.0E0*I_ESP_I3x3y_G4z_a-2*I_ESP_G3xy_G4z;
    abcd[iGrid*945+613] = 2.0E0*I_ESP_I3x2yz_G4z_a-1*I_ESP_G3xz_G4z;
    abcd[iGrid*945+614] = 2.0E0*I_ESP_I3xy2z_G4z_a;
    abcd[iGrid*945+615] = 2.0E0*I_ESP_I2x4y_G4z_a-3*I_ESP_G2x2y_G4z;
    abcd[iGrid*945+616] = 2.0E0*I_ESP_I2x3yz_G4z_a-2*I_ESP_G2xyz_G4z;
    abcd[iGrid*945+617] = 2.0E0*I_ESP_I2x2y2z_G4z_a-1*I_ESP_G2x2z_G4z;
    abcd[iGrid*945+618] = 2.0E0*I_ESP_I2xy3z_G4z_a;
    abcd[iGrid*945+619] = 2.0E0*I_ESP_Ix5y_G4z_a-4*I_ESP_Gx3y_G4z;
    abcd[iGrid*945+620] = 2.0E0*I_ESP_Ix4yz_G4z_a-3*I_ESP_Gx2yz_G4z;
    abcd[iGrid*945+621] = 2.0E0*I_ESP_Ix3y2z_G4z_a-2*I_ESP_Gxy2z_G4z;
    abcd[iGrid*945+622] = 2.0E0*I_ESP_Ix2y3z_G4z_a-1*I_ESP_Gx3z_G4z;
    abcd[iGrid*945+623] = 2.0E0*I_ESP_Ixy4z_G4z_a;
    abcd[iGrid*945+624] = 2.0E0*I_ESP_I6y_G4z_a-5*I_ESP_G4y_G4z;
    abcd[iGrid*945+625] = 2.0E0*I_ESP_I5yz_G4z_a-4*I_ESP_G3yz_G4z;
    abcd[iGrid*945+626] = 2.0E0*I_ESP_I4y2z_G4z_a-3*I_ESP_G2y2z_G4z;
    abcd[iGrid*945+627] = 2.0E0*I_ESP_I3y3z_G4z_a-2*I_ESP_Gy3z_G4z;
    abcd[iGrid*945+628] = 2.0E0*I_ESP_I2y4z_G4z_a-1*I_ESP_G4z_G4z;
    abcd[iGrid*945+629] = 2.0E0*I_ESP_Iy5z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G_a
     * RHS shell quartet name: SQ_ESP_G_G
     ************************************************************/
    abcd[iGrid*945+630] = 2.0E0*I_ESP_I5xz_G4x_a;
    abcd[iGrid*945+631] = 2.0E0*I_ESP_I4xyz_G4x_a;
    abcd[iGrid*945+632] = 2.0E0*I_ESP_I4x2z_G4x_a-1*I_ESP_G4x_G4x;
    abcd[iGrid*945+633] = 2.0E0*I_ESP_I3x2yz_G4x_a;
    abcd[iGrid*945+634] = 2.0E0*I_ESP_I3xy2z_G4x_a-1*I_ESP_G3xy_G4x;
    abcd[iGrid*945+635] = 2.0E0*I_ESP_I3x3z_G4x_a-2*I_ESP_G3xz_G4x;
    abcd[iGrid*945+636] = 2.0E0*I_ESP_I2x3yz_G4x_a;
    abcd[iGrid*945+637] = 2.0E0*I_ESP_I2x2y2z_G4x_a-1*I_ESP_G2x2y_G4x;
    abcd[iGrid*945+638] = 2.0E0*I_ESP_I2xy3z_G4x_a-2*I_ESP_G2xyz_G4x;
    abcd[iGrid*945+639] = 2.0E0*I_ESP_I2x4z_G4x_a-3*I_ESP_G2x2z_G4x;
    abcd[iGrid*945+640] = 2.0E0*I_ESP_Ix4yz_G4x_a;
    abcd[iGrid*945+641] = 2.0E0*I_ESP_Ix3y2z_G4x_a-1*I_ESP_Gx3y_G4x;
    abcd[iGrid*945+642] = 2.0E0*I_ESP_Ix2y3z_G4x_a-2*I_ESP_Gx2yz_G4x;
    abcd[iGrid*945+643] = 2.0E0*I_ESP_Ixy4z_G4x_a-3*I_ESP_Gxy2z_G4x;
    abcd[iGrid*945+644] = 2.0E0*I_ESP_Ix5z_G4x_a-4*I_ESP_Gx3z_G4x;
    abcd[iGrid*945+645] = 2.0E0*I_ESP_I5yz_G4x_a;
    abcd[iGrid*945+646] = 2.0E0*I_ESP_I4y2z_G4x_a-1*I_ESP_G4y_G4x;
    abcd[iGrid*945+647] = 2.0E0*I_ESP_I3y3z_G4x_a-2*I_ESP_G3yz_G4x;
    abcd[iGrid*945+648] = 2.0E0*I_ESP_I2y4z_G4x_a-3*I_ESP_G2y2z_G4x;
    abcd[iGrid*945+649] = 2.0E0*I_ESP_Iy5z_G4x_a-4*I_ESP_Gy3z_G4x;
    abcd[iGrid*945+650] = 2.0E0*I_ESP_I6z_G4x_a-5*I_ESP_G4z_G4x;
    abcd[iGrid*945+651] = 2.0E0*I_ESP_I5xz_G3xy_a;
    abcd[iGrid*945+652] = 2.0E0*I_ESP_I4xyz_G3xy_a;
    abcd[iGrid*945+653] = 2.0E0*I_ESP_I4x2z_G3xy_a-1*I_ESP_G4x_G3xy;
    abcd[iGrid*945+654] = 2.0E0*I_ESP_I3x2yz_G3xy_a;
    abcd[iGrid*945+655] = 2.0E0*I_ESP_I3xy2z_G3xy_a-1*I_ESP_G3xy_G3xy;
    abcd[iGrid*945+656] = 2.0E0*I_ESP_I3x3z_G3xy_a-2*I_ESP_G3xz_G3xy;
    abcd[iGrid*945+657] = 2.0E0*I_ESP_I2x3yz_G3xy_a;
    abcd[iGrid*945+658] = 2.0E0*I_ESP_I2x2y2z_G3xy_a-1*I_ESP_G2x2y_G3xy;
    abcd[iGrid*945+659] = 2.0E0*I_ESP_I2xy3z_G3xy_a-2*I_ESP_G2xyz_G3xy;
    abcd[iGrid*945+660] = 2.0E0*I_ESP_I2x4z_G3xy_a-3*I_ESP_G2x2z_G3xy;
    abcd[iGrid*945+661] = 2.0E0*I_ESP_Ix4yz_G3xy_a;
    abcd[iGrid*945+662] = 2.0E0*I_ESP_Ix3y2z_G3xy_a-1*I_ESP_Gx3y_G3xy;
    abcd[iGrid*945+663] = 2.0E0*I_ESP_Ix2y3z_G3xy_a-2*I_ESP_Gx2yz_G3xy;
    abcd[iGrid*945+664] = 2.0E0*I_ESP_Ixy4z_G3xy_a-3*I_ESP_Gxy2z_G3xy;
    abcd[iGrid*945+665] = 2.0E0*I_ESP_Ix5z_G3xy_a-4*I_ESP_Gx3z_G3xy;
    abcd[iGrid*945+666] = 2.0E0*I_ESP_I5yz_G3xy_a;
    abcd[iGrid*945+667] = 2.0E0*I_ESP_I4y2z_G3xy_a-1*I_ESP_G4y_G3xy;
    abcd[iGrid*945+668] = 2.0E0*I_ESP_I3y3z_G3xy_a-2*I_ESP_G3yz_G3xy;
    abcd[iGrid*945+669] = 2.0E0*I_ESP_I2y4z_G3xy_a-3*I_ESP_G2y2z_G3xy;
    abcd[iGrid*945+670] = 2.0E0*I_ESP_Iy5z_G3xy_a-4*I_ESP_Gy3z_G3xy;
    abcd[iGrid*945+671] = 2.0E0*I_ESP_I6z_G3xy_a-5*I_ESP_G4z_G3xy;
    abcd[iGrid*945+672] = 2.0E0*I_ESP_I5xz_G3xz_a;
    abcd[iGrid*945+673] = 2.0E0*I_ESP_I4xyz_G3xz_a;
    abcd[iGrid*945+674] = 2.0E0*I_ESP_I4x2z_G3xz_a-1*I_ESP_G4x_G3xz;
    abcd[iGrid*945+675] = 2.0E0*I_ESP_I3x2yz_G3xz_a;
    abcd[iGrid*945+676] = 2.0E0*I_ESP_I3xy2z_G3xz_a-1*I_ESP_G3xy_G3xz;
    abcd[iGrid*945+677] = 2.0E0*I_ESP_I3x3z_G3xz_a-2*I_ESP_G3xz_G3xz;
    abcd[iGrid*945+678] = 2.0E0*I_ESP_I2x3yz_G3xz_a;
    abcd[iGrid*945+679] = 2.0E0*I_ESP_I2x2y2z_G3xz_a-1*I_ESP_G2x2y_G3xz;
    abcd[iGrid*945+680] = 2.0E0*I_ESP_I2xy3z_G3xz_a-2*I_ESP_G2xyz_G3xz;
    abcd[iGrid*945+681] = 2.0E0*I_ESP_I2x4z_G3xz_a-3*I_ESP_G2x2z_G3xz;
    abcd[iGrid*945+682] = 2.0E0*I_ESP_Ix4yz_G3xz_a;
    abcd[iGrid*945+683] = 2.0E0*I_ESP_Ix3y2z_G3xz_a-1*I_ESP_Gx3y_G3xz;
    abcd[iGrid*945+684] = 2.0E0*I_ESP_Ix2y3z_G3xz_a-2*I_ESP_Gx2yz_G3xz;
    abcd[iGrid*945+685] = 2.0E0*I_ESP_Ixy4z_G3xz_a-3*I_ESP_Gxy2z_G3xz;
    abcd[iGrid*945+686] = 2.0E0*I_ESP_Ix5z_G3xz_a-4*I_ESP_Gx3z_G3xz;
    abcd[iGrid*945+687] = 2.0E0*I_ESP_I5yz_G3xz_a;
    abcd[iGrid*945+688] = 2.0E0*I_ESP_I4y2z_G3xz_a-1*I_ESP_G4y_G3xz;
    abcd[iGrid*945+689] = 2.0E0*I_ESP_I3y3z_G3xz_a-2*I_ESP_G3yz_G3xz;
    abcd[iGrid*945+690] = 2.0E0*I_ESP_I2y4z_G3xz_a-3*I_ESP_G2y2z_G3xz;
    abcd[iGrid*945+691] = 2.0E0*I_ESP_Iy5z_G3xz_a-4*I_ESP_Gy3z_G3xz;
    abcd[iGrid*945+692] = 2.0E0*I_ESP_I6z_G3xz_a-5*I_ESP_G4z_G3xz;
    abcd[iGrid*945+693] = 2.0E0*I_ESP_I5xz_G2x2y_a;
    abcd[iGrid*945+694] = 2.0E0*I_ESP_I4xyz_G2x2y_a;
    abcd[iGrid*945+695] = 2.0E0*I_ESP_I4x2z_G2x2y_a-1*I_ESP_G4x_G2x2y;
    abcd[iGrid*945+696] = 2.0E0*I_ESP_I3x2yz_G2x2y_a;
    abcd[iGrid*945+697] = 2.0E0*I_ESP_I3xy2z_G2x2y_a-1*I_ESP_G3xy_G2x2y;
    abcd[iGrid*945+698] = 2.0E0*I_ESP_I3x3z_G2x2y_a-2*I_ESP_G3xz_G2x2y;
    abcd[iGrid*945+699] = 2.0E0*I_ESP_I2x3yz_G2x2y_a;
    abcd[iGrid*945+700] = 2.0E0*I_ESP_I2x2y2z_G2x2y_a-1*I_ESP_G2x2y_G2x2y;
    abcd[iGrid*945+701] = 2.0E0*I_ESP_I2xy3z_G2x2y_a-2*I_ESP_G2xyz_G2x2y;
    abcd[iGrid*945+702] = 2.0E0*I_ESP_I2x4z_G2x2y_a-3*I_ESP_G2x2z_G2x2y;
    abcd[iGrid*945+703] = 2.0E0*I_ESP_Ix4yz_G2x2y_a;
    abcd[iGrid*945+704] = 2.0E0*I_ESP_Ix3y2z_G2x2y_a-1*I_ESP_Gx3y_G2x2y;
    abcd[iGrid*945+705] = 2.0E0*I_ESP_Ix2y3z_G2x2y_a-2*I_ESP_Gx2yz_G2x2y;
    abcd[iGrid*945+706] = 2.0E0*I_ESP_Ixy4z_G2x2y_a-3*I_ESP_Gxy2z_G2x2y;
    abcd[iGrid*945+707] = 2.0E0*I_ESP_Ix5z_G2x2y_a-4*I_ESP_Gx3z_G2x2y;
    abcd[iGrid*945+708] = 2.0E0*I_ESP_I5yz_G2x2y_a;
    abcd[iGrid*945+709] = 2.0E0*I_ESP_I4y2z_G2x2y_a-1*I_ESP_G4y_G2x2y;
    abcd[iGrid*945+710] = 2.0E0*I_ESP_I3y3z_G2x2y_a-2*I_ESP_G3yz_G2x2y;
    abcd[iGrid*945+711] = 2.0E0*I_ESP_I2y4z_G2x2y_a-3*I_ESP_G2y2z_G2x2y;
    abcd[iGrid*945+712] = 2.0E0*I_ESP_Iy5z_G2x2y_a-4*I_ESP_Gy3z_G2x2y;
    abcd[iGrid*945+713] = 2.0E0*I_ESP_I6z_G2x2y_a-5*I_ESP_G4z_G2x2y;
    abcd[iGrid*945+714] = 2.0E0*I_ESP_I5xz_G2xyz_a;
    abcd[iGrid*945+715] = 2.0E0*I_ESP_I4xyz_G2xyz_a;
    abcd[iGrid*945+716] = 2.0E0*I_ESP_I4x2z_G2xyz_a-1*I_ESP_G4x_G2xyz;
    abcd[iGrid*945+717] = 2.0E0*I_ESP_I3x2yz_G2xyz_a;
    abcd[iGrid*945+718] = 2.0E0*I_ESP_I3xy2z_G2xyz_a-1*I_ESP_G3xy_G2xyz;
    abcd[iGrid*945+719] = 2.0E0*I_ESP_I3x3z_G2xyz_a-2*I_ESP_G3xz_G2xyz;
    abcd[iGrid*945+720] = 2.0E0*I_ESP_I2x3yz_G2xyz_a;
    abcd[iGrid*945+721] = 2.0E0*I_ESP_I2x2y2z_G2xyz_a-1*I_ESP_G2x2y_G2xyz;
    abcd[iGrid*945+722] = 2.0E0*I_ESP_I2xy3z_G2xyz_a-2*I_ESP_G2xyz_G2xyz;
    abcd[iGrid*945+723] = 2.0E0*I_ESP_I2x4z_G2xyz_a-3*I_ESP_G2x2z_G2xyz;
    abcd[iGrid*945+724] = 2.0E0*I_ESP_Ix4yz_G2xyz_a;
    abcd[iGrid*945+725] = 2.0E0*I_ESP_Ix3y2z_G2xyz_a-1*I_ESP_Gx3y_G2xyz;
    abcd[iGrid*945+726] = 2.0E0*I_ESP_Ix2y3z_G2xyz_a-2*I_ESP_Gx2yz_G2xyz;
    abcd[iGrid*945+727] = 2.0E0*I_ESP_Ixy4z_G2xyz_a-3*I_ESP_Gxy2z_G2xyz;
    abcd[iGrid*945+728] = 2.0E0*I_ESP_Ix5z_G2xyz_a-4*I_ESP_Gx3z_G2xyz;
    abcd[iGrid*945+729] = 2.0E0*I_ESP_I5yz_G2xyz_a;
    abcd[iGrid*945+730] = 2.0E0*I_ESP_I4y2z_G2xyz_a-1*I_ESP_G4y_G2xyz;
    abcd[iGrid*945+731] = 2.0E0*I_ESP_I3y3z_G2xyz_a-2*I_ESP_G3yz_G2xyz;
    abcd[iGrid*945+732] = 2.0E0*I_ESP_I2y4z_G2xyz_a-3*I_ESP_G2y2z_G2xyz;
    abcd[iGrid*945+733] = 2.0E0*I_ESP_Iy5z_G2xyz_a-4*I_ESP_Gy3z_G2xyz;
    abcd[iGrid*945+734] = 2.0E0*I_ESP_I6z_G2xyz_a-5*I_ESP_G4z_G2xyz;
    abcd[iGrid*945+735] = 2.0E0*I_ESP_I5xz_G2x2z_a;
    abcd[iGrid*945+736] = 2.0E0*I_ESP_I4xyz_G2x2z_a;
    abcd[iGrid*945+737] = 2.0E0*I_ESP_I4x2z_G2x2z_a-1*I_ESP_G4x_G2x2z;
    abcd[iGrid*945+738] = 2.0E0*I_ESP_I3x2yz_G2x2z_a;
    abcd[iGrid*945+739] = 2.0E0*I_ESP_I3xy2z_G2x2z_a-1*I_ESP_G3xy_G2x2z;
    abcd[iGrid*945+740] = 2.0E0*I_ESP_I3x3z_G2x2z_a-2*I_ESP_G3xz_G2x2z;
    abcd[iGrid*945+741] = 2.0E0*I_ESP_I2x3yz_G2x2z_a;
    abcd[iGrid*945+742] = 2.0E0*I_ESP_I2x2y2z_G2x2z_a-1*I_ESP_G2x2y_G2x2z;
    abcd[iGrid*945+743] = 2.0E0*I_ESP_I2xy3z_G2x2z_a-2*I_ESP_G2xyz_G2x2z;
    abcd[iGrid*945+744] = 2.0E0*I_ESP_I2x4z_G2x2z_a-3*I_ESP_G2x2z_G2x2z;
    abcd[iGrid*945+745] = 2.0E0*I_ESP_Ix4yz_G2x2z_a;
    abcd[iGrid*945+746] = 2.0E0*I_ESP_Ix3y2z_G2x2z_a-1*I_ESP_Gx3y_G2x2z;
    abcd[iGrid*945+747] = 2.0E0*I_ESP_Ix2y3z_G2x2z_a-2*I_ESP_Gx2yz_G2x2z;
    abcd[iGrid*945+748] = 2.0E0*I_ESP_Ixy4z_G2x2z_a-3*I_ESP_Gxy2z_G2x2z;
    abcd[iGrid*945+749] = 2.0E0*I_ESP_Ix5z_G2x2z_a-4*I_ESP_Gx3z_G2x2z;
    abcd[iGrid*945+750] = 2.0E0*I_ESP_I5yz_G2x2z_a;
    abcd[iGrid*945+751] = 2.0E0*I_ESP_I4y2z_G2x2z_a-1*I_ESP_G4y_G2x2z;
    abcd[iGrid*945+752] = 2.0E0*I_ESP_I3y3z_G2x2z_a-2*I_ESP_G3yz_G2x2z;
    abcd[iGrid*945+753] = 2.0E0*I_ESP_I2y4z_G2x2z_a-3*I_ESP_G2y2z_G2x2z;
    abcd[iGrid*945+754] = 2.0E0*I_ESP_Iy5z_G2x2z_a-4*I_ESP_Gy3z_G2x2z;
    abcd[iGrid*945+755] = 2.0E0*I_ESP_I6z_G2x2z_a-5*I_ESP_G4z_G2x2z;
    abcd[iGrid*945+756] = 2.0E0*I_ESP_I5xz_Gx3y_a;
    abcd[iGrid*945+757] = 2.0E0*I_ESP_I4xyz_Gx3y_a;
    abcd[iGrid*945+758] = 2.0E0*I_ESP_I4x2z_Gx3y_a-1*I_ESP_G4x_Gx3y;
    abcd[iGrid*945+759] = 2.0E0*I_ESP_I3x2yz_Gx3y_a;
    abcd[iGrid*945+760] = 2.0E0*I_ESP_I3xy2z_Gx3y_a-1*I_ESP_G3xy_Gx3y;
    abcd[iGrid*945+761] = 2.0E0*I_ESP_I3x3z_Gx3y_a-2*I_ESP_G3xz_Gx3y;
    abcd[iGrid*945+762] = 2.0E0*I_ESP_I2x3yz_Gx3y_a;
    abcd[iGrid*945+763] = 2.0E0*I_ESP_I2x2y2z_Gx3y_a-1*I_ESP_G2x2y_Gx3y;
    abcd[iGrid*945+764] = 2.0E0*I_ESP_I2xy3z_Gx3y_a-2*I_ESP_G2xyz_Gx3y;
    abcd[iGrid*945+765] = 2.0E0*I_ESP_I2x4z_Gx3y_a-3*I_ESP_G2x2z_Gx3y;
    abcd[iGrid*945+766] = 2.0E0*I_ESP_Ix4yz_Gx3y_a;
    abcd[iGrid*945+767] = 2.0E0*I_ESP_Ix3y2z_Gx3y_a-1*I_ESP_Gx3y_Gx3y;
    abcd[iGrid*945+768] = 2.0E0*I_ESP_Ix2y3z_Gx3y_a-2*I_ESP_Gx2yz_Gx3y;
    abcd[iGrid*945+769] = 2.0E0*I_ESP_Ixy4z_Gx3y_a-3*I_ESP_Gxy2z_Gx3y;
    abcd[iGrid*945+770] = 2.0E0*I_ESP_Ix5z_Gx3y_a-4*I_ESP_Gx3z_Gx3y;
    abcd[iGrid*945+771] = 2.0E0*I_ESP_I5yz_Gx3y_a;
    abcd[iGrid*945+772] = 2.0E0*I_ESP_I4y2z_Gx3y_a-1*I_ESP_G4y_Gx3y;
    abcd[iGrid*945+773] = 2.0E0*I_ESP_I3y3z_Gx3y_a-2*I_ESP_G3yz_Gx3y;
    abcd[iGrid*945+774] = 2.0E0*I_ESP_I2y4z_Gx3y_a-3*I_ESP_G2y2z_Gx3y;
    abcd[iGrid*945+775] = 2.0E0*I_ESP_Iy5z_Gx3y_a-4*I_ESP_Gy3z_Gx3y;
    abcd[iGrid*945+776] = 2.0E0*I_ESP_I6z_Gx3y_a-5*I_ESP_G4z_Gx3y;
    abcd[iGrid*945+777] = 2.0E0*I_ESP_I5xz_Gx2yz_a;
    abcd[iGrid*945+778] = 2.0E0*I_ESP_I4xyz_Gx2yz_a;
    abcd[iGrid*945+779] = 2.0E0*I_ESP_I4x2z_Gx2yz_a-1*I_ESP_G4x_Gx2yz;
    abcd[iGrid*945+780] = 2.0E0*I_ESP_I3x2yz_Gx2yz_a;
    abcd[iGrid*945+781] = 2.0E0*I_ESP_I3xy2z_Gx2yz_a-1*I_ESP_G3xy_Gx2yz;
    abcd[iGrid*945+782] = 2.0E0*I_ESP_I3x3z_Gx2yz_a-2*I_ESP_G3xz_Gx2yz;
    abcd[iGrid*945+783] = 2.0E0*I_ESP_I2x3yz_Gx2yz_a;
    abcd[iGrid*945+784] = 2.0E0*I_ESP_I2x2y2z_Gx2yz_a-1*I_ESP_G2x2y_Gx2yz;
    abcd[iGrid*945+785] = 2.0E0*I_ESP_I2xy3z_Gx2yz_a-2*I_ESP_G2xyz_Gx2yz;
    abcd[iGrid*945+786] = 2.0E0*I_ESP_I2x4z_Gx2yz_a-3*I_ESP_G2x2z_Gx2yz;
    abcd[iGrid*945+787] = 2.0E0*I_ESP_Ix4yz_Gx2yz_a;
    abcd[iGrid*945+788] = 2.0E0*I_ESP_Ix3y2z_Gx2yz_a-1*I_ESP_Gx3y_Gx2yz;
    abcd[iGrid*945+789] = 2.0E0*I_ESP_Ix2y3z_Gx2yz_a-2*I_ESP_Gx2yz_Gx2yz;
    abcd[iGrid*945+790] = 2.0E0*I_ESP_Ixy4z_Gx2yz_a-3*I_ESP_Gxy2z_Gx2yz;
    abcd[iGrid*945+791] = 2.0E0*I_ESP_Ix5z_Gx2yz_a-4*I_ESP_Gx3z_Gx2yz;
    abcd[iGrid*945+792] = 2.0E0*I_ESP_I5yz_Gx2yz_a;
    abcd[iGrid*945+793] = 2.0E0*I_ESP_I4y2z_Gx2yz_a-1*I_ESP_G4y_Gx2yz;
    abcd[iGrid*945+794] = 2.0E0*I_ESP_I3y3z_Gx2yz_a-2*I_ESP_G3yz_Gx2yz;
    abcd[iGrid*945+795] = 2.0E0*I_ESP_I2y4z_Gx2yz_a-3*I_ESP_G2y2z_Gx2yz;
    abcd[iGrid*945+796] = 2.0E0*I_ESP_Iy5z_Gx2yz_a-4*I_ESP_Gy3z_Gx2yz;
    abcd[iGrid*945+797] = 2.0E0*I_ESP_I6z_Gx2yz_a-5*I_ESP_G4z_Gx2yz;
    abcd[iGrid*945+798] = 2.0E0*I_ESP_I5xz_Gxy2z_a;
    abcd[iGrid*945+799] = 2.0E0*I_ESP_I4xyz_Gxy2z_a;
    abcd[iGrid*945+800] = 2.0E0*I_ESP_I4x2z_Gxy2z_a-1*I_ESP_G4x_Gxy2z;
    abcd[iGrid*945+801] = 2.0E0*I_ESP_I3x2yz_Gxy2z_a;
    abcd[iGrid*945+802] = 2.0E0*I_ESP_I3xy2z_Gxy2z_a-1*I_ESP_G3xy_Gxy2z;
    abcd[iGrid*945+803] = 2.0E0*I_ESP_I3x3z_Gxy2z_a-2*I_ESP_G3xz_Gxy2z;
    abcd[iGrid*945+804] = 2.0E0*I_ESP_I2x3yz_Gxy2z_a;
    abcd[iGrid*945+805] = 2.0E0*I_ESP_I2x2y2z_Gxy2z_a-1*I_ESP_G2x2y_Gxy2z;
    abcd[iGrid*945+806] = 2.0E0*I_ESP_I2xy3z_Gxy2z_a-2*I_ESP_G2xyz_Gxy2z;
    abcd[iGrid*945+807] = 2.0E0*I_ESP_I2x4z_Gxy2z_a-3*I_ESP_G2x2z_Gxy2z;
    abcd[iGrid*945+808] = 2.0E0*I_ESP_Ix4yz_Gxy2z_a;
    abcd[iGrid*945+809] = 2.0E0*I_ESP_Ix3y2z_Gxy2z_a-1*I_ESP_Gx3y_Gxy2z;
    abcd[iGrid*945+810] = 2.0E0*I_ESP_Ix2y3z_Gxy2z_a-2*I_ESP_Gx2yz_Gxy2z;
    abcd[iGrid*945+811] = 2.0E0*I_ESP_Ixy4z_Gxy2z_a-3*I_ESP_Gxy2z_Gxy2z;
    abcd[iGrid*945+812] = 2.0E0*I_ESP_Ix5z_Gxy2z_a-4*I_ESP_Gx3z_Gxy2z;
    abcd[iGrid*945+813] = 2.0E0*I_ESP_I5yz_Gxy2z_a;
    abcd[iGrid*945+814] = 2.0E0*I_ESP_I4y2z_Gxy2z_a-1*I_ESP_G4y_Gxy2z;
    abcd[iGrid*945+815] = 2.0E0*I_ESP_I3y3z_Gxy2z_a-2*I_ESP_G3yz_Gxy2z;
    abcd[iGrid*945+816] = 2.0E0*I_ESP_I2y4z_Gxy2z_a-3*I_ESP_G2y2z_Gxy2z;
    abcd[iGrid*945+817] = 2.0E0*I_ESP_Iy5z_Gxy2z_a-4*I_ESP_Gy3z_Gxy2z;
    abcd[iGrid*945+818] = 2.0E0*I_ESP_I6z_Gxy2z_a-5*I_ESP_G4z_Gxy2z;
    abcd[iGrid*945+819] = 2.0E0*I_ESP_I5xz_Gx3z_a;
    abcd[iGrid*945+820] = 2.0E0*I_ESP_I4xyz_Gx3z_a;
    abcd[iGrid*945+821] = 2.0E0*I_ESP_I4x2z_Gx3z_a-1*I_ESP_G4x_Gx3z;
    abcd[iGrid*945+822] = 2.0E0*I_ESP_I3x2yz_Gx3z_a;
    abcd[iGrid*945+823] = 2.0E0*I_ESP_I3xy2z_Gx3z_a-1*I_ESP_G3xy_Gx3z;
    abcd[iGrid*945+824] = 2.0E0*I_ESP_I3x3z_Gx3z_a-2*I_ESP_G3xz_Gx3z;
    abcd[iGrid*945+825] = 2.0E0*I_ESP_I2x3yz_Gx3z_a;
    abcd[iGrid*945+826] = 2.0E0*I_ESP_I2x2y2z_Gx3z_a-1*I_ESP_G2x2y_Gx3z;
    abcd[iGrid*945+827] = 2.0E0*I_ESP_I2xy3z_Gx3z_a-2*I_ESP_G2xyz_Gx3z;
    abcd[iGrid*945+828] = 2.0E0*I_ESP_I2x4z_Gx3z_a-3*I_ESP_G2x2z_Gx3z;
    abcd[iGrid*945+829] = 2.0E0*I_ESP_Ix4yz_Gx3z_a;
    abcd[iGrid*945+830] = 2.0E0*I_ESP_Ix3y2z_Gx3z_a-1*I_ESP_Gx3y_Gx3z;
    abcd[iGrid*945+831] = 2.0E0*I_ESP_Ix2y3z_Gx3z_a-2*I_ESP_Gx2yz_Gx3z;
    abcd[iGrid*945+832] = 2.0E0*I_ESP_Ixy4z_Gx3z_a-3*I_ESP_Gxy2z_Gx3z;
    abcd[iGrid*945+833] = 2.0E0*I_ESP_Ix5z_Gx3z_a-4*I_ESP_Gx3z_Gx3z;
    abcd[iGrid*945+834] = 2.0E0*I_ESP_I5yz_Gx3z_a;
    abcd[iGrid*945+835] = 2.0E0*I_ESP_I4y2z_Gx3z_a-1*I_ESP_G4y_Gx3z;
    abcd[iGrid*945+836] = 2.0E0*I_ESP_I3y3z_Gx3z_a-2*I_ESP_G3yz_Gx3z;
    abcd[iGrid*945+837] = 2.0E0*I_ESP_I2y4z_Gx3z_a-3*I_ESP_G2y2z_Gx3z;
    abcd[iGrid*945+838] = 2.0E0*I_ESP_Iy5z_Gx3z_a-4*I_ESP_Gy3z_Gx3z;
    abcd[iGrid*945+839] = 2.0E0*I_ESP_I6z_Gx3z_a-5*I_ESP_G4z_Gx3z;
    abcd[iGrid*945+840] = 2.0E0*I_ESP_I5xz_G4y_a;
    abcd[iGrid*945+841] = 2.0E0*I_ESP_I4xyz_G4y_a;
    abcd[iGrid*945+842] = 2.0E0*I_ESP_I4x2z_G4y_a-1*I_ESP_G4x_G4y;
    abcd[iGrid*945+843] = 2.0E0*I_ESP_I3x2yz_G4y_a;
    abcd[iGrid*945+844] = 2.0E0*I_ESP_I3xy2z_G4y_a-1*I_ESP_G3xy_G4y;
    abcd[iGrid*945+845] = 2.0E0*I_ESP_I3x3z_G4y_a-2*I_ESP_G3xz_G4y;
    abcd[iGrid*945+846] = 2.0E0*I_ESP_I2x3yz_G4y_a;
    abcd[iGrid*945+847] = 2.0E0*I_ESP_I2x2y2z_G4y_a-1*I_ESP_G2x2y_G4y;
    abcd[iGrid*945+848] = 2.0E0*I_ESP_I2xy3z_G4y_a-2*I_ESP_G2xyz_G4y;
    abcd[iGrid*945+849] = 2.0E0*I_ESP_I2x4z_G4y_a-3*I_ESP_G2x2z_G4y;
    abcd[iGrid*945+850] = 2.0E0*I_ESP_Ix4yz_G4y_a;
    abcd[iGrid*945+851] = 2.0E0*I_ESP_Ix3y2z_G4y_a-1*I_ESP_Gx3y_G4y;
    abcd[iGrid*945+852] = 2.0E0*I_ESP_Ix2y3z_G4y_a-2*I_ESP_Gx2yz_G4y;
    abcd[iGrid*945+853] = 2.0E0*I_ESP_Ixy4z_G4y_a-3*I_ESP_Gxy2z_G4y;
    abcd[iGrid*945+854] = 2.0E0*I_ESP_Ix5z_G4y_a-4*I_ESP_Gx3z_G4y;
    abcd[iGrid*945+855] = 2.0E0*I_ESP_I5yz_G4y_a;
    abcd[iGrid*945+856] = 2.0E0*I_ESP_I4y2z_G4y_a-1*I_ESP_G4y_G4y;
    abcd[iGrid*945+857] = 2.0E0*I_ESP_I3y3z_G4y_a-2*I_ESP_G3yz_G4y;
    abcd[iGrid*945+858] = 2.0E0*I_ESP_I2y4z_G4y_a-3*I_ESP_G2y2z_G4y;
    abcd[iGrid*945+859] = 2.0E0*I_ESP_Iy5z_G4y_a-4*I_ESP_Gy3z_G4y;
    abcd[iGrid*945+860] = 2.0E0*I_ESP_I6z_G4y_a-5*I_ESP_G4z_G4y;
    abcd[iGrid*945+861] = 2.0E0*I_ESP_I5xz_G3yz_a;
    abcd[iGrid*945+862] = 2.0E0*I_ESP_I4xyz_G3yz_a;
    abcd[iGrid*945+863] = 2.0E0*I_ESP_I4x2z_G3yz_a-1*I_ESP_G4x_G3yz;
    abcd[iGrid*945+864] = 2.0E0*I_ESP_I3x2yz_G3yz_a;
    abcd[iGrid*945+865] = 2.0E0*I_ESP_I3xy2z_G3yz_a-1*I_ESP_G3xy_G3yz;
    abcd[iGrid*945+866] = 2.0E0*I_ESP_I3x3z_G3yz_a-2*I_ESP_G3xz_G3yz;
    abcd[iGrid*945+867] = 2.0E0*I_ESP_I2x3yz_G3yz_a;
    abcd[iGrid*945+868] = 2.0E0*I_ESP_I2x2y2z_G3yz_a-1*I_ESP_G2x2y_G3yz;
    abcd[iGrid*945+869] = 2.0E0*I_ESP_I2xy3z_G3yz_a-2*I_ESP_G2xyz_G3yz;
    abcd[iGrid*945+870] = 2.0E0*I_ESP_I2x4z_G3yz_a-3*I_ESP_G2x2z_G3yz;
    abcd[iGrid*945+871] = 2.0E0*I_ESP_Ix4yz_G3yz_a;
    abcd[iGrid*945+872] = 2.0E0*I_ESP_Ix3y2z_G3yz_a-1*I_ESP_Gx3y_G3yz;
    abcd[iGrid*945+873] = 2.0E0*I_ESP_Ix2y3z_G3yz_a-2*I_ESP_Gx2yz_G3yz;
    abcd[iGrid*945+874] = 2.0E0*I_ESP_Ixy4z_G3yz_a-3*I_ESP_Gxy2z_G3yz;
    abcd[iGrid*945+875] = 2.0E0*I_ESP_Ix5z_G3yz_a-4*I_ESP_Gx3z_G3yz;
    abcd[iGrid*945+876] = 2.0E0*I_ESP_I5yz_G3yz_a;
    abcd[iGrid*945+877] = 2.0E0*I_ESP_I4y2z_G3yz_a-1*I_ESP_G4y_G3yz;
    abcd[iGrid*945+878] = 2.0E0*I_ESP_I3y3z_G3yz_a-2*I_ESP_G3yz_G3yz;
    abcd[iGrid*945+879] = 2.0E0*I_ESP_I2y4z_G3yz_a-3*I_ESP_G2y2z_G3yz;
    abcd[iGrid*945+880] = 2.0E0*I_ESP_Iy5z_G3yz_a-4*I_ESP_Gy3z_G3yz;
    abcd[iGrid*945+881] = 2.0E0*I_ESP_I6z_G3yz_a-5*I_ESP_G4z_G3yz;
    abcd[iGrid*945+882] = 2.0E0*I_ESP_I5xz_G2y2z_a;
    abcd[iGrid*945+883] = 2.0E0*I_ESP_I4xyz_G2y2z_a;
    abcd[iGrid*945+884] = 2.0E0*I_ESP_I4x2z_G2y2z_a-1*I_ESP_G4x_G2y2z;
    abcd[iGrid*945+885] = 2.0E0*I_ESP_I3x2yz_G2y2z_a;
    abcd[iGrid*945+886] = 2.0E0*I_ESP_I3xy2z_G2y2z_a-1*I_ESP_G3xy_G2y2z;
    abcd[iGrid*945+887] = 2.0E0*I_ESP_I3x3z_G2y2z_a-2*I_ESP_G3xz_G2y2z;
    abcd[iGrid*945+888] = 2.0E0*I_ESP_I2x3yz_G2y2z_a;
    abcd[iGrid*945+889] = 2.0E0*I_ESP_I2x2y2z_G2y2z_a-1*I_ESP_G2x2y_G2y2z;
    abcd[iGrid*945+890] = 2.0E0*I_ESP_I2xy3z_G2y2z_a-2*I_ESP_G2xyz_G2y2z;
    abcd[iGrid*945+891] = 2.0E0*I_ESP_I2x4z_G2y2z_a-3*I_ESP_G2x2z_G2y2z;
    abcd[iGrid*945+892] = 2.0E0*I_ESP_Ix4yz_G2y2z_a;
    abcd[iGrid*945+893] = 2.0E0*I_ESP_Ix3y2z_G2y2z_a-1*I_ESP_Gx3y_G2y2z;
    abcd[iGrid*945+894] = 2.0E0*I_ESP_Ix2y3z_G2y2z_a-2*I_ESP_Gx2yz_G2y2z;
    abcd[iGrid*945+895] = 2.0E0*I_ESP_Ixy4z_G2y2z_a-3*I_ESP_Gxy2z_G2y2z;
    abcd[iGrid*945+896] = 2.0E0*I_ESP_Ix5z_G2y2z_a-4*I_ESP_Gx3z_G2y2z;
    abcd[iGrid*945+897] = 2.0E0*I_ESP_I5yz_G2y2z_a;
    abcd[iGrid*945+898] = 2.0E0*I_ESP_I4y2z_G2y2z_a-1*I_ESP_G4y_G2y2z;
    abcd[iGrid*945+899] = 2.0E0*I_ESP_I3y3z_G2y2z_a-2*I_ESP_G3yz_G2y2z;
    abcd[iGrid*945+900] = 2.0E0*I_ESP_I2y4z_G2y2z_a-3*I_ESP_G2y2z_G2y2z;
    abcd[iGrid*945+901] = 2.0E0*I_ESP_Iy5z_G2y2z_a-4*I_ESP_Gy3z_G2y2z;
    abcd[iGrid*945+902] = 2.0E0*I_ESP_I6z_G2y2z_a-5*I_ESP_G4z_G2y2z;
    abcd[iGrid*945+903] = 2.0E0*I_ESP_I5xz_Gy3z_a;
    abcd[iGrid*945+904] = 2.0E0*I_ESP_I4xyz_Gy3z_a;
    abcd[iGrid*945+905] = 2.0E0*I_ESP_I4x2z_Gy3z_a-1*I_ESP_G4x_Gy3z;
    abcd[iGrid*945+906] = 2.0E0*I_ESP_I3x2yz_Gy3z_a;
    abcd[iGrid*945+907] = 2.0E0*I_ESP_I3xy2z_Gy3z_a-1*I_ESP_G3xy_Gy3z;
    abcd[iGrid*945+908] = 2.0E0*I_ESP_I3x3z_Gy3z_a-2*I_ESP_G3xz_Gy3z;
    abcd[iGrid*945+909] = 2.0E0*I_ESP_I2x3yz_Gy3z_a;
    abcd[iGrid*945+910] = 2.0E0*I_ESP_I2x2y2z_Gy3z_a-1*I_ESP_G2x2y_Gy3z;
    abcd[iGrid*945+911] = 2.0E0*I_ESP_I2xy3z_Gy3z_a-2*I_ESP_G2xyz_Gy3z;
    abcd[iGrid*945+912] = 2.0E0*I_ESP_I2x4z_Gy3z_a-3*I_ESP_G2x2z_Gy3z;
    abcd[iGrid*945+913] = 2.0E0*I_ESP_Ix4yz_Gy3z_a;
    abcd[iGrid*945+914] = 2.0E0*I_ESP_Ix3y2z_Gy3z_a-1*I_ESP_Gx3y_Gy3z;
    abcd[iGrid*945+915] = 2.0E0*I_ESP_Ix2y3z_Gy3z_a-2*I_ESP_Gx2yz_Gy3z;
    abcd[iGrid*945+916] = 2.0E0*I_ESP_Ixy4z_Gy3z_a-3*I_ESP_Gxy2z_Gy3z;
    abcd[iGrid*945+917] = 2.0E0*I_ESP_Ix5z_Gy3z_a-4*I_ESP_Gx3z_Gy3z;
    abcd[iGrid*945+918] = 2.0E0*I_ESP_I5yz_Gy3z_a;
    abcd[iGrid*945+919] = 2.0E0*I_ESP_I4y2z_Gy3z_a-1*I_ESP_G4y_Gy3z;
    abcd[iGrid*945+920] = 2.0E0*I_ESP_I3y3z_Gy3z_a-2*I_ESP_G3yz_Gy3z;
    abcd[iGrid*945+921] = 2.0E0*I_ESP_I2y4z_Gy3z_a-3*I_ESP_G2y2z_Gy3z;
    abcd[iGrid*945+922] = 2.0E0*I_ESP_Iy5z_Gy3z_a-4*I_ESP_Gy3z_Gy3z;
    abcd[iGrid*945+923] = 2.0E0*I_ESP_I6z_Gy3z_a-5*I_ESP_G4z_Gy3z;
    abcd[iGrid*945+924] = 2.0E0*I_ESP_I5xz_G4z_a;
    abcd[iGrid*945+925] = 2.0E0*I_ESP_I4xyz_G4z_a;
    abcd[iGrid*945+926] = 2.0E0*I_ESP_I4x2z_G4z_a-1*I_ESP_G4x_G4z;
    abcd[iGrid*945+927] = 2.0E0*I_ESP_I3x2yz_G4z_a;
    abcd[iGrid*945+928] = 2.0E0*I_ESP_I3xy2z_G4z_a-1*I_ESP_G3xy_G4z;
    abcd[iGrid*945+929] = 2.0E0*I_ESP_I3x3z_G4z_a-2*I_ESP_G3xz_G4z;
    abcd[iGrid*945+930] = 2.0E0*I_ESP_I2x3yz_G4z_a;
    abcd[iGrid*945+931] = 2.0E0*I_ESP_I2x2y2z_G4z_a-1*I_ESP_G2x2y_G4z;
    abcd[iGrid*945+932] = 2.0E0*I_ESP_I2xy3z_G4z_a-2*I_ESP_G2xyz_G4z;
    abcd[iGrid*945+933] = 2.0E0*I_ESP_I2x4z_G4z_a-3*I_ESP_G2x2z_G4z;
    abcd[iGrid*945+934] = 2.0E0*I_ESP_Ix4yz_G4z_a;
    abcd[iGrid*945+935] = 2.0E0*I_ESP_Ix3y2z_G4z_a-1*I_ESP_Gx3y_G4z;
    abcd[iGrid*945+936] = 2.0E0*I_ESP_Ix2y3z_G4z_a-2*I_ESP_Gx2yz_G4z;
    abcd[iGrid*945+937] = 2.0E0*I_ESP_Ixy4z_G4z_a-3*I_ESP_Gxy2z_G4z;
    abcd[iGrid*945+938] = 2.0E0*I_ESP_Ix5z_G4z_a-4*I_ESP_Gx3z_G4z;
    abcd[iGrid*945+939] = 2.0E0*I_ESP_I5yz_G4z_a;
    abcd[iGrid*945+940] = 2.0E0*I_ESP_I4y2z_G4z_a-1*I_ESP_G4y_G4z;
    abcd[iGrid*945+941] = 2.0E0*I_ESP_I3y3z_G4z_a-2*I_ESP_G3yz_G4z;
    abcd[iGrid*945+942] = 2.0E0*I_ESP_I2y4z_G4z_a-3*I_ESP_G2y2z_G4z;
    abcd[iGrid*945+943] = 2.0E0*I_ESP_Iy5z_G4z_a-4*I_ESP_Gy3z_G4z;
    abcd[iGrid*945+944] = 2.0E0*I_ESP_I6z_G4z_a-5*I_ESP_G4z_G4z;
  }
}
