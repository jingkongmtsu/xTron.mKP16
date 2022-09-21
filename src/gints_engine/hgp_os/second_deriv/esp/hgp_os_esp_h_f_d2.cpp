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
// BRA1 as redundant position, total RHS integrals evaluated as: 11614
// BRA2 as redundant position, total RHS integrals evaluated as: 10070
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

void hgp_os_esp_h_f_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
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
     * totally 20 integrals are omitted 
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
     * totally 14 integrals are omitted 
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
    Double I_ESP_H4yz_Px = I_ESP_Ix4yz_S+ABX*I_ESP_H4yz_S;
    Double I_ESP_H3y2z_Px = I_ESP_Ix3y2z_S+ABX*I_ESP_H3y2z_S;
    Double I_ESP_H2y3z_Px = I_ESP_Ix2y3z_S+ABX*I_ESP_H2y3z_S;
    Double I_ESP_Hy4z_Px = I_ESP_Ixy4z_S+ABX*I_ESP_Hy4z_S;
    Double I_ESP_H4xy_Py = I_ESP_I4x2y_S+ABY*I_ESP_H4xy_S;
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
    Double I_ESP_H4xz_Pz = I_ESP_I4x2z_S+ABZ*I_ESP_H4xz_S;
    Double I_ESP_H3xyz_Pz = I_ESP_I3xy2z_S+ABZ*I_ESP_H3xyz_S;
    Double I_ESP_H3x2z_Pz = I_ESP_I3x3z_S+ABZ*I_ESP_H3x2z_S;
    Double I_ESP_H2x2yz_Pz = I_ESP_I2x2y2z_S+ABZ*I_ESP_H2x2yz_S;
    Double I_ESP_H2xy2z_Pz = I_ESP_I2xy3z_S+ABZ*I_ESP_H2xy2z_S;
    Double I_ESP_H2x3z_Pz = I_ESP_I2x4z_S+ABZ*I_ESP_H2x3z_S;
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
     * totally 35 integrals are omitted 
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
    Double I_ESP_G3xz_Dxy = I_ESP_H3xyz_Px+ABY*I_ESP_G3xz_Px;
    Double I_ESP_G2xyz_Dxy = I_ESP_H2x2yz_Px+ABY*I_ESP_G2xyz_Px;
    Double I_ESP_G2x2z_Dxy = I_ESP_H2xy2z_Px+ABY*I_ESP_G2x2z_Px;
    Double I_ESP_Gx2yz_Dxy = I_ESP_Hx3yz_Px+ABY*I_ESP_Gx2yz_Px;
    Double I_ESP_Gxy2z_Dxy = I_ESP_Hx2y2z_Px+ABY*I_ESP_Gxy2z_Px;
    Double I_ESP_Gx3z_Dxy = I_ESP_Hxy3z_Px+ABY*I_ESP_Gx3z_Px;
    Double I_ESP_G3yz_Dxy = I_ESP_H4yz_Px+ABY*I_ESP_G3yz_Px;
    Double I_ESP_G2y2z_Dxy = I_ESP_H3y2z_Px+ABY*I_ESP_G2y2z_Px;
    Double I_ESP_Gy3z_Dxy = I_ESP_H2y3z_Px+ABY*I_ESP_Gy3z_Px;
    Double I_ESP_G4z_Dxy = I_ESP_Hy4z_Px+ABY*I_ESP_G4z_Px;
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
     * totally 0 integrals are omitted 
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
    Double I_ESP_F3x_Fxyz = I_ESP_G3xz_Dxy+ABZ*I_ESP_F3x_Dxy;
    Double I_ESP_F2xy_Fxyz = I_ESP_G2xyz_Dxy+ABZ*I_ESP_F2xy_Dxy;
    Double I_ESP_F2xz_Fxyz = I_ESP_G2x2z_Dxy+ABZ*I_ESP_F2xz_Dxy;
    Double I_ESP_Fx2y_Fxyz = I_ESP_Gx2yz_Dxy+ABZ*I_ESP_Fx2y_Dxy;
    Double I_ESP_Fxyz_Fxyz = I_ESP_Gxy2z_Dxy+ABZ*I_ESP_Fxyz_Dxy;
    Double I_ESP_Fx2z_Fxyz = I_ESP_Gx3z_Dxy+ABZ*I_ESP_Fx2z_Dxy;
    Double I_ESP_F3y_Fxyz = I_ESP_G3yz_Dxy+ABZ*I_ESP_F3y_Dxy;
    Double I_ESP_F2yz_Fxyz = I_ESP_G2y2z_Dxy+ABZ*I_ESP_F2yz_Dxy;
    Double I_ESP_Fy2z_Fxyz = I_ESP_Gy3z_Dxy+ABZ*I_ESP_Fy2z_Dxy;
    Double I_ESP_F3z_Fxyz = I_ESP_G4z_Dxy+ABZ*I_ESP_F3z_Dxy;
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
    Double I_ESP_F3x_Fy2z = I_ESP_G3xy_D2z+ABY*I_ESP_F3x_D2z;
    Double I_ESP_F2xy_Fy2z = I_ESP_G2x2y_D2z+ABY*I_ESP_F2xy_D2z;
    Double I_ESP_F2xz_Fy2z = I_ESP_G2xyz_D2z+ABY*I_ESP_F2xz_D2z;
    Double I_ESP_Fx2y_Fy2z = I_ESP_Gx3y_D2z+ABY*I_ESP_Fx2y_D2z;
    Double I_ESP_Fxyz_Fy2z = I_ESP_Gx2yz_D2z+ABY*I_ESP_Fxyz_D2z;
    Double I_ESP_Fx2z_Fy2z = I_ESP_Gxy2z_D2z+ABY*I_ESP_Fx2z_D2z;
    Double I_ESP_F3y_Fy2z = I_ESP_G4y_D2z+ABY*I_ESP_F3y_D2z;
    Double I_ESP_F2yz_Fy2z = I_ESP_G3yz_D2z+ABY*I_ESP_F2yz_D2z;
    Double I_ESP_Fy2z_Fy2z = I_ESP_G2y2z_D2z+ABY*I_ESP_Fy2z_D2z;
    Double I_ESP_F3z_Fy2z = I_ESP_Gy3z_D2z+ABY*I_ESP_F3z_D2z;
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
     * totally 42 integrals are omitted 
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
     * totally 18 integrals are omitted 
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
    Double I_ESP_K6yz_Px_a = I_ESP_Lx6yz_S_a+ABX*I_ESP_K6yz_S_a;
    Double I_ESP_K5y2z_Px_a = I_ESP_Lx5y2z_S_a+ABX*I_ESP_K5y2z_S_a;
    Double I_ESP_K4y3z_Px_a = I_ESP_Lx4y3z_S_a+ABX*I_ESP_K4y3z_S_a;
    Double I_ESP_K3y4z_Px_a = I_ESP_Lx3y4z_S_a+ABX*I_ESP_K3y4z_S_a;
    Double I_ESP_K2y5z_Px_a = I_ESP_Lx2y5z_S_a+ABX*I_ESP_K2y5z_S_a;
    Double I_ESP_Ky6z_Px_a = I_ESP_Lxy6z_S_a+ABX*I_ESP_Ky6z_S_a;
    Double I_ESP_K6xy_Py_a = I_ESP_L6x2y_S_a+ABY*I_ESP_K6xy_S_a;
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
    Double I_ESP_K6xz_Pz_a = I_ESP_L6x2z_S_a+ABZ*I_ESP_K6xz_S_a;
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
     * totally 63 integrals are omitted 
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
    Double I_ESP_I5xz_Dxy_a = I_ESP_K5xyz_Px_a+ABY*I_ESP_I5xz_Px_a;
    Double I_ESP_I4xyz_Dxy_a = I_ESP_K4x2yz_Px_a+ABY*I_ESP_I4xyz_Px_a;
    Double I_ESP_I4x2z_Dxy_a = I_ESP_K4xy2z_Px_a+ABY*I_ESP_I4x2z_Px_a;
    Double I_ESP_I3x2yz_Dxy_a = I_ESP_K3x3yz_Px_a+ABY*I_ESP_I3x2yz_Px_a;
    Double I_ESP_I3xy2z_Dxy_a = I_ESP_K3x2y2z_Px_a+ABY*I_ESP_I3xy2z_Px_a;
    Double I_ESP_I3x3z_Dxy_a = I_ESP_K3xy3z_Px_a+ABY*I_ESP_I3x3z_Px_a;
    Double I_ESP_I2x3yz_Dxy_a = I_ESP_K2x4yz_Px_a+ABY*I_ESP_I2x3yz_Px_a;
    Double I_ESP_I2x2y2z_Dxy_a = I_ESP_K2x3y2z_Px_a+ABY*I_ESP_I2x2y2z_Px_a;
    Double I_ESP_I2xy3z_Dxy_a = I_ESP_K2x2y3z_Px_a+ABY*I_ESP_I2xy3z_Px_a;
    Double I_ESP_I2x4z_Dxy_a = I_ESP_K2xy4z_Px_a+ABY*I_ESP_I2x4z_Px_a;
    Double I_ESP_Ix4yz_Dxy_a = I_ESP_Kx5yz_Px_a+ABY*I_ESP_Ix4yz_Px_a;
    Double I_ESP_Ix3y2z_Dxy_a = I_ESP_Kx4y2z_Px_a+ABY*I_ESP_Ix3y2z_Px_a;
    Double I_ESP_Ix2y3z_Dxy_a = I_ESP_Kx3y3z_Px_a+ABY*I_ESP_Ix2y3z_Px_a;
    Double I_ESP_Ixy4z_Dxy_a = I_ESP_Kx2y4z_Px_a+ABY*I_ESP_Ixy4z_Px_a;
    Double I_ESP_Ix5z_Dxy_a = I_ESP_Kxy5z_Px_a+ABY*I_ESP_Ix5z_Px_a;
    Double I_ESP_I5yz_Dxy_a = I_ESP_K6yz_Px_a+ABY*I_ESP_I5yz_Px_a;
    Double I_ESP_I4y2z_Dxy_a = I_ESP_K5y2z_Px_a+ABY*I_ESP_I4y2z_Px_a;
    Double I_ESP_I3y3z_Dxy_a = I_ESP_K4y3z_Px_a+ABY*I_ESP_I3y3z_Px_a;
    Double I_ESP_I2y4z_Dxy_a = I_ESP_K3y4z_Px_a+ABY*I_ESP_I2y4z_Px_a;
    Double I_ESP_Iy5z_Dxy_a = I_ESP_K2y5z_Px_a+ABY*I_ESP_Iy5z_Px_a;
    Double I_ESP_I6z_Dxy_a = I_ESP_Ky6z_Px_a+ABY*I_ESP_I6z_Px_a;
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
     * totally 0 integrals are omitted 
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
    Double I_ESP_H5x_Fxyz_a = I_ESP_I5xz_Dxy_a+ABZ*I_ESP_H5x_Dxy_a;
    Double I_ESP_H4xy_Fxyz_a = I_ESP_I4xyz_Dxy_a+ABZ*I_ESP_H4xy_Dxy_a;
    Double I_ESP_H4xz_Fxyz_a = I_ESP_I4x2z_Dxy_a+ABZ*I_ESP_H4xz_Dxy_a;
    Double I_ESP_H3x2y_Fxyz_a = I_ESP_I3x2yz_Dxy_a+ABZ*I_ESP_H3x2y_Dxy_a;
    Double I_ESP_H3xyz_Fxyz_a = I_ESP_I3xy2z_Dxy_a+ABZ*I_ESP_H3xyz_Dxy_a;
    Double I_ESP_H3x2z_Fxyz_a = I_ESP_I3x3z_Dxy_a+ABZ*I_ESP_H3x2z_Dxy_a;
    Double I_ESP_H2x3y_Fxyz_a = I_ESP_I2x3yz_Dxy_a+ABZ*I_ESP_H2x3y_Dxy_a;
    Double I_ESP_H2x2yz_Fxyz_a = I_ESP_I2x2y2z_Dxy_a+ABZ*I_ESP_H2x2yz_Dxy_a;
    Double I_ESP_H2xy2z_Fxyz_a = I_ESP_I2xy3z_Dxy_a+ABZ*I_ESP_H2xy2z_Dxy_a;
    Double I_ESP_H2x3z_Fxyz_a = I_ESP_I2x4z_Dxy_a+ABZ*I_ESP_H2x3z_Dxy_a;
    Double I_ESP_Hx4y_Fxyz_a = I_ESP_Ix4yz_Dxy_a+ABZ*I_ESP_Hx4y_Dxy_a;
    Double I_ESP_Hx3yz_Fxyz_a = I_ESP_Ix3y2z_Dxy_a+ABZ*I_ESP_Hx3yz_Dxy_a;
    Double I_ESP_Hx2y2z_Fxyz_a = I_ESP_Ix2y3z_Dxy_a+ABZ*I_ESP_Hx2y2z_Dxy_a;
    Double I_ESP_Hxy3z_Fxyz_a = I_ESP_Ixy4z_Dxy_a+ABZ*I_ESP_Hxy3z_Dxy_a;
    Double I_ESP_Hx4z_Fxyz_a = I_ESP_Ix5z_Dxy_a+ABZ*I_ESP_Hx4z_Dxy_a;
    Double I_ESP_H5y_Fxyz_a = I_ESP_I5yz_Dxy_a+ABZ*I_ESP_H5y_Dxy_a;
    Double I_ESP_H4yz_Fxyz_a = I_ESP_I4y2z_Dxy_a+ABZ*I_ESP_H4yz_Dxy_a;
    Double I_ESP_H3y2z_Fxyz_a = I_ESP_I3y3z_Dxy_a+ABZ*I_ESP_H3y2z_Dxy_a;
    Double I_ESP_H2y3z_Fxyz_a = I_ESP_I2y4z_Dxy_a+ABZ*I_ESP_H2y3z_Dxy_a;
    Double I_ESP_Hy4z_Fxyz_a = I_ESP_Iy5z_Dxy_a+ABZ*I_ESP_Hy4z_Dxy_a;
    Double I_ESP_H5z_Fxyz_a = I_ESP_I6z_Dxy_a+ABZ*I_ESP_H5z_Dxy_a;
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
    Double I_ESP_H5x_Fy2z_a = I_ESP_I5xy_D2z_a+ABY*I_ESP_H5x_D2z_a;
    Double I_ESP_H4xy_Fy2z_a = I_ESP_I4x2y_D2z_a+ABY*I_ESP_H4xy_D2z_a;
    Double I_ESP_H4xz_Fy2z_a = I_ESP_I4xyz_D2z_a+ABY*I_ESP_H4xz_D2z_a;
    Double I_ESP_H3x2y_Fy2z_a = I_ESP_I3x3y_D2z_a+ABY*I_ESP_H3x2y_D2z_a;
    Double I_ESP_H3xyz_Fy2z_a = I_ESP_I3x2yz_D2z_a+ABY*I_ESP_H3xyz_D2z_a;
    Double I_ESP_H3x2z_Fy2z_a = I_ESP_I3xy2z_D2z_a+ABY*I_ESP_H3x2z_D2z_a;
    Double I_ESP_H2x3y_Fy2z_a = I_ESP_I2x4y_D2z_a+ABY*I_ESP_H2x3y_D2z_a;
    Double I_ESP_H2x2yz_Fy2z_a = I_ESP_I2x3yz_D2z_a+ABY*I_ESP_H2x2yz_D2z_a;
    Double I_ESP_H2xy2z_Fy2z_a = I_ESP_I2x2y2z_D2z_a+ABY*I_ESP_H2xy2z_D2z_a;
    Double I_ESP_H2x3z_Fy2z_a = I_ESP_I2xy3z_D2z_a+ABY*I_ESP_H2x3z_D2z_a;
    Double I_ESP_Hx4y_Fy2z_a = I_ESP_Ix5y_D2z_a+ABY*I_ESP_Hx4y_D2z_a;
    Double I_ESP_Hx3yz_Fy2z_a = I_ESP_Ix4yz_D2z_a+ABY*I_ESP_Hx3yz_D2z_a;
    Double I_ESP_Hx2y2z_Fy2z_a = I_ESP_Ix3y2z_D2z_a+ABY*I_ESP_Hx2y2z_D2z_a;
    Double I_ESP_Hxy3z_Fy2z_a = I_ESP_Ix2y3z_D2z_a+ABY*I_ESP_Hxy3z_D2z_a;
    Double I_ESP_Hx4z_Fy2z_a = I_ESP_Ixy4z_D2z_a+ABY*I_ESP_Hx4z_D2z_a;
    Double I_ESP_H5y_Fy2z_a = I_ESP_I6y_D2z_a+ABY*I_ESP_H5y_D2z_a;
    Double I_ESP_H4yz_Fy2z_a = I_ESP_I5yz_D2z_a+ABY*I_ESP_H4yz_D2z_a;
    Double I_ESP_H3y2z_Fy2z_a = I_ESP_I4y2z_D2z_a+ABY*I_ESP_H3y2z_D2z_a;
    Double I_ESP_H2y3z_Fy2z_a = I_ESP_I3y3z_D2z_a+ABY*I_ESP_H2y3z_D2z_a;
    Double I_ESP_Hy4z_Fy2z_a = I_ESP_I2y4z_D2z_a+ABY*I_ESP_Hy4z_D2z_a;
    Double I_ESP_H5z_Fy2z_a = I_ESP_Iy5z_D2z_a+ABY*I_ESP_H5z_D2z_a;
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
     * totally 0 integrals are omitted 
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
    Double I_ESP_L8x_Py_aa = I_ESP_M8xy_S_aa+ABY*I_ESP_L8x_S_aa;
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
    Double I_ESP_L8x_Pz_aa = I_ESP_M8xz_S_aa+ABZ*I_ESP_L8x_S_aa;
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
    Double I_ESP_L8y_Pz_aa = I_ESP_M8yz_S_aa+ABZ*I_ESP_L8y_S_aa;
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
     * totally 72 integrals are omitted 
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
     * shell quartet name: SQ_ESP_M_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 22 integrals are omitted 
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
    Double I_ESP_M8yz_Px_aa = I_ESP_Nx8yz_S_aa+ABX*I_ESP_M8yz_S_aa;
    Double I_ESP_M7y2z_Px_aa = I_ESP_Nx7y2z_S_aa+ABX*I_ESP_M7y2z_S_aa;
    Double I_ESP_M6y3z_Px_aa = I_ESP_Nx6y3z_S_aa+ABX*I_ESP_M6y3z_S_aa;
    Double I_ESP_M5y4z_Px_aa = I_ESP_Nx5y4z_S_aa+ABX*I_ESP_M5y4z_S_aa;
    Double I_ESP_M4y5z_Px_aa = I_ESP_Nx4y5z_S_aa+ABX*I_ESP_M4y5z_S_aa;
    Double I_ESP_M3y6z_Px_aa = I_ESP_Nx3y6z_S_aa+ABX*I_ESP_M3y6z_S_aa;
    Double I_ESP_M2y7z_Px_aa = I_ESP_Nx2y7z_S_aa+ABX*I_ESP_M2y7z_S_aa;
    Double I_ESP_My8z_Px_aa = I_ESP_Nxy8z_S_aa+ABX*I_ESP_My8z_S_aa;
    Double I_ESP_M8xy_Py_aa = I_ESP_N8x2y_S_aa+ABY*I_ESP_M8xy_S_aa;
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
    Double I_ESP_M8xz_Pz_aa = I_ESP_N8x2z_S_aa+ABZ*I_ESP_M8xz_S_aa;
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
    Double I_ESP_M8yz_Pz_aa = I_ESP_N8y2z_S_aa+ABZ*I_ESP_M8yz_S_aa;
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
     * totally 99 integrals are omitted 
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
    Double I_ESP_L7xz_Dxy_aa = I_ESP_M7xyz_Px_aa+ABY*I_ESP_L7xz_Px_aa;
    Double I_ESP_L6xyz_Dxy_aa = I_ESP_M6x2yz_Px_aa+ABY*I_ESP_L6xyz_Px_aa;
    Double I_ESP_L6x2z_Dxy_aa = I_ESP_M6xy2z_Px_aa+ABY*I_ESP_L6x2z_Px_aa;
    Double I_ESP_L5x2yz_Dxy_aa = I_ESP_M5x3yz_Px_aa+ABY*I_ESP_L5x2yz_Px_aa;
    Double I_ESP_L5xy2z_Dxy_aa = I_ESP_M5x2y2z_Px_aa+ABY*I_ESP_L5xy2z_Px_aa;
    Double I_ESP_L5x3z_Dxy_aa = I_ESP_M5xy3z_Px_aa+ABY*I_ESP_L5x3z_Px_aa;
    Double I_ESP_L4x3yz_Dxy_aa = I_ESP_M4x4yz_Px_aa+ABY*I_ESP_L4x3yz_Px_aa;
    Double I_ESP_L4x2y2z_Dxy_aa = I_ESP_M4x3y2z_Px_aa+ABY*I_ESP_L4x2y2z_Px_aa;
    Double I_ESP_L4xy3z_Dxy_aa = I_ESP_M4x2y3z_Px_aa+ABY*I_ESP_L4xy3z_Px_aa;
    Double I_ESP_L4x4z_Dxy_aa = I_ESP_M4xy4z_Px_aa+ABY*I_ESP_L4x4z_Px_aa;
    Double I_ESP_L3x4yz_Dxy_aa = I_ESP_M3x5yz_Px_aa+ABY*I_ESP_L3x4yz_Px_aa;
    Double I_ESP_L3x3y2z_Dxy_aa = I_ESP_M3x4y2z_Px_aa+ABY*I_ESP_L3x3y2z_Px_aa;
    Double I_ESP_L3x2y3z_Dxy_aa = I_ESP_M3x3y3z_Px_aa+ABY*I_ESP_L3x2y3z_Px_aa;
    Double I_ESP_L3xy4z_Dxy_aa = I_ESP_M3x2y4z_Px_aa+ABY*I_ESP_L3xy4z_Px_aa;
    Double I_ESP_L3x5z_Dxy_aa = I_ESP_M3xy5z_Px_aa+ABY*I_ESP_L3x5z_Px_aa;
    Double I_ESP_L2x5yz_Dxy_aa = I_ESP_M2x6yz_Px_aa+ABY*I_ESP_L2x5yz_Px_aa;
    Double I_ESP_L2x4y2z_Dxy_aa = I_ESP_M2x5y2z_Px_aa+ABY*I_ESP_L2x4y2z_Px_aa;
    Double I_ESP_L2x3y3z_Dxy_aa = I_ESP_M2x4y3z_Px_aa+ABY*I_ESP_L2x3y3z_Px_aa;
    Double I_ESP_L2x2y4z_Dxy_aa = I_ESP_M2x3y4z_Px_aa+ABY*I_ESP_L2x2y4z_Px_aa;
    Double I_ESP_L2xy5z_Dxy_aa = I_ESP_M2x2y5z_Px_aa+ABY*I_ESP_L2xy5z_Px_aa;
    Double I_ESP_L2x6z_Dxy_aa = I_ESP_M2xy6z_Px_aa+ABY*I_ESP_L2x6z_Px_aa;
    Double I_ESP_Lx6yz_Dxy_aa = I_ESP_Mx7yz_Px_aa+ABY*I_ESP_Lx6yz_Px_aa;
    Double I_ESP_Lx5y2z_Dxy_aa = I_ESP_Mx6y2z_Px_aa+ABY*I_ESP_Lx5y2z_Px_aa;
    Double I_ESP_Lx4y3z_Dxy_aa = I_ESP_Mx5y3z_Px_aa+ABY*I_ESP_Lx4y3z_Px_aa;
    Double I_ESP_Lx3y4z_Dxy_aa = I_ESP_Mx4y4z_Px_aa+ABY*I_ESP_Lx3y4z_Px_aa;
    Double I_ESP_Lx2y5z_Dxy_aa = I_ESP_Mx3y5z_Px_aa+ABY*I_ESP_Lx2y5z_Px_aa;
    Double I_ESP_Lxy6z_Dxy_aa = I_ESP_Mx2y6z_Px_aa+ABY*I_ESP_Lxy6z_Px_aa;
    Double I_ESP_Lx7z_Dxy_aa = I_ESP_Mxy7z_Px_aa+ABY*I_ESP_Lx7z_Px_aa;
    Double I_ESP_L7yz_Dxy_aa = I_ESP_M8yz_Px_aa+ABY*I_ESP_L7yz_Px_aa;
    Double I_ESP_L6y2z_Dxy_aa = I_ESP_M7y2z_Px_aa+ABY*I_ESP_L6y2z_Px_aa;
    Double I_ESP_L5y3z_Dxy_aa = I_ESP_M6y3z_Px_aa+ABY*I_ESP_L5y3z_Px_aa;
    Double I_ESP_L4y4z_Dxy_aa = I_ESP_M5y4z_Px_aa+ABY*I_ESP_L4y4z_Px_aa;
    Double I_ESP_L3y5z_Dxy_aa = I_ESP_M4y5z_Px_aa+ABY*I_ESP_L3y5z_Px_aa;
    Double I_ESP_L2y6z_Dxy_aa = I_ESP_M3y6z_Px_aa+ABY*I_ESP_L2y6z_Px_aa;
    Double I_ESP_Ly7z_Dxy_aa = I_ESP_M2y7z_Px_aa+ABY*I_ESP_Ly7z_Px_aa;
    Double I_ESP_L8z_Dxy_aa = I_ESP_My8z_Px_aa+ABY*I_ESP_L8z_Px_aa;
    Double I_ESP_L8x_D2y_aa = I_ESP_M8xy_Py_aa+ABY*I_ESP_L8x_Py_aa;
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
    Double I_ESP_L8x_D2z_aa = I_ESP_M8xz_Pz_aa+ABZ*I_ESP_L8x_Pz_aa;
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
    Double I_ESP_L8y_D2z_aa = I_ESP_M8yz_Pz_aa+ABZ*I_ESP_L8y_Pz_aa;
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
     * totally 0 integrals are omitted 
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
    Double I_ESP_K7x_F2xy_aa = I_ESP_L7xy_D2x_aa+ABY*I_ESP_K7x_D2x_aa;
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
    Double I_ESP_K7x_F2xz_aa = I_ESP_L7xz_D2x_aa+ABZ*I_ESP_K7x_D2x_aa;
    Double I_ESP_K6xy_F2xz_aa = I_ESP_L6xyz_D2x_aa+ABZ*I_ESP_K6xy_D2x_aa;
    Double I_ESP_K6xz_F2xz_aa = I_ESP_L6x2z_D2x_aa+ABZ*I_ESP_K6xz_D2x_aa;
    Double I_ESP_K5x2y_F2xz_aa = I_ESP_L5x2yz_D2x_aa+ABZ*I_ESP_K5x2y_D2x_aa;
    Double I_ESP_K5xyz_F2xz_aa = I_ESP_L5xy2z_D2x_aa+ABZ*I_ESP_K5xyz_D2x_aa;
    Double I_ESP_K5x2z_F2xz_aa = I_ESP_L5x3z_D2x_aa+ABZ*I_ESP_K5x2z_D2x_aa;
    Double I_ESP_K4x3y_F2xz_aa = I_ESP_L4x3yz_D2x_aa+ABZ*I_ESP_K4x3y_D2x_aa;
    Double I_ESP_K4x2yz_F2xz_aa = I_ESP_L4x2y2z_D2x_aa+ABZ*I_ESP_K4x2yz_D2x_aa;
    Double I_ESP_K4xy2z_F2xz_aa = I_ESP_L4xy3z_D2x_aa+ABZ*I_ESP_K4xy2z_D2x_aa;
    Double I_ESP_K4x3z_F2xz_aa = I_ESP_L4x4z_D2x_aa+ABZ*I_ESP_K4x3z_D2x_aa;
    Double I_ESP_K3x4y_F2xz_aa = I_ESP_L3x4yz_D2x_aa+ABZ*I_ESP_K3x4y_D2x_aa;
    Double I_ESP_K3x3yz_F2xz_aa = I_ESP_L3x3y2z_D2x_aa+ABZ*I_ESP_K3x3yz_D2x_aa;
    Double I_ESP_K3x2y2z_F2xz_aa = I_ESP_L3x2y3z_D2x_aa+ABZ*I_ESP_K3x2y2z_D2x_aa;
    Double I_ESP_K3xy3z_F2xz_aa = I_ESP_L3xy4z_D2x_aa+ABZ*I_ESP_K3xy3z_D2x_aa;
    Double I_ESP_K3x4z_F2xz_aa = I_ESP_L3x5z_D2x_aa+ABZ*I_ESP_K3x4z_D2x_aa;
    Double I_ESP_K2x5y_F2xz_aa = I_ESP_L2x5yz_D2x_aa+ABZ*I_ESP_K2x5y_D2x_aa;
    Double I_ESP_K2x4yz_F2xz_aa = I_ESP_L2x4y2z_D2x_aa+ABZ*I_ESP_K2x4yz_D2x_aa;
    Double I_ESP_K2x3y2z_F2xz_aa = I_ESP_L2x3y3z_D2x_aa+ABZ*I_ESP_K2x3y2z_D2x_aa;
    Double I_ESP_K2x2y3z_F2xz_aa = I_ESP_L2x2y4z_D2x_aa+ABZ*I_ESP_K2x2y3z_D2x_aa;
    Double I_ESP_K2xy4z_F2xz_aa = I_ESP_L2xy5z_D2x_aa+ABZ*I_ESP_K2xy4z_D2x_aa;
    Double I_ESP_K2x5z_F2xz_aa = I_ESP_L2x6z_D2x_aa+ABZ*I_ESP_K2x5z_D2x_aa;
    Double I_ESP_Kx6y_F2xz_aa = I_ESP_Lx6yz_D2x_aa+ABZ*I_ESP_Kx6y_D2x_aa;
    Double I_ESP_Kx5yz_F2xz_aa = I_ESP_Lx5y2z_D2x_aa+ABZ*I_ESP_Kx5yz_D2x_aa;
    Double I_ESP_Kx4y2z_F2xz_aa = I_ESP_Lx4y3z_D2x_aa+ABZ*I_ESP_Kx4y2z_D2x_aa;
    Double I_ESP_Kx3y3z_F2xz_aa = I_ESP_Lx3y4z_D2x_aa+ABZ*I_ESP_Kx3y3z_D2x_aa;
    Double I_ESP_Kx2y4z_F2xz_aa = I_ESP_Lx2y5z_D2x_aa+ABZ*I_ESP_Kx2y4z_D2x_aa;
    Double I_ESP_Kxy5z_F2xz_aa = I_ESP_Lxy6z_D2x_aa+ABZ*I_ESP_Kxy5z_D2x_aa;
    Double I_ESP_Kx6z_F2xz_aa = I_ESP_Lx7z_D2x_aa+ABZ*I_ESP_Kx6z_D2x_aa;
    Double I_ESP_K7y_F2xz_aa = I_ESP_L7yz_D2x_aa+ABZ*I_ESP_K7y_D2x_aa;
    Double I_ESP_K6yz_F2xz_aa = I_ESP_L6y2z_D2x_aa+ABZ*I_ESP_K6yz_D2x_aa;
    Double I_ESP_K5y2z_F2xz_aa = I_ESP_L5y3z_D2x_aa+ABZ*I_ESP_K5y2z_D2x_aa;
    Double I_ESP_K4y3z_F2xz_aa = I_ESP_L4y4z_D2x_aa+ABZ*I_ESP_K4y3z_D2x_aa;
    Double I_ESP_K3y4z_F2xz_aa = I_ESP_L3y5z_D2x_aa+ABZ*I_ESP_K3y4z_D2x_aa;
    Double I_ESP_K2y5z_F2xz_aa = I_ESP_L2y6z_D2x_aa+ABZ*I_ESP_K2y5z_D2x_aa;
    Double I_ESP_Ky6z_F2xz_aa = I_ESP_Ly7z_D2x_aa+ABZ*I_ESP_Ky6z_D2x_aa;
    Double I_ESP_K7z_F2xz_aa = I_ESP_L8z_D2x_aa+ABZ*I_ESP_K7z_D2x_aa;
    Double I_ESP_K7x_Fx2y_aa = I_ESP_L8x_D2y_aa+ABX*I_ESP_K7x_D2y_aa;
    Double I_ESP_K6xy_Fx2y_aa = I_ESP_L7xy_D2y_aa+ABX*I_ESP_K6xy_D2y_aa;
    Double I_ESP_K6xz_Fx2y_aa = I_ESP_L7xz_D2y_aa+ABX*I_ESP_K6xz_D2y_aa;
    Double I_ESP_K5x2y_Fx2y_aa = I_ESP_L6x2y_D2y_aa+ABX*I_ESP_K5x2y_D2y_aa;
    Double I_ESP_K5xyz_Fx2y_aa = I_ESP_L6xyz_D2y_aa+ABX*I_ESP_K5xyz_D2y_aa;
    Double I_ESP_K5x2z_Fx2y_aa = I_ESP_L6x2z_D2y_aa+ABX*I_ESP_K5x2z_D2y_aa;
    Double I_ESP_K4x3y_Fx2y_aa = I_ESP_L5x3y_D2y_aa+ABX*I_ESP_K4x3y_D2y_aa;
    Double I_ESP_K4x2yz_Fx2y_aa = I_ESP_L5x2yz_D2y_aa+ABX*I_ESP_K4x2yz_D2y_aa;
    Double I_ESP_K4xy2z_Fx2y_aa = I_ESP_L5xy2z_D2y_aa+ABX*I_ESP_K4xy2z_D2y_aa;
    Double I_ESP_K4x3z_Fx2y_aa = I_ESP_L5x3z_D2y_aa+ABX*I_ESP_K4x3z_D2y_aa;
    Double I_ESP_K3x4y_Fx2y_aa = I_ESP_L4x4y_D2y_aa+ABX*I_ESP_K3x4y_D2y_aa;
    Double I_ESP_K3x3yz_Fx2y_aa = I_ESP_L4x3yz_D2y_aa+ABX*I_ESP_K3x3yz_D2y_aa;
    Double I_ESP_K3x2y2z_Fx2y_aa = I_ESP_L4x2y2z_D2y_aa+ABX*I_ESP_K3x2y2z_D2y_aa;
    Double I_ESP_K3xy3z_Fx2y_aa = I_ESP_L4xy3z_D2y_aa+ABX*I_ESP_K3xy3z_D2y_aa;
    Double I_ESP_K3x4z_Fx2y_aa = I_ESP_L4x4z_D2y_aa+ABX*I_ESP_K3x4z_D2y_aa;
    Double I_ESP_K2x5y_Fx2y_aa = I_ESP_L3x5y_D2y_aa+ABX*I_ESP_K2x5y_D2y_aa;
    Double I_ESP_K2x4yz_Fx2y_aa = I_ESP_L3x4yz_D2y_aa+ABX*I_ESP_K2x4yz_D2y_aa;
    Double I_ESP_K2x3y2z_Fx2y_aa = I_ESP_L3x3y2z_D2y_aa+ABX*I_ESP_K2x3y2z_D2y_aa;
    Double I_ESP_K2x2y3z_Fx2y_aa = I_ESP_L3x2y3z_D2y_aa+ABX*I_ESP_K2x2y3z_D2y_aa;
    Double I_ESP_K2xy4z_Fx2y_aa = I_ESP_L3xy4z_D2y_aa+ABX*I_ESP_K2xy4z_D2y_aa;
    Double I_ESP_K2x5z_Fx2y_aa = I_ESP_L3x5z_D2y_aa+ABX*I_ESP_K2x5z_D2y_aa;
    Double I_ESP_Kx6y_Fx2y_aa = I_ESP_L2x6y_D2y_aa+ABX*I_ESP_Kx6y_D2y_aa;
    Double I_ESP_Kx5yz_Fx2y_aa = I_ESP_L2x5yz_D2y_aa+ABX*I_ESP_Kx5yz_D2y_aa;
    Double I_ESP_Kx4y2z_Fx2y_aa = I_ESP_L2x4y2z_D2y_aa+ABX*I_ESP_Kx4y2z_D2y_aa;
    Double I_ESP_Kx3y3z_Fx2y_aa = I_ESP_L2x3y3z_D2y_aa+ABX*I_ESP_Kx3y3z_D2y_aa;
    Double I_ESP_Kx2y4z_Fx2y_aa = I_ESP_L2x2y4z_D2y_aa+ABX*I_ESP_Kx2y4z_D2y_aa;
    Double I_ESP_Kxy5z_Fx2y_aa = I_ESP_L2xy5z_D2y_aa+ABX*I_ESP_Kxy5z_D2y_aa;
    Double I_ESP_Kx6z_Fx2y_aa = I_ESP_L2x6z_D2y_aa+ABX*I_ESP_Kx6z_D2y_aa;
    Double I_ESP_K7y_Fx2y_aa = I_ESP_Lx7y_D2y_aa+ABX*I_ESP_K7y_D2y_aa;
    Double I_ESP_K6yz_Fx2y_aa = I_ESP_Lx6yz_D2y_aa+ABX*I_ESP_K6yz_D2y_aa;
    Double I_ESP_K5y2z_Fx2y_aa = I_ESP_Lx5y2z_D2y_aa+ABX*I_ESP_K5y2z_D2y_aa;
    Double I_ESP_K4y3z_Fx2y_aa = I_ESP_Lx4y3z_D2y_aa+ABX*I_ESP_K4y3z_D2y_aa;
    Double I_ESP_K3y4z_Fx2y_aa = I_ESP_Lx3y4z_D2y_aa+ABX*I_ESP_K3y4z_D2y_aa;
    Double I_ESP_K2y5z_Fx2y_aa = I_ESP_Lx2y5z_D2y_aa+ABX*I_ESP_K2y5z_D2y_aa;
    Double I_ESP_Ky6z_Fx2y_aa = I_ESP_Lxy6z_D2y_aa+ABX*I_ESP_Ky6z_D2y_aa;
    Double I_ESP_K7z_Fx2y_aa = I_ESP_Lx7z_D2y_aa+ABX*I_ESP_K7z_D2y_aa;
    Double I_ESP_K7x_Fxyz_aa = I_ESP_L7xz_Dxy_aa+ABZ*I_ESP_K7x_Dxy_aa;
    Double I_ESP_K6xy_Fxyz_aa = I_ESP_L6xyz_Dxy_aa+ABZ*I_ESP_K6xy_Dxy_aa;
    Double I_ESP_K6xz_Fxyz_aa = I_ESP_L6x2z_Dxy_aa+ABZ*I_ESP_K6xz_Dxy_aa;
    Double I_ESP_K5x2y_Fxyz_aa = I_ESP_L5x2yz_Dxy_aa+ABZ*I_ESP_K5x2y_Dxy_aa;
    Double I_ESP_K5xyz_Fxyz_aa = I_ESP_L5xy2z_Dxy_aa+ABZ*I_ESP_K5xyz_Dxy_aa;
    Double I_ESP_K5x2z_Fxyz_aa = I_ESP_L5x3z_Dxy_aa+ABZ*I_ESP_K5x2z_Dxy_aa;
    Double I_ESP_K4x3y_Fxyz_aa = I_ESP_L4x3yz_Dxy_aa+ABZ*I_ESP_K4x3y_Dxy_aa;
    Double I_ESP_K4x2yz_Fxyz_aa = I_ESP_L4x2y2z_Dxy_aa+ABZ*I_ESP_K4x2yz_Dxy_aa;
    Double I_ESP_K4xy2z_Fxyz_aa = I_ESP_L4xy3z_Dxy_aa+ABZ*I_ESP_K4xy2z_Dxy_aa;
    Double I_ESP_K4x3z_Fxyz_aa = I_ESP_L4x4z_Dxy_aa+ABZ*I_ESP_K4x3z_Dxy_aa;
    Double I_ESP_K3x4y_Fxyz_aa = I_ESP_L3x4yz_Dxy_aa+ABZ*I_ESP_K3x4y_Dxy_aa;
    Double I_ESP_K3x3yz_Fxyz_aa = I_ESP_L3x3y2z_Dxy_aa+ABZ*I_ESP_K3x3yz_Dxy_aa;
    Double I_ESP_K3x2y2z_Fxyz_aa = I_ESP_L3x2y3z_Dxy_aa+ABZ*I_ESP_K3x2y2z_Dxy_aa;
    Double I_ESP_K3xy3z_Fxyz_aa = I_ESP_L3xy4z_Dxy_aa+ABZ*I_ESP_K3xy3z_Dxy_aa;
    Double I_ESP_K3x4z_Fxyz_aa = I_ESP_L3x5z_Dxy_aa+ABZ*I_ESP_K3x4z_Dxy_aa;
    Double I_ESP_K2x5y_Fxyz_aa = I_ESP_L2x5yz_Dxy_aa+ABZ*I_ESP_K2x5y_Dxy_aa;
    Double I_ESP_K2x4yz_Fxyz_aa = I_ESP_L2x4y2z_Dxy_aa+ABZ*I_ESP_K2x4yz_Dxy_aa;
    Double I_ESP_K2x3y2z_Fxyz_aa = I_ESP_L2x3y3z_Dxy_aa+ABZ*I_ESP_K2x3y2z_Dxy_aa;
    Double I_ESP_K2x2y3z_Fxyz_aa = I_ESP_L2x2y4z_Dxy_aa+ABZ*I_ESP_K2x2y3z_Dxy_aa;
    Double I_ESP_K2xy4z_Fxyz_aa = I_ESP_L2xy5z_Dxy_aa+ABZ*I_ESP_K2xy4z_Dxy_aa;
    Double I_ESP_K2x5z_Fxyz_aa = I_ESP_L2x6z_Dxy_aa+ABZ*I_ESP_K2x5z_Dxy_aa;
    Double I_ESP_Kx6y_Fxyz_aa = I_ESP_Lx6yz_Dxy_aa+ABZ*I_ESP_Kx6y_Dxy_aa;
    Double I_ESP_Kx5yz_Fxyz_aa = I_ESP_Lx5y2z_Dxy_aa+ABZ*I_ESP_Kx5yz_Dxy_aa;
    Double I_ESP_Kx4y2z_Fxyz_aa = I_ESP_Lx4y3z_Dxy_aa+ABZ*I_ESP_Kx4y2z_Dxy_aa;
    Double I_ESP_Kx3y3z_Fxyz_aa = I_ESP_Lx3y4z_Dxy_aa+ABZ*I_ESP_Kx3y3z_Dxy_aa;
    Double I_ESP_Kx2y4z_Fxyz_aa = I_ESP_Lx2y5z_Dxy_aa+ABZ*I_ESP_Kx2y4z_Dxy_aa;
    Double I_ESP_Kxy5z_Fxyz_aa = I_ESP_Lxy6z_Dxy_aa+ABZ*I_ESP_Kxy5z_Dxy_aa;
    Double I_ESP_Kx6z_Fxyz_aa = I_ESP_Lx7z_Dxy_aa+ABZ*I_ESP_Kx6z_Dxy_aa;
    Double I_ESP_K7y_Fxyz_aa = I_ESP_L7yz_Dxy_aa+ABZ*I_ESP_K7y_Dxy_aa;
    Double I_ESP_K6yz_Fxyz_aa = I_ESP_L6y2z_Dxy_aa+ABZ*I_ESP_K6yz_Dxy_aa;
    Double I_ESP_K5y2z_Fxyz_aa = I_ESP_L5y3z_Dxy_aa+ABZ*I_ESP_K5y2z_Dxy_aa;
    Double I_ESP_K4y3z_Fxyz_aa = I_ESP_L4y4z_Dxy_aa+ABZ*I_ESP_K4y3z_Dxy_aa;
    Double I_ESP_K3y4z_Fxyz_aa = I_ESP_L3y5z_Dxy_aa+ABZ*I_ESP_K3y4z_Dxy_aa;
    Double I_ESP_K2y5z_Fxyz_aa = I_ESP_L2y6z_Dxy_aa+ABZ*I_ESP_K2y5z_Dxy_aa;
    Double I_ESP_Ky6z_Fxyz_aa = I_ESP_Ly7z_Dxy_aa+ABZ*I_ESP_Ky6z_Dxy_aa;
    Double I_ESP_K7z_Fxyz_aa = I_ESP_L8z_Dxy_aa+ABZ*I_ESP_K7z_Dxy_aa;
    Double I_ESP_K7x_Fx2z_aa = I_ESP_L8x_D2z_aa+ABX*I_ESP_K7x_D2z_aa;
    Double I_ESP_K6xy_Fx2z_aa = I_ESP_L7xy_D2z_aa+ABX*I_ESP_K6xy_D2z_aa;
    Double I_ESP_K6xz_Fx2z_aa = I_ESP_L7xz_D2z_aa+ABX*I_ESP_K6xz_D2z_aa;
    Double I_ESP_K5x2y_Fx2z_aa = I_ESP_L6x2y_D2z_aa+ABX*I_ESP_K5x2y_D2z_aa;
    Double I_ESP_K5xyz_Fx2z_aa = I_ESP_L6xyz_D2z_aa+ABX*I_ESP_K5xyz_D2z_aa;
    Double I_ESP_K5x2z_Fx2z_aa = I_ESP_L6x2z_D2z_aa+ABX*I_ESP_K5x2z_D2z_aa;
    Double I_ESP_K4x3y_Fx2z_aa = I_ESP_L5x3y_D2z_aa+ABX*I_ESP_K4x3y_D2z_aa;
    Double I_ESP_K4x2yz_Fx2z_aa = I_ESP_L5x2yz_D2z_aa+ABX*I_ESP_K4x2yz_D2z_aa;
    Double I_ESP_K4xy2z_Fx2z_aa = I_ESP_L5xy2z_D2z_aa+ABX*I_ESP_K4xy2z_D2z_aa;
    Double I_ESP_K4x3z_Fx2z_aa = I_ESP_L5x3z_D2z_aa+ABX*I_ESP_K4x3z_D2z_aa;
    Double I_ESP_K3x4y_Fx2z_aa = I_ESP_L4x4y_D2z_aa+ABX*I_ESP_K3x4y_D2z_aa;
    Double I_ESP_K3x3yz_Fx2z_aa = I_ESP_L4x3yz_D2z_aa+ABX*I_ESP_K3x3yz_D2z_aa;
    Double I_ESP_K3x2y2z_Fx2z_aa = I_ESP_L4x2y2z_D2z_aa+ABX*I_ESP_K3x2y2z_D2z_aa;
    Double I_ESP_K3xy3z_Fx2z_aa = I_ESP_L4xy3z_D2z_aa+ABX*I_ESP_K3xy3z_D2z_aa;
    Double I_ESP_K3x4z_Fx2z_aa = I_ESP_L4x4z_D2z_aa+ABX*I_ESP_K3x4z_D2z_aa;
    Double I_ESP_K2x5y_Fx2z_aa = I_ESP_L3x5y_D2z_aa+ABX*I_ESP_K2x5y_D2z_aa;
    Double I_ESP_K2x4yz_Fx2z_aa = I_ESP_L3x4yz_D2z_aa+ABX*I_ESP_K2x4yz_D2z_aa;
    Double I_ESP_K2x3y2z_Fx2z_aa = I_ESP_L3x3y2z_D2z_aa+ABX*I_ESP_K2x3y2z_D2z_aa;
    Double I_ESP_K2x2y3z_Fx2z_aa = I_ESP_L3x2y3z_D2z_aa+ABX*I_ESP_K2x2y3z_D2z_aa;
    Double I_ESP_K2xy4z_Fx2z_aa = I_ESP_L3xy4z_D2z_aa+ABX*I_ESP_K2xy4z_D2z_aa;
    Double I_ESP_K2x5z_Fx2z_aa = I_ESP_L3x5z_D2z_aa+ABX*I_ESP_K2x5z_D2z_aa;
    Double I_ESP_Kx6y_Fx2z_aa = I_ESP_L2x6y_D2z_aa+ABX*I_ESP_Kx6y_D2z_aa;
    Double I_ESP_Kx5yz_Fx2z_aa = I_ESP_L2x5yz_D2z_aa+ABX*I_ESP_Kx5yz_D2z_aa;
    Double I_ESP_Kx4y2z_Fx2z_aa = I_ESP_L2x4y2z_D2z_aa+ABX*I_ESP_Kx4y2z_D2z_aa;
    Double I_ESP_Kx3y3z_Fx2z_aa = I_ESP_L2x3y3z_D2z_aa+ABX*I_ESP_Kx3y3z_D2z_aa;
    Double I_ESP_Kx2y4z_Fx2z_aa = I_ESP_L2x2y4z_D2z_aa+ABX*I_ESP_Kx2y4z_D2z_aa;
    Double I_ESP_Kxy5z_Fx2z_aa = I_ESP_L2xy5z_D2z_aa+ABX*I_ESP_Kxy5z_D2z_aa;
    Double I_ESP_Kx6z_Fx2z_aa = I_ESP_L2x6z_D2z_aa+ABX*I_ESP_Kx6z_D2z_aa;
    Double I_ESP_K7y_Fx2z_aa = I_ESP_Lx7y_D2z_aa+ABX*I_ESP_K7y_D2z_aa;
    Double I_ESP_K6yz_Fx2z_aa = I_ESP_Lx6yz_D2z_aa+ABX*I_ESP_K6yz_D2z_aa;
    Double I_ESP_K5y2z_Fx2z_aa = I_ESP_Lx5y2z_D2z_aa+ABX*I_ESP_K5y2z_D2z_aa;
    Double I_ESP_K4y3z_Fx2z_aa = I_ESP_Lx4y3z_D2z_aa+ABX*I_ESP_K4y3z_D2z_aa;
    Double I_ESP_K3y4z_Fx2z_aa = I_ESP_Lx3y4z_D2z_aa+ABX*I_ESP_K3y4z_D2z_aa;
    Double I_ESP_K2y5z_Fx2z_aa = I_ESP_Lx2y5z_D2z_aa+ABX*I_ESP_K2y5z_D2z_aa;
    Double I_ESP_Ky6z_Fx2z_aa = I_ESP_Lxy6z_D2z_aa+ABX*I_ESP_Ky6z_D2z_aa;
    Double I_ESP_K7z_Fx2z_aa = I_ESP_Lx7z_D2z_aa+ABX*I_ESP_K7z_D2z_aa;
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
    Double I_ESP_K7x_F2yz_aa = I_ESP_L7xz_D2y_aa+ABZ*I_ESP_K7x_D2y_aa;
    Double I_ESP_K6xy_F2yz_aa = I_ESP_L6xyz_D2y_aa+ABZ*I_ESP_K6xy_D2y_aa;
    Double I_ESP_K6xz_F2yz_aa = I_ESP_L6x2z_D2y_aa+ABZ*I_ESP_K6xz_D2y_aa;
    Double I_ESP_K5x2y_F2yz_aa = I_ESP_L5x2yz_D2y_aa+ABZ*I_ESP_K5x2y_D2y_aa;
    Double I_ESP_K5xyz_F2yz_aa = I_ESP_L5xy2z_D2y_aa+ABZ*I_ESP_K5xyz_D2y_aa;
    Double I_ESP_K5x2z_F2yz_aa = I_ESP_L5x3z_D2y_aa+ABZ*I_ESP_K5x2z_D2y_aa;
    Double I_ESP_K4x3y_F2yz_aa = I_ESP_L4x3yz_D2y_aa+ABZ*I_ESP_K4x3y_D2y_aa;
    Double I_ESP_K4x2yz_F2yz_aa = I_ESP_L4x2y2z_D2y_aa+ABZ*I_ESP_K4x2yz_D2y_aa;
    Double I_ESP_K4xy2z_F2yz_aa = I_ESP_L4xy3z_D2y_aa+ABZ*I_ESP_K4xy2z_D2y_aa;
    Double I_ESP_K4x3z_F2yz_aa = I_ESP_L4x4z_D2y_aa+ABZ*I_ESP_K4x3z_D2y_aa;
    Double I_ESP_K3x4y_F2yz_aa = I_ESP_L3x4yz_D2y_aa+ABZ*I_ESP_K3x4y_D2y_aa;
    Double I_ESP_K3x3yz_F2yz_aa = I_ESP_L3x3y2z_D2y_aa+ABZ*I_ESP_K3x3yz_D2y_aa;
    Double I_ESP_K3x2y2z_F2yz_aa = I_ESP_L3x2y3z_D2y_aa+ABZ*I_ESP_K3x2y2z_D2y_aa;
    Double I_ESP_K3xy3z_F2yz_aa = I_ESP_L3xy4z_D2y_aa+ABZ*I_ESP_K3xy3z_D2y_aa;
    Double I_ESP_K3x4z_F2yz_aa = I_ESP_L3x5z_D2y_aa+ABZ*I_ESP_K3x4z_D2y_aa;
    Double I_ESP_K2x5y_F2yz_aa = I_ESP_L2x5yz_D2y_aa+ABZ*I_ESP_K2x5y_D2y_aa;
    Double I_ESP_K2x4yz_F2yz_aa = I_ESP_L2x4y2z_D2y_aa+ABZ*I_ESP_K2x4yz_D2y_aa;
    Double I_ESP_K2x3y2z_F2yz_aa = I_ESP_L2x3y3z_D2y_aa+ABZ*I_ESP_K2x3y2z_D2y_aa;
    Double I_ESP_K2x2y3z_F2yz_aa = I_ESP_L2x2y4z_D2y_aa+ABZ*I_ESP_K2x2y3z_D2y_aa;
    Double I_ESP_K2xy4z_F2yz_aa = I_ESP_L2xy5z_D2y_aa+ABZ*I_ESP_K2xy4z_D2y_aa;
    Double I_ESP_K2x5z_F2yz_aa = I_ESP_L2x6z_D2y_aa+ABZ*I_ESP_K2x5z_D2y_aa;
    Double I_ESP_Kx6y_F2yz_aa = I_ESP_Lx6yz_D2y_aa+ABZ*I_ESP_Kx6y_D2y_aa;
    Double I_ESP_Kx5yz_F2yz_aa = I_ESP_Lx5y2z_D2y_aa+ABZ*I_ESP_Kx5yz_D2y_aa;
    Double I_ESP_Kx4y2z_F2yz_aa = I_ESP_Lx4y3z_D2y_aa+ABZ*I_ESP_Kx4y2z_D2y_aa;
    Double I_ESP_Kx3y3z_F2yz_aa = I_ESP_Lx3y4z_D2y_aa+ABZ*I_ESP_Kx3y3z_D2y_aa;
    Double I_ESP_Kx2y4z_F2yz_aa = I_ESP_Lx2y5z_D2y_aa+ABZ*I_ESP_Kx2y4z_D2y_aa;
    Double I_ESP_Kxy5z_F2yz_aa = I_ESP_Lxy6z_D2y_aa+ABZ*I_ESP_Kxy5z_D2y_aa;
    Double I_ESP_Kx6z_F2yz_aa = I_ESP_Lx7z_D2y_aa+ABZ*I_ESP_Kx6z_D2y_aa;
    Double I_ESP_K7y_F2yz_aa = I_ESP_L7yz_D2y_aa+ABZ*I_ESP_K7y_D2y_aa;
    Double I_ESP_K6yz_F2yz_aa = I_ESP_L6y2z_D2y_aa+ABZ*I_ESP_K6yz_D2y_aa;
    Double I_ESP_K5y2z_F2yz_aa = I_ESP_L5y3z_D2y_aa+ABZ*I_ESP_K5y2z_D2y_aa;
    Double I_ESP_K4y3z_F2yz_aa = I_ESP_L4y4z_D2y_aa+ABZ*I_ESP_K4y3z_D2y_aa;
    Double I_ESP_K3y4z_F2yz_aa = I_ESP_L3y5z_D2y_aa+ABZ*I_ESP_K3y4z_D2y_aa;
    Double I_ESP_K2y5z_F2yz_aa = I_ESP_L2y6z_D2y_aa+ABZ*I_ESP_K2y5z_D2y_aa;
    Double I_ESP_Ky6z_F2yz_aa = I_ESP_Ly7z_D2y_aa+ABZ*I_ESP_Ky6z_D2y_aa;
    Double I_ESP_K7z_F2yz_aa = I_ESP_L8z_D2y_aa+ABZ*I_ESP_K7z_D2y_aa;
    Double I_ESP_K7x_Fy2z_aa = I_ESP_L7xy_D2z_aa+ABY*I_ESP_K7x_D2z_aa;
    Double I_ESP_K6xy_Fy2z_aa = I_ESP_L6x2y_D2z_aa+ABY*I_ESP_K6xy_D2z_aa;
    Double I_ESP_K6xz_Fy2z_aa = I_ESP_L6xyz_D2z_aa+ABY*I_ESP_K6xz_D2z_aa;
    Double I_ESP_K5x2y_Fy2z_aa = I_ESP_L5x3y_D2z_aa+ABY*I_ESP_K5x2y_D2z_aa;
    Double I_ESP_K5xyz_Fy2z_aa = I_ESP_L5x2yz_D2z_aa+ABY*I_ESP_K5xyz_D2z_aa;
    Double I_ESP_K5x2z_Fy2z_aa = I_ESP_L5xy2z_D2z_aa+ABY*I_ESP_K5x2z_D2z_aa;
    Double I_ESP_K4x3y_Fy2z_aa = I_ESP_L4x4y_D2z_aa+ABY*I_ESP_K4x3y_D2z_aa;
    Double I_ESP_K4x2yz_Fy2z_aa = I_ESP_L4x3yz_D2z_aa+ABY*I_ESP_K4x2yz_D2z_aa;
    Double I_ESP_K4xy2z_Fy2z_aa = I_ESP_L4x2y2z_D2z_aa+ABY*I_ESP_K4xy2z_D2z_aa;
    Double I_ESP_K4x3z_Fy2z_aa = I_ESP_L4xy3z_D2z_aa+ABY*I_ESP_K4x3z_D2z_aa;
    Double I_ESP_K3x4y_Fy2z_aa = I_ESP_L3x5y_D2z_aa+ABY*I_ESP_K3x4y_D2z_aa;
    Double I_ESP_K3x3yz_Fy2z_aa = I_ESP_L3x4yz_D2z_aa+ABY*I_ESP_K3x3yz_D2z_aa;
    Double I_ESP_K3x2y2z_Fy2z_aa = I_ESP_L3x3y2z_D2z_aa+ABY*I_ESP_K3x2y2z_D2z_aa;
    Double I_ESP_K3xy3z_Fy2z_aa = I_ESP_L3x2y3z_D2z_aa+ABY*I_ESP_K3xy3z_D2z_aa;
    Double I_ESP_K3x4z_Fy2z_aa = I_ESP_L3xy4z_D2z_aa+ABY*I_ESP_K3x4z_D2z_aa;
    Double I_ESP_K2x5y_Fy2z_aa = I_ESP_L2x6y_D2z_aa+ABY*I_ESP_K2x5y_D2z_aa;
    Double I_ESP_K2x4yz_Fy2z_aa = I_ESP_L2x5yz_D2z_aa+ABY*I_ESP_K2x4yz_D2z_aa;
    Double I_ESP_K2x3y2z_Fy2z_aa = I_ESP_L2x4y2z_D2z_aa+ABY*I_ESP_K2x3y2z_D2z_aa;
    Double I_ESP_K2x2y3z_Fy2z_aa = I_ESP_L2x3y3z_D2z_aa+ABY*I_ESP_K2x2y3z_D2z_aa;
    Double I_ESP_K2xy4z_Fy2z_aa = I_ESP_L2x2y4z_D2z_aa+ABY*I_ESP_K2xy4z_D2z_aa;
    Double I_ESP_K2x5z_Fy2z_aa = I_ESP_L2xy5z_D2z_aa+ABY*I_ESP_K2x5z_D2z_aa;
    Double I_ESP_Kx6y_Fy2z_aa = I_ESP_Lx7y_D2z_aa+ABY*I_ESP_Kx6y_D2z_aa;
    Double I_ESP_Kx5yz_Fy2z_aa = I_ESP_Lx6yz_D2z_aa+ABY*I_ESP_Kx5yz_D2z_aa;
    Double I_ESP_Kx4y2z_Fy2z_aa = I_ESP_Lx5y2z_D2z_aa+ABY*I_ESP_Kx4y2z_D2z_aa;
    Double I_ESP_Kx3y3z_Fy2z_aa = I_ESP_Lx4y3z_D2z_aa+ABY*I_ESP_Kx3y3z_D2z_aa;
    Double I_ESP_Kx2y4z_Fy2z_aa = I_ESP_Lx3y4z_D2z_aa+ABY*I_ESP_Kx2y4z_D2z_aa;
    Double I_ESP_Kxy5z_Fy2z_aa = I_ESP_Lx2y5z_D2z_aa+ABY*I_ESP_Kxy5z_D2z_aa;
    Double I_ESP_Kx6z_Fy2z_aa = I_ESP_Lxy6z_D2z_aa+ABY*I_ESP_Kx6z_D2z_aa;
    Double I_ESP_K7y_Fy2z_aa = I_ESP_L8y_D2z_aa+ABY*I_ESP_K7y_D2z_aa;
    Double I_ESP_K6yz_Fy2z_aa = I_ESP_L7yz_D2z_aa+ABY*I_ESP_K6yz_D2z_aa;
    Double I_ESP_K5y2z_Fy2z_aa = I_ESP_L6y2z_D2z_aa+ABY*I_ESP_K5y2z_D2z_aa;
    Double I_ESP_K4y3z_Fy2z_aa = I_ESP_L5y3z_D2z_aa+ABY*I_ESP_K4y3z_D2z_aa;
    Double I_ESP_K3y4z_Fy2z_aa = I_ESP_L4y4z_D2z_aa+ABY*I_ESP_K3y4z_D2z_aa;
    Double I_ESP_K2y5z_Fy2z_aa = I_ESP_L3y5z_D2z_aa+ABY*I_ESP_K2y5z_D2z_aa;
    Double I_ESP_Ky6z_Fy2z_aa = I_ESP_L2y6z_D2z_aa+ABY*I_ESP_Ky6z_D2z_aa;
    Double I_ESP_K7z_Fy2z_aa = I_ESP_Ly7z_D2z_aa+ABY*I_ESP_K7z_D2z_aa;
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
     * shell quartet name: SQ_ESP_H_F_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_F_aa
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    abcd[iGrid*1260+0] = 4.0E0*I_ESP_K7x_F3x_aa-2.0E0*5*I_ESP_H5x_F3x_a-2.0E0*6*I_ESP_H5x_F3x_a+5*4*I_ESP_F3x_F3x;
    abcd[iGrid*1260+1] = 4.0E0*I_ESP_K6xy_F3x_aa-2.0E0*4*I_ESP_H4xy_F3x_a-2.0E0*5*I_ESP_H4xy_F3x_a+4*3*I_ESP_F2xy_F3x;
    abcd[iGrid*1260+2] = 4.0E0*I_ESP_K6xz_F3x_aa-2.0E0*4*I_ESP_H4xz_F3x_a-2.0E0*5*I_ESP_H4xz_F3x_a+4*3*I_ESP_F2xz_F3x;
    abcd[iGrid*1260+3] = 4.0E0*I_ESP_K5x2y_F3x_aa-2.0E0*3*I_ESP_H3x2y_F3x_a-2.0E0*4*I_ESP_H3x2y_F3x_a+3*2*I_ESP_Fx2y_F3x;
    abcd[iGrid*1260+4] = 4.0E0*I_ESP_K5xyz_F3x_aa-2.0E0*3*I_ESP_H3xyz_F3x_a-2.0E0*4*I_ESP_H3xyz_F3x_a+3*2*I_ESP_Fxyz_F3x;
    abcd[iGrid*1260+5] = 4.0E0*I_ESP_K5x2z_F3x_aa-2.0E0*3*I_ESP_H3x2z_F3x_a-2.0E0*4*I_ESP_H3x2z_F3x_a+3*2*I_ESP_Fx2z_F3x;
    abcd[iGrid*1260+6] = 4.0E0*I_ESP_K4x3y_F3x_aa-2.0E0*2*I_ESP_H2x3y_F3x_a-2.0E0*3*I_ESP_H2x3y_F3x_a+2*1*I_ESP_F3y_F3x;
    abcd[iGrid*1260+7] = 4.0E0*I_ESP_K4x2yz_F3x_aa-2.0E0*2*I_ESP_H2x2yz_F3x_a-2.0E0*3*I_ESP_H2x2yz_F3x_a+2*1*I_ESP_F2yz_F3x;
    abcd[iGrid*1260+8] = 4.0E0*I_ESP_K4xy2z_F3x_aa-2.0E0*2*I_ESP_H2xy2z_F3x_a-2.0E0*3*I_ESP_H2xy2z_F3x_a+2*1*I_ESP_Fy2z_F3x;
    abcd[iGrid*1260+9] = 4.0E0*I_ESP_K4x3z_F3x_aa-2.0E0*2*I_ESP_H2x3z_F3x_a-2.0E0*3*I_ESP_H2x3z_F3x_a+2*1*I_ESP_F3z_F3x;
    abcd[iGrid*1260+10] = 4.0E0*I_ESP_K3x4y_F3x_aa-2.0E0*1*I_ESP_Hx4y_F3x_a-2.0E0*2*I_ESP_Hx4y_F3x_a;
    abcd[iGrid*1260+11] = 4.0E0*I_ESP_K3x3yz_F3x_aa-2.0E0*1*I_ESP_Hx3yz_F3x_a-2.0E0*2*I_ESP_Hx3yz_F3x_a;
    abcd[iGrid*1260+12] = 4.0E0*I_ESP_K3x2y2z_F3x_aa-2.0E0*1*I_ESP_Hx2y2z_F3x_a-2.0E0*2*I_ESP_Hx2y2z_F3x_a;
    abcd[iGrid*1260+13] = 4.0E0*I_ESP_K3xy3z_F3x_aa-2.0E0*1*I_ESP_Hxy3z_F3x_a-2.0E0*2*I_ESP_Hxy3z_F3x_a;
    abcd[iGrid*1260+14] = 4.0E0*I_ESP_K3x4z_F3x_aa-2.0E0*1*I_ESP_Hx4z_F3x_a-2.0E0*2*I_ESP_Hx4z_F3x_a;
    abcd[iGrid*1260+15] = 4.0E0*I_ESP_K2x5y_F3x_aa-2.0E0*1*I_ESP_H5y_F3x_a;
    abcd[iGrid*1260+16] = 4.0E0*I_ESP_K2x4yz_F3x_aa-2.0E0*1*I_ESP_H4yz_F3x_a;
    abcd[iGrid*1260+17] = 4.0E0*I_ESP_K2x3y2z_F3x_aa-2.0E0*1*I_ESP_H3y2z_F3x_a;
    abcd[iGrid*1260+18] = 4.0E0*I_ESP_K2x2y3z_F3x_aa-2.0E0*1*I_ESP_H2y3z_F3x_a;
    abcd[iGrid*1260+19] = 4.0E0*I_ESP_K2xy4z_F3x_aa-2.0E0*1*I_ESP_Hy4z_F3x_a;
    abcd[iGrid*1260+20] = 4.0E0*I_ESP_K2x5z_F3x_aa-2.0E0*1*I_ESP_H5z_F3x_a;
    abcd[iGrid*1260+21] = 4.0E0*I_ESP_K7x_F2xy_aa-2.0E0*5*I_ESP_H5x_F2xy_a-2.0E0*6*I_ESP_H5x_F2xy_a+5*4*I_ESP_F3x_F2xy;
    abcd[iGrid*1260+22] = 4.0E0*I_ESP_K6xy_F2xy_aa-2.0E0*4*I_ESP_H4xy_F2xy_a-2.0E0*5*I_ESP_H4xy_F2xy_a+4*3*I_ESP_F2xy_F2xy;
    abcd[iGrid*1260+23] = 4.0E0*I_ESP_K6xz_F2xy_aa-2.0E0*4*I_ESP_H4xz_F2xy_a-2.0E0*5*I_ESP_H4xz_F2xy_a+4*3*I_ESP_F2xz_F2xy;
    abcd[iGrid*1260+24] = 4.0E0*I_ESP_K5x2y_F2xy_aa-2.0E0*3*I_ESP_H3x2y_F2xy_a-2.0E0*4*I_ESP_H3x2y_F2xy_a+3*2*I_ESP_Fx2y_F2xy;
    abcd[iGrid*1260+25] = 4.0E0*I_ESP_K5xyz_F2xy_aa-2.0E0*3*I_ESP_H3xyz_F2xy_a-2.0E0*4*I_ESP_H3xyz_F2xy_a+3*2*I_ESP_Fxyz_F2xy;
    abcd[iGrid*1260+26] = 4.0E0*I_ESP_K5x2z_F2xy_aa-2.0E0*3*I_ESP_H3x2z_F2xy_a-2.0E0*4*I_ESP_H3x2z_F2xy_a+3*2*I_ESP_Fx2z_F2xy;
    abcd[iGrid*1260+27] = 4.0E0*I_ESP_K4x3y_F2xy_aa-2.0E0*2*I_ESP_H2x3y_F2xy_a-2.0E0*3*I_ESP_H2x3y_F2xy_a+2*1*I_ESP_F3y_F2xy;
    abcd[iGrid*1260+28] = 4.0E0*I_ESP_K4x2yz_F2xy_aa-2.0E0*2*I_ESP_H2x2yz_F2xy_a-2.0E0*3*I_ESP_H2x2yz_F2xy_a+2*1*I_ESP_F2yz_F2xy;
    abcd[iGrid*1260+29] = 4.0E0*I_ESP_K4xy2z_F2xy_aa-2.0E0*2*I_ESP_H2xy2z_F2xy_a-2.0E0*3*I_ESP_H2xy2z_F2xy_a+2*1*I_ESP_Fy2z_F2xy;
    abcd[iGrid*1260+30] = 4.0E0*I_ESP_K4x3z_F2xy_aa-2.0E0*2*I_ESP_H2x3z_F2xy_a-2.0E0*3*I_ESP_H2x3z_F2xy_a+2*1*I_ESP_F3z_F2xy;
    abcd[iGrid*1260+31] = 4.0E0*I_ESP_K3x4y_F2xy_aa-2.0E0*1*I_ESP_Hx4y_F2xy_a-2.0E0*2*I_ESP_Hx4y_F2xy_a;
    abcd[iGrid*1260+32] = 4.0E0*I_ESP_K3x3yz_F2xy_aa-2.0E0*1*I_ESP_Hx3yz_F2xy_a-2.0E0*2*I_ESP_Hx3yz_F2xy_a;
    abcd[iGrid*1260+33] = 4.0E0*I_ESP_K3x2y2z_F2xy_aa-2.0E0*1*I_ESP_Hx2y2z_F2xy_a-2.0E0*2*I_ESP_Hx2y2z_F2xy_a;
    abcd[iGrid*1260+34] = 4.0E0*I_ESP_K3xy3z_F2xy_aa-2.0E0*1*I_ESP_Hxy3z_F2xy_a-2.0E0*2*I_ESP_Hxy3z_F2xy_a;
    abcd[iGrid*1260+35] = 4.0E0*I_ESP_K3x4z_F2xy_aa-2.0E0*1*I_ESP_Hx4z_F2xy_a-2.0E0*2*I_ESP_Hx4z_F2xy_a;
    abcd[iGrid*1260+36] = 4.0E0*I_ESP_K2x5y_F2xy_aa-2.0E0*1*I_ESP_H5y_F2xy_a;
    abcd[iGrid*1260+37] = 4.0E0*I_ESP_K2x4yz_F2xy_aa-2.0E0*1*I_ESP_H4yz_F2xy_a;
    abcd[iGrid*1260+38] = 4.0E0*I_ESP_K2x3y2z_F2xy_aa-2.0E0*1*I_ESP_H3y2z_F2xy_a;
    abcd[iGrid*1260+39] = 4.0E0*I_ESP_K2x2y3z_F2xy_aa-2.0E0*1*I_ESP_H2y3z_F2xy_a;
    abcd[iGrid*1260+40] = 4.0E0*I_ESP_K2xy4z_F2xy_aa-2.0E0*1*I_ESP_Hy4z_F2xy_a;
    abcd[iGrid*1260+41] = 4.0E0*I_ESP_K2x5z_F2xy_aa-2.0E0*1*I_ESP_H5z_F2xy_a;
    abcd[iGrid*1260+42] = 4.0E0*I_ESP_K7x_F2xz_aa-2.0E0*5*I_ESP_H5x_F2xz_a-2.0E0*6*I_ESP_H5x_F2xz_a+5*4*I_ESP_F3x_F2xz;
    abcd[iGrid*1260+43] = 4.0E0*I_ESP_K6xy_F2xz_aa-2.0E0*4*I_ESP_H4xy_F2xz_a-2.0E0*5*I_ESP_H4xy_F2xz_a+4*3*I_ESP_F2xy_F2xz;
    abcd[iGrid*1260+44] = 4.0E0*I_ESP_K6xz_F2xz_aa-2.0E0*4*I_ESP_H4xz_F2xz_a-2.0E0*5*I_ESP_H4xz_F2xz_a+4*3*I_ESP_F2xz_F2xz;
    abcd[iGrid*1260+45] = 4.0E0*I_ESP_K5x2y_F2xz_aa-2.0E0*3*I_ESP_H3x2y_F2xz_a-2.0E0*4*I_ESP_H3x2y_F2xz_a+3*2*I_ESP_Fx2y_F2xz;
    abcd[iGrid*1260+46] = 4.0E0*I_ESP_K5xyz_F2xz_aa-2.0E0*3*I_ESP_H3xyz_F2xz_a-2.0E0*4*I_ESP_H3xyz_F2xz_a+3*2*I_ESP_Fxyz_F2xz;
    abcd[iGrid*1260+47] = 4.0E0*I_ESP_K5x2z_F2xz_aa-2.0E0*3*I_ESP_H3x2z_F2xz_a-2.0E0*4*I_ESP_H3x2z_F2xz_a+3*2*I_ESP_Fx2z_F2xz;
    abcd[iGrid*1260+48] = 4.0E0*I_ESP_K4x3y_F2xz_aa-2.0E0*2*I_ESP_H2x3y_F2xz_a-2.0E0*3*I_ESP_H2x3y_F2xz_a+2*1*I_ESP_F3y_F2xz;
    abcd[iGrid*1260+49] = 4.0E0*I_ESP_K4x2yz_F2xz_aa-2.0E0*2*I_ESP_H2x2yz_F2xz_a-2.0E0*3*I_ESP_H2x2yz_F2xz_a+2*1*I_ESP_F2yz_F2xz;
    abcd[iGrid*1260+50] = 4.0E0*I_ESP_K4xy2z_F2xz_aa-2.0E0*2*I_ESP_H2xy2z_F2xz_a-2.0E0*3*I_ESP_H2xy2z_F2xz_a+2*1*I_ESP_Fy2z_F2xz;
    abcd[iGrid*1260+51] = 4.0E0*I_ESP_K4x3z_F2xz_aa-2.0E0*2*I_ESP_H2x3z_F2xz_a-2.0E0*3*I_ESP_H2x3z_F2xz_a+2*1*I_ESP_F3z_F2xz;
    abcd[iGrid*1260+52] = 4.0E0*I_ESP_K3x4y_F2xz_aa-2.0E0*1*I_ESP_Hx4y_F2xz_a-2.0E0*2*I_ESP_Hx4y_F2xz_a;
    abcd[iGrid*1260+53] = 4.0E0*I_ESP_K3x3yz_F2xz_aa-2.0E0*1*I_ESP_Hx3yz_F2xz_a-2.0E0*2*I_ESP_Hx3yz_F2xz_a;
    abcd[iGrid*1260+54] = 4.0E0*I_ESP_K3x2y2z_F2xz_aa-2.0E0*1*I_ESP_Hx2y2z_F2xz_a-2.0E0*2*I_ESP_Hx2y2z_F2xz_a;
    abcd[iGrid*1260+55] = 4.0E0*I_ESP_K3xy3z_F2xz_aa-2.0E0*1*I_ESP_Hxy3z_F2xz_a-2.0E0*2*I_ESP_Hxy3z_F2xz_a;
    abcd[iGrid*1260+56] = 4.0E0*I_ESP_K3x4z_F2xz_aa-2.0E0*1*I_ESP_Hx4z_F2xz_a-2.0E0*2*I_ESP_Hx4z_F2xz_a;
    abcd[iGrid*1260+57] = 4.0E0*I_ESP_K2x5y_F2xz_aa-2.0E0*1*I_ESP_H5y_F2xz_a;
    abcd[iGrid*1260+58] = 4.0E0*I_ESP_K2x4yz_F2xz_aa-2.0E0*1*I_ESP_H4yz_F2xz_a;
    abcd[iGrid*1260+59] = 4.0E0*I_ESP_K2x3y2z_F2xz_aa-2.0E0*1*I_ESP_H3y2z_F2xz_a;
    abcd[iGrid*1260+60] = 4.0E0*I_ESP_K2x2y3z_F2xz_aa-2.0E0*1*I_ESP_H2y3z_F2xz_a;
    abcd[iGrid*1260+61] = 4.0E0*I_ESP_K2xy4z_F2xz_aa-2.0E0*1*I_ESP_Hy4z_F2xz_a;
    abcd[iGrid*1260+62] = 4.0E0*I_ESP_K2x5z_F2xz_aa-2.0E0*1*I_ESP_H5z_F2xz_a;
    abcd[iGrid*1260+63] = 4.0E0*I_ESP_K7x_Fx2y_aa-2.0E0*5*I_ESP_H5x_Fx2y_a-2.0E0*6*I_ESP_H5x_Fx2y_a+5*4*I_ESP_F3x_Fx2y;
    abcd[iGrid*1260+64] = 4.0E0*I_ESP_K6xy_Fx2y_aa-2.0E0*4*I_ESP_H4xy_Fx2y_a-2.0E0*5*I_ESP_H4xy_Fx2y_a+4*3*I_ESP_F2xy_Fx2y;
    abcd[iGrid*1260+65] = 4.0E0*I_ESP_K6xz_Fx2y_aa-2.0E0*4*I_ESP_H4xz_Fx2y_a-2.0E0*5*I_ESP_H4xz_Fx2y_a+4*3*I_ESP_F2xz_Fx2y;
    abcd[iGrid*1260+66] = 4.0E0*I_ESP_K5x2y_Fx2y_aa-2.0E0*3*I_ESP_H3x2y_Fx2y_a-2.0E0*4*I_ESP_H3x2y_Fx2y_a+3*2*I_ESP_Fx2y_Fx2y;
    abcd[iGrid*1260+67] = 4.0E0*I_ESP_K5xyz_Fx2y_aa-2.0E0*3*I_ESP_H3xyz_Fx2y_a-2.0E0*4*I_ESP_H3xyz_Fx2y_a+3*2*I_ESP_Fxyz_Fx2y;
    abcd[iGrid*1260+68] = 4.0E0*I_ESP_K5x2z_Fx2y_aa-2.0E0*3*I_ESP_H3x2z_Fx2y_a-2.0E0*4*I_ESP_H3x2z_Fx2y_a+3*2*I_ESP_Fx2z_Fx2y;
    abcd[iGrid*1260+69] = 4.0E0*I_ESP_K4x3y_Fx2y_aa-2.0E0*2*I_ESP_H2x3y_Fx2y_a-2.0E0*3*I_ESP_H2x3y_Fx2y_a+2*1*I_ESP_F3y_Fx2y;
    abcd[iGrid*1260+70] = 4.0E0*I_ESP_K4x2yz_Fx2y_aa-2.0E0*2*I_ESP_H2x2yz_Fx2y_a-2.0E0*3*I_ESP_H2x2yz_Fx2y_a+2*1*I_ESP_F2yz_Fx2y;
    abcd[iGrid*1260+71] = 4.0E0*I_ESP_K4xy2z_Fx2y_aa-2.0E0*2*I_ESP_H2xy2z_Fx2y_a-2.0E0*3*I_ESP_H2xy2z_Fx2y_a+2*1*I_ESP_Fy2z_Fx2y;
    abcd[iGrid*1260+72] = 4.0E0*I_ESP_K4x3z_Fx2y_aa-2.0E0*2*I_ESP_H2x3z_Fx2y_a-2.0E0*3*I_ESP_H2x3z_Fx2y_a+2*1*I_ESP_F3z_Fx2y;
    abcd[iGrid*1260+73] = 4.0E0*I_ESP_K3x4y_Fx2y_aa-2.0E0*1*I_ESP_Hx4y_Fx2y_a-2.0E0*2*I_ESP_Hx4y_Fx2y_a;
    abcd[iGrid*1260+74] = 4.0E0*I_ESP_K3x3yz_Fx2y_aa-2.0E0*1*I_ESP_Hx3yz_Fx2y_a-2.0E0*2*I_ESP_Hx3yz_Fx2y_a;
    abcd[iGrid*1260+75] = 4.0E0*I_ESP_K3x2y2z_Fx2y_aa-2.0E0*1*I_ESP_Hx2y2z_Fx2y_a-2.0E0*2*I_ESP_Hx2y2z_Fx2y_a;
    abcd[iGrid*1260+76] = 4.0E0*I_ESP_K3xy3z_Fx2y_aa-2.0E0*1*I_ESP_Hxy3z_Fx2y_a-2.0E0*2*I_ESP_Hxy3z_Fx2y_a;
    abcd[iGrid*1260+77] = 4.0E0*I_ESP_K3x4z_Fx2y_aa-2.0E0*1*I_ESP_Hx4z_Fx2y_a-2.0E0*2*I_ESP_Hx4z_Fx2y_a;
    abcd[iGrid*1260+78] = 4.0E0*I_ESP_K2x5y_Fx2y_aa-2.0E0*1*I_ESP_H5y_Fx2y_a;
    abcd[iGrid*1260+79] = 4.0E0*I_ESP_K2x4yz_Fx2y_aa-2.0E0*1*I_ESP_H4yz_Fx2y_a;
    abcd[iGrid*1260+80] = 4.0E0*I_ESP_K2x3y2z_Fx2y_aa-2.0E0*1*I_ESP_H3y2z_Fx2y_a;
    abcd[iGrid*1260+81] = 4.0E0*I_ESP_K2x2y3z_Fx2y_aa-2.0E0*1*I_ESP_H2y3z_Fx2y_a;
    abcd[iGrid*1260+82] = 4.0E0*I_ESP_K2xy4z_Fx2y_aa-2.0E0*1*I_ESP_Hy4z_Fx2y_a;
    abcd[iGrid*1260+83] = 4.0E0*I_ESP_K2x5z_Fx2y_aa-2.0E0*1*I_ESP_H5z_Fx2y_a;
    abcd[iGrid*1260+84] = 4.0E0*I_ESP_K7x_Fxyz_aa-2.0E0*5*I_ESP_H5x_Fxyz_a-2.0E0*6*I_ESP_H5x_Fxyz_a+5*4*I_ESP_F3x_Fxyz;
    abcd[iGrid*1260+85] = 4.0E0*I_ESP_K6xy_Fxyz_aa-2.0E0*4*I_ESP_H4xy_Fxyz_a-2.0E0*5*I_ESP_H4xy_Fxyz_a+4*3*I_ESP_F2xy_Fxyz;
    abcd[iGrid*1260+86] = 4.0E0*I_ESP_K6xz_Fxyz_aa-2.0E0*4*I_ESP_H4xz_Fxyz_a-2.0E0*5*I_ESP_H4xz_Fxyz_a+4*3*I_ESP_F2xz_Fxyz;
    abcd[iGrid*1260+87] = 4.0E0*I_ESP_K5x2y_Fxyz_aa-2.0E0*3*I_ESP_H3x2y_Fxyz_a-2.0E0*4*I_ESP_H3x2y_Fxyz_a+3*2*I_ESP_Fx2y_Fxyz;
    abcd[iGrid*1260+88] = 4.0E0*I_ESP_K5xyz_Fxyz_aa-2.0E0*3*I_ESP_H3xyz_Fxyz_a-2.0E0*4*I_ESP_H3xyz_Fxyz_a+3*2*I_ESP_Fxyz_Fxyz;
    abcd[iGrid*1260+89] = 4.0E0*I_ESP_K5x2z_Fxyz_aa-2.0E0*3*I_ESP_H3x2z_Fxyz_a-2.0E0*4*I_ESP_H3x2z_Fxyz_a+3*2*I_ESP_Fx2z_Fxyz;
    abcd[iGrid*1260+90] = 4.0E0*I_ESP_K4x3y_Fxyz_aa-2.0E0*2*I_ESP_H2x3y_Fxyz_a-2.0E0*3*I_ESP_H2x3y_Fxyz_a+2*1*I_ESP_F3y_Fxyz;
    abcd[iGrid*1260+91] = 4.0E0*I_ESP_K4x2yz_Fxyz_aa-2.0E0*2*I_ESP_H2x2yz_Fxyz_a-2.0E0*3*I_ESP_H2x2yz_Fxyz_a+2*1*I_ESP_F2yz_Fxyz;
    abcd[iGrid*1260+92] = 4.0E0*I_ESP_K4xy2z_Fxyz_aa-2.0E0*2*I_ESP_H2xy2z_Fxyz_a-2.0E0*3*I_ESP_H2xy2z_Fxyz_a+2*1*I_ESP_Fy2z_Fxyz;
    abcd[iGrid*1260+93] = 4.0E0*I_ESP_K4x3z_Fxyz_aa-2.0E0*2*I_ESP_H2x3z_Fxyz_a-2.0E0*3*I_ESP_H2x3z_Fxyz_a+2*1*I_ESP_F3z_Fxyz;
    abcd[iGrid*1260+94] = 4.0E0*I_ESP_K3x4y_Fxyz_aa-2.0E0*1*I_ESP_Hx4y_Fxyz_a-2.0E0*2*I_ESP_Hx4y_Fxyz_a;
    abcd[iGrid*1260+95] = 4.0E0*I_ESP_K3x3yz_Fxyz_aa-2.0E0*1*I_ESP_Hx3yz_Fxyz_a-2.0E0*2*I_ESP_Hx3yz_Fxyz_a;
    abcd[iGrid*1260+96] = 4.0E0*I_ESP_K3x2y2z_Fxyz_aa-2.0E0*1*I_ESP_Hx2y2z_Fxyz_a-2.0E0*2*I_ESP_Hx2y2z_Fxyz_a;
    abcd[iGrid*1260+97] = 4.0E0*I_ESP_K3xy3z_Fxyz_aa-2.0E0*1*I_ESP_Hxy3z_Fxyz_a-2.0E0*2*I_ESP_Hxy3z_Fxyz_a;
    abcd[iGrid*1260+98] = 4.0E0*I_ESP_K3x4z_Fxyz_aa-2.0E0*1*I_ESP_Hx4z_Fxyz_a-2.0E0*2*I_ESP_Hx4z_Fxyz_a;
    abcd[iGrid*1260+99] = 4.0E0*I_ESP_K2x5y_Fxyz_aa-2.0E0*1*I_ESP_H5y_Fxyz_a;
    abcd[iGrid*1260+100] = 4.0E0*I_ESP_K2x4yz_Fxyz_aa-2.0E0*1*I_ESP_H4yz_Fxyz_a;
    abcd[iGrid*1260+101] = 4.0E0*I_ESP_K2x3y2z_Fxyz_aa-2.0E0*1*I_ESP_H3y2z_Fxyz_a;
    abcd[iGrid*1260+102] = 4.0E0*I_ESP_K2x2y3z_Fxyz_aa-2.0E0*1*I_ESP_H2y3z_Fxyz_a;
    abcd[iGrid*1260+103] = 4.0E0*I_ESP_K2xy4z_Fxyz_aa-2.0E0*1*I_ESP_Hy4z_Fxyz_a;
    abcd[iGrid*1260+104] = 4.0E0*I_ESP_K2x5z_Fxyz_aa-2.0E0*1*I_ESP_H5z_Fxyz_a;
    abcd[iGrid*1260+105] = 4.0E0*I_ESP_K7x_Fx2z_aa-2.0E0*5*I_ESP_H5x_Fx2z_a-2.0E0*6*I_ESP_H5x_Fx2z_a+5*4*I_ESP_F3x_Fx2z;
    abcd[iGrid*1260+106] = 4.0E0*I_ESP_K6xy_Fx2z_aa-2.0E0*4*I_ESP_H4xy_Fx2z_a-2.0E0*5*I_ESP_H4xy_Fx2z_a+4*3*I_ESP_F2xy_Fx2z;
    abcd[iGrid*1260+107] = 4.0E0*I_ESP_K6xz_Fx2z_aa-2.0E0*4*I_ESP_H4xz_Fx2z_a-2.0E0*5*I_ESP_H4xz_Fx2z_a+4*3*I_ESP_F2xz_Fx2z;
    abcd[iGrid*1260+108] = 4.0E0*I_ESP_K5x2y_Fx2z_aa-2.0E0*3*I_ESP_H3x2y_Fx2z_a-2.0E0*4*I_ESP_H3x2y_Fx2z_a+3*2*I_ESP_Fx2y_Fx2z;
    abcd[iGrid*1260+109] = 4.0E0*I_ESP_K5xyz_Fx2z_aa-2.0E0*3*I_ESP_H3xyz_Fx2z_a-2.0E0*4*I_ESP_H3xyz_Fx2z_a+3*2*I_ESP_Fxyz_Fx2z;
    abcd[iGrid*1260+110] = 4.0E0*I_ESP_K5x2z_Fx2z_aa-2.0E0*3*I_ESP_H3x2z_Fx2z_a-2.0E0*4*I_ESP_H3x2z_Fx2z_a+3*2*I_ESP_Fx2z_Fx2z;
    abcd[iGrid*1260+111] = 4.0E0*I_ESP_K4x3y_Fx2z_aa-2.0E0*2*I_ESP_H2x3y_Fx2z_a-2.0E0*3*I_ESP_H2x3y_Fx2z_a+2*1*I_ESP_F3y_Fx2z;
    abcd[iGrid*1260+112] = 4.0E0*I_ESP_K4x2yz_Fx2z_aa-2.0E0*2*I_ESP_H2x2yz_Fx2z_a-2.0E0*3*I_ESP_H2x2yz_Fx2z_a+2*1*I_ESP_F2yz_Fx2z;
    abcd[iGrid*1260+113] = 4.0E0*I_ESP_K4xy2z_Fx2z_aa-2.0E0*2*I_ESP_H2xy2z_Fx2z_a-2.0E0*3*I_ESP_H2xy2z_Fx2z_a+2*1*I_ESP_Fy2z_Fx2z;
    abcd[iGrid*1260+114] = 4.0E0*I_ESP_K4x3z_Fx2z_aa-2.0E0*2*I_ESP_H2x3z_Fx2z_a-2.0E0*3*I_ESP_H2x3z_Fx2z_a+2*1*I_ESP_F3z_Fx2z;
    abcd[iGrid*1260+115] = 4.0E0*I_ESP_K3x4y_Fx2z_aa-2.0E0*1*I_ESP_Hx4y_Fx2z_a-2.0E0*2*I_ESP_Hx4y_Fx2z_a;
    abcd[iGrid*1260+116] = 4.0E0*I_ESP_K3x3yz_Fx2z_aa-2.0E0*1*I_ESP_Hx3yz_Fx2z_a-2.0E0*2*I_ESP_Hx3yz_Fx2z_a;
    abcd[iGrid*1260+117] = 4.0E0*I_ESP_K3x2y2z_Fx2z_aa-2.0E0*1*I_ESP_Hx2y2z_Fx2z_a-2.0E0*2*I_ESP_Hx2y2z_Fx2z_a;
    abcd[iGrid*1260+118] = 4.0E0*I_ESP_K3xy3z_Fx2z_aa-2.0E0*1*I_ESP_Hxy3z_Fx2z_a-2.0E0*2*I_ESP_Hxy3z_Fx2z_a;
    abcd[iGrid*1260+119] = 4.0E0*I_ESP_K3x4z_Fx2z_aa-2.0E0*1*I_ESP_Hx4z_Fx2z_a-2.0E0*2*I_ESP_Hx4z_Fx2z_a;
    abcd[iGrid*1260+120] = 4.0E0*I_ESP_K2x5y_Fx2z_aa-2.0E0*1*I_ESP_H5y_Fx2z_a;
    abcd[iGrid*1260+121] = 4.0E0*I_ESP_K2x4yz_Fx2z_aa-2.0E0*1*I_ESP_H4yz_Fx2z_a;
    abcd[iGrid*1260+122] = 4.0E0*I_ESP_K2x3y2z_Fx2z_aa-2.0E0*1*I_ESP_H3y2z_Fx2z_a;
    abcd[iGrid*1260+123] = 4.0E0*I_ESP_K2x2y3z_Fx2z_aa-2.0E0*1*I_ESP_H2y3z_Fx2z_a;
    abcd[iGrid*1260+124] = 4.0E0*I_ESP_K2xy4z_Fx2z_aa-2.0E0*1*I_ESP_Hy4z_Fx2z_a;
    abcd[iGrid*1260+125] = 4.0E0*I_ESP_K2x5z_Fx2z_aa-2.0E0*1*I_ESP_H5z_Fx2z_a;
    abcd[iGrid*1260+126] = 4.0E0*I_ESP_K7x_F3y_aa-2.0E0*5*I_ESP_H5x_F3y_a-2.0E0*6*I_ESP_H5x_F3y_a+5*4*I_ESP_F3x_F3y;
    abcd[iGrid*1260+127] = 4.0E0*I_ESP_K6xy_F3y_aa-2.0E0*4*I_ESP_H4xy_F3y_a-2.0E0*5*I_ESP_H4xy_F3y_a+4*3*I_ESP_F2xy_F3y;
    abcd[iGrid*1260+128] = 4.0E0*I_ESP_K6xz_F3y_aa-2.0E0*4*I_ESP_H4xz_F3y_a-2.0E0*5*I_ESP_H4xz_F3y_a+4*3*I_ESP_F2xz_F3y;
    abcd[iGrid*1260+129] = 4.0E0*I_ESP_K5x2y_F3y_aa-2.0E0*3*I_ESP_H3x2y_F3y_a-2.0E0*4*I_ESP_H3x2y_F3y_a+3*2*I_ESP_Fx2y_F3y;
    abcd[iGrid*1260+130] = 4.0E0*I_ESP_K5xyz_F3y_aa-2.0E0*3*I_ESP_H3xyz_F3y_a-2.0E0*4*I_ESP_H3xyz_F3y_a+3*2*I_ESP_Fxyz_F3y;
    abcd[iGrid*1260+131] = 4.0E0*I_ESP_K5x2z_F3y_aa-2.0E0*3*I_ESP_H3x2z_F3y_a-2.0E0*4*I_ESP_H3x2z_F3y_a+3*2*I_ESP_Fx2z_F3y;
    abcd[iGrid*1260+132] = 4.0E0*I_ESP_K4x3y_F3y_aa-2.0E0*2*I_ESP_H2x3y_F3y_a-2.0E0*3*I_ESP_H2x3y_F3y_a+2*1*I_ESP_F3y_F3y;
    abcd[iGrid*1260+133] = 4.0E0*I_ESP_K4x2yz_F3y_aa-2.0E0*2*I_ESP_H2x2yz_F3y_a-2.0E0*3*I_ESP_H2x2yz_F3y_a+2*1*I_ESP_F2yz_F3y;
    abcd[iGrid*1260+134] = 4.0E0*I_ESP_K4xy2z_F3y_aa-2.0E0*2*I_ESP_H2xy2z_F3y_a-2.0E0*3*I_ESP_H2xy2z_F3y_a+2*1*I_ESP_Fy2z_F3y;
    abcd[iGrid*1260+135] = 4.0E0*I_ESP_K4x3z_F3y_aa-2.0E0*2*I_ESP_H2x3z_F3y_a-2.0E0*3*I_ESP_H2x3z_F3y_a+2*1*I_ESP_F3z_F3y;
    abcd[iGrid*1260+136] = 4.0E0*I_ESP_K3x4y_F3y_aa-2.0E0*1*I_ESP_Hx4y_F3y_a-2.0E0*2*I_ESP_Hx4y_F3y_a;
    abcd[iGrid*1260+137] = 4.0E0*I_ESP_K3x3yz_F3y_aa-2.0E0*1*I_ESP_Hx3yz_F3y_a-2.0E0*2*I_ESP_Hx3yz_F3y_a;
    abcd[iGrid*1260+138] = 4.0E0*I_ESP_K3x2y2z_F3y_aa-2.0E0*1*I_ESP_Hx2y2z_F3y_a-2.0E0*2*I_ESP_Hx2y2z_F3y_a;
    abcd[iGrid*1260+139] = 4.0E0*I_ESP_K3xy3z_F3y_aa-2.0E0*1*I_ESP_Hxy3z_F3y_a-2.0E0*2*I_ESP_Hxy3z_F3y_a;
    abcd[iGrid*1260+140] = 4.0E0*I_ESP_K3x4z_F3y_aa-2.0E0*1*I_ESP_Hx4z_F3y_a-2.0E0*2*I_ESP_Hx4z_F3y_a;
    abcd[iGrid*1260+141] = 4.0E0*I_ESP_K2x5y_F3y_aa-2.0E0*1*I_ESP_H5y_F3y_a;
    abcd[iGrid*1260+142] = 4.0E0*I_ESP_K2x4yz_F3y_aa-2.0E0*1*I_ESP_H4yz_F3y_a;
    abcd[iGrid*1260+143] = 4.0E0*I_ESP_K2x3y2z_F3y_aa-2.0E0*1*I_ESP_H3y2z_F3y_a;
    abcd[iGrid*1260+144] = 4.0E0*I_ESP_K2x2y3z_F3y_aa-2.0E0*1*I_ESP_H2y3z_F3y_a;
    abcd[iGrid*1260+145] = 4.0E0*I_ESP_K2xy4z_F3y_aa-2.0E0*1*I_ESP_Hy4z_F3y_a;
    abcd[iGrid*1260+146] = 4.0E0*I_ESP_K2x5z_F3y_aa-2.0E0*1*I_ESP_H5z_F3y_a;
    abcd[iGrid*1260+147] = 4.0E0*I_ESP_K7x_F2yz_aa-2.0E0*5*I_ESP_H5x_F2yz_a-2.0E0*6*I_ESP_H5x_F2yz_a+5*4*I_ESP_F3x_F2yz;
    abcd[iGrid*1260+148] = 4.0E0*I_ESP_K6xy_F2yz_aa-2.0E0*4*I_ESP_H4xy_F2yz_a-2.0E0*5*I_ESP_H4xy_F2yz_a+4*3*I_ESP_F2xy_F2yz;
    abcd[iGrid*1260+149] = 4.0E0*I_ESP_K6xz_F2yz_aa-2.0E0*4*I_ESP_H4xz_F2yz_a-2.0E0*5*I_ESP_H4xz_F2yz_a+4*3*I_ESP_F2xz_F2yz;
    abcd[iGrid*1260+150] = 4.0E0*I_ESP_K5x2y_F2yz_aa-2.0E0*3*I_ESP_H3x2y_F2yz_a-2.0E0*4*I_ESP_H3x2y_F2yz_a+3*2*I_ESP_Fx2y_F2yz;
    abcd[iGrid*1260+151] = 4.0E0*I_ESP_K5xyz_F2yz_aa-2.0E0*3*I_ESP_H3xyz_F2yz_a-2.0E0*4*I_ESP_H3xyz_F2yz_a+3*2*I_ESP_Fxyz_F2yz;
    abcd[iGrid*1260+152] = 4.0E0*I_ESP_K5x2z_F2yz_aa-2.0E0*3*I_ESP_H3x2z_F2yz_a-2.0E0*4*I_ESP_H3x2z_F2yz_a+3*2*I_ESP_Fx2z_F2yz;
    abcd[iGrid*1260+153] = 4.0E0*I_ESP_K4x3y_F2yz_aa-2.0E0*2*I_ESP_H2x3y_F2yz_a-2.0E0*3*I_ESP_H2x3y_F2yz_a+2*1*I_ESP_F3y_F2yz;
    abcd[iGrid*1260+154] = 4.0E0*I_ESP_K4x2yz_F2yz_aa-2.0E0*2*I_ESP_H2x2yz_F2yz_a-2.0E0*3*I_ESP_H2x2yz_F2yz_a+2*1*I_ESP_F2yz_F2yz;
    abcd[iGrid*1260+155] = 4.0E0*I_ESP_K4xy2z_F2yz_aa-2.0E0*2*I_ESP_H2xy2z_F2yz_a-2.0E0*3*I_ESP_H2xy2z_F2yz_a+2*1*I_ESP_Fy2z_F2yz;
    abcd[iGrid*1260+156] = 4.0E0*I_ESP_K4x3z_F2yz_aa-2.0E0*2*I_ESP_H2x3z_F2yz_a-2.0E0*3*I_ESP_H2x3z_F2yz_a+2*1*I_ESP_F3z_F2yz;
    abcd[iGrid*1260+157] = 4.0E0*I_ESP_K3x4y_F2yz_aa-2.0E0*1*I_ESP_Hx4y_F2yz_a-2.0E0*2*I_ESP_Hx4y_F2yz_a;
    abcd[iGrid*1260+158] = 4.0E0*I_ESP_K3x3yz_F2yz_aa-2.0E0*1*I_ESP_Hx3yz_F2yz_a-2.0E0*2*I_ESP_Hx3yz_F2yz_a;
    abcd[iGrid*1260+159] = 4.0E0*I_ESP_K3x2y2z_F2yz_aa-2.0E0*1*I_ESP_Hx2y2z_F2yz_a-2.0E0*2*I_ESP_Hx2y2z_F2yz_a;
    abcd[iGrid*1260+160] = 4.0E0*I_ESP_K3xy3z_F2yz_aa-2.0E0*1*I_ESP_Hxy3z_F2yz_a-2.0E0*2*I_ESP_Hxy3z_F2yz_a;
    abcd[iGrid*1260+161] = 4.0E0*I_ESP_K3x4z_F2yz_aa-2.0E0*1*I_ESP_Hx4z_F2yz_a-2.0E0*2*I_ESP_Hx4z_F2yz_a;
    abcd[iGrid*1260+162] = 4.0E0*I_ESP_K2x5y_F2yz_aa-2.0E0*1*I_ESP_H5y_F2yz_a;
    abcd[iGrid*1260+163] = 4.0E0*I_ESP_K2x4yz_F2yz_aa-2.0E0*1*I_ESP_H4yz_F2yz_a;
    abcd[iGrid*1260+164] = 4.0E0*I_ESP_K2x3y2z_F2yz_aa-2.0E0*1*I_ESP_H3y2z_F2yz_a;
    abcd[iGrid*1260+165] = 4.0E0*I_ESP_K2x2y3z_F2yz_aa-2.0E0*1*I_ESP_H2y3z_F2yz_a;
    abcd[iGrid*1260+166] = 4.0E0*I_ESP_K2xy4z_F2yz_aa-2.0E0*1*I_ESP_Hy4z_F2yz_a;
    abcd[iGrid*1260+167] = 4.0E0*I_ESP_K2x5z_F2yz_aa-2.0E0*1*I_ESP_H5z_F2yz_a;
    abcd[iGrid*1260+168] = 4.0E0*I_ESP_K7x_Fy2z_aa-2.0E0*5*I_ESP_H5x_Fy2z_a-2.0E0*6*I_ESP_H5x_Fy2z_a+5*4*I_ESP_F3x_Fy2z;
    abcd[iGrid*1260+169] = 4.0E0*I_ESP_K6xy_Fy2z_aa-2.0E0*4*I_ESP_H4xy_Fy2z_a-2.0E0*5*I_ESP_H4xy_Fy2z_a+4*3*I_ESP_F2xy_Fy2z;
    abcd[iGrid*1260+170] = 4.0E0*I_ESP_K6xz_Fy2z_aa-2.0E0*4*I_ESP_H4xz_Fy2z_a-2.0E0*5*I_ESP_H4xz_Fy2z_a+4*3*I_ESP_F2xz_Fy2z;
    abcd[iGrid*1260+171] = 4.0E0*I_ESP_K5x2y_Fy2z_aa-2.0E0*3*I_ESP_H3x2y_Fy2z_a-2.0E0*4*I_ESP_H3x2y_Fy2z_a+3*2*I_ESP_Fx2y_Fy2z;
    abcd[iGrid*1260+172] = 4.0E0*I_ESP_K5xyz_Fy2z_aa-2.0E0*3*I_ESP_H3xyz_Fy2z_a-2.0E0*4*I_ESP_H3xyz_Fy2z_a+3*2*I_ESP_Fxyz_Fy2z;
    abcd[iGrid*1260+173] = 4.0E0*I_ESP_K5x2z_Fy2z_aa-2.0E0*3*I_ESP_H3x2z_Fy2z_a-2.0E0*4*I_ESP_H3x2z_Fy2z_a+3*2*I_ESP_Fx2z_Fy2z;
    abcd[iGrid*1260+174] = 4.0E0*I_ESP_K4x3y_Fy2z_aa-2.0E0*2*I_ESP_H2x3y_Fy2z_a-2.0E0*3*I_ESP_H2x3y_Fy2z_a+2*1*I_ESP_F3y_Fy2z;
    abcd[iGrid*1260+175] = 4.0E0*I_ESP_K4x2yz_Fy2z_aa-2.0E0*2*I_ESP_H2x2yz_Fy2z_a-2.0E0*3*I_ESP_H2x2yz_Fy2z_a+2*1*I_ESP_F2yz_Fy2z;
    abcd[iGrid*1260+176] = 4.0E0*I_ESP_K4xy2z_Fy2z_aa-2.0E0*2*I_ESP_H2xy2z_Fy2z_a-2.0E0*3*I_ESP_H2xy2z_Fy2z_a+2*1*I_ESP_Fy2z_Fy2z;
    abcd[iGrid*1260+177] = 4.0E0*I_ESP_K4x3z_Fy2z_aa-2.0E0*2*I_ESP_H2x3z_Fy2z_a-2.0E0*3*I_ESP_H2x3z_Fy2z_a+2*1*I_ESP_F3z_Fy2z;
    abcd[iGrid*1260+178] = 4.0E0*I_ESP_K3x4y_Fy2z_aa-2.0E0*1*I_ESP_Hx4y_Fy2z_a-2.0E0*2*I_ESP_Hx4y_Fy2z_a;
    abcd[iGrid*1260+179] = 4.0E0*I_ESP_K3x3yz_Fy2z_aa-2.0E0*1*I_ESP_Hx3yz_Fy2z_a-2.0E0*2*I_ESP_Hx3yz_Fy2z_a;
    abcd[iGrid*1260+180] = 4.0E0*I_ESP_K3x2y2z_Fy2z_aa-2.0E0*1*I_ESP_Hx2y2z_Fy2z_a-2.0E0*2*I_ESP_Hx2y2z_Fy2z_a;
    abcd[iGrid*1260+181] = 4.0E0*I_ESP_K3xy3z_Fy2z_aa-2.0E0*1*I_ESP_Hxy3z_Fy2z_a-2.0E0*2*I_ESP_Hxy3z_Fy2z_a;
    abcd[iGrid*1260+182] = 4.0E0*I_ESP_K3x4z_Fy2z_aa-2.0E0*1*I_ESP_Hx4z_Fy2z_a-2.0E0*2*I_ESP_Hx4z_Fy2z_a;
    abcd[iGrid*1260+183] = 4.0E0*I_ESP_K2x5y_Fy2z_aa-2.0E0*1*I_ESP_H5y_Fy2z_a;
    abcd[iGrid*1260+184] = 4.0E0*I_ESP_K2x4yz_Fy2z_aa-2.0E0*1*I_ESP_H4yz_Fy2z_a;
    abcd[iGrid*1260+185] = 4.0E0*I_ESP_K2x3y2z_Fy2z_aa-2.0E0*1*I_ESP_H3y2z_Fy2z_a;
    abcd[iGrid*1260+186] = 4.0E0*I_ESP_K2x2y3z_Fy2z_aa-2.0E0*1*I_ESP_H2y3z_Fy2z_a;
    abcd[iGrid*1260+187] = 4.0E0*I_ESP_K2xy4z_Fy2z_aa-2.0E0*1*I_ESP_Hy4z_Fy2z_a;
    abcd[iGrid*1260+188] = 4.0E0*I_ESP_K2x5z_Fy2z_aa-2.0E0*1*I_ESP_H5z_Fy2z_a;
    abcd[iGrid*1260+189] = 4.0E0*I_ESP_K7x_F3z_aa-2.0E0*5*I_ESP_H5x_F3z_a-2.0E0*6*I_ESP_H5x_F3z_a+5*4*I_ESP_F3x_F3z;
    abcd[iGrid*1260+190] = 4.0E0*I_ESP_K6xy_F3z_aa-2.0E0*4*I_ESP_H4xy_F3z_a-2.0E0*5*I_ESP_H4xy_F3z_a+4*3*I_ESP_F2xy_F3z;
    abcd[iGrid*1260+191] = 4.0E0*I_ESP_K6xz_F3z_aa-2.0E0*4*I_ESP_H4xz_F3z_a-2.0E0*5*I_ESP_H4xz_F3z_a+4*3*I_ESP_F2xz_F3z;
    abcd[iGrid*1260+192] = 4.0E0*I_ESP_K5x2y_F3z_aa-2.0E0*3*I_ESP_H3x2y_F3z_a-2.0E0*4*I_ESP_H3x2y_F3z_a+3*2*I_ESP_Fx2y_F3z;
    abcd[iGrid*1260+193] = 4.0E0*I_ESP_K5xyz_F3z_aa-2.0E0*3*I_ESP_H3xyz_F3z_a-2.0E0*4*I_ESP_H3xyz_F3z_a+3*2*I_ESP_Fxyz_F3z;
    abcd[iGrid*1260+194] = 4.0E0*I_ESP_K5x2z_F3z_aa-2.0E0*3*I_ESP_H3x2z_F3z_a-2.0E0*4*I_ESP_H3x2z_F3z_a+3*2*I_ESP_Fx2z_F3z;
    abcd[iGrid*1260+195] = 4.0E0*I_ESP_K4x3y_F3z_aa-2.0E0*2*I_ESP_H2x3y_F3z_a-2.0E0*3*I_ESP_H2x3y_F3z_a+2*1*I_ESP_F3y_F3z;
    abcd[iGrid*1260+196] = 4.0E0*I_ESP_K4x2yz_F3z_aa-2.0E0*2*I_ESP_H2x2yz_F3z_a-2.0E0*3*I_ESP_H2x2yz_F3z_a+2*1*I_ESP_F2yz_F3z;
    abcd[iGrid*1260+197] = 4.0E0*I_ESP_K4xy2z_F3z_aa-2.0E0*2*I_ESP_H2xy2z_F3z_a-2.0E0*3*I_ESP_H2xy2z_F3z_a+2*1*I_ESP_Fy2z_F3z;
    abcd[iGrid*1260+198] = 4.0E0*I_ESP_K4x3z_F3z_aa-2.0E0*2*I_ESP_H2x3z_F3z_a-2.0E0*3*I_ESP_H2x3z_F3z_a+2*1*I_ESP_F3z_F3z;
    abcd[iGrid*1260+199] = 4.0E0*I_ESP_K3x4y_F3z_aa-2.0E0*1*I_ESP_Hx4y_F3z_a-2.0E0*2*I_ESP_Hx4y_F3z_a;
    abcd[iGrid*1260+200] = 4.0E0*I_ESP_K3x3yz_F3z_aa-2.0E0*1*I_ESP_Hx3yz_F3z_a-2.0E0*2*I_ESP_Hx3yz_F3z_a;
    abcd[iGrid*1260+201] = 4.0E0*I_ESP_K3x2y2z_F3z_aa-2.0E0*1*I_ESP_Hx2y2z_F3z_a-2.0E0*2*I_ESP_Hx2y2z_F3z_a;
    abcd[iGrid*1260+202] = 4.0E0*I_ESP_K3xy3z_F3z_aa-2.0E0*1*I_ESP_Hxy3z_F3z_a-2.0E0*2*I_ESP_Hxy3z_F3z_a;
    abcd[iGrid*1260+203] = 4.0E0*I_ESP_K3x4z_F3z_aa-2.0E0*1*I_ESP_Hx4z_F3z_a-2.0E0*2*I_ESP_Hx4z_F3z_a;
    abcd[iGrid*1260+204] = 4.0E0*I_ESP_K2x5y_F3z_aa-2.0E0*1*I_ESP_H5y_F3z_a;
    abcd[iGrid*1260+205] = 4.0E0*I_ESP_K2x4yz_F3z_aa-2.0E0*1*I_ESP_H4yz_F3z_a;
    abcd[iGrid*1260+206] = 4.0E0*I_ESP_K2x3y2z_F3z_aa-2.0E0*1*I_ESP_H3y2z_F3z_a;
    abcd[iGrid*1260+207] = 4.0E0*I_ESP_K2x2y3z_F3z_aa-2.0E0*1*I_ESP_H2y3z_F3z_a;
    abcd[iGrid*1260+208] = 4.0E0*I_ESP_K2xy4z_F3z_aa-2.0E0*1*I_ESP_Hy4z_F3z_a;
    abcd[iGrid*1260+209] = 4.0E0*I_ESP_K2x5z_F3z_aa-2.0E0*1*I_ESP_H5z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_F_aa
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    abcd[iGrid*1260+210] = 4.0E0*I_ESP_K6xy_F3x_aa-2.0E0*5*I_ESP_H4xy_F3x_a;
    abcd[iGrid*1260+211] = 4.0E0*I_ESP_K5x2y_F3x_aa-2.0E0*1*I_ESP_H5x_F3x_a-2.0E0*4*I_ESP_H3x2y_F3x_a+4*1*I_ESP_F3x_F3x;
    abcd[iGrid*1260+212] = 4.0E0*I_ESP_K5xyz_F3x_aa-2.0E0*4*I_ESP_H3xyz_F3x_a;
    abcd[iGrid*1260+213] = 4.0E0*I_ESP_K4x3y_F3x_aa-2.0E0*2*I_ESP_H4xy_F3x_a-2.0E0*3*I_ESP_H2x3y_F3x_a+3*2*I_ESP_F2xy_F3x;
    abcd[iGrid*1260+214] = 4.0E0*I_ESP_K4x2yz_F3x_aa-2.0E0*1*I_ESP_H4xz_F3x_a-2.0E0*3*I_ESP_H2x2yz_F3x_a+3*1*I_ESP_F2xz_F3x;
    abcd[iGrid*1260+215] = 4.0E0*I_ESP_K4xy2z_F3x_aa-2.0E0*3*I_ESP_H2xy2z_F3x_a;
    abcd[iGrid*1260+216] = 4.0E0*I_ESP_K3x4y_F3x_aa-2.0E0*3*I_ESP_H3x2y_F3x_a-2.0E0*2*I_ESP_Hx4y_F3x_a+2*3*I_ESP_Fx2y_F3x;
    abcd[iGrid*1260+217] = 4.0E0*I_ESP_K3x3yz_F3x_aa-2.0E0*2*I_ESP_H3xyz_F3x_a-2.0E0*2*I_ESP_Hx3yz_F3x_a+2*2*I_ESP_Fxyz_F3x;
    abcd[iGrid*1260+218] = 4.0E0*I_ESP_K3x2y2z_F3x_aa-2.0E0*1*I_ESP_H3x2z_F3x_a-2.0E0*2*I_ESP_Hx2y2z_F3x_a+2*1*I_ESP_Fx2z_F3x;
    abcd[iGrid*1260+219] = 4.0E0*I_ESP_K3xy3z_F3x_aa-2.0E0*2*I_ESP_Hxy3z_F3x_a;
    abcd[iGrid*1260+220] = 4.0E0*I_ESP_K2x5y_F3x_aa-2.0E0*4*I_ESP_H2x3y_F3x_a-2.0E0*1*I_ESP_H5y_F3x_a+4*I_ESP_F3y_F3x;
    abcd[iGrid*1260+221] = 4.0E0*I_ESP_K2x4yz_F3x_aa-2.0E0*3*I_ESP_H2x2yz_F3x_a-2.0E0*1*I_ESP_H4yz_F3x_a+3*I_ESP_F2yz_F3x;
    abcd[iGrid*1260+222] = 4.0E0*I_ESP_K2x3y2z_F3x_aa-2.0E0*2*I_ESP_H2xy2z_F3x_a-2.0E0*1*I_ESP_H3y2z_F3x_a+2*I_ESP_Fy2z_F3x;
    abcd[iGrid*1260+223] = 4.0E0*I_ESP_K2x2y3z_F3x_aa-2.0E0*1*I_ESP_H2x3z_F3x_a-2.0E0*1*I_ESP_H2y3z_F3x_a+1*I_ESP_F3z_F3x;
    abcd[iGrid*1260+224] = 4.0E0*I_ESP_K2xy4z_F3x_aa-2.0E0*1*I_ESP_Hy4z_F3x_a;
    abcd[iGrid*1260+225] = 4.0E0*I_ESP_Kx6y_F3x_aa-2.0E0*5*I_ESP_Hx4y_F3x_a;
    abcd[iGrid*1260+226] = 4.0E0*I_ESP_Kx5yz_F3x_aa-2.0E0*4*I_ESP_Hx3yz_F3x_a;
    abcd[iGrid*1260+227] = 4.0E0*I_ESP_Kx4y2z_F3x_aa-2.0E0*3*I_ESP_Hx2y2z_F3x_a;
    abcd[iGrid*1260+228] = 4.0E0*I_ESP_Kx3y3z_F3x_aa-2.0E0*2*I_ESP_Hxy3z_F3x_a;
    abcd[iGrid*1260+229] = 4.0E0*I_ESP_Kx2y4z_F3x_aa-2.0E0*1*I_ESP_Hx4z_F3x_a;
    abcd[iGrid*1260+230] = 4.0E0*I_ESP_Kxy5z_F3x_aa;
    abcd[iGrid*1260+231] = 4.0E0*I_ESP_K6xy_F2xy_aa-2.0E0*5*I_ESP_H4xy_F2xy_a;
    abcd[iGrid*1260+232] = 4.0E0*I_ESP_K5x2y_F2xy_aa-2.0E0*1*I_ESP_H5x_F2xy_a-2.0E0*4*I_ESP_H3x2y_F2xy_a+4*1*I_ESP_F3x_F2xy;
    abcd[iGrid*1260+233] = 4.0E0*I_ESP_K5xyz_F2xy_aa-2.0E0*4*I_ESP_H3xyz_F2xy_a;
    abcd[iGrid*1260+234] = 4.0E0*I_ESP_K4x3y_F2xy_aa-2.0E0*2*I_ESP_H4xy_F2xy_a-2.0E0*3*I_ESP_H2x3y_F2xy_a+3*2*I_ESP_F2xy_F2xy;
    abcd[iGrid*1260+235] = 4.0E0*I_ESP_K4x2yz_F2xy_aa-2.0E0*1*I_ESP_H4xz_F2xy_a-2.0E0*3*I_ESP_H2x2yz_F2xy_a+3*1*I_ESP_F2xz_F2xy;
    abcd[iGrid*1260+236] = 4.0E0*I_ESP_K4xy2z_F2xy_aa-2.0E0*3*I_ESP_H2xy2z_F2xy_a;
    abcd[iGrid*1260+237] = 4.0E0*I_ESP_K3x4y_F2xy_aa-2.0E0*3*I_ESP_H3x2y_F2xy_a-2.0E0*2*I_ESP_Hx4y_F2xy_a+2*3*I_ESP_Fx2y_F2xy;
    abcd[iGrid*1260+238] = 4.0E0*I_ESP_K3x3yz_F2xy_aa-2.0E0*2*I_ESP_H3xyz_F2xy_a-2.0E0*2*I_ESP_Hx3yz_F2xy_a+2*2*I_ESP_Fxyz_F2xy;
    abcd[iGrid*1260+239] = 4.0E0*I_ESP_K3x2y2z_F2xy_aa-2.0E0*1*I_ESP_H3x2z_F2xy_a-2.0E0*2*I_ESP_Hx2y2z_F2xy_a+2*1*I_ESP_Fx2z_F2xy;
    abcd[iGrid*1260+240] = 4.0E0*I_ESP_K3xy3z_F2xy_aa-2.0E0*2*I_ESP_Hxy3z_F2xy_a;
    abcd[iGrid*1260+241] = 4.0E0*I_ESP_K2x5y_F2xy_aa-2.0E0*4*I_ESP_H2x3y_F2xy_a-2.0E0*1*I_ESP_H5y_F2xy_a+4*I_ESP_F3y_F2xy;
    abcd[iGrid*1260+242] = 4.0E0*I_ESP_K2x4yz_F2xy_aa-2.0E0*3*I_ESP_H2x2yz_F2xy_a-2.0E0*1*I_ESP_H4yz_F2xy_a+3*I_ESP_F2yz_F2xy;
    abcd[iGrid*1260+243] = 4.0E0*I_ESP_K2x3y2z_F2xy_aa-2.0E0*2*I_ESP_H2xy2z_F2xy_a-2.0E0*1*I_ESP_H3y2z_F2xy_a+2*I_ESP_Fy2z_F2xy;
    abcd[iGrid*1260+244] = 4.0E0*I_ESP_K2x2y3z_F2xy_aa-2.0E0*1*I_ESP_H2x3z_F2xy_a-2.0E0*1*I_ESP_H2y3z_F2xy_a+1*I_ESP_F3z_F2xy;
    abcd[iGrid*1260+245] = 4.0E0*I_ESP_K2xy4z_F2xy_aa-2.0E0*1*I_ESP_Hy4z_F2xy_a;
    abcd[iGrid*1260+246] = 4.0E0*I_ESP_Kx6y_F2xy_aa-2.0E0*5*I_ESP_Hx4y_F2xy_a;
    abcd[iGrid*1260+247] = 4.0E0*I_ESP_Kx5yz_F2xy_aa-2.0E0*4*I_ESP_Hx3yz_F2xy_a;
    abcd[iGrid*1260+248] = 4.0E0*I_ESP_Kx4y2z_F2xy_aa-2.0E0*3*I_ESP_Hx2y2z_F2xy_a;
    abcd[iGrid*1260+249] = 4.0E0*I_ESP_Kx3y3z_F2xy_aa-2.0E0*2*I_ESP_Hxy3z_F2xy_a;
    abcd[iGrid*1260+250] = 4.0E0*I_ESP_Kx2y4z_F2xy_aa-2.0E0*1*I_ESP_Hx4z_F2xy_a;
    abcd[iGrid*1260+251] = 4.0E0*I_ESP_Kxy5z_F2xy_aa;
    abcd[iGrid*1260+252] = 4.0E0*I_ESP_K6xy_F2xz_aa-2.0E0*5*I_ESP_H4xy_F2xz_a;
    abcd[iGrid*1260+253] = 4.0E0*I_ESP_K5x2y_F2xz_aa-2.0E0*1*I_ESP_H5x_F2xz_a-2.0E0*4*I_ESP_H3x2y_F2xz_a+4*1*I_ESP_F3x_F2xz;
    abcd[iGrid*1260+254] = 4.0E0*I_ESP_K5xyz_F2xz_aa-2.0E0*4*I_ESP_H3xyz_F2xz_a;
    abcd[iGrid*1260+255] = 4.0E0*I_ESP_K4x3y_F2xz_aa-2.0E0*2*I_ESP_H4xy_F2xz_a-2.0E0*3*I_ESP_H2x3y_F2xz_a+3*2*I_ESP_F2xy_F2xz;
    abcd[iGrid*1260+256] = 4.0E0*I_ESP_K4x2yz_F2xz_aa-2.0E0*1*I_ESP_H4xz_F2xz_a-2.0E0*3*I_ESP_H2x2yz_F2xz_a+3*1*I_ESP_F2xz_F2xz;
    abcd[iGrid*1260+257] = 4.0E0*I_ESP_K4xy2z_F2xz_aa-2.0E0*3*I_ESP_H2xy2z_F2xz_a;
    abcd[iGrid*1260+258] = 4.0E0*I_ESP_K3x4y_F2xz_aa-2.0E0*3*I_ESP_H3x2y_F2xz_a-2.0E0*2*I_ESP_Hx4y_F2xz_a+2*3*I_ESP_Fx2y_F2xz;
    abcd[iGrid*1260+259] = 4.0E0*I_ESP_K3x3yz_F2xz_aa-2.0E0*2*I_ESP_H3xyz_F2xz_a-2.0E0*2*I_ESP_Hx3yz_F2xz_a+2*2*I_ESP_Fxyz_F2xz;
    abcd[iGrid*1260+260] = 4.0E0*I_ESP_K3x2y2z_F2xz_aa-2.0E0*1*I_ESP_H3x2z_F2xz_a-2.0E0*2*I_ESP_Hx2y2z_F2xz_a+2*1*I_ESP_Fx2z_F2xz;
    abcd[iGrid*1260+261] = 4.0E0*I_ESP_K3xy3z_F2xz_aa-2.0E0*2*I_ESP_Hxy3z_F2xz_a;
    abcd[iGrid*1260+262] = 4.0E0*I_ESP_K2x5y_F2xz_aa-2.0E0*4*I_ESP_H2x3y_F2xz_a-2.0E0*1*I_ESP_H5y_F2xz_a+4*I_ESP_F3y_F2xz;
    abcd[iGrid*1260+263] = 4.0E0*I_ESP_K2x4yz_F2xz_aa-2.0E0*3*I_ESP_H2x2yz_F2xz_a-2.0E0*1*I_ESP_H4yz_F2xz_a+3*I_ESP_F2yz_F2xz;
    abcd[iGrid*1260+264] = 4.0E0*I_ESP_K2x3y2z_F2xz_aa-2.0E0*2*I_ESP_H2xy2z_F2xz_a-2.0E0*1*I_ESP_H3y2z_F2xz_a+2*I_ESP_Fy2z_F2xz;
    abcd[iGrid*1260+265] = 4.0E0*I_ESP_K2x2y3z_F2xz_aa-2.0E0*1*I_ESP_H2x3z_F2xz_a-2.0E0*1*I_ESP_H2y3z_F2xz_a+1*I_ESP_F3z_F2xz;
    abcd[iGrid*1260+266] = 4.0E0*I_ESP_K2xy4z_F2xz_aa-2.0E0*1*I_ESP_Hy4z_F2xz_a;
    abcd[iGrid*1260+267] = 4.0E0*I_ESP_Kx6y_F2xz_aa-2.0E0*5*I_ESP_Hx4y_F2xz_a;
    abcd[iGrid*1260+268] = 4.0E0*I_ESP_Kx5yz_F2xz_aa-2.0E0*4*I_ESP_Hx3yz_F2xz_a;
    abcd[iGrid*1260+269] = 4.0E0*I_ESP_Kx4y2z_F2xz_aa-2.0E0*3*I_ESP_Hx2y2z_F2xz_a;
    abcd[iGrid*1260+270] = 4.0E0*I_ESP_Kx3y3z_F2xz_aa-2.0E0*2*I_ESP_Hxy3z_F2xz_a;
    abcd[iGrid*1260+271] = 4.0E0*I_ESP_Kx2y4z_F2xz_aa-2.0E0*1*I_ESP_Hx4z_F2xz_a;
    abcd[iGrid*1260+272] = 4.0E0*I_ESP_Kxy5z_F2xz_aa;
    abcd[iGrid*1260+273] = 4.0E0*I_ESP_K6xy_Fx2y_aa-2.0E0*5*I_ESP_H4xy_Fx2y_a;
    abcd[iGrid*1260+274] = 4.0E0*I_ESP_K5x2y_Fx2y_aa-2.0E0*1*I_ESP_H5x_Fx2y_a-2.0E0*4*I_ESP_H3x2y_Fx2y_a+4*1*I_ESP_F3x_Fx2y;
    abcd[iGrid*1260+275] = 4.0E0*I_ESP_K5xyz_Fx2y_aa-2.0E0*4*I_ESP_H3xyz_Fx2y_a;
    abcd[iGrid*1260+276] = 4.0E0*I_ESP_K4x3y_Fx2y_aa-2.0E0*2*I_ESP_H4xy_Fx2y_a-2.0E0*3*I_ESP_H2x3y_Fx2y_a+3*2*I_ESP_F2xy_Fx2y;
    abcd[iGrid*1260+277] = 4.0E0*I_ESP_K4x2yz_Fx2y_aa-2.0E0*1*I_ESP_H4xz_Fx2y_a-2.0E0*3*I_ESP_H2x2yz_Fx2y_a+3*1*I_ESP_F2xz_Fx2y;
    abcd[iGrid*1260+278] = 4.0E0*I_ESP_K4xy2z_Fx2y_aa-2.0E0*3*I_ESP_H2xy2z_Fx2y_a;
    abcd[iGrid*1260+279] = 4.0E0*I_ESP_K3x4y_Fx2y_aa-2.0E0*3*I_ESP_H3x2y_Fx2y_a-2.0E0*2*I_ESP_Hx4y_Fx2y_a+2*3*I_ESP_Fx2y_Fx2y;
    abcd[iGrid*1260+280] = 4.0E0*I_ESP_K3x3yz_Fx2y_aa-2.0E0*2*I_ESP_H3xyz_Fx2y_a-2.0E0*2*I_ESP_Hx3yz_Fx2y_a+2*2*I_ESP_Fxyz_Fx2y;
    abcd[iGrid*1260+281] = 4.0E0*I_ESP_K3x2y2z_Fx2y_aa-2.0E0*1*I_ESP_H3x2z_Fx2y_a-2.0E0*2*I_ESP_Hx2y2z_Fx2y_a+2*1*I_ESP_Fx2z_Fx2y;
    abcd[iGrid*1260+282] = 4.0E0*I_ESP_K3xy3z_Fx2y_aa-2.0E0*2*I_ESP_Hxy3z_Fx2y_a;
    abcd[iGrid*1260+283] = 4.0E0*I_ESP_K2x5y_Fx2y_aa-2.0E0*4*I_ESP_H2x3y_Fx2y_a-2.0E0*1*I_ESP_H5y_Fx2y_a+4*I_ESP_F3y_Fx2y;
    abcd[iGrid*1260+284] = 4.0E0*I_ESP_K2x4yz_Fx2y_aa-2.0E0*3*I_ESP_H2x2yz_Fx2y_a-2.0E0*1*I_ESP_H4yz_Fx2y_a+3*I_ESP_F2yz_Fx2y;
    abcd[iGrid*1260+285] = 4.0E0*I_ESP_K2x3y2z_Fx2y_aa-2.0E0*2*I_ESP_H2xy2z_Fx2y_a-2.0E0*1*I_ESP_H3y2z_Fx2y_a+2*I_ESP_Fy2z_Fx2y;
    abcd[iGrid*1260+286] = 4.0E0*I_ESP_K2x2y3z_Fx2y_aa-2.0E0*1*I_ESP_H2x3z_Fx2y_a-2.0E0*1*I_ESP_H2y3z_Fx2y_a+1*I_ESP_F3z_Fx2y;
    abcd[iGrid*1260+287] = 4.0E0*I_ESP_K2xy4z_Fx2y_aa-2.0E0*1*I_ESP_Hy4z_Fx2y_a;
    abcd[iGrid*1260+288] = 4.0E0*I_ESP_Kx6y_Fx2y_aa-2.0E0*5*I_ESP_Hx4y_Fx2y_a;
    abcd[iGrid*1260+289] = 4.0E0*I_ESP_Kx5yz_Fx2y_aa-2.0E0*4*I_ESP_Hx3yz_Fx2y_a;
    abcd[iGrid*1260+290] = 4.0E0*I_ESP_Kx4y2z_Fx2y_aa-2.0E0*3*I_ESP_Hx2y2z_Fx2y_a;
    abcd[iGrid*1260+291] = 4.0E0*I_ESP_Kx3y3z_Fx2y_aa-2.0E0*2*I_ESP_Hxy3z_Fx2y_a;
    abcd[iGrid*1260+292] = 4.0E0*I_ESP_Kx2y4z_Fx2y_aa-2.0E0*1*I_ESP_Hx4z_Fx2y_a;
    abcd[iGrid*1260+293] = 4.0E0*I_ESP_Kxy5z_Fx2y_aa;
    abcd[iGrid*1260+294] = 4.0E0*I_ESP_K6xy_Fxyz_aa-2.0E0*5*I_ESP_H4xy_Fxyz_a;
    abcd[iGrid*1260+295] = 4.0E0*I_ESP_K5x2y_Fxyz_aa-2.0E0*1*I_ESP_H5x_Fxyz_a-2.0E0*4*I_ESP_H3x2y_Fxyz_a+4*1*I_ESP_F3x_Fxyz;
    abcd[iGrid*1260+296] = 4.0E0*I_ESP_K5xyz_Fxyz_aa-2.0E0*4*I_ESP_H3xyz_Fxyz_a;
    abcd[iGrid*1260+297] = 4.0E0*I_ESP_K4x3y_Fxyz_aa-2.0E0*2*I_ESP_H4xy_Fxyz_a-2.0E0*3*I_ESP_H2x3y_Fxyz_a+3*2*I_ESP_F2xy_Fxyz;
    abcd[iGrid*1260+298] = 4.0E0*I_ESP_K4x2yz_Fxyz_aa-2.0E0*1*I_ESP_H4xz_Fxyz_a-2.0E0*3*I_ESP_H2x2yz_Fxyz_a+3*1*I_ESP_F2xz_Fxyz;
    abcd[iGrid*1260+299] = 4.0E0*I_ESP_K4xy2z_Fxyz_aa-2.0E0*3*I_ESP_H2xy2z_Fxyz_a;
    abcd[iGrid*1260+300] = 4.0E0*I_ESP_K3x4y_Fxyz_aa-2.0E0*3*I_ESP_H3x2y_Fxyz_a-2.0E0*2*I_ESP_Hx4y_Fxyz_a+2*3*I_ESP_Fx2y_Fxyz;
    abcd[iGrid*1260+301] = 4.0E0*I_ESP_K3x3yz_Fxyz_aa-2.0E0*2*I_ESP_H3xyz_Fxyz_a-2.0E0*2*I_ESP_Hx3yz_Fxyz_a+2*2*I_ESP_Fxyz_Fxyz;
    abcd[iGrid*1260+302] = 4.0E0*I_ESP_K3x2y2z_Fxyz_aa-2.0E0*1*I_ESP_H3x2z_Fxyz_a-2.0E0*2*I_ESP_Hx2y2z_Fxyz_a+2*1*I_ESP_Fx2z_Fxyz;
    abcd[iGrid*1260+303] = 4.0E0*I_ESP_K3xy3z_Fxyz_aa-2.0E0*2*I_ESP_Hxy3z_Fxyz_a;
    abcd[iGrid*1260+304] = 4.0E0*I_ESP_K2x5y_Fxyz_aa-2.0E0*4*I_ESP_H2x3y_Fxyz_a-2.0E0*1*I_ESP_H5y_Fxyz_a+4*I_ESP_F3y_Fxyz;
    abcd[iGrid*1260+305] = 4.0E0*I_ESP_K2x4yz_Fxyz_aa-2.0E0*3*I_ESP_H2x2yz_Fxyz_a-2.0E0*1*I_ESP_H4yz_Fxyz_a+3*I_ESP_F2yz_Fxyz;
    abcd[iGrid*1260+306] = 4.0E0*I_ESP_K2x3y2z_Fxyz_aa-2.0E0*2*I_ESP_H2xy2z_Fxyz_a-2.0E0*1*I_ESP_H3y2z_Fxyz_a+2*I_ESP_Fy2z_Fxyz;
    abcd[iGrid*1260+307] = 4.0E0*I_ESP_K2x2y3z_Fxyz_aa-2.0E0*1*I_ESP_H2x3z_Fxyz_a-2.0E0*1*I_ESP_H2y3z_Fxyz_a+1*I_ESP_F3z_Fxyz;
    abcd[iGrid*1260+308] = 4.0E0*I_ESP_K2xy4z_Fxyz_aa-2.0E0*1*I_ESP_Hy4z_Fxyz_a;
    abcd[iGrid*1260+309] = 4.0E0*I_ESP_Kx6y_Fxyz_aa-2.0E0*5*I_ESP_Hx4y_Fxyz_a;
    abcd[iGrid*1260+310] = 4.0E0*I_ESP_Kx5yz_Fxyz_aa-2.0E0*4*I_ESP_Hx3yz_Fxyz_a;
    abcd[iGrid*1260+311] = 4.0E0*I_ESP_Kx4y2z_Fxyz_aa-2.0E0*3*I_ESP_Hx2y2z_Fxyz_a;
    abcd[iGrid*1260+312] = 4.0E0*I_ESP_Kx3y3z_Fxyz_aa-2.0E0*2*I_ESP_Hxy3z_Fxyz_a;
    abcd[iGrid*1260+313] = 4.0E0*I_ESP_Kx2y4z_Fxyz_aa-2.0E0*1*I_ESP_Hx4z_Fxyz_a;
    abcd[iGrid*1260+314] = 4.0E0*I_ESP_Kxy5z_Fxyz_aa;
    abcd[iGrid*1260+315] = 4.0E0*I_ESP_K6xy_Fx2z_aa-2.0E0*5*I_ESP_H4xy_Fx2z_a;
    abcd[iGrid*1260+316] = 4.0E0*I_ESP_K5x2y_Fx2z_aa-2.0E0*1*I_ESP_H5x_Fx2z_a-2.0E0*4*I_ESP_H3x2y_Fx2z_a+4*1*I_ESP_F3x_Fx2z;
    abcd[iGrid*1260+317] = 4.0E0*I_ESP_K5xyz_Fx2z_aa-2.0E0*4*I_ESP_H3xyz_Fx2z_a;
    abcd[iGrid*1260+318] = 4.0E0*I_ESP_K4x3y_Fx2z_aa-2.0E0*2*I_ESP_H4xy_Fx2z_a-2.0E0*3*I_ESP_H2x3y_Fx2z_a+3*2*I_ESP_F2xy_Fx2z;
    abcd[iGrid*1260+319] = 4.0E0*I_ESP_K4x2yz_Fx2z_aa-2.0E0*1*I_ESP_H4xz_Fx2z_a-2.0E0*3*I_ESP_H2x2yz_Fx2z_a+3*1*I_ESP_F2xz_Fx2z;
    abcd[iGrid*1260+320] = 4.0E0*I_ESP_K4xy2z_Fx2z_aa-2.0E0*3*I_ESP_H2xy2z_Fx2z_a;
    abcd[iGrid*1260+321] = 4.0E0*I_ESP_K3x4y_Fx2z_aa-2.0E0*3*I_ESP_H3x2y_Fx2z_a-2.0E0*2*I_ESP_Hx4y_Fx2z_a+2*3*I_ESP_Fx2y_Fx2z;
    abcd[iGrid*1260+322] = 4.0E0*I_ESP_K3x3yz_Fx2z_aa-2.0E0*2*I_ESP_H3xyz_Fx2z_a-2.0E0*2*I_ESP_Hx3yz_Fx2z_a+2*2*I_ESP_Fxyz_Fx2z;
    abcd[iGrid*1260+323] = 4.0E0*I_ESP_K3x2y2z_Fx2z_aa-2.0E0*1*I_ESP_H3x2z_Fx2z_a-2.0E0*2*I_ESP_Hx2y2z_Fx2z_a+2*1*I_ESP_Fx2z_Fx2z;
    abcd[iGrid*1260+324] = 4.0E0*I_ESP_K3xy3z_Fx2z_aa-2.0E0*2*I_ESP_Hxy3z_Fx2z_a;
    abcd[iGrid*1260+325] = 4.0E0*I_ESP_K2x5y_Fx2z_aa-2.0E0*4*I_ESP_H2x3y_Fx2z_a-2.0E0*1*I_ESP_H5y_Fx2z_a+4*I_ESP_F3y_Fx2z;
    abcd[iGrid*1260+326] = 4.0E0*I_ESP_K2x4yz_Fx2z_aa-2.0E0*3*I_ESP_H2x2yz_Fx2z_a-2.0E0*1*I_ESP_H4yz_Fx2z_a+3*I_ESP_F2yz_Fx2z;
    abcd[iGrid*1260+327] = 4.0E0*I_ESP_K2x3y2z_Fx2z_aa-2.0E0*2*I_ESP_H2xy2z_Fx2z_a-2.0E0*1*I_ESP_H3y2z_Fx2z_a+2*I_ESP_Fy2z_Fx2z;
    abcd[iGrid*1260+328] = 4.0E0*I_ESP_K2x2y3z_Fx2z_aa-2.0E0*1*I_ESP_H2x3z_Fx2z_a-2.0E0*1*I_ESP_H2y3z_Fx2z_a+1*I_ESP_F3z_Fx2z;
    abcd[iGrid*1260+329] = 4.0E0*I_ESP_K2xy4z_Fx2z_aa-2.0E0*1*I_ESP_Hy4z_Fx2z_a;
    abcd[iGrid*1260+330] = 4.0E0*I_ESP_Kx6y_Fx2z_aa-2.0E0*5*I_ESP_Hx4y_Fx2z_a;
    abcd[iGrid*1260+331] = 4.0E0*I_ESP_Kx5yz_Fx2z_aa-2.0E0*4*I_ESP_Hx3yz_Fx2z_a;
    abcd[iGrid*1260+332] = 4.0E0*I_ESP_Kx4y2z_Fx2z_aa-2.0E0*3*I_ESP_Hx2y2z_Fx2z_a;
    abcd[iGrid*1260+333] = 4.0E0*I_ESP_Kx3y3z_Fx2z_aa-2.0E0*2*I_ESP_Hxy3z_Fx2z_a;
    abcd[iGrid*1260+334] = 4.0E0*I_ESP_Kx2y4z_Fx2z_aa-2.0E0*1*I_ESP_Hx4z_Fx2z_a;
    abcd[iGrid*1260+335] = 4.0E0*I_ESP_Kxy5z_Fx2z_aa;
    abcd[iGrid*1260+336] = 4.0E0*I_ESP_K6xy_F3y_aa-2.0E0*5*I_ESP_H4xy_F3y_a;
    abcd[iGrid*1260+337] = 4.0E0*I_ESP_K5x2y_F3y_aa-2.0E0*1*I_ESP_H5x_F3y_a-2.0E0*4*I_ESP_H3x2y_F3y_a+4*1*I_ESP_F3x_F3y;
    abcd[iGrid*1260+338] = 4.0E0*I_ESP_K5xyz_F3y_aa-2.0E0*4*I_ESP_H3xyz_F3y_a;
    abcd[iGrid*1260+339] = 4.0E0*I_ESP_K4x3y_F3y_aa-2.0E0*2*I_ESP_H4xy_F3y_a-2.0E0*3*I_ESP_H2x3y_F3y_a+3*2*I_ESP_F2xy_F3y;
    abcd[iGrid*1260+340] = 4.0E0*I_ESP_K4x2yz_F3y_aa-2.0E0*1*I_ESP_H4xz_F3y_a-2.0E0*3*I_ESP_H2x2yz_F3y_a+3*1*I_ESP_F2xz_F3y;
    abcd[iGrid*1260+341] = 4.0E0*I_ESP_K4xy2z_F3y_aa-2.0E0*3*I_ESP_H2xy2z_F3y_a;
    abcd[iGrid*1260+342] = 4.0E0*I_ESP_K3x4y_F3y_aa-2.0E0*3*I_ESP_H3x2y_F3y_a-2.0E0*2*I_ESP_Hx4y_F3y_a+2*3*I_ESP_Fx2y_F3y;
    abcd[iGrid*1260+343] = 4.0E0*I_ESP_K3x3yz_F3y_aa-2.0E0*2*I_ESP_H3xyz_F3y_a-2.0E0*2*I_ESP_Hx3yz_F3y_a+2*2*I_ESP_Fxyz_F3y;
    abcd[iGrid*1260+344] = 4.0E0*I_ESP_K3x2y2z_F3y_aa-2.0E0*1*I_ESP_H3x2z_F3y_a-2.0E0*2*I_ESP_Hx2y2z_F3y_a+2*1*I_ESP_Fx2z_F3y;
    abcd[iGrid*1260+345] = 4.0E0*I_ESP_K3xy3z_F3y_aa-2.0E0*2*I_ESP_Hxy3z_F3y_a;
    abcd[iGrid*1260+346] = 4.0E0*I_ESP_K2x5y_F3y_aa-2.0E0*4*I_ESP_H2x3y_F3y_a-2.0E0*1*I_ESP_H5y_F3y_a+4*I_ESP_F3y_F3y;
    abcd[iGrid*1260+347] = 4.0E0*I_ESP_K2x4yz_F3y_aa-2.0E0*3*I_ESP_H2x2yz_F3y_a-2.0E0*1*I_ESP_H4yz_F3y_a+3*I_ESP_F2yz_F3y;
    abcd[iGrid*1260+348] = 4.0E0*I_ESP_K2x3y2z_F3y_aa-2.0E0*2*I_ESP_H2xy2z_F3y_a-2.0E0*1*I_ESP_H3y2z_F3y_a+2*I_ESP_Fy2z_F3y;
    abcd[iGrid*1260+349] = 4.0E0*I_ESP_K2x2y3z_F3y_aa-2.0E0*1*I_ESP_H2x3z_F3y_a-2.0E0*1*I_ESP_H2y3z_F3y_a+1*I_ESP_F3z_F3y;
    abcd[iGrid*1260+350] = 4.0E0*I_ESP_K2xy4z_F3y_aa-2.0E0*1*I_ESP_Hy4z_F3y_a;
    abcd[iGrid*1260+351] = 4.0E0*I_ESP_Kx6y_F3y_aa-2.0E0*5*I_ESP_Hx4y_F3y_a;
    abcd[iGrid*1260+352] = 4.0E0*I_ESP_Kx5yz_F3y_aa-2.0E0*4*I_ESP_Hx3yz_F3y_a;
    abcd[iGrid*1260+353] = 4.0E0*I_ESP_Kx4y2z_F3y_aa-2.0E0*3*I_ESP_Hx2y2z_F3y_a;
    abcd[iGrid*1260+354] = 4.0E0*I_ESP_Kx3y3z_F3y_aa-2.0E0*2*I_ESP_Hxy3z_F3y_a;
    abcd[iGrid*1260+355] = 4.0E0*I_ESP_Kx2y4z_F3y_aa-2.0E0*1*I_ESP_Hx4z_F3y_a;
    abcd[iGrid*1260+356] = 4.0E0*I_ESP_Kxy5z_F3y_aa;
    abcd[iGrid*1260+357] = 4.0E0*I_ESP_K6xy_F2yz_aa-2.0E0*5*I_ESP_H4xy_F2yz_a;
    abcd[iGrid*1260+358] = 4.0E0*I_ESP_K5x2y_F2yz_aa-2.0E0*1*I_ESP_H5x_F2yz_a-2.0E0*4*I_ESP_H3x2y_F2yz_a+4*1*I_ESP_F3x_F2yz;
    abcd[iGrid*1260+359] = 4.0E0*I_ESP_K5xyz_F2yz_aa-2.0E0*4*I_ESP_H3xyz_F2yz_a;
    abcd[iGrid*1260+360] = 4.0E0*I_ESP_K4x3y_F2yz_aa-2.0E0*2*I_ESP_H4xy_F2yz_a-2.0E0*3*I_ESP_H2x3y_F2yz_a+3*2*I_ESP_F2xy_F2yz;
    abcd[iGrid*1260+361] = 4.0E0*I_ESP_K4x2yz_F2yz_aa-2.0E0*1*I_ESP_H4xz_F2yz_a-2.0E0*3*I_ESP_H2x2yz_F2yz_a+3*1*I_ESP_F2xz_F2yz;
    abcd[iGrid*1260+362] = 4.0E0*I_ESP_K4xy2z_F2yz_aa-2.0E0*3*I_ESP_H2xy2z_F2yz_a;
    abcd[iGrid*1260+363] = 4.0E0*I_ESP_K3x4y_F2yz_aa-2.0E0*3*I_ESP_H3x2y_F2yz_a-2.0E0*2*I_ESP_Hx4y_F2yz_a+2*3*I_ESP_Fx2y_F2yz;
    abcd[iGrid*1260+364] = 4.0E0*I_ESP_K3x3yz_F2yz_aa-2.0E0*2*I_ESP_H3xyz_F2yz_a-2.0E0*2*I_ESP_Hx3yz_F2yz_a+2*2*I_ESP_Fxyz_F2yz;
    abcd[iGrid*1260+365] = 4.0E0*I_ESP_K3x2y2z_F2yz_aa-2.0E0*1*I_ESP_H3x2z_F2yz_a-2.0E0*2*I_ESP_Hx2y2z_F2yz_a+2*1*I_ESP_Fx2z_F2yz;
    abcd[iGrid*1260+366] = 4.0E0*I_ESP_K3xy3z_F2yz_aa-2.0E0*2*I_ESP_Hxy3z_F2yz_a;
    abcd[iGrid*1260+367] = 4.0E0*I_ESP_K2x5y_F2yz_aa-2.0E0*4*I_ESP_H2x3y_F2yz_a-2.0E0*1*I_ESP_H5y_F2yz_a+4*I_ESP_F3y_F2yz;
    abcd[iGrid*1260+368] = 4.0E0*I_ESP_K2x4yz_F2yz_aa-2.0E0*3*I_ESP_H2x2yz_F2yz_a-2.0E0*1*I_ESP_H4yz_F2yz_a+3*I_ESP_F2yz_F2yz;
    abcd[iGrid*1260+369] = 4.0E0*I_ESP_K2x3y2z_F2yz_aa-2.0E0*2*I_ESP_H2xy2z_F2yz_a-2.0E0*1*I_ESP_H3y2z_F2yz_a+2*I_ESP_Fy2z_F2yz;
    abcd[iGrid*1260+370] = 4.0E0*I_ESP_K2x2y3z_F2yz_aa-2.0E0*1*I_ESP_H2x3z_F2yz_a-2.0E0*1*I_ESP_H2y3z_F2yz_a+1*I_ESP_F3z_F2yz;
    abcd[iGrid*1260+371] = 4.0E0*I_ESP_K2xy4z_F2yz_aa-2.0E0*1*I_ESP_Hy4z_F2yz_a;
    abcd[iGrid*1260+372] = 4.0E0*I_ESP_Kx6y_F2yz_aa-2.0E0*5*I_ESP_Hx4y_F2yz_a;
    abcd[iGrid*1260+373] = 4.0E0*I_ESP_Kx5yz_F2yz_aa-2.0E0*4*I_ESP_Hx3yz_F2yz_a;
    abcd[iGrid*1260+374] = 4.0E0*I_ESP_Kx4y2z_F2yz_aa-2.0E0*3*I_ESP_Hx2y2z_F2yz_a;
    abcd[iGrid*1260+375] = 4.0E0*I_ESP_Kx3y3z_F2yz_aa-2.0E0*2*I_ESP_Hxy3z_F2yz_a;
    abcd[iGrid*1260+376] = 4.0E0*I_ESP_Kx2y4z_F2yz_aa-2.0E0*1*I_ESP_Hx4z_F2yz_a;
    abcd[iGrid*1260+377] = 4.0E0*I_ESP_Kxy5z_F2yz_aa;
    abcd[iGrid*1260+378] = 4.0E0*I_ESP_K6xy_Fy2z_aa-2.0E0*5*I_ESP_H4xy_Fy2z_a;
    abcd[iGrid*1260+379] = 4.0E0*I_ESP_K5x2y_Fy2z_aa-2.0E0*1*I_ESP_H5x_Fy2z_a-2.0E0*4*I_ESP_H3x2y_Fy2z_a+4*1*I_ESP_F3x_Fy2z;
    abcd[iGrid*1260+380] = 4.0E0*I_ESP_K5xyz_Fy2z_aa-2.0E0*4*I_ESP_H3xyz_Fy2z_a;
    abcd[iGrid*1260+381] = 4.0E0*I_ESP_K4x3y_Fy2z_aa-2.0E0*2*I_ESP_H4xy_Fy2z_a-2.0E0*3*I_ESP_H2x3y_Fy2z_a+3*2*I_ESP_F2xy_Fy2z;
    abcd[iGrid*1260+382] = 4.0E0*I_ESP_K4x2yz_Fy2z_aa-2.0E0*1*I_ESP_H4xz_Fy2z_a-2.0E0*3*I_ESP_H2x2yz_Fy2z_a+3*1*I_ESP_F2xz_Fy2z;
    abcd[iGrid*1260+383] = 4.0E0*I_ESP_K4xy2z_Fy2z_aa-2.0E0*3*I_ESP_H2xy2z_Fy2z_a;
    abcd[iGrid*1260+384] = 4.0E0*I_ESP_K3x4y_Fy2z_aa-2.0E0*3*I_ESP_H3x2y_Fy2z_a-2.0E0*2*I_ESP_Hx4y_Fy2z_a+2*3*I_ESP_Fx2y_Fy2z;
    abcd[iGrid*1260+385] = 4.0E0*I_ESP_K3x3yz_Fy2z_aa-2.0E0*2*I_ESP_H3xyz_Fy2z_a-2.0E0*2*I_ESP_Hx3yz_Fy2z_a+2*2*I_ESP_Fxyz_Fy2z;
    abcd[iGrid*1260+386] = 4.0E0*I_ESP_K3x2y2z_Fy2z_aa-2.0E0*1*I_ESP_H3x2z_Fy2z_a-2.0E0*2*I_ESP_Hx2y2z_Fy2z_a+2*1*I_ESP_Fx2z_Fy2z;
    abcd[iGrid*1260+387] = 4.0E0*I_ESP_K3xy3z_Fy2z_aa-2.0E0*2*I_ESP_Hxy3z_Fy2z_a;
    abcd[iGrid*1260+388] = 4.0E0*I_ESP_K2x5y_Fy2z_aa-2.0E0*4*I_ESP_H2x3y_Fy2z_a-2.0E0*1*I_ESP_H5y_Fy2z_a+4*I_ESP_F3y_Fy2z;
    abcd[iGrid*1260+389] = 4.0E0*I_ESP_K2x4yz_Fy2z_aa-2.0E0*3*I_ESP_H2x2yz_Fy2z_a-2.0E0*1*I_ESP_H4yz_Fy2z_a+3*I_ESP_F2yz_Fy2z;
    abcd[iGrid*1260+390] = 4.0E0*I_ESP_K2x3y2z_Fy2z_aa-2.0E0*2*I_ESP_H2xy2z_Fy2z_a-2.0E0*1*I_ESP_H3y2z_Fy2z_a+2*I_ESP_Fy2z_Fy2z;
    abcd[iGrid*1260+391] = 4.0E0*I_ESP_K2x2y3z_Fy2z_aa-2.0E0*1*I_ESP_H2x3z_Fy2z_a-2.0E0*1*I_ESP_H2y3z_Fy2z_a+1*I_ESP_F3z_Fy2z;
    abcd[iGrid*1260+392] = 4.0E0*I_ESP_K2xy4z_Fy2z_aa-2.0E0*1*I_ESP_Hy4z_Fy2z_a;
    abcd[iGrid*1260+393] = 4.0E0*I_ESP_Kx6y_Fy2z_aa-2.0E0*5*I_ESP_Hx4y_Fy2z_a;
    abcd[iGrid*1260+394] = 4.0E0*I_ESP_Kx5yz_Fy2z_aa-2.0E0*4*I_ESP_Hx3yz_Fy2z_a;
    abcd[iGrid*1260+395] = 4.0E0*I_ESP_Kx4y2z_Fy2z_aa-2.0E0*3*I_ESP_Hx2y2z_Fy2z_a;
    abcd[iGrid*1260+396] = 4.0E0*I_ESP_Kx3y3z_Fy2z_aa-2.0E0*2*I_ESP_Hxy3z_Fy2z_a;
    abcd[iGrid*1260+397] = 4.0E0*I_ESP_Kx2y4z_Fy2z_aa-2.0E0*1*I_ESP_Hx4z_Fy2z_a;
    abcd[iGrid*1260+398] = 4.0E0*I_ESP_Kxy5z_Fy2z_aa;
    abcd[iGrid*1260+399] = 4.0E0*I_ESP_K6xy_F3z_aa-2.0E0*5*I_ESP_H4xy_F3z_a;
    abcd[iGrid*1260+400] = 4.0E0*I_ESP_K5x2y_F3z_aa-2.0E0*1*I_ESP_H5x_F3z_a-2.0E0*4*I_ESP_H3x2y_F3z_a+4*1*I_ESP_F3x_F3z;
    abcd[iGrid*1260+401] = 4.0E0*I_ESP_K5xyz_F3z_aa-2.0E0*4*I_ESP_H3xyz_F3z_a;
    abcd[iGrid*1260+402] = 4.0E0*I_ESP_K4x3y_F3z_aa-2.0E0*2*I_ESP_H4xy_F3z_a-2.0E0*3*I_ESP_H2x3y_F3z_a+3*2*I_ESP_F2xy_F3z;
    abcd[iGrid*1260+403] = 4.0E0*I_ESP_K4x2yz_F3z_aa-2.0E0*1*I_ESP_H4xz_F3z_a-2.0E0*3*I_ESP_H2x2yz_F3z_a+3*1*I_ESP_F2xz_F3z;
    abcd[iGrid*1260+404] = 4.0E0*I_ESP_K4xy2z_F3z_aa-2.0E0*3*I_ESP_H2xy2z_F3z_a;
    abcd[iGrid*1260+405] = 4.0E0*I_ESP_K3x4y_F3z_aa-2.0E0*3*I_ESP_H3x2y_F3z_a-2.0E0*2*I_ESP_Hx4y_F3z_a+2*3*I_ESP_Fx2y_F3z;
    abcd[iGrid*1260+406] = 4.0E0*I_ESP_K3x3yz_F3z_aa-2.0E0*2*I_ESP_H3xyz_F3z_a-2.0E0*2*I_ESP_Hx3yz_F3z_a+2*2*I_ESP_Fxyz_F3z;
    abcd[iGrid*1260+407] = 4.0E0*I_ESP_K3x2y2z_F3z_aa-2.0E0*1*I_ESP_H3x2z_F3z_a-2.0E0*2*I_ESP_Hx2y2z_F3z_a+2*1*I_ESP_Fx2z_F3z;
    abcd[iGrid*1260+408] = 4.0E0*I_ESP_K3xy3z_F3z_aa-2.0E0*2*I_ESP_Hxy3z_F3z_a;
    abcd[iGrid*1260+409] = 4.0E0*I_ESP_K2x5y_F3z_aa-2.0E0*4*I_ESP_H2x3y_F3z_a-2.0E0*1*I_ESP_H5y_F3z_a+4*I_ESP_F3y_F3z;
    abcd[iGrid*1260+410] = 4.0E0*I_ESP_K2x4yz_F3z_aa-2.0E0*3*I_ESP_H2x2yz_F3z_a-2.0E0*1*I_ESP_H4yz_F3z_a+3*I_ESP_F2yz_F3z;
    abcd[iGrid*1260+411] = 4.0E0*I_ESP_K2x3y2z_F3z_aa-2.0E0*2*I_ESP_H2xy2z_F3z_a-2.0E0*1*I_ESP_H3y2z_F3z_a+2*I_ESP_Fy2z_F3z;
    abcd[iGrid*1260+412] = 4.0E0*I_ESP_K2x2y3z_F3z_aa-2.0E0*1*I_ESP_H2x3z_F3z_a-2.0E0*1*I_ESP_H2y3z_F3z_a+1*I_ESP_F3z_F3z;
    abcd[iGrid*1260+413] = 4.0E0*I_ESP_K2xy4z_F3z_aa-2.0E0*1*I_ESP_Hy4z_F3z_a;
    abcd[iGrid*1260+414] = 4.0E0*I_ESP_Kx6y_F3z_aa-2.0E0*5*I_ESP_Hx4y_F3z_a;
    abcd[iGrid*1260+415] = 4.0E0*I_ESP_Kx5yz_F3z_aa-2.0E0*4*I_ESP_Hx3yz_F3z_a;
    abcd[iGrid*1260+416] = 4.0E0*I_ESP_Kx4y2z_F3z_aa-2.0E0*3*I_ESP_Hx2y2z_F3z_a;
    abcd[iGrid*1260+417] = 4.0E0*I_ESP_Kx3y3z_F3z_aa-2.0E0*2*I_ESP_Hxy3z_F3z_a;
    abcd[iGrid*1260+418] = 4.0E0*I_ESP_Kx2y4z_F3z_aa-2.0E0*1*I_ESP_Hx4z_F3z_a;
    abcd[iGrid*1260+419] = 4.0E0*I_ESP_Kxy5z_F3z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_F_aa
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    abcd[iGrid*1260+420] = 4.0E0*I_ESP_K6xz_F3x_aa-2.0E0*5*I_ESP_H4xz_F3x_a;
    abcd[iGrid*1260+421] = 4.0E0*I_ESP_K5xyz_F3x_aa-2.0E0*4*I_ESP_H3xyz_F3x_a;
    abcd[iGrid*1260+422] = 4.0E0*I_ESP_K5x2z_F3x_aa-2.0E0*1*I_ESP_H5x_F3x_a-2.0E0*4*I_ESP_H3x2z_F3x_a+4*1*I_ESP_F3x_F3x;
    abcd[iGrid*1260+423] = 4.0E0*I_ESP_K4x2yz_F3x_aa-2.0E0*3*I_ESP_H2x2yz_F3x_a;
    abcd[iGrid*1260+424] = 4.0E0*I_ESP_K4xy2z_F3x_aa-2.0E0*1*I_ESP_H4xy_F3x_a-2.0E0*3*I_ESP_H2xy2z_F3x_a+3*1*I_ESP_F2xy_F3x;
    abcd[iGrid*1260+425] = 4.0E0*I_ESP_K4x3z_F3x_aa-2.0E0*2*I_ESP_H4xz_F3x_a-2.0E0*3*I_ESP_H2x3z_F3x_a+3*2*I_ESP_F2xz_F3x;
    abcd[iGrid*1260+426] = 4.0E0*I_ESP_K3x3yz_F3x_aa-2.0E0*2*I_ESP_Hx3yz_F3x_a;
    abcd[iGrid*1260+427] = 4.0E0*I_ESP_K3x2y2z_F3x_aa-2.0E0*1*I_ESP_H3x2y_F3x_a-2.0E0*2*I_ESP_Hx2y2z_F3x_a+2*1*I_ESP_Fx2y_F3x;
    abcd[iGrid*1260+428] = 4.0E0*I_ESP_K3xy3z_F3x_aa-2.0E0*2*I_ESP_H3xyz_F3x_a-2.0E0*2*I_ESP_Hxy3z_F3x_a+2*2*I_ESP_Fxyz_F3x;
    abcd[iGrid*1260+429] = 4.0E0*I_ESP_K3x4z_F3x_aa-2.0E0*3*I_ESP_H3x2z_F3x_a-2.0E0*2*I_ESP_Hx4z_F3x_a+2*3*I_ESP_Fx2z_F3x;
    abcd[iGrid*1260+430] = 4.0E0*I_ESP_K2x4yz_F3x_aa-2.0E0*1*I_ESP_H4yz_F3x_a;
    abcd[iGrid*1260+431] = 4.0E0*I_ESP_K2x3y2z_F3x_aa-2.0E0*1*I_ESP_H2x3y_F3x_a-2.0E0*1*I_ESP_H3y2z_F3x_a+1*I_ESP_F3y_F3x;
    abcd[iGrid*1260+432] = 4.0E0*I_ESP_K2x2y3z_F3x_aa-2.0E0*2*I_ESP_H2x2yz_F3x_a-2.0E0*1*I_ESP_H2y3z_F3x_a+2*I_ESP_F2yz_F3x;
    abcd[iGrid*1260+433] = 4.0E0*I_ESP_K2xy4z_F3x_aa-2.0E0*3*I_ESP_H2xy2z_F3x_a-2.0E0*1*I_ESP_Hy4z_F3x_a+3*I_ESP_Fy2z_F3x;
    abcd[iGrid*1260+434] = 4.0E0*I_ESP_K2x5z_F3x_aa-2.0E0*4*I_ESP_H2x3z_F3x_a-2.0E0*1*I_ESP_H5z_F3x_a+4*I_ESP_F3z_F3x;
    abcd[iGrid*1260+435] = 4.0E0*I_ESP_Kx5yz_F3x_aa;
    abcd[iGrid*1260+436] = 4.0E0*I_ESP_Kx4y2z_F3x_aa-2.0E0*1*I_ESP_Hx4y_F3x_a;
    abcd[iGrid*1260+437] = 4.0E0*I_ESP_Kx3y3z_F3x_aa-2.0E0*2*I_ESP_Hx3yz_F3x_a;
    abcd[iGrid*1260+438] = 4.0E0*I_ESP_Kx2y4z_F3x_aa-2.0E0*3*I_ESP_Hx2y2z_F3x_a;
    abcd[iGrid*1260+439] = 4.0E0*I_ESP_Kxy5z_F3x_aa-2.0E0*4*I_ESP_Hxy3z_F3x_a;
    abcd[iGrid*1260+440] = 4.0E0*I_ESP_Kx6z_F3x_aa-2.0E0*5*I_ESP_Hx4z_F3x_a;
    abcd[iGrid*1260+441] = 4.0E0*I_ESP_K6xz_F2xy_aa-2.0E0*5*I_ESP_H4xz_F2xy_a;
    abcd[iGrid*1260+442] = 4.0E0*I_ESP_K5xyz_F2xy_aa-2.0E0*4*I_ESP_H3xyz_F2xy_a;
    abcd[iGrid*1260+443] = 4.0E0*I_ESP_K5x2z_F2xy_aa-2.0E0*1*I_ESP_H5x_F2xy_a-2.0E0*4*I_ESP_H3x2z_F2xy_a+4*1*I_ESP_F3x_F2xy;
    abcd[iGrid*1260+444] = 4.0E0*I_ESP_K4x2yz_F2xy_aa-2.0E0*3*I_ESP_H2x2yz_F2xy_a;
    abcd[iGrid*1260+445] = 4.0E0*I_ESP_K4xy2z_F2xy_aa-2.0E0*1*I_ESP_H4xy_F2xy_a-2.0E0*3*I_ESP_H2xy2z_F2xy_a+3*1*I_ESP_F2xy_F2xy;
    abcd[iGrid*1260+446] = 4.0E0*I_ESP_K4x3z_F2xy_aa-2.0E0*2*I_ESP_H4xz_F2xy_a-2.0E0*3*I_ESP_H2x3z_F2xy_a+3*2*I_ESP_F2xz_F2xy;
    abcd[iGrid*1260+447] = 4.0E0*I_ESP_K3x3yz_F2xy_aa-2.0E0*2*I_ESP_Hx3yz_F2xy_a;
    abcd[iGrid*1260+448] = 4.0E0*I_ESP_K3x2y2z_F2xy_aa-2.0E0*1*I_ESP_H3x2y_F2xy_a-2.0E0*2*I_ESP_Hx2y2z_F2xy_a+2*1*I_ESP_Fx2y_F2xy;
    abcd[iGrid*1260+449] = 4.0E0*I_ESP_K3xy3z_F2xy_aa-2.0E0*2*I_ESP_H3xyz_F2xy_a-2.0E0*2*I_ESP_Hxy3z_F2xy_a+2*2*I_ESP_Fxyz_F2xy;
    abcd[iGrid*1260+450] = 4.0E0*I_ESP_K3x4z_F2xy_aa-2.0E0*3*I_ESP_H3x2z_F2xy_a-2.0E0*2*I_ESP_Hx4z_F2xy_a+2*3*I_ESP_Fx2z_F2xy;
    abcd[iGrid*1260+451] = 4.0E0*I_ESP_K2x4yz_F2xy_aa-2.0E0*1*I_ESP_H4yz_F2xy_a;
    abcd[iGrid*1260+452] = 4.0E0*I_ESP_K2x3y2z_F2xy_aa-2.0E0*1*I_ESP_H2x3y_F2xy_a-2.0E0*1*I_ESP_H3y2z_F2xy_a+1*I_ESP_F3y_F2xy;
    abcd[iGrid*1260+453] = 4.0E0*I_ESP_K2x2y3z_F2xy_aa-2.0E0*2*I_ESP_H2x2yz_F2xy_a-2.0E0*1*I_ESP_H2y3z_F2xy_a+2*I_ESP_F2yz_F2xy;
    abcd[iGrid*1260+454] = 4.0E0*I_ESP_K2xy4z_F2xy_aa-2.0E0*3*I_ESP_H2xy2z_F2xy_a-2.0E0*1*I_ESP_Hy4z_F2xy_a+3*I_ESP_Fy2z_F2xy;
    abcd[iGrid*1260+455] = 4.0E0*I_ESP_K2x5z_F2xy_aa-2.0E0*4*I_ESP_H2x3z_F2xy_a-2.0E0*1*I_ESP_H5z_F2xy_a+4*I_ESP_F3z_F2xy;
    abcd[iGrid*1260+456] = 4.0E0*I_ESP_Kx5yz_F2xy_aa;
    abcd[iGrid*1260+457] = 4.0E0*I_ESP_Kx4y2z_F2xy_aa-2.0E0*1*I_ESP_Hx4y_F2xy_a;
    abcd[iGrid*1260+458] = 4.0E0*I_ESP_Kx3y3z_F2xy_aa-2.0E0*2*I_ESP_Hx3yz_F2xy_a;
    abcd[iGrid*1260+459] = 4.0E0*I_ESP_Kx2y4z_F2xy_aa-2.0E0*3*I_ESP_Hx2y2z_F2xy_a;
    abcd[iGrid*1260+460] = 4.0E0*I_ESP_Kxy5z_F2xy_aa-2.0E0*4*I_ESP_Hxy3z_F2xy_a;
    abcd[iGrid*1260+461] = 4.0E0*I_ESP_Kx6z_F2xy_aa-2.0E0*5*I_ESP_Hx4z_F2xy_a;
    abcd[iGrid*1260+462] = 4.0E0*I_ESP_K6xz_F2xz_aa-2.0E0*5*I_ESP_H4xz_F2xz_a;
    abcd[iGrid*1260+463] = 4.0E0*I_ESP_K5xyz_F2xz_aa-2.0E0*4*I_ESP_H3xyz_F2xz_a;
    abcd[iGrid*1260+464] = 4.0E0*I_ESP_K5x2z_F2xz_aa-2.0E0*1*I_ESP_H5x_F2xz_a-2.0E0*4*I_ESP_H3x2z_F2xz_a+4*1*I_ESP_F3x_F2xz;
    abcd[iGrid*1260+465] = 4.0E0*I_ESP_K4x2yz_F2xz_aa-2.0E0*3*I_ESP_H2x2yz_F2xz_a;
    abcd[iGrid*1260+466] = 4.0E0*I_ESP_K4xy2z_F2xz_aa-2.0E0*1*I_ESP_H4xy_F2xz_a-2.0E0*3*I_ESP_H2xy2z_F2xz_a+3*1*I_ESP_F2xy_F2xz;
    abcd[iGrid*1260+467] = 4.0E0*I_ESP_K4x3z_F2xz_aa-2.0E0*2*I_ESP_H4xz_F2xz_a-2.0E0*3*I_ESP_H2x3z_F2xz_a+3*2*I_ESP_F2xz_F2xz;
    abcd[iGrid*1260+468] = 4.0E0*I_ESP_K3x3yz_F2xz_aa-2.0E0*2*I_ESP_Hx3yz_F2xz_a;
    abcd[iGrid*1260+469] = 4.0E0*I_ESP_K3x2y2z_F2xz_aa-2.0E0*1*I_ESP_H3x2y_F2xz_a-2.0E0*2*I_ESP_Hx2y2z_F2xz_a+2*1*I_ESP_Fx2y_F2xz;
    abcd[iGrid*1260+470] = 4.0E0*I_ESP_K3xy3z_F2xz_aa-2.0E0*2*I_ESP_H3xyz_F2xz_a-2.0E0*2*I_ESP_Hxy3z_F2xz_a+2*2*I_ESP_Fxyz_F2xz;
    abcd[iGrid*1260+471] = 4.0E0*I_ESP_K3x4z_F2xz_aa-2.0E0*3*I_ESP_H3x2z_F2xz_a-2.0E0*2*I_ESP_Hx4z_F2xz_a+2*3*I_ESP_Fx2z_F2xz;
    abcd[iGrid*1260+472] = 4.0E0*I_ESP_K2x4yz_F2xz_aa-2.0E0*1*I_ESP_H4yz_F2xz_a;
    abcd[iGrid*1260+473] = 4.0E0*I_ESP_K2x3y2z_F2xz_aa-2.0E0*1*I_ESP_H2x3y_F2xz_a-2.0E0*1*I_ESP_H3y2z_F2xz_a+1*I_ESP_F3y_F2xz;
    abcd[iGrid*1260+474] = 4.0E0*I_ESP_K2x2y3z_F2xz_aa-2.0E0*2*I_ESP_H2x2yz_F2xz_a-2.0E0*1*I_ESP_H2y3z_F2xz_a+2*I_ESP_F2yz_F2xz;
    abcd[iGrid*1260+475] = 4.0E0*I_ESP_K2xy4z_F2xz_aa-2.0E0*3*I_ESP_H2xy2z_F2xz_a-2.0E0*1*I_ESP_Hy4z_F2xz_a+3*I_ESP_Fy2z_F2xz;
    abcd[iGrid*1260+476] = 4.0E0*I_ESP_K2x5z_F2xz_aa-2.0E0*4*I_ESP_H2x3z_F2xz_a-2.0E0*1*I_ESP_H5z_F2xz_a+4*I_ESP_F3z_F2xz;
    abcd[iGrid*1260+477] = 4.0E0*I_ESP_Kx5yz_F2xz_aa;
    abcd[iGrid*1260+478] = 4.0E0*I_ESP_Kx4y2z_F2xz_aa-2.0E0*1*I_ESP_Hx4y_F2xz_a;
    abcd[iGrid*1260+479] = 4.0E0*I_ESP_Kx3y3z_F2xz_aa-2.0E0*2*I_ESP_Hx3yz_F2xz_a;
    abcd[iGrid*1260+480] = 4.0E0*I_ESP_Kx2y4z_F2xz_aa-2.0E0*3*I_ESP_Hx2y2z_F2xz_a;
    abcd[iGrid*1260+481] = 4.0E0*I_ESP_Kxy5z_F2xz_aa-2.0E0*4*I_ESP_Hxy3z_F2xz_a;
    abcd[iGrid*1260+482] = 4.0E0*I_ESP_Kx6z_F2xz_aa-2.0E0*5*I_ESP_Hx4z_F2xz_a;
    abcd[iGrid*1260+483] = 4.0E0*I_ESP_K6xz_Fx2y_aa-2.0E0*5*I_ESP_H4xz_Fx2y_a;
    abcd[iGrid*1260+484] = 4.0E0*I_ESP_K5xyz_Fx2y_aa-2.0E0*4*I_ESP_H3xyz_Fx2y_a;
    abcd[iGrid*1260+485] = 4.0E0*I_ESP_K5x2z_Fx2y_aa-2.0E0*1*I_ESP_H5x_Fx2y_a-2.0E0*4*I_ESP_H3x2z_Fx2y_a+4*1*I_ESP_F3x_Fx2y;
    abcd[iGrid*1260+486] = 4.0E0*I_ESP_K4x2yz_Fx2y_aa-2.0E0*3*I_ESP_H2x2yz_Fx2y_a;
    abcd[iGrid*1260+487] = 4.0E0*I_ESP_K4xy2z_Fx2y_aa-2.0E0*1*I_ESP_H4xy_Fx2y_a-2.0E0*3*I_ESP_H2xy2z_Fx2y_a+3*1*I_ESP_F2xy_Fx2y;
    abcd[iGrid*1260+488] = 4.0E0*I_ESP_K4x3z_Fx2y_aa-2.0E0*2*I_ESP_H4xz_Fx2y_a-2.0E0*3*I_ESP_H2x3z_Fx2y_a+3*2*I_ESP_F2xz_Fx2y;
    abcd[iGrid*1260+489] = 4.0E0*I_ESP_K3x3yz_Fx2y_aa-2.0E0*2*I_ESP_Hx3yz_Fx2y_a;
    abcd[iGrid*1260+490] = 4.0E0*I_ESP_K3x2y2z_Fx2y_aa-2.0E0*1*I_ESP_H3x2y_Fx2y_a-2.0E0*2*I_ESP_Hx2y2z_Fx2y_a+2*1*I_ESP_Fx2y_Fx2y;
    abcd[iGrid*1260+491] = 4.0E0*I_ESP_K3xy3z_Fx2y_aa-2.0E0*2*I_ESP_H3xyz_Fx2y_a-2.0E0*2*I_ESP_Hxy3z_Fx2y_a+2*2*I_ESP_Fxyz_Fx2y;
    abcd[iGrid*1260+492] = 4.0E0*I_ESP_K3x4z_Fx2y_aa-2.0E0*3*I_ESP_H3x2z_Fx2y_a-2.0E0*2*I_ESP_Hx4z_Fx2y_a+2*3*I_ESP_Fx2z_Fx2y;
    abcd[iGrid*1260+493] = 4.0E0*I_ESP_K2x4yz_Fx2y_aa-2.0E0*1*I_ESP_H4yz_Fx2y_a;
    abcd[iGrid*1260+494] = 4.0E0*I_ESP_K2x3y2z_Fx2y_aa-2.0E0*1*I_ESP_H2x3y_Fx2y_a-2.0E0*1*I_ESP_H3y2z_Fx2y_a+1*I_ESP_F3y_Fx2y;
    abcd[iGrid*1260+495] = 4.0E0*I_ESP_K2x2y3z_Fx2y_aa-2.0E0*2*I_ESP_H2x2yz_Fx2y_a-2.0E0*1*I_ESP_H2y3z_Fx2y_a+2*I_ESP_F2yz_Fx2y;
    abcd[iGrid*1260+496] = 4.0E0*I_ESP_K2xy4z_Fx2y_aa-2.0E0*3*I_ESP_H2xy2z_Fx2y_a-2.0E0*1*I_ESP_Hy4z_Fx2y_a+3*I_ESP_Fy2z_Fx2y;
    abcd[iGrid*1260+497] = 4.0E0*I_ESP_K2x5z_Fx2y_aa-2.0E0*4*I_ESP_H2x3z_Fx2y_a-2.0E0*1*I_ESP_H5z_Fx2y_a+4*I_ESP_F3z_Fx2y;
    abcd[iGrid*1260+498] = 4.0E0*I_ESP_Kx5yz_Fx2y_aa;
    abcd[iGrid*1260+499] = 4.0E0*I_ESP_Kx4y2z_Fx2y_aa-2.0E0*1*I_ESP_Hx4y_Fx2y_a;
    abcd[iGrid*1260+500] = 4.0E0*I_ESP_Kx3y3z_Fx2y_aa-2.0E0*2*I_ESP_Hx3yz_Fx2y_a;
    abcd[iGrid*1260+501] = 4.0E0*I_ESP_Kx2y4z_Fx2y_aa-2.0E0*3*I_ESP_Hx2y2z_Fx2y_a;
    abcd[iGrid*1260+502] = 4.0E0*I_ESP_Kxy5z_Fx2y_aa-2.0E0*4*I_ESP_Hxy3z_Fx2y_a;
    abcd[iGrid*1260+503] = 4.0E0*I_ESP_Kx6z_Fx2y_aa-2.0E0*5*I_ESP_Hx4z_Fx2y_a;
    abcd[iGrid*1260+504] = 4.0E0*I_ESP_K6xz_Fxyz_aa-2.0E0*5*I_ESP_H4xz_Fxyz_a;
    abcd[iGrid*1260+505] = 4.0E0*I_ESP_K5xyz_Fxyz_aa-2.0E0*4*I_ESP_H3xyz_Fxyz_a;
    abcd[iGrid*1260+506] = 4.0E0*I_ESP_K5x2z_Fxyz_aa-2.0E0*1*I_ESP_H5x_Fxyz_a-2.0E0*4*I_ESP_H3x2z_Fxyz_a+4*1*I_ESP_F3x_Fxyz;
    abcd[iGrid*1260+507] = 4.0E0*I_ESP_K4x2yz_Fxyz_aa-2.0E0*3*I_ESP_H2x2yz_Fxyz_a;
    abcd[iGrid*1260+508] = 4.0E0*I_ESP_K4xy2z_Fxyz_aa-2.0E0*1*I_ESP_H4xy_Fxyz_a-2.0E0*3*I_ESP_H2xy2z_Fxyz_a+3*1*I_ESP_F2xy_Fxyz;
    abcd[iGrid*1260+509] = 4.0E0*I_ESP_K4x3z_Fxyz_aa-2.0E0*2*I_ESP_H4xz_Fxyz_a-2.0E0*3*I_ESP_H2x3z_Fxyz_a+3*2*I_ESP_F2xz_Fxyz;
    abcd[iGrid*1260+510] = 4.0E0*I_ESP_K3x3yz_Fxyz_aa-2.0E0*2*I_ESP_Hx3yz_Fxyz_a;
    abcd[iGrid*1260+511] = 4.0E0*I_ESP_K3x2y2z_Fxyz_aa-2.0E0*1*I_ESP_H3x2y_Fxyz_a-2.0E0*2*I_ESP_Hx2y2z_Fxyz_a+2*1*I_ESP_Fx2y_Fxyz;
    abcd[iGrid*1260+512] = 4.0E0*I_ESP_K3xy3z_Fxyz_aa-2.0E0*2*I_ESP_H3xyz_Fxyz_a-2.0E0*2*I_ESP_Hxy3z_Fxyz_a+2*2*I_ESP_Fxyz_Fxyz;
    abcd[iGrid*1260+513] = 4.0E0*I_ESP_K3x4z_Fxyz_aa-2.0E0*3*I_ESP_H3x2z_Fxyz_a-2.0E0*2*I_ESP_Hx4z_Fxyz_a+2*3*I_ESP_Fx2z_Fxyz;
    abcd[iGrid*1260+514] = 4.0E0*I_ESP_K2x4yz_Fxyz_aa-2.0E0*1*I_ESP_H4yz_Fxyz_a;
    abcd[iGrid*1260+515] = 4.0E0*I_ESP_K2x3y2z_Fxyz_aa-2.0E0*1*I_ESP_H2x3y_Fxyz_a-2.0E0*1*I_ESP_H3y2z_Fxyz_a+1*I_ESP_F3y_Fxyz;
    abcd[iGrid*1260+516] = 4.0E0*I_ESP_K2x2y3z_Fxyz_aa-2.0E0*2*I_ESP_H2x2yz_Fxyz_a-2.0E0*1*I_ESP_H2y3z_Fxyz_a+2*I_ESP_F2yz_Fxyz;
    abcd[iGrid*1260+517] = 4.0E0*I_ESP_K2xy4z_Fxyz_aa-2.0E0*3*I_ESP_H2xy2z_Fxyz_a-2.0E0*1*I_ESP_Hy4z_Fxyz_a+3*I_ESP_Fy2z_Fxyz;
    abcd[iGrid*1260+518] = 4.0E0*I_ESP_K2x5z_Fxyz_aa-2.0E0*4*I_ESP_H2x3z_Fxyz_a-2.0E0*1*I_ESP_H5z_Fxyz_a+4*I_ESP_F3z_Fxyz;
    abcd[iGrid*1260+519] = 4.0E0*I_ESP_Kx5yz_Fxyz_aa;
    abcd[iGrid*1260+520] = 4.0E0*I_ESP_Kx4y2z_Fxyz_aa-2.0E0*1*I_ESP_Hx4y_Fxyz_a;
    abcd[iGrid*1260+521] = 4.0E0*I_ESP_Kx3y3z_Fxyz_aa-2.0E0*2*I_ESP_Hx3yz_Fxyz_a;
    abcd[iGrid*1260+522] = 4.0E0*I_ESP_Kx2y4z_Fxyz_aa-2.0E0*3*I_ESP_Hx2y2z_Fxyz_a;
    abcd[iGrid*1260+523] = 4.0E0*I_ESP_Kxy5z_Fxyz_aa-2.0E0*4*I_ESP_Hxy3z_Fxyz_a;
    abcd[iGrid*1260+524] = 4.0E0*I_ESP_Kx6z_Fxyz_aa-2.0E0*5*I_ESP_Hx4z_Fxyz_a;
    abcd[iGrid*1260+525] = 4.0E0*I_ESP_K6xz_Fx2z_aa-2.0E0*5*I_ESP_H4xz_Fx2z_a;
    abcd[iGrid*1260+526] = 4.0E0*I_ESP_K5xyz_Fx2z_aa-2.0E0*4*I_ESP_H3xyz_Fx2z_a;
    abcd[iGrid*1260+527] = 4.0E0*I_ESP_K5x2z_Fx2z_aa-2.0E0*1*I_ESP_H5x_Fx2z_a-2.0E0*4*I_ESP_H3x2z_Fx2z_a+4*1*I_ESP_F3x_Fx2z;
    abcd[iGrid*1260+528] = 4.0E0*I_ESP_K4x2yz_Fx2z_aa-2.0E0*3*I_ESP_H2x2yz_Fx2z_a;
    abcd[iGrid*1260+529] = 4.0E0*I_ESP_K4xy2z_Fx2z_aa-2.0E0*1*I_ESP_H4xy_Fx2z_a-2.0E0*3*I_ESP_H2xy2z_Fx2z_a+3*1*I_ESP_F2xy_Fx2z;
    abcd[iGrid*1260+530] = 4.0E0*I_ESP_K4x3z_Fx2z_aa-2.0E0*2*I_ESP_H4xz_Fx2z_a-2.0E0*3*I_ESP_H2x3z_Fx2z_a+3*2*I_ESP_F2xz_Fx2z;
    abcd[iGrid*1260+531] = 4.0E0*I_ESP_K3x3yz_Fx2z_aa-2.0E0*2*I_ESP_Hx3yz_Fx2z_a;
    abcd[iGrid*1260+532] = 4.0E0*I_ESP_K3x2y2z_Fx2z_aa-2.0E0*1*I_ESP_H3x2y_Fx2z_a-2.0E0*2*I_ESP_Hx2y2z_Fx2z_a+2*1*I_ESP_Fx2y_Fx2z;
    abcd[iGrid*1260+533] = 4.0E0*I_ESP_K3xy3z_Fx2z_aa-2.0E0*2*I_ESP_H3xyz_Fx2z_a-2.0E0*2*I_ESP_Hxy3z_Fx2z_a+2*2*I_ESP_Fxyz_Fx2z;
    abcd[iGrid*1260+534] = 4.0E0*I_ESP_K3x4z_Fx2z_aa-2.0E0*3*I_ESP_H3x2z_Fx2z_a-2.0E0*2*I_ESP_Hx4z_Fx2z_a+2*3*I_ESP_Fx2z_Fx2z;
    abcd[iGrid*1260+535] = 4.0E0*I_ESP_K2x4yz_Fx2z_aa-2.0E0*1*I_ESP_H4yz_Fx2z_a;
    abcd[iGrid*1260+536] = 4.0E0*I_ESP_K2x3y2z_Fx2z_aa-2.0E0*1*I_ESP_H2x3y_Fx2z_a-2.0E0*1*I_ESP_H3y2z_Fx2z_a+1*I_ESP_F3y_Fx2z;
    abcd[iGrid*1260+537] = 4.0E0*I_ESP_K2x2y3z_Fx2z_aa-2.0E0*2*I_ESP_H2x2yz_Fx2z_a-2.0E0*1*I_ESP_H2y3z_Fx2z_a+2*I_ESP_F2yz_Fx2z;
    abcd[iGrid*1260+538] = 4.0E0*I_ESP_K2xy4z_Fx2z_aa-2.0E0*3*I_ESP_H2xy2z_Fx2z_a-2.0E0*1*I_ESP_Hy4z_Fx2z_a+3*I_ESP_Fy2z_Fx2z;
    abcd[iGrid*1260+539] = 4.0E0*I_ESP_K2x5z_Fx2z_aa-2.0E0*4*I_ESP_H2x3z_Fx2z_a-2.0E0*1*I_ESP_H5z_Fx2z_a+4*I_ESP_F3z_Fx2z;
    abcd[iGrid*1260+540] = 4.0E0*I_ESP_Kx5yz_Fx2z_aa;
    abcd[iGrid*1260+541] = 4.0E0*I_ESP_Kx4y2z_Fx2z_aa-2.0E0*1*I_ESP_Hx4y_Fx2z_a;
    abcd[iGrid*1260+542] = 4.0E0*I_ESP_Kx3y3z_Fx2z_aa-2.0E0*2*I_ESP_Hx3yz_Fx2z_a;
    abcd[iGrid*1260+543] = 4.0E0*I_ESP_Kx2y4z_Fx2z_aa-2.0E0*3*I_ESP_Hx2y2z_Fx2z_a;
    abcd[iGrid*1260+544] = 4.0E0*I_ESP_Kxy5z_Fx2z_aa-2.0E0*4*I_ESP_Hxy3z_Fx2z_a;
    abcd[iGrid*1260+545] = 4.0E0*I_ESP_Kx6z_Fx2z_aa-2.0E0*5*I_ESP_Hx4z_Fx2z_a;
    abcd[iGrid*1260+546] = 4.0E0*I_ESP_K6xz_F3y_aa-2.0E0*5*I_ESP_H4xz_F3y_a;
    abcd[iGrid*1260+547] = 4.0E0*I_ESP_K5xyz_F3y_aa-2.0E0*4*I_ESP_H3xyz_F3y_a;
    abcd[iGrid*1260+548] = 4.0E0*I_ESP_K5x2z_F3y_aa-2.0E0*1*I_ESP_H5x_F3y_a-2.0E0*4*I_ESP_H3x2z_F3y_a+4*1*I_ESP_F3x_F3y;
    abcd[iGrid*1260+549] = 4.0E0*I_ESP_K4x2yz_F3y_aa-2.0E0*3*I_ESP_H2x2yz_F3y_a;
    abcd[iGrid*1260+550] = 4.0E0*I_ESP_K4xy2z_F3y_aa-2.0E0*1*I_ESP_H4xy_F3y_a-2.0E0*3*I_ESP_H2xy2z_F3y_a+3*1*I_ESP_F2xy_F3y;
    abcd[iGrid*1260+551] = 4.0E0*I_ESP_K4x3z_F3y_aa-2.0E0*2*I_ESP_H4xz_F3y_a-2.0E0*3*I_ESP_H2x3z_F3y_a+3*2*I_ESP_F2xz_F3y;
    abcd[iGrid*1260+552] = 4.0E0*I_ESP_K3x3yz_F3y_aa-2.0E0*2*I_ESP_Hx3yz_F3y_a;
    abcd[iGrid*1260+553] = 4.0E0*I_ESP_K3x2y2z_F3y_aa-2.0E0*1*I_ESP_H3x2y_F3y_a-2.0E0*2*I_ESP_Hx2y2z_F3y_a+2*1*I_ESP_Fx2y_F3y;
    abcd[iGrid*1260+554] = 4.0E0*I_ESP_K3xy3z_F3y_aa-2.0E0*2*I_ESP_H3xyz_F3y_a-2.0E0*2*I_ESP_Hxy3z_F3y_a+2*2*I_ESP_Fxyz_F3y;
    abcd[iGrid*1260+555] = 4.0E0*I_ESP_K3x4z_F3y_aa-2.0E0*3*I_ESP_H3x2z_F3y_a-2.0E0*2*I_ESP_Hx4z_F3y_a+2*3*I_ESP_Fx2z_F3y;
    abcd[iGrid*1260+556] = 4.0E0*I_ESP_K2x4yz_F3y_aa-2.0E0*1*I_ESP_H4yz_F3y_a;
    abcd[iGrid*1260+557] = 4.0E0*I_ESP_K2x3y2z_F3y_aa-2.0E0*1*I_ESP_H2x3y_F3y_a-2.0E0*1*I_ESP_H3y2z_F3y_a+1*I_ESP_F3y_F3y;
    abcd[iGrid*1260+558] = 4.0E0*I_ESP_K2x2y3z_F3y_aa-2.0E0*2*I_ESP_H2x2yz_F3y_a-2.0E0*1*I_ESP_H2y3z_F3y_a+2*I_ESP_F2yz_F3y;
    abcd[iGrid*1260+559] = 4.0E0*I_ESP_K2xy4z_F3y_aa-2.0E0*3*I_ESP_H2xy2z_F3y_a-2.0E0*1*I_ESP_Hy4z_F3y_a+3*I_ESP_Fy2z_F3y;
    abcd[iGrid*1260+560] = 4.0E0*I_ESP_K2x5z_F3y_aa-2.0E0*4*I_ESP_H2x3z_F3y_a-2.0E0*1*I_ESP_H5z_F3y_a+4*I_ESP_F3z_F3y;
    abcd[iGrid*1260+561] = 4.0E0*I_ESP_Kx5yz_F3y_aa;
    abcd[iGrid*1260+562] = 4.0E0*I_ESP_Kx4y2z_F3y_aa-2.0E0*1*I_ESP_Hx4y_F3y_a;
    abcd[iGrid*1260+563] = 4.0E0*I_ESP_Kx3y3z_F3y_aa-2.0E0*2*I_ESP_Hx3yz_F3y_a;
    abcd[iGrid*1260+564] = 4.0E0*I_ESP_Kx2y4z_F3y_aa-2.0E0*3*I_ESP_Hx2y2z_F3y_a;
    abcd[iGrid*1260+565] = 4.0E0*I_ESP_Kxy5z_F3y_aa-2.0E0*4*I_ESP_Hxy3z_F3y_a;
    abcd[iGrid*1260+566] = 4.0E0*I_ESP_Kx6z_F3y_aa-2.0E0*5*I_ESP_Hx4z_F3y_a;
    abcd[iGrid*1260+567] = 4.0E0*I_ESP_K6xz_F2yz_aa-2.0E0*5*I_ESP_H4xz_F2yz_a;
    abcd[iGrid*1260+568] = 4.0E0*I_ESP_K5xyz_F2yz_aa-2.0E0*4*I_ESP_H3xyz_F2yz_a;
    abcd[iGrid*1260+569] = 4.0E0*I_ESP_K5x2z_F2yz_aa-2.0E0*1*I_ESP_H5x_F2yz_a-2.0E0*4*I_ESP_H3x2z_F2yz_a+4*1*I_ESP_F3x_F2yz;
    abcd[iGrid*1260+570] = 4.0E0*I_ESP_K4x2yz_F2yz_aa-2.0E0*3*I_ESP_H2x2yz_F2yz_a;
    abcd[iGrid*1260+571] = 4.0E0*I_ESP_K4xy2z_F2yz_aa-2.0E0*1*I_ESP_H4xy_F2yz_a-2.0E0*3*I_ESP_H2xy2z_F2yz_a+3*1*I_ESP_F2xy_F2yz;
    abcd[iGrid*1260+572] = 4.0E0*I_ESP_K4x3z_F2yz_aa-2.0E0*2*I_ESP_H4xz_F2yz_a-2.0E0*3*I_ESP_H2x3z_F2yz_a+3*2*I_ESP_F2xz_F2yz;
    abcd[iGrid*1260+573] = 4.0E0*I_ESP_K3x3yz_F2yz_aa-2.0E0*2*I_ESP_Hx3yz_F2yz_a;
    abcd[iGrid*1260+574] = 4.0E0*I_ESP_K3x2y2z_F2yz_aa-2.0E0*1*I_ESP_H3x2y_F2yz_a-2.0E0*2*I_ESP_Hx2y2z_F2yz_a+2*1*I_ESP_Fx2y_F2yz;
    abcd[iGrid*1260+575] = 4.0E0*I_ESP_K3xy3z_F2yz_aa-2.0E0*2*I_ESP_H3xyz_F2yz_a-2.0E0*2*I_ESP_Hxy3z_F2yz_a+2*2*I_ESP_Fxyz_F2yz;
    abcd[iGrid*1260+576] = 4.0E0*I_ESP_K3x4z_F2yz_aa-2.0E0*3*I_ESP_H3x2z_F2yz_a-2.0E0*2*I_ESP_Hx4z_F2yz_a+2*3*I_ESP_Fx2z_F2yz;
    abcd[iGrid*1260+577] = 4.0E0*I_ESP_K2x4yz_F2yz_aa-2.0E0*1*I_ESP_H4yz_F2yz_a;
    abcd[iGrid*1260+578] = 4.0E0*I_ESP_K2x3y2z_F2yz_aa-2.0E0*1*I_ESP_H2x3y_F2yz_a-2.0E0*1*I_ESP_H3y2z_F2yz_a+1*I_ESP_F3y_F2yz;
    abcd[iGrid*1260+579] = 4.0E0*I_ESP_K2x2y3z_F2yz_aa-2.0E0*2*I_ESP_H2x2yz_F2yz_a-2.0E0*1*I_ESP_H2y3z_F2yz_a+2*I_ESP_F2yz_F2yz;
    abcd[iGrid*1260+580] = 4.0E0*I_ESP_K2xy4z_F2yz_aa-2.0E0*3*I_ESP_H2xy2z_F2yz_a-2.0E0*1*I_ESP_Hy4z_F2yz_a+3*I_ESP_Fy2z_F2yz;
    abcd[iGrid*1260+581] = 4.0E0*I_ESP_K2x5z_F2yz_aa-2.0E0*4*I_ESP_H2x3z_F2yz_a-2.0E0*1*I_ESP_H5z_F2yz_a+4*I_ESP_F3z_F2yz;
    abcd[iGrid*1260+582] = 4.0E0*I_ESP_Kx5yz_F2yz_aa;
    abcd[iGrid*1260+583] = 4.0E0*I_ESP_Kx4y2z_F2yz_aa-2.0E0*1*I_ESP_Hx4y_F2yz_a;
    abcd[iGrid*1260+584] = 4.0E0*I_ESP_Kx3y3z_F2yz_aa-2.0E0*2*I_ESP_Hx3yz_F2yz_a;
    abcd[iGrid*1260+585] = 4.0E0*I_ESP_Kx2y4z_F2yz_aa-2.0E0*3*I_ESP_Hx2y2z_F2yz_a;
    abcd[iGrid*1260+586] = 4.0E0*I_ESP_Kxy5z_F2yz_aa-2.0E0*4*I_ESP_Hxy3z_F2yz_a;
    abcd[iGrid*1260+587] = 4.0E0*I_ESP_Kx6z_F2yz_aa-2.0E0*5*I_ESP_Hx4z_F2yz_a;
    abcd[iGrid*1260+588] = 4.0E0*I_ESP_K6xz_Fy2z_aa-2.0E0*5*I_ESP_H4xz_Fy2z_a;
    abcd[iGrid*1260+589] = 4.0E0*I_ESP_K5xyz_Fy2z_aa-2.0E0*4*I_ESP_H3xyz_Fy2z_a;
    abcd[iGrid*1260+590] = 4.0E0*I_ESP_K5x2z_Fy2z_aa-2.0E0*1*I_ESP_H5x_Fy2z_a-2.0E0*4*I_ESP_H3x2z_Fy2z_a+4*1*I_ESP_F3x_Fy2z;
    abcd[iGrid*1260+591] = 4.0E0*I_ESP_K4x2yz_Fy2z_aa-2.0E0*3*I_ESP_H2x2yz_Fy2z_a;
    abcd[iGrid*1260+592] = 4.0E0*I_ESP_K4xy2z_Fy2z_aa-2.0E0*1*I_ESP_H4xy_Fy2z_a-2.0E0*3*I_ESP_H2xy2z_Fy2z_a+3*1*I_ESP_F2xy_Fy2z;
    abcd[iGrid*1260+593] = 4.0E0*I_ESP_K4x3z_Fy2z_aa-2.0E0*2*I_ESP_H4xz_Fy2z_a-2.0E0*3*I_ESP_H2x3z_Fy2z_a+3*2*I_ESP_F2xz_Fy2z;
    abcd[iGrid*1260+594] = 4.0E0*I_ESP_K3x3yz_Fy2z_aa-2.0E0*2*I_ESP_Hx3yz_Fy2z_a;
    abcd[iGrid*1260+595] = 4.0E0*I_ESP_K3x2y2z_Fy2z_aa-2.0E0*1*I_ESP_H3x2y_Fy2z_a-2.0E0*2*I_ESP_Hx2y2z_Fy2z_a+2*1*I_ESP_Fx2y_Fy2z;
    abcd[iGrid*1260+596] = 4.0E0*I_ESP_K3xy3z_Fy2z_aa-2.0E0*2*I_ESP_H3xyz_Fy2z_a-2.0E0*2*I_ESP_Hxy3z_Fy2z_a+2*2*I_ESP_Fxyz_Fy2z;
    abcd[iGrid*1260+597] = 4.0E0*I_ESP_K3x4z_Fy2z_aa-2.0E0*3*I_ESP_H3x2z_Fy2z_a-2.0E0*2*I_ESP_Hx4z_Fy2z_a+2*3*I_ESP_Fx2z_Fy2z;
    abcd[iGrid*1260+598] = 4.0E0*I_ESP_K2x4yz_Fy2z_aa-2.0E0*1*I_ESP_H4yz_Fy2z_a;
    abcd[iGrid*1260+599] = 4.0E0*I_ESP_K2x3y2z_Fy2z_aa-2.0E0*1*I_ESP_H2x3y_Fy2z_a-2.0E0*1*I_ESP_H3y2z_Fy2z_a+1*I_ESP_F3y_Fy2z;
    abcd[iGrid*1260+600] = 4.0E0*I_ESP_K2x2y3z_Fy2z_aa-2.0E0*2*I_ESP_H2x2yz_Fy2z_a-2.0E0*1*I_ESP_H2y3z_Fy2z_a+2*I_ESP_F2yz_Fy2z;
    abcd[iGrid*1260+601] = 4.0E0*I_ESP_K2xy4z_Fy2z_aa-2.0E0*3*I_ESP_H2xy2z_Fy2z_a-2.0E0*1*I_ESP_Hy4z_Fy2z_a+3*I_ESP_Fy2z_Fy2z;
    abcd[iGrid*1260+602] = 4.0E0*I_ESP_K2x5z_Fy2z_aa-2.0E0*4*I_ESP_H2x3z_Fy2z_a-2.0E0*1*I_ESP_H5z_Fy2z_a+4*I_ESP_F3z_Fy2z;
    abcd[iGrid*1260+603] = 4.0E0*I_ESP_Kx5yz_Fy2z_aa;
    abcd[iGrid*1260+604] = 4.0E0*I_ESP_Kx4y2z_Fy2z_aa-2.0E0*1*I_ESP_Hx4y_Fy2z_a;
    abcd[iGrid*1260+605] = 4.0E0*I_ESP_Kx3y3z_Fy2z_aa-2.0E0*2*I_ESP_Hx3yz_Fy2z_a;
    abcd[iGrid*1260+606] = 4.0E0*I_ESP_Kx2y4z_Fy2z_aa-2.0E0*3*I_ESP_Hx2y2z_Fy2z_a;
    abcd[iGrid*1260+607] = 4.0E0*I_ESP_Kxy5z_Fy2z_aa-2.0E0*4*I_ESP_Hxy3z_Fy2z_a;
    abcd[iGrid*1260+608] = 4.0E0*I_ESP_Kx6z_Fy2z_aa-2.0E0*5*I_ESP_Hx4z_Fy2z_a;
    abcd[iGrid*1260+609] = 4.0E0*I_ESP_K6xz_F3z_aa-2.0E0*5*I_ESP_H4xz_F3z_a;
    abcd[iGrid*1260+610] = 4.0E0*I_ESP_K5xyz_F3z_aa-2.0E0*4*I_ESP_H3xyz_F3z_a;
    abcd[iGrid*1260+611] = 4.0E0*I_ESP_K5x2z_F3z_aa-2.0E0*1*I_ESP_H5x_F3z_a-2.0E0*4*I_ESP_H3x2z_F3z_a+4*1*I_ESP_F3x_F3z;
    abcd[iGrid*1260+612] = 4.0E0*I_ESP_K4x2yz_F3z_aa-2.0E0*3*I_ESP_H2x2yz_F3z_a;
    abcd[iGrid*1260+613] = 4.0E0*I_ESP_K4xy2z_F3z_aa-2.0E0*1*I_ESP_H4xy_F3z_a-2.0E0*3*I_ESP_H2xy2z_F3z_a+3*1*I_ESP_F2xy_F3z;
    abcd[iGrid*1260+614] = 4.0E0*I_ESP_K4x3z_F3z_aa-2.0E0*2*I_ESP_H4xz_F3z_a-2.0E0*3*I_ESP_H2x3z_F3z_a+3*2*I_ESP_F2xz_F3z;
    abcd[iGrid*1260+615] = 4.0E0*I_ESP_K3x3yz_F3z_aa-2.0E0*2*I_ESP_Hx3yz_F3z_a;
    abcd[iGrid*1260+616] = 4.0E0*I_ESP_K3x2y2z_F3z_aa-2.0E0*1*I_ESP_H3x2y_F3z_a-2.0E0*2*I_ESP_Hx2y2z_F3z_a+2*1*I_ESP_Fx2y_F3z;
    abcd[iGrid*1260+617] = 4.0E0*I_ESP_K3xy3z_F3z_aa-2.0E0*2*I_ESP_H3xyz_F3z_a-2.0E0*2*I_ESP_Hxy3z_F3z_a+2*2*I_ESP_Fxyz_F3z;
    abcd[iGrid*1260+618] = 4.0E0*I_ESP_K3x4z_F3z_aa-2.0E0*3*I_ESP_H3x2z_F3z_a-2.0E0*2*I_ESP_Hx4z_F3z_a+2*3*I_ESP_Fx2z_F3z;
    abcd[iGrid*1260+619] = 4.0E0*I_ESP_K2x4yz_F3z_aa-2.0E0*1*I_ESP_H4yz_F3z_a;
    abcd[iGrid*1260+620] = 4.0E0*I_ESP_K2x3y2z_F3z_aa-2.0E0*1*I_ESP_H2x3y_F3z_a-2.0E0*1*I_ESP_H3y2z_F3z_a+1*I_ESP_F3y_F3z;
    abcd[iGrid*1260+621] = 4.0E0*I_ESP_K2x2y3z_F3z_aa-2.0E0*2*I_ESP_H2x2yz_F3z_a-2.0E0*1*I_ESP_H2y3z_F3z_a+2*I_ESP_F2yz_F3z;
    abcd[iGrid*1260+622] = 4.0E0*I_ESP_K2xy4z_F3z_aa-2.0E0*3*I_ESP_H2xy2z_F3z_a-2.0E0*1*I_ESP_Hy4z_F3z_a+3*I_ESP_Fy2z_F3z;
    abcd[iGrid*1260+623] = 4.0E0*I_ESP_K2x5z_F3z_aa-2.0E0*4*I_ESP_H2x3z_F3z_a-2.0E0*1*I_ESP_H5z_F3z_a+4*I_ESP_F3z_F3z;
    abcd[iGrid*1260+624] = 4.0E0*I_ESP_Kx5yz_F3z_aa;
    abcd[iGrid*1260+625] = 4.0E0*I_ESP_Kx4y2z_F3z_aa-2.0E0*1*I_ESP_Hx4y_F3z_a;
    abcd[iGrid*1260+626] = 4.0E0*I_ESP_Kx3y3z_F3z_aa-2.0E0*2*I_ESP_Hx3yz_F3z_a;
    abcd[iGrid*1260+627] = 4.0E0*I_ESP_Kx2y4z_F3z_aa-2.0E0*3*I_ESP_Hx2y2z_F3z_a;
    abcd[iGrid*1260+628] = 4.0E0*I_ESP_Kxy5z_F3z_aa-2.0E0*4*I_ESP_Hxy3z_F3z_a;
    abcd[iGrid*1260+629] = 4.0E0*I_ESP_Kx6z_F3z_aa-2.0E0*5*I_ESP_Hx4z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_F_aa
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    abcd[iGrid*1260+630] = 4.0E0*I_ESP_K5x2y_F3x_aa-2.0E0*1*I_ESP_H5x_F3x_a;
    abcd[iGrid*1260+631] = 4.0E0*I_ESP_K4x3y_F3x_aa-2.0E0*1*I_ESP_H4xy_F3x_a-2.0E0*2*I_ESP_H4xy_F3x_a;
    abcd[iGrid*1260+632] = 4.0E0*I_ESP_K4x2yz_F3x_aa-2.0E0*1*I_ESP_H4xz_F3x_a;
    abcd[iGrid*1260+633] = 4.0E0*I_ESP_K3x4y_F3x_aa-2.0E0*2*I_ESP_H3x2y_F3x_a-2.0E0*3*I_ESP_H3x2y_F3x_a+2*1*I_ESP_F3x_F3x;
    abcd[iGrid*1260+634] = 4.0E0*I_ESP_K3x3yz_F3x_aa-2.0E0*1*I_ESP_H3xyz_F3x_a-2.0E0*2*I_ESP_H3xyz_F3x_a;
    abcd[iGrid*1260+635] = 4.0E0*I_ESP_K3x2y2z_F3x_aa-2.0E0*1*I_ESP_H3x2z_F3x_a;
    abcd[iGrid*1260+636] = 4.0E0*I_ESP_K2x5y_F3x_aa-2.0E0*3*I_ESP_H2x3y_F3x_a-2.0E0*4*I_ESP_H2x3y_F3x_a+3*2*I_ESP_F2xy_F3x;
    abcd[iGrid*1260+637] = 4.0E0*I_ESP_K2x4yz_F3x_aa-2.0E0*2*I_ESP_H2x2yz_F3x_a-2.0E0*3*I_ESP_H2x2yz_F3x_a+2*1*I_ESP_F2xz_F3x;
    abcd[iGrid*1260+638] = 4.0E0*I_ESP_K2x3y2z_F3x_aa-2.0E0*1*I_ESP_H2xy2z_F3x_a-2.0E0*2*I_ESP_H2xy2z_F3x_a;
    abcd[iGrid*1260+639] = 4.0E0*I_ESP_K2x2y3z_F3x_aa-2.0E0*1*I_ESP_H2x3z_F3x_a;
    abcd[iGrid*1260+640] = 4.0E0*I_ESP_Kx6y_F3x_aa-2.0E0*4*I_ESP_Hx4y_F3x_a-2.0E0*5*I_ESP_Hx4y_F3x_a+4*3*I_ESP_Fx2y_F3x;
    abcd[iGrid*1260+641] = 4.0E0*I_ESP_Kx5yz_F3x_aa-2.0E0*3*I_ESP_Hx3yz_F3x_a-2.0E0*4*I_ESP_Hx3yz_F3x_a+3*2*I_ESP_Fxyz_F3x;
    abcd[iGrid*1260+642] = 4.0E0*I_ESP_Kx4y2z_F3x_aa-2.0E0*2*I_ESP_Hx2y2z_F3x_a-2.0E0*3*I_ESP_Hx2y2z_F3x_a+2*1*I_ESP_Fx2z_F3x;
    abcd[iGrid*1260+643] = 4.0E0*I_ESP_Kx3y3z_F3x_aa-2.0E0*1*I_ESP_Hxy3z_F3x_a-2.0E0*2*I_ESP_Hxy3z_F3x_a;
    abcd[iGrid*1260+644] = 4.0E0*I_ESP_Kx2y4z_F3x_aa-2.0E0*1*I_ESP_Hx4z_F3x_a;
    abcd[iGrid*1260+645] = 4.0E0*I_ESP_K7y_F3x_aa-2.0E0*5*I_ESP_H5y_F3x_a-2.0E0*6*I_ESP_H5y_F3x_a+5*4*I_ESP_F3y_F3x;
    abcd[iGrid*1260+646] = 4.0E0*I_ESP_K6yz_F3x_aa-2.0E0*4*I_ESP_H4yz_F3x_a-2.0E0*5*I_ESP_H4yz_F3x_a+4*3*I_ESP_F2yz_F3x;
    abcd[iGrid*1260+647] = 4.0E0*I_ESP_K5y2z_F3x_aa-2.0E0*3*I_ESP_H3y2z_F3x_a-2.0E0*4*I_ESP_H3y2z_F3x_a+3*2*I_ESP_Fy2z_F3x;
    abcd[iGrid*1260+648] = 4.0E0*I_ESP_K4y3z_F3x_aa-2.0E0*2*I_ESP_H2y3z_F3x_a-2.0E0*3*I_ESP_H2y3z_F3x_a+2*1*I_ESP_F3z_F3x;
    abcd[iGrid*1260+649] = 4.0E0*I_ESP_K3y4z_F3x_aa-2.0E0*1*I_ESP_Hy4z_F3x_a-2.0E0*2*I_ESP_Hy4z_F3x_a;
    abcd[iGrid*1260+650] = 4.0E0*I_ESP_K2y5z_F3x_aa-2.0E0*1*I_ESP_H5z_F3x_a;
    abcd[iGrid*1260+651] = 4.0E0*I_ESP_K5x2y_F2xy_aa-2.0E0*1*I_ESP_H5x_F2xy_a;
    abcd[iGrid*1260+652] = 4.0E0*I_ESP_K4x3y_F2xy_aa-2.0E0*1*I_ESP_H4xy_F2xy_a-2.0E0*2*I_ESP_H4xy_F2xy_a;
    abcd[iGrid*1260+653] = 4.0E0*I_ESP_K4x2yz_F2xy_aa-2.0E0*1*I_ESP_H4xz_F2xy_a;
    abcd[iGrid*1260+654] = 4.0E0*I_ESP_K3x4y_F2xy_aa-2.0E0*2*I_ESP_H3x2y_F2xy_a-2.0E0*3*I_ESP_H3x2y_F2xy_a+2*1*I_ESP_F3x_F2xy;
    abcd[iGrid*1260+655] = 4.0E0*I_ESP_K3x3yz_F2xy_aa-2.0E0*1*I_ESP_H3xyz_F2xy_a-2.0E0*2*I_ESP_H3xyz_F2xy_a;
    abcd[iGrid*1260+656] = 4.0E0*I_ESP_K3x2y2z_F2xy_aa-2.0E0*1*I_ESP_H3x2z_F2xy_a;
    abcd[iGrid*1260+657] = 4.0E0*I_ESP_K2x5y_F2xy_aa-2.0E0*3*I_ESP_H2x3y_F2xy_a-2.0E0*4*I_ESP_H2x3y_F2xy_a+3*2*I_ESP_F2xy_F2xy;
    abcd[iGrid*1260+658] = 4.0E0*I_ESP_K2x4yz_F2xy_aa-2.0E0*2*I_ESP_H2x2yz_F2xy_a-2.0E0*3*I_ESP_H2x2yz_F2xy_a+2*1*I_ESP_F2xz_F2xy;
    abcd[iGrid*1260+659] = 4.0E0*I_ESP_K2x3y2z_F2xy_aa-2.0E0*1*I_ESP_H2xy2z_F2xy_a-2.0E0*2*I_ESP_H2xy2z_F2xy_a;
    abcd[iGrid*1260+660] = 4.0E0*I_ESP_K2x2y3z_F2xy_aa-2.0E0*1*I_ESP_H2x3z_F2xy_a;
    abcd[iGrid*1260+661] = 4.0E0*I_ESP_Kx6y_F2xy_aa-2.0E0*4*I_ESP_Hx4y_F2xy_a-2.0E0*5*I_ESP_Hx4y_F2xy_a+4*3*I_ESP_Fx2y_F2xy;
    abcd[iGrid*1260+662] = 4.0E0*I_ESP_Kx5yz_F2xy_aa-2.0E0*3*I_ESP_Hx3yz_F2xy_a-2.0E0*4*I_ESP_Hx3yz_F2xy_a+3*2*I_ESP_Fxyz_F2xy;
    abcd[iGrid*1260+663] = 4.0E0*I_ESP_Kx4y2z_F2xy_aa-2.0E0*2*I_ESP_Hx2y2z_F2xy_a-2.0E0*3*I_ESP_Hx2y2z_F2xy_a+2*1*I_ESP_Fx2z_F2xy;
    abcd[iGrid*1260+664] = 4.0E0*I_ESP_Kx3y3z_F2xy_aa-2.0E0*1*I_ESP_Hxy3z_F2xy_a-2.0E0*2*I_ESP_Hxy3z_F2xy_a;
    abcd[iGrid*1260+665] = 4.0E0*I_ESP_Kx2y4z_F2xy_aa-2.0E0*1*I_ESP_Hx4z_F2xy_a;
    abcd[iGrid*1260+666] = 4.0E0*I_ESP_K7y_F2xy_aa-2.0E0*5*I_ESP_H5y_F2xy_a-2.0E0*6*I_ESP_H5y_F2xy_a+5*4*I_ESP_F3y_F2xy;
    abcd[iGrid*1260+667] = 4.0E0*I_ESP_K6yz_F2xy_aa-2.0E0*4*I_ESP_H4yz_F2xy_a-2.0E0*5*I_ESP_H4yz_F2xy_a+4*3*I_ESP_F2yz_F2xy;
    abcd[iGrid*1260+668] = 4.0E0*I_ESP_K5y2z_F2xy_aa-2.0E0*3*I_ESP_H3y2z_F2xy_a-2.0E0*4*I_ESP_H3y2z_F2xy_a+3*2*I_ESP_Fy2z_F2xy;
    abcd[iGrid*1260+669] = 4.0E0*I_ESP_K4y3z_F2xy_aa-2.0E0*2*I_ESP_H2y3z_F2xy_a-2.0E0*3*I_ESP_H2y3z_F2xy_a+2*1*I_ESP_F3z_F2xy;
    abcd[iGrid*1260+670] = 4.0E0*I_ESP_K3y4z_F2xy_aa-2.0E0*1*I_ESP_Hy4z_F2xy_a-2.0E0*2*I_ESP_Hy4z_F2xy_a;
    abcd[iGrid*1260+671] = 4.0E0*I_ESP_K2y5z_F2xy_aa-2.0E0*1*I_ESP_H5z_F2xy_a;
    abcd[iGrid*1260+672] = 4.0E0*I_ESP_K5x2y_F2xz_aa-2.0E0*1*I_ESP_H5x_F2xz_a;
    abcd[iGrid*1260+673] = 4.0E0*I_ESP_K4x3y_F2xz_aa-2.0E0*1*I_ESP_H4xy_F2xz_a-2.0E0*2*I_ESP_H4xy_F2xz_a;
    abcd[iGrid*1260+674] = 4.0E0*I_ESP_K4x2yz_F2xz_aa-2.0E0*1*I_ESP_H4xz_F2xz_a;
    abcd[iGrid*1260+675] = 4.0E0*I_ESP_K3x4y_F2xz_aa-2.0E0*2*I_ESP_H3x2y_F2xz_a-2.0E0*3*I_ESP_H3x2y_F2xz_a+2*1*I_ESP_F3x_F2xz;
    abcd[iGrid*1260+676] = 4.0E0*I_ESP_K3x3yz_F2xz_aa-2.0E0*1*I_ESP_H3xyz_F2xz_a-2.0E0*2*I_ESP_H3xyz_F2xz_a;
    abcd[iGrid*1260+677] = 4.0E0*I_ESP_K3x2y2z_F2xz_aa-2.0E0*1*I_ESP_H3x2z_F2xz_a;
    abcd[iGrid*1260+678] = 4.0E0*I_ESP_K2x5y_F2xz_aa-2.0E0*3*I_ESP_H2x3y_F2xz_a-2.0E0*4*I_ESP_H2x3y_F2xz_a+3*2*I_ESP_F2xy_F2xz;
    abcd[iGrid*1260+679] = 4.0E0*I_ESP_K2x4yz_F2xz_aa-2.0E0*2*I_ESP_H2x2yz_F2xz_a-2.0E0*3*I_ESP_H2x2yz_F2xz_a+2*1*I_ESP_F2xz_F2xz;
    abcd[iGrid*1260+680] = 4.0E0*I_ESP_K2x3y2z_F2xz_aa-2.0E0*1*I_ESP_H2xy2z_F2xz_a-2.0E0*2*I_ESP_H2xy2z_F2xz_a;
    abcd[iGrid*1260+681] = 4.0E0*I_ESP_K2x2y3z_F2xz_aa-2.0E0*1*I_ESP_H2x3z_F2xz_a;
    abcd[iGrid*1260+682] = 4.0E0*I_ESP_Kx6y_F2xz_aa-2.0E0*4*I_ESP_Hx4y_F2xz_a-2.0E0*5*I_ESP_Hx4y_F2xz_a+4*3*I_ESP_Fx2y_F2xz;
    abcd[iGrid*1260+683] = 4.0E0*I_ESP_Kx5yz_F2xz_aa-2.0E0*3*I_ESP_Hx3yz_F2xz_a-2.0E0*4*I_ESP_Hx3yz_F2xz_a+3*2*I_ESP_Fxyz_F2xz;
    abcd[iGrid*1260+684] = 4.0E0*I_ESP_Kx4y2z_F2xz_aa-2.0E0*2*I_ESP_Hx2y2z_F2xz_a-2.0E0*3*I_ESP_Hx2y2z_F2xz_a+2*1*I_ESP_Fx2z_F2xz;
    abcd[iGrid*1260+685] = 4.0E0*I_ESP_Kx3y3z_F2xz_aa-2.0E0*1*I_ESP_Hxy3z_F2xz_a-2.0E0*2*I_ESP_Hxy3z_F2xz_a;
    abcd[iGrid*1260+686] = 4.0E0*I_ESP_Kx2y4z_F2xz_aa-2.0E0*1*I_ESP_Hx4z_F2xz_a;
    abcd[iGrid*1260+687] = 4.0E0*I_ESP_K7y_F2xz_aa-2.0E0*5*I_ESP_H5y_F2xz_a-2.0E0*6*I_ESP_H5y_F2xz_a+5*4*I_ESP_F3y_F2xz;
    abcd[iGrid*1260+688] = 4.0E0*I_ESP_K6yz_F2xz_aa-2.0E0*4*I_ESP_H4yz_F2xz_a-2.0E0*5*I_ESP_H4yz_F2xz_a+4*3*I_ESP_F2yz_F2xz;
    abcd[iGrid*1260+689] = 4.0E0*I_ESP_K5y2z_F2xz_aa-2.0E0*3*I_ESP_H3y2z_F2xz_a-2.0E0*4*I_ESP_H3y2z_F2xz_a+3*2*I_ESP_Fy2z_F2xz;
    abcd[iGrid*1260+690] = 4.0E0*I_ESP_K4y3z_F2xz_aa-2.0E0*2*I_ESP_H2y3z_F2xz_a-2.0E0*3*I_ESP_H2y3z_F2xz_a+2*1*I_ESP_F3z_F2xz;
    abcd[iGrid*1260+691] = 4.0E0*I_ESP_K3y4z_F2xz_aa-2.0E0*1*I_ESP_Hy4z_F2xz_a-2.0E0*2*I_ESP_Hy4z_F2xz_a;
    abcd[iGrid*1260+692] = 4.0E0*I_ESP_K2y5z_F2xz_aa-2.0E0*1*I_ESP_H5z_F2xz_a;
    abcd[iGrid*1260+693] = 4.0E0*I_ESP_K5x2y_Fx2y_aa-2.0E0*1*I_ESP_H5x_Fx2y_a;
    abcd[iGrid*1260+694] = 4.0E0*I_ESP_K4x3y_Fx2y_aa-2.0E0*1*I_ESP_H4xy_Fx2y_a-2.0E0*2*I_ESP_H4xy_Fx2y_a;
    abcd[iGrid*1260+695] = 4.0E0*I_ESP_K4x2yz_Fx2y_aa-2.0E0*1*I_ESP_H4xz_Fx2y_a;
    abcd[iGrid*1260+696] = 4.0E0*I_ESP_K3x4y_Fx2y_aa-2.0E0*2*I_ESP_H3x2y_Fx2y_a-2.0E0*3*I_ESP_H3x2y_Fx2y_a+2*1*I_ESP_F3x_Fx2y;
    abcd[iGrid*1260+697] = 4.0E0*I_ESP_K3x3yz_Fx2y_aa-2.0E0*1*I_ESP_H3xyz_Fx2y_a-2.0E0*2*I_ESP_H3xyz_Fx2y_a;
    abcd[iGrid*1260+698] = 4.0E0*I_ESP_K3x2y2z_Fx2y_aa-2.0E0*1*I_ESP_H3x2z_Fx2y_a;
    abcd[iGrid*1260+699] = 4.0E0*I_ESP_K2x5y_Fx2y_aa-2.0E0*3*I_ESP_H2x3y_Fx2y_a-2.0E0*4*I_ESP_H2x3y_Fx2y_a+3*2*I_ESP_F2xy_Fx2y;
    abcd[iGrid*1260+700] = 4.0E0*I_ESP_K2x4yz_Fx2y_aa-2.0E0*2*I_ESP_H2x2yz_Fx2y_a-2.0E0*3*I_ESP_H2x2yz_Fx2y_a+2*1*I_ESP_F2xz_Fx2y;
    abcd[iGrid*1260+701] = 4.0E0*I_ESP_K2x3y2z_Fx2y_aa-2.0E0*1*I_ESP_H2xy2z_Fx2y_a-2.0E0*2*I_ESP_H2xy2z_Fx2y_a;
    abcd[iGrid*1260+702] = 4.0E0*I_ESP_K2x2y3z_Fx2y_aa-2.0E0*1*I_ESP_H2x3z_Fx2y_a;
    abcd[iGrid*1260+703] = 4.0E0*I_ESP_Kx6y_Fx2y_aa-2.0E0*4*I_ESP_Hx4y_Fx2y_a-2.0E0*5*I_ESP_Hx4y_Fx2y_a+4*3*I_ESP_Fx2y_Fx2y;
    abcd[iGrid*1260+704] = 4.0E0*I_ESP_Kx5yz_Fx2y_aa-2.0E0*3*I_ESP_Hx3yz_Fx2y_a-2.0E0*4*I_ESP_Hx3yz_Fx2y_a+3*2*I_ESP_Fxyz_Fx2y;
    abcd[iGrid*1260+705] = 4.0E0*I_ESP_Kx4y2z_Fx2y_aa-2.0E0*2*I_ESP_Hx2y2z_Fx2y_a-2.0E0*3*I_ESP_Hx2y2z_Fx2y_a+2*1*I_ESP_Fx2z_Fx2y;
    abcd[iGrid*1260+706] = 4.0E0*I_ESP_Kx3y3z_Fx2y_aa-2.0E0*1*I_ESP_Hxy3z_Fx2y_a-2.0E0*2*I_ESP_Hxy3z_Fx2y_a;
    abcd[iGrid*1260+707] = 4.0E0*I_ESP_Kx2y4z_Fx2y_aa-2.0E0*1*I_ESP_Hx4z_Fx2y_a;
    abcd[iGrid*1260+708] = 4.0E0*I_ESP_K7y_Fx2y_aa-2.0E0*5*I_ESP_H5y_Fx2y_a-2.0E0*6*I_ESP_H5y_Fx2y_a+5*4*I_ESP_F3y_Fx2y;
    abcd[iGrid*1260+709] = 4.0E0*I_ESP_K6yz_Fx2y_aa-2.0E0*4*I_ESP_H4yz_Fx2y_a-2.0E0*5*I_ESP_H4yz_Fx2y_a+4*3*I_ESP_F2yz_Fx2y;
    abcd[iGrid*1260+710] = 4.0E0*I_ESP_K5y2z_Fx2y_aa-2.0E0*3*I_ESP_H3y2z_Fx2y_a-2.0E0*4*I_ESP_H3y2z_Fx2y_a+3*2*I_ESP_Fy2z_Fx2y;
    abcd[iGrid*1260+711] = 4.0E0*I_ESP_K4y3z_Fx2y_aa-2.0E0*2*I_ESP_H2y3z_Fx2y_a-2.0E0*3*I_ESP_H2y3z_Fx2y_a+2*1*I_ESP_F3z_Fx2y;
    abcd[iGrid*1260+712] = 4.0E0*I_ESP_K3y4z_Fx2y_aa-2.0E0*1*I_ESP_Hy4z_Fx2y_a-2.0E0*2*I_ESP_Hy4z_Fx2y_a;
    abcd[iGrid*1260+713] = 4.0E0*I_ESP_K2y5z_Fx2y_aa-2.0E0*1*I_ESP_H5z_Fx2y_a;
    abcd[iGrid*1260+714] = 4.0E0*I_ESP_K5x2y_Fxyz_aa-2.0E0*1*I_ESP_H5x_Fxyz_a;
    abcd[iGrid*1260+715] = 4.0E0*I_ESP_K4x3y_Fxyz_aa-2.0E0*1*I_ESP_H4xy_Fxyz_a-2.0E0*2*I_ESP_H4xy_Fxyz_a;
    abcd[iGrid*1260+716] = 4.0E0*I_ESP_K4x2yz_Fxyz_aa-2.0E0*1*I_ESP_H4xz_Fxyz_a;
    abcd[iGrid*1260+717] = 4.0E0*I_ESP_K3x4y_Fxyz_aa-2.0E0*2*I_ESP_H3x2y_Fxyz_a-2.0E0*3*I_ESP_H3x2y_Fxyz_a+2*1*I_ESP_F3x_Fxyz;
    abcd[iGrid*1260+718] = 4.0E0*I_ESP_K3x3yz_Fxyz_aa-2.0E0*1*I_ESP_H3xyz_Fxyz_a-2.0E0*2*I_ESP_H3xyz_Fxyz_a;
    abcd[iGrid*1260+719] = 4.0E0*I_ESP_K3x2y2z_Fxyz_aa-2.0E0*1*I_ESP_H3x2z_Fxyz_a;
    abcd[iGrid*1260+720] = 4.0E0*I_ESP_K2x5y_Fxyz_aa-2.0E0*3*I_ESP_H2x3y_Fxyz_a-2.0E0*4*I_ESP_H2x3y_Fxyz_a+3*2*I_ESP_F2xy_Fxyz;
    abcd[iGrid*1260+721] = 4.0E0*I_ESP_K2x4yz_Fxyz_aa-2.0E0*2*I_ESP_H2x2yz_Fxyz_a-2.0E0*3*I_ESP_H2x2yz_Fxyz_a+2*1*I_ESP_F2xz_Fxyz;
    abcd[iGrid*1260+722] = 4.0E0*I_ESP_K2x3y2z_Fxyz_aa-2.0E0*1*I_ESP_H2xy2z_Fxyz_a-2.0E0*2*I_ESP_H2xy2z_Fxyz_a;
    abcd[iGrid*1260+723] = 4.0E0*I_ESP_K2x2y3z_Fxyz_aa-2.0E0*1*I_ESP_H2x3z_Fxyz_a;
    abcd[iGrid*1260+724] = 4.0E0*I_ESP_Kx6y_Fxyz_aa-2.0E0*4*I_ESP_Hx4y_Fxyz_a-2.0E0*5*I_ESP_Hx4y_Fxyz_a+4*3*I_ESP_Fx2y_Fxyz;
    abcd[iGrid*1260+725] = 4.0E0*I_ESP_Kx5yz_Fxyz_aa-2.0E0*3*I_ESP_Hx3yz_Fxyz_a-2.0E0*4*I_ESP_Hx3yz_Fxyz_a+3*2*I_ESP_Fxyz_Fxyz;
    abcd[iGrid*1260+726] = 4.0E0*I_ESP_Kx4y2z_Fxyz_aa-2.0E0*2*I_ESP_Hx2y2z_Fxyz_a-2.0E0*3*I_ESP_Hx2y2z_Fxyz_a+2*1*I_ESP_Fx2z_Fxyz;
    abcd[iGrid*1260+727] = 4.0E0*I_ESP_Kx3y3z_Fxyz_aa-2.0E0*1*I_ESP_Hxy3z_Fxyz_a-2.0E0*2*I_ESP_Hxy3z_Fxyz_a;
    abcd[iGrid*1260+728] = 4.0E0*I_ESP_Kx2y4z_Fxyz_aa-2.0E0*1*I_ESP_Hx4z_Fxyz_a;
    abcd[iGrid*1260+729] = 4.0E0*I_ESP_K7y_Fxyz_aa-2.0E0*5*I_ESP_H5y_Fxyz_a-2.0E0*6*I_ESP_H5y_Fxyz_a+5*4*I_ESP_F3y_Fxyz;
    abcd[iGrid*1260+730] = 4.0E0*I_ESP_K6yz_Fxyz_aa-2.0E0*4*I_ESP_H4yz_Fxyz_a-2.0E0*5*I_ESP_H4yz_Fxyz_a+4*3*I_ESP_F2yz_Fxyz;
    abcd[iGrid*1260+731] = 4.0E0*I_ESP_K5y2z_Fxyz_aa-2.0E0*3*I_ESP_H3y2z_Fxyz_a-2.0E0*4*I_ESP_H3y2z_Fxyz_a+3*2*I_ESP_Fy2z_Fxyz;
    abcd[iGrid*1260+732] = 4.0E0*I_ESP_K4y3z_Fxyz_aa-2.0E0*2*I_ESP_H2y3z_Fxyz_a-2.0E0*3*I_ESP_H2y3z_Fxyz_a+2*1*I_ESP_F3z_Fxyz;
    abcd[iGrid*1260+733] = 4.0E0*I_ESP_K3y4z_Fxyz_aa-2.0E0*1*I_ESP_Hy4z_Fxyz_a-2.0E0*2*I_ESP_Hy4z_Fxyz_a;
    abcd[iGrid*1260+734] = 4.0E0*I_ESP_K2y5z_Fxyz_aa-2.0E0*1*I_ESP_H5z_Fxyz_a;
    abcd[iGrid*1260+735] = 4.0E0*I_ESP_K5x2y_Fx2z_aa-2.0E0*1*I_ESP_H5x_Fx2z_a;
    abcd[iGrid*1260+736] = 4.0E0*I_ESP_K4x3y_Fx2z_aa-2.0E0*1*I_ESP_H4xy_Fx2z_a-2.0E0*2*I_ESP_H4xy_Fx2z_a;
    abcd[iGrid*1260+737] = 4.0E0*I_ESP_K4x2yz_Fx2z_aa-2.0E0*1*I_ESP_H4xz_Fx2z_a;
    abcd[iGrid*1260+738] = 4.0E0*I_ESP_K3x4y_Fx2z_aa-2.0E0*2*I_ESP_H3x2y_Fx2z_a-2.0E0*3*I_ESP_H3x2y_Fx2z_a+2*1*I_ESP_F3x_Fx2z;
    abcd[iGrid*1260+739] = 4.0E0*I_ESP_K3x3yz_Fx2z_aa-2.0E0*1*I_ESP_H3xyz_Fx2z_a-2.0E0*2*I_ESP_H3xyz_Fx2z_a;
    abcd[iGrid*1260+740] = 4.0E0*I_ESP_K3x2y2z_Fx2z_aa-2.0E0*1*I_ESP_H3x2z_Fx2z_a;
    abcd[iGrid*1260+741] = 4.0E0*I_ESP_K2x5y_Fx2z_aa-2.0E0*3*I_ESP_H2x3y_Fx2z_a-2.0E0*4*I_ESP_H2x3y_Fx2z_a+3*2*I_ESP_F2xy_Fx2z;
    abcd[iGrid*1260+742] = 4.0E0*I_ESP_K2x4yz_Fx2z_aa-2.0E0*2*I_ESP_H2x2yz_Fx2z_a-2.0E0*3*I_ESP_H2x2yz_Fx2z_a+2*1*I_ESP_F2xz_Fx2z;
    abcd[iGrid*1260+743] = 4.0E0*I_ESP_K2x3y2z_Fx2z_aa-2.0E0*1*I_ESP_H2xy2z_Fx2z_a-2.0E0*2*I_ESP_H2xy2z_Fx2z_a;
    abcd[iGrid*1260+744] = 4.0E0*I_ESP_K2x2y3z_Fx2z_aa-2.0E0*1*I_ESP_H2x3z_Fx2z_a;
    abcd[iGrid*1260+745] = 4.0E0*I_ESP_Kx6y_Fx2z_aa-2.0E0*4*I_ESP_Hx4y_Fx2z_a-2.0E0*5*I_ESP_Hx4y_Fx2z_a+4*3*I_ESP_Fx2y_Fx2z;
    abcd[iGrid*1260+746] = 4.0E0*I_ESP_Kx5yz_Fx2z_aa-2.0E0*3*I_ESP_Hx3yz_Fx2z_a-2.0E0*4*I_ESP_Hx3yz_Fx2z_a+3*2*I_ESP_Fxyz_Fx2z;
    abcd[iGrid*1260+747] = 4.0E0*I_ESP_Kx4y2z_Fx2z_aa-2.0E0*2*I_ESP_Hx2y2z_Fx2z_a-2.0E0*3*I_ESP_Hx2y2z_Fx2z_a+2*1*I_ESP_Fx2z_Fx2z;
    abcd[iGrid*1260+748] = 4.0E0*I_ESP_Kx3y3z_Fx2z_aa-2.0E0*1*I_ESP_Hxy3z_Fx2z_a-2.0E0*2*I_ESP_Hxy3z_Fx2z_a;
    abcd[iGrid*1260+749] = 4.0E0*I_ESP_Kx2y4z_Fx2z_aa-2.0E0*1*I_ESP_Hx4z_Fx2z_a;
    abcd[iGrid*1260+750] = 4.0E0*I_ESP_K7y_Fx2z_aa-2.0E0*5*I_ESP_H5y_Fx2z_a-2.0E0*6*I_ESP_H5y_Fx2z_a+5*4*I_ESP_F3y_Fx2z;
    abcd[iGrid*1260+751] = 4.0E0*I_ESP_K6yz_Fx2z_aa-2.0E0*4*I_ESP_H4yz_Fx2z_a-2.0E0*5*I_ESP_H4yz_Fx2z_a+4*3*I_ESP_F2yz_Fx2z;
    abcd[iGrid*1260+752] = 4.0E0*I_ESP_K5y2z_Fx2z_aa-2.0E0*3*I_ESP_H3y2z_Fx2z_a-2.0E0*4*I_ESP_H3y2z_Fx2z_a+3*2*I_ESP_Fy2z_Fx2z;
    abcd[iGrid*1260+753] = 4.0E0*I_ESP_K4y3z_Fx2z_aa-2.0E0*2*I_ESP_H2y3z_Fx2z_a-2.0E0*3*I_ESP_H2y3z_Fx2z_a+2*1*I_ESP_F3z_Fx2z;
    abcd[iGrid*1260+754] = 4.0E0*I_ESP_K3y4z_Fx2z_aa-2.0E0*1*I_ESP_Hy4z_Fx2z_a-2.0E0*2*I_ESP_Hy4z_Fx2z_a;
    abcd[iGrid*1260+755] = 4.0E0*I_ESP_K2y5z_Fx2z_aa-2.0E0*1*I_ESP_H5z_Fx2z_a;
    abcd[iGrid*1260+756] = 4.0E0*I_ESP_K5x2y_F3y_aa-2.0E0*1*I_ESP_H5x_F3y_a;
    abcd[iGrid*1260+757] = 4.0E0*I_ESP_K4x3y_F3y_aa-2.0E0*1*I_ESP_H4xy_F3y_a-2.0E0*2*I_ESP_H4xy_F3y_a;
    abcd[iGrid*1260+758] = 4.0E0*I_ESP_K4x2yz_F3y_aa-2.0E0*1*I_ESP_H4xz_F3y_a;
    abcd[iGrid*1260+759] = 4.0E0*I_ESP_K3x4y_F3y_aa-2.0E0*2*I_ESP_H3x2y_F3y_a-2.0E0*3*I_ESP_H3x2y_F3y_a+2*1*I_ESP_F3x_F3y;
    abcd[iGrid*1260+760] = 4.0E0*I_ESP_K3x3yz_F3y_aa-2.0E0*1*I_ESP_H3xyz_F3y_a-2.0E0*2*I_ESP_H3xyz_F3y_a;
    abcd[iGrid*1260+761] = 4.0E0*I_ESP_K3x2y2z_F3y_aa-2.0E0*1*I_ESP_H3x2z_F3y_a;
    abcd[iGrid*1260+762] = 4.0E0*I_ESP_K2x5y_F3y_aa-2.0E0*3*I_ESP_H2x3y_F3y_a-2.0E0*4*I_ESP_H2x3y_F3y_a+3*2*I_ESP_F2xy_F3y;
    abcd[iGrid*1260+763] = 4.0E0*I_ESP_K2x4yz_F3y_aa-2.0E0*2*I_ESP_H2x2yz_F3y_a-2.0E0*3*I_ESP_H2x2yz_F3y_a+2*1*I_ESP_F2xz_F3y;
    abcd[iGrid*1260+764] = 4.0E0*I_ESP_K2x3y2z_F3y_aa-2.0E0*1*I_ESP_H2xy2z_F3y_a-2.0E0*2*I_ESP_H2xy2z_F3y_a;
    abcd[iGrid*1260+765] = 4.0E0*I_ESP_K2x2y3z_F3y_aa-2.0E0*1*I_ESP_H2x3z_F3y_a;
    abcd[iGrid*1260+766] = 4.0E0*I_ESP_Kx6y_F3y_aa-2.0E0*4*I_ESP_Hx4y_F3y_a-2.0E0*5*I_ESP_Hx4y_F3y_a+4*3*I_ESP_Fx2y_F3y;
    abcd[iGrid*1260+767] = 4.0E0*I_ESP_Kx5yz_F3y_aa-2.0E0*3*I_ESP_Hx3yz_F3y_a-2.0E0*4*I_ESP_Hx3yz_F3y_a+3*2*I_ESP_Fxyz_F3y;
    abcd[iGrid*1260+768] = 4.0E0*I_ESP_Kx4y2z_F3y_aa-2.0E0*2*I_ESP_Hx2y2z_F3y_a-2.0E0*3*I_ESP_Hx2y2z_F3y_a+2*1*I_ESP_Fx2z_F3y;
    abcd[iGrid*1260+769] = 4.0E0*I_ESP_Kx3y3z_F3y_aa-2.0E0*1*I_ESP_Hxy3z_F3y_a-2.0E0*2*I_ESP_Hxy3z_F3y_a;
    abcd[iGrid*1260+770] = 4.0E0*I_ESP_Kx2y4z_F3y_aa-2.0E0*1*I_ESP_Hx4z_F3y_a;
    abcd[iGrid*1260+771] = 4.0E0*I_ESP_K7y_F3y_aa-2.0E0*5*I_ESP_H5y_F3y_a-2.0E0*6*I_ESP_H5y_F3y_a+5*4*I_ESP_F3y_F3y;
    abcd[iGrid*1260+772] = 4.0E0*I_ESP_K6yz_F3y_aa-2.0E0*4*I_ESP_H4yz_F3y_a-2.0E0*5*I_ESP_H4yz_F3y_a+4*3*I_ESP_F2yz_F3y;
    abcd[iGrid*1260+773] = 4.0E0*I_ESP_K5y2z_F3y_aa-2.0E0*3*I_ESP_H3y2z_F3y_a-2.0E0*4*I_ESP_H3y2z_F3y_a+3*2*I_ESP_Fy2z_F3y;
    abcd[iGrid*1260+774] = 4.0E0*I_ESP_K4y3z_F3y_aa-2.0E0*2*I_ESP_H2y3z_F3y_a-2.0E0*3*I_ESP_H2y3z_F3y_a+2*1*I_ESP_F3z_F3y;
    abcd[iGrid*1260+775] = 4.0E0*I_ESP_K3y4z_F3y_aa-2.0E0*1*I_ESP_Hy4z_F3y_a-2.0E0*2*I_ESP_Hy4z_F3y_a;
    abcd[iGrid*1260+776] = 4.0E0*I_ESP_K2y5z_F3y_aa-2.0E0*1*I_ESP_H5z_F3y_a;
    abcd[iGrid*1260+777] = 4.0E0*I_ESP_K5x2y_F2yz_aa-2.0E0*1*I_ESP_H5x_F2yz_a;
    abcd[iGrid*1260+778] = 4.0E0*I_ESP_K4x3y_F2yz_aa-2.0E0*1*I_ESP_H4xy_F2yz_a-2.0E0*2*I_ESP_H4xy_F2yz_a;
    abcd[iGrid*1260+779] = 4.0E0*I_ESP_K4x2yz_F2yz_aa-2.0E0*1*I_ESP_H4xz_F2yz_a;
    abcd[iGrid*1260+780] = 4.0E0*I_ESP_K3x4y_F2yz_aa-2.0E0*2*I_ESP_H3x2y_F2yz_a-2.0E0*3*I_ESP_H3x2y_F2yz_a+2*1*I_ESP_F3x_F2yz;
    abcd[iGrid*1260+781] = 4.0E0*I_ESP_K3x3yz_F2yz_aa-2.0E0*1*I_ESP_H3xyz_F2yz_a-2.0E0*2*I_ESP_H3xyz_F2yz_a;
    abcd[iGrid*1260+782] = 4.0E0*I_ESP_K3x2y2z_F2yz_aa-2.0E0*1*I_ESP_H3x2z_F2yz_a;
    abcd[iGrid*1260+783] = 4.0E0*I_ESP_K2x5y_F2yz_aa-2.0E0*3*I_ESP_H2x3y_F2yz_a-2.0E0*4*I_ESP_H2x3y_F2yz_a+3*2*I_ESP_F2xy_F2yz;
    abcd[iGrid*1260+784] = 4.0E0*I_ESP_K2x4yz_F2yz_aa-2.0E0*2*I_ESP_H2x2yz_F2yz_a-2.0E0*3*I_ESP_H2x2yz_F2yz_a+2*1*I_ESP_F2xz_F2yz;
    abcd[iGrid*1260+785] = 4.0E0*I_ESP_K2x3y2z_F2yz_aa-2.0E0*1*I_ESP_H2xy2z_F2yz_a-2.0E0*2*I_ESP_H2xy2z_F2yz_a;
    abcd[iGrid*1260+786] = 4.0E0*I_ESP_K2x2y3z_F2yz_aa-2.0E0*1*I_ESP_H2x3z_F2yz_a;
    abcd[iGrid*1260+787] = 4.0E0*I_ESP_Kx6y_F2yz_aa-2.0E0*4*I_ESP_Hx4y_F2yz_a-2.0E0*5*I_ESP_Hx4y_F2yz_a+4*3*I_ESP_Fx2y_F2yz;
    abcd[iGrid*1260+788] = 4.0E0*I_ESP_Kx5yz_F2yz_aa-2.0E0*3*I_ESP_Hx3yz_F2yz_a-2.0E0*4*I_ESP_Hx3yz_F2yz_a+3*2*I_ESP_Fxyz_F2yz;
    abcd[iGrid*1260+789] = 4.0E0*I_ESP_Kx4y2z_F2yz_aa-2.0E0*2*I_ESP_Hx2y2z_F2yz_a-2.0E0*3*I_ESP_Hx2y2z_F2yz_a+2*1*I_ESP_Fx2z_F2yz;
    abcd[iGrid*1260+790] = 4.0E0*I_ESP_Kx3y3z_F2yz_aa-2.0E0*1*I_ESP_Hxy3z_F2yz_a-2.0E0*2*I_ESP_Hxy3z_F2yz_a;
    abcd[iGrid*1260+791] = 4.0E0*I_ESP_Kx2y4z_F2yz_aa-2.0E0*1*I_ESP_Hx4z_F2yz_a;
    abcd[iGrid*1260+792] = 4.0E0*I_ESP_K7y_F2yz_aa-2.0E0*5*I_ESP_H5y_F2yz_a-2.0E0*6*I_ESP_H5y_F2yz_a+5*4*I_ESP_F3y_F2yz;
    abcd[iGrid*1260+793] = 4.0E0*I_ESP_K6yz_F2yz_aa-2.0E0*4*I_ESP_H4yz_F2yz_a-2.0E0*5*I_ESP_H4yz_F2yz_a+4*3*I_ESP_F2yz_F2yz;
    abcd[iGrid*1260+794] = 4.0E0*I_ESP_K5y2z_F2yz_aa-2.0E0*3*I_ESP_H3y2z_F2yz_a-2.0E0*4*I_ESP_H3y2z_F2yz_a+3*2*I_ESP_Fy2z_F2yz;
    abcd[iGrid*1260+795] = 4.0E0*I_ESP_K4y3z_F2yz_aa-2.0E0*2*I_ESP_H2y3z_F2yz_a-2.0E0*3*I_ESP_H2y3z_F2yz_a+2*1*I_ESP_F3z_F2yz;
    abcd[iGrid*1260+796] = 4.0E0*I_ESP_K3y4z_F2yz_aa-2.0E0*1*I_ESP_Hy4z_F2yz_a-2.0E0*2*I_ESP_Hy4z_F2yz_a;
    abcd[iGrid*1260+797] = 4.0E0*I_ESP_K2y5z_F2yz_aa-2.0E0*1*I_ESP_H5z_F2yz_a;
    abcd[iGrid*1260+798] = 4.0E0*I_ESP_K5x2y_Fy2z_aa-2.0E0*1*I_ESP_H5x_Fy2z_a;
    abcd[iGrid*1260+799] = 4.0E0*I_ESP_K4x3y_Fy2z_aa-2.0E0*1*I_ESP_H4xy_Fy2z_a-2.0E0*2*I_ESP_H4xy_Fy2z_a;
    abcd[iGrid*1260+800] = 4.0E0*I_ESP_K4x2yz_Fy2z_aa-2.0E0*1*I_ESP_H4xz_Fy2z_a;
    abcd[iGrid*1260+801] = 4.0E0*I_ESP_K3x4y_Fy2z_aa-2.0E0*2*I_ESP_H3x2y_Fy2z_a-2.0E0*3*I_ESP_H3x2y_Fy2z_a+2*1*I_ESP_F3x_Fy2z;
    abcd[iGrid*1260+802] = 4.0E0*I_ESP_K3x3yz_Fy2z_aa-2.0E0*1*I_ESP_H3xyz_Fy2z_a-2.0E0*2*I_ESP_H3xyz_Fy2z_a;
    abcd[iGrid*1260+803] = 4.0E0*I_ESP_K3x2y2z_Fy2z_aa-2.0E0*1*I_ESP_H3x2z_Fy2z_a;
    abcd[iGrid*1260+804] = 4.0E0*I_ESP_K2x5y_Fy2z_aa-2.0E0*3*I_ESP_H2x3y_Fy2z_a-2.0E0*4*I_ESP_H2x3y_Fy2z_a+3*2*I_ESP_F2xy_Fy2z;
    abcd[iGrid*1260+805] = 4.0E0*I_ESP_K2x4yz_Fy2z_aa-2.0E0*2*I_ESP_H2x2yz_Fy2z_a-2.0E0*3*I_ESP_H2x2yz_Fy2z_a+2*1*I_ESP_F2xz_Fy2z;
    abcd[iGrid*1260+806] = 4.0E0*I_ESP_K2x3y2z_Fy2z_aa-2.0E0*1*I_ESP_H2xy2z_Fy2z_a-2.0E0*2*I_ESP_H2xy2z_Fy2z_a;
    abcd[iGrid*1260+807] = 4.0E0*I_ESP_K2x2y3z_Fy2z_aa-2.0E0*1*I_ESP_H2x3z_Fy2z_a;
    abcd[iGrid*1260+808] = 4.0E0*I_ESP_Kx6y_Fy2z_aa-2.0E0*4*I_ESP_Hx4y_Fy2z_a-2.0E0*5*I_ESP_Hx4y_Fy2z_a+4*3*I_ESP_Fx2y_Fy2z;
    abcd[iGrid*1260+809] = 4.0E0*I_ESP_Kx5yz_Fy2z_aa-2.0E0*3*I_ESP_Hx3yz_Fy2z_a-2.0E0*4*I_ESP_Hx3yz_Fy2z_a+3*2*I_ESP_Fxyz_Fy2z;
    abcd[iGrid*1260+810] = 4.0E0*I_ESP_Kx4y2z_Fy2z_aa-2.0E0*2*I_ESP_Hx2y2z_Fy2z_a-2.0E0*3*I_ESP_Hx2y2z_Fy2z_a+2*1*I_ESP_Fx2z_Fy2z;
    abcd[iGrid*1260+811] = 4.0E0*I_ESP_Kx3y3z_Fy2z_aa-2.0E0*1*I_ESP_Hxy3z_Fy2z_a-2.0E0*2*I_ESP_Hxy3z_Fy2z_a;
    abcd[iGrid*1260+812] = 4.0E0*I_ESP_Kx2y4z_Fy2z_aa-2.0E0*1*I_ESP_Hx4z_Fy2z_a;
    abcd[iGrid*1260+813] = 4.0E0*I_ESP_K7y_Fy2z_aa-2.0E0*5*I_ESP_H5y_Fy2z_a-2.0E0*6*I_ESP_H5y_Fy2z_a+5*4*I_ESP_F3y_Fy2z;
    abcd[iGrid*1260+814] = 4.0E0*I_ESP_K6yz_Fy2z_aa-2.0E0*4*I_ESP_H4yz_Fy2z_a-2.0E0*5*I_ESP_H4yz_Fy2z_a+4*3*I_ESP_F2yz_Fy2z;
    abcd[iGrid*1260+815] = 4.0E0*I_ESP_K5y2z_Fy2z_aa-2.0E0*3*I_ESP_H3y2z_Fy2z_a-2.0E0*4*I_ESP_H3y2z_Fy2z_a+3*2*I_ESP_Fy2z_Fy2z;
    abcd[iGrid*1260+816] = 4.0E0*I_ESP_K4y3z_Fy2z_aa-2.0E0*2*I_ESP_H2y3z_Fy2z_a-2.0E0*3*I_ESP_H2y3z_Fy2z_a+2*1*I_ESP_F3z_Fy2z;
    abcd[iGrid*1260+817] = 4.0E0*I_ESP_K3y4z_Fy2z_aa-2.0E0*1*I_ESP_Hy4z_Fy2z_a-2.0E0*2*I_ESP_Hy4z_Fy2z_a;
    abcd[iGrid*1260+818] = 4.0E0*I_ESP_K2y5z_Fy2z_aa-2.0E0*1*I_ESP_H5z_Fy2z_a;
    abcd[iGrid*1260+819] = 4.0E0*I_ESP_K5x2y_F3z_aa-2.0E0*1*I_ESP_H5x_F3z_a;
    abcd[iGrid*1260+820] = 4.0E0*I_ESP_K4x3y_F3z_aa-2.0E0*1*I_ESP_H4xy_F3z_a-2.0E0*2*I_ESP_H4xy_F3z_a;
    abcd[iGrid*1260+821] = 4.0E0*I_ESP_K4x2yz_F3z_aa-2.0E0*1*I_ESP_H4xz_F3z_a;
    abcd[iGrid*1260+822] = 4.0E0*I_ESP_K3x4y_F3z_aa-2.0E0*2*I_ESP_H3x2y_F3z_a-2.0E0*3*I_ESP_H3x2y_F3z_a+2*1*I_ESP_F3x_F3z;
    abcd[iGrid*1260+823] = 4.0E0*I_ESP_K3x3yz_F3z_aa-2.0E0*1*I_ESP_H3xyz_F3z_a-2.0E0*2*I_ESP_H3xyz_F3z_a;
    abcd[iGrid*1260+824] = 4.0E0*I_ESP_K3x2y2z_F3z_aa-2.0E0*1*I_ESP_H3x2z_F3z_a;
    abcd[iGrid*1260+825] = 4.0E0*I_ESP_K2x5y_F3z_aa-2.0E0*3*I_ESP_H2x3y_F3z_a-2.0E0*4*I_ESP_H2x3y_F3z_a+3*2*I_ESP_F2xy_F3z;
    abcd[iGrid*1260+826] = 4.0E0*I_ESP_K2x4yz_F3z_aa-2.0E0*2*I_ESP_H2x2yz_F3z_a-2.0E0*3*I_ESP_H2x2yz_F3z_a+2*1*I_ESP_F2xz_F3z;
    abcd[iGrid*1260+827] = 4.0E0*I_ESP_K2x3y2z_F3z_aa-2.0E0*1*I_ESP_H2xy2z_F3z_a-2.0E0*2*I_ESP_H2xy2z_F3z_a;
    abcd[iGrid*1260+828] = 4.0E0*I_ESP_K2x2y3z_F3z_aa-2.0E0*1*I_ESP_H2x3z_F3z_a;
    abcd[iGrid*1260+829] = 4.0E0*I_ESP_Kx6y_F3z_aa-2.0E0*4*I_ESP_Hx4y_F3z_a-2.0E0*5*I_ESP_Hx4y_F3z_a+4*3*I_ESP_Fx2y_F3z;
    abcd[iGrid*1260+830] = 4.0E0*I_ESP_Kx5yz_F3z_aa-2.0E0*3*I_ESP_Hx3yz_F3z_a-2.0E0*4*I_ESP_Hx3yz_F3z_a+3*2*I_ESP_Fxyz_F3z;
    abcd[iGrid*1260+831] = 4.0E0*I_ESP_Kx4y2z_F3z_aa-2.0E0*2*I_ESP_Hx2y2z_F3z_a-2.0E0*3*I_ESP_Hx2y2z_F3z_a+2*1*I_ESP_Fx2z_F3z;
    abcd[iGrid*1260+832] = 4.0E0*I_ESP_Kx3y3z_F3z_aa-2.0E0*1*I_ESP_Hxy3z_F3z_a-2.0E0*2*I_ESP_Hxy3z_F3z_a;
    abcd[iGrid*1260+833] = 4.0E0*I_ESP_Kx2y4z_F3z_aa-2.0E0*1*I_ESP_Hx4z_F3z_a;
    abcd[iGrid*1260+834] = 4.0E0*I_ESP_K7y_F3z_aa-2.0E0*5*I_ESP_H5y_F3z_a-2.0E0*6*I_ESP_H5y_F3z_a+5*4*I_ESP_F3y_F3z;
    abcd[iGrid*1260+835] = 4.0E0*I_ESP_K6yz_F3z_aa-2.0E0*4*I_ESP_H4yz_F3z_a-2.0E0*5*I_ESP_H4yz_F3z_a+4*3*I_ESP_F2yz_F3z;
    abcd[iGrid*1260+836] = 4.0E0*I_ESP_K5y2z_F3z_aa-2.0E0*3*I_ESP_H3y2z_F3z_a-2.0E0*4*I_ESP_H3y2z_F3z_a+3*2*I_ESP_Fy2z_F3z;
    abcd[iGrid*1260+837] = 4.0E0*I_ESP_K4y3z_F3z_aa-2.0E0*2*I_ESP_H2y3z_F3z_a-2.0E0*3*I_ESP_H2y3z_F3z_a+2*1*I_ESP_F3z_F3z;
    abcd[iGrid*1260+838] = 4.0E0*I_ESP_K3y4z_F3z_aa-2.0E0*1*I_ESP_Hy4z_F3z_a-2.0E0*2*I_ESP_Hy4z_F3z_a;
    abcd[iGrid*1260+839] = 4.0E0*I_ESP_K2y5z_F3z_aa-2.0E0*1*I_ESP_H5z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_F_aa
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    abcd[iGrid*1260+840] = 4.0E0*I_ESP_K5xyz_F3x_aa;
    abcd[iGrid*1260+841] = 4.0E0*I_ESP_K4x2yz_F3x_aa-2.0E0*1*I_ESP_H4xz_F3x_a;
    abcd[iGrid*1260+842] = 4.0E0*I_ESP_K4xy2z_F3x_aa-2.0E0*1*I_ESP_H4xy_F3x_a;
    abcd[iGrid*1260+843] = 4.0E0*I_ESP_K3x3yz_F3x_aa-2.0E0*2*I_ESP_H3xyz_F3x_a;
    abcd[iGrid*1260+844] = 4.0E0*I_ESP_K3x2y2z_F3x_aa-2.0E0*1*I_ESP_H3x2y_F3x_a-2.0E0*1*I_ESP_H3x2z_F3x_a+1*I_ESP_F3x_F3x;
    abcd[iGrid*1260+845] = 4.0E0*I_ESP_K3xy3z_F3x_aa-2.0E0*2*I_ESP_H3xyz_F3x_a;
    abcd[iGrid*1260+846] = 4.0E0*I_ESP_K2x4yz_F3x_aa-2.0E0*3*I_ESP_H2x2yz_F3x_a;
    abcd[iGrid*1260+847] = 4.0E0*I_ESP_K2x3y2z_F3x_aa-2.0E0*1*I_ESP_H2x3y_F3x_a-2.0E0*2*I_ESP_H2xy2z_F3x_a+2*1*I_ESP_F2xy_F3x;
    abcd[iGrid*1260+848] = 4.0E0*I_ESP_K2x2y3z_F3x_aa-2.0E0*2*I_ESP_H2x2yz_F3x_a-2.0E0*1*I_ESP_H2x3z_F3x_a+2*I_ESP_F2xz_F3x;
    abcd[iGrid*1260+849] = 4.0E0*I_ESP_K2xy4z_F3x_aa-2.0E0*3*I_ESP_H2xy2z_F3x_a;
    abcd[iGrid*1260+850] = 4.0E0*I_ESP_Kx5yz_F3x_aa-2.0E0*4*I_ESP_Hx3yz_F3x_a;
    abcd[iGrid*1260+851] = 4.0E0*I_ESP_Kx4y2z_F3x_aa-2.0E0*1*I_ESP_Hx4y_F3x_a-2.0E0*3*I_ESP_Hx2y2z_F3x_a+3*1*I_ESP_Fx2y_F3x;
    abcd[iGrid*1260+852] = 4.0E0*I_ESP_Kx3y3z_F3x_aa-2.0E0*2*I_ESP_Hx3yz_F3x_a-2.0E0*2*I_ESP_Hxy3z_F3x_a+2*2*I_ESP_Fxyz_F3x;
    abcd[iGrid*1260+853] = 4.0E0*I_ESP_Kx2y4z_F3x_aa-2.0E0*3*I_ESP_Hx2y2z_F3x_a-2.0E0*1*I_ESP_Hx4z_F3x_a+3*I_ESP_Fx2z_F3x;
    abcd[iGrid*1260+854] = 4.0E0*I_ESP_Kxy5z_F3x_aa-2.0E0*4*I_ESP_Hxy3z_F3x_a;
    abcd[iGrid*1260+855] = 4.0E0*I_ESP_K6yz_F3x_aa-2.0E0*5*I_ESP_H4yz_F3x_a;
    abcd[iGrid*1260+856] = 4.0E0*I_ESP_K5y2z_F3x_aa-2.0E0*1*I_ESP_H5y_F3x_a-2.0E0*4*I_ESP_H3y2z_F3x_a+4*1*I_ESP_F3y_F3x;
    abcd[iGrid*1260+857] = 4.0E0*I_ESP_K4y3z_F3x_aa-2.0E0*2*I_ESP_H4yz_F3x_a-2.0E0*3*I_ESP_H2y3z_F3x_a+3*2*I_ESP_F2yz_F3x;
    abcd[iGrid*1260+858] = 4.0E0*I_ESP_K3y4z_F3x_aa-2.0E0*3*I_ESP_H3y2z_F3x_a-2.0E0*2*I_ESP_Hy4z_F3x_a+2*3*I_ESP_Fy2z_F3x;
    abcd[iGrid*1260+859] = 4.0E0*I_ESP_K2y5z_F3x_aa-2.0E0*4*I_ESP_H2y3z_F3x_a-2.0E0*1*I_ESP_H5z_F3x_a+4*I_ESP_F3z_F3x;
    abcd[iGrid*1260+860] = 4.0E0*I_ESP_Ky6z_F3x_aa-2.0E0*5*I_ESP_Hy4z_F3x_a;
    abcd[iGrid*1260+861] = 4.0E0*I_ESP_K5xyz_F2xy_aa;
    abcd[iGrid*1260+862] = 4.0E0*I_ESP_K4x2yz_F2xy_aa-2.0E0*1*I_ESP_H4xz_F2xy_a;
    abcd[iGrid*1260+863] = 4.0E0*I_ESP_K4xy2z_F2xy_aa-2.0E0*1*I_ESP_H4xy_F2xy_a;
    abcd[iGrid*1260+864] = 4.0E0*I_ESP_K3x3yz_F2xy_aa-2.0E0*2*I_ESP_H3xyz_F2xy_a;
    abcd[iGrid*1260+865] = 4.0E0*I_ESP_K3x2y2z_F2xy_aa-2.0E0*1*I_ESP_H3x2y_F2xy_a-2.0E0*1*I_ESP_H3x2z_F2xy_a+1*I_ESP_F3x_F2xy;
    abcd[iGrid*1260+866] = 4.0E0*I_ESP_K3xy3z_F2xy_aa-2.0E0*2*I_ESP_H3xyz_F2xy_a;
    abcd[iGrid*1260+867] = 4.0E0*I_ESP_K2x4yz_F2xy_aa-2.0E0*3*I_ESP_H2x2yz_F2xy_a;
    abcd[iGrid*1260+868] = 4.0E0*I_ESP_K2x3y2z_F2xy_aa-2.0E0*1*I_ESP_H2x3y_F2xy_a-2.0E0*2*I_ESP_H2xy2z_F2xy_a+2*1*I_ESP_F2xy_F2xy;
    abcd[iGrid*1260+869] = 4.0E0*I_ESP_K2x2y3z_F2xy_aa-2.0E0*2*I_ESP_H2x2yz_F2xy_a-2.0E0*1*I_ESP_H2x3z_F2xy_a+2*I_ESP_F2xz_F2xy;
    abcd[iGrid*1260+870] = 4.0E0*I_ESP_K2xy4z_F2xy_aa-2.0E0*3*I_ESP_H2xy2z_F2xy_a;
    abcd[iGrid*1260+871] = 4.0E0*I_ESP_Kx5yz_F2xy_aa-2.0E0*4*I_ESP_Hx3yz_F2xy_a;
    abcd[iGrid*1260+872] = 4.0E0*I_ESP_Kx4y2z_F2xy_aa-2.0E0*1*I_ESP_Hx4y_F2xy_a-2.0E0*3*I_ESP_Hx2y2z_F2xy_a+3*1*I_ESP_Fx2y_F2xy;
    abcd[iGrid*1260+873] = 4.0E0*I_ESP_Kx3y3z_F2xy_aa-2.0E0*2*I_ESP_Hx3yz_F2xy_a-2.0E0*2*I_ESP_Hxy3z_F2xy_a+2*2*I_ESP_Fxyz_F2xy;
    abcd[iGrid*1260+874] = 4.0E0*I_ESP_Kx2y4z_F2xy_aa-2.0E0*3*I_ESP_Hx2y2z_F2xy_a-2.0E0*1*I_ESP_Hx4z_F2xy_a+3*I_ESP_Fx2z_F2xy;
    abcd[iGrid*1260+875] = 4.0E0*I_ESP_Kxy5z_F2xy_aa-2.0E0*4*I_ESP_Hxy3z_F2xy_a;
    abcd[iGrid*1260+876] = 4.0E0*I_ESP_K6yz_F2xy_aa-2.0E0*5*I_ESP_H4yz_F2xy_a;
    abcd[iGrid*1260+877] = 4.0E0*I_ESP_K5y2z_F2xy_aa-2.0E0*1*I_ESP_H5y_F2xy_a-2.0E0*4*I_ESP_H3y2z_F2xy_a+4*1*I_ESP_F3y_F2xy;
    abcd[iGrid*1260+878] = 4.0E0*I_ESP_K4y3z_F2xy_aa-2.0E0*2*I_ESP_H4yz_F2xy_a-2.0E0*3*I_ESP_H2y3z_F2xy_a+3*2*I_ESP_F2yz_F2xy;
    abcd[iGrid*1260+879] = 4.0E0*I_ESP_K3y4z_F2xy_aa-2.0E0*3*I_ESP_H3y2z_F2xy_a-2.0E0*2*I_ESP_Hy4z_F2xy_a+2*3*I_ESP_Fy2z_F2xy;
    abcd[iGrid*1260+880] = 4.0E0*I_ESP_K2y5z_F2xy_aa-2.0E0*4*I_ESP_H2y3z_F2xy_a-2.0E0*1*I_ESP_H5z_F2xy_a+4*I_ESP_F3z_F2xy;
    abcd[iGrid*1260+881] = 4.0E0*I_ESP_Ky6z_F2xy_aa-2.0E0*5*I_ESP_Hy4z_F2xy_a;
    abcd[iGrid*1260+882] = 4.0E0*I_ESP_K5xyz_F2xz_aa;
    abcd[iGrid*1260+883] = 4.0E0*I_ESP_K4x2yz_F2xz_aa-2.0E0*1*I_ESP_H4xz_F2xz_a;
    abcd[iGrid*1260+884] = 4.0E0*I_ESP_K4xy2z_F2xz_aa-2.0E0*1*I_ESP_H4xy_F2xz_a;
    abcd[iGrid*1260+885] = 4.0E0*I_ESP_K3x3yz_F2xz_aa-2.0E0*2*I_ESP_H3xyz_F2xz_a;
    abcd[iGrid*1260+886] = 4.0E0*I_ESP_K3x2y2z_F2xz_aa-2.0E0*1*I_ESP_H3x2y_F2xz_a-2.0E0*1*I_ESP_H3x2z_F2xz_a+1*I_ESP_F3x_F2xz;
    abcd[iGrid*1260+887] = 4.0E0*I_ESP_K3xy3z_F2xz_aa-2.0E0*2*I_ESP_H3xyz_F2xz_a;
    abcd[iGrid*1260+888] = 4.0E0*I_ESP_K2x4yz_F2xz_aa-2.0E0*3*I_ESP_H2x2yz_F2xz_a;
    abcd[iGrid*1260+889] = 4.0E0*I_ESP_K2x3y2z_F2xz_aa-2.0E0*1*I_ESP_H2x3y_F2xz_a-2.0E0*2*I_ESP_H2xy2z_F2xz_a+2*1*I_ESP_F2xy_F2xz;
    abcd[iGrid*1260+890] = 4.0E0*I_ESP_K2x2y3z_F2xz_aa-2.0E0*2*I_ESP_H2x2yz_F2xz_a-2.0E0*1*I_ESP_H2x3z_F2xz_a+2*I_ESP_F2xz_F2xz;
    abcd[iGrid*1260+891] = 4.0E0*I_ESP_K2xy4z_F2xz_aa-2.0E0*3*I_ESP_H2xy2z_F2xz_a;
    abcd[iGrid*1260+892] = 4.0E0*I_ESP_Kx5yz_F2xz_aa-2.0E0*4*I_ESP_Hx3yz_F2xz_a;
    abcd[iGrid*1260+893] = 4.0E0*I_ESP_Kx4y2z_F2xz_aa-2.0E0*1*I_ESP_Hx4y_F2xz_a-2.0E0*3*I_ESP_Hx2y2z_F2xz_a+3*1*I_ESP_Fx2y_F2xz;
    abcd[iGrid*1260+894] = 4.0E0*I_ESP_Kx3y3z_F2xz_aa-2.0E0*2*I_ESP_Hx3yz_F2xz_a-2.0E0*2*I_ESP_Hxy3z_F2xz_a+2*2*I_ESP_Fxyz_F2xz;
    abcd[iGrid*1260+895] = 4.0E0*I_ESP_Kx2y4z_F2xz_aa-2.0E0*3*I_ESP_Hx2y2z_F2xz_a-2.0E0*1*I_ESP_Hx4z_F2xz_a+3*I_ESP_Fx2z_F2xz;
    abcd[iGrid*1260+896] = 4.0E0*I_ESP_Kxy5z_F2xz_aa-2.0E0*4*I_ESP_Hxy3z_F2xz_a;
    abcd[iGrid*1260+897] = 4.0E0*I_ESP_K6yz_F2xz_aa-2.0E0*5*I_ESP_H4yz_F2xz_a;
    abcd[iGrid*1260+898] = 4.0E0*I_ESP_K5y2z_F2xz_aa-2.0E0*1*I_ESP_H5y_F2xz_a-2.0E0*4*I_ESP_H3y2z_F2xz_a+4*1*I_ESP_F3y_F2xz;
    abcd[iGrid*1260+899] = 4.0E0*I_ESP_K4y3z_F2xz_aa-2.0E0*2*I_ESP_H4yz_F2xz_a-2.0E0*3*I_ESP_H2y3z_F2xz_a+3*2*I_ESP_F2yz_F2xz;
    abcd[iGrid*1260+900] = 4.0E0*I_ESP_K3y4z_F2xz_aa-2.0E0*3*I_ESP_H3y2z_F2xz_a-2.0E0*2*I_ESP_Hy4z_F2xz_a+2*3*I_ESP_Fy2z_F2xz;
    abcd[iGrid*1260+901] = 4.0E0*I_ESP_K2y5z_F2xz_aa-2.0E0*4*I_ESP_H2y3z_F2xz_a-2.0E0*1*I_ESP_H5z_F2xz_a+4*I_ESP_F3z_F2xz;
    abcd[iGrid*1260+902] = 4.0E0*I_ESP_Ky6z_F2xz_aa-2.0E0*5*I_ESP_Hy4z_F2xz_a;
    abcd[iGrid*1260+903] = 4.0E0*I_ESP_K5xyz_Fx2y_aa;
    abcd[iGrid*1260+904] = 4.0E0*I_ESP_K4x2yz_Fx2y_aa-2.0E0*1*I_ESP_H4xz_Fx2y_a;
    abcd[iGrid*1260+905] = 4.0E0*I_ESP_K4xy2z_Fx2y_aa-2.0E0*1*I_ESP_H4xy_Fx2y_a;
    abcd[iGrid*1260+906] = 4.0E0*I_ESP_K3x3yz_Fx2y_aa-2.0E0*2*I_ESP_H3xyz_Fx2y_a;
    abcd[iGrid*1260+907] = 4.0E0*I_ESP_K3x2y2z_Fx2y_aa-2.0E0*1*I_ESP_H3x2y_Fx2y_a-2.0E0*1*I_ESP_H3x2z_Fx2y_a+1*I_ESP_F3x_Fx2y;
    abcd[iGrid*1260+908] = 4.0E0*I_ESP_K3xy3z_Fx2y_aa-2.0E0*2*I_ESP_H3xyz_Fx2y_a;
    abcd[iGrid*1260+909] = 4.0E0*I_ESP_K2x4yz_Fx2y_aa-2.0E0*3*I_ESP_H2x2yz_Fx2y_a;
    abcd[iGrid*1260+910] = 4.0E0*I_ESP_K2x3y2z_Fx2y_aa-2.0E0*1*I_ESP_H2x3y_Fx2y_a-2.0E0*2*I_ESP_H2xy2z_Fx2y_a+2*1*I_ESP_F2xy_Fx2y;
    abcd[iGrid*1260+911] = 4.0E0*I_ESP_K2x2y3z_Fx2y_aa-2.0E0*2*I_ESP_H2x2yz_Fx2y_a-2.0E0*1*I_ESP_H2x3z_Fx2y_a+2*I_ESP_F2xz_Fx2y;
    abcd[iGrid*1260+912] = 4.0E0*I_ESP_K2xy4z_Fx2y_aa-2.0E0*3*I_ESP_H2xy2z_Fx2y_a;
    abcd[iGrid*1260+913] = 4.0E0*I_ESP_Kx5yz_Fx2y_aa-2.0E0*4*I_ESP_Hx3yz_Fx2y_a;
    abcd[iGrid*1260+914] = 4.0E0*I_ESP_Kx4y2z_Fx2y_aa-2.0E0*1*I_ESP_Hx4y_Fx2y_a-2.0E0*3*I_ESP_Hx2y2z_Fx2y_a+3*1*I_ESP_Fx2y_Fx2y;
    abcd[iGrid*1260+915] = 4.0E0*I_ESP_Kx3y3z_Fx2y_aa-2.0E0*2*I_ESP_Hx3yz_Fx2y_a-2.0E0*2*I_ESP_Hxy3z_Fx2y_a+2*2*I_ESP_Fxyz_Fx2y;
    abcd[iGrid*1260+916] = 4.0E0*I_ESP_Kx2y4z_Fx2y_aa-2.0E0*3*I_ESP_Hx2y2z_Fx2y_a-2.0E0*1*I_ESP_Hx4z_Fx2y_a+3*I_ESP_Fx2z_Fx2y;
    abcd[iGrid*1260+917] = 4.0E0*I_ESP_Kxy5z_Fx2y_aa-2.0E0*4*I_ESP_Hxy3z_Fx2y_a;
    abcd[iGrid*1260+918] = 4.0E0*I_ESP_K6yz_Fx2y_aa-2.0E0*5*I_ESP_H4yz_Fx2y_a;
    abcd[iGrid*1260+919] = 4.0E0*I_ESP_K5y2z_Fx2y_aa-2.0E0*1*I_ESP_H5y_Fx2y_a-2.0E0*4*I_ESP_H3y2z_Fx2y_a+4*1*I_ESP_F3y_Fx2y;
    abcd[iGrid*1260+920] = 4.0E0*I_ESP_K4y3z_Fx2y_aa-2.0E0*2*I_ESP_H4yz_Fx2y_a-2.0E0*3*I_ESP_H2y3z_Fx2y_a+3*2*I_ESP_F2yz_Fx2y;
    abcd[iGrid*1260+921] = 4.0E0*I_ESP_K3y4z_Fx2y_aa-2.0E0*3*I_ESP_H3y2z_Fx2y_a-2.0E0*2*I_ESP_Hy4z_Fx2y_a+2*3*I_ESP_Fy2z_Fx2y;
    abcd[iGrid*1260+922] = 4.0E0*I_ESP_K2y5z_Fx2y_aa-2.0E0*4*I_ESP_H2y3z_Fx2y_a-2.0E0*1*I_ESP_H5z_Fx2y_a+4*I_ESP_F3z_Fx2y;
    abcd[iGrid*1260+923] = 4.0E0*I_ESP_Ky6z_Fx2y_aa-2.0E0*5*I_ESP_Hy4z_Fx2y_a;
    abcd[iGrid*1260+924] = 4.0E0*I_ESP_K5xyz_Fxyz_aa;
    abcd[iGrid*1260+925] = 4.0E0*I_ESP_K4x2yz_Fxyz_aa-2.0E0*1*I_ESP_H4xz_Fxyz_a;
    abcd[iGrid*1260+926] = 4.0E0*I_ESP_K4xy2z_Fxyz_aa-2.0E0*1*I_ESP_H4xy_Fxyz_a;
    abcd[iGrid*1260+927] = 4.0E0*I_ESP_K3x3yz_Fxyz_aa-2.0E0*2*I_ESP_H3xyz_Fxyz_a;
    abcd[iGrid*1260+928] = 4.0E0*I_ESP_K3x2y2z_Fxyz_aa-2.0E0*1*I_ESP_H3x2y_Fxyz_a-2.0E0*1*I_ESP_H3x2z_Fxyz_a+1*I_ESP_F3x_Fxyz;
    abcd[iGrid*1260+929] = 4.0E0*I_ESP_K3xy3z_Fxyz_aa-2.0E0*2*I_ESP_H3xyz_Fxyz_a;
    abcd[iGrid*1260+930] = 4.0E0*I_ESP_K2x4yz_Fxyz_aa-2.0E0*3*I_ESP_H2x2yz_Fxyz_a;
    abcd[iGrid*1260+931] = 4.0E0*I_ESP_K2x3y2z_Fxyz_aa-2.0E0*1*I_ESP_H2x3y_Fxyz_a-2.0E0*2*I_ESP_H2xy2z_Fxyz_a+2*1*I_ESP_F2xy_Fxyz;
    abcd[iGrid*1260+932] = 4.0E0*I_ESP_K2x2y3z_Fxyz_aa-2.0E0*2*I_ESP_H2x2yz_Fxyz_a-2.0E0*1*I_ESP_H2x3z_Fxyz_a+2*I_ESP_F2xz_Fxyz;
    abcd[iGrid*1260+933] = 4.0E0*I_ESP_K2xy4z_Fxyz_aa-2.0E0*3*I_ESP_H2xy2z_Fxyz_a;
    abcd[iGrid*1260+934] = 4.0E0*I_ESP_Kx5yz_Fxyz_aa-2.0E0*4*I_ESP_Hx3yz_Fxyz_a;
    abcd[iGrid*1260+935] = 4.0E0*I_ESP_Kx4y2z_Fxyz_aa-2.0E0*1*I_ESP_Hx4y_Fxyz_a-2.0E0*3*I_ESP_Hx2y2z_Fxyz_a+3*1*I_ESP_Fx2y_Fxyz;
    abcd[iGrid*1260+936] = 4.0E0*I_ESP_Kx3y3z_Fxyz_aa-2.0E0*2*I_ESP_Hx3yz_Fxyz_a-2.0E0*2*I_ESP_Hxy3z_Fxyz_a+2*2*I_ESP_Fxyz_Fxyz;
    abcd[iGrid*1260+937] = 4.0E0*I_ESP_Kx2y4z_Fxyz_aa-2.0E0*3*I_ESP_Hx2y2z_Fxyz_a-2.0E0*1*I_ESP_Hx4z_Fxyz_a+3*I_ESP_Fx2z_Fxyz;
    abcd[iGrid*1260+938] = 4.0E0*I_ESP_Kxy5z_Fxyz_aa-2.0E0*4*I_ESP_Hxy3z_Fxyz_a;
    abcd[iGrid*1260+939] = 4.0E0*I_ESP_K6yz_Fxyz_aa-2.0E0*5*I_ESP_H4yz_Fxyz_a;
    abcd[iGrid*1260+940] = 4.0E0*I_ESP_K5y2z_Fxyz_aa-2.0E0*1*I_ESP_H5y_Fxyz_a-2.0E0*4*I_ESP_H3y2z_Fxyz_a+4*1*I_ESP_F3y_Fxyz;
    abcd[iGrid*1260+941] = 4.0E0*I_ESP_K4y3z_Fxyz_aa-2.0E0*2*I_ESP_H4yz_Fxyz_a-2.0E0*3*I_ESP_H2y3z_Fxyz_a+3*2*I_ESP_F2yz_Fxyz;
    abcd[iGrid*1260+942] = 4.0E0*I_ESP_K3y4z_Fxyz_aa-2.0E0*3*I_ESP_H3y2z_Fxyz_a-2.0E0*2*I_ESP_Hy4z_Fxyz_a+2*3*I_ESP_Fy2z_Fxyz;
    abcd[iGrid*1260+943] = 4.0E0*I_ESP_K2y5z_Fxyz_aa-2.0E0*4*I_ESP_H2y3z_Fxyz_a-2.0E0*1*I_ESP_H5z_Fxyz_a+4*I_ESP_F3z_Fxyz;
    abcd[iGrid*1260+944] = 4.0E0*I_ESP_Ky6z_Fxyz_aa-2.0E0*5*I_ESP_Hy4z_Fxyz_a;
    abcd[iGrid*1260+945] = 4.0E0*I_ESP_K5xyz_Fx2z_aa;
    abcd[iGrid*1260+946] = 4.0E0*I_ESP_K4x2yz_Fx2z_aa-2.0E0*1*I_ESP_H4xz_Fx2z_a;
    abcd[iGrid*1260+947] = 4.0E0*I_ESP_K4xy2z_Fx2z_aa-2.0E0*1*I_ESP_H4xy_Fx2z_a;
    abcd[iGrid*1260+948] = 4.0E0*I_ESP_K3x3yz_Fx2z_aa-2.0E0*2*I_ESP_H3xyz_Fx2z_a;
    abcd[iGrid*1260+949] = 4.0E0*I_ESP_K3x2y2z_Fx2z_aa-2.0E0*1*I_ESP_H3x2y_Fx2z_a-2.0E0*1*I_ESP_H3x2z_Fx2z_a+1*I_ESP_F3x_Fx2z;
    abcd[iGrid*1260+950] = 4.0E0*I_ESP_K3xy3z_Fx2z_aa-2.0E0*2*I_ESP_H3xyz_Fx2z_a;
    abcd[iGrid*1260+951] = 4.0E0*I_ESP_K2x4yz_Fx2z_aa-2.0E0*3*I_ESP_H2x2yz_Fx2z_a;
    abcd[iGrid*1260+952] = 4.0E0*I_ESP_K2x3y2z_Fx2z_aa-2.0E0*1*I_ESP_H2x3y_Fx2z_a-2.0E0*2*I_ESP_H2xy2z_Fx2z_a+2*1*I_ESP_F2xy_Fx2z;
    abcd[iGrid*1260+953] = 4.0E0*I_ESP_K2x2y3z_Fx2z_aa-2.0E0*2*I_ESP_H2x2yz_Fx2z_a-2.0E0*1*I_ESP_H2x3z_Fx2z_a+2*I_ESP_F2xz_Fx2z;
    abcd[iGrid*1260+954] = 4.0E0*I_ESP_K2xy4z_Fx2z_aa-2.0E0*3*I_ESP_H2xy2z_Fx2z_a;
    abcd[iGrid*1260+955] = 4.0E0*I_ESP_Kx5yz_Fx2z_aa-2.0E0*4*I_ESP_Hx3yz_Fx2z_a;
    abcd[iGrid*1260+956] = 4.0E0*I_ESP_Kx4y2z_Fx2z_aa-2.0E0*1*I_ESP_Hx4y_Fx2z_a-2.0E0*3*I_ESP_Hx2y2z_Fx2z_a+3*1*I_ESP_Fx2y_Fx2z;
    abcd[iGrid*1260+957] = 4.0E0*I_ESP_Kx3y3z_Fx2z_aa-2.0E0*2*I_ESP_Hx3yz_Fx2z_a-2.0E0*2*I_ESP_Hxy3z_Fx2z_a+2*2*I_ESP_Fxyz_Fx2z;
    abcd[iGrid*1260+958] = 4.0E0*I_ESP_Kx2y4z_Fx2z_aa-2.0E0*3*I_ESP_Hx2y2z_Fx2z_a-2.0E0*1*I_ESP_Hx4z_Fx2z_a+3*I_ESP_Fx2z_Fx2z;
    abcd[iGrid*1260+959] = 4.0E0*I_ESP_Kxy5z_Fx2z_aa-2.0E0*4*I_ESP_Hxy3z_Fx2z_a;
    abcd[iGrid*1260+960] = 4.0E0*I_ESP_K6yz_Fx2z_aa-2.0E0*5*I_ESP_H4yz_Fx2z_a;
    abcd[iGrid*1260+961] = 4.0E0*I_ESP_K5y2z_Fx2z_aa-2.0E0*1*I_ESP_H5y_Fx2z_a-2.0E0*4*I_ESP_H3y2z_Fx2z_a+4*1*I_ESP_F3y_Fx2z;
    abcd[iGrid*1260+962] = 4.0E0*I_ESP_K4y3z_Fx2z_aa-2.0E0*2*I_ESP_H4yz_Fx2z_a-2.0E0*3*I_ESP_H2y3z_Fx2z_a+3*2*I_ESP_F2yz_Fx2z;
    abcd[iGrid*1260+963] = 4.0E0*I_ESP_K3y4z_Fx2z_aa-2.0E0*3*I_ESP_H3y2z_Fx2z_a-2.0E0*2*I_ESP_Hy4z_Fx2z_a+2*3*I_ESP_Fy2z_Fx2z;
    abcd[iGrid*1260+964] = 4.0E0*I_ESP_K2y5z_Fx2z_aa-2.0E0*4*I_ESP_H2y3z_Fx2z_a-2.0E0*1*I_ESP_H5z_Fx2z_a+4*I_ESP_F3z_Fx2z;
    abcd[iGrid*1260+965] = 4.0E0*I_ESP_Ky6z_Fx2z_aa-2.0E0*5*I_ESP_Hy4z_Fx2z_a;
    abcd[iGrid*1260+966] = 4.0E0*I_ESP_K5xyz_F3y_aa;
    abcd[iGrid*1260+967] = 4.0E0*I_ESP_K4x2yz_F3y_aa-2.0E0*1*I_ESP_H4xz_F3y_a;
    abcd[iGrid*1260+968] = 4.0E0*I_ESP_K4xy2z_F3y_aa-2.0E0*1*I_ESP_H4xy_F3y_a;
    abcd[iGrid*1260+969] = 4.0E0*I_ESP_K3x3yz_F3y_aa-2.0E0*2*I_ESP_H3xyz_F3y_a;
    abcd[iGrid*1260+970] = 4.0E0*I_ESP_K3x2y2z_F3y_aa-2.0E0*1*I_ESP_H3x2y_F3y_a-2.0E0*1*I_ESP_H3x2z_F3y_a+1*I_ESP_F3x_F3y;
    abcd[iGrid*1260+971] = 4.0E0*I_ESP_K3xy3z_F3y_aa-2.0E0*2*I_ESP_H3xyz_F3y_a;
    abcd[iGrid*1260+972] = 4.0E0*I_ESP_K2x4yz_F3y_aa-2.0E0*3*I_ESP_H2x2yz_F3y_a;
    abcd[iGrid*1260+973] = 4.0E0*I_ESP_K2x3y2z_F3y_aa-2.0E0*1*I_ESP_H2x3y_F3y_a-2.0E0*2*I_ESP_H2xy2z_F3y_a+2*1*I_ESP_F2xy_F3y;
    abcd[iGrid*1260+974] = 4.0E0*I_ESP_K2x2y3z_F3y_aa-2.0E0*2*I_ESP_H2x2yz_F3y_a-2.0E0*1*I_ESP_H2x3z_F3y_a+2*I_ESP_F2xz_F3y;
    abcd[iGrid*1260+975] = 4.0E0*I_ESP_K2xy4z_F3y_aa-2.0E0*3*I_ESP_H2xy2z_F3y_a;
    abcd[iGrid*1260+976] = 4.0E0*I_ESP_Kx5yz_F3y_aa-2.0E0*4*I_ESP_Hx3yz_F3y_a;
    abcd[iGrid*1260+977] = 4.0E0*I_ESP_Kx4y2z_F3y_aa-2.0E0*1*I_ESP_Hx4y_F3y_a-2.0E0*3*I_ESP_Hx2y2z_F3y_a+3*1*I_ESP_Fx2y_F3y;
    abcd[iGrid*1260+978] = 4.0E0*I_ESP_Kx3y3z_F3y_aa-2.0E0*2*I_ESP_Hx3yz_F3y_a-2.0E0*2*I_ESP_Hxy3z_F3y_a+2*2*I_ESP_Fxyz_F3y;
    abcd[iGrid*1260+979] = 4.0E0*I_ESP_Kx2y4z_F3y_aa-2.0E0*3*I_ESP_Hx2y2z_F3y_a-2.0E0*1*I_ESP_Hx4z_F3y_a+3*I_ESP_Fx2z_F3y;
    abcd[iGrid*1260+980] = 4.0E0*I_ESP_Kxy5z_F3y_aa-2.0E0*4*I_ESP_Hxy3z_F3y_a;
    abcd[iGrid*1260+981] = 4.0E0*I_ESP_K6yz_F3y_aa-2.0E0*5*I_ESP_H4yz_F3y_a;
    abcd[iGrid*1260+982] = 4.0E0*I_ESP_K5y2z_F3y_aa-2.0E0*1*I_ESP_H5y_F3y_a-2.0E0*4*I_ESP_H3y2z_F3y_a+4*1*I_ESP_F3y_F3y;
    abcd[iGrid*1260+983] = 4.0E0*I_ESP_K4y3z_F3y_aa-2.0E0*2*I_ESP_H4yz_F3y_a-2.0E0*3*I_ESP_H2y3z_F3y_a+3*2*I_ESP_F2yz_F3y;
    abcd[iGrid*1260+984] = 4.0E0*I_ESP_K3y4z_F3y_aa-2.0E0*3*I_ESP_H3y2z_F3y_a-2.0E0*2*I_ESP_Hy4z_F3y_a+2*3*I_ESP_Fy2z_F3y;
    abcd[iGrid*1260+985] = 4.0E0*I_ESP_K2y5z_F3y_aa-2.0E0*4*I_ESP_H2y3z_F3y_a-2.0E0*1*I_ESP_H5z_F3y_a+4*I_ESP_F3z_F3y;
    abcd[iGrid*1260+986] = 4.0E0*I_ESP_Ky6z_F3y_aa-2.0E0*5*I_ESP_Hy4z_F3y_a;
    abcd[iGrid*1260+987] = 4.0E0*I_ESP_K5xyz_F2yz_aa;
    abcd[iGrid*1260+988] = 4.0E0*I_ESP_K4x2yz_F2yz_aa-2.0E0*1*I_ESP_H4xz_F2yz_a;
    abcd[iGrid*1260+989] = 4.0E0*I_ESP_K4xy2z_F2yz_aa-2.0E0*1*I_ESP_H4xy_F2yz_a;
    abcd[iGrid*1260+990] = 4.0E0*I_ESP_K3x3yz_F2yz_aa-2.0E0*2*I_ESP_H3xyz_F2yz_a;
    abcd[iGrid*1260+991] = 4.0E0*I_ESP_K3x2y2z_F2yz_aa-2.0E0*1*I_ESP_H3x2y_F2yz_a-2.0E0*1*I_ESP_H3x2z_F2yz_a+1*I_ESP_F3x_F2yz;
    abcd[iGrid*1260+992] = 4.0E0*I_ESP_K3xy3z_F2yz_aa-2.0E0*2*I_ESP_H3xyz_F2yz_a;
    abcd[iGrid*1260+993] = 4.0E0*I_ESP_K2x4yz_F2yz_aa-2.0E0*3*I_ESP_H2x2yz_F2yz_a;
    abcd[iGrid*1260+994] = 4.0E0*I_ESP_K2x3y2z_F2yz_aa-2.0E0*1*I_ESP_H2x3y_F2yz_a-2.0E0*2*I_ESP_H2xy2z_F2yz_a+2*1*I_ESP_F2xy_F2yz;
    abcd[iGrid*1260+995] = 4.0E0*I_ESP_K2x2y3z_F2yz_aa-2.0E0*2*I_ESP_H2x2yz_F2yz_a-2.0E0*1*I_ESP_H2x3z_F2yz_a+2*I_ESP_F2xz_F2yz;
    abcd[iGrid*1260+996] = 4.0E0*I_ESP_K2xy4z_F2yz_aa-2.0E0*3*I_ESP_H2xy2z_F2yz_a;
    abcd[iGrid*1260+997] = 4.0E0*I_ESP_Kx5yz_F2yz_aa-2.0E0*4*I_ESP_Hx3yz_F2yz_a;
    abcd[iGrid*1260+998] = 4.0E0*I_ESP_Kx4y2z_F2yz_aa-2.0E0*1*I_ESP_Hx4y_F2yz_a-2.0E0*3*I_ESP_Hx2y2z_F2yz_a+3*1*I_ESP_Fx2y_F2yz;
    abcd[iGrid*1260+999] = 4.0E0*I_ESP_Kx3y3z_F2yz_aa-2.0E0*2*I_ESP_Hx3yz_F2yz_a-2.0E0*2*I_ESP_Hxy3z_F2yz_a+2*2*I_ESP_Fxyz_F2yz;
    abcd[iGrid*1260+1000] = 4.0E0*I_ESP_Kx2y4z_F2yz_aa-2.0E0*3*I_ESP_Hx2y2z_F2yz_a-2.0E0*1*I_ESP_Hx4z_F2yz_a+3*I_ESP_Fx2z_F2yz;
    abcd[iGrid*1260+1001] = 4.0E0*I_ESP_Kxy5z_F2yz_aa-2.0E0*4*I_ESP_Hxy3z_F2yz_a;
    abcd[iGrid*1260+1002] = 4.0E0*I_ESP_K6yz_F2yz_aa-2.0E0*5*I_ESP_H4yz_F2yz_a;
    abcd[iGrid*1260+1003] = 4.0E0*I_ESP_K5y2z_F2yz_aa-2.0E0*1*I_ESP_H5y_F2yz_a-2.0E0*4*I_ESP_H3y2z_F2yz_a+4*1*I_ESP_F3y_F2yz;
    abcd[iGrid*1260+1004] = 4.0E0*I_ESP_K4y3z_F2yz_aa-2.0E0*2*I_ESP_H4yz_F2yz_a-2.0E0*3*I_ESP_H2y3z_F2yz_a+3*2*I_ESP_F2yz_F2yz;
    abcd[iGrid*1260+1005] = 4.0E0*I_ESP_K3y4z_F2yz_aa-2.0E0*3*I_ESP_H3y2z_F2yz_a-2.0E0*2*I_ESP_Hy4z_F2yz_a+2*3*I_ESP_Fy2z_F2yz;
    abcd[iGrid*1260+1006] = 4.0E0*I_ESP_K2y5z_F2yz_aa-2.0E0*4*I_ESP_H2y3z_F2yz_a-2.0E0*1*I_ESP_H5z_F2yz_a+4*I_ESP_F3z_F2yz;
    abcd[iGrid*1260+1007] = 4.0E0*I_ESP_Ky6z_F2yz_aa-2.0E0*5*I_ESP_Hy4z_F2yz_a;
    abcd[iGrid*1260+1008] = 4.0E0*I_ESP_K5xyz_Fy2z_aa;
    abcd[iGrid*1260+1009] = 4.0E0*I_ESP_K4x2yz_Fy2z_aa-2.0E0*1*I_ESP_H4xz_Fy2z_a;
    abcd[iGrid*1260+1010] = 4.0E0*I_ESP_K4xy2z_Fy2z_aa-2.0E0*1*I_ESP_H4xy_Fy2z_a;
    abcd[iGrid*1260+1011] = 4.0E0*I_ESP_K3x3yz_Fy2z_aa-2.0E0*2*I_ESP_H3xyz_Fy2z_a;
    abcd[iGrid*1260+1012] = 4.0E0*I_ESP_K3x2y2z_Fy2z_aa-2.0E0*1*I_ESP_H3x2y_Fy2z_a-2.0E0*1*I_ESP_H3x2z_Fy2z_a+1*I_ESP_F3x_Fy2z;
    abcd[iGrid*1260+1013] = 4.0E0*I_ESP_K3xy3z_Fy2z_aa-2.0E0*2*I_ESP_H3xyz_Fy2z_a;
    abcd[iGrid*1260+1014] = 4.0E0*I_ESP_K2x4yz_Fy2z_aa-2.0E0*3*I_ESP_H2x2yz_Fy2z_a;
    abcd[iGrid*1260+1015] = 4.0E0*I_ESP_K2x3y2z_Fy2z_aa-2.0E0*1*I_ESP_H2x3y_Fy2z_a-2.0E0*2*I_ESP_H2xy2z_Fy2z_a+2*1*I_ESP_F2xy_Fy2z;
    abcd[iGrid*1260+1016] = 4.0E0*I_ESP_K2x2y3z_Fy2z_aa-2.0E0*2*I_ESP_H2x2yz_Fy2z_a-2.0E0*1*I_ESP_H2x3z_Fy2z_a+2*I_ESP_F2xz_Fy2z;
    abcd[iGrid*1260+1017] = 4.0E0*I_ESP_K2xy4z_Fy2z_aa-2.0E0*3*I_ESP_H2xy2z_Fy2z_a;
    abcd[iGrid*1260+1018] = 4.0E0*I_ESP_Kx5yz_Fy2z_aa-2.0E0*4*I_ESP_Hx3yz_Fy2z_a;
    abcd[iGrid*1260+1019] = 4.0E0*I_ESP_Kx4y2z_Fy2z_aa-2.0E0*1*I_ESP_Hx4y_Fy2z_a-2.0E0*3*I_ESP_Hx2y2z_Fy2z_a+3*1*I_ESP_Fx2y_Fy2z;
    abcd[iGrid*1260+1020] = 4.0E0*I_ESP_Kx3y3z_Fy2z_aa-2.0E0*2*I_ESP_Hx3yz_Fy2z_a-2.0E0*2*I_ESP_Hxy3z_Fy2z_a+2*2*I_ESP_Fxyz_Fy2z;
    abcd[iGrid*1260+1021] = 4.0E0*I_ESP_Kx2y4z_Fy2z_aa-2.0E0*3*I_ESP_Hx2y2z_Fy2z_a-2.0E0*1*I_ESP_Hx4z_Fy2z_a+3*I_ESP_Fx2z_Fy2z;
    abcd[iGrid*1260+1022] = 4.0E0*I_ESP_Kxy5z_Fy2z_aa-2.0E0*4*I_ESP_Hxy3z_Fy2z_a;
    abcd[iGrid*1260+1023] = 4.0E0*I_ESP_K6yz_Fy2z_aa-2.0E0*5*I_ESP_H4yz_Fy2z_a;
    abcd[iGrid*1260+1024] = 4.0E0*I_ESP_K5y2z_Fy2z_aa-2.0E0*1*I_ESP_H5y_Fy2z_a-2.0E0*4*I_ESP_H3y2z_Fy2z_a+4*1*I_ESP_F3y_Fy2z;
    abcd[iGrid*1260+1025] = 4.0E0*I_ESP_K4y3z_Fy2z_aa-2.0E0*2*I_ESP_H4yz_Fy2z_a-2.0E0*3*I_ESP_H2y3z_Fy2z_a+3*2*I_ESP_F2yz_Fy2z;
    abcd[iGrid*1260+1026] = 4.0E0*I_ESP_K3y4z_Fy2z_aa-2.0E0*3*I_ESP_H3y2z_Fy2z_a-2.0E0*2*I_ESP_Hy4z_Fy2z_a+2*3*I_ESP_Fy2z_Fy2z;
    abcd[iGrid*1260+1027] = 4.0E0*I_ESP_K2y5z_Fy2z_aa-2.0E0*4*I_ESP_H2y3z_Fy2z_a-2.0E0*1*I_ESP_H5z_Fy2z_a+4*I_ESP_F3z_Fy2z;
    abcd[iGrid*1260+1028] = 4.0E0*I_ESP_Ky6z_Fy2z_aa-2.0E0*5*I_ESP_Hy4z_Fy2z_a;
    abcd[iGrid*1260+1029] = 4.0E0*I_ESP_K5xyz_F3z_aa;
    abcd[iGrid*1260+1030] = 4.0E0*I_ESP_K4x2yz_F3z_aa-2.0E0*1*I_ESP_H4xz_F3z_a;
    abcd[iGrid*1260+1031] = 4.0E0*I_ESP_K4xy2z_F3z_aa-2.0E0*1*I_ESP_H4xy_F3z_a;
    abcd[iGrid*1260+1032] = 4.0E0*I_ESP_K3x3yz_F3z_aa-2.0E0*2*I_ESP_H3xyz_F3z_a;
    abcd[iGrid*1260+1033] = 4.0E0*I_ESP_K3x2y2z_F3z_aa-2.0E0*1*I_ESP_H3x2y_F3z_a-2.0E0*1*I_ESP_H3x2z_F3z_a+1*I_ESP_F3x_F3z;
    abcd[iGrid*1260+1034] = 4.0E0*I_ESP_K3xy3z_F3z_aa-2.0E0*2*I_ESP_H3xyz_F3z_a;
    abcd[iGrid*1260+1035] = 4.0E0*I_ESP_K2x4yz_F3z_aa-2.0E0*3*I_ESP_H2x2yz_F3z_a;
    abcd[iGrid*1260+1036] = 4.0E0*I_ESP_K2x3y2z_F3z_aa-2.0E0*1*I_ESP_H2x3y_F3z_a-2.0E0*2*I_ESP_H2xy2z_F3z_a+2*1*I_ESP_F2xy_F3z;
    abcd[iGrid*1260+1037] = 4.0E0*I_ESP_K2x2y3z_F3z_aa-2.0E0*2*I_ESP_H2x2yz_F3z_a-2.0E0*1*I_ESP_H2x3z_F3z_a+2*I_ESP_F2xz_F3z;
    abcd[iGrid*1260+1038] = 4.0E0*I_ESP_K2xy4z_F3z_aa-2.0E0*3*I_ESP_H2xy2z_F3z_a;
    abcd[iGrid*1260+1039] = 4.0E0*I_ESP_Kx5yz_F3z_aa-2.0E0*4*I_ESP_Hx3yz_F3z_a;
    abcd[iGrid*1260+1040] = 4.0E0*I_ESP_Kx4y2z_F3z_aa-2.0E0*1*I_ESP_Hx4y_F3z_a-2.0E0*3*I_ESP_Hx2y2z_F3z_a+3*1*I_ESP_Fx2y_F3z;
    abcd[iGrid*1260+1041] = 4.0E0*I_ESP_Kx3y3z_F3z_aa-2.0E0*2*I_ESP_Hx3yz_F3z_a-2.0E0*2*I_ESP_Hxy3z_F3z_a+2*2*I_ESP_Fxyz_F3z;
    abcd[iGrid*1260+1042] = 4.0E0*I_ESP_Kx2y4z_F3z_aa-2.0E0*3*I_ESP_Hx2y2z_F3z_a-2.0E0*1*I_ESP_Hx4z_F3z_a+3*I_ESP_Fx2z_F3z;
    abcd[iGrid*1260+1043] = 4.0E0*I_ESP_Kxy5z_F3z_aa-2.0E0*4*I_ESP_Hxy3z_F3z_a;
    abcd[iGrid*1260+1044] = 4.0E0*I_ESP_K6yz_F3z_aa-2.0E0*5*I_ESP_H4yz_F3z_a;
    abcd[iGrid*1260+1045] = 4.0E0*I_ESP_K5y2z_F3z_aa-2.0E0*1*I_ESP_H5y_F3z_a-2.0E0*4*I_ESP_H3y2z_F3z_a+4*1*I_ESP_F3y_F3z;
    abcd[iGrid*1260+1046] = 4.0E0*I_ESP_K4y3z_F3z_aa-2.0E0*2*I_ESP_H4yz_F3z_a-2.0E0*3*I_ESP_H2y3z_F3z_a+3*2*I_ESP_F2yz_F3z;
    abcd[iGrid*1260+1047] = 4.0E0*I_ESP_K3y4z_F3z_aa-2.0E0*3*I_ESP_H3y2z_F3z_a-2.0E0*2*I_ESP_Hy4z_F3z_a+2*3*I_ESP_Fy2z_F3z;
    abcd[iGrid*1260+1048] = 4.0E0*I_ESP_K2y5z_F3z_aa-2.0E0*4*I_ESP_H2y3z_F3z_a-2.0E0*1*I_ESP_H5z_F3z_a+4*I_ESP_F3z_F3z;
    abcd[iGrid*1260+1049] = 4.0E0*I_ESP_Ky6z_F3z_aa-2.0E0*5*I_ESP_Hy4z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_F_aa
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_H_F_a
     * RHS shell quartet name: SQ_ESP_F_F
     ************************************************************/
    abcd[iGrid*1260+1050] = 4.0E0*I_ESP_K5x2z_F3x_aa-2.0E0*1*I_ESP_H5x_F3x_a;
    abcd[iGrid*1260+1051] = 4.0E0*I_ESP_K4xy2z_F3x_aa-2.0E0*1*I_ESP_H4xy_F3x_a;
    abcd[iGrid*1260+1052] = 4.0E0*I_ESP_K4x3z_F3x_aa-2.0E0*1*I_ESP_H4xz_F3x_a-2.0E0*2*I_ESP_H4xz_F3x_a;
    abcd[iGrid*1260+1053] = 4.0E0*I_ESP_K3x2y2z_F3x_aa-2.0E0*1*I_ESP_H3x2y_F3x_a;
    abcd[iGrid*1260+1054] = 4.0E0*I_ESP_K3xy3z_F3x_aa-2.0E0*1*I_ESP_H3xyz_F3x_a-2.0E0*2*I_ESP_H3xyz_F3x_a;
    abcd[iGrid*1260+1055] = 4.0E0*I_ESP_K3x4z_F3x_aa-2.0E0*2*I_ESP_H3x2z_F3x_a-2.0E0*3*I_ESP_H3x2z_F3x_a+2*1*I_ESP_F3x_F3x;
    abcd[iGrid*1260+1056] = 4.0E0*I_ESP_K2x3y2z_F3x_aa-2.0E0*1*I_ESP_H2x3y_F3x_a;
    abcd[iGrid*1260+1057] = 4.0E0*I_ESP_K2x2y3z_F3x_aa-2.0E0*1*I_ESP_H2x2yz_F3x_a-2.0E0*2*I_ESP_H2x2yz_F3x_a;
    abcd[iGrid*1260+1058] = 4.0E0*I_ESP_K2xy4z_F3x_aa-2.0E0*2*I_ESP_H2xy2z_F3x_a-2.0E0*3*I_ESP_H2xy2z_F3x_a+2*1*I_ESP_F2xy_F3x;
    abcd[iGrid*1260+1059] = 4.0E0*I_ESP_K2x5z_F3x_aa-2.0E0*3*I_ESP_H2x3z_F3x_a-2.0E0*4*I_ESP_H2x3z_F3x_a+3*2*I_ESP_F2xz_F3x;
    abcd[iGrid*1260+1060] = 4.0E0*I_ESP_Kx4y2z_F3x_aa-2.0E0*1*I_ESP_Hx4y_F3x_a;
    abcd[iGrid*1260+1061] = 4.0E0*I_ESP_Kx3y3z_F3x_aa-2.0E0*1*I_ESP_Hx3yz_F3x_a-2.0E0*2*I_ESP_Hx3yz_F3x_a;
    abcd[iGrid*1260+1062] = 4.0E0*I_ESP_Kx2y4z_F3x_aa-2.0E0*2*I_ESP_Hx2y2z_F3x_a-2.0E0*3*I_ESP_Hx2y2z_F3x_a+2*1*I_ESP_Fx2y_F3x;
    abcd[iGrid*1260+1063] = 4.0E0*I_ESP_Kxy5z_F3x_aa-2.0E0*3*I_ESP_Hxy3z_F3x_a-2.0E0*4*I_ESP_Hxy3z_F3x_a+3*2*I_ESP_Fxyz_F3x;
    abcd[iGrid*1260+1064] = 4.0E0*I_ESP_Kx6z_F3x_aa-2.0E0*4*I_ESP_Hx4z_F3x_a-2.0E0*5*I_ESP_Hx4z_F3x_a+4*3*I_ESP_Fx2z_F3x;
    abcd[iGrid*1260+1065] = 4.0E0*I_ESP_K5y2z_F3x_aa-2.0E0*1*I_ESP_H5y_F3x_a;
    abcd[iGrid*1260+1066] = 4.0E0*I_ESP_K4y3z_F3x_aa-2.0E0*1*I_ESP_H4yz_F3x_a-2.0E0*2*I_ESP_H4yz_F3x_a;
    abcd[iGrid*1260+1067] = 4.0E0*I_ESP_K3y4z_F3x_aa-2.0E0*2*I_ESP_H3y2z_F3x_a-2.0E0*3*I_ESP_H3y2z_F3x_a+2*1*I_ESP_F3y_F3x;
    abcd[iGrid*1260+1068] = 4.0E0*I_ESP_K2y5z_F3x_aa-2.0E0*3*I_ESP_H2y3z_F3x_a-2.0E0*4*I_ESP_H2y3z_F3x_a+3*2*I_ESP_F2yz_F3x;
    abcd[iGrid*1260+1069] = 4.0E0*I_ESP_Ky6z_F3x_aa-2.0E0*4*I_ESP_Hy4z_F3x_a-2.0E0*5*I_ESP_Hy4z_F3x_a+4*3*I_ESP_Fy2z_F3x;
    abcd[iGrid*1260+1070] = 4.0E0*I_ESP_K7z_F3x_aa-2.0E0*5*I_ESP_H5z_F3x_a-2.0E0*6*I_ESP_H5z_F3x_a+5*4*I_ESP_F3z_F3x;
    abcd[iGrid*1260+1071] = 4.0E0*I_ESP_K5x2z_F2xy_aa-2.0E0*1*I_ESP_H5x_F2xy_a;
    abcd[iGrid*1260+1072] = 4.0E0*I_ESP_K4xy2z_F2xy_aa-2.0E0*1*I_ESP_H4xy_F2xy_a;
    abcd[iGrid*1260+1073] = 4.0E0*I_ESP_K4x3z_F2xy_aa-2.0E0*1*I_ESP_H4xz_F2xy_a-2.0E0*2*I_ESP_H4xz_F2xy_a;
    abcd[iGrid*1260+1074] = 4.0E0*I_ESP_K3x2y2z_F2xy_aa-2.0E0*1*I_ESP_H3x2y_F2xy_a;
    abcd[iGrid*1260+1075] = 4.0E0*I_ESP_K3xy3z_F2xy_aa-2.0E0*1*I_ESP_H3xyz_F2xy_a-2.0E0*2*I_ESP_H3xyz_F2xy_a;
    abcd[iGrid*1260+1076] = 4.0E0*I_ESP_K3x4z_F2xy_aa-2.0E0*2*I_ESP_H3x2z_F2xy_a-2.0E0*3*I_ESP_H3x2z_F2xy_a+2*1*I_ESP_F3x_F2xy;
    abcd[iGrid*1260+1077] = 4.0E0*I_ESP_K2x3y2z_F2xy_aa-2.0E0*1*I_ESP_H2x3y_F2xy_a;
    abcd[iGrid*1260+1078] = 4.0E0*I_ESP_K2x2y3z_F2xy_aa-2.0E0*1*I_ESP_H2x2yz_F2xy_a-2.0E0*2*I_ESP_H2x2yz_F2xy_a;
    abcd[iGrid*1260+1079] = 4.0E0*I_ESP_K2xy4z_F2xy_aa-2.0E0*2*I_ESP_H2xy2z_F2xy_a-2.0E0*3*I_ESP_H2xy2z_F2xy_a+2*1*I_ESP_F2xy_F2xy;
    abcd[iGrid*1260+1080] = 4.0E0*I_ESP_K2x5z_F2xy_aa-2.0E0*3*I_ESP_H2x3z_F2xy_a-2.0E0*4*I_ESP_H2x3z_F2xy_a+3*2*I_ESP_F2xz_F2xy;
    abcd[iGrid*1260+1081] = 4.0E0*I_ESP_Kx4y2z_F2xy_aa-2.0E0*1*I_ESP_Hx4y_F2xy_a;
    abcd[iGrid*1260+1082] = 4.0E0*I_ESP_Kx3y3z_F2xy_aa-2.0E0*1*I_ESP_Hx3yz_F2xy_a-2.0E0*2*I_ESP_Hx3yz_F2xy_a;
    abcd[iGrid*1260+1083] = 4.0E0*I_ESP_Kx2y4z_F2xy_aa-2.0E0*2*I_ESP_Hx2y2z_F2xy_a-2.0E0*3*I_ESP_Hx2y2z_F2xy_a+2*1*I_ESP_Fx2y_F2xy;
    abcd[iGrid*1260+1084] = 4.0E0*I_ESP_Kxy5z_F2xy_aa-2.0E0*3*I_ESP_Hxy3z_F2xy_a-2.0E0*4*I_ESP_Hxy3z_F2xy_a+3*2*I_ESP_Fxyz_F2xy;
    abcd[iGrid*1260+1085] = 4.0E0*I_ESP_Kx6z_F2xy_aa-2.0E0*4*I_ESP_Hx4z_F2xy_a-2.0E0*5*I_ESP_Hx4z_F2xy_a+4*3*I_ESP_Fx2z_F2xy;
    abcd[iGrid*1260+1086] = 4.0E0*I_ESP_K5y2z_F2xy_aa-2.0E0*1*I_ESP_H5y_F2xy_a;
    abcd[iGrid*1260+1087] = 4.0E0*I_ESP_K4y3z_F2xy_aa-2.0E0*1*I_ESP_H4yz_F2xy_a-2.0E0*2*I_ESP_H4yz_F2xy_a;
    abcd[iGrid*1260+1088] = 4.0E0*I_ESP_K3y4z_F2xy_aa-2.0E0*2*I_ESP_H3y2z_F2xy_a-2.0E0*3*I_ESP_H3y2z_F2xy_a+2*1*I_ESP_F3y_F2xy;
    abcd[iGrid*1260+1089] = 4.0E0*I_ESP_K2y5z_F2xy_aa-2.0E0*3*I_ESP_H2y3z_F2xy_a-2.0E0*4*I_ESP_H2y3z_F2xy_a+3*2*I_ESP_F2yz_F2xy;
    abcd[iGrid*1260+1090] = 4.0E0*I_ESP_Ky6z_F2xy_aa-2.0E0*4*I_ESP_Hy4z_F2xy_a-2.0E0*5*I_ESP_Hy4z_F2xy_a+4*3*I_ESP_Fy2z_F2xy;
    abcd[iGrid*1260+1091] = 4.0E0*I_ESP_K7z_F2xy_aa-2.0E0*5*I_ESP_H5z_F2xy_a-2.0E0*6*I_ESP_H5z_F2xy_a+5*4*I_ESP_F3z_F2xy;
    abcd[iGrid*1260+1092] = 4.0E0*I_ESP_K5x2z_F2xz_aa-2.0E0*1*I_ESP_H5x_F2xz_a;
    abcd[iGrid*1260+1093] = 4.0E0*I_ESP_K4xy2z_F2xz_aa-2.0E0*1*I_ESP_H4xy_F2xz_a;
    abcd[iGrid*1260+1094] = 4.0E0*I_ESP_K4x3z_F2xz_aa-2.0E0*1*I_ESP_H4xz_F2xz_a-2.0E0*2*I_ESP_H4xz_F2xz_a;
    abcd[iGrid*1260+1095] = 4.0E0*I_ESP_K3x2y2z_F2xz_aa-2.0E0*1*I_ESP_H3x2y_F2xz_a;
    abcd[iGrid*1260+1096] = 4.0E0*I_ESP_K3xy3z_F2xz_aa-2.0E0*1*I_ESP_H3xyz_F2xz_a-2.0E0*2*I_ESP_H3xyz_F2xz_a;
    abcd[iGrid*1260+1097] = 4.0E0*I_ESP_K3x4z_F2xz_aa-2.0E0*2*I_ESP_H3x2z_F2xz_a-2.0E0*3*I_ESP_H3x2z_F2xz_a+2*1*I_ESP_F3x_F2xz;
    abcd[iGrid*1260+1098] = 4.0E0*I_ESP_K2x3y2z_F2xz_aa-2.0E0*1*I_ESP_H2x3y_F2xz_a;
    abcd[iGrid*1260+1099] = 4.0E0*I_ESP_K2x2y3z_F2xz_aa-2.0E0*1*I_ESP_H2x2yz_F2xz_a-2.0E0*2*I_ESP_H2x2yz_F2xz_a;
    abcd[iGrid*1260+1100] = 4.0E0*I_ESP_K2xy4z_F2xz_aa-2.0E0*2*I_ESP_H2xy2z_F2xz_a-2.0E0*3*I_ESP_H2xy2z_F2xz_a+2*1*I_ESP_F2xy_F2xz;
    abcd[iGrid*1260+1101] = 4.0E0*I_ESP_K2x5z_F2xz_aa-2.0E0*3*I_ESP_H2x3z_F2xz_a-2.0E0*4*I_ESP_H2x3z_F2xz_a+3*2*I_ESP_F2xz_F2xz;
    abcd[iGrid*1260+1102] = 4.0E0*I_ESP_Kx4y2z_F2xz_aa-2.0E0*1*I_ESP_Hx4y_F2xz_a;
    abcd[iGrid*1260+1103] = 4.0E0*I_ESP_Kx3y3z_F2xz_aa-2.0E0*1*I_ESP_Hx3yz_F2xz_a-2.0E0*2*I_ESP_Hx3yz_F2xz_a;
    abcd[iGrid*1260+1104] = 4.0E0*I_ESP_Kx2y4z_F2xz_aa-2.0E0*2*I_ESP_Hx2y2z_F2xz_a-2.0E0*3*I_ESP_Hx2y2z_F2xz_a+2*1*I_ESP_Fx2y_F2xz;
    abcd[iGrid*1260+1105] = 4.0E0*I_ESP_Kxy5z_F2xz_aa-2.0E0*3*I_ESP_Hxy3z_F2xz_a-2.0E0*4*I_ESP_Hxy3z_F2xz_a+3*2*I_ESP_Fxyz_F2xz;
    abcd[iGrid*1260+1106] = 4.0E0*I_ESP_Kx6z_F2xz_aa-2.0E0*4*I_ESP_Hx4z_F2xz_a-2.0E0*5*I_ESP_Hx4z_F2xz_a+4*3*I_ESP_Fx2z_F2xz;
    abcd[iGrid*1260+1107] = 4.0E0*I_ESP_K5y2z_F2xz_aa-2.0E0*1*I_ESP_H5y_F2xz_a;
    abcd[iGrid*1260+1108] = 4.0E0*I_ESP_K4y3z_F2xz_aa-2.0E0*1*I_ESP_H4yz_F2xz_a-2.0E0*2*I_ESP_H4yz_F2xz_a;
    abcd[iGrid*1260+1109] = 4.0E0*I_ESP_K3y4z_F2xz_aa-2.0E0*2*I_ESP_H3y2z_F2xz_a-2.0E0*3*I_ESP_H3y2z_F2xz_a+2*1*I_ESP_F3y_F2xz;
    abcd[iGrid*1260+1110] = 4.0E0*I_ESP_K2y5z_F2xz_aa-2.0E0*3*I_ESP_H2y3z_F2xz_a-2.0E0*4*I_ESP_H2y3z_F2xz_a+3*2*I_ESP_F2yz_F2xz;
    abcd[iGrid*1260+1111] = 4.0E0*I_ESP_Ky6z_F2xz_aa-2.0E0*4*I_ESP_Hy4z_F2xz_a-2.0E0*5*I_ESP_Hy4z_F2xz_a+4*3*I_ESP_Fy2z_F2xz;
    abcd[iGrid*1260+1112] = 4.0E0*I_ESP_K7z_F2xz_aa-2.0E0*5*I_ESP_H5z_F2xz_a-2.0E0*6*I_ESP_H5z_F2xz_a+5*4*I_ESP_F3z_F2xz;
    abcd[iGrid*1260+1113] = 4.0E0*I_ESP_K5x2z_Fx2y_aa-2.0E0*1*I_ESP_H5x_Fx2y_a;
    abcd[iGrid*1260+1114] = 4.0E0*I_ESP_K4xy2z_Fx2y_aa-2.0E0*1*I_ESP_H4xy_Fx2y_a;
    abcd[iGrid*1260+1115] = 4.0E0*I_ESP_K4x3z_Fx2y_aa-2.0E0*1*I_ESP_H4xz_Fx2y_a-2.0E0*2*I_ESP_H4xz_Fx2y_a;
    abcd[iGrid*1260+1116] = 4.0E0*I_ESP_K3x2y2z_Fx2y_aa-2.0E0*1*I_ESP_H3x2y_Fx2y_a;
    abcd[iGrid*1260+1117] = 4.0E0*I_ESP_K3xy3z_Fx2y_aa-2.0E0*1*I_ESP_H3xyz_Fx2y_a-2.0E0*2*I_ESP_H3xyz_Fx2y_a;
    abcd[iGrid*1260+1118] = 4.0E0*I_ESP_K3x4z_Fx2y_aa-2.0E0*2*I_ESP_H3x2z_Fx2y_a-2.0E0*3*I_ESP_H3x2z_Fx2y_a+2*1*I_ESP_F3x_Fx2y;
    abcd[iGrid*1260+1119] = 4.0E0*I_ESP_K2x3y2z_Fx2y_aa-2.0E0*1*I_ESP_H2x3y_Fx2y_a;
    abcd[iGrid*1260+1120] = 4.0E0*I_ESP_K2x2y3z_Fx2y_aa-2.0E0*1*I_ESP_H2x2yz_Fx2y_a-2.0E0*2*I_ESP_H2x2yz_Fx2y_a;
    abcd[iGrid*1260+1121] = 4.0E0*I_ESP_K2xy4z_Fx2y_aa-2.0E0*2*I_ESP_H2xy2z_Fx2y_a-2.0E0*3*I_ESP_H2xy2z_Fx2y_a+2*1*I_ESP_F2xy_Fx2y;
    abcd[iGrid*1260+1122] = 4.0E0*I_ESP_K2x5z_Fx2y_aa-2.0E0*3*I_ESP_H2x3z_Fx2y_a-2.0E0*4*I_ESP_H2x3z_Fx2y_a+3*2*I_ESP_F2xz_Fx2y;
    abcd[iGrid*1260+1123] = 4.0E0*I_ESP_Kx4y2z_Fx2y_aa-2.0E0*1*I_ESP_Hx4y_Fx2y_a;
    abcd[iGrid*1260+1124] = 4.0E0*I_ESP_Kx3y3z_Fx2y_aa-2.0E0*1*I_ESP_Hx3yz_Fx2y_a-2.0E0*2*I_ESP_Hx3yz_Fx2y_a;
    abcd[iGrid*1260+1125] = 4.0E0*I_ESP_Kx2y4z_Fx2y_aa-2.0E0*2*I_ESP_Hx2y2z_Fx2y_a-2.0E0*3*I_ESP_Hx2y2z_Fx2y_a+2*1*I_ESP_Fx2y_Fx2y;
    abcd[iGrid*1260+1126] = 4.0E0*I_ESP_Kxy5z_Fx2y_aa-2.0E0*3*I_ESP_Hxy3z_Fx2y_a-2.0E0*4*I_ESP_Hxy3z_Fx2y_a+3*2*I_ESP_Fxyz_Fx2y;
    abcd[iGrid*1260+1127] = 4.0E0*I_ESP_Kx6z_Fx2y_aa-2.0E0*4*I_ESP_Hx4z_Fx2y_a-2.0E0*5*I_ESP_Hx4z_Fx2y_a+4*3*I_ESP_Fx2z_Fx2y;
    abcd[iGrid*1260+1128] = 4.0E0*I_ESP_K5y2z_Fx2y_aa-2.0E0*1*I_ESP_H5y_Fx2y_a;
    abcd[iGrid*1260+1129] = 4.0E0*I_ESP_K4y3z_Fx2y_aa-2.0E0*1*I_ESP_H4yz_Fx2y_a-2.0E0*2*I_ESP_H4yz_Fx2y_a;
    abcd[iGrid*1260+1130] = 4.0E0*I_ESP_K3y4z_Fx2y_aa-2.0E0*2*I_ESP_H3y2z_Fx2y_a-2.0E0*3*I_ESP_H3y2z_Fx2y_a+2*1*I_ESP_F3y_Fx2y;
    abcd[iGrid*1260+1131] = 4.0E0*I_ESP_K2y5z_Fx2y_aa-2.0E0*3*I_ESP_H2y3z_Fx2y_a-2.0E0*4*I_ESP_H2y3z_Fx2y_a+3*2*I_ESP_F2yz_Fx2y;
    abcd[iGrid*1260+1132] = 4.0E0*I_ESP_Ky6z_Fx2y_aa-2.0E0*4*I_ESP_Hy4z_Fx2y_a-2.0E0*5*I_ESP_Hy4z_Fx2y_a+4*3*I_ESP_Fy2z_Fx2y;
    abcd[iGrid*1260+1133] = 4.0E0*I_ESP_K7z_Fx2y_aa-2.0E0*5*I_ESP_H5z_Fx2y_a-2.0E0*6*I_ESP_H5z_Fx2y_a+5*4*I_ESP_F3z_Fx2y;
    abcd[iGrid*1260+1134] = 4.0E0*I_ESP_K5x2z_Fxyz_aa-2.0E0*1*I_ESP_H5x_Fxyz_a;
    abcd[iGrid*1260+1135] = 4.0E0*I_ESP_K4xy2z_Fxyz_aa-2.0E0*1*I_ESP_H4xy_Fxyz_a;
    abcd[iGrid*1260+1136] = 4.0E0*I_ESP_K4x3z_Fxyz_aa-2.0E0*1*I_ESP_H4xz_Fxyz_a-2.0E0*2*I_ESP_H4xz_Fxyz_a;
    abcd[iGrid*1260+1137] = 4.0E0*I_ESP_K3x2y2z_Fxyz_aa-2.0E0*1*I_ESP_H3x2y_Fxyz_a;
    abcd[iGrid*1260+1138] = 4.0E0*I_ESP_K3xy3z_Fxyz_aa-2.0E0*1*I_ESP_H3xyz_Fxyz_a-2.0E0*2*I_ESP_H3xyz_Fxyz_a;
    abcd[iGrid*1260+1139] = 4.0E0*I_ESP_K3x4z_Fxyz_aa-2.0E0*2*I_ESP_H3x2z_Fxyz_a-2.0E0*3*I_ESP_H3x2z_Fxyz_a+2*1*I_ESP_F3x_Fxyz;
    abcd[iGrid*1260+1140] = 4.0E0*I_ESP_K2x3y2z_Fxyz_aa-2.0E0*1*I_ESP_H2x3y_Fxyz_a;
    abcd[iGrid*1260+1141] = 4.0E0*I_ESP_K2x2y3z_Fxyz_aa-2.0E0*1*I_ESP_H2x2yz_Fxyz_a-2.0E0*2*I_ESP_H2x2yz_Fxyz_a;
    abcd[iGrid*1260+1142] = 4.0E0*I_ESP_K2xy4z_Fxyz_aa-2.0E0*2*I_ESP_H2xy2z_Fxyz_a-2.0E0*3*I_ESP_H2xy2z_Fxyz_a+2*1*I_ESP_F2xy_Fxyz;
    abcd[iGrid*1260+1143] = 4.0E0*I_ESP_K2x5z_Fxyz_aa-2.0E0*3*I_ESP_H2x3z_Fxyz_a-2.0E0*4*I_ESP_H2x3z_Fxyz_a+3*2*I_ESP_F2xz_Fxyz;
    abcd[iGrid*1260+1144] = 4.0E0*I_ESP_Kx4y2z_Fxyz_aa-2.0E0*1*I_ESP_Hx4y_Fxyz_a;
    abcd[iGrid*1260+1145] = 4.0E0*I_ESP_Kx3y3z_Fxyz_aa-2.0E0*1*I_ESP_Hx3yz_Fxyz_a-2.0E0*2*I_ESP_Hx3yz_Fxyz_a;
    abcd[iGrid*1260+1146] = 4.0E0*I_ESP_Kx2y4z_Fxyz_aa-2.0E0*2*I_ESP_Hx2y2z_Fxyz_a-2.0E0*3*I_ESP_Hx2y2z_Fxyz_a+2*1*I_ESP_Fx2y_Fxyz;
    abcd[iGrid*1260+1147] = 4.0E0*I_ESP_Kxy5z_Fxyz_aa-2.0E0*3*I_ESP_Hxy3z_Fxyz_a-2.0E0*4*I_ESP_Hxy3z_Fxyz_a+3*2*I_ESP_Fxyz_Fxyz;
    abcd[iGrid*1260+1148] = 4.0E0*I_ESP_Kx6z_Fxyz_aa-2.0E0*4*I_ESP_Hx4z_Fxyz_a-2.0E0*5*I_ESP_Hx4z_Fxyz_a+4*3*I_ESP_Fx2z_Fxyz;
    abcd[iGrid*1260+1149] = 4.0E0*I_ESP_K5y2z_Fxyz_aa-2.0E0*1*I_ESP_H5y_Fxyz_a;
    abcd[iGrid*1260+1150] = 4.0E0*I_ESP_K4y3z_Fxyz_aa-2.0E0*1*I_ESP_H4yz_Fxyz_a-2.0E0*2*I_ESP_H4yz_Fxyz_a;
    abcd[iGrid*1260+1151] = 4.0E0*I_ESP_K3y4z_Fxyz_aa-2.0E0*2*I_ESP_H3y2z_Fxyz_a-2.0E0*3*I_ESP_H3y2z_Fxyz_a+2*1*I_ESP_F3y_Fxyz;
    abcd[iGrid*1260+1152] = 4.0E0*I_ESP_K2y5z_Fxyz_aa-2.0E0*3*I_ESP_H2y3z_Fxyz_a-2.0E0*4*I_ESP_H2y3z_Fxyz_a+3*2*I_ESP_F2yz_Fxyz;
    abcd[iGrid*1260+1153] = 4.0E0*I_ESP_Ky6z_Fxyz_aa-2.0E0*4*I_ESP_Hy4z_Fxyz_a-2.0E0*5*I_ESP_Hy4z_Fxyz_a+4*3*I_ESP_Fy2z_Fxyz;
    abcd[iGrid*1260+1154] = 4.0E0*I_ESP_K7z_Fxyz_aa-2.0E0*5*I_ESP_H5z_Fxyz_a-2.0E0*6*I_ESP_H5z_Fxyz_a+5*4*I_ESP_F3z_Fxyz;
    abcd[iGrid*1260+1155] = 4.0E0*I_ESP_K5x2z_Fx2z_aa-2.0E0*1*I_ESP_H5x_Fx2z_a;
    abcd[iGrid*1260+1156] = 4.0E0*I_ESP_K4xy2z_Fx2z_aa-2.0E0*1*I_ESP_H4xy_Fx2z_a;
    abcd[iGrid*1260+1157] = 4.0E0*I_ESP_K4x3z_Fx2z_aa-2.0E0*1*I_ESP_H4xz_Fx2z_a-2.0E0*2*I_ESP_H4xz_Fx2z_a;
    abcd[iGrid*1260+1158] = 4.0E0*I_ESP_K3x2y2z_Fx2z_aa-2.0E0*1*I_ESP_H3x2y_Fx2z_a;
    abcd[iGrid*1260+1159] = 4.0E0*I_ESP_K3xy3z_Fx2z_aa-2.0E0*1*I_ESP_H3xyz_Fx2z_a-2.0E0*2*I_ESP_H3xyz_Fx2z_a;
    abcd[iGrid*1260+1160] = 4.0E0*I_ESP_K3x4z_Fx2z_aa-2.0E0*2*I_ESP_H3x2z_Fx2z_a-2.0E0*3*I_ESP_H3x2z_Fx2z_a+2*1*I_ESP_F3x_Fx2z;
    abcd[iGrid*1260+1161] = 4.0E0*I_ESP_K2x3y2z_Fx2z_aa-2.0E0*1*I_ESP_H2x3y_Fx2z_a;
    abcd[iGrid*1260+1162] = 4.0E0*I_ESP_K2x2y3z_Fx2z_aa-2.0E0*1*I_ESP_H2x2yz_Fx2z_a-2.0E0*2*I_ESP_H2x2yz_Fx2z_a;
    abcd[iGrid*1260+1163] = 4.0E0*I_ESP_K2xy4z_Fx2z_aa-2.0E0*2*I_ESP_H2xy2z_Fx2z_a-2.0E0*3*I_ESP_H2xy2z_Fx2z_a+2*1*I_ESP_F2xy_Fx2z;
    abcd[iGrid*1260+1164] = 4.0E0*I_ESP_K2x5z_Fx2z_aa-2.0E0*3*I_ESP_H2x3z_Fx2z_a-2.0E0*4*I_ESP_H2x3z_Fx2z_a+3*2*I_ESP_F2xz_Fx2z;
    abcd[iGrid*1260+1165] = 4.0E0*I_ESP_Kx4y2z_Fx2z_aa-2.0E0*1*I_ESP_Hx4y_Fx2z_a;
    abcd[iGrid*1260+1166] = 4.0E0*I_ESP_Kx3y3z_Fx2z_aa-2.0E0*1*I_ESP_Hx3yz_Fx2z_a-2.0E0*2*I_ESP_Hx3yz_Fx2z_a;
    abcd[iGrid*1260+1167] = 4.0E0*I_ESP_Kx2y4z_Fx2z_aa-2.0E0*2*I_ESP_Hx2y2z_Fx2z_a-2.0E0*3*I_ESP_Hx2y2z_Fx2z_a+2*1*I_ESP_Fx2y_Fx2z;
    abcd[iGrid*1260+1168] = 4.0E0*I_ESP_Kxy5z_Fx2z_aa-2.0E0*3*I_ESP_Hxy3z_Fx2z_a-2.0E0*4*I_ESP_Hxy3z_Fx2z_a+3*2*I_ESP_Fxyz_Fx2z;
    abcd[iGrid*1260+1169] = 4.0E0*I_ESP_Kx6z_Fx2z_aa-2.0E0*4*I_ESP_Hx4z_Fx2z_a-2.0E0*5*I_ESP_Hx4z_Fx2z_a+4*3*I_ESP_Fx2z_Fx2z;
    abcd[iGrid*1260+1170] = 4.0E0*I_ESP_K5y2z_Fx2z_aa-2.0E0*1*I_ESP_H5y_Fx2z_a;
    abcd[iGrid*1260+1171] = 4.0E0*I_ESP_K4y3z_Fx2z_aa-2.0E0*1*I_ESP_H4yz_Fx2z_a-2.0E0*2*I_ESP_H4yz_Fx2z_a;
    abcd[iGrid*1260+1172] = 4.0E0*I_ESP_K3y4z_Fx2z_aa-2.0E0*2*I_ESP_H3y2z_Fx2z_a-2.0E0*3*I_ESP_H3y2z_Fx2z_a+2*1*I_ESP_F3y_Fx2z;
    abcd[iGrid*1260+1173] = 4.0E0*I_ESP_K2y5z_Fx2z_aa-2.0E0*3*I_ESP_H2y3z_Fx2z_a-2.0E0*4*I_ESP_H2y3z_Fx2z_a+3*2*I_ESP_F2yz_Fx2z;
    abcd[iGrid*1260+1174] = 4.0E0*I_ESP_Ky6z_Fx2z_aa-2.0E0*4*I_ESP_Hy4z_Fx2z_a-2.0E0*5*I_ESP_Hy4z_Fx2z_a+4*3*I_ESP_Fy2z_Fx2z;
    abcd[iGrid*1260+1175] = 4.0E0*I_ESP_K7z_Fx2z_aa-2.0E0*5*I_ESP_H5z_Fx2z_a-2.0E0*6*I_ESP_H5z_Fx2z_a+5*4*I_ESP_F3z_Fx2z;
    abcd[iGrid*1260+1176] = 4.0E0*I_ESP_K5x2z_F3y_aa-2.0E0*1*I_ESP_H5x_F3y_a;
    abcd[iGrid*1260+1177] = 4.0E0*I_ESP_K4xy2z_F3y_aa-2.0E0*1*I_ESP_H4xy_F3y_a;
    abcd[iGrid*1260+1178] = 4.0E0*I_ESP_K4x3z_F3y_aa-2.0E0*1*I_ESP_H4xz_F3y_a-2.0E0*2*I_ESP_H4xz_F3y_a;
    abcd[iGrid*1260+1179] = 4.0E0*I_ESP_K3x2y2z_F3y_aa-2.0E0*1*I_ESP_H3x2y_F3y_a;
    abcd[iGrid*1260+1180] = 4.0E0*I_ESP_K3xy3z_F3y_aa-2.0E0*1*I_ESP_H3xyz_F3y_a-2.0E0*2*I_ESP_H3xyz_F3y_a;
    abcd[iGrid*1260+1181] = 4.0E0*I_ESP_K3x4z_F3y_aa-2.0E0*2*I_ESP_H3x2z_F3y_a-2.0E0*3*I_ESP_H3x2z_F3y_a+2*1*I_ESP_F3x_F3y;
    abcd[iGrid*1260+1182] = 4.0E0*I_ESP_K2x3y2z_F3y_aa-2.0E0*1*I_ESP_H2x3y_F3y_a;
    abcd[iGrid*1260+1183] = 4.0E0*I_ESP_K2x2y3z_F3y_aa-2.0E0*1*I_ESP_H2x2yz_F3y_a-2.0E0*2*I_ESP_H2x2yz_F3y_a;
    abcd[iGrid*1260+1184] = 4.0E0*I_ESP_K2xy4z_F3y_aa-2.0E0*2*I_ESP_H2xy2z_F3y_a-2.0E0*3*I_ESP_H2xy2z_F3y_a+2*1*I_ESP_F2xy_F3y;
    abcd[iGrid*1260+1185] = 4.0E0*I_ESP_K2x5z_F3y_aa-2.0E0*3*I_ESP_H2x3z_F3y_a-2.0E0*4*I_ESP_H2x3z_F3y_a+3*2*I_ESP_F2xz_F3y;
    abcd[iGrid*1260+1186] = 4.0E0*I_ESP_Kx4y2z_F3y_aa-2.0E0*1*I_ESP_Hx4y_F3y_a;
    abcd[iGrid*1260+1187] = 4.0E0*I_ESP_Kx3y3z_F3y_aa-2.0E0*1*I_ESP_Hx3yz_F3y_a-2.0E0*2*I_ESP_Hx3yz_F3y_a;
    abcd[iGrid*1260+1188] = 4.0E0*I_ESP_Kx2y4z_F3y_aa-2.0E0*2*I_ESP_Hx2y2z_F3y_a-2.0E0*3*I_ESP_Hx2y2z_F3y_a+2*1*I_ESP_Fx2y_F3y;
    abcd[iGrid*1260+1189] = 4.0E0*I_ESP_Kxy5z_F3y_aa-2.0E0*3*I_ESP_Hxy3z_F3y_a-2.0E0*4*I_ESP_Hxy3z_F3y_a+3*2*I_ESP_Fxyz_F3y;
    abcd[iGrid*1260+1190] = 4.0E0*I_ESP_Kx6z_F3y_aa-2.0E0*4*I_ESP_Hx4z_F3y_a-2.0E0*5*I_ESP_Hx4z_F3y_a+4*3*I_ESP_Fx2z_F3y;
    abcd[iGrid*1260+1191] = 4.0E0*I_ESP_K5y2z_F3y_aa-2.0E0*1*I_ESP_H5y_F3y_a;
    abcd[iGrid*1260+1192] = 4.0E0*I_ESP_K4y3z_F3y_aa-2.0E0*1*I_ESP_H4yz_F3y_a-2.0E0*2*I_ESP_H4yz_F3y_a;
    abcd[iGrid*1260+1193] = 4.0E0*I_ESP_K3y4z_F3y_aa-2.0E0*2*I_ESP_H3y2z_F3y_a-2.0E0*3*I_ESP_H3y2z_F3y_a+2*1*I_ESP_F3y_F3y;
    abcd[iGrid*1260+1194] = 4.0E0*I_ESP_K2y5z_F3y_aa-2.0E0*3*I_ESP_H2y3z_F3y_a-2.0E0*4*I_ESP_H2y3z_F3y_a+3*2*I_ESP_F2yz_F3y;
    abcd[iGrid*1260+1195] = 4.0E0*I_ESP_Ky6z_F3y_aa-2.0E0*4*I_ESP_Hy4z_F3y_a-2.0E0*5*I_ESP_Hy4z_F3y_a+4*3*I_ESP_Fy2z_F3y;
    abcd[iGrid*1260+1196] = 4.0E0*I_ESP_K7z_F3y_aa-2.0E0*5*I_ESP_H5z_F3y_a-2.0E0*6*I_ESP_H5z_F3y_a+5*4*I_ESP_F3z_F3y;
    abcd[iGrid*1260+1197] = 4.0E0*I_ESP_K5x2z_F2yz_aa-2.0E0*1*I_ESP_H5x_F2yz_a;
    abcd[iGrid*1260+1198] = 4.0E0*I_ESP_K4xy2z_F2yz_aa-2.0E0*1*I_ESP_H4xy_F2yz_a;
    abcd[iGrid*1260+1199] = 4.0E0*I_ESP_K4x3z_F2yz_aa-2.0E0*1*I_ESP_H4xz_F2yz_a-2.0E0*2*I_ESP_H4xz_F2yz_a;
    abcd[iGrid*1260+1200] = 4.0E0*I_ESP_K3x2y2z_F2yz_aa-2.0E0*1*I_ESP_H3x2y_F2yz_a;
    abcd[iGrid*1260+1201] = 4.0E0*I_ESP_K3xy3z_F2yz_aa-2.0E0*1*I_ESP_H3xyz_F2yz_a-2.0E0*2*I_ESP_H3xyz_F2yz_a;
    abcd[iGrid*1260+1202] = 4.0E0*I_ESP_K3x4z_F2yz_aa-2.0E0*2*I_ESP_H3x2z_F2yz_a-2.0E0*3*I_ESP_H3x2z_F2yz_a+2*1*I_ESP_F3x_F2yz;
    abcd[iGrid*1260+1203] = 4.0E0*I_ESP_K2x3y2z_F2yz_aa-2.0E0*1*I_ESP_H2x3y_F2yz_a;
    abcd[iGrid*1260+1204] = 4.0E0*I_ESP_K2x2y3z_F2yz_aa-2.0E0*1*I_ESP_H2x2yz_F2yz_a-2.0E0*2*I_ESP_H2x2yz_F2yz_a;
    abcd[iGrid*1260+1205] = 4.0E0*I_ESP_K2xy4z_F2yz_aa-2.0E0*2*I_ESP_H2xy2z_F2yz_a-2.0E0*3*I_ESP_H2xy2z_F2yz_a+2*1*I_ESP_F2xy_F2yz;
    abcd[iGrid*1260+1206] = 4.0E0*I_ESP_K2x5z_F2yz_aa-2.0E0*3*I_ESP_H2x3z_F2yz_a-2.0E0*4*I_ESP_H2x3z_F2yz_a+3*2*I_ESP_F2xz_F2yz;
    abcd[iGrid*1260+1207] = 4.0E0*I_ESP_Kx4y2z_F2yz_aa-2.0E0*1*I_ESP_Hx4y_F2yz_a;
    abcd[iGrid*1260+1208] = 4.0E0*I_ESP_Kx3y3z_F2yz_aa-2.0E0*1*I_ESP_Hx3yz_F2yz_a-2.0E0*2*I_ESP_Hx3yz_F2yz_a;
    abcd[iGrid*1260+1209] = 4.0E0*I_ESP_Kx2y4z_F2yz_aa-2.0E0*2*I_ESP_Hx2y2z_F2yz_a-2.0E0*3*I_ESP_Hx2y2z_F2yz_a+2*1*I_ESP_Fx2y_F2yz;
    abcd[iGrid*1260+1210] = 4.0E0*I_ESP_Kxy5z_F2yz_aa-2.0E0*3*I_ESP_Hxy3z_F2yz_a-2.0E0*4*I_ESP_Hxy3z_F2yz_a+3*2*I_ESP_Fxyz_F2yz;
    abcd[iGrid*1260+1211] = 4.0E0*I_ESP_Kx6z_F2yz_aa-2.0E0*4*I_ESP_Hx4z_F2yz_a-2.0E0*5*I_ESP_Hx4z_F2yz_a+4*3*I_ESP_Fx2z_F2yz;
    abcd[iGrid*1260+1212] = 4.0E0*I_ESP_K5y2z_F2yz_aa-2.0E0*1*I_ESP_H5y_F2yz_a;
    abcd[iGrid*1260+1213] = 4.0E0*I_ESP_K4y3z_F2yz_aa-2.0E0*1*I_ESP_H4yz_F2yz_a-2.0E0*2*I_ESP_H4yz_F2yz_a;
    abcd[iGrid*1260+1214] = 4.0E0*I_ESP_K3y4z_F2yz_aa-2.0E0*2*I_ESP_H3y2z_F2yz_a-2.0E0*3*I_ESP_H3y2z_F2yz_a+2*1*I_ESP_F3y_F2yz;
    abcd[iGrid*1260+1215] = 4.0E0*I_ESP_K2y5z_F2yz_aa-2.0E0*3*I_ESP_H2y3z_F2yz_a-2.0E0*4*I_ESP_H2y3z_F2yz_a+3*2*I_ESP_F2yz_F2yz;
    abcd[iGrid*1260+1216] = 4.0E0*I_ESP_Ky6z_F2yz_aa-2.0E0*4*I_ESP_Hy4z_F2yz_a-2.0E0*5*I_ESP_Hy4z_F2yz_a+4*3*I_ESP_Fy2z_F2yz;
    abcd[iGrid*1260+1217] = 4.0E0*I_ESP_K7z_F2yz_aa-2.0E0*5*I_ESP_H5z_F2yz_a-2.0E0*6*I_ESP_H5z_F2yz_a+5*4*I_ESP_F3z_F2yz;
    abcd[iGrid*1260+1218] = 4.0E0*I_ESP_K5x2z_Fy2z_aa-2.0E0*1*I_ESP_H5x_Fy2z_a;
    abcd[iGrid*1260+1219] = 4.0E0*I_ESP_K4xy2z_Fy2z_aa-2.0E0*1*I_ESP_H4xy_Fy2z_a;
    abcd[iGrid*1260+1220] = 4.0E0*I_ESP_K4x3z_Fy2z_aa-2.0E0*1*I_ESP_H4xz_Fy2z_a-2.0E0*2*I_ESP_H4xz_Fy2z_a;
    abcd[iGrid*1260+1221] = 4.0E0*I_ESP_K3x2y2z_Fy2z_aa-2.0E0*1*I_ESP_H3x2y_Fy2z_a;
    abcd[iGrid*1260+1222] = 4.0E0*I_ESP_K3xy3z_Fy2z_aa-2.0E0*1*I_ESP_H3xyz_Fy2z_a-2.0E0*2*I_ESP_H3xyz_Fy2z_a;
    abcd[iGrid*1260+1223] = 4.0E0*I_ESP_K3x4z_Fy2z_aa-2.0E0*2*I_ESP_H3x2z_Fy2z_a-2.0E0*3*I_ESP_H3x2z_Fy2z_a+2*1*I_ESP_F3x_Fy2z;
    abcd[iGrid*1260+1224] = 4.0E0*I_ESP_K2x3y2z_Fy2z_aa-2.0E0*1*I_ESP_H2x3y_Fy2z_a;
    abcd[iGrid*1260+1225] = 4.0E0*I_ESP_K2x2y3z_Fy2z_aa-2.0E0*1*I_ESP_H2x2yz_Fy2z_a-2.0E0*2*I_ESP_H2x2yz_Fy2z_a;
    abcd[iGrid*1260+1226] = 4.0E0*I_ESP_K2xy4z_Fy2z_aa-2.0E0*2*I_ESP_H2xy2z_Fy2z_a-2.0E0*3*I_ESP_H2xy2z_Fy2z_a+2*1*I_ESP_F2xy_Fy2z;
    abcd[iGrid*1260+1227] = 4.0E0*I_ESP_K2x5z_Fy2z_aa-2.0E0*3*I_ESP_H2x3z_Fy2z_a-2.0E0*4*I_ESP_H2x3z_Fy2z_a+3*2*I_ESP_F2xz_Fy2z;
    abcd[iGrid*1260+1228] = 4.0E0*I_ESP_Kx4y2z_Fy2z_aa-2.0E0*1*I_ESP_Hx4y_Fy2z_a;
    abcd[iGrid*1260+1229] = 4.0E0*I_ESP_Kx3y3z_Fy2z_aa-2.0E0*1*I_ESP_Hx3yz_Fy2z_a-2.0E0*2*I_ESP_Hx3yz_Fy2z_a;
    abcd[iGrid*1260+1230] = 4.0E0*I_ESP_Kx2y4z_Fy2z_aa-2.0E0*2*I_ESP_Hx2y2z_Fy2z_a-2.0E0*3*I_ESP_Hx2y2z_Fy2z_a+2*1*I_ESP_Fx2y_Fy2z;
    abcd[iGrid*1260+1231] = 4.0E0*I_ESP_Kxy5z_Fy2z_aa-2.0E0*3*I_ESP_Hxy3z_Fy2z_a-2.0E0*4*I_ESP_Hxy3z_Fy2z_a+3*2*I_ESP_Fxyz_Fy2z;
    abcd[iGrid*1260+1232] = 4.0E0*I_ESP_Kx6z_Fy2z_aa-2.0E0*4*I_ESP_Hx4z_Fy2z_a-2.0E0*5*I_ESP_Hx4z_Fy2z_a+4*3*I_ESP_Fx2z_Fy2z;
    abcd[iGrid*1260+1233] = 4.0E0*I_ESP_K5y2z_Fy2z_aa-2.0E0*1*I_ESP_H5y_Fy2z_a;
    abcd[iGrid*1260+1234] = 4.0E0*I_ESP_K4y3z_Fy2z_aa-2.0E0*1*I_ESP_H4yz_Fy2z_a-2.0E0*2*I_ESP_H4yz_Fy2z_a;
    abcd[iGrid*1260+1235] = 4.0E0*I_ESP_K3y4z_Fy2z_aa-2.0E0*2*I_ESP_H3y2z_Fy2z_a-2.0E0*3*I_ESP_H3y2z_Fy2z_a+2*1*I_ESP_F3y_Fy2z;
    abcd[iGrid*1260+1236] = 4.0E0*I_ESP_K2y5z_Fy2z_aa-2.0E0*3*I_ESP_H2y3z_Fy2z_a-2.0E0*4*I_ESP_H2y3z_Fy2z_a+3*2*I_ESP_F2yz_Fy2z;
    abcd[iGrid*1260+1237] = 4.0E0*I_ESP_Ky6z_Fy2z_aa-2.0E0*4*I_ESP_Hy4z_Fy2z_a-2.0E0*5*I_ESP_Hy4z_Fy2z_a+4*3*I_ESP_Fy2z_Fy2z;
    abcd[iGrid*1260+1238] = 4.0E0*I_ESP_K7z_Fy2z_aa-2.0E0*5*I_ESP_H5z_Fy2z_a-2.0E0*6*I_ESP_H5z_Fy2z_a+5*4*I_ESP_F3z_Fy2z;
    abcd[iGrid*1260+1239] = 4.0E0*I_ESP_K5x2z_F3z_aa-2.0E0*1*I_ESP_H5x_F3z_a;
    abcd[iGrid*1260+1240] = 4.0E0*I_ESP_K4xy2z_F3z_aa-2.0E0*1*I_ESP_H4xy_F3z_a;
    abcd[iGrid*1260+1241] = 4.0E0*I_ESP_K4x3z_F3z_aa-2.0E0*1*I_ESP_H4xz_F3z_a-2.0E0*2*I_ESP_H4xz_F3z_a;
    abcd[iGrid*1260+1242] = 4.0E0*I_ESP_K3x2y2z_F3z_aa-2.0E0*1*I_ESP_H3x2y_F3z_a;
    abcd[iGrid*1260+1243] = 4.0E0*I_ESP_K3xy3z_F3z_aa-2.0E0*1*I_ESP_H3xyz_F3z_a-2.0E0*2*I_ESP_H3xyz_F3z_a;
    abcd[iGrid*1260+1244] = 4.0E0*I_ESP_K3x4z_F3z_aa-2.0E0*2*I_ESP_H3x2z_F3z_a-2.0E0*3*I_ESP_H3x2z_F3z_a+2*1*I_ESP_F3x_F3z;
    abcd[iGrid*1260+1245] = 4.0E0*I_ESP_K2x3y2z_F3z_aa-2.0E0*1*I_ESP_H2x3y_F3z_a;
    abcd[iGrid*1260+1246] = 4.0E0*I_ESP_K2x2y3z_F3z_aa-2.0E0*1*I_ESP_H2x2yz_F3z_a-2.0E0*2*I_ESP_H2x2yz_F3z_a;
    abcd[iGrid*1260+1247] = 4.0E0*I_ESP_K2xy4z_F3z_aa-2.0E0*2*I_ESP_H2xy2z_F3z_a-2.0E0*3*I_ESP_H2xy2z_F3z_a+2*1*I_ESP_F2xy_F3z;
    abcd[iGrid*1260+1248] = 4.0E0*I_ESP_K2x5z_F3z_aa-2.0E0*3*I_ESP_H2x3z_F3z_a-2.0E0*4*I_ESP_H2x3z_F3z_a+3*2*I_ESP_F2xz_F3z;
    abcd[iGrid*1260+1249] = 4.0E0*I_ESP_Kx4y2z_F3z_aa-2.0E0*1*I_ESP_Hx4y_F3z_a;
    abcd[iGrid*1260+1250] = 4.0E0*I_ESP_Kx3y3z_F3z_aa-2.0E0*1*I_ESP_Hx3yz_F3z_a-2.0E0*2*I_ESP_Hx3yz_F3z_a;
    abcd[iGrid*1260+1251] = 4.0E0*I_ESP_Kx2y4z_F3z_aa-2.0E0*2*I_ESP_Hx2y2z_F3z_a-2.0E0*3*I_ESP_Hx2y2z_F3z_a+2*1*I_ESP_Fx2y_F3z;
    abcd[iGrid*1260+1252] = 4.0E0*I_ESP_Kxy5z_F3z_aa-2.0E0*3*I_ESP_Hxy3z_F3z_a-2.0E0*4*I_ESP_Hxy3z_F3z_a+3*2*I_ESP_Fxyz_F3z;
    abcd[iGrid*1260+1253] = 4.0E0*I_ESP_Kx6z_F3z_aa-2.0E0*4*I_ESP_Hx4z_F3z_a-2.0E0*5*I_ESP_Hx4z_F3z_a+4*3*I_ESP_Fx2z_F3z;
    abcd[iGrid*1260+1254] = 4.0E0*I_ESP_K5y2z_F3z_aa-2.0E0*1*I_ESP_H5y_F3z_a;
    abcd[iGrid*1260+1255] = 4.0E0*I_ESP_K4y3z_F3z_aa-2.0E0*1*I_ESP_H4yz_F3z_a-2.0E0*2*I_ESP_H4yz_F3z_a;
    abcd[iGrid*1260+1256] = 4.0E0*I_ESP_K3y4z_F3z_aa-2.0E0*2*I_ESP_H3y2z_F3z_a-2.0E0*3*I_ESP_H3y2z_F3z_a+2*1*I_ESP_F3y_F3z;
    abcd[iGrid*1260+1257] = 4.0E0*I_ESP_K2y5z_F3z_aa-2.0E0*3*I_ESP_H2y3z_F3z_a-2.0E0*4*I_ESP_H2y3z_F3z_a+3*2*I_ESP_F2yz_F3z;
    abcd[iGrid*1260+1258] = 4.0E0*I_ESP_Ky6z_F3z_aa-2.0E0*4*I_ESP_Hy4z_F3z_a-2.0E0*5*I_ESP_Hy4z_F3z_a+4*3*I_ESP_Fy2z_F3z;
    abcd[iGrid*1260+1259] = 4.0E0*I_ESP_K7z_F3z_aa-2.0E0*5*I_ESP_H5z_F3z_a-2.0E0*6*I_ESP_H5z_F3z_a+5*4*I_ESP_F3z_F3z;
  }
}
