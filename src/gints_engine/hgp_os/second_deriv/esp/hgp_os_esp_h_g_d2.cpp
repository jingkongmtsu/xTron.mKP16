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
#include <boost/math/special_functions/gamma.hpp>

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
// BRA1 as redundant position, total RHS integrals evaluated as: 18999
// BRA2 as redundant position, total RHS integrals evaluated as: 16858
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

void hgp_os_esp_h_g_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_O11x_S_aa = 0.0E0;
    Double I_ESP_O10xy_S_aa = 0.0E0;
    Double I_ESP_O10xz_S_aa = 0.0E0;
    Double I_ESP_O9x2y_S_aa = 0.0E0;
    Double I_ESP_O9xyz_S_aa = 0.0E0;
    Double I_ESP_O9x2z_S_aa = 0.0E0;
    Double I_ESP_O8x3y_S_aa = 0.0E0;
    Double I_ESP_O8x2yz_S_aa = 0.0E0;
    Double I_ESP_O8xy2z_S_aa = 0.0E0;
    Double I_ESP_O8x3z_S_aa = 0.0E0;
    Double I_ESP_O7x4y_S_aa = 0.0E0;
    Double I_ESP_O7x3yz_S_aa = 0.0E0;
    Double I_ESP_O7x2y2z_S_aa = 0.0E0;
    Double I_ESP_O7xy3z_S_aa = 0.0E0;
    Double I_ESP_O7x4z_S_aa = 0.0E0;
    Double I_ESP_O6x5y_S_aa = 0.0E0;
    Double I_ESP_O6x4yz_S_aa = 0.0E0;
    Double I_ESP_O6x3y2z_S_aa = 0.0E0;
    Double I_ESP_O6x2y3z_S_aa = 0.0E0;
    Double I_ESP_O6xy4z_S_aa = 0.0E0;
    Double I_ESP_O6x5z_S_aa = 0.0E0;
    Double I_ESP_O5x6y_S_aa = 0.0E0;
    Double I_ESP_O5x5yz_S_aa = 0.0E0;
    Double I_ESP_O5x4y2z_S_aa = 0.0E0;
    Double I_ESP_O5x3y3z_S_aa = 0.0E0;
    Double I_ESP_O5x2y4z_S_aa = 0.0E0;
    Double I_ESP_O5xy5z_S_aa = 0.0E0;
    Double I_ESP_O5x6z_S_aa = 0.0E0;
    Double I_ESP_O4x7y_S_aa = 0.0E0;
    Double I_ESP_O4x6yz_S_aa = 0.0E0;
    Double I_ESP_O4x5y2z_S_aa = 0.0E0;
    Double I_ESP_O4x4y3z_S_aa = 0.0E0;
    Double I_ESP_O4x3y4z_S_aa = 0.0E0;
    Double I_ESP_O4x2y5z_S_aa = 0.0E0;
    Double I_ESP_O4xy6z_S_aa = 0.0E0;
    Double I_ESP_O4x7z_S_aa = 0.0E0;
    Double I_ESP_O3x8y_S_aa = 0.0E0;
    Double I_ESP_O3x7yz_S_aa = 0.0E0;
    Double I_ESP_O3x6y2z_S_aa = 0.0E0;
    Double I_ESP_O3x5y3z_S_aa = 0.0E0;
    Double I_ESP_O3x4y4z_S_aa = 0.0E0;
    Double I_ESP_O3x3y5z_S_aa = 0.0E0;
    Double I_ESP_O3x2y6z_S_aa = 0.0E0;
    Double I_ESP_O3xy7z_S_aa = 0.0E0;
    Double I_ESP_O3x8z_S_aa = 0.0E0;
    Double I_ESP_O2x9y_S_aa = 0.0E0;
    Double I_ESP_O2x8yz_S_aa = 0.0E0;
    Double I_ESP_O2x7y2z_S_aa = 0.0E0;
    Double I_ESP_O2x6y3z_S_aa = 0.0E0;
    Double I_ESP_O2x5y4z_S_aa = 0.0E0;
    Double I_ESP_O2x4y5z_S_aa = 0.0E0;
    Double I_ESP_O2x3y6z_S_aa = 0.0E0;
    Double I_ESP_O2x2y7z_S_aa = 0.0E0;
    Double I_ESP_O2xy8z_S_aa = 0.0E0;
    Double I_ESP_O2x9z_S_aa = 0.0E0;
    Double I_ESP_Ox10y_S_aa = 0.0E0;
    Double I_ESP_Ox9yz_S_aa = 0.0E0;
    Double I_ESP_Ox8y2z_S_aa = 0.0E0;
    Double I_ESP_Ox7y3z_S_aa = 0.0E0;
    Double I_ESP_Ox6y4z_S_aa = 0.0E0;
    Double I_ESP_Ox5y5z_S_aa = 0.0E0;
    Double I_ESP_Ox4y6z_S_aa = 0.0E0;
    Double I_ESP_Ox3y7z_S_aa = 0.0E0;
    Double I_ESP_Ox2y8z_S_aa = 0.0E0;
    Double I_ESP_Oxy9z_S_aa = 0.0E0;
    Double I_ESP_Ox10z_S_aa = 0.0E0;
    Double I_ESP_O11y_S_aa = 0.0E0;
    Double I_ESP_O10yz_S_aa = 0.0E0;
    Double I_ESP_O9y2z_S_aa = 0.0E0;
    Double I_ESP_O8y3z_S_aa = 0.0E0;
    Double I_ESP_O7y4z_S_aa = 0.0E0;
    Double I_ESP_O6y5z_S_aa = 0.0E0;
    Double I_ESP_O5y6z_S_aa = 0.0E0;
    Double I_ESP_O4y7z_S_aa = 0.0E0;
    Double I_ESP_O3y8z_S_aa = 0.0E0;
    Double I_ESP_O2y9z_S_aa = 0.0E0;
    Double I_ESP_Oy10z_S_aa = 0.0E0;
    Double I_ESP_O11z_S_aa = 0.0E0;
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
      Double I_ESP_S_S_M10_vrr  = 0.0E0;
      Double I_ESP_S_S_M11_vrr  = 0.0E0;

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER57;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER55*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER53*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER51*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER49*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER47*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER45*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER43*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER41*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER39*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER37*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER35*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER33*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER31*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER29*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER27*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = 1.0E0+u2*ONEOVER25*I_ESP_S_S_M11_vrr;
        I_ESP_S_S_M11_vrr = ONEOVER23*I_ESP_S_S_M11_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M11_vrr  = f*I_ESP_S_S_M11_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ESP_S_S_M10_vrr  = ONEOVER21*(u2*I_ESP_S_S_M11_vrr+f);
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

        //
        // now here for maxM>M_limit
        // use external function to calculate f_{Mmax}(t)
        // then use down recursive relation to get others
        //

        // calculate (SS|SS)^{Mmax} with incomplete gamma function
        // currently we use boost library for calculation
        if (fabs(u)<THRESHOLD_MATH) {
          I_ESP_S_S_M11_vrr = 1.0E0/(2.0E0*11+1.0E0);
        }else{
          I_ESP_S_S_M11_vrr = (0.5E0/squ)*boost::math::tgamma_lower(11+0.5E0,u);
          Double oneOveru = 1.0E0/u;
          for(UInt i=0; i<11; i++) {
            I_ESP_S_S_M11_vrr = I_ESP_S_S_M11_vrr*oneOveru;
          }
        }
        Double f = TWOOVERSQRTPI*prefactor*sqrho;
        I_ESP_S_S_M11_vrr *= f;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        Double u2 = 2.0E0*u;
        Double eu = exp(-u);
        f = f*eu;
        I_ESP_S_S_M10_vrr  = ONEOVER21*(u2*I_ESP_S_S_M11_vrr+f);
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
      }


        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M10
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M10
         * RHS shell quartet name: SQ_ESP_S_S_M11
         ************************************************************/
        Double I_ESP_Px_S_M10_vrr = PAX*I_ESP_S_S_M10_vrr-PRX*I_ESP_S_S_M11_vrr;
        Double I_ESP_Py_S_M10_vrr = PAY*I_ESP_S_S_M10_vrr-PRY*I_ESP_S_S_M11_vrr;
        Double I_ESP_Pz_S_M10_vrr = PAZ*I_ESP_S_S_M10_vrr-PRZ*I_ESP_S_S_M11_vrr;

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
         * shell quartet name: SQ_ESP_D_S_M9
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M9
         * RHS shell quartet name: SQ_ESP_P_S_M10
         * RHS shell quartet name: SQ_ESP_S_S_M9
         * RHS shell quartet name: SQ_ESP_S_S_M10
         ************************************************************/
        Double I_ESP_D2x_S_M9_vrr = PAX*I_ESP_Px_S_M9_vrr-PRX*I_ESP_Px_S_M10_vrr+oned2z*I_ESP_S_S_M9_vrr-oned2z*I_ESP_S_S_M10_vrr;
        Double I_ESP_D2y_S_M9_vrr = PAY*I_ESP_Py_S_M9_vrr-PRY*I_ESP_Py_S_M10_vrr+oned2z*I_ESP_S_S_M9_vrr-oned2z*I_ESP_S_S_M10_vrr;
        Double I_ESP_D2z_S_M9_vrr = PAZ*I_ESP_Pz_S_M9_vrr-PRZ*I_ESP_Pz_S_M10_vrr+oned2z*I_ESP_S_S_M9_vrr-oned2z*I_ESP_S_S_M10_vrr;

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
         * shell quartet name: SQ_ESP_F_S_M8
         * expanding position: BRA1
         * code section is: VRR
         * totally 7 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M8
         * RHS shell quartet name: SQ_ESP_D_S_M9
         * RHS shell quartet name: SQ_ESP_P_S_M8
         * RHS shell quartet name: SQ_ESP_P_S_M9
         ************************************************************/
        Double I_ESP_F3x_S_M8_vrr = PAX*I_ESP_D2x_S_M8_vrr-PRX*I_ESP_D2x_S_M9_vrr+2*oned2z*I_ESP_Px_S_M8_vrr-2*oned2z*I_ESP_Px_S_M9_vrr;
        Double I_ESP_F3y_S_M8_vrr = PAY*I_ESP_D2y_S_M8_vrr-PRY*I_ESP_D2y_S_M9_vrr+2*oned2z*I_ESP_Py_S_M8_vrr-2*oned2z*I_ESP_Py_S_M9_vrr;
        Double I_ESP_F3z_S_M8_vrr = PAZ*I_ESP_D2z_S_M8_vrr-PRZ*I_ESP_D2z_S_M9_vrr+2*oned2z*I_ESP_Pz_S_M8_vrr-2*oned2z*I_ESP_Pz_S_M9_vrr;

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
         * shell quartet name: SQ_ESP_G_S_M7
         * expanding position: BRA1
         * code section is: VRR
         * totally 12 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M7
         * RHS shell quartet name: SQ_ESP_F_S_M8
         * RHS shell quartet name: SQ_ESP_D_S_M7
         * RHS shell quartet name: SQ_ESP_D_S_M8
         ************************************************************/
        Double I_ESP_G4x_S_M7_vrr = PAX*I_ESP_F3x_S_M7_vrr-PRX*I_ESP_F3x_S_M8_vrr+3*oned2z*I_ESP_D2x_S_M7_vrr-3*oned2z*I_ESP_D2x_S_M8_vrr;
        Double I_ESP_G4y_S_M7_vrr = PAY*I_ESP_F3y_S_M7_vrr-PRY*I_ESP_F3y_S_M8_vrr+3*oned2z*I_ESP_D2y_S_M7_vrr-3*oned2z*I_ESP_D2y_S_M8_vrr;
        Double I_ESP_G4z_S_M7_vrr = PAZ*I_ESP_F3z_S_M7_vrr-PRZ*I_ESP_F3z_S_M8_vrr+3*oned2z*I_ESP_D2z_S_M7_vrr-3*oned2z*I_ESP_D2z_S_M8_vrr;

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
         * shell quartet name: SQ_ESP_H_S_M6
         * expanding position: BRA1
         * code section is: VRR
         * totally 15 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S_M6
         * RHS shell quartet name: SQ_ESP_G_S_M7
         * RHS shell quartet name: SQ_ESP_F_S_M6
         * RHS shell quartet name: SQ_ESP_F_S_M7
         ************************************************************/
        Double I_ESP_H5x_S_M6_vrr = PAX*I_ESP_G4x_S_M6_vrr-PRX*I_ESP_G4x_S_M7_vrr+4*oned2z*I_ESP_F3x_S_M6_vrr-4*oned2z*I_ESP_F3x_S_M7_vrr;
        Double I_ESP_H4xy_S_M6_vrr = PAY*I_ESP_G4x_S_M6_vrr-PRY*I_ESP_G4x_S_M7_vrr;
        Double I_ESP_H4xz_S_M6_vrr = PAZ*I_ESP_G4x_S_M6_vrr-PRZ*I_ESP_G4x_S_M7_vrr;
        Double I_ESP_H5y_S_M6_vrr = PAY*I_ESP_G4y_S_M6_vrr-PRY*I_ESP_G4y_S_M7_vrr+4*oned2z*I_ESP_F3y_S_M6_vrr-4*oned2z*I_ESP_F3y_S_M7_vrr;
        Double I_ESP_H4yz_S_M6_vrr = PAZ*I_ESP_G4y_S_M6_vrr-PRZ*I_ESP_G4y_S_M7_vrr;
        Double I_ESP_H5z_S_M6_vrr = PAZ*I_ESP_G4z_S_M6_vrr-PRZ*I_ESP_G4z_S_M7_vrr+4*oned2z*I_ESP_F3z_S_M6_vrr-4*oned2z*I_ESP_F3z_S_M7_vrr;

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
         * shell quartet name: SQ_ESP_I_S_M5
         * expanding position: BRA1
         * code section is: VRR
         * totally 16 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_H_S_M5
         * RHS shell quartet name: SQ_ESP_H_S_M6
         * RHS shell quartet name: SQ_ESP_G_S_M5
         * RHS shell quartet name: SQ_ESP_G_S_M6
         ************************************************************/
        Double I_ESP_I6x_S_M5_vrr = PAX*I_ESP_H5x_S_M5_vrr-PRX*I_ESP_H5x_S_M6_vrr+5*oned2z*I_ESP_G4x_S_M5_vrr-5*oned2z*I_ESP_G4x_S_M6_vrr;
        Double I_ESP_I5xy_S_M5_vrr = PAY*I_ESP_H5x_S_M5_vrr-PRY*I_ESP_H5x_S_M6_vrr;
        Double I_ESP_I5xz_S_M5_vrr = PAZ*I_ESP_H5x_S_M5_vrr-PRZ*I_ESP_H5x_S_M6_vrr;
        Double I_ESP_I4x2y_S_M5_vrr = PAY*I_ESP_H4xy_S_M5_vrr-PRY*I_ESP_H4xy_S_M6_vrr+oned2z*I_ESP_G4x_S_M5_vrr-oned2z*I_ESP_G4x_S_M6_vrr;
        Double I_ESP_I4x2z_S_M5_vrr = PAZ*I_ESP_H4xz_S_M5_vrr-PRZ*I_ESP_H4xz_S_M6_vrr+oned2z*I_ESP_G4x_S_M5_vrr-oned2z*I_ESP_G4x_S_M6_vrr;
        Double I_ESP_Ix5y_S_M5_vrr = PAX*I_ESP_H5y_S_M5_vrr-PRX*I_ESP_H5y_S_M6_vrr;
        Double I_ESP_Ix5z_S_M5_vrr = PAX*I_ESP_H5z_S_M5_vrr-PRX*I_ESP_H5z_S_M6_vrr;
        Double I_ESP_I6y_S_M5_vrr = PAY*I_ESP_H5y_S_M5_vrr-PRY*I_ESP_H5y_S_M6_vrr+5*oned2z*I_ESP_G4y_S_M5_vrr-5*oned2z*I_ESP_G4y_S_M6_vrr;
        Double I_ESP_I5yz_S_M5_vrr = PAZ*I_ESP_H5y_S_M5_vrr-PRZ*I_ESP_H5y_S_M6_vrr;
        Double I_ESP_I4y2z_S_M5_vrr = PAZ*I_ESP_H4yz_S_M5_vrr-PRZ*I_ESP_H4yz_S_M6_vrr+oned2z*I_ESP_G4y_S_M5_vrr-oned2z*I_ESP_G4y_S_M6_vrr;
        Double I_ESP_Iy5z_S_M5_vrr = PAY*I_ESP_H5z_S_M5_vrr-PRY*I_ESP_H5z_S_M6_vrr;
        Double I_ESP_I6z_S_M5_vrr = PAZ*I_ESP_H5z_S_M5_vrr-PRZ*I_ESP_H5z_S_M6_vrr+5*oned2z*I_ESP_G4z_S_M5_vrr-5*oned2z*I_ESP_G4z_S_M6_vrr;

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
         * shell quartet name: SQ_ESP_K_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 18 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_I_S_M4
         * RHS shell quartet name: SQ_ESP_I_S_M5
         * RHS shell quartet name: SQ_ESP_H_S_M4
         * RHS shell quartet name: SQ_ESP_H_S_M5
         ************************************************************/
        Double I_ESP_K7x_S_M4_vrr = PAX*I_ESP_I6x_S_M4_vrr-PRX*I_ESP_I6x_S_M5_vrr+6*oned2z*I_ESP_H5x_S_M4_vrr-6*oned2z*I_ESP_H5x_S_M5_vrr;
        Double I_ESP_K6xy_S_M4_vrr = PAY*I_ESP_I6x_S_M4_vrr-PRY*I_ESP_I6x_S_M5_vrr;
        Double I_ESP_K6xz_S_M4_vrr = PAZ*I_ESP_I6x_S_M4_vrr-PRZ*I_ESP_I6x_S_M5_vrr;
        Double I_ESP_K5x2y_S_M4_vrr = PAY*I_ESP_I5xy_S_M4_vrr-PRY*I_ESP_I5xy_S_M5_vrr+oned2z*I_ESP_H5x_S_M4_vrr-oned2z*I_ESP_H5x_S_M5_vrr;
        Double I_ESP_K5x2z_S_M4_vrr = PAZ*I_ESP_I5xz_S_M4_vrr-PRZ*I_ESP_I5xz_S_M5_vrr+oned2z*I_ESP_H5x_S_M4_vrr-oned2z*I_ESP_H5x_S_M5_vrr;
        Double I_ESP_K4x3y_S_M4_vrr = PAY*I_ESP_I4x2y_S_M4_vrr-PRY*I_ESP_I4x2y_S_M5_vrr+2*oned2z*I_ESP_H4xy_S_M4_vrr-2*oned2z*I_ESP_H4xy_S_M5_vrr;
        Double I_ESP_K4x3z_S_M4_vrr = PAZ*I_ESP_I4x2z_S_M4_vrr-PRZ*I_ESP_I4x2z_S_M5_vrr+2*oned2z*I_ESP_H4xz_S_M4_vrr-2*oned2z*I_ESP_H4xz_S_M5_vrr;
        Double I_ESP_K2x5y_S_M4_vrr = PAX*I_ESP_Ix5y_S_M4_vrr-PRX*I_ESP_Ix5y_S_M5_vrr+oned2z*I_ESP_H5y_S_M4_vrr-oned2z*I_ESP_H5y_S_M5_vrr;
        Double I_ESP_K2x5z_S_M4_vrr = PAX*I_ESP_Ix5z_S_M4_vrr-PRX*I_ESP_Ix5z_S_M5_vrr+oned2z*I_ESP_H5z_S_M4_vrr-oned2z*I_ESP_H5z_S_M5_vrr;
        Double I_ESP_Kx6y_S_M4_vrr = PAX*I_ESP_I6y_S_M4_vrr-PRX*I_ESP_I6y_S_M5_vrr;
        Double I_ESP_Kx6z_S_M4_vrr = PAX*I_ESP_I6z_S_M4_vrr-PRX*I_ESP_I6z_S_M5_vrr;
        Double I_ESP_K7y_S_M4_vrr = PAY*I_ESP_I6y_S_M4_vrr-PRY*I_ESP_I6y_S_M5_vrr+6*oned2z*I_ESP_H5y_S_M4_vrr-6*oned2z*I_ESP_H5y_S_M5_vrr;
        Double I_ESP_K6yz_S_M4_vrr = PAZ*I_ESP_I6y_S_M4_vrr-PRZ*I_ESP_I6y_S_M5_vrr;
        Double I_ESP_K5y2z_S_M4_vrr = PAZ*I_ESP_I5yz_S_M4_vrr-PRZ*I_ESP_I5yz_S_M5_vrr+oned2z*I_ESP_H5y_S_M4_vrr-oned2z*I_ESP_H5y_S_M5_vrr;
        Double I_ESP_K4y3z_S_M4_vrr = PAZ*I_ESP_I4y2z_S_M4_vrr-PRZ*I_ESP_I4y2z_S_M5_vrr+2*oned2z*I_ESP_H4yz_S_M4_vrr-2*oned2z*I_ESP_H4yz_S_M5_vrr;
        Double I_ESP_K2y5z_S_M4_vrr = PAY*I_ESP_Iy5z_S_M4_vrr-PRY*I_ESP_Iy5z_S_M5_vrr+oned2z*I_ESP_H5z_S_M4_vrr-oned2z*I_ESP_H5z_S_M5_vrr;
        Double I_ESP_Ky6z_S_M4_vrr = PAY*I_ESP_I6z_S_M4_vrr-PRY*I_ESP_I6z_S_M5_vrr;
        Double I_ESP_K7z_S_M4_vrr = PAZ*I_ESP_I6z_S_M4_vrr-PRZ*I_ESP_I6z_S_M5_vrr+6*oned2z*I_ESP_H5z_S_M4_vrr-6*oned2z*I_ESP_H5z_S_M5_vrr;

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
         * shell quartet name: SQ_ESP_L_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 21 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_K_S_M3
         * RHS shell quartet name: SQ_ESP_K_S_M4
         * RHS shell quartet name: SQ_ESP_I_S_M3
         * RHS shell quartet name: SQ_ESP_I_S_M4
         ************************************************************/
        Double I_ESP_L8x_S_M3_vrr = PAX*I_ESP_K7x_S_M3_vrr-PRX*I_ESP_K7x_S_M4_vrr+7*oned2z*I_ESP_I6x_S_M3_vrr-7*oned2z*I_ESP_I6x_S_M4_vrr;
        Double I_ESP_L7xy_S_M3_vrr = PAY*I_ESP_K7x_S_M3_vrr-PRY*I_ESP_K7x_S_M4_vrr;
        Double I_ESP_L7xz_S_M3_vrr = PAZ*I_ESP_K7x_S_M3_vrr-PRZ*I_ESP_K7x_S_M4_vrr;
        Double I_ESP_L6x2y_S_M3_vrr = PAY*I_ESP_K6xy_S_M3_vrr-PRY*I_ESP_K6xy_S_M4_vrr+oned2z*I_ESP_I6x_S_M3_vrr-oned2z*I_ESP_I6x_S_M4_vrr;
        Double I_ESP_L6x2z_S_M3_vrr = PAZ*I_ESP_K6xz_S_M3_vrr-PRZ*I_ESP_K6xz_S_M4_vrr+oned2z*I_ESP_I6x_S_M3_vrr-oned2z*I_ESP_I6x_S_M4_vrr;
        Double I_ESP_L5x3y_S_M3_vrr = PAY*I_ESP_K5x2y_S_M3_vrr-PRY*I_ESP_K5x2y_S_M4_vrr+2*oned2z*I_ESP_I5xy_S_M3_vrr-2*oned2z*I_ESP_I5xy_S_M4_vrr;
        Double I_ESP_L5x3z_S_M3_vrr = PAZ*I_ESP_K5x2z_S_M3_vrr-PRZ*I_ESP_K5x2z_S_M4_vrr+2*oned2z*I_ESP_I5xz_S_M3_vrr-2*oned2z*I_ESP_I5xz_S_M4_vrr;
        Double I_ESP_L4x4y_S_M3_vrr = PAY*I_ESP_K4x3y_S_M3_vrr-PRY*I_ESP_K4x3y_S_M4_vrr+3*oned2z*I_ESP_I4x2y_S_M3_vrr-3*oned2z*I_ESP_I4x2y_S_M4_vrr;
        Double I_ESP_L4x4z_S_M3_vrr = PAZ*I_ESP_K4x3z_S_M3_vrr-PRZ*I_ESP_K4x3z_S_M4_vrr+3*oned2z*I_ESP_I4x2z_S_M3_vrr-3*oned2z*I_ESP_I4x2z_S_M4_vrr;
        Double I_ESP_L3x5y_S_M3_vrr = PAX*I_ESP_K2x5y_S_M3_vrr-PRX*I_ESP_K2x5y_S_M4_vrr+2*oned2z*I_ESP_Ix5y_S_M3_vrr-2*oned2z*I_ESP_Ix5y_S_M4_vrr;
        Double I_ESP_L3x5z_S_M3_vrr = PAX*I_ESP_K2x5z_S_M3_vrr-PRX*I_ESP_K2x5z_S_M4_vrr+2*oned2z*I_ESP_Ix5z_S_M3_vrr-2*oned2z*I_ESP_Ix5z_S_M4_vrr;
        Double I_ESP_L2x6y_S_M3_vrr = PAX*I_ESP_Kx6y_S_M3_vrr-PRX*I_ESP_Kx6y_S_M4_vrr+oned2z*I_ESP_I6y_S_M3_vrr-oned2z*I_ESP_I6y_S_M4_vrr;
        Double I_ESP_L2x6z_S_M3_vrr = PAX*I_ESP_Kx6z_S_M3_vrr-PRX*I_ESP_Kx6z_S_M4_vrr+oned2z*I_ESP_I6z_S_M3_vrr-oned2z*I_ESP_I6z_S_M4_vrr;
        Double I_ESP_Lx7y_S_M3_vrr = PAX*I_ESP_K7y_S_M3_vrr-PRX*I_ESP_K7y_S_M4_vrr;
        Double I_ESP_Lx7z_S_M3_vrr = PAX*I_ESP_K7z_S_M3_vrr-PRX*I_ESP_K7z_S_M4_vrr;
        Double I_ESP_L8y_S_M3_vrr = PAY*I_ESP_K7y_S_M3_vrr-PRY*I_ESP_K7y_S_M4_vrr+7*oned2z*I_ESP_I6y_S_M3_vrr-7*oned2z*I_ESP_I6y_S_M4_vrr;
        Double I_ESP_L7yz_S_M3_vrr = PAZ*I_ESP_K7y_S_M3_vrr-PRZ*I_ESP_K7y_S_M4_vrr;
        Double I_ESP_L6y2z_S_M3_vrr = PAZ*I_ESP_K6yz_S_M3_vrr-PRZ*I_ESP_K6yz_S_M4_vrr+oned2z*I_ESP_I6y_S_M3_vrr-oned2z*I_ESP_I6y_S_M4_vrr;
        Double I_ESP_L5y3z_S_M3_vrr = PAZ*I_ESP_K5y2z_S_M3_vrr-PRZ*I_ESP_K5y2z_S_M4_vrr+2*oned2z*I_ESP_I5yz_S_M3_vrr-2*oned2z*I_ESP_I5yz_S_M4_vrr;
        Double I_ESP_L4y4z_S_M3_vrr = PAZ*I_ESP_K4y3z_S_M3_vrr-PRZ*I_ESP_K4y3z_S_M4_vrr+3*oned2z*I_ESP_I4y2z_S_M3_vrr-3*oned2z*I_ESP_I4y2z_S_M4_vrr;
        Double I_ESP_L3y5z_S_M3_vrr = PAY*I_ESP_K2y5z_S_M3_vrr-PRY*I_ESP_K2y5z_S_M4_vrr+2*oned2z*I_ESP_Iy5z_S_M3_vrr-2*oned2z*I_ESP_Iy5z_S_M4_vrr;
        Double I_ESP_L2y6z_S_M3_vrr = PAY*I_ESP_Ky6z_S_M3_vrr-PRY*I_ESP_Ky6z_S_M4_vrr+oned2z*I_ESP_I6z_S_M3_vrr-oned2z*I_ESP_I6z_S_M4_vrr;
        Double I_ESP_Ly7z_S_M3_vrr = PAY*I_ESP_K7z_S_M3_vrr-PRY*I_ESP_K7z_S_M4_vrr;
        Double I_ESP_L8z_S_M3_vrr = PAZ*I_ESP_K7z_S_M3_vrr-PRZ*I_ESP_K7z_S_M4_vrr+7*oned2z*I_ESP_I6z_S_M3_vrr-7*oned2z*I_ESP_I6z_S_M4_vrr;

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
         * shell quartet name: SQ_ESP_M_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 22 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_L_S_M2
         * RHS shell quartet name: SQ_ESP_L_S_M3
         * RHS shell quartet name: SQ_ESP_K_S_M2
         * RHS shell quartet name: SQ_ESP_K_S_M3
         ************************************************************/
        Double I_ESP_M9x_S_M2_vrr = PAX*I_ESP_L8x_S_M2_vrr-PRX*I_ESP_L8x_S_M3_vrr+8*oned2z*I_ESP_K7x_S_M2_vrr-8*oned2z*I_ESP_K7x_S_M3_vrr;
        Double I_ESP_M8xy_S_M2_vrr = PAY*I_ESP_L8x_S_M2_vrr-PRY*I_ESP_L8x_S_M3_vrr;
        Double I_ESP_M8xz_S_M2_vrr = PAZ*I_ESP_L8x_S_M2_vrr-PRZ*I_ESP_L8x_S_M3_vrr;
        Double I_ESP_M7x2y_S_M2_vrr = PAY*I_ESP_L7xy_S_M2_vrr-PRY*I_ESP_L7xy_S_M3_vrr+oned2z*I_ESP_K7x_S_M2_vrr-oned2z*I_ESP_K7x_S_M3_vrr;
        Double I_ESP_M7x2z_S_M2_vrr = PAZ*I_ESP_L7xz_S_M2_vrr-PRZ*I_ESP_L7xz_S_M3_vrr+oned2z*I_ESP_K7x_S_M2_vrr-oned2z*I_ESP_K7x_S_M3_vrr;
        Double I_ESP_M6x3y_S_M2_vrr = PAY*I_ESP_L6x2y_S_M2_vrr-PRY*I_ESP_L6x2y_S_M3_vrr+2*oned2z*I_ESP_K6xy_S_M2_vrr-2*oned2z*I_ESP_K6xy_S_M3_vrr;
        Double I_ESP_M6x3z_S_M2_vrr = PAZ*I_ESP_L6x2z_S_M2_vrr-PRZ*I_ESP_L6x2z_S_M3_vrr+2*oned2z*I_ESP_K6xz_S_M2_vrr-2*oned2z*I_ESP_K6xz_S_M3_vrr;
        Double I_ESP_M5x4y_S_M2_vrr = PAY*I_ESP_L5x3y_S_M2_vrr-PRY*I_ESP_L5x3y_S_M3_vrr+3*oned2z*I_ESP_K5x2y_S_M2_vrr-3*oned2z*I_ESP_K5x2y_S_M3_vrr;
        Double I_ESP_M5x3yz_S_M2_vrr = PAZ*I_ESP_L5x3y_S_M2_vrr-PRZ*I_ESP_L5x3y_S_M3_vrr;
        Double I_ESP_M5x4z_S_M2_vrr = PAZ*I_ESP_L5x3z_S_M2_vrr-PRZ*I_ESP_L5x3z_S_M3_vrr+3*oned2z*I_ESP_K5x2z_S_M2_vrr-3*oned2z*I_ESP_K5x2z_S_M3_vrr;
        Double I_ESP_M4x5y_S_M2_vrr = PAX*I_ESP_L3x5y_S_M2_vrr-PRX*I_ESP_L3x5y_S_M3_vrr+3*oned2z*I_ESP_K2x5y_S_M2_vrr-3*oned2z*I_ESP_K2x5y_S_M3_vrr;
        Double I_ESP_M4x4yz_S_M2_vrr = PAZ*I_ESP_L4x4y_S_M2_vrr-PRZ*I_ESP_L4x4y_S_M3_vrr;
        Double I_ESP_M4xy4z_S_M2_vrr = PAY*I_ESP_L4x4z_S_M2_vrr-PRY*I_ESP_L4x4z_S_M3_vrr;
        Double I_ESP_M4x5z_S_M2_vrr = PAX*I_ESP_L3x5z_S_M2_vrr-PRX*I_ESP_L3x5z_S_M3_vrr+3*oned2z*I_ESP_K2x5z_S_M2_vrr-3*oned2z*I_ESP_K2x5z_S_M3_vrr;
        Double I_ESP_M3x6y_S_M2_vrr = PAX*I_ESP_L2x6y_S_M2_vrr-PRX*I_ESP_L2x6y_S_M3_vrr+2*oned2z*I_ESP_Kx6y_S_M2_vrr-2*oned2z*I_ESP_Kx6y_S_M3_vrr;
        Double I_ESP_M3x5yz_S_M2_vrr = PAZ*I_ESP_L3x5y_S_M2_vrr-PRZ*I_ESP_L3x5y_S_M3_vrr;
        Double I_ESP_M3xy5z_S_M2_vrr = PAY*I_ESP_L3x5z_S_M2_vrr-PRY*I_ESP_L3x5z_S_M3_vrr;
        Double I_ESP_M3x6z_S_M2_vrr = PAX*I_ESP_L2x6z_S_M2_vrr-PRX*I_ESP_L2x6z_S_M3_vrr+2*oned2z*I_ESP_Kx6z_S_M2_vrr-2*oned2z*I_ESP_Kx6z_S_M3_vrr;
        Double I_ESP_M2x7y_S_M2_vrr = PAX*I_ESP_Lx7y_S_M2_vrr-PRX*I_ESP_Lx7y_S_M3_vrr+oned2z*I_ESP_K7y_S_M2_vrr-oned2z*I_ESP_K7y_S_M3_vrr;
        Double I_ESP_M2x7z_S_M2_vrr = PAX*I_ESP_Lx7z_S_M2_vrr-PRX*I_ESP_Lx7z_S_M3_vrr+oned2z*I_ESP_K7z_S_M2_vrr-oned2z*I_ESP_K7z_S_M3_vrr;
        Double I_ESP_Mx8y_S_M2_vrr = PAX*I_ESP_L8y_S_M2_vrr-PRX*I_ESP_L8y_S_M3_vrr;
        Double I_ESP_Mx4y4z_S_M2_vrr = PAX*I_ESP_L4y4z_S_M2_vrr-PRX*I_ESP_L4y4z_S_M3_vrr;
        Double I_ESP_Mx8z_S_M2_vrr = PAX*I_ESP_L8z_S_M2_vrr-PRX*I_ESP_L8z_S_M3_vrr;
        Double I_ESP_M9y_S_M2_vrr = PAY*I_ESP_L8y_S_M2_vrr-PRY*I_ESP_L8y_S_M3_vrr+8*oned2z*I_ESP_K7y_S_M2_vrr-8*oned2z*I_ESP_K7y_S_M3_vrr;
        Double I_ESP_M8yz_S_M2_vrr = PAZ*I_ESP_L8y_S_M2_vrr-PRZ*I_ESP_L8y_S_M3_vrr;
        Double I_ESP_M7y2z_S_M2_vrr = PAZ*I_ESP_L7yz_S_M2_vrr-PRZ*I_ESP_L7yz_S_M3_vrr+oned2z*I_ESP_K7y_S_M2_vrr-oned2z*I_ESP_K7y_S_M3_vrr;
        Double I_ESP_M6y3z_S_M2_vrr = PAZ*I_ESP_L6y2z_S_M2_vrr-PRZ*I_ESP_L6y2z_S_M3_vrr+2*oned2z*I_ESP_K6yz_S_M2_vrr-2*oned2z*I_ESP_K6yz_S_M3_vrr;
        Double I_ESP_M5y4z_S_M2_vrr = PAZ*I_ESP_L5y3z_S_M2_vrr-PRZ*I_ESP_L5y3z_S_M3_vrr+3*oned2z*I_ESP_K5y2z_S_M2_vrr-3*oned2z*I_ESP_K5y2z_S_M3_vrr;
        Double I_ESP_M4y5z_S_M2_vrr = PAY*I_ESP_L3y5z_S_M2_vrr-PRY*I_ESP_L3y5z_S_M3_vrr+3*oned2z*I_ESP_K2y5z_S_M2_vrr-3*oned2z*I_ESP_K2y5z_S_M3_vrr;
        Double I_ESP_M3y6z_S_M2_vrr = PAY*I_ESP_L2y6z_S_M2_vrr-PRY*I_ESP_L2y6z_S_M3_vrr+2*oned2z*I_ESP_Ky6z_S_M2_vrr-2*oned2z*I_ESP_Ky6z_S_M3_vrr;
        Double I_ESP_M2y7z_S_M2_vrr = PAY*I_ESP_Ly7z_S_M2_vrr-PRY*I_ESP_Ly7z_S_M3_vrr+oned2z*I_ESP_K7z_S_M2_vrr-oned2z*I_ESP_K7z_S_M3_vrr;
        Double I_ESP_My8z_S_M2_vrr = PAY*I_ESP_L8z_S_M2_vrr-PRY*I_ESP_L8z_S_M3_vrr;
        Double I_ESP_M9z_S_M2_vrr = PAZ*I_ESP_L8z_S_M2_vrr-PRZ*I_ESP_L8z_S_M3_vrr+8*oned2z*I_ESP_K7z_S_M2_vrr-8*oned2z*I_ESP_K7z_S_M3_vrr;

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
         * shell quartet name: SQ_ESP_N_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 15 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_M_S_M1
         * RHS shell quartet name: SQ_ESP_M_S_M2
         * RHS shell quartet name: SQ_ESP_L_S_M1
         * RHS shell quartet name: SQ_ESP_L_S_M2
         ************************************************************/
        Double I_ESP_N10x_S_M1_vrr = PAX*I_ESP_M9x_S_M1_vrr-PRX*I_ESP_M9x_S_M2_vrr+9*oned2z*I_ESP_L8x_S_M1_vrr-9*oned2z*I_ESP_L8x_S_M2_vrr;
        Double I_ESP_N9xy_S_M1_vrr = PAY*I_ESP_M9x_S_M1_vrr-PRY*I_ESP_M9x_S_M2_vrr;
        Double I_ESP_N9xz_S_M1_vrr = PAZ*I_ESP_M9x_S_M1_vrr-PRZ*I_ESP_M9x_S_M2_vrr;
        Double I_ESP_N8x2y_S_M1_vrr = PAY*I_ESP_M8xy_S_M1_vrr-PRY*I_ESP_M8xy_S_M2_vrr+oned2z*I_ESP_L8x_S_M1_vrr-oned2z*I_ESP_L8x_S_M2_vrr;
        Double I_ESP_N8x2z_S_M1_vrr = PAZ*I_ESP_M8xz_S_M1_vrr-PRZ*I_ESP_M8xz_S_M2_vrr+oned2z*I_ESP_L8x_S_M1_vrr-oned2z*I_ESP_L8x_S_M2_vrr;
        Double I_ESP_N7x3y_S_M1_vrr = PAY*I_ESP_M7x2y_S_M1_vrr-PRY*I_ESP_M7x2y_S_M2_vrr+2*oned2z*I_ESP_L7xy_S_M1_vrr-2*oned2z*I_ESP_L7xy_S_M2_vrr;
        Double I_ESP_N7x2yz_S_M1_vrr = PAZ*I_ESP_M7x2y_S_M1_vrr-PRZ*I_ESP_M7x2y_S_M2_vrr;
        Double I_ESP_N7x3z_S_M1_vrr = PAZ*I_ESP_M7x2z_S_M1_vrr-PRZ*I_ESP_M7x2z_S_M2_vrr+2*oned2z*I_ESP_L7xz_S_M1_vrr-2*oned2z*I_ESP_L7xz_S_M2_vrr;
        Double I_ESP_N6x4y_S_M1_vrr = PAY*I_ESP_M6x3y_S_M1_vrr-PRY*I_ESP_M6x3y_S_M2_vrr+3*oned2z*I_ESP_L6x2y_S_M1_vrr-3*oned2z*I_ESP_L6x2y_S_M2_vrr;
        Double I_ESP_N6x3yz_S_M1_vrr = PAZ*I_ESP_M6x3y_S_M1_vrr-PRZ*I_ESP_M6x3y_S_M2_vrr;
        Double I_ESP_N6xy3z_S_M1_vrr = PAY*I_ESP_M6x3z_S_M1_vrr-PRY*I_ESP_M6x3z_S_M2_vrr;
        Double I_ESP_N6x4z_S_M1_vrr = PAZ*I_ESP_M6x3z_S_M1_vrr-PRZ*I_ESP_M6x3z_S_M2_vrr+3*oned2z*I_ESP_L6x2z_S_M1_vrr-3*oned2z*I_ESP_L6x2z_S_M2_vrr;
        Double I_ESP_N5x5y_S_M1_vrr = PAY*I_ESP_M5x4y_S_M1_vrr-PRY*I_ESP_M5x4y_S_M2_vrr+4*oned2z*I_ESP_L5x3y_S_M1_vrr-4*oned2z*I_ESP_L5x3y_S_M2_vrr;
        Double I_ESP_N5x4yz_S_M1_vrr = PAZ*I_ESP_M5x4y_S_M1_vrr-PRZ*I_ESP_M5x4y_S_M2_vrr;
        Double I_ESP_N5x3y2z_S_M1_vrr = PAZ*I_ESP_M5x3yz_S_M1_vrr-PRZ*I_ESP_M5x3yz_S_M2_vrr+oned2z*I_ESP_L5x3y_S_M1_vrr-oned2z*I_ESP_L5x3y_S_M2_vrr;
        Double I_ESP_N5xy4z_S_M1_vrr = PAY*I_ESP_M5x4z_S_M1_vrr-PRY*I_ESP_M5x4z_S_M2_vrr;
        Double I_ESP_N5x5z_S_M1_vrr = PAZ*I_ESP_M5x4z_S_M1_vrr-PRZ*I_ESP_M5x4z_S_M2_vrr+4*oned2z*I_ESP_L5x3z_S_M1_vrr-4*oned2z*I_ESP_L5x3z_S_M2_vrr;
        Double I_ESP_N4x6y_S_M1_vrr = PAX*I_ESP_M3x6y_S_M1_vrr-PRX*I_ESP_M3x6y_S_M2_vrr+3*oned2z*I_ESP_L2x6y_S_M1_vrr-3*oned2z*I_ESP_L2x6y_S_M2_vrr;
        Double I_ESP_N4x5yz_S_M1_vrr = PAZ*I_ESP_M4x5y_S_M1_vrr-PRZ*I_ESP_M4x5y_S_M2_vrr;
        Double I_ESP_N4x4y2z_S_M1_vrr = PAZ*I_ESP_M4x4yz_S_M1_vrr-PRZ*I_ESP_M4x4yz_S_M2_vrr+oned2z*I_ESP_L4x4y_S_M1_vrr-oned2z*I_ESP_L4x4y_S_M2_vrr;
        Double I_ESP_N4x2y4z_S_M1_vrr = PAY*I_ESP_M4xy4z_S_M1_vrr-PRY*I_ESP_M4xy4z_S_M2_vrr+oned2z*I_ESP_L4x4z_S_M1_vrr-oned2z*I_ESP_L4x4z_S_M2_vrr;
        Double I_ESP_N4xy5z_S_M1_vrr = PAY*I_ESP_M4x5z_S_M1_vrr-PRY*I_ESP_M4x5z_S_M2_vrr;
        Double I_ESP_N4x6z_S_M1_vrr = PAX*I_ESP_M3x6z_S_M1_vrr-PRX*I_ESP_M3x6z_S_M2_vrr+3*oned2z*I_ESP_L2x6z_S_M1_vrr-3*oned2z*I_ESP_L2x6z_S_M2_vrr;
        Double I_ESP_N3x7y_S_M1_vrr = PAX*I_ESP_M2x7y_S_M1_vrr-PRX*I_ESP_M2x7y_S_M2_vrr+2*oned2z*I_ESP_Lx7y_S_M1_vrr-2*oned2z*I_ESP_Lx7y_S_M2_vrr;
        Double I_ESP_N3x6yz_S_M1_vrr = PAZ*I_ESP_M3x6y_S_M1_vrr-PRZ*I_ESP_M3x6y_S_M2_vrr;
        Double I_ESP_N3x5y2z_S_M1_vrr = PAZ*I_ESP_M3x5yz_S_M1_vrr-PRZ*I_ESP_M3x5yz_S_M2_vrr+oned2z*I_ESP_L3x5y_S_M1_vrr-oned2z*I_ESP_L3x5y_S_M2_vrr;
        Double I_ESP_N3x2y5z_S_M1_vrr = PAY*I_ESP_M3xy5z_S_M1_vrr-PRY*I_ESP_M3xy5z_S_M2_vrr+oned2z*I_ESP_L3x5z_S_M1_vrr-oned2z*I_ESP_L3x5z_S_M2_vrr;
        Double I_ESP_N3xy6z_S_M1_vrr = PAY*I_ESP_M3x6z_S_M1_vrr-PRY*I_ESP_M3x6z_S_M2_vrr;
        Double I_ESP_N3x7z_S_M1_vrr = PAX*I_ESP_M2x7z_S_M1_vrr-PRX*I_ESP_M2x7z_S_M2_vrr+2*oned2z*I_ESP_Lx7z_S_M1_vrr-2*oned2z*I_ESP_Lx7z_S_M2_vrr;
        Double I_ESP_N2x8y_S_M1_vrr = PAX*I_ESP_Mx8y_S_M1_vrr-PRX*I_ESP_Mx8y_S_M2_vrr+oned2z*I_ESP_L8y_S_M1_vrr-oned2z*I_ESP_L8y_S_M2_vrr;
        Double I_ESP_N2x7yz_S_M1_vrr = PAZ*I_ESP_M2x7y_S_M1_vrr-PRZ*I_ESP_M2x7y_S_M2_vrr;
        Double I_ESP_N2x4y4z_S_M1_vrr = PAX*I_ESP_Mx4y4z_S_M1_vrr-PRX*I_ESP_Mx4y4z_S_M2_vrr+oned2z*I_ESP_L4y4z_S_M1_vrr-oned2z*I_ESP_L4y4z_S_M2_vrr;
        Double I_ESP_N2xy7z_S_M1_vrr = PAY*I_ESP_M2x7z_S_M1_vrr-PRY*I_ESP_M2x7z_S_M2_vrr;
        Double I_ESP_N2x8z_S_M1_vrr = PAX*I_ESP_Mx8z_S_M1_vrr-PRX*I_ESP_Mx8z_S_M2_vrr+oned2z*I_ESP_L8z_S_M1_vrr-oned2z*I_ESP_L8z_S_M2_vrr;
        Double I_ESP_Nx9y_S_M1_vrr = PAX*I_ESP_M9y_S_M1_vrr-PRX*I_ESP_M9y_S_M2_vrr;
        Double I_ESP_Nx6y3z_S_M1_vrr = PAX*I_ESP_M6y3z_S_M1_vrr-PRX*I_ESP_M6y3z_S_M2_vrr;
        Double I_ESP_Nx5y4z_S_M1_vrr = PAX*I_ESP_M5y4z_S_M1_vrr-PRX*I_ESP_M5y4z_S_M2_vrr;
        Double I_ESP_Nx4y5z_S_M1_vrr = PAX*I_ESP_M4y5z_S_M1_vrr-PRX*I_ESP_M4y5z_S_M2_vrr;
        Double I_ESP_Nx3y6z_S_M1_vrr = PAX*I_ESP_M3y6z_S_M1_vrr-PRX*I_ESP_M3y6z_S_M2_vrr;
        Double I_ESP_Nx9z_S_M1_vrr = PAX*I_ESP_M9z_S_M1_vrr-PRX*I_ESP_M9z_S_M2_vrr;
        Double I_ESP_N10y_S_M1_vrr = PAY*I_ESP_M9y_S_M1_vrr-PRY*I_ESP_M9y_S_M2_vrr+9*oned2z*I_ESP_L8y_S_M1_vrr-9*oned2z*I_ESP_L8y_S_M2_vrr;
        Double I_ESP_N9yz_S_M1_vrr = PAZ*I_ESP_M9y_S_M1_vrr-PRZ*I_ESP_M9y_S_M2_vrr;
        Double I_ESP_N8y2z_S_M1_vrr = PAZ*I_ESP_M8yz_S_M1_vrr-PRZ*I_ESP_M8yz_S_M2_vrr+oned2z*I_ESP_L8y_S_M1_vrr-oned2z*I_ESP_L8y_S_M2_vrr;
        Double I_ESP_N7y3z_S_M1_vrr = PAZ*I_ESP_M7y2z_S_M1_vrr-PRZ*I_ESP_M7y2z_S_M2_vrr+2*oned2z*I_ESP_L7yz_S_M1_vrr-2*oned2z*I_ESP_L7yz_S_M2_vrr;
        Double I_ESP_N6y4z_S_M1_vrr = PAZ*I_ESP_M6y3z_S_M1_vrr-PRZ*I_ESP_M6y3z_S_M2_vrr+3*oned2z*I_ESP_L6y2z_S_M1_vrr-3*oned2z*I_ESP_L6y2z_S_M2_vrr;
        Double I_ESP_N5y5z_S_M1_vrr = PAZ*I_ESP_M5y4z_S_M1_vrr-PRZ*I_ESP_M5y4z_S_M2_vrr+4*oned2z*I_ESP_L5y3z_S_M1_vrr-4*oned2z*I_ESP_L5y3z_S_M2_vrr;
        Double I_ESP_N4y6z_S_M1_vrr = PAY*I_ESP_M3y6z_S_M1_vrr-PRY*I_ESP_M3y6z_S_M2_vrr+3*oned2z*I_ESP_L2y6z_S_M1_vrr-3*oned2z*I_ESP_L2y6z_S_M2_vrr;
        Double I_ESP_N3y7z_S_M1_vrr = PAY*I_ESP_M2y7z_S_M1_vrr-PRY*I_ESP_M2y7z_S_M2_vrr+2*oned2z*I_ESP_Ly7z_S_M1_vrr-2*oned2z*I_ESP_Ly7z_S_M2_vrr;
        Double I_ESP_N2y8z_S_M1_vrr = PAY*I_ESP_My8z_S_M1_vrr-PRY*I_ESP_My8z_S_M2_vrr+oned2z*I_ESP_L8z_S_M1_vrr-oned2z*I_ESP_L8z_S_M2_vrr;
        Double I_ESP_Ny9z_S_M1_vrr = PAY*I_ESP_M9z_S_M1_vrr-PRY*I_ESP_M9z_S_M2_vrr;
        Double I_ESP_N10z_S_M1_vrr = PAZ*I_ESP_M9z_S_M1_vrr-PRZ*I_ESP_M9z_S_M2_vrr+9*oned2z*I_ESP_L8z_S_M1_vrr-9*oned2z*I_ESP_L8z_S_M2_vrr;

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
         * shell quartet name: SQ_ESP_O_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_N_S
         * RHS shell quartet name: SQ_ESP_N_S_M1
         * RHS shell quartet name: SQ_ESP_M_S
         * RHS shell quartet name: SQ_ESP_M_S_M1
         ************************************************************/
        Double I_ESP_O11x_S_vrr = PAX*I_ESP_N10x_S_vrr-PRX*I_ESP_N10x_S_M1_vrr+10*oned2z*I_ESP_M9x_S_vrr-10*oned2z*I_ESP_M9x_S_M1_vrr;
        Double I_ESP_O10xy_S_vrr = PAY*I_ESP_N10x_S_vrr-PRY*I_ESP_N10x_S_M1_vrr;
        Double I_ESP_O10xz_S_vrr = PAZ*I_ESP_N10x_S_vrr-PRZ*I_ESP_N10x_S_M1_vrr;
        Double I_ESP_O9x2y_S_vrr = PAY*I_ESP_N9xy_S_vrr-PRY*I_ESP_N9xy_S_M1_vrr+oned2z*I_ESP_M9x_S_vrr-oned2z*I_ESP_M9x_S_M1_vrr;
        Double I_ESP_O9xyz_S_vrr = PAZ*I_ESP_N9xy_S_vrr-PRZ*I_ESP_N9xy_S_M1_vrr;
        Double I_ESP_O9x2z_S_vrr = PAZ*I_ESP_N9xz_S_vrr-PRZ*I_ESP_N9xz_S_M1_vrr+oned2z*I_ESP_M9x_S_vrr-oned2z*I_ESP_M9x_S_M1_vrr;
        Double I_ESP_O8x3y_S_vrr = PAY*I_ESP_N8x2y_S_vrr-PRY*I_ESP_N8x2y_S_M1_vrr+2*oned2z*I_ESP_M8xy_S_vrr-2*oned2z*I_ESP_M8xy_S_M1_vrr;
        Double I_ESP_O8x2yz_S_vrr = PAZ*I_ESP_N8x2y_S_vrr-PRZ*I_ESP_N8x2y_S_M1_vrr;
        Double I_ESP_O8xy2z_S_vrr = PAY*I_ESP_N8x2z_S_vrr-PRY*I_ESP_N8x2z_S_M1_vrr;
        Double I_ESP_O8x3z_S_vrr = PAZ*I_ESP_N8x2z_S_vrr-PRZ*I_ESP_N8x2z_S_M1_vrr+2*oned2z*I_ESP_M8xz_S_vrr-2*oned2z*I_ESP_M8xz_S_M1_vrr;
        Double I_ESP_O7x4y_S_vrr = PAY*I_ESP_N7x3y_S_vrr-PRY*I_ESP_N7x3y_S_M1_vrr+3*oned2z*I_ESP_M7x2y_S_vrr-3*oned2z*I_ESP_M7x2y_S_M1_vrr;
        Double I_ESP_O7x3yz_S_vrr = PAZ*I_ESP_N7x3y_S_vrr-PRZ*I_ESP_N7x3y_S_M1_vrr;
        Double I_ESP_O7x2y2z_S_vrr = PAZ*I_ESP_N7x2yz_S_vrr-PRZ*I_ESP_N7x2yz_S_M1_vrr+oned2z*I_ESP_M7x2y_S_vrr-oned2z*I_ESP_M7x2y_S_M1_vrr;
        Double I_ESP_O7xy3z_S_vrr = PAY*I_ESP_N7x3z_S_vrr-PRY*I_ESP_N7x3z_S_M1_vrr;
        Double I_ESP_O7x4z_S_vrr = PAZ*I_ESP_N7x3z_S_vrr-PRZ*I_ESP_N7x3z_S_M1_vrr+3*oned2z*I_ESP_M7x2z_S_vrr-3*oned2z*I_ESP_M7x2z_S_M1_vrr;
        Double I_ESP_O6x5y_S_vrr = PAY*I_ESP_N6x4y_S_vrr-PRY*I_ESP_N6x4y_S_M1_vrr+4*oned2z*I_ESP_M6x3y_S_vrr-4*oned2z*I_ESP_M6x3y_S_M1_vrr;
        Double I_ESP_O6x4yz_S_vrr = PAZ*I_ESP_N6x4y_S_vrr-PRZ*I_ESP_N6x4y_S_M1_vrr;
        Double I_ESP_O6x3y2z_S_vrr = PAZ*I_ESP_N6x3yz_S_vrr-PRZ*I_ESP_N6x3yz_S_M1_vrr+oned2z*I_ESP_M6x3y_S_vrr-oned2z*I_ESP_M6x3y_S_M1_vrr;
        Double I_ESP_O6x2y3z_S_vrr = PAY*I_ESP_N6xy3z_S_vrr-PRY*I_ESP_N6xy3z_S_M1_vrr+oned2z*I_ESP_M6x3z_S_vrr-oned2z*I_ESP_M6x3z_S_M1_vrr;
        Double I_ESP_O6xy4z_S_vrr = PAY*I_ESP_N6x4z_S_vrr-PRY*I_ESP_N6x4z_S_M1_vrr;
        Double I_ESP_O6x5z_S_vrr = PAZ*I_ESP_N6x4z_S_vrr-PRZ*I_ESP_N6x4z_S_M1_vrr+4*oned2z*I_ESP_M6x3z_S_vrr-4*oned2z*I_ESP_M6x3z_S_M1_vrr;
        Double I_ESP_O5x6y_S_vrr = PAX*I_ESP_N4x6y_S_vrr-PRX*I_ESP_N4x6y_S_M1_vrr+4*oned2z*I_ESP_M3x6y_S_vrr-4*oned2z*I_ESP_M3x6y_S_M1_vrr;
        Double I_ESP_O5x5yz_S_vrr = PAZ*I_ESP_N5x5y_S_vrr-PRZ*I_ESP_N5x5y_S_M1_vrr;
        Double I_ESP_O5x4y2z_S_vrr = PAZ*I_ESP_N5x4yz_S_vrr-PRZ*I_ESP_N5x4yz_S_M1_vrr+oned2z*I_ESP_M5x4y_S_vrr-oned2z*I_ESP_M5x4y_S_M1_vrr;
        Double I_ESP_O5x3y3z_S_vrr = PAZ*I_ESP_N5x3y2z_S_vrr-PRZ*I_ESP_N5x3y2z_S_M1_vrr+2*oned2z*I_ESP_M5x3yz_S_vrr-2*oned2z*I_ESP_M5x3yz_S_M1_vrr;
        Double I_ESP_O5x2y4z_S_vrr = PAY*I_ESP_N5xy4z_S_vrr-PRY*I_ESP_N5xy4z_S_M1_vrr+oned2z*I_ESP_M5x4z_S_vrr-oned2z*I_ESP_M5x4z_S_M1_vrr;
        Double I_ESP_O5xy5z_S_vrr = PAY*I_ESP_N5x5z_S_vrr-PRY*I_ESP_N5x5z_S_M1_vrr;
        Double I_ESP_O5x6z_S_vrr = PAX*I_ESP_N4x6z_S_vrr-PRX*I_ESP_N4x6z_S_M1_vrr+4*oned2z*I_ESP_M3x6z_S_vrr-4*oned2z*I_ESP_M3x6z_S_M1_vrr;
        Double I_ESP_O4x7y_S_vrr = PAX*I_ESP_N3x7y_S_vrr-PRX*I_ESP_N3x7y_S_M1_vrr+3*oned2z*I_ESP_M2x7y_S_vrr-3*oned2z*I_ESP_M2x7y_S_M1_vrr;
        Double I_ESP_O4x6yz_S_vrr = PAZ*I_ESP_N4x6y_S_vrr-PRZ*I_ESP_N4x6y_S_M1_vrr;
        Double I_ESP_O4x5y2z_S_vrr = PAZ*I_ESP_N4x5yz_S_vrr-PRZ*I_ESP_N4x5yz_S_M1_vrr+oned2z*I_ESP_M4x5y_S_vrr-oned2z*I_ESP_M4x5y_S_M1_vrr;
        Double I_ESP_O4x4y3z_S_vrr = PAZ*I_ESP_N4x4y2z_S_vrr-PRZ*I_ESP_N4x4y2z_S_M1_vrr+2*oned2z*I_ESP_M4x4yz_S_vrr-2*oned2z*I_ESP_M4x4yz_S_M1_vrr;
        Double I_ESP_O4x3y4z_S_vrr = PAY*I_ESP_N4x2y4z_S_vrr-PRY*I_ESP_N4x2y4z_S_M1_vrr+2*oned2z*I_ESP_M4xy4z_S_vrr-2*oned2z*I_ESP_M4xy4z_S_M1_vrr;
        Double I_ESP_O4x2y5z_S_vrr = PAY*I_ESP_N4xy5z_S_vrr-PRY*I_ESP_N4xy5z_S_M1_vrr+oned2z*I_ESP_M4x5z_S_vrr-oned2z*I_ESP_M4x5z_S_M1_vrr;
        Double I_ESP_O4xy6z_S_vrr = PAY*I_ESP_N4x6z_S_vrr-PRY*I_ESP_N4x6z_S_M1_vrr;
        Double I_ESP_O4x7z_S_vrr = PAX*I_ESP_N3x7z_S_vrr-PRX*I_ESP_N3x7z_S_M1_vrr+3*oned2z*I_ESP_M2x7z_S_vrr-3*oned2z*I_ESP_M2x7z_S_M1_vrr;
        Double I_ESP_O3x8y_S_vrr = PAX*I_ESP_N2x8y_S_vrr-PRX*I_ESP_N2x8y_S_M1_vrr+2*oned2z*I_ESP_Mx8y_S_vrr-2*oned2z*I_ESP_Mx8y_S_M1_vrr;
        Double I_ESP_O3x7yz_S_vrr = PAZ*I_ESP_N3x7y_S_vrr-PRZ*I_ESP_N3x7y_S_M1_vrr;
        Double I_ESP_O3x6y2z_S_vrr = PAZ*I_ESP_N3x6yz_S_vrr-PRZ*I_ESP_N3x6yz_S_M1_vrr+oned2z*I_ESP_M3x6y_S_vrr-oned2z*I_ESP_M3x6y_S_M1_vrr;
        Double I_ESP_O3x5y3z_S_vrr = PAZ*I_ESP_N3x5y2z_S_vrr-PRZ*I_ESP_N3x5y2z_S_M1_vrr+2*oned2z*I_ESP_M3x5yz_S_vrr-2*oned2z*I_ESP_M3x5yz_S_M1_vrr;
        Double I_ESP_O3x4y4z_S_vrr = PAX*I_ESP_N2x4y4z_S_vrr-PRX*I_ESP_N2x4y4z_S_M1_vrr+2*oned2z*I_ESP_Mx4y4z_S_vrr-2*oned2z*I_ESP_Mx4y4z_S_M1_vrr;
        Double I_ESP_O3x3y5z_S_vrr = PAY*I_ESP_N3x2y5z_S_vrr-PRY*I_ESP_N3x2y5z_S_M1_vrr+2*oned2z*I_ESP_M3xy5z_S_vrr-2*oned2z*I_ESP_M3xy5z_S_M1_vrr;
        Double I_ESP_O3x2y6z_S_vrr = PAY*I_ESP_N3xy6z_S_vrr-PRY*I_ESP_N3xy6z_S_M1_vrr+oned2z*I_ESP_M3x6z_S_vrr-oned2z*I_ESP_M3x6z_S_M1_vrr;
        Double I_ESP_O3xy7z_S_vrr = PAY*I_ESP_N3x7z_S_vrr-PRY*I_ESP_N3x7z_S_M1_vrr;
        Double I_ESP_O3x8z_S_vrr = PAX*I_ESP_N2x8z_S_vrr-PRX*I_ESP_N2x8z_S_M1_vrr+2*oned2z*I_ESP_Mx8z_S_vrr-2*oned2z*I_ESP_Mx8z_S_M1_vrr;
        Double I_ESP_O2x9y_S_vrr = PAX*I_ESP_Nx9y_S_vrr-PRX*I_ESP_Nx9y_S_M1_vrr+oned2z*I_ESP_M9y_S_vrr-oned2z*I_ESP_M9y_S_M1_vrr;
        Double I_ESP_O2x8yz_S_vrr = PAZ*I_ESP_N2x8y_S_vrr-PRZ*I_ESP_N2x8y_S_M1_vrr;
        Double I_ESP_O2x7y2z_S_vrr = PAZ*I_ESP_N2x7yz_S_vrr-PRZ*I_ESP_N2x7yz_S_M1_vrr+oned2z*I_ESP_M2x7y_S_vrr-oned2z*I_ESP_M2x7y_S_M1_vrr;
        Double I_ESP_O2x6y3z_S_vrr = PAX*I_ESP_Nx6y3z_S_vrr-PRX*I_ESP_Nx6y3z_S_M1_vrr+oned2z*I_ESP_M6y3z_S_vrr-oned2z*I_ESP_M6y3z_S_M1_vrr;
        Double I_ESP_O2x5y4z_S_vrr = PAX*I_ESP_Nx5y4z_S_vrr-PRX*I_ESP_Nx5y4z_S_M1_vrr+oned2z*I_ESP_M5y4z_S_vrr-oned2z*I_ESP_M5y4z_S_M1_vrr;
        Double I_ESP_O2x4y5z_S_vrr = PAX*I_ESP_Nx4y5z_S_vrr-PRX*I_ESP_Nx4y5z_S_M1_vrr+oned2z*I_ESP_M4y5z_S_vrr-oned2z*I_ESP_M4y5z_S_M1_vrr;
        Double I_ESP_O2x3y6z_S_vrr = PAX*I_ESP_Nx3y6z_S_vrr-PRX*I_ESP_Nx3y6z_S_M1_vrr+oned2z*I_ESP_M3y6z_S_vrr-oned2z*I_ESP_M3y6z_S_M1_vrr;
        Double I_ESP_O2x2y7z_S_vrr = PAY*I_ESP_N2xy7z_S_vrr-PRY*I_ESP_N2xy7z_S_M1_vrr+oned2z*I_ESP_M2x7z_S_vrr-oned2z*I_ESP_M2x7z_S_M1_vrr;
        Double I_ESP_O2xy8z_S_vrr = PAY*I_ESP_N2x8z_S_vrr-PRY*I_ESP_N2x8z_S_M1_vrr;
        Double I_ESP_O2x9z_S_vrr = PAX*I_ESP_Nx9z_S_vrr-PRX*I_ESP_Nx9z_S_M1_vrr+oned2z*I_ESP_M9z_S_vrr-oned2z*I_ESP_M9z_S_M1_vrr;
        Double I_ESP_Ox10y_S_vrr = PAX*I_ESP_N10y_S_vrr-PRX*I_ESP_N10y_S_M1_vrr;
        Double I_ESP_Ox9yz_S_vrr = PAZ*I_ESP_Nx9y_S_vrr-PRZ*I_ESP_Nx9y_S_M1_vrr;
        Double I_ESP_Ox8y2z_S_vrr = PAX*I_ESP_N8y2z_S_vrr-PRX*I_ESP_N8y2z_S_M1_vrr;
        Double I_ESP_Ox7y3z_S_vrr = PAX*I_ESP_N7y3z_S_vrr-PRX*I_ESP_N7y3z_S_M1_vrr;
        Double I_ESP_Ox6y4z_S_vrr = PAX*I_ESP_N6y4z_S_vrr-PRX*I_ESP_N6y4z_S_M1_vrr;
        Double I_ESP_Ox5y5z_S_vrr = PAX*I_ESP_N5y5z_S_vrr-PRX*I_ESP_N5y5z_S_M1_vrr;
        Double I_ESP_Ox4y6z_S_vrr = PAX*I_ESP_N4y6z_S_vrr-PRX*I_ESP_N4y6z_S_M1_vrr;
        Double I_ESP_Ox3y7z_S_vrr = PAX*I_ESP_N3y7z_S_vrr-PRX*I_ESP_N3y7z_S_M1_vrr;
        Double I_ESP_Ox2y8z_S_vrr = PAX*I_ESP_N2y8z_S_vrr-PRX*I_ESP_N2y8z_S_M1_vrr;
        Double I_ESP_Oxy9z_S_vrr = PAY*I_ESP_Nx9z_S_vrr-PRY*I_ESP_Nx9z_S_M1_vrr;
        Double I_ESP_Ox10z_S_vrr = PAX*I_ESP_N10z_S_vrr-PRX*I_ESP_N10z_S_M1_vrr;
        Double I_ESP_O11y_S_vrr = PAY*I_ESP_N10y_S_vrr-PRY*I_ESP_N10y_S_M1_vrr+10*oned2z*I_ESP_M9y_S_vrr-10*oned2z*I_ESP_M9y_S_M1_vrr;
        Double I_ESP_O10yz_S_vrr = PAZ*I_ESP_N10y_S_vrr-PRZ*I_ESP_N10y_S_M1_vrr;
        Double I_ESP_O9y2z_S_vrr = PAZ*I_ESP_N9yz_S_vrr-PRZ*I_ESP_N9yz_S_M1_vrr+oned2z*I_ESP_M9y_S_vrr-oned2z*I_ESP_M9y_S_M1_vrr;
        Double I_ESP_O8y3z_S_vrr = PAZ*I_ESP_N8y2z_S_vrr-PRZ*I_ESP_N8y2z_S_M1_vrr+2*oned2z*I_ESP_M8yz_S_vrr-2*oned2z*I_ESP_M8yz_S_M1_vrr;
        Double I_ESP_O7y4z_S_vrr = PAZ*I_ESP_N7y3z_S_vrr-PRZ*I_ESP_N7y3z_S_M1_vrr+3*oned2z*I_ESP_M7y2z_S_vrr-3*oned2z*I_ESP_M7y2z_S_M1_vrr;
        Double I_ESP_O6y5z_S_vrr = PAZ*I_ESP_N6y4z_S_vrr-PRZ*I_ESP_N6y4z_S_M1_vrr+4*oned2z*I_ESP_M6y3z_S_vrr-4*oned2z*I_ESP_M6y3z_S_M1_vrr;
        Double I_ESP_O5y6z_S_vrr = PAY*I_ESP_N4y6z_S_vrr-PRY*I_ESP_N4y6z_S_M1_vrr+4*oned2z*I_ESP_M3y6z_S_vrr-4*oned2z*I_ESP_M3y6z_S_M1_vrr;
        Double I_ESP_O4y7z_S_vrr = PAY*I_ESP_N3y7z_S_vrr-PRY*I_ESP_N3y7z_S_M1_vrr+3*oned2z*I_ESP_M2y7z_S_vrr-3*oned2z*I_ESP_M2y7z_S_M1_vrr;
        Double I_ESP_O3y8z_S_vrr = PAY*I_ESP_N2y8z_S_vrr-PRY*I_ESP_N2y8z_S_M1_vrr+2*oned2z*I_ESP_My8z_S_vrr-2*oned2z*I_ESP_My8z_S_M1_vrr;
        Double I_ESP_O2y9z_S_vrr = PAY*I_ESP_Ny9z_S_vrr-PRY*I_ESP_Ny9z_S_M1_vrr+oned2z*I_ESP_M9z_S_vrr-oned2z*I_ESP_M9z_S_M1_vrr;
        Double I_ESP_Oy10z_S_vrr = PAY*I_ESP_N10z_S_vrr-PRY*I_ESP_N10z_S_M1_vrr;
        Double I_ESP_O11z_S_vrr = PAZ*I_ESP_N10z_S_vrr-PRZ*I_ESP_N10z_S_M1_vrr+10*oned2z*I_ESP_M9z_S_vrr-10*oned2z*I_ESP_M9z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_O_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_O_S_aa_coefs = alpha*alpha;
        I_ESP_O11x_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O11x_S_vrr;
        I_ESP_O10xy_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O10xy_S_vrr;
        I_ESP_O10xz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O10xz_S_vrr;
        I_ESP_O9x2y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O9x2y_S_vrr;
        I_ESP_O9xyz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O9xyz_S_vrr;
        I_ESP_O9x2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O9x2z_S_vrr;
        I_ESP_O8x3y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O8x3y_S_vrr;
        I_ESP_O8x2yz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O8x2yz_S_vrr;
        I_ESP_O8xy2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O8xy2z_S_vrr;
        I_ESP_O8x3z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O8x3z_S_vrr;
        I_ESP_O7x4y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O7x4y_S_vrr;
        I_ESP_O7x3yz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O7x3yz_S_vrr;
        I_ESP_O7x2y2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O7x2y2z_S_vrr;
        I_ESP_O7xy3z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O7xy3z_S_vrr;
        I_ESP_O7x4z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O7x4z_S_vrr;
        I_ESP_O6x5y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O6x5y_S_vrr;
        I_ESP_O6x4yz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O6x4yz_S_vrr;
        I_ESP_O6x3y2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O6x3y2z_S_vrr;
        I_ESP_O6x2y3z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O6x2y3z_S_vrr;
        I_ESP_O6xy4z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O6xy4z_S_vrr;
        I_ESP_O6x5z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O6x5z_S_vrr;
        I_ESP_O5x6y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O5x6y_S_vrr;
        I_ESP_O5x5yz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O5x5yz_S_vrr;
        I_ESP_O5x4y2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O5x4y2z_S_vrr;
        I_ESP_O5x3y3z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O5x3y3z_S_vrr;
        I_ESP_O5x2y4z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O5x2y4z_S_vrr;
        I_ESP_O5xy5z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O5xy5z_S_vrr;
        I_ESP_O5x6z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O5x6z_S_vrr;
        I_ESP_O4x7y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O4x7y_S_vrr;
        I_ESP_O4x6yz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O4x6yz_S_vrr;
        I_ESP_O4x5y2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O4x5y2z_S_vrr;
        I_ESP_O4x4y3z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O4x4y3z_S_vrr;
        I_ESP_O4x3y4z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O4x3y4z_S_vrr;
        I_ESP_O4x2y5z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O4x2y5z_S_vrr;
        I_ESP_O4xy6z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O4xy6z_S_vrr;
        I_ESP_O4x7z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O4x7z_S_vrr;
        I_ESP_O3x8y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3x8y_S_vrr;
        I_ESP_O3x7yz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3x7yz_S_vrr;
        I_ESP_O3x6y2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3x6y2z_S_vrr;
        I_ESP_O3x5y3z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3x5y3z_S_vrr;
        I_ESP_O3x4y4z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3x4y4z_S_vrr;
        I_ESP_O3x3y5z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3x3y5z_S_vrr;
        I_ESP_O3x2y6z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3x2y6z_S_vrr;
        I_ESP_O3xy7z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3xy7z_S_vrr;
        I_ESP_O3x8z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3x8z_S_vrr;
        I_ESP_O2x9y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2x9y_S_vrr;
        I_ESP_O2x8yz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2x8yz_S_vrr;
        I_ESP_O2x7y2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2x7y2z_S_vrr;
        I_ESP_O2x6y3z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2x6y3z_S_vrr;
        I_ESP_O2x5y4z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2x5y4z_S_vrr;
        I_ESP_O2x4y5z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2x4y5z_S_vrr;
        I_ESP_O2x3y6z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2x3y6z_S_vrr;
        I_ESP_O2x2y7z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2x2y7z_S_vrr;
        I_ESP_O2xy8z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2xy8z_S_vrr;
        I_ESP_O2x9z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2x9z_S_vrr;
        I_ESP_Ox10y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox10y_S_vrr;
        I_ESP_Ox9yz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox9yz_S_vrr;
        I_ESP_Ox8y2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox8y2z_S_vrr;
        I_ESP_Ox7y3z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox7y3z_S_vrr;
        I_ESP_Ox6y4z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox6y4z_S_vrr;
        I_ESP_Ox5y5z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox5y5z_S_vrr;
        I_ESP_Ox4y6z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox4y6z_S_vrr;
        I_ESP_Ox3y7z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox3y7z_S_vrr;
        I_ESP_Ox2y8z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox2y8z_S_vrr;
        I_ESP_Oxy9z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Oxy9z_S_vrr;
        I_ESP_Ox10z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Ox10z_S_vrr;
        I_ESP_O11y_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O11y_S_vrr;
        I_ESP_O10yz_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O10yz_S_vrr;
        I_ESP_O9y2z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O9y2z_S_vrr;
        I_ESP_O8y3z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O8y3z_S_vrr;
        I_ESP_O7y4z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O7y4z_S_vrr;
        I_ESP_O6y5z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O6y5z_S_vrr;
        I_ESP_O5y6z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O5y6z_S_vrr;
        I_ESP_O4y7z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O4y7z_S_vrr;
        I_ESP_O3y8z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O3y8z_S_vrr;
        I_ESP_O2y9z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O2y9z_S_vrr;
        I_ESP_Oy10z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_Oy10z_S_vrr;
        I_ESP_O11z_S_aa += SQ_ESP_O_S_aa_coefs*I_ESP_O11z_S_vrr;

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
     * shell quartet name: SQ_ESP_M_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 3 integrals are omitted 
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
    Double I_ESP_M9y_Px_aa = I_ESP_Nx9y_S_aa+ABX*I_ESP_M9y_S_aa;
    Double I_ESP_M8yz_Px_aa = I_ESP_Nx8yz_S_aa+ABX*I_ESP_M8yz_S_aa;
    Double I_ESP_M7y2z_Px_aa = I_ESP_Nx7y2z_S_aa+ABX*I_ESP_M7y2z_S_aa;
    Double I_ESP_M6y3z_Px_aa = I_ESP_Nx6y3z_S_aa+ABX*I_ESP_M6y3z_S_aa;
    Double I_ESP_M5y4z_Px_aa = I_ESP_Nx5y4z_S_aa+ABX*I_ESP_M5y4z_S_aa;
    Double I_ESP_M4y5z_Px_aa = I_ESP_Nx4y5z_S_aa+ABX*I_ESP_M4y5z_S_aa;
    Double I_ESP_M3y6z_Px_aa = I_ESP_Nx3y6z_S_aa+ABX*I_ESP_M3y6z_S_aa;
    Double I_ESP_M2y7z_Px_aa = I_ESP_Nx2y7z_S_aa+ABX*I_ESP_M2y7z_S_aa;
    Double I_ESP_My8z_Px_aa = I_ESP_Nxy8z_S_aa+ABX*I_ESP_My8z_S_aa;
    Double I_ESP_M9z_Px_aa = I_ESP_Nx9z_S_aa+ABX*I_ESP_M9z_S_aa;
    Double I_ESP_M8xy_Py_aa = I_ESP_N8x2y_S_aa+ABY*I_ESP_M8xy_S_aa;
    Double I_ESP_M8xz_Py_aa = I_ESP_N8xyz_S_aa+ABY*I_ESP_M8xz_S_aa;
    Double I_ESP_M7x2y_Py_aa = I_ESP_N7x3y_S_aa+ABY*I_ESP_M7x2y_S_aa;
    Double I_ESP_M7xyz_Py_aa = I_ESP_N7x2yz_S_aa+ABY*I_ESP_M7xyz_S_aa;
    Double I_ESP_M7x2z_Py_aa = I_ESP_N7xy2z_S_aa+ABY*I_ESP_M7x2z_S_aa;
    Double I_ESP_M6x3y_Py_aa = I_ESP_N6x4y_S_aa+ABY*I_ESP_M6x3y_S_aa;
    Double I_ESP_M6x2yz_Py_aa = I_ESP_N6x3yz_S_aa+ABY*I_ESP_M6x2yz_S_aa;
    Double I_ESP_M6xy2z_Py_aa = I_ESP_N6x2y2z_S_aa+ABY*I_ESP_M6xy2z_S_aa;
    Double I_ESP_M6x3z_Py_aa = I_ESP_N6xy3z_S_aa+ABY*I_ESP_M6x3z_S_aa;
    Double I_ESP_M5x4y_Py_aa = I_ESP_N5x5y_S_aa+ABY*I_ESP_M5x4y_S_aa;
    Double I_ESP_M5x3yz_Py_aa = I_ESP_N5x4yz_S_aa+ABY*I_ESP_M5x3yz_S_aa;
    Double I_ESP_M5x2y2z_Py_aa = I_ESP_N5x3y2z_S_aa+ABY*I_ESP_M5x2y2z_S_aa;
    Double I_ESP_M5xy3z_Py_aa = I_ESP_N5x2y3z_S_aa+ABY*I_ESP_M5xy3z_S_aa;
    Double I_ESP_M5x4z_Py_aa = I_ESP_N5xy4z_S_aa+ABY*I_ESP_M5x4z_S_aa;
    Double I_ESP_M4x5y_Py_aa = I_ESP_N4x6y_S_aa+ABY*I_ESP_M4x5y_S_aa;
    Double I_ESP_M4x4yz_Py_aa = I_ESP_N4x5yz_S_aa+ABY*I_ESP_M4x4yz_S_aa;
    Double I_ESP_M4x3y2z_Py_aa = I_ESP_N4x4y2z_S_aa+ABY*I_ESP_M4x3y2z_S_aa;
    Double I_ESP_M4x2y3z_Py_aa = I_ESP_N4x3y3z_S_aa+ABY*I_ESP_M4x2y3z_S_aa;
    Double I_ESP_M4xy4z_Py_aa = I_ESP_N4x2y4z_S_aa+ABY*I_ESP_M4xy4z_S_aa;
    Double I_ESP_M4x5z_Py_aa = I_ESP_N4xy5z_S_aa+ABY*I_ESP_M4x5z_S_aa;
    Double I_ESP_M3x6y_Py_aa = I_ESP_N3x7y_S_aa+ABY*I_ESP_M3x6y_S_aa;
    Double I_ESP_M3x5yz_Py_aa = I_ESP_N3x6yz_S_aa+ABY*I_ESP_M3x5yz_S_aa;
    Double I_ESP_M3x4y2z_Py_aa = I_ESP_N3x5y2z_S_aa+ABY*I_ESP_M3x4y2z_S_aa;
    Double I_ESP_M3x3y3z_Py_aa = I_ESP_N3x4y3z_S_aa+ABY*I_ESP_M3x3y3z_S_aa;
    Double I_ESP_M3x2y4z_Py_aa = I_ESP_N3x3y4z_S_aa+ABY*I_ESP_M3x2y4z_S_aa;
    Double I_ESP_M3xy5z_Py_aa = I_ESP_N3x2y5z_S_aa+ABY*I_ESP_M3xy5z_S_aa;
    Double I_ESP_M3x6z_Py_aa = I_ESP_N3xy6z_S_aa+ABY*I_ESP_M3x6z_S_aa;
    Double I_ESP_M2x7y_Py_aa = I_ESP_N2x8y_S_aa+ABY*I_ESP_M2x7y_S_aa;
    Double I_ESP_M2x6yz_Py_aa = I_ESP_N2x7yz_S_aa+ABY*I_ESP_M2x6yz_S_aa;
    Double I_ESP_M2x5y2z_Py_aa = I_ESP_N2x6y2z_S_aa+ABY*I_ESP_M2x5y2z_S_aa;
    Double I_ESP_M2x4y3z_Py_aa = I_ESP_N2x5y3z_S_aa+ABY*I_ESP_M2x4y3z_S_aa;
    Double I_ESP_M2x3y4z_Py_aa = I_ESP_N2x4y4z_S_aa+ABY*I_ESP_M2x3y4z_S_aa;
    Double I_ESP_M2x2y5z_Py_aa = I_ESP_N2x3y5z_S_aa+ABY*I_ESP_M2x2y5z_S_aa;
    Double I_ESP_M2xy6z_Py_aa = I_ESP_N2x2y6z_S_aa+ABY*I_ESP_M2xy6z_S_aa;
    Double I_ESP_M2x7z_Py_aa = I_ESP_N2xy7z_S_aa+ABY*I_ESP_M2x7z_S_aa;
    Double I_ESP_Mx8y_Py_aa = I_ESP_Nx9y_S_aa+ABY*I_ESP_Mx8y_S_aa;
    Double I_ESP_Mx7yz_Py_aa = I_ESP_Nx8yz_S_aa+ABY*I_ESP_Mx7yz_S_aa;
    Double I_ESP_Mx6y2z_Py_aa = I_ESP_Nx7y2z_S_aa+ABY*I_ESP_Mx6y2z_S_aa;
    Double I_ESP_Mx5y3z_Py_aa = I_ESP_Nx6y3z_S_aa+ABY*I_ESP_Mx5y3z_S_aa;
    Double I_ESP_Mx4y4z_Py_aa = I_ESP_Nx5y4z_S_aa+ABY*I_ESP_Mx4y4z_S_aa;
    Double I_ESP_Mx3y5z_Py_aa = I_ESP_Nx4y5z_S_aa+ABY*I_ESP_Mx3y5z_S_aa;
    Double I_ESP_Mx2y6z_Py_aa = I_ESP_Nx3y6z_S_aa+ABY*I_ESP_Mx2y6z_S_aa;
    Double I_ESP_Mxy7z_Py_aa = I_ESP_Nx2y7z_S_aa+ABY*I_ESP_Mxy7z_S_aa;
    Double I_ESP_Mx8z_Py_aa = I_ESP_Nxy8z_S_aa+ABY*I_ESP_Mx8z_S_aa;
    Double I_ESP_M9y_Py_aa = I_ESP_N10y_S_aa+ABY*I_ESP_M9y_S_aa;
    Double I_ESP_M8yz_Py_aa = I_ESP_N9yz_S_aa+ABY*I_ESP_M8yz_S_aa;
    Double I_ESP_M7y2z_Py_aa = I_ESP_N8y2z_S_aa+ABY*I_ESP_M7y2z_S_aa;
    Double I_ESP_M6y3z_Py_aa = I_ESP_N7y3z_S_aa+ABY*I_ESP_M6y3z_S_aa;
    Double I_ESP_M5y4z_Py_aa = I_ESP_N6y4z_S_aa+ABY*I_ESP_M5y4z_S_aa;
    Double I_ESP_M4y5z_Py_aa = I_ESP_N5y5z_S_aa+ABY*I_ESP_M4y5z_S_aa;
    Double I_ESP_M3y6z_Py_aa = I_ESP_N4y6z_S_aa+ABY*I_ESP_M3y6z_S_aa;
    Double I_ESP_M2y7z_Py_aa = I_ESP_N3y7z_S_aa+ABY*I_ESP_M2y7z_S_aa;
    Double I_ESP_My8z_Py_aa = I_ESP_N2y8z_S_aa+ABY*I_ESP_My8z_S_aa;
    Double I_ESP_M9z_Py_aa = I_ESP_Ny9z_S_aa+ABY*I_ESP_M9z_S_aa;
    Double I_ESP_M8xy_Pz_aa = I_ESP_N8xyz_S_aa+ABZ*I_ESP_M8xy_S_aa;
    Double I_ESP_M8xz_Pz_aa = I_ESP_N8x2z_S_aa+ABZ*I_ESP_M8xz_S_aa;
    Double I_ESP_M7x2y_Pz_aa = I_ESP_N7x2yz_S_aa+ABZ*I_ESP_M7x2y_S_aa;
    Double I_ESP_M7xyz_Pz_aa = I_ESP_N7xy2z_S_aa+ABZ*I_ESP_M7xyz_S_aa;
    Double I_ESP_M7x2z_Pz_aa = I_ESP_N7x3z_S_aa+ABZ*I_ESP_M7x2z_S_aa;
    Double I_ESP_M6x3y_Pz_aa = I_ESP_N6x3yz_S_aa+ABZ*I_ESP_M6x3y_S_aa;
    Double I_ESP_M6x2yz_Pz_aa = I_ESP_N6x2y2z_S_aa+ABZ*I_ESP_M6x2yz_S_aa;
    Double I_ESP_M6xy2z_Pz_aa = I_ESP_N6xy3z_S_aa+ABZ*I_ESP_M6xy2z_S_aa;
    Double I_ESP_M6x3z_Pz_aa = I_ESP_N6x4z_S_aa+ABZ*I_ESP_M6x3z_S_aa;
    Double I_ESP_M5x4y_Pz_aa = I_ESP_N5x4yz_S_aa+ABZ*I_ESP_M5x4y_S_aa;
    Double I_ESP_M5x3yz_Pz_aa = I_ESP_N5x3y2z_S_aa+ABZ*I_ESP_M5x3yz_S_aa;
    Double I_ESP_M5x2y2z_Pz_aa = I_ESP_N5x2y3z_S_aa+ABZ*I_ESP_M5x2y2z_S_aa;
    Double I_ESP_M5xy3z_Pz_aa = I_ESP_N5xy4z_S_aa+ABZ*I_ESP_M5xy3z_S_aa;
    Double I_ESP_M5x4z_Pz_aa = I_ESP_N5x5z_S_aa+ABZ*I_ESP_M5x4z_S_aa;
    Double I_ESP_M4x5y_Pz_aa = I_ESP_N4x5yz_S_aa+ABZ*I_ESP_M4x5y_S_aa;
    Double I_ESP_M4x4yz_Pz_aa = I_ESP_N4x4y2z_S_aa+ABZ*I_ESP_M4x4yz_S_aa;
    Double I_ESP_M4x3y2z_Pz_aa = I_ESP_N4x3y3z_S_aa+ABZ*I_ESP_M4x3y2z_S_aa;
    Double I_ESP_M4x2y3z_Pz_aa = I_ESP_N4x2y4z_S_aa+ABZ*I_ESP_M4x2y3z_S_aa;
    Double I_ESP_M4xy4z_Pz_aa = I_ESP_N4xy5z_S_aa+ABZ*I_ESP_M4xy4z_S_aa;
    Double I_ESP_M4x5z_Pz_aa = I_ESP_N4x6z_S_aa+ABZ*I_ESP_M4x5z_S_aa;
    Double I_ESP_M3x6y_Pz_aa = I_ESP_N3x6yz_S_aa+ABZ*I_ESP_M3x6y_S_aa;
    Double I_ESP_M3x5yz_Pz_aa = I_ESP_N3x5y2z_S_aa+ABZ*I_ESP_M3x5yz_S_aa;
    Double I_ESP_M3x4y2z_Pz_aa = I_ESP_N3x4y3z_S_aa+ABZ*I_ESP_M3x4y2z_S_aa;
    Double I_ESP_M3x3y3z_Pz_aa = I_ESP_N3x3y4z_S_aa+ABZ*I_ESP_M3x3y3z_S_aa;
    Double I_ESP_M3x2y4z_Pz_aa = I_ESP_N3x2y5z_S_aa+ABZ*I_ESP_M3x2y4z_S_aa;
    Double I_ESP_M3xy5z_Pz_aa = I_ESP_N3xy6z_S_aa+ABZ*I_ESP_M3xy5z_S_aa;
    Double I_ESP_M3x6z_Pz_aa = I_ESP_N3x7z_S_aa+ABZ*I_ESP_M3x6z_S_aa;
    Double I_ESP_M2x7y_Pz_aa = I_ESP_N2x7yz_S_aa+ABZ*I_ESP_M2x7y_S_aa;
    Double I_ESP_M2x6yz_Pz_aa = I_ESP_N2x6y2z_S_aa+ABZ*I_ESP_M2x6yz_S_aa;
    Double I_ESP_M2x5y2z_Pz_aa = I_ESP_N2x5y3z_S_aa+ABZ*I_ESP_M2x5y2z_S_aa;
    Double I_ESP_M2x4y3z_Pz_aa = I_ESP_N2x4y4z_S_aa+ABZ*I_ESP_M2x4y3z_S_aa;
    Double I_ESP_M2x3y4z_Pz_aa = I_ESP_N2x3y5z_S_aa+ABZ*I_ESP_M2x3y4z_S_aa;
    Double I_ESP_M2x2y5z_Pz_aa = I_ESP_N2x2y6z_S_aa+ABZ*I_ESP_M2x2y5z_S_aa;
    Double I_ESP_M2xy6z_Pz_aa = I_ESP_N2xy7z_S_aa+ABZ*I_ESP_M2xy6z_S_aa;
    Double I_ESP_M2x7z_Pz_aa = I_ESP_N2x8z_S_aa+ABZ*I_ESP_M2x7z_S_aa;
    Double I_ESP_Mx8y_Pz_aa = I_ESP_Nx8yz_S_aa+ABZ*I_ESP_Mx8y_S_aa;
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
     * totally 135 integrals are omitted 
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
     * totally 72 integrals are omitted 
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
     * shell quartet name: SQ_ESP_N_P_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 36 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_O_S_aa
     * RHS shell quartet name: SQ_ESP_N_S_aa
     ************************************************************/
    Double I_ESP_N10x_Px_aa = I_ESP_O11x_S_aa+ABX*I_ESP_N10x_S_aa;
    Double I_ESP_N9xy_Px_aa = I_ESP_O10xy_S_aa+ABX*I_ESP_N9xy_S_aa;
    Double I_ESP_N9xz_Px_aa = I_ESP_O10xz_S_aa+ABX*I_ESP_N9xz_S_aa;
    Double I_ESP_N8x2y_Px_aa = I_ESP_O9x2y_S_aa+ABX*I_ESP_N8x2y_S_aa;
    Double I_ESP_N8xyz_Px_aa = I_ESP_O9xyz_S_aa+ABX*I_ESP_N8xyz_S_aa;
    Double I_ESP_N8x2z_Px_aa = I_ESP_O9x2z_S_aa+ABX*I_ESP_N8x2z_S_aa;
    Double I_ESP_N7x3y_Px_aa = I_ESP_O8x3y_S_aa+ABX*I_ESP_N7x3y_S_aa;
    Double I_ESP_N7x2yz_Px_aa = I_ESP_O8x2yz_S_aa+ABX*I_ESP_N7x2yz_S_aa;
    Double I_ESP_N7xy2z_Px_aa = I_ESP_O8xy2z_S_aa+ABX*I_ESP_N7xy2z_S_aa;
    Double I_ESP_N7x3z_Px_aa = I_ESP_O8x3z_S_aa+ABX*I_ESP_N7x3z_S_aa;
    Double I_ESP_N6x4y_Px_aa = I_ESP_O7x4y_S_aa+ABX*I_ESP_N6x4y_S_aa;
    Double I_ESP_N6x3yz_Px_aa = I_ESP_O7x3yz_S_aa+ABX*I_ESP_N6x3yz_S_aa;
    Double I_ESP_N6x2y2z_Px_aa = I_ESP_O7x2y2z_S_aa+ABX*I_ESP_N6x2y2z_S_aa;
    Double I_ESP_N6xy3z_Px_aa = I_ESP_O7xy3z_S_aa+ABX*I_ESP_N6xy3z_S_aa;
    Double I_ESP_N6x4z_Px_aa = I_ESP_O7x4z_S_aa+ABX*I_ESP_N6x4z_S_aa;
    Double I_ESP_N5x5y_Px_aa = I_ESP_O6x5y_S_aa+ABX*I_ESP_N5x5y_S_aa;
    Double I_ESP_N5x4yz_Px_aa = I_ESP_O6x4yz_S_aa+ABX*I_ESP_N5x4yz_S_aa;
    Double I_ESP_N5x3y2z_Px_aa = I_ESP_O6x3y2z_S_aa+ABX*I_ESP_N5x3y2z_S_aa;
    Double I_ESP_N5x2y3z_Px_aa = I_ESP_O6x2y3z_S_aa+ABX*I_ESP_N5x2y3z_S_aa;
    Double I_ESP_N5xy4z_Px_aa = I_ESP_O6xy4z_S_aa+ABX*I_ESP_N5xy4z_S_aa;
    Double I_ESP_N5x5z_Px_aa = I_ESP_O6x5z_S_aa+ABX*I_ESP_N5x5z_S_aa;
    Double I_ESP_N4x6y_Px_aa = I_ESP_O5x6y_S_aa+ABX*I_ESP_N4x6y_S_aa;
    Double I_ESP_N4x5yz_Px_aa = I_ESP_O5x5yz_S_aa+ABX*I_ESP_N4x5yz_S_aa;
    Double I_ESP_N4x4y2z_Px_aa = I_ESP_O5x4y2z_S_aa+ABX*I_ESP_N4x4y2z_S_aa;
    Double I_ESP_N4x3y3z_Px_aa = I_ESP_O5x3y3z_S_aa+ABX*I_ESP_N4x3y3z_S_aa;
    Double I_ESP_N4x2y4z_Px_aa = I_ESP_O5x2y4z_S_aa+ABX*I_ESP_N4x2y4z_S_aa;
    Double I_ESP_N4xy5z_Px_aa = I_ESP_O5xy5z_S_aa+ABX*I_ESP_N4xy5z_S_aa;
    Double I_ESP_N4x6z_Px_aa = I_ESP_O5x6z_S_aa+ABX*I_ESP_N4x6z_S_aa;
    Double I_ESP_N3x7y_Px_aa = I_ESP_O4x7y_S_aa+ABX*I_ESP_N3x7y_S_aa;
    Double I_ESP_N3x6yz_Px_aa = I_ESP_O4x6yz_S_aa+ABX*I_ESP_N3x6yz_S_aa;
    Double I_ESP_N3x5y2z_Px_aa = I_ESP_O4x5y2z_S_aa+ABX*I_ESP_N3x5y2z_S_aa;
    Double I_ESP_N3x4y3z_Px_aa = I_ESP_O4x4y3z_S_aa+ABX*I_ESP_N3x4y3z_S_aa;
    Double I_ESP_N3x3y4z_Px_aa = I_ESP_O4x3y4z_S_aa+ABX*I_ESP_N3x3y4z_S_aa;
    Double I_ESP_N3x2y5z_Px_aa = I_ESP_O4x2y5z_S_aa+ABX*I_ESP_N3x2y5z_S_aa;
    Double I_ESP_N3xy6z_Px_aa = I_ESP_O4xy6z_S_aa+ABX*I_ESP_N3xy6z_S_aa;
    Double I_ESP_N3x7z_Px_aa = I_ESP_O4x7z_S_aa+ABX*I_ESP_N3x7z_S_aa;
    Double I_ESP_N2x8y_Px_aa = I_ESP_O3x8y_S_aa+ABX*I_ESP_N2x8y_S_aa;
    Double I_ESP_N2x7yz_Px_aa = I_ESP_O3x7yz_S_aa+ABX*I_ESP_N2x7yz_S_aa;
    Double I_ESP_N2x6y2z_Px_aa = I_ESP_O3x6y2z_S_aa+ABX*I_ESP_N2x6y2z_S_aa;
    Double I_ESP_N2x5y3z_Px_aa = I_ESP_O3x5y3z_S_aa+ABX*I_ESP_N2x5y3z_S_aa;
    Double I_ESP_N2x4y4z_Px_aa = I_ESP_O3x4y4z_S_aa+ABX*I_ESP_N2x4y4z_S_aa;
    Double I_ESP_N2x3y5z_Px_aa = I_ESP_O3x3y5z_S_aa+ABX*I_ESP_N2x3y5z_S_aa;
    Double I_ESP_N2x2y6z_Px_aa = I_ESP_O3x2y6z_S_aa+ABX*I_ESP_N2x2y6z_S_aa;
    Double I_ESP_N2xy7z_Px_aa = I_ESP_O3xy7z_S_aa+ABX*I_ESP_N2xy7z_S_aa;
    Double I_ESP_N2x8z_Px_aa = I_ESP_O3x8z_S_aa+ABX*I_ESP_N2x8z_S_aa;
    Double I_ESP_Nx9y_Px_aa = I_ESP_O2x9y_S_aa+ABX*I_ESP_Nx9y_S_aa;
    Double I_ESP_Nx8yz_Px_aa = I_ESP_O2x8yz_S_aa+ABX*I_ESP_Nx8yz_S_aa;
    Double I_ESP_Nx7y2z_Px_aa = I_ESP_O2x7y2z_S_aa+ABX*I_ESP_Nx7y2z_S_aa;
    Double I_ESP_Nx6y3z_Px_aa = I_ESP_O2x6y3z_S_aa+ABX*I_ESP_Nx6y3z_S_aa;
    Double I_ESP_Nx5y4z_Px_aa = I_ESP_O2x5y4z_S_aa+ABX*I_ESP_Nx5y4z_S_aa;
    Double I_ESP_Nx4y5z_Px_aa = I_ESP_O2x4y5z_S_aa+ABX*I_ESP_Nx4y5z_S_aa;
    Double I_ESP_Nx3y6z_Px_aa = I_ESP_O2x3y6z_S_aa+ABX*I_ESP_Nx3y6z_S_aa;
    Double I_ESP_Nx2y7z_Px_aa = I_ESP_O2x2y7z_S_aa+ABX*I_ESP_Nx2y7z_S_aa;
    Double I_ESP_Nxy8z_Px_aa = I_ESP_O2xy8z_S_aa+ABX*I_ESP_Nxy8z_S_aa;
    Double I_ESP_Nx9z_Px_aa = I_ESP_O2x9z_S_aa+ABX*I_ESP_Nx9z_S_aa;
    Double I_ESP_N8x2y_Py_aa = I_ESP_O8x3y_S_aa+ABY*I_ESP_N8x2y_S_aa;
    Double I_ESP_N8xyz_Py_aa = I_ESP_O8x2yz_S_aa+ABY*I_ESP_N8xyz_S_aa;
    Double I_ESP_N7x3y_Py_aa = I_ESP_O7x4y_S_aa+ABY*I_ESP_N7x3y_S_aa;
    Double I_ESP_N7x2yz_Py_aa = I_ESP_O7x3yz_S_aa+ABY*I_ESP_N7x2yz_S_aa;
    Double I_ESP_N7xy2z_Py_aa = I_ESP_O7x2y2z_S_aa+ABY*I_ESP_N7xy2z_S_aa;
    Double I_ESP_N6x4y_Py_aa = I_ESP_O6x5y_S_aa+ABY*I_ESP_N6x4y_S_aa;
    Double I_ESP_N6x3yz_Py_aa = I_ESP_O6x4yz_S_aa+ABY*I_ESP_N6x3yz_S_aa;
    Double I_ESP_N6x2y2z_Py_aa = I_ESP_O6x3y2z_S_aa+ABY*I_ESP_N6x2y2z_S_aa;
    Double I_ESP_N6xy3z_Py_aa = I_ESP_O6x2y3z_S_aa+ABY*I_ESP_N6xy3z_S_aa;
    Double I_ESP_N5x5y_Py_aa = I_ESP_O5x6y_S_aa+ABY*I_ESP_N5x5y_S_aa;
    Double I_ESP_N5x4yz_Py_aa = I_ESP_O5x5yz_S_aa+ABY*I_ESP_N5x4yz_S_aa;
    Double I_ESP_N5x3y2z_Py_aa = I_ESP_O5x4y2z_S_aa+ABY*I_ESP_N5x3y2z_S_aa;
    Double I_ESP_N5x2y3z_Py_aa = I_ESP_O5x3y3z_S_aa+ABY*I_ESP_N5x2y3z_S_aa;
    Double I_ESP_N5xy4z_Py_aa = I_ESP_O5x2y4z_S_aa+ABY*I_ESP_N5xy4z_S_aa;
    Double I_ESP_N4x6y_Py_aa = I_ESP_O4x7y_S_aa+ABY*I_ESP_N4x6y_S_aa;
    Double I_ESP_N4x5yz_Py_aa = I_ESP_O4x6yz_S_aa+ABY*I_ESP_N4x5yz_S_aa;
    Double I_ESP_N4x4y2z_Py_aa = I_ESP_O4x5y2z_S_aa+ABY*I_ESP_N4x4y2z_S_aa;
    Double I_ESP_N4x3y3z_Py_aa = I_ESP_O4x4y3z_S_aa+ABY*I_ESP_N4x3y3z_S_aa;
    Double I_ESP_N4x2y4z_Py_aa = I_ESP_O4x3y4z_S_aa+ABY*I_ESP_N4x2y4z_S_aa;
    Double I_ESP_N4xy5z_Py_aa = I_ESP_O4x2y5z_S_aa+ABY*I_ESP_N4xy5z_S_aa;
    Double I_ESP_N3x7y_Py_aa = I_ESP_O3x8y_S_aa+ABY*I_ESP_N3x7y_S_aa;
    Double I_ESP_N3x6yz_Py_aa = I_ESP_O3x7yz_S_aa+ABY*I_ESP_N3x6yz_S_aa;
    Double I_ESP_N3x5y2z_Py_aa = I_ESP_O3x6y2z_S_aa+ABY*I_ESP_N3x5y2z_S_aa;
    Double I_ESP_N3x4y3z_Py_aa = I_ESP_O3x5y3z_S_aa+ABY*I_ESP_N3x4y3z_S_aa;
    Double I_ESP_N3x3y4z_Py_aa = I_ESP_O3x4y4z_S_aa+ABY*I_ESP_N3x3y4z_S_aa;
    Double I_ESP_N3x2y5z_Py_aa = I_ESP_O3x3y5z_S_aa+ABY*I_ESP_N3x2y5z_S_aa;
    Double I_ESP_N3xy6z_Py_aa = I_ESP_O3x2y6z_S_aa+ABY*I_ESP_N3xy6z_S_aa;
    Double I_ESP_N2x8y_Py_aa = I_ESP_O2x9y_S_aa+ABY*I_ESP_N2x8y_S_aa;
    Double I_ESP_N2x7yz_Py_aa = I_ESP_O2x8yz_S_aa+ABY*I_ESP_N2x7yz_S_aa;
    Double I_ESP_N2x6y2z_Py_aa = I_ESP_O2x7y2z_S_aa+ABY*I_ESP_N2x6y2z_S_aa;
    Double I_ESP_N2x5y3z_Py_aa = I_ESP_O2x6y3z_S_aa+ABY*I_ESP_N2x5y3z_S_aa;
    Double I_ESP_N2x4y4z_Py_aa = I_ESP_O2x5y4z_S_aa+ABY*I_ESP_N2x4y4z_S_aa;
    Double I_ESP_N2x3y5z_Py_aa = I_ESP_O2x4y5z_S_aa+ABY*I_ESP_N2x3y5z_S_aa;
    Double I_ESP_N2x2y6z_Py_aa = I_ESP_O2x3y6z_S_aa+ABY*I_ESP_N2x2y6z_S_aa;
    Double I_ESP_N2xy7z_Py_aa = I_ESP_O2x2y7z_S_aa+ABY*I_ESP_N2xy7z_S_aa;
    Double I_ESP_Nx9y_Py_aa = I_ESP_Ox10y_S_aa+ABY*I_ESP_Nx9y_S_aa;
    Double I_ESP_Nx8yz_Py_aa = I_ESP_Ox9yz_S_aa+ABY*I_ESP_Nx8yz_S_aa;
    Double I_ESP_Nx7y2z_Py_aa = I_ESP_Ox8y2z_S_aa+ABY*I_ESP_Nx7y2z_S_aa;
    Double I_ESP_Nx6y3z_Py_aa = I_ESP_Ox7y3z_S_aa+ABY*I_ESP_Nx6y3z_S_aa;
    Double I_ESP_Nx5y4z_Py_aa = I_ESP_Ox6y4z_S_aa+ABY*I_ESP_Nx5y4z_S_aa;
    Double I_ESP_Nx4y5z_Py_aa = I_ESP_Ox5y5z_S_aa+ABY*I_ESP_Nx4y5z_S_aa;
    Double I_ESP_Nx3y6z_Py_aa = I_ESP_Ox4y6z_S_aa+ABY*I_ESP_Nx3y6z_S_aa;
    Double I_ESP_Nx2y7z_Py_aa = I_ESP_Ox3y7z_S_aa+ABY*I_ESP_Nx2y7z_S_aa;
    Double I_ESP_Nxy8z_Py_aa = I_ESP_Ox2y8z_S_aa+ABY*I_ESP_Nxy8z_S_aa;
    Double I_ESP_N10y_Py_aa = I_ESP_O11y_S_aa+ABY*I_ESP_N10y_S_aa;
    Double I_ESP_N9yz_Py_aa = I_ESP_O10yz_S_aa+ABY*I_ESP_N9yz_S_aa;
    Double I_ESP_N8y2z_Py_aa = I_ESP_O9y2z_S_aa+ABY*I_ESP_N8y2z_S_aa;
    Double I_ESP_N7y3z_Py_aa = I_ESP_O8y3z_S_aa+ABY*I_ESP_N7y3z_S_aa;
    Double I_ESP_N6y4z_Py_aa = I_ESP_O7y4z_S_aa+ABY*I_ESP_N6y4z_S_aa;
    Double I_ESP_N5y5z_Py_aa = I_ESP_O6y5z_S_aa+ABY*I_ESP_N5y5z_S_aa;
    Double I_ESP_N4y6z_Py_aa = I_ESP_O5y6z_S_aa+ABY*I_ESP_N4y6z_S_aa;
    Double I_ESP_N3y7z_Py_aa = I_ESP_O4y7z_S_aa+ABY*I_ESP_N3y7z_S_aa;
    Double I_ESP_N2y8z_Py_aa = I_ESP_O3y8z_S_aa+ABY*I_ESP_N2y8z_S_aa;
    Double I_ESP_Ny9z_Py_aa = I_ESP_O2y9z_S_aa+ABY*I_ESP_Ny9z_S_aa;
    Double I_ESP_N8xyz_Pz_aa = I_ESP_O8xy2z_S_aa+ABZ*I_ESP_N8xyz_S_aa;
    Double I_ESP_N8x2z_Pz_aa = I_ESP_O8x3z_S_aa+ABZ*I_ESP_N8x2z_S_aa;
    Double I_ESP_N7x2yz_Pz_aa = I_ESP_O7x2y2z_S_aa+ABZ*I_ESP_N7x2yz_S_aa;
    Double I_ESP_N7xy2z_Pz_aa = I_ESP_O7xy3z_S_aa+ABZ*I_ESP_N7xy2z_S_aa;
    Double I_ESP_N7x3z_Pz_aa = I_ESP_O7x4z_S_aa+ABZ*I_ESP_N7x3z_S_aa;
    Double I_ESP_N6x3yz_Pz_aa = I_ESP_O6x3y2z_S_aa+ABZ*I_ESP_N6x3yz_S_aa;
    Double I_ESP_N6x2y2z_Pz_aa = I_ESP_O6x2y3z_S_aa+ABZ*I_ESP_N6x2y2z_S_aa;
    Double I_ESP_N6xy3z_Pz_aa = I_ESP_O6xy4z_S_aa+ABZ*I_ESP_N6xy3z_S_aa;
    Double I_ESP_N6x4z_Pz_aa = I_ESP_O6x5z_S_aa+ABZ*I_ESP_N6x4z_S_aa;
    Double I_ESP_N5x4yz_Pz_aa = I_ESP_O5x4y2z_S_aa+ABZ*I_ESP_N5x4yz_S_aa;
    Double I_ESP_N5x3y2z_Pz_aa = I_ESP_O5x3y3z_S_aa+ABZ*I_ESP_N5x3y2z_S_aa;
    Double I_ESP_N5x2y3z_Pz_aa = I_ESP_O5x2y4z_S_aa+ABZ*I_ESP_N5x2y3z_S_aa;
    Double I_ESP_N5xy4z_Pz_aa = I_ESP_O5xy5z_S_aa+ABZ*I_ESP_N5xy4z_S_aa;
    Double I_ESP_N5x5z_Pz_aa = I_ESP_O5x6z_S_aa+ABZ*I_ESP_N5x5z_S_aa;
    Double I_ESP_N4x5yz_Pz_aa = I_ESP_O4x5y2z_S_aa+ABZ*I_ESP_N4x5yz_S_aa;
    Double I_ESP_N4x4y2z_Pz_aa = I_ESP_O4x4y3z_S_aa+ABZ*I_ESP_N4x4y2z_S_aa;
    Double I_ESP_N4x3y3z_Pz_aa = I_ESP_O4x3y4z_S_aa+ABZ*I_ESP_N4x3y3z_S_aa;
    Double I_ESP_N4x2y4z_Pz_aa = I_ESP_O4x2y5z_S_aa+ABZ*I_ESP_N4x2y4z_S_aa;
    Double I_ESP_N4xy5z_Pz_aa = I_ESP_O4xy6z_S_aa+ABZ*I_ESP_N4xy5z_S_aa;
    Double I_ESP_N4x6z_Pz_aa = I_ESP_O4x7z_S_aa+ABZ*I_ESP_N4x6z_S_aa;
    Double I_ESP_N3x6yz_Pz_aa = I_ESP_O3x6y2z_S_aa+ABZ*I_ESP_N3x6yz_S_aa;
    Double I_ESP_N3x5y2z_Pz_aa = I_ESP_O3x5y3z_S_aa+ABZ*I_ESP_N3x5y2z_S_aa;
    Double I_ESP_N3x4y3z_Pz_aa = I_ESP_O3x4y4z_S_aa+ABZ*I_ESP_N3x4y3z_S_aa;
    Double I_ESP_N3x3y4z_Pz_aa = I_ESP_O3x3y5z_S_aa+ABZ*I_ESP_N3x3y4z_S_aa;
    Double I_ESP_N3x2y5z_Pz_aa = I_ESP_O3x2y6z_S_aa+ABZ*I_ESP_N3x2y5z_S_aa;
    Double I_ESP_N3xy6z_Pz_aa = I_ESP_O3xy7z_S_aa+ABZ*I_ESP_N3xy6z_S_aa;
    Double I_ESP_N3x7z_Pz_aa = I_ESP_O3x8z_S_aa+ABZ*I_ESP_N3x7z_S_aa;
    Double I_ESP_N2x7yz_Pz_aa = I_ESP_O2x7y2z_S_aa+ABZ*I_ESP_N2x7yz_S_aa;
    Double I_ESP_N2x6y2z_Pz_aa = I_ESP_O2x6y3z_S_aa+ABZ*I_ESP_N2x6y2z_S_aa;
    Double I_ESP_N2x5y3z_Pz_aa = I_ESP_O2x5y4z_S_aa+ABZ*I_ESP_N2x5y3z_S_aa;
    Double I_ESP_N2x4y4z_Pz_aa = I_ESP_O2x4y5z_S_aa+ABZ*I_ESP_N2x4y4z_S_aa;
    Double I_ESP_N2x3y5z_Pz_aa = I_ESP_O2x3y6z_S_aa+ABZ*I_ESP_N2x3y5z_S_aa;
    Double I_ESP_N2x2y6z_Pz_aa = I_ESP_O2x2y7z_S_aa+ABZ*I_ESP_N2x2y6z_S_aa;
    Double I_ESP_N2xy7z_Pz_aa = I_ESP_O2xy8z_S_aa+ABZ*I_ESP_N2xy7z_S_aa;
    Double I_ESP_N2x8z_Pz_aa = I_ESP_O2x9z_S_aa+ABZ*I_ESP_N2x8z_S_aa;
    Double I_ESP_Nx8yz_Pz_aa = I_ESP_Ox8y2z_S_aa+ABZ*I_ESP_Nx8yz_S_aa;
    Double I_ESP_Nx7y2z_Pz_aa = I_ESP_Ox7y3z_S_aa+ABZ*I_ESP_Nx7y2z_S_aa;
    Double I_ESP_Nx6y3z_Pz_aa = I_ESP_Ox6y4z_S_aa+ABZ*I_ESP_Nx6y3z_S_aa;
    Double I_ESP_Nx5y4z_Pz_aa = I_ESP_Ox5y5z_S_aa+ABZ*I_ESP_Nx5y4z_S_aa;
    Double I_ESP_Nx4y5z_Pz_aa = I_ESP_Ox4y6z_S_aa+ABZ*I_ESP_Nx4y5z_S_aa;
    Double I_ESP_Nx3y6z_Pz_aa = I_ESP_Ox3y7z_S_aa+ABZ*I_ESP_Nx3y6z_S_aa;
    Double I_ESP_Nx2y7z_Pz_aa = I_ESP_Ox2y8z_S_aa+ABZ*I_ESP_Nx2y7z_S_aa;
    Double I_ESP_Nxy8z_Pz_aa = I_ESP_Oxy9z_S_aa+ABZ*I_ESP_Nxy8z_S_aa;
    Double I_ESP_Nx9z_Pz_aa = I_ESP_Ox10z_S_aa+ABZ*I_ESP_Nx9z_S_aa;
    Double I_ESP_N8y2z_Pz_aa = I_ESP_O8y3z_S_aa+ABZ*I_ESP_N8y2z_S_aa;
    Double I_ESP_N7y3z_Pz_aa = I_ESP_O7y4z_S_aa+ABZ*I_ESP_N7y3z_S_aa;
    Double I_ESP_N6y4z_Pz_aa = I_ESP_O6y5z_S_aa+ABZ*I_ESP_N6y4z_S_aa;
    Double I_ESP_N5y5z_Pz_aa = I_ESP_O5y6z_S_aa+ABZ*I_ESP_N5y5z_S_aa;
    Double I_ESP_N4y6z_Pz_aa = I_ESP_O4y7z_S_aa+ABZ*I_ESP_N4y6z_S_aa;
    Double I_ESP_N3y7z_Pz_aa = I_ESP_O3y8z_S_aa+ABZ*I_ESP_N3y7z_S_aa;
    Double I_ESP_N2y8z_Pz_aa = I_ESP_O2y9z_S_aa+ABZ*I_ESP_N2y8z_S_aa;
    Double I_ESP_Ny9z_Pz_aa = I_ESP_Oy10z_S_aa+ABZ*I_ESP_Ny9z_S_aa;
    Double I_ESP_N10z_Pz_aa = I_ESP_O11z_S_aa+ABZ*I_ESP_N10z_S_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_M_D_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 168 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_N_P_aa
     * RHS shell quartet name: SQ_ESP_M_P_aa
     ************************************************************/
    Double I_ESP_M9x_D2x_aa = I_ESP_N10x_Px_aa+ABX*I_ESP_M9x_Px_aa;
    Double I_ESP_M8xy_D2x_aa = I_ESP_N9xy_Px_aa+ABX*I_ESP_M8xy_Px_aa;
    Double I_ESP_M8xz_D2x_aa = I_ESP_N9xz_Px_aa+ABX*I_ESP_M8xz_Px_aa;
    Double I_ESP_M7x2y_D2x_aa = I_ESP_N8x2y_Px_aa+ABX*I_ESP_M7x2y_Px_aa;
    Double I_ESP_M7xyz_D2x_aa = I_ESP_N8xyz_Px_aa+ABX*I_ESP_M7xyz_Px_aa;
    Double I_ESP_M7x2z_D2x_aa = I_ESP_N8x2z_Px_aa+ABX*I_ESP_M7x2z_Px_aa;
    Double I_ESP_M6x3y_D2x_aa = I_ESP_N7x3y_Px_aa+ABX*I_ESP_M6x3y_Px_aa;
    Double I_ESP_M6x2yz_D2x_aa = I_ESP_N7x2yz_Px_aa+ABX*I_ESP_M6x2yz_Px_aa;
    Double I_ESP_M6xy2z_D2x_aa = I_ESP_N7xy2z_Px_aa+ABX*I_ESP_M6xy2z_Px_aa;
    Double I_ESP_M6x3z_D2x_aa = I_ESP_N7x3z_Px_aa+ABX*I_ESP_M6x3z_Px_aa;
    Double I_ESP_M5x4y_D2x_aa = I_ESP_N6x4y_Px_aa+ABX*I_ESP_M5x4y_Px_aa;
    Double I_ESP_M5x3yz_D2x_aa = I_ESP_N6x3yz_Px_aa+ABX*I_ESP_M5x3yz_Px_aa;
    Double I_ESP_M5x2y2z_D2x_aa = I_ESP_N6x2y2z_Px_aa+ABX*I_ESP_M5x2y2z_Px_aa;
    Double I_ESP_M5xy3z_D2x_aa = I_ESP_N6xy3z_Px_aa+ABX*I_ESP_M5xy3z_Px_aa;
    Double I_ESP_M5x4z_D2x_aa = I_ESP_N6x4z_Px_aa+ABX*I_ESP_M5x4z_Px_aa;
    Double I_ESP_M4x5y_D2x_aa = I_ESP_N5x5y_Px_aa+ABX*I_ESP_M4x5y_Px_aa;
    Double I_ESP_M4x4yz_D2x_aa = I_ESP_N5x4yz_Px_aa+ABX*I_ESP_M4x4yz_Px_aa;
    Double I_ESP_M4x3y2z_D2x_aa = I_ESP_N5x3y2z_Px_aa+ABX*I_ESP_M4x3y2z_Px_aa;
    Double I_ESP_M4x2y3z_D2x_aa = I_ESP_N5x2y3z_Px_aa+ABX*I_ESP_M4x2y3z_Px_aa;
    Double I_ESP_M4xy4z_D2x_aa = I_ESP_N5xy4z_Px_aa+ABX*I_ESP_M4xy4z_Px_aa;
    Double I_ESP_M4x5z_D2x_aa = I_ESP_N5x5z_Px_aa+ABX*I_ESP_M4x5z_Px_aa;
    Double I_ESP_M3x6y_D2x_aa = I_ESP_N4x6y_Px_aa+ABX*I_ESP_M3x6y_Px_aa;
    Double I_ESP_M3x5yz_D2x_aa = I_ESP_N4x5yz_Px_aa+ABX*I_ESP_M3x5yz_Px_aa;
    Double I_ESP_M3x4y2z_D2x_aa = I_ESP_N4x4y2z_Px_aa+ABX*I_ESP_M3x4y2z_Px_aa;
    Double I_ESP_M3x3y3z_D2x_aa = I_ESP_N4x3y3z_Px_aa+ABX*I_ESP_M3x3y3z_Px_aa;
    Double I_ESP_M3x2y4z_D2x_aa = I_ESP_N4x2y4z_Px_aa+ABX*I_ESP_M3x2y4z_Px_aa;
    Double I_ESP_M3xy5z_D2x_aa = I_ESP_N4xy5z_Px_aa+ABX*I_ESP_M3xy5z_Px_aa;
    Double I_ESP_M3x6z_D2x_aa = I_ESP_N4x6z_Px_aa+ABX*I_ESP_M3x6z_Px_aa;
    Double I_ESP_M2x7y_D2x_aa = I_ESP_N3x7y_Px_aa+ABX*I_ESP_M2x7y_Px_aa;
    Double I_ESP_M2x6yz_D2x_aa = I_ESP_N3x6yz_Px_aa+ABX*I_ESP_M2x6yz_Px_aa;
    Double I_ESP_M2x5y2z_D2x_aa = I_ESP_N3x5y2z_Px_aa+ABX*I_ESP_M2x5y2z_Px_aa;
    Double I_ESP_M2x4y3z_D2x_aa = I_ESP_N3x4y3z_Px_aa+ABX*I_ESP_M2x4y3z_Px_aa;
    Double I_ESP_M2x3y4z_D2x_aa = I_ESP_N3x3y4z_Px_aa+ABX*I_ESP_M2x3y4z_Px_aa;
    Double I_ESP_M2x2y5z_D2x_aa = I_ESP_N3x2y5z_Px_aa+ABX*I_ESP_M2x2y5z_Px_aa;
    Double I_ESP_M2xy6z_D2x_aa = I_ESP_N3xy6z_Px_aa+ABX*I_ESP_M2xy6z_Px_aa;
    Double I_ESP_M2x7z_D2x_aa = I_ESP_N3x7z_Px_aa+ABX*I_ESP_M2x7z_Px_aa;
    Double I_ESP_Mx8y_D2x_aa = I_ESP_N2x8y_Px_aa+ABX*I_ESP_Mx8y_Px_aa;
    Double I_ESP_Mx7yz_D2x_aa = I_ESP_N2x7yz_Px_aa+ABX*I_ESP_Mx7yz_Px_aa;
    Double I_ESP_Mx6y2z_D2x_aa = I_ESP_N2x6y2z_Px_aa+ABX*I_ESP_Mx6y2z_Px_aa;
    Double I_ESP_Mx5y3z_D2x_aa = I_ESP_N2x5y3z_Px_aa+ABX*I_ESP_Mx5y3z_Px_aa;
    Double I_ESP_Mx4y4z_D2x_aa = I_ESP_N2x4y4z_Px_aa+ABX*I_ESP_Mx4y4z_Px_aa;
    Double I_ESP_Mx3y5z_D2x_aa = I_ESP_N2x3y5z_Px_aa+ABX*I_ESP_Mx3y5z_Px_aa;
    Double I_ESP_Mx2y6z_D2x_aa = I_ESP_N2x2y6z_Px_aa+ABX*I_ESP_Mx2y6z_Px_aa;
    Double I_ESP_Mxy7z_D2x_aa = I_ESP_N2xy7z_Px_aa+ABX*I_ESP_Mxy7z_Px_aa;
    Double I_ESP_Mx8z_D2x_aa = I_ESP_N2x8z_Px_aa+ABX*I_ESP_Mx8z_Px_aa;
    Double I_ESP_M9y_D2x_aa = I_ESP_Nx9y_Px_aa+ABX*I_ESP_M9y_Px_aa;
    Double I_ESP_M8yz_D2x_aa = I_ESP_Nx8yz_Px_aa+ABX*I_ESP_M8yz_Px_aa;
    Double I_ESP_M7y2z_D2x_aa = I_ESP_Nx7y2z_Px_aa+ABX*I_ESP_M7y2z_Px_aa;
    Double I_ESP_M6y3z_D2x_aa = I_ESP_Nx6y3z_Px_aa+ABX*I_ESP_M6y3z_Px_aa;
    Double I_ESP_M5y4z_D2x_aa = I_ESP_Nx5y4z_Px_aa+ABX*I_ESP_M5y4z_Px_aa;
    Double I_ESP_M4y5z_D2x_aa = I_ESP_Nx4y5z_Px_aa+ABX*I_ESP_M4y5z_Px_aa;
    Double I_ESP_M3y6z_D2x_aa = I_ESP_Nx3y6z_Px_aa+ABX*I_ESP_M3y6z_Px_aa;
    Double I_ESP_M2y7z_D2x_aa = I_ESP_Nx2y7z_Px_aa+ABX*I_ESP_M2y7z_Px_aa;
    Double I_ESP_My8z_D2x_aa = I_ESP_Nxy8z_Px_aa+ABX*I_ESP_My8z_Px_aa;
    Double I_ESP_M9z_D2x_aa = I_ESP_Nx9z_Px_aa+ABX*I_ESP_M9z_Px_aa;
    Double I_ESP_M8xy_D2y_aa = I_ESP_N8x2y_Py_aa+ABY*I_ESP_M8xy_Py_aa;
    Double I_ESP_M8xz_D2y_aa = I_ESP_N8xyz_Py_aa+ABY*I_ESP_M8xz_Py_aa;
    Double I_ESP_M7x2y_D2y_aa = I_ESP_N7x3y_Py_aa+ABY*I_ESP_M7x2y_Py_aa;
    Double I_ESP_M7xyz_D2y_aa = I_ESP_N7x2yz_Py_aa+ABY*I_ESP_M7xyz_Py_aa;
    Double I_ESP_M7x2z_D2y_aa = I_ESP_N7xy2z_Py_aa+ABY*I_ESP_M7x2z_Py_aa;
    Double I_ESP_M6x3y_D2y_aa = I_ESP_N6x4y_Py_aa+ABY*I_ESP_M6x3y_Py_aa;
    Double I_ESP_M6x2yz_D2y_aa = I_ESP_N6x3yz_Py_aa+ABY*I_ESP_M6x2yz_Py_aa;
    Double I_ESP_M6xy2z_D2y_aa = I_ESP_N6x2y2z_Py_aa+ABY*I_ESP_M6xy2z_Py_aa;
    Double I_ESP_M6x3z_D2y_aa = I_ESP_N6xy3z_Py_aa+ABY*I_ESP_M6x3z_Py_aa;
    Double I_ESP_M5x4y_D2y_aa = I_ESP_N5x5y_Py_aa+ABY*I_ESP_M5x4y_Py_aa;
    Double I_ESP_M5x3yz_D2y_aa = I_ESP_N5x4yz_Py_aa+ABY*I_ESP_M5x3yz_Py_aa;
    Double I_ESP_M5x2y2z_D2y_aa = I_ESP_N5x3y2z_Py_aa+ABY*I_ESP_M5x2y2z_Py_aa;
    Double I_ESP_M5xy3z_D2y_aa = I_ESP_N5x2y3z_Py_aa+ABY*I_ESP_M5xy3z_Py_aa;
    Double I_ESP_M5x4z_D2y_aa = I_ESP_N5xy4z_Py_aa+ABY*I_ESP_M5x4z_Py_aa;
    Double I_ESP_M4x5y_D2y_aa = I_ESP_N4x6y_Py_aa+ABY*I_ESP_M4x5y_Py_aa;
    Double I_ESP_M4x4yz_D2y_aa = I_ESP_N4x5yz_Py_aa+ABY*I_ESP_M4x4yz_Py_aa;
    Double I_ESP_M4x3y2z_D2y_aa = I_ESP_N4x4y2z_Py_aa+ABY*I_ESP_M4x3y2z_Py_aa;
    Double I_ESP_M4x2y3z_D2y_aa = I_ESP_N4x3y3z_Py_aa+ABY*I_ESP_M4x2y3z_Py_aa;
    Double I_ESP_M4xy4z_D2y_aa = I_ESP_N4x2y4z_Py_aa+ABY*I_ESP_M4xy4z_Py_aa;
    Double I_ESP_M4x5z_D2y_aa = I_ESP_N4xy5z_Py_aa+ABY*I_ESP_M4x5z_Py_aa;
    Double I_ESP_M3x6y_D2y_aa = I_ESP_N3x7y_Py_aa+ABY*I_ESP_M3x6y_Py_aa;
    Double I_ESP_M3x5yz_D2y_aa = I_ESP_N3x6yz_Py_aa+ABY*I_ESP_M3x5yz_Py_aa;
    Double I_ESP_M3x4y2z_D2y_aa = I_ESP_N3x5y2z_Py_aa+ABY*I_ESP_M3x4y2z_Py_aa;
    Double I_ESP_M3x3y3z_D2y_aa = I_ESP_N3x4y3z_Py_aa+ABY*I_ESP_M3x3y3z_Py_aa;
    Double I_ESP_M3x2y4z_D2y_aa = I_ESP_N3x3y4z_Py_aa+ABY*I_ESP_M3x2y4z_Py_aa;
    Double I_ESP_M3xy5z_D2y_aa = I_ESP_N3x2y5z_Py_aa+ABY*I_ESP_M3xy5z_Py_aa;
    Double I_ESP_M3x6z_D2y_aa = I_ESP_N3xy6z_Py_aa+ABY*I_ESP_M3x6z_Py_aa;
    Double I_ESP_M2x7y_D2y_aa = I_ESP_N2x8y_Py_aa+ABY*I_ESP_M2x7y_Py_aa;
    Double I_ESP_M2x6yz_D2y_aa = I_ESP_N2x7yz_Py_aa+ABY*I_ESP_M2x6yz_Py_aa;
    Double I_ESP_M2x5y2z_D2y_aa = I_ESP_N2x6y2z_Py_aa+ABY*I_ESP_M2x5y2z_Py_aa;
    Double I_ESP_M2x4y3z_D2y_aa = I_ESP_N2x5y3z_Py_aa+ABY*I_ESP_M2x4y3z_Py_aa;
    Double I_ESP_M2x3y4z_D2y_aa = I_ESP_N2x4y4z_Py_aa+ABY*I_ESP_M2x3y4z_Py_aa;
    Double I_ESP_M2x2y5z_D2y_aa = I_ESP_N2x3y5z_Py_aa+ABY*I_ESP_M2x2y5z_Py_aa;
    Double I_ESP_M2xy6z_D2y_aa = I_ESP_N2x2y6z_Py_aa+ABY*I_ESP_M2xy6z_Py_aa;
    Double I_ESP_M2x7z_D2y_aa = I_ESP_N2xy7z_Py_aa+ABY*I_ESP_M2x7z_Py_aa;
    Double I_ESP_Mx8y_D2y_aa = I_ESP_Nx9y_Py_aa+ABY*I_ESP_Mx8y_Py_aa;
    Double I_ESP_Mx7yz_D2y_aa = I_ESP_Nx8yz_Py_aa+ABY*I_ESP_Mx7yz_Py_aa;
    Double I_ESP_Mx6y2z_D2y_aa = I_ESP_Nx7y2z_Py_aa+ABY*I_ESP_Mx6y2z_Py_aa;
    Double I_ESP_Mx5y3z_D2y_aa = I_ESP_Nx6y3z_Py_aa+ABY*I_ESP_Mx5y3z_Py_aa;
    Double I_ESP_Mx4y4z_D2y_aa = I_ESP_Nx5y4z_Py_aa+ABY*I_ESP_Mx4y4z_Py_aa;
    Double I_ESP_Mx3y5z_D2y_aa = I_ESP_Nx4y5z_Py_aa+ABY*I_ESP_Mx3y5z_Py_aa;
    Double I_ESP_Mx2y6z_D2y_aa = I_ESP_Nx3y6z_Py_aa+ABY*I_ESP_Mx2y6z_Py_aa;
    Double I_ESP_Mxy7z_D2y_aa = I_ESP_Nx2y7z_Py_aa+ABY*I_ESP_Mxy7z_Py_aa;
    Double I_ESP_Mx8z_D2y_aa = I_ESP_Nxy8z_Py_aa+ABY*I_ESP_Mx8z_Py_aa;
    Double I_ESP_M9y_D2y_aa = I_ESP_N10y_Py_aa+ABY*I_ESP_M9y_Py_aa;
    Double I_ESP_M8yz_D2y_aa = I_ESP_N9yz_Py_aa+ABY*I_ESP_M8yz_Py_aa;
    Double I_ESP_M7y2z_D2y_aa = I_ESP_N8y2z_Py_aa+ABY*I_ESP_M7y2z_Py_aa;
    Double I_ESP_M6y3z_D2y_aa = I_ESP_N7y3z_Py_aa+ABY*I_ESP_M6y3z_Py_aa;
    Double I_ESP_M5y4z_D2y_aa = I_ESP_N6y4z_Py_aa+ABY*I_ESP_M5y4z_Py_aa;
    Double I_ESP_M4y5z_D2y_aa = I_ESP_N5y5z_Py_aa+ABY*I_ESP_M4y5z_Py_aa;
    Double I_ESP_M3y6z_D2y_aa = I_ESP_N4y6z_Py_aa+ABY*I_ESP_M3y6z_Py_aa;
    Double I_ESP_M2y7z_D2y_aa = I_ESP_N3y7z_Py_aa+ABY*I_ESP_M2y7z_Py_aa;
    Double I_ESP_My8z_D2y_aa = I_ESP_N2y8z_Py_aa+ABY*I_ESP_My8z_Py_aa;
    Double I_ESP_M9z_D2y_aa = I_ESP_Ny9z_Py_aa+ABY*I_ESP_M9z_Py_aa;
    Double I_ESP_M8xy_D2z_aa = I_ESP_N8xyz_Pz_aa+ABZ*I_ESP_M8xy_Pz_aa;
    Double I_ESP_M8xz_D2z_aa = I_ESP_N8x2z_Pz_aa+ABZ*I_ESP_M8xz_Pz_aa;
    Double I_ESP_M7x2y_D2z_aa = I_ESP_N7x2yz_Pz_aa+ABZ*I_ESP_M7x2y_Pz_aa;
    Double I_ESP_M7xyz_D2z_aa = I_ESP_N7xy2z_Pz_aa+ABZ*I_ESP_M7xyz_Pz_aa;
    Double I_ESP_M7x2z_D2z_aa = I_ESP_N7x3z_Pz_aa+ABZ*I_ESP_M7x2z_Pz_aa;
    Double I_ESP_M6x3y_D2z_aa = I_ESP_N6x3yz_Pz_aa+ABZ*I_ESP_M6x3y_Pz_aa;
    Double I_ESP_M6x2yz_D2z_aa = I_ESP_N6x2y2z_Pz_aa+ABZ*I_ESP_M6x2yz_Pz_aa;
    Double I_ESP_M6xy2z_D2z_aa = I_ESP_N6xy3z_Pz_aa+ABZ*I_ESP_M6xy2z_Pz_aa;
    Double I_ESP_M6x3z_D2z_aa = I_ESP_N6x4z_Pz_aa+ABZ*I_ESP_M6x3z_Pz_aa;
    Double I_ESP_M5x4y_D2z_aa = I_ESP_N5x4yz_Pz_aa+ABZ*I_ESP_M5x4y_Pz_aa;
    Double I_ESP_M5x3yz_D2z_aa = I_ESP_N5x3y2z_Pz_aa+ABZ*I_ESP_M5x3yz_Pz_aa;
    Double I_ESP_M5x2y2z_D2z_aa = I_ESP_N5x2y3z_Pz_aa+ABZ*I_ESP_M5x2y2z_Pz_aa;
    Double I_ESP_M5xy3z_D2z_aa = I_ESP_N5xy4z_Pz_aa+ABZ*I_ESP_M5xy3z_Pz_aa;
    Double I_ESP_M5x4z_D2z_aa = I_ESP_N5x5z_Pz_aa+ABZ*I_ESP_M5x4z_Pz_aa;
    Double I_ESP_M4x5y_D2z_aa = I_ESP_N4x5yz_Pz_aa+ABZ*I_ESP_M4x5y_Pz_aa;
    Double I_ESP_M4x4yz_D2z_aa = I_ESP_N4x4y2z_Pz_aa+ABZ*I_ESP_M4x4yz_Pz_aa;
    Double I_ESP_M4x3y2z_D2z_aa = I_ESP_N4x3y3z_Pz_aa+ABZ*I_ESP_M4x3y2z_Pz_aa;
    Double I_ESP_M4x2y3z_D2z_aa = I_ESP_N4x2y4z_Pz_aa+ABZ*I_ESP_M4x2y3z_Pz_aa;
    Double I_ESP_M4xy4z_D2z_aa = I_ESP_N4xy5z_Pz_aa+ABZ*I_ESP_M4xy4z_Pz_aa;
    Double I_ESP_M4x5z_D2z_aa = I_ESP_N4x6z_Pz_aa+ABZ*I_ESP_M4x5z_Pz_aa;
    Double I_ESP_M3x6y_D2z_aa = I_ESP_N3x6yz_Pz_aa+ABZ*I_ESP_M3x6y_Pz_aa;
    Double I_ESP_M3x5yz_D2z_aa = I_ESP_N3x5y2z_Pz_aa+ABZ*I_ESP_M3x5yz_Pz_aa;
    Double I_ESP_M3x4y2z_D2z_aa = I_ESP_N3x4y3z_Pz_aa+ABZ*I_ESP_M3x4y2z_Pz_aa;
    Double I_ESP_M3x3y3z_D2z_aa = I_ESP_N3x3y4z_Pz_aa+ABZ*I_ESP_M3x3y3z_Pz_aa;
    Double I_ESP_M3x2y4z_D2z_aa = I_ESP_N3x2y5z_Pz_aa+ABZ*I_ESP_M3x2y4z_Pz_aa;
    Double I_ESP_M3xy5z_D2z_aa = I_ESP_N3xy6z_Pz_aa+ABZ*I_ESP_M3xy5z_Pz_aa;
    Double I_ESP_M3x6z_D2z_aa = I_ESP_N3x7z_Pz_aa+ABZ*I_ESP_M3x6z_Pz_aa;
    Double I_ESP_M2x7y_D2z_aa = I_ESP_N2x7yz_Pz_aa+ABZ*I_ESP_M2x7y_Pz_aa;
    Double I_ESP_M2x6yz_D2z_aa = I_ESP_N2x6y2z_Pz_aa+ABZ*I_ESP_M2x6yz_Pz_aa;
    Double I_ESP_M2x5y2z_D2z_aa = I_ESP_N2x5y3z_Pz_aa+ABZ*I_ESP_M2x5y2z_Pz_aa;
    Double I_ESP_M2x4y3z_D2z_aa = I_ESP_N2x4y4z_Pz_aa+ABZ*I_ESP_M2x4y3z_Pz_aa;
    Double I_ESP_M2x3y4z_D2z_aa = I_ESP_N2x3y5z_Pz_aa+ABZ*I_ESP_M2x3y4z_Pz_aa;
    Double I_ESP_M2x2y5z_D2z_aa = I_ESP_N2x2y6z_Pz_aa+ABZ*I_ESP_M2x2y5z_Pz_aa;
    Double I_ESP_M2xy6z_D2z_aa = I_ESP_N2xy7z_Pz_aa+ABZ*I_ESP_M2xy6z_Pz_aa;
    Double I_ESP_M2x7z_D2z_aa = I_ESP_N2x8z_Pz_aa+ABZ*I_ESP_M2x7z_Pz_aa;
    Double I_ESP_Mx8y_D2z_aa = I_ESP_Nx8yz_Pz_aa+ABZ*I_ESP_Mx8y_Pz_aa;
    Double I_ESP_Mx7yz_D2z_aa = I_ESP_Nx7y2z_Pz_aa+ABZ*I_ESP_Mx7yz_Pz_aa;
    Double I_ESP_Mx6y2z_D2z_aa = I_ESP_Nx6y3z_Pz_aa+ABZ*I_ESP_Mx6y2z_Pz_aa;
    Double I_ESP_Mx5y3z_D2z_aa = I_ESP_Nx5y4z_Pz_aa+ABZ*I_ESP_Mx5y3z_Pz_aa;
    Double I_ESP_Mx4y4z_D2z_aa = I_ESP_Nx4y5z_Pz_aa+ABZ*I_ESP_Mx4y4z_Pz_aa;
    Double I_ESP_Mx3y5z_D2z_aa = I_ESP_Nx3y6z_Pz_aa+ABZ*I_ESP_Mx3y5z_Pz_aa;
    Double I_ESP_Mx2y6z_D2z_aa = I_ESP_Nx2y7z_Pz_aa+ABZ*I_ESP_Mx2y6z_Pz_aa;
    Double I_ESP_Mxy7z_D2z_aa = I_ESP_Nxy8z_Pz_aa+ABZ*I_ESP_Mxy7z_Pz_aa;
    Double I_ESP_Mx8z_D2z_aa = I_ESP_Nx9z_Pz_aa+ABZ*I_ESP_Mx8z_Pz_aa;
    Double I_ESP_M8yz_D2z_aa = I_ESP_N8y2z_Pz_aa+ABZ*I_ESP_M8yz_Pz_aa;
    Double I_ESP_M7y2z_D2z_aa = I_ESP_N7y3z_Pz_aa+ABZ*I_ESP_M7y2z_Pz_aa;
    Double I_ESP_M6y3z_D2z_aa = I_ESP_N6y4z_Pz_aa+ABZ*I_ESP_M6y3z_Pz_aa;
    Double I_ESP_M5y4z_D2z_aa = I_ESP_N5y5z_Pz_aa+ABZ*I_ESP_M5y4z_Pz_aa;
    Double I_ESP_M4y5z_D2z_aa = I_ESP_N4y6z_Pz_aa+ABZ*I_ESP_M4y5z_Pz_aa;
    Double I_ESP_M3y6z_D2z_aa = I_ESP_N3y7z_Pz_aa+ABZ*I_ESP_M3y6z_Pz_aa;
    Double I_ESP_M2y7z_D2z_aa = I_ESP_N2y8z_Pz_aa+ABZ*I_ESP_M2y7z_Pz_aa;
    Double I_ESP_My8z_D2z_aa = I_ESP_Ny9z_Pz_aa+ABZ*I_ESP_My8z_Pz_aa;
    Double I_ESP_M9z_D2z_aa = I_ESP_N10z_Pz_aa+ABZ*I_ESP_M9z_Pz_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_L_F_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 127 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_M_D_aa
     * RHS shell quartet name: SQ_ESP_L_D_aa
     ************************************************************/
    Double I_ESP_L8x_F3x_aa = I_ESP_M9x_D2x_aa+ABX*I_ESP_L8x_D2x_aa;
    Double I_ESP_L7xy_F3x_aa = I_ESP_M8xy_D2x_aa+ABX*I_ESP_L7xy_D2x_aa;
    Double I_ESP_L7xz_F3x_aa = I_ESP_M8xz_D2x_aa+ABX*I_ESP_L7xz_D2x_aa;
    Double I_ESP_L6x2y_F3x_aa = I_ESP_M7x2y_D2x_aa+ABX*I_ESP_L6x2y_D2x_aa;
    Double I_ESP_L6xyz_F3x_aa = I_ESP_M7xyz_D2x_aa+ABX*I_ESP_L6xyz_D2x_aa;
    Double I_ESP_L6x2z_F3x_aa = I_ESP_M7x2z_D2x_aa+ABX*I_ESP_L6x2z_D2x_aa;
    Double I_ESP_L5x3y_F3x_aa = I_ESP_M6x3y_D2x_aa+ABX*I_ESP_L5x3y_D2x_aa;
    Double I_ESP_L5x2yz_F3x_aa = I_ESP_M6x2yz_D2x_aa+ABX*I_ESP_L5x2yz_D2x_aa;
    Double I_ESP_L5xy2z_F3x_aa = I_ESP_M6xy2z_D2x_aa+ABX*I_ESP_L5xy2z_D2x_aa;
    Double I_ESP_L5x3z_F3x_aa = I_ESP_M6x3z_D2x_aa+ABX*I_ESP_L5x3z_D2x_aa;
    Double I_ESP_L4x4y_F3x_aa = I_ESP_M5x4y_D2x_aa+ABX*I_ESP_L4x4y_D2x_aa;
    Double I_ESP_L4x3yz_F3x_aa = I_ESP_M5x3yz_D2x_aa+ABX*I_ESP_L4x3yz_D2x_aa;
    Double I_ESP_L4x2y2z_F3x_aa = I_ESP_M5x2y2z_D2x_aa+ABX*I_ESP_L4x2y2z_D2x_aa;
    Double I_ESP_L4xy3z_F3x_aa = I_ESP_M5xy3z_D2x_aa+ABX*I_ESP_L4xy3z_D2x_aa;
    Double I_ESP_L4x4z_F3x_aa = I_ESP_M5x4z_D2x_aa+ABX*I_ESP_L4x4z_D2x_aa;
    Double I_ESP_L3x5y_F3x_aa = I_ESP_M4x5y_D2x_aa+ABX*I_ESP_L3x5y_D2x_aa;
    Double I_ESP_L3x4yz_F3x_aa = I_ESP_M4x4yz_D2x_aa+ABX*I_ESP_L3x4yz_D2x_aa;
    Double I_ESP_L3x3y2z_F3x_aa = I_ESP_M4x3y2z_D2x_aa+ABX*I_ESP_L3x3y2z_D2x_aa;
    Double I_ESP_L3x2y3z_F3x_aa = I_ESP_M4x2y3z_D2x_aa+ABX*I_ESP_L3x2y3z_D2x_aa;
    Double I_ESP_L3xy4z_F3x_aa = I_ESP_M4xy4z_D2x_aa+ABX*I_ESP_L3xy4z_D2x_aa;
    Double I_ESP_L3x5z_F3x_aa = I_ESP_M4x5z_D2x_aa+ABX*I_ESP_L3x5z_D2x_aa;
    Double I_ESP_L2x6y_F3x_aa = I_ESP_M3x6y_D2x_aa+ABX*I_ESP_L2x6y_D2x_aa;
    Double I_ESP_L2x5yz_F3x_aa = I_ESP_M3x5yz_D2x_aa+ABX*I_ESP_L2x5yz_D2x_aa;
    Double I_ESP_L2x4y2z_F3x_aa = I_ESP_M3x4y2z_D2x_aa+ABX*I_ESP_L2x4y2z_D2x_aa;
    Double I_ESP_L2x3y3z_F3x_aa = I_ESP_M3x3y3z_D2x_aa+ABX*I_ESP_L2x3y3z_D2x_aa;
    Double I_ESP_L2x2y4z_F3x_aa = I_ESP_M3x2y4z_D2x_aa+ABX*I_ESP_L2x2y4z_D2x_aa;
    Double I_ESP_L2xy5z_F3x_aa = I_ESP_M3xy5z_D2x_aa+ABX*I_ESP_L2xy5z_D2x_aa;
    Double I_ESP_L2x6z_F3x_aa = I_ESP_M3x6z_D2x_aa+ABX*I_ESP_L2x6z_D2x_aa;
    Double I_ESP_Lx7y_F3x_aa = I_ESP_M2x7y_D2x_aa+ABX*I_ESP_Lx7y_D2x_aa;
    Double I_ESP_Lx6yz_F3x_aa = I_ESP_M2x6yz_D2x_aa+ABX*I_ESP_Lx6yz_D2x_aa;
    Double I_ESP_Lx5y2z_F3x_aa = I_ESP_M2x5y2z_D2x_aa+ABX*I_ESP_Lx5y2z_D2x_aa;
    Double I_ESP_Lx4y3z_F3x_aa = I_ESP_M2x4y3z_D2x_aa+ABX*I_ESP_Lx4y3z_D2x_aa;
    Double I_ESP_Lx3y4z_F3x_aa = I_ESP_M2x3y4z_D2x_aa+ABX*I_ESP_Lx3y4z_D2x_aa;
    Double I_ESP_Lx2y5z_F3x_aa = I_ESP_M2x2y5z_D2x_aa+ABX*I_ESP_Lx2y5z_D2x_aa;
    Double I_ESP_Lxy6z_F3x_aa = I_ESP_M2xy6z_D2x_aa+ABX*I_ESP_Lxy6z_D2x_aa;
    Double I_ESP_Lx7z_F3x_aa = I_ESP_M2x7z_D2x_aa+ABX*I_ESP_Lx7z_D2x_aa;
    Double I_ESP_L8y_F3x_aa = I_ESP_Mx8y_D2x_aa+ABX*I_ESP_L8y_D2x_aa;
    Double I_ESP_L7yz_F3x_aa = I_ESP_Mx7yz_D2x_aa+ABX*I_ESP_L7yz_D2x_aa;
    Double I_ESP_L6y2z_F3x_aa = I_ESP_Mx6y2z_D2x_aa+ABX*I_ESP_L6y2z_D2x_aa;
    Double I_ESP_L5y3z_F3x_aa = I_ESP_Mx5y3z_D2x_aa+ABX*I_ESP_L5y3z_D2x_aa;
    Double I_ESP_L4y4z_F3x_aa = I_ESP_Mx4y4z_D2x_aa+ABX*I_ESP_L4y4z_D2x_aa;
    Double I_ESP_L3y5z_F3x_aa = I_ESP_Mx3y5z_D2x_aa+ABX*I_ESP_L3y5z_D2x_aa;
    Double I_ESP_L2y6z_F3x_aa = I_ESP_Mx2y6z_D2x_aa+ABX*I_ESP_L2y6z_D2x_aa;
    Double I_ESP_Ly7z_F3x_aa = I_ESP_Mxy7z_D2x_aa+ABX*I_ESP_Ly7z_D2x_aa;
    Double I_ESP_L8z_F3x_aa = I_ESP_Mx8z_D2x_aa+ABX*I_ESP_L8z_D2x_aa;
    Double I_ESP_L7xy_F2xy_aa = I_ESP_M7x2y_D2x_aa+ABY*I_ESP_L7xy_D2x_aa;
    Double I_ESP_L7xz_F2xy_aa = I_ESP_M7xyz_D2x_aa+ABY*I_ESP_L7xz_D2x_aa;
    Double I_ESP_L6x2y_F2xy_aa = I_ESP_M6x3y_D2x_aa+ABY*I_ESP_L6x2y_D2x_aa;
    Double I_ESP_L6xyz_F2xy_aa = I_ESP_M6x2yz_D2x_aa+ABY*I_ESP_L6xyz_D2x_aa;
    Double I_ESP_L6x2z_F2xy_aa = I_ESP_M6xy2z_D2x_aa+ABY*I_ESP_L6x2z_D2x_aa;
    Double I_ESP_L5x3y_F2xy_aa = I_ESP_M5x4y_D2x_aa+ABY*I_ESP_L5x3y_D2x_aa;
    Double I_ESP_L5x2yz_F2xy_aa = I_ESP_M5x3yz_D2x_aa+ABY*I_ESP_L5x2yz_D2x_aa;
    Double I_ESP_L5xy2z_F2xy_aa = I_ESP_M5x2y2z_D2x_aa+ABY*I_ESP_L5xy2z_D2x_aa;
    Double I_ESP_L5x3z_F2xy_aa = I_ESP_M5xy3z_D2x_aa+ABY*I_ESP_L5x3z_D2x_aa;
    Double I_ESP_L4x4y_F2xy_aa = I_ESP_M4x5y_D2x_aa+ABY*I_ESP_L4x4y_D2x_aa;
    Double I_ESP_L4x3yz_F2xy_aa = I_ESP_M4x4yz_D2x_aa+ABY*I_ESP_L4x3yz_D2x_aa;
    Double I_ESP_L4x2y2z_F2xy_aa = I_ESP_M4x3y2z_D2x_aa+ABY*I_ESP_L4x2y2z_D2x_aa;
    Double I_ESP_L4xy3z_F2xy_aa = I_ESP_M4x2y3z_D2x_aa+ABY*I_ESP_L4xy3z_D2x_aa;
    Double I_ESP_L4x4z_F2xy_aa = I_ESP_M4xy4z_D2x_aa+ABY*I_ESP_L4x4z_D2x_aa;
    Double I_ESP_L3x5y_F2xy_aa = I_ESP_M3x6y_D2x_aa+ABY*I_ESP_L3x5y_D2x_aa;
    Double I_ESP_L3x4yz_F2xy_aa = I_ESP_M3x5yz_D2x_aa+ABY*I_ESP_L3x4yz_D2x_aa;
    Double I_ESP_L3x3y2z_F2xy_aa = I_ESP_M3x4y2z_D2x_aa+ABY*I_ESP_L3x3y2z_D2x_aa;
    Double I_ESP_L3x2y3z_F2xy_aa = I_ESP_M3x3y3z_D2x_aa+ABY*I_ESP_L3x2y3z_D2x_aa;
    Double I_ESP_L3xy4z_F2xy_aa = I_ESP_M3x2y4z_D2x_aa+ABY*I_ESP_L3xy4z_D2x_aa;
    Double I_ESP_L3x5z_F2xy_aa = I_ESP_M3xy5z_D2x_aa+ABY*I_ESP_L3x5z_D2x_aa;
    Double I_ESP_L2x6y_F2xy_aa = I_ESP_M2x7y_D2x_aa+ABY*I_ESP_L2x6y_D2x_aa;
    Double I_ESP_L2x5yz_F2xy_aa = I_ESP_M2x6yz_D2x_aa+ABY*I_ESP_L2x5yz_D2x_aa;
    Double I_ESP_L2x4y2z_F2xy_aa = I_ESP_M2x5y2z_D2x_aa+ABY*I_ESP_L2x4y2z_D2x_aa;
    Double I_ESP_L2x3y3z_F2xy_aa = I_ESP_M2x4y3z_D2x_aa+ABY*I_ESP_L2x3y3z_D2x_aa;
    Double I_ESP_L2x2y4z_F2xy_aa = I_ESP_M2x3y4z_D2x_aa+ABY*I_ESP_L2x2y4z_D2x_aa;
    Double I_ESP_L2xy5z_F2xy_aa = I_ESP_M2x2y5z_D2x_aa+ABY*I_ESP_L2xy5z_D2x_aa;
    Double I_ESP_L2x6z_F2xy_aa = I_ESP_M2xy6z_D2x_aa+ABY*I_ESP_L2x6z_D2x_aa;
    Double I_ESP_Lx7y_F2xy_aa = I_ESP_Mx8y_D2x_aa+ABY*I_ESP_Lx7y_D2x_aa;
    Double I_ESP_Lx6yz_F2xy_aa = I_ESP_Mx7yz_D2x_aa+ABY*I_ESP_Lx6yz_D2x_aa;
    Double I_ESP_Lx5y2z_F2xy_aa = I_ESP_Mx6y2z_D2x_aa+ABY*I_ESP_Lx5y2z_D2x_aa;
    Double I_ESP_Lx4y3z_F2xy_aa = I_ESP_Mx5y3z_D2x_aa+ABY*I_ESP_Lx4y3z_D2x_aa;
    Double I_ESP_Lx3y4z_F2xy_aa = I_ESP_Mx4y4z_D2x_aa+ABY*I_ESP_Lx3y4z_D2x_aa;
    Double I_ESP_Lx2y5z_F2xy_aa = I_ESP_Mx3y5z_D2x_aa+ABY*I_ESP_Lx2y5z_D2x_aa;
    Double I_ESP_Lxy6z_F2xy_aa = I_ESP_Mx2y6z_D2x_aa+ABY*I_ESP_Lxy6z_D2x_aa;
    Double I_ESP_Lx7z_F2xy_aa = I_ESP_Mxy7z_D2x_aa+ABY*I_ESP_Lx7z_D2x_aa;
    Double I_ESP_L8y_F2xy_aa = I_ESP_M9y_D2x_aa+ABY*I_ESP_L8y_D2x_aa;
    Double I_ESP_L7yz_F2xy_aa = I_ESP_M8yz_D2x_aa+ABY*I_ESP_L7yz_D2x_aa;
    Double I_ESP_L6y2z_F2xy_aa = I_ESP_M7y2z_D2x_aa+ABY*I_ESP_L6y2z_D2x_aa;
    Double I_ESP_L5y3z_F2xy_aa = I_ESP_M6y3z_D2x_aa+ABY*I_ESP_L5y3z_D2x_aa;
    Double I_ESP_L4y4z_F2xy_aa = I_ESP_M5y4z_D2x_aa+ABY*I_ESP_L4y4z_D2x_aa;
    Double I_ESP_L3y5z_F2xy_aa = I_ESP_M4y5z_D2x_aa+ABY*I_ESP_L3y5z_D2x_aa;
    Double I_ESP_L2y6z_F2xy_aa = I_ESP_M3y6z_D2x_aa+ABY*I_ESP_L2y6z_D2x_aa;
    Double I_ESP_Ly7z_F2xy_aa = I_ESP_M2y7z_D2x_aa+ABY*I_ESP_Ly7z_D2x_aa;
    Double I_ESP_L8z_F2xy_aa = I_ESP_My8z_D2x_aa+ABY*I_ESP_L8z_D2x_aa;
    Double I_ESP_L7xz_F2xz_aa = I_ESP_M7x2z_D2x_aa+ABZ*I_ESP_L7xz_D2x_aa;
    Double I_ESP_L6xyz_F2xz_aa = I_ESP_M6xy2z_D2x_aa+ABZ*I_ESP_L6xyz_D2x_aa;
    Double I_ESP_L6x2z_F2xz_aa = I_ESP_M6x3z_D2x_aa+ABZ*I_ESP_L6x2z_D2x_aa;
    Double I_ESP_L5x2yz_F2xz_aa = I_ESP_M5x2y2z_D2x_aa+ABZ*I_ESP_L5x2yz_D2x_aa;
    Double I_ESP_L5xy2z_F2xz_aa = I_ESP_M5xy3z_D2x_aa+ABZ*I_ESP_L5xy2z_D2x_aa;
    Double I_ESP_L5x3z_F2xz_aa = I_ESP_M5x4z_D2x_aa+ABZ*I_ESP_L5x3z_D2x_aa;
    Double I_ESP_L4x3yz_F2xz_aa = I_ESP_M4x3y2z_D2x_aa+ABZ*I_ESP_L4x3yz_D2x_aa;
    Double I_ESP_L4x2y2z_F2xz_aa = I_ESP_M4x2y3z_D2x_aa+ABZ*I_ESP_L4x2y2z_D2x_aa;
    Double I_ESP_L4xy3z_F2xz_aa = I_ESP_M4xy4z_D2x_aa+ABZ*I_ESP_L4xy3z_D2x_aa;
    Double I_ESP_L4x4z_F2xz_aa = I_ESP_M4x5z_D2x_aa+ABZ*I_ESP_L4x4z_D2x_aa;
    Double I_ESP_L3x4yz_F2xz_aa = I_ESP_M3x4y2z_D2x_aa+ABZ*I_ESP_L3x4yz_D2x_aa;
    Double I_ESP_L3x3y2z_F2xz_aa = I_ESP_M3x3y3z_D2x_aa+ABZ*I_ESP_L3x3y2z_D2x_aa;
    Double I_ESP_L3x2y3z_F2xz_aa = I_ESP_M3x2y4z_D2x_aa+ABZ*I_ESP_L3x2y3z_D2x_aa;
    Double I_ESP_L3xy4z_F2xz_aa = I_ESP_M3xy5z_D2x_aa+ABZ*I_ESP_L3xy4z_D2x_aa;
    Double I_ESP_L3x5z_F2xz_aa = I_ESP_M3x6z_D2x_aa+ABZ*I_ESP_L3x5z_D2x_aa;
    Double I_ESP_L2x5yz_F2xz_aa = I_ESP_M2x5y2z_D2x_aa+ABZ*I_ESP_L2x5yz_D2x_aa;
    Double I_ESP_L2x4y2z_F2xz_aa = I_ESP_M2x4y3z_D2x_aa+ABZ*I_ESP_L2x4y2z_D2x_aa;
    Double I_ESP_L2x3y3z_F2xz_aa = I_ESP_M2x3y4z_D2x_aa+ABZ*I_ESP_L2x3y3z_D2x_aa;
    Double I_ESP_L2x2y4z_F2xz_aa = I_ESP_M2x2y5z_D2x_aa+ABZ*I_ESP_L2x2y4z_D2x_aa;
    Double I_ESP_L2xy5z_F2xz_aa = I_ESP_M2xy6z_D2x_aa+ABZ*I_ESP_L2xy5z_D2x_aa;
    Double I_ESP_L2x6z_F2xz_aa = I_ESP_M2x7z_D2x_aa+ABZ*I_ESP_L2x6z_D2x_aa;
    Double I_ESP_Lx6yz_F2xz_aa = I_ESP_Mx6y2z_D2x_aa+ABZ*I_ESP_Lx6yz_D2x_aa;
    Double I_ESP_Lx5y2z_F2xz_aa = I_ESP_Mx5y3z_D2x_aa+ABZ*I_ESP_Lx5y2z_D2x_aa;
    Double I_ESP_Lx4y3z_F2xz_aa = I_ESP_Mx4y4z_D2x_aa+ABZ*I_ESP_Lx4y3z_D2x_aa;
    Double I_ESP_Lx3y4z_F2xz_aa = I_ESP_Mx3y5z_D2x_aa+ABZ*I_ESP_Lx3y4z_D2x_aa;
    Double I_ESP_Lx2y5z_F2xz_aa = I_ESP_Mx2y6z_D2x_aa+ABZ*I_ESP_Lx2y5z_D2x_aa;
    Double I_ESP_Lxy6z_F2xz_aa = I_ESP_Mxy7z_D2x_aa+ABZ*I_ESP_Lxy6z_D2x_aa;
    Double I_ESP_Lx7z_F2xz_aa = I_ESP_Mx8z_D2x_aa+ABZ*I_ESP_Lx7z_D2x_aa;
    Double I_ESP_L7yz_F2xz_aa = I_ESP_M7y2z_D2x_aa+ABZ*I_ESP_L7yz_D2x_aa;
    Double I_ESP_L6y2z_F2xz_aa = I_ESP_M6y3z_D2x_aa+ABZ*I_ESP_L6y2z_D2x_aa;
    Double I_ESP_L5y3z_F2xz_aa = I_ESP_M5y4z_D2x_aa+ABZ*I_ESP_L5y3z_D2x_aa;
    Double I_ESP_L4y4z_F2xz_aa = I_ESP_M4y5z_D2x_aa+ABZ*I_ESP_L4y4z_D2x_aa;
    Double I_ESP_L3y5z_F2xz_aa = I_ESP_M3y6z_D2x_aa+ABZ*I_ESP_L3y5z_D2x_aa;
    Double I_ESP_L2y6z_F2xz_aa = I_ESP_M2y7z_D2x_aa+ABZ*I_ESP_L2y6z_D2x_aa;
    Double I_ESP_Ly7z_F2xz_aa = I_ESP_My8z_D2x_aa+ABZ*I_ESP_Ly7z_D2x_aa;
    Double I_ESP_L8z_F2xz_aa = I_ESP_M9z_D2x_aa+ABZ*I_ESP_L8z_D2x_aa;
    Double I_ESP_L7xz_Fx2y_aa = I_ESP_M8xz_D2y_aa+ABX*I_ESP_L7xz_D2y_aa;
    Double I_ESP_L6xyz_Fx2y_aa = I_ESP_M7xyz_D2y_aa+ABX*I_ESP_L6xyz_D2y_aa;
    Double I_ESP_L6x2z_Fx2y_aa = I_ESP_M7x2z_D2y_aa+ABX*I_ESP_L6x2z_D2y_aa;
    Double I_ESP_L5x2yz_Fx2y_aa = I_ESP_M6x2yz_D2y_aa+ABX*I_ESP_L5x2yz_D2y_aa;
    Double I_ESP_L5xy2z_Fx2y_aa = I_ESP_M6xy2z_D2y_aa+ABX*I_ESP_L5xy2z_D2y_aa;
    Double I_ESP_L5x3z_Fx2y_aa = I_ESP_M6x3z_D2y_aa+ABX*I_ESP_L5x3z_D2y_aa;
    Double I_ESP_L4x3yz_Fx2y_aa = I_ESP_M5x3yz_D2y_aa+ABX*I_ESP_L4x3yz_D2y_aa;
    Double I_ESP_L4x2y2z_Fx2y_aa = I_ESP_M5x2y2z_D2y_aa+ABX*I_ESP_L4x2y2z_D2y_aa;
    Double I_ESP_L4xy3z_Fx2y_aa = I_ESP_M5xy3z_D2y_aa+ABX*I_ESP_L4xy3z_D2y_aa;
    Double I_ESP_L4x4z_Fx2y_aa = I_ESP_M5x4z_D2y_aa+ABX*I_ESP_L4x4z_D2y_aa;
    Double I_ESP_L3x4yz_Fx2y_aa = I_ESP_M4x4yz_D2y_aa+ABX*I_ESP_L3x4yz_D2y_aa;
    Double I_ESP_L3x3y2z_Fx2y_aa = I_ESP_M4x3y2z_D2y_aa+ABX*I_ESP_L3x3y2z_D2y_aa;
    Double I_ESP_L3x2y3z_Fx2y_aa = I_ESP_M4x2y3z_D2y_aa+ABX*I_ESP_L3x2y3z_D2y_aa;
    Double I_ESP_L3xy4z_Fx2y_aa = I_ESP_M4xy4z_D2y_aa+ABX*I_ESP_L3xy4z_D2y_aa;
    Double I_ESP_L3x5z_Fx2y_aa = I_ESP_M4x5z_D2y_aa+ABX*I_ESP_L3x5z_D2y_aa;
    Double I_ESP_L2x5yz_Fx2y_aa = I_ESP_M3x5yz_D2y_aa+ABX*I_ESP_L2x5yz_D2y_aa;
    Double I_ESP_L2x4y2z_Fx2y_aa = I_ESP_M3x4y2z_D2y_aa+ABX*I_ESP_L2x4y2z_D2y_aa;
    Double I_ESP_L2x3y3z_Fx2y_aa = I_ESP_M3x3y3z_D2y_aa+ABX*I_ESP_L2x3y3z_D2y_aa;
    Double I_ESP_L2x2y4z_Fx2y_aa = I_ESP_M3x2y4z_D2y_aa+ABX*I_ESP_L2x2y4z_D2y_aa;
    Double I_ESP_L2xy5z_Fx2y_aa = I_ESP_M3xy5z_D2y_aa+ABX*I_ESP_L2xy5z_D2y_aa;
    Double I_ESP_L2x6z_Fx2y_aa = I_ESP_M3x6z_D2y_aa+ABX*I_ESP_L2x6z_D2y_aa;
    Double I_ESP_Lx6yz_Fx2y_aa = I_ESP_M2x6yz_D2y_aa+ABX*I_ESP_Lx6yz_D2y_aa;
    Double I_ESP_Lx5y2z_Fx2y_aa = I_ESP_M2x5y2z_D2y_aa+ABX*I_ESP_Lx5y2z_D2y_aa;
    Double I_ESP_Lx4y3z_Fx2y_aa = I_ESP_M2x4y3z_D2y_aa+ABX*I_ESP_Lx4y3z_D2y_aa;
    Double I_ESP_Lx3y4z_Fx2y_aa = I_ESP_M2x3y4z_D2y_aa+ABX*I_ESP_Lx3y4z_D2y_aa;
    Double I_ESP_Lx2y5z_Fx2y_aa = I_ESP_M2x2y5z_D2y_aa+ABX*I_ESP_Lx2y5z_D2y_aa;
    Double I_ESP_Lxy6z_Fx2y_aa = I_ESP_M2xy6z_D2y_aa+ABX*I_ESP_Lxy6z_D2y_aa;
    Double I_ESP_Lx7z_Fx2y_aa = I_ESP_M2x7z_D2y_aa+ABX*I_ESP_Lx7z_D2y_aa;
    Double I_ESP_L7yz_Fx2y_aa = I_ESP_Mx7yz_D2y_aa+ABX*I_ESP_L7yz_D2y_aa;
    Double I_ESP_L6y2z_Fx2y_aa = I_ESP_Mx6y2z_D2y_aa+ABX*I_ESP_L6y2z_D2y_aa;
    Double I_ESP_L5y3z_Fx2y_aa = I_ESP_Mx5y3z_D2y_aa+ABX*I_ESP_L5y3z_D2y_aa;
    Double I_ESP_L4y4z_Fx2y_aa = I_ESP_Mx4y4z_D2y_aa+ABX*I_ESP_L4y4z_D2y_aa;
    Double I_ESP_L3y5z_Fx2y_aa = I_ESP_Mx3y5z_D2y_aa+ABX*I_ESP_L3y5z_D2y_aa;
    Double I_ESP_L2y6z_Fx2y_aa = I_ESP_Mx2y6z_D2y_aa+ABX*I_ESP_L2y6z_D2y_aa;
    Double I_ESP_Ly7z_Fx2y_aa = I_ESP_Mxy7z_D2y_aa+ABX*I_ESP_Ly7z_D2y_aa;
    Double I_ESP_L8z_Fx2y_aa = I_ESP_Mx8z_D2y_aa+ABX*I_ESP_L8z_D2y_aa;
    Double I_ESP_L7xy_Fx2z_aa = I_ESP_M8xy_D2z_aa+ABX*I_ESP_L7xy_D2z_aa;
    Double I_ESP_L6x2y_Fx2z_aa = I_ESP_M7x2y_D2z_aa+ABX*I_ESP_L6x2y_D2z_aa;
    Double I_ESP_L6xyz_Fx2z_aa = I_ESP_M7xyz_D2z_aa+ABX*I_ESP_L6xyz_D2z_aa;
    Double I_ESP_L5x3y_Fx2z_aa = I_ESP_M6x3y_D2z_aa+ABX*I_ESP_L5x3y_D2z_aa;
    Double I_ESP_L5x2yz_Fx2z_aa = I_ESP_M6x2yz_D2z_aa+ABX*I_ESP_L5x2yz_D2z_aa;
    Double I_ESP_L5xy2z_Fx2z_aa = I_ESP_M6xy2z_D2z_aa+ABX*I_ESP_L5xy2z_D2z_aa;
    Double I_ESP_L4x4y_Fx2z_aa = I_ESP_M5x4y_D2z_aa+ABX*I_ESP_L4x4y_D2z_aa;
    Double I_ESP_L4x3yz_Fx2z_aa = I_ESP_M5x3yz_D2z_aa+ABX*I_ESP_L4x3yz_D2z_aa;
    Double I_ESP_L4x2y2z_Fx2z_aa = I_ESP_M5x2y2z_D2z_aa+ABX*I_ESP_L4x2y2z_D2z_aa;
    Double I_ESP_L4xy3z_Fx2z_aa = I_ESP_M5xy3z_D2z_aa+ABX*I_ESP_L4xy3z_D2z_aa;
    Double I_ESP_L3x5y_Fx2z_aa = I_ESP_M4x5y_D2z_aa+ABX*I_ESP_L3x5y_D2z_aa;
    Double I_ESP_L3x4yz_Fx2z_aa = I_ESP_M4x4yz_D2z_aa+ABX*I_ESP_L3x4yz_D2z_aa;
    Double I_ESP_L3x3y2z_Fx2z_aa = I_ESP_M4x3y2z_D2z_aa+ABX*I_ESP_L3x3y2z_D2z_aa;
    Double I_ESP_L3x2y3z_Fx2z_aa = I_ESP_M4x2y3z_D2z_aa+ABX*I_ESP_L3x2y3z_D2z_aa;
    Double I_ESP_L3xy4z_Fx2z_aa = I_ESP_M4xy4z_D2z_aa+ABX*I_ESP_L3xy4z_D2z_aa;
    Double I_ESP_L2x6y_Fx2z_aa = I_ESP_M3x6y_D2z_aa+ABX*I_ESP_L2x6y_D2z_aa;
    Double I_ESP_L2x5yz_Fx2z_aa = I_ESP_M3x5yz_D2z_aa+ABX*I_ESP_L2x5yz_D2z_aa;
    Double I_ESP_L2x4y2z_Fx2z_aa = I_ESP_M3x4y2z_D2z_aa+ABX*I_ESP_L2x4y2z_D2z_aa;
    Double I_ESP_L2x3y3z_Fx2z_aa = I_ESP_M3x3y3z_D2z_aa+ABX*I_ESP_L2x3y3z_D2z_aa;
    Double I_ESP_L2x2y4z_Fx2z_aa = I_ESP_M3x2y4z_D2z_aa+ABX*I_ESP_L2x2y4z_D2z_aa;
    Double I_ESP_L2xy5z_Fx2z_aa = I_ESP_M3xy5z_D2z_aa+ABX*I_ESP_L2xy5z_D2z_aa;
    Double I_ESP_Lx7y_Fx2z_aa = I_ESP_M2x7y_D2z_aa+ABX*I_ESP_Lx7y_D2z_aa;
    Double I_ESP_Lx6yz_Fx2z_aa = I_ESP_M2x6yz_D2z_aa+ABX*I_ESP_Lx6yz_D2z_aa;
    Double I_ESP_Lx5y2z_Fx2z_aa = I_ESP_M2x5y2z_D2z_aa+ABX*I_ESP_Lx5y2z_D2z_aa;
    Double I_ESP_Lx4y3z_Fx2z_aa = I_ESP_M2x4y3z_D2z_aa+ABX*I_ESP_Lx4y3z_D2z_aa;
    Double I_ESP_Lx3y4z_Fx2z_aa = I_ESP_M2x3y4z_D2z_aa+ABX*I_ESP_Lx3y4z_D2z_aa;
    Double I_ESP_Lx2y5z_Fx2z_aa = I_ESP_M2x2y5z_D2z_aa+ABX*I_ESP_Lx2y5z_D2z_aa;
    Double I_ESP_Lxy6z_Fx2z_aa = I_ESP_M2xy6z_D2z_aa+ABX*I_ESP_Lxy6z_D2z_aa;
    Double I_ESP_L8y_Fx2z_aa = I_ESP_Mx8y_D2z_aa+ABX*I_ESP_L8y_D2z_aa;
    Double I_ESP_L7yz_Fx2z_aa = I_ESP_Mx7yz_D2z_aa+ABX*I_ESP_L7yz_D2z_aa;
    Double I_ESP_L6y2z_Fx2z_aa = I_ESP_Mx6y2z_D2z_aa+ABX*I_ESP_L6y2z_D2z_aa;
    Double I_ESP_L5y3z_Fx2z_aa = I_ESP_Mx5y3z_D2z_aa+ABX*I_ESP_L5y3z_D2z_aa;
    Double I_ESP_L4y4z_Fx2z_aa = I_ESP_Mx4y4z_D2z_aa+ABX*I_ESP_L4y4z_D2z_aa;
    Double I_ESP_L3y5z_Fx2z_aa = I_ESP_Mx3y5z_D2z_aa+ABX*I_ESP_L3y5z_D2z_aa;
    Double I_ESP_L2y6z_Fx2z_aa = I_ESP_Mx2y6z_D2z_aa+ABX*I_ESP_L2y6z_D2z_aa;
    Double I_ESP_Ly7z_Fx2z_aa = I_ESP_Mxy7z_D2z_aa+ABX*I_ESP_Ly7z_D2z_aa;
    Double I_ESP_L8x_F3y_aa = I_ESP_M8xy_D2y_aa+ABY*I_ESP_L8x_D2y_aa;
    Double I_ESP_L7xy_F3y_aa = I_ESP_M7x2y_D2y_aa+ABY*I_ESP_L7xy_D2y_aa;
    Double I_ESP_L7xz_F3y_aa = I_ESP_M7xyz_D2y_aa+ABY*I_ESP_L7xz_D2y_aa;
    Double I_ESP_L6x2y_F3y_aa = I_ESP_M6x3y_D2y_aa+ABY*I_ESP_L6x2y_D2y_aa;
    Double I_ESP_L6xyz_F3y_aa = I_ESP_M6x2yz_D2y_aa+ABY*I_ESP_L6xyz_D2y_aa;
    Double I_ESP_L6x2z_F3y_aa = I_ESP_M6xy2z_D2y_aa+ABY*I_ESP_L6x2z_D2y_aa;
    Double I_ESP_L5x3y_F3y_aa = I_ESP_M5x4y_D2y_aa+ABY*I_ESP_L5x3y_D2y_aa;
    Double I_ESP_L5x2yz_F3y_aa = I_ESP_M5x3yz_D2y_aa+ABY*I_ESP_L5x2yz_D2y_aa;
    Double I_ESP_L5xy2z_F3y_aa = I_ESP_M5x2y2z_D2y_aa+ABY*I_ESP_L5xy2z_D2y_aa;
    Double I_ESP_L5x3z_F3y_aa = I_ESP_M5xy3z_D2y_aa+ABY*I_ESP_L5x3z_D2y_aa;
    Double I_ESP_L4x4y_F3y_aa = I_ESP_M4x5y_D2y_aa+ABY*I_ESP_L4x4y_D2y_aa;
    Double I_ESP_L4x3yz_F3y_aa = I_ESP_M4x4yz_D2y_aa+ABY*I_ESP_L4x3yz_D2y_aa;
    Double I_ESP_L4x2y2z_F3y_aa = I_ESP_M4x3y2z_D2y_aa+ABY*I_ESP_L4x2y2z_D2y_aa;
    Double I_ESP_L4xy3z_F3y_aa = I_ESP_M4x2y3z_D2y_aa+ABY*I_ESP_L4xy3z_D2y_aa;
    Double I_ESP_L4x4z_F3y_aa = I_ESP_M4xy4z_D2y_aa+ABY*I_ESP_L4x4z_D2y_aa;
    Double I_ESP_L3x5y_F3y_aa = I_ESP_M3x6y_D2y_aa+ABY*I_ESP_L3x5y_D2y_aa;
    Double I_ESP_L3x4yz_F3y_aa = I_ESP_M3x5yz_D2y_aa+ABY*I_ESP_L3x4yz_D2y_aa;
    Double I_ESP_L3x3y2z_F3y_aa = I_ESP_M3x4y2z_D2y_aa+ABY*I_ESP_L3x3y2z_D2y_aa;
    Double I_ESP_L3x2y3z_F3y_aa = I_ESP_M3x3y3z_D2y_aa+ABY*I_ESP_L3x2y3z_D2y_aa;
    Double I_ESP_L3xy4z_F3y_aa = I_ESP_M3x2y4z_D2y_aa+ABY*I_ESP_L3xy4z_D2y_aa;
    Double I_ESP_L3x5z_F3y_aa = I_ESP_M3xy5z_D2y_aa+ABY*I_ESP_L3x5z_D2y_aa;
    Double I_ESP_L2x6y_F3y_aa = I_ESP_M2x7y_D2y_aa+ABY*I_ESP_L2x6y_D2y_aa;
    Double I_ESP_L2x5yz_F3y_aa = I_ESP_M2x6yz_D2y_aa+ABY*I_ESP_L2x5yz_D2y_aa;
    Double I_ESP_L2x4y2z_F3y_aa = I_ESP_M2x5y2z_D2y_aa+ABY*I_ESP_L2x4y2z_D2y_aa;
    Double I_ESP_L2x3y3z_F3y_aa = I_ESP_M2x4y3z_D2y_aa+ABY*I_ESP_L2x3y3z_D2y_aa;
    Double I_ESP_L2x2y4z_F3y_aa = I_ESP_M2x3y4z_D2y_aa+ABY*I_ESP_L2x2y4z_D2y_aa;
    Double I_ESP_L2xy5z_F3y_aa = I_ESP_M2x2y5z_D2y_aa+ABY*I_ESP_L2xy5z_D2y_aa;
    Double I_ESP_L2x6z_F3y_aa = I_ESP_M2xy6z_D2y_aa+ABY*I_ESP_L2x6z_D2y_aa;
    Double I_ESP_Lx7y_F3y_aa = I_ESP_Mx8y_D2y_aa+ABY*I_ESP_Lx7y_D2y_aa;
    Double I_ESP_Lx6yz_F3y_aa = I_ESP_Mx7yz_D2y_aa+ABY*I_ESP_Lx6yz_D2y_aa;
    Double I_ESP_Lx5y2z_F3y_aa = I_ESP_Mx6y2z_D2y_aa+ABY*I_ESP_Lx5y2z_D2y_aa;
    Double I_ESP_Lx4y3z_F3y_aa = I_ESP_Mx5y3z_D2y_aa+ABY*I_ESP_Lx4y3z_D2y_aa;
    Double I_ESP_Lx3y4z_F3y_aa = I_ESP_Mx4y4z_D2y_aa+ABY*I_ESP_Lx3y4z_D2y_aa;
    Double I_ESP_Lx2y5z_F3y_aa = I_ESP_Mx3y5z_D2y_aa+ABY*I_ESP_Lx2y5z_D2y_aa;
    Double I_ESP_Lxy6z_F3y_aa = I_ESP_Mx2y6z_D2y_aa+ABY*I_ESP_Lxy6z_D2y_aa;
    Double I_ESP_Lx7z_F3y_aa = I_ESP_Mxy7z_D2y_aa+ABY*I_ESP_Lx7z_D2y_aa;
    Double I_ESP_L8y_F3y_aa = I_ESP_M9y_D2y_aa+ABY*I_ESP_L8y_D2y_aa;
    Double I_ESP_L7yz_F3y_aa = I_ESP_M8yz_D2y_aa+ABY*I_ESP_L7yz_D2y_aa;
    Double I_ESP_L6y2z_F3y_aa = I_ESP_M7y2z_D2y_aa+ABY*I_ESP_L6y2z_D2y_aa;
    Double I_ESP_L5y3z_F3y_aa = I_ESP_M6y3z_D2y_aa+ABY*I_ESP_L5y3z_D2y_aa;
    Double I_ESP_L4y4z_F3y_aa = I_ESP_M5y4z_D2y_aa+ABY*I_ESP_L4y4z_D2y_aa;
    Double I_ESP_L3y5z_F3y_aa = I_ESP_M4y5z_D2y_aa+ABY*I_ESP_L3y5z_D2y_aa;
    Double I_ESP_L2y6z_F3y_aa = I_ESP_M3y6z_D2y_aa+ABY*I_ESP_L2y6z_D2y_aa;
    Double I_ESP_Ly7z_F3y_aa = I_ESP_M2y7z_D2y_aa+ABY*I_ESP_Ly7z_D2y_aa;
    Double I_ESP_L8z_F3y_aa = I_ESP_My8z_D2y_aa+ABY*I_ESP_L8z_D2y_aa;
    Double I_ESP_L7xz_F2yz_aa = I_ESP_M7x2z_D2y_aa+ABZ*I_ESP_L7xz_D2y_aa;
    Double I_ESP_L6xyz_F2yz_aa = I_ESP_M6xy2z_D2y_aa+ABZ*I_ESP_L6xyz_D2y_aa;
    Double I_ESP_L6x2z_F2yz_aa = I_ESP_M6x3z_D2y_aa+ABZ*I_ESP_L6x2z_D2y_aa;
    Double I_ESP_L5x2yz_F2yz_aa = I_ESP_M5x2y2z_D2y_aa+ABZ*I_ESP_L5x2yz_D2y_aa;
    Double I_ESP_L5xy2z_F2yz_aa = I_ESP_M5xy3z_D2y_aa+ABZ*I_ESP_L5xy2z_D2y_aa;
    Double I_ESP_L5x3z_F2yz_aa = I_ESP_M5x4z_D2y_aa+ABZ*I_ESP_L5x3z_D2y_aa;
    Double I_ESP_L4x3yz_F2yz_aa = I_ESP_M4x3y2z_D2y_aa+ABZ*I_ESP_L4x3yz_D2y_aa;
    Double I_ESP_L4x2y2z_F2yz_aa = I_ESP_M4x2y3z_D2y_aa+ABZ*I_ESP_L4x2y2z_D2y_aa;
    Double I_ESP_L4xy3z_F2yz_aa = I_ESP_M4xy4z_D2y_aa+ABZ*I_ESP_L4xy3z_D2y_aa;
    Double I_ESP_L4x4z_F2yz_aa = I_ESP_M4x5z_D2y_aa+ABZ*I_ESP_L4x4z_D2y_aa;
    Double I_ESP_L3x4yz_F2yz_aa = I_ESP_M3x4y2z_D2y_aa+ABZ*I_ESP_L3x4yz_D2y_aa;
    Double I_ESP_L3x3y2z_F2yz_aa = I_ESP_M3x3y3z_D2y_aa+ABZ*I_ESP_L3x3y2z_D2y_aa;
    Double I_ESP_L3x2y3z_F2yz_aa = I_ESP_M3x2y4z_D2y_aa+ABZ*I_ESP_L3x2y3z_D2y_aa;
    Double I_ESP_L3xy4z_F2yz_aa = I_ESP_M3xy5z_D2y_aa+ABZ*I_ESP_L3xy4z_D2y_aa;
    Double I_ESP_L3x5z_F2yz_aa = I_ESP_M3x6z_D2y_aa+ABZ*I_ESP_L3x5z_D2y_aa;
    Double I_ESP_L2x5yz_F2yz_aa = I_ESP_M2x5y2z_D2y_aa+ABZ*I_ESP_L2x5yz_D2y_aa;
    Double I_ESP_L2x4y2z_F2yz_aa = I_ESP_M2x4y3z_D2y_aa+ABZ*I_ESP_L2x4y2z_D2y_aa;
    Double I_ESP_L2x3y3z_F2yz_aa = I_ESP_M2x3y4z_D2y_aa+ABZ*I_ESP_L2x3y3z_D2y_aa;
    Double I_ESP_L2x2y4z_F2yz_aa = I_ESP_M2x2y5z_D2y_aa+ABZ*I_ESP_L2x2y4z_D2y_aa;
    Double I_ESP_L2xy5z_F2yz_aa = I_ESP_M2xy6z_D2y_aa+ABZ*I_ESP_L2xy5z_D2y_aa;
    Double I_ESP_L2x6z_F2yz_aa = I_ESP_M2x7z_D2y_aa+ABZ*I_ESP_L2x6z_D2y_aa;
    Double I_ESP_Lx6yz_F2yz_aa = I_ESP_Mx6y2z_D2y_aa+ABZ*I_ESP_Lx6yz_D2y_aa;
    Double I_ESP_Lx5y2z_F2yz_aa = I_ESP_Mx5y3z_D2y_aa+ABZ*I_ESP_Lx5y2z_D2y_aa;
    Double I_ESP_Lx4y3z_F2yz_aa = I_ESP_Mx4y4z_D2y_aa+ABZ*I_ESP_Lx4y3z_D2y_aa;
    Double I_ESP_Lx3y4z_F2yz_aa = I_ESP_Mx3y5z_D2y_aa+ABZ*I_ESP_Lx3y4z_D2y_aa;
    Double I_ESP_Lx2y5z_F2yz_aa = I_ESP_Mx2y6z_D2y_aa+ABZ*I_ESP_Lx2y5z_D2y_aa;
    Double I_ESP_Lxy6z_F2yz_aa = I_ESP_Mxy7z_D2y_aa+ABZ*I_ESP_Lxy6z_D2y_aa;
    Double I_ESP_Lx7z_F2yz_aa = I_ESP_Mx8z_D2y_aa+ABZ*I_ESP_Lx7z_D2y_aa;
    Double I_ESP_L7yz_F2yz_aa = I_ESP_M7y2z_D2y_aa+ABZ*I_ESP_L7yz_D2y_aa;
    Double I_ESP_L6y2z_F2yz_aa = I_ESP_M6y3z_D2y_aa+ABZ*I_ESP_L6y2z_D2y_aa;
    Double I_ESP_L5y3z_F2yz_aa = I_ESP_M5y4z_D2y_aa+ABZ*I_ESP_L5y3z_D2y_aa;
    Double I_ESP_L4y4z_F2yz_aa = I_ESP_M4y5z_D2y_aa+ABZ*I_ESP_L4y4z_D2y_aa;
    Double I_ESP_L3y5z_F2yz_aa = I_ESP_M3y6z_D2y_aa+ABZ*I_ESP_L3y5z_D2y_aa;
    Double I_ESP_L2y6z_F2yz_aa = I_ESP_M2y7z_D2y_aa+ABZ*I_ESP_L2y6z_D2y_aa;
    Double I_ESP_Ly7z_F2yz_aa = I_ESP_My8z_D2y_aa+ABZ*I_ESP_Ly7z_D2y_aa;
    Double I_ESP_L8z_F2yz_aa = I_ESP_M9z_D2y_aa+ABZ*I_ESP_L8z_D2y_aa;
    Double I_ESP_L8x_F3z_aa = I_ESP_M8xz_D2z_aa+ABZ*I_ESP_L8x_D2z_aa;
    Double I_ESP_L7xy_F3z_aa = I_ESP_M7xyz_D2z_aa+ABZ*I_ESP_L7xy_D2z_aa;
    Double I_ESP_L7xz_F3z_aa = I_ESP_M7x2z_D2z_aa+ABZ*I_ESP_L7xz_D2z_aa;
    Double I_ESP_L6x2y_F3z_aa = I_ESP_M6x2yz_D2z_aa+ABZ*I_ESP_L6x2y_D2z_aa;
    Double I_ESP_L6xyz_F3z_aa = I_ESP_M6xy2z_D2z_aa+ABZ*I_ESP_L6xyz_D2z_aa;
    Double I_ESP_L6x2z_F3z_aa = I_ESP_M6x3z_D2z_aa+ABZ*I_ESP_L6x2z_D2z_aa;
    Double I_ESP_L5x3y_F3z_aa = I_ESP_M5x3yz_D2z_aa+ABZ*I_ESP_L5x3y_D2z_aa;
    Double I_ESP_L5x2yz_F3z_aa = I_ESP_M5x2y2z_D2z_aa+ABZ*I_ESP_L5x2yz_D2z_aa;
    Double I_ESP_L5xy2z_F3z_aa = I_ESP_M5xy3z_D2z_aa+ABZ*I_ESP_L5xy2z_D2z_aa;
    Double I_ESP_L5x3z_F3z_aa = I_ESP_M5x4z_D2z_aa+ABZ*I_ESP_L5x3z_D2z_aa;
    Double I_ESP_L4x4y_F3z_aa = I_ESP_M4x4yz_D2z_aa+ABZ*I_ESP_L4x4y_D2z_aa;
    Double I_ESP_L4x3yz_F3z_aa = I_ESP_M4x3y2z_D2z_aa+ABZ*I_ESP_L4x3yz_D2z_aa;
    Double I_ESP_L4x2y2z_F3z_aa = I_ESP_M4x2y3z_D2z_aa+ABZ*I_ESP_L4x2y2z_D2z_aa;
    Double I_ESP_L4xy3z_F3z_aa = I_ESP_M4xy4z_D2z_aa+ABZ*I_ESP_L4xy3z_D2z_aa;
    Double I_ESP_L4x4z_F3z_aa = I_ESP_M4x5z_D2z_aa+ABZ*I_ESP_L4x4z_D2z_aa;
    Double I_ESP_L3x5y_F3z_aa = I_ESP_M3x5yz_D2z_aa+ABZ*I_ESP_L3x5y_D2z_aa;
    Double I_ESP_L3x4yz_F3z_aa = I_ESP_M3x4y2z_D2z_aa+ABZ*I_ESP_L3x4yz_D2z_aa;
    Double I_ESP_L3x3y2z_F3z_aa = I_ESP_M3x3y3z_D2z_aa+ABZ*I_ESP_L3x3y2z_D2z_aa;
    Double I_ESP_L3x2y3z_F3z_aa = I_ESP_M3x2y4z_D2z_aa+ABZ*I_ESP_L3x2y3z_D2z_aa;
    Double I_ESP_L3xy4z_F3z_aa = I_ESP_M3xy5z_D2z_aa+ABZ*I_ESP_L3xy4z_D2z_aa;
    Double I_ESP_L3x5z_F3z_aa = I_ESP_M3x6z_D2z_aa+ABZ*I_ESP_L3x5z_D2z_aa;
    Double I_ESP_L2x6y_F3z_aa = I_ESP_M2x6yz_D2z_aa+ABZ*I_ESP_L2x6y_D2z_aa;
    Double I_ESP_L2x5yz_F3z_aa = I_ESP_M2x5y2z_D2z_aa+ABZ*I_ESP_L2x5yz_D2z_aa;
    Double I_ESP_L2x4y2z_F3z_aa = I_ESP_M2x4y3z_D2z_aa+ABZ*I_ESP_L2x4y2z_D2z_aa;
    Double I_ESP_L2x3y3z_F3z_aa = I_ESP_M2x3y4z_D2z_aa+ABZ*I_ESP_L2x3y3z_D2z_aa;
    Double I_ESP_L2x2y4z_F3z_aa = I_ESP_M2x2y5z_D2z_aa+ABZ*I_ESP_L2x2y4z_D2z_aa;
    Double I_ESP_L2xy5z_F3z_aa = I_ESP_M2xy6z_D2z_aa+ABZ*I_ESP_L2xy5z_D2z_aa;
    Double I_ESP_L2x6z_F3z_aa = I_ESP_M2x7z_D2z_aa+ABZ*I_ESP_L2x6z_D2z_aa;
    Double I_ESP_Lx7y_F3z_aa = I_ESP_Mx7yz_D2z_aa+ABZ*I_ESP_Lx7y_D2z_aa;
    Double I_ESP_Lx6yz_F3z_aa = I_ESP_Mx6y2z_D2z_aa+ABZ*I_ESP_Lx6yz_D2z_aa;
    Double I_ESP_Lx5y2z_F3z_aa = I_ESP_Mx5y3z_D2z_aa+ABZ*I_ESP_Lx5y2z_D2z_aa;
    Double I_ESP_Lx4y3z_F3z_aa = I_ESP_Mx4y4z_D2z_aa+ABZ*I_ESP_Lx4y3z_D2z_aa;
    Double I_ESP_Lx3y4z_F3z_aa = I_ESP_Mx3y5z_D2z_aa+ABZ*I_ESP_Lx3y4z_D2z_aa;
    Double I_ESP_Lx2y5z_F3z_aa = I_ESP_Mx2y6z_D2z_aa+ABZ*I_ESP_Lx2y5z_D2z_aa;
    Double I_ESP_Lxy6z_F3z_aa = I_ESP_Mxy7z_D2z_aa+ABZ*I_ESP_Lxy6z_D2z_aa;
    Double I_ESP_Lx7z_F3z_aa = I_ESP_Mx8z_D2z_aa+ABZ*I_ESP_Lx7z_D2z_aa;
    Double I_ESP_L8y_F3z_aa = I_ESP_M8yz_D2z_aa+ABZ*I_ESP_L8y_D2z_aa;
    Double I_ESP_L7yz_F3z_aa = I_ESP_M7y2z_D2z_aa+ABZ*I_ESP_L7yz_D2z_aa;
    Double I_ESP_L6y2z_F3z_aa = I_ESP_M6y3z_D2z_aa+ABZ*I_ESP_L6y2z_D2z_aa;
    Double I_ESP_L5y3z_F3z_aa = I_ESP_M5y4z_D2z_aa+ABZ*I_ESP_L5y3z_D2z_aa;
    Double I_ESP_L4y4z_F3z_aa = I_ESP_M4y5z_D2z_aa+ABZ*I_ESP_L4y4z_D2z_aa;
    Double I_ESP_L3y5z_F3z_aa = I_ESP_M3y6z_D2z_aa+ABZ*I_ESP_L3y5z_D2z_aa;
    Double I_ESP_L2y6z_F3z_aa = I_ESP_M2y7z_D2z_aa+ABZ*I_ESP_L2y6z_D2z_aa;
    Double I_ESP_Ly7z_F3z_aa = I_ESP_My8z_D2z_aa+ABZ*I_ESP_Ly7z_D2z_aa;
    Double I_ESP_L8z_F3z_aa = I_ESP_M9z_D2z_aa+ABZ*I_ESP_L8z_D2z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_K_G_aa
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_F_aa
     * RHS shell quartet name: SQ_ESP_K_F_aa
     ************************************************************/
    Double I_ESP_K7x_G4x_aa = I_ESP_L8x_F3x_aa+ABX*I_ESP_K7x_F3x_aa;
    Double I_ESP_K6xy_G4x_aa = I_ESP_L7xy_F3x_aa+ABX*I_ESP_K6xy_F3x_aa;
    Double I_ESP_K6xz_G4x_aa = I_ESP_L7xz_F3x_aa+ABX*I_ESP_K6xz_F3x_aa;
    Double I_ESP_K5x2y_G4x_aa = I_ESP_L6x2y_F3x_aa+ABX*I_ESP_K5x2y_F3x_aa;
    Double I_ESP_K5xyz_G4x_aa = I_ESP_L6xyz_F3x_aa+ABX*I_ESP_K5xyz_F3x_aa;
    Double I_ESP_K5x2z_G4x_aa = I_ESP_L6x2z_F3x_aa+ABX*I_ESP_K5x2z_F3x_aa;
    Double I_ESP_K4x3y_G4x_aa = I_ESP_L5x3y_F3x_aa+ABX*I_ESP_K4x3y_F3x_aa;
    Double I_ESP_K4x2yz_G4x_aa = I_ESP_L5x2yz_F3x_aa+ABX*I_ESP_K4x2yz_F3x_aa;
    Double I_ESP_K4xy2z_G4x_aa = I_ESP_L5xy2z_F3x_aa+ABX*I_ESP_K4xy2z_F3x_aa;
    Double I_ESP_K4x3z_G4x_aa = I_ESP_L5x3z_F3x_aa+ABX*I_ESP_K4x3z_F3x_aa;
    Double I_ESP_K3x4y_G4x_aa = I_ESP_L4x4y_F3x_aa+ABX*I_ESP_K3x4y_F3x_aa;
    Double I_ESP_K3x3yz_G4x_aa = I_ESP_L4x3yz_F3x_aa+ABX*I_ESP_K3x3yz_F3x_aa;
    Double I_ESP_K3x2y2z_G4x_aa = I_ESP_L4x2y2z_F3x_aa+ABX*I_ESP_K3x2y2z_F3x_aa;
    Double I_ESP_K3xy3z_G4x_aa = I_ESP_L4xy3z_F3x_aa+ABX*I_ESP_K3xy3z_F3x_aa;
    Double I_ESP_K3x4z_G4x_aa = I_ESP_L4x4z_F3x_aa+ABX*I_ESP_K3x4z_F3x_aa;
    Double I_ESP_K2x5y_G4x_aa = I_ESP_L3x5y_F3x_aa+ABX*I_ESP_K2x5y_F3x_aa;
    Double I_ESP_K2x4yz_G4x_aa = I_ESP_L3x4yz_F3x_aa+ABX*I_ESP_K2x4yz_F3x_aa;
    Double I_ESP_K2x3y2z_G4x_aa = I_ESP_L3x3y2z_F3x_aa+ABX*I_ESP_K2x3y2z_F3x_aa;
    Double I_ESP_K2x2y3z_G4x_aa = I_ESP_L3x2y3z_F3x_aa+ABX*I_ESP_K2x2y3z_F3x_aa;
    Double I_ESP_K2xy4z_G4x_aa = I_ESP_L3xy4z_F3x_aa+ABX*I_ESP_K2xy4z_F3x_aa;
    Double I_ESP_K2x5z_G4x_aa = I_ESP_L3x5z_F3x_aa+ABX*I_ESP_K2x5z_F3x_aa;
    Double I_ESP_Kx6y_G4x_aa = I_ESP_L2x6y_F3x_aa+ABX*I_ESP_Kx6y_F3x_aa;
    Double I_ESP_Kx5yz_G4x_aa = I_ESP_L2x5yz_F3x_aa+ABX*I_ESP_Kx5yz_F3x_aa;
    Double I_ESP_Kx4y2z_G4x_aa = I_ESP_L2x4y2z_F3x_aa+ABX*I_ESP_Kx4y2z_F3x_aa;
    Double I_ESP_Kx3y3z_G4x_aa = I_ESP_L2x3y3z_F3x_aa+ABX*I_ESP_Kx3y3z_F3x_aa;
    Double I_ESP_Kx2y4z_G4x_aa = I_ESP_L2x2y4z_F3x_aa+ABX*I_ESP_Kx2y4z_F3x_aa;
    Double I_ESP_Kxy5z_G4x_aa = I_ESP_L2xy5z_F3x_aa+ABX*I_ESP_Kxy5z_F3x_aa;
    Double I_ESP_Kx6z_G4x_aa = I_ESP_L2x6z_F3x_aa+ABX*I_ESP_Kx6z_F3x_aa;
    Double I_ESP_K7y_G4x_aa = I_ESP_Lx7y_F3x_aa+ABX*I_ESP_K7y_F3x_aa;
    Double I_ESP_K6yz_G4x_aa = I_ESP_Lx6yz_F3x_aa+ABX*I_ESP_K6yz_F3x_aa;
    Double I_ESP_K5y2z_G4x_aa = I_ESP_Lx5y2z_F3x_aa+ABX*I_ESP_K5y2z_F3x_aa;
    Double I_ESP_K4y3z_G4x_aa = I_ESP_Lx4y3z_F3x_aa+ABX*I_ESP_K4y3z_F3x_aa;
    Double I_ESP_K3y4z_G4x_aa = I_ESP_Lx3y4z_F3x_aa+ABX*I_ESP_K3y4z_F3x_aa;
    Double I_ESP_K2y5z_G4x_aa = I_ESP_Lx2y5z_F3x_aa+ABX*I_ESP_K2y5z_F3x_aa;
    Double I_ESP_Ky6z_G4x_aa = I_ESP_Lxy6z_F3x_aa+ABX*I_ESP_Ky6z_F3x_aa;
    Double I_ESP_K7z_G4x_aa = I_ESP_Lx7z_F3x_aa+ABX*I_ESP_K7z_F3x_aa;
    Double I_ESP_K7x_G3xy_aa = I_ESP_L7xy_F3x_aa+ABY*I_ESP_K7x_F3x_aa;
    Double I_ESP_K6xy_G3xy_aa = I_ESP_L6x2y_F3x_aa+ABY*I_ESP_K6xy_F3x_aa;
    Double I_ESP_K6xz_G3xy_aa = I_ESP_L6xyz_F3x_aa+ABY*I_ESP_K6xz_F3x_aa;
    Double I_ESP_K5x2y_G3xy_aa = I_ESP_L5x3y_F3x_aa+ABY*I_ESP_K5x2y_F3x_aa;
    Double I_ESP_K5xyz_G3xy_aa = I_ESP_L5x2yz_F3x_aa+ABY*I_ESP_K5xyz_F3x_aa;
    Double I_ESP_K5x2z_G3xy_aa = I_ESP_L5xy2z_F3x_aa+ABY*I_ESP_K5x2z_F3x_aa;
    Double I_ESP_K4x3y_G3xy_aa = I_ESP_L4x4y_F3x_aa+ABY*I_ESP_K4x3y_F3x_aa;
    Double I_ESP_K4x2yz_G3xy_aa = I_ESP_L4x3yz_F3x_aa+ABY*I_ESP_K4x2yz_F3x_aa;
    Double I_ESP_K4xy2z_G3xy_aa = I_ESP_L4x2y2z_F3x_aa+ABY*I_ESP_K4xy2z_F3x_aa;
    Double I_ESP_K4x3z_G3xy_aa = I_ESP_L4xy3z_F3x_aa+ABY*I_ESP_K4x3z_F3x_aa;
    Double I_ESP_K3x4y_G3xy_aa = I_ESP_L3x5y_F3x_aa+ABY*I_ESP_K3x4y_F3x_aa;
    Double I_ESP_K3x3yz_G3xy_aa = I_ESP_L3x4yz_F3x_aa+ABY*I_ESP_K3x3yz_F3x_aa;
    Double I_ESP_K3x2y2z_G3xy_aa = I_ESP_L3x3y2z_F3x_aa+ABY*I_ESP_K3x2y2z_F3x_aa;
    Double I_ESP_K3xy3z_G3xy_aa = I_ESP_L3x2y3z_F3x_aa+ABY*I_ESP_K3xy3z_F3x_aa;
    Double I_ESP_K3x4z_G3xy_aa = I_ESP_L3xy4z_F3x_aa+ABY*I_ESP_K3x4z_F3x_aa;
    Double I_ESP_K2x5y_G3xy_aa = I_ESP_L2x6y_F3x_aa+ABY*I_ESP_K2x5y_F3x_aa;
    Double I_ESP_K2x4yz_G3xy_aa = I_ESP_L2x5yz_F3x_aa+ABY*I_ESP_K2x4yz_F3x_aa;
    Double I_ESP_K2x3y2z_G3xy_aa = I_ESP_L2x4y2z_F3x_aa+ABY*I_ESP_K2x3y2z_F3x_aa;
    Double I_ESP_K2x2y3z_G3xy_aa = I_ESP_L2x3y3z_F3x_aa+ABY*I_ESP_K2x2y3z_F3x_aa;
    Double I_ESP_K2xy4z_G3xy_aa = I_ESP_L2x2y4z_F3x_aa+ABY*I_ESP_K2xy4z_F3x_aa;
    Double I_ESP_K2x5z_G3xy_aa = I_ESP_L2xy5z_F3x_aa+ABY*I_ESP_K2x5z_F3x_aa;
    Double I_ESP_Kx6y_G3xy_aa = I_ESP_Lx7y_F3x_aa+ABY*I_ESP_Kx6y_F3x_aa;
    Double I_ESP_Kx5yz_G3xy_aa = I_ESP_Lx6yz_F3x_aa+ABY*I_ESP_Kx5yz_F3x_aa;
    Double I_ESP_Kx4y2z_G3xy_aa = I_ESP_Lx5y2z_F3x_aa+ABY*I_ESP_Kx4y2z_F3x_aa;
    Double I_ESP_Kx3y3z_G3xy_aa = I_ESP_Lx4y3z_F3x_aa+ABY*I_ESP_Kx3y3z_F3x_aa;
    Double I_ESP_Kx2y4z_G3xy_aa = I_ESP_Lx3y4z_F3x_aa+ABY*I_ESP_Kx2y4z_F3x_aa;
    Double I_ESP_Kxy5z_G3xy_aa = I_ESP_Lx2y5z_F3x_aa+ABY*I_ESP_Kxy5z_F3x_aa;
    Double I_ESP_Kx6z_G3xy_aa = I_ESP_Lxy6z_F3x_aa+ABY*I_ESP_Kx6z_F3x_aa;
    Double I_ESP_K7y_G3xy_aa = I_ESP_L8y_F3x_aa+ABY*I_ESP_K7y_F3x_aa;
    Double I_ESP_K6yz_G3xy_aa = I_ESP_L7yz_F3x_aa+ABY*I_ESP_K6yz_F3x_aa;
    Double I_ESP_K5y2z_G3xy_aa = I_ESP_L6y2z_F3x_aa+ABY*I_ESP_K5y2z_F3x_aa;
    Double I_ESP_K4y3z_G3xy_aa = I_ESP_L5y3z_F3x_aa+ABY*I_ESP_K4y3z_F3x_aa;
    Double I_ESP_K3y4z_G3xy_aa = I_ESP_L4y4z_F3x_aa+ABY*I_ESP_K3y4z_F3x_aa;
    Double I_ESP_K2y5z_G3xy_aa = I_ESP_L3y5z_F3x_aa+ABY*I_ESP_K2y5z_F3x_aa;
    Double I_ESP_Ky6z_G3xy_aa = I_ESP_L2y6z_F3x_aa+ABY*I_ESP_Ky6z_F3x_aa;
    Double I_ESP_K7z_G3xy_aa = I_ESP_Ly7z_F3x_aa+ABY*I_ESP_K7z_F3x_aa;
    Double I_ESP_K7x_G3xz_aa = I_ESP_L7xz_F3x_aa+ABZ*I_ESP_K7x_F3x_aa;
    Double I_ESP_K6xy_G3xz_aa = I_ESP_L6xyz_F3x_aa+ABZ*I_ESP_K6xy_F3x_aa;
    Double I_ESP_K6xz_G3xz_aa = I_ESP_L6x2z_F3x_aa+ABZ*I_ESP_K6xz_F3x_aa;
    Double I_ESP_K5x2y_G3xz_aa = I_ESP_L5x2yz_F3x_aa+ABZ*I_ESP_K5x2y_F3x_aa;
    Double I_ESP_K5xyz_G3xz_aa = I_ESP_L5xy2z_F3x_aa+ABZ*I_ESP_K5xyz_F3x_aa;
    Double I_ESP_K5x2z_G3xz_aa = I_ESP_L5x3z_F3x_aa+ABZ*I_ESP_K5x2z_F3x_aa;
    Double I_ESP_K4x3y_G3xz_aa = I_ESP_L4x3yz_F3x_aa+ABZ*I_ESP_K4x3y_F3x_aa;
    Double I_ESP_K4x2yz_G3xz_aa = I_ESP_L4x2y2z_F3x_aa+ABZ*I_ESP_K4x2yz_F3x_aa;
    Double I_ESP_K4xy2z_G3xz_aa = I_ESP_L4xy3z_F3x_aa+ABZ*I_ESP_K4xy2z_F3x_aa;
    Double I_ESP_K4x3z_G3xz_aa = I_ESP_L4x4z_F3x_aa+ABZ*I_ESP_K4x3z_F3x_aa;
    Double I_ESP_K3x4y_G3xz_aa = I_ESP_L3x4yz_F3x_aa+ABZ*I_ESP_K3x4y_F3x_aa;
    Double I_ESP_K3x3yz_G3xz_aa = I_ESP_L3x3y2z_F3x_aa+ABZ*I_ESP_K3x3yz_F3x_aa;
    Double I_ESP_K3x2y2z_G3xz_aa = I_ESP_L3x2y3z_F3x_aa+ABZ*I_ESP_K3x2y2z_F3x_aa;
    Double I_ESP_K3xy3z_G3xz_aa = I_ESP_L3xy4z_F3x_aa+ABZ*I_ESP_K3xy3z_F3x_aa;
    Double I_ESP_K3x4z_G3xz_aa = I_ESP_L3x5z_F3x_aa+ABZ*I_ESP_K3x4z_F3x_aa;
    Double I_ESP_K2x5y_G3xz_aa = I_ESP_L2x5yz_F3x_aa+ABZ*I_ESP_K2x5y_F3x_aa;
    Double I_ESP_K2x4yz_G3xz_aa = I_ESP_L2x4y2z_F3x_aa+ABZ*I_ESP_K2x4yz_F3x_aa;
    Double I_ESP_K2x3y2z_G3xz_aa = I_ESP_L2x3y3z_F3x_aa+ABZ*I_ESP_K2x3y2z_F3x_aa;
    Double I_ESP_K2x2y3z_G3xz_aa = I_ESP_L2x2y4z_F3x_aa+ABZ*I_ESP_K2x2y3z_F3x_aa;
    Double I_ESP_K2xy4z_G3xz_aa = I_ESP_L2xy5z_F3x_aa+ABZ*I_ESP_K2xy4z_F3x_aa;
    Double I_ESP_K2x5z_G3xz_aa = I_ESP_L2x6z_F3x_aa+ABZ*I_ESP_K2x5z_F3x_aa;
    Double I_ESP_Kx6y_G3xz_aa = I_ESP_Lx6yz_F3x_aa+ABZ*I_ESP_Kx6y_F3x_aa;
    Double I_ESP_Kx5yz_G3xz_aa = I_ESP_Lx5y2z_F3x_aa+ABZ*I_ESP_Kx5yz_F3x_aa;
    Double I_ESP_Kx4y2z_G3xz_aa = I_ESP_Lx4y3z_F3x_aa+ABZ*I_ESP_Kx4y2z_F3x_aa;
    Double I_ESP_Kx3y3z_G3xz_aa = I_ESP_Lx3y4z_F3x_aa+ABZ*I_ESP_Kx3y3z_F3x_aa;
    Double I_ESP_Kx2y4z_G3xz_aa = I_ESP_Lx2y5z_F3x_aa+ABZ*I_ESP_Kx2y4z_F3x_aa;
    Double I_ESP_Kxy5z_G3xz_aa = I_ESP_Lxy6z_F3x_aa+ABZ*I_ESP_Kxy5z_F3x_aa;
    Double I_ESP_Kx6z_G3xz_aa = I_ESP_Lx7z_F3x_aa+ABZ*I_ESP_Kx6z_F3x_aa;
    Double I_ESP_K7y_G3xz_aa = I_ESP_L7yz_F3x_aa+ABZ*I_ESP_K7y_F3x_aa;
    Double I_ESP_K6yz_G3xz_aa = I_ESP_L6y2z_F3x_aa+ABZ*I_ESP_K6yz_F3x_aa;
    Double I_ESP_K5y2z_G3xz_aa = I_ESP_L5y3z_F3x_aa+ABZ*I_ESP_K5y2z_F3x_aa;
    Double I_ESP_K4y3z_G3xz_aa = I_ESP_L4y4z_F3x_aa+ABZ*I_ESP_K4y3z_F3x_aa;
    Double I_ESP_K3y4z_G3xz_aa = I_ESP_L3y5z_F3x_aa+ABZ*I_ESP_K3y4z_F3x_aa;
    Double I_ESP_K2y5z_G3xz_aa = I_ESP_L2y6z_F3x_aa+ABZ*I_ESP_K2y5z_F3x_aa;
    Double I_ESP_Ky6z_G3xz_aa = I_ESP_Ly7z_F3x_aa+ABZ*I_ESP_Ky6z_F3x_aa;
    Double I_ESP_K7z_G3xz_aa = I_ESP_L8z_F3x_aa+ABZ*I_ESP_K7z_F3x_aa;
    Double I_ESP_K7x_G2x2y_aa = I_ESP_L7xy_F2xy_aa+ABY*I_ESP_K7x_F2xy_aa;
    Double I_ESP_K6xy_G2x2y_aa = I_ESP_L6x2y_F2xy_aa+ABY*I_ESP_K6xy_F2xy_aa;
    Double I_ESP_K6xz_G2x2y_aa = I_ESP_L6xyz_F2xy_aa+ABY*I_ESP_K6xz_F2xy_aa;
    Double I_ESP_K5x2y_G2x2y_aa = I_ESP_L5x3y_F2xy_aa+ABY*I_ESP_K5x2y_F2xy_aa;
    Double I_ESP_K5xyz_G2x2y_aa = I_ESP_L5x2yz_F2xy_aa+ABY*I_ESP_K5xyz_F2xy_aa;
    Double I_ESP_K5x2z_G2x2y_aa = I_ESP_L5xy2z_F2xy_aa+ABY*I_ESP_K5x2z_F2xy_aa;
    Double I_ESP_K4x3y_G2x2y_aa = I_ESP_L4x4y_F2xy_aa+ABY*I_ESP_K4x3y_F2xy_aa;
    Double I_ESP_K4x2yz_G2x2y_aa = I_ESP_L4x3yz_F2xy_aa+ABY*I_ESP_K4x2yz_F2xy_aa;
    Double I_ESP_K4xy2z_G2x2y_aa = I_ESP_L4x2y2z_F2xy_aa+ABY*I_ESP_K4xy2z_F2xy_aa;
    Double I_ESP_K4x3z_G2x2y_aa = I_ESP_L4xy3z_F2xy_aa+ABY*I_ESP_K4x3z_F2xy_aa;
    Double I_ESP_K3x4y_G2x2y_aa = I_ESP_L3x5y_F2xy_aa+ABY*I_ESP_K3x4y_F2xy_aa;
    Double I_ESP_K3x3yz_G2x2y_aa = I_ESP_L3x4yz_F2xy_aa+ABY*I_ESP_K3x3yz_F2xy_aa;
    Double I_ESP_K3x2y2z_G2x2y_aa = I_ESP_L3x3y2z_F2xy_aa+ABY*I_ESP_K3x2y2z_F2xy_aa;
    Double I_ESP_K3xy3z_G2x2y_aa = I_ESP_L3x2y3z_F2xy_aa+ABY*I_ESP_K3xy3z_F2xy_aa;
    Double I_ESP_K3x4z_G2x2y_aa = I_ESP_L3xy4z_F2xy_aa+ABY*I_ESP_K3x4z_F2xy_aa;
    Double I_ESP_K2x5y_G2x2y_aa = I_ESP_L2x6y_F2xy_aa+ABY*I_ESP_K2x5y_F2xy_aa;
    Double I_ESP_K2x4yz_G2x2y_aa = I_ESP_L2x5yz_F2xy_aa+ABY*I_ESP_K2x4yz_F2xy_aa;
    Double I_ESP_K2x3y2z_G2x2y_aa = I_ESP_L2x4y2z_F2xy_aa+ABY*I_ESP_K2x3y2z_F2xy_aa;
    Double I_ESP_K2x2y3z_G2x2y_aa = I_ESP_L2x3y3z_F2xy_aa+ABY*I_ESP_K2x2y3z_F2xy_aa;
    Double I_ESP_K2xy4z_G2x2y_aa = I_ESP_L2x2y4z_F2xy_aa+ABY*I_ESP_K2xy4z_F2xy_aa;
    Double I_ESP_K2x5z_G2x2y_aa = I_ESP_L2xy5z_F2xy_aa+ABY*I_ESP_K2x5z_F2xy_aa;
    Double I_ESP_Kx6y_G2x2y_aa = I_ESP_Lx7y_F2xy_aa+ABY*I_ESP_Kx6y_F2xy_aa;
    Double I_ESP_Kx5yz_G2x2y_aa = I_ESP_Lx6yz_F2xy_aa+ABY*I_ESP_Kx5yz_F2xy_aa;
    Double I_ESP_Kx4y2z_G2x2y_aa = I_ESP_Lx5y2z_F2xy_aa+ABY*I_ESP_Kx4y2z_F2xy_aa;
    Double I_ESP_Kx3y3z_G2x2y_aa = I_ESP_Lx4y3z_F2xy_aa+ABY*I_ESP_Kx3y3z_F2xy_aa;
    Double I_ESP_Kx2y4z_G2x2y_aa = I_ESP_Lx3y4z_F2xy_aa+ABY*I_ESP_Kx2y4z_F2xy_aa;
    Double I_ESP_Kxy5z_G2x2y_aa = I_ESP_Lx2y5z_F2xy_aa+ABY*I_ESP_Kxy5z_F2xy_aa;
    Double I_ESP_Kx6z_G2x2y_aa = I_ESP_Lxy6z_F2xy_aa+ABY*I_ESP_Kx6z_F2xy_aa;
    Double I_ESP_K7y_G2x2y_aa = I_ESP_L8y_F2xy_aa+ABY*I_ESP_K7y_F2xy_aa;
    Double I_ESP_K6yz_G2x2y_aa = I_ESP_L7yz_F2xy_aa+ABY*I_ESP_K6yz_F2xy_aa;
    Double I_ESP_K5y2z_G2x2y_aa = I_ESP_L6y2z_F2xy_aa+ABY*I_ESP_K5y2z_F2xy_aa;
    Double I_ESP_K4y3z_G2x2y_aa = I_ESP_L5y3z_F2xy_aa+ABY*I_ESP_K4y3z_F2xy_aa;
    Double I_ESP_K3y4z_G2x2y_aa = I_ESP_L4y4z_F2xy_aa+ABY*I_ESP_K3y4z_F2xy_aa;
    Double I_ESP_K2y5z_G2x2y_aa = I_ESP_L3y5z_F2xy_aa+ABY*I_ESP_K2y5z_F2xy_aa;
    Double I_ESP_Ky6z_G2x2y_aa = I_ESP_L2y6z_F2xy_aa+ABY*I_ESP_Ky6z_F2xy_aa;
    Double I_ESP_K7z_G2x2y_aa = I_ESP_Ly7z_F2xy_aa+ABY*I_ESP_K7z_F2xy_aa;
    Double I_ESP_K7x_G2xyz_aa = I_ESP_L7xz_F2xy_aa+ABZ*I_ESP_K7x_F2xy_aa;
    Double I_ESP_K6xy_G2xyz_aa = I_ESP_L6xyz_F2xy_aa+ABZ*I_ESP_K6xy_F2xy_aa;
    Double I_ESP_K6xz_G2xyz_aa = I_ESP_L6x2z_F2xy_aa+ABZ*I_ESP_K6xz_F2xy_aa;
    Double I_ESP_K5x2y_G2xyz_aa = I_ESP_L5x2yz_F2xy_aa+ABZ*I_ESP_K5x2y_F2xy_aa;
    Double I_ESP_K5xyz_G2xyz_aa = I_ESP_L5xy2z_F2xy_aa+ABZ*I_ESP_K5xyz_F2xy_aa;
    Double I_ESP_K5x2z_G2xyz_aa = I_ESP_L5x3z_F2xy_aa+ABZ*I_ESP_K5x2z_F2xy_aa;
    Double I_ESP_K4x3y_G2xyz_aa = I_ESP_L4x3yz_F2xy_aa+ABZ*I_ESP_K4x3y_F2xy_aa;
    Double I_ESP_K4x2yz_G2xyz_aa = I_ESP_L4x2y2z_F2xy_aa+ABZ*I_ESP_K4x2yz_F2xy_aa;
    Double I_ESP_K4xy2z_G2xyz_aa = I_ESP_L4xy3z_F2xy_aa+ABZ*I_ESP_K4xy2z_F2xy_aa;
    Double I_ESP_K4x3z_G2xyz_aa = I_ESP_L4x4z_F2xy_aa+ABZ*I_ESP_K4x3z_F2xy_aa;
    Double I_ESP_K3x4y_G2xyz_aa = I_ESP_L3x4yz_F2xy_aa+ABZ*I_ESP_K3x4y_F2xy_aa;
    Double I_ESP_K3x3yz_G2xyz_aa = I_ESP_L3x3y2z_F2xy_aa+ABZ*I_ESP_K3x3yz_F2xy_aa;
    Double I_ESP_K3x2y2z_G2xyz_aa = I_ESP_L3x2y3z_F2xy_aa+ABZ*I_ESP_K3x2y2z_F2xy_aa;
    Double I_ESP_K3xy3z_G2xyz_aa = I_ESP_L3xy4z_F2xy_aa+ABZ*I_ESP_K3xy3z_F2xy_aa;
    Double I_ESP_K3x4z_G2xyz_aa = I_ESP_L3x5z_F2xy_aa+ABZ*I_ESP_K3x4z_F2xy_aa;
    Double I_ESP_K2x5y_G2xyz_aa = I_ESP_L2x5yz_F2xy_aa+ABZ*I_ESP_K2x5y_F2xy_aa;
    Double I_ESP_K2x4yz_G2xyz_aa = I_ESP_L2x4y2z_F2xy_aa+ABZ*I_ESP_K2x4yz_F2xy_aa;
    Double I_ESP_K2x3y2z_G2xyz_aa = I_ESP_L2x3y3z_F2xy_aa+ABZ*I_ESP_K2x3y2z_F2xy_aa;
    Double I_ESP_K2x2y3z_G2xyz_aa = I_ESP_L2x2y4z_F2xy_aa+ABZ*I_ESP_K2x2y3z_F2xy_aa;
    Double I_ESP_K2xy4z_G2xyz_aa = I_ESP_L2xy5z_F2xy_aa+ABZ*I_ESP_K2xy4z_F2xy_aa;
    Double I_ESP_K2x5z_G2xyz_aa = I_ESP_L2x6z_F2xy_aa+ABZ*I_ESP_K2x5z_F2xy_aa;
    Double I_ESP_Kx6y_G2xyz_aa = I_ESP_Lx6yz_F2xy_aa+ABZ*I_ESP_Kx6y_F2xy_aa;
    Double I_ESP_Kx5yz_G2xyz_aa = I_ESP_Lx5y2z_F2xy_aa+ABZ*I_ESP_Kx5yz_F2xy_aa;
    Double I_ESP_Kx4y2z_G2xyz_aa = I_ESP_Lx4y3z_F2xy_aa+ABZ*I_ESP_Kx4y2z_F2xy_aa;
    Double I_ESP_Kx3y3z_G2xyz_aa = I_ESP_Lx3y4z_F2xy_aa+ABZ*I_ESP_Kx3y3z_F2xy_aa;
    Double I_ESP_Kx2y4z_G2xyz_aa = I_ESP_Lx2y5z_F2xy_aa+ABZ*I_ESP_Kx2y4z_F2xy_aa;
    Double I_ESP_Kxy5z_G2xyz_aa = I_ESP_Lxy6z_F2xy_aa+ABZ*I_ESP_Kxy5z_F2xy_aa;
    Double I_ESP_Kx6z_G2xyz_aa = I_ESP_Lx7z_F2xy_aa+ABZ*I_ESP_Kx6z_F2xy_aa;
    Double I_ESP_K7y_G2xyz_aa = I_ESP_L7yz_F2xy_aa+ABZ*I_ESP_K7y_F2xy_aa;
    Double I_ESP_K6yz_G2xyz_aa = I_ESP_L6y2z_F2xy_aa+ABZ*I_ESP_K6yz_F2xy_aa;
    Double I_ESP_K5y2z_G2xyz_aa = I_ESP_L5y3z_F2xy_aa+ABZ*I_ESP_K5y2z_F2xy_aa;
    Double I_ESP_K4y3z_G2xyz_aa = I_ESP_L4y4z_F2xy_aa+ABZ*I_ESP_K4y3z_F2xy_aa;
    Double I_ESP_K3y4z_G2xyz_aa = I_ESP_L3y5z_F2xy_aa+ABZ*I_ESP_K3y4z_F2xy_aa;
    Double I_ESP_K2y5z_G2xyz_aa = I_ESP_L2y6z_F2xy_aa+ABZ*I_ESP_K2y5z_F2xy_aa;
    Double I_ESP_Ky6z_G2xyz_aa = I_ESP_Ly7z_F2xy_aa+ABZ*I_ESP_Ky6z_F2xy_aa;
    Double I_ESP_K7z_G2xyz_aa = I_ESP_L8z_F2xy_aa+ABZ*I_ESP_K7z_F2xy_aa;
    Double I_ESP_K7x_G2x2z_aa = I_ESP_L7xz_F2xz_aa+ABZ*I_ESP_K7x_F2xz_aa;
    Double I_ESP_K6xy_G2x2z_aa = I_ESP_L6xyz_F2xz_aa+ABZ*I_ESP_K6xy_F2xz_aa;
    Double I_ESP_K6xz_G2x2z_aa = I_ESP_L6x2z_F2xz_aa+ABZ*I_ESP_K6xz_F2xz_aa;
    Double I_ESP_K5x2y_G2x2z_aa = I_ESP_L5x2yz_F2xz_aa+ABZ*I_ESP_K5x2y_F2xz_aa;
    Double I_ESP_K5xyz_G2x2z_aa = I_ESP_L5xy2z_F2xz_aa+ABZ*I_ESP_K5xyz_F2xz_aa;
    Double I_ESP_K5x2z_G2x2z_aa = I_ESP_L5x3z_F2xz_aa+ABZ*I_ESP_K5x2z_F2xz_aa;
    Double I_ESP_K4x3y_G2x2z_aa = I_ESP_L4x3yz_F2xz_aa+ABZ*I_ESP_K4x3y_F2xz_aa;
    Double I_ESP_K4x2yz_G2x2z_aa = I_ESP_L4x2y2z_F2xz_aa+ABZ*I_ESP_K4x2yz_F2xz_aa;
    Double I_ESP_K4xy2z_G2x2z_aa = I_ESP_L4xy3z_F2xz_aa+ABZ*I_ESP_K4xy2z_F2xz_aa;
    Double I_ESP_K4x3z_G2x2z_aa = I_ESP_L4x4z_F2xz_aa+ABZ*I_ESP_K4x3z_F2xz_aa;
    Double I_ESP_K3x4y_G2x2z_aa = I_ESP_L3x4yz_F2xz_aa+ABZ*I_ESP_K3x4y_F2xz_aa;
    Double I_ESP_K3x3yz_G2x2z_aa = I_ESP_L3x3y2z_F2xz_aa+ABZ*I_ESP_K3x3yz_F2xz_aa;
    Double I_ESP_K3x2y2z_G2x2z_aa = I_ESP_L3x2y3z_F2xz_aa+ABZ*I_ESP_K3x2y2z_F2xz_aa;
    Double I_ESP_K3xy3z_G2x2z_aa = I_ESP_L3xy4z_F2xz_aa+ABZ*I_ESP_K3xy3z_F2xz_aa;
    Double I_ESP_K3x4z_G2x2z_aa = I_ESP_L3x5z_F2xz_aa+ABZ*I_ESP_K3x4z_F2xz_aa;
    Double I_ESP_K2x5y_G2x2z_aa = I_ESP_L2x5yz_F2xz_aa+ABZ*I_ESP_K2x5y_F2xz_aa;
    Double I_ESP_K2x4yz_G2x2z_aa = I_ESP_L2x4y2z_F2xz_aa+ABZ*I_ESP_K2x4yz_F2xz_aa;
    Double I_ESP_K2x3y2z_G2x2z_aa = I_ESP_L2x3y3z_F2xz_aa+ABZ*I_ESP_K2x3y2z_F2xz_aa;
    Double I_ESP_K2x2y3z_G2x2z_aa = I_ESP_L2x2y4z_F2xz_aa+ABZ*I_ESP_K2x2y3z_F2xz_aa;
    Double I_ESP_K2xy4z_G2x2z_aa = I_ESP_L2xy5z_F2xz_aa+ABZ*I_ESP_K2xy4z_F2xz_aa;
    Double I_ESP_K2x5z_G2x2z_aa = I_ESP_L2x6z_F2xz_aa+ABZ*I_ESP_K2x5z_F2xz_aa;
    Double I_ESP_Kx6y_G2x2z_aa = I_ESP_Lx6yz_F2xz_aa+ABZ*I_ESP_Kx6y_F2xz_aa;
    Double I_ESP_Kx5yz_G2x2z_aa = I_ESP_Lx5y2z_F2xz_aa+ABZ*I_ESP_Kx5yz_F2xz_aa;
    Double I_ESP_Kx4y2z_G2x2z_aa = I_ESP_Lx4y3z_F2xz_aa+ABZ*I_ESP_Kx4y2z_F2xz_aa;
    Double I_ESP_Kx3y3z_G2x2z_aa = I_ESP_Lx3y4z_F2xz_aa+ABZ*I_ESP_Kx3y3z_F2xz_aa;
    Double I_ESP_Kx2y4z_G2x2z_aa = I_ESP_Lx2y5z_F2xz_aa+ABZ*I_ESP_Kx2y4z_F2xz_aa;
    Double I_ESP_Kxy5z_G2x2z_aa = I_ESP_Lxy6z_F2xz_aa+ABZ*I_ESP_Kxy5z_F2xz_aa;
    Double I_ESP_Kx6z_G2x2z_aa = I_ESP_Lx7z_F2xz_aa+ABZ*I_ESP_Kx6z_F2xz_aa;
    Double I_ESP_K7y_G2x2z_aa = I_ESP_L7yz_F2xz_aa+ABZ*I_ESP_K7y_F2xz_aa;
    Double I_ESP_K6yz_G2x2z_aa = I_ESP_L6y2z_F2xz_aa+ABZ*I_ESP_K6yz_F2xz_aa;
    Double I_ESP_K5y2z_G2x2z_aa = I_ESP_L5y3z_F2xz_aa+ABZ*I_ESP_K5y2z_F2xz_aa;
    Double I_ESP_K4y3z_G2x2z_aa = I_ESP_L4y4z_F2xz_aa+ABZ*I_ESP_K4y3z_F2xz_aa;
    Double I_ESP_K3y4z_G2x2z_aa = I_ESP_L3y5z_F2xz_aa+ABZ*I_ESP_K3y4z_F2xz_aa;
    Double I_ESP_K2y5z_G2x2z_aa = I_ESP_L2y6z_F2xz_aa+ABZ*I_ESP_K2y5z_F2xz_aa;
    Double I_ESP_Ky6z_G2x2z_aa = I_ESP_Ly7z_F2xz_aa+ABZ*I_ESP_Ky6z_F2xz_aa;
    Double I_ESP_K7z_G2x2z_aa = I_ESP_L8z_F2xz_aa+ABZ*I_ESP_K7z_F2xz_aa;
    Double I_ESP_K7x_Gx3y_aa = I_ESP_L8x_F3y_aa+ABX*I_ESP_K7x_F3y_aa;
    Double I_ESP_K6xy_Gx3y_aa = I_ESP_L7xy_F3y_aa+ABX*I_ESP_K6xy_F3y_aa;
    Double I_ESP_K6xz_Gx3y_aa = I_ESP_L7xz_F3y_aa+ABX*I_ESP_K6xz_F3y_aa;
    Double I_ESP_K5x2y_Gx3y_aa = I_ESP_L6x2y_F3y_aa+ABX*I_ESP_K5x2y_F3y_aa;
    Double I_ESP_K5xyz_Gx3y_aa = I_ESP_L6xyz_F3y_aa+ABX*I_ESP_K5xyz_F3y_aa;
    Double I_ESP_K5x2z_Gx3y_aa = I_ESP_L6x2z_F3y_aa+ABX*I_ESP_K5x2z_F3y_aa;
    Double I_ESP_K4x3y_Gx3y_aa = I_ESP_L5x3y_F3y_aa+ABX*I_ESP_K4x3y_F3y_aa;
    Double I_ESP_K4x2yz_Gx3y_aa = I_ESP_L5x2yz_F3y_aa+ABX*I_ESP_K4x2yz_F3y_aa;
    Double I_ESP_K4xy2z_Gx3y_aa = I_ESP_L5xy2z_F3y_aa+ABX*I_ESP_K4xy2z_F3y_aa;
    Double I_ESP_K4x3z_Gx3y_aa = I_ESP_L5x3z_F3y_aa+ABX*I_ESP_K4x3z_F3y_aa;
    Double I_ESP_K3x4y_Gx3y_aa = I_ESP_L4x4y_F3y_aa+ABX*I_ESP_K3x4y_F3y_aa;
    Double I_ESP_K3x3yz_Gx3y_aa = I_ESP_L4x3yz_F3y_aa+ABX*I_ESP_K3x3yz_F3y_aa;
    Double I_ESP_K3x2y2z_Gx3y_aa = I_ESP_L4x2y2z_F3y_aa+ABX*I_ESP_K3x2y2z_F3y_aa;
    Double I_ESP_K3xy3z_Gx3y_aa = I_ESP_L4xy3z_F3y_aa+ABX*I_ESP_K3xy3z_F3y_aa;
    Double I_ESP_K3x4z_Gx3y_aa = I_ESP_L4x4z_F3y_aa+ABX*I_ESP_K3x4z_F3y_aa;
    Double I_ESP_K2x5y_Gx3y_aa = I_ESP_L3x5y_F3y_aa+ABX*I_ESP_K2x5y_F3y_aa;
    Double I_ESP_K2x4yz_Gx3y_aa = I_ESP_L3x4yz_F3y_aa+ABX*I_ESP_K2x4yz_F3y_aa;
    Double I_ESP_K2x3y2z_Gx3y_aa = I_ESP_L3x3y2z_F3y_aa+ABX*I_ESP_K2x3y2z_F3y_aa;
    Double I_ESP_K2x2y3z_Gx3y_aa = I_ESP_L3x2y3z_F3y_aa+ABX*I_ESP_K2x2y3z_F3y_aa;
    Double I_ESP_K2xy4z_Gx3y_aa = I_ESP_L3xy4z_F3y_aa+ABX*I_ESP_K2xy4z_F3y_aa;
    Double I_ESP_K2x5z_Gx3y_aa = I_ESP_L3x5z_F3y_aa+ABX*I_ESP_K2x5z_F3y_aa;
    Double I_ESP_Kx6y_Gx3y_aa = I_ESP_L2x6y_F3y_aa+ABX*I_ESP_Kx6y_F3y_aa;
    Double I_ESP_Kx5yz_Gx3y_aa = I_ESP_L2x5yz_F3y_aa+ABX*I_ESP_Kx5yz_F3y_aa;
    Double I_ESP_Kx4y2z_Gx3y_aa = I_ESP_L2x4y2z_F3y_aa+ABX*I_ESP_Kx4y2z_F3y_aa;
    Double I_ESP_Kx3y3z_Gx3y_aa = I_ESP_L2x3y3z_F3y_aa+ABX*I_ESP_Kx3y3z_F3y_aa;
    Double I_ESP_Kx2y4z_Gx3y_aa = I_ESP_L2x2y4z_F3y_aa+ABX*I_ESP_Kx2y4z_F3y_aa;
    Double I_ESP_Kxy5z_Gx3y_aa = I_ESP_L2xy5z_F3y_aa+ABX*I_ESP_Kxy5z_F3y_aa;
    Double I_ESP_Kx6z_Gx3y_aa = I_ESP_L2x6z_F3y_aa+ABX*I_ESP_Kx6z_F3y_aa;
    Double I_ESP_K7y_Gx3y_aa = I_ESP_Lx7y_F3y_aa+ABX*I_ESP_K7y_F3y_aa;
    Double I_ESP_K6yz_Gx3y_aa = I_ESP_Lx6yz_F3y_aa+ABX*I_ESP_K6yz_F3y_aa;
    Double I_ESP_K5y2z_Gx3y_aa = I_ESP_Lx5y2z_F3y_aa+ABX*I_ESP_K5y2z_F3y_aa;
    Double I_ESP_K4y3z_Gx3y_aa = I_ESP_Lx4y3z_F3y_aa+ABX*I_ESP_K4y3z_F3y_aa;
    Double I_ESP_K3y4z_Gx3y_aa = I_ESP_Lx3y4z_F3y_aa+ABX*I_ESP_K3y4z_F3y_aa;
    Double I_ESP_K2y5z_Gx3y_aa = I_ESP_Lx2y5z_F3y_aa+ABX*I_ESP_K2y5z_F3y_aa;
    Double I_ESP_Ky6z_Gx3y_aa = I_ESP_Lxy6z_F3y_aa+ABX*I_ESP_Ky6z_F3y_aa;
    Double I_ESP_K7z_Gx3y_aa = I_ESP_Lx7z_F3y_aa+ABX*I_ESP_K7z_F3y_aa;
    Double I_ESP_K7x_Gx2yz_aa = I_ESP_L7xz_Fx2y_aa+ABZ*I_ESP_K7x_Fx2y_aa;
    Double I_ESP_K6xy_Gx2yz_aa = I_ESP_L6xyz_Fx2y_aa+ABZ*I_ESP_K6xy_Fx2y_aa;
    Double I_ESP_K6xz_Gx2yz_aa = I_ESP_L6x2z_Fx2y_aa+ABZ*I_ESP_K6xz_Fx2y_aa;
    Double I_ESP_K5x2y_Gx2yz_aa = I_ESP_L5x2yz_Fx2y_aa+ABZ*I_ESP_K5x2y_Fx2y_aa;
    Double I_ESP_K5xyz_Gx2yz_aa = I_ESP_L5xy2z_Fx2y_aa+ABZ*I_ESP_K5xyz_Fx2y_aa;
    Double I_ESP_K5x2z_Gx2yz_aa = I_ESP_L5x3z_Fx2y_aa+ABZ*I_ESP_K5x2z_Fx2y_aa;
    Double I_ESP_K4x3y_Gx2yz_aa = I_ESP_L4x3yz_Fx2y_aa+ABZ*I_ESP_K4x3y_Fx2y_aa;
    Double I_ESP_K4x2yz_Gx2yz_aa = I_ESP_L4x2y2z_Fx2y_aa+ABZ*I_ESP_K4x2yz_Fx2y_aa;
    Double I_ESP_K4xy2z_Gx2yz_aa = I_ESP_L4xy3z_Fx2y_aa+ABZ*I_ESP_K4xy2z_Fx2y_aa;
    Double I_ESP_K4x3z_Gx2yz_aa = I_ESP_L4x4z_Fx2y_aa+ABZ*I_ESP_K4x3z_Fx2y_aa;
    Double I_ESP_K3x4y_Gx2yz_aa = I_ESP_L3x4yz_Fx2y_aa+ABZ*I_ESP_K3x4y_Fx2y_aa;
    Double I_ESP_K3x3yz_Gx2yz_aa = I_ESP_L3x3y2z_Fx2y_aa+ABZ*I_ESP_K3x3yz_Fx2y_aa;
    Double I_ESP_K3x2y2z_Gx2yz_aa = I_ESP_L3x2y3z_Fx2y_aa+ABZ*I_ESP_K3x2y2z_Fx2y_aa;
    Double I_ESP_K3xy3z_Gx2yz_aa = I_ESP_L3xy4z_Fx2y_aa+ABZ*I_ESP_K3xy3z_Fx2y_aa;
    Double I_ESP_K3x4z_Gx2yz_aa = I_ESP_L3x5z_Fx2y_aa+ABZ*I_ESP_K3x4z_Fx2y_aa;
    Double I_ESP_K2x5y_Gx2yz_aa = I_ESP_L2x5yz_Fx2y_aa+ABZ*I_ESP_K2x5y_Fx2y_aa;
    Double I_ESP_K2x4yz_Gx2yz_aa = I_ESP_L2x4y2z_Fx2y_aa+ABZ*I_ESP_K2x4yz_Fx2y_aa;
    Double I_ESP_K2x3y2z_Gx2yz_aa = I_ESP_L2x3y3z_Fx2y_aa+ABZ*I_ESP_K2x3y2z_Fx2y_aa;
    Double I_ESP_K2x2y3z_Gx2yz_aa = I_ESP_L2x2y4z_Fx2y_aa+ABZ*I_ESP_K2x2y3z_Fx2y_aa;
    Double I_ESP_K2xy4z_Gx2yz_aa = I_ESP_L2xy5z_Fx2y_aa+ABZ*I_ESP_K2xy4z_Fx2y_aa;
    Double I_ESP_K2x5z_Gx2yz_aa = I_ESP_L2x6z_Fx2y_aa+ABZ*I_ESP_K2x5z_Fx2y_aa;
    Double I_ESP_Kx6y_Gx2yz_aa = I_ESP_Lx6yz_Fx2y_aa+ABZ*I_ESP_Kx6y_Fx2y_aa;
    Double I_ESP_Kx5yz_Gx2yz_aa = I_ESP_Lx5y2z_Fx2y_aa+ABZ*I_ESP_Kx5yz_Fx2y_aa;
    Double I_ESP_Kx4y2z_Gx2yz_aa = I_ESP_Lx4y3z_Fx2y_aa+ABZ*I_ESP_Kx4y2z_Fx2y_aa;
    Double I_ESP_Kx3y3z_Gx2yz_aa = I_ESP_Lx3y4z_Fx2y_aa+ABZ*I_ESP_Kx3y3z_Fx2y_aa;
    Double I_ESP_Kx2y4z_Gx2yz_aa = I_ESP_Lx2y5z_Fx2y_aa+ABZ*I_ESP_Kx2y4z_Fx2y_aa;
    Double I_ESP_Kxy5z_Gx2yz_aa = I_ESP_Lxy6z_Fx2y_aa+ABZ*I_ESP_Kxy5z_Fx2y_aa;
    Double I_ESP_Kx6z_Gx2yz_aa = I_ESP_Lx7z_Fx2y_aa+ABZ*I_ESP_Kx6z_Fx2y_aa;
    Double I_ESP_K7y_Gx2yz_aa = I_ESP_L7yz_Fx2y_aa+ABZ*I_ESP_K7y_Fx2y_aa;
    Double I_ESP_K6yz_Gx2yz_aa = I_ESP_L6y2z_Fx2y_aa+ABZ*I_ESP_K6yz_Fx2y_aa;
    Double I_ESP_K5y2z_Gx2yz_aa = I_ESP_L5y3z_Fx2y_aa+ABZ*I_ESP_K5y2z_Fx2y_aa;
    Double I_ESP_K4y3z_Gx2yz_aa = I_ESP_L4y4z_Fx2y_aa+ABZ*I_ESP_K4y3z_Fx2y_aa;
    Double I_ESP_K3y4z_Gx2yz_aa = I_ESP_L3y5z_Fx2y_aa+ABZ*I_ESP_K3y4z_Fx2y_aa;
    Double I_ESP_K2y5z_Gx2yz_aa = I_ESP_L2y6z_Fx2y_aa+ABZ*I_ESP_K2y5z_Fx2y_aa;
    Double I_ESP_Ky6z_Gx2yz_aa = I_ESP_Ly7z_Fx2y_aa+ABZ*I_ESP_Ky6z_Fx2y_aa;
    Double I_ESP_K7z_Gx2yz_aa = I_ESP_L8z_Fx2y_aa+ABZ*I_ESP_K7z_Fx2y_aa;
    Double I_ESP_K7x_Gxy2z_aa = I_ESP_L7xy_Fx2z_aa+ABY*I_ESP_K7x_Fx2z_aa;
    Double I_ESP_K6xy_Gxy2z_aa = I_ESP_L6x2y_Fx2z_aa+ABY*I_ESP_K6xy_Fx2z_aa;
    Double I_ESP_K6xz_Gxy2z_aa = I_ESP_L6xyz_Fx2z_aa+ABY*I_ESP_K6xz_Fx2z_aa;
    Double I_ESP_K5x2y_Gxy2z_aa = I_ESP_L5x3y_Fx2z_aa+ABY*I_ESP_K5x2y_Fx2z_aa;
    Double I_ESP_K5xyz_Gxy2z_aa = I_ESP_L5x2yz_Fx2z_aa+ABY*I_ESP_K5xyz_Fx2z_aa;
    Double I_ESP_K5x2z_Gxy2z_aa = I_ESP_L5xy2z_Fx2z_aa+ABY*I_ESP_K5x2z_Fx2z_aa;
    Double I_ESP_K4x3y_Gxy2z_aa = I_ESP_L4x4y_Fx2z_aa+ABY*I_ESP_K4x3y_Fx2z_aa;
    Double I_ESP_K4x2yz_Gxy2z_aa = I_ESP_L4x3yz_Fx2z_aa+ABY*I_ESP_K4x2yz_Fx2z_aa;
    Double I_ESP_K4xy2z_Gxy2z_aa = I_ESP_L4x2y2z_Fx2z_aa+ABY*I_ESP_K4xy2z_Fx2z_aa;
    Double I_ESP_K4x3z_Gxy2z_aa = I_ESP_L4xy3z_Fx2z_aa+ABY*I_ESP_K4x3z_Fx2z_aa;
    Double I_ESP_K3x4y_Gxy2z_aa = I_ESP_L3x5y_Fx2z_aa+ABY*I_ESP_K3x4y_Fx2z_aa;
    Double I_ESP_K3x3yz_Gxy2z_aa = I_ESP_L3x4yz_Fx2z_aa+ABY*I_ESP_K3x3yz_Fx2z_aa;
    Double I_ESP_K3x2y2z_Gxy2z_aa = I_ESP_L3x3y2z_Fx2z_aa+ABY*I_ESP_K3x2y2z_Fx2z_aa;
    Double I_ESP_K3xy3z_Gxy2z_aa = I_ESP_L3x2y3z_Fx2z_aa+ABY*I_ESP_K3xy3z_Fx2z_aa;
    Double I_ESP_K3x4z_Gxy2z_aa = I_ESP_L3xy4z_Fx2z_aa+ABY*I_ESP_K3x4z_Fx2z_aa;
    Double I_ESP_K2x5y_Gxy2z_aa = I_ESP_L2x6y_Fx2z_aa+ABY*I_ESP_K2x5y_Fx2z_aa;
    Double I_ESP_K2x4yz_Gxy2z_aa = I_ESP_L2x5yz_Fx2z_aa+ABY*I_ESP_K2x4yz_Fx2z_aa;
    Double I_ESP_K2x3y2z_Gxy2z_aa = I_ESP_L2x4y2z_Fx2z_aa+ABY*I_ESP_K2x3y2z_Fx2z_aa;
    Double I_ESP_K2x2y3z_Gxy2z_aa = I_ESP_L2x3y3z_Fx2z_aa+ABY*I_ESP_K2x2y3z_Fx2z_aa;
    Double I_ESP_K2xy4z_Gxy2z_aa = I_ESP_L2x2y4z_Fx2z_aa+ABY*I_ESP_K2xy4z_Fx2z_aa;
    Double I_ESP_K2x5z_Gxy2z_aa = I_ESP_L2xy5z_Fx2z_aa+ABY*I_ESP_K2x5z_Fx2z_aa;
    Double I_ESP_Kx6y_Gxy2z_aa = I_ESP_Lx7y_Fx2z_aa+ABY*I_ESP_Kx6y_Fx2z_aa;
    Double I_ESP_Kx5yz_Gxy2z_aa = I_ESP_Lx6yz_Fx2z_aa+ABY*I_ESP_Kx5yz_Fx2z_aa;
    Double I_ESP_Kx4y2z_Gxy2z_aa = I_ESP_Lx5y2z_Fx2z_aa+ABY*I_ESP_Kx4y2z_Fx2z_aa;
    Double I_ESP_Kx3y3z_Gxy2z_aa = I_ESP_Lx4y3z_Fx2z_aa+ABY*I_ESP_Kx3y3z_Fx2z_aa;
    Double I_ESP_Kx2y4z_Gxy2z_aa = I_ESP_Lx3y4z_Fx2z_aa+ABY*I_ESP_Kx2y4z_Fx2z_aa;
    Double I_ESP_Kxy5z_Gxy2z_aa = I_ESP_Lx2y5z_Fx2z_aa+ABY*I_ESP_Kxy5z_Fx2z_aa;
    Double I_ESP_Kx6z_Gxy2z_aa = I_ESP_Lxy6z_Fx2z_aa+ABY*I_ESP_Kx6z_Fx2z_aa;
    Double I_ESP_K7y_Gxy2z_aa = I_ESP_L8y_Fx2z_aa+ABY*I_ESP_K7y_Fx2z_aa;
    Double I_ESP_K6yz_Gxy2z_aa = I_ESP_L7yz_Fx2z_aa+ABY*I_ESP_K6yz_Fx2z_aa;
    Double I_ESP_K5y2z_Gxy2z_aa = I_ESP_L6y2z_Fx2z_aa+ABY*I_ESP_K5y2z_Fx2z_aa;
    Double I_ESP_K4y3z_Gxy2z_aa = I_ESP_L5y3z_Fx2z_aa+ABY*I_ESP_K4y3z_Fx2z_aa;
    Double I_ESP_K3y4z_Gxy2z_aa = I_ESP_L4y4z_Fx2z_aa+ABY*I_ESP_K3y4z_Fx2z_aa;
    Double I_ESP_K2y5z_Gxy2z_aa = I_ESP_L3y5z_Fx2z_aa+ABY*I_ESP_K2y5z_Fx2z_aa;
    Double I_ESP_Ky6z_Gxy2z_aa = I_ESP_L2y6z_Fx2z_aa+ABY*I_ESP_Ky6z_Fx2z_aa;
    Double I_ESP_K7z_Gxy2z_aa = I_ESP_Ly7z_Fx2z_aa+ABY*I_ESP_K7z_Fx2z_aa;
    Double I_ESP_K7x_Gx3z_aa = I_ESP_L8x_F3z_aa+ABX*I_ESP_K7x_F3z_aa;
    Double I_ESP_K6xy_Gx3z_aa = I_ESP_L7xy_F3z_aa+ABX*I_ESP_K6xy_F3z_aa;
    Double I_ESP_K6xz_Gx3z_aa = I_ESP_L7xz_F3z_aa+ABX*I_ESP_K6xz_F3z_aa;
    Double I_ESP_K5x2y_Gx3z_aa = I_ESP_L6x2y_F3z_aa+ABX*I_ESP_K5x2y_F3z_aa;
    Double I_ESP_K5xyz_Gx3z_aa = I_ESP_L6xyz_F3z_aa+ABX*I_ESP_K5xyz_F3z_aa;
    Double I_ESP_K5x2z_Gx3z_aa = I_ESP_L6x2z_F3z_aa+ABX*I_ESP_K5x2z_F3z_aa;
    Double I_ESP_K4x3y_Gx3z_aa = I_ESP_L5x3y_F3z_aa+ABX*I_ESP_K4x3y_F3z_aa;
    Double I_ESP_K4x2yz_Gx3z_aa = I_ESP_L5x2yz_F3z_aa+ABX*I_ESP_K4x2yz_F3z_aa;
    Double I_ESP_K4xy2z_Gx3z_aa = I_ESP_L5xy2z_F3z_aa+ABX*I_ESP_K4xy2z_F3z_aa;
    Double I_ESP_K4x3z_Gx3z_aa = I_ESP_L5x3z_F3z_aa+ABX*I_ESP_K4x3z_F3z_aa;
    Double I_ESP_K3x4y_Gx3z_aa = I_ESP_L4x4y_F3z_aa+ABX*I_ESP_K3x4y_F3z_aa;
    Double I_ESP_K3x3yz_Gx3z_aa = I_ESP_L4x3yz_F3z_aa+ABX*I_ESP_K3x3yz_F3z_aa;
    Double I_ESP_K3x2y2z_Gx3z_aa = I_ESP_L4x2y2z_F3z_aa+ABX*I_ESP_K3x2y2z_F3z_aa;
    Double I_ESP_K3xy3z_Gx3z_aa = I_ESP_L4xy3z_F3z_aa+ABX*I_ESP_K3xy3z_F3z_aa;
    Double I_ESP_K3x4z_Gx3z_aa = I_ESP_L4x4z_F3z_aa+ABX*I_ESP_K3x4z_F3z_aa;
    Double I_ESP_K2x5y_Gx3z_aa = I_ESP_L3x5y_F3z_aa+ABX*I_ESP_K2x5y_F3z_aa;
    Double I_ESP_K2x4yz_Gx3z_aa = I_ESP_L3x4yz_F3z_aa+ABX*I_ESP_K2x4yz_F3z_aa;
    Double I_ESP_K2x3y2z_Gx3z_aa = I_ESP_L3x3y2z_F3z_aa+ABX*I_ESP_K2x3y2z_F3z_aa;
    Double I_ESP_K2x2y3z_Gx3z_aa = I_ESP_L3x2y3z_F3z_aa+ABX*I_ESP_K2x2y3z_F3z_aa;
    Double I_ESP_K2xy4z_Gx3z_aa = I_ESP_L3xy4z_F3z_aa+ABX*I_ESP_K2xy4z_F3z_aa;
    Double I_ESP_K2x5z_Gx3z_aa = I_ESP_L3x5z_F3z_aa+ABX*I_ESP_K2x5z_F3z_aa;
    Double I_ESP_Kx6y_Gx3z_aa = I_ESP_L2x6y_F3z_aa+ABX*I_ESP_Kx6y_F3z_aa;
    Double I_ESP_Kx5yz_Gx3z_aa = I_ESP_L2x5yz_F3z_aa+ABX*I_ESP_Kx5yz_F3z_aa;
    Double I_ESP_Kx4y2z_Gx3z_aa = I_ESP_L2x4y2z_F3z_aa+ABX*I_ESP_Kx4y2z_F3z_aa;
    Double I_ESP_Kx3y3z_Gx3z_aa = I_ESP_L2x3y3z_F3z_aa+ABX*I_ESP_Kx3y3z_F3z_aa;
    Double I_ESP_Kx2y4z_Gx3z_aa = I_ESP_L2x2y4z_F3z_aa+ABX*I_ESP_Kx2y4z_F3z_aa;
    Double I_ESP_Kxy5z_Gx3z_aa = I_ESP_L2xy5z_F3z_aa+ABX*I_ESP_Kxy5z_F3z_aa;
    Double I_ESP_Kx6z_Gx3z_aa = I_ESP_L2x6z_F3z_aa+ABX*I_ESP_Kx6z_F3z_aa;
    Double I_ESP_K7y_Gx3z_aa = I_ESP_Lx7y_F3z_aa+ABX*I_ESP_K7y_F3z_aa;
    Double I_ESP_K6yz_Gx3z_aa = I_ESP_Lx6yz_F3z_aa+ABX*I_ESP_K6yz_F3z_aa;
    Double I_ESP_K5y2z_Gx3z_aa = I_ESP_Lx5y2z_F3z_aa+ABX*I_ESP_K5y2z_F3z_aa;
    Double I_ESP_K4y3z_Gx3z_aa = I_ESP_Lx4y3z_F3z_aa+ABX*I_ESP_K4y3z_F3z_aa;
    Double I_ESP_K3y4z_Gx3z_aa = I_ESP_Lx3y4z_F3z_aa+ABX*I_ESP_K3y4z_F3z_aa;
    Double I_ESP_K2y5z_Gx3z_aa = I_ESP_Lx2y5z_F3z_aa+ABX*I_ESP_K2y5z_F3z_aa;
    Double I_ESP_Ky6z_Gx3z_aa = I_ESP_Lxy6z_F3z_aa+ABX*I_ESP_Ky6z_F3z_aa;
    Double I_ESP_K7z_Gx3z_aa = I_ESP_Lx7z_F3z_aa+ABX*I_ESP_K7z_F3z_aa;
    Double I_ESP_K7x_G4y_aa = I_ESP_L7xy_F3y_aa+ABY*I_ESP_K7x_F3y_aa;
    Double I_ESP_K6xy_G4y_aa = I_ESP_L6x2y_F3y_aa+ABY*I_ESP_K6xy_F3y_aa;
    Double I_ESP_K6xz_G4y_aa = I_ESP_L6xyz_F3y_aa+ABY*I_ESP_K6xz_F3y_aa;
    Double I_ESP_K5x2y_G4y_aa = I_ESP_L5x3y_F3y_aa+ABY*I_ESP_K5x2y_F3y_aa;
    Double I_ESP_K5xyz_G4y_aa = I_ESP_L5x2yz_F3y_aa+ABY*I_ESP_K5xyz_F3y_aa;
    Double I_ESP_K5x2z_G4y_aa = I_ESP_L5xy2z_F3y_aa+ABY*I_ESP_K5x2z_F3y_aa;
    Double I_ESP_K4x3y_G4y_aa = I_ESP_L4x4y_F3y_aa+ABY*I_ESP_K4x3y_F3y_aa;
    Double I_ESP_K4x2yz_G4y_aa = I_ESP_L4x3yz_F3y_aa+ABY*I_ESP_K4x2yz_F3y_aa;
    Double I_ESP_K4xy2z_G4y_aa = I_ESP_L4x2y2z_F3y_aa+ABY*I_ESP_K4xy2z_F3y_aa;
    Double I_ESP_K4x3z_G4y_aa = I_ESP_L4xy3z_F3y_aa+ABY*I_ESP_K4x3z_F3y_aa;
    Double I_ESP_K3x4y_G4y_aa = I_ESP_L3x5y_F3y_aa+ABY*I_ESP_K3x4y_F3y_aa;
    Double I_ESP_K3x3yz_G4y_aa = I_ESP_L3x4yz_F3y_aa+ABY*I_ESP_K3x3yz_F3y_aa;
    Double I_ESP_K3x2y2z_G4y_aa = I_ESP_L3x3y2z_F3y_aa+ABY*I_ESP_K3x2y2z_F3y_aa;
    Double I_ESP_K3xy3z_G4y_aa = I_ESP_L3x2y3z_F3y_aa+ABY*I_ESP_K3xy3z_F3y_aa;
    Double I_ESP_K3x4z_G4y_aa = I_ESP_L3xy4z_F3y_aa+ABY*I_ESP_K3x4z_F3y_aa;
    Double I_ESP_K2x5y_G4y_aa = I_ESP_L2x6y_F3y_aa+ABY*I_ESP_K2x5y_F3y_aa;
    Double I_ESP_K2x4yz_G4y_aa = I_ESP_L2x5yz_F3y_aa+ABY*I_ESP_K2x4yz_F3y_aa;
    Double I_ESP_K2x3y2z_G4y_aa = I_ESP_L2x4y2z_F3y_aa+ABY*I_ESP_K2x3y2z_F3y_aa;
    Double I_ESP_K2x2y3z_G4y_aa = I_ESP_L2x3y3z_F3y_aa+ABY*I_ESP_K2x2y3z_F3y_aa;
    Double I_ESP_K2xy4z_G4y_aa = I_ESP_L2x2y4z_F3y_aa+ABY*I_ESP_K2xy4z_F3y_aa;
    Double I_ESP_K2x5z_G4y_aa = I_ESP_L2xy5z_F3y_aa+ABY*I_ESP_K2x5z_F3y_aa;
    Double I_ESP_Kx6y_G4y_aa = I_ESP_Lx7y_F3y_aa+ABY*I_ESP_Kx6y_F3y_aa;
    Double I_ESP_Kx5yz_G4y_aa = I_ESP_Lx6yz_F3y_aa+ABY*I_ESP_Kx5yz_F3y_aa;
    Double I_ESP_Kx4y2z_G4y_aa = I_ESP_Lx5y2z_F3y_aa+ABY*I_ESP_Kx4y2z_F3y_aa;
    Double I_ESP_Kx3y3z_G4y_aa = I_ESP_Lx4y3z_F3y_aa+ABY*I_ESP_Kx3y3z_F3y_aa;
    Double I_ESP_Kx2y4z_G4y_aa = I_ESP_Lx3y4z_F3y_aa+ABY*I_ESP_Kx2y4z_F3y_aa;
    Double I_ESP_Kxy5z_G4y_aa = I_ESP_Lx2y5z_F3y_aa+ABY*I_ESP_Kxy5z_F3y_aa;
    Double I_ESP_Kx6z_G4y_aa = I_ESP_Lxy6z_F3y_aa+ABY*I_ESP_Kx6z_F3y_aa;
    Double I_ESP_K7y_G4y_aa = I_ESP_L8y_F3y_aa+ABY*I_ESP_K7y_F3y_aa;
    Double I_ESP_K6yz_G4y_aa = I_ESP_L7yz_F3y_aa+ABY*I_ESP_K6yz_F3y_aa;
    Double I_ESP_K5y2z_G4y_aa = I_ESP_L6y2z_F3y_aa+ABY*I_ESP_K5y2z_F3y_aa;
    Double I_ESP_K4y3z_G4y_aa = I_ESP_L5y3z_F3y_aa+ABY*I_ESP_K4y3z_F3y_aa;
    Double I_ESP_K3y4z_G4y_aa = I_ESP_L4y4z_F3y_aa+ABY*I_ESP_K3y4z_F3y_aa;
    Double I_ESP_K2y5z_G4y_aa = I_ESP_L3y5z_F3y_aa+ABY*I_ESP_K2y5z_F3y_aa;
    Double I_ESP_Ky6z_G4y_aa = I_ESP_L2y6z_F3y_aa+ABY*I_ESP_Ky6z_F3y_aa;
    Double I_ESP_K7z_G4y_aa = I_ESP_Ly7z_F3y_aa+ABY*I_ESP_K7z_F3y_aa;
    Double I_ESP_K7x_G3yz_aa = I_ESP_L7xz_F3y_aa+ABZ*I_ESP_K7x_F3y_aa;
    Double I_ESP_K6xy_G3yz_aa = I_ESP_L6xyz_F3y_aa+ABZ*I_ESP_K6xy_F3y_aa;
    Double I_ESP_K6xz_G3yz_aa = I_ESP_L6x2z_F3y_aa+ABZ*I_ESP_K6xz_F3y_aa;
    Double I_ESP_K5x2y_G3yz_aa = I_ESP_L5x2yz_F3y_aa+ABZ*I_ESP_K5x2y_F3y_aa;
    Double I_ESP_K5xyz_G3yz_aa = I_ESP_L5xy2z_F3y_aa+ABZ*I_ESP_K5xyz_F3y_aa;
    Double I_ESP_K5x2z_G3yz_aa = I_ESP_L5x3z_F3y_aa+ABZ*I_ESP_K5x2z_F3y_aa;
    Double I_ESP_K4x3y_G3yz_aa = I_ESP_L4x3yz_F3y_aa+ABZ*I_ESP_K4x3y_F3y_aa;
    Double I_ESP_K4x2yz_G3yz_aa = I_ESP_L4x2y2z_F3y_aa+ABZ*I_ESP_K4x2yz_F3y_aa;
    Double I_ESP_K4xy2z_G3yz_aa = I_ESP_L4xy3z_F3y_aa+ABZ*I_ESP_K4xy2z_F3y_aa;
    Double I_ESP_K4x3z_G3yz_aa = I_ESP_L4x4z_F3y_aa+ABZ*I_ESP_K4x3z_F3y_aa;
    Double I_ESP_K3x4y_G3yz_aa = I_ESP_L3x4yz_F3y_aa+ABZ*I_ESP_K3x4y_F3y_aa;
    Double I_ESP_K3x3yz_G3yz_aa = I_ESP_L3x3y2z_F3y_aa+ABZ*I_ESP_K3x3yz_F3y_aa;
    Double I_ESP_K3x2y2z_G3yz_aa = I_ESP_L3x2y3z_F3y_aa+ABZ*I_ESP_K3x2y2z_F3y_aa;
    Double I_ESP_K3xy3z_G3yz_aa = I_ESP_L3xy4z_F3y_aa+ABZ*I_ESP_K3xy3z_F3y_aa;
    Double I_ESP_K3x4z_G3yz_aa = I_ESP_L3x5z_F3y_aa+ABZ*I_ESP_K3x4z_F3y_aa;
    Double I_ESP_K2x5y_G3yz_aa = I_ESP_L2x5yz_F3y_aa+ABZ*I_ESP_K2x5y_F3y_aa;
    Double I_ESP_K2x4yz_G3yz_aa = I_ESP_L2x4y2z_F3y_aa+ABZ*I_ESP_K2x4yz_F3y_aa;
    Double I_ESP_K2x3y2z_G3yz_aa = I_ESP_L2x3y3z_F3y_aa+ABZ*I_ESP_K2x3y2z_F3y_aa;
    Double I_ESP_K2x2y3z_G3yz_aa = I_ESP_L2x2y4z_F3y_aa+ABZ*I_ESP_K2x2y3z_F3y_aa;
    Double I_ESP_K2xy4z_G3yz_aa = I_ESP_L2xy5z_F3y_aa+ABZ*I_ESP_K2xy4z_F3y_aa;
    Double I_ESP_K2x5z_G3yz_aa = I_ESP_L2x6z_F3y_aa+ABZ*I_ESP_K2x5z_F3y_aa;
    Double I_ESP_Kx6y_G3yz_aa = I_ESP_Lx6yz_F3y_aa+ABZ*I_ESP_Kx6y_F3y_aa;
    Double I_ESP_Kx5yz_G3yz_aa = I_ESP_Lx5y2z_F3y_aa+ABZ*I_ESP_Kx5yz_F3y_aa;
    Double I_ESP_Kx4y2z_G3yz_aa = I_ESP_Lx4y3z_F3y_aa+ABZ*I_ESP_Kx4y2z_F3y_aa;
    Double I_ESP_Kx3y3z_G3yz_aa = I_ESP_Lx3y4z_F3y_aa+ABZ*I_ESP_Kx3y3z_F3y_aa;
    Double I_ESP_Kx2y4z_G3yz_aa = I_ESP_Lx2y5z_F3y_aa+ABZ*I_ESP_Kx2y4z_F3y_aa;
    Double I_ESP_Kxy5z_G3yz_aa = I_ESP_Lxy6z_F3y_aa+ABZ*I_ESP_Kxy5z_F3y_aa;
    Double I_ESP_Kx6z_G3yz_aa = I_ESP_Lx7z_F3y_aa+ABZ*I_ESP_Kx6z_F3y_aa;
    Double I_ESP_K7y_G3yz_aa = I_ESP_L7yz_F3y_aa+ABZ*I_ESP_K7y_F3y_aa;
    Double I_ESP_K6yz_G3yz_aa = I_ESP_L6y2z_F3y_aa+ABZ*I_ESP_K6yz_F3y_aa;
    Double I_ESP_K5y2z_G3yz_aa = I_ESP_L5y3z_F3y_aa+ABZ*I_ESP_K5y2z_F3y_aa;
    Double I_ESP_K4y3z_G3yz_aa = I_ESP_L4y4z_F3y_aa+ABZ*I_ESP_K4y3z_F3y_aa;
    Double I_ESP_K3y4z_G3yz_aa = I_ESP_L3y5z_F3y_aa+ABZ*I_ESP_K3y4z_F3y_aa;
    Double I_ESP_K2y5z_G3yz_aa = I_ESP_L2y6z_F3y_aa+ABZ*I_ESP_K2y5z_F3y_aa;
    Double I_ESP_Ky6z_G3yz_aa = I_ESP_Ly7z_F3y_aa+ABZ*I_ESP_Ky6z_F3y_aa;
    Double I_ESP_K7z_G3yz_aa = I_ESP_L8z_F3y_aa+ABZ*I_ESP_K7z_F3y_aa;
    Double I_ESP_K7x_G2y2z_aa = I_ESP_L7xz_F2yz_aa+ABZ*I_ESP_K7x_F2yz_aa;
    Double I_ESP_K6xy_G2y2z_aa = I_ESP_L6xyz_F2yz_aa+ABZ*I_ESP_K6xy_F2yz_aa;
    Double I_ESP_K6xz_G2y2z_aa = I_ESP_L6x2z_F2yz_aa+ABZ*I_ESP_K6xz_F2yz_aa;
    Double I_ESP_K5x2y_G2y2z_aa = I_ESP_L5x2yz_F2yz_aa+ABZ*I_ESP_K5x2y_F2yz_aa;
    Double I_ESP_K5xyz_G2y2z_aa = I_ESP_L5xy2z_F2yz_aa+ABZ*I_ESP_K5xyz_F2yz_aa;
    Double I_ESP_K5x2z_G2y2z_aa = I_ESP_L5x3z_F2yz_aa+ABZ*I_ESP_K5x2z_F2yz_aa;
    Double I_ESP_K4x3y_G2y2z_aa = I_ESP_L4x3yz_F2yz_aa+ABZ*I_ESP_K4x3y_F2yz_aa;
    Double I_ESP_K4x2yz_G2y2z_aa = I_ESP_L4x2y2z_F2yz_aa+ABZ*I_ESP_K4x2yz_F2yz_aa;
    Double I_ESP_K4xy2z_G2y2z_aa = I_ESP_L4xy3z_F2yz_aa+ABZ*I_ESP_K4xy2z_F2yz_aa;
    Double I_ESP_K4x3z_G2y2z_aa = I_ESP_L4x4z_F2yz_aa+ABZ*I_ESP_K4x3z_F2yz_aa;
    Double I_ESP_K3x4y_G2y2z_aa = I_ESP_L3x4yz_F2yz_aa+ABZ*I_ESP_K3x4y_F2yz_aa;
    Double I_ESP_K3x3yz_G2y2z_aa = I_ESP_L3x3y2z_F2yz_aa+ABZ*I_ESP_K3x3yz_F2yz_aa;
    Double I_ESP_K3x2y2z_G2y2z_aa = I_ESP_L3x2y3z_F2yz_aa+ABZ*I_ESP_K3x2y2z_F2yz_aa;
    Double I_ESP_K3xy3z_G2y2z_aa = I_ESP_L3xy4z_F2yz_aa+ABZ*I_ESP_K3xy3z_F2yz_aa;
    Double I_ESP_K3x4z_G2y2z_aa = I_ESP_L3x5z_F2yz_aa+ABZ*I_ESP_K3x4z_F2yz_aa;
    Double I_ESP_K2x5y_G2y2z_aa = I_ESP_L2x5yz_F2yz_aa+ABZ*I_ESP_K2x5y_F2yz_aa;
    Double I_ESP_K2x4yz_G2y2z_aa = I_ESP_L2x4y2z_F2yz_aa+ABZ*I_ESP_K2x4yz_F2yz_aa;
    Double I_ESP_K2x3y2z_G2y2z_aa = I_ESP_L2x3y3z_F2yz_aa+ABZ*I_ESP_K2x3y2z_F2yz_aa;
    Double I_ESP_K2x2y3z_G2y2z_aa = I_ESP_L2x2y4z_F2yz_aa+ABZ*I_ESP_K2x2y3z_F2yz_aa;
    Double I_ESP_K2xy4z_G2y2z_aa = I_ESP_L2xy5z_F2yz_aa+ABZ*I_ESP_K2xy4z_F2yz_aa;
    Double I_ESP_K2x5z_G2y2z_aa = I_ESP_L2x6z_F2yz_aa+ABZ*I_ESP_K2x5z_F2yz_aa;
    Double I_ESP_Kx6y_G2y2z_aa = I_ESP_Lx6yz_F2yz_aa+ABZ*I_ESP_Kx6y_F2yz_aa;
    Double I_ESP_Kx5yz_G2y2z_aa = I_ESP_Lx5y2z_F2yz_aa+ABZ*I_ESP_Kx5yz_F2yz_aa;
    Double I_ESP_Kx4y2z_G2y2z_aa = I_ESP_Lx4y3z_F2yz_aa+ABZ*I_ESP_Kx4y2z_F2yz_aa;
    Double I_ESP_Kx3y3z_G2y2z_aa = I_ESP_Lx3y4z_F2yz_aa+ABZ*I_ESP_Kx3y3z_F2yz_aa;
    Double I_ESP_Kx2y4z_G2y2z_aa = I_ESP_Lx2y5z_F2yz_aa+ABZ*I_ESP_Kx2y4z_F2yz_aa;
    Double I_ESP_Kxy5z_G2y2z_aa = I_ESP_Lxy6z_F2yz_aa+ABZ*I_ESP_Kxy5z_F2yz_aa;
    Double I_ESP_Kx6z_G2y2z_aa = I_ESP_Lx7z_F2yz_aa+ABZ*I_ESP_Kx6z_F2yz_aa;
    Double I_ESP_K7y_G2y2z_aa = I_ESP_L7yz_F2yz_aa+ABZ*I_ESP_K7y_F2yz_aa;
    Double I_ESP_K6yz_G2y2z_aa = I_ESP_L6y2z_F2yz_aa+ABZ*I_ESP_K6yz_F2yz_aa;
    Double I_ESP_K5y2z_G2y2z_aa = I_ESP_L5y3z_F2yz_aa+ABZ*I_ESP_K5y2z_F2yz_aa;
    Double I_ESP_K4y3z_G2y2z_aa = I_ESP_L4y4z_F2yz_aa+ABZ*I_ESP_K4y3z_F2yz_aa;
    Double I_ESP_K3y4z_G2y2z_aa = I_ESP_L3y5z_F2yz_aa+ABZ*I_ESP_K3y4z_F2yz_aa;
    Double I_ESP_K2y5z_G2y2z_aa = I_ESP_L2y6z_F2yz_aa+ABZ*I_ESP_K2y5z_F2yz_aa;
    Double I_ESP_Ky6z_G2y2z_aa = I_ESP_Ly7z_F2yz_aa+ABZ*I_ESP_Ky6z_F2yz_aa;
    Double I_ESP_K7z_G2y2z_aa = I_ESP_L8z_F2yz_aa+ABZ*I_ESP_K7z_F2yz_aa;
    Double I_ESP_K7x_Gy3z_aa = I_ESP_L7xy_F3z_aa+ABY*I_ESP_K7x_F3z_aa;
    Double I_ESP_K6xy_Gy3z_aa = I_ESP_L6x2y_F3z_aa+ABY*I_ESP_K6xy_F3z_aa;
    Double I_ESP_K6xz_Gy3z_aa = I_ESP_L6xyz_F3z_aa+ABY*I_ESP_K6xz_F3z_aa;
    Double I_ESP_K5x2y_Gy3z_aa = I_ESP_L5x3y_F3z_aa+ABY*I_ESP_K5x2y_F3z_aa;
    Double I_ESP_K5xyz_Gy3z_aa = I_ESP_L5x2yz_F3z_aa+ABY*I_ESP_K5xyz_F3z_aa;
    Double I_ESP_K5x2z_Gy3z_aa = I_ESP_L5xy2z_F3z_aa+ABY*I_ESP_K5x2z_F3z_aa;
    Double I_ESP_K4x3y_Gy3z_aa = I_ESP_L4x4y_F3z_aa+ABY*I_ESP_K4x3y_F3z_aa;
    Double I_ESP_K4x2yz_Gy3z_aa = I_ESP_L4x3yz_F3z_aa+ABY*I_ESP_K4x2yz_F3z_aa;
    Double I_ESP_K4xy2z_Gy3z_aa = I_ESP_L4x2y2z_F3z_aa+ABY*I_ESP_K4xy2z_F3z_aa;
    Double I_ESP_K4x3z_Gy3z_aa = I_ESP_L4xy3z_F3z_aa+ABY*I_ESP_K4x3z_F3z_aa;
    Double I_ESP_K3x4y_Gy3z_aa = I_ESP_L3x5y_F3z_aa+ABY*I_ESP_K3x4y_F3z_aa;
    Double I_ESP_K3x3yz_Gy3z_aa = I_ESP_L3x4yz_F3z_aa+ABY*I_ESP_K3x3yz_F3z_aa;
    Double I_ESP_K3x2y2z_Gy3z_aa = I_ESP_L3x3y2z_F3z_aa+ABY*I_ESP_K3x2y2z_F3z_aa;
    Double I_ESP_K3xy3z_Gy3z_aa = I_ESP_L3x2y3z_F3z_aa+ABY*I_ESP_K3xy3z_F3z_aa;
    Double I_ESP_K3x4z_Gy3z_aa = I_ESP_L3xy4z_F3z_aa+ABY*I_ESP_K3x4z_F3z_aa;
    Double I_ESP_K2x5y_Gy3z_aa = I_ESP_L2x6y_F3z_aa+ABY*I_ESP_K2x5y_F3z_aa;
    Double I_ESP_K2x4yz_Gy3z_aa = I_ESP_L2x5yz_F3z_aa+ABY*I_ESP_K2x4yz_F3z_aa;
    Double I_ESP_K2x3y2z_Gy3z_aa = I_ESP_L2x4y2z_F3z_aa+ABY*I_ESP_K2x3y2z_F3z_aa;
    Double I_ESP_K2x2y3z_Gy3z_aa = I_ESP_L2x3y3z_F3z_aa+ABY*I_ESP_K2x2y3z_F3z_aa;
    Double I_ESP_K2xy4z_Gy3z_aa = I_ESP_L2x2y4z_F3z_aa+ABY*I_ESP_K2xy4z_F3z_aa;
    Double I_ESP_K2x5z_Gy3z_aa = I_ESP_L2xy5z_F3z_aa+ABY*I_ESP_K2x5z_F3z_aa;
    Double I_ESP_Kx6y_Gy3z_aa = I_ESP_Lx7y_F3z_aa+ABY*I_ESP_Kx6y_F3z_aa;
    Double I_ESP_Kx5yz_Gy3z_aa = I_ESP_Lx6yz_F3z_aa+ABY*I_ESP_Kx5yz_F3z_aa;
    Double I_ESP_Kx4y2z_Gy3z_aa = I_ESP_Lx5y2z_F3z_aa+ABY*I_ESP_Kx4y2z_F3z_aa;
    Double I_ESP_Kx3y3z_Gy3z_aa = I_ESP_Lx4y3z_F3z_aa+ABY*I_ESP_Kx3y3z_F3z_aa;
    Double I_ESP_Kx2y4z_Gy3z_aa = I_ESP_Lx3y4z_F3z_aa+ABY*I_ESP_Kx2y4z_F3z_aa;
    Double I_ESP_Kxy5z_Gy3z_aa = I_ESP_Lx2y5z_F3z_aa+ABY*I_ESP_Kxy5z_F3z_aa;
    Double I_ESP_Kx6z_Gy3z_aa = I_ESP_Lxy6z_F3z_aa+ABY*I_ESP_Kx6z_F3z_aa;
    Double I_ESP_K7y_Gy3z_aa = I_ESP_L8y_F3z_aa+ABY*I_ESP_K7y_F3z_aa;
    Double I_ESP_K6yz_Gy3z_aa = I_ESP_L7yz_F3z_aa+ABY*I_ESP_K6yz_F3z_aa;
    Double I_ESP_K5y2z_Gy3z_aa = I_ESP_L6y2z_F3z_aa+ABY*I_ESP_K5y2z_F3z_aa;
    Double I_ESP_K4y3z_Gy3z_aa = I_ESP_L5y3z_F3z_aa+ABY*I_ESP_K4y3z_F3z_aa;
    Double I_ESP_K3y4z_Gy3z_aa = I_ESP_L4y4z_F3z_aa+ABY*I_ESP_K3y4z_F3z_aa;
    Double I_ESP_K2y5z_Gy3z_aa = I_ESP_L3y5z_F3z_aa+ABY*I_ESP_K2y5z_F3z_aa;
    Double I_ESP_Ky6z_Gy3z_aa = I_ESP_L2y6z_F3z_aa+ABY*I_ESP_Ky6z_F3z_aa;
    Double I_ESP_K7z_Gy3z_aa = I_ESP_Ly7z_F3z_aa+ABY*I_ESP_K7z_F3z_aa;
    Double I_ESP_K7x_G4z_aa = I_ESP_L7xz_F3z_aa+ABZ*I_ESP_K7x_F3z_aa;
    Double I_ESP_K6xy_G4z_aa = I_ESP_L6xyz_F3z_aa+ABZ*I_ESP_K6xy_F3z_aa;
    Double I_ESP_K6xz_G4z_aa = I_ESP_L6x2z_F3z_aa+ABZ*I_ESP_K6xz_F3z_aa;
    Double I_ESP_K5x2y_G4z_aa = I_ESP_L5x2yz_F3z_aa+ABZ*I_ESP_K5x2y_F3z_aa;
    Double I_ESP_K5xyz_G4z_aa = I_ESP_L5xy2z_F3z_aa+ABZ*I_ESP_K5xyz_F3z_aa;
    Double I_ESP_K5x2z_G4z_aa = I_ESP_L5x3z_F3z_aa+ABZ*I_ESP_K5x2z_F3z_aa;
    Double I_ESP_K4x3y_G4z_aa = I_ESP_L4x3yz_F3z_aa+ABZ*I_ESP_K4x3y_F3z_aa;
    Double I_ESP_K4x2yz_G4z_aa = I_ESP_L4x2y2z_F3z_aa+ABZ*I_ESP_K4x2yz_F3z_aa;
    Double I_ESP_K4xy2z_G4z_aa = I_ESP_L4xy3z_F3z_aa+ABZ*I_ESP_K4xy2z_F3z_aa;
    Double I_ESP_K4x3z_G4z_aa = I_ESP_L4x4z_F3z_aa+ABZ*I_ESP_K4x3z_F3z_aa;
    Double I_ESP_K3x4y_G4z_aa = I_ESP_L3x4yz_F3z_aa+ABZ*I_ESP_K3x4y_F3z_aa;
    Double I_ESP_K3x3yz_G4z_aa = I_ESP_L3x3y2z_F3z_aa+ABZ*I_ESP_K3x3yz_F3z_aa;
    Double I_ESP_K3x2y2z_G4z_aa = I_ESP_L3x2y3z_F3z_aa+ABZ*I_ESP_K3x2y2z_F3z_aa;
    Double I_ESP_K3xy3z_G4z_aa = I_ESP_L3xy4z_F3z_aa+ABZ*I_ESP_K3xy3z_F3z_aa;
    Double I_ESP_K3x4z_G4z_aa = I_ESP_L3x5z_F3z_aa+ABZ*I_ESP_K3x4z_F3z_aa;
    Double I_ESP_K2x5y_G4z_aa = I_ESP_L2x5yz_F3z_aa+ABZ*I_ESP_K2x5y_F3z_aa;
    Double I_ESP_K2x4yz_G4z_aa = I_ESP_L2x4y2z_F3z_aa+ABZ*I_ESP_K2x4yz_F3z_aa;
    Double I_ESP_K2x3y2z_G4z_aa = I_ESP_L2x3y3z_F3z_aa+ABZ*I_ESP_K2x3y2z_F3z_aa;
    Double I_ESP_K2x2y3z_G4z_aa = I_ESP_L2x2y4z_F3z_aa+ABZ*I_ESP_K2x2y3z_F3z_aa;
    Double I_ESP_K2xy4z_G4z_aa = I_ESP_L2xy5z_F3z_aa+ABZ*I_ESP_K2xy4z_F3z_aa;
    Double I_ESP_K2x5z_G4z_aa = I_ESP_L2x6z_F3z_aa+ABZ*I_ESP_K2x5z_F3z_aa;
    Double I_ESP_Kx6y_G4z_aa = I_ESP_Lx6yz_F3z_aa+ABZ*I_ESP_Kx6y_F3z_aa;
    Double I_ESP_Kx5yz_G4z_aa = I_ESP_Lx5y2z_F3z_aa+ABZ*I_ESP_Kx5yz_F3z_aa;
    Double I_ESP_Kx4y2z_G4z_aa = I_ESP_Lx4y3z_F3z_aa+ABZ*I_ESP_Kx4y2z_F3z_aa;
    Double I_ESP_Kx3y3z_G4z_aa = I_ESP_Lx3y4z_F3z_aa+ABZ*I_ESP_Kx3y3z_F3z_aa;
    Double I_ESP_Kx2y4z_G4z_aa = I_ESP_Lx2y5z_F3z_aa+ABZ*I_ESP_Kx2y4z_F3z_aa;
    Double I_ESP_Kxy5z_G4z_aa = I_ESP_Lxy6z_F3z_aa+ABZ*I_ESP_Kxy5z_F3z_aa;
    Double I_ESP_Kx6z_G4z_aa = I_ESP_Lx7z_F3z_aa+ABZ*I_ESP_Kx6z_F3z_aa;
    Double I_ESP_K7y_G4z_aa = I_ESP_L7yz_F3z_aa+ABZ*I_ESP_K7y_F3z_aa;
    Double I_ESP_K6yz_G4z_aa = I_ESP_L6y2z_F3z_aa+ABZ*I_ESP_K6yz_F3z_aa;
    Double I_ESP_K5y2z_G4z_aa = I_ESP_L5y3z_F3z_aa+ABZ*I_ESP_K5y2z_F3z_aa;
    Double I_ESP_K4y3z_G4z_aa = I_ESP_L4y4z_F3z_aa+ABZ*I_ESP_K4y3z_F3z_aa;
    Double I_ESP_K3y4z_G4z_aa = I_ESP_L3y5z_F3z_aa+ABZ*I_ESP_K3y4z_F3z_aa;
    Double I_ESP_K2y5z_G4z_aa = I_ESP_L2y6z_F3z_aa+ABZ*I_ESP_K2y5z_F3z_aa;
    Double I_ESP_Ky6z_G4z_aa = I_ESP_Ly7z_F3z_aa+ABZ*I_ESP_Ky6z_F3z_aa;
    Double I_ESP_K7z_G4z_aa = I_ESP_L8z_F3z_aa+ABZ*I_ESP_K7z_F3z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_G_aa
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_F_G
     ************************************************************/
    abcd[iGrid*1890+0] = 4.0E0*I_ESP_K7x_G4x_aa-2.0E0*5*I_ESP_H5x_G4x_a-2.0E0*6*I_ESP_H5x_G4x_a+5*4*I_ESP_F3x_G4x;
    abcd[iGrid*1890+1] = 4.0E0*I_ESP_K6xy_G4x_aa-2.0E0*4*I_ESP_H4xy_G4x_a-2.0E0*5*I_ESP_H4xy_G4x_a+4*3*I_ESP_F2xy_G4x;
    abcd[iGrid*1890+2] = 4.0E0*I_ESP_K6xz_G4x_aa-2.0E0*4*I_ESP_H4xz_G4x_a-2.0E0*5*I_ESP_H4xz_G4x_a+4*3*I_ESP_F2xz_G4x;
    abcd[iGrid*1890+3] = 4.0E0*I_ESP_K5x2y_G4x_aa-2.0E0*3*I_ESP_H3x2y_G4x_a-2.0E0*4*I_ESP_H3x2y_G4x_a+3*2*I_ESP_Fx2y_G4x;
    abcd[iGrid*1890+4] = 4.0E0*I_ESP_K5xyz_G4x_aa-2.0E0*3*I_ESP_H3xyz_G4x_a-2.0E0*4*I_ESP_H3xyz_G4x_a+3*2*I_ESP_Fxyz_G4x;
    abcd[iGrid*1890+5] = 4.0E0*I_ESP_K5x2z_G4x_aa-2.0E0*3*I_ESP_H3x2z_G4x_a-2.0E0*4*I_ESP_H3x2z_G4x_a+3*2*I_ESP_Fx2z_G4x;
    abcd[iGrid*1890+6] = 4.0E0*I_ESP_K4x3y_G4x_aa-2.0E0*2*I_ESP_H2x3y_G4x_a-2.0E0*3*I_ESP_H2x3y_G4x_a+2*1*I_ESP_F3y_G4x;
    abcd[iGrid*1890+7] = 4.0E0*I_ESP_K4x2yz_G4x_aa-2.0E0*2*I_ESP_H2x2yz_G4x_a-2.0E0*3*I_ESP_H2x2yz_G4x_a+2*1*I_ESP_F2yz_G4x;
    abcd[iGrid*1890+8] = 4.0E0*I_ESP_K4xy2z_G4x_aa-2.0E0*2*I_ESP_H2xy2z_G4x_a-2.0E0*3*I_ESP_H2xy2z_G4x_a+2*1*I_ESP_Fy2z_G4x;
    abcd[iGrid*1890+9] = 4.0E0*I_ESP_K4x3z_G4x_aa-2.0E0*2*I_ESP_H2x3z_G4x_a-2.0E0*3*I_ESP_H2x3z_G4x_a+2*1*I_ESP_F3z_G4x;
    abcd[iGrid*1890+10] = 4.0E0*I_ESP_K3x4y_G4x_aa-2.0E0*1*I_ESP_Hx4y_G4x_a-2.0E0*2*I_ESP_Hx4y_G4x_a;
    abcd[iGrid*1890+11] = 4.0E0*I_ESP_K3x3yz_G4x_aa-2.0E0*1*I_ESP_Hx3yz_G4x_a-2.0E0*2*I_ESP_Hx3yz_G4x_a;
    abcd[iGrid*1890+12] = 4.0E0*I_ESP_K3x2y2z_G4x_aa-2.0E0*1*I_ESP_Hx2y2z_G4x_a-2.0E0*2*I_ESP_Hx2y2z_G4x_a;
    abcd[iGrid*1890+13] = 4.0E0*I_ESP_K3xy3z_G4x_aa-2.0E0*1*I_ESP_Hxy3z_G4x_a-2.0E0*2*I_ESP_Hxy3z_G4x_a;
    abcd[iGrid*1890+14] = 4.0E0*I_ESP_K3x4z_G4x_aa-2.0E0*1*I_ESP_Hx4z_G4x_a-2.0E0*2*I_ESP_Hx4z_G4x_a;
    abcd[iGrid*1890+15] = 4.0E0*I_ESP_K2x5y_G4x_aa-2.0E0*1*I_ESP_H5y_G4x_a;
    abcd[iGrid*1890+16] = 4.0E0*I_ESP_K2x4yz_G4x_aa-2.0E0*1*I_ESP_H4yz_G4x_a;
    abcd[iGrid*1890+17] = 4.0E0*I_ESP_K2x3y2z_G4x_aa-2.0E0*1*I_ESP_H3y2z_G4x_a;
    abcd[iGrid*1890+18] = 4.0E0*I_ESP_K2x2y3z_G4x_aa-2.0E0*1*I_ESP_H2y3z_G4x_a;
    abcd[iGrid*1890+19] = 4.0E0*I_ESP_K2xy4z_G4x_aa-2.0E0*1*I_ESP_Hy4z_G4x_a;
    abcd[iGrid*1890+20] = 4.0E0*I_ESP_K2x5z_G4x_aa-2.0E0*1*I_ESP_H5z_G4x_a;
    abcd[iGrid*1890+21] = 4.0E0*I_ESP_K7x_G3xy_aa-2.0E0*5*I_ESP_H5x_G3xy_a-2.0E0*6*I_ESP_H5x_G3xy_a+5*4*I_ESP_F3x_G3xy;
    abcd[iGrid*1890+22] = 4.0E0*I_ESP_K6xy_G3xy_aa-2.0E0*4*I_ESP_H4xy_G3xy_a-2.0E0*5*I_ESP_H4xy_G3xy_a+4*3*I_ESP_F2xy_G3xy;
    abcd[iGrid*1890+23] = 4.0E0*I_ESP_K6xz_G3xy_aa-2.0E0*4*I_ESP_H4xz_G3xy_a-2.0E0*5*I_ESP_H4xz_G3xy_a+4*3*I_ESP_F2xz_G3xy;
    abcd[iGrid*1890+24] = 4.0E0*I_ESP_K5x2y_G3xy_aa-2.0E0*3*I_ESP_H3x2y_G3xy_a-2.0E0*4*I_ESP_H3x2y_G3xy_a+3*2*I_ESP_Fx2y_G3xy;
    abcd[iGrid*1890+25] = 4.0E0*I_ESP_K5xyz_G3xy_aa-2.0E0*3*I_ESP_H3xyz_G3xy_a-2.0E0*4*I_ESP_H3xyz_G3xy_a+3*2*I_ESP_Fxyz_G3xy;
    abcd[iGrid*1890+26] = 4.0E0*I_ESP_K5x2z_G3xy_aa-2.0E0*3*I_ESP_H3x2z_G3xy_a-2.0E0*4*I_ESP_H3x2z_G3xy_a+3*2*I_ESP_Fx2z_G3xy;
    abcd[iGrid*1890+27] = 4.0E0*I_ESP_K4x3y_G3xy_aa-2.0E0*2*I_ESP_H2x3y_G3xy_a-2.0E0*3*I_ESP_H2x3y_G3xy_a+2*1*I_ESP_F3y_G3xy;
    abcd[iGrid*1890+28] = 4.0E0*I_ESP_K4x2yz_G3xy_aa-2.0E0*2*I_ESP_H2x2yz_G3xy_a-2.0E0*3*I_ESP_H2x2yz_G3xy_a+2*1*I_ESP_F2yz_G3xy;
    abcd[iGrid*1890+29] = 4.0E0*I_ESP_K4xy2z_G3xy_aa-2.0E0*2*I_ESP_H2xy2z_G3xy_a-2.0E0*3*I_ESP_H2xy2z_G3xy_a+2*1*I_ESP_Fy2z_G3xy;
    abcd[iGrid*1890+30] = 4.0E0*I_ESP_K4x3z_G3xy_aa-2.0E0*2*I_ESP_H2x3z_G3xy_a-2.0E0*3*I_ESP_H2x3z_G3xy_a+2*1*I_ESP_F3z_G3xy;
    abcd[iGrid*1890+31] = 4.0E0*I_ESP_K3x4y_G3xy_aa-2.0E0*1*I_ESP_Hx4y_G3xy_a-2.0E0*2*I_ESP_Hx4y_G3xy_a;
    abcd[iGrid*1890+32] = 4.0E0*I_ESP_K3x3yz_G3xy_aa-2.0E0*1*I_ESP_Hx3yz_G3xy_a-2.0E0*2*I_ESP_Hx3yz_G3xy_a;
    abcd[iGrid*1890+33] = 4.0E0*I_ESP_K3x2y2z_G3xy_aa-2.0E0*1*I_ESP_Hx2y2z_G3xy_a-2.0E0*2*I_ESP_Hx2y2z_G3xy_a;
    abcd[iGrid*1890+34] = 4.0E0*I_ESP_K3xy3z_G3xy_aa-2.0E0*1*I_ESP_Hxy3z_G3xy_a-2.0E0*2*I_ESP_Hxy3z_G3xy_a;
    abcd[iGrid*1890+35] = 4.0E0*I_ESP_K3x4z_G3xy_aa-2.0E0*1*I_ESP_Hx4z_G3xy_a-2.0E0*2*I_ESP_Hx4z_G3xy_a;
    abcd[iGrid*1890+36] = 4.0E0*I_ESP_K2x5y_G3xy_aa-2.0E0*1*I_ESP_H5y_G3xy_a;
    abcd[iGrid*1890+37] = 4.0E0*I_ESP_K2x4yz_G3xy_aa-2.0E0*1*I_ESP_H4yz_G3xy_a;
    abcd[iGrid*1890+38] = 4.0E0*I_ESP_K2x3y2z_G3xy_aa-2.0E0*1*I_ESP_H3y2z_G3xy_a;
    abcd[iGrid*1890+39] = 4.0E0*I_ESP_K2x2y3z_G3xy_aa-2.0E0*1*I_ESP_H2y3z_G3xy_a;
    abcd[iGrid*1890+40] = 4.0E0*I_ESP_K2xy4z_G3xy_aa-2.0E0*1*I_ESP_Hy4z_G3xy_a;
    abcd[iGrid*1890+41] = 4.0E0*I_ESP_K2x5z_G3xy_aa-2.0E0*1*I_ESP_H5z_G3xy_a;
    abcd[iGrid*1890+42] = 4.0E0*I_ESP_K7x_G3xz_aa-2.0E0*5*I_ESP_H5x_G3xz_a-2.0E0*6*I_ESP_H5x_G3xz_a+5*4*I_ESP_F3x_G3xz;
    abcd[iGrid*1890+43] = 4.0E0*I_ESP_K6xy_G3xz_aa-2.0E0*4*I_ESP_H4xy_G3xz_a-2.0E0*5*I_ESP_H4xy_G3xz_a+4*3*I_ESP_F2xy_G3xz;
    abcd[iGrid*1890+44] = 4.0E0*I_ESP_K6xz_G3xz_aa-2.0E0*4*I_ESP_H4xz_G3xz_a-2.0E0*5*I_ESP_H4xz_G3xz_a+4*3*I_ESP_F2xz_G3xz;
    abcd[iGrid*1890+45] = 4.0E0*I_ESP_K5x2y_G3xz_aa-2.0E0*3*I_ESP_H3x2y_G3xz_a-2.0E0*4*I_ESP_H3x2y_G3xz_a+3*2*I_ESP_Fx2y_G3xz;
    abcd[iGrid*1890+46] = 4.0E0*I_ESP_K5xyz_G3xz_aa-2.0E0*3*I_ESP_H3xyz_G3xz_a-2.0E0*4*I_ESP_H3xyz_G3xz_a+3*2*I_ESP_Fxyz_G3xz;
    abcd[iGrid*1890+47] = 4.0E0*I_ESP_K5x2z_G3xz_aa-2.0E0*3*I_ESP_H3x2z_G3xz_a-2.0E0*4*I_ESP_H3x2z_G3xz_a+3*2*I_ESP_Fx2z_G3xz;
    abcd[iGrid*1890+48] = 4.0E0*I_ESP_K4x3y_G3xz_aa-2.0E0*2*I_ESP_H2x3y_G3xz_a-2.0E0*3*I_ESP_H2x3y_G3xz_a+2*1*I_ESP_F3y_G3xz;
    abcd[iGrid*1890+49] = 4.0E0*I_ESP_K4x2yz_G3xz_aa-2.0E0*2*I_ESP_H2x2yz_G3xz_a-2.0E0*3*I_ESP_H2x2yz_G3xz_a+2*1*I_ESP_F2yz_G3xz;
    abcd[iGrid*1890+50] = 4.0E0*I_ESP_K4xy2z_G3xz_aa-2.0E0*2*I_ESP_H2xy2z_G3xz_a-2.0E0*3*I_ESP_H2xy2z_G3xz_a+2*1*I_ESP_Fy2z_G3xz;
    abcd[iGrid*1890+51] = 4.0E0*I_ESP_K4x3z_G3xz_aa-2.0E0*2*I_ESP_H2x3z_G3xz_a-2.0E0*3*I_ESP_H2x3z_G3xz_a+2*1*I_ESP_F3z_G3xz;
    abcd[iGrid*1890+52] = 4.0E0*I_ESP_K3x4y_G3xz_aa-2.0E0*1*I_ESP_Hx4y_G3xz_a-2.0E0*2*I_ESP_Hx4y_G3xz_a;
    abcd[iGrid*1890+53] = 4.0E0*I_ESP_K3x3yz_G3xz_aa-2.0E0*1*I_ESP_Hx3yz_G3xz_a-2.0E0*2*I_ESP_Hx3yz_G3xz_a;
    abcd[iGrid*1890+54] = 4.0E0*I_ESP_K3x2y2z_G3xz_aa-2.0E0*1*I_ESP_Hx2y2z_G3xz_a-2.0E0*2*I_ESP_Hx2y2z_G3xz_a;
    abcd[iGrid*1890+55] = 4.0E0*I_ESP_K3xy3z_G3xz_aa-2.0E0*1*I_ESP_Hxy3z_G3xz_a-2.0E0*2*I_ESP_Hxy3z_G3xz_a;
    abcd[iGrid*1890+56] = 4.0E0*I_ESP_K3x4z_G3xz_aa-2.0E0*1*I_ESP_Hx4z_G3xz_a-2.0E0*2*I_ESP_Hx4z_G3xz_a;
    abcd[iGrid*1890+57] = 4.0E0*I_ESP_K2x5y_G3xz_aa-2.0E0*1*I_ESP_H5y_G3xz_a;
    abcd[iGrid*1890+58] = 4.0E0*I_ESP_K2x4yz_G3xz_aa-2.0E0*1*I_ESP_H4yz_G3xz_a;
    abcd[iGrid*1890+59] = 4.0E0*I_ESP_K2x3y2z_G3xz_aa-2.0E0*1*I_ESP_H3y2z_G3xz_a;
    abcd[iGrid*1890+60] = 4.0E0*I_ESP_K2x2y3z_G3xz_aa-2.0E0*1*I_ESP_H2y3z_G3xz_a;
    abcd[iGrid*1890+61] = 4.0E0*I_ESP_K2xy4z_G3xz_aa-2.0E0*1*I_ESP_Hy4z_G3xz_a;
    abcd[iGrid*1890+62] = 4.0E0*I_ESP_K2x5z_G3xz_aa-2.0E0*1*I_ESP_H5z_G3xz_a;
    abcd[iGrid*1890+63] = 4.0E0*I_ESP_K7x_G2x2y_aa-2.0E0*5*I_ESP_H5x_G2x2y_a-2.0E0*6*I_ESP_H5x_G2x2y_a+5*4*I_ESP_F3x_G2x2y;
    abcd[iGrid*1890+64] = 4.0E0*I_ESP_K6xy_G2x2y_aa-2.0E0*4*I_ESP_H4xy_G2x2y_a-2.0E0*5*I_ESP_H4xy_G2x2y_a+4*3*I_ESP_F2xy_G2x2y;
    abcd[iGrid*1890+65] = 4.0E0*I_ESP_K6xz_G2x2y_aa-2.0E0*4*I_ESP_H4xz_G2x2y_a-2.0E0*5*I_ESP_H4xz_G2x2y_a+4*3*I_ESP_F2xz_G2x2y;
    abcd[iGrid*1890+66] = 4.0E0*I_ESP_K5x2y_G2x2y_aa-2.0E0*3*I_ESP_H3x2y_G2x2y_a-2.0E0*4*I_ESP_H3x2y_G2x2y_a+3*2*I_ESP_Fx2y_G2x2y;
    abcd[iGrid*1890+67] = 4.0E0*I_ESP_K5xyz_G2x2y_aa-2.0E0*3*I_ESP_H3xyz_G2x2y_a-2.0E0*4*I_ESP_H3xyz_G2x2y_a+3*2*I_ESP_Fxyz_G2x2y;
    abcd[iGrid*1890+68] = 4.0E0*I_ESP_K5x2z_G2x2y_aa-2.0E0*3*I_ESP_H3x2z_G2x2y_a-2.0E0*4*I_ESP_H3x2z_G2x2y_a+3*2*I_ESP_Fx2z_G2x2y;
    abcd[iGrid*1890+69] = 4.0E0*I_ESP_K4x3y_G2x2y_aa-2.0E0*2*I_ESP_H2x3y_G2x2y_a-2.0E0*3*I_ESP_H2x3y_G2x2y_a+2*1*I_ESP_F3y_G2x2y;
    abcd[iGrid*1890+70] = 4.0E0*I_ESP_K4x2yz_G2x2y_aa-2.0E0*2*I_ESP_H2x2yz_G2x2y_a-2.0E0*3*I_ESP_H2x2yz_G2x2y_a+2*1*I_ESP_F2yz_G2x2y;
    abcd[iGrid*1890+71] = 4.0E0*I_ESP_K4xy2z_G2x2y_aa-2.0E0*2*I_ESP_H2xy2z_G2x2y_a-2.0E0*3*I_ESP_H2xy2z_G2x2y_a+2*1*I_ESP_Fy2z_G2x2y;
    abcd[iGrid*1890+72] = 4.0E0*I_ESP_K4x3z_G2x2y_aa-2.0E0*2*I_ESP_H2x3z_G2x2y_a-2.0E0*3*I_ESP_H2x3z_G2x2y_a+2*1*I_ESP_F3z_G2x2y;
    abcd[iGrid*1890+73] = 4.0E0*I_ESP_K3x4y_G2x2y_aa-2.0E0*1*I_ESP_Hx4y_G2x2y_a-2.0E0*2*I_ESP_Hx4y_G2x2y_a;
    abcd[iGrid*1890+74] = 4.0E0*I_ESP_K3x3yz_G2x2y_aa-2.0E0*1*I_ESP_Hx3yz_G2x2y_a-2.0E0*2*I_ESP_Hx3yz_G2x2y_a;
    abcd[iGrid*1890+75] = 4.0E0*I_ESP_K3x2y2z_G2x2y_aa-2.0E0*1*I_ESP_Hx2y2z_G2x2y_a-2.0E0*2*I_ESP_Hx2y2z_G2x2y_a;
    abcd[iGrid*1890+76] = 4.0E0*I_ESP_K3xy3z_G2x2y_aa-2.0E0*1*I_ESP_Hxy3z_G2x2y_a-2.0E0*2*I_ESP_Hxy3z_G2x2y_a;
    abcd[iGrid*1890+77] = 4.0E0*I_ESP_K3x4z_G2x2y_aa-2.0E0*1*I_ESP_Hx4z_G2x2y_a-2.0E0*2*I_ESP_Hx4z_G2x2y_a;
    abcd[iGrid*1890+78] = 4.0E0*I_ESP_K2x5y_G2x2y_aa-2.0E0*1*I_ESP_H5y_G2x2y_a;
    abcd[iGrid*1890+79] = 4.0E0*I_ESP_K2x4yz_G2x2y_aa-2.0E0*1*I_ESP_H4yz_G2x2y_a;
    abcd[iGrid*1890+80] = 4.0E0*I_ESP_K2x3y2z_G2x2y_aa-2.0E0*1*I_ESP_H3y2z_G2x2y_a;
    abcd[iGrid*1890+81] = 4.0E0*I_ESP_K2x2y3z_G2x2y_aa-2.0E0*1*I_ESP_H2y3z_G2x2y_a;
    abcd[iGrid*1890+82] = 4.0E0*I_ESP_K2xy4z_G2x2y_aa-2.0E0*1*I_ESP_Hy4z_G2x2y_a;
    abcd[iGrid*1890+83] = 4.0E0*I_ESP_K2x5z_G2x2y_aa-2.0E0*1*I_ESP_H5z_G2x2y_a;
    abcd[iGrid*1890+84] = 4.0E0*I_ESP_K7x_G2xyz_aa-2.0E0*5*I_ESP_H5x_G2xyz_a-2.0E0*6*I_ESP_H5x_G2xyz_a+5*4*I_ESP_F3x_G2xyz;
    abcd[iGrid*1890+85] = 4.0E0*I_ESP_K6xy_G2xyz_aa-2.0E0*4*I_ESP_H4xy_G2xyz_a-2.0E0*5*I_ESP_H4xy_G2xyz_a+4*3*I_ESP_F2xy_G2xyz;
    abcd[iGrid*1890+86] = 4.0E0*I_ESP_K6xz_G2xyz_aa-2.0E0*4*I_ESP_H4xz_G2xyz_a-2.0E0*5*I_ESP_H4xz_G2xyz_a+4*3*I_ESP_F2xz_G2xyz;
    abcd[iGrid*1890+87] = 4.0E0*I_ESP_K5x2y_G2xyz_aa-2.0E0*3*I_ESP_H3x2y_G2xyz_a-2.0E0*4*I_ESP_H3x2y_G2xyz_a+3*2*I_ESP_Fx2y_G2xyz;
    abcd[iGrid*1890+88] = 4.0E0*I_ESP_K5xyz_G2xyz_aa-2.0E0*3*I_ESP_H3xyz_G2xyz_a-2.0E0*4*I_ESP_H3xyz_G2xyz_a+3*2*I_ESP_Fxyz_G2xyz;
    abcd[iGrid*1890+89] = 4.0E0*I_ESP_K5x2z_G2xyz_aa-2.0E0*3*I_ESP_H3x2z_G2xyz_a-2.0E0*4*I_ESP_H3x2z_G2xyz_a+3*2*I_ESP_Fx2z_G2xyz;
    abcd[iGrid*1890+90] = 4.0E0*I_ESP_K4x3y_G2xyz_aa-2.0E0*2*I_ESP_H2x3y_G2xyz_a-2.0E0*3*I_ESP_H2x3y_G2xyz_a+2*1*I_ESP_F3y_G2xyz;
    abcd[iGrid*1890+91] = 4.0E0*I_ESP_K4x2yz_G2xyz_aa-2.0E0*2*I_ESP_H2x2yz_G2xyz_a-2.0E0*3*I_ESP_H2x2yz_G2xyz_a+2*1*I_ESP_F2yz_G2xyz;
    abcd[iGrid*1890+92] = 4.0E0*I_ESP_K4xy2z_G2xyz_aa-2.0E0*2*I_ESP_H2xy2z_G2xyz_a-2.0E0*3*I_ESP_H2xy2z_G2xyz_a+2*1*I_ESP_Fy2z_G2xyz;
    abcd[iGrid*1890+93] = 4.0E0*I_ESP_K4x3z_G2xyz_aa-2.0E0*2*I_ESP_H2x3z_G2xyz_a-2.0E0*3*I_ESP_H2x3z_G2xyz_a+2*1*I_ESP_F3z_G2xyz;
    abcd[iGrid*1890+94] = 4.0E0*I_ESP_K3x4y_G2xyz_aa-2.0E0*1*I_ESP_Hx4y_G2xyz_a-2.0E0*2*I_ESP_Hx4y_G2xyz_a;
    abcd[iGrid*1890+95] = 4.0E0*I_ESP_K3x3yz_G2xyz_aa-2.0E0*1*I_ESP_Hx3yz_G2xyz_a-2.0E0*2*I_ESP_Hx3yz_G2xyz_a;
    abcd[iGrid*1890+96] = 4.0E0*I_ESP_K3x2y2z_G2xyz_aa-2.0E0*1*I_ESP_Hx2y2z_G2xyz_a-2.0E0*2*I_ESP_Hx2y2z_G2xyz_a;
    abcd[iGrid*1890+97] = 4.0E0*I_ESP_K3xy3z_G2xyz_aa-2.0E0*1*I_ESP_Hxy3z_G2xyz_a-2.0E0*2*I_ESP_Hxy3z_G2xyz_a;
    abcd[iGrid*1890+98] = 4.0E0*I_ESP_K3x4z_G2xyz_aa-2.0E0*1*I_ESP_Hx4z_G2xyz_a-2.0E0*2*I_ESP_Hx4z_G2xyz_a;
    abcd[iGrid*1890+99] = 4.0E0*I_ESP_K2x5y_G2xyz_aa-2.0E0*1*I_ESP_H5y_G2xyz_a;
    abcd[iGrid*1890+100] = 4.0E0*I_ESP_K2x4yz_G2xyz_aa-2.0E0*1*I_ESP_H4yz_G2xyz_a;
    abcd[iGrid*1890+101] = 4.0E0*I_ESP_K2x3y2z_G2xyz_aa-2.0E0*1*I_ESP_H3y2z_G2xyz_a;
    abcd[iGrid*1890+102] = 4.0E0*I_ESP_K2x2y3z_G2xyz_aa-2.0E0*1*I_ESP_H2y3z_G2xyz_a;
    abcd[iGrid*1890+103] = 4.0E0*I_ESP_K2xy4z_G2xyz_aa-2.0E0*1*I_ESP_Hy4z_G2xyz_a;
    abcd[iGrid*1890+104] = 4.0E0*I_ESP_K2x5z_G2xyz_aa-2.0E0*1*I_ESP_H5z_G2xyz_a;
    abcd[iGrid*1890+105] = 4.0E0*I_ESP_K7x_G2x2z_aa-2.0E0*5*I_ESP_H5x_G2x2z_a-2.0E0*6*I_ESP_H5x_G2x2z_a+5*4*I_ESP_F3x_G2x2z;
    abcd[iGrid*1890+106] = 4.0E0*I_ESP_K6xy_G2x2z_aa-2.0E0*4*I_ESP_H4xy_G2x2z_a-2.0E0*5*I_ESP_H4xy_G2x2z_a+4*3*I_ESP_F2xy_G2x2z;
    abcd[iGrid*1890+107] = 4.0E0*I_ESP_K6xz_G2x2z_aa-2.0E0*4*I_ESP_H4xz_G2x2z_a-2.0E0*5*I_ESP_H4xz_G2x2z_a+4*3*I_ESP_F2xz_G2x2z;
    abcd[iGrid*1890+108] = 4.0E0*I_ESP_K5x2y_G2x2z_aa-2.0E0*3*I_ESP_H3x2y_G2x2z_a-2.0E0*4*I_ESP_H3x2y_G2x2z_a+3*2*I_ESP_Fx2y_G2x2z;
    abcd[iGrid*1890+109] = 4.0E0*I_ESP_K5xyz_G2x2z_aa-2.0E0*3*I_ESP_H3xyz_G2x2z_a-2.0E0*4*I_ESP_H3xyz_G2x2z_a+3*2*I_ESP_Fxyz_G2x2z;
    abcd[iGrid*1890+110] = 4.0E0*I_ESP_K5x2z_G2x2z_aa-2.0E0*3*I_ESP_H3x2z_G2x2z_a-2.0E0*4*I_ESP_H3x2z_G2x2z_a+3*2*I_ESP_Fx2z_G2x2z;
    abcd[iGrid*1890+111] = 4.0E0*I_ESP_K4x3y_G2x2z_aa-2.0E0*2*I_ESP_H2x3y_G2x2z_a-2.0E0*3*I_ESP_H2x3y_G2x2z_a+2*1*I_ESP_F3y_G2x2z;
    abcd[iGrid*1890+112] = 4.0E0*I_ESP_K4x2yz_G2x2z_aa-2.0E0*2*I_ESP_H2x2yz_G2x2z_a-2.0E0*3*I_ESP_H2x2yz_G2x2z_a+2*1*I_ESP_F2yz_G2x2z;
    abcd[iGrid*1890+113] = 4.0E0*I_ESP_K4xy2z_G2x2z_aa-2.0E0*2*I_ESP_H2xy2z_G2x2z_a-2.0E0*3*I_ESP_H2xy2z_G2x2z_a+2*1*I_ESP_Fy2z_G2x2z;
    abcd[iGrid*1890+114] = 4.0E0*I_ESP_K4x3z_G2x2z_aa-2.0E0*2*I_ESP_H2x3z_G2x2z_a-2.0E0*3*I_ESP_H2x3z_G2x2z_a+2*1*I_ESP_F3z_G2x2z;
    abcd[iGrid*1890+115] = 4.0E0*I_ESP_K3x4y_G2x2z_aa-2.0E0*1*I_ESP_Hx4y_G2x2z_a-2.0E0*2*I_ESP_Hx4y_G2x2z_a;
    abcd[iGrid*1890+116] = 4.0E0*I_ESP_K3x3yz_G2x2z_aa-2.0E0*1*I_ESP_Hx3yz_G2x2z_a-2.0E0*2*I_ESP_Hx3yz_G2x2z_a;
    abcd[iGrid*1890+117] = 4.0E0*I_ESP_K3x2y2z_G2x2z_aa-2.0E0*1*I_ESP_Hx2y2z_G2x2z_a-2.0E0*2*I_ESP_Hx2y2z_G2x2z_a;
    abcd[iGrid*1890+118] = 4.0E0*I_ESP_K3xy3z_G2x2z_aa-2.0E0*1*I_ESP_Hxy3z_G2x2z_a-2.0E0*2*I_ESP_Hxy3z_G2x2z_a;
    abcd[iGrid*1890+119] = 4.0E0*I_ESP_K3x4z_G2x2z_aa-2.0E0*1*I_ESP_Hx4z_G2x2z_a-2.0E0*2*I_ESP_Hx4z_G2x2z_a;
    abcd[iGrid*1890+120] = 4.0E0*I_ESP_K2x5y_G2x2z_aa-2.0E0*1*I_ESP_H5y_G2x2z_a;
    abcd[iGrid*1890+121] = 4.0E0*I_ESP_K2x4yz_G2x2z_aa-2.0E0*1*I_ESP_H4yz_G2x2z_a;
    abcd[iGrid*1890+122] = 4.0E0*I_ESP_K2x3y2z_G2x2z_aa-2.0E0*1*I_ESP_H3y2z_G2x2z_a;
    abcd[iGrid*1890+123] = 4.0E0*I_ESP_K2x2y3z_G2x2z_aa-2.0E0*1*I_ESP_H2y3z_G2x2z_a;
    abcd[iGrid*1890+124] = 4.0E0*I_ESP_K2xy4z_G2x2z_aa-2.0E0*1*I_ESP_Hy4z_G2x2z_a;
    abcd[iGrid*1890+125] = 4.0E0*I_ESP_K2x5z_G2x2z_aa-2.0E0*1*I_ESP_H5z_G2x2z_a;
    abcd[iGrid*1890+126] = 4.0E0*I_ESP_K7x_Gx3y_aa-2.0E0*5*I_ESP_H5x_Gx3y_a-2.0E0*6*I_ESP_H5x_Gx3y_a+5*4*I_ESP_F3x_Gx3y;
    abcd[iGrid*1890+127] = 4.0E0*I_ESP_K6xy_Gx3y_aa-2.0E0*4*I_ESP_H4xy_Gx3y_a-2.0E0*5*I_ESP_H4xy_Gx3y_a+4*3*I_ESP_F2xy_Gx3y;
    abcd[iGrid*1890+128] = 4.0E0*I_ESP_K6xz_Gx3y_aa-2.0E0*4*I_ESP_H4xz_Gx3y_a-2.0E0*5*I_ESP_H4xz_Gx3y_a+4*3*I_ESP_F2xz_Gx3y;
    abcd[iGrid*1890+129] = 4.0E0*I_ESP_K5x2y_Gx3y_aa-2.0E0*3*I_ESP_H3x2y_Gx3y_a-2.0E0*4*I_ESP_H3x2y_Gx3y_a+3*2*I_ESP_Fx2y_Gx3y;
    abcd[iGrid*1890+130] = 4.0E0*I_ESP_K5xyz_Gx3y_aa-2.0E0*3*I_ESP_H3xyz_Gx3y_a-2.0E0*4*I_ESP_H3xyz_Gx3y_a+3*2*I_ESP_Fxyz_Gx3y;
    abcd[iGrid*1890+131] = 4.0E0*I_ESP_K5x2z_Gx3y_aa-2.0E0*3*I_ESP_H3x2z_Gx3y_a-2.0E0*4*I_ESP_H3x2z_Gx3y_a+3*2*I_ESP_Fx2z_Gx3y;
    abcd[iGrid*1890+132] = 4.0E0*I_ESP_K4x3y_Gx3y_aa-2.0E0*2*I_ESP_H2x3y_Gx3y_a-2.0E0*3*I_ESP_H2x3y_Gx3y_a+2*1*I_ESP_F3y_Gx3y;
    abcd[iGrid*1890+133] = 4.0E0*I_ESP_K4x2yz_Gx3y_aa-2.0E0*2*I_ESP_H2x2yz_Gx3y_a-2.0E0*3*I_ESP_H2x2yz_Gx3y_a+2*1*I_ESP_F2yz_Gx3y;
    abcd[iGrid*1890+134] = 4.0E0*I_ESP_K4xy2z_Gx3y_aa-2.0E0*2*I_ESP_H2xy2z_Gx3y_a-2.0E0*3*I_ESP_H2xy2z_Gx3y_a+2*1*I_ESP_Fy2z_Gx3y;
    abcd[iGrid*1890+135] = 4.0E0*I_ESP_K4x3z_Gx3y_aa-2.0E0*2*I_ESP_H2x3z_Gx3y_a-2.0E0*3*I_ESP_H2x3z_Gx3y_a+2*1*I_ESP_F3z_Gx3y;
    abcd[iGrid*1890+136] = 4.0E0*I_ESP_K3x4y_Gx3y_aa-2.0E0*1*I_ESP_Hx4y_Gx3y_a-2.0E0*2*I_ESP_Hx4y_Gx3y_a;
    abcd[iGrid*1890+137] = 4.0E0*I_ESP_K3x3yz_Gx3y_aa-2.0E0*1*I_ESP_Hx3yz_Gx3y_a-2.0E0*2*I_ESP_Hx3yz_Gx3y_a;
    abcd[iGrid*1890+138] = 4.0E0*I_ESP_K3x2y2z_Gx3y_aa-2.0E0*1*I_ESP_Hx2y2z_Gx3y_a-2.0E0*2*I_ESP_Hx2y2z_Gx3y_a;
    abcd[iGrid*1890+139] = 4.0E0*I_ESP_K3xy3z_Gx3y_aa-2.0E0*1*I_ESP_Hxy3z_Gx3y_a-2.0E0*2*I_ESP_Hxy3z_Gx3y_a;
    abcd[iGrid*1890+140] = 4.0E0*I_ESP_K3x4z_Gx3y_aa-2.0E0*1*I_ESP_Hx4z_Gx3y_a-2.0E0*2*I_ESP_Hx4z_Gx3y_a;
    abcd[iGrid*1890+141] = 4.0E0*I_ESP_K2x5y_Gx3y_aa-2.0E0*1*I_ESP_H5y_Gx3y_a;
    abcd[iGrid*1890+142] = 4.0E0*I_ESP_K2x4yz_Gx3y_aa-2.0E0*1*I_ESP_H4yz_Gx3y_a;
    abcd[iGrid*1890+143] = 4.0E0*I_ESP_K2x3y2z_Gx3y_aa-2.0E0*1*I_ESP_H3y2z_Gx3y_a;
    abcd[iGrid*1890+144] = 4.0E0*I_ESP_K2x2y3z_Gx3y_aa-2.0E0*1*I_ESP_H2y3z_Gx3y_a;
    abcd[iGrid*1890+145] = 4.0E0*I_ESP_K2xy4z_Gx3y_aa-2.0E0*1*I_ESP_Hy4z_Gx3y_a;
    abcd[iGrid*1890+146] = 4.0E0*I_ESP_K2x5z_Gx3y_aa-2.0E0*1*I_ESP_H5z_Gx3y_a;
    abcd[iGrid*1890+147] = 4.0E0*I_ESP_K7x_Gx2yz_aa-2.0E0*5*I_ESP_H5x_Gx2yz_a-2.0E0*6*I_ESP_H5x_Gx2yz_a+5*4*I_ESP_F3x_Gx2yz;
    abcd[iGrid*1890+148] = 4.0E0*I_ESP_K6xy_Gx2yz_aa-2.0E0*4*I_ESP_H4xy_Gx2yz_a-2.0E0*5*I_ESP_H4xy_Gx2yz_a+4*3*I_ESP_F2xy_Gx2yz;
    abcd[iGrid*1890+149] = 4.0E0*I_ESP_K6xz_Gx2yz_aa-2.0E0*4*I_ESP_H4xz_Gx2yz_a-2.0E0*5*I_ESP_H4xz_Gx2yz_a+4*3*I_ESP_F2xz_Gx2yz;
    abcd[iGrid*1890+150] = 4.0E0*I_ESP_K5x2y_Gx2yz_aa-2.0E0*3*I_ESP_H3x2y_Gx2yz_a-2.0E0*4*I_ESP_H3x2y_Gx2yz_a+3*2*I_ESP_Fx2y_Gx2yz;
    abcd[iGrid*1890+151] = 4.0E0*I_ESP_K5xyz_Gx2yz_aa-2.0E0*3*I_ESP_H3xyz_Gx2yz_a-2.0E0*4*I_ESP_H3xyz_Gx2yz_a+3*2*I_ESP_Fxyz_Gx2yz;
    abcd[iGrid*1890+152] = 4.0E0*I_ESP_K5x2z_Gx2yz_aa-2.0E0*3*I_ESP_H3x2z_Gx2yz_a-2.0E0*4*I_ESP_H3x2z_Gx2yz_a+3*2*I_ESP_Fx2z_Gx2yz;
    abcd[iGrid*1890+153] = 4.0E0*I_ESP_K4x3y_Gx2yz_aa-2.0E0*2*I_ESP_H2x3y_Gx2yz_a-2.0E0*3*I_ESP_H2x3y_Gx2yz_a+2*1*I_ESP_F3y_Gx2yz;
    abcd[iGrid*1890+154] = 4.0E0*I_ESP_K4x2yz_Gx2yz_aa-2.0E0*2*I_ESP_H2x2yz_Gx2yz_a-2.0E0*3*I_ESP_H2x2yz_Gx2yz_a+2*1*I_ESP_F2yz_Gx2yz;
    abcd[iGrid*1890+155] = 4.0E0*I_ESP_K4xy2z_Gx2yz_aa-2.0E0*2*I_ESP_H2xy2z_Gx2yz_a-2.0E0*3*I_ESP_H2xy2z_Gx2yz_a+2*1*I_ESP_Fy2z_Gx2yz;
    abcd[iGrid*1890+156] = 4.0E0*I_ESP_K4x3z_Gx2yz_aa-2.0E0*2*I_ESP_H2x3z_Gx2yz_a-2.0E0*3*I_ESP_H2x3z_Gx2yz_a+2*1*I_ESP_F3z_Gx2yz;
    abcd[iGrid*1890+157] = 4.0E0*I_ESP_K3x4y_Gx2yz_aa-2.0E0*1*I_ESP_Hx4y_Gx2yz_a-2.0E0*2*I_ESP_Hx4y_Gx2yz_a;
    abcd[iGrid*1890+158] = 4.0E0*I_ESP_K3x3yz_Gx2yz_aa-2.0E0*1*I_ESP_Hx3yz_Gx2yz_a-2.0E0*2*I_ESP_Hx3yz_Gx2yz_a;
    abcd[iGrid*1890+159] = 4.0E0*I_ESP_K3x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_Hx2y2z_Gx2yz_a-2.0E0*2*I_ESP_Hx2y2z_Gx2yz_a;
    abcd[iGrid*1890+160] = 4.0E0*I_ESP_K3xy3z_Gx2yz_aa-2.0E0*1*I_ESP_Hxy3z_Gx2yz_a-2.0E0*2*I_ESP_Hxy3z_Gx2yz_a;
    abcd[iGrid*1890+161] = 4.0E0*I_ESP_K3x4z_Gx2yz_aa-2.0E0*1*I_ESP_Hx4z_Gx2yz_a-2.0E0*2*I_ESP_Hx4z_Gx2yz_a;
    abcd[iGrid*1890+162] = 4.0E0*I_ESP_K2x5y_Gx2yz_aa-2.0E0*1*I_ESP_H5y_Gx2yz_a;
    abcd[iGrid*1890+163] = 4.0E0*I_ESP_K2x4yz_Gx2yz_aa-2.0E0*1*I_ESP_H4yz_Gx2yz_a;
    abcd[iGrid*1890+164] = 4.0E0*I_ESP_K2x3y2z_Gx2yz_aa-2.0E0*1*I_ESP_H3y2z_Gx2yz_a;
    abcd[iGrid*1890+165] = 4.0E0*I_ESP_K2x2y3z_Gx2yz_aa-2.0E0*1*I_ESP_H2y3z_Gx2yz_a;
    abcd[iGrid*1890+166] = 4.0E0*I_ESP_K2xy4z_Gx2yz_aa-2.0E0*1*I_ESP_Hy4z_Gx2yz_a;
    abcd[iGrid*1890+167] = 4.0E0*I_ESP_K2x5z_Gx2yz_aa-2.0E0*1*I_ESP_H5z_Gx2yz_a;
    abcd[iGrid*1890+168] = 4.0E0*I_ESP_K7x_Gxy2z_aa-2.0E0*5*I_ESP_H5x_Gxy2z_a-2.0E0*6*I_ESP_H5x_Gxy2z_a+5*4*I_ESP_F3x_Gxy2z;
    abcd[iGrid*1890+169] = 4.0E0*I_ESP_K6xy_Gxy2z_aa-2.0E0*4*I_ESP_H4xy_Gxy2z_a-2.0E0*5*I_ESP_H4xy_Gxy2z_a+4*3*I_ESP_F2xy_Gxy2z;
    abcd[iGrid*1890+170] = 4.0E0*I_ESP_K6xz_Gxy2z_aa-2.0E0*4*I_ESP_H4xz_Gxy2z_a-2.0E0*5*I_ESP_H4xz_Gxy2z_a+4*3*I_ESP_F2xz_Gxy2z;
    abcd[iGrid*1890+171] = 4.0E0*I_ESP_K5x2y_Gxy2z_aa-2.0E0*3*I_ESP_H3x2y_Gxy2z_a-2.0E0*4*I_ESP_H3x2y_Gxy2z_a+3*2*I_ESP_Fx2y_Gxy2z;
    abcd[iGrid*1890+172] = 4.0E0*I_ESP_K5xyz_Gxy2z_aa-2.0E0*3*I_ESP_H3xyz_Gxy2z_a-2.0E0*4*I_ESP_H3xyz_Gxy2z_a+3*2*I_ESP_Fxyz_Gxy2z;
    abcd[iGrid*1890+173] = 4.0E0*I_ESP_K5x2z_Gxy2z_aa-2.0E0*3*I_ESP_H3x2z_Gxy2z_a-2.0E0*4*I_ESP_H3x2z_Gxy2z_a+3*2*I_ESP_Fx2z_Gxy2z;
    abcd[iGrid*1890+174] = 4.0E0*I_ESP_K4x3y_Gxy2z_aa-2.0E0*2*I_ESP_H2x3y_Gxy2z_a-2.0E0*3*I_ESP_H2x3y_Gxy2z_a+2*1*I_ESP_F3y_Gxy2z;
    abcd[iGrid*1890+175] = 4.0E0*I_ESP_K4x2yz_Gxy2z_aa-2.0E0*2*I_ESP_H2x2yz_Gxy2z_a-2.0E0*3*I_ESP_H2x2yz_Gxy2z_a+2*1*I_ESP_F2yz_Gxy2z;
    abcd[iGrid*1890+176] = 4.0E0*I_ESP_K4xy2z_Gxy2z_aa-2.0E0*2*I_ESP_H2xy2z_Gxy2z_a-2.0E0*3*I_ESP_H2xy2z_Gxy2z_a+2*1*I_ESP_Fy2z_Gxy2z;
    abcd[iGrid*1890+177] = 4.0E0*I_ESP_K4x3z_Gxy2z_aa-2.0E0*2*I_ESP_H2x3z_Gxy2z_a-2.0E0*3*I_ESP_H2x3z_Gxy2z_a+2*1*I_ESP_F3z_Gxy2z;
    abcd[iGrid*1890+178] = 4.0E0*I_ESP_K3x4y_Gxy2z_aa-2.0E0*1*I_ESP_Hx4y_Gxy2z_a-2.0E0*2*I_ESP_Hx4y_Gxy2z_a;
    abcd[iGrid*1890+179] = 4.0E0*I_ESP_K3x3yz_Gxy2z_aa-2.0E0*1*I_ESP_Hx3yz_Gxy2z_a-2.0E0*2*I_ESP_Hx3yz_Gxy2z_a;
    abcd[iGrid*1890+180] = 4.0E0*I_ESP_K3x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_Hx2y2z_Gxy2z_a-2.0E0*2*I_ESP_Hx2y2z_Gxy2z_a;
    abcd[iGrid*1890+181] = 4.0E0*I_ESP_K3xy3z_Gxy2z_aa-2.0E0*1*I_ESP_Hxy3z_Gxy2z_a-2.0E0*2*I_ESP_Hxy3z_Gxy2z_a;
    abcd[iGrid*1890+182] = 4.0E0*I_ESP_K3x4z_Gxy2z_aa-2.0E0*1*I_ESP_Hx4z_Gxy2z_a-2.0E0*2*I_ESP_Hx4z_Gxy2z_a;
    abcd[iGrid*1890+183] = 4.0E0*I_ESP_K2x5y_Gxy2z_aa-2.0E0*1*I_ESP_H5y_Gxy2z_a;
    abcd[iGrid*1890+184] = 4.0E0*I_ESP_K2x4yz_Gxy2z_aa-2.0E0*1*I_ESP_H4yz_Gxy2z_a;
    abcd[iGrid*1890+185] = 4.0E0*I_ESP_K2x3y2z_Gxy2z_aa-2.0E0*1*I_ESP_H3y2z_Gxy2z_a;
    abcd[iGrid*1890+186] = 4.0E0*I_ESP_K2x2y3z_Gxy2z_aa-2.0E0*1*I_ESP_H2y3z_Gxy2z_a;
    abcd[iGrid*1890+187] = 4.0E0*I_ESP_K2xy4z_Gxy2z_aa-2.0E0*1*I_ESP_Hy4z_Gxy2z_a;
    abcd[iGrid*1890+188] = 4.0E0*I_ESP_K2x5z_Gxy2z_aa-2.0E0*1*I_ESP_H5z_Gxy2z_a;
    abcd[iGrid*1890+189] = 4.0E0*I_ESP_K7x_Gx3z_aa-2.0E0*5*I_ESP_H5x_Gx3z_a-2.0E0*6*I_ESP_H5x_Gx3z_a+5*4*I_ESP_F3x_Gx3z;
    abcd[iGrid*1890+190] = 4.0E0*I_ESP_K6xy_Gx3z_aa-2.0E0*4*I_ESP_H4xy_Gx3z_a-2.0E0*5*I_ESP_H4xy_Gx3z_a+4*3*I_ESP_F2xy_Gx3z;
    abcd[iGrid*1890+191] = 4.0E0*I_ESP_K6xz_Gx3z_aa-2.0E0*4*I_ESP_H4xz_Gx3z_a-2.0E0*5*I_ESP_H4xz_Gx3z_a+4*3*I_ESP_F2xz_Gx3z;
    abcd[iGrid*1890+192] = 4.0E0*I_ESP_K5x2y_Gx3z_aa-2.0E0*3*I_ESP_H3x2y_Gx3z_a-2.0E0*4*I_ESP_H3x2y_Gx3z_a+3*2*I_ESP_Fx2y_Gx3z;
    abcd[iGrid*1890+193] = 4.0E0*I_ESP_K5xyz_Gx3z_aa-2.0E0*3*I_ESP_H3xyz_Gx3z_a-2.0E0*4*I_ESP_H3xyz_Gx3z_a+3*2*I_ESP_Fxyz_Gx3z;
    abcd[iGrid*1890+194] = 4.0E0*I_ESP_K5x2z_Gx3z_aa-2.0E0*3*I_ESP_H3x2z_Gx3z_a-2.0E0*4*I_ESP_H3x2z_Gx3z_a+3*2*I_ESP_Fx2z_Gx3z;
    abcd[iGrid*1890+195] = 4.0E0*I_ESP_K4x3y_Gx3z_aa-2.0E0*2*I_ESP_H2x3y_Gx3z_a-2.0E0*3*I_ESP_H2x3y_Gx3z_a+2*1*I_ESP_F3y_Gx3z;
    abcd[iGrid*1890+196] = 4.0E0*I_ESP_K4x2yz_Gx3z_aa-2.0E0*2*I_ESP_H2x2yz_Gx3z_a-2.0E0*3*I_ESP_H2x2yz_Gx3z_a+2*1*I_ESP_F2yz_Gx3z;
    abcd[iGrid*1890+197] = 4.0E0*I_ESP_K4xy2z_Gx3z_aa-2.0E0*2*I_ESP_H2xy2z_Gx3z_a-2.0E0*3*I_ESP_H2xy2z_Gx3z_a+2*1*I_ESP_Fy2z_Gx3z;
    abcd[iGrid*1890+198] = 4.0E0*I_ESP_K4x3z_Gx3z_aa-2.0E0*2*I_ESP_H2x3z_Gx3z_a-2.0E0*3*I_ESP_H2x3z_Gx3z_a+2*1*I_ESP_F3z_Gx3z;
    abcd[iGrid*1890+199] = 4.0E0*I_ESP_K3x4y_Gx3z_aa-2.0E0*1*I_ESP_Hx4y_Gx3z_a-2.0E0*2*I_ESP_Hx4y_Gx3z_a;
    abcd[iGrid*1890+200] = 4.0E0*I_ESP_K3x3yz_Gx3z_aa-2.0E0*1*I_ESP_Hx3yz_Gx3z_a-2.0E0*2*I_ESP_Hx3yz_Gx3z_a;
    abcd[iGrid*1890+201] = 4.0E0*I_ESP_K3x2y2z_Gx3z_aa-2.0E0*1*I_ESP_Hx2y2z_Gx3z_a-2.0E0*2*I_ESP_Hx2y2z_Gx3z_a;
    abcd[iGrid*1890+202] = 4.0E0*I_ESP_K3xy3z_Gx3z_aa-2.0E0*1*I_ESP_Hxy3z_Gx3z_a-2.0E0*2*I_ESP_Hxy3z_Gx3z_a;
    abcd[iGrid*1890+203] = 4.0E0*I_ESP_K3x4z_Gx3z_aa-2.0E0*1*I_ESP_Hx4z_Gx3z_a-2.0E0*2*I_ESP_Hx4z_Gx3z_a;
    abcd[iGrid*1890+204] = 4.0E0*I_ESP_K2x5y_Gx3z_aa-2.0E0*1*I_ESP_H5y_Gx3z_a;
    abcd[iGrid*1890+205] = 4.0E0*I_ESP_K2x4yz_Gx3z_aa-2.0E0*1*I_ESP_H4yz_Gx3z_a;
    abcd[iGrid*1890+206] = 4.0E0*I_ESP_K2x3y2z_Gx3z_aa-2.0E0*1*I_ESP_H3y2z_Gx3z_a;
    abcd[iGrid*1890+207] = 4.0E0*I_ESP_K2x2y3z_Gx3z_aa-2.0E0*1*I_ESP_H2y3z_Gx3z_a;
    abcd[iGrid*1890+208] = 4.0E0*I_ESP_K2xy4z_Gx3z_aa-2.0E0*1*I_ESP_Hy4z_Gx3z_a;
    abcd[iGrid*1890+209] = 4.0E0*I_ESP_K2x5z_Gx3z_aa-2.0E0*1*I_ESP_H5z_Gx3z_a;
    abcd[iGrid*1890+210] = 4.0E0*I_ESP_K7x_G4y_aa-2.0E0*5*I_ESP_H5x_G4y_a-2.0E0*6*I_ESP_H5x_G4y_a+5*4*I_ESP_F3x_G4y;
    abcd[iGrid*1890+211] = 4.0E0*I_ESP_K6xy_G4y_aa-2.0E0*4*I_ESP_H4xy_G4y_a-2.0E0*5*I_ESP_H4xy_G4y_a+4*3*I_ESP_F2xy_G4y;
    abcd[iGrid*1890+212] = 4.0E0*I_ESP_K6xz_G4y_aa-2.0E0*4*I_ESP_H4xz_G4y_a-2.0E0*5*I_ESP_H4xz_G4y_a+4*3*I_ESP_F2xz_G4y;
    abcd[iGrid*1890+213] = 4.0E0*I_ESP_K5x2y_G4y_aa-2.0E0*3*I_ESP_H3x2y_G4y_a-2.0E0*4*I_ESP_H3x2y_G4y_a+3*2*I_ESP_Fx2y_G4y;
    abcd[iGrid*1890+214] = 4.0E0*I_ESP_K5xyz_G4y_aa-2.0E0*3*I_ESP_H3xyz_G4y_a-2.0E0*4*I_ESP_H3xyz_G4y_a+3*2*I_ESP_Fxyz_G4y;
    abcd[iGrid*1890+215] = 4.0E0*I_ESP_K5x2z_G4y_aa-2.0E0*3*I_ESP_H3x2z_G4y_a-2.0E0*4*I_ESP_H3x2z_G4y_a+3*2*I_ESP_Fx2z_G4y;
    abcd[iGrid*1890+216] = 4.0E0*I_ESP_K4x3y_G4y_aa-2.0E0*2*I_ESP_H2x3y_G4y_a-2.0E0*3*I_ESP_H2x3y_G4y_a+2*1*I_ESP_F3y_G4y;
    abcd[iGrid*1890+217] = 4.0E0*I_ESP_K4x2yz_G4y_aa-2.0E0*2*I_ESP_H2x2yz_G4y_a-2.0E0*3*I_ESP_H2x2yz_G4y_a+2*1*I_ESP_F2yz_G4y;
    abcd[iGrid*1890+218] = 4.0E0*I_ESP_K4xy2z_G4y_aa-2.0E0*2*I_ESP_H2xy2z_G4y_a-2.0E0*3*I_ESP_H2xy2z_G4y_a+2*1*I_ESP_Fy2z_G4y;
    abcd[iGrid*1890+219] = 4.0E0*I_ESP_K4x3z_G4y_aa-2.0E0*2*I_ESP_H2x3z_G4y_a-2.0E0*3*I_ESP_H2x3z_G4y_a+2*1*I_ESP_F3z_G4y;
    abcd[iGrid*1890+220] = 4.0E0*I_ESP_K3x4y_G4y_aa-2.0E0*1*I_ESP_Hx4y_G4y_a-2.0E0*2*I_ESP_Hx4y_G4y_a;
    abcd[iGrid*1890+221] = 4.0E0*I_ESP_K3x3yz_G4y_aa-2.0E0*1*I_ESP_Hx3yz_G4y_a-2.0E0*2*I_ESP_Hx3yz_G4y_a;
    abcd[iGrid*1890+222] = 4.0E0*I_ESP_K3x2y2z_G4y_aa-2.0E0*1*I_ESP_Hx2y2z_G4y_a-2.0E0*2*I_ESP_Hx2y2z_G4y_a;
    abcd[iGrid*1890+223] = 4.0E0*I_ESP_K3xy3z_G4y_aa-2.0E0*1*I_ESP_Hxy3z_G4y_a-2.0E0*2*I_ESP_Hxy3z_G4y_a;
    abcd[iGrid*1890+224] = 4.0E0*I_ESP_K3x4z_G4y_aa-2.0E0*1*I_ESP_Hx4z_G4y_a-2.0E0*2*I_ESP_Hx4z_G4y_a;
    abcd[iGrid*1890+225] = 4.0E0*I_ESP_K2x5y_G4y_aa-2.0E0*1*I_ESP_H5y_G4y_a;
    abcd[iGrid*1890+226] = 4.0E0*I_ESP_K2x4yz_G4y_aa-2.0E0*1*I_ESP_H4yz_G4y_a;
    abcd[iGrid*1890+227] = 4.0E0*I_ESP_K2x3y2z_G4y_aa-2.0E0*1*I_ESP_H3y2z_G4y_a;
    abcd[iGrid*1890+228] = 4.0E0*I_ESP_K2x2y3z_G4y_aa-2.0E0*1*I_ESP_H2y3z_G4y_a;
    abcd[iGrid*1890+229] = 4.0E0*I_ESP_K2xy4z_G4y_aa-2.0E0*1*I_ESP_Hy4z_G4y_a;
    abcd[iGrid*1890+230] = 4.0E0*I_ESP_K2x5z_G4y_aa-2.0E0*1*I_ESP_H5z_G4y_a;
    abcd[iGrid*1890+231] = 4.0E0*I_ESP_K7x_G3yz_aa-2.0E0*5*I_ESP_H5x_G3yz_a-2.0E0*6*I_ESP_H5x_G3yz_a+5*4*I_ESP_F3x_G3yz;
    abcd[iGrid*1890+232] = 4.0E0*I_ESP_K6xy_G3yz_aa-2.0E0*4*I_ESP_H4xy_G3yz_a-2.0E0*5*I_ESP_H4xy_G3yz_a+4*3*I_ESP_F2xy_G3yz;
    abcd[iGrid*1890+233] = 4.0E0*I_ESP_K6xz_G3yz_aa-2.0E0*4*I_ESP_H4xz_G3yz_a-2.0E0*5*I_ESP_H4xz_G3yz_a+4*3*I_ESP_F2xz_G3yz;
    abcd[iGrid*1890+234] = 4.0E0*I_ESP_K5x2y_G3yz_aa-2.0E0*3*I_ESP_H3x2y_G3yz_a-2.0E0*4*I_ESP_H3x2y_G3yz_a+3*2*I_ESP_Fx2y_G3yz;
    abcd[iGrid*1890+235] = 4.0E0*I_ESP_K5xyz_G3yz_aa-2.0E0*3*I_ESP_H3xyz_G3yz_a-2.0E0*4*I_ESP_H3xyz_G3yz_a+3*2*I_ESP_Fxyz_G3yz;
    abcd[iGrid*1890+236] = 4.0E0*I_ESP_K5x2z_G3yz_aa-2.0E0*3*I_ESP_H3x2z_G3yz_a-2.0E0*4*I_ESP_H3x2z_G3yz_a+3*2*I_ESP_Fx2z_G3yz;
    abcd[iGrid*1890+237] = 4.0E0*I_ESP_K4x3y_G3yz_aa-2.0E0*2*I_ESP_H2x3y_G3yz_a-2.0E0*3*I_ESP_H2x3y_G3yz_a+2*1*I_ESP_F3y_G3yz;
    abcd[iGrid*1890+238] = 4.0E0*I_ESP_K4x2yz_G3yz_aa-2.0E0*2*I_ESP_H2x2yz_G3yz_a-2.0E0*3*I_ESP_H2x2yz_G3yz_a+2*1*I_ESP_F2yz_G3yz;
    abcd[iGrid*1890+239] = 4.0E0*I_ESP_K4xy2z_G3yz_aa-2.0E0*2*I_ESP_H2xy2z_G3yz_a-2.0E0*3*I_ESP_H2xy2z_G3yz_a+2*1*I_ESP_Fy2z_G3yz;
    abcd[iGrid*1890+240] = 4.0E0*I_ESP_K4x3z_G3yz_aa-2.0E0*2*I_ESP_H2x3z_G3yz_a-2.0E0*3*I_ESP_H2x3z_G3yz_a+2*1*I_ESP_F3z_G3yz;
    abcd[iGrid*1890+241] = 4.0E0*I_ESP_K3x4y_G3yz_aa-2.0E0*1*I_ESP_Hx4y_G3yz_a-2.0E0*2*I_ESP_Hx4y_G3yz_a;
    abcd[iGrid*1890+242] = 4.0E0*I_ESP_K3x3yz_G3yz_aa-2.0E0*1*I_ESP_Hx3yz_G3yz_a-2.0E0*2*I_ESP_Hx3yz_G3yz_a;
    abcd[iGrid*1890+243] = 4.0E0*I_ESP_K3x2y2z_G3yz_aa-2.0E0*1*I_ESP_Hx2y2z_G3yz_a-2.0E0*2*I_ESP_Hx2y2z_G3yz_a;
    abcd[iGrid*1890+244] = 4.0E0*I_ESP_K3xy3z_G3yz_aa-2.0E0*1*I_ESP_Hxy3z_G3yz_a-2.0E0*2*I_ESP_Hxy3z_G3yz_a;
    abcd[iGrid*1890+245] = 4.0E0*I_ESP_K3x4z_G3yz_aa-2.0E0*1*I_ESP_Hx4z_G3yz_a-2.0E0*2*I_ESP_Hx4z_G3yz_a;
    abcd[iGrid*1890+246] = 4.0E0*I_ESP_K2x5y_G3yz_aa-2.0E0*1*I_ESP_H5y_G3yz_a;
    abcd[iGrid*1890+247] = 4.0E0*I_ESP_K2x4yz_G3yz_aa-2.0E0*1*I_ESP_H4yz_G3yz_a;
    abcd[iGrid*1890+248] = 4.0E0*I_ESP_K2x3y2z_G3yz_aa-2.0E0*1*I_ESP_H3y2z_G3yz_a;
    abcd[iGrid*1890+249] = 4.0E0*I_ESP_K2x2y3z_G3yz_aa-2.0E0*1*I_ESP_H2y3z_G3yz_a;
    abcd[iGrid*1890+250] = 4.0E0*I_ESP_K2xy4z_G3yz_aa-2.0E0*1*I_ESP_Hy4z_G3yz_a;
    abcd[iGrid*1890+251] = 4.0E0*I_ESP_K2x5z_G3yz_aa-2.0E0*1*I_ESP_H5z_G3yz_a;
    abcd[iGrid*1890+252] = 4.0E0*I_ESP_K7x_G2y2z_aa-2.0E0*5*I_ESP_H5x_G2y2z_a-2.0E0*6*I_ESP_H5x_G2y2z_a+5*4*I_ESP_F3x_G2y2z;
    abcd[iGrid*1890+253] = 4.0E0*I_ESP_K6xy_G2y2z_aa-2.0E0*4*I_ESP_H4xy_G2y2z_a-2.0E0*5*I_ESP_H4xy_G2y2z_a+4*3*I_ESP_F2xy_G2y2z;
    abcd[iGrid*1890+254] = 4.0E0*I_ESP_K6xz_G2y2z_aa-2.0E0*4*I_ESP_H4xz_G2y2z_a-2.0E0*5*I_ESP_H4xz_G2y2z_a+4*3*I_ESP_F2xz_G2y2z;
    abcd[iGrid*1890+255] = 4.0E0*I_ESP_K5x2y_G2y2z_aa-2.0E0*3*I_ESP_H3x2y_G2y2z_a-2.0E0*4*I_ESP_H3x2y_G2y2z_a+3*2*I_ESP_Fx2y_G2y2z;
    abcd[iGrid*1890+256] = 4.0E0*I_ESP_K5xyz_G2y2z_aa-2.0E0*3*I_ESP_H3xyz_G2y2z_a-2.0E0*4*I_ESP_H3xyz_G2y2z_a+3*2*I_ESP_Fxyz_G2y2z;
    abcd[iGrid*1890+257] = 4.0E0*I_ESP_K5x2z_G2y2z_aa-2.0E0*3*I_ESP_H3x2z_G2y2z_a-2.0E0*4*I_ESP_H3x2z_G2y2z_a+3*2*I_ESP_Fx2z_G2y2z;
    abcd[iGrid*1890+258] = 4.0E0*I_ESP_K4x3y_G2y2z_aa-2.0E0*2*I_ESP_H2x3y_G2y2z_a-2.0E0*3*I_ESP_H2x3y_G2y2z_a+2*1*I_ESP_F3y_G2y2z;
    abcd[iGrid*1890+259] = 4.0E0*I_ESP_K4x2yz_G2y2z_aa-2.0E0*2*I_ESP_H2x2yz_G2y2z_a-2.0E0*3*I_ESP_H2x2yz_G2y2z_a+2*1*I_ESP_F2yz_G2y2z;
    abcd[iGrid*1890+260] = 4.0E0*I_ESP_K4xy2z_G2y2z_aa-2.0E0*2*I_ESP_H2xy2z_G2y2z_a-2.0E0*3*I_ESP_H2xy2z_G2y2z_a+2*1*I_ESP_Fy2z_G2y2z;
    abcd[iGrid*1890+261] = 4.0E0*I_ESP_K4x3z_G2y2z_aa-2.0E0*2*I_ESP_H2x3z_G2y2z_a-2.0E0*3*I_ESP_H2x3z_G2y2z_a+2*1*I_ESP_F3z_G2y2z;
    abcd[iGrid*1890+262] = 4.0E0*I_ESP_K3x4y_G2y2z_aa-2.0E0*1*I_ESP_Hx4y_G2y2z_a-2.0E0*2*I_ESP_Hx4y_G2y2z_a;
    abcd[iGrid*1890+263] = 4.0E0*I_ESP_K3x3yz_G2y2z_aa-2.0E0*1*I_ESP_Hx3yz_G2y2z_a-2.0E0*2*I_ESP_Hx3yz_G2y2z_a;
    abcd[iGrid*1890+264] = 4.0E0*I_ESP_K3x2y2z_G2y2z_aa-2.0E0*1*I_ESP_Hx2y2z_G2y2z_a-2.0E0*2*I_ESP_Hx2y2z_G2y2z_a;
    abcd[iGrid*1890+265] = 4.0E0*I_ESP_K3xy3z_G2y2z_aa-2.0E0*1*I_ESP_Hxy3z_G2y2z_a-2.0E0*2*I_ESP_Hxy3z_G2y2z_a;
    abcd[iGrid*1890+266] = 4.0E0*I_ESP_K3x4z_G2y2z_aa-2.0E0*1*I_ESP_Hx4z_G2y2z_a-2.0E0*2*I_ESP_Hx4z_G2y2z_a;
    abcd[iGrid*1890+267] = 4.0E0*I_ESP_K2x5y_G2y2z_aa-2.0E0*1*I_ESP_H5y_G2y2z_a;
    abcd[iGrid*1890+268] = 4.0E0*I_ESP_K2x4yz_G2y2z_aa-2.0E0*1*I_ESP_H4yz_G2y2z_a;
    abcd[iGrid*1890+269] = 4.0E0*I_ESP_K2x3y2z_G2y2z_aa-2.0E0*1*I_ESP_H3y2z_G2y2z_a;
    abcd[iGrid*1890+270] = 4.0E0*I_ESP_K2x2y3z_G2y2z_aa-2.0E0*1*I_ESP_H2y3z_G2y2z_a;
    abcd[iGrid*1890+271] = 4.0E0*I_ESP_K2xy4z_G2y2z_aa-2.0E0*1*I_ESP_Hy4z_G2y2z_a;
    abcd[iGrid*1890+272] = 4.0E0*I_ESP_K2x5z_G2y2z_aa-2.0E0*1*I_ESP_H5z_G2y2z_a;
    abcd[iGrid*1890+273] = 4.0E0*I_ESP_K7x_Gy3z_aa-2.0E0*5*I_ESP_H5x_Gy3z_a-2.0E0*6*I_ESP_H5x_Gy3z_a+5*4*I_ESP_F3x_Gy3z;
    abcd[iGrid*1890+274] = 4.0E0*I_ESP_K6xy_Gy3z_aa-2.0E0*4*I_ESP_H4xy_Gy3z_a-2.0E0*5*I_ESP_H4xy_Gy3z_a+4*3*I_ESP_F2xy_Gy3z;
    abcd[iGrid*1890+275] = 4.0E0*I_ESP_K6xz_Gy3z_aa-2.0E0*4*I_ESP_H4xz_Gy3z_a-2.0E0*5*I_ESP_H4xz_Gy3z_a+4*3*I_ESP_F2xz_Gy3z;
    abcd[iGrid*1890+276] = 4.0E0*I_ESP_K5x2y_Gy3z_aa-2.0E0*3*I_ESP_H3x2y_Gy3z_a-2.0E0*4*I_ESP_H3x2y_Gy3z_a+3*2*I_ESP_Fx2y_Gy3z;
    abcd[iGrid*1890+277] = 4.0E0*I_ESP_K5xyz_Gy3z_aa-2.0E0*3*I_ESP_H3xyz_Gy3z_a-2.0E0*4*I_ESP_H3xyz_Gy3z_a+3*2*I_ESP_Fxyz_Gy3z;
    abcd[iGrid*1890+278] = 4.0E0*I_ESP_K5x2z_Gy3z_aa-2.0E0*3*I_ESP_H3x2z_Gy3z_a-2.0E0*4*I_ESP_H3x2z_Gy3z_a+3*2*I_ESP_Fx2z_Gy3z;
    abcd[iGrid*1890+279] = 4.0E0*I_ESP_K4x3y_Gy3z_aa-2.0E0*2*I_ESP_H2x3y_Gy3z_a-2.0E0*3*I_ESP_H2x3y_Gy3z_a+2*1*I_ESP_F3y_Gy3z;
    abcd[iGrid*1890+280] = 4.0E0*I_ESP_K4x2yz_Gy3z_aa-2.0E0*2*I_ESP_H2x2yz_Gy3z_a-2.0E0*3*I_ESP_H2x2yz_Gy3z_a+2*1*I_ESP_F2yz_Gy3z;
    abcd[iGrid*1890+281] = 4.0E0*I_ESP_K4xy2z_Gy3z_aa-2.0E0*2*I_ESP_H2xy2z_Gy3z_a-2.0E0*3*I_ESP_H2xy2z_Gy3z_a+2*1*I_ESP_Fy2z_Gy3z;
    abcd[iGrid*1890+282] = 4.0E0*I_ESP_K4x3z_Gy3z_aa-2.0E0*2*I_ESP_H2x3z_Gy3z_a-2.0E0*3*I_ESP_H2x3z_Gy3z_a+2*1*I_ESP_F3z_Gy3z;
    abcd[iGrid*1890+283] = 4.0E0*I_ESP_K3x4y_Gy3z_aa-2.0E0*1*I_ESP_Hx4y_Gy3z_a-2.0E0*2*I_ESP_Hx4y_Gy3z_a;
    abcd[iGrid*1890+284] = 4.0E0*I_ESP_K3x3yz_Gy3z_aa-2.0E0*1*I_ESP_Hx3yz_Gy3z_a-2.0E0*2*I_ESP_Hx3yz_Gy3z_a;
    abcd[iGrid*1890+285] = 4.0E0*I_ESP_K3x2y2z_Gy3z_aa-2.0E0*1*I_ESP_Hx2y2z_Gy3z_a-2.0E0*2*I_ESP_Hx2y2z_Gy3z_a;
    abcd[iGrid*1890+286] = 4.0E0*I_ESP_K3xy3z_Gy3z_aa-2.0E0*1*I_ESP_Hxy3z_Gy3z_a-2.0E0*2*I_ESP_Hxy3z_Gy3z_a;
    abcd[iGrid*1890+287] = 4.0E0*I_ESP_K3x4z_Gy3z_aa-2.0E0*1*I_ESP_Hx4z_Gy3z_a-2.0E0*2*I_ESP_Hx4z_Gy3z_a;
    abcd[iGrid*1890+288] = 4.0E0*I_ESP_K2x5y_Gy3z_aa-2.0E0*1*I_ESP_H5y_Gy3z_a;
    abcd[iGrid*1890+289] = 4.0E0*I_ESP_K2x4yz_Gy3z_aa-2.0E0*1*I_ESP_H4yz_Gy3z_a;
    abcd[iGrid*1890+290] = 4.0E0*I_ESP_K2x3y2z_Gy3z_aa-2.0E0*1*I_ESP_H3y2z_Gy3z_a;
    abcd[iGrid*1890+291] = 4.0E0*I_ESP_K2x2y3z_Gy3z_aa-2.0E0*1*I_ESP_H2y3z_Gy3z_a;
    abcd[iGrid*1890+292] = 4.0E0*I_ESP_K2xy4z_Gy3z_aa-2.0E0*1*I_ESP_Hy4z_Gy3z_a;
    abcd[iGrid*1890+293] = 4.0E0*I_ESP_K2x5z_Gy3z_aa-2.0E0*1*I_ESP_H5z_Gy3z_a;
    abcd[iGrid*1890+294] = 4.0E0*I_ESP_K7x_G4z_aa-2.0E0*5*I_ESP_H5x_G4z_a-2.0E0*6*I_ESP_H5x_G4z_a+5*4*I_ESP_F3x_G4z;
    abcd[iGrid*1890+295] = 4.0E0*I_ESP_K6xy_G4z_aa-2.0E0*4*I_ESP_H4xy_G4z_a-2.0E0*5*I_ESP_H4xy_G4z_a+4*3*I_ESP_F2xy_G4z;
    abcd[iGrid*1890+296] = 4.0E0*I_ESP_K6xz_G4z_aa-2.0E0*4*I_ESP_H4xz_G4z_a-2.0E0*5*I_ESP_H4xz_G4z_a+4*3*I_ESP_F2xz_G4z;
    abcd[iGrid*1890+297] = 4.0E0*I_ESP_K5x2y_G4z_aa-2.0E0*3*I_ESP_H3x2y_G4z_a-2.0E0*4*I_ESP_H3x2y_G4z_a+3*2*I_ESP_Fx2y_G4z;
    abcd[iGrid*1890+298] = 4.0E0*I_ESP_K5xyz_G4z_aa-2.0E0*3*I_ESP_H3xyz_G4z_a-2.0E0*4*I_ESP_H3xyz_G4z_a+3*2*I_ESP_Fxyz_G4z;
    abcd[iGrid*1890+299] = 4.0E0*I_ESP_K5x2z_G4z_aa-2.0E0*3*I_ESP_H3x2z_G4z_a-2.0E0*4*I_ESP_H3x2z_G4z_a+3*2*I_ESP_Fx2z_G4z;
    abcd[iGrid*1890+300] = 4.0E0*I_ESP_K4x3y_G4z_aa-2.0E0*2*I_ESP_H2x3y_G4z_a-2.0E0*3*I_ESP_H2x3y_G4z_a+2*1*I_ESP_F3y_G4z;
    abcd[iGrid*1890+301] = 4.0E0*I_ESP_K4x2yz_G4z_aa-2.0E0*2*I_ESP_H2x2yz_G4z_a-2.0E0*3*I_ESP_H2x2yz_G4z_a+2*1*I_ESP_F2yz_G4z;
    abcd[iGrid*1890+302] = 4.0E0*I_ESP_K4xy2z_G4z_aa-2.0E0*2*I_ESP_H2xy2z_G4z_a-2.0E0*3*I_ESP_H2xy2z_G4z_a+2*1*I_ESP_Fy2z_G4z;
    abcd[iGrid*1890+303] = 4.0E0*I_ESP_K4x3z_G4z_aa-2.0E0*2*I_ESP_H2x3z_G4z_a-2.0E0*3*I_ESP_H2x3z_G4z_a+2*1*I_ESP_F3z_G4z;
    abcd[iGrid*1890+304] = 4.0E0*I_ESP_K3x4y_G4z_aa-2.0E0*1*I_ESP_Hx4y_G4z_a-2.0E0*2*I_ESP_Hx4y_G4z_a;
    abcd[iGrid*1890+305] = 4.0E0*I_ESP_K3x3yz_G4z_aa-2.0E0*1*I_ESP_Hx3yz_G4z_a-2.0E0*2*I_ESP_Hx3yz_G4z_a;
    abcd[iGrid*1890+306] = 4.0E0*I_ESP_K3x2y2z_G4z_aa-2.0E0*1*I_ESP_Hx2y2z_G4z_a-2.0E0*2*I_ESP_Hx2y2z_G4z_a;
    abcd[iGrid*1890+307] = 4.0E0*I_ESP_K3xy3z_G4z_aa-2.0E0*1*I_ESP_Hxy3z_G4z_a-2.0E0*2*I_ESP_Hxy3z_G4z_a;
    abcd[iGrid*1890+308] = 4.0E0*I_ESP_K3x4z_G4z_aa-2.0E0*1*I_ESP_Hx4z_G4z_a-2.0E0*2*I_ESP_Hx4z_G4z_a;
    abcd[iGrid*1890+309] = 4.0E0*I_ESP_K2x5y_G4z_aa-2.0E0*1*I_ESP_H5y_G4z_a;
    abcd[iGrid*1890+310] = 4.0E0*I_ESP_K2x4yz_G4z_aa-2.0E0*1*I_ESP_H4yz_G4z_a;
    abcd[iGrid*1890+311] = 4.0E0*I_ESP_K2x3y2z_G4z_aa-2.0E0*1*I_ESP_H3y2z_G4z_a;
    abcd[iGrid*1890+312] = 4.0E0*I_ESP_K2x2y3z_G4z_aa-2.0E0*1*I_ESP_H2y3z_G4z_a;
    abcd[iGrid*1890+313] = 4.0E0*I_ESP_K2xy4z_G4z_aa-2.0E0*1*I_ESP_Hy4z_G4z_a;
    abcd[iGrid*1890+314] = 4.0E0*I_ESP_K2x5z_G4z_aa-2.0E0*1*I_ESP_H5z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_G_aa
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_F_G
     ************************************************************/
    abcd[iGrid*1890+315] = 4.0E0*I_ESP_K6xy_G4x_aa-2.0E0*5*I_ESP_H4xy_G4x_a;
    abcd[iGrid*1890+316] = 4.0E0*I_ESP_K5x2y_G4x_aa-2.0E0*1*I_ESP_H5x_G4x_a-2.0E0*4*I_ESP_H3x2y_G4x_a+4*1*I_ESP_F3x_G4x;
    abcd[iGrid*1890+317] = 4.0E0*I_ESP_K5xyz_G4x_aa-2.0E0*4*I_ESP_H3xyz_G4x_a;
    abcd[iGrid*1890+318] = 4.0E0*I_ESP_K4x3y_G4x_aa-2.0E0*2*I_ESP_H4xy_G4x_a-2.0E0*3*I_ESP_H2x3y_G4x_a+3*2*I_ESP_F2xy_G4x;
    abcd[iGrid*1890+319] = 4.0E0*I_ESP_K4x2yz_G4x_aa-2.0E0*1*I_ESP_H4xz_G4x_a-2.0E0*3*I_ESP_H2x2yz_G4x_a+3*1*I_ESP_F2xz_G4x;
    abcd[iGrid*1890+320] = 4.0E0*I_ESP_K4xy2z_G4x_aa-2.0E0*3*I_ESP_H2xy2z_G4x_a;
    abcd[iGrid*1890+321] = 4.0E0*I_ESP_K3x4y_G4x_aa-2.0E0*3*I_ESP_H3x2y_G4x_a-2.0E0*2*I_ESP_Hx4y_G4x_a+2*3*I_ESP_Fx2y_G4x;
    abcd[iGrid*1890+322] = 4.0E0*I_ESP_K3x3yz_G4x_aa-2.0E0*2*I_ESP_H3xyz_G4x_a-2.0E0*2*I_ESP_Hx3yz_G4x_a+2*2*I_ESP_Fxyz_G4x;
    abcd[iGrid*1890+323] = 4.0E0*I_ESP_K3x2y2z_G4x_aa-2.0E0*1*I_ESP_H3x2z_G4x_a-2.0E0*2*I_ESP_Hx2y2z_G4x_a+2*1*I_ESP_Fx2z_G4x;
    abcd[iGrid*1890+324] = 4.0E0*I_ESP_K3xy3z_G4x_aa-2.0E0*2*I_ESP_Hxy3z_G4x_a;
    abcd[iGrid*1890+325] = 4.0E0*I_ESP_K2x5y_G4x_aa-2.0E0*4*I_ESP_H2x3y_G4x_a-2.0E0*1*I_ESP_H5y_G4x_a+4*I_ESP_F3y_G4x;
    abcd[iGrid*1890+326] = 4.0E0*I_ESP_K2x4yz_G4x_aa-2.0E0*3*I_ESP_H2x2yz_G4x_a-2.0E0*1*I_ESP_H4yz_G4x_a+3*I_ESP_F2yz_G4x;
    abcd[iGrid*1890+327] = 4.0E0*I_ESP_K2x3y2z_G4x_aa-2.0E0*2*I_ESP_H2xy2z_G4x_a-2.0E0*1*I_ESP_H3y2z_G4x_a+2*I_ESP_Fy2z_G4x;
    abcd[iGrid*1890+328] = 4.0E0*I_ESP_K2x2y3z_G4x_aa-2.0E0*1*I_ESP_H2x3z_G4x_a-2.0E0*1*I_ESP_H2y3z_G4x_a+1*I_ESP_F3z_G4x;
    abcd[iGrid*1890+329] = 4.0E0*I_ESP_K2xy4z_G4x_aa-2.0E0*1*I_ESP_Hy4z_G4x_a;
    abcd[iGrid*1890+330] = 4.0E0*I_ESP_Kx6y_G4x_aa-2.0E0*5*I_ESP_Hx4y_G4x_a;
    abcd[iGrid*1890+331] = 4.0E0*I_ESP_Kx5yz_G4x_aa-2.0E0*4*I_ESP_Hx3yz_G4x_a;
    abcd[iGrid*1890+332] = 4.0E0*I_ESP_Kx4y2z_G4x_aa-2.0E0*3*I_ESP_Hx2y2z_G4x_a;
    abcd[iGrid*1890+333] = 4.0E0*I_ESP_Kx3y3z_G4x_aa-2.0E0*2*I_ESP_Hxy3z_G4x_a;
    abcd[iGrid*1890+334] = 4.0E0*I_ESP_Kx2y4z_G4x_aa-2.0E0*1*I_ESP_Hx4z_G4x_a;
    abcd[iGrid*1890+335] = 4.0E0*I_ESP_Kxy5z_G4x_aa;
    abcd[iGrid*1890+336] = 4.0E0*I_ESP_K6xy_G3xy_aa-2.0E0*5*I_ESP_H4xy_G3xy_a;
    abcd[iGrid*1890+337] = 4.0E0*I_ESP_K5x2y_G3xy_aa-2.0E0*1*I_ESP_H5x_G3xy_a-2.0E0*4*I_ESP_H3x2y_G3xy_a+4*1*I_ESP_F3x_G3xy;
    abcd[iGrid*1890+338] = 4.0E0*I_ESP_K5xyz_G3xy_aa-2.0E0*4*I_ESP_H3xyz_G3xy_a;
    abcd[iGrid*1890+339] = 4.0E0*I_ESP_K4x3y_G3xy_aa-2.0E0*2*I_ESP_H4xy_G3xy_a-2.0E0*3*I_ESP_H2x3y_G3xy_a+3*2*I_ESP_F2xy_G3xy;
    abcd[iGrid*1890+340] = 4.0E0*I_ESP_K4x2yz_G3xy_aa-2.0E0*1*I_ESP_H4xz_G3xy_a-2.0E0*3*I_ESP_H2x2yz_G3xy_a+3*1*I_ESP_F2xz_G3xy;
    abcd[iGrid*1890+341] = 4.0E0*I_ESP_K4xy2z_G3xy_aa-2.0E0*3*I_ESP_H2xy2z_G3xy_a;
    abcd[iGrid*1890+342] = 4.0E0*I_ESP_K3x4y_G3xy_aa-2.0E0*3*I_ESP_H3x2y_G3xy_a-2.0E0*2*I_ESP_Hx4y_G3xy_a+2*3*I_ESP_Fx2y_G3xy;
    abcd[iGrid*1890+343] = 4.0E0*I_ESP_K3x3yz_G3xy_aa-2.0E0*2*I_ESP_H3xyz_G3xy_a-2.0E0*2*I_ESP_Hx3yz_G3xy_a+2*2*I_ESP_Fxyz_G3xy;
    abcd[iGrid*1890+344] = 4.0E0*I_ESP_K3x2y2z_G3xy_aa-2.0E0*1*I_ESP_H3x2z_G3xy_a-2.0E0*2*I_ESP_Hx2y2z_G3xy_a+2*1*I_ESP_Fx2z_G3xy;
    abcd[iGrid*1890+345] = 4.0E0*I_ESP_K3xy3z_G3xy_aa-2.0E0*2*I_ESP_Hxy3z_G3xy_a;
    abcd[iGrid*1890+346] = 4.0E0*I_ESP_K2x5y_G3xy_aa-2.0E0*4*I_ESP_H2x3y_G3xy_a-2.0E0*1*I_ESP_H5y_G3xy_a+4*I_ESP_F3y_G3xy;
    abcd[iGrid*1890+347] = 4.0E0*I_ESP_K2x4yz_G3xy_aa-2.0E0*3*I_ESP_H2x2yz_G3xy_a-2.0E0*1*I_ESP_H4yz_G3xy_a+3*I_ESP_F2yz_G3xy;
    abcd[iGrid*1890+348] = 4.0E0*I_ESP_K2x3y2z_G3xy_aa-2.0E0*2*I_ESP_H2xy2z_G3xy_a-2.0E0*1*I_ESP_H3y2z_G3xy_a+2*I_ESP_Fy2z_G3xy;
    abcd[iGrid*1890+349] = 4.0E0*I_ESP_K2x2y3z_G3xy_aa-2.0E0*1*I_ESP_H2x3z_G3xy_a-2.0E0*1*I_ESP_H2y3z_G3xy_a+1*I_ESP_F3z_G3xy;
    abcd[iGrid*1890+350] = 4.0E0*I_ESP_K2xy4z_G3xy_aa-2.0E0*1*I_ESP_Hy4z_G3xy_a;
    abcd[iGrid*1890+351] = 4.0E0*I_ESP_Kx6y_G3xy_aa-2.0E0*5*I_ESP_Hx4y_G3xy_a;
    abcd[iGrid*1890+352] = 4.0E0*I_ESP_Kx5yz_G3xy_aa-2.0E0*4*I_ESP_Hx3yz_G3xy_a;
    abcd[iGrid*1890+353] = 4.0E0*I_ESP_Kx4y2z_G3xy_aa-2.0E0*3*I_ESP_Hx2y2z_G3xy_a;
    abcd[iGrid*1890+354] = 4.0E0*I_ESP_Kx3y3z_G3xy_aa-2.0E0*2*I_ESP_Hxy3z_G3xy_a;
    abcd[iGrid*1890+355] = 4.0E0*I_ESP_Kx2y4z_G3xy_aa-2.0E0*1*I_ESP_Hx4z_G3xy_a;
    abcd[iGrid*1890+356] = 4.0E0*I_ESP_Kxy5z_G3xy_aa;
    abcd[iGrid*1890+357] = 4.0E0*I_ESP_K6xy_G3xz_aa-2.0E0*5*I_ESP_H4xy_G3xz_a;
    abcd[iGrid*1890+358] = 4.0E0*I_ESP_K5x2y_G3xz_aa-2.0E0*1*I_ESP_H5x_G3xz_a-2.0E0*4*I_ESP_H3x2y_G3xz_a+4*1*I_ESP_F3x_G3xz;
    abcd[iGrid*1890+359] = 4.0E0*I_ESP_K5xyz_G3xz_aa-2.0E0*4*I_ESP_H3xyz_G3xz_a;
    abcd[iGrid*1890+360] = 4.0E0*I_ESP_K4x3y_G3xz_aa-2.0E0*2*I_ESP_H4xy_G3xz_a-2.0E0*3*I_ESP_H2x3y_G3xz_a+3*2*I_ESP_F2xy_G3xz;
    abcd[iGrid*1890+361] = 4.0E0*I_ESP_K4x2yz_G3xz_aa-2.0E0*1*I_ESP_H4xz_G3xz_a-2.0E0*3*I_ESP_H2x2yz_G3xz_a+3*1*I_ESP_F2xz_G3xz;
    abcd[iGrid*1890+362] = 4.0E0*I_ESP_K4xy2z_G3xz_aa-2.0E0*3*I_ESP_H2xy2z_G3xz_a;
    abcd[iGrid*1890+363] = 4.0E0*I_ESP_K3x4y_G3xz_aa-2.0E0*3*I_ESP_H3x2y_G3xz_a-2.0E0*2*I_ESP_Hx4y_G3xz_a+2*3*I_ESP_Fx2y_G3xz;
    abcd[iGrid*1890+364] = 4.0E0*I_ESP_K3x3yz_G3xz_aa-2.0E0*2*I_ESP_H3xyz_G3xz_a-2.0E0*2*I_ESP_Hx3yz_G3xz_a+2*2*I_ESP_Fxyz_G3xz;
    abcd[iGrid*1890+365] = 4.0E0*I_ESP_K3x2y2z_G3xz_aa-2.0E0*1*I_ESP_H3x2z_G3xz_a-2.0E0*2*I_ESP_Hx2y2z_G3xz_a+2*1*I_ESP_Fx2z_G3xz;
    abcd[iGrid*1890+366] = 4.0E0*I_ESP_K3xy3z_G3xz_aa-2.0E0*2*I_ESP_Hxy3z_G3xz_a;
    abcd[iGrid*1890+367] = 4.0E0*I_ESP_K2x5y_G3xz_aa-2.0E0*4*I_ESP_H2x3y_G3xz_a-2.0E0*1*I_ESP_H5y_G3xz_a+4*I_ESP_F3y_G3xz;
    abcd[iGrid*1890+368] = 4.0E0*I_ESP_K2x4yz_G3xz_aa-2.0E0*3*I_ESP_H2x2yz_G3xz_a-2.0E0*1*I_ESP_H4yz_G3xz_a+3*I_ESP_F2yz_G3xz;
    abcd[iGrid*1890+369] = 4.0E0*I_ESP_K2x3y2z_G3xz_aa-2.0E0*2*I_ESP_H2xy2z_G3xz_a-2.0E0*1*I_ESP_H3y2z_G3xz_a+2*I_ESP_Fy2z_G3xz;
    abcd[iGrid*1890+370] = 4.0E0*I_ESP_K2x2y3z_G3xz_aa-2.0E0*1*I_ESP_H2x3z_G3xz_a-2.0E0*1*I_ESP_H2y3z_G3xz_a+1*I_ESP_F3z_G3xz;
    abcd[iGrid*1890+371] = 4.0E0*I_ESP_K2xy4z_G3xz_aa-2.0E0*1*I_ESP_Hy4z_G3xz_a;
    abcd[iGrid*1890+372] = 4.0E0*I_ESP_Kx6y_G3xz_aa-2.0E0*5*I_ESP_Hx4y_G3xz_a;
    abcd[iGrid*1890+373] = 4.0E0*I_ESP_Kx5yz_G3xz_aa-2.0E0*4*I_ESP_Hx3yz_G3xz_a;
    abcd[iGrid*1890+374] = 4.0E0*I_ESP_Kx4y2z_G3xz_aa-2.0E0*3*I_ESP_Hx2y2z_G3xz_a;
    abcd[iGrid*1890+375] = 4.0E0*I_ESP_Kx3y3z_G3xz_aa-2.0E0*2*I_ESP_Hxy3z_G3xz_a;
    abcd[iGrid*1890+376] = 4.0E0*I_ESP_Kx2y4z_G3xz_aa-2.0E0*1*I_ESP_Hx4z_G3xz_a;
    abcd[iGrid*1890+377] = 4.0E0*I_ESP_Kxy5z_G3xz_aa;
    abcd[iGrid*1890+378] = 4.0E0*I_ESP_K6xy_G2x2y_aa-2.0E0*5*I_ESP_H4xy_G2x2y_a;
    abcd[iGrid*1890+379] = 4.0E0*I_ESP_K5x2y_G2x2y_aa-2.0E0*1*I_ESP_H5x_G2x2y_a-2.0E0*4*I_ESP_H3x2y_G2x2y_a+4*1*I_ESP_F3x_G2x2y;
    abcd[iGrid*1890+380] = 4.0E0*I_ESP_K5xyz_G2x2y_aa-2.0E0*4*I_ESP_H3xyz_G2x2y_a;
    abcd[iGrid*1890+381] = 4.0E0*I_ESP_K4x3y_G2x2y_aa-2.0E0*2*I_ESP_H4xy_G2x2y_a-2.0E0*3*I_ESP_H2x3y_G2x2y_a+3*2*I_ESP_F2xy_G2x2y;
    abcd[iGrid*1890+382] = 4.0E0*I_ESP_K4x2yz_G2x2y_aa-2.0E0*1*I_ESP_H4xz_G2x2y_a-2.0E0*3*I_ESP_H2x2yz_G2x2y_a+3*1*I_ESP_F2xz_G2x2y;
    abcd[iGrid*1890+383] = 4.0E0*I_ESP_K4xy2z_G2x2y_aa-2.0E0*3*I_ESP_H2xy2z_G2x2y_a;
    abcd[iGrid*1890+384] = 4.0E0*I_ESP_K3x4y_G2x2y_aa-2.0E0*3*I_ESP_H3x2y_G2x2y_a-2.0E0*2*I_ESP_Hx4y_G2x2y_a+2*3*I_ESP_Fx2y_G2x2y;
    abcd[iGrid*1890+385] = 4.0E0*I_ESP_K3x3yz_G2x2y_aa-2.0E0*2*I_ESP_H3xyz_G2x2y_a-2.0E0*2*I_ESP_Hx3yz_G2x2y_a+2*2*I_ESP_Fxyz_G2x2y;
    abcd[iGrid*1890+386] = 4.0E0*I_ESP_K3x2y2z_G2x2y_aa-2.0E0*1*I_ESP_H3x2z_G2x2y_a-2.0E0*2*I_ESP_Hx2y2z_G2x2y_a+2*1*I_ESP_Fx2z_G2x2y;
    abcd[iGrid*1890+387] = 4.0E0*I_ESP_K3xy3z_G2x2y_aa-2.0E0*2*I_ESP_Hxy3z_G2x2y_a;
    abcd[iGrid*1890+388] = 4.0E0*I_ESP_K2x5y_G2x2y_aa-2.0E0*4*I_ESP_H2x3y_G2x2y_a-2.0E0*1*I_ESP_H5y_G2x2y_a+4*I_ESP_F3y_G2x2y;
    abcd[iGrid*1890+389] = 4.0E0*I_ESP_K2x4yz_G2x2y_aa-2.0E0*3*I_ESP_H2x2yz_G2x2y_a-2.0E0*1*I_ESP_H4yz_G2x2y_a+3*I_ESP_F2yz_G2x2y;
    abcd[iGrid*1890+390] = 4.0E0*I_ESP_K2x3y2z_G2x2y_aa-2.0E0*2*I_ESP_H2xy2z_G2x2y_a-2.0E0*1*I_ESP_H3y2z_G2x2y_a+2*I_ESP_Fy2z_G2x2y;
    abcd[iGrid*1890+391] = 4.0E0*I_ESP_K2x2y3z_G2x2y_aa-2.0E0*1*I_ESP_H2x3z_G2x2y_a-2.0E0*1*I_ESP_H2y3z_G2x2y_a+1*I_ESP_F3z_G2x2y;
    abcd[iGrid*1890+392] = 4.0E0*I_ESP_K2xy4z_G2x2y_aa-2.0E0*1*I_ESP_Hy4z_G2x2y_a;
    abcd[iGrid*1890+393] = 4.0E0*I_ESP_Kx6y_G2x2y_aa-2.0E0*5*I_ESP_Hx4y_G2x2y_a;
    abcd[iGrid*1890+394] = 4.0E0*I_ESP_Kx5yz_G2x2y_aa-2.0E0*4*I_ESP_Hx3yz_G2x2y_a;
    abcd[iGrid*1890+395] = 4.0E0*I_ESP_Kx4y2z_G2x2y_aa-2.0E0*3*I_ESP_Hx2y2z_G2x2y_a;
    abcd[iGrid*1890+396] = 4.0E0*I_ESP_Kx3y3z_G2x2y_aa-2.0E0*2*I_ESP_Hxy3z_G2x2y_a;
    abcd[iGrid*1890+397] = 4.0E0*I_ESP_Kx2y4z_G2x2y_aa-2.0E0*1*I_ESP_Hx4z_G2x2y_a;
    abcd[iGrid*1890+398] = 4.0E0*I_ESP_Kxy5z_G2x2y_aa;
    abcd[iGrid*1890+399] = 4.0E0*I_ESP_K6xy_G2xyz_aa-2.0E0*5*I_ESP_H4xy_G2xyz_a;
    abcd[iGrid*1890+400] = 4.0E0*I_ESP_K5x2y_G2xyz_aa-2.0E0*1*I_ESP_H5x_G2xyz_a-2.0E0*4*I_ESP_H3x2y_G2xyz_a+4*1*I_ESP_F3x_G2xyz;
    abcd[iGrid*1890+401] = 4.0E0*I_ESP_K5xyz_G2xyz_aa-2.0E0*4*I_ESP_H3xyz_G2xyz_a;
    abcd[iGrid*1890+402] = 4.0E0*I_ESP_K4x3y_G2xyz_aa-2.0E0*2*I_ESP_H4xy_G2xyz_a-2.0E0*3*I_ESP_H2x3y_G2xyz_a+3*2*I_ESP_F2xy_G2xyz;
    abcd[iGrid*1890+403] = 4.0E0*I_ESP_K4x2yz_G2xyz_aa-2.0E0*1*I_ESP_H4xz_G2xyz_a-2.0E0*3*I_ESP_H2x2yz_G2xyz_a+3*1*I_ESP_F2xz_G2xyz;
    abcd[iGrid*1890+404] = 4.0E0*I_ESP_K4xy2z_G2xyz_aa-2.0E0*3*I_ESP_H2xy2z_G2xyz_a;
    abcd[iGrid*1890+405] = 4.0E0*I_ESP_K3x4y_G2xyz_aa-2.0E0*3*I_ESP_H3x2y_G2xyz_a-2.0E0*2*I_ESP_Hx4y_G2xyz_a+2*3*I_ESP_Fx2y_G2xyz;
    abcd[iGrid*1890+406] = 4.0E0*I_ESP_K3x3yz_G2xyz_aa-2.0E0*2*I_ESP_H3xyz_G2xyz_a-2.0E0*2*I_ESP_Hx3yz_G2xyz_a+2*2*I_ESP_Fxyz_G2xyz;
    abcd[iGrid*1890+407] = 4.0E0*I_ESP_K3x2y2z_G2xyz_aa-2.0E0*1*I_ESP_H3x2z_G2xyz_a-2.0E0*2*I_ESP_Hx2y2z_G2xyz_a+2*1*I_ESP_Fx2z_G2xyz;
    abcd[iGrid*1890+408] = 4.0E0*I_ESP_K3xy3z_G2xyz_aa-2.0E0*2*I_ESP_Hxy3z_G2xyz_a;
    abcd[iGrid*1890+409] = 4.0E0*I_ESP_K2x5y_G2xyz_aa-2.0E0*4*I_ESP_H2x3y_G2xyz_a-2.0E0*1*I_ESP_H5y_G2xyz_a+4*I_ESP_F3y_G2xyz;
    abcd[iGrid*1890+410] = 4.0E0*I_ESP_K2x4yz_G2xyz_aa-2.0E0*3*I_ESP_H2x2yz_G2xyz_a-2.0E0*1*I_ESP_H4yz_G2xyz_a+3*I_ESP_F2yz_G2xyz;
    abcd[iGrid*1890+411] = 4.0E0*I_ESP_K2x3y2z_G2xyz_aa-2.0E0*2*I_ESP_H2xy2z_G2xyz_a-2.0E0*1*I_ESP_H3y2z_G2xyz_a+2*I_ESP_Fy2z_G2xyz;
    abcd[iGrid*1890+412] = 4.0E0*I_ESP_K2x2y3z_G2xyz_aa-2.0E0*1*I_ESP_H2x3z_G2xyz_a-2.0E0*1*I_ESP_H2y3z_G2xyz_a+1*I_ESP_F3z_G2xyz;
    abcd[iGrid*1890+413] = 4.0E0*I_ESP_K2xy4z_G2xyz_aa-2.0E0*1*I_ESP_Hy4z_G2xyz_a;
    abcd[iGrid*1890+414] = 4.0E0*I_ESP_Kx6y_G2xyz_aa-2.0E0*5*I_ESP_Hx4y_G2xyz_a;
    abcd[iGrid*1890+415] = 4.0E0*I_ESP_Kx5yz_G2xyz_aa-2.0E0*4*I_ESP_Hx3yz_G2xyz_a;
    abcd[iGrid*1890+416] = 4.0E0*I_ESP_Kx4y2z_G2xyz_aa-2.0E0*3*I_ESP_Hx2y2z_G2xyz_a;
    abcd[iGrid*1890+417] = 4.0E0*I_ESP_Kx3y3z_G2xyz_aa-2.0E0*2*I_ESP_Hxy3z_G2xyz_a;
    abcd[iGrid*1890+418] = 4.0E0*I_ESP_Kx2y4z_G2xyz_aa-2.0E0*1*I_ESP_Hx4z_G2xyz_a;
    abcd[iGrid*1890+419] = 4.0E0*I_ESP_Kxy5z_G2xyz_aa;
    abcd[iGrid*1890+420] = 4.0E0*I_ESP_K6xy_G2x2z_aa-2.0E0*5*I_ESP_H4xy_G2x2z_a;
    abcd[iGrid*1890+421] = 4.0E0*I_ESP_K5x2y_G2x2z_aa-2.0E0*1*I_ESP_H5x_G2x2z_a-2.0E0*4*I_ESP_H3x2y_G2x2z_a+4*1*I_ESP_F3x_G2x2z;
    abcd[iGrid*1890+422] = 4.0E0*I_ESP_K5xyz_G2x2z_aa-2.0E0*4*I_ESP_H3xyz_G2x2z_a;
    abcd[iGrid*1890+423] = 4.0E0*I_ESP_K4x3y_G2x2z_aa-2.0E0*2*I_ESP_H4xy_G2x2z_a-2.0E0*3*I_ESP_H2x3y_G2x2z_a+3*2*I_ESP_F2xy_G2x2z;
    abcd[iGrid*1890+424] = 4.0E0*I_ESP_K4x2yz_G2x2z_aa-2.0E0*1*I_ESP_H4xz_G2x2z_a-2.0E0*3*I_ESP_H2x2yz_G2x2z_a+3*1*I_ESP_F2xz_G2x2z;
    abcd[iGrid*1890+425] = 4.0E0*I_ESP_K4xy2z_G2x2z_aa-2.0E0*3*I_ESP_H2xy2z_G2x2z_a;
    abcd[iGrid*1890+426] = 4.0E0*I_ESP_K3x4y_G2x2z_aa-2.0E0*3*I_ESP_H3x2y_G2x2z_a-2.0E0*2*I_ESP_Hx4y_G2x2z_a+2*3*I_ESP_Fx2y_G2x2z;
    abcd[iGrid*1890+427] = 4.0E0*I_ESP_K3x3yz_G2x2z_aa-2.0E0*2*I_ESP_H3xyz_G2x2z_a-2.0E0*2*I_ESP_Hx3yz_G2x2z_a+2*2*I_ESP_Fxyz_G2x2z;
    abcd[iGrid*1890+428] = 4.0E0*I_ESP_K3x2y2z_G2x2z_aa-2.0E0*1*I_ESP_H3x2z_G2x2z_a-2.0E0*2*I_ESP_Hx2y2z_G2x2z_a+2*1*I_ESP_Fx2z_G2x2z;
    abcd[iGrid*1890+429] = 4.0E0*I_ESP_K3xy3z_G2x2z_aa-2.0E0*2*I_ESP_Hxy3z_G2x2z_a;
    abcd[iGrid*1890+430] = 4.0E0*I_ESP_K2x5y_G2x2z_aa-2.0E0*4*I_ESP_H2x3y_G2x2z_a-2.0E0*1*I_ESP_H5y_G2x2z_a+4*I_ESP_F3y_G2x2z;
    abcd[iGrid*1890+431] = 4.0E0*I_ESP_K2x4yz_G2x2z_aa-2.0E0*3*I_ESP_H2x2yz_G2x2z_a-2.0E0*1*I_ESP_H4yz_G2x2z_a+3*I_ESP_F2yz_G2x2z;
    abcd[iGrid*1890+432] = 4.0E0*I_ESP_K2x3y2z_G2x2z_aa-2.0E0*2*I_ESP_H2xy2z_G2x2z_a-2.0E0*1*I_ESP_H3y2z_G2x2z_a+2*I_ESP_Fy2z_G2x2z;
    abcd[iGrid*1890+433] = 4.0E0*I_ESP_K2x2y3z_G2x2z_aa-2.0E0*1*I_ESP_H2x3z_G2x2z_a-2.0E0*1*I_ESP_H2y3z_G2x2z_a+1*I_ESP_F3z_G2x2z;
    abcd[iGrid*1890+434] = 4.0E0*I_ESP_K2xy4z_G2x2z_aa-2.0E0*1*I_ESP_Hy4z_G2x2z_a;
    abcd[iGrid*1890+435] = 4.0E0*I_ESP_Kx6y_G2x2z_aa-2.0E0*5*I_ESP_Hx4y_G2x2z_a;
    abcd[iGrid*1890+436] = 4.0E0*I_ESP_Kx5yz_G2x2z_aa-2.0E0*4*I_ESP_Hx3yz_G2x2z_a;
    abcd[iGrid*1890+437] = 4.0E0*I_ESP_Kx4y2z_G2x2z_aa-2.0E0*3*I_ESP_Hx2y2z_G2x2z_a;
    abcd[iGrid*1890+438] = 4.0E0*I_ESP_Kx3y3z_G2x2z_aa-2.0E0*2*I_ESP_Hxy3z_G2x2z_a;
    abcd[iGrid*1890+439] = 4.0E0*I_ESP_Kx2y4z_G2x2z_aa-2.0E0*1*I_ESP_Hx4z_G2x2z_a;
    abcd[iGrid*1890+440] = 4.0E0*I_ESP_Kxy5z_G2x2z_aa;
    abcd[iGrid*1890+441] = 4.0E0*I_ESP_K6xy_Gx3y_aa-2.0E0*5*I_ESP_H4xy_Gx3y_a;
    abcd[iGrid*1890+442] = 4.0E0*I_ESP_K5x2y_Gx3y_aa-2.0E0*1*I_ESP_H5x_Gx3y_a-2.0E0*4*I_ESP_H3x2y_Gx3y_a+4*1*I_ESP_F3x_Gx3y;
    abcd[iGrid*1890+443] = 4.0E0*I_ESP_K5xyz_Gx3y_aa-2.0E0*4*I_ESP_H3xyz_Gx3y_a;
    abcd[iGrid*1890+444] = 4.0E0*I_ESP_K4x3y_Gx3y_aa-2.0E0*2*I_ESP_H4xy_Gx3y_a-2.0E0*3*I_ESP_H2x3y_Gx3y_a+3*2*I_ESP_F2xy_Gx3y;
    abcd[iGrid*1890+445] = 4.0E0*I_ESP_K4x2yz_Gx3y_aa-2.0E0*1*I_ESP_H4xz_Gx3y_a-2.0E0*3*I_ESP_H2x2yz_Gx3y_a+3*1*I_ESP_F2xz_Gx3y;
    abcd[iGrid*1890+446] = 4.0E0*I_ESP_K4xy2z_Gx3y_aa-2.0E0*3*I_ESP_H2xy2z_Gx3y_a;
    abcd[iGrid*1890+447] = 4.0E0*I_ESP_K3x4y_Gx3y_aa-2.0E0*3*I_ESP_H3x2y_Gx3y_a-2.0E0*2*I_ESP_Hx4y_Gx3y_a+2*3*I_ESP_Fx2y_Gx3y;
    abcd[iGrid*1890+448] = 4.0E0*I_ESP_K3x3yz_Gx3y_aa-2.0E0*2*I_ESP_H3xyz_Gx3y_a-2.0E0*2*I_ESP_Hx3yz_Gx3y_a+2*2*I_ESP_Fxyz_Gx3y;
    abcd[iGrid*1890+449] = 4.0E0*I_ESP_K3x2y2z_Gx3y_aa-2.0E0*1*I_ESP_H3x2z_Gx3y_a-2.0E0*2*I_ESP_Hx2y2z_Gx3y_a+2*1*I_ESP_Fx2z_Gx3y;
    abcd[iGrid*1890+450] = 4.0E0*I_ESP_K3xy3z_Gx3y_aa-2.0E0*2*I_ESP_Hxy3z_Gx3y_a;
    abcd[iGrid*1890+451] = 4.0E0*I_ESP_K2x5y_Gx3y_aa-2.0E0*4*I_ESP_H2x3y_Gx3y_a-2.0E0*1*I_ESP_H5y_Gx3y_a+4*I_ESP_F3y_Gx3y;
    abcd[iGrid*1890+452] = 4.0E0*I_ESP_K2x4yz_Gx3y_aa-2.0E0*3*I_ESP_H2x2yz_Gx3y_a-2.0E0*1*I_ESP_H4yz_Gx3y_a+3*I_ESP_F2yz_Gx3y;
    abcd[iGrid*1890+453] = 4.0E0*I_ESP_K2x3y2z_Gx3y_aa-2.0E0*2*I_ESP_H2xy2z_Gx3y_a-2.0E0*1*I_ESP_H3y2z_Gx3y_a+2*I_ESP_Fy2z_Gx3y;
    abcd[iGrid*1890+454] = 4.0E0*I_ESP_K2x2y3z_Gx3y_aa-2.0E0*1*I_ESP_H2x3z_Gx3y_a-2.0E0*1*I_ESP_H2y3z_Gx3y_a+1*I_ESP_F3z_Gx3y;
    abcd[iGrid*1890+455] = 4.0E0*I_ESP_K2xy4z_Gx3y_aa-2.0E0*1*I_ESP_Hy4z_Gx3y_a;
    abcd[iGrid*1890+456] = 4.0E0*I_ESP_Kx6y_Gx3y_aa-2.0E0*5*I_ESP_Hx4y_Gx3y_a;
    abcd[iGrid*1890+457] = 4.0E0*I_ESP_Kx5yz_Gx3y_aa-2.0E0*4*I_ESP_Hx3yz_Gx3y_a;
    abcd[iGrid*1890+458] = 4.0E0*I_ESP_Kx4y2z_Gx3y_aa-2.0E0*3*I_ESP_Hx2y2z_Gx3y_a;
    abcd[iGrid*1890+459] = 4.0E0*I_ESP_Kx3y3z_Gx3y_aa-2.0E0*2*I_ESP_Hxy3z_Gx3y_a;
    abcd[iGrid*1890+460] = 4.0E0*I_ESP_Kx2y4z_Gx3y_aa-2.0E0*1*I_ESP_Hx4z_Gx3y_a;
    abcd[iGrid*1890+461] = 4.0E0*I_ESP_Kxy5z_Gx3y_aa;
    abcd[iGrid*1890+462] = 4.0E0*I_ESP_K6xy_Gx2yz_aa-2.0E0*5*I_ESP_H4xy_Gx2yz_a;
    abcd[iGrid*1890+463] = 4.0E0*I_ESP_K5x2y_Gx2yz_aa-2.0E0*1*I_ESP_H5x_Gx2yz_a-2.0E0*4*I_ESP_H3x2y_Gx2yz_a+4*1*I_ESP_F3x_Gx2yz;
    abcd[iGrid*1890+464] = 4.0E0*I_ESP_K5xyz_Gx2yz_aa-2.0E0*4*I_ESP_H3xyz_Gx2yz_a;
    abcd[iGrid*1890+465] = 4.0E0*I_ESP_K4x3y_Gx2yz_aa-2.0E0*2*I_ESP_H4xy_Gx2yz_a-2.0E0*3*I_ESP_H2x3y_Gx2yz_a+3*2*I_ESP_F2xy_Gx2yz;
    abcd[iGrid*1890+466] = 4.0E0*I_ESP_K4x2yz_Gx2yz_aa-2.0E0*1*I_ESP_H4xz_Gx2yz_a-2.0E0*3*I_ESP_H2x2yz_Gx2yz_a+3*1*I_ESP_F2xz_Gx2yz;
    abcd[iGrid*1890+467] = 4.0E0*I_ESP_K4xy2z_Gx2yz_aa-2.0E0*3*I_ESP_H2xy2z_Gx2yz_a;
    abcd[iGrid*1890+468] = 4.0E0*I_ESP_K3x4y_Gx2yz_aa-2.0E0*3*I_ESP_H3x2y_Gx2yz_a-2.0E0*2*I_ESP_Hx4y_Gx2yz_a+2*3*I_ESP_Fx2y_Gx2yz;
    abcd[iGrid*1890+469] = 4.0E0*I_ESP_K3x3yz_Gx2yz_aa-2.0E0*2*I_ESP_H3xyz_Gx2yz_a-2.0E0*2*I_ESP_Hx3yz_Gx2yz_a+2*2*I_ESP_Fxyz_Gx2yz;
    abcd[iGrid*1890+470] = 4.0E0*I_ESP_K3x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_H3x2z_Gx2yz_a-2.0E0*2*I_ESP_Hx2y2z_Gx2yz_a+2*1*I_ESP_Fx2z_Gx2yz;
    abcd[iGrid*1890+471] = 4.0E0*I_ESP_K3xy3z_Gx2yz_aa-2.0E0*2*I_ESP_Hxy3z_Gx2yz_a;
    abcd[iGrid*1890+472] = 4.0E0*I_ESP_K2x5y_Gx2yz_aa-2.0E0*4*I_ESP_H2x3y_Gx2yz_a-2.0E0*1*I_ESP_H5y_Gx2yz_a+4*I_ESP_F3y_Gx2yz;
    abcd[iGrid*1890+473] = 4.0E0*I_ESP_K2x4yz_Gx2yz_aa-2.0E0*3*I_ESP_H2x2yz_Gx2yz_a-2.0E0*1*I_ESP_H4yz_Gx2yz_a+3*I_ESP_F2yz_Gx2yz;
    abcd[iGrid*1890+474] = 4.0E0*I_ESP_K2x3y2z_Gx2yz_aa-2.0E0*2*I_ESP_H2xy2z_Gx2yz_a-2.0E0*1*I_ESP_H3y2z_Gx2yz_a+2*I_ESP_Fy2z_Gx2yz;
    abcd[iGrid*1890+475] = 4.0E0*I_ESP_K2x2y3z_Gx2yz_aa-2.0E0*1*I_ESP_H2x3z_Gx2yz_a-2.0E0*1*I_ESP_H2y3z_Gx2yz_a+1*I_ESP_F3z_Gx2yz;
    abcd[iGrid*1890+476] = 4.0E0*I_ESP_K2xy4z_Gx2yz_aa-2.0E0*1*I_ESP_Hy4z_Gx2yz_a;
    abcd[iGrid*1890+477] = 4.0E0*I_ESP_Kx6y_Gx2yz_aa-2.0E0*5*I_ESP_Hx4y_Gx2yz_a;
    abcd[iGrid*1890+478] = 4.0E0*I_ESP_Kx5yz_Gx2yz_aa-2.0E0*4*I_ESP_Hx3yz_Gx2yz_a;
    abcd[iGrid*1890+479] = 4.0E0*I_ESP_Kx4y2z_Gx2yz_aa-2.0E0*3*I_ESP_Hx2y2z_Gx2yz_a;
    abcd[iGrid*1890+480] = 4.0E0*I_ESP_Kx3y3z_Gx2yz_aa-2.0E0*2*I_ESP_Hxy3z_Gx2yz_a;
    abcd[iGrid*1890+481] = 4.0E0*I_ESP_Kx2y4z_Gx2yz_aa-2.0E0*1*I_ESP_Hx4z_Gx2yz_a;
    abcd[iGrid*1890+482] = 4.0E0*I_ESP_Kxy5z_Gx2yz_aa;
    abcd[iGrid*1890+483] = 4.0E0*I_ESP_K6xy_Gxy2z_aa-2.0E0*5*I_ESP_H4xy_Gxy2z_a;
    abcd[iGrid*1890+484] = 4.0E0*I_ESP_K5x2y_Gxy2z_aa-2.0E0*1*I_ESP_H5x_Gxy2z_a-2.0E0*4*I_ESP_H3x2y_Gxy2z_a+4*1*I_ESP_F3x_Gxy2z;
    abcd[iGrid*1890+485] = 4.0E0*I_ESP_K5xyz_Gxy2z_aa-2.0E0*4*I_ESP_H3xyz_Gxy2z_a;
    abcd[iGrid*1890+486] = 4.0E0*I_ESP_K4x3y_Gxy2z_aa-2.0E0*2*I_ESP_H4xy_Gxy2z_a-2.0E0*3*I_ESP_H2x3y_Gxy2z_a+3*2*I_ESP_F2xy_Gxy2z;
    abcd[iGrid*1890+487] = 4.0E0*I_ESP_K4x2yz_Gxy2z_aa-2.0E0*1*I_ESP_H4xz_Gxy2z_a-2.0E0*3*I_ESP_H2x2yz_Gxy2z_a+3*1*I_ESP_F2xz_Gxy2z;
    abcd[iGrid*1890+488] = 4.0E0*I_ESP_K4xy2z_Gxy2z_aa-2.0E0*3*I_ESP_H2xy2z_Gxy2z_a;
    abcd[iGrid*1890+489] = 4.0E0*I_ESP_K3x4y_Gxy2z_aa-2.0E0*3*I_ESP_H3x2y_Gxy2z_a-2.0E0*2*I_ESP_Hx4y_Gxy2z_a+2*3*I_ESP_Fx2y_Gxy2z;
    abcd[iGrid*1890+490] = 4.0E0*I_ESP_K3x3yz_Gxy2z_aa-2.0E0*2*I_ESP_H3xyz_Gxy2z_a-2.0E0*2*I_ESP_Hx3yz_Gxy2z_a+2*2*I_ESP_Fxyz_Gxy2z;
    abcd[iGrid*1890+491] = 4.0E0*I_ESP_K3x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_H3x2z_Gxy2z_a-2.0E0*2*I_ESP_Hx2y2z_Gxy2z_a+2*1*I_ESP_Fx2z_Gxy2z;
    abcd[iGrid*1890+492] = 4.0E0*I_ESP_K3xy3z_Gxy2z_aa-2.0E0*2*I_ESP_Hxy3z_Gxy2z_a;
    abcd[iGrid*1890+493] = 4.0E0*I_ESP_K2x5y_Gxy2z_aa-2.0E0*4*I_ESP_H2x3y_Gxy2z_a-2.0E0*1*I_ESP_H5y_Gxy2z_a+4*I_ESP_F3y_Gxy2z;
    abcd[iGrid*1890+494] = 4.0E0*I_ESP_K2x4yz_Gxy2z_aa-2.0E0*3*I_ESP_H2x2yz_Gxy2z_a-2.0E0*1*I_ESP_H4yz_Gxy2z_a+3*I_ESP_F2yz_Gxy2z;
    abcd[iGrid*1890+495] = 4.0E0*I_ESP_K2x3y2z_Gxy2z_aa-2.0E0*2*I_ESP_H2xy2z_Gxy2z_a-2.0E0*1*I_ESP_H3y2z_Gxy2z_a+2*I_ESP_Fy2z_Gxy2z;
    abcd[iGrid*1890+496] = 4.0E0*I_ESP_K2x2y3z_Gxy2z_aa-2.0E0*1*I_ESP_H2x3z_Gxy2z_a-2.0E0*1*I_ESP_H2y3z_Gxy2z_a+1*I_ESP_F3z_Gxy2z;
    abcd[iGrid*1890+497] = 4.0E0*I_ESP_K2xy4z_Gxy2z_aa-2.0E0*1*I_ESP_Hy4z_Gxy2z_a;
    abcd[iGrid*1890+498] = 4.0E0*I_ESP_Kx6y_Gxy2z_aa-2.0E0*5*I_ESP_Hx4y_Gxy2z_a;
    abcd[iGrid*1890+499] = 4.0E0*I_ESP_Kx5yz_Gxy2z_aa-2.0E0*4*I_ESP_Hx3yz_Gxy2z_a;
    abcd[iGrid*1890+500] = 4.0E0*I_ESP_Kx4y2z_Gxy2z_aa-2.0E0*3*I_ESP_Hx2y2z_Gxy2z_a;
    abcd[iGrid*1890+501] = 4.0E0*I_ESP_Kx3y3z_Gxy2z_aa-2.0E0*2*I_ESP_Hxy3z_Gxy2z_a;
    abcd[iGrid*1890+502] = 4.0E0*I_ESP_Kx2y4z_Gxy2z_aa-2.0E0*1*I_ESP_Hx4z_Gxy2z_a;
    abcd[iGrid*1890+503] = 4.0E0*I_ESP_Kxy5z_Gxy2z_aa;
    abcd[iGrid*1890+504] = 4.0E0*I_ESP_K6xy_Gx3z_aa-2.0E0*5*I_ESP_H4xy_Gx3z_a;
    abcd[iGrid*1890+505] = 4.0E0*I_ESP_K5x2y_Gx3z_aa-2.0E0*1*I_ESP_H5x_Gx3z_a-2.0E0*4*I_ESP_H3x2y_Gx3z_a+4*1*I_ESP_F3x_Gx3z;
    abcd[iGrid*1890+506] = 4.0E0*I_ESP_K5xyz_Gx3z_aa-2.0E0*4*I_ESP_H3xyz_Gx3z_a;
    abcd[iGrid*1890+507] = 4.0E0*I_ESP_K4x3y_Gx3z_aa-2.0E0*2*I_ESP_H4xy_Gx3z_a-2.0E0*3*I_ESP_H2x3y_Gx3z_a+3*2*I_ESP_F2xy_Gx3z;
    abcd[iGrid*1890+508] = 4.0E0*I_ESP_K4x2yz_Gx3z_aa-2.0E0*1*I_ESP_H4xz_Gx3z_a-2.0E0*3*I_ESP_H2x2yz_Gx3z_a+3*1*I_ESP_F2xz_Gx3z;
    abcd[iGrid*1890+509] = 4.0E0*I_ESP_K4xy2z_Gx3z_aa-2.0E0*3*I_ESP_H2xy2z_Gx3z_a;
    abcd[iGrid*1890+510] = 4.0E0*I_ESP_K3x4y_Gx3z_aa-2.0E0*3*I_ESP_H3x2y_Gx3z_a-2.0E0*2*I_ESP_Hx4y_Gx3z_a+2*3*I_ESP_Fx2y_Gx3z;
    abcd[iGrid*1890+511] = 4.0E0*I_ESP_K3x3yz_Gx3z_aa-2.0E0*2*I_ESP_H3xyz_Gx3z_a-2.0E0*2*I_ESP_Hx3yz_Gx3z_a+2*2*I_ESP_Fxyz_Gx3z;
    abcd[iGrid*1890+512] = 4.0E0*I_ESP_K3x2y2z_Gx3z_aa-2.0E0*1*I_ESP_H3x2z_Gx3z_a-2.0E0*2*I_ESP_Hx2y2z_Gx3z_a+2*1*I_ESP_Fx2z_Gx3z;
    abcd[iGrid*1890+513] = 4.0E0*I_ESP_K3xy3z_Gx3z_aa-2.0E0*2*I_ESP_Hxy3z_Gx3z_a;
    abcd[iGrid*1890+514] = 4.0E0*I_ESP_K2x5y_Gx3z_aa-2.0E0*4*I_ESP_H2x3y_Gx3z_a-2.0E0*1*I_ESP_H5y_Gx3z_a+4*I_ESP_F3y_Gx3z;
    abcd[iGrid*1890+515] = 4.0E0*I_ESP_K2x4yz_Gx3z_aa-2.0E0*3*I_ESP_H2x2yz_Gx3z_a-2.0E0*1*I_ESP_H4yz_Gx3z_a+3*I_ESP_F2yz_Gx3z;
    abcd[iGrid*1890+516] = 4.0E0*I_ESP_K2x3y2z_Gx3z_aa-2.0E0*2*I_ESP_H2xy2z_Gx3z_a-2.0E0*1*I_ESP_H3y2z_Gx3z_a+2*I_ESP_Fy2z_Gx3z;
    abcd[iGrid*1890+517] = 4.0E0*I_ESP_K2x2y3z_Gx3z_aa-2.0E0*1*I_ESP_H2x3z_Gx3z_a-2.0E0*1*I_ESP_H2y3z_Gx3z_a+1*I_ESP_F3z_Gx3z;
    abcd[iGrid*1890+518] = 4.0E0*I_ESP_K2xy4z_Gx3z_aa-2.0E0*1*I_ESP_Hy4z_Gx3z_a;
    abcd[iGrid*1890+519] = 4.0E0*I_ESP_Kx6y_Gx3z_aa-2.0E0*5*I_ESP_Hx4y_Gx3z_a;
    abcd[iGrid*1890+520] = 4.0E0*I_ESP_Kx5yz_Gx3z_aa-2.0E0*4*I_ESP_Hx3yz_Gx3z_a;
    abcd[iGrid*1890+521] = 4.0E0*I_ESP_Kx4y2z_Gx3z_aa-2.0E0*3*I_ESP_Hx2y2z_Gx3z_a;
    abcd[iGrid*1890+522] = 4.0E0*I_ESP_Kx3y3z_Gx3z_aa-2.0E0*2*I_ESP_Hxy3z_Gx3z_a;
    abcd[iGrid*1890+523] = 4.0E0*I_ESP_Kx2y4z_Gx3z_aa-2.0E0*1*I_ESP_Hx4z_Gx3z_a;
    abcd[iGrid*1890+524] = 4.0E0*I_ESP_Kxy5z_Gx3z_aa;
    abcd[iGrid*1890+525] = 4.0E0*I_ESP_K6xy_G4y_aa-2.0E0*5*I_ESP_H4xy_G4y_a;
    abcd[iGrid*1890+526] = 4.0E0*I_ESP_K5x2y_G4y_aa-2.0E0*1*I_ESP_H5x_G4y_a-2.0E0*4*I_ESP_H3x2y_G4y_a+4*1*I_ESP_F3x_G4y;
    abcd[iGrid*1890+527] = 4.0E0*I_ESP_K5xyz_G4y_aa-2.0E0*4*I_ESP_H3xyz_G4y_a;
    abcd[iGrid*1890+528] = 4.0E0*I_ESP_K4x3y_G4y_aa-2.0E0*2*I_ESP_H4xy_G4y_a-2.0E0*3*I_ESP_H2x3y_G4y_a+3*2*I_ESP_F2xy_G4y;
    abcd[iGrid*1890+529] = 4.0E0*I_ESP_K4x2yz_G4y_aa-2.0E0*1*I_ESP_H4xz_G4y_a-2.0E0*3*I_ESP_H2x2yz_G4y_a+3*1*I_ESP_F2xz_G4y;
    abcd[iGrid*1890+530] = 4.0E0*I_ESP_K4xy2z_G4y_aa-2.0E0*3*I_ESP_H2xy2z_G4y_a;
    abcd[iGrid*1890+531] = 4.0E0*I_ESP_K3x4y_G4y_aa-2.0E0*3*I_ESP_H3x2y_G4y_a-2.0E0*2*I_ESP_Hx4y_G4y_a+2*3*I_ESP_Fx2y_G4y;
    abcd[iGrid*1890+532] = 4.0E0*I_ESP_K3x3yz_G4y_aa-2.0E0*2*I_ESP_H3xyz_G4y_a-2.0E0*2*I_ESP_Hx3yz_G4y_a+2*2*I_ESP_Fxyz_G4y;
    abcd[iGrid*1890+533] = 4.0E0*I_ESP_K3x2y2z_G4y_aa-2.0E0*1*I_ESP_H3x2z_G4y_a-2.0E0*2*I_ESP_Hx2y2z_G4y_a+2*1*I_ESP_Fx2z_G4y;
    abcd[iGrid*1890+534] = 4.0E0*I_ESP_K3xy3z_G4y_aa-2.0E0*2*I_ESP_Hxy3z_G4y_a;
    abcd[iGrid*1890+535] = 4.0E0*I_ESP_K2x5y_G4y_aa-2.0E0*4*I_ESP_H2x3y_G4y_a-2.0E0*1*I_ESP_H5y_G4y_a+4*I_ESP_F3y_G4y;
    abcd[iGrid*1890+536] = 4.0E0*I_ESP_K2x4yz_G4y_aa-2.0E0*3*I_ESP_H2x2yz_G4y_a-2.0E0*1*I_ESP_H4yz_G4y_a+3*I_ESP_F2yz_G4y;
    abcd[iGrid*1890+537] = 4.0E0*I_ESP_K2x3y2z_G4y_aa-2.0E0*2*I_ESP_H2xy2z_G4y_a-2.0E0*1*I_ESP_H3y2z_G4y_a+2*I_ESP_Fy2z_G4y;
    abcd[iGrid*1890+538] = 4.0E0*I_ESP_K2x2y3z_G4y_aa-2.0E0*1*I_ESP_H2x3z_G4y_a-2.0E0*1*I_ESP_H2y3z_G4y_a+1*I_ESP_F3z_G4y;
    abcd[iGrid*1890+539] = 4.0E0*I_ESP_K2xy4z_G4y_aa-2.0E0*1*I_ESP_Hy4z_G4y_a;
    abcd[iGrid*1890+540] = 4.0E0*I_ESP_Kx6y_G4y_aa-2.0E0*5*I_ESP_Hx4y_G4y_a;
    abcd[iGrid*1890+541] = 4.0E0*I_ESP_Kx5yz_G4y_aa-2.0E0*4*I_ESP_Hx3yz_G4y_a;
    abcd[iGrid*1890+542] = 4.0E0*I_ESP_Kx4y2z_G4y_aa-2.0E0*3*I_ESP_Hx2y2z_G4y_a;
    abcd[iGrid*1890+543] = 4.0E0*I_ESP_Kx3y3z_G4y_aa-2.0E0*2*I_ESP_Hxy3z_G4y_a;
    abcd[iGrid*1890+544] = 4.0E0*I_ESP_Kx2y4z_G4y_aa-2.0E0*1*I_ESP_Hx4z_G4y_a;
    abcd[iGrid*1890+545] = 4.0E0*I_ESP_Kxy5z_G4y_aa;
    abcd[iGrid*1890+546] = 4.0E0*I_ESP_K6xy_G3yz_aa-2.0E0*5*I_ESP_H4xy_G3yz_a;
    abcd[iGrid*1890+547] = 4.0E0*I_ESP_K5x2y_G3yz_aa-2.0E0*1*I_ESP_H5x_G3yz_a-2.0E0*4*I_ESP_H3x2y_G3yz_a+4*1*I_ESP_F3x_G3yz;
    abcd[iGrid*1890+548] = 4.0E0*I_ESP_K5xyz_G3yz_aa-2.0E0*4*I_ESP_H3xyz_G3yz_a;
    abcd[iGrid*1890+549] = 4.0E0*I_ESP_K4x3y_G3yz_aa-2.0E0*2*I_ESP_H4xy_G3yz_a-2.0E0*3*I_ESP_H2x3y_G3yz_a+3*2*I_ESP_F2xy_G3yz;
    abcd[iGrid*1890+550] = 4.0E0*I_ESP_K4x2yz_G3yz_aa-2.0E0*1*I_ESP_H4xz_G3yz_a-2.0E0*3*I_ESP_H2x2yz_G3yz_a+3*1*I_ESP_F2xz_G3yz;
    abcd[iGrid*1890+551] = 4.0E0*I_ESP_K4xy2z_G3yz_aa-2.0E0*3*I_ESP_H2xy2z_G3yz_a;
    abcd[iGrid*1890+552] = 4.0E0*I_ESP_K3x4y_G3yz_aa-2.0E0*3*I_ESP_H3x2y_G3yz_a-2.0E0*2*I_ESP_Hx4y_G3yz_a+2*3*I_ESP_Fx2y_G3yz;
    abcd[iGrid*1890+553] = 4.0E0*I_ESP_K3x3yz_G3yz_aa-2.0E0*2*I_ESP_H3xyz_G3yz_a-2.0E0*2*I_ESP_Hx3yz_G3yz_a+2*2*I_ESP_Fxyz_G3yz;
    abcd[iGrid*1890+554] = 4.0E0*I_ESP_K3x2y2z_G3yz_aa-2.0E0*1*I_ESP_H3x2z_G3yz_a-2.0E0*2*I_ESP_Hx2y2z_G3yz_a+2*1*I_ESP_Fx2z_G3yz;
    abcd[iGrid*1890+555] = 4.0E0*I_ESP_K3xy3z_G3yz_aa-2.0E0*2*I_ESP_Hxy3z_G3yz_a;
    abcd[iGrid*1890+556] = 4.0E0*I_ESP_K2x5y_G3yz_aa-2.0E0*4*I_ESP_H2x3y_G3yz_a-2.0E0*1*I_ESP_H5y_G3yz_a+4*I_ESP_F3y_G3yz;
    abcd[iGrid*1890+557] = 4.0E0*I_ESP_K2x4yz_G3yz_aa-2.0E0*3*I_ESP_H2x2yz_G3yz_a-2.0E0*1*I_ESP_H4yz_G3yz_a+3*I_ESP_F2yz_G3yz;
    abcd[iGrid*1890+558] = 4.0E0*I_ESP_K2x3y2z_G3yz_aa-2.0E0*2*I_ESP_H2xy2z_G3yz_a-2.0E0*1*I_ESP_H3y2z_G3yz_a+2*I_ESP_Fy2z_G3yz;
    abcd[iGrid*1890+559] = 4.0E0*I_ESP_K2x2y3z_G3yz_aa-2.0E0*1*I_ESP_H2x3z_G3yz_a-2.0E0*1*I_ESP_H2y3z_G3yz_a+1*I_ESP_F3z_G3yz;
    abcd[iGrid*1890+560] = 4.0E0*I_ESP_K2xy4z_G3yz_aa-2.0E0*1*I_ESP_Hy4z_G3yz_a;
    abcd[iGrid*1890+561] = 4.0E0*I_ESP_Kx6y_G3yz_aa-2.0E0*5*I_ESP_Hx4y_G3yz_a;
    abcd[iGrid*1890+562] = 4.0E0*I_ESP_Kx5yz_G3yz_aa-2.0E0*4*I_ESP_Hx3yz_G3yz_a;
    abcd[iGrid*1890+563] = 4.0E0*I_ESP_Kx4y2z_G3yz_aa-2.0E0*3*I_ESP_Hx2y2z_G3yz_a;
    abcd[iGrid*1890+564] = 4.0E0*I_ESP_Kx3y3z_G3yz_aa-2.0E0*2*I_ESP_Hxy3z_G3yz_a;
    abcd[iGrid*1890+565] = 4.0E0*I_ESP_Kx2y4z_G3yz_aa-2.0E0*1*I_ESP_Hx4z_G3yz_a;
    abcd[iGrid*1890+566] = 4.0E0*I_ESP_Kxy5z_G3yz_aa;
    abcd[iGrid*1890+567] = 4.0E0*I_ESP_K6xy_G2y2z_aa-2.0E0*5*I_ESP_H4xy_G2y2z_a;
    abcd[iGrid*1890+568] = 4.0E0*I_ESP_K5x2y_G2y2z_aa-2.0E0*1*I_ESP_H5x_G2y2z_a-2.0E0*4*I_ESP_H3x2y_G2y2z_a+4*1*I_ESP_F3x_G2y2z;
    abcd[iGrid*1890+569] = 4.0E0*I_ESP_K5xyz_G2y2z_aa-2.0E0*4*I_ESP_H3xyz_G2y2z_a;
    abcd[iGrid*1890+570] = 4.0E0*I_ESP_K4x3y_G2y2z_aa-2.0E0*2*I_ESP_H4xy_G2y2z_a-2.0E0*3*I_ESP_H2x3y_G2y2z_a+3*2*I_ESP_F2xy_G2y2z;
    abcd[iGrid*1890+571] = 4.0E0*I_ESP_K4x2yz_G2y2z_aa-2.0E0*1*I_ESP_H4xz_G2y2z_a-2.0E0*3*I_ESP_H2x2yz_G2y2z_a+3*1*I_ESP_F2xz_G2y2z;
    abcd[iGrid*1890+572] = 4.0E0*I_ESP_K4xy2z_G2y2z_aa-2.0E0*3*I_ESP_H2xy2z_G2y2z_a;
    abcd[iGrid*1890+573] = 4.0E0*I_ESP_K3x4y_G2y2z_aa-2.0E0*3*I_ESP_H3x2y_G2y2z_a-2.0E0*2*I_ESP_Hx4y_G2y2z_a+2*3*I_ESP_Fx2y_G2y2z;
    abcd[iGrid*1890+574] = 4.0E0*I_ESP_K3x3yz_G2y2z_aa-2.0E0*2*I_ESP_H3xyz_G2y2z_a-2.0E0*2*I_ESP_Hx3yz_G2y2z_a+2*2*I_ESP_Fxyz_G2y2z;
    abcd[iGrid*1890+575] = 4.0E0*I_ESP_K3x2y2z_G2y2z_aa-2.0E0*1*I_ESP_H3x2z_G2y2z_a-2.0E0*2*I_ESP_Hx2y2z_G2y2z_a+2*1*I_ESP_Fx2z_G2y2z;
    abcd[iGrid*1890+576] = 4.0E0*I_ESP_K3xy3z_G2y2z_aa-2.0E0*2*I_ESP_Hxy3z_G2y2z_a;
    abcd[iGrid*1890+577] = 4.0E0*I_ESP_K2x5y_G2y2z_aa-2.0E0*4*I_ESP_H2x3y_G2y2z_a-2.0E0*1*I_ESP_H5y_G2y2z_a+4*I_ESP_F3y_G2y2z;
    abcd[iGrid*1890+578] = 4.0E0*I_ESP_K2x4yz_G2y2z_aa-2.0E0*3*I_ESP_H2x2yz_G2y2z_a-2.0E0*1*I_ESP_H4yz_G2y2z_a+3*I_ESP_F2yz_G2y2z;
    abcd[iGrid*1890+579] = 4.0E0*I_ESP_K2x3y2z_G2y2z_aa-2.0E0*2*I_ESP_H2xy2z_G2y2z_a-2.0E0*1*I_ESP_H3y2z_G2y2z_a+2*I_ESP_Fy2z_G2y2z;
    abcd[iGrid*1890+580] = 4.0E0*I_ESP_K2x2y3z_G2y2z_aa-2.0E0*1*I_ESP_H2x3z_G2y2z_a-2.0E0*1*I_ESP_H2y3z_G2y2z_a+1*I_ESP_F3z_G2y2z;
    abcd[iGrid*1890+581] = 4.0E0*I_ESP_K2xy4z_G2y2z_aa-2.0E0*1*I_ESP_Hy4z_G2y2z_a;
    abcd[iGrid*1890+582] = 4.0E0*I_ESP_Kx6y_G2y2z_aa-2.0E0*5*I_ESP_Hx4y_G2y2z_a;
    abcd[iGrid*1890+583] = 4.0E0*I_ESP_Kx5yz_G2y2z_aa-2.0E0*4*I_ESP_Hx3yz_G2y2z_a;
    abcd[iGrid*1890+584] = 4.0E0*I_ESP_Kx4y2z_G2y2z_aa-2.0E0*3*I_ESP_Hx2y2z_G2y2z_a;
    abcd[iGrid*1890+585] = 4.0E0*I_ESP_Kx3y3z_G2y2z_aa-2.0E0*2*I_ESP_Hxy3z_G2y2z_a;
    abcd[iGrid*1890+586] = 4.0E0*I_ESP_Kx2y4z_G2y2z_aa-2.0E0*1*I_ESP_Hx4z_G2y2z_a;
    abcd[iGrid*1890+587] = 4.0E0*I_ESP_Kxy5z_G2y2z_aa;
    abcd[iGrid*1890+588] = 4.0E0*I_ESP_K6xy_Gy3z_aa-2.0E0*5*I_ESP_H4xy_Gy3z_a;
    abcd[iGrid*1890+589] = 4.0E0*I_ESP_K5x2y_Gy3z_aa-2.0E0*1*I_ESP_H5x_Gy3z_a-2.0E0*4*I_ESP_H3x2y_Gy3z_a+4*1*I_ESP_F3x_Gy3z;
    abcd[iGrid*1890+590] = 4.0E0*I_ESP_K5xyz_Gy3z_aa-2.0E0*4*I_ESP_H3xyz_Gy3z_a;
    abcd[iGrid*1890+591] = 4.0E0*I_ESP_K4x3y_Gy3z_aa-2.0E0*2*I_ESP_H4xy_Gy3z_a-2.0E0*3*I_ESP_H2x3y_Gy3z_a+3*2*I_ESP_F2xy_Gy3z;
    abcd[iGrid*1890+592] = 4.0E0*I_ESP_K4x2yz_Gy3z_aa-2.0E0*1*I_ESP_H4xz_Gy3z_a-2.0E0*3*I_ESP_H2x2yz_Gy3z_a+3*1*I_ESP_F2xz_Gy3z;
    abcd[iGrid*1890+593] = 4.0E0*I_ESP_K4xy2z_Gy3z_aa-2.0E0*3*I_ESP_H2xy2z_Gy3z_a;
    abcd[iGrid*1890+594] = 4.0E0*I_ESP_K3x4y_Gy3z_aa-2.0E0*3*I_ESP_H3x2y_Gy3z_a-2.0E0*2*I_ESP_Hx4y_Gy3z_a+2*3*I_ESP_Fx2y_Gy3z;
    abcd[iGrid*1890+595] = 4.0E0*I_ESP_K3x3yz_Gy3z_aa-2.0E0*2*I_ESP_H3xyz_Gy3z_a-2.0E0*2*I_ESP_Hx3yz_Gy3z_a+2*2*I_ESP_Fxyz_Gy3z;
    abcd[iGrid*1890+596] = 4.0E0*I_ESP_K3x2y2z_Gy3z_aa-2.0E0*1*I_ESP_H3x2z_Gy3z_a-2.0E0*2*I_ESP_Hx2y2z_Gy3z_a+2*1*I_ESP_Fx2z_Gy3z;
    abcd[iGrid*1890+597] = 4.0E0*I_ESP_K3xy3z_Gy3z_aa-2.0E0*2*I_ESP_Hxy3z_Gy3z_a;
    abcd[iGrid*1890+598] = 4.0E0*I_ESP_K2x5y_Gy3z_aa-2.0E0*4*I_ESP_H2x3y_Gy3z_a-2.0E0*1*I_ESP_H5y_Gy3z_a+4*I_ESP_F3y_Gy3z;
    abcd[iGrid*1890+599] = 4.0E0*I_ESP_K2x4yz_Gy3z_aa-2.0E0*3*I_ESP_H2x2yz_Gy3z_a-2.0E0*1*I_ESP_H4yz_Gy3z_a+3*I_ESP_F2yz_Gy3z;
    abcd[iGrid*1890+600] = 4.0E0*I_ESP_K2x3y2z_Gy3z_aa-2.0E0*2*I_ESP_H2xy2z_Gy3z_a-2.0E0*1*I_ESP_H3y2z_Gy3z_a+2*I_ESP_Fy2z_Gy3z;
    abcd[iGrid*1890+601] = 4.0E0*I_ESP_K2x2y3z_Gy3z_aa-2.0E0*1*I_ESP_H2x3z_Gy3z_a-2.0E0*1*I_ESP_H2y3z_Gy3z_a+1*I_ESP_F3z_Gy3z;
    abcd[iGrid*1890+602] = 4.0E0*I_ESP_K2xy4z_Gy3z_aa-2.0E0*1*I_ESP_Hy4z_Gy3z_a;
    abcd[iGrid*1890+603] = 4.0E0*I_ESP_Kx6y_Gy3z_aa-2.0E0*5*I_ESP_Hx4y_Gy3z_a;
    abcd[iGrid*1890+604] = 4.0E0*I_ESP_Kx5yz_Gy3z_aa-2.0E0*4*I_ESP_Hx3yz_Gy3z_a;
    abcd[iGrid*1890+605] = 4.0E0*I_ESP_Kx4y2z_Gy3z_aa-2.0E0*3*I_ESP_Hx2y2z_Gy3z_a;
    abcd[iGrid*1890+606] = 4.0E0*I_ESP_Kx3y3z_Gy3z_aa-2.0E0*2*I_ESP_Hxy3z_Gy3z_a;
    abcd[iGrid*1890+607] = 4.0E0*I_ESP_Kx2y4z_Gy3z_aa-2.0E0*1*I_ESP_Hx4z_Gy3z_a;
    abcd[iGrid*1890+608] = 4.0E0*I_ESP_Kxy5z_Gy3z_aa;
    abcd[iGrid*1890+609] = 4.0E0*I_ESP_K6xy_G4z_aa-2.0E0*5*I_ESP_H4xy_G4z_a;
    abcd[iGrid*1890+610] = 4.0E0*I_ESP_K5x2y_G4z_aa-2.0E0*1*I_ESP_H5x_G4z_a-2.0E0*4*I_ESP_H3x2y_G4z_a+4*1*I_ESP_F3x_G4z;
    abcd[iGrid*1890+611] = 4.0E0*I_ESP_K5xyz_G4z_aa-2.0E0*4*I_ESP_H3xyz_G4z_a;
    abcd[iGrid*1890+612] = 4.0E0*I_ESP_K4x3y_G4z_aa-2.0E0*2*I_ESP_H4xy_G4z_a-2.0E0*3*I_ESP_H2x3y_G4z_a+3*2*I_ESP_F2xy_G4z;
    abcd[iGrid*1890+613] = 4.0E0*I_ESP_K4x2yz_G4z_aa-2.0E0*1*I_ESP_H4xz_G4z_a-2.0E0*3*I_ESP_H2x2yz_G4z_a+3*1*I_ESP_F2xz_G4z;
    abcd[iGrid*1890+614] = 4.0E0*I_ESP_K4xy2z_G4z_aa-2.0E0*3*I_ESP_H2xy2z_G4z_a;
    abcd[iGrid*1890+615] = 4.0E0*I_ESP_K3x4y_G4z_aa-2.0E0*3*I_ESP_H3x2y_G4z_a-2.0E0*2*I_ESP_Hx4y_G4z_a+2*3*I_ESP_Fx2y_G4z;
    abcd[iGrid*1890+616] = 4.0E0*I_ESP_K3x3yz_G4z_aa-2.0E0*2*I_ESP_H3xyz_G4z_a-2.0E0*2*I_ESP_Hx3yz_G4z_a+2*2*I_ESP_Fxyz_G4z;
    abcd[iGrid*1890+617] = 4.0E0*I_ESP_K3x2y2z_G4z_aa-2.0E0*1*I_ESP_H3x2z_G4z_a-2.0E0*2*I_ESP_Hx2y2z_G4z_a+2*1*I_ESP_Fx2z_G4z;
    abcd[iGrid*1890+618] = 4.0E0*I_ESP_K3xy3z_G4z_aa-2.0E0*2*I_ESP_Hxy3z_G4z_a;
    abcd[iGrid*1890+619] = 4.0E0*I_ESP_K2x5y_G4z_aa-2.0E0*4*I_ESP_H2x3y_G4z_a-2.0E0*1*I_ESP_H5y_G4z_a+4*I_ESP_F3y_G4z;
    abcd[iGrid*1890+620] = 4.0E0*I_ESP_K2x4yz_G4z_aa-2.0E0*3*I_ESP_H2x2yz_G4z_a-2.0E0*1*I_ESP_H4yz_G4z_a+3*I_ESP_F2yz_G4z;
    abcd[iGrid*1890+621] = 4.0E0*I_ESP_K2x3y2z_G4z_aa-2.0E0*2*I_ESP_H2xy2z_G4z_a-2.0E0*1*I_ESP_H3y2z_G4z_a+2*I_ESP_Fy2z_G4z;
    abcd[iGrid*1890+622] = 4.0E0*I_ESP_K2x2y3z_G4z_aa-2.0E0*1*I_ESP_H2x3z_G4z_a-2.0E0*1*I_ESP_H2y3z_G4z_a+1*I_ESP_F3z_G4z;
    abcd[iGrid*1890+623] = 4.0E0*I_ESP_K2xy4z_G4z_aa-2.0E0*1*I_ESP_Hy4z_G4z_a;
    abcd[iGrid*1890+624] = 4.0E0*I_ESP_Kx6y_G4z_aa-2.0E0*5*I_ESP_Hx4y_G4z_a;
    abcd[iGrid*1890+625] = 4.0E0*I_ESP_Kx5yz_G4z_aa-2.0E0*4*I_ESP_Hx3yz_G4z_a;
    abcd[iGrid*1890+626] = 4.0E0*I_ESP_Kx4y2z_G4z_aa-2.0E0*3*I_ESP_Hx2y2z_G4z_a;
    abcd[iGrid*1890+627] = 4.0E0*I_ESP_Kx3y3z_G4z_aa-2.0E0*2*I_ESP_Hxy3z_G4z_a;
    abcd[iGrid*1890+628] = 4.0E0*I_ESP_Kx2y4z_G4z_aa-2.0E0*1*I_ESP_Hx4z_G4z_a;
    abcd[iGrid*1890+629] = 4.0E0*I_ESP_Kxy5z_G4z_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_G_aa
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_F_G
     ************************************************************/
    abcd[iGrid*1890+630] = 4.0E0*I_ESP_K6xz_G4x_aa-2.0E0*5*I_ESP_H4xz_G4x_a;
    abcd[iGrid*1890+631] = 4.0E0*I_ESP_K5xyz_G4x_aa-2.0E0*4*I_ESP_H3xyz_G4x_a;
    abcd[iGrid*1890+632] = 4.0E0*I_ESP_K5x2z_G4x_aa-2.0E0*1*I_ESP_H5x_G4x_a-2.0E0*4*I_ESP_H3x2z_G4x_a+4*1*I_ESP_F3x_G4x;
    abcd[iGrid*1890+633] = 4.0E0*I_ESP_K4x2yz_G4x_aa-2.0E0*3*I_ESP_H2x2yz_G4x_a;
    abcd[iGrid*1890+634] = 4.0E0*I_ESP_K4xy2z_G4x_aa-2.0E0*1*I_ESP_H4xy_G4x_a-2.0E0*3*I_ESP_H2xy2z_G4x_a+3*1*I_ESP_F2xy_G4x;
    abcd[iGrid*1890+635] = 4.0E0*I_ESP_K4x3z_G4x_aa-2.0E0*2*I_ESP_H4xz_G4x_a-2.0E0*3*I_ESP_H2x3z_G4x_a+3*2*I_ESP_F2xz_G4x;
    abcd[iGrid*1890+636] = 4.0E0*I_ESP_K3x3yz_G4x_aa-2.0E0*2*I_ESP_Hx3yz_G4x_a;
    abcd[iGrid*1890+637] = 4.0E0*I_ESP_K3x2y2z_G4x_aa-2.0E0*1*I_ESP_H3x2y_G4x_a-2.0E0*2*I_ESP_Hx2y2z_G4x_a+2*1*I_ESP_Fx2y_G4x;
    abcd[iGrid*1890+638] = 4.0E0*I_ESP_K3xy3z_G4x_aa-2.0E0*2*I_ESP_H3xyz_G4x_a-2.0E0*2*I_ESP_Hxy3z_G4x_a+2*2*I_ESP_Fxyz_G4x;
    abcd[iGrid*1890+639] = 4.0E0*I_ESP_K3x4z_G4x_aa-2.0E0*3*I_ESP_H3x2z_G4x_a-2.0E0*2*I_ESP_Hx4z_G4x_a+2*3*I_ESP_Fx2z_G4x;
    abcd[iGrid*1890+640] = 4.0E0*I_ESP_K2x4yz_G4x_aa-2.0E0*1*I_ESP_H4yz_G4x_a;
    abcd[iGrid*1890+641] = 4.0E0*I_ESP_K2x3y2z_G4x_aa-2.0E0*1*I_ESP_H2x3y_G4x_a-2.0E0*1*I_ESP_H3y2z_G4x_a+1*I_ESP_F3y_G4x;
    abcd[iGrid*1890+642] = 4.0E0*I_ESP_K2x2y3z_G4x_aa-2.0E0*2*I_ESP_H2x2yz_G4x_a-2.0E0*1*I_ESP_H2y3z_G4x_a+2*I_ESP_F2yz_G4x;
    abcd[iGrid*1890+643] = 4.0E0*I_ESP_K2xy4z_G4x_aa-2.0E0*3*I_ESP_H2xy2z_G4x_a-2.0E0*1*I_ESP_Hy4z_G4x_a+3*I_ESP_Fy2z_G4x;
    abcd[iGrid*1890+644] = 4.0E0*I_ESP_K2x5z_G4x_aa-2.0E0*4*I_ESP_H2x3z_G4x_a-2.0E0*1*I_ESP_H5z_G4x_a+4*I_ESP_F3z_G4x;
    abcd[iGrid*1890+645] = 4.0E0*I_ESP_Kx5yz_G4x_aa;
    abcd[iGrid*1890+646] = 4.0E0*I_ESP_Kx4y2z_G4x_aa-2.0E0*1*I_ESP_Hx4y_G4x_a;
    abcd[iGrid*1890+647] = 4.0E0*I_ESP_Kx3y3z_G4x_aa-2.0E0*2*I_ESP_Hx3yz_G4x_a;
    abcd[iGrid*1890+648] = 4.0E0*I_ESP_Kx2y4z_G4x_aa-2.0E0*3*I_ESP_Hx2y2z_G4x_a;
    abcd[iGrid*1890+649] = 4.0E0*I_ESP_Kxy5z_G4x_aa-2.0E0*4*I_ESP_Hxy3z_G4x_a;
    abcd[iGrid*1890+650] = 4.0E0*I_ESP_Kx6z_G4x_aa-2.0E0*5*I_ESP_Hx4z_G4x_a;
    abcd[iGrid*1890+651] = 4.0E0*I_ESP_K6xz_G3xy_aa-2.0E0*5*I_ESP_H4xz_G3xy_a;
    abcd[iGrid*1890+652] = 4.0E0*I_ESP_K5xyz_G3xy_aa-2.0E0*4*I_ESP_H3xyz_G3xy_a;
    abcd[iGrid*1890+653] = 4.0E0*I_ESP_K5x2z_G3xy_aa-2.0E0*1*I_ESP_H5x_G3xy_a-2.0E0*4*I_ESP_H3x2z_G3xy_a+4*1*I_ESP_F3x_G3xy;
    abcd[iGrid*1890+654] = 4.0E0*I_ESP_K4x2yz_G3xy_aa-2.0E0*3*I_ESP_H2x2yz_G3xy_a;
    abcd[iGrid*1890+655] = 4.0E0*I_ESP_K4xy2z_G3xy_aa-2.0E0*1*I_ESP_H4xy_G3xy_a-2.0E0*3*I_ESP_H2xy2z_G3xy_a+3*1*I_ESP_F2xy_G3xy;
    abcd[iGrid*1890+656] = 4.0E0*I_ESP_K4x3z_G3xy_aa-2.0E0*2*I_ESP_H4xz_G3xy_a-2.0E0*3*I_ESP_H2x3z_G3xy_a+3*2*I_ESP_F2xz_G3xy;
    abcd[iGrid*1890+657] = 4.0E0*I_ESP_K3x3yz_G3xy_aa-2.0E0*2*I_ESP_Hx3yz_G3xy_a;
    abcd[iGrid*1890+658] = 4.0E0*I_ESP_K3x2y2z_G3xy_aa-2.0E0*1*I_ESP_H3x2y_G3xy_a-2.0E0*2*I_ESP_Hx2y2z_G3xy_a+2*1*I_ESP_Fx2y_G3xy;
    abcd[iGrid*1890+659] = 4.0E0*I_ESP_K3xy3z_G3xy_aa-2.0E0*2*I_ESP_H3xyz_G3xy_a-2.0E0*2*I_ESP_Hxy3z_G3xy_a+2*2*I_ESP_Fxyz_G3xy;
    abcd[iGrid*1890+660] = 4.0E0*I_ESP_K3x4z_G3xy_aa-2.0E0*3*I_ESP_H3x2z_G3xy_a-2.0E0*2*I_ESP_Hx4z_G3xy_a+2*3*I_ESP_Fx2z_G3xy;
    abcd[iGrid*1890+661] = 4.0E0*I_ESP_K2x4yz_G3xy_aa-2.0E0*1*I_ESP_H4yz_G3xy_a;
    abcd[iGrid*1890+662] = 4.0E0*I_ESP_K2x3y2z_G3xy_aa-2.0E0*1*I_ESP_H2x3y_G3xy_a-2.0E0*1*I_ESP_H3y2z_G3xy_a+1*I_ESP_F3y_G3xy;
    abcd[iGrid*1890+663] = 4.0E0*I_ESP_K2x2y3z_G3xy_aa-2.0E0*2*I_ESP_H2x2yz_G3xy_a-2.0E0*1*I_ESP_H2y3z_G3xy_a+2*I_ESP_F2yz_G3xy;
    abcd[iGrid*1890+664] = 4.0E0*I_ESP_K2xy4z_G3xy_aa-2.0E0*3*I_ESP_H2xy2z_G3xy_a-2.0E0*1*I_ESP_Hy4z_G3xy_a+3*I_ESP_Fy2z_G3xy;
    abcd[iGrid*1890+665] = 4.0E0*I_ESP_K2x5z_G3xy_aa-2.0E0*4*I_ESP_H2x3z_G3xy_a-2.0E0*1*I_ESP_H5z_G3xy_a+4*I_ESP_F3z_G3xy;
    abcd[iGrid*1890+666] = 4.0E0*I_ESP_Kx5yz_G3xy_aa;
    abcd[iGrid*1890+667] = 4.0E0*I_ESP_Kx4y2z_G3xy_aa-2.0E0*1*I_ESP_Hx4y_G3xy_a;
    abcd[iGrid*1890+668] = 4.0E0*I_ESP_Kx3y3z_G3xy_aa-2.0E0*2*I_ESP_Hx3yz_G3xy_a;
    abcd[iGrid*1890+669] = 4.0E0*I_ESP_Kx2y4z_G3xy_aa-2.0E0*3*I_ESP_Hx2y2z_G3xy_a;
    abcd[iGrid*1890+670] = 4.0E0*I_ESP_Kxy5z_G3xy_aa-2.0E0*4*I_ESP_Hxy3z_G3xy_a;
    abcd[iGrid*1890+671] = 4.0E0*I_ESP_Kx6z_G3xy_aa-2.0E0*5*I_ESP_Hx4z_G3xy_a;
    abcd[iGrid*1890+672] = 4.0E0*I_ESP_K6xz_G3xz_aa-2.0E0*5*I_ESP_H4xz_G3xz_a;
    abcd[iGrid*1890+673] = 4.0E0*I_ESP_K5xyz_G3xz_aa-2.0E0*4*I_ESP_H3xyz_G3xz_a;
    abcd[iGrid*1890+674] = 4.0E0*I_ESP_K5x2z_G3xz_aa-2.0E0*1*I_ESP_H5x_G3xz_a-2.0E0*4*I_ESP_H3x2z_G3xz_a+4*1*I_ESP_F3x_G3xz;
    abcd[iGrid*1890+675] = 4.0E0*I_ESP_K4x2yz_G3xz_aa-2.0E0*3*I_ESP_H2x2yz_G3xz_a;
    abcd[iGrid*1890+676] = 4.0E0*I_ESP_K4xy2z_G3xz_aa-2.0E0*1*I_ESP_H4xy_G3xz_a-2.0E0*3*I_ESP_H2xy2z_G3xz_a+3*1*I_ESP_F2xy_G3xz;
    abcd[iGrid*1890+677] = 4.0E0*I_ESP_K4x3z_G3xz_aa-2.0E0*2*I_ESP_H4xz_G3xz_a-2.0E0*3*I_ESP_H2x3z_G3xz_a+3*2*I_ESP_F2xz_G3xz;
    abcd[iGrid*1890+678] = 4.0E0*I_ESP_K3x3yz_G3xz_aa-2.0E0*2*I_ESP_Hx3yz_G3xz_a;
    abcd[iGrid*1890+679] = 4.0E0*I_ESP_K3x2y2z_G3xz_aa-2.0E0*1*I_ESP_H3x2y_G3xz_a-2.0E0*2*I_ESP_Hx2y2z_G3xz_a+2*1*I_ESP_Fx2y_G3xz;
    abcd[iGrid*1890+680] = 4.0E0*I_ESP_K3xy3z_G3xz_aa-2.0E0*2*I_ESP_H3xyz_G3xz_a-2.0E0*2*I_ESP_Hxy3z_G3xz_a+2*2*I_ESP_Fxyz_G3xz;
    abcd[iGrid*1890+681] = 4.0E0*I_ESP_K3x4z_G3xz_aa-2.0E0*3*I_ESP_H3x2z_G3xz_a-2.0E0*2*I_ESP_Hx4z_G3xz_a+2*3*I_ESP_Fx2z_G3xz;
    abcd[iGrid*1890+682] = 4.0E0*I_ESP_K2x4yz_G3xz_aa-2.0E0*1*I_ESP_H4yz_G3xz_a;
    abcd[iGrid*1890+683] = 4.0E0*I_ESP_K2x3y2z_G3xz_aa-2.0E0*1*I_ESP_H2x3y_G3xz_a-2.0E0*1*I_ESP_H3y2z_G3xz_a+1*I_ESP_F3y_G3xz;
    abcd[iGrid*1890+684] = 4.0E0*I_ESP_K2x2y3z_G3xz_aa-2.0E0*2*I_ESP_H2x2yz_G3xz_a-2.0E0*1*I_ESP_H2y3z_G3xz_a+2*I_ESP_F2yz_G3xz;
    abcd[iGrid*1890+685] = 4.0E0*I_ESP_K2xy4z_G3xz_aa-2.0E0*3*I_ESP_H2xy2z_G3xz_a-2.0E0*1*I_ESP_Hy4z_G3xz_a+3*I_ESP_Fy2z_G3xz;
    abcd[iGrid*1890+686] = 4.0E0*I_ESP_K2x5z_G3xz_aa-2.0E0*4*I_ESP_H2x3z_G3xz_a-2.0E0*1*I_ESP_H5z_G3xz_a+4*I_ESP_F3z_G3xz;
    abcd[iGrid*1890+687] = 4.0E0*I_ESP_Kx5yz_G3xz_aa;
    abcd[iGrid*1890+688] = 4.0E0*I_ESP_Kx4y2z_G3xz_aa-2.0E0*1*I_ESP_Hx4y_G3xz_a;
    abcd[iGrid*1890+689] = 4.0E0*I_ESP_Kx3y3z_G3xz_aa-2.0E0*2*I_ESP_Hx3yz_G3xz_a;
    abcd[iGrid*1890+690] = 4.0E0*I_ESP_Kx2y4z_G3xz_aa-2.0E0*3*I_ESP_Hx2y2z_G3xz_a;
    abcd[iGrid*1890+691] = 4.0E0*I_ESP_Kxy5z_G3xz_aa-2.0E0*4*I_ESP_Hxy3z_G3xz_a;
    abcd[iGrid*1890+692] = 4.0E0*I_ESP_Kx6z_G3xz_aa-2.0E0*5*I_ESP_Hx4z_G3xz_a;
    abcd[iGrid*1890+693] = 4.0E0*I_ESP_K6xz_G2x2y_aa-2.0E0*5*I_ESP_H4xz_G2x2y_a;
    abcd[iGrid*1890+694] = 4.0E0*I_ESP_K5xyz_G2x2y_aa-2.0E0*4*I_ESP_H3xyz_G2x2y_a;
    abcd[iGrid*1890+695] = 4.0E0*I_ESP_K5x2z_G2x2y_aa-2.0E0*1*I_ESP_H5x_G2x2y_a-2.0E0*4*I_ESP_H3x2z_G2x2y_a+4*1*I_ESP_F3x_G2x2y;
    abcd[iGrid*1890+696] = 4.0E0*I_ESP_K4x2yz_G2x2y_aa-2.0E0*3*I_ESP_H2x2yz_G2x2y_a;
    abcd[iGrid*1890+697] = 4.0E0*I_ESP_K4xy2z_G2x2y_aa-2.0E0*1*I_ESP_H4xy_G2x2y_a-2.0E0*3*I_ESP_H2xy2z_G2x2y_a+3*1*I_ESP_F2xy_G2x2y;
    abcd[iGrid*1890+698] = 4.0E0*I_ESP_K4x3z_G2x2y_aa-2.0E0*2*I_ESP_H4xz_G2x2y_a-2.0E0*3*I_ESP_H2x3z_G2x2y_a+3*2*I_ESP_F2xz_G2x2y;
    abcd[iGrid*1890+699] = 4.0E0*I_ESP_K3x3yz_G2x2y_aa-2.0E0*2*I_ESP_Hx3yz_G2x2y_a;
    abcd[iGrid*1890+700] = 4.0E0*I_ESP_K3x2y2z_G2x2y_aa-2.0E0*1*I_ESP_H3x2y_G2x2y_a-2.0E0*2*I_ESP_Hx2y2z_G2x2y_a+2*1*I_ESP_Fx2y_G2x2y;
    abcd[iGrid*1890+701] = 4.0E0*I_ESP_K3xy3z_G2x2y_aa-2.0E0*2*I_ESP_H3xyz_G2x2y_a-2.0E0*2*I_ESP_Hxy3z_G2x2y_a+2*2*I_ESP_Fxyz_G2x2y;
    abcd[iGrid*1890+702] = 4.0E0*I_ESP_K3x4z_G2x2y_aa-2.0E0*3*I_ESP_H3x2z_G2x2y_a-2.0E0*2*I_ESP_Hx4z_G2x2y_a+2*3*I_ESP_Fx2z_G2x2y;
    abcd[iGrid*1890+703] = 4.0E0*I_ESP_K2x4yz_G2x2y_aa-2.0E0*1*I_ESP_H4yz_G2x2y_a;
    abcd[iGrid*1890+704] = 4.0E0*I_ESP_K2x3y2z_G2x2y_aa-2.0E0*1*I_ESP_H2x3y_G2x2y_a-2.0E0*1*I_ESP_H3y2z_G2x2y_a+1*I_ESP_F3y_G2x2y;
    abcd[iGrid*1890+705] = 4.0E0*I_ESP_K2x2y3z_G2x2y_aa-2.0E0*2*I_ESP_H2x2yz_G2x2y_a-2.0E0*1*I_ESP_H2y3z_G2x2y_a+2*I_ESP_F2yz_G2x2y;
    abcd[iGrid*1890+706] = 4.0E0*I_ESP_K2xy4z_G2x2y_aa-2.0E0*3*I_ESP_H2xy2z_G2x2y_a-2.0E0*1*I_ESP_Hy4z_G2x2y_a+3*I_ESP_Fy2z_G2x2y;
    abcd[iGrid*1890+707] = 4.0E0*I_ESP_K2x5z_G2x2y_aa-2.0E0*4*I_ESP_H2x3z_G2x2y_a-2.0E0*1*I_ESP_H5z_G2x2y_a+4*I_ESP_F3z_G2x2y;
    abcd[iGrid*1890+708] = 4.0E0*I_ESP_Kx5yz_G2x2y_aa;
    abcd[iGrid*1890+709] = 4.0E0*I_ESP_Kx4y2z_G2x2y_aa-2.0E0*1*I_ESP_Hx4y_G2x2y_a;
    abcd[iGrid*1890+710] = 4.0E0*I_ESP_Kx3y3z_G2x2y_aa-2.0E0*2*I_ESP_Hx3yz_G2x2y_a;
    abcd[iGrid*1890+711] = 4.0E0*I_ESP_Kx2y4z_G2x2y_aa-2.0E0*3*I_ESP_Hx2y2z_G2x2y_a;
    abcd[iGrid*1890+712] = 4.0E0*I_ESP_Kxy5z_G2x2y_aa-2.0E0*4*I_ESP_Hxy3z_G2x2y_a;
    abcd[iGrid*1890+713] = 4.0E0*I_ESP_Kx6z_G2x2y_aa-2.0E0*5*I_ESP_Hx4z_G2x2y_a;
    abcd[iGrid*1890+714] = 4.0E0*I_ESP_K6xz_G2xyz_aa-2.0E0*5*I_ESP_H4xz_G2xyz_a;
    abcd[iGrid*1890+715] = 4.0E0*I_ESP_K5xyz_G2xyz_aa-2.0E0*4*I_ESP_H3xyz_G2xyz_a;
    abcd[iGrid*1890+716] = 4.0E0*I_ESP_K5x2z_G2xyz_aa-2.0E0*1*I_ESP_H5x_G2xyz_a-2.0E0*4*I_ESP_H3x2z_G2xyz_a+4*1*I_ESP_F3x_G2xyz;
    abcd[iGrid*1890+717] = 4.0E0*I_ESP_K4x2yz_G2xyz_aa-2.0E0*3*I_ESP_H2x2yz_G2xyz_a;
    abcd[iGrid*1890+718] = 4.0E0*I_ESP_K4xy2z_G2xyz_aa-2.0E0*1*I_ESP_H4xy_G2xyz_a-2.0E0*3*I_ESP_H2xy2z_G2xyz_a+3*1*I_ESP_F2xy_G2xyz;
    abcd[iGrid*1890+719] = 4.0E0*I_ESP_K4x3z_G2xyz_aa-2.0E0*2*I_ESP_H4xz_G2xyz_a-2.0E0*3*I_ESP_H2x3z_G2xyz_a+3*2*I_ESP_F2xz_G2xyz;
    abcd[iGrid*1890+720] = 4.0E0*I_ESP_K3x3yz_G2xyz_aa-2.0E0*2*I_ESP_Hx3yz_G2xyz_a;
    abcd[iGrid*1890+721] = 4.0E0*I_ESP_K3x2y2z_G2xyz_aa-2.0E0*1*I_ESP_H3x2y_G2xyz_a-2.0E0*2*I_ESP_Hx2y2z_G2xyz_a+2*1*I_ESP_Fx2y_G2xyz;
    abcd[iGrid*1890+722] = 4.0E0*I_ESP_K3xy3z_G2xyz_aa-2.0E0*2*I_ESP_H3xyz_G2xyz_a-2.0E0*2*I_ESP_Hxy3z_G2xyz_a+2*2*I_ESP_Fxyz_G2xyz;
    abcd[iGrid*1890+723] = 4.0E0*I_ESP_K3x4z_G2xyz_aa-2.0E0*3*I_ESP_H3x2z_G2xyz_a-2.0E0*2*I_ESP_Hx4z_G2xyz_a+2*3*I_ESP_Fx2z_G2xyz;
    abcd[iGrid*1890+724] = 4.0E0*I_ESP_K2x4yz_G2xyz_aa-2.0E0*1*I_ESP_H4yz_G2xyz_a;
    abcd[iGrid*1890+725] = 4.0E0*I_ESP_K2x3y2z_G2xyz_aa-2.0E0*1*I_ESP_H2x3y_G2xyz_a-2.0E0*1*I_ESP_H3y2z_G2xyz_a+1*I_ESP_F3y_G2xyz;
    abcd[iGrid*1890+726] = 4.0E0*I_ESP_K2x2y3z_G2xyz_aa-2.0E0*2*I_ESP_H2x2yz_G2xyz_a-2.0E0*1*I_ESP_H2y3z_G2xyz_a+2*I_ESP_F2yz_G2xyz;
    abcd[iGrid*1890+727] = 4.0E0*I_ESP_K2xy4z_G2xyz_aa-2.0E0*3*I_ESP_H2xy2z_G2xyz_a-2.0E0*1*I_ESP_Hy4z_G2xyz_a+3*I_ESP_Fy2z_G2xyz;
    abcd[iGrid*1890+728] = 4.0E0*I_ESP_K2x5z_G2xyz_aa-2.0E0*4*I_ESP_H2x3z_G2xyz_a-2.0E0*1*I_ESP_H5z_G2xyz_a+4*I_ESP_F3z_G2xyz;
    abcd[iGrid*1890+729] = 4.0E0*I_ESP_Kx5yz_G2xyz_aa;
    abcd[iGrid*1890+730] = 4.0E0*I_ESP_Kx4y2z_G2xyz_aa-2.0E0*1*I_ESP_Hx4y_G2xyz_a;
    abcd[iGrid*1890+731] = 4.0E0*I_ESP_Kx3y3z_G2xyz_aa-2.0E0*2*I_ESP_Hx3yz_G2xyz_a;
    abcd[iGrid*1890+732] = 4.0E0*I_ESP_Kx2y4z_G2xyz_aa-2.0E0*3*I_ESP_Hx2y2z_G2xyz_a;
    abcd[iGrid*1890+733] = 4.0E0*I_ESP_Kxy5z_G2xyz_aa-2.0E0*4*I_ESP_Hxy3z_G2xyz_a;
    abcd[iGrid*1890+734] = 4.0E0*I_ESP_Kx6z_G2xyz_aa-2.0E0*5*I_ESP_Hx4z_G2xyz_a;
    abcd[iGrid*1890+735] = 4.0E0*I_ESP_K6xz_G2x2z_aa-2.0E0*5*I_ESP_H4xz_G2x2z_a;
    abcd[iGrid*1890+736] = 4.0E0*I_ESP_K5xyz_G2x2z_aa-2.0E0*4*I_ESP_H3xyz_G2x2z_a;
    abcd[iGrid*1890+737] = 4.0E0*I_ESP_K5x2z_G2x2z_aa-2.0E0*1*I_ESP_H5x_G2x2z_a-2.0E0*4*I_ESP_H3x2z_G2x2z_a+4*1*I_ESP_F3x_G2x2z;
    abcd[iGrid*1890+738] = 4.0E0*I_ESP_K4x2yz_G2x2z_aa-2.0E0*3*I_ESP_H2x2yz_G2x2z_a;
    abcd[iGrid*1890+739] = 4.0E0*I_ESP_K4xy2z_G2x2z_aa-2.0E0*1*I_ESP_H4xy_G2x2z_a-2.0E0*3*I_ESP_H2xy2z_G2x2z_a+3*1*I_ESP_F2xy_G2x2z;
    abcd[iGrid*1890+740] = 4.0E0*I_ESP_K4x3z_G2x2z_aa-2.0E0*2*I_ESP_H4xz_G2x2z_a-2.0E0*3*I_ESP_H2x3z_G2x2z_a+3*2*I_ESP_F2xz_G2x2z;
    abcd[iGrid*1890+741] = 4.0E0*I_ESP_K3x3yz_G2x2z_aa-2.0E0*2*I_ESP_Hx3yz_G2x2z_a;
    abcd[iGrid*1890+742] = 4.0E0*I_ESP_K3x2y2z_G2x2z_aa-2.0E0*1*I_ESP_H3x2y_G2x2z_a-2.0E0*2*I_ESP_Hx2y2z_G2x2z_a+2*1*I_ESP_Fx2y_G2x2z;
    abcd[iGrid*1890+743] = 4.0E0*I_ESP_K3xy3z_G2x2z_aa-2.0E0*2*I_ESP_H3xyz_G2x2z_a-2.0E0*2*I_ESP_Hxy3z_G2x2z_a+2*2*I_ESP_Fxyz_G2x2z;
    abcd[iGrid*1890+744] = 4.0E0*I_ESP_K3x4z_G2x2z_aa-2.0E0*3*I_ESP_H3x2z_G2x2z_a-2.0E0*2*I_ESP_Hx4z_G2x2z_a+2*3*I_ESP_Fx2z_G2x2z;
    abcd[iGrid*1890+745] = 4.0E0*I_ESP_K2x4yz_G2x2z_aa-2.0E0*1*I_ESP_H4yz_G2x2z_a;
    abcd[iGrid*1890+746] = 4.0E0*I_ESP_K2x3y2z_G2x2z_aa-2.0E0*1*I_ESP_H2x3y_G2x2z_a-2.0E0*1*I_ESP_H3y2z_G2x2z_a+1*I_ESP_F3y_G2x2z;
    abcd[iGrid*1890+747] = 4.0E0*I_ESP_K2x2y3z_G2x2z_aa-2.0E0*2*I_ESP_H2x2yz_G2x2z_a-2.0E0*1*I_ESP_H2y3z_G2x2z_a+2*I_ESP_F2yz_G2x2z;
    abcd[iGrid*1890+748] = 4.0E0*I_ESP_K2xy4z_G2x2z_aa-2.0E0*3*I_ESP_H2xy2z_G2x2z_a-2.0E0*1*I_ESP_Hy4z_G2x2z_a+3*I_ESP_Fy2z_G2x2z;
    abcd[iGrid*1890+749] = 4.0E0*I_ESP_K2x5z_G2x2z_aa-2.0E0*4*I_ESP_H2x3z_G2x2z_a-2.0E0*1*I_ESP_H5z_G2x2z_a+4*I_ESP_F3z_G2x2z;
    abcd[iGrid*1890+750] = 4.0E0*I_ESP_Kx5yz_G2x2z_aa;
    abcd[iGrid*1890+751] = 4.0E0*I_ESP_Kx4y2z_G2x2z_aa-2.0E0*1*I_ESP_Hx4y_G2x2z_a;
    abcd[iGrid*1890+752] = 4.0E0*I_ESP_Kx3y3z_G2x2z_aa-2.0E0*2*I_ESP_Hx3yz_G2x2z_a;
    abcd[iGrid*1890+753] = 4.0E0*I_ESP_Kx2y4z_G2x2z_aa-2.0E0*3*I_ESP_Hx2y2z_G2x2z_a;
    abcd[iGrid*1890+754] = 4.0E0*I_ESP_Kxy5z_G2x2z_aa-2.0E0*4*I_ESP_Hxy3z_G2x2z_a;
    abcd[iGrid*1890+755] = 4.0E0*I_ESP_Kx6z_G2x2z_aa-2.0E0*5*I_ESP_Hx4z_G2x2z_a;
    abcd[iGrid*1890+756] = 4.0E0*I_ESP_K6xz_Gx3y_aa-2.0E0*5*I_ESP_H4xz_Gx3y_a;
    abcd[iGrid*1890+757] = 4.0E0*I_ESP_K5xyz_Gx3y_aa-2.0E0*4*I_ESP_H3xyz_Gx3y_a;
    abcd[iGrid*1890+758] = 4.0E0*I_ESP_K5x2z_Gx3y_aa-2.0E0*1*I_ESP_H5x_Gx3y_a-2.0E0*4*I_ESP_H3x2z_Gx3y_a+4*1*I_ESP_F3x_Gx3y;
    abcd[iGrid*1890+759] = 4.0E0*I_ESP_K4x2yz_Gx3y_aa-2.0E0*3*I_ESP_H2x2yz_Gx3y_a;
    abcd[iGrid*1890+760] = 4.0E0*I_ESP_K4xy2z_Gx3y_aa-2.0E0*1*I_ESP_H4xy_Gx3y_a-2.0E0*3*I_ESP_H2xy2z_Gx3y_a+3*1*I_ESP_F2xy_Gx3y;
    abcd[iGrid*1890+761] = 4.0E0*I_ESP_K4x3z_Gx3y_aa-2.0E0*2*I_ESP_H4xz_Gx3y_a-2.0E0*3*I_ESP_H2x3z_Gx3y_a+3*2*I_ESP_F2xz_Gx3y;
    abcd[iGrid*1890+762] = 4.0E0*I_ESP_K3x3yz_Gx3y_aa-2.0E0*2*I_ESP_Hx3yz_Gx3y_a;
    abcd[iGrid*1890+763] = 4.0E0*I_ESP_K3x2y2z_Gx3y_aa-2.0E0*1*I_ESP_H3x2y_Gx3y_a-2.0E0*2*I_ESP_Hx2y2z_Gx3y_a+2*1*I_ESP_Fx2y_Gx3y;
    abcd[iGrid*1890+764] = 4.0E0*I_ESP_K3xy3z_Gx3y_aa-2.0E0*2*I_ESP_H3xyz_Gx3y_a-2.0E0*2*I_ESP_Hxy3z_Gx3y_a+2*2*I_ESP_Fxyz_Gx3y;
    abcd[iGrid*1890+765] = 4.0E0*I_ESP_K3x4z_Gx3y_aa-2.0E0*3*I_ESP_H3x2z_Gx3y_a-2.0E0*2*I_ESP_Hx4z_Gx3y_a+2*3*I_ESP_Fx2z_Gx3y;
    abcd[iGrid*1890+766] = 4.0E0*I_ESP_K2x4yz_Gx3y_aa-2.0E0*1*I_ESP_H4yz_Gx3y_a;
    abcd[iGrid*1890+767] = 4.0E0*I_ESP_K2x3y2z_Gx3y_aa-2.0E0*1*I_ESP_H2x3y_Gx3y_a-2.0E0*1*I_ESP_H3y2z_Gx3y_a+1*I_ESP_F3y_Gx3y;
    abcd[iGrid*1890+768] = 4.0E0*I_ESP_K2x2y3z_Gx3y_aa-2.0E0*2*I_ESP_H2x2yz_Gx3y_a-2.0E0*1*I_ESP_H2y3z_Gx3y_a+2*I_ESP_F2yz_Gx3y;
    abcd[iGrid*1890+769] = 4.0E0*I_ESP_K2xy4z_Gx3y_aa-2.0E0*3*I_ESP_H2xy2z_Gx3y_a-2.0E0*1*I_ESP_Hy4z_Gx3y_a+3*I_ESP_Fy2z_Gx3y;
    abcd[iGrid*1890+770] = 4.0E0*I_ESP_K2x5z_Gx3y_aa-2.0E0*4*I_ESP_H2x3z_Gx3y_a-2.0E0*1*I_ESP_H5z_Gx3y_a+4*I_ESP_F3z_Gx3y;
    abcd[iGrid*1890+771] = 4.0E0*I_ESP_Kx5yz_Gx3y_aa;
    abcd[iGrid*1890+772] = 4.0E0*I_ESP_Kx4y2z_Gx3y_aa-2.0E0*1*I_ESP_Hx4y_Gx3y_a;
    abcd[iGrid*1890+773] = 4.0E0*I_ESP_Kx3y3z_Gx3y_aa-2.0E0*2*I_ESP_Hx3yz_Gx3y_a;
    abcd[iGrid*1890+774] = 4.0E0*I_ESP_Kx2y4z_Gx3y_aa-2.0E0*3*I_ESP_Hx2y2z_Gx3y_a;
    abcd[iGrid*1890+775] = 4.0E0*I_ESP_Kxy5z_Gx3y_aa-2.0E0*4*I_ESP_Hxy3z_Gx3y_a;
    abcd[iGrid*1890+776] = 4.0E0*I_ESP_Kx6z_Gx3y_aa-2.0E0*5*I_ESP_Hx4z_Gx3y_a;
    abcd[iGrid*1890+777] = 4.0E0*I_ESP_K6xz_Gx2yz_aa-2.0E0*5*I_ESP_H4xz_Gx2yz_a;
    abcd[iGrid*1890+778] = 4.0E0*I_ESP_K5xyz_Gx2yz_aa-2.0E0*4*I_ESP_H3xyz_Gx2yz_a;
    abcd[iGrid*1890+779] = 4.0E0*I_ESP_K5x2z_Gx2yz_aa-2.0E0*1*I_ESP_H5x_Gx2yz_a-2.0E0*4*I_ESP_H3x2z_Gx2yz_a+4*1*I_ESP_F3x_Gx2yz;
    abcd[iGrid*1890+780] = 4.0E0*I_ESP_K4x2yz_Gx2yz_aa-2.0E0*3*I_ESP_H2x2yz_Gx2yz_a;
    abcd[iGrid*1890+781] = 4.0E0*I_ESP_K4xy2z_Gx2yz_aa-2.0E0*1*I_ESP_H4xy_Gx2yz_a-2.0E0*3*I_ESP_H2xy2z_Gx2yz_a+3*1*I_ESP_F2xy_Gx2yz;
    abcd[iGrid*1890+782] = 4.0E0*I_ESP_K4x3z_Gx2yz_aa-2.0E0*2*I_ESP_H4xz_Gx2yz_a-2.0E0*3*I_ESP_H2x3z_Gx2yz_a+3*2*I_ESP_F2xz_Gx2yz;
    abcd[iGrid*1890+783] = 4.0E0*I_ESP_K3x3yz_Gx2yz_aa-2.0E0*2*I_ESP_Hx3yz_Gx2yz_a;
    abcd[iGrid*1890+784] = 4.0E0*I_ESP_K3x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_H3x2y_Gx2yz_a-2.0E0*2*I_ESP_Hx2y2z_Gx2yz_a+2*1*I_ESP_Fx2y_Gx2yz;
    abcd[iGrid*1890+785] = 4.0E0*I_ESP_K3xy3z_Gx2yz_aa-2.0E0*2*I_ESP_H3xyz_Gx2yz_a-2.0E0*2*I_ESP_Hxy3z_Gx2yz_a+2*2*I_ESP_Fxyz_Gx2yz;
    abcd[iGrid*1890+786] = 4.0E0*I_ESP_K3x4z_Gx2yz_aa-2.0E0*3*I_ESP_H3x2z_Gx2yz_a-2.0E0*2*I_ESP_Hx4z_Gx2yz_a+2*3*I_ESP_Fx2z_Gx2yz;
    abcd[iGrid*1890+787] = 4.0E0*I_ESP_K2x4yz_Gx2yz_aa-2.0E0*1*I_ESP_H4yz_Gx2yz_a;
    abcd[iGrid*1890+788] = 4.0E0*I_ESP_K2x3y2z_Gx2yz_aa-2.0E0*1*I_ESP_H2x3y_Gx2yz_a-2.0E0*1*I_ESP_H3y2z_Gx2yz_a+1*I_ESP_F3y_Gx2yz;
    abcd[iGrid*1890+789] = 4.0E0*I_ESP_K2x2y3z_Gx2yz_aa-2.0E0*2*I_ESP_H2x2yz_Gx2yz_a-2.0E0*1*I_ESP_H2y3z_Gx2yz_a+2*I_ESP_F2yz_Gx2yz;
    abcd[iGrid*1890+790] = 4.0E0*I_ESP_K2xy4z_Gx2yz_aa-2.0E0*3*I_ESP_H2xy2z_Gx2yz_a-2.0E0*1*I_ESP_Hy4z_Gx2yz_a+3*I_ESP_Fy2z_Gx2yz;
    abcd[iGrid*1890+791] = 4.0E0*I_ESP_K2x5z_Gx2yz_aa-2.0E0*4*I_ESP_H2x3z_Gx2yz_a-2.0E0*1*I_ESP_H5z_Gx2yz_a+4*I_ESP_F3z_Gx2yz;
    abcd[iGrid*1890+792] = 4.0E0*I_ESP_Kx5yz_Gx2yz_aa;
    abcd[iGrid*1890+793] = 4.0E0*I_ESP_Kx4y2z_Gx2yz_aa-2.0E0*1*I_ESP_Hx4y_Gx2yz_a;
    abcd[iGrid*1890+794] = 4.0E0*I_ESP_Kx3y3z_Gx2yz_aa-2.0E0*2*I_ESP_Hx3yz_Gx2yz_a;
    abcd[iGrid*1890+795] = 4.0E0*I_ESP_Kx2y4z_Gx2yz_aa-2.0E0*3*I_ESP_Hx2y2z_Gx2yz_a;
    abcd[iGrid*1890+796] = 4.0E0*I_ESP_Kxy5z_Gx2yz_aa-2.0E0*4*I_ESP_Hxy3z_Gx2yz_a;
    abcd[iGrid*1890+797] = 4.0E0*I_ESP_Kx6z_Gx2yz_aa-2.0E0*5*I_ESP_Hx4z_Gx2yz_a;
    abcd[iGrid*1890+798] = 4.0E0*I_ESP_K6xz_Gxy2z_aa-2.0E0*5*I_ESP_H4xz_Gxy2z_a;
    abcd[iGrid*1890+799] = 4.0E0*I_ESP_K5xyz_Gxy2z_aa-2.0E0*4*I_ESP_H3xyz_Gxy2z_a;
    abcd[iGrid*1890+800] = 4.0E0*I_ESP_K5x2z_Gxy2z_aa-2.0E0*1*I_ESP_H5x_Gxy2z_a-2.0E0*4*I_ESP_H3x2z_Gxy2z_a+4*1*I_ESP_F3x_Gxy2z;
    abcd[iGrid*1890+801] = 4.0E0*I_ESP_K4x2yz_Gxy2z_aa-2.0E0*3*I_ESP_H2x2yz_Gxy2z_a;
    abcd[iGrid*1890+802] = 4.0E0*I_ESP_K4xy2z_Gxy2z_aa-2.0E0*1*I_ESP_H4xy_Gxy2z_a-2.0E0*3*I_ESP_H2xy2z_Gxy2z_a+3*1*I_ESP_F2xy_Gxy2z;
    abcd[iGrid*1890+803] = 4.0E0*I_ESP_K4x3z_Gxy2z_aa-2.0E0*2*I_ESP_H4xz_Gxy2z_a-2.0E0*3*I_ESP_H2x3z_Gxy2z_a+3*2*I_ESP_F2xz_Gxy2z;
    abcd[iGrid*1890+804] = 4.0E0*I_ESP_K3x3yz_Gxy2z_aa-2.0E0*2*I_ESP_Hx3yz_Gxy2z_a;
    abcd[iGrid*1890+805] = 4.0E0*I_ESP_K3x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_H3x2y_Gxy2z_a-2.0E0*2*I_ESP_Hx2y2z_Gxy2z_a+2*1*I_ESP_Fx2y_Gxy2z;
    abcd[iGrid*1890+806] = 4.0E0*I_ESP_K3xy3z_Gxy2z_aa-2.0E0*2*I_ESP_H3xyz_Gxy2z_a-2.0E0*2*I_ESP_Hxy3z_Gxy2z_a+2*2*I_ESP_Fxyz_Gxy2z;
    abcd[iGrid*1890+807] = 4.0E0*I_ESP_K3x4z_Gxy2z_aa-2.0E0*3*I_ESP_H3x2z_Gxy2z_a-2.0E0*2*I_ESP_Hx4z_Gxy2z_a+2*3*I_ESP_Fx2z_Gxy2z;
    abcd[iGrid*1890+808] = 4.0E0*I_ESP_K2x4yz_Gxy2z_aa-2.0E0*1*I_ESP_H4yz_Gxy2z_a;
    abcd[iGrid*1890+809] = 4.0E0*I_ESP_K2x3y2z_Gxy2z_aa-2.0E0*1*I_ESP_H2x3y_Gxy2z_a-2.0E0*1*I_ESP_H3y2z_Gxy2z_a+1*I_ESP_F3y_Gxy2z;
    abcd[iGrid*1890+810] = 4.0E0*I_ESP_K2x2y3z_Gxy2z_aa-2.0E0*2*I_ESP_H2x2yz_Gxy2z_a-2.0E0*1*I_ESP_H2y3z_Gxy2z_a+2*I_ESP_F2yz_Gxy2z;
    abcd[iGrid*1890+811] = 4.0E0*I_ESP_K2xy4z_Gxy2z_aa-2.0E0*3*I_ESP_H2xy2z_Gxy2z_a-2.0E0*1*I_ESP_Hy4z_Gxy2z_a+3*I_ESP_Fy2z_Gxy2z;
    abcd[iGrid*1890+812] = 4.0E0*I_ESP_K2x5z_Gxy2z_aa-2.0E0*4*I_ESP_H2x3z_Gxy2z_a-2.0E0*1*I_ESP_H5z_Gxy2z_a+4*I_ESP_F3z_Gxy2z;
    abcd[iGrid*1890+813] = 4.0E0*I_ESP_Kx5yz_Gxy2z_aa;
    abcd[iGrid*1890+814] = 4.0E0*I_ESP_Kx4y2z_Gxy2z_aa-2.0E0*1*I_ESP_Hx4y_Gxy2z_a;
    abcd[iGrid*1890+815] = 4.0E0*I_ESP_Kx3y3z_Gxy2z_aa-2.0E0*2*I_ESP_Hx3yz_Gxy2z_a;
    abcd[iGrid*1890+816] = 4.0E0*I_ESP_Kx2y4z_Gxy2z_aa-2.0E0*3*I_ESP_Hx2y2z_Gxy2z_a;
    abcd[iGrid*1890+817] = 4.0E0*I_ESP_Kxy5z_Gxy2z_aa-2.0E0*4*I_ESP_Hxy3z_Gxy2z_a;
    abcd[iGrid*1890+818] = 4.0E0*I_ESP_Kx6z_Gxy2z_aa-2.0E0*5*I_ESP_Hx4z_Gxy2z_a;
    abcd[iGrid*1890+819] = 4.0E0*I_ESP_K6xz_Gx3z_aa-2.0E0*5*I_ESP_H4xz_Gx3z_a;
    abcd[iGrid*1890+820] = 4.0E0*I_ESP_K5xyz_Gx3z_aa-2.0E0*4*I_ESP_H3xyz_Gx3z_a;
    abcd[iGrid*1890+821] = 4.0E0*I_ESP_K5x2z_Gx3z_aa-2.0E0*1*I_ESP_H5x_Gx3z_a-2.0E0*4*I_ESP_H3x2z_Gx3z_a+4*1*I_ESP_F3x_Gx3z;
    abcd[iGrid*1890+822] = 4.0E0*I_ESP_K4x2yz_Gx3z_aa-2.0E0*3*I_ESP_H2x2yz_Gx3z_a;
    abcd[iGrid*1890+823] = 4.0E0*I_ESP_K4xy2z_Gx3z_aa-2.0E0*1*I_ESP_H4xy_Gx3z_a-2.0E0*3*I_ESP_H2xy2z_Gx3z_a+3*1*I_ESP_F2xy_Gx3z;
    abcd[iGrid*1890+824] = 4.0E0*I_ESP_K4x3z_Gx3z_aa-2.0E0*2*I_ESP_H4xz_Gx3z_a-2.0E0*3*I_ESP_H2x3z_Gx3z_a+3*2*I_ESP_F2xz_Gx3z;
    abcd[iGrid*1890+825] = 4.0E0*I_ESP_K3x3yz_Gx3z_aa-2.0E0*2*I_ESP_Hx3yz_Gx3z_a;
    abcd[iGrid*1890+826] = 4.0E0*I_ESP_K3x2y2z_Gx3z_aa-2.0E0*1*I_ESP_H3x2y_Gx3z_a-2.0E0*2*I_ESP_Hx2y2z_Gx3z_a+2*1*I_ESP_Fx2y_Gx3z;
    abcd[iGrid*1890+827] = 4.0E0*I_ESP_K3xy3z_Gx3z_aa-2.0E0*2*I_ESP_H3xyz_Gx3z_a-2.0E0*2*I_ESP_Hxy3z_Gx3z_a+2*2*I_ESP_Fxyz_Gx3z;
    abcd[iGrid*1890+828] = 4.0E0*I_ESP_K3x4z_Gx3z_aa-2.0E0*3*I_ESP_H3x2z_Gx3z_a-2.0E0*2*I_ESP_Hx4z_Gx3z_a+2*3*I_ESP_Fx2z_Gx3z;
    abcd[iGrid*1890+829] = 4.0E0*I_ESP_K2x4yz_Gx3z_aa-2.0E0*1*I_ESP_H4yz_Gx3z_a;
    abcd[iGrid*1890+830] = 4.0E0*I_ESP_K2x3y2z_Gx3z_aa-2.0E0*1*I_ESP_H2x3y_Gx3z_a-2.0E0*1*I_ESP_H3y2z_Gx3z_a+1*I_ESP_F3y_Gx3z;
    abcd[iGrid*1890+831] = 4.0E0*I_ESP_K2x2y3z_Gx3z_aa-2.0E0*2*I_ESP_H2x2yz_Gx3z_a-2.0E0*1*I_ESP_H2y3z_Gx3z_a+2*I_ESP_F2yz_Gx3z;
    abcd[iGrid*1890+832] = 4.0E0*I_ESP_K2xy4z_Gx3z_aa-2.0E0*3*I_ESP_H2xy2z_Gx3z_a-2.0E0*1*I_ESP_Hy4z_Gx3z_a+3*I_ESP_Fy2z_Gx3z;
    abcd[iGrid*1890+833] = 4.0E0*I_ESP_K2x5z_Gx3z_aa-2.0E0*4*I_ESP_H2x3z_Gx3z_a-2.0E0*1*I_ESP_H5z_Gx3z_a+4*I_ESP_F3z_Gx3z;
    abcd[iGrid*1890+834] = 4.0E0*I_ESP_Kx5yz_Gx3z_aa;
    abcd[iGrid*1890+835] = 4.0E0*I_ESP_Kx4y2z_Gx3z_aa-2.0E0*1*I_ESP_Hx4y_Gx3z_a;
    abcd[iGrid*1890+836] = 4.0E0*I_ESP_Kx3y3z_Gx3z_aa-2.0E0*2*I_ESP_Hx3yz_Gx3z_a;
    abcd[iGrid*1890+837] = 4.0E0*I_ESP_Kx2y4z_Gx3z_aa-2.0E0*3*I_ESP_Hx2y2z_Gx3z_a;
    abcd[iGrid*1890+838] = 4.0E0*I_ESP_Kxy5z_Gx3z_aa-2.0E0*4*I_ESP_Hxy3z_Gx3z_a;
    abcd[iGrid*1890+839] = 4.0E0*I_ESP_Kx6z_Gx3z_aa-2.0E0*5*I_ESP_Hx4z_Gx3z_a;
    abcd[iGrid*1890+840] = 4.0E0*I_ESP_K6xz_G4y_aa-2.0E0*5*I_ESP_H4xz_G4y_a;
    abcd[iGrid*1890+841] = 4.0E0*I_ESP_K5xyz_G4y_aa-2.0E0*4*I_ESP_H3xyz_G4y_a;
    abcd[iGrid*1890+842] = 4.0E0*I_ESP_K5x2z_G4y_aa-2.0E0*1*I_ESP_H5x_G4y_a-2.0E0*4*I_ESP_H3x2z_G4y_a+4*1*I_ESP_F3x_G4y;
    abcd[iGrid*1890+843] = 4.0E0*I_ESP_K4x2yz_G4y_aa-2.0E0*3*I_ESP_H2x2yz_G4y_a;
    abcd[iGrid*1890+844] = 4.0E0*I_ESP_K4xy2z_G4y_aa-2.0E0*1*I_ESP_H4xy_G4y_a-2.0E0*3*I_ESP_H2xy2z_G4y_a+3*1*I_ESP_F2xy_G4y;
    abcd[iGrid*1890+845] = 4.0E0*I_ESP_K4x3z_G4y_aa-2.0E0*2*I_ESP_H4xz_G4y_a-2.0E0*3*I_ESP_H2x3z_G4y_a+3*2*I_ESP_F2xz_G4y;
    abcd[iGrid*1890+846] = 4.0E0*I_ESP_K3x3yz_G4y_aa-2.0E0*2*I_ESP_Hx3yz_G4y_a;
    abcd[iGrid*1890+847] = 4.0E0*I_ESP_K3x2y2z_G4y_aa-2.0E0*1*I_ESP_H3x2y_G4y_a-2.0E0*2*I_ESP_Hx2y2z_G4y_a+2*1*I_ESP_Fx2y_G4y;
    abcd[iGrid*1890+848] = 4.0E0*I_ESP_K3xy3z_G4y_aa-2.0E0*2*I_ESP_H3xyz_G4y_a-2.0E0*2*I_ESP_Hxy3z_G4y_a+2*2*I_ESP_Fxyz_G4y;
    abcd[iGrid*1890+849] = 4.0E0*I_ESP_K3x4z_G4y_aa-2.0E0*3*I_ESP_H3x2z_G4y_a-2.0E0*2*I_ESP_Hx4z_G4y_a+2*3*I_ESP_Fx2z_G4y;
    abcd[iGrid*1890+850] = 4.0E0*I_ESP_K2x4yz_G4y_aa-2.0E0*1*I_ESP_H4yz_G4y_a;
    abcd[iGrid*1890+851] = 4.0E0*I_ESP_K2x3y2z_G4y_aa-2.0E0*1*I_ESP_H2x3y_G4y_a-2.0E0*1*I_ESP_H3y2z_G4y_a+1*I_ESP_F3y_G4y;
    abcd[iGrid*1890+852] = 4.0E0*I_ESP_K2x2y3z_G4y_aa-2.0E0*2*I_ESP_H2x2yz_G4y_a-2.0E0*1*I_ESP_H2y3z_G4y_a+2*I_ESP_F2yz_G4y;
    abcd[iGrid*1890+853] = 4.0E0*I_ESP_K2xy4z_G4y_aa-2.0E0*3*I_ESP_H2xy2z_G4y_a-2.0E0*1*I_ESP_Hy4z_G4y_a+3*I_ESP_Fy2z_G4y;
    abcd[iGrid*1890+854] = 4.0E0*I_ESP_K2x5z_G4y_aa-2.0E0*4*I_ESP_H2x3z_G4y_a-2.0E0*1*I_ESP_H5z_G4y_a+4*I_ESP_F3z_G4y;
    abcd[iGrid*1890+855] = 4.0E0*I_ESP_Kx5yz_G4y_aa;
    abcd[iGrid*1890+856] = 4.0E0*I_ESP_Kx4y2z_G4y_aa-2.0E0*1*I_ESP_Hx4y_G4y_a;
    abcd[iGrid*1890+857] = 4.0E0*I_ESP_Kx3y3z_G4y_aa-2.0E0*2*I_ESP_Hx3yz_G4y_a;
    abcd[iGrid*1890+858] = 4.0E0*I_ESP_Kx2y4z_G4y_aa-2.0E0*3*I_ESP_Hx2y2z_G4y_a;
    abcd[iGrid*1890+859] = 4.0E0*I_ESP_Kxy5z_G4y_aa-2.0E0*4*I_ESP_Hxy3z_G4y_a;
    abcd[iGrid*1890+860] = 4.0E0*I_ESP_Kx6z_G4y_aa-2.0E0*5*I_ESP_Hx4z_G4y_a;
    abcd[iGrid*1890+861] = 4.0E0*I_ESP_K6xz_G3yz_aa-2.0E0*5*I_ESP_H4xz_G3yz_a;
    abcd[iGrid*1890+862] = 4.0E0*I_ESP_K5xyz_G3yz_aa-2.0E0*4*I_ESP_H3xyz_G3yz_a;
    abcd[iGrid*1890+863] = 4.0E0*I_ESP_K5x2z_G3yz_aa-2.0E0*1*I_ESP_H5x_G3yz_a-2.0E0*4*I_ESP_H3x2z_G3yz_a+4*1*I_ESP_F3x_G3yz;
    abcd[iGrid*1890+864] = 4.0E0*I_ESP_K4x2yz_G3yz_aa-2.0E0*3*I_ESP_H2x2yz_G3yz_a;
    abcd[iGrid*1890+865] = 4.0E0*I_ESP_K4xy2z_G3yz_aa-2.0E0*1*I_ESP_H4xy_G3yz_a-2.0E0*3*I_ESP_H2xy2z_G3yz_a+3*1*I_ESP_F2xy_G3yz;
    abcd[iGrid*1890+866] = 4.0E0*I_ESP_K4x3z_G3yz_aa-2.0E0*2*I_ESP_H4xz_G3yz_a-2.0E0*3*I_ESP_H2x3z_G3yz_a+3*2*I_ESP_F2xz_G3yz;
    abcd[iGrid*1890+867] = 4.0E0*I_ESP_K3x3yz_G3yz_aa-2.0E0*2*I_ESP_Hx3yz_G3yz_a;
    abcd[iGrid*1890+868] = 4.0E0*I_ESP_K3x2y2z_G3yz_aa-2.0E0*1*I_ESP_H3x2y_G3yz_a-2.0E0*2*I_ESP_Hx2y2z_G3yz_a+2*1*I_ESP_Fx2y_G3yz;
    abcd[iGrid*1890+869] = 4.0E0*I_ESP_K3xy3z_G3yz_aa-2.0E0*2*I_ESP_H3xyz_G3yz_a-2.0E0*2*I_ESP_Hxy3z_G3yz_a+2*2*I_ESP_Fxyz_G3yz;
    abcd[iGrid*1890+870] = 4.0E0*I_ESP_K3x4z_G3yz_aa-2.0E0*3*I_ESP_H3x2z_G3yz_a-2.0E0*2*I_ESP_Hx4z_G3yz_a+2*3*I_ESP_Fx2z_G3yz;
    abcd[iGrid*1890+871] = 4.0E0*I_ESP_K2x4yz_G3yz_aa-2.0E0*1*I_ESP_H4yz_G3yz_a;
    abcd[iGrid*1890+872] = 4.0E0*I_ESP_K2x3y2z_G3yz_aa-2.0E0*1*I_ESP_H2x3y_G3yz_a-2.0E0*1*I_ESP_H3y2z_G3yz_a+1*I_ESP_F3y_G3yz;
    abcd[iGrid*1890+873] = 4.0E0*I_ESP_K2x2y3z_G3yz_aa-2.0E0*2*I_ESP_H2x2yz_G3yz_a-2.0E0*1*I_ESP_H2y3z_G3yz_a+2*I_ESP_F2yz_G3yz;
    abcd[iGrid*1890+874] = 4.0E0*I_ESP_K2xy4z_G3yz_aa-2.0E0*3*I_ESP_H2xy2z_G3yz_a-2.0E0*1*I_ESP_Hy4z_G3yz_a+3*I_ESP_Fy2z_G3yz;
    abcd[iGrid*1890+875] = 4.0E0*I_ESP_K2x5z_G3yz_aa-2.0E0*4*I_ESP_H2x3z_G3yz_a-2.0E0*1*I_ESP_H5z_G3yz_a+4*I_ESP_F3z_G3yz;
    abcd[iGrid*1890+876] = 4.0E0*I_ESP_Kx5yz_G3yz_aa;
    abcd[iGrid*1890+877] = 4.0E0*I_ESP_Kx4y2z_G3yz_aa-2.0E0*1*I_ESP_Hx4y_G3yz_a;
    abcd[iGrid*1890+878] = 4.0E0*I_ESP_Kx3y3z_G3yz_aa-2.0E0*2*I_ESP_Hx3yz_G3yz_a;
    abcd[iGrid*1890+879] = 4.0E0*I_ESP_Kx2y4z_G3yz_aa-2.0E0*3*I_ESP_Hx2y2z_G3yz_a;
    abcd[iGrid*1890+880] = 4.0E0*I_ESP_Kxy5z_G3yz_aa-2.0E0*4*I_ESP_Hxy3z_G3yz_a;
    abcd[iGrid*1890+881] = 4.0E0*I_ESP_Kx6z_G3yz_aa-2.0E0*5*I_ESP_Hx4z_G3yz_a;
    abcd[iGrid*1890+882] = 4.0E0*I_ESP_K6xz_G2y2z_aa-2.0E0*5*I_ESP_H4xz_G2y2z_a;
    abcd[iGrid*1890+883] = 4.0E0*I_ESP_K5xyz_G2y2z_aa-2.0E0*4*I_ESP_H3xyz_G2y2z_a;
    abcd[iGrid*1890+884] = 4.0E0*I_ESP_K5x2z_G2y2z_aa-2.0E0*1*I_ESP_H5x_G2y2z_a-2.0E0*4*I_ESP_H3x2z_G2y2z_a+4*1*I_ESP_F3x_G2y2z;
    abcd[iGrid*1890+885] = 4.0E0*I_ESP_K4x2yz_G2y2z_aa-2.0E0*3*I_ESP_H2x2yz_G2y2z_a;
    abcd[iGrid*1890+886] = 4.0E0*I_ESP_K4xy2z_G2y2z_aa-2.0E0*1*I_ESP_H4xy_G2y2z_a-2.0E0*3*I_ESP_H2xy2z_G2y2z_a+3*1*I_ESP_F2xy_G2y2z;
    abcd[iGrid*1890+887] = 4.0E0*I_ESP_K4x3z_G2y2z_aa-2.0E0*2*I_ESP_H4xz_G2y2z_a-2.0E0*3*I_ESP_H2x3z_G2y2z_a+3*2*I_ESP_F2xz_G2y2z;
    abcd[iGrid*1890+888] = 4.0E0*I_ESP_K3x3yz_G2y2z_aa-2.0E0*2*I_ESP_Hx3yz_G2y2z_a;
    abcd[iGrid*1890+889] = 4.0E0*I_ESP_K3x2y2z_G2y2z_aa-2.0E0*1*I_ESP_H3x2y_G2y2z_a-2.0E0*2*I_ESP_Hx2y2z_G2y2z_a+2*1*I_ESP_Fx2y_G2y2z;
    abcd[iGrid*1890+890] = 4.0E0*I_ESP_K3xy3z_G2y2z_aa-2.0E0*2*I_ESP_H3xyz_G2y2z_a-2.0E0*2*I_ESP_Hxy3z_G2y2z_a+2*2*I_ESP_Fxyz_G2y2z;
    abcd[iGrid*1890+891] = 4.0E0*I_ESP_K3x4z_G2y2z_aa-2.0E0*3*I_ESP_H3x2z_G2y2z_a-2.0E0*2*I_ESP_Hx4z_G2y2z_a+2*3*I_ESP_Fx2z_G2y2z;
    abcd[iGrid*1890+892] = 4.0E0*I_ESP_K2x4yz_G2y2z_aa-2.0E0*1*I_ESP_H4yz_G2y2z_a;
    abcd[iGrid*1890+893] = 4.0E0*I_ESP_K2x3y2z_G2y2z_aa-2.0E0*1*I_ESP_H2x3y_G2y2z_a-2.0E0*1*I_ESP_H3y2z_G2y2z_a+1*I_ESP_F3y_G2y2z;
    abcd[iGrid*1890+894] = 4.0E0*I_ESP_K2x2y3z_G2y2z_aa-2.0E0*2*I_ESP_H2x2yz_G2y2z_a-2.0E0*1*I_ESP_H2y3z_G2y2z_a+2*I_ESP_F2yz_G2y2z;
    abcd[iGrid*1890+895] = 4.0E0*I_ESP_K2xy4z_G2y2z_aa-2.0E0*3*I_ESP_H2xy2z_G2y2z_a-2.0E0*1*I_ESP_Hy4z_G2y2z_a+3*I_ESP_Fy2z_G2y2z;
    abcd[iGrid*1890+896] = 4.0E0*I_ESP_K2x5z_G2y2z_aa-2.0E0*4*I_ESP_H2x3z_G2y2z_a-2.0E0*1*I_ESP_H5z_G2y2z_a+4*I_ESP_F3z_G2y2z;
    abcd[iGrid*1890+897] = 4.0E0*I_ESP_Kx5yz_G2y2z_aa;
    abcd[iGrid*1890+898] = 4.0E0*I_ESP_Kx4y2z_G2y2z_aa-2.0E0*1*I_ESP_Hx4y_G2y2z_a;
    abcd[iGrid*1890+899] = 4.0E0*I_ESP_Kx3y3z_G2y2z_aa-2.0E0*2*I_ESP_Hx3yz_G2y2z_a;
    abcd[iGrid*1890+900] = 4.0E0*I_ESP_Kx2y4z_G2y2z_aa-2.0E0*3*I_ESP_Hx2y2z_G2y2z_a;
    abcd[iGrid*1890+901] = 4.0E0*I_ESP_Kxy5z_G2y2z_aa-2.0E0*4*I_ESP_Hxy3z_G2y2z_a;
    abcd[iGrid*1890+902] = 4.0E0*I_ESP_Kx6z_G2y2z_aa-2.0E0*5*I_ESP_Hx4z_G2y2z_a;
    abcd[iGrid*1890+903] = 4.0E0*I_ESP_K6xz_Gy3z_aa-2.0E0*5*I_ESP_H4xz_Gy3z_a;
    abcd[iGrid*1890+904] = 4.0E0*I_ESP_K5xyz_Gy3z_aa-2.0E0*4*I_ESP_H3xyz_Gy3z_a;
    abcd[iGrid*1890+905] = 4.0E0*I_ESP_K5x2z_Gy3z_aa-2.0E0*1*I_ESP_H5x_Gy3z_a-2.0E0*4*I_ESP_H3x2z_Gy3z_a+4*1*I_ESP_F3x_Gy3z;
    abcd[iGrid*1890+906] = 4.0E0*I_ESP_K4x2yz_Gy3z_aa-2.0E0*3*I_ESP_H2x2yz_Gy3z_a;
    abcd[iGrid*1890+907] = 4.0E0*I_ESP_K4xy2z_Gy3z_aa-2.0E0*1*I_ESP_H4xy_Gy3z_a-2.0E0*3*I_ESP_H2xy2z_Gy3z_a+3*1*I_ESP_F2xy_Gy3z;
    abcd[iGrid*1890+908] = 4.0E0*I_ESP_K4x3z_Gy3z_aa-2.0E0*2*I_ESP_H4xz_Gy3z_a-2.0E0*3*I_ESP_H2x3z_Gy3z_a+3*2*I_ESP_F2xz_Gy3z;
    abcd[iGrid*1890+909] = 4.0E0*I_ESP_K3x3yz_Gy3z_aa-2.0E0*2*I_ESP_Hx3yz_Gy3z_a;
    abcd[iGrid*1890+910] = 4.0E0*I_ESP_K3x2y2z_Gy3z_aa-2.0E0*1*I_ESP_H3x2y_Gy3z_a-2.0E0*2*I_ESP_Hx2y2z_Gy3z_a+2*1*I_ESP_Fx2y_Gy3z;
    abcd[iGrid*1890+911] = 4.0E0*I_ESP_K3xy3z_Gy3z_aa-2.0E0*2*I_ESP_H3xyz_Gy3z_a-2.0E0*2*I_ESP_Hxy3z_Gy3z_a+2*2*I_ESP_Fxyz_Gy3z;
    abcd[iGrid*1890+912] = 4.0E0*I_ESP_K3x4z_Gy3z_aa-2.0E0*3*I_ESP_H3x2z_Gy3z_a-2.0E0*2*I_ESP_Hx4z_Gy3z_a+2*3*I_ESP_Fx2z_Gy3z;
    abcd[iGrid*1890+913] = 4.0E0*I_ESP_K2x4yz_Gy3z_aa-2.0E0*1*I_ESP_H4yz_Gy3z_a;
    abcd[iGrid*1890+914] = 4.0E0*I_ESP_K2x3y2z_Gy3z_aa-2.0E0*1*I_ESP_H2x3y_Gy3z_a-2.0E0*1*I_ESP_H3y2z_Gy3z_a+1*I_ESP_F3y_Gy3z;
    abcd[iGrid*1890+915] = 4.0E0*I_ESP_K2x2y3z_Gy3z_aa-2.0E0*2*I_ESP_H2x2yz_Gy3z_a-2.0E0*1*I_ESP_H2y3z_Gy3z_a+2*I_ESP_F2yz_Gy3z;
    abcd[iGrid*1890+916] = 4.0E0*I_ESP_K2xy4z_Gy3z_aa-2.0E0*3*I_ESP_H2xy2z_Gy3z_a-2.0E0*1*I_ESP_Hy4z_Gy3z_a+3*I_ESP_Fy2z_Gy3z;
    abcd[iGrid*1890+917] = 4.0E0*I_ESP_K2x5z_Gy3z_aa-2.0E0*4*I_ESP_H2x3z_Gy3z_a-2.0E0*1*I_ESP_H5z_Gy3z_a+4*I_ESP_F3z_Gy3z;
    abcd[iGrid*1890+918] = 4.0E0*I_ESP_Kx5yz_Gy3z_aa;
    abcd[iGrid*1890+919] = 4.0E0*I_ESP_Kx4y2z_Gy3z_aa-2.0E0*1*I_ESP_Hx4y_Gy3z_a;
    abcd[iGrid*1890+920] = 4.0E0*I_ESP_Kx3y3z_Gy3z_aa-2.0E0*2*I_ESP_Hx3yz_Gy3z_a;
    abcd[iGrid*1890+921] = 4.0E0*I_ESP_Kx2y4z_Gy3z_aa-2.0E0*3*I_ESP_Hx2y2z_Gy3z_a;
    abcd[iGrid*1890+922] = 4.0E0*I_ESP_Kxy5z_Gy3z_aa-2.0E0*4*I_ESP_Hxy3z_Gy3z_a;
    abcd[iGrid*1890+923] = 4.0E0*I_ESP_Kx6z_Gy3z_aa-2.0E0*5*I_ESP_Hx4z_Gy3z_a;
    abcd[iGrid*1890+924] = 4.0E0*I_ESP_K6xz_G4z_aa-2.0E0*5*I_ESP_H4xz_G4z_a;
    abcd[iGrid*1890+925] = 4.0E0*I_ESP_K5xyz_G4z_aa-2.0E0*4*I_ESP_H3xyz_G4z_a;
    abcd[iGrid*1890+926] = 4.0E0*I_ESP_K5x2z_G4z_aa-2.0E0*1*I_ESP_H5x_G4z_a-2.0E0*4*I_ESP_H3x2z_G4z_a+4*1*I_ESP_F3x_G4z;
    abcd[iGrid*1890+927] = 4.0E0*I_ESP_K4x2yz_G4z_aa-2.0E0*3*I_ESP_H2x2yz_G4z_a;
    abcd[iGrid*1890+928] = 4.0E0*I_ESP_K4xy2z_G4z_aa-2.0E0*1*I_ESP_H4xy_G4z_a-2.0E0*3*I_ESP_H2xy2z_G4z_a+3*1*I_ESP_F2xy_G4z;
    abcd[iGrid*1890+929] = 4.0E0*I_ESP_K4x3z_G4z_aa-2.0E0*2*I_ESP_H4xz_G4z_a-2.0E0*3*I_ESP_H2x3z_G4z_a+3*2*I_ESP_F2xz_G4z;
    abcd[iGrid*1890+930] = 4.0E0*I_ESP_K3x3yz_G4z_aa-2.0E0*2*I_ESP_Hx3yz_G4z_a;
    abcd[iGrid*1890+931] = 4.0E0*I_ESP_K3x2y2z_G4z_aa-2.0E0*1*I_ESP_H3x2y_G4z_a-2.0E0*2*I_ESP_Hx2y2z_G4z_a+2*1*I_ESP_Fx2y_G4z;
    abcd[iGrid*1890+932] = 4.0E0*I_ESP_K3xy3z_G4z_aa-2.0E0*2*I_ESP_H3xyz_G4z_a-2.0E0*2*I_ESP_Hxy3z_G4z_a+2*2*I_ESP_Fxyz_G4z;
    abcd[iGrid*1890+933] = 4.0E0*I_ESP_K3x4z_G4z_aa-2.0E0*3*I_ESP_H3x2z_G4z_a-2.0E0*2*I_ESP_Hx4z_G4z_a+2*3*I_ESP_Fx2z_G4z;
    abcd[iGrid*1890+934] = 4.0E0*I_ESP_K2x4yz_G4z_aa-2.0E0*1*I_ESP_H4yz_G4z_a;
    abcd[iGrid*1890+935] = 4.0E0*I_ESP_K2x3y2z_G4z_aa-2.0E0*1*I_ESP_H2x3y_G4z_a-2.0E0*1*I_ESP_H3y2z_G4z_a+1*I_ESP_F3y_G4z;
    abcd[iGrid*1890+936] = 4.0E0*I_ESP_K2x2y3z_G4z_aa-2.0E0*2*I_ESP_H2x2yz_G4z_a-2.0E0*1*I_ESP_H2y3z_G4z_a+2*I_ESP_F2yz_G4z;
    abcd[iGrid*1890+937] = 4.0E0*I_ESP_K2xy4z_G4z_aa-2.0E0*3*I_ESP_H2xy2z_G4z_a-2.0E0*1*I_ESP_Hy4z_G4z_a+3*I_ESP_Fy2z_G4z;
    abcd[iGrid*1890+938] = 4.0E0*I_ESP_K2x5z_G4z_aa-2.0E0*4*I_ESP_H2x3z_G4z_a-2.0E0*1*I_ESP_H5z_G4z_a+4*I_ESP_F3z_G4z;
    abcd[iGrid*1890+939] = 4.0E0*I_ESP_Kx5yz_G4z_aa;
    abcd[iGrid*1890+940] = 4.0E0*I_ESP_Kx4y2z_G4z_aa-2.0E0*1*I_ESP_Hx4y_G4z_a;
    abcd[iGrid*1890+941] = 4.0E0*I_ESP_Kx3y3z_G4z_aa-2.0E0*2*I_ESP_Hx3yz_G4z_a;
    abcd[iGrid*1890+942] = 4.0E0*I_ESP_Kx2y4z_G4z_aa-2.0E0*3*I_ESP_Hx2y2z_G4z_a;
    abcd[iGrid*1890+943] = 4.0E0*I_ESP_Kxy5z_G4z_aa-2.0E0*4*I_ESP_Hxy3z_G4z_a;
    abcd[iGrid*1890+944] = 4.0E0*I_ESP_Kx6z_G4z_aa-2.0E0*5*I_ESP_Hx4z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_G_aa
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_F_G
     ************************************************************/
    abcd[iGrid*1890+945] = 4.0E0*I_ESP_K5x2y_G4x_aa-2.0E0*1*I_ESP_H5x_G4x_a;
    abcd[iGrid*1890+946] = 4.0E0*I_ESP_K4x3y_G4x_aa-2.0E0*1*I_ESP_H4xy_G4x_a-2.0E0*2*I_ESP_H4xy_G4x_a;
    abcd[iGrid*1890+947] = 4.0E0*I_ESP_K4x2yz_G4x_aa-2.0E0*1*I_ESP_H4xz_G4x_a;
    abcd[iGrid*1890+948] = 4.0E0*I_ESP_K3x4y_G4x_aa-2.0E0*2*I_ESP_H3x2y_G4x_a-2.0E0*3*I_ESP_H3x2y_G4x_a+2*1*I_ESP_F3x_G4x;
    abcd[iGrid*1890+949] = 4.0E0*I_ESP_K3x3yz_G4x_aa-2.0E0*1*I_ESP_H3xyz_G4x_a-2.0E0*2*I_ESP_H3xyz_G4x_a;
    abcd[iGrid*1890+950] = 4.0E0*I_ESP_K3x2y2z_G4x_aa-2.0E0*1*I_ESP_H3x2z_G4x_a;
    abcd[iGrid*1890+951] = 4.0E0*I_ESP_K2x5y_G4x_aa-2.0E0*3*I_ESP_H2x3y_G4x_a-2.0E0*4*I_ESP_H2x3y_G4x_a+3*2*I_ESP_F2xy_G4x;
    abcd[iGrid*1890+952] = 4.0E0*I_ESP_K2x4yz_G4x_aa-2.0E0*2*I_ESP_H2x2yz_G4x_a-2.0E0*3*I_ESP_H2x2yz_G4x_a+2*1*I_ESP_F2xz_G4x;
    abcd[iGrid*1890+953] = 4.0E0*I_ESP_K2x3y2z_G4x_aa-2.0E0*1*I_ESP_H2xy2z_G4x_a-2.0E0*2*I_ESP_H2xy2z_G4x_a;
    abcd[iGrid*1890+954] = 4.0E0*I_ESP_K2x2y3z_G4x_aa-2.0E0*1*I_ESP_H2x3z_G4x_a;
    abcd[iGrid*1890+955] = 4.0E0*I_ESP_Kx6y_G4x_aa-2.0E0*4*I_ESP_Hx4y_G4x_a-2.0E0*5*I_ESP_Hx4y_G4x_a+4*3*I_ESP_Fx2y_G4x;
    abcd[iGrid*1890+956] = 4.0E0*I_ESP_Kx5yz_G4x_aa-2.0E0*3*I_ESP_Hx3yz_G4x_a-2.0E0*4*I_ESP_Hx3yz_G4x_a+3*2*I_ESP_Fxyz_G4x;
    abcd[iGrid*1890+957] = 4.0E0*I_ESP_Kx4y2z_G4x_aa-2.0E0*2*I_ESP_Hx2y2z_G4x_a-2.0E0*3*I_ESP_Hx2y2z_G4x_a+2*1*I_ESP_Fx2z_G4x;
    abcd[iGrid*1890+958] = 4.0E0*I_ESP_Kx3y3z_G4x_aa-2.0E0*1*I_ESP_Hxy3z_G4x_a-2.0E0*2*I_ESP_Hxy3z_G4x_a;
    abcd[iGrid*1890+959] = 4.0E0*I_ESP_Kx2y4z_G4x_aa-2.0E0*1*I_ESP_Hx4z_G4x_a;
    abcd[iGrid*1890+960] = 4.0E0*I_ESP_K7y_G4x_aa-2.0E0*5*I_ESP_H5y_G4x_a-2.0E0*6*I_ESP_H5y_G4x_a+5*4*I_ESP_F3y_G4x;
    abcd[iGrid*1890+961] = 4.0E0*I_ESP_K6yz_G4x_aa-2.0E0*4*I_ESP_H4yz_G4x_a-2.0E0*5*I_ESP_H4yz_G4x_a+4*3*I_ESP_F2yz_G4x;
    abcd[iGrid*1890+962] = 4.0E0*I_ESP_K5y2z_G4x_aa-2.0E0*3*I_ESP_H3y2z_G4x_a-2.0E0*4*I_ESP_H3y2z_G4x_a+3*2*I_ESP_Fy2z_G4x;
    abcd[iGrid*1890+963] = 4.0E0*I_ESP_K4y3z_G4x_aa-2.0E0*2*I_ESP_H2y3z_G4x_a-2.0E0*3*I_ESP_H2y3z_G4x_a+2*1*I_ESP_F3z_G4x;
    abcd[iGrid*1890+964] = 4.0E0*I_ESP_K3y4z_G4x_aa-2.0E0*1*I_ESP_Hy4z_G4x_a-2.0E0*2*I_ESP_Hy4z_G4x_a;
    abcd[iGrid*1890+965] = 4.0E0*I_ESP_K2y5z_G4x_aa-2.0E0*1*I_ESP_H5z_G4x_a;
    abcd[iGrid*1890+966] = 4.0E0*I_ESP_K5x2y_G3xy_aa-2.0E0*1*I_ESP_H5x_G3xy_a;
    abcd[iGrid*1890+967] = 4.0E0*I_ESP_K4x3y_G3xy_aa-2.0E0*1*I_ESP_H4xy_G3xy_a-2.0E0*2*I_ESP_H4xy_G3xy_a;
    abcd[iGrid*1890+968] = 4.0E0*I_ESP_K4x2yz_G3xy_aa-2.0E0*1*I_ESP_H4xz_G3xy_a;
    abcd[iGrid*1890+969] = 4.0E0*I_ESP_K3x4y_G3xy_aa-2.0E0*2*I_ESP_H3x2y_G3xy_a-2.0E0*3*I_ESP_H3x2y_G3xy_a+2*1*I_ESP_F3x_G3xy;
    abcd[iGrid*1890+970] = 4.0E0*I_ESP_K3x3yz_G3xy_aa-2.0E0*1*I_ESP_H3xyz_G3xy_a-2.0E0*2*I_ESP_H3xyz_G3xy_a;
    abcd[iGrid*1890+971] = 4.0E0*I_ESP_K3x2y2z_G3xy_aa-2.0E0*1*I_ESP_H3x2z_G3xy_a;
    abcd[iGrid*1890+972] = 4.0E0*I_ESP_K2x5y_G3xy_aa-2.0E0*3*I_ESP_H2x3y_G3xy_a-2.0E0*4*I_ESP_H2x3y_G3xy_a+3*2*I_ESP_F2xy_G3xy;
    abcd[iGrid*1890+973] = 4.0E0*I_ESP_K2x4yz_G3xy_aa-2.0E0*2*I_ESP_H2x2yz_G3xy_a-2.0E0*3*I_ESP_H2x2yz_G3xy_a+2*1*I_ESP_F2xz_G3xy;
    abcd[iGrid*1890+974] = 4.0E0*I_ESP_K2x3y2z_G3xy_aa-2.0E0*1*I_ESP_H2xy2z_G3xy_a-2.0E0*2*I_ESP_H2xy2z_G3xy_a;
    abcd[iGrid*1890+975] = 4.0E0*I_ESP_K2x2y3z_G3xy_aa-2.0E0*1*I_ESP_H2x3z_G3xy_a;
    abcd[iGrid*1890+976] = 4.0E0*I_ESP_Kx6y_G3xy_aa-2.0E0*4*I_ESP_Hx4y_G3xy_a-2.0E0*5*I_ESP_Hx4y_G3xy_a+4*3*I_ESP_Fx2y_G3xy;
    abcd[iGrid*1890+977] = 4.0E0*I_ESP_Kx5yz_G3xy_aa-2.0E0*3*I_ESP_Hx3yz_G3xy_a-2.0E0*4*I_ESP_Hx3yz_G3xy_a+3*2*I_ESP_Fxyz_G3xy;
    abcd[iGrid*1890+978] = 4.0E0*I_ESP_Kx4y2z_G3xy_aa-2.0E0*2*I_ESP_Hx2y2z_G3xy_a-2.0E0*3*I_ESP_Hx2y2z_G3xy_a+2*1*I_ESP_Fx2z_G3xy;
    abcd[iGrid*1890+979] = 4.0E0*I_ESP_Kx3y3z_G3xy_aa-2.0E0*1*I_ESP_Hxy3z_G3xy_a-2.0E0*2*I_ESP_Hxy3z_G3xy_a;
    abcd[iGrid*1890+980] = 4.0E0*I_ESP_Kx2y4z_G3xy_aa-2.0E0*1*I_ESP_Hx4z_G3xy_a;
    abcd[iGrid*1890+981] = 4.0E0*I_ESP_K7y_G3xy_aa-2.0E0*5*I_ESP_H5y_G3xy_a-2.0E0*6*I_ESP_H5y_G3xy_a+5*4*I_ESP_F3y_G3xy;
    abcd[iGrid*1890+982] = 4.0E0*I_ESP_K6yz_G3xy_aa-2.0E0*4*I_ESP_H4yz_G3xy_a-2.0E0*5*I_ESP_H4yz_G3xy_a+4*3*I_ESP_F2yz_G3xy;
    abcd[iGrid*1890+983] = 4.0E0*I_ESP_K5y2z_G3xy_aa-2.0E0*3*I_ESP_H3y2z_G3xy_a-2.0E0*4*I_ESP_H3y2z_G3xy_a+3*2*I_ESP_Fy2z_G3xy;
    abcd[iGrid*1890+984] = 4.0E0*I_ESP_K4y3z_G3xy_aa-2.0E0*2*I_ESP_H2y3z_G3xy_a-2.0E0*3*I_ESP_H2y3z_G3xy_a+2*1*I_ESP_F3z_G3xy;
    abcd[iGrid*1890+985] = 4.0E0*I_ESP_K3y4z_G3xy_aa-2.0E0*1*I_ESP_Hy4z_G3xy_a-2.0E0*2*I_ESP_Hy4z_G3xy_a;
    abcd[iGrid*1890+986] = 4.0E0*I_ESP_K2y5z_G3xy_aa-2.0E0*1*I_ESP_H5z_G3xy_a;
    abcd[iGrid*1890+987] = 4.0E0*I_ESP_K5x2y_G3xz_aa-2.0E0*1*I_ESP_H5x_G3xz_a;
    abcd[iGrid*1890+988] = 4.0E0*I_ESP_K4x3y_G3xz_aa-2.0E0*1*I_ESP_H4xy_G3xz_a-2.0E0*2*I_ESP_H4xy_G3xz_a;
    abcd[iGrid*1890+989] = 4.0E0*I_ESP_K4x2yz_G3xz_aa-2.0E0*1*I_ESP_H4xz_G3xz_a;
    abcd[iGrid*1890+990] = 4.0E0*I_ESP_K3x4y_G3xz_aa-2.0E0*2*I_ESP_H3x2y_G3xz_a-2.0E0*3*I_ESP_H3x2y_G3xz_a+2*1*I_ESP_F3x_G3xz;
    abcd[iGrid*1890+991] = 4.0E0*I_ESP_K3x3yz_G3xz_aa-2.0E0*1*I_ESP_H3xyz_G3xz_a-2.0E0*2*I_ESP_H3xyz_G3xz_a;
    abcd[iGrid*1890+992] = 4.0E0*I_ESP_K3x2y2z_G3xz_aa-2.0E0*1*I_ESP_H3x2z_G3xz_a;
    abcd[iGrid*1890+993] = 4.0E0*I_ESP_K2x5y_G3xz_aa-2.0E0*3*I_ESP_H2x3y_G3xz_a-2.0E0*4*I_ESP_H2x3y_G3xz_a+3*2*I_ESP_F2xy_G3xz;
    abcd[iGrid*1890+994] = 4.0E0*I_ESP_K2x4yz_G3xz_aa-2.0E0*2*I_ESP_H2x2yz_G3xz_a-2.0E0*3*I_ESP_H2x2yz_G3xz_a+2*1*I_ESP_F2xz_G3xz;
    abcd[iGrid*1890+995] = 4.0E0*I_ESP_K2x3y2z_G3xz_aa-2.0E0*1*I_ESP_H2xy2z_G3xz_a-2.0E0*2*I_ESP_H2xy2z_G3xz_a;
    abcd[iGrid*1890+996] = 4.0E0*I_ESP_K2x2y3z_G3xz_aa-2.0E0*1*I_ESP_H2x3z_G3xz_a;
    abcd[iGrid*1890+997] = 4.0E0*I_ESP_Kx6y_G3xz_aa-2.0E0*4*I_ESP_Hx4y_G3xz_a-2.0E0*5*I_ESP_Hx4y_G3xz_a+4*3*I_ESP_Fx2y_G3xz;
    abcd[iGrid*1890+998] = 4.0E0*I_ESP_Kx5yz_G3xz_aa-2.0E0*3*I_ESP_Hx3yz_G3xz_a-2.0E0*4*I_ESP_Hx3yz_G3xz_a+3*2*I_ESP_Fxyz_G3xz;
    abcd[iGrid*1890+999] = 4.0E0*I_ESP_Kx4y2z_G3xz_aa-2.0E0*2*I_ESP_Hx2y2z_G3xz_a-2.0E0*3*I_ESP_Hx2y2z_G3xz_a+2*1*I_ESP_Fx2z_G3xz;
    abcd[iGrid*1890+1000] = 4.0E0*I_ESP_Kx3y3z_G3xz_aa-2.0E0*1*I_ESP_Hxy3z_G3xz_a-2.0E0*2*I_ESP_Hxy3z_G3xz_a;
    abcd[iGrid*1890+1001] = 4.0E0*I_ESP_Kx2y4z_G3xz_aa-2.0E0*1*I_ESP_Hx4z_G3xz_a;
    abcd[iGrid*1890+1002] = 4.0E0*I_ESP_K7y_G3xz_aa-2.0E0*5*I_ESP_H5y_G3xz_a-2.0E0*6*I_ESP_H5y_G3xz_a+5*4*I_ESP_F3y_G3xz;
    abcd[iGrid*1890+1003] = 4.0E0*I_ESP_K6yz_G3xz_aa-2.0E0*4*I_ESP_H4yz_G3xz_a-2.0E0*5*I_ESP_H4yz_G3xz_a+4*3*I_ESP_F2yz_G3xz;
    abcd[iGrid*1890+1004] = 4.0E0*I_ESP_K5y2z_G3xz_aa-2.0E0*3*I_ESP_H3y2z_G3xz_a-2.0E0*4*I_ESP_H3y2z_G3xz_a+3*2*I_ESP_Fy2z_G3xz;
    abcd[iGrid*1890+1005] = 4.0E0*I_ESP_K4y3z_G3xz_aa-2.0E0*2*I_ESP_H2y3z_G3xz_a-2.0E0*3*I_ESP_H2y3z_G3xz_a+2*1*I_ESP_F3z_G3xz;
    abcd[iGrid*1890+1006] = 4.0E0*I_ESP_K3y4z_G3xz_aa-2.0E0*1*I_ESP_Hy4z_G3xz_a-2.0E0*2*I_ESP_Hy4z_G3xz_a;
    abcd[iGrid*1890+1007] = 4.0E0*I_ESP_K2y5z_G3xz_aa-2.0E0*1*I_ESP_H5z_G3xz_a;
    abcd[iGrid*1890+1008] = 4.0E0*I_ESP_K5x2y_G2x2y_aa-2.0E0*1*I_ESP_H5x_G2x2y_a;
    abcd[iGrid*1890+1009] = 4.0E0*I_ESP_K4x3y_G2x2y_aa-2.0E0*1*I_ESP_H4xy_G2x2y_a-2.0E0*2*I_ESP_H4xy_G2x2y_a;
    abcd[iGrid*1890+1010] = 4.0E0*I_ESP_K4x2yz_G2x2y_aa-2.0E0*1*I_ESP_H4xz_G2x2y_a;
    abcd[iGrid*1890+1011] = 4.0E0*I_ESP_K3x4y_G2x2y_aa-2.0E0*2*I_ESP_H3x2y_G2x2y_a-2.0E0*3*I_ESP_H3x2y_G2x2y_a+2*1*I_ESP_F3x_G2x2y;
    abcd[iGrid*1890+1012] = 4.0E0*I_ESP_K3x3yz_G2x2y_aa-2.0E0*1*I_ESP_H3xyz_G2x2y_a-2.0E0*2*I_ESP_H3xyz_G2x2y_a;
    abcd[iGrid*1890+1013] = 4.0E0*I_ESP_K3x2y2z_G2x2y_aa-2.0E0*1*I_ESP_H3x2z_G2x2y_a;
    abcd[iGrid*1890+1014] = 4.0E0*I_ESP_K2x5y_G2x2y_aa-2.0E0*3*I_ESP_H2x3y_G2x2y_a-2.0E0*4*I_ESP_H2x3y_G2x2y_a+3*2*I_ESP_F2xy_G2x2y;
    abcd[iGrid*1890+1015] = 4.0E0*I_ESP_K2x4yz_G2x2y_aa-2.0E0*2*I_ESP_H2x2yz_G2x2y_a-2.0E0*3*I_ESP_H2x2yz_G2x2y_a+2*1*I_ESP_F2xz_G2x2y;
    abcd[iGrid*1890+1016] = 4.0E0*I_ESP_K2x3y2z_G2x2y_aa-2.0E0*1*I_ESP_H2xy2z_G2x2y_a-2.0E0*2*I_ESP_H2xy2z_G2x2y_a;
    abcd[iGrid*1890+1017] = 4.0E0*I_ESP_K2x2y3z_G2x2y_aa-2.0E0*1*I_ESP_H2x3z_G2x2y_a;
    abcd[iGrid*1890+1018] = 4.0E0*I_ESP_Kx6y_G2x2y_aa-2.0E0*4*I_ESP_Hx4y_G2x2y_a-2.0E0*5*I_ESP_Hx4y_G2x2y_a+4*3*I_ESP_Fx2y_G2x2y;
    abcd[iGrid*1890+1019] = 4.0E0*I_ESP_Kx5yz_G2x2y_aa-2.0E0*3*I_ESP_Hx3yz_G2x2y_a-2.0E0*4*I_ESP_Hx3yz_G2x2y_a+3*2*I_ESP_Fxyz_G2x2y;
    abcd[iGrid*1890+1020] = 4.0E0*I_ESP_Kx4y2z_G2x2y_aa-2.0E0*2*I_ESP_Hx2y2z_G2x2y_a-2.0E0*3*I_ESP_Hx2y2z_G2x2y_a+2*1*I_ESP_Fx2z_G2x2y;
    abcd[iGrid*1890+1021] = 4.0E0*I_ESP_Kx3y3z_G2x2y_aa-2.0E0*1*I_ESP_Hxy3z_G2x2y_a-2.0E0*2*I_ESP_Hxy3z_G2x2y_a;
    abcd[iGrid*1890+1022] = 4.0E0*I_ESP_Kx2y4z_G2x2y_aa-2.0E0*1*I_ESP_Hx4z_G2x2y_a;
    abcd[iGrid*1890+1023] = 4.0E0*I_ESP_K7y_G2x2y_aa-2.0E0*5*I_ESP_H5y_G2x2y_a-2.0E0*6*I_ESP_H5y_G2x2y_a+5*4*I_ESP_F3y_G2x2y;
    abcd[iGrid*1890+1024] = 4.0E0*I_ESP_K6yz_G2x2y_aa-2.0E0*4*I_ESP_H4yz_G2x2y_a-2.0E0*5*I_ESP_H4yz_G2x2y_a+4*3*I_ESP_F2yz_G2x2y;
    abcd[iGrid*1890+1025] = 4.0E0*I_ESP_K5y2z_G2x2y_aa-2.0E0*3*I_ESP_H3y2z_G2x2y_a-2.0E0*4*I_ESP_H3y2z_G2x2y_a+3*2*I_ESP_Fy2z_G2x2y;
    abcd[iGrid*1890+1026] = 4.0E0*I_ESP_K4y3z_G2x2y_aa-2.0E0*2*I_ESP_H2y3z_G2x2y_a-2.0E0*3*I_ESP_H2y3z_G2x2y_a+2*1*I_ESP_F3z_G2x2y;
    abcd[iGrid*1890+1027] = 4.0E0*I_ESP_K3y4z_G2x2y_aa-2.0E0*1*I_ESP_Hy4z_G2x2y_a-2.0E0*2*I_ESP_Hy4z_G2x2y_a;
    abcd[iGrid*1890+1028] = 4.0E0*I_ESP_K2y5z_G2x2y_aa-2.0E0*1*I_ESP_H5z_G2x2y_a;
    abcd[iGrid*1890+1029] = 4.0E0*I_ESP_K5x2y_G2xyz_aa-2.0E0*1*I_ESP_H5x_G2xyz_a;
    abcd[iGrid*1890+1030] = 4.0E0*I_ESP_K4x3y_G2xyz_aa-2.0E0*1*I_ESP_H4xy_G2xyz_a-2.0E0*2*I_ESP_H4xy_G2xyz_a;
    abcd[iGrid*1890+1031] = 4.0E0*I_ESP_K4x2yz_G2xyz_aa-2.0E0*1*I_ESP_H4xz_G2xyz_a;
    abcd[iGrid*1890+1032] = 4.0E0*I_ESP_K3x4y_G2xyz_aa-2.0E0*2*I_ESP_H3x2y_G2xyz_a-2.0E0*3*I_ESP_H3x2y_G2xyz_a+2*1*I_ESP_F3x_G2xyz;
    abcd[iGrid*1890+1033] = 4.0E0*I_ESP_K3x3yz_G2xyz_aa-2.0E0*1*I_ESP_H3xyz_G2xyz_a-2.0E0*2*I_ESP_H3xyz_G2xyz_a;
    abcd[iGrid*1890+1034] = 4.0E0*I_ESP_K3x2y2z_G2xyz_aa-2.0E0*1*I_ESP_H3x2z_G2xyz_a;
    abcd[iGrid*1890+1035] = 4.0E0*I_ESP_K2x5y_G2xyz_aa-2.0E0*3*I_ESP_H2x3y_G2xyz_a-2.0E0*4*I_ESP_H2x3y_G2xyz_a+3*2*I_ESP_F2xy_G2xyz;
    abcd[iGrid*1890+1036] = 4.0E0*I_ESP_K2x4yz_G2xyz_aa-2.0E0*2*I_ESP_H2x2yz_G2xyz_a-2.0E0*3*I_ESP_H2x2yz_G2xyz_a+2*1*I_ESP_F2xz_G2xyz;
    abcd[iGrid*1890+1037] = 4.0E0*I_ESP_K2x3y2z_G2xyz_aa-2.0E0*1*I_ESP_H2xy2z_G2xyz_a-2.0E0*2*I_ESP_H2xy2z_G2xyz_a;
    abcd[iGrid*1890+1038] = 4.0E0*I_ESP_K2x2y3z_G2xyz_aa-2.0E0*1*I_ESP_H2x3z_G2xyz_a;
    abcd[iGrid*1890+1039] = 4.0E0*I_ESP_Kx6y_G2xyz_aa-2.0E0*4*I_ESP_Hx4y_G2xyz_a-2.0E0*5*I_ESP_Hx4y_G2xyz_a+4*3*I_ESP_Fx2y_G2xyz;
    abcd[iGrid*1890+1040] = 4.0E0*I_ESP_Kx5yz_G2xyz_aa-2.0E0*3*I_ESP_Hx3yz_G2xyz_a-2.0E0*4*I_ESP_Hx3yz_G2xyz_a+3*2*I_ESP_Fxyz_G2xyz;
    abcd[iGrid*1890+1041] = 4.0E0*I_ESP_Kx4y2z_G2xyz_aa-2.0E0*2*I_ESP_Hx2y2z_G2xyz_a-2.0E0*3*I_ESP_Hx2y2z_G2xyz_a+2*1*I_ESP_Fx2z_G2xyz;
    abcd[iGrid*1890+1042] = 4.0E0*I_ESP_Kx3y3z_G2xyz_aa-2.0E0*1*I_ESP_Hxy3z_G2xyz_a-2.0E0*2*I_ESP_Hxy3z_G2xyz_a;
    abcd[iGrid*1890+1043] = 4.0E0*I_ESP_Kx2y4z_G2xyz_aa-2.0E0*1*I_ESP_Hx4z_G2xyz_a;
    abcd[iGrid*1890+1044] = 4.0E0*I_ESP_K7y_G2xyz_aa-2.0E0*5*I_ESP_H5y_G2xyz_a-2.0E0*6*I_ESP_H5y_G2xyz_a+5*4*I_ESP_F3y_G2xyz;
    abcd[iGrid*1890+1045] = 4.0E0*I_ESP_K6yz_G2xyz_aa-2.0E0*4*I_ESP_H4yz_G2xyz_a-2.0E0*5*I_ESP_H4yz_G2xyz_a+4*3*I_ESP_F2yz_G2xyz;
    abcd[iGrid*1890+1046] = 4.0E0*I_ESP_K5y2z_G2xyz_aa-2.0E0*3*I_ESP_H3y2z_G2xyz_a-2.0E0*4*I_ESP_H3y2z_G2xyz_a+3*2*I_ESP_Fy2z_G2xyz;
    abcd[iGrid*1890+1047] = 4.0E0*I_ESP_K4y3z_G2xyz_aa-2.0E0*2*I_ESP_H2y3z_G2xyz_a-2.0E0*3*I_ESP_H2y3z_G2xyz_a+2*1*I_ESP_F3z_G2xyz;
    abcd[iGrid*1890+1048] = 4.0E0*I_ESP_K3y4z_G2xyz_aa-2.0E0*1*I_ESP_Hy4z_G2xyz_a-2.0E0*2*I_ESP_Hy4z_G2xyz_a;
    abcd[iGrid*1890+1049] = 4.0E0*I_ESP_K2y5z_G2xyz_aa-2.0E0*1*I_ESP_H5z_G2xyz_a;
    abcd[iGrid*1890+1050] = 4.0E0*I_ESP_K5x2y_G2x2z_aa-2.0E0*1*I_ESP_H5x_G2x2z_a;
    abcd[iGrid*1890+1051] = 4.0E0*I_ESP_K4x3y_G2x2z_aa-2.0E0*1*I_ESP_H4xy_G2x2z_a-2.0E0*2*I_ESP_H4xy_G2x2z_a;
    abcd[iGrid*1890+1052] = 4.0E0*I_ESP_K4x2yz_G2x2z_aa-2.0E0*1*I_ESP_H4xz_G2x2z_a;
    abcd[iGrid*1890+1053] = 4.0E0*I_ESP_K3x4y_G2x2z_aa-2.0E0*2*I_ESP_H3x2y_G2x2z_a-2.0E0*3*I_ESP_H3x2y_G2x2z_a+2*1*I_ESP_F3x_G2x2z;
    abcd[iGrid*1890+1054] = 4.0E0*I_ESP_K3x3yz_G2x2z_aa-2.0E0*1*I_ESP_H3xyz_G2x2z_a-2.0E0*2*I_ESP_H3xyz_G2x2z_a;
    abcd[iGrid*1890+1055] = 4.0E0*I_ESP_K3x2y2z_G2x2z_aa-2.0E0*1*I_ESP_H3x2z_G2x2z_a;
    abcd[iGrid*1890+1056] = 4.0E0*I_ESP_K2x5y_G2x2z_aa-2.0E0*3*I_ESP_H2x3y_G2x2z_a-2.0E0*4*I_ESP_H2x3y_G2x2z_a+3*2*I_ESP_F2xy_G2x2z;
    abcd[iGrid*1890+1057] = 4.0E0*I_ESP_K2x4yz_G2x2z_aa-2.0E0*2*I_ESP_H2x2yz_G2x2z_a-2.0E0*3*I_ESP_H2x2yz_G2x2z_a+2*1*I_ESP_F2xz_G2x2z;
    abcd[iGrid*1890+1058] = 4.0E0*I_ESP_K2x3y2z_G2x2z_aa-2.0E0*1*I_ESP_H2xy2z_G2x2z_a-2.0E0*2*I_ESP_H2xy2z_G2x2z_a;
    abcd[iGrid*1890+1059] = 4.0E0*I_ESP_K2x2y3z_G2x2z_aa-2.0E0*1*I_ESP_H2x3z_G2x2z_a;
    abcd[iGrid*1890+1060] = 4.0E0*I_ESP_Kx6y_G2x2z_aa-2.0E0*4*I_ESP_Hx4y_G2x2z_a-2.0E0*5*I_ESP_Hx4y_G2x2z_a+4*3*I_ESP_Fx2y_G2x2z;
    abcd[iGrid*1890+1061] = 4.0E0*I_ESP_Kx5yz_G2x2z_aa-2.0E0*3*I_ESP_Hx3yz_G2x2z_a-2.0E0*4*I_ESP_Hx3yz_G2x2z_a+3*2*I_ESP_Fxyz_G2x2z;
    abcd[iGrid*1890+1062] = 4.0E0*I_ESP_Kx4y2z_G2x2z_aa-2.0E0*2*I_ESP_Hx2y2z_G2x2z_a-2.0E0*3*I_ESP_Hx2y2z_G2x2z_a+2*1*I_ESP_Fx2z_G2x2z;
    abcd[iGrid*1890+1063] = 4.0E0*I_ESP_Kx3y3z_G2x2z_aa-2.0E0*1*I_ESP_Hxy3z_G2x2z_a-2.0E0*2*I_ESP_Hxy3z_G2x2z_a;
    abcd[iGrid*1890+1064] = 4.0E0*I_ESP_Kx2y4z_G2x2z_aa-2.0E0*1*I_ESP_Hx4z_G2x2z_a;
    abcd[iGrid*1890+1065] = 4.0E0*I_ESP_K7y_G2x2z_aa-2.0E0*5*I_ESP_H5y_G2x2z_a-2.0E0*6*I_ESP_H5y_G2x2z_a+5*4*I_ESP_F3y_G2x2z;
    abcd[iGrid*1890+1066] = 4.0E0*I_ESP_K6yz_G2x2z_aa-2.0E0*4*I_ESP_H4yz_G2x2z_a-2.0E0*5*I_ESP_H4yz_G2x2z_a+4*3*I_ESP_F2yz_G2x2z;
    abcd[iGrid*1890+1067] = 4.0E0*I_ESP_K5y2z_G2x2z_aa-2.0E0*3*I_ESP_H3y2z_G2x2z_a-2.0E0*4*I_ESP_H3y2z_G2x2z_a+3*2*I_ESP_Fy2z_G2x2z;
    abcd[iGrid*1890+1068] = 4.0E0*I_ESP_K4y3z_G2x2z_aa-2.0E0*2*I_ESP_H2y3z_G2x2z_a-2.0E0*3*I_ESP_H2y3z_G2x2z_a+2*1*I_ESP_F3z_G2x2z;
    abcd[iGrid*1890+1069] = 4.0E0*I_ESP_K3y4z_G2x2z_aa-2.0E0*1*I_ESP_Hy4z_G2x2z_a-2.0E0*2*I_ESP_Hy4z_G2x2z_a;
    abcd[iGrid*1890+1070] = 4.0E0*I_ESP_K2y5z_G2x2z_aa-2.0E0*1*I_ESP_H5z_G2x2z_a;
    abcd[iGrid*1890+1071] = 4.0E0*I_ESP_K5x2y_Gx3y_aa-2.0E0*1*I_ESP_H5x_Gx3y_a;
    abcd[iGrid*1890+1072] = 4.0E0*I_ESP_K4x3y_Gx3y_aa-2.0E0*1*I_ESP_H4xy_Gx3y_a-2.0E0*2*I_ESP_H4xy_Gx3y_a;
    abcd[iGrid*1890+1073] = 4.0E0*I_ESP_K4x2yz_Gx3y_aa-2.0E0*1*I_ESP_H4xz_Gx3y_a;
    abcd[iGrid*1890+1074] = 4.0E0*I_ESP_K3x4y_Gx3y_aa-2.0E0*2*I_ESP_H3x2y_Gx3y_a-2.0E0*3*I_ESP_H3x2y_Gx3y_a+2*1*I_ESP_F3x_Gx3y;
    abcd[iGrid*1890+1075] = 4.0E0*I_ESP_K3x3yz_Gx3y_aa-2.0E0*1*I_ESP_H3xyz_Gx3y_a-2.0E0*2*I_ESP_H3xyz_Gx3y_a;
    abcd[iGrid*1890+1076] = 4.0E0*I_ESP_K3x2y2z_Gx3y_aa-2.0E0*1*I_ESP_H3x2z_Gx3y_a;
    abcd[iGrid*1890+1077] = 4.0E0*I_ESP_K2x5y_Gx3y_aa-2.0E0*3*I_ESP_H2x3y_Gx3y_a-2.0E0*4*I_ESP_H2x3y_Gx3y_a+3*2*I_ESP_F2xy_Gx3y;
    abcd[iGrid*1890+1078] = 4.0E0*I_ESP_K2x4yz_Gx3y_aa-2.0E0*2*I_ESP_H2x2yz_Gx3y_a-2.0E0*3*I_ESP_H2x2yz_Gx3y_a+2*1*I_ESP_F2xz_Gx3y;
    abcd[iGrid*1890+1079] = 4.0E0*I_ESP_K2x3y2z_Gx3y_aa-2.0E0*1*I_ESP_H2xy2z_Gx3y_a-2.0E0*2*I_ESP_H2xy2z_Gx3y_a;
    abcd[iGrid*1890+1080] = 4.0E0*I_ESP_K2x2y3z_Gx3y_aa-2.0E0*1*I_ESP_H2x3z_Gx3y_a;
    abcd[iGrid*1890+1081] = 4.0E0*I_ESP_Kx6y_Gx3y_aa-2.0E0*4*I_ESP_Hx4y_Gx3y_a-2.0E0*5*I_ESP_Hx4y_Gx3y_a+4*3*I_ESP_Fx2y_Gx3y;
    abcd[iGrid*1890+1082] = 4.0E0*I_ESP_Kx5yz_Gx3y_aa-2.0E0*3*I_ESP_Hx3yz_Gx3y_a-2.0E0*4*I_ESP_Hx3yz_Gx3y_a+3*2*I_ESP_Fxyz_Gx3y;
    abcd[iGrid*1890+1083] = 4.0E0*I_ESP_Kx4y2z_Gx3y_aa-2.0E0*2*I_ESP_Hx2y2z_Gx3y_a-2.0E0*3*I_ESP_Hx2y2z_Gx3y_a+2*1*I_ESP_Fx2z_Gx3y;
    abcd[iGrid*1890+1084] = 4.0E0*I_ESP_Kx3y3z_Gx3y_aa-2.0E0*1*I_ESP_Hxy3z_Gx3y_a-2.0E0*2*I_ESP_Hxy3z_Gx3y_a;
    abcd[iGrid*1890+1085] = 4.0E0*I_ESP_Kx2y4z_Gx3y_aa-2.0E0*1*I_ESP_Hx4z_Gx3y_a;
    abcd[iGrid*1890+1086] = 4.0E0*I_ESP_K7y_Gx3y_aa-2.0E0*5*I_ESP_H5y_Gx3y_a-2.0E0*6*I_ESP_H5y_Gx3y_a+5*4*I_ESP_F3y_Gx3y;
    abcd[iGrid*1890+1087] = 4.0E0*I_ESP_K6yz_Gx3y_aa-2.0E0*4*I_ESP_H4yz_Gx3y_a-2.0E0*5*I_ESP_H4yz_Gx3y_a+4*3*I_ESP_F2yz_Gx3y;
    abcd[iGrid*1890+1088] = 4.0E0*I_ESP_K5y2z_Gx3y_aa-2.0E0*3*I_ESP_H3y2z_Gx3y_a-2.0E0*4*I_ESP_H3y2z_Gx3y_a+3*2*I_ESP_Fy2z_Gx3y;
    abcd[iGrid*1890+1089] = 4.0E0*I_ESP_K4y3z_Gx3y_aa-2.0E0*2*I_ESP_H2y3z_Gx3y_a-2.0E0*3*I_ESP_H2y3z_Gx3y_a+2*1*I_ESP_F3z_Gx3y;
    abcd[iGrid*1890+1090] = 4.0E0*I_ESP_K3y4z_Gx3y_aa-2.0E0*1*I_ESP_Hy4z_Gx3y_a-2.0E0*2*I_ESP_Hy4z_Gx3y_a;
    abcd[iGrid*1890+1091] = 4.0E0*I_ESP_K2y5z_Gx3y_aa-2.0E0*1*I_ESP_H5z_Gx3y_a;
    abcd[iGrid*1890+1092] = 4.0E0*I_ESP_K5x2y_Gx2yz_aa-2.0E0*1*I_ESP_H5x_Gx2yz_a;
    abcd[iGrid*1890+1093] = 4.0E0*I_ESP_K4x3y_Gx2yz_aa-2.0E0*1*I_ESP_H4xy_Gx2yz_a-2.0E0*2*I_ESP_H4xy_Gx2yz_a;
    abcd[iGrid*1890+1094] = 4.0E0*I_ESP_K4x2yz_Gx2yz_aa-2.0E0*1*I_ESP_H4xz_Gx2yz_a;
    abcd[iGrid*1890+1095] = 4.0E0*I_ESP_K3x4y_Gx2yz_aa-2.0E0*2*I_ESP_H3x2y_Gx2yz_a-2.0E0*3*I_ESP_H3x2y_Gx2yz_a+2*1*I_ESP_F3x_Gx2yz;
    abcd[iGrid*1890+1096] = 4.0E0*I_ESP_K3x3yz_Gx2yz_aa-2.0E0*1*I_ESP_H3xyz_Gx2yz_a-2.0E0*2*I_ESP_H3xyz_Gx2yz_a;
    abcd[iGrid*1890+1097] = 4.0E0*I_ESP_K3x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_H3x2z_Gx2yz_a;
    abcd[iGrid*1890+1098] = 4.0E0*I_ESP_K2x5y_Gx2yz_aa-2.0E0*3*I_ESP_H2x3y_Gx2yz_a-2.0E0*4*I_ESP_H2x3y_Gx2yz_a+3*2*I_ESP_F2xy_Gx2yz;
    abcd[iGrid*1890+1099] = 4.0E0*I_ESP_K2x4yz_Gx2yz_aa-2.0E0*2*I_ESP_H2x2yz_Gx2yz_a-2.0E0*3*I_ESP_H2x2yz_Gx2yz_a+2*1*I_ESP_F2xz_Gx2yz;
    abcd[iGrid*1890+1100] = 4.0E0*I_ESP_K2x3y2z_Gx2yz_aa-2.0E0*1*I_ESP_H2xy2z_Gx2yz_a-2.0E0*2*I_ESP_H2xy2z_Gx2yz_a;
    abcd[iGrid*1890+1101] = 4.0E0*I_ESP_K2x2y3z_Gx2yz_aa-2.0E0*1*I_ESP_H2x3z_Gx2yz_a;
    abcd[iGrid*1890+1102] = 4.0E0*I_ESP_Kx6y_Gx2yz_aa-2.0E0*4*I_ESP_Hx4y_Gx2yz_a-2.0E0*5*I_ESP_Hx4y_Gx2yz_a+4*3*I_ESP_Fx2y_Gx2yz;
    abcd[iGrid*1890+1103] = 4.0E0*I_ESP_Kx5yz_Gx2yz_aa-2.0E0*3*I_ESP_Hx3yz_Gx2yz_a-2.0E0*4*I_ESP_Hx3yz_Gx2yz_a+3*2*I_ESP_Fxyz_Gx2yz;
    abcd[iGrid*1890+1104] = 4.0E0*I_ESP_Kx4y2z_Gx2yz_aa-2.0E0*2*I_ESP_Hx2y2z_Gx2yz_a-2.0E0*3*I_ESP_Hx2y2z_Gx2yz_a+2*1*I_ESP_Fx2z_Gx2yz;
    abcd[iGrid*1890+1105] = 4.0E0*I_ESP_Kx3y3z_Gx2yz_aa-2.0E0*1*I_ESP_Hxy3z_Gx2yz_a-2.0E0*2*I_ESP_Hxy3z_Gx2yz_a;
    abcd[iGrid*1890+1106] = 4.0E0*I_ESP_Kx2y4z_Gx2yz_aa-2.0E0*1*I_ESP_Hx4z_Gx2yz_a;
    abcd[iGrid*1890+1107] = 4.0E0*I_ESP_K7y_Gx2yz_aa-2.0E0*5*I_ESP_H5y_Gx2yz_a-2.0E0*6*I_ESP_H5y_Gx2yz_a+5*4*I_ESP_F3y_Gx2yz;
    abcd[iGrid*1890+1108] = 4.0E0*I_ESP_K6yz_Gx2yz_aa-2.0E0*4*I_ESP_H4yz_Gx2yz_a-2.0E0*5*I_ESP_H4yz_Gx2yz_a+4*3*I_ESP_F2yz_Gx2yz;
    abcd[iGrid*1890+1109] = 4.0E0*I_ESP_K5y2z_Gx2yz_aa-2.0E0*3*I_ESP_H3y2z_Gx2yz_a-2.0E0*4*I_ESP_H3y2z_Gx2yz_a+3*2*I_ESP_Fy2z_Gx2yz;
    abcd[iGrid*1890+1110] = 4.0E0*I_ESP_K4y3z_Gx2yz_aa-2.0E0*2*I_ESP_H2y3z_Gx2yz_a-2.0E0*3*I_ESP_H2y3z_Gx2yz_a+2*1*I_ESP_F3z_Gx2yz;
    abcd[iGrid*1890+1111] = 4.0E0*I_ESP_K3y4z_Gx2yz_aa-2.0E0*1*I_ESP_Hy4z_Gx2yz_a-2.0E0*2*I_ESP_Hy4z_Gx2yz_a;
    abcd[iGrid*1890+1112] = 4.0E0*I_ESP_K2y5z_Gx2yz_aa-2.0E0*1*I_ESP_H5z_Gx2yz_a;
    abcd[iGrid*1890+1113] = 4.0E0*I_ESP_K5x2y_Gxy2z_aa-2.0E0*1*I_ESP_H5x_Gxy2z_a;
    abcd[iGrid*1890+1114] = 4.0E0*I_ESP_K4x3y_Gxy2z_aa-2.0E0*1*I_ESP_H4xy_Gxy2z_a-2.0E0*2*I_ESP_H4xy_Gxy2z_a;
    abcd[iGrid*1890+1115] = 4.0E0*I_ESP_K4x2yz_Gxy2z_aa-2.0E0*1*I_ESP_H4xz_Gxy2z_a;
    abcd[iGrid*1890+1116] = 4.0E0*I_ESP_K3x4y_Gxy2z_aa-2.0E0*2*I_ESP_H3x2y_Gxy2z_a-2.0E0*3*I_ESP_H3x2y_Gxy2z_a+2*1*I_ESP_F3x_Gxy2z;
    abcd[iGrid*1890+1117] = 4.0E0*I_ESP_K3x3yz_Gxy2z_aa-2.0E0*1*I_ESP_H3xyz_Gxy2z_a-2.0E0*2*I_ESP_H3xyz_Gxy2z_a;
    abcd[iGrid*1890+1118] = 4.0E0*I_ESP_K3x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_H3x2z_Gxy2z_a;
    abcd[iGrid*1890+1119] = 4.0E0*I_ESP_K2x5y_Gxy2z_aa-2.0E0*3*I_ESP_H2x3y_Gxy2z_a-2.0E0*4*I_ESP_H2x3y_Gxy2z_a+3*2*I_ESP_F2xy_Gxy2z;
    abcd[iGrid*1890+1120] = 4.0E0*I_ESP_K2x4yz_Gxy2z_aa-2.0E0*2*I_ESP_H2x2yz_Gxy2z_a-2.0E0*3*I_ESP_H2x2yz_Gxy2z_a+2*1*I_ESP_F2xz_Gxy2z;
    abcd[iGrid*1890+1121] = 4.0E0*I_ESP_K2x3y2z_Gxy2z_aa-2.0E0*1*I_ESP_H2xy2z_Gxy2z_a-2.0E0*2*I_ESP_H2xy2z_Gxy2z_a;
    abcd[iGrid*1890+1122] = 4.0E0*I_ESP_K2x2y3z_Gxy2z_aa-2.0E0*1*I_ESP_H2x3z_Gxy2z_a;
    abcd[iGrid*1890+1123] = 4.0E0*I_ESP_Kx6y_Gxy2z_aa-2.0E0*4*I_ESP_Hx4y_Gxy2z_a-2.0E0*5*I_ESP_Hx4y_Gxy2z_a+4*3*I_ESP_Fx2y_Gxy2z;
    abcd[iGrid*1890+1124] = 4.0E0*I_ESP_Kx5yz_Gxy2z_aa-2.0E0*3*I_ESP_Hx3yz_Gxy2z_a-2.0E0*4*I_ESP_Hx3yz_Gxy2z_a+3*2*I_ESP_Fxyz_Gxy2z;
    abcd[iGrid*1890+1125] = 4.0E0*I_ESP_Kx4y2z_Gxy2z_aa-2.0E0*2*I_ESP_Hx2y2z_Gxy2z_a-2.0E0*3*I_ESP_Hx2y2z_Gxy2z_a+2*1*I_ESP_Fx2z_Gxy2z;
    abcd[iGrid*1890+1126] = 4.0E0*I_ESP_Kx3y3z_Gxy2z_aa-2.0E0*1*I_ESP_Hxy3z_Gxy2z_a-2.0E0*2*I_ESP_Hxy3z_Gxy2z_a;
    abcd[iGrid*1890+1127] = 4.0E0*I_ESP_Kx2y4z_Gxy2z_aa-2.0E0*1*I_ESP_Hx4z_Gxy2z_a;
    abcd[iGrid*1890+1128] = 4.0E0*I_ESP_K7y_Gxy2z_aa-2.0E0*5*I_ESP_H5y_Gxy2z_a-2.0E0*6*I_ESP_H5y_Gxy2z_a+5*4*I_ESP_F3y_Gxy2z;
    abcd[iGrid*1890+1129] = 4.0E0*I_ESP_K6yz_Gxy2z_aa-2.0E0*4*I_ESP_H4yz_Gxy2z_a-2.0E0*5*I_ESP_H4yz_Gxy2z_a+4*3*I_ESP_F2yz_Gxy2z;
    abcd[iGrid*1890+1130] = 4.0E0*I_ESP_K5y2z_Gxy2z_aa-2.0E0*3*I_ESP_H3y2z_Gxy2z_a-2.0E0*4*I_ESP_H3y2z_Gxy2z_a+3*2*I_ESP_Fy2z_Gxy2z;
    abcd[iGrid*1890+1131] = 4.0E0*I_ESP_K4y3z_Gxy2z_aa-2.0E0*2*I_ESP_H2y3z_Gxy2z_a-2.0E0*3*I_ESP_H2y3z_Gxy2z_a+2*1*I_ESP_F3z_Gxy2z;
    abcd[iGrid*1890+1132] = 4.0E0*I_ESP_K3y4z_Gxy2z_aa-2.0E0*1*I_ESP_Hy4z_Gxy2z_a-2.0E0*2*I_ESP_Hy4z_Gxy2z_a;
    abcd[iGrid*1890+1133] = 4.0E0*I_ESP_K2y5z_Gxy2z_aa-2.0E0*1*I_ESP_H5z_Gxy2z_a;
    abcd[iGrid*1890+1134] = 4.0E0*I_ESP_K5x2y_Gx3z_aa-2.0E0*1*I_ESP_H5x_Gx3z_a;
    abcd[iGrid*1890+1135] = 4.0E0*I_ESP_K4x3y_Gx3z_aa-2.0E0*1*I_ESP_H4xy_Gx3z_a-2.0E0*2*I_ESP_H4xy_Gx3z_a;
    abcd[iGrid*1890+1136] = 4.0E0*I_ESP_K4x2yz_Gx3z_aa-2.0E0*1*I_ESP_H4xz_Gx3z_a;
    abcd[iGrid*1890+1137] = 4.0E0*I_ESP_K3x4y_Gx3z_aa-2.0E0*2*I_ESP_H3x2y_Gx3z_a-2.0E0*3*I_ESP_H3x2y_Gx3z_a+2*1*I_ESP_F3x_Gx3z;
    abcd[iGrid*1890+1138] = 4.0E0*I_ESP_K3x3yz_Gx3z_aa-2.0E0*1*I_ESP_H3xyz_Gx3z_a-2.0E0*2*I_ESP_H3xyz_Gx3z_a;
    abcd[iGrid*1890+1139] = 4.0E0*I_ESP_K3x2y2z_Gx3z_aa-2.0E0*1*I_ESP_H3x2z_Gx3z_a;
    abcd[iGrid*1890+1140] = 4.0E0*I_ESP_K2x5y_Gx3z_aa-2.0E0*3*I_ESP_H2x3y_Gx3z_a-2.0E0*4*I_ESP_H2x3y_Gx3z_a+3*2*I_ESP_F2xy_Gx3z;
    abcd[iGrid*1890+1141] = 4.0E0*I_ESP_K2x4yz_Gx3z_aa-2.0E0*2*I_ESP_H2x2yz_Gx3z_a-2.0E0*3*I_ESP_H2x2yz_Gx3z_a+2*1*I_ESP_F2xz_Gx3z;
    abcd[iGrid*1890+1142] = 4.0E0*I_ESP_K2x3y2z_Gx3z_aa-2.0E0*1*I_ESP_H2xy2z_Gx3z_a-2.0E0*2*I_ESP_H2xy2z_Gx3z_a;
    abcd[iGrid*1890+1143] = 4.0E0*I_ESP_K2x2y3z_Gx3z_aa-2.0E0*1*I_ESP_H2x3z_Gx3z_a;
    abcd[iGrid*1890+1144] = 4.0E0*I_ESP_Kx6y_Gx3z_aa-2.0E0*4*I_ESP_Hx4y_Gx3z_a-2.0E0*5*I_ESP_Hx4y_Gx3z_a+4*3*I_ESP_Fx2y_Gx3z;
    abcd[iGrid*1890+1145] = 4.0E0*I_ESP_Kx5yz_Gx3z_aa-2.0E0*3*I_ESP_Hx3yz_Gx3z_a-2.0E0*4*I_ESP_Hx3yz_Gx3z_a+3*2*I_ESP_Fxyz_Gx3z;
    abcd[iGrid*1890+1146] = 4.0E0*I_ESP_Kx4y2z_Gx3z_aa-2.0E0*2*I_ESP_Hx2y2z_Gx3z_a-2.0E0*3*I_ESP_Hx2y2z_Gx3z_a+2*1*I_ESP_Fx2z_Gx3z;
    abcd[iGrid*1890+1147] = 4.0E0*I_ESP_Kx3y3z_Gx3z_aa-2.0E0*1*I_ESP_Hxy3z_Gx3z_a-2.0E0*2*I_ESP_Hxy3z_Gx3z_a;
    abcd[iGrid*1890+1148] = 4.0E0*I_ESP_Kx2y4z_Gx3z_aa-2.0E0*1*I_ESP_Hx4z_Gx3z_a;
    abcd[iGrid*1890+1149] = 4.0E0*I_ESP_K7y_Gx3z_aa-2.0E0*5*I_ESP_H5y_Gx3z_a-2.0E0*6*I_ESP_H5y_Gx3z_a+5*4*I_ESP_F3y_Gx3z;
    abcd[iGrid*1890+1150] = 4.0E0*I_ESP_K6yz_Gx3z_aa-2.0E0*4*I_ESP_H4yz_Gx3z_a-2.0E0*5*I_ESP_H4yz_Gx3z_a+4*3*I_ESP_F2yz_Gx3z;
    abcd[iGrid*1890+1151] = 4.0E0*I_ESP_K5y2z_Gx3z_aa-2.0E0*3*I_ESP_H3y2z_Gx3z_a-2.0E0*4*I_ESP_H3y2z_Gx3z_a+3*2*I_ESP_Fy2z_Gx3z;
    abcd[iGrid*1890+1152] = 4.0E0*I_ESP_K4y3z_Gx3z_aa-2.0E0*2*I_ESP_H2y3z_Gx3z_a-2.0E0*3*I_ESP_H2y3z_Gx3z_a+2*1*I_ESP_F3z_Gx3z;
    abcd[iGrid*1890+1153] = 4.0E0*I_ESP_K3y4z_Gx3z_aa-2.0E0*1*I_ESP_Hy4z_Gx3z_a-2.0E0*2*I_ESP_Hy4z_Gx3z_a;
    abcd[iGrid*1890+1154] = 4.0E0*I_ESP_K2y5z_Gx3z_aa-2.0E0*1*I_ESP_H5z_Gx3z_a;
    abcd[iGrid*1890+1155] = 4.0E0*I_ESP_K5x2y_G4y_aa-2.0E0*1*I_ESP_H5x_G4y_a;
    abcd[iGrid*1890+1156] = 4.0E0*I_ESP_K4x3y_G4y_aa-2.0E0*1*I_ESP_H4xy_G4y_a-2.0E0*2*I_ESP_H4xy_G4y_a;
    abcd[iGrid*1890+1157] = 4.0E0*I_ESP_K4x2yz_G4y_aa-2.0E0*1*I_ESP_H4xz_G4y_a;
    abcd[iGrid*1890+1158] = 4.0E0*I_ESP_K3x4y_G4y_aa-2.0E0*2*I_ESP_H3x2y_G4y_a-2.0E0*3*I_ESP_H3x2y_G4y_a+2*1*I_ESP_F3x_G4y;
    abcd[iGrid*1890+1159] = 4.0E0*I_ESP_K3x3yz_G4y_aa-2.0E0*1*I_ESP_H3xyz_G4y_a-2.0E0*2*I_ESP_H3xyz_G4y_a;
    abcd[iGrid*1890+1160] = 4.0E0*I_ESP_K3x2y2z_G4y_aa-2.0E0*1*I_ESP_H3x2z_G4y_a;
    abcd[iGrid*1890+1161] = 4.0E0*I_ESP_K2x5y_G4y_aa-2.0E0*3*I_ESP_H2x3y_G4y_a-2.0E0*4*I_ESP_H2x3y_G4y_a+3*2*I_ESP_F2xy_G4y;
    abcd[iGrid*1890+1162] = 4.0E0*I_ESP_K2x4yz_G4y_aa-2.0E0*2*I_ESP_H2x2yz_G4y_a-2.0E0*3*I_ESP_H2x2yz_G4y_a+2*1*I_ESP_F2xz_G4y;
    abcd[iGrid*1890+1163] = 4.0E0*I_ESP_K2x3y2z_G4y_aa-2.0E0*1*I_ESP_H2xy2z_G4y_a-2.0E0*2*I_ESP_H2xy2z_G4y_a;
    abcd[iGrid*1890+1164] = 4.0E0*I_ESP_K2x2y3z_G4y_aa-2.0E0*1*I_ESP_H2x3z_G4y_a;
    abcd[iGrid*1890+1165] = 4.0E0*I_ESP_Kx6y_G4y_aa-2.0E0*4*I_ESP_Hx4y_G4y_a-2.0E0*5*I_ESP_Hx4y_G4y_a+4*3*I_ESP_Fx2y_G4y;
    abcd[iGrid*1890+1166] = 4.0E0*I_ESP_Kx5yz_G4y_aa-2.0E0*3*I_ESP_Hx3yz_G4y_a-2.0E0*4*I_ESP_Hx3yz_G4y_a+3*2*I_ESP_Fxyz_G4y;
    abcd[iGrid*1890+1167] = 4.0E0*I_ESP_Kx4y2z_G4y_aa-2.0E0*2*I_ESP_Hx2y2z_G4y_a-2.0E0*3*I_ESP_Hx2y2z_G4y_a+2*1*I_ESP_Fx2z_G4y;
    abcd[iGrid*1890+1168] = 4.0E0*I_ESP_Kx3y3z_G4y_aa-2.0E0*1*I_ESP_Hxy3z_G4y_a-2.0E0*2*I_ESP_Hxy3z_G4y_a;
    abcd[iGrid*1890+1169] = 4.0E0*I_ESP_Kx2y4z_G4y_aa-2.0E0*1*I_ESP_Hx4z_G4y_a;
    abcd[iGrid*1890+1170] = 4.0E0*I_ESP_K7y_G4y_aa-2.0E0*5*I_ESP_H5y_G4y_a-2.0E0*6*I_ESP_H5y_G4y_a+5*4*I_ESP_F3y_G4y;
    abcd[iGrid*1890+1171] = 4.0E0*I_ESP_K6yz_G4y_aa-2.0E0*4*I_ESP_H4yz_G4y_a-2.0E0*5*I_ESP_H4yz_G4y_a+4*3*I_ESP_F2yz_G4y;
    abcd[iGrid*1890+1172] = 4.0E0*I_ESP_K5y2z_G4y_aa-2.0E0*3*I_ESP_H3y2z_G4y_a-2.0E0*4*I_ESP_H3y2z_G4y_a+3*2*I_ESP_Fy2z_G4y;
    abcd[iGrid*1890+1173] = 4.0E0*I_ESP_K4y3z_G4y_aa-2.0E0*2*I_ESP_H2y3z_G4y_a-2.0E0*3*I_ESP_H2y3z_G4y_a+2*1*I_ESP_F3z_G4y;
    abcd[iGrid*1890+1174] = 4.0E0*I_ESP_K3y4z_G4y_aa-2.0E0*1*I_ESP_Hy4z_G4y_a-2.0E0*2*I_ESP_Hy4z_G4y_a;
    abcd[iGrid*1890+1175] = 4.0E0*I_ESP_K2y5z_G4y_aa-2.0E0*1*I_ESP_H5z_G4y_a;
    abcd[iGrid*1890+1176] = 4.0E0*I_ESP_K5x2y_G3yz_aa-2.0E0*1*I_ESP_H5x_G3yz_a;
    abcd[iGrid*1890+1177] = 4.0E0*I_ESP_K4x3y_G3yz_aa-2.0E0*1*I_ESP_H4xy_G3yz_a-2.0E0*2*I_ESP_H4xy_G3yz_a;
    abcd[iGrid*1890+1178] = 4.0E0*I_ESP_K4x2yz_G3yz_aa-2.0E0*1*I_ESP_H4xz_G3yz_a;
    abcd[iGrid*1890+1179] = 4.0E0*I_ESP_K3x4y_G3yz_aa-2.0E0*2*I_ESP_H3x2y_G3yz_a-2.0E0*3*I_ESP_H3x2y_G3yz_a+2*1*I_ESP_F3x_G3yz;
    abcd[iGrid*1890+1180] = 4.0E0*I_ESP_K3x3yz_G3yz_aa-2.0E0*1*I_ESP_H3xyz_G3yz_a-2.0E0*2*I_ESP_H3xyz_G3yz_a;
    abcd[iGrid*1890+1181] = 4.0E0*I_ESP_K3x2y2z_G3yz_aa-2.0E0*1*I_ESP_H3x2z_G3yz_a;
    abcd[iGrid*1890+1182] = 4.0E0*I_ESP_K2x5y_G3yz_aa-2.0E0*3*I_ESP_H2x3y_G3yz_a-2.0E0*4*I_ESP_H2x3y_G3yz_a+3*2*I_ESP_F2xy_G3yz;
    abcd[iGrid*1890+1183] = 4.0E0*I_ESP_K2x4yz_G3yz_aa-2.0E0*2*I_ESP_H2x2yz_G3yz_a-2.0E0*3*I_ESP_H2x2yz_G3yz_a+2*1*I_ESP_F2xz_G3yz;
    abcd[iGrid*1890+1184] = 4.0E0*I_ESP_K2x3y2z_G3yz_aa-2.0E0*1*I_ESP_H2xy2z_G3yz_a-2.0E0*2*I_ESP_H2xy2z_G3yz_a;
    abcd[iGrid*1890+1185] = 4.0E0*I_ESP_K2x2y3z_G3yz_aa-2.0E0*1*I_ESP_H2x3z_G3yz_a;
    abcd[iGrid*1890+1186] = 4.0E0*I_ESP_Kx6y_G3yz_aa-2.0E0*4*I_ESP_Hx4y_G3yz_a-2.0E0*5*I_ESP_Hx4y_G3yz_a+4*3*I_ESP_Fx2y_G3yz;
    abcd[iGrid*1890+1187] = 4.0E0*I_ESP_Kx5yz_G3yz_aa-2.0E0*3*I_ESP_Hx3yz_G3yz_a-2.0E0*4*I_ESP_Hx3yz_G3yz_a+3*2*I_ESP_Fxyz_G3yz;
    abcd[iGrid*1890+1188] = 4.0E0*I_ESP_Kx4y2z_G3yz_aa-2.0E0*2*I_ESP_Hx2y2z_G3yz_a-2.0E0*3*I_ESP_Hx2y2z_G3yz_a+2*1*I_ESP_Fx2z_G3yz;
    abcd[iGrid*1890+1189] = 4.0E0*I_ESP_Kx3y3z_G3yz_aa-2.0E0*1*I_ESP_Hxy3z_G3yz_a-2.0E0*2*I_ESP_Hxy3z_G3yz_a;
    abcd[iGrid*1890+1190] = 4.0E0*I_ESP_Kx2y4z_G3yz_aa-2.0E0*1*I_ESP_Hx4z_G3yz_a;
    abcd[iGrid*1890+1191] = 4.0E0*I_ESP_K7y_G3yz_aa-2.0E0*5*I_ESP_H5y_G3yz_a-2.0E0*6*I_ESP_H5y_G3yz_a+5*4*I_ESP_F3y_G3yz;
    abcd[iGrid*1890+1192] = 4.0E0*I_ESP_K6yz_G3yz_aa-2.0E0*4*I_ESP_H4yz_G3yz_a-2.0E0*5*I_ESP_H4yz_G3yz_a+4*3*I_ESP_F2yz_G3yz;
    abcd[iGrid*1890+1193] = 4.0E0*I_ESP_K5y2z_G3yz_aa-2.0E0*3*I_ESP_H3y2z_G3yz_a-2.0E0*4*I_ESP_H3y2z_G3yz_a+3*2*I_ESP_Fy2z_G3yz;
    abcd[iGrid*1890+1194] = 4.0E0*I_ESP_K4y3z_G3yz_aa-2.0E0*2*I_ESP_H2y3z_G3yz_a-2.0E0*3*I_ESP_H2y3z_G3yz_a+2*1*I_ESP_F3z_G3yz;
    abcd[iGrid*1890+1195] = 4.0E0*I_ESP_K3y4z_G3yz_aa-2.0E0*1*I_ESP_Hy4z_G3yz_a-2.0E0*2*I_ESP_Hy4z_G3yz_a;
    abcd[iGrid*1890+1196] = 4.0E0*I_ESP_K2y5z_G3yz_aa-2.0E0*1*I_ESP_H5z_G3yz_a;
    abcd[iGrid*1890+1197] = 4.0E0*I_ESP_K5x2y_G2y2z_aa-2.0E0*1*I_ESP_H5x_G2y2z_a;
    abcd[iGrid*1890+1198] = 4.0E0*I_ESP_K4x3y_G2y2z_aa-2.0E0*1*I_ESP_H4xy_G2y2z_a-2.0E0*2*I_ESP_H4xy_G2y2z_a;
    abcd[iGrid*1890+1199] = 4.0E0*I_ESP_K4x2yz_G2y2z_aa-2.0E0*1*I_ESP_H4xz_G2y2z_a;
    abcd[iGrid*1890+1200] = 4.0E0*I_ESP_K3x4y_G2y2z_aa-2.0E0*2*I_ESP_H3x2y_G2y2z_a-2.0E0*3*I_ESP_H3x2y_G2y2z_a+2*1*I_ESP_F3x_G2y2z;
    abcd[iGrid*1890+1201] = 4.0E0*I_ESP_K3x3yz_G2y2z_aa-2.0E0*1*I_ESP_H3xyz_G2y2z_a-2.0E0*2*I_ESP_H3xyz_G2y2z_a;
    abcd[iGrid*1890+1202] = 4.0E0*I_ESP_K3x2y2z_G2y2z_aa-2.0E0*1*I_ESP_H3x2z_G2y2z_a;
    abcd[iGrid*1890+1203] = 4.0E0*I_ESP_K2x5y_G2y2z_aa-2.0E0*3*I_ESP_H2x3y_G2y2z_a-2.0E0*4*I_ESP_H2x3y_G2y2z_a+3*2*I_ESP_F2xy_G2y2z;
    abcd[iGrid*1890+1204] = 4.0E0*I_ESP_K2x4yz_G2y2z_aa-2.0E0*2*I_ESP_H2x2yz_G2y2z_a-2.0E0*3*I_ESP_H2x2yz_G2y2z_a+2*1*I_ESP_F2xz_G2y2z;
    abcd[iGrid*1890+1205] = 4.0E0*I_ESP_K2x3y2z_G2y2z_aa-2.0E0*1*I_ESP_H2xy2z_G2y2z_a-2.0E0*2*I_ESP_H2xy2z_G2y2z_a;
    abcd[iGrid*1890+1206] = 4.0E0*I_ESP_K2x2y3z_G2y2z_aa-2.0E0*1*I_ESP_H2x3z_G2y2z_a;
    abcd[iGrid*1890+1207] = 4.0E0*I_ESP_Kx6y_G2y2z_aa-2.0E0*4*I_ESP_Hx4y_G2y2z_a-2.0E0*5*I_ESP_Hx4y_G2y2z_a+4*3*I_ESP_Fx2y_G2y2z;
    abcd[iGrid*1890+1208] = 4.0E0*I_ESP_Kx5yz_G2y2z_aa-2.0E0*3*I_ESP_Hx3yz_G2y2z_a-2.0E0*4*I_ESP_Hx3yz_G2y2z_a+3*2*I_ESP_Fxyz_G2y2z;
    abcd[iGrid*1890+1209] = 4.0E0*I_ESP_Kx4y2z_G2y2z_aa-2.0E0*2*I_ESP_Hx2y2z_G2y2z_a-2.0E0*3*I_ESP_Hx2y2z_G2y2z_a+2*1*I_ESP_Fx2z_G2y2z;
    abcd[iGrid*1890+1210] = 4.0E0*I_ESP_Kx3y3z_G2y2z_aa-2.0E0*1*I_ESP_Hxy3z_G2y2z_a-2.0E0*2*I_ESP_Hxy3z_G2y2z_a;
    abcd[iGrid*1890+1211] = 4.0E0*I_ESP_Kx2y4z_G2y2z_aa-2.0E0*1*I_ESP_Hx4z_G2y2z_a;
    abcd[iGrid*1890+1212] = 4.0E0*I_ESP_K7y_G2y2z_aa-2.0E0*5*I_ESP_H5y_G2y2z_a-2.0E0*6*I_ESP_H5y_G2y2z_a+5*4*I_ESP_F3y_G2y2z;
    abcd[iGrid*1890+1213] = 4.0E0*I_ESP_K6yz_G2y2z_aa-2.0E0*4*I_ESP_H4yz_G2y2z_a-2.0E0*5*I_ESP_H4yz_G2y2z_a+4*3*I_ESP_F2yz_G2y2z;
    abcd[iGrid*1890+1214] = 4.0E0*I_ESP_K5y2z_G2y2z_aa-2.0E0*3*I_ESP_H3y2z_G2y2z_a-2.0E0*4*I_ESP_H3y2z_G2y2z_a+3*2*I_ESP_Fy2z_G2y2z;
    abcd[iGrid*1890+1215] = 4.0E0*I_ESP_K4y3z_G2y2z_aa-2.0E0*2*I_ESP_H2y3z_G2y2z_a-2.0E0*3*I_ESP_H2y3z_G2y2z_a+2*1*I_ESP_F3z_G2y2z;
    abcd[iGrid*1890+1216] = 4.0E0*I_ESP_K3y4z_G2y2z_aa-2.0E0*1*I_ESP_Hy4z_G2y2z_a-2.0E0*2*I_ESP_Hy4z_G2y2z_a;
    abcd[iGrid*1890+1217] = 4.0E0*I_ESP_K2y5z_G2y2z_aa-2.0E0*1*I_ESP_H5z_G2y2z_a;
    abcd[iGrid*1890+1218] = 4.0E0*I_ESP_K5x2y_Gy3z_aa-2.0E0*1*I_ESP_H5x_Gy3z_a;
    abcd[iGrid*1890+1219] = 4.0E0*I_ESP_K4x3y_Gy3z_aa-2.0E0*1*I_ESP_H4xy_Gy3z_a-2.0E0*2*I_ESP_H4xy_Gy3z_a;
    abcd[iGrid*1890+1220] = 4.0E0*I_ESP_K4x2yz_Gy3z_aa-2.0E0*1*I_ESP_H4xz_Gy3z_a;
    abcd[iGrid*1890+1221] = 4.0E0*I_ESP_K3x4y_Gy3z_aa-2.0E0*2*I_ESP_H3x2y_Gy3z_a-2.0E0*3*I_ESP_H3x2y_Gy3z_a+2*1*I_ESP_F3x_Gy3z;
    abcd[iGrid*1890+1222] = 4.0E0*I_ESP_K3x3yz_Gy3z_aa-2.0E0*1*I_ESP_H3xyz_Gy3z_a-2.0E0*2*I_ESP_H3xyz_Gy3z_a;
    abcd[iGrid*1890+1223] = 4.0E0*I_ESP_K3x2y2z_Gy3z_aa-2.0E0*1*I_ESP_H3x2z_Gy3z_a;
    abcd[iGrid*1890+1224] = 4.0E0*I_ESP_K2x5y_Gy3z_aa-2.0E0*3*I_ESP_H2x3y_Gy3z_a-2.0E0*4*I_ESP_H2x3y_Gy3z_a+3*2*I_ESP_F2xy_Gy3z;
    abcd[iGrid*1890+1225] = 4.0E0*I_ESP_K2x4yz_Gy3z_aa-2.0E0*2*I_ESP_H2x2yz_Gy3z_a-2.0E0*3*I_ESP_H2x2yz_Gy3z_a+2*1*I_ESP_F2xz_Gy3z;
    abcd[iGrid*1890+1226] = 4.0E0*I_ESP_K2x3y2z_Gy3z_aa-2.0E0*1*I_ESP_H2xy2z_Gy3z_a-2.0E0*2*I_ESP_H2xy2z_Gy3z_a;
    abcd[iGrid*1890+1227] = 4.0E0*I_ESP_K2x2y3z_Gy3z_aa-2.0E0*1*I_ESP_H2x3z_Gy3z_a;
    abcd[iGrid*1890+1228] = 4.0E0*I_ESP_Kx6y_Gy3z_aa-2.0E0*4*I_ESP_Hx4y_Gy3z_a-2.0E0*5*I_ESP_Hx4y_Gy3z_a+4*3*I_ESP_Fx2y_Gy3z;
    abcd[iGrid*1890+1229] = 4.0E0*I_ESP_Kx5yz_Gy3z_aa-2.0E0*3*I_ESP_Hx3yz_Gy3z_a-2.0E0*4*I_ESP_Hx3yz_Gy3z_a+3*2*I_ESP_Fxyz_Gy3z;
    abcd[iGrid*1890+1230] = 4.0E0*I_ESP_Kx4y2z_Gy3z_aa-2.0E0*2*I_ESP_Hx2y2z_Gy3z_a-2.0E0*3*I_ESP_Hx2y2z_Gy3z_a+2*1*I_ESP_Fx2z_Gy3z;
    abcd[iGrid*1890+1231] = 4.0E0*I_ESP_Kx3y3z_Gy3z_aa-2.0E0*1*I_ESP_Hxy3z_Gy3z_a-2.0E0*2*I_ESP_Hxy3z_Gy3z_a;
    abcd[iGrid*1890+1232] = 4.0E0*I_ESP_Kx2y4z_Gy3z_aa-2.0E0*1*I_ESP_Hx4z_Gy3z_a;
    abcd[iGrid*1890+1233] = 4.0E0*I_ESP_K7y_Gy3z_aa-2.0E0*5*I_ESP_H5y_Gy3z_a-2.0E0*6*I_ESP_H5y_Gy3z_a+5*4*I_ESP_F3y_Gy3z;
    abcd[iGrid*1890+1234] = 4.0E0*I_ESP_K6yz_Gy3z_aa-2.0E0*4*I_ESP_H4yz_Gy3z_a-2.0E0*5*I_ESP_H4yz_Gy3z_a+4*3*I_ESP_F2yz_Gy3z;
    abcd[iGrid*1890+1235] = 4.0E0*I_ESP_K5y2z_Gy3z_aa-2.0E0*3*I_ESP_H3y2z_Gy3z_a-2.0E0*4*I_ESP_H3y2z_Gy3z_a+3*2*I_ESP_Fy2z_Gy3z;
    abcd[iGrid*1890+1236] = 4.0E0*I_ESP_K4y3z_Gy3z_aa-2.0E0*2*I_ESP_H2y3z_Gy3z_a-2.0E0*3*I_ESP_H2y3z_Gy3z_a+2*1*I_ESP_F3z_Gy3z;
    abcd[iGrid*1890+1237] = 4.0E0*I_ESP_K3y4z_Gy3z_aa-2.0E0*1*I_ESP_Hy4z_Gy3z_a-2.0E0*2*I_ESP_Hy4z_Gy3z_a;
    abcd[iGrid*1890+1238] = 4.0E0*I_ESP_K2y5z_Gy3z_aa-2.0E0*1*I_ESP_H5z_Gy3z_a;
    abcd[iGrid*1890+1239] = 4.0E0*I_ESP_K5x2y_G4z_aa-2.0E0*1*I_ESP_H5x_G4z_a;
    abcd[iGrid*1890+1240] = 4.0E0*I_ESP_K4x3y_G4z_aa-2.0E0*1*I_ESP_H4xy_G4z_a-2.0E0*2*I_ESP_H4xy_G4z_a;
    abcd[iGrid*1890+1241] = 4.0E0*I_ESP_K4x2yz_G4z_aa-2.0E0*1*I_ESP_H4xz_G4z_a;
    abcd[iGrid*1890+1242] = 4.0E0*I_ESP_K3x4y_G4z_aa-2.0E0*2*I_ESP_H3x2y_G4z_a-2.0E0*3*I_ESP_H3x2y_G4z_a+2*1*I_ESP_F3x_G4z;
    abcd[iGrid*1890+1243] = 4.0E0*I_ESP_K3x3yz_G4z_aa-2.0E0*1*I_ESP_H3xyz_G4z_a-2.0E0*2*I_ESP_H3xyz_G4z_a;
    abcd[iGrid*1890+1244] = 4.0E0*I_ESP_K3x2y2z_G4z_aa-2.0E0*1*I_ESP_H3x2z_G4z_a;
    abcd[iGrid*1890+1245] = 4.0E0*I_ESP_K2x5y_G4z_aa-2.0E0*3*I_ESP_H2x3y_G4z_a-2.0E0*4*I_ESP_H2x3y_G4z_a+3*2*I_ESP_F2xy_G4z;
    abcd[iGrid*1890+1246] = 4.0E0*I_ESP_K2x4yz_G4z_aa-2.0E0*2*I_ESP_H2x2yz_G4z_a-2.0E0*3*I_ESP_H2x2yz_G4z_a+2*1*I_ESP_F2xz_G4z;
    abcd[iGrid*1890+1247] = 4.0E0*I_ESP_K2x3y2z_G4z_aa-2.0E0*1*I_ESP_H2xy2z_G4z_a-2.0E0*2*I_ESP_H2xy2z_G4z_a;
    abcd[iGrid*1890+1248] = 4.0E0*I_ESP_K2x2y3z_G4z_aa-2.0E0*1*I_ESP_H2x3z_G4z_a;
    abcd[iGrid*1890+1249] = 4.0E0*I_ESP_Kx6y_G4z_aa-2.0E0*4*I_ESP_Hx4y_G4z_a-2.0E0*5*I_ESP_Hx4y_G4z_a+4*3*I_ESP_Fx2y_G4z;
    abcd[iGrid*1890+1250] = 4.0E0*I_ESP_Kx5yz_G4z_aa-2.0E0*3*I_ESP_Hx3yz_G4z_a-2.0E0*4*I_ESP_Hx3yz_G4z_a+3*2*I_ESP_Fxyz_G4z;
    abcd[iGrid*1890+1251] = 4.0E0*I_ESP_Kx4y2z_G4z_aa-2.0E0*2*I_ESP_Hx2y2z_G4z_a-2.0E0*3*I_ESP_Hx2y2z_G4z_a+2*1*I_ESP_Fx2z_G4z;
    abcd[iGrid*1890+1252] = 4.0E0*I_ESP_Kx3y3z_G4z_aa-2.0E0*1*I_ESP_Hxy3z_G4z_a-2.0E0*2*I_ESP_Hxy3z_G4z_a;
    abcd[iGrid*1890+1253] = 4.0E0*I_ESP_Kx2y4z_G4z_aa-2.0E0*1*I_ESP_Hx4z_G4z_a;
    abcd[iGrid*1890+1254] = 4.0E0*I_ESP_K7y_G4z_aa-2.0E0*5*I_ESP_H5y_G4z_a-2.0E0*6*I_ESP_H5y_G4z_a+5*4*I_ESP_F3y_G4z;
    abcd[iGrid*1890+1255] = 4.0E0*I_ESP_K6yz_G4z_aa-2.0E0*4*I_ESP_H4yz_G4z_a-2.0E0*5*I_ESP_H4yz_G4z_a+4*3*I_ESP_F2yz_G4z;
    abcd[iGrid*1890+1256] = 4.0E0*I_ESP_K5y2z_G4z_aa-2.0E0*3*I_ESP_H3y2z_G4z_a-2.0E0*4*I_ESP_H3y2z_G4z_a+3*2*I_ESP_Fy2z_G4z;
    abcd[iGrid*1890+1257] = 4.0E0*I_ESP_K4y3z_G4z_aa-2.0E0*2*I_ESP_H2y3z_G4z_a-2.0E0*3*I_ESP_H2y3z_G4z_a+2*1*I_ESP_F3z_G4z;
    abcd[iGrid*1890+1258] = 4.0E0*I_ESP_K3y4z_G4z_aa-2.0E0*1*I_ESP_Hy4z_G4z_a-2.0E0*2*I_ESP_Hy4z_G4z_a;
    abcd[iGrid*1890+1259] = 4.0E0*I_ESP_K2y5z_G4z_aa-2.0E0*1*I_ESP_H5z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_G_aa
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_F_G
     ************************************************************/
    abcd[iGrid*1890+1260] = 4.0E0*I_ESP_K5xyz_G4x_aa;
    abcd[iGrid*1890+1261] = 4.0E0*I_ESP_K4x2yz_G4x_aa-2.0E0*1*I_ESP_H4xz_G4x_a;
    abcd[iGrid*1890+1262] = 4.0E0*I_ESP_K4xy2z_G4x_aa-2.0E0*1*I_ESP_H4xy_G4x_a;
    abcd[iGrid*1890+1263] = 4.0E0*I_ESP_K3x3yz_G4x_aa-2.0E0*2*I_ESP_H3xyz_G4x_a;
    abcd[iGrid*1890+1264] = 4.0E0*I_ESP_K3x2y2z_G4x_aa-2.0E0*1*I_ESP_H3x2y_G4x_a-2.0E0*1*I_ESP_H3x2z_G4x_a+1*I_ESP_F3x_G4x;
    abcd[iGrid*1890+1265] = 4.0E0*I_ESP_K3xy3z_G4x_aa-2.0E0*2*I_ESP_H3xyz_G4x_a;
    abcd[iGrid*1890+1266] = 4.0E0*I_ESP_K2x4yz_G4x_aa-2.0E0*3*I_ESP_H2x2yz_G4x_a;
    abcd[iGrid*1890+1267] = 4.0E0*I_ESP_K2x3y2z_G4x_aa-2.0E0*1*I_ESP_H2x3y_G4x_a-2.0E0*2*I_ESP_H2xy2z_G4x_a+2*1*I_ESP_F2xy_G4x;
    abcd[iGrid*1890+1268] = 4.0E0*I_ESP_K2x2y3z_G4x_aa-2.0E0*2*I_ESP_H2x2yz_G4x_a-2.0E0*1*I_ESP_H2x3z_G4x_a+2*I_ESP_F2xz_G4x;
    abcd[iGrid*1890+1269] = 4.0E0*I_ESP_K2xy4z_G4x_aa-2.0E0*3*I_ESP_H2xy2z_G4x_a;
    abcd[iGrid*1890+1270] = 4.0E0*I_ESP_Kx5yz_G4x_aa-2.0E0*4*I_ESP_Hx3yz_G4x_a;
    abcd[iGrid*1890+1271] = 4.0E0*I_ESP_Kx4y2z_G4x_aa-2.0E0*1*I_ESP_Hx4y_G4x_a-2.0E0*3*I_ESP_Hx2y2z_G4x_a+3*1*I_ESP_Fx2y_G4x;
    abcd[iGrid*1890+1272] = 4.0E0*I_ESP_Kx3y3z_G4x_aa-2.0E0*2*I_ESP_Hx3yz_G4x_a-2.0E0*2*I_ESP_Hxy3z_G4x_a+2*2*I_ESP_Fxyz_G4x;
    abcd[iGrid*1890+1273] = 4.0E0*I_ESP_Kx2y4z_G4x_aa-2.0E0*3*I_ESP_Hx2y2z_G4x_a-2.0E0*1*I_ESP_Hx4z_G4x_a+3*I_ESP_Fx2z_G4x;
    abcd[iGrid*1890+1274] = 4.0E0*I_ESP_Kxy5z_G4x_aa-2.0E0*4*I_ESP_Hxy3z_G4x_a;
    abcd[iGrid*1890+1275] = 4.0E0*I_ESP_K6yz_G4x_aa-2.0E0*5*I_ESP_H4yz_G4x_a;
    abcd[iGrid*1890+1276] = 4.0E0*I_ESP_K5y2z_G4x_aa-2.0E0*1*I_ESP_H5y_G4x_a-2.0E0*4*I_ESP_H3y2z_G4x_a+4*1*I_ESP_F3y_G4x;
    abcd[iGrid*1890+1277] = 4.0E0*I_ESP_K4y3z_G4x_aa-2.0E0*2*I_ESP_H4yz_G4x_a-2.0E0*3*I_ESP_H2y3z_G4x_a+3*2*I_ESP_F2yz_G4x;
    abcd[iGrid*1890+1278] = 4.0E0*I_ESP_K3y4z_G4x_aa-2.0E0*3*I_ESP_H3y2z_G4x_a-2.0E0*2*I_ESP_Hy4z_G4x_a+2*3*I_ESP_Fy2z_G4x;
    abcd[iGrid*1890+1279] = 4.0E0*I_ESP_K2y5z_G4x_aa-2.0E0*4*I_ESP_H2y3z_G4x_a-2.0E0*1*I_ESP_H5z_G4x_a+4*I_ESP_F3z_G4x;
    abcd[iGrid*1890+1280] = 4.0E0*I_ESP_Ky6z_G4x_aa-2.0E0*5*I_ESP_Hy4z_G4x_a;
    abcd[iGrid*1890+1281] = 4.0E0*I_ESP_K5xyz_G3xy_aa;
    abcd[iGrid*1890+1282] = 4.0E0*I_ESP_K4x2yz_G3xy_aa-2.0E0*1*I_ESP_H4xz_G3xy_a;
    abcd[iGrid*1890+1283] = 4.0E0*I_ESP_K4xy2z_G3xy_aa-2.0E0*1*I_ESP_H4xy_G3xy_a;
    abcd[iGrid*1890+1284] = 4.0E0*I_ESP_K3x3yz_G3xy_aa-2.0E0*2*I_ESP_H3xyz_G3xy_a;
    abcd[iGrid*1890+1285] = 4.0E0*I_ESP_K3x2y2z_G3xy_aa-2.0E0*1*I_ESP_H3x2y_G3xy_a-2.0E0*1*I_ESP_H3x2z_G3xy_a+1*I_ESP_F3x_G3xy;
    abcd[iGrid*1890+1286] = 4.0E0*I_ESP_K3xy3z_G3xy_aa-2.0E0*2*I_ESP_H3xyz_G3xy_a;
    abcd[iGrid*1890+1287] = 4.0E0*I_ESP_K2x4yz_G3xy_aa-2.0E0*3*I_ESP_H2x2yz_G3xy_a;
    abcd[iGrid*1890+1288] = 4.0E0*I_ESP_K2x3y2z_G3xy_aa-2.0E0*1*I_ESP_H2x3y_G3xy_a-2.0E0*2*I_ESP_H2xy2z_G3xy_a+2*1*I_ESP_F2xy_G3xy;
    abcd[iGrid*1890+1289] = 4.0E0*I_ESP_K2x2y3z_G3xy_aa-2.0E0*2*I_ESP_H2x2yz_G3xy_a-2.0E0*1*I_ESP_H2x3z_G3xy_a+2*I_ESP_F2xz_G3xy;
    abcd[iGrid*1890+1290] = 4.0E0*I_ESP_K2xy4z_G3xy_aa-2.0E0*3*I_ESP_H2xy2z_G3xy_a;
    abcd[iGrid*1890+1291] = 4.0E0*I_ESP_Kx5yz_G3xy_aa-2.0E0*4*I_ESP_Hx3yz_G3xy_a;
    abcd[iGrid*1890+1292] = 4.0E0*I_ESP_Kx4y2z_G3xy_aa-2.0E0*1*I_ESP_Hx4y_G3xy_a-2.0E0*3*I_ESP_Hx2y2z_G3xy_a+3*1*I_ESP_Fx2y_G3xy;
    abcd[iGrid*1890+1293] = 4.0E0*I_ESP_Kx3y3z_G3xy_aa-2.0E0*2*I_ESP_Hx3yz_G3xy_a-2.0E0*2*I_ESP_Hxy3z_G3xy_a+2*2*I_ESP_Fxyz_G3xy;
    abcd[iGrid*1890+1294] = 4.0E0*I_ESP_Kx2y4z_G3xy_aa-2.0E0*3*I_ESP_Hx2y2z_G3xy_a-2.0E0*1*I_ESP_Hx4z_G3xy_a+3*I_ESP_Fx2z_G3xy;
    abcd[iGrid*1890+1295] = 4.0E0*I_ESP_Kxy5z_G3xy_aa-2.0E0*4*I_ESP_Hxy3z_G3xy_a;
    abcd[iGrid*1890+1296] = 4.0E0*I_ESP_K6yz_G3xy_aa-2.0E0*5*I_ESP_H4yz_G3xy_a;
    abcd[iGrid*1890+1297] = 4.0E0*I_ESP_K5y2z_G3xy_aa-2.0E0*1*I_ESP_H5y_G3xy_a-2.0E0*4*I_ESP_H3y2z_G3xy_a+4*1*I_ESP_F3y_G3xy;
    abcd[iGrid*1890+1298] = 4.0E0*I_ESP_K4y3z_G3xy_aa-2.0E0*2*I_ESP_H4yz_G3xy_a-2.0E0*3*I_ESP_H2y3z_G3xy_a+3*2*I_ESP_F2yz_G3xy;
    abcd[iGrid*1890+1299] = 4.0E0*I_ESP_K3y4z_G3xy_aa-2.0E0*3*I_ESP_H3y2z_G3xy_a-2.0E0*2*I_ESP_Hy4z_G3xy_a+2*3*I_ESP_Fy2z_G3xy;
    abcd[iGrid*1890+1300] = 4.0E0*I_ESP_K2y5z_G3xy_aa-2.0E0*4*I_ESP_H2y3z_G3xy_a-2.0E0*1*I_ESP_H5z_G3xy_a+4*I_ESP_F3z_G3xy;
    abcd[iGrid*1890+1301] = 4.0E0*I_ESP_Ky6z_G3xy_aa-2.0E0*5*I_ESP_Hy4z_G3xy_a;
    abcd[iGrid*1890+1302] = 4.0E0*I_ESP_K5xyz_G3xz_aa;
    abcd[iGrid*1890+1303] = 4.0E0*I_ESP_K4x2yz_G3xz_aa-2.0E0*1*I_ESP_H4xz_G3xz_a;
    abcd[iGrid*1890+1304] = 4.0E0*I_ESP_K4xy2z_G3xz_aa-2.0E0*1*I_ESP_H4xy_G3xz_a;
    abcd[iGrid*1890+1305] = 4.0E0*I_ESP_K3x3yz_G3xz_aa-2.0E0*2*I_ESP_H3xyz_G3xz_a;
    abcd[iGrid*1890+1306] = 4.0E0*I_ESP_K3x2y2z_G3xz_aa-2.0E0*1*I_ESP_H3x2y_G3xz_a-2.0E0*1*I_ESP_H3x2z_G3xz_a+1*I_ESP_F3x_G3xz;
    abcd[iGrid*1890+1307] = 4.0E0*I_ESP_K3xy3z_G3xz_aa-2.0E0*2*I_ESP_H3xyz_G3xz_a;
    abcd[iGrid*1890+1308] = 4.0E0*I_ESP_K2x4yz_G3xz_aa-2.0E0*3*I_ESP_H2x2yz_G3xz_a;
    abcd[iGrid*1890+1309] = 4.0E0*I_ESP_K2x3y2z_G3xz_aa-2.0E0*1*I_ESP_H2x3y_G3xz_a-2.0E0*2*I_ESP_H2xy2z_G3xz_a+2*1*I_ESP_F2xy_G3xz;
    abcd[iGrid*1890+1310] = 4.0E0*I_ESP_K2x2y3z_G3xz_aa-2.0E0*2*I_ESP_H2x2yz_G3xz_a-2.0E0*1*I_ESP_H2x3z_G3xz_a+2*I_ESP_F2xz_G3xz;
    abcd[iGrid*1890+1311] = 4.0E0*I_ESP_K2xy4z_G3xz_aa-2.0E0*3*I_ESP_H2xy2z_G3xz_a;
    abcd[iGrid*1890+1312] = 4.0E0*I_ESP_Kx5yz_G3xz_aa-2.0E0*4*I_ESP_Hx3yz_G3xz_a;
    abcd[iGrid*1890+1313] = 4.0E0*I_ESP_Kx4y2z_G3xz_aa-2.0E0*1*I_ESP_Hx4y_G3xz_a-2.0E0*3*I_ESP_Hx2y2z_G3xz_a+3*1*I_ESP_Fx2y_G3xz;
    abcd[iGrid*1890+1314] = 4.0E0*I_ESP_Kx3y3z_G3xz_aa-2.0E0*2*I_ESP_Hx3yz_G3xz_a-2.0E0*2*I_ESP_Hxy3z_G3xz_a+2*2*I_ESP_Fxyz_G3xz;
    abcd[iGrid*1890+1315] = 4.0E0*I_ESP_Kx2y4z_G3xz_aa-2.0E0*3*I_ESP_Hx2y2z_G3xz_a-2.0E0*1*I_ESP_Hx4z_G3xz_a+3*I_ESP_Fx2z_G3xz;
    abcd[iGrid*1890+1316] = 4.0E0*I_ESP_Kxy5z_G3xz_aa-2.0E0*4*I_ESP_Hxy3z_G3xz_a;
    abcd[iGrid*1890+1317] = 4.0E0*I_ESP_K6yz_G3xz_aa-2.0E0*5*I_ESP_H4yz_G3xz_a;
    abcd[iGrid*1890+1318] = 4.0E0*I_ESP_K5y2z_G3xz_aa-2.0E0*1*I_ESP_H5y_G3xz_a-2.0E0*4*I_ESP_H3y2z_G3xz_a+4*1*I_ESP_F3y_G3xz;
    abcd[iGrid*1890+1319] = 4.0E0*I_ESP_K4y3z_G3xz_aa-2.0E0*2*I_ESP_H4yz_G3xz_a-2.0E0*3*I_ESP_H2y3z_G3xz_a+3*2*I_ESP_F2yz_G3xz;
    abcd[iGrid*1890+1320] = 4.0E0*I_ESP_K3y4z_G3xz_aa-2.0E0*3*I_ESP_H3y2z_G3xz_a-2.0E0*2*I_ESP_Hy4z_G3xz_a+2*3*I_ESP_Fy2z_G3xz;
    abcd[iGrid*1890+1321] = 4.0E0*I_ESP_K2y5z_G3xz_aa-2.0E0*4*I_ESP_H2y3z_G3xz_a-2.0E0*1*I_ESP_H5z_G3xz_a+4*I_ESP_F3z_G3xz;
    abcd[iGrid*1890+1322] = 4.0E0*I_ESP_Ky6z_G3xz_aa-2.0E0*5*I_ESP_Hy4z_G3xz_a;
    abcd[iGrid*1890+1323] = 4.0E0*I_ESP_K5xyz_G2x2y_aa;
    abcd[iGrid*1890+1324] = 4.0E0*I_ESP_K4x2yz_G2x2y_aa-2.0E0*1*I_ESP_H4xz_G2x2y_a;
    abcd[iGrid*1890+1325] = 4.0E0*I_ESP_K4xy2z_G2x2y_aa-2.0E0*1*I_ESP_H4xy_G2x2y_a;
    abcd[iGrid*1890+1326] = 4.0E0*I_ESP_K3x3yz_G2x2y_aa-2.0E0*2*I_ESP_H3xyz_G2x2y_a;
    abcd[iGrid*1890+1327] = 4.0E0*I_ESP_K3x2y2z_G2x2y_aa-2.0E0*1*I_ESP_H3x2y_G2x2y_a-2.0E0*1*I_ESP_H3x2z_G2x2y_a+1*I_ESP_F3x_G2x2y;
    abcd[iGrid*1890+1328] = 4.0E0*I_ESP_K3xy3z_G2x2y_aa-2.0E0*2*I_ESP_H3xyz_G2x2y_a;
    abcd[iGrid*1890+1329] = 4.0E0*I_ESP_K2x4yz_G2x2y_aa-2.0E0*3*I_ESP_H2x2yz_G2x2y_a;
    abcd[iGrid*1890+1330] = 4.0E0*I_ESP_K2x3y2z_G2x2y_aa-2.0E0*1*I_ESP_H2x3y_G2x2y_a-2.0E0*2*I_ESP_H2xy2z_G2x2y_a+2*1*I_ESP_F2xy_G2x2y;
    abcd[iGrid*1890+1331] = 4.0E0*I_ESP_K2x2y3z_G2x2y_aa-2.0E0*2*I_ESP_H2x2yz_G2x2y_a-2.0E0*1*I_ESP_H2x3z_G2x2y_a+2*I_ESP_F2xz_G2x2y;
    abcd[iGrid*1890+1332] = 4.0E0*I_ESP_K2xy4z_G2x2y_aa-2.0E0*3*I_ESP_H2xy2z_G2x2y_a;
    abcd[iGrid*1890+1333] = 4.0E0*I_ESP_Kx5yz_G2x2y_aa-2.0E0*4*I_ESP_Hx3yz_G2x2y_a;
    abcd[iGrid*1890+1334] = 4.0E0*I_ESP_Kx4y2z_G2x2y_aa-2.0E0*1*I_ESP_Hx4y_G2x2y_a-2.0E0*3*I_ESP_Hx2y2z_G2x2y_a+3*1*I_ESP_Fx2y_G2x2y;
    abcd[iGrid*1890+1335] = 4.0E0*I_ESP_Kx3y3z_G2x2y_aa-2.0E0*2*I_ESP_Hx3yz_G2x2y_a-2.0E0*2*I_ESP_Hxy3z_G2x2y_a+2*2*I_ESP_Fxyz_G2x2y;
    abcd[iGrid*1890+1336] = 4.0E0*I_ESP_Kx2y4z_G2x2y_aa-2.0E0*3*I_ESP_Hx2y2z_G2x2y_a-2.0E0*1*I_ESP_Hx4z_G2x2y_a+3*I_ESP_Fx2z_G2x2y;
    abcd[iGrid*1890+1337] = 4.0E0*I_ESP_Kxy5z_G2x2y_aa-2.0E0*4*I_ESP_Hxy3z_G2x2y_a;
    abcd[iGrid*1890+1338] = 4.0E0*I_ESP_K6yz_G2x2y_aa-2.0E0*5*I_ESP_H4yz_G2x2y_a;
    abcd[iGrid*1890+1339] = 4.0E0*I_ESP_K5y2z_G2x2y_aa-2.0E0*1*I_ESP_H5y_G2x2y_a-2.0E0*4*I_ESP_H3y2z_G2x2y_a+4*1*I_ESP_F3y_G2x2y;
    abcd[iGrid*1890+1340] = 4.0E0*I_ESP_K4y3z_G2x2y_aa-2.0E0*2*I_ESP_H4yz_G2x2y_a-2.0E0*3*I_ESP_H2y3z_G2x2y_a+3*2*I_ESP_F2yz_G2x2y;
    abcd[iGrid*1890+1341] = 4.0E0*I_ESP_K3y4z_G2x2y_aa-2.0E0*3*I_ESP_H3y2z_G2x2y_a-2.0E0*2*I_ESP_Hy4z_G2x2y_a+2*3*I_ESP_Fy2z_G2x2y;
    abcd[iGrid*1890+1342] = 4.0E0*I_ESP_K2y5z_G2x2y_aa-2.0E0*4*I_ESP_H2y3z_G2x2y_a-2.0E0*1*I_ESP_H5z_G2x2y_a+4*I_ESP_F3z_G2x2y;
    abcd[iGrid*1890+1343] = 4.0E0*I_ESP_Ky6z_G2x2y_aa-2.0E0*5*I_ESP_Hy4z_G2x2y_a;
    abcd[iGrid*1890+1344] = 4.0E0*I_ESP_K5xyz_G2xyz_aa;
    abcd[iGrid*1890+1345] = 4.0E0*I_ESP_K4x2yz_G2xyz_aa-2.0E0*1*I_ESP_H4xz_G2xyz_a;
    abcd[iGrid*1890+1346] = 4.0E0*I_ESP_K4xy2z_G2xyz_aa-2.0E0*1*I_ESP_H4xy_G2xyz_a;
    abcd[iGrid*1890+1347] = 4.0E0*I_ESP_K3x3yz_G2xyz_aa-2.0E0*2*I_ESP_H3xyz_G2xyz_a;
    abcd[iGrid*1890+1348] = 4.0E0*I_ESP_K3x2y2z_G2xyz_aa-2.0E0*1*I_ESP_H3x2y_G2xyz_a-2.0E0*1*I_ESP_H3x2z_G2xyz_a+1*I_ESP_F3x_G2xyz;
    abcd[iGrid*1890+1349] = 4.0E0*I_ESP_K3xy3z_G2xyz_aa-2.0E0*2*I_ESP_H3xyz_G2xyz_a;
    abcd[iGrid*1890+1350] = 4.0E0*I_ESP_K2x4yz_G2xyz_aa-2.0E0*3*I_ESP_H2x2yz_G2xyz_a;
    abcd[iGrid*1890+1351] = 4.0E0*I_ESP_K2x3y2z_G2xyz_aa-2.0E0*1*I_ESP_H2x3y_G2xyz_a-2.0E0*2*I_ESP_H2xy2z_G2xyz_a+2*1*I_ESP_F2xy_G2xyz;
    abcd[iGrid*1890+1352] = 4.0E0*I_ESP_K2x2y3z_G2xyz_aa-2.0E0*2*I_ESP_H2x2yz_G2xyz_a-2.0E0*1*I_ESP_H2x3z_G2xyz_a+2*I_ESP_F2xz_G2xyz;
    abcd[iGrid*1890+1353] = 4.0E0*I_ESP_K2xy4z_G2xyz_aa-2.0E0*3*I_ESP_H2xy2z_G2xyz_a;
    abcd[iGrid*1890+1354] = 4.0E0*I_ESP_Kx5yz_G2xyz_aa-2.0E0*4*I_ESP_Hx3yz_G2xyz_a;
    abcd[iGrid*1890+1355] = 4.0E0*I_ESP_Kx4y2z_G2xyz_aa-2.0E0*1*I_ESP_Hx4y_G2xyz_a-2.0E0*3*I_ESP_Hx2y2z_G2xyz_a+3*1*I_ESP_Fx2y_G2xyz;
    abcd[iGrid*1890+1356] = 4.0E0*I_ESP_Kx3y3z_G2xyz_aa-2.0E0*2*I_ESP_Hx3yz_G2xyz_a-2.0E0*2*I_ESP_Hxy3z_G2xyz_a+2*2*I_ESP_Fxyz_G2xyz;
    abcd[iGrid*1890+1357] = 4.0E0*I_ESP_Kx2y4z_G2xyz_aa-2.0E0*3*I_ESP_Hx2y2z_G2xyz_a-2.0E0*1*I_ESP_Hx4z_G2xyz_a+3*I_ESP_Fx2z_G2xyz;
    abcd[iGrid*1890+1358] = 4.0E0*I_ESP_Kxy5z_G2xyz_aa-2.0E0*4*I_ESP_Hxy3z_G2xyz_a;
    abcd[iGrid*1890+1359] = 4.0E0*I_ESP_K6yz_G2xyz_aa-2.0E0*5*I_ESP_H4yz_G2xyz_a;
    abcd[iGrid*1890+1360] = 4.0E0*I_ESP_K5y2z_G2xyz_aa-2.0E0*1*I_ESP_H5y_G2xyz_a-2.0E0*4*I_ESP_H3y2z_G2xyz_a+4*1*I_ESP_F3y_G2xyz;
    abcd[iGrid*1890+1361] = 4.0E0*I_ESP_K4y3z_G2xyz_aa-2.0E0*2*I_ESP_H4yz_G2xyz_a-2.0E0*3*I_ESP_H2y3z_G2xyz_a+3*2*I_ESP_F2yz_G2xyz;
    abcd[iGrid*1890+1362] = 4.0E0*I_ESP_K3y4z_G2xyz_aa-2.0E0*3*I_ESP_H3y2z_G2xyz_a-2.0E0*2*I_ESP_Hy4z_G2xyz_a+2*3*I_ESP_Fy2z_G2xyz;
    abcd[iGrid*1890+1363] = 4.0E0*I_ESP_K2y5z_G2xyz_aa-2.0E0*4*I_ESP_H2y3z_G2xyz_a-2.0E0*1*I_ESP_H5z_G2xyz_a+4*I_ESP_F3z_G2xyz;
    abcd[iGrid*1890+1364] = 4.0E0*I_ESP_Ky6z_G2xyz_aa-2.0E0*5*I_ESP_Hy4z_G2xyz_a;
    abcd[iGrid*1890+1365] = 4.0E0*I_ESP_K5xyz_G2x2z_aa;
    abcd[iGrid*1890+1366] = 4.0E0*I_ESP_K4x2yz_G2x2z_aa-2.0E0*1*I_ESP_H4xz_G2x2z_a;
    abcd[iGrid*1890+1367] = 4.0E0*I_ESP_K4xy2z_G2x2z_aa-2.0E0*1*I_ESP_H4xy_G2x2z_a;
    abcd[iGrid*1890+1368] = 4.0E0*I_ESP_K3x3yz_G2x2z_aa-2.0E0*2*I_ESP_H3xyz_G2x2z_a;
    abcd[iGrid*1890+1369] = 4.0E0*I_ESP_K3x2y2z_G2x2z_aa-2.0E0*1*I_ESP_H3x2y_G2x2z_a-2.0E0*1*I_ESP_H3x2z_G2x2z_a+1*I_ESP_F3x_G2x2z;
    abcd[iGrid*1890+1370] = 4.0E0*I_ESP_K3xy3z_G2x2z_aa-2.0E0*2*I_ESP_H3xyz_G2x2z_a;
    abcd[iGrid*1890+1371] = 4.0E0*I_ESP_K2x4yz_G2x2z_aa-2.0E0*3*I_ESP_H2x2yz_G2x2z_a;
    abcd[iGrid*1890+1372] = 4.0E0*I_ESP_K2x3y2z_G2x2z_aa-2.0E0*1*I_ESP_H2x3y_G2x2z_a-2.0E0*2*I_ESP_H2xy2z_G2x2z_a+2*1*I_ESP_F2xy_G2x2z;
    abcd[iGrid*1890+1373] = 4.0E0*I_ESP_K2x2y3z_G2x2z_aa-2.0E0*2*I_ESP_H2x2yz_G2x2z_a-2.0E0*1*I_ESP_H2x3z_G2x2z_a+2*I_ESP_F2xz_G2x2z;
    abcd[iGrid*1890+1374] = 4.0E0*I_ESP_K2xy4z_G2x2z_aa-2.0E0*3*I_ESP_H2xy2z_G2x2z_a;
    abcd[iGrid*1890+1375] = 4.0E0*I_ESP_Kx5yz_G2x2z_aa-2.0E0*4*I_ESP_Hx3yz_G2x2z_a;
    abcd[iGrid*1890+1376] = 4.0E0*I_ESP_Kx4y2z_G2x2z_aa-2.0E0*1*I_ESP_Hx4y_G2x2z_a-2.0E0*3*I_ESP_Hx2y2z_G2x2z_a+3*1*I_ESP_Fx2y_G2x2z;
    abcd[iGrid*1890+1377] = 4.0E0*I_ESP_Kx3y3z_G2x2z_aa-2.0E0*2*I_ESP_Hx3yz_G2x2z_a-2.0E0*2*I_ESP_Hxy3z_G2x2z_a+2*2*I_ESP_Fxyz_G2x2z;
    abcd[iGrid*1890+1378] = 4.0E0*I_ESP_Kx2y4z_G2x2z_aa-2.0E0*3*I_ESP_Hx2y2z_G2x2z_a-2.0E0*1*I_ESP_Hx4z_G2x2z_a+3*I_ESP_Fx2z_G2x2z;
    abcd[iGrid*1890+1379] = 4.0E0*I_ESP_Kxy5z_G2x2z_aa-2.0E0*4*I_ESP_Hxy3z_G2x2z_a;
    abcd[iGrid*1890+1380] = 4.0E0*I_ESP_K6yz_G2x2z_aa-2.0E0*5*I_ESP_H4yz_G2x2z_a;
    abcd[iGrid*1890+1381] = 4.0E0*I_ESP_K5y2z_G2x2z_aa-2.0E0*1*I_ESP_H5y_G2x2z_a-2.0E0*4*I_ESP_H3y2z_G2x2z_a+4*1*I_ESP_F3y_G2x2z;
    abcd[iGrid*1890+1382] = 4.0E0*I_ESP_K4y3z_G2x2z_aa-2.0E0*2*I_ESP_H4yz_G2x2z_a-2.0E0*3*I_ESP_H2y3z_G2x2z_a+3*2*I_ESP_F2yz_G2x2z;
    abcd[iGrid*1890+1383] = 4.0E0*I_ESP_K3y4z_G2x2z_aa-2.0E0*3*I_ESP_H3y2z_G2x2z_a-2.0E0*2*I_ESP_Hy4z_G2x2z_a+2*3*I_ESP_Fy2z_G2x2z;
    abcd[iGrid*1890+1384] = 4.0E0*I_ESP_K2y5z_G2x2z_aa-2.0E0*4*I_ESP_H2y3z_G2x2z_a-2.0E0*1*I_ESP_H5z_G2x2z_a+4*I_ESP_F3z_G2x2z;
    abcd[iGrid*1890+1385] = 4.0E0*I_ESP_Ky6z_G2x2z_aa-2.0E0*5*I_ESP_Hy4z_G2x2z_a;
    abcd[iGrid*1890+1386] = 4.0E0*I_ESP_K5xyz_Gx3y_aa;
    abcd[iGrid*1890+1387] = 4.0E0*I_ESP_K4x2yz_Gx3y_aa-2.0E0*1*I_ESP_H4xz_Gx3y_a;
    abcd[iGrid*1890+1388] = 4.0E0*I_ESP_K4xy2z_Gx3y_aa-2.0E0*1*I_ESP_H4xy_Gx3y_a;
    abcd[iGrid*1890+1389] = 4.0E0*I_ESP_K3x3yz_Gx3y_aa-2.0E0*2*I_ESP_H3xyz_Gx3y_a;
    abcd[iGrid*1890+1390] = 4.0E0*I_ESP_K3x2y2z_Gx3y_aa-2.0E0*1*I_ESP_H3x2y_Gx3y_a-2.0E0*1*I_ESP_H3x2z_Gx3y_a+1*I_ESP_F3x_Gx3y;
    abcd[iGrid*1890+1391] = 4.0E0*I_ESP_K3xy3z_Gx3y_aa-2.0E0*2*I_ESP_H3xyz_Gx3y_a;
    abcd[iGrid*1890+1392] = 4.0E0*I_ESP_K2x4yz_Gx3y_aa-2.0E0*3*I_ESP_H2x2yz_Gx3y_a;
    abcd[iGrid*1890+1393] = 4.0E0*I_ESP_K2x3y2z_Gx3y_aa-2.0E0*1*I_ESP_H2x3y_Gx3y_a-2.0E0*2*I_ESP_H2xy2z_Gx3y_a+2*1*I_ESP_F2xy_Gx3y;
    abcd[iGrid*1890+1394] = 4.0E0*I_ESP_K2x2y3z_Gx3y_aa-2.0E0*2*I_ESP_H2x2yz_Gx3y_a-2.0E0*1*I_ESP_H2x3z_Gx3y_a+2*I_ESP_F2xz_Gx3y;
    abcd[iGrid*1890+1395] = 4.0E0*I_ESP_K2xy4z_Gx3y_aa-2.0E0*3*I_ESP_H2xy2z_Gx3y_a;
    abcd[iGrid*1890+1396] = 4.0E0*I_ESP_Kx5yz_Gx3y_aa-2.0E0*4*I_ESP_Hx3yz_Gx3y_a;
    abcd[iGrid*1890+1397] = 4.0E0*I_ESP_Kx4y2z_Gx3y_aa-2.0E0*1*I_ESP_Hx4y_Gx3y_a-2.0E0*3*I_ESP_Hx2y2z_Gx3y_a+3*1*I_ESP_Fx2y_Gx3y;
    abcd[iGrid*1890+1398] = 4.0E0*I_ESP_Kx3y3z_Gx3y_aa-2.0E0*2*I_ESP_Hx3yz_Gx3y_a-2.0E0*2*I_ESP_Hxy3z_Gx3y_a+2*2*I_ESP_Fxyz_Gx3y;
    abcd[iGrid*1890+1399] = 4.0E0*I_ESP_Kx2y4z_Gx3y_aa-2.0E0*3*I_ESP_Hx2y2z_Gx3y_a-2.0E0*1*I_ESP_Hx4z_Gx3y_a+3*I_ESP_Fx2z_Gx3y;
    abcd[iGrid*1890+1400] = 4.0E0*I_ESP_Kxy5z_Gx3y_aa-2.0E0*4*I_ESP_Hxy3z_Gx3y_a;
    abcd[iGrid*1890+1401] = 4.0E0*I_ESP_K6yz_Gx3y_aa-2.0E0*5*I_ESP_H4yz_Gx3y_a;
    abcd[iGrid*1890+1402] = 4.0E0*I_ESP_K5y2z_Gx3y_aa-2.0E0*1*I_ESP_H5y_Gx3y_a-2.0E0*4*I_ESP_H3y2z_Gx3y_a+4*1*I_ESP_F3y_Gx3y;
    abcd[iGrid*1890+1403] = 4.0E0*I_ESP_K4y3z_Gx3y_aa-2.0E0*2*I_ESP_H4yz_Gx3y_a-2.0E0*3*I_ESP_H2y3z_Gx3y_a+3*2*I_ESP_F2yz_Gx3y;
    abcd[iGrid*1890+1404] = 4.0E0*I_ESP_K3y4z_Gx3y_aa-2.0E0*3*I_ESP_H3y2z_Gx3y_a-2.0E0*2*I_ESP_Hy4z_Gx3y_a+2*3*I_ESP_Fy2z_Gx3y;
    abcd[iGrid*1890+1405] = 4.0E0*I_ESP_K2y5z_Gx3y_aa-2.0E0*4*I_ESP_H2y3z_Gx3y_a-2.0E0*1*I_ESP_H5z_Gx3y_a+4*I_ESP_F3z_Gx3y;
    abcd[iGrid*1890+1406] = 4.0E0*I_ESP_Ky6z_Gx3y_aa-2.0E0*5*I_ESP_Hy4z_Gx3y_a;
    abcd[iGrid*1890+1407] = 4.0E0*I_ESP_K5xyz_Gx2yz_aa;
    abcd[iGrid*1890+1408] = 4.0E0*I_ESP_K4x2yz_Gx2yz_aa-2.0E0*1*I_ESP_H4xz_Gx2yz_a;
    abcd[iGrid*1890+1409] = 4.0E0*I_ESP_K4xy2z_Gx2yz_aa-2.0E0*1*I_ESP_H4xy_Gx2yz_a;
    abcd[iGrid*1890+1410] = 4.0E0*I_ESP_K3x3yz_Gx2yz_aa-2.0E0*2*I_ESP_H3xyz_Gx2yz_a;
    abcd[iGrid*1890+1411] = 4.0E0*I_ESP_K3x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_H3x2y_Gx2yz_a-2.0E0*1*I_ESP_H3x2z_Gx2yz_a+1*I_ESP_F3x_Gx2yz;
    abcd[iGrid*1890+1412] = 4.0E0*I_ESP_K3xy3z_Gx2yz_aa-2.0E0*2*I_ESP_H3xyz_Gx2yz_a;
    abcd[iGrid*1890+1413] = 4.0E0*I_ESP_K2x4yz_Gx2yz_aa-2.0E0*3*I_ESP_H2x2yz_Gx2yz_a;
    abcd[iGrid*1890+1414] = 4.0E0*I_ESP_K2x3y2z_Gx2yz_aa-2.0E0*1*I_ESP_H2x3y_Gx2yz_a-2.0E0*2*I_ESP_H2xy2z_Gx2yz_a+2*1*I_ESP_F2xy_Gx2yz;
    abcd[iGrid*1890+1415] = 4.0E0*I_ESP_K2x2y3z_Gx2yz_aa-2.0E0*2*I_ESP_H2x2yz_Gx2yz_a-2.0E0*1*I_ESP_H2x3z_Gx2yz_a+2*I_ESP_F2xz_Gx2yz;
    abcd[iGrid*1890+1416] = 4.0E0*I_ESP_K2xy4z_Gx2yz_aa-2.0E0*3*I_ESP_H2xy2z_Gx2yz_a;
    abcd[iGrid*1890+1417] = 4.0E0*I_ESP_Kx5yz_Gx2yz_aa-2.0E0*4*I_ESP_Hx3yz_Gx2yz_a;
    abcd[iGrid*1890+1418] = 4.0E0*I_ESP_Kx4y2z_Gx2yz_aa-2.0E0*1*I_ESP_Hx4y_Gx2yz_a-2.0E0*3*I_ESP_Hx2y2z_Gx2yz_a+3*1*I_ESP_Fx2y_Gx2yz;
    abcd[iGrid*1890+1419] = 4.0E0*I_ESP_Kx3y3z_Gx2yz_aa-2.0E0*2*I_ESP_Hx3yz_Gx2yz_a-2.0E0*2*I_ESP_Hxy3z_Gx2yz_a+2*2*I_ESP_Fxyz_Gx2yz;
    abcd[iGrid*1890+1420] = 4.0E0*I_ESP_Kx2y4z_Gx2yz_aa-2.0E0*3*I_ESP_Hx2y2z_Gx2yz_a-2.0E0*1*I_ESP_Hx4z_Gx2yz_a+3*I_ESP_Fx2z_Gx2yz;
    abcd[iGrid*1890+1421] = 4.0E0*I_ESP_Kxy5z_Gx2yz_aa-2.0E0*4*I_ESP_Hxy3z_Gx2yz_a;
    abcd[iGrid*1890+1422] = 4.0E0*I_ESP_K6yz_Gx2yz_aa-2.0E0*5*I_ESP_H4yz_Gx2yz_a;
    abcd[iGrid*1890+1423] = 4.0E0*I_ESP_K5y2z_Gx2yz_aa-2.0E0*1*I_ESP_H5y_Gx2yz_a-2.0E0*4*I_ESP_H3y2z_Gx2yz_a+4*1*I_ESP_F3y_Gx2yz;
    abcd[iGrid*1890+1424] = 4.0E0*I_ESP_K4y3z_Gx2yz_aa-2.0E0*2*I_ESP_H4yz_Gx2yz_a-2.0E0*3*I_ESP_H2y3z_Gx2yz_a+3*2*I_ESP_F2yz_Gx2yz;
    abcd[iGrid*1890+1425] = 4.0E0*I_ESP_K3y4z_Gx2yz_aa-2.0E0*3*I_ESP_H3y2z_Gx2yz_a-2.0E0*2*I_ESP_Hy4z_Gx2yz_a+2*3*I_ESP_Fy2z_Gx2yz;
    abcd[iGrid*1890+1426] = 4.0E0*I_ESP_K2y5z_Gx2yz_aa-2.0E0*4*I_ESP_H2y3z_Gx2yz_a-2.0E0*1*I_ESP_H5z_Gx2yz_a+4*I_ESP_F3z_Gx2yz;
    abcd[iGrid*1890+1427] = 4.0E0*I_ESP_Ky6z_Gx2yz_aa-2.0E0*5*I_ESP_Hy4z_Gx2yz_a;
    abcd[iGrid*1890+1428] = 4.0E0*I_ESP_K5xyz_Gxy2z_aa;
    abcd[iGrid*1890+1429] = 4.0E0*I_ESP_K4x2yz_Gxy2z_aa-2.0E0*1*I_ESP_H4xz_Gxy2z_a;
    abcd[iGrid*1890+1430] = 4.0E0*I_ESP_K4xy2z_Gxy2z_aa-2.0E0*1*I_ESP_H4xy_Gxy2z_a;
    abcd[iGrid*1890+1431] = 4.0E0*I_ESP_K3x3yz_Gxy2z_aa-2.0E0*2*I_ESP_H3xyz_Gxy2z_a;
    abcd[iGrid*1890+1432] = 4.0E0*I_ESP_K3x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_H3x2y_Gxy2z_a-2.0E0*1*I_ESP_H3x2z_Gxy2z_a+1*I_ESP_F3x_Gxy2z;
    abcd[iGrid*1890+1433] = 4.0E0*I_ESP_K3xy3z_Gxy2z_aa-2.0E0*2*I_ESP_H3xyz_Gxy2z_a;
    abcd[iGrid*1890+1434] = 4.0E0*I_ESP_K2x4yz_Gxy2z_aa-2.0E0*3*I_ESP_H2x2yz_Gxy2z_a;
    abcd[iGrid*1890+1435] = 4.0E0*I_ESP_K2x3y2z_Gxy2z_aa-2.0E0*1*I_ESP_H2x3y_Gxy2z_a-2.0E0*2*I_ESP_H2xy2z_Gxy2z_a+2*1*I_ESP_F2xy_Gxy2z;
    abcd[iGrid*1890+1436] = 4.0E0*I_ESP_K2x2y3z_Gxy2z_aa-2.0E0*2*I_ESP_H2x2yz_Gxy2z_a-2.0E0*1*I_ESP_H2x3z_Gxy2z_a+2*I_ESP_F2xz_Gxy2z;
    abcd[iGrid*1890+1437] = 4.0E0*I_ESP_K2xy4z_Gxy2z_aa-2.0E0*3*I_ESP_H2xy2z_Gxy2z_a;
    abcd[iGrid*1890+1438] = 4.0E0*I_ESP_Kx5yz_Gxy2z_aa-2.0E0*4*I_ESP_Hx3yz_Gxy2z_a;
    abcd[iGrid*1890+1439] = 4.0E0*I_ESP_Kx4y2z_Gxy2z_aa-2.0E0*1*I_ESP_Hx4y_Gxy2z_a-2.0E0*3*I_ESP_Hx2y2z_Gxy2z_a+3*1*I_ESP_Fx2y_Gxy2z;
    abcd[iGrid*1890+1440] = 4.0E0*I_ESP_Kx3y3z_Gxy2z_aa-2.0E0*2*I_ESP_Hx3yz_Gxy2z_a-2.0E0*2*I_ESP_Hxy3z_Gxy2z_a+2*2*I_ESP_Fxyz_Gxy2z;
    abcd[iGrid*1890+1441] = 4.0E0*I_ESP_Kx2y4z_Gxy2z_aa-2.0E0*3*I_ESP_Hx2y2z_Gxy2z_a-2.0E0*1*I_ESP_Hx4z_Gxy2z_a+3*I_ESP_Fx2z_Gxy2z;
    abcd[iGrid*1890+1442] = 4.0E0*I_ESP_Kxy5z_Gxy2z_aa-2.0E0*4*I_ESP_Hxy3z_Gxy2z_a;
    abcd[iGrid*1890+1443] = 4.0E0*I_ESP_K6yz_Gxy2z_aa-2.0E0*5*I_ESP_H4yz_Gxy2z_a;
    abcd[iGrid*1890+1444] = 4.0E0*I_ESP_K5y2z_Gxy2z_aa-2.0E0*1*I_ESP_H5y_Gxy2z_a-2.0E0*4*I_ESP_H3y2z_Gxy2z_a+4*1*I_ESP_F3y_Gxy2z;
    abcd[iGrid*1890+1445] = 4.0E0*I_ESP_K4y3z_Gxy2z_aa-2.0E0*2*I_ESP_H4yz_Gxy2z_a-2.0E0*3*I_ESP_H2y3z_Gxy2z_a+3*2*I_ESP_F2yz_Gxy2z;
    abcd[iGrid*1890+1446] = 4.0E0*I_ESP_K3y4z_Gxy2z_aa-2.0E0*3*I_ESP_H3y2z_Gxy2z_a-2.0E0*2*I_ESP_Hy4z_Gxy2z_a+2*3*I_ESP_Fy2z_Gxy2z;
    abcd[iGrid*1890+1447] = 4.0E0*I_ESP_K2y5z_Gxy2z_aa-2.0E0*4*I_ESP_H2y3z_Gxy2z_a-2.0E0*1*I_ESP_H5z_Gxy2z_a+4*I_ESP_F3z_Gxy2z;
    abcd[iGrid*1890+1448] = 4.0E0*I_ESP_Ky6z_Gxy2z_aa-2.0E0*5*I_ESP_Hy4z_Gxy2z_a;
    abcd[iGrid*1890+1449] = 4.0E0*I_ESP_K5xyz_Gx3z_aa;
    abcd[iGrid*1890+1450] = 4.0E0*I_ESP_K4x2yz_Gx3z_aa-2.0E0*1*I_ESP_H4xz_Gx3z_a;
    abcd[iGrid*1890+1451] = 4.0E0*I_ESP_K4xy2z_Gx3z_aa-2.0E0*1*I_ESP_H4xy_Gx3z_a;
    abcd[iGrid*1890+1452] = 4.0E0*I_ESP_K3x3yz_Gx3z_aa-2.0E0*2*I_ESP_H3xyz_Gx3z_a;
    abcd[iGrid*1890+1453] = 4.0E0*I_ESP_K3x2y2z_Gx3z_aa-2.0E0*1*I_ESP_H3x2y_Gx3z_a-2.0E0*1*I_ESP_H3x2z_Gx3z_a+1*I_ESP_F3x_Gx3z;
    abcd[iGrid*1890+1454] = 4.0E0*I_ESP_K3xy3z_Gx3z_aa-2.0E0*2*I_ESP_H3xyz_Gx3z_a;
    abcd[iGrid*1890+1455] = 4.0E0*I_ESP_K2x4yz_Gx3z_aa-2.0E0*3*I_ESP_H2x2yz_Gx3z_a;
    abcd[iGrid*1890+1456] = 4.0E0*I_ESP_K2x3y2z_Gx3z_aa-2.0E0*1*I_ESP_H2x3y_Gx3z_a-2.0E0*2*I_ESP_H2xy2z_Gx3z_a+2*1*I_ESP_F2xy_Gx3z;
    abcd[iGrid*1890+1457] = 4.0E0*I_ESP_K2x2y3z_Gx3z_aa-2.0E0*2*I_ESP_H2x2yz_Gx3z_a-2.0E0*1*I_ESP_H2x3z_Gx3z_a+2*I_ESP_F2xz_Gx3z;
    abcd[iGrid*1890+1458] = 4.0E0*I_ESP_K2xy4z_Gx3z_aa-2.0E0*3*I_ESP_H2xy2z_Gx3z_a;
    abcd[iGrid*1890+1459] = 4.0E0*I_ESP_Kx5yz_Gx3z_aa-2.0E0*4*I_ESP_Hx3yz_Gx3z_a;
    abcd[iGrid*1890+1460] = 4.0E0*I_ESP_Kx4y2z_Gx3z_aa-2.0E0*1*I_ESP_Hx4y_Gx3z_a-2.0E0*3*I_ESP_Hx2y2z_Gx3z_a+3*1*I_ESP_Fx2y_Gx3z;
    abcd[iGrid*1890+1461] = 4.0E0*I_ESP_Kx3y3z_Gx3z_aa-2.0E0*2*I_ESP_Hx3yz_Gx3z_a-2.0E0*2*I_ESP_Hxy3z_Gx3z_a+2*2*I_ESP_Fxyz_Gx3z;
    abcd[iGrid*1890+1462] = 4.0E0*I_ESP_Kx2y4z_Gx3z_aa-2.0E0*3*I_ESP_Hx2y2z_Gx3z_a-2.0E0*1*I_ESP_Hx4z_Gx3z_a+3*I_ESP_Fx2z_Gx3z;
    abcd[iGrid*1890+1463] = 4.0E0*I_ESP_Kxy5z_Gx3z_aa-2.0E0*4*I_ESP_Hxy3z_Gx3z_a;
    abcd[iGrid*1890+1464] = 4.0E0*I_ESP_K6yz_Gx3z_aa-2.0E0*5*I_ESP_H4yz_Gx3z_a;
    abcd[iGrid*1890+1465] = 4.0E0*I_ESP_K5y2z_Gx3z_aa-2.0E0*1*I_ESP_H5y_Gx3z_a-2.0E0*4*I_ESP_H3y2z_Gx3z_a+4*1*I_ESP_F3y_Gx3z;
    abcd[iGrid*1890+1466] = 4.0E0*I_ESP_K4y3z_Gx3z_aa-2.0E0*2*I_ESP_H4yz_Gx3z_a-2.0E0*3*I_ESP_H2y3z_Gx3z_a+3*2*I_ESP_F2yz_Gx3z;
    abcd[iGrid*1890+1467] = 4.0E0*I_ESP_K3y4z_Gx3z_aa-2.0E0*3*I_ESP_H3y2z_Gx3z_a-2.0E0*2*I_ESP_Hy4z_Gx3z_a+2*3*I_ESP_Fy2z_Gx3z;
    abcd[iGrid*1890+1468] = 4.0E0*I_ESP_K2y5z_Gx3z_aa-2.0E0*4*I_ESP_H2y3z_Gx3z_a-2.0E0*1*I_ESP_H5z_Gx3z_a+4*I_ESP_F3z_Gx3z;
    abcd[iGrid*1890+1469] = 4.0E0*I_ESP_Ky6z_Gx3z_aa-2.0E0*5*I_ESP_Hy4z_Gx3z_a;
    abcd[iGrid*1890+1470] = 4.0E0*I_ESP_K5xyz_G4y_aa;
    abcd[iGrid*1890+1471] = 4.0E0*I_ESP_K4x2yz_G4y_aa-2.0E0*1*I_ESP_H4xz_G4y_a;
    abcd[iGrid*1890+1472] = 4.0E0*I_ESP_K4xy2z_G4y_aa-2.0E0*1*I_ESP_H4xy_G4y_a;
    abcd[iGrid*1890+1473] = 4.0E0*I_ESP_K3x3yz_G4y_aa-2.0E0*2*I_ESP_H3xyz_G4y_a;
    abcd[iGrid*1890+1474] = 4.0E0*I_ESP_K3x2y2z_G4y_aa-2.0E0*1*I_ESP_H3x2y_G4y_a-2.0E0*1*I_ESP_H3x2z_G4y_a+1*I_ESP_F3x_G4y;
    abcd[iGrid*1890+1475] = 4.0E0*I_ESP_K3xy3z_G4y_aa-2.0E0*2*I_ESP_H3xyz_G4y_a;
    abcd[iGrid*1890+1476] = 4.0E0*I_ESP_K2x4yz_G4y_aa-2.0E0*3*I_ESP_H2x2yz_G4y_a;
    abcd[iGrid*1890+1477] = 4.0E0*I_ESP_K2x3y2z_G4y_aa-2.0E0*1*I_ESP_H2x3y_G4y_a-2.0E0*2*I_ESP_H2xy2z_G4y_a+2*1*I_ESP_F2xy_G4y;
    abcd[iGrid*1890+1478] = 4.0E0*I_ESP_K2x2y3z_G4y_aa-2.0E0*2*I_ESP_H2x2yz_G4y_a-2.0E0*1*I_ESP_H2x3z_G4y_a+2*I_ESP_F2xz_G4y;
    abcd[iGrid*1890+1479] = 4.0E0*I_ESP_K2xy4z_G4y_aa-2.0E0*3*I_ESP_H2xy2z_G4y_a;
    abcd[iGrid*1890+1480] = 4.0E0*I_ESP_Kx5yz_G4y_aa-2.0E0*4*I_ESP_Hx3yz_G4y_a;
    abcd[iGrid*1890+1481] = 4.0E0*I_ESP_Kx4y2z_G4y_aa-2.0E0*1*I_ESP_Hx4y_G4y_a-2.0E0*3*I_ESP_Hx2y2z_G4y_a+3*1*I_ESP_Fx2y_G4y;
    abcd[iGrid*1890+1482] = 4.0E0*I_ESP_Kx3y3z_G4y_aa-2.0E0*2*I_ESP_Hx3yz_G4y_a-2.0E0*2*I_ESP_Hxy3z_G4y_a+2*2*I_ESP_Fxyz_G4y;
    abcd[iGrid*1890+1483] = 4.0E0*I_ESP_Kx2y4z_G4y_aa-2.0E0*3*I_ESP_Hx2y2z_G4y_a-2.0E0*1*I_ESP_Hx4z_G4y_a+3*I_ESP_Fx2z_G4y;
    abcd[iGrid*1890+1484] = 4.0E0*I_ESP_Kxy5z_G4y_aa-2.0E0*4*I_ESP_Hxy3z_G4y_a;
    abcd[iGrid*1890+1485] = 4.0E0*I_ESP_K6yz_G4y_aa-2.0E0*5*I_ESP_H4yz_G4y_a;
    abcd[iGrid*1890+1486] = 4.0E0*I_ESP_K5y2z_G4y_aa-2.0E0*1*I_ESP_H5y_G4y_a-2.0E0*4*I_ESP_H3y2z_G4y_a+4*1*I_ESP_F3y_G4y;
    abcd[iGrid*1890+1487] = 4.0E0*I_ESP_K4y3z_G4y_aa-2.0E0*2*I_ESP_H4yz_G4y_a-2.0E0*3*I_ESP_H2y3z_G4y_a+3*2*I_ESP_F2yz_G4y;
    abcd[iGrid*1890+1488] = 4.0E0*I_ESP_K3y4z_G4y_aa-2.0E0*3*I_ESP_H3y2z_G4y_a-2.0E0*2*I_ESP_Hy4z_G4y_a+2*3*I_ESP_Fy2z_G4y;
    abcd[iGrid*1890+1489] = 4.0E0*I_ESP_K2y5z_G4y_aa-2.0E0*4*I_ESP_H2y3z_G4y_a-2.0E0*1*I_ESP_H5z_G4y_a+4*I_ESP_F3z_G4y;
    abcd[iGrid*1890+1490] = 4.0E0*I_ESP_Ky6z_G4y_aa-2.0E0*5*I_ESP_Hy4z_G4y_a;
    abcd[iGrid*1890+1491] = 4.0E0*I_ESP_K5xyz_G3yz_aa;
    abcd[iGrid*1890+1492] = 4.0E0*I_ESP_K4x2yz_G3yz_aa-2.0E0*1*I_ESP_H4xz_G3yz_a;
    abcd[iGrid*1890+1493] = 4.0E0*I_ESP_K4xy2z_G3yz_aa-2.0E0*1*I_ESP_H4xy_G3yz_a;
    abcd[iGrid*1890+1494] = 4.0E0*I_ESP_K3x3yz_G3yz_aa-2.0E0*2*I_ESP_H3xyz_G3yz_a;
    abcd[iGrid*1890+1495] = 4.0E0*I_ESP_K3x2y2z_G3yz_aa-2.0E0*1*I_ESP_H3x2y_G3yz_a-2.0E0*1*I_ESP_H3x2z_G3yz_a+1*I_ESP_F3x_G3yz;
    abcd[iGrid*1890+1496] = 4.0E0*I_ESP_K3xy3z_G3yz_aa-2.0E0*2*I_ESP_H3xyz_G3yz_a;
    abcd[iGrid*1890+1497] = 4.0E0*I_ESP_K2x4yz_G3yz_aa-2.0E0*3*I_ESP_H2x2yz_G3yz_a;
    abcd[iGrid*1890+1498] = 4.0E0*I_ESP_K2x3y2z_G3yz_aa-2.0E0*1*I_ESP_H2x3y_G3yz_a-2.0E0*2*I_ESP_H2xy2z_G3yz_a+2*1*I_ESP_F2xy_G3yz;
    abcd[iGrid*1890+1499] = 4.0E0*I_ESP_K2x2y3z_G3yz_aa-2.0E0*2*I_ESP_H2x2yz_G3yz_a-2.0E0*1*I_ESP_H2x3z_G3yz_a+2*I_ESP_F2xz_G3yz;
    abcd[iGrid*1890+1500] = 4.0E0*I_ESP_K2xy4z_G3yz_aa-2.0E0*3*I_ESP_H2xy2z_G3yz_a;
    abcd[iGrid*1890+1501] = 4.0E0*I_ESP_Kx5yz_G3yz_aa-2.0E0*4*I_ESP_Hx3yz_G3yz_a;
    abcd[iGrid*1890+1502] = 4.0E0*I_ESP_Kx4y2z_G3yz_aa-2.0E0*1*I_ESP_Hx4y_G3yz_a-2.0E0*3*I_ESP_Hx2y2z_G3yz_a+3*1*I_ESP_Fx2y_G3yz;
    abcd[iGrid*1890+1503] = 4.0E0*I_ESP_Kx3y3z_G3yz_aa-2.0E0*2*I_ESP_Hx3yz_G3yz_a-2.0E0*2*I_ESP_Hxy3z_G3yz_a+2*2*I_ESP_Fxyz_G3yz;
    abcd[iGrid*1890+1504] = 4.0E0*I_ESP_Kx2y4z_G3yz_aa-2.0E0*3*I_ESP_Hx2y2z_G3yz_a-2.0E0*1*I_ESP_Hx4z_G3yz_a+3*I_ESP_Fx2z_G3yz;
    abcd[iGrid*1890+1505] = 4.0E0*I_ESP_Kxy5z_G3yz_aa-2.0E0*4*I_ESP_Hxy3z_G3yz_a;
    abcd[iGrid*1890+1506] = 4.0E0*I_ESP_K6yz_G3yz_aa-2.0E0*5*I_ESP_H4yz_G3yz_a;
    abcd[iGrid*1890+1507] = 4.0E0*I_ESP_K5y2z_G3yz_aa-2.0E0*1*I_ESP_H5y_G3yz_a-2.0E0*4*I_ESP_H3y2z_G3yz_a+4*1*I_ESP_F3y_G3yz;
    abcd[iGrid*1890+1508] = 4.0E0*I_ESP_K4y3z_G3yz_aa-2.0E0*2*I_ESP_H4yz_G3yz_a-2.0E0*3*I_ESP_H2y3z_G3yz_a+3*2*I_ESP_F2yz_G3yz;
    abcd[iGrid*1890+1509] = 4.0E0*I_ESP_K3y4z_G3yz_aa-2.0E0*3*I_ESP_H3y2z_G3yz_a-2.0E0*2*I_ESP_Hy4z_G3yz_a+2*3*I_ESP_Fy2z_G3yz;
    abcd[iGrid*1890+1510] = 4.0E0*I_ESP_K2y5z_G3yz_aa-2.0E0*4*I_ESP_H2y3z_G3yz_a-2.0E0*1*I_ESP_H5z_G3yz_a+4*I_ESP_F3z_G3yz;
    abcd[iGrid*1890+1511] = 4.0E0*I_ESP_Ky6z_G3yz_aa-2.0E0*5*I_ESP_Hy4z_G3yz_a;
    abcd[iGrid*1890+1512] = 4.0E0*I_ESP_K5xyz_G2y2z_aa;
    abcd[iGrid*1890+1513] = 4.0E0*I_ESP_K4x2yz_G2y2z_aa-2.0E0*1*I_ESP_H4xz_G2y2z_a;
    abcd[iGrid*1890+1514] = 4.0E0*I_ESP_K4xy2z_G2y2z_aa-2.0E0*1*I_ESP_H4xy_G2y2z_a;
    abcd[iGrid*1890+1515] = 4.0E0*I_ESP_K3x3yz_G2y2z_aa-2.0E0*2*I_ESP_H3xyz_G2y2z_a;
    abcd[iGrid*1890+1516] = 4.0E0*I_ESP_K3x2y2z_G2y2z_aa-2.0E0*1*I_ESP_H3x2y_G2y2z_a-2.0E0*1*I_ESP_H3x2z_G2y2z_a+1*I_ESP_F3x_G2y2z;
    abcd[iGrid*1890+1517] = 4.0E0*I_ESP_K3xy3z_G2y2z_aa-2.0E0*2*I_ESP_H3xyz_G2y2z_a;
    abcd[iGrid*1890+1518] = 4.0E0*I_ESP_K2x4yz_G2y2z_aa-2.0E0*3*I_ESP_H2x2yz_G2y2z_a;
    abcd[iGrid*1890+1519] = 4.0E0*I_ESP_K2x3y2z_G2y2z_aa-2.0E0*1*I_ESP_H2x3y_G2y2z_a-2.0E0*2*I_ESP_H2xy2z_G2y2z_a+2*1*I_ESP_F2xy_G2y2z;
    abcd[iGrid*1890+1520] = 4.0E0*I_ESP_K2x2y3z_G2y2z_aa-2.0E0*2*I_ESP_H2x2yz_G2y2z_a-2.0E0*1*I_ESP_H2x3z_G2y2z_a+2*I_ESP_F2xz_G2y2z;
    abcd[iGrid*1890+1521] = 4.0E0*I_ESP_K2xy4z_G2y2z_aa-2.0E0*3*I_ESP_H2xy2z_G2y2z_a;
    abcd[iGrid*1890+1522] = 4.0E0*I_ESP_Kx5yz_G2y2z_aa-2.0E0*4*I_ESP_Hx3yz_G2y2z_a;
    abcd[iGrid*1890+1523] = 4.0E0*I_ESP_Kx4y2z_G2y2z_aa-2.0E0*1*I_ESP_Hx4y_G2y2z_a-2.0E0*3*I_ESP_Hx2y2z_G2y2z_a+3*1*I_ESP_Fx2y_G2y2z;
    abcd[iGrid*1890+1524] = 4.0E0*I_ESP_Kx3y3z_G2y2z_aa-2.0E0*2*I_ESP_Hx3yz_G2y2z_a-2.0E0*2*I_ESP_Hxy3z_G2y2z_a+2*2*I_ESP_Fxyz_G2y2z;
    abcd[iGrid*1890+1525] = 4.0E0*I_ESP_Kx2y4z_G2y2z_aa-2.0E0*3*I_ESP_Hx2y2z_G2y2z_a-2.0E0*1*I_ESP_Hx4z_G2y2z_a+3*I_ESP_Fx2z_G2y2z;
    abcd[iGrid*1890+1526] = 4.0E0*I_ESP_Kxy5z_G2y2z_aa-2.0E0*4*I_ESP_Hxy3z_G2y2z_a;
    abcd[iGrid*1890+1527] = 4.0E0*I_ESP_K6yz_G2y2z_aa-2.0E0*5*I_ESP_H4yz_G2y2z_a;
    abcd[iGrid*1890+1528] = 4.0E0*I_ESP_K5y2z_G2y2z_aa-2.0E0*1*I_ESP_H5y_G2y2z_a-2.0E0*4*I_ESP_H3y2z_G2y2z_a+4*1*I_ESP_F3y_G2y2z;
    abcd[iGrid*1890+1529] = 4.0E0*I_ESP_K4y3z_G2y2z_aa-2.0E0*2*I_ESP_H4yz_G2y2z_a-2.0E0*3*I_ESP_H2y3z_G2y2z_a+3*2*I_ESP_F2yz_G2y2z;
    abcd[iGrid*1890+1530] = 4.0E0*I_ESP_K3y4z_G2y2z_aa-2.0E0*3*I_ESP_H3y2z_G2y2z_a-2.0E0*2*I_ESP_Hy4z_G2y2z_a+2*3*I_ESP_Fy2z_G2y2z;
    abcd[iGrid*1890+1531] = 4.0E0*I_ESP_K2y5z_G2y2z_aa-2.0E0*4*I_ESP_H2y3z_G2y2z_a-2.0E0*1*I_ESP_H5z_G2y2z_a+4*I_ESP_F3z_G2y2z;
    abcd[iGrid*1890+1532] = 4.0E0*I_ESP_Ky6z_G2y2z_aa-2.0E0*5*I_ESP_Hy4z_G2y2z_a;
    abcd[iGrid*1890+1533] = 4.0E0*I_ESP_K5xyz_Gy3z_aa;
    abcd[iGrid*1890+1534] = 4.0E0*I_ESP_K4x2yz_Gy3z_aa-2.0E0*1*I_ESP_H4xz_Gy3z_a;
    abcd[iGrid*1890+1535] = 4.0E0*I_ESP_K4xy2z_Gy3z_aa-2.0E0*1*I_ESP_H4xy_Gy3z_a;
    abcd[iGrid*1890+1536] = 4.0E0*I_ESP_K3x3yz_Gy3z_aa-2.0E0*2*I_ESP_H3xyz_Gy3z_a;
    abcd[iGrid*1890+1537] = 4.0E0*I_ESP_K3x2y2z_Gy3z_aa-2.0E0*1*I_ESP_H3x2y_Gy3z_a-2.0E0*1*I_ESP_H3x2z_Gy3z_a+1*I_ESP_F3x_Gy3z;
    abcd[iGrid*1890+1538] = 4.0E0*I_ESP_K3xy3z_Gy3z_aa-2.0E0*2*I_ESP_H3xyz_Gy3z_a;
    abcd[iGrid*1890+1539] = 4.0E0*I_ESP_K2x4yz_Gy3z_aa-2.0E0*3*I_ESP_H2x2yz_Gy3z_a;
    abcd[iGrid*1890+1540] = 4.0E0*I_ESP_K2x3y2z_Gy3z_aa-2.0E0*1*I_ESP_H2x3y_Gy3z_a-2.0E0*2*I_ESP_H2xy2z_Gy3z_a+2*1*I_ESP_F2xy_Gy3z;
    abcd[iGrid*1890+1541] = 4.0E0*I_ESP_K2x2y3z_Gy3z_aa-2.0E0*2*I_ESP_H2x2yz_Gy3z_a-2.0E0*1*I_ESP_H2x3z_Gy3z_a+2*I_ESP_F2xz_Gy3z;
    abcd[iGrid*1890+1542] = 4.0E0*I_ESP_K2xy4z_Gy3z_aa-2.0E0*3*I_ESP_H2xy2z_Gy3z_a;
    abcd[iGrid*1890+1543] = 4.0E0*I_ESP_Kx5yz_Gy3z_aa-2.0E0*4*I_ESP_Hx3yz_Gy3z_a;
    abcd[iGrid*1890+1544] = 4.0E0*I_ESP_Kx4y2z_Gy3z_aa-2.0E0*1*I_ESP_Hx4y_Gy3z_a-2.0E0*3*I_ESP_Hx2y2z_Gy3z_a+3*1*I_ESP_Fx2y_Gy3z;
    abcd[iGrid*1890+1545] = 4.0E0*I_ESP_Kx3y3z_Gy3z_aa-2.0E0*2*I_ESP_Hx3yz_Gy3z_a-2.0E0*2*I_ESP_Hxy3z_Gy3z_a+2*2*I_ESP_Fxyz_Gy3z;
    abcd[iGrid*1890+1546] = 4.0E0*I_ESP_Kx2y4z_Gy3z_aa-2.0E0*3*I_ESP_Hx2y2z_Gy3z_a-2.0E0*1*I_ESP_Hx4z_Gy3z_a+3*I_ESP_Fx2z_Gy3z;
    abcd[iGrid*1890+1547] = 4.0E0*I_ESP_Kxy5z_Gy3z_aa-2.0E0*4*I_ESP_Hxy3z_Gy3z_a;
    abcd[iGrid*1890+1548] = 4.0E0*I_ESP_K6yz_Gy3z_aa-2.0E0*5*I_ESP_H4yz_Gy3z_a;
    abcd[iGrid*1890+1549] = 4.0E0*I_ESP_K5y2z_Gy3z_aa-2.0E0*1*I_ESP_H5y_Gy3z_a-2.0E0*4*I_ESP_H3y2z_Gy3z_a+4*1*I_ESP_F3y_Gy3z;
    abcd[iGrid*1890+1550] = 4.0E0*I_ESP_K4y3z_Gy3z_aa-2.0E0*2*I_ESP_H4yz_Gy3z_a-2.0E0*3*I_ESP_H2y3z_Gy3z_a+3*2*I_ESP_F2yz_Gy3z;
    abcd[iGrid*1890+1551] = 4.0E0*I_ESP_K3y4z_Gy3z_aa-2.0E0*3*I_ESP_H3y2z_Gy3z_a-2.0E0*2*I_ESP_Hy4z_Gy3z_a+2*3*I_ESP_Fy2z_Gy3z;
    abcd[iGrid*1890+1552] = 4.0E0*I_ESP_K2y5z_Gy3z_aa-2.0E0*4*I_ESP_H2y3z_Gy3z_a-2.0E0*1*I_ESP_H5z_Gy3z_a+4*I_ESP_F3z_Gy3z;
    abcd[iGrid*1890+1553] = 4.0E0*I_ESP_Ky6z_Gy3z_aa-2.0E0*5*I_ESP_Hy4z_Gy3z_a;
    abcd[iGrid*1890+1554] = 4.0E0*I_ESP_K5xyz_G4z_aa;
    abcd[iGrid*1890+1555] = 4.0E0*I_ESP_K4x2yz_G4z_aa-2.0E0*1*I_ESP_H4xz_G4z_a;
    abcd[iGrid*1890+1556] = 4.0E0*I_ESP_K4xy2z_G4z_aa-2.0E0*1*I_ESP_H4xy_G4z_a;
    abcd[iGrid*1890+1557] = 4.0E0*I_ESP_K3x3yz_G4z_aa-2.0E0*2*I_ESP_H3xyz_G4z_a;
    abcd[iGrid*1890+1558] = 4.0E0*I_ESP_K3x2y2z_G4z_aa-2.0E0*1*I_ESP_H3x2y_G4z_a-2.0E0*1*I_ESP_H3x2z_G4z_a+1*I_ESP_F3x_G4z;
    abcd[iGrid*1890+1559] = 4.0E0*I_ESP_K3xy3z_G4z_aa-2.0E0*2*I_ESP_H3xyz_G4z_a;
    abcd[iGrid*1890+1560] = 4.0E0*I_ESP_K2x4yz_G4z_aa-2.0E0*3*I_ESP_H2x2yz_G4z_a;
    abcd[iGrid*1890+1561] = 4.0E0*I_ESP_K2x3y2z_G4z_aa-2.0E0*1*I_ESP_H2x3y_G4z_a-2.0E0*2*I_ESP_H2xy2z_G4z_a+2*1*I_ESP_F2xy_G4z;
    abcd[iGrid*1890+1562] = 4.0E0*I_ESP_K2x2y3z_G4z_aa-2.0E0*2*I_ESP_H2x2yz_G4z_a-2.0E0*1*I_ESP_H2x3z_G4z_a+2*I_ESP_F2xz_G4z;
    abcd[iGrid*1890+1563] = 4.0E0*I_ESP_K2xy4z_G4z_aa-2.0E0*3*I_ESP_H2xy2z_G4z_a;
    abcd[iGrid*1890+1564] = 4.0E0*I_ESP_Kx5yz_G4z_aa-2.0E0*4*I_ESP_Hx3yz_G4z_a;
    abcd[iGrid*1890+1565] = 4.0E0*I_ESP_Kx4y2z_G4z_aa-2.0E0*1*I_ESP_Hx4y_G4z_a-2.0E0*3*I_ESP_Hx2y2z_G4z_a+3*1*I_ESP_Fx2y_G4z;
    abcd[iGrid*1890+1566] = 4.0E0*I_ESP_Kx3y3z_G4z_aa-2.0E0*2*I_ESP_Hx3yz_G4z_a-2.0E0*2*I_ESP_Hxy3z_G4z_a+2*2*I_ESP_Fxyz_G4z;
    abcd[iGrid*1890+1567] = 4.0E0*I_ESP_Kx2y4z_G4z_aa-2.0E0*3*I_ESP_Hx2y2z_G4z_a-2.0E0*1*I_ESP_Hx4z_G4z_a+3*I_ESP_Fx2z_G4z;
    abcd[iGrid*1890+1568] = 4.0E0*I_ESP_Kxy5z_G4z_aa-2.0E0*4*I_ESP_Hxy3z_G4z_a;
    abcd[iGrid*1890+1569] = 4.0E0*I_ESP_K6yz_G4z_aa-2.0E0*5*I_ESP_H4yz_G4z_a;
    abcd[iGrid*1890+1570] = 4.0E0*I_ESP_K5y2z_G4z_aa-2.0E0*1*I_ESP_H5y_G4z_a-2.0E0*4*I_ESP_H3y2z_G4z_a+4*1*I_ESP_F3y_G4z;
    abcd[iGrid*1890+1571] = 4.0E0*I_ESP_K4y3z_G4z_aa-2.0E0*2*I_ESP_H4yz_G4z_a-2.0E0*3*I_ESP_H2y3z_G4z_a+3*2*I_ESP_F2yz_G4z;
    abcd[iGrid*1890+1572] = 4.0E0*I_ESP_K3y4z_G4z_aa-2.0E0*3*I_ESP_H3y2z_G4z_a-2.0E0*2*I_ESP_Hy4z_G4z_a+2*3*I_ESP_Fy2z_G4z;
    abcd[iGrid*1890+1573] = 4.0E0*I_ESP_K2y5z_G4z_aa-2.0E0*4*I_ESP_H2y3z_G4z_a-2.0E0*1*I_ESP_H5z_G4z_a+4*I_ESP_F3z_G4z;
    abcd[iGrid*1890+1574] = 4.0E0*I_ESP_Ky6z_G4z_aa-2.0E0*5*I_ESP_Hy4z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_G_aa
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_H_G_a
     * RHS shell quartet name: SQ_ESP_F_G
     ************************************************************/
    abcd[iGrid*1890+1575] = 4.0E0*I_ESP_K5x2z_G4x_aa-2.0E0*1*I_ESP_H5x_G4x_a;
    abcd[iGrid*1890+1576] = 4.0E0*I_ESP_K4xy2z_G4x_aa-2.0E0*1*I_ESP_H4xy_G4x_a;
    abcd[iGrid*1890+1577] = 4.0E0*I_ESP_K4x3z_G4x_aa-2.0E0*1*I_ESP_H4xz_G4x_a-2.0E0*2*I_ESP_H4xz_G4x_a;
    abcd[iGrid*1890+1578] = 4.0E0*I_ESP_K3x2y2z_G4x_aa-2.0E0*1*I_ESP_H3x2y_G4x_a;
    abcd[iGrid*1890+1579] = 4.0E0*I_ESP_K3xy3z_G4x_aa-2.0E0*1*I_ESP_H3xyz_G4x_a-2.0E0*2*I_ESP_H3xyz_G4x_a;
    abcd[iGrid*1890+1580] = 4.0E0*I_ESP_K3x4z_G4x_aa-2.0E0*2*I_ESP_H3x2z_G4x_a-2.0E0*3*I_ESP_H3x2z_G4x_a+2*1*I_ESP_F3x_G4x;
    abcd[iGrid*1890+1581] = 4.0E0*I_ESP_K2x3y2z_G4x_aa-2.0E0*1*I_ESP_H2x3y_G4x_a;
    abcd[iGrid*1890+1582] = 4.0E0*I_ESP_K2x2y3z_G4x_aa-2.0E0*1*I_ESP_H2x2yz_G4x_a-2.0E0*2*I_ESP_H2x2yz_G4x_a;
    abcd[iGrid*1890+1583] = 4.0E0*I_ESP_K2xy4z_G4x_aa-2.0E0*2*I_ESP_H2xy2z_G4x_a-2.0E0*3*I_ESP_H2xy2z_G4x_a+2*1*I_ESP_F2xy_G4x;
    abcd[iGrid*1890+1584] = 4.0E0*I_ESP_K2x5z_G4x_aa-2.0E0*3*I_ESP_H2x3z_G4x_a-2.0E0*4*I_ESP_H2x3z_G4x_a+3*2*I_ESP_F2xz_G4x;
    abcd[iGrid*1890+1585] = 4.0E0*I_ESP_Kx4y2z_G4x_aa-2.0E0*1*I_ESP_Hx4y_G4x_a;
    abcd[iGrid*1890+1586] = 4.0E0*I_ESP_Kx3y3z_G4x_aa-2.0E0*1*I_ESP_Hx3yz_G4x_a-2.0E0*2*I_ESP_Hx3yz_G4x_a;
    abcd[iGrid*1890+1587] = 4.0E0*I_ESP_Kx2y4z_G4x_aa-2.0E0*2*I_ESP_Hx2y2z_G4x_a-2.0E0*3*I_ESP_Hx2y2z_G4x_a+2*1*I_ESP_Fx2y_G4x;
    abcd[iGrid*1890+1588] = 4.0E0*I_ESP_Kxy5z_G4x_aa-2.0E0*3*I_ESP_Hxy3z_G4x_a-2.0E0*4*I_ESP_Hxy3z_G4x_a+3*2*I_ESP_Fxyz_G4x;
    abcd[iGrid*1890+1589] = 4.0E0*I_ESP_Kx6z_G4x_aa-2.0E0*4*I_ESP_Hx4z_G4x_a-2.0E0*5*I_ESP_Hx4z_G4x_a+4*3*I_ESP_Fx2z_G4x;
    abcd[iGrid*1890+1590] = 4.0E0*I_ESP_K5y2z_G4x_aa-2.0E0*1*I_ESP_H5y_G4x_a;
    abcd[iGrid*1890+1591] = 4.0E0*I_ESP_K4y3z_G4x_aa-2.0E0*1*I_ESP_H4yz_G4x_a-2.0E0*2*I_ESP_H4yz_G4x_a;
    abcd[iGrid*1890+1592] = 4.0E0*I_ESP_K3y4z_G4x_aa-2.0E0*2*I_ESP_H3y2z_G4x_a-2.0E0*3*I_ESP_H3y2z_G4x_a+2*1*I_ESP_F3y_G4x;
    abcd[iGrid*1890+1593] = 4.0E0*I_ESP_K2y5z_G4x_aa-2.0E0*3*I_ESP_H2y3z_G4x_a-2.0E0*4*I_ESP_H2y3z_G4x_a+3*2*I_ESP_F2yz_G4x;
    abcd[iGrid*1890+1594] = 4.0E0*I_ESP_Ky6z_G4x_aa-2.0E0*4*I_ESP_Hy4z_G4x_a-2.0E0*5*I_ESP_Hy4z_G4x_a+4*3*I_ESP_Fy2z_G4x;
    abcd[iGrid*1890+1595] = 4.0E0*I_ESP_K7z_G4x_aa-2.0E0*5*I_ESP_H5z_G4x_a-2.0E0*6*I_ESP_H5z_G4x_a+5*4*I_ESP_F3z_G4x;
    abcd[iGrid*1890+1596] = 4.0E0*I_ESP_K5x2z_G3xy_aa-2.0E0*1*I_ESP_H5x_G3xy_a;
    abcd[iGrid*1890+1597] = 4.0E0*I_ESP_K4xy2z_G3xy_aa-2.0E0*1*I_ESP_H4xy_G3xy_a;
    abcd[iGrid*1890+1598] = 4.0E0*I_ESP_K4x3z_G3xy_aa-2.0E0*1*I_ESP_H4xz_G3xy_a-2.0E0*2*I_ESP_H4xz_G3xy_a;
    abcd[iGrid*1890+1599] = 4.0E0*I_ESP_K3x2y2z_G3xy_aa-2.0E0*1*I_ESP_H3x2y_G3xy_a;
    abcd[iGrid*1890+1600] = 4.0E0*I_ESP_K3xy3z_G3xy_aa-2.0E0*1*I_ESP_H3xyz_G3xy_a-2.0E0*2*I_ESP_H3xyz_G3xy_a;
    abcd[iGrid*1890+1601] = 4.0E0*I_ESP_K3x4z_G3xy_aa-2.0E0*2*I_ESP_H3x2z_G3xy_a-2.0E0*3*I_ESP_H3x2z_G3xy_a+2*1*I_ESP_F3x_G3xy;
    abcd[iGrid*1890+1602] = 4.0E0*I_ESP_K2x3y2z_G3xy_aa-2.0E0*1*I_ESP_H2x3y_G3xy_a;
    abcd[iGrid*1890+1603] = 4.0E0*I_ESP_K2x2y3z_G3xy_aa-2.0E0*1*I_ESP_H2x2yz_G3xy_a-2.0E0*2*I_ESP_H2x2yz_G3xy_a;
    abcd[iGrid*1890+1604] = 4.0E0*I_ESP_K2xy4z_G3xy_aa-2.0E0*2*I_ESP_H2xy2z_G3xy_a-2.0E0*3*I_ESP_H2xy2z_G3xy_a+2*1*I_ESP_F2xy_G3xy;
    abcd[iGrid*1890+1605] = 4.0E0*I_ESP_K2x5z_G3xy_aa-2.0E0*3*I_ESP_H2x3z_G3xy_a-2.0E0*4*I_ESP_H2x3z_G3xy_a+3*2*I_ESP_F2xz_G3xy;
    abcd[iGrid*1890+1606] = 4.0E0*I_ESP_Kx4y2z_G3xy_aa-2.0E0*1*I_ESP_Hx4y_G3xy_a;
    abcd[iGrid*1890+1607] = 4.0E0*I_ESP_Kx3y3z_G3xy_aa-2.0E0*1*I_ESP_Hx3yz_G3xy_a-2.0E0*2*I_ESP_Hx3yz_G3xy_a;
    abcd[iGrid*1890+1608] = 4.0E0*I_ESP_Kx2y4z_G3xy_aa-2.0E0*2*I_ESP_Hx2y2z_G3xy_a-2.0E0*3*I_ESP_Hx2y2z_G3xy_a+2*1*I_ESP_Fx2y_G3xy;
    abcd[iGrid*1890+1609] = 4.0E0*I_ESP_Kxy5z_G3xy_aa-2.0E0*3*I_ESP_Hxy3z_G3xy_a-2.0E0*4*I_ESP_Hxy3z_G3xy_a+3*2*I_ESP_Fxyz_G3xy;
    abcd[iGrid*1890+1610] = 4.0E0*I_ESP_Kx6z_G3xy_aa-2.0E0*4*I_ESP_Hx4z_G3xy_a-2.0E0*5*I_ESP_Hx4z_G3xy_a+4*3*I_ESP_Fx2z_G3xy;
    abcd[iGrid*1890+1611] = 4.0E0*I_ESP_K5y2z_G3xy_aa-2.0E0*1*I_ESP_H5y_G3xy_a;
    abcd[iGrid*1890+1612] = 4.0E0*I_ESP_K4y3z_G3xy_aa-2.0E0*1*I_ESP_H4yz_G3xy_a-2.0E0*2*I_ESP_H4yz_G3xy_a;
    abcd[iGrid*1890+1613] = 4.0E0*I_ESP_K3y4z_G3xy_aa-2.0E0*2*I_ESP_H3y2z_G3xy_a-2.0E0*3*I_ESP_H3y2z_G3xy_a+2*1*I_ESP_F3y_G3xy;
    abcd[iGrid*1890+1614] = 4.0E0*I_ESP_K2y5z_G3xy_aa-2.0E0*3*I_ESP_H2y3z_G3xy_a-2.0E0*4*I_ESP_H2y3z_G3xy_a+3*2*I_ESP_F2yz_G3xy;
    abcd[iGrid*1890+1615] = 4.0E0*I_ESP_Ky6z_G3xy_aa-2.0E0*4*I_ESP_Hy4z_G3xy_a-2.0E0*5*I_ESP_Hy4z_G3xy_a+4*3*I_ESP_Fy2z_G3xy;
    abcd[iGrid*1890+1616] = 4.0E0*I_ESP_K7z_G3xy_aa-2.0E0*5*I_ESP_H5z_G3xy_a-2.0E0*6*I_ESP_H5z_G3xy_a+5*4*I_ESP_F3z_G3xy;
    abcd[iGrid*1890+1617] = 4.0E0*I_ESP_K5x2z_G3xz_aa-2.0E0*1*I_ESP_H5x_G3xz_a;
    abcd[iGrid*1890+1618] = 4.0E0*I_ESP_K4xy2z_G3xz_aa-2.0E0*1*I_ESP_H4xy_G3xz_a;
    abcd[iGrid*1890+1619] = 4.0E0*I_ESP_K4x3z_G3xz_aa-2.0E0*1*I_ESP_H4xz_G3xz_a-2.0E0*2*I_ESP_H4xz_G3xz_a;
    abcd[iGrid*1890+1620] = 4.0E0*I_ESP_K3x2y2z_G3xz_aa-2.0E0*1*I_ESP_H3x2y_G3xz_a;
    abcd[iGrid*1890+1621] = 4.0E0*I_ESP_K3xy3z_G3xz_aa-2.0E0*1*I_ESP_H3xyz_G3xz_a-2.0E0*2*I_ESP_H3xyz_G3xz_a;
    abcd[iGrid*1890+1622] = 4.0E0*I_ESP_K3x4z_G3xz_aa-2.0E0*2*I_ESP_H3x2z_G3xz_a-2.0E0*3*I_ESP_H3x2z_G3xz_a+2*1*I_ESP_F3x_G3xz;
    abcd[iGrid*1890+1623] = 4.0E0*I_ESP_K2x3y2z_G3xz_aa-2.0E0*1*I_ESP_H2x3y_G3xz_a;
    abcd[iGrid*1890+1624] = 4.0E0*I_ESP_K2x2y3z_G3xz_aa-2.0E0*1*I_ESP_H2x2yz_G3xz_a-2.0E0*2*I_ESP_H2x2yz_G3xz_a;
    abcd[iGrid*1890+1625] = 4.0E0*I_ESP_K2xy4z_G3xz_aa-2.0E0*2*I_ESP_H2xy2z_G3xz_a-2.0E0*3*I_ESP_H2xy2z_G3xz_a+2*1*I_ESP_F2xy_G3xz;
    abcd[iGrid*1890+1626] = 4.0E0*I_ESP_K2x5z_G3xz_aa-2.0E0*3*I_ESP_H2x3z_G3xz_a-2.0E0*4*I_ESP_H2x3z_G3xz_a+3*2*I_ESP_F2xz_G3xz;
    abcd[iGrid*1890+1627] = 4.0E0*I_ESP_Kx4y2z_G3xz_aa-2.0E0*1*I_ESP_Hx4y_G3xz_a;
    abcd[iGrid*1890+1628] = 4.0E0*I_ESP_Kx3y3z_G3xz_aa-2.0E0*1*I_ESP_Hx3yz_G3xz_a-2.0E0*2*I_ESP_Hx3yz_G3xz_a;
    abcd[iGrid*1890+1629] = 4.0E0*I_ESP_Kx2y4z_G3xz_aa-2.0E0*2*I_ESP_Hx2y2z_G3xz_a-2.0E0*3*I_ESP_Hx2y2z_G3xz_a+2*1*I_ESP_Fx2y_G3xz;
    abcd[iGrid*1890+1630] = 4.0E0*I_ESP_Kxy5z_G3xz_aa-2.0E0*3*I_ESP_Hxy3z_G3xz_a-2.0E0*4*I_ESP_Hxy3z_G3xz_a+3*2*I_ESP_Fxyz_G3xz;
    abcd[iGrid*1890+1631] = 4.0E0*I_ESP_Kx6z_G3xz_aa-2.0E0*4*I_ESP_Hx4z_G3xz_a-2.0E0*5*I_ESP_Hx4z_G3xz_a+4*3*I_ESP_Fx2z_G3xz;
    abcd[iGrid*1890+1632] = 4.0E0*I_ESP_K5y2z_G3xz_aa-2.0E0*1*I_ESP_H5y_G3xz_a;
    abcd[iGrid*1890+1633] = 4.0E0*I_ESP_K4y3z_G3xz_aa-2.0E0*1*I_ESP_H4yz_G3xz_a-2.0E0*2*I_ESP_H4yz_G3xz_a;
    abcd[iGrid*1890+1634] = 4.0E0*I_ESP_K3y4z_G3xz_aa-2.0E0*2*I_ESP_H3y2z_G3xz_a-2.0E0*3*I_ESP_H3y2z_G3xz_a+2*1*I_ESP_F3y_G3xz;
    abcd[iGrid*1890+1635] = 4.0E0*I_ESP_K2y5z_G3xz_aa-2.0E0*3*I_ESP_H2y3z_G3xz_a-2.0E0*4*I_ESP_H2y3z_G3xz_a+3*2*I_ESP_F2yz_G3xz;
    abcd[iGrid*1890+1636] = 4.0E0*I_ESP_Ky6z_G3xz_aa-2.0E0*4*I_ESP_Hy4z_G3xz_a-2.0E0*5*I_ESP_Hy4z_G3xz_a+4*3*I_ESP_Fy2z_G3xz;
    abcd[iGrid*1890+1637] = 4.0E0*I_ESP_K7z_G3xz_aa-2.0E0*5*I_ESP_H5z_G3xz_a-2.0E0*6*I_ESP_H5z_G3xz_a+5*4*I_ESP_F3z_G3xz;
    abcd[iGrid*1890+1638] = 4.0E0*I_ESP_K5x2z_G2x2y_aa-2.0E0*1*I_ESP_H5x_G2x2y_a;
    abcd[iGrid*1890+1639] = 4.0E0*I_ESP_K4xy2z_G2x2y_aa-2.0E0*1*I_ESP_H4xy_G2x2y_a;
    abcd[iGrid*1890+1640] = 4.0E0*I_ESP_K4x3z_G2x2y_aa-2.0E0*1*I_ESP_H4xz_G2x2y_a-2.0E0*2*I_ESP_H4xz_G2x2y_a;
    abcd[iGrid*1890+1641] = 4.0E0*I_ESP_K3x2y2z_G2x2y_aa-2.0E0*1*I_ESP_H3x2y_G2x2y_a;
    abcd[iGrid*1890+1642] = 4.0E0*I_ESP_K3xy3z_G2x2y_aa-2.0E0*1*I_ESP_H3xyz_G2x2y_a-2.0E0*2*I_ESP_H3xyz_G2x2y_a;
    abcd[iGrid*1890+1643] = 4.0E0*I_ESP_K3x4z_G2x2y_aa-2.0E0*2*I_ESP_H3x2z_G2x2y_a-2.0E0*3*I_ESP_H3x2z_G2x2y_a+2*1*I_ESP_F3x_G2x2y;
    abcd[iGrid*1890+1644] = 4.0E0*I_ESP_K2x3y2z_G2x2y_aa-2.0E0*1*I_ESP_H2x3y_G2x2y_a;
    abcd[iGrid*1890+1645] = 4.0E0*I_ESP_K2x2y3z_G2x2y_aa-2.0E0*1*I_ESP_H2x2yz_G2x2y_a-2.0E0*2*I_ESP_H2x2yz_G2x2y_a;
    abcd[iGrid*1890+1646] = 4.0E0*I_ESP_K2xy4z_G2x2y_aa-2.0E0*2*I_ESP_H2xy2z_G2x2y_a-2.0E0*3*I_ESP_H2xy2z_G2x2y_a+2*1*I_ESP_F2xy_G2x2y;
    abcd[iGrid*1890+1647] = 4.0E0*I_ESP_K2x5z_G2x2y_aa-2.0E0*3*I_ESP_H2x3z_G2x2y_a-2.0E0*4*I_ESP_H2x3z_G2x2y_a+3*2*I_ESP_F2xz_G2x2y;
    abcd[iGrid*1890+1648] = 4.0E0*I_ESP_Kx4y2z_G2x2y_aa-2.0E0*1*I_ESP_Hx4y_G2x2y_a;
    abcd[iGrid*1890+1649] = 4.0E0*I_ESP_Kx3y3z_G2x2y_aa-2.0E0*1*I_ESP_Hx3yz_G2x2y_a-2.0E0*2*I_ESP_Hx3yz_G2x2y_a;
    abcd[iGrid*1890+1650] = 4.0E0*I_ESP_Kx2y4z_G2x2y_aa-2.0E0*2*I_ESP_Hx2y2z_G2x2y_a-2.0E0*3*I_ESP_Hx2y2z_G2x2y_a+2*1*I_ESP_Fx2y_G2x2y;
    abcd[iGrid*1890+1651] = 4.0E0*I_ESP_Kxy5z_G2x2y_aa-2.0E0*3*I_ESP_Hxy3z_G2x2y_a-2.0E0*4*I_ESP_Hxy3z_G2x2y_a+3*2*I_ESP_Fxyz_G2x2y;
    abcd[iGrid*1890+1652] = 4.0E0*I_ESP_Kx6z_G2x2y_aa-2.0E0*4*I_ESP_Hx4z_G2x2y_a-2.0E0*5*I_ESP_Hx4z_G2x2y_a+4*3*I_ESP_Fx2z_G2x2y;
    abcd[iGrid*1890+1653] = 4.0E0*I_ESP_K5y2z_G2x2y_aa-2.0E0*1*I_ESP_H5y_G2x2y_a;
    abcd[iGrid*1890+1654] = 4.0E0*I_ESP_K4y3z_G2x2y_aa-2.0E0*1*I_ESP_H4yz_G2x2y_a-2.0E0*2*I_ESP_H4yz_G2x2y_a;
    abcd[iGrid*1890+1655] = 4.0E0*I_ESP_K3y4z_G2x2y_aa-2.0E0*2*I_ESP_H3y2z_G2x2y_a-2.0E0*3*I_ESP_H3y2z_G2x2y_a+2*1*I_ESP_F3y_G2x2y;
    abcd[iGrid*1890+1656] = 4.0E0*I_ESP_K2y5z_G2x2y_aa-2.0E0*3*I_ESP_H2y3z_G2x2y_a-2.0E0*4*I_ESP_H2y3z_G2x2y_a+3*2*I_ESP_F2yz_G2x2y;
    abcd[iGrid*1890+1657] = 4.0E0*I_ESP_Ky6z_G2x2y_aa-2.0E0*4*I_ESP_Hy4z_G2x2y_a-2.0E0*5*I_ESP_Hy4z_G2x2y_a+4*3*I_ESP_Fy2z_G2x2y;
    abcd[iGrid*1890+1658] = 4.0E0*I_ESP_K7z_G2x2y_aa-2.0E0*5*I_ESP_H5z_G2x2y_a-2.0E0*6*I_ESP_H5z_G2x2y_a+5*4*I_ESP_F3z_G2x2y;
    abcd[iGrid*1890+1659] = 4.0E0*I_ESP_K5x2z_G2xyz_aa-2.0E0*1*I_ESP_H5x_G2xyz_a;
    abcd[iGrid*1890+1660] = 4.0E0*I_ESP_K4xy2z_G2xyz_aa-2.0E0*1*I_ESP_H4xy_G2xyz_a;
    abcd[iGrid*1890+1661] = 4.0E0*I_ESP_K4x3z_G2xyz_aa-2.0E0*1*I_ESP_H4xz_G2xyz_a-2.0E0*2*I_ESP_H4xz_G2xyz_a;
    abcd[iGrid*1890+1662] = 4.0E0*I_ESP_K3x2y2z_G2xyz_aa-2.0E0*1*I_ESP_H3x2y_G2xyz_a;
    abcd[iGrid*1890+1663] = 4.0E0*I_ESP_K3xy3z_G2xyz_aa-2.0E0*1*I_ESP_H3xyz_G2xyz_a-2.0E0*2*I_ESP_H3xyz_G2xyz_a;
    abcd[iGrid*1890+1664] = 4.0E0*I_ESP_K3x4z_G2xyz_aa-2.0E0*2*I_ESP_H3x2z_G2xyz_a-2.0E0*3*I_ESP_H3x2z_G2xyz_a+2*1*I_ESP_F3x_G2xyz;
    abcd[iGrid*1890+1665] = 4.0E0*I_ESP_K2x3y2z_G2xyz_aa-2.0E0*1*I_ESP_H2x3y_G2xyz_a;
    abcd[iGrid*1890+1666] = 4.0E0*I_ESP_K2x2y3z_G2xyz_aa-2.0E0*1*I_ESP_H2x2yz_G2xyz_a-2.0E0*2*I_ESP_H2x2yz_G2xyz_a;
    abcd[iGrid*1890+1667] = 4.0E0*I_ESP_K2xy4z_G2xyz_aa-2.0E0*2*I_ESP_H2xy2z_G2xyz_a-2.0E0*3*I_ESP_H2xy2z_G2xyz_a+2*1*I_ESP_F2xy_G2xyz;
    abcd[iGrid*1890+1668] = 4.0E0*I_ESP_K2x5z_G2xyz_aa-2.0E0*3*I_ESP_H2x3z_G2xyz_a-2.0E0*4*I_ESP_H2x3z_G2xyz_a+3*2*I_ESP_F2xz_G2xyz;
    abcd[iGrid*1890+1669] = 4.0E0*I_ESP_Kx4y2z_G2xyz_aa-2.0E0*1*I_ESP_Hx4y_G2xyz_a;
    abcd[iGrid*1890+1670] = 4.0E0*I_ESP_Kx3y3z_G2xyz_aa-2.0E0*1*I_ESP_Hx3yz_G2xyz_a-2.0E0*2*I_ESP_Hx3yz_G2xyz_a;
    abcd[iGrid*1890+1671] = 4.0E0*I_ESP_Kx2y4z_G2xyz_aa-2.0E0*2*I_ESP_Hx2y2z_G2xyz_a-2.0E0*3*I_ESP_Hx2y2z_G2xyz_a+2*1*I_ESP_Fx2y_G2xyz;
    abcd[iGrid*1890+1672] = 4.0E0*I_ESP_Kxy5z_G2xyz_aa-2.0E0*3*I_ESP_Hxy3z_G2xyz_a-2.0E0*4*I_ESP_Hxy3z_G2xyz_a+3*2*I_ESP_Fxyz_G2xyz;
    abcd[iGrid*1890+1673] = 4.0E0*I_ESP_Kx6z_G2xyz_aa-2.0E0*4*I_ESP_Hx4z_G2xyz_a-2.0E0*5*I_ESP_Hx4z_G2xyz_a+4*3*I_ESP_Fx2z_G2xyz;
    abcd[iGrid*1890+1674] = 4.0E0*I_ESP_K5y2z_G2xyz_aa-2.0E0*1*I_ESP_H5y_G2xyz_a;
    abcd[iGrid*1890+1675] = 4.0E0*I_ESP_K4y3z_G2xyz_aa-2.0E0*1*I_ESP_H4yz_G2xyz_a-2.0E0*2*I_ESP_H4yz_G2xyz_a;
    abcd[iGrid*1890+1676] = 4.0E0*I_ESP_K3y4z_G2xyz_aa-2.0E0*2*I_ESP_H3y2z_G2xyz_a-2.0E0*3*I_ESP_H3y2z_G2xyz_a+2*1*I_ESP_F3y_G2xyz;
    abcd[iGrid*1890+1677] = 4.0E0*I_ESP_K2y5z_G2xyz_aa-2.0E0*3*I_ESP_H2y3z_G2xyz_a-2.0E0*4*I_ESP_H2y3z_G2xyz_a+3*2*I_ESP_F2yz_G2xyz;
    abcd[iGrid*1890+1678] = 4.0E0*I_ESP_Ky6z_G2xyz_aa-2.0E0*4*I_ESP_Hy4z_G2xyz_a-2.0E0*5*I_ESP_Hy4z_G2xyz_a+4*3*I_ESP_Fy2z_G2xyz;
    abcd[iGrid*1890+1679] = 4.0E0*I_ESP_K7z_G2xyz_aa-2.0E0*5*I_ESP_H5z_G2xyz_a-2.0E0*6*I_ESP_H5z_G2xyz_a+5*4*I_ESP_F3z_G2xyz;
    abcd[iGrid*1890+1680] = 4.0E0*I_ESP_K5x2z_G2x2z_aa-2.0E0*1*I_ESP_H5x_G2x2z_a;
    abcd[iGrid*1890+1681] = 4.0E0*I_ESP_K4xy2z_G2x2z_aa-2.0E0*1*I_ESP_H4xy_G2x2z_a;
    abcd[iGrid*1890+1682] = 4.0E0*I_ESP_K4x3z_G2x2z_aa-2.0E0*1*I_ESP_H4xz_G2x2z_a-2.0E0*2*I_ESP_H4xz_G2x2z_a;
    abcd[iGrid*1890+1683] = 4.0E0*I_ESP_K3x2y2z_G2x2z_aa-2.0E0*1*I_ESP_H3x2y_G2x2z_a;
    abcd[iGrid*1890+1684] = 4.0E0*I_ESP_K3xy3z_G2x2z_aa-2.0E0*1*I_ESP_H3xyz_G2x2z_a-2.0E0*2*I_ESP_H3xyz_G2x2z_a;
    abcd[iGrid*1890+1685] = 4.0E0*I_ESP_K3x4z_G2x2z_aa-2.0E0*2*I_ESP_H3x2z_G2x2z_a-2.0E0*3*I_ESP_H3x2z_G2x2z_a+2*1*I_ESP_F3x_G2x2z;
    abcd[iGrid*1890+1686] = 4.0E0*I_ESP_K2x3y2z_G2x2z_aa-2.0E0*1*I_ESP_H2x3y_G2x2z_a;
    abcd[iGrid*1890+1687] = 4.0E0*I_ESP_K2x2y3z_G2x2z_aa-2.0E0*1*I_ESP_H2x2yz_G2x2z_a-2.0E0*2*I_ESP_H2x2yz_G2x2z_a;
    abcd[iGrid*1890+1688] = 4.0E0*I_ESP_K2xy4z_G2x2z_aa-2.0E0*2*I_ESP_H2xy2z_G2x2z_a-2.0E0*3*I_ESP_H2xy2z_G2x2z_a+2*1*I_ESP_F2xy_G2x2z;
    abcd[iGrid*1890+1689] = 4.0E0*I_ESP_K2x5z_G2x2z_aa-2.0E0*3*I_ESP_H2x3z_G2x2z_a-2.0E0*4*I_ESP_H2x3z_G2x2z_a+3*2*I_ESP_F2xz_G2x2z;
    abcd[iGrid*1890+1690] = 4.0E0*I_ESP_Kx4y2z_G2x2z_aa-2.0E0*1*I_ESP_Hx4y_G2x2z_a;
    abcd[iGrid*1890+1691] = 4.0E0*I_ESP_Kx3y3z_G2x2z_aa-2.0E0*1*I_ESP_Hx3yz_G2x2z_a-2.0E0*2*I_ESP_Hx3yz_G2x2z_a;
    abcd[iGrid*1890+1692] = 4.0E0*I_ESP_Kx2y4z_G2x2z_aa-2.0E0*2*I_ESP_Hx2y2z_G2x2z_a-2.0E0*3*I_ESP_Hx2y2z_G2x2z_a+2*1*I_ESP_Fx2y_G2x2z;
    abcd[iGrid*1890+1693] = 4.0E0*I_ESP_Kxy5z_G2x2z_aa-2.0E0*3*I_ESP_Hxy3z_G2x2z_a-2.0E0*4*I_ESP_Hxy3z_G2x2z_a+3*2*I_ESP_Fxyz_G2x2z;
    abcd[iGrid*1890+1694] = 4.0E0*I_ESP_Kx6z_G2x2z_aa-2.0E0*4*I_ESP_Hx4z_G2x2z_a-2.0E0*5*I_ESP_Hx4z_G2x2z_a+4*3*I_ESP_Fx2z_G2x2z;
    abcd[iGrid*1890+1695] = 4.0E0*I_ESP_K5y2z_G2x2z_aa-2.0E0*1*I_ESP_H5y_G2x2z_a;
    abcd[iGrid*1890+1696] = 4.0E0*I_ESP_K4y3z_G2x2z_aa-2.0E0*1*I_ESP_H4yz_G2x2z_a-2.0E0*2*I_ESP_H4yz_G2x2z_a;
    abcd[iGrid*1890+1697] = 4.0E0*I_ESP_K3y4z_G2x2z_aa-2.0E0*2*I_ESP_H3y2z_G2x2z_a-2.0E0*3*I_ESP_H3y2z_G2x2z_a+2*1*I_ESP_F3y_G2x2z;
    abcd[iGrid*1890+1698] = 4.0E0*I_ESP_K2y5z_G2x2z_aa-2.0E0*3*I_ESP_H2y3z_G2x2z_a-2.0E0*4*I_ESP_H2y3z_G2x2z_a+3*2*I_ESP_F2yz_G2x2z;
    abcd[iGrid*1890+1699] = 4.0E0*I_ESP_Ky6z_G2x2z_aa-2.0E0*4*I_ESP_Hy4z_G2x2z_a-2.0E0*5*I_ESP_Hy4z_G2x2z_a+4*3*I_ESP_Fy2z_G2x2z;
    abcd[iGrid*1890+1700] = 4.0E0*I_ESP_K7z_G2x2z_aa-2.0E0*5*I_ESP_H5z_G2x2z_a-2.0E0*6*I_ESP_H5z_G2x2z_a+5*4*I_ESP_F3z_G2x2z;
    abcd[iGrid*1890+1701] = 4.0E0*I_ESP_K5x2z_Gx3y_aa-2.0E0*1*I_ESP_H5x_Gx3y_a;
    abcd[iGrid*1890+1702] = 4.0E0*I_ESP_K4xy2z_Gx3y_aa-2.0E0*1*I_ESP_H4xy_Gx3y_a;
    abcd[iGrid*1890+1703] = 4.0E0*I_ESP_K4x3z_Gx3y_aa-2.0E0*1*I_ESP_H4xz_Gx3y_a-2.0E0*2*I_ESP_H4xz_Gx3y_a;
    abcd[iGrid*1890+1704] = 4.0E0*I_ESP_K3x2y2z_Gx3y_aa-2.0E0*1*I_ESP_H3x2y_Gx3y_a;
    abcd[iGrid*1890+1705] = 4.0E0*I_ESP_K3xy3z_Gx3y_aa-2.0E0*1*I_ESP_H3xyz_Gx3y_a-2.0E0*2*I_ESP_H3xyz_Gx3y_a;
    abcd[iGrid*1890+1706] = 4.0E0*I_ESP_K3x4z_Gx3y_aa-2.0E0*2*I_ESP_H3x2z_Gx3y_a-2.0E0*3*I_ESP_H3x2z_Gx3y_a+2*1*I_ESP_F3x_Gx3y;
    abcd[iGrid*1890+1707] = 4.0E0*I_ESP_K2x3y2z_Gx3y_aa-2.0E0*1*I_ESP_H2x3y_Gx3y_a;
    abcd[iGrid*1890+1708] = 4.0E0*I_ESP_K2x2y3z_Gx3y_aa-2.0E0*1*I_ESP_H2x2yz_Gx3y_a-2.0E0*2*I_ESP_H2x2yz_Gx3y_a;
    abcd[iGrid*1890+1709] = 4.0E0*I_ESP_K2xy4z_Gx3y_aa-2.0E0*2*I_ESP_H2xy2z_Gx3y_a-2.0E0*3*I_ESP_H2xy2z_Gx3y_a+2*1*I_ESP_F2xy_Gx3y;
    abcd[iGrid*1890+1710] = 4.0E0*I_ESP_K2x5z_Gx3y_aa-2.0E0*3*I_ESP_H2x3z_Gx3y_a-2.0E0*4*I_ESP_H2x3z_Gx3y_a+3*2*I_ESP_F2xz_Gx3y;
    abcd[iGrid*1890+1711] = 4.0E0*I_ESP_Kx4y2z_Gx3y_aa-2.0E0*1*I_ESP_Hx4y_Gx3y_a;
    abcd[iGrid*1890+1712] = 4.0E0*I_ESP_Kx3y3z_Gx3y_aa-2.0E0*1*I_ESP_Hx3yz_Gx3y_a-2.0E0*2*I_ESP_Hx3yz_Gx3y_a;
    abcd[iGrid*1890+1713] = 4.0E0*I_ESP_Kx2y4z_Gx3y_aa-2.0E0*2*I_ESP_Hx2y2z_Gx3y_a-2.0E0*3*I_ESP_Hx2y2z_Gx3y_a+2*1*I_ESP_Fx2y_Gx3y;
    abcd[iGrid*1890+1714] = 4.0E0*I_ESP_Kxy5z_Gx3y_aa-2.0E0*3*I_ESP_Hxy3z_Gx3y_a-2.0E0*4*I_ESP_Hxy3z_Gx3y_a+3*2*I_ESP_Fxyz_Gx3y;
    abcd[iGrid*1890+1715] = 4.0E0*I_ESP_Kx6z_Gx3y_aa-2.0E0*4*I_ESP_Hx4z_Gx3y_a-2.0E0*5*I_ESP_Hx4z_Gx3y_a+4*3*I_ESP_Fx2z_Gx3y;
    abcd[iGrid*1890+1716] = 4.0E0*I_ESP_K5y2z_Gx3y_aa-2.0E0*1*I_ESP_H5y_Gx3y_a;
    abcd[iGrid*1890+1717] = 4.0E0*I_ESP_K4y3z_Gx3y_aa-2.0E0*1*I_ESP_H4yz_Gx3y_a-2.0E0*2*I_ESP_H4yz_Gx3y_a;
    abcd[iGrid*1890+1718] = 4.0E0*I_ESP_K3y4z_Gx3y_aa-2.0E0*2*I_ESP_H3y2z_Gx3y_a-2.0E0*3*I_ESP_H3y2z_Gx3y_a+2*1*I_ESP_F3y_Gx3y;
    abcd[iGrid*1890+1719] = 4.0E0*I_ESP_K2y5z_Gx3y_aa-2.0E0*3*I_ESP_H2y3z_Gx3y_a-2.0E0*4*I_ESP_H2y3z_Gx3y_a+3*2*I_ESP_F2yz_Gx3y;
    abcd[iGrid*1890+1720] = 4.0E0*I_ESP_Ky6z_Gx3y_aa-2.0E0*4*I_ESP_Hy4z_Gx3y_a-2.0E0*5*I_ESP_Hy4z_Gx3y_a+4*3*I_ESP_Fy2z_Gx3y;
    abcd[iGrid*1890+1721] = 4.0E0*I_ESP_K7z_Gx3y_aa-2.0E0*5*I_ESP_H5z_Gx3y_a-2.0E0*6*I_ESP_H5z_Gx3y_a+5*4*I_ESP_F3z_Gx3y;
    abcd[iGrid*1890+1722] = 4.0E0*I_ESP_K5x2z_Gx2yz_aa-2.0E0*1*I_ESP_H5x_Gx2yz_a;
    abcd[iGrid*1890+1723] = 4.0E0*I_ESP_K4xy2z_Gx2yz_aa-2.0E0*1*I_ESP_H4xy_Gx2yz_a;
    abcd[iGrid*1890+1724] = 4.0E0*I_ESP_K4x3z_Gx2yz_aa-2.0E0*1*I_ESP_H4xz_Gx2yz_a-2.0E0*2*I_ESP_H4xz_Gx2yz_a;
    abcd[iGrid*1890+1725] = 4.0E0*I_ESP_K3x2y2z_Gx2yz_aa-2.0E0*1*I_ESP_H3x2y_Gx2yz_a;
    abcd[iGrid*1890+1726] = 4.0E0*I_ESP_K3xy3z_Gx2yz_aa-2.0E0*1*I_ESP_H3xyz_Gx2yz_a-2.0E0*2*I_ESP_H3xyz_Gx2yz_a;
    abcd[iGrid*1890+1727] = 4.0E0*I_ESP_K3x4z_Gx2yz_aa-2.0E0*2*I_ESP_H3x2z_Gx2yz_a-2.0E0*3*I_ESP_H3x2z_Gx2yz_a+2*1*I_ESP_F3x_Gx2yz;
    abcd[iGrid*1890+1728] = 4.0E0*I_ESP_K2x3y2z_Gx2yz_aa-2.0E0*1*I_ESP_H2x3y_Gx2yz_a;
    abcd[iGrid*1890+1729] = 4.0E0*I_ESP_K2x2y3z_Gx2yz_aa-2.0E0*1*I_ESP_H2x2yz_Gx2yz_a-2.0E0*2*I_ESP_H2x2yz_Gx2yz_a;
    abcd[iGrid*1890+1730] = 4.0E0*I_ESP_K2xy4z_Gx2yz_aa-2.0E0*2*I_ESP_H2xy2z_Gx2yz_a-2.0E0*3*I_ESP_H2xy2z_Gx2yz_a+2*1*I_ESP_F2xy_Gx2yz;
    abcd[iGrid*1890+1731] = 4.0E0*I_ESP_K2x5z_Gx2yz_aa-2.0E0*3*I_ESP_H2x3z_Gx2yz_a-2.0E0*4*I_ESP_H2x3z_Gx2yz_a+3*2*I_ESP_F2xz_Gx2yz;
    abcd[iGrid*1890+1732] = 4.0E0*I_ESP_Kx4y2z_Gx2yz_aa-2.0E0*1*I_ESP_Hx4y_Gx2yz_a;
    abcd[iGrid*1890+1733] = 4.0E0*I_ESP_Kx3y3z_Gx2yz_aa-2.0E0*1*I_ESP_Hx3yz_Gx2yz_a-2.0E0*2*I_ESP_Hx3yz_Gx2yz_a;
    abcd[iGrid*1890+1734] = 4.0E0*I_ESP_Kx2y4z_Gx2yz_aa-2.0E0*2*I_ESP_Hx2y2z_Gx2yz_a-2.0E0*3*I_ESP_Hx2y2z_Gx2yz_a+2*1*I_ESP_Fx2y_Gx2yz;
    abcd[iGrid*1890+1735] = 4.0E0*I_ESP_Kxy5z_Gx2yz_aa-2.0E0*3*I_ESP_Hxy3z_Gx2yz_a-2.0E0*4*I_ESP_Hxy3z_Gx2yz_a+3*2*I_ESP_Fxyz_Gx2yz;
    abcd[iGrid*1890+1736] = 4.0E0*I_ESP_Kx6z_Gx2yz_aa-2.0E0*4*I_ESP_Hx4z_Gx2yz_a-2.0E0*5*I_ESP_Hx4z_Gx2yz_a+4*3*I_ESP_Fx2z_Gx2yz;
    abcd[iGrid*1890+1737] = 4.0E0*I_ESP_K5y2z_Gx2yz_aa-2.0E0*1*I_ESP_H5y_Gx2yz_a;
    abcd[iGrid*1890+1738] = 4.0E0*I_ESP_K4y3z_Gx2yz_aa-2.0E0*1*I_ESP_H4yz_Gx2yz_a-2.0E0*2*I_ESP_H4yz_Gx2yz_a;
    abcd[iGrid*1890+1739] = 4.0E0*I_ESP_K3y4z_Gx2yz_aa-2.0E0*2*I_ESP_H3y2z_Gx2yz_a-2.0E0*3*I_ESP_H3y2z_Gx2yz_a+2*1*I_ESP_F3y_Gx2yz;
    abcd[iGrid*1890+1740] = 4.0E0*I_ESP_K2y5z_Gx2yz_aa-2.0E0*3*I_ESP_H2y3z_Gx2yz_a-2.0E0*4*I_ESP_H2y3z_Gx2yz_a+3*2*I_ESP_F2yz_Gx2yz;
    abcd[iGrid*1890+1741] = 4.0E0*I_ESP_Ky6z_Gx2yz_aa-2.0E0*4*I_ESP_Hy4z_Gx2yz_a-2.0E0*5*I_ESP_Hy4z_Gx2yz_a+4*3*I_ESP_Fy2z_Gx2yz;
    abcd[iGrid*1890+1742] = 4.0E0*I_ESP_K7z_Gx2yz_aa-2.0E0*5*I_ESP_H5z_Gx2yz_a-2.0E0*6*I_ESP_H5z_Gx2yz_a+5*4*I_ESP_F3z_Gx2yz;
    abcd[iGrid*1890+1743] = 4.0E0*I_ESP_K5x2z_Gxy2z_aa-2.0E0*1*I_ESP_H5x_Gxy2z_a;
    abcd[iGrid*1890+1744] = 4.0E0*I_ESP_K4xy2z_Gxy2z_aa-2.0E0*1*I_ESP_H4xy_Gxy2z_a;
    abcd[iGrid*1890+1745] = 4.0E0*I_ESP_K4x3z_Gxy2z_aa-2.0E0*1*I_ESP_H4xz_Gxy2z_a-2.0E0*2*I_ESP_H4xz_Gxy2z_a;
    abcd[iGrid*1890+1746] = 4.0E0*I_ESP_K3x2y2z_Gxy2z_aa-2.0E0*1*I_ESP_H3x2y_Gxy2z_a;
    abcd[iGrid*1890+1747] = 4.0E0*I_ESP_K3xy3z_Gxy2z_aa-2.0E0*1*I_ESP_H3xyz_Gxy2z_a-2.0E0*2*I_ESP_H3xyz_Gxy2z_a;
    abcd[iGrid*1890+1748] = 4.0E0*I_ESP_K3x4z_Gxy2z_aa-2.0E0*2*I_ESP_H3x2z_Gxy2z_a-2.0E0*3*I_ESP_H3x2z_Gxy2z_a+2*1*I_ESP_F3x_Gxy2z;
    abcd[iGrid*1890+1749] = 4.0E0*I_ESP_K2x3y2z_Gxy2z_aa-2.0E0*1*I_ESP_H2x3y_Gxy2z_a;
    abcd[iGrid*1890+1750] = 4.0E0*I_ESP_K2x2y3z_Gxy2z_aa-2.0E0*1*I_ESP_H2x2yz_Gxy2z_a-2.0E0*2*I_ESP_H2x2yz_Gxy2z_a;
    abcd[iGrid*1890+1751] = 4.0E0*I_ESP_K2xy4z_Gxy2z_aa-2.0E0*2*I_ESP_H2xy2z_Gxy2z_a-2.0E0*3*I_ESP_H2xy2z_Gxy2z_a+2*1*I_ESP_F2xy_Gxy2z;
    abcd[iGrid*1890+1752] = 4.0E0*I_ESP_K2x5z_Gxy2z_aa-2.0E0*3*I_ESP_H2x3z_Gxy2z_a-2.0E0*4*I_ESP_H2x3z_Gxy2z_a+3*2*I_ESP_F2xz_Gxy2z;
    abcd[iGrid*1890+1753] = 4.0E0*I_ESP_Kx4y2z_Gxy2z_aa-2.0E0*1*I_ESP_Hx4y_Gxy2z_a;
    abcd[iGrid*1890+1754] = 4.0E0*I_ESP_Kx3y3z_Gxy2z_aa-2.0E0*1*I_ESP_Hx3yz_Gxy2z_a-2.0E0*2*I_ESP_Hx3yz_Gxy2z_a;
    abcd[iGrid*1890+1755] = 4.0E0*I_ESP_Kx2y4z_Gxy2z_aa-2.0E0*2*I_ESP_Hx2y2z_Gxy2z_a-2.0E0*3*I_ESP_Hx2y2z_Gxy2z_a+2*1*I_ESP_Fx2y_Gxy2z;
    abcd[iGrid*1890+1756] = 4.0E0*I_ESP_Kxy5z_Gxy2z_aa-2.0E0*3*I_ESP_Hxy3z_Gxy2z_a-2.0E0*4*I_ESP_Hxy3z_Gxy2z_a+3*2*I_ESP_Fxyz_Gxy2z;
    abcd[iGrid*1890+1757] = 4.0E0*I_ESP_Kx6z_Gxy2z_aa-2.0E0*4*I_ESP_Hx4z_Gxy2z_a-2.0E0*5*I_ESP_Hx4z_Gxy2z_a+4*3*I_ESP_Fx2z_Gxy2z;
    abcd[iGrid*1890+1758] = 4.0E0*I_ESP_K5y2z_Gxy2z_aa-2.0E0*1*I_ESP_H5y_Gxy2z_a;
    abcd[iGrid*1890+1759] = 4.0E0*I_ESP_K4y3z_Gxy2z_aa-2.0E0*1*I_ESP_H4yz_Gxy2z_a-2.0E0*2*I_ESP_H4yz_Gxy2z_a;
    abcd[iGrid*1890+1760] = 4.0E0*I_ESP_K3y4z_Gxy2z_aa-2.0E0*2*I_ESP_H3y2z_Gxy2z_a-2.0E0*3*I_ESP_H3y2z_Gxy2z_a+2*1*I_ESP_F3y_Gxy2z;
    abcd[iGrid*1890+1761] = 4.0E0*I_ESP_K2y5z_Gxy2z_aa-2.0E0*3*I_ESP_H2y3z_Gxy2z_a-2.0E0*4*I_ESP_H2y3z_Gxy2z_a+3*2*I_ESP_F2yz_Gxy2z;
    abcd[iGrid*1890+1762] = 4.0E0*I_ESP_Ky6z_Gxy2z_aa-2.0E0*4*I_ESP_Hy4z_Gxy2z_a-2.0E0*5*I_ESP_Hy4z_Gxy2z_a+4*3*I_ESP_Fy2z_Gxy2z;
    abcd[iGrid*1890+1763] = 4.0E0*I_ESP_K7z_Gxy2z_aa-2.0E0*5*I_ESP_H5z_Gxy2z_a-2.0E0*6*I_ESP_H5z_Gxy2z_a+5*4*I_ESP_F3z_Gxy2z;
    abcd[iGrid*1890+1764] = 4.0E0*I_ESP_K5x2z_Gx3z_aa-2.0E0*1*I_ESP_H5x_Gx3z_a;
    abcd[iGrid*1890+1765] = 4.0E0*I_ESP_K4xy2z_Gx3z_aa-2.0E0*1*I_ESP_H4xy_Gx3z_a;
    abcd[iGrid*1890+1766] = 4.0E0*I_ESP_K4x3z_Gx3z_aa-2.0E0*1*I_ESP_H4xz_Gx3z_a-2.0E0*2*I_ESP_H4xz_Gx3z_a;
    abcd[iGrid*1890+1767] = 4.0E0*I_ESP_K3x2y2z_Gx3z_aa-2.0E0*1*I_ESP_H3x2y_Gx3z_a;
    abcd[iGrid*1890+1768] = 4.0E0*I_ESP_K3xy3z_Gx3z_aa-2.0E0*1*I_ESP_H3xyz_Gx3z_a-2.0E0*2*I_ESP_H3xyz_Gx3z_a;
    abcd[iGrid*1890+1769] = 4.0E0*I_ESP_K3x4z_Gx3z_aa-2.0E0*2*I_ESP_H3x2z_Gx3z_a-2.0E0*3*I_ESP_H3x2z_Gx3z_a+2*1*I_ESP_F3x_Gx3z;
    abcd[iGrid*1890+1770] = 4.0E0*I_ESP_K2x3y2z_Gx3z_aa-2.0E0*1*I_ESP_H2x3y_Gx3z_a;
    abcd[iGrid*1890+1771] = 4.0E0*I_ESP_K2x2y3z_Gx3z_aa-2.0E0*1*I_ESP_H2x2yz_Gx3z_a-2.0E0*2*I_ESP_H2x2yz_Gx3z_a;
    abcd[iGrid*1890+1772] = 4.0E0*I_ESP_K2xy4z_Gx3z_aa-2.0E0*2*I_ESP_H2xy2z_Gx3z_a-2.0E0*3*I_ESP_H2xy2z_Gx3z_a+2*1*I_ESP_F2xy_Gx3z;
    abcd[iGrid*1890+1773] = 4.0E0*I_ESP_K2x5z_Gx3z_aa-2.0E0*3*I_ESP_H2x3z_Gx3z_a-2.0E0*4*I_ESP_H2x3z_Gx3z_a+3*2*I_ESP_F2xz_Gx3z;
    abcd[iGrid*1890+1774] = 4.0E0*I_ESP_Kx4y2z_Gx3z_aa-2.0E0*1*I_ESP_Hx4y_Gx3z_a;
    abcd[iGrid*1890+1775] = 4.0E0*I_ESP_Kx3y3z_Gx3z_aa-2.0E0*1*I_ESP_Hx3yz_Gx3z_a-2.0E0*2*I_ESP_Hx3yz_Gx3z_a;
    abcd[iGrid*1890+1776] = 4.0E0*I_ESP_Kx2y4z_Gx3z_aa-2.0E0*2*I_ESP_Hx2y2z_Gx3z_a-2.0E0*3*I_ESP_Hx2y2z_Gx3z_a+2*1*I_ESP_Fx2y_Gx3z;
    abcd[iGrid*1890+1777] = 4.0E0*I_ESP_Kxy5z_Gx3z_aa-2.0E0*3*I_ESP_Hxy3z_Gx3z_a-2.0E0*4*I_ESP_Hxy3z_Gx3z_a+3*2*I_ESP_Fxyz_Gx3z;
    abcd[iGrid*1890+1778] = 4.0E0*I_ESP_Kx6z_Gx3z_aa-2.0E0*4*I_ESP_Hx4z_Gx3z_a-2.0E0*5*I_ESP_Hx4z_Gx3z_a+4*3*I_ESP_Fx2z_Gx3z;
    abcd[iGrid*1890+1779] = 4.0E0*I_ESP_K5y2z_Gx3z_aa-2.0E0*1*I_ESP_H5y_Gx3z_a;
    abcd[iGrid*1890+1780] = 4.0E0*I_ESP_K4y3z_Gx3z_aa-2.0E0*1*I_ESP_H4yz_Gx3z_a-2.0E0*2*I_ESP_H4yz_Gx3z_a;
    abcd[iGrid*1890+1781] = 4.0E0*I_ESP_K3y4z_Gx3z_aa-2.0E0*2*I_ESP_H3y2z_Gx3z_a-2.0E0*3*I_ESP_H3y2z_Gx3z_a+2*1*I_ESP_F3y_Gx3z;
    abcd[iGrid*1890+1782] = 4.0E0*I_ESP_K2y5z_Gx3z_aa-2.0E0*3*I_ESP_H2y3z_Gx3z_a-2.0E0*4*I_ESP_H2y3z_Gx3z_a+3*2*I_ESP_F2yz_Gx3z;
    abcd[iGrid*1890+1783] = 4.0E0*I_ESP_Ky6z_Gx3z_aa-2.0E0*4*I_ESP_Hy4z_Gx3z_a-2.0E0*5*I_ESP_Hy4z_Gx3z_a+4*3*I_ESP_Fy2z_Gx3z;
    abcd[iGrid*1890+1784] = 4.0E0*I_ESP_K7z_Gx3z_aa-2.0E0*5*I_ESP_H5z_Gx3z_a-2.0E0*6*I_ESP_H5z_Gx3z_a+5*4*I_ESP_F3z_Gx3z;
    abcd[iGrid*1890+1785] = 4.0E0*I_ESP_K5x2z_G4y_aa-2.0E0*1*I_ESP_H5x_G4y_a;
    abcd[iGrid*1890+1786] = 4.0E0*I_ESP_K4xy2z_G4y_aa-2.0E0*1*I_ESP_H4xy_G4y_a;
    abcd[iGrid*1890+1787] = 4.0E0*I_ESP_K4x3z_G4y_aa-2.0E0*1*I_ESP_H4xz_G4y_a-2.0E0*2*I_ESP_H4xz_G4y_a;
    abcd[iGrid*1890+1788] = 4.0E0*I_ESP_K3x2y2z_G4y_aa-2.0E0*1*I_ESP_H3x2y_G4y_a;
    abcd[iGrid*1890+1789] = 4.0E0*I_ESP_K3xy3z_G4y_aa-2.0E0*1*I_ESP_H3xyz_G4y_a-2.0E0*2*I_ESP_H3xyz_G4y_a;
    abcd[iGrid*1890+1790] = 4.0E0*I_ESP_K3x4z_G4y_aa-2.0E0*2*I_ESP_H3x2z_G4y_a-2.0E0*3*I_ESP_H3x2z_G4y_a+2*1*I_ESP_F3x_G4y;
    abcd[iGrid*1890+1791] = 4.0E0*I_ESP_K2x3y2z_G4y_aa-2.0E0*1*I_ESP_H2x3y_G4y_a;
    abcd[iGrid*1890+1792] = 4.0E0*I_ESP_K2x2y3z_G4y_aa-2.0E0*1*I_ESP_H2x2yz_G4y_a-2.0E0*2*I_ESP_H2x2yz_G4y_a;
    abcd[iGrid*1890+1793] = 4.0E0*I_ESP_K2xy4z_G4y_aa-2.0E0*2*I_ESP_H2xy2z_G4y_a-2.0E0*3*I_ESP_H2xy2z_G4y_a+2*1*I_ESP_F2xy_G4y;
    abcd[iGrid*1890+1794] = 4.0E0*I_ESP_K2x5z_G4y_aa-2.0E0*3*I_ESP_H2x3z_G4y_a-2.0E0*4*I_ESP_H2x3z_G4y_a+3*2*I_ESP_F2xz_G4y;
    abcd[iGrid*1890+1795] = 4.0E0*I_ESP_Kx4y2z_G4y_aa-2.0E0*1*I_ESP_Hx4y_G4y_a;
    abcd[iGrid*1890+1796] = 4.0E0*I_ESP_Kx3y3z_G4y_aa-2.0E0*1*I_ESP_Hx3yz_G4y_a-2.0E0*2*I_ESP_Hx3yz_G4y_a;
    abcd[iGrid*1890+1797] = 4.0E0*I_ESP_Kx2y4z_G4y_aa-2.0E0*2*I_ESP_Hx2y2z_G4y_a-2.0E0*3*I_ESP_Hx2y2z_G4y_a+2*1*I_ESP_Fx2y_G4y;
    abcd[iGrid*1890+1798] = 4.0E0*I_ESP_Kxy5z_G4y_aa-2.0E0*3*I_ESP_Hxy3z_G4y_a-2.0E0*4*I_ESP_Hxy3z_G4y_a+3*2*I_ESP_Fxyz_G4y;
    abcd[iGrid*1890+1799] = 4.0E0*I_ESP_Kx6z_G4y_aa-2.0E0*4*I_ESP_Hx4z_G4y_a-2.0E0*5*I_ESP_Hx4z_G4y_a+4*3*I_ESP_Fx2z_G4y;
    abcd[iGrid*1890+1800] = 4.0E0*I_ESP_K5y2z_G4y_aa-2.0E0*1*I_ESP_H5y_G4y_a;
    abcd[iGrid*1890+1801] = 4.0E0*I_ESP_K4y3z_G4y_aa-2.0E0*1*I_ESP_H4yz_G4y_a-2.0E0*2*I_ESP_H4yz_G4y_a;
    abcd[iGrid*1890+1802] = 4.0E0*I_ESP_K3y4z_G4y_aa-2.0E0*2*I_ESP_H3y2z_G4y_a-2.0E0*3*I_ESP_H3y2z_G4y_a+2*1*I_ESP_F3y_G4y;
    abcd[iGrid*1890+1803] = 4.0E0*I_ESP_K2y5z_G4y_aa-2.0E0*3*I_ESP_H2y3z_G4y_a-2.0E0*4*I_ESP_H2y3z_G4y_a+3*2*I_ESP_F2yz_G4y;
    abcd[iGrid*1890+1804] = 4.0E0*I_ESP_Ky6z_G4y_aa-2.0E0*4*I_ESP_Hy4z_G4y_a-2.0E0*5*I_ESP_Hy4z_G4y_a+4*3*I_ESP_Fy2z_G4y;
    abcd[iGrid*1890+1805] = 4.0E0*I_ESP_K7z_G4y_aa-2.0E0*5*I_ESP_H5z_G4y_a-2.0E0*6*I_ESP_H5z_G4y_a+5*4*I_ESP_F3z_G4y;
    abcd[iGrid*1890+1806] = 4.0E0*I_ESP_K5x2z_G3yz_aa-2.0E0*1*I_ESP_H5x_G3yz_a;
    abcd[iGrid*1890+1807] = 4.0E0*I_ESP_K4xy2z_G3yz_aa-2.0E0*1*I_ESP_H4xy_G3yz_a;
    abcd[iGrid*1890+1808] = 4.0E0*I_ESP_K4x3z_G3yz_aa-2.0E0*1*I_ESP_H4xz_G3yz_a-2.0E0*2*I_ESP_H4xz_G3yz_a;
    abcd[iGrid*1890+1809] = 4.0E0*I_ESP_K3x2y2z_G3yz_aa-2.0E0*1*I_ESP_H3x2y_G3yz_a;
    abcd[iGrid*1890+1810] = 4.0E0*I_ESP_K3xy3z_G3yz_aa-2.0E0*1*I_ESP_H3xyz_G3yz_a-2.0E0*2*I_ESP_H3xyz_G3yz_a;
    abcd[iGrid*1890+1811] = 4.0E0*I_ESP_K3x4z_G3yz_aa-2.0E0*2*I_ESP_H3x2z_G3yz_a-2.0E0*3*I_ESP_H3x2z_G3yz_a+2*1*I_ESP_F3x_G3yz;
    abcd[iGrid*1890+1812] = 4.0E0*I_ESP_K2x3y2z_G3yz_aa-2.0E0*1*I_ESP_H2x3y_G3yz_a;
    abcd[iGrid*1890+1813] = 4.0E0*I_ESP_K2x2y3z_G3yz_aa-2.0E0*1*I_ESP_H2x2yz_G3yz_a-2.0E0*2*I_ESP_H2x2yz_G3yz_a;
    abcd[iGrid*1890+1814] = 4.0E0*I_ESP_K2xy4z_G3yz_aa-2.0E0*2*I_ESP_H2xy2z_G3yz_a-2.0E0*3*I_ESP_H2xy2z_G3yz_a+2*1*I_ESP_F2xy_G3yz;
    abcd[iGrid*1890+1815] = 4.0E0*I_ESP_K2x5z_G3yz_aa-2.0E0*3*I_ESP_H2x3z_G3yz_a-2.0E0*4*I_ESP_H2x3z_G3yz_a+3*2*I_ESP_F2xz_G3yz;
    abcd[iGrid*1890+1816] = 4.0E0*I_ESP_Kx4y2z_G3yz_aa-2.0E0*1*I_ESP_Hx4y_G3yz_a;
    abcd[iGrid*1890+1817] = 4.0E0*I_ESP_Kx3y3z_G3yz_aa-2.0E0*1*I_ESP_Hx3yz_G3yz_a-2.0E0*2*I_ESP_Hx3yz_G3yz_a;
    abcd[iGrid*1890+1818] = 4.0E0*I_ESP_Kx2y4z_G3yz_aa-2.0E0*2*I_ESP_Hx2y2z_G3yz_a-2.0E0*3*I_ESP_Hx2y2z_G3yz_a+2*1*I_ESP_Fx2y_G3yz;
    abcd[iGrid*1890+1819] = 4.0E0*I_ESP_Kxy5z_G3yz_aa-2.0E0*3*I_ESP_Hxy3z_G3yz_a-2.0E0*4*I_ESP_Hxy3z_G3yz_a+3*2*I_ESP_Fxyz_G3yz;
    abcd[iGrid*1890+1820] = 4.0E0*I_ESP_Kx6z_G3yz_aa-2.0E0*4*I_ESP_Hx4z_G3yz_a-2.0E0*5*I_ESP_Hx4z_G3yz_a+4*3*I_ESP_Fx2z_G3yz;
    abcd[iGrid*1890+1821] = 4.0E0*I_ESP_K5y2z_G3yz_aa-2.0E0*1*I_ESP_H5y_G3yz_a;
    abcd[iGrid*1890+1822] = 4.0E0*I_ESP_K4y3z_G3yz_aa-2.0E0*1*I_ESP_H4yz_G3yz_a-2.0E0*2*I_ESP_H4yz_G3yz_a;
    abcd[iGrid*1890+1823] = 4.0E0*I_ESP_K3y4z_G3yz_aa-2.0E0*2*I_ESP_H3y2z_G3yz_a-2.0E0*3*I_ESP_H3y2z_G3yz_a+2*1*I_ESP_F3y_G3yz;
    abcd[iGrid*1890+1824] = 4.0E0*I_ESP_K2y5z_G3yz_aa-2.0E0*3*I_ESP_H2y3z_G3yz_a-2.0E0*4*I_ESP_H2y3z_G3yz_a+3*2*I_ESP_F2yz_G3yz;
    abcd[iGrid*1890+1825] = 4.0E0*I_ESP_Ky6z_G3yz_aa-2.0E0*4*I_ESP_Hy4z_G3yz_a-2.0E0*5*I_ESP_Hy4z_G3yz_a+4*3*I_ESP_Fy2z_G3yz;
    abcd[iGrid*1890+1826] = 4.0E0*I_ESP_K7z_G3yz_aa-2.0E0*5*I_ESP_H5z_G3yz_a-2.0E0*6*I_ESP_H5z_G3yz_a+5*4*I_ESP_F3z_G3yz;
    abcd[iGrid*1890+1827] = 4.0E0*I_ESP_K5x2z_G2y2z_aa-2.0E0*1*I_ESP_H5x_G2y2z_a;
    abcd[iGrid*1890+1828] = 4.0E0*I_ESP_K4xy2z_G2y2z_aa-2.0E0*1*I_ESP_H4xy_G2y2z_a;
    abcd[iGrid*1890+1829] = 4.0E0*I_ESP_K4x3z_G2y2z_aa-2.0E0*1*I_ESP_H4xz_G2y2z_a-2.0E0*2*I_ESP_H4xz_G2y2z_a;
    abcd[iGrid*1890+1830] = 4.0E0*I_ESP_K3x2y2z_G2y2z_aa-2.0E0*1*I_ESP_H3x2y_G2y2z_a;
    abcd[iGrid*1890+1831] = 4.0E0*I_ESP_K3xy3z_G2y2z_aa-2.0E0*1*I_ESP_H3xyz_G2y2z_a-2.0E0*2*I_ESP_H3xyz_G2y2z_a;
    abcd[iGrid*1890+1832] = 4.0E0*I_ESP_K3x4z_G2y2z_aa-2.0E0*2*I_ESP_H3x2z_G2y2z_a-2.0E0*3*I_ESP_H3x2z_G2y2z_a+2*1*I_ESP_F3x_G2y2z;
    abcd[iGrid*1890+1833] = 4.0E0*I_ESP_K2x3y2z_G2y2z_aa-2.0E0*1*I_ESP_H2x3y_G2y2z_a;
    abcd[iGrid*1890+1834] = 4.0E0*I_ESP_K2x2y3z_G2y2z_aa-2.0E0*1*I_ESP_H2x2yz_G2y2z_a-2.0E0*2*I_ESP_H2x2yz_G2y2z_a;
    abcd[iGrid*1890+1835] = 4.0E0*I_ESP_K2xy4z_G2y2z_aa-2.0E0*2*I_ESP_H2xy2z_G2y2z_a-2.0E0*3*I_ESP_H2xy2z_G2y2z_a+2*1*I_ESP_F2xy_G2y2z;
    abcd[iGrid*1890+1836] = 4.0E0*I_ESP_K2x5z_G2y2z_aa-2.0E0*3*I_ESP_H2x3z_G2y2z_a-2.0E0*4*I_ESP_H2x3z_G2y2z_a+3*2*I_ESP_F2xz_G2y2z;
    abcd[iGrid*1890+1837] = 4.0E0*I_ESP_Kx4y2z_G2y2z_aa-2.0E0*1*I_ESP_Hx4y_G2y2z_a;
    abcd[iGrid*1890+1838] = 4.0E0*I_ESP_Kx3y3z_G2y2z_aa-2.0E0*1*I_ESP_Hx3yz_G2y2z_a-2.0E0*2*I_ESP_Hx3yz_G2y2z_a;
    abcd[iGrid*1890+1839] = 4.0E0*I_ESP_Kx2y4z_G2y2z_aa-2.0E0*2*I_ESP_Hx2y2z_G2y2z_a-2.0E0*3*I_ESP_Hx2y2z_G2y2z_a+2*1*I_ESP_Fx2y_G2y2z;
    abcd[iGrid*1890+1840] = 4.0E0*I_ESP_Kxy5z_G2y2z_aa-2.0E0*3*I_ESP_Hxy3z_G2y2z_a-2.0E0*4*I_ESP_Hxy3z_G2y2z_a+3*2*I_ESP_Fxyz_G2y2z;
    abcd[iGrid*1890+1841] = 4.0E0*I_ESP_Kx6z_G2y2z_aa-2.0E0*4*I_ESP_Hx4z_G2y2z_a-2.0E0*5*I_ESP_Hx4z_G2y2z_a+4*3*I_ESP_Fx2z_G2y2z;
    abcd[iGrid*1890+1842] = 4.0E0*I_ESP_K5y2z_G2y2z_aa-2.0E0*1*I_ESP_H5y_G2y2z_a;
    abcd[iGrid*1890+1843] = 4.0E0*I_ESP_K4y3z_G2y2z_aa-2.0E0*1*I_ESP_H4yz_G2y2z_a-2.0E0*2*I_ESP_H4yz_G2y2z_a;
    abcd[iGrid*1890+1844] = 4.0E0*I_ESP_K3y4z_G2y2z_aa-2.0E0*2*I_ESP_H3y2z_G2y2z_a-2.0E0*3*I_ESP_H3y2z_G2y2z_a+2*1*I_ESP_F3y_G2y2z;
    abcd[iGrid*1890+1845] = 4.0E0*I_ESP_K2y5z_G2y2z_aa-2.0E0*3*I_ESP_H2y3z_G2y2z_a-2.0E0*4*I_ESP_H2y3z_G2y2z_a+3*2*I_ESP_F2yz_G2y2z;
    abcd[iGrid*1890+1846] = 4.0E0*I_ESP_Ky6z_G2y2z_aa-2.0E0*4*I_ESP_Hy4z_G2y2z_a-2.0E0*5*I_ESP_Hy4z_G2y2z_a+4*3*I_ESP_Fy2z_G2y2z;
    abcd[iGrid*1890+1847] = 4.0E0*I_ESP_K7z_G2y2z_aa-2.0E0*5*I_ESP_H5z_G2y2z_a-2.0E0*6*I_ESP_H5z_G2y2z_a+5*4*I_ESP_F3z_G2y2z;
    abcd[iGrid*1890+1848] = 4.0E0*I_ESP_K5x2z_Gy3z_aa-2.0E0*1*I_ESP_H5x_Gy3z_a;
    abcd[iGrid*1890+1849] = 4.0E0*I_ESP_K4xy2z_Gy3z_aa-2.0E0*1*I_ESP_H4xy_Gy3z_a;
    abcd[iGrid*1890+1850] = 4.0E0*I_ESP_K4x3z_Gy3z_aa-2.0E0*1*I_ESP_H4xz_Gy3z_a-2.0E0*2*I_ESP_H4xz_Gy3z_a;
    abcd[iGrid*1890+1851] = 4.0E0*I_ESP_K3x2y2z_Gy3z_aa-2.0E0*1*I_ESP_H3x2y_Gy3z_a;
    abcd[iGrid*1890+1852] = 4.0E0*I_ESP_K3xy3z_Gy3z_aa-2.0E0*1*I_ESP_H3xyz_Gy3z_a-2.0E0*2*I_ESP_H3xyz_Gy3z_a;
    abcd[iGrid*1890+1853] = 4.0E0*I_ESP_K3x4z_Gy3z_aa-2.0E0*2*I_ESP_H3x2z_Gy3z_a-2.0E0*3*I_ESP_H3x2z_Gy3z_a+2*1*I_ESP_F3x_Gy3z;
    abcd[iGrid*1890+1854] = 4.0E0*I_ESP_K2x3y2z_Gy3z_aa-2.0E0*1*I_ESP_H2x3y_Gy3z_a;
    abcd[iGrid*1890+1855] = 4.0E0*I_ESP_K2x2y3z_Gy3z_aa-2.0E0*1*I_ESP_H2x2yz_Gy3z_a-2.0E0*2*I_ESP_H2x2yz_Gy3z_a;
    abcd[iGrid*1890+1856] = 4.0E0*I_ESP_K2xy4z_Gy3z_aa-2.0E0*2*I_ESP_H2xy2z_Gy3z_a-2.0E0*3*I_ESP_H2xy2z_Gy3z_a+2*1*I_ESP_F2xy_Gy3z;
    abcd[iGrid*1890+1857] = 4.0E0*I_ESP_K2x5z_Gy3z_aa-2.0E0*3*I_ESP_H2x3z_Gy3z_a-2.0E0*4*I_ESP_H2x3z_Gy3z_a+3*2*I_ESP_F2xz_Gy3z;
    abcd[iGrid*1890+1858] = 4.0E0*I_ESP_Kx4y2z_Gy3z_aa-2.0E0*1*I_ESP_Hx4y_Gy3z_a;
    abcd[iGrid*1890+1859] = 4.0E0*I_ESP_Kx3y3z_Gy3z_aa-2.0E0*1*I_ESP_Hx3yz_Gy3z_a-2.0E0*2*I_ESP_Hx3yz_Gy3z_a;
    abcd[iGrid*1890+1860] = 4.0E0*I_ESP_Kx2y4z_Gy3z_aa-2.0E0*2*I_ESP_Hx2y2z_Gy3z_a-2.0E0*3*I_ESP_Hx2y2z_Gy3z_a+2*1*I_ESP_Fx2y_Gy3z;
    abcd[iGrid*1890+1861] = 4.0E0*I_ESP_Kxy5z_Gy3z_aa-2.0E0*3*I_ESP_Hxy3z_Gy3z_a-2.0E0*4*I_ESP_Hxy3z_Gy3z_a+3*2*I_ESP_Fxyz_Gy3z;
    abcd[iGrid*1890+1862] = 4.0E0*I_ESP_Kx6z_Gy3z_aa-2.0E0*4*I_ESP_Hx4z_Gy3z_a-2.0E0*5*I_ESP_Hx4z_Gy3z_a+4*3*I_ESP_Fx2z_Gy3z;
    abcd[iGrid*1890+1863] = 4.0E0*I_ESP_K5y2z_Gy3z_aa-2.0E0*1*I_ESP_H5y_Gy3z_a;
    abcd[iGrid*1890+1864] = 4.0E0*I_ESP_K4y3z_Gy3z_aa-2.0E0*1*I_ESP_H4yz_Gy3z_a-2.0E0*2*I_ESP_H4yz_Gy3z_a;
    abcd[iGrid*1890+1865] = 4.0E0*I_ESP_K3y4z_Gy3z_aa-2.0E0*2*I_ESP_H3y2z_Gy3z_a-2.0E0*3*I_ESP_H3y2z_Gy3z_a+2*1*I_ESP_F3y_Gy3z;
    abcd[iGrid*1890+1866] = 4.0E0*I_ESP_K2y5z_Gy3z_aa-2.0E0*3*I_ESP_H2y3z_Gy3z_a-2.0E0*4*I_ESP_H2y3z_Gy3z_a+3*2*I_ESP_F2yz_Gy3z;
    abcd[iGrid*1890+1867] = 4.0E0*I_ESP_Ky6z_Gy3z_aa-2.0E0*4*I_ESP_Hy4z_Gy3z_a-2.0E0*5*I_ESP_Hy4z_Gy3z_a+4*3*I_ESP_Fy2z_Gy3z;
    abcd[iGrid*1890+1868] = 4.0E0*I_ESP_K7z_Gy3z_aa-2.0E0*5*I_ESP_H5z_Gy3z_a-2.0E0*6*I_ESP_H5z_Gy3z_a+5*4*I_ESP_F3z_Gy3z;
    abcd[iGrid*1890+1869] = 4.0E0*I_ESP_K5x2z_G4z_aa-2.0E0*1*I_ESP_H5x_G4z_a;
    abcd[iGrid*1890+1870] = 4.0E0*I_ESP_K4xy2z_G4z_aa-2.0E0*1*I_ESP_H4xy_G4z_a;
    abcd[iGrid*1890+1871] = 4.0E0*I_ESP_K4x3z_G4z_aa-2.0E0*1*I_ESP_H4xz_G4z_a-2.0E0*2*I_ESP_H4xz_G4z_a;
    abcd[iGrid*1890+1872] = 4.0E0*I_ESP_K3x2y2z_G4z_aa-2.0E0*1*I_ESP_H3x2y_G4z_a;
    abcd[iGrid*1890+1873] = 4.0E0*I_ESP_K3xy3z_G4z_aa-2.0E0*1*I_ESP_H3xyz_G4z_a-2.0E0*2*I_ESP_H3xyz_G4z_a;
    abcd[iGrid*1890+1874] = 4.0E0*I_ESP_K3x4z_G4z_aa-2.0E0*2*I_ESP_H3x2z_G4z_a-2.0E0*3*I_ESP_H3x2z_G4z_a+2*1*I_ESP_F3x_G4z;
    abcd[iGrid*1890+1875] = 4.0E0*I_ESP_K2x3y2z_G4z_aa-2.0E0*1*I_ESP_H2x3y_G4z_a;
    abcd[iGrid*1890+1876] = 4.0E0*I_ESP_K2x2y3z_G4z_aa-2.0E0*1*I_ESP_H2x2yz_G4z_a-2.0E0*2*I_ESP_H2x2yz_G4z_a;
    abcd[iGrid*1890+1877] = 4.0E0*I_ESP_K2xy4z_G4z_aa-2.0E0*2*I_ESP_H2xy2z_G4z_a-2.0E0*3*I_ESP_H2xy2z_G4z_a+2*1*I_ESP_F2xy_G4z;
    abcd[iGrid*1890+1878] = 4.0E0*I_ESP_K2x5z_G4z_aa-2.0E0*3*I_ESP_H2x3z_G4z_a-2.0E0*4*I_ESP_H2x3z_G4z_a+3*2*I_ESP_F2xz_G4z;
    abcd[iGrid*1890+1879] = 4.0E0*I_ESP_Kx4y2z_G4z_aa-2.0E0*1*I_ESP_Hx4y_G4z_a;
    abcd[iGrid*1890+1880] = 4.0E0*I_ESP_Kx3y3z_G4z_aa-2.0E0*1*I_ESP_Hx3yz_G4z_a-2.0E0*2*I_ESP_Hx3yz_G4z_a;
    abcd[iGrid*1890+1881] = 4.0E0*I_ESP_Kx2y4z_G4z_aa-2.0E0*2*I_ESP_Hx2y2z_G4z_a-2.0E0*3*I_ESP_Hx2y2z_G4z_a+2*1*I_ESP_Fx2y_G4z;
    abcd[iGrid*1890+1882] = 4.0E0*I_ESP_Kxy5z_G4z_aa-2.0E0*3*I_ESP_Hxy3z_G4z_a-2.0E0*4*I_ESP_Hxy3z_G4z_a+3*2*I_ESP_Fxyz_G4z;
    abcd[iGrid*1890+1883] = 4.0E0*I_ESP_Kx6z_G4z_aa-2.0E0*4*I_ESP_Hx4z_G4z_a-2.0E0*5*I_ESP_Hx4z_G4z_a+4*3*I_ESP_Fx2z_G4z;
    abcd[iGrid*1890+1884] = 4.0E0*I_ESP_K5y2z_G4z_aa-2.0E0*1*I_ESP_H5y_G4z_a;
    abcd[iGrid*1890+1885] = 4.0E0*I_ESP_K4y3z_G4z_aa-2.0E0*1*I_ESP_H4yz_G4z_a-2.0E0*2*I_ESP_H4yz_G4z_a;
    abcd[iGrid*1890+1886] = 4.0E0*I_ESP_K3y4z_G4z_aa-2.0E0*2*I_ESP_H3y2z_G4z_a-2.0E0*3*I_ESP_H3y2z_G4z_a+2*1*I_ESP_F3y_G4z;
    abcd[iGrid*1890+1887] = 4.0E0*I_ESP_K2y5z_G4z_aa-2.0E0*3*I_ESP_H2y3z_G4z_a-2.0E0*4*I_ESP_H2y3z_G4z_a+3*2*I_ESP_F2yz_G4z;
    abcd[iGrid*1890+1888] = 4.0E0*I_ESP_Ky6z_G4z_aa-2.0E0*4*I_ESP_Hy4z_G4z_a-2.0E0*5*I_ESP_Hy4z_G4z_a+4*3*I_ESP_Fy2z_G4z;
    abcd[iGrid*1890+1889] = 4.0E0*I_ESP_K7z_G4z_aa-2.0E0*5*I_ESP_H5z_G4z_a-2.0E0*6*I_ESP_H5z_G4z_a+5*4*I_ESP_F3z_G4z;
  }
}
