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

void hgp_os_esp_h_h_d1(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_O11x_S_a = 0.0E0;
    Double I_ESP_O10xy_S_a = 0.0E0;
    Double I_ESP_O10xz_S_a = 0.0E0;
    Double I_ESP_O9x2y_S_a = 0.0E0;
    Double I_ESP_O9xyz_S_a = 0.0E0;
    Double I_ESP_O9x2z_S_a = 0.0E0;
    Double I_ESP_O8x3y_S_a = 0.0E0;
    Double I_ESP_O8x2yz_S_a = 0.0E0;
    Double I_ESP_O8xy2z_S_a = 0.0E0;
    Double I_ESP_O8x3z_S_a = 0.0E0;
    Double I_ESP_O7x4y_S_a = 0.0E0;
    Double I_ESP_O7x3yz_S_a = 0.0E0;
    Double I_ESP_O7x2y2z_S_a = 0.0E0;
    Double I_ESP_O7xy3z_S_a = 0.0E0;
    Double I_ESP_O7x4z_S_a = 0.0E0;
    Double I_ESP_O6x5y_S_a = 0.0E0;
    Double I_ESP_O6x4yz_S_a = 0.0E0;
    Double I_ESP_O6x3y2z_S_a = 0.0E0;
    Double I_ESP_O6x2y3z_S_a = 0.0E0;
    Double I_ESP_O6xy4z_S_a = 0.0E0;
    Double I_ESP_O6x5z_S_a = 0.0E0;
    Double I_ESP_O5x6y_S_a = 0.0E0;
    Double I_ESP_O5x5yz_S_a = 0.0E0;
    Double I_ESP_O5x4y2z_S_a = 0.0E0;
    Double I_ESP_O5x3y3z_S_a = 0.0E0;
    Double I_ESP_O5x2y4z_S_a = 0.0E0;
    Double I_ESP_O5xy5z_S_a = 0.0E0;
    Double I_ESP_O5x6z_S_a = 0.0E0;
    Double I_ESP_O4x7y_S_a = 0.0E0;
    Double I_ESP_O4x6yz_S_a = 0.0E0;
    Double I_ESP_O4x5y2z_S_a = 0.0E0;
    Double I_ESP_O4x4y3z_S_a = 0.0E0;
    Double I_ESP_O4x3y4z_S_a = 0.0E0;
    Double I_ESP_O4x2y5z_S_a = 0.0E0;
    Double I_ESP_O4xy6z_S_a = 0.0E0;
    Double I_ESP_O4x7z_S_a = 0.0E0;
    Double I_ESP_O3x8y_S_a = 0.0E0;
    Double I_ESP_O3x7yz_S_a = 0.0E0;
    Double I_ESP_O3x6y2z_S_a = 0.0E0;
    Double I_ESP_O3x5y3z_S_a = 0.0E0;
    Double I_ESP_O3x4y4z_S_a = 0.0E0;
    Double I_ESP_O3x3y5z_S_a = 0.0E0;
    Double I_ESP_O3x2y6z_S_a = 0.0E0;
    Double I_ESP_O3xy7z_S_a = 0.0E0;
    Double I_ESP_O3x8z_S_a = 0.0E0;
    Double I_ESP_O2x9y_S_a = 0.0E0;
    Double I_ESP_O2x8yz_S_a = 0.0E0;
    Double I_ESP_O2x7y2z_S_a = 0.0E0;
    Double I_ESP_O2x6y3z_S_a = 0.0E0;
    Double I_ESP_O2x5y4z_S_a = 0.0E0;
    Double I_ESP_O2x4y5z_S_a = 0.0E0;
    Double I_ESP_O2x3y6z_S_a = 0.0E0;
    Double I_ESP_O2x2y7z_S_a = 0.0E0;
    Double I_ESP_O2xy8z_S_a = 0.0E0;
    Double I_ESP_O2x9z_S_a = 0.0E0;
    Double I_ESP_Ox10y_S_a = 0.0E0;
    Double I_ESP_Ox9yz_S_a = 0.0E0;
    Double I_ESP_Ox8y2z_S_a = 0.0E0;
    Double I_ESP_Ox7y3z_S_a = 0.0E0;
    Double I_ESP_Ox6y4z_S_a = 0.0E0;
    Double I_ESP_Ox5y5z_S_a = 0.0E0;
    Double I_ESP_Ox4y6z_S_a = 0.0E0;
    Double I_ESP_Ox3y7z_S_a = 0.0E0;
    Double I_ESP_Ox2y8z_S_a = 0.0E0;
    Double I_ESP_Oxy9z_S_a = 0.0E0;
    Double I_ESP_Ox10z_S_a = 0.0E0;
    Double I_ESP_O11y_S_a = 0.0E0;
    Double I_ESP_O10yz_S_a = 0.0E0;
    Double I_ESP_O9y2z_S_a = 0.0E0;
    Double I_ESP_O8y3z_S_a = 0.0E0;
    Double I_ESP_O7y4z_S_a = 0.0E0;
    Double I_ESP_O6y5z_S_a = 0.0E0;
    Double I_ESP_O5y6z_S_a = 0.0E0;
    Double I_ESP_O4y7z_S_a = 0.0E0;
    Double I_ESP_O3y8z_S_a = 0.0E0;
    Double I_ESP_O2y9z_S_a = 0.0E0;
    Double I_ESP_Oy10z_S_a = 0.0E0;
    Double I_ESP_O11z_S_a = 0.0E0;
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
    Double I_ESP_M9x_S = 0.0E0;
    Double I_ESP_M8xy_S = 0.0E0;
    Double I_ESP_M8xz_S = 0.0E0;
    Double I_ESP_M7x2y_S = 0.0E0;
    Double I_ESP_M7xyz_S = 0.0E0;
    Double I_ESP_M7x2z_S = 0.0E0;
    Double I_ESP_M6x3y_S = 0.0E0;
    Double I_ESP_M6x2yz_S = 0.0E0;
    Double I_ESP_M6xy2z_S = 0.0E0;
    Double I_ESP_M6x3z_S = 0.0E0;
    Double I_ESP_M5x4y_S = 0.0E0;
    Double I_ESP_M5x3yz_S = 0.0E0;
    Double I_ESP_M5x2y2z_S = 0.0E0;
    Double I_ESP_M5xy3z_S = 0.0E0;
    Double I_ESP_M5x4z_S = 0.0E0;
    Double I_ESP_M4x5y_S = 0.0E0;
    Double I_ESP_M4x4yz_S = 0.0E0;
    Double I_ESP_M4x3y2z_S = 0.0E0;
    Double I_ESP_M4x2y3z_S = 0.0E0;
    Double I_ESP_M4xy4z_S = 0.0E0;
    Double I_ESP_M4x5z_S = 0.0E0;
    Double I_ESP_M3x6y_S = 0.0E0;
    Double I_ESP_M3x5yz_S = 0.0E0;
    Double I_ESP_M3x4y2z_S = 0.0E0;
    Double I_ESP_M3x3y3z_S = 0.0E0;
    Double I_ESP_M3x2y4z_S = 0.0E0;
    Double I_ESP_M3xy5z_S = 0.0E0;
    Double I_ESP_M3x6z_S = 0.0E0;
    Double I_ESP_M2x7y_S = 0.0E0;
    Double I_ESP_M2x6yz_S = 0.0E0;
    Double I_ESP_M2x5y2z_S = 0.0E0;
    Double I_ESP_M2x4y3z_S = 0.0E0;
    Double I_ESP_M2x3y4z_S = 0.0E0;
    Double I_ESP_M2x2y5z_S = 0.0E0;
    Double I_ESP_M2xy6z_S = 0.0E0;
    Double I_ESP_M2x7z_S = 0.0E0;
    Double I_ESP_Mx8y_S = 0.0E0;
    Double I_ESP_Mx7yz_S = 0.0E0;
    Double I_ESP_Mx6y2z_S = 0.0E0;
    Double I_ESP_Mx5y3z_S = 0.0E0;
    Double I_ESP_Mx4y4z_S = 0.0E0;
    Double I_ESP_Mx3y5z_S = 0.0E0;
    Double I_ESP_Mx2y6z_S = 0.0E0;
    Double I_ESP_Mxy7z_S = 0.0E0;
    Double I_ESP_Mx8z_S = 0.0E0;
    Double I_ESP_M9y_S = 0.0E0;
    Double I_ESP_M8yz_S = 0.0E0;
    Double I_ESP_M7y2z_S = 0.0E0;
    Double I_ESP_M6y3z_S = 0.0E0;
    Double I_ESP_M5y4z_S = 0.0E0;
    Double I_ESP_M4y5z_S = 0.0E0;
    Double I_ESP_M3y6z_S = 0.0E0;
    Double I_ESP_M2y7z_S = 0.0E0;
    Double I_ESP_My8z_S = 0.0E0;
    Double I_ESP_M9z_S = 0.0E0;
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
         * shell quartet name: SQ_ESP_O_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_O_S_a_coefs = alpha;
        I_ESP_O11x_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O11x_S_vrr;
        I_ESP_O10xy_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O10xy_S_vrr;
        I_ESP_O10xz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O10xz_S_vrr;
        I_ESP_O9x2y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O9x2y_S_vrr;
        I_ESP_O9xyz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O9xyz_S_vrr;
        I_ESP_O9x2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O9x2z_S_vrr;
        I_ESP_O8x3y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O8x3y_S_vrr;
        I_ESP_O8x2yz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O8x2yz_S_vrr;
        I_ESP_O8xy2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O8xy2z_S_vrr;
        I_ESP_O8x3z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O8x3z_S_vrr;
        I_ESP_O7x4y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O7x4y_S_vrr;
        I_ESP_O7x3yz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O7x3yz_S_vrr;
        I_ESP_O7x2y2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O7x2y2z_S_vrr;
        I_ESP_O7xy3z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O7xy3z_S_vrr;
        I_ESP_O7x4z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O7x4z_S_vrr;
        I_ESP_O6x5y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O6x5y_S_vrr;
        I_ESP_O6x4yz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O6x4yz_S_vrr;
        I_ESP_O6x3y2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O6x3y2z_S_vrr;
        I_ESP_O6x2y3z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O6x2y3z_S_vrr;
        I_ESP_O6xy4z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O6xy4z_S_vrr;
        I_ESP_O6x5z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O6x5z_S_vrr;
        I_ESP_O5x6y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O5x6y_S_vrr;
        I_ESP_O5x5yz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O5x5yz_S_vrr;
        I_ESP_O5x4y2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O5x4y2z_S_vrr;
        I_ESP_O5x3y3z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O5x3y3z_S_vrr;
        I_ESP_O5x2y4z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O5x2y4z_S_vrr;
        I_ESP_O5xy5z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O5xy5z_S_vrr;
        I_ESP_O5x6z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O5x6z_S_vrr;
        I_ESP_O4x7y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O4x7y_S_vrr;
        I_ESP_O4x6yz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O4x6yz_S_vrr;
        I_ESP_O4x5y2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O4x5y2z_S_vrr;
        I_ESP_O4x4y3z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O4x4y3z_S_vrr;
        I_ESP_O4x3y4z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O4x3y4z_S_vrr;
        I_ESP_O4x2y5z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O4x2y5z_S_vrr;
        I_ESP_O4xy6z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O4xy6z_S_vrr;
        I_ESP_O4x7z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O4x7z_S_vrr;
        I_ESP_O3x8y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3x8y_S_vrr;
        I_ESP_O3x7yz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3x7yz_S_vrr;
        I_ESP_O3x6y2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3x6y2z_S_vrr;
        I_ESP_O3x5y3z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3x5y3z_S_vrr;
        I_ESP_O3x4y4z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3x4y4z_S_vrr;
        I_ESP_O3x3y5z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3x3y5z_S_vrr;
        I_ESP_O3x2y6z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3x2y6z_S_vrr;
        I_ESP_O3xy7z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3xy7z_S_vrr;
        I_ESP_O3x8z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3x8z_S_vrr;
        I_ESP_O2x9y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2x9y_S_vrr;
        I_ESP_O2x8yz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2x8yz_S_vrr;
        I_ESP_O2x7y2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2x7y2z_S_vrr;
        I_ESP_O2x6y3z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2x6y3z_S_vrr;
        I_ESP_O2x5y4z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2x5y4z_S_vrr;
        I_ESP_O2x4y5z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2x4y5z_S_vrr;
        I_ESP_O2x3y6z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2x3y6z_S_vrr;
        I_ESP_O2x2y7z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2x2y7z_S_vrr;
        I_ESP_O2xy8z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2xy8z_S_vrr;
        I_ESP_O2x9z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2x9z_S_vrr;
        I_ESP_Ox10y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox10y_S_vrr;
        I_ESP_Ox9yz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox9yz_S_vrr;
        I_ESP_Ox8y2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox8y2z_S_vrr;
        I_ESP_Ox7y3z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox7y3z_S_vrr;
        I_ESP_Ox6y4z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox6y4z_S_vrr;
        I_ESP_Ox5y5z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox5y5z_S_vrr;
        I_ESP_Ox4y6z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox4y6z_S_vrr;
        I_ESP_Ox3y7z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox3y7z_S_vrr;
        I_ESP_Ox2y8z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox2y8z_S_vrr;
        I_ESP_Oxy9z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Oxy9z_S_vrr;
        I_ESP_Ox10z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Ox10z_S_vrr;
        I_ESP_O11y_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O11y_S_vrr;
        I_ESP_O10yz_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O10yz_S_vrr;
        I_ESP_O9y2z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O9y2z_S_vrr;
        I_ESP_O8y3z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O8y3z_S_vrr;
        I_ESP_O7y4z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O7y4z_S_vrr;
        I_ESP_O6y5z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O6y5z_S_vrr;
        I_ESP_O5y6z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O5y6z_S_vrr;
        I_ESP_O4y7z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O4y7z_S_vrr;
        I_ESP_O3y8z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O3y8z_S_vrr;
        I_ESP_O2y9z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O2y9z_S_vrr;
        I_ESP_Oy10z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_Oy10z_S_vrr;
        I_ESP_O11z_S_a += SQ_ESP_O_S_a_coefs*I_ESP_O11z_S_vrr;

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
         * shell quartet name: SQ_ESP_M_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_M9x_S += I_ESP_M9x_S_vrr;
        I_ESP_M8xy_S += I_ESP_M8xy_S_vrr;
        I_ESP_M8xz_S += I_ESP_M8xz_S_vrr;
        I_ESP_M7x2y_S += I_ESP_M7x2y_S_vrr;
        I_ESP_M7xyz_S += I_ESP_M7xyz_S_vrr;
        I_ESP_M7x2z_S += I_ESP_M7x2z_S_vrr;
        I_ESP_M6x3y_S += I_ESP_M6x3y_S_vrr;
        I_ESP_M6x2yz_S += I_ESP_M6x2yz_S_vrr;
        I_ESP_M6xy2z_S += I_ESP_M6xy2z_S_vrr;
        I_ESP_M6x3z_S += I_ESP_M6x3z_S_vrr;
        I_ESP_M5x4y_S += I_ESP_M5x4y_S_vrr;
        I_ESP_M5x3yz_S += I_ESP_M5x3yz_S_vrr;
        I_ESP_M5x2y2z_S += I_ESP_M5x2y2z_S_vrr;
        I_ESP_M5xy3z_S += I_ESP_M5xy3z_S_vrr;
        I_ESP_M5x4z_S += I_ESP_M5x4z_S_vrr;
        I_ESP_M4x5y_S += I_ESP_M4x5y_S_vrr;
        I_ESP_M4x4yz_S += I_ESP_M4x4yz_S_vrr;
        I_ESP_M4x3y2z_S += I_ESP_M4x3y2z_S_vrr;
        I_ESP_M4x2y3z_S += I_ESP_M4x2y3z_S_vrr;
        I_ESP_M4xy4z_S += I_ESP_M4xy4z_S_vrr;
        I_ESP_M4x5z_S += I_ESP_M4x5z_S_vrr;
        I_ESP_M3x6y_S += I_ESP_M3x6y_S_vrr;
        I_ESP_M3x5yz_S += I_ESP_M3x5yz_S_vrr;
        I_ESP_M3x4y2z_S += I_ESP_M3x4y2z_S_vrr;
        I_ESP_M3x3y3z_S += I_ESP_M3x3y3z_S_vrr;
        I_ESP_M3x2y4z_S += I_ESP_M3x2y4z_S_vrr;
        I_ESP_M3xy5z_S += I_ESP_M3xy5z_S_vrr;
        I_ESP_M3x6z_S += I_ESP_M3x6z_S_vrr;
        I_ESP_M2x7y_S += I_ESP_M2x7y_S_vrr;
        I_ESP_M2x6yz_S += I_ESP_M2x6yz_S_vrr;
        I_ESP_M2x5y2z_S += I_ESP_M2x5y2z_S_vrr;
        I_ESP_M2x4y3z_S += I_ESP_M2x4y3z_S_vrr;
        I_ESP_M2x3y4z_S += I_ESP_M2x3y4z_S_vrr;
        I_ESP_M2x2y5z_S += I_ESP_M2x2y5z_S_vrr;
        I_ESP_M2xy6z_S += I_ESP_M2xy6z_S_vrr;
        I_ESP_M2x7z_S += I_ESP_M2x7z_S_vrr;
        I_ESP_Mx8y_S += I_ESP_Mx8y_S_vrr;
        I_ESP_Mx7yz_S += I_ESP_Mx7yz_S_vrr;
        I_ESP_Mx6y2z_S += I_ESP_Mx6y2z_S_vrr;
        I_ESP_Mx5y3z_S += I_ESP_Mx5y3z_S_vrr;
        I_ESP_Mx4y4z_S += I_ESP_Mx4y4z_S_vrr;
        I_ESP_Mx3y5z_S += I_ESP_Mx3y5z_S_vrr;
        I_ESP_Mx2y6z_S += I_ESP_Mx2y6z_S_vrr;
        I_ESP_Mxy7z_S += I_ESP_Mxy7z_S_vrr;
        I_ESP_Mx8z_S += I_ESP_Mx8z_S_vrr;
        I_ESP_M9y_S += I_ESP_M9y_S_vrr;
        I_ESP_M8yz_S += I_ESP_M8yz_S_vrr;
        I_ESP_M7y2z_S += I_ESP_M7y2z_S_vrr;
        I_ESP_M6y3z_S += I_ESP_M6y3z_S_vrr;
        I_ESP_M5y4z_S += I_ESP_M5y4z_S_vrr;
        I_ESP_M4y5z_S += I_ESP_M4y5z_S_vrr;
        I_ESP_M3y6z_S += I_ESP_M3y6z_S_vrr;
        I_ESP_M2y7z_S += I_ESP_M2y7z_S_vrr;
        I_ESP_My8z_S += I_ESP_My8z_S_vrr;
        I_ESP_M9z_S += I_ESP_M9z_S_vrr;

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
     * totally 0 integrals are omitted 
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
    Double I_ESP_I6x_Py = I_ESP_K6xy_S+ABY*I_ESP_I6x_S;
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
    Double I_ESP_I6x_Pz = I_ESP_K6xz_S+ABZ*I_ESP_I6x_S;
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
    Double I_ESP_I6y_Pz = I_ESP_K6yz_S+ABZ*I_ESP_I6y_S;
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
     * totally 60 integrals are omitted 
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
     * totally 13 integrals are omitted 
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
    Double I_ESP_K6yz_Px = I_ESP_Lx6yz_S+ABX*I_ESP_K6yz_S;
    Double I_ESP_K5y2z_Px = I_ESP_Lx5y2z_S+ABX*I_ESP_K5y2z_S;
    Double I_ESP_K4y3z_Px = I_ESP_Lx4y3z_S+ABX*I_ESP_K4y3z_S;
    Double I_ESP_K3y4z_Px = I_ESP_Lx3y4z_S+ABX*I_ESP_K3y4z_S;
    Double I_ESP_K2y5z_Px = I_ESP_Lx2y5z_S+ABX*I_ESP_K2y5z_S;
    Double I_ESP_Ky6z_Px = I_ESP_Lxy6z_S+ABX*I_ESP_Ky6z_S;
    Double I_ESP_K6xy_Py = I_ESP_L6x2y_S+ABY*I_ESP_K6xy_S;
    Double I_ESP_K5x2y_Py = I_ESP_L5x3y_S+ABY*I_ESP_K5x2y_S;
    Double I_ESP_K5xyz_Py = I_ESP_L5x2yz_S+ABY*I_ESP_K5xyz_S;
    Double I_ESP_K5x2z_Py = I_ESP_L5xy2z_S+ABY*I_ESP_K5x2z_S;
    Double I_ESP_K4x3y_Py = I_ESP_L4x4y_S+ABY*I_ESP_K4x3y_S;
    Double I_ESP_K4x2yz_Py = I_ESP_L4x3yz_S+ABY*I_ESP_K4x2yz_S;
    Double I_ESP_K4xy2z_Py = I_ESP_L4x2y2z_S+ABY*I_ESP_K4xy2z_S;
    Double I_ESP_K4x3z_Py = I_ESP_L4xy3z_S+ABY*I_ESP_K4x3z_S;
    Double I_ESP_K3x4y_Py = I_ESP_L3x5y_S+ABY*I_ESP_K3x4y_S;
    Double I_ESP_K3x3yz_Py = I_ESP_L3x4yz_S+ABY*I_ESP_K3x3yz_S;
    Double I_ESP_K3x2y2z_Py = I_ESP_L3x3y2z_S+ABY*I_ESP_K3x2y2z_S;
    Double I_ESP_K3xy3z_Py = I_ESP_L3x2y3z_S+ABY*I_ESP_K3xy3z_S;
    Double I_ESP_K3x4z_Py = I_ESP_L3xy4z_S+ABY*I_ESP_K3x4z_S;
    Double I_ESP_K2x5y_Py = I_ESP_L2x6y_S+ABY*I_ESP_K2x5y_S;
    Double I_ESP_K2x4yz_Py = I_ESP_L2x5yz_S+ABY*I_ESP_K2x4yz_S;
    Double I_ESP_K2x3y2z_Py = I_ESP_L2x4y2z_S+ABY*I_ESP_K2x3y2z_S;
    Double I_ESP_K2x2y3z_Py = I_ESP_L2x3y3z_S+ABY*I_ESP_K2x2y3z_S;
    Double I_ESP_K2xy4z_Py = I_ESP_L2x2y4z_S+ABY*I_ESP_K2xy4z_S;
    Double I_ESP_K2x5z_Py = I_ESP_L2xy5z_S+ABY*I_ESP_K2x5z_S;
    Double I_ESP_Kx6y_Py = I_ESP_Lx7y_S+ABY*I_ESP_Kx6y_S;
    Double I_ESP_Kx5yz_Py = I_ESP_Lx6yz_S+ABY*I_ESP_Kx5yz_S;
    Double I_ESP_Kx4y2z_Py = I_ESP_Lx5y2z_S+ABY*I_ESP_Kx4y2z_S;
    Double I_ESP_Kx3y3z_Py = I_ESP_Lx4y3z_S+ABY*I_ESP_Kx3y3z_S;
    Double I_ESP_Kx2y4z_Py = I_ESP_Lx3y4z_S+ABY*I_ESP_Kx2y4z_S;
    Double I_ESP_Kxy5z_Py = I_ESP_Lx2y5z_S+ABY*I_ESP_Kxy5z_S;
    Double I_ESP_Kx6z_Py = I_ESP_Lxy6z_S+ABY*I_ESP_Kx6z_S;
    Double I_ESP_K7y_Py = I_ESP_L8y_S+ABY*I_ESP_K7y_S;
    Double I_ESP_K6yz_Py = I_ESP_L7yz_S+ABY*I_ESP_K6yz_S;
    Double I_ESP_K5y2z_Py = I_ESP_L6y2z_S+ABY*I_ESP_K5y2z_S;
    Double I_ESP_K4y3z_Py = I_ESP_L5y3z_S+ABY*I_ESP_K4y3z_S;
    Double I_ESP_K3y4z_Py = I_ESP_L4y4z_S+ABY*I_ESP_K3y4z_S;
    Double I_ESP_K2y5z_Py = I_ESP_L3y5z_S+ABY*I_ESP_K2y5z_S;
    Double I_ESP_Ky6z_Py = I_ESP_L2y6z_S+ABY*I_ESP_Ky6z_S;
    Double I_ESP_K6xz_Pz = I_ESP_L6x2z_S+ABZ*I_ESP_K6xz_S;
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
    Double I_ESP_K6yz_Pz = I_ESP_L6y2z_S+ABZ*I_ESP_K6yz_S;
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
     * totally 84 integrals are omitted 
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
    Double I_ESP_I6x_D2y = I_ESP_K6xy_Py+ABY*I_ESP_I6x_Py;
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
    Double I_ESP_I6x_D2z = I_ESP_K6xz_Pz+ABZ*I_ESP_I6x_Pz;
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
    Double I_ESP_I6y_D2z = I_ESP_K6yz_Pz+ABZ*I_ESP_I6y_Pz;
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
     * totally 87 integrals are omitted 
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
    Double I_ESP_H4xy_F2xz = I_ESP_I4xyz_D2x+ABZ*I_ESP_H4xy_D2x;
    Double I_ESP_H4xz_F2xz = I_ESP_I4x2z_D2x+ABZ*I_ESP_H4xz_D2x;
    Double I_ESP_H3x2y_F2xz = I_ESP_I3x2yz_D2x+ABZ*I_ESP_H3x2y_D2x;
    Double I_ESP_H3xyz_F2xz = I_ESP_I3xy2z_D2x+ABZ*I_ESP_H3xyz_D2x;
    Double I_ESP_H3x2z_F2xz = I_ESP_I3x3z_D2x+ABZ*I_ESP_H3x2z_D2x;
    Double I_ESP_H2x3y_F2xz = I_ESP_I2x3yz_D2x+ABZ*I_ESP_H2x3y_D2x;
    Double I_ESP_H2x2yz_F2xz = I_ESP_I2x2y2z_D2x+ABZ*I_ESP_H2x2yz_D2x;
    Double I_ESP_H2xy2z_F2xz = I_ESP_I2xy3z_D2x+ABZ*I_ESP_H2xy2z_D2x;
    Double I_ESP_H2x3z_F2xz = I_ESP_I2x4z_D2x+ABZ*I_ESP_H2x3z_D2x;
    Double I_ESP_Hx4y_F2xz = I_ESP_Ix4yz_D2x+ABZ*I_ESP_Hx4y_D2x;
    Double I_ESP_Hx3yz_F2xz = I_ESP_Ix3y2z_D2x+ABZ*I_ESP_Hx3yz_D2x;
    Double I_ESP_Hx2y2z_F2xz = I_ESP_Ix2y3z_D2x+ABZ*I_ESP_Hx2y2z_D2x;
    Double I_ESP_Hxy3z_F2xz = I_ESP_Ixy4z_D2x+ABZ*I_ESP_Hxy3z_D2x;
    Double I_ESP_Hx4z_F2xz = I_ESP_Ix5z_D2x+ABZ*I_ESP_Hx4z_D2x;
    Double I_ESP_H5y_F2xz = I_ESP_I5yz_D2x+ABZ*I_ESP_H5y_D2x;
    Double I_ESP_H4yz_F2xz = I_ESP_I4y2z_D2x+ABZ*I_ESP_H4yz_D2x;
    Double I_ESP_H3y2z_F2xz = I_ESP_I3y3z_D2x+ABZ*I_ESP_H3y2z_D2x;
    Double I_ESP_H2y3z_F2xz = I_ESP_I2y4z_D2x+ABZ*I_ESP_H2y3z_D2x;
    Double I_ESP_Hy4z_F2xz = I_ESP_Iy5z_D2x+ABZ*I_ESP_Hy4z_D2x;
    Double I_ESP_H5z_F2xz = I_ESP_I6z_D2x+ABZ*I_ESP_H5z_D2x;
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
    Double I_ESP_H5x_F2yz = I_ESP_I5xz_D2y+ABZ*I_ESP_H5x_D2y;
    Double I_ESP_H4xy_F2yz = I_ESP_I4xyz_D2y+ABZ*I_ESP_H4xy_D2y;
    Double I_ESP_H4xz_F2yz = I_ESP_I4x2z_D2y+ABZ*I_ESP_H4xz_D2y;
    Double I_ESP_H3x2y_F2yz = I_ESP_I3x2yz_D2y+ABZ*I_ESP_H3x2y_D2y;
    Double I_ESP_H3xyz_F2yz = I_ESP_I3xy2z_D2y+ABZ*I_ESP_H3xyz_D2y;
    Double I_ESP_H3x2z_F2yz = I_ESP_I3x3z_D2y+ABZ*I_ESP_H3x2z_D2y;
    Double I_ESP_H2x3y_F2yz = I_ESP_I2x3yz_D2y+ABZ*I_ESP_H2x3y_D2y;
    Double I_ESP_H2x2yz_F2yz = I_ESP_I2x2y2z_D2y+ABZ*I_ESP_H2x2yz_D2y;
    Double I_ESP_H2xy2z_F2yz = I_ESP_I2xy3z_D2y+ABZ*I_ESP_H2xy2z_D2y;
    Double I_ESP_H2x3z_F2yz = I_ESP_I2x4z_D2y+ABZ*I_ESP_H2x3z_D2y;
    Double I_ESP_Hx4y_F2yz = I_ESP_Ix4yz_D2y+ABZ*I_ESP_Hx4y_D2y;
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
     * totally 45 integrals are omitted 
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
     * shell quartet name: SQ_ESP_L_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 40 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_M_S
     * RHS shell quartet name: SQ_ESP_L_S
     ************************************************************/
    Double I_ESP_L8x_Px = I_ESP_M9x_S+ABX*I_ESP_L8x_S;
    Double I_ESP_L7xy_Px = I_ESP_M8xy_S+ABX*I_ESP_L7xy_S;
    Double I_ESP_L7xz_Px = I_ESP_M8xz_S+ABX*I_ESP_L7xz_S;
    Double I_ESP_L6x2y_Px = I_ESP_M7x2y_S+ABX*I_ESP_L6x2y_S;
    Double I_ESP_L6xyz_Px = I_ESP_M7xyz_S+ABX*I_ESP_L6xyz_S;
    Double I_ESP_L6x2z_Px = I_ESP_M7x2z_S+ABX*I_ESP_L6x2z_S;
    Double I_ESP_L5x3y_Px = I_ESP_M6x3y_S+ABX*I_ESP_L5x3y_S;
    Double I_ESP_L5x2yz_Px = I_ESP_M6x2yz_S+ABX*I_ESP_L5x2yz_S;
    Double I_ESP_L5xy2z_Px = I_ESP_M6xy2z_S+ABX*I_ESP_L5xy2z_S;
    Double I_ESP_L5x3z_Px = I_ESP_M6x3z_S+ABX*I_ESP_L5x3z_S;
    Double I_ESP_L4x4y_Px = I_ESP_M5x4y_S+ABX*I_ESP_L4x4y_S;
    Double I_ESP_L4x3yz_Px = I_ESP_M5x3yz_S+ABX*I_ESP_L4x3yz_S;
    Double I_ESP_L4x2y2z_Px = I_ESP_M5x2y2z_S+ABX*I_ESP_L4x2y2z_S;
    Double I_ESP_L4xy3z_Px = I_ESP_M5xy3z_S+ABX*I_ESP_L4xy3z_S;
    Double I_ESP_L4x4z_Px = I_ESP_M5x4z_S+ABX*I_ESP_L4x4z_S;
    Double I_ESP_L3x5y_Px = I_ESP_M4x5y_S+ABX*I_ESP_L3x5y_S;
    Double I_ESP_L3x4yz_Px = I_ESP_M4x4yz_S+ABX*I_ESP_L3x4yz_S;
    Double I_ESP_L3x3y2z_Px = I_ESP_M4x3y2z_S+ABX*I_ESP_L3x3y2z_S;
    Double I_ESP_L3x2y3z_Px = I_ESP_M4x2y3z_S+ABX*I_ESP_L3x2y3z_S;
    Double I_ESP_L3xy4z_Px = I_ESP_M4xy4z_S+ABX*I_ESP_L3xy4z_S;
    Double I_ESP_L3x5z_Px = I_ESP_M4x5z_S+ABX*I_ESP_L3x5z_S;
    Double I_ESP_L2x6y_Px = I_ESP_M3x6y_S+ABX*I_ESP_L2x6y_S;
    Double I_ESP_L2x5yz_Px = I_ESP_M3x5yz_S+ABX*I_ESP_L2x5yz_S;
    Double I_ESP_L2x4y2z_Px = I_ESP_M3x4y2z_S+ABX*I_ESP_L2x4y2z_S;
    Double I_ESP_L2x3y3z_Px = I_ESP_M3x3y3z_S+ABX*I_ESP_L2x3y3z_S;
    Double I_ESP_L2x2y4z_Px = I_ESP_M3x2y4z_S+ABX*I_ESP_L2x2y4z_S;
    Double I_ESP_L2xy5z_Px = I_ESP_M3xy5z_S+ABX*I_ESP_L2xy5z_S;
    Double I_ESP_L2x6z_Px = I_ESP_M3x6z_S+ABX*I_ESP_L2x6z_S;
    Double I_ESP_Lx6yz_Px = I_ESP_M2x6yz_S+ABX*I_ESP_Lx6yz_S;
    Double I_ESP_Lx5y2z_Px = I_ESP_M2x5y2z_S+ABX*I_ESP_Lx5y2z_S;
    Double I_ESP_Lx4y3z_Px = I_ESP_M2x4y3z_S+ABX*I_ESP_Lx4y3z_S;
    Double I_ESP_Lx3y4z_Px = I_ESP_M2x3y4z_S+ABX*I_ESP_Lx3y4z_S;
    Double I_ESP_Lx2y5z_Px = I_ESP_M2x2y5z_S+ABX*I_ESP_Lx2y5z_S;
    Double I_ESP_Lxy6z_Px = I_ESP_M2xy6z_S+ABX*I_ESP_Lxy6z_S;
    Double I_ESP_L6x2y_Py = I_ESP_M6x3y_S+ABY*I_ESP_L6x2y_S;
    Double I_ESP_L5x3y_Py = I_ESP_M5x4y_S+ABY*I_ESP_L5x3y_S;
    Double I_ESP_L5x2yz_Py = I_ESP_M5x3yz_S+ABY*I_ESP_L5x2yz_S;
    Double I_ESP_L5xy2z_Py = I_ESP_M5x2y2z_S+ABY*I_ESP_L5xy2z_S;
    Double I_ESP_L4x4y_Py = I_ESP_M4x5y_S+ABY*I_ESP_L4x4y_S;
    Double I_ESP_L4x3yz_Py = I_ESP_M4x4yz_S+ABY*I_ESP_L4x3yz_S;
    Double I_ESP_L4x2y2z_Py = I_ESP_M4x3y2z_S+ABY*I_ESP_L4x2y2z_S;
    Double I_ESP_L4xy3z_Py = I_ESP_M4x2y3z_S+ABY*I_ESP_L4xy3z_S;
    Double I_ESP_L3x5y_Py = I_ESP_M3x6y_S+ABY*I_ESP_L3x5y_S;
    Double I_ESP_L3x4yz_Py = I_ESP_M3x5yz_S+ABY*I_ESP_L3x4yz_S;
    Double I_ESP_L3x3y2z_Py = I_ESP_M3x4y2z_S+ABY*I_ESP_L3x3y2z_S;
    Double I_ESP_L3x2y3z_Py = I_ESP_M3x3y3z_S+ABY*I_ESP_L3x2y3z_S;
    Double I_ESP_L3xy4z_Py = I_ESP_M3x2y4z_S+ABY*I_ESP_L3xy4z_S;
    Double I_ESP_L2x6y_Py = I_ESP_M2x7y_S+ABY*I_ESP_L2x6y_S;
    Double I_ESP_L2x5yz_Py = I_ESP_M2x6yz_S+ABY*I_ESP_L2x5yz_S;
    Double I_ESP_L2x4y2z_Py = I_ESP_M2x5y2z_S+ABY*I_ESP_L2x4y2z_S;
    Double I_ESP_L2x3y3z_Py = I_ESP_M2x4y3z_S+ABY*I_ESP_L2x3y3z_S;
    Double I_ESP_L2x2y4z_Py = I_ESP_M2x3y4z_S+ABY*I_ESP_L2x2y4z_S;
    Double I_ESP_L2xy5z_Py = I_ESP_M2x2y5z_S+ABY*I_ESP_L2xy5z_S;
    Double I_ESP_Lx7y_Py = I_ESP_Mx8y_S+ABY*I_ESP_Lx7y_S;
    Double I_ESP_Lx6yz_Py = I_ESP_Mx7yz_S+ABY*I_ESP_Lx6yz_S;
    Double I_ESP_Lx5y2z_Py = I_ESP_Mx6y2z_S+ABY*I_ESP_Lx5y2z_S;
    Double I_ESP_Lx4y3z_Py = I_ESP_Mx5y3z_S+ABY*I_ESP_Lx4y3z_S;
    Double I_ESP_Lx3y4z_Py = I_ESP_Mx4y4z_S+ABY*I_ESP_Lx3y4z_S;
    Double I_ESP_Lx2y5z_Py = I_ESP_Mx3y5z_S+ABY*I_ESP_Lx2y5z_S;
    Double I_ESP_Lxy6z_Py = I_ESP_Mx2y6z_S+ABY*I_ESP_Lxy6z_S;
    Double I_ESP_L8y_Py = I_ESP_M9y_S+ABY*I_ESP_L8y_S;
    Double I_ESP_L7yz_Py = I_ESP_M8yz_S+ABY*I_ESP_L7yz_S;
    Double I_ESP_L6y2z_Py = I_ESP_M7y2z_S+ABY*I_ESP_L6y2z_S;
    Double I_ESP_L5y3z_Py = I_ESP_M6y3z_S+ABY*I_ESP_L5y3z_S;
    Double I_ESP_L4y4z_Py = I_ESP_M5y4z_S+ABY*I_ESP_L4y4z_S;
    Double I_ESP_L3y5z_Py = I_ESP_M4y5z_S+ABY*I_ESP_L3y5z_S;
    Double I_ESP_L2y6z_Py = I_ESP_M3y6z_S+ABY*I_ESP_L2y6z_S;
    Double I_ESP_L6x2z_Pz = I_ESP_M6x3z_S+ABZ*I_ESP_L6x2z_S;
    Double I_ESP_L5xy2z_Pz = I_ESP_M5xy3z_S+ABZ*I_ESP_L5xy2z_S;
    Double I_ESP_L5x3z_Pz = I_ESP_M5x4z_S+ABZ*I_ESP_L5x3z_S;
    Double I_ESP_L4x2y2z_Pz = I_ESP_M4x2y3z_S+ABZ*I_ESP_L4x2y2z_S;
    Double I_ESP_L4xy3z_Pz = I_ESP_M4xy4z_S+ABZ*I_ESP_L4xy3z_S;
    Double I_ESP_L4x4z_Pz = I_ESP_M4x5z_S+ABZ*I_ESP_L4x4z_S;
    Double I_ESP_L3x3y2z_Pz = I_ESP_M3x3y3z_S+ABZ*I_ESP_L3x3y2z_S;
    Double I_ESP_L3x2y3z_Pz = I_ESP_M3x2y4z_S+ABZ*I_ESP_L3x2y3z_S;
    Double I_ESP_L3xy4z_Pz = I_ESP_M3xy5z_S+ABZ*I_ESP_L3xy4z_S;
    Double I_ESP_L3x5z_Pz = I_ESP_M3x6z_S+ABZ*I_ESP_L3x5z_S;
    Double I_ESP_L2x4y2z_Pz = I_ESP_M2x4y3z_S+ABZ*I_ESP_L2x4y2z_S;
    Double I_ESP_L2x3y3z_Pz = I_ESP_M2x3y4z_S+ABZ*I_ESP_L2x3y3z_S;
    Double I_ESP_L2x2y4z_Pz = I_ESP_M2x2y5z_S+ABZ*I_ESP_L2x2y4z_S;
    Double I_ESP_L2xy5z_Pz = I_ESP_M2xy6z_S+ABZ*I_ESP_L2xy5z_S;
    Double I_ESP_L2x6z_Pz = I_ESP_M2x7z_S+ABZ*I_ESP_L2x6z_S;
    Double I_ESP_Lx5y2z_Pz = I_ESP_Mx5y3z_S+ABZ*I_ESP_Lx5y2z_S;
    Double I_ESP_Lx4y3z_Pz = I_ESP_Mx4y4z_S+ABZ*I_ESP_Lx4y3z_S;
    Double I_ESP_Lx3y4z_Pz = I_ESP_Mx3y5z_S+ABZ*I_ESP_Lx3y4z_S;
    Double I_ESP_Lx2y5z_Pz = I_ESP_Mx2y6z_S+ABZ*I_ESP_Lx2y5z_S;
    Double I_ESP_Lxy6z_Pz = I_ESP_Mxy7z_S+ABZ*I_ESP_Lxy6z_S;
    Double I_ESP_Lx7z_Pz = I_ESP_Mx8z_S+ABZ*I_ESP_Lx7z_S;
    Double I_ESP_L6y2z_Pz = I_ESP_M6y3z_S+ABZ*I_ESP_L6y2z_S;
    Double I_ESP_L5y3z_Pz = I_ESP_M5y4z_S+ABZ*I_ESP_L5y3z_S;
    Double I_ESP_L4y4z_Pz = I_ESP_M4y5z_S+ABZ*I_ESP_L4y4z_S;
    Double I_ESP_L3y5z_Pz = I_ESP_M3y6z_S+ABZ*I_ESP_L3y5z_S;
    Double I_ESP_L2y6z_Pz = I_ESP_M2y7z_S+ABZ*I_ESP_L2y6z_S;
    Double I_ESP_Ly7z_Pz = I_ESP_My8z_S+ABZ*I_ESP_Ly7z_S;
    Double I_ESP_L8z_Pz = I_ESP_M9z_S+ABZ*I_ESP_L8z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_K_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 121 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_P
     * RHS shell quartet name: SQ_ESP_K_P
     ************************************************************/
    Double I_ESP_K7x_D2x = I_ESP_L8x_Px+ABX*I_ESP_K7x_Px;
    Double I_ESP_K6xy_D2x = I_ESP_L7xy_Px+ABX*I_ESP_K6xy_Px;
    Double I_ESP_K6xz_D2x = I_ESP_L7xz_Px+ABX*I_ESP_K6xz_Px;
    Double I_ESP_K5x2y_D2x = I_ESP_L6x2y_Px+ABX*I_ESP_K5x2y_Px;
    Double I_ESP_K5xyz_D2x = I_ESP_L6xyz_Px+ABX*I_ESP_K5xyz_Px;
    Double I_ESP_K5x2z_D2x = I_ESP_L6x2z_Px+ABX*I_ESP_K5x2z_Px;
    Double I_ESP_K4x3y_D2x = I_ESP_L5x3y_Px+ABX*I_ESP_K4x3y_Px;
    Double I_ESP_K4x2yz_D2x = I_ESP_L5x2yz_Px+ABX*I_ESP_K4x2yz_Px;
    Double I_ESP_K4xy2z_D2x = I_ESP_L5xy2z_Px+ABX*I_ESP_K4xy2z_Px;
    Double I_ESP_K4x3z_D2x = I_ESP_L5x3z_Px+ABX*I_ESP_K4x3z_Px;
    Double I_ESP_K3x4y_D2x = I_ESP_L4x4y_Px+ABX*I_ESP_K3x4y_Px;
    Double I_ESP_K3x3yz_D2x = I_ESP_L4x3yz_Px+ABX*I_ESP_K3x3yz_Px;
    Double I_ESP_K3x2y2z_D2x = I_ESP_L4x2y2z_Px+ABX*I_ESP_K3x2y2z_Px;
    Double I_ESP_K3xy3z_D2x = I_ESP_L4xy3z_Px+ABX*I_ESP_K3xy3z_Px;
    Double I_ESP_K3x4z_D2x = I_ESP_L4x4z_Px+ABX*I_ESP_K3x4z_Px;
    Double I_ESP_K2x5y_D2x = I_ESP_L3x5y_Px+ABX*I_ESP_K2x5y_Px;
    Double I_ESP_K2x4yz_D2x = I_ESP_L3x4yz_Px+ABX*I_ESP_K2x4yz_Px;
    Double I_ESP_K2x3y2z_D2x = I_ESP_L3x3y2z_Px+ABX*I_ESP_K2x3y2z_Px;
    Double I_ESP_K2x2y3z_D2x = I_ESP_L3x2y3z_Px+ABX*I_ESP_K2x2y3z_Px;
    Double I_ESP_K2xy4z_D2x = I_ESP_L3xy4z_Px+ABX*I_ESP_K2xy4z_Px;
    Double I_ESP_K2x5z_D2x = I_ESP_L3x5z_Px+ABX*I_ESP_K2x5z_Px;
    Double I_ESP_Kx6y_D2x = I_ESP_L2x6y_Px+ABX*I_ESP_Kx6y_Px;
    Double I_ESP_Kx5yz_D2x = I_ESP_L2x5yz_Px+ABX*I_ESP_Kx5yz_Px;
    Double I_ESP_Kx4y2z_D2x = I_ESP_L2x4y2z_Px+ABX*I_ESP_Kx4y2z_Px;
    Double I_ESP_Kx3y3z_D2x = I_ESP_L2x3y3z_Px+ABX*I_ESP_Kx3y3z_Px;
    Double I_ESP_Kx2y4z_D2x = I_ESP_L2x2y4z_Px+ABX*I_ESP_Kx2y4z_Px;
    Double I_ESP_Kxy5z_D2x = I_ESP_L2xy5z_Px+ABX*I_ESP_Kxy5z_Px;
    Double I_ESP_Kx6z_D2x = I_ESP_L2x6z_Px+ABX*I_ESP_Kx6z_Px;
    Double I_ESP_K6yz_D2x = I_ESP_Lx6yz_Px+ABX*I_ESP_K6yz_Px;
    Double I_ESP_K5y2z_D2x = I_ESP_Lx5y2z_Px+ABX*I_ESP_K5y2z_Px;
    Double I_ESP_K4y3z_D2x = I_ESP_Lx4y3z_Px+ABX*I_ESP_K4y3z_Px;
    Double I_ESP_K3y4z_D2x = I_ESP_Lx3y4z_Px+ABX*I_ESP_K3y4z_Px;
    Double I_ESP_K2y5z_D2x = I_ESP_Lx2y5z_Px+ABX*I_ESP_K2y5z_Px;
    Double I_ESP_Ky6z_D2x = I_ESP_Lxy6z_Px+ABX*I_ESP_Ky6z_Px;
    Double I_ESP_K6xy_D2y = I_ESP_L6x2y_Py+ABY*I_ESP_K6xy_Py;
    Double I_ESP_K5x2y_D2y = I_ESP_L5x3y_Py+ABY*I_ESP_K5x2y_Py;
    Double I_ESP_K5xyz_D2y = I_ESP_L5x2yz_Py+ABY*I_ESP_K5xyz_Py;
    Double I_ESP_K5x2z_D2y = I_ESP_L5xy2z_Py+ABY*I_ESP_K5x2z_Py;
    Double I_ESP_K4x3y_D2y = I_ESP_L4x4y_Py+ABY*I_ESP_K4x3y_Py;
    Double I_ESP_K4x2yz_D2y = I_ESP_L4x3yz_Py+ABY*I_ESP_K4x2yz_Py;
    Double I_ESP_K4xy2z_D2y = I_ESP_L4x2y2z_Py+ABY*I_ESP_K4xy2z_Py;
    Double I_ESP_K4x3z_D2y = I_ESP_L4xy3z_Py+ABY*I_ESP_K4x3z_Py;
    Double I_ESP_K3x4y_D2y = I_ESP_L3x5y_Py+ABY*I_ESP_K3x4y_Py;
    Double I_ESP_K3x3yz_D2y = I_ESP_L3x4yz_Py+ABY*I_ESP_K3x3yz_Py;
    Double I_ESP_K3x2y2z_D2y = I_ESP_L3x3y2z_Py+ABY*I_ESP_K3x2y2z_Py;
    Double I_ESP_K3xy3z_D2y = I_ESP_L3x2y3z_Py+ABY*I_ESP_K3xy3z_Py;
    Double I_ESP_K3x4z_D2y = I_ESP_L3xy4z_Py+ABY*I_ESP_K3x4z_Py;
    Double I_ESP_K2x5y_D2y = I_ESP_L2x6y_Py+ABY*I_ESP_K2x5y_Py;
    Double I_ESP_K2x4yz_D2y = I_ESP_L2x5yz_Py+ABY*I_ESP_K2x4yz_Py;
    Double I_ESP_K2x3y2z_D2y = I_ESP_L2x4y2z_Py+ABY*I_ESP_K2x3y2z_Py;
    Double I_ESP_K2x2y3z_D2y = I_ESP_L2x3y3z_Py+ABY*I_ESP_K2x2y3z_Py;
    Double I_ESP_K2xy4z_D2y = I_ESP_L2x2y4z_Py+ABY*I_ESP_K2xy4z_Py;
    Double I_ESP_K2x5z_D2y = I_ESP_L2xy5z_Py+ABY*I_ESP_K2x5z_Py;
    Double I_ESP_Kx6y_D2y = I_ESP_Lx7y_Py+ABY*I_ESP_Kx6y_Py;
    Double I_ESP_Kx5yz_D2y = I_ESP_Lx6yz_Py+ABY*I_ESP_Kx5yz_Py;
    Double I_ESP_Kx4y2z_D2y = I_ESP_Lx5y2z_Py+ABY*I_ESP_Kx4y2z_Py;
    Double I_ESP_Kx3y3z_D2y = I_ESP_Lx4y3z_Py+ABY*I_ESP_Kx3y3z_Py;
    Double I_ESP_Kx2y4z_D2y = I_ESP_Lx3y4z_Py+ABY*I_ESP_Kx2y4z_Py;
    Double I_ESP_Kxy5z_D2y = I_ESP_Lx2y5z_Py+ABY*I_ESP_Kxy5z_Py;
    Double I_ESP_Kx6z_D2y = I_ESP_Lxy6z_Py+ABY*I_ESP_Kx6z_Py;
    Double I_ESP_K7y_D2y = I_ESP_L8y_Py+ABY*I_ESP_K7y_Py;
    Double I_ESP_K6yz_D2y = I_ESP_L7yz_Py+ABY*I_ESP_K6yz_Py;
    Double I_ESP_K5y2z_D2y = I_ESP_L6y2z_Py+ABY*I_ESP_K5y2z_Py;
    Double I_ESP_K4y3z_D2y = I_ESP_L5y3z_Py+ABY*I_ESP_K4y3z_Py;
    Double I_ESP_K3y4z_D2y = I_ESP_L4y4z_Py+ABY*I_ESP_K3y4z_Py;
    Double I_ESP_K2y5z_D2y = I_ESP_L3y5z_Py+ABY*I_ESP_K2y5z_Py;
    Double I_ESP_Ky6z_D2y = I_ESP_L2y6z_Py+ABY*I_ESP_Ky6z_Py;
    Double I_ESP_K6xz_D2z = I_ESP_L6x2z_Pz+ABZ*I_ESP_K6xz_Pz;
    Double I_ESP_K5xyz_D2z = I_ESP_L5xy2z_Pz+ABZ*I_ESP_K5xyz_Pz;
    Double I_ESP_K5x2z_D2z = I_ESP_L5x3z_Pz+ABZ*I_ESP_K5x2z_Pz;
    Double I_ESP_K4x2yz_D2z = I_ESP_L4x2y2z_Pz+ABZ*I_ESP_K4x2yz_Pz;
    Double I_ESP_K4xy2z_D2z = I_ESP_L4xy3z_Pz+ABZ*I_ESP_K4xy2z_Pz;
    Double I_ESP_K4x3z_D2z = I_ESP_L4x4z_Pz+ABZ*I_ESP_K4x3z_Pz;
    Double I_ESP_K3x3yz_D2z = I_ESP_L3x3y2z_Pz+ABZ*I_ESP_K3x3yz_Pz;
    Double I_ESP_K3x2y2z_D2z = I_ESP_L3x2y3z_Pz+ABZ*I_ESP_K3x2y2z_Pz;
    Double I_ESP_K3xy3z_D2z = I_ESP_L3xy4z_Pz+ABZ*I_ESP_K3xy3z_Pz;
    Double I_ESP_K3x4z_D2z = I_ESP_L3x5z_Pz+ABZ*I_ESP_K3x4z_Pz;
    Double I_ESP_K2x4yz_D2z = I_ESP_L2x4y2z_Pz+ABZ*I_ESP_K2x4yz_Pz;
    Double I_ESP_K2x3y2z_D2z = I_ESP_L2x3y3z_Pz+ABZ*I_ESP_K2x3y2z_Pz;
    Double I_ESP_K2x2y3z_D2z = I_ESP_L2x2y4z_Pz+ABZ*I_ESP_K2x2y3z_Pz;
    Double I_ESP_K2xy4z_D2z = I_ESP_L2xy5z_Pz+ABZ*I_ESP_K2xy4z_Pz;
    Double I_ESP_K2x5z_D2z = I_ESP_L2x6z_Pz+ABZ*I_ESP_K2x5z_Pz;
    Double I_ESP_Kx5yz_D2z = I_ESP_Lx5y2z_Pz+ABZ*I_ESP_Kx5yz_Pz;
    Double I_ESP_Kx4y2z_D2z = I_ESP_Lx4y3z_Pz+ABZ*I_ESP_Kx4y2z_Pz;
    Double I_ESP_Kx3y3z_D2z = I_ESP_Lx3y4z_Pz+ABZ*I_ESP_Kx3y3z_Pz;
    Double I_ESP_Kx2y4z_D2z = I_ESP_Lx2y5z_Pz+ABZ*I_ESP_Kx2y4z_Pz;
    Double I_ESP_Kxy5z_D2z = I_ESP_Lxy6z_Pz+ABZ*I_ESP_Kxy5z_Pz;
    Double I_ESP_Kx6z_D2z = I_ESP_Lx7z_Pz+ABZ*I_ESP_Kx6z_Pz;
    Double I_ESP_K6yz_D2z = I_ESP_L6y2z_Pz+ABZ*I_ESP_K6yz_Pz;
    Double I_ESP_K5y2z_D2z = I_ESP_L5y3z_Pz+ABZ*I_ESP_K5y2z_Pz;
    Double I_ESP_K4y3z_D2z = I_ESP_L4y4z_Pz+ABZ*I_ESP_K4y3z_Pz;
    Double I_ESP_K3y4z_D2z = I_ESP_L3y5z_Pz+ABZ*I_ESP_K3y4z_Pz;
    Double I_ESP_K2y5z_D2z = I_ESP_L2y6z_Pz+ABZ*I_ESP_K2y5z_Pz;
    Double I_ESP_Ky6z_D2z = I_ESP_Ly7z_Pz+ABZ*I_ESP_Ky6z_Pz;
    Double I_ESP_K7z_D2z = I_ESP_L8z_Pz+ABZ*I_ESP_K7z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_I_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 151 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_D
     * RHS shell quartet name: SQ_ESP_I_D
     ************************************************************/
    Double I_ESP_I6x_F3x = I_ESP_K7x_D2x+ABX*I_ESP_I6x_D2x;
    Double I_ESP_I5xy_F3x = I_ESP_K6xy_D2x+ABX*I_ESP_I5xy_D2x;
    Double I_ESP_I5xz_F3x = I_ESP_K6xz_D2x+ABX*I_ESP_I5xz_D2x;
    Double I_ESP_I4x2y_F3x = I_ESP_K5x2y_D2x+ABX*I_ESP_I4x2y_D2x;
    Double I_ESP_I4xyz_F3x = I_ESP_K5xyz_D2x+ABX*I_ESP_I4xyz_D2x;
    Double I_ESP_I4x2z_F3x = I_ESP_K5x2z_D2x+ABX*I_ESP_I4x2z_D2x;
    Double I_ESP_I3x3y_F3x = I_ESP_K4x3y_D2x+ABX*I_ESP_I3x3y_D2x;
    Double I_ESP_I3x2yz_F3x = I_ESP_K4x2yz_D2x+ABX*I_ESP_I3x2yz_D2x;
    Double I_ESP_I3xy2z_F3x = I_ESP_K4xy2z_D2x+ABX*I_ESP_I3xy2z_D2x;
    Double I_ESP_I3x3z_F3x = I_ESP_K4x3z_D2x+ABX*I_ESP_I3x3z_D2x;
    Double I_ESP_I2x4y_F3x = I_ESP_K3x4y_D2x+ABX*I_ESP_I2x4y_D2x;
    Double I_ESP_I2x3yz_F3x = I_ESP_K3x3yz_D2x+ABX*I_ESP_I2x3yz_D2x;
    Double I_ESP_I2x2y2z_F3x = I_ESP_K3x2y2z_D2x+ABX*I_ESP_I2x2y2z_D2x;
    Double I_ESP_I2xy3z_F3x = I_ESP_K3xy3z_D2x+ABX*I_ESP_I2xy3z_D2x;
    Double I_ESP_I2x4z_F3x = I_ESP_K3x4z_D2x+ABX*I_ESP_I2x4z_D2x;
    Double I_ESP_Ix5y_F3x = I_ESP_K2x5y_D2x+ABX*I_ESP_Ix5y_D2x;
    Double I_ESP_Ix4yz_F3x = I_ESP_K2x4yz_D2x+ABX*I_ESP_Ix4yz_D2x;
    Double I_ESP_Ix3y2z_F3x = I_ESP_K2x3y2z_D2x+ABX*I_ESP_Ix3y2z_D2x;
    Double I_ESP_Ix2y3z_F3x = I_ESP_K2x2y3z_D2x+ABX*I_ESP_Ix2y3z_D2x;
    Double I_ESP_Ixy4z_F3x = I_ESP_K2xy4z_D2x+ABX*I_ESP_Ixy4z_D2x;
    Double I_ESP_Ix5z_F3x = I_ESP_K2x5z_D2x+ABX*I_ESP_Ix5z_D2x;
    Double I_ESP_I6y_F3x = I_ESP_Kx6y_D2x+ABX*I_ESP_I6y_D2x;
    Double I_ESP_I5yz_F3x = I_ESP_Kx5yz_D2x+ABX*I_ESP_I5yz_D2x;
    Double I_ESP_I4y2z_F3x = I_ESP_Kx4y2z_D2x+ABX*I_ESP_I4y2z_D2x;
    Double I_ESP_I3y3z_F3x = I_ESP_Kx3y3z_D2x+ABX*I_ESP_I3y3z_D2x;
    Double I_ESP_I2y4z_F3x = I_ESP_Kx2y4z_D2x+ABX*I_ESP_I2y4z_D2x;
    Double I_ESP_Iy5z_F3x = I_ESP_Kxy5z_D2x+ABX*I_ESP_Iy5z_D2x;
    Double I_ESP_I6z_F3x = I_ESP_Kx6z_D2x+ABX*I_ESP_I6z_D2x;
    Double I_ESP_I4xyz_F2xy = I_ESP_K4x2yz_D2x+ABY*I_ESP_I4xyz_D2x;
    Double I_ESP_I3x2yz_F2xy = I_ESP_K3x3yz_D2x+ABY*I_ESP_I3x2yz_D2x;
    Double I_ESP_I3xy2z_F2xy = I_ESP_K3x2y2z_D2x+ABY*I_ESP_I3xy2z_D2x;
    Double I_ESP_I2x3yz_F2xy = I_ESP_K2x4yz_D2x+ABY*I_ESP_I2x3yz_D2x;
    Double I_ESP_I2x2y2z_F2xy = I_ESP_K2x3y2z_D2x+ABY*I_ESP_I2x2y2z_D2x;
    Double I_ESP_I2xy3z_F2xy = I_ESP_K2x2y3z_D2x+ABY*I_ESP_I2xy3z_D2x;
    Double I_ESP_Ix4yz_F2xy = I_ESP_Kx5yz_D2x+ABY*I_ESP_Ix4yz_D2x;
    Double I_ESP_Ix3y2z_F2xy = I_ESP_Kx4y2z_D2x+ABY*I_ESP_Ix3y2z_D2x;
    Double I_ESP_Ix2y3z_F2xy = I_ESP_Kx3y3z_D2x+ABY*I_ESP_Ix2y3z_D2x;
    Double I_ESP_Ixy4z_F2xy = I_ESP_Kx2y4z_D2x+ABY*I_ESP_Ixy4z_D2x;
    Double I_ESP_I5yz_F2xy = I_ESP_K6yz_D2x+ABY*I_ESP_I5yz_D2x;
    Double I_ESP_I4y2z_F2xy = I_ESP_K5y2z_D2x+ABY*I_ESP_I4y2z_D2x;
    Double I_ESP_I3y3z_F2xy = I_ESP_K4y3z_D2x+ABY*I_ESP_I3y3z_D2x;
    Double I_ESP_I2y4z_F2xy = I_ESP_K3y4z_D2x+ABY*I_ESP_I2y4z_D2x;
    Double I_ESP_Iy5z_F2xy = I_ESP_K2y5z_D2x+ABY*I_ESP_Iy5z_D2x;
    Double I_ESP_I4xyz_F2xz = I_ESP_K4xy2z_D2x+ABZ*I_ESP_I4xyz_D2x;
    Double I_ESP_I3x2yz_F2xz = I_ESP_K3x2y2z_D2x+ABZ*I_ESP_I3x2yz_D2x;
    Double I_ESP_I3xy2z_F2xz = I_ESP_K3xy3z_D2x+ABZ*I_ESP_I3xy2z_D2x;
    Double I_ESP_I2x3yz_F2xz = I_ESP_K2x3y2z_D2x+ABZ*I_ESP_I2x3yz_D2x;
    Double I_ESP_I2x2y2z_F2xz = I_ESP_K2x2y3z_D2x+ABZ*I_ESP_I2x2y2z_D2x;
    Double I_ESP_I2xy3z_F2xz = I_ESP_K2xy4z_D2x+ABZ*I_ESP_I2xy3z_D2x;
    Double I_ESP_Ix4yz_F2xz = I_ESP_Kx4y2z_D2x+ABZ*I_ESP_Ix4yz_D2x;
    Double I_ESP_Ix3y2z_F2xz = I_ESP_Kx3y3z_D2x+ABZ*I_ESP_Ix3y2z_D2x;
    Double I_ESP_Ix2y3z_F2xz = I_ESP_Kx2y4z_D2x+ABZ*I_ESP_Ix2y3z_D2x;
    Double I_ESP_Ixy4z_F2xz = I_ESP_Kxy5z_D2x+ABZ*I_ESP_Ixy4z_D2x;
    Double I_ESP_I5yz_F2xz = I_ESP_K5y2z_D2x+ABZ*I_ESP_I5yz_D2x;
    Double I_ESP_I4y2z_F2xz = I_ESP_K4y3z_D2x+ABZ*I_ESP_I4y2z_D2x;
    Double I_ESP_I3y3z_F2xz = I_ESP_K3y4z_D2x+ABZ*I_ESP_I3y3z_D2x;
    Double I_ESP_I2y4z_F2xz = I_ESP_K2y5z_D2x+ABZ*I_ESP_I2y4z_D2x;
    Double I_ESP_Iy5z_F2xz = I_ESP_Ky6z_D2x+ABZ*I_ESP_Iy5z_D2x;
    Double I_ESP_I6x_F3y = I_ESP_K6xy_D2y+ABY*I_ESP_I6x_D2y;
    Double I_ESP_I5xy_F3y = I_ESP_K5x2y_D2y+ABY*I_ESP_I5xy_D2y;
    Double I_ESP_I5xz_F3y = I_ESP_K5xyz_D2y+ABY*I_ESP_I5xz_D2y;
    Double I_ESP_I4x2y_F3y = I_ESP_K4x3y_D2y+ABY*I_ESP_I4x2y_D2y;
    Double I_ESP_I4xyz_F3y = I_ESP_K4x2yz_D2y+ABY*I_ESP_I4xyz_D2y;
    Double I_ESP_I4x2z_F3y = I_ESP_K4xy2z_D2y+ABY*I_ESP_I4x2z_D2y;
    Double I_ESP_I3x3y_F3y = I_ESP_K3x4y_D2y+ABY*I_ESP_I3x3y_D2y;
    Double I_ESP_I3x2yz_F3y = I_ESP_K3x3yz_D2y+ABY*I_ESP_I3x2yz_D2y;
    Double I_ESP_I3xy2z_F3y = I_ESP_K3x2y2z_D2y+ABY*I_ESP_I3xy2z_D2y;
    Double I_ESP_I3x3z_F3y = I_ESP_K3xy3z_D2y+ABY*I_ESP_I3x3z_D2y;
    Double I_ESP_I2x4y_F3y = I_ESP_K2x5y_D2y+ABY*I_ESP_I2x4y_D2y;
    Double I_ESP_I2x3yz_F3y = I_ESP_K2x4yz_D2y+ABY*I_ESP_I2x3yz_D2y;
    Double I_ESP_I2x2y2z_F3y = I_ESP_K2x3y2z_D2y+ABY*I_ESP_I2x2y2z_D2y;
    Double I_ESP_I2xy3z_F3y = I_ESP_K2x2y3z_D2y+ABY*I_ESP_I2xy3z_D2y;
    Double I_ESP_I2x4z_F3y = I_ESP_K2xy4z_D2y+ABY*I_ESP_I2x4z_D2y;
    Double I_ESP_Ix5y_F3y = I_ESP_Kx6y_D2y+ABY*I_ESP_Ix5y_D2y;
    Double I_ESP_Ix4yz_F3y = I_ESP_Kx5yz_D2y+ABY*I_ESP_Ix4yz_D2y;
    Double I_ESP_Ix3y2z_F3y = I_ESP_Kx4y2z_D2y+ABY*I_ESP_Ix3y2z_D2y;
    Double I_ESP_Ix2y3z_F3y = I_ESP_Kx3y3z_D2y+ABY*I_ESP_Ix2y3z_D2y;
    Double I_ESP_Ixy4z_F3y = I_ESP_Kx2y4z_D2y+ABY*I_ESP_Ixy4z_D2y;
    Double I_ESP_Ix5z_F3y = I_ESP_Kxy5z_D2y+ABY*I_ESP_Ix5z_D2y;
    Double I_ESP_I6y_F3y = I_ESP_K7y_D2y+ABY*I_ESP_I6y_D2y;
    Double I_ESP_I5yz_F3y = I_ESP_K6yz_D2y+ABY*I_ESP_I5yz_D2y;
    Double I_ESP_I4y2z_F3y = I_ESP_K5y2z_D2y+ABY*I_ESP_I4y2z_D2y;
    Double I_ESP_I3y3z_F3y = I_ESP_K4y3z_D2y+ABY*I_ESP_I3y3z_D2y;
    Double I_ESP_I2y4z_F3y = I_ESP_K3y4z_D2y+ABY*I_ESP_I2y4z_D2y;
    Double I_ESP_Iy5z_F3y = I_ESP_K2y5z_D2y+ABY*I_ESP_Iy5z_D2y;
    Double I_ESP_I6z_F3y = I_ESP_Ky6z_D2y+ABY*I_ESP_I6z_D2y;
    Double I_ESP_I5xz_F2yz = I_ESP_K5x2z_D2y+ABZ*I_ESP_I5xz_D2y;
    Double I_ESP_I4xyz_F2yz = I_ESP_K4xy2z_D2y+ABZ*I_ESP_I4xyz_D2y;
    Double I_ESP_I4x2z_F2yz = I_ESP_K4x3z_D2y+ABZ*I_ESP_I4x2z_D2y;
    Double I_ESP_I3x2yz_F2yz = I_ESP_K3x2y2z_D2y+ABZ*I_ESP_I3x2yz_D2y;
    Double I_ESP_I3xy2z_F2yz = I_ESP_K3xy3z_D2y+ABZ*I_ESP_I3xy2z_D2y;
    Double I_ESP_I3x3z_F2yz = I_ESP_K3x4z_D2y+ABZ*I_ESP_I3x3z_D2y;
    Double I_ESP_I2x3yz_F2yz = I_ESP_K2x3y2z_D2y+ABZ*I_ESP_I2x3yz_D2y;
    Double I_ESP_I2x2y2z_F2yz = I_ESP_K2x2y3z_D2y+ABZ*I_ESP_I2x2y2z_D2y;
    Double I_ESP_I2xy3z_F2yz = I_ESP_K2xy4z_D2y+ABZ*I_ESP_I2xy3z_D2y;
    Double I_ESP_I2x4z_F2yz = I_ESP_K2x5z_D2y+ABZ*I_ESP_I2x4z_D2y;
    Double I_ESP_Ix4yz_F2yz = I_ESP_Kx4y2z_D2y+ABZ*I_ESP_Ix4yz_D2y;
    Double I_ESP_Ix3y2z_F2yz = I_ESP_Kx3y3z_D2y+ABZ*I_ESP_Ix3y2z_D2y;
    Double I_ESP_Ix2y3z_F2yz = I_ESP_Kx2y4z_D2y+ABZ*I_ESP_Ix2y3z_D2y;
    Double I_ESP_Ixy4z_F2yz = I_ESP_Kxy5z_D2y+ABZ*I_ESP_Ixy4z_D2y;
    Double I_ESP_Ix5z_F2yz = I_ESP_Kx6z_D2y+ABZ*I_ESP_Ix5z_D2y;
    Double I_ESP_I6x_F3z = I_ESP_K6xz_D2z+ABZ*I_ESP_I6x_D2z;
    Double I_ESP_I5xy_F3z = I_ESP_K5xyz_D2z+ABZ*I_ESP_I5xy_D2z;
    Double I_ESP_I5xz_F3z = I_ESP_K5x2z_D2z+ABZ*I_ESP_I5xz_D2z;
    Double I_ESP_I4x2y_F3z = I_ESP_K4x2yz_D2z+ABZ*I_ESP_I4x2y_D2z;
    Double I_ESP_I4xyz_F3z = I_ESP_K4xy2z_D2z+ABZ*I_ESP_I4xyz_D2z;
    Double I_ESP_I4x2z_F3z = I_ESP_K4x3z_D2z+ABZ*I_ESP_I4x2z_D2z;
    Double I_ESP_I3x3y_F3z = I_ESP_K3x3yz_D2z+ABZ*I_ESP_I3x3y_D2z;
    Double I_ESP_I3x2yz_F3z = I_ESP_K3x2y2z_D2z+ABZ*I_ESP_I3x2yz_D2z;
    Double I_ESP_I3xy2z_F3z = I_ESP_K3xy3z_D2z+ABZ*I_ESP_I3xy2z_D2z;
    Double I_ESP_I3x3z_F3z = I_ESP_K3x4z_D2z+ABZ*I_ESP_I3x3z_D2z;
    Double I_ESP_I2x4y_F3z = I_ESP_K2x4yz_D2z+ABZ*I_ESP_I2x4y_D2z;
    Double I_ESP_I2x3yz_F3z = I_ESP_K2x3y2z_D2z+ABZ*I_ESP_I2x3yz_D2z;
    Double I_ESP_I2x2y2z_F3z = I_ESP_K2x2y3z_D2z+ABZ*I_ESP_I2x2y2z_D2z;
    Double I_ESP_I2xy3z_F3z = I_ESP_K2xy4z_D2z+ABZ*I_ESP_I2xy3z_D2z;
    Double I_ESP_I2x4z_F3z = I_ESP_K2x5z_D2z+ABZ*I_ESP_I2x4z_D2z;
    Double I_ESP_Ix5y_F3z = I_ESP_Kx5yz_D2z+ABZ*I_ESP_Ix5y_D2z;
    Double I_ESP_Ix4yz_F3z = I_ESP_Kx4y2z_D2z+ABZ*I_ESP_Ix4yz_D2z;
    Double I_ESP_Ix3y2z_F3z = I_ESP_Kx3y3z_D2z+ABZ*I_ESP_Ix3y2z_D2z;
    Double I_ESP_Ix2y3z_F3z = I_ESP_Kx2y4z_D2z+ABZ*I_ESP_Ix2y3z_D2z;
    Double I_ESP_Ixy4z_F3z = I_ESP_Kxy5z_D2z+ABZ*I_ESP_Ixy4z_D2z;
    Double I_ESP_Ix5z_F3z = I_ESP_Kx6z_D2z+ABZ*I_ESP_Ix5z_D2z;
    Double I_ESP_I6y_F3z = I_ESP_K6yz_D2z+ABZ*I_ESP_I6y_D2z;
    Double I_ESP_I5yz_F3z = I_ESP_K5y2z_D2z+ABZ*I_ESP_I5yz_D2z;
    Double I_ESP_I4y2z_F3z = I_ESP_K4y3z_D2z+ABZ*I_ESP_I4y2z_D2z;
    Double I_ESP_I3y3z_F3z = I_ESP_K3y4z_D2z+ABZ*I_ESP_I3y3z_D2z;
    Double I_ESP_I2y4z_F3z = I_ESP_K2y5z_D2z+ABZ*I_ESP_I2y4z_D2z;
    Double I_ESP_Iy5z_F3z = I_ESP_Ky6z_D2z+ABZ*I_ESP_Iy5z_D2z;
    Double I_ESP_I6z_F3z = I_ESP_K7z_D2z+ABZ*I_ESP_I6z_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_H_G
     * expanding position: BRA2
     * code section is: HRR
     * totally 102 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F
     * RHS shell quartet name: SQ_ESP_H_F
     ************************************************************/
    Double I_ESP_H5x_G4x = I_ESP_I6x_F3x+ABX*I_ESP_H5x_F3x;
    Double I_ESP_H4xy_G4x = I_ESP_I5xy_F3x+ABX*I_ESP_H4xy_F3x;
    Double I_ESP_H4xz_G4x = I_ESP_I5xz_F3x+ABX*I_ESP_H4xz_F3x;
    Double I_ESP_H3x2y_G4x = I_ESP_I4x2y_F3x+ABX*I_ESP_H3x2y_F3x;
    Double I_ESP_H3xyz_G4x = I_ESP_I4xyz_F3x+ABX*I_ESP_H3xyz_F3x;
    Double I_ESP_H3x2z_G4x = I_ESP_I4x2z_F3x+ABX*I_ESP_H3x2z_F3x;
    Double I_ESP_H2x3y_G4x = I_ESP_I3x3y_F3x+ABX*I_ESP_H2x3y_F3x;
    Double I_ESP_H2x2yz_G4x = I_ESP_I3x2yz_F3x+ABX*I_ESP_H2x2yz_F3x;
    Double I_ESP_H2xy2z_G4x = I_ESP_I3xy2z_F3x+ABX*I_ESP_H2xy2z_F3x;
    Double I_ESP_H2x3z_G4x = I_ESP_I3x3z_F3x+ABX*I_ESP_H2x3z_F3x;
    Double I_ESP_Hx4y_G4x = I_ESP_I2x4y_F3x+ABX*I_ESP_Hx4y_F3x;
    Double I_ESP_Hx3yz_G4x = I_ESP_I2x3yz_F3x+ABX*I_ESP_Hx3yz_F3x;
    Double I_ESP_Hx2y2z_G4x = I_ESP_I2x2y2z_F3x+ABX*I_ESP_Hx2y2z_F3x;
    Double I_ESP_Hxy3z_G4x = I_ESP_I2xy3z_F3x+ABX*I_ESP_Hxy3z_F3x;
    Double I_ESP_Hx4z_G4x = I_ESP_I2x4z_F3x+ABX*I_ESP_Hx4z_F3x;
    Double I_ESP_H5y_G4x = I_ESP_Ix5y_F3x+ABX*I_ESP_H5y_F3x;
    Double I_ESP_H4yz_G4x = I_ESP_Ix4yz_F3x+ABX*I_ESP_H4yz_F3x;
    Double I_ESP_H3y2z_G4x = I_ESP_Ix3y2z_F3x+ABX*I_ESP_H3y2z_F3x;
    Double I_ESP_H2y3z_G4x = I_ESP_Ix2y3z_F3x+ABX*I_ESP_H2y3z_F3x;
    Double I_ESP_Hy4z_G4x = I_ESP_Ixy4z_F3x+ABX*I_ESP_Hy4z_F3x;
    Double I_ESP_H5z_G4x = I_ESP_Ix5z_F3x+ABX*I_ESP_H5z_F3x;
    Double I_ESP_H4xy_G3xy = I_ESP_I4x2y_F3x+ABY*I_ESP_H4xy_F3x;
    Double I_ESP_H4xz_G3xy = I_ESP_I4xyz_F3x+ABY*I_ESP_H4xz_F3x;
    Double I_ESP_H3x2y_G3xy = I_ESP_I3x3y_F3x+ABY*I_ESP_H3x2y_F3x;
    Double I_ESP_H3xyz_G3xy = I_ESP_I3x2yz_F3x+ABY*I_ESP_H3xyz_F3x;
    Double I_ESP_H3x2z_G3xy = I_ESP_I3xy2z_F3x+ABY*I_ESP_H3x2z_F3x;
    Double I_ESP_H2x3y_G3xy = I_ESP_I2x4y_F3x+ABY*I_ESP_H2x3y_F3x;
    Double I_ESP_H2x2yz_G3xy = I_ESP_I2x3yz_F3x+ABY*I_ESP_H2x2yz_F3x;
    Double I_ESP_H2xy2z_G3xy = I_ESP_I2x2y2z_F3x+ABY*I_ESP_H2xy2z_F3x;
    Double I_ESP_H2x3z_G3xy = I_ESP_I2xy3z_F3x+ABY*I_ESP_H2x3z_F3x;
    Double I_ESP_Hx4y_G3xy = I_ESP_Ix5y_F3x+ABY*I_ESP_Hx4y_F3x;
    Double I_ESP_Hx3yz_G3xy = I_ESP_Ix4yz_F3x+ABY*I_ESP_Hx3yz_F3x;
    Double I_ESP_Hx2y2z_G3xy = I_ESP_Ix3y2z_F3x+ABY*I_ESP_Hx2y2z_F3x;
    Double I_ESP_Hxy3z_G3xy = I_ESP_Ix2y3z_F3x+ABY*I_ESP_Hxy3z_F3x;
    Double I_ESP_Hx4z_G3xy = I_ESP_Ixy4z_F3x+ABY*I_ESP_Hx4z_F3x;
    Double I_ESP_H5y_G3xy = I_ESP_I6y_F3x+ABY*I_ESP_H5y_F3x;
    Double I_ESP_H4yz_G3xy = I_ESP_I5yz_F3x+ABY*I_ESP_H4yz_F3x;
    Double I_ESP_H3y2z_G3xy = I_ESP_I4y2z_F3x+ABY*I_ESP_H3y2z_F3x;
    Double I_ESP_H2y3z_G3xy = I_ESP_I3y3z_F3x+ABY*I_ESP_H2y3z_F3x;
    Double I_ESP_Hy4z_G3xy = I_ESP_I2y4z_F3x+ABY*I_ESP_Hy4z_F3x;
    Double I_ESP_H5z_G3xy = I_ESP_Iy5z_F3x+ABY*I_ESP_H5z_F3x;
    Double I_ESP_H4xz_G3xz = I_ESP_I4x2z_F3x+ABZ*I_ESP_H4xz_F3x;
    Double I_ESP_H3xyz_G3xz = I_ESP_I3xy2z_F3x+ABZ*I_ESP_H3xyz_F3x;
    Double I_ESP_H3x2z_G3xz = I_ESP_I3x3z_F3x+ABZ*I_ESP_H3x2z_F3x;
    Double I_ESP_H2x2yz_G3xz = I_ESP_I2x2y2z_F3x+ABZ*I_ESP_H2x2yz_F3x;
    Double I_ESP_H2xy2z_G3xz = I_ESP_I2xy3z_F3x+ABZ*I_ESP_H2xy2z_F3x;
    Double I_ESP_H2x3z_G3xz = I_ESP_I2x4z_F3x+ABZ*I_ESP_H2x3z_F3x;
    Double I_ESP_Hx3yz_G3xz = I_ESP_Ix3y2z_F3x+ABZ*I_ESP_Hx3yz_F3x;
    Double I_ESP_Hx2y2z_G3xz = I_ESP_Ix2y3z_F3x+ABZ*I_ESP_Hx2y2z_F3x;
    Double I_ESP_Hxy3z_G3xz = I_ESP_Ixy4z_F3x+ABZ*I_ESP_Hxy3z_F3x;
    Double I_ESP_Hx4z_G3xz = I_ESP_Ix5z_F3x+ABZ*I_ESP_Hx4z_F3x;
    Double I_ESP_H4yz_G3xz = I_ESP_I4y2z_F3x+ABZ*I_ESP_H4yz_F3x;
    Double I_ESP_H3y2z_G3xz = I_ESP_I3y3z_F3x+ABZ*I_ESP_H3y2z_F3x;
    Double I_ESP_H2y3z_G3xz = I_ESP_I2y4z_F3x+ABZ*I_ESP_H2y3z_F3x;
    Double I_ESP_Hy4z_G3xz = I_ESP_Iy5z_F3x+ABZ*I_ESP_Hy4z_F3x;
    Double I_ESP_H5z_G3xz = I_ESP_I6z_F3x+ABZ*I_ESP_H5z_F3x;
    Double I_ESP_H4xz_G2x2y = I_ESP_I4xyz_F2xy+ABY*I_ESP_H4xz_F2xy;
    Double I_ESP_H3xyz_G2x2y = I_ESP_I3x2yz_F2xy+ABY*I_ESP_H3xyz_F2xy;
    Double I_ESP_H3x2z_G2x2y = I_ESP_I3xy2z_F2xy+ABY*I_ESP_H3x2z_F2xy;
    Double I_ESP_H2x2yz_G2x2y = I_ESP_I2x3yz_F2xy+ABY*I_ESP_H2x2yz_F2xy;
    Double I_ESP_H2xy2z_G2x2y = I_ESP_I2x2y2z_F2xy+ABY*I_ESP_H2xy2z_F2xy;
    Double I_ESP_H2x3z_G2x2y = I_ESP_I2xy3z_F2xy+ABY*I_ESP_H2x3z_F2xy;
    Double I_ESP_Hx3yz_G2x2y = I_ESP_Ix4yz_F2xy+ABY*I_ESP_Hx3yz_F2xy;
    Double I_ESP_Hx2y2z_G2x2y = I_ESP_Ix3y2z_F2xy+ABY*I_ESP_Hx2y2z_F2xy;
    Double I_ESP_Hxy3z_G2x2y = I_ESP_Ix2y3z_F2xy+ABY*I_ESP_Hxy3z_F2xy;
    Double I_ESP_Hx4z_G2x2y = I_ESP_Ixy4z_F2xy+ABY*I_ESP_Hx4z_F2xy;
    Double I_ESP_H4yz_G2x2y = I_ESP_I5yz_F2xy+ABY*I_ESP_H4yz_F2xy;
    Double I_ESP_H3y2z_G2x2y = I_ESP_I4y2z_F2xy+ABY*I_ESP_H3y2z_F2xy;
    Double I_ESP_H2y3z_G2x2y = I_ESP_I3y3z_F2xy+ABY*I_ESP_H2y3z_F2xy;
    Double I_ESP_Hy4z_G2x2y = I_ESP_I2y4z_F2xy+ABY*I_ESP_Hy4z_F2xy;
    Double I_ESP_H5z_G2x2y = I_ESP_Iy5z_F2xy+ABY*I_ESP_H5z_F2xy;
    Double I_ESP_H4xy_G2x2z = I_ESP_I4xyz_F2xz+ABZ*I_ESP_H4xy_F2xz;
    Double I_ESP_H3x2y_G2x2z = I_ESP_I3x2yz_F2xz+ABZ*I_ESP_H3x2y_F2xz;
    Double I_ESP_H3xyz_G2x2z = I_ESP_I3xy2z_F2xz+ABZ*I_ESP_H3xyz_F2xz;
    Double I_ESP_H2x3y_G2x2z = I_ESP_I2x3yz_F2xz+ABZ*I_ESP_H2x3y_F2xz;
    Double I_ESP_H2x2yz_G2x2z = I_ESP_I2x2y2z_F2xz+ABZ*I_ESP_H2x2yz_F2xz;
    Double I_ESP_H2xy2z_G2x2z = I_ESP_I2xy3z_F2xz+ABZ*I_ESP_H2xy2z_F2xz;
    Double I_ESP_Hx4y_G2x2z = I_ESP_Ix4yz_F2xz+ABZ*I_ESP_Hx4y_F2xz;
    Double I_ESP_Hx3yz_G2x2z = I_ESP_Ix3y2z_F2xz+ABZ*I_ESP_Hx3yz_F2xz;
    Double I_ESP_Hx2y2z_G2x2z = I_ESP_Ix2y3z_F2xz+ABZ*I_ESP_Hx2y2z_F2xz;
    Double I_ESP_Hxy3z_G2x2z = I_ESP_Ixy4z_F2xz+ABZ*I_ESP_Hxy3z_F2xz;
    Double I_ESP_H5y_G2x2z = I_ESP_I5yz_F2xz+ABZ*I_ESP_H5y_F2xz;
    Double I_ESP_H4yz_G2x2z = I_ESP_I4y2z_F2xz+ABZ*I_ESP_H4yz_F2xz;
    Double I_ESP_H3y2z_G2x2z = I_ESP_I3y3z_F2xz+ABZ*I_ESP_H3y2z_F2xz;
    Double I_ESP_H2y3z_G2x2z = I_ESP_I2y4z_F2xz+ABZ*I_ESP_H2y3z_F2xz;
    Double I_ESP_Hy4z_G2x2z = I_ESP_Iy5z_F2xz+ABZ*I_ESP_Hy4z_F2xz;
    Double I_ESP_H5x_Gx3y = I_ESP_I6x_F3y+ABX*I_ESP_H5x_F3y;
    Double I_ESP_H4xy_Gx3y = I_ESP_I5xy_F3y+ABX*I_ESP_H4xy_F3y;
    Double I_ESP_H4xz_Gx3y = I_ESP_I5xz_F3y+ABX*I_ESP_H4xz_F3y;
    Double I_ESP_H3x2y_Gx3y = I_ESP_I4x2y_F3y+ABX*I_ESP_H3x2y_F3y;
    Double I_ESP_H3xyz_Gx3y = I_ESP_I4xyz_F3y+ABX*I_ESP_H3xyz_F3y;
    Double I_ESP_H3x2z_Gx3y = I_ESP_I4x2z_F3y+ABX*I_ESP_H3x2z_F3y;
    Double I_ESP_H2x3y_Gx3y = I_ESP_I3x3y_F3y+ABX*I_ESP_H2x3y_F3y;
    Double I_ESP_H2x2yz_Gx3y = I_ESP_I3x2yz_F3y+ABX*I_ESP_H2x2yz_F3y;
    Double I_ESP_H2xy2z_Gx3y = I_ESP_I3xy2z_F3y+ABX*I_ESP_H2xy2z_F3y;
    Double I_ESP_H2x3z_Gx3y = I_ESP_I3x3z_F3y+ABX*I_ESP_H2x3z_F3y;
    Double I_ESP_Hx4y_Gx3y = I_ESP_I2x4y_F3y+ABX*I_ESP_Hx4y_F3y;
    Double I_ESP_Hx3yz_Gx3y = I_ESP_I2x3yz_F3y+ABX*I_ESP_Hx3yz_F3y;
    Double I_ESP_Hx2y2z_Gx3y = I_ESP_I2x2y2z_F3y+ABX*I_ESP_Hx2y2z_F3y;
    Double I_ESP_Hxy3z_Gx3y = I_ESP_I2xy3z_F3y+ABX*I_ESP_Hxy3z_F3y;
    Double I_ESP_Hx4z_Gx3y = I_ESP_I2x4z_F3y+ABX*I_ESP_Hx4z_F3y;
    Double I_ESP_H4yz_Gx3y = I_ESP_Ix4yz_F3y+ABX*I_ESP_H4yz_F3y;
    Double I_ESP_H3y2z_Gx3y = I_ESP_Ix3y2z_F3y+ABX*I_ESP_H3y2z_F3y;
    Double I_ESP_H2y3z_Gx3y = I_ESP_Ix2y3z_F3y+ABX*I_ESP_H2y3z_F3y;
    Double I_ESP_Hy4z_Gx3y = I_ESP_Ixy4z_F3y+ABX*I_ESP_Hy4z_F3y;
    Double I_ESP_H5z_Gx3y = I_ESP_Ix5z_F3y+ABX*I_ESP_H5z_F3y;
    Double I_ESP_H5x_Gx3z = I_ESP_I6x_F3z+ABX*I_ESP_H5x_F3z;
    Double I_ESP_H4xy_Gx3z = I_ESP_I5xy_F3z+ABX*I_ESP_H4xy_F3z;
    Double I_ESP_H4xz_Gx3z = I_ESP_I5xz_F3z+ABX*I_ESP_H4xz_F3z;
    Double I_ESP_H3x2y_Gx3z = I_ESP_I4x2y_F3z+ABX*I_ESP_H3x2y_F3z;
    Double I_ESP_H3xyz_Gx3z = I_ESP_I4xyz_F3z+ABX*I_ESP_H3xyz_F3z;
    Double I_ESP_H3x2z_Gx3z = I_ESP_I4x2z_F3z+ABX*I_ESP_H3x2z_F3z;
    Double I_ESP_H2x3y_Gx3z = I_ESP_I3x3y_F3z+ABX*I_ESP_H2x3y_F3z;
    Double I_ESP_H2x2yz_Gx3z = I_ESP_I3x2yz_F3z+ABX*I_ESP_H2x2yz_F3z;
    Double I_ESP_H2xy2z_Gx3z = I_ESP_I3xy2z_F3z+ABX*I_ESP_H2xy2z_F3z;
    Double I_ESP_H2x3z_Gx3z = I_ESP_I3x3z_F3z+ABX*I_ESP_H2x3z_F3z;
    Double I_ESP_Hx4y_Gx3z = I_ESP_I2x4y_F3z+ABX*I_ESP_Hx4y_F3z;
    Double I_ESP_Hx3yz_Gx3z = I_ESP_I2x3yz_F3z+ABX*I_ESP_Hx3yz_F3z;
    Double I_ESP_Hx2y2z_Gx3z = I_ESP_I2x2y2z_F3z+ABX*I_ESP_Hx2y2z_F3z;
    Double I_ESP_Hxy3z_Gx3z = I_ESP_I2xy3z_F3z+ABX*I_ESP_Hxy3z_F3z;
    Double I_ESP_Hx4z_Gx3z = I_ESP_I2x4z_F3z+ABX*I_ESP_Hx4z_F3z;
    Double I_ESP_H5y_Gx3z = I_ESP_Ix5y_F3z+ABX*I_ESP_H5y_F3z;
    Double I_ESP_H4yz_Gx3z = I_ESP_Ix4yz_F3z+ABX*I_ESP_H4yz_F3z;
    Double I_ESP_H3y2z_Gx3z = I_ESP_Ix3y2z_F3z+ABX*I_ESP_H3y2z_F3z;
    Double I_ESP_H2y3z_Gx3z = I_ESP_Ix2y3z_F3z+ABX*I_ESP_H2y3z_F3z;
    Double I_ESP_Hy4z_Gx3z = I_ESP_Ixy4z_F3z+ABX*I_ESP_Hy4z_F3z;
    Double I_ESP_H5x_G4y = I_ESP_I5xy_F3y+ABY*I_ESP_H5x_F3y;
    Double I_ESP_H4xy_G4y = I_ESP_I4x2y_F3y+ABY*I_ESP_H4xy_F3y;
    Double I_ESP_H4xz_G4y = I_ESP_I4xyz_F3y+ABY*I_ESP_H4xz_F3y;
    Double I_ESP_H3x2y_G4y = I_ESP_I3x3y_F3y+ABY*I_ESP_H3x2y_F3y;
    Double I_ESP_H3xyz_G4y = I_ESP_I3x2yz_F3y+ABY*I_ESP_H3xyz_F3y;
    Double I_ESP_H3x2z_G4y = I_ESP_I3xy2z_F3y+ABY*I_ESP_H3x2z_F3y;
    Double I_ESP_H2x3y_G4y = I_ESP_I2x4y_F3y+ABY*I_ESP_H2x3y_F3y;
    Double I_ESP_H2x2yz_G4y = I_ESP_I2x3yz_F3y+ABY*I_ESP_H2x2yz_F3y;
    Double I_ESP_H2xy2z_G4y = I_ESP_I2x2y2z_F3y+ABY*I_ESP_H2xy2z_F3y;
    Double I_ESP_H2x3z_G4y = I_ESP_I2xy3z_F3y+ABY*I_ESP_H2x3z_F3y;
    Double I_ESP_Hx4y_G4y = I_ESP_Ix5y_F3y+ABY*I_ESP_Hx4y_F3y;
    Double I_ESP_Hx3yz_G4y = I_ESP_Ix4yz_F3y+ABY*I_ESP_Hx3yz_F3y;
    Double I_ESP_Hx2y2z_G4y = I_ESP_Ix3y2z_F3y+ABY*I_ESP_Hx2y2z_F3y;
    Double I_ESP_Hxy3z_G4y = I_ESP_Ix2y3z_F3y+ABY*I_ESP_Hxy3z_F3y;
    Double I_ESP_Hx4z_G4y = I_ESP_Ixy4z_F3y+ABY*I_ESP_Hx4z_F3y;
    Double I_ESP_H5y_G4y = I_ESP_I6y_F3y+ABY*I_ESP_H5y_F3y;
    Double I_ESP_H4yz_G4y = I_ESP_I5yz_F3y+ABY*I_ESP_H4yz_F3y;
    Double I_ESP_H3y2z_G4y = I_ESP_I4y2z_F3y+ABY*I_ESP_H3y2z_F3y;
    Double I_ESP_H2y3z_G4y = I_ESP_I3y3z_F3y+ABY*I_ESP_H2y3z_F3y;
    Double I_ESP_Hy4z_G4y = I_ESP_I2y4z_F3y+ABY*I_ESP_Hy4z_F3y;
    Double I_ESP_H5z_G4y = I_ESP_Iy5z_F3y+ABY*I_ESP_H5z_F3y;
    Double I_ESP_H4xz_G3yz = I_ESP_I4x2z_F3y+ABZ*I_ESP_H4xz_F3y;
    Double I_ESP_H3xyz_G3yz = I_ESP_I3xy2z_F3y+ABZ*I_ESP_H3xyz_F3y;
    Double I_ESP_H3x2z_G3yz = I_ESP_I3x3z_F3y+ABZ*I_ESP_H3x2z_F3y;
    Double I_ESP_H2x2yz_G3yz = I_ESP_I2x2y2z_F3y+ABZ*I_ESP_H2x2yz_F3y;
    Double I_ESP_H2xy2z_G3yz = I_ESP_I2xy3z_F3y+ABZ*I_ESP_H2xy2z_F3y;
    Double I_ESP_H2x3z_G3yz = I_ESP_I2x4z_F3y+ABZ*I_ESP_H2x3z_F3y;
    Double I_ESP_Hx3yz_G3yz = I_ESP_Ix3y2z_F3y+ABZ*I_ESP_Hx3yz_F3y;
    Double I_ESP_Hx2y2z_G3yz = I_ESP_Ix2y3z_F3y+ABZ*I_ESP_Hx2y2z_F3y;
    Double I_ESP_Hxy3z_G3yz = I_ESP_Ixy4z_F3y+ABZ*I_ESP_Hxy3z_F3y;
    Double I_ESP_Hx4z_G3yz = I_ESP_Ix5z_F3y+ABZ*I_ESP_Hx4z_F3y;
    Double I_ESP_H4yz_G3yz = I_ESP_I4y2z_F3y+ABZ*I_ESP_H4yz_F3y;
    Double I_ESP_H3y2z_G3yz = I_ESP_I3y3z_F3y+ABZ*I_ESP_H3y2z_F3y;
    Double I_ESP_H2y3z_G3yz = I_ESP_I2y4z_F3y+ABZ*I_ESP_H2y3z_F3y;
    Double I_ESP_Hy4z_G3yz = I_ESP_Iy5z_F3y+ABZ*I_ESP_Hy4z_F3y;
    Double I_ESP_H5z_G3yz = I_ESP_I6z_F3y+ABZ*I_ESP_H5z_F3y;
    Double I_ESP_H5x_G2y2z = I_ESP_I5xz_F2yz+ABZ*I_ESP_H5x_F2yz;
    Double I_ESP_H4xy_G2y2z = I_ESP_I4xyz_F2yz+ABZ*I_ESP_H4xy_F2yz;
    Double I_ESP_H4xz_G2y2z = I_ESP_I4x2z_F2yz+ABZ*I_ESP_H4xz_F2yz;
    Double I_ESP_H3x2y_G2y2z = I_ESP_I3x2yz_F2yz+ABZ*I_ESP_H3x2y_F2yz;
    Double I_ESP_H3xyz_G2y2z = I_ESP_I3xy2z_F2yz+ABZ*I_ESP_H3xyz_F2yz;
    Double I_ESP_H3x2z_G2y2z = I_ESP_I3x3z_F2yz+ABZ*I_ESP_H3x2z_F2yz;
    Double I_ESP_H2x3y_G2y2z = I_ESP_I2x3yz_F2yz+ABZ*I_ESP_H2x3y_F2yz;
    Double I_ESP_H2x2yz_G2y2z = I_ESP_I2x2y2z_F2yz+ABZ*I_ESP_H2x2yz_F2yz;
    Double I_ESP_H2xy2z_G2y2z = I_ESP_I2xy3z_F2yz+ABZ*I_ESP_H2xy2z_F2yz;
    Double I_ESP_H2x3z_G2y2z = I_ESP_I2x4z_F2yz+ABZ*I_ESP_H2x3z_F2yz;
    Double I_ESP_Hx4y_G2y2z = I_ESP_Ix4yz_F2yz+ABZ*I_ESP_Hx4y_F2yz;
    Double I_ESP_Hx3yz_G2y2z = I_ESP_Ix3y2z_F2yz+ABZ*I_ESP_Hx3yz_F2yz;
    Double I_ESP_Hx2y2z_G2y2z = I_ESP_Ix2y3z_F2yz+ABZ*I_ESP_Hx2y2z_F2yz;
    Double I_ESP_Hxy3z_G2y2z = I_ESP_Ixy4z_F2yz+ABZ*I_ESP_Hxy3z_F2yz;
    Double I_ESP_Hx4z_G2y2z = I_ESP_Ix5z_F2yz+ABZ*I_ESP_Hx4z_F2yz;
    Double I_ESP_H4xy_Gy3z = I_ESP_I4x2y_F3z+ABY*I_ESP_H4xy_F3z;
    Double I_ESP_H3x2y_Gy3z = I_ESP_I3x3y_F3z+ABY*I_ESP_H3x2y_F3z;
    Double I_ESP_H3xyz_Gy3z = I_ESP_I3x2yz_F3z+ABY*I_ESP_H3xyz_F3z;
    Double I_ESP_H2x3y_Gy3z = I_ESP_I2x4y_F3z+ABY*I_ESP_H2x3y_F3z;
    Double I_ESP_H2x2yz_Gy3z = I_ESP_I2x3yz_F3z+ABY*I_ESP_H2x2yz_F3z;
    Double I_ESP_H2xy2z_Gy3z = I_ESP_I2x2y2z_F3z+ABY*I_ESP_H2xy2z_F3z;
    Double I_ESP_Hx4y_Gy3z = I_ESP_Ix5y_F3z+ABY*I_ESP_Hx4y_F3z;
    Double I_ESP_Hx3yz_Gy3z = I_ESP_Ix4yz_F3z+ABY*I_ESP_Hx3yz_F3z;
    Double I_ESP_Hx2y2z_Gy3z = I_ESP_Ix3y2z_F3z+ABY*I_ESP_Hx2y2z_F3z;
    Double I_ESP_Hxy3z_Gy3z = I_ESP_Ix2y3z_F3z+ABY*I_ESP_Hxy3z_F3z;
    Double I_ESP_H5y_Gy3z = I_ESP_I6y_F3z+ABY*I_ESP_H5y_F3z;
    Double I_ESP_H4yz_Gy3z = I_ESP_I5yz_F3z+ABY*I_ESP_H4yz_F3z;
    Double I_ESP_H3y2z_Gy3z = I_ESP_I4y2z_F3z+ABY*I_ESP_H3y2z_F3z;
    Double I_ESP_H2y3z_Gy3z = I_ESP_I3y3z_F3z+ABY*I_ESP_H2y3z_F3z;
    Double I_ESP_Hy4z_Gy3z = I_ESP_I2y4z_F3z+ABY*I_ESP_Hy4z_F3z;
    Double I_ESP_H5x_G4z = I_ESP_I5xz_F3z+ABZ*I_ESP_H5x_F3z;
    Double I_ESP_H4xy_G4z = I_ESP_I4xyz_F3z+ABZ*I_ESP_H4xy_F3z;
    Double I_ESP_H4xz_G4z = I_ESP_I4x2z_F3z+ABZ*I_ESP_H4xz_F3z;
    Double I_ESP_H3x2y_G4z = I_ESP_I3x2yz_F3z+ABZ*I_ESP_H3x2y_F3z;
    Double I_ESP_H3xyz_G4z = I_ESP_I3xy2z_F3z+ABZ*I_ESP_H3xyz_F3z;
    Double I_ESP_H3x2z_G4z = I_ESP_I3x3z_F3z+ABZ*I_ESP_H3x2z_F3z;
    Double I_ESP_H2x3y_G4z = I_ESP_I2x3yz_F3z+ABZ*I_ESP_H2x3y_F3z;
    Double I_ESP_H2x2yz_G4z = I_ESP_I2x2y2z_F3z+ABZ*I_ESP_H2x2yz_F3z;
    Double I_ESP_H2xy2z_G4z = I_ESP_I2xy3z_F3z+ABZ*I_ESP_H2xy2z_F3z;
    Double I_ESP_H2x3z_G4z = I_ESP_I2x4z_F3z+ABZ*I_ESP_H2x3z_F3z;
    Double I_ESP_Hx4y_G4z = I_ESP_Ix4yz_F3z+ABZ*I_ESP_Hx4y_F3z;
    Double I_ESP_Hx3yz_G4z = I_ESP_Ix3y2z_F3z+ABZ*I_ESP_Hx3yz_F3z;
    Double I_ESP_Hx2y2z_G4z = I_ESP_Ix2y3z_F3z+ABZ*I_ESP_Hx2y2z_F3z;
    Double I_ESP_Hxy3z_G4z = I_ESP_Ixy4z_F3z+ABZ*I_ESP_Hxy3z_F3z;
    Double I_ESP_Hx4z_G4z = I_ESP_Ix5z_F3z+ABZ*I_ESP_Hx4z_F3z;
    Double I_ESP_H5y_G4z = I_ESP_I5yz_F3z+ABZ*I_ESP_H5y_F3z;
    Double I_ESP_H4yz_G4z = I_ESP_I4y2z_F3z+ABZ*I_ESP_H4yz_F3z;
    Double I_ESP_H3y2z_G4z = I_ESP_I3y3z_F3z+ABZ*I_ESP_H3y2z_F3z;
    Double I_ESP_H2y3z_G4z = I_ESP_I2y4z_F3z+ABZ*I_ESP_H2y3z_F3z;
    Double I_ESP_Hy4z_G4z = I_ESP_Iy5z_F3z+ABZ*I_ESP_Hy4z_F3z;
    Double I_ESP_H5z_G4z = I_ESP_I6z_F3z+ABZ*I_ESP_H5z_F3z;

    /************************************************************
     * shell quartet name: SQ_ESP_G_H
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_G
     * RHS shell quartet name: SQ_ESP_G_G
     ************************************************************/
    Double I_ESP_G4x_H5x = I_ESP_H5x_G4x+ABX*I_ESP_G4x_G4x;
    Double I_ESP_G3xy_H5x = I_ESP_H4xy_G4x+ABX*I_ESP_G3xy_G4x;
    Double I_ESP_G3xz_H5x = I_ESP_H4xz_G4x+ABX*I_ESP_G3xz_G4x;
    Double I_ESP_G2x2y_H5x = I_ESP_H3x2y_G4x+ABX*I_ESP_G2x2y_G4x;
    Double I_ESP_G2xyz_H5x = I_ESP_H3xyz_G4x+ABX*I_ESP_G2xyz_G4x;
    Double I_ESP_G2x2z_H5x = I_ESP_H3x2z_G4x+ABX*I_ESP_G2x2z_G4x;
    Double I_ESP_Gx3y_H5x = I_ESP_H2x3y_G4x+ABX*I_ESP_Gx3y_G4x;
    Double I_ESP_Gx2yz_H5x = I_ESP_H2x2yz_G4x+ABX*I_ESP_Gx2yz_G4x;
    Double I_ESP_Gxy2z_H5x = I_ESP_H2xy2z_G4x+ABX*I_ESP_Gxy2z_G4x;
    Double I_ESP_Gx3z_H5x = I_ESP_H2x3z_G4x+ABX*I_ESP_Gx3z_G4x;
    Double I_ESP_G4y_H5x = I_ESP_Hx4y_G4x+ABX*I_ESP_G4y_G4x;
    Double I_ESP_G3yz_H5x = I_ESP_Hx3yz_G4x+ABX*I_ESP_G3yz_G4x;
    Double I_ESP_G2y2z_H5x = I_ESP_Hx2y2z_G4x+ABX*I_ESP_G2y2z_G4x;
    Double I_ESP_Gy3z_H5x = I_ESP_Hxy3z_G4x+ABX*I_ESP_Gy3z_G4x;
    Double I_ESP_G4z_H5x = I_ESP_Hx4z_G4x+ABX*I_ESP_G4z_G4x;
    Double I_ESP_G4x_H4xy = I_ESP_H4xy_G4x+ABY*I_ESP_G4x_G4x;
    Double I_ESP_G3xy_H4xy = I_ESP_H3x2y_G4x+ABY*I_ESP_G3xy_G4x;
    Double I_ESP_G3xz_H4xy = I_ESP_H3xyz_G4x+ABY*I_ESP_G3xz_G4x;
    Double I_ESP_G2x2y_H4xy = I_ESP_H2x3y_G4x+ABY*I_ESP_G2x2y_G4x;
    Double I_ESP_G2xyz_H4xy = I_ESP_H2x2yz_G4x+ABY*I_ESP_G2xyz_G4x;
    Double I_ESP_G2x2z_H4xy = I_ESP_H2xy2z_G4x+ABY*I_ESP_G2x2z_G4x;
    Double I_ESP_Gx3y_H4xy = I_ESP_Hx4y_G4x+ABY*I_ESP_Gx3y_G4x;
    Double I_ESP_Gx2yz_H4xy = I_ESP_Hx3yz_G4x+ABY*I_ESP_Gx2yz_G4x;
    Double I_ESP_Gxy2z_H4xy = I_ESP_Hx2y2z_G4x+ABY*I_ESP_Gxy2z_G4x;
    Double I_ESP_Gx3z_H4xy = I_ESP_Hxy3z_G4x+ABY*I_ESP_Gx3z_G4x;
    Double I_ESP_G4y_H4xy = I_ESP_H5y_G4x+ABY*I_ESP_G4y_G4x;
    Double I_ESP_G3yz_H4xy = I_ESP_H4yz_G4x+ABY*I_ESP_G3yz_G4x;
    Double I_ESP_G2y2z_H4xy = I_ESP_H3y2z_G4x+ABY*I_ESP_G2y2z_G4x;
    Double I_ESP_Gy3z_H4xy = I_ESP_H2y3z_G4x+ABY*I_ESP_Gy3z_G4x;
    Double I_ESP_G4z_H4xy = I_ESP_Hy4z_G4x+ABY*I_ESP_G4z_G4x;
    Double I_ESP_G4x_H4xz = I_ESP_H4xz_G4x+ABZ*I_ESP_G4x_G4x;
    Double I_ESP_G3xy_H4xz = I_ESP_H3xyz_G4x+ABZ*I_ESP_G3xy_G4x;
    Double I_ESP_G3xz_H4xz = I_ESP_H3x2z_G4x+ABZ*I_ESP_G3xz_G4x;
    Double I_ESP_G2x2y_H4xz = I_ESP_H2x2yz_G4x+ABZ*I_ESP_G2x2y_G4x;
    Double I_ESP_G2xyz_H4xz = I_ESP_H2xy2z_G4x+ABZ*I_ESP_G2xyz_G4x;
    Double I_ESP_G2x2z_H4xz = I_ESP_H2x3z_G4x+ABZ*I_ESP_G2x2z_G4x;
    Double I_ESP_Gx3y_H4xz = I_ESP_Hx3yz_G4x+ABZ*I_ESP_Gx3y_G4x;
    Double I_ESP_Gx2yz_H4xz = I_ESP_Hx2y2z_G4x+ABZ*I_ESP_Gx2yz_G4x;
    Double I_ESP_Gxy2z_H4xz = I_ESP_Hxy3z_G4x+ABZ*I_ESP_Gxy2z_G4x;
    Double I_ESP_Gx3z_H4xz = I_ESP_Hx4z_G4x+ABZ*I_ESP_Gx3z_G4x;
    Double I_ESP_G4y_H4xz = I_ESP_H4yz_G4x+ABZ*I_ESP_G4y_G4x;
    Double I_ESP_G3yz_H4xz = I_ESP_H3y2z_G4x+ABZ*I_ESP_G3yz_G4x;
    Double I_ESP_G2y2z_H4xz = I_ESP_H2y3z_G4x+ABZ*I_ESP_G2y2z_G4x;
    Double I_ESP_Gy3z_H4xz = I_ESP_Hy4z_G4x+ABZ*I_ESP_Gy3z_G4x;
    Double I_ESP_G4z_H4xz = I_ESP_H5z_G4x+ABZ*I_ESP_G4z_G4x;
    Double I_ESP_G4x_H3x2y = I_ESP_H4xy_G3xy+ABY*I_ESP_G4x_G3xy;
    Double I_ESP_G3xy_H3x2y = I_ESP_H3x2y_G3xy+ABY*I_ESP_G3xy_G3xy;
    Double I_ESP_G3xz_H3x2y = I_ESP_H3xyz_G3xy+ABY*I_ESP_G3xz_G3xy;
    Double I_ESP_G2x2y_H3x2y = I_ESP_H2x3y_G3xy+ABY*I_ESP_G2x2y_G3xy;
    Double I_ESP_G2xyz_H3x2y = I_ESP_H2x2yz_G3xy+ABY*I_ESP_G2xyz_G3xy;
    Double I_ESP_G2x2z_H3x2y = I_ESP_H2xy2z_G3xy+ABY*I_ESP_G2x2z_G3xy;
    Double I_ESP_Gx3y_H3x2y = I_ESP_Hx4y_G3xy+ABY*I_ESP_Gx3y_G3xy;
    Double I_ESP_Gx2yz_H3x2y = I_ESP_Hx3yz_G3xy+ABY*I_ESP_Gx2yz_G3xy;
    Double I_ESP_Gxy2z_H3x2y = I_ESP_Hx2y2z_G3xy+ABY*I_ESP_Gxy2z_G3xy;
    Double I_ESP_Gx3z_H3x2y = I_ESP_Hxy3z_G3xy+ABY*I_ESP_Gx3z_G3xy;
    Double I_ESP_G4y_H3x2y = I_ESP_H5y_G3xy+ABY*I_ESP_G4y_G3xy;
    Double I_ESP_G3yz_H3x2y = I_ESP_H4yz_G3xy+ABY*I_ESP_G3yz_G3xy;
    Double I_ESP_G2y2z_H3x2y = I_ESP_H3y2z_G3xy+ABY*I_ESP_G2y2z_G3xy;
    Double I_ESP_Gy3z_H3x2y = I_ESP_H2y3z_G3xy+ABY*I_ESP_Gy3z_G3xy;
    Double I_ESP_G4z_H3x2y = I_ESP_Hy4z_G3xy+ABY*I_ESP_G4z_G3xy;
    Double I_ESP_G4x_H3xyz = I_ESP_H4xz_G3xy+ABZ*I_ESP_G4x_G3xy;
    Double I_ESP_G3xy_H3xyz = I_ESP_H3xyz_G3xy+ABZ*I_ESP_G3xy_G3xy;
    Double I_ESP_G3xz_H3xyz = I_ESP_H3x2z_G3xy+ABZ*I_ESP_G3xz_G3xy;
    Double I_ESP_G2x2y_H3xyz = I_ESP_H2x2yz_G3xy+ABZ*I_ESP_G2x2y_G3xy;
    Double I_ESP_G2xyz_H3xyz = I_ESP_H2xy2z_G3xy+ABZ*I_ESP_G2xyz_G3xy;
    Double I_ESP_G2x2z_H3xyz = I_ESP_H2x3z_G3xy+ABZ*I_ESP_G2x2z_G3xy;
    Double I_ESP_Gx3y_H3xyz = I_ESP_Hx3yz_G3xy+ABZ*I_ESP_Gx3y_G3xy;
    Double I_ESP_Gx2yz_H3xyz = I_ESP_Hx2y2z_G3xy+ABZ*I_ESP_Gx2yz_G3xy;
    Double I_ESP_Gxy2z_H3xyz = I_ESP_Hxy3z_G3xy+ABZ*I_ESP_Gxy2z_G3xy;
    Double I_ESP_Gx3z_H3xyz = I_ESP_Hx4z_G3xy+ABZ*I_ESP_Gx3z_G3xy;
    Double I_ESP_G4y_H3xyz = I_ESP_H4yz_G3xy+ABZ*I_ESP_G4y_G3xy;
    Double I_ESP_G3yz_H3xyz = I_ESP_H3y2z_G3xy+ABZ*I_ESP_G3yz_G3xy;
    Double I_ESP_G2y2z_H3xyz = I_ESP_H2y3z_G3xy+ABZ*I_ESP_G2y2z_G3xy;
    Double I_ESP_Gy3z_H3xyz = I_ESP_Hy4z_G3xy+ABZ*I_ESP_Gy3z_G3xy;
    Double I_ESP_G4z_H3xyz = I_ESP_H5z_G3xy+ABZ*I_ESP_G4z_G3xy;
    Double I_ESP_G4x_H3x2z = I_ESP_H4xz_G3xz+ABZ*I_ESP_G4x_G3xz;
    Double I_ESP_G3xy_H3x2z = I_ESP_H3xyz_G3xz+ABZ*I_ESP_G3xy_G3xz;
    Double I_ESP_G3xz_H3x2z = I_ESP_H3x2z_G3xz+ABZ*I_ESP_G3xz_G3xz;
    Double I_ESP_G2x2y_H3x2z = I_ESP_H2x2yz_G3xz+ABZ*I_ESP_G2x2y_G3xz;
    Double I_ESP_G2xyz_H3x2z = I_ESP_H2xy2z_G3xz+ABZ*I_ESP_G2xyz_G3xz;
    Double I_ESP_G2x2z_H3x2z = I_ESP_H2x3z_G3xz+ABZ*I_ESP_G2x2z_G3xz;
    Double I_ESP_Gx3y_H3x2z = I_ESP_Hx3yz_G3xz+ABZ*I_ESP_Gx3y_G3xz;
    Double I_ESP_Gx2yz_H3x2z = I_ESP_Hx2y2z_G3xz+ABZ*I_ESP_Gx2yz_G3xz;
    Double I_ESP_Gxy2z_H3x2z = I_ESP_Hxy3z_G3xz+ABZ*I_ESP_Gxy2z_G3xz;
    Double I_ESP_Gx3z_H3x2z = I_ESP_Hx4z_G3xz+ABZ*I_ESP_Gx3z_G3xz;
    Double I_ESP_G4y_H3x2z = I_ESP_H4yz_G3xz+ABZ*I_ESP_G4y_G3xz;
    Double I_ESP_G3yz_H3x2z = I_ESP_H3y2z_G3xz+ABZ*I_ESP_G3yz_G3xz;
    Double I_ESP_G2y2z_H3x2z = I_ESP_H2y3z_G3xz+ABZ*I_ESP_G2y2z_G3xz;
    Double I_ESP_Gy3z_H3x2z = I_ESP_Hy4z_G3xz+ABZ*I_ESP_Gy3z_G3xz;
    Double I_ESP_G4z_H3x2z = I_ESP_H5z_G3xz+ABZ*I_ESP_G4z_G3xz;
    Double I_ESP_G4x_H2x3y = I_ESP_H5x_Gx3y+ABX*I_ESP_G4x_Gx3y;
    Double I_ESP_G3xy_H2x3y = I_ESP_H4xy_Gx3y+ABX*I_ESP_G3xy_Gx3y;
    Double I_ESP_G3xz_H2x3y = I_ESP_H4xz_Gx3y+ABX*I_ESP_G3xz_Gx3y;
    Double I_ESP_G2x2y_H2x3y = I_ESP_H3x2y_Gx3y+ABX*I_ESP_G2x2y_Gx3y;
    Double I_ESP_G2xyz_H2x3y = I_ESP_H3xyz_Gx3y+ABX*I_ESP_G2xyz_Gx3y;
    Double I_ESP_G2x2z_H2x3y = I_ESP_H3x2z_Gx3y+ABX*I_ESP_G2x2z_Gx3y;
    Double I_ESP_Gx3y_H2x3y = I_ESP_H2x3y_Gx3y+ABX*I_ESP_Gx3y_Gx3y;
    Double I_ESP_Gx2yz_H2x3y = I_ESP_H2x2yz_Gx3y+ABX*I_ESP_Gx2yz_Gx3y;
    Double I_ESP_Gxy2z_H2x3y = I_ESP_H2xy2z_Gx3y+ABX*I_ESP_Gxy2z_Gx3y;
    Double I_ESP_Gx3z_H2x3y = I_ESP_H2x3z_Gx3y+ABX*I_ESP_Gx3z_Gx3y;
    Double I_ESP_G4y_H2x3y = I_ESP_Hx4y_Gx3y+ABX*I_ESP_G4y_Gx3y;
    Double I_ESP_G3yz_H2x3y = I_ESP_Hx3yz_Gx3y+ABX*I_ESP_G3yz_Gx3y;
    Double I_ESP_G2y2z_H2x3y = I_ESP_Hx2y2z_Gx3y+ABX*I_ESP_G2y2z_Gx3y;
    Double I_ESP_Gy3z_H2x3y = I_ESP_Hxy3z_Gx3y+ABX*I_ESP_Gy3z_Gx3y;
    Double I_ESP_G4z_H2x3y = I_ESP_Hx4z_Gx3y+ABX*I_ESP_G4z_Gx3y;
    Double I_ESP_G4x_H2x2yz = I_ESP_H4xz_G2x2y+ABZ*I_ESP_G4x_G2x2y;
    Double I_ESP_G3xy_H2x2yz = I_ESP_H3xyz_G2x2y+ABZ*I_ESP_G3xy_G2x2y;
    Double I_ESP_G3xz_H2x2yz = I_ESP_H3x2z_G2x2y+ABZ*I_ESP_G3xz_G2x2y;
    Double I_ESP_G2x2y_H2x2yz = I_ESP_H2x2yz_G2x2y+ABZ*I_ESP_G2x2y_G2x2y;
    Double I_ESP_G2xyz_H2x2yz = I_ESP_H2xy2z_G2x2y+ABZ*I_ESP_G2xyz_G2x2y;
    Double I_ESP_G2x2z_H2x2yz = I_ESP_H2x3z_G2x2y+ABZ*I_ESP_G2x2z_G2x2y;
    Double I_ESP_Gx3y_H2x2yz = I_ESP_Hx3yz_G2x2y+ABZ*I_ESP_Gx3y_G2x2y;
    Double I_ESP_Gx2yz_H2x2yz = I_ESP_Hx2y2z_G2x2y+ABZ*I_ESP_Gx2yz_G2x2y;
    Double I_ESP_Gxy2z_H2x2yz = I_ESP_Hxy3z_G2x2y+ABZ*I_ESP_Gxy2z_G2x2y;
    Double I_ESP_Gx3z_H2x2yz = I_ESP_Hx4z_G2x2y+ABZ*I_ESP_Gx3z_G2x2y;
    Double I_ESP_G4y_H2x2yz = I_ESP_H4yz_G2x2y+ABZ*I_ESP_G4y_G2x2y;
    Double I_ESP_G3yz_H2x2yz = I_ESP_H3y2z_G2x2y+ABZ*I_ESP_G3yz_G2x2y;
    Double I_ESP_G2y2z_H2x2yz = I_ESP_H2y3z_G2x2y+ABZ*I_ESP_G2y2z_G2x2y;
    Double I_ESP_Gy3z_H2x2yz = I_ESP_Hy4z_G2x2y+ABZ*I_ESP_Gy3z_G2x2y;
    Double I_ESP_G4z_H2x2yz = I_ESP_H5z_G2x2y+ABZ*I_ESP_G4z_G2x2y;
    Double I_ESP_G4x_H2xy2z = I_ESP_H4xy_G2x2z+ABY*I_ESP_G4x_G2x2z;
    Double I_ESP_G3xy_H2xy2z = I_ESP_H3x2y_G2x2z+ABY*I_ESP_G3xy_G2x2z;
    Double I_ESP_G3xz_H2xy2z = I_ESP_H3xyz_G2x2z+ABY*I_ESP_G3xz_G2x2z;
    Double I_ESP_G2x2y_H2xy2z = I_ESP_H2x3y_G2x2z+ABY*I_ESP_G2x2y_G2x2z;
    Double I_ESP_G2xyz_H2xy2z = I_ESP_H2x2yz_G2x2z+ABY*I_ESP_G2xyz_G2x2z;
    Double I_ESP_G2x2z_H2xy2z = I_ESP_H2xy2z_G2x2z+ABY*I_ESP_G2x2z_G2x2z;
    Double I_ESP_Gx3y_H2xy2z = I_ESP_Hx4y_G2x2z+ABY*I_ESP_Gx3y_G2x2z;
    Double I_ESP_Gx2yz_H2xy2z = I_ESP_Hx3yz_G2x2z+ABY*I_ESP_Gx2yz_G2x2z;
    Double I_ESP_Gxy2z_H2xy2z = I_ESP_Hx2y2z_G2x2z+ABY*I_ESP_Gxy2z_G2x2z;
    Double I_ESP_Gx3z_H2xy2z = I_ESP_Hxy3z_G2x2z+ABY*I_ESP_Gx3z_G2x2z;
    Double I_ESP_G4y_H2xy2z = I_ESP_H5y_G2x2z+ABY*I_ESP_G4y_G2x2z;
    Double I_ESP_G3yz_H2xy2z = I_ESP_H4yz_G2x2z+ABY*I_ESP_G3yz_G2x2z;
    Double I_ESP_G2y2z_H2xy2z = I_ESP_H3y2z_G2x2z+ABY*I_ESP_G2y2z_G2x2z;
    Double I_ESP_Gy3z_H2xy2z = I_ESP_H2y3z_G2x2z+ABY*I_ESP_Gy3z_G2x2z;
    Double I_ESP_G4z_H2xy2z = I_ESP_Hy4z_G2x2z+ABY*I_ESP_G4z_G2x2z;
    Double I_ESP_G4x_H2x3z = I_ESP_H5x_Gx3z+ABX*I_ESP_G4x_Gx3z;
    Double I_ESP_G3xy_H2x3z = I_ESP_H4xy_Gx3z+ABX*I_ESP_G3xy_Gx3z;
    Double I_ESP_G3xz_H2x3z = I_ESP_H4xz_Gx3z+ABX*I_ESP_G3xz_Gx3z;
    Double I_ESP_G2x2y_H2x3z = I_ESP_H3x2y_Gx3z+ABX*I_ESP_G2x2y_Gx3z;
    Double I_ESP_G2xyz_H2x3z = I_ESP_H3xyz_Gx3z+ABX*I_ESP_G2xyz_Gx3z;
    Double I_ESP_G2x2z_H2x3z = I_ESP_H3x2z_Gx3z+ABX*I_ESP_G2x2z_Gx3z;
    Double I_ESP_Gx3y_H2x3z = I_ESP_H2x3y_Gx3z+ABX*I_ESP_Gx3y_Gx3z;
    Double I_ESP_Gx2yz_H2x3z = I_ESP_H2x2yz_Gx3z+ABX*I_ESP_Gx2yz_Gx3z;
    Double I_ESP_Gxy2z_H2x3z = I_ESP_H2xy2z_Gx3z+ABX*I_ESP_Gxy2z_Gx3z;
    Double I_ESP_Gx3z_H2x3z = I_ESP_H2x3z_Gx3z+ABX*I_ESP_Gx3z_Gx3z;
    Double I_ESP_G4y_H2x3z = I_ESP_Hx4y_Gx3z+ABX*I_ESP_G4y_Gx3z;
    Double I_ESP_G3yz_H2x3z = I_ESP_Hx3yz_Gx3z+ABX*I_ESP_G3yz_Gx3z;
    Double I_ESP_G2y2z_H2x3z = I_ESP_Hx2y2z_Gx3z+ABX*I_ESP_G2y2z_Gx3z;
    Double I_ESP_Gy3z_H2x3z = I_ESP_Hxy3z_Gx3z+ABX*I_ESP_Gy3z_Gx3z;
    Double I_ESP_G4z_H2x3z = I_ESP_Hx4z_Gx3z+ABX*I_ESP_G4z_Gx3z;
    Double I_ESP_G4x_Hx4y = I_ESP_H5x_G4y+ABX*I_ESP_G4x_G4y;
    Double I_ESP_G3xy_Hx4y = I_ESP_H4xy_G4y+ABX*I_ESP_G3xy_G4y;
    Double I_ESP_G3xz_Hx4y = I_ESP_H4xz_G4y+ABX*I_ESP_G3xz_G4y;
    Double I_ESP_G2x2y_Hx4y = I_ESP_H3x2y_G4y+ABX*I_ESP_G2x2y_G4y;
    Double I_ESP_G2xyz_Hx4y = I_ESP_H3xyz_G4y+ABX*I_ESP_G2xyz_G4y;
    Double I_ESP_G2x2z_Hx4y = I_ESP_H3x2z_G4y+ABX*I_ESP_G2x2z_G4y;
    Double I_ESP_Gx3y_Hx4y = I_ESP_H2x3y_G4y+ABX*I_ESP_Gx3y_G4y;
    Double I_ESP_Gx2yz_Hx4y = I_ESP_H2x2yz_G4y+ABX*I_ESP_Gx2yz_G4y;
    Double I_ESP_Gxy2z_Hx4y = I_ESP_H2xy2z_G4y+ABX*I_ESP_Gxy2z_G4y;
    Double I_ESP_Gx3z_Hx4y = I_ESP_H2x3z_G4y+ABX*I_ESP_Gx3z_G4y;
    Double I_ESP_G4y_Hx4y = I_ESP_Hx4y_G4y+ABX*I_ESP_G4y_G4y;
    Double I_ESP_G3yz_Hx4y = I_ESP_Hx3yz_G4y+ABX*I_ESP_G3yz_G4y;
    Double I_ESP_G2y2z_Hx4y = I_ESP_Hx2y2z_G4y+ABX*I_ESP_G2y2z_G4y;
    Double I_ESP_Gy3z_Hx4y = I_ESP_Hxy3z_G4y+ABX*I_ESP_Gy3z_G4y;
    Double I_ESP_G4z_Hx4y = I_ESP_Hx4z_G4y+ABX*I_ESP_G4z_G4y;
    Double I_ESP_G4x_Hx3yz = I_ESP_H4xz_Gx3y+ABZ*I_ESP_G4x_Gx3y;
    Double I_ESP_G3xy_Hx3yz = I_ESP_H3xyz_Gx3y+ABZ*I_ESP_G3xy_Gx3y;
    Double I_ESP_G3xz_Hx3yz = I_ESP_H3x2z_Gx3y+ABZ*I_ESP_G3xz_Gx3y;
    Double I_ESP_G2x2y_Hx3yz = I_ESP_H2x2yz_Gx3y+ABZ*I_ESP_G2x2y_Gx3y;
    Double I_ESP_G2xyz_Hx3yz = I_ESP_H2xy2z_Gx3y+ABZ*I_ESP_G2xyz_Gx3y;
    Double I_ESP_G2x2z_Hx3yz = I_ESP_H2x3z_Gx3y+ABZ*I_ESP_G2x2z_Gx3y;
    Double I_ESP_Gx3y_Hx3yz = I_ESP_Hx3yz_Gx3y+ABZ*I_ESP_Gx3y_Gx3y;
    Double I_ESP_Gx2yz_Hx3yz = I_ESP_Hx2y2z_Gx3y+ABZ*I_ESP_Gx2yz_Gx3y;
    Double I_ESP_Gxy2z_Hx3yz = I_ESP_Hxy3z_Gx3y+ABZ*I_ESP_Gxy2z_Gx3y;
    Double I_ESP_Gx3z_Hx3yz = I_ESP_Hx4z_Gx3y+ABZ*I_ESP_Gx3z_Gx3y;
    Double I_ESP_G4y_Hx3yz = I_ESP_H4yz_Gx3y+ABZ*I_ESP_G4y_Gx3y;
    Double I_ESP_G3yz_Hx3yz = I_ESP_H3y2z_Gx3y+ABZ*I_ESP_G3yz_Gx3y;
    Double I_ESP_G2y2z_Hx3yz = I_ESP_H2y3z_Gx3y+ABZ*I_ESP_G2y2z_Gx3y;
    Double I_ESP_Gy3z_Hx3yz = I_ESP_Hy4z_Gx3y+ABZ*I_ESP_Gy3z_Gx3y;
    Double I_ESP_G4z_Hx3yz = I_ESP_H5z_Gx3y+ABZ*I_ESP_G4z_Gx3y;
    Double I_ESP_G4x_Hx2y2z = I_ESP_H5x_G2y2z+ABX*I_ESP_G4x_G2y2z;
    Double I_ESP_G3xy_Hx2y2z = I_ESP_H4xy_G2y2z+ABX*I_ESP_G3xy_G2y2z;
    Double I_ESP_G3xz_Hx2y2z = I_ESP_H4xz_G2y2z+ABX*I_ESP_G3xz_G2y2z;
    Double I_ESP_G2x2y_Hx2y2z = I_ESP_H3x2y_G2y2z+ABX*I_ESP_G2x2y_G2y2z;
    Double I_ESP_G2xyz_Hx2y2z = I_ESP_H3xyz_G2y2z+ABX*I_ESP_G2xyz_G2y2z;
    Double I_ESP_G2x2z_Hx2y2z = I_ESP_H3x2z_G2y2z+ABX*I_ESP_G2x2z_G2y2z;
    Double I_ESP_Gx3y_Hx2y2z = I_ESP_H2x3y_G2y2z+ABX*I_ESP_Gx3y_G2y2z;
    Double I_ESP_Gx2yz_Hx2y2z = I_ESP_H2x2yz_G2y2z+ABX*I_ESP_Gx2yz_G2y2z;
    Double I_ESP_Gxy2z_Hx2y2z = I_ESP_H2xy2z_G2y2z+ABX*I_ESP_Gxy2z_G2y2z;
    Double I_ESP_Gx3z_Hx2y2z = I_ESP_H2x3z_G2y2z+ABX*I_ESP_Gx3z_G2y2z;
    Double I_ESP_G4y_Hx2y2z = I_ESP_Hx4y_G2y2z+ABX*I_ESP_G4y_G2y2z;
    Double I_ESP_G3yz_Hx2y2z = I_ESP_Hx3yz_G2y2z+ABX*I_ESP_G3yz_G2y2z;
    Double I_ESP_G2y2z_Hx2y2z = I_ESP_Hx2y2z_G2y2z+ABX*I_ESP_G2y2z_G2y2z;
    Double I_ESP_Gy3z_Hx2y2z = I_ESP_Hxy3z_G2y2z+ABX*I_ESP_Gy3z_G2y2z;
    Double I_ESP_G4z_Hx2y2z = I_ESP_Hx4z_G2y2z+ABX*I_ESP_G4z_G2y2z;
    Double I_ESP_G4x_Hxy3z = I_ESP_H4xy_Gx3z+ABY*I_ESP_G4x_Gx3z;
    Double I_ESP_G3xy_Hxy3z = I_ESP_H3x2y_Gx3z+ABY*I_ESP_G3xy_Gx3z;
    Double I_ESP_G3xz_Hxy3z = I_ESP_H3xyz_Gx3z+ABY*I_ESP_G3xz_Gx3z;
    Double I_ESP_G2x2y_Hxy3z = I_ESP_H2x3y_Gx3z+ABY*I_ESP_G2x2y_Gx3z;
    Double I_ESP_G2xyz_Hxy3z = I_ESP_H2x2yz_Gx3z+ABY*I_ESP_G2xyz_Gx3z;
    Double I_ESP_G2x2z_Hxy3z = I_ESP_H2xy2z_Gx3z+ABY*I_ESP_G2x2z_Gx3z;
    Double I_ESP_Gx3y_Hxy3z = I_ESP_Hx4y_Gx3z+ABY*I_ESP_Gx3y_Gx3z;
    Double I_ESP_Gx2yz_Hxy3z = I_ESP_Hx3yz_Gx3z+ABY*I_ESP_Gx2yz_Gx3z;
    Double I_ESP_Gxy2z_Hxy3z = I_ESP_Hx2y2z_Gx3z+ABY*I_ESP_Gxy2z_Gx3z;
    Double I_ESP_Gx3z_Hxy3z = I_ESP_Hxy3z_Gx3z+ABY*I_ESP_Gx3z_Gx3z;
    Double I_ESP_G4y_Hxy3z = I_ESP_H5y_Gx3z+ABY*I_ESP_G4y_Gx3z;
    Double I_ESP_G3yz_Hxy3z = I_ESP_H4yz_Gx3z+ABY*I_ESP_G3yz_Gx3z;
    Double I_ESP_G2y2z_Hxy3z = I_ESP_H3y2z_Gx3z+ABY*I_ESP_G2y2z_Gx3z;
    Double I_ESP_Gy3z_Hxy3z = I_ESP_H2y3z_Gx3z+ABY*I_ESP_Gy3z_Gx3z;
    Double I_ESP_G4z_Hxy3z = I_ESP_Hy4z_Gx3z+ABY*I_ESP_G4z_Gx3z;
    Double I_ESP_G4x_Hx4z = I_ESP_H5x_G4z+ABX*I_ESP_G4x_G4z;
    Double I_ESP_G3xy_Hx4z = I_ESP_H4xy_G4z+ABX*I_ESP_G3xy_G4z;
    Double I_ESP_G3xz_Hx4z = I_ESP_H4xz_G4z+ABX*I_ESP_G3xz_G4z;
    Double I_ESP_G2x2y_Hx4z = I_ESP_H3x2y_G4z+ABX*I_ESP_G2x2y_G4z;
    Double I_ESP_G2xyz_Hx4z = I_ESP_H3xyz_G4z+ABX*I_ESP_G2xyz_G4z;
    Double I_ESP_G2x2z_Hx4z = I_ESP_H3x2z_G4z+ABX*I_ESP_G2x2z_G4z;
    Double I_ESP_Gx3y_Hx4z = I_ESP_H2x3y_G4z+ABX*I_ESP_Gx3y_G4z;
    Double I_ESP_Gx2yz_Hx4z = I_ESP_H2x2yz_G4z+ABX*I_ESP_Gx2yz_G4z;
    Double I_ESP_Gxy2z_Hx4z = I_ESP_H2xy2z_G4z+ABX*I_ESP_Gxy2z_G4z;
    Double I_ESP_Gx3z_Hx4z = I_ESP_H2x3z_G4z+ABX*I_ESP_Gx3z_G4z;
    Double I_ESP_G4y_Hx4z = I_ESP_Hx4y_G4z+ABX*I_ESP_G4y_G4z;
    Double I_ESP_G3yz_Hx4z = I_ESP_Hx3yz_G4z+ABX*I_ESP_G3yz_G4z;
    Double I_ESP_G2y2z_Hx4z = I_ESP_Hx2y2z_G4z+ABX*I_ESP_G2y2z_G4z;
    Double I_ESP_Gy3z_Hx4z = I_ESP_Hxy3z_G4z+ABX*I_ESP_Gy3z_G4z;
    Double I_ESP_G4z_Hx4z = I_ESP_Hx4z_G4z+ABX*I_ESP_G4z_G4z;
    Double I_ESP_G4x_H5y = I_ESP_H4xy_G4y+ABY*I_ESP_G4x_G4y;
    Double I_ESP_G3xy_H5y = I_ESP_H3x2y_G4y+ABY*I_ESP_G3xy_G4y;
    Double I_ESP_G3xz_H5y = I_ESP_H3xyz_G4y+ABY*I_ESP_G3xz_G4y;
    Double I_ESP_G2x2y_H5y = I_ESP_H2x3y_G4y+ABY*I_ESP_G2x2y_G4y;
    Double I_ESP_G2xyz_H5y = I_ESP_H2x2yz_G4y+ABY*I_ESP_G2xyz_G4y;
    Double I_ESP_G2x2z_H5y = I_ESP_H2xy2z_G4y+ABY*I_ESP_G2x2z_G4y;
    Double I_ESP_Gx3y_H5y = I_ESP_Hx4y_G4y+ABY*I_ESP_Gx3y_G4y;
    Double I_ESP_Gx2yz_H5y = I_ESP_Hx3yz_G4y+ABY*I_ESP_Gx2yz_G4y;
    Double I_ESP_Gxy2z_H5y = I_ESP_Hx2y2z_G4y+ABY*I_ESP_Gxy2z_G4y;
    Double I_ESP_Gx3z_H5y = I_ESP_Hxy3z_G4y+ABY*I_ESP_Gx3z_G4y;
    Double I_ESP_G4y_H5y = I_ESP_H5y_G4y+ABY*I_ESP_G4y_G4y;
    Double I_ESP_G3yz_H5y = I_ESP_H4yz_G4y+ABY*I_ESP_G3yz_G4y;
    Double I_ESP_G2y2z_H5y = I_ESP_H3y2z_G4y+ABY*I_ESP_G2y2z_G4y;
    Double I_ESP_Gy3z_H5y = I_ESP_H2y3z_G4y+ABY*I_ESP_Gy3z_G4y;
    Double I_ESP_G4z_H5y = I_ESP_Hy4z_G4y+ABY*I_ESP_G4z_G4y;
    Double I_ESP_G4x_H4yz = I_ESP_H4xz_G4y+ABZ*I_ESP_G4x_G4y;
    Double I_ESP_G3xy_H4yz = I_ESP_H3xyz_G4y+ABZ*I_ESP_G3xy_G4y;
    Double I_ESP_G3xz_H4yz = I_ESP_H3x2z_G4y+ABZ*I_ESP_G3xz_G4y;
    Double I_ESP_G2x2y_H4yz = I_ESP_H2x2yz_G4y+ABZ*I_ESP_G2x2y_G4y;
    Double I_ESP_G2xyz_H4yz = I_ESP_H2xy2z_G4y+ABZ*I_ESP_G2xyz_G4y;
    Double I_ESP_G2x2z_H4yz = I_ESP_H2x3z_G4y+ABZ*I_ESP_G2x2z_G4y;
    Double I_ESP_Gx3y_H4yz = I_ESP_Hx3yz_G4y+ABZ*I_ESP_Gx3y_G4y;
    Double I_ESP_Gx2yz_H4yz = I_ESP_Hx2y2z_G4y+ABZ*I_ESP_Gx2yz_G4y;
    Double I_ESP_Gxy2z_H4yz = I_ESP_Hxy3z_G4y+ABZ*I_ESP_Gxy2z_G4y;
    Double I_ESP_Gx3z_H4yz = I_ESP_Hx4z_G4y+ABZ*I_ESP_Gx3z_G4y;
    Double I_ESP_G4y_H4yz = I_ESP_H4yz_G4y+ABZ*I_ESP_G4y_G4y;
    Double I_ESP_G3yz_H4yz = I_ESP_H3y2z_G4y+ABZ*I_ESP_G3yz_G4y;
    Double I_ESP_G2y2z_H4yz = I_ESP_H2y3z_G4y+ABZ*I_ESP_G2y2z_G4y;
    Double I_ESP_Gy3z_H4yz = I_ESP_Hy4z_G4y+ABZ*I_ESP_Gy3z_G4y;
    Double I_ESP_G4z_H4yz = I_ESP_H5z_G4y+ABZ*I_ESP_G4z_G4y;
    Double I_ESP_G4x_H3y2z = I_ESP_H4xz_G3yz+ABZ*I_ESP_G4x_G3yz;
    Double I_ESP_G3xy_H3y2z = I_ESP_H3xyz_G3yz+ABZ*I_ESP_G3xy_G3yz;
    Double I_ESP_G3xz_H3y2z = I_ESP_H3x2z_G3yz+ABZ*I_ESP_G3xz_G3yz;
    Double I_ESP_G2x2y_H3y2z = I_ESP_H2x2yz_G3yz+ABZ*I_ESP_G2x2y_G3yz;
    Double I_ESP_G2xyz_H3y2z = I_ESP_H2xy2z_G3yz+ABZ*I_ESP_G2xyz_G3yz;
    Double I_ESP_G2x2z_H3y2z = I_ESP_H2x3z_G3yz+ABZ*I_ESP_G2x2z_G3yz;
    Double I_ESP_Gx3y_H3y2z = I_ESP_Hx3yz_G3yz+ABZ*I_ESP_Gx3y_G3yz;
    Double I_ESP_Gx2yz_H3y2z = I_ESP_Hx2y2z_G3yz+ABZ*I_ESP_Gx2yz_G3yz;
    Double I_ESP_Gxy2z_H3y2z = I_ESP_Hxy3z_G3yz+ABZ*I_ESP_Gxy2z_G3yz;
    Double I_ESP_Gx3z_H3y2z = I_ESP_Hx4z_G3yz+ABZ*I_ESP_Gx3z_G3yz;
    Double I_ESP_G4y_H3y2z = I_ESP_H4yz_G3yz+ABZ*I_ESP_G4y_G3yz;
    Double I_ESP_G3yz_H3y2z = I_ESP_H3y2z_G3yz+ABZ*I_ESP_G3yz_G3yz;
    Double I_ESP_G2y2z_H3y2z = I_ESP_H2y3z_G3yz+ABZ*I_ESP_G2y2z_G3yz;
    Double I_ESP_Gy3z_H3y2z = I_ESP_Hy4z_G3yz+ABZ*I_ESP_Gy3z_G3yz;
    Double I_ESP_G4z_H3y2z = I_ESP_H5z_G3yz+ABZ*I_ESP_G4z_G3yz;
    Double I_ESP_G4x_H2y3z = I_ESP_H4xy_Gy3z+ABY*I_ESP_G4x_Gy3z;
    Double I_ESP_G3xy_H2y3z = I_ESP_H3x2y_Gy3z+ABY*I_ESP_G3xy_Gy3z;
    Double I_ESP_G3xz_H2y3z = I_ESP_H3xyz_Gy3z+ABY*I_ESP_G3xz_Gy3z;
    Double I_ESP_G2x2y_H2y3z = I_ESP_H2x3y_Gy3z+ABY*I_ESP_G2x2y_Gy3z;
    Double I_ESP_G2xyz_H2y3z = I_ESP_H2x2yz_Gy3z+ABY*I_ESP_G2xyz_Gy3z;
    Double I_ESP_G2x2z_H2y3z = I_ESP_H2xy2z_Gy3z+ABY*I_ESP_G2x2z_Gy3z;
    Double I_ESP_Gx3y_H2y3z = I_ESP_Hx4y_Gy3z+ABY*I_ESP_Gx3y_Gy3z;
    Double I_ESP_Gx2yz_H2y3z = I_ESP_Hx3yz_Gy3z+ABY*I_ESP_Gx2yz_Gy3z;
    Double I_ESP_Gxy2z_H2y3z = I_ESP_Hx2y2z_Gy3z+ABY*I_ESP_Gxy2z_Gy3z;
    Double I_ESP_Gx3z_H2y3z = I_ESP_Hxy3z_Gy3z+ABY*I_ESP_Gx3z_Gy3z;
    Double I_ESP_G4y_H2y3z = I_ESP_H5y_Gy3z+ABY*I_ESP_G4y_Gy3z;
    Double I_ESP_G3yz_H2y3z = I_ESP_H4yz_Gy3z+ABY*I_ESP_G3yz_Gy3z;
    Double I_ESP_G2y2z_H2y3z = I_ESP_H3y2z_Gy3z+ABY*I_ESP_G2y2z_Gy3z;
    Double I_ESP_Gy3z_H2y3z = I_ESP_H2y3z_Gy3z+ABY*I_ESP_Gy3z_Gy3z;
    Double I_ESP_G4z_H2y3z = I_ESP_Hy4z_Gy3z+ABY*I_ESP_G4z_Gy3z;
    Double I_ESP_G4x_Hy4z = I_ESP_H4xy_G4z+ABY*I_ESP_G4x_G4z;
    Double I_ESP_G3xy_Hy4z = I_ESP_H3x2y_G4z+ABY*I_ESP_G3xy_G4z;
    Double I_ESP_G3xz_Hy4z = I_ESP_H3xyz_G4z+ABY*I_ESP_G3xz_G4z;
    Double I_ESP_G2x2y_Hy4z = I_ESP_H2x3y_G4z+ABY*I_ESP_G2x2y_G4z;
    Double I_ESP_G2xyz_Hy4z = I_ESP_H2x2yz_G4z+ABY*I_ESP_G2xyz_G4z;
    Double I_ESP_G2x2z_Hy4z = I_ESP_H2xy2z_G4z+ABY*I_ESP_G2x2z_G4z;
    Double I_ESP_Gx3y_Hy4z = I_ESP_Hx4y_G4z+ABY*I_ESP_Gx3y_G4z;
    Double I_ESP_Gx2yz_Hy4z = I_ESP_Hx3yz_G4z+ABY*I_ESP_Gx2yz_G4z;
    Double I_ESP_Gxy2z_Hy4z = I_ESP_Hx2y2z_G4z+ABY*I_ESP_Gxy2z_G4z;
    Double I_ESP_Gx3z_Hy4z = I_ESP_Hxy3z_G4z+ABY*I_ESP_Gx3z_G4z;
    Double I_ESP_G4y_Hy4z = I_ESP_H5y_G4z+ABY*I_ESP_G4y_G4z;
    Double I_ESP_G3yz_Hy4z = I_ESP_H4yz_G4z+ABY*I_ESP_G3yz_G4z;
    Double I_ESP_G2y2z_Hy4z = I_ESP_H3y2z_G4z+ABY*I_ESP_G2y2z_G4z;
    Double I_ESP_Gy3z_Hy4z = I_ESP_H2y3z_G4z+ABY*I_ESP_Gy3z_G4z;
    Double I_ESP_G4z_Hy4z = I_ESP_Hy4z_G4z+ABY*I_ESP_G4z_G4z;
    Double I_ESP_G4x_H5z = I_ESP_H4xz_G4z+ABZ*I_ESP_G4x_G4z;
    Double I_ESP_G3xy_H5z = I_ESP_H3xyz_G4z+ABZ*I_ESP_G3xy_G4z;
    Double I_ESP_G3xz_H5z = I_ESP_H3x2z_G4z+ABZ*I_ESP_G3xz_G4z;
    Double I_ESP_G2x2y_H5z = I_ESP_H2x2yz_G4z+ABZ*I_ESP_G2x2y_G4z;
    Double I_ESP_G2xyz_H5z = I_ESP_H2xy2z_G4z+ABZ*I_ESP_G2xyz_G4z;
    Double I_ESP_G2x2z_H5z = I_ESP_H2x3z_G4z+ABZ*I_ESP_G2x2z_G4z;
    Double I_ESP_Gx3y_H5z = I_ESP_Hx3yz_G4z+ABZ*I_ESP_Gx3y_G4z;
    Double I_ESP_Gx2yz_H5z = I_ESP_Hx2y2z_G4z+ABZ*I_ESP_Gx2yz_G4z;
    Double I_ESP_Gxy2z_H5z = I_ESP_Hxy3z_G4z+ABZ*I_ESP_Gxy2z_G4z;
    Double I_ESP_Gx3z_H5z = I_ESP_Hx4z_G4z+ABZ*I_ESP_Gx3z_G4z;
    Double I_ESP_G4y_H5z = I_ESP_H4yz_G4z+ABZ*I_ESP_G4y_G4z;
    Double I_ESP_G3yz_H5z = I_ESP_H3y2z_G4z+ABZ*I_ESP_G3yz_G4z;
    Double I_ESP_G2y2z_H5z = I_ESP_H2y3z_G4z+ABZ*I_ESP_G2y2z_G4z;
    Double I_ESP_Gy3z_H5z = I_ESP_Hy4z_G4z+ABZ*I_ESP_Gy3z_G4z;
    Double I_ESP_G4z_H5z = I_ESP_H5z_G4z+ABZ*I_ESP_G4z_G4z;

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
     * totally 0 integrals are omitted 
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
    Double I_ESP_L8x_Py_a = I_ESP_M8xy_S_a+ABY*I_ESP_L8x_S_a;
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
    Double I_ESP_L8x_Pz_a = I_ESP_M8xz_S_a+ABZ*I_ESP_L8x_S_a;
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
    Double I_ESP_L8y_Pz_a = I_ESP_M8yz_S_a+ABZ*I_ESP_L8y_S_a;
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
     * totally 112 integrals are omitted 
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
     * totally 15 integrals are omitted 
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
    Double I_ESP_M8yz_Px_a = I_ESP_Nx8yz_S_a+ABX*I_ESP_M8yz_S_a;
    Double I_ESP_M7y2z_Px_a = I_ESP_Nx7y2z_S_a+ABX*I_ESP_M7y2z_S_a;
    Double I_ESP_M6y3z_Px_a = I_ESP_Nx6y3z_S_a+ABX*I_ESP_M6y3z_S_a;
    Double I_ESP_M5y4z_Px_a = I_ESP_Nx5y4z_S_a+ABX*I_ESP_M5y4z_S_a;
    Double I_ESP_M4y5z_Px_a = I_ESP_Nx4y5z_S_a+ABX*I_ESP_M4y5z_S_a;
    Double I_ESP_M3y6z_Px_a = I_ESP_Nx3y6z_S_a+ABX*I_ESP_M3y6z_S_a;
    Double I_ESP_M2y7z_Px_a = I_ESP_Nx2y7z_S_a+ABX*I_ESP_M2y7z_S_a;
    Double I_ESP_My8z_Px_a = I_ESP_Nxy8z_S_a+ABX*I_ESP_My8z_S_a;
    Double I_ESP_M8xy_Py_a = I_ESP_N8x2y_S_a+ABY*I_ESP_M8xy_S_a;
    Double I_ESP_M7x2y_Py_a = I_ESP_N7x3y_S_a+ABY*I_ESP_M7x2y_S_a;
    Double I_ESP_M7xyz_Py_a = I_ESP_N7x2yz_S_a+ABY*I_ESP_M7xyz_S_a;
    Double I_ESP_M7x2z_Py_a = I_ESP_N7xy2z_S_a+ABY*I_ESP_M7x2z_S_a;
    Double I_ESP_M6x3y_Py_a = I_ESP_N6x4y_S_a+ABY*I_ESP_M6x3y_S_a;
    Double I_ESP_M6x2yz_Py_a = I_ESP_N6x3yz_S_a+ABY*I_ESP_M6x2yz_S_a;
    Double I_ESP_M6xy2z_Py_a = I_ESP_N6x2y2z_S_a+ABY*I_ESP_M6xy2z_S_a;
    Double I_ESP_M6x3z_Py_a = I_ESP_N6xy3z_S_a+ABY*I_ESP_M6x3z_S_a;
    Double I_ESP_M5x4y_Py_a = I_ESP_N5x5y_S_a+ABY*I_ESP_M5x4y_S_a;
    Double I_ESP_M5x3yz_Py_a = I_ESP_N5x4yz_S_a+ABY*I_ESP_M5x3yz_S_a;
    Double I_ESP_M5x2y2z_Py_a = I_ESP_N5x3y2z_S_a+ABY*I_ESP_M5x2y2z_S_a;
    Double I_ESP_M5xy3z_Py_a = I_ESP_N5x2y3z_S_a+ABY*I_ESP_M5xy3z_S_a;
    Double I_ESP_M5x4z_Py_a = I_ESP_N5xy4z_S_a+ABY*I_ESP_M5x4z_S_a;
    Double I_ESP_M4x5y_Py_a = I_ESP_N4x6y_S_a+ABY*I_ESP_M4x5y_S_a;
    Double I_ESP_M4x4yz_Py_a = I_ESP_N4x5yz_S_a+ABY*I_ESP_M4x4yz_S_a;
    Double I_ESP_M4x3y2z_Py_a = I_ESP_N4x4y2z_S_a+ABY*I_ESP_M4x3y2z_S_a;
    Double I_ESP_M4x2y3z_Py_a = I_ESP_N4x3y3z_S_a+ABY*I_ESP_M4x2y3z_S_a;
    Double I_ESP_M4xy4z_Py_a = I_ESP_N4x2y4z_S_a+ABY*I_ESP_M4xy4z_S_a;
    Double I_ESP_M4x5z_Py_a = I_ESP_N4xy5z_S_a+ABY*I_ESP_M4x5z_S_a;
    Double I_ESP_M3x6y_Py_a = I_ESP_N3x7y_S_a+ABY*I_ESP_M3x6y_S_a;
    Double I_ESP_M3x5yz_Py_a = I_ESP_N3x6yz_S_a+ABY*I_ESP_M3x5yz_S_a;
    Double I_ESP_M3x4y2z_Py_a = I_ESP_N3x5y2z_S_a+ABY*I_ESP_M3x4y2z_S_a;
    Double I_ESP_M3x3y3z_Py_a = I_ESP_N3x4y3z_S_a+ABY*I_ESP_M3x3y3z_S_a;
    Double I_ESP_M3x2y4z_Py_a = I_ESP_N3x3y4z_S_a+ABY*I_ESP_M3x2y4z_S_a;
    Double I_ESP_M3xy5z_Py_a = I_ESP_N3x2y5z_S_a+ABY*I_ESP_M3xy5z_S_a;
    Double I_ESP_M3x6z_Py_a = I_ESP_N3xy6z_S_a+ABY*I_ESP_M3x6z_S_a;
    Double I_ESP_M2x7y_Py_a = I_ESP_N2x8y_S_a+ABY*I_ESP_M2x7y_S_a;
    Double I_ESP_M2x6yz_Py_a = I_ESP_N2x7yz_S_a+ABY*I_ESP_M2x6yz_S_a;
    Double I_ESP_M2x5y2z_Py_a = I_ESP_N2x6y2z_S_a+ABY*I_ESP_M2x5y2z_S_a;
    Double I_ESP_M2x4y3z_Py_a = I_ESP_N2x5y3z_S_a+ABY*I_ESP_M2x4y3z_S_a;
    Double I_ESP_M2x3y4z_Py_a = I_ESP_N2x4y4z_S_a+ABY*I_ESP_M2x3y4z_S_a;
    Double I_ESP_M2x2y5z_Py_a = I_ESP_N2x3y5z_S_a+ABY*I_ESP_M2x2y5z_S_a;
    Double I_ESP_M2xy6z_Py_a = I_ESP_N2x2y6z_S_a+ABY*I_ESP_M2xy6z_S_a;
    Double I_ESP_M2x7z_Py_a = I_ESP_N2xy7z_S_a+ABY*I_ESP_M2x7z_S_a;
    Double I_ESP_Mx8y_Py_a = I_ESP_Nx9y_S_a+ABY*I_ESP_Mx8y_S_a;
    Double I_ESP_Mx7yz_Py_a = I_ESP_Nx8yz_S_a+ABY*I_ESP_Mx7yz_S_a;
    Double I_ESP_Mx6y2z_Py_a = I_ESP_Nx7y2z_S_a+ABY*I_ESP_Mx6y2z_S_a;
    Double I_ESP_Mx5y3z_Py_a = I_ESP_Nx6y3z_S_a+ABY*I_ESP_Mx5y3z_S_a;
    Double I_ESP_Mx4y4z_Py_a = I_ESP_Nx5y4z_S_a+ABY*I_ESP_Mx4y4z_S_a;
    Double I_ESP_Mx3y5z_Py_a = I_ESP_Nx4y5z_S_a+ABY*I_ESP_Mx3y5z_S_a;
    Double I_ESP_Mx2y6z_Py_a = I_ESP_Nx3y6z_S_a+ABY*I_ESP_Mx2y6z_S_a;
    Double I_ESP_Mxy7z_Py_a = I_ESP_Nx2y7z_S_a+ABY*I_ESP_Mxy7z_S_a;
    Double I_ESP_Mx8z_Py_a = I_ESP_Nxy8z_S_a+ABY*I_ESP_Mx8z_S_a;
    Double I_ESP_M9y_Py_a = I_ESP_N10y_S_a+ABY*I_ESP_M9y_S_a;
    Double I_ESP_M8yz_Py_a = I_ESP_N9yz_S_a+ABY*I_ESP_M8yz_S_a;
    Double I_ESP_M7y2z_Py_a = I_ESP_N8y2z_S_a+ABY*I_ESP_M7y2z_S_a;
    Double I_ESP_M6y3z_Py_a = I_ESP_N7y3z_S_a+ABY*I_ESP_M6y3z_S_a;
    Double I_ESP_M5y4z_Py_a = I_ESP_N6y4z_S_a+ABY*I_ESP_M5y4z_S_a;
    Double I_ESP_M4y5z_Py_a = I_ESP_N5y5z_S_a+ABY*I_ESP_M4y5z_S_a;
    Double I_ESP_M3y6z_Py_a = I_ESP_N4y6z_S_a+ABY*I_ESP_M3y6z_S_a;
    Double I_ESP_M2y7z_Py_a = I_ESP_N3y7z_S_a+ABY*I_ESP_M2y7z_S_a;
    Double I_ESP_My8z_Py_a = I_ESP_N2y8z_S_a+ABY*I_ESP_My8z_S_a;
    Double I_ESP_M8xz_Pz_a = I_ESP_N8x2z_S_a+ABZ*I_ESP_M8xz_S_a;
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
    Double I_ESP_M8yz_Pz_a = I_ESP_N8y2z_S_a+ABZ*I_ESP_M8yz_S_a;
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
     * totally 135 integrals are omitted 
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
    Double I_ESP_L8x_D2y_a = I_ESP_M8xy_Py_a+ABY*I_ESP_L8x_Py_a;
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
    Double I_ESP_L8x_D2z_a = I_ESP_M8xz_Pz_a+ABZ*I_ESP_L8x_Pz_a;
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
    Double I_ESP_L8y_D2z_a = I_ESP_M8yz_Pz_a+ABZ*I_ESP_L8y_Pz_a;
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
     * totally 147 integrals are omitted 
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
    Double I_ESP_K6xy_F2xz_a = I_ESP_L6xyz_D2x_a+ABZ*I_ESP_K6xy_D2x_a;
    Double I_ESP_K6xz_F2xz_a = I_ESP_L6x2z_D2x_a+ABZ*I_ESP_K6xz_D2x_a;
    Double I_ESP_K5x2y_F2xz_a = I_ESP_L5x2yz_D2x_a+ABZ*I_ESP_K5x2y_D2x_a;
    Double I_ESP_K5xyz_F2xz_a = I_ESP_L5xy2z_D2x_a+ABZ*I_ESP_K5xyz_D2x_a;
    Double I_ESP_K5x2z_F2xz_a = I_ESP_L5x3z_D2x_a+ABZ*I_ESP_K5x2z_D2x_a;
    Double I_ESP_K4x3y_F2xz_a = I_ESP_L4x3yz_D2x_a+ABZ*I_ESP_K4x3y_D2x_a;
    Double I_ESP_K4x2yz_F2xz_a = I_ESP_L4x2y2z_D2x_a+ABZ*I_ESP_K4x2yz_D2x_a;
    Double I_ESP_K4xy2z_F2xz_a = I_ESP_L4xy3z_D2x_a+ABZ*I_ESP_K4xy2z_D2x_a;
    Double I_ESP_K4x3z_F2xz_a = I_ESP_L4x4z_D2x_a+ABZ*I_ESP_K4x3z_D2x_a;
    Double I_ESP_K3x4y_F2xz_a = I_ESP_L3x4yz_D2x_a+ABZ*I_ESP_K3x4y_D2x_a;
    Double I_ESP_K3x3yz_F2xz_a = I_ESP_L3x3y2z_D2x_a+ABZ*I_ESP_K3x3yz_D2x_a;
    Double I_ESP_K3x2y2z_F2xz_a = I_ESP_L3x2y3z_D2x_a+ABZ*I_ESP_K3x2y2z_D2x_a;
    Double I_ESP_K3xy3z_F2xz_a = I_ESP_L3xy4z_D2x_a+ABZ*I_ESP_K3xy3z_D2x_a;
    Double I_ESP_K3x4z_F2xz_a = I_ESP_L3x5z_D2x_a+ABZ*I_ESP_K3x4z_D2x_a;
    Double I_ESP_K2x5y_F2xz_a = I_ESP_L2x5yz_D2x_a+ABZ*I_ESP_K2x5y_D2x_a;
    Double I_ESP_K2x4yz_F2xz_a = I_ESP_L2x4y2z_D2x_a+ABZ*I_ESP_K2x4yz_D2x_a;
    Double I_ESP_K2x3y2z_F2xz_a = I_ESP_L2x3y3z_D2x_a+ABZ*I_ESP_K2x3y2z_D2x_a;
    Double I_ESP_K2x2y3z_F2xz_a = I_ESP_L2x2y4z_D2x_a+ABZ*I_ESP_K2x2y3z_D2x_a;
    Double I_ESP_K2xy4z_F2xz_a = I_ESP_L2xy5z_D2x_a+ABZ*I_ESP_K2xy4z_D2x_a;
    Double I_ESP_K2x5z_F2xz_a = I_ESP_L2x6z_D2x_a+ABZ*I_ESP_K2x5z_D2x_a;
    Double I_ESP_Kx6y_F2xz_a = I_ESP_Lx6yz_D2x_a+ABZ*I_ESP_Kx6y_D2x_a;
    Double I_ESP_Kx5yz_F2xz_a = I_ESP_Lx5y2z_D2x_a+ABZ*I_ESP_Kx5yz_D2x_a;
    Double I_ESP_Kx4y2z_F2xz_a = I_ESP_Lx4y3z_D2x_a+ABZ*I_ESP_Kx4y2z_D2x_a;
    Double I_ESP_Kx3y3z_F2xz_a = I_ESP_Lx3y4z_D2x_a+ABZ*I_ESP_Kx3y3z_D2x_a;
    Double I_ESP_Kx2y4z_F2xz_a = I_ESP_Lx2y5z_D2x_a+ABZ*I_ESP_Kx2y4z_D2x_a;
    Double I_ESP_Kxy5z_F2xz_a = I_ESP_Lxy6z_D2x_a+ABZ*I_ESP_Kxy5z_D2x_a;
    Double I_ESP_Kx6z_F2xz_a = I_ESP_Lx7z_D2x_a+ABZ*I_ESP_Kx6z_D2x_a;
    Double I_ESP_K7y_F2xz_a = I_ESP_L7yz_D2x_a+ABZ*I_ESP_K7y_D2x_a;
    Double I_ESP_K6yz_F2xz_a = I_ESP_L6y2z_D2x_a+ABZ*I_ESP_K6yz_D2x_a;
    Double I_ESP_K5y2z_F2xz_a = I_ESP_L5y3z_D2x_a+ABZ*I_ESP_K5y2z_D2x_a;
    Double I_ESP_K4y3z_F2xz_a = I_ESP_L4y4z_D2x_a+ABZ*I_ESP_K4y3z_D2x_a;
    Double I_ESP_K3y4z_F2xz_a = I_ESP_L3y5z_D2x_a+ABZ*I_ESP_K3y4z_D2x_a;
    Double I_ESP_K2y5z_F2xz_a = I_ESP_L2y6z_D2x_a+ABZ*I_ESP_K2y5z_D2x_a;
    Double I_ESP_Ky6z_F2xz_a = I_ESP_Ly7z_D2x_a+ABZ*I_ESP_Ky6z_D2x_a;
    Double I_ESP_K7z_F2xz_a = I_ESP_L8z_D2x_a+ABZ*I_ESP_K7z_D2x_a;
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
    Double I_ESP_K7x_F2yz_a = I_ESP_L7xz_D2y_a+ABZ*I_ESP_K7x_D2y_a;
    Double I_ESP_K6xy_F2yz_a = I_ESP_L6xyz_D2y_a+ABZ*I_ESP_K6xy_D2y_a;
    Double I_ESP_K6xz_F2yz_a = I_ESP_L6x2z_D2y_a+ABZ*I_ESP_K6xz_D2y_a;
    Double I_ESP_K5x2y_F2yz_a = I_ESP_L5x2yz_D2y_a+ABZ*I_ESP_K5x2y_D2y_a;
    Double I_ESP_K5xyz_F2yz_a = I_ESP_L5xy2z_D2y_a+ABZ*I_ESP_K5xyz_D2y_a;
    Double I_ESP_K5x2z_F2yz_a = I_ESP_L5x3z_D2y_a+ABZ*I_ESP_K5x2z_D2y_a;
    Double I_ESP_K4x3y_F2yz_a = I_ESP_L4x3yz_D2y_a+ABZ*I_ESP_K4x3y_D2y_a;
    Double I_ESP_K4x2yz_F2yz_a = I_ESP_L4x2y2z_D2y_a+ABZ*I_ESP_K4x2yz_D2y_a;
    Double I_ESP_K4xy2z_F2yz_a = I_ESP_L4xy3z_D2y_a+ABZ*I_ESP_K4xy2z_D2y_a;
    Double I_ESP_K4x3z_F2yz_a = I_ESP_L4x4z_D2y_a+ABZ*I_ESP_K4x3z_D2y_a;
    Double I_ESP_K3x4y_F2yz_a = I_ESP_L3x4yz_D2y_a+ABZ*I_ESP_K3x4y_D2y_a;
    Double I_ESP_K3x3yz_F2yz_a = I_ESP_L3x3y2z_D2y_a+ABZ*I_ESP_K3x3yz_D2y_a;
    Double I_ESP_K3x2y2z_F2yz_a = I_ESP_L3x2y3z_D2y_a+ABZ*I_ESP_K3x2y2z_D2y_a;
    Double I_ESP_K3xy3z_F2yz_a = I_ESP_L3xy4z_D2y_a+ABZ*I_ESP_K3xy3z_D2y_a;
    Double I_ESP_K3x4z_F2yz_a = I_ESP_L3x5z_D2y_a+ABZ*I_ESP_K3x4z_D2y_a;
    Double I_ESP_K2x5y_F2yz_a = I_ESP_L2x5yz_D2y_a+ABZ*I_ESP_K2x5y_D2y_a;
    Double I_ESP_K2x4yz_F2yz_a = I_ESP_L2x4y2z_D2y_a+ABZ*I_ESP_K2x4yz_D2y_a;
    Double I_ESP_K2x3y2z_F2yz_a = I_ESP_L2x3y3z_D2y_a+ABZ*I_ESP_K2x3y2z_D2y_a;
    Double I_ESP_K2x2y3z_F2yz_a = I_ESP_L2x2y4z_D2y_a+ABZ*I_ESP_K2x2y3z_D2y_a;
    Double I_ESP_K2xy4z_F2yz_a = I_ESP_L2xy5z_D2y_a+ABZ*I_ESP_K2xy4z_D2y_a;
    Double I_ESP_K2x5z_F2yz_a = I_ESP_L2x6z_D2y_a+ABZ*I_ESP_K2x5z_D2y_a;
    Double I_ESP_Kx6y_F2yz_a = I_ESP_Lx6yz_D2y_a+ABZ*I_ESP_Kx6y_D2y_a;
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
     * totally 84 integrals are omitted 
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
     * shell quartet name: SQ_ESP_N_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 48 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_O_S_a
     * RHS shell quartet name: SQ_ESP_N_S_a
     ************************************************************/
    Double I_ESP_N10x_Px_a = I_ESP_O11x_S_a+ABX*I_ESP_N10x_S_a;
    Double I_ESP_N9xy_Px_a = I_ESP_O10xy_S_a+ABX*I_ESP_N9xy_S_a;
    Double I_ESP_N9xz_Px_a = I_ESP_O10xz_S_a+ABX*I_ESP_N9xz_S_a;
    Double I_ESP_N8x2y_Px_a = I_ESP_O9x2y_S_a+ABX*I_ESP_N8x2y_S_a;
    Double I_ESP_N8xyz_Px_a = I_ESP_O9xyz_S_a+ABX*I_ESP_N8xyz_S_a;
    Double I_ESP_N8x2z_Px_a = I_ESP_O9x2z_S_a+ABX*I_ESP_N8x2z_S_a;
    Double I_ESP_N7x3y_Px_a = I_ESP_O8x3y_S_a+ABX*I_ESP_N7x3y_S_a;
    Double I_ESP_N7x2yz_Px_a = I_ESP_O8x2yz_S_a+ABX*I_ESP_N7x2yz_S_a;
    Double I_ESP_N7xy2z_Px_a = I_ESP_O8xy2z_S_a+ABX*I_ESP_N7xy2z_S_a;
    Double I_ESP_N7x3z_Px_a = I_ESP_O8x3z_S_a+ABX*I_ESP_N7x3z_S_a;
    Double I_ESP_N6x4y_Px_a = I_ESP_O7x4y_S_a+ABX*I_ESP_N6x4y_S_a;
    Double I_ESP_N6x3yz_Px_a = I_ESP_O7x3yz_S_a+ABX*I_ESP_N6x3yz_S_a;
    Double I_ESP_N6x2y2z_Px_a = I_ESP_O7x2y2z_S_a+ABX*I_ESP_N6x2y2z_S_a;
    Double I_ESP_N6xy3z_Px_a = I_ESP_O7xy3z_S_a+ABX*I_ESP_N6xy3z_S_a;
    Double I_ESP_N6x4z_Px_a = I_ESP_O7x4z_S_a+ABX*I_ESP_N6x4z_S_a;
    Double I_ESP_N5x5y_Px_a = I_ESP_O6x5y_S_a+ABX*I_ESP_N5x5y_S_a;
    Double I_ESP_N5x4yz_Px_a = I_ESP_O6x4yz_S_a+ABX*I_ESP_N5x4yz_S_a;
    Double I_ESP_N5x3y2z_Px_a = I_ESP_O6x3y2z_S_a+ABX*I_ESP_N5x3y2z_S_a;
    Double I_ESP_N5x2y3z_Px_a = I_ESP_O6x2y3z_S_a+ABX*I_ESP_N5x2y3z_S_a;
    Double I_ESP_N5xy4z_Px_a = I_ESP_O6xy4z_S_a+ABX*I_ESP_N5xy4z_S_a;
    Double I_ESP_N5x5z_Px_a = I_ESP_O6x5z_S_a+ABX*I_ESP_N5x5z_S_a;
    Double I_ESP_N4x6y_Px_a = I_ESP_O5x6y_S_a+ABX*I_ESP_N4x6y_S_a;
    Double I_ESP_N4x5yz_Px_a = I_ESP_O5x5yz_S_a+ABX*I_ESP_N4x5yz_S_a;
    Double I_ESP_N4x4y2z_Px_a = I_ESP_O5x4y2z_S_a+ABX*I_ESP_N4x4y2z_S_a;
    Double I_ESP_N4x3y3z_Px_a = I_ESP_O5x3y3z_S_a+ABX*I_ESP_N4x3y3z_S_a;
    Double I_ESP_N4x2y4z_Px_a = I_ESP_O5x2y4z_S_a+ABX*I_ESP_N4x2y4z_S_a;
    Double I_ESP_N4xy5z_Px_a = I_ESP_O5xy5z_S_a+ABX*I_ESP_N4xy5z_S_a;
    Double I_ESP_N4x6z_Px_a = I_ESP_O5x6z_S_a+ABX*I_ESP_N4x6z_S_a;
    Double I_ESP_N3x7y_Px_a = I_ESP_O4x7y_S_a+ABX*I_ESP_N3x7y_S_a;
    Double I_ESP_N3x6yz_Px_a = I_ESP_O4x6yz_S_a+ABX*I_ESP_N3x6yz_S_a;
    Double I_ESP_N3x5y2z_Px_a = I_ESP_O4x5y2z_S_a+ABX*I_ESP_N3x5y2z_S_a;
    Double I_ESP_N3x4y3z_Px_a = I_ESP_O4x4y3z_S_a+ABX*I_ESP_N3x4y3z_S_a;
    Double I_ESP_N3x3y4z_Px_a = I_ESP_O4x3y4z_S_a+ABX*I_ESP_N3x3y4z_S_a;
    Double I_ESP_N3x2y5z_Px_a = I_ESP_O4x2y5z_S_a+ABX*I_ESP_N3x2y5z_S_a;
    Double I_ESP_N3xy6z_Px_a = I_ESP_O4xy6z_S_a+ABX*I_ESP_N3xy6z_S_a;
    Double I_ESP_N3x7z_Px_a = I_ESP_O4x7z_S_a+ABX*I_ESP_N3x7z_S_a;
    Double I_ESP_N2x8y_Px_a = I_ESP_O3x8y_S_a+ABX*I_ESP_N2x8y_S_a;
    Double I_ESP_N2x7yz_Px_a = I_ESP_O3x7yz_S_a+ABX*I_ESP_N2x7yz_S_a;
    Double I_ESP_N2x6y2z_Px_a = I_ESP_O3x6y2z_S_a+ABX*I_ESP_N2x6y2z_S_a;
    Double I_ESP_N2x5y3z_Px_a = I_ESP_O3x5y3z_S_a+ABX*I_ESP_N2x5y3z_S_a;
    Double I_ESP_N2x4y4z_Px_a = I_ESP_O3x4y4z_S_a+ABX*I_ESP_N2x4y4z_S_a;
    Double I_ESP_N2x3y5z_Px_a = I_ESP_O3x3y5z_S_a+ABX*I_ESP_N2x3y5z_S_a;
    Double I_ESP_N2x2y6z_Px_a = I_ESP_O3x2y6z_S_a+ABX*I_ESP_N2x2y6z_S_a;
    Double I_ESP_N2xy7z_Px_a = I_ESP_O3xy7z_S_a+ABX*I_ESP_N2xy7z_S_a;
    Double I_ESP_N2x8z_Px_a = I_ESP_O3x8z_S_a+ABX*I_ESP_N2x8z_S_a;
    Double I_ESP_Nx8yz_Px_a = I_ESP_O2x8yz_S_a+ABX*I_ESP_Nx8yz_S_a;
    Double I_ESP_Nx7y2z_Px_a = I_ESP_O2x7y2z_S_a+ABX*I_ESP_Nx7y2z_S_a;
    Double I_ESP_Nx6y3z_Px_a = I_ESP_O2x6y3z_S_a+ABX*I_ESP_Nx6y3z_S_a;
    Double I_ESP_Nx5y4z_Px_a = I_ESP_O2x5y4z_S_a+ABX*I_ESP_Nx5y4z_S_a;
    Double I_ESP_Nx4y5z_Px_a = I_ESP_O2x4y5z_S_a+ABX*I_ESP_Nx4y5z_S_a;
    Double I_ESP_Nx3y6z_Px_a = I_ESP_O2x3y6z_S_a+ABX*I_ESP_Nx3y6z_S_a;
    Double I_ESP_Nx2y7z_Px_a = I_ESP_O2x2y7z_S_a+ABX*I_ESP_Nx2y7z_S_a;
    Double I_ESP_Nxy8z_Px_a = I_ESP_O2xy8z_S_a+ABX*I_ESP_Nxy8z_S_a;
    Double I_ESP_N8x2y_Py_a = I_ESP_O8x3y_S_a+ABY*I_ESP_N8x2y_S_a;
    Double I_ESP_N7x3y_Py_a = I_ESP_O7x4y_S_a+ABY*I_ESP_N7x3y_S_a;
    Double I_ESP_N7x2yz_Py_a = I_ESP_O7x3yz_S_a+ABY*I_ESP_N7x2yz_S_a;
    Double I_ESP_N7xy2z_Py_a = I_ESP_O7x2y2z_S_a+ABY*I_ESP_N7xy2z_S_a;
    Double I_ESP_N6x4y_Py_a = I_ESP_O6x5y_S_a+ABY*I_ESP_N6x4y_S_a;
    Double I_ESP_N6x3yz_Py_a = I_ESP_O6x4yz_S_a+ABY*I_ESP_N6x3yz_S_a;
    Double I_ESP_N6x2y2z_Py_a = I_ESP_O6x3y2z_S_a+ABY*I_ESP_N6x2y2z_S_a;
    Double I_ESP_N6xy3z_Py_a = I_ESP_O6x2y3z_S_a+ABY*I_ESP_N6xy3z_S_a;
    Double I_ESP_N5x5y_Py_a = I_ESP_O5x6y_S_a+ABY*I_ESP_N5x5y_S_a;
    Double I_ESP_N5x4yz_Py_a = I_ESP_O5x5yz_S_a+ABY*I_ESP_N5x4yz_S_a;
    Double I_ESP_N5x3y2z_Py_a = I_ESP_O5x4y2z_S_a+ABY*I_ESP_N5x3y2z_S_a;
    Double I_ESP_N5x2y3z_Py_a = I_ESP_O5x3y3z_S_a+ABY*I_ESP_N5x2y3z_S_a;
    Double I_ESP_N5xy4z_Py_a = I_ESP_O5x2y4z_S_a+ABY*I_ESP_N5xy4z_S_a;
    Double I_ESP_N4x6y_Py_a = I_ESP_O4x7y_S_a+ABY*I_ESP_N4x6y_S_a;
    Double I_ESP_N4x5yz_Py_a = I_ESP_O4x6yz_S_a+ABY*I_ESP_N4x5yz_S_a;
    Double I_ESP_N4x4y2z_Py_a = I_ESP_O4x5y2z_S_a+ABY*I_ESP_N4x4y2z_S_a;
    Double I_ESP_N4x3y3z_Py_a = I_ESP_O4x4y3z_S_a+ABY*I_ESP_N4x3y3z_S_a;
    Double I_ESP_N4x2y4z_Py_a = I_ESP_O4x3y4z_S_a+ABY*I_ESP_N4x2y4z_S_a;
    Double I_ESP_N4xy5z_Py_a = I_ESP_O4x2y5z_S_a+ABY*I_ESP_N4xy5z_S_a;
    Double I_ESP_N3x7y_Py_a = I_ESP_O3x8y_S_a+ABY*I_ESP_N3x7y_S_a;
    Double I_ESP_N3x6yz_Py_a = I_ESP_O3x7yz_S_a+ABY*I_ESP_N3x6yz_S_a;
    Double I_ESP_N3x5y2z_Py_a = I_ESP_O3x6y2z_S_a+ABY*I_ESP_N3x5y2z_S_a;
    Double I_ESP_N3x4y3z_Py_a = I_ESP_O3x5y3z_S_a+ABY*I_ESP_N3x4y3z_S_a;
    Double I_ESP_N3x3y4z_Py_a = I_ESP_O3x4y4z_S_a+ABY*I_ESP_N3x3y4z_S_a;
    Double I_ESP_N3x2y5z_Py_a = I_ESP_O3x3y5z_S_a+ABY*I_ESP_N3x2y5z_S_a;
    Double I_ESP_N3xy6z_Py_a = I_ESP_O3x2y6z_S_a+ABY*I_ESP_N3xy6z_S_a;
    Double I_ESP_N2x8y_Py_a = I_ESP_O2x9y_S_a+ABY*I_ESP_N2x8y_S_a;
    Double I_ESP_N2x7yz_Py_a = I_ESP_O2x8yz_S_a+ABY*I_ESP_N2x7yz_S_a;
    Double I_ESP_N2x6y2z_Py_a = I_ESP_O2x7y2z_S_a+ABY*I_ESP_N2x6y2z_S_a;
    Double I_ESP_N2x5y3z_Py_a = I_ESP_O2x6y3z_S_a+ABY*I_ESP_N2x5y3z_S_a;
    Double I_ESP_N2x4y4z_Py_a = I_ESP_O2x5y4z_S_a+ABY*I_ESP_N2x4y4z_S_a;
    Double I_ESP_N2x3y5z_Py_a = I_ESP_O2x4y5z_S_a+ABY*I_ESP_N2x3y5z_S_a;
    Double I_ESP_N2x2y6z_Py_a = I_ESP_O2x3y6z_S_a+ABY*I_ESP_N2x2y6z_S_a;
    Double I_ESP_N2xy7z_Py_a = I_ESP_O2x2y7z_S_a+ABY*I_ESP_N2xy7z_S_a;
    Double I_ESP_Nx9y_Py_a = I_ESP_Ox10y_S_a+ABY*I_ESP_Nx9y_S_a;
    Double I_ESP_Nx8yz_Py_a = I_ESP_Ox9yz_S_a+ABY*I_ESP_Nx8yz_S_a;
    Double I_ESP_Nx7y2z_Py_a = I_ESP_Ox8y2z_S_a+ABY*I_ESP_Nx7y2z_S_a;
    Double I_ESP_Nx6y3z_Py_a = I_ESP_Ox7y3z_S_a+ABY*I_ESP_Nx6y3z_S_a;
    Double I_ESP_Nx5y4z_Py_a = I_ESP_Ox6y4z_S_a+ABY*I_ESP_Nx5y4z_S_a;
    Double I_ESP_Nx4y5z_Py_a = I_ESP_Ox5y5z_S_a+ABY*I_ESP_Nx4y5z_S_a;
    Double I_ESP_Nx3y6z_Py_a = I_ESP_Ox4y6z_S_a+ABY*I_ESP_Nx3y6z_S_a;
    Double I_ESP_Nx2y7z_Py_a = I_ESP_Ox3y7z_S_a+ABY*I_ESP_Nx2y7z_S_a;
    Double I_ESP_Nxy8z_Py_a = I_ESP_Ox2y8z_S_a+ABY*I_ESP_Nxy8z_S_a;
    Double I_ESP_N10y_Py_a = I_ESP_O11y_S_a+ABY*I_ESP_N10y_S_a;
    Double I_ESP_N9yz_Py_a = I_ESP_O10yz_S_a+ABY*I_ESP_N9yz_S_a;
    Double I_ESP_N8y2z_Py_a = I_ESP_O9y2z_S_a+ABY*I_ESP_N8y2z_S_a;
    Double I_ESP_N7y3z_Py_a = I_ESP_O8y3z_S_a+ABY*I_ESP_N7y3z_S_a;
    Double I_ESP_N6y4z_Py_a = I_ESP_O7y4z_S_a+ABY*I_ESP_N6y4z_S_a;
    Double I_ESP_N5y5z_Py_a = I_ESP_O6y5z_S_a+ABY*I_ESP_N5y5z_S_a;
    Double I_ESP_N4y6z_Py_a = I_ESP_O5y6z_S_a+ABY*I_ESP_N4y6z_S_a;
    Double I_ESP_N3y7z_Py_a = I_ESP_O4y7z_S_a+ABY*I_ESP_N3y7z_S_a;
    Double I_ESP_N2y8z_Py_a = I_ESP_O3y8z_S_a+ABY*I_ESP_N2y8z_S_a;
    Double I_ESP_N8x2z_Pz_a = I_ESP_O8x3z_S_a+ABZ*I_ESP_N8x2z_S_a;
    Double I_ESP_N7xy2z_Pz_a = I_ESP_O7xy3z_S_a+ABZ*I_ESP_N7xy2z_S_a;
    Double I_ESP_N7x3z_Pz_a = I_ESP_O7x4z_S_a+ABZ*I_ESP_N7x3z_S_a;
    Double I_ESP_N6x2y2z_Pz_a = I_ESP_O6x2y3z_S_a+ABZ*I_ESP_N6x2y2z_S_a;
    Double I_ESP_N6xy3z_Pz_a = I_ESP_O6xy4z_S_a+ABZ*I_ESP_N6xy3z_S_a;
    Double I_ESP_N6x4z_Pz_a = I_ESP_O6x5z_S_a+ABZ*I_ESP_N6x4z_S_a;
    Double I_ESP_N5x3y2z_Pz_a = I_ESP_O5x3y3z_S_a+ABZ*I_ESP_N5x3y2z_S_a;
    Double I_ESP_N5x2y3z_Pz_a = I_ESP_O5x2y4z_S_a+ABZ*I_ESP_N5x2y3z_S_a;
    Double I_ESP_N5xy4z_Pz_a = I_ESP_O5xy5z_S_a+ABZ*I_ESP_N5xy4z_S_a;
    Double I_ESP_N5x5z_Pz_a = I_ESP_O5x6z_S_a+ABZ*I_ESP_N5x5z_S_a;
    Double I_ESP_N4x4y2z_Pz_a = I_ESP_O4x4y3z_S_a+ABZ*I_ESP_N4x4y2z_S_a;
    Double I_ESP_N4x3y3z_Pz_a = I_ESP_O4x3y4z_S_a+ABZ*I_ESP_N4x3y3z_S_a;
    Double I_ESP_N4x2y4z_Pz_a = I_ESP_O4x2y5z_S_a+ABZ*I_ESP_N4x2y4z_S_a;
    Double I_ESP_N4xy5z_Pz_a = I_ESP_O4xy6z_S_a+ABZ*I_ESP_N4xy5z_S_a;
    Double I_ESP_N4x6z_Pz_a = I_ESP_O4x7z_S_a+ABZ*I_ESP_N4x6z_S_a;
    Double I_ESP_N3x5y2z_Pz_a = I_ESP_O3x5y3z_S_a+ABZ*I_ESP_N3x5y2z_S_a;
    Double I_ESP_N3x4y3z_Pz_a = I_ESP_O3x4y4z_S_a+ABZ*I_ESP_N3x4y3z_S_a;
    Double I_ESP_N3x3y4z_Pz_a = I_ESP_O3x3y5z_S_a+ABZ*I_ESP_N3x3y4z_S_a;
    Double I_ESP_N3x2y5z_Pz_a = I_ESP_O3x2y6z_S_a+ABZ*I_ESP_N3x2y5z_S_a;
    Double I_ESP_N3xy6z_Pz_a = I_ESP_O3xy7z_S_a+ABZ*I_ESP_N3xy6z_S_a;
    Double I_ESP_N3x7z_Pz_a = I_ESP_O3x8z_S_a+ABZ*I_ESP_N3x7z_S_a;
    Double I_ESP_N2x6y2z_Pz_a = I_ESP_O2x6y3z_S_a+ABZ*I_ESP_N2x6y2z_S_a;
    Double I_ESP_N2x5y3z_Pz_a = I_ESP_O2x5y4z_S_a+ABZ*I_ESP_N2x5y3z_S_a;
    Double I_ESP_N2x4y4z_Pz_a = I_ESP_O2x4y5z_S_a+ABZ*I_ESP_N2x4y4z_S_a;
    Double I_ESP_N2x3y5z_Pz_a = I_ESP_O2x3y6z_S_a+ABZ*I_ESP_N2x3y5z_S_a;
    Double I_ESP_N2x2y6z_Pz_a = I_ESP_O2x2y7z_S_a+ABZ*I_ESP_N2x2y6z_S_a;
    Double I_ESP_N2xy7z_Pz_a = I_ESP_O2xy8z_S_a+ABZ*I_ESP_N2xy7z_S_a;
    Double I_ESP_N2x8z_Pz_a = I_ESP_O2x9z_S_a+ABZ*I_ESP_N2x8z_S_a;
    Double I_ESP_Nx7y2z_Pz_a = I_ESP_Ox7y3z_S_a+ABZ*I_ESP_Nx7y2z_S_a;
    Double I_ESP_Nx6y3z_Pz_a = I_ESP_Ox6y4z_S_a+ABZ*I_ESP_Nx6y3z_S_a;
    Double I_ESP_Nx5y4z_Pz_a = I_ESP_Ox5y5z_S_a+ABZ*I_ESP_Nx5y4z_S_a;
    Double I_ESP_Nx4y5z_Pz_a = I_ESP_Ox4y6z_S_a+ABZ*I_ESP_Nx4y5z_S_a;
    Double I_ESP_Nx3y6z_Pz_a = I_ESP_Ox3y7z_S_a+ABZ*I_ESP_Nx3y6z_S_a;
    Double I_ESP_Nx2y7z_Pz_a = I_ESP_Ox2y8z_S_a+ABZ*I_ESP_Nx2y7z_S_a;
    Double I_ESP_Nxy8z_Pz_a = I_ESP_Oxy9z_S_a+ABZ*I_ESP_Nxy8z_S_a;
    Double I_ESP_Nx9z_Pz_a = I_ESP_Ox10z_S_a+ABZ*I_ESP_Nx9z_S_a;
    Double I_ESP_N8y2z_Pz_a = I_ESP_O8y3z_S_a+ABZ*I_ESP_N8y2z_S_a;
    Double I_ESP_N7y3z_Pz_a = I_ESP_O7y4z_S_a+ABZ*I_ESP_N7y3z_S_a;
    Double I_ESP_N6y4z_Pz_a = I_ESP_O6y5z_S_a+ABZ*I_ESP_N6y4z_S_a;
    Double I_ESP_N5y5z_Pz_a = I_ESP_O5y6z_S_a+ABZ*I_ESP_N5y5z_S_a;
    Double I_ESP_N4y6z_Pz_a = I_ESP_O4y7z_S_a+ABZ*I_ESP_N4y6z_S_a;
    Double I_ESP_N3y7z_Pz_a = I_ESP_O3y8z_S_a+ABZ*I_ESP_N3y7z_S_a;
    Double I_ESP_N2y8z_Pz_a = I_ESP_O2y9z_S_a+ABZ*I_ESP_N2y8z_S_a;
    Double I_ESP_Ny9z_Pz_a = I_ESP_Oy10z_S_a+ABZ*I_ESP_Ny9z_S_a;
    Double I_ESP_N10z_Pz_a = I_ESP_O11z_S_a+ABZ*I_ESP_N10z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_M_D_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 180 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_N_P_a
     * RHS shell quartet name: SQ_ESP_M_P_a
     ************************************************************/
    Double I_ESP_M9x_D2x_a = I_ESP_N10x_Px_a+ABX*I_ESP_M9x_Px_a;
    Double I_ESP_M8xy_D2x_a = I_ESP_N9xy_Px_a+ABX*I_ESP_M8xy_Px_a;
    Double I_ESP_M8xz_D2x_a = I_ESP_N9xz_Px_a+ABX*I_ESP_M8xz_Px_a;
    Double I_ESP_M7x2y_D2x_a = I_ESP_N8x2y_Px_a+ABX*I_ESP_M7x2y_Px_a;
    Double I_ESP_M7xyz_D2x_a = I_ESP_N8xyz_Px_a+ABX*I_ESP_M7xyz_Px_a;
    Double I_ESP_M7x2z_D2x_a = I_ESP_N8x2z_Px_a+ABX*I_ESP_M7x2z_Px_a;
    Double I_ESP_M6x3y_D2x_a = I_ESP_N7x3y_Px_a+ABX*I_ESP_M6x3y_Px_a;
    Double I_ESP_M6x2yz_D2x_a = I_ESP_N7x2yz_Px_a+ABX*I_ESP_M6x2yz_Px_a;
    Double I_ESP_M6xy2z_D2x_a = I_ESP_N7xy2z_Px_a+ABX*I_ESP_M6xy2z_Px_a;
    Double I_ESP_M6x3z_D2x_a = I_ESP_N7x3z_Px_a+ABX*I_ESP_M6x3z_Px_a;
    Double I_ESP_M5x4y_D2x_a = I_ESP_N6x4y_Px_a+ABX*I_ESP_M5x4y_Px_a;
    Double I_ESP_M5x3yz_D2x_a = I_ESP_N6x3yz_Px_a+ABX*I_ESP_M5x3yz_Px_a;
    Double I_ESP_M5x2y2z_D2x_a = I_ESP_N6x2y2z_Px_a+ABX*I_ESP_M5x2y2z_Px_a;
    Double I_ESP_M5xy3z_D2x_a = I_ESP_N6xy3z_Px_a+ABX*I_ESP_M5xy3z_Px_a;
    Double I_ESP_M5x4z_D2x_a = I_ESP_N6x4z_Px_a+ABX*I_ESP_M5x4z_Px_a;
    Double I_ESP_M4x5y_D2x_a = I_ESP_N5x5y_Px_a+ABX*I_ESP_M4x5y_Px_a;
    Double I_ESP_M4x4yz_D2x_a = I_ESP_N5x4yz_Px_a+ABX*I_ESP_M4x4yz_Px_a;
    Double I_ESP_M4x3y2z_D2x_a = I_ESP_N5x3y2z_Px_a+ABX*I_ESP_M4x3y2z_Px_a;
    Double I_ESP_M4x2y3z_D2x_a = I_ESP_N5x2y3z_Px_a+ABX*I_ESP_M4x2y3z_Px_a;
    Double I_ESP_M4xy4z_D2x_a = I_ESP_N5xy4z_Px_a+ABX*I_ESP_M4xy4z_Px_a;
    Double I_ESP_M4x5z_D2x_a = I_ESP_N5x5z_Px_a+ABX*I_ESP_M4x5z_Px_a;
    Double I_ESP_M3x6y_D2x_a = I_ESP_N4x6y_Px_a+ABX*I_ESP_M3x6y_Px_a;
    Double I_ESP_M3x5yz_D2x_a = I_ESP_N4x5yz_Px_a+ABX*I_ESP_M3x5yz_Px_a;
    Double I_ESP_M3x4y2z_D2x_a = I_ESP_N4x4y2z_Px_a+ABX*I_ESP_M3x4y2z_Px_a;
    Double I_ESP_M3x3y3z_D2x_a = I_ESP_N4x3y3z_Px_a+ABX*I_ESP_M3x3y3z_Px_a;
    Double I_ESP_M3x2y4z_D2x_a = I_ESP_N4x2y4z_Px_a+ABX*I_ESP_M3x2y4z_Px_a;
    Double I_ESP_M3xy5z_D2x_a = I_ESP_N4xy5z_Px_a+ABX*I_ESP_M3xy5z_Px_a;
    Double I_ESP_M3x6z_D2x_a = I_ESP_N4x6z_Px_a+ABX*I_ESP_M3x6z_Px_a;
    Double I_ESP_M2x7y_D2x_a = I_ESP_N3x7y_Px_a+ABX*I_ESP_M2x7y_Px_a;
    Double I_ESP_M2x6yz_D2x_a = I_ESP_N3x6yz_Px_a+ABX*I_ESP_M2x6yz_Px_a;
    Double I_ESP_M2x5y2z_D2x_a = I_ESP_N3x5y2z_Px_a+ABX*I_ESP_M2x5y2z_Px_a;
    Double I_ESP_M2x4y3z_D2x_a = I_ESP_N3x4y3z_Px_a+ABX*I_ESP_M2x4y3z_Px_a;
    Double I_ESP_M2x3y4z_D2x_a = I_ESP_N3x3y4z_Px_a+ABX*I_ESP_M2x3y4z_Px_a;
    Double I_ESP_M2x2y5z_D2x_a = I_ESP_N3x2y5z_Px_a+ABX*I_ESP_M2x2y5z_Px_a;
    Double I_ESP_M2xy6z_D2x_a = I_ESP_N3xy6z_Px_a+ABX*I_ESP_M2xy6z_Px_a;
    Double I_ESP_M2x7z_D2x_a = I_ESP_N3x7z_Px_a+ABX*I_ESP_M2x7z_Px_a;
    Double I_ESP_Mx8y_D2x_a = I_ESP_N2x8y_Px_a+ABX*I_ESP_Mx8y_Px_a;
    Double I_ESP_Mx7yz_D2x_a = I_ESP_N2x7yz_Px_a+ABX*I_ESP_Mx7yz_Px_a;
    Double I_ESP_Mx6y2z_D2x_a = I_ESP_N2x6y2z_Px_a+ABX*I_ESP_Mx6y2z_Px_a;
    Double I_ESP_Mx5y3z_D2x_a = I_ESP_N2x5y3z_Px_a+ABX*I_ESP_Mx5y3z_Px_a;
    Double I_ESP_Mx4y4z_D2x_a = I_ESP_N2x4y4z_Px_a+ABX*I_ESP_Mx4y4z_Px_a;
    Double I_ESP_Mx3y5z_D2x_a = I_ESP_N2x3y5z_Px_a+ABX*I_ESP_Mx3y5z_Px_a;
    Double I_ESP_Mx2y6z_D2x_a = I_ESP_N2x2y6z_Px_a+ABX*I_ESP_Mx2y6z_Px_a;
    Double I_ESP_Mxy7z_D2x_a = I_ESP_N2xy7z_Px_a+ABX*I_ESP_Mxy7z_Px_a;
    Double I_ESP_Mx8z_D2x_a = I_ESP_N2x8z_Px_a+ABX*I_ESP_Mx8z_Px_a;
    Double I_ESP_M8yz_D2x_a = I_ESP_Nx8yz_Px_a+ABX*I_ESP_M8yz_Px_a;
    Double I_ESP_M7y2z_D2x_a = I_ESP_Nx7y2z_Px_a+ABX*I_ESP_M7y2z_Px_a;
    Double I_ESP_M6y3z_D2x_a = I_ESP_Nx6y3z_Px_a+ABX*I_ESP_M6y3z_Px_a;
    Double I_ESP_M5y4z_D2x_a = I_ESP_Nx5y4z_Px_a+ABX*I_ESP_M5y4z_Px_a;
    Double I_ESP_M4y5z_D2x_a = I_ESP_Nx4y5z_Px_a+ABX*I_ESP_M4y5z_Px_a;
    Double I_ESP_M3y6z_D2x_a = I_ESP_Nx3y6z_Px_a+ABX*I_ESP_M3y6z_Px_a;
    Double I_ESP_M2y7z_D2x_a = I_ESP_Nx2y7z_Px_a+ABX*I_ESP_M2y7z_Px_a;
    Double I_ESP_My8z_D2x_a = I_ESP_Nxy8z_Px_a+ABX*I_ESP_My8z_Px_a;
    Double I_ESP_M8xy_D2y_a = I_ESP_N8x2y_Py_a+ABY*I_ESP_M8xy_Py_a;
    Double I_ESP_M7x2y_D2y_a = I_ESP_N7x3y_Py_a+ABY*I_ESP_M7x2y_Py_a;
    Double I_ESP_M7xyz_D2y_a = I_ESP_N7x2yz_Py_a+ABY*I_ESP_M7xyz_Py_a;
    Double I_ESP_M7x2z_D2y_a = I_ESP_N7xy2z_Py_a+ABY*I_ESP_M7x2z_Py_a;
    Double I_ESP_M6x3y_D2y_a = I_ESP_N6x4y_Py_a+ABY*I_ESP_M6x3y_Py_a;
    Double I_ESP_M6x2yz_D2y_a = I_ESP_N6x3yz_Py_a+ABY*I_ESP_M6x2yz_Py_a;
    Double I_ESP_M6xy2z_D2y_a = I_ESP_N6x2y2z_Py_a+ABY*I_ESP_M6xy2z_Py_a;
    Double I_ESP_M6x3z_D2y_a = I_ESP_N6xy3z_Py_a+ABY*I_ESP_M6x3z_Py_a;
    Double I_ESP_M5x4y_D2y_a = I_ESP_N5x5y_Py_a+ABY*I_ESP_M5x4y_Py_a;
    Double I_ESP_M5x3yz_D2y_a = I_ESP_N5x4yz_Py_a+ABY*I_ESP_M5x3yz_Py_a;
    Double I_ESP_M5x2y2z_D2y_a = I_ESP_N5x3y2z_Py_a+ABY*I_ESP_M5x2y2z_Py_a;
    Double I_ESP_M5xy3z_D2y_a = I_ESP_N5x2y3z_Py_a+ABY*I_ESP_M5xy3z_Py_a;
    Double I_ESP_M5x4z_D2y_a = I_ESP_N5xy4z_Py_a+ABY*I_ESP_M5x4z_Py_a;
    Double I_ESP_M4x5y_D2y_a = I_ESP_N4x6y_Py_a+ABY*I_ESP_M4x5y_Py_a;
    Double I_ESP_M4x4yz_D2y_a = I_ESP_N4x5yz_Py_a+ABY*I_ESP_M4x4yz_Py_a;
    Double I_ESP_M4x3y2z_D2y_a = I_ESP_N4x4y2z_Py_a+ABY*I_ESP_M4x3y2z_Py_a;
    Double I_ESP_M4x2y3z_D2y_a = I_ESP_N4x3y3z_Py_a+ABY*I_ESP_M4x2y3z_Py_a;
    Double I_ESP_M4xy4z_D2y_a = I_ESP_N4x2y4z_Py_a+ABY*I_ESP_M4xy4z_Py_a;
    Double I_ESP_M4x5z_D2y_a = I_ESP_N4xy5z_Py_a+ABY*I_ESP_M4x5z_Py_a;
    Double I_ESP_M3x6y_D2y_a = I_ESP_N3x7y_Py_a+ABY*I_ESP_M3x6y_Py_a;
    Double I_ESP_M3x5yz_D2y_a = I_ESP_N3x6yz_Py_a+ABY*I_ESP_M3x5yz_Py_a;
    Double I_ESP_M3x4y2z_D2y_a = I_ESP_N3x5y2z_Py_a+ABY*I_ESP_M3x4y2z_Py_a;
    Double I_ESP_M3x3y3z_D2y_a = I_ESP_N3x4y3z_Py_a+ABY*I_ESP_M3x3y3z_Py_a;
    Double I_ESP_M3x2y4z_D2y_a = I_ESP_N3x3y4z_Py_a+ABY*I_ESP_M3x2y4z_Py_a;
    Double I_ESP_M3xy5z_D2y_a = I_ESP_N3x2y5z_Py_a+ABY*I_ESP_M3xy5z_Py_a;
    Double I_ESP_M3x6z_D2y_a = I_ESP_N3xy6z_Py_a+ABY*I_ESP_M3x6z_Py_a;
    Double I_ESP_M2x7y_D2y_a = I_ESP_N2x8y_Py_a+ABY*I_ESP_M2x7y_Py_a;
    Double I_ESP_M2x6yz_D2y_a = I_ESP_N2x7yz_Py_a+ABY*I_ESP_M2x6yz_Py_a;
    Double I_ESP_M2x5y2z_D2y_a = I_ESP_N2x6y2z_Py_a+ABY*I_ESP_M2x5y2z_Py_a;
    Double I_ESP_M2x4y3z_D2y_a = I_ESP_N2x5y3z_Py_a+ABY*I_ESP_M2x4y3z_Py_a;
    Double I_ESP_M2x3y4z_D2y_a = I_ESP_N2x4y4z_Py_a+ABY*I_ESP_M2x3y4z_Py_a;
    Double I_ESP_M2x2y5z_D2y_a = I_ESP_N2x3y5z_Py_a+ABY*I_ESP_M2x2y5z_Py_a;
    Double I_ESP_M2xy6z_D2y_a = I_ESP_N2x2y6z_Py_a+ABY*I_ESP_M2xy6z_Py_a;
    Double I_ESP_M2x7z_D2y_a = I_ESP_N2xy7z_Py_a+ABY*I_ESP_M2x7z_Py_a;
    Double I_ESP_Mx8y_D2y_a = I_ESP_Nx9y_Py_a+ABY*I_ESP_Mx8y_Py_a;
    Double I_ESP_Mx7yz_D2y_a = I_ESP_Nx8yz_Py_a+ABY*I_ESP_Mx7yz_Py_a;
    Double I_ESP_Mx6y2z_D2y_a = I_ESP_Nx7y2z_Py_a+ABY*I_ESP_Mx6y2z_Py_a;
    Double I_ESP_Mx5y3z_D2y_a = I_ESP_Nx6y3z_Py_a+ABY*I_ESP_Mx5y3z_Py_a;
    Double I_ESP_Mx4y4z_D2y_a = I_ESP_Nx5y4z_Py_a+ABY*I_ESP_Mx4y4z_Py_a;
    Double I_ESP_Mx3y5z_D2y_a = I_ESP_Nx4y5z_Py_a+ABY*I_ESP_Mx3y5z_Py_a;
    Double I_ESP_Mx2y6z_D2y_a = I_ESP_Nx3y6z_Py_a+ABY*I_ESP_Mx2y6z_Py_a;
    Double I_ESP_Mxy7z_D2y_a = I_ESP_Nx2y7z_Py_a+ABY*I_ESP_Mxy7z_Py_a;
    Double I_ESP_Mx8z_D2y_a = I_ESP_Nxy8z_Py_a+ABY*I_ESP_Mx8z_Py_a;
    Double I_ESP_M9y_D2y_a = I_ESP_N10y_Py_a+ABY*I_ESP_M9y_Py_a;
    Double I_ESP_M8yz_D2y_a = I_ESP_N9yz_Py_a+ABY*I_ESP_M8yz_Py_a;
    Double I_ESP_M7y2z_D2y_a = I_ESP_N8y2z_Py_a+ABY*I_ESP_M7y2z_Py_a;
    Double I_ESP_M6y3z_D2y_a = I_ESP_N7y3z_Py_a+ABY*I_ESP_M6y3z_Py_a;
    Double I_ESP_M5y4z_D2y_a = I_ESP_N6y4z_Py_a+ABY*I_ESP_M5y4z_Py_a;
    Double I_ESP_M4y5z_D2y_a = I_ESP_N5y5z_Py_a+ABY*I_ESP_M4y5z_Py_a;
    Double I_ESP_M3y6z_D2y_a = I_ESP_N4y6z_Py_a+ABY*I_ESP_M3y6z_Py_a;
    Double I_ESP_M2y7z_D2y_a = I_ESP_N3y7z_Py_a+ABY*I_ESP_M2y7z_Py_a;
    Double I_ESP_My8z_D2y_a = I_ESP_N2y8z_Py_a+ABY*I_ESP_My8z_Py_a;
    Double I_ESP_M8xz_D2z_a = I_ESP_N8x2z_Pz_a+ABZ*I_ESP_M8xz_Pz_a;
    Double I_ESP_M7xyz_D2z_a = I_ESP_N7xy2z_Pz_a+ABZ*I_ESP_M7xyz_Pz_a;
    Double I_ESP_M7x2z_D2z_a = I_ESP_N7x3z_Pz_a+ABZ*I_ESP_M7x2z_Pz_a;
    Double I_ESP_M6x2yz_D2z_a = I_ESP_N6x2y2z_Pz_a+ABZ*I_ESP_M6x2yz_Pz_a;
    Double I_ESP_M6xy2z_D2z_a = I_ESP_N6xy3z_Pz_a+ABZ*I_ESP_M6xy2z_Pz_a;
    Double I_ESP_M6x3z_D2z_a = I_ESP_N6x4z_Pz_a+ABZ*I_ESP_M6x3z_Pz_a;
    Double I_ESP_M5x3yz_D2z_a = I_ESP_N5x3y2z_Pz_a+ABZ*I_ESP_M5x3yz_Pz_a;
    Double I_ESP_M5x2y2z_D2z_a = I_ESP_N5x2y3z_Pz_a+ABZ*I_ESP_M5x2y2z_Pz_a;
    Double I_ESP_M5xy3z_D2z_a = I_ESP_N5xy4z_Pz_a+ABZ*I_ESP_M5xy3z_Pz_a;
    Double I_ESP_M5x4z_D2z_a = I_ESP_N5x5z_Pz_a+ABZ*I_ESP_M5x4z_Pz_a;
    Double I_ESP_M4x4yz_D2z_a = I_ESP_N4x4y2z_Pz_a+ABZ*I_ESP_M4x4yz_Pz_a;
    Double I_ESP_M4x3y2z_D2z_a = I_ESP_N4x3y3z_Pz_a+ABZ*I_ESP_M4x3y2z_Pz_a;
    Double I_ESP_M4x2y3z_D2z_a = I_ESP_N4x2y4z_Pz_a+ABZ*I_ESP_M4x2y3z_Pz_a;
    Double I_ESP_M4xy4z_D2z_a = I_ESP_N4xy5z_Pz_a+ABZ*I_ESP_M4xy4z_Pz_a;
    Double I_ESP_M4x5z_D2z_a = I_ESP_N4x6z_Pz_a+ABZ*I_ESP_M4x5z_Pz_a;
    Double I_ESP_M3x5yz_D2z_a = I_ESP_N3x5y2z_Pz_a+ABZ*I_ESP_M3x5yz_Pz_a;
    Double I_ESP_M3x4y2z_D2z_a = I_ESP_N3x4y3z_Pz_a+ABZ*I_ESP_M3x4y2z_Pz_a;
    Double I_ESP_M3x3y3z_D2z_a = I_ESP_N3x3y4z_Pz_a+ABZ*I_ESP_M3x3y3z_Pz_a;
    Double I_ESP_M3x2y4z_D2z_a = I_ESP_N3x2y5z_Pz_a+ABZ*I_ESP_M3x2y4z_Pz_a;
    Double I_ESP_M3xy5z_D2z_a = I_ESP_N3xy6z_Pz_a+ABZ*I_ESP_M3xy5z_Pz_a;
    Double I_ESP_M3x6z_D2z_a = I_ESP_N3x7z_Pz_a+ABZ*I_ESP_M3x6z_Pz_a;
    Double I_ESP_M2x6yz_D2z_a = I_ESP_N2x6y2z_Pz_a+ABZ*I_ESP_M2x6yz_Pz_a;
    Double I_ESP_M2x5y2z_D2z_a = I_ESP_N2x5y3z_Pz_a+ABZ*I_ESP_M2x5y2z_Pz_a;
    Double I_ESP_M2x4y3z_D2z_a = I_ESP_N2x4y4z_Pz_a+ABZ*I_ESP_M2x4y3z_Pz_a;
    Double I_ESP_M2x3y4z_D2z_a = I_ESP_N2x3y5z_Pz_a+ABZ*I_ESP_M2x3y4z_Pz_a;
    Double I_ESP_M2x2y5z_D2z_a = I_ESP_N2x2y6z_Pz_a+ABZ*I_ESP_M2x2y5z_Pz_a;
    Double I_ESP_M2xy6z_D2z_a = I_ESP_N2xy7z_Pz_a+ABZ*I_ESP_M2xy6z_Pz_a;
    Double I_ESP_M2x7z_D2z_a = I_ESP_N2x8z_Pz_a+ABZ*I_ESP_M2x7z_Pz_a;
    Double I_ESP_Mx7yz_D2z_a = I_ESP_Nx7y2z_Pz_a+ABZ*I_ESP_Mx7yz_Pz_a;
    Double I_ESP_Mx6y2z_D2z_a = I_ESP_Nx6y3z_Pz_a+ABZ*I_ESP_Mx6y2z_Pz_a;
    Double I_ESP_Mx5y3z_D2z_a = I_ESP_Nx5y4z_Pz_a+ABZ*I_ESP_Mx5y3z_Pz_a;
    Double I_ESP_Mx4y4z_D2z_a = I_ESP_Nx4y5z_Pz_a+ABZ*I_ESP_Mx4y4z_Pz_a;
    Double I_ESP_Mx3y5z_D2z_a = I_ESP_Nx3y6z_Pz_a+ABZ*I_ESP_Mx3y5z_Pz_a;
    Double I_ESP_Mx2y6z_D2z_a = I_ESP_Nx2y7z_Pz_a+ABZ*I_ESP_Mx2y6z_Pz_a;
    Double I_ESP_Mxy7z_D2z_a = I_ESP_Nxy8z_Pz_a+ABZ*I_ESP_Mxy7z_Pz_a;
    Double I_ESP_Mx8z_D2z_a = I_ESP_Nx9z_Pz_a+ABZ*I_ESP_Mx8z_Pz_a;
    Double I_ESP_M8yz_D2z_a = I_ESP_N8y2z_Pz_a+ABZ*I_ESP_M8yz_Pz_a;
    Double I_ESP_M7y2z_D2z_a = I_ESP_N7y3z_Pz_a+ABZ*I_ESP_M7y2z_Pz_a;
    Double I_ESP_M6y3z_D2z_a = I_ESP_N6y4z_Pz_a+ABZ*I_ESP_M6y3z_Pz_a;
    Double I_ESP_M5y4z_D2z_a = I_ESP_N5y5z_Pz_a+ABZ*I_ESP_M5y4z_Pz_a;
    Double I_ESP_M4y5z_D2z_a = I_ESP_N4y6z_Pz_a+ABZ*I_ESP_M4y5z_Pz_a;
    Double I_ESP_M3y6z_D2z_a = I_ESP_N3y7z_Pz_a+ABZ*I_ESP_M3y6z_Pz_a;
    Double I_ESP_M2y7z_D2z_a = I_ESP_N2y8z_Pz_a+ABZ*I_ESP_M2y7z_Pz_a;
    Double I_ESP_My8z_D2z_a = I_ESP_Ny9z_Pz_a+ABZ*I_ESP_My8z_Pz_a;
    Double I_ESP_M9z_D2z_a = I_ESP_N10z_Pz_a+ABZ*I_ESP_M9z_Pz_a;

    /************************************************************
     * shell quartet name: SQ_ESP_L_F_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 231 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_M_D_a
     * RHS shell quartet name: SQ_ESP_L_D_a
     ************************************************************/
    Double I_ESP_L8x_F3x_a = I_ESP_M9x_D2x_a+ABX*I_ESP_L8x_D2x_a;
    Double I_ESP_L7xy_F3x_a = I_ESP_M8xy_D2x_a+ABX*I_ESP_L7xy_D2x_a;
    Double I_ESP_L7xz_F3x_a = I_ESP_M8xz_D2x_a+ABX*I_ESP_L7xz_D2x_a;
    Double I_ESP_L6x2y_F3x_a = I_ESP_M7x2y_D2x_a+ABX*I_ESP_L6x2y_D2x_a;
    Double I_ESP_L6xyz_F3x_a = I_ESP_M7xyz_D2x_a+ABX*I_ESP_L6xyz_D2x_a;
    Double I_ESP_L6x2z_F3x_a = I_ESP_M7x2z_D2x_a+ABX*I_ESP_L6x2z_D2x_a;
    Double I_ESP_L5x3y_F3x_a = I_ESP_M6x3y_D2x_a+ABX*I_ESP_L5x3y_D2x_a;
    Double I_ESP_L5x2yz_F3x_a = I_ESP_M6x2yz_D2x_a+ABX*I_ESP_L5x2yz_D2x_a;
    Double I_ESP_L5xy2z_F3x_a = I_ESP_M6xy2z_D2x_a+ABX*I_ESP_L5xy2z_D2x_a;
    Double I_ESP_L5x3z_F3x_a = I_ESP_M6x3z_D2x_a+ABX*I_ESP_L5x3z_D2x_a;
    Double I_ESP_L4x4y_F3x_a = I_ESP_M5x4y_D2x_a+ABX*I_ESP_L4x4y_D2x_a;
    Double I_ESP_L4x3yz_F3x_a = I_ESP_M5x3yz_D2x_a+ABX*I_ESP_L4x3yz_D2x_a;
    Double I_ESP_L4x2y2z_F3x_a = I_ESP_M5x2y2z_D2x_a+ABX*I_ESP_L4x2y2z_D2x_a;
    Double I_ESP_L4xy3z_F3x_a = I_ESP_M5xy3z_D2x_a+ABX*I_ESP_L4xy3z_D2x_a;
    Double I_ESP_L4x4z_F3x_a = I_ESP_M5x4z_D2x_a+ABX*I_ESP_L4x4z_D2x_a;
    Double I_ESP_L3x5y_F3x_a = I_ESP_M4x5y_D2x_a+ABX*I_ESP_L3x5y_D2x_a;
    Double I_ESP_L3x4yz_F3x_a = I_ESP_M4x4yz_D2x_a+ABX*I_ESP_L3x4yz_D2x_a;
    Double I_ESP_L3x3y2z_F3x_a = I_ESP_M4x3y2z_D2x_a+ABX*I_ESP_L3x3y2z_D2x_a;
    Double I_ESP_L3x2y3z_F3x_a = I_ESP_M4x2y3z_D2x_a+ABX*I_ESP_L3x2y3z_D2x_a;
    Double I_ESP_L3xy4z_F3x_a = I_ESP_M4xy4z_D2x_a+ABX*I_ESP_L3xy4z_D2x_a;
    Double I_ESP_L3x5z_F3x_a = I_ESP_M4x5z_D2x_a+ABX*I_ESP_L3x5z_D2x_a;
    Double I_ESP_L2x6y_F3x_a = I_ESP_M3x6y_D2x_a+ABX*I_ESP_L2x6y_D2x_a;
    Double I_ESP_L2x5yz_F3x_a = I_ESP_M3x5yz_D2x_a+ABX*I_ESP_L2x5yz_D2x_a;
    Double I_ESP_L2x4y2z_F3x_a = I_ESP_M3x4y2z_D2x_a+ABX*I_ESP_L2x4y2z_D2x_a;
    Double I_ESP_L2x3y3z_F3x_a = I_ESP_M3x3y3z_D2x_a+ABX*I_ESP_L2x3y3z_D2x_a;
    Double I_ESP_L2x2y4z_F3x_a = I_ESP_M3x2y4z_D2x_a+ABX*I_ESP_L2x2y4z_D2x_a;
    Double I_ESP_L2xy5z_F3x_a = I_ESP_M3xy5z_D2x_a+ABX*I_ESP_L2xy5z_D2x_a;
    Double I_ESP_L2x6z_F3x_a = I_ESP_M3x6z_D2x_a+ABX*I_ESP_L2x6z_D2x_a;
    Double I_ESP_Lx7y_F3x_a = I_ESP_M2x7y_D2x_a+ABX*I_ESP_Lx7y_D2x_a;
    Double I_ESP_Lx6yz_F3x_a = I_ESP_M2x6yz_D2x_a+ABX*I_ESP_Lx6yz_D2x_a;
    Double I_ESP_Lx5y2z_F3x_a = I_ESP_M2x5y2z_D2x_a+ABX*I_ESP_Lx5y2z_D2x_a;
    Double I_ESP_Lx4y3z_F3x_a = I_ESP_M2x4y3z_D2x_a+ABX*I_ESP_Lx4y3z_D2x_a;
    Double I_ESP_Lx3y4z_F3x_a = I_ESP_M2x3y4z_D2x_a+ABX*I_ESP_Lx3y4z_D2x_a;
    Double I_ESP_Lx2y5z_F3x_a = I_ESP_M2x2y5z_D2x_a+ABX*I_ESP_Lx2y5z_D2x_a;
    Double I_ESP_Lxy6z_F3x_a = I_ESP_M2xy6z_D2x_a+ABX*I_ESP_Lxy6z_D2x_a;
    Double I_ESP_Lx7z_F3x_a = I_ESP_M2x7z_D2x_a+ABX*I_ESP_Lx7z_D2x_a;
    Double I_ESP_L8y_F3x_a = I_ESP_Mx8y_D2x_a+ABX*I_ESP_L8y_D2x_a;
    Double I_ESP_L7yz_F3x_a = I_ESP_Mx7yz_D2x_a+ABX*I_ESP_L7yz_D2x_a;
    Double I_ESP_L6y2z_F3x_a = I_ESP_Mx6y2z_D2x_a+ABX*I_ESP_L6y2z_D2x_a;
    Double I_ESP_L5y3z_F3x_a = I_ESP_Mx5y3z_D2x_a+ABX*I_ESP_L5y3z_D2x_a;
    Double I_ESP_L4y4z_F3x_a = I_ESP_Mx4y4z_D2x_a+ABX*I_ESP_L4y4z_D2x_a;
    Double I_ESP_L3y5z_F3x_a = I_ESP_Mx3y5z_D2x_a+ABX*I_ESP_L3y5z_D2x_a;
    Double I_ESP_L2y6z_F3x_a = I_ESP_Mx2y6z_D2x_a+ABX*I_ESP_L2y6z_D2x_a;
    Double I_ESP_Ly7z_F3x_a = I_ESP_Mxy7z_D2x_a+ABX*I_ESP_Ly7z_D2x_a;
    Double I_ESP_L8z_F3x_a = I_ESP_Mx8z_D2x_a+ABX*I_ESP_L8z_D2x_a;
    Double I_ESP_L6xyz_F2xy_a = I_ESP_M6x2yz_D2x_a+ABY*I_ESP_L6xyz_D2x_a;
    Double I_ESP_L5x2yz_F2xy_a = I_ESP_M5x3yz_D2x_a+ABY*I_ESP_L5x2yz_D2x_a;
    Double I_ESP_L5xy2z_F2xy_a = I_ESP_M5x2y2z_D2x_a+ABY*I_ESP_L5xy2z_D2x_a;
    Double I_ESP_L4x3yz_F2xy_a = I_ESP_M4x4yz_D2x_a+ABY*I_ESP_L4x3yz_D2x_a;
    Double I_ESP_L4x2y2z_F2xy_a = I_ESP_M4x3y2z_D2x_a+ABY*I_ESP_L4x2y2z_D2x_a;
    Double I_ESP_L4xy3z_F2xy_a = I_ESP_M4x2y3z_D2x_a+ABY*I_ESP_L4xy3z_D2x_a;
    Double I_ESP_L3x4yz_F2xy_a = I_ESP_M3x5yz_D2x_a+ABY*I_ESP_L3x4yz_D2x_a;
    Double I_ESP_L3x3y2z_F2xy_a = I_ESP_M3x4y2z_D2x_a+ABY*I_ESP_L3x3y2z_D2x_a;
    Double I_ESP_L3x2y3z_F2xy_a = I_ESP_M3x3y3z_D2x_a+ABY*I_ESP_L3x2y3z_D2x_a;
    Double I_ESP_L3xy4z_F2xy_a = I_ESP_M3x2y4z_D2x_a+ABY*I_ESP_L3xy4z_D2x_a;
    Double I_ESP_L2x5yz_F2xy_a = I_ESP_M2x6yz_D2x_a+ABY*I_ESP_L2x5yz_D2x_a;
    Double I_ESP_L2x4y2z_F2xy_a = I_ESP_M2x5y2z_D2x_a+ABY*I_ESP_L2x4y2z_D2x_a;
    Double I_ESP_L2x3y3z_F2xy_a = I_ESP_M2x4y3z_D2x_a+ABY*I_ESP_L2x3y3z_D2x_a;
    Double I_ESP_L2x2y4z_F2xy_a = I_ESP_M2x3y4z_D2x_a+ABY*I_ESP_L2x2y4z_D2x_a;
    Double I_ESP_L2xy5z_F2xy_a = I_ESP_M2x2y5z_D2x_a+ABY*I_ESP_L2xy5z_D2x_a;
    Double I_ESP_Lx6yz_F2xy_a = I_ESP_Mx7yz_D2x_a+ABY*I_ESP_Lx6yz_D2x_a;
    Double I_ESP_Lx5y2z_F2xy_a = I_ESP_Mx6y2z_D2x_a+ABY*I_ESP_Lx5y2z_D2x_a;
    Double I_ESP_Lx4y3z_F2xy_a = I_ESP_Mx5y3z_D2x_a+ABY*I_ESP_Lx4y3z_D2x_a;
    Double I_ESP_Lx3y4z_F2xy_a = I_ESP_Mx4y4z_D2x_a+ABY*I_ESP_Lx3y4z_D2x_a;
    Double I_ESP_Lx2y5z_F2xy_a = I_ESP_Mx3y5z_D2x_a+ABY*I_ESP_Lx2y5z_D2x_a;
    Double I_ESP_Lxy6z_F2xy_a = I_ESP_Mx2y6z_D2x_a+ABY*I_ESP_Lxy6z_D2x_a;
    Double I_ESP_L7yz_F2xy_a = I_ESP_M8yz_D2x_a+ABY*I_ESP_L7yz_D2x_a;
    Double I_ESP_L6y2z_F2xy_a = I_ESP_M7y2z_D2x_a+ABY*I_ESP_L6y2z_D2x_a;
    Double I_ESP_L5y3z_F2xy_a = I_ESP_M6y3z_D2x_a+ABY*I_ESP_L5y3z_D2x_a;
    Double I_ESP_L4y4z_F2xy_a = I_ESP_M5y4z_D2x_a+ABY*I_ESP_L4y4z_D2x_a;
    Double I_ESP_L3y5z_F2xy_a = I_ESP_M4y5z_D2x_a+ABY*I_ESP_L3y5z_D2x_a;
    Double I_ESP_L2y6z_F2xy_a = I_ESP_M3y6z_D2x_a+ABY*I_ESP_L2y6z_D2x_a;
    Double I_ESP_Ly7z_F2xy_a = I_ESP_M2y7z_D2x_a+ABY*I_ESP_Ly7z_D2x_a;
    Double I_ESP_L6xyz_F2xz_a = I_ESP_M6xy2z_D2x_a+ABZ*I_ESP_L6xyz_D2x_a;
    Double I_ESP_L5x2yz_F2xz_a = I_ESP_M5x2y2z_D2x_a+ABZ*I_ESP_L5x2yz_D2x_a;
    Double I_ESP_L5xy2z_F2xz_a = I_ESP_M5xy3z_D2x_a+ABZ*I_ESP_L5xy2z_D2x_a;
    Double I_ESP_L4x3yz_F2xz_a = I_ESP_M4x3y2z_D2x_a+ABZ*I_ESP_L4x3yz_D2x_a;
    Double I_ESP_L4x2y2z_F2xz_a = I_ESP_M4x2y3z_D2x_a+ABZ*I_ESP_L4x2y2z_D2x_a;
    Double I_ESP_L4xy3z_F2xz_a = I_ESP_M4xy4z_D2x_a+ABZ*I_ESP_L4xy3z_D2x_a;
    Double I_ESP_L3x4yz_F2xz_a = I_ESP_M3x4y2z_D2x_a+ABZ*I_ESP_L3x4yz_D2x_a;
    Double I_ESP_L3x3y2z_F2xz_a = I_ESP_M3x3y3z_D2x_a+ABZ*I_ESP_L3x3y2z_D2x_a;
    Double I_ESP_L3x2y3z_F2xz_a = I_ESP_M3x2y4z_D2x_a+ABZ*I_ESP_L3x2y3z_D2x_a;
    Double I_ESP_L3xy4z_F2xz_a = I_ESP_M3xy5z_D2x_a+ABZ*I_ESP_L3xy4z_D2x_a;
    Double I_ESP_L2x5yz_F2xz_a = I_ESP_M2x5y2z_D2x_a+ABZ*I_ESP_L2x5yz_D2x_a;
    Double I_ESP_L2x4y2z_F2xz_a = I_ESP_M2x4y3z_D2x_a+ABZ*I_ESP_L2x4y2z_D2x_a;
    Double I_ESP_L2x3y3z_F2xz_a = I_ESP_M2x3y4z_D2x_a+ABZ*I_ESP_L2x3y3z_D2x_a;
    Double I_ESP_L2x2y4z_F2xz_a = I_ESP_M2x2y5z_D2x_a+ABZ*I_ESP_L2x2y4z_D2x_a;
    Double I_ESP_L2xy5z_F2xz_a = I_ESP_M2xy6z_D2x_a+ABZ*I_ESP_L2xy5z_D2x_a;
    Double I_ESP_Lx6yz_F2xz_a = I_ESP_Mx6y2z_D2x_a+ABZ*I_ESP_Lx6yz_D2x_a;
    Double I_ESP_Lx5y2z_F2xz_a = I_ESP_Mx5y3z_D2x_a+ABZ*I_ESP_Lx5y2z_D2x_a;
    Double I_ESP_Lx4y3z_F2xz_a = I_ESP_Mx4y4z_D2x_a+ABZ*I_ESP_Lx4y3z_D2x_a;
    Double I_ESP_Lx3y4z_F2xz_a = I_ESP_Mx3y5z_D2x_a+ABZ*I_ESP_Lx3y4z_D2x_a;
    Double I_ESP_Lx2y5z_F2xz_a = I_ESP_Mx2y6z_D2x_a+ABZ*I_ESP_Lx2y5z_D2x_a;
    Double I_ESP_Lxy6z_F2xz_a = I_ESP_Mxy7z_D2x_a+ABZ*I_ESP_Lxy6z_D2x_a;
    Double I_ESP_L7yz_F2xz_a = I_ESP_M7y2z_D2x_a+ABZ*I_ESP_L7yz_D2x_a;
    Double I_ESP_L6y2z_F2xz_a = I_ESP_M6y3z_D2x_a+ABZ*I_ESP_L6y2z_D2x_a;
    Double I_ESP_L5y3z_F2xz_a = I_ESP_M5y4z_D2x_a+ABZ*I_ESP_L5y3z_D2x_a;
    Double I_ESP_L4y4z_F2xz_a = I_ESP_M4y5z_D2x_a+ABZ*I_ESP_L4y4z_D2x_a;
    Double I_ESP_L3y5z_F2xz_a = I_ESP_M3y6z_D2x_a+ABZ*I_ESP_L3y5z_D2x_a;
    Double I_ESP_L2y6z_F2xz_a = I_ESP_M2y7z_D2x_a+ABZ*I_ESP_L2y6z_D2x_a;
    Double I_ESP_Ly7z_F2xz_a = I_ESP_My8z_D2x_a+ABZ*I_ESP_Ly7z_D2x_a;
    Double I_ESP_L8x_F3y_a = I_ESP_M8xy_D2y_a+ABY*I_ESP_L8x_D2y_a;
    Double I_ESP_L7xy_F3y_a = I_ESP_M7x2y_D2y_a+ABY*I_ESP_L7xy_D2y_a;
    Double I_ESP_L7xz_F3y_a = I_ESP_M7xyz_D2y_a+ABY*I_ESP_L7xz_D2y_a;
    Double I_ESP_L6x2y_F3y_a = I_ESP_M6x3y_D2y_a+ABY*I_ESP_L6x2y_D2y_a;
    Double I_ESP_L6xyz_F3y_a = I_ESP_M6x2yz_D2y_a+ABY*I_ESP_L6xyz_D2y_a;
    Double I_ESP_L6x2z_F3y_a = I_ESP_M6xy2z_D2y_a+ABY*I_ESP_L6x2z_D2y_a;
    Double I_ESP_L5x3y_F3y_a = I_ESP_M5x4y_D2y_a+ABY*I_ESP_L5x3y_D2y_a;
    Double I_ESP_L5x2yz_F3y_a = I_ESP_M5x3yz_D2y_a+ABY*I_ESP_L5x2yz_D2y_a;
    Double I_ESP_L5xy2z_F3y_a = I_ESP_M5x2y2z_D2y_a+ABY*I_ESP_L5xy2z_D2y_a;
    Double I_ESP_L5x3z_F3y_a = I_ESP_M5xy3z_D2y_a+ABY*I_ESP_L5x3z_D2y_a;
    Double I_ESP_L4x4y_F3y_a = I_ESP_M4x5y_D2y_a+ABY*I_ESP_L4x4y_D2y_a;
    Double I_ESP_L4x3yz_F3y_a = I_ESP_M4x4yz_D2y_a+ABY*I_ESP_L4x3yz_D2y_a;
    Double I_ESP_L4x2y2z_F3y_a = I_ESP_M4x3y2z_D2y_a+ABY*I_ESP_L4x2y2z_D2y_a;
    Double I_ESP_L4xy3z_F3y_a = I_ESP_M4x2y3z_D2y_a+ABY*I_ESP_L4xy3z_D2y_a;
    Double I_ESP_L4x4z_F3y_a = I_ESP_M4xy4z_D2y_a+ABY*I_ESP_L4x4z_D2y_a;
    Double I_ESP_L3x5y_F3y_a = I_ESP_M3x6y_D2y_a+ABY*I_ESP_L3x5y_D2y_a;
    Double I_ESP_L3x4yz_F3y_a = I_ESP_M3x5yz_D2y_a+ABY*I_ESP_L3x4yz_D2y_a;
    Double I_ESP_L3x3y2z_F3y_a = I_ESP_M3x4y2z_D2y_a+ABY*I_ESP_L3x3y2z_D2y_a;
    Double I_ESP_L3x2y3z_F3y_a = I_ESP_M3x3y3z_D2y_a+ABY*I_ESP_L3x2y3z_D2y_a;
    Double I_ESP_L3xy4z_F3y_a = I_ESP_M3x2y4z_D2y_a+ABY*I_ESP_L3xy4z_D2y_a;
    Double I_ESP_L3x5z_F3y_a = I_ESP_M3xy5z_D2y_a+ABY*I_ESP_L3x5z_D2y_a;
    Double I_ESP_L2x6y_F3y_a = I_ESP_M2x7y_D2y_a+ABY*I_ESP_L2x6y_D2y_a;
    Double I_ESP_L2x5yz_F3y_a = I_ESP_M2x6yz_D2y_a+ABY*I_ESP_L2x5yz_D2y_a;
    Double I_ESP_L2x4y2z_F3y_a = I_ESP_M2x5y2z_D2y_a+ABY*I_ESP_L2x4y2z_D2y_a;
    Double I_ESP_L2x3y3z_F3y_a = I_ESP_M2x4y3z_D2y_a+ABY*I_ESP_L2x3y3z_D2y_a;
    Double I_ESP_L2x2y4z_F3y_a = I_ESP_M2x3y4z_D2y_a+ABY*I_ESP_L2x2y4z_D2y_a;
    Double I_ESP_L2xy5z_F3y_a = I_ESP_M2x2y5z_D2y_a+ABY*I_ESP_L2xy5z_D2y_a;
    Double I_ESP_L2x6z_F3y_a = I_ESP_M2xy6z_D2y_a+ABY*I_ESP_L2x6z_D2y_a;
    Double I_ESP_Lx7y_F3y_a = I_ESP_Mx8y_D2y_a+ABY*I_ESP_Lx7y_D2y_a;
    Double I_ESP_Lx6yz_F3y_a = I_ESP_Mx7yz_D2y_a+ABY*I_ESP_Lx6yz_D2y_a;
    Double I_ESP_Lx5y2z_F3y_a = I_ESP_Mx6y2z_D2y_a+ABY*I_ESP_Lx5y2z_D2y_a;
    Double I_ESP_Lx4y3z_F3y_a = I_ESP_Mx5y3z_D2y_a+ABY*I_ESP_Lx4y3z_D2y_a;
    Double I_ESP_Lx3y4z_F3y_a = I_ESP_Mx4y4z_D2y_a+ABY*I_ESP_Lx3y4z_D2y_a;
    Double I_ESP_Lx2y5z_F3y_a = I_ESP_Mx3y5z_D2y_a+ABY*I_ESP_Lx2y5z_D2y_a;
    Double I_ESP_Lxy6z_F3y_a = I_ESP_Mx2y6z_D2y_a+ABY*I_ESP_Lxy6z_D2y_a;
    Double I_ESP_Lx7z_F3y_a = I_ESP_Mxy7z_D2y_a+ABY*I_ESP_Lx7z_D2y_a;
    Double I_ESP_L8y_F3y_a = I_ESP_M9y_D2y_a+ABY*I_ESP_L8y_D2y_a;
    Double I_ESP_L7yz_F3y_a = I_ESP_M8yz_D2y_a+ABY*I_ESP_L7yz_D2y_a;
    Double I_ESP_L6y2z_F3y_a = I_ESP_M7y2z_D2y_a+ABY*I_ESP_L6y2z_D2y_a;
    Double I_ESP_L5y3z_F3y_a = I_ESP_M6y3z_D2y_a+ABY*I_ESP_L5y3z_D2y_a;
    Double I_ESP_L4y4z_F3y_a = I_ESP_M5y4z_D2y_a+ABY*I_ESP_L4y4z_D2y_a;
    Double I_ESP_L3y5z_F3y_a = I_ESP_M4y5z_D2y_a+ABY*I_ESP_L3y5z_D2y_a;
    Double I_ESP_L2y6z_F3y_a = I_ESP_M3y6z_D2y_a+ABY*I_ESP_L2y6z_D2y_a;
    Double I_ESP_Ly7z_F3y_a = I_ESP_M2y7z_D2y_a+ABY*I_ESP_Ly7z_D2y_a;
    Double I_ESP_L8z_F3y_a = I_ESP_My8z_D2y_a+ABY*I_ESP_L8z_D2y_a;
    Double I_ESP_L7xz_F2yz_a = I_ESP_M7x2z_D2y_a+ABZ*I_ESP_L7xz_D2y_a;
    Double I_ESP_L6xyz_F2yz_a = I_ESP_M6xy2z_D2y_a+ABZ*I_ESP_L6xyz_D2y_a;
    Double I_ESP_L6x2z_F2yz_a = I_ESP_M6x3z_D2y_a+ABZ*I_ESP_L6x2z_D2y_a;
    Double I_ESP_L5x2yz_F2yz_a = I_ESP_M5x2y2z_D2y_a+ABZ*I_ESP_L5x2yz_D2y_a;
    Double I_ESP_L5xy2z_F2yz_a = I_ESP_M5xy3z_D2y_a+ABZ*I_ESP_L5xy2z_D2y_a;
    Double I_ESP_L5x3z_F2yz_a = I_ESP_M5x4z_D2y_a+ABZ*I_ESP_L5x3z_D2y_a;
    Double I_ESP_L4x3yz_F2yz_a = I_ESP_M4x3y2z_D2y_a+ABZ*I_ESP_L4x3yz_D2y_a;
    Double I_ESP_L4x2y2z_F2yz_a = I_ESP_M4x2y3z_D2y_a+ABZ*I_ESP_L4x2y2z_D2y_a;
    Double I_ESP_L4xy3z_F2yz_a = I_ESP_M4xy4z_D2y_a+ABZ*I_ESP_L4xy3z_D2y_a;
    Double I_ESP_L4x4z_F2yz_a = I_ESP_M4x5z_D2y_a+ABZ*I_ESP_L4x4z_D2y_a;
    Double I_ESP_L3x4yz_F2yz_a = I_ESP_M3x4y2z_D2y_a+ABZ*I_ESP_L3x4yz_D2y_a;
    Double I_ESP_L3x3y2z_F2yz_a = I_ESP_M3x3y3z_D2y_a+ABZ*I_ESP_L3x3y2z_D2y_a;
    Double I_ESP_L3x2y3z_F2yz_a = I_ESP_M3x2y4z_D2y_a+ABZ*I_ESP_L3x2y3z_D2y_a;
    Double I_ESP_L3xy4z_F2yz_a = I_ESP_M3xy5z_D2y_a+ABZ*I_ESP_L3xy4z_D2y_a;
    Double I_ESP_L3x5z_F2yz_a = I_ESP_M3x6z_D2y_a+ABZ*I_ESP_L3x5z_D2y_a;
    Double I_ESP_L2x5yz_F2yz_a = I_ESP_M2x5y2z_D2y_a+ABZ*I_ESP_L2x5yz_D2y_a;
    Double I_ESP_L2x4y2z_F2yz_a = I_ESP_M2x4y3z_D2y_a+ABZ*I_ESP_L2x4y2z_D2y_a;
    Double I_ESP_L2x3y3z_F2yz_a = I_ESP_M2x3y4z_D2y_a+ABZ*I_ESP_L2x3y3z_D2y_a;
    Double I_ESP_L2x2y4z_F2yz_a = I_ESP_M2x2y5z_D2y_a+ABZ*I_ESP_L2x2y4z_D2y_a;
    Double I_ESP_L2xy5z_F2yz_a = I_ESP_M2xy6z_D2y_a+ABZ*I_ESP_L2xy5z_D2y_a;
    Double I_ESP_L2x6z_F2yz_a = I_ESP_M2x7z_D2y_a+ABZ*I_ESP_L2x6z_D2y_a;
    Double I_ESP_Lx6yz_F2yz_a = I_ESP_Mx6y2z_D2y_a+ABZ*I_ESP_Lx6yz_D2y_a;
    Double I_ESP_Lx5y2z_F2yz_a = I_ESP_Mx5y3z_D2y_a+ABZ*I_ESP_Lx5y2z_D2y_a;
    Double I_ESP_Lx4y3z_F2yz_a = I_ESP_Mx4y4z_D2y_a+ABZ*I_ESP_Lx4y3z_D2y_a;
    Double I_ESP_Lx3y4z_F2yz_a = I_ESP_Mx3y5z_D2y_a+ABZ*I_ESP_Lx3y4z_D2y_a;
    Double I_ESP_Lx2y5z_F2yz_a = I_ESP_Mx2y6z_D2y_a+ABZ*I_ESP_Lx2y5z_D2y_a;
    Double I_ESP_Lxy6z_F2yz_a = I_ESP_Mxy7z_D2y_a+ABZ*I_ESP_Lxy6z_D2y_a;
    Double I_ESP_Lx7z_F2yz_a = I_ESP_Mx8z_D2y_a+ABZ*I_ESP_Lx7z_D2y_a;
    Double I_ESP_L8x_F3z_a = I_ESP_M8xz_D2z_a+ABZ*I_ESP_L8x_D2z_a;
    Double I_ESP_L7xy_F3z_a = I_ESP_M7xyz_D2z_a+ABZ*I_ESP_L7xy_D2z_a;
    Double I_ESP_L7xz_F3z_a = I_ESP_M7x2z_D2z_a+ABZ*I_ESP_L7xz_D2z_a;
    Double I_ESP_L6x2y_F3z_a = I_ESP_M6x2yz_D2z_a+ABZ*I_ESP_L6x2y_D2z_a;
    Double I_ESP_L6xyz_F3z_a = I_ESP_M6xy2z_D2z_a+ABZ*I_ESP_L6xyz_D2z_a;
    Double I_ESP_L6x2z_F3z_a = I_ESP_M6x3z_D2z_a+ABZ*I_ESP_L6x2z_D2z_a;
    Double I_ESP_L5x3y_F3z_a = I_ESP_M5x3yz_D2z_a+ABZ*I_ESP_L5x3y_D2z_a;
    Double I_ESP_L5x2yz_F3z_a = I_ESP_M5x2y2z_D2z_a+ABZ*I_ESP_L5x2yz_D2z_a;
    Double I_ESP_L5xy2z_F3z_a = I_ESP_M5xy3z_D2z_a+ABZ*I_ESP_L5xy2z_D2z_a;
    Double I_ESP_L5x3z_F3z_a = I_ESP_M5x4z_D2z_a+ABZ*I_ESP_L5x3z_D2z_a;
    Double I_ESP_L4x4y_F3z_a = I_ESP_M4x4yz_D2z_a+ABZ*I_ESP_L4x4y_D2z_a;
    Double I_ESP_L4x3yz_F3z_a = I_ESP_M4x3y2z_D2z_a+ABZ*I_ESP_L4x3yz_D2z_a;
    Double I_ESP_L4x2y2z_F3z_a = I_ESP_M4x2y3z_D2z_a+ABZ*I_ESP_L4x2y2z_D2z_a;
    Double I_ESP_L4xy3z_F3z_a = I_ESP_M4xy4z_D2z_a+ABZ*I_ESP_L4xy3z_D2z_a;
    Double I_ESP_L4x4z_F3z_a = I_ESP_M4x5z_D2z_a+ABZ*I_ESP_L4x4z_D2z_a;
    Double I_ESP_L3x5y_F3z_a = I_ESP_M3x5yz_D2z_a+ABZ*I_ESP_L3x5y_D2z_a;
    Double I_ESP_L3x4yz_F3z_a = I_ESP_M3x4y2z_D2z_a+ABZ*I_ESP_L3x4yz_D2z_a;
    Double I_ESP_L3x3y2z_F3z_a = I_ESP_M3x3y3z_D2z_a+ABZ*I_ESP_L3x3y2z_D2z_a;
    Double I_ESP_L3x2y3z_F3z_a = I_ESP_M3x2y4z_D2z_a+ABZ*I_ESP_L3x2y3z_D2z_a;
    Double I_ESP_L3xy4z_F3z_a = I_ESP_M3xy5z_D2z_a+ABZ*I_ESP_L3xy4z_D2z_a;
    Double I_ESP_L3x5z_F3z_a = I_ESP_M3x6z_D2z_a+ABZ*I_ESP_L3x5z_D2z_a;
    Double I_ESP_L2x6y_F3z_a = I_ESP_M2x6yz_D2z_a+ABZ*I_ESP_L2x6y_D2z_a;
    Double I_ESP_L2x5yz_F3z_a = I_ESP_M2x5y2z_D2z_a+ABZ*I_ESP_L2x5yz_D2z_a;
    Double I_ESP_L2x4y2z_F3z_a = I_ESP_M2x4y3z_D2z_a+ABZ*I_ESP_L2x4y2z_D2z_a;
    Double I_ESP_L2x3y3z_F3z_a = I_ESP_M2x3y4z_D2z_a+ABZ*I_ESP_L2x3y3z_D2z_a;
    Double I_ESP_L2x2y4z_F3z_a = I_ESP_M2x2y5z_D2z_a+ABZ*I_ESP_L2x2y4z_D2z_a;
    Double I_ESP_L2xy5z_F3z_a = I_ESP_M2xy6z_D2z_a+ABZ*I_ESP_L2xy5z_D2z_a;
    Double I_ESP_L2x6z_F3z_a = I_ESP_M2x7z_D2z_a+ABZ*I_ESP_L2x6z_D2z_a;
    Double I_ESP_Lx7y_F3z_a = I_ESP_Mx7yz_D2z_a+ABZ*I_ESP_Lx7y_D2z_a;
    Double I_ESP_Lx6yz_F3z_a = I_ESP_Mx6y2z_D2z_a+ABZ*I_ESP_Lx6yz_D2z_a;
    Double I_ESP_Lx5y2z_F3z_a = I_ESP_Mx5y3z_D2z_a+ABZ*I_ESP_Lx5y2z_D2z_a;
    Double I_ESP_Lx4y3z_F3z_a = I_ESP_Mx4y4z_D2z_a+ABZ*I_ESP_Lx4y3z_D2z_a;
    Double I_ESP_Lx3y4z_F3z_a = I_ESP_Mx3y5z_D2z_a+ABZ*I_ESP_Lx3y4z_D2z_a;
    Double I_ESP_Lx2y5z_F3z_a = I_ESP_Mx2y6z_D2z_a+ABZ*I_ESP_Lx2y5z_D2z_a;
    Double I_ESP_Lxy6z_F3z_a = I_ESP_Mxy7z_D2z_a+ABZ*I_ESP_Lxy6z_D2z_a;
    Double I_ESP_Lx7z_F3z_a = I_ESP_Mx8z_D2z_a+ABZ*I_ESP_Lx7z_D2z_a;
    Double I_ESP_L8y_F3z_a = I_ESP_M8yz_D2z_a+ABZ*I_ESP_L8y_D2z_a;
    Double I_ESP_L7yz_F3z_a = I_ESP_M7y2z_D2z_a+ABZ*I_ESP_L7yz_D2z_a;
    Double I_ESP_L6y2z_F3z_a = I_ESP_M6y3z_D2z_a+ABZ*I_ESP_L6y2z_D2z_a;
    Double I_ESP_L5y3z_F3z_a = I_ESP_M5y4z_D2z_a+ABZ*I_ESP_L5y3z_D2z_a;
    Double I_ESP_L4y4z_F3z_a = I_ESP_M4y5z_D2z_a+ABZ*I_ESP_L4y4z_D2z_a;
    Double I_ESP_L3y5z_F3z_a = I_ESP_M3y6z_D2z_a+ABZ*I_ESP_L3y5z_D2z_a;
    Double I_ESP_L2y6z_F3z_a = I_ESP_M2y7z_D2z_a+ABZ*I_ESP_L2y6z_D2z_a;
    Double I_ESP_Ly7z_F3z_a = I_ESP_My8z_D2z_a+ABZ*I_ESP_Ly7z_D2z_a;
    Double I_ESP_L8z_F3z_a = I_ESP_M9z_D2z_a+ABZ*I_ESP_L8z_D2z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_K_G_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 159 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_F_a
     * RHS shell quartet name: SQ_ESP_K_F_a
     ************************************************************/
    Double I_ESP_K7x_G4x_a = I_ESP_L8x_F3x_a+ABX*I_ESP_K7x_F3x_a;
    Double I_ESP_K6xy_G4x_a = I_ESP_L7xy_F3x_a+ABX*I_ESP_K6xy_F3x_a;
    Double I_ESP_K6xz_G4x_a = I_ESP_L7xz_F3x_a+ABX*I_ESP_K6xz_F3x_a;
    Double I_ESP_K5x2y_G4x_a = I_ESP_L6x2y_F3x_a+ABX*I_ESP_K5x2y_F3x_a;
    Double I_ESP_K5xyz_G4x_a = I_ESP_L6xyz_F3x_a+ABX*I_ESP_K5xyz_F3x_a;
    Double I_ESP_K5x2z_G4x_a = I_ESP_L6x2z_F3x_a+ABX*I_ESP_K5x2z_F3x_a;
    Double I_ESP_K4x3y_G4x_a = I_ESP_L5x3y_F3x_a+ABX*I_ESP_K4x3y_F3x_a;
    Double I_ESP_K4x2yz_G4x_a = I_ESP_L5x2yz_F3x_a+ABX*I_ESP_K4x2yz_F3x_a;
    Double I_ESP_K4xy2z_G4x_a = I_ESP_L5xy2z_F3x_a+ABX*I_ESP_K4xy2z_F3x_a;
    Double I_ESP_K4x3z_G4x_a = I_ESP_L5x3z_F3x_a+ABX*I_ESP_K4x3z_F3x_a;
    Double I_ESP_K3x4y_G4x_a = I_ESP_L4x4y_F3x_a+ABX*I_ESP_K3x4y_F3x_a;
    Double I_ESP_K3x3yz_G4x_a = I_ESP_L4x3yz_F3x_a+ABX*I_ESP_K3x3yz_F3x_a;
    Double I_ESP_K3x2y2z_G4x_a = I_ESP_L4x2y2z_F3x_a+ABX*I_ESP_K3x2y2z_F3x_a;
    Double I_ESP_K3xy3z_G4x_a = I_ESP_L4xy3z_F3x_a+ABX*I_ESP_K3xy3z_F3x_a;
    Double I_ESP_K3x4z_G4x_a = I_ESP_L4x4z_F3x_a+ABX*I_ESP_K3x4z_F3x_a;
    Double I_ESP_K2x5y_G4x_a = I_ESP_L3x5y_F3x_a+ABX*I_ESP_K2x5y_F3x_a;
    Double I_ESP_K2x4yz_G4x_a = I_ESP_L3x4yz_F3x_a+ABX*I_ESP_K2x4yz_F3x_a;
    Double I_ESP_K2x3y2z_G4x_a = I_ESP_L3x3y2z_F3x_a+ABX*I_ESP_K2x3y2z_F3x_a;
    Double I_ESP_K2x2y3z_G4x_a = I_ESP_L3x2y3z_F3x_a+ABX*I_ESP_K2x2y3z_F3x_a;
    Double I_ESP_K2xy4z_G4x_a = I_ESP_L3xy4z_F3x_a+ABX*I_ESP_K2xy4z_F3x_a;
    Double I_ESP_K2x5z_G4x_a = I_ESP_L3x5z_F3x_a+ABX*I_ESP_K2x5z_F3x_a;
    Double I_ESP_Kx6y_G4x_a = I_ESP_L2x6y_F3x_a+ABX*I_ESP_Kx6y_F3x_a;
    Double I_ESP_Kx5yz_G4x_a = I_ESP_L2x5yz_F3x_a+ABX*I_ESP_Kx5yz_F3x_a;
    Double I_ESP_Kx4y2z_G4x_a = I_ESP_L2x4y2z_F3x_a+ABX*I_ESP_Kx4y2z_F3x_a;
    Double I_ESP_Kx3y3z_G4x_a = I_ESP_L2x3y3z_F3x_a+ABX*I_ESP_Kx3y3z_F3x_a;
    Double I_ESP_Kx2y4z_G4x_a = I_ESP_L2x2y4z_F3x_a+ABX*I_ESP_Kx2y4z_F3x_a;
    Double I_ESP_Kxy5z_G4x_a = I_ESP_L2xy5z_F3x_a+ABX*I_ESP_Kxy5z_F3x_a;
    Double I_ESP_Kx6z_G4x_a = I_ESP_L2x6z_F3x_a+ABX*I_ESP_Kx6z_F3x_a;
    Double I_ESP_K7y_G4x_a = I_ESP_Lx7y_F3x_a+ABX*I_ESP_K7y_F3x_a;
    Double I_ESP_K6yz_G4x_a = I_ESP_Lx6yz_F3x_a+ABX*I_ESP_K6yz_F3x_a;
    Double I_ESP_K5y2z_G4x_a = I_ESP_Lx5y2z_F3x_a+ABX*I_ESP_K5y2z_F3x_a;
    Double I_ESP_K4y3z_G4x_a = I_ESP_Lx4y3z_F3x_a+ABX*I_ESP_K4y3z_F3x_a;
    Double I_ESP_K3y4z_G4x_a = I_ESP_Lx3y4z_F3x_a+ABX*I_ESP_K3y4z_F3x_a;
    Double I_ESP_K2y5z_G4x_a = I_ESP_Lx2y5z_F3x_a+ABX*I_ESP_K2y5z_F3x_a;
    Double I_ESP_Ky6z_G4x_a = I_ESP_Lxy6z_F3x_a+ABX*I_ESP_Ky6z_F3x_a;
    Double I_ESP_K7z_G4x_a = I_ESP_Lx7z_F3x_a+ABX*I_ESP_K7z_F3x_a;
    Double I_ESP_K6xy_G3xy_a = I_ESP_L6x2y_F3x_a+ABY*I_ESP_K6xy_F3x_a;
    Double I_ESP_K6xz_G3xy_a = I_ESP_L6xyz_F3x_a+ABY*I_ESP_K6xz_F3x_a;
    Double I_ESP_K5x2y_G3xy_a = I_ESP_L5x3y_F3x_a+ABY*I_ESP_K5x2y_F3x_a;
    Double I_ESP_K5xyz_G3xy_a = I_ESP_L5x2yz_F3x_a+ABY*I_ESP_K5xyz_F3x_a;
    Double I_ESP_K5x2z_G3xy_a = I_ESP_L5xy2z_F3x_a+ABY*I_ESP_K5x2z_F3x_a;
    Double I_ESP_K4x3y_G3xy_a = I_ESP_L4x4y_F3x_a+ABY*I_ESP_K4x3y_F3x_a;
    Double I_ESP_K4x2yz_G3xy_a = I_ESP_L4x3yz_F3x_a+ABY*I_ESP_K4x2yz_F3x_a;
    Double I_ESP_K4xy2z_G3xy_a = I_ESP_L4x2y2z_F3x_a+ABY*I_ESP_K4xy2z_F3x_a;
    Double I_ESP_K4x3z_G3xy_a = I_ESP_L4xy3z_F3x_a+ABY*I_ESP_K4x3z_F3x_a;
    Double I_ESP_K3x4y_G3xy_a = I_ESP_L3x5y_F3x_a+ABY*I_ESP_K3x4y_F3x_a;
    Double I_ESP_K3x3yz_G3xy_a = I_ESP_L3x4yz_F3x_a+ABY*I_ESP_K3x3yz_F3x_a;
    Double I_ESP_K3x2y2z_G3xy_a = I_ESP_L3x3y2z_F3x_a+ABY*I_ESP_K3x2y2z_F3x_a;
    Double I_ESP_K3xy3z_G3xy_a = I_ESP_L3x2y3z_F3x_a+ABY*I_ESP_K3xy3z_F3x_a;
    Double I_ESP_K3x4z_G3xy_a = I_ESP_L3xy4z_F3x_a+ABY*I_ESP_K3x4z_F3x_a;
    Double I_ESP_K2x5y_G3xy_a = I_ESP_L2x6y_F3x_a+ABY*I_ESP_K2x5y_F3x_a;
    Double I_ESP_K2x4yz_G3xy_a = I_ESP_L2x5yz_F3x_a+ABY*I_ESP_K2x4yz_F3x_a;
    Double I_ESP_K2x3y2z_G3xy_a = I_ESP_L2x4y2z_F3x_a+ABY*I_ESP_K2x3y2z_F3x_a;
    Double I_ESP_K2x2y3z_G3xy_a = I_ESP_L2x3y3z_F3x_a+ABY*I_ESP_K2x2y3z_F3x_a;
    Double I_ESP_K2xy4z_G3xy_a = I_ESP_L2x2y4z_F3x_a+ABY*I_ESP_K2xy4z_F3x_a;
    Double I_ESP_K2x5z_G3xy_a = I_ESP_L2xy5z_F3x_a+ABY*I_ESP_K2x5z_F3x_a;
    Double I_ESP_Kx6y_G3xy_a = I_ESP_Lx7y_F3x_a+ABY*I_ESP_Kx6y_F3x_a;
    Double I_ESP_Kx5yz_G3xy_a = I_ESP_Lx6yz_F3x_a+ABY*I_ESP_Kx5yz_F3x_a;
    Double I_ESP_Kx4y2z_G3xy_a = I_ESP_Lx5y2z_F3x_a+ABY*I_ESP_Kx4y2z_F3x_a;
    Double I_ESP_Kx3y3z_G3xy_a = I_ESP_Lx4y3z_F3x_a+ABY*I_ESP_Kx3y3z_F3x_a;
    Double I_ESP_Kx2y4z_G3xy_a = I_ESP_Lx3y4z_F3x_a+ABY*I_ESP_Kx2y4z_F3x_a;
    Double I_ESP_Kxy5z_G3xy_a = I_ESP_Lx2y5z_F3x_a+ABY*I_ESP_Kxy5z_F3x_a;
    Double I_ESP_Kx6z_G3xy_a = I_ESP_Lxy6z_F3x_a+ABY*I_ESP_Kx6z_F3x_a;
    Double I_ESP_K7y_G3xy_a = I_ESP_L8y_F3x_a+ABY*I_ESP_K7y_F3x_a;
    Double I_ESP_K6yz_G3xy_a = I_ESP_L7yz_F3x_a+ABY*I_ESP_K6yz_F3x_a;
    Double I_ESP_K5y2z_G3xy_a = I_ESP_L6y2z_F3x_a+ABY*I_ESP_K5y2z_F3x_a;
    Double I_ESP_K4y3z_G3xy_a = I_ESP_L5y3z_F3x_a+ABY*I_ESP_K4y3z_F3x_a;
    Double I_ESP_K3y4z_G3xy_a = I_ESP_L4y4z_F3x_a+ABY*I_ESP_K3y4z_F3x_a;
    Double I_ESP_K2y5z_G3xy_a = I_ESP_L3y5z_F3x_a+ABY*I_ESP_K2y5z_F3x_a;
    Double I_ESP_Ky6z_G3xy_a = I_ESP_L2y6z_F3x_a+ABY*I_ESP_Ky6z_F3x_a;
    Double I_ESP_K7z_G3xy_a = I_ESP_Ly7z_F3x_a+ABY*I_ESP_K7z_F3x_a;
    Double I_ESP_K6xz_G3xz_a = I_ESP_L6x2z_F3x_a+ABZ*I_ESP_K6xz_F3x_a;
    Double I_ESP_K5xyz_G3xz_a = I_ESP_L5xy2z_F3x_a+ABZ*I_ESP_K5xyz_F3x_a;
    Double I_ESP_K5x2z_G3xz_a = I_ESP_L5x3z_F3x_a+ABZ*I_ESP_K5x2z_F3x_a;
    Double I_ESP_K4x2yz_G3xz_a = I_ESP_L4x2y2z_F3x_a+ABZ*I_ESP_K4x2yz_F3x_a;
    Double I_ESP_K4xy2z_G3xz_a = I_ESP_L4xy3z_F3x_a+ABZ*I_ESP_K4xy2z_F3x_a;
    Double I_ESP_K4x3z_G3xz_a = I_ESP_L4x4z_F3x_a+ABZ*I_ESP_K4x3z_F3x_a;
    Double I_ESP_K3x3yz_G3xz_a = I_ESP_L3x3y2z_F3x_a+ABZ*I_ESP_K3x3yz_F3x_a;
    Double I_ESP_K3x2y2z_G3xz_a = I_ESP_L3x2y3z_F3x_a+ABZ*I_ESP_K3x2y2z_F3x_a;
    Double I_ESP_K3xy3z_G3xz_a = I_ESP_L3xy4z_F3x_a+ABZ*I_ESP_K3xy3z_F3x_a;
    Double I_ESP_K3x4z_G3xz_a = I_ESP_L3x5z_F3x_a+ABZ*I_ESP_K3x4z_F3x_a;
    Double I_ESP_K2x4yz_G3xz_a = I_ESP_L2x4y2z_F3x_a+ABZ*I_ESP_K2x4yz_F3x_a;
    Double I_ESP_K2x3y2z_G3xz_a = I_ESP_L2x3y3z_F3x_a+ABZ*I_ESP_K2x3y2z_F3x_a;
    Double I_ESP_K2x2y3z_G3xz_a = I_ESP_L2x2y4z_F3x_a+ABZ*I_ESP_K2x2y3z_F3x_a;
    Double I_ESP_K2xy4z_G3xz_a = I_ESP_L2xy5z_F3x_a+ABZ*I_ESP_K2xy4z_F3x_a;
    Double I_ESP_K2x5z_G3xz_a = I_ESP_L2x6z_F3x_a+ABZ*I_ESP_K2x5z_F3x_a;
    Double I_ESP_Kx5yz_G3xz_a = I_ESP_Lx5y2z_F3x_a+ABZ*I_ESP_Kx5yz_F3x_a;
    Double I_ESP_Kx4y2z_G3xz_a = I_ESP_Lx4y3z_F3x_a+ABZ*I_ESP_Kx4y2z_F3x_a;
    Double I_ESP_Kx3y3z_G3xz_a = I_ESP_Lx3y4z_F3x_a+ABZ*I_ESP_Kx3y3z_F3x_a;
    Double I_ESP_Kx2y4z_G3xz_a = I_ESP_Lx2y5z_F3x_a+ABZ*I_ESP_Kx2y4z_F3x_a;
    Double I_ESP_Kxy5z_G3xz_a = I_ESP_Lxy6z_F3x_a+ABZ*I_ESP_Kxy5z_F3x_a;
    Double I_ESP_Kx6z_G3xz_a = I_ESP_Lx7z_F3x_a+ABZ*I_ESP_Kx6z_F3x_a;
    Double I_ESP_K6yz_G3xz_a = I_ESP_L6y2z_F3x_a+ABZ*I_ESP_K6yz_F3x_a;
    Double I_ESP_K5y2z_G3xz_a = I_ESP_L5y3z_F3x_a+ABZ*I_ESP_K5y2z_F3x_a;
    Double I_ESP_K4y3z_G3xz_a = I_ESP_L4y4z_F3x_a+ABZ*I_ESP_K4y3z_F3x_a;
    Double I_ESP_K3y4z_G3xz_a = I_ESP_L3y5z_F3x_a+ABZ*I_ESP_K3y4z_F3x_a;
    Double I_ESP_K2y5z_G3xz_a = I_ESP_L2y6z_F3x_a+ABZ*I_ESP_K2y5z_F3x_a;
    Double I_ESP_Ky6z_G3xz_a = I_ESP_Ly7z_F3x_a+ABZ*I_ESP_Ky6z_F3x_a;
    Double I_ESP_K7z_G3xz_a = I_ESP_L8z_F3x_a+ABZ*I_ESP_K7z_F3x_a;
    Double I_ESP_K6xz_G2x2y_a = I_ESP_L6xyz_F2xy_a+ABY*I_ESP_K6xz_F2xy_a;
    Double I_ESP_K5xyz_G2x2y_a = I_ESP_L5x2yz_F2xy_a+ABY*I_ESP_K5xyz_F2xy_a;
    Double I_ESP_K5x2z_G2x2y_a = I_ESP_L5xy2z_F2xy_a+ABY*I_ESP_K5x2z_F2xy_a;
    Double I_ESP_K4x2yz_G2x2y_a = I_ESP_L4x3yz_F2xy_a+ABY*I_ESP_K4x2yz_F2xy_a;
    Double I_ESP_K4xy2z_G2x2y_a = I_ESP_L4x2y2z_F2xy_a+ABY*I_ESP_K4xy2z_F2xy_a;
    Double I_ESP_K4x3z_G2x2y_a = I_ESP_L4xy3z_F2xy_a+ABY*I_ESP_K4x3z_F2xy_a;
    Double I_ESP_K3x3yz_G2x2y_a = I_ESP_L3x4yz_F2xy_a+ABY*I_ESP_K3x3yz_F2xy_a;
    Double I_ESP_K3x2y2z_G2x2y_a = I_ESP_L3x3y2z_F2xy_a+ABY*I_ESP_K3x2y2z_F2xy_a;
    Double I_ESP_K3xy3z_G2x2y_a = I_ESP_L3x2y3z_F2xy_a+ABY*I_ESP_K3xy3z_F2xy_a;
    Double I_ESP_K3x4z_G2x2y_a = I_ESP_L3xy4z_F2xy_a+ABY*I_ESP_K3x4z_F2xy_a;
    Double I_ESP_K2x4yz_G2x2y_a = I_ESP_L2x5yz_F2xy_a+ABY*I_ESP_K2x4yz_F2xy_a;
    Double I_ESP_K2x3y2z_G2x2y_a = I_ESP_L2x4y2z_F2xy_a+ABY*I_ESP_K2x3y2z_F2xy_a;
    Double I_ESP_K2x2y3z_G2x2y_a = I_ESP_L2x3y3z_F2xy_a+ABY*I_ESP_K2x2y3z_F2xy_a;
    Double I_ESP_K2xy4z_G2x2y_a = I_ESP_L2x2y4z_F2xy_a+ABY*I_ESP_K2xy4z_F2xy_a;
    Double I_ESP_K2x5z_G2x2y_a = I_ESP_L2xy5z_F2xy_a+ABY*I_ESP_K2x5z_F2xy_a;
    Double I_ESP_Kx5yz_G2x2y_a = I_ESP_Lx6yz_F2xy_a+ABY*I_ESP_Kx5yz_F2xy_a;
    Double I_ESP_Kx4y2z_G2x2y_a = I_ESP_Lx5y2z_F2xy_a+ABY*I_ESP_Kx4y2z_F2xy_a;
    Double I_ESP_Kx3y3z_G2x2y_a = I_ESP_Lx4y3z_F2xy_a+ABY*I_ESP_Kx3y3z_F2xy_a;
    Double I_ESP_Kx2y4z_G2x2y_a = I_ESP_Lx3y4z_F2xy_a+ABY*I_ESP_Kx2y4z_F2xy_a;
    Double I_ESP_Kxy5z_G2x2y_a = I_ESP_Lx2y5z_F2xy_a+ABY*I_ESP_Kxy5z_F2xy_a;
    Double I_ESP_Kx6z_G2x2y_a = I_ESP_Lxy6z_F2xy_a+ABY*I_ESP_Kx6z_F2xy_a;
    Double I_ESP_K6yz_G2x2y_a = I_ESP_L7yz_F2xy_a+ABY*I_ESP_K6yz_F2xy_a;
    Double I_ESP_K5y2z_G2x2y_a = I_ESP_L6y2z_F2xy_a+ABY*I_ESP_K5y2z_F2xy_a;
    Double I_ESP_K4y3z_G2x2y_a = I_ESP_L5y3z_F2xy_a+ABY*I_ESP_K4y3z_F2xy_a;
    Double I_ESP_K3y4z_G2x2y_a = I_ESP_L4y4z_F2xy_a+ABY*I_ESP_K3y4z_F2xy_a;
    Double I_ESP_K2y5z_G2x2y_a = I_ESP_L3y5z_F2xy_a+ABY*I_ESP_K2y5z_F2xy_a;
    Double I_ESP_Ky6z_G2x2y_a = I_ESP_L2y6z_F2xy_a+ABY*I_ESP_Ky6z_F2xy_a;
    Double I_ESP_K7z_G2x2y_a = I_ESP_Ly7z_F2xy_a+ABY*I_ESP_K7z_F2xy_a;
    Double I_ESP_K6xy_G2x2z_a = I_ESP_L6xyz_F2xz_a+ABZ*I_ESP_K6xy_F2xz_a;
    Double I_ESP_K5x2y_G2x2z_a = I_ESP_L5x2yz_F2xz_a+ABZ*I_ESP_K5x2y_F2xz_a;
    Double I_ESP_K5xyz_G2x2z_a = I_ESP_L5xy2z_F2xz_a+ABZ*I_ESP_K5xyz_F2xz_a;
    Double I_ESP_K4x3y_G2x2z_a = I_ESP_L4x3yz_F2xz_a+ABZ*I_ESP_K4x3y_F2xz_a;
    Double I_ESP_K4x2yz_G2x2z_a = I_ESP_L4x2y2z_F2xz_a+ABZ*I_ESP_K4x2yz_F2xz_a;
    Double I_ESP_K4xy2z_G2x2z_a = I_ESP_L4xy3z_F2xz_a+ABZ*I_ESP_K4xy2z_F2xz_a;
    Double I_ESP_K3x4y_G2x2z_a = I_ESP_L3x4yz_F2xz_a+ABZ*I_ESP_K3x4y_F2xz_a;
    Double I_ESP_K3x3yz_G2x2z_a = I_ESP_L3x3y2z_F2xz_a+ABZ*I_ESP_K3x3yz_F2xz_a;
    Double I_ESP_K3x2y2z_G2x2z_a = I_ESP_L3x2y3z_F2xz_a+ABZ*I_ESP_K3x2y2z_F2xz_a;
    Double I_ESP_K3xy3z_G2x2z_a = I_ESP_L3xy4z_F2xz_a+ABZ*I_ESP_K3xy3z_F2xz_a;
    Double I_ESP_K2x5y_G2x2z_a = I_ESP_L2x5yz_F2xz_a+ABZ*I_ESP_K2x5y_F2xz_a;
    Double I_ESP_K2x4yz_G2x2z_a = I_ESP_L2x4y2z_F2xz_a+ABZ*I_ESP_K2x4yz_F2xz_a;
    Double I_ESP_K2x3y2z_G2x2z_a = I_ESP_L2x3y3z_F2xz_a+ABZ*I_ESP_K2x3y2z_F2xz_a;
    Double I_ESP_K2x2y3z_G2x2z_a = I_ESP_L2x2y4z_F2xz_a+ABZ*I_ESP_K2x2y3z_F2xz_a;
    Double I_ESP_K2xy4z_G2x2z_a = I_ESP_L2xy5z_F2xz_a+ABZ*I_ESP_K2xy4z_F2xz_a;
    Double I_ESP_Kx6y_G2x2z_a = I_ESP_Lx6yz_F2xz_a+ABZ*I_ESP_Kx6y_F2xz_a;
    Double I_ESP_Kx5yz_G2x2z_a = I_ESP_Lx5y2z_F2xz_a+ABZ*I_ESP_Kx5yz_F2xz_a;
    Double I_ESP_Kx4y2z_G2x2z_a = I_ESP_Lx4y3z_F2xz_a+ABZ*I_ESP_Kx4y2z_F2xz_a;
    Double I_ESP_Kx3y3z_G2x2z_a = I_ESP_Lx3y4z_F2xz_a+ABZ*I_ESP_Kx3y3z_F2xz_a;
    Double I_ESP_Kx2y4z_G2x2z_a = I_ESP_Lx2y5z_F2xz_a+ABZ*I_ESP_Kx2y4z_F2xz_a;
    Double I_ESP_Kxy5z_G2x2z_a = I_ESP_Lxy6z_F2xz_a+ABZ*I_ESP_Kxy5z_F2xz_a;
    Double I_ESP_K7y_G2x2z_a = I_ESP_L7yz_F2xz_a+ABZ*I_ESP_K7y_F2xz_a;
    Double I_ESP_K6yz_G2x2z_a = I_ESP_L6y2z_F2xz_a+ABZ*I_ESP_K6yz_F2xz_a;
    Double I_ESP_K5y2z_G2x2z_a = I_ESP_L5y3z_F2xz_a+ABZ*I_ESP_K5y2z_F2xz_a;
    Double I_ESP_K4y3z_G2x2z_a = I_ESP_L4y4z_F2xz_a+ABZ*I_ESP_K4y3z_F2xz_a;
    Double I_ESP_K3y4z_G2x2z_a = I_ESP_L3y5z_F2xz_a+ABZ*I_ESP_K3y4z_F2xz_a;
    Double I_ESP_K2y5z_G2x2z_a = I_ESP_L2y6z_F2xz_a+ABZ*I_ESP_K2y5z_F2xz_a;
    Double I_ESP_Ky6z_G2x2z_a = I_ESP_Ly7z_F2xz_a+ABZ*I_ESP_Ky6z_F2xz_a;
    Double I_ESP_K7x_Gx3y_a = I_ESP_L8x_F3y_a+ABX*I_ESP_K7x_F3y_a;
    Double I_ESP_K6xy_Gx3y_a = I_ESP_L7xy_F3y_a+ABX*I_ESP_K6xy_F3y_a;
    Double I_ESP_K6xz_Gx3y_a = I_ESP_L7xz_F3y_a+ABX*I_ESP_K6xz_F3y_a;
    Double I_ESP_K5x2y_Gx3y_a = I_ESP_L6x2y_F3y_a+ABX*I_ESP_K5x2y_F3y_a;
    Double I_ESP_K5xyz_Gx3y_a = I_ESP_L6xyz_F3y_a+ABX*I_ESP_K5xyz_F3y_a;
    Double I_ESP_K5x2z_Gx3y_a = I_ESP_L6x2z_F3y_a+ABX*I_ESP_K5x2z_F3y_a;
    Double I_ESP_K4x3y_Gx3y_a = I_ESP_L5x3y_F3y_a+ABX*I_ESP_K4x3y_F3y_a;
    Double I_ESP_K4x2yz_Gx3y_a = I_ESP_L5x2yz_F3y_a+ABX*I_ESP_K4x2yz_F3y_a;
    Double I_ESP_K4xy2z_Gx3y_a = I_ESP_L5xy2z_F3y_a+ABX*I_ESP_K4xy2z_F3y_a;
    Double I_ESP_K4x3z_Gx3y_a = I_ESP_L5x3z_F3y_a+ABX*I_ESP_K4x3z_F3y_a;
    Double I_ESP_K3x4y_Gx3y_a = I_ESP_L4x4y_F3y_a+ABX*I_ESP_K3x4y_F3y_a;
    Double I_ESP_K3x3yz_Gx3y_a = I_ESP_L4x3yz_F3y_a+ABX*I_ESP_K3x3yz_F3y_a;
    Double I_ESP_K3x2y2z_Gx3y_a = I_ESP_L4x2y2z_F3y_a+ABX*I_ESP_K3x2y2z_F3y_a;
    Double I_ESP_K3xy3z_Gx3y_a = I_ESP_L4xy3z_F3y_a+ABX*I_ESP_K3xy3z_F3y_a;
    Double I_ESP_K3x4z_Gx3y_a = I_ESP_L4x4z_F3y_a+ABX*I_ESP_K3x4z_F3y_a;
    Double I_ESP_K2x5y_Gx3y_a = I_ESP_L3x5y_F3y_a+ABX*I_ESP_K2x5y_F3y_a;
    Double I_ESP_K2x4yz_Gx3y_a = I_ESP_L3x4yz_F3y_a+ABX*I_ESP_K2x4yz_F3y_a;
    Double I_ESP_K2x3y2z_Gx3y_a = I_ESP_L3x3y2z_F3y_a+ABX*I_ESP_K2x3y2z_F3y_a;
    Double I_ESP_K2x2y3z_Gx3y_a = I_ESP_L3x2y3z_F3y_a+ABX*I_ESP_K2x2y3z_F3y_a;
    Double I_ESP_K2xy4z_Gx3y_a = I_ESP_L3xy4z_F3y_a+ABX*I_ESP_K2xy4z_F3y_a;
    Double I_ESP_K2x5z_Gx3y_a = I_ESP_L3x5z_F3y_a+ABX*I_ESP_K2x5z_F3y_a;
    Double I_ESP_Kx6y_Gx3y_a = I_ESP_L2x6y_F3y_a+ABX*I_ESP_Kx6y_F3y_a;
    Double I_ESP_Kx5yz_Gx3y_a = I_ESP_L2x5yz_F3y_a+ABX*I_ESP_Kx5yz_F3y_a;
    Double I_ESP_Kx4y2z_Gx3y_a = I_ESP_L2x4y2z_F3y_a+ABX*I_ESP_Kx4y2z_F3y_a;
    Double I_ESP_Kx3y3z_Gx3y_a = I_ESP_L2x3y3z_F3y_a+ABX*I_ESP_Kx3y3z_F3y_a;
    Double I_ESP_Kx2y4z_Gx3y_a = I_ESP_L2x2y4z_F3y_a+ABX*I_ESP_Kx2y4z_F3y_a;
    Double I_ESP_Kxy5z_Gx3y_a = I_ESP_L2xy5z_F3y_a+ABX*I_ESP_Kxy5z_F3y_a;
    Double I_ESP_Kx6z_Gx3y_a = I_ESP_L2x6z_F3y_a+ABX*I_ESP_Kx6z_F3y_a;
    Double I_ESP_K6yz_Gx3y_a = I_ESP_Lx6yz_F3y_a+ABX*I_ESP_K6yz_F3y_a;
    Double I_ESP_K5y2z_Gx3y_a = I_ESP_Lx5y2z_F3y_a+ABX*I_ESP_K5y2z_F3y_a;
    Double I_ESP_K4y3z_Gx3y_a = I_ESP_Lx4y3z_F3y_a+ABX*I_ESP_K4y3z_F3y_a;
    Double I_ESP_K3y4z_Gx3y_a = I_ESP_Lx3y4z_F3y_a+ABX*I_ESP_K3y4z_F3y_a;
    Double I_ESP_K2y5z_Gx3y_a = I_ESP_Lx2y5z_F3y_a+ABX*I_ESP_K2y5z_F3y_a;
    Double I_ESP_Ky6z_Gx3y_a = I_ESP_Lxy6z_F3y_a+ABX*I_ESP_Ky6z_F3y_a;
    Double I_ESP_K7z_Gx3y_a = I_ESP_Lx7z_F3y_a+ABX*I_ESP_K7z_F3y_a;
    Double I_ESP_K7x_Gx3z_a = I_ESP_L8x_F3z_a+ABX*I_ESP_K7x_F3z_a;
    Double I_ESP_K6xy_Gx3z_a = I_ESP_L7xy_F3z_a+ABX*I_ESP_K6xy_F3z_a;
    Double I_ESP_K6xz_Gx3z_a = I_ESP_L7xz_F3z_a+ABX*I_ESP_K6xz_F3z_a;
    Double I_ESP_K5x2y_Gx3z_a = I_ESP_L6x2y_F3z_a+ABX*I_ESP_K5x2y_F3z_a;
    Double I_ESP_K5xyz_Gx3z_a = I_ESP_L6xyz_F3z_a+ABX*I_ESP_K5xyz_F3z_a;
    Double I_ESP_K5x2z_Gx3z_a = I_ESP_L6x2z_F3z_a+ABX*I_ESP_K5x2z_F3z_a;
    Double I_ESP_K4x3y_Gx3z_a = I_ESP_L5x3y_F3z_a+ABX*I_ESP_K4x3y_F3z_a;
    Double I_ESP_K4x2yz_Gx3z_a = I_ESP_L5x2yz_F3z_a+ABX*I_ESP_K4x2yz_F3z_a;
    Double I_ESP_K4xy2z_Gx3z_a = I_ESP_L5xy2z_F3z_a+ABX*I_ESP_K4xy2z_F3z_a;
    Double I_ESP_K4x3z_Gx3z_a = I_ESP_L5x3z_F3z_a+ABX*I_ESP_K4x3z_F3z_a;
    Double I_ESP_K3x4y_Gx3z_a = I_ESP_L4x4y_F3z_a+ABX*I_ESP_K3x4y_F3z_a;
    Double I_ESP_K3x3yz_Gx3z_a = I_ESP_L4x3yz_F3z_a+ABX*I_ESP_K3x3yz_F3z_a;
    Double I_ESP_K3x2y2z_Gx3z_a = I_ESP_L4x2y2z_F3z_a+ABX*I_ESP_K3x2y2z_F3z_a;
    Double I_ESP_K3xy3z_Gx3z_a = I_ESP_L4xy3z_F3z_a+ABX*I_ESP_K3xy3z_F3z_a;
    Double I_ESP_K3x4z_Gx3z_a = I_ESP_L4x4z_F3z_a+ABX*I_ESP_K3x4z_F3z_a;
    Double I_ESP_K2x5y_Gx3z_a = I_ESP_L3x5y_F3z_a+ABX*I_ESP_K2x5y_F3z_a;
    Double I_ESP_K2x4yz_Gx3z_a = I_ESP_L3x4yz_F3z_a+ABX*I_ESP_K2x4yz_F3z_a;
    Double I_ESP_K2x3y2z_Gx3z_a = I_ESP_L3x3y2z_F3z_a+ABX*I_ESP_K2x3y2z_F3z_a;
    Double I_ESP_K2x2y3z_Gx3z_a = I_ESP_L3x2y3z_F3z_a+ABX*I_ESP_K2x2y3z_F3z_a;
    Double I_ESP_K2xy4z_Gx3z_a = I_ESP_L3xy4z_F3z_a+ABX*I_ESP_K2xy4z_F3z_a;
    Double I_ESP_K2x5z_Gx3z_a = I_ESP_L3x5z_F3z_a+ABX*I_ESP_K2x5z_F3z_a;
    Double I_ESP_Kx6y_Gx3z_a = I_ESP_L2x6y_F3z_a+ABX*I_ESP_Kx6y_F3z_a;
    Double I_ESP_Kx5yz_Gx3z_a = I_ESP_L2x5yz_F3z_a+ABX*I_ESP_Kx5yz_F3z_a;
    Double I_ESP_Kx4y2z_Gx3z_a = I_ESP_L2x4y2z_F3z_a+ABX*I_ESP_Kx4y2z_F3z_a;
    Double I_ESP_Kx3y3z_Gx3z_a = I_ESP_L2x3y3z_F3z_a+ABX*I_ESP_Kx3y3z_F3z_a;
    Double I_ESP_Kx2y4z_Gx3z_a = I_ESP_L2x2y4z_F3z_a+ABX*I_ESP_Kx2y4z_F3z_a;
    Double I_ESP_Kxy5z_Gx3z_a = I_ESP_L2xy5z_F3z_a+ABX*I_ESP_Kxy5z_F3z_a;
    Double I_ESP_Kx6z_Gx3z_a = I_ESP_L2x6z_F3z_a+ABX*I_ESP_Kx6z_F3z_a;
    Double I_ESP_K7y_Gx3z_a = I_ESP_Lx7y_F3z_a+ABX*I_ESP_K7y_F3z_a;
    Double I_ESP_K6yz_Gx3z_a = I_ESP_Lx6yz_F3z_a+ABX*I_ESP_K6yz_F3z_a;
    Double I_ESP_K5y2z_Gx3z_a = I_ESP_Lx5y2z_F3z_a+ABX*I_ESP_K5y2z_F3z_a;
    Double I_ESP_K4y3z_Gx3z_a = I_ESP_Lx4y3z_F3z_a+ABX*I_ESP_K4y3z_F3z_a;
    Double I_ESP_K3y4z_Gx3z_a = I_ESP_Lx3y4z_F3z_a+ABX*I_ESP_K3y4z_F3z_a;
    Double I_ESP_K2y5z_Gx3z_a = I_ESP_Lx2y5z_F3z_a+ABX*I_ESP_K2y5z_F3z_a;
    Double I_ESP_Ky6z_Gx3z_a = I_ESP_Lxy6z_F3z_a+ABX*I_ESP_Ky6z_F3z_a;
    Double I_ESP_K7x_G4y_a = I_ESP_L7xy_F3y_a+ABY*I_ESP_K7x_F3y_a;
    Double I_ESP_K6xy_G4y_a = I_ESP_L6x2y_F3y_a+ABY*I_ESP_K6xy_F3y_a;
    Double I_ESP_K6xz_G4y_a = I_ESP_L6xyz_F3y_a+ABY*I_ESP_K6xz_F3y_a;
    Double I_ESP_K5x2y_G4y_a = I_ESP_L5x3y_F3y_a+ABY*I_ESP_K5x2y_F3y_a;
    Double I_ESP_K5xyz_G4y_a = I_ESP_L5x2yz_F3y_a+ABY*I_ESP_K5xyz_F3y_a;
    Double I_ESP_K5x2z_G4y_a = I_ESP_L5xy2z_F3y_a+ABY*I_ESP_K5x2z_F3y_a;
    Double I_ESP_K4x3y_G4y_a = I_ESP_L4x4y_F3y_a+ABY*I_ESP_K4x3y_F3y_a;
    Double I_ESP_K4x2yz_G4y_a = I_ESP_L4x3yz_F3y_a+ABY*I_ESP_K4x2yz_F3y_a;
    Double I_ESP_K4xy2z_G4y_a = I_ESP_L4x2y2z_F3y_a+ABY*I_ESP_K4xy2z_F3y_a;
    Double I_ESP_K4x3z_G4y_a = I_ESP_L4xy3z_F3y_a+ABY*I_ESP_K4x3z_F3y_a;
    Double I_ESP_K3x4y_G4y_a = I_ESP_L3x5y_F3y_a+ABY*I_ESP_K3x4y_F3y_a;
    Double I_ESP_K3x3yz_G4y_a = I_ESP_L3x4yz_F3y_a+ABY*I_ESP_K3x3yz_F3y_a;
    Double I_ESP_K3x2y2z_G4y_a = I_ESP_L3x3y2z_F3y_a+ABY*I_ESP_K3x2y2z_F3y_a;
    Double I_ESP_K3xy3z_G4y_a = I_ESP_L3x2y3z_F3y_a+ABY*I_ESP_K3xy3z_F3y_a;
    Double I_ESP_K3x4z_G4y_a = I_ESP_L3xy4z_F3y_a+ABY*I_ESP_K3x4z_F3y_a;
    Double I_ESP_K2x5y_G4y_a = I_ESP_L2x6y_F3y_a+ABY*I_ESP_K2x5y_F3y_a;
    Double I_ESP_K2x4yz_G4y_a = I_ESP_L2x5yz_F3y_a+ABY*I_ESP_K2x4yz_F3y_a;
    Double I_ESP_K2x3y2z_G4y_a = I_ESP_L2x4y2z_F3y_a+ABY*I_ESP_K2x3y2z_F3y_a;
    Double I_ESP_K2x2y3z_G4y_a = I_ESP_L2x3y3z_F3y_a+ABY*I_ESP_K2x2y3z_F3y_a;
    Double I_ESP_K2xy4z_G4y_a = I_ESP_L2x2y4z_F3y_a+ABY*I_ESP_K2xy4z_F3y_a;
    Double I_ESP_K2x5z_G4y_a = I_ESP_L2xy5z_F3y_a+ABY*I_ESP_K2x5z_F3y_a;
    Double I_ESP_Kx6y_G4y_a = I_ESP_Lx7y_F3y_a+ABY*I_ESP_Kx6y_F3y_a;
    Double I_ESP_Kx5yz_G4y_a = I_ESP_Lx6yz_F3y_a+ABY*I_ESP_Kx5yz_F3y_a;
    Double I_ESP_Kx4y2z_G4y_a = I_ESP_Lx5y2z_F3y_a+ABY*I_ESP_Kx4y2z_F3y_a;
    Double I_ESP_Kx3y3z_G4y_a = I_ESP_Lx4y3z_F3y_a+ABY*I_ESP_Kx3y3z_F3y_a;
    Double I_ESP_Kx2y4z_G4y_a = I_ESP_Lx3y4z_F3y_a+ABY*I_ESP_Kx2y4z_F3y_a;
    Double I_ESP_Kxy5z_G4y_a = I_ESP_Lx2y5z_F3y_a+ABY*I_ESP_Kxy5z_F3y_a;
    Double I_ESP_Kx6z_G4y_a = I_ESP_Lxy6z_F3y_a+ABY*I_ESP_Kx6z_F3y_a;
    Double I_ESP_K7y_G4y_a = I_ESP_L8y_F3y_a+ABY*I_ESP_K7y_F3y_a;
    Double I_ESP_K6yz_G4y_a = I_ESP_L7yz_F3y_a+ABY*I_ESP_K6yz_F3y_a;
    Double I_ESP_K5y2z_G4y_a = I_ESP_L6y2z_F3y_a+ABY*I_ESP_K5y2z_F3y_a;
    Double I_ESP_K4y3z_G4y_a = I_ESP_L5y3z_F3y_a+ABY*I_ESP_K4y3z_F3y_a;
    Double I_ESP_K3y4z_G4y_a = I_ESP_L4y4z_F3y_a+ABY*I_ESP_K3y4z_F3y_a;
    Double I_ESP_K2y5z_G4y_a = I_ESP_L3y5z_F3y_a+ABY*I_ESP_K2y5z_F3y_a;
    Double I_ESP_Ky6z_G4y_a = I_ESP_L2y6z_F3y_a+ABY*I_ESP_Ky6z_F3y_a;
    Double I_ESP_K7z_G4y_a = I_ESP_Ly7z_F3y_a+ABY*I_ESP_K7z_F3y_a;
    Double I_ESP_K6xz_G3yz_a = I_ESP_L6x2z_F3y_a+ABZ*I_ESP_K6xz_F3y_a;
    Double I_ESP_K5xyz_G3yz_a = I_ESP_L5xy2z_F3y_a+ABZ*I_ESP_K5xyz_F3y_a;
    Double I_ESP_K5x2z_G3yz_a = I_ESP_L5x3z_F3y_a+ABZ*I_ESP_K5x2z_F3y_a;
    Double I_ESP_K4x2yz_G3yz_a = I_ESP_L4x2y2z_F3y_a+ABZ*I_ESP_K4x2yz_F3y_a;
    Double I_ESP_K4xy2z_G3yz_a = I_ESP_L4xy3z_F3y_a+ABZ*I_ESP_K4xy2z_F3y_a;
    Double I_ESP_K4x3z_G3yz_a = I_ESP_L4x4z_F3y_a+ABZ*I_ESP_K4x3z_F3y_a;
    Double I_ESP_K3x3yz_G3yz_a = I_ESP_L3x3y2z_F3y_a+ABZ*I_ESP_K3x3yz_F3y_a;
    Double I_ESP_K3x2y2z_G3yz_a = I_ESP_L3x2y3z_F3y_a+ABZ*I_ESP_K3x2y2z_F3y_a;
    Double I_ESP_K3xy3z_G3yz_a = I_ESP_L3xy4z_F3y_a+ABZ*I_ESP_K3xy3z_F3y_a;
    Double I_ESP_K3x4z_G3yz_a = I_ESP_L3x5z_F3y_a+ABZ*I_ESP_K3x4z_F3y_a;
    Double I_ESP_K2x4yz_G3yz_a = I_ESP_L2x4y2z_F3y_a+ABZ*I_ESP_K2x4yz_F3y_a;
    Double I_ESP_K2x3y2z_G3yz_a = I_ESP_L2x3y3z_F3y_a+ABZ*I_ESP_K2x3y2z_F3y_a;
    Double I_ESP_K2x2y3z_G3yz_a = I_ESP_L2x2y4z_F3y_a+ABZ*I_ESP_K2x2y3z_F3y_a;
    Double I_ESP_K2xy4z_G3yz_a = I_ESP_L2xy5z_F3y_a+ABZ*I_ESP_K2xy4z_F3y_a;
    Double I_ESP_K2x5z_G3yz_a = I_ESP_L2x6z_F3y_a+ABZ*I_ESP_K2x5z_F3y_a;
    Double I_ESP_Kx5yz_G3yz_a = I_ESP_Lx5y2z_F3y_a+ABZ*I_ESP_Kx5yz_F3y_a;
    Double I_ESP_Kx4y2z_G3yz_a = I_ESP_Lx4y3z_F3y_a+ABZ*I_ESP_Kx4y2z_F3y_a;
    Double I_ESP_Kx3y3z_G3yz_a = I_ESP_Lx3y4z_F3y_a+ABZ*I_ESP_Kx3y3z_F3y_a;
    Double I_ESP_Kx2y4z_G3yz_a = I_ESP_Lx2y5z_F3y_a+ABZ*I_ESP_Kx2y4z_F3y_a;
    Double I_ESP_Kxy5z_G3yz_a = I_ESP_Lxy6z_F3y_a+ABZ*I_ESP_Kxy5z_F3y_a;
    Double I_ESP_Kx6z_G3yz_a = I_ESP_Lx7z_F3y_a+ABZ*I_ESP_Kx6z_F3y_a;
    Double I_ESP_K6yz_G3yz_a = I_ESP_L6y2z_F3y_a+ABZ*I_ESP_K6yz_F3y_a;
    Double I_ESP_K5y2z_G3yz_a = I_ESP_L5y3z_F3y_a+ABZ*I_ESP_K5y2z_F3y_a;
    Double I_ESP_K4y3z_G3yz_a = I_ESP_L4y4z_F3y_a+ABZ*I_ESP_K4y3z_F3y_a;
    Double I_ESP_K3y4z_G3yz_a = I_ESP_L3y5z_F3y_a+ABZ*I_ESP_K3y4z_F3y_a;
    Double I_ESP_K2y5z_G3yz_a = I_ESP_L2y6z_F3y_a+ABZ*I_ESP_K2y5z_F3y_a;
    Double I_ESP_Ky6z_G3yz_a = I_ESP_Ly7z_F3y_a+ABZ*I_ESP_Ky6z_F3y_a;
    Double I_ESP_K7z_G3yz_a = I_ESP_L8z_F3y_a+ABZ*I_ESP_K7z_F3y_a;
    Double I_ESP_K7x_G2y2z_a = I_ESP_L7xz_F2yz_a+ABZ*I_ESP_K7x_F2yz_a;
    Double I_ESP_K6xy_G2y2z_a = I_ESP_L6xyz_F2yz_a+ABZ*I_ESP_K6xy_F2yz_a;
    Double I_ESP_K6xz_G2y2z_a = I_ESP_L6x2z_F2yz_a+ABZ*I_ESP_K6xz_F2yz_a;
    Double I_ESP_K5x2y_G2y2z_a = I_ESP_L5x2yz_F2yz_a+ABZ*I_ESP_K5x2y_F2yz_a;
    Double I_ESP_K5xyz_G2y2z_a = I_ESP_L5xy2z_F2yz_a+ABZ*I_ESP_K5xyz_F2yz_a;
    Double I_ESP_K5x2z_G2y2z_a = I_ESP_L5x3z_F2yz_a+ABZ*I_ESP_K5x2z_F2yz_a;
    Double I_ESP_K4x3y_G2y2z_a = I_ESP_L4x3yz_F2yz_a+ABZ*I_ESP_K4x3y_F2yz_a;
    Double I_ESP_K4x2yz_G2y2z_a = I_ESP_L4x2y2z_F2yz_a+ABZ*I_ESP_K4x2yz_F2yz_a;
    Double I_ESP_K4xy2z_G2y2z_a = I_ESP_L4xy3z_F2yz_a+ABZ*I_ESP_K4xy2z_F2yz_a;
    Double I_ESP_K4x3z_G2y2z_a = I_ESP_L4x4z_F2yz_a+ABZ*I_ESP_K4x3z_F2yz_a;
    Double I_ESP_K3x4y_G2y2z_a = I_ESP_L3x4yz_F2yz_a+ABZ*I_ESP_K3x4y_F2yz_a;
    Double I_ESP_K3x3yz_G2y2z_a = I_ESP_L3x3y2z_F2yz_a+ABZ*I_ESP_K3x3yz_F2yz_a;
    Double I_ESP_K3x2y2z_G2y2z_a = I_ESP_L3x2y3z_F2yz_a+ABZ*I_ESP_K3x2y2z_F2yz_a;
    Double I_ESP_K3xy3z_G2y2z_a = I_ESP_L3xy4z_F2yz_a+ABZ*I_ESP_K3xy3z_F2yz_a;
    Double I_ESP_K3x4z_G2y2z_a = I_ESP_L3x5z_F2yz_a+ABZ*I_ESP_K3x4z_F2yz_a;
    Double I_ESP_K2x5y_G2y2z_a = I_ESP_L2x5yz_F2yz_a+ABZ*I_ESP_K2x5y_F2yz_a;
    Double I_ESP_K2x4yz_G2y2z_a = I_ESP_L2x4y2z_F2yz_a+ABZ*I_ESP_K2x4yz_F2yz_a;
    Double I_ESP_K2x3y2z_G2y2z_a = I_ESP_L2x3y3z_F2yz_a+ABZ*I_ESP_K2x3y2z_F2yz_a;
    Double I_ESP_K2x2y3z_G2y2z_a = I_ESP_L2x2y4z_F2yz_a+ABZ*I_ESP_K2x2y3z_F2yz_a;
    Double I_ESP_K2xy4z_G2y2z_a = I_ESP_L2xy5z_F2yz_a+ABZ*I_ESP_K2xy4z_F2yz_a;
    Double I_ESP_K2x5z_G2y2z_a = I_ESP_L2x6z_F2yz_a+ABZ*I_ESP_K2x5z_F2yz_a;
    Double I_ESP_Kx6y_G2y2z_a = I_ESP_Lx6yz_F2yz_a+ABZ*I_ESP_Kx6y_F2yz_a;
    Double I_ESP_Kx5yz_G2y2z_a = I_ESP_Lx5y2z_F2yz_a+ABZ*I_ESP_Kx5yz_F2yz_a;
    Double I_ESP_Kx4y2z_G2y2z_a = I_ESP_Lx4y3z_F2yz_a+ABZ*I_ESP_Kx4y2z_F2yz_a;
    Double I_ESP_Kx3y3z_G2y2z_a = I_ESP_Lx3y4z_F2yz_a+ABZ*I_ESP_Kx3y3z_F2yz_a;
    Double I_ESP_Kx2y4z_G2y2z_a = I_ESP_Lx2y5z_F2yz_a+ABZ*I_ESP_Kx2y4z_F2yz_a;
    Double I_ESP_Kxy5z_G2y2z_a = I_ESP_Lxy6z_F2yz_a+ABZ*I_ESP_Kxy5z_F2yz_a;
    Double I_ESP_Kx6z_G2y2z_a = I_ESP_Lx7z_F2yz_a+ABZ*I_ESP_Kx6z_F2yz_a;
    Double I_ESP_K6xy_Gy3z_a = I_ESP_L6x2y_F3z_a+ABY*I_ESP_K6xy_F3z_a;
    Double I_ESP_K5x2y_Gy3z_a = I_ESP_L5x3y_F3z_a+ABY*I_ESP_K5x2y_F3z_a;
    Double I_ESP_K5xyz_Gy3z_a = I_ESP_L5x2yz_F3z_a+ABY*I_ESP_K5xyz_F3z_a;
    Double I_ESP_K4x3y_Gy3z_a = I_ESP_L4x4y_F3z_a+ABY*I_ESP_K4x3y_F3z_a;
    Double I_ESP_K4x2yz_Gy3z_a = I_ESP_L4x3yz_F3z_a+ABY*I_ESP_K4x2yz_F3z_a;
    Double I_ESP_K4xy2z_Gy3z_a = I_ESP_L4x2y2z_F3z_a+ABY*I_ESP_K4xy2z_F3z_a;
    Double I_ESP_K3x4y_Gy3z_a = I_ESP_L3x5y_F3z_a+ABY*I_ESP_K3x4y_F3z_a;
    Double I_ESP_K3x3yz_Gy3z_a = I_ESP_L3x4yz_F3z_a+ABY*I_ESP_K3x3yz_F3z_a;
    Double I_ESP_K3x2y2z_Gy3z_a = I_ESP_L3x3y2z_F3z_a+ABY*I_ESP_K3x2y2z_F3z_a;
    Double I_ESP_K3xy3z_Gy3z_a = I_ESP_L3x2y3z_F3z_a+ABY*I_ESP_K3xy3z_F3z_a;
    Double I_ESP_K2x5y_Gy3z_a = I_ESP_L2x6y_F3z_a+ABY*I_ESP_K2x5y_F3z_a;
    Double I_ESP_K2x4yz_Gy3z_a = I_ESP_L2x5yz_F3z_a+ABY*I_ESP_K2x4yz_F3z_a;
    Double I_ESP_K2x3y2z_Gy3z_a = I_ESP_L2x4y2z_F3z_a+ABY*I_ESP_K2x3y2z_F3z_a;
    Double I_ESP_K2x2y3z_Gy3z_a = I_ESP_L2x3y3z_F3z_a+ABY*I_ESP_K2x2y3z_F3z_a;
    Double I_ESP_K2xy4z_Gy3z_a = I_ESP_L2x2y4z_F3z_a+ABY*I_ESP_K2xy4z_F3z_a;
    Double I_ESP_Kx6y_Gy3z_a = I_ESP_Lx7y_F3z_a+ABY*I_ESP_Kx6y_F3z_a;
    Double I_ESP_Kx5yz_Gy3z_a = I_ESP_Lx6yz_F3z_a+ABY*I_ESP_Kx5yz_F3z_a;
    Double I_ESP_Kx4y2z_Gy3z_a = I_ESP_Lx5y2z_F3z_a+ABY*I_ESP_Kx4y2z_F3z_a;
    Double I_ESP_Kx3y3z_Gy3z_a = I_ESP_Lx4y3z_F3z_a+ABY*I_ESP_Kx3y3z_F3z_a;
    Double I_ESP_Kx2y4z_Gy3z_a = I_ESP_Lx3y4z_F3z_a+ABY*I_ESP_Kx2y4z_F3z_a;
    Double I_ESP_Kxy5z_Gy3z_a = I_ESP_Lx2y5z_F3z_a+ABY*I_ESP_Kxy5z_F3z_a;
    Double I_ESP_K7y_Gy3z_a = I_ESP_L8y_F3z_a+ABY*I_ESP_K7y_F3z_a;
    Double I_ESP_K6yz_Gy3z_a = I_ESP_L7yz_F3z_a+ABY*I_ESP_K6yz_F3z_a;
    Double I_ESP_K5y2z_Gy3z_a = I_ESP_L6y2z_F3z_a+ABY*I_ESP_K5y2z_F3z_a;
    Double I_ESP_K4y3z_Gy3z_a = I_ESP_L5y3z_F3z_a+ABY*I_ESP_K4y3z_F3z_a;
    Double I_ESP_K3y4z_Gy3z_a = I_ESP_L4y4z_F3z_a+ABY*I_ESP_K3y4z_F3z_a;
    Double I_ESP_K2y5z_Gy3z_a = I_ESP_L3y5z_F3z_a+ABY*I_ESP_K2y5z_F3z_a;
    Double I_ESP_Ky6z_Gy3z_a = I_ESP_L2y6z_F3z_a+ABY*I_ESP_Ky6z_F3z_a;
    Double I_ESP_K7x_G4z_a = I_ESP_L7xz_F3z_a+ABZ*I_ESP_K7x_F3z_a;
    Double I_ESP_K6xy_G4z_a = I_ESP_L6xyz_F3z_a+ABZ*I_ESP_K6xy_F3z_a;
    Double I_ESP_K6xz_G4z_a = I_ESP_L6x2z_F3z_a+ABZ*I_ESP_K6xz_F3z_a;
    Double I_ESP_K5x2y_G4z_a = I_ESP_L5x2yz_F3z_a+ABZ*I_ESP_K5x2y_F3z_a;
    Double I_ESP_K5xyz_G4z_a = I_ESP_L5xy2z_F3z_a+ABZ*I_ESP_K5xyz_F3z_a;
    Double I_ESP_K5x2z_G4z_a = I_ESP_L5x3z_F3z_a+ABZ*I_ESP_K5x2z_F3z_a;
    Double I_ESP_K4x3y_G4z_a = I_ESP_L4x3yz_F3z_a+ABZ*I_ESP_K4x3y_F3z_a;
    Double I_ESP_K4x2yz_G4z_a = I_ESP_L4x2y2z_F3z_a+ABZ*I_ESP_K4x2yz_F3z_a;
    Double I_ESP_K4xy2z_G4z_a = I_ESP_L4xy3z_F3z_a+ABZ*I_ESP_K4xy2z_F3z_a;
    Double I_ESP_K4x3z_G4z_a = I_ESP_L4x4z_F3z_a+ABZ*I_ESP_K4x3z_F3z_a;
    Double I_ESP_K3x4y_G4z_a = I_ESP_L3x4yz_F3z_a+ABZ*I_ESP_K3x4y_F3z_a;
    Double I_ESP_K3x3yz_G4z_a = I_ESP_L3x3y2z_F3z_a+ABZ*I_ESP_K3x3yz_F3z_a;
    Double I_ESP_K3x2y2z_G4z_a = I_ESP_L3x2y3z_F3z_a+ABZ*I_ESP_K3x2y2z_F3z_a;
    Double I_ESP_K3xy3z_G4z_a = I_ESP_L3xy4z_F3z_a+ABZ*I_ESP_K3xy3z_F3z_a;
    Double I_ESP_K3x4z_G4z_a = I_ESP_L3x5z_F3z_a+ABZ*I_ESP_K3x4z_F3z_a;
    Double I_ESP_K2x5y_G4z_a = I_ESP_L2x5yz_F3z_a+ABZ*I_ESP_K2x5y_F3z_a;
    Double I_ESP_K2x4yz_G4z_a = I_ESP_L2x4y2z_F3z_a+ABZ*I_ESP_K2x4yz_F3z_a;
    Double I_ESP_K2x3y2z_G4z_a = I_ESP_L2x3y3z_F3z_a+ABZ*I_ESP_K2x3y2z_F3z_a;
    Double I_ESP_K2x2y3z_G4z_a = I_ESP_L2x2y4z_F3z_a+ABZ*I_ESP_K2x2y3z_F3z_a;
    Double I_ESP_K2xy4z_G4z_a = I_ESP_L2xy5z_F3z_a+ABZ*I_ESP_K2xy4z_F3z_a;
    Double I_ESP_K2x5z_G4z_a = I_ESP_L2x6z_F3z_a+ABZ*I_ESP_K2x5z_F3z_a;
    Double I_ESP_Kx6y_G4z_a = I_ESP_Lx6yz_F3z_a+ABZ*I_ESP_Kx6y_F3z_a;
    Double I_ESP_Kx5yz_G4z_a = I_ESP_Lx5y2z_F3z_a+ABZ*I_ESP_Kx5yz_F3z_a;
    Double I_ESP_Kx4y2z_G4z_a = I_ESP_Lx4y3z_F3z_a+ABZ*I_ESP_Kx4y2z_F3z_a;
    Double I_ESP_Kx3y3z_G4z_a = I_ESP_Lx3y4z_F3z_a+ABZ*I_ESP_Kx3y3z_F3z_a;
    Double I_ESP_Kx2y4z_G4z_a = I_ESP_Lx2y5z_F3z_a+ABZ*I_ESP_Kx2y4z_F3z_a;
    Double I_ESP_Kxy5z_G4z_a = I_ESP_Lxy6z_F3z_a+ABZ*I_ESP_Kxy5z_F3z_a;
    Double I_ESP_Kx6z_G4z_a = I_ESP_Lx7z_F3z_a+ABZ*I_ESP_Kx6z_F3z_a;
    Double I_ESP_K7y_G4z_a = I_ESP_L7yz_F3z_a+ABZ*I_ESP_K7y_F3z_a;
    Double I_ESP_K6yz_G4z_a = I_ESP_L6y2z_F3z_a+ABZ*I_ESP_K6yz_F3z_a;
    Double I_ESP_K5y2z_G4z_a = I_ESP_L5y3z_F3z_a+ABZ*I_ESP_K5y2z_F3z_a;
    Double I_ESP_K4y3z_G4z_a = I_ESP_L4y4z_F3z_a+ABZ*I_ESP_K4y3z_F3z_a;
    Double I_ESP_K3y4z_G4z_a = I_ESP_L3y5z_F3z_a+ABZ*I_ESP_K3y4z_F3z_a;
    Double I_ESP_K2y5z_G4z_a = I_ESP_L2y6z_F3z_a+ABZ*I_ESP_K2y5z_F3z_a;
    Double I_ESP_Ky6z_G4z_a = I_ESP_Ly7z_F3z_a+ABZ*I_ESP_Ky6z_F3z_a;
    Double I_ESP_K7z_G4z_a = I_ESP_L8z_F3z_a+ABZ*I_ESP_K7z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_I_H_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_G_a
     * RHS shell quartet name: SQ_ESP_I_G_a
     ************************************************************/
    Double I_ESP_I6x_H5x_a = I_ESP_K7x_G4x_a+ABX*I_ESP_I6x_G4x_a;
    Double I_ESP_I5xy_H5x_a = I_ESP_K6xy_G4x_a+ABX*I_ESP_I5xy_G4x_a;
    Double I_ESP_I5xz_H5x_a = I_ESP_K6xz_G4x_a+ABX*I_ESP_I5xz_G4x_a;
    Double I_ESP_I4x2y_H5x_a = I_ESP_K5x2y_G4x_a+ABX*I_ESP_I4x2y_G4x_a;
    Double I_ESP_I4xyz_H5x_a = I_ESP_K5xyz_G4x_a+ABX*I_ESP_I4xyz_G4x_a;
    Double I_ESP_I4x2z_H5x_a = I_ESP_K5x2z_G4x_a+ABX*I_ESP_I4x2z_G4x_a;
    Double I_ESP_I3x3y_H5x_a = I_ESP_K4x3y_G4x_a+ABX*I_ESP_I3x3y_G4x_a;
    Double I_ESP_I3x2yz_H5x_a = I_ESP_K4x2yz_G4x_a+ABX*I_ESP_I3x2yz_G4x_a;
    Double I_ESP_I3xy2z_H5x_a = I_ESP_K4xy2z_G4x_a+ABX*I_ESP_I3xy2z_G4x_a;
    Double I_ESP_I3x3z_H5x_a = I_ESP_K4x3z_G4x_a+ABX*I_ESP_I3x3z_G4x_a;
    Double I_ESP_I2x4y_H5x_a = I_ESP_K3x4y_G4x_a+ABX*I_ESP_I2x4y_G4x_a;
    Double I_ESP_I2x3yz_H5x_a = I_ESP_K3x3yz_G4x_a+ABX*I_ESP_I2x3yz_G4x_a;
    Double I_ESP_I2x2y2z_H5x_a = I_ESP_K3x2y2z_G4x_a+ABX*I_ESP_I2x2y2z_G4x_a;
    Double I_ESP_I2xy3z_H5x_a = I_ESP_K3xy3z_G4x_a+ABX*I_ESP_I2xy3z_G4x_a;
    Double I_ESP_I2x4z_H5x_a = I_ESP_K3x4z_G4x_a+ABX*I_ESP_I2x4z_G4x_a;
    Double I_ESP_Ix5y_H5x_a = I_ESP_K2x5y_G4x_a+ABX*I_ESP_Ix5y_G4x_a;
    Double I_ESP_Ix4yz_H5x_a = I_ESP_K2x4yz_G4x_a+ABX*I_ESP_Ix4yz_G4x_a;
    Double I_ESP_Ix3y2z_H5x_a = I_ESP_K2x3y2z_G4x_a+ABX*I_ESP_Ix3y2z_G4x_a;
    Double I_ESP_Ix2y3z_H5x_a = I_ESP_K2x2y3z_G4x_a+ABX*I_ESP_Ix2y3z_G4x_a;
    Double I_ESP_Ixy4z_H5x_a = I_ESP_K2xy4z_G4x_a+ABX*I_ESP_Ixy4z_G4x_a;
    Double I_ESP_Ix5z_H5x_a = I_ESP_K2x5z_G4x_a+ABX*I_ESP_Ix5z_G4x_a;
    Double I_ESP_I6y_H5x_a = I_ESP_Kx6y_G4x_a+ABX*I_ESP_I6y_G4x_a;
    Double I_ESP_I5yz_H5x_a = I_ESP_Kx5yz_G4x_a+ABX*I_ESP_I5yz_G4x_a;
    Double I_ESP_I4y2z_H5x_a = I_ESP_Kx4y2z_G4x_a+ABX*I_ESP_I4y2z_G4x_a;
    Double I_ESP_I3y3z_H5x_a = I_ESP_Kx3y3z_G4x_a+ABX*I_ESP_I3y3z_G4x_a;
    Double I_ESP_I2y4z_H5x_a = I_ESP_Kx2y4z_G4x_a+ABX*I_ESP_I2y4z_G4x_a;
    Double I_ESP_Iy5z_H5x_a = I_ESP_Kxy5z_G4x_a+ABX*I_ESP_Iy5z_G4x_a;
    Double I_ESP_I6z_H5x_a = I_ESP_Kx6z_G4x_a+ABX*I_ESP_I6z_G4x_a;
    Double I_ESP_I6x_H4xy_a = I_ESP_K6xy_G4x_a+ABY*I_ESP_I6x_G4x_a;
    Double I_ESP_I5xy_H4xy_a = I_ESP_K5x2y_G4x_a+ABY*I_ESP_I5xy_G4x_a;
    Double I_ESP_I5xz_H4xy_a = I_ESP_K5xyz_G4x_a+ABY*I_ESP_I5xz_G4x_a;
    Double I_ESP_I4x2y_H4xy_a = I_ESP_K4x3y_G4x_a+ABY*I_ESP_I4x2y_G4x_a;
    Double I_ESP_I4xyz_H4xy_a = I_ESP_K4x2yz_G4x_a+ABY*I_ESP_I4xyz_G4x_a;
    Double I_ESP_I4x2z_H4xy_a = I_ESP_K4xy2z_G4x_a+ABY*I_ESP_I4x2z_G4x_a;
    Double I_ESP_I3x3y_H4xy_a = I_ESP_K3x4y_G4x_a+ABY*I_ESP_I3x3y_G4x_a;
    Double I_ESP_I3x2yz_H4xy_a = I_ESP_K3x3yz_G4x_a+ABY*I_ESP_I3x2yz_G4x_a;
    Double I_ESP_I3xy2z_H4xy_a = I_ESP_K3x2y2z_G4x_a+ABY*I_ESP_I3xy2z_G4x_a;
    Double I_ESP_I3x3z_H4xy_a = I_ESP_K3xy3z_G4x_a+ABY*I_ESP_I3x3z_G4x_a;
    Double I_ESP_I2x4y_H4xy_a = I_ESP_K2x5y_G4x_a+ABY*I_ESP_I2x4y_G4x_a;
    Double I_ESP_I2x3yz_H4xy_a = I_ESP_K2x4yz_G4x_a+ABY*I_ESP_I2x3yz_G4x_a;
    Double I_ESP_I2x2y2z_H4xy_a = I_ESP_K2x3y2z_G4x_a+ABY*I_ESP_I2x2y2z_G4x_a;
    Double I_ESP_I2xy3z_H4xy_a = I_ESP_K2x2y3z_G4x_a+ABY*I_ESP_I2xy3z_G4x_a;
    Double I_ESP_I2x4z_H4xy_a = I_ESP_K2xy4z_G4x_a+ABY*I_ESP_I2x4z_G4x_a;
    Double I_ESP_Ix5y_H4xy_a = I_ESP_Kx6y_G4x_a+ABY*I_ESP_Ix5y_G4x_a;
    Double I_ESP_Ix4yz_H4xy_a = I_ESP_Kx5yz_G4x_a+ABY*I_ESP_Ix4yz_G4x_a;
    Double I_ESP_Ix3y2z_H4xy_a = I_ESP_Kx4y2z_G4x_a+ABY*I_ESP_Ix3y2z_G4x_a;
    Double I_ESP_Ix2y3z_H4xy_a = I_ESP_Kx3y3z_G4x_a+ABY*I_ESP_Ix2y3z_G4x_a;
    Double I_ESP_Ixy4z_H4xy_a = I_ESP_Kx2y4z_G4x_a+ABY*I_ESP_Ixy4z_G4x_a;
    Double I_ESP_Ix5z_H4xy_a = I_ESP_Kxy5z_G4x_a+ABY*I_ESP_Ix5z_G4x_a;
    Double I_ESP_I6y_H4xy_a = I_ESP_K7y_G4x_a+ABY*I_ESP_I6y_G4x_a;
    Double I_ESP_I5yz_H4xy_a = I_ESP_K6yz_G4x_a+ABY*I_ESP_I5yz_G4x_a;
    Double I_ESP_I4y2z_H4xy_a = I_ESP_K5y2z_G4x_a+ABY*I_ESP_I4y2z_G4x_a;
    Double I_ESP_I3y3z_H4xy_a = I_ESP_K4y3z_G4x_a+ABY*I_ESP_I3y3z_G4x_a;
    Double I_ESP_I2y4z_H4xy_a = I_ESP_K3y4z_G4x_a+ABY*I_ESP_I2y4z_G4x_a;
    Double I_ESP_Iy5z_H4xy_a = I_ESP_K2y5z_G4x_a+ABY*I_ESP_Iy5z_G4x_a;
    Double I_ESP_I6z_H4xy_a = I_ESP_Ky6z_G4x_a+ABY*I_ESP_I6z_G4x_a;
    Double I_ESP_I6x_H4xz_a = I_ESP_K6xz_G4x_a+ABZ*I_ESP_I6x_G4x_a;
    Double I_ESP_I5xy_H4xz_a = I_ESP_K5xyz_G4x_a+ABZ*I_ESP_I5xy_G4x_a;
    Double I_ESP_I5xz_H4xz_a = I_ESP_K5x2z_G4x_a+ABZ*I_ESP_I5xz_G4x_a;
    Double I_ESP_I4x2y_H4xz_a = I_ESP_K4x2yz_G4x_a+ABZ*I_ESP_I4x2y_G4x_a;
    Double I_ESP_I4xyz_H4xz_a = I_ESP_K4xy2z_G4x_a+ABZ*I_ESP_I4xyz_G4x_a;
    Double I_ESP_I4x2z_H4xz_a = I_ESP_K4x3z_G4x_a+ABZ*I_ESP_I4x2z_G4x_a;
    Double I_ESP_I3x3y_H4xz_a = I_ESP_K3x3yz_G4x_a+ABZ*I_ESP_I3x3y_G4x_a;
    Double I_ESP_I3x2yz_H4xz_a = I_ESP_K3x2y2z_G4x_a+ABZ*I_ESP_I3x2yz_G4x_a;
    Double I_ESP_I3xy2z_H4xz_a = I_ESP_K3xy3z_G4x_a+ABZ*I_ESP_I3xy2z_G4x_a;
    Double I_ESP_I3x3z_H4xz_a = I_ESP_K3x4z_G4x_a+ABZ*I_ESP_I3x3z_G4x_a;
    Double I_ESP_I2x4y_H4xz_a = I_ESP_K2x4yz_G4x_a+ABZ*I_ESP_I2x4y_G4x_a;
    Double I_ESP_I2x3yz_H4xz_a = I_ESP_K2x3y2z_G4x_a+ABZ*I_ESP_I2x3yz_G4x_a;
    Double I_ESP_I2x2y2z_H4xz_a = I_ESP_K2x2y3z_G4x_a+ABZ*I_ESP_I2x2y2z_G4x_a;
    Double I_ESP_I2xy3z_H4xz_a = I_ESP_K2xy4z_G4x_a+ABZ*I_ESP_I2xy3z_G4x_a;
    Double I_ESP_I2x4z_H4xz_a = I_ESP_K2x5z_G4x_a+ABZ*I_ESP_I2x4z_G4x_a;
    Double I_ESP_Ix5y_H4xz_a = I_ESP_Kx5yz_G4x_a+ABZ*I_ESP_Ix5y_G4x_a;
    Double I_ESP_Ix4yz_H4xz_a = I_ESP_Kx4y2z_G4x_a+ABZ*I_ESP_Ix4yz_G4x_a;
    Double I_ESP_Ix3y2z_H4xz_a = I_ESP_Kx3y3z_G4x_a+ABZ*I_ESP_Ix3y2z_G4x_a;
    Double I_ESP_Ix2y3z_H4xz_a = I_ESP_Kx2y4z_G4x_a+ABZ*I_ESP_Ix2y3z_G4x_a;
    Double I_ESP_Ixy4z_H4xz_a = I_ESP_Kxy5z_G4x_a+ABZ*I_ESP_Ixy4z_G4x_a;
    Double I_ESP_Ix5z_H4xz_a = I_ESP_Kx6z_G4x_a+ABZ*I_ESP_Ix5z_G4x_a;
    Double I_ESP_I6y_H4xz_a = I_ESP_K6yz_G4x_a+ABZ*I_ESP_I6y_G4x_a;
    Double I_ESP_I5yz_H4xz_a = I_ESP_K5y2z_G4x_a+ABZ*I_ESP_I5yz_G4x_a;
    Double I_ESP_I4y2z_H4xz_a = I_ESP_K4y3z_G4x_a+ABZ*I_ESP_I4y2z_G4x_a;
    Double I_ESP_I3y3z_H4xz_a = I_ESP_K3y4z_G4x_a+ABZ*I_ESP_I3y3z_G4x_a;
    Double I_ESP_I2y4z_H4xz_a = I_ESP_K2y5z_G4x_a+ABZ*I_ESP_I2y4z_G4x_a;
    Double I_ESP_Iy5z_H4xz_a = I_ESP_Ky6z_G4x_a+ABZ*I_ESP_Iy5z_G4x_a;
    Double I_ESP_I6z_H4xz_a = I_ESP_K7z_G4x_a+ABZ*I_ESP_I6z_G4x_a;
    Double I_ESP_I6x_H3x2y_a = I_ESP_K6xy_G3xy_a+ABY*I_ESP_I6x_G3xy_a;
    Double I_ESP_I5xy_H3x2y_a = I_ESP_K5x2y_G3xy_a+ABY*I_ESP_I5xy_G3xy_a;
    Double I_ESP_I5xz_H3x2y_a = I_ESP_K5xyz_G3xy_a+ABY*I_ESP_I5xz_G3xy_a;
    Double I_ESP_I4x2y_H3x2y_a = I_ESP_K4x3y_G3xy_a+ABY*I_ESP_I4x2y_G3xy_a;
    Double I_ESP_I4xyz_H3x2y_a = I_ESP_K4x2yz_G3xy_a+ABY*I_ESP_I4xyz_G3xy_a;
    Double I_ESP_I4x2z_H3x2y_a = I_ESP_K4xy2z_G3xy_a+ABY*I_ESP_I4x2z_G3xy_a;
    Double I_ESP_I3x3y_H3x2y_a = I_ESP_K3x4y_G3xy_a+ABY*I_ESP_I3x3y_G3xy_a;
    Double I_ESP_I3x2yz_H3x2y_a = I_ESP_K3x3yz_G3xy_a+ABY*I_ESP_I3x2yz_G3xy_a;
    Double I_ESP_I3xy2z_H3x2y_a = I_ESP_K3x2y2z_G3xy_a+ABY*I_ESP_I3xy2z_G3xy_a;
    Double I_ESP_I3x3z_H3x2y_a = I_ESP_K3xy3z_G3xy_a+ABY*I_ESP_I3x3z_G3xy_a;
    Double I_ESP_I2x4y_H3x2y_a = I_ESP_K2x5y_G3xy_a+ABY*I_ESP_I2x4y_G3xy_a;
    Double I_ESP_I2x3yz_H3x2y_a = I_ESP_K2x4yz_G3xy_a+ABY*I_ESP_I2x3yz_G3xy_a;
    Double I_ESP_I2x2y2z_H3x2y_a = I_ESP_K2x3y2z_G3xy_a+ABY*I_ESP_I2x2y2z_G3xy_a;
    Double I_ESP_I2xy3z_H3x2y_a = I_ESP_K2x2y3z_G3xy_a+ABY*I_ESP_I2xy3z_G3xy_a;
    Double I_ESP_I2x4z_H3x2y_a = I_ESP_K2xy4z_G3xy_a+ABY*I_ESP_I2x4z_G3xy_a;
    Double I_ESP_Ix5y_H3x2y_a = I_ESP_Kx6y_G3xy_a+ABY*I_ESP_Ix5y_G3xy_a;
    Double I_ESP_Ix4yz_H3x2y_a = I_ESP_Kx5yz_G3xy_a+ABY*I_ESP_Ix4yz_G3xy_a;
    Double I_ESP_Ix3y2z_H3x2y_a = I_ESP_Kx4y2z_G3xy_a+ABY*I_ESP_Ix3y2z_G3xy_a;
    Double I_ESP_Ix2y3z_H3x2y_a = I_ESP_Kx3y3z_G3xy_a+ABY*I_ESP_Ix2y3z_G3xy_a;
    Double I_ESP_Ixy4z_H3x2y_a = I_ESP_Kx2y4z_G3xy_a+ABY*I_ESP_Ixy4z_G3xy_a;
    Double I_ESP_Ix5z_H3x2y_a = I_ESP_Kxy5z_G3xy_a+ABY*I_ESP_Ix5z_G3xy_a;
    Double I_ESP_I6y_H3x2y_a = I_ESP_K7y_G3xy_a+ABY*I_ESP_I6y_G3xy_a;
    Double I_ESP_I5yz_H3x2y_a = I_ESP_K6yz_G3xy_a+ABY*I_ESP_I5yz_G3xy_a;
    Double I_ESP_I4y2z_H3x2y_a = I_ESP_K5y2z_G3xy_a+ABY*I_ESP_I4y2z_G3xy_a;
    Double I_ESP_I3y3z_H3x2y_a = I_ESP_K4y3z_G3xy_a+ABY*I_ESP_I3y3z_G3xy_a;
    Double I_ESP_I2y4z_H3x2y_a = I_ESP_K3y4z_G3xy_a+ABY*I_ESP_I2y4z_G3xy_a;
    Double I_ESP_Iy5z_H3x2y_a = I_ESP_K2y5z_G3xy_a+ABY*I_ESP_Iy5z_G3xy_a;
    Double I_ESP_I6z_H3x2y_a = I_ESP_Ky6z_G3xy_a+ABY*I_ESP_I6z_G3xy_a;
    Double I_ESP_I6x_H3xyz_a = I_ESP_K6xz_G3xy_a+ABZ*I_ESP_I6x_G3xy_a;
    Double I_ESP_I5xy_H3xyz_a = I_ESP_K5xyz_G3xy_a+ABZ*I_ESP_I5xy_G3xy_a;
    Double I_ESP_I5xz_H3xyz_a = I_ESP_K5x2z_G3xy_a+ABZ*I_ESP_I5xz_G3xy_a;
    Double I_ESP_I4x2y_H3xyz_a = I_ESP_K4x2yz_G3xy_a+ABZ*I_ESP_I4x2y_G3xy_a;
    Double I_ESP_I4xyz_H3xyz_a = I_ESP_K4xy2z_G3xy_a+ABZ*I_ESP_I4xyz_G3xy_a;
    Double I_ESP_I4x2z_H3xyz_a = I_ESP_K4x3z_G3xy_a+ABZ*I_ESP_I4x2z_G3xy_a;
    Double I_ESP_I3x3y_H3xyz_a = I_ESP_K3x3yz_G3xy_a+ABZ*I_ESP_I3x3y_G3xy_a;
    Double I_ESP_I3x2yz_H3xyz_a = I_ESP_K3x2y2z_G3xy_a+ABZ*I_ESP_I3x2yz_G3xy_a;
    Double I_ESP_I3xy2z_H3xyz_a = I_ESP_K3xy3z_G3xy_a+ABZ*I_ESP_I3xy2z_G3xy_a;
    Double I_ESP_I3x3z_H3xyz_a = I_ESP_K3x4z_G3xy_a+ABZ*I_ESP_I3x3z_G3xy_a;
    Double I_ESP_I2x4y_H3xyz_a = I_ESP_K2x4yz_G3xy_a+ABZ*I_ESP_I2x4y_G3xy_a;
    Double I_ESP_I2x3yz_H3xyz_a = I_ESP_K2x3y2z_G3xy_a+ABZ*I_ESP_I2x3yz_G3xy_a;
    Double I_ESP_I2x2y2z_H3xyz_a = I_ESP_K2x2y3z_G3xy_a+ABZ*I_ESP_I2x2y2z_G3xy_a;
    Double I_ESP_I2xy3z_H3xyz_a = I_ESP_K2xy4z_G3xy_a+ABZ*I_ESP_I2xy3z_G3xy_a;
    Double I_ESP_I2x4z_H3xyz_a = I_ESP_K2x5z_G3xy_a+ABZ*I_ESP_I2x4z_G3xy_a;
    Double I_ESP_Ix5y_H3xyz_a = I_ESP_Kx5yz_G3xy_a+ABZ*I_ESP_Ix5y_G3xy_a;
    Double I_ESP_Ix4yz_H3xyz_a = I_ESP_Kx4y2z_G3xy_a+ABZ*I_ESP_Ix4yz_G3xy_a;
    Double I_ESP_Ix3y2z_H3xyz_a = I_ESP_Kx3y3z_G3xy_a+ABZ*I_ESP_Ix3y2z_G3xy_a;
    Double I_ESP_Ix2y3z_H3xyz_a = I_ESP_Kx2y4z_G3xy_a+ABZ*I_ESP_Ix2y3z_G3xy_a;
    Double I_ESP_Ixy4z_H3xyz_a = I_ESP_Kxy5z_G3xy_a+ABZ*I_ESP_Ixy4z_G3xy_a;
    Double I_ESP_Ix5z_H3xyz_a = I_ESP_Kx6z_G3xy_a+ABZ*I_ESP_Ix5z_G3xy_a;
    Double I_ESP_I6y_H3xyz_a = I_ESP_K6yz_G3xy_a+ABZ*I_ESP_I6y_G3xy_a;
    Double I_ESP_I5yz_H3xyz_a = I_ESP_K5y2z_G3xy_a+ABZ*I_ESP_I5yz_G3xy_a;
    Double I_ESP_I4y2z_H3xyz_a = I_ESP_K4y3z_G3xy_a+ABZ*I_ESP_I4y2z_G3xy_a;
    Double I_ESP_I3y3z_H3xyz_a = I_ESP_K3y4z_G3xy_a+ABZ*I_ESP_I3y3z_G3xy_a;
    Double I_ESP_I2y4z_H3xyz_a = I_ESP_K2y5z_G3xy_a+ABZ*I_ESP_I2y4z_G3xy_a;
    Double I_ESP_Iy5z_H3xyz_a = I_ESP_Ky6z_G3xy_a+ABZ*I_ESP_Iy5z_G3xy_a;
    Double I_ESP_I6z_H3xyz_a = I_ESP_K7z_G3xy_a+ABZ*I_ESP_I6z_G3xy_a;
    Double I_ESP_I6x_H3x2z_a = I_ESP_K6xz_G3xz_a+ABZ*I_ESP_I6x_G3xz_a;
    Double I_ESP_I5xy_H3x2z_a = I_ESP_K5xyz_G3xz_a+ABZ*I_ESP_I5xy_G3xz_a;
    Double I_ESP_I5xz_H3x2z_a = I_ESP_K5x2z_G3xz_a+ABZ*I_ESP_I5xz_G3xz_a;
    Double I_ESP_I4x2y_H3x2z_a = I_ESP_K4x2yz_G3xz_a+ABZ*I_ESP_I4x2y_G3xz_a;
    Double I_ESP_I4xyz_H3x2z_a = I_ESP_K4xy2z_G3xz_a+ABZ*I_ESP_I4xyz_G3xz_a;
    Double I_ESP_I4x2z_H3x2z_a = I_ESP_K4x3z_G3xz_a+ABZ*I_ESP_I4x2z_G3xz_a;
    Double I_ESP_I3x3y_H3x2z_a = I_ESP_K3x3yz_G3xz_a+ABZ*I_ESP_I3x3y_G3xz_a;
    Double I_ESP_I3x2yz_H3x2z_a = I_ESP_K3x2y2z_G3xz_a+ABZ*I_ESP_I3x2yz_G3xz_a;
    Double I_ESP_I3xy2z_H3x2z_a = I_ESP_K3xy3z_G3xz_a+ABZ*I_ESP_I3xy2z_G3xz_a;
    Double I_ESP_I3x3z_H3x2z_a = I_ESP_K3x4z_G3xz_a+ABZ*I_ESP_I3x3z_G3xz_a;
    Double I_ESP_I2x4y_H3x2z_a = I_ESP_K2x4yz_G3xz_a+ABZ*I_ESP_I2x4y_G3xz_a;
    Double I_ESP_I2x3yz_H3x2z_a = I_ESP_K2x3y2z_G3xz_a+ABZ*I_ESP_I2x3yz_G3xz_a;
    Double I_ESP_I2x2y2z_H3x2z_a = I_ESP_K2x2y3z_G3xz_a+ABZ*I_ESP_I2x2y2z_G3xz_a;
    Double I_ESP_I2xy3z_H3x2z_a = I_ESP_K2xy4z_G3xz_a+ABZ*I_ESP_I2xy3z_G3xz_a;
    Double I_ESP_I2x4z_H3x2z_a = I_ESP_K2x5z_G3xz_a+ABZ*I_ESP_I2x4z_G3xz_a;
    Double I_ESP_Ix5y_H3x2z_a = I_ESP_Kx5yz_G3xz_a+ABZ*I_ESP_Ix5y_G3xz_a;
    Double I_ESP_Ix4yz_H3x2z_a = I_ESP_Kx4y2z_G3xz_a+ABZ*I_ESP_Ix4yz_G3xz_a;
    Double I_ESP_Ix3y2z_H3x2z_a = I_ESP_Kx3y3z_G3xz_a+ABZ*I_ESP_Ix3y2z_G3xz_a;
    Double I_ESP_Ix2y3z_H3x2z_a = I_ESP_Kx2y4z_G3xz_a+ABZ*I_ESP_Ix2y3z_G3xz_a;
    Double I_ESP_Ixy4z_H3x2z_a = I_ESP_Kxy5z_G3xz_a+ABZ*I_ESP_Ixy4z_G3xz_a;
    Double I_ESP_Ix5z_H3x2z_a = I_ESP_Kx6z_G3xz_a+ABZ*I_ESP_Ix5z_G3xz_a;
    Double I_ESP_I6y_H3x2z_a = I_ESP_K6yz_G3xz_a+ABZ*I_ESP_I6y_G3xz_a;
    Double I_ESP_I5yz_H3x2z_a = I_ESP_K5y2z_G3xz_a+ABZ*I_ESP_I5yz_G3xz_a;
    Double I_ESP_I4y2z_H3x2z_a = I_ESP_K4y3z_G3xz_a+ABZ*I_ESP_I4y2z_G3xz_a;
    Double I_ESP_I3y3z_H3x2z_a = I_ESP_K3y4z_G3xz_a+ABZ*I_ESP_I3y3z_G3xz_a;
    Double I_ESP_I2y4z_H3x2z_a = I_ESP_K2y5z_G3xz_a+ABZ*I_ESP_I2y4z_G3xz_a;
    Double I_ESP_Iy5z_H3x2z_a = I_ESP_Ky6z_G3xz_a+ABZ*I_ESP_Iy5z_G3xz_a;
    Double I_ESP_I6z_H3x2z_a = I_ESP_K7z_G3xz_a+ABZ*I_ESP_I6z_G3xz_a;
    Double I_ESP_I6x_H2x3y_a = I_ESP_K7x_Gx3y_a+ABX*I_ESP_I6x_Gx3y_a;
    Double I_ESP_I5xy_H2x3y_a = I_ESP_K6xy_Gx3y_a+ABX*I_ESP_I5xy_Gx3y_a;
    Double I_ESP_I5xz_H2x3y_a = I_ESP_K6xz_Gx3y_a+ABX*I_ESP_I5xz_Gx3y_a;
    Double I_ESP_I4x2y_H2x3y_a = I_ESP_K5x2y_Gx3y_a+ABX*I_ESP_I4x2y_Gx3y_a;
    Double I_ESP_I4xyz_H2x3y_a = I_ESP_K5xyz_Gx3y_a+ABX*I_ESP_I4xyz_Gx3y_a;
    Double I_ESP_I4x2z_H2x3y_a = I_ESP_K5x2z_Gx3y_a+ABX*I_ESP_I4x2z_Gx3y_a;
    Double I_ESP_I3x3y_H2x3y_a = I_ESP_K4x3y_Gx3y_a+ABX*I_ESP_I3x3y_Gx3y_a;
    Double I_ESP_I3x2yz_H2x3y_a = I_ESP_K4x2yz_Gx3y_a+ABX*I_ESP_I3x2yz_Gx3y_a;
    Double I_ESP_I3xy2z_H2x3y_a = I_ESP_K4xy2z_Gx3y_a+ABX*I_ESP_I3xy2z_Gx3y_a;
    Double I_ESP_I3x3z_H2x3y_a = I_ESP_K4x3z_Gx3y_a+ABX*I_ESP_I3x3z_Gx3y_a;
    Double I_ESP_I2x4y_H2x3y_a = I_ESP_K3x4y_Gx3y_a+ABX*I_ESP_I2x4y_Gx3y_a;
    Double I_ESP_I2x3yz_H2x3y_a = I_ESP_K3x3yz_Gx3y_a+ABX*I_ESP_I2x3yz_Gx3y_a;
    Double I_ESP_I2x2y2z_H2x3y_a = I_ESP_K3x2y2z_Gx3y_a+ABX*I_ESP_I2x2y2z_Gx3y_a;
    Double I_ESP_I2xy3z_H2x3y_a = I_ESP_K3xy3z_Gx3y_a+ABX*I_ESP_I2xy3z_Gx3y_a;
    Double I_ESP_I2x4z_H2x3y_a = I_ESP_K3x4z_Gx3y_a+ABX*I_ESP_I2x4z_Gx3y_a;
    Double I_ESP_Ix5y_H2x3y_a = I_ESP_K2x5y_Gx3y_a+ABX*I_ESP_Ix5y_Gx3y_a;
    Double I_ESP_Ix4yz_H2x3y_a = I_ESP_K2x4yz_Gx3y_a+ABX*I_ESP_Ix4yz_Gx3y_a;
    Double I_ESP_Ix3y2z_H2x3y_a = I_ESP_K2x3y2z_Gx3y_a+ABX*I_ESP_Ix3y2z_Gx3y_a;
    Double I_ESP_Ix2y3z_H2x3y_a = I_ESP_K2x2y3z_Gx3y_a+ABX*I_ESP_Ix2y3z_Gx3y_a;
    Double I_ESP_Ixy4z_H2x3y_a = I_ESP_K2xy4z_Gx3y_a+ABX*I_ESP_Ixy4z_Gx3y_a;
    Double I_ESP_Ix5z_H2x3y_a = I_ESP_K2x5z_Gx3y_a+ABX*I_ESP_Ix5z_Gx3y_a;
    Double I_ESP_I6y_H2x3y_a = I_ESP_Kx6y_Gx3y_a+ABX*I_ESP_I6y_Gx3y_a;
    Double I_ESP_I5yz_H2x3y_a = I_ESP_Kx5yz_Gx3y_a+ABX*I_ESP_I5yz_Gx3y_a;
    Double I_ESP_I4y2z_H2x3y_a = I_ESP_Kx4y2z_Gx3y_a+ABX*I_ESP_I4y2z_Gx3y_a;
    Double I_ESP_I3y3z_H2x3y_a = I_ESP_Kx3y3z_Gx3y_a+ABX*I_ESP_I3y3z_Gx3y_a;
    Double I_ESP_I2y4z_H2x3y_a = I_ESP_Kx2y4z_Gx3y_a+ABX*I_ESP_I2y4z_Gx3y_a;
    Double I_ESP_Iy5z_H2x3y_a = I_ESP_Kxy5z_Gx3y_a+ABX*I_ESP_Iy5z_Gx3y_a;
    Double I_ESP_I6z_H2x3y_a = I_ESP_Kx6z_Gx3y_a+ABX*I_ESP_I6z_Gx3y_a;
    Double I_ESP_I6x_H2x2yz_a = I_ESP_K6xz_G2x2y_a+ABZ*I_ESP_I6x_G2x2y_a;
    Double I_ESP_I5xy_H2x2yz_a = I_ESP_K5xyz_G2x2y_a+ABZ*I_ESP_I5xy_G2x2y_a;
    Double I_ESP_I5xz_H2x2yz_a = I_ESP_K5x2z_G2x2y_a+ABZ*I_ESP_I5xz_G2x2y_a;
    Double I_ESP_I4x2y_H2x2yz_a = I_ESP_K4x2yz_G2x2y_a+ABZ*I_ESP_I4x2y_G2x2y_a;
    Double I_ESP_I4xyz_H2x2yz_a = I_ESP_K4xy2z_G2x2y_a+ABZ*I_ESP_I4xyz_G2x2y_a;
    Double I_ESP_I4x2z_H2x2yz_a = I_ESP_K4x3z_G2x2y_a+ABZ*I_ESP_I4x2z_G2x2y_a;
    Double I_ESP_I3x3y_H2x2yz_a = I_ESP_K3x3yz_G2x2y_a+ABZ*I_ESP_I3x3y_G2x2y_a;
    Double I_ESP_I3x2yz_H2x2yz_a = I_ESP_K3x2y2z_G2x2y_a+ABZ*I_ESP_I3x2yz_G2x2y_a;
    Double I_ESP_I3xy2z_H2x2yz_a = I_ESP_K3xy3z_G2x2y_a+ABZ*I_ESP_I3xy2z_G2x2y_a;
    Double I_ESP_I3x3z_H2x2yz_a = I_ESP_K3x4z_G2x2y_a+ABZ*I_ESP_I3x3z_G2x2y_a;
    Double I_ESP_I2x4y_H2x2yz_a = I_ESP_K2x4yz_G2x2y_a+ABZ*I_ESP_I2x4y_G2x2y_a;
    Double I_ESP_I2x3yz_H2x2yz_a = I_ESP_K2x3y2z_G2x2y_a+ABZ*I_ESP_I2x3yz_G2x2y_a;
    Double I_ESP_I2x2y2z_H2x2yz_a = I_ESP_K2x2y3z_G2x2y_a+ABZ*I_ESP_I2x2y2z_G2x2y_a;
    Double I_ESP_I2xy3z_H2x2yz_a = I_ESP_K2xy4z_G2x2y_a+ABZ*I_ESP_I2xy3z_G2x2y_a;
    Double I_ESP_I2x4z_H2x2yz_a = I_ESP_K2x5z_G2x2y_a+ABZ*I_ESP_I2x4z_G2x2y_a;
    Double I_ESP_Ix5y_H2x2yz_a = I_ESP_Kx5yz_G2x2y_a+ABZ*I_ESP_Ix5y_G2x2y_a;
    Double I_ESP_Ix4yz_H2x2yz_a = I_ESP_Kx4y2z_G2x2y_a+ABZ*I_ESP_Ix4yz_G2x2y_a;
    Double I_ESP_Ix3y2z_H2x2yz_a = I_ESP_Kx3y3z_G2x2y_a+ABZ*I_ESP_Ix3y2z_G2x2y_a;
    Double I_ESP_Ix2y3z_H2x2yz_a = I_ESP_Kx2y4z_G2x2y_a+ABZ*I_ESP_Ix2y3z_G2x2y_a;
    Double I_ESP_Ixy4z_H2x2yz_a = I_ESP_Kxy5z_G2x2y_a+ABZ*I_ESP_Ixy4z_G2x2y_a;
    Double I_ESP_Ix5z_H2x2yz_a = I_ESP_Kx6z_G2x2y_a+ABZ*I_ESP_Ix5z_G2x2y_a;
    Double I_ESP_I6y_H2x2yz_a = I_ESP_K6yz_G2x2y_a+ABZ*I_ESP_I6y_G2x2y_a;
    Double I_ESP_I5yz_H2x2yz_a = I_ESP_K5y2z_G2x2y_a+ABZ*I_ESP_I5yz_G2x2y_a;
    Double I_ESP_I4y2z_H2x2yz_a = I_ESP_K4y3z_G2x2y_a+ABZ*I_ESP_I4y2z_G2x2y_a;
    Double I_ESP_I3y3z_H2x2yz_a = I_ESP_K3y4z_G2x2y_a+ABZ*I_ESP_I3y3z_G2x2y_a;
    Double I_ESP_I2y4z_H2x2yz_a = I_ESP_K2y5z_G2x2y_a+ABZ*I_ESP_I2y4z_G2x2y_a;
    Double I_ESP_Iy5z_H2x2yz_a = I_ESP_Ky6z_G2x2y_a+ABZ*I_ESP_Iy5z_G2x2y_a;
    Double I_ESP_I6z_H2x2yz_a = I_ESP_K7z_G2x2y_a+ABZ*I_ESP_I6z_G2x2y_a;
    Double I_ESP_I6x_H2xy2z_a = I_ESP_K6xy_G2x2z_a+ABY*I_ESP_I6x_G2x2z_a;
    Double I_ESP_I5xy_H2xy2z_a = I_ESP_K5x2y_G2x2z_a+ABY*I_ESP_I5xy_G2x2z_a;
    Double I_ESP_I5xz_H2xy2z_a = I_ESP_K5xyz_G2x2z_a+ABY*I_ESP_I5xz_G2x2z_a;
    Double I_ESP_I4x2y_H2xy2z_a = I_ESP_K4x3y_G2x2z_a+ABY*I_ESP_I4x2y_G2x2z_a;
    Double I_ESP_I4xyz_H2xy2z_a = I_ESP_K4x2yz_G2x2z_a+ABY*I_ESP_I4xyz_G2x2z_a;
    Double I_ESP_I4x2z_H2xy2z_a = I_ESP_K4xy2z_G2x2z_a+ABY*I_ESP_I4x2z_G2x2z_a;
    Double I_ESP_I3x3y_H2xy2z_a = I_ESP_K3x4y_G2x2z_a+ABY*I_ESP_I3x3y_G2x2z_a;
    Double I_ESP_I3x2yz_H2xy2z_a = I_ESP_K3x3yz_G2x2z_a+ABY*I_ESP_I3x2yz_G2x2z_a;
    Double I_ESP_I3xy2z_H2xy2z_a = I_ESP_K3x2y2z_G2x2z_a+ABY*I_ESP_I3xy2z_G2x2z_a;
    Double I_ESP_I3x3z_H2xy2z_a = I_ESP_K3xy3z_G2x2z_a+ABY*I_ESP_I3x3z_G2x2z_a;
    Double I_ESP_I2x4y_H2xy2z_a = I_ESP_K2x5y_G2x2z_a+ABY*I_ESP_I2x4y_G2x2z_a;
    Double I_ESP_I2x3yz_H2xy2z_a = I_ESP_K2x4yz_G2x2z_a+ABY*I_ESP_I2x3yz_G2x2z_a;
    Double I_ESP_I2x2y2z_H2xy2z_a = I_ESP_K2x3y2z_G2x2z_a+ABY*I_ESP_I2x2y2z_G2x2z_a;
    Double I_ESP_I2xy3z_H2xy2z_a = I_ESP_K2x2y3z_G2x2z_a+ABY*I_ESP_I2xy3z_G2x2z_a;
    Double I_ESP_I2x4z_H2xy2z_a = I_ESP_K2xy4z_G2x2z_a+ABY*I_ESP_I2x4z_G2x2z_a;
    Double I_ESP_Ix5y_H2xy2z_a = I_ESP_Kx6y_G2x2z_a+ABY*I_ESP_Ix5y_G2x2z_a;
    Double I_ESP_Ix4yz_H2xy2z_a = I_ESP_Kx5yz_G2x2z_a+ABY*I_ESP_Ix4yz_G2x2z_a;
    Double I_ESP_Ix3y2z_H2xy2z_a = I_ESP_Kx4y2z_G2x2z_a+ABY*I_ESP_Ix3y2z_G2x2z_a;
    Double I_ESP_Ix2y3z_H2xy2z_a = I_ESP_Kx3y3z_G2x2z_a+ABY*I_ESP_Ix2y3z_G2x2z_a;
    Double I_ESP_Ixy4z_H2xy2z_a = I_ESP_Kx2y4z_G2x2z_a+ABY*I_ESP_Ixy4z_G2x2z_a;
    Double I_ESP_Ix5z_H2xy2z_a = I_ESP_Kxy5z_G2x2z_a+ABY*I_ESP_Ix5z_G2x2z_a;
    Double I_ESP_I6y_H2xy2z_a = I_ESP_K7y_G2x2z_a+ABY*I_ESP_I6y_G2x2z_a;
    Double I_ESP_I5yz_H2xy2z_a = I_ESP_K6yz_G2x2z_a+ABY*I_ESP_I5yz_G2x2z_a;
    Double I_ESP_I4y2z_H2xy2z_a = I_ESP_K5y2z_G2x2z_a+ABY*I_ESP_I4y2z_G2x2z_a;
    Double I_ESP_I3y3z_H2xy2z_a = I_ESP_K4y3z_G2x2z_a+ABY*I_ESP_I3y3z_G2x2z_a;
    Double I_ESP_I2y4z_H2xy2z_a = I_ESP_K3y4z_G2x2z_a+ABY*I_ESP_I2y4z_G2x2z_a;
    Double I_ESP_Iy5z_H2xy2z_a = I_ESP_K2y5z_G2x2z_a+ABY*I_ESP_Iy5z_G2x2z_a;
    Double I_ESP_I6z_H2xy2z_a = I_ESP_Ky6z_G2x2z_a+ABY*I_ESP_I6z_G2x2z_a;
    Double I_ESP_I6x_H2x3z_a = I_ESP_K7x_Gx3z_a+ABX*I_ESP_I6x_Gx3z_a;
    Double I_ESP_I5xy_H2x3z_a = I_ESP_K6xy_Gx3z_a+ABX*I_ESP_I5xy_Gx3z_a;
    Double I_ESP_I5xz_H2x3z_a = I_ESP_K6xz_Gx3z_a+ABX*I_ESP_I5xz_Gx3z_a;
    Double I_ESP_I4x2y_H2x3z_a = I_ESP_K5x2y_Gx3z_a+ABX*I_ESP_I4x2y_Gx3z_a;
    Double I_ESP_I4xyz_H2x3z_a = I_ESP_K5xyz_Gx3z_a+ABX*I_ESP_I4xyz_Gx3z_a;
    Double I_ESP_I4x2z_H2x3z_a = I_ESP_K5x2z_Gx3z_a+ABX*I_ESP_I4x2z_Gx3z_a;
    Double I_ESP_I3x3y_H2x3z_a = I_ESP_K4x3y_Gx3z_a+ABX*I_ESP_I3x3y_Gx3z_a;
    Double I_ESP_I3x2yz_H2x3z_a = I_ESP_K4x2yz_Gx3z_a+ABX*I_ESP_I3x2yz_Gx3z_a;
    Double I_ESP_I3xy2z_H2x3z_a = I_ESP_K4xy2z_Gx3z_a+ABX*I_ESP_I3xy2z_Gx3z_a;
    Double I_ESP_I3x3z_H2x3z_a = I_ESP_K4x3z_Gx3z_a+ABX*I_ESP_I3x3z_Gx3z_a;
    Double I_ESP_I2x4y_H2x3z_a = I_ESP_K3x4y_Gx3z_a+ABX*I_ESP_I2x4y_Gx3z_a;
    Double I_ESP_I2x3yz_H2x3z_a = I_ESP_K3x3yz_Gx3z_a+ABX*I_ESP_I2x3yz_Gx3z_a;
    Double I_ESP_I2x2y2z_H2x3z_a = I_ESP_K3x2y2z_Gx3z_a+ABX*I_ESP_I2x2y2z_Gx3z_a;
    Double I_ESP_I2xy3z_H2x3z_a = I_ESP_K3xy3z_Gx3z_a+ABX*I_ESP_I2xy3z_Gx3z_a;
    Double I_ESP_I2x4z_H2x3z_a = I_ESP_K3x4z_Gx3z_a+ABX*I_ESP_I2x4z_Gx3z_a;
    Double I_ESP_Ix5y_H2x3z_a = I_ESP_K2x5y_Gx3z_a+ABX*I_ESP_Ix5y_Gx3z_a;
    Double I_ESP_Ix4yz_H2x3z_a = I_ESP_K2x4yz_Gx3z_a+ABX*I_ESP_Ix4yz_Gx3z_a;
    Double I_ESP_Ix3y2z_H2x3z_a = I_ESP_K2x3y2z_Gx3z_a+ABX*I_ESP_Ix3y2z_Gx3z_a;
    Double I_ESP_Ix2y3z_H2x3z_a = I_ESP_K2x2y3z_Gx3z_a+ABX*I_ESP_Ix2y3z_Gx3z_a;
    Double I_ESP_Ixy4z_H2x3z_a = I_ESP_K2xy4z_Gx3z_a+ABX*I_ESP_Ixy4z_Gx3z_a;
    Double I_ESP_Ix5z_H2x3z_a = I_ESP_K2x5z_Gx3z_a+ABX*I_ESP_Ix5z_Gx3z_a;
    Double I_ESP_I6y_H2x3z_a = I_ESP_Kx6y_Gx3z_a+ABX*I_ESP_I6y_Gx3z_a;
    Double I_ESP_I5yz_H2x3z_a = I_ESP_Kx5yz_Gx3z_a+ABX*I_ESP_I5yz_Gx3z_a;
    Double I_ESP_I4y2z_H2x3z_a = I_ESP_Kx4y2z_Gx3z_a+ABX*I_ESP_I4y2z_Gx3z_a;
    Double I_ESP_I3y3z_H2x3z_a = I_ESP_Kx3y3z_Gx3z_a+ABX*I_ESP_I3y3z_Gx3z_a;
    Double I_ESP_I2y4z_H2x3z_a = I_ESP_Kx2y4z_Gx3z_a+ABX*I_ESP_I2y4z_Gx3z_a;
    Double I_ESP_Iy5z_H2x3z_a = I_ESP_Kxy5z_Gx3z_a+ABX*I_ESP_Iy5z_Gx3z_a;
    Double I_ESP_I6z_H2x3z_a = I_ESP_Kx6z_Gx3z_a+ABX*I_ESP_I6z_Gx3z_a;
    Double I_ESP_I6x_Hx4y_a = I_ESP_K7x_G4y_a+ABX*I_ESP_I6x_G4y_a;
    Double I_ESP_I5xy_Hx4y_a = I_ESP_K6xy_G4y_a+ABX*I_ESP_I5xy_G4y_a;
    Double I_ESP_I5xz_Hx4y_a = I_ESP_K6xz_G4y_a+ABX*I_ESP_I5xz_G4y_a;
    Double I_ESP_I4x2y_Hx4y_a = I_ESP_K5x2y_G4y_a+ABX*I_ESP_I4x2y_G4y_a;
    Double I_ESP_I4xyz_Hx4y_a = I_ESP_K5xyz_G4y_a+ABX*I_ESP_I4xyz_G4y_a;
    Double I_ESP_I4x2z_Hx4y_a = I_ESP_K5x2z_G4y_a+ABX*I_ESP_I4x2z_G4y_a;
    Double I_ESP_I3x3y_Hx4y_a = I_ESP_K4x3y_G4y_a+ABX*I_ESP_I3x3y_G4y_a;
    Double I_ESP_I3x2yz_Hx4y_a = I_ESP_K4x2yz_G4y_a+ABX*I_ESP_I3x2yz_G4y_a;
    Double I_ESP_I3xy2z_Hx4y_a = I_ESP_K4xy2z_G4y_a+ABX*I_ESP_I3xy2z_G4y_a;
    Double I_ESP_I3x3z_Hx4y_a = I_ESP_K4x3z_G4y_a+ABX*I_ESP_I3x3z_G4y_a;
    Double I_ESP_I2x4y_Hx4y_a = I_ESP_K3x4y_G4y_a+ABX*I_ESP_I2x4y_G4y_a;
    Double I_ESP_I2x3yz_Hx4y_a = I_ESP_K3x3yz_G4y_a+ABX*I_ESP_I2x3yz_G4y_a;
    Double I_ESP_I2x2y2z_Hx4y_a = I_ESP_K3x2y2z_G4y_a+ABX*I_ESP_I2x2y2z_G4y_a;
    Double I_ESP_I2xy3z_Hx4y_a = I_ESP_K3xy3z_G4y_a+ABX*I_ESP_I2xy3z_G4y_a;
    Double I_ESP_I2x4z_Hx4y_a = I_ESP_K3x4z_G4y_a+ABX*I_ESP_I2x4z_G4y_a;
    Double I_ESP_Ix5y_Hx4y_a = I_ESP_K2x5y_G4y_a+ABX*I_ESP_Ix5y_G4y_a;
    Double I_ESP_Ix4yz_Hx4y_a = I_ESP_K2x4yz_G4y_a+ABX*I_ESP_Ix4yz_G4y_a;
    Double I_ESP_Ix3y2z_Hx4y_a = I_ESP_K2x3y2z_G4y_a+ABX*I_ESP_Ix3y2z_G4y_a;
    Double I_ESP_Ix2y3z_Hx4y_a = I_ESP_K2x2y3z_G4y_a+ABX*I_ESP_Ix2y3z_G4y_a;
    Double I_ESP_Ixy4z_Hx4y_a = I_ESP_K2xy4z_G4y_a+ABX*I_ESP_Ixy4z_G4y_a;
    Double I_ESP_Ix5z_Hx4y_a = I_ESP_K2x5z_G4y_a+ABX*I_ESP_Ix5z_G4y_a;
    Double I_ESP_I6y_Hx4y_a = I_ESP_Kx6y_G4y_a+ABX*I_ESP_I6y_G4y_a;
    Double I_ESP_I5yz_Hx4y_a = I_ESP_Kx5yz_G4y_a+ABX*I_ESP_I5yz_G4y_a;
    Double I_ESP_I4y2z_Hx4y_a = I_ESP_Kx4y2z_G4y_a+ABX*I_ESP_I4y2z_G4y_a;
    Double I_ESP_I3y3z_Hx4y_a = I_ESP_Kx3y3z_G4y_a+ABX*I_ESP_I3y3z_G4y_a;
    Double I_ESP_I2y4z_Hx4y_a = I_ESP_Kx2y4z_G4y_a+ABX*I_ESP_I2y4z_G4y_a;
    Double I_ESP_Iy5z_Hx4y_a = I_ESP_Kxy5z_G4y_a+ABX*I_ESP_Iy5z_G4y_a;
    Double I_ESP_I6z_Hx4y_a = I_ESP_Kx6z_G4y_a+ABX*I_ESP_I6z_G4y_a;
    Double I_ESP_I6x_Hx3yz_a = I_ESP_K6xz_Gx3y_a+ABZ*I_ESP_I6x_Gx3y_a;
    Double I_ESP_I5xy_Hx3yz_a = I_ESP_K5xyz_Gx3y_a+ABZ*I_ESP_I5xy_Gx3y_a;
    Double I_ESP_I5xz_Hx3yz_a = I_ESP_K5x2z_Gx3y_a+ABZ*I_ESP_I5xz_Gx3y_a;
    Double I_ESP_I4x2y_Hx3yz_a = I_ESP_K4x2yz_Gx3y_a+ABZ*I_ESP_I4x2y_Gx3y_a;
    Double I_ESP_I4xyz_Hx3yz_a = I_ESP_K4xy2z_Gx3y_a+ABZ*I_ESP_I4xyz_Gx3y_a;
    Double I_ESP_I4x2z_Hx3yz_a = I_ESP_K4x3z_Gx3y_a+ABZ*I_ESP_I4x2z_Gx3y_a;
    Double I_ESP_I3x3y_Hx3yz_a = I_ESP_K3x3yz_Gx3y_a+ABZ*I_ESP_I3x3y_Gx3y_a;
    Double I_ESP_I3x2yz_Hx3yz_a = I_ESP_K3x2y2z_Gx3y_a+ABZ*I_ESP_I3x2yz_Gx3y_a;
    Double I_ESP_I3xy2z_Hx3yz_a = I_ESP_K3xy3z_Gx3y_a+ABZ*I_ESP_I3xy2z_Gx3y_a;
    Double I_ESP_I3x3z_Hx3yz_a = I_ESP_K3x4z_Gx3y_a+ABZ*I_ESP_I3x3z_Gx3y_a;
    Double I_ESP_I2x4y_Hx3yz_a = I_ESP_K2x4yz_Gx3y_a+ABZ*I_ESP_I2x4y_Gx3y_a;
    Double I_ESP_I2x3yz_Hx3yz_a = I_ESP_K2x3y2z_Gx3y_a+ABZ*I_ESP_I2x3yz_Gx3y_a;
    Double I_ESP_I2x2y2z_Hx3yz_a = I_ESP_K2x2y3z_Gx3y_a+ABZ*I_ESP_I2x2y2z_Gx3y_a;
    Double I_ESP_I2xy3z_Hx3yz_a = I_ESP_K2xy4z_Gx3y_a+ABZ*I_ESP_I2xy3z_Gx3y_a;
    Double I_ESP_I2x4z_Hx3yz_a = I_ESP_K2x5z_Gx3y_a+ABZ*I_ESP_I2x4z_Gx3y_a;
    Double I_ESP_Ix5y_Hx3yz_a = I_ESP_Kx5yz_Gx3y_a+ABZ*I_ESP_Ix5y_Gx3y_a;
    Double I_ESP_Ix4yz_Hx3yz_a = I_ESP_Kx4y2z_Gx3y_a+ABZ*I_ESP_Ix4yz_Gx3y_a;
    Double I_ESP_Ix3y2z_Hx3yz_a = I_ESP_Kx3y3z_Gx3y_a+ABZ*I_ESP_Ix3y2z_Gx3y_a;
    Double I_ESP_Ix2y3z_Hx3yz_a = I_ESP_Kx2y4z_Gx3y_a+ABZ*I_ESP_Ix2y3z_Gx3y_a;
    Double I_ESP_Ixy4z_Hx3yz_a = I_ESP_Kxy5z_Gx3y_a+ABZ*I_ESP_Ixy4z_Gx3y_a;
    Double I_ESP_Ix5z_Hx3yz_a = I_ESP_Kx6z_Gx3y_a+ABZ*I_ESP_Ix5z_Gx3y_a;
    Double I_ESP_I6y_Hx3yz_a = I_ESP_K6yz_Gx3y_a+ABZ*I_ESP_I6y_Gx3y_a;
    Double I_ESP_I5yz_Hx3yz_a = I_ESP_K5y2z_Gx3y_a+ABZ*I_ESP_I5yz_Gx3y_a;
    Double I_ESP_I4y2z_Hx3yz_a = I_ESP_K4y3z_Gx3y_a+ABZ*I_ESP_I4y2z_Gx3y_a;
    Double I_ESP_I3y3z_Hx3yz_a = I_ESP_K3y4z_Gx3y_a+ABZ*I_ESP_I3y3z_Gx3y_a;
    Double I_ESP_I2y4z_Hx3yz_a = I_ESP_K2y5z_Gx3y_a+ABZ*I_ESP_I2y4z_Gx3y_a;
    Double I_ESP_Iy5z_Hx3yz_a = I_ESP_Ky6z_Gx3y_a+ABZ*I_ESP_Iy5z_Gx3y_a;
    Double I_ESP_I6z_Hx3yz_a = I_ESP_K7z_Gx3y_a+ABZ*I_ESP_I6z_Gx3y_a;
    Double I_ESP_I6x_Hx2y2z_a = I_ESP_K7x_G2y2z_a+ABX*I_ESP_I6x_G2y2z_a;
    Double I_ESP_I5xy_Hx2y2z_a = I_ESP_K6xy_G2y2z_a+ABX*I_ESP_I5xy_G2y2z_a;
    Double I_ESP_I5xz_Hx2y2z_a = I_ESP_K6xz_G2y2z_a+ABX*I_ESP_I5xz_G2y2z_a;
    Double I_ESP_I4x2y_Hx2y2z_a = I_ESP_K5x2y_G2y2z_a+ABX*I_ESP_I4x2y_G2y2z_a;
    Double I_ESP_I4xyz_Hx2y2z_a = I_ESP_K5xyz_G2y2z_a+ABX*I_ESP_I4xyz_G2y2z_a;
    Double I_ESP_I4x2z_Hx2y2z_a = I_ESP_K5x2z_G2y2z_a+ABX*I_ESP_I4x2z_G2y2z_a;
    Double I_ESP_I3x3y_Hx2y2z_a = I_ESP_K4x3y_G2y2z_a+ABX*I_ESP_I3x3y_G2y2z_a;
    Double I_ESP_I3x2yz_Hx2y2z_a = I_ESP_K4x2yz_G2y2z_a+ABX*I_ESP_I3x2yz_G2y2z_a;
    Double I_ESP_I3xy2z_Hx2y2z_a = I_ESP_K4xy2z_G2y2z_a+ABX*I_ESP_I3xy2z_G2y2z_a;
    Double I_ESP_I3x3z_Hx2y2z_a = I_ESP_K4x3z_G2y2z_a+ABX*I_ESP_I3x3z_G2y2z_a;
    Double I_ESP_I2x4y_Hx2y2z_a = I_ESP_K3x4y_G2y2z_a+ABX*I_ESP_I2x4y_G2y2z_a;
    Double I_ESP_I2x3yz_Hx2y2z_a = I_ESP_K3x3yz_G2y2z_a+ABX*I_ESP_I2x3yz_G2y2z_a;
    Double I_ESP_I2x2y2z_Hx2y2z_a = I_ESP_K3x2y2z_G2y2z_a+ABX*I_ESP_I2x2y2z_G2y2z_a;
    Double I_ESP_I2xy3z_Hx2y2z_a = I_ESP_K3xy3z_G2y2z_a+ABX*I_ESP_I2xy3z_G2y2z_a;
    Double I_ESP_I2x4z_Hx2y2z_a = I_ESP_K3x4z_G2y2z_a+ABX*I_ESP_I2x4z_G2y2z_a;
    Double I_ESP_Ix5y_Hx2y2z_a = I_ESP_K2x5y_G2y2z_a+ABX*I_ESP_Ix5y_G2y2z_a;
    Double I_ESP_Ix4yz_Hx2y2z_a = I_ESP_K2x4yz_G2y2z_a+ABX*I_ESP_Ix4yz_G2y2z_a;
    Double I_ESP_Ix3y2z_Hx2y2z_a = I_ESP_K2x3y2z_G2y2z_a+ABX*I_ESP_Ix3y2z_G2y2z_a;
    Double I_ESP_Ix2y3z_Hx2y2z_a = I_ESP_K2x2y3z_G2y2z_a+ABX*I_ESP_Ix2y3z_G2y2z_a;
    Double I_ESP_Ixy4z_Hx2y2z_a = I_ESP_K2xy4z_G2y2z_a+ABX*I_ESP_Ixy4z_G2y2z_a;
    Double I_ESP_Ix5z_Hx2y2z_a = I_ESP_K2x5z_G2y2z_a+ABX*I_ESP_Ix5z_G2y2z_a;
    Double I_ESP_I6y_Hx2y2z_a = I_ESP_Kx6y_G2y2z_a+ABX*I_ESP_I6y_G2y2z_a;
    Double I_ESP_I5yz_Hx2y2z_a = I_ESP_Kx5yz_G2y2z_a+ABX*I_ESP_I5yz_G2y2z_a;
    Double I_ESP_I4y2z_Hx2y2z_a = I_ESP_Kx4y2z_G2y2z_a+ABX*I_ESP_I4y2z_G2y2z_a;
    Double I_ESP_I3y3z_Hx2y2z_a = I_ESP_Kx3y3z_G2y2z_a+ABX*I_ESP_I3y3z_G2y2z_a;
    Double I_ESP_I2y4z_Hx2y2z_a = I_ESP_Kx2y4z_G2y2z_a+ABX*I_ESP_I2y4z_G2y2z_a;
    Double I_ESP_Iy5z_Hx2y2z_a = I_ESP_Kxy5z_G2y2z_a+ABX*I_ESP_Iy5z_G2y2z_a;
    Double I_ESP_I6z_Hx2y2z_a = I_ESP_Kx6z_G2y2z_a+ABX*I_ESP_I6z_G2y2z_a;
    Double I_ESP_I6x_Hxy3z_a = I_ESP_K6xy_Gx3z_a+ABY*I_ESP_I6x_Gx3z_a;
    Double I_ESP_I5xy_Hxy3z_a = I_ESP_K5x2y_Gx3z_a+ABY*I_ESP_I5xy_Gx3z_a;
    Double I_ESP_I5xz_Hxy3z_a = I_ESP_K5xyz_Gx3z_a+ABY*I_ESP_I5xz_Gx3z_a;
    Double I_ESP_I4x2y_Hxy3z_a = I_ESP_K4x3y_Gx3z_a+ABY*I_ESP_I4x2y_Gx3z_a;
    Double I_ESP_I4xyz_Hxy3z_a = I_ESP_K4x2yz_Gx3z_a+ABY*I_ESP_I4xyz_Gx3z_a;
    Double I_ESP_I4x2z_Hxy3z_a = I_ESP_K4xy2z_Gx3z_a+ABY*I_ESP_I4x2z_Gx3z_a;
    Double I_ESP_I3x3y_Hxy3z_a = I_ESP_K3x4y_Gx3z_a+ABY*I_ESP_I3x3y_Gx3z_a;
    Double I_ESP_I3x2yz_Hxy3z_a = I_ESP_K3x3yz_Gx3z_a+ABY*I_ESP_I3x2yz_Gx3z_a;
    Double I_ESP_I3xy2z_Hxy3z_a = I_ESP_K3x2y2z_Gx3z_a+ABY*I_ESP_I3xy2z_Gx3z_a;
    Double I_ESP_I3x3z_Hxy3z_a = I_ESP_K3xy3z_Gx3z_a+ABY*I_ESP_I3x3z_Gx3z_a;
    Double I_ESP_I2x4y_Hxy3z_a = I_ESP_K2x5y_Gx3z_a+ABY*I_ESP_I2x4y_Gx3z_a;
    Double I_ESP_I2x3yz_Hxy3z_a = I_ESP_K2x4yz_Gx3z_a+ABY*I_ESP_I2x3yz_Gx3z_a;
    Double I_ESP_I2x2y2z_Hxy3z_a = I_ESP_K2x3y2z_Gx3z_a+ABY*I_ESP_I2x2y2z_Gx3z_a;
    Double I_ESP_I2xy3z_Hxy3z_a = I_ESP_K2x2y3z_Gx3z_a+ABY*I_ESP_I2xy3z_Gx3z_a;
    Double I_ESP_I2x4z_Hxy3z_a = I_ESP_K2xy4z_Gx3z_a+ABY*I_ESP_I2x4z_Gx3z_a;
    Double I_ESP_Ix5y_Hxy3z_a = I_ESP_Kx6y_Gx3z_a+ABY*I_ESP_Ix5y_Gx3z_a;
    Double I_ESP_Ix4yz_Hxy3z_a = I_ESP_Kx5yz_Gx3z_a+ABY*I_ESP_Ix4yz_Gx3z_a;
    Double I_ESP_Ix3y2z_Hxy3z_a = I_ESP_Kx4y2z_Gx3z_a+ABY*I_ESP_Ix3y2z_Gx3z_a;
    Double I_ESP_Ix2y3z_Hxy3z_a = I_ESP_Kx3y3z_Gx3z_a+ABY*I_ESP_Ix2y3z_Gx3z_a;
    Double I_ESP_Ixy4z_Hxy3z_a = I_ESP_Kx2y4z_Gx3z_a+ABY*I_ESP_Ixy4z_Gx3z_a;
    Double I_ESP_Ix5z_Hxy3z_a = I_ESP_Kxy5z_Gx3z_a+ABY*I_ESP_Ix5z_Gx3z_a;
    Double I_ESP_I6y_Hxy3z_a = I_ESP_K7y_Gx3z_a+ABY*I_ESP_I6y_Gx3z_a;
    Double I_ESP_I5yz_Hxy3z_a = I_ESP_K6yz_Gx3z_a+ABY*I_ESP_I5yz_Gx3z_a;
    Double I_ESP_I4y2z_Hxy3z_a = I_ESP_K5y2z_Gx3z_a+ABY*I_ESP_I4y2z_Gx3z_a;
    Double I_ESP_I3y3z_Hxy3z_a = I_ESP_K4y3z_Gx3z_a+ABY*I_ESP_I3y3z_Gx3z_a;
    Double I_ESP_I2y4z_Hxy3z_a = I_ESP_K3y4z_Gx3z_a+ABY*I_ESP_I2y4z_Gx3z_a;
    Double I_ESP_Iy5z_Hxy3z_a = I_ESP_K2y5z_Gx3z_a+ABY*I_ESP_Iy5z_Gx3z_a;
    Double I_ESP_I6z_Hxy3z_a = I_ESP_Ky6z_Gx3z_a+ABY*I_ESP_I6z_Gx3z_a;
    Double I_ESP_I6x_Hx4z_a = I_ESP_K7x_G4z_a+ABX*I_ESP_I6x_G4z_a;
    Double I_ESP_I5xy_Hx4z_a = I_ESP_K6xy_G4z_a+ABX*I_ESP_I5xy_G4z_a;
    Double I_ESP_I5xz_Hx4z_a = I_ESP_K6xz_G4z_a+ABX*I_ESP_I5xz_G4z_a;
    Double I_ESP_I4x2y_Hx4z_a = I_ESP_K5x2y_G4z_a+ABX*I_ESP_I4x2y_G4z_a;
    Double I_ESP_I4xyz_Hx4z_a = I_ESP_K5xyz_G4z_a+ABX*I_ESP_I4xyz_G4z_a;
    Double I_ESP_I4x2z_Hx4z_a = I_ESP_K5x2z_G4z_a+ABX*I_ESP_I4x2z_G4z_a;
    Double I_ESP_I3x3y_Hx4z_a = I_ESP_K4x3y_G4z_a+ABX*I_ESP_I3x3y_G4z_a;
    Double I_ESP_I3x2yz_Hx4z_a = I_ESP_K4x2yz_G4z_a+ABX*I_ESP_I3x2yz_G4z_a;
    Double I_ESP_I3xy2z_Hx4z_a = I_ESP_K4xy2z_G4z_a+ABX*I_ESP_I3xy2z_G4z_a;
    Double I_ESP_I3x3z_Hx4z_a = I_ESP_K4x3z_G4z_a+ABX*I_ESP_I3x3z_G4z_a;
    Double I_ESP_I2x4y_Hx4z_a = I_ESP_K3x4y_G4z_a+ABX*I_ESP_I2x4y_G4z_a;
    Double I_ESP_I2x3yz_Hx4z_a = I_ESP_K3x3yz_G4z_a+ABX*I_ESP_I2x3yz_G4z_a;
    Double I_ESP_I2x2y2z_Hx4z_a = I_ESP_K3x2y2z_G4z_a+ABX*I_ESP_I2x2y2z_G4z_a;
    Double I_ESP_I2xy3z_Hx4z_a = I_ESP_K3xy3z_G4z_a+ABX*I_ESP_I2xy3z_G4z_a;
    Double I_ESP_I2x4z_Hx4z_a = I_ESP_K3x4z_G4z_a+ABX*I_ESP_I2x4z_G4z_a;
    Double I_ESP_Ix5y_Hx4z_a = I_ESP_K2x5y_G4z_a+ABX*I_ESP_Ix5y_G4z_a;
    Double I_ESP_Ix4yz_Hx4z_a = I_ESP_K2x4yz_G4z_a+ABX*I_ESP_Ix4yz_G4z_a;
    Double I_ESP_Ix3y2z_Hx4z_a = I_ESP_K2x3y2z_G4z_a+ABX*I_ESP_Ix3y2z_G4z_a;
    Double I_ESP_Ix2y3z_Hx4z_a = I_ESP_K2x2y3z_G4z_a+ABX*I_ESP_Ix2y3z_G4z_a;
    Double I_ESP_Ixy4z_Hx4z_a = I_ESP_K2xy4z_G4z_a+ABX*I_ESP_Ixy4z_G4z_a;
    Double I_ESP_Ix5z_Hx4z_a = I_ESP_K2x5z_G4z_a+ABX*I_ESP_Ix5z_G4z_a;
    Double I_ESP_I6y_Hx4z_a = I_ESP_Kx6y_G4z_a+ABX*I_ESP_I6y_G4z_a;
    Double I_ESP_I5yz_Hx4z_a = I_ESP_Kx5yz_G4z_a+ABX*I_ESP_I5yz_G4z_a;
    Double I_ESP_I4y2z_Hx4z_a = I_ESP_Kx4y2z_G4z_a+ABX*I_ESP_I4y2z_G4z_a;
    Double I_ESP_I3y3z_Hx4z_a = I_ESP_Kx3y3z_G4z_a+ABX*I_ESP_I3y3z_G4z_a;
    Double I_ESP_I2y4z_Hx4z_a = I_ESP_Kx2y4z_G4z_a+ABX*I_ESP_I2y4z_G4z_a;
    Double I_ESP_Iy5z_Hx4z_a = I_ESP_Kxy5z_G4z_a+ABX*I_ESP_Iy5z_G4z_a;
    Double I_ESP_I6z_Hx4z_a = I_ESP_Kx6z_G4z_a+ABX*I_ESP_I6z_G4z_a;
    Double I_ESP_I6x_H5y_a = I_ESP_K6xy_G4y_a+ABY*I_ESP_I6x_G4y_a;
    Double I_ESP_I5xy_H5y_a = I_ESP_K5x2y_G4y_a+ABY*I_ESP_I5xy_G4y_a;
    Double I_ESP_I5xz_H5y_a = I_ESP_K5xyz_G4y_a+ABY*I_ESP_I5xz_G4y_a;
    Double I_ESP_I4x2y_H5y_a = I_ESP_K4x3y_G4y_a+ABY*I_ESP_I4x2y_G4y_a;
    Double I_ESP_I4xyz_H5y_a = I_ESP_K4x2yz_G4y_a+ABY*I_ESP_I4xyz_G4y_a;
    Double I_ESP_I4x2z_H5y_a = I_ESP_K4xy2z_G4y_a+ABY*I_ESP_I4x2z_G4y_a;
    Double I_ESP_I3x3y_H5y_a = I_ESP_K3x4y_G4y_a+ABY*I_ESP_I3x3y_G4y_a;
    Double I_ESP_I3x2yz_H5y_a = I_ESP_K3x3yz_G4y_a+ABY*I_ESP_I3x2yz_G4y_a;
    Double I_ESP_I3xy2z_H5y_a = I_ESP_K3x2y2z_G4y_a+ABY*I_ESP_I3xy2z_G4y_a;
    Double I_ESP_I3x3z_H5y_a = I_ESP_K3xy3z_G4y_a+ABY*I_ESP_I3x3z_G4y_a;
    Double I_ESP_I2x4y_H5y_a = I_ESP_K2x5y_G4y_a+ABY*I_ESP_I2x4y_G4y_a;
    Double I_ESP_I2x3yz_H5y_a = I_ESP_K2x4yz_G4y_a+ABY*I_ESP_I2x3yz_G4y_a;
    Double I_ESP_I2x2y2z_H5y_a = I_ESP_K2x3y2z_G4y_a+ABY*I_ESP_I2x2y2z_G4y_a;
    Double I_ESP_I2xy3z_H5y_a = I_ESP_K2x2y3z_G4y_a+ABY*I_ESP_I2xy3z_G4y_a;
    Double I_ESP_I2x4z_H5y_a = I_ESP_K2xy4z_G4y_a+ABY*I_ESP_I2x4z_G4y_a;
    Double I_ESP_Ix5y_H5y_a = I_ESP_Kx6y_G4y_a+ABY*I_ESP_Ix5y_G4y_a;
    Double I_ESP_Ix4yz_H5y_a = I_ESP_Kx5yz_G4y_a+ABY*I_ESP_Ix4yz_G4y_a;
    Double I_ESP_Ix3y2z_H5y_a = I_ESP_Kx4y2z_G4y_a+ABY*I_ESP_Ix3y2z_G4y_a;
    Double I_ESP_Ix2y3z_H5y_a = I_ESP_Kx3y3z_G4y_a+ABY*I_ESP_Ix2y3z_G4y_a;
    Double I_ESP_Ixy4z_H5y_a = I_ESP_Kx2y4z_G4y_a+ABY*I_ESP_Ixy4z_G4y_a;
    Double I_ESP_Ix5z_H5y_a = I_ESP_Kxy5z_G4y_a+ABY*I_ESP_Ix5z_G4y_a;
    Double I_ESP_I6y_H5y_a = I_ESP_K7y_G4y_a+ABY*I_ESP_I6y_G4y_a;
    Double I_ESP_I5yz_H5y_a = I_ESP_K6yz_G4y_a+ABY*I_ESP_I5yz_G4y_a;
    Double I_ESP_I4y2z_H5y_a = I_ESP_K5y2z_G4y_a+ABY*I_ESP_I4y2z_G4y_a;
    Double I_ESP_I3y3z_H5y_a = I_ESP_K4y3z_G4y_a+ABY*I_ESP_I3y3z_G4y_a;
    Double I_ESP_I2y4z_H5y_a = I_ESP_K3y4z_G4y_a+ABY*I_ESP_I2y4z_G4y_a;
    Double I_ESP_Iy5z_H5y_a = I_ESP_K2y5z_G4y_a+ABY*I_ESP_Iy5z_G4y_a;
    Double I_ESP_I6z_H5y_a = I_ESP_Ky6z_G4y_a+ABY*I_ESP_I6z_G4y_a;
    Double I_ESP_I6x_H4yz_a = I_ESP_K6xz_G4y_a+ABZ*I_ESP_I6x_G4y_a;
    Double I_ESP_I5xy_H4yz_a = I_ESP_K5xyz_G4y_a+ABZ*I_ESP_I5xy_G4y_a;
    Double I_ESP_I5xz_H4yz_a = I_ESP_K5x2z_G4y_a+ABZ*I_ESP_I5xz_G4y_a;
    Double I_ESP_I4x2y_H4yz_a = I_ESP_K4x2yz_G4y_a+ABZ*I_ESP_I4x2y_G4y_a;
    Double I_ESP_I4xyz_H4yz_a = I_ESP_K4xy2z_G4y_a+ABZ*I_ESP_I4xyz_G4y_a;
    Double I_ESP_I4x2z_H4yz_a = I_ESP_K4x3z_G4y_a+ABZ*I_ESP_I4x2z_G4y_a;
    Double I_ESP_I3x3y_H4yz_a = I_ESP_K3x3yz_G4y_a+ABZ*I_ESP_I3x3y_G4y_a;
    Double I_ESP_I3x2yz_H4yz_a = I_ESP_K3x2y2z_G4y_a+ABZ*I_ESP_I3x2yz_G4y_a;
    Double I_ESP_I3xy2z_H4yz_a = I_ESP_K3xy3z_G4y_a+ABZ*I_ESP_I3xy2z_G4y_a;
    Double I_ESP_I3x3z_H4yz_a = I_ESP_K3x4z_G4y_a+ABZ*I_ESP_I3x3z_G4y_a;
    Double I_ESP_I2x4y_H4yz_a = I_ESP_K2x4yz_G4y_a+ABZ*I_ESP_I2x4y_G4y_a;
    Double I_ESP_I2x3yz_H4yz_a = I_ESP_K2x3y2z_G4y_a+ABZ*I_ESP_I2x3yz_G4y_a;
    Double I_ESP_I2x2y2z_H4yz_a = I_ESP_K2x2y3z_G4y_a+ABZ*I_ESP_I2x2y2z_G4y_a;
    Double I_ESP_I2xy3z_H4yz_a = I_ESP_K2xy4z_G4y_a+ABZ*I_ESP_I2xy3z_G4y_a;
    Double I_ESP_I2x4z_H4yz_a = I_ESP_K2x5z_G4y_a+ABZ*I_ESP_I2x4z_G4y_a;
    Double I_ESP_Ix5y_H4yz_a = I_ESP_Kx5yz_G4y_a+ABZ*I_ESP_Ix5y_G4y_a;
    Double I_ESP_Ix4yz_H4yz_a = I_ESP_Kx4y2z_G4y_a+ABZ*I_ESP_Ix4yz_G4y_a;
    Double I_ESP_Ix3y2z_H4yz_a = I_ESP_Kx3y3z_G4y_a+ABZ*I_ESP_Ix3y2z_G4y_a;
    Double I_ESP_Ix2y3z_H4yz_a = I_ESP_Kx2y4z_G4y_a+ABZ*I_ESP_Ix2y3z_G4y_a;
    Double I_ESP_Ixy4z_H4yz_a = I_ESP_Kxy5z_G4y_a+ABZ*I_ESP_Ixy4z_G4y_a;
    Double I_ESP_Ix5z_H4yz_a = I_ESP_Kx6z_G4y_a+ABZ*I_ESP_Ix5z_G4y_a;
    Double I_ESP_I6y_H4yz_a = I_ESP_K6yz_G4y_a+ABZ*I_ESP_I6y_G4y_a;
    Double I_ESP_I5yz_H4yz_a = I_ESP_K5y2z_G4y_a+ABZ*I_ESP_I5yz_G4y_a;
    Double I_ESP_I4y2z_H4yz_a = I_ESP_K4y3z_G4y_a+ABZ*I_ESP_I4y2z_G4y_a;
    Double I_ESP_I3y3z_H4yz_a = I_ESP_K3y4z_G4y_a+ABZ*I_ESP_I3y3z_G4y_a;
    Double I_ESP_I2y4z_H4yz_a = I_ESP_K2y5z_G4y_a+ABZ*I_ESP_I2y4z_G4y_a;
    Double I_ESP_Iy5z_H4yz_a = I_ESP_Ky6z_G4y_a+ABZ*I_ESP_Iy5z_G4y_a;
    Double I_ESP_I6z_H4yz_a = I_ESP_K7z_G4y_a+ABZ*I_ESP_I6z_G4y_a;
    Double I_ESP_I6x_H3y2z_a = I_ESP_K6xz_G3yz_a+ABZ*I_ESP_I6x_G3yz_a;
    Double I_ESP_I5xy_H3y2z_a = I_ESP_K5xyz_G3yz_a+ABZ*I_ESP_I5xy_G3yz_a;
    Double I_ESP_I5xz_H3y2z_a = I_ESP_K5x2z_G3yz_a+ABZ*I_ESP_I5xz_G3yz_a;
    Double I_ESP_I4x2y_H3y2z_a = I_ESP_K4x2yz_G3yz_a+ABZ*I_ESP_I4x2y_G3yz_a;
    Double I_ESP_I4xyz_H3y2z_a = I_ESP_K4xy2z_G3yz_a+ABZ*I_ESP_I4xyz_G3yz_a;
    Double I_ESP_I4x2z_H3y2z_a = I_ESP_K4x3z_G3yz_a+ABZ*I_ESP_I4x2z_G3yz_a;
    Double I_ESP_I3x3y_H3y2z_a = I_ESP_K3x3yz_G3yz_a+ABZ*I_ESP_I3x3y_G3yz_a;
    Double I_ESP_I3x2yz_H3y2z_a = I_ESP_K3x2y2z_G3yz_a+ABZ*I_ESP_I3x2yz_G3yz_a;
    Double I_ESP_I3xy2z_H3y2z_a = I_ESP_K3xy3z_G3yz_a+ABZ*I_ESP_I3xy2z_G3yz_a;
    Double I_ESP_I3x3z_H3y2z_a = I_ESP_K3x4z_G3yz_a+ABZ*I_ESP_I3x3z_G3yz_a;
    Double I_ESP_I2x4y_H3y2z_a = I_ESP_K2x4yz_G3yz_a+ABZ*I_ESP_I2x4y_G3yz_a;
    Double I_ESP_I2x3yz_H3y2z_a = I_ESP_K2x3y2z_G3yz_a+ABZ*I_ESP_I2x3yz_G3yz_a;
    Double I_ESP_I2x2y2z_H3y2z_a = I_ESP_K2x2y3z_G3yz_a+ABZ*I_ESP_I2x2y2z_G3yz_a;
    Double I_ESP_I2xy3z_H3y2z_a = I_ESP_K2xy4z_G3yz_a+ABZ*I_ESP_I2xy3z_G3yz_a;
    Double I_ESP_I2x4z_H3y2z_a = I_ESP_K2x5z_G3yz_a+ABZ*I_ESP_I2x4z_G3yz_a;
    Double I_ESP_Ix5y_H3y2z_a = I_ESP_Kx5yz_G3yz_a+ABZ*I_ESP_Ix5y_G3yz_a;
    Double I_ESP_Ix4yz_H3y2z_a = I_ESP_Kx4y2z_G3yz_a+ABZ*I_ESP_Ix4yz_G3yz_a;
    Double I_ESP_Ix3y2z_H3y2z_a = I_ESP_Kx3y3z_G3yz_a+ABZ*I_ESP_Ix3y2z_G3yz_a;
    Double I_ESP_Ix2y3z_H3y2z_a = I_ESP_Kx2y4z_G3yz_a+ABZ*I_ESP_Ix2y3z_G3yz_a;
    Double I_ESP_Ixy4z_H3y2z_a = I_ESP_Kxy5z_G3yz_a+ABZ*I_ESP_Ixy4z_G3yz_a;
    Double I_ESP_Ix5z_H3y2z_a = I_ESP_Kx6z_G3yz_a+ABZ*I_ESP_Ix5z_G3yz_a;
    Double I_ESP_I6y_H3y2z_a = I_ESP_K6yz_G3yz_a+ABZ*I_ESP_I6y_G3yz_a;
    Double I_ESP_I5yz_H3y2z_a = I_ESP_K5y2z_G3yz_a+ABZ*I_ESP_I5yz_G3yz_a;
    Double I_ESP_I4y2z_H3y2z_a = I_ESP_K4y3z_G3yz_a+ABZ*I_ESP_I4y2z_G3yz_a;
    Double I_ESP_I3y3z_H3y2z_a = I_ESP_K3y4z_G3yz_a+ABZ*I_ESP_I3y3z_G3yz_a;
    Double I_ESP_I2y4z_H3y2z_a = I_ESP_K2y5z_G3yz_a+ABZ*I_ESP_I2y4z_G3yz_a;
    Double I_ESP_Iy5z_H3y2z_a = I_ESP_Ky6z_G3yz_a+ABZ*I_ESP_Iy5z_G3yz_a;
    Double I_ESP_I6z_H3y2z_a = I_ESP_K7z_G3yz_a+ABZ*I_ESP_I6z_G3yz_a;
    Double I_ESP_I6x_H2y3z_a = I_ESP_K6xy_Gy3z_a+ABY*I_ESP_I6x_Gy3z_a;
    Double I_ESP_I5xy_H2y3z_a = I_ESP_K5x2y_Gy3z_a+ABY*I_ESP_I5xy_Gy3z_a;
    Double I_ESP_I5xz_H2y3z_a = I_ESP_K5xyz_Gy3z_a+ABY*I_ESP_I5xz_Gy3z_a;
    Double I_ESP_I4x2y_H2y3z_a = I_ESP_K4x3y_Gy3z_a+ABY*I_ESP_I4x2y_Gy3z_a;
    Double I_ESP_I4xyz_H2y3z_a = I_ESP_K4x2yz_Gy3z_a+ABY*I_ESP_I4xyz_Gy3z_a;
    Double I_ESP_I4x2z_H2y3z_a = I_ESP_K4xy2z_Gy3z_a+ABY*I_ESP_I4x2z_Gy3z_a;
    Double I_ESP_I3x3y_H2y3z_a = I_ESP_K3x4y_Gy3z_a+ABY*I_ESP_I3x3y_Gy3z_a;
    Double I_ESP_I3x2yz_H2y3z_a = I_ESP_K3x3yz_Gy3z_a+ABY*I_ESP_I3x2yz_Gy3z_a;
    Double I_ESP_I3xy2z_H2y3z_a = I_ESP_K3x2y2z_Gy3z_a+ABY*I_ESP_I3xy2z_Gy3z_a;
    Double I_ESP_I3x3z_H2y3z_a = I_ESP_K3xy3z_Gy3z_a+ABY*I_ESP_I3x3z_Gy3z_a;
    Double I_ESP_I2x4y_H2y3z_a = I_ESP_K2x5y_Gy3z_a+ABY*I_ESP_I2x4y_Gy3z_a;
    Double I_ESP_I2x3yz_H2y3z_a = I_ESP_K2x4yz_Gy3z_a+ABY*I_ESP_I2x3yz_Gy3z_a;
    Double I_ESP_I2x2y2z_H2y3z_a = I_ESP_K2x3y2z_Gy3z_a+ABY*I_ESP_I2x2y2z_Gy3z_a;
    Double I_ESP_I2xy3z_H2y3z_a = I_ESP_K2x2y3z_Gy3z_a+ABY*I_ESP_I2xy3z_Gy3z_a;
    Double I_ESP_I2x4z_H2y3z_a = I_ESP_K2xy4z_Gy3z_a+ABY*I_ESP_I2x4z_Gy3z_a;
    Double I_ESP_Ix5y_H2y3z_a = I_ESP_Kx6y_Gy3z_a+ABY*I_ESP_Ix5y_Gy3z_a;
    Double I_ESP_Ix4yz_H2y3z_a = I_ESP_Kx5yz_Gy3z_a+ABY*I_ESP_Ix4yz_Gy3z_a;
    Double I_ESP_Ix3y2z_H2y3z_a = I_ESP_Kx4y2z_Gy3z_a+ABY*I_ESP_Ix3y2z_Gy3z_a;
    Double I_ESP_Ix2y3z_H2y3z_a = I_ESP_Kx3y3z_Gy3z_a+ABY*I_ESP_Ix2y3z_Gy3z_a;
    Double I_ESP_Ixy4z_H2y3z_a = I_ESP_Kx2y4z_Gy3z_a+ABY*I_ESP_Ixy4z_Gy3z_a;
    Double I_ESP_Ix5z_H2y3z_a = I_ESP_Kxy5z_Gy3z_a+ABY*I_ESP_Ix5z_Gy3z_a;
    Double I_ESP_I6y_H2y3z_a = I_ESP_K7y_Gy3z_a+ABY*I_ESP_I6y_Gy3z_a;
    Double I_ESP_I5yz_H2y3z_a = I_ESP_K6yz_Gy3z_a+ABY*I_ESP_I5yz_Gy3z_a;
    Double I_ESP_I4y2z_H2y3z_a = I_ESP_K5y2z_Gy3z_a+ABY*I_ESP_I4y2z_Gy3z_a;
    Double I_ESP_I3y3z_H2y3z_a = I_ESP_K4y3z_Gy3z_a+ABY*I_ESP_I3y3z_Gy3z_a;
    Double I_ESP_I2y4z_H2y3z_a = I_ESP_K3y4z_Gy3z_a+ABY*I_ESP_I2y4z_Gy3z_a;
    Double I_ESP_Iy5z_H2y3z_a = I_ESP_K2y5z_Gy3z_a+ABY*I_ESP_Iy5z_Gy3z_a;
    Double I_ESP_I6z_H2y3z_a = I_ESP_Ky6z_Gy3z_a+ABY*I_ESP_I6z_Gy3z_a;
    Double I_ESP_I6x_Hy4z_a = I_ESP_K6xy_G4z_a+ABY*I_ESP_I6x_G4z_a;
    Double I_ESP_I5xy_Hy4z_a = I_ESP_K5x2y_G4z_a+ABY*I_ESP_I5xy_G4z_a;
    Double I_ESP_I5xz_Hy4z_a = I_ESP_K5xyz_G4z_a+ABY*I_ESP_I5xz_G4z_a;
    Double I_ESP_I4x2y_Hy4z_a = I_ESP_K4x3y_G4z_a+ABY*I_ESP_I4x2y_G4z_a;
    Double I_ESP_I4xyz_Hy4z_a = I_ESP_K4x2yz_G4z_a+ABY*I_ESP_I4xyz_G4z_a;
    Double I_ESP_I4x2z_Hy4z_a = I_ESP_K4xy2z_G4z_a+ABY*I_ESP_I4x2z_G4z_a;
    Double I_ESP_I3x3y_Hy4z_a = I_ESP_K3x4y_G4z_a+ABY*I_ESP_I3x3y_G4z_a;
    Double I_ESP_I3x2yz_Hy4z_a = I_ESP_K3x3yz_G4z_a+ABY*I_ESP_I3x2yz_G4z_a;
    Double I_ESP_I3xy2z_Hy4z_a = I_ESP_K3x2y2z_G4z_a+ABY*I_ESP_I3xy2z_G4z_a;
    Double I_ESP_I3x3z_Hy4z_a = I_ESP_K3xy3z_G4z_a+ABY*I_ESP_I3x3z_G4z_a;
    Double I_ESP_I2x4y_Hy4z_a = I_ESP_K2x5y_G4z_a+ABY*I_ESP_I2x4y_G4z_a;
    Double I_ESP_I2x3yz_Hy4z_a = I_ESP_K2x4yz_G4z_a+ABY*I_ESP_I2x3yz_G4z_a;
    Double I_ESP_I2x2y2z_Hy4z_a = I_ESP_K2x3y2z_G4z_a+ABY*I_ESP_I2x2y2z_G4z_a;
    Double I_ESP_I2xy3z_Hy4z_a = I_ESP_K2x2y3z_G4z_a+ABY*I_ESP_I2xy3z_G4z_a;
    Double I_ESP_I2x4z_Hy4z_a = I_ESP_K2xy4z_G4z_a+ABY*I_ESP_I2x4z_G4z_a;
    Double I_ESP_Ix5y_Hy4z_a = I_ESP_Kx6y_G4z_a+ABY*I_ESP_Ix5y_G4z_a;
    Double I_ESP_Ix4yz_Hy4z_a = I_ESP_Kx5yz_G4z_a+ABY*I_ESP_Ix4yz_G4z_a;
    Double I_ESP_Ix3y2z_Hy4z_a = I_ESP_Kx4y2z_G4z_a+ABY*I_ESP_Ix3y2z_G4z_a;
    Double I_ESP_Ix2y3z_Hy4z_a = I_ESP_Kx3y3z_G4z_a+ABY*I_ESP_Ix2y3z_G4z_a;
    Double I_ESP_Ixy4z_Hy4z_a = I_ESP_Kx2y4z_G4z_a+ABY*I_ESP_Ixy4z_G4z_a;
    Double I_ESP_Ix5z_Hy4z_a = I_ESP_Kxy5z_G4z_a+ABY*I_ESP_Ix5z_G4z_a;
    Double I_ESP_I6y_Hy4z_a = I_ESP_K7y_G4z_a+ABY*I_ESP_I6y_G4z_a;
    Double I_ESP_I5yz_Hy4z_a = I_ESP_K6yz_G4z_a+ABY*I_ESP_I5yz_G4z_a;
    Double I_ESP_I4y2z_Hy4z_a = I_ESP_K5y2z_G4z_a+ABY*I_ESP_I4y2z_G4z_a;
    Double I_ESP_I3y3z_Hy4z_a = I_ESP_K4y3z_G4z_a+ABY*I_ESP_I3y3z_G4z_a;
    Double I_ESP_I2y4z_Hy4z_a = I_ESP_K3y4z_G4z_a+ABY*I_ESP_I2y4z_G4z_a;
    Double I_ESP_Iy5z_Hy4z_a = I_ESP_K2y5z_G4z_a+ABY*I_ESP_Iy5z_G4z_a;
    Double I_ESP_I6z_Hy4z_a = I_ESP_Ky6z_G4z_a+ABY*I_ESP_I6z_G4z_a;
    Double I_ESP_I6x_H5z_a = I_ESP_K6xz_G4z_a+ABZ*I_ESP_I6x_G4z_a;
    Double I_ESP_I5xy_H5z_a = I_ESP_K5xyz_G4z_a+ABZ*I_ESP_I5xy_G4z_a;
    Double I_ESP_I5xz_H5z_a = I_ESP_K5x2z_G4z_a+ABZ*I_ESP_I5xz_G4z_a;
    Double I_ESP_I4x2y_H5z_a = I_ESP_K4x2yz_G4z_a+ABZ*I_ESP_I4x2y_G4z_a;
    Double I_ESP_I4xyz_H5z_a = I_ESP_K4xy2z_G4z_a+ABZ*I_ESP_I4xyz_G4z_a;
    Double I_ESP_I4x2z_H5z_a = I_ESP_K4x3z_G4z_a+ABZ*I_ESP_I4x2z_G4z_a;
    Double I_ESP_I3x3y_H5z_a = I_ESP_K3x3yz_G4z_a+ABZ*I_ESP_I3x3y_G4z_a;
    Double I_ESP_I3x2yz_H5z_a = I_ESP_K3x2y2z_G4z_a+ABZ*I_ESP_I3x2yz_G4z_a;
    Double I_ESP_I3xy2z_H5z_a = I_ESP_K3xy3z_G4z_a+ABZ*I_ESP_I3xy2z_G4z_a;
    Double I_ESP_I3x3z_H5z_a = I_ESP_K3x4z_G4z_a+ABZ*I_ESP_I3x3z_G4z_a;
    Double I_ESP_I2x4y_H5z_a = I_ESP_K2x4yz_G4z_a+ABZ*I_ESP_I2x4y_G4z_a;
    Double I_ESP_I2x3yz_H5z_a = I_ESP_K2x3y2z_G4z_a+ABZ*I_ESP_I2x3yz_G4z_a;
    Double I_ESP_I2x2y2z_H5z_a = I_ESP_K2x2y3z_G4z_a+ABZ*I_ESP_I2x2y2z_G4z_a;
    Double I_ESP_I2xy3z_H5z_a = I_ESP_K2xy4z_G4z_a+ABZ*I_ESP_I2xy3z_G4z_a;
    Double I_ESP_I2x4z_H5z_a = I_ESP_K2x5z_G4z_a+ABZ*I_ESP_I2x4z_G4z_a;
    Double I_ESP_Ix5y_H5z_a = I_ESP_Kx5yz_G4z_a+ABZ*I_ESP_Ix5y_G4z_a;
    Double I_ESP_Ix4yz_H5z_a = I_ESP_Kx4y2z_G4z_a+ABZ*I_ESP_Ix4yz_G4z_a;
    Double I_ESP_Ix3y2z_H5z_a = I_ESP_Kx3y3z_G4z_a+ABZ*I_ESP_Ix3y2z_G4z_a;
    Double I_ESP_Ix2y3z_H5z_a = I_ESP_Kx2y4z_G4z_a+ABZ*I_ESP_Ix2y3z_G4z_a;
    Double I_ESP_Ixy4z_H5z_a = I_ESP_Kxy5z_G4z_a+ABZ*I_ESP_Ixy4z_G4z_a;
    Double I_ESP_Ix5z_H5z_a = I_ESP_Kx6z_G4z_a+ABZ*I_ESP_Ix5z_G4z_a;
    Double I_ESP_I6y_H5z_a = I_ESP_K6yz_G4z_a+ABZ*I_ESP_I6y_G4z_a;
    Double I_ESP_I5yz_H5z_a = I_ESP_K5y2z_G4z_a+ABZ*I_ESP_I5yz_G4z_a;
    Double I_ESP_I4y2z_H5z_a = I_ESP_K4y3z_G4z_a+ABZ*I_ESP_I4y2z_G4z_a;
    Double I_ESP_I3y3z_H5z_a = I_ESP_K3y4z_G4z_a+ABZ*I_ESP_I3y3z_G4z_a;
    Double I_ESP_I2y4z_H5z_a = I_ESP_K2y5z_G4z_a+ABZ*I_ESP_I2y4z_G4z_a;
    Double I_ESP_Iy5z_H5z_a = I_ESP_Ky6z_G4z_a+ABZ*I_ESP_Iy5z_G4z_a;
    Double I_ESP_I6z_H5z_a = I_ESP_K7z_G4z_a+ABZ*I_ESP_I6z_G4z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_H_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_H_a
     * RHS shell quartet name: SQ_ESP_G_H
     ************************************************************/
    abcd[iGrid*1323+0] = 2.0E0*I_ESP_I6x_H5x_a-5*I_ESP_G4x_H5x;
    abcd[iGrid*1323+1] = 2.0E0*I_ESP_I5xy_H5x_a-4*I_ESP_G3xy_H5x;
    abcd[iGrid*1323+2] = 2.0E0*I_ESP_I5xz_H5x_a-4*I_ESP_G3xz_H5x;
    abcd[iGrid*1323+3] = 2.0E0*I_ESP_I4x2y_H5x_a-3*I_ESP_G2x2y_H5x;
    abcd[iGrid*1323+4] = 2.0E0*I_ESP_I4xyz_H5x_a-3*I_ESP_G2xyz_H5x;
    abcd[iGrid*1323+5] = 2.0E0*I_ESP_I4x2z_H5x_a-3*I_ESP_G2x2z_H5x;
    abcd[iGrid*1323+6] = 2.0E0*I_ESP_I3x3y_H5x_a-2*I_ESP_Gx3y_H5x;
    abcd[iGrid*1323+7] = 2.0E0*I_ESP_I3x2yz_H5x_a-2*I_ESP_Gx2yz_H5x;
    abcd[iGrid*1323+8] = 2.0E0*I_ESP_I3xy2z_H5x_a-2*I_ESP_Gxy2z_H5x;
    abcd[iGrid*1323+9] = 2.0E0*I_ESP_I3x3z_H5x_a-2*I_ESP_Gx3z_H5x;
    abcd[iGrid*1323+10] = 2.0E0*I_ESP_I2x4y_H5x_a-1*I_ESP_G4y_H5x;
    abcd[iGrid*1323+11] = 2.0E0*I_ESP_I2x3yz_H5x_a-1*I_ESP_G3yz_H5x;
    abcd[iGrid*1323+12] = 2.0E0*I_ESP_I2x2y2z_H5x_a-1*I_ESP_G2y2z_H5x;
    abcd[iGrid*1323+13] = 2.0E0*I_ESP_I2xy3z_H5x_a-1*I_ESP_Gy3z_H5x;
    abcd[iGrid*1323+14] = 2.0E0*I_ESP_I2x4z_H5x_a-1*I_ESP_G4z_H5x;
    abcd[iGrid*1323+15] = 2.0E0*I_ESP_Ix5y_H5x_a;
    abcd[iGrid*1323+16] = 2.0E0*I_ESP_Ix4yz_H5x_a;
    abcd[iGrid*1323+17] = 2.0E0*I_ESP_Ix3y2z_H5x_a;
    abcd[iGrid*1323+18] = 2.0E0*I_ESP_Ix2y3z_H5x_a;
    abcd[iGrid*1323+19] = 2.0E0*I_ESP_Ixy4z_H5x_a;
    abcd[iGrid*1323+20] = 2.0E0*I_ESP_Ix5z_H5x_a;
    abcd[iGrid*1323+21] = 2.0E0*I_ESP_I6x_H4xy_a-5*I_ESP_G4x_H4xy;
    abcd[iGrid*1323+22] = 2.0E0*I_ESP_I5xy_H4xy_a-4*I_ESP_G3xy_H4xy;
    abcd[iGrid*1323+23] = 2.0E0*I_ESP_I5xz_H4xy_a-4*I_ESP_G3xz_H4xy;
    abcd[iGrid*1323+24] = 2.0E0*I_ESP_I4x2y_H4xy_a-3*I_ESP_G2x2y_H4xy;
    abcd[iGrid*1323+25] = 2.0E0*I_ESP_I4xyz_H4xy_a-3*I_ESP_G2xyz_H4xy;
    abcd[iGrid*1323+26] = 2.0E0*I_ESP_I4x2z_H4xy_a-3*I_ESP_G2x2z_H4xy;
    abcd[iGrid*1323+27] = 2.0E0*I_ESP_I3x3y_H4xy_a-2*I_ESP_Gx3y_H4xy;
    abcd[iGrid*1323+28] = 2.0E0*I_ESP_I3x2yz_H4xy_a-2*I_ESP_Gx2yz_H4xy;
    abcd[iGrid*1323+29] = 2.0E0*I_ESP_I3xy2z_H4xy_a-2*I_ESP_Gxy2z_H4xy;
    abcd[iGrid*1323+30] = 2.0E0*I_ESP_I3x3z_H4xy_a-2*I_ESP_Gx3z_H4xy;
    abcd[iGrid*1323+31] = 2.0E0*I_ESP_I2x4y_H4xy_a-1*I_ESP_G4y_H4xy;
    abcd[iGrid*1323+32] = 2.0E0*I_ESP_I2x3yz_H4xy_a-1*I_ESP_G3yz_H4xy;
    abcd[iGrid*1323+33] = 2.0E0*I_ESP_I2x2y2z_H4xy_a-1*I_ESP_G2y2z_H4xy;
    abcd[iGrid*1323+34] = 2.0E0*I_ESP_I2xy3z_H4xy_a-1*I_ESP_Gy3z_H4xy;
    abcd[iGrid*1323+35] = 2.0E0*I_ESP_I2x4z_H4xy_a-1*I_ESP_G4z_H4xy;
    abcd[iGrid*1323+36] = 2.0E0*I_ESP_Ix5y_H4xy_a;
    abcd[iGrid*1323+37] = 2.0E0*I_ESP_Ix4yz_H4xy_a;
    abcd[iGrid*1323+38] = 2.0E0*I_ESP_Ix3y2z_H4xy_a;
    abcd[iGrid*1323+39] = 2.0E0*I_ESP_Ix2y3z_H4xy_a;
    abcd[iGrid*1323+40] = 2.0E0*I_ESP_Ixy4z_H4xy_a;
    abcd[iGrid*1323+41] = 2.0E0*I_ESP_Ix5z_H4xy_a;
    abcd[iGrid*1323+42] = 2.0E0*I_ESP_I6x_H4xz_a-5*I_ESP_G4x_H4xz;
    abcd[iGrid*1323+43] = 2.0E0*I_ESP_I5xy_H4xz_a-4*I_ESP_G3xy_H4xz;
    abcd[iGrid*1323+44] = 2.0E0*I_ESP_I5xz_H4xz_a-4*I_ESP_G3xz_H4xz;
    abcd[iGrid*1323+45] = 2.0E0*I_ESP_I4x2y_H4xz_a-3*I_ESP_G2x2y_H4xz;
    abcd[iGrid*1323+46] = 2.0E0*I_ESP_I4xyz_H4xz_a-3*I_ESP_G2xyz_H4xz;
    abcd[iGrid*1323+47] = 2.0E0*I_ESP_I4x2z_H4xz_a-3*I_ESP_G2x2z_H4xz;
    abcd[iGrid*1323+48] = 2.0E0*I_ESP_I3x3y_H4xz_a-2*I_ESP_Gx3y_H4xz;
    abcd[iGrid*1323+49] = 2.0E0*I_ESP_I3x2yz_H4xz_a-2*I_ESP_Gx2yz_H4xz;
    abcd[iGrid*1323+50] = 2.0E0*I_ESP_I3xy2z_H4xz_a-2*I_ESP_Gxy2z_H4xz;
    abcd[iGrid*1323+51] = 2.0E0*I_ESP_I3x3z_H4xz_a-2*I_ESP_Gx3z_H4xz;
    abcd[iGrid*1323+52] = 2.0E0*I_ESP_I2x4y_H4xz_a-1*I_ESP_G4y_H4xz;
    abcd[iGrid*1323+53] = 2.0E0*I_ESP_I2x3yz_H4xz_a-1*I_ESP_G3yz_H4xz;
    abcd[iGrid*1323+54] = 2.0E0*I_ESP_I2x2y2z_H4xz_a-1*I_ESP_G2y2z_H4xz;
    abcd[iGrid*1323+55] = 2.0E0*I_ESP_I2xy3z_H4xz_a-1*I_ESP_Gy3z_H4xz;
    abcd[iGrid*1323+56] = 2.0E0*I_ESP_I2x4z_H4xz_a-1*I_ESP_G4z_H4xz;
    abcd[iGrid*1323+57] = 2.0E0*I_ESP_Ix5y_H4xz_a;
    abcd[iGrid*1323+58] = 2.0E0*I_ESP_Ix4yz_H4xz_a;
    abcd[iGrid*1323+59] = 2.0E0*I_ESP_Ix3y2z_H4xz_a;
    abcd[iGrid*1323+60] = 2.0E0*I_ESP_Ix2y3z_H4xz_a;
    abcd[iGrid*1323+61] = 2.0E0*I_ESP_Ixy4z_H4xz_a;
    abcd[iGrid*1323+62] = 2.0E0*I_ESP_Ix5z_H4xz_a;
    abcd[iGrid*1323+63] = 2.0E0*I_ESP_I6x_H3x2y_a-5*I_ESP_G4x_H3x2y;
    abcd[iGrid*1323+64] = 2.0E0*I_ESP_I5xy_H3x2y_a-4*I_ESP_G3xy_H3x2y;
    abcd[iGrid*1323+65] = 2.0E0*I_ESP_I5xz_H3x2y_a-4*I_ESP_G3xz_H3x2y;
    abcd[iGrid*1323+66] = 2.0E0*I_ESP_I4x2y_H3x2y_a-3*I_ESP_G2x2y_H3x2y;
    abcd[iGrid*1323+67] = 2.0E0*I_ESP_I4xyz_H3x2y_a-3*I_ESP_G2xyz_H3x2y;
    abcd[iGrid*1323+68] = 2.0E0*I_ESP_I4x2z_H3x2y_a-3*I_ESP_G2x2z_H3x2y;
    abcd[iGrid*1323+69] = 2.0E0*I_ESP_I3x3y_H3x2y_a-2*I_ESP_Gx3y_H3x2y;
    abcd[iGrid*1323+70] = 2.0E0*I_ESP_I3x2yz_H3x2y_a-2*I_ESP_Gx2yz_H3x2y;
    abcd[iGrid*1323+71] = 2.0E0*I_ESP_I3xy2z_H3x2y_a-2*I_ESP_Gxy2z_H3x2y;
    abcd[iGrid*1323+72] = 2.0E0*I_ESP_I3x3z_H3x2y_a-2*I_ESP_Gx3z_H3x2y;
    abcd[iGrid*1323+73] = 2.0E0*I_ESP_I2x4y_H3x2y_a-1*I_ESP_G4y_H3x2y;
    abcd[iGrid*1323+74] = 2.0E0*I_ESP_I2x3yz_H3x2y_a-1*I_ESP_G3yz_H3x2y;
    abcd[iGrid*1323+75] = 2.0E0*I_ESP_I2x2y2z_H3x2y_a-1*I_ESP_G2y2z_H3x2y;
    abcd[iGrid*1323+76] = 2.0E0*I_ESP_I2xy3z_H3x2y_a-1*I_ESP_Gy3z_H3x2y;
    abcd[iGrid*1323+77] = 2.0E0*I_ESP_I2x4z_H3x2y_a-1*I_ESP_G4z_H3x2y;
    abcd[iGrid*1323+78] = 2.0E0*I_ESP_Ix5y_H3x2y_a;
    abcd[iGrid*1323+79] = 2.0E0*I_ESP_Ix4yz_H3x2y_a;
    abcd[iGrid*1323+80] = 2.0E0*I_ESP_Ix3y2z_H3x2y_a;
    abcd[iGrid*1323+81] = 2.0E0*I_ESP_Ix2y3z_H3x2y_a;
    abcd[iGrid*1323+82] = 2.0E0*I_ESP_Ixy4z_H3x2y_a;
    abcd[iGrid*1323+83] = 2.0E0*I_ESP_Ix5z_H3x2y_a;
    abcd[iGrid*1323+84] = 2.0E0*I_ESP_I6x_H3xyz_a-5*I_ESP_G4x_H3xyz;
    abcd[iGrid*1323+85] = 2.0E0*I_ESP_I5xy_H3xyz_a-4*I_ESP_G3xy_H3xyz;
    abcd[iGrid*1323+86] = 2.0E0*I_ESP_I5xz_H3xyz_a-4*I_ESP_G3xz_H3xyz;
    abcd[iGrid*1323+87] = 2.0E0*I_ESP_I4x2y_H3xyz_a-3*I_ESP_G2x2y_H3xyz;
    abcd[iGrid*1323+88] = 2.0E0*I_ESP_I4xyz_H3xyz_a-3*I_ESP_G2xyz_H3xyz;
    abcd[iGrid*1323+89] = 2.0E0*I_ESP_I4x2z_H3xyz_a-3*I_ESP_G2x2z_H3xyz;
    abcd[iGrid*1323+90] = 2.0E0*I_ESP_I3x3y_H3xyz_a-2*I_ESP_Gx3y_H3xyz;
    abcd[iGrid*1323+91] = 2.0E0*I_ESP_I3x2yz_H3xyz_a-2*I_ESP_Gx2yz_H3xyz;
    abcd[iGrid*1323+92] = 2.0E0*I_ESP_I3xy2z_H3xyz_a-2*I_ESP_Gxy2z_H3xyz;
    abcd[iGrid*1323+93] = 2.0E0*I_ESP_I3x3z_H3xyz_a-2*I_ESP_Gx3z_H3xyz;
    abcd[iGrid*1323+94] = 2.0E0*I_ESP_I2x4y_H3xyz_a-1*I_ESP_G4y_H3xyz;
    abcd[iGrid*1323+95] = 2.0E0*I_ESP_I2x3yz_H3xyz_a-1*I_ESP_G3yz_H3xyz;
    abcd[iGrid*1323+96] = 2.0E0*I_ESP_I2x2y2z_H3xyz_a-1*I_ESP_G2y2z_H3xyz;
    abcd[iGrid*1323+97] = 2.0E0*I_ESP_I2xy3z_H3xyz_a-1*I_ESP_Gy3z_H3xyz;
    abcd[iGrid*1323+98] = 2.0E0*I_ESP_I2x4z_H3xyz_a-1*I_ESP_G4z_H3xyz;
    abcd[iGrid*1323+99] = 2.0E0*I_ESP_Ix5y_H3xyz_a;
    abcd[iGrid*1323+100] = 2.0E0*I_ESP_Ix4yz_H3xyz_a;
    abcd[iGrid*1323+101] = 2.0E0*I_ESP_Ix3y2z_H3xyz_a;
    abcd[iGrid*1323+102] = 2.0E0*I_ESP_Ix2y3z_H3xyz_a;
    abcd[iGrid*1323+103] = 2.0E0*I_ESP_Ixy4z_H3xyz_a;
    abcd[iGrid*1323+104] = 2.0E0*I_ESP_Ix5z_H3xyz_a;
    abcd[iGrid*1323+105] = 2.0E0*I_ESP_I6x_H3x2z_a-5*I_ESP_G4x_H3x2z;
    abcd[iGrid*1323+106] = 2.0E0*I_ESP_I5xy_H3x2z_a-4*I_ESP_G3xy_H3x2z;
    abcd[iGrid*1323+107] = 2.0E0*I_ESP_I5xz_H3x2z_a-4*I_ESP_G3xz_H3x2z;
    abcd[iGrid*1323+108] = 2.0E0*I_ESP_I4x2y_H3x2z_a-3*I_ESP_G2x2y_H3x2z;
    abcd[iGrid*1323+109] = 2.0E0*I_ESP_I4xyz_H3x2z_a-3*I_ESP_G2xyz_H3x2z;
    abcd[iGrid*1323+110] = 2.0E0*I_ESP_I4x2z_H3x2z_a-3*I_ESP_G2x2z_H3x2z;
    abcd[iGrid*1323+111] = 2.0E0*I_ESP_I3x3y_H3x2z_a-2*I_ESP_Gx3y_H3x2z;
    abcd[iGrid*1323+112] = 2.0E0*I_ESP_I3x2yz_H3x2z_a-2*I_ESP_Gx2yz_H3x2z;
    abcd[iGrid*1323+113] = 2.0E0*I_ESP_I3xy2z_H3x2z_a-2*I_ESP_Gxy2z_H3x2z;
    abcd[iGrid*1323+114] = 2.0E0*I_ESP_I3x3z_H3x2z_a-2*I_ESP_Gx3z_H3x2z;
    abcd[iGrid*1323+115] = 2.0E0*I_ESP_I2x4y_H3x2z_a-1*I_ESP_G4y_H3x2z;
    abcd[iGrid*1323+116] = 2.0E0*I_ESP_I2x3yz_H3x2z_a-1*I_ESP_G3yz_H3x2z;
    abcd[iGrid*1323+117] = 2.0E0*I_ESP_I2x2y2z_H3x2z_a-1*I_ESP_G2y2z_H3x2z;
    abcd[iGrid*1323+118] = 2.0E0*I_ESP_I2xy3z_H3x2z_a-1*I_ESP_Gy3z_H3x2z;
    abcd[iGrid*1323+119] = 2.0E0*I_ESP_I2x4z_H3x2z_a-1*I_ESP_G4z_H3x2z;
    abcd[iGrid*1323+120] = 2.0E0*I_ESP_Ix5y_H3x2z_a;
    abcd[iGrid*1323+121] = 2.0E0*I_ESP_Ix4yz_H3x2z_a;
    abcd[iGrid*1323+122] = 2.0E0*I_ESP_Ix3y2z_H3x2z_a;
    abcd[iGrid*1323+123] = 2.0E0*I_ESP_Ix2y3z_H3x2z_a;
    abcd[iGrid*1323+124] = 2.0E0*I_ESP_Ixy4z_H3x2z_a;
    abcd[iGrid*1323+125] = 2.0E0*I_ESP_Ix5z_H3x2z_a;
    abcd[iGrid*1323+126] = 2.0E0*I_ESP_I6x_H2x3y_a-5*I_ESP_G4x_H2x3y;
    abcd[iGrid*1323+127] = 2.0E0*I_ESP_I5xy_H2x3y_a-4*I_ESP_G3xy_H2x3y;
    abcd[iGrid*1323+128] = 2.0E0*I_ESP_I5xz_H2x3y_a-4*I_ESP_G3xz_H2x3y;
    abcd[iGrid*1323+129] = 2.0E0*I_ESP_I4x2y_H2x3y_a-3*I_ESP_G2x2y_H2x3y;
    abcd[iGrid*1323+130] = 2.0E0*I_ESP_I4xyz_H2x3y_a-3*I_ESP_G2xyz_H2x3y;
    abcd[iGrid*1323+131] = 2.0E0*I_ESP_I4x2z_H2x3y_a-3*I_ESP_G2x2z_H2x3y;
    abcd[iGrid*1323+132] = 2.0E0*I_ESP_I3x3y_H2x3y_a-2*I_ESP_Gx3y_H2x3y;
    abcd[iGrid*1323+133] = 2.0E0*I_ESP_I3x2yz_H2x3y_a-2*I_ESP_Gx2yz_H2x3y;
    abcd[iGrid*1323+134] = 2.0E0*I_ESP_I3xy2z_H2x3y_a-2*I_ESP_Gxy2z_H2x3y;
    abcd[iGrid*1323+135] = 2.0E0*I_ESP_I3x3z_H2x3y_a-2*I_ESP_Gx3z_H2x3y;
    abcd[iGrid*1323+136] = 2.0E0*I_ESP_I2x4y_H2x3y_a-1*I_ESP_G4y_H2x3y;
    abcd[iGrid*1323+137] = 2.0E0*I_ESP_I2x3yz_H2x3y_a-1*I_ESP_G3yz_H2x3y;
    abcd[iGrid*1323+138] = 2.0E0*I_ESP_I2x2y2z_H2x3y_a-1*I_ESP_G2y2z_H2x3y;
    abcd[iGrid*1323+139] = 2.0E0*I_ESP_I2xy3z_H2x3y_a-1*I_ESP_Gy3z_H2x3y;
    abcd[iGrid*1323+140] = 2.0E0*I_ESP_I2x4z_H2x3y_a-1*I_ESP_G4z_H2x3y;
    abcd[iGrid*1323+141] = 2.0E0*I_ESP_Ix5y_H2x3y_a;
    abcd[iGrid*1323+142] = 2.0E0*I_ESP_Ix4yz_H2x3y_a;
    abcd[iGrid*1323+143] = 2.0E0*I_ESP_Ix3y2z_H2x3y_a;
    abcd[iGrid*1323+144] = 2.0E0*I_ESP_Ix2y3z_H2x3y_a;
    abcd[iGrid*1323+145] = 2.0E0*I_ESP_Ixy4z_H2x3y_a;
    abcd[iGrid*1323+146] = 2.0E0*I_ESP_Ix5z_H2x3y_a;
    abcd[iGrid*1323+147] = 2.0E0*I_ESP_I6x_H2x2yz_a-5*I_ESP_G4x_H2x2yz;
    abcd[iGrid*1323+148] = 2.0E0*I_ESP_I5xy_H2x2yz_a-4*I_ESP_G3xy_H2x2yz;
    abcd[iGrid*1323+149] = 2.0E0*I_ESP_I5xz_H2x2yz_a-4*I_ESP_G3xz_H2x2yz;
    abcd[iGrid*1323+150] = 2.0E0*I_ESP_I4x2y_H2x2yz_a-3*I_ESP_G2x2y_H2x2yz;
    abcd[iGrid*1323+151] = 2.0E0*I_ESP_I4xyz_H2x2yz_a-3*I_ESP_G2xyz_H2x2yz;
    abcd[iGrid*1323+152] = 2.0E0*I_ESP_I4x2z_H2x2yz_a-3*I_ESP_G2x2z_H2x2yz;
    abcd[iGrid*1323+153] = 2.0E0*I_ESP_I3x3y_H2x2yz_a-2*I_ESP_Gx3y_H2x2yz;
    abcd[iGrid*1323+154] = 2.0E0*I_ESP_I3x2yz_H2x2yz_a-2*I_ESP_Gx2yz_H2x2yz;
    abcd[iGrid*1323+155] = 2.0E0*I_ESP_I3xy2z_H2x2yz_a-2*I_ESP_Gxy2z_H2x2yz;
    abcd[iGrid*1323+156] = 2.0E0*I_ESP_I3x3z_H2x2yz_a-2*I_ESP_Gx3z_H2x2yz;
    abcd[iGrid*1323+157] = 2.0E0*I_ESP_I2x4y_H2x2yz_a-1*I_ESP_G4y_H2x2yz;
    abcd[iGrid*1323+158] = 2.0E0*I_ESP_I2x3yz_H2x2yz_a-1*I_ESP_G3yz_H2x2yz;
    abcd[iGrid*1323+159] = 2.0E0*I_ESP_I2x2y2z_H2x2yz_a-1*I_ESP_G2y2z_H2x2yz;
    abcd[iGrid*1323+160] = 2.0E0*I_ESP_I2xy3z_H2x2yz_a-1*I_ESP_Gy3z_H2x2yz;
    abcd[iGrid*1323+161] = 2.0E0*I_ESP_I2x4z_H2x2yz_a-1*I_ESP_G4z_H2x2yz;
    abcd[iGrid*1323+162] = 2.0E0*I_ESP_Ix5y_H2x2yz_a;
    abcd[iGrid*1323+163] = 2.0E0*I_ESP_Ix4yz_H2x2yz_a;
    abcd[iGrid*1323+164] = 2.0E0*I_ESP_Ix3y2z_H2x2yz_a;
    abcd[iGrid*1323+165] = 2.0E0*I_ESP_Ix2y3z_H2x2yz_a;
    abcd[iGrid*1323+166] = 2.0E0*I_ESP_Ixy4z_H2x2yz_a;
    abcd[iGrid*1323+167] = 2.0E0*I_ESP_Ix5z_H2x2yz_a;
    abcd[iGrid*1323+168] = 2.0E0*I_ESP_I6x_H2xy2z_a-5*I_ESP_G4x_H2xy2z;
    abcd[iGrid*1323+169] = 2.0E0*I_ESP_I5xy_H2xy2z_a-4*I_ESP_G3xy_H2xy2z;
    abcd[iGrid*1323+170] = 2.0E0*I_ESP_I5xz_H2xy2z_a-4*I_ESP_G3xz_H2xy2z;
    abcd[iGrid*1323+171] = 2.0E0*I_ESP_I4x2y_H2xy2z_a-3*I_ESP_G2x2y_H2xy2z;
    abcd[iGrid*1323+172] = 2.0E0*I_ESP_I4xyz_H2xy2z_a-3*I_ESP_G2xyz_H2xy2z;
    abcd[iGrid*1323+173] = 2.0E0*I_ESP_I4x2z_H2xy2z_a-3*I_ESP_G2x2z_H2xy2z;
    abcd[iGrid*1323+174] = 2.0E0*I_ESP_I3x3y_H2xy2z_a-2*I_ESP_Gx3y_H2xy2z;
    abcd[iGrid*1323+175] = 2.0E0*I_ESP_I3x2yz_H2xy2z_a-2*I_ESP_Gx2yz_H2xy2z;
    abcd[iGrid*1323+176] = 2.0E0*I_ESP_I3xy2z_H2xy2z_a-2*I_ESP_Gxy2z_H2xy2z;
    abcd[iGrid*1323+177] = 2.0E0*I_ESP_I3x3z_H2xy2z_a-2*I_ESP_Gx3z_H2xy2z;
    abcd[iGrid*1323+178] = 2.0E0*I_ESP_I2x4y_H2xy2z_a-1*I_ESP_G4y_H2xy2z;
    abcd[iGrid*1323+179] = 2.0E0*I_ESP_I2x3yz_H2xy2z_a-1*I_ESP_G3yz_H2xy2z;
    abcd[iGrid*1323+180] = 2.0E0*I_ESP_I2x2y2z_H2xy2z_a-1*I_ESP_G2y2z_H2xy2z;
    abcd[iGrid*1323+181] = 2.0E0*I_ESP_I2xy3z_H2xy2z_a-1*I_ESP_Gy3z_H2xy2z;
    abcd[iGrid*1323+182] = 2.0E0*I_ESP_I2x4z_H2xy2z_a-1*I_ESP_G4z_H2xy2z;
    abcd[iGrid*1323+183] = 2.0E0*I_ESP_Ix5y_H2xy2z_a;
    abcd[iGrid*1323+184] = 2.0E0*I_ESP_Ix4yz_H2xy2z_a;
    abcd[iGrid*1323+185] = 2.0E0*I_ESP_Ix3y2z_H2xy2z_a;
    abcd[iGrid*1323+186] = 2.0E0*I_ESP_Ix2y3z_H2xy2z_a;
    abcd[iGrid*1323+187] = 2.0E0*I_ESP_Ixy4z_H2xy2z_a;
    abcd[iGrid*1323+188] = 2.0E0*I_ESP_Ix5z_H2xy2z_a;
    abcd[iGrid*1323+189] = 2.0E0*I_ESP_I6x_H2x3z_a-5*I_ESP_G4x_H2x3z;
    abcd[iGrid*1323+190] = 2.0E0*I_ESP_I5xy_H2x3z_a-4*I_ESP_G3xy_H2x3z;
    abcd[iGrid*1323+191] = 2.0E0*I_ESP_I5xz_H2x3z_a-4*I_ESP_G3xz_H2x3z;
    abcd[iGrid*1323+192] = 2.0E0*I_ESP_I4x2y_H2x3z_a-3*I_ESP_G2x2y_H2x3z;
    abcd[iGrid*1323+193] = 2.0E0*I_ESP_I4xyz_H2x3z_a-3*I_ESP_G2xyz_H2x3z;
    abcd[iGrid*1323+194] = 2.0E0*I_ESP_I4x2z_H2x3z_a-3*I_ESP_G2x2z_H2x3z;
    abcd[iGrid*1323+195] = 2.0E0*I_ESP_I3x3y_H2x3z_a-2*I_ESP_Gx3y_H2x3z;
    abcd[iGrid*1323+196] = 2.0E0*I_ESP_I3x2yz_H2x3z_a-2*I_ESP_Gx2yz_H2x3z;
    abcd[iGrid*1323+197] = 2.0E0*I_ESP_I3xy2z_H2x3z_a-2*I_ESP_Gxy2z_H2x3z;
    abcd[iGrid*1323+198] = 2.0E0*I_ESP_I3x3z_H2x3z_a-2*I_ESP_Gx3z_H2x3z;
    abcd[iGrid*1323+199] = 2.0E0*I_ESP_I2x4y_H2x3z_a-1*I_ESP_G4y_H2x3z;
    abcd[iGrid*1323+200] = 2.0E0*I_ESP_I2x3yz_H2x3z_a-1*I_ESP_G3yz_H2x3z;
    abcd[iGrid*1323+201] = 2.0E0*I_ESP_I2x2y2z_H2x3z_a-1*I_ESP_G2y2z_H2x3z;
    abcd[iGrid*1323+202] = 2.0E0*I_ESP_I2xy3z_H2x3z_a-1*I_ESP_Gy3z_H2x3z;
    abcd[iGrid*1323+203] = 2.0E0*I_ESP_I2x4z_H2x3z_a-1*I_ESP_G4z_H2x3z;
    abcd[iGrid*1323+204] = 2.0E0*I_ESP_Ix5y_H2x3z_a;
    abcd[iGrid*1323+205] = 2.0E0*I_ESP_Ix4yz_H2x3z_a;
    abcd[iGrid*1323+206] = 2.0E0*I_ESP_Ix3y2z_H2x3z_a;
    abcd[iGrid*1323+207] = 2.0E0*I_ESP_Ix2y3z_H2x3z_a;
    abcd[iGrid*1323+208] = 2.0E0*I_ESP_Ixy4z_H2x3z_a;
    abcd[iGrid*1323+209] = 2.0E0*I_ESP_Ix5z_H2x3z_a;
    abcd[iGrid*1323+210] = 2.0E0*I_ESP_I6x_Hx4y_a-5*I_ESP_G4x_Hx4y;
    abcd[iGrid*1323+211] = 2.0E0*I_ESP_I5xy_Hx4y_a-4*I_ESP_G3xy_Hx4y;
    abcd[iGrid*1323+212] = 2.0E0*I_ESP_I5xz_Hx4y_a-4*I_ESP_G3xz_Hx4y;
    abcd[iGrid*1323+213] = 2.0E0*I_ESP_I4x2y_Hx4y_a-3*I_ESP_G2x2y_Hx4y;
    abcd[iGrid*1323+214] = 2.0E0*I_ESP_I4xyz_Hx4y_a-3*I_ESP_G2xyz_Hx4y;
    abcd[iGrid*1323+215] = 2.0E0*I_ESP_I4x2z_Hx4y_a-3*I_ESP_G2x2z_Hx4y;
    abcd[iGrid*1323+216] = 2.0E0*I_ESP_I3x3y_Hx4y_a-2*I_ESP_Gx3y_Hx4y;
    abcd[iGrid*1323+217] = 2.0E0*I_ESP_I3x2yz_Hx4y_a-2*I_ESP_Gx2yz_Hx4y;
    abcd[iGrid*1323+218] = 2.0E0*I_ESP_I3xy2z_Hx4y_a-2*I_ESP_Gxy2z_Hx4y;
    abcd[iGrid*1323+219] = 2.0E0*I_ESP_I3x3z_Hx4y_a-2*I_ESP_Gx3z_Hx4y;
    abcd[iGrid*1323+220] = 2.0E0*I_ESP_I2x4y_Hx4y_a-1*I_ESP_G4y_Hx4y;
    abcd[iGrid*1323+221] = 2.0E0*I_ESP_I2x3yz_Hx4y_a-1*I_ESP_G3yz_Hx4y;
    abcd[iGrid*1323+222] = 2.0E0*I_ESP_I2x2y2z_Hx4y_a-1*I_ESP_G2y2z_Hx4y;
    abcd[iGrid*1323+223] = 2.0E0*I_ESP_I2xy3z_Hx4y_a-1*I_ESP_Gy3z_Hx4y;
    abcd[iGrid*1323+224] = 2.0E0*I_ESP_I2x4z_Hx4y_a-1*I_ESP_G4z_Hx4y;
    abcd[iGrid*1323+225] = 2.0E0*I_ESP_Ix5y_Hx4y_a;
    abcd[iGrid*1323+226] = 2.0E0*I_ESP_Ix4yz_Hx4y_a;
    abcd[iGrid*1323+227] = 2.0E0*I_ESP_Ix3y2z_Hx4y_a;
    abcd[iGrid*1323+228] = 2.0E0*I_ESP_Ix2y3z_Hx4y_a;
    abcd[iGrid*1323+229] = 2.0E0*I_ESP_Ixy4z_Hx4y_a;
    abcd[iGrid*1323+230] = 2.0E0*I_ESP_Ix5z_Hx4y_a;
    abcd[iGrid*1323+231] = 2.0E0*I_ESP_I6x_Hx3yz_a-5*I_ESP_G4x_Hx3yz;
    abcd[iGrid*1323+232] = 2.0E0*I_ESP_I5xy_Hx3yz_a-4*I_ESP_G3xy_Hx3yz;
    abcd[iGrid*1323+233] = 2.0E0*I_ESP_I5xz_Hx3yz_a-4*I_ESP_G3xz_Hx3yz;
    abcd[iGrid*1323+234] = 2.0E0*I_ESP_I4x2y_Hx3yz_a-3*I_ESP_G2x2y_Hx3yz;
    abcd[iGrid*1323+235] = 2.0E0*I_ESP_I4xyz_Hx3yz_a-3*I_ESP_G2xyz_Hx3yz;
    abcd[iGrid*1323+236] = 2.0E0*I_ESP_I4x2z_Hx3yz_a-3*I_ESP_G2x2z_Hx3yz;
    abcd[iGrid*1323+237] = 2.0E0*I_ESP_I3x3y_Hx3yz_a-2*I_ESP_Gx3y_Hx3yz;
    abcd[iGrid*1323+238] = 2.0E0*I_ESP_I3x2yz_Hx3yz_a-2*I_ESP_Gx2yz_Hx3yz;
    abcd[iGrid*1323+239] = 2.0E0*I_ESP_I3xy2z_Hx3yz_a-2*I_ESP_Gxy2z_Hx3yz;
    abcd[iGrid*1323+240] = 2.0E0*I_ESP_I3x3z_Hx3yz_a-2*I_ESP_Gx3z_Hx3yz;
    abcd[iGrid*1323+241] = 2.0E0*I_ESP_I2x4y_Hx3yz_a-1*I_ESP_G4y_Hx3yz;
    abcd[iGrid*1323+242] = 2.0E0*I_ESP_I2x3yz_Hx3yz_a-1*I_ESP_G3yz_Hx3yz;
    abcd[iGrid*1323+243] = 2.0E0*I_ESP_I2x2y2z_Hx3yz_a-1*I_ESP_G2y2z_Hx3yz;
    abcd[iGrid*1323+244] = 2.0E0*I_ESP_I2xy3z_Hx3yz_a-1*I_ESP_Gy3z_Hx3yz;
    abcd[iGrid*1323+245] = 2.0E0*I_ESP_I2x4z_Hx3yz_a-1*I_ESP_G4z_Hx3yz;
    abcd[iGrid*1323+246] = 2.0E0*I_ESP_Ix5y_Hx3yz_a;
    abcd[iGrid*1323+247] = 2.0E0*I_ESP_Ix4yz_Hx3yz_a;
    abcd[iGrid*1323+248] = 2.0E0*I_ESP_Ix3y2z_Hx3yz_a;
    abcd[iGrid*1323+249] = 2.0E0*I_ESP_Ix2y3z_Hx3yz_a;
    abcd[iGrid*1323+250] = 2.0E0*I_ESP_Ixy4z_Hx3yz_a;
    abcd[iGrid*1323+251] = 2.0E0*I_ESP_Ix5z_Hx3yz_a;
    abcd[iGrid*1323+252] = 2.0E0*I_ESP_I6x_Hx2y2z_a-5*I_ESP_G4x_Hx2y2z;
    abcd[iGrid*1323+253] = 2.0E0*I_ESP_I5xy_Hx2y2z_a-4*I_ESP_G3xy_Hx2y2z;
    abcd[iGrid*1323+254] = 2.0E0*I_ESP_I5xz_Hx2y2z_a-4*I_ESP_G3xz_Hx2y2z;
    abcd[iGrid*1323+255] = 2.0E0*I_ESP_I4x2y_Hx2y2z_a-3*I_ESP_G2x2y_Hx2y2z;
    abcd[iGrid*1323+256] = 2.0E0*I_ESP_I4xyz_Hx2y2z_a-3*I_ESP_G2xyz_Hx2y2z;
    abcd[iGrid*1323+257] = 2.0E0*I_ESP_I4x2z_Hx2y2z_a-3*I_ESP_G2x2z_Hx2y2z;
    abcd[iGrid*1323+258] = 2.0E0*I_ESP_I3x3y_Hx2y2z_a-2*I_ESP_Gx3y_Hx2y2z;
    abcd[iGrid*1323+259] = 2.0E0*I_ESP_I3x2yz_Hx2y2z_a-2*I_ESP_Gx2yz_Hx2y2z;
    abcd[iGrid*1323+260] = 2.0E0*I_ESP_I3xy2z_Hx2y2z_a-2*I_ESP_Gxy2z_Hx2y2z;
    abcd[iGrid*1323+261] = 2.0E0*I_ESP_I3x3z_Hx2y2z_a-2*I_ESP_Gx3z_Hx2y2z;
    abcd[iGrid*1323+262] = 2.0E0*I_ESP_I2x4y_Hx2y2z_a-1*I_ESP_G4y_Hx2y2z;
    abcd[iGrid*1323+263] = 2.0E0*I_ESP_I2x3yz_Hx2y2z_a-1*I_ESP_G3yz_Hx2y2z;
    abcd[iGrid*1323+264] = 2.0E0*I_ESP_I2x2y2z_Hx2y2z_a-1*I_ESP_G2y2z_Hx2y2z;
    abcd[iGrid*1323+265] = 2.0E0*I_ESP_I2xy3z_Hx2y2z_a-1*I_ESP_Gy3z_Hx2y2z;
    abcd[iGrid*1323+266] = 2.0E0*I_ESP_I2x4z_Hx2y2z_a-1*I_ESP_G4z_Hx2y2z;
    abcd[iGrid*1323+267] = 2.0E0*I_ESP_Ix5y_Hx2y2z_a;
    abcd[iGrid*1323+268] = 2.0E0*I_ESP_Ix4yz_Hx2y2z_a;
    abcd[iGrid*1323+269] = 2.0E0*I_ESP_Ix3y2z_Hx2y2z_a;
    abcd[iGrid*1323+270] = 2.0E0*I_ESP_Ix2y3z_Hx2y2z_a;
    abcd[iGrid*1323+271] = 2.0E0*I_ESP_Ixy4z_Hx2y2z_a;
    abcd[iGrid*1323+272] = 2.0E0*I_ESP_Ix5z_Hx2y2z_a;
    abcd[iGrid*1323+273] = 2.0E0*I_ESP_I6x_Hxy3z_a-5*I_ESP_G4x_Hxy3z;
    abcd[iGrid*1323+274] = 2.0E0*I_ESP_I5xy_Hxy3z_a-4*I_ESP_G3xy_Hxy3z;
    abcd[iGrid*1323+275] = 2.0E0*I_ESP_I5xz_Hxy3z_a-4*I_ESP_G3xz_Hxy3z;
    abcd[iGrid*1323+276] = 2.0E0*I_ESP_I4x2y_Hxy3z_a-3*I_ESP_G2x2y_Hxy3z;
    abcd[iGrid*1323+277] = 2.0E0*I_ESP_I4xyz_Hxy3z_a-3*I_ESP_G2xyz_Hxy3z;
    abcd[iGrid*1323+278] = 2.0E0*I_ESP_I4x2z_Hxy3z_a-3*I_ESP_G2x2z_Hxy3z;
    abcd[iGrid*1323+279] = 2.0E0*I_ESP_I3x3y_Hxy3z_a-2*I_ESP_Gx3y_Hxy3z;
    abcd[iGrid*1323+280] = 2.0E0*I_ESP_I3x2yz_Hxy3z_a-2*I_ESP_Gx2yz_Hxy3z;
    abcd[iGrid*1323+281] = 2.0E0*I_ESP_I3xy2z_Hxy3z_a-2*I_ESP_Gxy2z_Hxy3z;
    abcd[iGrid*1323+282] = 2.0E0*I_ESP_I3x3z_Hxy3z_a-2*I_ESP_Gx3z_Hxy3z;
    abcd[iGrid*1323+283] = 2.0E0*I_ESP_I2x4y_Hxy3z_a-1*I_ESP_G4y_Hxy3z;
    abcd[iGrid*1323+284] = 2.0E0*I_ESP_I2x3yz_Hxy3z_a-1*I_ESP_G3yz_Hxy3z;
    abcd[iGrid*1323+285] = 2.0E0*I_ESP_I2x2y2z_Hxy3z_a-1*I_ESP_G2y2z_Hxy3z;
    abcd[iGrid*1323+286] = 2.0E0*I_ESP_I2xy3z_Hxy3z_a-1*I_ESP_Gy3z_Hxy3z;
    abcd[iGrid*1323+287] = 2.0E0*I_ESP_I2x4z_Hxy3z_a-1*I_ESP_G4z_Hxy3z;
    abcd[iGrid*1323+288] = 2.0E0*I_ESP_Ix5y_Hxy3z_a;
    abcd[iGrid*1323+289] = 2.0E0*I_ESP_Ix4yz_Hxy3z_a;
    abcd[iGrid*1323+290] = 2.0E0*I_ESP_Ix3y2z_Hxy3z_a;
    abcd[iGrid*1323+291] = 2.0E0*I_ESP_Ix2y3z_Hxy3z_a;
    abcd[iGrid*1323+292] = 2.0E0*I_ESP_Ixy4z_Hxy3z_a;
    abcd[iGrid*1323+293] = 2.0E0*I_ESP_Ix5z_Hxy3z_a;
    abcd[iGrid*1323+294] = 2.0E0*I_ESP_I6x_Hx4z_a-5*I_ESP_G4x_Hx4z;
    abcd[iGrid*1323+295] = 2.0E0*I_ESP_I5xy_Hx4z_a-4*I_ESP_G3xy_Hx4z;
    abcd[iGrid*1323+296] = 2.0E0*I_ESP_I5xz_Hx4z_a-4*I_ESP_G3xz_Hx4z;
    abcd[iGrid*1323+297] = 2.0E0*I_ESP_I4x2y_Hx4z_a-3*I_ESP_G2x2y_Hx4z;
    abcd[iGrid*1323+298] = 2.0E0*I_ESP_I4xyz_Hx4z_a-3*I_ESP_G2xyz_Hx4z;
    abcd[iGrid*1323+299] = 2.0E0*I_ESP_I4x2z_Hx4z_a-3*I_ESP_G2x2z_Hx4z;
    abcd[iGrid*1323+300] = 2.0E0*I_ESP_I3x3y_Hx4z_a-2*I_ESP_Gx3y_Hx4z;
    abcd[iGrid*1323+301] = 2.0E0*I_ESP_I3x2yz_Hx4z_a-2*I_ESP_Gx2yz_Hx4z;
    abcd[iGrid*1323+302] = 2.0E0*I_ESP_I3xy2z_Hx4z_a-2*I_ESP_Gxy2z_Hx4z;
    abcd[iGrid*1323+303] = 2.0E0*I_ESP_I3x3z_Hx4z_a-2*I_ESP_Gx3z_Hx4z;
    abcd[iGrid*1323+304] = 2.0E0*I_ESP_I2x4y_Hx4z_a-1*I_ESP_G4y_Hx4z;
    abcd[iGrid*1323+305] = 2.0E0*I_ESP_I2x3yz_Hx4z_a-1*I_ESP_G3yz_Hx4z;
    abcd[iGrid*1323+306] = 2.0E0*I_ESP_I2x2y2z_Hx4z_a-1*I_ESP_G2y2z_Hx4z;
    abcd[iGrid*1323+307] = 2.0E0*I_ESP_I2xy3z_Hx4z_a-1*I_ESP_Gy3z_Hx4z;
    abcd[iGrid*1323+308] = 2.0E0*I_ESP_I2x4z_Hx4z_a-1*I_ESP_G4z_Hx4z;
    abcd[iGrid*1323+309] = 2.0E0*I_ESP_Ix5y_Hx4z_a;
    abcd[iGrid*1323+310] = 2.0E0*I_ESP_Ix4yz_Hx4z_a;
    abcd[iGrid*1323+311] = 2.0E0*I_ESP_Ix3y2z_Hx4z_a;
    abcd[iGrid*1323+312] = 2.0E0*I_ESP_Ix2y3z_Hx4z_a;
    abcd[iGrid*1323+313] = 2.0E0*I_ESP_Ixy4z_Hx4z_a;
    abcd[iGrid*1323+314] = 2.0E0*I_ESP_Ix5z_Hx4z_a;
    abcd[iGrid*1323+315] = 2.0E0*I_ESP_I6x_H5y_a-5*I_ESP_G4x_H5y;
    abcd[iGrid*1323+316] = 2.0E0*I_ESP_I5xy_H5y_a-4*I_ESP_G3xy_H5y;
    abcd[iGrid*1323+317] = 2.0E0*I_ESP_I5xz_H5y_a-4*I_ESP_G3xz_H5y;
    abcd[iGrid*1323+318] = 2.0E0*I_ESP_I4x2y_H5y_a-3*I_ESP_G2x2y_H5y;
    abcd[iGrid*1323+319] = 2.0E0*I_ESP_I4xyz_H5y_a-3*I_ESP_G2xyz_H5y;
    abcd[iGrid*1323+320] = 2.0E0*I_ESP_I4x2z_H5y_a-3*I_ESP_G2x2z_H5y;
    abcd[iGrid*1323+321] = 2.0E0*I_ESP_I3x3y_H5y_a-2*I_ESP_Gx3y_H5y;
    abcd[iGrid*1323+322] = 2.0E0*I_ESP_I3x2yz_H5y_a-2*I_ESP_Gx2yz_H5y;
    abcd[iGrid*1323+323] = 2.0E0*I_ESP_I3xy2z_H5y_a-2*I_ESP_Gxy2z_H5y;
    abcd[iGrid*1323+324] = 2.0E0*I_ESP_I3x3z_H5y_a-2*I_ESP_Gx3z_H5y;
    abcd[iGrid*1323+325] = 2.0E0*I_ESP_I2x4y_H5y_a-1*I_ESP_G4y_H5y;
    abcd[iGrid*1323+326] = 2.0E0*I_ESP_I2x3yz_H5y_a-1*I_ESP_G3yz_H5y;
    abcd[iGrid*1323+327] = 2.0E0*I_ESP_I2x2y2z_H5y_a-1*I_ESP_G2y2z_H5y;
    abcd[iGrid*1323+328] = 2.0E0*I_ESP_I2xy3z_H5y_a-1*I_ESP_Gy3z_H5y;
    abcd[iGrid*1323+329] = 2.0E0*I_ESP_I2x4z_H5y_a-1*I_ESP_G4z_H5y;
    abcd[iGrid*1323+330] = 2.0E0*I_ESP_Ix5y_H5y_a;
    abcd[iGrid*1323+331] = 2.0E0*I_ESP_Ix4yz_H5y_a;
    abcd[iGrid*1323+332] = 2.0E0*I_ESP_Ix3y2z_H5y_a;
    abcd[iGrid*1323+333] = 2.0E0*I_ESP_Ix2y3z_H5y_a;
    abcd[iGrid*1323+334] = 2.0E0*I_ESP_Ixy4z_H5y_a;
    abcd[iGrid*1323+335] = 2.0E0*I_ESP_Ix5z_H5y_a;
    abcd[iGrid*1323+336] = 2.0E0*I_ESP_I6x_H4yz_a-5*I_ESP_G4x_H4yz;
    abcd[iGrid*1323+337] = 2.0E0*I_ESP_I5xy_H4yz_a-4*I_ESP_G3xy_H4yz;
    abcd[iGrid*1323+338] = 2.0E0*I_ESP_I5xz_H4yz_a-4*I_ESP_G3xz_H4yz;
    abcd[iGrid*1323+339] = 2.0E0*I_ESP_I4x2y_H4yz_a-3*I_ESP_G2x2y_H4yz;
    abcd[iGrid*1323+340] = 2.0E0*I_ESP_I4xyz_H4yz_a-3*I_ESP_G2xyz_H4yz;
    abcd[iGrid*1323+341] = 2.0E0*I_ESP_I4x2z_H4yz_a-3*I_ESP_G2x2z_H4yz;
    abcd[iGrid*1323+342] = 2.0E0*I_ESP_I3x3y_H4yz_a-2*I_ESP_Gx3y_H4yz;
    abcd[iGrid*1323+343] = 2.0E0*I_ESP_I3x2yz_H4yz_a-2*I_ESP_Gx2yz_H4yz;
    abcd[iGrid*1323+344] = 2.0E0*I_ESP_I3xy2z_H4yz_a-2*I_ESP_Gxy2z_H4yz;
    abcd[iGrid*1323+345] = 2.0E0*I_ESP_I3x3z_H4yz_a-2*I_ESP_Gx3z_H4yz;
    abcd[iGrid*1323+346] = 2.0E0*I_ESP_I2x4y_H4yz_a-1*I_ESP_G4y_H4yz;
    abcd[iGrid*1323+347] = 2.0E0*I_ESP_I2x3yz_H4yz_a-1*I_ESP_G3yz_H4yz;
    abcd[iGrid*1323+348] = 2.0E0*I_ESP_I2x2y2z_H4yz_a-1*I_ESP_G2y2z_H4yz;
    abcd[iGrid*1323+349] = 2.0E0*I_ESP_I2xy3z_H4yz_a-1*I_ESP_Gy3z_H4yz;
    abcd[iGrid*1323+350] = 2.0E0*I_ESP_I2x4z_H4yz_a-1*I_ESP_G4z_H4yz;
    abcd[iGrid*1323+351] = 2.0E0*I_ESP_Ix5y_H4yz_a;
    abcd[iGrid*1323+352] = 2.0E0*I_ESP_Ix4yz_H4yz_a;
    abcd[iGrid*1323+353] = 2.0E0*I_ESP_Ix3y2z_H4yz_a;
    abcd[iGrid*1323+354] = 2.0E0*I_ESP_Ix2y3z_H4yz_a;
    abcd[iGrid*1323+355] = 2.0E0*I_ESP_Ixy4z_H4yz_a;
    abcd[iGrid*1323+356] = 2.0E0*I_ESP_Ix5z_H4yz_a;
    abcd[iGrid*1323+357] = 2.0E0*I_ESP_I6x_H3y2z_a-5*I_ESP_G4x_H3y2z;
    abcd[iGrid*1323+358] = 2.0E0*I_ESP_I5xy_H3y2z_a-4*I_ESP_G3xy_H3y2z;
    abcd[iGrid*1323+359] = 2.0E0*I_ESP_I5xz_H3y2z_a-4*I_ESP_G3xz_H3y2z;
    abcd[iGrid*1323+360] = 2.0E0*I_ESP_I4x2y_H3y2z_a-3*I_ESP_G2x2y_H3y2z;
    abcd[iGrid*1323+361] = 2.0E0*I_ESP_I4xyz_H3y2z_a-3*I_ESP_G2xyz_H3y2z;
    abcd[iGrid*1323+362] = 2.0E0*I_ESP_I4x2z_H3y2z_a-3*I_ESP_G2x2z_H3y2z;
    abcd[iGrid*1323+363] = 2.0E0*I_ESP_I3x3y_H3y2z_a-2*I_ESP_Gx3y_H3y2z;
    abcd[iGrid*1323+364] = 2.0E0*I_ESP_I3x2yz_H3y2z_a-2*I_ESP_Gx2yz_H3y2z;
    abcd[iGrid*1323+365] = 2.0E0*I_ESP_I3xy2z_H3y2z_a-2*I_ESP_Gxy2z_H3y2z;
    abcd[iGrid*1323+366] = 2.0E0*I_ESP_I3x3z_H3y2z_a-2*I_ESP_Gx3z_H3y2z;
    abcd[iGrid*1323+367] = 2.0E0*I_ESP_I2x4y_H3y2z_a-1*I_ESP_G4y_H3y2z;
    abcd[iGrid*1323+368] = 2.0E0*I_ESP_I2x3yz_H3y2z_a-1*I_ESP_G3yz_H3y2z;
    abcd[iGrid*1323+369] = 2.0E0*I_ESP_I2x2y2z_H3y2z_a-1*I_ESP_G2y2z_H3y2z;
    abcd[iGrid*1323+370] = 2.0E0*I_ESP_I2xy3z_H3y2z_a-1*I_ESP_Gy3z_H3y2z;
    abcd[iGrid*1323+371] = 2.0E0*I_ESP_I2x4z_H3y2z_a-1*I_ESP_G4z_H3y2z;
    abcd[iGrid*1323+372] = 2.0E0*I_ESP_Ix5y_H3y2z_a;
    abcd[iGrid*1323+373] = 2.0E0*I_ESP_Ix4yz_H3y2z_a;
    abcd[iGrid*1323+374] = 2.0E0*I_ESP_Ix3y2z_H3y2z_a;
    abcd[iGrid*1323+375] = 2.0E0*I_ESP_Ix2y3z_H3y2z_a;
    abcd[iGrid*1323+376] = 2.0E0*I_ESP_Ixy4z_H3y2z_a;
    abcd[iGrid*1323+377] = 2.0E0*I_ESP_Ix5z_H3y2z_a;
    abcd[iGrid*1323+378] = 2.0E0*I_ESP_I6x_H2y3z_a-5*I_ESP_G4x_H2y3z;
    abcd[iGrid*1323+379] = 2.0E0*I_ESP_I5xy_H2y3z_a-4*I_ESP_G3xy_H2y3z;
    abcd[iGrid*1323+380] = 2.0E0*I_ESP_I5xz_H2y3z_a-4*I_ESP_G3xz_H2y3z;
    abcd[iGrid*1323+381] = 2.0E0*I_ESP_I4x2y_H2y3z_a-3*I_ESP_G2x2y_H2y3z;
    abcd[iGrid*1323+382] = 2.0E0*I_ESP_I4xyz_H2y3z_a-3*I_ESP_G2xyz_H2y3z;
    abcd[iGrid*1323+383] = 2.0E0*I_ESP_I4x2z_H2y3z_a-3*I_ESP_G2x2z_H2y3z;
    abcd[iGrid*1323+384] = 2.0E0*I_ESP_I3x3y_H2y3z_a-2*I_ESP_Gx3y_H2y3z;
    abcd[iGrid*1323+385] = 2.0E0*I_ESP_I3x2yz_H2y3z_a-2*I_ESP_Gx2yz_H2y3z;
    abcd[iGrid*1323+386] = 2.0E0*I_ESP_I3xy2z_H2y3z_a-2*I_ESP_Gxy2z_H2y3z;
    abcd[iGrid*1323+387] = 2.0E0*I_ESP_I3x3z_H2y3z_a-2*I_ESP_Gx3z_H2y3z;
    abcd[iGrid*1323+388] = 2.0E0*I_ESP_I2x4y_H2y3z_a-1*I_ESP_G4y_H2y3z;
    abcd[iGrid*1323+389] = 2.0E0*I_ESP_I2x3yz_H2y3z_a-1*I_ESP_G3yz_H2y3z;
    abcd[iGrid*1323+390] = 2.0E0*I_ESP_I2x2y2z_H2y3z_a-1*I_ESP_G2y2z_H2y3z;
    abcd[iGrid*1323+391] = 2.0E0*I_ESP_I2xy3z_H2y3z_a-1*I_ESP_Gy3z_H2y3z;
    abcd[iGrid*1323+392] = 2.0E0*I_ESP_I2x4z_H2y3z_a-1*I_ESP_G4z_H2y3z;
    abcd[iGrid*1323+393] = 2.0E0*I_ESP_Ix5y_H2y3z_a;
    abcd[iGrid*1323+394] = 2.0E0*I_ESP_Ix4yz_H2y3z_a;
    abcd[iGrid*1323+395] = 2.0E0*I_ESP_Ix3y2z_H2y3z_a;
    abcd[iGrid*1323+396] = 2.0E0*I_ESP_Ix2y3z_H2y3z_a;
    abcd[iGrid*1323+397] = 2.0E0*I_ESP_Ixy4z_H2y3z_a;
    abcd[iGrid*1323+398] = 2.0E0*I_ESP_Ix5z_H2y3z_a;
    abcd[iGrid*1323+399] = 2.0E0*I_ESP_I6x_Hy4z_a-5*I_ESP_G4x_Hy4z;
    abcd[iGrid*1323+400] = 2.0E0*I_ESP_I5xy_Hy4z_a-4*I_ESP_G3xy_Hy4z;
    abcd[iGrid*1323+401] = 2.0E0*I_ESP_I5xz_Hy4z_a-4*I_ESP_G3xz_Hy4z;
    abcd[iGrid*1323+402] = 2.0E0*I_ESP_I4x2y_Hy4z_a-3*I_ESP_G2x2y_Hy4z;
    abcd[iGrid*1323+403] = 2.0E0*I_ESP_I4xyz_Hy4z_a-3*I_ESP_G2xyz_Hy4z;
    abcd[iGrid*1323+404] = 2.0E0*I_ESP_I4x2z_Hy4z_a-3*I_ESP_G2x2z_Hy4z;
    abcd[iGrid*1323+405] = 2.0E0*I_ESP_I3x3y_Hy4z_a-2*I_ESP_Gx3y_Hy4z;
    abcd[iGrid*1323+406] = 2.0E0*I_ESP_I3x2yz_Hy4z_a-2*I_ESP_Gx2yz_Hy4z;
    abcd[iGrid*1323+407] = 2.0E0*I_ESP_I3xy2z_Hy4z_a-2*I_ESP_Gxy2z_Hy4z;
    abcd[iGrid*1323+408] = 2.0E0*I_ESP_I3x3z_Hy4z_a-2*I_ESP_Gx3z_Hy4z;
    abcd[iGrid*1323+409] = 2.0E0*I_ESP_I2x4y_Hy4z_a-1*I_ESP_G4y_Hy4z;
    abcd[iGrid*1323+410] = 2.0E0*I_ESP_I2x3yz_Hy4z_a-1*I_ESP_G3yz_Hy4z;
    abcd[iGrid*1323+411] = 2.0E0*I_ESP_I2x2y2z_Hy4z_a-1*I_ESP_G2y2z_Hy4z;
    abcd[iGrid*1323+412] = 2.0E0*I_ESP_I2xy3z_Hy4z_a-1*I_ESP_Gy3z_Hy4z;
    abcd[iGrid*1323+413] = 2.0E0*I_ESP_I2x4z_Hy4z_a-1*I_ESP_G4z_Hy4z;
    abcd[iGrid*1323+414] = 2.0E0*I_ESP_Ix5y_Hy4z_a;
    abcd[iGrid*1323+415] = 2.0E0*I_ESP_Ix4yz_Hy4z_a;
    abcd[iGrid*1323+416] = 2.0E0*I_ESP_Ix3y2z_Hy4z_a;
    abcd[iGrid*1323+417] = 2.0E0*I_ESP_Ix2y3z_Hy4z_a;
    abcd[iGrid*1323+418] = 2.0E0*I_ESP_Ixy4z_Hy4z_a;
    abcd[iGrid*1323+419] = 2.0E0*I_ESP_Ix5z_Hy4z_a;
    abcd[iGrid*1323+420] = 2.0E0*I_ESP_I6x_H5z_a-5*I_ESP_G4x_H5z;
    abcd[iGrid*1323+421] = 2.0E0*I_ESP_I5xy_H5z_a-4*I_ESP_G3xy_H5z;
    abcd[iGrid*1323+422] = 2.0E0*I_ESP_I5xz_H5z_a-4*I_ESP_G3xz_H5z;
    abcd[iGrid*1323+423] = 2.0E0*I_ESP_I4x2y_H5z_a-3*I_ESP_G2x2y_H5z;
    abcd[iGrid*1323+424] = 2.0E0*I_ESP_I4xyz_H5z_a-3*I_ESP_G2xyz_H5z;
    abcd[iGrid*1323+425] = 2.0E0*I_ESP_I4x2z_H5z_a-3*I_ESP_G2x2z_H5z;
    abcd[iGrid*1323+426] = 2.0E0*I_ESP_I3x3y_H5z_a-2*I_ESP_Gx3y_H5z;
    abcd[iGrid*1323+427] = 2.0E0*I_ESP_I3x2yz_H5z_a-2*I_ESP_Gx2yz_H5z;
    abcd[iGrid*1323+428] = 2.0E0*I_ESP_I3xy2z_H5z_a-2*I_ESP_Gxy2z_H5z;
    abcd[iGrid*1323+429] = 2.0E0*I_ESP_I3x3z_H5z_a-2*I_ESP_Gx3z_H5z;
    abcd[iGrid*1323+430] = 2.0E0*I_ESP_I2x4y_H5z_a-1*I_ESP_G4y_H5z;
    abcd[iGrid*1323+431] = 2.0E0*I_ESP_I2x3yz_H5z_a-1*I_ESP_G3yz_H5z;
    abcd[iGrid*1323+432] = 2.0E0*I_ESP_I2x2y2z_H5z_a-1*I_ESP_G2y2z_H5z;
    abcd[iGrid*1323+433] = 2.0E0*I_ESP_I2xy3z_H5z_a-1*I_ESP_Gy3z_H5z;
    abcd[iGrid*1323+434] = 2.0E0*I_ESP_I2x4z_H5z_a-1*I_ESP_G4z_H5z;
    abcd[iGrid*1323+435] = 2.0E0*I_ESP_Ix5y_H5z_a;
    abcd[iGrid*1323+436] = 2.0E0*I_ESP_Ix4yz_H5z_a;
    abcd[iGrid*1323+437] = 2.0E0*I_ESP_Ix3y2z_H5z_a;
    abcd[iGrid*1323+438] = 2.0E0*I_ESP_Ix2y3z_H5z_a;
    abcd[iGrid*1323+439] = 2.0E0*I_ESP_Ixy4z_H5z_a;
    abcd[iGrid*1323+440] = 2.0E0*I_ESP_Ix5z_H5z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_H_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_H_a
     * RHS shell quartet name: SQ_ESP_G_H
     ************************************************************/
    abcd[iGrid*1323+441] = 2.0E0*I_ESP_I5xy_H5x_a;
    abcd[iGrid*1323+442] = 2.0E0*I_ESP_I4x2y_H5x_a-1*I_ESP_G4x_H5x;
    abcd[iGrid*1323+443] = 2.0E0*I_ESP_I4xyz_H5x_a;
    abcd[iGrid*1323+444] = 2.0E0*I_ESP_I3x3y_H5x_a-2*I_ESP_G3xy_H5x;
    abcd[iGrid*1323+445] = 2.0E0*I_ESP_I3x2yz_H5x_a-1*I_ESP_G3xz_H5x;
    abcd[iGrid*1323+446] = 2.0E0*I_ESP_I3xy2z_H5x_a;
    abcd[iGrid*1323+447] = 2.0E0*I_ESP_I2x4y_H5x_a-3*I_ESP_G2x2y_H5x;
    abcd[iGrid*1323+448] = 2.0E0*I_ESP_I2x3yz_H5x_a-2*I_ESP_G2xyz_H5x;
    abcd[iGrid*1323+449] = 2.0E0*I_ESP_I2x2y2z_H5x_a-1*I_ESP_G2x2z_H5x;
    abcd[iGrid*1323+450] = 2.0E0*I_ESP_I2xy3z_H5x_a;
    abcd[iGrid*1323+451] = 2.0E0*I_ESP_Ix5y_H5x_a-4*I_ESP_Gx3y_H5x;
    abcd[iGrid*1323+452] = 2.0E0*I_ESP_Ix4yz_H5x_a-3*I_ESP_Gx2yz_H5x;
    abcd[iGrid*1323+453] = 2.0E0*I_ESP_Ix3y2z_H5x_a-2*I_ESP_Gxy2z_H5x;
    abcd[iGrid*1323+454] = 2.0E0*I_ESP_Ix2y3z_H5x_a-1*I_ESP_Gx3z_H5x;
    abcd[iGrid*1323+455] = 2.0E0*I_ESP_Ixy4z_H5x_a;
    abcd[iGrid*1323+456] = 2.0E0*I_ESP_I6y_H5x_a-5*I_ESP_G4y_H5x;
    abcd[iGrid*1323+457] = 2.0E0*I_ESP_I5yz_H5x_a-4*I_ESP_G3yz_H5x;
    abcd[iGrid*1323+458] = 2.0E0*I_ESP_I4y2z_H5x_a-3*I_ESP_G2y2z_H5x;
    abcd[iGrid*1323+459] = 2.0E0*I_ESP_I3y3z_H5x_a-2*I_ESP_Gy3z_H5x;
    abcd[iGrid*1323+460] = 2.0E0*I_ESP_I2y4z_H5x_a-1*I_ESP_G4z_H5x;
    abcd[iGrid*1323+461] = 2.0E0*I_ESP_Iy5z_H5x_a;
    abcd[iGrid*1323+462] = 2.0E0*I_ESP_I5xy_H4xy_a;
    abcd[iGrid*1323+463] = 2.0E0*I_ESP_I4x2y_H4xy_a-1*I_ESP_G4x_H4xy;
    abcd[iGrid*1323+464] = 2.0E0*I_ESP_I4xyz_H4xy_a;
    abcd[iGrid*1323+465] = 2.0E0*I_ESP_I3x3y_H4xy_a-2*I_ESP_G3xy_H4xy;
    abcd[iGrid*1323+466] = 2.0E0*I_ESP_I3x2yz_H4xy_a-1*I_ESP_G3xz_H4xy;
    abcd[iGrid*1323+467] = 2.0E0*I_ESP_I3xy2z_H4xy_a;
    abcd[iGrid*1323+468] = 2.0E0*I_ESP_I2x4y_H4xy_a-3*I_ESP_G2x2y_H4xy;
    abcd[iGrid*1323+469] = 2.0E0*I_ESP_I2x3yz_H4xy_a-2*I_ESP_G2xyz_H4xy;
    abcd[iGrid*1323+470] = 2.0E0*I_ESP_I2x2y2z_H4xy_a-1*I_ESP_G2x2z_H4xy;
    abcd[iGrid*1323+471] = 2.0E0*I_ESP_I2xy3z_H4xy_a;
    abcd[iGrid*1323+472] = 2.0E0*I_ESP_Ix5y_H4xy_a-4*I_ESP_Gx3y_H4xy;
    abcd[iGrid*1323+473] = 2.0E0*I_ESP_Ix4yz_H4xy_a-3*I_ESP_Gx2yz_H4xy;
    abcd[iGrid*1323+474] = 2.0E0*I_ESP_Ix3y2z_H4xy_a-2*I_ESP_Gxy2z_H4xy;
    abcd[iGrid*1323+475] = 2.0E0*I_ESP_Ix2y3z_H4xy_a-1*I_ESP_Gx3z_H4xy;
    abcd[iGrid*1323+476] = 2.0E0*I_ESP_Ixy4z_H4xy_a;
    abcd[iGrid*1323+477] = 2.0E0*I_ESP_I6y_H4xy_a-5*I_ESP_G4y_H4xy;
    abcd[iGrid*1323+478] = 2.0E0*I_ESP_I5yz_H4xy_a-4*I_ESP_G3yz_H4xy;
    abcd[iGrid*1323+479] = 2.0E0*I_ESP_I4y2z_H4xy_a-3*I_ESP_G2y2z_H4xy;
    abcd[iGrid*1323+480] = 2.0E0*I_ESP_I3y3z_H4xy_a-2*I_ESP_Gy3z_H4xy;
    abcd[iGrid*1323+481] = 2.0E0*I_ESP_I2y4z_H4xy_a-1*I_ESP_G4z_H4xy;
    abcd[iGrid*1323+482] = 2.0E0*I_ESP_Iy5z_H4xy_a;
    abcd[iGrid*1323+483] = 2.0E0*I_ESP_I5xy_H4xz_a;
    abcd[iGrid*1323+484] = 2.0E0*I_ESP_I4x2y_H4xz_a-1*I_ESP_G4x_H4xz;
    abcd[iGrid*1323+485] = 2.0E0*I_ESP_I4xyz_H4xz_a;
    abcd[iGrid*1323+486] = 2.0E0*I_ESP_I3x3y_H4xz_a-2*I_ESP_G3xy_H4xz;
    abcd[iGrid*1323+487] = 2.0E0*I_ESP_I3x2yz_H4xz_a-1*I_ESP_G3xz_H4xz;
    abcd[iGrid*1323+488] = 2.0E0*I_ESP_I3xy2z_H4xz_a;
    abcd[iGrid*1323+489] = 2.0E0*I_ESP_I2x4y_H4xz_a-3*I_ESP_G2x2y_H4xz;
    abcd[iGrid*1323+490] = 2.0E0*I_ESP_I2x3yz_H4xz_a-2*I_ESP_G2xyz_H4xz;
    abcd[iGrid*1323+491] = 2.0E0*I_ESP_I2x2y2z_H4xz_a-1*I_ESP_G2x2z_H4xz;
    abcd[iGrid*1323+492] = 2.0E0*I_ESP_I2xy3z_H4xz_a;
    abcd[iGrid*1323+493] = 2.0E0*I_ESP_Ix5y_H4xz_a-4*I_ESP_Gx3y_H4xz;
    abcd[iGrid*1323+494] = 2.0E0*I_ESP_Ix4yz_H4xz_a-3*I_ESP_Gx2yz_H4xz;
    abcd[iGrid*1323+495] = 2.0E0*I_ESP_Ix3y2z_H4xz_a-2*I_ESP_Gxy2z_H4xz;
    abcd[iGrid*1323+496] = 2.0E0*I_ESP_Ix2y3z_H4xz_a-1*I_ESP_Gx3z_H4xz;
    abcd[iGrid*1323+497] = 2.0E0*I_ESP_Ixy4z_H4xz_a;
    abcd[iGrid*1323+498] = 2.0E0*I_ESP_I6y_H4xz_a-5*I_ESP_G4y_H4xz;
    abcd[iGrid*1323+499] = 2.0E0*I_ESP_I5yz_H4xz_a-4*I_ESP_G3yz_H4xz;
    abcd[iGrid*1323+500] = 2.0E0*I_ESP_I4y2z_H4xz_a-3*I_ESP_G2y2z_H4xz;
    abcd[iGrid*1323+501] = 2.0E0*I_ESP_I3y3z_H4xz_a-2*I_ESP_Gy3z_H4xz;
    abcd[iGrid*1323+502] = 2.0E0*I_ESP_I2y4z_H4xz_a-1*I_ESP_G4z_H4xz;
    abcd[iGrid*1323+503] = 2.0E0*I_ESP_Iy5z_H4xz_a;
    abcd[iGrid*1323+504] = 2.0E0*I_ESP_I5xy_H3x2y_a;
    abcd[iGrid*1323+505] = 2.0E0*I_ESP_I4x2y_H3x2y_a-1*I_ESP_G4x_H3x2y;
    abcd[iGrid*1323+506] = 2.0E0*I_ESP_I4xyz_H3x2y_a;
    abcd[iGrid*1323+507] = 2.0E0*I_ESP_I3x3y_H3x2y_a-2*I_ESP_G3xy_H3x2y;
    abcd[iGrid*1323+508] = 2.0E0*I_ESP_I3x2yz_H3x2y_a-1*I_ESP_G3xz_H3x2y;
    abcd[iGrid*1323+509] = 2.0E0*I_ESP_I3xy2z_H3x2y_a;
    abcd[iGrid*1323+510] = 2.0E0*I_ESP_I2x4y_H3x2y_a-3*I_ESP_G2x2y_H3x2y;
    abcd[iGrid*1323+511] = 2.0E0*I_ESP_I2x3yz_H3x2y_a-2*I_ESP_G2xyz_H3x2y;
    abcd[iGrid*1323+512] = 2.0E0*I_ESP_I2x2y2z_H3x2y_a-1*I_ESP_G2x2z_H3x2y;
    abcd[iGrid*1323+513] = 2.0E0*I_ESP_I2xy3z_H3x2y_a;
    abcd[iGrid*1323+514] = 2.0E0*I_ESP_Ix5y_H3x2y_a-4*I_ESP_Gx3y_H3x2y;
    abcd[iGrid*1323+515] = 2.0E0*I_ESP_Ix4yz_H3x2y_a-3*I_ESP_Gx2yz_H3x2y;
    abcd[iGrid*1323+516] = 2.0E0*I_ESP_Ix3y2z_H3x2y_a-2*I_ESP_Gxy2z_H3x2y;
    abcd[iGrid*1323+517] = 2.0E0*I_ESP_Ix2y3z_H3x2y_a-1*I_ESP_Gx3z_H3x2y;
    abcd[iGrid*1323+518] = 2.0E0*I_ESP_Ixy4z_H3x2y_a;
    abcd[iGrid*1323+519] = 2.0E0*I_ESP_I6y_H3x2y_a-5*I_ESP_G4y_H3x2y;
    abcd[iGrid*1323+520] = 2.0E0*I_ESP_I5yz_H3x2y_a-4*I_ESP_G3yz_H3x2y;
    abcd[iGrid*1323+521] = 2.0E0*I_ESP_I4y2z_H3x2y_a-3*I_ESP_G2y2z_H3x2y;
    abcd[iGrid*1323+522] = 2.0E0*I_ESP_I3y3z_H3x2y_a-2*I_ESP_Gy3z_H3x2y;
    abcd[iGrid*1323+523] = 2.0E0*I_ESP_I2y4z_H3x2y_a-1*I_ESP_G4z_H3x2y;
    abcd[iGrid*1323+524] = 2.0E0*I_ESP_Iy5z_H3x2y_a;
    abcd[iGrid*1323+525] = 2.0E0*I_ESP_I5xy_H3xyz_a;
    abcd[iGrid*1323+526] = 2.0E0*I_ESP_I4x2y_H3xyz_a-1*I_ESP_G4x_H3xyz;
    abcd[iGrid*1323+527] = 2.0E0*I_ESP_I4xyz_H3xyz_a;
    abcd[iGrid*1323+528] = 2.0E0*I_ESP_I3x3y_H3xyz_a-2*I_ESP_G3xy_H3xyz;
    abcd[iGrid*1323+529] = 2.0E0*I_ESP_I3x2yz_H3xyz_a-1*I_ESP_G3xz_H3xyz;
    abcd[iGrid*1323+530] = 2.0E0*I_ESP_I3xy2z_H3xyz_a;
    abcd[iGrid*1323+531] = 2.0E0*I_ESP_I2x4y_H3xyz_a-3*I_ESP_G2x2y_H3xyz;
    abcd[iGrid*1323+532] = 2.0E0*I_ESP_I2x3yz_H3xyz_a-2*I_ESP_G2xyz_H3xyz;
    abcd[iGrid*1323+533] = 2.0E0*I_ESP_I2x2y2z_H3xyz_a-1*I_ESP_G2x2z_H3xyz;
    abcd[iGrid*1323+534] = 2.0E0*I_ESP_I2xy3z_H3xyz_a;
    abcd[iGrid*1323+535] = 2.0E0*I_ESP_Ix5y_H3xyz_a-4*I_ESP_Gx3y_H3xyz;
    abcd[iGrid*1323+536] = 2.0E0*I_ESP_Ix4yz_H3xyz_a-3*I_ESP_Gx2yz_H3xyz;
    abcd[iGrid*1323+537] = 2.0E0*I_ESP_Ix3y2z_H3xyz_a-2*I_ESP_Gxy2z_H3xyz;
    abcd[iGrid*1323+538] = 2.0E0*I_ESP_Ix2y3z_H3xyz_a-1*I_ESP_Gx3z_H3xyz;
    abcd[iGrid*1323+539] = 2.0E0*I_ESP_Ixy4z_H3xyz_a;
    abcd[iGrid*1323+540] = 2.0E0*I_ESP_I6y_H3xyz_a-5*I_ESP_G4y_H3xyz;
    abcd[iGrid*1323+541] = 2.0E0*I_ESP_I5yz_H3xyz_a-4*I_ESP_G3yz_H3xyz;
    abcd[iGrid*1323+542] = 2.0E0*I_ESP_I4y2z_H3xyz_a-3*I_ESP_G2y2z_H3xyz;
    abcd[iGrid*1323+543] = 2.0E0*I_ESP_I3y3z_H3xyz_a-2*I_ESP_Gy3z_H3xyz;
    abcd[iGrid*1323+544] = 2.0E0*I_ESP_I2y4z_H3xyz_a-1*I_ESP_G4z_H3xyz;
    abcd[iGrid*1323+545] = 2.0E0*I_ESP_Iy5z_H3xyz_a;
    abcd[iGrid*1323+546] = 2.0E0*I_ESP_I5xy_H3x2z_a;
    abcd[iGrid*1323+547] = 2.0E0*I_ESP_I4x2y_H3x2z_a-1*I_ESP_G4x_H3x2z;
    abcd[iGrid*1323+548] = 2.0E0*I_ESP_I4xyz_H3x2z_a;
    abcd[iGrid*1323+549] = 2.0E0*I_ESP_I3x3y_H3x2z_a-2*I_ESP_G3xy_H3x2z;
    abcd[iGrid*1323+550] = 2.0E0*I_ESP_I3x2yz_H3x2z_a-1*I_ESP_G3xz_H3x2z;
    abcd[iGrid*1323+551] = 2.0E0*I_ESP_I3xy2z_H3x2z_a;
    abcd[iGrid*1323+552] = 2.0E0*I_ESP_I2x4y_H3x2z_a-3*I_ESP_G2x2y_H3x2z;
    abcd[iGrid*1323+553] = 2.0E0*I_ESP_I2x3yz_H3x2z_a-2*I_ESP_G2xyz_H3x2z;
    abcd[iGrid*1323+554] = 2.0E0*I_ESP_I2x2y2z_H3x2z_a-1*I_ESP_G2x2z_H3x2z;
    abcd[iGrid*1323+555] = 2.0E0*I_ESP_I2xy3z_H3x2z_a;
    abcd[iGrid*1323+556] = 2.0E0*I_ESP_Ix5y_H3x2z_a-4*I_ESP_Gx3y_H3x2z;
    abcd[iGrid*1323+557] = 2.0E0*I_ESP_Ix4yz_H3x2z_a-3*I_ESP_Gx2yz_H3x2z;
    abcd[iGrid*1323+558] = 2.0E0*I_ESP_Ix3y2z_H3x2z_a-2*I_ESP_Gxy2z_H3x2z;
    abcd[iGrid*1323+559] = 2.0E0*I_ESP_Ix2y3z_H3x2z_a-1*I_ESP_Gx3z_H3x2z;
    abcd[iGrid*1323+560] = 2.0E0*I_ESP_Ixy4z_H3x2z_a;
    abcd[iGrid*1323+561] = 2.0E0*I_ESP_I6y_H3x2z_a-5*I_ESP_G4y_H3x2z;
    abcd[iGrid*1323+562] = 2.0E0*I_ESP_I5yz_H3x2z_a-4*I_ESP_G3yz_H3x2z;
    abcd[iGrid*1323+563] = 2.0E0*I_ESP_I4y2z_H3x2z_a-3*I_ESP_G2y2z_H3x2z;
    abcd[iGrid*1323+564] = 2.0E0*I_ESP_I3y3z_H3x2z_a-2*I_ESP_Gy3z_H3x2z;
    abcd[iGrid*1323+565] = 2.0E0*I_ESP_I2y4z_H3x2z_a-1*I_ESP_G4z_H3x2z;
    abcd[iGrid*1323+566] = 2.0E0*I_ESP_Iy5z_H3x2z_a;
    abcd[iGrid*1323+567] = 2.0E0*I_ESP_I5xy_H2x3y_a;
    abcd[iGrid*1323+568] = 2.0E0*I_ESP_I4x2y_H2x3y_a-1*I_ESP_G4x_H2x3y;
    abcd[iGrid*1323+569] = 2.0E0*I_ESP_I4xyz_H2x3y_a;
    abcd[iGrid*1323+570] = 2.0E0*I_ESP_I3x3y_H2x3y_a-2*I_ESP_G3xy_H2x3y;
    abcd[iGrid*1323+571] = 2.0E0*I_ESP_I3x2yz_H2x3y_a-1*I_ESP_G3xz_H2x3y;
    abcd[iGrid*1323+572] = 2.0E0*I_ESP_I3xy2z_H2x3y_a;
    abcd[iGrid*1323+573] = 2.0E0*I_ESP_I2x4y_H2x3y_a-3*I_ESP_G2x2y_H2x3y;
    abcd[iGrid*1323+574] = 2.0E0*I_ESP_I2x3yz_H2x3y_a-2*I_ESP_G2xyz_H2x3y;
    abcd[iGrid*1323+575] = 2.0E0*I_ESP_I2x2y2z_H2x3y_a-1*I_ESP_G2x2z_H2x3y;
    abcd[iGrid*1323+576] = 2.0E0*I_ESP_I2xy3z_H2x3y_a;
    abcd[iGrid*1323+577] = 2.0E0*I_ESP_Ix5y_H2x3y_a-4*I_ESP_Gx3y_H2x3y;
    abcd[iGrid*1323+578] = 2.0E0*I_ESP_Ix4yz_H2x3y_a-3*I_ESP_Gx2yz_H2x3y;
    abcd[iGrid*1323+579] = 2.0E0*I_ESP_Ix3y2z_H2x3y_a-2*I_ESP_Gxy2z_H2x3y;
    abcd[iGrid*1323+580] = 2.0E0*I_ESP_Ix2y3z_H2x3y_a-1*I_ESP_Gx3z_H2x3y;
    abcd[iGrid*1323+581] = 2.0E0*I_ESP_Ixy4z_H2x3y_a;
    abcd[iGrid*1323+582] = 2.0E0*I_ESP_I6y_H2x3y_a-5*I_ESP_G4y_H2x3y;
    abcd[iGrid*1323+583] = 2.0E0*I_ESP_I5yz_H2x3y_a-4*I_ESP_G3yz_H2x3y;
    abcd[iGrid*1323+584] = 2.0E0*I_ESP_I4y2z_H2x3y_a-3*I_ESP_G2y2z_H2x3y;
    abcd[iGrid*1323+585] = 2.0E0*I_ESP_I3y3z_H2x3y_a-2*I_ESP_Gy3z_H2x3y;
    abcd[iGrid*1323+586] = 2.0E0*I_ESP_I2y4z_H2x3y_a-1*I_ESP_G4z_H2x3y;
    abcd[iGrid*1323+587] = 2.0E0*I_ESP_Iy5z_H2x3y_a;
    abcd[iGrid*1323+588] = 2.0E0*I_ESP_I5xy_H2x2yz_a;
    abcd[iGrid*1323+589] = 2.0E0*I_ESP_I4x2y_H2x2yz_a-1*I_ESP_G4x_H2x2yz;
    abcd[iGrid*1323+590] = 2.0E0*I_ESP_I4xyz_H2x2yz_a;
    abcd[iGrid*1323+591] = 2.0E0*I_ESP_I3x3y_H2x2yz_a-2*I_ESP_G3xy_H2x2yz;
    abcd[iGrid*1323+592] = 2.0E0*I_ESP_I3x2yz_H2x2yz_a-1*I_ESP_G3xz_H2x2yz;
    abcd[iGrid*1323+593] = 2.0E0*I_ESP_I3xy2z_H2x2yz_a;
    abcd[iGrid*1323+594] = 2.0E0*I_ESP_I2x4y_H2x2yz_a-3*I_ESP_G2x2y_H2x2yz;
    abcd[iGrid*1323+595] = 2.0E0*I_ESP_I2x3yz_H2x2yz_a-2*I_ESP_G2xyz_H2x2yz;
    abcd[iGrid*1323+596] = 2.0E0*I_ESP_I2x2y2z_H2x2yz_a-1*I_ESP_G2x2z_H2x2yz;
    abcd[iGrid*1323+597] = 2.0E0*I_ESP_I2xy3z_H2x2yz_a;
    abcd[iGrid*1323+598] = 2.0E0*I_ESP_Ix5y_H2x2yz_a-4*I_ESP_Gx3y_H2x2yz;
    abcd[iGrid*1323+599] = 2.0E0*I_ESP_Ix4yz_H2x2yz_a-3*I_ESP_Gx2yz_H2x2yz;
    abcd[iGrid*1323+600] = 2.0E0*I_ESP_Ix3y2z_H2x2yz_a-2*I_ESP_Gxy2z_H2x2yz;
    abcd[iGrid*1323+601] = 2.0E0*I_ESP_Ix2y3z_H2x2yz_a-1*I_ESP_Gx3z_H2x2yz;
    abcd[iGrid*1323+602] = 2.0E0*I_ESP_Ixy4z_H2x2yz_a;
    abcd[iGrid*1323+603] = 2.0E0*I_ESP_I6y_H2x2yz_a-5*I_ESP_G4y_H2x2yz;
    abcd[iGrid*1323+604] = 2.0E0*I_ESP_I5yz_H2x2yz_a-4*I_ESP_G3yz_H2x2yz;
    abcd[iGrid*1323+605] = 2.0E0*I_ESP_I4y2z_H2x2yz_a-3*I_ESP_G2y2z_H2x2yz;
    abcd[iGrid*1323+606] = 2.0E0*I_ESP_I3y3z_H2x2yz_a-2*I_ESP_Gy3z_H2x2yz;
    abcd[iGrid*1323+607] = 2.0E0*I_ESP_I2y4z_H2x2yz_a-1*I_ESP_G4z_H2x2yz;
    abcd[iGrid*1323+608] = 2.0E0*I_ESP_Iy5z_H2x2yz_a;
    abcd[iGrid*1323+609] = 2.0E0*I_ESP_I5xy_H2xy2z_a;
    abcd[iGrid*1323+610] = 2.0E0*I_ESP_I4x2y_H2xy2z_a-1*I_ESP_G4x_H2xy2z;
    abcd[iGrid*1323+611] = 2.0E0*I_ESP_I4xyz_H2xy2z_a;
    abcd[iGrid*1323+612] = 2.0E0*I_ESP_I3x3y_H2xy2z_a-2*I_ESP_G3xy_H2xy2z;
    abcd[iGrid*1323+613] = 2.0E0*I_ESP_I3x2yz_H2xy2z_a-1*I_ESP_G3xz_H2xy2z;
    abcd[iGrid*1323+614] = 2.0E0*I_ESP_I3xy2z_H2xy2z_a;
    abcd[iGrid*1323+615] = 2.0E0*I_ESP_I2x4y_H2xy2z_a-3*I_ESP_G2x2y_H2xy2z;
    abcd[iGrid*1323+616] = 2.0E0*I_ESP_I2x3yz_H2xy2z_a-2*I_ESP_G2xyz_H2xy2z;
    abcd[iGrid*1323+617] = 2.0E0*I_ESP_I2x2y2z_H2xy2z_a-1*I_ESP_G2x2z_H2xy2z;
    abcd[iGrid*1323+618] = 2.0E0*I_ESP_I2xy3z_H2xy2z_a;
    abcd[iGrid*1323+619] = 2.0E0*I_ESP_Ix5y_H2xy2z_a-4*I_ESP_Gx3y_H2xy2z;
    abcd[iGrid*1323+620] = 2.0E0*I_ESP_Ix4yz_H2xy2z_a-3*I_ESP_Gx2yz_H2xy2z;
    abcd[iGrid*1323+621] = 2.0E0*I_ESP_Ix3y2z_H2xy2z_a-2*I_ESP_Gxy2z_H2xy2z;
    abcd[iGrid*1323+622] = 2.0E0*I_ESP_Ix2y3z_H2xy2z_a-1*I_ESP_Gx3z_H2xy2z;
    abcd[iGrid*1323+623] = 2.0E0*I_ESP_Ixy4z_H2xy2z_a;
    abcd[iGrid*1323+624] = 2.0E0*I_ESP_I6y_H2xy2z_a-5*I_ESP_G4y_H2xy2z;
    abcd[iGrid*1323+625] = 2.0E0*I_ESP_I5yz_H2xy2z_a-4*I_ESP_G3yz_H2xy2z;
    abcd[iGrid*1323+626] = 2.0E0*I_ESP_I4y2z_H2xy2z_a-3*I_ESP_G2y2z_H2xy2z;
    abcd[iGrid*1323+627] = 2.0E0*I_ESP_I3y3z_H2xy2z_a-2*I_ESP_Gy3z_H2xy2z;
    abcd[iGrid*1323+628] = 2.0E0*I_ESP_I2y4z_H2xy2z_a-1*I_ESP_G4z_H2xy2z;
    abcd[iGrid*1323+629] = 2.0E0*I_ESP_Iy5z_H2xy2z_a;
    abcd[iGrid*1323+630] = 2.0E0*I_ESP_I5xy_H2x3z_a;
    abcd[iGrid*1323+631] = 2.0E0*I_ESP_I4x2y_H2x3z_a-1*I_ESP_G4x_H2x3z;
    abcd[iGrid*1323+632] = 2.0E0*I_ESP_I4xyz_H2x3z_a;
    abcd[iGrid*1323+633] = 2.0E0*I_ESP_I3x3y_H2x3z_a-2*I_ESP_G3xy_H2x3z;
    abcd[iGrid*1323+634] = 2.0E0*I_ESP_I3x2yz_H2x3z_a-1*I_ESP_G3xz_H2x3z;
    abcd[iGrid*1323+635] = 2.0E0*I_ESP_I3xy2z_H2x3z_a;
    abcd[iGrid*1323+636] = 2.0E0*I_ESP_I2x4y_H2x3z_a-3*I_ESP_G2x2y_H2x3z;
    abcd[iGrid*1323+637] = 2.0E0*I_ESP_I2x3yz_H2x3z_a-2*I_ESP_G2xyz_H2x3z;
    abcd[iGrid*1323+638] = 2.0E0*I_ESP_I2x2y2z_H2x3z_a-1*I_ESP_G2x2z_H2x3z;
    abcd[iGrid*1323+639] = 2.0E0*I_ESP_I2xy3z_H2x3z_a;
    abcd[iGrid*1323+640] = 2.0E0*I_ESP_Ix5y_H2x3z_a-4*I_ESP_Gx3y_H2x3z;
    abcd[iGrid*1323+641] = 2.0E0*I_ESP_Ix4yz_H2x3z_a-3*I_ESP_Gx2yz_H2x3z;
    abcd[iGrid*1323+642] = 2.0E0*I_ESP_Ix3y2z_H2x3z_a-2*I_ESP_Gxy2z_H2x3z;
    abcd[iGrid*1323+643] = 2.0E0*I_ESP_Ix2y3z_H2x3z_a-1*I_ESP_Gx3z_H2x3z;
    abcd[iGrid*1323+644] = 2.0E0*I_ESP_Ixy4z_H2x3z_a;
    abcd[iGrid*1323+645] = 2.0E0*I_ESP_I6y_H2x3z_a-5*I_ESP_G4y_H2x3z;
    abcd[iGrid*1323+646] = 2.0E0*I_ESP_I5yz_H2x3z_a-4*I_ESP_G3yz_H2x3z;
    abcd[iGrid*1323+647] = 2.0E0*I_ESP_I4y2z_H2x3z_a-3*I_ESP_G2y2z_H2x3z;
    abcd[iGrid*1323+648] = 2.0E0*I_ESP_I3y3z_H2x3z_a-2*I_ESP_Gy3z_H2x3z;
    abcd[iGrid*1323+649] = 2.0E0*I_ESP_I2y4z_H2x3z_a-1*I_ESP_G4z_H2x3z;
    abcd[iGrid*1323+650] = 2.0E0*I_ESP_Iy5z_H2x3z_a;
    abcd[iGrid*1323+651] = 2.0E0*I_ESP_I5xy_Hx4y_a;
    abcd[iGrid*1323+652] = 2.0E0*I_ESP_I4x2y_Hx4y_a-1*I_ESP_G4x_Hx4y;
    abcd[iGrid*1323+653] = 2.0E0*I_ESP_I4xyz_Hx4y_a;
    abcd[iGrid*1323+654] = 2.0E0*I_ESP_I3x3y_Hx4y_a-2*I_ESP_G3xy_Hx4y;
    abcd[iGrid*1323+655] = 2.0E0*I_ESP_I3x2yz_Hx4y_a-1*I_ESP_G3xz_Hx4y;
    abcd[iGrid*1323+656] = 2.0E0*I_ESP_I3xy2z_Hx4y_a;
    abcd[iGrid*1323+657] = 2.0E0*I_ESP_I2x4y_Hx4y_a-3*I_ESP_G2x2y_Hx4y;
    abcd[iGrid*1323+658] = 2.0E0*I_ESP_I2x3yz_Hx4y_a-2*I_ESP_G2xyz_Hx4y;
    abcd[iGrid*1323+659] = 2.0E0*I_ESP_I2x2y2z_Hx4y_a-1*I_ESP_G2x2z_Hx4y;
    abcd[iGrid*1323+660] = 2.0E0*I_ESP_I2xy3z_Hx4y_a;
    abcd[iGrid*1323+661] = 2.0E0*I_ESP_Ix5y_Hx4y_a-4*I_ESP_Gx3y_Hx4y;
    abcd[iGrid*1323+662] = 2.0E0*I_ESP_Ix4yz_Hx4y_a-3*I_ESP_Gx2yz_Hx4y;
    abcd[iGrid*1323+663] = 2.0E0*I_ESP_Ix3y2z_Hx4y_a-2*I_ESP_Gxy2z_Hx4y;
    abcd[iGrid*1323+664] = 2.0E0*I_ESP_Ix2y3z_Hx4y_a-1*I_ESP_Gx3z_Hx4y;
    abcd[iGrid*1323+665] = 2.0E0*I_ESP_Ixy4z_Hx4y_a;
    abcd[iGrid*1323+666] = 2.0E0*I_ESP_I6y_Hx4y_a-5*I_ESP_G4y_Hx4y;
    abcd[iGrid*1323+667] = 2.0E0*I_ESP_I5yz_Hx4y_a-4*I_ESP_G3yz_Hx4y;
    abcd[iGrid*1323+668] = 2.0E0*I_ESP_I4y2z_Hx4y_a-3*I_ESP_G2y2z_Hx4y;
    abcd[iGrid*1323+669] = 2.0E0*I_ESP_I3y3z_Hx4y_a-2*I_ESP_Gy3z_Hx4y;
    abcd[iGrid*1323+670] = 2.0E0*I_ESP_I2y4z_Hx4y_a-1*I_ESP_G4z_Hx4y;
    abcd[iGrid*1323+671] = 2.0E0*I_ESP_Iy5z_Hx4y_a;
    abcd[iGrid*1323+672] = 2.0E0*I_ESP_I5xy_Hx3yz_a;
    abcd[iGrid*1323+673] = 2.0E0*I_ESP_I4x2y_Hx3yz_a-1*I_ESP_G4x_Hx3yz;
    abcd[iGrid*1323+674] = 2.0E0*I_ESP_I4xyz_Hx3yz_a;
    abcd[iGrid*1323+675] = 2.0E0*I_ESP_I3x3y_Hx3yz_a-2*I_ESP_G3xy_Hx3yz;
    abcd[iGrid*1323+676] = 2.0E0*I_ESP_I3x2yz_Hx3yz_a-1*I_ESP_G3xz_Hx3yz;
    abcd[iGrid*1323+677] = 2.0E0*I_ESP_I3xy2z_Hx3yz_a;
    abcd[iGrid*1323+678] = 2.0E0*I_ESP_I2x4y_Hx3yz_a-3*I_ESP_G2x2y_Hx3yz;
    abcd[iGrid*1323+679] = 2.0E0*I_ESP_I2x3yz_Hx3yz_a-2*I_ESP_G2xyz_Hx3yz;
    abcd[iGrid*1323+680] = 2.0E0*I_ESP_I2x2y2z_Hx3yz_a-1*I_ESP_G2x2z_Hx3yz;
    abcd[iGrid*1323+681] = 2.0E0*I_ESP_I2xy3z_Hx3yz_a;
    abcd[iGrid*1323+682] = 2.0E0*I_ESP_Ix5y_Hx3yz_a-4*I_ESP_Gx3y_Hx3yz;
    abcd[iGrid*1323+683] = 2.0E0*I_ESP_Ix4yz_Hx3yz_a-3*I_ESP_Gx2yz_Hx3yz;
    abcd[iGrid*1323+684] = 2.0E0*I_ESP_Ix3y2z_Hx3yz_a-2*I_ESP_Gxy2z_Hx3yz;
    abcd[iGrid*1323+685] = 2.0E0*I_ESP_Ix2y3z_Hx3yz_a-1*I_ESP_Gx3z_Hx3yz;
    abcd[iGrid*1323+686] = 2.0E0*I_ESP_Ixy4z_Hx3yz_a;
    abcd[iGrid*1323+687] = 2.0E0*I_ESP_I6y_Hx3yz_a-5*I_ESP_G4y_Hx3yz;
    abcd[iGrid*1323+688] = 2.0E0*I_ESP_I5yz_Hx3yz_a-4*I_ESP_G3yz_Hx3yz;
    abcd[iGrid*1323+689] = 2.0E0*I_ESP_I4y2z_Hx3yz_a-3*I_ESP_G2y2z_Hx3yz;
    abcd[iGrid*1323+690] = 2.0E0*I_ESP_I3y3z_Hx3yz_a-2*I_ESP_Gy3z_Hx3yz;
    abcd[iGrid*1323+691] = 2.0E0*I_ESP_I2y4z_Hx3yz_a-1*I_ESP_G4z_Hx3yz;
    abcd[iGrid*1323+692] = 2.0E0*I_ESP_Iy5z_Hx3yz_a;
    abcd[iGrid*1323+693] = 2.0E0*I_ESP_I5xy_Hx2y2z_a;
    abcd[iGrid*1323+694] = 2.0E0*I_ESP_I4x2y_Hx2y2z_a-1*I_ESP_G4x_Hx2y2z;
    abcd[iGrid*1323+695] = 2.0E0*I_ESP_I4xyz_Hx2y2z_a;
    abcd[iGrid*1323+696] = 2.0E0*I_ESP_I3x3y_Hx2y2z_a-2*I_ESP_G3xy_Hx2y2z;
    abcd[iGrid*1323+697] = 2.0E0*I_ESP_I3x2yz_Hx2y2z_a-1*I_ESP_G3xz_Hx2y2z;
    abcd[iGrid*1323+698] = 2.0E0*I_ESP_I3xy2z_Hx2y2z_a;
    abcd[iGrid*1323+699] = 2.0E0*I_ESP_I2x4y_Hx2y2z_a-3*I_ESP_G2x2y_Hx2y2z;
    abcd[iGrid*1323+700] = 2.0E0*I_ESP_I2x3yz_Hx2y2z_a-2*I_ESP_G2xyz_Hx2y2z;
    abcd[iGrid*1323+701] = 2.0E0*I_ESP_I2x2y2z_Hx2y2z_a-1*I_ESP_G2x2z_Hx2y2z;
    abcd[iGrid*1323+702] = 2.0E0*I_ESP_I2xy3z_Hx2y2z_a;
    abcd[iGrid*1323+703] = 2.0E0*I_ESP_Ix5y_Hx2y2z_a-4*I_ESP_Gx3y_Hx2y2z;
    abcd[iGrid*1323+704] = 2.0E0*I_ESP_Ix4yz_Hx2y2z_a-3*I_ESP_Gx2yz_Hx2y2z;
    abcd[iGrid*1323+705] = 2.0E0*I_ESP_Ix3y2z_Hx2y2z_a-2*I_ESP_Gxy2z_Hx2y2z;
    abcd[iGrid*1323+706] = 2.0E0*I_ESP_Ix2y3z_Hx2y2z_a-1*I_ESP_Gx3z_Hx2y2z;
    abcd[iGrid*1323+707] = 2.0E0*I_ESP_Ixy4z_Hx2y2z_a;
    abcd[iGrid*1323+708] = 2.0E0*I_ESP_I6y_Hx2y2z_a-5*I_ESP_G4y_Hx2y2z;
    abcd[iGrid*1323+709] = 2.0E0*I_ESP_I5yz_Hx2y2z_a-4*I_ESP_G3yz_Hx2y2z;
    abcd[iGrid*1323+710] = 2.0E0*I_ESP_I4y2z_Hx2y2z_a-3*I_ESP_G2y2z_Hx2y2z;
    abcd[iGrid*1323+711] = 2.0E0*I_ESP_I3y3z_Hx2y2z_a-2*I_ESP_Gy3z_Hx2y2z;
    abcd[iGrid*1323+712] = 2.0E0*I_ESP_I2y4z_Hx2y2z_a-1*I_ESP_G4z_Hx2y2z;
    abcd[iGrid*1323+713] = 2.0E0*I_ESP_Iy5z_Hx2y2z_a;
    abcd[iGrid*1323+714] = 2.0E0*I_ESP_I5xy_Hxy3z_a;
    abcd[iGrid*1323+715] = 2.0E0*I_ESP_I4x2y_Hxy3z_a-1*I_ESP_G4x_Hxy3z;
    abcd[iGrid*1323+716] = 2.0E0*I_ESP_I4xyz_Hxy3z_a;
    abcd[iGrid*1323+717] = 2.0E0*I_ESP_I3x3y_Hxy3z_a-2*I_ESP_G3xy_Hxy3z;
    abcd[iGrid*1323+718] = 2.0E0*I_ESP_I3x2yz_Hxy3z_a-1*I_ESP_G3xz_Hxy3z;
    abcd[iGrid*1323+719] = 2.0E0*I_ESP_I3xy2z_Hxy3z_a;
    abcd[iGrid*1323+720] = 2.0E0*I_ESP_I2x4y_Hxy3z_a-3*I_ESP_G2x2y_Hxy3z;
    abcd[iGrid*1323+721] = 2.0E0*I_ESP_I2x3yz_Hxy3z_a-2*I_ESP_G2xyz_Hxy3z;
    abcd[iGrid*1323+722] = 2.0E0*I_ESP_I2x2y2z_Hxy3z_a-1*I_ESP_G2x2z_Hxy3z;
    abcd[iGrid*1323+723] = 2.0E0*I_ESP_I2xy3z_Hxy3z_a;
    abcd[iGrid*1323+724] = 2.0E0*I_ESP_Ix5y_Hxy3z_a-4*I_ESP_Gx3y_Hxy3z;
    abcd[iGrid*1323+725] = 2.0E0*I_ESP_Ix4yz_Hxy3z_a-3*I_ESP_Gx2yz_Hxy3z;
    abcd[iGrid*1323+726] = 2.0E0*I_ESP_Ix3y2z_Hxy3z_a-2*I_ESP_Gxy2z_Hxy3z;
    abcd[iGrid*1323+727] = 2.0E0*I_ESP_Ix2y3z_Hxy3z_a-1*I_ESP_Gx3z_Hxy3z;
    abcd[iGrid*1323+728] = 2.0E0*I_ESP_Ixy4z_Hxy3z_a;
    abcd[iGrid*1323+729] = 2.0E0*I_ESP_I6y_Hxy3z_a-5*I_ESP_G4y_Hxy3z;
    abcd[iGrid*1323+730] = 2.0E0*I_ESP_I5yz_Hxy3z_a-4*I_ESP_G3yz_Hxy3z;
    abcd[iGrid*1323+731] = 2.0E0*I_ESP_I4y2z_Hxy3z_a-3*I_ESP_G2y2z_Hxy3z;
    abcd[iGrid*1323+732] = 2.0E0*I_ESP_I3y3z_Hxy3z_a-2*I_ESP_Gy3z_Hxy3z;
    abcd[iGrid*1323+733] = 2.0E0*I_ESP_I2y4z_Hxy3z_a-1*I_ESP_G4z_Hxy3z;
    abcd[iGrid*1323+734] = 2.0E0*I_ESP_Iy5z_Hxy3z_a;
    abcd[iGrid*1323+735] = 2.0E0*I_ESP_I5xy_Hx4z_a;
    abcd[iGrid*1323+736] = 2.0E0*I_ESP_I4x2y_Hx4z_a-1*I_ESP_G4x_Hx4z;
    abcd[iGrid*1323+737] = 2.0E0*I_ESP_I4xyz_Hx4z_a;
    abcd[iGrid*1323+738] = 2.0E0*I_ESP_I3x3y_Hx4z_a-2*I_ESP_G3xy_Hx4z;
    abcd[iGrid*1323+739] = 2.0E0*I_ESP_I3x2yz_Hx4z_a-1*I_ESP_G3xz_Hx4z;
    abcd[iGrid*1323+740] = 2.0E0*I_ESP_I3xy2z_Hx4z_a;
    abcd[iGrid*1323+741] = 2.0E0*I_ESP_I2x4y_Hx4z_a-3*I_ESP_G2x2y_Hx4z;
    abcd[iGrid*1323+742] = 2.0E0*I_ESP_I2x3yz_Hx4z_a-2*I_ESP_G2xyz_Hx4z;
    abcd[iGrid*1323+743] = 2.0E0*I_ESP_I2x2y2z_Hx4z_a-1*I_ESP_G2x2z_Hx4z;
    abcd[iGrid*1323+744] = 2.0E0*I_ESP_I2xy3z_Hx4z_a;
    abcd[iGrid*1323+745] = 2.0E0*I_ESP_Ix5y_Hx4z_a-4*I_ESP_Gx3y_Hx4z;
    abcd[iGrid*1323+746] = 2.0E0*I_ESP_Ix4yz_Hx4z_a-3*I_ESP_Gx2yz_Hx4z;
    abcd[iGrid*1323+747] = 2.0E0*I_ESP_Ix3y2z_Hx4z_a-2*I_ESP_Gxy2z_Hx4z;
    abcd[iGrid*1323+748] = 2.0E0*I_ESP_Ix2y3z_Hx4z_a-1*I_ESP_Gx3z_Hx4z;
    abcd[iGrid*1323+749] = 2.0E0*I_ESP_Ixy4z_Hx4z_a;
    abcd[iGrid*1323+750] = 2.0E0*I_ESP_I6y_Hx4z_a-5*I_ESP_G4y_Hx4z;
    abcd[iGrid*1323+751] = 2.0E0*I_ESP_I5yz_Hx4z_a-4*I_ESP_G3yz_Hx4z;
    abcd[iGrid*1323+752] = 2.0E0*I_ESP_I4y2z_Hx4z_a-3*I_ESP_G2y2z_Hx4z;
    abcd[iGrid*1323+753] = 2.0E0*I_ESP_I3y3z_Hx4z_a-2*I_ESP_Gy3z_Hx4z;
    abcd[iGrid*1323+754] = 2.0E0*I_ESP_I2y4z_Hx4z_a-1*I_ESP_G4z_Hx4z;
    abcd[iGrid*1323+755] = 2.0E0*I_ESP_Iy5z_Hx4z_a;
    abcd[iGrid*1323+756] = 2.0E0*I_ESP_I5xy_H5y_a;
    abcd[iGrid*1323+757] = 2.0E0*I_ESP_I4x2y_H5y_a-1*I_ESP_G4x_H5y;
    abcd[iGrid*1323+758] = 2.0E0*I_ESP_I4xyz_H5y_a;
    abcd[iGrid*1323+759] = 2.0E0*I_ESP_I3x3y_H5y_a-2*I_ESP_G3xy_H5y;
    abcd[iGrid*1323+760] = 2.0E0*I_ESP_I3x2yz_H5y_a-1*I_ESP_G3xz_H5y;
    abcd[iGrid*1323+761] = 2.0E0*I_ESP_I3xy2z_H5y_a;
    abcd[iGrid*1323+762] = 2.0E0*I_ESP_I2x4y_H5y_a-3*I_ESP_G2x2y_H5y;
    abcd[iGrid*1323+763] = 2.0E0*I_ESP_I2x3yz_H5y_a-2*I_ESP_G2xyz_H5y;
    abcd[iGrid*1323+764] = 2.0E0*I_ESP_I2x2y2z_H5y_a-1*I_ESP_G2x2z_H5y;
    abcd[iGrid*1323+765] = 2.0E0*I_ESP_I2xy3z_H5y_a;
    abcd[iGrid*1323+766] = 2.0E0*I_ESP_Ix5y_H5y_a-4*I_ESP_Gx3y_H5y;
    abcd[iGrid*1323+767] = 2.0E0*I_ESP_Ix4yz_H5y_a-3*I_ESP_Gx2yz_H5y;
    abcd[iGrid*1323+768] = 2.0E0*I_ESP_Ix3y2z_H5y_a-2*I_ESP_Gxy2z_H5y;
    abcd[iGrid*1323+769] = 2.0E0*I_ESP_Ix2y3z_H5y_a-1*I_ESP_Gx3z_H5y;
    abcd[iGrid*1323+770] = 2.0E0*I_ESP_Ixy4z_H5y_a;
    abcd[iGrid*1323+771] = 2.0E0*I_ESP_I6y_H5y_a-5*I_ESP_G4y_H5y;
    abcd[iGrid*1323+772] = 2.0E0*I_ESP_I5yz_H5y_a-4*I_ESP_G3yz_H5y;
    abcd[iGrid*1323+773] = 2.0E0*I_ESP_I4y2z_H5y_a-3*I_ESP_G2y2z_H5y;
    abcd[iGrid*1323+774] = 2.0E0*I_ESP_I3y3z_H5y_a-2*I_ESP_Gy3z_H5y;
    abcd[iGrid*1323+775] = 2.0E0*I_ESP_I2y4z_H5y_a-1*I_ESP_G4z_H5y;
    abcd[iGrid*1323+776] = 2.0E0*I_ESP_Iy5z_H5y_a;
    abcd[iGrid*1323+777] = 2.0E0*I_ESP_I5xy_H4yz_a;
    abcd[iGrid*1323+778] = 2.0E0*I_ESP_I4x2y_H4yz_a-1*I_ESP_G4x_H4yz;
    abcd[iGrid*1323+779] = 2.0E0*I_ESP_I4xyz_H4yz_a;
    abcd[iGrid*1323+780] = 2.0E0*I_ESP_I3x3y_H4yz_a-2*I_ESP_G3xy_H4yz;
    abcd[iGrid*1323+781] = 2.0E0*I_ESP_I3x2yz_H4yz_a-1*I_ESP_G3xz_H4yz;
    abcd[iGrid*1323+782] = 2.0E0*I_ESP_I3xy2z_H4yz_a;
    abcd[iGrid*1323+783] = 2.0E0*I_ESP_I2x4y_H4yz_a-3*I_ESP_G2x2y_H4yz;
    abcd[iGrid*1323+784] = 2.0E0*I_ESP_I2x3yz_H4yz_a-2*I_ESP_G2xyz_H4yz;
    abcd[iGrid*1323+785] = 2.0E0*I_ESP_I2x2y2z_H4yz_a-1*I_ESP_G2x2z_H4yz;
    abcd[iGrid*1323+786] = 2.0E0*I_ESP_I2xy3z_H4yz_a;
    abcd[iGrid*1323+787] = 2.0E0*I_ESP_Ix5y_H4yz_a-4*I_ESP_Gx3y_H4yz;
    abcd[iGrid*1323+788] = 2.0E0*I_ESP_Ix4yz_H4yz_a-3*I_ESP_Gx2yz_H4yz;
    abcd[iGrid*1323+789] = 2.0E0*I_ESP_Ix3y2z_H4yz_a-2*I_ESP_Gxy2z_H4yz;
    abcd[iGrid*1323+790] = 2.0E0*I_ESP_Ix2y3z_H4yz_a-1*I_ESP_Gx3z_H4yz;
    abcd[iGrid*1323+791] = 2.0E0*I_ESP_Ixy4z_H4yz_a;
    abcd[iGrid*1323+792] = 2.0E0*I_ESP_I6y_H4yz_a-5*I_ESP_G4y_H4yz;
    abcd[iGrid*1323+793] = 2.0E0*I_ESP_I5yz_H4yz_a-4*I_ESP_G3yz_H4yz;
    abcd[iGrid*1323+794] = 2.0E0*I_ESP_I4y2z_H4yz_a-3*I_ESP_G2y2z_H4yz;
    abcd[iGrid*1323+795] = 2.0E0*I_ESP_I3y3z_H4yz_a-2*I_ESP_Gy3z_H4yz;
    abcd[iGrid*1323+796] = 2.0E0*I_ESP_I2y4z_H4yz_a-1*I_ESP_G4z_H4yz;
    abcd[iGrid*1323+797] = 2.0E0*I_ESP_Iy5z_H4yz_a;
    abcd[iGrid*1323+798] = 2.0E0*I_ESP_I5xy_H3y2z_a;
    abcd[iGrid*1323+799] = 2.0E0*I_ESP_I4x2y_H3y2z_a-1*I_ESP_G4x_H3y2z;
    abcd[iGrid*1323+800] = 2.0E0*I_ESP_I4xyz_H3y2z_a;
    abcd[iGrid*1323+801] = 2.0E0*I_ESP_I3x3y_H3y2z_a-2*I_ESP_G3xy_H3y2z;
    abcd[iGrid*1323+802] = 2.0E0*I_ESP_I3x2yz_H3y2z_a-1*I_ESP_G3xz_H3y2z;
    abcd[iGrid*1323+803] = 2.0E0*I_ESP_I3xy2z_H3y2z_a;
    abcd[iGrid*1323+804] = 2.0E0*I_ESP_I2x4y_H3y2z_a-3*I_ESP_G2x2y_H3y2z;
    abcd[iGrid*1323+805] = 2.0E0*I_ESP_I2x3yz_H3y2z_a-2*I_ESP_G2xyz_H3y2z;
    abcd[iGrid*1323+806] = 2.0E0*I_ESP_I2x2y2z_H3y2z_a-1*I_ESP_G2x2z_H3y2z;
    abcd[iGrid*1323+807] = 2.0E0*I_ESP_I2xy3z_H3y2z_a;
    abcd[iGrid*1323+808] = 2.0E0*I_ESP_Ix5y_H3y2z_a-4*I_ESP_Gx3y_H3y2z;
    abcd[iGrid*1323+809] = 2.0E0*I_ESP_Ix4yz_H3y2z_a-3*I_ESP_Gx2yz_H3y2z;
    abcd[iGrid*1323+810] = 2.0E0*I_ESP_Ix3y2z_H3y2z_a-2*I_ESP_Gxy2z_H3y2z;
    abcd[iGrid*1323+811] = 2.0E0*I_ESP_Ix2y3z_H3y2z_a-1*I_ESP_Gx3z_H3y2z;
    abcd[iGrid*1323+812] = 2.0E0*I_ESP_Ixy4z_H3y2z_a;
    abcd[iGrid*1323+813] = 2.0E0*I_ESP_I6y_H3y2z_a-5*I_ESP_G4y_H3y2z;
    abcd[iGrid*1323+814] = 2.0E0*I_ESP_I5yz_H3y2z_a-4*I_ESP_G3yz_H3y2z;
    abcd[iGrid*1323+815] = 2.0E0*I_ESP_I4y2z_H3y2z_a-3*I_ESP_G2y2z_H3y2z;
    abcd[iGrid*1323+816] = 2.0E0*I_ESP_I3y3z_H3y2z_a-2*I_ESP_Gy3z_H3y2z;
    abcd[iGrid*1323+817] = 2.0E0*I_ESP_I2y4z_H3y2z_a-1*I_ESP_G4z_H3y2z;
    abcd[iGrid*1323+818] = 2.0E0*I_ESP_Iy5z_H3y2z_a;
    abcd[iGrid*1323+819] = 2.0E0*I_ESP_I5xy_H2y3z_a;
    abcd[iGrid*1323+820] = 2.0E0*I_ESP_I4x2y_H2y3z_a-1*I_ESP_G4x_H2y3z;
    abcd[iGrid*1323+821] = 2.0E0*I_ESP_I4xyz_H2y3z_a;
    abcd[iGrid*1323+822] = 2.0E0*I_ESP_I3x3y_H2y3z_a-2*I_ESP_G3xy_H2y3z;
    abcd[iGrid*1323+823] = 2.0E0*I_ESP_I3x2yz_H2y3z_a-1*I_ESP_G3xz_H2y3z;
    abcd[iGrid*1323+824] = 2.0E0*I_ESP_I3xy2z_H2y3z_a;
    abcd[iGrid*1323+825] = 2.0E0*I_ESP_I2x4y_H2y3z_a-3*I_ESP_G2x2y_H2y3z;
    abcd[iGrid*1323+826] = 2.0E0*I_ESP_I2x3yz_H2y3z_a-2*I_ESP_G2xyz_H2y3z;
    abcd[iGrid*1323+827] = 2.0E0*I_ESP_I2x2y2z_H2y3z_a-1*I_ESP_G2x2z_H2y3z;
    abcd[iGrid*1323+828] = 2.0E0*I_ESP_I2xy3z_H2y3z_a;
    abcd[iGrid*1323+829] = 2.0E0*I_ESP_Ix5y_H2y3z_a-4*I_ESP_Gx3y_H2y3z;
    abcd[iGrid*1323+830] = 2.0E0*I_ESP_Ix4yz_H2y3z_a-3*I_ESP_Gx2yz_H2y3z;
    abcd[iGrid*1323+831] = 2.0E0*I_ESP_Ix3y2z_H2y3z_a-2*I_ESP_Gxy2z_H2y3z;
    abcd[iGrid*1323+832] = 2.0E0*I_ESP_Ix2y3z_H2y3z_a-1*I_ESP_Gx3z_H2y3z;
    abcd[iGrid*1323+833] = 2.0E0*I_ESP_Ixy4z_H2y3z_a;
    abcd[iGrid*1323+834] = 2.0E0*I_ESP_I6y_H2y3z_a-5*I_ESP_G4y_H2y3z;
    abcd[iGrid*1323+835] = 2.0E0*I_ESP_I5yz_H2y3z_a-4*I_ESP_G3yz_H2y3z;
    abcd[iGrid*1323+836] = 2.0E0*I_ESP_I4y2z_H2y3z_a-3*I_ESP_G2y2z_H2y3z;
    abcd[iGrid*1323+837] = 2.0E0*I_ESP_I3y3z_H2y3z_a-2*I_ESP_Gy3z_H2y3z;
    abcd[iGrid*1323+838] = 2.0E0*I_ESP_I2y4z_H2y3z_a-1*I_ESP_G4z_H2y3z;
    abcd[iGrid*1323+839] = 2.0E0*I_ESP_Iy5z_H2y3z_a;
    abcd[iGrid*1323+840] = 2.0E0*I_ESP_I5xy_Hy4z_a;
    abcd[iGrid*1323+841] = 2.0E0*I_ESP_I4x2y_Hy4z_a-1*I_ESP_G4x_Hy4z;
    abcd[iGrid*1323+842] = 2.0E0*I_ESP_I4xyz_Hy4z_a;
    abcd[iGrid*1323+843] = 2.0E0*I_ESP_I3x3y_Hy4z_a-2*I_ESP_G3xy_Hy4z;
    abcd[iGrid*1323+844] = 2.0E0*I_ESP_I3x2yz_Hy4z_a-1*I_ESP_G3xz_Hy4z;
    abcd[iGrid*1323+845] = 2.0E0*I_ESP_I3xy2z_Hy4z_a;
    abcd[iGrid*1323+846] = 2.0E0*I_ESP_I2x4y_Hy4z_a-3*I_ESP_G2x2y_Hy4z;
    abcd[iGrid*1323+847] = 2.0E0*I_ESP_I2x3yz_Hy4z_a-2*I_ESP_G2xyz_Hy4z;
    abcd[iGrid*1323+848] = 2.0E0*I_ESP_I2x2y2z_Hy4z_a-1*I_ESP_G2x2z_Hy4z;
    abcd[iGrid*1323+849] = 2.0E0*I_ESP_I2xy3z_Hy4z_a;
    abcd[iGrid*1323+850] = 2.0E0*I_ESP_Ix5y_Hy4z_a-4*I_ESP_Gx3y_Hy4z;
    abcd[iGrid*1323+851] = 2.0E0*I_ESP_Ix4yz_Hy4z_a-3*I_ESP_Gx2yz_Hy4z;
    abcd[iGrid*1323+852] = 2.0E0*I_ESP_Ix3y2z_Hy4z_a-2*I_ESP_Gxy2z_Hy4z;
    abcd[iGrid*1323+853] = 2.0E0*I_ESP_Ix2y3z_Hy4z_a-1*I_ESP_Gx3z_Hy4z;
    abcd[iGrid*1323+854] = 2.0E0*I_ESP_Ixy4z_Hy4z_a;
    abcd[iGrid*1323+855] = 2.0E0*I_ESP_I6y_Hy4z_a-5*I_ESP_G4y_Hy4z;
    abcd[iGrid*1323+856] = 2.0E0*I_ESP_I5yz_Hy4z_a-4*I_ESP_G3yz_Hy4z;
    abcd[iGrid*1323+857] = 2.0E0*I_ESP_I4y2z_Hy4z_a-3*I_ESP_G2y2z_Hy4z;
    abcd[iGrid*1323+858] = 2.0E0*I_ESP_I3y3z_Hy4z_a-2*I_ESP_Gy3z_Hy4z;
    abcd[iGrid*1323+859] = 2.0E0*I_ESP_I2y4z_Hy4z_a-1*I_ESP_G4z_Hy4z;
    abcd[iGrid*1323+860] = 2.0E0*I_ESP_Iy5z_Hy4z_a;
    abcd[iGrid*1323+861] = 2.0E0*I_ESP_I5xy_H5z_a;
    abcd[iGrid*1323+862] = 2.0E0*I_ESP_I4x2y_H5z_a-1*I_ESP_G4x_H5z;
    abcd[iGrid*1323+863] = 2.0E0*I_ESP_I4xyz_H5z_a;
    abcd[iGrid*1323+864] = 2.0E0*I_ESP_I3x3y_H5z_a-2*I_ESP_G3xy_H5z;
    abcd[iGrid*1323+865] = 2.0E0*I_ESP_I3x2yz_H5z_a-1*I_ESP_G3xz_H5z;
    abcd[iGrid*1323+866] = 2.0E0*I_ESP_I3xy2z_H5z_a;
    abcd[iGrid*1323+867] = 2.0E0*I_ESP_I2x4y_H5z_a-3*I_ESP_G2x2y_H5z;
    abcd[iGrid*1323+868] = 2.0E0*I_ESP_I2x3yz_H5z_a-2*I_ESP_G2xyz_H5z;
    abcd[iGrid*1323+869] = 2.0E0*I_ESP_I2x2y2z_H5z_a-1*I_ESP_G2x2z_H5z;
    abcd[iGrid*1323+870] = 2.0E0*I_ESP_I2xy3z_H5z_a;
    abcd[iGrid*1323+871] = 2.0E0*I_ESP_Ix5y_H5z_a-4*I_ESP_Gx3y_H5z;
    abcd[iGrid*1323+872] = 2.0E0*I_ESP_Ix4yz_H5z_a-3*I_ESP_Gx2yz_H5z;
    abcd[iGrid*1323+873] = 2.0E0*I_ESP_Ix3y2z_H5z_a-2*I_ESP_Gxy2z_H5z;
    abcd[iGrid*1323+874] = 2.0E0*I_ESP_Ix2y3z_H5z_a-1*I_ESP_Gx3z_H5z;
    abcd[iGrid*1323+875] = 2.0E0*I_ESP_Ixy4z_H5z_a;
    abcd[iGrid*1323+876] = 2.0E0*I_ESP_I6y_H5z_a-5*I_ESP_G4y_H5z;
    abcd[iGrid*1323+877] = 2.0E0*I_ESP_I5yz_H5z_a-4*I_ESP_G3yz_H5z;
    abcd[iGrid*1323+878] = 2.0E0*I_ESP_I4y2z_H5z_a-3*I_ESP_G2y2z_H5z;
    abcd[iGrid*1323+879] = 2.0E0*I_ESP_I3y3z_H5z_a-2*I_ESP_Gy3z_H5z;
    abcd[iGrid*1323+880] = 2.0E0*I_ESP_I2y4z_H5z_a-1*I_ESP_G4z_H5z;
    abcd[iGrid*1323+881] = 2.0E0*I_ESP_Iy5z_H5z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_H_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_H_a
     * RHS shell quartet name: SQ_ESP_G_H
     ************************************************************/
    abcd[iGrid*1323+882] = 2.0E0*I_ESP_I5xz_H5x_a;
    abcd[iGrid*1323+883] = 2.0E0*I_ESP_I4xyz_H5x_a;
    abcd[iGrid*1323+884] = 2.0E0*I_ESP_I4x2z_H5x_a-1*I_ESP_G4x_H5x;
    abcd[iGrid*1323+885] = 2.0E0*I_ESP_I3x2yz_H5x_a;
    abcd[iGrid*1323+886] = 2.0E0*I_ESP_I3xy2z_H5x_a-1*I_ESP_G3xy_H5x;
    abcd[iGrid*1323+887] = 2.0E0*I_ESP_I3x3z_H5x_a-2*I_ESP_G3xz_H5x;
    abcd[iGrid*1323+888] = 2.0E0*I_ESP_I2x3yz_H5x_a;
    abcd[iGrid*1323+889] = 2.0E0*I_ESP_I2x2y2z_H5x_a-1*I_ESP_G2x2y_H5x;
    abcd[iGrid*1323+890] = 2.0E0*I_ESP_I2xy3z_H5x_a-2*I_ESP_G2xyz_H5x;
    abcd[iGrid*1323+891] = 2.0E0*I_ESP_I2x4z_H5x_a-3*I_ESP_G2x2z_H5x;
    abcd[iGrid*1323+892] = 2.0E0*I_ESP_Ix4yz_H5x_a;
    abcd[iGrid*1323+893] = 2.0E0*I_ESP_Ix3y2z_H5x_a-1*I_ESP_Gx3y_H5x;
    abcd[iGrid*1323+894] = 2.0E0*I_ESP_Ix2y3z_H5x_a-2*I_ESP_Gx2yz_H5x;
    abcd[iGrid*1323+895] = 2.0E0*I_ESP_Ixy4z_H5x_a-3*I_ESP_Gxy2z_H5x;
    abcd[iGrid*1323+896] = 2.0E0*I_ESP_Ix5z_H5x_a-4*I_ESP_Gx3z_H5x;
    abcd[iGrid*1323+897] = 2.0E0*I_ESP_I5yz_H5x_a;
    abcd[iGrid*1323+898] = 2.0E0*I_ESP_I4y2z_H5x_a-1*I_ESP_G4y_H5x;
    abcd[iGrid*1323+899] = 2.0E0*I_ESP_I3y3z_H5x_a-2*I_ESP_G3yz_H5x;
    abcd[iGrid*1323+900] = 2.0E0*I_ESP_I2y4z_H5x_a-3*I_ESP_G2y2z_H5x;
    abcd[iGrid*1323+901] = 2.0E0*I_ESP_Iy5z_H5x_a-4*I_ESP_Gy3z_H5x;
    abcd[iGrid*1323+902] = 2.0E0*I_ESP_I6z_H5x_a-5*I_ESP_G4z_H5x;
    abcd[iGrid*1323+903] = 2.0E0*I_ESP_I5xz_H4xy_a;
    abcd[iGrid*1323+904] = 2.0E0*I_ESP_I4xyz_H4xy_a;
    abcd[iGrid*1323+905] = 2.0E0*I_ESP_I4x2z_H4xy_a-1*I_ESP_G4x_H4xy;
    abcd[iGrid*1323+906] = 2.0E0*I_ESP_I3x2yz_H4xy_a;
    abcd[iGrid*1323+907] = 2.0E0*I_ESP_I3xy2z_H4xy_a-1*I_ESP_G3xy_H4xy;
    abcd[iGrid*1323+908] = 2.0E0*I_ESP_I3x3z_H4xy_a-2*I_ESP_G3xz_H4xy;
    abcd[iGrid*1323+909] = 2.0E0*I_ESP_I2x3yz_H4xy_a;
    abcd[iGrid*1323+910] = 2.0E0*I_ESP_I2x2y2z_H4xy_a-1*I_ESP_G2x2y_H4xy;
    abcd[iGrid*1323+911] = 2.0E0*I_ESP_I2xy3z_H4xy_a-2*I_ESP_G2xyz_H4xy;
    abcd[iGrid*1323+912] = 2.0E0*I_ESP_I2x4z_H4xy_a-3*I_ESP_G2x2z_H4xy;
    abcd[iGrid*1323+913] = 2.0E0*I_ESP_Ix4yz_H4xy_a;
    abcd[iGrid*1323+914] = 2.0E0*I_ESP_Ix3y2z_H4xy_a-1*I_ESP_Gx3y_H4xy;
    abcd[iGrid*1323+915] = 2.0E0*I_ESP_Ix2y3z_H4xy_a-2*I_ESP_Gx2yz_H4xy;
    abcd[iGrid*1323+916] = 2.0E0*I_ESP_Ixy4z_H4xy_a-3*I_ESP_Gxy2z_H4xy;
    abcd[iGrid*1323+917] = 2.0E0*I_ESP_Ix5z_H4xy_a-4*I_ESP_Gx3z_H4xy;
    abcd[iGrid*1323+918] = 2.0E0*I_ESP_I5yz_H4xy_a;
    abcd[iGrid*1323+919] = 2.0E0*I_ESP_I4y2z_H4xy_a-1*I_ESP_G4y_H4xy;
    abcd[iGrid*1323+920] = 2.0E0*I_ESP_I3y3z_H4xy_a-2*I_ESP_G3yz_H4xy;
    abcd[iGrid*1323+921] = 2.0E0*I_ESP_I2y4z_H4xy_a-3*I_ESP_G2y2z_H4xy;
    abcd[iGrid*1323+922] = 2.0E0*I_ESP_Iy5z_H4xy_a-4*I_ESP_Gy3z_H4xy;
    abcd[iGrid*1323+923] = 2.0E0*I_ESP_I6z_H4xy_a-5*I_ESP_G4z_H4xy;
    abcd[iGrid*1323+924] = 2.0E0*I_ESP_I5xz_H4xz_a;
    abcd[iGrid*1323+925] = 2.0E0*I_ESP_I4xyz_H4xz_a;
    abcd[iGrid*1323+926] = 2.0E0*I_ESP_I4x2z_H4xz_a-1*I_ESP_G4x_H4xz;
    abcd[iGrid*1323+927] = 2.0E0*I_ESP_I3x2yz_H4xz_a;
    abcd[iGrid*1323+928] = 2.0E0*I_ESP_I3xy2z_H4xz_a-1*I_ESP_G3xy_H4xz;
    abcd[iGrid*1323+929] = 2.0E0*I_ESP_I3x3z_H4xz_a-2*I_ESP_G3xz_H4xz;
    abcd[iGrid*1323+930] = 2.0E0*I_ESP_I2x3yz_H4xz_a;
    abcd[iGrid*1323+931] = 2.0E0*I_ESP_I2x2y2z_H4xz_a-1*I_ESP_G2x2y_H4xz;
    abcd[iGrid*1323+932] = 2.0E0*I_ESP_I2xy3z_H4xz_a-2*I_ESP_G2xyz_H4xz;
    abcd[iGrid*1323+933] = 2.0E0*I_ESP_I2x4z_H4xz_a-3*I_ESP_G2x2z_H4xz;
    abcd[iGrid*1323+934] = 2.0E0*I_ESP_Ix4yz_H4xz_a;
    abcd[iGrid*1323+935] = 2.0E0*I_ESP_Ix3y2z_H4xz_a-1*I_ESP_Gx3y_H4xz;
    abcd[iGrid*1323+936] = 2.0E0*I_ESP_Ix2y3z_H4xz_a-2*I_ESP_Gx2yz_H4xz;
    abcd[iGrid*1323+937] = 2.0E0*I_ESP_Ixy4z_H4xz_a-3*I_ESP_Gxy2z_H4xz;
    abcd[iGrid*1323+938] = 2.0E0*I_ESP_Ix5z_H4xz_a-4*I_ESP_Gx3z_H4xz;
    abcd[iGrid*1323+939] = 2.0E0*I_ESP_I5yz_H4xz_a;
    abcd[iGrid*1323+940] = 2.0E0*I_ESP_I4y2z_H4xz_a-1*I_ESP_G4y_H4xz;
    abcd[iGrid*1323+941] = 2.0E0*I_ESP_I3y3z_H4xz_a-2*I_ESP_G3yz_H4xz;
    abcd[iGrid*1323+942] = 2.0E0*I_ESP_I2y4z_H4xz_a-3*I_ESP_G2y2z_H4xz;
    abcd[iGrid*1323+943] = 2.0E0*I_ESP_Iy5z_H4xz_a-4*I_ESP_Gy3z_H4xz;
    abcd[iGrid*1323+944] = 2.0E0*I_ESP_I6z_H4xz_a-5*I_ESP_G4z_H4xz;
    abcd[iGrid*1323+945] = 2.0E0*I_ESP_I5xz_H3x2y_a;
    abcd[iGrid*1323+946] = 2.0E0*I_ESP_I4xyz_H3x2y_a;
    abcd[iGrid*1323+947] = 2.0E0*I_ESP_I4x2z_H3x2y_a-1*I_ESP_G4x_H3x2y;
    abcd[iGrid*1323+948] = 2.0E0*I_ESP_I3x2yz_H3x2y_a;
    abcd[iGrid*1323+949] = 2.0E0*I_ESP_I3xy2z_H3x2y_a-1*I_ESP_G3xy_H3x2y;
    abcd[iGrid*1323+950] = 2.0E0*I_ESP_I3x3z_H3x2y_a-2*I_ESP_G3xz_H3x2y;
    abcd[iGrid*1323+951] = 2.0E0*I_ESP_I2x3yz_H3x2y_a;
    abcd[iGrid*1323+952] = 2.0E0*I_ESP_I2x2y2z_H3x2y_a-1*I_ESP_G2x2y_H3x2y;
    abcd[iGrid*1323+953] = 2.0E0*I_ESP_I2xy3z_H3x2y_a-2*I_ESP_G2xyz_H3x2y;
    abcd[iGrid*1323+954] = 2.0E0*I_ESP_I2x4z_H3x2y_a-3*I_ESP_G2x2z_H3x2y;
    abcd[iGrid*1323+955] = 2.0E0*I_ESP_Ix4yz_H3x2y_a;
    abcd[iGrid*1323+956] = 2.0E0*I_ESP_Ix3y2z_H3x2y_a-1*I_ESP_Gx3y_H3x2y;
    abcd[iGrid*1323+957] = 2.0E0*I_ESP_Ix2y3z_H3x2y_a-2*I_ESP_Gx2yz_H3x2y;
    abcd[iGrid*1323+958] = 2.0E0*I_ESP_Ixy4z_H3x2y_a-3*I_ESP_Gxy2z_H3x2y;
    abcd[iGrid*1323+959] = 2.0E0*I_ESP_Ix5z_H3x2y_a-4*I_ESP_Gx3z_H3x2y;
    abcd[iGrid*1323+960] = 2.0E0*I_ESP_I5yz_H3x2y_a;
    abcd[iGrid*1323+961] = 2.0E0*I_ESP_I4y2z_H3x2y_a-1*I_ESP_G4y_H3x2y;
    abcd[iGrid*1323+962] = 2.0E0*I_ESP_I3y3z_H3x2y_a-2*I_ESP_G3yz_H3x2y;
    abcd[iGrid*1323+963] = 2.0E0*I_ESP_I2y4z_H3x2y_a-3*I_ESP_G2y2z_H3x2y;
    abcd[iGrid*1323+964] = 2.0E0*I_ESP_Iy5z_H3x2y_a-4*I_ESP_Gy3z_H3x2y;
    abcd[iGrid*1323+965] = 2.0E0*I_ESP_I6z_H3x2y_a-5*I_ESP_G4z_H3x2y;
    abcd[iGrid*1323+966] = 2.0E0*I_ESP_I5xz_H3xyz_a;
    abcd[iGrid*1323+967] = 2.0E0*I_ESP_I4xyz_H3xyz_a;
    abcd[iGrid*1323+968] = 2.0E0*I_ESP_I4x2z_H3xyz_a-1*I_ESP_G4x_H3xyz;
    abcd[iGrid*1323+969] = 2.0E0*I_ESP_I3x2yz_H3xyz_a;
    abcd[iGrid*1323+970] = 2.0E0*I_ESP_I3xy2z_H3xyz_a-1*I_ESP_G3xy_H3xyz;
    abcd[iGrid*1323+971] = 2.0E0*I_ESP_I3x3z_H3xyz_a-2*I_ESP_G3xz_H3xyz;
    abcd[iGrid*1323+972] = 2.0E0*I_ESP_I2x3yz_H3xyz_a;
    abcd[iGrid*1323+973] = 2.0E0*I_ESP_I2x2y2z_H3xyz_a-1*I_ESP_G2x2y_H3xyz;
    abcd[iGrid*1323+974] = 2.0E0*I_ESP_I2xy3z_H3xyz_a-2*I_ESP_G2xyz_H3xyz;
    abcd[iGrid*1323+975] = 2.0E0*I_ESP_I2x4z_H3xyz_a-3*I_ESP_G2x2z_H3xyz;
    abcd[iGrid*1323+976] = 2.0E0*I_ESP_Ix4yz_H3xyz_a;
    abcd[iGrid*1323+977] = 2.0E0*I_ESP_Ix3y2z_H3xyz_a-1*I_ESP_Gx3y_H3xyz;
    abcd[iGrid*1323+978] = 2.0E0*I_ESP_Ix2y3z_H3xyz_a-2*I_ESP_Gx2yz_H3xyz;
    abcd[iGrid*1323+979] = 2.0E0*I_ESP_Ixy4z_H3xyz_a-3*I_ESP_Gxy2z_H3xyz;
    abcd[iGrid*1323+980] = 2.0E0*I_ESP_Ix5z_H3xyz_a-4*I_ESP_Gx3z_H3xyz;
    abcd[iGrid*1323+981] = 2.0E0*I_ESP_I5yz_H3xyz_a;
    abcd[iGrid*1323+982] = 2.0E0*I_ESP_I4y2z_H3xyz_a-1*I_ESP_G4y_H3xyz;
    abcd[iGrid*1323+983] = 2.0E0*I_ESP_I3y3z_H3xyz_a-2*I_ESP_G3yz_H3xyz;
    abcd[iGrid*1323+984] = 2.0E0*I_ESP_I2y4z_H3xyz_a-3*I_ESP_G2y2z_H3xyz;
    abcd[iGrid*1323+985] = 2.0E0*I_ESP_Iy5z_H3xyz_a-4*I_ESP_Gy3z_H3xyz;
    abcd[iGrid*1323+986] = 2.0E0*I_ESP_I6z_H3xyz_a-5*I_ESP_G4z_H3xyz;
    abcd[iGrid*1323+987] = 2.0E0*I_ESP_I5xz_H3x2z_a;
    abcd[iGrid*1323+988] = 2.0E0*I_ESP_I4xyz_H3x2z_a;
    abcd[iGrid*1323+989] = 2.0E0*I_ESP_I4x2z_H3x2z_a-1*I_ESP_G4x_H3x2z;
    abcd[iGrid*1323+990] = 2.0E0*I_ESP_I3x2yz_H3x2z_a;
    abcd[iGrid*1323+991] = 2.0E0*I_ESP_I3xy2z_H3x2z_a-1*I_ESP_G3xy_H3x2z;
    abcd[iGrid*1323+992] = 2.0E0*I_ESP_I3x3z_H3x2z_a-2*I_ESP_G3xz_H3x2z;
    abcd[iGrid*1323+993] = 2.0E0*I_ESP_I2x3yz_H3x2z_a;
    abcd[iGrid*1323+994] = 2.0E0*I_ESP_I2x2y2z_H3x2z_a-1*I_ESP_G2x2y_H3x2z;
    abcd[iGrid*1323+995] = 2.0E0*I_ESP_I2xy3z_H3x2z_a-2*I_ESP_G2xyz_H3x2z;
    abcd[iGrid*1323+996] = 2.0E0*I_ESP_I2x4z_H3x2z_a-3*I_ESP_G2x2z_H3x2z;
    abcd[iGrid*1323+997] = 2.0E0*I_ESP_Ix4yz_H3x2z_a;
    abcd[iGrid*1323+998] = 2.0E0*I_ESP_Ix3y2z_H3x2z_a-1*I_ESP_Gx3y_H3x2z;
    abcd[iGrid*1323+999] = 2.0E0*I_ESP_Ix2y3z_H3x2z_a-2*I_ESP_Gx2yz_H3x2z;
    abcd[iGrid*1323+1000] = 2.0E0*I_ESP_Ixy4z_H3x2z_a-3*I_ESP_Gxy2z_H3x2z;
    abcd[iGrid*1323+1001] = 2.0E0*I_ESP_Ix5z_H3x2z_a-4*I_ESP_Gx3z_H3x2z;
    abcd[iGrid*1323+1002] = 2.0E0*I_ESP_I5yz_H3x2z_a;
    abcd[iGrid*1323+1003] = 2.0E0*I_ESP_I4y2z_H3x2z_a-1*I_ESP_G4y_H3x2z;
    abcd[iGrid*1323+1004] = 2.0E0*I_ESP_I3y3z_H3x2z_a-2*I_ESP_G3yz_H3x2z;
    abcd[iGrid*1323+1005] = 2.0E0*I_ESP_I2y4z_H3x2z_a-3*I_ESP_G2y2z_H3x2z;
    abcd[iGrid*1323+1006] = 2.0E0*I_ESP_Iy5z_H3x2z_a-4*I_ESP_Gy3z_H3x2z;
    abcd[iGrid*1323+1007] = 2.0E0*I_ESP_I6z_H3x2z_a-5*I_ESP_G4z_H3x2z;
    abcd[iGrid*1323+1008] = 2.0E0*I_ESP_I5xz_H2x3y_a;
    abcd[iGrid*1323+1009] = 2.0E0*I_ESP_I4xyz_H2x3y_a;
    abcd[iGrid*1323+1010] = 2.0E0*I_ESP_I4x2z_H2x3y_a-1*I_ESP_G4x_H2x3y;
    abcd[iGrid*1323+1011] = 2.0E0*I_ESP_I3x2yz_H2x3y_a;
    abcd[iGrid*1323+1012] = 2.0E0*I_ESP_I3xy2z_H2x3y_a-1*I_ESP_G3xy_H2x3y;
    abcd[iGrid*1323+1013] = 2.0E0*I_ESP_I3x3z_H2x3y_a-2*I_ESP_G3xz_H2x3y;
    abcd[iGrid*1323+1014] = 2.0E0*I_ESP_I2x3yz_H2x3y_a;
    abcd[iGrid*1323+1015] = 2.0E0*I_ESP_I2x2y2z_H2x3y_a-1*I_ESP_G2x2y_H2x3y;
    abcd[iGrid*1323+1016] = 2.0E0*I_ESP_I2xy3z_H2x3y_a-2*I_ESP_G2xyz_H2x3y;
    abcd[iGrid*1323+1017] = 2.0E0*I_ESP_I2x4z_H2x3y_a-3*I_ESP_G2x2z_H2x3y;
    abcd[iGrid*1323+1018] = 2.0E0*I_ESP_Ix4yz_H2x3y_a;
    abcd[iGrid*1323+1019] = 2.0E0*I_ESP_Ix3y2z_H2x3y_a-1*I_ESP_Gx3y_H2x3y;
    abcd[iGrid*1323+1020] = 2.0E0*I_ESP_Ix2y3z_H2x3y_a-2*I_ESP_Gx2yz_H2x3y;
    abcd[iGrid*1323+1021] = 2.0E0*I_ESP_Ixy4z_H2x3y_a-3*I_ESP_Gxy2z_H2x3y;
    abcd[iGrid*1323+1022] = 2.0E0*I_ESP_Ix5z_H2x3y_a-4*I_ESP_Gx3z_H2x3y;
    abcd[iGrid*1323+1023] = 2.0E0*I_ESP_I5yz_H2x3y_a;
    abcd[iGrid*1323+1024] = 2.0E0*I_ESP_I4y2z_H2x3y_a-1*I_ESP_G4y_H2x3y;
    abcd[iGrid*1323+1025] = 2.0E0*I_ESP_I3y3z_H2x3y_a-2*I_ESP_G3yz_H2x3y;
    abcd[iGrid*1323+1026] = 2.0E0*I_ESP_I2y4z_H2x3y_a-3*I_ESP_G2y2z_H2x3y;
    abcd[iGrid*1323+1027] = 2.0E0*I_ESP_Iy5z_H2x3y_a-4*I_ESP_Gy3z_H2x3y;
    abcd[iGrid*1323+1028] = 2.0E0*I_ESP_I6z_H2x3y_a-5*I_ESP_G4z_H2x3y;
    abcd[iGrid*1323+1029] = 2.0E0*I_ESP_I5xz_H2x2yz_a;
    abcd[iGrid*1323+1030] = 2.0E0*I_ESP_I4xyz_H2x2yz_a;
    abcd[iGrid*1323+1031] = 2.0E0*I_ESP_I4x2z_H2x2yz_a-1*I_ESP_G4x_H2x2yz;
    abcd[iGrid*1323+1032] = 2.0E0*I_ESP_I3x2yz_H2x2yz_a;
    abcd[iGrid*1323+1033] = 2.0E0*I_ESP_I3xy2z_H2x2yz_a-1*I_ESP_G3xy_H2x2yz;
    abcd[iGrid*1323+1034] = 2.0E0*I_ESP_I3x3z_H2x2yz_a-2*I_ESP_G3xz_H2x2yz;
    abcd[iGrid*1323+1035] = 2.0E0*I_ESP_I2x3yz_H2x2yz_a;
    abcd[iGrid*1323+1036] = 2.0E0*I_ESP_I2x2y2z_H2x2yz_a-1*I_ESP_G2x2y_H2x2yz;
    abcd[iGrid*1323+1037] = 2.0E0*I_ESP_I2xy3z_H2x2yz_a-2*I_ESP_G2xyz_H2x2yz;
    abcd[iGrid*1323+1038] = 2.0E0*I_ESP_I2x4z_H2x2yz_a-3*I_ESP_G2x2z_H2x2yz;
    abcd[iGrid*1323+1039] = 2.0E0*I_ESP_Ix4yz_H2x2yz_a;
    abcd[iGrid*1323+1040] = 2.0E0*I_ESP_Ix3y2z_H2x2yz_a-1*I_ESP_Gx3y_H2x2yz;
    abcd[iGrid*1323+1041] = 2.0E0*I_ESP_Ix2y3z_H2x2yz_a-2*I_ESP_Gx2yz_H2x2yz;
    abcd[iGrid*1323+1042] = 2.0E0*I_ESP_Ixy4z_H2x2yz_a-3*I_ESP_Gxy2z_H2x2yz;
    abcd[iGrid*1323+1043] = 2.0E0*I_ESP_Ix5z_H2x2yz_a-4*I_ESP_Gx3z_H2x2yz;
    abcd[iGrid*1323+1044] = 2.0E0*I_ESP_I5yz_H2x2yz_a;
    abcd[iGrid*1323+1045] = 2.0E0*I_ESP_I4y2z_H2x2yz_a-1*I_ESP_G4y_H2x2yz;
    abcd[iGrid*1323+1046] = 2.0E0*I_ESP_I3y3z_H2x2yz_a-2*I_ESP_G3yz_H2x2yz;
    abcd[iGrid*1323+1047] = 2.0E0*I_ESP_I2y4z_H2x2yz_a-3*I_ESP_G2y2z_H2x2yz;
    abcd[iGrid*1323+1048] = 2.0E0*I_ESP_Iy5z_H2x2yz_a-4*I_ESP_Gy3z_H2x2yz;
    abcd[iGrid*1323+1049] = 2.0E0*I_ESP_I6z_H2x2yz_a-5*I_ESP_G4z_H2x2yz;
    abcd[iGrid*1323+1050] = 2.0E0*I_ESP_I5xz_H2xy2z_a;
    abcd[iGrid*1323+1051] = 2.0E0*I_ESP_I4xyz_H2xy2z_a;
    abcd[iGrid*1323+1052] = 2.0E0*I_ESP_I4x2z_H2xy2z_a-1*I_ESP_G4x_H2xy2z;
    abcd[iGrid*1323+1053] = 2.0E0*I_ESP_I3x2yz_H2xy2z_a;
    abcd[iGrid*1323+1054] = 2.0E0*I_ESP_I3xy2z_H2xy2z_a-1*I_ESP_G3xy_H2xy2z;
    abcd[iGrid*1323+1055] = 2.0E0*I_ESP_I3x3z_H2xy2z_a-2*I_ESP_G3xz_H2xy2z;
    abcd[iGrid*1323+1056] = 2.0E0*I_ESP_I2x3yz_H2xy2z_a;
    abcd[iGrid*1323+1057] = 2.0E0*I_ESP_I2x2y2z_H2xy2z_a-1*I_ESP_G2x2y_H2xy2z;
    abcd[iGrid*1323+1058] = 2.0E0*I_ESP_I2xy3z_H2xy2z_a-2*I_ESP_G2xyz_H2xy2z;
    abcd[iGrid*1323+1059] = 2.0E0*I_ESP_I2x4z_H2xy2z_a-3*I_ESP_G2x2z_H2xy2z;
    abcd[iGrid*1323+1060] = 2.0E0*I_ESP_Ix4yz_H2xy2z_a;
    abcd[iGrid*1323+1061] = 2.0E0*I_ESP_Ix3y2z_H2xy2z_a-1*I_ESP_Gx3y_H2xy2z;
    abcd[iGrid*1323+1062] = 2.0E0*I_ESP_Ix2y3z_H2xy2z_a-2*I_ESP_Gx2yz_H2xy2z;
    abcd[iGrid*1323+1063] = 2.0E0*I_ESP_Ixy4z_H2xy2z_a-3*I_ESP_Gxy2z_H2xy2z;
    abcd[iGrid*1323+1064] = 2.0E0*I_ESP_Ix5z_H2xy2z_a-4*I_ESP_Gx3z_H2xy2z;
    abcd[iGrid*1323+1065] = 2.0E0*I_ESP_I5yz_H2xy2z_a;
    abcd[iGrid*1323+1066] = 2.0E0*I_ESP_I4y2z_H2xy2z_a-1*I_ESP_G4y_H2xy2z;
    abcd[iGrid*1323+1067] = 2.0E0*I_ESP_I3y3z_H2xy2z_a-2*I_ESP_G3yz_H2xy2z;
    abcd[iGrid*1323+1068] = 2.0E0*I_ESP_I2y4z_H2xy2z_a-3*I_ESP_G2y2z_H2xy2z;
    abcd[iGrid*1323+1069] = 2.0E0*I_ESP_Iy5z_H2xy2z_a-4*I_ESP_Gy3z_H2xy2z;
    abcd[iGrid*1323+1070] = 2.0E0*I_ESP_I6z_H2xy2z_a-5*I_ESP_G4z_H2xy2z;
    abcd[iGrid*1323+1071] = 2.0E0*I_ESP_I5xz_H2x3z_a;
    abcd[iGrid*1323+1072] = 2.0E0*I_ESP_I4xyz_H2x3z_a;
    abcd[iGrid*1323+1073] = 2.0E0*I_ESP_I4x2z_H2x3z_a-1*I_ESP_G4x_H2x3z;
    abcd[iGrid*1323+1074] = 2.0E0*I_ESP_I3x2yz_H2x3z_a;
    abcd[iGrid*1323+1075] = 2.0E0*I_ESP_I3xy2z_H2x3z_a-1*I_ESP_G3xy_H2x3z;
    abcd[iGrid*1323+1076] = 2.0E0*I_ESP_I3x3z_H2x3z_a-2*I_ESP_G3xz_H2x3z;
    abcd[iGrid*1323+1077] = 2.0E0*I_ESP_I2x3yz_H2x3z_a;
    abcd[iGrid*1323+1078] = 2.0E0*I_ESP_I2x2y2z_H2x3z_a-1*I_ESP_G2x2y_H2x3z;
    abcd[iGrid*1323+1079] = 2.0E0*I_ESP_I2xy3z_H2x3z_a-2*I_ESP_G2xyz_H2x3z;
    abcd[iGrid*1323+1080] = 2.0E0*I_ESP_I2x4z_H2x3z_a-3*I_ESP_G2x2z_H2x3z;
    abcd[iGrid*1323+1081] = 2.0E0*I_ESP_Ix4yz_H2x3z_a;
    abcd[iGrid*1323+1082] = 2.0E0*I_ESP_Ix3y2z_H2x3z_a-1*I_ESP_Gx3y_H2x3z;
    abcd[iGrid*1323+1083] = 2.0E0*I_ESP_Ix2y3z_H2x3z_a-2*I_ESP_Gx2yz_H2x3z;
    abcd[iGrid*1323+1084] = 2.0E0*I_ESP_Ixy4z_H2x3z_a-3*I_ESP_Gxy2z_H2x3z;
    abcd[iGrid*1323+1085] = 2.0E0*I_ESP_Ix5z_H2x3z_a-4*I_ESP_Gx3z_H2x3z;
    abcd[iGrid*1323+1086] = 2.0E0*I_ESP_I5yz_H2x3z_a;
    abcd[iGrid*1323+1087] = 2.0E0*I_ESP_I4y2z_H2x3z_a-1*I_ESP_G4y_H2x3z;
    abcd[iGrid*1323+1088] = 2.0E0*I_ESP_I3y3z_H2x3z_a-2*I_ESP_G3yz_H2x3z;
    abcd[iGrid*1323+1089] = 2.0E0*I_ESP_I2y4z_H2x3z_a-3*I_ESP_G2y2z_H2x3z;
    abcd[iGrid*1323+1090] = 2.0E0*I_ESP_Iy5z_H2x3z_a-4*I_ESP_Gy3z_H2x3z;
    abcd[iGrid*1323+1091] = 2.0E0*I_ESP_I6z_H2x3z_a-5*I_ESP_G4z_H2x3z;
    abcd[iGrid*1323+1092] = 2.0E0*I_ESP_I5xz_Hx4y_a;
    abcd[iGrid*1323+1093] = 2.0E0*I_ESP_I4xyz_Hx4y_a;
    abcd[iGrid*1323+1094] = 2.0E0*I_ESP_I4x2z_Hx4y_a-1*I_ESP_G4x_Hx4y;
    abcd[iGrid*1323+1095] = 2.0E0*I_ESP_I3x2yz_Hx4y_a;
    abcd[iGrid*1323+1096] = 2.0E0*I_ESP_I3xy2z_Hx4y_a-1*I_ESP_G3xy_Hx4y;
    abcd[iGrid*1323+1097] = 2.0E0*I_ESP_I3x3z_Hx4y_a-2*I_ESP_G3xz_Hx4y;
    abcd[iGrid*1323+1098] = 2.0E0*I_ESP_I2x3yz_Hx4y_a;
    abcd[iGrid*1323+1099] = 2.0E0*I_ESP_I2x2y2z_Hx4y_a-1*I_ESP_G2x2y_Hx4y;
    abcd[iGrid*1323+1100] = 2.0E0*I_ESP_I2xy3z_Hx4y_a-2*I_ESP_G2xyz_Hx4y;
    abcd[iGrid*1323+1101] = 2.0E0*I_ESP_I2x4z_Hx4y_a-3*I_ESP_G2x2z_Hx4y;
    abcd[iGrid*1323+1102] = 2.0E0*I_ESP_Ix4yz_Hx4y_a;
    abcd[iGrid*1323+1103] = 2.0E0*I_ESP_Ix3y2z_Hx4y_a-1*I_ESP_Gx3y_Hx4y;
    abcd[iGrid*1323+1104] = 2.0E0*I_ESP_Ix2y3z_Hx4y_a-2*I_ESP_Gx2yz_Hx4y;
    abcd[iGrid*1323+1105] = 2.0E0*I_ESP_Ixy4z_Hx4y_a-3*I_ESP_Gxy2z_Hx4y;
    abcd[iGrid*1323+1106] = 2.0E0*I_ESP_Ix5z_Hx4y_a-4*I_ESP_Gx3z_Hx4y;
    abcd[iGrid*1323+1107] = 2.0E0*I_ESP_I5yz_Hx4y_a;
    abcd[iGrid*1323+1108] = 2.0E0*I_ESP_I4y2z_Hx4y_a-1*I_ESP_G4y_Hx4y;
    abcd[iGrid*1323+1109] = 2.0E0*I_ESP_I3y3z_Hx4y_a-2*I_ESP_G3yz_Hx4y;
    abcd[iGrid*1323+1110] = 2.0E0*I_ESP_I2y4z_Hx4y_a-3*I_ESP_G2y2z_Hx4y;
    abcd[iGrid*1323+1111] = 2.0E0*I_ESP_Iy5z_Hx4y_a-4*I_ESP_Gy3z_Hx4y;
    abcd[iGrid*1323+1112] = 2.0E0*I_ESP_I6z_Hx4y_a-5*I_ESP_G4z_Hx4y;
    abcd[iGrid*1323+1113] = 2.0E0*I_ESP_I5xz_Hx3yz_a;
    abcd[iGrid*1323+1114] = 2.0E0*I_ESP_I4xyz_Hx3yz_a;
    abcd[iGrid*1323+1115] = 2.0E0*I_ESP_I4x2z_Hx3yz_a-1*I_ESP_G4x_Hx3yz;
    abcd[iGrid*1323+1116] = 2.0E0*I_ESP_I3x2yz_Hx3yz_a;
    abcd[iGrid*1323+1117] = 2.0E0*I_ESP_I3xy2z_Hx3yz_a-1*I_ESP_G3xy_Hx3yz;
    abcd[iGrid*1323+1118] = 2.0E0*I_ESP_I3x3z_Hx3yz_a-2*I_ESP_G3xz_Hx3yz;
    abcd[iGrid*1323+1119] = 2.0E0*I_ESP_I2x3yz_Hx3yz_a;
    abcd[iGrid*1323+1120] = 2.0E0*I_ESP_I2x2y2z_Hx3yz_a-1*I_ESP_G2x2y_Hx3yz;
    abcd[iGrid*1323+1121] = 2.0E0*I_ESP_I2xy3z_Hx3yz_a-2*I_ESP_G2xyz_Hx3yz;
    abcd[iGrid*1323+1122] = 2.0E0*I_ESP_I2x4z_Hx3yz_a-3*I_ESP_G2x2z_Hx3yz;
    abcd[iGrid*1323+1123] = 2.0E0*I_ESP_Ix4yz_Hx3yz_a;
    abcd[iGrid*1323+1124] = 2.0E0*I_ESP_Ix3y2z_Hx3yz_a-1*I_ESP_Gx3y_Hx3yz;
    abcd[iGrid*1323+1125] = 2.0E0*I_ESP_Ix2y3z_Hx3yz_a-2*I_ESP_Gx2yz_Hx3yz;
    abcd[iGrid*1323+1126] = 2.0E0*I_ESP_Ixy4z_Hx3yz_a-3*I_ESP_Gxy2z_Hx3yz;
    abcd[iGrid*1323+1127] = 2.0E0*I_ESP_Ix5z_Hx3yz_a-4*I_ESP_Gx3z_Hx3yz;
    abcd[iGrid*1323+1128] = 2.0E0*I_ESP_I5yz_Hx3yz_a;
    abcd[iGrid*1323+1129] = 2.0E0*I_ESP_I4y2z_Hx3yz_a-1*I_ESP_G4y_Hx3yz;
    abcd[iGrid*1323+1130] = 2.0E0*I_ESP_I3y3z_Hx3yz_a-2*I_ESP_G3yz_Hx3yz;
    abcd[iGrid*1323+1131] = 2.0E0*I_ESP_I2y4z_Hx3yz_a-3*I_ESP_G2y2z_Hx3yz;
    abcd[iGrid*1323+1132] = 2.0E0*I_ESP_Iy5z_Hx3yz_a-4*I_ESP_Gy3z_Hx3yz;
    abcd[iGrid*1323+1133] = 2.0E0*I_ESP_I6z_Hx3yz_a-5*I_ESP_G4z_Hx3yz;
    abcd[iGrid*1323+1134] = 2.0E0*I_ESP_I5xz_Hx2y2z_a;
    abcd[iGrid*1323+1135] = 2.0E0*I_ESP_I4xyz_Hx2y2z_a;
    abcd[iGrid*1323+1136] = 2.0E0*I_ESP_I4x2z_Hx2y2z_a-1*I_ESP_G4x_Hx2y2z;
    abcd[iGrid*1323+1137] = 2.0E0*I_ESP_I3x2yz_Hx2y2z_a;
    abcd[iGrid*1323+1138] = 2.0E0*I_ESP_I3xy2z_Hx2y2z_a-1*I_ESP_G3xy_Hx2y2z;
    abcd[iGrid*1323+1139] = 2.0E0*I_ESP_I3x3z_Hx2y2z_a-2*I_ESP_G3xz_Hx2y2z;
    abcd[iGrid*1323+1140] = 2.0E0*I_ESP_I2x3yz_Hx2y2z_a;
    abcd[iGrid*1323+1141] = 2.0E0*I_ESP_I2x2y2z_Hx2y2z_a-1*I_ESP_G2x2y_Hx2y2z;
    abcd[iGrid*1323+1142] = 2.0E0*I_ESP_I2xy3z_Hx2y2z_a-2*I_ESP_G2xyz_Hx2y2z;
    abcd[iGrid*1323+1143] = 2.0E0*I_ESP_I2x4z_Hx2y2z_a-3*I_ESP_G2x2z_Hx2y2z;
    abcd[iGrid*1323+1144] = 2.0E0*I_ESP_Ix4yz_Hx2y2z_a;
    abcd[iGrid*1323+1145] = 2.0E0*I_ESP_Ix3y2z_Hx2y2z_a-1*I_ESP_Gx3y_Hx2y2z;
    abcd[iGrid*1323+1146] = 2.0E0*I_ESP_Ix2y3z_Hx2y2z_a-2*I_ESP_Gx2yz_Hx2y2z;
    abcd[iGrid*1323+1147] = 2.0E0*I_ESP_Ixy4z_Hx2y2z_a-3*I_ESP_Gxy2z_Hx2y2z;
    abcd[iGrid*1323+1148] = 2.0E0*I_ESP_Ix5z_Hx2y2z_a-4*I_ESP_Gx3z_Hx2y2z;
    abcd[iGrid*1323+1149] = 2.0E0*I_ESP_I5yz_Hx2y2z_a;
    abcd[iGrid*1323+1150] = 2.0E0*I_ESP_I4y2z_Hx2y2z_a-1*I_ESP_G4y_Hx2y2z;
    abcd[iGrid*1323+1151] = 2.0E0*I_ESP_I3y3z_Hx2y2z_a-2*I_ESP_G3yz_Hx2y2z;
    abcd[iGrid*1323+1152] = 2.0E0*I_ESP_I2y4z_Hx2y2z_a-3*I_ESP_G2y2z_Hx2y2z;
    abcd[iGrid*1323+1153] = 2.0E0*I_ESP_Iy5z_Hx2y2z_a-4*I_ESP_Gy3z_Hx2y2z;
    abcd[iGrid*1323+1154] = 2.0E0*I_ESP_I6z_Hx2y2z_a-5*I_ESP_G4z_Hx2y2z;
    abcd[iGrid*1323+1155] = 2.0E0*I_ESP_I5xz_Hxy3z_a;
    abcd[iGrid*1323+1156] = 2.0E0*I_ESP_I4xyz_Hxy3z_a;
    abcd[iGrid*1323+1157] = 2.0E0*I_ESP_I4x2z_Hxy3z_a-1*I_ESP_G4x_Hxy3z;
    abcd[iGrid*1323+1158] = 2.0E0*I_ESP_I3x2yz_Hxy3z_a;
    abcd[iGrid*1323+1159] = 2.0E0*I_ESP_I3xy2z_Hxy3z_a-1*I_ESP_G3xy_Hxy3z;
    abcd[iGrid*1323+1160] = 2.0E0*I_ESP_I3x3z_Hxy3z_a-2*I_ESP_G3xz_Hxy3z;
    abcd[iGrid*1323+1161] = 2.0E0*I_ESP_I2x3yz_Hxy3z_a;
    abcd[iGrid*1323+1162] = 2.0E0*I_ESP_I2x2y2z_Hxy3z_a-1*I_ESP_G2x2y_Hxy3z;
    abcd[iGrid*1323+1163] = 2.0E0*I_ESP_I2xy3z_Hxy3z_a-2*I_ESP_G2xyz_Hxy3z;
    abcd[iGrid*1323+1164] = 2.0E0*I_ESP_I2x4z_Hxy3z_a-3*I_ESP_G2x2z_Hxy3z;
    abcd[iGrid*1323+1165] = 2.0E0*I_ESP_Ix4yz_Hxy3z_a;
    abcd[iGrid*1323+1166] = 2.0E0*I_ESP_Ix3y2z_Hxy3z_a-1*I_ESP_Gx3y_Hxy3z;
    abcd[iGrid*1323+1167] = 2.0E0*I_ESP_Ix2y3z_Hxy3z_a-2*I_ESP_Gx2yz_Hxy3z;
    abcd[iGrid*1323+1168] = 2.0E0*I_ESP_Ixy4z_Hxy3z_a-3*I_ESP_Gxy2z_Hxy3z;
    abcd[iGrid*1323+1169] = 2.0E0*I_ESP_Ix5z_Hxy3z_a-4*I_ESP_Gx3z_Hxy3z;
    abcd[iGrid*1323+1170] = 2.0E0*I_ESP_I5yz_Hxy3z_a;
    abcd[iGrid*1323+1171] = 2.0E0*I_ESP_I4y2z_Hxy3z_a-1*I_ESP_G4y_Hxy3z;
    abcd[iGrid*1323+1172] = 2.0E0*I_ESP_I3y3z_Hxy3z_a-2*I_ESP_G3yz_Hxy3z;
    abcd[iGrid*1323+1173] = 2.0E0*I_ESP_I2y4z_Hxy3z_a-3*I_ESP_G2y2z_Hxy3z;
    abcd[iGrid*1323+1174] = 2.0E0*I_ESP_Iy5z_Hxy3z_a-4*I_ESP_Gy3z_Hxy3z;
    abcd[iGrid*1323+1175] = 2.0E0*I_ESP_I6z_Hxy3z_a-5*I_ESP_G4z_Hxy3z;
    abcd[iGrid*1323+1176] = 2.0E0*I_ESP_I5xz_Hx4z_a;
    abcd[iGrid*1323+1177] = 2.0E0*I_ESP_I4xyz_Hx4z_a;
    abcd[iGrid*1323+1178] = 2.0E0*I_ESP_I4x2z_Hx4z_a-1*I_ESP_G4x_Hx4z;
    abcd[iGrid*1323+1179] = 2.0E0*I_ESP_I3x2yz_Hx4z_a;
    abcd[iGrid*1323+1180] = 2.0E0*I_ESP_I3xy2z_Hx4z_a-1*I_ESP_G3xy_Hx4z;
    abcd[iGrid*1323+1181] = 2.0E0*I_ESP_I3x3z_Hx4z_a-2*I_ESP_G3xz_Hx4z;
    abcd[iGrid*1323+1182] = 2.0E0*I_ESP_I2x3yz_Hx4z_a;
    abcd[iGrid*1323+1183] = 2.0E0*I_ESP_I2x2y2z_Hx4z_a-1*I_ESP_G2x2y_Hx4z;
    abcd[iGrid*1323+1184] = 2.0E0*I_ESP_I2xy3z_Hx4z_a-2*I_ESP_G2xyz_Hx4z;
    abcd[iGrid*1323+1185] = 2.0E0*I_ESP_I2x4z_Hx4z_a-3*I_ESP_G2x2z_Hx4z;
    abcd[iGrid*1323+1186] = 2.0E0*I_ESP_Ix4yz_Hx4z_a;
    abcd[iGrid*1323+1187] = 2.0E0*I_ESP_Ix3y2z_Hx4z_a-1*I_ESP_Gx3y_Hx4z;
    abcd[iGrid*1323+1188] = 2.0E0*I_ESP_Ix2y3z_Hx4z_a-2*I_ESP_Gx2yz_Hx4z;
    abcd[iGrid*1323+1189] = 2.0E0*I_ESP_Ixy4z_Hx4z_a-3*I_ESP_Gxy2z_Hx4z;
    abcd[iGrid*1323+1190] = 2.0E0*I_ESP_Ix5z_Hx4z_a-4*I_ESP_Gx3z_Hx4z;
    abcd[iGrid*1323+1191] = 2.0E0*I_ESP_I5yz_Hx4z_a;
    abcd[iGrid*1323+1192] = 2.0E0*I_ESP_I4y2z_Hx4z_a-1*I_ESP_G4y_Hx4z;
    abcd[iGrid*1323+1193] = 2.0E0*I_ESP_I3y3z_Hx4z_a-2*I_ESP_G3yz_Hx4z;
    abcd[iGrid*1323+1194] = 2.0E0*I_ESP_I2y4z_Hx4z_a-3*I_ESP_G2y2z_Hx4z;
    abcd[iGrid*1323+1195] = 2.0E0*I_ESP_Iy5z_Hx4z_a-4*I_ESP_Gy3z_Hx4z;
    abcd[iGrid*1323+1196] = 2.0E0*I_ESP_I6z_Hx4z_a-5*I_ESP_G4z_Hx4z;
    abcd[iGrid*1323+1197] = 2.0E0*I_ESP_I5xz_H5y_a;
    abcd[iGrid*1323+1198] = 2.0E0*I_ESP_I4xyz_H5y_a;
    abcd[iGrid*1323+1199] = 2.0E0*I_ESP_I4x2z_H5y_a-1*I_ESP_G4x_H5y;
    abcd[iGrid*1323+1200] = 2.0E0*I_ESP_I3x2yz_H5y_a;
    abcd[iGrid*1323+1201] = 2.0E0*I_ESP_I3xy2z_H5y_a-1*I_ESP_G3xy_H5y;
    abcd[iGrid*1323+1202] = 2.0E0*I_ESP_I3x3z_H5y_a-2*I_ESP_G3xz_H5y;
    abcd[iGrid*1323+1203] = 2.0E0*I_ESP_I2x3yz_H5y_a;
    abcd[iGrid*1323+1204] = 2.0E0*I_ESP_I2x2y2z_H5y_a-1*I_ESP_G2x2y_H5y;
    abcd[iGrid*1323+1205] = 2.0E0*I_ESP_I2xy3z_H5y_a-2*I_ESP_G2xyz_H5y;
    abcd[iGrid*1323+1206] = 2.0E0*I_ESP_I2x4z_H5y_a-3*I_ESP_G2x2z_H5y;
    abcd[iGrid*1323+1207] = 2.0E0*I_ESP_Ix4yz_H5y_a;
    abcd[iGrid*1323+1208] = 2.0E0*I_ESP_Ix3y2z_H5y_a-1*I_ESP_Gx3y_H5y;
    abcd[iGrid*1323+1209] = 2.0E0*I_ESP_Ix2y3z_H5y_a-2*I_ESP_Gx2yz_H5y;
    abcd[iGrid*1323+1210] = 2.0E0*I_ESP_Ixy4z_H5y_a-3*I_ESP_Gxy2z_H5y;
    abcd[iGrid*1323+1211] = 2.0E0*I_ESP_Ix5z_H5y_a-4*I_ESP_Gx3z_H5y;
    abcd[iGrid*1323+1212] = 2.0E0*I_ESP_I5yz_H5y_a;
    abcd[iGrid*1323+1213] = 2.0E0*I_ESP_I4y2z_H5y_a-1*I_ESP_G4y_H5y;
    abcd[iGrid*1323+1214] = 2.0E0*I_ESP_I3y3z_H5y_a-2*I_ESP_G3yz_H5y;
    abcd[iGrid*1323+1215] = 2.0E0*I_ESP_I2y4z_H5y_a-3*I_ESP_G2y2z_H5y;
    abcd[iGrid*1323+1216] = 2.0E0*I_ESP_Iy5z_H5y_a-4*I_ESP_Gy3z_H5y;
    abcd[iGrid*1323+1217] = 2.0E0*I_ESP_I6z_H5y_a-5*I_ESP_G4z_H5y;
    abcd[iGrid*1323+1218] = 2.0E0*I_ESP_I5xz_H4yz_a;
    abcd[iGrid*1323+1219] = 2.0E0*I_ESP_I4xyz_H4yz_a;
    abcd[iGrid*1323+1220] = 2.0E0*I_ESP_I4x2z_H4yz_a-1*I_ESP_G4x_H4yz;
    abcd[iGrid*1323+1221] = 2.0E0*I_ESP_I3x2yz_H4yz_a;
    abcd[iGrid*1323+1222] = 2.0E0*I_ESP_I3xy2z_H4yz_a-1*I_ESP_G3xy_H4yz;
    abcd[iGrid*1323+1223] = 2.0E0*I_ESP_I3x3z_H4yz_a-2*I_ESP_G3xz_H4yz;
    abcd[iGrid*1323+1224] = 2.0E0*I_ESP_I2x3yz_H4yz_a;
    abcd[iGrid*1323+1225] = 2.0E0*I_ESP_I2x2y2z_H4yz_a-1*I_ESP_G2x2y_H4yz;
    abcd[iGrid*1323+1226] = 2.0E0*I_ESP_I2xy3z_H4yz_a-2*I_ESP_G2xyz_H4yz;
    abcd[iGrid*1323+1227] = 2.0E0*I_ESP_I2x4z_H4yz_a-3*I_ESP_G2x2z_H4yz;
    abcd[iGrid*1323+1228] = 2.0E0*I_ESP_Ix4yz_H4yz_a;
    abcd[iGrid*1323+1229] = 2.0E0*I_ESP_Ix3y2z_H4yz_a-1*I_ESP_Gx3y_H4yz;
    abcd[iGrid*1323+1230] = 2.0E0*I_ESP_Ix2y3z_H4yz_a-2*I_ESP_Gx2yz_H4yz;
    abcd[iGrid*1323+1231] = 2.0E0*I_ESP_Ixy4z_H4yz_a-3*I_ESP_Gxy2z_H4yz;
    abcd[iGrid*1323+1232] = 2.0E0*I_ESP_Ix5z_H4yz_a-4*I_ESP_Gx3z_H4yz;
    abcd[iGrid*1323+1233] = 2.0E0*I_ESP_I5yz_H4yz_a;
    abcd[iGrid*1323+1234] = 2.0E0*I_ESP_I4y2z_H4yz_a-1*I_ESP_G4y_H4yz;
    abcd[iGrid*1323+1235] = 2.0E0*I_ESP_I3y3z_H4yz_a-2*I_ESP_G3yz_H4yz;
    abcd[iGrid*1323+1236] = 2.0E0*I_ESP_I2y4z_H4yz_a-3*I_ESP_G2y2z_H4yz;
    abcd[iGrid*1323+1237] = 2.0E0*I_ESP_Iy5z_H4yz_a-4*I_ESP_Gy3z_H4yz;
    abcd[iGrid*1323+1238] = 2.0E0*I_ESP_I6z_H4yz_a-5*I_ESP_G4z_H4yz;
    abcd[iGrid*1323+1239] = 2.0E0*I_ESP_I5xz_H3y2z_a;
    abcd[iGrid*1323+1240] = 2.0E0*I_ESP_I4xyz_H3y2z_a;
    abcd[iGrid*1323+1241] = 2.0E0*I_ESP_I4x2z_H3y2z_a-1*I_ESP_G4x_H3y2z;
    abcd[iGrid*1323+1242] = 2.0E0*I_ESP_I3x2yz_H3y2z_a;
    abcd[iGrid*1323+1243] = 2.0E0*I_ESP_I3xy2z_H3y2z_a-1*I_ESP_G3xy_H3y2z;
    abcd[iGrid*1323+1244] = 2.0E0*I_ESP_I3x3z_H3y2z_a-2*I_ESP_G3xz_H3y2z;
    abcd[iGrid*1323+1245] = 2.0E0*I_ESP_I2x3yz_H3y2z_a;
    abcd[iGrid*1323+1246] = 2.0E0*I_ESP_I2x2y2z_H3y2z_a-1*I_ESP_G2x2y_H3y2z;
    abcd[iGrid*1323+1247] = 2.0E0*I_ESP_I2xy3z_H3y2z_a-2*I_ESP_G2xyz_H3y2z;
    abcd[iGrid*1323+1248] = 2.0E0*I_ESP_I2x4z_H3y2z_a-3*I_ESP_G2x2z_H3y2z;
    abcd[iGrid*1323+1249] = 2.0E0*I_ESP_Ix4yz_H3y2z_a;
    abcd[iGrid*1323+1250] = 2.0E0*I_ESP_Ix3y2z_H3y2z_a-1*I_ESP_Gx3y_H3y2z;
    abcd[iGrid*1323+1251] = 2.0E0*I_ESP_Ix2y3z_H3y2z_a-2*I_ESP_Gx2yz_H3y2z;
    abcd[iGrid*1323+1252] = 2.0E0*I_ESP_Ixy4z_H3y2z_a-3*I_ESP_Gxy2z_H3y2z;
    abcd[iGrid*1323+1253] = 2.0E0*I_ESP_Ix5z_H3y2z_a-4*I_ESP_Gx3z_H3y2z;
    abcd[iGrid*1323+1254] = 2.0E0*I_ESP_I5yz_H3y2z_a;
    abcd[iGrid*1323+1255] = 2.0E0*I_ESP_I4y2z_H3y2z_a-1*I_ESP_G4y_H3y2z;
    abcd[iGrid*1323+1256] = 2.0E0*I_ESP_I3y3z_H3y2z_a-2*I_ESP_G3yz_H3y2z;
    abcd[iGrid*1323+1257] = 2.0E0*I_ESP_I2y4z_H3y2z_a-3*I_ESP_G2y2z_H3y2z;
    abcd[iGrid*1323+1258] = 2.0E0*I_ESP_Iy5z_H3y2z_a-4*I_ESP_Gy3z_H3y2z;
    abcd[iGrid*1323+1259] = 2.0E0*I_ESP_I6z_H3y2z_a-5*I_ESP_G4z_H3y2z;
    abcd[iGrid*1323+1260] = 2.0E0*I_ESP_I5xz_H2y3z_a;
    abcd[iGrid*1323+1261] = 2.0E0*I_ESP_I4xyz_H2y3z_a;
    abcd[iGrid*1323+1262] = 2.0E0*I_ESP_I4x2z_H2y3z_a-1*I_ESP_G4x_H2y3z;
    abcd[iGrid*1323+1263] = 2.0E0*I_ESP_I3x2yz_H2y3z_a;
    abcd[iGrid*1323+1264] = 2.0E0*I_ESP_I3xy2z_H2y3z_a-1*I_ESP_G3xy_H2y3z;
    abcd[iGrid*1323+1265] = 2.0E0*I_ESP_I3x3z_H2y3z_a-2*I_ESP_G3xz_H2y3z;
    abcd[iGrid*1323+1266] = 2.0E0*I_ESP_I2x3yz_H2y3z_a;
    abcd[iGrid*1323+1267] = 2.0E0*I_ESP_I2x2y2z_H2y3z_a-1*I_ESP_G2x2y_H2y3z;
    abcd[iGrid*1323+1268] = 2.0E0*I_ESP_I2xy3z_H2y3z_a-2*I_ESP_G2xyz_H2y3z;
    abcd[iGrid*1323+1269] = 2.0E0*I_ESP_I2x4z_H2y3z_a-3*I_ESP_G2x2z_H2y3z;
    abcd[iGrid*1323+1270] = 2.0E0*I_ESP_Ix4yz_H2y3z_a;
    abcd[iGrid*1323+1271] = 2.0E0*I_ESP_Ix3y2z_H2y3z_a-1*I_ESP_Gx3y_H2y3z;
    abcd[iGrid*1323+1272] = 2.0E0*I_ESP_Ix2y3z_H2y3z_a-2*I_ESP_Gx2yz_H2y3z;
    abcd[iGrid*1323+1273] = 2.0E0*I_ESP_Ixy4z_H2y3z_a-3*I_ESP_Gxy2z_H2y3z;
    abcd[iGrid*1323+1274] = 2.0E0*I_ESP_Ix5z_H2y3z_a-4*I_ESP_Gx3z_H2y3z;
    abcd[iGrid*1323+1275] = 2.0E0*I_ESP_I5yz_H2y3z_a;
    abcd[iGrid*1323+1276] = 2.0E0*I_ESP_I4y2z_H2y3z_a-1*I_ESP_G4y_H2y3z;
    abcd[iGrid*1323+1277] = 2.0E0*I_ESP_I3y3z_H2y3z_a-2*I_ESP_G3yz_H2y3z;
    abcd[iGrid*1323+1278] = 2.0E0*I_ESP_I2y4z_H2y3z_a-3*I_ESP_G2y2z_H2y3z;
    abcd[iGrid*1323+1279] = 2.0E0*I_ESP_Iy5z_H2y3z_a-4*I_ESP_Gy3z_H2y3z;
    abcd[iGrid*1323+1280] = 2.0E0*I_ESP_I6z_H2y3z_a-5*I_ESP_G4z_H2y3z;
    abcd[iGrid*1323+1281] = 2.0E0*I_ESP_I5xz_Hy4z_a;
    abcd[iGrid*1323+1282] = 2.0E0*I_ESP_I4xyz_Hy4z_a;
    abcd[iGrid*1323+1283] = 2.0E0*I_ESP_I4x2z_Hy4z_a-1*I_ESP_G4x_Hy4z;
    abcd[iGrid*1323+1284] = 2.0E0*I_ESP_I3x2yz_Hy4z_a;
    abcd[iGrid*1323+1285] = 2.0E0*I_ESP_I3xy2z_Hy4z_a-1*I_ESP_G3xy_Hy4z;
    abcd[iGrid*1323+1286] = 2.0E0*I_ESP_I3x3z_Hy4z_a-2*I_ESP_G3xz_Hy4z;
    abcd[iGrid*1323+1287] = 2.0E0*I_ESP_I2x3yz_Hy4z_a;
    abcd[iGrid*1323+1288] = 2.0E0*I_ESP_I2x2y2z_Hy4z_a-1*I_ESP_G2x2y_Hy4z;
    abcd[iGrid*1323+1289] = 2.0E0*I_ESP_I2xy3z_Hy4z_a-2*I_ESP_G2xyz_Hy4z;
    abcd[iGrid*1323+1290] = 2.0E0*I_ESP_I2x4z_Hy4z_a-3*I_ESP_G2x2z_Hy4z;
    abcd[iGrid*1323+1291] = 2.0E0*I_ESP_Ix4yz_Hy4z_a;
    abcd[iGrid*1323+1292] = 2.0E0*I_ESP_Ix3y2z_Hy4z_a-1*I_ESP_Gx3y_Hy4z;
    abcd[iGrid*1323+1293] = 2.0E0*I_ESP_Ix2y3z_Hy4z_a-2*I_ESP_Gx2yz_Hy4z;
    abcd[iGrid*1323+1294] = 2.0E0*I_ESP_Ixy4z_Hy4z_a-3*I_ESP_Gxy2z_Hy4z;
    abcd[iGrid*1323+1295] = 2.0E0*I_ESP_Ix5z_Hy4z_a-4*I_ESP_Gx3z_Hy4z;
    abcd[iGrid*1323+1296] = 2.0E0*I_ESP_I5yz_Hy4z_a;
    abcd[iGrid*1323+1297] = 2.0E0*I_ESP_I4y2z_Hy4z_a-1*I_ESP_G4y_Hy4z;
    abcd[iGrid*1323+1298] = 2.0E0*I_ESP_I3y3z_Hy4z_a-2*I_ESP_G3yz_Hy4z;
    abcd[iGrid*1323+1299] = 2.0E0*I_ESP_I2y4z_Hy4z_a-3*I_ESP_G2y2z_Hy4z;
    abcd[iGrid*1323+1300] = 2.0E0*I_ESP_Iy5z_Hy4z_a-4*I_ESP_Gy3z_Hy4z;
    abcd[iGrid*1323+1301] = 2.0E0*I_ESP_I6z_Hy4z_a-5*I_ESP_G4z_Hy4z;
    abcd[iGrid*1323+1302] = 2.0E0*I_ESP_I5xz_H5z_a;
    abcd[iGrid*1323+1303] = 2.0E0*I_ESP_I4xyz_H5z_a;
    abcd[iGrid*1323+1304] = 2.0E0*I_ESP_I4x2z_H5z_a-1*I_ESP_G4x_H5z;
    abcd[iGrid*1323+1305] = 2.0E0*I_ESP_I3x2yz_H5z_a;
    abcd[iGrid*1323+1306] = 2.0E0*I_ESP_I3xy2z_H5z_a-1*I_ESP_G3xy_H5z;
    abcd[iGrid*1323+1307] = 2.0E0*I_ESP_I3x3z_H5z_a-2*I_ESP_G3xz_H5z;
    abcd[iGrid*1323+1308] = 2.0E0*I_ESP_I2x3yz_H5z_a;
    abcd[iGrid*1323+1309] = 2.0E0*I_ESP_I2x2y2z_H5z_a-1*I_ESP_G2x2y_H5z;
    abcd[iGrid*1323+1310] = 2.0E0*I_ESP_I2xy3z_H5z_a-2*I_ESP_G2xyz_H5z;
    abcd[iGrid*1323+1311] = 2.0E0*I_ESP_I2x4z_H5z_a-3*I_ESP_G2x2z_H5z;
    abcd[iGrid*1323+1312] = 2.0E0*I_ESP_Ix4yz_H5z_a;
    abcd[iGrid*1323+1313] = 2.0E0*I_ESP_Ix3y2z_H5z_a-1*I_ESP_Gx3y_H5z;
    abcd[iGrid*1323+1314] = 2.0E0*I_ESP_Ix2y3z_H5z_a-2*I_ESP_Gx2yz_H5z;
    abcd[iGrid*1323+1315] = 2.0E0*I_ESP_Ixy4z_H5z_a-3*I_ESP_Gxy2z_H5z;
    abcd[iGrid*1323+1316] = 2.0E0*I_ESP_Ix5z_H5z_a-4*I_ESP_Gx3z_H5z;
    abcd[iGrid*1323+1317] = 2.0E0*I_ESP_I5yz_H5z_a;
    abcd[iGrid*1323+1318] = 2.0E0*I_ESP_I4y2z_H5z_a-1*I_ESP_G4y_H5z;
    abcd[iGrid*1323+1319] = 2.0E0*I_ESP_I3y3z_H5z_a-2*I_ESP_G3yz_H5z;
    abcd[iGrid*1323+1320] = 2.0E0*I_ESP_I2y4z_H5z_a-3*I_ESP_G2y2z_H5z;
    abcd[iGrid*1323+1321] = 2.0E0*I_ESP_Iy5z_H5z_a-4*I_ESP_Gy3z_H5z;
    abcd[iGrid*1323+1322] = 2.0E0*I_ESP_I6z_H5z_a-5*I_ESP_G4z_H5z;
  }
}
