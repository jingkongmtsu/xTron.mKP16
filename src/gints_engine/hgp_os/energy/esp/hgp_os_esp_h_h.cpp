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

void hgp_os_esp_h_h(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_N10x_S = 0.0E0;
    Double I_ESP_N9xy_S = 0.0E0;
    Double I_ESP_N9xz_S = 0.0E0;
    Double I_ESP_N8x2y_S = 0.0E0;
    Double I_ESP_N8xyz_S = 0.0E0;
    Double I_ESP_N8x2z_S = 0.0E0;
    Double I_ESP_N7x3y_S = 0.0E0;
    Double I_ESP_N7x2yz_S = 0.0E0;
    Double I_ESP_N7xy2z_S = 0.0E0;
    Double I_ESP_N7x3z_S = 0.0E0;
    Double I_ESP_N6x4y_S = 0.0E0;
    Double I_ESP_N6x3yz_S = 0.0E0;
    Double I_ESP_N6x2y2z_S = 0.0E0;
    Double I_ESP_N6xy3z_S = 0.0E0;
    Double I_ESP_N6x4z_S = 0.0E0;
    Double I_ESP_N5x5y_S = 0.0E0;
    Double I_ESP_N5x4yz_S = 0.0E0;
    Double I_ESP_N5x3y2z_S = 0.0E0;
    Double I_ESP_N5x2y3z_S = 0.0E0;
    Double I_ESP_N5xy4z_S = 0.0E0;
    Double I_ESP_N5x5z_S = 0.0E0;
    Double I_ESP_N4x6y_S = 0.0E0;
    Double I_ESP_N4x5yz_S = 0.0E0;
    Double I_ESP_N4x4y2z_S = 0.0E0;
    Double I_ESP_N4x3y3z_S = 0.0E0;
    Double I_ESP_N4x2y4z_S = 0.0E0;
    Double I_ESP_N4xy5z_S = 0.0E0;
    Double I_ESP_N4x6z_S = 0.0E0;
    Double I_ESP_N3x7y_S = 0.0E0;
    Double I_ESP_N3x6yz_S = 0.0E0;
    Double I_ESP_N3x5y2z_S = 0.0E0;
    Double I_ESP_N3x4y3z_S = 0.0E0;
    Double I_ESP_N3x3y4z_S = 0.0E0;
    Double I_ESP_N3x2y5z_S = 0.0E0;
    Double I_ESP_N3xy6z_S = 0.0E0;
    Double I_ESP_N3x7z_S = 0.0E0;
    Double I_ESP_N2x8y_S = 0.0E0;
    Double I_ESP_N2x7yz_S = 0.0E0;
    Double I_ESP_N2x6y2z_S = 0.0E0;
    Double I_ESP_N2x5y3z_S = 0.0E0;
    Double I_ESP_N2x4y4z_S = 0.0E0;
    Double I_ESP_N2x3y5z_S = 0.0E0;
    Double I_ESP_N2x2y6z_S = 0.0E0;
    Double I_ESP_N2xy7z_S = 0.0E0;
    Double I_ESP_N2x8z_S = 0.0E0;
    Double I_ESP_Nx9y_S = 0.0E0;
    Double I_ESP_Nx8yz_S = 0.0E0;
    Double I_ESP_Nx7y2z_S = 0.0E0;
    Double I_ESP_Nx6y3z_S = 0.0E0;
    Double I_ESP_Nx5y4z_S = 0.0E0;
    Double I_ESP_Nx4y5z_S = 0.0E0;
    Double I_ESP_Nx3y6z_S = 0.0E0;
    Double I_ESP_Nx2y7z_S = 0.0E0;
    Double I_ESP_Nxy8z_S = 0.0E0;
    Double I_ESP_Nx9z_S = 0.0E0;
    Double I_ESP_N10y_S = 0.0E0;
    Double I_ESP_N9yz_S = 0.0E0;
    Double I_ESP_N8y2z_S = 0.0E0;
    Double I_ESP_N7y3z_S = 0.0E0;
    Double I_ESP_N6y4z_S = 0.0E0;
    Double I_ESP_N5y5z_S = 0.0E0;
    Double I_ESP_N4y6z_S = 0.0E0;
    Double I_ESP_N3y7z_S = 0.0E0;
    Double I_ESP_N2y8z_S = 0.0E0;
    Double I_ESP_Ny9z_S = 0.0E0;
    Double I_ESP_N10z_S = 0.0E0;
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

    for(UInt ip2=0; ip2<inp2; ip2++) {
      Double ic2   = icoe[ip2];
      Double onedz = iexp[ip2];
      Double rho   = 1.0E0/onedz;
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
         * totally 4 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         ************************************************************/
        Double I_ESP_F3x_S_M1_vrr = PAX*I_ESP_D2x_S_M1_vrr-PRX*I_ESP_D2x_S_M2_vrr+2*oned2z*I_ESP_Px_S_M1_vrr-2*oned2z*I_ESP_Px_S_M2_vrr;
        Double I_ESP_F2xy_S_M1_vrr = PAY*I_ESP_D2x_S_M1_vrr-PRY*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_F2xz_S_M1_vrr = PAZ*I_ESP_D2x_S_M1_vrr-PRZ*I_ESP_D2x_S_M2_vrr;
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
         * totally 4 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         ************************************************************/
        Double I_ESP_F3x_S_vrr = PAX*I_ESP_D2x_S_vrr-PRX*I_ESP_D2x_S_M1_vrr+2*oned2z*I_ESP_Px_S_vrr-2*oned2z*I_ESP_Px_S_M1_vrr;
        Double I_ESP_F2xy_S_vrr = PAY*I_ESP_D2x_S_vrr-PRY*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F2xz_S_vrr = PAZ*I_ESP_D2x_S_vrr-PRZ*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F3y_S_vrr = PAY*I_ESP_D2y_S_vrr-PRY*I_ESP_D2y_S_M1_vrr+2*oned2z*I_ESP_Py_S_vrr-2*oned2z*I_ESP_Py_S_M1_vrr;
        Double I_ESP_F2yz_S_vrr = PAZ*I_ESP_D2y_S_vrr-PRZ*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_F3z_S_vrr = PAZ*I_ESP_D2z_S_vrr-PRZ*I_ESP_D2z_S_M1_vrr+2*oned2z*I_ESP_Pz_S_vrr-2*oned2z*I_ESP_Pz_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         ************************************************************/
        Double I_ESP_G4x_S_vrr = PAX*I_ESP_F3x_S_vrr-PRX*I_ESP_F3x_S_M1_vrr+3*oned2z*I_ESP_D2x_S_vrr-3*oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G3xy_S_vrr = PAY*I_ESP_F3x_S_vrr-PRY*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G3xz_S_vrr = PAZ*I_ESP_F3x_S_vrr-PRZ*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G2x2y_S_vrr = PAY*I_ESP_F2xy_S_vrr-PRY*I_ESP_F2xy_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G2x2z_S_vrr = PAZ*I_ESP_F2xz_S_vrr-PRZ*I_ESP_F2xz_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Gx3y_S_vrr = PAX*I_ESP_F3y_S_vrr-PRX*I_ESP_F3y_S_M1_vrr;
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
         * shell quartet name: SQ_ESP_N_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_N10x_S += I_ESP_N10x_S_vrr;
        I_ESP_N9xy_S += I_ESP_N9xy_S_vrr;
        I_ESP_N9xz_S += I_ESP_N9xz_S_vrr;
        I_ESP_N8x2y_S += I_ESP_N8x2y_S_vrr;
        I_ESP_N8xyz_S += I_ESP_N8xyz_S_vrr;
        I_ESP_N8x2z_S += I_ESP_N8x2z_S_vrr;
        I_ESP_N7x3y_S += I_ESP_N7x3y_S_vrr;
        I_ESP_N7x2yz_S += I_ESP_N7x2yz_S_vrr;
        I_ESP_N7xy2z_S += I_ESP_N7xy2z_S_vrr;
        I_ESP_N7x3z_S += I_ESP_N7x3z_S_vrr;
        I_ESP_N6x4y_S += I_ESP_N6x4y_S_vrr;
        I_ESP_N6x3yz_S += I_ESP_N6x3yz_S_vrr;
        I_ESP_N6x2y2z_S += I_ESP_N6x2y2z_S_vrr;
        I_ESP_N6xy3z_S += I_ESP_N6xy3z_S_vrr;
        I_ESP_N6x4z_S += I_ESP_N6x4z_S_vrr;
        I_ESP_N5x5y_S += I_ESP_N5x5y_S_vrr;
        I_ESP_N5x4yz_S += I_ESP_N5x4yz_S_vrr;
        I_ESP_N5x3y2z_S += I_ESP_N5x3y2z_S_vrr;
        I_ESP_N5x2y3z_S += I_ESP_N5x2y3z_S_vrr;
        I_ESP_N5xy4z_S += I_ESP_N5xy4z_S_vrr;
        I_ESP_N5x5z_S += I_ESP_N5x5z_S_vrr;
        I_ESP_N4x6y_S += I_ESP_N4x6y_S_vrr;
        I_ESP_N4x5yz_S += I_ESP_N4x5yz_S_vrr;
        I_ESP_N4x4y2z_S += I_ESP_N4x4y2z_S_vrr;
        I_ESP_N4x3y3z_S += I_ESP_N4x3y3z_S_vrr;
        I_ESP_N4x2y4z_S += I_ESP_N4x2y4z_S_vrr;
        I_ESP_N4xy5z_S += I_ESP_N4xy5z_S_vrr;
        I_ESP_N4x6z_S += I_ESP_N4x6z_S_vrr;
        I_ESP_N3x7y_S += I_ESP_N3x7y_S_vrr;
        I_ESP_N3x6yz_S += I_ESP_N3x6yz_S_vrr;
        I_ESP_N3x5y2z_S += I_ESP_N3x5y2z_S_vrr;
        I_ESP_N3x4y3z_S += I_ESP_N3x4y3z_S_vrr;
        I_ESP_N3x3y4z_S += I_ESP_N3x3y4z_S_vrr;
        I_ESP_N3x2y5z_S += I_ESP_N3x2y5z_S_vrr;
        I_ESP_N3xy6z_S += I_ESP_N3xy6z_S_vrr;
        I_ESP_N3x7z_S += I_ESP_N3x7z_S_vrr;
        I_ESP_N2x8y_S += I_ESP_N2x8y_S_vrr;
        I_ESP_N2x7yz_S += I_ESP_N2x7yz_S_vrr;
        I_ESP_N2x6y2z_S += I_ESP_N2x6y2z_S_vrr;
        I_ESP_N2x5y3z_S += I_ESP_N2x5y3z_S_vrr;
        I_ESP_N2x4y4z_S += I_ESP_N2x4y4z_S_vrr;
        I_ESP_N2x3y5z_S += I_ESP_N2x3y5z_S_vrr;
        I_ESP_N2x2y6z_S += I_ESP_N2x2y6z_S_vrr;
        I_ESP_N2xy7z_S += I_ESP_N2xy7z_S_vrr;
        I_ESP_N2x8z_S += I_ESP_N2x8z_S_vrr;
        I_ESP_Nx9y_S += I_ESP_Nx9y_S_vrr;
        I_ESP_Nx8yz_S += I_ESP_Nx8yz_S_vrr;
        I_ESP_Nx7y2z_S += I_ESP_Nx7y2z_S_vrr;
        I_ESP_Nx6y3z_S += I_ESP_Nx6y3z_S_vrr;
        I_ESP_Nx5y4z_S += I_ESP_Nx5y4z_S_vrr;
        I_ESP_Nx4y5z_S += I_ESP_Nx4y5z_S_vrr;
        I_ESP_Nx3y6z_S += I_ESP_Nx3y6z_S_vrr;
        I_ESP_Nx2y7z_S += I_ESP_Nx2y7z_S_vrr;
        I_ESP_Nxy8z_S += I_ESP_Nxy8z_S_vrr;
        I_ESP_Nx9z_S += I_ESP_Nx9z_S_vrr;
        I_ESP_N10y_S += I_ESP_N10y_S_vrr;
        I_ESP_N9yz_S += I_ESP_N9yz_S_vrr;
        I_ESP_N8y2z_S += I_ESP_N8y2z_S_vrr;
        I_ESP_N7y3z_S += I_ESP_N7y3z_S_vrr;
        I_ESP_N6y4z_S += I_ESP_N6y4z_S_vrr;
        I_ESP_N5y5z_S += I_ESP_N5y5z_S_vrr;
        I_ESP_N4y6z_S += I_ESP_N4y6z_S_vrr;
        I_ESP_N3y7z_S += I_ESP_N3y7z_S_vrr;
        I_ESP_N2y8z_S += I_ESP_N2y8z_S_vrr;
        I_ESP_Ny9z_S += I_ESP_Ny9z_S_vrr;
        I_ESP_N10z_S += I_ESP_N10z_S_vrr;

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
     * shell quartet name: SQ_ESP_K_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
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
    Double I_ESP_K7y_Px = I_ESP_Lx7y_S+ABX*I_ESP_K7y_S;
    Double I_ESP_K6yz_Px = I_ESP_Lx6yz_S+ABX*I_ESP_K6yz_S;
    Double I_ESP_K5y2z_Px = I_ESP_Lx5y2z_S+ABX*I_ESP_K5y2z_S;
    Double I_ESP_K4y3z_Px = I_ESP_Lx4y3z_S+ABX*I_ESP_K4y3z_S;
    Double I_ESP_K3y4z_Px = I_ESP_Lx3y4z_S+ABX*I_ESP_K3y4z_S;
    Double I_ESP_K2y5z_Px = I_ESP_Lx2y5z_S+ABX*I_ESP_K2y5z_S;
    Double I_ESP_Ky6z_Px = I_ESP_Lxy6z_S+ABX*I_ESP_Ky6z_S;
    Double I_ESP_K7z_Px = I_ESP_Lx7z_S+ABX*I_ESP_K7z_S;
    Double I_ESP_K7x_Py = I_ESP_L7xy_S+ABY*I_ESP_K7x_S;
    Double I_ESP_K6xy_Py = I_ESP_L6x2y_S+ABY*I_ESP_K6xy_S;
    Double I_ESP_K6xz_Py = I_ESP_L6xyz_S+ABY*I_ESP_K6xz_S;
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
    Double I_ESP_K7z_Py = I_ESP_Ly7z_S+ABY*I_ESP_K7z_S;
    Double I_ESP_K7x_Pz = I_ESP_L7xz_S+ABZ*I_ESP_K7x_S;
    Double I_ESP_K6xy_Pz = I_ESP_L6xyz_S+ABZ*I_ESP_K6xy_S;
    Double I_ESP_K6xz_Pz = I_ESP_L6x2z_S+ABZ*I_ESP_K6xz_S;
    Double I_ESP_K5x2y_Pz = I_ESP_L5x2yz_S+ABZ*I_ESP_K5x2y_S;
    Double I_ESP_K5xyz_Pz = I_ESP_L5xy2z_S+ABZ*I_ESP_K5xyz_S;
    Double I_ESP_K5x2z_Pz = I_ESP_L5x3z_S+ABZ*I_ESP_K5x2z_S;
    Double I_ESP_K4x3y_Pz = I_ESP_L4x3yz_S+ABZ*I_ESP_K4x3y_S;
    Double I_ESP_K4x2yz_Pz = I_ESP_L4x2y2z_S+ABZ*I_ESP_K4x2yz_S;
    Double I_ESP_K4xy2z_Pz = I_ESP_L4xy3z_S+ABZ*I_ESP_K4xy2z_S;
    Double I_ESP_K4x3z_Pz = I_ESP_L4x4z_S+ABZ*I_ESP_K4x3z_S;
    Double I_ESP_K3x4y_Pz = I_ESP_L3x4yz_S+ABZ*I_ESP_K3x4y_S;
    Double I_ESP_K3x3yz_Pz = I_ESP_L3x3y2z_S+ABZ*I_ESP_K3x3yz_S;
    Double I_ESP_K3x2y2z_Pz = I_ESP_L3x2y3z_S+ABZ*I_ESP_K3x2y2z_S;
    Double I_ESP_K3xy3z_Pz = I_ESP_L3xy4z_S+ABZ*I_ESP_K3xy3z_S;
    Double I_ESP_K3x4z_Pz = I_ESP_L3x5z_S+ABZ*I_ESP_K3x4z_S;
    Double I_ESP_K2x5y_Pz = I_ESP_L2x5yz_S+ABZ*I_ESP_K2x5y_S;
    Double I_ESP_K2x4yz_Pz = I_ESP_L2x4y2z_S+ABZ*I_ESP_K2x4yz_S;
    Double I_ESP_K2x3y2z_Pz = I_ESP_L2x3y3z_S+ABZ*I_ESP_K2x3y2z_S;
    Double I_ESP_K2x2y3z_Pz = I_ESP_L2x2y4z_S+ABZ*I_ESP_K2x2y3z_S;
    Double I_ESP_K2xy4z_Pz = I_ESP_L2xy5z_S+ABZ*I_ESP_K2xy4z_S;
    Double I_ESP_K2x5z_Pz = I_ESP_L2x6z_S+ABZ*I_ESP_K2x5z_S;
    Double I_ESP_Kx6y_Pz = I_ESP_Lx6yz_S+ABZ*I_ESP_Kx6y_S;
    Double I_ESP_Kx5yz_Pz = I_ESP_Lx5y2z_S+ABZ*I_ESP_Kx5yz_S;
    Double I_ESP_Kx4y2z_Pz = I_ESP_Lx4y3z_S+ABZ*I_ESP_Kx4y2z_S;
    Double I_ESP_Kx3y3z_Pz = I_ESP_Lx3y4z_S+ABZ*I_ESP_Kx3y3z_S;
    Double I_ESP_Kx2y4z_Pz = I_ESP_Lx2y5z_S+ABZ*I_ESP_Kx2y4z_S;
    Double I_ESP_Kxy5z_Pz = I_ESP_Lxy6z_S+ABZ*I_ESP_Kxy5z_S;
    Double I_ESP_Kx6z_Pz = I_ESP_Lx7z_S+ABZ*I_ESP_Kx6z_S;
    Double I_ESP_K7y_Pz = I_ESP_L7yz_S+ABZ*I_ESP_K7y_S;
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
     * totally 84 integrals are omitted 
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
    Double I_ESP_H5x_F2xy = I_ESP_I5xy_D2x+ABY*I_ESP_H5x_D2x;
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
    Double I_ESP_H5x_F2xz = I_ESP_I5xz_D2x+ABZ*I_ESP_H5x_D2x;
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
    Double I_ESP_H5y_F2yz = I_ESP_I5yz_D2y+ABZ*I_ESP_H5y_D2y;
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
     * shell quartet name: SQ_ESP_L_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 14 integrals are omitted 
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
    Double I_ESP_Lx7y_Px = I_ESP_M2x7y_S+ABX*I_ESP_Lx7y_S;
    Double I_ESP_Lx6yz_Px = I_ESP_M2x6yz_S+ABX*I_ESP_Lx6yz_S;
    Double I_ESP_Lx5y2z_Px = I_ESP_M2x5y2z_S+ABX*I_ESP_Lx5y2z_S;
    Double I_ESP_Lx4y3z_Px = I_ESP_M2x4y3z_S+ABX*I_ESP_Lx4y3z_S;
    Double I_ESP_Lx3y4z_Px = I_ESP_M2x3y4z_S+ABX*I_ESP_Lx3y4z_S;
    Double I_ESP_Lx2y5z_Px = I_ESP_M2x2y5z_S+ABX*I_ESP_Lx2y5z_S;
    Double I_ESP_Lxy6z_Px = I_ESP_M2xy6z_S+ABX*I_ESP_Lxy6z_S;
    Double I_ESP_Lx7z_Px = I_ESP_M2x7z_S+ABX*I_ESP_Lx7z_S;
    Double I_ESP_L7yz_Px = I_ESP_Mx7yz_S+ABX*I_ESP_L7yz_S;
    Double I_ESP_L6y2z_Px = I_ESP_Mx6y2z_S+ABX*I_ESP_L6y2z_S;
    Double I_ESP_L5y3z_Px = I_ESP_Mx5y3z_S+ABX*I_ESP_L5y3z_S;
    Double I_ESP_L4y4z_Px = I_ESP_Mx4y4z_S+ABX*I_ESP_L4y4z_S;
    Double I_ESP_L3y5z_Px = I_ESP_Mx3y5z_S+ABX*I_ESP_L3y5z_S;
    Double I_ESP_L2y6z_Px = I_ESP_Mx2y6z_S+ABX*I_ESP_L2y6z_S;
    Double I_ESP_Ly7z_Px = I_ESP_Mxy7z_S+ABX*I_ESP_Ly7z_S;
    Double I_ESP_L7xy_Py = I_ESP_M7x2y_S+ABY*I_ESP_L7xy_S;
    Double I_ESP_L6x2y_Py = I_ESP_M6x3y_S+ABY*I_ESP_L6x2y_S;
    Double I_ESP_L6xyz_Py = I_ESP_M6x2yz_S+ABY*I_ESP_L6xyz_S;
    Double I_ESP_L6x2z_Py = I_ESP_M6xy2z_S+ABY*I_ESP_L6x2z_S;
    Double I_ESP_L5x3y_Py = I_ESP_M5x4y_S+ABY*I_ESP_L5x3y_S;
    Double I_ESP_L5x2yz_Py = I_ESP_M5x3yz_S+ABY*I_ESP_L5x2yz_S;
    Double I_ESP_L5xy2z_Py = I_ESP_M5x2y2z_S+ABY*I_ESP_L5xy2z_S;
    Double I_ESP_L5x3z_Py = I_ESP_M5xy3z_S+ABY*I_ESP_L5x3z_S;
    Double I_ESP_L4x4y_Py = I_ESP_M4x5y_S+ABY*I_ESP_L4x4y_S;
    Double I_ESP_L4x3yz_Py = I_ESP_M4x4yz_S+ABY*I_ESP_L4x3yz_S;
    Double I_ESP_L4x2y2z_Py = I_ESP_M4x3y2z_S+ABY*I_ESP_L4x2y2z_S;
    Double I_ESP_L4xy3z_Py = I_ESP_M4x2y3z_S+ABY*I_ESP_L4xy3z_S;
    Double I_ESP_L4x4z_Py = I_ESP_M4xy4z_S+ABY*I_ESP_L4x4z_S;
    Double I_ESP_L3x5y_Py = I_ESP_M3x6y_S+ABY*I_ESP_L3x5y_S;
    Double I_ESP_L3x4yz_Py = I_ESP_M3x5yz_S+ABY*I_ESP_L3x4yz_S;
    Double I_ESP_L3x3y2z_Py = I_ESP_M3x4y2z_S+ABY*I_ESP_L3x3y2z_S;
    Double I_ESP_L3x2y3z_Py = I_ESP_M3x3y3z_S+ABY*I_ESP_L3x2y3z_S;
    Double I_ESP_L3xy4z_Py = I_ESP_M3x2y4z_S+ABY*I_ESP_L3xy4z_S;
    Double I_ESP_L3x5z_Py = I_ESP_M3xy5z_S+ABY*I_ESP_L3x5z_S;
    Double I_ESP_L2x6y_Py = I_ESP_M2x7y_S+ABY*I_ESP_L2x6y_S;
    Double I_ESP_L2x5yz_Py = I_ESP_M2x6yz_S+ABY*I_ESP_L2x5yz_S;
    Double I_ESP_L2x4y2z_Py = I_ESP_M2x5y2z_S+ABY*I_ESP_L2x4y2z_S;
    Double I_ESP_L2x3y3z_Py = I_ESP_M2x4y3z_S+ABY*I_ESP_L2x3y3z_S;
    Double I_ESP_L2x2y4z_Py = I_ESP_M2x3y4z_S+ABY*I_ESP_L2x2y4z_S;
    Double I_ESP_L2xy5z_Py = I_ESP_M2x2y5z_S+ABY*I_ESP_L2xy5z_S;
    Double I_ESP_L2x6z_Py = I_ESP_M2xy6z_S+ABY*I_ESP_L2x6z_S;
    Double I_ESP_Lx7y_Py = I_ESP_Mx8y_S+ABY*I_ESP_Lx7y_S;
    Double I_ESP_Lx6yz_Py = I_ESP_Mx7yz_S+ABY*I_ESP_Lx6yz_S;
    Double I_ESP_Lx5y2z_Py = I_ESP_Mx6y2z_S+ABY*I_ESP_Lx5y2z_S;
    Double I_ESP_Lx4y3z_Py = I_ESP_Mx5y3z_S+ABY*I_ESP_Lx4y3z_S;
    Double I_ESP_Lx3y4z_Py = I_ESP_Mx4y4z_S+ABY*I_ESP_Lx3y4z_S;
    Double I_ESP_Lx2y5z_Py = I_ESP_Mx3y5z_S+ABY*I_ESP_Lx2y5z_S;
    Double I_ESP_Lxy6z_Py = I_ESP_Mx2y6z_S+ABY*I_ESP_Lxy6z_S;
    Double I_ESP_Lx7z_Py = I_ESP_Mxy7z_S+ABY*I_ESP_Lx7z_S;
    Double I_ESP_L8y_Py = I_ESP_M9y_S+ABY*I_ESP_L8y_S;
    Double I_ESP_L7yz_Py = I_ESP_M8yz_S+ABY*I_ESP_L7yz_S;
    Double I_ESP_L6y2z_Py = I_ESP_M7y2z_S+ABY*I_ESP_L6y2z_S;
    Double I_ESP_L5y3z_Py = I_ESP_M6y3z_S+ABY*I_ESP_L5y3z_S;
    Double I_ESP_L4y4z_Py = I_ESP_M5y4z_S+ABY*I_ESP_L4y4z_S;
    Double I_ESP_L3y5z_Py = I_ESP_M4y5z_S+ABY*I_ESP_L3y5z_S;
    Double I_ESP_L2y6z_Py = I_ESP_M3y6z_S+ABY*I_ESP_L2y6z_S;
    Double I_ESP_Ly7z_Py = I_ESP_M2y7z_S+ABY*I_ESP_Ly7z_S;
    Double I_ESP_L7xz_Pz = I_ESP_M7x2z_S+ABZ*I_ESP_L7xz_S;
    Double I_ESP_L6xyz_Pz = I_ESP_M6xy2z_S+ABZ*I_ESP_L6xyz_S;
    Double I_ESP_L6x2z_Pz = I_ESP_M6x3z_S+ABZ*I_ESP_L6x2z_S;
    Double I_ESP_L5x2yz_Pz = I_ESP_M5x2y2z_S+ABZ*I_ESP_L5x2yz_S;
    Double I_ESP_L5xy2z_Pz = I_ESP_M5xy3z_S+ABZ*I_ESP_L5xy2z_S;
    Double I_ESP_L5x3z_Pz = I_ESP_M5x4z_S+ABZ*I_ESP_L5x3z_S;
    Double I_ESP_L4x3yz_Pz = I_ESP_M4x3y2z_S+ABZ*I_ESP_L4x3yz_S;
    Double I_ESP_L4x2y2z_Pz = I_ESP_M4x2y3z_S+ABZ*I_ESP_L4x2y2z_S;
    Double I_ESP_L4xy3z_Pz = I_ESP_M4xy4z_S+ABZ*I_ESP_L4xy3z_S;
    Double I_ESP_L4x4z_Pz = I_ESP_M4x5z_S+ABZ*I_ESP_L4x4z_S;
    Double I_ESP_L3x4yz_Pz = I_ESP_M3x4y2z_S+ABZ*I_ESP_L3x4yz_S;
    Double I_ESP_L3x3y2z_Pz = I_ESP_M3x3y3z_S+ABZ*I_ESP_L3x3y2z_S;
    Double I_ESP_L3x2y3z_Pz = I_ESP_M3x2y4z_S+ABZ*I_ESP_L3x2y3z_S;
    Double I_ESP_L3xy4z_Pz = I_ESP_M3xy5z_S+ABZ*I_ESP_L3xy4z_S;
    Double I_ESP_L3x5z_Pz = I_ESP_M3x6z_S+ABZ*I_ESP_L3x5z_S;
    Double I_ESP_L2x5yz_Pz = I_ESP_M2x5y2z_S+ABZ*I_ESP_L2x5yz_S;
    Double I_ESP_L2x4y2z_Pz = I_ESP_M2x4y3z_S+ABZ*I_ESP_L2x4y2z_S;
    Double I_ESP_L2x3y3z_Pz = I_ESP_M2x3y4z_S+ABZ*I_ESP_L2x3y3z_S;
    Double I_ESP_L2x2y4z_Pz = I_ESP_M2x2y5z_S+ABZ*I_ESP_L2x2y4z_S;
    Double I_ESP_L2xy5z_Pz = I_ESP_M2xy6z_S+ABZ*I_ESP_L2xy5z_S;
    Double I_ESP_L2x6z_Pz = I_ESP_M2x7z_S+ABZ*I_ESP_L2x6z_S;
    Double I_ESP_Lx6yz_Pz = I_ESP_Mx6y2z_S+ABZ*I_ESP_Lx6yz_S;
    Double I_ESP_Lx5y2z_Pz = I_ESP_Mx5y3z_S+ABZ*I_ESP_Lx5y2z_S;
    Double I_ESP_Lx4y3z_Pz = I_ESP_Mx4y4z_S+ABZ*I_ESP_Lx4y3z_S;
    Double I_ESP_Lx3y4z_Pz = I_ESP_Mx3y5z_S+ABZ*I_ESP_Lx3y4z_S;
    Double I_ESP_Lx2y5z_Pz = I_ESP_Mx2y6z_S+ABZ*I_ESP_Lx2y5z_S;
    Double I_ESP_Lxy6z_Pz = I_ESP_Mxy7z_S+ABZ*I_ESP_Lxy6z_S;
    Double I_ESP_Lx7z_Pz = I_ESP_Mx8z_S+ABZ*I_ESP_Lx7z_S;
    Double I_ESP_L7yz_Pz = I_ESP_M7y2z_S+ABZ*I_ESP_L7yz_S;
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
     * totally 108 integrals are omitted 
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
    Double I_ESP_K7y_D2x = I_ESP_Lx7y_Px+ABX*I_ESP_K7y_Px;
    Double I_ESP_K6yz_D2x = I_ESP_Lx6yz_Px+ABX*I_ESP_K6yz_Px;
    Double I_ESP_K5y2z_D2x = I_ESP_Lx5y2z_Px+ABX*I_ESP_K5y2z_Px;
    Double I_ESP_K4y3z_D2x = I_ESP_Lx4y3z_Px+ABX*I_ESP_K4y3z_Px;
    Double I_ESP_K3y4z_D2x = I_ESP_Lx3y4z_Px+ABX*I_ESP_K3y4z_Px;
    Double I_ESP_K2y5z_D2x = I_ESP_Lx2y5z_Px+ABX*I_ESP_K2y5z_Px;
    Double I_ESP_Ky6z_D2x = I_ESP_Lxy6z_Px+ABX*I_ESP_Ky6z_Px;
    Double I_ESP_K7z_D2x = I_ESP_Lx7z_Px+ABX*I_ESP_K7z_Px;
    Double I_ESP_K7x_D2y = I_ESP_L7xy_Py+ABY*I_ESP_K7x_Py;
    Double I_ESP_K6xy_D2y = I_ESP_L6x2y_Py+ABY*I_ESP_K6xy_Py;
    Double I_ESP_K6xz_D2y = I_ESP_L6xyz_Py+ABY*I_ESP_K6xz_Py;
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
    Double I_ESP_K7z_D2y = I_ESP_Ly7z_Py+ABY*I_ESP_K7z_Py;
    Double I_ESP_K7x_D2z = I_ESP_L7xz_Pz+ABZ*I_ESP_K7x_Pz;
    Double I_ESP_K6xy_D2z = I_ESP_L6xyz_Pz+ABZ*I_ESP_K6xy_Pz;
    Double I_ESP_K6xz_D2z = I_ESP_L6x2z_Pz+ABZ*I_ESP_K6xz_Pz;
    Double I_ESP_K5x2y_D2z = I_ESP_L5x2yz_Pz+ABZ*I_ESP_K5x2y_Pz;
    Double I_ESP_K5xyz_D2z = I_ESP_L5xy2z_Pz+ABZ*I_ESP_K5xyz_Pz;
    Double I_ESP_K5x2z_D2z = I_ESP_L5x3z_Pz+ABZ*I_ESP_K5x2z_Pz;
    Double I_ESP_K4x3y_D2z = I_ESP_L4x3yz_Pz+ABZ*I_ESP_K4x3y_Pz;
    Double I_ESP_K4x2yz_D2z = I_ESP_L4x2y2z_Pz+ABZ*I_ESP_K4x2yz_Pz;
    Double I_ESP_K4xy2z_D2z = I_ESP_L4xy3z_Pz+ABZ*I_ESP_K4xy2z_Pz;
    Double I_ESP_K4x3z_D2z = I_ESP_L4x4z_Pz+ABZ*I_ESP_K4x3z_Pz;
    Double I_ESP_K3x4y_D2z = I_ESP_L3x4yz_Pz+ABZ*I_ESP_K3x4y_Pz;
    Double I_ESP_K3x3yz_D2z = I_ESP_L3x3y2z_Pz+ABZ*I_ESP_K3x3yz_Pz;
    Double I_ESP_K3x2y2z_D2z = I_ESP_L3x2y3z_Pz+ABZ*I_ESP_K3x2y2z_Pz;
    Double I_ESP_K3xy3z_D2z = I_ESP_L3xy4z_Pz+ABZ*I_ESP_K3xy3z_Pz;
    Double I_ESP_K3x4z_D2z = I_ESP_L3x5z_Pz+ABZ*I_ESP_K3x4z_Pz;
    Double I_ESP_K2x5y_D2z = I_ESP_L2x5yz_Pz+ABZ*I_ESP_K2x5y_Pz;
    Double I_ESP_K2x4yz_D2z = I_ESP_L2x4y2z_Pz+ABZ*I_ESP_K2x4yz_Pz;
    Double I_ESP_K2x3y2z_D2z = I_ESP_L2x3y3z_Pz+ABZ*I_ESP_K2x3y2z_Pz;
    Double I_ESP_K2x2y3z_D2z = I_ESP_L2x2y4z_Pz+ABZ*I_ESP_K2x2y3z_Pz;
    Double I_ESP_K2xy4z_D2z = I_ESP_L2xy5z_Pz+ABZ*I_ESP_K2xy4z_Pz;
    Double I_ESP_K2x5z_D2z = I_ESP_L2x6z_Pz+ABZ*I_ESP_K2x5z_Pz;
    Double I_ESP_Kx6y_D2z = I_ESP_Lx6yz_Pz+ABZ*I_ESP_Kx6y_Pz;
    Double I_ESP_Kx5yz_D2z = I_ESP_Lx5y2z_Pz+ABZ*I_ESP_Kx5yz_Pz;
    Double I_ESP_Kx4y2z_D2z = I_ESP_Lx4y3z_Pz+ABZ*I_ESP_Kx4y2z_Pz;
    Double I_ESP_Kx3y3z_D2z = I_ESP_Lx3y4z_Pz+ABZ*I_ESP_Kx3y3z_Pz;
    Double I_ESP_Kx2y4z_D2z = I_ESP_Lx2y5z_Pz+ABZ*I_ESP_Kx2y4z_Pz;
    Double I_ESP_Kxy5z_D2z = I_ESP_Lxy6z_Pz+ABZ*I_ESP_Kxy5z_Pz;
    Double I_ESP_Kx6z_D2z = I_ESP_Lx7z_Pz+ABZ*I_ESP_Kx6z_Pz;
    Double I_ESP_K7y_D2z = I_ESP_L7yz_Pz+ABZ*I_ESP_K7y_Pz;
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
     * totally 115 integrals are omitted 
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
    Double I_ESP_I5xy_F2xy = I_ESP_K5x2y_D2x+ABY*I_ESP_I5xy_D2x;
    Double I_ESP_I5xz_F2xy = I_ESP_K5xyz_D2x+ABY*I_ESP_I5xz_D2x;
    Double I_ESP_I4x2y_F2xy = I_ESP_K4x3y_D2x+ABY*I_ESP_I4x2y_D2x;
    Double I_ESP_I4xyz_F2xy = I_ESP_K4x2yz_D2x+ABY*I_ESP_I4xyz_D2x;
    Double I_ESP_I4x2z_F2xy = I_ESP_K4xy2z_D2x+ABY*I_ESP_I4x2z_D2x;
    Double I_ESP_I3x3y_F2xy = I_ESP_K3x4y_D2x+ABY*I_ESP_I3x3y_D2x;
    Double I_ESP_I3x2yz_F2xy = I_ESP_K3x3yz_D2x+ABY*I_ESP_I3x2yz_D2x;
    Double I_ESP_I3xy2z_F2xy = I_ESP_K3x2y2z_D2x+ABY*I_ESP_I3xy2z_D2x;
    Double I_ESP_I3x3z_F2xy = I_ESP_K3xy3z_D2x+ABY*I_ESP_I3x3z_D2x;
    Double I_ESP_I2x4y_F2xy = I_ESP_K2x5y_D2x+ABY*I_ESP_I2x4y_D2x;
    Double I_ESP_I2x3yz_F2xy = I_ESP_K2x4yz_D2x+ABY*I_ESP_I2x3yz_D2x;
    Double I_ESP_I2x2y2z_F2xy = I_ESP_K2x3y2z_D2x+ABY*I_ESP_I2x2y2z_D2x;
    Double I_ESP_I2xy3z_F2xy = I_ESP_K2x2y3z_D2x+ABY*I_ESP_I2xy3z_D2x;
    Double I_ESP_I2x4z_F2xy = I_ESP_K2xy4z_D2x+ABY*I_ESP_I2x4z_D2x;
    Double I_ESP_Ix5y_F2xy = I_ESP_Kx6y_D2x+ABY*I_ESP_Ix5y_D2x;
    Double I_ESP_Ix4yz_F2xy = I_ESP_Kx5yz_D2x+ABY*I_ESP_Ix4yz_D2x;
    Double I_ESP_Ix3y2z_F2xy = I_ESP_Kx4y2z_D2x+ABY*I_ESP_Ix3y2z_D2x;
    Double I_ESP_Ix2y3z_F2xy = I_ESP_Kx3y3z_D2x+ABY*I_ESP_Ix2y3z_D2x;
    Double I_ESP_Ixy4z_F2xy = I_ESP_Kx2y4z_D2x+ABY*I_ESP_Ixy4z_D2x;
    Double I_ESP_Ix5z_F2xy = I_ESP_Kxy5z_D2x+ABY*I_ESP_Ix5z_D2x;
    Double I_ESP_I6y_F2xy = I_ESP_K7y_D2x+ABY*I_ESP_I6y_D2x;
    Double I_ESP_I5yz_F2xy = I_ESP_K6yz_D2x+ABY*I_ESP_I5yz_D2x;
    Double I_ESP_I4y2z_F2xy = I_ESP_K5y2z_D2x+ABY*I_ESP_I4y2z_D2x;
    Double I_ESP_I3y3z_F2xy = I_ESP_K4y3z_D2x+ABY*I_ESP_I3y3z_D2x;
    Double I_ESP_I2y4z_F2xy = I_ESP_K3y4z_D2x+ABY*I_ESP_I2y4z_D2x;
    Double I_ESP_Iy5z_F2xy = I_ESP_K2y5z_D2x+ABY*I_ESP_Iy5z_D2x;
    Double I_ESP_I6z_F2xy = I_ESP_Ky6z_D2x+ABY*I_ESP_I6z_D2x;
    Double I_ESP_I5xy_F2xz = I_ESP_K5xyz_D2x+ABZ*I_ESP_I5xy_D2x;
    Double I_ESP_I5xz_F2xz = I_ESP_K5x2z_D2x+ABZ*I_ESP_I5xz_D2x;
    Double I_ESP_I4x2y_F2xz = I_ESP_K4x2yz_D2x+ABZ*I_ESP_I4x2y_D2x;
    Double I_ESP_I4xyz_F2xz = I_ESP_K4xy2z_D2x+ABZ*I_ESP_I4xyz_D2x;
    Double I_ESP_I4x2z_F2xz = I_ESP_K4x3z_D2x+ABZ*I_ESP_I4x2z_D2x;
    Double I_ESP_I3x3y_F2xz = I_ESP_K3x3yz_D2x+ABZ*I_ESP_I3x3y_D2x;
    Double I_ESP_I3x2yz_F2xz = I_ESP_K3x2y2z_D2x+ABZ*I_ESP_I3x2yz_D2x;
    Double I_ESP_I3xy2z_F2xz = I_ESP_K3xy3z_D2x+ABZ*I_ESP_I3xy2z_D2x;
    Double I_ESP_I3x3z_F2xz = I_ESP_K3x4z_D2x+ABZ*I_ESP_I3x3z_D2x;
    Double I_ESP_I2x4y_F2xz = I_ESP_K2x4yz_D2x+ABZ*I_ESP_I2x4y_D2x;
    Double I_ESP_I2x3yz_F2xz = I_ESP_K2x3y2z_D2x+ABZ*I_ESP_I2x3yz_D2x;
    Double I_ESP_I2x2y2z_F2xz = I_ESP_K2x2y3z_D2x+ABZ*I_ESP_I2x2y2z_D2x;
    Double I_ESP_I2xy3z_F2xz = I_ESP_K2xy4z_D2x+ABZ*I_ESP_I2xy3z_D2x;
    Double I_ESP_I2x4z_F2xz = I_ESP_K2x5z_D2x+ABZ*I_ESP_I2x4z_D2x;
    Double I_ESP_Ix5y_F2xz = I_ESP_Kx5yz_D2x+ABZ*I_ESP_Ix5y_D2x;
    Double I_ESP_Ix4yz_F2xz = I_ESP_Kx4y2z_D2x+ABZ*I_ESP_Ix4yz_D2x;
    Double I_ESP_Ix3y2z_F2xz = I_ESP_Kx3y3z_D2x+ABZ*I_ESP_Ix3y2z_D2x;
    Double I_ESP_Ix2y3z_F2xz = I_ESP_Kx2y4z_D2x+ABZ*I_ESP_Ix2y3z_D2x;
    Double I_ESP_Ixy4z_F2xz = I_ESP_Kxy5z_D2x+ABZ*I_ESP_Ixy4z_D2x;
    Double I_ESP_Ix5z_F2xz = I_ESP_Kx6z_D2x+ABZ*I_ESP_Ix5z_D2x;
    Double I_ESP_I6y_F2xz = I_ESP_K6yz_D2x+ABZ*I_ESP_I6y_D2x;
    Double I_ESP_I5yz_F2xz = I_ESP_K5y2z_D2x+ABZ*I_ESP_I5yz_D2x;
    Double I_ESP_I4y2z_F2xz = I_ESP_K4y3z_D2x+ABZ*I_ESP_I4y2z_D2x;
    Double I_ESP_I3y3z_F2xz = I_ESP_K3y4z_D2x+ABZ*I_ESP_I3y3z_D2x;
    Double I_ESP_I2y4z_F2xz = I_ESP_K2y5z_D2x+ABZ*I_ESP_I2y4z_D2x;
    Double I_ESP_Iy5z_F2xz = I_ESP_Ky6z_D2x+ABZ*I_ESP_Iy5z_D2x;
    Double I_ESP_I6z_F2xz = I_ESP_K7z_D2x+ABZ*I_ESP_I6z_D2x;
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
    Double I_ESP_I6x_F2yz = I_ESP_K6xz_D2y+ABZ*I_ESP_I6x_D2y;
    Double I_ESP_I5xy_F2yz = I_ESP_K5xyz_D2y+ABZ*I_ESP_I5xy_D2y;
    Double I_ESP_I5xz_F2yz = I_ESP_K5x2z_D2y+ABZ*I_ESP_I5xz_D2y;
    Double I_ESP_I4x2y_F2yz = I_ESP_K4x2yz_D2y+ABZ*I_ESP_I4x2y_D2y;
    Double I_ESP_I4xyz_F2yz = I_ESP_K4xy2z_D2y+ABZ*I_ESP_I4xyz_D2y;
    Double I_ESP_I4x2z_F2yz = I_ESP_K4x3z_D2y+ABZ*I_ESP_I4x2z_D2y;
    Double I_ESP_I3x3y_F2yz = I_ESP_K3x3yz_D2y+ABZ*I_ESP_I3x3y_D2y;
    Double I_ESP_I3x2yz_F2yz = I_ESP_K3x2y2z_D2y+ABZ*I_ESP_I3x2yz_D2y;
    Double I_ESP_I3xy2z_F2yz = I_ESP_K3xy3z_D2y+ABZ*I_ESP_I3xy2z_D2y;
    Double I_ESP_I3x3z_F2yz = I_ESP_K3x4z_D2y+ABZ*I_ESP_I3x3z_D2y;
    Double I_ESP_I2x4y_F2yz = I_ESP_K2x4yz_D2y+ABZ*I_ESP_I2x4y_D2y;
    Double I_ESP_I2x3yz_F2yz = I_ESP_K2x3y2z_D2y+ABZ*I_ESP_I2x3yz_D2y;
    Double I_ESP_I2x2y2z_F2yz = I_ESP_K2x2y3z_D2y+ABZ*I_ESP_I2x2y2z_D2y;
    Double I_ESP_I2xy3z_F2yz = I_ESP_K2xy4z_D2y+ABZ*I_ESP_I2xy3z_D2y;
    Double I_ESP_I2x4z_F2yz = I_ESP_K2x5z_D2y+ABZ*I_ESP_I2x4z_D2y;
    Double I_ESP_Ix5y_F2yz = I_ESP_Kx5yz_D2y+ABZ*I_ESP_Ix5y_D2y;
    Double I_ESP_Ix4yz_F2yz = I_ESP_Kx4y2z_D2y+ABZ*I_ESP_Ix4yz_D2y;
    Double I_ESP_Ix3y2z_F2yz = I_ESP_Kx3y3z_D2y+ABZ*I_ESP_Ix3y2z_D2y;
    Double I_ESP_Ix2y3z_F2yz = I_ESP_Kx2y4z_D2y+ABZ*I_ESP_Ix2y3z_D2y;
    Double I_ESP_Ixy4z_F2yz = I_ESP_Kxy5z_D2y+ABZ*I_ESP_Ixy4z_D2y;
    Double I_ESP_Ix5z_F2yz = I_ESP_Kx6z_D2y+ABZ*I_ESP_Ix5z_D2y;
    Double I_ESP_I5yz_F2yz = I_ESP_K5y2z_D2y+ABZ*I_ESP_I5yz_D2y;
    Double I_ESP_I4y2z_F2yz = I_ESP_K4y3z_D2y+ABZ*I_ESP_I4y2z_D2y;
    Double I_ESP_I3y3z_F2yz = I_ESP_K3y4z_D2y+ABZ*I_ESP_I3y3z_D2y;
    Double I_ESP_I2y4z_F2yz = I_ESP_K2y5z_D2y+ABZ*I_ESP_I2y4z_D2y;
    Double I_ESP_Iy5z_F2yz = I_ESP_Ky6z_D2y+ABZ*I_ESP_Iy5z_D2y;
    Double I_ESP_I6z_F2yz = I_ESP_K7z_D2y+ABZ*I_ESP_I6z_D2y;
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
     * totally 63 integrals are omitted 
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
    Double I_ESP_H5x_G3xy = I_ESP_I5xy_F3x+ABY*I_ESP_H5x_F3x;
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
    Double I_ESP_H5x_G3xz = I_ESP_I5xz_F3x+ABZ*I_ESP_H5x_F3x;
    Double I_ESP_H4xy_G3xz = I_ESP_I4xyz_F3x+ABZ*I_ESP_H4xy_F3x;
    Double I_ESP_H4xz_G3xz = I_ESP_I4x2z_F3x+ABZ*I_ESP_H4xz_F3x;
    Double I_ESP_H3x2y_G3xz = I_ESP_I3x2yz_F3x+ABZ*I_ESP_H3x2y_F3x;
    Double I_ESP_H3xyz_G3xz = I_ESP_I3xy2z_F3x+ABZ*I_ESP_H3xyz_F3x;
    Double I_ESP_H3x2z_G3xz = I_ESP_I3x3z_F3x+ABZ*I_ESP_H3x2z_F3x;
    Double I_ESP_H2x3y_G3xz = I_ESP_I2x3yz_F3x+ABZ*I_ESP_H2x3y_F3x;
    Double I_ESP_H2x2yz_G3xz = I_ESP_I2x2y2z_F3x+ABZ*I_ESP_H2x2yz_F3x;
    Double I_ESP_H2xy2z_G3xz = I_ESP_I2xy3z_F3x+ABZ*I_ESP_H2xy2z_F3x;
    Double I_ESP_H2x3z_G3xz = I_ESP_I2x4z_F3x+ABZ*I_ESP_H2x3z_F3x;
    Double I_ESP_Hx4y_G3xz = I_ESP_Ix4yz_F3x+ABZ*I_ESP_Hx4y_F3x;
    Double I_ESP_Hx3yz_G3xz = I_ESP_Ix3y2z_F3x+ABZ*I_ESP_Hx3yz_F3x;
    Double I_ESP_Hx2y2z_G3xz = I_ESP_Ix2y3z_F3x+ABZ*I_ESP_Hx2y2z_F3x;
    Double I_ESP_Hxy3z_G3xz = I_ESP_Ixy4z_F3x+ABZ*I_ESP_Hxy3z_F3x;
    Double I_ESP_Hx4z_G3xz = I_ESP_Ix5z_F3x+ABZ*I_ESP_Hx4z_F3x;
    Double I_ESP_H5y_G3xz = I_ESP_I5yz_F3x+ABZ*I_ESP_H5y_F3x;
    Double I_ESP_H4yz_G3xz = I_ESP_I4y2z_F3x+ABZ*I_ESP_H4yz_F3x;
    Double I_ESP_H3y2z_G3xz = I_ESP_I3y3z_F3x+ABZ*I_ESP_H3y2z_F3x;
    Double I_ESP_H2y3z_G3xz = I_ESP_I2y4z_F3x+ABZ*I_ESP_H2y3z_F3x;
    Double I_ESP_Hy4z_G3xz = I_ESP_Iy5z_F3x+ABZ*I_ESP_Hy4z_F3x;
    Double I_ESP_H5z_G3xz = I_ESP_I6z_F3x+ABZ*I_ESP_H5z_F3x;
    Double I_ESP_H5x_G2x2y = I_ESP_I5xy_F2xy+ABY*I_ESP_H5x_F2xy;
    Double I_ESP_H4xy_G2x2y = I_ESP_I4x2y_F2xy+ABY*I_ESP_H4xy_F2xy;
    Double I_ESP_H4xz_G2x2y = I_ESP_I4xyz_F2xy+ABY*I_ESP_H4xz_F2xy;
    Double I_ESP_H3x2y_G2x2y = I_ESP_I3x3y_F2xy+ABY*I_ESP_H3x2y_F2xy;
    Double I_ESP_H3xyz_G2x2y = I_ESP_I3x2yz_F2xy+ABY*I_ESP_H3xyz_F2xy;
    Double I_ESP_H3x2z_G2x2y = I_ESP_I3xy2z_F2xy+ABY*I_ESP_H3x2z_F2xy;
    Double I_ESP_H2x3y_G2x2y = I_ESP_I2x4y_F2xy+ABY*I_ESP_H2x3y_F2xy;
    Double I_ESP_H2x2yz_G2x2y = I_ESP_I2x3yz_F2xy+ABY*I_ESP_H2x2yz_F2xy;
    Double I_ESP_H2xy2z_G2x2y = I_ESP_I2x2y2z_F2xy+ABY*I_ESP_H2xy2z_F2xy;
    Double I_ESP_H2x3z_G2x2y = I_ESP_I2xy3z_F2xy+ABY*I_ESP_H2x3z_F2xy;
    Double I_ESP_Hx4y_G2x2y = I_ESP_Ix5y_F2xy+ABY*I_ESP_Hx4y_F2xy;
    Double I_ESP_Hx3yz_G2x2y = I_ESP_Ix4yz_F2xy+ABY*I_ESP_Hx3yz_F2xy;
    Double I_ESP_Hx2y2z_G2x2y = I_ESP_Ix3y2z_F2xy+ABY*I_ESP_Hx2y2z_F2xy;
    Double I_ESP_Hxy3z_G2x2y = I_ESP_Ix2y3z_F2xy+ABY*I_ESP_Hxy3z_F2xy;
    Double I_ESP_Hx4z_G2x2y = I_ESP_Ixy4z_F2xy+ABY*I_ESP_Hx4z_F2xy;
    Double I_ESP_H5y_G2x2y = I_ESP_I6y_F2xy+ABY*I_ESP_H5y_F2xy;
    Double I_ESP_H4yz_G2x2y = I_ESP_I5yz_F2xy+ABY*I_ESP_H4yz_F2xy;
    Double I_ESP_H3y2z_G2x2y = I_ESP_I4y2z_F2xy+ABY*I_ESP_H3y2z_F2xy;
    Double I_ESP_H2y3z_G2x2y = I_ESP_I3y3z_F2xy+ABY*I_ESP_H2y3z_F2xy;
    Double I_ESP_Hy4z_G2x2y = I_ESP_I2y4z_F2xy+ABY*I_ESP_Hy4z_F2xy;
    Double I_ESP_H5z_G2x2y = I_ESP_Iy5z_F2xy+ABY*I_ESP_H5z_F2xy;
    Double I_ESP_H5x_G2x2z = I_ESP_I5xz_F2xz+ABZ*I_ESP_H5x_F2xz;
    Double I_ESP_H4xy_G2x2z = I_ESP_I4xyz_F2xz+ABZ*I_ESP_H4xy_F2xz;
    Double I_ESP_H4xz_G2x2z = I_ESP_I4x2z_F2xz+ABZ*I_ESP_H4xz_F2xz;
    Double I_ESP_H3x2y_G2x2z = I_ESP_I3x2yz_F2xz+ABZ*I_ESP_H3x2y_F2xz;
    Double I_ESP_H3xyz_G2x2z = I_ESP_I3xy2z_F2xz+ABZ*I_ESP_H3xyz_F2xz;
    Double I_ESP_H3x2z_G2x2z = I_ESP_I3x3z_F2xz+ABZ*I_ESP_H3x2z_F2xz;
    Double I_ESP_H2x3y_G2x2z = I_ESP_I2x3yz_F2xz+ABZ*I_ESP_H2x3y_F2xz;
    Double I_ESP_H2x2yz_G2x2z = I_ESP_I2x2y2z_F2xz+ABZ*I_ESP_H2x2yz_F2xz;
    Double I_ESP_H2xy2z_G2x2z = I_ESP_I2xy3z_F2xz+ABZ*I_ESP_H2xy2z_F2xz;
    Double I_ESP_H2x3z_G2x2z = I_ESP_I2x4z_F2xz+ABZ*I_ESP_H2x3z_F2xz;
    Double I_ESP_Hx4y_G2x2z = I_ESP_Ix4yz_F2xz+ABZ*I_ESP_Hx4y_F2xz;
    Double I_ESP_Hx3yz_G2x2z = I_ESP_Ix3y2z_F2xz+ABZ*I_ESP_Hx3yz_F2xz;
    Double I_ESP_Hx2y2z_G2x2z = I_ESP_Ix2y3z_F2xz+ABZ*I_ESP_Hx2y2z_F2xz;
    Double I_ESP_Hxy3z_G2x2z = I_ESP_Ixy4z_F2xz+ABZ*I_ESP_Hxy3z_F2xz;
    Double I_ESP_Hx4z_G2x2z = I_ESP_Ix5z_F2xz+ABZ*I_ESP_Hx4z_F2xz;
    Double I_ESP_H5y_G2x2z = I_ESP_I5yz_F2xz+ABZ*I_ESP_H5y_F2xz;
    Double I_ESP_H4yz_G2x2z = I_ESP_I4y2z_F2xz+ABZ*I_ESP_H4yz_F2xz;
    Double I_ESP_H3y2z_G2x2z = I_ESP_I3y3z_F2xz+ABZ*I_ESP_H3y2z_F2xz;
    Double I_ESP_H2y3z_G2x2z = I_ESP_I2y4z_F2xz+ABZ*I_ESP_H2y3z_F2xz;
    Double I_ESP_Hy4z_G2x2z = I_ESP_Iy5z_F2xz+ABZ*I_ESP_Hy4z_F2xz;
    Double I_ESP_H5z_G2x2z = I_ESP_I6z_F2xz+ABZ*I_ESP_H5z_F2xz;
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
    Double I_ESP_H5y_Gx3y = I_ESP_Ix5y_F3y+ABX*I_ESP_H5y_F3y;
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
    Double I_ESP_H5z_Gx3z = I_ESP_Ix5z_F3z+ABX*I_ESP_H5z_F3z;
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
    Double I_ESP_H5x_G3yz = I_ESP_I5xz_F3y+ABZ*I_ESP_H5x_F3y;
    Double I_ESP_H4xy_G3yz = I_ESP_I4xyz_F3y+ABZ*I_ESP_H4xy_F3y;
    Double I_ESP_H4xz_G3yz = I_ESP_I4x2z_F3y+ABZ*I_ESP_H4xz_F3y;
    Double I_ESP_H3x2y_G3yz = I_ESP_I3x2yz_F3y+ABZ*I_ESP_H3x2y_F3y;
    Double I_ESP_H3xyz_G3yz = I_ESP_I3xy2z_F3y+ABZ*I_ESP_H3xyz_F3y;
    Double I_ESP_H3x2z_G3yz = I_ESP_I3x3z_F3y+ABZ*I_ESP_H3x2z_F3y;
    Double I_ESP_H2x3y_G3yz = I_ESP_I2x3yz_F3y+ABZ*I_ESP_H2x3y_F3y;
    Double I_ESP_H2x2yz_G3yz = I_ESP_I2x2y2z_F3y+ABZ*I_ESP_H2x2yz_F3y;
    Double I_ESP_H2xy2z_G3yz = I_ESP_I2xy3z_F3y+ABZ*I_ESP_H2xy2z_F3y;
    Double I_ESP_H2x3z_G3yz = I_ESP_I2x4z_F3y+ABZ*I_ESP_H2x3z_F3y;
    Double I_ESP_Hx4y_G3yz = I_ESP_Ix4yz_F3y+ABZ*I_ESP_Hx4y_F3y;
    Double I_ESP_Hx3yz_G3yz = I_ESP_Ix3y2z_F3y+ABZ*I_ESP_Hx3yz_F3y;
    Double I_ESP_Hx2y2z_G3yz = I_ESP_Ix2y3z_F3y+ABZ*I_ESP_Hx2y2z_F3y;
    Double I_ESP_Hxy3z_G3yz = I_ESP_Ixy4z_F3y+ABZ*I_ESP_Hxy3z_F3y;
    Double I_ESP_Hx4z_G3yz = I_ESP_Ix5z_F3y+ABZ*I_ESP_Hx4z_F3y;
    Double I_ESP_H5y_G3yz = I_ESP_I5yz_F3y+ABZ*I_ESP_H5y_F3y;
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
    Double I_ESP_H5y_G2y2z = I_ESP_I5yz_F2yz+ABZ*I_ESP_H5y_F2yz;
    Double I_ESP_H4yz_G2y2z = I_ESP_I4y2z_F2yz+ABZ*I_ESP_H4yz_F2yz;
    Double I_ESP_H3y2z_G2y2z = I_ESP_I3y3z_F2yz+ABZ*I_ESP_H3y2z_F2yz;
    Double I_ESP_H2y3z_G2y2z = I_ESP_I2y4z_F2yz+ABZ*I_ESP_H2y3z_F2yz;
    Double I_ESP_Hy4z_G2y2z = I_ESP_Iy5z_F2yz+ABZ*I_ESP_Hy4z_F2yz;
    Double I_ESP_H5z_G2y2z = I_ESP_I6z_F2yz+ABZ*I_ESP_H5z_F2yz;
    Double I_ESP_H5x_Gy3z = I_ESP_I5xy_F3z+ABY*I_ESP_H5x_F3z;
    Double I_ESP_H4xy_Gy3z = I_ESP_I4x2y_F3z+ABY*I_ESP_H4xy_F3z;
    Double I_ESP_H4xz_Gy3z = I_ESP_I4xyz_F3z+ABY*I_ESP_H4xz_F3z;
    Double I_ESP_H3x2y_Gy3z = I_ESP_I3x3y_F3z+ABY*I_ESP_H3x2y_F3z;
    Double I_ESP_H3xyz_Gy3z = I_ESP_I3x2yz_F3z+ABY*I_ESP_H3xyz_F3z;
    Double I_ESP_H3x2z_Gy3z = I_ESP_I3xy2z_F3z+ABY*I_ESP_H3x2z_F3z;
    Double I_ESP_H2x3y_Gy3z = I_ESP_I2x4y_F3z+ABY*I_ESP_H2x3y_F3z;
    Double I_ESP_H2x2yz_Gy3z = I_ESP_I2x3yz_F3z+ABY*I_ESP_H2x2yz_F3z;
    Double I_ESP_H2xy2z_Gy3z = I_ESP_I2x2y2z_F3z+ABY*I_ESP_H2xy2z_F3z;
    Double I_ESP_H2x3z_Gy3z = I_ESP_I2xy3z_F3z+ABY*I_ESP_H2x3z_F3z;
    Double I_ESP_Hx4y_Gy3z = I_ESP_Ix5y_F3z+ABY*I_ESP_Hx4y_F3z;
    Double I_ESP_Hx3yz_Gy3z = I_ESP_Ix4yz_F3z+ABY*I_ESP_Hx3yz_F3z;
    Double I_ESP_Hx2y2z_Gy3z = I_ESP_Ix3y2z_F3z+ABY*I_ESP_Hx2y2z_F3z;
    Double I_ESP_Hxy3z_Gy3z = I_ESP_Ix2y3z_F3z+ABY*I_ESP_Hxy3z_F3z;
    Double I_ESP_Hx4z_Gy3z = I_ESP_Ixy4z_F3z+ABY*I_ESP_Hx4z_F3z;
    Double I_ESP_H5y_Gy3z = I_ESP_I6y_F3z+ABY*I_ESP_H5y_F3z;
    Double I_ESP_H4yz_Gy3z = I_ESP_I5yz_F3z+ABY*I_ESP_H4yz_F3z;
    Double I_ESP_H3y2z_Gy3z = I_ESP_I4y2z_F3z+ABY*I_ESP_H3y2z_F3z;
    Double I_ESP_H2y3z_Gy3z = I_ESP_I3y3z_F3z+ABY*I_ESP_H2y3z_F3z;
    Double I_ESP_Hy4z_Gy3z = I_ESP_I2y4z_F3z+ABY*I_ESP_Hy4z_F3z;
    Double I_ESP_H5z_Gy3z = I_ESP_Iy5z_F3z+ABY*I_ESP_H5z_F3z;
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
     * shell quartet name: SQ_ESP_M_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 44 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_N_S
     * RHS shell quartet name: SQ_ESP_M_S
     ************************************************************/
    Double I_ESP_M9x_Px = I_ESP_N10x_S+ABX*I_ESP_M9x_S;
    Double I_ESP_M8xy_Px = I_ESP_N9xy_S+ABX*I_ESP_M8xy_S;
    Double I_ESP_M8xz_Px = I_ESP_N9xz_S+ABX*I_ESP_M8xz_S;
    Double I_ESP_M7x2y_Px = I_ESP_N8x2y_S+ABX*I_ESP_M7x2y_S;
    Double I_ESP_M7xyz_Px = I_ESP_N8xyz_S+ABX*I_ESP_M7xyz_S;
    Double I_ESP_M7x2z_Px = I_ESP_N8x2z_S+ABX*I_ESP_M7x2z_S;
    Double I_ESP_M6x3y_Px = I_ESP_N7x3y_S+ABX*I_ESP_M6x3y_S;
    Double I_ESP_M6x2yz_Px = I_ESP_N7x2yz_S+ABX*I_ESP_M6x2yz_S;
    Double I_ESP_M6xy2z_Px = I_ESP_N7xy2z_S+ABX*I_ESP_M6xy2z_S;
    Double I_ESP_M6x3z_Px = I_ESP_N7x3z_S+ABX*I_ESP_M6x3z_S;
    Double I_ESP_M5x4y_Px = I_ESP_N6x4y_S+ABX*I_ESP_M5x4y_S;
    Double I_ESP_M5x3yz_Px = I_ESP_N6x3yz_S+ABX*I_ESP_M5x3yz_S;
    Double I_ESP_M5x2y2z_Px = I_ESP_N6x2y2z_S+ABX*I_ESP_M5x2y2z_S;
    Double I_ESP_M5xy3z_Px = I_ESP_N6xy3z_S+ABX*I_ESP_M5xy3z_S;
    Double I_ESP_M5x4z_Px = I_ESP_N6x4z_S+ABX*I_ESP_M5x4z_S;
    Double I_ESP_M4x5y_Px = I_ESP_N5x5y_S+ABX*I_ESP_M4x5y_S;
    Double I_ESP_M4x4yz_Px = I_ESP_N5x4yz_S+ABX*I_ESP_M4x4yz_S;
    Double I_ESP_M4x3y2z_Px = I_ESP_N5x3y2z_S+ABX*I_ESP_M4x3y2z_S;
    Double I_ESP_M4x2y3z_Px = I_ESP_N5x2y3z_S+ABX*I_ESP_M4x2y3z_S;
    Double I_ESP_M4xy4z_Px = I_ESP_N5xy4z_S+ABX*I_ESP_M4xy4z_S;
    Double I_ESP_M4x5z_Px = I_ESP_N5x5z_S+ABX*I_ESP_M4x5z_S;
    Double I_ESP_M3x6y_Px = I_ESP_N4x6y_S+ABX*I_ESP_M3x6y_S;
    Double I_ESP_M3x5yz_Px = I_ESP_N4x5yz_S+ABX*I_ESP_M3x5yz_S;
    Double I_ESP_M3x4y2z_Px = I_ESP_N4x4y2z_S+ABX*I_ESP_M3x4y2z_S;
    Double I_ESP_M3x3y3z_Px = I_ESP_N4x3y3z_S+ABX*I_ESP_M3x3y3z_S;
    Double I_ESP_M3x2y4z_Px = I_ESP_N4x2y4z_S+ABX*I_ESP_M3x2y4z_S;
    Double I_ESP_M3xy5z_Px = I_ESP_N4xy5z_S+ABX*I_ESP_M3xy5z_S;
    Double I_ESP_M3x6z_Px = I_ESP_N4x6z_S+ABX*I_ESP_M3x6z_S;
    Double I_ESP_M2x7y_Px = I_ESP_N3x7y_S+ABX*I_ESP_M2x7y_S;
    Double I_ESP_M2x6yz_Px = I_ESP_N3x6yz_S+ABX*I_ESP_M2x6yz_S;
    Double I_ESP_M2x5y2z_Px = I_ESP_N3x5y2z_S+ABX*I_ESP_M2x5y2z_S;
    Double I_ESP_M2x4y3z_Px = I_ESP_N3x4y3z_S+ABX*I_ESP_M2x4y3z_S;
    Double I_ESP_M2x3y4z_Px = I_ESP_N3x3y4z_S+ABX*I_ESP_M2x3y4z_S;
    Double I_ESP_M2x2y5z_Px = I_ESP_N3x2y5z_S+ABX*I_ESP_M2x2y5z_S;
    Double I_ESP_M2xy6z_Px = I_ESP_N3xy6z_S+ABX*I_ESP_M2xy6z_S;
    Double I_ESP_M2x7z_Px = I_ESP_N3x7z_S+ABX*I_ESP_M2x7z_S;
    Double I_ESP_Mx7yz_Px = I_ESP_N2x7yz_S+ABX*I_ESP_Mx7yz_S;
    Double I_ESP_Mx6y2z_Px = I_ESP_N2x6y2z_S+ABX*I_ESP_Mx6y2z_S;
    Double I_ESP_Mx5y3z_Px = I_ESP_N2x5y3z_S+ABX*I_ESP_Mx5y3z_S;
    Double I_ESP_Mx4y4z_Px = I_ESP_N2x4y4z_S+ABX*I_ESP_Mx4y4z_S;
    Double I_ESP_Mx3y5z_Px = I_ESP_N2x3y5z_S+ABX*I_ESP_Mx3y5z_S;
    Double I_ESP_Mx2y6z_Px = I_ESP_N2x2y6z_S+ABX*I_ESP_Mx2y6z_S;
    Double I_ESP_Mxy7z_Px = I_ESP_N2xy7z_S+ABX*I_ESP_Mxy7z_S;
    Double I_ESP_M7x2y_Py = I_ESP_N7x3y_S+ABY*I_ESP_M7x2y_S;
    Double I_ESP_M6x3y_Py = I_ESP_N6x4y_S+ABY*I_ESP_M6x3y_S;
    Double I_ESP_M6x2yz_Py = I_ESP_N6x3yz_S+ABY*I_ESP_M6x2yz_S;
    Double I_ESP_M6xy2z_Py = I_ESP_N6x2y2z_S+ABY*I_ESP_M6xy2z_S;
    Double I_ESP_M5x4y_Py = I_ESP_N5x5y_S+ABY*I_ESP_M5x4y_S;
    Double I_ESP_M5x3yz_Py = I_ESP_N5x4yz_S+ABY*I_ESP_M5x3yz_S;
    Double I_ESP_M5x2y2z_Py = I_ESP_N5x3y2z_S+ABY*I_ESP_M5x2y2z_S;
    Double I_ESP_M5xy3z_Py = I_ESP_N5x2y3z_S+ABY*I_ESP_M5xy3z_S;
    Double I_ESP_M4x5y_Py = I_ESP_N4x6y_S+ABY*I_ESP_M4x5y_S;
    Double I_ESP_M4x4yz_Py = I_ESP_N4x5yz_S+ABY*I_ESP_M4x4yz_S;
    Double I_ESP_M4x3y2z_Py = I_ESP_N4x4y2z_S+ABY*I_ESP_M4x3y2z_S;
    Double I_ESP_M4x2y3z_Py = I_ESP_N4x3y3z_S+ABY*I_ESP_M4x2y3z_S;
    Double I_ESP_M4xy4z_Py = I_ESP_N4x2y4z_S+ABY*I_ESP_M4xy4z_S;
    Double I_ESP_M3x6y_Py = I_ESP_N3x7y_S+ABY*I_ESP_M3x6y_S;
    Double I_ESP_M3x5yz_Py = I_ESP_N3x6yz_S+ABY*I_ESP_M3x5yz_S;
    Double I_ESP_M3x4y2z_Py = I_ESP_N3x5y2z_S+ABY*I_ESP_M3x4y2z_S;
    Double I_ESP_M3x3y3z_Py = I_ESP_N3x4y3z_S+ABY*I_ESP_M3x3y3z_S;
    Double I_ESP_M3x2y4z_Py = I_ESP_N3x3y4z_S+ABY*I_ESP_M3x2y4z_S;
    Double I_ESP_M3xy5z_Py = I_ESP_N3x2y5z_S+ABY*I_ESP_M3xy5z_S;
    Double I_ESP_M2x7y_Py = I_ESP_N2x8y_S+ABY*I_ESP_M2x7y_S;
    Double I_ESP_M2x6yz_Py = I_ESP_N2x7yz_S+ABY*I_ESP_M2x6yz_S;
    Double I_ESP_M2x5y2z_Py = I_ESP_N2x6y2z_S+ABY*I_ESP_M2x5y2z_S;
    Double I_ESP_M2x4y3z_Py = I_ESP_N2x5y3z_S+ABY*I_ESP_M2x4y3z_S;
    Double I_ESP_M2x3y4z_Py = I_ESP_N2x4y4z_S+ABY*I_ESP_M2x3y4z_S;
    Double I_ESP_M2x2y5z_Py = I_ESP_N2x3y5z_S+ABY*I_ESP_M2x2y5z_S;
    Double I_ESP_M2xy6z_Py = I_ESP_N2x2y6z_S+ABY*I_ESP_M2xy6z_S;
    Double I_ESP_Mx8y_Py = I_ESP_Nx9y_S+ABY*I_ESP_Mx8y_S;
    Double I_ESP_Mx7yz_Py = I_ESP_Nx8yz_S+ABY*I_ESP_Mx7yz_S;
    Double I_ESP_Mx6y2z_Py = I_ESP_Nx7y2z_S+ABY*I_ESP_Mx6y2z_S;
    Double I_ESP_Mx5y3z_Py = I_ESP_Nx6y3z_S+ABY*I_ESP_Mx5y3z_S;
    Double I_ESP_Mx4y4z_Py = I_ESP_Nx5y4z_S+ABY*I_ESP_Mx4y4z_S;
    Double I_ESP_Mx3y5z_Py = I_ESP_Nx4y5z_S+ABY*I_ESP_Mx3y5z_S;
    Double I_ESP_Mx2y6z_Py = I_ESP_Nx3y6z_S+ABY*I_ESP_Mx2y6z_S;
    Double I_ESP_Mxy7z_Py = I_ESP_Nx2y7z_S+ABY*I_ESP_Mxy7z_S;
    Double I_ESP_M9y_Py = I_ESP_N10y_S+ABY*I_ESP_M9y_S;
    Double I_ESP_M8yz_Py = I_ESP_N9yz_S+ABY*I_ESP_M8yz_S;
    Double I_ESP_M7y2z_Py = I_ESP_N8y2z_S+ABY*I_ESP_M7y2z_S;
    Double I_ESP_M6y3z_Py = I_ESP_N7y3z_S+ABY*I_ESP_M6y3z_S;
    Double I_ESP_M5y4z_Py = I_ESP_N6y4z_S+ABY*I_ESP_M5y4z_S;
    Double I_ESP_M4y5z_Py = I_ESP_N5y5z_S+ABY*I_ESP_M4y5z_S;
    Double I_ESP_M3y6z_Py = I_ESP_N4y6z_S+ABY*I_ESP_M3y6z_S;
    Double I_ESP_M2y7z_Py = I_ESP_N3y7z_S+ABY*I_ESP_M2y7z_S;
    Double I_ESP_M7x2z_Pz = I_ESP_N7x3z_S+ABZ*I_ESP_M7x2z_S;
    Double I_ESP_M6xy2z_Pz = I_ESP_N6xy3z_S+ABZ*I_ESP_M6xy2z_S;
    Double I_ESP_M6x3z_Pz = I_ESP_N6x4z_S+ABZ*I_ESP_M6x3z_S;
    Double I_ESP_M5x2y2z_Pz = I_ESP_N5x2y3z_S+ABZ*I_ESP_M5x2y2z_S;
    Double I_ESP_M5xy3z_Pz = I_ESP_N5xy4z_S+ABZ*I_ESP_M5xy3z_S;
    Double I_ESP_M5x4z_Pz = I_ESP_N5x5z_S+ABZ*I_ESP_M5x4z_S;
    Double I_ESP_M4x3y2z_Pz = I_ESP_N4x3y3z_S+ABZ*I_ESP_M4x3y2z_S;
    Double I_ESP_M4x2y3z_Pz = I_ESP_N4x2y4z_S+ABZ*I_ESP_M4x2y3z_S;
    Double I_ESP_M4xy4z_Pz = I_ESP_N4xy5z_S+ABZ*I_ESP_M4xy4z_S;
    Double I_ESP_M4x5z_Pz = I_ESP_N4x6z_S+ABZ*I_ESP_M4x5z_S;
    Double I_ESP_M3x4y2z_Pz = I_ESP_N3x4y3z_S+ABZ*I_ESP_M3x4y2z_S;
    Double I_ESP_M3x3y3z_Pz = I_ESP_N3x3y4z_S+ABZ*I_ESP_M3x3y3z_S;
    Double I_ESP_M3x2y4z_Pz = I_ESP_N3x2y5z_S+ABZ*I_ESP_M3x2y4z_S;
    Double I_ESP_M3xy5z_Pz = I_ESP_N3xy6z_S+ABZ*I_ESP_M3xy5z_S;
    Double I_ESP_M3x6z_Pz = I_ESP_N3x7z_S+ABZ*I_ESP_M3x6z_S;
    Double I_ESP_M2x5y2z_Pz = I_ESP_N2x5y3z_S+ABZ*I_ESP_M2x5y2z_S;
    Double I_ESP_M2x4y3z_Pz = I_ESP_N2x4y4z_S+ABZ*I_ESP_M2x4y3z_S;
    Double I_ESP_M2x3y4z_Pz = I_ESP_N2x3y5z_S+ABZ*I_ESP_M2x3y4z_S;
    Double I_ESP_M2x2y5z_Pz = I_ESP_N2x2y6z_S+ABZ*I_ESP_M2x2y5z_S;
    Double I_ESP_M2xy6z_Pz = I_ESP_N2xy7z_S+ABZ*I_ESP_M2xy6z_S;
    Double I_ESP_M2x7z_Pz = I_ESP_N2x8z_S+ABZ*I_ESP_M2x7z_S;
    Double I_ESP_Mx6y2z_Pz = I_ESP_Nx6y3z_S+ABZ*I_ESP_Mx6y2z_S;
    Double I_ESP_Mx5y3z_Pz = I_ESP_Nx5y4z_S+ABZ*I_ESP_Mx5y3z_S;
    Double I_ESP_Mx4y4z_Pz = I_ESP_Nx4y5z_S+ABZ*I_ESP_Mx4y4z_S;
    Double I_ESP_Mx3y5z_Pz = I_ESP_Nx3y6z_S+ABZ*I_ESP_Mx3y5z_S;
    Double I_ESP_Mx2y6z_Pz = I_ESP_Nx2y7z_S+ABZ*I_ESP_Mx2y6z_S;
    Double I_ESP_Mxy7z_Pz = I_ESP_Nxy8z_S+ABZ*I_ESP_Mxy7z_S;
    Double I_ESP_Mx8z_Pz = I_ESP_Nx9z_S+ABZ*I_ESP_Mx8z_S;
    Double I_ESP_M7y2z_Pz = I_ESP_N7y3z_S+ABZ*I_ESP_M7y2z_S;
    Double I_ESP_M6y3z_Pz = I_ESP_N6y4z_S+ABZ*I_ESP_M6y3z_S;
    Double I_ESP_M5y4z_Pz = I_ESP_N5y5z_S+ABZ*I_ESP_M5y4z_S;
    Double I_ESP_M4y5z_Pz = I_ESP_N4y6z_S+ABZ*I_ESP_M4y5z_S;
    Double I_ESP_M3y6z_Pz = I_ESP_N3y7z_S+ABZ*I_ESP_M3y6z_S;
    Double I_ESP_M2y7z_Pz = I_ESP_N2y8z_S+ABZ*I_ESP_M2y7z_S;
    Double I_ESP_My8z_Pz = I_ESP_Ny9z_S+ABZ*I_ESP_My8z_S;
    Double I_ESP_M9z_Pz = I_ESP_N10z_S+ABZ*I_ESP_M9z_S;

    /************************************************************
     * shell quartet name: SQ_ESP_L_D
     * expanding position: BRA2
     * code section is: HRR
     * totally 149 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_M_P
     * RHS shell quartet name: SQ_ESP_L_P
     ************************************************************/
    Double I_ESP_L8x_D2x = I_ESP_M9x_Px+ABX*I_ESP_L8x_Px;
    Double I_ESP_L7xy_D2x = I_ESP_M8xy_Px+ABX*I_ESP_L7xy_Px;
    Double I_ESP_L7xz_D2x = I_ESP_M8xz_Px+ABX*I_ESP_L7xz_Px;
    Double I_ESP_L6x2y_D2x = I_ESP_M7x2y_Px+ABX*I_ESP_L6x2y_Px;
    Double I_ESP_L6xyz_D2x = I_ESP_M7xyz_Px+ABX*I_ESP_L6xyz_Px;
    Double I_ESP_L6x2z_D2x = I_ESP_M7x2z_Px+ABX*I_ESP_L6x2z_Px;
    Double I_ESP_L5x3y_D2x = I_ESP_M6x3y_Px+ABX*I_ESP_L5x3y_Px;
    Double I_ESP_L5x2yz_D2x = I_ESP_M6x2yz_Px+ABX*I_ESP_L5x2yz_Px;
    Double I_ESP_L5xy2z_D2x = I_ESP_M6xy2z_Px+ABX*I_ESP_L5xy2z_Px;
    Double I_ESP_L5x3z_D2x = I_ESP_M6x3z_Px+ABX*I_ESP_L5x3z_Px;
    Double I_ESP_L4x4y_D2x = I_ESP_M5x4y_Px+ABX*I_ESP_L4x4y_Px;
    Double I_ESP_L4x3yz_D2x = I_ESP_M5x3yz_Px+ABX*I_ESP_L4x3yz_Px;
    Double I_ESP_L4x2y2z_D2x = I_ESP_M5x2y2z_Px+ABX*I_ESP_L4x2y2z_Px;
    Double I_ESP_L4xy3z_D2x = I_ESP_M5xy3z_Px+ABX*I_ESP_L4xy3z_Px;
    Double I_ESP_L4x4z_D2x = I_ESP_M5x4z_Px+ABX*I_ESP_L4x4z_Px;
    Double I_ESP_L3x5y_D2x = I_ESP_M4x5y_Px+ABX*I_ESP_L3x5y_Px;
    Double I_ESP_L3x4yz_D2x = I_ESP_M4x4yz_Px+ABX*I_ESP_L3x4yz_Px;
    Double I_ESP_L3x3y2z_D2x = I_ESP_M4x3y2z_Px+ABX*I_ESP_L3x3y2z_Px;
    Double I_ESP_L3x2y3z_D2x = I_ESP_M4x2y3z_Px+ABX*I_ESP_L3x2y3z_Px;
    Double I_ESP_L3xy4z_D2x = I_ESP_M4xy4z_Px+ABX*I_ESP_L3xy4z_Px;
    Double I_ESP_L3x5z_D2x = I_ESP_M4x5z_Px+ABX*I_ESP_L3x5z_Px;
    Double I_ESP_L2x6y_D2x = I_ESP_M3x6y_Px+ABX*I_ESP_L2x6y_Px;
    Double I_ESP_L2x5yz_D2x = I_ESP_M3x5yz_Px+ABX*I_ESP_L2x5yz_Px;
    Double I_ESP_L2x4y2z_D2x = I_ESP_M3x4y2z_Px+ABX*I_ESP_L2x4y2z_Px;
    Double I_ESP_L2x3y3z_D2x = I_ESP_M3x3y3z_Px+ABX*I_ESP_L2x3y3z_Px;
    Double I_ESP_L2x2y4z_D2x = I_ESP_M3x2y4z_Px+ABX*I_ESP_L2x2y4z_Px;
    Double I_ESP_L2xy5z_D2x = I_ESP_M3xy5z_Px+ABX*I_ESP_L2xy5z_Px;
    Double I_ESP_L2x6z_D2x = I_ESP_M3x6z_Px+ABX*I_ESP_L2x6z_Px;
    Double I_ESP_Lx7y_D2x = I_ESP_M2x7y_Px+ABX*I_ESP_Lx7y_Px;
    Double I_ESP_Lx6yz_D2x = I_ESP_M2x6yz_Px+ABX*I_ESP_Lx6yz_Px;
    Double I_ESP_Lx5y2z_D2x = I_ESP_M2x5y2z_Px+ABX*I_ESP_Lx5y2z_Px;
    Double I_ESP_Lx4y3z_D2x = I_ESP_M2x4y3z_Px+ABX*I_ESP_Lx4y3z_Px;
    Double I_ESP_Lx3y4z_D2x = I_ESP_M2x3y4z_Px+ABX*I_ESP_Lx3y4z_Px;
    Double I_ESP_Lx2y5z_D2x = I_ESP_M2x2y5z_Px+ABX*I_ESP_Lx2y5z_Px;
    Double I_ESP_Lxy6z_D2x = I_ESP_M2xy6z_Px+ABX*I_ESP_Lxy6z_Px;
    Double I_ESP_Lx7z_D2x = I_ESP_M2x7z_Px+ABX*I_ESP_Lx7z_Px;
    Double I_ESP_L7yz_D2x = I_ESP_Mx7yz_Px+ABX*I_ESP_L7yz_Px;
    Double I_ESP_L6y2z_D2x = I_ESP_Mx6y2z_Px+ABX*I_ESP_L6y2z_Px;
    Double I_ESP_L5y3z_D2x = I_ESP_Mx5y3z_Px+ABX*I_ESP_L5y3z_Px;
    Double I_ESP_L4y4z_D2x = I_ESP_Mx4y4z_Px+ABX*I_ESP_L4y4z_Px;
    Double I_ESP_L3y5z_D2x = I_ESP_Mx3y5z_Px+ABX*I_ESP_L3y5z_Px;
    Double I_ESP_L2y6z_D2x = I_ESP_Mx2y6z_Px+ABX*I_ESP_L2y6z_Px;
    Double I_ESP_Ly7z_D2x = I_ESP_Mxy7z_Px+ABX*I_ESP_Ly7z_Px;
    Double I_ESP_L7xy_D2y = I_ESP_M7x2y_Py+ABY*I_ESP_L7xy_Py;
    Double I_ESP_L6x2y_D2y = I_ESP_M6x3y_Py+ABY*I_ESP_L6x2y_Py;
    Double I_ESP_L6xyz_D2y = I_ESP_M6x2yz_Py+ABY*I_ESP_L6xyz_Py;
    Double I_ESP_L6x2z_D2y = I_ESP_M6xy2z_Py+ABY*I_ESP_L6x2z_Py;
    Double I_ESP_L5x3y_D2y = I_ESP_M5x4y_Py+ABY*I_ESP_L5x3y_Py;
    Double I_ESP_L5x2yz_D2y = I_ESP_M5x3yz_Py+ABY*I_ESP_L5x2yz_Py;
    Double I_ESP_L5xy2z_D2y = I_ESP_M5x2y2z_Py+ABY*I_ESP_L5xy2z_Py;
    Double I_ESP_L5x3z_D2y = I_ESP_M5xy3z_Py+ABY*I_ESP_L5x3z_Py;
    Double I_ESP_L4x4y_D2y = I_ESP_M4x5y_Py+ABY*I_ESP_L4x4y_Py;
    Double I_ESP_L4x3yz_D2y = I_ESP_M4x4yz_Py+ABY*I_ESP_L4x3yz_Py;
    Double I_ESP_L4x2y2z_D2y = I_ESP_M4x3y2z_Py+ABY*I_ESP_L4x2y2z_Py;
    Double I_ESP_L4xy3z_D2y = I_ESP_M4x2y3z_Py+ABY*I_ESP_L4xy3z_Py;
    Double I_ESP_L4x4z_D2y = I_ESP_M4xy4z_Py+ABY*I_ESP_L4x4z_Py;
    Double I_ESP_L3x5y_D2y = I_ESP_M3x6y_Py+ABY*I_ESP_L3x5y_Py;
    Double I_ESP_L3x4yz_D2y = I_ESP_M3x5yz_Py+ABY*I_ESP_L3x4yz_Py;
    Double I_ESP_L3x3y2z_D2y = I_ESP_M3x4y2z_Py+ABY*I_ESP_L3x3y2z_Py;
    Double I_ESP_L3x2y3z_D2y = I_ESP_M3x3y3z_Py+ABY*I_ESP_L3x2y3z_Py;
    Double I_ESP_L3xy4z_D2y = I_ESP_M3x2y4z_Py+ABY*I_ESP_L3xy4z_Py;
    Double I_ESP_L3x5z_D2y = I_ESP_M3xy5z_Py+ABY*I_ESP_L3x5z_Py;
    Double I_ESP_L2x6y_D2y = I_ESP_M2x7y_Py+ABY*I_ESP_L2x6y_Py;
    Double I_ESP_L2x5yz_D2y = I_ESP_M2x6yz_Py+ABY*I_ESP_L2x5yz_Py;
    Double I_ESP_L2x4y2z_D2y = I_ESP_M2x5y2z_Py+ABY*I_ESP_L2x4y2z_Py;
    Double I_ESP_L2x3y3z_D2y = I_ESP_M2x4y3z_Py+ABY*I_ESP_L2x3y3z_Py;
    Double I_ESP_L2x2y4z_D2y = I_ESP_M2x3y4z_Py+ABY*I_ESP_L2x2y4z_Py;
    Double I_ESP_L2xy5z_D2y = I_ESP_M2x2y5z_Py+ABY*I_ESP_L2xy5z_Py;
    Double I_ESP_L2x6z_D2y = I_ESP_M2xy6z_Py+ABY*I_ESP_L2x6z_Py;
    Double I_ESP_Lx7y_D2y = I_ESP_Mx8y_Py+ABY*I_ESP_Lx7y_Py;
    Double I_ESP_Lx6yz_D2y = I_ESP_Mx7yz_Py+ABY*I_ESP_Lx6yz_Py;
    Double I_ESP_Lx5y2z_D2y = I_ESP_Mx6y2z_Py+ABY*I_ESP_Lx5y2z_Py;
    Double I_ESP_Lx4y3z_D2y = I_ESP_Mx5y3z_Py+ABY*I_ESP_Lx4y3z_Py;
    Double I_ESP_Lx3y4z_D2y = I_ESP_Mx4y4z_Py+ABY*I_ESP_Lx3y4z_Py;
    Double I_ESP_Lx2y5z_D2y = I_ESP_Mx3y5z_Py+ABY*I_ESP_Lx2y5z_Py;
    Double I_ESP_Lxy6z_D2y = I_ESP_Mx2y6z_Py+ABY*I_ESP_Lxy6z_Py;
    Double I_ESP_Lx7z_D2y = I_ESP_Mxy7z_Py+ABY*I_ESP_Lx7z_Py;
    Double I_ESP_L8y_D2y = I_ESP_M9y_Py+ABY*I_ESP_L8y_Py;
    Double I_ESP_L7yz_D2y = I_ESP_M8yz_Py+ABY*I_ESP_L7yz_Py;
    Double I_ESP_L6y2z_D2y = I_ESP_M7y2z_Py+ABY*I_ESP_L6y2z_Py;
    Double I_ESP_L5y3z_D2y = I_ESP_M6y3z_Py+ABY*I_ESP_L5y3z_Py;
    Double I_ESP_L4y4z_D2y = I_ESP_M5y4z_Py+ABY*I_ESP_L4y4z_Py;
    Double I_ESP_L3y5z_D2y = I_ESP_M4y5z_Py+ABY*I_ESP_L3y5z_Py;
    Double I_ESP_L2y6z_D2y = I_ESP_M3y6z_Py+ABY*I_ESP_L2y6z_Py;
    Double I_ESP_Ly7z_D2y = I_ESP_M2y7z_Py+ABY*I_ESP_Ly7z_Py;
    Double I_ESP_L7xz_D2z = I_ESP_M7x2z_Pz+ABZ*I_ESP_L7xz_Pz;
    Double I_ESP_L6xyz_D2z = I_ESP_M6xy2z_Pz+ABZ*I_ESP_L6xyz_Pz;
    Double I_ESP_L6x2z_D2z = I_ESP_M6x3z_Pz+ABZ*I_ESP_L6x2z_Pz;
    Double I_ESP_L5x2yz_D2z = I_ESP_M5x2y2z_Pz+ABZ*I_ESP_L5x2yz_Pz;
    Double I_ESP_L5xy2z_D2z = I_ESP_M5xy3z_Pz+ABZ*I_ESP_L5xy2z_Pz;
    Double I_ESP_L5x3z_D2z = I_ESP_M5x4z_Pz+ABZ*I_ESP_L5x3z_Pz;
    Double I_ESP_L4x3yz_D2z = I_ESP_M4x3y2z_Pz+ABZ*I_ESP_L4x3yz_Pz;
    Double I_ESP_L4x2y2z_D2z = I_ESP_M4x2y3z_Pz+ABZ*I_ESP_L4x2y2z_Pz;
    Double I_ESP_L4xy3z_D2z = I_ESP_M4xy4z_Pz+ABZ*I_ESP_L4xy3z_Pz;
    Double I_ESP_L4x4z_D2z = I_ESP_M4x5z_Pz+ABZ*I_ESP_L4x4z_Pz;
    Double I_ESP_L3x4yz_D2z = I_ESP_M3x4y2z_Pz+ABZ*I_ESP_L3x4yz_Pz;
    Double I_ESP_L3x3y2z_D2z = I_ESP_M3x3y3z_Pz+ABZ*I_ESP_L3x3y2z_Pz;
    Double I_ESP_L3x2y3z_D2z = I_ESP_M3x2y4z_Pz+ABZ*I_ESP_L3x2y3z_Pz;
    Double I_ESP_L3xy4z_D2z = I_ESP_M3xy5z_Pz+ABZ*I_ESP_L3xy4z_Pz;
    Double I_ESP_L3x5z_D2z = I_ESP_M3x6z_Pz+ABZ*I_ESP_L3x5z_Pz;
    Double I_ESP_L2x5yz_D2z = I_ESP_M2x5y2z_Pz+ABZ*I_ESP_L2x5yz_Pz;
    Double I_ESP_L2x4y2z_D2z = I_ESP_M2x4y3z_Pz+ABZ*I_ESP_L2x4y2z_Pz;
    Double I_ESP_L2x3y3z_D2z = I_ESP_M2x3y4z_Pz+ABZ*I_ESP_L2x3y3z_Pz;
    Double I_ESP_L2x2y4z_D2z = I_ESP_M2x2y5z_Pz+ABZ*I_ESP_L2x2y4z_Pz;
    Double I_ESP_L2xy5z_D2z = I_ESP_M2xy6z_Pz+ABZ*I_ESP_L2xy5z_Pz;
    Double I_ESP_L2x6z_D2z = I_ESP_M2x7z_Pz+ABZ*I_ESP_L2x6z_Pz;
    Double I_ESP_Lx6yz_D2z = I_ESP_Mx6y2z_Pz+ABZ*I_ESP_Lx6yz_Pz;
    Double I_ESP_Lx5y2z_D2z = I_ESP_Mx5y3z_Pz+ABZ*I_ESP_Lx5y2z_Pz;
    Double I_ESP_Lx4y3z_D2z = I_ESP_Mx4y4z_Pz+ABZ*I_ESP_Lx4y3z_Pz;
    Double I_ESP_Lx3y4z_D2z = I_ESP_Mx3y5z_Pz+ABZ*I_ESP_Lx3y4z_Pz;
    Double I_ESP_Lx2y5z_D2z = I_ESP_Mx2y6z_Pz+ABZ*I_ESP_Lx2y5z_Pz;
    Double I_ESP_Lxy6z_D2z = I_ESP_Mxy7z_Pz+ABZ*I_ESP_Lxy6z_Pz;
    Double I_ESP_Lx7z_D2z = I_ESP_Mx8z_Pz+ABZ*I_ESP_Lx7z_Pz;
    Double I_ESP_L7yz_D2z = I_ESP_M7y2z_Pz+ABZ*I_ESP_L7yz_Pz;
    Double I_ESP_L6y2z_D2z = I_ESP_M6y3z_Pz+ABZ*I_ESP_L6y2z_Pz;
    Double I_ESP_L5y3z_D2z = I_ESP_M5y4z_Pz+ABZ*I_ESP_L5y3z_Pz;
    Double I_ESP_L4y4z_D2z = I_ESP_M4y5z_Pz+ABZ*I_ESP_L4y4z_Pz;
    Double I_ESP_L3y5z_D2z = I_ESP_M3y6z_Pz+ABZ*I_ESP_L3y5z_Pz;
    Double I_ESP_L2y6z_D2z = I_ESP_M2y7z_Pz+ABZ*I_ESP_L2y6z_Pz;
    Double I_ESP_Ly7z_D2z = I_ESP_My8z_Pz+ABZ*I_ESP_Ly7z_Pz;
    Double I_ESP_L8z_D2z = I_ESP_M9z_Pz+ABZ*I_ESP_L8z_Pz;

    /************************************************************
     * shell quartet name: SQ_ESP_K_F
     * expanding position: BRA2
     * code section is: HRR
     * totally 189 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_L_D
     * RHS shell quartet name: SQ_ESP_K_D
     ************************************************************/
    Double I_ESP_K7x_F3x = I_ESP_L8x_D2x+ABX*I_ESP_K7x_D2x;
    Double I_ESP_K6xy_F3x = I_ESP_L7xy_D2x+ABX*I_ESP_K6xy_D2x;
    Double I_ESP_K6xz_F3x = I_ESP_L7xz_D2x+ABX*I_ESP_K6xz_D2x;
    Double I_ESP_K5x2y_F3x = I_ESP_L6x2y_D2x+ABX*I_ESP_K5x2y_D2x;
    Double I_ESP_K5xyz_F3x = I_ESP_L6xyz_D2x+ABX*I_ESP_K5xyz_D2x;
    Double I_ESP_K5x2z_F3x = I_ESP_L6x2z_D2x+ABX*I_ESP_K5x2z_D2x;
    Double I_ESP_K4x3y_F3x = I_ESP_L5x3y_D2x+ABX*I_ESP_K4x3y_D2x;
    Double I_ESP_K4x2yz_F3x = I_ESP_L5x2yz_D2x+ABX*I_ESP_K4x2yz_D2x;
    Double I_ESP_K4xy2z_F3x = I_ESP_L5xy2z_D2x+ABX*I_ESP_K4xy2z_D2x;
    Double I_ESP_K4x3z_F3x = I_ESP_L5x3z_D2x+ABX*I_ESP_K4x3z_D2x;
    Double I_ESP_K3x4y_F3x = I_ESP_L4x4y_D2x+ABX*I_ESP_K3x4y_D2x;
    Double I_ESP_K3x3yz_F3x = I_ESP_L4x3yz_D2x+ABX*I_ESP_K3x3yz_D2x;
    Double I_ESP_K3x2y2z_F3x = I_ESP_L4x2y2z_D2x+ABX*I_ESP_K3x2y2z_D2x;
    Double I_ESP_K3xy3z_F3x = I_ESP_L4xy3z_D2x+ABX*I_ESP_K3xy3z_D2x;
    Double I_ESP_K3x4z_F3x = I_ESP_L4x4z_D2x+ABX*I_ESP_K3x4z_D2x;
    Double I_ESP_K2x5y_F3x = I_ESP_L3x5y_D2x+ABX*I_ESP_K2x5y_D2x;
    Double I_ESP_K2x4yz_F3x = I_ESP_L3x4yz_D2x+ABX*I_ESP_K2x4yz_D2x;
    Double I_ESP_K2x3y2z_F3x = I_ESP_L3x3y2z_D2x+ABX*I_ESP_K2x3y2z_D2x;
    Double I_ESP_K2x2y3z_F3x = I_ESP_L3x2y3z_D2x+ABX*I_ESP_K2x2y3z_D2x;
    Double I_ESP_K2xy4z_F3x = I_ESP_L3xy4z_D2x+ABX*I_ESP_K2xy4z_D2x;
    Double I_ESP_K2x5z_F3x = I_ESP_L3x5z_D2x+ABX*I_ESP_K2x5z_D2x;
    Double I_ESP_Kx6y_F3x = I_ESP_L2x6y_D2x+ABX*I_ESP_Kx6y_D2x;
    Double I_ESP_Kx5yz_F3x = I_ESP_L2x5yz_D2x+ABX*I_ESP_Kx5yz_D2x;
    Double I_ESP_Kx4y2z_F3x = I_ESP_L2x4y2z_D2x+ABX*I_ESP_Kx4y2z_D2x;
    Double I_ESP_Kx3y3z_F3x = I_ESP_L2x3y3z_D2x+ABX*I_ESP_Kx3y3z_D2x;
    Double I_ESP_Kx2y4z_F3x = I_ESP_L2x2y4z_D2x+ABX*I_ESP_Kx2y4z_D2x;
    Double I_ESP_Kxy5z_F3x = I_ESP_L2xy5z_D2x+ABX*I_ESP_Kxy5z_D2x;
    Double I_ESP_Kx6z_F3x = I_ESP_L2x6z_D2x+ABX*I_ESP_Kx6z_D2x;
    Double I_ESP_K7y_F3x = I_ESP_Lx7y_D2x+ABX*I_ESP_K7y_D2x;
    Double I_ESP_K6yz_F3x = I_ESP_Lx6yz_D2x+ABX*I_ESP_K6yz_D2x;
    Double I_ESP_K5y2z_F3x = I_ESP_Lx5y2z_D2x+ABX*I_ESP_K5y2z_D2x;
    Double I_ESP_K4y3z_F3x = I_ESP_Lx4y3z_D2x+ABX*I_ESP_K4y3z_D2x;
    Double I_ESP_K3y4z_F3x = I_ESP_Lx3y4z_D2x+ABX*I_ESP_K3y4z_D2x;
    Double I_ESP_K2y5z_F3x = I_ESP_Lx2y5z_D2x+ABX*I_ESP_K2y5z_D2x;
    Double I_ESP_Ky6z_F3x = I_ESP_Lxy6z_D2x+ABX*I_ESP_Ky6z_D2x;
    Double I_ESP_K7z_F3x = I_ESP_Lx7z_D2x+ABX*I_ESP_K7z_D2x;
    Double I_ESP_K5xyz_F2xy = I_ESP_L5x2yz_D2x+ABY*I_ESP_K5xyz_D2x;
    Double I_ESP_K4x2yz_F2xy = I_ESP_L4x3yz_D2x+ABY*I_ESP_K4x2yz_D2x;
    Double I_ESP_K4xy2z_F2xy = I_ESP_L4x2y2z_D2x+ABY*I_ESP_K4xy2z_D2x;
    Double I_ESP_K3x3yz_F2xy = I_ESP_L3x4yz_D2x+ABY*I_ESP_K3x3yz_D2x;
    Double I_ESP_K3x2y2z_F2xy = I_ESP_L3x3y2z_D2x+ABY*I_ESP_K3x2y2z_D2x;
    Double I_ESP_K3xy3z_F2xy = I_ESP_L3x2y3z_D2x+ABY*I_ESP_K3xy3z_D2x;
    Double I_ESP_K2x4yz_F2xy = I_ESP_L2x5yz_D2x+ABY*I_ESP_K2x4yz_D2x;
    Double I_ESP_K2x3y2z_F2xy = I_ESP_L2x4y2z_D2x+ABY*I_ESP_K2x3y2z_D2x;
    Double I_ESP_K2x2y3z_F2xy = I_ESP_L2x3y3z_D2x+ABY*I_ESP_K2x2y3z_D2x;
    Double I_ESP_K2xy4z_F2xy = I_ESP_L2x2y4z_D2x+ABY*I_ESP_K2xy4z_D2x;
    Double I_ESP_Kx5yz_F2xy = I_ESP_Lx6yz_D2x+ABY*I_ESP_Kx5yz_D2x;
    Double I_ESP_Kx4y2z_F2xy = I_ESP_Lx5y2z_D2x+ABY*I_ESP_Kx4y2z_D2x;
    Double I_ESP_Kx3y3z_F2xy = I_ESP_Lx4y3z_D2x+ABY*I_ESP_Kx3y3z_D2x;
    Double I_ESP_Kx2y4z_F2xy = I_ESP_Lx3y4z_D2x+ABY*I_ESP_Kx2y4z_D2x;
    Double I_ESP_Kxy5z_F2xy = I_ESP_Lx2y5z_D2x+ABY*I_ESP_Kxy5z_D2x;
    Double I_ESP_K6yz_F2xy = I_ESP_L7yz_D2x+ABY*I_ESP_K6yz_D2x;
    Double I_ESP_K5y2z_F2xy = I_ESP_L6y2z_D2x+ABY*I_ESP_K5y2z_D2x;
    Double I_ESP_K4y3z_F2xy = I_ESP_L5y3z_D2x+ABY*I_ESP_K4y3z_D2x;
    Double I_ESP_K3y4z_F2xy = I_ESP_L4y4z_D2x+ABY*I_ESP_K3y4z_D2x;
    Double I_ESP_K2y5z_F2xy = I_ESP_L3y5z_D2x+ABY*I_ESP_K2y5z_D2x;
    Double I_ESP_Ky6z_F2xy = I_ESP_L2y6z_D2x+ABY*I_ESP_Ky6z_D2x;
    Double I_ESP_K5xyz_F2xz = I_ESP_L5xy2z_D2x+ABZ*I_ESP_K5xyz_D2x;
    Double I_ESP_K4x2yz_F2xz = I_ESP_L4x2y2z_D2x+ABZ*I_ESP_K4x2yz_D2x;
    Double I_ESP_K4xy2z_F2xz = I_ESP_L4xy3z_D2x+ABZ*I_ESP_K4xy2z_D2x;
    Double I_ESP_K3x3yz_F2xz = I_ESP_L3x3y2z_D2x+ABZ*I_ESP_K3x3yz_D2x;
    Double I_ESP_K3x2y2z_F2xz = I_ESP_L3x2y3z_D2x+ABZ*I_ESP_K3x2y2z_D2x;
    Double I_ESP_K3xy3z_F2xz = I_ESP_L3xy4z_D2x+ABZ*I_ESP_K3xy3z_D2x;
    Double I_ESP_K2x4yz_F2xz = I_ESP_L2x4y2z_D2x+ABZ*I_ESP_K2x4yz_D2x;
    Double I_ESP_K2x3y2z_F2xz = I_ESP_L2x3y3z_D2x+ABZ*I_ESP_K2x3y2z_D2x;
    Double I_ESP_K2x2y3z_F2xz = I_ESP_L2x2y4z_D2x+ABZ*I_ESP_K2x2y3z_D2x;
    Double I_ESP_K2xy4z_F2xz = I_ESP_L2xy5z_D2x+ABZ*I_ESP_K2xy4z_D2x;
    Double I_ESP_Kx5yz_F2xz = I_ESP_Lx5y2z_D2x+ABZ*I_ESP_Kx5yz_D2x;
    Double I_ESP_Kx4y2z_F2xz = I_ESP_Lx4y3z_D2x+ABZ*I_ESP_Kx4y2z_D2x;
    Double I_ESP_Kx3y3z_F2xz = I_ESP_Lx3y4z_D2x+ABZ*I_ESP_Kx3y3z_D2x;
    Double I_ESP_Kx2y4z_F2xz = I_ESP_Lx2y5z_D2x+ABZ*I_ESP_Kx2y4z_D2x;
    Double I_ESP_Kxy5z_F2xz = I_ESP_Lxy6z_D2x+ABZ*I_ESP_Kxy5z_D2x;
    Double I_ESP_K6yz_F2xz = I_ESP_L6y2z_D2x+ABZ*I_ESP_K6yz_D2x;
    Double I_ESP_K5y2z_F2xz = I_ESP_L5y3z_D2x+ABZ*I_ESP_K5y2z_D2x;
    Double I_ESP_K4y3z_F2xz = I_ESP_L4y4z_D2x+ABZ*I_ESP_K4y3z_D2x;
    Double I_ESP_K3y4z_F2xz = I_ESP_L3y5z_D2x+ABZ*I_ESP_K3y4z_D2x;
    Double I_ESP_K2y5z_F2xz = I_ESP_L2y6z_D2x+ABZ*I_ESP_K2y5z_D2x;
    Double I_ESP_Ky6z_F2xz = I_ESP_Ly7z_D2x+ABZ*I_ESP_Ky6z_D2x;
    Double I_ESP_K7x_F3y = I_ESP_L7xy_D2y+ABY*I_ESP_K7x_D2y;
    Double I_ESP_K6xy_F3y = I_ESP_L6x2y_D2y+ABY*I_ESP_K6xy_D2y;
    Double I_ESP_K6xz_F3y = I_ESP_L6xyz_D2y+ABY*I_ESP_K6xz_D2y;
    Double I_ESP_K5x2y_F3y = I_ESP_L5x3y_D2y+ABY*I_ESP_K5x2y_D2y;
    Double I_ESP_K5xyz_F3y = I_ESP_L5x2yz_D2y+ABY*I_ESP_K5xyz_D2y;
    Double I_ESP_K5x2z_F3y = I_ESP_L5xy2z_D2y+ABY*I_ESP_K5x2z_D2y;
    Double I_ESP_K4x3y_F3y = I_ESP_L4x4y_D2y+ABY*I_ESP_K4x3y_D2y;
    Double I_ESP_K4x2yz_F3y = I_ESP_L4x3yz_D2y+ABY*I_ESP_K4x2yz_D2y;
    Double I_ESP_K4xy2z_F3y = I_ESP_L4x2y2z_D2y+ABY*I_ESP_K4xy2z_D2y;
    Double I_ESP_K4x3z_F3y = I_ESP_L4xy3z_D2y+ABY*I_ESP_K4x3z_D2y;
    Double I_ESP_K3x4y_F3y = I_ESP_L3x5y_D2y+ABY*I_ESP_K3x4y_D2y;
    Double I_ESP_K3x3yz_F3y = I_ESP_L3x4yz_D2y+ABY*I_ESP_K3x3yz_D2y;
    Double I_ESP_K3x2y2z_F3y = I_ESP_L3x3y2z_D2y+ABY*I_ESP_K3x2y2z_D2y;
    Double I_ESP_K3xy3z_F3y = I_ESP_L3x2y3z_D2y+ABY*I_ESP_K3xy3z_D2y;
    Double I_ESP_K3x4z_F3y = I_ESP_L3xy4z_D2y+ABY*I_ESP_K3x4z_D2y;
    Double I_ESP_K2x5y_F3y = I_ESP_L2x6y_D2y+ABY*I_ESP_K2x5y_D2y;
    Double I_ESP_K2x4yz_F3y = I_ESP_L2x5yz_D2y+ABY*I_ESP_K2x4yz_D2y;
    Double I_ESP_K2x3y2z_F3y = I_ESP_L2x4y2z_D2y+ABY*I_ESP_K2x3y2z_D2y;
    Double I_ESP_K2x2y3z_F3y = I_ESP_L2x3y3z_D2y+ABY*I_ESP_K2x2y3z_D2y;
    Double I_ESP_K2xy4z_F3y = I_ESP_L2x2y4z_D2y+ABY*I_ESP_K2xy4z_D2y;
    Double I_ESP_K2x5z_F3y = I_ESP_L2xy5z_D2y+ABY*I_ESP_K2x5z_D2y;
    Double I_ESP_Kx6y_F3y = I_ESP_Lx7y_D2y+ABY*I_ESP_Kx6y_D2y;
    Double I_ESP_Kx5yz_F3y = I_ESP_Lx6yz_D2y+ABY*I_ESP_Kx5yz_D2y;
    Double I_ESP_Kx4y2z_F3y = I_ESP_Lx5y2z_D2y+ABY*I_ESP_Kx4y2z_D2y;
    Double I_ESP_Kx3y3z_F3y = I_ESP_Lx4y3z_D2y+ABY*I_ESP_Kx3y3z_D2y;
    Double I_ESP_Kx2y4z_F3y = I_ESP_Lx3y4z_D2y+ABY*I_ESP_Kx2y4z_D2y;
    Double I_ESP_Kxy5z_F3y = I_ESP_Lx2y5z_D2y+ABY*I_ESP_Kxy5z_D2y;
    Double I_ESP_Kx6z_F3y = I_ESP_Lxy6z_D2y+ABY*I_ESP_Kx6z_D2y;
    Double I_ESP_K7y_F3y = I_ESP_L8y_D2y+ABY*I_ESP_K7y_D2y;
    Double I_ESP_K6yz_F3y = I_ESP_L7yz_D2y+ABY*I_ESP_K6yz_D2y;
    Double I_ESP_K5y2z_F3y = I_ESP_L6y2z_D2y+ABY*I_ESP_K5y2z_D2y;
    Double I_ESP_K4y3z_F3y = I_ESP_L5y3z_D2y+ABY*I_ESP_K4y3z_D2y;
    Double I_ESP_K3y4z_F3y = I_ESP_L4y4z_D2y+ABY*I_ESP_K3y4z_D2y;
    Double I_ESP_K2y5z_F3y = I_ESP_L3y5z_D2y+ABY*I_ESP_K2y5z_D2y;
    Double I_ESP_Ky6z_F3y = I_ESP_L2y6z_D2y+ABY*I_ESP_Ky6z_D2y;
    Double I_ESP_K7z_F3y = I_ESP_Ly7z_D2y+ABY*I_ESP_K7z_D2y;
    Double I_ESP_K6xz_F2yz = I_ESP_L6x2z_D2y+ABZ*I_ESP_K6xz_D2y;
    Double I_ESP_K5xyz_F2yz = I_ESP_L5xy2z_D2y+ABZ*I_ESP_K5xyz_D2y;
    Double I_ESP_K5x2z_F2yz = I_ESP_L5x3z_D2y+ABZ*I_ESP_K5x2z_D2y;
    Double I_ESP_K4x2yz_F2yz = I_ESP_L4x2y2z_D2y+ABZ*I_ESP_K4x2yz_D2y;
    Double I_ESP_K4xy2z_F2yz = I_ESP_L4xy3z_D2y+ABZ*I_ESP_K4xy2z_D2y;
    Double I_ESP_K4x3z_F2yz = I_ESP_L4x4z_D2y+ABZ*I_ESP_K4x3z_D2y;
    Double I_ESP_K3x3yz_F2yz = I_ESP_L3x3y2z_D2y+ABZ*I_ESP_K3x3yz_D2y;
    Double I_ESP_K3x2y2z_F2yz = I_ESP_L3x2y3z_D2y+ABZ*I_ESP_K3x2y2z_D2y;
    Double I_ESP_K3xy3z_F2yz = I_ESP_L3xy4z_D2y+ABZ*I_ESP_K3xy3z_D2y;
    Double I_ESP_K3x4z_F2yz = I_ESP_L3x5z_D2y+ABZ*I_ESP_K3x4z_D2y;
    Double I_ESP_K2x4yz_F2yz = I_ESP_L2x4y2z_D2y+ABZ*I_ESP_K2x4yz_D2y;
    Double I_ESP_K2x3y2z_F2yz = I_ESP_L2x3y3z_D2y+ABZ*I_ESP_K2x3y2z_D2y;
    Double I_ESP_K2x2y3z_F2yz = I_ESP_L2x2y4z_D2y+ABZ*I_ESP_K2x2y3z_D2y;
    Double I_ESP_K2xy4z_F2yz = I_ESP_L2xy5z_D2y+ABZ*I_ESP_K2xy4z_D2y;
    Double I_ESP_K2x5z_F2yz = I_ESP_L2x6z_D2y+ABZ*I_ESP_K2x5z_D2y;
    Double I_ESP_Kx5yz_F2yz = I_ESP_Lx5y2z_D2y+ABZ*I_ESP_Kx5yz_D2y;
    Double I_ESP_Kx4y2z_F2yz = I_ESP_Lx4y3z_D2y+ABZ*I_ESP_Kx4y2z_D2y;
    Double I_ESP_Kx3y3z_F2yz = I_ESP_Lx3y4z_D2y+ABZ*I_ESP_Kx3y3z_D2y;
    Double I_ESP_Kx2y4z_F2yz = I_ESP_Lx2y5z_D2y+ABZ*I_ESP_Kx2y4z_D2y;
    Double I_ESP_Kxy5z_F2yz = I_ESP_Lxy6z_D2y+ABZ*I_ESP_Kxy5z_D2y;
    Double I_ESP_Kx6z_F2yz = I_ESP_Lx7z_D2y+ABZ*I_ESP_Kx6z_D2y;
    Double I_ESP_K7x_F3z = I_ESP_L7xz_D2z+ABZ*I_ESP_K7x_D2z;
    Double I_ESP_K6xy_F3z = I_ESP_L6xyz_D2z+ABZ*I_ESP_K6xy_D2z;
    Double I_ESP_K6xz_F3z = I_ESP_L6x2z_D2z+ABZ*I_ESP_K6xz_D2z;
    Double I_ESP_K5x2y_F3z = I_ESP_L5x2yz_D2z+ABZ*I_ESP_K5x2y_D2z;
    Double I_ESP_K5xyz_F3z = I_ESP_L5xy2z_D2z+ABZ*I_ESP_K5xyz_D2z;
    Double I_ESP_K5x2z_F3z = I_ESP_L5x3z_D2z+ABZ*I_ESP_K5x2z_D2z;
    Double I_ESP_K4x3y_F3z = I_ESP_L4x3yz_D2z+ABZ*I_ESP_K4x3y_D2z;
    Double I_ESP_K4x2yz_F3z = I_ESP_L4x2y2z_D2z+ABZ*I_ESP_K4x2yz_D2z;
    Double I_ESP_K4xy2z_F3z = I_ESP_L4xy3z_D2z+ABZ*I_ESP_K4xy2z_D2z;
    Double I_ESP_K4x3z_F3z = I_ESP_L4x4z_D2z+ABZ*I_ESP_K4x3z_D2z;
    Double I_ESP_K3x4y_F3z = I_ESP_L3x4yz_D2z+ABZ*I_ESP_K3x4y_D2z;
    Double I_ESP_K3x3yz_F3z = I_ESP_L3x3y2z_D2z+ABZ*I_ESP_K3x3yz_D2z;
    Double I_ESP_K3x2y2z_F3z = I_ESP_L3x2y3z_D2z+ABZ*I_ESP_K3x2y2z_D2z;
    Double I_ESP_K3xy3z_F3z = I_ESP_L3xy4z_D2z+ABZ*I_ESP_K3xy3z_D2z;
    Double I_ESP_K3x4z_F3z = I_ESP_L3x5z_D2z+ABZ*I_ESP_K3x4z_D2z;
    Double I_ESP_K2x5y_F3z = I_ESP_L2x5yz_D2z+ABZ*I_ESP_K2x5y_D2z;
    Double I_ESP_K2x4yz_F3z = I_ESP_L2x4y2z_D2z+ABZ*I_ESP_K2x4yz_D2z;
    Double I_ESP_K2x3y2z_F3z = I_ESP_L2x3y3z_D2z+ABZ*I_ESP_K2x3y2z_D2z;
    Double I_ESP_K2x2y3z_F3z = I_ESP_L2x2y4z_D2z+ABZ*I_ESP_K2x2y3z_D2z;
    Double I_ESP_K2xy4z_F3z = I_ESP_L2xy5z_D2z+ABZ*I_ESP_K2xy4z_D2z;
    Double I_ESP_K2x5z_F3z = I_ESP_L2x6z_D2z+ABZ*I_ESP_K2x5z_D2z;
    Double I_ESP_Kx6y_F3z = I_ESP_Lx6yz_D2z+ABZ*I_ESP_Kx6y_D2z;
    Double I_ESP_Kx5yz_F3z = I_ESP_Lx5y2z_D2z+ABZ*I_ESP_Kx5yz_D2z;
    Double I_ESP_Kx4y2z_F3z = I_ESP_Lx4y3z_D2z+ABZ*I_ESP_Kx4y2z_D2z;
    Double I_ESP_Kx3y3z_F3z = I_ESP_Lx3y4z_D2z+ABZ*I_ESP_Kx3y3z_D2z;
    Double I_ESP_Kx2y4z_F3z = I_ESP_Lx2y5z_D2z+ABZ*I_ESP_Kx2y4z_D2z;
    Double I_ESP_Kxy5z_F3z = I_ESP_Lxy6z_D2z+ABZ*I_ESP_Kxy5z_D2z;
    Double I_ESP_Kx6z_F3z = I_ESP_Lx7z_D2z+ABZ*I_ESP_Kx6z_D2z;
    Double I_ESP_K7y_F3z = I_ESP_L7yz_D2z+ABZ*I_ESP_K7y_D2z;
    Double I_ESP_K6yz_F3z = I_ESP_L6y2z_D2z+ABZ*I_ESP_K6yz_D2z;
    Double I_ESP_K5y2z_F3z = I_ESP_L5y3z_D2z+ABZ*I_ESP_K5y2z_D2z;
    Double I_ESP_K4y3z_F3z = I_ESP_L4y4z_D2z+ABZ*I_ESP_K4y3z_D2z;
    Double I_ESP_K3y4z_F3z = I_ESP_L3y5z_D2z+ABZ*I_ESP_K3y4z_D2z;
    Double I_ESP_K2y5z_F3z = I_ESP_L2y6z_D2z+ABZ*I_ESP_K2y5z_D2z;
    Double I_ESP_Ky6z_F3z = I_ESP_Ly7z_D2z+ABZ*I_ESP_Ky6z_D2z;
    Double I_ESP_K7z_F3z = I_ESP_L8z_D2z+ABZ*I_ESP_K7z_D2z;

    /************************************************************
     * shell quartet name: SQ_ESP_I_G
     * expanding position: BRA2
     * code section is: HRR
     * totally 129 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_K_F
     * RHS shell quartet name: SQ_ESP_I_F
     ************************************************************/
    Double I_ESP_I6x_G4x = I_ESP_K7x_F3x+ABX*I_ESP_I6x_F3x;
    Double I_ESP_I5xy_G4x = I_ESP_K6xy_F3x+ABX*I_ESP_I5xy_F3x;
    Double I_ESP_I5xz_G4x = I_ESP_K6xz_F3x+ABX*I_ESP_I5xz_F3x;
    Double I_ESP_I4x2y_G4x = I_ESP_K5x2y_F3x+ABX*I_ESP_I4x2y_F3x;
    Double I_ESP_I4xyz_G4x = I_ESP_K5xyz_F3x+ABX*I_ESP_I4xyz_F3x;
    Double I_ESP_I4x2z_G4x = I_ESP_K5x2z_F3x+ABX*I_ESP_I4x2z_F3x;
    Double I_ESP_I3x3y_G4x = I_ESP_K4x3y_F3x+ABX*I_ESP_I3x3y_F3x;
    Double I_ESP_I3x2yz_G4x = I_ESP_K4x2yz_F3x+ABX*I_ESP_I3x2yz_F3x;
    Double I_ESP_I3xy2z_G4x = I_ESP_K4xy2z_F3x+ABX*I_ESP_I3xy2z_F3x;
    Double I_ESP_I3x3z_G4x = I_ESP_K4x3z_F3x+ABX*I_ESP_I3x3z_F3x;
    Double I_ESP_I2x4y_G4x = I_ESP_K3x4y_F3x+ABX*I_ESP_I2x4y_F3x;
    Double I_ESP_I2x3yz_G4x = I_ESP_K3x3yz_F3x+ABX*I_ESP_I2x3yz_F3x;
    Double I_ESP_I2x2y2z_G4x = I_ESP_K3x2y2z_F3x+ABX*I_ESP_I2x2y2z_F3x;
    Double I_ESP_I2xy3z_G4x = I_ESP_K3xy3z_F3x+ABX*I_ESP_I2xy3z_F3x;
    Double I_ESP_I2x4z_G4x = I_ESP_K3x4z_F3x+ABX*I_ESP_I2x4z_F3x;
    Double I_ESP_Ix5y_G4x = I_ESP_K2x5y_F3x+ABX*I_ESP_Ix5y_F3x;
    Double I_ESP_Ix4yz_G4x = I_ESP_K2x4yz_F3x+ABX*I_ESP_Ix4yz_F3x;
    Double I_ESP_Ix3y2z_G4x = I_ESP_K2x3y2z_F3x+ABX*I_ESP_Ix3y2z_F3x;
    Double I_ESP_Ix2y3z_G4x = I_ESP_K2x2y3z_F3x+ABX*I_ESP_Ix2y3z_F3x;
    Double I_ESP_Ixy4z_G4x = I_ESP_K2xy4z_F3x+ABX*I_ESP_Ixy4z_F3x;
    Double I_ESP_Ix5z_G4x = I_ESP_K2x5z_F3x+ABX*I_ESP_Ix5z_F3x;
    Double I_ESP_I6y_G4x = I_ESP_Kx6y_F3x+ABX*I_ESP_I6y_F3x;
    Double I_ESP_I5yz_G4x = I_ESP_Kx5yz_F3x+ABX*I_ESP_I5yz_F3x;
    Double I_ESP_I4y2z_G4x = I_ESP_Kx4y2z_F3x+ABX*I_ESP_I4y2z_F3x;
    Double I_ESP_I3y3z_G4x = I_ESP_Kx3y3z_F3x+ABX*I_ESP_I3y3z_F3x;
    Double I_ESP_I2y4z_G4x = I_ESP_Kx2y4z_F3x+ABX*I_ESP_I2y4z_F3x;
    Double I_ESP_Iy5z_G4x = I_ESP_Kxy5z_F3x+ABX*I_ESP_Iy5z_F3x;
    Double I_ESP_I6z_G4x = I_ESP_Kx6z_F3x+ABX*I_ESP_I6z_F3x;
    Double I_ESP_I5xy_G3xy = I_ESP_K5x2y_F3x+ABY*I_ESP_I5xy_F3x;
    Double I_ESP_I5xz_G3xy = I_ESP_K5xyz_F3x+ABY*I_ESP_I5xz_F3x;
    Double I_ESP_I4x2y_G3xy = I_ESP_K4x3y_F3x+ABY*I_ESP_I4x2y_F3x;
    Double I_ESP_I4xyz_G3xy = I_ESP_K4x2yz_F3x+ABY*I_ESP_I4xyz_F3x;
    Double I_ESP_I4x2z_G3xy = I_ESP_K4xy2z_F3x+ABY*I_ESP_I4x2z_F3x;
    Double I_ESP_I3x3y_G3xy = I_ESP_K3x4y_F3x+ABY*I_ESP_I3x3y_F3x;
    Double I_ESP_I3x2yz_G3xy = I_ESP_K3x3yz_F3x+ABY*I_ESP_I3x2yz_F3x;
    Double I_ESP_I3xy2z_G3xy = I_ESP_K3x2y2z_F3x+ABY*I_ESP_I3xy2z_F3x;
    Double I_ESP_I3x3z_G3xy = I_ESP_K3xy3z_F3x+ABY*I_ESP_I3x3z_F3x;
    Double I_ESP_I2x4y_G3xy = I_ESP_K2x5y_F3x+ABY*I_ESP_I2x4y_F3x;
    Double I_ESP_I2x3yz_G3xy = I_ESP_K2x4yz_F3x+ABY*I_ESP_I2x3yz_F3x;
    Double I_ESP_I2x2y2z_G3xy = I_ESP_K2x3y2z_F3x+ABY*I_ESP_I2x2y2z_F3x;
    Double I_ESP_I2xy3z_G3xy = I_ESP_K2x2y3z_F3x+ABY*I_ESP_I2xy3z_F3x;
    Double I_ESP_I2x4z_G3xy = I_ESP_K2xy4z_F3x+ABY*I_ESP_I2x4z_F3x;
    Double I_ESP_Ix5y_G3xy = I_ESP_Kx6y_F3x+ABY*I_ESP_Ix5y_F3x;
    Double I_ESP_Ix4yz_G3xy = I_ESP_Kx5yz_F3x+ABY*I_ESP_Ix4yz_F3x;
    Double I_ESP_Ix3y2z_G3xy = I_ESP_Kx4y2z_F3x+ABY*I_ESP_Ix3y2z_F3x;
    Double I_ESP_Ix2y3z_G3xy = I_ESP_Kx3y3z_F3x+ABY*I_ESP_Ix2y3z_F3x;
    Double I_ESP_Ixy4z_G3xy = I_ESP_Kx2y4z_F3x+ABY*I_ESP_Ixy4z_F3x;
    Double I_ESP_Ix5z_G3xy = I_ESP_Kxy5z_F3x+ABY*I_ESP_Ix5z_F3x;
    Double I_ESP_I6y_G3xy = I_ESP_K7y_F3x+ABY*I_ESP_I6y_F3x;
    Double I_ESP_I5yz_G3xy = I_ESP_K6yz_F3x+ABY*I_ESP_I5yz_F3x;
    Double I_ESP_I4y2z_G3xy = I_ESP_K5y2z_F3x+ABY*I_ESP_I4y2z_F3x;
    Double I_ESP_I3y3z_G3xy = I_ESP_K4y3z_F3x+ABY*I_ESP_I3y3z_F3x;
    Double I_ESP_I2y4z_G3xy = I_ESP_K3y4z_F3x+ABY*I_ESP_I2y4z_F3x;
    Double I_ESP_Iy5z_G3xy = I_ESP_K2y5z_F3x+ABY*I_ESP_Iy5z_F3x;
    Double I_ESP_I6z_G3xy = I_ESP_Ky6z_F3x+ABY*I_ESP_I6z_F3x;
    Double I_ESP_I5xz_G3xz = I_ESP_K5x2z_F3x+ABZ*I_ESP_I5xz_F3x;
    Double I_ESP_I4xyz_G3xz = I_ESP_K4xy2z_F3x+ABZ*I_ESP_I4xyz_F3x;
    Double I_ESP_I4x2z_G3xz = I_ESP_K4x3z_F3x+ABZ*I_ESP_I4x2z_F3x;
    Double I_ESP_I3x2yz_G3xz = I_ESP_K3x2y2z_F3x+ABZ*I_ESP_I3x2yz_F3x;
    Double I_ESP_I3xy2z_G3xz = I_ESP_K3xy3z_F3x+ABZ*I_ESP_I3xy2z_F3x;
    Double I_ESP_I3x3z_G3xz = I_ESP_K3x4z_F3x+ABZ*I_ESP_I3x3z_F3x;
    Double I_ESP_I2x3yz_G3xz = I_ESP_K2x3y2z_F3x+ABZ*I_ESP_I2x3yz_F3x;
    Double I_ESP_I2x2y2z_G3xz = I_ESP_K2x2y3z_F3x+ABZ*I_ESP_I2x2y2z_F3x;
    Double I_ESP_I2xy3z_G3xz = I_ESP_K2xy4z_F3x+ABZ*I_ESP_I2xy3z_F3x;
    Double I_ESP_I2x4z_G3xz = I_ESP_K2x5z_F3x+ABZ*I_ESP_I2x4z_F3x;
    Double I_ESP_Ix4yz_G3xz = I_ESP_Kx4y2z_F3x+ABZ*I_ESP_Ix4yz_F3x;
    Double I_ESP_Ix3y2z_G3xz = I_ESP_Kx3y3z_F3x+ABZ*I_ESP_Ix3y2z_F3x;
    Double I_ESP_Ix2y3z_G3xz = I_ESP_Kx2y4z_F3x+ABZ*I_ESP_Ix2y3z_F3x;
    Double I_ESP_Ixy4z_G3xz = I_ESP_Kxy5z_F3x+ABZ*I_ESP_Ixy4z_F3x;
    Double I_ESP_Ix5z_G3xz = I_ESP_Kx6z_F3x+ABZ*I_ESP_Ix5z_F3x;
    Double I_ESP_I5yz_G3xz = I_ESP_K5y2z_F3x+ABZ*I_ESP_I5yz_F3x;
    Double I_ESP_I4y2z_G3xz = I_ESP_K4y3z_F3x+ABZ*I_ESP_I4y2z_F3x;
    Double I_ESP_I3y3z_G3xz = I_ESP_K3y4z_F3x+ABZ*I_ESP_I3y3z_F3x;
    Double I_ESP_I2y4z_G3xz = I_ESP_K2y5z_F3x+ABZ*I_ESP_I2y4z_F3x;
    Double I_ESP_Iy5z_G3xz = I_ESP_Ky6z_F3x+ABZ*I_ESP_Iy5z_F3x;
    Double I_ESP_I6z_G3xz = I_ESP_K7z_F3x+ABZ*I_ESP_I6z_F3x;
    Double I_ESP_I5xz_G2x2y = I_ESP_K5xyz_F2xy+ABY*I_ESP_I5xz_F2xy;
    Double I_ESP_I4xyz_G2x2y = I_ESP_K4x2yz_F2xy+ABY*I_ESP_I4xyz_F2xy;
    Double I_ESP_I4x2z_G2x2y = I_ESP_K4xy2z_F2xy+ABY*I_ESP_I4x2z_F2xy;
    Double I_ESP_I3x2yz_G2x2y = I_ESP_K3x3yz_F2xy+ABY*I_ESP_I3x2yz_F2xy;
    Double I_ESP_I3xy2z_G2x2y = I_ESP_K3x2y2z_F2xy+ABY*I_ESP_I3xy2z_F2xy;
    Double I_ESP_I3x3z_G2x2y = I_ESP_K3xy3z_F2xy+ABY*I_ESP_I3x3z_F2xy;
    Double I_ESP_I2x3yz_G2x2y = I_ESP_K2x4yz_F2xy+ABY*I_ESP_I2x3yz_F2xy;
    Double I_ESP_I2x2y2z_G2x2y = I_ESP_K2x3y2z_F2xy+ABY*I_ESP_I2x2y2z_F2xy;
    Double I_ESP_I2xy3z_G2x2y = I_ESP_K2x2y3z_F2xy+ABY*I_ESP_I2xy3z_F2xy;
    Double I_ESP_I2x4z_G2x2y = I_ESP_K2xy4z_F2xy+ABY*I_ESP_I2x4z_F2xy;
    Double I_ESP_Ix4yz_G2x2y = I_ESP_Kx5yz_F2xy+ABY*I_ESP_Ix4yz_F2xy;
    Double I_ESP_Ix3y2z_G2x2y = I_ESP_Kx4y2z_F2xy+ABY*I_ESP_Ix3y2z_F2xy;
    Double I_ESP_Ix2y3z_G2x2y = I_ESP_Kx3y3z_F2xy+ABY*I_ESP_Ix2y3z_F2xy;
    Double I_ESP_Ixy4z_G2x2y = I_ESP_Kx2y4z_F2xy+ABY*I_ESP_Ixy4z_F2xy;
    Double I_ESP_Ix5z_G2x2y = I_ESP_Kxy5z_F2xy+ABY*I_ESP_Ix5z_F2xy;
    Double I_ESP_I5yz_G2x2y = I_ESP_K6yz_F2xy+ABY*I_ESP_I5yz_F2xy;
    Double I_ESP_I4y2z_G2x2y = I_ESP_K5y2z_F2xy+ABY*I_ESP_I4y2z_F2xy;
    Double I_ESP_I3y3z_G2x2y = I_ESP_K4y3z_F2xy+ABY*I_ESP_I3y3z_F2xy;
    Double I_ESP_I2y4z_G2x2y = I_ESP_K3y4z_F2xy+ABY*I_ESP_I2y4z_F2xy;
    Double I_ESP_Iy5z_G2x2y = I_ESP_K2y5z_F2xy+ABY*I_ESP_Iy5z_F2xy;
    Double I_ESP_I6z_G2x2y = I_ESP_Ky6z_F2xy+ABY*I_ESP_I6z_F2xy;
    Double I_ESP_I5xy_G2x2z = I_ESP_K5xyz_F2xz+ABZ*I_ESP_I5xy_F2xz;
    Double I_ESP_I4x2y_G2x2z = I_ESP_K4x2yz_F2xz+ABZ*I_ESP_I4x2y_F2xz;
    Double I_ESP_I4xyz_G2x2z = I_ESP_K4xy2z_F2xz+ABZ*I_ESP_I4xyz_F2xz;
    Double I_ESP_I3x3y_G2x2z = I_ESP_K3x3yz_F2xz+ABZ*I_ESP_I3x3y_F2xz;
    Double I_ESP_I3x2yz_G2x2z = I_ESP_K3x2y2z_F2xz+ABZ*I_ESP_I3x2yz_F2xz;
    Double I_ESP_I3xy2z_G2x2z = I_ESP_K3xy3z_F2xz+ABZ*I_ESP_I3xy2z_F2xz;
    Double I_ESP_I2x4y_G2x2z = I_ESP_K2x4yz_F2xz+ABZ*I_ESP_I2x4y_F2xz;
    Double I_ESP_I2x3yz_G2x2z = I_ESP_K2x3y2z_F2xz+ABZ*I_ESP_I2x3yz_F2xz;
    Double I_ESP_I2x2y2z_G2x2z = I_ESP_K2x2y3z_F2xz+ABZ*I_ESP_I2x2y2z_F2xz;
    Double I_ESP_I2xy3z_G2x2z = I_ESP_K2xy4z_F2xz+ABZ*I_ESP_I2xy3z_F2xz;
    Double I_ESP_Ix5y_G2x2z = I_ESP_Kx5yz_F2xz+ABZ*I_ESP_Ix5y_F2xz;
    Double I_ESP_Ix4yz_G2x2z = I_ESP_Kx4y2z_F2xz+ABZ*I_ESP_Ix4yz_F2xz;
    Double I_ESP_Ix3y2z_G2x2z = I_ESP_Kx3y3z_F2xz+ABZ*I_ESP_Ix3y2z_F2xz;
    Double I_ESP_Ix2y3z_G2x2z = I_ESP_Kx2y4z_F2xz+ABZ*I_ESP_Ix2y3z_F2xz;
    Double I_ESP_Ixy4z_G2x2z = I_ESP_Kxy5z_F2xz+ABZ*I_ESP_Ixy4z_F2xz;
    Double I_ESP_I6y_G2x2z = I_ESP_K6yz_F2xz+ABZ*I_ESP_I6y_F2xz;
    Double I_ESP_I5yz_G2x2z = I_ESP_K5y2z_F2xz+ABZ*I_ESP_I5yz_F2xz;
    Double I_ESP_I4y2z_G2x2z = I_ESP_K4y3z_F2xz+ABZ*I_ESP_I4y2z_F2xz;
    Double I_ESP_I3y3z_G2x2z = I_ESP_K3y4z_F2xz+ABZ*I_ESP_I3y3z_F2xz;
    Double I_ESP_I2y4z_G2x2z = I_ESP_K2y5z_F2xz+ABZ*I_ESP_I2y4z_F2xz;
    Double I_ESP_Iy5z_G2x2z = I_ESP_Ky6z_F2xz+ABZ*I_ESP_Iy5z_F2xz;
    Double I_ESP_I6x_Gx3y = I_ESP_K7x_F3y+ABX*I_ESP_I6x_F3y;
    Double I_ESP_I5xy_Gx3y = I_ESP_K6xy_F3y+ABX*I_ESP_I5xy_F3y;
    Double I_ESP_I5xz_Gx3y = I_ESP_K6xz_F3y+ABX*I_ESP_I5xz_F3y;
    Double I_ESP_I4x2y_Gx3y = I_ESP_K5x2y_F3y+ABX*I_ESP_I4x2y_F3y;
    Double I_ESP_I4xyz_Gx3y = I_ESP_K5xyz_F3y+ABX*I_ESP_I4xyz_F3y;
    Double I_ESP_I4x2z_Gx3y = I_ESP_K5x2z_F3y+ABX*I_ESP_I4x2z_F3y;
    Double I_ESP_I3x3y_Gx3y = I_ESP_K4x3y_F3y+ABX*I_ESP_I3x3y_F3y;
    Double I_ESP_I3x2yz_Gx3y = I_ESP_K4x2yz_F3y+ABX*I_ESP_I3x2yz_F3y;
    Double I_ESP_I3xy2z_Gx3y = I_ESP_K4xy2z_F3y+ABX*I_ESP_I3xy2z_F3y;
    Double I_ESP_I3x3z_Gx3y = I_ESP_K4x3z_F3y+ABX*I_ESP_I3x3z_F3y;
    Double I_ESP_I2x4y_Gx3y = I_ESP_K3x4y_F3y+ABX*I_ESP_I2x4y_F3y;
    Double I_ESP_I2x3yz_Gx3y = I_ESP_K3x3yz_F3y+ABX*I_ESP_I2x3yz_F3y;
    Double I_ESP_I2x2y2z_Gx3y = I_ESP_K3x2y2z_F3y+ABX*I_ESP_I2x2y2z_F3y;
    Double I_ESP_I2xy3z_Gx3y = I_ESP_K3xy3z_F3y+ABX*I_ESP_I2xy3z_F3y;
    Double I_ESP_I2x4z_Gx3y = I_ESP_K3x4z_F3y+ABX*I_ESP_I2x4z_F3y;
    Double I_ESP_Ix5y_Gx3y = I_ESP_K2x5y_F3y+ABX*I_ESP_Ix5y_F3y;
    Double I_ESP_Ix4yz_Gx3y = I_ESP_K2x4yz_F3y+ABX*I_ESP_Ix4yz_F3y;
    Double I_ESP_Ix3y2z_Gx3y = I_ESP_K2x3y2z_F3y+ABX*I_ESP_Ix3y2z_F3y;
    Double I_ESP_Ix2y3z_Gx3y = I_ESP_K2x2y3z_F3y+ABX*I_ESP_Ix2y3z_F3y;
    Double I_ESP_Ixy4z_Gx3y = I_ESP_K2xy4z_F3y+ABX*I_ESP_Ixy4z_F3y;
    Double I_ESP_Ix5z_Gx3y = I_ESP_K2x5z_F3y+ABX*I_ESP_Ix5z_F3y;
    Double I_ESP_I5yz_Gx3y = I_ESP_Kx5yz_F3y+ABX*I_ESP_I5yz_F3y;
    Double I_ESP_I4y2z_Gx3y = I_ESP_Kx4y2z_F3y+ABX*I_ESP_I4y2z_F3y;
    Double I_ESP_I3y3z_Gx3y = I_ESP_Kx3y3z_F3y+ABX*I_ESP_I3y3z_F3y;
    Double I_ESP_I2y4z_Gx3y = I_ESP_Kx2y4z_F3y+ABX*I_ESP_I2y4z_F3y;
    Double I_ESP_Iy5z_Gx3y = I_ESP_Kxy5z_F3y+ABX*I_ESP_Iy5z_F3y;
    Double I_ESP_I6z_Gx3y = I_ESP_Kx6z_F3y+ABX*I_ESP_I6z_F3y;
    Double I_ESP_I6x_Gx3z = I_ESP_K7x_F3z+ABX*I_ESP_I6x_F3z;
    Double I_ESP_I5xy_Gx3z = I_ESP_K6xy_F3z+ABX*I_ESP_I5xy_F3z;
    Double I_ESP_I5xz_Gx3z = I_ESP_K6xz_F3z+ABX*I_ESP_I5xz_F3z;
    Double I_ESP_I4x2y_Gx3z = I_ESP_K5x2y_F3z+ABX*I_ESP_I4x2y_F3z;
    Double I_ESP_I4xyz_Gx3z = I_ESP_K5xyz_F3z+ABX*I_ESP_I4xyz_F3z;
    Double I_ESP_I4x2z_Gx3z = I_ESP_K5x2z_F3z+ABX*I_ESP_I4x2z_F3z;
    Double I_ESP_I3x3y_Gx3z = I_ESP_K4x3y_F3z+ABX*I_ESP_I3x3y_F3z;
    Double I_ESP_I3x2yz_Gx3z = I_ESP_K4x2yz_F3z+ABX*I_ESP_I3x2yz_F3z;
    Double I_ESP_I3xy2z_Gx3z = I_ESP_K4xy2z_F3z+ABX*I_ESP_I3xy2z_F3z;
    Double I_ESP_I3x3z_Gx3z = I_ESP_K4x3z_F3z+ABX*I_ESP_I3x3z_F3z;
    Double I_ESP_I2x4y_Gx3z = I_ESP_K3x4y_F3z+ABX*I_ESP_I2x4y_F3z;
    Double I_ESP_I2x3yz_Gx3z = I_ESP_K3x3yz_F3z+ABX*I_ESP_I2x3yz_F3z;
    Double I_ESP_I2x2y2z_Gx3z = I_ESP_K3x2y2z_F3z+ABX*I_ESP_I2x2y2z_F3z;
    Double I_ESP_I2xy3z_Gx3z = I_ESP_K3xy3z_F3z+ABX*I_ESP_I2xy3z_F3z;
    Double I_ESP_I2x4z_Gx3z = I_ESP_K3x4z_F3z+ABX*I_ESP_I2x4z_F3z;
    Double I_ESP_Ix5y_Gx3z = I_ESP_K2x5y_F3z+ABX*I_ESP_Ix5y_F3z;
    Double I_ESP_Ix4yz_Gx3z = I_ESP_K2x4yz_F3z+ABX*I_ESP_Ix4yz_F3z;
    Double I_ESP_Ix3y2z_Gx3z = I_ESP_K2x3y2z_F3z+ABX*I_ESP_Ix3y2z_F3z;
    Double I_ESP_Ix2y3z_Gx3z = I_ESP_K2x2y3z_F3z+ABX*I_ESP_Ix2y3z_F3z;
    Double I_ESP_Ixy4z_Gx3z = I_ESP_K2xy4z_F3z+ABX*I_ESP_Ixy4z_F3z;
    Double I_ESP_Ix5z_Gx3z = I_ESP_K2x5z_F3z+ABX*I_ESP_Ix5z_F3z;
    Double I_ESP_I6y_Gx3z = I_ESP_Kx6y_F3z+ABX*I_ESP_I6y_F3z;
    Double I_ESP_I5yz_Gx3z = I_ESP_Kx5yz_F3z+ABX*I_ESP_I5yz_F3z;
    Double I_ESP_I4y2z_Gx3z = I_ESP_Kx4y2z_F3z+ABX*I_ESP_I4y2z_F3z;
    Double I_ESP_I3y3z_Gx3z = I_ESP_Kx3y3z_F3z+ABX*I_ESP_I3y3z_F3z;
    Double I_ESP_I2y4z_Gx3z = I_ESP_Kx2y4z_F3z+ABX*I_ESP_I2y4z_F3z;
    Double I_ESP_Iy5z_Gx3z = I_ESP_Kxy5z_F3z+ABX*I_ESP_Iy5z_F3z;
    Double I_ESP_I6x_G4y = I_ESP_K6xy_F3y+ABY*I_ESP_I6x_F3y;
    Double I_ESP_I5xy_G4y = I_ESP_K5x2y_F3y+ABY*I_ESP_I5xy_F3y;
    Double I_ESP_I5xz_G4y = I_ESP_K5xyz_F3y+ABY*I_ESP_I5xz_F3y;
    Double I_ESP_I4x2y_G4y = I_ESP_K4x3y_F3y+ABY*I_ESP_I4x2y_F3y;
    Double I_ESP_I4xyz_G4y = I_ESP_K4x2yz_F3y+ABY*I_ESP_I4xyz_F3y;
    Double I_ESP_I4x2z_G4y = I_ESP_K4xy2z_F3y+ABY*I_ESP_I4x2z_F3y;
    Double I_ESP_I3x3y_G4y = I_ESP_K3x4y_F3y+ABY*I_ESP_I3x3y_F3y;
    Double I_ESP_I3x2yz_G4y = I_ESP_K3x3yz_F3y+ABY*I_ESP_I3x2yz_F3y;
    Double I_ESP_I3xy2z_G4y = I_ESP_K3x2y2z_F3y+ABY*I_ESP_I3xy2z_F3y;
    Double I_ESP_I3x3z_G4y = I_ESP_K3xy3z_F3y+ABY*I_ESP_I3x3z_F3y;
    Double I_ESP_I2x4y_G4y = I_ESP_K2x5y_F3y+ABY*I_ESP_I2x4y_F3y;
    Double I_ESP_I2x3yz_G4y = I_ESP_K2x4yz_F3y+ABY*I_ESP_I2x3yz_F3y;
    Double I_ESP_I2x2y2z_G4y = I_ESP_K2x3y2z_F3y+ABY*I_ESP_I2x2y2z_F3y;
    Double I_ESP_I2xy3z_G4y = I_ESP_K2x2y3z_F3y+ABY*I_ESP_I2xy3z_F3y;
    Double I_ESP_I2x4z_G4y = I_ESP_K2xy4z_F3y+ABY*I_ESP_I2x4z_F3y;
    Double I_ESP_Ix5y_G4y = I_ESP_Kx6y_F3y+ABY*I_ESP_Ix5y_F3y;
    Double I_ESP_Ix4yz_G4y = I_ESP_Kx5yz_F3y+ABY*I_ESP_Ix4yz_F3y;
    Double I_ESP_Ix3y2z_G4y = I_ESP_Kx4y2z_F3y+ABY*I_ESP_Ix3y2z_F3y;
    Double I_ESP_Ix2y3z_G4y = I_ESP_Kx3y3z_F3y+ABY*I_ESP_Ix2y3z_F3y;
    Double I_ESP_Ixy4z_G4y = I_ESP_Kx2y4z_F3y+ABY*I_ESP_Ixy4z_F3y;
    Double I_ESP_Ix5z_G4y = I_ESP_Kxy5z_F3y+ABY*I_ESP_Ix5z_F3y;
    Double I_ESP_I6y_G4y = I_ESP_K7y_F3y+ABY*I_ESP_I6y_F3y;
    Double I_ESP_I5yz_G4y = I_ESP_K6yz_F3y+ABY*I_ESP_I5yz_F3y;
    Double I_ESP_I4y2z_G4y = I_ESP_K5y2z_F3y+ABY*I_ESP_I4y2z_F3y;
    Double I_ESP_I3y3z_G4y = I_ESP_K4y3z_F3y+ABY*I_ESP_I3y3z_F3y;
    Double I_ESP_I2y4z_G4y = I_ESP_K3y4z_F3y+ABY*I_ESP_I2y4z_F3y;
    Double I_ESP_Iy5z_G4y = I_ESP_K2y5z_F3y+ABY*I_ESP_Iy5z_F3y;
    Double I_ESP_I6z_G4y = I_ESP_Ky6z_F3y+ABY*I_ESP_I6z_F3y;
    Double I_ESP_I5xz_G3yz = I_ESP_K5x2z_F3y+ABZ*I_ESP_I5xz_F3y;
    Double I_ESP_I4xyz_G3yz = I_ESP_K4xy2z_F3y+ABZ*I_ESP_I4xyz_F3y;
    Double I_ESP_I4x2z_G3yz = I_ESP_K4x3z_F3y+ABZ*I_ESP_I4x2z_F3y;
    Double I_ESP_I3x2yz_G3yz = I_ESP_K3x2y2z_F3y+ABZ*I_ESP_I3x2yz_F3y;
    Double I_ESP_I3xy2z_G3yz = I_ESP_K3xy3z_F3y+ABZ*I_ESP_I3xy2z_F3y;
    Double I_ESP_I3x3z_G3yz = I_ESP_K3x4z_F3y+ABZ*I_ESP_I3x3z_F3y;
    Double I_ESP_I2x3yz_G3yz = I_ESP_K2x3y2z_F3y+ABZ*I_ESP_I2x3yz_F3y;
    Double I_ESP_I2x2y2z_G3yz = I_ESP_K2x2y3z_F3y+ABZ*I_ESP_I2x2y2z_F3y;
    Double I_ESP_I2xy3z_G3yz = I_ESP_K2xy4z_F3y+ABZ*I_ESP_I2xy3z_F3y;
    Double I_ESP_I2x4z_G3yz = I_ESP_K2x5z_F3y+ABZ*I_ESP_I2x4z_F3y;
    Double I_ESP_Ix4yz_G3yz = I_ESP_Kx4y2z_F3y+ABZ*I_ESP_Ix4yz_F3y;
    Double I_ESP_Ix3y2z_G3yz = I_ESP_Kx3y3z_F3y+ABZ*I_ESP_Ix3y2z_F3y;
    Double I_ESP_Ix2y3z_G3yz = I_ESP_Kx2y4z_F3y+ABZ*I_ESP_Ix2y3z_F3y;
    Double I_ESP_Ixy4z_G3yz = I_ESP_Kxy5z_F3y+ABZ*I_ESP_Ixy4z_F3y;
    Double I_ESP_Ix5z_G3yz = I_ESP_Kx6z_F3y+ABZ*I_ESP_Ix5z_F3y;
    Double I_ESP_I5yz_G3yz = I_ESP_K5y2z_F3y+ABZ*I_ESP_I5yz_F3y;
    Double I_ESP_I4y2z_G3yz = I_ESP_K4y3z_F3y+ABZ*I_ESP_I4y2z_F3y;
    Double I_ESP_I3y3z_G3yz = I_ESP_K3y4z_F3y+ABZ*I_ESP_I3y3z_F3y;
    Double I_ESP_I2y4z_G3yz = I_ESP_K2y5z_F3y+ABZ*I_ESP_I2y4z_F3y;
    Double I_ESP_Iy5z_G3yz = I_ESP_Ky6z_F3y+ABZ*I_ESP_Iy5z_F3y;
    Double I_ESP_I6z_G3yz = I_ESP_K7z_F3y+ABZ*I_ESP_I6z_F3y;
    Double I_ESP_I6x_G2y2z = I_ESP_K6xz_F2yz+ABZ*I_ESP_I6x_F2yz;
    Double I_ESP_I5xy_G2y2z = I_ESP_K5xyz_F2yz+ABZ*I_ESP_I5xy_F2yz;
    Double I_ESP_I5xz_G2y2z = I_ESP_K5x2z_F2yz+ABZ*I_ESP_I5xz_F2yz;
    Double I_ESP_I4x2y_G2y2z = I_ESP_K4x2yz_F2yz+ABZ*I_ESP_I4x2y_F2yz;
    Double I_ESP_I4xyz_G2y2z = I_ESP_K4xy2z_F2yz+ABZ*I_ESP_I4xyz_F2yz;
    Double I_ESP_I4x2z_G2y2z = I_ESP_K4x3z_F2yz+ABZ*I_ESP_I4x2z_F2yz;
    Double I_ESP_I3x3y_G2y2z = I_ESP_K3x3yz_F2yz+ABZ*I_ESP_I3x3y_F2yz;
    Double I_ESP_I3x2yz_G2y2z = I_ESP_K3x2y2z_F2yz+ABZ*I_ESP_I3x2yz_F2yz;
    Double I_ESP_I3xy2z_G2y2z = I_ESP_K3xy3z_F2yz+ABZ*I_ESP_I3xy2z_F2yz;
    Double I_ESP_I3x3z_G2y2z = I_ESP_K3x4z_F2yz+ABZ*I_ESP_I3x3z_F2yz;
    Double I_ESP_I2x4y_G2y2z = I_ESP_K2x4yz_F2yz+ABZ*I_ESP_I2x4y_F2yz;
    Double I_ESP_I2x3yz_G2y2z = I_ESP_K2x3y2z_F2yz+ABZ*I_ESP_I2x3yz_F2yz;
    Double I_ESP_I2x2y2z_G2y2z = I_ESP_K2x2y3z_F2yz+ABZ*I_ESP_I2x2y2z_F2yz;
    Double I_ESP_I2xy3z_G2y2z = I_ESP_K2xy4z_F2yz+ABZ*I_ESP_I2xy3z_F2yz;
    Double I_ESP_I2x4z_G2y2z = I_ESP_K2x5z_F2yz+ABZ*I_ESP_I2x4z_F2yz;
    Double I_ESP_Ix5y_G2y2z = I_ESP_Kx5yz_F2yz+ABZ*I_ESP_Ix5y_F2yz;
    Double I_ESP_Ix4yz_G2y2z = I_ESP_Kx4y2z_F2yz+ABZ*I_ESP_Ix4yz_F2yz;
    Double I_ESP_Ix3y2z_G2y2z = I_ESP_Kx3y3z_F2yz+ABZ*I_ESP_Ix3y2z_F2yz;
    Double I_ESP_Ix2y3z_G2y2z = I_ESP_Kx2y4z_F2yz+ABZ*I_ESP_Ix2y3z_F2yz;
    Double I_ESP_Ixy4z_G2y2z = I_ESP_Kxy5z_F2yz+ABZ*I_ESP_Ixy4z_F2yz;
    Double I_ESP_Ix5z_G2y2z = I_ESP_Kx6z_F2yz+ABZ*I_ESP_Ix5z_F2yz;
    Double I_ESP_I5xy_Gy3z = I_ESP_K5x2y_F3z+ABY*I_ESP_I5xy_F3z;
    Double I_ESP_I4x2y_Gy3z = I_ESP_K4x3y_F3z+ABY*I_ESP_I4x2y_F3z;
    Double I_ESP_I4xyz_Gy3z = I_ESP_K4x2yz_F3z+ABY*I_ESP_I4xyz_F3z;
    Double I_ESP_I3x3y_Gy3z = I_ESP_K3x4y_F3z+ABY*I_ESP_I3x3y_F3z;
    Double I_ESP_I3x2yz_Gy3z = I_ESP_K3x3yz_F3z+ABY*I_ESP_I3x2yz_F3z;
    Double I_ESP_I3xy2z_Gy3z = I_ESP_K3x2y2z_F3z+ABY*I_ESP_I3xy2z_F3z;
    Double I_ESP_I2x4y_Gy3z = I_ESP_K2x5y_F3z+ABY*I_ESP_I2x4y_F3z;
    Double I_ESP_I2x3yz_Gy3z = I_ESP_K2x4yz_F3z+ABY*I_ESP_I2x3yz_F3z;
    Double I_ESP_I2x2y2z_Gy3z = I_ESP_K2x3y2z_F3z+ABY*I_ESP_I2x2y2z_F3z;
    Double I_ESP_I2xy3z_Gy3z = I_ESP_K2x2y3z_F3z+ABY*I_ESP_I2xy3z_F3z;
    Double I_ESP_Ix5y_Gy3z = I_ESP_Kx6y_F3z+ABY*I_ESP_Ix5y_F3z;
    Double I_ESP_Ix4yz_Gy3z = I_ESP_Kx5yz_F3z+ABY*I_ESP_Ix4yz_F3z;
    Double I_ESP_Ix3y2z_Gy3z = I_ESP_Kx4y2z_F3z+ABY*I_ESP_Ix3y2z_F3z;
    Double I_ESP_Ix2y3z_Gy3z = I_ESP_Kx3y3z_F3z+ABY*I_ESP_Ix2y3z_F3z;
    Double I_ESP_Ixy4z_Gy3z = I_ESP_Kx2y4z_F3z+ABY*I_ESP_Ixy4z_F3z;
    Double I_ESP_I6y_Gy3z = I_ESP_K7y_F3z+ABY*I_ESP_I6y_F3z;
    Double I_ESP_I5yz_Gy3z = I_ESP_K6yz_F3z+ABY*I_ESP_I5yz_F3z;
    Double I_ESP_I4y2z_Gy3z = I_ESP_K5y2z_F3z+ABY*I_ESP_I4y2z_F3z;
    Double I_ESP_I3y3z_Gy3z = I_ESP_K4y3z_F3z+ABY*I_ESP_I3y3z_F3z;
    Double I_ESP_I2y4z_Gy3z = I_ESP_K3y4z_F3z+ABY*I_ESP_I2y4z_F3z;
    Double I_ESP_Iy5z_Gy3z = I_ESP_K2y5z_F3z+ABY*I_ESP_Iy5z_F3z;
    Double I_ESP_I6x_G4z = I_ESP_K6xz_F3z+ABZ*I_ESP_I6x_F3z;
    Double I_ESP_I5xy_G4z = I_ESP_K5xyz_F3z+ABZ*I_ESP_I5xy_F3z;
    Double I_ESP_I5xz_G4z = I_ESP_K5x2z_F3z+ABZ*I_ESP_I5xz_F3z;
    Double I_ESP_I4x2y_G4z = I_ESP_K4x2yz_F3z+ABZ*I_ESP_I4x2y_F3z;
    Double I_ESP_I4xyz_G4z = I_ESP_K4xy2z_F3z+ABZ*I_ESP_I4xyz_F3z;
    Double I_ESP_I4x2z_G4z = I_ESP_K4x3z_F3z+ABZ*I_ESP_I4x2z_F3z;
    Double I_ESP_I3x3y_G4z = I_ESP_K3x3yz_F3z+ABZ*I_ESP_I3x3y_F3z;
    Double I_ESP_I3x2yz_G4z = I_ESP_K3x2y2z_F3z+ABZ*I_ESP_I3x2yz_F3z;
    Double I_ESP_I3xy2z_G4z = I_ESP_K3xy3z_F3z+ABZ*I_ESP_I3xy2z_F3z;
    Double I_ESP_I3x3z_G4z = I_ESP_K3x4z_F3z+ABZ*I_ESP_I3x3z_F3z;
    Double I_ESP_I2x4y_G4z = I_ESP_K2x4yz_F3z+ABZ*I_ESP_I2x4y_F3z;
    Double I_ESP_I2x3yz_G4z = I_ESP_K2x3y2z_F3z+ABZ*I_ESP_I2x3yz_F3z;
    Double I_ESP_I2x2y2z_G4z = I_ESP_K2x2y3z_F3z+ABZ*I_ESP_I2x2y2z_F3z;
    Double I_ESP_I2xy3z_G4z = I_ESP_K2xy4z_F3z+ABZ*I_ESP_I2xy3z_F3z;
    Double I_ESP_I2x4z_G4z = I_ESP_K2x5z_F3z+ABZ*I_ESP_I2x4z_F3z;
    Double I_ESP_Ix5y_G4z = I_ESP_Kx5yz_F3z+ABZ*I_ESP_Ix5y_F3z;
    Double I_ESP_Ix4yz_G4z = I_ESP_Kx4y2z_F3z+ABZ*I_ESP_Ix4yz_F3z;
    Double I_ESP_Ix3y2z_G4z = I_ESP_Kx3y3z_F3z+ABZ*I_ESP_Ix3y2z_F3z;
    Double I_ESP_Ix2y3z_G4z = I_ESP_Kx2y4z_F3z+ABZ*I_ESP_Ix2y3z_F3z;
    Double I_ESP_Ixy4z_G4z = I_ESP_Kxy5z_F3z+ABZ*I_ESP_Ixy4z_F3z;
    Double I_ESP_Ix5z_G4z = I_ESP_Kx6z_F3z+ABZ*I_ESP_Ix5z_F3z;
    Double I_ESP_I6y_G4z = I_ESP_K6yz_F3z+ABZ*I_ESP_I6y_F3z;
    Double I_ESP_I5yz_G4z = I_ESP_K5y2z_F3z+ABZ*I_ESP_I5yz_F3z;
    Double I_ESP_I4y2z_G4z = I_ESP_K4y3z_F3z+ABZ*I_ESP_I4y2z_F3z;
    Double I_ESP_I3y3z_G4z = I_ESP_K3y4z_F3z+ABZ*I_ESP_I3y3z_F3z;
    Double I_ESP_I2y4z_G4z = I_ESP_K2y5z_F3z+ABZ*I_ESP_I2y4z_F3z;
    Double I_ESP_Iy5z_G4z = I_ESP_Ky6z_F3z+ABZ*I_ESP_Iy5z_F3z;
    Double I_ESP_I6z_G4z = I_ESP_K7z_F3z+ABZ*I_ESP_I6z_F3z;

    /************************************************************
     * shell quartet name: SQ_ESP_H_H
     * expanding position: BRA2
     * code section is: HRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_G
     * RHS shell quartet name: SQ_ESP_H_G
     ************************************************************/
    abcd[iGrid*441+0] = I_ESP_I6x_G4x+ABX*I_ESP_H5x_G4x;
    abcd[iGrid*441+1] = I_ESP_I5xy_G4x+ABX*I_ESP_H4xy_G4x;
    abcd[iGrid*441+2] = I_ESP_I5xz_G4x+ABX*I_ESP_H4xz_G4x;
    abcd[iGrid*441+3] = I_ESP_I4x2y_G4x+ABX*I_ESP_H3x2y_G4x;
    abcd[iGrid*441+4] = I_ESP_I4xyz_G4x+ABX*I_ESP_H3xyz_G4x;
    abcd[iGrid*441+5] = I_ESP_I4x2z_G4x+ABX*I_ESP_H3x2z_G4x;
    abcd[iGrid*441+6] = I_ESP_I3x3y_G4x+ABX*I_ESP_H2x3y_G4x;
    abcd[iGrid*441+7] = I_ESP_I3x2yz_G4x+ABX*I_ESP_H2x2yz_G4x;
    abcd[iGrid*441+8] = I_ESP_I3xy2z_G4x+ABX*I_ESP_H2xy2z_G4x;
    abcd[iGrid*441+9] = I_ESP_I3x3z_G4x+ABX*I_ESP_H2x3z_G4x;
    abcd[iGrid*441+10] = I_ESP_I2x4y_G4x+ABX*I_ESP_Hx4y_G4x;
    abcd[iGrid*441+11] = I_ESP_I2x3yz_G4x+ABX*I_ESP_Hx3yz_G4x;
    abcd[iGrid*441+12] = I_ESP_I2x2y2z_G4x+ABX*I_ESP_Hx2y2z_G4x;
    abcd[iGrid*441+13] = I_ESP_I2xy3z_G4x+ABX*I_ESP_Hxy3z_G4x;
    abcd[iGrid*441+14] = I_ESP_I2x4z_G4x+ABX*I_ESP_Hx4z_G4x;
    abcd[iGrid*441+15] = I_ESP_Ix5y_G4x+ABX*I_ESP_H5y_G4x;
    abcd[iGrid*441+16] = I_ESP_Ix4yz_G4x+ABX*I_ESP_H4yz_G4x;
    abcd[iGrid*441+17] = I_ESP_Ix3y2z_G4x+ABX*I_ESP_H3y2z_G4x;
    abcd[iGrid*441+18] = I_ESP_Ix2y3z_G4x+ABX*I_ESP_H2y3z_G4x;
    abcd[iGrid*441+19] = I_ESP_Ixy4z_G4x+ABX*I_ESP_Hy4z_G4x;
    abcd[iGrid*441+20] = I_ESP_Ix5z_G4x+ABX*I_ESP_H5z_G4x;
    abcd[iGrid*441+21] = I_ESP_I5xy_G4x+ABY*I_ESP_H5x_G4x;
    abcd[iGrid*441+22] = I_ESP_I4x2y_G4x+ABY*I_ESP_H4xy_G4x;
    abcd[iGrid*441+23] = I_ESP_I4xyz_G4x+ABY*I_ESP_H4xz_G4x;
    abcd[iGrid*441+24] = I_ESP_I3x3y_G4x+ABY*I_ESP_H3x2y_G4x;
    abcd[iGrid*441+25] = I_ESP_I3x2yz_G4x+ABY*I_ESP_H3xyz_G4x;
    abcd[iGrid*441+26] = I_ESP_I3xy2z_G4x+ABY*I_ESP_H3x2z_G4x;
    abcd[iGrid*441+27] = I_ESP_I2x4y_G4x+ABY*I_ESP_H2x3y_G4x;
    abcd[iGrid*441+28] = I_ESP_I2x3yz_G4x+ABY*I_ESP_H2x2yz_G4x;
    abcd[iGrid*441+29] = I_ESP_I2x2y2z_G4x+ABY*I_ESP_H2xy2z_G4x;
    abcd[iGrid*441+30] = I_ESP_I2xy3z_G4x+ABY*I_ESP_H2x3z_G4x;
    abcd[iGrid*441+31] = I_ESP_Ix5y_G4x+ABY*I_ESP_Hx4y_G4x;
    abcd[iGrid*441+32] = I_ESP_Ix4yz_G4x+ABY*I_ESP_Hx3yz_G4x;
    abcd[iGrid*441+33] = I_ESP_Ix3y2z_G4x+ABY*I_ESP_Hx2y2z_G4x;
    abcd[iGrid*441+34] = I_ESP_Ix2y3z_G4x+ABY*I_ESP_Hxy3z_G4x;
    abcd[iGrid*441+35] = I_ESP_Ixy4z_G4x+ABY*I_ESP_Hx4z_G4x;
    abcd[iGrid*441+36] = I_ESP_I6y_G4x+ABY*I_ESP_H5y_G4x;
    abcd[iGrid*441+37] = I_ESP_I5yz_G4x+ABY*I_ESP_H4yz_G4x;
    abcd[iGrid*441+38] = I_ESP_I4y2z_G4x+ABY*I_ESP_H3y2z_G4x;
    abcd[iGrid*441+39] = I_ESP_I3y3z_G4x+ABY*I_ESP_H2y3z_G4x;
    abcd[iGrid*441+40] = I_ESP_I2y4z_G4x+ABY*I_ESP_Hy4z_G4x;
    abcd[iGrid*441+41] = I_ESP_Iy5z_G4x+ABY*I_ESP_H5z_G4x;
    abcd[iGrid*441+42] = I_ESP_I5xz_G4x+ABZ*I_ESP_H5x_G4x;
    abcd[iGrid*441+43] = I_ESP_I4xyz_G4x+ABZ*I_ESP_H4xy_G4x;
    abcd[iGrid*441+44] = I_ESP_I4x2z_G4x+ABZ*I_ESP_H4xz_G4x;
    abcd[iGrid*441+45] = I_ESP_I3x2yz_G4x+ABZ*I_ESP_H3x2y_G4x;
    abcd[iGrid*441+46] = I_ESP_I3xy2z_G4x+ABZ*I_ESP_H3xyz_G4x;
    abcd[iGrid*441+47] = I_ESP_I3x3z_G4x+ABZ*I_ESP_H3x2z_G4x;
    abcd[iGrid*441+48] = I_ESP_I2x3yz_G4x+ABZ*I_ESP_H2x3y_G4x;
    abcd[iGrid*441+49] = I_ESP_I2x2y2z_G4x+ABZ*I_ESP_H2x2yz_G4x;
    abcd[iGrid*441+50] = I_ESP_I2xy3z_G4x+ABZ*I_ESP_H2xy2z_G4x;
    abcd[iGrid*441+51] = I_ESP_I2x4z_G4x+ABZ*I_ESP_H2x3z_G4x;
    abcd[iGrid*441+52] = I_ESP_Ix4yz_G4x+ABZ*I_ESP_Hx4y_G4x;
    abcd[iGrid*441+53] = I_ESP_Ix3y2z_G4x+ABZ*I_ESP_Hx3yz_G4x;
    abcd[iGrid*441+54] = I_ESP_Ix2y3z_G4x+ABZ*I_ESP_Hx2y2z_G4x;
    abcd[iGrid*441+55] = I_ESP_Ixy4z_G4x+ABZ*I_ESP_Hxy3z_G4x;
    abcd[iGrid*441+56] = I_ESP_Ix5z_G4x+ABZ*I_ESP_Hx4z_G4x;
    abcd[iGrid*441+57] = I_ESP_I5yz_G4x+ABZ*I_ESP_H5y_G4x;
    abcd[iGrid*441+58] = I_ESP_I4y2z_G4x+ABZ*I_ESP_H4yz_G4x;
    abcd[iGrid*441+59] = I_ESP_I3y3z_G4x+ABZ*I_ESP_H3y2z_G4x;
    abcd[iGrid*441+60] = I_ESP_I2y4z_G4x+ABZ*I_ESP_H2y3z_G4x;
    abcd[iGrid*441+61] = I_ESP_Iy5z_G4x+ABZ*I_ESP_Hy4z_G4x;
    abcd[iGrid*441+62] = I_ESP_I6z_G4x+ABZ*I_ESP_H5z_G4x;
    abcd[iGrid*441+63] = I_ESP_I5xy_G3xy+ABY*I_ESP_H5x_G3xy;
    abcd[iGrid*441+64] = I_ESP_I4x2y_G3xy+ABY*I_ESP_H4xy_G3xy;
    abcd[iGrid*441+65] = I_ESP_I4xyz_G3xy+ABY*I_ESP_H4xz_G3xy;
    abcd[iGrid*441+66] = I_ESP_I3x3y_G3xy+ABY*I_ESP_H3x2y_G3xy;
    abcd[iGrid*441+67] = I_ESP_I3x2yz_G3xy+ABY*I_ESP_H3xyz_G3xy;
    abcd[iGrid*441+68] = I_ESP_I3xy2z_G3xy+ABY*I_ESP_H3x2z_G3xy;
    abcd[iGrid*441+69] = I_ESP_I2x4y_G3xy+ABY*I_ESP_H2x3y_G3xy;
    abcd[iGrid*441+70] = I_ESP_I2x3yz_G3xy+ABY*I_ESP_H2x2yz_G3xy;
    abcd[iGrid*441+71] = I_ESP_I2x2y2z_G3xy+ABY*I_ESP_H2xy2z_G3xy;
    abcd[iGrid*441+72] = I_ESP_I2xy3z_G3xy+ABY*I_ESP_H2x3z_G3xy;
    abcd[iGrid*441+73] = I_ESP_Ix5y_G3xy+ABY*I_ESP_Hx4y_G3xy;
    abcd[iGrid*441+74] = I_ESP_Ix4yz_G3xy+ABY*I_ESP_Hx3yz_G3xy;
    abcd[iGrid*441+75] = I_ESP_Ix3y2z_G3xy+ABY*I_ESP_Hx2y2z_G3xy;
    abcd[iGrid*441+76] = I_ESP_Ix2y3z_G3xy+ABY*I_ESP_Hxy3z_G3xy;
    abcd[iGrid*441+77] = I_ESP_Ixy4z_G3xy+ABY*I_ESP_Hx4z_G3xy;
    abcd[iGrid*441+78] = I_ESP_I6y_G3xy+ABY*I_ESP_H5y_G3xy;
    abcd[iGrid*441+79] = I_ESP_I5yz_G3xy+ABY*I_ESP_H4yz_G3xy;
    abcd[iGrid*441+80] = I_ESP_I4y2z_G3xy+ABY*I_ESP_H3y2z_G3xy;
    abcd[iGrid*441+81] = I_ESP_I3y3z_G3xy+ABY*I_ESP_H2y3z_G3xy;
    abcd[iGrid*441+82] = I_ESP_I2y4z_G3xy+ABY*I_ESP_Hy4z_G3xy;
    abcd[iGrid*441+83] = I_ESP_Iy5z_G3xy+ABY*I_ESP_H5z_G3xy;
    abcd[iGrid*441+84] = I_ESP_I5xz_G3xy+ABZ*I_ESP_H5x_G3xy;
    abcd[iGrid*441+85] = I_ESP_I4xyz_G3xy+ABZ*I_ESP_H4xy_G3xy;
    abcd[iGrid*441+86] = I_ESP_I4x2z_G3xy+ABZ*I_ESP_H4xz_G3xy;
    abcd[iGrid*441+87] = I_ESP_I3x2yz_G3xy+ABZ*I_ESP_H3x2y_G3xy;
    abcd[iGrid*441+88] = I_ESP_I3xy2z_G3xy+ABZ*I_ESP_H3xyz_G3xy;
    abcd[iGrid*441+89] = I_ESP_I3x3z_G3xy+ABZ*I_ESP_H3x2z_G3xy;
    abcd[iGrid*441+90] = I_ESP_I2x3yz_G3xy+ABZ*I_ESP_H2x3y_G3xy;
    abcd[iGrid*441+91] = I_ESP_I2x2y2z_G3xy+ABZ*I_ESP_H2x2yz_G3xy;
    abcd[iGrid*441+92] = I_ESP_I2xy3z_G3xy+ABZ*I_ESP_H2xy2z_G3xy;
    abcd[iGrid*441+93] = I_ESP_I2x4z_G3xy+ABZ*I_ESP_H2x3z_G3xy;
    abcd[iGrid*441+94] = I_ESP_Ix4yz_G3xy+ABZ*I_ESP_Hx4y_G3xy;
    abcd[iGrid*441+95] = I_ESP_Ix3y2z_G3xy+ABZ*I_ESP_Hx3yz_G3xy;
    abcd[iGrid*441+96] = I_ESP_Ix2y3z_G3xy+ABZ*I_ESP_Hx2y2z_G3xy;
    abcd[iGrid*441+97] = I_ESP_Ixy4z_G3xy+ABZ*I_ESP_Hxy3z_G3xy;
    abcd[iGrid*441+98] = I_ESP_Ix5z_G3xy+ABZ*I_ESP_Hx4z_G3xy;
    abcd[iGrid*441+99] = I_ESP_I5yz_G3xy+ABZ*I_ESP_H5y_G3xy;
    abcd[iGrid*441+100] = I_ESP_I4y2z_G3xy+ABZ*I_ESP_H4yz_G3xy;
    abcd[iGrid*441+101] = I_ESP_I3y3z_G3xy+ABZ*I_ESP_H3y2z_G3xy;
    abcd[iGrid*441+102] = I_ESP_I2y4z_G3xy+ABZ*I_ESP_H2y3z_G3xy;
    abcd[iGrid*441+103] = I_ESP_Iy5z_G3xy+ABZ*I_ESP_Hy4z_G3xy;
    abcd[iGrid*441+104] = I_ESP_I6z_G3xy+ABZ*I_ESP_H5z_G3xy;
    abcd[iGrid*441+105] = I_ESP_I5xz_G3xz+ABZ*I_ESP_H5x_G3xz;
    abcd[iGrid*441+106] = I_ESP_I4xyz_G3xz+ABZ*I_ESP_H4xy_G3xz;
    abcd[iGrid*441+107] = I_ESP_I4x2z_G3xz+ABZ*I_ESP_H4xz_G3xz;
    abcd[iGrid*441+108] = I_ESP_I3x2yz_G3xz+ABZ*I_ESP_H3x2y_G3xz;
    abcd[iGrid*441+109] = I_ESP_I3xy2z_G3xz+ABZ*I_ESP_H3xyz_G3xz;
    abcd[iGrid*441+110] = I_ESP_I3x3z_G3xz+ABZ*I_ESP_H3x2z_G3xz;
    abcd[iGrid*441+111] = I_ESP_I2x3yz_G3xz+ABZ*I_ESP_H2x3y_G3xz;
    abcd[iGrid*441+112] = I_ESP_I2x2y2z_G3xz+ABZ*I_ESP_H2x2yz_G3xz;
    abcd[iGrid*441+113] = I_ESP_I2xy3z_G3xz+ABZ*I_ESP_H2xy2z_G3xz;
    abcd[iGrid*441+114] = I_ESP_I2x4z_G3xz+ABZ*I_ESP_H2x3z_G3xz;
    abcd[iGrid*441+115] = I_ESP_Ix4yz_G3xz+ABZ*I_ESP_Hx4y_G3xz;
    abcd[iGrid*441+116] = I_ESP_Ix3y2z_G3xz+ABZ*I_ESP_Hx3yz_G3xz;
    abcd[iGrid*441+117] = I_ESP_Ix2y3z_G3xz+ABZ*I_ESP_Hx2y2z_G3xz;
    abcd[iGrid*441+118] = I_ESP_Ixy4z_G3xz+ABZ*I_ESP_Hxy3z_G3xz;
    abcd[iGrid*441+119] = I_ESP_Ix5z_G3xz+ABZ*I_ESP_Hx4z_G3xz;
    abcd[iGrid*441+120] = I_ESP_I5yz_G3xz+ABZ*I_ESP_H5y_G3xz;
    abcd[iGrid*441+121] = I_ESP_I4y2z_G3xz+ABZ*I_ESP_H4yz_G3xz;
    abcd[iGrid*441+122] = I_ESP_I3y3z_G3xz+ABZ*I_ESP_H3y2z_G3xz;
    abcd[iGrid*441+123] = I_ESP_I2y4z_G3xz+ABZ*I_ESP_H2y3z_G3xz;
    abcd[iGrid*441+124] = I_ESP_Iy5z_G3xz+ABZ*I_ESP_Hy4z_G3xz;
    abcd[iGrid*441+125] = I_ESP_I6z_G3xz+ABZ*I_ESP_H5z_G3xz;
    abcd[iGrid*441+126] = I_ESP_I6x_Gx3y+ABX*I_ESP_H5x_Gx3y;
    abcd[iGrid*441+127] = I_ESP_I5xy_Gx3y+ABX*I_ESP_H4xy_Gx3y;
    abcd[iGrid*441+128] = I_ESP_I5xz_Gx3y+ABX*I_ESP_H4xz_Gx3y;
    abcd[iGrid*441+129] = I_ESP_I4x2y_Gx3y+ABX*I_ESP_H3x2y_Gx3y;
    abcd[iGrid*441+130] = I_ESP_I4xyz_Gx3y+ABX*I_ESP_H3xyz_Gx3y;
    abcd[iGrid*441+131] = I_ESP_I4x2z_Gx3y+ABX*I_ESP_H3x2z_Gx3y;
    abcd[iGrid*441+132] = I_ESP_I3x3y_Gx3y+ABX*I_ESP_H2x3y_Gx3y;
    abcd[iGrid*441+133] = I_ESP_I3x2yz_Gx3y+ABX*I_ESP_H2x2yz_Gx3y;
    abcd[iGrid*441+134] = I_ESP_I3xy2z_Gx3y+ABX*I_ESP_H2xy2z_Gx3y;
    abcd[iGrid*441+135] = I_ESP_I3x3z_Gx3y+ABX*I_ESP_H2x3z_Gx3y;
    abcd[iGrid*441+136] = I_ESP_I2x4y_Gx3y+ABX*I_ESP_Hx4y_Gx3y;
    abcd[iGrid*441+137] = I_ESP_I2x3yz_Gx3y+ABX*I_ESP_Hx3yz_Gx3y;
    abcd[iGrid*441+138] = I_ESP_I2x2y2z_Gx3y+ABX*I_ESP_Hx2y2z_Gx3y;
    abcd[iGrid*441+139] = I_ESP_I2xy3z_Gx3y+ABX*I_ESP_Hxy3z_Gx3y;
    abcd[iGrid*441+140] = I_ESP_I2x4z_Gx3y+ABX*I_ESP_Hx4z_Gx3y;
    abcd[iGrid*441+141] = I_ESP_Ix5y_Gx3y+ABX*I_ESP_H5y_Gx3y;
    abcd[iGrid*441+142] = I_ESP_Ix4yz_Gx3y+ABX*I_ESP_H4yz_Gx3y;
    abcd[iGrid*441+143] = I_ESP_Ix3y2z_Gx3y+ABX*I_ESP_H3y2z_Gx3y;
    abcd[iGrid*441+144] = I_ESP_Ix2y3z_Gx3y+ABX*I_ESP_H2y3z_Gx3y;
    abcd[iGrid*441+145] = I_ESP_Ixy4z_Gx3y+ABX*I_ESP_Hy4z_Gx3y;
    abcd[iGrid*441+146] = I_ESP_Ix5z_Gx3y+ABX*I_ESP_H5z_Gx3y;
    abcd[iGrid*441+147] = I_ESP_I5xz_G2x2y+ABZ*I_ESP_H5x_G2x2y;
    abcd[iGrid*441+148] = I_ESP_I4xyz_G2x2y+ABZ*I_ESP_H4xy_G2x2y;
    abcd[iGrid*441+149] = I_ESP_I4x2z_G2x2y+ABZ*I_ESP_H4xz_G2x2y;
    abcd[iGrid*441+150] = I_ESP_I3x2yz_G2x2y+ABZ*I_ESP_H3x2y_G2x2y;
    abcd[iGrid*441+151] = I_ESP_I3xy2z_G2x2y+ABZ*I_ESP_H3xyz_G2x2y;
    abcd[iGrid*441+152] = I_ESP_I3x3z_G2x2y+ABZ*I_ESP_H3x2z_G2x2y;
    abcd[iGrid*441+153] = I_ESP_I2x3yz_G2x2y+ABZ*I_ESP_H2x3y_G2x2y;
    abcd[iGrid*441+154] = I_ESP_I2x2y2z_G2x2y+ABZ*I_ESP_H2x2yz_G2x2y;
    abcd[iGrid*441+155] = I_ESP_I2xy3z_G2x2y+ABZ*I_ESP_H2xy2z_G2x2y;
    abcd[iGrid*441+156] = I_ESP_I2x4z_G2x2y+ABZ*I_ESP_H2x3z_G2x2y;
    abcd[iGrid*441+157] = I_ESP_Ix4yz_G2x2y+ABZ*I_ESP_Hx4y_G2x2y;
    abcd[iGrid*441+158] = I_ESP_Ix3y2z_G2x2y+ABZ*I_ESP_Hx3yz_G2x2y;
    abcd[iGrid*441+159] = I_ESP_Ix2y3z_G2x2y+ABZ*I_ESP_Hx2y2z_G2x2y;
    abcd[iGrid*441+160] = I_ESP_Ixy4z_G2x2y+ABZ*I_ESP_Hxy3z_G2x2y;
    abcd[iGrid*441+161] = I_ESP_Ix5z_G2x2y+ABZ*I_ESP_Hx4z_G2x2y;
    abcd[iGrid*441+162] = I_ESP_I5yz_G2x2y+ABZ*I_ESP_H5y_G2x2y;
    abcd[iGrid*441+163] = I_ESP_I4y2z_G2x2y+ABZ*I_ESP_H4yz_G2x2y;
    abcd[iGrid*441+164] = I_ESP_I3y3z_G2x2y+ABZ*I_ESP_H3y2z_G2x2y;
    abcd[iGrid*441+165] = I_ESP_I2y4z_G2x2y+ABZ*I_ESP_H2y3z_G2x2y;
    abcd[iGrid*441+166] = I_ESP_Iy5z_G2x2y+ABZ*I_ESP_Hy4z_G2x2y;
    abcd[iGrid*441+167] = I_ESP_I6z_G2x2y+ABZ*I_ESP_H5z_G2x2y;
    abcd[iGrid*441+168] = I_ESP_I5xy_G2x2z+ABY*I_ESP_H5x_G2x2z;
    abcd[iGrid*441+169] = I_ESP_I4x2y_G2x2z+ABY*I_ESP_H4xy_G2x2z;
    abcd[iGrid*441+170] = I_ESP_I4xyz_G2x2z+ABY*I_ESP_H4xz_G2x2z;
    abcd[iGrid*441+171] = I_ESP_I3x3y_G2x2z+ABY*I_ESP_H3x2y_G2x2z;
    abcd[iGrid*441+172] = I_ESP_I3x2yz_G2x2z+ABY*I_ESP_H3xyz_G2x2z;
    abcd[iGrid*441+173] = I_ESP_I3xy2z_G2x2z+ABY*I_ESP_H3x2z_G2x2z;
    abcd[iGrid*441+174] = I_ESP_I2x4y_G2x2z+ABY*I_ESP_H2x3y_G2x2z;
    abcd[iGrid*441+175] = I_ESP_I2x3yz_G2x2z+ABY*I_ESP_H2x2yz_G2x2z;
    abcd[iGrid*441+176] = I_ESP_I2x2y2z_G2x2z+ABY*I_ESP_H2xy2z_G2x2z;
    abcd[iGrid*441+177] = I_ESP_I2xy3z_G2x2z+ABY*I_ESP_H2x3z_G2x2z;
    abcd[iGrid*441+178] = I_ESP_Ix5y_G2x2z+ABY*I_ESP_Hx4y_G2x2z;
    abcd[iGrid*441+179] = I_ESP_Ix4yz_G2x2z+ABY*I_ESP_Hx3yz_G2x2z;
    abcd[iGrid*441+180] = I_ESP_Ix3y2z_G2x2z+ABY*I_ESP_Hx2y2z_G2x2z;
    abcd[iGrid*441+181] = I_ESP_Ix2y3z_G2x2z+ABY*I_ESP_Hxy3z_G2x2z;
    abcd[iGrid*441+182] = I_ESP_Ixy4z_G2x2z+ABY*I_ESP_Hx4z_G2x2z;
    abcd[iGrid*441+183] = I_ESP_I6y_G2x2z+ABY*I_ESP_H5y_G2x2z;
    abcd[iGrid*441+184] = I_ESP_I5yz_G2x2z+ABY*I_ESP_H4yz_G2x2z;
    abcd[iGrid*441+185] = I_ESP_I4y2z_G2x2z+ABY*I_ESP_H3y2z_G2x2z;
    abcd[iGrid*441+186] = I_ESP_I3y3z_G2x2z+ABY*I_ESP_H2y3z_G2x2z;
    abcd[iGrid*441+187] = I_ESP_I2y4z_G2x2z+ABY*I_ESP_Hy4z_G2x2z;
    abcd[iGrid*441+188] = I_ESP_Iy5z_G2x2z+ABY*I_ESP_H5z_G2x2z;
    abcd[iGrid*441+189] = I_ESP_I6x_Gx3z+ABX*I_ESP_H5x_Gx3z;
    abcd[iGrid*441+190] = I_ESP_I5xy_Gx3z+ABX*I_ESP_H4xy_Gx3z;
    abcd[iGrid*441+191] = I_ESP_I5xz_Gx3z+ABX*I_ESP_H4xz_Gx3z;
    abcd[iGrid*441+192] = I_ESP_I4x2y_Gx3z+ABX*I_ESP_H3x2y_Gx3z;
    abcd[iGrid*441+193] = I_ESP_I4xyz_Gx3z+ABX*I_ESP_H3xyz_Gx3z;
    abcd[iGrid*441+194] = I_ESP_I4x2z_Gx3z+ABX*I_ESP_H3x2z_Gx3z;
    abcd[iGrid*441+195] = I_ESP_I3x3y_Gx3z+ABX*I_ESP_H2x3y_Gx3z;
    abcd[iGrid*441+196] = I_ESP_I3x2yz_Gx3z+ABX*I_ESP_H2x2yz_Gx3z;
    abcd[iGrid*441+197] = I_ESP_I3xy2z_Gx3z+ABX*I_ESP_H2xy2z_Gx3z;
    abcd[iGrid*441+198] = I_ESP_I3x3z_Gx3z+ABX*I_ESP_H2x3z_Gx3z;
    abcd[iGrid*441+199] = I_ESP_I2x4y_Gx3z+ABX*I_ESP_Hx4y_Gx3z;
    abcd[iGrid*441+200] = I_ESP_I2x3yz_Gx3z+ABX*I_ESP_Hx3yz_Gx3z;
    abcd[iGrid*441+201] = I_ESP_I2x2y2z_Gx3z+ABX*I_ESP_Hx2y2z_Gx3z;
    abcd[iGrid*441+202] = I_ESP_I2xy3z_Gx3z+ABX*I_ESP_Hxy3z_Gx3z;
    abcd[iGrid*441+203] = I_ESP_I2x4z_Gx3z+ABX*I_ESP_Hx4z_Gx3z;
    abcd[iGrid*441+204] = I_ESP_Ix5y_Gx3z+ABX*I_ESP_H5y_Gx3z;
    abcd[iGrid*441+205] = I_ESP_Ix4yz_Gx3z+ABX*I_ESP_H4yz_Gx3z;
    abcd[iGrid*441+206] = I_ESP_Ix3y2z_Gx3z+ABX*I_ESP_H3y2z_Gx3z;
    abcd[iGrid*441+207] = I_ESP_Ix2y3z_Gx3z+ABX*I_ESP_H2y3z_Gx3z;
    abcd[iGrid*441+208] = I_ESP_Ixy4z_Gx3z+ABX*I_ESP_Hy4z_Gx3z;
    abcd[iGrid*441+209] = I_ESP_Ix5z_Gx3z+ABX*I_ESP_H5z_Gx3z;
    abcd[iGrid*441+210] = I_ESP_I6x_G4y+ABX*I_ESP_H5x_G4y;
    abcd[iGrid*441+211] = I_ESP_I5xy_G4y+ABX*I_ESP_H4xy_G4y;
    abcd[iGrid*441+212] = I_ESP_I5xz_G4y+ABX*I_ESP_H4xz_G4y;
    abcd[iGrid*441+213] = I_ESP_I4x2y_G4y+ABX*I_ESP_H3x2y_G4y;
    abcd[iGrid*441+214] = I_ESP_I4xyz_G4y+ABX*I_ESP_H3xyz_G4y;
    abcd[iGrid*441+215] = I_ESP_I4x2z_G4y+ABX*I_ESP_H3x2z_G4y;
    abcd[iGrid*441+216] = I_ESP_I3x3y_G4y+ABX*I_ESP_H2x3y_G4y;
    abcd[iGrid*441+217] = I_ESP_I3x2yz_G4y+ABX*I_ESP_H2x2yz_G4y;
    abcd[iGrid*441+218] = I_ESP_I3xy2z_G4y+ABX*I_ESP_H2xy2z_G4y;
    abcd[iGrid*441+219] = I_ESP_I3x3z_G4y+ABX*I_ESP_H2x3z_G4y;
    abcd[iGrid*441+220] = I_ESP_I2x4y_G4y+ABX*I_ESP_Hx4y_G4y;
    abcd[iGrid*441+221] = I_ESP_I2x3yz_G4y+ABX*I_ESP_Hx3yz_G4y;
    abcd[iGrid*441+222] = I_ESP_I2x2y2z_G4y+ABX*I_ESP_Hx2y2z_G4y;
    abcd[iGrid*441+223] = I_ESP_I2xy3z_G4y+ABX*I_ESP_Hxy3z_G4y;
    abcd[iGrid*441+224] = I_ESP_I2x4z_G4y+ABX*I_ESP_Hx4z_G4y;
    abcd[iGrid*441+225] = I_ESP_Ix5y_G4y+ABX*I_ESP_H5y_G4y;
    abcd[iGrid*441+226] = I_ESP_Ix4yz_G4y+ABX*I_ESP_H4yz_G4y;
    abcd[iGrid*441+227] = I_ESP_Ix3y2z_G4y+ABX*I_ESP_H3y2z_G4y;
    abcd[iGrid*441+228] = I_ESP_Ix2y3z_G4y+ABX*I_ESP_H2y3z_G4y;
    abcd[iGrid*441+229] = I_ESP_Ixy4z_G4y+ABX*I_ESP_Hy4z_G4y;
    abcd[iGrid*441+230] = I_ESP_Ix5z_G4y+ABX*I_ESP_H5z_G4y;
    abcd[iGrid*441+231] = I_ESP_I5xz_Gx3y+ABZ*I_ESP_H5x_Gx3y;
    abcd[iGrid*441+232] = I_ESP_I4xyz_Gx3y+ABZ*I_ESP_H4xy_Gx3y;
    abcd[iGrid*441+233] = I_ESP_I4x2z_Gx3y+ABZ*I_ESP_H4xz_Gx3y;
    abcd[iGrid*441+234] = I_ESP_I3x2yz_Gx3y+ABZ*I_ESP_H3x2y_Gx3y;
    abcd[iGrid*441+235] = I_ESP_I3xy2z_Gx3y+ABZ*I_ESP_H3xyz_Gx3y;
    abcd[iGrid*441+236] = I_ESP_I3x3z_Gx3y+ABZ*I_ESP_H3x2z_Gx3y;
    abcd[iGrid*441+237] = I_ESP_I2x3yz_Gx3y+ABZ*I_ESP_H2x3y_Gx3y;
    abcd[iGrid*441+238] = I_ESP_I2x2y2z_Gx3y+ABZ*I_ESP_H2x2yz_Gx3y;
    abcd[iGrid*441+239] = I_ESP_I2xy3z_Gx3y+ABZ*I_ESP_H2xy2z_Gx3y;
    abcd[iGrid*441+240] = I_ESP_I2x4z_Gx3y+ABZ*I_ESP_H2x3z_Gx3y;
    abcd[iGrid*441+241] = I_ESP_Ix4yz_Gx3y+ABZ*I_ESP_Hx4y_Gx3y;
    abcd[iGrid*441+242] = I_ESP_Ix3y2z_Gx3y+ABZ*I_ESP_Hx3yz_Gx3y;
    abcd[iGrid*441+243] = I_ESP_Ix2y3z_Gx3y+ABZ*I_ESP_Hx2y2z_Gx3y;
    abcd[iGrid*441+244] = I_ESP_Ixy4z_Gx3y+ABZ*I_ESP_Hxy3z_Gx3y;
    abcd[iGrid*441+245] = I_ESP_Ix5z_Gx3y+ABZ*I_ESP_Hx4z_Gx3y;
    abcd[iGrid*441+246] = I_ESP_I5yz_Gx3y+ABZ*I_ESP_H5y_Gx3y;
    abcd[iGrid*441+247] = I_ESP_I4y2z_Gx3y+ABZ*I_ESP_H4yz_Gx3y;
    abcd[iGrid*441+248] = I_ESP_I3y3z_Gx3y+ABZ*I_ESP_H3y2z_Gx3y;
    abcd[iGrid*441+249] = I_ESP_I2y4z_Gx3y+ABZ*I_ESP_H2y3z_Gx3y;
    abcd[iGrid*441+250] = I_ESP_Iy5z_Gx3y+ABZ*I_ESP_Hy4z_Gx3y;
    abcd[iGrid*441+251] = I_ESP_I6z_Gx3y+ABZ*I_ESP_H5z_Gx3y;
    abcd[iGrid*441+252] = I_ESP_I6x_G2y2z+ABX*I_ESP_H5x_G2y2z;
    abcd[iGrid*441+253] = I_ESP_I5xy_G2y2z+ABX*I_ESP_H4xy_G2y2z;
    abcd[iGrid*441+254] = I_ESP_I5xz_G2y2z+ABX*I_ESP_H4xz_G2y2z;
    abcd[iGrid*441+255] = I_ESP_I4x2y_G2y2z+ABX*I_ESP_H3x2y_G2y2z;
    abcd[iGrid*441+256] = I_ESP_I4xyz_G2y2z+ABX*I_ESP_H3xyz_G2y2z;
    abcd[iGrid*441+257] = I_ESP_I4x2z_G2y2z+ABX*I_ESP_H3x2z_G2y2z;
    abcd[iGrid*441+258] = I_ESP_I3x3y_G2y2z+ABX*I_ESP_H2x3y_G2y2z;
    abcd[iGrid*441+259] = I_ESP_I3x2yz_G2y2z+ABX*I_ESP_H2x2yz_G2y2z;
    abcd[iGrid*441+260] = I_ESP_I3xy2z_G2y2z+ABX*I_ESP_H2xy2z_G2y2z;
    abcd[iGrid*441+261] = I_ESP_I3x3z_G2y2z+ABX*I_ESP_H2x3z_G2y2z;
    abcd[iGrid*441+262] = I_ESP_I2x4y_G2y2z+ABX*I_ESP_Hx4y_G2y2z;
    abcd[iGrid*441+263] = I_ESP_I2x3yz_G2y2z+ABX*I_ESP_Hx3yz_G2y2z;
    abcd[iGrid*441+264] = I_ESP_I2x2y2z_G2y2z+ABX*I_ESP_Hx2y2z_G2y2z;
    abcd[iGrid*441+265] = I_ESP_I2xy3z_G2y2z+ABX*I_ESP_Hxy3z_G2y2z;
    abcd[iGrid*441+266] = I_ESP_I2x4z_G2y2z+ABX*I_ESP_Hx4z_G2y2z;
    abcd[iGrid*441+267] = I_ESP_Ix5y_G2y2z+ABX*I_ESP_H5y_G2y2z;
    abcd[iGrid*441+268] = I_ESP_Ix4yz_G2y2z+ABX*I_ESP_H4yz_G2y2z;
    abcd[iGrid*441+269] = I_ESP_Ix3y2z_G2y2z+ABX*I_ESP_H3y2z_G2y2z;
    abcd[iGrid*441+270] = I_ESP_Ix2y3z_G2y2z+ABX*I_ESP_H2y3z_G2y2z;
    abcd[iGrid*441+271] = I_ESP_Ixy4z_G2y2z+ABX*I_ESP_Hy4z_G2y2z;
    abcd[iGrid*441+272] = I_ESP_Ix5z_G2y2z+ABX*I_ESP_H5z_G2y2z;
    abcd[iGrid*441+273] = I_ESP_I5xy_Gx3z+ABY*I_ESP_H5x_Gx3z;
    abcd[iGrid*441+274] = I_ESP_I4x2y_Gx3z+ABY*I_ESP_H4xy_Gx3z;
    abcd[iGrid*441+275] = I_ESP_I4xyz_Gx3z+ABY*I_ESP_H4xz_Gx3z;
    abcd[iGrid*441+276] = I_ESP_I3x3y_Gx3z+ABY*I_ESP_H3x2y_Gx3z;
    abcd[iGrid*441+277] = I_ESP_I3x2yz_Gx3z+ABY*I_ESP_H3xyz_Gx3z;
    abcd[iGrid*441+278] = I_ESP_I3xy2z_Gx3z+ABY*I_ESP_H3x2z_Gx3z;
    abcd[iGrid*441+279] = I_ESP_I2x4y_Gx3z+ABY*I_ESP_H2x3y_Gx3z;
    abcd[iGrid*441+280] = I_ESP_I2x3yz_Gx3z+ABY*I_ESP_H2x2yz_Gx3z;
    abcd[iGrid*441+281] = I_ESP_I2x2y2z_Gx3z+ABY*I_ESP_H2xy2z_Gx3z;
    abcd[iGrid*441+282] = I_ESP_I2xy3z_Gx3z+ABY*I_ESP_H2x3z_Gx3z;
    abcd[iGrid*441+283] = I_ESP_Ix5y_Gx3z+ABY*I_ESP_Hx4y_Gx3z;
    abcd[iGrid*441+284] = I_ESP_Ix4yz_Gx3z+ABY*I_ESP_Hx3yz_Gx3z;
    abcd[iGrid*441+285] = I_ESP_Ix3y2z_Gx3z+ABY*I_ESP_Hx2y2z_Gx3z;
    abcd[iGrid*441+286] = I_ESP_Ix2y3z_Gx3z+ABY*I_ESP_Hxy3z_Gx3z;
    abcd[iGrid*441+287] = I_ESP_Ixy4z_Gx3z+ABY*I_ESP_Hx4z_Gx3z;
    abcd[iGrid*441+288] = I_ESP_I6y_Gx3z+ABY*I_ESP_H5y_Gx3z;
    abcd[iGrid*441+289] = I_ESP_I5yz_Gx3z+ABY*I_ESP_H4yz_Gx3z;
    abcd[iGrid*441+290] = I_ESP_I4y2z_Gx3z+ABY*I_ESP_H3y2z_Gx3z;
    abcd[iGrid*441+291] = I_ESP_I3y3z_Gx3z+ABY*I_ESP_H2y3z_Gx3z;
    abcd[iGrid*441+292] = I_ESP_I2y4z_Gx3z+ABY*I_ESP_Hy4z_Gx3z;
    abcd[iGrid*441+293] = I_ESP_Iy5z_Gx3z+ABY*I_ESP_H5z_Gx3z;
    abcd[iGrid*441+294] = I_ESP_I6x_G4z+ABX*I_ESP_H5x_G4z;
    abcd[iGrid*441+295] = I_ESP_I5xy_G4z+ABX*I_ESP_H4xy_G4z;
    abcd[iGrid*441+296] = I_ESP_I5xz_G4z+ABX*I_ESP_H4xz_G4z;
    abcd[iGrid*441+297] = I_ESP_I4x2y_G4z+ABX*I_ESP_H3x2y_G4z;
    abcd[iGrid*441+298] = I_ESP_I4xyz_G4z+ABX*I_ESP_H3xyz_G4z;
    abcd[iGrid*441+299] = I_ESP_I4x2z_G4z+ABX*I_ESP_H3x2z_G4z;
    abcd[iGrid*441+300] = I_ESP_I3x3y_G4z+ABX*I_ESP_H2x3y_G4z;
    abcd[iGrid*441+301] = I_ESP_I3x2yz_G4z+ABX*I_ESP_H2x2yz_G4z;
    abcd[iGrid*441+302] = I_ESP_I3xy2z_G4z+ABX*I_ESP_H2xy2z_G4z;
    abcd[iGrid*441+303] = I_ESP_I3x3z_G4z+ABX*I_ESP_H2x3z_G4z;
    abcd[iGrid*441+304] = I_ESP_I2x4y_G4z+ABX*I_ESP_Hx4y_G4z;
    abcd[iGrid*441+305] = I_ESP_I2x3yz_G4z+ABX*I_ESP_Hx3yz_G4z;
    abcd[iGrid*441+306] = I_ESP_I2x2y2z_G4z+ABX*I_ESP_Hx2y2z_G4z;
    abcd[iGrid*441+307] = I_ESP_I2xy3z_G4z+ABX*I_ESP_Hxy3z_G4z;
    abcd[iGrid*441+308] = I_ESP_I2x4z_G4z+ABX*I_ESP_Hx4z_G4z;
    abcd[iGrid*441+309] = I_ESP_Ix5y_G4z+ABX*I_ESP_H5y_G4z;
    abcd[iGrid*441+310] = I_ESP_Ix4yz_G4z+ABX*I_ESP_H4yz_G4z;
    abcd[iGrid*441+311] = I_ESP_Ix3y2z_G4z+ABX*I_ESP_H3y2z_G4z;
    abcd[iGrid*441+312] = I_ESP_Ix2y3z_G4z+ABX*I_ESP_H2y3z_G4z;
    abcd[iGrid*441+313] = I_ESP_Ixy4z_G4z+ABX*I_ESP_Hy4z_G4z;
    abcd[iGrid*441+314] = I_ESP_Ix5z_G4z+ABX*I_ESP_H5z_G4z;
    abcd[iGrid*441+315] = I_ESP_I5xy_G4y+ABY*I_ESP_H5x_G4y;
    abcd[iGrid*441+316] = I_ESP_I4x2y_G4y+ABY*I_ESP_H4xy_G4y;
    abcd[iGrid*441+317] = I_ESP_I4xyz_G4y+ABY*I_ESP_H4xz_G4y;
    abcd[iGrid*441+318] = I_ESP_I3x3y_G4y+ABY*I_ESP_H3x2y_G4y;
    abcd[iGrid*441+319] = I_ESP_I3x2yz_G4y+ABY*I_ESP_H3xyz_G4y;
    abcd[iGrid*441+320] = I_ESP_I3xy2z_G4y+ABY*I_ESP_H3x2z_G4y;
    abcd[iGrid*441+321] = I_ESP_I2x4y_G4y+ABY*I_ESP_H2x3y_G4y;
    abcd[iGrid*441+322] = I_ESP_I2x3yz_G4y+ABY*I_ESP_H2x2yz_G4y;
    abcd[iGrid*441+323] = I_ESP_I2x2y2z_G4y+ABY*I_ESP_H2xy2z_G4y;
    abcd[iGrid*441+324] = I_ESP_I2xy3z_G4y+ABY*I_ESP_H2x3z_G4y;
    abcd[iGrid*441+325] = I_ESP_Ix5y_G4y+ABY*I_ESP_Hx4y_G4y;
    abcd[iGrid*441+326] = I_ESP_Ix4yz_G4y+ABY*I_ESP_Hx3yz_G4y;
    abcd[iGrid*441+327] = I_ESP_Ix3y2z_G4y+ABY*I_ESP_Hx2y2z_G4y;
    abcd[iGrid*441+328] = I_ESP_Ix2y3z_G4y+ABY*I_ESP_Hxy3z_G4y;
    abcd[iGrid*441+329] = I_ESP_Ixy4z_G4y+ABY*I_ESP_Hx4z_G4y;
    abcd[iGrid*441+330] = I_ESP_I6y_G4y+ABY*I_ESP_H5y_G4y;
    abcd[iGrid*441+331] = I_ESP_I5yz_G4y+ABY*I_ESP_H4yz_G4y;
    abcd[iGrid*441+332] = I_ESP_I4y2z_G4y+ABY*I_ESP_H3y2z_G4y;
    abcd[iGrid*441+333] = I_ESP_I3y3z_G4y+ABY*I_ESP_H2y3z_G4y;
    abcd[iGrid*441+334] = I_ESP_I2y4z_G4y+ABY*I_ESP_Hy4z_G4y;
    abcd[iGrid*441+335] = I_ESP_Iy5z_G4y+ABY*I_ESP_H5z_G4y;
    abcd[iGrid*441+336] = I_ESP_I5xz_G4y+ABZ*I_ESP_H5x_G4y;
    abcd[iGrid*441+337] = I_ESP_I4xyz_G4y+ABZ*I_ESP_H4xy_G4y;
    abcd[iGrid*441+338] = I_ESP_I4x2z_G4y+ABZ*I_ESP_H4xz_G4y;
    abcd[iGrid*441+339] = I_ESP_I3x2yz_G4y+ABZ*I_ESP_H3x2y_G4y;
    abcd[iGrid*441+340] = I_ESP_I3xy2z_G4y+ABZ*I_ESP_H3xyz_G4y;
    abcd[iGrid*441+341] = I_ESP_I3x3z_G4y+ABZ*I_ESP_H3x2z_G4y;
    abcd[iGrid*441+342] = I_ESP_I2x3yz_G4y+ABZ*I_ESP_H2x3y_G4y;
    abcd[iGrid*441+343] = I_ESP_I2x2y2z_G4y+ABZ*I_ESP_H2x2yz_G4y;
    abcd[iGrid*441+344] = I_ESP_I2xy3z_G4y+ABZ*I_ESP_H2xy2z_G4y;
    abcd[iGrid*441+345] = I_ESP_I2x4z_G4y+ABZ*I_ESP_H2x3z_G4y;
    abcd[iGrid*441+346] = I_ESP_Ix4yz_G4y+ABZ*I_ESP_Hx4y_G4y;
    abcd[iGrid*441+347] = I_ESP_Ix3y2z_G4y+ABZ*I_ESP_Hx3yz_G4y;
    abcd[iGrid*441+348] = I_ESP_Ix2y3z_G4y+ABZ*I_ESP_Hx2y2z_G4y;
    abcd[iGrid*441+349] = I_ESP_Ixy4z_G4y+ABZ*I_ESP_Hxy3z_G4y;
    abcd[iGrid*441+350] = I_ESP_Ix5z_G4y+ABZ*I_ESP_Hx4z_G4y;
    abcd[iGrid*441+351] = I_ESP_I5yz_G4y+ABZ*I_ESP_H5y_G4y;
    abcd[iGrid*441+352] = I_ESP_I4y2z_G4y+ABZ*I_ESP_H4yz_G4y;
    abcd[iGrid*441+353] = I_ESP_I3y3z_G4y+ABZ*I_ESP_H3y2z_G4y;
    abcd[iGrid*441+354] = I_ESP_I2y4z_G4y+ABZ*I_ESP_H2y3z_G4y;
    abcd[iGrid*441+355] = I_ESP_Iy5z_G4y+ABZ*I_ESP_Hy4z_G4y;
    abcd[iGrid*441+356] = I_ESP_I6z_G4y+ABZ*I_ESP_H5z_G4y;
    abcd[iGrid*441+357] = I_ESP_I5xz_G3yz+ABZ*I_ESP_H5x_G3yz;
    abcd[iGrid*441+358] = I_ESP_I4xyz_G3yz+ABZ*I_ESP_H4xy_G3yz;
    abcd[iGrid*441+359] = I_ESP_I4x2z_G3yz+ABZ*I_ESP_H4xz_G3yz;
    abcd[iGrid*441+360] = I_ESP_I3x2yz_G3yz+ABZ*I_ESP_H3x2y_G3yz;
    abcd[iGrid*441+361] = I_ESP_I3xy2z_G3yz+ABZ*I_ESP_H3xyz_G3yz;
    abcd[iGrid*441+362] = I_ESP_I3x3z_G3yz+ABZ*I_ESP_H3x2z_G3yz;
    abcd[iGrid*441+363] = I_ESP_I2x3yz_G3yz+ABZ*I_ESP_H2x3y_G3yz;
    abcd[iGrid*441+364] = I_ESP_I2x2y2z_G3yz+ABZ*I_ESP_H2x2yz_G3yz;
    abcd[iGrid*441+365] = I_ESP_I2xy3z_G3yz+ABZ*I_ESP_H2xy2z_G3yz;
    abcd[iGrid*441+366] = I_ESP_I2x4z_G3yz+ABZ*I_ESP_H2x3z_G3yz;
    abcd[iGrid*441+367] = I_ESP_Ix4yz_G3yz+ABZ*I_ESP_Hx4y_G3yz;
    abcd[iGrid*441+368] = I_ESP_Ix3y2z_G3yz+ABZ*I_ESP_Hx3yz_G3yz;
    abcd[iGrid*441+369] = I_ESP_Ix2y3z_G3yz+ABZ*I_ESP_Hx2y2z_G3yz;
    abcd[iGrid*441+370] = I_ESP_Ixy4z_G3yz+ABZ*I_ESP_Hxy3z_G3yz;
    abcd[iGrid*441+371] = I_ESP_Ix5z_G3yz+ABZ*I_ESP_Hx4z_G3yz;
    abcd[iGrid*441+372] = I_ESP_I5yz_G3yz+ABZ*I_ESP_H5y_G3yz;
    abcd[iGrid*441+373] = I_ESP_I4y2z_G3yz+ABZ*I_ESP_H4yz_G3yz;
    abcd[iGrid*441+374] = I_ESP_I3y3z_G3yz+ABZ*I_ESP_H3y2z_G3yz;
    abcd[iGrid*441+375] = I_ESP_I2y4z_G3yz+ABZ*I_ESP_H2y3z_G3yz;
    abcd[iGrid*441+376] = I_ESP_Iy5z_G3yz+ABZ*I_ESP_Hy4z_G3yz;
    abcd[iGrid*441+377] = I_ESP_I6z_G3yz+ABZ*I_ESP_H5z_G3yz;
    abcd[iGrid*441+378] = I_ESP_I5xy_Gy3z+ABY*I_ESP_H5x_Gy3z;
    abcd[iGrid*441+379] = I_ESP_I4x2y_Gy3z+ABY*I_ESP_H4xy_Gy3z;
    abcd[iGrid*441+380] = I_ESP_I4xyz_Gy3z+ABY*I_ESP_H4xz_Gy3z;
    abcd[iGrid*441+381] = I_ESP_I3x3y_Gy3z+ABY*I_ESP_H3x2y_Gy3z;
    abcd[iGrid*441+382] = I_ESP_I3x2yz_Gy3z+ABY*I_ESP_H3xyz_Gy3z;
    abcd[iGrid*441+383] = I_ESP_I3xy2z_Gy3z+ABY*I_ESP_H3x2z_Gy3z;
    abcd[iGrid*441+384] = I_ESP_I2x4y_Gy3z+ABY*I_ESP_H2x3y_Gy3z;
    abcd[iGrid*441+385] = I_ESP_I2x3yz_Gy3z+ABY*I_ESP_H2x2yz_Gy3z;
    abcd[iGrid*441+386] = I_ESP_I2x2y2z_Gy3z+ABY*I_ESP_H2xy2z_Gy3z;
    abcd[iGrid*441+387] = I_ESP_I2xy3z_Gy3z+ABY*I_ESP_H2x3z_Gy3z;
    abcd[iGrid*441+388] = I_ESP_Ix5y_Gy3z+ABY*I_ESP_Hx4y_Gy3z;
    abcd[iGrid*441+389] = I_ESP_Ix4yz_Gy3z+ABY*I_ESP_Hx3yz_Gy3z;
    abcd[iGrid*441+390] = I_ESP_Ix3y2z_Gy3z+ABY*I_ESP_Hx2y2z_Gy3z;
    abcd[iGrid*441+391] = I_ESP_Ix2y3z_Gy3z+ABY*I_ESP_Hxy3z_Gy3z;
    abcd[iGrid*441+392] = I_ESP_Ixy4z_Gy3z+ABY*I_ESP_Hx4z_Gy3z;
    abcd[iGrid*441+393] = I_ESP_I6y_Gy3z+ABY*I_ESP_H5y_Gy3z;
    abcd[iGrid*441+394] = I_ESP_I5yz_Gy3z+ABY*I_ESP_H4yz_Gy3z;
    abcd[iGrid*441+395] = I_ESP_I4y2z_Gy3z+ABY*I_ESP_H3y2z_Gy3z;
    abcd[iGrid*441+396] = I_ESP_I3y3z_Gy3z+ABY*I_ESP_H2y3z_Gy3z;
    abcd[iGrid*441+397] = I_ESP_I2y4z_Gy3z+ABY*I_ESP_Hy4z_Gy3z;
    abcd[iGrid*441+398] = I_ESP_Iy5z_Gy3z+ABY*I_ESP_H5z_Gy3z;
    abcd[iGrid*441+399] = I_ESP_I5xy_G4z+ABY*I_ESP_H5x_G4z;
    abcd[iGrid*441+400] = I_ESP_I4x2y_G4z+ABY*I_ESP_H4xy_G4z;
    abcd[iGrid*441+401] = I_ESP_I4xyz_G4z+ABY*I_ESP_H4xz_G4z;
    abcd[iGrid*441+402] = I_ESP_I3x3y_G4z+ABY*I_ESP_H3x2y_G4z;
    abcd[iGrid*441+403] = I_ESP_I3x2yz_G4z+ABY*I_ESP_H3xyz_G4z;
    abcd[iGrid*441+404] = I_ESP_I3xy2z_G4z+ABY*I_ESP_H3x2z_G4z;
    abcd[iGrid*441+405] = I_ESP_I2x4y_G4z+ABY*I_ESP_H2x3y_G4z;
    abcd[iGrid*441+406] = I_ESP_I2x3yz_G4z+ABY*I_ESP_H2x2yz_G4z;
    abcd[iGrid*441+407] = I_ESP_I2x2y2z_G4z+ABY*I_ESP_H2xy2z_G4z;
    abcd[iGrid*441+408] = I_ESP_I2xy3z_G4z+ABY*I_ESP_H2x3z_G4z;
    abcd[iGrid*441+409] = I_ESP_Ix5y_G4z+ABY*I_ESP_Hx4y_G4z;
    abcd[iGrid*441+410] = I_ESP_Ix4yz_G4z+ABY*I_ESP_Hx3yz_G4z;
    abcd[iGrid*441+411] = I_ESP_Ix3y2z_G4z+ABY*I_ESP_Hx2y2z_G4z;
    abcd[iGrid*441+412] = I_ESP_Ix2y3z_G4z+ABY*I_ESP_Hxy3z_G4z;
    abcd[iGrid*441+413] = I_ESP_Ixy4z_G4z+ABY*I_ESP_Hx4z_G4z;
    abcd[iGrid*441+414] = I_ESP_I6y_G4z+ABY*I_ESP_H5y_G4z;
    abcd[iGrid*441+415] = I_ESP_I5yz_G4z+ABY*I_ESP_H4yz_G4z;
    abcd[iGrid*441+416] = I_ESP_I4y2z_G4z+ABY*I_ESP_H3y2z_G4z;
    abcd[iGrid*441+417] = I_ESP_I3y3z_G4z+ABY*I_ESP_H2y3z_G4z;
    abcd[iGrid*441+418] = I_ESP_I2y4z_G4z+ABY*I_ESP_Hy4z_G4z;
    abcd[iGrid*441+419] = I_ESP_Iy5z_G4z+ABY*I_ESP_H5z_G4z;
    abcd[iGrid*441+420] = I_ESP_I5xz_G4z+ABZ*I_ESP_H5x_G4z;
    abcd[iGrid*441+421] = I_ESP_I4xyz_G4z+ABZ*I_ESP_H4xy_G4z;
    abcd[iGrid*441+422] = I_ESP_I4x2z_G4z+ABZ*I_ESP_H4xz_G4z;
    abcd[iGrid*441+423] = I_ESP_I3x2yz_G4z+ABZ*I_ESP_H3x2y_G4z;
    abcd[iGrid*441+424] = I_ESP_I3xy2z_G4z+ABZ*I_ESP_H3xyz_G4z;
    abcd[iGrid*441+425] = I_ESP_I3x3z_G4z+ABZ*I_ESP_H3x2z_G4z;
    abcd[iGrid*441+426] = I_ESP_I2x3yz_G4z+ABZ*I_ESP_H2x3y_G4z;
    abcd[iGrid*441+427] = I_ESP_I2x2y2z_G4z+ABZ*I_ESP_H2x2yz_G4z;
    abcd[iGrid*441+428] = I_ESP_I2xy3z_G4z+ABZ*I_ESP_H2xy2z_G4z;
    abcd[iGrid*441+429] = I_ESP_I2x4z_G4z+ABZ*I_ESP_H2x3z_G4z;
    abcd[iGrid*441+430] = I_ESP_Ix4yz_G4z+ABZ*I_ESP_Hx4y_G4z;
    abcd[iGrid*441+431] = I_ESP_Ix3y2z_G4z+ABZ*I_ESP_Hx3yz_G4z;
    abcd[iGrid*441+432] = I_ESP_Ix2y3z_G4z+ABZ*I_ESP_Hx2y2z_G4z;
    abcd[iGrid*441+433] = I_ESP_Ixy4z_G4z+ABZ*I_ESP_Hxy3z_G4z;
    abcd[iGrid*441+434] = I_ESP_Ix5z_G4z+ABZ*I_ESP_Hx4z_G4z;
    abcd[iGrid*441+435] = I_ESP_I5yz_G4z+ABZ*I_ESP_H5y_G4z;
    abcd[iGrid*441+436] = I_ESP_I4y2z_G4z+ABZ*I_ESP_H4yz_G4z;
    abcd[iGrid*441+437] = I_ESP_I3y3z_G4z+ABZ*I_ESP_H3y2z_G4z;
    abcd[iGrid*441+438] = I_ESP_I2y4z_G4z+ABZ*I_ESP_H2y3z_G4z;
    abcd[iGrid*441+439] = I_ESP_Iy5z_G4z+ABZ*I_ESP_Hy4z_G4z;
    abcd[iGrid*441+440] = I_ESP_I6z_G4z+ABZ*I_ESP_H5z_G4z;
  }
}
