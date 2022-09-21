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
// BRA1 as redundant position, total RHS integrals evaluated as: 5680
// BRA2 as redundant position, total RHS integrals evaluated as: 5260
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

void hgp_os_esp_h_f_d1(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
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
     * totally 30 integrals are omitted 
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
    Double I_ESP_G4x_Dxy = I_ESP_H4xy_Px+ABY*I_ESP_G4x_Px;
    Double I_ESP_G3xy_Dxy = I_ESP_H3x2y_Px+ABY*I_ESP_G3xy_Px;
    Double I_ESP_G3xz_Dxy = I_ESP_H3xyz_Px+ABY*I_ESP_G3xz_Px;
    Double I_ESP_G2x2y_Dxy = I_ESP_H2x3y_Px+ABY*I_ESP_G2x2y_Px;
    Double I_ESP_G2xyz_Dxy = I_ESP_H2x2yz_Px+ABY*I_ESP_G2xyz_Px;
    Double I_ESP_G2x2z_Dxy = I_ESP_H2xy2z_Px+ABY*I_ESP_G2x2z_Px;
    Double I_ESP_Gx3y_Dxy = I_ESP_Hx4y_Px+ABY*I_ESP_Gx3y_Px;
    Double I_ESP_Gx2yz_Dxy = I_ESP_Hx3yz_Px+ABY*I_ESP_Gx2yz_Px;
    Double I_ESP_Gxy2z_Dxy = I_ESP_Hx2y2z_Px+ABY*I_ESP_Gxy2z_Px;
    Double I_ESP_Gx3z_Dxy = I_ESP_Hxy3z_Px+ABY*I_ESP_Gx3z_Px;
    Double I_ESP_G4y_Dxy = I_ESP_H5y_Px+ABY*I_ESP_G4y_Px;
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
     * shell quartet name: SQ_ESP_I_P
     * expanding position: BRA2
     * code section is: HRR
     * totally 16 integrals are omitted 
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
    Double I_ESP_I5yz_Px = I_ESP_Kx5yz_S+ABX*I_ESP_I5yz_S;
    Double I_ESP_I4y2z_Px = I_ESP_Kx4y2z_S+ABX*I_ESP_I4y2z_S;
    Double I_ESP_I3y3z_Px = I_ESP_Kx3y3z_S+ABX*I_ESP_I3y3z_S;
    Double I_ESP_I2y4z_Px = I_ESP_Kx2y4z_S+ABX*I_ESP_I2y4z_S;
    Double I_ESP_Iy5z_Px = I_ESP_Kxy5z_S+ABX*I_ESP_Iy5z_S;
    Double I_ESP_I5xy_Py = I_ESP_K5x2y_S+ABY*I_ESP_I5xy_S;
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
    Double I_ESP_I5xz_Pz = I_ESP_K5x2z_S+ABZ*I_ESP_I5xz_S;
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
     * totally 48 integrals are omitted 
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
    Double I_ESP_H4xz_Dxy = I_ESP_I4xyz_Px+ABY*I_ESP_H4xz_Px;
    Double I_ESP_H3xyz_Dxy = I_ESP_I3x2yz_Px+ABY*I_ESP_H3xyz_Px;
    Double I_ESP_H3x2z_Dxy = I_ESP_I3xy2z_Px+ABY*I_ESP_H3x2z_Px;
    Double I_ESP_H2x2yz_Dxy = I_ESP_I2x3yz_Px+ABY*I_ESP_H2x2yz_Px;
    Double I_ESP_H2xy2z_Dxy = I_ESP_I2x2y2z_Px+ABY*I_ESP_H2xy2z_Px;
    Double I_ESP_H2x3z_Dxy = I_ESP_I2xy3z_Px+ABY*I_ESP_H2x3z_Px;
    Double I_ESP_Hx3yz_Dxy = I_ESP_Ix4yz_Px+ABY*I_ESP_Hx3yz_Px;
    Double I_ESP_Hx2y2z_Dxy = I_ESP_Ix3y2z_Px+ABY*I_ESP_Hx2y2z_Px;
    Double I_ESP_Hxy3z_Dxy = I_ESP_Ix2y3z_Px+ABY*I_ESP_Hxy3z_Px;
    Double I_ESP_Hx4z_Dxy = I_ESP_Ixy4z_Px+ABY*I_ESP_Hx4z_Px;
    Double I_ESP_H4yz_Dxy = I_ESP_I5yz_Px+ABY*I_ESP_H4yz_Px;
    Double I_ESP_H3y2z_Dxy = I_ESP_I4y2z_Px+ABY*I_ESP_H3y2z_Px;
    Double I_ESP_H2y3z_Dxy = I_ESP_I3y3z_Px+ABY*I_ESP_H2y3z_Px;
    Double I_ESP_Hy4z_Dxy = I_ESP_I2y4z_Px+ABY*I_ESP_Hy4z_Px;
    Double I_ESP_H5z_Dxy = I_ESP_Iy5z_Px+ABY*I_ESP_H5z_Px;
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
     * totally 0 integrals are omitted 
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
    Double I_ESP_G4x_Fxyz = I_ESP_H4xz_Dxy+ABZ*I_ESP_G4x_Dxy;
    Double I_ESP_G3xy_Fxyz = I_ESP_H3xyz_Dxy+ABZ*I_ESP_G3xy_Dxy;
    Double I_ESP_G3xz_Fxyz = I_ESP_H3x2z_Dxy+ABZ*I_ESP_G3xz_Dxy;
    Double I_ESP_G2x2y_Fxyz = I_ESP_H2x2yz_Dxy+ABZ*I_ESP_G2x2y_Dxy;
    Double I_ESP_G2xyz_Fxyz = I_ESP_H2xy2z_Dxy+ABZ*I_ESP_G2xyz_Dxy;
    Double I_ESP_G2x2z_Fxyz = I_ESP_H2x3z_Dxy+ABZ*I_ESP_G2x2z_Dxy;
    Double I_ESP_Gx3y_Fxyz = I_ESP_Hx3yz_Dxy+ABZ*I_ESP_Gx3y_Dxy;
    Double I_ESP_Gx2yz_Fxyz = I_ESP_Hx2y2z_Dxy+ABZ*I_ESP_Gx2yz_Dxy;
    Double I_ESP_Gxy2z_Fxyz = I_ESP_Hxy3z_Dxy+ABZ*I_ESP_Gxy2z_Dxy;
    Double I_ESP_Gx3z_Fxyz = I_ESP_Hx4z_Dxy+ABZ*I_ESP_Gx3z_Dxy;
    Double I_ESP_G4y_Fxyz = I_ESP_H4yz_Dxy+ABZ*I_ESP_G4y_Dxy;
    Double I_ESP_G3yz_Fxyz = I_ESP_H3y2z_Dxy+ABZ*I_ESP_G3yz_Dxy;
    Double I_ESP_G2y2z_Fxyz = I_ESP_H2y3z_Dxy+ABZ*I_ESP_G2y2z_Dxy;
    Double I_ESP_Gy3z_Fxyz = I_ESP_Hy4z_Dxy+ABZ*I_ESP_Gy3z_Dxy;
    Double I_ESP_G4z_Fxyz = I_ESP_H5z_Dxy+ABZ*I_ESP_G4z_Dxy;
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
    Double I_ESP_G4x_Fy2z = I_ESP_H4xy_D2z+ABY*I_ESP_G4x_D2z;
    Double I_ESP_G3xy_Fy2z = I_ESP_H3x2y_D2z+ABY*I_ESP_G3xy_D2z;
    Double I_ESP_G3xz_Fy2z = I_ESP_H3xyz_D2z+ABY*I_ESP_G3xz_D2z;
    Double I_ESP_G2x2y_Fy2z = I_ESP_H2x3y_D2z+ABY*I_ESP_G2x2y_D2z;
    Double I_ESP_G2xyz_Fy2z = I_ESP_H2x2yz_D2z+ABY*I_ESP_G2xyz_D2z;
    Double I_ESP_G2x2z_Fy2z = I_ESP_H2xy2z_D2z+ABY*I_ESP_G2x2z_D2z;
    Double I_ESP_Gx3y_Fy2z = I_ESP_Hx4y_D2z+ABY*I_ESP_Gx3y_D2z;
    Double I_ESP_Gx2yz_Fy2z = I_ESP_Hx3yz_D2z+ABY*I_ESP_Gx2yz_D2z;
    Double I_ESP_Gxy2z_Fy2z = I_ESP_Hx2y2z_D2z+ABY*I_ESP_Gxy2z_D2z;
    Double I_ESP_Gx3z_Fy2z = I_ESP_Hxy3z_D2z+ABY*I_ESP_Gx3z_D2z;
    Double I_ESP_G4y_Fy2z = I_ESP_H5y_D2z+ABY*I_ESP_G4y_D2z;
    Double I_ESP_G3yz_Fy2z = I_ESP_H4yz_D2z+ABY*I_ESP_G3yz_D2z;
    Double I_ESP_G2y2z_Fy2z = I_ESP_H3y2z_D2z+ABY*I_ESP_G2y2z_D2z;
    Double I_ESP_Gy3z_Fy2z = I_ESP_H2y3z_D2z+ABY*I_ESP_Gy3z_D2z;
    Double I_ESP_G4z_Fy2z = I_ESP_Hy4z_D2z+ABY*I_ESP_G4z_D2z;
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
     * totally 56 integrals are omitted 
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
    Double I_ESP_I6x_Dxy_a = I_ESP_K6xy_Px_a+ABY*I_ESP_I6x_Px_a;
    Double I_ESP_I5xy_Dxy_a = I_ESP_K5x2y_Px_a+ABY*I_ESP_I5xy_Px_a;
    Double I_ESP_I5xz_Dxy_a = I_ESP_K5xyz_Px_a+ABY*I_ESP_I5xz_Px_a;
    Double I_ESP_I4x2y_Dxy_a = I_ESP_K4x3y_Px_a+ABY*I_ESP_I4x2y_Px_a;
    Double I_ESP_I4xyz_Dxy_a = I_ESP_K4x2yz_Px_a+ABY*I_ESP_I4xyz_Px_a;
    Double I_ESP_I4x2z_Dxy_a = I_ESP_K4xy2z_Px_a+ABY*I_ESP_I4x2z_Px_a;
    Double I_ESP_I3x3y_Dxy_a = I_ESP_K3x4y_Px_a+ABY*I_ESP_I3x3y_Px_a;
    Double I_ESP_I3x2yz_Dxy_a = I_ESP_K3x3yz_Px_a+ABY*I_ESP_I3x2yz_Px_a;
    Double I_ESP_I3xy2z_Dxy_a = I_ESP_K3x2y2z_Px_a+ABY*I_ESP_I3xy2z_Px_a;
    Double I_ESP_I3x3z_Dxy_a = I_ESP_K3xy3z_Px_a+ABY*I_ESP_I3x3z_Px_a;
    Double I_ESP_I2x4y_Dxy_a = I_ESP_K2x5y_Px_a+ABY*I_ESP_I2x4y_Px_a;
    Double I_ESP_I2x3yz_Dxy_a = I_ESP_K2x4yz_Px_a+ABY*I_ESP_I2x3yz_Px_a;
    Double I_ESP_I2x2y2z_Dxy_a = I_ESP_K2x3y2z_Px_a+ABY*I_ESP_I2x2y2z_Px_a;
    Double I_ESP_I2xy3z_Dxy_a = I_ESP_K2x2y3z_Px_a+ABY*I_ESP_I2xy3z_Px_a;
    Double I_ESP_I2x4z_Dxy_a = I_ESP_K2xy4z_Px_a+ABY*I_ESP_I2x4z_Px_a;
    Double I_ESP_Ix5y_Dxy_a = I_ESP_Kx6y_Px_a+ABY*I_ESP_Ix5y_Px_a;
    Double I_ESP_Ix4yz_Dxy_a = I_ESP_Kx5yz_Px_a+ABY*I_ESP_Ix4yz_Px_a;
    Double I_ESP_Ix3y2z_Dxy_a = I_ESP_Kx4y2z_Px_a+ABY*I_ESP_Ix3y2z_Px_a;
    Double I_ESP_Ix2y3z_Dxy_a = I_ESP_Kx3y3z_Px_a+ABY*I_ESP_Ix2y3z_Px_a;
    Double I_ESP_Ixy4z_Dxy_a = I_ESP_Kx2y4z_Px_a+ABY*I_ESP_Ixy4z_Px_a;
    Double I_ESP_Ix5z_Dxy_a = I_ESP_Kxy5z_Px_a+ABY*I_ESP_Ix5z_Px_a;
    Double I_ESP_I6y_Dxy_a = I_ESP_K7y_Px_a+ABY*I_ESP_I6y_Px_a;
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
     * shell quartet name: SQ_ESP_L_P_a
     * expanding position: BRA2
     * code section is: HRR
     * totally 20 integrals are omitted 
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
    Double I_ESP_L7yz_Px_a = I_ESP_Mx7yz_S_a+ABX*I_ESP_L7yz_S_a;
    Double I_ESP_L6y2z_Px_a = I_ESP_Mx6y2z_S_a+ABX*I_ESP_L6y2z_S_a;
    Double I_ESP_L5y3z_Px_a = I_ESP_Mx5y3z_S_a+ABX*I_ESP_L5y3z_S_a;
    Double I_ESP_L4y4z_Px_a = I_ESP_Mx4y4z_S_a+ABX*I_ESP_L4y4z_S_a;
    Double I_ESP_L3y5z_Px_a = I_ESP_Mx3y5z_S_a+ABX*I_ESP_L3y5z_S_a;
    Double I_ESP_L2y6z_Px_a = I_ESP_Mx2y6z_S_a+ABX*I_ESP_L2y6z_S_a;
    Double I_ESP_Ly7z_Px_a = I_ESP_Mxy7z_S_a+ABX*I_ESP_Ly7z_S_a;
    Double I_ESP_L7xy_Py_a = I_ESP_M7x2y_S_a+ABY*I_ESP_L7xy_S_a;
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
    Double I_ESP_L7xz_Pz_a = I_ESP_M7x2z_S_a+ABZ*I_ESP_L7xz_S_a;
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
     * totally 80 integrals are omitted 
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
    Double I_ESP_K6xz_Dxy_a = I_ESP_L6xyz_Px_a+ABY*I_ESP_K6xz_Px_a;
    Double I_ESP_K5xyz_Dxy_a = I_ESP_L5x2yz_Px_a+ABY*I_ESP_K5xyz_Px_a;
    Double I_ESP_K5x2z_Dxy_a = I_ESP_L5xy2z_Px_a+ABY*I_ESP_K5x2z_Px_a;
    Double I_ESP_K4x2yz_Dxy_a = I_ESP_L4x3yz_Px_a+ABY*I_ESP_K4x2yz_Px_a;
    Double I_ESP_K4xy2z_Dxy_a = I_ESP_L4x2y2z_Px_a+ABY*I_ESP_K4xy2z_Px_a;
    Double I_ESP_K4x3z_Dxy_a = I_ESP_L4xy3z_Px_a+ABY*I_ESP_K4x3z_Px_a;
    Double I_ESP_K3x3yz_Dxy_a = I_ESP_L3x4yz_Px_a+ABY*I_ESP_K3x3yz_Px_a;
    Double I_ESP_K3x2y2z_Dxy_a = I_ESP_L3x3y2z_Px_a+ABY*I_ESP_K3x2y2z_Px_a;
    Double I_ESP_K3xy3z_Dxy_a = I_ESP_L3x2y3z_Px_a+ABY*I_ESP_K3xy3z_Px_a;
    Double I_ESP_K3x4z_Dxy_a = I_ESP_L3xy4z_Px_a+ABY*I_ESP_K3x4z_Px_a;
    Double I_ESP_K2x4yz_Dxy_a = I_ESP_L2x5yz_Px_a+ABY*I_ESP_K2x4yz_Px_a;
    Double I_ESP_K2x3y2z_Dxy_a = I_ESP_L2x4y2z_Px_a+ABY*I_ESP_K2x3y2z_Px_a;
    Double I_ESP_K2x2y3z_Dxy_a = I_ESP_L2x3y3z_Px_a+ABY*I_ESP_K2x2y3z_Px_a;
    Double I_ESP_K2xy4z_Dxy_a = I_ESP_L2x2y4z_Px_a+ABY*I_ESP_K2xy4z_Px_a;
    Double I_ESP_K2x5z_Dxy_a = I_ESP_L2xy5z_Px_a+ABY*I_ESP_K2x5z_Px_a;
    Double I_ESP_Kx5yz_Dxy_a = I_ESP_Lx6yz_Px_a+ABY*I_ESP_Kx5yz_Px_a;
    Double I_ESP_Kx4y2z_Dxy_a = I_ESP_Lx5y2z_Px_a+ABY*I_ESP_Kx4y2z_Px_a;
    Double I_ESP_Kx3y3z_Dxy_a = I_ESP_Lx4y3z_Px_a+ABY*I_ESP_Kx3y3z_Px_a;
    Double I_ESP_Kx2y4z_Dxy_a = I_ESP_Lx3y4z_Px_a+ABY*I_ESP_Kx2y4z_Px_a;
    Double I_ESP_Kxy5z_Dxy_a = I_ESP_Lx2y5z_Px_a+ABY*I_ESP_Kxy5z_Px_a;
    Double I_ESP_Kx6z_Dxy_a = I_ESP_Lxy6z_Px_a+ABY*I_ESP_Kx6z_Px_a;
    Double I_ESP_K6yz_Dxy_a = I_ESP_L7yz_Px_a+ABY*I_ESP_K6yz_Px_a;
    Double I_ESP_K5y2z_Dxy_a = I_ESP_L6y2z_Px_a+ABY*I_ESP_K5y2z_Px_a;
    Double I_ESP_K4y3z_Dxy_a = I_ESP_L5y3z_Px_a+ABY*I_ESP_K4y3z_Px_a;
    Double I_ESP_K3y4z_Dxy_a = I_ESP_L4y4z_Px_a+ABY*I_ESP_K3y4z_Px_a;
    Double I_ESP_K2y5z_Dxy_a = I_ESP_L3y5z_Px_a+ABY*I_ESP_K2y5z_Px_a;
    Double I_ESP_Ky6z_Dxy_a = I_ESP_L2y6z_Px_a+ABY*I_ESP_Ky6z_Px_a;
    Double I_ESP_K7z_Dxy_a = I_ESP_Ly7z_Px_a+ABY*I_ESP_K7z_Px_a;
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
     * totally 0 integrals are omitted 
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
    Double I_ESP_I6x_Fxyz_a = I_ESP_K6xz_Dxy_a+ABZ*I_ESP_I6x_Dxy_a;
    Double I_ESP_I5xy_Fxyz_a = I_ESP_K5xyz_Dxy_a+ABZ*I_ESP_I5xy_Dxy_a;
    Double I_ESP_I5xz_Fxyz_a = I_ESP_K5x2z_Dxy_a+ABZ*I_ESP_I5xz_Dxy_a;
    Double I_ESP_I4x2y_Fxyz_a = I_ESP_K4x2yz_Dxy_a+ABZ*I_ESP_I4x2y_Dxy_a;
    Double I_ESP_I4xyz_Fxyz_a = I_ESP_K4xy2z_Dxy_a+ABZ*I_ESP_I4xyz_Dxy_a;
    Double I_ESP_I4x2z_Fxyz_a = I_ESP_K4x3z_Dxy_a+ABZ*I_ESP_I4x2z_Dxy_a;
    Double I_ESP_I3x3y_Fxyz_a = I_ESP_K3x3yz_Dxy_a+ABZ*I_ESP_I3x3y_Dxy_a;
    Double I_ESP_I3x2yz_Fxyz_a = I_ESP_K3x2y2z_Dxy_a+ABZ*I_ESP_I3x2yz_Dxy_a;
    Double I_ESP_I3xy2z_Fxyz_a = I_ESP_K3xy3z_Dxy_a+ABZ*I_ESP_I3xy2z_Dxy_a;
    Double I_ESP_I3x3z_Fxyz_a = I_ESP_K3x4z_Dxy_a+ABZ*I_ESP_I3x3z_Dxy_a;
    Double I_ESP_I2x4y_Fxyz_a = I_ESP_K2x4yz_Dxy_a+ABZ*I_ESP_I2x4y_Dxy_a;
    Double I_ESP_I2x3yz_Fxyz_a = I_ESP_K2x3y2z_Dxy_a+ABZ*I_ESP_I2x3yz_Dxy_a;
    Double I_ESP_I2x2y2z_Fxyz_a = I_ESP_K2x2y3z_Dxy_a+ABZ*I_ESP_I2x2y2z_Dxy_a;
    Double I_ESP_I2xy3z_Fxyz_a = I_ESP_K2xy4z_Dxy_a+ABZ*I_ESP_I2xy3z_Dxy_a;
    Double I_ESP_I2x4z_Fxyz_a = I_ESP_K2x5z_Dxy_a+ABZ*I_ESP_I2x4z_Dxy_a;
    Double I_ESP_Ix5y_Fxyz_a = I_ESP_Kx5yz_Dxy_a+ABZ*I_ESP_Ix5y_Dxy_a;
    Double I_ESP_Ix4yz_Fxyz_a = I_ESP_Kx4y2z_Dxy_a+ABZ*I_ESP_Ix4yz_Dxy_a;
    Double I_ESP_Ix3y2z_Fxyz_a = I_ESP_Kx3y3z_Dxy_a+ABZ*I_ESP_Ix3y2z_Dxy_a;
    Double I_ESP_Ix2y3z_Fxyz_a = I_ESP_Kx2y4z_Dxy_a+ABZ*I_ESP_Ix2y3z_Dxy_a;
    Double I_ESP_Ixy4z_Fxyz_a = I_ESP_Kxy5z_Dxy_a+ABZ*I_ESP_Ixy4z_Dxy_a;
    Double I_ESP_Ix5z_Fxyz_a = I_ESP_Kx6z_Dxy_a+ABZ*I_ESP_Ix5z_Dxy_a;
    Double I_ESP_I6y_Fxyz_a = I_ESP_K6yz_Dxy_a+ABZ*I_ESP_I6y_Dxy_a;
    Double I_ESP_I5yz_Fxyz_a = I_ESP_K5y2z_Dxy_a+ABZ*I_ESP_I5yz_Dxy_a;
    Double I_ESP_I4y2z_Fxyz_a = I_ESP_K4y3z_Dxy_a+ABZ*I_ESP_I4y2z_Dxy_a;
    Double I_ESP_I3y3z_Fxyz_a = I_ESP_K3y4z_Dxy_a+ABZ*I_ESP_I3y3z_Dxy_a;
    Double I_ESP_I2y4z_Fxyz_a = I_ESP_K2y5z_Dxy_a+ABZ*I_ESP_I2y4z_Dxy_a;
    Double I_ESP_Iy5z_Fxyz_a = I_ESP_Ky6z_Dxy_a+ABZ*I_ESP_Iy5z_Dxy_a;
    Double I_ESP_I6z_Fxyz_a = I_ESP_K7z_Dxy_a+ABZ*I_ESP_I6z_Dxy_a;
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
    Double I_ESP_I6x_Fy2z_a = I_ESP_K6xy_D2z_a+ABY*I_ESP_I6x_D2z_a;
    Double I_ESP_I5xy_Fy2z_a = I_ESP_K5x2y_D2z_a+ABY*I_ESP_I5xy_D2z_a;
    Double I_ESP_I5xz_Fy2z_a = I_ESP_K5xyz_D2z_a+ABY*I_ESP_I5xz_D2z_a;
    Double I_ESP_I4x2y_Fy2z_a = I_ESP_K4x3y_D2z_a+ABY*I_ESP_I4x2y_D2z_a;
    Double I_ESP_I4xyz_Fy2z_a = I_ESP_K4x2yz_D2z_a+ABY*I_ESP_I4xyz_D2z_a;
    Double I_ESP_I4x2z_Fy2z_a = I_ESP_K4xy2z_D2z_a+ABY*I_ESP_I4x2z_D2z_a;
    Double I_ESP_I3x3y_Fy2z_a = I_ESP_K3x4y_D2z_a+ABY*I_ESP_I3x3y_D2z_a;
    Double I_ESP_I3x2yz_Fy2z_a = I_ESP_K3x3yz_D2z_a+ABY*I_ESP_I3x2yz_D2z_a;
    Double I_ESP_I3xy2z_Fy2z_a = I_ESP_K3x2y2z_D2z_a+ABY*I_ESP_I3xy2z_D2z_a;
    Double I_ESP_I3x3z_Fy2z_a = I_ESP_K3xy3z_D2z_a+ABY*I_ESP_I3x3z_D2z_a;
    Double I_ESP_I2x4y_Fy2z_a = I_ESP_K2x5y_D2z_a+ABY*I_ESP_I2x4y_D2z_a;
    Double I_ESP_I2x3yz_Fy2z_a = I_ESP_K2x4yz_D2z_a+ABY*I_ESP_I2x3yz_D2z_a;
    Double I_ESP_I2x2y2z_Fy2z_a = I_ESP_K2x3y2z_D2z_a+ABY*I_ESP_I2x2y2z_D2z_a;
    Double I_ESP_I2xy3z_Fy2z_a = I_ESP_K2x2y3z_D2z_a+ABY*I_ESP_I2xy3z_D2z_a;
    Double I_ESP_I2x4z_Fy2z_a = I_ESP_K2xy4z_D2z_a+ABY*I_ESP_I2x4z_D2z_a;
    Double I_ESP_Ix5y_Fy2z_a = I_ESP_Kx6y_D2z_a+ABY*I_ESP_Ix5y_D2z_a;
    Double I_ESP_Ix4yz_Fy2z_a = I_ESP_Kx5yz_D2z_a+ABY*I_ESP_Ix4yz_D2z_a;
    Double I_ESP_Ix3y2z_Fy2z_a = I_ESP_Kx4y2z_D2z_a+ABY*I_ESP_Ix3y2z_D2z_a;
    Double I_ESP_Ix2y3z_Fy2z_a = I_ESP_Kx3y3z_D2z_a+ABY*I_ESP_Ix2y3z_D2z_a;
    Double I_ESP_Ixy4z_Fy2z_a = I_ESP_Kx2y4z_D2z_a+ABY*I_ESP_Ixy4z_D2z_a;
    Double I_ESP_Ix5z_Fy2z_a = I_ESP_Kxy5z_D2z_a+ABY*I_ESP_Ix5z_D2z_a;
    Double I_ESP_I6y_Fy2z_a = I_ESP_K7y_D2z_a+ABY*I_ESP_I6y_D2z_a;
    Double I_ESP_I5yz_Fy2z_a = I_ESP_K6yz_D2z_a+ABY*I_ESP_I5yz_D2z_a;
    Double I_ESP_I4y2z_Fy2z_a = I_ESP_K5y2z_D2z_a+ABY*I_ESP_I4y2z_D2z_a;
    Double I_ESP_I3y3z_Fy2z_a = I_ESP_K4y3z_D2z_a+ABY*I_ESP_I3y3z_D2z_a;
    Double I_ESP_I2y4z_Fy2z_a = I_ESP_K3y4z_D2z_a+ABY*I_ESP_I2y4z_D2z_a;
    Double I_ESP_Iy5z_Fy2z_a = I_ESP_K2y5z_D2z_a+ABY*I_ESP_Iy5z_D2z_a;
    Double I_ESP_I6z_Fy2z_a = I_ESP_Ky6z_D2z_a+ABY*I_ESP_I6z_D2z_a;
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
     * shell quartet name: SQ_ESP_H_F_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_a
     * RHS shell quartet name: SQ_ESP_G_F
     ************************************************************/
    abcd[iGrid*630+0] = 2.0E0*I_ESP_I6x_F3x_a-5*I_ESP_G4x_F3x;
    abcd[iGrid*630+1] = 2.0E0*I_ESP_I5xy_F3x_a-4*I_ESP_G3xy_F3x;
    abcd[iGrid*630+2] = 2.0E0*I_ESP_I5xz_F3x_a-4*I_ESP_G3xz_F3x;
    abcd[iGrid*630+3] = 2.0E0*I_ESP_I4x2y_F3x_a-3*I_ESP_G2x2y_F3x;
    abcd[iGrid*630+4] = 2.0E0*I_ESP_I4xyz_F3x_a-3*I_ESP_G2xyz_F3x;
    abcd[iGrid*630+5] = 2.0E0*I_ESP_I4x2z_F3x_a-3*I_ESP_G2x2z_F3x;
    abcd[iGrid*630+6] = 2.0E0*I_ESP_I3x3y_F3x_a-2*I_ESP_Gx3y_F3x;
    abcd[iGrid*630+7] = 2.0E0*I_ESP_I3x2yz_F3x_a-2*I_ESP_Gx2yz_F3x;
    abcd[iGrid*630+8] = 2.0E0*I_ESP_I3xy2z_F3x_a-2*I_ESP_Gxy2z_F3x;
    abcd[iGrid*630+9] = 2.0E0*I_ESP_I3x3z_F3x_a-2*I_ESP_Gx3z_F3x;
    abcd[iGrid*630+10] = 2.0E0*I_ESP_I2x4y_F3x_a-1*I_ESP_G4y_F3x;
    abcd[iGrid*630+11] = 2.0E0*I_ESP_I2x3yz_F3x_a-1*I_ESP_G3yz_F3x;
    abcd[iGrid*630+12] = 2.0E0*I_ESP_I2x2y2z_F3x_a-1*I_ESP_G2y2z_F3x;
    abcd[iGrid*630+13] = 2.0E0*I_ESP_I2xy3z_F3x_a-1*I_ESP_Gy3z_F3x;
    abcd[iGrid*630+14] = 2.0E0*I_ESP_I2x4z_F3x_a-1*I_ESP_G4z_F3x;
    abcd[iGrid*630+15] = 2.0E0*I_ESP_Ix5y_F3x_a;
    abcd[iGrid*630+16] = 2.0E0*I_ESP_Ix4yz_F3x_a;
    abcd[iGrid*630+17] = 2.0E0*I_ESP_Ix3y2z_F3x_a;
    abcd[iGrid*630+18] = 2.0E0*I_ESP_Ix2y3z_F3x_a;
    abcd[iGrid*630+19] = 2.0E0*I_ESP_Ixy4z_F3x_a;
    abcd[iGrid*630+20] = 2.0E0*I_ESP_Ix5z_F3x_a;
    abcd[iGrid*630+21] = 2.0E0*I_ESP_I6x_F2xy_a-5*I_ESP_G4x_F2xy;
    abcd[iGrid*630+22] = 2.0E0*I_ESP_I5xy_F2xy_a-4*I_ESP_G3xy_F2xy;
    abcd[iGrid*630+23] = 2.0E0*I_ESP_I5xz_F2xy_a-4*I_ESP_G3xz_F2xy;
    abcd[iGrid*630+24] = 2.0E0*I_ESP_I4x2y_F2xy_a-3*I_ESP_G2x2y_F2xy;
    abcd[iGrid*630+25] = 2.0E0*I_ESP_I4xyz_F2xy_a-3*I_ESP_G2xyz_F2xy;
    abcd[iGrid*630+26] = 2.0E0*I_ESP_I4x2z_F2xy_a-3*I_ESP_G2x2z_F2xy;
    abcd[iGrid*630+27] = 2.0E0*I_ESP_I3x3y_F2xy_a-2*I_ESP_Gx3y_F2xy;
    abcd[iGrid*630+28] = 2.0E0*I_ESP_I3x2yz_F2xy_a-2*I_ESP_Gx2yz_F2xy;
    abcd[iGrid*630+29] = 2.0E0*I_ESP_I3xy2z_F2xy_a-2*I_ESP_Gxy2z_F2xy;
    abcd[iGrid*630+30] = 2.0E0*I_ESP_I3x3z_F2xy_a-2*I_ESP_Gx3z_F2xy;
    abcd[iGrid*630+31] = 2.0E0*I_ESP_I2x4y_F2xy_a-1*I_ESP_G4y_F2xy;
    abcd[iGrid*630+32] = 2.0E0*I_ESP_I2x3yz_F2xy_a-1*I_ESP_G3yz_F2xy;
    abcd[iGrid*630+33] = 2.0E0*I_ESP_I2x2y2z_F2xy_a-1*I_ESP_G2y2z_F2xy;
    abcd[iGrid*630+34] = 2.0E0*I_ESP_I2xy3z_F2xy_a-1*I_ESP_Gy3z_F2xy;
    abcd[iGrid*630+35] = 2.0E0*I_ESP_I2x4z_F2xy_a-1*I_ESP_G4z_F2xy;
    abcd[iGrid*630+36] = 2.0E0*I_ESP_Ix5y_F2xy_a;
    abcd[iGrid*630+37] = 2.0E0*I_ESP_Ix4yz_F2xy_a;
    abcd[iGrid*630+38] = 2.0E0*I_ESP_Ix3y2z_F2xy_a;
    abcd[iGrid*630+39] = 2.0E0*I_ESP_Ix2y3z_F2xy_a;
    abcd[iGrid*630+40] = 2.0E0*I_ESP_Ixy4z_F2xy_a;
    abcd[iGrid*630+41] = 2.0E0*I_ESP_Ix5z_F2xy_a;
    abcd[iGrid*630+42] = 2.0E0*I_ESP_I6x_F2xz_a-5*I_ESP_G4x_F2xz;
    abcd[iGrid*630+43] = 2.0E0*I_ESP_I5xy_F2xz_a-4*I_ESP_G3xy_F2xz;
    abcd[iGrid*630+44] = 2.0E0*I_ESP_I5xz_F2xz_a-4*I_ESP_G3xz_F2xz;
    abcd[iGrid*630+45] = 2.0E0*I_ESP_I4x2y_F2xz_a-3*I_ESP_G2x2y_F2xz;
    abcd[iGrid*630+46] = 2.0E0*I_ESP_I4xyz_F2xz_a-3*I_ESP_G2xyz_F2xz;
    abcd[iGrid*630+47] = 2.0E0*I_ESP_I4x2z_F2xz_a-3*I_ESP_G2x2z_F2xz;
    abcd[iGrid*630+48] = 2.0E0*I_ESP_I3x3y_F2xz_a-2*I_ESP_Gx3y_F2xz;
    abcd[iGrid*630+49] = 2.0E0*I_ESP_I3x2yz_F2xz_a-2*I_ESP_Gx2yz_F2xz;
    abcd[iGrid*630+50] = 2.0E0*I_ESP_I3xy2z_F2xz_a-2*I_ESP_Gxy2z_F2xz;
    abcd[iGrid*630+51] = 2.0E0*I_ESP_I3x3z_F2xz_a-2*I_ESP_Gx3z_F2xz;
    abcd[iGrid*630+52] = 2.0E0*I_ESP_I2x4y_F2xz_a-1*I_ESP_G4y_F2xz;
    abcd[iGrid*630+53] = 2.0E0*I_ESP_I2x3yz_F2xz_a-1*I_ESP_G3yz_F2xz;
    abcd[iGrid*630+54] = 2.0E0*I_ESP_I2x2y2z_F2xz_a-1*I_ESP_G2y2z_F2xz;
    abcd[iGrid*630+55] = 2.0E0*I_ESP_I2xy3z_F2xz_a-1*I_ESP_Gy3z_F2xz;
    abcd[iGrid*630+56] = 2.0E0*I_ESP_I2x4z_F2xz_a-1*I_ESP_G4z_F2xz;
    abcd[iGrid*630+57] = 2.0E0*I_ESP_Ix5y_F2xz_a;
    abcd[iGrid*630+58] = 2.0E0*I_ESP_Ix4yz_F2xz_a;
    abcd[iGrid*630+59] = 2.0E0*I_ESP_Ix3y2z_F2xz_a;
    abcd[iGrid*630+60] = 2.0E0*I_ESP_Ix2y3z_F2xz_a;
    abcd[iGrid*630+61] = 2.0E0*I_ESP_Ixy4z_F2xz_a;
    abcd[iGrid*630+62] = 2.0E0*I_ESP_Ix5z_F2xz_a;
    abcd[iGrid*630+63] = 2.0E0*I_ESP_I6x_Fx2y_a-5*I_ESP_G4x_Fx2y;
    abcd[iGrid*630+64] = 2.0E0*I_ESP_I5xy_Fx2y_a-4*I_ESP_G3xy_Fx2y;
    abcd[iGrid*630+65] = 2.0E0*I_ESP_I5xz_Fx2y_a-4*I_ESP_G3xz_Fx2y;
    abcd[iGrid*630+66] = 2.0E0*I_ESP_I4x2y_Fx2y_a-3*I_ESP_G2x2y_Fx2y;
    abcd[iGrid*630+67] = 2.0E0*I_ESP_I4xyz_Fx2y_a-3*I_ESP_G2xyz_Fx2y;
    abcd[iGrid*630+68] = 2.0E0*I_ESP_I4x2z_Fx2y_a-3*I_ESP_G2x2z_Fx2y;
    abcd[iGrid*630+69] = 2.0E0*I_ESP_I3x3y_Fx2y_a-2*I_ESP_Gx3y_Fx2y;
    abcd[iGrid*630+70] = 2.0E0*I_ESP_I3x2yz_Fx2y_a-2*I_ESP_Gx2yz_Fx2y;
    abcd[iGrid*630+71] = 2.0E0*I_ESP_I3xy2z_Fx2y_a-2*I_ESP_Gxy2z_Fx2y;
    abcd[iGrid*630+72] = 2.0E0*I_ESP_I3x3z_Fx2y_a-2*I_ESP_Gx3z_Fx2y;
    abcd[iGrid*630+73] = 2.0E0*I_ESP_I2x4y_Fx2y_a-1*I_ESP_G4y_Fx2y;
    abcd[iGrid*630+74] = 2.0E0*I_ESP_I2x3yz_Fx2y_a-1*I_ESP_G3yz_Fx2y;
    abcd[iGrid*630+75] = 2.0E0*I_ESP_I2x2y2z_Fx2y_a-1*I_ESP_G2y2z_Fx2y;
    abcd[iGrid*630+76] = 2.0E0*I_ESP_I2xy3z_Fx2y_a-1*I_ESP_Gy3z_Fx2y;
    abcd[iGrid*630+77] = 2.0E0*I_ESP_I2x4z_Fx2y_a-1*I_ESP_G4z_Fx2y;
    abcd[iGrid*630+78] = 2.0E0*I_ESP_Ix5y_Fx2y_a;
    abcd[iGrid*630+79] = 2.0E0*I_ESP_Ix4yz_Fx2y_a;
    abcd[iGrid*630+80] = 2.0E0*I_ESP_Ix3y2z_Fx2y_a;
    abcd[iGrid*630+81] = 2.0E0*I_ESP_Ix2y3z_Fx2y_a;
    abcd[iGrid*630+82] = 2.0E0*I_ESP_Ixy4z_Fx2y_a;
    abcd[iGrid*630+83] = 2.0E0*I_ESP_Ix5z_Fx2y_a;
    abcd[iGrid*630+84] = 2.0E0*I_ESP_I6x_Fxyz_a-5*I_ESP_G4x_Fxyz;
    abcd[iGrid*630+85] = 2.0E0*I_ESP_I5xy_Fxyz_a-4*I_ESP_G3xy_Fxyz;
    abcd[iGrid*630+86] = 2.0E0*I_ESP_I5xz_Fxyz_a-4*I_ESP_G3xz_Fxyz;
    abcd[iGrid*630+87] = 2.0E0*I_ESP_I4x2y_Fxyz_a-3*I_ESP_G2x2y_Fxyz;
    abcd[iGrid*630+88] = 2.0E0*I_ESP_I4xyz_Fxyz_a-3*I_ESP_G2xyz_Fxyz;
    abcd[iGrid*630+89] = 2.0E0*I_ESP_I4x2z_Fxyz_a-3*I_ESP_G2x2z_Fxyz;
    abcd[iGrid*630+90] = 2.0E0*I_ESP_I3x3y_Fxyz_a-2*I_ESP_Gx3y_Fxyz;
    abcd[iGrid*630+91] = 2.0E0*I_ESP_I3x2yz_Fxyz_a-2*I_ESP_Gx2yz_Fxyz;
    abcd[iGrid*630+92] = 2.0E0*I_ESP_I3xy2z_Fxyz_a-2*I_ESP_Gxy2z_Fxyz;
    abcd[iGrid*630+93] = 2.0E0*I_ESP_I3x3z_Fxyz_a-2*I_ESP_Gx3z_Fxyz;
    abcd[iGrid*630+94] = 2.0E0*I_ESP_I2x4y_Fxyz_a-1*I_ESP_G4y_Fxyz;
    abcd[iGrid*630+95] = 2.0E0*I_ESP_I2x3yz_Fxyz_a-1*I_ESP_G3yz_Fxyz;
    abcd[iGrid*630+96] = 2.0E0*I_ESP_I2x2y2z_Fxyz_a-1*I_ESP_G2y2z_Fxyz;
    abcd[iGrid*630+97] = 2.0E0*I_ESP_I2xy3z_Fxyz_a-1*I_ESP_Gy3z_Fxyz;
    abcd[iGrid*630+98] = 2.0E0*I_ESP_I2x4z_Fxyz_a-1*I_ESP_G4z_Fxyz;
    abcd[iGrid*630+99] = 2.0E0*I_ESP_Ix5y_Fxyz_a;
    abcd[iGrid*630+100] = 2.0E0*I_ESP_Ix4yz_Fxyz_a;
    abcd[iGrid*630+101] = 2.0E0*I_ESP_Ix3y2z_Fxyz_a;
    abcd[iGrid*630+102] = 2.0E0*I_ESP_Ix2y3z_Fxyz_a;
    abcd[iGrid*630+103] = 2.0E0*I_ESP_Ixy4z_Fxyz_a;
    abcd[iGrid*630+104] = 2.0E0*I_ESP_Ix5z_Fxyz_a;
    abcd[iGrid*630+105] = 2.0E0*I_ESP_I6x_Fx2z_a-5*I_ESP_G4x_Fx2z;
    abcd[iGrid*630+106] = 2.0E0*I_ESP_I5xy_Fx2z_a-4*I_ESP_G3xy_Fx2z;
    abcd[iGrid*630+107] = 2.0E0*I_ESP_I5xz_Fx2z_a-4*I_ESP_G3xz_Fx2z;
    abcd[iGrid*630+108] = 2.0E0*I_ESP_I4x2y_Fx2z_a-3*I_ESP_G2x2y_Fx2z;
    abcd[iGrid*630+109] = 2.0E0*I_ESP_I4xyz_Fx2z_a-3*I_ESP_G2xyz_Fx2z;
    abcd[iGrid*630+110] = 2.0E0*I_ESP_I4x2z_Fx2z_a-3*I_ESP_G2x2z_Fx2z;
    abcd[iGrid*630+111] = 2.0E0*I_ESP_I3x3y_Fx2z_a-2*I_ESP_Gx3y_Fx2z;
    abcd[iGrid*630+112] = 2.0E0*I_ESP_I3x2yz_Fx2z_a-2*I_ESP_Gx2yz_Fx2z;
    abcd[iGrid*630+113] = 2.0E0*I_ESP_I3xy2z_Fx2z_a-2*I_ESP_Gxy2z_Fx2z;
    abcd[iGrid*630+114] = 2.0E0*I_ESP_I3x3z_Fx2z_a-2*I_ESP_Gx3z_Fx2z;
    abcd[iGrid*630+115] = 2.0E0*I_ESP_I2x4y_Fx2z_a-1*I_ESP_G4y_Fx2z;
    abcd[iGrid*630+116] = 2.0E0*I_ESP_I2x3yz_Fx2z_a-1*I_ESP_G3yz_Fx2z;
    abcd[iGrid*630+117] = 2.0E0*I_ESP_I2x2y2z_Fx2z_a-1*I_ESP_G2y2z_Fx2z;
    abcd[iGrid*630+118] = 2.0E0*I_ESP_I2xy3z_Fx2z_a-1*I_ESP_Gy3z_Fx2z;
    abcd[iGrid*630+119] = 2.0E0*I_ESP_I2x4z_Fx2z_a-1*I_ESP_G4z_Fx2z;
    abcd[iGrid*630+120] = 2.0E0*I_ESP_Ix5y_Fx2z_a;
    abcd[iGrid*630+121] = 2.0E0*I_ESP_Ix4yz_Fx2z_a;
    abcd[iGrid*630+122] = 2.0E0*I_ESP_Ix3y2z_Fx2z_a;
    abcd[iGrid*630+123] = 2.0E0*I_ESP_Ix2y3z_Fx2z_a;
    abcd[iGrid*630+124] = 2.0E0*I_ESP_Ixy4z_Fx2z_a;
    abcd[iGrid*630+125] = 2.0E0*I_ESP_Ix5z_Fx2z_a;
    abcd[iGrid*630+126] = 2.0E0*I_ESP_I6x_F3y_a-5*I_ESP_G4x_F3y;
    abcd[iGrid*630+127] = 2.0E0*I_ESP_I5xy_F3y_a-4*I_ESP_G3xy_F3y;
    abcd[iGrid*630+128] = 2.0E0*I_ESP_I5xz_F3y_a-4*I_ESP_G3xz_F3y;
    abcd[iGrid*630+129] = 2.0E0*I_ESP_I4x2y_F3y_a-3*I_ESP_G2x2y_F3y;
    abcd[iGrid*630+130] = 2.0E0*I_ESP_I4xyz_F3y_a-3*I_ESP_G2xyz_F3y;
    abcd[iGrid*630+131] = 2.0E0*I_ESP_I4x2z_F3y_a-3*I_ESP_G2x2z_F3y;
    abcd[iGrid*630+132] = 2.0E0*I_ESP_I3x3y_F3y_a-2*I_ESP_Gx3y_F3y;
    abcd[iGrid*630+133] = 2.0E0*I_ESP_I3x2yz_F3y_a-2*I_ESP_Gx2yz_F3y;
    abcd[iGrid*630+134] = 2.0E0*I_ESP_I3xy2z_F3y_a-2*I_ESP_Gxy2z_F3y;
    abcd[iGrid*630+135] = 2.0E0*I_ESP_I3x3z_F3y_a-2*I_ESP_Gx3z_F3y;
    abcd[iGrid*630+136] = 2.0E0*I_ESP_I2x4y_F3y_a-1*I_ESP_G4y_F3y;
    abcd[iGrid*630+137] = 2.0E0*I_ESP_I2x3yz_F3y_a-1*I_ESP_G3yz_F3y;
    abcd[iGrid*630+138] = 2.0E0*I_ESP_I2x2y2z_F3y_a-1*I_ESP_G2y2z_F3y;
    abcd[iGrid*630+139] = 2.0E0*I_ESP_I2xy3z_F3y_a-1*I_ESP_Gy3z_F3y;
    abcd[iGrid*630+140] = 2.0E0*I_ESP_I2x4z_F3y_a-1*I_ESP_G4z_F3y;
    abcd[iGrid*630+141] = 2.0E0*I_ESP_Ix5y_F3y_a;
    abcd[iGrid*630+142] = 2.0E0*I_ESP_Ix4yz_F3y_a;
    abcd[iGrid*630+143] = 2.0E0*I_ESP_Ix3y2z_F3y_a;
    abcd[iGrid*630+144] = 2.0E0*I_ESP_Ix2y3z_F3y_a;
    abcd[iGrid*630+145] = 2.0E0*I_ESP_Ixy4z_F3y_a;
    abcd[iGrid*630+146] = 2.0E0*I_ESP_Ix5z_F3y_a;
    abcd[iGrid*630+147] = 2.0E0*I_ESP_I6x_F2yz_a-5*I_ESP_G4x_F2yz;
    abcd[iGrid*630+148] = 2.0E0*I_ESP_I5xy_F2yz_a-4*I_ESP_G3xy_F2yz;
    abcd[iGrid*630+149] = 2.0E0*I_ESP_I5xz_F2yz_a-4*I_ESP_G3xz_F2yz;
    abcd[iGrid*630+150] = 2.0E0*I_ESP_I4x2y_F2yz_a-3*I_ESP_G2x2y_F2yz;
    abcd[iGrid*630+151] = 2.0E0*I_ESP_I4xyz_F2yz_a-3*I_ESP_G2xyz_F2yz;
    abcd[iGrid*630+152] = 2.0E0*I_ESP_I4x2z_F2yz_a-3*I_ESP_G2x2z_F2yz;
    abcd[iGrid*630+153] = 2.0E0*I_ESP_I3x3y_F2yz_a-2*I_ESP_Gx3y_F2yz;
    abcd[iGrid*630+154] = 2.0E0*I_ESP_I3x2yz_F2yz_a-2*I_ESP_Gx2yz_F2yz;
    abcd[iGrid*630+155] = 2.0E0*I_ESP_I3xy2z_F2yz_a-2*I_ESP_Gxy2z_F2yz;
    abcd[iGrid*630+156] = 2.0E0*I_ESP_I3x3z_F2yz_a-2*I_ESP_Gx3z_F2yz;
    abcd[iGrid*630+157] = 2.0E0*I_ESP_I2x4y_F2yz_a-1*I_ESP_G4y_F2yz;
    abcd[iGrid*630+158] = 2.0E0*I_ESP_I2x3yz_F2yz_a-1*I_ESP_G3yz_F2yz;
    abcd[iGrid*630+159] = 2.0E0*I_ESP_I2x2y2z_F2yz_a-1*I_ESP_G2y2z_F2yz;
    abcd[iGrid*630+160] = 2.0E0*I_ESP_I2xy3z_F2yz_a-1*I_ESP_Gy3z_F2yz;
    abcd[iGrid*630+161] = 2.0E0*I_ESP_I2x4z_F2yz_a-1*I_ESP_G4z_F2yz;
    abcd[iGrid*630+162] = 2.0E0*I_ESP_Ix5y_F2yz_a;
    abcd[iGrid*630+163] = 2.0E0*I_ESP_Ix4yz_F2yz_a;
    abcd[iGrid*630+164] = 2.0E0*I_ESP_Ix3y2z_F2yz_a;
    abcd[iGrid*630+165] = 2.0E0*I_ESP_Ix2y3z_F2yz_a;
    abcd[iGrid*630+166] = 2.0E0*I_ESP_Ixy4z_F2yz_a;
    abcd[iGrid*630+167] = 2.0E0*I_ESP_Ix5z_F2yz_a;
    abcd[iGrid*630+168] = 2.0E0*I_ESP_I6x_Fy2z_a-5*I_ESP_G4x_Fy2z;
    abcd[iGrid*630+169] = 2.0E0*I_ESP_I5xy_Fy2z_a-4*I_ESP_G3xy_Fy2z;
    abcd[iGrid*630+170] = 2.0E0*I_ESP_I5xz_Fy2z_a-4*I_ESP_G3xz_Fy2z;
    abcd[iGrid*630+171] = 2.0E0*I_ESP_I4x2y_Fy2z_a-3*I_ESP_G2x2y_Fy2z;
    abcd[iGrid*630+172] = 2.0E0*I_ESP_I4xyz_Fy2z_a-3*I_ESP_G2xyz_Fy2z;
    abcd[iGrid*630+173] = 2.0E0*I_ESP_I4x2z_Fy2z_a-3*I_ESP_G2x2z_Fy2z;
    abcd[iGrid*630+174] = 2.0E0*I_ESP_I3x3y_Fy2z_a-2*I_ESP_Gx3y_Fy2z;
    abcd[iGrid*630+175] = 2.0E0*I_ESP_I3x2yz_Fy2z_a-2*I_ESP_Gx2yz_Fy2z;
    abcd[iGrid*630+176] = 2.0E0*I_ESP_I3xy2z_Fy2z_a-2*I_ESP_Gxy2z_Fy2z;
    abcd[iGrid*630+177] = 2.0E0*I_ESP_I3x3z_Fy2z_a-2*I_ESP_Gx3z_Fy2z;
    abcd[iGrid*630+178] = 2.0E0*I_ESP_I2x4y_Fy2z_a-1*I_ESP_G4y_Fy2z;
    abcd[iGrid*630+179] = 2.0E0*I_ESP_I2x3yz_Fy2z_a-1*I_ESP_G3yz_Fy2z;
    abcd[iGrid*630+180] = 2.0E0*I_ESP_I2x2y2z_Fy2z_a-1*I_ESP_G2y2z_Fy2z;
    abcd[iGrid*630+181] = 2.0E0*I_ESP_I2xy3z_Fy2z_a-1*I_ESP_Gy3z_Fy2z;
    abcd[iGrid*630+182] = 2.0E0*I_ESP_I2x4z_Fy2z_a-1*I_ESP_G4z_Fy2z;
    abcd[iGrid*630+183] = 2.0E0*I_ESP_Ix5y_Fy2z_a;
    abcd[iGrid*630+184] = 2.0E0*I_ESP_Ix4yz_Fy2z_a;
    abcd[iGrid*630+185] = 2.0E0*I_ESP_Ix3y2z_Fy2z_a;
    abcd[iGrid*630+186] = 2.0E0*I_ESP_Ix2y3z_Fy2z_a;
    abcd[iGrid*630+187] = 2.0E0*I_ESP_Ixy4z_Fy2z_a;
    abcd[iGrid*630+188] = 2.0E0*I_ESP_Ix5z_Fy2z_a;
    abcd[iGrid*630+189] = 2.0E0*I_ESP_I6x_F3z_a-5*I_ESP_G4x_F3z;
    abcd[iGrid*630+190] = 2.0E0*I_ESP_I5xy_F3z_a-4*I_ESP_G3xy_F3z;
    abcd[iGrid*630+191] = 2.0E0*I_ESP_I5xz_F3z_a-4*I_ESP_G3xz_F3z;
    abcd[iGrid*630+192] = 2.0E0*I_ESP_I4x2y_F3z_a-3*I_ESP_G2x2y_F3z;
    abcd[iGrid*630+193] = 2.0E0*I_ESP_I4xyz_F3z_a-3*I_ESP_G2xyz_F3z;
    abcd[iGrid*630+194] = 2.0E0*I_ESP_I4x2z_F3z_a-3*I_ESP_G2x2z_F3z;
    abcd[iGrid*630+195] = 2.0E0*I_ESP_I3x3y_F3z_a-2*I_ESP_Gx3y_F3z;
    abcd[iGrid*630+196] = 2.0E0*I_ESP_I3x2yz_F3z_a-2*I_ESP_Gx2yz_F3z;
    abcd[iGrid*630+197] = 2.0E0*I_ESP_I3xy2z_F3z_a-2*I_ESP_Gxy2z_F3z;
    abcd[iGrid*630+198] = 2.0E0*I_ESP_I3x3z_F3z_a-2*I_ESP_Gx3z_F3z;
    abcd[iGrid*630+199] = 2.0E0*I_ESP_I2x4y_F3z_a-1*I_ESP_G4y_F3z;
    abcd[iGrid*630+200] = 2.0E0*I_ESP_I2x3yz_F3z_a-1*I_ESP_G3yz_F3z;
    abcd[iGrid*630+201] = 2.0E0*I_ESP_I2x2y2z_F3z_a-1*I_ESP_G2y2z_F3z;
    abcd[iGrid*630+202] = 2.0E0*I_ESP_I2xy3z_F3z_a-1*I_ESP_Gy3z_F3z;
    abcd[iGrid*630+203] = 2.0E0*I_ESP_I2x4z_F3z_a-1*I_ESP_G4z_F3z;
    abcd[iGrid*630+204] = 2.0E0*I_ESP_Ix5y_F3z_a;
    abcd[iGrid*630+205] = 2.0E0*I_ESP_Ix4yz_F3z_a;
    abcd[iGrid*630+206] = 2.0E0*I_ESP_Ix3y2z_F3z_a;
    abcd[iGrid*630+207] = 2.0E0*I_ESP_Ix2y3z_F3z_a;
    abcd[iGrid*630+208] = 2.0E0*I_ESP_Ixy4z_F3z_a;
    abcd[iGrid*630+209] = 2.0E0*I_ESP_Ix5z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_a
     * RHS shell quartet name: SQ_ESP_G_F
     ************************************************************/
    abcd[iGrid*630+210] = 2.0E0*I_ESP_I5xy_F3x_a;
    abcd[iGrid*630+211] = 2.0E0*I_ESP_I4x2y_F3x_a-1*I_ESP_G4x_F3x;
    abcd[iGrid*630+212] = 2.0E0*I_ESP_I4xyz_F3x_a;
    abcd[iGrid*630+213] = 2.0E0*I_ESP_I3x3y_F3x_a-2*I_ESP_G3xy_F3x;
    abcd[iGrid*630+214] = 2.0E0*I_ESP_I3x2yz_F3x_a-1*I_ESP_G3xz_F3x;
    abcd[iGrid*630+215] = 2.0E0*I_ESP_I3xy2z_F3x_a;
    abcd[iGrid*630+216] = 2.0E0*I_ESP_I2x4y_F3x_a-3*I_ESP_G2x2y_F3x;
    abcd[iGrid*630+217] = 2.0E0*I_ESP_I2x3yz_F3x_a-2*I_ESP_G2xyz_F3x;
    abcd[iGrid*630+218] = 2.0E0*I_ESP_I2x2y2z_F3x_a-1*I_ESP_G2x2z_F3x;
    abcd[iGrid*630+219] = 2.0E0*I_ESP_I2xy3z_F3x_a;
    abcd[iGrid*630+220] = 2.0E0*I_ESP_Ix5y_F3x_a-4*I_ESP_Gx3y_F3x;
    abcd[iGrid*630+221] = 2.0E0*I_ESP_Ix4yz_F3x_a-3*I_ESP_Gx2yz_F3x;
    abcd[iGrid*630+222] = 2.0E0*I_ESP_Ix3y2z_F3x_a-2*I_ESP_Gxy2z_F3x;
    abcd[iGrid*630+223] = 2.0E0*I_ESP_Ix2y3z_F3x_a-1*I_ESP_Gx3z_F3x;
    abcd[iGrid*630+224] = 2.0E0*I_ESP_Ixy4z_F3x_a;
    abcd[iGrid*630+225] = 2.0E0*I_ESP_I6y_F3x_a-5*I_ESP_G4y_F3x;
    abcd[iGrid*630+226] = 2.0E0*I_ESP_I5yz_F3x_a-4*I_ESP_G3yz_F3x;
    abcd[iGrid*630+227] = 2.0E0*I_ESP_I4y2z_F3x_a-3*I_ESP_G2y2z_F3x;
    abcd[iGrid*630+228] = 2.0E0*I_ESP_I3y3z_F3x_a-2*I_ESP_Gy3z_F3x;
    abcd[iGrid*630+229] = 2.0E0*I_ESP_I2y4z_F3x_a-1*I_ESP_G4z_F3x;
    abcd[iGrid*630+230] = 2.0E0*I_ESP_Iy5z_F3x_a;
    abcd[iGrid*630+231] = 2.0E0*I_ESP_I5xy_F2xy_a;
    abcd[iGrid*630+232] = 2.0E0*I_ESP_I4x2y_F2xy_a-1*I_ESP_G4x_F2xy;
    abcd[iGrid*630+233] = 2.0E0*I_ESP_I4xyz_F2xy_a;
    abcd[iGrid*630+234] = 2.0E0*I_ESP_I3x3y_F2xy_a-2*I_ESP_G3xy_F2xy;
    abcd[iGrid*630+235] = 2.0E0*I_ESP_I3x2yz_F2xy_a-1*I_ESP_G3xz_F2xy;
    abcd[iGrid*630+236] = 2.0E0*I_ESP_I3xy2z_F2xy_a;
    abcd[iGrid*630+237] = 2.0E0*I_ESP_I2x4y_F2xy_a-3*I_ESP_G2x2y_F2xy;
    abcd[iGrid*630+238] = 2.0E0*I_ESP_I2x3yz_F2xy_a-2*I_ESP_G2xyz_F2xy;
    abcd[iGrid*630+239] = 2.0E0*I_ESP_I2x2y2z_F2xy_a-1*I_ESP_G2x2z_F2xy;
    abcd[iGrid*630+240] = 2.0E0*I_ESP_I2xy3z_F2xy_a;
    abcd[iGrid*630+241] = 2.0E0*I_ESP_Ix5y_F2xy_a-4*I_ESP_Gx3y_F2xy;
    abcd[iGrid*630+242] = 2.0E0*I_ESP_Ix4yz_F2xy_a-3*I_ESP_Gx2yz_F2xy;
    abcd[iGrid*630+243] = 2.0E0*I_ESP_Ix3y2z_F2xy_a-2*I_ESP_Gxy2z_F2xy;
    abcd[iGrid*630+244] = 2.0E0*I_ESP_Ix2y3z_F2xy_a-1*I_ESP_Gx3z_F2xy;
    abcd[iGrid*630+245] = 2.0E0*I_ESP_Ixy4z_F2xy_a;
    abcd[iGrid*630+246] = 2.0E0*I_ESP_I6y_F2xy_a-5*I_ESP_G4y_F2xy;
    abcd[iGrid*630+247] = 2.0E0*I_ESP_I5yz_F2xy_a-4*I_ESP_G3yz_F2xy;
    abcd[iGrid*630+248] = 2.0E0*I_ESP_I4y2z_F2xy_a-3*I_ESP_G2y2z_F2xy;
    abcd[iGrid*630+249] = 2.0E0*I_ESP_I3y3z_F2xy_a-2*I_ESP_Gy3z_F2xy;
    abcd[iGrid*630+250] = 2.0E0*I_ESP_I2y4z_F2xy_a-1*I_ESP_G4z_F2xy;
    abcd[iGrid*630+251] = 2.0E0*I_ESP_Iy5z_F2xy_a;
    abcd[iGrid*630+252] = 2.0E0*I_ESP_I5xy_F2xz_a;
    abcd[iGrid*630+253] = 2.0E0*I_ESP_I4x2y_F2xz_a-1*I_ESP_G4x_F2xz;
    abcd[iGrid*630+254] = 2.0E0*I_ESP_I4xyz_F2xz_a;
    abcd[iGrid*630+255] = 2.0E0*I_ESP_I3x3y_F2xz_a-2*I_ESP_G3xy_F2xz;
    abcd[iGrid*630+256] = 2.0E0*I_ESP_I3x2yz_F2xz_a-1*I_ESP_G3xz_F2xz;
    abcd[iGrid*630+257] = 2.0E0*I_ESP_I3xy2z_F2xz_a;
    abcd[iGrid*630+258] = 2.0E0*I_ESP_I2x4y_F2xz_a-3*I_ESP_G2x2y_F2xz;
    abcd[iGrid*630+259] = 2.0E0*I_ESP_I2x3yz_F2xz_a-2*I_ESP_G2xyz_F2xz;
    abcd[iGrid*630+260] = 2.0E0*I_ESP_I2x2y2z_F2xz_a-1*I_ESP_G2x2z_F2xz;
    abcd[iGrid*630+261] = 2.0E0*I_ESP_I2xy3z_F2xz_a;
    abcd[iGrid*630+262] = 2.0E0*I_ESP_Ix5y_F2xz_a-4*I_ESP_Gx3y_F2xz;
    abcd[iGrid*630+263] = 2.0E0*I_ESP_Ix4yz_F2xz_a-3*I_ESP_Gx2yz_F2xz;
    abcd[iGrid*630+264] = 2.0E0*I_ESP_Ix3y2z_F2xz_a-2*I_ESP_Gxy2z_F2xz;
    abcd[iGrid*630+265] = 2.0E0*I_ESP_Ix2y3z_F2xz_a-1*I_ESP_Gx3z_F2xz;
    abcd[iGrid*630+266] = 2.0E0*I_ESP_Ixy4z_F2xz_a;
    abcd[iGrid*630+267] = 2.0E0*I_ESP_I6y_F2xz_a-5*I_ESP_G4y_F2xz;
    abcd[iGrid*630+268] = 2.0E0*I_ESP_I5yz_F2xz_a-4*I_ESP_G3yz_F2xz;
    abcd[iGrid*630+269] = 2.0E0*I_ESP_I4y2z_F2xz_a-3*I_ESP_G2y2z_F2xz;
    abcd[iGrid*630+270] = 2.0E0*I_ESP_I3y3z_F2xz_a-2*I_ESP_Gy3z_F2xz;
    abcd[iGrid*630+271] = 2.0E0*I_ESP_I2y4z_F2xz_a-1*I_ESP_G4z_F2xz;
    abcd[iGrid*630+272] = 2.0E0*I_ESP_Iy5z_F2xz_a;
    abcd[iGrid*630+273] = 2.0E0*I_ESP_I5xy_Fx2y_a;
    abcd[iGrid*630+274] = 2.0E0*I_ESP_I4x2y_Fx2y_a-1*I_ESP_G4x_Fx2y;
    abcd[iGrid*630+275] = 2.0E0*I_ESP_I4xyz_Fx2y_a;
    abcd[iGrid*630+276] = 2.0E0*I_ESP_I3x3y_Fx2y_a-2*I_ESP_G3xy_Fx2y;
    abcd[iGrid*630+277] = 2.0E0*I_ESP_I3x2yz_Fx2y_a-1*I_ESP_G3xz_Fx2y;
    abcd[iGrid*630+278] = 2.0E0*I_ESP_I3xy2z_Fx2y_a;
    abcd[iGrid*630+279] = 2.0E0*I_ESP_I2x4y_Fx2y_a-3*I_ESP_G2x2y_Fx2y;
    abcd[iGrid*630+280] = 2.0E0*I_ESP_I2x3yz_Fx2y_a-2*I_ESP_G2xyz_Fx2y;
    abcd[iGrid*630+281] = 2.0E0*I_ESP_I2x2y2z_Fx2y_a-1*I_ESP_G2x2z_Fx2y;
    abcd[iGrid*630+282] = 2.0E0*I_ESP_I2xy3z_Fx2y_a;
    abcd[iGrid*630+283] = 2.0E0*I_ESP_Ix5y_Fx2y_a-4*I_ESP_Gx3y_Fx2y;
    abcd[iGrid*630+284] = 2.0E0*I_ESP_Ix4yz_Fx2y_a-3*I_ESP_Gx2yz_Fx2y;
    abcd[iGrid*630+285] = 2.0E0*I_ESP_Ix3y2z_Fx2y_a-2*I_ESP_Gxy2z_Fx2y;
    abcd[iGrid*630+286] = 2.0E0*I_ESP_Ix2y3z_Fx2y_a-1*I_ESP_Gx3z_Fx2y;
    abcd[iGrid*630+287] = 2.0E0*I_ESP_Ixy4z_Fx2y_a;
    abcd[iGrid*630+288] = 2.0E0*I_ESP_I6y_Fx2y_a-5*I_ESP_G4y_Fx2y;
    abcd[iGrid*630+289] = 2.0E0*I_ESP_I5yz_Fx2y_a-4*I_ESP_G3yz_Fx2y;
    abcd[iGrid*630+290] = 2.0E0*I_ESP_I4y2z_Fx2y_a-3*I_ESP_G2y2z_Fx2y;
    abcd[iGrid*630+291] = 2.0E0*I_ESP_I3y3z_Fx2y_a-2*I_ESP_Gy3z_Fx2y;
    abcd[iGrid*630+292] = 2.0E0*I_ESP_I2y4z_Fx2y_a-1*I_ESP_G4z_Fx2y;
    abcd[iGrid*630+293] = 2.0E0*I_ESP_Iy5z_Fx2y_a;
    abcd[iGrid*630+294] = 2.0E0*I_ESP_I5xy_Fxyz_a;
    abcd[iGrid*630+295] = 2.0E0*I_ESP_I4x2y_Fxyz_a-1*I_ESP_G4x_Fxyz;
    abcd[iGrid*630+296] = 2.0E0*I_ESP_I4xyz_Fxyz_a;
    abcd[iGrid*630+297] = 2.0E0*I_ESP_I3x3y_Fxyz_a-2*I_ESP_G3xy_Fxyz;
    abcd[iGrid*630+298] = 2.0E0*I_ESP_I3x2yz_Fxyz_a-1*I_ESP_G3xz_Fxyz;
    abcd[iGrid*630+299] = 2.0E0*I_ESP_I3xy2z_Fxyz_a;
    abcd[iGrid*630+300] = 2.0E0*I_ESP_I2x4y_Fxyz_a-3*I_ESP_G2x2y_Fxyz;
    abcd[iGrid*630+301] = 2.0E0*I_ESP_I2x3yz_Fxyz_a-2*I_ESP_G2xyz_Fxyz;
    abcd[iGrid*630+302] = 2.0E0*I_ESP_I2x2y2z_Fxyz_a-1*I_ESP_G2x2z_Fxyz;
    abcd[iGrid*630+303] = 2.0E0*I_ESP_I2xy3z_Fxyz_a;
    abcd[iGrid*630+304] = 2.0E0*I_ESP_Ix5y_Fxyz_a-4*I_ESP_Gx3y_Fxyz;
    abcd[iGrid*630+305] = 2.0E0*I_ESP_Ix4yz_Fxyz_a-3*I_ESP_Gx2yz_Fxyz;
    abcd[iGrid*630+306] = 2.0E0*I_ESP_Ix3y2z_Fxyz_a-2*I_ESP_Gxy2z_Fxyz;
    abcd[iGrid*630+307] = 2.0E0*I_ESP_Ix2y3z_Fxyz_a-1*I_ESP_Gx3z_Fxyz;
    abcd[iGrid*630+308] = 2.0E0*I_ESP_Ixy4z_Fxyz_a;
    abcd[iGrid*630+309] = 2.0E0*I_ESP_I6y_Fxyz_a-5*I_ESP_G4y_Fxyz;
    abcd[iGrid*630+310] = 2.0E0*I_ESP_I5yz_Fxyz_a-4*I_ESP_G3yz_Fxyz;
    abcd[iGrid*630+311] = 2.0E0*I_ESP_I4y2z_Fxyz_a-3*I_ESP_G2y2z_Fxyz;
    abcd[iGrid*630+312] = 2.0E0*I_ESP_I3y3z_Fxyz_a-2*I_ESP_Gy3z_Fxyz;
    abcd[iGrid*630+313] = 2.0E0*I_ESP_I2y4z_Fxyz_a-1*I_ESP_G4z_Fxyz;
    abcd[iGrid*630+314] = 2.0E0*I_ESP_Iy5z_Fxyz_a;
    abcd[iGrid*630+315] = 2.0E0*I_ESP_I5xy_Fx2z_a;
    abcd[iGrid*630+316] = 2.0E0*I_ESP_I4x2y_Fx2z_a-1*I_ESP_G4x_Fx2z;
    abcd[iGrid*630+317] = 2.0E0*I_ESP_I4xyz_Fx2z_a;
    abcd[iGrid*630+318] = 2.0E0*I_ESP_I3x3y_Fx2z_a-2*I_ESP_G3xy_Fx2z;
    abcd[iGrid*630+319] = 2.0E0*I_ESP_I3x2yz_Fx2z_a-1*I_ESP_G3xz_Fx2z;
    abcd[iGrid*630+320] = 2.0E0*I_ESP_I3xy2z_Fx2z_a;
    abcd[iGrid*630+321] = 2.0E0*I_ESP_I2x4y_Fx2z_a-3*I_ESP_G2x2y_Fx2z;
    abcd[iGrid*630+322] = 2.0E0*I_ESP_I2x3yz_Fx2z_a-2*I_ESP_G2xyz_Fx2z;
    abcd[iGrid*630+323] = 2.0E0*I_ESP_I2x2y2z_Fx2z_a-1*I_ESP_G2x2z_Fx2z;
    abcd[iGrid*630+324] = 2.0E0*I_ESP_I2xy3z_Fx2z_a;
    abcd[iGrid*630+325] = 2.0E0*I_ESP_Ix5y_Fx2z_a-4*I_ESP_Gx3y_Fx2z;
    abcd[iGrid*630+326] = 2.0E0*I_ESP_Ix4yz_Fx2z_a-3*I_ESP_Gx2yz_Fx2z;
    abcd[iGrid*630+327] = 2.0E0*I_ESP_Ix3y2z_Fx2z_a-2*I_ESP_Gxy2z_Fx2z;
    abcd[iGrid*630+328] = 2.0E0*I_ESP_Ix2y3z_Fx2z_a-1*I_ESP_Gx3z_Fx2z;
    abcd[iGrid*630+329] = 2.0E0*I_ESP_Ixy4z_Fx2z_a;
    abcd[iGrid*630+330] = 2.0E0*I_ESP_I6y_Fx2z_a-5*I_ESP_G4y_Fx2z;
    abcd[iGrid*630+331] = 2.0E0*I_ESP_I5yz_Fx2z_a-4*I_ESP_G3yz_Fx2z;
    abcd[iGrid*630+332] = 2.0E0*I_ESP_I4y2z_Fx2z_a-3*I_ESP_G2y2z_Fx2z;
    abcd[iGrid*630+333] = 2.0E0*I_ESP_I3y3z_Fx2z_a-2*I_ESP_Gy3z_Fx2z;
    abcd[iGrid*630+334] = 2.0E0*I_ESP_I2y4z_Fx2z_a-1*I_ESP_G4z_Fx2z;
    abcd[iGrid*630+335] = 2.0E0*I_ESP_Iy5z_Fx2z_a;
    abcd[iGrid*630+336] = 2.0E0*I_ESP_I5xy_F3y_a;
    abcd[iGrid*630+337] = 2.0E0*I_ESP_I4x2y_F3y_a-1*I_ESP_G4x_F3y;
    abcd[iGrid*630+338] = 2.0E0*I_ESP_I4xyz_F3y_a;
    abcd[iGrid*630+339] = 2.0E0*I_ESP_I3x3y_F3y_a-2*I_ESP_G3xy_F3y;
    abcd[iGrid*630+340] = 2.0E0*I_ESP_I3x2yz_F3y_a-1*I_ESP_G3xz_F3y;
    abcd[iGrid*630+341] = 2.0E0*I_ESP_I3xy2z_F3y_a;
    abcd[iGrid*630+342] = 2.0E0*I_ESP_I2x4y_F3y_a-3*I_ESP_G2x2y_F3y;
    abcd[iGrid*630+343] = 2.0E0*I_ESP_I2x3yz_F3y_a-2*I_ESP_G2xyz_F3y;
    abcd[iGrid*630+344] = 2.0E0*I_ESP_I2x2y2z_F3y_a-1*I_ESP_G2x2z_F3y;
    abcd[iGrid*630+345] = 2.0E0*I_ESP_I2xy3z_F3y_a;
    abcd[iGrid*630+346] = 2.0E0*I_ESP_Ix5y_F3y_a-4*I_ESP_Gx3y_F3y;
    abcd[iGrid*630+347] = 2.0E0*I_ESP_Ix4yz_F3y_a-3*I_ESP_Gx2yz_F3y;
    abcd[iGrid*630+348] = 2.0E0*I_ESP_Ix3y2z_F3y_a-2*I_ESP_Gxy2z_F3y;
    abcd[iGrid*630+349] = 2.0E0*I_ESP_Ix2y3z_F3y_a-1*I_ESP_Gx3z_F3y;
    abcd[iGrid*630+350] = 2.0E0*I_ESP_Ixy4z_F3y_a;
    abcd[iGrid*630+351] = 2.0E0*I_ESP_I6y_F3y_a-5*I_ESP_G4y_F3y;
    abcd[iGrid*630+352] = 2.0E0*I_ESP_I5yz_F3y_a-4*I_ESP_G3yz_F3y;
    abcd[iGrid*630+353] = 2.0E0*I_ESP_I4y2z_F3y_a-3*I_ESP_G2y2z_F3y;
    abcd[iGrid*630+354] = 2.0E0*I_ESP_I3y3z_F3y_a-2*I_ESP_Gy3z_F3y;
    abcd[iGrid*630+355] = 2.0E0*I_ESP_I2y4z_F3y_a-1*I_ESP_G4z_F3y;
    abcd[iGrid*630+356] = 2.0E0*I_ESP_Iy5z_F3y_a;
    abcd[iGrid*630+357] = 2.0E0*I_ESP_I5xy_F2yz_a;
    abcd[iGrid*630+358] = 2.0E0*I_ESP_I4x2y_F2yz_a-1*I_ESP_G4x_F2yz;
    abcd[iGrid*630+359] = 2.0E0*I_ESP_I4xyz_F2yz_a;
    abcd[iGrid*630+360] = 2.0E0*I_ESP_I3x3y_F2yz_a-2*I_ESP_G3xy_F2yz;
    abcd[iGrid*630+361] = 2.0E0*I_ESP_I3x2yz_F2yz_a-1*I_ESP_G3xz_F2yz;
    abcd[iGrid*630+362] = 2.0E0*I_ESP_I3xy2z_F2yz_a;
    abcd[iGrid*630+363] = 2.0E0*I_ESP_I2x4y_F2yz_a-3*I_ESP_G2x2y_F2yz;
    abcd[iGrid*630+364] = 2.0E0*I_ESP_I2x3yz_F2yz_a-2*I_ESP_G2xyz_F2yz;
    abcd[iGrid*630+365] = 2.0E0*I_ESP_I2x2y2z_F2yz_a-1*I_ESP_G2x2z_F2yz;
    abcd[iGrid*630+366] = 2.0E0*I_ESP_I2xy3z_F2yz_a;
    abcd[iGrid*630+367] = 2.0E0*I_ESP_Ix5y_F2yz_a-4*I_ESP_Gx3y_F2yz;
    abcd[iGrid*630+368] = 2.0E0*I_ESP_Ix4yz_F2yz_a-3*I_ESP_Gx2yz_F2yz;
    abcd[iGrid*630+369] = 2.0E0*I_ESP_Ix3y2z_F2yz_a-2*I_ESP_Gxy2z_F2yz;
    abcd[iGrid*630+370] = 2.0E0*I_ESP_Ix2y3z_F2yz_a-1*I_ESP_Gx3z_F2yz;
    abcd[iGrid*630+371] = 2.0E0*I_ESP_Ixy4z_F2yz_a;
    abcd[iGrid*630+372] = 2.0E0*I_ESP_I6y_F2yz_a-5*I_ESP_G4y_F2yz;
    abcd[iGrid*630+373] = 2.0E0*I_ESP_I5yz_F2yz_a-4*I_ESP_G3yz_F2yz;
    abcd[iGrid*630+374] = 2.0E0*I_ESP_I4y2z_F2yz_a-3*I_ESP_G2y2z_F2yz;
    abcd[iGrid*630+375] = 2.0E0*I_ESP_I3y3z_F2yz_a-2*I_ESP_Gy3z_F2yz;
    abcd[iGrid*630+376] = 2.0E0*I_ESP_I2y4z_F2yz_a-1*I_ESP_G4z_F2yz;
    abcd[iGrid*630+377] = 2.0E0*I_ESP_Iy5z_F2yz_a;
    abcd[iGrid*630+378] = 2.0E0*I_ESP_I5xy_Fy2z_a;
    abcd[iGrid*630+379] = 2.0E0*I_ESP_I4x2y_Fy2z_a-1*I_ESP_G4x_Fy2z;
    abcd[iGrid*630+380] = 2.0E0*I_ESP_I4xyz_Fy2z_a;
    abcd[iGrid*630+381] = 2.0E0*I_ESP_I3x3y_Fy2z_a-2*I_ESP_G3xy_Fy2z;
    abcd[iGrid*630+382] = 2.0E0*I_ESP_I3x2yz_Fy2z_a-1*I_ESP_G3xz_Fy2z;
    abcd[iGrid*630+383] = 2.0E0*I_ESP_I3xy2z_Fy2z_a;
    abcd[iGrid*630+384] = 2.0E0*I_ESP_I2x4y_Fy2z_a-3*I_ESP_G2x2y_Fy2z;
    abcd[iGrid*630+385] = 2.0E0*I_ESP_I2x3yz_Fy2z_a-2*I_ESP_G2xyz_Fy2z;
    abcd[iGrid*630+386] = 2.0E0*I_ESP_I2x2y2z_Fy2z_a-1*I_ESP_G2x2z_Fy2z;
    abcd[iGrid*630+387] = 2.0E0*I_ESP_I2xy3z_Fy2z_a;
    abcd[iGrid*630+388] = 2.0E0*I_ESP_Ix5y_Fy2z_a-4*I_ESP_Gx3y_Fy2z;
    abcd[iGrid*630+389] = 2.0E0*I_ESP_Ix4yz_Fy2z_a-3*I_ESP_Gx2yz_Fy2z;
    abcd[iGrid*630+390] = 2.0E0*I_ESP_Ix3y2z_Fy2z_a-2*I_ESP_Gxy2z_Fy2z;
    abcd[iGrid*630+391] = 2.0E0*I_ESP_Ix2y3z_Fy2z_a-1*I_ESP_Gx3z_Fy2z;
    abcd[iGrid*630+392] = 2.0E0*I_ESP_Ixy4z_Fy2z_a;
    abcd[iGrid*630+393] = 2.0E0*I_ESP_I6y_Fy2z_a-5*I_ESP_G4y_Fy2z;
    abcd[iGrid*630+394] = 2.0E0*I_ESP_I5yz_Fy2z_a-4*I_ESP_G3yz_Fy2z;
    abcd[iGrid*630+395] = 2.0E0*I_ESP_I4y2z_Fy2z_a-3*I_ESP_G2y2z_Fy2z;
    abcd[iGrid*630+396] = 2.0E0*I_ESP_I3y3z_Fy2z_a-2*I_ESP_Gy3z_Fy2z;
    abcd[iGrid*630+397] = 2.0E0*I_ESP_I2y4z_Fy2z_a-1*I_ESP_G4z_Fy2z;
    abcd[iGrid*630+398] = 2.0E0*I_ESP_Iy5z_Fy2z_a;
    abcd[iGrid*630+399] = 2.0E0*I_ESP_I5xy_F3z_a;
    abcd[iGrid*630+400] = 2.0E0*I_ESP_I4x2y_F3z_a-1*I_ESP_G4x_F3z;
    abcd[iGrid*630+401] = 2.0E0*I_ESP_I4xyz_F3z_a;
    abcd[iGrid*630+402] = 2.0E0*I_ESP_I3x3y_F3z_a-2*I_ESP_G3xy_F3z;
    abcd[iGrid*630+403] = 2.0E0*I_ESP_I3x2yz_F3z_a-1*I_ESP_G3xz_F3z;
    abcd[iGrid*630+404] = 2.0E0*I_ESP_I3xy2z_F3z_a;
    abcd[iGrid*630+405] = 2.0E0*I_ESP_I2x4y_F3z_a-3*I_ESP_G2x2y_F3z;
    abcd[iGrid*630+406] = 2.0E0*I_ESP_I2x3yz_F3z_a-2*I_ESP_G2xyz_F3z;
    abcd[iGrid*630+407] = 2.0E0*I_ESP_I2x2y2z_F3z_a-1*I_ESP_G2x2z_F3z;
    abcd[iGrid*630+408] = 2.0E0*I_ESP_I2xy3z_F3z_a;
    abcd[iGrid*630+409] = 2.0E0*I_ESP_Ix5y_F3z_a-4*I_ESP_Gx3y_F3z;
    abcd[iGrid*630+410] = 2.0E0*I_ESP_Ix4yz_F3z_a-3*I_ESP_Gx2yz_F3z;
    abcd[iGrid*630+411] = 2.0E0*I_ESP_Ix3y2z_F3z_a-2*I_ESP_Gxy2z_F3z;
    abcd[iGrid*630+412] = 2.0E0*I_ESP_Ix2y3z_F3z_a-1*I_ESP_Gx3z_F3z;
    abcd[iGrid*630+413] = 2.0E0*I_ESP_Ixy4z_F3z_a;
    abcd[iGrid*630+414] = 2.0E0*I_ESP_I6y_F3z_a-5*I_ESP_G4y_F3z;
    abcd[iGrid*630+415] = 2.0E0*I_ESP_I5yz_F3z_a-4*I_ESP_G3yz_F3z;
    abcd[iGrid*630+416] = 2.0E0*I_ESP_I4y2z_F3z_a-3*I_ESP_G2y2z_F3z;
    abcd[iGrid*630+417] = 2.0E0*I_ESP_I3y3z_F3z_a-2*I_ESP_Gy3z_F3z;
    abcd[iGrid*630+418] = 2.0E0*I_ESP_I2y4z_F3z_a-1*I_ESP_G4z_F3z;
    abcd[iGrid*630+419] = 2.0E0*I_ESP_Iy5z_F3z_a;

    /************************************************************
     * shell quartet name: SQ_ESP_H_F_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_I_F_a
     * RHS shell quartet name: SQ_ESP_G_F
     ************************************************************/
    abcd[iGrid*630+420] = 2.0E0*I_ESP_I5xz_F3x_a;
    abcd[iGrid*630+421] = 2.0E0*I_ESP_I4xyz_F3x_a;
    abcd[iGrid*630+422] = 2.0E0*I_ESP_I4x2z_F3x_a-1*I_ESP_G4x_F3x;
    abcd[iGrid*630+423] = 2.0E0*I_ESP_I3x2yz_F3x_a;
    abcd[iGrid*630+424] = 2.0E0*I_ESP_I3xy2z_F3x_a-1*I_ESP_G3xy_F3x;
    abcd[iGrid*630+425] = 2.0E0*I_ESP_I3x3z_F3x_a-2*I_ESP_G3xz_F3x;
    abcd[iGrid*630+426] = 2.0E0*I_ESP_I2x3yz_F3x_a;
    abcd[iGrid*630+427] = 2.0E0*I_ESP_I2x2y2z_F3x_a-1*I_ESP_G2x2y_F3x;
    abcd[iGrid*630+428] = 2.0E0*I_ESP_I2xy3z_F3x_a-2*I_ESP_G2xyz_F3x;
    abcd[iGrid*630+429] = 2.0E0*I_ESP_I2x4z_F3x_a-3*I_ESP_G2x2z_F3x;
    abcd[iGrid*630+430] = 2.0E0*I_ESP_Ix4yz_F3x_a;
    abcd[iGrid*630+431] = 2.0E0*I_ESP_Ix3y2z_F3x_a-1*I_ESP_Gx3y_F3x;
    abcd[iGrid*630+432] = 2.0E0*I_ESP_Ix2y3z_F3x_a-2*I_ESP_Gx2yz_F3x;
    abcd[iGrid*630+433] = 2.0E0*I_ESP_Ixy4z_F3x_a-3*I_ESP_Gxy2z_F3x;
    abcd[iGrid*630+434] = 2.0E0*I_ESP_Ix5z_F3x_a-4*I_ESP_Gx3z_F3x;
    abcd[iGrid*630+435] = 2.0E0*I_ESP_I5yz_F3x_a;
    abcd[iGrid*630+436] = 2.0E0*I_ESP_I4y2z_F3x_a-1*I_ESP_G4y_F3x;
    abcd[iGrid*630+437] = 2.0E0*I_ESP_I3y3z_F3x_a-2*I_ESP_G3yz_F3x;
    abcd[iGrid*630+438] = 2.0E0*I_ESP_I2y4z_F3x_a-3*I_ESP_G2y2z_F3x;
    abcd[iGrid*630+439] = 2.0E0*I_ESP_Iy5z_F3x_a-4*I_ESP_Gy3z_F3x;
    abcd[iGrid*630+440] = 2.0E0*I_ESP_I6z_F3x_a-5*I_ESP_G4z_F3x;
    abcd[iGrid*630+441] = 2.0E0*I_ESP_I5xz_F2xy_a;
    abcd[iGrid*630+442] = 2.0E0*I_ESP_I4xyz_F2xy_a;
    abcd[iGrid*630+443] = 2.0E0*I_ESP_I4x2z_F2xy_a-1*I_ESP_G4x_F2xy;
    abcd[iGrid*630+444] = 2.0E0*I_ESP_I3x2yz_F2xy_a;
    abcd[iGrid*630+445] = 2.0E0*I_ESP_I3xy2z_F2xy_a-1*I_ESP_G3xy_F2xy;
    abcd[iGrid*630+446] = 2.0E0*I_ESP_I3x3z_F2xy_a-2*I_ESP_G3xz_F2xy;
    abcd[iGrid*630+447] = 2.0E0*I_ESP_I2x3yz_F2xy_a;
    abcd[iGrid*630+448] = 2.0E0*I_ESP_I2x2y2z_F2xy_a-1*I_ESP_G2x2y_F2xy;
    abcd[iGrid*630+449] = 2.0E0*I_ESP_I2xy3z_F2xy_a-2*I_ESP_G2xyz_F2xy;
    abcd[iGrid*630+450] = 2.0E0*I_ESP_I2x4z_F2xy_a-3*I_ESP_G2x2z_F2xy;
    abcd[iGrid*630+451] = 2.0E0*I_ESP_Ix4yz_F2xy_a;
    abcd[iGrid*630+452] = 2.0E0*I_ESP_Ix3y2z_F2xy_a-1*I_ESP_Gx3y_F2xy;
    abcd[iGrid*630+453] = 2.0E0*I_ESP_Ix2y3z_F2xy_a-2*I_ESP_Gx2yz_F2xy;
    abcd[iGrid*630+454] = 2.0E0*I_ESP_Ixy4z_F2xy_a-3*I_ESP_Gxy2z_F2xy;
    abcd[iGrid*630+455] = 2.0E0*I_ESP_Ix5z_F2xy_a-4*I_ESP_Gx3z_F2xy;
    abcd[iGrid*630+456] = 2.0E0*I_ESP_I5yz_F2xy_a;
    abcd[iGrid*630+457] = 2.0E0*I_ESP_I4y2z_F2xy_a-1*I_ESP_G4y_F2xy;
    abcd[iGrid*630+458] = 2.0E0*I_ESP_I3y3z_F2xy_a-2*I_ESP_G3yz_F2xy;
    abcd[iGrid*630+459] = 2.0E0*I_ESP_I2y4z_F2xy_a-3*I_ESP_G2y2z_F2xy;
    abcd[iGrid*630+460] = 2.0E0*I_ESP_Iy5z_F2xy_a-4*I_ESP_Gy3z_F2xy;
    abcd[iGrid*630+461] = 2.0E0*I_ESP_I6z_F2xy_a-5*I_ESP_G4z_F2xy;
    abcd[iGrid*630+462] = 2.0E0*I_ESP_I5xz_F2xz_a;
    abcd[iGrid*630+463] = 2.0E0*I_ESP_I4xyz_F2xz_a;
    abcd[iGrid*630+464] = 2.0E0*I_ESP_I4x2z_F2xz_a-1*I_ESP_G4x_F2xz;
    abcd[iGrid*630+465] = 2.0E0*I_ESP_I3x2yz_F2xz_a;
    abcd[iGrid*630+466] = 2.0E0*I_ESP_I3xy2z_F2xz_a-1*I_ESP_G3xy_F2xz;
    abcd[iGrid*630+467] = 2.0E0*I_ESP_I3x3z_F2xz_a-2*I_ESP_G3xz_F2xz;
    abcd[iGrid*630+468] = 2.0E0*I_ESP_I2x3yz_F2xz_a;
    abcd[iGrid*630+469] = 2.0E0*I_ESP_I2x2y2z_F2xz_a-1*I_ESP_G2x2y_F2xz;
    abcd[iGrid*630+470] = 2.0E0*I_ESP_I2xy3z_F2xz_a-2*I_ESP_G2xyz_F2xz;
    abcd[iGrid*630+471] = 2.0E0*I_ESP_I2x4z_F2xz_a-3*I_ESP_G2x2z_F2xz;
    abcd[iGrid*630+472] = 2.0E0*I_ESP_Ix4yz_F2xz_a;
    abcd[iGrid*630+473] = 2.0E0*I_ESP_Ix3y2z_F2xz_a-1*I_ESP_Gx3y_F2xz;
    abcd[iGrid*630+474] = 2.0E0*I_ESP_Ix2y3z_F2xz_a-2*I_ESP_Gx2yz_F2xz;
    abcd[iGrid*630+475] = 2.0E0*I_ESP_Ixy4z_F2xz_a-3*I_ESP_Gxy2z_F2xz;
    abcd[iGrid*630+476] = 2.0E0*I_ESP_Ix5z_F2xz_a-4*I_ESP_Gx3z_F2xz;
    abcd[iGrid*630+477] = 2.0E0*I_ESP_I5yz_F2xz_a;
    abcd[iGrid*630+478] = 2.0E0*I_ESP_I4y2z_F2xz_a-1*I_ESP_G4y_F2xz;
    abcd[iGrid*630+479] = 2.0E0*I_ESP_I3y3z_F2xz_a-2*I_ESP_G3yz_F2xz;
    abcd[iGrid*630+480] = 2.0E0*I_ESP_I2y4z_F2xz_a-3*I_ESP_G2y2z_F2xz;
    abcd[iGrid*630+481] = 2.0E0*I_ESP_Iy5z_F2xz_a-4*I_ESP_Gy3z_F2xz;
    abcd[iGrid*630+482] = 2.0E0*I_ESP_I6z_F2xz_a-5*I_ESP_G4z_F2xz;
    abcd[iGrid*630+483] = 2.0E0*I_ESP_I5xz_Fx2y_a;
    abcd[iGrid*630+484] = 2.0E0*I_ESP_I4xyz_Fx2y_a;
    abcd[iGrid*630+485] = 2.0E0*I_ESP_I4x2z_Fx2y_a-1*I_ESP_G4x_Fx2y;
    abcd[iGrid*630+486] = 2.0E0*I_ESP_I3x2yz_Fx2y_a;
    abcd[iGrid*630+487] = 2.0E0*I_ESP_I3xy2z_Fx2y_a-1*I_ESP_G3xy_Fx2y;
    abcd[iGrid*630+488] = 2.0E0*I_ESP_I3x3z_Fx2y_a-2*I_ESP_G3xz_Fx2y;
    abcd[iGrid*630+489] = 2.0E0*I_ESP_I2x3yz_Fx2y_a;
    abcd[iGrid*630+490] = 2.0E0*I_ESP_I2x2y2z_Fx2y_a-1*I_ESP_G2x2y_Fx2y;
    abcd[iGrid*630+491] = 2.0E0*I_ESP_I2xy3z_Fx2y_a-2*I_ESP_G2xyz_Fx2y;
    abcd[iGrid*630+492] = 2.0E0*I_ESP_I2x4z_Fx2y_a-3*I_ESP_G2x2z_Fx2y;
    abcd[iGrid*630+493] = 2.0E0*I_ESP_Ix4yz_Fx2y_a;
    abcd[iGrid*630+494] = 2.0E0*I_ESP_Ix3y2z_Fx2y_a-1*I_ESP_Gx3y_Fx2y;
    abcd[iGrid*630+495] = 2.0E0*I_ESP_Ix2y3z_Fx2y_a-2*I_ESP_Gx2yz_Fx2y;
    abcd[iGrid*630+496] = 2.0E0*I_ESP_Ixy4z_Fx2y_a-3*I_ESP_Gxy2z_Fx2y;
    abcd[iGrid*630+497] = 2.0E0*I_ESP_Ix5z_Fx2y_a-4*I_ESP_Gx3z_Fx2y;
    abcd[iGrid*630+498] = 2.0E0*I_ESP_I5yz_Fx2y_a;
    abcd[iGrid*630+499] = 2.0E0*I_ESP_I4y2z_Fx2y_a-1*I_ESP_G4y_Fx2y;
    abcd[iGrid*630+500] = 2.0E0*I_ESP_I3y3z_Fx2y_a-2*I_ESP_G3yz_Fx2y;
    abcd[iGrid*630+501] = 2.0E0*I_ESP_I2y4z_Fx2y_a-3*I_ESP_G2y2z_Fx2y;
    abcd[iGrid*630+502] = 2.0E0*I_ESP_Iy5z_Fx2y_a-4*I_ESP_Gy3z_Fx2y;
    abcd[iGrid*630+503] = 2.0E0*I_ESP_I6z_Fx2y_a-5*I_ESP_G4z_Fx2y;
    abcd[iGrid*630+504] = 2.0E0*I_ESP_I5xz_Fxyz_a;
    abcd[iGrid*630+505] = 2.0E0*I_ESP_I4xyz_Fxyz_a;
    abcd[iGrid*630+506] = 2.0E0*I_ESP_I4x2z_Fxyz_a-1*I_ESP_G4x_Fxyz;
    abcd[iGrid*630+507] = 2.0E0*I_ESP_I3x2yz_Fxyz_a;
    abcd[iGrid*630+508] = 2.0E0*I_ESP_I3xy2z_Fxyz_a-1*I_ESP_G3xy_Fxyz;
    abcd[iGrid*630+509] = 2.0E0*I_ESP_I3x3z_Fxyz_a-2*I_ESP_G3xz_Fxyz;
    abcd[iGrid*630+510] = 2.0E0*I_ESP_I2x3yz_Fxyz_a;
    abcd[iGrid*630+511] = 2.0E0*I_ESP_I2x2y2z_Fxyz_a-1*I_ESP_G2x2y_Fxyz;
    abcd[iGrid*630+512] = 2.0E0*I_ESP_I2xy3z_Fxyz_a-2*I_ESP_G2xyz_Fxyz;
    abcd[iGrid*630+513] = 2.0E0*I_ESP_I2x4z_Fxyz_a-3*I_ESP_G2x2z_Fxyz;
    abcd[iGrid*630+514] = 2.0E0*I_ESP_Ix4yz_Fxyz_a;
    abcd[iGrid*630+515] = 2.0E0*I_ESP_Ix3y2z_Fxyz_a-1*I_ESP_Gx3y_Fxyz;
    abcd[iGrid*630+516] = 2.0E0*I_ESP_Ix2y3z_Fxyz_a-2*I_ESP_Gx2yz_Fxyz;
    abcd[iGrid*630+517] = 2.0E0*I_ESP_Ixy4z_Fxyz_a-3*I_ESP_Gxy2z_Fxyz;
    abcd[iGrid*630+518] = 2.0E0*I_ESP_Ix5z_Fxyz_a-4*I_ESP_Gx3z_Fxyz;
    abcd[iGrid*630+519] = 2.0E0*I_ESP_I5yz_Fxyz_a;
    abcd[iGrid*630+520] = 2.0E0*I_ESP_I4y2z_Fxyz_a-1*I_ESP_G4y_Fxyz;
    abcd[iGrid*630+521] = 2.0E0*I_ESP_I3y3z_Fxyz_a-2*I_ESP_G3yz_Fxyz;
    abcd[iGrid*630+522] = 2.0E0*I_ESP_I2y4z_Fxyz_a-3*I_ESP_G2y2z_Fxyz;
    abcd[iGrid*630+523] = 2.0E0*I_ESP_Iy5z_Fxyz_a-4*I_ESP_Gy3z_Fxyz;
    abcd[iGrid*630+524] = 2.0E0*I_ESP_I6z_Fxyz_a-5*I_ESP_G4z_Fxyz;
    abcd[iGrid*630+525] = 2.0E0*I_ESP_I5xz_Fx2z_a;
    abcd[iGrid*630+526] = 2.0E0*I_ESP_I4xyz_Fx2z_a;
    abcd[iGrid*630+527] = 2.0E0*I_ESP_I4x2z_Fx2z_a-1*I_ESP_G4x_Fx2z;
    abcd[iGrid*630+528] = 2.0E0*I_ESP_I3x2yz_Fx2z_a;
    abcd[iGrid*630+529] = 2.0E0*I_ESP_I3xy2z_Fx2z_a-1*I_ESP_G3xy_Fx2z;
    abcd[iGrid*630+530] = 2.0E0*I_ESP_I3x3z_Fx2z_a-2*I_ESP_G3xz_Fx2z;
    abcd[iGrid*630+531] = 2.0E0*I_ESP_I2x3yz_Fx2z_a;
    abcd[iGrid*630+532] = 2.0E0*I_ESP_I2x2y2z_Fx2z_a-1*I_ESP_G2x2y_Fx2z;
    abcd[iGrid*630+533] = 2.0E0*I_ESP_I2xy3z_Fx2z_a-2*I_ESP_G2xyz_Fx2z;
    abcd[iGrid*630+534] = 2.0E0*I_ESP_I2x4z_Fx2z_a-3*I_ESP_G2x2z_Fx2z;
    abcd[iGrid*630+535] = 2.0E0*I_ESP_Ix4yz_Fx2z_a;
    abcd[iGrid*630+536] = 2.0E0*I_ESP_Ix3y2z_Fx2z_a-1*I_ESP_Gx3y_Fx2z;
    abcd[iGrid*630+537] = 2.0E0*I_ESP_Ix2y3z_Fx2z_a-2*I_ESP_Gx2yz_Fx2z;
    abcd[iGrid*630+538] = 2.0E0*I_ESP_Ixy4z_Fx2z_a-3*I_ESP_Gxy2z_Fx2z;
    abcd[iGrid*630+539] = 2.0E0*I_ESP_Ix5z_Fx2z_a-4*I_ESP_Gx3z_Fx2z;
    abcd[iGrid*630+540] = 2.0E0*I_ESP_I5yz_Fx2z_a;
    abcd[iGrid*630+541] = 2.0E0*I_ESP_I4y2z_Fx2z_a-1*I_ESP_G4y_Fx2z;
    abcd[iGrid*630+542] = 2.0E0*I_ESP_I3y3z_Fx2z_a-2*I_ESP_G3yz_Fx2z;
    abcd[iGrid*630+543] = 2.0E0*I_ESP_I2y4z_Fx2z_a-3*I_ESP_G2y2z_Fx2z;
    abcd[iGrid*630+544] = 2.0E0*I_ESP_Iy5z_Fx2z_a-4*I_ESP_Gy3z_Fx2z;
    abcd[iGrid*630+545] = 2.0E0*I_ESP_I6z_Fx2z_a-5*I_ESP_G4z_Fx2z;
    abcd[iGrid*630+546] = 2.0E0*I_ESP_I5xz_F3y_a;
    abcd[iGrid*630+547] = 2.0E0*I_ESP_I4xyz_F3y_a;
    abcd[iGrid*630+548] = 2.0E0*I_ESP_I4x2z_F3y_a-1*I_ESP_G4x_F3y;
    abcd[iGrid*630+549] = 2.0E0*I_ESP_I3x2yz_F3y_a;
    abcd[iGrid*630+550] = 2.0E0*I_ESP_I3xy2z_F3y_a-1*I_ESP_G3xy_F3y;
    abcd[iGrid*630+551] = 2.0E0*I_ESP_I3x3z_F3y_a-2*I_ESP_G3xz_F3y;
    abcd[iGrid*630+552] = 2.0E0*I_ESP_I2x3yz_F3y_a;
    abcd[iGrid*630+553] = 2.0E0*I_ESP_I2x2y2z_F3y_a-1*I_ESP_G2x2y_F3y;
    abcd[iGrid*630+554] = 2.0E0*I_ESP_I2xy3z_F3y_a-2*I_ESP_G2xyz_F3y;
    abcd[iGrid*630+555] = 2.0E0*I_ESP_I2x4z_F3y_a-3*I_ESP_G2x2z_F3y;
    abcd[iGrid*630+556] = 2.0E0*I_ESP_Ix4yz_F3y_a;
    abcd[iGrid*630+557] = 2.0E0*I_ESP_Ix3y2z_F3y_a-1*I_ESP_Gx3y_F3y;
    abcd[iGrid*630+558] = 2.0E0*I_ESP_Ix2y3z_F3y_a-2*I_ESP_Gx2yz_F3y;
    abcd[iGrid*630+559] = 2.0E0*I_ESP_Ixy4z_F3y_a-3*I_ESP_Gxy2z_F3y;
    abcd[iGrid*630+560] = 2.0E0*I_ESP_Ix5z_F3y_a-4*I_ESP_Gx3z_F3y;
    abcd[iGrid*630+561] = 2.0E0*I_ESP_I5yz_F3y_a;
    abcd[iGrid*630+562] = 2.0E0*I_ESP_I4y2z_F3y_a-1*I_ESP_G4y_F3y;
    abcd[iGrid*630+563] = 2.0E0*I_ESP_I3y3z_F3y_a-2*I_ESP_G3yz_F3y;
    abcd[iGrid*630+564] = 2.0E0*I_ESP_I2y4z_F3y_a-3*I_ESP_G2y2z_F3y;
    abcd[iGrid*630+565] = 2.0E0*I_ESP_Iy5z_F3y_a-4*I_ESP_Gy3z_F3y;
    abcd[iGrid*630+566] = 2.0E0*I_ESP_I6z_F3y_a-5*I_ESP_G4z_F3y;
    abcd[iGrid*630+567] = 2.0E0*I_ESP_I5xz_F2yz_a;
    abcd[iGrid*630+568] = 2.0E0*I_ESP_I4xyz_F2yz_a;
    abcd[iGrid*630+569] = 2.0E0*I_ESP_I4x2z_F2yz_a-1*I_ESP_G4x_F2yz;
    abcd[iGrid*630+570] = 2.0E0*I_ESP_I3x2yz_F2yz_a;
    abcd[iGrid*630+571] = 2.0E0*I_ESP_I3xy2z_F2yz_a-1*I_ESP_G3xy_F2yz;
    abcd[iGrid*630+572] = 2.0E0*I_ESP_I3x3z_F2yz_a-2*I_ESP_G3xz_F2yz;
    abcd[iGrid*630+573] = 2.0E0*I_ESP_I2x3yz_F2yz_a;
    abcd[iGrid*630+574] = 2.0E0*I_ESP_I2x2y2z_F2yz_a-1*I_ESP_G2x2y_F2yz;
    abcd[iGrid*630+575] = 2.0E0*I_ESP_I2xy3z_F2yz_a-2*I_ESP_G2xyz_F2yz;
    abcd[iGrid*630+576] = 2.0E0*I_ESP_I2x4z_F2yz_a-3*I_ESP_G2x2z_F2yz;
    abcd[iGrid*630+577] = 2.0E0*I_ESP_Ix4yz_F2yz_a;
    abcd[iGrid*630+578] = 2.0E0*I_ESP_Ix3y2z_F2yz_a-1*I_ESP_Gx3y_F2yz;
    abcd[iGrid*630+579] = 2.0E0*I_ESP_Ix2y3z_F2yz_a-2*I_ESP_Gx2yz_F2yz;
    abcd[iGrid*630+580] = 2.0E0*I_ESP_Ixy4z_F2yz_a-3*I_ESP_Gxy2z_F2yz;
    abcd[iGrid*630+581] = 2.0E0*I_ESP_Ix5z_F2yz_a-4*I_ESP_Gx3z_F2yz;
    abcd[iGrid*630+582] = 2.0E0*I_ESP_I5yz_F2yz_a;
    abcd[iGrid*630+583] = 2.0E0*I_ESP_I4y2z_F2yz_a-1*I_ESP_G4y_F2yz;
    abcd[iGrid*630+584] = 2.0E0*I_ESP_I3y3z_F2yz_a-2*I_ESP_G3yz_F2yz;
    abcd[iGrid*630+585] = 2.0E0*I_ESP_I2y4z_F2yz_a-3*I_ESP_G2y2z_F2yz;
    abcd[iGrid*630+586] = 2.0E0*I_ESP_Iy5z_F2yz_a-4*I_ESP_Gy3z_F2yz;
    abcd[iGrid*630+587] = 2.0E0*I_ESP_I6z_F2yz_a-5*I_ESP_G4z_F2yz;
    abcd[iGrid*630+588] = 2.0E0*I_ESP_I5xz_Fy2z_a;
    abcd[iGrid*630+589] = 2.0E0*I_ESP_I4xyz_Fy2z_a;
    abcd[iGrid*630+590] = 2.0E0*I_ESP_I4x2z_Fy2z_a-1*I_ESP_G4x_Fy2z;
    abcd[iGrid*630+591] = 2.0E0*I_ESP_I3x2yz_Fy2z_a;
    abcd[iGrid*630+592] = 2.0E0*I_ESP_I3xy2z_Fy2z_a-1*I_ESP_G3xy_Fy2z;
    abcd[iGrid*630+593] = 2.0E0*I_ESP_I3x3z_Fy2z_a-2*I_ESP_G3xz_Fy2z;
    abcd[iGrid*630+594] = 2.0E0*I_ESP_I2x3yz_Fy2z_a;
    abcd[iGrid*630+595] = 2.0E0*I_ESP_I2x2y2z_Fy2z_a-1*I_ESP_G2x2y_Fy2z;
    abcd[iGrid*630+596] = 2.0E0*I_ESP_I2xy3z_Fy2z_a-2*I_ESP_G2xyz_Fy2z;
    abcd[iGrid*630+597] = 2.0E0*I_ESP_I2x4z_Fy2z_a-3*I_ESP_G2x2z_Fy2z;
    abcd[iGrid*630+598] = 2.0E0*I_ESP_Ix4yz_Fy2z_a;
    abcd[iGrid*630+599] = 2.0E0*I_ESP_Ix3y2z_Fy2z_a-1*I_ESP_Gx3y_Fy2z;
    abcd[iGrid*630+600] = 2.0E0*I_ESP_Ix2y3z_Fy2z_a-2*I_ESP_Gx2yz_Fy2z;
    abcd[iGrid*630+601] = 2.0E0*I_ESP_Ixy4z_Fy2z_a-3*I_ESP_Gxy2z_Fy2z;
    abcd[iGrid*630+602] = 2.0E0*I_ESP_Ix5z_Fy2z_a-4*I_ESP_Gx3z_Fy2z;
    abcd[iGrid*630+603] = 2.0E0*I_ESP_I5yz_Fy2z_a;
    abcd[iGrid*630+604] = 2.0E0*I_ESP_I4y2z_Fy2z_a-1*I_ESP_G4y_Fy2z;
    abcd[iGrid*630+605] = 2.0E0*I_ESP_I3y3z_Fy2z_a-2*I_ESP_G3yz_Fy2z;
    abcd[iGrid*630+606] = 2.0E0*I_ESP_I2y4z_Fy2z_a-3*I_ESP_G2y2z_Fy2z;
    abcd[iGrid*630+607] = 2.0E0*I_ESP_Iy5z_Fy2z_a-4*I_ESP_Gy3z_Fy2z;
    abcd[iGrid*630+608] = 2.0E0*I_ESP_I6z_Fy2z_a-5*I_ESP_G4z_Fy2z;
    abcd[iGrid*630+609] = 2.0E0*I_ESP_I5xz_F3z_a;
    abcd[iGrid*630+610] = 2.0E0*I_ESP_I4xyz_F3z_a;
    abcd[iGrid*630+611] = 2.0E0*I_ESP_I4x2z_F3z_a-1*I_ESP_G4x_F3z;
    abcd[iGrid*630+612] = 2.0E0*I_ESP_I3x2yz_F3z_a;
    abcd[iGrid*630+613] = 2.0E0*I_ESP_I3xy2z_F3z_a-1*I_ESP_G3xy_F3z;
    abcd[iGrid*630+614] = 2.0E0*I_ESP_I3x3z_F3z_a-2*I_ESP_G3xz_F3z;
    abcd[iGrid*630+615] = 2.0E0*I_ESP_I2x3yz_F3z_a;
    abcd[iGrid*630+616] = 2.0E0*I_ESP_I2x2y2z_F3z_a-1*I_ESP_G2x2y_F3z;
    abcd[iGrid*630+617] = 2.0E0*I_ESP_I2xy3z_F3z_a-2*I_ESP_G2xyz_F3z;
    abcd[iGrid*630+618] = 2.0E0*I_ESP_I2x4z_F3z_a-3*I_ESP_G2x2z_F3z;
    abcd[iGrid*630+619] = 2.0E0*I_ESP_Ix4yz_F3z_a;
    abcd[iGrid*630+620] = 2.0E0*I_ESP_Ix3y2z_F3z_a-1*I_ESP_Gx3y_F3z;
    abcd[iGrid*630+621] = 2.0E0*I_ESP_Ix2y3z_F3z_a-2*I_ESP_Gx2yz_F3z;
    abcd[iGrid*630+622] = 2.0E0*I_ESP_Ixy4z_F3z_a-3*I_ESP_Gxy2z_F3z;
    abcd[iGrid*630+623] = 2.0E0*I_ESP_Ix5z_F3z_a-4*I_ESP_Gx3z_F3z;
    abcd[iGrid*630+624] = 2.0E0*I_ESP_I5yz_F3z_a;
    abcd[iGrid*630+625] = 2.0E0*I_ESP_I4y2z_F3z_a-1*I_ESP_G4y_F3z;
    abcd[iGrid*630+626] = 2.0E0*I_ESP_I3y3z_F3z_a-2*I_ESP_G3yz_F3z;
    abcd[iGrid*630+627] = 2.0E0*I_ESP_I2y4z_F3z_a-3*I_ESP_G2y2z_F3z;
    abcd[iGrid*630+628] = 2.0E0*I_ESP_Iy5z_F3z_a-4*I_ESP_Gy3z_F3z;
    abcd[iGrid*630+629] = 2.0E0*I_ESP_I6z_F3z_a-5*I_ESP_G4z_F3z;
  }
}
