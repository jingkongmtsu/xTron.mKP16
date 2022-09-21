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
// BRA1 as redundant position, total RHS integrals evaluated as: 3816
// BRA2 as redundant position, total RHS integrals evaluated as: 3048
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

void hgp_os_twobodyoverlap_g_d_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_L8x_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xy_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6xyz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x2yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5xy2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x3yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x2y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4xy3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x4yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x3y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x2y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3xy4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x5yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x4y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x3y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x2y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2xy5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx6yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx5y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx4y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx3y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx2y5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lxy6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L8y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L7yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3y5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2y6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ly7z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L8z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7x_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xy_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5xyz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x2yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4xy2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x3yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x2y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3xy3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x4yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x3y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x2y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2xy4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx5yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx4y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx3y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx2y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kxy5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2y5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ky6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6x_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F3x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Fxyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Fy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_D2x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_D2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Dyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_D2z_S = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double fbra  = ifac[ip2];
    Double onedz = iexp[ip2];
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;
    if(fabs(I_TWOBODYOVERLAP_S_S_vrr)<THRESHOLD_MATH) continue;


    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_S_vrr = PAX*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Py_S_vrr = PAY*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_S_vrr = PAZ*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_S_vrr = PAY*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_S_vrr = PAZ*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_S_vrr = PAZ*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_S_vrr = PAZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_S_vrr = PAX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_S_vrr = PAY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_S_vrr = PAX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_S_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_S_vrr = PAX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_S_vrr = PAY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_S_vrr = PAY*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_S_vrr = PAZ*I_TWOBODYOVERLAP_D2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_S_vrr = PAX*I_TWOBODYOVERLAP_F3x_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_S_vrr = PAY*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_S_vrr = PAY*I_TWOBODYOVERLAP_F2xy_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_S_vrr = PAX*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_S_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_S_vrr = PAX*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_S_vrr = PAY*I_TWOBODYOVERLAP_F3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_S_vrr = PAY*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_S_vrr = PAZ*I_TWOBODYOVERLAP_F3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_S_vrr = PAX*I_TWOBODYOVERLAP_G4x_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_S_vrr = PAY*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_S_vrr = PAZ*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_S_vrr = PAY*I_TWOBODYOVERLAP_G3xy_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_S_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_S_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_S_vrr = PAX*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_S_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_S_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_S_vrr = PAX*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_S_vrr = PAY*I_TWOBODYOVERLAP_G4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_S_vrr = PAY*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_S_vrr = PAZ*I_TWOBODYOVERLAP_G4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_S_vrr = PAX*I_TWOBODYOVERLAP_H5x_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_S_vrr = PAY*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_S_vrr = PAZ*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_S_vrr = PAY*I_TWOBODYOVERLAP_H4xy_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_S_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_S_vrr = PAX*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_S_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_S_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_S_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_S_vrr = PAX*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_S_vrr = PAY*I_TWOBODYOVERLAP_H5y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_S_vrr = PAY*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_S_vrr = PAZ*I_TWOBODYOVERLAP_H5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_K7x_S_vrr = PAX*I_TWOBODYOVERLAP_I6x_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xy_S_vrr = PAY*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xz_S_vrr = PAZ*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_S_vrr = PAY*I_TWOBODYOVERLAP_I5xy_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I5xz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_S_vrr = PAY*I_TWOBODYOVERLAP_I4x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_I4x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_S_vrr = PAX*I_TWOBODYOVERLAP_I2x4y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I3x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_S_vrr = PAX*I_TWOBODYOVERLAP_I2x4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_S_vrr = PAX*I_TWOBODYOVERLAP_Ix5y_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I2x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_I2xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_S_vrr = PAY*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_S_vrr = PAX*I_TWOBODYOVERLAP_Ix5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_S_vrr = PAX*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_S_vrr = PAX*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_S_vrr = PAX*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_S_vrr = PAX*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_S_vrr = PAY*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_S_vrr = PAX*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7y_S_vrr = PAY*I_TWOBODYOVERLAP_I6y_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_K6yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I5yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_I4y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_S_vrr = PAY*I_TWOBODYOVERLAP_I2y4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_S_vrr = PAY*I_TWOBODYOVERLAP_Iy5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_S_vrr = PAY*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7z_S_vrr = PAZ*I_TWOBODYOVERLAP_I6z_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_L_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_L8x_S_vrr = PAX*I_TWOBODYOVERLAP_K7x_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L7xy_S_vrr = PAY*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_L7xz_S_vrr = PAZ*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_L6x2y_S_vrr = PAY*I_TWOBODYOVERLAP_K6xy_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L6xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_K6xy_S_vrr;
    Double I_TWOBODYOVERLAP_L6x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K6xz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L5x3y_S_vrr = PAY*I_TWOBODYOVERLAP_K5x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_L5x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L5xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    Double I_TWOBODYOVERLAP_L5x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_K5x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_L4x4y_S_vrr = PAY*I_TWOBODYOVERLAP_K4x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L4x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    Double I_TWOBODYOVERLAP_L4x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L4xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    Double I_TWOBODYOVERLAP_L4x4z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_L3x5y_S_vrr = PAX*I_TWOBODYOVERLAP_K2x5y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K3x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_K3xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_L3xy4z_S_vrr = PAY*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    Double I_TWOBODYOVERLAP_L3x5z_S_vrr = PAX*I_TWOBODYOVERLAP_K2x5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x6y_S_vrr = PAX*I_TWOBODYOVERLAP_Kx6y_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K2x4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x3y3z_S_vrr = PAX*I_TWOBODYOVERLAP_Kx3y3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_K2xy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_L2xy5z_S_vrr = PAY*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x6z_S_vrr = PAX*I_TWOBODYOVERLAP_Kx6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx7y_S_vrr = PAX*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_Lx6yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    Double I_TWOBODYOVERLAP_Lx5y2z_S_vrr = PAX*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx4y3z_S_vrr = PAX*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx3y4z_S_vrr = PAX*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx2y5z_S_vrr = PAX*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    Double I_TWOBODYOVERLAP_Lxy6z_S_vrr = PAY*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx7z_S_vrr = PAX*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_L8y_S_vrr = PAY*I_TWOBODYOVERLAP_K7y_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L7yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_L6y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K6yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L5y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_K5y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_L4y4z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_L3y5z_S_vrr = PAY*I_TWOBODYOVERLAP_K2y5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2y6z_S_vrr = PAY*I_TWOBODYOVERLAP_Ky6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_Ly7z_S_vrr = PAY*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_L8z_S_vrr = PAZ*I_TWOBODYOVERLAP_K7z_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_L_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_L_S_aa_coefs = alpha*alpha;
    I_TWOBODYOVERLAP_L8x_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L8x_S_vrr;
    I_TWOBODYOVERLAP_L7xy_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L7xy_S_vrr;
    I_TWOBODYOVERLAP_L7xz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L7xz_S_vrr;
    I_TWOBODYOVERLAP_L6x2y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L6x2y_S_vrr;
    I_TWOBODYOVERLAP_L6xyz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L6xyz_S_vrr;
    I_TWOBODYOVERLAP_L6x2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L6x2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5x3y_S_vrr;
    I_TWOBODYOVERLAP_L5x2yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5x2yz_S_vrr;
    I_TWOBODYOVERLAP_L5xy2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5xy2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5x3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4x4y_S_vrr;
    I_TWOBODYOVERLAP_L4x3yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4x3yz_S_vrr;
    I_TWOBODYOVERLAP_L4x2y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4x2y2z_S_vrr;
    I_TWOBODYOVERLAP_L4xy3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4xy3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4x4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x5y_S_vrr;
    I_TWOBODYOVERLAP_L3x4yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x4yz_S_vrr;
    I_TWOBODYOVERLAP_L3x3y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x3y2z_S_vrr;
    I_TWOBODYOVERLAP_L3x2y3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x2y3z_S_vrr;
    I_TWOBODYOVERLAP_L3xy4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3xy4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x6y_S_vrr;
    I_TWOBODYOVERLAP_L2x5yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x5yz_S_vrr;
    I_TWOBODYOVERLAP_L2x4y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x4y2z_S_vrr;
    I_TWOBODYOVERLAP_L2x3y3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x3y3z_S_vrr;
    I_TWOBODYOVERLAP_L2x2y4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x2y4z_S_vrr;
    I_TWOBODYOVERLAP_L2xy5z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2xy5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx7y_S_vrr;
    I_TWOBODYOVERLAP_Lx6yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx6yz_S_vrr;
    I_TWOBODYOVERLAP_Lx5y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx5y2z_S_vrr;
    I_TWOBODYOVERLAP_Lx4y3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx4y3z_S_vrr;
    I_TWOBODYOVERLAP_Lx3y4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx3y4z_S_vrr;
    I_TWOBODYOVERLAP_Lx2y5z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx2y5z_S_vrr;
    I_TWOBODYOVERLAP_Lxy6z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lxy6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx7z_S_vrr;
    I_TWOBODYOVERLAP_L8y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L8y_S_vrr;
    I_TWOBODYOVERLAP_L7yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L7yz_S_vrr;
    I_TWOBODYOVERLAP_L6y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L6y2z_S_vrr;
    I_TWOBODYOVERLAP_L5y3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5y3z_S_vrr;
    I_TWOBODYOVERLAP_L4y4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4y4z_S_vrr;
    I_TWOBODYOVERLAP_L3y5z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3y5z_S_vrr;
    I_TWOBODYOVERLAP_L2y6z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2y6z_S_vrr;
    I_TWOBODYOVERLAP_Ly7z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Ly7z_S_vrr;
    I_TWOBODYOVERLAP_L8z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L8z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_K_S_aa_coefs = alpha*alpha;
    I_TWOBODYOVERLAP_K7x_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K7x_S_vrr;
    I_TWOBODYOVERLAP_K6xy_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K6xy_S_vrr;
    I_TWOBODYOVERLAP_K6xz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K6xz_S_vrr;
    I_TWOBODYOVERLAP_K5x2y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    I_TWOBODYOVERLAP_K5xyz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K5xyz_S_vrr;
    I_TWOBODYOVERLAP_K5x2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    I_TWOBODYOVERLAP_K4x2yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4x2yz_S_vrr;
    I_TWOBODYOVERLAP_K4xy2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4xy2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    I_TWOBODYOVERLAP_K3x3yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    I_TWOBODYOVERLAP_K3x2y2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3x2y2z_S_vrr;
    I_TWOBODYOVERLAP_K3xy3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3xy3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    I_TWOBODYOVERLAP_K2x4yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x4yz_S_vrr;
    I_TWOBODYOVERLAP_K2x3y2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x3y2z_S_vrr;
    I_TWOBODYOVERLAP_K2x2y3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x2y3z_S_vrr;
    I_TWOBODYOVERLAP_K2xy4z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2xy4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    I_TWOBODYOVERLAP_Kx5yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx5yz_S_vrr;
    I_TWOBODYOVERLAP_Kx4y2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx4y2z_S_vrr;
    I_TWOBODYOVERLAP_Kx3y3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx3y3z_S_vrr;
    I_TWOBODYOVERLAP_Kx2y4z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx2y4z_S_vrr;
    I_TWOBODYOVERLAP_Kxy5z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kxy5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    I_TWOBODYOVERLAP_K7y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K7y_S_vrr;
    I_TWOBODYOVERLAP_K6yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K6yz_S_vrr;
    I_TWOBODYOVERLAP_K5y2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    I_TWOBODYOVERLAP_K4y3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    I_TWOBODYOVERLAP_K3y4z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    I_TWOBODYOVERLAP_K2y5z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    I_TWOBODYOVERLAP_Ky6z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    I_TWOBODYOVERLAP_K7z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_I_S_aa_coefs = alpha*alpha;
    I_TWOBODYOVERLAP_I6x_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_I_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_I6x_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_H_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_H5x_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_G_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_G4x_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_G4x_S += I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S += I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S += I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S += I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S += I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S += I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S += I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S += I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S += I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S += I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S += I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S += I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S += I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S += I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S += I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_F3x_S += I_TWOBODYOVERLAP_F3x_S_vrr;
    I_TWOBODYOVERLAP_F2xy_S += I_TWOBODYOVERLAP_F2xy_S_vrr;
    I_TWOBODYOVERLAP_F2xz_S += I_TWOBODYOVERLAP_F2xz_S_vrr;
    I_TWOBODYOVERLAP_Fx2y_S += I_TWOBODYOVERLAP_Fx2y_S_vrr;
    I_TWOBODYOVERLAP_Fxyz_S += I_TWOBODYOVERLAP_Fxyz_S_vrr;
    I_TWOBODYOVERLAP_Fx2z_S += I_TWOBODYOVERLAP_Fx2z_S_vrr;
    I_TWOBODYOVERLAP_F3y_S += I_TWOBODYOVERLAP_F3y_S_vrr;
    I_TWOBODYOVERLAP_F2yz_S += I_TWOBODYOVERLAP_F2yz_S_vrr;
    I_TWOBODYOVERLAP_Fy2z_S += I_TWOBODYOVERLAP_Fy2z_S_vrr;
    I_TWOBODYOVERLAP_F3z_S += I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_D2x_S += I_TWOBODYOVERLAP_D2x_S_vrr;
    I_TWOBODYOVERLAP_Dxy_S += I_TWOBODYOVERLAP_Dxy_S_vrr;
    I_TWOBODYOVERLAP_Dxz_S += I_TWOBODYOVERLAP_Dxz_S_vrr;
    I_TWOBODYOVERLAP_D2y_S += I_TWOBODYOVERLAP_D2y_S_vrr;
    I_TWOBODYOVERLAP_Dyz_S += I_TWOBODYOVERLAP_Dyz_S_vrr;
    I_TWOBODYOVERLAP_D2z_S += I_TWOBODYOVERLAP_D2z_S_vrr;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_D_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_D2x_Px = I_TWOBODYOVERLAP_F3x_S+ABX*I_TWOBODYOVERLAP_D2x_S;
  Double I_TWOBODYOVERLAP_Dxy_Px = I_TWOBODYOVERLAP_F2xy_S+ABX*I_TWOBODYOVERLAP_Dxy_S;
  Double I_TWOBODYOVERLAP_Dxz_Px = I_TWOBODYOVERLAP_F2xz_S+ABX*I_TWOBODYOVERLAP_Dxz_S;
  Double I_TWOBODYOVERLAP_D2y_Px = I_TWOBODYOVERLAP_Fx2y_S+ABX*I_TWOBODYOVERLAP_D2y_S;
  Double I_TWOBODYOVERLAP_Dyz_Px = I_TWOBODYOVERLAP_Fxyz_S+ABX*I_TWOBODYOVERLAP_Dyz_S;
  Double I_TWOBODYOVERLAP_D2z_Px = I_TWOBODYOVERLAP_Fx2z_S+ABX*I_TWOBODYOVERLAP_D2z_S;
  Double I_TWOBODYOVERLAP_D2x_Py = I_TWOBODYOVERLAP_F2xy_S+ABY*I_TWOBODYOVERLAP_D2x_S;
  Double I_TWOBODYOVERLAP_Dxy_Py = I_TWOBODYOVERLAP_Fx2y_S+ABY*I_TWOBODYOVERLAP_Dxy_S;
  Double I_TWOBODYOVERLAP_Dxz_Py = I_TWOBODYOVERLAP_Fxyz_S+ABY*I_TWOBODYOVERLAP_Dxz_S;
  Double I_TWOBODYOVERLAP_D2y_Py = I_TWOBODYOVERLAP_F3y_S+ABY*I_TWOBODYOVERLAP_D2y_S;
  Double I_TWOBODYOVERLAP_Dyz_Py = I_TWOBODYOVERLAP_F2yz_S+ABY*I_TWOBODYOVERLAP_Dyz_S;
  Double I_TWOBODYOVERLAP_D2z_Py = I_TWOBODYOVERLAP_Fy2z_S+ABY*I_TWOBODYOVERLAP_D2z_S;
  Double I_TWOBODYOVERLAP_D2x_Pz = I_TWOBODYOVERLAP_F2xz_S+ABZ*I_TWOBODYOVERLAP_D2x_S;
  Double I_TWOBODYOVERLAP_Dxy_Pz = I_TWOBODYOVERLAP_Fxyz_S+ABZ*I_TWOBODYOVERLAP_Dxy_S;
  Double I_TWOBODYOVERLAP_Dxz_Pz = I_TWOBODYOVERLAP_Fx2z_S+ABZ*I_TWOBODYOVERLAP_Dxz_S;
  Double I_TWOBODYOVERLAP_D2y_Pz = I_TWOBODYOVERLAP_F2yz_S+ABZ*I_TWOBODYOVERLAP_D2y_S;
  Double I_TWOBODYOVERLAP_Dyz_Pz = I_TWOBODYOVERLAP_Fy2z_S+ABZ*I_TWOBODYOVERLAP_Dyz_S;
  Double I_TWOBODYOVERLAP_D2z_Pz = I_TWOBODYOVERLAP_F3z_S+ABZ*I_TWOBODYOVERLAP_D2z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 5 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_F3x_Px = I_TWOBODYOVERLAP_G4x_S+ABX*I_TWOBODYOVERLAP_F3x_S;
  Double I_TWOBODYOVERLAP_F2xy_Px = I_TWOBODYOVERLAP_G3xy_S+ABX*I_TWOBODYOVERLAP_F2xy_S;
  Double I_TWOBODYOVERLAP_F2xz_Px = I_TWOBODYOVERLAP_G3xz_S+ABX*I_TWOBODYOVERLAP_F2xz_S;
  Double I_TWOBODYOVERLAP_Fx2y_Px = I_TWOBODYOVERLAP_G2x2y_S+ABX*I_TWOBODYOVERLAP_Fx2y_S;
  Double I_TWOBODYOVERLAP_Fxyz_Px = I_TWOBODYOVERLAP_G2xyz_S+ABX*I_TWOBODYOVERLAP_Fxyz_S;
  Double I_TWOBODYOVERLAP_Fx2z_Px = I_TWOBODYOVERLAP_G2x2z_S+ABX*I_TWOBODYOVERLAP_Fx2z_S;
  Double I_TWOBODYOVERLAP_F3y_Px = I_TWOBODYOVERLAP_Gx3y_S+ABX*I_TWOBODYOVERLAP_F3y_S;
  Double I_TWOBODYOVERLAP_F2yz_Px = I_TWOBODYOVERLAP_Gx2yz_S+ABX*I_TWOBODYOVERLAP_F2yz_S;
  Double I_TWOBODYOVERLAP_Fy2z_Px = I_TWOBODYOVERLAP_Gxy2z_S+ABX*I_TWOBODYOVERLAP_Fy2z_S;
  Double I_TWOBODYOVERLAP_F3z_Px = I_TWOBODYOVERLAP_Gx3z_S+ABX*I_TWOBODYOVERLAP_F3z_S;
  Double I_TWOBODYOVERLAP_F2xy_Py = I_TWOBODYOVERLAP_G2x2y_S+ABY*I_TWOBODYOVERLAP_F2xy_S;
  Double I_TWOBODYOVERLAP_F2xz_Py = I_TWOBODYOVERLAP_G2xyz_S+ABY*I_TWOBODYOVERLAP_F2xz_S;
  Double I_TWOBODYOVERLAP_Fx2y_Py = I_TWOBODYOVERLAP_Gx3y_S+ABY*I_TWOBODYOVERLAP_Fx2y_S;
  Double I_TWOBODYOVERLAP_Fxyz_Py = I_TWOBODYOVERLAP_Gx2yz_S+ABY*I_TWOBODYOVERLAP_Fxyz_S;
  Double I_TWOBODYOVERLAP_Fx2z_Py = I_TWOBODYOVERLAP_Gxy2z_S+ABY*I_TWOBODYOVERLAP_Fx2z_S;
  Double I_TWOBODYOVERLAP_F3y_Py = I_TWOBODYOVERLAP_G4y_S+ABY*I_TWOBODYOVERLAP_F3y_S;
  Double I_TWOBODYOVERLAP_F2yz_Py = I_TWOBODYOVERLAP_G3yz_S+ABY*I_TWOBODYOVERLAP_F2yz_S;
  Double I_TWOBODYOVERLAP_Fy2z_Py = I_TWOBODYOVERLAP_G2y2z_S+ABY*I_TWOBODYOVERLAP_Fy2z_S;
  Double I_TWOBODYOVERLAP_F3z_Py = I_TWOBODYOVERLAP_Gy3z_S+ABY*I_TWOBODYOVERLAP_F3z_S;
  Double I_TWOBODYOVERLAP_F2xz_Pz = I_TWOBODYOVERLAP_G2x2z_S+ABZ*I_TWOBODYOVERLAP_F2xz_S;
  Double I_TWOBODYOVERLAP_Fxyz_Pz = I_TWOBODYOVERLAP_Gxy2z_S+ABZ*I_TWOBODYOVERLAP_Fxyz_S;
  Double I_TWOBODYOVERLAP_Fx2z_Pz = I_TWOBODYOVERLAP_Gx3z_S+ABZ*I_TWOBODYOVERLAP_Fx2z_S;
  Double I_TWOBODYOVERLAP_F2yz_Pz = I_TWOBODYOVERLAP_G2y2z_S+ABZ*I_TWOBODYOVERLAP_F2yz_S;
  Double I_TWOBODYOVERLAP_Fy2z_Pz = I_TWOBODYOVERLAP_Gy3z_S+ABZ*I_TWOBODYOVERLAP_Fy2z_S;
  Double I_TWOBODYOVERLAP_F3z_Pz = I_TWOBODYOVERLAP_G4z_S+ABZ*I_TWOBODYOVERLAP_F3z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_D_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
   ************************************************************/
  Double I_TWOBODYOVERLAP_D2x_D2x = I_TWOBODYOVERLAP_F3x_Px+ABX*I_TWOBODYOVERLAP_D2x_Px;
  Double I_TWOBODYOVERLAP_Dxy_D2x = I_TWOBODYOVERLAP_F2xy_Px+ABX*I_TWOBODYOVERLAP_Dxy_Px;
  Double I_TWOBODYOVERLAP_Dxz_D2x = I_TWOBODYOVERLAP_F2xz_Px+ABX*I_TWOBODYOVERLAP_Dxz_Px;
  Double I_TWOBODYOVERLAP_D2y_D2x = I_TWOBODYOVERLAP_Fx2y_Px+ABX*I_TWOBODYOVERLAP_D2y_Px;
  Double I_TWOBODYOVERLAP_Dyz_D2x = I_TWOBODYOVERLAP_Fxyz_Px+ABX*I_TWOBODYOVERLAP_Dyz_Px;
  Double I_TWOBODYOVERLAP_D2z_D2x = I_TWOBODYOVERLAP_Fx2z_Px+ABX*I_TWOBODYOVERLAP_D2z_Px;
  Double I_TWOBODYOVERLAP_D2x_Dxy = I_TWOBODYOVERLAP_F2xy_Px+ABY*I_TWOBODYOVERLAP_D2x_Px;
  Double I_TWOBODYOVERLAP_Dxy_Dxy = I_TWOBODYOVERLAP_Fx2y_Px+ABY*I_TWOBODYOVERLAP_Dxy_Px;
  Double I_TWOBODYOVERLAP_Dxz_Dxy = I_TWOBODYOVERLAP_Fxyz_Px+ABY*I_TWOBODYOVERLAP_Dxz_Px;
  Double I_TWOBODYOVERLAP_D2y_Dxy = I_TWOBODYOVERLAP_F3y_Px+ABY*I_TWOBODYOVERLAP_D2y_Px;
  Double I_TWOBODYOVERLAP_Dyz_Dxy = I_TWOBODYOVERLAP_F2yz_Px+ABY*I_TWOBODYOVERLAP_Dyz_Px;
  Double I_TWOBODYOVERLAP_D2z_Dxy = I_TWOBODYOVERLAP_Fy2z_Px+ABY*I_TWOBODYOVERLAP_D2z_Px;
  Double I_TWOBODYOVERLAP_D2x_Dxz = I_TWOBODYOVERLAP_F2xz_Px+ABZ*I_TWOBODYOVERLAP_D2x_Px;
  Double I_TWOBODYOVERLAP_Dxy_Dxz = I_TWOBODYOVERLAP_Fxyz_Px+ABZ*I_TWOBODYOVERLAP_Dxy_Px;
  Double I_TWOBODYOVERLAP_Dxz_Dxz = I_TWOBODYOVERLAP_Fx2z_Px+ABZ*I_TWOBODYOVERLAP_Dxz_Px;
  Double I_TWOBODYOVERLAP_D2y_Dxz = I_TWOBODYOVERLAP_F2yz_Px+ABZ*I_TWOBODYOVERLAP_D2y_Px;
  Double I_TWOBODYOVERLAP_Dyz_Dxz = I_TWOBODYOVERLAP_Fy2z_Px+ABZ*I_TWOBODYOVERLAP_Dyz_Px;
  Double I_TWOBODYOVERLAP_D2z_Dxz = I_TWOBODYOVERLAP_F3z_Px+ABZ*I_TWOBODYOVERLAP_D2z_Px;
  Double I_TWOBODYOVERLAP_D2x_D2y = I_TWOBODYOVERLAP_F2xy_Py+ABY*I_TWOBODYOVERLAP_D2x_Py;
  Double I_TWOBODYOVERLAP_Dxy_D2y = I_TWOBODYOVERLAP_Fx2y_Py+ABY*I_TWOBODYOVERLAP_Dxy_Py;
  Double I_TWOBODYOVERLAP_Dxz_D2y = I_TWOBODYOVERLAP_Fxyz_Py+ABY*I_TWOBODYOVERLAP_Dxz_Py;
  Double I_TWOBODYOVERLAP_D2y_D2y = I_TWOBODYOVERLAP_F3y_Py+ABY*I_TWOBODYOVERLAP_D2y_Py;
  Double I_TWOBODYOVERLAP_Dyz_D2y = I_TWOBODYOVERLAP_F2yz_Py+ABY*I_TWOBODYOVERLAP_Dyz_Py;
  Double I_TWOBODYOVERLAP_D2z_D2y = I_TWOBODYOVERLAP_Fy2z_Py+ABY*I_TWOBODYOVERLAP_D2z_Py;
  Double I_TWOBODYOVERLAP_D2x_Dyz = I_TWOBODYOVERLAP_F2xz_Py+ABZ*I_TWOBODYOVERLAP_D2x_Py;
  Double I_TWOBODYOVERLAP_Dxy_Dyz = I_TWOBODYOVERLAP_Fxyz_Py+ABZ*I_TWOBODYOVERLAP_Dxy_Py;
  Double I_TWOBODYOVERLAP_Dxz_Dyz = I_TWOBODYOVERLAP_Fx2z_Py+ABZ*I_TWOBODYOVERLAP_Dxz_Py;
  Double I_TWOBODYOVERLAP_D2y_Dyz = I_TWOBODYOVERLAP_F2yz_Py+ABZ*I_TWOBODYOVERLAP_D2y_Py;
  Double I_TWOBODYOVERLAP_Dyz_Dyz = I_TWOBODYOVERLAP_Fy2z_Py+ABZ*I_TWOBODYOVERLAP_Dyz_Py;
  Double I_TWOBODYOVERLAP_D2z_Dyz = I_TWOBODYOVERLAP_F3z_Py+ABZ*I_TWOBODYOVERLAP_D2z_Py;
  Double I_TWOBODYOVERLAP_D2x_D2z = I_TWOBODYOVERLAP_F2xz_Pz+ABZ*I_TWOBODYOVERLAP_D2x_Pz;
  Double I_TWOBODYOVERLAP_Dxy_D2z = I_TWOBODYOVERLAP_Fxyz_Pz+ABZ*I_TWOBODYOVERLAP_Dxy_Pz;
  Double I_TWOBODYOVERLAP_Dxz_D2z = I_TWOBODYOVERLAP_Fx2z_Pz+ABZ*I_TWOBODYOVERLAP_Dxz_Pz;
  Double I_TWOBODYOVERLAP_D2y_D2z = I_TWOBODYOVERLAP_F2yz_Pz+ABZ*I_TWOBODYOVERLAP_D2y_Pz;
  Double I_TWOBODYOVERLAP_Dyz_D2z = I_TWOBODYOVERLAP_Fy2z_Pz+ABZ*I_TWOBODYOVERLAP_Dyz_Pz;
  Double I_TWOBODYOVERLAP_D2z_D2z = I_TWOBODYOVERLAP_F3z_Pz+ABZ*I_TWOBODYOVERLAP_D2z_Pz;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_Px_a = I_TWOBODYOVERLAP_H5x_S_a+ABX*I_TWOBODYOVERLAP_G4x_S_a;
  Double I_TWOBODYOVERLAP_G3xy_Px_a = I_TWOBODYOVERLAP_H4xy_S_a+ABX*I_TWOBODYOVERLAP_G3xy_S_a;
  Double I_TWOBODYOVERLAP_G3xz_Px_a = I_TWOBODYOVERLAP_H4xz_S_a+ABX*I_TWOBODYOVERLAP_G3xz_S_a;
  Double I_TWOBODYOVERLAP_G2x2y_Px_a = I_TWOBODYOVERLAP_H3x2y_S_a+ABX*I_TWOBODYOVERLAP_G2x2y_S_a;
  Double I_TWOBODYOVERLAP_G2xyz_Px_a = I_TWOBODYOVERLAP_H3xyz_S_a+ABX*I_TWOBODYOVERLAP_G2xyz_S_a;
  Double I_TWOBODYOVERLAP_G2x2z_Px_a = I_TWOBODYOVERLAP_H3x2z_S_a+ABX*I_TWOBODYOVERLAP_G2x2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3y_Px_a = I_TWOBODYOVERLAP_H2x3y_S_a+ABX*I_TWOBODYOVERLAP_Gx3y_S_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Px_a = I_TWOBODYOVERLAP_H2x2yz_S_a+ABX*I_TWOBODYOVERLAP_Gx2yz_S_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Px_a = I_TWOBODYOVERLAP_H2xy2z_S_a+ABX*I_TWOBODYOVERLAP_Gxy2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3z_Px_a = I_TWOBODYOVERLAP_H2x3z_S_a+ABX*I_TWOBODYOVERLAP_Gx3z_S_a;
  Double I_TWOBODYOVERLAP_G4y_Px_a = I_TWOBODYOVERLAP_Hx4y_S_a+ABX*I_TWOBODYOVERLAP_G4y_S_a;
  Double I_TWOBODYOVERLAP_G3yz_Px_a = I_TWOBODYOVERLAP_Hx3yz_S_a+ABX*I_TWOBODYOVERLAP_G3yz_S_a;
  Double I_TWOBODYOVERLAP_G2y2z_Px_a = I_TWOBODYOVERLAP_Hx2y2z_S_a+ABX*I_TWOBODYOVERLAP_G2y2z_S_a;
  Double I_TWOBODYOVERLAP_Gy3z_Px_a = I_TWOBODYOVERLAP_Hxy3z_S_a+ABX*I_TWOBODYOVERLAP_Gy3z_S_a;
  Double I_TWOBODYOVERLAP_G4z_Px_a = I_TWOBODYOVERLAP_Hx4z_S_a+ABX*I_TWOBODYOVERLAP_G4z_S_a;
  Double I_TWOBODYOVERLAP_G4x_Py_a = I_TWOBODYOVERLAP_H4xy_S_a+ABY*I_TWOBODYOVERLAP_G4x_S_a;
  Double I_TWOBODYOVERLAP_G3xy_Py_a = I_TWOBODYOVERLAP_H3x2y_S_a+ABY*I_TWOBODYOVERLAP_G3xy_S_a;
  Double I_TWOBODYOVERLAP_G3xz_Py_a = I_TWOBODYOVERLAP_H3xyz_S_a+ABY*I_TWOBODYOVERLAP_G3xz_S_a;
  Double I_TWOBODYOVERLAP_G2x2y_Py_a = I_TWOBODYOVERLAP_H2x3y_S_a+ABY*I_TWOBODYOVERLAP_G2x2y_S_a;
  Double I_TWOBODYOVERLAP_G2xyz_Py_a = I_TWOBODYOVERLAP_H2x2yz_S_a+ABY*I_TWOBODYOVERLAP_G2xyz_S_a;
  Double I_TWOBODYOVERLAP_G2x2z_Py_a = I_TWOBODYOVERLAP_H2xy2z_S_a+ABY*I_TWOBODYOVERLAP_G2x2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3y_Py_a = I_TWOBODYOVERLAP_Hx4y_S_a+ABY*I_TWOBODYOVERLAP_Gx3y_S_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Py_a = I_TWOBODYOVERLAP_Hx3yz_S_a+ABY*I_TWOBODYOVERLAP_Gx2yz_S_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Py_a = I_TWOBODYOVERLAP_Hx2y2z_S_a+ABY*I_TWOBODYOVERLAP_Gxy2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3z_Py_a = I_TWOBODYOVERLAP_Hxy3z_S_a+ABY*I_TWOBODYOVERLAP_Gx3z_S_a;
  Double I_TWOBODYOVERLAP_G4y_Py_a = I_TWOBODYOVERLAP_H5y_S_a+ABY*I_TWOBODYOVERLAP_G4y_S_a;
  Double I_TWOBODYOVERLAP_G3yz_Py_a = I_TWOBODYOVERLAP_H4yz_S_a+ABY*I_TWOBODYOVERLAP_G3yz_S_a;
  Double I_TWOBODYOVERLAP_G2y2z_Py_a = I_TWOBODYOVERLAP_H3y2z_S_a+ABY*I_TWOBODYOVERLAP_G2y2z_S_a;
  Double I_TWOBODYOVERLAP_Gy3z_Py_a = I_TWOBODYOVERLAP_H2y3z_S_a+ABY*I_TWOBODYOVERLAP_Gy3z_S_a;
  Double I_TWOBODYOVERLAP_G4z_Py_a = I_TWOBODYOVERLAP_Hy4z_S_a+ABY*I_TWOBODYOVERLAP_G4z_S_a;
  Double I_TWOBODYOVERLAP_G4x_Pz_a = I_TWOBODYOVERLAP_H4xz_S_a+ABZ*I_TWOBODYOVERLAP_G4x_S_a;
  Double I_TWOBODYOVERLAP_G3xy_Pz_a = I_TWOBODYOVERLAP_H3xyz_S_a+ABZ*I_TWOBODYOVERLAP_G3xy_S_a;
  Double I_TWOBODYOVERLAP_G3xz_Pz_a = I_TWOBODYOVERLAP_H3x2z_S_a+ABZ*I_TWOBODYOVERLAP_G3xz_S_a;
  Double I_TWOBODYOVERLAP_G2x2y_Pz_a = I_TWOBODYOVERLAP_H2x2yz_S_a+ABZ*I_TWOBODYOVERLAP_G2x2y_S_a;
  Double I_TWOBODYOVERLAP_G2xyz_Pz_a = I_TWOBODYOVERLAP_H2xy2z_S_a+ABZ*I_TWOBODYOVERLAP_G2xyz_S_a;
  Double I_TWOBODYOVERLAP_G2x2z_Pz_a = I_TWOBODYOVERLAP_H2x3z_S_a+ABZ*I_TWOBODYOVERLAP_G2x2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3y_Pz_a = I_TWOBODYOVERLAP_Hx3yz_S_a+ABZ*I_TWOBODYOVERLAP_Gx3y_S_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Pz_a = I_TWOBODYOVERLAP_Hx2y2z_S_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_S_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Pz_a = I_TWOBODYOVERLAP_Hxy3z_S_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3z_Pz_a = I_TWOBODYOVERLAP_Hx4z_S_a+ABZ*I_TWOBODYOVERLAP_Gx3z_S_a;
  Double I_TWOBODYOVERLAP_G4y_Pz_a = I_TWOBODYOVERLAP_H4yz_S_a+ABZ*I_TWOBODYOVERLAP_G4y_S_a;
  Double I_TWOBODYOVERLAP_G3yz_Pz_a = I_TWOBODYOVERLAP_H3y2z_S_a+ABZ*I_TWOBODYOVERLAP_G3yz_S_a;
  Double I_TWOBODYOVERLAP_G2y2z_Pz_a = I_TWOBODYOVERLAP_H2y3z_S_a+ABZ*I_TWOBODYOVERLAP_G2y2z_S_a;
  Double I_TWOBODYOVERLAP_Gy3z_Pz_a = I_TWOBODYOVERLAP_Hy4z_S_a+ABZ*I_TWOBODYOVERLAP_Gy3z_S_a;
  Double I_TWOBODYOVERLAP_G4z_Pz_a = I_TWOBODYOVERLAP_H5z_S_a+ABZ*I_TWOBODYOVERLAP_G4z_S_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 7 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_Px_a = I_TWOBODYOVERLAP_I6x_S_a+ABX*I_TWOBODYOVERLAP_H5x_S_a;
  Double I_TWOBODYOVERLAP_H4xy_Px_a = I_TWOBODYOVERLAP_I5xy_S_a+ABX*I_TWOBODYOVERLAP_H4xy_S_a;
  Double I_TWOBODYOVERLAP_H4xz_Px_a = I_TWOBODYOVERLAP_I5xz_S_a+ABX*I_TWOBODYOVERLAP_H4xz_S_a;
  Double I_TWOBODYOVERLAP_H3x2y_Px_a = I_TWOBODYOVERLAP_I4x2y_S_a+ABX*I_TWOBODYOVERLAP_H3x2y_S_a;
  Double I_TWOBODYOVERLAP_H3xyz_Px_a = I_TWOBODYOVERLAP_I4xyz_S_a+ABX*I_TWOBODYOVERLAP_H3xyz_S_a;
  Double I_TWOBODYOVERLAP_H3x2z_Px_a = I_TWOBODYOVERLAP_I4x2z_S_a+ABX*I_TWOBODYOVERLAP_H3x2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3y_Px_a = I_TWOBODYOVERLAP_I3x3y_S_a+ABX*I_TWOBODYOVERLAP_H2x3y_S_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Px_a = I_TWOBODYOVERLAP_I3x2yz_S_a+ABX*I_TWOBODYOVERLAP_H2x2yz_S_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Px_a = I_TWOBODYOVERLAP_I3xy2z_S_a+ABX*I_TWOBODYOVERLAP_H2xy2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3z_Px_a = I_TWOBODYOVERLAP_I3x3z_S_a+ABX*I_TWOBODYOVERLAP_H2x3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4y_Px_a = I_TWOBODYOVERLAP_I2x4y_S_a+ABX*I_TWOBODYOVERLAP_Hx4y_S_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Px_a = I_TWOBODYOVERLAP_I2x3yz_S_a+ABX*I_TWOBODYOVERLAP_Hx3yz_S_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Px_a = I_TWOBODYOVERLAP_I2x2y2z_S_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_S_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Px_a = I_TWOBODYOVERLAP_I2xy3z_S_a+ABX*I_TWOBODYOVERLAP_Hxy3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4z_Px_a = I_TWOBODYOVERLAP_I2x4z_S_a+ABX*I_TWOBODYOVERLAP_Hx4z_S_a;
  Double I_TWOBODYOVERLAP_H5y_Px_a = I_TWOBODYOVERLAP_Ix5y_S_a+ABX*I_TWOBODYOVERLAP_H5y_S_a;
  Double I_TWOBODYOVERLAP_H4yz_Px_a = I_TWOBODYOVERLAP_Ix4yz_S_a+ABX*I_TWOBODYOVERLAP_H4yz_S_a;
  Double I_TWOBODYOVERLAP_H3y2z_Px_a = I_TWOBODYOVERLAP_Ix3y2z_S_a+ABX*I_TWOBODYOVERLAP_H3y2z_S_a;
  Double I_TWOBODYOVERLAP_H2y3z_Px_a = I_TWOBODYOVERLAP_Ix2y3z_S_a+ABX*I_TWOBODYOVERLAP_H2y3z_S_a;
  Double I_TWOBODYOVERLAP_Hy4z_Px_a = I_TWOBODYOVERLAP_Ixy4z_S_a+ABX*I_TWOBODYOVERLAP_Hy4z_S_a;
  Double I_TWOBODYOVERLAP_H5z_Px_a = I_TWOBODYOVERLAP_Ix5z_S_a+ABX*I_TWOBODYOVERLAP_H5z_S_a;
  Double I_TWOBODYOVERLAP_H4xy_Py_a = I_TWOBODYOVERLAP_I4x2y_S_a+ABY*I_TWOBODYOVERLAP_H4xy_S_a;
  Double I_TWOBODYOVERLAP_H4xz_Py_a = I_TWOBODYOVERLAP_I4xyz_S_a+ABY*I_TWOBODYOVERLAP_H4xz_S_a;
  Double I_TWOBODYOVERLAP_H3x2y_Py_a = I_TWOBODYOVERLAP_I3x3y_S_a+ABY*I_TWOBODYOVERLAP_H3x2y_S_a;
  Double I_TWOBODYOVERLAP_H3xyz_Py_a = I_TWOBODYOVERLAP_I3x2yz_S_a+ABY*I_TWOBODYOVERLAP_H3xyz_S_a;
  Double I_TWOBODYOVERLAP_H3x2z_Py_a = I_TWOBODYOVERLAP_I3xy2z_S_a+ABY*I_TWOBODYOVERLAP_H3x2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3y_Py_a = I_TWOBODYOVERLAP_I2x4y_S_a+ABY*I_TWOBODYOVERLAP_H2x3y_S_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Py_a = I_TWOBODYOVERLAP_I2x3yz_S_a+ABY*I_TWOBODYOVERLAP_H2x2yz_S_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Py_a = I_TWOBODYOVERLAP_I2x2y2z_S_a+ABY*I_TWOBODYOVERLAP_H2xy2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3z_Py_a = I_TWOBODYOVERLAP_I2xy3z_S_a+ABY*I_TWOBODYOVERLAP_H2x3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4y_Py_a = I_TWOBODYOVERLAP_Ix5y_S_a+ABY*I_TWOBODYOVERLAP_Hx4y_S_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Py_a = I_TWOBODYOVERLAP_Ix4yz_S_a+ABY*I_TWOBODYOVERLAP_Hx3yz_S_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Py_a = I_TWOBODYOVERLAP_Ix3y2z_S_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_S_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Py_a = I_TWOBODYOVERLAP_Ix2y3z_S_a+ABY*I_TWOBODYOVERLAP_Hxy3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4z_Py_a = I_TWOBODYOVERLAP_Ixy4z_S_a+ABY*I_TWOBODYOVERLAP_Hx4z_S_a;
  Double I_TWOBODYOVERLAP_H5y_Py_a = I_TWOBODYOVERLAP_I6y_S_a+ABY*I_TWOBODYOVERLAP_H5y_S_a;
  Double I_TWOBODYOVERLAP_H4yz_Py_a = I_TWOBODYOVERLAP_I5yz_S_a+ABY*I_TWOBODYOVERLAP_H4yz_S_a;
  Double I_TWOBODYOVERLAP_H3y2z_Py_a = I_TWOBODYOVERLAP_I4y2z_S_a+ABY*I_TWOBODYOVERLAP_H3y2z_S_a;
  Double I_TWOBODYOVERLAP_H2y3z_Py_a = I_TWOBODYOVERLAP_I3y3z_S_a+ABY*I_TWOBODYOVERLAP_H2y3z_S_a;
  Double I_TWOBODYOVERLAP_Hy4z_Py_a = I_TWOBODYOVERLAP_I2y4z_S_a+ABY*I_TWOBODYOVERLAP_Hy4z_S_a;
  Double I_TWOBODYOVERLAP_H5z_Py_a = I_TWOBODYOVERLAP_Iy5z_S_a+ABY*I_TWOBODYOVERLAP_H5z_S_a;
  Double I_TWOBODYOVERLAP_H4xz_Pz_a = I_TWOBODYOVERLAP_I4x2z_S_a+ABZ*I_TWOBODYOVERLAP_H4xz_S_a;
  Double I_TWOBODYOVERLAP_H3xyz_Pz_a = I_TWOBODYOVERLAP_I3xy2z_S_a+ABZ*I_TWOBODYOVERLAP_H3xyz_S_a;
  Double I_TWOBODYOVERLAP_H3x2z_Pz_a = I_TWOBODYOVERLAP_I3x3z_S_a+ABZ*I_TWOBODYOVERLAP_H3x2z_S_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Pz_a = I_TWOBODYOVERLAP_I2x2y2z_S_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_S_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Pz_a = I_TWOBODYOVERLAP_I2xy3z_S_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3z_Pz_a = I_TWOBODYOVERLAP_I2x4z_S_a+ABZ*I_TWOBODYOVERLAP_H2x3z_S_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Pz_a = I_TWOBODYOVERLAP_Ix3y2z_S_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_S_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Pz_a = I_TWOBODYOVERLAP_Ix2y3z_S_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_S_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Pz_a = I_TWOBODYOVERLAP_Ixy4z_S_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4z_Pz_a = I_TWOBODYOVERLAP_Ix5z_S_a+ABZ*I_TWOBODYOVERLAP_Hx4z_S_a;
  Double I_TWOBODYOVERLAP_H4yz_Pz_a = I_TWOBODYOVERLAP_I4y2z_S_a+ABZ*I_TWOBODYOVERLAP_H4yz_S_a;
  Double I_TWOBODYOVERLAP_H3y2z_Pz_a = I_TWOBODYOVERLAP_I3y3z_S_a+ABZ*I_TWOBODYOVERLAP_H3y2z_S_a;
  Double I_TWOBODYOVERLAP_H2y3z_Pz_a = I_TWOBODYOVERLAP_I2y4z_S_a+ABZ*I_TWOBODYOVERLAP_H2y3z_S_a;
  Double I_TWOBODYOVERLAP_Hy4z_Pz_a = I_TWOBODYOVERLAP_Iy5z_S_a+ABZ*I_TWOBODYOVERLAP_Hy4z_S_a;
  Double I_TWOBODYOVERLAP_H5z_Pz_a = I_TWOBODYOVERLAP_I6z_S_a+ABZ*I_TWOBODYOVERLAP_H5z_S_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_D2x_a = I_TWOBODYOVERLAP_H5x_Px_a+ABX*I_TWOBODYOVERLAP_G4x_Px_a;
  Double I_TWOBODYOVERLAP_G3xy_D2x_a = I_TWOBODYOVERLAP_H4xy_Px_a+ABX*I_TWOBODYOVERLAP_G3xy_Px_a;
  Double I_TWOBODYOVERLAP_G3xz_D2x_a = I_TWOBODYOVERLAP_H4xz_Px_a+ABX*I_TWOBODYOVERLAP_G3xz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2y_D2x_a = I_TWOBODYOVERLAP_H3x2y_Px_a+ABX*I_TWOBODYOVERLAP_G2x2y_Px_a;
  Double I_TWOBODYOVERLAP_G2xyz_D2x_a = I_TWOBODYOVERLAP_H3xyz_Px_a+ABX*I_TWOBODYOVERLAP_G2xyz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2z_D2x_a = I_TWOBODYOVERLAP_H3x2z_Px_a+ABX*I_TWOBODYOVERLAP_G2x2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3y_D2x_a = I_TWOBODYOVERLAP_H2x3y_Px_a+ABX*I_TWOBODYOVERLAP_Gx3y_Px_a;
  Double I_TWOBODYOVERLAP_Gx2yz_D2x_a = I_TWOBODYOVERLAP_H2x2yz_Px_a+ABX*I_TWOBODYOVERLAP_Gx2yz_Px_a;
  Double I_TWOBODYOVERLAP_Gxy2z_D2x_a = I_TWOBODYOVERLAP_H2xy2z_Px_a+ABX*I_TWOBODYOVERLAP_Gxy2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3z_D2x_a = I_TWOBODYOVERLAP_H2x3z_Px_a+ABX*I_TWOBODYOVERLAP_Gx3z_Px_a;
  Double I_TWOBODYOVERLAP_G4y_D2x_a = I_TWOBODYOVERLAP_Hx4y_Px_a+ABX*I_TWOBODYOVERLAP_G4y_Px_a;
  Double I_TWOBODYOVERLAP_G3yz_D2x_a = I_TWOBODYOVERLAP_Hx3yz_Px_a+ABX*I_TWOBODYOVERLAP_G3yz_Px_a;
  Double I_TWOBODYOVERLAP_G2y2z_D2x_a = I_TWOBODYOVERLAP_Hx2y2z_Px_a+ABX*I_TWOBODYOVERLAP_G2y2z_Px_a;
  Double I_TWOBODYOVERLAP_Gy3z_D2x_a = I_TWOBODYOVERLAP_Hxy3z_Px_a+ABX*I_TWOBODYOVERLAP_Gy3z_Px_a;
  Double I_TWOBODYOVERLAP_G4z_D2x_a = I_TWOBODYOVERLAP_Hx4z_Px_a+ABX*I_TWOBODYOVERLAP_G4z_Px_a;
  Double I_TWOBODYOVERLAP_G4x_Dxy_a = I_TWOBODYOVERLAP_H4xy_Px_a+ABY*I_TWOBODYOVERLAP_G4x_Px_a;
  Double I_TWOBODYOVERLAP_G3xy_Dxy_a = I_TWOBODYOVERLAP_H3x2y_Px_a+ABY*I_TWOBODYOVERLAP_G3xy_Px_a;
  Double I_TWOBODYOVERLAP_G3xz_Dxy_a = I_TWOBODYOVERLAP_H3xyz_Px_a+ABY*I_TWOBODYOVERLAP_G3xz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2y_Dxy_a = I_TWOBODYOVERLAP_H2x3y_Px_a+ABY*I_TWOBODYOVERLAP_G2x2y_Px_a;
  Double I_TWOBODYOVERLAP_G2xyz_Dxy_a = I_TWOBODYOVERLAP_H2x2yz_Px_a+ABY*I_TWOBODYOVERLAP_G2xyz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2z_Dxy_a = I_TWOBODYOVERLAP_H2xy2z_Px_a+ABY*I_TWOBODYOVERLAP_G2x2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3y_Dxy_a = I_TWOBODYOVERLAP_Hx4y_Px_a+ABY*I_TWOBODYOVERLAP_Gx3y_Px_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Dxy_a = I_TWOBODYOVERLAP_Hx3yz_Px_a+ABY*I_TWOBODYOVERLAP_Gx2yz_Px_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Dxy_a = I_TWOBODYOVERLAP_Hx2y2z_Px_a+ABY*I_TWOBODYOVERLAP_Gxy2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3z_Dxy_a = I_TWOBODYOVERLAP_Hxy3z_Px_a+ABY*I_TWOBODYOVERLAP_Gx3z_Px_a;
  Double I_TWOBODYOVERLAP_G4y_Dxy_a = I_TWOBODYOVERLAP_H5y_Px_a+ABY*I_TWOBODYOVERLAP_G4y_Px_a;
  Double I_TWOBODYOVERLAP_G3yz_Dxy_a = I_TWOBODYOVERLAP_H4yz_Px_a+ABY*I_TWOBODYOVERLAP_G3yz_Px_a;
  Double I_TWOBODYOVERLAP_G2y2z_Dxy_a = I_TWOBODYOVERLAP_H3y2z_Px_a+ABY*I_TWOBODYOVERLAP_G2y2z_Px_a;
  Double I_TWOBODYOVERLAP_Gy3z_Dxy_a = I_TWOBODYOVERLAP_H2y3z_Px_a+ABY*I_TWOBODYOVERLAP_Gy3z_Px_a;
  Double I_TWOBODYOVERLAP_G4z_Dxy_a = I_TWOBODYOVERLAP_Hy4z_Px_a+ABY*I_TWOBODYOVERLAP_G4z_Px_a;
  Double I_TWOBODYOVERLAP_G4x_Dxz_a = I_TWOBODYOVERLAP_H4xz_Px_a+ABZ*I_TWOBODYOVERLAP_G4x_Px_a;
  Double I_TWOBODYOVERLAP_G3xy_Dxz_a = I_TWOBODYOVERLAP_H3xyz_Px_a+ABZ*I_TWOBODYOVERLAP_G3xy_Px_a;
  Double I_TWOBODYOVERLAP_G3xz_Dxz_a = I_TWOBODYOVERLAP_H3x2z_Px_a+ABZ*I_TWOBODYOVERLAP_G3xz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2y_Dxz_a = I_TWOBODYOVERLAP_H2x2yz_Px_a+ABZ*I_TWOBODYOVERLAP_G2x2y_Px_a;
  Double I_TWOBODYOVERLAP_G2xyz_Dxz_a = I_TWOBODYOVERLAP_H2xy2z_Px_a+ABZ*I_TWOBODYOVERLAP_G2xyz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2z_Dxz_a = I_TWOBODYOVERLAP_H2x3z_Px_a+ABZ*I_TWOBODYOVERLAP_G2x2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3y_Dxz_a = I_TWOBODYOVERLAP_Hx3yz_Px_a+ABZ*I_TWOBODYOVERLAP_Gx3y_Px_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Dxz_a = I_TWOBODYOVERLAP_Hx2y2z_Px_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_Px_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Dxz_a = I_TWOBODYOVERLAP_Hxy3z_Px_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3z_Dxz_a = I_TWOBODYOVERLAP_Hx4z_Px_a+ABZ*I_TWOBODYOVERLAP_Gx3z_Px_a;
  Double I_TWOBODYOVERLAP_G4y_Dxz_a = I_TWOBODYOVERLAP_H4yz_Px_a+ABZ*I_TWOBODYOVERLAP_G4y_Px_a;
  Double I_TWOBODYOVERLAP_G3yz_Dxz_a = I_TWOBODYOVERLAP_H3y2z_Px_a+ABZ*I_TWOBODYOVERLAP_G3yz_Px_a;
  Double I_TWOBODYOVERLAP_G2y2z_Dxz_a = I_TWOBODYOVERLAP_H2y3z_Px_a+ABZ*I_TWOBODYOVERLAP_G2y2z_Px_a;
  Double I_TWOBODYOVERLAP_Gy3z_Dxz_a = I_TWOBODYOVERLAP_Hy4z_Px_a+ABZ*I_TWOBODYOVERLAP_Gy3z_Px_a;
  Double I_TWOBODYOVERLAP_G4z_Dxz_a = I_TWOBODYOVERLAP_H5z_Px_a+ABZ*I_TWOBODYOVERLAP_G4z_Px_a;
  Double I_TWOBODYOVERLAP_G4x_D2y_a = I_TWOBODYOVERLAP_H4xy_Py_a+ABY*I_TWOBODYOVERLAP_G4x_Py_a;
  Double I_TWOBODYOVERLAP_G3xy_D2y_a = I_TWOBODYOVERLAP_H3x2y_Py_a+ABY*I_TWOBODYOVERLAP_G3xy_Py_a;
  Double I_TWOBODYOVERLAP_G3xz_D2y_a = I_TWOBODYOVERLAP_H3xyz_Py_a+ABY*I_TWOBODYOVERLAP_G3xz_Py_a;
  Double I_TWOBODYOVERLAP_G2x2y_D2y_a = I_TWOBODYOVERLAP_H2x3y_Py_a+ABY*I_TWOBODYOVERLAP_G2x2y_Py_a;
  Double I_TWOBODYOVERLAP_G2xyz_D2y_a = I_TWOBODYOVERLAP_H2x2yz_Py_a+ABY*I_TWOBODYOVERLAP_G2xyz_Py_a;
  Double I_TWOBODYOVERLAP_G2x2z_D2y_a = I_TWOBODYOVERLAP_H2xy2z_Py_a+ABY*I_TWOBODYOVERLAP_G2x2z_Py_a;
  Double I_TWOBODYOVERLAP_Gx3y_D2y_a = I_TWOBODYOVERLAP_Hx4y_Py_a+ABY*I_TWOBODYOVERLAP_Gx3y_Py_a;
  Double I_TWOBODYOVERLAP_Gx2yz_D2y_a = I_TWOBODYOVERLAP_Hx3yz_Py_a+ABY*I_TWOBODYOVERLAP_Gx2yz_Py_a;
  Double I_TWOBODYOVERLAP_Gxy2z_D2y_a = I_TWOBODYOVERLAP_Hx2y2z_Py_a+ABY*I_TWOBODYOVERLAP_Gxy2z_Py_a;
  Double I_TWOBODYOVERLAP_Gx3z_D2y_a = I_TWOBODYOVERLAP_Hxy3z_Py_a+ABY*I_TWOBODYOVERLAP_Gx3z_Py_a;
  Double I_TWOBODYOVERLAP_G4y_D2y_a = I_TWOBODYOVERLAP_H5y_Py_a+ABY*I_TWOBODYOVERLAP_G4y_Py_a;
  Double I_TWOBODYOVERLAP_G3yz_D2y_a = I_TWOBODYOVERLAP_H4yz_Py_a+ABY*I_TWOBODYOVERLAP_G3yz_Py_a;
  Double I_TWOBODYOVERLAP_G2y2z_D2y_a = I_TWOBODYOVERLAP_H3y2z_Py_a+ABY*I_TWOBODYOVERLAP_G2y2z_Py_a;
  Double I_TWOBODYOVERLAP_Gy3z_D2y_a = I_TWOBODYOVERLAP_H2y3z_Py_a+ABY*I_TWOBODYOVERLAP_Gy3z_Py_a;
  Double I_TWOBODYOVERLAP_G4z_D2y_a = I_TWOBODYOVERLAP_Hy4z_Py_a+ABY*I_TWOBODYOVERLAP_G4z_Py_a;
  Double I_TWOBODYOVERLAP_G4x_Dyz_a = I_TWOBODYOVERLAP_H4xz_Py_a+ABZ*I_TWOBODYOVERLAP_G4x_Py_a;
  Double I_TWOBODYOVERLAP_G3xy_Dyz_a = I_TWOBODYOVERLAP_H3xyz_Py_a+ABZ*I_TWOBODYOVERLAP_G3xy_Py_a;
  Double I_TWOBODYOVERLAP_G3xz_Dyz_a = I_TWOBODYOVERLAP_H3x2z_Py_a+ABZ*I_TWOBODYOVERLAP_G3xz_Py_a;
  Double I_TWOBODYOVERLAP_G2x2y_Dyz_a = I_TWOBODYOVERLAP_H2x2yz_Py_a+ABZ*I_TWOBODYOVERLAP_G2x2y_Py_a;
  Double I_TWOBODYOVERLAP_G2xyz_Dyz_a = I_TWOBODYOVERLAP_H2xy2z_Py_a+ABZ*I_TWOBODYOVERLAP_G2xyz_Py_a;
  Double I_TWOBODYOVERLAP_G2x2z_Dyz_a = I_TWOBODYOVERLAP_H2x3z_Py_a+ABZ*I_TWOBODYOVERLAP_G2x2z_Py_a;
  Double I_TWOBODYOVERLAP_Gx3y_Dyz_a = I_TWOBODYOVERLAP_Hx3yz_Py_a+ABZ*I_TWOBODYOVERLAP_Gx3y_Py_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Dyz_a = I_TWOBODYOVERLAP_Hx2y2z_Py_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_Py_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Dyz_a = I_TWOBODYOVERLAP_Hxy3z_Py_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_Py_a;
  Double I_TWOBODYOVERLAP_Gx3z_Dyz_a = I_TWOBODYOVERLAP_Hx4z_Py_a+ABZ*I_TWOBODYOVERLAP_Gx3z_Py_a;
  Double I_TWOBODYOVERLAP_G4y_Dyz_a = I_TWOBODYOVERLAP_H4yz_Py_a+ABZ*I_TWOBODYOVERLAP_G4y_Py_a;
  Double I_TWOBODYOVERLAP_G3yz_Dyz_a = I_TWOBODYOVERLAP_H3y2z_Py_a+ABZ*I_TWOBODYOVERLAP_G3yz_Py_a;
  Double I_TWOBODYOVERLAP_G2y2z_Dyz_a = I_TWOBODYOVERLAP_H2y3z_Py_a+ABZ*I_TWOBODYOVERLAP_G2y2z_Py_a;
  Double I_TWOBODYOVERLAP_Gy3z_Dyz_a = I_TWOBODYOVERLAP_Hy4z_Py_a+ABZ*I_TWOBODYOVERLAP_Gy3z_Py_a;
  Double I_TWOBODYOVERLAP_G4z_Dyz_a = I_TWOBODYOVERLAP_H5z_Py_a+ABZ*I_TWOBODYOVERLAP_G4z_Py_a;
  Double I_TWOBODYOVERLAP_G4x_D2z_a = I_TWOBODYOVERLAP_H4xz_Pz_a+ABZ*I_TWOBODYOVERLAP_G4x_Pz_a;
  Double I_TWOBODYOVERLAP_G3xy_D2z_a = I_TWOBODYOVERLAP_H3xyz_Pz_a+ABZ*I_TWOBODYOVERLAP_G3xy_Pz_a;
  Double I_TWOBODYOVERLAP_G3xz_D2z_a = I_TWOBODYOVERLAP_H3x2z_Pz_a+ABZ*I_TWOBODYOVERLAP_G3xz_Pz_a;
  Double I_TWOBODYOVERLAP_G2x2y_D2z_a = I_TWOBODYOVERLAP_H2x2yz_Pz_a+ABZ*I_TWOBODYOVERLAP_G2x2y_Pz_a;
  Double I_TWOBODYOVERLAP_G2xyz_D2z_a = I_TWOBODYOVERLAP_H2xy2z_Pz_a+ABZ*I_TWOBODYOVERLAP_G2xyz_Pz_a;
  Double I_TWOBODYOVERLAP_G2x2z_D2z_a = I_TWOBODYOVERLAP_H2x3z_Pz_a+ABZ*I_TWOBODYOVERLAP_G2x2z_Pz_a;
  Double I_TWOBODYOVERLAP_Gx3y_D2z_a = I_TWOBODYOVERLAP_Hx3yz_Pz_a+ABZ*I_TWOBODYOVERLAP_Gx3y_Pz_a;
  Double I_TWOBODYOVERLAP_Gx2yz_D2z_a = I_TWOBODYOVERLAP_Hx2y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_Pz_a;
  Double I_TWOBODYOVERLAP_Gxy2z_D2z_a = I_TWOBODYOVERLAP_Hxy3z_Pz_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_Pz_a;
  Double I_TWOBODYOVERLAP_Gx3z_D2z_a = I_TWOBODYOVERLAP_Hx4z_Pz_a+ABZ*I_TWOBODYOVERLAP_Gx3z_Pz_a;
  Double I_TWOBODYOVERLAP_G4y_D2z_a = I_TWOBODYOVERLAP_H4yz_Pz_a+ABZ*I_TWOBODYOVERLAP_G4y_Pz_a;
  Double I_TWOBODYOVERLAP_G3yz_D2z_a = I_TWOBODYOVERLAP_H3y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_G3yz_Pz_a;
  Double I_TWOBODYOVERLAP_G2y2z_D2z_a = I_TWOBODYOVERLAP_H2y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_G2y2z_Pz_a;
  Double I_TWOBODYOVERLAP_Gy3z_D2z_a = I_TWOBODYOVERLAP_Hy4z_Pz_a+ABZ*I_TWOBODYOVERLAP_Gy3z_Pz_a;
  Double I_TWOBODYOVERLAP_G4z_D2z_a = I_TWOBODYOVERLAP_H5z_Pz_a+ABZ*I_TWOBODYOVERLAP_G4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_Px_aa = I_TWOBODYOVERLAP_K7x_S_aa+ABX*I_TWOBODYOVERLAP_I6x_S_aa;
  Double I_TWOBODYOVERLAP_I5xy_Px_aa = I_TWOBODYOVERLAP_K6xy_S_aa+ABX*I_TWOBODYOVERLAP_I5xy_S_aa;
  Double I_TWOBODYOVERLAP_I5xz_Px_aa = I_TWOBODYOVERLAP_K6xz_S_aa+ABX*I_TWOBODYOVERLAP_I5xz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Px_aa = I_TWOBODYOVERLAP_K5x2y_S_aa+ABX*I_TWOBODYOVERLAP_I4x2y_S_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Px_aa = I_TWOBODYOVERLAP_K5xyz_S_aa+ABX*I_TWOBODYOVERLAP_I4xyz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Px_aa = I_TWOBODYOVERLAP_K5x2z_S_aa+ABX*I_TWOBODYOVERLAP_I4x2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Px_aa = I_TWOBODYOVERLAP_K4x3y_S_aa+ABX*I_TWOBODYOVERLAP_I3x3y_S_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Px_aa = I_TWOBODYOVERLAP_K4x2yz_S_aa+ABX*I_TWOBODYOVERLAP_I3x2yz_S_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Px_aa = I_TWOBODYOVERLAP_K4xy2z_S_aa+ABX*I_TWOBODYOVERLAP_I3xy2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Px_aa = I_TWOBODYOVERLAP_K4x3z_S_aa+ABX*I_TWOBODYOVERLAP_I3x3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Px_aa = I_TWOBODYOVERLAP_K3x4y_S_aa+ABX*I_TWOBODYOVERLAP_I2x4y_S_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Px_aa = I_TWOBODYOVERLAP_K3x3yz_S_aa+ABX*I_TWOBODYOVERLAP_I2x3yz_S_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Px_aa = I_TWOBODYOVERLAP_K3x2y2z_S_aa+ABX*I_TWOBODYOVERLAP_I2x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Px_aa = I_TWOBODYOVERLAP_K3xy3z_S_aa+ABX*I_TWOBODYOVERLAP_I2xy3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Px_aa = I_TWOBODYOVERLAP_K3x4z_S_aa+ABX*I_TWOBODYOVERLAP_I2x4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Px_aa = I_TWOBODYOVERLAP_K2x5y_S_aa+ABX*I_TWOBODYOVERLAP_Ix5y_S_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Px_aa = I_TWOBODYOVERLAP_K2x4yz_S_aa+ABX*I_TWOBODYOVERLAP_Ix4yz_S_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Px_aa = I_TWOBODYOVERLAP_K2x3y2z_S_aa+ABX*I_TWOBODYOVERLAP_Ix3y2z_S_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Px_aa = I_TWOBODYOVERLAP_K2x2y3z_S_aa+ABX*I_TWOBODYOVERLAP_Ix2y3z_S_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Px_aa = I_TWOBODYOVERLAP_K2xy4z_S_aa+ABX*I_TWOBODYOVERLAP_Ixy4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Px_aa = I_TWOBODYOVERLAP_K2x5z_S_aa+ABX*I_TWOBODYOVERLAP_Ix5z_S_aa;
  Double I_TWOBODYOVERLAP_I6y_Px_aa = I_TWOBODYOVERLAP_Kx6y_S_aa+ABX*I_TWOBODYOVERLAP_I6y_S_aa;
  Double I_TWOBODYOVERLAP_I5yz_Px_aa = I_TWOBODYOVERLAP_Kx5yz_S_aa+ABX*I_TWOBODYOVERLAP_I5yz_S_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Px_aa = I_TWOBODYOVERLAP_Kx4y2z_S_aa+ABX*I_TWOBODYOVERLAP_I4y2z_S_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Px_aa = I_TWOBODYOVERLAP_Kx3y3z_S_aa+ABX*I_TWOBODYOVERLAP_I3y3z_S_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Px_aa = I_TWOBODYOVERLAP_Kx2y4z_S_aa+ABX*I_TWOBODYOVERLAP_I2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Px_aa = I_TWOBODYOVERLAP_Kxy5z_S_aa+ABX*I_TWOBODYOVERLAP_Iy5z_S_aa;
  Double I_TWOBODYOVERLAP_I6z_Px_aa = I_TWOBODYOVERLAP_Kx6z_S_aa+ABX*I_TWOBODYOVERLAP_I6z_S_aa;
  Double I_TWOBODYOVERLAP_I6x_Py_aa = I_TWOBODYOVERLAP_K6xy_S_aa+ABY*I_TWOBODYOVERLAP_I6x_S_aa;
  Double I_TWOBODYOVERLAP_I5xy_Py_aa = I_TWOBODYOVERLAP_K5x2y_S_aa+ABY*I_TWOBODYOVERLAP_I5xy_S_aa;
  Double I_TWOBODYOVERLAP_I5xz_Py_aa = I_TWOBODYOVERLAP_K5xyz_S_aa+ABY*I_TWOBODYOVERLAP_I5xz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Py_aa = I_TWOBODYOVERLAP_K4x3y_S_aa+ABY*I_TWOBODYOVERLAP_I4x2y_S_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Py_aa = I_TWOBODYOVERLAP_K4x2yz_S_aa+ABY*I_TWOBODYOVERLAP_I4xyz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Py_aa = I_TWOBODYOVERLAP_K4xy2z_S_aa+ABY*I_TWOBODYOVERLAP_I4x2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Py_aa = I_TWOBODYOVERLAP_K3x4y_S_aa+ABY*I_TWOBODYOVERLAP_I3x3y_S_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Py_aa = I_TWOBODYOVERLAP_K3x3yz_S_aa+ABY*I_TWOBODYOVERLAP_I3x2yz_S_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Py_aa = I_TWOBODYOVERLAP_K3x2y2z_S_aa+ABY*I_TWOBODYOVERLAP_I3xy2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Py_aa = I_TWOBODYOVERLAP_K3xy3z_S_aa+ABY*I_TWOBODYOVERLAP_I3x3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Py_aa = I_TWOBODYOVERLAP_K2x5y_S_aa+ABY*I_TWOBODYOVERLAP_I2x4y_S_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Py_aa = I_TWOBODYOVERLAP_K2x4yz_S_aa+ABY*I_TWOBODYOVERLAP_I2x3yz_S_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Py_aa = I_TWOBODYOVERLAP_K2x3y2z_S_aa+ABY*I_TWOBODYOVERLAP_I2x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Py_aa = I_TWOBODYOVERLAP_K2x2y3z_S_aa+ABY*I_TWOBODYOVERLAP_I2xy3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Py_aa = I_TWOBODYOVERLAP_K2xy4z_S_aa+ABY*I_TWOBODYOVERLAP_I2x4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Py_aa = I_TWOBODYOVERLAP_Kx6y_S_aa+ABY*I_TWOBODYOVERLAP_Ix5y_S_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Py_aa = I_TWOBODYOVERLAP_Kx5yz_S_aa+ABY*I_TWOBODYOVERLAP_Ix4yz_S_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Py_aa = I_TWOBODYOVERLAP_Kx4y2z_S_aa+ABY*I_TWOBODYOVERLAP_Ix3y2z_S_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Py_aa = I_TWOBODYOVERLAP_Kx3y3z_S_aa+ABY*I_TWOBODYOVERLAP_Ix2y3z_S_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Py_aa = I_TWOBODYOVERLAP_Kx2y4z_S_aa+ABY*I_TWOBODYOVERLAP_Ixy4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Py_aa = I_TWOBODYOVERLAP_Kxy5z_S_aa+ABY*I_TWOBODYOVERLAP_Ix5z_S_aa;
  Double I_TWOBODYOVERLAP_I6y_Py_aa = I_TWOBODYOVERLAP_K7y_S_aa+ABY*I_TWOBODYOVERLAP_I6y_S_aa;
  Double I_TWOBODYOVERLAP_I5yz_Py_aa = I_TWOBODYOVERLAP_K6yz_S_aa+ABY*I_TWOBODYOVERLAP_I5yz_S_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Py_aa = I_TWOBODYOVERLAP_K5y2z_S_aa+ABY*I_TWOBODYOVERLAP_I4y2z_S_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Py_aa = I_TWOBODYOVERLAP_K4y3z_S_aa+ABY*I_TWOBODYOVERLAP_I3y3z_S_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Py_aa = I_TWOBODYOVERLAP_K3y4z_S_aa+ABY*I_TWOBODYOVERLAP_I2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Py_aa = I_TWOBODYOVERLAP_K2y5z_S_aa+ABY*I_TWOBODYOVERLAP_Iy5z_S_aa;
  Double I_TWOBODYOVERLAP_I6z_Py_aa = I_TWOBODYOVERLAP_Ky6z_S_aa+ABY*I_TWOBODYOVERLAP_I6z_S_aa;
  Double I_TWOBODYOVERLAP_I6x_Pz_aa = I_TWOBODYOVERLAP_K6xz_S_aa+ABZ*I_TWOBODYOVERLAP_I6x_S_aa;
  Double I_TWOBODYOVERLAP_I5xy_Pz_aa = I_TWOBODYOVERLAP_K5xyz_S_aa+ABZ*I_TWOBODYOVERLAP_I5xy_S_aa;
  Double I_TWOBODYOVERLAP_I5xz_Pz_aa = I_TWOBODYOVERLAP_K5x2z_S_aa+ABZ*I_TWOBODYOVERLAP_I5xz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Pz_aa = I_TWOBODYOVERLAP_K4x2yz_S_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_S_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Pz_aa = I_TWOBODYOVERLAP_K4xy2z_S_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Pz_aa = I_TWOBODYOVERLAP_K4x3z_S_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Pz_aa = I_TWOBODYOVERLAP_K3x3yz_S_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_S_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Pz_aa = I_TWOBODYOVERLAP_K3x2y2z_S_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_S_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Pz_aa = I_TWOBODYOVERLAP_K3xy3z_S_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Pz_aa = I_TWOBODYOVERLAP_K3x4z_S_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Pz_aa = I_TWOBODYOVERLAP_K2x4yz_S_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_S_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Pz_aa = I_TWOBODYOVERLAP_K2x3y2z_S_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_S_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Pz_aa = I_TWOBODYOVERLAP_K2x2y3z_S_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Pz_aa = I_TWOBODYOVERLAP_K2xy4z_S_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Pz_aa = I_TWOBODYOVERLAP_K2x5z_S_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Pz_aa = I_TWOBODYOVERLAP_Kx5yz_S_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_S_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Pz_aa = I_TWOBODYOVERLAP_Kx4y2z_S_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_S_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Pz_aa = I_TWOBODYOVERLAP_Kx3y3z_S_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_S_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Pz_aa = I_TWOBODYOVERLAP_Kx2y4z_S_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_S_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Pz_aa = I_TWOBODYOVERLAP_Kxy5z_S_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Pz_aa = I_TWOBODYOVERLAP_Kx6z_S_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_S_aa;
  Double I_TWOBODYOVERLAP_I6y_Pz_aa = I_TWOBODYOVERLAP_K6yz_S_aa+ABZ*I_TWOBODYOVERLAP_I6y_S_aa;
  Double I_TWOBODYOVERLAP_I5yz_Pz_aa = I_TWOBODYOVERLAP_K5y2z_S_aa+ABZ*I_TWOBODYOVERLAP_I5yz_S_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Pz_aa = I_TWOBODYOVERLAP_K4y3z_S_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_S_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Pz_aa = I_TWOBODYOVERLAP_K3y4z_S_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_S_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Pz_aa = I_TWOBODYOVERLAP_K2y5z_S_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Pz_aa = I_TWOBODYOVERLAP_Ky6z_S_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_S_aa;
  Double I_TWOBODYOVERLAP_I6z_Pz_aa = I_TWOBODYOVERLAP_K7z_S_aa+ABZ*I_TWOBODYOVERLAP_I6z_S_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_K_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 9 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_K7x_Px_aa = I_TWOBODYOVERLAP_L8x_S_aa+ABX*I_TWOBODYOVERLAP_K7x_S_aa;
  Double I_TWOBODYOVERLAP_K6xy_Px_aa = I_TWOBODYOVERLAP_L7xy_S_aa+ABX*I_TWOBODYOVERLAP_K6xy_S_aa;
  Double I_TWOBODYOVERLAP_K6xz_Px_aa = I_TWOBODYOVERLAP_L7xz_S_aa+ABX*I_TWOBODYOVERLAP_K6xz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2y_Px_aa = I_TWOBODYOVERLAP_L6x2y_S_aa+ABX*I_TWOBODYOVERLAP_K5x2y_S_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Px_aa = I_TWOBODYOVERLAP_L6xyz_S_aa+ABX*I_TWOBODYOVERLAP_K5xyz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Px_aa = I_TWOBODYOVERLAP_L6x2z_S_aa+ABX*I_TWOBODYOVERLAP_K5x2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3y_Px_aa = I_TWOBODYOVERLAP_L5x3y_S_aa+ABX*I_TWOBODYOVERLAP_K4x3y_S_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Px_aa = I_TWOBODYOVERLAP_L5x2yz_S_aa+ABX*I_TWOBODYOVERLAP_K4x2yz_S_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Px_aa = I_TWOBODYOVERLAP_L5xy2z_S_aa+ABX*I_TWOBODYOVERLAP_K4xy2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Px_aa = I_TWOBODYOVERLAP_L5x3z_S_aa+ABX*I_TWOBODYOVERLAP_K4x3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4y_Px_aa = I_TWOBODYOVERLAP_L4x4y_S_aa+ABX*I_TWOBODYOVERLAP_K3x4y_S_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Px_aa = I_TWOBODYOVERLAP_L4x3yz_S_aa+ABX*I_TWOBODYOVERLAP_K3x3yz_S_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Px_aa = I_TWOBODYOVERLAP_L4x2y2z_S_aa+ABX*I_TWOBODYOVERLAP_K3x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Px_aa = I_TWOBODYOVERLAP_L4xy3z_S_aa+ABX*I_TWOBODYOVERLAP_K3xy3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Px_aa = I_TWOBODYOVERLAP_L4x4z_S_aa+ABX*I_TWOBODYOVERLAP_K3x4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5y_Px_aa = I_TWOBODYOVERLAP_L3x5y_S_aa+ABX*I_TWOBODYOVERLAP_K2x5y_S_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Px_aa = I_TWOBODYOVERLAP_L3x4yz_S_aa+ABX*I_TWOBODYOVERLAP_K2x4yz_S_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Px_aa = I_TWOBODYOVERLAP_L3x3y2z_S_aa+ABX*I_TWOBODYOVERLAP_K2x3y2z_S_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Px_aa = I_TWOBODYOVERLAP_L3x2y3z_S_aa+ABX*I_TWOBODYOVERLAP_K2x2y3z_S_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Px_aa = I_TWOBODYOVERLAP_L3xy4z_S_aa+ABX*I_TWOBODYOVERLAP_K2xy4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Px_aa = I_TWOBODYOVERLAP_L3x5z_S_aa+ABX*I_TWOBODYOVERLAP_K2x5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6y_Px_aa = I_TWOBODYOVERLAP_L2x6y_S_aa+ABX*I_TWOBODYOVERLAP_Kx6y_S_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Px_aa = I_TWOBODYOVERLAP_L2x5yz_S_aa+ABX*I_TWOBODYOVERLAP_Kx5yz_S_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Px_aa = I_TWOBODYOVERLAP_L2x4y2z_S_aa+ABX*I_TWOBODYOVERLAP_Kx4y2z_S_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Px_aa = I_TWOBODYOVERLAP_L2x3y3z_S_aa+ABX*I_TWOBODYOVERLAP_Kx3y3z_S_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Px_aa = I_TWOBODYOVERLAP_L2x2y4z_S_aa+ABX*I_TWOBODYOVERLAP_Kx2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Px_aa = I_TWOBODYOVERLAP_L2xy5z_S_aa+ABX*I_TWOBODYOVERLAP_Kxy5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Px_aa = I_TWOBODYOVERLAP_L2x6z_S_aa+ABX*I_TWOBODYOVERLAP_Kx6z_S_aa;
  Double I_TWOBODYOVERLAP_K7y_Px_aa = I_TWOBODYOVERLAP_Lx7y_S_aa+ABX*I_TWOBODYOVERLAP_K7y_S_aa;
  Double I_TWOBODYOVERLAP_K6yz_Px_aa = I_TWOBODYOVERLAP_Lx6yz_S_aa+ABX*I_TWOBODYOVERLAP_K6yz_S_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Px_aa = I_TWOBODYOVERLAP_Lx5y2z_S_aa+ABX*I_TWOBODYOVERLAP_K5y2z_S_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Px_aa = I_TWOBODYOVERLAP_Lx4y3z_S_aa+ABX*I_TWOBODYOVERLAP_K4y3z_S_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Px_aa = I_TWOBODYOVERLAP_Lx3y4z_S_aa+ABX*I_TWOBODYOVERLAP_K3y4z_S_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Px_aa = I_TWOBODYOVERLAP_Lx2y5z_S_aa+ABX*I_TWOBODYOVERLAP_K2y5z_S_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Px_aa = I_TWOBODYOVERLAP_Lxy6z_S_aa+ABX*I_TWOBODYOVERLAP_Ky6z_S_aa;
  Double I_TWOBODYOVERLAP_K7z_Px_aa = I_TWOBODYOVERLAP_Lx7z_S_aa+ABX*I_TWOBODYOVERLAP_K7z_S_aa;
  Double I_TWOBODYOVERLAP_K6xy_Py_aa = I_TWOBODYOVERLAP_L6x2y_S_aa+ABY*I_TWOBODYOVERLAP_K6xy_S_aa;
  Double I_TWOBODYOVERLAP_K6xz_Py_aa = I_TWOBODYOVERLAP_L6xyz_S_aa+ABY*I_TWOBODYOVERLAP_K6xz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2y_Py_aa = I_TWOBODYOVERLAP_L5x3y_S_aa+ABY*I_TWOBODYOVERLAP_K5x2y_S_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Py_aa = I_TWOBODYOVERLAP_L5x2yz_S_aa+ABY*I_TWOBODYOVERLAP_K5xyz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Py_aa = I_TWOBODYOVERLAP_L5xy2z_S_aa+ABY*I_TWOBODYOVERLAP_K5x2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3y_Py_aa = I_TWOBODYOVERLAP_L4x4y_S_aa+ABY*I_TWOBODYOVERLAP_K4x3y_S_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Py_aa = I_TWOBODYOVERLAP_L4x3yz_S_aa+ABY*I_TWOBODYOVERLAP_K4x2yz_S_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Py_aa = I_TWOBODYOVERLAP_L4x2y2z_S_aa+ABY*I_TWOBODYOVERLAP_K4xy2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Py_aa = I_TWOBODYOVERLAP_L4xy3z_S_aa+ABY*I_TWOBODYOVERLAP_K4x3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4y_Py_aa = I_TWOBODYOVERLAP_L3x5y_S_aa+ABY*I_TWOBODYOVERLAP_K3x4y_S_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Py_aa = I_TWOBODYOVERLAP_L3x4yz_S_aa+ABY*I_TWOBODYOVERLAP_K3x3yz_S_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Py_aa = I_TWOBODYOVERLAP_L3x3y2z_S_aa+ABY*I_TWOBODYOVERLAP_K3x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Py_aa = I_TWOBODYOVERLAP_L3x2y3z_S_aa+ABY*I_TWOBODYOVERLAP_K3xy3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Py_aa = I_TWOBODYOVERLAP_L3xy4z_S_aa+ABY*I_TWOBODYOVERLAP_K3x4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5y_Py_aa = I_TWOBODYOVERLAP_L2x6y_S_aa+ABY*I_TWOBODYOVERLAP_K2x5y_S_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Py_aa = I_TWOBODYOVERLAP_L2x5yz_S_aa+ABY*I_TWOBODYOVERLAP_K2x4yz_S_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Py_aa = I_TWOBODYOVERLAP_L2x4y2z_S_aa+ABY*I_TWOBODYOVERLAP_K2x3y2z_S_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Py_aa = I_TWOBODYOVERLAP_L2x3y3z_S_aa+ABY*I_TWOBODYOVERLAP_K2x2y3z_S_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Py_aa = I_TWOBODYOVERLAP_L2x2y4z_S_aa+ABY*I_TWOBODYOVERLAP_K2xy4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Py_aa = I_TWOBODYOVERLAP_L2xy5z_S_aa+ABY*I_TWOBODYOVERLAP_K2x5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6y_Py_aa = I_TWOBODYOVERLAP_Lx7y_S_aa+ABY*I_TWOBODYOVERLAP_Kx6y_S_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Py_aa = I_TWOBODYOVERLAP_Lx6yz_S_aa+ABY*I_TWOBODYOVERLAP_Kx5yz_S_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Py_aa = I_TWOBODYOVERLAP_Lx5y2z_S_aa+ABY*I_TWOBODYOVERLAP_Kx4y2z_S_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Py_aa = I_TWOBODYOVERLAP_Lx4y3z_S_aa+ABY*I_TWOBODYOVERLAP_Kx3y3z_S_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Py_aa = I_TWOBODYOVERLAP_Lx3y4z_S_aa+ABY*I_TWOBODYOVERLAP_Kx2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Py_aa = I_TWOBODYOVERLAP_Lx2y5z_S_aa+ABY*I_TWOBODYOVERLAP_Kxy5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Py_aa = I_TWOBODYOVERLAP_Lxy6z_S_aa+ABY*I_TWOBODYOVERLAP_Kx6z_S_aa;
  Double I_TWOBODYOVERLAP_K7y_Py_aa = I_TWOBODYOVERLAP_L8y_S_aa+ABY*I_TWOBODYOVERLAP_K7y_S_aa;
  Double I_TWOBODYOVERLAP_K6yz_Py_aa = I_TWOBODYOVERLAP_L7yz_S_aa+ABY*I_TWOBODYOVERLAP_K6yz_S_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Py_aa = I_TWOBODYOVERLAP_L6y2z_S_aa+ABY*I_TWOBODYOVERLAP_K5y2z_S_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Py_aa = I_TWOBODYOVERLAP_L5y3z_S_aa+ABY*I_TWOBODYOVERLAP_K4y3z_S_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Py_aa = I_TWOBODYOVERLAP_L4y4z_S_aa+ABY*I_TWOBODYOVERLAP_K3y4z_S_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Py_aa = I_TWOBODYOVERLAP_L3y5z_S_aa+ABY*I_TWOBODYOVERLAP_K2y5z_S_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Py_aa = I_TWOBODYOVERLAP_L2y6z_S_aa+ABY*I_TWOBODYOVERLAP_Ky6z_S_aa;
  Double I_TWOBODYOVERLAP_K7z_Py_aa = I_TWOBODYOVERLAP_Ly7z_S_aa+ABY*I_TWOBODYOVERLAP_K7z_S_aa;
  Double I_TWOBODYOVERLAP_K6xz_Pz_aa = I_TWOBODYOVERLAP_L6x2z_S_aa+ABZ*I_TWOBODYOVERLAP_K6xz_S_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Pz_aa = I_TWOBODYOVERLAP_L5xy2z_S_aa+ABZ*I_TWOBODYOVERLAP_K5xyz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Pz_aa = I_TWOBODYOVERLAP_L5x3z_S_aa+ABZ*I_TWOBODYOVERLAP_K5x2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Pz_aa = I_TWOBODYOVERLAP_L4x2y2z_S_aa+ABZ*I_TWOBODYOVERLAP_K4x2yz_S_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Pz_aa = I_TWOBODYOVERLAP_L4xy3z_S_aa+ABZ*I_TWOBODYOVERLAP_K4xy2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Pz_aa = I_TWOBODYOVERLAP_L4x4z_S_aa+ABZ*I_TWOBODYOVERLAP_K4x3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Pz_aa = I_TWOBODYOVERLAP_L3x3y2z_S_aa+ABZ*I_TWOBODYOVERLAP_K3x3yz_S_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Pz_aa = I_TWOBODYOVERLAP_L3x2y3z_S_aa+ABZ*I_TWOBODYOVERLAP_K3x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Pz_aa = I_TWOBODYOVERLAP_L3xy4z_S_aa+ABZ*I_TWOBODYOVERLAP_K3xy3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Pz_aa = I_TWOBODYOVERLAP_L3x5z_S_aa+ABZ*I_TWOBODYOVERLAP_K3x4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Pz_aa = I_TWOBODYOVERLAP_L2x4y2z_S_aa+ABZ*I_TWOBODYOVERLAP_K2x4yz_S_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Pz_aa = I_TWOBODYOVERLAP_L2x3y3z_S_aa+ABZ*I_TWOBODYOVERLAP_K2x3y2z_S_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Pz_aa = I_TWOBODYOVERLAP_L2x2y4z_S_aa+ABZ*I_TWOBODYOVERLAP_K2x2y3z_S_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Pz_aa = I_TWOBODYOVERLAP_L2xy5z_S_aa+ABZ*I_TWOBODYOVERLAP_K2xy4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Pz_aa = I_TWOBODYOVERLAP_L2x6z_S_aa+ABZ*I_TWOBODYOVERLAP_K2x5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Pz_aa = I_TWOBODYOVERLAP_Lx5y2z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx5yz_S_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Pz_aa = I_TWOBODYOVERLAP_Lx4y3z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx4y2z_S_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Pz_aa = I_TWOBODYOVERLAP_Lx3y4z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx3y3z_S_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Pz_aa = I_TWOBODYOVERLAP_Lx2y5z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Pz_aa = I_TWOBODYOVERLAP_Lxy6z_S_aa+ABZ*I_TWOBODYOVERLAP_Kxy5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Pz_aa = I_TWOBODYOVERLAP_Lx7z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx6z_S_aa;
  Double I_TWOBODYOVERLAP_K6yz_Pz_aa = I_TWOBODYOVERLAP_L6y2z_S_aa+ABZ*I_TWOBODYOVERLAP_K6yz_S_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Pz_aa = I_TWOBODYOVERLAP_L5y3z_S_aa+ABZ*I_TWOBODYOVERLAP_K5y2z_S_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Pz_aa = I_TWOBODYOVERLAP_L4y4z_S_aa+ABZ*I_TWOBODYOVERLAP_K4y3z_S_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Pz_aa = I_TWOBODYOVERLAP_L3y5z_S_aa+ABZ*I_TWOBODYOVERLAP_K3y4z_S_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Pz_aa = I_TWOBODYOVERLAP_L2y6z_S_aa+ABZ*I_TWOBODYOVERLAP_K2y5z_S_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Pz_aa = I_TWOBODYOVERLAP_Ly7z_S_aa+ABZ*I_TWOBODYOVERLAP_Ky6z_S_aa;
  Double I_TWOBODYOVERLAP_K7z_Pz_aa = I_TWOBODYOVERLAP_L8z_S_aa+ABZ*I_TWOBODYOVERLAP_K7z_S_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_D_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_D2x_aa = I_TWOBODYOVERLAP_K7x_Px_aa+ABX*I_TWOBODYOVERLAP_I6x_Px_aa;
  Double I_TWOBODYOVERLAP_I5xy_D2x_aa = I_TWOBODYOVERLAP_K6xy_Px_aa+ABX*I_TWOBODYOVERLAP_I5xy_Px_aa;
  Double I_TWOBODYOVERLAP_I5xz_D2x_aa = I_TWOBODYOVERLAP_K6xz_Px_aa+ABX*I_TWOBODYOVERLAP_I5xz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2y_D2x_aa = I_TWOBODYOVERLAP_K5x2y_Px_aa+ABX*I_TWOBODYOVERLAP_I4x2y_Px_aa;
  Double I_TWOBODYOVERLAP_I4xyz_D2x_aa = I_TWOBODYOVERLAP_K5xyz_Px_aa+ABX*I_TWOBODYOVERLAP_I4xyz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2z_D2x_aa = I_TWOBODYOVERLAP_K5x2z_Px_aa+ABX*I_TWOBODYOVERLAP_I4x2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3y_D2x_aa = I_TWOBODYOVERLAP_K4x3y_Px_aa+ABX*I_TWOBODYOVERLAP_I3x3y_Px_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_D2x_aa = I_TWOBODYOVERLAP_K4x2yz_Px_aa+ABX*I_TWOBODYOVERLAP_I3x2yz_Px_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_D2x_aa = I_TWOBODYOVERLAP_K4xy2z_Px_aa+ABX*I_TWOBODYOVERLAP_I3xy2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3z_D2x_aa = I_TWOBODYOVERLAP_K4x3z_Px_aa+ABX*I_TWOBODYOVERLAP_I3x3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4y_D2x_aa = I_TWOBODYOVERLAP_K3x4y_Px_aa+ABX*I_TWOBODYOVERLAP_I2x4y_Px_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_D2x_aa = I_TWOBODYOVERLAP_K3x3yz_Px_aa+ABX*I_TWOBODYOVERLAP_I2x3yz_Px_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2x_aa = I_TWOBODYOVERLAP_K3x2y2z_Px_aa+ABX*I_TWOBODYOVERLAP_I2x2y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_D2x_aa = I_TWOBODYOVERLAP_K3xy3z_Px_aa+ABX*I_TWOBODYOVERLAP_I2xy3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4z_D2x_aa = I_TWOBODYOVERLAP_K3x4z_Px_aa+ABX*I_TWOBODYOVERLAP_I2x4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5y_D2x_aa = I_TWOBODYOVERLAP_K2x5y_Px_aa+ABX*I_TWOBODYOVERLAP_Ix5y_Px_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_D2x_aa = I_TWOBODYOVERLAP_K2x4yz_Px_aa+ABX*I_TWOBODYOVERLAP_Ix4yz_Px_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2x_aa = I_TWOBODYOVERLAP_K2x3y2z_Px_aa+ABX*I_TWOBODYOVERLAP_Ix3y2z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2x_aa = I_TWOBODYOVERLAP_K2x2y3z_Px_aa+ABX*I_TWOBODYOVERLAP_Ix2y3z_Px_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_D2x_aa = I_TWOBODYOVERLAP_K2xy4z_Px_aa+ABX*I_TWOBODYOVERLAP_Ixy4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5z_D2x_aa = I_TWOBODYOVERLAP_K2x5z_Px_aa+ABX*I_TWOBODYOVERLAP_Ix5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6y_D2x_aa = I_TWOBODYOVERLAP_Kx6y_Px_aa+ABX*I_TWOBODYOVERLAP_I6y_Px_aa;
  Double I_TWOBODYOVERLAP_I5yz_D2x_aa = I_TWOBODYOVERLAP_Kx5yz_Px_aa+ABX*I_TWOBODYOVERLAP_I5yz_Px_aa;
  Double I_TWOBODYOVERLAP_I4y2z_D2x_aa = I_TWOBODYOVERLAP_Kx4y2z_Px_aa+ABX*I_TWOBODYOVERLAP_I4y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3y3z_D2x_aa = I_TWOBODYOVERLAP_Kx3y3z_Px_aa+ABX*I_TWOBODYOVERLAP_I3y3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2y4z_D2x_aa = I_TWOBODYOVERLAP_Kx2y4z_Px_aa+ABX*I_TWOBODYOVERLAP_I2y4z_Px_aa;
  Double I_TWOBODYOVERLAP_Iy5z_D2x_aa = I_TWOBODYOVERLAP_Kxy5z_Px_aa+ABX*I_TWOBODYOVERLAP_Iy5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6z_D2x_aa = I_TWOBODYOVERLAP_Kx6z_Px_aa+ABX*I_TWOBODYOVERLAP_I6z_Px_aa;
  Double I_TWOBODYOVERLAP_I6x_Dxy_aa = I_TWOBODYOVERLAP_K6xy_Px_aa+ABY*I_TWOBODYOVERLAP_I6x_Px_aa;
  Double I_TWOBODYOVERLAP_I5xy_Dxy_aa = I_TWOBODYOVERLAP_K5x2y_Px_aa+ABY*I_TWOBODYOVERLAP_I5xy_Px_aa;
  Double I_TWOBODYOVERLAP_I5xz_Dxy_aa = I_TWOBODYOVERLAP_K5xyz_Px_aa+ABY*I_TWOBODYOVERLAP_I5xz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Dxy_aa = I_TWOBODYOVERLAP_K4x3y_Px_aa+ABY*I_TWOBODYOVERLAP_I4x2y_Px_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Dxy_aa = I_TWOBODYOVERLAP_K4x2yz_Px_aa+ABY*I_TWOBODYOVERLAP_I4xyz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Dxy_aa = I_TWOBODYOVERLAP_K4xy2z_Px_aa+ABY*I_TWOBODYOVERLAP_I4x2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Dxy_aa = I_TWOBODYOVERLAP_K3x4y_Px_aa+ABY*I_TWOBODYOVERLAP_I3x3y_Px_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Dxy_aa = I_TWOBODYOVERLAP_K3x3yz_Px_aa+ABY*I_TWOBODYOVERLAP_I3x2yz_Px_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Dxy_aa = I_TWOBODYOVERLAP_K3x2y2z_Px_aa+ABY*I_TWOBODYOVERLAP_I3xy2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Dxy_aa = I_TWOBODYOVERLAP_K3xy3z_Px_aa+ABY*I_TWOBODYOVERLAP_I3x3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Dxy_aa = I_TWOBODYOVERLAP_K2x5y_Px_aa+ABY*I_TWOBODYOVERLAP_I2x4y_Px_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Dxy_aa = I_TWOBODYOVERLAP_K2x4yz_Px_aa+ABY*I_TWOBODYOVERLAP_I2x3yz_Px_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Dxy_aa = I_TWOBODYOVERLAP_K2x3y2z_Px_aa+ABY*I_TWOBODYOVERLAP_I2x2y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Dxy_aa = I_TWOBODYOVERLAP_K2x2y3z_Px_aa+ABY*I_TWOBODYOVERLAP_I2xy3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Dxy_aa = I_TWOBODYOVERLAP_K2xy4z_Px_aa+ABY*I_TWOBODYOVERLAP_I2x4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Dxy_aa = I_TWOBODYOVERLAP_Kx6y_Px_aa+ABY*I_TWOBODYOVERLAP_Ix5y_Px_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Dxy_aa = I_TWOBODYOVERLAP_Kx5yz_Px_aa+ABY*I_TWOBODYOVERLAP_Ix4yz_Px_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Dxy_aa = I_TWOBODYOVERLAP_Kx4y2z_Px_aa+ABY*I_TWOBODYOVERLAP_Ix3y2z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Dxy_aa = I_TWOBODYOVERLAP_Kx3y3z_Px_aa+ABY*I_TWOBODYOVERLAP_Ix2y3z_Px_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Dxy_aa = I_TWOBODYOVERLAP_Kx2y4z_Px_aa+ABY*I_TWOBODYOVERLAP_Ixy4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Dxy_aa = I_TWOBODYOVERLAP_Kxy5z_Px_aa+ABY*I_TWOBODYOVERLAP_Ix5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6y_Dxy_aa = I_TWOBODYOVERLAP_K7y_Px_aa+ABY*I_TWOBODYOVERLAP_I6y_Px_aa;
  Double I_TWOBODYOVERLAP_I5yz_Dxy_aa = I_TWOBODYOVERLAP_K6yz_Px_aa+ABY*I_TWOBODYOVERLAP_I5yz_Px_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Dxy_aa = I_TWOBODYOVERLAP_K5y2z_Px_aa+ABY*I_TWOBODYOVERLAP_I4y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Dxy_aa = I_TWOBODYOVERLAP_K4y3z_Px_aa+ABY*I_TWOBODYOVERLAP_I3y3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Dxy_aa = I_TWOBODYOVERLAP_K3y4z_Px_aa+ABY*I_TWOBODYOVERLAP_I2y4z_Px_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Dxy_aa = I_TWOBODYOVERLAP_K2y5z_Px_aa+ABY*I_TWOBODYOVERLAP_Iy5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6z_Dxy_aa = I_TWOBODYOVERLAP_Ky6z_Px_aa+ABY*I_TWOBODYOVERLAP_I6z_Px_aa;
  Double I_TWOBODYOVERLAP_I6x_Dxz_aa = I_TWOBODYOVERLAP_K6xz_Px_aa+ABZ*I_TWOBODYOVERLAP_I6x_Px_aa;
  Double I_TWOBODYOVERLAP_I5xy_Dxz_aa = I_TWOBODYOVERLAP_K5xyz_Px_aa+ABZ*I_TWOBODYOVERLAP_I5xy_Px_aa;
  Double I_TWOBODYOVERLAP_I5xz_Dxz_aa = I_TWOBODYOVERLAP_K5x2z_Px_aa+ABZ*I_TWOBODYOVERLAP_I5xz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Dxz_aa = I_TWOBODYOVERLAP_K4x2yz_Px_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_Px_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Dxz_aa = I_TWOBODYOVERLAP_K4xy2z_Px_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Dxz_aa = I_TWOBODYOVERLAP_K4x3z_Px_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Dxz_aa = I_TWOBODYOVERLAP_K3x3yz_Px_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_Px_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Dxz_aa = I_TWOBODYOVERLAP_K3x2y2z_Px_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_Px_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Dxz_aa = I_TWOBODYOVERLAP_K3xy3z_Px_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Dxz_aa = I_TWOBODYOVERLAP_K3x4z_Px_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Dxz_aa = I_TWOBODYOVERLAP_K2x4yz_Px_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_Px_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Dxz_aa = I_TWOBODYOVERLAP_K2x3y2z_Px_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_Px_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Dxz_aa = I_TWOBODYOVERLAP_K2x2y3z_Px_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Dxz_aa = I_TWOBODYOVERLAP_K2xy4z_Px_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Dxz_aa = I_TWOBODYOVERLAP_K2x5z_Px_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Dxz_aa = I_TWOBODYOVERLAP_Kx5yz_Px_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_Px_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Dxz_aa = I_TWOBODYOVERLAP_Kx4y2z_Px_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_Px_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Dxz_aa = I_TWOBODYOVERLAP_Kx3y3z_Px_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Dxz_aa = I_TWOBODYOVERLAP_Kx2y4z_Px_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_Px_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Dxz_aa = I_TWOBODYOVERLAP_Kxy5z_Px_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Dxz_aa = I_TWOBODYOVERLAP_Kx6z_Px_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6y_Dxz_aa = I_TWOBODYOVERLAP_K6yz_Px_aa+ABZ*I_TWOBODYOVERLAP_I6y_Px_aa;
  Double I_TWOBODYOVERLAP_I5yz_Dxz_aa = I_TWOBODYOVERLAP_K5y2z_Px_aa+ABZ*I_TWOBODYOVERLAP_I5yz_Px_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Dxz_aa = I_TWOBODYOVERLAP_K4y3z_Px_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Dxz_aa = I_TWOBODYOVERLAP_K3y4z_Px_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Dxz_aa = I_TWOBODYOVERLAP_K2y5z_Px_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_Px_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Dxz_aa = I_TWOBODYOVERLAP_Ky6z_Px_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6z_Dxz_aa = I_TWOBODYOVERLAP_K7z_Px_aa+ABZ*I_TWOBODYOVERLAP_I6z_Px_aa;
  Double I_TWOBODYOVERLAP_I6x_D2y_aa = I_TWOBODYOVERLAP_K6xy_Py_aa+ABY*I_TWOBODYOVERLAP_I6x_Py_aa;
  Double I_TWOBODYOVERLAP_I5xy_D2y_aa = I_TWOBODYOVERLAP_K5x2y_Py_aa+ABY*I_TWOBODYOVERLAP_I5xy_Py_aa;
  Double I_TWOBODYOVERLAP_I5xz_D2y_aa = I_TWOBODYOVERLAP_K5xyz_Py_aa+ABY*I_TWOBODYOVERLAP_I5xz_Py_aa;
  Double I_TWOBODYOVERLAP_I4x2y_D2y_aa = I_TWOBODYOVERLAP_K4x3y_Py_aa+ABY*I_TWOBODYOVERLAP_I4x2y_Py_aa;
  Double I_TWOBODYOVERLAP_I4xyz_D2y_aa = I_TWOBODYOVERLAP_K4x2yz_Py_aa+ABY*I_TWOBODYOVERLAP_I4xyz_Py_aa;
  Double I_TWOBODYOVERLAP_I4x2z_D2y_aa = I_TWOBODYOVERLAP_K4xy2z_Py_aa+ABY*I_TWOBODYOVERLAP_I4x2z_Py_aa;
  Double I_TWOBODYOVERLAP_I3x3y_D2y_aa = I_TWOBODYOVERLAP_K3x4y_Py_aa+ABY*I_TWOBODYOVERLAP_I3x3y_Py_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_D2y_aa = I_TWOBODYOVERLAP_K3x3yz_Py_aa+ABY*I_TWOBODYOVERLAP_I3x2yz_Py_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_D2y_aa = I_TWOBODYOVERLAP_K3x2y2z_Py_aa+ABY*I_TWOBODYOVERLAP_I3xy2z_Py_aa;
  Double I_TWOBODYOVERLAP_I3x3z_D2y_aa = I_TWOBODYOVERLAP_K3xy3z_Py_aa+ABY*I_TWOBODYOVERLAP_I3x3z_Py_aa;
  Double I_TWOBODYOVERLAP_I2x4y_D2y_aa = I_TWOBODYOVERLAP_K2x5y_Py_aa+ABY*I_TWOBODYOVERLAP_I2x4y_Py_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_D2y_aa = I_TWOBODYOVERLAP_K2x4yz_Py_aa+ABY*I_TWOBODYOVERLAP_I2x3yz_Py_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2y_aa = I_TWOBODYOVERLAP_K2x3y2z_Py_aa+ABY*I_TWOBODYOVERLAP_I2x2y2z_Py_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_D2y_aa = I_TWOBODYOVERLAP_K2x2y3z_Py_aa+ABY*I_TWOBODYOVERLAP_I2xy3z_Py_aa;
  Double I_TWOBODYOVERLAP_I2x4z_D2y_aa = I_TWOBODYOVERLAP_K2xy4z_Py_aa+ABY*I_TWOBODYOVERLAP_I2x4z_Py_aa;
  Double I_TWOBODYOVERLAP_Ix5y_D2y_aa = I_TWOBODYOVERLAP_Kx6y_Py_aa+ABY*I_TWOBODYOVERLAP_Ix5y_Py_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_D2y_aa = I_TWOBODYOVERLAP_Kx5yz_Py_aa+ABY*I_TWOBODYOVERLAP_Ix4yz_Py_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2y_aa = I_TWOBODYOVERLAP_Kx4y2z_Py_aa+ABY*I_TWOBODYOVERLAP_Ix3y2z_Py_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2y_aa = I_TWOBODYOVERLAP_Kx3y3z_Py_aa+ABY*I_TWOBODYOVERLAP_Ix2y3z_Py_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_D2y_aa = I_TWOBODYOVERLAP_Kx2y4z_Py_aa+ABY*I_TWOBODYOVERLAP_Ixy4z_Py_aa;
  Double I_TWOBODYOVERLAP_Ix5z_D2y_aa = I_TWOBODYOVERLAP_Kxy5z_Py_aa+ABY*I_TWOBODYOVERLAP_Ix5z_Py_aa;
  Double I_TWOBODYOVERLAP_I6y_D2y_aa = I_TWOBODYOVERLAP_K7y_Py_aa+ABY*I_TWOBODYOVERLAP_I6y_Py_aa;
  Double I_TWOBODYOVERLAP_I5yz_D2y_aa = I_TWOBODYOVERLAP_K6yz_Py_aa+ABY*I_TWOBODYOVERLAP_I5yz_Py_aa;
  Double I_TWOBODYOVERLAP_I4y2z_D2y_aa = I_TWOBODYOVERLAP_K5y2z_Py_aa+ABY*I_TWOBODYOVERLAP_I4y2z_Py_aa;
  Double I_TWOBODYOVERLAP_I3y3z_D2y_aa = I_TWOBODYOVERLAP_K4y3z_Py_aa+ABY*I_TWOBODYOVERLAP_I3y3z_Py_aa;
  Double I_TWOBODYOVERLAP_I2y4z_D2y_aa = I_TWOBODYOVERLAP_K3y4z_Py_aa+ABY*I_TWOBODYOVERLAP_I2y4z_Py_aa;
  Double I_TWOBODYOVERLAP_Iy5z_D2y_aa = I_TWOBODYOVERLAP_K2y5z_Py_aa+ABY*I_TWOBODYOVERLAP_Iy5z_Py_aa;
  Double I_TWOBODYOVERLAP_I6z_D2y_aa = I_TWOBODYOVERLAP_Ky6z_Py_aa+ABY*I_TWOBODYOVERLAP_I6z_Py_aa;
  Double I_TWOBODYOVERLAP_I6x_Dyz_aa = I_TWOBODYOVERLAP_K6xz_Py_aa+ABZ*I_TWOBODYOVERLAP_I6x_Py_aa;
  Double I_TWOBODYOVERLAP_I5xy_Dyz_aa = I_TWOBODYOVERLAP_K5xyz_Py_aa+ABZ*I_TWOBODYOVERLAP_I5xy_Py_aa;
  Double I_TWOBODYOVERLAP_I5xz_Dyz_aa = I_TWOBODYOVERLAP_K5x2z_Py_aa+ABZ*I_TWOBODYOVERLAP_I5xz_Py_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Dyz_aa = I_TWOBODYOVERLAP_K4x2yz_Py_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_Py_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Dyz_aa = I_TWOBODYOVERLAP_K4xy2z_Py_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_Py_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Dyz_aa = I_TWOBODYOVERLAP_K4x3z_Py_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_Py_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Dyz_aa = I_TWOBODYOVERLAP_K3x3yz_Py_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_Py_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Dyz_aa = I_TWOBODYOVERLAP_K3x2y2z_Py_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_Py_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Dyz_aa = I_TWOBODYOVERLAP_K3xy3z_Py_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_Py_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Dyz_aa = I_TWOBODYOVERLAP_K3x4z_Py_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_Py_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Dyz_aa = I_TWOBODYOVERLAP_K2x4yz_Py_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_Py_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Dyz_aa = I_TWOBODYOVERLAP_K2x3y2z_Py_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_Py_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Dyz_aa = I_TWOBODYOVERLAP_K2x2y3z_Py_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_Py_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Dyz_aa = I_TWOBODYOVERLAP_K2xy4z_Py_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_Py_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Dyz_aa = I_TWOBODYOVERLAP_K2x5z_Py_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_Py_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Dyz_aa = I_TWOBODYOVERLAP_Kx5yz_Py_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_Py_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Dyz_aa = I_TWOBODYOVERLAP_Kx4y2z_Py_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_Py_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Dyz_aa = I_TWOBODYOVERLAP_Kx3y3z_Py_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_Py_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Dyz_aa = I_TWOBODYOVERLAP_Kx2y4z_Py_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_Py_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Dyz_aa = I_TWOBODYOVERLAP_Kxy5z_Py_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_Py_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Dyz_aa = I_TWOBODYOVERLAP_Kx6z_Py_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_Py_aa;
  Double I_TWOBODYOVERLAP_I6y_Dyz_aa = I_TWOBODYOVERLAP_K6yz_Py_aa+ABZ*I_TWOBODYOVERLAP_I6y_Py_aa;
  Double I_TWOBODYOVERLAP_I5yz_Dyz_aa = I_TWOBODYOVERLAP_K5y2z_Py_aa+ABZ*I_TWOBODYOVERLAP_I5yz_Py_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Dyz_aa = I_TWOBODYOVERLAP_K4y3z_Py_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_Py_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Dyz_aa = I_TWOBODYOVERLAP_K3y4z_Py_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_Py_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Dyz_aa = I_TWOBODYOVERLAP_K2y5z_Py_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_Py_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Dyz_aa = I_TWOBODYOVERLAP_Ky6z_Py_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_Py_aa;
  Double I_TWOBODYOVERLAP_I6z_Dyz_aa = I_TWOBODYOVERLAP_K7z_Py_aa+ABZ*I_TWOBODYOVERLAP_I6z_Py_aa;
  Double I_TWOBODYOVERLAP_I6x_D2z_aa = I_TWOBODYOVERLAP_K6xz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I6x_Pz_aa;
  Double I_TWOBODYOVERLAP_I5xy_D2z_aa = I_TWOBODYOVERLAP_K5xyz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I5xy_Pz_aa;
  Double I_TWOBODYOVERLAP_I5xz_D2z_aa = I_TWOBODYOVERLAP_K5x2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I5xz_Pz_aa;
  Double I_TWOBODYOVERLAP_I4x2y_D2z_aa = I_TWOBODYOVERLAP_K4x2yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_Pz_aa;
  Double I_TWOBODYOVERLAP_I4xyz_D2z_aa = I_TWOBODYOVERLAP_K4xy2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_Pz_aa;
  Double I_TWOBODYOVERLAP_I4x2z_D2z_aa = I_TWOBODYOVERLAP_K4x3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_Pz_aa;
  Double I_TWOBODYOVERLAP_I3x3y_D2z_aa = I_TWOBODYOVERLAP_K3x3yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_Pz_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_D2z_aa = I_TWOBODYOVERLAP_K3x2y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_Pz_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_D2z_aa = I_TWOBODYOVERLAP_K3xy3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_Pz_aa;
  Double I_TWOBODYOVERLAP_I3x3z_D2z_aa = I_TWOBODYOVERLAP_K3x4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_Pz_aa;
  Double I_TWOBODYOVERLAP_I2x4y_D2z_aa = I_TWOBODYOVERLAP_K2x4yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_Pz_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_D2z_aa = I_TWOBODYOVERLAP_K2x3y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_Pz_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2z_aa = I_TWOBODYOVERLAP_K2x2y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_D2z_aa = I_TWOBODYOVERLAP_K2xy4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_Pz_aa;
  Double I_TWOBODYOVERLAP_I2x4z_D2z_aa = I_TWOBODYOVERLAP_K2x5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix5y_D2z_aa = I_TWOBODYOVERLAP_Kx5yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_D2z_aa = I_TWOBODYOVERLAP_Kx4y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2z_aa = I_TWOBODYOVERLAP_Kx3y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2z_aa = I_TWOBODYOVERLAP_Kx2y4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_Pz_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_D2z_aa = I_TWOBODYOVERLAP_Kxy5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix5z_D2z_aa = I_TWOBODYOVERLAP_Kx6z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_Pz_aa;
  Double I_TWOBODYOVERLAP_I6y_D2z_aa = I_TWOBODYOVERLAP_K6yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I6y_Pz_aa;
  Double I_TWOBODYOVERLAP_I5yz_D2z_aa = I_TWOBODYOVERLAP_K5y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I5yz_Pz_aa;
  Double I_TWOBODYOVERLAP_I4y2z_D2z_aa = I_TWOBODYOVERLAP_K4y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_I3y3z_D2z_aa = I_TWOBODYOVERLAP_K3y4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_Pz_aa;
  Double I_TWOBODYOVERLAP_I2y4z_D2z_aa = I_TWOBODYOVERLAP_K2y5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_Pz_aa;
  Double I_TWOBODYOVERLAP_Iy5z_D2z_aa = I_TWOBODYOVERLAP_Ky6z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_Pz_aa;
  Double I_TWOBODYOVERLAP_I6z_D2z_aa = I_TWOBODYOVERLAP_K7z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I6z_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
   ************************************************************/
  abcd[0] = 4.0E0*I_TWOBODYOVERLAP_I6x_D2x_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_D2x_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_D2x_a+4*3*I_TWOBODYOVERLAP_D2x_D2x;
  abcd[1] = 4.0E0*I_TWOBODYOVERLAP_I5xy_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_D2x_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_D2x_a+3*2*I_TWOBODYOVERLAP_Dxy_D2x;
  abcd[2] = 4.0E0*I_TWOBODYOVERLAP_I5xz_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_D2x_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_D2x_a+3*2*I_TWOBODYOVERLAP_Dxz_D2x;
  abcd[3] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2x_a+2*1*I_TWOBODYOVERLAP_D2y_D2x;
  abcd[4] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_D2x_a+2*1*I_TWOBODYOVERLAP_Dyz_D2x;
  abcd[5] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2x_a+2*1*I_TWOBODYOVERLAP_D2z_D2x;
  abcd[6] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_D2x_a;
  abcd[7] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2x_a;
  abcd[8] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2x_a;
  abcd[9] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_D2x_a;
  abcd[10] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2x_a;
  abcd[11] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2x_a;
  abcd[12] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_D2x_a;
  abcd[13] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2x_a;
  abcd[14] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2x_a;
  abcd[15] = 4.0E0*I_TWOBODYOVERLAP_I6x_Dxy_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_Dxy_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_Dxy_a+4*3*I_TWOBODYOVERLAP_D2x_Dxy;
  abcd[16] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_Dxy_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Dxy_a+3*2*I_TWOBODYOVERLAP_Dxy_Dxy;
  abcd[17] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_Dxy_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Dxy_a+3*2*I_TWOBODYOVERLAP_Dxz_Dxy;
  abcd[18] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dxy_a+2*1*I_TWOBODYOVERLAP_D2y_Dxy;
  abcd[19] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Dxy_a+2*1*I_TWOBODYOVERLAP_Dyz_Dxy;
  abcd[20] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dxy_a+2*1*I_TWOBODYOVERLAP_D2z_Dxy;
  abcd[21] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Dxy_a;
  abcd[22] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxy_a;
  abcd[23] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxy_a;
  abcd[24] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Dxy_a;
  abcd[25] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dxy_a;
  abcd[26] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dxy_a;
  abcd[27] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Dxy_a;
  abcd[28] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dxy_a;
  abcd[29] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dxy_a;
  abcd[30] = 4.0E0*I_TWOBODYOVERLAP_I6x_Dxz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_Dxz_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_Dxz_a+4*3*I_TWOBODYOVERLAP_D2x_Dxz;
  abcd[31] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_Dxz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Dxz_a+3*2*I_TWOBODYOVERLAP_Dxy_Dxz;
  abcd[32] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_Dxz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Dxz_a+3*2*I_TWOBODYOVERLAP_Dxz_Dxz;
  abcd[33] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dxz_a+2*1*I_TWOBODYOVERLAP_D2y_Dxz;
  abcd[34] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Dxz_a+2*1*I_TWOBODYOVERLAP_Dyz_Dxz;
  abcd[35] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dxz_a+2*1*I_TWOBODYOVERLAP_D2z_Dxz;
  abcd[36] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Dxz_a;
  abcd[37] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxz_a;
  abcd[38] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxz_a;
  abcd[39] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Dxz_a;
  abcd[40] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dxz_a;
  abcd[41] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dxz_a;
  abcd[42] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Dxz_a;
  abcd[43] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dxz_a;
  abcd[44] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dxz_a;
  abcd[45] = 4.0E0*I_TWOBODYOVERLAP_I6x_D2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_D2y_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_D2y_a+4*3*I_TWOBODYOVERLAP_D2x_D2y;
  abcd[46] = 4.0E0*I_TWOBODYOVERLAP_I5xy_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_D2y_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_D2y_a+3*2*I_TWOBODYOVERLAP_Dxy_D2y;
  abcd[47] = 4.0E0*I_TWOBODYOVERLAP_I5xz_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_D2y_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_D2y_a+3*2*I_TWOBODYOVERLAP_Dxz_D2y;
  abcd[48] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2y_a+2*1*I_TWOBODYOVERLAP_D2y_D2y;
  abcd[49] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_D2y_a+2*1*I_TWOBODYOVERLAP_Dyz_D2y;
  abcd[50] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2y_a+2*1*I_TWOBODYOVERLAP_D2z_D2y;
  abcd[51] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_D2y_a;
  abcd[52] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2y_a;
  abcd[53] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2y_a;
  abcd[54] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_D2y_a;
  abcd[55] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2y_a;
  abcd[56] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2y_a;
  abcd[57] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_D2y_a;
  abcd[58] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2y_a;
  abcd[59] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2y_a;
  abcd[60] = 4.0E0*I_TWOBODYOVERLAP_I6x_Dyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_Dyz_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_Dyz_a+4*3*I_TWOBODYOVERLAP_D2x_Dyz;
  abcd[61] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_Dyz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Dyz_a+3*2*I_TWOBODYOVERLAP_Dxy_Dyz;
  abcd[62] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_Dyz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Dyz_a+3*2*I_TWOBODYOVERLAP_Dxz_Dyz;
  abcd[63] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dyz_a+2*1*I_TWOBODYOVERLAP_D2y_Dyz;
  abcd[64] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Dyz_a+2*1*I_TWOBODYOVERLAP_Dyz_Dyz;
  abcd[65] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dyz_a+2*1*I_TWOBODYOVERLAP_D2z_Dyz;
  abcd[66] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Dyz_a;
  abcd[67] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dyz_a;
  abcd[68] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dyz_a;
  abcd[69] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Dyz_a;
  abcd[70] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dyz_a;
  abcd[71] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dyz_a;
  abcd[72] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Dyz_a;
  abcd[73] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dyz_a;
  abcd[74] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dyz_a;
  abcd[75] = 4.0E0*I_TWOBODYOVERLAP_I6x_D2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_D2z_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_D2z_a+4*3*I_TWOBODYOVERLAP_D2x_D2z;
  abcd[76] = 4.0E0*I_TWOBODYOVERLAP_I5xy_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_D2z_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_D2z_a+3*2*I_TWOBODYOVERLAP_Dxy_D2z;
  abcd[77] = 4.0E0*I_TWOBODYOVERLAP_I5xz_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_D2z_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_D2z_a+3*2*I_TWOBODYOVERLAP_Dxz_D2z;
  abcd[78] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2z_a+2*1*I_TWOBODYOVERLAP_D2y_D2z;
  abcd[79] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_D2z_a+2*1*I_TWOBODYOVERLAP_Dyz_D2z;
  abcd[80] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2z_a+2*1*I_TWOBODYOVERLAP_D2z_D2z;
  abcd[81] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_D2z_a;
  abcd[82] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2z_a;
  abcd[83] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2z_a;
  abcd[84] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_D2z_a;
  abcd[85] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2z_a;
  abcd[86] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2z_a;
  abcd[87] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_D2z_a;
  abcd[88] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2z_a;
  abcd[89] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
   ************************************************************/
  abcd[90] = 4.0E0*I_TWOBODYOVERLAP_I5xy_D2x_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_D2x_a;
  abcd[91] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2x_a+3*1*I_TWOBODYOVERLAP_D2x_D2x;
  abcd[92] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_D2x_a;
  abcd[93] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_D2x_a+2*2*I_TWOBODYOVERLAP_Dxy_D2x;
  abcd[94] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2x_a+2*1*I_TWOBODYOVERLAP_Dxz_D2x;
  abcd[95] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2x_a;
  abcd[96] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2x_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2x_a+3*I_TWOBODYOVERLAP_D2y_D2x;
  abcd[97] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2x_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2x_a+2*I_TWOBODYOVERLAP_Dyz_D2x;
  abcd[98] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_D2x_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_D2x_a+1*I_TWOBODYOVERLAP_D2z_D2x;
  abcd[99] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2x_a;
  abcd[100] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_D2x_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_D2x_a;
  abcd[101] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_D2x_a;
  abcd[102] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2x_a;
  abcd[103] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2x_a;
  abcd[104] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2x_aa;
  abcd[105] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Dxy_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Dxy_a;
  abcd[106] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dxy_a+3*1*I_TWOBODYOVERLAP_D2x_Dxy;
  abcd[107] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Dxy_a;
  abcd[108] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Dxy_a+2*2*I_TWOBODYOVERLAP_Dxy_Dxy;
  abcd[109] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxy_a+2*1*I_TWOBODYOVERLAP_Dxz_Dxy;
  abcd[110] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxy_a;
  abcd[111] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dxy_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dxy_a+3*I_TWOBODYOVERLAP_D2y_Dxy;
  abcd[112] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxy_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dxy_a+2*I_TWOBODYOVERLAP_Dyz_Dxy;
  abcd[113] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Dxy_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Dxy_a+1*I_TWOBODYOVERLAP_D2z_Dxy;
  abcd[114] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dxy_a;
  abcd[115] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Dxy_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Dxy_a;
  abcd[116] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Dxy_a;
  abcd[117] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxy_a;
  abcd[118] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dxy_a;
  abcd[119] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dxy_aa;
  abcd[120] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Dxz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Dxz_a;
  abcd[121] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dxz_a+3*1*I_TWOBODYOVERLAP_D2x_Dxz;
  abcd[122] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Dxz_a;
  abcd[123] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Dxz_a+2*2*I_TWOBODYOVERLAP_Dxy_Dxz;
  abcd[124] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxz_a+2*1*I_TWOBODYOVERLAP_Dxz_Dxz;
  abcd[125] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxz_a;
  abcd[126] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dxz_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dxz_a+3*I_TWOBODYOVERLAP_D2y_Dxz;
  abcd[127] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxz_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dxz_a+2*I_TWOBODYOVERLAP_Dyz_Dxz;
  abcd[128] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Dxz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Dxz_a+1*I_TWOBODYOVERLAP_D2z_Dxz;
  abcd[129] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dxz_a;
  abcd[130] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Dxz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Dxz_a;
  abcd[131] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Dxz_a;
  abcd[132] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxz_a;
  abcd[133] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dxz_a;
  abcd[134] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dxz_aa;
  abcd[135] = 4.0E0*I_TWOBODYOVERLAP_I5xy_D2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_D2y_a;
  abcd[136] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2y_a+3*1*I_TWOBODYOVERLAP_D2x_D2y;
  abcd[137] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_D2y_a;
  abcd[138] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_D2y_a+2*2*I_TWOBODYOVERLAP_Dxy_D2y;
  abcd[139] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2y_a+2*1*I_TWOBODYOVERLAP_Dxz_D2y;
  abcd[140] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2y_a;
  abcd[141] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2y_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2y_a+3*I_TWOBODYOVERLAP_D2y_D2y;
  abcd[142] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2y_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2y_a+2*I_TWOBODYOVERLAP_Dyz_D2y;
  abcd[143] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_D2y_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_D2y_a+1*I_TWOBODYOVERLAP_D2z_D2y;
  abcd[144] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2y_a;
  abcd[145] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_D2y_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_D2y_a;
  abcd[146] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_D2y_a;
  abcd[147] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2y_a;
  abcd[148] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2y_a;
  abcd[149] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2y_aa;
  abcd[150] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Dyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Dyz_a;
  abcd[151] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dyz_a+3*1*I_TWOBODYOVERLAP_D2x_Dyz;
  abcd[152] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Dyz_a;
  abcd[153] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Dyz_a+2*2*I_TWOBODYOVERLAP_Dxy_Dyz;
  abcd[154] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dyz_a+2*1*I_TWOBODYOVERLAP_Dxz_Dyz;
  abcd[155] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dyz_a;
  abcd[156] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dyz_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dyz_a+3*I_TWOBODYOVERLAP_D2y_Dyz;
  abcd[157] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dyz_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dyz_a+2*I_TWOBODYOVERLAP_Dyz_Dyz;
  abcd[158] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Dyz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Dyz_a+1*I_TWOBODYOVERLAP_D2z_Dyz;
  abcd[159] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dyz_a;
  abcd[160] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Dyz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Dyz_a;
  abcd[161] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Dyz_a;
  abcd[162] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dyz_a;
  abcd[163] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dyz_a;
  abcd[164] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dyz_aa;
  abcd[165] = 4.0E0*I_TWOBODYOVERLAP_I5xy_D2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_D2z_a;
  abcd[166] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2z_a+3*1*I_TWOBODYOVERLAP_D2x_D2z;
  abcd[167] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_D2z_a;
  abcd[168] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_D2z_a+2*2*I_TWOBODYOVERLAP_Dxy_D2z;
  abcd[169] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2z_a+2*1*I_TWOBODYOVERLAP_Dxz_D2z;
  abcd[170] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2z_a;
  abcd[171] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2z_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2z_a+3*I_TWOBODYOVERLAP_D2y_D2z;
  abcd[172] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2z_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2z_a+2*I_TWOBODYOVERLAP_Dyz_D2z;
  abcd[173] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_D2z_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_D2z_a+1*I_TWOBODYOVERLAP_D2z_D2z;
  abcd[174] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2z_a;
  abcd[175] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_D2z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_D2z_a;
  abcd[176] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_D2z_a;
  abcd[177] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2z_a;
  abcd[178] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2z_a;
  abcd[179] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2z_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
   ************************************************************/
  abcd[180] = 4.0E0*I_TWOBODYOVERLAP_I5xz_D2x_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_D2x_a;
  abcd[181] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_D2x_a;
  abcd[182] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2x_a+3*1*I_TWOBODYOVERLAP_D2x_D2x;
  abcd[183] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2x_a;
  abcd[184] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2x_a+2*1*I_TWOBODYOVERLAP_Dxy_D2x;
  abcd[185] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_D2x_a+2*2*I_TWOBODYOVERLAP_Dxz_D2x;
  abcd[186] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2x_a;
  abcd[187] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_D2x_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_D2x_a+1*I_TWOBODYOVERLAP_D2y_D2x;
  abcd[188] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2x_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2x_a+2*I_TWOBODYOVERLAP_Dyz_D2x;
  abcd[189] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2x_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2x_a+3*I_TWOBODYOVERLAP_D2z_D2x;
  abcd[190] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2x_aa;
  abcd[191] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2x_a;
  abcd[192] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2x_a;
  abcd[193] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_D2x_a;
  abcd[194] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_D2x_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_D2x_a;
  abcd[195] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Dxy_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Dxy_a;
  abcd[196] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Dxy_a;
  abcd[197] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dxy_a+3*1*I_TWOBODYOVERLAP_D2x_Dxy;
  abcd[198] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxy_a;
  abcd[199] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxy_a+2*1*I_TWOBODYOVERLAP_Dxy_Dxy;
  abcd[200] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Dxy_a+2*2*I_TWOBODYOVERLAP_Dxz_Dxy;
  abcd[201] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dxy_a;
  abcd[202] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Dxy_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Dxy_a+1*I_TWOBODYOVERLAP_D2y_Dxy;
  abcd[203] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxy_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dxy_a+2*I_TWOBODYOVERLAP_Dyz_Dxy;
  abcd[204] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dxy_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dxy_a+3*I_TWOBODYOVERLAP_D2z_Dxy;
  abcd[205] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dxy_aa;
  abcd[206] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dxy_a;
  abcd[207] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxy_a;
  abcd[208] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Dxy_a;
  abcd[209] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Dxy_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Dxy_a;
  abcd[210] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Dxz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Dxz_a;
  abcd[211] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Dxz_a;
  abcd[212] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dxz_a+3*1*I_TWOBODYOVERLAP_D2x_Dxz;
  abcd[213] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxz_a;
  abcd[214] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxz_a+2*1*I_TWOBODYOVERLAP_Dxy_Dxz;
  abcd[215] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Dxz_a+2*2*I_TWOBODYOVERLAP_Dxz_Dxz;
  abcd[216] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dxz_a;
  abcd[217] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Dxz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Dxz_a+1*I_TWOBODYOVERLAP_D2y_Dxz;
  abcd[218] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxz_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dxz_a+2*I_TWOBODYOVERLAP_Dyz_Dxz;
  abcd[219] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dxz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dxz_a+3*I_TWOBODYOVERLAP_D2z_Dxz;
  abcd[220] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dxz_aa;
  abcd[221] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dxz_a;
  abcd[222] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxz_a;
  abcd[223] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Dxz_a;
  abcd[224] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Dxz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Dxz_a;
  abcd[225] = 4.0E0*I_TWOBODYOVERLAP_I5xz_D2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_D2y_a;
  abcd[226] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_D2y_a;
  abcd[227] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2y_a+3*1*I_TWOBODYOVERLAP_D2x_D2y;
  abcd[228] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2y_a;
  abcd[229] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2y_a+2*1*I_TWOBODYOVERLAP_Dxy_D2y;
  abcd[230] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_D2y_a+2*2*I_TWOBODYOVERLAP_Dxz_D2y;
  abcd[231] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2y_a;
  abcd[232] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_D2y_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_D2y_a+1*I_TWOBODYOVERLAP_D2y_D2y;
  abcd[233] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2y_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2y_a+2*I_TWOBODYOVERLAP_Dyz_D2y;
  abcd[234] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2y_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2y_a+3*I_TWOBODYOVERLAP_D2z_D2y;
  abcd[235] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2y_aa;
  abcd[236] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2y_a;
  abcd[237] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2y_a;
  abcd[238] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_D2y_a;
  abcd[239] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_D2y_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_D2y_a;
  abcd[240] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Dyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Dyz_a;
  abcd[241] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Dyz_a;
  abcd[242] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dyz_a+3*1*I_TWOBODYOVERLAP_D2x_Dyz;
  abcd[243] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dyz_a;
  abcd[244] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dyz_a+2*1*I_TWOBODYOVERLAP_Dxy_Dyz;
  abcd[245] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Dyz_a+2*2*I_TWOBODYOVERLAP_Dxz_Dyz;
  abcd[246] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dyz_a;
  abcd[247] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Dyz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Dyz_a+1*I_TWOBODYOVERLAP_D2y_Dyz;
  abcd[248] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dyz_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dyz_a+2*I_TWOBODYOVERLAP_Dyz_Dyz;
  abcd[249] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dyz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dyz_a+3*I_TWOBODYOVERLAP_D2z_Dyz;
  abcd[250] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dyz_aa;
  abcd[251] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dyz_a;
  abcd[252] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dyz_a;
  abcd[253] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Dyz_a;
  abcd[254] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Dyz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Dyz_a;
  abcd[255] = 4.0E0*I_TWOBODYOVERLAP_I5xz_D2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_D2z_a;
  abcd[256] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_D2z_a;
  abcd[257] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2z_a+3*1*I_TWOBODYOVERLAP_D2x_D2z;
  abcd[258] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2z_a;
  abcd[259] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2z_a+2*1*I_TWOBODYOVERLAP_Dxy_D2z;
  abcd[260] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_D2z_a+2*2*I_TWOBODYOVERLAP_Dxz_D2z;
  abcd[261] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2z_a;
  abcd[262] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_D2z_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_D2z_a+1*I_TWOBODYOVERLAP_D2y_D2z;
  abcd[263] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2z_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2z_a+2*I_TWOBODYOVERLAP_Dyz_D2z;
  abcd[264] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2z_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2z_a+3*I_TWOBODYOVERLAP_D2z_D2z;
  abcd[265] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2z_aa;
  abcd[266] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2z_a;
  abcd[267] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2z_a;
  abcd[268] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_D2z_a;
  abcd[269] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_D2z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
   ************************************************************/
  abcd[270] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2x_a;
  abcd[271] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_D2x_a;
  abcd[272] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2x_a;
  abcd[273] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2x_a+2*1*I_TWOBODYOVERLAP_D2x_D2x;
  abcd[274] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2x_a;
  abcd[275] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_D2x_a;
  abcd[276] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_D2x_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_D2x_a+3*2*I_TWOBODYOVERLAP_Dxy_D2x;
  abcd[277] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_D2x_a+2*1*I_TWOBODYOVERLAP_Dxz_D2x;
  abcd[278] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2x_a;
  abcd[279] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2x_a;
  abcd[280] = 4.0E0*I_TWOBODYOVERLAP_I6y_D2x_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_D2x_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_D2x_a+4*3*I_TWOBODYOVERLAP_D2y_D2x;
  abcd[281] = 4.0E0*I_TWOBODYOVERLAP_I5yz_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_D2x_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_D2x_a+3*2*I_TWOBODYOVERLAP_Dyz_D2x;
  abcd[282] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2x_a+2*1*I_TWOBODYOVERLAP_D2z_D2x;
  abcd[283] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_D2x_a;
  abcd[284] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2x_a;
  abcd[285] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dxy_a;
  abcd[286] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Dxy_a;
  abcd[287] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dxy_a;
  abcd[288] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dxy_a+2*1*I_TWOBODYOVERLAP_D2x_Dxy;
  abcd[289] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxy_a;
  abcd[290] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Dxy_a;
  abcd[291] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_Dxy_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Dxy_a+3*2*I_TWOBODYOVERLAP_Dxy_Dxy;
  abcd[292] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Dxy_a+2*1*I_TWOBODYOVERLAP_Dxz_Dxy;
  abcd[293] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxy_a;
  abcd[294] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dxy_a;
  abcd[295] = 4.0E0*I_TWOBODYOVERLAP_I6y_Dxy_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_Dxy_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_Dxy_a+4*3*I_TWOBODYOVERLAP_D2y_Dxy;
  abcd[296] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_Dxy_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Dxy_a+3*2*I_TWOBODYOVERLAP_Dyz_Dxy;
  abcd[297] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dxy_a+2*1*I_TWOBODYOVERLAP_D2z_Dxy;
  abcd[298] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Dxy_a;
  abcd[299] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dxy_a;
  abcd[300] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dxz_a;
  abcd[301] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Dxz_a;
  abcd[302] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dxz_a;
  abcd[303] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dxz_a+2*1*I_TWOBODYOVERLAP_D2x_Dxz;
  abcd[304] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxz_a;
  abcd[305] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Dxz_a;
  abcd[306] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_Dxz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Dxz_a+3*2*I_TWOBODYOVERLAP_Dxy_Dxz;
  abcd[307] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Dxz_a+2*1*I_TWOBODYOVERLAP_Dxz_Dxz;
  abcd[308] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxz_a;
  abcd[309] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dxz_a;
  abcd[310] = 4.0E0*I_TWOBODYOVERLAP_I6y_Dxz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_Dxz_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_Dxz_a+4*3*I_TWOBODYOVERLAP_D2y_Dxz;
  abcd[311] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_Dxz_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Dxz_a+3*2*I_TWOBODYOVERLAP_Dyz_Dxz;
  abcd[312] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dxz_a+2*1*I_TWOBODYOVERLAP_D2z_Dxz;
  abcd[313] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Dxz_a;
  abcd[314] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dxz_a;
  abcd[315] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2y_a;
  abcd[316] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_D2y_a;
  abcd[317] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2y_a;
  abcd[318] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2y_a+2*1*I_TWOBODYOVERLAP_D2x_D2y;
  abcd[319] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2y_a;
  abcd[320] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_D2y_a;
  abcd[321] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_D2y_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_D2y_a+3*2*I_TWOBODYOVERLAP_Dxy_D2y;
  abcd[322] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_D2y_a+2*1*I_TWOBODYOVERLAP_Dxz_D2y;
  abcd[323] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2y_a;
  abcd[324] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2y_a;
  abcd[325] = 4.0E0*I_TWOBODYOVERLAP_I6y_D2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_D2y_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_D2y_a+4*3*I_TWOBODYOVERLAP_D2y_D2y;
  abcd[326] = 4.0E0*I_TWOBODYOVERLAP_I5yz_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_D2y_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_D2y_a+3*2*I_TWOBODYOVERLAP_Dyz_D2y;
  abcd[327] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2y_a+2*1*I_TWOBODYOVERLAP_D2z_D2y;
  abcd[328] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_D2y_a;
  abcd[329] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2y_a;
  abcd[330] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dyz_a;
  abcd[331] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Dyz_a;
  abcd[332] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dyz_a;
  abcd[333] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Dyz_a+2*1*I_TWOBODYOVERLAP_D2x_Dyz;
  abcd[334] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dyz_a;
  abcd[335] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Dyz_a;
  abcd[336] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_Dyz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Dyz_a+3*2*I_TWOBODYOVERLAP_Dxy_Dyz;
  abcd[337] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Dyz_a+2*1*I_TWOBODYOVERLAP_Dxz_Dyz;
  abcd[338] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dyz_a;
  abcd[339] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dyz_a;
  abcd[340] = 4.0E0*I_TWOBODYOVERLAP_I6y_Dyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_Dyz_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_Dyz_a+4*3*I_TWOBODYOVERLAP_D2y_Dyz;
  abcd[341] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_Dyz_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Dyz_a+3*2*I_TWOBODYOVERLAP_Dyz_Dyz;
  abcd[342] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dyz_a+2*1*I_TWOBODYOVERLAP_D2z_Dyz;
  abcd[343] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Dyz_a;
  abcd[344] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dyz_a;
  abcd[345] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2z_a;
  abcd[346] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_D2z_a;
  abcd[347] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2z_a;
  abcd[348] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_D2z_a+2*1*I_TWOBODYOVERLAP_D2x_D2z;
  abcd[349] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2z_a;
  abcd[350] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_D2z_a;
  abcd[351] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_D2z_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_D2z_a+3*2*I_TWOBODYOVERLAP_Dxy_D2z;
  abcd[352] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_D2z_a+2*1*I_TWOBODYOVERLAP_Dxz_D2z;
  abcd[353] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2z_a;
  abcd[354] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2z_a;
  abcd[355] = 4.0E0*I_TWOBODYOVERLAP_I6y_D2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_D2z_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_D2z_a+4*3*I_TWOBODYOVERLAP_D2y_D2z;
  abcd[356] = 4.0E0*I_TWOBODYOVERLAP_I5yz_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_D2z_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_D2z_a+3*2*I_TWOBODYOVERLAP_Dyz_D2z;
  abcd[357] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2z_a+2*1*I_TWOBODYOVERLAP_D2z_D2z;
  abcd[358] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_D2z_a;
  abcd[359] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
   ************************************************************/
  abcd[360] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2x_aa;
  abcd[361] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2x_a;
  abcd[362] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2x_a;
  abcd[363] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2x_a;
  abcd[364] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_D2x_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_D2x_a+1*I_TWOBODYOVERLAP_D2x_D2x;
  abcd[365] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2x_a;
  abcd[366] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_D2x_a;
  abcd[367] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2x_a+2*1*I_TWOBODYOVERLAP_Dxy_D2x;
  abcd[368] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2x_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2x_a+2*I_TWOBODYOVERLAP_Dxz_D2x;
  abcd[369] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_D2x_a;
  abcd[370] = 4.0E0*I_TWOBODYOVERLAP_I5yz_D2x_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_D2x_a;
  abcd[371] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2x_a+3*1*I_TWOBODYOVERLAP_D2y_D2x;
  abcd[372] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_D2x_a+2*2*I_TWOBODYOVERLAP_Dyz_D2x;
  abcd[373] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2x_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2x_a+3*I_TWOBODYOVERLAP_D2z_D2x;
  abcd[374] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_D2x_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_D2x_a;
  abcd[375] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dxy_aa;
  abcd[376] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dxy_a;
  abcd[377] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dxy_a;
  abcd[378] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxy_a;
  abcd[379] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Dxy_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Dxy_a+1*I_TWOBODYOVERLAP_D2x_Dxy;
  abcd[380] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxy_a;
  abcd[381] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Dxy_a;
  abcd[382] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxy_a+2*1*I_TWOBODYOVERLAP_Dxy_Dxy;
  abcd[383] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxy_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dxy_a+2*I_TWOBODYOVERLAP_Dxz_Dxy;
  abcd[384] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Dxy_a;
  abcd[385] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Dxy_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Dxy_a;
  abcd[386] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dxy_a+3*1*I_TWOBODYOVERLAP_D2y_Dxy;
  abcd[387] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Dxy_a+2*2*I_TWOBODYOVERLAP_Dyz_Dxy;
  abcd[388] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dxy_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dxy_a+3*I_TWOBODYOVERLAP_D2z_Dxy;
  abcd[389] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Dxy_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Dxy_a;
  abcd[390] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dxz_aa;
  abcd[391] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dxz_a;
  abcd[392] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dxz_a;
  abcd[393] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxz_a;
  abcd[394] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Dxz_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Dxz_a+1*I_TWOBODYOVERLAP_D2x_Dxz;
  abcd[395] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxz_a;
  abcd[396] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Dxz_a;
  abcd[397] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxz_a+2*1*I_TWOBODYOVERLAP_Dxy_Dxz;
  abcd[398] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxz_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dxz_a+2*I_TWOBODYOVERLAP_Dxz_Dxz;
  abcd[399] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Dxz_a;
  abcd[400] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Dxz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Dxz_a;
  abcd[401] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dxz_a+3*1*I_TWOBODYOVERLAP_D2y_Dxz;
  abcd[402] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Dxz_a+2*2*I_TWOBODYOVERLAP_Dyz_Dxz;
  abcd[403] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dxz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dxz_a+3*I_TWOBODYOVERLAP_D2z_Dxz;
  abcd[404] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Dxz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Dxz_a;
  abcd[405] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2y_aa;
  abcd[406] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2y_a;
  abcd[407] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2y_a;
  abcd[408] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2y_a;
  abcd[409] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_D2y_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_D2y_a+1*I_TWOBODYOVERLAP_D2x_D2y;
  abcd[410] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2y_a;
  abcd[411] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_D2y_a;
  abcd[412] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2y_a+2*1*I_TWOBODYOVERLAP_Dxy_D2y;
  abcd[413] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2y_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2y_a+2*I_TWOBODYOVERLAP_Dxz_D2y;
  abcd[414] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_D2y_a;
  abcd[415] = 4.0E0*I_TWOBODYOVERLAP_I5yz_D2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_D2y_a;
  abcd[416] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2y_a+3*1*I_TWOBODYOVERLAP_D2y_D2y;
  abcd[417] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_D2y_a+2*2*I_TWOBODYOVERLAP_Dyz_D2y;
  abcd[418] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2y_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2y_a+3*I_TWOBODYOVERLAP_D2z_D2y;
  abcd[419] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_D2y_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_D2y_a;
  abcd[420] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Dyz_aa;
  abcd[421] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dyz_a;
  abcd[422] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dyz_a;
  abcd[423] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dyz_a;
  abcd[424] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Dyz_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Dyz_a+1*I_TWOBODYOVERLAP_D2x_Dyz;
  abcd[425] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dyz_a;
  abcd[426] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Dyz_a;
  abcd[427] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dyz_a+2*1*I_TWOBODYOVERLAP_Dxy_Dyz;
  abcd[428] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dyz_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Dyz_a+2*I_TWOBODYOVERLAP_Dxz_Dyz;
  abcd[429] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Dyz_a;
  abcd[430] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Dyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Dyz_a;
  abcd[431] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dyz_a+3*1*I_TWOBODYOVERLAP_D2y_Dyz;
  abcd[432] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Dyz_a+2*2*I_TWOBODYOVERLAP_Dyz_Dyz;
  abcd[433] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dyz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Dyz_a+3*I_TWOBODYOVERLAP_D2z_Dyz;
  abcd[434] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Dyz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Dyz_a;
  abcd[435] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_D2z_aa;
  abcd[436] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2z_a;
  abcd[437] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2z_a;
  abcd[438] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2z_a;
  abcd[439] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_D2z_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_D2z_a+1*I_TWOBODYOVERLAP_D2x_D2z;
  abcd[440] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2z_a;
  abcd[441] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_D2z_a;
  abcd[442] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2z_a+2*1*I_TWOBODYOVERLAP_Dxy_D2z;
  abcd[443] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2z_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_D2z_a+2*I_TWOBODYOVERLAP_Dxz_D2z;
  abcd[444] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_D2z_a;
  abcd[445] = 4.0E0*I_TWOBODYOVERLAP_I5yz_D2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_D2z_a;
  abcd[446] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2z_a+3*1*I_TWOBODYOVERLAP_D2y_D2z;
  abcd[447] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_D2z_a+2*2*I_TWOBODYOVERLAP_Dyz_D2z;
  abcd[448] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2z_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_D2z_a+3*I_TWOBODYOVERLAP_D2z_D2z;
  abcd[449] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_D2z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
   ************************************************************/
  abcd[450] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2x_a;
  abcd[451] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2x_a;
  abcd[452] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_D2x_a;
  abcd[453] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_D2x_a;
  abcd[454] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2x_a;
  abcd[455] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2x_a+2*1*I_TWOBODYOVERLAP_D2x_D2x;
  abcd[456] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2x_a;
  abcd[457] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2x_a;
  abcd[458] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_D2x_a+2*1*I_TWOBODYOVERLAP_Dxy_D2x;
  abcd[459] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_D2x_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_D2x_a+3*2*I_TWOBODYOVERLAP_Dxz_D2x;
  abcd[460] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2x_a;
  abcd[461] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_D2x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2x_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_D2x_a;
  abcd[462] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_D2x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_D2x_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2x_a+2*1*I_TWOBODYOVERLAP_D2y_D2x;
  abcd[463] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_D2x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_D2x_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_D2x_a+3*2*I_TWOBODYOVERLAP_Dyz_D2x;
  abcd[464] = 4.0E0*I_TWOBODYOVERLAP_I6z_D2x_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_D2x_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_D2x_a+4*3*I_TWOBODYOVERLAP_D2z_D2x;
  abcd[465] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dxy_a;
  abcd[466] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dxy_a;
  abcd[467] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Dxy_a;
  abcd[468] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Dxy_a;
  abcd[469] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxy_a;
  abcd[470] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dxy_a+2*1*I_TWOBODYOVERLAP_D2x_Dxy;
  abcd[471] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dxy_a;
  abcd[472] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxy_a;
  abcd[473] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Dxy_a+2*1*I_TWOBODYOVERLAP_Dxy_Dxy;
  abcd[474] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_Dxy_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Dxy_a+3*2*I_TWOBODYOVERLAP_Dxz_Dxy;
  abcd[475] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dxy_a;
  abcd[476] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Dxy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dxy_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Dxy_a;
  abcd[477] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Dxy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Dxy_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dxy_a+2*1*I_TWOBODYOVERLAP_D2y_Dxy;
  abcd[478] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Dxy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_Dxy_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Dxy_a+3*2*I_TWOBODYOVERLAP_Dyz_Dxy;
  abcd[479] = 4.0E0*I_TWOBODYOVERLAP_I6z_Dxy_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_Dxy_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_Dxy_a+4*3*I_TWOBODYOVERLAP_D2z_Dxy;
  abcd[480] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dxz_a;
  abcd[481] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dxz_a;
  abcd[482] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Dxz_a;
  abcd[483] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Dxz_a;
  abcd[484] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dxz_a;
  abcd[485] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dxz_a+2*1*I_TWOBODYOVERLAP_D2x_Dxz;
  abcd[486] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dxz_a;
  abcd[487] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dxz_a;
  abcd[488] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Dxz_a+2*1*I_TWOBODYOVERLAP_Dxy_Dxz;
  abcd[489] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_Dxz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Dxz_a+3*2*I_TWOBODYOVERLAP_Dxz_Dxz;
  abcd[490] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dxz_a;
  abcd[491] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Dxz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dxz_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Dxz_a;
  abcd[492] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Dxz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Dxz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dxz_a+2*1*I_TWOBODYOVERLAP_D2y_Dxz;
  abcd[493] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Dxz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_Dxz_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Dxz_a+3*2*I_TWOBODYOVERLAP_Dyz_Dxz;
  abcd[494] = 4.0E0*I_TWOBODYOVERLAP_I6z_Dxz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_Dxz_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_Dxz_a+4*3*I_TWOBODYOVERLAP_D2z_Dxz;
  abcd[495] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2y_a;
  abcd[496] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2y_a;
  abcd[497] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_D2y_a;
  abcd[498] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_D2y_a;
  abcd[499] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2y_a;
  abcd[500] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2y_a+2*1*I_TWOBODYOVERLAP_D2x_D2y;
  abcd[501] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2y_a;
  abcd[502] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2y_a;
  abcd[503] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_D2y_a+2*1*I_TWOBODYOVERLAP_Dxy_D2y;
  abcd[504] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_D2y_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_D2y_a+3*2*I_TWOBODYOVERLAP_Dxz_D2y;
  abcd[505] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2y_a;
  abcd[506] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_D2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2y_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_D2y_a;
  abcd[507] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_D2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_D2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2y_a+2*1*I_TWOBODYOVERLAP_D2y_D2y;
  abcd[508] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_D2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_D2y_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_D2y_a+3*2*I_TWOBODYOVERLAP_Dyz_D2y;
  abcd[509] = 4.0E0*I_TWOBODYOVERLAP_I6z_D2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_D2y_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_D2y_a+4*3*I_TWOBODYOVERLAP_D2z_D2y;
  abcd[510] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Dyz_a;
  abcd[511] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Dyz_a;
  abcd[512] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Dyz_a;
  abcd[513] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Dyz_a;
  abcd[514] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Dyz_a;
  abcd[515] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Dyz_a+2*1*I_TWOBODYOVERLAP_D2x_Dyz;
  abcd[516] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Dyz_a;
  abcd[517] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Dyz_a;
  abcd[518] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Dyz_a+2*1*I_TWOBODYOVERLAP_Dxy_Dyz;
  abcd[519] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_Dyz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Dyz_a+3*2*I_TWOBODYOVERLAP_Dxz_Dyz;
  abcd[520] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Dyz_a;
  abcd[521] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Dyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Dyz_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Dyz_a;
  abcd[522] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Dyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Dyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Dyz_a+2*1*I_TWOBODYOVERLAP_D2y_Dyz;
  abcd[523] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Dyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_Dyz_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Dyz_a+3*2*I_TWOBODYOVERLAP_Dyz_Dyz;
  abcd[524] = 4.0E0*I_TWOBODYOVERLAP_I6z_Dyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_Dyz_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_Dyz_a+4*3*I_TWOBODYOVERLAP_D2z_Dyz;
  abcd[525] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_D2z_a;
  abcd[526] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_D2z_a;
  abcd[527] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_D2z_a;
  abcd[528] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_D2z_a;
  abcd[529] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_D2z_a;
  abcd[530] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_D2z_a+2*1*I_TWOBODYOVERLAP_D2x_D2z;
  abcd[531] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_D2z_a;
  abcd[532] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_D2z_a;
  abcd[533] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_D2z_a+2*1*I_TWOBODYOVERLAP_Dxy_D2z;
  abcd[534] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_D2z_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_D2z_a+3*2*I_TWOBODYOVERLAP_Dxz_D2z;
  abcd[535] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_D2z_a;
  abcd[536] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_D2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_D2z_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_D2z_a;
  abcd[537] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_D2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_D2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_D2z_a+2*1*I_TWOBODYOVERLAP_D2y_D2z;
  abcd[538] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_D2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_D2z_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_D2z_a+3*2*I_TWOBODYOVERLAP_Dyz_D2z;
  abcd[539] = 4.0E0*I_TWOBODYOVERLAP_I6z_D2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_D2z_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_D2z_a+4*3*I_TWOBODYOVERLAP_D2z_D2z;
}
