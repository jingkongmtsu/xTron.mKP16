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
// BRA1 as redundant position, total RHS integrals evaluated as: 1454
// BRA2 as redundant position, total RHS integrals evaluated as: 1188
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

void hgp_os_twobodyoverlap_h_sp_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_I6x_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S_C5_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4x_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S_C5 = 0.0E0;
  Double I_TWOBODYOVERLAP_K7x_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xy_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K5xyz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x2yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4xy2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x3yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x2y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3xy3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x4yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x3y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x2y3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2xy4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx5yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx4y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx3y3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx2y4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kxy5z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K7y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K6yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K5y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4y3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3y4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2y5z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ky6z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K7z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6x_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S_C1005_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5x_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G4x_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S_C1005 = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double ic2_1 = icoe[ip2+1*inp2];
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
    Double I_TWOBODYOVERLAP_S_S_vrr = fbra;
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
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_S_vrr = PAZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_S_vrr = PAX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_S_vrr = PAY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_S_vrr = PAX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_S_vrr = PAX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_S_vrr = PAY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2y_S_vrr;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S_C5_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_I_S_C5_a_coefs = ic2*alpha;
    I_TWOBODYOVERLAP_I6x_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S_C5_a += SQ_TWOBODYOVERLAP_I_S_C5_a_coefs*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S_C5
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_G_S_C5_coefs = ic2;
    I_TWOBODYOVERLAP_G4x_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S_C5 += SQ_TWOBODYOVERLAP_G_S_C5_coefs*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S_C1005_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs = ic2_1*alpha;
    I_TWOBODYOVERLAP_K7x_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K7x_S_vrr;
    I_TWOBODYOVERLAP_K6xy_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K6xy_S_vrr;
    I_TWOBODYOVERLAP_K6xz_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K6xz_S_vrr;
    I_TWOBODYOVERLAP_K5x2y_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    I_TWOBODYOVERLAP_K5xyz_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K5xyz_S_vrr;
    I_TWOBODYOVERLAP_K5x2z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3y_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    I_TWOBODYOVERLAP_K4x2yz_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K4x2yz_S_vrr;
    I_TWOBODYOVERLAP_K4xy2z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K4xy2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4y_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    I_TWOBODYOVERLAP_K3x3yz_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    I_TWOBODYOVERLAP_K3x2y2z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K3x2y2z_S_vrr;
    I_TWOBODYOVERLAP_K3xy3z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K3xy3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5y_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    I_TWOBODYOVERLAP_K2x4yz_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K2x4yz_S_vrr;
    I_TWOBODYOVERLAP_K2x3y2z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K2x3y2z_S_vrr;
    I_TWOBODYOVERLAP_K2x2y3z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K2x2y3z_S_vrr;
    I_TWOBODYOVERLAP_K2xy4z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K2xy4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6y_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    I_TWOBODYOVERLAP_Kx5yz_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_Kx5yz_S_vrr;
    I_TWOBODYOVERLAP_Kx4y2z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_Kx4y2z_S_vrr;
    I_TWOBODYOVERLAP_Kx3y3z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_Kx3y3z_S_vrr;
    I_TWOBODYOVERLAP_Kx2y4z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_Kx2y4z_S_vrr;
    I_TWOBODYOVERLAP_Kxy5z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_Kxy5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    I_TWOBODYOVERLAP_K7y_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K7y_S_vrr;
    I_TWOBODYOVERLAP_K6yz_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K6yz_S_vrr;
    I_TWOBODYOVERLAP_K5y2z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    I_TWOBODYOVERLAP_K4y3z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    I_TWOBODYOVERLAP_K3y4z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    I_TWOBODYOVERLAP_K2y5z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    I_TWOBODYOVERLAP_Ky6z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    I_TWOBODYOVERLAP_K7z_S_C1005_a += SQ_TWOBODYOVERLAP_K_S_C1005_a_coefs*I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S_C1005_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs = ic2_1*alpha;
    I_TWOBODYOVERLAP_I6x_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S_C1005_a += SQ_TWOBODYOVERLAP_I_S_C1005_a_coefs*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C1005
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_H_S_C1005_coefs = ic2_1;
    I_TWOBODYOVERLAP_H5x_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S_C1005
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_G_S_C1005_coefs = ic2_1;
    I_TWOBODYOVERLAP_G4x_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S_C1005 += SQ_TWOBODYOVERLAP_G_S_C1005_coefs*I_TWOBODYOVERLAP_G4z_S_vrr;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_G_P_C1005
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C1005
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_C1005
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_Px_C1005 = I_TWOBODYOVERLAP_H5x_S_C1005+ABX*I_TWOBODYOVERLAP_G4x_S_C1005;
  Double I_TWOBODYOVERLAP_G3xy_Px_C1005 = I_TWOBODYOVERLAP_H4xy_S_C1005+ABX*I_TWOBODYOVERLAP_G3xy_S_C1005;
  Double I_TWOBODYOVERLAP_G3xz_Px_C1005 = I_TWOBODYOVERLAP_H4xz_S_C1005+ABX*I_TWOBODYOVERLAP_G3xz_S_C1005;
  Double I_TWOBODYOVERLAP_G2x2y_Px_C1005 = I_TWOBODYOVERLAP_H3x2y_S_C1005+ABX*I_TWOBODYOVERLAP_G2x2y_S_C1005;
  Double I_TWOBODYOVERLAP_G2xyz_Px_C1005 = I_TWOBODYOVERLAP_H3xyz_S_C1005+ABX*I_TWOBODYOVERLAP_G2xyz_S_C1005;
  Double I_TWOBODYOVERLAP_G2x2z_Px_C1005 = I_TWOBODYOVERLAP_H3x2z_S_C1005+ABX*I_TWOBODYOVERLAP_G2x2z_S_C1005;
  Double I_TWOBODYOVERLAP_Gx3y_Px_C1005 = I_TWOBODYOVERLAP_H2x3y_S_C1005+ABX*I_TWOBODYOVERLAP_Gx3y_S_C1005;
  Double I_TWOBODYOVERLAP_Gx2yz_Px_C1005 = I_TWOBODYOVERLAP_H2x2yz_S_C1005+ABX*I_TWOBODYOVERLAP_Gx2yz_S_C1005;
  Double I_TWOBODYOVERLAP_Gxy2z_Px_C1005 = I_TWOBODYOVERLAP_H2xy2z_S_C1005+ABX*I_TWOBODYOVERLAP_Gxy2z_S_C1005;
  Double I_TWOBODYOVERLAP_Gx3z_Px_C1005 = I_TWOBODYOVERLAP_H2x3z_S_C1005+ABX*I_TWOBODYOVERLAP_Gx3z_S_C1005;
  Double I_TWOBODYOVERLAP_G4y_Px_C1005 = I_TWOBODYOVERLAP_Hx4y_S_C1005+ABX*I_TWOBODYOVERLAP_G4y_S_C1005;
  Double I_TWOBODYOVERLAP_G3yz_Px_C1005 = I_TWOBODYOVERLAP_Hx3yz_S_C1005+ABX*I_TWOBODYOVERLAP_G3yz_S_C1005;
  Double I_TWOBODYOVERLAP_G2y2z_Px_C1005 = I_TWOBODYOVERLAP_Hx2y2z_S_C1005+ABX*I_TWOBODYOVERLAP_G2y2z_S_C1005;
  Double I_TWOBODYOVERLAP_Gy3z_Px_C1005 = I_TWOBODYOVERLAP_Hxy3z_S_C1005+ABX*I_TWOBODYOVERLAP_Gy3z_S_C1005;
  Double I_TWOBODYOVERLAP_G4z_Px_C1005 = I_TWOBODYOVERLAP_Hx4z_S_C1005+ABX*I_TWOBODYOVERLAP_G4z_S_C1005;
  Double I_TWOBODYOVERLAP_G4x_Py_C1005 = I_TWOBODYOVERLAP_H4xy_S_C1005+ABY*I_TWOBODYOVERLAP_G4x_S_C1005;
  Double I_TWOBODYOVERLAP_G3xy_Py_C1005 = I_TWOBODYOVERLAP_H3x2y_S_C1005+ABY*I_TWOBODYOVERLAP_G3xy_S_C1005;
  Double I_TWOBODYOVERLAP_G3xz_Py_C1005 = I_TWOBODYOVERLAP_H3xyz_S_C1005+ABY*I_TWOBODYOVERLAP_G3xz_S_C1005;
  Double I_TWOBODYOVERLAP_G2x2y_Py_C1005 = I_TWOBODYOVERLAP_H2x3y_S_C1005+ABY*I_TWOBODYOVERLAP_G2x2y_S_C1005;
  Double I_TWOBODYOVERLAP_G2xyz_Py_C1005 = I_TWOBODYOVERLAP_H2x2yz_S_C1005+ABY*I_TWOBODYOVERLAP_G2xyz_S_C1005;
  Double I_TWOBODYOVERLAP_G2x2z_Py_C1005 = I_TWOBODYOVERLAP_H2xy2z_S_C1005+ABY*I_TWOBODYOVERLAP_G2x2z_S_C1005;
  Double I_TWOBODYOVERLAP_Gx3y_Py_C1005 = I_TWOBODYOVERLAP_Hx4y_S_C1005+ABY*I_TWOBODYOVERLAP_Gx3y_S_C1005;
  Double I_TWOBODYOVERLAP_Gx2yz_Py_C1005 = I_TWOBODYOVERLAP_Hx3yz_S_C1005+ABY*I_TWOBODYOVERLAP_Gx2yz_S_C1005;
  Double I_TWOBODYOVERLAP_Gxy2z_Py_C1005 = I_TWOBODYOVERLAP_Hx2y2z_S_C1005+ABY*I_TWOBODYOVERLAP_Gxy2z_S_C1005;
  Double I_TWOBODYOVERLAP_Gx3z_Py_C1005 = I_TWOBODYOVERLAP_Hxy3z_S_C1005+ABY*I_TWOBODYOVERLAP_Gx3z_S_C1005;
  Double I_TWOBODYOVERLAP_G4y_Py_C1005 = I_TWOBODYOVERLAP_H5y_S_C1005+ABY*I_TWOBODYOVERLAP_G4y_S_C1005;
  Double I_TWOBODYOVERLAP_G3yz_Py_C1005 = I_TWOBODYOVERLAP_H4yz_S_C1005+ABY*I_TWOBODYOVERLAP_G3yz_S_C1005;
  Double I_TWOBODYOVERLAP_G2y2z_Py_C1005 = I_TWOBODYOVERLAP_H3y2z_S_C1005+ABY*I_TWOBODYOVERLAP_G2y2z_S_C1005;
  Double I_TWOBODYOVERLAP_Gy3z_Py_C1005 = I_TWOBODYOVERLAP_H2y3z_S_C1005+ABY*I_TWOBODYOVERLAP_Gy3z_S_C1005;
  Double I_TWOBODYOVERLAP_G4z_Py_C1005 = I_TWOBODYOVERLAP_Hy4z_S_C1005+ABY*I_TWOBODYOVERLAP_G4z_S_C1005;
  Double I_TWOBODYOVERLAP_G4x_Pz_C1005 = I_TWOBODYOVERLAP_H4xz_S_C1005+ABZ*I_TWOBODYOVERLAP_G4x_S_C1005;
  Double I_TWOBODYOVERLAP_G3xy_Pz_C1005 = I_TWOBODYOVERLAP_H3xyz_S_C1005+ABZ*I_TWOBODYOVERLAP_G3xy_S_C1005;
  Double I_TWOBODYOVERLAP_G3xz_Pz_C1005 = I_TWOBODYOVERLAP_H3x2z_S_C1005+ABZ*I_TWOBODYOVERLAP_G3xz_S_C1005;
  Double I_TWOBODYOVERLAP_G2x2y_Pz_C1005 = I_TWOBODYOVERLAP_H2x2yz_S_C1005+ABZ*I_TWOBODYOVERLAP_G2x2y_S_C1005;
  Double I_TWOBODYOVERLAP_G2xyz_Pz_C1005 = I_TWOBODYOVERLAP_H2xy2z_S_C1005+ABZ*I_TWOBODYOVERLAP_G2xyz_S_C1005;
  Double I_TWOBODYOVERLAP_G2x2z_Pz_C1005 = I_TWOBODYOVERLAP_H2x3z_S_C1005+ABZ*I_TWOBODYOVERLAP_G2x2z_S_C1005;
  Double I_TWOBODYOVERLAP_Gx3y_Pz_C1005 = I_TWOBODYOVERLAP_Hx3yz_S_C1005+ABZ*I_TWOBODYOVERLAP_Gx3y_S_C1005;
  Double I_TWOBODYOVERLAP_Gx2yz_Pz_C1005 = I_TWOBODYOVERLAP_Hx2y2z_S_C1005+ABZ*I_TWOBODYOVERLAP_Gx2yz_S_C1005;
  Double I_TWOBODYOVERLAP_Gxy2z_Pz_C1005 = I_TWOBODYOVERLAP_Hxy3z_S_C1005+ABZ*I_TWOBODYOVERLAP_Gxy2z_S_C1005;
  Double I_TWOBODYOVERLAP_Gx3z_Pz_C1005 = I_TWOBODYOVERLAP_Hx4z_S_C1005+ABZ*I_TWOBODYOVERLAP_Gx3z_S_C1005;
  Double I_TWOBODYOVERLAP_G4y_Pz_C1005 = I_TWOBODYOVERLAP_H4yz_S_C1005+ABZ*I_TWOBODYOVERLAP_G4y_S_C1005;
  Double I_TWOBODYOVERLAP_G3yz_Pz_C1005 = I_TWOBODYOVERLAP_H3y2z_S_C1005+ABZ*I_TWOBODYOVERLAP_G3yz_S_C1005;
  Double I_TWOBODYOVERLAP_G2y2z_Pz_C1005 = I_TWOBODYOVERLAP_H2y3z_S_C1005+ABZ*I_TWOBODYOVERLAP_G2y2z_S_C1005;
  Double I_TWOBODYOVERLAP_Gy3z_Pz_C1005 = I_TWOBODYOVERLAP_Hy4z_S_C1005+ABZ*I_TWOBODYOVERLAP_Gy3z_S_C1005;
  Double I_TWOBODYOVERLAP_G4z_Pz_C1005 = I_TWOBODYOVERLAP_H5z_S_C1005+ABZ*I_TWOBODYOVERLAP_G4z_S_C1005;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_P_C1005_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_C1005_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_Px_C1005_a = I_TWOBODYOVERLAP_K7x_S_C1005_a+ABX*I_TWOBODYOVERLAP_I6x_S_C1005_a;
  Double I_TWOBODYOVERLAP_I5xy_Px_C1005_a = I_TWOBODYOVERLAP_K6xy_S_C1005_a+ABX*I_TWOBODYOVERLAP_I5xy_S_C1005_a;
  Double I_TWOBODYOVERLAP_I5xz_Px_C1005_a = I_TWOBODYOVERLAP_K6xz_S_C1005_a+ABX*I_TWOBODYOVERLAP_I5xz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4x2y_Px_C1005_a = I_TWOBODYOVERLAP_K5x2y_S_C1005_a+ABX*I_TWOBODYOVERLAP_I4x2y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4xyz_Px_C1005_a = I_TWOBODYOVERLAP_K5xyz_S_C1005_a+ABX*I_TWOBODYOVERLAP_I4xyz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4x2z_Px_C1005_a = I_TWOBODYOVERLAP_K5x2z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I4x2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3x3y_Px_C1005_a = I_TWOBODYOVERLAP_K4x3y_S_C1005_a+ABX*I_TWOBODYOVERLAP_I3x3y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Px_C1005_a = I_TWOBODYOVERLAP_K4x2yz_S_C1005_a+ABX*I_TWOBODYOVERLAP_I3x2yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Px_C1005_a = I_TWOBODYOVERLAP_K4xy2z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I3xy2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3x3z_Px_C1005_a = I_TWOBODYOVERLAP_K4x3z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I3x3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x4y_Px_C1005_a = I_TWOBODYOVERLAP_K3x4y_S_C1005_a+ABX*I_TWOBODYOVERLAP_I2x4y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Px_C1005_a = I_TWOBODYOVERLAP_K3x3yz_S_C1005_a+ABX*I_TWOBODYOVERLAP_I2x3yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Px_C1005_a = I_TWOBODYOVERLAP_K3x2y2z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Px_C1005_a = I_TWOBODYOVERLAP_K3xy3z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I2xy3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x4z_Px_C1005_a = I_TWOBODYOVERLAP_K3x4z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I2x4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix5y_Px_C1005_a = I_TWOBODYOVERLAP_K2x5y_S_C1005_a+ABX*I_TWOBODYOVERLAP_Ix5y_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Px_C1005_a = I_TWOBODYOVERLAP_K2x4yz_S_C1005_a+ABX*I_TWOBODYOVERLAP_Ix4yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Px_C1005_a = I_TWOBODYOVERLAP_K2x3y2z_S_C1005_a+ABX*I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Px_C1005_a = I_TWOBODYOVERLAP_K2x2y3z_S_C1005_a+ABX*I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Px_C1005_a = I_TWOBODYOVERLAP_K2xy4z_S_C1005_a+ABX*I_TWOBODYOVERLAP_Ixy4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix5z_Px_C1005_a = I_TWOBODYOVERLAP_K2x5z_S_C1005_a+ABX*I_TWOBODYOVERLAP_Ix5z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I6y_Px_C1005_a = I_TWOBODYOVERLAP_Kx6y_S_C1005_a+ABX*I_TWOBODYOVERLAP_I6y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I5yz_Px_C1005_a = I_TWOBODYOVERLAP_Kx5yz_S_C1005_a+ABX*I_TWOBODYOVERLAP_I5yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4y2z_Px_C1005_a = I_TWOBODYOVERLAP_Kx4y2z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I4y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3y3z_Px_C1005_a = I_TWOBODYOVERLAP_Kx3y3z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I3y3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2y4z_Px_C1005_a = I_TWOBODYOVERLAP_Kx2y4z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I2y4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Iy5z_Px_C1005_a = I_TWOBODYOVERLAP_Kxy5z_S_C1005_a+ABX*I_TWOBODYOVERLAP_Iy5z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I6z_Px_C1005_a = I_TWOBODYOVERLAP_Kx6z_S_C1005_a+ABX*I_TWOBODYOVERLAP_I6z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I6x_Py_C1005_a = I_TWOBODYOVERLAP_K6xy_S_C1005_a+ABY*I_TWOBODYOVERLAP_I6x_S_C1005_a;
  Double I_TWOBODYOVERLAP_I5xy_Py_C1005_a = I_TWOBODYOVERLAP_K5x2y_S_C1005_a+ABY*I_TWOBODYOVERLAP_I5xy_S_C1005_a;
  Double I_TWOBODYOVERLAP_I5xz_Py_C1005_a = I_TWOBODYOVERLAP_K5xyz_S_C1005_a+ABY*I_TWOBODYOVERLAP_I5xz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4x2y_Py_C1005_a = I_TWOBODYOVERLAP_K4x3y_S_C1005_a+ABY*I_TWOBODYOVERLAP_I4x2y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4xyz_Py_C1005_a = I_TWOBODYOVERLAP_K4x2yz_S_C1005_a+ABY*I_TWOBODYOVERLAP_I4xyz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4x2z_Py_C1005_a = I_TWOBODYOVERLAP_K4xy2z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I4x2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3x3y_Py_C1005_a = I_TWOBODYOVERLAP_K3x4y_S_C1005_a+ABY*I_TWOBODYOVERLAP_I3x3y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Py_C1005_a = I_TWOBODYOVERLAP_K3x3yz_S_C1005_a+ABY*I_TWOBODYOVERLAP_I3x2yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Py_C1005_a = I_TWOBODYOVERLAP_K3x2y2z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I3xy2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3x3z_Py_C1005_a = I_TWOBODYOVERLAP_K3xy3z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I3x3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x4y_Py_C1005_a = I_TWOBODYOVERLAP_K2x5y_S_C1005_a+ABY*I_TWOBODYOVERLAP_I2x4y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Py_C1005_a = I_TWOBODYOVERLAP_K2x4yz_S_C1005_a+ABY*I_TWOBODYOVERLAP_I2x3yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Py_C1005_a = I_TWOBODYOVERLAP_K2x3y2z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Py_C1005_a = I_TWOBODYOVERLAP_K2x2y3z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I2xy3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x4z_Py_C1005_a = I_TWOBODYOVERLAP_K2xy4z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I2x4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix5y_Py_C1005_a = I_TWOBODYOVERLAP_Kx6y_S_C1005_a+ABY*I_TWOBODYOVERLAP_Ix5y_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Py_C1005_a = I_TWOBODYOVERLAP_Kx5yz_S_C1005_a+ABY*I_TWOBODYOVERLAP_Ix4yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Py_C1005_a = I_TWOBODYOVERLAP_Kx4y2z_S_C1005_a+ABY*I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Py_C1005_a = I_TWOBODYOVERLAP_Kx3y3z_S_C1005_a+ABY*I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Py_C1005_a = I_TWOBODYOVERLAP_Kx2y4z_S_C1005_a+ABY*I_TWOBODYOVERLAP_Ixy4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix5z_Py_C1005_a = I_TWOBODYOVERLAP_Kxy5z_S_C1005_a+ABY*I_TWOBODYOVERLAP_Ix5z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I6y_Py_C1005_a = I_TWOBODYOVERLAP_K7y_S_C1005_a+ABY*I_TWOBODYOVERLAP_I6y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I5yz_Py_C1005_a = I_TWOBODYOVERLAP_K6yz_S_C1005_a+ABY*I_TWOBODYOVERLAP_I5yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4y2z_Py_C1005_a = I_TWOBODYOVERLAP_K5y2z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I4y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3y3z_Py_C1005_a = I_TWOBODYOVERLAP_K4y3z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I3y3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2y4z_Py_C1005_a = I_TWOBODYOVERLAP_K3y4z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I2y4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Iy5z_Py_C1005_a = I_TWOBODYOVERLAP_K2y5z_S_C1005_a+ABY*I_TWOBODYOVERLAP_Iy5z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I6z_Py_C1005_a = I_TWOBODYOVERLAP_Ky6z_S_C1005_a+ABY*I_TWOBODYOVERLAP_I6z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I6x_Pz_C1005_a = I_TWOBODYOVERLAP_K6xz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I6x_S_C1005_a;
  Double I_TWOBODYOVERLAP_I5xy_Pz_C1005_a = I_TWOBODYOVERLAP_K5xyz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I5xy_S_C1005_a;
  Double I_TWOBODYOVERLAP_I5xz_Pz_C1005_a = I_TWOBODYOVERLAP_K5x2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I5xz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4x2y_Pz_C1005_a = I_TWOBODYOVERLAP_K4x2yz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I4x2y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4xyz_Pz_C1005_a = I_TWOBODYOVERLAP_K4xy2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I4xyz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4x2z_Pz_C1005_a = I_TWOBODYOVERLAP_K4x3z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I4x2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3x3y_Pz_C1005_a = I_TWOBODYOVERLAP_K3x3yz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I3x3y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Pz_C1005_a = I_TWOBODYOVERLAP_K3x2y2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I3x2yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Pz_C1005_a = I_TWOBODYOVERLAP_K3xy3z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I3xy2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3x3z_Pz_C1005_a = I_TWOBODYOVERLAP_K3x4z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I3x3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x4y_Pz_C1005_a = I_TWOBODYOVERLAP_K2x4yz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I2x4y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Pz_C1005_a = I_TWOBODYOVERLAP_K2x3y2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I2x3yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Pz_C1005_a = I_TWOBODYOVERLAP_K2x2y3z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I2x2y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Pz_C1005_a = I_TWOBODYOVERLAP_K2xy4z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I2xy3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2x4z_Pz_C1005_a = I_TWOBODYOVERLAP_K2x5z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I2x4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix5y_Pz_C1005_a = I_TWOBODYOVERLAP_Kx5yz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Ix5y_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Pz_C1005_a = I_TWOBODYOVERLAP_Kx4y2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Ix4yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Pz_C1005_a = I_TWOBODYOVERLAP_Kx3y3z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Ix3y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Pz_C1005_a = I_TWOBODYOVERLAP_Kx2y4z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Ix2y3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Pz_C1005_a = I_TWOBODYOVERLAP_Kxy5z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Ixy4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Ix5z_Pz_C1005_a = I_TWOBODYOVERLAP_Kx6z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Ix5z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I6y_Pz_C1005_a = I_TWOBODYOVERLAP_K6yz_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I6y_S_C1005_a;
  Double I_TWOBODYOVERLAP_I5yz_Pz_C1005_a = I_TWOBODYOVERLAP_K5y2z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I5yz_S_C1005_a;
  Double I_TWOBODYOVERLAP_I4y2z_Pz_C1005_a = I_TWOBODYOVERLAP_K4y3z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I4y2z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I3y3z_Pz_C1005_a = I_TWOBODYOVERLAP_K3y4z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I3y3z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I2y4z_Pz_C1005_a = I_TWOBODYOVERLAP_K2y5z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I2y4z_S_C1005_a;
  Double I_TWOBODYOVERLAP_Iy5z_Pz_C1005_a = I_TWOBODYOVERLAP_Ky6z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_Iy5z_S_C1005_a;
  Double I_TWOBODYOVERLAP_I6z_Pz_C1005_a = I_TWOBODYOVERLAP_K7z_S_C1005_a+ABZ*I_TWOBODYOVERLAP_I6z_S_C1005_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_C5
   ************************************************************/
  abcd[0] = 2.0E0*I_TWOBODYOVERLAP_I6x_S_C5_a-5*I_TWOBODYOVERLAP_G4x_S_C5;
  abcd[1] = 2.0E0*I_TWOBODYOVERLAP_I5xy_S_C5_a-4*I_TWOBODYOVERLAP_G3xy_S_C5;
  abcd[2] = 2.0E0*I_TWOBODYOVERLAP_I5xz_S_C5_a-4*I_TWOBODYOVERLAP_G3xz_S_C5;
  abcd[3] = 2.0E0*I_TWOBODYOVERLAP_I4x2y_S_C5_a-3*I_TWOBODYOVERLAP_G2x2y_S_C5;
  abcd[4] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_S_C5_a-3*I_TWOBODYOVERLAP_G2xyz_S_C5;
  abcd[5] = 2.0E0*I_TWOBODYOVERLAP_I4x2z_S_C5_a-3*I_TWOBODYOVERLAP_G2x2z_S_C5;
  abcd[6] = 2.0E0*I_TWOBODYOVERLAP_I3x3y_S_C5_a-2*I_TWOBODYOVERLAP_Gx3y_S_C5;
  abcd[7] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_S_C5_a-2*I_TWOBODYOVERLAP_Gx2yz_S_C5;
  abcd[8] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_S_C5_a-2*I_TWOBODYOVERLAP_Gxy2z_S_C5;
  abcd[9] = 2.0E0*I_TWOBODYOVERLAP_I3x3z_S_C5_a-2*I_TWOBODYOVERLAP_Gx3z_S_C5;
  abcd[10] = 2.0E0*I_TWOBODYOVERLAP_I2x4y_S_C5_a-1*I_TWOBODYOVERLAP_G4y_S_C5;
  abcd[11] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_S_C5_a-1*I_TWOBODYOVERLAP_G3yz_S_C5;
  abcd[12] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_S_C5_a-1*I_TWOBODYOVERLAP_G2y2z_S_C5;
  abcd[13] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_S_C5_a-1*I_TWOBODYOVERLAP_Gy3z_S_C5;
  abcd[14] = 2.0E0*I_TWOBODYOVERLAP_I2x4z_S_C5_a-1*I_TWOBODYOVERLAP_G4z_S_C5;
  abcd[15] = 2.0E0*I_TWOBODYOVERLAP_Ix5y_S_C5_a;
  abcd[16] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_S_C5_a;
  abcd[17] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_S_C5_a;
  abcd[18] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_S_C5_a;
  abcd[19] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_S_C5_a;
  abcd[20] = 2.0E0*I_TWOBODYOVERLAP_Ix5z_S_C5_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P_C1005
   ************************************************************/
  abcd[21] = 2.0E0*I_TWOBODYOVERLAP_I6x_Px_C1005_a-5*I_TWOBODYOVERLAP_G4x_Px_C1005;
  abcd[22] = 2.0E0*I_TWOBODYOVERLAP_I5xy_Px_C1005_a-4*I_TWOBODYOVERLAP_G3xy_Px_C1005;
  abcd[23] = 2.0E0*I_TWOBODYOVERLAP_I5xz_Px_C1005_a-4*I_TWOBODYOVERLAP_G3xz_Px_C1005;
  abcd[24] = 2.0E0*I_TWOBODYOVERLAP_I4x2y_Px_C1005_a-3*I_TWOBODYOVERLAP_G2x2y_Px_C1005;
  abcd[25] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_Px_C1005_a-3*I_TWOBODYOVERLAP_G2xyz_Px_C1005;
  abcd[26] = 2.0E0*I_TWOBODYOVERLAP_I4x2z_Px_C1005_a-3*I_TWOBODYOVERLAP_G2x2z_Px_C1005;
  abcd[27] = 2.0E0*I_TWOBODYOVERLAP_I3x3y_Px_C1005_a-2*I_TWOBODYOVERLAP_Gx3y_Px_C1005;
  abcd[28] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_Px_C1005_a-2*I_TWOBODYOVERLAP_Gx2yz_Px_C1005;
  abcd[29] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_Px_C1005_a-2*I_TWOBODYOVERLAP_Gxy2z_Px_C1005;
  abcd[30] = 2.0E0*I_TWOBODYOVERLAP_I3x3z_Px_C1005_a-2*I_TWOBODYOVERLAP_Gx3z_Px_C1005;
  abcd[31] = 2.0E0*I_TWOBODYOVERLAP_I2x4y_Px_C1005_a-1*I_TWOBODYOVERLAP_G4y_Px_C1005;
  abcd[32] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_Px_C1005_a-1*I_TWOBODYOVERLAP_G3yz_Px_C1005;
  abcd[33] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_Px_C1005_a-1*I_TWOBODYOVERLAP_G2y2z_Px_C1005;
  abcd[34] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_Px_C1005_a-1*I_TWOBODYOVERLAP_Gy3z_Px_C1005;
  abcd[35] = 2.0E0*I_TWOBODYOVERLAP_I2x4z_Px_C1005_a-1*I_TWOBODYOVERLAP_G4z_Px_C1005;
  abcd[36] = 2.0E0*I_TWOBODYOVERLAP_Ix5y_Px_C1005_a;
  abcd[37] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_Px_C1005_a;
  abcd[38] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_Px_C1005_a;
  abcd[39] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_Px_C1005_a;
  abcd[40] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_Px_C1005_a;
  abcd[41] = 2.0E0*I_TWOBODYOVERLAP_Ix5z_Px_C1005_a;
  abcd[42] = 2.0E0*I_TWOBODYOVERLAP_I6x_Py_C1005_a-5*I_TWOBODYOVERLAP_G4x_Py_C1005;
  abcd[43] = 2.0E0*I_TWOBODYOVERLAP_I5xy_Py_C1005_a-4*I_TWOBODYOVERLAP_G3xy_Py_C1005;
  abcd[44] = 2.0E0*I_TWOBODYOVERLAP_I5xz_Py_C1005_a-4*I_TWOBODYOVERLAP_G3xz_Py_C1005;
  abcd[45] = 2.0E0*I_TWOBODYOVERLAP_I4x2y_Py_C1005_a-3*I_TWOBODYOVERLAP_G2x2y_Py_C1005;
  abcd[46] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_Py_C1005_a-3*I_TWOBODYOVERLAP_G2xyz_Py_C1005;
  abcd[47] = 2.0E0*I_TWOBODYOVERLAP_I4x2z_Py_C1005_a-3*I_TWOBODYOVERLAP_G2x2z_Py_C1005;
  abcd[48] = 2.0E0*I_TWOBODYOVERLAP_I3x3y_Py_C1005_a-2*I_TWOBODYOVERLAP_Gx3y_Py_C1005;
  abcd[49] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_Py_C1005_a-2*I_TWOBODYOVERLAP_Gx2yz_Py_C1005;
  abcd[50] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_Py_C1005_a-2*I_TWOBODYOVERLAP_Gxy2z_Py_C1005;
  abcd[51] = 2.0E0*I_TWOBODYOVERLAP_I3x3z_Py_C1005_a-2*I_TWOBODYOVERLAP_Gx3z_Py_C1005;
  abcd[52] = 2.0E0*I_TWOBODYOVERLAP_I2x4y_Py_C1005_a-1*I_TWOBODYOVERLAP_G4y_Py_C1005;
  abcd[53] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_Py_C1005_a-1*I_TWOBODYOVERLAP_G3yz_Py_C1005;
  abcd[54] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_Py_C1005_a-1*I_TWOBODYOVERLAP_G2y2z_Py_C1005;
  abcd[55] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_Py_C1005_a-1*I_TWOBODYOVERLAP_Gy3z_Py_C1005;
  abcd[56] = 2.0E0*I_TWOBODYOVERLAP_I2x4z_Py_C1005_a-1*I_TWOBODYOVERLAP_G4z_Py_C1005;
  abcd[57] = 2.0E0*I_TWOBODYOVERLAP_Ix5y_Py_C1005_a;
  abcd[58] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_Py_C1005_a;
  abcd[59] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_Py_C1005_a;
  abcd[60] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_Py_C1005_a;
  abcd[61] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_Py_C1005_a;
  abcd[62] = 2.0E0*I_TWOBODYOVERLAP_Ix5z_Py_C1005_a;
  abcd[63] = 2.0E0*I_TWOBODYOVERLAP_I6x_Pz_C1005_a-5*I_TWOBODYOVERLAP_G4x_Pz_C1005;
  abcd[64] = 2.0E0*I_TWOBODYOVERLAP_I5xy_Pz_C1005_a-4*I_TWOBODYOVERLAP_G3xy_Pz_C1005;
  abcd[65] = 2.0E0*I_TWOBODYOVERLAP_I5xz_Pz_C1005_a-4*I_TWOBODYOVERLAP_G3xz_Pz_C1005;
  abcd[66] = 2.0E0*I_TWOBODYOVERLAP_I4x2y_Pz_C1005_a-3*I_TWOBODYOVERLAP_G2x2y_Pz_C1005;
  abcd[67] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_Pz_C1005_a-3*I_TWOBODYOVERLAP_G2xyz_Pz_C1005;
  abcd[68] = 2.0E0*I_TWOBODYOVERLAP_I4x2z_Pz_C1005_a-3*I_TWOBODYOVERLAP_G2x2z_Pz_C1005;
  abcd[69] = 2.0E0*I_TWOBODYOVERLAP_I3x3y_Pz_C1005_a-2*I_TWOBODYOVERLAP_Gx3y_Pz_C1005;
  abcd[70] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_Pz_C1005_a-2*I_TWOBODYOVERLAP_Gx2yz_Pz_C1005;
  abcd[71] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_Pz_C1005_a-2*I_TWOBODYOVERLAP_Gxy2z_Pz_C1005;
  abcd[72] = 2.0E0*I_TWOBODYOVERLAP_I3x3z_Pz_C1005_a-2*I_TWOBODYOVERLAP_Gx3z_Pz_C1005;
  abcd[73] = 2.0E0*I_TWOBODYOVERLAP_I2x4y_Pz_C1005_a-1*I_TWOBODYOVERLAP_G4y_Pz_C1005;
  abcd[74] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_Pz_C1005_a-1*I_TWOBODYOVERLAP_G3yz_Pz_C1005;
  abcd[75] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_Pz_C1005_a-1*I_TWOBODYOVERLAP_G2y2z_Pz_C1005;
  abcd[76] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_Pz_C1005_a-1*I_TWOBODYOVERLAP_Gy3z_Pz_C1005;
  abcd[77] = 2.0E0*I_TWOBODYOVERLAP_I2x4z_Pz_C1005_a-1*I_TWOBODYOVERLAP_G4z_Pz_C1005;
  abcd[78] = 2.0E0*I_TWOBODYOVERLAP_Ix5y_Pz_C1005_a;
  abcd[79] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_Pz_C1005_a;
  abcd[80] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_Pz_C1005_a;
  abcd[81] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_Pz_C1005_a;
  abcd[82] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_Pz_C1005_a;
  abcd[83] = 2.0E0*I_TWOBODYOVERLAP_Ix5z_Pz_C1005_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_C5
   ************************************************************/
  abcd[84] = 2.0E0*I_TWOBODYOVERLAP_I5xy_S_C5_a;
  abcd[85] = 2.0E0*I_TWOBODYOVERLAP_I4x2y_S_C5_a-1*I_TWOBODYOVERLAP_G4x_S_C5;
  abcd[86] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_S_C5_a;
  abcd[87] = 2.0E0*I_TWOBODYOVERLAP_I3x3y_S_C5_a-2*I_TWOBODYOVERLAP_G3xy_S_C5;
  abcd[88] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_S_C5_a-1*I_TWOBODYOVERLAP_G3xz_S_C5;
  abcd[89] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_S_C5_a;
  abcd[90] = 2.0E0*I_TWOBODYOVERLAP_I2x4y_S_C5_a-3*I_TWOBODYOVERLAP_G2x2y_S_C5;
  abcd[91] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_S_C5_a-2*I_TWOBODYOVERLAP_G2xyz_S_C5;
  abcd[92] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_S_C5_a-1*I_TWOBODYOVERLAP_G2x2z_S_C5;
  abcd[93] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_S_C5_a;
  abcd[94] = 2.0E0*I_TWOBODYOVERLAP_Ix5y_S_C5_a-4*I_TWOBODYOVERLAP_Gx3y_S_C5;
  abcd[95] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_S_C5_a-3*I_TWOBODYOVERLAP_Gx2yz_S_C5;
  abcd[96] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_S_C5_a-2*I_TWOBODYOVERLAP_Gxy2z_S_C5;
  abcd[97] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_S_C5_a-1*I_TWOBODYOVERLAP_Gx3z_S_C5;
  abcd[98] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_S_C5_a;
  abcd[99] = 2.0E0*I_TWOBODYOVERLAP_I6y_S_C5_a-5*I_TWOBODYOVERLAP_G4y_S_C5;
  abcd[100] = 2.0E0*I_TWOBODYOVERLAP_I5yz_S_C5_a-4*I_TWOBODYOVERLAP_G3yz_S_C5;
  abcd[101] = 2.0E0*I_TWOBODYOVERLAP_I4y2z_S_C5_a-3*I_TWOBODYOVERLAP_G2y2z_S_C5;
  abcd[102] = 2.0E0*I_TWOBODYOVERLAP_I3y3z_S_C5_a-2*I_TWOBODYOVERLAP_Gy3z_S_C5;
  abcd[103] = 2.0E0*I_TWOBODYOVERLAP_I2y4z_S_C5_a-1*I_TWOBODYOVERLAP_G4z_S_C5;
  abcd[104] = 2.0E0*I_TWOBODYOVERLAP_Iy5z_S_C5_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P_C1005
   ************************************************************/
  abcd[105] = 2.0E0*I_TWOBODYOVERLAP_I5xy_Px_C1005_a;
  abcd[106] = 2.0E0*I_TWOBODYOVERLAP_I4x2y_Px_C1005_a-1*I_TWOBODYOVERLAP_G4x_Px_C1005;
  abcd[107] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_Px_C1005_a;
  abcd[108] = 2.0E0*I_TWOBODYOVERLAP_I3x3y_Px_C1005_a-2*I_TWOBODYOVERLAP_G3xy_Px_C1005;
  abcd[109] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_Px_C1005_a-1*I_TWOBODYOVERLAP_G3xz_Px_C1005;
  abcd[110] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_Px_C1005_a;
  abcd[111] = 2.0E0*I_TWOBODYOVERLAP_I2x4y_Px_C1005_a-3*I_TWOBODYOVERLAP_G2x2y_Px_C1005;
  abcd[112] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_Px_C1005_a-2*I_TWOBODYOVERLAP_G2xyz_Px_C1005;
  abcd[113] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_Px_C1005_a-1*I_TWOBODYOVERLAP_G2x2z_Px_C1005;
  abcd[114] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_Px_C1005_a;
  abcd[115] = 2.0E0*I_TWOBODYOVERLAP_Ix5y_Px_C1005_a-4*I_TWOBODYOVERLAP_Gx3y_Px_C1005;
  abcd[116] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_Px_C1005_a-3*I_TWOBODYOVERLAP_Gx2yz_Px_C1005;
  abcd[117] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_Px_C1005_a-2*I_TWOBODYOVERLAP_Gxy2z_Px_C1005;
  abcd[118] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_Px_C1005_a-1*I_TWOBODYOVERLAP_Gx3z_Px_C1005;
  abcd[119] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_Px_C1005_a;
  abcd[120] = 2.0E0*I_TWOBODYOVERLAP_I6y_Px_C1005_a-5*I_TWOBODYOVERLAP_G4y_Px_C1005;
  abcd[121] = 2.0E0*I_TWOBODYOVERLAP_I5yz_Px_C1005_a-4*I_TWOBODYOVERLAP_G3yz_Px_C1005;
  abcd[122] = 2.0E0*I_TWOBODYOVERLAP_I4y2z_Px_C1005_a-3*I_TWOBODYOVERLAP_G2y2z_Px_C1005;
  abcd[123] = 2.0E0*I_TWOBODYOVERLAP_I3y3z_Px_C1005_a-2*I_TWOBODYOVERLAP_Gy3z_Px_C1005;
  abcd[124] = 2.0E0*I_TWOBODYOVERLAP_I2y4z_Px_C1005_a-1*I_TWOBODYOVERLAP_G4z_Px_C1005;
  abcd[125] = 2.0E0*I_TWOBODYOVERLAP_Iy5z_Px_C1005_a;
  abcd[126] = 2.0E0*I_TWOBODYOVERLAP_I5xy_Py_C1005_a;
  abcd[127] = 2.0E0*I_TWOBODYOVERLAP_I4x2y_Py_C1005_a-1*I_TWOBODYOVERLAP_G4x_Py_C1005;
  abcd[128] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_Py_C1005_a;
  abcd[129] = 2.0E0*I_TWOBODYOVERLAP_I3x3y_Py_C1005_a-2*I_TWOBODYOVERLAP_G3xy_Py_C1005;
  abcd[130] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_Py_C1005_a-1*I_TWOBODYOVERLAP_G3xz_Py_C1005;
  abcd[131] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_Py_C1005_a;
  abcd[132] = 2.0E0*I_TWOBODYOVERLAP_I2x4y_Py_C1005_a-3*I_TWOBODYOVERLAP_G2x2y_Py_C1005;
  abcd[133] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_Py_C1005_a-2*I_TWOBODYOVERLAP_G2xyz_Py_C1005;
  abcd[134] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_Py_C1005_a-1*I_TWOBODYOVERLAP_G2x2z_Py_C1005;
  abcd[135] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_Py_C1005_a;
  abcd[136] = 2.0E0*I_TWOBODYOVERLAP_Ix5y_Py_C1005_a-4*I_TWOBODYOVERLAP_Gx3y_Py_C1005;
  abcd[137] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_Py_C1005_a-3*I_TWOBODYOVERLAP_Gx2yz_Py_C1005;
  abcd[138] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_Py_C1005_a-2*I_TWOBODYOVERLAP_Gxy2z_Py_C1005;
  abcd[139] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_Py_C1005_a-1*I_TWOBODYOVERLAP_Gx3z_Py_C1005;
  abcd[140] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_Py_C1005_a;
  abcd[141] = 2.0E0*I_TWOBODYOVERLAP_I6y_Py_C1005_a-5*I_TWOBODYOVERLAP_G4y_Py_C1005;
  abcd[142] = 2.0E0*I_TWOBODYOVERLAP_I5yz_Py_C1005_a-4*I_TWOBODYOVERLAP_G3yz_Py_C1005;
  abcd[143] = 2.0E0*I_TWOBODYOVERLAP_I4y2z_Py_C1005_a-3*I_TWOBODYOVERLAP_G2y2z_Py_C1005;
  abcd[144] = 2.0E0*I_TWOBODYOVERLAP_I3y3z_Py_C1005_a-2*I_TWOBODYOVERLAP_Gy3z_Py_C1005;
  abcd[145] = 2.0E0*I_TWOBODYOVERLAP_I2y4z_Py_C1005_a-1*I_TWOBODYOVERLAP_G4z_Py_C1005;
  abcd[146] = 2.0E0*I_TWOBODYOVERLAP_Iy5z_Py_C1005_a;
  abcd[147] = 2.0E0*I_TWOBODYOVERLAP_I5xy_Pz_C1005_a;
  abcd[148] = 2.0E0*I_TWOBODYOVERLAP_I4x2y_Pz_C1005_a-1*I_TWOBODYOVERLAP_G4x_Pz_C1005;
  abcd[149] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_Pz_C1005_a;
  abcd[150] = 2.0E0*I_TWOBODYOVERLAP_I3x3y_Pz_C1005_a-2*I_TWOBODYOVERLAP_G3xy_Pz_C1005;
  abcd[151] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_Pz_C1005_a-1*I_TWOBODYOVERLAP_G3xz_Pz_C1005;
  abcd[152] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_Pz_C1005_a;
  abcd[153] = 2.0E0*I_TWOBODYOVERLAP_I2x4y_Pz_C1005_a-3*I_TWOBODYOVERLAP_G2x2y_Pz_C1005;
  abcd[154] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_Pz_C1005_a-2*I_TWOBODYOVERLAP_G2xyz_Pz_C1005;
  abcd[155] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_Pz_C1005_a-1*I_TWOBODYOVERLAP_G2x2z_Pz_C1005;
  abcd[156] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_Pz_C1005_a;
  abcd[157] = 2.0E0*I_TWOBODYOVERLAP_Ix5y_Pz_C1005_a-4*I_TWOBODYOVERLAP_Gx3y_Pz_C1005;
  abcd[158] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_Pz_C1005_a-3*I_TWOBODYOVERLAP_Gx2yz_Pz_C1005;
  abcd[159] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_Pz_C1005_a-2*I_TWOBODYOVERLAP_Gxy2z_Pz_C1005;
  abcd[160] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_Pz_C1005_a-1*I_TWOBODYOVERLAP_Gx3z_Pz_C1005;
  abcd[161] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_Pz_C1005_a;
  abcd[162] = 2.0E0*I_TWOBODYOVERLAP_I6y_Pz_C1005_a-5*I_TWOBODYOVERLAP_G4y_Pz_C1005;
  abcd[163] = 2.0E0*I_TWOBODYOVERLAP_I5yz_Pz_C1005_a-4*I_TWOBODYOVERLAP_G3yz_Pz_C1005;
  abcd[164] = 2.0E0*I_TWOBODYOVERLAP_I4y2z_Pz_C1005_a-3*I_TWOBODYOVERLAP_G2y2z_Pz_C1005;
  abcd[165] = 2.0E0*I_TWOBODYOVERLAP_I3y3z_Pz_C1005_a-2*I_TWOBODYOVERLAP_Gy3z_Pz_C1005;
  abcd[166] = 2.0E0*I_TWOBODYOVERLAP_I2y4z_Pz_C1005_a-1*I_TWOBODYOVERLAP_G4z_Pz_C1005;
  abcd[167] = 2.0E0*I_TWOBODYOVERLAP_Iy5z_Pz_C1005_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_C5_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_C5
   ************************************************************/
  abcd[168] = 2.0E0*I_TWOBODYOVERLAP_I5xz_S_C5_a;
  abcd[169] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_S_C5_a;
  abcd[170] = 2.0E0*I_TWOBODYOVERLAP_I4x2z_S_C5_a-1*I_TWOBODYOVERLAP_G4x_S_C5;
  abcd[171] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_S_C5_a;
  abcd[172] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_S_C5_a-1*I_TWOBODYOVERLAP_G3xy_S_C5;
  abcd[173] = 2.0E0*I_TWOBODYOVERLAP_I3x3z_S_C5_a-2*I_TWOBODYOVERLAP_G3xz_S_C5;
  abcd[174] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_S_C5_a;
  abcd[175] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_S_C5_a-1*I_TWOBODYOVERLAP_G2x2y_S_C5;
  abcd[176] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_S_C5_a-2*I_TWOBODYOVERLAP_G2xyz_S_C5;
  abcd[177] = 2.0E0*I_TWOBODYOVERLAP_I2x4z_S_C5_a-3*I_TWOBODYOVERLAP_G2x2z_S_C5;
  abcd[178] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_S_C5_a;
  abcd[179] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_S_C5_a-1*I_TWOBODYOVERLAP_Gx3y_S_C5;
  abcd[180] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_S_C5_a-2*I_TWOBODYOVERLAP_Gx2yz_S_C5;
  abcd[181] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_S_C5_a-3*I_TWOBODYOVERLAP_Gxy2z_S_C5;
  abcd[182] = 2.0E0*I_TWOBODYOVERLAP_Ix5z_S_C5_a-4*I_TWOBODYOVERLAP_Gx3z_S_C5;
  abcd[183] = 2.0E0*I_TWOBODYOVERLAP_I5yz_S_C5_a;
  abcd[184] = 2.0E0*I_TWOBODYOVERLAP_I4y2z_S_C5_a-1*I_TWOBODYOVERLAP_G4y_S_C5;
  abcd[185] = 2.0E0*I_TWOBODYOVERLAP_I3y3z_S_C5_a-2*I_TWOBODYOVERLAP_G3yz_S_C5;
  abcd[186] = 2.0E0*I_TWOBODYOVERLAP_I2y4z_S_C5_a-3*I_TWOBODYOVERLAP_G2y2z_S_C5;
  abcd[187] = 2.0E0*I_TWOBODYOVERLAP_Iy5z_S_C5_a-4*I_TWOBODYOVERLAP_Gy3z_S_C5;
  abcd[188] = 2.0E0*I_TWOBODYOVERLAP_I6z_S_C5_a-5*I_TWOBODYOVERLAP_G4z_S_C5;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P_C1005_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P_C1005
   ************************************************************/
  abcd[189] = 2.0E0*I_TWOBODYOVERLAP_I5xz_Px_C1005_a;
  abcd[190] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_Px_C1005_a;
  abcd[191] = 2.0E0*I_TWOBODYOVERLAP_I4x2z_Px_C1005_a-1*I_TWOBODYOVERLAP_G4x_Px_C1005;
  abcd[192] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_Px_C1005_a;
  abcd[193] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_Px_C1005_a-1*I_TWOBODYOVERLAP_G3xy_Px_C1005;
  abcd[194] = 2.0E0*I_TWOBODYOVERLAP_I3x3z_Px_C1005_a-2*I_TWOBODYOVERLAP_G3xz_Px_C1005;
  abcd[195] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_Px_C1005_a;
  abcd[196] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_Px_C1005_a-1*I_TWOBODYOVERLAP_G2x2y_Px_C1005;
  abcd[197] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_Px_C1005_a-2*I_TWOBODYOVERLAP_G2xyz_Px_C1005;
  abcd[198] = 2.0E0*I_TWOBODYOVERLAP_I2x4z_Px_C1005_a-3*I_TWOBODYOVERLAP_G2x2z_Px_C1005;
  abcd[199] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_Px_C1005_a;
  abcd[200] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_Px_C1005_a-1*I_TWOBODYOVERLAP_Gx3y_Px_C1005;
  abcd[201] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_Px_C1005_a-2*I_TWOBODYOVERLAP_Gx2yz_Px_C1005;
  abcd[202] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_Px_C1005_a-3*I_TWOBODYOVERLAP_Gxy2z_Px_C1005;
  abcd[203] = 2.0E0*I_TWOBODYOVERLAP_Ix5z_Px_C1005_a-4*I_TWOBODYOVERLAP_Gx3z_Px_C1005;
  abcd[204] = 2.0E0*I_TWOBODYOVERLAP_I5yz_Px_C1005_a;
  abcd[205] = 2.0E0*I_TWOBODYOVERLAP_I4y2z_Px_C1005_a-1*I_TWOBODYOVERLAP_G4y_Px_C1005;
  abcd[206] = 2.0E0*I_TWOBODYOVERLAP_I3y3z_Px_C1005_a-2*I_TWOBODYOVERLAP_G3yz_Px_C1005;
  abcd[207] = 2.0E0*I_TWOBODYOVERLAP_I2y4z_Px_C1005_a-3*I_TWOBODYOVERLAP_G2y2z_Px_C1005;
  abcd[208] = 2.0E0*I_TWOBODYOVERLAP_Iy5z_Px_C1005_a-4*I_TWOBODYOVERLAP_Gy3z_Px_C1005;
  abcd[209] = 2.0E0*I_TWOBODYOVERLAP_I6z_Px_C1005_a-5*I_TWOBODYOVERLAP_G4z_Px_C1005;
  abcd[210] = 2.0E0*I_TWOBODYOVERLAP_I5xz_Py_C1005_a;
  abcd[211] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_Py_C1005_a;
  abcd[212] = 2.0E0*I_TWOBODYOVERLAP_I4x2z_Py_C1005_a-1*I_TWOBODYOVERLAP_G4x_Py_C1005;
  abcd[213] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_Py_C1005_a;
  abcd[214] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_Py_C1005_a-1*I_TWOBODYOVERLAP_G3xy_Py_C1005;
  abcd[215] = 2.0E0*I_TWOBODYOVERLAP_I3x3z_Py_C1005_a-2*I_TWOBODYOVERLAP_G3xz_Py_C1005;
  abcd[216] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_Py_C1005_a;
  abcd[217] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_Py_C1005_a-1*I_TWOBODYOVERLAP_G2x2y_Py_C1005;
  abcd[218] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_Py_C1005_a-2*I_TWOBODYOVERLAP_G2xyz_Py_C1005;
  abcd[219] = 2.0E0*I_TWOBODYOVERLAP_I2x4z_Py_C1005_a-3*I_TWOBODYOVERLAP_G2x2z_Py_C1005;
  abcd[220] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_Py_C1005_a;
  abcd[221] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_Py_C1005_a-1*I_TWOBODYOVERLAP_Gx3y_Py_C1005;
  abcd[222] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_Py_C1005_a-2*I_TWOBODYOVERLAP_Gx2yz_Py_C1005;
  abcd[223] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_Py_C1005_a-3*I_TWOBODYOVERLAP_Gxy2z_Py_C1005;
  abcd[224] = 2.0E0*I_TWOBODYOVERLAP_Ix5z_Py_C1005_a-4*I_TWOBODYOVERLAP_Gx3z_Py_C1005;
  abcd[225] = 2.0E0*I_TWOBODYOVERLAP_I5yz_Py_C1005_a;
  abcd[226] = 2.0E0*I_TWOBODYOVERLAP_I4y2z_Py_C1005_a-1*I_TWOBODYOVERLAP_G4y_Py_C1005;
  abcd[227] = 2.0E0*I_TWOBODYOVERLAP_I3y3z_Py_C1005_a-2*I_TWOBODYOVERLAP_G3yz_Py_C1005;
  abcd[228] = 2.0E0*I_TWOBODYOVERLAP_I2y4z_Py_C1005_a-3*I_TWOBODYOVERLAP_G2y2z_Py_C1005;
  abcd[229] = 2.0E0*I_TWOBODYOVERLAP_Iy5z_Py_C1005_a-4*I_TWOBODYOVERLAP_Gy3z_Py_C1005;
  abcd[230] = 2.0E0*I_TWOBODYOVERLAP_I6z_Py_C1005_a-5*I_TWOBODYOVERLAP_G4z_Py_C1005;
  abcd[231] = 2.0E0*I_TWOBODYOVERLAP_I5xz_Pz_C1005_a;
  abcd[232] = 2.0E0*I_TWOBODYOVERLAP_I4xyz_Pz_C1005_a;
  abcd[233] = 2.0E0*I_TWOBODYOVERLAP_I4x2z_Pz_C1005_a-1*I_TWOBODYOVERLAP_G4x_Pz_C1005;
  abcd[234] = 2.0E0*I_TWOBODYOVERLAP_I3x2yz_Pz_C1005_a;
  abcd[235] = 2.0E0*I_TWOBODYOVERLAP_I3xy2z_Pz_C1005_a-1*I_TWOBODYOVERLAP_G3xy_Pz_C1005;
  abcd[236] = 2.0E0*I_TWOBODYOVERLAP_I3x3z_Pz_C1005_a-2*I_TWOBODYOVERLAP_G3xz_Pz_C1005;
  abcd[237] = 2.0E0*I_TWOBODYOVERLAP_I2x3yz_Pz_C1005_a;
  abcd[238] = 2.0E0*I_TWOBODYOVERLAP_I2x2y2z_Pz_C1005_a-1*I_TWOBODYOVERLAP_G2x2y_Pz_C1005;
  abcd[239] = 2.0E0*I_TWOBODYOVERLAP_I2xy3z_Pz_C1005_a-2*I_TWOBODYOVERLAP_G2xyz_Pz_C1005;
  abcd[240] = 2.0E0*I_TWOBODYOVERLAP_I2x4z_Pz_C1005_a-3*I_TWOBODYOVERLAP_G2x2z_Pz_C1005;
  abcd[241] = 2.0E0*I_TWOBODYOVERLAP_Ix4yz_Pz_C1005_a;
  abcd[242] = 2.0E0*I_TWOBODYOVERLAP_Ix3y2z_Pz_C1005_a-1*I_TWOBODYOVERLAP_Gx3y_Pz_C1005;
  abcd[243] = 2.0E0*I_TWOBODYOVERLAP_Ix2y3z_Pz_C1005_a-2*I_TWOBODYOVERLAP_Gx2yz_Pz_C1005;
  abcd[244] = 2.0E0*I_TWOBODYOVERLAP_Ixy4z_Pz_C1005_a-3*I_TWOBODYOVERLAP_Gxy2z_Pz_C1005;
  abcd[245] = 2.0E0*I_TWOBODYOVERLAP_Ix5z_Pz_C1005_a-4*I_TWOBODYOVERLAP_Gx3z_Pz_C1005;
  abcd[246] = 2.0E0*I_TWOBODYOVERLAP_I5yz_Pz_C1005_a;
  abcd[247] = 2.0E0*I_TWOBODYOVERLAP_I4y2z_Pz_C1005_a-1*I_TWOBODYOVERLAP_G4y_Pz_C1005;
  abcd[248] = 2.0E0*I_TWOBODYOVERLAP_I3y3z_Pz_C1005_a-2*I_TWOBODYOVERLAP_G3yz_Pz_C1005;
  abcd[249] = 2.0E0*I_TWOBODYOVERLAP_I2y4z_Pz_C1005_a-3*I_TWOBODYOVERLAP_G2y2z_Pz_C1005;
  abcd[250] = 2.0E0*I_TWOBODYOVERLAP_Iy5z_Pz_C1005_a-4*I_TWOBODYOVERLAP_Gy3z_Pz_C1005;
  abcd[251] = 2.0E0*I_TWOBODYOVERLAP_I6z_Pz_C1005_a-5*I_TWOBODYOVERLAP_G4z_Pz_C1005;
}
