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
// BRA1 as redundant position, total RHS integrals evaluated as: 1751
// BRA2 as redundant position, total RHS integrals evaluated as: 1203
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

void hgp_os_twobodyoverlap_f_sp_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_H5x_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S_C3_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_F3x_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xy_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xz_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2y_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fxyz_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2z_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3y_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2yz_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fy2z_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3z_S_C3_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Px_S_C3 = 0.0E0;
  Double I_TWOBODYOVERLAP_Py_S_C3 = 0.0E0;
  Double I_TWOBODYOVERLAP_Pz_S_C3 = 0.0E0;
  Double I_TWOBODYOVERLAP_I6x_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H5x_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S_C1003_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G4x_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3x_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xy_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xz_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2y_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fxyz_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2z_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3y_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2yz_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fy2z_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3z_S_C1003_a = 0.0E0;
  Double I_TWOBODYOVERLAP_D2x_S_C1003 = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxy_S_C1003 = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxz_S_C1003 = 0.0E0;
  Double I_TWOBODYOVERLAP_D2y_S_C1003 = 0.0E0;
  Double I_TWOBODYOVERLAP_Dyz_S_C1003 = 0.0E0;
  Double I_TWOBODYOVERLAP_D2z_S_C1003 = 0.0E0;
  Double I_TWOBODYOVERLAP_Px_S_C1003 = 0.0E0;
  Double I_TWOBODYOVERLAP_Py_S_C1003 = 0.0E0;
  Double I_TWOBODYOVERLAP_Pz_S_C1003 = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C3_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs = ic2*alpha*alpha;
    I_TWOBODYOVERLAP_H5x_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S_C3_aa += SQ_TWOBODYOVERLAP_H_S_C3_aa_coefs*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_F_S_C3_a_coefs = ic2*alpha;
    I_TWOBODYOVERLAP_F3x_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_F3x_S_vrr;
    I_TWOBODYOVERLAP_F2xy_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_F2xy_S_vrr;
    I_TWOBODYOVERLAP_F2xz_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_F2xz_S_vrr;
    I_TWOBODYOVERLAP_Fx2y_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    I_TWOBODYOVERLAP_Fxyz_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    I_TWOBODYOVERLAP_Fx2z_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    I_TWOBODYOVERLAP_F3y_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_F3y_S_vrr;
    I_TWOBODYOVERLAP_F2yz_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_F2yz_S_vrr;
    I_TWOBODYOVERLAP_Fy2z_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    I_TWOBODYOVERLAP_F3z_S_C3_a += SQ_TWOBODYOVERLAP_F_S_C3_a_coefs*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S_C3
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_P_S_C3_coefs = ic2;
    I_TWOBODYOVERLAP_Px_S_C3 += SQ_TWOBODYOVERLAP_P_S_C3_coefs*I_TWOBODYOVERLAP_Px_S_vrr;
    I_TWOBODYOVERLAP_Py_S_C3 += SQ_TWOBODYOVERLAP_P_S_C3_coefs*I_TWOBODYOVERLAP_Py_S_vrr;
    I_TWOBODYOVERLAP_Pz_S_C3 += SQ_TWOBODYOVERLAP_P_S_C3_coefs*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S_C1003_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs = ic2_1*alpha*alpha;
    I_TWOBODYOVERLAP_I6x_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S_C1003_aa += SQ_TWOBODYOVERLAP_I_S_C1003_aa_coefs*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C1003_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs = ic2_1*alpha*alpha;
    I_TWOBODYOVERLAP_H5x_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S_C1003_aa += SQ_TWOBODYOVERLAP_H_S_C1003_aa_coefs*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S_C1003_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs = ic2_1*alpha;
    I_TWOBODYOVERLAP_G4x_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S_C1003_a += SQ_TWOBODYOVERLAP_G_S_C1003_a_coefs*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C1003_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs = ic2_1*alpha;
    I_TWOBODYOVERLAP_F3x_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_F3x_S_vrr;
    I_TWOBODYOVERLAP_F2xy_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_F2xy_S_vrr;
    I_TWOBODYOVERLAP_F2xz_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_F2xz_S_vrr;
    I_TWOBODYOVERLAP_Fx2y_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    I_TWOBODYOVERLAP_Fxyz_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    I_TWOBODYOVERLAP_Fx2z_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    I_TWOBODYOVERLAP_F3y_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_F3y_S_vrr;
    I_TWOBODYOVERLAP_F2yz_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_F2yz_S_vrr;
    I_TWOBODYOVERLAP_Fy2z_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    I_TWOBODYOVERLAP_F3z_S_C1003_a += SQ_TWOBODYOVERLAP_F_S_C1003_a_coefs*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S_C1003
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_D_S_C1003_coefs = ic2_1;
    I_TWOBODYOVERLAP_D2x_S_C1003 += SQ_TWOBODYOVERLAP_D_S_C1003_coefs*I_TWOBODYOVERLAP_D2x_S_vrr;
    I_TWOBODYOVERLAP_Dxy_S_C1003 += SQ_TWOBODYOVERLAP_D_S_C1003_coefs*I_TWOBODYOVERLAP_Dxy_S_vrr;
    I_TWOBODYOVERLAP_Dxz_S_C1003 += SQ_TWOBODYOVERLAP_D_S_C1003_coefs*I_TWOBODYOVERLAP_Dxz_S_vrr;
    I_TWOBODYOVERLAP_D2y_S_C1003 += SQ_TWOBODYOVERLAP_D_S_C1003_coefs*I_TWOBODYOVERLAP_D2y_S_vrr;
    I_TWOBODYOVERLAP_Dyz_S_C1003 += SQ_TWOBODYOVERLAP_D_S_C1003_coefs*I_TWOBODYOVERLAP_Dyz_S_vrr;
    I_TWOBODYOVERLAP_D2z_S_C1003 += SQ_TWOBODYOVERLAP_D_S_C1003_coefs*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S_C1003
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_P_S_C1003_coefs = ic2_1;
    I_TWOBODYOVERLAP_Px_S_C1003 += SQ_TWOBODYOVERLAP_P_S_C1003_coefs*I_TWOBODYOVERLAP_Px_S_vrr;
    I_TWOBODYOVERLAP_Py_S_C1003 += SQ_TWOBODYOVERLAP_P_S_C1003_coefs*I_TWOBODYOVERLAP_Py_S_vrr;
    I_TWOBODYOVERLAP_Pz_S_C1003 += SQ_TWOBODYOVERLAP_P_S_C1003_coefs*I_TWOBODYOVERLAP_Pz_S_vrr;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1003
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S_C1003
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C1003
   ************************************************************/
  Double I_TWOBODYOVERLAP_Px_Px_C1003 = I_TWOBODYOVERLAP_D2x_S_C1003+ABX*I_TWOBODYOVERLAP_Px_S_C1003;
  Double I_TWOBODYOVERLAP_Py_Px_C1003 = I_TWOBODYOVERLAP_Dxy_S_C1003+ABX*I_TWOBODYOVERLAP_Py_S_C1003;
  Double I_TWOBODYOVERLAP_Pz_Px_C1003 = I_TWOBODYOVERLAP_Dxz_S_C1003+ABX*I_TWOBODYOVERLAP_Pz_S_C1003;
  Double I_TWOBODYOVERLAP_Px_Py_C1003 = I_TWOBODYOVERLAP_Dxy_S_C1003+ABY*I_TWOBODYOVERLAP_Px_S_C1003;
  Double I_TWOBODYOVERLAP_Py_Py_C1003 = I_TWOBODYOVERLAP_D2y_S_C1003+ABY*I_TWOBODYOVERLAP_Py_S_C1003;
  Double I_TWOBODYOVERLAP_Pz_Py_C1003 = I_TWOBODYOVERLAP_Dyz_S_C1003+ABY*I_TWOBODYOVERLAP_Pz_S_C1003;
  Double I_TWOBODYOVERLAP_Px_Pz_C1003 = I_TWOBODYOVERLAP_Dxz_S_C1003+ABZ*I_TWOBODYOVERLAP_Px_S_C1003;
  Double I_TWOBODYOVERLAP_Py_Pz_C1003 = I_TWOBODYOVERLAP_Dyz_S_C1003+ABZ*I_TWOBODYOVERLAP_Py_S_C1003;
  Double I_TWOBODYOVERLAP_Pz_Pz_C1003 = I_TWOBODYOVERLAP_D2z_S_C1003+ABZ*I_TWOBODYOVERLAP_Pz_S_C1003;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C1003_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_F3x_Px_C1003_a = I_TWOBODYOVERLAP_G4x_S_C1003_a+ABX*I_TWOBODYOVERLAP_F3x_S_C1003_a;
  Double I_TWOBODYOVERLAP_F2xy_Px_C1003_a = I_TWOBODYOVERLAP_G3xy_S_C1003_a+ABX*I_TWOBODYOVERLAP_F2xy_S_C1003_a;
  Double I_TWOBODYOVERLAP_F2xz_Px_C1003_a = I_TWOBODYOVERLAP_G3xz_S_C1003_a+ABX*I_TWOBODYOVERLAP_F2xz_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fx2y_Px_C1003_a = I_TWOBODYOVERLAP_G2x2y_S_C1003_a+ABX*I_TWOBODYOVERLAP_Fx2y_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fxyz_Px_C1003_a = I_TWOBODYOVERLAP_G2xyz_S_C1003_a+ABX*I_TWOBODYOVERLAP_Fxyz_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fx2z_Px_C1003_a = I_TWOBODYOVERLAP_G2x2z_S_C1003_a+ABX*I_TWOBODYOVERLAP_Fx2z_S_C1003_a;
  Double I_TWOBODYOVERLAP_F3y_Px_C1003_a = I_TWOBODYOVERLAP_Gx3y_S_C1003_a+ABX*I_TWOBODYOVERLAP_F3y_S_C1003_a;
  Double I_TWOBODYOVERLAP_F2yz_Px_C1003_a = I_TWOBODYOVERLAP_Gx2yz_S_C1003_a+ABX*I_TWOBODYOVERLAP_F2yz_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fy2z_Px_C1003_a = I_TWOBODYOVERLAP_Gxy2z_S_C1003_a+ABX*I_TWOBODYOVERLAP_Fy2z_S_C1003_a;
  Double I_TWOBODYOVERLAP_F3z_Px_C1003_a = I_TWOBODYOVERLAP_Gx3z_S_C1003_a+ABX*I_TWOBODYOVERLAP_F3z_S_C1003_a;
  Double I_TWOBODYOVERLAP_F3x_Py_C1003_a = I_TWOBODYOVERLAP_G3xy_S_C1003_a+ABY*I_TWOBODYOVERLAP_F3x_S_C1003_a;
  Double I_TWOBODYOVERLAP_F2xy_Py_C1003_a = I_TWOBODYOVERLAP_G2x2y_S_C1003_a+ABY*I_TWOBODYOVERLAP_F2xy_S_C1003_a;
  Double I_TWOBODYOVERLAP_F2xz_Py_C1003_a = I_TWOBODYOVERLAP_G2xyz_S_C1003_a+ABY*I_TWOBODYOVERLAP_F2xz_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fx2y_Py_C1003_a = I_TWOBODYOVERLAP_Gx3y_S_C1003_a+ABY*I_TWOBODYOVERLAP_Fx2y_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fxyz_Py_C1003_a = I_TWOBODYOVERLAP_Gx2yz_S_C1003_a+ABY*I_TWOBODYOVERLAP_Fxyz_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fx2z_Py_C1003_a = I_TWOBODYOVERLAP_Gxy2z_S_C1003_a+ABY*I_TWOBODYOVERLAP_Fx2z_S_C1003_a;
  Double I_TWOBODYOVERLAP_F3y_Py_C1003_a = I_TWOBODYOVERLAP_G4y_S_C1003_a+ABY*I_TWOBODYOVERLAP_F3y_S_C1003_a;
  Double I_TWOBODYOVERLAP_F2yz_Py_C1003_a = I_TWOBODYOVERLAP_G3yz_S_C1003_a+ABY*I_TWOBODYOVERLAP_F2yz_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fy2z_Py_C1003_a = I_TWOBODYOVERLAP_G2y2z_S_C1003_a+ABY*I_TWOBODYOVERLAP_Fy2z_S_C1003_a;
  Double I_TWOBODYOVERLAP_F3z_Py_C1003_a = I_TWOBODYOVERLAP_Gy3z_S_C1003_a+ABY*I_TWOBODYOVERLAP_F3z_S_C1003_a;
  Double I_TWOBODYOVERLAP_F3x_Pz_C1003_a = I_TWOBODYOVERLAP_G3xz_S_C1003_a+ABZ*I_TWOBODYOVERLAP_F3x_S_C1003_a;
  Double I_TWOBODYOVERLAP_F2xy_Pz_C1003_a = I_TWOBODYOVERLAP_G2xyz_S_C1003_a+ABZ*I_TWOBODYOVERLAP_F2xy_S_C1003_a;
  Double I_TWOBODYOVERLAP_F2xz_Pz_C1003_a = I_TWOBODYOVERLAP_G2x2z_S_C1003_a+ABZ*I_TWOBODYOVERLAP_F2xz_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a = I_TWOBODYOVERLAP_Gx2yz_S_C1003_a+ABZ*I_TWOBODYOVERLAP_Fx2y_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a = I_TWOBODYOVERLAP_Gxy2z_S_C1003_a+ABZ*I_TWOBODYOVERLAP_Fxyz_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a = I_TWOBODYOVERLAP_Gx3z_S_C1003_a+ABZ*I_TWOBODYOVERLAP_Fx2z_S_C1003_a;
  Double I_TWOBODYOVERLAP_F3y_Pz_C1003_a = I_TWOBODYOVERLAP_G3yz_S_C1003_a+ABZ*I_TWOBODYOVERLAP_F3y_S_C1003_a;
  Double I_TWOBODYOVERLAP_F2yz_Pz_C1003_a = I_TWOBODYOVERLAP_G2y2z_S_C1003_a+ABZ*I_TWOBODYOVERLAP_F2yz_S_C1003_a;
  Double I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a = I_TWOBODYOVERLAP_Gy3z_S_C1003_a+ABZ*I_TWOBODYOVERLAP_Fy2z_S_C1003_a;
  Double I_TWOBODYOVERLAP_F3z_Pz_C1003_a = I_TWOBODYOVERLAP_G4z_S_C1003_a+ABZ*I_TWOBODYOVERLAP_F3z_S_C1003_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1003_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_C1003_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C1003_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_Px_C1003_aa = I_TWOBODYOVERLAP_I6x_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H5x_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H4xy_Px_C1003_aa = I_TWOBODYOVERLAP_I5xy_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H4xy_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H4xz_Px_C1003_aa = I_TWOBODYOVERLAP_I5xz_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H4xz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3x2y_Px_C1003_aa = I_TWOBODYOVERLAP_I4x2y_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H3x2y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3xyz_Px_C1003_aa = I_TWOBODYOVERLAP_I4xyz_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H3xyz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3x2z_Px_C1003_aa = I_TWOBODYOVERLAP_I4x2z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H3x2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2x3y_Px_C1003_aa = I_TWOBODYOVERLAP_I3x3y_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H2x3y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2x2yz_Px_C1003_aa = I_TWOBODYOVERLAP_I3x2yz_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H2x2yz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2xy2z_Px_C1003_aa = I_TWOBODYOVERLAP_I3xy2z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H2xy2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2x3z_Px_C1003_aa = I_TWOBODYOVERLAP_I3x3z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H2x3z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx4y_Px_C1003_aa = I_TWOBODYOVERLAP_I2x4y_S_C1003_aa+ABX*I_TWOBODYOVERLAP_Hx4y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx3yz_Px_C1003_aa = I_TWOBODYOVERLAP_I2x3yz_S_C1003_aa+ABX*I_TWOBODYOVERLAP_Hx3yz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx2y2z_Px_C1003_aa = I_TWOBODYOVERLAP_I2x2y2z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_Hx2y2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hxy3z_Px_C1003_aa = I_TWOBODYOVERLAP_I2xy3z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_Hxy3z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx4z_Px_C1003_aa = I_TWOBODYOVERLAP_I2x4z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_Hx4z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H5y_Px_C1003_aa = I_TWOBODYOVERLAP_Ix5y_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H5y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H4yz_Px_C1003_aa = I_TWOBODYOVERLAP_Ix4yz_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H4yz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3y2z_Px_C1003_aa = I_TWOBODYOVERLAP_Ix3y2z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H3y2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2y3z_Px_C1003_aa = I_TWOBODYOVERLAP_Ix2y3z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H2y3z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hy4z_Px_C1003_aa = I_TWOBODYOVERLAP_Ixy4z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_Hy4z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H5z_Px_C1003_aa = I_TWOBODYOVERLAP_Ix5z_S_C1003_aa+ABX*I_TWOBODYOVERLAP_H5z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H5x_Py_C1003_aa = I_TWOBODYOVERLAP_I5xy_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H5x_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H4xy_Py_C1003_aa = I_TWOBODYOVERLAP_I4x2y_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H4xy_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H4xz_Py_C1003_aa = I_TWOBODYOVERLAP_I4xyz_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H4xz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3x2y_Py_C1003_aa = I_TWOBODYOVERLAP_I3x3y_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H3x2y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3xyz_Py_C1003_aa = I_TWOBODYOVERLAP_I3x2yz_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H3xyz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3x2z_Py_C1003_aa = I_TWOBODYOVERLAP_I3xy2z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H3x2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2x3y_Py_C1003_aa = I_TWOBODYOVERLAP_I2x4y_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H2x3y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2x2yz_Py_C1003_aa = I_TWOBODYOVERLAP_I2x3yz_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H2x2yz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2xy2z_Py_C1003_aa = I_TWOBODYOVERLAP_I2x2y2z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H2xy2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2x3z_Py_C1003_aa = I_TWOBODYOVERLAP_I2xy3z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H2x3z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx4y_Py_C1003_aa = I_TWOBODYOVERLAP_Ix5y_S_C1003_aa+ABY*I_TWOBODYOVERLAP_Hx4y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx3yz_Py_C1003_aa = I_TWOBODYOVERLAP_Ix4yz_S_C1003_aa+ABY*I_TWOBODYOVERLAP_Hx3yz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx2y2z_Py_C1003_aa = I_TWOBODYOVERLAP_Ix3y2z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_Hx2y2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hxy3z_Py_C1003_aa = I_TWOBODYOVERLAP_Ix2y3z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_Hxy3z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx4z_Py_C1003_aa = I_TWOBODYOVERLAP_Ixy4z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_Hx4z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H5y_Py_C1003_aa = I_TWOBODYOVERLAP_I6y_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H5y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H4yz_Py_C1003_aa = I_TWOBODYOVERLAP_I5yz_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H4yz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3y2z_Py_C1003_aa = I_TWOBODYOVERLAP_I4y2z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H3y2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2y3z_Py_C1003_aa = I_TWOBODYOVERLAP_I3y3z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H2y3z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hy4z_Py_C1003_aa = I_TWOBODYOVERLAP_I2y4z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_Hy4z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H5z_Py_C1003_aa = I_TWOBODYOVERLAP_Iy5z_S_C1003_aa+ABY*I_TWOBODYOVERLAP_H5z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H5x_Pz_C1003_aa = I_TWOBODYOVERLAP_I5xz_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H5x_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H4xy_Pz_C1003_aa = I_TWOBODYOVERLAP_I4xyz_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H4xy_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H4xz_Pz_C1003_aa = I_TWOBODYOVERLAP_I4x2z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H4xz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3x2y_Pz_C1003_aa = I_TWOBODYOVERLAP_I3x2yz_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H3x2y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3xyz_Pz_C1003_aa = I_TWOBODYOVERLAP_I3xy2z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H3xyz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3x2z_Pz_C1003_aa = I_TWOBODYOVERLAP_I3x3z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H3x2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2x3y_Pz_C1003_aa = I_TWOBODYOVERLAP_I2x3yz_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H2x3y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2x2yz_Pz_C1003_aa = I_TWOBODYOVERLAP_I2x2y2z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H2x2yz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2xy2z_Pz_C1003_aa = I_TWOBODYOVERLAP_I2xy3z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H2xy2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2x3z_Pz_C1003_aa = I_TWOBODYOVERLAP_I2x4z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H2x3z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx4y_Pz_C1003_aa = I_TWOBODYOVERLAP_Ix4yz_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_Hx4y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx3yz_Pz_C1003_aa = I_TWOBODYOVERLAP_Ix3y2z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_Hx3yz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx2y2z_Pz_C1003_aa = I_TWOBODYOVERLAP_Ix2y3z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_Hx2y2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hxy3z_Pz_C1003_aa = I_TWOBODYOVERLAP_Ixy4z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_Hxy3z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hx4z_Pz_C1003_aa = I_TWOBODYOVERLAP_Ix5z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_Hx4z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H5y_Pz_C1003_aa = I_TWOBODYOVERLAP_I5yz_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H5y_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H4yz_Pz_C1003_aa = I_TWOBODYOVERLAP_I4y2z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H4yz_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H3y2z_Pz_C1003_aa = I_TWOBODYOVERLAP_I3y3z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H3y2z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H2y3z_Pz_C1003_aa = I_TWOBODYOVERLAP_I2y4z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H2y3z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_Hy4z_Pz_C1003_aa = I_TWOBODYOVERLAP_Iy5z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_Hy4z_S_C1003_aa;
  Double I_TWOBODYOVERLAP_H5z_Pz_C1003_aa = I_TWOBODYOVERLAP_I6z_S_C1003_aa+ABZ*I_TWOBODYOVERLAP_H5z_S_C1003_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C3_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C3
   ************************************************************/
  abcd[0] = 4.0E0*I_TWOBODYOVERLAP_H5x_S_C3_aa-2.0E0*3*I_TWOBODYOVERLAP_F3x_S_C3_a-2.0E0*4*I_TWOBODYOVERLAP_F3x_S_C3_a+3*2*I_TWOBODYOVERLAP_Px_S_C3;
  abcd[1] = 4.0E0*I_TWOBODYOVERLAP_H4xy_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xy_S_C3_a-2.0E0*3*I_TWOBODYOVERLAP_F2xy_S_C3_a+2*1*I_TWOBODYOVERLAP_Py_S_C3;
  abcd[2] = 4.0E0*I_TWOBODYOVERLAP_H4xz_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xz_S_C3_a-2.0E0*3*I_TWOBODYOVERLAP_F2xz_S_C3_a+2*1*I_TWOBODYOVERLAP_Pz_S_C3;
  abcd[3] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_S_C3_a;
  abcd[4] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_S_C3_a;
  abcd[5] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_S_C3_a;
  abcd[6] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_S_C3_a;
  abcd[7] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_S_C3_a;
  abcd[8] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_S_C3_a;
  abcd[9] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F3z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1003_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1003
   ************************************************************/
  abcd[10] = 4.0E0*I_TWOBODYOVERLAP_H5x_Px_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F3x_Px_C1003_a-2.0E0*4*I_TWOBODYOVERLAP_F3x_Px_C1003_a+3*2*I_TWOBODYOVERLAP_Px_Px_C1003;
  abcd[11] = 4.0E0*I_TWOBODYOVERLAP_H4xy_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xy_Px_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_F2xy_Px_C1003_a+2*1*I_TWOBODYOVERLAP_Py_Px_C1003;
  abcd[12] = 4.0E0*I_TWOBODYOVERLAP_H4xz_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xz_Px_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_F2xz_Px_C1003_a+2*1*I_TWOBODYOVERLAP_Pz_Px_C1003;
  abcd[13] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_Px_C1003_a;
  abcd[14] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a;
  abcd[15] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_Px_C1003_a;
  abcd[16] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_Px_C1003_a;
  abcd[17] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Px_C1003_a;
  abcd[18] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Px_C1003_a;
  abcd[19] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3z_Px_C1003_a;
  abcd[20] = 4.0E0*I_TWOBODYOVERLAP_H5x_Py_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F3x_Py_C1003_a-2.0E0*4*I_TWOBODYOVERLAP_F3x_Py_C1003_a+3*2*I_TWOBODYOVERLAP_Px_Py_C1003;
  abcd[21] = 4.0E0*I_TWOBODYOVERLAP_H4xy_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xy_Py_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_F2xy_Py_C1003_a+2*1*I_TWOBODYOVERLAP_Py_Py_C1003;
  abcd[22] = 4.0E0*I_TWOBODYOVERLAP_H4xz_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xz_Py_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_F2xz_Py_C1003_a+2*1*I_TWOBODYOVERLAP_Pz_Py_C1003;
  abcd[23] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_Py_C1003_a;
  abcd[24] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a;
  abcd[25] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_Py_C1003_a;
  abcd[26] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_Py_C1003_a;
  abcd[27] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Py_C1003_a;
  abcd[28] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Py_C1003_a;
  abcd[29] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3z_Py_C1003_a;
  abcd[30] = 4.0E0*I_TWOBODYOVERLAP_H5x_Pz_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F3x_Pz_C1003_a-2.0E0*4*I_TWOBODYOVERLAP_F3x_Pz_C1003_a+3*2*I_TWOBODYOVERLAP_Px_Pz_C1003;
  abcd[31] = 4.0E0*I_TWOBODYOVERLAP_H4xy_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xy_Pz_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_F2xy_Pz_C1003_a+2*1*I_TWOBODYOVERLAP_Py_Pz_C1003;
  abcd[32] = 4.0E0*I_TWOBODYOVERLAP_H4xz_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xz_Pz_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_F2xz_Pz_C1003_a+2*1*I_TWOBODYOVERLAP_Pz_Pz_C1003;
  abcd[33] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a;
  abcd[34] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a;
  abcd[35] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a;
  abcd[36] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_Pz_C1003_a;
  abcd[37] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Pz_C1003_a;
  abcd[38] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a;
  abcd[39] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C3_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C3
   ************************************************************/
  abcd[40] = 4.0E0*I_TWOBODYOVERLAP_H4xy_S_C3_aa-2.0E0*3*I_TWOBODYOVERLAP_F2xy_S_C3_a;
  abcd[41] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_S_C3_a+2*1*I_TWOBODYOVERLAP_Px_S_C3;
  abcd[42] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_S_C3_a;
  abcd[43] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xy_S_C3_a-2.0E0*1*I_TWOBODYOVERLAP_F3y_S_C3_a+2*I_TWOBODYOVERLAP_Py_S_C3;
  abcd[44] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_S_C3_a-2.0E0*1*I_TWOBODYOVERLAP_F2yz_S_C3_a+1*I_TWOBODYOVERLAP_Pz_S_C3;
  abcd[45] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_S_C3_a;
  abcd[46] = 4.0E0*I_TWOBODYOVERLAP_Hx4y_S_C3_aa-2.0E0*3*I_TWOBODYOVERLAP_Fx2y_S_C3_a;
  abcd[47] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_S_C3_a;
  abcd[48] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_S_C3_a;
  abcd[49] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_S_C3_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1003_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1003
   ************************************************************/
  abcd[50] = 4.0E0*I_TWOBODYOVERLAP_H4xy_Px_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F2xy_Px_C1003_a;
  abcd[51] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_Px_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Px_C1003;
  abcd[52] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a;
  abcd[53] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xy_Px_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F3y_Px_C1003_a+2*I_TWOBODYOVERLAP_Py_Px_C1003;
  abcd[54] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Px_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Px_C1003_a+1*I_TWOBODYOVERLAP_Pz_Px_C1003;
  abcd[55] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Px_C1003_a;
  abcd[56] = 4.0E0*I_TWOBODYOVERLAP_Hx4y_Px_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_Fx2y_Px_C1003_a;
  abcd[57] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a;
  abcd[58] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Px_C1003_a;
  abcd[59] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Px_C1003_aa;
  abcd[60] = 4.0E0*I_TWOBODYOVERLAP_H4xy_Py_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F2xy_Py_C1003_a;
  abcd[61] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_Py_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Py_C1003;
  abcd[62] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a;
  abcd[63] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xy_Py_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F3y_Py_C1003_a+2*I_TWOBODYOVERLAP_Py_Py_C1003;
  abcd[64] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Py_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Py_C1003_a+1*I_TWOBODYOVERLAP_Pz_Py_C1003;
  abcd[65] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Py_C1003_a;
  abcd[66] = 4.0E0*I_TWOBODYOVERLAP_Hx4y_Py_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_Fx2y_Py_C1003_a;
  abcd[67] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a;
  abcd[68] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Py_C1003_a;
  abcd[69] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Py_C1003_aa;
  abcd[70] = 4.0E0*I_TWOBODYOVERLAP_H4xy_Pz_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F2xy_Pz_C1003_a;
  abcd[71] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Pz_C1003;
  abcd[72] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a;
  abcd[73] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xy_Pz_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F3y_Pz_C1003_a+2*I_TWOBODYOVERLAP_Py_Pz_C1003;
  abcd[74] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Pz_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Pz_C1003_a+1*I_TWOBODYOVERLAP_Pz_Pz_C1003;
  abcd[75] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a;
  abcd[76] = 4.0E0*I_TWOBODYOVERLAP_Hx4y_Pz_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a;
  abcd[77] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a;
  abcd[78] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a;
  abcd[79] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Pz_C1003_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C3_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C3
   ************************************************************/
  abcd[80] = 4.0E0*I_TWOBODYOVERLAP_H4xz_S_C3_aa-2.0E0*3*I_TWOBODYOVERLAP_F2xz_S_C3_a;
  abcd[81] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_S_C3_a;
  abcd[82] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_S_C3_a+2*1*I_TWOBODYOVERLAP_Px_S_C3;
  abcd[83] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_S_C3_a;
  abcd[84] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_S_C3_a-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_S_C3_a+1*I_TWOBODYOVERLAP_Py_S_C3;
  abcd[85] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xz_S_C3_a-2.0E0*1*I_TWOBODYOVERLAP_F3z_S_C3_a+2*I_TWOBODYOVERLAP_Pz_S_C3;
  abcd[86] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_S_C3_aa;
  abcd[87] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_S_C3_a;
  abcd[88] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_S_C3_a;
  abcd[89] = 4.0E0*I_TWOBODYOVERLAP_Hx4z_S_C3_aa-2.0E0*3*I_TWOBODYOVERLAP_Fx2z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1003_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1003
   ************************************************************/
  abcd[90] = 4.0E0*I_TWOBODYOVERLAP_H4xz_Px_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F2xz_Px_C1003_a;
  abcd[91] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a;
  abcd[92] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_Px_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Px_C1003;
  abcd[93] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Px_C1003_a;
  abcd[94] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Px_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Px_C1003_a+1*I_TWOBODYOVERLAP_Py_Px_C1003;
  abcd[95] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xz_Px_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F3z_Px_C1003_a+2*I_TWOBODYOVERLAP_Pz_Px_C1003;
  abcd[96] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Px_C1003_aa;
  abcd[97] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Px_C1003_a;
  abcd[98] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a;
  abcd[99] = 4.0E0*I_TWOBODYOVERLAP_Hx4z_Px_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_Fx2z_Px_C1003_a;
  abcd[100] = 4.0E0*I_TWOBODYOVERLAP_H4xz_Py_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F2xz_Py_C1003_a;
  abcd[101] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a;
  abcd[102] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_Py_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Py_C1003;
  abcd[103] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Py_C1003_a;
  abcd[104] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Py_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Py_C1003_a+1*I_TWOBODYOVERLAP_Py_Py_C1003;
  abcd[105] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xz_Py_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F3z_Py_C1003_a+2*I_TWOBODYOVERLAP_Pz_Py_C1003;
  abcd[106] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Py_C1003_aa;
  abcd[107] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Py_C1003_a;
  abcd[108] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a;
  abcd[109] = 4.0E0*I_TWOBODYOVERLAP_Hx4z_Py_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_Fx2z_Py_C1003_a;
  abcd[110] = 4.0E0*I_TWOBODYOVERLAP_H4xz_Pz_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F2xz_Pz_C1003_a;
  abcd[111] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a;
  abcd[112] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Pz_C1003;
  abcd[113] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Pz_C1003_a;
  abcd[114] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Pz_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a+1*I_TWOBODYOVERLAP_Py_Pz_C1003;
  abcd[115] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2xz_Pz_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F3z_Pz_C1003_a+2*I_TWOBODYOVERLAP_Pz_Pz_C1003;
  abcd[116] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Pz_C1003_aa;
  abcd[117] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a;
  abcd[118] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a;
  abcd[119] = 4.0E0*I_TWOBODYOVERLAP_Hx4z_Pz_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C3_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C3
   ************************************************************/
  abcd[120] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_S_C3_a;
  abcd[121] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_F2xy_S_C3_a;
  abcd[122] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_S_C3_a;
  abcd[123] = 4.0E0*I_TWOBODYOVERLAP_Hx4y_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_S_C3_a-2.0E0*3*I_TWOBODYOVERLAP_Fx2y_S_C3_a+2*1*I_TWOBODYOVERLAP_Px_S_C3;
  abcd[124] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_S_C3_a;
  abcd[125] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_S_C3_a;
  abcd[126] = 4.0E0*I_TWOBODYOVERLAP_H5y_S_C3_aa-2.0E0*3*I_TWOBODYOVERLAP_F3y_S_C3_a-2.0E0*4*I_TWOBODYOVERLAP_F3y_S_C3_a+3*2*I_TWOBODYOVERLAP_Py_S_C3;
  abcd[127] = 4.0E0*I_TWOBODYOVERLAP_H4yz_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_F2yz_S_C3_a-2.0E0*3*I_TWOBODYOVERLAP_F2yz_S_C3_a+2*1*I_TWOBODYOVERLAP_Pz_S_C3;
  abcd[128] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_S_C3_a;
  abcd[129] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F3z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1003_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1003
   ************************************************************/
  abcd[130] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Px_C1003_a;
  abcd[131] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_F2xy_Px_C1003_a;
  abcd[132] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Px_C1003_a;
  abcd[133] = 4.0E0*I_TWOBODYOVERLAP_Hx4y_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_Px_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_Fx2y_Px_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Px_C1003;
  abcd[134] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a;
  abcd[135] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Px_C1003_a;
  abcd[136] = 4.0E0*I_TWOBODYOVERLAP_H5y_Px_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F3y_Px_C1003_a-2.0E0*4*I_TWOBODYOVERLAP_F3y_Px_C1003_a+3*2*I_TWOBODYOVERLAP_Py_Px_C1003;
  abcd[137] = 4.0E0*I_TWOBODYOVERLAP_H4yz_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2yz_Px_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_F2yz_Px_C1003_a+2*1*I_TWOBODYOVERLAP_Pz_Px_C1003;
  abcd[138] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_Px_C1003_a;
  abcd[139] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3z_Px_C1003_a;
  abcd[140] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Py_C1003_a;
  abcd[141] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_F2xy_Py_C1003_a;
  abcd[142] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Py_C1003_a;
  abcd[143] = 4.0E0*I_TWOBODYOVERLAP_Hx4y_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_Py_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_Fx2y_Py_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Py_C1003;
  abcd[144] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a;
  abcd[145] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Py_C1003_a;
  abcd[146] = 4.0E0*I_TWOBODYOVERLAP_H5y_Py_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F3y_Py_C1003_a-2.0E0*4*I_TWOBODYOVERLAP_F3y_Py_C1003_a+3*2*I_TWOBODYOVERLAP_Py_Py_C1003;
  abcd[147] = 4.0E0*I_TWOBODYOVERLAP_H4yz_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2yz_Py_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_F2yz_Py_C1003_a+2*1*I_TWOBODYOVERLAP_Pz_Py_C1003;
  abcd[148] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_Py_C1003_a;
  abcd[149] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3z_Py_C1003_a;
  abcd[150] = 4.0E0*I_TWOBODYOVERLAP_H3x2y_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Pz_C1003_a;
  abcd[151] = 4.0E0*I_TWOBODYOVERLAP_H2x3y_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_F2xy_Pz_C1003_a;
  abcd[152] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Pz_C1003_a;
  abcd[153] = 4.0E0*I_TWOBODYOVERLAP_Hx4y_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Pz_C1003;
  abcd[154] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a;
  abcd[155] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a;
  abcd[156] = 4.0E0*I_TWOBODYOVERLAP_H5y_Pz_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F3y_Pz_C1003_a-2.0E0*4*I_TWOBODYOVERLAP_F3y_Pz_C1003_a+3*2*I_TWOBODYOVERLAP_Py_Pz_C1003;
  abcd[157] = 4.0E0*I_TWOBODYOVERLAP_H4yz_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2yz_Pz_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_F2yz_Pz_C1003_a+2*1*I_TWOBODYOVERLAP_Pz_Pz_C1003;
  abcd[158] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a;
  abcd[159] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C3_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C3
   ************************************************************/
  abcd[160] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_S_C3_aa;
  abcd[161] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_S_C3_a;
  abcd[162] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_S_C3_a;
  abcd[163] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_S_C3_a;
  abcd[164] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_S_C3_a-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_S_C3_a+1*I_TWOBODYOVERLAP_Px_S_C3;
  abcd[165] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_S_C3_a;
  abcd[166] = 4.0E0*I_TWOBODYOVERLAP_H4yz_S_C3_aa-2.0E0*3*I_TWOBODYOVERLAP_F2yz_S_C3_a;
  abcd[167] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_S_C3_a+2*1*I_TWOBODYOVERLAP_Py_S_C3;
  abcd[168] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_F2yz_S_C3_a-2.0E0*1*I_TWOBODYOVERLAP_F3z_S_C3_a+2*I_TWOBODYOVERLAP_Pz_S_C3;
  abcd[169] = 4.0E0*I_TWOBODYOVERLAP_Hy4z_S_C3_aa-2.0E0*3*I_TWOBODYOVERLAP_Fy2z_S_C3_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1003_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1003
   ************************************************************/
  abcd[170] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Px_C1003_aa;
  abcd[171] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Px_C1003_a;
  abcd[172] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Px_C1003_a;
  abcd[173] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a;
  abcd[174] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Px_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Px_C1003_a+1*I_TWOBODYOVERLAP_Px_Px_C1003;
  abcd[175] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a;
  abcd[176] = 4.0E0*I_TWOBODYOVERLAP_H4yz_Px_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F2yz_Px_C1003_a;
  abcd[177] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_Px_C1003_a+2*1*I_TWOBODYOVERLAP_Py_Px_C1003;
  abcd[178] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2yz_Px_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F3z_Px_C1003_a+2*I_TWOBODYOVERLAP_Pz_Px_C1003;
  abcd[179] = 4.0E0*I_TWOBODYOVERLAP_Hy4z_Px_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_Fy2z_Px_C1003_a;
  abcd[180] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Py_C1003_aa;
  abcd[181] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Py_C1003_a;
  abcd[182] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Py_C1003_a;
  abcd[183] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a;
  abcd[184] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Py_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Py_C1003_a+1*I_TWOBODYOVERLAP_Px_Py_C1003;
  abcd[185] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a;
  abcd[186] = 4.0E0*I_TWOBODYOVERLAP_H4yz_Py_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F2yz_Py_C1003_a;
  abcd[187] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_Py_C1003_a+2*1*I_TWOBODYOVERLAP_Py_Py_C1003;
  abcd[188] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2yz_Py_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F3z_Py_C1003_a+2*I_TWOBODYOVERLAP_Pz_Py_C1003;
  abcd[189] = 4.0E0*I_TWOBODYOVERLAP_Hy4z_Py_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_Fy2z_Py_C1003_a;
  abcd[190] = 4.0E0*I_TWOBODYOVERLAP_H3xyz_Pz_C1003_aa;
  abcd[191] = 4.0E0*I_TWOBODYOVERLAP_H2x2yz_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Pz_C1003_a;
  abcd[192] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Pz_C1003_a;
  abcd[193] = 4.0E0*I_TWOBODYOVERLAP_Hx3yz_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a;
  abcd[194] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a+1*I_TWOBODYOVERLAP_Px_Pz_C1003;
  abcd[195] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a;
  abcd[196] = 4.0E0*I_TWOBODYOVERLAP_H4yz_Pz_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F2yz_Pz_C1003_a;
  abcd[197] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a+2*1*I_TWOBODYOVERLAP_Py_Pz_C1003;
  abcd[198] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_F2yz_Pz_C1003_a-2.0E0*1*I_TWOBODYOVERLAP_F3z_Pz_C1003_a+2*I_TWOBODYOVERLAP_Pz_Pz_C1003;
  abcd[199] = 4.0E0*I_TWOBODYOVERLAP_Hy4z_Pz_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C3_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C3_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C3
   ************************************************************/
  abcd[200] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_S_C3_a;
  abcd[201] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_S_C3_a;
  abcd[202] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_F2xz_S_C3_a;
  abcd[203] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_S_C3_a;
  abcd[204] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_S_C3_a;
  abcd[205] = 4.0E0*I_TWOBODYOVERLAP_Hx4z_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_S_C3_a-2.0E0*3*I_TWOBODYOVERLAP_Fx2z_S_C3_a+2*1*I_TWOBODYOVERLAP_Px_S_C3;
  abcd[206] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_S_C3_a;
  abcd[207] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_S_C3_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_S_C3_a-2.0E0*2*I_TWOBODYOVERLAP_F2yz_S_C3_a;
  abcd[208] = 4.0E0*I_TWOBODYOVERLAP_Hy4z_S_C3_aa-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_S_C3_a-2.0E0*3*I_TWOBODYOVERLAP_Fy2z_S_C3_a+2*1*I_TWOBODYOVERLAP_Py_S_C3;
  abcd[209] = 4.0E0*I_TWOBODYOVERLAP_H5z_S_C3_aa-2.0E0*3*I_TWOBODYOVERLAP_F3z_S_C3_a-2.0E0*4*I_TWOBODYOVERLAP_F3z_S_C3_a+3*2*I_TWOBODYOVERLAP_Pz_S_C3;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1003_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1003_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1003
   ************************************************************/
  abcd[210] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Px_C1003_a;
  abcd[211] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Px_C1003_a;
  abcd[212] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_F2xz_Px_C1003_a;
  abcd[213] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Px_C1003_a;
  abcd[214] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Px_C1003_a;
  abcd[215] = 4.0E0*I_TWOBODYOVERLAP_Hx4z_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_Px_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_Fx2z_Px_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Px_C1003;
  abcd[216] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_Px_C1003_a;
  abcd[217] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_Px_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Px_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_F2yz_Px_C1003_a;
  abcd[218] = 4.0E0*I_TWOBODYOVERLAP_Hy4z_Px_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_Px_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_Fy2z_Px_C1003_a+2*1*I_TWOBODYOVERLAP_Py_Px_C1003;
  abcd[219] = 4.0E0*I_TWOBODYOVERLAP_H5z_Px_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F3z_Px_C1003_a-2.0E0*4*I_TWOBODYOVERLAP_F3z_Px_C1003_a+3*2*I_TWOBODYOVERLAP_Pz_Px_C1003;
  abcd[220] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Py_C1003_a;
  abcd[221] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Py_C1003_a;
  abcd[222] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_F2xz_Py_C1003_a;
  abcd[223] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Py_C1003_a;
  abcd[224] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Py_C1003_a;
  abcd[225] = 4.0E0*I_TWOBODYOVERLAP_Hx4z_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_Py_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_Fx2z_Py_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Py_C1003;
  abcd[226] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_Py_C1003_a;
  abcd[227] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_Py_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Py_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_F2yz_Py_C1003_a;
  abcd[228] = 4.0E0*I_TWOBODYOVERLAP_Hy4z_Py_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_Py_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_Fy2z_Py_C1003_a+2*1*I_TWOBODYOVERLAP_Py_Py_C1003;
  abcd[229] = 4.0E0*I_TWOBODYOVERLAP_H5z_Py_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F3z_Py_C1003_a-2.0E0*4*I_TWOBODYOVERLAP_F3z_Py_C1003_a+3*2*I_TWOBODYOVERLAP_Pz_Py_C1003;
  abcd[230] = 4.0E0*I_TWOBODYOVERLAP_H3x2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3x_Pz_C1003_a;
  abcd[231] = 4.0E0*I_TWOBODYOVERLAP_H2xy2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xy_Pz_C1003_a;
  abcd[232] = 4.0E0*I_TWOBODYOVERLAP_H2x3z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2xz_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_F2xz_Pz_C1003_a;
  abcd[233] = 4.0E0*I_TWOBODYOVERLAP_Hx2y2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fx2y_Pz_C1003_a;
  abcd[234] = 4.0E0*I_TWOBODYOVERLAP_Hxy3z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_Fxyz_Pz_C1003_a;
  abcd[235] = 4.0E0*I_TWOBODYOVERLAP_Hx4z_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_Fx2z_Pz_C1003_a+2*1*I_TWOBODYOVERLAP_Px_Pz_C1003;
  abcd[236] = 4.0E0*I_TWOBODYOVERLAP_H3y2z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F3y_Pz_C1003_a;
  abcd[237] = 4.0E0*I_TWOBODYOVERLAP_H2y3z_Pz_C1003_aa-2.0E0*1*I_TWOBODYOVERLAP_F2yz_Pz_C1003_a-2.0E0*2*I_TWOBODYOVERLAP_F2yz_Pz_C1003_a;
  abcd[238] = 4.0E0*I_TWOBODYOVERLAP_Hy4z_Pz_C1003_aa-2.0E0*2*I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a-2.0E0*3*I_TWOBODYOVERLAP_Fy2z_Pz_C1003_a+2*1*I_TWOBODYOVERLAP_Py_Pz_C1003;
  abcd[239] = 4.0E0*I_TWOBODYOVERLAP_H5z_Pz_C1003_aa-2.0E0*3*I_TWOBODYOVERLAP_F3z_Pz_C1003_a-2.0E0*4*I_TWOBODYOVERLAP_F3z_Pz_C1003_a+3*2*I_TWOBODYOVERLAP_Pz_Pz_C1003;
}
