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
// BRA1 as redundant position, total RHS integrals evaluated as: 422
// BRA2 as redundant position, total RHS integrals evaluated as: 324
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

void hgp_os_twobodyoverlap_d_sp_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_F3x_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xy_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xz_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2y_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fxyz_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2z_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3y_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2yz_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fy2z_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3z_S_C2_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Px_S_C2 = 0.0E0;
  Double I_TWOBODYOVERLAP_Py_S_C2 = 0.0E0;
  Double I_TWOBODYOVERLAP_Pz_S_C2 = 0.0E0;
  Double I_TWOBODYOVERLAP_G4x_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3x_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xy_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xz_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2y_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fxyz_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2z_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3y_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F2yz_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Fy2z_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_F3z_S_C1002_a = 0.0E0;
  Double I_TWOBODYOVERLAP_D2x_S_C1002 = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxy_S_C1002 = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxz_S_C1002 = 0.0E0;
  Double I_TWOBODYOVERLAP_D2y_S_C1002 = 0.0E0;
  Double I_TWOBODYOVERLAP_Dyz_S_C1002 = 0.0E0;
  Double I_TWOBODYOVERLAP_D2z_S_C1002 = 0.0E0;
  Double I_TWOBODYOVERLAP_Px_S_C1002 = 0.0E0;
  Double I_TWOBODYOVERLAP_Py_S_C1002 = 0.0E0;
  Double I_TWOBODYOVERLAP_Pz_S_C1002 = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C2_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_F_S_C2_a_coefs = ic2*alpha;
    I_TWOBODYOVERLAP_F3x_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_F3x_S_vrr;
    I_TWOBODYOVERLAP_F2xy_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_F2xy_S_vrr;
    I_TWOBODYOVERLAP_F2xz_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_F2xz_S_vrr;
    I_TWOBODYOVERLAP_Fx2y_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    I_TWOBODYOVERLAP_Fxyz_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    I_TWOBODYOVERLAP_Fx2z_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    I_TWOBODYOVERLAP_F3y_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_F3y_S_vrr;
    I_TWOBODYOVERLAP_F2yz_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_F2yz_S_vrr;
    I_TWOBODYOVERLAP_Fy2z_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    I_TWOBODYOVERLAP_F3z_S_C2_a += SQ_TWOBODYOVERLAP_F_S_C2_a_coefs*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S_C2
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_P_S_C2_coefs = ic2;
    I_TWOBODYOVERLAP_Px_S_C2 += SQ_TWOBODYOVERLAP_P_S_C2_coefs*I_TWOBODYOVERLAP_Px_S_vrr;
    I_TWOBODYOVERLAP_Py_S_C2 += SQ_TWOBODYOVERLAP_P_S_C2_coefs*I_TWOBODYOVERLAP_Py_S_vrr;
    I_TWOBODYOVERLAP_Pz_S_C2 += SQ_TWOBODYOVERLAP_P_S_C2_coefs*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S_C1002_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs = ic2_1*alpha;
    I_TWOBODYOVERLAP_G4x_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S_C1002_a += SQ_TWOBODYOVERLAP_G_S_C1002_a_coefs*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S_C1002_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs = ic2_1*alpha;
    I_TWOBODYOVERLAP_F3x_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_F3x_S_vrr;
    I_TWOBODYOVERLAP_F2xy_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_F2xy_S_vrr;
    I_TWOBODYOVERLAP_F2xz_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_F2xz_S_vrr;
    I_TWOBODYOVERLAP_Fx2y_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    I_TWOBODYOVERLAP_Fxyz_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    I_TWOBODYOVERLAP_Fx2z_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    I_TWOBODYOVERLAP_F3y_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_F3y_S_vrr;
    I_TWOBODYOVERLAP_F2yz_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_F2yz_S_vrr;
    I_TWOBODYOVERLAP_Fy2z_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    I_TWOBODYOVERLAP_F3z_S_C1002_a += SQ_TWOBODYOVERLAP_F_S_C1002_a_coefs*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S_C1002
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_D_S_C1002_coefs = ic2_1;
    I_TWOBODYOVERLAP_D2x_S_C1002 += SQ_TWOBODYOVERLAP_D_S_C1002_coefs*I_TWOBODYOVERLAP_D2x_S_vrr;
    I_TWOBODYOVERLAP_Dxy_S_C1002 += SQ_TWOBODYOVERLAP_D_S_C1002_coefs*I_TWOBODYOVERLAP_Dxy_S_vrr;
    I_TWOBODYOVERLAP_Dxz_S_C1002 += SQ_TWOBODYOVERLAP_D_S_C1002_coefs*I_TWOBODYOVERLAP_Dxz_S_vrr;
    I_TWOBODYOVERLAP_D2y_S_C1002 += SQ_TWOBODYOVERLAP_D_S_C1002_coefs*I_TWOBODYOVERLAP_D2y_S_vrr;
    I_TWOBODYOVERLAP_Dyz_S_C1002 += SQ_TWOBODYOVERLAP_D_S_C1002_coefs*I_TWOBODYOVERLAP_Dyz_S_vrr;
    I_TWOBODYOVERLAP_D2z_S_C1002 += SQ_TWOBODYOVERLAP_D_S_C1002_coefs*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S_C1002
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_P_S_C1002_coefs = ic2_1;
    I_TWOBODYOVERLAP_Px_S_C1002 += SQ_TWOBODYOVERLAP_P_S_C1002_coefs*I_TWOBODYOVERLAP_Px_S_vrr;
    I_TWOBODYOVERLAP_Py_S_C1002 += SQ_TWOBODYOVERLAP_P_S_C1002_coefs*I_TWOBODYOVERLAP_Py_S_vrr;
    I_TWOBODYOVERLAP_Pz_S_C1002 += SQ_TWOBODYOVERLAP_P_S_C1002_coefs*I_TWOBODYOVERLAP_Pz_S_vrr;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1002
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S_C1002
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C1002
   ************************************************************/
  Double I_TWOBODYOVERLAP_Px_Px_C1002 = I_TWOBODYOVERLAP_D2x_S_C1002+ABX*I_TWOBODYOVERLAP_Px_S_C1002;
  Double I_TWOBODYOVERLAP_Py_Px_C1002 = I_TWOBODYOVERLAP_Dxy_S_C1002+ABX*I_TWOBODYOVERLAP_Py_S_C1002;
  Double I_TWOBODYOVERLAP_Pz_Px_C1002 = I_TWOBODYOVERLAP_Dxz_S_C1002+ABX*I_TWOBODYOVERLAP_Pz_S_C1002;
  Double I_TWOBODYOVERLAP_Px_Py_C1002 = I_TWOBODYOVERLAP_Dxy_S_C1002+ABY*I_TWOBODYOVERLAP_Px_S_C1002;
  Double I_TWOBODYOVERLAP_Py_Py_C1002 = I_TWOBODYOVERLAP_D2y_S_C1002+ABY*I_TWOBODYOVERLAP_Py_S_C1002;
  Double I_TWOBODYOVERLAP_Pz_Py_C1002 = I_TWOBODYOVERLAP_Dyz_S_C1002+ABY*I_TWOBODYOVERLAP_Pz_S_C1002;
  Double I_TWOBODYOVERLAP_Px_Pz_C1002 = I_TWOBODYOVERLAP_Dxz_S_C1002+ABZ*I_TWOBODYOVERLAP_Px_S_C1002;
  Double I_TWOBODYOVERLAP_Py_Pz_C1002 = I_TWOBODYOVERLAP_Dyz_S_C1002+ABZ*I_TWOBODYOVERLAP_Py_S_C1002;
  Double I_TWOBODYOVERLAP_Pz_Pz_C1002 = I_TWOBODYOVERLAP_D2z_S_C1002+ABZ*I_TWOBODYOVERLAP_Pz_S_C1002;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1002_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_C1002_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C1002_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_F3x_Px_C1002_a = I_TWOBODYOVERLAP_G4x_S_C1002_a+ABX*I_TWOBODYOVERLAP_F3x_S_C1002_a;
  Double I_TWOBODYOVERLAP_F2xy_Px_C1002_a = I_TWOBODYOVERLAP_G3xy_S_C1002_a+ABX*I_TWOBODYOVERLAP_F2xy_S_C1002_a;
  Double I_TWOBODYOVERLAP_F2xz_Px_C1002_a = I_TWOBODYOVERLAP_G3xz_S_C1002_a+ABX*I_TWOBODYOVERLAP_F2xz_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fx2y_Px_C1002_a = I_TWOBODYOVERLAP_G2x2y_S_C1002_a+ABX*I_TWOBODYOVERLAP_Fx2y_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fxyz_Px_C1002_a = I_TWOBODYOVERLAP_G2xyz_S_C1002_a+ABX*I_TWOBODYOVERLAP_Fxyz_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fx2z_Px_C1002_a = I_TWOBODYOVERLAP_G2x2z_S_C1002_a+ABX*I_TWOBODYOVERLAP_Fx2z_S_C1002_a;
  Double I_TWOBODYOVERLAP_F3y_Px_C1002_a = I_TWOBODYOVERLAP_Gx3y_S_C1002_a+ABX*I_TWOBODYOVERLAP_F3y_S_C1002_a;
  Double I_TWOBODYOVERLAP_F2yz_Px_C1002_a = I_TWOBODYOVERLAP_Gx2yz_S_C1002_a+ABX*I_TWOBODYOVERLAP_F2yz_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fy2z_Px_C1002_a = I_TWOBODYOVERLAP_Gxy2z_S_C1002_a+ABX*I_TWOBODYOVERLAP_Fy2z_S_C1002_a;
  Double I_TWOBODYOVERLAP_F3z_Px_C1002_a = I_TWOBODYOVERLAP_Gx3z_S_C1002_a+ABX*I_TWOBODYOVERLAP_F3z_S_C1002_a;
  Double I_TWOBODYOVERLAP_F3x_Py_C1002_a = I_TWOBODYOVERLAP_G3xy_S_C1002_a+ABY*I_TWOBODYOVERLAP_F3x_S_C1002_a;
  Double I_TWOBODYOVERLAP_F2xy_Py_C1002_a = I_TWOBODYOVERLAP_G2x2y_S_C1002_a+ABY*I_TWOBODYOVERLAP_F2xy_S_C1002_a;
  Double I_TWOBODYOVERLAP_F2xz_Py_C1002_a = I_TWOBODYOVERLAP_G2xyz_S_C1002_a+ABY*I_TWOBODYOVERLAP_F2xz_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fx2y_Py_C1002_a = I_TWOBODYOVERLAP_Gx3y_S_C1002_a+ABY*I_TWOBODYOVERLAP_Fx2y_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fxyz_Py_C1002_a = I_TWOBODYOVERLAP_Gx2yz_S_C1002_a+ABY*I_TWOBODYOVERLAP_Fxyz_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fx2z_Py_C1002_a = I_TWOBODYOVERLAP_Gxy2z_S_C1002_a+ABY*I_TWOBODYOVERLAP_Fx2z_S_C1002_a;
  Double I_TWOBODYOVERLAP_F3y_Py_C1002_a = I_TWOBODYOVERLAP_G4y_S_C1002_a+ABY*I_TWOBODYOVERLAP_F3y_S_C1002_a;
  Double I_TWOBODYOVERLAP_F2yz_Py_C1002_a = I_TWOBODYOVERLAP_G3yz_S_C1002_a+ABY*I_TWOBODYOVERLAP_F2yz_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fy2z_Py_C1002_a = I_TWOBODYOVERLAP_G2y2z_S_C1002_a+ABY*I_TWOBODYOVERLAP_Fy2z_S_C1002_a;
  Double I_TWOBODYOVERLAP_F3z_Py_C1002_a = I_TWOBODYOVERLAP_Gy3z_S_C1002_a+ABY*I_TWOBODYOVERLAP_F3z_S_C1002_a;
  Double I_TWOBODYOVERLAP_F3x_Pz_C1002_a = I_TWOBODYOVERLAP_G3xz_S_C1002_a+ABZ*I_TWOBODYOVERLAP_F3x_S_C1002_a;
  Double I_TWOBODYOVERLAP_F2xy_Pz_C1002_a = I_TWOBODYOVERLAP_G2xyz_S_C1002_a+ABZ*I_TWOBODYOVERLAP_F2xy_S_C1002_a;
  Double I_TWOBODYOVERLAP_F2xz_Pz_C1002_a = I_TWOBODYOVERLAP_G2x2z_S_C1002_a+ABZ*I_TWOBODYOVERLAP_F2xz_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fx2y_Pz_C1002_a = I_TWOBODYOVERLAP_Gx2yz_S_C1002_a+ABZ*I_TWOBODYOVERLAP_Fx2y_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fxyz_Pz_C1002_a = I_TWOBODYOVERLAP_Gxy2z_S_C1002_a+ABZ*I_TWOBODYOVERLAP_Fxyz_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fx2z_Pz_C1002_a = I_TWOBODYOVERLAP_Gx3z_S_C1002_a+ABZ*I_TWOBODYOVERLAP_Fx2z_S_C1002_a;
  Double I_TWOBODYOVERLAP_F3y_Pz_C1002_a = I_TWOBODYOVERLAP_G3yz_S_C1002_a+ABZ*I_TWOBODYOVERLAP_F3y_S_C1002_a;
  Double I_TWOBODYOVERLAP_F2yz_Pz_C1002_a = I_TWOBODYOVERLAP_G2y2z_S_C1002_a+ABZ*I_TWOBODYOVERLAP_F2yz_S_C1002_a;
  Double I_TWOBODYOVERLAP_Fy2z_Pz_C1002_a = I_TWOBODYOVERLAP_Gy3z_S_C1002_a+ABZ*I_TWOBODYOVERLAP_Fy2z_S_C1002_a;
  Double I_TWOBODYOVERLAP_F3z_Pz_C1002_a = I_TWOBODYOVERLAP_G4z_S_C1002_a+ABZ*I_TWOBODYOVERLAP_F3z_S_C1002_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_D_S_C2_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C2_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C2
   ************************************************************/
  abcd[0] = 2.0E0*I_TWOBODYOVERLAP_F3x_S_C2_a-2*I_TWOBODYOVERLAP_Px_S_C2;
  abcd[1] = 2.0E0*I_TWOBODYOVERLAP_F2xy_S_C2_a-1*I_TWOBODYOVERLAP_Py_S_C2;
  abcd[2] = 2.0E0*I_TWOBODYOVERLAP_F2xz_S_C2_a-1*I_TWOBODYOVERLAP_Pz_S_C2;
  abcd[3] = 2.0E0*I_TWOBODYOVERLAP_Fx2y_S_C2_a;
  abcd[4] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_S_C2_a;
  abcd[5] = 2.0E0*I_TWOBODYOVERLAP_Fx2z_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_D_P_C1002_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1002_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1002
   ************************************************************/
  abcd[6] = 2.0E0*I_TWOBODYOVERLAP_F3x_Px_C1002_a-2*I_TWOBODYOVERLAP_Px_Px_C1002;
  abcd[7] = 2.0E0*I_TWOBODYOVERLAP_F2xy_Px_C1002_a-1*I_TWOBODYOVERLAP_Py_Px_C1002;
  abcd[8] = 2.0E0*I_TWOBODYOVERLAP_F2xz_Px_C1002_a-1*I_TWOBODYOVERLAP_Pz_Px_C1002;
  abcd[9] = 2.0E0*I_TWOBODYOVERLAP_Fx2y_Px_C1002_a;
  abcd[10] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_Px_C1002_a;
  abcd[11] = 2.0E0*I_TWOBODYOVERLAP_Fx2z_Px_C1002_a;
  abcd[12] = 2.0E0*I_TWOBODYOVERLAP_F3x_Py_C1002_a-2*I_TWOBODYOVERLAP_Px_Py_C1002;
  abcd[13] = 2.0E0*I_TWOBODYOVERLAP_F2xy_Py_C1002_a-1*I_TWOBODYOVERLAP_Py_Py_C1002;
  abcd[14] = 2.0E0*I_TWOBODYOVERLAP_F2xz_Py_C1002_a-1*I_TWOBODYOVERLAP_Pz_Py_C1002;
  abcd[15] = 2.0E0*I_TWOBODYOVERLAP_Fx2y_Py_C1002_a;
  abcd[16] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_Py_C1002_a;
  abcd[17] = 2.0E0*I_TWOBODYOVERLAP_Fx2z_Py_C1002_a;
  abcd[18] = 2.0E0*I_TWOBODYOVERLAP_F3x_Pz_C1002_a-2*I_TWOBODYOVERLAP_Px_Pz_C1002;
  abcd[19] = 2.0E0*I_TWOBODYOVERLAP_F2xy_Pz_C1002_a-1*I_TWOBODYOVERLAP_Py_Pz_C1002;
  abcd[20] = 2.0E0*I_TWOBODYOVERLAP_F2xz_Pz_C1002_a-1*I_TWOBODYOVERLAP_Pz_Pz_C1002;
  abcd[21] = 2.0E0*I_TWOBODYOVERLAP_Fx2y_Pz_C1002_a;
  abcd[22] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_Pz_C1002_a;
  abcd[23] = 2.0E0*I_TWOBODYOVERLAP_Fx2z_Pz_C1002_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_D_S_C2_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C2_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C2
   ************************************************************/
  abcd[24] = 2.0E0*I_TWOBODYOVERLAP_F2xy_S_C2_a;
  abcd[25] = 2.0E0*I_TWOBODYOVERLAP_Fx2y_S_C2_a-1*I_TWOBODYOVERLAP_Px_S_C2;
  abcd[26] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_S_C2_a;
  abcd[27] = 2.0E0*I_TWOBODYOVERLAP_F3y_S_C2_a-2*I_TWOBODYOVERLAP_Py_S_C2;
  abcd[28] = 2.0E0*I_TWOBODYOVERLAP_F2yz_S_C2_a-1*I_TWOBODYOVERLAP_Pz_S_C2;
  abcd[29] = 2.0E0*I_TWOBODYOVERLAP_Fy2z_S_C2_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_D_P_C1002_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1002_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1002
   ************************************************************/
  abcd[30] = 2.0E0*I_TWOBODYOVERLAP_F2xy_Px_C1002_a;
  abcd[31] = 2.0E0*I_TWOBODYOVERLAP_Fx2y_Px_C1002_a-1*I_TWOBODYOVERLAP_Px_Px_C1002;
  abcd[32] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_Px_C1002_a;
  abcd[33] = 2.0E0*I_TWOBODYOVERLAP_F3y_Px_C1002_a-2*I_TWOBODYOVERLAP_Py_Px_C1002;
  abcd[34] = 2.0E0*I_TWOBODYOVERLAP_F2yz_Px_C1002_a-1*I_TWOBODYOVERLAP_Pz_Px_C1002;
  abcd[35] = 2.0E0*I_TWOBODYOVERLAP_Fy2z_Px_C1002_a;
  abcd[36] = 2.0E0*I_TWOBODYOVERLAP_F2xy_Py_C1002_a;
  abcd[37] = 2.0E0*I_TWOBODYOVERLAP_Fx2y_Py_C1002_a-1*I_TWOBODYOVERLAP_Px_Py_C1002;
  abcd[38] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_Py_C1002_a;
  abcd[39] = 2.0E0*I_TWOBODYOVERLAP_F3y_Py_C1002_a-2*I_TWOBODYOVERLAP_Py_Py_C1002;
  abcd[40] = 2.0E0*I_TWOBODYOVERLAP_F2yz_Py_C1002_a-1*I_TWOBODYOVERLAP_Pz_Py_C1002;
  abcd[41] = 2.0E0*I_TWOBODYOVERLAP_Fy2z_Py_C1002_a;
  abcd[42] = 2.0E0*I_TWOBODYOVERLAP_F2xy_Pz_C1002_a;
  abcd[43] = 2.0E0*I_TWOBODYOVERLAP_Fx2y_Pz_C1002_a-1*I_TWOBODYOVERLAP_Px_Pz_C1002;
  abcd[44] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_Pz_C1002_a;
  abcd[45] = 2.0E0*I_TWOBODYOVERLAP_F3y_Pz_C1002_a-2*I_TWOBODYOVERLAP_Py_Pz_C1002;
  abcd[46] = 2.0E0*I_TWOBODYOVERLAP_F2yz_Pz_C1002_a-1*I_TWOBODYOVERLAP_Pz_Pz_C1002;
  abcd[47] = 2.0E0*I_TWOBODYOVERLAP_Fy2z_Pz_C1002_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_D_S_C2_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_C2_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C2
   ************************************************************/
  abcd[48] = 2.0E0*I_TWOBODYOVERLAP_F2xz_S_C2_a;
  abcd[49] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_S_C2_a;
  abcd[50] = 2.0E0*I_TWOBODYOVERLAP_Fx2z_S_C2_a-1*I_TWOBODYOVERLAP_Px_S_C2;
  abcd[51] = 2.0E0*I_TWOBODYOVERLAP_F2yz_S_C2_a;
  abcd[52] = 2.0E0*I_TWOBODYOVERLAP_Fy2z_S_C2_a-1*I_TWOBODYOVERLAP_Py_S_C2;
  abcd[53] = 2.0E0*I_TWOBODYOVERLAP_F3z_S_C2_a-2*I_TWOBODYOVERLAP_Pz_S_C2;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_D_P_C1002_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_C1002_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_C1002
   ************************************************************/
  abcd[54] = 2.0E0*I_TWOBODYOVERLAP_F2xz_Px_C1002_a;
  abcd[55] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_Px_C1002_a;
  abcd[56] = 2.0E0*I_TWOBODYOVERLAP_Fx2z_Px_C1002_a-1*I_TWOBODYOVERLAP_Px_Px_C1002;
  abcd[57] = 2.0E0*I_TWOBODYOVERLAP_F2yz_Px_C1002_a;
  abcd[58] = 2.0E0*I_TWOBODYOVERLAP_Fy2z_Px_C1002_a-1*I_TWOBODYOVERLAP_Py_Px_C1002;
  abcd[59] = 2.0E0*I_TWOBODYOVERLAP_F3z_Px_C1002_a-2*I_TWOBODYOVERLAP_Pz_Px_C1002;
  abcd[60] = 2.0E0*I_TWOBODYOVERLAP_F2xz_Py_C1002_a;
  abcd[61] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_Py_C1002_a;
  abcd[62] = 2.0E0*I_TWOBODYOVERLAP_Fx2z_Py_C1002_a-1*I_TWOBODYOVERLAP_Px_Py_C1002;
  abcd[63] = 2.0E0*I_TWOBODYOVERLAP_F2yz_Py_C1002_a;
  abcd[64] = 2.0E0*I_TWOBODYOVERLAP_Fy2z_Py_C1002_a-1*I_TWOBODYOVERLAP_Py_Py_C1002;
  abcd[65] = 2.0E0*I_TWOBODYOVERLAP_F3z_Py_C1002_a-2*I_TWOBODYOVERLAP_Pz_Py_C1002;
  abcd[66] = 2.0E0*I_TWOBODYOVERLAP_F2xz_Pz_C1002_a;
  abcd[67] = 2.0E0*I_TWOBODYOVERLAP_Fxyz_Pz_C1002_a;
  abcd[68] = 2.0E0*I_TWOBODYOVERLAP_Fx2z_Pz_C1002_a-1*I_TWOBODYOVERLAP_Px_Pz_C1002;
  abcd[69] = 2.0E0*I_TWOBODYOVERLAP_F2yz_Pz_C1002_a;
  abcd[70] = 2.0E0*I_TWOBODYOVERLAP_Fy2z_Pz_C1002_a-1*I_TWOBODYOVERLAP_Py_Pz_C1002;
  abcd[71] = 2.0E0*I_TWOBODYOVERLAP_F3z_Pz_C1002_a-2*I_TWOBODYOVERLAP_Pz_Pz_C1002;
}
