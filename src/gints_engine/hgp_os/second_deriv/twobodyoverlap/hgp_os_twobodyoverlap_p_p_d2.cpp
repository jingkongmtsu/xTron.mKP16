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
// BRA1 as redundant position, total RHS integrals evaluated as: 0
// BRA2 as redundant position, total RHS integrals evaluated as: 0
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: NOT AVIALABLE
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

void hgp_os_twobodyoverlap_p_p_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_G4x_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_F3x_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xy_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Fxyz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_F3y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_F2yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Fy2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_F3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_D2x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_D2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Dyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_D2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Px_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Py_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Pz_S_a = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_G_S_aa_coefs = alpha*alpha;
    I_TWOBODYOVERLAP_G4x_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S_aa += SQ_TWOBODYOVERLAP_G_S_aa_coefs*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_F_S_aa_coefs = alpha*alpha;
    I_TWOBODYOVERLAP_F3x_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_F3x_S_vrr;
    I_TWOBODYOVERLAP_F2xy_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_F2xy_S_vrr;
    I_TWOBODYOVERLAP_F2xz_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_F2xz_S_vrr;
    I_TWOBODYOVERLAP_Fx2y_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    I_TWOBODYOVERLAP_Fxyz_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    I_TWOBODYOVERLAP_Fx2z_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    I_TWOBODYOVERLAP_F3y_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_F3y_S_vrr;
    I_TWOBODYOVERLAP_F2yz_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_F2yz_S_vrr;
    I_TWOBODYOVERLAP_Fy2z_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    I_TWOBODYOVERLAP_F3z_S_aa += SQ_TWOBODYOVERLAP_F_S_aa_coefs*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_D_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_D2x_S_a += SQ_TWOBODYOVERLAP_D_S_a_coefs*I_TWOBODYOVERLAP_D2x_S_vrr;
    I_TWOBODYOVERLAP_Dxy_S_a += SQ_TWOBODYOVERLAP_D_S_a_coefs*I_TWOBODYOVERLAP_Dxy_S_vrr;
    I_TWOBODYOVERLAP_Dxz_S_a += SQ_TWOBODYOVERLAP_D_S_a_coefs*I_TWOBODYOVERLAP_Dxz_S_vrr;
    I_TWOBODYOVERLAP_D2y_S_a += SQ_TWOBODYOVERLAP_D_S_a_coefs*I_TWOBODYOVERLAP_D2y_S_vrr;
    I_TWOBODYOVERLAP_Dyz_S_a += SQ_TWOBODYOVERLAP_D_S_a_coefs*I_TWOBODYOVERLAP_Dyz_S_vrr;
    I_TWOBODYOVERLAP_D2z_S_a += SQ_TWOBODYOVERLAP_D_S_a_coefs*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_P_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_Px_S_a += SQ_TWOBODYOVERLAP_P_S_a_coefs*I_TWOBODYOVERLAP_Px_S_vrr;
    I_TWOBODYOVERLAP_Py_S_a += SQ_TWOBODYOVERLAP_P_S_a_coefs*I_TWOBODYOVERLAP_Py_S_vrr;
    I_TWOBODYOVERLAP_Pz_S_a += SQ_TWOBODYOVERLAP_P_S_a_coefs*I_TWOBODYOVERLAP_Pz_S_vrr;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_Px_Px_a = I_TWOBODYOVERLAP_D2x_S_a+ABX*I_TWOBODYOVERLAP_Px_S_a;
  Double I_TWOBODYOVERLAP_Py_Px_a = I_TWOBODYOVERLAP_Dxy_S_a+ABX*I_TWOBODYOVERLAP_Py_S_a;
  Double I_TWOBODYOVERLAP_Pz_Px_a = I_TWOBODYOVERLAP_Dxz_S_a+ABX*I_TWOBODYOVERLAP_Pz_S_a;
  Double I_TWOBODYOVERLAP_Px_Py_a = I_TWOBODYOVERLAP_Dxy_S_a+ABY*I_TWOBODYOVERLAP_Px_S_a;
  Double I_TWOBODYOVERLAP_Py_Py_a = I_TWOBODYOVERLAP_D2y_S_a+ABY*I_TWOBODYOVERLAP_Py_S_a;
  Double I_TWOBODYOVERLAP_Pz_Py_a = I_TWOBODYOVERLAP_Dyz_S_a+ABY*I_TWOBODYOVERLAP_Pz_S_a;
  Double I_TWOBODYOVERLAP_Px_Pz_a = I_TWOBODYOVERLAP_Dxz_S_a+ABZ*I_TWOBODYOVERLAP_Px_S_a;
  Double I_TWOBODYOVERLAP_Py_Pz_a = I_TWOBODYOVERLAP_Dyz_S_a+ABZ*I_TWOBODYOVERLAP_Py_S_a;
  Double I_TWOBODYOVERLAP_Pz_Pz_a = I_TWOBODYOVERLAP_D2z_S_a+ABZ*I_TWOBODYOVERLAP_Pz_S_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_F3x_Px_aa = I_TWOBODYOVERLAP_G4x_S_aa+ABX*I_TWOBODYOVERLAP_F3x_S_aa;
  Double I_TWOBODYOVERLAP_F2xy_Px_aa = I_TWOBODYOVERLAP_G3xy_S_aa+ABX*I_TWOBODYOVERLAP_F2xy_S_aa;
  Double I_TWOBODYOVERLAP_F2xz_Px_aa = I_TWOBODYOVERLAP_G3xz_S_aa+ABX*I_TWOBODYOVERLAP_F2xz_S_aa;
  Double I_TWOBODYOVERLAP_Fx2y_Px_aa = I_TWOBODYOVERLAP_G2x2y_S_aa+ABX*I_TWOBODYOVERLAP_Fx2y_S_aa;
  Double I_TWOBODYOVERLAP_Fxyz_Px_aa = I_TWOBODYOVERLAP_G2xyz_S_aa+ABX*I_TWOBODYOVERLAP_Fxyz_S_aa;
  Double I_TWOBODYOVERLAP_Fx2z_Px_aa = I_TWOBODYOVERLAP_G2x2z_S_aa+ABX*I_TWOBODYOVERLAP_Fx2z_S_aa;
  Double I_TWOBODYOVERLAP_F3y_Px_aa = I_TWOBODYOVERLAP_Gx3y_S_aa+ABX*I_TWOBODYOVERLAP_F3y_S_aa;
  Double I_TWOBODYOVERLAP_F2yz_Px_aa = I_TWOBODYOVERLAP_Gx2yz_S_aa+ABX*I_TWOBODYOVERLAP_F2yz_S_aa;
  Double I_TWOBODYOVERLAP_Fy2z_Px_aa = I_TWOBODYOVERLAP_Gxy2z_S_aa+ABX*I_TWOBODYOVERLAP_Fy2z_S_aa;
  Double I_TWOBODYOVERLAP_F3z_Px_aa = I_TWOBODYOVERLAP_Gx3z_S_aa+ABX*I_TWOBODYOVERLAP_F3z_S_aa;
  Double I_TWOBODYOVERLAP_F3x_Py_aa = I_TWOBODYOVERLAP_G3xy_S_aa+ABY*I_TWOBODYOVERLAP_F3x_S_aa;
  Double I_TWOBODYOVERLAP_F2xy_Py_aa = I_TWOBODYOVERLAP_G2x2y_S_aa+ABY*I_TWOBODYOVERLAP_F2xy_S_aa;
  Double I_TWOBODYOVERLAP_F2xz_Py_aa = I_TWOBODYOVERLAP_G2xyz_S_aa+ABY*I_TWOBODYOVERLAP_F2xz_S_aa;
  Double I_TWOBODYOVERLAP_Fx2y_Py_aa = I_TWOBODYOVERLAP_Gx3y_S_aa+ABY*I_TWOBODYOVERLAP_Fx2y_S_aa;
  Double I_TWOBODYOVERLAP_Fxyz_Py_aa = I_TWOBODYOVERLAP_Gx2yz_S_aa+ABY*I_TWOBODYOVERLAP_Fxyz_S_aa;
  Double I_TWOBODYOVERLAP_Fx2z_Py_aa = I_TWOBODYOVERLAP_Gxy2z_S_aa+ABY*I_TWOBODYOVERLAP_Fx2z_S_aa;
  Double I_TWOBODYOVERLAP_F3y_Py_aa = I_TWOBODYOVERLAP_G4y_S_aa+ABY*I_TWOBODYOVERLAP_F3y_S_aa;
  Double I_TWOBODYOVERLAP_F2yz_Py_aa = I_TWOBODYOVERLAP_G3yz_S_aa+ABY*I_TWOBODYOVERLAP_F2yz_S_aa;
  Double I_TWOBODYOVERLAP_Fy2z_Py_aa = I_TWOBODYOVERLAP_G2y2z_S_aa+ABY*I_TWOBODYOVERLAP_Fy2z_S_aa;
  Double I_TWOBODYOVERLAP_F3z_Py_aa = I_TWOBODYOVERLAP_Gy3z_S_aa+ABY*I_TWOBODYOVERLAP_F3z_S_aa;
  Double I_TWOBODYOVERLAP_F3x_Pz_aa = I_TWOBODYOVERLAP_G3xz_S_aa+ABZ*I_TWOBODYOVERLAP_F3x_S_aa;
  Double I_TWOBODYOVERLAP_F2xy_Pz_aa = I_TWOBODYOVERLAP_G2xyz_S_aa+ABZ*I_TWOBODYOVERLAP_F2xy_S_aa;
  Double I_TWOBODYOVERLAP_F2xz_Pz_aa = I_TWOBODYOVERLAP_G2x2z_S_aa+ABZ*I_TWOBODYOVERLAP_F2xz_S_aa;
  Double I_TWOBODYOVERLAP_Fx2y_Pz_aa = I_TWOBODYOVERLAP_Gx2yz_S_aa+ABZ*I_TWOBODYOVERLAP_Fx2y_S_aa;
  Double I_TWOBODYOVERLAP_Fxyz_Pz_aa = I_TWOBODYOVERLAP_Gxy2z_S_aa+ABZ*I_TWOBODYOVERLAP_Fxyz_S_aa;
  Double I_TWOBODYOVERLAP_Fx2z_Pz_aa = I_TWOBODYOVERLAP_Gx3z_S_aa+ABZ*I_TWOBODYOVERLAP_Fx2z_S_aa;
  Double I_TWOBODYOVERLAP_F3y_Pz_aa = I_TWOBODYOVERLAP_G3yz_S_aa+ABZ*I_TWOBODYOVERLAP_F3y_S_aa;
  Double I_TWOBODYOVERLAP_F2yz_Pz_aa = I_TWOBODYOVERLAP_G2y2z_S_aa+ABZ*I_TWOBODYOVERLAP_F2yz_S_aa;
  Double I_TWOBODYOVERLAP_Fy2z_Pz_aa = I_TWOBODYOVERLAP_Gy3z_S_aa+ABZ*I_TWOBODYOVERLAP_Fy2z_S_aa;
  Double I_TWOBODYOVERLAP_F3z_Pz_aa = I_TWOBODYOVERLAP_G4z_S_aa+ABZ*I_TWOBODYOVERLAP_F3z_S_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_P_P_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   ************************************************************/
  abcd[0] = 4.0E0*I_TWOBODYOVERLAP_F3x_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Px_a-2.0E0*2*I_TWOBODYOVERLAP_Px_Px_a;
  abcd[1] = 4.0E0*I_TWOBODYOVERLAP_F2xy_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Px_a;
  abcd[2] = 4.0E0*I_TWOBODYOVERLAP_F2xz_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Px_a;
  abcd[3] = 4.0E0*I_TWOBODYOVERLAP_F3x_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Py_a-2.0E0*2*I_TWOBODYOVERLAP_Px_Py_a;
  abcd[4] = 4.0E0*I_TWOBODYOVERLAP_F2xy_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Py_a;
  abcd[5] = 4.0E0*I_TWOBODYOVERLAP_F2xz_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Py_a;
  abcd[6] = 4.0E0*I_TWOBODYOVERLAP_F3x_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Pz_a-2.0E0*2*I_TWOBODYOVERLAP_Px_Pz_a;
  abcd[7] = 4.0E0*I_TWOBODYOVERLAP_F2xy_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Pz_a;
  abcd[8] = 4.0E0*I_TWOBODYOVERLAP_F2xz_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Pz_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_P_P_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   ************************************************************/
  abcd[9] = 4.0E0*I_TWOBODYOVERLAP_F2xy_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Px_a;
  abcd[10] = 4.0E0*I_TWOBODYOVERLAP_Fx2y_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Px_a;
  abcd[11] = 4.0E0*I_TWOBODYOVERLAP_Fxyz_Px_aa;
  abcd[12] = 4.0E0*I_TWOBODYOVERLAP_F2xy_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Py_a;
  abcd[13] = 4.0E0*I_TWOBODYOVERLAP_Fx2y_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Py_a;
  abcd[14] = 4.0E0*I_TWOBODYOVERLAP_Fxyz_Py_aa;
  abcd[15] = 4.0E0*I_TWOBODYOVERLAP_F2xy_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Pz_a;
  abcd[16] = 4.0E0*I_TWOBODYOVERLAP_Fx2y_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Pz_a;
  abcd[17] = 4.0E0*I_TWOBODYOVERLAP_Fxyz_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_P_P_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   ************************************************************/
  abcd[18] = 4.0E0*I_TWOBODYOVERLAP_F2xz_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Px_a;
  abcd[19] = 4.0E0*I_TWOBODYOVERLAP_Fxyz_Px_aa;
  abcd[20] = 4.0E0*I_TWOBODYOVERLAP_Fx2z_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Px_a;
  abcd[21] = 4.0E0*I_TWOBODYOVERLAP_F2xz_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Py_a;
  abcd[22] = 4.0E0*I_TWOBODYOVERLAP_Fxyz_Py_aa;
  abcd[23] = 4.0E0*I_TWOBODYOVERLAP_Fx2z_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Py_a;
  abcd[24] = 4.0E0*I_TWOBODYOVERLAP_F2xz_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Pz_a;
  abcd[25] = 4.0E0*I_TWOBODYOVERLAP_Fxyz_Pz_aa;
  abcd[26] = 4.0E0*I_TWOBODYOVERLAP_Fx2z_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Pz_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_P_P_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   ************************************************************/
  abcd[27] = 4.0E0*I_TWOBODYOVERLAP_Fx2y_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Px_a;
  abcd[28] = 4.0E0*I_TWOBODYOVERLAP_F3y_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Px_a-2.0E0*2*I_TWOBODYOVERLAP_Py_Px_a;
  abcd[29] = 4.0E0*I_TWOBODYOVERLAP_F2yz_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Px_a;
  abcd[30] = 4.0E0*I_TWOBODYOVERLAP_Fx2y_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Py_a;
  abcd[31] = 4.0E0*I_TWOBODYOVERLAP_F3y_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Py_a-2.0E0*2*I_TWOBODYOVERLAP_Py_Py_a;
  abcd[32] = 4.0E0*I_TWOBODYOVERLAP_F2yz_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Py_a;
  abcd[33] = 4.0E0*I_TWOBODYOVERLAP_Fx2y_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Pz_a;
  abcd[34] = 4.0E0*I_TWOBODYOVERLAP_F3y_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Pz_a-2.0E0*2*I_TWOBODYOVERLAP_Py_Pz_a;
  abcd[35] = 4.0E0*I_TWOBODYOVERLAP_F2yz_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Pz_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_P_P_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   ************************************************************/
  abcd[36] = 4.0E0*I_TWOBODYOVERLAP_Fxyz_Px_aa;
  abcd[37] = 4.0E0*I_TWOBODYOVERLAP_F2yz_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Px_a;
  abcd[38] = 4.0E0*I_TWOBODYOVERLAP_Fy2z_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Px_a;
  abcd[39] = 4.0E0*I_TWOBODYOVERLAP_Fxyz_Py_aa;
  abcd[40] = 4.0E0*I_TWOBODYOVERLAP_F2yz_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Py_a;
  abcd[41] = 4.0E0*I_TWOBODYOVERLAP_Fy2z_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Py_a;
  abcd[42] = 4.0E0*I_TWOBODYOVERLAP_Fxyz_Pz_aa;
  abcd[43] = 4.0E0*I_TWOBODYOVERLAP_F2yz_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Pz_a;
  abcd[44] = 4.0E0*I_TWOBODYOVERLAP_Fy2z_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Pz_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_P_P_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P_a
   ************************************************************/
  abcd[45] = 4.0E0*I_TWOBODYOVERLAP_Fx2z_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Px_a;
  abcd[46] = 4.0E0*I_TWOBODYOVERLAP_Fy2z_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Px_a;
  abcd[47] = 4.0E0*I_TWOBODYOVERLAP_F3z_Px_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Px_a-2.0E0*2*I_TWOBODYOVERLAP_Pz_Px_a;
  abcd[48] = 4.0E0*I_TWOBODYOVERLAP_Fx2z_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Py_a;
  abcd[49] = 4.0E0*I_TWOBODYOVERLAP_Fy2z_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Py_a;
  abcd[50] = 4.0E0*I_TWOBODYOVERLAP_F3z_Py_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Py_a-2.0E0*2*I_TWOBODYOVERLAP_Pz_Py_a;
  abcd[51] = 4.0E0*I_TWOBODYOVERLAP_Fx2z_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Px_Pz_a;
  abcd[52] = 4.0E0*I_TWOBODYOVERLAP_Fy2z_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Py_Pz_a;
  abcd[53] = 4.0E0*I_TWOBODYOVERLAP_F3z_Pz_aa-2.0E0*1*I_TWOBODYOVERLAP_Pz_Pz_a-2.0E0*2*I_TWOBODYOVERLAP_Pz_Pz_a;
}
