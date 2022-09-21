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
// BRA1 as redundant position, total RHS integrals evaluated as: 165
// BRA2 as redundant position, total RHS integrals evaluated as: 123
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

void hgp_os_twobodyoverlap_sp_s_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_Px_S_C0_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Py_S_C0_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Pz_S_C0_a = 0.0E0;
  Double I_TWOBODYOVERLAP_D2x_S_C1_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxy_S_C1_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxz_S_C1_a = 0.0E0;
  Double I_TWOBODYOVERLAP_D2y_S_C1_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Dyz_S_C1_a = 0.0E0;
  Double I_TWOBODYOVERLAP_D2z_S_C1_a = 0.0E0;
  Double I_TWOBODYOVERLAP_S_S_C1 = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S_C0_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_P_S_C0_a_coefs = ic2*alpha;
    I_TWOBODYOVERLAP_Px_S_C0_a += SQ_TWOBODYOVERLAP_P_S_C0_a_coefs*I_TWOBODYOVERLAP_Px_S_vrr;
    I_TWOBODYOVERLAP_Py_S_C0_a += SQ_TWOBODYOVERLAP_P_S_C0_a_coefs*I_TWOBODYOVERLAP_Py_S_vrr;
    I_TWOBODYOVERLAP_Pz_S_C0_a += SQ_TWOBODYOVERLAP_P_S_C0_a_coefs*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S_C1_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_D_S_C1_a_coefs = ic2_1*alpha;
    I_TWOBODYOVERLAP_D2x_S_C1_a += SQ_TWOBODYOVERLAP_D_S_C1_a_coefs*I_TWOBODYOVERLAP_D2x_S_vrr;
    I_TWOBODYOVERLAP_Dxy_S_C1_a += SQ_TWOBODYOVERLAP_D_S_C1_a_coefs*I_TWOBODYOVERLAP_Dxy_S_vrr;
    I_TWOBODYOVERLAP_Dxz_S_C1_a += SQ_TWOBODYOVERLAP_D_S_C1_a_coefs*I_TWOBODYOVERLAP_Dxz_S_vrr;
    I_TWOBODYOVERLAP_D2y_S_C1_a += SQ_TWOBODYOVERLAP_D_S_C1_a_coefs*I_TWOBODYOVERLAP_D2y_S_vrr;
    I_TWOBODYOVERLAP_Dyz_S_C1_a += SQ_TWOBODYOVERLAP_D_S_C1_a_coefs*I_TWOBODYOVERLAP_Dyz_S_vrr;
    I_TWOBODYOVERLAP_D2z_S_C1_a += SQ_TWOBODYOVERLAP_D_S_C1_a_coefs*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_S_S_C1
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_S_S_C1_coefs = ic2_1;
    I_TWOBODYOVERLAP_S_S_C1 += SQ_TWOBODYOVERLAP_S_S_C1_coefs*I_TWOBODYOVERLAP_S_S_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_S_S_C0_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C0_a
   ************************************************************/
  abcd[0] = 2.0E0*I_TWOBODYOVERLAP_Px_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_P_S_C1_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S_C1_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S_C1
   ************************************************************/
  abcd[1] = 2.0E0*I_TWOBODYOVERLAP_D2x_S_C1_a-1*I_TWOBODYOVERLAP_S_S_C1;
  abcd[2] = 2.0E0*I_TWOBODYOVERLAP_Dxy_S_C1_a;
  abcd[3] = 2.0E0*I_TWOBODYOVERLAP_Dxz_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_S_S_C0_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C0_a
   ************************************************************/
  abcd[4] = 2.0E0*I_TWOBODYOVERLAP_Py_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_P_S_C1_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S_C1_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S_C1
   ************************************************************/
  abcd[5] = 2.0E0*I_TWOBODYOVERLAP_Dxy_S_C1_a;
  abcd[6] = 2.0E0*I_TWOBODYOVERLAP_D2y_S_C1_a-1*I_TWOBODYOVERLAP_S_S_C1;
  abcd[7] = 2.0E0*I_TWOBODYOVERLAP_Dyz_S_C1_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_S_S_C0_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S_C0_a
   ************************************************************/
  abcd[8] = 2.0E0*I_TWOBODYOVERLAP_Pz_S_C0_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_P_S_C1_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S_C1_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S_C1
   ************************************************************/
  abcd[9] = 2.0E0*I_TWOBODYOVERLAP_Dxz_S_C1_a;
  abcd[10] = 2.0E0*I_TWOBODYOVERLAP_Dyz_S_C1_a;
  abcd[11] = 2.0E0*I_TWOBODYOVERLAP_D2z_S_C1_a-1*I_TWOBODYOVERLAP_S_S_C1;
}
