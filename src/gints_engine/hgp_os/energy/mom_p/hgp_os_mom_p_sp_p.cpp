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

void hgp_os_mom_p_sp_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_MOM_D2x_S_Px_C1001001 = 0.0E0;
  Double I_MOM_Dxy_S_Px_C1001001 = 0.0E0;
  Double I_MOM_Dxz_S_Px_C1001001 = 0.0E0;
  Double I_MOM_D2y_S_Px_C1001001 = 0.0E0;
  Double I_MOM_Dyz_S_Px_C1001001 = 0.0E0;
  Double I_MOM_D2z_S_Px_C1001001 = 0.0E0;
  Double I_MOM_D2x_S_Py_C1001001 = 0.0E0;
  Double I_MOM_Dxy_S_Py_C1001001 = 0.0E0;
  Double I_MOM_Dxz_S_Py_C1001001 = 0.0E0;
  Double I_MOM_D2y_S_Py_C1001001 = 0.0E0;
  Double I_MOM_Dyz_S_Py_C1001001 = 0.0E0;
  Double I_MOM_D2z_S_Py_C1001001 = 0.0E0;
  Double I_MOM_D2x_S_Pz_C1001001 = 0.0E0;
  Double I_MOM_Dxy_S_Pz_C1001001 = 0.0E0;
  Double I_MOM_Dxz_S_Pz_C1001001 = 0.0E0;
  Double I_MOM_D2y_S_Pz_C1001001 = 0.0E0;
  Double I_MOM_Dyz_S_Pz_C1001001 = 0.0E0;
  Double I_MOM_D2z_S_Pz_C1001001 = 0.0E0;
  Double I_MOM_Px_S_Px_C1001001 = 0.0E0;
  Double I_MOM_Py_S_Px_C1001001 = 0.0E0;
  Double I_MOM_Pz_S_Px_C1001001 = 0.0E0;
  Double I_MOM_Px_S_Py_C1001001 = 0.0E0;
  Double I_MOM_Py_S_Py_C1001001 = 0.0E0;
  Double I_MOM_Pz_S_Py_C1001001 = 0.0E0;
  Double I_MOM_Px_S_Pz_C1001001 = 0.0E0;
  Double I_MOM_Py_S_Pz_C1001001 = 0.0E0;
  Double I_MOM_Pz_S_Pz_C1001001 = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double ic2_1 = icoe[ip2+1*inp2];
    Double fbra  = ifac[ip2];
    Double onedz = iexp[ip2];
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PCX   = PX - C[0];
    Double PCY   = PY - C[1];
    Double PCZ   = PZ - C[2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    Double I_TWOBODYOVERLAP_S_S_vrr = fbra;
    if(fabs(I_TWOBODYOVERLAP_S_S_vrr)<THRESHOLD_MATH) continue;


    // 
    // now create the bottom integrals for momentum
    // 
    Double I_MOM_S_S_S_vrr   = I_TWOBODYOVERLAP_S_S_vrr;
    Double I_MOM_S_S_Px_vrr = PCX*I_MOM_S_S_S_vrr;
    Double I_MOM_S_S_Py_vrr = PCY*I_MOM_S_S_S_vrr;
    Double I_MOM_S_S_Pz_vrr = PCZ*I_MOM_S_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_P_S_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_MOM_S_S_S
     ************************************************************/
    Double I_MOM_Px_S_S_vrr = PAX*I_MOM_S_S_S_vrr;
    Double I_MOM_Py_S_S_vrr = PAY*I_MOM_S_S_S_vrr;
    Double I_MOM_Pz_S_S_vrr = PAZ*I_MOM_S_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_P_S_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_MOM_S_S_P
     * RHS shell quartet name: SQ_MOM_S_S_S
     ************************************************************/
    Double I_MOM_Px_S_Px_vrr = PAX*I_MOM_S_S_Px_vrr+oned2z*I_MOM_S_S_S_vrr;
    Double I_MOM_Py_S_Px_vrr = PAY*I_MOM_S_S_Px_vrr;
    Double I_MOM_Pz_S_Px_vrr = PAZ*I_MOM_S_S_Px_vrr;
    Double I_MOM_Px_S_Py_vrr = PAX*I_MOM_S_S_Py_vrr;
    Double I_MOM_Py_S_Py_vrr = PAY*I_MOM_S_S_Py_vrr+oned2z*I_MOM_S_S_S_vrr;
    Double I_MOM_Pz_S_Py_vrr = PAZ*I_MOM_S_S_Py_vrr;
    Double I_MOM_Px_S_Pz_vrr = PAX*I_MOM_S_S_Pz_vrr;
    Double I_MOM_Py_S_Pz_vrr = PAY*I_MOM_S_S_Pz_vrr;
    Double I_MOM_Pz_S_Pz_vrr = PAZ*I_MOM_S_S_Pz_vrr+oned2z*I_MOM_S_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_D_S_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_MOM_P_S_P
     * RHS shell quartet name: SQ_MOM_S_S_P
     * RHS shell quartet name: SQ_MOM_P_S_S
     ************************************************************/
    Double I_MOM_D2x_S_Px_vrr = PAX*I_MOM_Px_S_Px_vrr+oned2z*I_MOM_S_S_Px_vrr+oned2z*I_MOM_Px_S_S_vrr;
    Double I_MOM_Dxy_S_Px_vrr = PAY*I_MOM_Px_S_Px_vrr;
    Double I_MOM_Dxz_S_Px_vrr = PAZ*I_MOM_Px_S_Px_vrr;
    Double I_MOM_D2y_S_Px_vrr = PAY*I_MOM_Py_S_Px_vrr+oned2z*I_MOM_S_S_Px_vrr;
    Double I_MOM_Dyz_S_Px_vrr = PAZ*I_MOM_Py_S_Px_vrr;
    Double I_MOM_D2z_S_Px_vrr = PAZ*I_MOM_Pz_S_Px_vrr+oned2z*I_MOM_S_S_Px_vrr;
    Double I_MOM_D2x_S_Py_vrr = PAX*I_MOM_Px_S_Py_vrr+oned2z*I_MOM_S_S_Py_vrr;
    Double I_MOM_Dxy_S_Py_vrr = PAY*I_MOM_Px_S_Py_vrr+oned2z*I_MOM_Px_S_S_vrr;
    Double I_MOM_Dxz_S_Py_vrr = PAZ*I_MOM_Px_S_Py_vrr;
    Double I_MOM_D2y_S_Py_vrr = PAY*I_MOM_Py_S_Py_vrr+oned2z*I_MOM_S_S_Py_vrr+oned2z*I_MOM_Py_S_S_vrr;
    Double I_MOM_Dyz_S_Py_vrr = PAZ*I_MOM_Py_S_Py_vrr;
    Double I_MOM_D2z_S_Py_vrr = PAZ*I_MOM_Pz_S_Py_vrr+oned2z*I_MOM_S_S_Py_vrr;
    Double I_MOM_D2x_S_Pz_vrr = PAX*I_MOM_Px_S_Pz_vrr+oned2z*I_MOM_S_S_Pz_vrr;
    Double I_MOM_Dxy_S_Pz_vrr = PAY*I_MOM_Px_S_Pz_vrr;
    Double I_MOM_Dxz_S_Pz_vrr = PAZ*I_MOM_Px_S_Pz_vrr+oned2z*I_MOM_Px_S_S_vrr;
    Double I_MOM_D2y_S_Pz_vrr = PAY*I_MOM_Py_S_Pz_vrr+oned2z*I_MOM_S_S_Pz_vrr;
    Double I_MOM_Dyz_S_Pz_vrr = PAZ*I_MOM_Py_S_Pz_vrr+oned2z*I_MOM_Py_S_S_vrr;
    Double I_MOM_D2z_S_Pz_vrr = PAZ*I_MOM_Pz_S_Pz_vrr+oned2z*I_MOM_S_S_Pz_vrr+oned2z*I_MOM_Pz_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_P_S_P_C1000001
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_MOM_P_S_P_C1000001_coefs = ic2;
    abcd[0] += SQ_MOM_P_S_P_C1000001_coefs*I_MOM_Px_S_Px_vrr;
    abcd[1] += SQ_MOM_P_S_P_C1000001_coefs*I_MOM_Py_S_Px_vrr;
    abcd[2] += SQ_MOM_P_S_P_C1000001_coefs*I_MOM_Pz_S_Px_vrr;
    abcd[12] += SQ_MOM_P_S_P_C1000001_coefs*I_MOM_Px_S_Py_vrr;
    abcd[13] += SQ_MOM_P_S_P_C1000001_coefs*I_MOM_Py_S_Py_vrr;
    abcd[14] += SQ_MOM_P_S_P_C1000001_coefs*I_MOM_Pz_S_Py_vrr;
    abcd[24] += SQ_MOM_P_S_P_C1000001_coefs*I_MOM_Px_S_Pz_vrr;
    abcd[25] += SQ_MOM_P_S_P_C1000001_coefs*I_MOM_Py_S_Pz_vrr;
    abcd[26] += SQ_MOM_P_S_P_C1000001_coefs*I_MOM_Pz_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_D_S_P_C1001001
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_MOM_D_S_P_C1001001_coefs = ic2_1;
    I_MOM_D2x_S_Px_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_D2x_S_Px_vrr;
    I_MOM_Dxy_S_Px_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_Dxy_S_Px_vrr;
    I_MOM_Dxz_S_Px_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_Dxz_S_Px_vrr;
    I_MOM_D2y_S_Px_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_D2y_S_Px_vrr;
    I_MOM_Dyz_S_Px_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_Dyz_S_Px_vrr;
    I_MOM_D2z_S_Px_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_D2z_S_Px_vrr;
    I_MOM_D2x_S_Py_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_D2x_S_Py_vrr;
    I_MOM_Dxy_S_Py_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_Dxy_S_Py_vrr;
    I_MOM_Dxz_S_Py_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_Dxz_S_Py_vrr;
    I_MOM_D2y_S_Py_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_D2y_S_Py_vrr;
    I_MOM_Dyz_S_Py_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_Dyz_S_Py_vrr;
    I_MOM_D2z_S_Py_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_D2z_S_Py_vrr;
    I_MOM_D2x_S_Pz_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_D2x_S_Pz_vrr;
    I_MOM_Dxy_S_Pz_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_Dxy_S_Pz_vrr;
    I_MOM_Dxz_S_Pz_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_Dxz_S_Pz_vrr;
    I_MOM_D2y_S_Pz_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_D2y_S_Pz_vrr;
    I_MOM_Dyz_S_Pz_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_Dyz_S_Pz_vrr;
    I_MOM_D2z_S_Pz_C1001001 += SQ_MOM_D_S_P_C1001001_coefs*I_MOM_D2z_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_P_S_P_C1001001
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_MOM_P_S_P_C1001001_coefs = ic2_1;
    I_MOM_Px_S_Px_C1001001 += SQ_MOM_P_S_P_C1001001_coefs*I_MOM_Px_S_Px_vrr;
    I_MOM_Py_S_Px_C1001001 += SQ_MOM_P_S_P_C1001001_coefs*I_MOM_Py_S_Px_vrr;
    I_MOM_Pz_S_Px_C1001001 += SQ_MOM_P_S_P_C1001001_coefs*I_MOM_Pz_S_Px_vrr;
    I_MOM_Px_S_Py_C1001001 += SQ_MOM_P_S_P_C1001001_coefs*I_MOM_Px_S_Py_vrr;
    I_MOM_Py_S_Py_C1001001 += SQ_MOM_P_S_P_C1001001_coefs*I_MOM_Py_S_Py_vrr;
    I_MOM_Pz_S_Py_C1001001 += SQ_MOM_P_S_P_C1001001_coefs*I_MOM_Pz_S_Py_vrr;
    I_MOM_Px_S_Pz_C1001001 += SQ_MOM_P_S_P_C1001001_coefs*I_MOM_Px_S_Pz_vrr;
    I_MOM_Py_S_Pz_C1001001 += SQ_MOM_P_S_P_C1001001_coefs*I_MOM_Py_S_Pz_vrr;
    I_MOM_Pz_S_Pz_C1001001 += SQ_MOM_P_S_P_C1001001_coefs*I_MOM_Pz_S_Pz_vrr;
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
   * shell quartet name: SQ_MOM_P_P_P_C1001001
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_MOM_D_S_P_C1001001
   * RHS shell quartet name: SQ_MOM_P_S_P_C1001001
   ************************************************************/
  abcd[3] = I_MOM_D2x_S_Px_C1001001+ABX*I_MOM_Px_S_Px_C1001001;
  abcd[4] = I_MOM_Dxy_S_Px_C1001001+ABX*I_MOM_Py_S_Px_C1001001;
  abcd[5] = I_MOM_Dxz_S_Px_C1001001+ABX*I_MOM_Pz_S_Px_C1001001;
  abcd[6] = I_MOM_Dxy_S_Px_C1001001+ABY*I_MOM_Px_S_Px_C1001001;
  abcd[7] = I_MOM_D2y_S_Px_C1001001+ABY*I_MOM_Py_S_Px_C1001001;
  abcd[8] = I_MOM_Dyz_S_Px_C1001001+ABY*I_MOM_Pz_S_Px_C1001001;
  abcd[9] = I_MOM_Dxz_S_Px_C1001001+ABZ*I_MOM_Px_S_Px_C1001001;
  abcd[10] = I_MOM_Dyz_S_Px_C1001001+ABZ*I_MOM_Py_S_Px_C1001001;
  abcd[11] = I_MOM_D2z_S_Px_C1001001+ABZ*I_MOM_Pz_S_Px_C1001001;
  abcd[15] = I_MOM_D2x_S_Py_C1001001+ABX*I_MOM_Px_S_Py_C1001001;
  abcd[16] = I_MOM_Dxy_S_Py_C1001001+ABX*I_MOM_Py_S_Py_C1001001;
  abcd[17] = I_MOM_Dxz_S_Py_C1001001+ABX*I_MOM_Pz_S_Py_C1001001;
  abcd[18] = I_MOM_Dxy_S_Py_C1001001+ABY*I_MOM_Px_S_Py_C1001001;
  abcd[19] = I_MOM_D2y_S_Py_C1001001+ABY*I_MOM_Py_S_Py_C1001001;
  abcd[20] = I_MOM_Dyz_S_Py_C1001001+ABY*I_MOM_Pz_S_Py_C1001001;
  abcd[21] = I_MOM_Dxz_S_Py_C1001001+ABZ*I_MOM_Px_S_Py_C1001001;
  abcd[22] = I_MOM_Dyz_S_Py_C1001001+ABZ*I_MOM_Py_S_Py_C1001001;
  abcd[23] = I_MOM_D2z_S_Py_C1001001+ABZ*I_MOM_Pz_S_Py_C1001001;
  abcd[27] = I_MOM_D2x_S_Pz_C1001001+ABX*I_MOM_Px_S_Pz_C1001001;
  abcd[28] = I_MOM_Dxy_S_Pz_C1001001+ABX*I_MOM_Py_S_Pz_C1001001;
  abcd[29] = I_MOM_Dxz_S_Pz_C1001001+ABX*I_MOM_Pz_S_Pz_C1001001;
  abcd[30] = I_MOM_Dxy_S_Pz_C1001001+ABY*I_MOM_Px_S_Pz_C1001001;
  abcd[31] = I_MOM_D2y_S_Pz_C1001001+ABY*I_MOM_Py_S_Pz_C1001001;
  abcd[32] = I_MOM_Dyz_S_Pz_C1001001+ABY*I_MOM_Pz_S_Pz_C1001001;
  abcd[33] = I_MOM_Dxz_S_Pz_C1001001+ABZ*I_MOM_Px_S_Pz_C1001001;
  abcd[34] = I_MOM_Dyz_S_Pz_C1001001+ABZ*I_MOM_Py_S_Pz_C1001001;
  abcd[35] = I_MOM_D2z_S_Pz_C1001001+ABZ*I_MOM_Pz_S_Pz_C1001001;
}
