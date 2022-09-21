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
// BRA1
// X
// Y
// Z
// ####

void hgp_os_kinetic_d_d_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_F3x_D2x_a = 0.0E0;
  Double I_KINETIC_F2xy_D2x_a = 0.0E0;
  Double I_KINETIC_F2xz_D2x_a = 0.0E0;
  Double I_KINETIC_Fx2y_D2x_a = 0.0E0;
  Double I_KINETIC_Fxyz_D2x_a = 0.0E0;
  Double I_KINETIC_Fx2z_D2x_a = 0.0E0;
  Double I_KINETIC_F3y_D2x_a = 0.0E0;
  Double I_KINETIC_F2yz_D2x_a = 0.0E0;
  Double I_KINETIC_Fy2z_D2x_a = 0.0E0;
  Double I_KINETIC_F3z_D2x_a = 0.0E0;
  Double I_KINETIC_F3x_Dxy_a = 0.0E0;
  Double I_KINETIC_F2xy_Dxy_a = 0.0E0;
  Double I_KINETIC_F2xz_Dxy_a = 0.0E0;
  Double I_KINETIC_Fx2y_Dxy_a = 0.0E0;
  Double I_KINETIC_Fxyz_Dxy_a = 0.0E0;
  Double I_KINETIC_Fx2z_Dxy_a = 0.0E0;
  Double I_KINETIC_F3y_Dxy_a = 0.0E0;
  Double I_KINETIC_F2yz_Dxy_a = 0.0E0;
  Double I_KINETIC_Fy2z_Dxy_a = 0.0E0;
  Double I_KINETIC_F3z_Dxy_a = 0.0E0;
  Double I_KINETIC_F3x_Dxz_a = 0.0E0;
  Double I_KINETIC_F2xy_Dxz_a = 0.0E0;
  Double I_KINETIC_F2xz_Dxz_a = 0.0E0;
  Double I_KINETIC_Fx2y_Dxz_a = 0.0E0;
  Double I_KINETIC_Fxyz_Dxz_a = 0.0E0;
  Double I_KINETIC_Fx2z_Dxz_a = 0.0E0;
  Double I_KINETIC_F3y_Dxz_a = 0.0E0;
  Double I_KINETIC_F2yz_Dxz_a = 0.0E0;
  Double I_KINETIC_Fy2z_Dxz_a = 0.0E0;
  Double I_KINETIC_F3z_Dxz_a = 0.0E0;
  Double I_KINETIC_F3x_D2y_a = 0.0E0;
  Double I_KINETIC_F2xy_D2y_a = 0.0E0;
  Double I_KINETIC_F2xz_D2y_a = 0.0E0;
  Double I_KINETIC_Fx2y_D2y_a = 0.0E0;
  Double I_KINETIC_Fxyz_D2y_a = 0.0E0;
  Double I_KINETIC_Fx2z_D2y_a = 0.0E0;
  Double I_KINETIC_F3y_D2y_a = 0.0E0;
  Double I_KINETIC_F2yz_D2y_a = 0.0E0;
  Double I_KINETIC_Fy2z_D2y_a = 0.0E0;
  Double I_KINETIC_F3z_D2y_a = 0.0E0;
  Double I_KINETIC_F3x_Dyz_a = 0.0E0;
  Double I_KINETIC_F2xy_Dyz_a = 0.0E0;
  Double I_KINETIC_F2xz_Dyz_a = 0.0E0;
  Double I_KINETIC_Fx2y_Dyz_a = 0.0E0;
  Double I_KINETIC_Fxyz_Dyz_a = 0.0E0;
  Double I_KINETIC_Fx2z_Dyz_a = 0.0E0;
  Double I_KINETIC_F3y_Dyz_a = 0.0E0;
  Double I_KINETIC_F2yz_Dyz_a = 0.0E0;
  Double I_KINETIC_Fy2z_Dyz_a = 0.0E0;
  Double I_KINETIC_F3z_Dyz_a = 0.0E0;
  Double I_KINETIC_F3x_D2z_a = 0.0E0;
  Double I_KINETIC_F2xy_D2z_a = 0.0E0;
  Double I_KINETIC_F2xz_D2z_a = 0.0E0;
  Double I_KINETIC_Fx2y_D2z_a = 0.0E0;
  Double I_KINETIC_Fxyz_D2z_a = 0.0E0;
  Double I_KINETIC_Fx2z_D2z_a = 0.0E0;
  Double I_KINETIC_F3y_D2z_a = 0.0E0;
  Double I_KINETIC_F2yz_D2z_a = 0.0E0;
  Double I_KINETIC_Fy2z_D2z_a = 0.0E0;
  Double I_KINETIC_F3z_D2z_a = 0.0E0;
  Double I_KINETIC_Px_D2x = 0.0E0;
  Double I_KINETIC_Py_D2x = 0.0E0;
  Double I_KINETIC_Pz_D2x = 0.0E0;
  Double I_KINETIC_Px_Dxy = 0.0E0;
  Double I_KINETIC_Py_Dxy = 0.0E0;
  Double I_KINETIC_Pz_Dxy = 0.0E0;
  Double I_KINETIC_Px_Dxz = 0.0E0;
  Double I_KINETIC_Py_Dxz = 0.0E0;
  Double I_KINETIC_Pz_Dxz = 0.0E0;
  Double I_KINETIC_Px_D2y = 0.0E0;
  Double I_KINETIC_Py_D2y = 0.0E0;
  Double I_KINETIC_Pz_D2y = 0.0E0;
  Double I_KINETIC_Px_Dyz = 0.0E0;
  Double I_KINETIC_Py_Dyz = 0.0E0;
  Double I_KINETIC_Pz_Dyz = 0.0E0;
  Double I_KINETIC_Px_D2z = 0.0E0;
  Double I_KINETIC_Py_D2z = 0.0E0;
  Double I_KINETIC_Pz_D2z = 0.0E0;

  Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double onedz = iexp[ip2];
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double xi    = alpha*beta*onedz;
    Double twoxi = 2.0E0*xi;
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    Double adz   = alpha*onedz;
    Double bdz   = beta*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    Double PBX   = PX - B[0];
    Double PBY   = PY - B[1];
    Double PBZ   = PZ - B[2];
    Double I_KINETIC_S_S_vrr = ic2*fbra*xi*(3.0E0-twoxi*AB2);
    Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;
    if(fabs(I_KINETIC_S_S_vrr)<THRESHOLD_MATH) continue;


    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_S_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_S_Px_vrr = PBX*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_S_Py_vrr = PBY*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_S_Pz_vrr = PBZ*I_TWOBODYOVERLAP_S_S_vrr;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_Px_vrr = PBX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Py_Px_vrr = PBX*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_Px_vrr = PBX*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_Px_Py_vrr = PBY*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Py_Py_vrr = PBY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_Py_vrr = PBY*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_Px_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Py_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_S_D2x_vrr = PBX*I_TWOBODYOVERLAP_S_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_S_Dxy_vrr = PBY*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_S_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_S_D2y_vrr = PBY*I_TWOBODYOVERLAP_S_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_S_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_S_D2z_vrr = PBZ*I_TWOBODYOVERLAP_S_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_D2x_vrr = PAX*I_TWOBODYOVERLAP_S_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_Py_D2x_vrr = PAY*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Pz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Px_Dxy_vrr = PAX*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Py_Dxy_vrr = PAY*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_Pz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Px_Dxz_vrr = PAX*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_Py_Dxz_vrr = PAY*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Pz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_Px_D2y_vrr = PAX*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Py_D2y_vrr = PAY*I_TWOBODYOVERLAP_S_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Pz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Px_Dyz_vrr = PAX*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Py_Dyz_vrr = PAY*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_Pz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Px_D2z_vrr = PAX*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Py_D2z_vrr = PAY*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Pz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_S_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_Px_vrr = PAX*I_TWOBODYOVERLAP_Px_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Px_vrr = PAY*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_Px_vrr = PAY*I_TWOBODYOVERLAP_Py_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_Py_vrr = PAX*I_TWOBODYOVERLAP_Px_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Py_vrr = PAY*I_TWOBODYOVERLAP_Px_Py_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Py_vrr = PAY*I_TWOBODYOVERLAP_Py_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_Pz_vrr = PAX*I_TWOBODYOVERLAP_Px_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Pz_vrr = PAY*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_Py_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_D2x_vrr = PAX*I_TWOBODYOVERLAP_Px_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2x_vrr = PAY*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2x_vrr = PAY*I_TWOBODYOVERLAP_Py_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Py_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Px_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Py_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2y_vrr = PAX*I_TWOBODYOVERLAP_Px_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2y_vrr = PAY*I_TWOBODYOVERLAP_Px_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2y_vrr = PAY*I_TWOBODYOVERLAP_Py_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Px_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Px_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Py_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2z_vrr = PAX*I_TWOBODYOVERLAP_Px_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxy_D2z_vrr = PAY*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2z_vrr = PAY*I_TWOBODYOVERLAP_Py_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_D2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_D2z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_D2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2y_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2y_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2y_vrr = PAX*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2y_vrr = PAX*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2y_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2y_vrr = PAY*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_D2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_D2z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2z_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2z_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2z_vrr = PAX*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2z_vrr = PAX*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2z_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2z_vrr = PAY*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_S_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_KINETIC_S_Px_vrr = PBX*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_S_Py_vrr = PBY*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_S_Pz_vrr = PBZ*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_Px_S_vrr = PAX*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_Py_S_vrr = PAY*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_Pz_S_vrr = PAZ*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_KINETIC_Px_Px_vrr = PBX*I_KINETIC_Px_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_KINETIC_Py_Px_vrr = PBX*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_KINETIC_Pz_Px_vrr = PBX*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_KINETIC_Px_Py_vrr = PBY*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_KINETIC_Py_Py_vrr = PBY*I_KINETIC_Py_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_KINETIC_Pz_Py_vrr = PBY*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_KINETIC_Px_Pz_vrr = PBZ*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_KINETIC_Py_Pz_vrr = PBZ*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_KINETIC_Pz_Pz_vrr = PBZ*I_KINETIC_Pz_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_S_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_KINETIC_S_D2x_vrr = PBX*I_KINETIC_S_Px_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_D2x_vrr-adz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_S_Dxy_vrr = PBY*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_S_Dxz_vrr = PBZ*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_S_D2y_vrr = PBY*I_KINETIC_S_Py_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_D2y_vrr-adz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_S_Dyz_vrr = PBZ*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_S_D2z_vrr = PBZ*I_KINETIC_S_Pz_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_S_D2z_vrr-adz*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_D
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     ************************************************************/
    Double I_KINETIC_Px_D2x_vrr = PAX*I_KINETIC_S_D2x_vrr+2*oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_KINETIC_Py_D2x_vrr = PAY*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_KINETIC_Pz_D2x_vrr = PAZ*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_KINETIC_Px_Dxy_vrr = PAX*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_KINETIC_Py_Dxy_vrr = PAY*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_KINETIC_Pz_Dxy_vrr = PAZ*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_KINETIC_Px_Dxz_vrr = PAX*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_KINETIC_Py_Dxz_vrr = PAY*I_KINETIC_S_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_KINETIC_Pz_Dxz_vrr = PAZ*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Dxz_vrr;
    Double I_KINETIC_Px_D2y_vrr = PAX*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_KINETIC_Py_D2y_vrr = PAY*I_KINETIC_S_D2y_vrr+2*oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_KINETIC_Pz_D2y_vrr = PAZ*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_KINETIC_Px_Dyz_vrr = PAX*I_KINETIC_S_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_KINETIC_Py_Dyz_vrr = PAY*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Py_Dyz_vrr;
    Double I_KINETIC_Pz_Dyz_vrr = PAZ*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Dyz_vrr;
    Double I_KINETIC_Px_D2z_vrr = PAX*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_KINETIC_Py_D2z_vrr = PAY*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_KINETIC_Pz_D2z_vrr = PAZ*I_KINETIC_S_D2z_vrr+2*oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Pz_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 6 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_KINETIC_D2x_Px_vrr = PAX*I_KINETIC_Px_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_Dxy_Px_vrr = PAY*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_KINETIC_D2y_Px_vrr = PAY*I_KINETIC_Py_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_D2z_Px_vrr = PAZ*I_KINETIC_Pz_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_D2x_Py_vrr = PAX*I_KINETIC_Px_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_Dxy_Py_vrr = PAY*I_KINETIC_Px_Py_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_KINETIC_D2y_Py_vrr = PAY*I_KINETIC_Py_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_D2z_Py_vrr = PAZ*I_KINETIC_Pz_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_D2x_Pz_vrr = PAX*I_KINETIC_Px_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_KINETIC_Dxy_Pz_vrr = PAY*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_KINETIC_D2y_Pz_vrr = PAY*I_KINETIC_Py_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_KINETIC_D2z_Pz_vrr = PAZ*I_KINETIC_Pz_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_S_D
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     ************************************************************/
    Double I_KINETIC_D2x_D2x_vrr = PAX*I_KINETIC_Px_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+2*oned2z*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_Dxy_D2x_vrr = PAY*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2x_vrr;
    Double I_KINETIC_D2y_D2x_vrr = PAY*I_KINETIC_Py_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_D2z_D2x_vrr = PAZ*I_KINETIC_Pz_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_D2x_Dxy_vrr = PAX*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_Dxy_Dxy_vrr = PAY*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Dxy_vrr;
    Double I_KINETIC_D2y_Dxy_vrr = PAY*I_KINETIC_Py_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2z_Dxy_vrr = PAZ*I_KINETIC_Pz_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2x_Dxz_vrr = PAX*I_KINETIC_Px_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_Dxy_Dxz_vrr = PAY*I_KINETIC_Px_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Dxz_vrr;
    Double I_KINETIC_D2y_Dxz_vrr = PAY*I_KINETIC_Py_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_D2z_Dxz_vrr = PAZ*I_KINETIC_Pz_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_D2x_D2y_vrr = PAX*I_KINETIC_Px_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_Dxy_D2y_vrr = PAY*I_KINETIC_Px_D2y_vrr+2*oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2y_vrr;
    Double I_KINETIC_D2y_D2y_vrr = PAY*I_KINETIC_Py_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+2*oned2z*I_KINETIC_Py_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_D2z_D2y_vrr = PAZ*I_KINETIC_Pz_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_D2x_Dyz_vrr = PAX*I_KINETIC_Px_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_Dxy_Dyz_vrr = PAY*I_KINETIC_Px_Dyz_vrr+oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Dyz_vrr;
    Double I_KINETIC_D2y_Dyz_vrr = PAY*I_KINETIC_Py_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_Py_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2z_Dyz_vrr = PAZ*I_KINETIC_Pz_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_Pz_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2x_D2z_vrr = PAX*I_KINETIC_Px_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_KINETIC_Dxy_D2z_vrr = PAY*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_D2z_vrr;
    Double I_KINETIC_D2y_D2z_vrr = PAY*I_KINETIC_Py_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_KINETIC_D2z_D2z_vrr = PAZ*I_KINETIC_Pz_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+2*oned2z*I_KINETIC_Pz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     ************************************************************/
    Double I_KINETIC_F3x_D2x_vrr = PAX*I_KINETIC_D2x_D2x_vrr+2*oned2z*I_KINETIC_Px_D2x_vrr+2*oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_KINETIC_F2xy_D2x_vrr = PAY*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_KINETIC_F2xz_D2x_vrr = PAZ*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_KINETIC_Fx2y_D2x_vrr = PAX*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_KINETIC_Fxyz_D2x_vrr = PAZ*I_KINETIC_Dxy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_KINETIC_Fx2z_D2x_vrr = PAX*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_KINETIC_F3y_D2x_vrr = PAY*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_KINETIC_F2yz_D2x_vrr = PAZ*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_KINETIC_Fy2z_D2x_vrr = PAY*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_KINETIC_F3z_D2x_vrr = PAZ*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_KINETIC_F3x_Dxy_vrr = PAX*I_KINETIC_D2x_Dxy_vrr+2*oned2z*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_KINETIC_F2xy_Dxy_vrr = PAY*I_KINETIC_D2x_Dxy_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_KINETIC_F2xz_Dxy_vrr = PAZ*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_KINETIC_Fx2y_Dxy_vrr = PAX*I_KINETIC_D2y_Dxy_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_KINETIC_Fxyz_Dxy_vrr = PAZ*I_KINETIC_Dxy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dxy_vrr;
    Double I_KINETIC_Fx2z_Dxy_vrr = PAX*I_KINETIC_D2z_Dxy_vrr+oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr;
    Double I_KINETIC_F3y_Dxy_vrr = PAY*I_KINETIC_D2y_Dxy_vrr+2*oned2z*I_KINETIC_Py_Dxy_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_KINETIC_F2yz_Dxy_vrr = PAZ*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_KINETIC_Fy2z_Dxy_vrr = PAY*I_KINETIC_D2z_Dxy_vrr+oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dxy_vrr;
    Double I_KINETIC_F3z_Dxy_vrr = PAZ*I_KINETIC_D2z_Dxy_vrr+2*oned2z*I_KINETIC_Pz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_KINETIC_F3x_Dxz_vrr = PAX*I_KINETIC_D2x_Dxz_vrr+2*oned2z*I_KINETIC_Px_Dxz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_KINETIC_F2xy_Dxz_vrr = PAY*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxz_vrr;
    Double I_KINETIC_F2xz_Dxz_vrr = PAZ*I_KINETIC_D2x_Dxz_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxz_vrr;
    Double I_KINETIC_Fx2y_Dxz_vrr = PAX*I_KINETIC_D2y_Dxz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxz_vrr;
    Double I_KINETIC_Fxyz_Dxz_vrr = PAZ*I_KINETIC_Dxy_Dxz_vrr+oned2z*I_KINETIC_Dxy_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dxz_vrr;
    Double I_KINETIC_Fx2z_Dxz_vrr = PAX*I_KINETIC_D2z_Dxz_vrr+oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxz_vrr;
    Double I_KINETIC_F3y_Dxz_vrr = PAY*I_KINETIC_D2y_Dxz_vrr+2*oned2z*I_KINETIC_Py_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_KINETIC_F2yz_Dxz_vrr = PAZ*I_KINETIC_D2y_Dxz_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxz_vrr;
    Double I_KINETIC_Fy2z_Dxz_vrr = PAY*I_KINETIC_D2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dxz_vrr;
    Double I_KINETIC_F3z_Dxz_vrr = PAZ*I_KINETIC_D2z_Dxz_vrr+2*oned2z*I_KINETIC_Pz_Dxz_vrr+oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Dxz_vrr;
    Double I_KINETIC_F3x_D2y_vrr = PAX*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_KINETIC_F2xy_D2y_vrr = PAY*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_KINETIC_F2xz_D2y_vrr = PAZ*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_KINETIC_Fx2y_D2y_vrr = PAX*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_KINETIC_Fxyz_D2y_vrr = PAZ*I_KINETIC_Dxy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_KINETIC_Fx2z_D2y_vrr = PAX*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_KINETIC_F3y_D2y_vrr = PAY*I_KINETIC_D2y_D2y_vrr+2*oned2z*I_KINETIC_Py_D2y_vrr+2*oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_KINETIC_F2yz_D2y_vrr = PAZ*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_KINETIC_Fy2z_D2y_vrr = PAY*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_KINETIC_F3z_D2y_vrr = PAZ*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_Pz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_KINETIC_F3x_Dyz_vrr = PAX*I_KINETIC_D2x_Dyz_vrr+2*oned2z*I_KINETIC_Px_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_KINETIC_F2xy_Dyz_vrr = PAY*I_KINETIC_D2x_Dyz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dyz_vrr;
    Double I_KINETIC_F2xz_Dyz_vrr = PAZ*I_KINETIC_D2x_Dyz_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dyz_vrr;
    Double I_KINETIC_Fx2y_Dyz_vrr = PAX*I_KINETIC_D2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dyz_vrr;
    Double I_KINETIC_Fxyz_Dyz_vrr = PAZ*I_KINETIC_Dxy_Dyz_vrr+oned2z*I_KINETIC_Dxy_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dyz_vrr;
    Double I_KINETIC_Fx2z_Dyz_vrr = PAX*I_KINETIC_D2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dyz_vrr;
    Double I_KINETIC_F3y_Dyz_vrr = PAY*I_KINETIC_D2y_Dyz_vrr+2*oned2z*I_KINETIC_Py_Dyz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dyz_vrr;
    Double I_KINETIC_F2yz_Dyz_vrr = PAZ*I_KINETIC_D2y_Dyz_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dyz_vrr;
    Double I_KINETIC_Fy2z_Dyz_vrr = PAY*I_KINETIC_D2z_Dyz_vrr+oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dyz_vrr;
    Double I_KINETIC_F3z_Dyz_vrr = PAZ*I_KINETIC_D2z_Dyz_vrr+2*oned2z*I_KINETIC_Pz_Dyz_vrr+oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Dyz_vrr;
    Double I_KINETIC_F3x_D2z_vrr = PAX*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_KINETIC_F2xy_D2z_vrr = PAY*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_KINETIC_F2xz_D2z_vrr = PAZ*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_KINETIC_Fx2y_D2z_vrr = PAX*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_KINETIC_Fxyz_D2z_vrr = PAZ*I_KINETIC_Dxy_D2z_vrr+2*oned2z*I_KINETIC_Dxy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2z_vrr;
    Double I_KINETIC_Fx2z_D2z_vrr = PAX*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_KINETIC_F3y_D2z_vrr = PAY*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_KINETIC_F2yz_D2z_vrr = PAZ*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_KINETIC_Fy2z_D2z_vrr = PAY*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2z_vrr;
    Double I_KINETIC_F3z_D2z_vrr = PAZ*I_KINETIC_D2z_D2z_vrr+2*oned2z*I_KINETIC_Pz_D2z_vrr+2*oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_F_D_a_coefs = alpha;
    I_KINETIC_F3x_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3x_D2x_vrr;
    I_KINETIC_F2xy_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xy_D2x_vrr;
    I_KINETIC_F2xz_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xz_D2x_vrr;
    I_KINETIC_Fx2y_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2y_D2x_vrr;
    I_KINETIC_Fxyz_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fxyz_D2x_vrr;
    I_KINETIC_Fx2z_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2z_D2x_vrr;
    I_KINETIC_F3y_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3y_D2x_vrr;
    I_KINETIC_F2yz_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2yz_D2x_vrr;
    I_KINETIC_Fy2z_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fy2z_D2x_vrr;
    I_KINETIC_F3z_D2x_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3z_D2x_vrr;
    I_KINETIC_F3x_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3x_Dxy_vrr;
    I_KINETIC_F2xy_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xy_Dxy_vrr;
    I_KINETIC_F2xz_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xz_Dxy_vrr;
    I_KINETIC_Fx2y_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2y_Dxy_vrr;
    I_KINETIC_Fxyz_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fxyz_Dxy_vrr;
    I_KINETIC_Fx2z_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2z_Dxy_vrr;
    I_KINETIC_F3y_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3y_Dxy_vrr;
    I_KINETIC_F2yz_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2yz_Dxy_vrr;
    I_KINETIC_Fy2z_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fy2z_Dxy_vrr;
    I_KINETIC_F3z_Dxy_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3z_Dxy_vrr;
    I_KINETIC_F3x_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3x_Dxz_vrr;
    I_KINETIC_F2xy_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xy_Dxz_vrr;
    I_KINETIC_F2xz_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xz_Dxz_vrr;
    I_KINETIC_Fx2y_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2y_Dxz_vrr;
    I_KINETIC_Fxyz_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fxyz_Dxz_vrr;
    I_KINETIC_Fx2z_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2z_Dxz_vrr;
    I_KINETIC_F3y_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3y_Dxz_vrr;
    I_KINETIC_F2yz_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2yz_Dxz_vrr;
    I_KINETIC_Fy2z_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fy2z_Dxz_vrr;
    I_KINETIC_F3z_Dxz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3z_Dxz_vrr;
    I_KINETIC_F3x_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3x_D2y_vrr;
    I_KINETIC_F2xy_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xy_D2y_vrr;
    I_KINETIC_F2xz_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xz_D2y_vrr;
    I_KINETIC_Fx2y_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2y_D2y_vrr;
    I_KINETIC_Fxyz_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fxyz_D2y_vrr;
    I_KINETIC_Fx2z_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2z_D2y_vrr;
    I_KINETIC_F3y_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3y_D2y_vrr;
    I_KINETIC_F2yz_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2yz_D2y_vrr;
    I_KINETIC_Fy2z_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fy2z_D2y_vrr;
    I_KINETIC_F3z_D2y_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3z_D2y_vrr;
    I_KINETIC_F3x_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3x_Dyz_vrr;
    I_KINETIC_F2xy_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xy_Dyz_vrr;
    I_KINETIC_F2xz_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xz_Dyz_vrr;
    I_KINETIC_Fx2y_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2y_Dyz_vrr;
    I_KINETIC_Fxyz_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fxyz_Dyz_vrr;
    I_KINETIC_Fx2z_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2z_Dyz_vrr;
    I_KINETIC_F3y_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3y_Dyz_vrr;
    I_KINETIC_F2yz_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2yz_Dyz_vrr;
    I_KINETIC_Fy2z_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fy2z_Dyz_vrr;
    I_KINETIC_F3z_Dyz_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3z_Dyz_vrr;
    I_KINETIC_F3x_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3x_D2z_vrr;
    I_KINETIC_F2xy_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xy_D2z_vrr;
    I_KINETIC_F2xz_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2xz_D2z_vrr;
    I_KINETIC_Fx2y_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2y_D2z_vrr;
    I_KINETIC_Fxyz_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fxyz_D2z_vrr;
    I_KINETIC_Fx2z_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fx2z_D2z_vrr;
    I_KINETIC_F3y_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3y_D2z_vrr;
    I_KINETIC_F2yz_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F2yz_D2z_vrr;
    I_KINETIC_Fy2z_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_Fy2z_D2z_vrr;
    I_KINETIC_F3z_D2z_a += SQ_KINETIC_F_D_a_coefs*I_KINETIC_F3z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_D
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_Px_D2x += I_KINETIC_Px_D2x_vrr;
    I_KINETIC_Py_D2x += I_KINETIC_Py_D2x_vrr;
    I_KINETIC_Pz_D2x += I_KINETIC_Pz_D2x_vrr;
    I_KINETIC_Px_Dxy += I_KINETIC_Px_Dxy_vrr;
    I_KINETIC_Py_Dxy += I_KINETIC_Py_Dxy_vrr;
    I_KINETIC_Pz_Dxy += I_KINETIC_Pz_Dxy_vrr;
    I_KINETIC_Px_Dxz += I_KINETIC_Px_Dxz_vrr;
    I_KINETIC_Py_Dxz += I_KINETIC_Py_Dxz_vrr;
    I_KINETIC_Pz_Dxz += I_KINETIC_Pz_Dxz_vrr;
    I_KINETIC_Px_D2y += I_KINETIC_Px_D2y_vrr;
    I_KINETIC_Py_D2y += I_KINETIC_Py_D2y_vrr;
    I_KINETIC_Pz_D2y += I_KINETIC_Pz_D2y_vrr;
    I_KINETIC_Px_Dyz += I_KINETIC_Px_Dyz_vrr;
    I_KINETIC_Py_Dyz += I_KINETIC_Py_Dyz_vrr;
    I_KINETIC_Pz_Dyz += I_KINETIC_Pz_Dyz_vrr;
    I_KINETIC_Px_D2z += I_KINETIC_Px_D2z_vrr;
    I_KINETIC_Py_D2z += I_KINETIC_Py_D2z_vrr;
    I_KINETIC_Pz_D2z += I_KINETIC_Pz_D2z_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_D_D_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_P_D
   ************************************************************/
  abcd[0] = 2.0E0*I_KINETIC_F3x_D2x_a-2*I_KINETIC_Px_D2x;
  abcd[1] = 2.0E0*I_KINETIC_F2xy_D2x_a-1*I_KINETIC_Py_D2x;
  abcd[2] = 2.0E0*I_KINETIC_F2xz_D2x_a-1*I_KINETIC_Pz_D2x;
  abcd[3] = 2.0E0*I_KINETIC_Fx2y_D2x_a;
  abcd[4] = 2.0E0*I_KINETIC_Fxyz_D2x_a;
  abcd[5] = 2.0E0*I_KINETIC_Fx2z_D2x_a;
  abcd[6] = 2.0E0*I_KINETIC_F3x_Dxy_a-2*I_KINETIC_Px_Dxy;
  abcd[7] = 2.0E0*I_KINETIC_F2xy_Dxy_a-1*I_KINETIC_Py_Dxy;
  abcd[8] = 2.0E0*I_KINETIC_F2xz_Dxy_a-1*I_KINETIC_Pz_Dxy;
  abcd[9] = 2.0E0*I_KINETIC_Fx2y_Dxy_a;
  abcd[10] = 2.0E0*I_KINETIC_Fxyz_Dxy_a;
  abcd[11] = 2.0E0*I_KINETIC_Fx2z_Dxy_a;
  abcd[12] = 2.0E0*I_KINETIC_F3x_Dxz_a-2*I_KINETIC_Px_Dxz;
  abcd[13] = 2.0E0*I_KINETIC_F2xy_Dxz_a-1*I_KINETIC_Py_Dxz;
  abcd[14] = 2.0E0*I_KINETIC_F2xz_Dxz_a-1*I_KINETIC_Pz_Dxz;
  abcd[15] = 2.0E0*I_KINETIC_Fx2y_Dxz_a;
  abcd[16] = 2.0E0*I_KINETIC_Fxyz_Dxz_a;
  abcd[17] = 2.0E0*I_KINETIC_Fx2z_Dxz_a;
  abcd[18] = 2.0E0*I_KINETIC_F3x_D2y_a-2*I_KINETIC_Px_D2y;
  abcd[19] = 2.0E0*I_KINETIC_F2xy_D2y_a-1*I_KINETIC_Py_D2y;
  abcd[20] = 2.0E0*I_KINETIC_F2xz_D2y_a-1*I_KINETIC_Pz_D2y;
  abcd[21] = 2.0E0*I_KINETIC_Fx2y_D2y_a;
  abcd[22] = 2.0E0*I_KINETIC_Fxyz_D2y_a;
  abcd[23] = 2.0E0*I_KINETIC_Fx2z_D2y_a;
  abcd[24] = 2.0E0*I_KINETIC_F3x_Dyz_a-2*I_KINETIC_Px_Dyz;
  abcd[25] = 2.0E0*I_KINETIC_F2xy_Dyz_a-1*I_KINETIC_Py_Dyz;
  abcd[26] = 2.0E0*I_KINETIC_F2xz_Dyz_a-1*I_KINETIC_Pz_Dyz;
  abcd[27] = 2.0E0*I_KINETIC_Fx2y_Dyz_a;
  abcd[28] = 2.0E0*I_KINETIC_Fxyz_Dyz_a;
  abcd[29] = 2.0E0*I_KINETIC_Fx2z_Dyz_a;
  abcd[30] = 2.0E0*I_KINETIC_F3x_D2z_a-2*I_KINETIC_Px_D2z;
  abcd[31] = 2.0E0*I_KINETIC_F2xy_D2z_a-1*I_KINETIC_Py_D2z;
  abcd[32] = 2.0E0*I_KINETIC_F2xz_D2z_a-1*I_KINETIC_Pz_D2z;
  abcd[33] = 2.0E0*I_KINETIC_Fx2y_D2z_a;
  abcd[34] = 2.0E0*I_KINETIC_Fxyz_D2z_a;
  abcd[35] = 2.0E0*I_KINETIC_Fx2z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_D_D_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_P_D
   ************************************************************/
  abcd[36] = 2.0E0*I_KINETIC_F2xy_D2x_a;
  abcd[37] = 2.0E0*I_KINETIC_Fx2y_D2x_a-1*I_KINETIC_Px_D2x;
  abcd[38] = 2.0E0*I_KINETIC_Fxyz_D2x_a;
  abcd[39] = 2.0E0*I_KINETIC_F3y_D2x_a-2*I_KINETIC_Py_D2x;
  abcd[40] = 2.0E0*I_KINETIC_F2yz_D2x_a-1*I_KINETIC_Pz_D2x;
  abcd[41] = 2.0E0*I_KINETIC_Fy2z_D2x_a;
  abcd[42] = 2.0E0*I_KINETIC_F2xy_Dxy_a;
  abcd[43] = 2.0E0*I_KINETIC_Fx2y_Dxy_a-1*I_KINETIC_Px_Dxy;
  abcd[44] = 2.0E0*I_KINETIC_Fxyz_Dxy_a;
  abcd[45] = 2.0E0*I_KINETIC_F3y_Dxy_a-2*I_KINETIC_Py_Dxy;
  abcd[46] = 2.0E0*I_KINETIC_F2yz_Dxy_a-1*I_KINETIC_Pz_Dxy;
  abcd[47] = 2.0E0*I_KINETIC_Fy2z_Dxy_a;
  abcd[48] = 2.0E0*I_KINETIC_F2xy_Dxz_a;
  abcd[49] = 2.0E0*I_KINETIC_Fx2y_Dxz_a-1*I_KINETIC_Px_Dxz;
  abcd[50] = 2.0E0*I_KINETIC_Fxyz_Dxz_a;
  abcd[51] = 2.0E0*I_KINETIC_F3y_Dxz_a-2*I_KINETIC_Py_Dxz;
  abcd[52] = 2.0E0*I_KINETIC_F2yz_Dxz_a-1*I_KINETIC_Pz_Dxz;
  abcd[53] = 2.0E0*I_KINETIC_Fy2z_Dxz_a;
  abcd[54] = 2.0E0*I_KINETIC_F2xy_D2y_a;
  abcd[55] = 2.0E0*I_KINETIC_Fx2y_D2y_a-1*I_KINETIC_Px_D2y;
  abcd[56] = 2.0E0*I_KINETIC_Fxyz_D2y_a;
  abcd[57] = 2.0E0*I_KINETIC_F3y_D2y_a-2*I_KINETIC_Py_D2y;
  abcd[58] = 2.0E0*I_KINETIC_F2yz_D2y_a-1*I_KINETIC_Pz_D2y;
  abcd[59] = 2.0E0*I_KINETIC_Fy2z_D2y_a;
  abcd[60] = 2.0E0*I_KINETIC_F2xy_Dyz_a;
  abcd[61] = 2.0E0*I_KINETIC_Fx2y_Dyz_a-1*I_KINETIC_Px_Dyz;
  abcd[62] = 2.0E0*I_KINETIC_Fxyz_Dyz_a;
  abcd[63] = 2.0E0*I_KINETIC_F3y_Dyz_a-2*I_KINETIC_Py_Dyz;
  abcd[64] = 2.0E0*I_KINETIC_F2yz_Dyz_a-1*I_KINETIC_Pz_Dyz;
  abcd[65] = 2.0E0*I_KINETIC_Fy2z_Dyz_a;
  abcd[66] = 2.0E0*I_KINETIC_F2xy_D2z_a;
  abcd[67] = 2.0E0*I_KINETIC_Fx2y_D2z_a-1*I_KINETIC_Px_D2z;
  abcd[68] = 2.0E0*I_KINETIC_Fxyz_D2z_a;
  abcd[69] = 2.0E0*I_KINETIC_F3y_D2z_a-2*I_KINETIC_Py_D2z;
  abcd[70] = 2.0E0*I_KINETIC_F2yz_D2z_a-1*I_KINETIC_Pz_D2z;
  abcd[71] = 2.0E0*I_KINETIC_Fy2z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_D_D_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_P_D
   ************************************************************/
  abcd[72] = 2.0E0*I_KINETIC_F2xz_D2x_a;
  abcd[73] = 2.0E0*I_KINETIC_Fxyz_D2x_a;
  abcd[74] = 2.0E0*I_KINETIC_Fx2z_D2x_a-1*I_KINETIC_Px_D2x;
  abcd[75] = 2.0E0*I_KINETIC_F2yz_D2x_a;
  abcd[76] = 2.0E0*I_KINETIC_Fy2z_D2x_a-1*I_KINETIC_Py_D2x;
  abcd[77] = 2.0E0*I_KINETIC_F3z_D2x_a-2*I_KINETIC_Pz_D2x;
  abcd[78] = 2.0E0*I_KINETIC_F2xz_Dxy_a;
  abcd[79] = 2.0E0*I_KINETIC_Fxyz_Dxy_a;
  abcd[80] = 2.0E0*I_KINETIC_Fx2z_Dxy_a-1*I_KINETIC_Px_Dxy;
  abcd[81] = 2.0E0*I_KINETIC_F2yz_Dxy_a;
  abcd[82] = 2.0E0*I_KINETIC_Fy2z_Dxy_a-1*I_KINETIC_Py_Dxy;
  abcd[83] = 2.0E0*I_KINETIC_F3z_Dxy_a-2*I_KINETIC_Pz_Dxy;
  abcd[84] = 2.0E0*I_KINETIC_F2xz_Dxz_a;
  abcd[85] = 2.0E0*I_KINETIC_Fxyz_Dxz_a;
  abcd[86] = 2.0E0*I_KINETIC_Fx2z_Dxz_a-1*I_KINETIC_Px_Dxz;
  abcd[87] = 2.0E0*I_KINETIC_F2yz_Dxz_a;
  abcd[88] = 2.0E0*I_KINETIC_Fy2z_Dxz_a-1*I_KINETIC_Py_Dxz;
  abcd[89] = 2.0E0*I_KINETIC_F3z_Dxz_a-2*I_KINETIC_Pz_Dxz;
  abcd[90] = 2.0E0*I_KINETIC_F2xz_D2y_a;
  abcd[91] = 2.0E0*I_KINETIC_Fxyz_D2y_a;
  abcd[92] = 2.0E0*I_KINETIC_Fx2z_D2y_a-1*I_KINETIC_Px_D2y;
  abcd[93] = 2.0E0*I_KINETIC_F2yz_D2y_a;
  abcd[94] = 2.0E0*I_KINETIC_Fy2z_D2y_a-1*I_KINETIC_Py_D2y;
  abcd[95] = 2.0E0*I_KINETIC_F3z_D2y_a-2*I_KINETIC_Pz_D2y;
  abcd[96] = 2.0E0*I_KINETIC_F2xz_Dyz_a;
  abcd[97] = 2.0E0*I_KINETIC_Fxyz_Dyz_a;
  abcd[98] = 2.0E0*I_KINETIC_Fx2z_Dyz_a-1*I_KINETIC_Px_Dyz;
  abcd[99] = 2.0E0*I_KINETIC_F2yz_Dyz_a;
  abcd[100] = 2.0E0*I_KINETIC_Fy2z_Dyz_a-1*I_KINETIC_Py_Dyz;
  abcd[101] = 2.0E0*I_KINETIC_F3z_Dyz_a-2*I_KINETIC_Pz_Dyz;
  abcd[102] = 2.0E0*I_KINETIC_F2xz_D2z_a;
  abcd[103] = 2.0E0*I_KINETIC_Fxyz_D2z_a;
  abcd[104] = 2.0E0*I_KINETIC_Fx2z_D2z_a-1*I_KINETIC_Px_D2z;
  abcd[105] = 2.0E0*I_KINETIC_F2yz_D2z_a;
  abcd[106] = 2.0E0*I_KINETIC_Fy2z_D2z_a-1*I_KINETIC_Py_D2z;
  abcd[107] = 2.0E0*I_KINETIC_F3z_D2z_a-2*I_KINETIC_Pz_D2z;
}
