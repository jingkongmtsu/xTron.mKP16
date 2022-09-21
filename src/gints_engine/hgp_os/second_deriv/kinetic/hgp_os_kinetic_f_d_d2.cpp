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
// BRA1 as redundant position, total RHS integrals evaluated as: 3206
// BRA2 as redundant position, total RHS integrals evaluated as: 3009
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

void hgp_os_kinetic_f_d_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_H5x_D2x_aa = 0.0E0;
  Double I_KINETIC_H4xy_D2x_aa = 0.0E0;
  Double I_KINETIC_H4xz_D2x_aa = 0.0E0;
  Double I_KINETIC_H3x2y_D2x_aa = 0.0E0;
  Double I_KINETIC_H3xyz_D2x_aa = 0.0E0;
  Double I_KINETIC_H3x2z_D2x_aa = 0.0E0;
  Double I_KINETIC_H2x3y_D2x_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_D2x_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_D2x_aa = 0.0E0;
  Double I_KINETIC_H2x3z_D2x_aa = 0.0E0;
  Double I_KINETIC_Hx4y_D2x_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_D2x_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_D2x_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_D2x_aa = 0.0E0;
  Double I_KINETIC_Hx4z_D2x_aa = 0.0E0;
  Double I_KINETIC_H5y_D2x_aa = 0.0E0;
  Double I_KINETIC_H4yz_D2x_aa = 0.0E0;
  Double I_KINETIC_H3y2z_D2x_aa = 0.0E0;
  Double I_KINETIC_H2y3z_D2x_aa = 0.0E0;
  Double I_KINETIC_Hy4z_D2x_aa = 0.0E0;
  Double I_KINETIC_H5z_D2x_aa = 0.0E0;
  Double I_KINETIC_H5x_Dxy_aa = 0.0E0;
  Double I_KINETIC_H4xy_Dxy_aa = 0.0E0;
  Double I_KINETIC_H4xz_Dxy_aa = 0.0E0;
  Double I_KINETIC_H3x2y_Dxy_aa = 0.0E0;
  Double I_KINETIC_H3xyz_Dxy_aa = 0.0E0;
  Double I_KINETIC_H3x2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_H2x3y_Dxy_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_Dxy_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_H2x3z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Hx4y_Dxy_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_Dxy_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Hx4z_Dxy_aa = 0.0E0;
  Double I_KINETIC_H5y_Dxy_aa = 0.0E0;
  Double I_KINETIC_H4yz_Dxy_aa = 0.0E0;
  Double I_KINETIC_H3y2z_Dxy_aa = 0.0E0;
  Double I_KINETIC_H2y3z_Dxy_aa = 0.0E0;
  Double I_KINETIC_Hy4z_Dxy_aa = 0.0E0;
  Double I_KINETIC_H5z_Dxy_aa = 0.0E0;
  Double I_KINETIC_H5x_Dxz_aa = 0.0E0;
  Double I_KINETIC_H4xy_Dxz_aa = 0.0E0;
  Double I_KINETIC_H4xz_Dxz_aa = 0.0E0;
  Double I_KINETIC_H3x2y_Dxz_aa = 0.0E0;
  Double I_KINETIC_H3xyz_Dxz_aa = 0.0E0;
  Double I_KINETIC_H3x2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_H2x3y_Dxz_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_Dxz_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_H2x3z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Hx4y_Dxz_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_Dxz_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Hx4z_Dxz_aa = 0.0E0;
  Double I_KINETIC_H5y_Dxz_aa = 0.0E0;
  Double I_KINETIC_H4yz_Dxz_aa = 0.0E0;
  Double I_KINETIC_H3y2z_Dxz_aa = 0.0E0;
  Double I_KINETIC_H2y3z_Dxz_aa = 0.0E0;
  Double I_KINETIC_Hy4z_Dxz_aa = 0.0E0;
  Double I_KINETIC_H5z_Dxz_aa = 0.0E0;
  Double I_KINETIC_H5x_D2y_aa = 0.0E0;
  Double I_KINETIC_H4xy_D2y_aa = 0.0E0;
  Double I_KINETIC_H4xz_D2y_aa = 0.0E0;
  Double I_KINETIC_H3x2y_D2y_aa = 0.0E0;
  Double I_KINETIC_H3xyz_D2y_aa = 0.0E0;
  Double I_KINETIC_H3x2z_D2y_aa = 0.0E0;
  Double I_KINETIC_H2x3y_D2y_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_D2y_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_D2y_aa = 0.0E0;
  Double I_KINETIC_H2x3z_D2y_aa = 0.0E0;
  Double I_KINETIC_Hx4y_D2y_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_D2y_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_D2y_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_D2y_aa = 0.0E0;
  Double I_KINETIC_Hx4z_D2y_aa = 0.0E0;
  Double I_KINETIC_H5y_D2y_aa = 0.0E0;
  Double I_KINETIC_H4yz_D2y_aa = 0.0E0;
  Double I_KINETIC_H3y2z_D2y_aa = 0.0E0;
  Double I_KINETIC_H2y3z_D2y_aa = 0.0E0;
  Double I_KINETIC_Hy4z_D2y_aa = 0.0E0;
  Double I_KINETIC_H5z_D2y_aa = 0.0E0;
  Double I_KINETIC_H5x_Dyz_aa = 0.0E0;
  Double I_KINETIC_H4xy_Dyz_aa = 0.0E0;
  Double I_KINETIC_H4xz_Dyz_aa = 0.0E0;
  Double I_KINETIC_H3x2y_Dyz_aa = 0.0E0;
  Double I_KINETIC_H3xyz_Dyz_aa = 0.0E0;
  Double I_KINETIC_H3x2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_H2x3y_Dyz_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_Dyz_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_H2x3z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Hx4y_Dyz_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_Dyz_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Hx4z_Dyz_aa = 0.0E0;
  Double I_KINETIC_H5y_Dyz_aa = 0.0E0;
  Double I_KINETIC_H4yz_Dyz_aa = 0.0E0;
  Double I_KINETIC_H3y2z_Dyz_aa = 0.0E0;
  Double I_KINETIC_H2y3z_Dyz_aa = 0.0E0;
  Double I_KINETIC_Hy4z_Dyz_aa = 0.0E0;
  Double I_KINETIC_H5z_Dyz_aa = 0.0E0;
  Double I_KINETIC_H5x_D2z_aa = 0.0E0;
  Double I_KINETIC_H4xy_D2z_aa = 0.0E0;
  Double I_KINETIC_H4xz_D2z_aa = 0.0E0;
  Double I_KINETIC_H3x2y_D2z_aa = 0.0E0;
  Double I_KINETIC_H3xyz_D2z_aa = 0.0E0;
  Double I_KINETIC_H3x2z_D2z_aa = 0.0E0;
  Double I_KINETIC_H2x3y_D2z_aa = 0.0E0;
  Double I_KINETIC_H2x2yz_D2z_aa = 0.0E0;
  Double I_KINETIC_H2xy2z_D2z_aa = 0.0E0;
  Double I_KINETIC_H2x3z_D2z_aa = 0.0E0;
  Double I_KINETIC_Hx4y_D2z_aa = 0.0E0;
  Double I_KINETIC_Hx3yz_D2z_aa = 0.0E0;
  Double I_KINETIC_Hx2y2z_D2z_aa = 0.0E0;
  Double I_KINETIC_Hxy3z_D2z_aa = 0.0E0;
  Double I_KINETIC_Hx4z_D2z_aa = 0.0E0;
  Double I_KINETIC_H5y_D2z_aa = 0.0E0;
  Double I_KINETIC_H4yz_D2z_aa = 0.0E0;
  Double I_KINETIC_H3y2z_D2z_aa = 0.0E0;
  Double I_KINETIC_H2y3z_D2z_aa = 0.0E0;
  Double I_KINETIC_Hy4z_D2z_aa = 0.0E0;
  Double I_KINETIC_H5z_D2z_aa = 0.0E0;
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
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_Px_vrr = PBX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Px_vrr = PBX*I_TWOBODYOVERLAP_Dxy_S_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Px_vrr = PBX*I_TWOBODYOVERLAP_Dxz_S_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Px_vrr = PBX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Px_vrr = PBX*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Px_vrr = PBX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_D2x_Py_vrr = PBY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Py_vrr = PBY*I_TWOBODYOVERLAP_Dxy_S_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Py_vrr = PBY*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Py_vrr = PBY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Py_vrr = PBY*I_TWOBODYOVERLAP_Dyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Py_vrr = PBY*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_D2x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Dxz_S_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Dyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_D2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_Px_vrr = PBX*I_TWOBODYOVERLAP_F3x_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Px_vrr = PBX*I_TWOBODYOVERLAP_F2xy_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Px_vrr = PBX*I_TWOBODYOVERLAP_F2xz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Px_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Px_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_S_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Px_vrr = PBX*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Px_vrr = PBX*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Py_vrr = PBY*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Py_vrr = PBY*I_TWOBODYOVERLAP_F2xy_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Py_vrr = PBY*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Py_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Py_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Py_vrr = PBY*I_TWOBODYOVERLAP_F3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_F2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Py_vrr = PBY*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_S_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_S_vrr = PAX*I_TWOBODYOVERLAP_F3x_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_S_vrr = PAY*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_S_vrr = PAY*I_TWOBODYOVERLAP_F2xy_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_S_vrr = PAX*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_S_vrr = PAX*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_S_vrr = PAY*I_TWOBODYOVERLAP_F3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_S_vrr = PAY*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_S_vrr = PAZ*I_TWOBODYOVERLAP_F3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_D2x_vrr = PBX*I_TWOBODYOVERLAP_F3x_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2x_vrr = PBX*I_TWOBODYOVERLAP_F2xy_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2x_vrr = PBX*I_TWOBODYOVERLAP_F2xz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fx2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fx2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2x_vrr = PBX*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Fy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2y_vrr = PBY*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2y_vrr = PBY*I_TWOBODYOVERLAP_F2xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2y_vrr = PBY*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fx2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxy_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fxyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Dxz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fx2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2y_vrr = PBY*I_TWOBODYOVERLAP_F3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_F2yz_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Fy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Py_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fxyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Dxy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fx2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dxz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Fy2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Dyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_F3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 9 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_Px_vrr = PAX*I_TWOBODYOVERLAP_F3x_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Px_vrr = PAY*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Px_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Px_vrr = PAX*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Px_vrr = PAX*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Px_vrr = PAY*I_TWOBODYOVERLAP_F3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Px_vrr = PAY*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_Py_vrr = PAX*I_TWOBODYOVERLAP_F3x_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Py_vrr = PAY*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Py_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Py_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Py_vrr = PAX*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Py_vrr = PAX*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4y_Py_vrr = PAY*I_TWOBODYOVERLAP_F3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Py_vrr = PAY*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4x_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3x_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Pz_vrr = PAY*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4y_Pz_vrr = PAY*I_TWOBODYOVERLAP_F3y_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Pz_vrr = PAY*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 18 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_D2x_vrr = PBX*I_TWOBODYOVERLAP_G4x_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2x_vrr = PBX*I_TWOBODYOVERLAP_G3xy_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2x_vrr = PBX*I_TWOBODYOVERLAP_G3xz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2x_vrr = PBX*I_TWOBODYOVERLAP_G4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_G3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2x_vrr = PBX*I_TWOBODYOVERLAP_G4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G3xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G4y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G3yz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dxz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_D2y_vrr = PBY*I_TWOBODYOVERLAP_G4x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2y_vrr = PBY*I_TWOBODYOVERLAP_G3xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2y_vrr = PBY*I_TWOBODYOVERLAP_G3xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2y_vrr = PBY*I_TWOBODYOVERLAP_G4y_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_G3yz_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2y_vrr = PBY*I_TWOBODYOVERLAP_G4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dyz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G4x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_D2x_vrr = PAX*I_TWOBODYOVERLAP_G4x_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2x_vrr = PAY*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2x_vrr = PAY*I_TWOBODYOVERLAP_G3xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2x_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_D2x_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2x_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2x_vrr = PAX*I_TWOBODYOVERLAP_G4y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_D2x_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2x_vrr = PAX*I_TWOBODYOVERLAP_G4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2x_vrr = PAY*I_TWOBODYOVERLAP_G4y_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2x_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2x_vrr = PAY*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_G4z_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_G4x_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_G4x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_G4y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_G4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_G4y_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dxy_vrr = PAY*I_TWOBODYOVERLAP_G4z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_G4x_Dxz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_G4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_G4z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_G4y_Dxz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dxz_vrr = PAY*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Dxz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_D2y_vrr = PAX*I_TWOBODYOVERLAP_G4x_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2y_vrr = PAY*I_TWOBODYOVERLAP_G4x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2y_vrr = PAY*I_TWOBODYOVERLAP_G3xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2y_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_D2y_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2y_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2y_vrr = PAX*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_D2y_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2y_vrr = PAX*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2y_vrr = PAY*I_TWOBODYOVERLAP_G4y_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2y_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2y_vrr = PAY*I_TWOBODYOVERLAP_G4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_G4z_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_G4x_Dyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_G4x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_H5y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_G4y_Dyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Dyz_vrr = PAY*I_TWOBODYOVERLAP_G4z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Dyz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5x_D2z_vrr = PAX*I_TWOBODYOVERLAP_G4x_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_D2z_vrr = PAY*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G4x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_D2z_vrr = PAY*I_TWOBODYOVERLAP_G3xy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_D2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_D2z_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_D2z_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_D2z_vrr = PAX*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_D2z_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_D2z_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_D2z_vrr = PAX*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_D2z_vrr = PAY*I_TWOBODYOVERLAP_G4y_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_D2z_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_D2z_vrr = PAY*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_G4z_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;

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
     * shell quartet name: SQ_KINETIC_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_KINETIC_D2x_S_vrr = PAX*I_KINETIC_Px_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_Dxy_S_vrr = PAY*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_KINETIC_Dxz_S_vrr = PAZ*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_S_vrr;
    Double I_KINETIC_D2y_S_vrr = PAY*I_KINETIC_Py_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_Dyz_S_vrr = PAZ*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_S_vrr;
    Double I_KINETIC_D2z_S_vrr = PAZ*I_KINETIC_Pz_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;

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
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_KINETIC_D2x_Px_vrr = PBX*I_KINETIC_D2x_S_vrr+2*oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_Dxy_Px_vrr = PBX*I_KINETIC_Dxy_S_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Px_vrr;
    Double I_KINETIC_Dxz_Px_vrr = PBX*I_KINETIC_Dxz_S_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Px_vrr;
    Double I_KINETIC_D2y_Px_vrr = PBX*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_Dyz_Px_vrr = PBX*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Px_vrr;
    Double I_KINETIC_D2z_Px_vrr = PBX*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_KINETIC_D2x_Py_vrr = PBY*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_Dxy_Py_vrr = PBY*I_KINETIC_Dxy_S_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Py_vrr;
    Double I_KINETIC_Dxz_Py_vrr = PBY*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Py_vrr;
    Double I_KINETIC_D2y_Py_vrr = PBY*I_KINETIC_D2y_S_vrr+2*oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_Dyz_Py_vrr = PBY*I_KINETIC_Dyz_S_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Py_vrr;
    Double I_KINETIC_D2z_Py_vrr = PBY*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_KINETIC_D2x_Pz_vrr = PBZ*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_Dxy_Pz_vrr = PBZ*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Pz_vrr;
    Double I_KINETIC_Dxz_Pz_vrr = PBZ*I_KINETIC_Dxz_S_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Pz_vrr;
    Double I_KINETIC_D2y_Pz_vrr = PBZ*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_Dyz_Pz_vrr = PBZ*I_KINETIC_Dyz_S_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Pz_vrr;
    Double I_KINETIC_D2z_Pz_vrr = PBZ*I_KINETIC_D2z_S_vrr+2*oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_F3x_S_vrr = PAX*I_KINETIC_D2x_S_vrr+2*oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_S_vrr-2*bdz*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_F2xy_S_vrr = PAY*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_S_vrr = PAZ*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_S_vrr = PAX*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_S_vrr = PAZ*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_S_vrr = PAX*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_S_vrr = PAY*I_KINETIC_D2y_S_vrr+2*oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_F2yz_S_vrr = PAZ*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_S_vrr = PAY*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_S_vrr = PAZ*I_KINETIC_D2z_S_vrr+2*oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     ************************************************************/
    Double I_KINETIC_F3x_Px_vrr = PBX*I_KINETIC_F3x_S_vrr+3*oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_KINETIC_F2xy_Px_vrr = PBX*I_KINETIC_F2xy_S_vrr+2*oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_KINETIC_F2xz_Px_vrr = PBX*I_KINETIC_F2xz_S_vrr+2*oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_KINETIC_Fx2y_Px_vrr = PBX*I_KINETIC_Fx2y_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_KINETIC_Fxyz_Px_vrr = PBX*I_KINETIC_Fxyz_S_vrr+oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Px_vrr;
    Double I_KINETIC_Fx2z_Px_vrr = PBX*I_KINETIC_Fx2z_S_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_KINETIC_F3y_Px_vrr = PBX*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_KINETIC_F2yz_Px_vrr = PBX*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_KINETIC_Fy2z_Px_vrr = PBX*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_KINETIC_F3z_Px_vrr = PBX*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_KINETIC_F3x_Py_vrr = PBY*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_KINETIC_F2xy_Py_vrr = PBY*I_KINETIC_F2xy_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_KINETIC_F2xz_Py_vrr = PBY*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Py_vrr;
    Double I_KINETIC_Fx2y_Py_vrr = PBY*I_KINETIC_Fx2y_S_vrr+2*oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_KINETIC_Fxyz_Py_vrr = PBY*I_KINETIC_Fxyz_S_vrr+oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Py_vrr;
    Double I_KINETIC_Fx2z_Py_vrr = PBY*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Py_vrr;
    Double I_KINETIC_F3y_Py_vrr = PBY*I_KINETIC_F3y_S_vrr+3*oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_KINETIC_F2yz_Py_vrr = PBY*I_KINETIC_F2yz_S_vrr+2*oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Py_vrr;
    Double I_KINETIC_Fy2z_Py_vrr = PBY*I_KINETIC_Fy2z_S_vrr+oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Py_vrr;
    Double I_KINETIC_F3z_Py_vrr = PBY*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_KINETIC_F3x_Pz_vrr = PBZ*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_KINETIC_F2xy_Pz_vrr = PBZ*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Pz_vrr;
    Double I_KINETIC_F2xz_Pz_vrr = PBZ*I_KINETIC_F2xz_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Pz_vrr;
    Double I_KINETIC_Fx2y_Pz_vrr = PBZ*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Pz_vrr;
    Double I_KINETIC_Fxyz_Pz_vrr = PBZ*I_KINETIC_Fxyz_S_vrr+oned2z*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Pz_vrr;
    Double I_KINETIC_Fx2z_Pz_vrr = PBZ*I_KINETIC_Fx2z_S_vrr+2*oned2z*I_KINETIC_Dxz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Pz_vrr;
    Double I_KINETIC_F3y_Pz_vrr = PBZ*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_KINETIC_F2yz_Pz_vrr = PBZ*I_KINETIC_F2yz_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Pz_vrr;
    Double I_KINETIC_Fy2z_Pz_vrr = PBZ*I_KINETIC_Fy2z_S_vrr+2*oned2z*I_KINETIC_Dyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Pz_vrr;
    Double I_KINETIC_F3z_Pz_vrr = PBZ*I_KINETIC_F3z_S_vrr+3*oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_KINETIC_G4x_S_vrr = PAX*I_KINETIC_F3x_S_vrr+3*oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G3xy_S_vrr = PAY*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_S_vrr = PAZ*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_S_vrr = PAY*I_KINETIC_F2xy_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G2x2z_S_vrr = PAZ*I_KINETIC_F2xz_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_Gx3y_S_vrr = PAX*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx3z_S_vrr = PAX*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_S_vrr = PAY*I_KINETIC_F3y_S_vrr+3*oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_G3yz_S_vrr = PAZ*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_S_vrr = PAZ*I_KINETIC_F2yz_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_Gy3z_S_vrr = PAY*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_S_vrr = PAZ*I_KINETIC_F3z_S_vrr+3*oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_KINETIC_F3x_D2x_vrr = PBX*I_KINETIC_F3x_Px_vrr+3*oned2z*I_KINETIC_D2x_Px_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2x_vrr-adz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_F2xy_D2x_vrr = PBX*I_KINETIC_F2xy_Px_vrr+2*oned2z*I_KINETIC_Dxy_Px_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2x_vrr-adz*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_D2x_vrr = PBX*I_KINETIC_F2xz_Px_vrr+2*oned2z*I_KINETIC_Dxz_Px_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2x_vrr-adz*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_D2x_vrr = PBX*I_KINETIC_Fx2y_Px_vrr+oned2z*I_KINETIC_D2y_Px_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2x_vrr-adz*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_D2x_vrr = PBX*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_Dyz_Px_vrr+oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2x_vrr-adz*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_D2x_vrr = PBX*I_KINETIC_Fx2z_Px_vrr+oned2z*I_KINETIC_D2z_Px_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_D2x_vrr = PBX*I_KINETIC_F3y_Px_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2x_vrr-adz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_F2yz_D2x_vrr = PBX*I_KINETIC_F2yz_Px_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2x_vrr-adz*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_D2x_vrr = PBX*I_KINETIC_Fy2z_Px_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_D2x_vrr = PBX*I_KINETIC_F3z_Px_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2x_vrr-adz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_F3x_Dxy_vrr = PBY*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_KINETIC_F2xy_Dxy_vrr = PBY*I_KINETIC_F2xy_Px_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_KINETIC_F2xz_Dxy_vrr = PBY*I_KINETIC_F2xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_KINETIC_Fx2y_Dxy_vrr = PBY*I_KINETIC_Fx2y_Px_vrr+2*oned2z*I_KINETIC_Dxy_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_KINETIC_Fxyz_Dxy_vrr = PBY*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_Dxz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dxy_vrr;
    Double I_KINETIC_Fx2z_Dxy_vrr = PBY*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr;
    Double I_KINETIC_F3y_Dxy_vrr = PBY*I_KINETIC_F3y_Px_vrr+3*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_KINETIC_F2yz_Dxy_vrr = PBY*I_KINETIC_F2yz_Px_vrr+2*oned2z*I_KINETIC_Dyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_KINETIC_Fy2z_Dxy_vrr = PBY*I_KINETIC_Fy2z_Px_vrr+oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dxy_vrr;
    Double I_KINETIC_F3z_Dxy_vrr = PBY*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_KINETIC_F3x_Dxz_vrr = PBZ*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_KINETIC_F2xy_Dxz_vrr = PBZ*I_KINETIC_F2xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxz_vrr;
    Double I_KINETIC_F2xz_Dxz_vrr = PBZ*I_KINETIC_F2xz_Px_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxz_vrr;
    Double I_KINETIC_Fx2y_Dxz_vrr = PBZ*I_KINETIC_Fx2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxz_vrr;
    Double I_KINETIC_Fxyz_Dxz_vrr = PBZ*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_Dxy_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dxz_vrr;
    Double I_KINETIC_Fx2z_Dxz_vrr = PBZ*I_KINETIC_Fx2z_Px_vrr+2*oned2z*I_KINETIC_Dxz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxz_vrr;
    Double I_KINETIC_F3y_Dxz_vrr = PBZ*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_KINETIC_F2yz_Dxz_vrr = PBZ*I_KINETIC_F2yz_Px_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxz_vrr;
    Double I_KINETIC_Fy2z_Dxz_vrr = PBZ*I_KINETIC_Fy2z_Px_vrr+2*oned2z*I_KINETIC_Dyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dxz_vrr;
    Double I_KINETIC_F3z_Dxz_vrr = PBZ*I_KINETIC_F3z_Px_vrr+3*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_KINETIC_F3x_D2y_vrr = PBY*I_KINETIC_F3x_Py_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2y_vrr-adz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_F2xy_D2y_vrr = PBY*I_KINETIC_F2xy_Py_vrr+oned2z*I_KINETIC_D2x_Py_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2y_vrr-adz*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_D2y_vrr = PBY*I_KINETIC_F2xz_Py_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2y_vrr-adz*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_D2y_vrr = PBY*I_KINETIC_Fx2y_Py_vrr+2*oned2z*I_KINETIC_Dxy_Py_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2y_vrr-adz*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_D2y_vrr = PBY*I_KINETIC_Fxyz_Py_vrr+oned2z*I_KINETIC_Dxz_Py_vrr+oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2y_vrr-adz*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_D2y_vrr = PBY*I_KINETIC_Fx2z_Py_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_D2y_vrr = PBY*I_KINETIC_F3y_Py_vrr+3*oned2z*I_KINETIC_D2y_Py_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2y_vrr-adz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_F2yz_D2y_vrr = PBY*I_KINETIC_F2yz_Py_vrr+2*oned2z*I_KINETIC_Dyz_Py_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2y_vrr-adz*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_D2y_vrr = PBY*I_KINETIC_Fy2z_Py_vrr+oned2z*I_KINETIC_D2z_Py_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_D2y_vrr = PBY*I_KINETIC_F3z_Py_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2y_vrr-adz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_F3x_Dyz_vrr = PBZ*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_KINETIC_F2xy_Dyz_vrr = PBZ*I_KINETIC_F2xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dyz_vrr;
    Double I_KINETIC_F2xz_Dyz_vrr = PBZ*I_KINETIC_F2xz_Py_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dyz_vrr;
    Double I_KINETIC_Fx2y_Dyz_vrr = PBZ*I_KINETIC_Fx2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dyz_vrr;
    Double I_KINETIC_Fxyz_Dyz_vrr = PBZ*I_KINETIC_Fxyz_Py_vrr+oned2z*I_KINETIC_Dxy_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_Dyz_vrr;
    Double I_KINETIC_Fx2z_Dyz_vrr = PBZ*I_KINETIC_Fx2z_Py_vrr+2*oned2z*I_KINETIC_Dxz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dyz_vrr;
    Double I_KINETIC_F3y_Dyz_vrr = PBZ*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_KINETIC_F2yz_Dyz_vrr = PBZ*I_KINETIC_F2yz_Py_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dyz_vrr;
    Double I_KINETIC_Fy2z_Dyz_vrr = PBZ*I_KINETIC_Fy2z_Py_vrr+2*oned2z*I_KINETIC_Dyz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_Dyz_vrr;
    Double I_KINETIC_F3z_Dyz_vrr = PBZ*I_KINETIC_F3z_Py_vrr+3*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_KINETIC_F3x_D2z_vrr = PBZ*I_KINETIC_F3x_Pz_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2z_vrr-adz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_F2xy_D2z_vrr = PBZ*I_KINETIC_F2xy_Pz_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2z_vrr-adz*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_D2z_vrr = PBZ*I_KINETIC_F2xz_Pz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2z_vrr-adz*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_D2z_vrr = PBZ*I_KINETIC_Fx2y_Pz_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2z_vrr-adz*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_D2z_vrr = PBZ*I_KINETIC_Fxyz_Pz_vrr+oned2z*I_KINETIC_Dxy_Pz_vrr+oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_D2z_vrr-adz*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_D2z_vrr = PBZ*I_KINETIC_Fx2z_Pz_vrr+2*oned2z*I_KINETIC_Dxz_Pz_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_D2z_vrr = PBZ*I_KINETIC_F3y_Pz_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2z_vrr-adz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_F2yz_D2z_vrr = PBZ*I_KINETIC_F2yz_Pz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2z_vrr-adz*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_D2z_vrr = PBZ*I_KINETIC_Fy2z_Pz_vrr+2*oned2z*I_KINETIC_Dyz_Pz_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_D2z_vrr = PBZ*I_KINETIC_F3z_Pz_vrr+3*oned2z*I_KINETIC_D2z_Pz_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2z_vrr-adz*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 9 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_KINETIC_D_P
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_KINETIC_G4x_Px_vrr = PAX*I_KINETIC_F3x_Px_vrr+3*oned2z*I_KINETIC_D2x_Px_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Px_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_G3xy_Px_vrr = PAY*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_KINETIC_G3xz_Px_vrr = PAZ*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_KINETIC_G2x2y_Px_vrr = PAY*I_KINETIC_F2xy_Px_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Px_vrr-bdz*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_G2x2z_Px_vrr = PAZ*I_KINETIC_F2xz_Px_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Px_vrr-bdz*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_KINETIC_Gx3y_Px_vrr = PAX*I_KINETIC_F3y_Px_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_KINETIC_Gx3z_Px_vrr = PAX*I_KINETIC_F3z_Px_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_KINETIC_G4y_Px_vrr = PAY*I_KINETIC_F3y_Px_vrr+3*oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Px_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_G3yz_Px_vrr = PAZ*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_KINETIC_G2y2z_Px_vrr = PAZ*I_KINETIC_F2yz_Px_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Px_vrr-bdz*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_KINETIC_Gy3z_Px_vrr = PAY*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_KINETIC_G4z_Px_vrr = PAZ*I_KINETIC_F3z_Px_vrr+3*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Px_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_KINETIC_G4x_Py_vrr = PAX*I_KINETIC_F3x_Py_vrr+3*oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Py_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_G3xy_Py_vrr = PAY*I_KINETIC_F3x_Py_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_KINETIC_G3xz_Py_vrr = PAZ*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_KINETIC_G2x2y_Py_vrr = PAY*I_KINETIC_F2xy_Py_vrr+oned2z*I_KINETIC_D2x_Py_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Py_vrr-bdz*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_G2x2z_Py_vrr = PAZ*I_KINETIC_F2xz_Py_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Py_vrr-bdz*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_KINETIC_Gx3y_Py_vrr = PAX*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_KINETIC_Gx3z_Py_vrr = PAX*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_KINETIC_G4y_Py_vrr = PAY*I_KINETIC_F3y_Py_vrr+3*oned2z*I_KINETIC_D2y_Py_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Py_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_G3yz_Py_vrr = PAZ*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_KINETIC_G2y2z_Py_vrr = PAZ*I_KINETIC_F2yz_Py_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Py_vrr-bdz*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_KINETIC_Gy3z_Py_vrr = PAY*I_KINETIC_F3z_Py_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_KINETIC_G4z_Py_vrr = PAZ*I_KINETIC_F3z_Py_vrr+3*oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Py_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_KINETIC_G4x_Pz_vrr = PAX*I_KINETIC_F3x_Pz_vrr+3*oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Pz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_G3xy_Pz_vrr = PAY*I_KINETIC_F3x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_KINETIC_G3xz_Pz_vrr = PAZ*I_KINETIC_F3x_Pz_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_KINETIC_G2x2y_Pz_vrr = PAY*I_KINETIC_F2xy_Pz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Pz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_G2x2z_Pz_vrr = PAZ*I_KINETIC_F2xz_Pz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_KINETIC_Gx3y_Pz_vrr = PAX*I_KINETIC_F3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_KINETIC_Gx3z_Pz_vrr = PAX*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_KINETIC_G4y_Pz_vrr = PAY*I_KINETIC_F3y_Pz_vrr+3*oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Pz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_G3yz_Pz_vrr = PAZ*I_KINETIC_F3y_Pz_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_KINETIC_G2y2z_Pz_vrr = PAZ*I_KINETIC_F2yz_Pz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_KINETIC_Gy3z_Pz_vrr = PAY*I_KINETIC_F3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_KINETIC_G4z_Pz_vrr = PAZ*I_KINETIC_F3z_Pz_vrr+3*oned2z*I_KINETIC_D2z_Pz_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Pz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 18 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_KINETIC_F_P
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_KINETIC_G4x_D2x_vrr = PBX*I_KINETIC_G4x_Px_vrr+4*oned2z*I_KINETIC_F3x_Px_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2x_vrr-adz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_G3xy_D2x_vrr = PBX*I_KINETIC_G3xy_Px_vrr+3*oned2z*I_KINETIC_F2xy_Px_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2x_vrr-adz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_D2x_vrr = PBX*I_KINETIC_G3xz_Px_vrr+3*oned2z*I_KINETIC_F2xz_Px_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2x_vrr-adz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_D2x_vrr = PBX*I_KINETIC_G2x2y_Px_vrr+2*oned2z*I_KINETIC_Fx2y_Px_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2x_vrr-adz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_G2x2z_D2x_vrr = PBX*I_KINETIC_G2x2z_Px_vrr+2*oned2z*I_KINETIC_Fx2z_Px_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2x_vrr-adz*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_KINETIC_Gx3y_D2x_vrr = PBX*I_KINETIC_Gx3y_Px_vrr+oned2z*I_KINETIC_F3y_Px_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2x_vrr-adz*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx3z_D2x_vrr = PBX*I_KINETIC_Gx3z_Px_vrr+oned2z*I_KINETIC_F3z_Px_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2x_vrr-adz*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_D2x_vrr = PBX*I_KINETIC_G4y_Px_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2x_vrr-adz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_G3yz_D2x_vrr = PBX*I_KINETIC_G3yz_Px_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2x_vrr-adz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_D2x_vrr = PBX*I_KINETIC_G2y2z_Px_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2x_vrr-adz*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_KINETIC_Gy3z_D2x_vrr = PBX*I_KINETIC_Gy3z_Px_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2x_vrr-adz*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_D2x_vrr = PBX*I_KINETIC_G4z_Px_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2x_vrr-adz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_G4x_Dxy_vrr = PBY*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_KINETIC_G3xy_Dxy_vrr = PBY*I_KINETIC_G3xy_Px_vrr+oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Dxy_vrr;
    Double I_KINETIC_G3xz_Dxy_vrr = PBY*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Dxy_vrr;
    Double I_KINETIC_G2x2y_Dxy_vrr = PBY*I_KINETIC_G2x2y_Px_vrr+2*oned2z*I_KINETIC_F2xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_KINETIC_G2x2z_Dxy_vrr = PBY*I_KINETIC_G2x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dxy_vrr;
    Double I_KINETIC_Gx3y_Dxy_vrr = PBY*I_KINETIC_Gx3y_Px_vrr+3*oned2z*I_KINETIC_Fx2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr;
    Double I_KINETIC_Gx3z_Dxy_vrr = PBY*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr;
    Double I_KINETIC_G4y_Dxy_vrr = PBY*I_KINETIC_G4y_Px_vrr+4*oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_KINETIC_G3yz_Dxy_vrr = PBY*I_KINETIC_G3yz_Px_vrr+3*oned2z*I_KINETIC_F2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dxy_vrr;
    Double I_KINETIC_G2y2z_Dxy_vrr = PBY*I_KINETIC_G2y2z_Px_vrr+2*oned2z*I_KINETIC_Fy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dxy_vrr;
    Double I_KINETIC_Gy3z_Dxy_vrr = PBY*I_KINETIC_Gy3z_Px_vrr+oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr;
    Double I_KINETIC_G4z_Dxy_vrr = PBY*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_KINETIC_G4x_Dxz_vrr = PBZ*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Dxz_vrr;
    Double I_KINETIC_G3xy_Dxz_vrr = PBZ*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Dxz_vrr;
    Double I_KINETIC_G3xz_Dxz_vrr = PBZ*I_KINETIC_G3xz_Px_vrr+oned2z*I_KINETIC_F3x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Dxz_vrr;
    Double I_KINETIC_G2x2y_Dxz_vrr = PBZ*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Dxz_vrr;
    Double I_KINETIC_G2x2z_Dxz_vrr = PBZ*I_KINETIC_G2x2z_Px_vrr+2*oned2z*I_KINETIC_F2xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dxz_vrr;
    Double I_KINETIC_Gx3y_Dxz_vrr = PBZ*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dxz_vrr;
    Double I_KINETIC_Gx3z_Dxz_vrr = PBZ*I_KINETIC_Gx3z_Px_vrr+3*oned2z*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dxz_vrr;
    Double I_KINETIC_G4y_Dxz_vrr = PBZ*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dxz_vrr;
    Double I_KINETIC_G3yz_Dxz_vrr = PBZ*I_KINETIC_G3yz_Px_vrr+oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dxz_vrr;
    Double I_KINETIC_G2y2z_Dxz_vrr = PBZ*I_KINETIC_G2y2z_Px_vrr+2*oned2z*I_KINETIC_F2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dxz_vrr;
    Double I_KINETIC_Gy3z_Dxz_vrr = PBZ*I_KINETIC_Gy3z_Px_vrr+3*oned2z*I_KINETIC_Fy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dxz_vrr;
    Double I_KINETIC_G4z_Dxz_vrr = PBZ*I_KINETIC_G4z_Px_vrr+4*oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dxz_vrr;
    Double I_KINETIC_G4x_D2y_vrr = PBY*I_KINETIC_G4x_Py_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2y_vrr-adz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_G3xy_D2y_vrr = PBY*I_KINETIC_G3xy_Py_vrr+oned2z*I_KINETIC_F3x_Py_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2y_vrr-adz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_D2y_vrr = PBY*I_KINETIC_G3xz_Py_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2y_vrr-adz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_D2y_vrr = PBY*I_KINETIC_G2x2y_Py_vrr+2*oned2z*I_KINETIC_F2xy_Py_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2y_vrr-adz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_G2x2z_D2y_vrr = PBY*I_KINETIC_G2x2z_Py_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2y_vrr-adz*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_KINETIC_Gx3y_D2y_vrr = PBY*I_KINETIC_Gx3y_Py_vrr+3*oned2z*I_KINETIC_Fx2y_Py_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2y_vrr-adz*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx3z_D2y_vrr = PBY*I_KINETIC_Gx3z_Py_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2y_vrr-adz*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_D2y_vrr = PBY*I_KINETIC_G4y_Py_vrr+4*oned2z*I_KINETIC_F3y_Py_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2y_vrr-adz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_G3yz_D2y_vrr = PBY*I_KINETIC_G3yz_Py_vrr+3*oned2z*I_KINETIC_F2yz_Py_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2y_vrr-adz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_D2y_vrr = PBY*I_KINETIC_G2y2z_Py_vrr+2*oned2z*I_KINETIC_Fy2z_Py_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2y_vrr-adz*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_KINETIC_Gy3z_D2y_vrr = PBY*I_KINETIC_Gy3z_Py_vrr+oned2z*I_KINETIC_F3z_Py_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2y_vrr-adz*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_D2y_vrr = PBY*I_KINETIC_G4z_Py_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2y_vrr-adz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_G4x_Dyz_vrr = PBZ*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Dyz_vrr;
    Double I_KINETIC_G3xy_Dyz_vrr = PBZ*I_KINETIC_G3xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Dyz_vrr;
    Double I_KINETIC_G3xz_Dyz_vrr = PBZ*I_KINETIC_G3xz_Py_vrr+oned2z*I_KINETIC_F3x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Dyz_vrr;
    Double I_KINETIC_G2x2y_Dyz_vrr = PBZ*I_KINETIC_G2x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Dyz_vrr;
    Double I_KINETIC_G2x2z_Dyz_vrr = PBZ*I_KINETIC_G2x2z_Py_vrr+2*oned2z*I_KINETIC_F2xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dyz_vrr;
    Double I_KINETIC_Gx3y_Dyz_vrr = PBZ*I_KINETIC_Gx3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dyz_vrr;
    Double I_KINETIC_Gx3z_Dyz_vrr = PBZ*I_KINETIC_Gx3z_Py_vrr+3*oned2z*I_KINETIC_Fx2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dyz_vrr;
    Double I_KINETIC_G4y_Dyz_vrr = PBZ*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dyz_vrr;
    Double I_KINETIC_G3yz_Dyz_vrr = PBZ*I_KINETIC_G3yz_Py_vrr+oned2z*I_KINETIC_F3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dyz_vrr;
    Double I_KINETIC_G2y2z_Dyz_vrr = PBZ*I_KINETIC_G2y2z_Py_vrr+2*oned2z*I_KINETIC_F2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dyz_vrr;
    Double I_KINETIC_Gy3z_Dyz_vrr = PBZ*I_KINETIC_Gy3z_Py_vrr+3*oned2z*I_KINETIC_Fy2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dyz_vrr;
    Double I_KINETIC_G4z_Dyz_vrr = PBZ*I_KINETIC_G4z_Py_vrr+4*oned2z*I_KINETIC_F3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dyz_vrr;
    Double I_KINETIC_G4x_D2z_vrr = PBZ*I_KINETIC_G4x_Pz_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2z_vrr-adz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_G3xy_D2z_vrr = PBZ*I_KINETIC_G3xy_Pz_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2z_vrr-adz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_D2z_vrr = PBZ*I_KINETIC_G3xz_Pz_vrr+oned2z*I_KINETIC_F3x_Pz_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2z_vrr-adz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_D2z_vrr = PBZ*I_KINETIC_G2x2y_Pz_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2z_vrr-adz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_G2x2z_D2z_vrr = PBZ*I_KINETIC_G2x2z_Pz_vrr+2*oned2z*I_KINETIC_F2xz_Pz_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2z_vrr-adz*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_KINETIC_Gx3y_D2z_vrr = PBZ*I_KINETIC_Gx3y_Pz_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2z_vrr-adz*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx3z_D2z_vrr = PBZ*I_KINETIC_Gx3z_Pz_vrr+3*oned2z*I_KINETIC_Fx2z_Pz_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2z_vrr-adz*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_D2z_vrr = PBZ*I_KINETIC_G4y_Pz_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2z_vrr-adz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_G3yz_D2z_vrr = PBZ*I_KINETIC_G3yz_Pz_vrr+oned2z*I_KINETIC_F3y_Pz_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2z_vrr-adz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_D2z_vrr = PBZ*I_KINETIC_G2y2z_Pz_vrr+2*oned2z*I_KINETIC_F2yz_Pz_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2z_vrr-adz*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_KINETIC_Gy3z_D2z_vrr = PBZ*I_KINETIC_Gy3z_Pz_vrr+3*oned2z*I_KINETIC_Fy2z_Pz_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2z_vrr-adz*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_D2z_vrr = PBZ*I_KINETIC_G4z_Pz_vrr+4*oned2z*I_KINETIC_F3z_Pz_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2z_vrr-adz*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_D
     * RHS shell quartet name: SQ_KINETIC_F_D
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     ************************************************************/
    Double I_KINETIC_H5x_D2x_vrr = PAX*I_KINETIC_G4x_D2x_vrr+4*oned2z*I_KINETIC_F3x_D2x_vrr+2*oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2x_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_KINETIC_H4xy_D2x_vrr = PAY*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_KINETIC_H4xz_D2x_vrr = PAZ*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_KINETIC_H3x2y_D2x_vrr = PAY*I_KINETIC_G3xy_D2x_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_KINETIC_H3xyz_D2x_vrr = PAZ*I_KINETIC_G3xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_D2x_vrr;
    Double I_KINETIC_H3x2z_D2x_vrr = PAZ*I_KINETIC_G3xz_D2x_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_KINETIC_H2x3y_D2x_vrr = PAX*I_KINETIC_Gx3y_D2x_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+2*oned2z*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_KINETIC_H2x2yz_D2x_vrr = PAZ*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr;
    Double I_KINETIC_H2xy2z_D2x_vrr = PAY*I_KINETIC_G2x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr;
    Double I_KINETIC_H2x3z_D2x_vrr = PAX*I_KINETIC_Gx3z_D2x_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+2*oned2z*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_KINETIC_Hx4y_D2x_vrr = PAX*I_KINETIC_G4y_D2x_vrr+2*oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_KINETIC_Hx3yz_D2x_vrr = PAZ*I_KINETIC_Gx3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr;
    Double I_KINETIC_Hx2y2z_D2x_vrr = PAX*I_KINETIC_G2y2z_D2x_vrr+2*oned2z*I_KINETIC_G2y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr;
    Double I_KINETIC_Hxy3z_D2x_vrr = PAY*I_KINETIC_Gx3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr;
    Double I_KINETIC_Hx4z_D2x_vrr = PAX*I_KINETIC_G4z_D2x_vrr+2*oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_KINETIC_H5y_D2x_vrr = PAY*I_KINETIC_G4y_D2x_vrr+4*oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2x_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_KINETIC_H4yz_D2x_vrr = PAZ*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_KINETIC_H3y2z_D2x_vrr = PAZ*I_KINETIC_G3yz_D2x_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_KINETIC_H2y3z_D2x_vrr = PAY*I_KINETIC_Gy3z_D2x_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2x_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_KINETIC_Hy4z_D2x_vrr = PAY*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_KINETIC_H5z_D2x_vrr = PAZ*I_KINETIC_G4z_D2x_vrr+4*oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2x_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_KINETIC_H5x_Dxy_vrr = PAX*I_KINETIC_G4x_Dxy_vrr+4*oned2z*I_KINETIC_F3x_Dxy_vrr+oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dxy_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_KINETIC_H4xy_Dxy_vrr = PAY*I_KINETIC_G4x_Dxy_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dxy_vrr;
    Double I_KINETIC_H4xz_Dxy_vrr = PAZ*I_KINETIC_G4x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dxy_vrr;
    Double I_KINETIC_H3x2y_Dxy_vrr = PAY*I_KINETIC_G3xy_Dxy_vrr+oned2z*I_KINETIC_F3x_Dxy_vrr+oned2z*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_KINETIC_H3xyz_Dxy_vrr = PAZ*I_KINETIC_G3xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Dxy_vrr;
    Double I_KINETIC_H3x2z_Dxy_vrr = PAZ*I_KINETIC_G3xz_Dxy_vrr+oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_KINETIC_H2x3y_Dxy_vrr = PAX*I_KINETIC_Gx3y_Dxy_vrr+oned2z*I_KINETIC_F3y_Dxy_vrr+oned2z*I_KINETIC_Gx3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_KINETIC_H2x2yz_Dxy_vrr = PAZ*I_KINETIC_G2x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dxy_vrr;
    Double I_KINETIC_H2xy2z_Dxy_vrr = PAY*I_KINETIC_G2x2z_Dxy_vrr+oned2z*I_KINETIC_G2x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Dxy_vrr;
    Double I_KINETIC_H2x3z_Dxy_vrr = PAX*I_KINETIC_Gx3z_Dxy_vrr+oned2z*I_KINETIC_F3z_Dxy_vrr+oned2z*I_KINETIC_Gx3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_KINETIC_Hx4y_Dxy_vrr = PAX*I_KINETIC_G4y_Dxy_vrr+oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dxy_vrr;
    Double I_KINETIC_Hx3yz_Dxy_vrr = PAZ*I_KINETIC_Gx3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Dxy_vrr;
    Double I_KINETIC_Hx2y2z_Dxy_vrr = PAX*I_KINETIC_G2y2z_Dxy_vrr+oned2z*I_KINETIC_G2y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Dxy_vrr;
    Double I_KINETIC_Hxy3z_Dxy_vrr = PAY*I_KINETIC_Gx3z_Dxy_vrr+oned2z*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Dxy_vrr;
    Double I_KINETIC_Hx4z_Dxy_vrr = PAX*I_KINETIC_G4z_Dxy_vrr+oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dxy_vrr;
    Double I_KINETIC_H5y_Dxy_vrr = PAY*I_KINETIC_G4y_Dxy_vrr+4*oned2z*I_KINETIC_F3y_Dxy_vrr+oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dxy_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_KINETIC_H4yz_Dxy_vrr = PAZ*I_KINETIC_G4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dxy_vrr;
    Double I_KINETIC_H3y2z_Dxy_vrr = PAZ*I_KINETIC_G3yz_Dxy_vrr+oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_KINETIC_H2y3z_Dxy_vrr = PAY*I_KINETIC_Gy3z_Dxy_vrr+oned2z*I_KINETIC_F3z_Dxy_vrr+oned2z*I_KINETIC_Gy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_KINETIC_Hy4z_Dxy_vrr = PAY*I_KINETIC_G4z_Dxy_vrr+oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dxy_vrr;
    Double I_KINETIC_H5z_Dxy_vrr = PAZ*I_KINETIC_G4z_Dxy_vrr+4*oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dxy_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_KINETIC_H5x_Dxz_vrr = PAX*I_KINETIC_G4x_Dxz_vrr+4*oned2z*I_KINETIC_F3x_Dxz_vrr+oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dxz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_KINETIC_H4xy_Dxz_vrr = PAY*I_KINETIC_G4x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dxz_vrr;
    Double I_KINETIC_H4xz_Dxz_vrr = PAZ*I_KINETIC_G4x_Dxz_vrr+oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dxz_vrr;
    Double I_KINETIC_H3x2y_Dxz_vrr = PAY*I_KINETIC_G3xy_Dxz_vrr+oned2z*I_KINETIC_F3x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_KINETIC_H3xyz_Dxz_vrr = PAZ*I_KINETIC_G3xy_Dxz_vrr+oned2z*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Dxz_vrr;
    Double I_KINETIC_H3x2z_Dxz_vrr = PAZ*I_KINETIC_G3xz_Dxz_vrr+oned2z*I_KINETIC_F3x_Dxz_vrr+oned2z*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_KINETIC_H2x3y_Dxz_vrr = PAX*I_KINETIC_Gx3y_Dxz_vrr+oned2z*I_KINETIC_F3y_Dxz_vrr+oned2z*I_KINETIC_Gx3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_KINETIC_H2x2yz_Dxz_vrr = PAZ*I_KINETIC_G2x2y_Dxz_vrr+oned2z*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dxz_vrr;
    Double I_KINETIC_H2xy2z_Dxz_vrr = PAY*I_KINETIC_G2x2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Dxz_vrr;
    Double I_KINETIC_H2x3z_Dxz_vrr = PAX*I_KINETIC_Gx3z_Dxz_vrr+oned2z*I_KINETIC_F3z_Dxz_vrr+oned2z*I_KINETIC_Gx3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_KINETIC_Hx4y_Dxz_vrr = PAX*I_KINETIC_G4y_Dxz_vrr+oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dxz_vrr;
    Double I_KINETIC_Hx3yz_Dxz_vrr = PAZ*I_KINETIC_Gx3y_Dxz_vrr+oned2z*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Dxz_vrr;
    Double I_KINETIC_Hx2y2z_Dxz_vrr = PAX*I_KINETIC_G2y2z_Dxz_vrr+oned2z*I_KINETIC_G2y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Dxz_vrr;
    Double I_KINETIC_Hxy3z_Dxz_vrr = PAY*I_KINETIC_Gx3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Dxz_vrr;
    Double I_KINETIC_Hx4z_Dxz_vrr = PAX*I_KINETIC_G4z_Dxz_vrr+oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dxz_vrr;
    Double I_KINETIC_H5y_Dxz_vrr = PAY*I_KINETIC_G4y_Dxz_vrr+4*oned2z*I_KINETIC_F3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dxz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_KINETIC_H4yz_Dxz_vrr = PAZ*I_KINETIC_G4y_Dxz_vrr+oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dxz_vrr;
    Double I_KINETIC_H3y2z_Dxz_vrr = PAZ*I_KINETIC_G3yz_Dxz_vrr+oned2z*I_KINETIC_F3y_Dxz_vrr+oned2z*I_KINETIC_G3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_KINETIC_H2y3z_Dxz_vrr = PAY*I_KINETIC_Gy3z_Dxz_vrr+oned2z*I_KINETIC_F3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_KINETIC_Hy4z_Dxz_vrr = PAY*I_KINETIC_G4z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dxz_vrr;
    Double I_KINETIC_H5z_Dxz_vrr = PAZ*I_KINETIC_G4z_Dxz_vrr+4*oned2z*I_KINETIC_F3z_Dxz_vrr+oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dxz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_KINETIC_H5x_D2y_vrr = PAX*I_KINETIC_G4x_D2y_vrr+4*oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_KINETIC_H4xy_D2y_vrr = PAY*I_KINETIC_G4x_D2y_vrr+2*oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_KINETIC_H4xz_D2y_vrr = PAZ*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2y_vrr;
    Double I_KINETIC_H3x2y_D2y_vrr = PAY*I_KINETIC_G3xy_D2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+2*oned2z*I_KINETIC_G3xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_KINETIC_H3xyz_D2y_vrr = PAZ*I_KINETIC_G3xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_D2y_vrr;
    Double I_KINETIC_H3x2z_D2y_vrr = PAZ*I_KINETIC_G3xz_D2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_KINETIC_H2x3y_D2y_vrr = PAX*I_KINETIC_Gx3y_D2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_KINETIC_H2x2yz_D2y_vrr = PAZ*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr;
    Double I_KINETIC_H2xy2z_D2y_vrr = PAY*I_KINETIC_G2x2z_D2y_vrr+2*oned2z*I_KINETIC_G2x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr;
    Double I_KINETIC_H2x3z_D2y_vrr = PAX*I_KINETIC_Gx3z_D2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_KINETIC_Hx4y_D2y_vrr = PAX*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_KINETIC_Hx3yz_D2y_vrr = PAZ*I_KINETIC_Gx3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr;
    Double I_KINETIC_Hx2y2z_D2y_vrr = PAX*I_KINETIC_G2y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr;
    Double I_KINETIC_Hxy3z_D2y_vrr = PAY*I_KINETIC_Gx3z_D2y_vrr+2*oned2z*I_KINETIC_Gx3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr;
    Double I_KINETIC_Hx4z_D2y_vrr = PAX*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_KINETIC_H5y_D2y_vrr = PAY*I_KINETIC_G4y_D2y_vrr+4*oned2z*I_KINETIC_F3y_D2y_vrr+2*oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_KINETIC_H4yz_D2y_vrr = PAZ*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_KINETIC_H3y2z_D2y_vrr = PAZ*I_KINETIC_G3yz_D2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_KINETIC_H2y3z_D2y_vrr = PAY*I_KINETIC_Gy3z_D2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+2*oned2z*I_KINETIC_Gy3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2y_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_KINETIC_Hy4z_D2y_vrr = PAY*I_KINETIC_G4z_D2y_vrr+2*oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_KINETIC_H5z_D2y_vrr = PAZ*I_KINETIC_G4z_D2y_vrr+4*oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2y_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_KINETIC_H5x_Dyz_vrr = PAX*I_KINETIC_G4x_Dyz_vrr+4*oned2z*I_KINETIC_F3x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Dyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_KINETIC_H4xy_Dyz_vrr = PAY*I_KINETIC_G4x_Dyz_vrr+oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Dyz_vrr;
    Double I_KINETIC_H4xz_Dyz_vrr = PAZ*I_KINETIC_G4x_Dyz_vrr+oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Dyz_vrr;
    Double I_KINETIC_H3x2y_Dyz_vrr = PAY*I_KINETIC_G3xy_Dyz_vrr+oned2z*I_KINETIC_F3x_Dyz_vrr+oned2z*I_KINETIC_G3xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_KINETIC_H3xyz_Dyz_vrr = PAZ*I_KINETIC_G3xy_Dyz_vrr+oned2z*I_KINETIC_G3xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Dyz_vrr;
    Double I_KINETIC_H3x2z_Dyz_vrr = PAZ*I_KINETIC_G3xz_Dyz_vrr+oned2z*I_KINETIC_F3x_Dyz_vrr+oned2z*I_KINETIC_G3xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_KINETIC_H2x3y_Dyz_vrr = PAX*I_KINETIC_Gx3y_Dyz_vrr+oned2z*I_KINETIC_F3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_KINETIC_H2x2yz_Dyz_vrr = PAZ*I_KINETIC_G2x2y_Dyz_vrr+oned2z*I_KINETIC_G2x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Dyz_vrr;
    Double I_KINETIC_H2xy2z_Dyz_vrr = PAY*I_KINETIC_G2x2z_Dyz_vrr+oned2z*I_KINETIC_G2x2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Dyz_vrr;
    Double I_KINETIC_H2x3z_Dyz_vrr = PAX*I_KINETIC_Gx3z_Dyz_vrr+oned2z*I_KINETIC_F3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_KINETIC_Hx4y_Dyz_vrr = PAX*I_KINETIC_G4y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Dyz_vrr;
    Double I_KINETIC_Hx3yz_Dyz_vrr = PAZ*I_KINETIC_Gx3y_Dyz_vrr+oned2z*I_KINETIC_Gx3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Dyz_vrr;
    Double I_KINETIC_Hx2y2z_Dyz_vrr = PAX*I_KINETIC_G2y2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Dyz_vrr;
    Double I_KINETIC_Hxy3z_Dyz_vrr = PAY*I_KINETIC_Gx3z_Dyz_vrr+oned2z*I_KINETIC_Gx3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Dyz_vrr;
    Double I_KINETIC_Hx4z_Dyz_vrr = PAX*I_KINETIC_G4z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Dyz_vrr;
    Double I_KINETIC_H5y_Dyz_vrr = PAY*I_KINETIC_G4y_Dyz_vrr+4*oned2z*I_KINETIC_F3y_Dyz_vrr+oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Dyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_KINETIC_H4yz_Dyz_vrr = PAZ*I_KINETIC_G4y_Dyz_vrr+oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Dyz_vrr;
    Double I_KINETIC_H3y2z_Dyz_vrr = PAZ*I_KINETIC_G3yz_Dyz_vrr+oned2z*I_KINETIC_F3y_Dyz_vrr+oned2z*I_KINETIC_G3yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_KINETIC_H2y3z_Dyz_vrr = PAY*I_KINETIC_Gy3z_Dyz_vrr+oned2z*I_KINETIC_F3z_Dyz_vrr+oned2z*I_KINETIC_Gy3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_KINETIC_Hy4z_Dyz_vrr = PAY*I_KINETIC_G4z_Dyz_vrr+oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Dyz_vrr;
    Double I_KINETIC_H5z_Dyz_vrr = PAZ*I_KINETIC_G4z_Dyz_vrr+4*oned2z*I_KINETIC_F3z_Dyz_vrr+oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Dyz_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_KINETIC_H5x_D2z_vrr = PAX*I_KINETIC_G4x_D2z_vrr+4*oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_KINETIC_H4xy_D2z_vrr = PAY*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2z_vrr;
    Double I_KINETIC_H4xz_D2z_vrr = PAZ*I_KINETIC_G4x_D2z_vrr+2*oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2z_vrr;
    Double I_KINETIC_H3x2y_D2z_vrr = PAY*I_KINETIC_G3xy_D2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_KINETIC_H3xyz_D2z_vrr = PAZ*I_KINETIC_G3xy_D2z_vrr+2*oned2z*I_KINETIC_G3xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_D2z_vrr;
    Double I_KINETIC_H3x2z_D2z_vrr = PAZ*I_KINETIC_G3xz_D2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+2*oned2z*I_KINETIC_G3xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_KINETIC_H2x3y_D2z_vrr = PAX*I_KINETIC_Gx3y_D2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_KINETIC_H2x2yz_D2z_vrr = PAZ*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_G2x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr;
    Double I_KINETIC_H2xy2z_D2z_vrr = PAY*I_KINETIC_G2x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_D2z_vrr;
    Double I_KINETIC_H2x3z_D2z_vrr = PAX*I_KINETIC_Gx3z_D2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_KINETIC_Hx4y_D2z_vrr = PAX*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2z_vrr;
    Double I_KINETIC_Hx3yz_D2z_vrr = PAZ*I_KINETIC_Gx3y_D2z_vrr+2*oned2z*I_KINETIC_Gx3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_D2z_vrr;
    Double I_KINETIC_Hx2y2z_D2z_vrr = PAX*I_KINETIC_G2y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_D2z_vrr;
    Double I_KINETIC_Hxy3z_D2z_vrr = PAY*I_KINETIC_Gx3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_D2z_vrr;
    Double I_KINETIC_Hx4z_D2z_vrr = PAX*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_KINETIC_H5y_D2z_vrr = PAY*I_KINETIC_G4y_D2z_vrr+4*oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_KINETIC_H4yz_D2z_vrr = PAZ*I_KINETIC_G4y_D2z_vrr+2*oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2z_vrr;
    Double I_KINETIC_H3y2z_D2z_vrr = PAZ*I_KINETIC_G3yz_D2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+2*oned2z*I_KINETIC_G3yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_KINETIC_H2y3z_D2z_vrr = PAY*I_KINETIC_Gy3z_D2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2z_vrr-bdz*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_KINETIC_Hy4z_D2z_vrr = PAY*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2z_vrr;
    Double I_KINETIC_H5z_D2z_vrr = PAZ*I_KINETIC_G4z_D2z_vrr+4*oned2z*I_KINETIC_F3z_D2z_vrr+2*oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2z_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_D_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_H_D_aa_coefs = alpha*alpha;
    I_KINETIC_H5x_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5x_D2x_vrr;
    I_KINETIC_H4xy_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xy_D2x_vrr;
    I_KINETIC_H4xz_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xz_D2x_vrr;
    I_KINETIC_H3x2y_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2y_D2x_vrr;
    I_KINETIC_H3xyz_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3xyz_D2x_vrr;
    I_KINETIC_H3x2z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2z_D2x_vrr;
    I_KINETIC_H2x3y_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3y_D2x_vrr;
    I_KINETIC_H2x2yz_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x2yz_D2x_vrr;
    I_KINETIC_H2xy2z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2xy2z_D2x_vrr;
    I_KINETIC_H2x3z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3z_D2x_vrr;
    I_KINETIC_Hx4y_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4y_D2x_vrr;
    I_KINETIC_Hx3yz_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx3yz_D2x_vrr;
    I_KINETIC_Hx2y2z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx2y2z_D2x_vrr;
    I_KINETIC_Hxy3z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hxy3z_D2x_vrr;
    I_KINETIC_Hx4z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4z_D2x_vrr;
    I_KINETIC_H5y_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5y_D2x_vrr;
    I_KINETIC_H4yz_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4yz_D2x_vrr;
    I_KINETIC_H3y2z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3y2z_D2x_vrr;
    I_KINETIC_H2y3z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2y3z_D2x_vrr;
    I_KINETIC_Hy4z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hy4z_D2x_vrr;
    I_KINETIC_H5z_D2x_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5z_D2x_vrr;
    I_KINETIC_H5x_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5x_Dxy_vrr;
    I_KINETIC_H4xy_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xy_Dxy_vrr;
    I_KINETIC_H4xz_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xz_Dxy_vrr;
    I_KINETIC_H3x2y_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2y_Dxy_vrr;
    I_KINETIC_H3xyz_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3xyz_Dxy_vrr;
    I_KINETIC_H3x2z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2z_Dxy_vrr;
    I_KINETIC_H2x3y_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3y_Dxy_vrr;
    I_KINETIC_H2x2yz_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x2yz_Dxy_vrr;
    I_KINETIC_H2xy2z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2xy2z_Dxy_vrr;
    I_KINETIC_H2x3z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3z_Dxy_vrr;
    I_KINETIC_Hx4y_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4y_Dxy_vrr;
    I_KINETIC_Hx3yz_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx3yz_Dxy_vrr;
    I_KINETIC_Hx2y2z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx2y2z_Dxy_vrr;
    I_KINETIC_Hxy3z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hxy3z_Dxy_vrr;
    I_KINETIC_Hx4z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4z_Dxy_vrr;
    I_KINETIC_H5y_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5y_Dxy_vrr;
    I_KINETIC_H4yz_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4yz_Dxy_vrr;
    I_KINETIC_H3y2z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3y2z_Dxy_vrr;
    I_KINETIC_H2y3z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2y3z_Dxy_vrr;
    I_KINETIC_Hy4z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hy4z_Dxy_vrr;
    I_KINETIC_H5z_Dxy_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5z_Dxy_vrr;
    I_KINETIC_H5x_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5x_Dxz_vrr;
    I_KINETIC_H4xy_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xy_Dxz_vrr;
    I_KINETIC_H4xz_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xz_Dxz_vrr;
    I_KINETIC_H3x2y_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2y_Dxz_vrr;
    I_KINETIC_H3xyz_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3xyz_Dxz_vrr;
    I_KINETIC_H3x2z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2z_Dxz_vrr;
    I_KINETIC_H2x3y_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3y_Dxz_vrr;
    I_KINETIC_H2x2yz_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x2yz_Dxz_vrr;
    I_KINETIC_H2xy2z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2xy2z_Dxz_vrr;
    I_KINETIC_H2x3z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3z_Dxz_vrr;
    I_KINETIC_Hx4y_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4y_Dxz_vrr;
    I_KINETIC_Hx3yz_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx3yz_Dxz_vrr;
    I_KINETIC_Hx2y2z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx2y2z_Dxz_vrr;
    I_KINETIC_Hxy3z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hxy3z_Dxz_vrr;
    I_KINETIC_Hx4z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4z_Dxz_vrr;
    I_KINETIC_H5y_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5y_Dxz_vrr;
    I_KINETIC_H4yz_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4yz_Dxz_vrr;
    I_KINETIC_H3y2z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3y2z_Dxz_vrr;
    I_KINETIC_H2y3z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2y3z_Dxz_vrr;
    I_KINETIC_Hy4z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hy4z_Dxz_vrr;
    I_KINETIC_H5z_Dxz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5z_Dxz_vrr;
    I_KINETIC_H5x_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5x_D2y_vrr;
    I_KINETIC_H4xy_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xy_D2y_vrr;
    I_KINETIC_H4xz_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xz_D2y_vrr;
    I_KINETIC_H3x2y_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2y_D2y_vrr;
    I_KINETIC_H3xyz_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3xyz_D2y_vrr;
    I_KINETIC_H3x2z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2z_D2y_vrr;
    I_KINETIC_H2x3y_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3y_D2y_vrr;
    I_KINETIC_H2x2yz_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x2yz_D2y_vrr;
    I_KINETIC_H2xy2z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2xy2z_D2y_vrr;
    I_KINETIC_H2x3z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3z_D2y_vrr;
    I_KINETIC_Hx4y_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4y_D2y_vrr;
    I_KINETIC_Hx3yz_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx3yz_D2y_vrr;
    I_KINETIC_Hx2y2z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx2y2z_D2y_vrr;
    I_KINETIC_Hxy3z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hxy3z_D2y_vrr;
    I_KINETIC_Hx4z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4z_D2y_vrr;
    I_KINETIC_H5y_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5y_D2y_vrr;
    I_KINETIC_H4yz_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4yz_D2y_vrr;
    I_KINETIC_H3y2z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3y2z_D2y_vrr;
    I_KINETIC_H2y3z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2y3z_D2y_vrr;
    I_KINETIC_Hy4z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hy4z_D2y_vrr;
    I_KINETIC_H5z_D2y_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5z_D2y_vrr;
    I_KINETIC_H5x_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5x_Dyz_vrr;
    I_KINETIC_H4xy_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xy_Dyz_vrr;
    I_KINETIC_H4xz_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xz_Dyz_vrr;
    I_KINETIC_H3x2y_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2y_Dyz_vrr;
    I_KINETIC_H3xyz_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3xyz_Dyz_vrr;
    I_KINETIC_H3x2z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2z_Dyz_vrr;
    I_KINETIC_H2x3y_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3y_Dyz_vrr;
    I_KINETIC_H2x2yz_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x2yz_Dyz_vrr;
    I_KINETIC_H2xy2z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2xy2z_Dyz_vrr;
    I_KINETIC_H2x3z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3z_Dyz_vrr;
    I_KINETIC_Hx4y_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4y_Dyz_vrr;
    I_KINETIC_Hx3yz_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx3yz_Dyz_vrr;
    I_KINETIC_Hx2y2z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx2y2z_Dyz_vrr;
    I_KINETIC_Hxy3z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hxy3z_Dyz_vrr;
    I_KINETIC_Hx4z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4z_Dyz_vrr;
    I_KINETIC_H5y_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5y_Dyz_vrr;
    I_KINETIC_H4yz_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4yz_Dyz_vrr;
    I_KINETIC_H3y2z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3y2z_Dyz_vrr;
    I_KINETIC_H2y3z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2y3z_Dyz_vrr;
    I_KINETIC_Hy4z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hy4z_Dyz_vrr;
    I_KINETIC_H5z_Dyz_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5z_Dyz_vrr;
    I_KINETIC_H5x_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5x_D2z_vrr;
    I_KINETIC_H4xy_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xy_D2z_vrr;
    I_KINETIC_H4xz_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4xz_D2z_vrr;
    I_KINETIC_H3x2y_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2y_D2z_vrr;
    I_KINETIC_H3xyz_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3xyz_D2z_vrr;
    I_KINETIC_H3x2z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3x2z_D2z_vrr;
    I_KINETIC_H2x3y_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3y_D2z_vrr;
    I_KINETIC_H2x2yz_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x2yz_D2z_vrr;
    I_KINETIC_H2xy2z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2xy2z_D2z_vrr;
    I_KINETIC_H2x3z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2x3z_D2z_vrr;
    I_KINETIC_Hx4y_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4y_D2z_vrr;
    I_KINETIC_Hx3yz_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx3yz_D2z_vrr;
    I_KINETIC_Hx2y2z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx2y2z_D2z_vrr;
    I_KINETIC_Hxy3z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hxy3z_D2z_vrr;
    I_KINETIC_Hx4z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hx4z_D2z_vrr;
    I_KINETIC_H5y_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5y_D2z_vrr;
    I_KINETIC_H4yz_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H4yz_D2z_vrr;
    I_KINETIC_H3y2z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H3y2z_D2z_vrr;
    I_KINETIC_H2y3z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H2y3z_D2z_vrr;
    I_KINETIC_Hy4z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_Hy4z_D2z_vrr;
    I_KINETIC_H5z_D2z_aa += SQ_KINETIC_H_D_aa_coefs*I_KINETIC_H5z_D2z_vrr;

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
   * shell quartet name: SQ_KINETIC_F_D_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_D_aa
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_P_D
   ************************************************************/
  abcd[0] = 4.0E0*I_KINETIC_H5x_D2x_aa-2.0E0*3*I_KINETIC_F3x_D2x_a-2.0E0*4*I_KINETIC_F3x_D2x_a+3*2*I_KINETIC_Px_D2x;
  abcd[1] = 4.0E0*I_KINETIC_H4xy_D2x_aa-2.0E0*2*I_KINETIC_F2xy_D2x_a-2.0E0*3*I_KINETIC_F2xy_D2x_a+2*1*I_KINETIC_Py_D2x;
  abcd[2] = 4.0E0*I_KINETIC_H4xz_D2x_aa-2.0E0*2*I_KINETIC_F2xz_D2x_a-2.0E0*3*I_KINETIC_F2xz_D2x_a+2*1*I_KINETIC_Pz_D2x;
  abcd[3] = 4.0E0*I_KINETIC_H3x2y_D2x_aa-2.0E0*1*I_KINETIC_Fx2y_D2x_a-2.0E0*2*I_KINETIC_Fx2y_D2x_a;
  abcd[4] = 4.0E0*I_KINETIC_H3xyz_D2x_aa-2.0E0*1*I_KINETIC_Fxyz_D2x_a-2.0E0*2*I_KINETIC_Fxyz_D2x_a;
  abcd[5] = 4.0E0*I_KINETIC_H3x2z_D2x_aa-2.0E0*1*I_KINETIC_Fx2z_D2x_a-2.0E0*2*I_KINETIC_Fx2z_D2x_a;
  abcd[6] = 4.0E0*I_KINETIC_H2x3y_D2x_aa-2.0E0*1*I_KINETIC_F3y_D2x_a;
  abcd[7] = 4.0E0*I_KINETIC_H2x2yz_D2x_aa-2.0E0*1*I_KINETIC_F2yz_D2x_a;
  abcd[8] = 4.0E0*I_KINETIC_H2xy2z_D2x_aa-2.0E0*1*I_KINETIC_Fy2z_D2x_a;
  abcd[9] = 4.0E0*I_KINETIC_H2x3z_D2x_aa-2.0E0*1*I_KINETIC_F3z_D2x_a;
  abcd[10] = 4.0E0*I_KINETIC_H5x_Dxy_aa-2.0E0*3*I_KINETIC_F3x_Dxy_a-2.0E0*4*I_KINETIC_F3x_Dxy_a+3*2*I_KINETIC_Px_Dxy;
  abcd[11] = 4.0E0*I_KINETIC_H4xy_Dxy_aa-2.0E0*2*I_KINETIC_F2xy_Dxy_a-2.0E0*3*I_KINETIC_F2xy_Dxy_a+2*1*I_KINETIC_Py_Dxy;
  abcd[12] = 4.0E0*I_KINETIC_H4xz_Dxy_aa-2.0E0*2*I_KINETIC_F2xz_Dxy_a-2.0E0*3*I_KINETIC_F2xz_Dxy_a+2*1*I_KINETIC_Pz_Dxy;
  abcd[13] = 4.0E0*I_KINETIC_H3x2y_Dxy_aa-2.0E0*1*I_KINETIC_Fx2y_Dxy_a-2.0E0*2*I_KINETIC_Fx2y_Dxy_a;
  abcd[14] = 4.0E0*I_KINETIC_H3xyz_Dxy_aa-2.0E0*1*I_KINETIC_Fxyz_Dxy_a-2.0E0*2*I_KINETIC_Fxyz_Dxy_a;
  abcd[15] = 4.0E0*I_KINETIC_H3x2z_Dxy_aa-2.0E0*1*I_KINETIC_Fx2z_Dxy_a-2.0E0*2*I_KINETIC_Fx2z_Dxy_a;
  abcd[16] = 4.0E0*I_KINETIC_H2x3y_Dxy_aa-2.0E0*1*I_KINETIC_F3y_Dxy_a;
  abcd[17] = 4.0E0*I_KINETIC_H2x2yz_Dxy_aa-2.0E0*1*I_KINETIC_F2yz_Dxy_a;
  abcd[18] = 4.0E0*I_KINETIC_H2xy2z_Dxy_aa-2.0E0*1*I_KINETIC_Fy2z_Dxy_a;
  abcd[19] = 4.0E0*I_KINETIC_H2x3z_Dxy_aa-2.0E0*1*I_KINETIC_F3z_Dxy_a;
  abcd[20] = 4.0E0*I_KINETIC_H5x_Dxz_aa-2.0E0*3*I_KINETIC_F3x_Dxz_a-2.0E0*4*I_KINETIC_F3x_Dxz_a+3*2*I_KINETIC_Px_Dxz;
  abcd[21] = 4.0E0*I_KINETIC_H4xy_Dxz_aa-2.0E0*2*I_KINETIC_F2xy_Dxz_a-2.0E0*3*I_KINETIC_F2xy_Dxz_a+2*1*I_KINETIC_Py_Dxz;
  abcd[22] = 4.0E0*I_KINETIC_H4xz_Dxz_aa-2.0E0*2*I_KINETIC_F2xz_Dxz_a-2.0E0*3*I_KINETIC_F2xz_Dxz_a+2*1*I_KINETIC_Pz_Dxz;
  abcd[23] = 4.0E0*I_KINETIC_H3x2y_Dxz_aa-2.0E0*1*I_KINETIC_Fx2y_Dxz_a-2.0E0*2*I_KINETIC_Fx2y_Dxz_a;
  abcd[24] = 4.0E0*I_KINETIC_H3xyz_Dxz_aa-2.0E0*1*I_KINETIC_Fxyz_Dxz_a-2.0E0*2*I_KINETIC_Fxyz_Dxz_a;
  abcd[25] = 4.0E0*I_KINETIC_H3x2z_Dxz_aa-2.0E0*1*I_KINETIC_Fx2z_Dxz_a-2.0E0*2*I_KINETIC_Fx2z_Dxz_a;
  abcd[26] = 4.0E0*I_KINETIC_H2x3y_Dxz_aa-2.0E0*1*I_KINETIC_F3y_Dxz_a;
  abcd[27] = 4.0E0*I_KINETIC_H2x2yz_Dxz_aa-2.0E0*1*I_KINETIC_F2yz_Dxz_a;
  abcd[28] = 4.0E0*I_KINETIC_H2xy2z_Dxz_aa-2.0E0*1*I_KINETIC_Fy2z_Dxz_a;
  abcd[29] = 4.0E0*I_KINETIC_H2x3z_Dxz_aa-2.0E0*1*I_KINETIC_F3z_Dxz_a;
  abcd[30] = 4.0E0*I_KINETIC_H5x_D2y_aa-2.0E0*3*I_KINETIC_F3x_D2y_a-2.0E0*4*I_KINETIC_F3x_D2y_a+3*2*I_KINETIC_Px_D2y;
  abcd[31] = 4.0E0*I_KINETIC_H4xy_D2y_aa-2.0E0*2*I_KINETIC_F2xy_D2y_a-2.0E0*3*I_KINETIC_F2xy_D2y_a+2*1*I_KINETIC_Py_D2y;
  abcd[32] = 4.0E0*I_KINETIC_H4xz_D2y_aa-2.0E0*2*I_KINETIC_F2xz_D2y_a-2.0E0*3*I_KINETIC_F2xz_D2y_a+2*1*I_KINETIC_Pz_D2y;
  abcd[33] = 4.0E0*I_KINETIC_H3x2y_D2y_aa-2.0E0*1*I_KINETIC_Fx2y_D2y_a-2.0E0*2*I_KINETIC_Fx2y_D2y_a;
  abcd[34] = 4.0E0*I_KINETIC_H3xyz_D2y_aa-2.0E0*1*I_KINETIC_Fxyz_D2y_a-2.0E0*2*I_KINETIC_Fxyz_D2y_a;
  abcd[35] = 4.0E0*I_KINETIC_H3x2z_D2y_aa-2.0E0*1*I_KINETIC_Fx2z_D2y_a-2.0E0*2*I_KINETIC_Fx2z_D2y_a;
  abcd[36] = 4.0E0*I_KINETIC_H2x3y_D2y_aa-2.0E0*1*I_KINETIC_F3y_D2y_a;
  abcd[37] = 4.0E0*I_KINETIC_H2x2yz_D2y_aa-2.0E0*1*I_KINETIC_F2yz_D2y_a;
  abcd[38] = 4.0E0*I_KINETIC_H2xy2z_D2y_aa-2.0E0*1*I_KINETIC_Fy2z_D2y_a;
  abcd[39] = 4.0E0*I_KINETIC_H2x3z_D2y_aa-2.0E0*1*I_KINETIC_F3z_D2y_a;
  abcd[40] = 4.0E0*I_KINETIC_H5x_Dyz_aa-2.0E0*3*I_KINETIC_F3x_Dyz_a-2.0E0*4*I_KINETIC_F3x_Dyz_a+3*2*I_KINETIC_Px_Dyz;
  abcd[41] = 4.0E0*I_KINETIC_H4xy_Dyz_aa-2.0E0*2*I_KINETIC_F2xy_Dyz_a-2.0E0*3*I_KINETIC_F2xy_Dyz_a+2*1*I_KINETIC_Py_Dyz;
  abcd[42] = 4.0E0*I_KINETIC_H4xz_Dyz_aa-2.0E0*2*I_KINETIC_F2xz_Dyz_a-2.0E0*3*I_KINETIC_F2xz_Dyz_a+2*1*I_KINETIC_Pz_Dyz;
  abcd[43] = 4.0E0*I_KINETIC_H3x2y_Dyz_aa-2.0E0*1*I_KINETIC_Fx2y_Dyz_a-2.0E0*2*I_KINETIC_Fx2y_Dyz_a;
  abcd[44] = 4.0E0*I_KINETIC_H3xyz_Dyz_aa-2.0E0*1*I_KINETIC_Fxyz_Dyz_a-2.0E0*2*I_KINETIC_Fxyz_Dyz_a;
  abcd[45] = 4.0E0*I_KINETIC_H3x2z_Dyz_aa-2.0E0*1*I_KINETIC_Fx2z_Dyz_a-2.0E0*2*I_KINETIC_Fx2z_Dyz_a;
  abcd[46] = 4.0E0*I_KINETIC_H2x3y_Dyz_aa-2.0E0*1*I_KINETIC_F3y_Dyz_a;
  abcd[47] = 4.0E0*I_KINETIC_H2x2yz_Dyz_aa-2.0E0*1*I_KINETIC_F2yz_Dyz_a;
  abcd[48] = 4.0E0*I_KINETIC_H2xy2z_Dyz_aa-2.0E0*1*I_KINETIC_Fy2z_Dyz_a;
  abcd[49] = 4.0E0*I_KINETIC_H2x3z_Dyz_aa-2.0E0*1*I_KINETIC_F3z_Dyz_a;
  abcd[50] = 4.0E0*I_KINETIC_H5x_D2z_aa-2.0E0*3*I_KINETIC_F3x_D2z_a-2.0E0*4*I_KINETIC_F3x_D2z_a+3*2*I_KINETIC_Px_D2z;
  abcd[51] = 4.0E0*I_KINETIC_H4xy_D2z_aa-2.0E0*2*I_KINETIC_F2xy_D2z_a-2.0E0*3*I_KINETIC_F2xy_D2z_a+2*1*I_KINETIC_Py_D2z;
  abcd[52] = 4.0E0*I_KINETIC_H4xz_D2z_aa-2.0E0*2*I_KINETIC_F2xz_D2z_a-2.0E0*3*I_KINETIC_F2xz_D2z_a+2*1*I_KINETIC_Pz_D2z;
  abcd[53] = 4.0E0*I_KINETIC_H3x2y_D2z_aa-2.0E0*1*I_KINETIC_Fx2y_D2z_a-2.0E0*2*I_KINETIC_Fx2y_D2z_a;
  abcd[54] = 4.0E0*I_KINETIC_H3xyz_D2z_aa-2.0E0*1*I_KINETIC_Fxyz_D2z_a-2.0E0*2*I_KINETIC_Fxyz_D2z_a;
  abcd[55] = 4.0E0*I_KINETIC_H3x2z_D2z_aa-2.0E0*1*I_KINETIC_Fx2z_D2z_a-2.0E0*2*I_KINETIC_Fx2z_D2z_a;
  abcd[56] = 4.0E0*I_KINETIC_H2x3y_D2z_aa-2.0E0*1*I_KINETIC_F3y_D2z_a;
  abcd[57] = 4.0E0*I_KINETIC_H2x2yz_D2z_aa-2.0E0*1*I_KINETIC_F2yz_D2z_a;
  abcd[58] = 4.0E0*I_KINETIC_H2xy2z_D2z_aa-2.0E0*1*I_KINETIC_Fy2z_D2z_a;
  abcd[59] = 4.0E0*I_KINETIC_H2x3z_D2z_aa-2.0E0*1*I_KINETIC_F3z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_D_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_D_aa
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_P_D
   ************************************************************/
  abcd[60] = 4.0E0*I_KINETIC_H4xy_D2x_aa-2.0E0*3*I_KINETIC_F2xy_D2x_a;
  abcd[61] = 4.0E0*I_KINETIC_H3x2y_D2x_aa-2.0E0*1*I_KINETIC_F3x_D2x_a-2.0E0*2*I_KINETIC_Fx2y_D2x_a+2*1*I_KINETIC_Px_D2x;
  abcd[62] = 4.0E0*I_KINETIC_H3xyz_D2x_aa-2.0E0*2*I_KINETIC_Fxyz_D2x_a;
  abcd[63] = 4.0E0*I_KINETIC_H2x3y_D2x_aa-2.0E0*2*I_KINETIC_F2xy_D2x_a-2.0E0*1*I_KINETIC_F3y_D2x_a+2*I_KINETIC_Py_D2x;
  abcd[64] = 4.0E0*I_KINETIC_H2x2yz_D2x_aa-2.0E0*1*I_KINETIC_F2xz_D2x_a-2.0E0*1*I_KINETIC_F2yz_D2x_a+1*I_KINETIC_Pz_D2x;
  abcd[65] = 4.0E0*I_KINETIC_H2xy2z_D2x_aa-2.0E0*1*I_KINETIC_Fy2z_D2x_a;
  abcd[66] = 4.0E0*I_KINETIC_Hx4y_D2x_aa-2.0E0*3*I_KINETIC_Fx2y_D2x_a;
  abcd[67] = 4.0E0*I_KINETIC_Hx3yz_D2x_aa-2.0E0*2*I_KINETIC_Fxyz_D2x_a;
  abcd[68] = 4.0E0*I_KINETIC_Hx2y2z_D2x_aa-2.0E0*1*I_KINETIC_Fx2z_D2x_a;
  abcd[69] = 4.0E0*I_KINETIC_Hxy3z_D2x_aa;
  abcd[70] = 4.0E0*I_KINETIC_H4xy_Dxy_aa-2.0E0*3*I_KINETIC_F2xy_Dxy_a;
  abcd[71] = 4.0E0*I_KINETIC_H3x2y_Dxy_aa-2.0E0*1*I_KINETIC_F3x_Dxy_a-2.0E0*2*I_KINETIC_Fx2y_Dxy_a+2*1*I_KINETIC_Px_Dxy;
  abcd[72] = 4.0E0*I_KINETIC_H3xyz_Dxy_aa-2.0E0*2*I_KINETIC_Fxyz_Dxy_a;
  abcd[73] = 4.0E0*I_KINETIC_H2x3y_Dxy_aa-2.0E0*2*I_KINETIC_F2xy_Dxy_a-2.0E0*1*I_KINETIC_F3y_Dxy_a+2*I_KINETIC_Py_Dxy;
  abcd[74] = 4.0E0*I_KINETIC_H2x2yz_Dxy_aa-2.0E0*1*I_KINETIC_F2xz_Dxy_a-2.0E0*1*I_KINETIC_F2yz_Dxy_a+1*I_KINETIC_Pz_Dxy;
  abcd[75] = 4.0E0*I_KINETIC_H2xy2z_Dxy_aa-2.0E0*1*I_KINETIC_Fy2z_Dxy_a;
  abcd[76] = 4.0E0*I_KINETIC_Hx4y_Dxy_aa-2.0E0*3*I_KINETIC_Fx2y_Dxy_a;
  abcd[77] = 4.0E0*I_KINETIC_Hx3yz_Dxy_aa-2.0E0*2*I_KINETIC_Fxyz_Dxy_a;
  abcd[78] = 4.0E0*I_KINETIC_Hx2y2z_Dxy_aa-2.0E0*1*I_KINETIC_Fx2z_Dxy_a;
  abcd[79] = 4.0E0*I_KINETIC_Hxy3z_Dxy_aa;
  abcd[80] = 4.0E0*I_KINETIC_H4xy_Dxz_aa-2.0E0*3*I_KINETIC_F2xy_Dxz_a;
  abcd[81] = 4.0E0*I_KINETIC_H3x2y_Dxz_aa-2.0E0*1*I_KINETIC_F3x_Dxz_a-2.0E0*2*I_KINETIC_Fx2y_Dxz_a+2*1*I_KINETIC_Px_Dxz;
  abcd[82] = 4.0E0*I_KINETIC_H3xyz_Dxz_aa-2.0E0*2*I_KINETIC_Fxyz_Dxz_a;
  abcd[83] = 4.0E0*I_KINETIC_H2x3y_Dxz_aa-2.0E0*2*I_KINETIC_F2xy_Dxz_a-2.0E0*1*I_KINETIC_F3y_Dxz_a+2*I_KINETIC_Py_Dxz;
  abcd[84] = 4.0E0*I_KINETIC_H2x2yz_Dxz_aa-2.0E0*1*I_KINETIC_F2xz_Dxz_a-2.0E0*1*I_KINETIC_F2yz_Dxz_a+1*I_KINETIC_Pz_Dxz;
  abcd[85] = 4.0E0*I_KINETIC_H2xy2z_Dxz_aa-2.0E0*1*I_KINETIC_Fy2z_Dxz_a;
  abcd[86] = 4.0E0*I_KINETIC_Hx4y_Dxz_aa-2.0E0*3*I_KINETIC_Fx2y_Dxz_a;
  abcd[87] = 4.0E0*I_KINETIC_Hx3yz_Dxz_aa-2.0E0*2*I_KINETIC_Fxyz_Dxz_a;
  abcd[88] = 4.0E0*I_KINETIC_Hx2y2z_Dxz_aa-2.0E0*1*I_KINETIC_Fx2z_Dxz_a;
  abcd[89] = 4.0E0*I_KINETIC_Hxy3z_Dxz_aa;
  abcd[90] = 4.0E0*I_KINETIC_H4xy_D2y_aa-2.0E0*3*I_KINETIC_F2xy_D2y_a;
  abcd[91] = 4.0E0*I_KINETIC_H3x2y_D2y_aa-2.0E0*1*I_KINETIC_F3x_D2y_a-2.0E0*2*I_KINETIC_Fx2y_D2y_a+2*1*I_KINETIC_Px_D2y;
  abcd[92] = 4.0E0*I_KINETIC_H3xyz_D2y_aa-2.0E0*2*I_KINETIC_Fxyz_D2y_a;
  abcd[93] = 4.0E0*I_KINETIC_H2x3y_D2y_aa-2.0E0*2*I_KINETIC_F2xy_D2y_a-2.0E0*1*I_KINETIC_F3y_D2y_a+2*I_KINETIC_Py_D2y;
  abcd[94] = 4.0E0*I_KINETIC_H2x2yz_D2y_aa-2.0E0*1*I_KINETIC_F2xz_D2y_a-2.0E0*1*I_KINETIC_F2yz_D2y_a+1*I_KINETIC_Pz_D2y;
  abcd[95] = 4.0E0*I_KINETIC_H2xy2z_D2y_aa-2.0E0*1*I_KINETIC_Fy2z_D2y_a;
  abcd[96] = 4.0E0*I_KINETIC_Hx4y_D2y_aa-2.0E0*3*I_KINETIC_Fx2y_D2y_a;
  abcd[97] = 4.0E0*I_KINETIC_Hx3yz_D2y_aa-2.0E0*2*I_KINETIC_Fxyz_D2y_a;
  abcd[98] = 4.0E0*I_KINETIC_Hx2y2z_D2y_aa-2.0E0*1*I_KINETIC_Fx2z_D2y_a;
  abcd[99] = 4.0E0*I_KINETIC_Hxy3z_D2y_aa;
  abcd[100] = 4.0E0*I_KINETIC_H4xy_Dyz_aa-2.0E0*3*I_KINETIC_F2xy_Dyz_a;
  abcd[101] = 4.0E0*I_KINETIC_H3x2y_Dyz_aa-2.0E0*1*I_KINETIC_F3x_Dyz_a-2.0E0*2*I_KINETIC_Fx2y_Dyz_a+2*1*I_KINETIC_Px_Dyz;
  abcd[102] = 4.0E0*I_KINETIC_H3xyz_Dyz_aa-2.0E0*2*I_KINETIC_Fxyz_Dyz_a;
  abcd[103] = 4.0E0*I_KINETIC_H2x3y_Dyz_aa-2.0E0*2*I_KINETIC_F2xy_Dyz_a-2.0E0*1*I_KINETIC_F3y_Dyz_a+2*I_KINETIC_Py_Dyz;
  abcd[104] = 4.0E0*I_KINETIC_H2x2yz_Dyz_aa-2.0E0*1*I_KINETIC_F2xz_Dyz_a-2.0E0*1*I_KINETIC_F2yz_Dyz_a+1*I_KINETIC_Pz_Dyz;
  abcd[105] = 4.0E0*I_KINETIC_H2xy2z_Dyz_aa-2.0E0*1*I_KINETIC_Fy2z_Dyz_a;
  abcd[106] = 4.0E0*I_KINETIC_Hx4y_Dyz_aa-2.0E0*3*I_KINETIC_Fx2y_Dyz_a;
  abcd[107] = 4.0E0*I_KINETIC_Hx3yz_Dyz_aa-2.0E0*2*I_KINETIC_Fxyz_Dyz_a;
  abcd[108] = 4.0E0*I_KINETIC_Hx2y2z_Dyz_aa-2.0E0*1*I_KINETIC_Fx2z_Dyz_a;
  abcd[109] = 4.0E0*I_KINETIC_Hxy3z_Dyz_aa;
  abcd[110] = 4.0E0*I_KINETIC_H4xy_D2z_aa-2.0E0*3*I_KINETIC_F2xy_D2z_a;
  abcd[111] = 4.0E0*I_KINETIC_H3x2y_D2z_aa-2.0E0*1*I_KINETIC_F3x_D2z_a-2.0E0*2*I_KINETIC_Fx2y_D2z_a+2*1*I_KINETIC_Px_D2z;
  abcd[112] = 4.0E0*I_KINETIC_H3xyz_D2z_aa-2.0E0*2*I_KINETIC_Fxyz_D2z_a;
  abcd[113] = 4.0E0*I_KINETIC_H2x3y_D2z_aa-2.0E0*2*I_KINETIC_F2xy_D2z_a-2.0E0*1*I_KINETIC_F3y_D2z_a+2*I_KINETIC_Py_D2z;
  abcd[114] = 4.0E0*I_KINETIC_H2x2yz_D2z_aa-2.0E0*1*I_KINETIC_F2xz_D2z_a-2.0E0*1*I_KINETIC_F2yz_D2z_a+1*I_KINETIC_Pz_D2z;
  abcd[115] = 4.0E0*I_KINETIC_H2xy2z_D2z_aa-2.0E0*1*I_KINETIC_Fy2z_D2z_a;
  abcd[116] = 4.0E0*I_KINETIC_Hx4y_D2z_aa-2.0E0*3*I_KINETIC_Fx2y_D2z_a;
  abcd[117] = 4.0E0*I_KINETIC_Hx3yz_D2z_aa-2.0E0*2*I_KINETIC_Fxyz_D2z_a;
  abcd[118] = 4.0E0*I_KINETIC_Hx2y2z_D2z_aa-2.0E0*1*I_KINETIC_Fx2z_D2z_a;
  abcd[119] = 4.0E0*I_KINETIC_Hxy3z_D2z_aa;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_D_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_D_aa
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_P_D
   ************************************************************/
  abcd[120] = 4.0E0*I_KINETIC_H4xz_D2x_aa-2.0E0*3*I_KINETIC_F2xz_D2x_a;
  abcd[121] = 4.0E0*I_KINETIC_H3xyz_D2x_aa-2.0E0*2*I_KINETIC_Fxyz_D2x_a;
  abcd[122] = 4.0E0*I_KINETIC_H3x2z_D2x_aa-2.0E0*1*I_KINETIC_F3x_D2x_a-2.0E0*2*I_KINETIC_Fx2z_D2x_a+2*1*I_KINETIC_Px_D2x;
  abcd[123] = 4.0E0*I_KINETIC_H2x2yz_D2x_aa-2.0E0*1*I_KINETIC_F2yz_D2x_a;
  abcd[124] = 4.0E0*I_KINETIC_H2xy2z_D2x_aa-2.0E0*1*I_KINETIC_F2xy_D2x_a-2.0E0*1*I_KINETIC_Fy2z_D2x_a+1*I_KINETIC_Py_D2x;
  abcd[125] = 4.0E0*I_KINETIC_H2x3z_D2x_aa-2.0E0*2*I_KINETIC_F2xz_D2x_a-2.0E0*1*I_KINETIC_F3z_D2x_a+2*I_KINETIC_Pz_D2x;
  abcd[126] = 4.0E0*I_KINETIC_Hx3yz_D2x_aa;
  abcd[127] = 4.0E0*I_KINETIC_Hx2y2z_D2x_aa-2.0E0*1*I_KINETIC_Fx2y_D2x_a;
  abcd[128] = 4.0E0*I_KINETIC_Hxy3z_D2x_aa-2.0E0*2*I_KINETIC_Fxyz_D2x_a;
  abcd[129] = 4.0E0*I_KINETIC_Hx4z_D2x_aa-2.0E0*3*I_KINETIC_Fx2z_D2x_a;
  abcd[130] = 4.0E0*I_KINETIC_H4xz_Dxy_aa-2.0E0*3*I_KINETIC_F2xz_Dxy_a;
  abcd[131] = 4.0E0*I_KINETIC_H3xyz_Dxy_aa-2.0E0*2*I_KINETIC_Fxyz_Dxy_a;
  abcd[132] = 4.0E0*I_KINETIC_H3x2z_Dxy_aa-2.0E0*1*I_KINETIC_F3x_Dxy_a-2.0E0*2*I_KINETIC_Fx2z_Dxy_a+2*1*I_KINETIC_Px_Dxy;
  abcd[133] = 4.0E0*I_KINETIC_H2x2yz_Dxy_aa-2.0E0*1*I_KINETIC_F2yz_Dxy_a;
  abcd[134] = 4.0E0*I_KINETIC_H2xy2z_Dxy_aa-2.0E0*1*I_KINETIC_F2xy_Dxy_a-2.0E0*1*I_KINETIC_Fy2z_Dxy_a+1*I_KINETIC_Py_Dxy;
  abcd[135] = 4.0E0*I_KINETIC_H2x3z_Dxy_aa-2.0E0*2*I_KINETIC_F2xz_Dxy_a-2.0E0*1*I_KINETIC_F3z_Dxy_a+2*I_KINETIC_Pz_Dxy;
  abcd[136] = 4.0E0*I_KINETIC_Hx3yz_Dxy_aa;
  abcd[137] = 4.0E0*I_KINETIC_Hx2y2z_Dxy_aa-2.0E0*1*I_KINETIC_Fx2y_Dxy_a;
  abcd[138] = 4.0E0*I_KINETIC_Hxy3z_Dxy_aa-2.0E0*2*I_KINETIC_Fxyz_Dxy_a;
  abcd[139] = 4.0E0*I_KINETIC_Hx4z_Dxy_aa-2.0E0*3*I_KINETIC_Fx2z_Dxy_a;
  abcd[140] = 4.0E0*I_KINETIC_H4xz_Dxz_aa-2.0E0*3*I_KINETIC_F2xz_Dxz_a;
  abcd[141] = 4.0E0*I_KINETIC_H3xyz_Dxz_aa-2.0E0*2*I_KINETIC_Fxyz_Dxz_a;
  abcd[142] = 4.0E0*I_KINETIC_H3x2z_Dxz_aa-2.0E0*1*I_KINETIC_F3x_Dxz_a-2.0E0*2*I_KINETIC_Fx2z_Dxz_a+2*1*I_KINETIC_Px_Dxz;
  abcd[143] = 4.0E0*I_KINETIC_H2x2yz_Dxz_aa-2.0E0*1*I_KINETIC_F2yz_Dxz_a;
  abcd[144] = 4.0E0*I_KINETIC_H2xy2z_Dxz_aa-2.0E0*1*I_KINETIC_F2xy_Dxz_a-2.0E0*1*I_KINETIC_Fy2z_Dxz_a+1*I_KINETIC_Py_Dxz;
  abcd[145] = 4.0E0*I_KINETIC_H2x3z_Dxz_aa-2.0E0*2*I_KINETIC_F2xz_Dxz_a-2.0E0*1*I_KINETIC_F3z_Dxz_a+2*I_KINETIC_Pz_Dxz;
  abcd[146] = 4.0E0*I_KINETIC_Hx3yz_Dxz_aa;
  abcd[147] = 4.0E0*I_KINETIC_Hx2y2z_Dxz_aa-2.0E0*1*I_KINETIC_Fx2y_Dxz_a;
  abcd[148] = 4.0E0*I_KINETIC_Hxy3z_Dxz_aa-2.0E0*2*I_KINETIC_Fxyz_Dxz_a;
  abcd[149] = 4.0E0*I_KINETIC_Hx4z_Dxz_aa-2.0E0*3*I_KINETIC_Fx2z_Dxz_a;
  abcd[150] = 4.0E0*I_KINETIC_H4xz_D2y_aa-2.0E0*3*I_KINETIC_F2xz_D2y_a;
  abcd[151] = 4.0E0*I_KINETIC_H3xyz_D2y_aa-2.0E0*2*I_KINETIC_Fxyz_D2y_a;
  abcd[152] = 4.0E0*I_KINETIC_H3x2z_D2y_aa-2.0E0*1*I_KINETIC_F3x_D2y_a-2.0E0*2*I_KINETIC_Fx2z_D2y_a+2*1*I_KINETIC_Px_D2y;
  abcd[153] = 4.0E0*I_KINETIC_H2x2yz_D2y_aa-2.0E0*1*I_KINETIC_F2yz_D2y_a;
  abcd[154] = 4.0E0*I_KINETIC_H2xy2z_D2y_aa-2.0E0*1*I_KINETIC_F2xy_D2y_a-2.0E0*1*I_KINETIC_Fy2z_D2y_a+1*I_KINETIC_Py_D2y;
  abcd[155] = 4.0E0*I_KINETIC_H2x3z_D2y_aa-2.0E0*2*I_KINETIC_F2xz_D2y_a-2.0E0*1*I_KINETIC_F3z_D2y_a+2*I_KINETIC_Pz_D2y;
  abcd[156] = 4.0E0*I_KINETIC_Hx3yz_D2y_aa;
  abcd[157] = 4.0E0*I_KINETIC_Hx2y2z_D2y_aa-2.0E0*1*I_KINETIC_Fx2y_D2y_a;
  abcd[158] = 4.0E0*I_KINETIC_Hxy3z_D2y_aa-2.0E0*2*I_KINETIC_Fxyz_D2y_a;
  abcd[159] = 4.0E0*I_KINETIC_Hx4z_D2y_aa-2.0E0*3*I_KINETIC_Fx2z_D2y_a;
  abcd[160] = 4.0E0*I_KINETIC_H4xz_Dyz_aa-2.0E0*3*I_KINETIC_F2xz_Dyz_a;
  abcd[161] = 4.0E0*I_KINETIC_H3xyz_Dyz_aa-2.0E0*2*I_KINETIC_Fxyz_Dyz_a;
  abcd[162] = 4.0E0*I_KINETIC_H3x2z_Dyz_aa-2.0E0*1*I_KINETIC_F3x_Dyz_a-2.0E0*2*I_KINETIC_Fx2z_Dyz_a+2*1*I_KINETIC_Px_Dyz;
  abcd[163] = 4.0E0*I_KINETIC_H2x2yz_Dyz_aa-2.0E0*1*I_KINETIC_F2yz_Dyz_a;
  abcd[164] = 4.0E0*I_KINETIC_H2xy2z_Dyz_aa-2.0E0*1*I_KINETIC_F2xy_Dyz_a-2.0E0*1*I_KINETIC_Fy2z_Dyz_a+1*I_KINETIC_Py_Dyz;
  abcd[165] = 4.0E0*I_KINETIC_H2x3z_Dyz_aa-2.0E0*2*I_KINETIC_F2xz_Dyz_a-2.0E0*1*I_KINETIC_F3z_Dyz_a+2*I_KINETIC_Pz_Dyz;
  abcd[166] = 4.0E0*I_KINETIC_Hx3yz_Dyz_aa;
  abcd[167] = 4.0E0*I_KINETIC_Hx2y2z_Dyz_aa-2.0E0*1*I_KINETIC_Fx2y_Dyz_a;
  abcd[168] = 4.0E0*I_KINETIC_Hxy3z_Dyz_aa-2.0E0*2*I_KINETIC_Fxyz_Dyz_a;
  abcd[169] = 4.0E0*I_KINETIC_Hx4z_Dyz_aa-2.0E0*3*I_KINETIC_Fx2z_Dyz_a;
  abcd[170] = 4.0E0*I_KINETIC_H4xz_D2z_aa-2.0E0*3*I_KINETIC_F2xz_D2z_a;
  abcd[171] = 4.0E0*I_KINETIC_H3xyz_D2z_aa-2.0E0*2*I_KINETIC_Fxyz_D2z_a;
  abcd[172] = 4.0E0*I_KINETIC_H3x2z_D2z_aa-2.0E0*1*I_KINETIC_F3x_D2z_a-2.0E0*2*I_KINETIC_Fx2z_D2z_a+2*1*I_KINETIC_Px_D2z;
  abcd[173] = 4.0E0*I_KINETIC_H2x2yz_D2z_aa-2.0E0*1*I_KINETIC_F2yz_D2z_a;
  abcd[174] = 4.0E0*I_KINETIC_H2xy2z_D2z_aa-2.0E0*1*I_KINETIC_F2xy_D2z_a-2.0E0*1*I_KINETIC_Fy2z_D2z_a+1*I_KINETIC_Py_D2z;
  abcd[175] = 4.0E0*I_KINETIC_H2x3z_D2z_aa-2.0E0*2*I_KINETIC_F2xz_D2z_a-2.0E0*1*I_KINETIC_F3z_D2z_a+2*I_KINETIC_Pz_D2z;
  abcd[176] = 4.0E0*I_KINETIC_Hx3yz_D2z_aa;
  abcd[177] = 4.0E0*I_KINETIC_Hx2y2z_D2z_aa-2.0E0*1*I_KINETIC_Fx2y_D2z_a;
  abcd[178] = 4.0E0*I_KINETIC_Hxy3z_D2z_aa-2.0E0*2*I_KINETIC_Fxyz_D2z_a;
  abcd[179] = 4.0E0*I_KINETIC_Hx4z_D2z_aa-2.0E0*3*I_KINETIC_Fx2z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_D_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_D_aa
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_P_D
   ************************************************************/
  abcd[180] = 4.0E0*I_KINETIC_H3x2y_D2x_aa-2.0E0*1*I_KINETIC_F3x_D2x_a;
  abcd[181] = 4.0E0*I_KINETIC_H2x3y_D2x_aa-2.0E0*1*I_KINETIC_F2xy_D2x_a-2.0E0*2*I_KINETIC_F2xy_D2x_a;
  abcd[182] = 4.0E0*I_KINETIC_H2x2yz_D2x_aa-2.0E0*1*I_KINETIC_F2xz_D2x_a;
  abcd[183] = 4.0E0*I_KINETIC_Hx4y_D2x_aa-2.0E0*2*I_KINETIC_Fx2y_D2x_a-2.0E0*3*I_KINETIC_Fx2y_D2x_a+2*1*I_KINETIC_Px_D2x;
  abcd[184] = 4.0E0*I_KINETIC_Hx3yz_D2x_aa-2.0E0*1*I_KINETIC_Fxyz_D2x_a-2.0E0*2*I_KINETIC_Fxyz_D2x_a;
  abcd[185] = 4.0E0*I_KINETIC_Hx2y2z_D2x_aa-2.0E0*1*I_KINETIC_Fx2z_D2x_a;
  abcd[186] = 4.0E0*I_KINETIC_H5y_D2x_aa-2.0E0*3*I_KINETIC_F3y_D2x_a-2.0E0*4*I_KINETIC_F3y_D2x_a+3*2*I_KINETIC_Py_D2x;
  abcd[187] = 4.0E0*I_KINETIC_H4yz_D2x_aa-2.0E0*2*I_KINETIC_F2yz_D2x_a-2.0E0*3*I_KINETIC_F2yz_D2x_a+2*1*I_KINETIC_Pz_D2x;
  abcd[188] = 4.0E0*I_KINETIC_H3y2z_D2x_aa-2.0E0*1*I_KINETIC_Fy2z_D2x_a-2.0E0*2*I_KINETIC_Fy2z_D2x_a;
  abcd[189] = 4.0E0*I_KINETIC_H2y3z_D2x_aa-2.0E0*1*I_KINETIC_F3z_D2x_a;
  abcd[190] = 4.0E0*I_KINETIC_H3x2y_Dxy_aa-2.0E0*1*I_KINETIC_F3x_Dxy_a;
  abcd[191] = 4.0E0*I_KINETIC_H2x3y_Dxy_aa-2.0E0*1*I_KINETIC_F2xy_Dxy_a-2.0E0*2*I_KINETIC_F2xy_Dxy_a;
  abcd[192] = 4.0E0*I_KINETIC_H2x2yz_Dxy_aa-2.0E0*1*I_KINETIC_F2xz_Dxy_a;
  abcd[193] = 4.0E0*I_KINETIC_Hx4y_Dxy_aa-2.0E0*2*I_KINETIC_Fx2y_Dxy_a-2.0E0*3*I_KINETIC_Fx2y_Dxy_a+2*1*I_KINETIC_Px_Dxy;
  abcd[194] = 4.0E0*I_KINETIC_Hx3yz_Dxy_aa-2.0E0*1*I_KINETIC_Fxyz_Dxy_a-2.0E0*2*I_KINETIC_Fxyz_Dxy_a;
  abcd[195] = 4.0E0*I_KINETIC_Hx2y2z_Dxy_aa-2.0E0*1*I_KINETIC_Fx2z_Dxy_a;
  abcd[196] = 4.0E0*I_KINETIC_H5y_Dxy_aa-2.0E0*3*I_KINETIC_F3y_Dxy_a-2.0E0*4*I_KINETIC_F3y_Dxy_a+3*2*I_KINETIC_Py_Dxy;
  abcd[197] = 4.0E0*I_KINETIC_H4yz_Dxy_aa-2.0E0*2*I_KINETIC_F2yz_Dxy_a-2.0E0*3*I_KINETIC_F2yz_Dxy_a+2*1*I_KINETIC_Pz_Dxy;
  abcd[198] = 4.0E0*I_KINETIC_H3y2z_Dxy_aa-2.0E0*1*I_KINETIC_Fy2z_Dxy_a-2.0E0*2*I_KINETIC_Fy2z_Dxy_a;
  abcd[199] = 4.0E0*I_KINETIC_H2y3z_Dxy_aa-2.0E0*1*I_KINETIC_F3z_Dxy_a;
  abcd[200] = 4.0E0*I_KINETIC_H3x2y_Dxz_aa-2.0E0*1*I_KINETIC_F3x_Dxz_a;
  abcd[201] = 4.0E0*I_KINETIC_H2x3y_Dxz_aa-2.0E0*1*I_KINETIC_F2xy_Dxz_a-2.0E0*2*I_KINETIC_F2xy_Dxz_a;
  abcd[202] = 4.0E0*I_KINETIC_H2x2yz_Dxz_aa-2.0E0*1*I_KINETIC_F2xz_Dxz_a;
  abcd[203] = 4.0E0*I_KINETIC_Hx4y_Dxz_aa-2.0E0*2*I_KINETIC_Fx2y_Dxz_a-2.0E0*3*I_KINETIC_Fx2y_Dxz_a+2*1*I_KINETIC_Px_Dxz;
  abcd[204] = 4.0E0*I_KINETIC_Hx3yz_Dxz_aa-2.0E0*1*I_KINETIC_Fxyz_Dxz_a-2.0E0*2*I_KINETIC_Fxyz_Dxz_a;
  abcd[205] = 4.0E0*I_KINETIC_Hx2y2z_Dxz_aa-2.0E0*1*I_KINETIC_Fx2z_Dxz_a;
  abcd[206] = 4.0E0*I_KINETIC_H5y_Dxz_aa-2.0E0*3*I_KINETIC_F3y_Dxz_a-2.0E0*4*I_KINETIC_F3y_Dxz_a+3*2*I_KINETIC_Py_Dxz;
  abcd[207] = 4.0E0*I_KINETIC_H4yz_Dxz_aa-2.0E0*2*I_KINETIC_F2yz_Dxz_a-2.0E0*3*I_KINETIC_F2yz_Dxz_a+2*1*I_KINETIC_Pz_Dxz;
  abcd[208] = 4.0E0*I_KINETIC_H3y2z_Dxz_aa-2.0E0*1*I_KINETIC_Fy2z_Dxz_a-2.0E0*2*I_KINETIC_Fy2z_Dxz_a;
  abcd[209] = 4.0E0*I_KINETIC_H2y3z_Dxz_aa-2.0E0*1*I_KINETIC_F3z_Dxz_a;
  abcd[210] = 4.0E0*I_KINETIC_H3x2y_D2y_aa-2.0E0*1*I_KINETIC_F3x_D2y_a;
  abcd[211] = 4.0E0*I_KINETIC_H2x3y_D2y_aa-2.0E0*1*I_KINETIC_F2xy_D2y_a-2.0E0*2*I_KINETIC_F2xy_D2y_a;
  abcd[212] = 4.0E0*I_KINETIC_H2x2yz_D2y_aa-2.0E0*1*I_KINETIC_F2xz_D2y_a;
  abcd[213] = 4.0E0*I_KINETIC_Hx4y_D2y_aa-2.0E0*2*I_KINETIC_Fx2y_D2y_a-2.0E0*3*I_KINETIC_Fx2y_D2y_a+2*1*I_KINETIC_Px_D2y;
  abcd[214] = 4.0E0*I_KINETIC_Hx3yz_D2y_aa-2.0E0*1*I_KINETIC_Fxyz_D2y_a-2.0E0*2*I_KINETIC_Fxyz_D2y_a;
  abcd[215] = 4.0E0*I_KINETIC_Hx2y2z_D2y_aa-2.0E0*1*I_KINETIC_Fx2z_D2y_a;
  abcd[216] = 4.0E0*I_KINETIC_H5y_D2y_aa-2.0E0*3*I_KINETIC_F3y_D2y_a-2.0E0*4*I_KINETIC_F3y_D2y_a+3*2*I_KINETIC_Py_D2y;
  abcd[217] = 4.0E0*I_KINETIC_H4yz_D2y_aa-2.0E0*2*I_KINETIC_F2yz_D2y_a-2.0E0*3*I_KINETIC_F2yz_D2y_a+2*1*I_KINETIC_Pz_D2y;
  abcd[218] = 4.0E0*I_KINETIC_H3y2z_D2y_aa-2.0E0*1*I_KINETIC_Fy2z_D2y_a-2.0E0*2*I_KINETIC_Fy2z_D2y_a;
  abcd[219] = 4.0E0*I_KINETIC_H2y3z_D2y_aa-2.0E0*1*I_KINETIC_F3z_D2y_a;
  abcd[220] = 4.0E0*I_KINETIC_H3x2y_Dyz_aa-2.0E0*1*I_KINETIC_F3x_Dyz_a;
  abcd[221] = 4.0E0*I_KINETIC_H2x3y_Dyz_aa-2.0E0*1*I_KINETIC_F2xy_Dyz_a-2.0E0*2*I_KINETIC_F2xy_Dyz_a;
  abcd[222] = 4.0E0*I_KINETIC_H2x2yz_Dyz_aa-2.0E0*1*I_KINETIC_F2xz_Dyz_a;
  abcd[223] = 4.0E0*I_KINETIC_Hx4y_Dyz_aa-2.0E0*2*I_KINETIC_Fx2y_Dyz_a-2.0E0*3*I_KINETIC_Fx2y_Dyz_a+2*1*I_KINETIC_Px_Dyz;
  abcd[224] = 4.0E0*I_KINETIC_Hx3yz_Dyz_aa-2.0E0*1*I_KINETIC_Fxyz_Dyz_a-2.0E0*2*I_KINETIC_Fxyz_Dyz_a;
  abcd[225] = 4.0E0*I_KINETIC_Hx2y2z_Dyz_aa-2.0E0*1*I_KINETIC_Fx2z_Dyz_a;
  abcd[226] = 4.0E0*I_KINETIC_H5y_Dyz_aa-2.0E0*3*I_KINETIC_F3y_Dyz_a-2.0E0*4*I_KINETIC_F3y_Dyz_a+3*2*I_KINETIC_Py_Dyz;
  abcd[227] = 4.0E0*I_KINETIC_H4yz_Dyz_aa-2.0E0*2*I_KINETIC_F2yz_Dyz_a-2.0E0*3*I_KINETIC_F2yz_Dyz_a+2*1*I_KINETIC_Pz_Dyz;
  abcd[228] = 4.0E0*I_KINETIC_H3y2z_Dyz_aa-2.0E0*1*I_KINETIC_Fy2z_Dyz_a-2.0E0*2*I_KINETIC_Fy2z_Dyz_a;
  abcd[229] = 4.0E0*I_KINETIC_H2y3z_Dyz_aa-2.0E0*1*I_KINETIC_F3z_Dyz_a;
  abcd[230] = 4.0E0*I_KINETIC_H3x2y_D2z_aa-2.0E0*1*I_KINETIC_F3x_D2z_a;
  abcd[231] = 4.0E0*I_KINETIC_H2x3y_D2z_aa-2.0E0*1*I_KINETIC_F2xy_D2z_a-2.0E0*2*I_KINETIC_F2xy_D2z_a;
  abcd[232] = 4.0E0*I_KINETIC_H2x2yz_D2z_aa-2.0E0*1*I_KINETIC_F2xz_D2z_a;
  abcd[233] = 4.0E0*I_KINETIC_Hx4y_D2z_aa-2.0E0*2*I_KINETIC_Fx2y_D2z_a-2.0E0*3*I_KINETIC_Fx2y_D2z_a+2*1*I_KINETIC_Px_D2z;
  abcd[234] = 4.0E0*I_KINETIC_Hx3yz_D2z_aa-2.0E0*1*I_KINETIC_Fxyz_D2z_a-2.0E0*2*I_KINETIC_Fxyz_D2z_a;
  abcd[235] = 4.0E0*I_KINETIC_Hx2y2z_D2z_aa-2.0E0*1*I_KINETIC_Fx2z_D2z_a;
  abcd[236] = 4.0E0*I_KINETIC_H5y_D2z_aa-2.0E0*3*I_KINETIC_F3y_D2z_a-2.0E0*4*I_KINETIC_F3y_D2z_a+3*2*I_KINETIC_Py_D2z;
  abcd[237] = 4.0E0*I_KINETIC_H4yz_D2z_aa-2.0E0*2*I_KINETIC_F2yz_D2z_a-2.0E0*3*I_KINETIC_F2yz_D2z_a+2*1*I_KINETIC_Pz_D2z;
  abcd[238] = 4.0E0*I_KINETIC_H3y2z_D2z_aa-2.0E0*1*I_KINETIC_Fy2z_D2z_a-2.0E0*2*I_KINETIC_Fy2z_D2z_a;
  abcd[239] = 4.0E0*I_KINETIC_H2y3z_D2z_aa-2.0E0*1*I_KINETIC_F3z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_D_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_D_aa
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_P_D
   ************************************************************/
  abcd[240] = 4.0E0*I_KINETIC_H3xyz_D2x_aa;
  abcd[241] = 4.0E0*I_KINETIC_H2x2yz_D2x_aa-2.0E0*1*I_KINETIC_F2xz_D2x_a;
  abcd[242] = 4.0E0*I_KINETIC_H2xy2z_D2x_aa-2.0E0*1*I_KINETIC_F2xy_D2x_a;
  abcd[243] = 4.0E0*I_KINETIC_Hx3yz_D2x_aa-2.0E0*2*I_KINETIC_Fxyz_D2x_a;
  abcd[244] = 4.0E0*I_KINETIC_Hx2y2z_D2x_aa-2.0E0*1*I_KINETIC_Fx2y_D2x_a-2.0E0*1*I_KINETIC_Fx2z_D2x_a+1*I_KINETIC_Px_D2x;
  abcd[245] = 4.0E0*I_KINETIC_Hxy3z_D2x_aa-2.0E0*2*I_KINETIC_Fxyz_D2x_a;
  abcd[246] = 4.0E0*I_KINETIC_H4yz_D2x_aa-2.0E0*3*I_KINETIC_F2yz_D2x_a;
  abcd[247] = 4.0E0*I_KINETIC_H3y2z_D2x_aa-2.0E0*1*I_KINETIC_F3y_D2x_a-2.0E0*2*I_KINETIC_Fy2z_D2x_a+2*1*I_KINETIC_Py_D2x;
  abcd[248] = 4.0E0*I_KINETIC_H2y3z_D2x_aa-2.0E0*2*I_KINETIC_F2yz_D2x_a-2.0E0*1*I_KINETIC_F3z_D2x_a+2*I_KINETIC_Pz_D2x;
  abcd[249] = 4.0E0*I_KINETIC_Hy4z_D2x_aa-2.0E0*3*I_KINETIC_Fy2z_D2x_a;
  abcd[250] = 4.0E0*I_KINETIC_H3xyz_Dxy_aa;
  abcd[251] = 4.0E0*I_KINETIC_H2x2yz_Dxy_aa-2.0E0*1*I_KINETIC_F2xz_Dxy_a;
  abcd[252] = 4.0E0*I_KINETIC_H2xy2z_Dxy_aa-2.0E0*1*I_KINETIC_F2xy_Dxy_a;
  abcd[253] = 4.0E0*I_KINETIC_Hx3yz_Dxy_aa-2.0E0*2*I_KINETIC_Fxyz_Dxy_a;
  abcd[254] = 4.0E0*I_KINETIC_Hx2y2z_Dxy_aa-2.0E0*1*I_KINETIC_Fx2y_Dxy_a-2.0E0*1*I_KINETIC_Fx2z_Dxy_a+1*I_KINETIC_Px_Dxy;
  abcd[255] = 4.0E0*I_KINETIC_Hxy3z_Dxy_aa-2.0E0*2*I_KINETIC_Fxyz_Dxy_a;
  abcd[256] = 4.0E0*I_KINETIC_H4yz_Dxy_aa-2.0E0*3*I_KINETIC_F2yz_Dxy_a;
  abcd[257] = 4.0E0*I_KINETIC_H3y2z_Dxy_aa-2.0E0*1*I_KINETIC_F3y_Dxy_a-2.0E0*2*I_KINETIC_Fy2z_Dxy_a+2*1*I_KINETIC_Py_Dxy;
  abcd[258] = 4.0E0*I_KINETIC_H2y3z_Dxy_aa-2.0E0*2*I_KINETIC_F2yz_Dxy_a-2.0E0*1*I_KINETIC_F3z_Dxy_a+2*I_KINETIC_Pz_Dxy;
  abcd[259] = 4.0E0*I_KINETIC_Hy4z_Dxy_aa-2.0E0*3*I_KINETIC_Fy2z_Dxy_a;
  abcd[260] = 4.0E0*I_KINETIC_H3xyz_Dxz_aa;
  abcd[261] = 4.0E0*I_KINETIC_H2x2yz_Dxz_aa-2.0E0*1*I_KINETIC_F2xz_Dxz_a;
  abcd[262] = 4.0E0*I_KINETIC_H2xy2z_Dxz_aa-2.0E0*1*I_KINETIC_F2xy_Dxz_a;
  abcd[263] = 4.0E0*I_KINETIC_Hx3yz_Dxz_aa-2.0E0*2*I_KINETIC_Fxyz_Dxz_a;
  abcd[264] = 4.0E0*I_KINETIC_Hx2y2z_Dxz_aa-2.0E0*1*I_KINETIC_Fx2y_Dxz_a-2.0E0*1*I_KINETIC_Fx2z_Dxz_a+1*I_KINETIC_Px_Dxz;
  abcd[265] = 4.0E0*I_KINETIC_Hxy3z_Dxz_aa-2.0E0*2*I_KINETIC_Fxyz_Dxz_a;
  abcd[266] = 4.0E0*I_KINETIC_H4yz_Dxz_aa-2.0E0*3*I_KINETIC_F2yz_Dxz_a;
  abcd[267] = 4.0E0*I_KINETIC_H3y2z_Dxz_aa-2.0E0*1*I_KINETIC_F3y_Dxz_a-2.0E0*2*I_KINETIC_Fy2z_Dxz_a+2*1*I_KINETIC_Py_Dxz;
  abcd[268] = 4.0E0*I_KINETIC_H2y3z_Dxz_aa-2.0E0*2*I_KINETIC_F2yz_Dxz_a-2.0E0*1*I_KINETIC_F3z_Dxz_a+2*I_KINETIC_Pz_Dxz;
  abcd[269] = 4.0E0*I_KINETIC_Hy4z_Dxz_aa-2.0E0*3*I_KINETIC_Fy2z_Dxz_a;
  abcd[270] = 4.0E0*I_KINETIC_H3xyz_D2y_aa;
  abcd[271] = 4.0E0*I_KINETIC_H2x2yz_D2y_aa-2.0E0*1*I_KINETIC_F2xz_D2y_a;
  abcd[272] = 4.0E0*I_KINETIC_H2xy2z_D2y_aa-2.0E0*1*I_KINETIC_F2xy_D2y_a;
  abcd[273] = 4.0E0*I_KINETIC_Hx3yz_D2y_aa-2.0E0*2*I_KINETIC_Fxyz_D2y_a;
  abcd[274] = 4.0E0*I_KINETIC_Hx2y2z_D2y_aa-2.0E0*1*I_KINETIC_Fx2y_D2y_a-2.0E0*1*I_KINETIC_Fx2z_D2y_a+1*I_KINETIC_Px_D2y;
  abcd[275] = 4.0E0*I_KINETIC_Hxy3z_D2y_aa-2.0E0*2*I_KINETIC_Fxyz_D2y_a;
  abcd[276] = 4.0E0*I_KINETIC_H4yz_D2y_aa-2.0E0*3*I_KINETIC_F2yz_D2y_a;
  abcd[277] = 4.0E0*I_KINETIC_H3y2z_D2y_aa-2.0E0*1*I_KINETIC_F3y_D2y_a-2.0E0*2*I_KINETIC_Fy2z_D2y_a+2*1*I_KINETIC_Py_D2y;
  abcd[278] = 4.0E0*I_KINETIC_H2y3z_D2y_aa-2.0E0*2*I_KINETIC_F2yz_D2y_a-2.0E0*1*I_KINETIC_F3z_D2y_a+2*I_KINETIC_Pz_D2y;
  abcd[279] = 4.0E0*I_KINETIC_Hy4z_D2y_aa-2.0E0*3*I_KINETIC_Fy2z_D2y_a;
  abcd[280] = 4.0E0*I_KINETIC_H3xyz_Dyz_aa;
  abcd[281] = 4.0E0*I_KINETIC_H2x2yz_Dyz_aa-2.0E0*1*I_KINETIC_F2xz_Dyz_a;
  abcd[282] = 4.0E0*I_KINETIC_H2xy2z_Dyz_aa-2.0E0*1*I_KINETIC_F2xy_Dyz_a;
  abcd[283] = 4.0E0*I_KINETIC_Hx3yz_Dyz_aa-2.0E0*2*I_KINETIC_Fxyz_Dyz_a;
  abcd[284] = 4.0E0*I_KINETIC_Hx2y2z_Dyz_aa-2.0E0*1*I_KINETIC_Fx2y_Dyz_a-2.0E0*1*I_KINETIC_Fx2z_Dyz_a+1*I_KINETIC_Px_Dyz;
  abcd[285] = 4.0E0*I_KINETIC_Hxy3z_Dyz_aa-2.0E0*2*I_KINETIC_Fxyz_Dyz_a;
  abcd[286] = 4.0E0*I_KINETIC_H4yz_Dyz_aa-2.0E0*3*I_KINETIC_F2yz_Dyz_a;
  abcd[287] = 4.0E0*I_KINETIC_H3y2z_Dyz_aa-2.0E0*1*I_KINETIC_F3y_Dyz_a-2.0E0*2*I_KINETIC_Fy2z_Dyz_a+2*1*I_KINETIC_Py_Dyz;
  abcd[288] = 4.0E0*I_KINETIC_H2y3z_Dyz_aa-2.0E0*2*I_KINETIC_F2yz_Dyz_a-2.0E0*1*I_KINETIC_F3z_Dyz_a+2*I_KINETIC_Pz_Dyz;
  abcd[289] = 4.0E0*I_KINETIC_Hy4z_Dyz_aa-2.0E0*3*I_KINETIC_Fy2z_Dyz_a;
  abcd[290] = 4.0E0*I_KINETIC_H3xyz_D2z_aa;
  abcd[291] = 4.0E0*I_KINETIC_H2x2yz_D2z_aa-2.0E0*1*I_KINETIC_F2xz_D2z_a;
  abcd[292] = 4.0E0*I_KINETIC_H2xy2z_D2z_aa-2.0E0*1*I_KINETIC_F2xy_D2z_a;
  abcd[293] = 4.0E0*I_KINETIC_Hx3yz_D2z_aa-2.0E0*2*I_KINETIC_Fxyz_D2z_a;
  abcd[294] = 4.0E0*I_KINETIC_Hx2y2z_D2z_aa-2.0E0*1*I_KINETIC_Fx2y_D2z_a-2.0E0*1*I_KINETIC_Fx2z_D2z_a+1*I_KINETIC_Px_D2z;
  abcd[295] = 4.0E0*I_KINETIC_Hxy3z_D2z_aa-2.0E0*2*I_KINETIC_Fxyz_D2z_a;
  abcd[296] = 4.0E0*I_KINETIC_H4yz_D2z_aa-2.0E0*3*I_KINETIC_F2yz_D2z_a;
  abcd[297] = 4.0E0*I_KINETIC_H3y2z_D2z_aa-2.0E0*1*I_KINETIC_F3y_D2z_a-2.0E0*2*I_KINETIC_Fy2z_D2z_a+2*1*I_KINETIC_Py_D2z;
  abcd[298] = 4.0E0*I_KINETIC_H2y3z_D2z_aa-2.0E0*2*I_KINETIC_F2yz_D2z_a-2.0E0*1*I_KINETIC_F3z_D2z_a+2*I_KINETIC_Pz_D2z;
  abcd[299] = 4.0E0*I_KINETIC_Hy4z_D2z_aa-2.0E0*3*I_KINETIC_Fy2z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_D_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_H_D_aa
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_F_D_a
   * RHS shell quartet name: SQ_KINETIC_P_D
   ************************************************************/
  abcd[300] = 4.0E0*I_KINETIC_H3x2z_D2x_aa-2.0E0*1*I_KINETIC_F3x_D2x_a;
  abcd[301] = 4.0E0*I_KINETIC_H2xy2z_D2x_aa-2.0E0*1*I_KINETIC_F2xy_D2x_a;
  abcd[302] = 4.0E0*I_KINETIC_H2x3z_D2x_aa-2.0E0*1*I_KINETIC_F2xz_D2x_a-2.0E0*2*I_KINETIC_F2xz_D2x_a;
  abcd[303] = 4.0E0*I_KINETIC_Hx2y2z_D2x_aa-2.0E0*1*I_KINETIC_Fx2y_D2x_a;
  abcd[304] = 4.0E0*I_KINETIC_Hxy3z_D2x_aa-2.0E0*1*I_KINETIC_Fxyz_D2x_a-2.0E0*2*I_KINETIC_Fxyz_D2x_a;
  abcd[305] = 4.0E0*I_KINETIC_Hx4z_D2x_aa-2.0E0*2*I_KINETIC_Fx2z_D2x_a-2.0E0*3*I_KINETIC_Fx2z_D2x_a+2*1*I_KINETIC_Px_D2x;
  abcd[306] = 4.0E0*I_KINETIC_H3y2z_D2x_aa-2.0E0*1*I_KINETIC_F3y_D2x_a;
  abcd[307] = 4.0E0*I_KINETIC_H2y3z_D2x_aa-2.0E0*1*I_KINETIC_F2yz_D2x_a-2.0E0*2*I_KINETIC_F2yz_D2x_a;
  abcd[308] = 4.0E0*I_KINETIC_Hy4z_D2x_aa-2.0E0*2*I_KINETIC_Fy2z_D2x_a-2.0E0*3*I_KINETIC_Fy2z_D2x_a+2*1*I_KINETIC_Py_D2x;
  abcd[309] = 4.0E0*I_KINETIC_H5z_D2x_aa-2.0E0*3*I_KINETIC_F3z_D2x_a-2.0E0*4*I_KINETIC_F3z_D2x_a+3*2*I_KINETIC_Pz_D2x;
  abcd[310] = 4.0E0*I_KINETIC_H3x2z_Dxy_aa-2.0E0*1*I_KINETIC_F3x_Dxy_a;
  abcd[311] = 4.0E0*I_KINETIC_H2xy2z_Dxy_aa-2.0E0*1*I_KINETIC_F2xy_Dxy_a;
  abcd[312] = 4.0E0*I_KINETIC_H2x3z_Dxy_aa-2.0E0*1*I_KINETIC_F2xz_Dxy_a-2.0E0*2*I_KINETIC_F2xz_Dxy_a;
  abcd[313] = 4.0E0*I_KINETIC_Hx2y2z_Dxy_aa-2.0E0*1*I_KINETIC_Fx2y_Dxy_a;
  abcd[314] = 4.0E0*I_KINETIC_Hxy3z_Dxy_aa-2.0E0*1*I_KINETIC_Fxyz_Dxy_a-2.0E0*2*I_KINETIC_Fxyz_Dxy_a;
  abcd[315] = 4.0E0*I_KINETIC_Hx4z_Dxy_aa-2.0E0*2*I_KINETIC_Fx2z_Dxy_a-2.0E0*3*I_KINETIC_Fx2z_Dxy_a+2*1*I_KINETIC_Px_Dxy;
  abcd[316] = 4.0E0*I_KINETIC_H3y2z_Dxy_aa-2.0E0*1*I_KINETIC_F3y_Dxy_a;
  abcd[317] = 4.0E0*I_KINETIC_H2y3z_Dxy_aa-2.0E0*1*I_KINETIC_F2yz_Dxy_a-2.0E0*2*I_KINETIC_F2yz_Dxy_a;
  abcd[318] = 4.0E0*I_KINETIC_Hy4z_Dxy_aa-2.0E0*2*I_KINETIC_Fy2z_Dxy_a-2.0E0*3*I_KINETIC_Fy2z_Dxy_a+2*1*I_KINETIC_Py_Dxy;
  abcd[319] = 4.0E0*I_KINETIC_H5z_Dxy_aa-2.0E0*3*I_KINETIC_F3z_Dxy_a-2.0E0*4*I_KINETIC_F3z_Dxy_a+3*2*I_KINETIC_Pz_Dxy;
  abcd[320] = 4.0E0*I_KINETIC_H3x2z_Dxz_aa-2.0E0*1*I_KINETIC_F3x_Dxz_a;
  abcd[321] = 4.0E0*I_KINETIC_H2xy2z_Dxz_aa-2.0E0*1*I_KINETIC_F2xy_Dxz_a;
  abcd[322] = 4.0E0*I_KINETIC_H2x3z_Dxz_aa-2.0E0*1*I_KINETIC_F2xz_Dxz_a-2.0E0*2*I_KINETIC_F2xz_Dxz_a;
  abcd[323] = 4.0E0*I_KINETIC_Hx2y2z_Dxz_aa-2.0E0*1*I_KINETIC_Fx2y_Dxz_a;
  abcd[324] = 4.0E0*I_KINETIC_Hxy3z_Dxz_aa-2.0E0*1*I_KINETIC_Fxyz_Dxz_a-2.0E0*2*I_KINETIC_Fxyz_Dxz_a;
  abcd[325] = 4.0E0*I_KINETIC_Hx4z_Dxz_aa-2.0E0*2*I_KINETIC_Fx2z_Dxz_a-2.0E0*3*I_KINETIC_Fx2z_Dxz_a+2*1*I_KINETIC_Px_Dxz;
  abcd[326] = 4.0E0*I_KINETIC_H3y2z_Dxz_aa-2.0E0*1*I_KINETIC_F3y_Dxz_a;
  abcd[327] = 4.0E0*I_KINETIC_H2y3z_Dxz_aa-2.0E0*1*I_KINETIC_F2yz_Dxz_a-2.0E0*2*I_KINETIC_F2yz_Dxz_a;
  abcd[328] = 4.0E0*I_KINETIC_Hy4z_Dxz_aa-2.0E0*2*I_KINETIC_Fy2z_Dxz_a-2.0E0*3*I_KINETIC_Fy2z_Dxz_a+2*1*I_KINETIC_Py_Dxz;
  abcd[329] = 4.0E0*I_KINETIC_H5z_Dxz_aa-2.0E0*3*I_KINETIC_F3z_Dxz_a-2.0E0*4*I_KINETIC_F3z_Dxz_a+3*2*I_KINETIC_Pz_Dxz;
  abcd[330] = 4.0E0*I_KINETIC_H3x2z_D2y_aa-2.0E0*1*I_KINETIC_F3x_D2y_a;
  abcd[331] = 4.0E0*I_KINETIC_H2xy2z_D2y_aa-2.0E0*1*I_KINETIC_F2xy_D2y_a;
  abcd[332] = 4.0E0*I_KINETIC_H2x3z_D2y_aa-2.0E0*1*I_KINETIC_F2xz_D2y_a-2.0E0*2*I_KINETIC_F2xz_D2y_a;
  abcd[333] = 4.0E0*I_KINETIC_Hx2y2z_D2y_aa-2.0E0*1*I_KINETIC_Fx2y_D2y_a;
  abcd[334] = 4.0E0*I_KINETIC_Hxy3z_D2y_aa-2.0E0*1*I_KINETIC_Fxyz_D2y_a-2.0E0*2*I_KINETIC_Fxyz_D2y_a;
  abcd[335] = 4.0E0*I_KINETIC_Hx4z_D2y_aa-2.0E0*2*I_KINETIC_Fx2z_D2y_a-2.0E0*3*I_KINETIC_Fx2z_D2y_a+2*1*I_KINETIC_Px_D2y;
  abcd[336] = 4.0E0*I_KINETIC_H3y2z_D2y_aa-2.0E0*1*I_KINETIC_F3y_D2y_a;
  abcd[337] = 4.0E0*I_KINETIC_H2y3z_D2y_aa-2.0E0*1*I_KINETIC_F2yz_D2y_a-2.0E0*2*I_KINETIC_F2yz_D2y_a;
  abcd[338] = 4.0E0*I_KINETIC_Hy4z_D2y_aa-2.0E0*2*I_KINETIC_Fy2z_D2y_a-2.0E0*3*I_KINETIC_Fy2z_D2y_a+2*1*I_KINETIC_Py_D2y;
  abcd[339] = 4.0E0*I_KINETIC_H5z_D2y_aa-2.0E0*3*I_KINETIC_F3z_D2y_a-2.0E0*4*I_KINETIC_F3z_D2y_a+3*2*I_KINETIC_Pz_D2y;
  abcd[340] = 4.0E0*I_KINETIC_H3x2z_Dyz_aa-2.0E0*1*I_KINETIC_F3x_Dyz_a;
  abcd[341] = 4.0E0*I_KINETIC_H2xy2z_Dyz_aa-2.0E0*1*I_KINETIC_F2xy_Dyz_a;
  abcd[342] = 4.0E0*I_KINETIC_H2x3z_Dyz_aa-2.0E0*1*I_KINETIC_F2xz_Dyz_a-2.0E0*2*I_KINETIC_F2xz_Dyz_a;
  abcd[343] = 4.0E0*I_KINETIC_Hx2y2z_Dyz_aa-2.0E0*1*I_KINETIC_Fx2y_Dyz_a;
  abcd[344] = 4.0E0*I_KINETIC_Hxy3z_Dyz_aa-2.0E0*1*I_KINETIC_Fxyz_Dyz_a-2.0E0*2*I_KINETIC_Fxyz_Dyz_a;
  abcd[345] = 4.0E0*I_KINETIC_Hx4z_Dyz_aa-2.0E0*2*I_KINETIC_Fx2z_Dyz_a-2.0E0*3*I_KINETIC_Fx2z_Dyz_a+2*1*I_KINETIC_Px_Dyz;
  abcd[346] = 4.0E0*I_KINETIC_H3y2z_Dyz_aa-2.0E0*1*I_KINETIC_F3y_Dyz_a;
  abcd[347] = 4.0E0*I_KINETIC_H2y3z_Dyz_aa-2.0E0*1*I_KINETIC_F2yz_Dyz_a-2.0E0*2*I_KINETIC_F2yz_Dyz_a;
  abcd[348] = 4.0E0*I_KINETIC_Hy4z_Dyz_aa-2.0E0*2*I_KINETIC_Fy2z_Dyz_a-2.0E0*3*I_KINETIC_Fy2z_Dyz_a+2*1*I_KINETIC_Py_Dyz;
  abcd[349] = 4.0E0*I_KINETIC_H5z_Dyz_aa-2.0E0*3*I_KINETIC_F3z_Dyz_a-2.0E0*4*I_KINETIC_F3z_Dyz_a+3*2*I_KINETIC_Pz_Dyz;
  abcd[350] = 4.0E0*I_KINETIC_H3x2z_D2z_aa-2.0E0*1*I_KINETIC_F3x_D2z_a;
  abcd[351] = 4.0E0*I_KINETIC_H2xy2z_D2z_aa-2.0E0*1*I_KINETIC_F2xy_D2z_a;
  abcd[352] = 4.0E0*I_KINETIC_H2x3z_D2z_aa-2.0E0*1*I_KINETIC_F2xz_D2z_a-2.0E0*2*I_KINETIC_F2xz_D2z_a;
  abcd[353] = 4.0E0*I_KINETIC_Hx2y2z_D2z_aa-2.0E0*1*I_KINETIC_Fx2y_D2z_a;
  abcd[354] = 4.0E0*I_KINETIC_Hxy3z_D2z_aa-2.0E0*1*I_KINETIC_Fxyz_D2z_a-2.0E0*2*I_KINETIC_Fxyz_D2z_a;
  abcd[355] = 4.0E0*I_KINETIC_Hx4z_D2z_aa-2.0E0*2*I_KINETIC_Fx2z_D2z_a-2.0E0*3*I_KINETIC_Fx2z_D2z_a+2*1*I_KINETIC_Px_D2z;
  abcd[356] = 4.0E0*I_KINETIC_H3y2z_D2z_aa-2.0E0*1*I_KINETIC_F3y_D2z_a;
  abcd[357] = 4.0E0*I_KINETIC_H2y3z_D2z_aa-2.0E0*1*I_KINETIC_F2yz_D2z_a-2.0E0*2*I_KINETIC_F2yz_D2z_a;
  abcd[358] = 4.0E0*I_KINETIC_Hy4z_D2z_aa-2.0E0*2*I_KINETIC_Fy2z_D2z_a-2.0E0*3*I_KINETIC_Fy2z_D2z_a+2*1*I_KINETIC_Py_D2z;
  abcd[359] = 4.0E0*I_KINETIC_H5z_D2z_aa-2.0E0*3*I_KINETIC_F3z_D2z_a-2.0E0*4*I_KINETIC_F3z_D2z_a+3*2*I_KINETIC_Pz_D2z;
}
