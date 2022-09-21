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
// BRA1 as redundant position, total RHS integrals evaluated as: 9474
// BRA2 as redundant position, total RHS integrals evaluated as: 5382
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

void hgp_os_kinetic_h_p_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_K7x_Px_aa = 0.0E0;
  Double I_KINETIC_K6xy_Px_aa = 0.0E0;
  Double I_KINETIC_K6xz_Px_aa = 0.0E0;
  Double I_KINETIC_K5x2y_Px_aa = 0.0E0;
  Double I_KINETIC_K5xyz_Px_aa = 0.0E0;
  Double I_KINETIC_K5x2z_Px_aa = 0.0E0;
  Double I_KINETIC_K4x3y_Px_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_Px_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_Px_aa = 0.0E0;
  Double I_KINETIC_K4x3z_Px_aa = 0.0E0;
  Double I_KINETIC_K3x4y_Px_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_Px_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_Px_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_Px_aa = 0.0E0;
  Double I_KINETIC_K3x4z_Px_aa = 0.0E0;
  Double I_KINETIC_K2x5y_Px_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_Px_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_Px_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_Px_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_Px_aa = 0.0E0;
  Double I_KINETIC_K2x5z_Px_aa = 0.0E0;
  Double I_KINETIC_Kx6y_Px_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_Px_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_Px_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_Px_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_Px_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_Px_aa = 0.0E0;
  Double I_KINETIC_Kx6z_Px_aa = 0.0E0;
  Double I_KINETIC_K7y_Px_aa = 0.0E0;
  Double I_KINETIC_K6yz_Px_aa = 0.0E0;
  Double I_KINETIC_K5y2z_Px_aa = 0.0E0;
  Double I_KINETIC_K4y3z_Px_aa = 0.0E0;
  Double I_KINETIC_K3y4z_Px_aa = 0.0E0;
  Double I_KINETIC_K2y5z_Px_aa = 0.0E0;
  Double I_KINETIC_Ky6z_Px_aa = 0.0E0;
  Double I_KINETIC_K7z_Px_aa = 0.0E0;
  Double I_KINETIC_K7x_Py_aa = 0.0E0;
  Double I_KINETIC_K6xy_Py_aa = 0.0E0;
  Double I_KINETIC_K6xz_Py_aa = 0.0E0;
  Double I_KINETIC_K5x2y_Py_aa = 0.0E0;
  Double I_KINETIC_K5xyz_Py_aa = 0.0E0;
  Double I_KINETIC_K5x2z_Py_aa = 0.0E0;
  Double I_KINETIC_K4x3y_Py_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_Py_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_Py_aa = 0.0E0;
  Double I_KINETIC_K4x3z_Py_aa = 0.0E0;
  Double I_KINETIC_K3x4y_Py_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_Py_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_Py_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_Py_aa = 0.0E0;
  Double I_KINETIC_K3x4z_Py_aa = 0.0E0;
  Double I_KINETIC_K2x5y_Py_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_Py_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_Py_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_Py_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_Py_aa = 0.0E0;
  Double I_KINETIC_K2x5z_Py_aa = 0.0E0;
  Double I_KINETIC_Kx6y_Py_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_Py_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_Py_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_Py_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_Py_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_Py_aa = 0.0E0;
  Double I_KINETIC_Kx6z_Py_aa = 0.0E0;
  Double I_KINETIC_K7y_Py_aa = 0.0E0;
  Double I_KINETIC_K6yz_Py_aa = 0.0E0;
  Double I_KINETIC_K5y2z_Py_aa = 0.0E0;
  Double I_KINETIC_K4y3z_Py_aa = 0.0E0;
  Double I_KINETIC_K3y4z_Py_aa = 0.0E0;
  Double I_KINETIC_K2y5z_Py_aa = 0.0E0;
  Double I_KINETIC_Ky6z_Py_aa = 0.0E0;
  Double I_KINETIC_K7z_Py_aa = 0.0E0;
  Double I_KINETIC_K7x_Pz_aa = 0.0E0;
  Double I_KINETIC_K6xy_Pz_aa = 0.0E0;
  Double I_KINETIC_K6xz_Pz_aa = 0.0E0;
  Double I_KINETIC_K5x2y_Pz_aa = 0.0E0;
  Double I_KINETIC_K5xyz_Pz_aa = 0.0E0;
  Double I_KINETIC_K5x2z_Pz_aa = 0.0E0;
  Double I_KINETIC_K4x3y_Pz_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_Pz_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_Pz_aa = 0.0E0;
  Double I_KINETIC_K4x3z_Pz_aa = 0.0E0;
  Double I_KINETIC_K3x4y_Pz_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_Pz_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_Pz_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_Pz_aa = 0.0E0;
  Double I_KINETIC_K3x4z_Pz_aa = 0.0E0;
  Double I_KINETIC_K2x5y_Pz_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_Pz_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_Pz_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_Pz_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_Pz_aa = 0.0E0;
  Double I_KINETIC_K2x5z_Pz_aa = 0.0E0;
  Double I_KINETIC_Kx6y_Pz_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_Pz_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_Pz_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_Pz_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_Pz_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_Pz_aa = 0.0E0;
  Double I_KINETIC_Kx6z_Pz_aa = 0.0E0;
  Double I_KINETIC_K7y_Pz_aa = 0.0E0;
  Double I_KINETIC_K6yz_Pz_aa = 0.0E0;
  Double I_KINETIC_K5y2z_Pz_aa = 0.0E0;
  Double I_KINETIC_K4y3z_Pz_aa = 0.0E0;
  Double I_KINETIC_K3y4z_Pz_aa = 0.0E0;
  Double I_KINETIC_K2y5z_Pz_aa = 0.0E0;
  Double I_KINETIC_Ky6z_Pz_aa = 0.0E0;
  Double I_KINETIC_K7z_Pz_aa = 0.0E0;
  Double I_KINETIC_H5x_Px_a = 0.0E0;
  Double I_KINETIC_H4xy_Px_a = 0.0E0;
  Double I_KINETIC_H4xz_Px_a = 0.0E0;
  Double I_KINETIC_H3x2y_Px_a = 0.0E0;
  Double I_KINETIC_H3xyz_Px_a = 0.0E0;
  Double I_KINETIC_H3x2z_Px_a = 0.0E0;
  Double I_KINETIC_H2x3y_Px_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Px_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Px_a = 0.0E0;
  Double I_KINETIC_H2x3z_Px_a = 0.0E0;
  Double I_KINETIC_Hx4y_Px_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Px_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Px_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Px_a = 0.0E0;
  Double I_KINETIC_Hx4z_Px_a = 0.0E0;
  Double I_KINETIC_H5y_Px_a = 0.0E0;
  Double I_KINETIC_H4yz_Px_a = 0.0E0;
  Double I_KINETIC_H3y2z_Px_a = 0.0E0;
  Double I_KINETIC_H2y3z_Px_a = 0.0E0;
  Double I_KINETIC_Hy4z_Px_a = 0.0E0;
  Double I_KINETIC_H5z_Px_a = 0.0E0;
  Double I_KINETIC_H5x_Py_a = 0.0E0;
  Double I_KINETIC_H4xy_Py_a = 0.0E0;
  Double I_KINETIC_H4xz_Py_a = 0.0E0;
  Double I_KINETIC_H3x2y_Py_a = 0.0E0;
  Double I_KINETIC_H3xyz_Py_a = 0.0E0;
  Double I_KINETIC_H3x2z_Py_a = 0.0E0;
  Double I_KINETIC_H2x3y_Py_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Py_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Py_a = 0.0E0;
  Double I_KINETIC_H2x3z_Py_a = 0.0E0;
  Double I_KINETIC_Hx4y_Py_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Py_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Py_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Py_a = 0.0E0;
  Double I_KINETIC_Hx4z_Py_a = 0.0E0;
  Double I_KINETIC_H5y_Py_a = 0.0E0;
  Double I_KINETIC_H4yz_Py_a = 0.0E0;
  Double I_KINETIC_H3y2z_Py_a = 0.0E0;
  Double I_KINETIC_H2y3z_Py_a = 0.0E0;
  Double I_KINETIC_Hy4z_Py_a = 0.0E0;
  Double I_KINETIC_H5z_Py_a = 0.0E0;
  Double I_KINETIC_H5x_Pz_a = 0.0E0;
  Double I_KINETIC_H4xy_Pz_a = 0.0E0;
  Double I_KINETIC_H4xz_Pz_a = 0.0E0;
  Double I_KINETIC_H3x2y_Pz_a = 0.0E0;
  Double I_KINETIC_H3xyz_Pz_a = 0.0E0;
  Double I_KINETIC_H3x2z_Pz_a = 0.0E0;
  Double I_KINETIC_H2x3y_Pz_a = 0.0E0;
  Double I_KINETIC_H2x2yz_Pz_a = 0.0E0;
  Double I_KINETIC_H2xy2z_Pz_a = 0.0E0;
  Double I_KINETIC_H2x3z_Pz_a = 0.0E0;
  Double I_KINETIC_Hx4y_Pz_a = 0.0E0;
  Double I_KINETIC_Hx3yz_Pz_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_Pz_a = 0.0E0;
  Double I_KINETIC_Hxy3z_Pz_a = 0.0E0;
  Double I_KINETIC_Hx4z_Pz_a = 0.0E0;
  Double I_KINETIC_H5y_Pz_a = 0.0E0;
  Double I_KINETIC_H4yz_Pz_a = 0.0E0;
  Double I_KINETIC_H3y2z_Pz_a = 0.0E0;
  Double I_KINETIC_H2y3z_Pz_a = 0.0E0;
  Double I_KINETIC_Hy4z_Pz_a = 0.0E0;
  Double I_KINETIC_H5z_Pz_a = 0.0E0;
  Double I_KINETIC_F3x_Px = 0.0E0;
  Double I_KINETIC_F2xy_Px = 0.0E0;
  Double I_KINETIC_F2xz_Px = 0.0E0;
  Double I_KINETIC_Fx2y_Px = 0.0E0;
  Double I_KINETIC_Fxyz_Px = 0.0E0;
  Double I_KINETIC_Fx2z_Px = 0.0E0;
  Double I_KINETIC_F3y_Px = 0.0E0;
  Double I_KINETIC_F2yz_Px = 0.0E0;
  Double I_KINETIC_Fy2z_Px = 0.0E0;
  Double I_KINETIC_F3z_Px = 0.0E0;
  Double I_KINETIC_F3x_Py = 0.0E0;
  Double I_KINETIC_F2xy_Py = 0.0E0;
  Double I_KINETIC_F2xz_Py = 0.0E0;
  Double I_KINETIC_Fx2y_Py = 0.0E0;
  Double I_KINETIC_Fxyz_Py = 0.0E0;
  Double I_KINETIC_Fx2z_Py = 0.0E0;
  Double I_KINETIC_F3y_Py = 0.0E0;
  Double I_KINETIC_F2yz_Py = 0.0E0;
  Double I_KINETIC_Fy2z_Py = 0.0E0;
  Double I_KINETIC_F3z_Py = 0.0E0;
  Double I_KINETIC_F3x_Pz = 0.0E0;
  Double I_KINETIC_F2xy_Pz = 0.0E0;
  Double I_KINETIC_F2xz_Pz = 0.0E0;
  Double I_KINETIC_Fx2y_Pz = 0.0E0;
  Double I_KINETIC_Fxyz_Pz = 0.0E0;
  Double I_KINETIC_Fx2z_Pz = 0.0E0;
  Double I_KINETIC_F3y_Pz = 0.0E0;
  Double I_KINETIC_F2yz_Pz = 0.0E0;
  Double I_KINETIC_Fy2z_Pz = 0.0E0;
  Double I_KINETIC_F3z_Pz = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_Px_vrr = PBX*I_TWOBODYOVERLAP_H5x_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Px_vrr = PBX*I_TWOBODYOVERLAP_H4xy_S_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Px_vrr = PBX*I_TWOBODYOVERLAP_H4xz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Px_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Px_vrr = PBX*I_TWOBODYOVERLAP_H3xyz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Px_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Px_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_H2xy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Px_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Px_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Px_vrr = PBX*I_TWOBODYOVERLAP_Hx3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Hx2y2z_S_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Hxy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Px_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_Px_vrr = PBX*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Px_vrr = PBX*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Px_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Px_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_Px_vrr = PBX*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_H5x_Py_vrr = PBY*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Py_vrr = PBY*I_TWOBODYOVERLAP_H4xy_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Py_vrr = PBY*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Py_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Py_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Py_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Py_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Py_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Py_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Py_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Py_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_Py_vrr = PBY*I_TWOBODYOVERLAP_H5y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Py_vrr = PBY*I_TWOBODYOVERLAP_H4yz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Py_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Py_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_Py_vrr = PBY*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_H5x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_K_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_K7x_Px_vrr = PBX*I_TWOBODYOVERLAP_K7x_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xy_Px_vrr = PBX*I_TWOBODYOVERLAP_K6xy_S_vrr+6*oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_K6xz_Px_vrr = PBX*I_TWOBODYOVERLAP_K6xz_S_vrr+6*oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_Px_vrr = PBX*I_TWOBODYOVERLAP_K5x2y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_Px_vrr = PBX*I_TWOBODYOVERLAP_K5xyz_S_vrr+5*oned2z*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_Px_vrr = PBX*I_TWOBODYOVERLAP_K5x2z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_Px_vrr = PBX*I_TWOBODYOVERLAP_K4x3y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_K4x2yz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_K4xy2z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_Px_vrr = PBX*I_TWOBODYOVERLAP_K4x3z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_Px_vrr = PBX*I_TWOBODYOVERLAP_K3x4y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_Px_vrr = PBX*I_TWOBODYOVERLAP_K3x3yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_K3x2y2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_Px_vrr = PBX*I_TWOBODYOVERLAP_K3xy3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_Px_vrr = PBX*I_TWOBODYOVERLAP_K3x4z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_Px_vrr = PBX*I_TWOBODYOVERLAP_K2x5y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_Px_vrr = PBX*I_TWOBODYOVERLAP_K2x4yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_K2x3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_Px_vrr = PBX*I_TWOBODYOVERLAP_K2x2y3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_Px_vrr = PBX*I_TWOBODYOVERLAP_K2xy4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_Px_vrr = PBX*I_TWOBODYOVERLAP_K2x5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_Px_vrr = PBX*I_TWOBODYOVERLAP_Kx6y_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_Px_vrr = PBX*I_TWOBODYOVERLAP_Kx5yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Kx4y2z_S_vrr+oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Kx3y3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_Px_vrr = PBX*I_TWOBODYOVERLAP_Kx2y4z_S_vrr+oned2z*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_Px_vrr = PBX*I_TWOBODYOVERLAP_Kxy5z_S_vrr+oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_Px_vrr = PBX*I_TWOBODYOVERLAP_Kx6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7y_Px_vrr = PBX*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_K6yz_Px_vrr = PBX*I_TWOBODYOVERLAP_K6yz_S_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_Px_vrr = PBX*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_Px_vrr = PBX*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_Px_vrr = PBX*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_Px_vrr = PBX*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7z_Px_vrr = PBX*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_K7x_Py_vrr = PBY*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xy_Py_vrr = PBY*I_TWOBODYOVERLAP_K6xy_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xz_Py_vrr = PBY*I_TWOBODYOVERLAP_K6xz_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_Py_vrr = PBY*I_TWOBODYOVERLAP_K5x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_Py_vrr = PBY*I_TWOBODYOVERLAP_K5xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_Py_vrr = PBY*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_Py_vrr = PBY*I_TWOBODYOVERLAP_K4x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_K4x2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_K4xy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_Py_vrr = PBY*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_Py_vrr = PBY*I_TWOBODYOVERLAP_K3x4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_Py_vrr = PBY*I_TWOBODYOVERLAP_K3x3yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_K3x2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_Py_vrr = PBY*I_TWOBODYOVERLAP_K3xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_Py_vrr = PBY*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_Py_vrr = PBY*I_TWOBODYOVERLAP_K2x5y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_Py_vrr = PBY*I_TWOBODYOVERLAP_K2x4yz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_K2x3y2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_Py_vrr = PBY*I_TWOBODYOVERLAP_K2x2y3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_Py_vrr = PBY*I_TWOBODYOVERLAP_K2xy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_Py_vrr = PBY*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_Py_vrr = PBY*I_TWOBODYOVERLAP_Kx6y_S_vrr+6*oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_Py_vrr = PBY*I_TWOBODYOVERLAP_Kx5yz_S_vrr+5*oned2z*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Kx4y2z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Kx3y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_Py_vrr = PBY*I_TWOBODYOVERLAP_Kx2y4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_Py_vrr = PBY*I_TWOBODYOVERLAP_Kxy5z_S_vrr+oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_Py_vrr = PBY*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7y_Py_vrr = PBY*I_TWOBODYOVERLAP_K7y_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_K6yz_Py_vrr = PBY*I_TWOBODYOVERLAP_K6yz_S_vrr+6*oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_K5y2z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_Py_vrr = PBY*I_TWOBODYOVERLAP_K4y3z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_Py_vrr = PBY*I_TWOBODYOVERLAP_K3y4z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_Py_vrr = PBY*I_TWOBODYOVERLAP_K2y5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_Py_vrr = PBY*I_TWOBODYOVERLAP_Ky6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7z_Py_vrr = PBY*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_K7x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K6xy_S_vrr;
    Double I_TWOBODYOVERLAP_K6xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K6xz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K5xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K5x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K4x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K4xy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K4x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K3x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K3x2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K3xy3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K3x4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K2x4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K2x3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K2x2y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K2xy4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K2x5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Kx5yz_S_vrr+oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Kx4y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Kx3y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Kx2y4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Kxy5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Kx6z_S_vrr+6*oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_K7y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_K6yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K6yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K5y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K4y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K3y4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K2y5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ky6z_S_vrr+6*oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_K7z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_K7z_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;

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
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_KINETIC_G4x_S_vrr = PAX*I_KINETIC_F3x_S_vrr+3*oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G3xy_S_vrr = PAY*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_S_vrr = PAZ*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_S_vrr = PAY*I_KINETIC_F2xy_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G2xyz_S_vrr = PAZ*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_S_vrr = PAZ*I_KINETIC_F2xz_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_Gx3y_S_vrr = PAX*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_S_vrr = PAZ*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_S_vrr = PAY*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_KINETIC_Gx3z_S_vrr = PAX*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_S_vrr = PAY*I_KINETIC_F3y_S_vrr+3*oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_G3yz_S_vrr = PAZ*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_S_vrr = PAZ*I_KINETIC_F2yz_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_Gy3z_S_vrr = PAY*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_S_vrr = PAZ*I_KINETIC_F3z_S_vrr+3*oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_KINETIC_H5x_S_vrr = PAX*I_KINETIC_G4x_S_vrr+4*oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H4xy_S_vrr = PAY*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_S_vrr = PAZ*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_S_vrr = PAY*I_KINETIC_G3xy_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_S_vrr-bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H3xyz_S_vrr = PAZ*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_KINETIC_H3x2z_S_vrr = PAZ*I_KINETIC_G3xz_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_S_vrr-bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H2x3y_S_vrr = PAX*I_KINETIC_Gx3y_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_S_vrr-bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H2x2yz_S_vrr = PAZ*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2xy2z_S_vrr = PAY*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_KINETIC_H2x3z_S_vrr = PAX*I_KINETIC_Gx3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_S_vrr-bdz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_Hx4y_S_vrr = PAX*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx3yz_S_vrr = PAZ*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_KINETIC_Hx2y2z_S_vrr = PAX*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_KINETIC_Hxy3z_S_vrr = PAY*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_KINETIC_Hx4z_S_vrr = PAX*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_S_vrr = PAY*I_KINETIC_G4y_S_vrr+4*oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H4yz_S_vrr = PAZ*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_S_vrr = PAZ*I_KINETIC_G3yz_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_S_vrr-bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H2y3z_S_vrr = PAY*I_KINETIC_Gy3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_S_vrr-bdz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_Hy4z_S_vrr = PAY*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_S_vrr = PAZ*I_KINETIC_G4z_S_vrr+4*oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     ************************************************************/
    Double I_KINETIC_H5x_Px_vrr = PBX*I_KINETIC_H5x_S_vrr+5*oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_KINETIC_H4xy_Px_vrr = PBX*I_KINETIC_H4xy_S_vrr+4*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_KINETIC_H4xz_Px_vrr = PBX*I_KINETIC_H4xz_S_vrr+4*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_KINETIC_H3x2y_Px_vrr = PBX*I_KINETIC_H3x2y_S_vrr+3*oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_KINETIC_H3xyz_Px_vrr = PBX*I_KINETIC_H3xyz_S_vrr+3*oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Px_vrr;
    Double I_KINETIC_H3x2z_Px_vrr = PBX*I_KINETIC_H3x2z_S_vrr+3*oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_KINETIC_H2x3y_Px_vrr = PBX*I_KINETIC_H2x3y_S_vrr+2*oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_KINETIC_H2x2yz_Px_vrr = PBX*I_KINETIC_H2x2yz_S_vrr+2*oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_KINETIC_H2xy2z_Px_vrr = PBX*I_KINETIC_H2xy2z_S_vrr+2*oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Px_vrr;
    Double I_KINETIC_H2x3z_Px_vrr = PBX*I_KINETIC_H2x3z_S_vrr+2*oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_KINETIC_Hx4y_Px_vrr = PBX*I_KINETIC_Hx4y_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_KINETIC_Hx3yz_Px_vrr = PBX*I_KINETIC_Hx3yz_S_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Px_vrr;
    Double I_KINETIC_Hx2y2z_Px_vrr = PBX*I_KINETIC_Hx2y2z_S_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Px_vrr;
    Double I_KINETIC_Hxy3z_Px_vrr = PBX*I_KINETIC_Hxy3z_S_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Px_vrr;
    Double I_KINETIC_Hx4z_Px_vrr = PBX*I_KINETIC_Hx4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_KINETIC_H5y_Px_vrr = PBX*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_KINETIC_H4yz_Px_vrr = PBX*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_KINETIC_H3y2z_Px_vrr = PBX*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_KINETIC_H2y3z_Px_vrr = PBX*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_KINETIC_Hy4z_Px_vrr = PBX*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_KINETIC_H5z_Px_vrr = PBX*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_KINETIC_H5x_Py_vrr = PBY*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_KINETIC_H4xy_Py_vrr = PBY*I_KINETIC_H4xy_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_KINETIC_H4xz_Py_vrr = PBY*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_KINETIC_H3x2y_Py_vrr = PBY*I_KINETIC_H3x2y_S_vrr+2*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_KINETIC_H3xyz_Py_vrr = PBY*I_KINETIC_H3xyz_S_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Py_vrr;
    Double I_KINETIC_H3x2z_Py_vrr = PBY*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_KINETIC_H2x3y_Py_vrr = PBY*I_KINETIC_H2x3y_S_vrr+3*oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_KINETIC_H2x2yz_Py_vrr = PBY*I_KINETIC_H2x2yz_S_vrr+2*oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Py_vrr;
    Double I_KINETIC_H2xy2z_Py_vrr = PBY*I_KINETIC_H2xy2z_S_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Py_vrr;
    Double I_KINETIC_H2x3z_Py_vrr = PBY*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Py_vrr;
    Double I_KINETIC_Hx4y_Py_vrr = PBY*I_KINETIC_Hx4y_S_vrr+4*oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_KINETIC_Hx3yz_Py_vrr = PBY*I_KINETIC_Hx3yz_S_vrr+3*oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Py_vrr;
    Double I_KINETIC_Hx2y2z_Py_vrr = PBY*I_KINETIC_Hx2y2z_S_vrr+2*oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Py_vrr;
    Double I_KINETIC_Hxy3z_Py_vrr = PBY*I_KINETIC_Hxy3z_S_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Py_vrr;
    Double I_KINETIC_Hx4z_Py_vrr = PBY*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_KINETIC_H5y_Py_vrr = PBY*I_KINETIC_H5y_S_vrr+5*oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_KINETIC_H4yz_Py_vrr = PBY*I_KINETIC_H4yz_S_vrr+4*oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_KINETIC_H3y2z_Py_vrr = PBY*I_KINETIC_H3y2z_S_vrr+3*oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_KINETIC_H2y3z_Py_vrr = PBY*I_KINETIC_H2y3z_S_vrr+2*oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Py_vrr;
    Double I_KINETIC_Hy4z_Py_vrr = PBY*I_KINETIC_Hy4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_KINETIC_H5z_Py_vrr = PBY*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_KINETIC_H5x_Pz_vrr = PBZ*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_KINETIC_H4xy_Pz_vrr = PBZ*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_KINETIC_H4xz_Pz_vrr = PBZ*I_KINETIC_H4xz_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Pz_vrr;
    Double I_KINETIC_H3x2y_Pz_vrr = PBZ*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Pz_vrr;
    Double I_KINETIC_H3xyz_Pz_vrr = PBZ*I_KINETIC_H3xyz_S_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Pz_vrr;
    Double I_KINETIC_H3x2z_Pz_vrr = PBZ*I_KINETIC_H3x2z_S_vrr+2*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Pz_vrr;
    Double I_KINETIC_H2x3y_Pz_vrr = PBZ*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Pz_vrr;
    Double I_KINETIC_H2x2yz_Pz_vrr = PBZ*I_KINETIC_H2x2yz_S_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr;
    Double I_KINETIC_H2xy2z_Pz_vrr = PBZ*I_KINETIC_H2xy2z_S_vrr+2*oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Pz_vrr;
    Double I_KINETIC_H2x3z_Pz_vrr = PBZ*I_KINETIC_H2x3z_S_vrr+3*oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Pz_vrr;
    Double I_KINETIC_Hx4y_Pz_vrr = PBZ*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_KINETIC_Hx3yz_Pz_vrr = PBZ*I_KINETIC_Hx3yz_S_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Pz_vrr;
    Double I_KINETIC_Hx2y2z_Pz_vrr = PBZ*I_KINETIC_Hx2y2z_S_vrr+2*oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr;
    Double I_KINETIC_Hxy3z_Pz_vrr = PBZ*I_KINETIC_Hxy3z_S_vrr+3*oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Pz_vrr;
    Double I_KINETIC_Hx4z_Pz_vrr = PBZ*I_KINETIC_Hx4z_S_vrr+4*oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_KINETIC_H5y_Pz_vrr = PBZ*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_KINETIC_H4yz_Pz_vrr = PBZ*I_KINETIC_H4yz_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Pz_vrr;
    Double I_KINETIC_H3y2z_Pz_vrr = PBZ*I_KINETIC_H3y2z_S_vrr+2*oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Pz_vrr;
    Double I_KINETIC_H2y3z_Pz_vrr = PBZ*I_KINETIC_H2y3z_S_vrr+3*oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Pz_vrr;
    Double I_KINETIC_Hy4z_Pz_vrr = PBZ*I_KINETIC_Hy4z_S_vrr+4*oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Pz_vrr;
    Double I_KINETIC_H5z_Pz_vrr = PBZ*I_KINETIC_H5z_S_vrr+5*oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_KINETIC_I6x_S_vrr = PAX*I_KINETIC_H5x_S_vrr+5*oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I5xy_S_vrr = PAY*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_KINETIC_I5xz_S_vrr = PAZ*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_KINETIC_I4x2y_S_vrr = PAY*I_KINETIC_H4xy_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_S_vrr-bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I4xyz_S_vrr = PAZ*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    Double I_KINETIC_I4x2z_S_vrr = PAZ*I_KINETIC_H4xz_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_S_vrr-bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I3x3y_S_vrr = PAY*I_KINETIC_H3x2y_S_vrr+2*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_I3x2yz_S_vrr = PAZ*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_KINETIC_I3xy2z_S_vrr = PAY*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    Double I_KINETIC_I3x3z_S_vrr = PAZ*I_KINETIC_H3x2z_S_vrr+2*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_I2x4y_S_vrr = PAX*I_KINETIC_Hx4y_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_S_vrr-bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I2x3yz_S_vrr = PAZ*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_KINETIC_I2x2y2z_S_vrr = PAZ*I_KINETIC_H2x2yz_S_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_S_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_I2xy3z_S_vrr = PAY*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_KINETIC_I2x4z_S_vrr = PAX*I_KINETIC_Hx4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_S_vrr-bdz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_Ix5y_S_vrr = PAX*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_KINETIC_Ix4yz_S_vrr = PAZ*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    Double I_KINETIC_Ix3y2z_S_vrr = PAX*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    Double I_KINETIC_Ix2y3z_S_vrr = PAX*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    Double I_KINETIC_Ixy4z_S_vrr = PAY*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    Double I_KINETIC_Ix5z_S_vrr = PAX*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_KINETIC_I6y_S_vrr = PAY*I_KINETIC_H5y_S_vrr+5*oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I5yz_S_vrr = PAZ*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_KINETIC_I4y2z_S_vrr = PAZ*I_KINETIC_H4yz_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_S_vrr-bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I3y3z_S_vrr = PAZ*I_KINETIC_H3y2z_S_vrr+2*oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_I2y4z_S_vrr = PAY*I_KINETIC_Hy4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_S_vrr-bdz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_Iy5z_S_vrr = PAY*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_KINETIC_I6z_S_vrr = PAZ*I_KINETIC_H5z_S_vrr+5*oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_K_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_I_S
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_KINETIC_K7x_S_vrr = PAX*I_KINETIC_I6x_S_vrr+6*oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_K7x_S_vrr-6*bdz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_K6xy_S_vrr = PAY*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_S_vrr;
    Double I_KINETIC_K6xz_S_vrr = PAZ*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_S_vrr;
    Double I_KINETIC_K5x2y_S_vrr = PAY*I_KINETIC_I5xy_S_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_S_vrr-bdz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_K5xyz_S_vrr = PAZ*I_KINETIC_I5xy_S_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_S_vrr;
    Double I_KINETIC_K5x2z_S_vrr = PAZ*I_KINETIC_I5xz_S_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_S_vrr-bdz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_K4x3y_S_vrr = PAY*I_KINETIC_I4x2y_S_vrr+2*oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_K4x2yz_S_vrr = PAZ*I_KINETIC_I4x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_S_vrr;
    Double I_KINETIC_K4xy2z_S_vrr = PAY*I_KINETIC_I4x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_S_vrr;
    Double I_KINETIC_K4x3z_S_vrr = PAZ*I_KINETIC_I4x2z_S_vrr+2*oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_K3x4y_S_vrr = PAX*I_KINETIC_I2x4y_S_vrr+2*oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_S_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_K3x3yz_S_vrr = PAZ*I_KINETIC_I3x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    Double I_KINETIC_K3x2y2z_S_vrr = PAZ*I_KINETIC_I3x2yz_S_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_S_vrr-bdz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_K3xy3z_S_vrr = PAY*I_KINETIC_I3x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_S_vrr;
    Double I_KINETIC_K3x4z_S_vrr = PAX*I_KINETIC_I2x4z_S_vrr+2*oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_S_vrr-2*bdz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_K2x5y_S_vrr = PAX*I_KINETIC_Ix5y_S_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_S_vrr-bdz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_K2x4yz_S_vrr = PAZ*I_KINETIC_I2x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_S_vrr;
    Double I_KINETIC_K2x3y2z_S_vrr = PAZ*I_KINETIC_I2x3yz_S_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_S_vrr-bdz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_K2x2y3z_S_vrr = PAY*I_KINETIC_I2xy3z_S_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_S_vrr-bdz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_K2xy4z_S_vrr = PAY*I_KINETIC_I2x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_S_vrr;
    Double I_KINETIC_K2x5z_S_vrr = PAX*I_KINETIC_Ix5z_S_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_S_vrr-bdz*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_KINETIC_Kx6y_S_vrr = PAX*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    Double I_KINETIC_Kx5yz_S_vrr = PAZ*I_KINETIC_Ix5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_S_vrr;
    Double I_KINETIC_Kx4y2z_S_vrr = PAX*I_KINETIC_I4y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_S_vrr;
    Double I_KINETIC_Kx3y3z_S_vrr = PAX*I_KINETIC_I3y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_S_vrr;
    Double I_KINETIC_Kx2y4z_S_vrr = PAX*I_KINETIC_I2y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_S_vrr;
    Double I_KINETIC_Kxy5z_S_vrr = PAY*I_KINETIC_Ix5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_S_vrr;
    Double I_KINETIC_Kx6z_S_vrr = PAX*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    Double I_KINETIC_K7y_S_vrr = PAY*I_KINETIC_I6y_S_vrr+6*oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_K7y_S_vrr-6*bdz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_K6yz_S_vrr = PAZ*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_S_vrr;
    Double I_KINETIC_K5y2z_S_vrr = PAZ*I_KINETIC_I5yz_S_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_S_vrr-bdz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_K4y3z_S_vrr = PAZ*I_KINETIC_I4y2z_S_vrr+2*oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_K3y4z_S_vrr = PAY*I_KINETIC_I2y4z_S_vrr+2*oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_S_vrr-2*bdz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_K2y5z_S_vrr = PAY*I_KINETIC_Iy5z_S_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_S_vrr-bdz*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_KINETIC_Ky6z_S_vrr = PAY*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    Double I_KINETIC_K7z_S_vrr = PAZ*I_KINETIC_I6z_S_vrr+6*oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_K7z_S_vrr-6*bdz*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_K_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_K_S
     * RHS shell quartet name: SQ_KINETIC_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P
     ************************************************************/
    Double I_KINETIC_K7x_Px_vrr = PBX*I_KINETIC_K7x_S_vrr+7*oned2z*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_K7x_Px_vrr;
    Double I_KINETIC_K6xy_Px_vrr = PBX*I_KINETIC_K6xy_S_vrr+6*oned2z*I_KINETIC_I5xy_S_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_Px_vrr;
    Double I_KINETIC_K6xz_Px_vrr = PBX*I_KINETIC_K6xz_S_vrr+6*oned2z*I_KINETIC_I5xz_S_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_Px_vrr;
    Double I_KINETIC_K5x2y_Px_vrr = PBX*I_KINETIC_K5x2y_S_vrr+5*oned2z*I_KINETIC_I4x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_Px_vrr;
    Double I_KINETIC_K5xyz_Px_vrr = PBX*I_KINETIC_K5xyz_S_vrr+5*oned2z*I_KINETIC_I4xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_Px_vrr;
    Double I_KINETIC_K5x2z_Px_vrr = PBX*I_KINETIC_K5x2z_S_vrr+5*oned2z*I_KINETIC_I4x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_Px_vrr;
    Double I_KINETIC_K4x3y_Px_vrr = PBX*I_KINETIC_K4x3y_S_vrr+4*oned2z*I_KINETIC_I3x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_Px_vrr;
    Double I_KINETIC_K4x2yz_Px_vrr = PBX*I_KINETIC_K4x2yz_S_vrr+4*oned2z*I_KINETIC_I3x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_Px_vrr;
    Double I_KINETIC_K4xy2z_Px_vrr = PBX*I_KINETIC_K4xy2z_S_vrr+4*oned2z*I_KINETIC_I3xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_Px_vrr;
    Double I_KINETIC_K4x3z_Px_vrr = PBX*I_KINETIC_K4x3z_S_vrr+4*oned2z*I_KINETIC_I3x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_Px_vrr;
    Double I_KINETIC_K3x4y_Px_vrr = PBX*I_KINETIC_K3x4y_S_vrr+3*oned2z*I_KINETIC_I2x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_Px_vrr;
    Double I_KINETIC_K3x3yz_Px_vrr = PBX*I_KINETIC_K3x3yz_S_vrr+3*oned2z*I_KINETIC_I2x3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_Px_vrr;
    Double I_KINETIC_K3x2y2z_Px_vrr = PBX*I_KINETIC_K3x2y2z_S_vrr+3*oned2z*I_KINETIC_I2x2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_Px_vrr;
    Double I_KINETIC_K3xy3z_Px_vrr = PBX*I_KINETIC_K3xy3z_S_vrr+3*oned2z*I_KINETIC_I2xy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_Px_vrr;
    Double I_KINETIC_K3x4z_Px_vrr = PBX*I_KINETIC_K3x4z_S_vrr+3*oned2z*I_KINETIC_I2x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_Px_vrr;
    Double I_KINETIC_K2x5y_Px_vrr = PBX*I_KINETIC_K2x5y_S_vrr+2*oned2z*I_KINETIC_Ix5y_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_Px_vrr;
    Double I_KINETIC_K2x4yz_Px_vrr = PBX*I_KINETIC_K2x4yz_S_vrr+2*oned2z*I_KINETIC_Ix4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_Px_vrr;
    Double I_KINETIC_K2x3y2z_Px_vrr = PBX*I_KINETIC_K2x3y2z_S_vrr+2*oned2z*I_KINETIC_Ix3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_Px_vrr;
    Double I_KINETIC_K2x2y3z_Px_vrr = PBX*I_KINETIC_K2x2y3z_S_vrr+2*oned2z*I_KINETIC_Ix2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_Px_vrr;
    Double I_KINETIC_K2xy4z_Px_vrr = PBX*I_KINETIC_K2xy4z_S_vrr+2*oned2z*I_KINETIC_Ixy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_Px_vrr;
    Double I_KINETIC_K2x5z_Px_vrr = PBX*I_KINETIC_K2x5z_S_vrr+2*oned2z*I_KINETIC_Ix5z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_Px_vrr;
    Double I_KINETIC_Kx6y_Px_vrr = PBX*I_KINETIC_Kx6y_S_vrr+oned2z*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_Px_vrr;
    Double I_KINETIC_Kx5yz_Px_vrr = PBX*I_KINETIC_Kx5yz_S_vrr+oned2z*I_KINETIC_I5yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_Px_vrr;
    Double I_KINETIC_Kx4y2z_Px_vrr = PBX*I_KINETIC_Kx4y2z_S_vrr+oned2z*I_KINETIC_I4y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_Px_vrr;
    Double I_KINETIC_Kx3y3z_Px_vrr = PBX*I_KINETIC_Kx3y3z_S_vrr+oned2z*I_KINETIC_I3y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_Px_vrr;
    Double I_KINETIC_Kx2y4z_Px_vrr = PBX*I_KINETIC_Kx2y4z_S_vrr+oned2z*I_KINETIC_I2y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_Px_vrr;
    Double I_KINETIC_Kxy5z_Px_vrr = PBX*I_KINETIC_Kxy5z_S_vrr+oned2z*I_KINETIC_Iy5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_Px_vrr;
    Double I_KINETIC_Kx6z_Px_vrr = PBX*I_KINETIC_Kx6z_S_vrr+oned2z*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_Px_vrr;
    Double I_KINETIC_K7y_Px_vrr = PBX*I_KINETIC_K7y_S_vrr+twoxi*I_TWOBODYOVERLAP_K7y_Px_vrr;
    Double I_KINETIC_K6yz_Px_vrr = PBX*I_KINETIC_K6yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_Px_vrr;
    Double I_KINETIC_K5y2z_Px_vrr = PBX*I_KINETIC_K5y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_Px_vrr;
    Double I_KINETIC_K4y3z_Px_vrr = PBX*I_KINETIC_K4y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_Px_vrr;
    Double I_KINETIC_K3y4z_Px_vrr = PBX*I_KINETIC_K3y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_Px_vrr;
    Double I_KINETIC_K2y5z_Px_vrr = PBX*I_KINETIC_K2y5z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_Px_vrr;
    Double I_KINETIC_Ky6z_Px_vrr = PBX*I_KINETIC_Ky6z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_Px_vrr;
    Double I_KINETIC_K7z_Px_vrr = PBX*I_KINETIC_K7z_S_vrr+twoxi*I_TWOBODYOVERLAP_K7z_Px_vrr;
    Double I_KINETIC_K7x_Py_vrr = PBY*I_KINETIC_K7x_S_vrr+twoxi*I_TWOBODYOVERLAP_K7x_Py_vrr;
    Double I_KINETIC_K6xy_Py_vrr = PBY*I_KINETIC_K6xy_S_vrr+oned2z*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_Py_vrr;
    Double I_KINETIC_K6xz_Py_vrr = PBY*I_KINETIC_K6xz_S_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_Py_vrr;
    Double I_KINETIC_K5x2y_Py_vrr = PBY*I_KINETIC_K5x2y_S_vrr+2*oned2z*I_KINETIC_I5xy_S_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_Py_vrr;
    Double I_KINETIC_K5xyz_Py_vrr = PBY*I_KINETIC_K5xyz_S_vrr+oned2z*I_KINETIC_I5xz_S_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_Py_vrr;
    Double I_KINETIC_K5x2z_Py_vrr = PBY*I_KINETIC_K5x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_Py_vrr;
    Double I_KINETIC_K4x3y_Py_vrr = PBY*I_KINETIC_K4x3y_S_vrr+3*oned2z*I_KINETIC_I4x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_Py_vrr;
    Double I_KINETIC_K4x2yz_Py_vrr = PBY*I_KINETIC_K4x2yz_S_vrr+2*oned2z*I_KINETIC_I4xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_Py_vrr;
    Double I_KINETIC_K4xy2z_Py_vrr = PBY*I_KINETIC_K4xy2z_S_vrr+oned2z*I_KINETIC_I4x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_Py_vrr;
    Double I_KINETIC_K4x3z_Py_vrr = PBY*I_KINETIC_K4x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_Py_vrr;
    Double I_KINETIC_K3x4y_Py_vrr = PBY*I_KINETIC_K3x4y_S_vrr+4*oned2z*I_KINETIC_I3x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_Py_vrr;
    Double I_KINETIC_K3x3yz_Py_vrr = PBY*I_KINETIC_K3x3yz_S_vrr+3*oned2z*I_KINETIC_I3x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_Py_vrr;
    Double I_KINETIC_K3x2y2z_Py_vrr = PBY*I_KINETIC_K3x2y2z_S_vrr+2*oned2z*I_KINETIC_I3xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_Py_vrr;
    Double I_KINETIC_K3xy3z_Py_vrr = PBY*I_KINETIC_K3xy3z_S_vrr+oned2z*I_KINETIC_I3x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_Py_vrr;
    Double I_KINETIC_K3x4z_Py_vrr = PBY*I_KINETIC_K3x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_Py_vrr;
    Double I_KINETIC_K2x5y_Py_vrr = PBY*I_KINETIC_K2x5y_S_vrr+5*oned2z*I_KINETIC_I2x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_Py_vrr;
    Double I_KINETIC_K2x4yz_Py_vrr = PBY*I_KINETIC_K2x4yz_S_vrr+4*oned2z*I_KINETIC_I2x3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_Py_vrr;
    Double I_KINETIC_K2x3y2z_Py_vrr = PBY*I_KINETIC_K2x3y2z_S_vrr+3*oned2z*I_KINETIC_I2x2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_Py_vrr;
    Double I_KINETIC_K2x2y3z_Py_vrr = PBY*I_KINETIC_K2x2y3z_S_vrr+2*oned2z*I_KINETIC_I2xy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_Py_vrr;
    Double I_KINETIC_K2xy4z_Py_vrr = PBY*I_KINETIC_K2xy4z_S_vrr+oned2z*I_KINETIC_I2x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_Py_vrr;
    Double I_KINETIC_K2x5z_Py_vrr = PBY*I_KINETIC_K2x5z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_Py_vrr;
    Double I_KINETIC_Kx6y_Py_vrr = PBY*I_KINETIC_Kx6y_S_vrr+6*oned2z*I_KINETIC_Ix5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_Py_vrr;
    Double I_KINETIC_Kx5yz_Py_vrr = PBY*I_KINETIC_Kx5yz_S_vrr+5*oned2z*I_KINETIC_Ix4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_Py_vrr;
    Double I_KINETIC_Kx4y2z_Py_vrr = PBY*I_KINETIC_Kx4y2z_S_vrr+4*oned2z*I_KINETIC_Ix3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_Py_vrr;
    Double I_KINETIC_Kx3y3z_Py_vrr = PBY*I_KINETIC_Kx3y3z_S_vrr+3*oned2z*I_KINETIC_Ix2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_Py_vrr;
    Double I_KINETIC_Kx2y4z_Py_vrr = PBY*I_KINETIC_Kx2y4z_S_vrr+2*oned2z*I_KINETIC_Ixy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_Py_vrr;
    Double I_KINETIC_Kxy5z_Py_vrr = PBY*I_KINETIC_Kxy5z_S_vrr+oned2z*I_KINETIC_Ix5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_Py_vrr;
    Double I_KINETIC_Kx6z_Py_vrr = PBY*I_KINETIC_Kx6z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_Py_vrr;
    Double I_KINETIC_K7y_Py_vrr = PBY*I_KINETIC_K7y_S_vrr+7*oned2z*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_K7y_Py_vrr;
    Double I_KINETIC_K6yz_Py_vrr = PBY*I_KINETIC_K6yz_S_vrr+6*oned2z*I_KINETIC_I5yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_Py_vrr;
    Double I_KINETIC_K5y2z_Py_vrr = PBY*I_KINETIC_K5y2z_S_vrr+5*oned2z*I_KINETIC_I4y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_Py_vrr;
    Double I_KINETIC_K4y3z_Py_vrr = PBY*I_KINETIC_K4y3z_S_vrr+4*oned2z*I_KINETIC_I3y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_Py_vrr;
    Double I_KINETIC_K3y4z_Py_vrr = PBY*I_KINETIC_K3y4z_S_vrr+3*oned2z*I_KINETIC_I2y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_Py_vrr;
    Double I_KINETIC_K2y5z_Py_vrr = PBY*I_KINETIC_K2y5z_S_vrr+2*oned2z*I_KINETIC_Iy5z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_Py_vrr;
    Double I_KINETIC_Ky6z_Py_vrr = PBY*I_KINETIC_Ky6z_S_vrr+oned2z*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_Py_vrr;
    Double I_KINETIC_K7z_Py_vrr = PBY*I_KINETIC_K7z_S_vrr+twoxi*I_TWOBODYOVERLAP_K7z_Py_vrr;
    Double I_KINETIC_K7x_Pz_vrr = PBZ*I_KINETIC_K7x_S_vrr+twoxi*I_TWOBODYOVERLAP_K7x_Pz_vrr;
    Double I_KINETIC_K6xy_Pz_vrr = PBZ*I_KINETIC_K6xy_S_vrr+twoxi*I_TWOBODYOVERLAP_K6xy_Pz_vrr;
    Double I_KINETIC_K6xz_Pz_vrr = PBZ*I_KINETIC_K6xz_S_vrr+oned2z*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_K6xz_Pz_vrr;
    Double I_KINETIC_K5x2y_Pz_vrr = PBZ*I_KINETIC_K5x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_K5x2y_Pz_vrr;
    Double I_KINETIC_K5xyz_Pz_vrr = PBZ*I_KINETIC_K5xyz_S_vrr+oned2z*I_KINETIC_I5xy_S_vrr+twoxi*I_TWOBODYOVERLAP_K5xyz_Pz_vrr;
    Double I_KINETIC_K5x2z_Pz_vrr = PBZ*I_KINETIC_K5x2z_S_vrr+2*oned2z*I_KINETIC_I5xz_S_vrr+twoxi*I_TWOBODYOVERLAP_K5x2z_Pz_vrr;
    Double I_KINETIC_K4x3y_Pz_vrr = PBZ*I_KINETIC_K4x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x3y_Pz_vrr;
    Double I_KINETIC_K4x2yz_Pz_vrr = PBZ*I_KINETIC_K4x2yz_S_vrr+oned2z*I_KINETIC_I4x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x2yz_Pz_vrr;
    Double I_KINETIC_K4xy2z_Pz_vrr = PBZ*I_KINETIC_K4xy2z_S_vrr+2*oned2z*I_KINETIC_I4xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_K4xy2z_Pz_vrr;
    Double I_KINETIC_K4x3z_Pz_vrr = PBZ*I_KINETIC_K4x3z_S_vrr+3*oned2z*I_KINETIC_I4x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K4x3z_Pz_vrr;
    Double I_KINETIC_K3x4y_Pz_vrr = PBZ*I_KINETIC_K3x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x4y_Pz_vrr;
    Double I_KINETIC_K3x3yz_Pz_vrr = PBZ*I_KINETIC_K3x3yz_S_vrr+oned2z*I_KINETIC_I3x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x3yz_Pz_vrr;
    Double I_KINETIC_K3x2y2z_Pz_vrr = PBZ*I_KINETIC_K3x2y2z_S_vrr+2*oned2z*I_KINETIC_I3x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x2y2z_Pz_vrr;
    Double I_KINETIC_K3xy3z_Pz_vrr = PBZ*I_KINETIC_K3xy3z_S_vrr+3*oned2z*I_KINETIC_I3xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3xy3z_Pz_vrr;
    Double I_KINETIC_K3x4z_Pz_vrr = PBZ*I_KINETIC_K3x4z_S_vrr+4*oned2z*I_KINETIC_I3x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3x4z_Pz_vrr;
    Double I_KINETIC_K2x5y_Pz_vrr = PBZ*I_KINETIC_K2x5y_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x5y_Pz_vrr;
    Double I_KINETIC_K2x4yz_Pz_vrr = PBZ*I_KINETIC_K2x4yz_S_vrr+oned2z*I_KINETIC_I2x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x4yz_Pz_vrr;
    Double I_KINETIC_K2x3y2z_Pz_vrr = PBZ*I_KINETIC_K2x3y2z_S_vrr+2*oned2z*I_KINETIC_I2x3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x3y2z_Pz_vrr;
    Double I_KINETIC_K2x2y3z_Pz_vrr = PBZ*I_KINETIC_K2x2y3z_S_vrr+3*oned2z*I_KINETIC_I2x2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x2y3z_Pz_vrr;
    Double I_KINETIC_K2xy4z_Pz_vrr = PBZ*I_KINETIC_K2xy4z_S_vrr+4*oned2z*I_KINETIC_I2xy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2xy4z_Pz_vrr;
    Double I_KINETIC_K2x5z_Pz_vrr = PBZ*I_KINETIC_K2x5z_S_vrr+5*oned2z*I_KINETIC_I2x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2x5z_Pz_vrr;
    Double I_KINETIC_Kx6y_Pz_vrr = PBZ*I_KINETIC_Kx6y_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx6y_Pz_vrr;
    Double I_KINETIC_Kx5yz_Pz_vrr = PBZ*I_KINETIC_Kx5yz_S_vrr+oned2z*I_KINETIC_Ix5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx5yz_Pz_vrr;
    Double I_KINETIC_Kx4y2z_Pz_vrr = PBZ*I_KINETIC_Kx4y2z_S_vrr+2*oned2z*I_KINETIC_Ix4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx4y2z_Pz_vrr;
    Double I_KINETIC_Kx3y3z_Pz_vrr = PBZ*I_KINETIC_Kx3y3z_S_vrr+3*oned2z*I_KINETIC_Ix3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx3y3z_Pz_vrr;
    Double I_KINETIC_Kx2y4z_Pz_vrr = PBZ*I_KINETIC_Kx2y4z_S_vrr+4*oned2z*I_KINETIC_Ix2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx2y4z_Pz_vrr;
    Double I_KINETIC_Kxy5z_Pz_vrr = PBZ*I_KINETIC_Kxy5z_S_vrr+5*oned2z*I_KINETIC_Ixy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kxy5z_Pz_vrr;
    Double I_KINETIC_Kx6z_Pz_vrr = PBZ*I_KINETIC_Kx6z_S_vrr+6*oned2z*I_KINETIC_Ix5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Kx6z_Pz_vrr;
    Double I_KINETIC_K7y_Pz_vrr = PBZ*I_KINETIC_K7y_S_vrr+twoxi*I_TWOBODYOVERLAP_K7y_Pz_vrr;
    Double I_KINETIC_K6yz_Pz_vrr = PBZ*I_KINETIC_K6yz_S_vrr+oned2z*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_K6yz_Pz_vrr;
    Double I_KINETIC_K5y2z_Pz_vrr = PBZ*I_KINETIC_K5y2z_S_vrr+2*oned2z*I_KINETIC_I5yz_S_vrr+twoxi*I_TWOBODYOVERLAP_K5y2z_Pz_vrr;
    Double I_KINETIC_K4y3z_Pz_vrr = PBZ*I_KINETIC_K4y3z_S_vrr+3*oned2z*I_KINETIC_I4y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_K4y3z_Pz_vrr;
    Double I_KINETIC_K3y4z_Pz_vrr = PBZ*I_KINETIC_K3y4z_S_vrr+4*oned2z*I_KINETIC_I3y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_K3y4z_Pz_vrr;
    Double I_KINETIC_K2y5z_Pz_vrr = PBZ*I_KINETIC_K2y5z_S_vrr+5*oned2z*I_KINETIC_I2y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_K2y5z_Pz_vrr;
    Double I_KINETIC_Ky6z_Pz_vrr = PBZ*I_KINETIC_Ky6z_S_vrr+6*oned2z*I_KINETIC_Iy5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ky6z_Pz_vrr;
    Double I_KINETIC_K7z_Pz_vrr = PBZ*I_KINETIC_K7z_S_vrr+7*oned2z*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_K7z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_K_P_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_K_P_aa_coefs = alpha*alpha;
    I_KINETIC_K7x_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K7x_Px_vrr;
    I_KINETIC_K6xy_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K6xy_Px_vrr;
    I_KINETIC_K6xz_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K6xz_Px_vrr;
    I_KINETIC_K5x2y_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5x2y_Px_vrr;
    I_KINETIC_K5xyz_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5xyz_Px_vrr;
    I_KINETIC_K5x2z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5x2z_Px_vrr;
    I_KINETIC_K4x3y_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4x3y_Px_vrr;
    I_KINETIC_K4x2yz_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4x2yz_Px_vrr;
    I_KINETIC_K4xy2z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4xy2z_Px_vrr;
    I_KINETIC_K4x3z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4x3z_Px_vrr;
    I_KINETIC_K3x4y_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x4y_Px_vrr;
    I_KINETIC_K3x3yz_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x3yz_Px_vrr;
    I_KINETIC_K3x2y2z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x2y2z_Px_vrr;
    I_KINETIC_K3xy3z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3xy3z_Px_vrr;
    I_KINETIC_K3x4z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x4z_Px_vrr;
    I_KINETIC_K2x5y_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x5y_Px_vrr;
    I_KINETIC_K2x4yz_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x4yz_Px_vrr;
    I_KINETIC_K2x3y2z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x3y2z_Px_vrr;
    I_KINETIC_K2x2y3z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x2y3z_Px_vrr;
    I_KINETIC_K2xy4z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2xy4z_Px_vrr;
    I_KINETIC_K2x5z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x5z_Px_vrr;
    I_KINETIC_Kx6y_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx6y_Px_vrr;
    I_KINETIC_Kx5yz_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx5yz_Px_vrr;
    I_KINETIC_Kx4y2z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx4y2z_Px_vrr;
    I_KINETIC_Kx3y3z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx3y3z_Px_vrr;
    I_KINETIC_Kx2y4z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx2y4z_Px_vrr;
    I_KINETIC_Kxy5z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kxy5z_Px_vrr;
    I_KINETIC_Kx6z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx6z_Px_vrr;
    I_KINETIC_K7y_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K7y_Px_vrr;
    I_KINETIC_K6yz_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K6yz_Px_vrr;
    I_KINETIC_K5y2z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5y2z_Px_vrr;
    I_KINETIC_K4y3z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4y3z_Px_vrr;
    I_KINETIC_K3y4z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3y4z_Px_vrr;
    I_KINETIC_K2y5z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2y5z_Px_vrr;
    I_KINETIC_Ky6z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Ky6z_Px_vrr;
    I_KINETIC_K7z_Px_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K7z_Px_vrr;
    I_KINETIC_K7x_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K7x_Py_vrr;
    I_KINETIC_K6xy_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K6xy_Py_vrr;
    I_KINETIC_K6xz_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K6xz_Py_vrr;
    I_KINETIC_K5x2y_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5x2y_Py_vrr;
    I_KINETIC_K5xyz_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5xyz_Py_vrr;
    I_KINETIC_K5x2z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5x2z_Py_vrr;
    I_KINETIC_K4x3y_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4x3y_Py_vrr;
    I_KINETIC_K4x2yz_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4x2yz_Py_vrr;
    I_KINETIC_K4xy2z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4xy2z_Py_vrr;
    I_KINETIC_K4x3z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4x3z_Py_vrr;
    I_KINETIC_K3x4y_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x4y_Py_vrr;
    I_KINETIC_K3x3yz_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x3yz_Py_vrr;
    I_KINETIC_K3x2y2z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x2y2z_Py_vrr;
    I_KINETIC_K3xy3z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3xy3z_Py_vrr;
    I_KINETIC_K3x4z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x4z_Py_vrr;
    I_KINETIC_K2x5y_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x5y_Py_vrr;
    I_KINETIC_K2x4yz_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x4yz_Py_vrr;
    I_KINETIC_K2x3y2z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x3y2z_Py_vrr;
    I_KINETIC_K2x2y3z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x2y3z_Py_vrr;
    I_KINETIC_K2xy4z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2xy4z_Py_vrr;
    I_KINETIC_K2x5z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x5z_Py_vrr;
    I_KINETIC_Kx6y_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx6y_Py_vrr;
    I_KINETIC_Kx5yz_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx5yz_Py_vrr;
    I_KINETIC_Kx4y2z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx4y2z_Py_vrr;
    I_KINETIC_Kx3y3z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx3y3z_Py_vrr;
    I_KINETIC_Kx2y4z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx2y4z_Py_vrr;
    I_KINETIC_Kxy5z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kxy5z_Py_vrr;
    I_KINETIC_Kx6z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx6z_Py_vrr;
    I_KINETIC_K7y_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K7y_Py_vrr;
    I_KINETIC_K6yz_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K6yz_Py_vrr;
    I_KINETIC_K5y2z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5y2z_Py_vrr;
    I_KINETIC_K4y3z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4y3z_Py_vrr;
    I_KINETIC_K3y4z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3y4z_Py_vrr;
    I_KINETIC_K2y5z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2y5z_Py_vrr;
    I_KINETIC_Ky6z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Ky6z_Py_vrr;
    I_KINETIC_K7z_Py_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K7z_Py_vrr;
    I_KINETIC_K7x_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K7x_Pz_vrr;
    I_KINETIC_K6xy_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K6xy_Pz_vrr;
    I_KINETIC_K6xz_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K6xz_Pz_vrr;
    I_KINETIC_K5x2y_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5x2y_Pz_vrr;
    I_KINETIC_K5xyz_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5xyz_Pz_vrr;
    I_KINETIC_K5x2z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5x2z_Pz_vrr;
    I_KINETIC_K4x3y_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4x3y_Pz_vrr;
    I_KINETIC_K4x2yz_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4x2yz_Pz_vrr;
    I_KINETIC_K4xy2z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4xy2z_Pz_vrr;
    I_KINETIC_K4x3z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4x3z_Pz_vrr;
    I_KINETIC_K3x4y_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x4y_Pz_vrr;
    I_KINETIC_K3x3yz_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x3yz_Pz_vrr;
    I_KINETIC_K3x2y2z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x2y2z_Pz_vrr;
    I_KINETIC_K3xy3z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3xy3z_Pz_vrr;
    I_KINETIC_K3x4z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3x4z_Pz_vrr;
    I_KINETIC_K2x5y_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x5y_Pz_vrr;
    I_KINETIC_K2x4yz_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x4yz_Pz_vrr;
    I_KINETIC_K2x3y2z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x3y2z_Pz_vrr;
    I_KINETIC_K2x2y3z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x2y3z_Pz_vrr;
    I_KINETIC_K2xy4z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2xy4z_Pz_vrr;
    I_KINETIC_K2x5z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2x5z_Pz_vrr;
    I_KINETIC_Kx6y_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx6y_Pz_vrr;
    I_KINETIC_Kx5yz_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx5yz_Pz_vrr;
    I_KINETIC_Kx4y2z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx4y2z_Pz_vrr;
    I_KINETIC_Kx3y3z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx3y3z_Pz_vrr;
    I_KINETIC_Kx2y4z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx2y4z_Pz_vrr;
    I_KINETIC_Kxy5z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kxy5z_Pz_vrr;
    I_KINETIC_Kx6z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Kx6z_Pz_vrr;
    I_KINETIC_K7y_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K7y_Pz_vrr;
    I_KINETIC_K6yz_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K6yz_Pz_vrr;
    I_KINETIC_K5y2z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K5y2z_Pz_vrr;
    I_KINETIC_K4y3z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K4y3z_Pz_vrr;
    I_KINETIC_K3y4z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K3y4z_Pz_vrr;
    I_KINETIC_K2y5z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K2y5z_Pz_vrr;
    I_KINETIC_Ky6z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_Ky6z_Pz_vrr;
    I_KINETIC_K7z_Pz_aa += SQ_KINETIC_K_P_aa_coefs*I_KINETIC_K7z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_P_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_H_P_a_coefs = alpha;
    I_KINETIC_H5x_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H5x_Px_vrr;
    I_KINETIC_H4xy_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H4xy_Px_vrr;
    I_KINETIC_H4xz_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H4xz_Px_vrr;
    I_KINETIC_H3x2y_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3x2y_Px_vrr;
    I_KINETIC_H3xyz_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3xyz_Px_vrr;
    I_KINETIC_H3x2z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3x2z_Px_vrr;
    I_KINETIC_H2x3y_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2x3y_Px_vrr;
    I_KINETIC_H2x2yz_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2x2yz_Px_vrr;
    I_KINETIC_H2xy2z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2xy2z_Px_vrr;
    I_KINETIC_H2x3z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2x3z_Px_vrr;
    I_KINETIC_Hx4y_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx4y_Px_vrr;
    I_KINETIC_Hx3yz_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx3yz_Px_vrr;
    I_KINETIC_Hx2y2z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx2y2z_Px_vrr;
    I_KINETIC_Hxy3z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hxy3z_Px_vrr;
    I_KINETIC_Hx4z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx4z_Px_vrr;
    I_KINETIC_H5y_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H5y_Px_vrr;
    I_KINETIC_H4yz_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H4yz_Px_vrr;
    I_KINETIC_H3y2z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3y2z_Px_vrr;
    I_KINETIC_H2y3z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2y3z_Px_vrr;
    I_KINETIC_Hy4z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hy4z_Px_vrr;
    I_KINETIC_H5z_Px_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H5z_Px_vrr;
    I_KINETIC_H5x_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H5x_Py_vrr;
    I_KINETIC_H4xy_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H4xy_Py_vrr;
    I_KINETIC_H4xz_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H4xz_Py_vrr;
    I_KINETIC_H3x2y_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3x2y_Py_vrr;
    I_KINETIC_H3xyz_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3xyz_Py_vrr;
    I_KINETIC_H3x2z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3x2z_Py_vrr;
    I_KINETIC_H2x3y_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2x3y_Py_vrr;
    I_KINETIC_H2x2yz_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2x2yz_Py_vrr;
    I_KINETIC_H2xy2z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2xy2z_Py_vrr;
    I_KINETIC_H2x3z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2x3z_Py_vrr;
    I_KINETIC_Hx4y_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx4y_Py_vrr;
    I_KINETIC_Hx3yz_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx3yz_Py_vrr;
    I_KINETIC_Hx2y2z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx2y2z_Py_vrr;
    I_KINETIC_Hxy3z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hxy3z_Py_vrr;
    I_KINETIC_Hx4z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx4z_Py_vrr;
    I_KINETIC_H5y_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H5y_Py_vrr;
    I_KINETIC_H4yz_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H4yz_Py_vrr;
    I_KINETIC_H3y2z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3y2z_Py_vrr;
    I_KINETIC_H2y3z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2y3z_Py_vrr;
    I_KINETIC_Hy4z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hy4z_Py_vrr;
    I_KINETIC_H5z_Py_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H5z_Py_vrr;
    I_KINETIC_H5x_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H5x_Pz_vrr;
    I_KINETIC_H4xy_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H4xy_Pz_vrr;
    I_KINETIC_H4xz_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H4xz_Pz_vrr;
    I_KINETIC_H3x2y_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3x2y_Pz_vrr;
    I_KINETIC_H3xyz_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3xyz_Pz_vrr;
    I_KINETIC_H3x2z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3x2z_Pz_vrr;
    I_KINETIC_H2x3y_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2x3y_Pz_vrr;
    I_KINETIC_H2x2yz_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2x2yz_Pz_vrr;
    I_KINETIC_H2xy2z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2xy2z_Pz_vrr;
    I_KINETIC_H2x3z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2x3z_Pz_vrr;
    I_KINETIC_Hx4y_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx4y_Pz_vrr;
    I_KINETIC_Hx3yz_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx3yz_Pz_vrr;
    I_KINETIC_Hx2y2z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx2y2z_Pz_vrr;
    I_KINETIC_Hxy3z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hxy3z_Pz_vrr;
    I_KINETIC_Hx4z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hx4z_Pz_vrr;
    I_KINETIC_H5y_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H5y_Pz_vrr;
    I_KINETIC_H4yz_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H4yz_Pz_vrr;
    I_KINETIC_H3y2z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H3y2z_Pz_vrr;
    I_KINETIC_H2y3z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H2y3z_Pz_vrr;
    I_KINETIC_Hy4z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_Hy4z_Pz_vrr;
    I_KINETIC_H5z_Pz_a += SQ_KINETIC_H_P_a_coefs*I_KINETIC_H5z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_P
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_F3x_Px += I_KINETIC_F3x_Px_vrr;
    I_KINETIC_F2xy_Px += I_KINETIC_F2xy_Px_vrr;
    I_KINETIC_F2xz_Px += I_KINETIC_F2xz_Px_vrr;
    I_KINETIC_Fx2y_Px += I_KINETIC_Fx2y_Px_vrr;
    I_KINETIC_Fxyz_Px += I_KINETIC_Fxyz_Px_vrr;
    I_KINETIC_Fx2z_Px += I_KINETIC_Fx2z_Px_vrr;
    I_KINETIC_F3y_Px += I_KINETIC_F3y_Px_vrr;
    I_KINETIC_F2yz_Px += I_KINETIC_F2yz_Px_vrr;
    I_KINETIC_Fy2z_Px += I_KINETIC_Fy2z_Px_vrr;
    I_KINETIC_F3z_Px += I_KINETIC_F3z_Px_vrr;
    I_KINETIC_F3x_Py += I_KINETIC_F3x_Py_vrr;
    I_KINETIC_F2xy_Py += I_KINETIC_F2xy_Py_vrr;
    I_KINETIC_F2xz_Py += I_KINETIC_F2xz_Py_vrr;
    I_KINETIC_Fx2y_Py += I_KINETIC_Fx2y_Py_vrr;
    I_KINETIC_Fxyz_Py += I_KINETIC_Fxyz_Py_vrr;
    I_KINETIC_Fx2z_Py += I_KINETIC_Fx2z_Py_vrr;
    I_KINETIC_F3y_Py += I_KINETIC_F3y_Py_vrr;
    I_KINETIC_F2yz_Py += I_KINETIC_F2yz_Py_vrr;
    I_KINETIC_Fy2z_Py += I_KINETIC_Fy2z_Py_vrr;
    I_KINETIC_F3z_Py += I_KINETIC_F3z_Py_vrr;
    I_KINETIC_F3x_Pz += I_KINETIC_F3x_Pz_vrr;
    I_KINETIC_F2xy_Pz += I_KINETIC_F2xy_Pz_vrr;
    I_KINETIC_F2xz_Pz += I_KINETIC_F2xz_Pz_vrr;
    I_KINETIC_Fx2y_Pz += I_KINETIC_Fx2y_Pz_vrr;
    I_KINETIC_Fxyz_Pz += I_KINETIC_Fxyz_Pz_vrr;
    I_KINETIC_Fx2z_Pz += I_KINETIC_Fx2z_Pz_vrr;
    I_KINETIC_F3y_Pz += I_KINETIC_F3y_Pz_vrr;
    I_KINETIC_F2yz_Pz += I_KINETIC_F2yz_Pz_vrr;
    I_KINETIC_Fy2z_Pz += I_KINETIC_Fy2z_Pz_vrr;
    I_KINETIC_F3z_Pz += I_KINETIC_F3z_Pz_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_P_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_P_aa
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_F_P
   ************************************************************/
  abcd[0] = 4.0E0*I_KINETIC_K7x_Px_aa-2.0E0*5*I_KINETIC_H5x_Px_a-2.0E0*6*I_KINETIC_H5x_Px_a+5*4*I_KINETIC_F3x_Px;
  abcd[1] = 4.0E0*I_KINETIC_K6xy_Px_aa-2.0E0*4*I_KINETIC_H4xy_Px_a-2.0E0*5*I_KINETIC_H4xy_Px_a+4*3*I_KINETIC_F2xy_Px;
  abcd[2] = 4.0E0*I_KINETIC_K6xz_Px_aa-2.0E0*4*I_KINETIC_H4xz_Px_a-2.0E0*5*I_KINETIC_H4xz_Px_a+4*3*I_KINETIC_F2xz_Px;
  abcd[3] = 4.0E0*I_KINETIC_K5x2y_Px_aa-2.0E0*3*I_KINETIC_H3x2y_Px_a-2.0E0*4*I_KINETIC_H3x2y_Px_a+3*2*I_KINETIC_Fx2y_Px;
  abcd[4] = 4.0E0*I_KINETIC_K5xyz_Px_aa-2.0E0*3*I_KINETIC_H3xyz_Px_a-2.0E0*4*I_KINETIC_H3xyz_Px_a+3*2*I_KINETIC_Fxyz_Px;
  abcd[5] = 4.0E0*I_KINETIC_K5x2z_Px_aa-2.0E0*3*I_KINETIC_H3x2z_Px_a-2.0E0*4*I_KINETIC_H3x2z_Px_a+3*2*I_KINETIC_Fx2z_Px;
  abcd[6] = 4.0E0*I_KINETIC_K4x3y_Px_aa-2.0E0*2*I_KINETIC_H2x3y_Px_a-2.0E0*3*I_KINETIC_H2x3y_Px_a+2*1*I_KINETIC_F3y_Px;
  abcd[7] = 4.0E0*I_KINETIC_K4x2yz_Px_aa-2.0E0*2*I_KINETIC_H2x2yz_Px_a-2.0E0*3*I_KINETIC_H2x2yz_Px_a+2*1*I_KINETIC_F2yz_Px;
  abcd[8] = 4.0E0*I_KINETIC_K4xy2z_Px_aa-2.0E0*2*I_KINETIC_H2xy2z_Px_a-2.0E0*3*I_KINETIC_H2xy2z_Px_a+2*1*I_KINETIC_Fy2z_Px;
  abcd[9] = 4.0E0*I_KINETIC_K4x3z_Px_aa-2.0E0*2*I_KINETIC_H2x3z_Px_a-2.0E0*3*I_KINETIC_H2x3z_Px_a+2*1*I_KINETIC_F3z_Px;
  abcd[10] = 4.0E0*I_KINETIC_K3x4y_Px_aa-2.0E0*1*I_KINETIC_Hx4y_Px_a-2.0E0*2*I_KINETIC_Hx4y_Px_a;
  abcd[11] = 4.0E0*I_KINETIC_K3x3yz_Px_aa-2.0E0*1*I_KINETIC_Hx3yz_Px_a-2.0E0*2*I_KINETIC_Hx3yz_Px_a;
  abcd[12] = 4.0E0*I_KINETIC_K3x2y2z_Px_aa-2.0E0*1*I_KINETIC_Hx2y2z_Px_a-2.0E0*2*I_KINETIC_Hx2y2z_Px_a;
  abcd[13] = 4.0E0*I_KINETIC_K3xy3z_Px_aa-2.0E0*1*I_KINETIC_Hxy3z_Px_a-2.0E0*2*I_KINETIC_Hxy3z_Px_a;
  abcd[14] = 4.0E0*I_KINETIC_K3x4z_Px_aa-2.0E0*1*I_KINETIC_Hx4z_Px_a-2.0E0*2*I_KINETIC_Hx4z_Px_a;
  abcd[15] = 4.0E0*I_KINETIC_K2x5y_Px_aa-2.0E0*1*I_KINETIC_H5y_Px_a;
  abcd[16] = 4.0E0*I_KINETIC_K2x4yz_Px_aa-2.0E0*1*I_KINETIC_H4yz_Px_a;
  abcd[17] = 4.0E0*I_KINETIC_K2x3y2z_Px_aa-2.0E0*1*I_KINETIC_H3y2z_Px_a;
  abcd[18] = 4.0E0*I_KINETIC_K2x2y3z_Px_aa-2.0E0*1*I_KINETIC_H2y3z_Px_a;
  abcd[19] = 4.0E0*I_KINETIC_K2xy4z_Px_aa-2.0E0*1*I_KINETIC_Hy4z_Px_a;
  abcd[20] = 4.0E0*I_KINETIC_K2x5z_Px_aa-2.0E0*1*I_KINETIC_H5z_Px_a;
  abcd[21] = 4.0E0*I_KINETIC_K7x_Py_aa-2.0E0*5*I_KINETIC_H5x_Py_a-2.0E0*6*I_KINETIC_H5x_Py_a+5*4*I_KINETIC_F3x_Py;
  abcd[22] = 4.0E0*I_KINETIC_K6xy_Py_aa-2.0E0*4*I_KINETIC_H4xy_Py_a-2.0E0*5*I_KINETIC_H4xy_Py_a+4*3*I_KINETIC_F2xy_Py;
  abcd[23] = 4.0E0*I_KINETIC_K6xz_Py_aa-2.0E0*4*I_KINETIC_H4xz_Py_a-2.0E0*5*I_KINETIC_H4xz_Py_a+4*3*I_KINETIC_F2xz_Py;
  abcd[24] = 4.0E0*I_KINETIC_K5x2y_Py_aa-2.0E0*3*I_KINETIC_H3x2y_Py_a-2.0E0*4*I_KINETIC_H3x2y_Py_a+3*2*I_KINETIC_Fx2y_Py;
  abcd[25] = 4.0E0*I_KINETIC_K5xyz_Py_aa-2.0E0*3*I_KINETIC_H3xyz_Py_a-2.0E0*4*I_KINETIC_H3xyz_Py_a+3*2*I_KINETIC_Fxyz_Py;
  abcd[26] = 4.0E0*I_KINETIC_K5x2z_Py_aa-2.0E0*3*I_KINETIC_H3x2z_Py_a-2.0E0*4*I_KINETIC_H3x2z_Py_a+3*2*I_KINETIC_Fx2z_Py;
  abcd[27] = 4.0E0*I_KINETIC_K4x3y_Py_aa-2.0E0*2*I_KINETIC_H2x3y_Py_a-2.0E0*3*I_KINETIC_H2x3y_Py_a+2*1*I_KINETIC_F3y_Py;
  abcd[28] = 4.0E0*I_KINETIC_K4x2yz_Py_aa-2.0E0*2*I_KINETIC_H2x2yz_Py_a-2.0E0*3*I_KINETIC_H2x2yz_Py_a+2*1*I_KINETIC_F2yz_Py;
  abcd[29] = 4.0E0*I_KINETIC_K4xy2z_Py_aa-2.0E0*2*I_KINETIC_H2xy2z_Py_a-2.0E0*3*I_KINETIC_H2xy2z_Py_a+2*1*I_KINETIC_Fy2z_Py;
  abcd[30] = 4.0E0*I_KINETIC_K4x3z_Py_aa-2.0E0*2*I_KINETIC_H2x3z_Py_a-2.0E0*3*I_KINETIC_H2x3z_Py_a+2*1*I_KINETIC_F3z_Py;
  abcd[31] = 4.0E0*I_KINETIC_K3x4y_Py_aa-2.0E0*1*I_KINETIC_Hx4y_Py_a-2.0E0*2*I_KINETIC_Hx4y_Py_a;
  abcd[32] = 4.0E0*I_KINETIC_K3x3yz_Py_aa-2.0E0*1*I_KINETIC_Hx3yz_Py_a-2.0E0*2*I_KINETIC_Hx3yz_Py_a;
  abcd[33] = 4.0E0*I_KINETIC_K3x2y2z_Py_aa-2.0E0*1*I_KINETIC_Hx2y2z_Py_a-2.0E0*2*I_KINETIC_Hx2y2z_Py_a;
  abcd[34] = 4.0E0*I_KINETIC_K3xy3z_Py_aa-2.0E0*1*I_KINETIC_Hxy3z_Py_a-2.0E0*2*I_KINETIC_Hxy3z_Py_a;
  abcd[35] = 4.0E0*I_KINETIC_K3x4z_Py_aa-2.0E0*1*I_KINETIC_Hx4z_Py_a-2.0E0*2*I_KINETIC_Hx4z_Py_a;
  abcd[36] = 4.0E0*I_KINETIC_K2x5y_Py_aa-2.0E0*1*I_KINETIC_H5y_Py_a;
  abcd[37] = 4.0E0*I_KINETIC_K2x4yz_Py_aa-2.0E0*1*I_KINETIC_H4yz_Py_a;
  abcd[38] = 4.0E0*I_KINETIC_K2x3y2z_Py_aa-2.0E0*1*I_KINETIC_H3y2z_Py_a;
  abcd[39] = 4.0E0*I_KINETIC_K2x2y3z_Py_aa-2.0E0*1*I_KINETIC_H2y3z_Py_a;
  abcd[40] = 4.0E0*I_KINETIC_K2xy4z_Py_aa-2.0E0*1*I_KINETIC_Hy4z_Py_a;
  abcd[41] = 4.0E0*I_KINETIC_K2x5z_Py_aa-2.0E0*1*I_KINETIC_H5z_Py_a;
  abcd[42] = 4.0E0*I_KINETIC_K7x_Pz_aa-2.0E0*5*I_KINETIC_H5x_Pz_a-2.0E0*6*I_KINETIC_H5x_Pz_a+5*4*I_KINETIC_F3x_Pz;
  abcd[43] = 4.0E0*I_KINETIC_K6xy_Pz_aa-2.0E0*4*I_KINETIC_H4xy_Pz_a-2.0E0*5*I_KINETIC_H4xy_Pz_a+4*3*I_KINETIC_F2xy_Pz;
  abcd[44] = 4.0E0*I_KINETIC_K6xz_Pz_aa-2.0E0*4*I_KINETIC_H4xz_Pz_a-2.0E0*5*I_KINETIC_H4xz_Pz_a+4*3*I_KINETIC_F2xz_Pz;
  abcd[45] = 4.0E0*I_KINETIC_K5x2y_Pz_aa-2.0E0*3*I_KINETIC_H3x2y_Pz_a-2.0E0*4*I_KINETIC_H3x2y_Pz_a+3*2*I_KINETIC_Fx2y_Pz;
  abcd[46] = 4.0E0*I_KINETIC_K5xyz_Pz_aa-2.0E0*3*I_KINETIC_H3xyz_Pz_a-2.0E0*4*I_KINETIC_H3xyz_Pz_a+3*2*I_KINETIC_Fxyz_Pz;
  abcd[47] = 4.0E0*I_KINETIC_K5x2z_Pz_aa-2.0E0*3*I_KINETIC_H3x2z_Pz_a-2.0E0*4*I_KINETIC_H3x2z_Pz_a+3*2*I_KINETIC_Fx2z_Pz;
  abcd[48] = 4.0E0*I_KINETIC_K4x3y_Pz_aa-2.0E0*2*I_KINETIC_H2x3y_Pz_a-2.0E0*3*I_KINETIC_H2x3y_Pz_a+2*1*I_KINETIC_F3y_Pz;
  abcd[49] = 4.0E0*I_KINETIC_K4x2yz_Pz_aa-2.0E0*2*I_KINETIC_H2x2yz_Pz_a-2.0E0*3*I_KINETIC_H2x2yz_Pz_a+2*1*I_KINETIC_F2yz_Pz;
  abcd[50] = 4.0E0*I_KINETIC_K4xy2z_Pz_aa-2.0E0*2*I_KINETIC_H2xy2z_Pz_a-2.0E0*3*I_KINETIC_H2xy2z_Pz_a+2*1*I_KINETIC_Fy2z_Pz;
  abcd[51] = 4.0E0*I_KINETIC_K4x3z_Pz_aa-2.0E0*2*I_KINETIC_H2x3z_Pz_a-2.0E0*3*I_KINETIC_H2x3z_Pz_a+2*1*I_KINETIC_F3z_Pz;
  abcd[52] = 4.0E0*I_KINETIC_K3x4y_Pz_aa-2.0E0*1*I_KINETIC_Hx4y_Pz_a-2.0E0*2*I_KINETIC_Hx4y_Pz_a;
  abcd[53] = 4.0E0*I_KINETIC_K3x3yz_Pz_aa-2.0E0*1*I_KINETIC_Hx3yz_Pz_a-2.0E0*2*I_KINETIC_Hx3yz_Pz_a;
  abcd[54] = 4.0E0*I_KINETIC_K3x2y2z_Pz_aa-2.0E0*1*I_KINETIC_Hx2y2z_Pz_a-2.0E0*2*I_KINETIC_Hx2y2z_Pz_a;
  abcd[55] = 4.0E0*I_KINETIC_K3xy3z_Pz_aa-2.0E0*1*I_KINETIC_Hxy3z_Pz_a-2.0E0*2*I_KINETIC_Hxy3z_Pz_a;
  abcd[56] = 4.0E0*I_KINETIC_K3x4z_Pz_aa-2.0E0*1*I_KINETIC_Hx4z_Pz_a-2.0E0*2*I_KINETIC_Hx4z_Pz_a;
  abcd[57] = 4.0E0*I_KINETIC_K2x5y_Pz_aa-2.0E0*1*I_KINETIC_H5y_Pz_a;
  abcd[58] = 4.0E0*I_KINETIC_K2x4yz_Pz_aa-2.0E0*1*I_KINETIC_H4yz_Pz_a;
  abcd[59] = 4.0E0*I_KINETIC_K2x3y2z_Pz_aa-2.0E0*1*I_KINETIC_H3y2z_Pz_a;
  abcd[60] = 4.0E0*I_KINETIC_K2x2y3z_Pz_aa-2.0E0*1*I_KINETIC_H2y3z_Pz_a;
  abcd[61] = 4.0E0*I_KINETIC_K2xy4z_Pz_aa-2.0E0*1*I_KINETIC_Hy4z_Pz_a;
  abcd[62] = 4.0E0*I_KINETIC_K2x5z_Pz_aa-2.0E0*1*I_KINETIC_H5z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_P_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_P_aa
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_F_P
   ************************************************************/
  abcd[63] = 4.0E0*I_KINETIC_K6xy_Px_aa-2.0E0*5*I_KINETIC_H4xy_Px_a;
  abcd[64] = 4.0E0*I_KINETIC_K5x2y_Px_aa-2.0E0*1*I_KINETIC_H5x_Px_a-2.0E0*4*I_KINETIC_H3x2y_Px_a+4*1*I_KINETIC_F3x_Px;
  abcd[65] = 4.0E0*I_KINETIC_K5xyz_Px_aa-2.0E0*4*I_KINETIC_H3xyz_Px_a;
  abcd[66] = 4.0E0*I_KINETIC_K4x3y_Px_aa-2.0E0*2*I_KINETIC_H4xy_Px_a-2.0E0*3*I_KINETIC_H2x3y_Px_a+3*2*I_KINETIC_F2xy_Px;
  abcd[67] = 4.0E0*I_KINETIC_K4x2yz_Px_aa-2.0E0*1*I_KINETIC_H4xz_Px_a-2.0E0*3*I_KINETIC_H2x2yz_Px_a+3*1*I_KINETIC_F2xz_Px;
  abcd[68] = 4.0E0*I_KINETIC_K4xy2z_Px_aa-2.0E0*3*I_KINETIC_H2xy2z_Px_a;
  abcd[69] = 4.0E0*I_KINETIC_K3x4y_Px_aa-2.0E0*3*I_KINETIC_H3x2y_Px_a-2.0E0*2*I_KINETIC_Hx4y_Px_a+2*3*I_KINETIC_Fx2y_Px;
  abcd[70] = 4.0E0*I_KINETIC_K3x3yz_Px_aa-2.0E0*2*I_KINETIC_H3xyz_Px_a-2.0E0*2*I_KINETIC_Hx3yz_Px_a+2*2*I_KINETIC_Fxyz_Px;
  abcd[71] = 4.0E0*I_KINETIC_K3x2y2z_Px_aa-2.0E0*1*I_KINETIC_H3x2z_Px_a-2.0E0*2*I_KINETIC_Hx2y2z_Px_a+2*1*I_KINETIC_Fx2z_Px;
  abcd[72] = 4.0E0*I_KINETIC_K3xy3z_Px_aa-2.0E0*2*I_KINETIC_Hxy3z_Px_a;
  abcd[73] = 4.0E0*I_KINETIC_K2x5y_Px_aa-2.0E0*4*I_KINETIC_H2x3y_Px_a-2.0E0*1*I_KINETIC_H5y_Px_a+4*I_KINETIC_F3y_Px;
  abcd[74] = 4.0E0*I_KINETIC_K2x4yz_Px_aa-2.0E0*3*I_KINETIC_H2x2yz_Px_a-2.0E0*1*I_KINETIC_H4yz_Px_a+3*I_KINETIC_F2yz_Px;
  abcd[75] = 4.0E0*I_KINETIC_K2x3y2z_Px_aa-2.0E0*2*I_KINETIC_H2xy2z_Px_a-2.0E0*1*I_KINETIC_H3y2z_Px_a+2*I_KINETIC_Fy2z_Px;
  abcd[76] = 4.0E0*I_KINETIC_K2x2y3z_Px_aa-2.0E0*1*I_KINETIC_H2x3z_Px_a-2.0E0*1*I_KINETIC_H2y3z_Px_a+1*I_KINETIC_F3z_Px;
  abcd[77] = 4.0E0*I_KINETIC_K2xy4z_Px_aa-2.0E0*1*I_KINETIC_Hy4z_Px_a;
  abcd[78] = 4.0E0*I_KINETIC_Kx6y_Px_aa-2.0E0*5*I_KINETIC_Hx4y_Px_a;
  abcd[79] = 4.0E0*I_KINETIC_Kx5yz_Px_aa-2.0E0*4*I_KINETIC_Hx3yz_Px_a;
  abcd[80] = 4.0E0*I_KINETIC_Kx4y2z_Px_aa-2.0E0*3*I_KINETIC_Hx2y2z_Px_a;
  abcd[81] = 4.0E0*I_KINETIC_Kx3y3z_Px_aa-2.0E0*2*I_KINETIC_Hxy3z_Px_a;
  abcd[82] = 4.0E0*I_KINETIC_Kx2y4z_Px_aa-2.0E0*1*I_KINETIC_Hx4z_Px_a;
  abcd[83] = 4.0E0*I_KINETIC_Kxy5z_Px_aa;
  abcd[84] = 4.0E0*I_KINETIC_K6xy_Py_aa-2.0E0*5*I_KINETIC_H4xy_Py_a;
  abcd[85] = 4.0E0*I_KINETIC_K5x2y_Py_aa-2.0E0*1*I_KINETIC_H5x_Py_a-2.0E0*4*I_KINETIC_H3x2y_Py_a+4*1*I_KINETIC_F3x_Py;
  abcd[86] = 4.0E0*I_KINETIC_K5xyz_Py_aa-2.0E0*4*I_KINETIC_H3xyz_Py_a;
  abcd[87] = 4.0E0*I_KINETIC_K4x3y_Py_aa-2.0E0*2*I_KINETIC_H4xy_Py_a-2.0E0*3*I_KINETIC_H2x3y_Py_a+3*2*I_KINETIC_F2xy_Py;
  abcd[88] = 4.0E0*I_KINETIC_K4x2yz_Py_aa-2.0E0*1*I_KINETIC_H4xz_Py_a-2.0E0*3*I_KINETIC_H2x2yz_Py_a+3*1*I_KINETIC_F2xz_Py;
  abcd[89] = 4.0E0*I_KINETIC_K4xy2z_Py_aa-2.0E0*3*I_KINETIC_H2xy2z_Py_a;
  abcd[90] = 4.0E0*I_KINETIC_K3x4y_Py_aa-2.0E0*3*I_KINETIC_H3x2y_Py_a-2.0E0*2*I_KINETIC_Hx4y_Py_a+2*3*I_KINETIC_Fx2y_Py;
  abcd[91] = 4.0E0*I_KINETIC_K3x3yz_Py_aa-2.0E0*2*I_KINETIC_H3xyz_Py_a-2.0E0*2*I_KINETIC_Hx3yz_Py_a+2*2*I_KINETIC_Fxyz_Py;
  abcd[92] = 4.0E0*I_KINETIC_K3x2y2z_Py_aa-2.0E0*1*I_KINETIC_H3x2z_Py_a-2.0E0*2*I_KINETIC_Hx2y2z_Py_a+2*1*I_KINETIC_Fx2z_Py;
  abcd[93] = 4.0E0*I_KINETIC_K3xy3z_Py_aa-2.0E0*2*I_KINETIC_Hxy3z_Py_a;
  abcd[94] = 4.0E0*I_KINETIC_K2x5y_Py_aa-2.0E0*4*I_KINETIC_H2x3y_Py_a-2.0E0*1*I_KINETIC_H5y_Py_a+4*I_KINETIC_F3y_Py;
  abcd[95] = 4.0E0*I_KINETIC_K2x4yz_Py_aa-2.0E0*3*I_KINETIC_H2x2yz_Py_a-2.0E0*1*I_KINETIC_H4yz_Py_a+3*I_KINETIC_F2yz_Py;
  abcd[96] = 4.0E0*I_KINETIC_K2x3y2z_Py_aa-2.0E0*2*I_KINETIC_H2xy2z_Py_a-2.0E0*1*I_KINETIC_H3y2z_Py_a+2*I_KINETIC_Fy2z_Py;
  abcd[97] = 4.0E0*I_KINETIC_K2x2y3z_Py_aa-2.0E0*1*I_KINETIC_H2x3z_Py_a-2.0E0*1*I_KINETIC_H2y3z_Py_a+1*I_KINETIC_F3z_Py;
  abcd[98] = 4.0E0*I_KINETIC_K2xy4z_Py_aa-2.0E0*1*I_KINETIC_Hy4z_Py_a;
  abcd[99] = 4.0E0*I_KINETIC_Kx6y_Py_aa-2.0E0*5*I_KINETIC_Hx4y_Py_a;
  abcd[100] = 4.0E0*I_KINETIC_Kx5yz_Py_aa-2.0E0*4*I_KINETIC_Hx3yz_Py_a;
  abcd[101] = 4.0E0*I_KINETIC_Kx4y2z_Py_aa-2.0E0*3*I_KINETIC_Hx2y2z_Py_a;
  abcd[102] = 4.0E0*I_KINETIC_Kx3y3z_Py_aa-2.0E0*2*I_KINETIC_Hxy3z_Py_a;
  abcd[103] = 4.0E0*I_KINETIC_Kx2y4z_Py_aa-2.0E0*1*I_KINETIC_Hx4z_Py_a;
  abcd[104] = 4.0E0*I_KINETIC_Kxy5z_Py_aa;
  abcd[105] = 4.0E0*I_KINETIC_K6xy_Pz_aa-2.0E0*5*I_KINETIC_H4xy_Pz_a;
  abcd[106] = 4.0E0*I_KINETIC_K5x2y_Pz_aa-2.0E0*1*I_KINETIC_H5x_Pz_a-2.0E0*4*I_KINETIC_H3x2y_Pz_a+4*1*I_KINETIC_F3x_Pz;
  abcd[107] = 4.0E0*I_KINETIC_K5xyz_Pz_aa-2.0E0*4*I_KINETIC_H3xyz_Pz_a;
  abcd[108] = 4.0E0*I_KINETIC_K4x3y_Pz_aa-2.0E0*2*I_KINETIC_H4xy_Pz_a-2.0E0*3*I_KINETIC_H2x3y_Pz_a+3*2*I_KINETIC_F2xy_Pz;
  abcd[109] = 4.0E0*I_KINETIC_K4x2yz_Pz_aa-2.0E0*1*I_KINETIC_H4xz_Pz_a-2.0E0*3*I_KINETIC_H2x2yz_Pz_a+3*1*I_KINETIC_F2xz_Pz;
  abcd[110] = 4.0E0*I_KINETIC_K4xy2z_Pz_aa-2.0E0*3*I_KINETIC_H2xy2z_Pz_a;
  abcd[111] = 4.0E0*I_KINETIC_K3x4y_Pz_aa-2.0E0*3*I_KINETIC_H3x2y_Pz_a-2.0E0*2*I_KINETIC_Hx4y_Pz_a+2*3*I_KINETIC_Fx2y_Pz;
  abcd[112] = 4.0E0*I_KINETIC_K3x3yz_Pz_aa-2.0E0*2*I_KINETIC_H3xyz_Pz_a-2.0E0*2*I_KINETIC_Hx3yz_Pz_a+2*2*I_KINETIC_Fxyz_Pz;
  abcd[113] = 4.0E0*I_KINETIC_K3x2y2z_Pz_aa-2.0E0*1*I_KINETIC_H3x2z_Pz_a-2.0E0*2*I_KINETIC_Hx2y2z_Pz_a+2*1*I_KINETIC_Fx2z_Pz;
  abcd[114] = 4.0E0*I_KINETIC_K3xy3z_Pz_aa-2.0E0*2*I_KINETIC_Hxy3z_Pz_a;
  abcd[115] = 4.0E0*I_KINETIC_K2x5y_Pz_aa-2.0E0*4*I_KINETIC_H2x3y_Pz_a-2.0E0*1*I_KINETIC_H5y_Pz_a+4*I_KINETIC_F3y_Pz;
  abcd[116] = 4.0E0*I_KINETIC_K2x4yz_Pz_aa-2.0E0*3*I_KINETIC_H2x2yz_Pz_a-2.0E0*1*I_KINETIC_H4yz_Pz_a+3*I_KINETIC_F2yz_Pz;
  abcd[117] = 4.0E0*I_KINETIC_K2x3y2z_Pz_aa-2.0E0*2*I_KINETIC_H2xy2z_Pz_a-2.0E0*1*I_KINETIC_H3y2z_Pz_a+2*I_KINETIC_Fy2z_Pz;
  abcd[118] = 4.0E0*I_KINETIC_K2x2y3z_Pz_aa-2.0E0*1*I_KINETIC_H2x3z_Pz_a-2.0E0*1*I_KINETIC_H2y3z_Pz_a+1*I_KINETIC_F3z_Pz;
  abcd[119] = 4.0E0*I_KINETIC_K2xy4z_Pz_aa-2.0E0*1*I_KINETIC_Hy4z_Pz_a;
  abcd[120] = 4.0E0*I_KINETIC_Kx6y_Pz_aa-2.0E0*5*I_KINETIC_Hx4y_Pz_a;
  abcd[121] = 4.0E0*I_KINETIC_Kx5yz_Pz_aa-2.0E0*4*I_KINETIC_Hx3yz_Pz_a;
  abcd[122] = 4.0E0*I_KINETIC_Kx4y2z_Pz_aa-2.0E0*3*I_KINETIC_Hx2y2z_Pz_a;
  abcd[123] = 4.0E0*I_KINETIC_Kx3y3z_Pz_aa-2.0E0*2*I_KINETIC_Hxy3z_Pz_a;
  abcd[124] = 4.0E0*I_KINETIC_Kx2y4z_Pz_aa-2.0E0*1*I_KINETIC_Hx4z_Pz_a;
  abcd[125] = 4.0E0*I_KINETIC_Kxy5z_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_P_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_P_aa
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_F_P
   ************************************************************/
  abcd[126] = 4.0E0*I_KINETIC_K6xz_Px_aa-2.0E0*5*I_KINETIC_H4xz_Px_a;
  abcd[127] = 4.0E0*I_KINETIC_K5xyz_Px_aa-2.0E0*4*I_KINETIC_H3xyz_Px_a;
  abcd[128] = 4.0E0*I_KINETIC_K5x2z_Px_aa-2.0E0*1*I_KINETIC_H5x_Px_a-2.0E0*4*I_KINETIC_H3x2z_Px_a+4*1*I_KINETIC_F3x_Px;
  abcd[129] = 4.0E0*I_KINETIC_K4x2yz_Px_aa-2.0E0*3*I_KINETIC_H2x2yz_Px_a;
  abcd[130] = 4.0E0*I_KINETIC_K4xy2z_Px_aa-2.0E0*1*I_KINETIC_H4xy_Px_a-2.0E0*3*I_KINETIC_H2xy2z_Px_a+3*1*I_KINETIC_F2xy_Px;
  abcd[131] = 4.0E0*I_KINETIC_K4x3z_Px_aa-2.0E0*2*I_KINETIC_H4xz_Px_a-2.0E0*3*I_KINETIC_H2x3z_Px_a+3*2*I_KINETIC_F2xz_Px;
  abcd[132] = 4.0E0*I_KINETIC_K3x3yz_Px_aa-2.0E0*2*I_KINETIC_Hx3yz_Px_a;
  abcd[133] = 4.0E0*I_KINETIC_K3x2y2z_Px_aa-2.0E0*1*I_KINETIC_H3x2y_Px_a-2.0E0*2*I_KINETIC_Hx2y2z_Px_a+2*1*I_KINETIC_Fx2y_Px;
  abcd[134] = 4.0E0*I_KINETIC_K3xy3z_Px_aa-2.0E0*2*I_KINETIC_H3xyz_Px_a-2.0E0*2*I_KINETIC_Hxy3z_Px_a+2*2*I_KINETIC_Fxyz_Px;
  abcd[135] = 4.0E0*I_KINETIC_K3x4z_Px_aa-2.0E0*3*I_KINETIC_H3x2z_Px_a-2.0E0*2*I_KINETIC_Hx4z_Px_a+2*3*I_KINETIC_Fx2z_Px;
  abcd[136] = 4.0E0*I_KINETIC_K2x4yz_Px_aa-2.0E0*1*I_KINETIC_H4yz_Px_a;
  abcd[137] = 4.0E0*I_KINETIC_K2x3y2z_Px_aa-2.0E0*1*I_KINETIC_H2x3y_Px_a-2.0E0*1*I_KINETIC_H3y2z_Px_a+1*I_KINETIC_F3y_Px;
  abcd[138] = 4.0E0*I_KINETIC_K2x2y3z_Px_aa-2.0E0*2*I_KINETIC_H2x2yz_Px_a-2.0E0*1*I_KINETIC_H2y3z_Px_a+2*I_KINETIC_F2yz_Px;
  abcd[139] = 4.0E0*I_KINETIC_K2xy4z_Px_aa-2.0E0*3*I_KINETIC_H2xy2z_Px_a-2.0E0*1*I_KINETIC_Hy4z_Px_a+3*I_KINETIC_Fy2z_Px;
  abcd[140] = 4.0E0*I_KINETIC_K2x5z_Px_aa-2.0E0*4*I_KINETIC_H2x3z_Px_a-2.0E0*1*I_KINETIC_H5z_Px_a+4*I_KINETIC_F3z_Px;
  abcd[141] = 4.0E0*I_KINETIC_Kx5yz_Px_aa;
  abcd[142] = 4.0E0*I_KINETIC_Kx4y2z_Px_aa-2.0E0*1*I_KINETIC_Hx4y_Px_a;
  abcd[143] = 4.0E0*I_KINETIC_Kx3y3z_Px_aa-2.0E0*2*I_KINETIC_Hx3yz_Px_a;
  abcd[144] = 4.0E0*I_KINETIC_Kx2y4z_Px_aa-2.0E0*3*I_KINETIC_Hx2y2z_Px_a;
  abcd[145] = 4.0E0*I_KINETIC_Kxy5z_Px_aa-2.0E0*4*I_KINETIC_Hxy3z_Px_a;
  abcd[146] = 4.0E0*I_KINETIC_Kx6z_Px_aa-2.0E0*5*I_KINETIC_Hx4z_Px_a;
  abcd[147] = 4.0E0*I_KINETIC_K6xz_Py_aa-2.0E0*5*I_KINETIC_H4xz_Py_a;
  abcd[148] = 4.0E0*I_KINETIC_K5xyz_Py_aa-2.0E0*4*I_KINETIC_H3xyz_Py_a;
  abcd[149] = 4.0E0*I_KINETIC_K5x2z_Py_aa-2.0E0*1*I_KINETIC_H5x_Py_a-2.0E0*4*I_KINETIC_H3x2z_Py_a+4*1*I_KINETIC_F3x_Py;
  abcd[150] = 4.0E0*I_KINETIC_K4x2yz_Py_aa-2.0E0*3*I_KINETIC_H2x2yz_Py_a;
  abcd[151] = 4.0E0*I_KINETIC_K4xy2z_Py_aa-2.0E0*1*I_KINETIC_H4xy_Py_a-2.0E0*3*I_KINETIC_H2xy2z_Py_a+3*1*I_KINETIC_F2xy_Py;
  abcd[152] = 4.0E0*I_KINETIC_K4x3z_Py_aa-2.0E0*2*I_KINETIC_H4xz_Py_a-2.0E0*3*I_KINETIC_H2x3z_Py_a+3*2*I_KINETIC_F2xz_Py;
  abcd[153] = 4.0E0*I_KINETIC_K3x3yz_Py_aa-2.0E0*2*I_KINETIC_Hx3yz_Py_a;
  abcd[154] = 4.0E0*I_KINETIC_K3x2y2z_Py_aa-2.0E0*1*I_KINETIC_H3x2y_Py_a-2.0E0*2*I_KINETIC_Hx2y2z_Py_a+2*1*I_KINETIC_Fx2y_Py;
  abcd[155] = 4.0E0*I_KINETIC_K3xy3z_Py_aa-2.0E0*2*I_KINETIC_H3xyz_Py_a-2.0E0*2*I_KINETIC_Hxy3z_Py_a+2*2*I_KINETIC_Fxyz_Py;
  abcd[156] = 4.0E0*I_KINETIC_K3x4z_Py_aa-2.0E0*3*I_KINETIC_H3x2z_Py_a-2.0E0*2*I_KINETIC_Hx4z_Py_a+2*3*I_KINETIC_Fx2z_Py;
  abcd[157] = 4.0E0*I_KINETIC_K2x4yz_Py_aa-2.0E0*1*I_KINETIC_H4yz_Py_a;
  abcd[158] = 4.0E0*I_KINETIC_K2x3y2z_Py_aa-2.0E0*1*I_KINETIC_H2x3y_Py_a-2.0E0*1*I_KINETIC_H3y2z_Py_a+1*I_KINETIC_F3y_Py;
  abcd[159] = 4.0E0*I_KINETIC_K2x2y3z_Py_aa-2.0E0*2*I_KINETIC_H2x2yz_Py_a-2.0E0*1*I_KINETIC_H2y3z_Py_a+2*I_KINETIC_F2yz_Py;
  abcd[160] = 4.0E0*I_KINETIC_K2xy4z_Py_aa-2.0E0*3*I_KINETIC_H2xy2z_Py_a-2.0E0*1*I_KINETIC_Hy4z_Py_a+3*I_KINETIC_Fy2z_Py;
  abcd[161] = 4.0E0*I_KINETIC_K2x5z_Py_aa-2.0E0*4*I_KINETIC_H2x3z_Py_a-2.0E0*1*I_KINETIC_H5z_Py_a+4*I_KINETIC_F3z_Py;
  abcd[162] = 4.0E0*I_KINETIC_Kx5yz_Py_aa;
  abcd[163] = 4.0E0*I_KINETIC_Kx4y2z_Py_aa-2.0E0*1*I_KINETIC_Hx4y_Py_a;
  abcd[164] = 4.0E0*I_KINETIC_Kx3y3z_Py_aa-2.0E0*2*I_KINETIC_Hx3yz_Py_a;
  abcd[165] = 4.0E0*I_KINETIC_Kx2y4z_Py_aa-2.0E0*3*I_KINETIC_Hx2y2z_Py_a;
  abcd[166] = 4.0E0*I_KINETIC_Kxy5z_Py_aa-2.0E0*4*I_KINETIC_Hxy3z_Py_a;
  abcd[167] = 4.0E0*I_KINETIC_Kx6z_Py_aa-2.0E0*5*I_KINETIC_Hx4z_Py_a;
  abcd[168] = 4.0E0*I_KINETIC_K6xz_Pz_aa-2.0E0*5*I_KINETIC_H4xz_Pz_a;
  abcd[169] = 4.0E0*I_KINETIC_K5xyz_Pz_aa-2.0E0*4*I_KINETIC_H3xyz_Pz_a;
  abcd[170] = 4.0E0*I_KINETIC_K5x2z_Pz_aa-2.0E0*1*I_KINETIC_H5x_Pz_a-2.0E0*4*I_KINETIC_H3x2z_Pz_a+4*1*I_KINETIC_F3x_Pz;
  abcd[171] = 4.0E0*I_KINETIC_K4x2yz_Pz_aa-2.0E0*3*I_KINETIC_H2x2yz_Pz_a;
  abcd[172] = 4.0E0*I_KINETIC_K4xy2z_Pz_aa-2.0E0*1*I_KINETIC_H4xy_Pz_a-2.0E0*3*I_KINETIC_H2xy2z_Pz_a+3*1*I_KINETIC_F2xy_Pz;
  abcd[173] = 4.0E0*I_KINETIC_K4x3z_Pz_aa-2.0E0*2*I_KINETIC_H4xz_Pz_a-2.0E0*3*I_KINETIC_H2x3z_Pz_a+3*2*I_KINETIC_F2xz_Pz;
  abcd[174] = 4.0E0*I_KINETIC_K3x3yz_Pz_aa-2.0E0*2*I_KINETIC_Hx3yz_Pz_a;
  abcd[175] = 4.0E0*I_KINETIC_K3x2y2z_Pz_aa-2.0E0*1*I_KINETIC_H3x2y_Pz_a-2.0E0*2*I_KINETIC_Hx2y2z_Pz_a+2*1*I_KINETIC_Fx2y_Pz;
  abcd[176] = 4.0E0*I_KINETIC_K3xy3z_Pz_aa-2.0E0*2*I_KINETIC_H3xyz_Pz_a-2.0E0*2*I_KINETIC_Hxy3z_Pz_a+2*2*I_KINETIC_Fxyz_Pz;
  abcd[177] = 4.0E0*I_KINETIC_K3x4z_Pz_aa-2.0E0*3*I_KINETIC_H3x2z_Pz_a-2.0E0*2*I_KINETIC_Hx4z_Pz_a+2*3*I_KINETIC_Fx2z_Pz;
  abcd[178] = 4.0E0*I_KINETIC_K2x4yz_Pz_aa-2.0E0*1*I_KINETIC_H4yz_Pz_a;
  abcd[179] = 4.0E0*I_KINETIC_K2x3y2z_Pz_aa-2.0E0*1*I_KINETIC_H2x3y_Pz_a-2.0E0*1*I_KINETIC_H3y2z_Pz_a+1*I_KINETIC_F3y_Pz;
  abcd[180] = 4.0E0*I_KINETIC_K2x2y3z_Pz_aa-2.0E0*2*I_KINETIC_H2x2yz_Pz_a-2.0E0*1*I_KINETIC_H2y3z_Pz_a+2*I_KINETIC_F2yz_Pz;
  abcd[181] = 4.0E0*I_KINETIC_K2xy4z_Pz_aa-2.0E0*3*I_KINETIC_H2xy2z_Pz_a-2.0E0*1*I_KINETIC_Hy4z_Pz_a+3*I_KINETIC_Fy2z_Pz;
  abcd[182] = 4.0E0*I_KINETIC_K2x5z_Pz_aa-2.0E0*4*I_KINETIC_H2x3z_Pz_a-2.0E0*1*I_KINETIC_H5z_Pz_a+4*I_KINETIC_F3z_Pz;
  abcd[183] = 4.0E0*I_KINETIC_Kx5yz_Pz_aa;
  abcd[184] = 4.0E0*I_KINETIC_Kx4y2z_Pz_aa-2.0E0*1*I_KINETIC_Hx4y_Pz_a;
  abcd[185] = 4.0E0*I_KINETIC_Kx3y3z_Pz_aa-2.0E0*2*I_KINETIC_Hx3yz_Pz_a;
  abcd[186] = 4.0E0*I_KINETIC_Kx2y4z_Pz_aa-2.0E0*3*I_KINETIC_Hx2y2z_Pz_a;
  abcd[187] = 4.0E0*I_KINETIC_Kxy5z_Pz_aa-2.0E0*4*I_KINETIC_Hxy3z_Pz_a;
  abcd[188] = 4.0E0*I_KINETIC_Kx6z_Pz_aa-2.0E0*5*I_KINETIC_Hx4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_P_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_P_aa
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_F_P
   ************************************************************/
  abcd[189] = 4.0E0*I_KINETIC_K5x2y_Px_aa-2.0E0*1*I_KINETIC_H5x_Px_a;
  abcd[190] = 4.0E0*I_KINETIC_K4x3y_Px_aa-2.0E0*1*I_KINETIC_H4xy_Px_a-2.0E0*2*I_KINETIC_H4xy_Px_a;
  abcd[191] = 4.0E0*I_KINETIC_K4x2yz_Px_aa-2.0E0*1*I_KINETIC_H4xz_Px_a;
  abcd[192] = 4.0E0*I_KINETIC_K3x4y_Px_aa-2.0E0*2*I_KINETIC_H3x2y_Px_a-2.0E0*3*I_KINETIC_H3x2y_Px_a+2*1*I_KINETIC_F3x_Px;
  abcd[193] = 4.0E0*I_KINETIC_K3x3yz_Px_aa-2.0E0*1*I_KINETIC_H3xyz_Px_a-2.0E0*2*I_KINETIC_H3xyz_Px_a;
  abcd[194] = 4.0E0*I_KINETIC_K3x2y2z_Px_aa-2.0E0*1*I_KINETIC_H3x2z_Px_a;
  abcd[195] = 4.0E0*I_KINETIC_K2x5y_Px_aa-2.0E0*3*I_KINETIC_H2x3y_Px_a-2.0E0*4*I_KINETIC_H2x3y_Px_a+3*2*I_KINETIC_F2xy_Px;
  abcd[196] = 4.0E0*I_KINETIC_K2x4yz_Px_aa-2.0E0*2*I_KINETIC_H2x2yz_Px_a-2.0E0*3*I_KINETIC_H2x2yz_Px_a+2*1*I_KINETIC_F2xz_Px;
  abcd[197] = 4.0E0*I_KINETIC_K2x3y2z_Px_aa-2.0E0*1*I_KINETIC_H2xy2z_Px_a-2.0E0*2*I_KINETIC_H2xy2z_Px_a;
  abcd[198] = 4.0E0*I_KINETIC_K2x2y3z_Px_aa-2.0E0*1*I_KINETIC_H2x3z_Px_a;
  abcd[199] = 4.0E0*I_KINETIC_Kx6y_Px_aa-2.0E0*4*I_KINETIC_Hx4y_Px_a-2.0E0*5*I_KINETIC_Hx4y_Px_a+4*3*I_KINETIC_Fx2y_Px;
  abcd[200] = 4.0E0*I_KINETIC_Kx5yz_Px_aa-2.0E0*3*I_KINETIC_Hx3yz_Px_a-2.0E0*4*I_KINETIC_Hx3yz_Px_a+3*2*I_KINETIC_Fxyz_Px;
  abcd[201] = 4.0E0*I_KINETIC_Kx4y2z_Px_aa-2.0E0*2*I_KINETIC_Hx2y2z_Px_a-2.0E0*3*I_KINETIC_Hx2y2z_Px_a+2*1*I_KINETIC_Fx2z_Px;
  abcd[202] = 4.0E0*I_KINETIC_Kx3y3z_Px_aa-2.0E0*1*I_KINETIC_Hxy3z_Px_a-2.0E0*2*I_KINETIC_Hxy3z_Px_a;
  abcd[203] = 4.0E0*I_KINETIC_Kx2y4z_Px_aa-2.0E0*1*I_KINETIC_Hx4z_Px_a;
  abcd[204] = 4.0E0*I_KINETIC_K7y_Px_aa-2.0E0*5*I_KINETIC_H5y_Px_a-2.0E0*6*I_KINETIC_H5y_Px_a+5*4*I_KINETIC_F3y_Px;
  abcd[205] = 4.0E0*I_KINETIC_K6yz_Px_aa-2.0E0*4*I_KINETIC_H4yz_Px_a-2.0E0*5*I_KINETIC_H4yz_Px_a+4*3*I_KINETIC_F2yz_Px;
  abcd[206] = 4.0E0*I_KINETIC_K5y2z_Px_aa-2.0E0*3*I_KINETIC_H3y2z_Px_a-2.0E0*4*I_KINETIC_H3y2z_Px_a+3*2*I_KINETIC_Fy2z_Px;
  abcd[207] = 4.0E0*I_KINETIC_K4y3z_Px_aa-2.0E0*2*I_KINETIC_H2y3z_Px_a-2.0E0*3*I_KINETIC_H2y3z_Px_a+2*1*I_KINETIC_F3z_Px;
  abcd[208] = 4.0E0*I_KINETIC_K3y4z_Px_aa-2.0E0*1*I_KINETIC_Hy4z_Px_a-2.0E0*2*I_KINETIC_Hy4z_Px_a;
  abcd[209] = 4.0E0*I_KINETIC_K2y5z_Px_aa-2.0E0*1*I_KINETIC_H5z_Px_a;
  abcd[210] = 4.0E0*I_KINETIC_K5x2y_Py_aa-2.0E0*1*I_KINETIC_H5x_Py_a;
  abcd[211] = 4.0E0*I_KINETIC_K4x3y_Py_aa-2.0E0*1*I_KINETIC_H4xy_Py_a-2.0E0*2*I_KINETIC_H4xy_Py_a;
  abcd[212] = 4.0E0*I_KINETIC_K4x2yz_Py_aa-2.0E0*1*I_KINETIC_H4xz_Py_a;
  abcd[213] = 4.0E0*I_KINETIC_K3x4y_Py_aa-2.0E0*2*I_KINETIC_H3x2y_Py_a-2.0E0*3*I_KINETIC_H3x2y_Py_a+2*1*I_KINETIC_F3x_Py;
  abcd[214] = 4.0E0*I_KINETIC_K3x3yz_Py_aa-2.0E0*1*I_KINETIC_H3xyz_Py_a-2.0E0*2*I_KINETIC_H3xyz_Py_a;
  abcd[215] = 4.0E0*I_KINETIC_K3x2y2z_Py_aa-2.0E0*1*I_KINETIC_H3x2z_Py_a;
  abcd[216] = 4.0E0*I_KINETIC_K2x5y_Py_aa-2.0E0*3*I_KINETIC_H2x3y_Py_a-2.0E0*4*I_KINETIC_H2x3y_Py_a+3*2*I_KINETIC_F2xy_Py;
  abcd[217] = 4.0E0*I_KINETIC_K2x4yz_Py_aa-2.0E0*2*I_KINETIC_H2x2yz_Py_a-2.0E0*3*I_KINETIC_H2x2yz_Py_a+2*1*I_KINETIC_F2xz_Py;
  abcd[218] = 4.0E0*I_KINETIC_K2x3y2z_Py_aa-2.0E0*1*I_KINETIC_H2xy2z_Py_a-2.0E0*2*I_KINETIC_H2xy2z_Py_a;
  abcd[219] = 4.0E0*I_KINETIC_K2x2y3z_Py_aa-2.0E0*1*I_KINETIC_H2x3z_Py_a;
  abcd[220] = 4.0E0*I_KINETIC_Kx6y_Py_aa-2.0E0*4*I_KINETIC_Hx4y_Py_a-2.0E0*5*I_KINETIC_Hx4y_Py_a+4*3*I_KINETIC_Fx2y_Py;
  abcd[221] = 4.0E0*I_KINETIC_Kx5yz_Py_aa-2.0E0*3*I_KINETIC_Hx3yz_Py_a-2.0E0*4*I_KINETIC_Hx3yz_Py_a+3*2*I_KINETIC_Fxyz_Py;
  abcd[222] = 4.0E0*I_KINETIC_Kx4y2z_Py_aa-2.0E0*2*I_KINETIC_Hx2y2z_Py_a-2.0E0*3*I_KINETIC_Hx2y2z_Py_a+2*1*I_KINETIC_Fx2z_Py;
  abcd[223] = 4.0E0*I_KINETIC_Kx3y3z_Py_aa-2.0E0*1*I_KINETIC_Hxy3z_Py_a-2.0E0*2*I_KINETIC_Hxy3z_Py_a;
  abcd[224] = 4.0E0*I_KINETIC_Kx2y4z_Py_aa-2.0E0*1*I_KINETIC_Hx4z_Py_a;
  abcd[225] = 4.0E0*I_KINETIC_K7y_Py_aa-2.0E0*5*I_KINETIC_H5y_Py_a-2.0E0*6*I_KINETIC_H5y_Py_a+5*4*I_KINETIC_F3y_Py;
  abcd[226] = 4.0E0*I_KINETIC_K6yz_Py_aa-2.0E0*4*I_KINETIC_H4yz_Py_a-2.0E0*5*I_KINETIC_H4yz_Py_a+4*3*I_KINETIC_F2yz_Py;
  abcd[227] = 4.0E0*I_KINETIC_K5y2z_Py_aa-2.0E0*3*I_KINETIC_H3y2z_Py_a-2.0E0*4*I_KINETIC_H3y2z_Py_a+3*2*I_KINETIC_Fy2z_Py;
  abcd[228] = 4.0E0*I_KINETIC_K4y3z_Py_aa-2.0E0*2*I_KINETIC_H2y3z_Py_a-2.0E0*3*I_KINETIC_H2y3z_Py_a+2*1*I_KINETIC_F3z_Py;
  abcd[229] = 4.0E0*I_KINETIC_K3y4z_Py_aa-2.0E0*1*I_KINETIC_Hy4z_Py_a-2.0E0*2*I_KINETIC_Hy4z_Py_a;
  abcd[230] = 4.0E0*I_KINETIC_K2y5z_Py_aa-2.0E0*1*I_KINETIC_H5z_Py_a;
  abcd[231] = 4.0E0*I_KINETIC_K5x2y_Pz_aa-2.0E0*1*I_KINETIC_H5x_Pz_a;
  abcd[232] = 4.0E0*I_KINETIC_K4x3y_Pz_aa-2.0E0*1*I_KINETIC_H4xy_Pz_a-2.0E0*2*I_KINETIC_H4xy_Pz_a;
  abcd[233] = 4.0E0*I_KINETIC_K4x2yz_Pz_aa-2.0E0*1*I_KINETIC_H4xz_Pz_a;
  abcd[234] = 4.0E0*I_KINETIC_K3x4y_Pz_aa-2.0E0*2*I_KINETIC_H3x2y_Pz_a-2.0E0*3*I_KINETIC_H3x2y_Pz_a+2*1*I_KINETIC_F3x_Pz;
  abcd[235] = 4.0E0*I_KINETIC_K3x3yz_Pz_aa-2.0E0*1*I_KINETIC_H3xyz_Pz_a-2.0E0*2*I_KINETIC_H3xyz_Pz_a;
  abcd[236] = 4.0E0*I_KINETIC_K3x2y2z_Pz_aa-2.0E0*1*I_KINETIC_H3x2z_Pz_a;
  abcd[237] = 4.0E0*I_KINETIC_K2x5y_Pz_aa-2.0E0*3*I_KINETIC_H2x3y_Pz_a-2.0E0*4*I_KINETIC_H2x3y_Pz_a+3*2*I_KINETIC_F2xy_Pz;
  abcd[238] = 4.0E0*I_KINETIC_K2x4yz_Pz_aa-2.0E0*2*I_KINETIC_H2x2yz_Pz_a-2.0E0*3*I_KINETIC_H2x2yz_Pz_a+2*1*I_KINETIC_F2xz_Pz;
  abcd[239] = 4.0E0*I_KINETIC_K2x3y2z_Pz_aa-2.0E0*1*I_KINETIC_H2xy2z_Pz_a-2.0E0*2*I_KINETIC_H2xy2z_Pz_a;
  abcd[240] = 4.0E0*I_KINETIC_K2x2y3z_Pz_aa-2.0E0*1*I_KINETIC_H2x3z_Pz_a;
  abcd[241] = 4.0E0*I_KINETIC_Kx6y_Pz_aa-2.0E0*4*I_KINETIC_Hx4y_Pz_a-2.0E0*5*I_KINETIC_Hx4y_Pz_a+4*3*I_KINETIC_Fx2y_Pz;
  abcd[242] = 4.0E0*I_KINETIC_Kx5yz_Pz_aa-2.0E0*3*I_KINETIC_Hx3yz_Pz_a-2.0E0*4*I_KINETIC_Hx3yz_Pz_a+3*2*I_KINETIC_Fxyz_Pz;
  abcd[243] = 4.0E0*I_KINETIC_Kx4y2z_Pz_aa-2.0E0*2*I_KINETIC_Hx2y2z_Pz_a-2.0E0*3*I_KINETIC_Hx2y2z_Pz_a+2*1*I_KINETIC_Fx2z_Pz;
  abcd[244] = 4.0E0*I_KINETIC_Kx3y3z_Pz_aa-2.0E0*1*I_KINETIC_Hxy3z_Pz_a-2.0E0*2*I_KINETIC_Hxy3z_Pz_a;
  abcd[245] = 4.0E0*I_KINETIC_Kx2y4z_Pz_aa-2.0E0*1*I_KINETIC_Hx4z_Pz_a;
  abcd[246] = 4.0E0*I_KINETIC_K7y_Pz_aa-2.0E0*5*I_KINETIC_H5y_Pz_a-2.0E0*6*I_KINETIC_H5y_Pz_a+5*4*I_KINETIC_F3y_Pz;
  abcd[247] = 4.0E0*I_KINETIC_K6yz_Pz_aa-2.0E0*4*I_KINETIC_H4yz_Pz_a-2.0E0*5*I_KINETIC_H4yz_Pz_a+4*3*I_KINETIC_F2yz_Pz;
  abcd[248] = 4.0E0*I_KINETIC_K5y2z_Pz_aa-2.0E0*3*I_KINETIC_H3y2z_Pz_a-2.0E0*4*I_KINETIC_H3y2z_Pz_a+3*2*I_KINETIC_Fy2z_Pz;
  abcd[249] = 4.0E0*I_KINETIC_K4y3z_Pz_aa-2.0E0*2*I_KINETIC_H2y3z_Pz_a-2.0E0*3*I_KINETIC_H2y3z_Pz_a+2*1*I_KINETIC_F3z_Pz;
  abcd[250] = 4.0E0*I_KINETIC_K3y4z_Pz_aa-2.0E0*1*I_KINETIC_Hy4z_Pz_a-2.0E0*2*I_KINETIC_Hy4z_Pz_a;
  abcd[251] = 4.0E0*I_KINETIC_K2y5z_Pz_aa-2.0E0*1*I_KINETIC_H5z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_P_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_P_aa
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_F_P
   ************************************************************/
  abcd[252] = 4.0E0*I_KINETIC_K5xyz_Px_aa;
  abcd[253] = 4.0E0*I_KINETIC_K4x2yz_Px_aa-2.0E0*1*I_KINETIC_H4xz_Px_a;
  abcd[254] = 4.0E0*I_KINETIC_K4xy2z_Px_aa-2.0E0*1*I_KINETIC_H4xy_Px_a;
  abcd[255] = 4.0E0*I_KINETIC_K3x3yz_Px_aa-2.0E0*2*I_KINETIC_H3xyz_Px_a;
  abcd[256] = 4.0E0*I_KINETIC_K3x2y2z_Px_aa-2.0E0*1*I_KINETIC_H3x2y_Px_a-2.0E0*1*I_KINETIC_H3x2z_Px_a+1*I_KINETIC_F3x_Px;
  abcd[257] = 4.0E0*I_KINETIC_K3xy3z_Px_aa-2.0E0*2*I_KINETIC_H3xyz_Px_a;
  abcd[258] = 4.0E0*I_KINETIC_K2x4yz_Px_aa-2.0E0*3*I_KINETIC_H2x2yz_Px_a;
  abcd[259] = 4.0E0*I_KINETIC_K2x3y2z_Px_aa-2.0E0*1*I_KINETIC_H2x3y_Px_a-2.0E0*2*I_KINETIC_H2xy2z_Px_a+2*1*I_KINETIC_F2xy_Px;
  abcd[260] = 4.0E0*I_KINETIC_K2x2y3z_Px_aa-2.0E0*2*I_KINETIC_H2x2yz_Px_a-2.0E0*1*I_KINETIC_H2x3z_Px_a+2*I_KINETIC_F2xz_Px;
  abcd[261] = 4.0E0*I_KINETIC_K2xy4z_Px_aa-2.0E0*3*I_KINETIC_H2xy2z_Px_a;
  abcd[262] = 4.0E0*I_KINETIC_Kx5yz_Px_aa-2.0E0*4*I_KINETIC_Hx3yz_Px_a;
  abcd[263] = 4.0E0*I_KINETIC_Kx4y2z_Px_aa-2.0E0*1*I_KINETIC_Hx4y_Px_a-2.0E0*3*I_KINETIC_Hx2y2z_Px_a+3*1*I_KINETIC_Fx2y_Px;
  abcd[264] = 4.0E0*I_KINETIC_Kx3y3z_Px_aa-2.0E0*2*I_KINETIC_Hx3yz_Px_a-2.0E0*2*I_KINETIC_Hxy3z_Px_a+2*2*I_KINETIC_Fxyz_Px;
  abcd[265] = 4.0E0*I_KINETIC_Kx2y4z_Px_aa-2.0E0*3*I_KINETIC_Hx2y2z_Px_a-2.0E0*1*I_KINETIC_Hx4z_Px_a+3*I_KINETIC_Fx2z_Px;
  abcd[266] = 4.0E0*I_KINETIC_Kxy5z_Px_aa-2.0E0*4*I_KINETIC_Hxy3z_Px_a;
  abcd[267] = 4.0E0*I_KINETIC_K6yz_Px_aa-2.0E0*5*I_KINETIC_H4yz_Px_a;
  abcd[268] = 4.0E0*I_KINETIC_K5y2z_Px_aa-2.0E0*1*I_KINETIC_H5y_Px_a-2.0E0*4*I_KINETIC_H3y2z_Px_a+4*1*I_KINETIC_F3y_Px;
  abcd[269] = 4.0E0*I_KINETIC_K4y3z_Px_aa-2.0E0*2*I_KINETIC_H4yz_Px_a-2.0E0*3*I_KINETIC_H2y3z_Px_a+3*2*I_KINETIC_F2yz_Px;
  abcd[270] = 4.0E0*I_KINETIC_K3y4z_Px_aa-2.0E0*3*I_KINETIC_H3y2z_Px_a-2.0E0*2*I_KINETIC_Hy4z_Px_a+2*3*I_KINETIC_Fy2z_Px;
  abcd[271] = 4.0E0*I_KINETIC_K2y5z_Px_aa-2.0E0*4*I_KINETIC_H2y3z_Px_a-2.0E0*1*I_KINETIC_H5z_Px_a+4*I_KINETIC_F3z_Px;
  abcd[272] = 4.0E0*I_KINETIC_Ky6z_Px_aa-2.0E0*5*I_KINETIC_Hy4z_Px_a;
  abcd[273] = 4.0E0*I_KINETIC_K5xyz_Py_aa;
  abcd[274] = 4.0E0*I_KINETIC_K4x2yz_Py_aa-2.0E0*1*I_KINETIC_H4xz_Py_a;
  abcd[275] = 4.0E0*I_KINETIC_K4xy2z_Py_aa-2.0E0*1*I_KINETIC_H4xy_Py_a;
  abcd[276] = 4.0E0*I_KINETIC_K3x3yz_Py_aa-2.0E0*2*I_KINETIC_H3xyz_Py_a;
  abcd[277] = 4.0E0*I_KINETIC_K3x2y2z_Py_aa-2.0E0*1*I_KINETIC_H3x2y_Py_a-2.0E0*1*I_KINETIC_H3x2z_Py_a+1*I_KINETIC_F3x_Py;
  abcd[278] = 4.0E0*I_KINETIC_K3xy3z_Py_aa-2.0E0*2*I_KINETIC_H3xyz_Py_a;
  abcd[279] = 4.0E0*I_KINETIC_K2x4yz_Py_aa-2.0E0*3*I_KINETIC_H2x2yz_Py_a;
  abcd[280] = 4.0E0*I_KINETIC_K2x3y2z_Py_aa-2.0E0*1*I_KINETIC_H2x3y_Py_a-2.0E0*2*I_KINETIC_H2xy2z_Py_a+2*1*I_KINETIC_F2xy_Py;
  abcd[281] = 4.0E0*I_KINETIC_K2x2y3z_Py_aa-2.0E0*2*I_KINETIC_H2x2yz_Py_a-2.0E0*1*I_KINETIC_H2x3z_Py_a+2*I_KINETIC_F2xz_Py;
  abcd[282] = 4.0E0*I_KINETIC_K2xy4z_Py_aa-2.0E0*3*I_KINETIC_H2xy2z_Py_a;
  abcd[283] = 4.0E0*I_KINETIC_Kx5yz_Py_aa-2.0E0*4*I_KINETIC_Hx3yz_Py_a;
  abcd[284] = 4.0E0*I_KINETIC_Kx4y2z_Py_aa-2.0E0*1*I_KINETIC_Hx4y_Py_a-2.0E0*3*I_KINETIC_Hx2y2z_Py_a+3*1*I_KINETIC_Fx2y_Py;
  abcd[285] = 4.0E0*I_KINETIC_Kx3y3z_Py_aa-2.0E0*2*I_KINETIC_Hx3yz_Py_a-2.0E0*2*I_KINETIC_Hxy3z_Py_a+2*2*I_KINETIC_Fxyz_Py;
  abcd[286] = 4.0E0*I_KINETIC_Kx2y4z_Py_aa-2.0E0*3*I_KINETIC_Hx2y2z_Py_a-2.0E0*1*I_KINETIC_Hx4z_Py_a+3*I_KINETIC_Fx2z_Py;
  abcd[287] = 4.0E0*I_KINETIC_Kxy5z_Py_aa-2.0E0*4*I_KINETIC_Hxy3z_Py_a;
  abcd[288] = 4.0E0*I_KINETIC_K6yz_Py_aa-2.0E0*5*I_KINETIC_H4yz_Py_a;
  abcd[289] = 4.0E0*I_KINETIC_K5y2z_Py_aa-2.0E0*1*I_KINETIC_H5y_Py_a-2.0E0*4*I_KINETIC_H3y2z_Py_a+4*1*I_KINETIC_F3y_Py;
  abcd[290] = 4.0E0*I_KINETIC_K4y3z_Py_aa-2.0E0*2*I_KINETIC_H4yz_Py_a-2.0E0*3*I_KINETIC_H2y3z_Py_a+3*2*I_KINETIC_F2yz_Py;
  abcd[291] = 4.0E0*I_KINETIC_K3y4z_Py_aa-2.0E0*3*I_KINETIC_H3y2z_Py_a-2.0E0*2*I_KINETIC_Hy4z_Py_a+2*3*I_KINETIC_Fy2z_Py;
  abcd[292] = 4.0E0*I_KINETIC_K2y5z_Py_aa-2.0E0*4*I_KINETIC_H2y3z_Py_a-2.0E0*1*I_KINETIC_H5z_Py_a+4*I_KINETIC_F3z_Py;
  abcd[293] = 4.0E0*I_KINETIC_Ky6z_Py_aa-2.0E0*5*I_KINETIC_Hy4z_Py_a;
  abcd[294] = 4.0E0*I_KINETIC_K5xyz_Pz_aa;
  abcd[295] = 4.0E0*I_KINETIC_K4x2yz_Pz_aa-2.0E0*1*I_KINETIC_H4xz_Pz_a;
  abcd[296] = 4.0E0*I_KINETIC_K4xy2z_Pz_aa-2.0E0*1*I_KINETIC_H4xy_Pz_a;
  abcd[297] = 4.0E0*I_KINETIC_K3x3yz_Pz_aa-2.0E0*2*I_KINETIC_H3xyz_Pz_a;
  abcd[298] = 4.0E0*I_KINETIC_K3x2y2z_Pz_aa-2.0E0*1*I_KINETIC_H3x2y_Pz_a-2.0E0*1*I_KINETIC_H3x2z_Pz_a+1*I_KINETIC_F3x_Pz;
  abcd[299] = 4.0E0*I_KINETIC_K3xy3z_Pz_aa-2.0E0*2*I_KINETIC_H3xyz_Pz_a;
  abcd[300] = 4.0E0*I_KINETIC_K2x4yz_Pz_aa-2.0E0*3*I_KINETIC_H2x2yz_Pz_a;
  abcd[301] = 4.0E0*I_KINETIC_K2x3y2z_Pz_aa-2.0E0*1*I_KINETIC_H2x3y_Pz_a-2.0E0*2*I_KINETIC_H2xy2z_Pz_a+2*1*I_KINETIC_F2xy_Pz;
  abcd[302] = 4.0E0*I_KINETIC_K2x2y3z_Pz_aa-2.0E0*2*I_KINETIC_H2x2yz_Pz_a-2.0E0*1*I_KINETIC_H2x3z_Pz_a+2*I_KINETIC_F2xz_Pz;
  abcd[303] = 4.0E0*I_KINETIC_K2xy4z_Pz_aa-2.0E0*3*I_KINETIC_H2xy2z_Pz_a;
  abcd[304] = 4.0E0*I_KINETIC_Kx5yz_Pz_aa-2.0E0*4*I_KINETIC_Hx3yz_Pz_a;
  abcd[305] = 4.0E0*I_KINETIC_Kx4y2z_Pz_aa-2.0E0*1*I_KINETIC_Hx4y_Pz_a-2.0E0*3*I_KINETIC_Hx2y2z_Pz_a+3*1*I_KINETIC_Fx2y_Pz;
  abcd[306] = 4.0E0*I_KINETIC_Kx3y3z_Pz_aa-2.0E0*2*I_KINETIC_Hx3yz_Pz_a-2.0E0*2*I_KINETIC_Hxy3z_Pz_a+2*2*I_KINETIC_Fxyz_Pz;
  abcd[307] = 4.0E0*I_KINETIC_Kx2y4z_Pz_aa-2.0E0*3*I_KINETIC_Hx2y2z_Pz_a-2.0E0*1*I_KINETIC_Hx4z_Pz_a+3*I_KINETIC_Fx2z_Pz;
  abcd[308] = 4.0E0*I_KINETIC_Kxy5z_Pz_aa-2.0E0*4*I_KINETIC_Hxy3z_Pz_a;
  abcd[309] = 4.0E0*I_KINETIC_K6yz_Pz_aa-2.0E0*5*I_KINETIC_H4yz_Pz_a;
  abcd[310] = 4.0E0*I_KINETIC_K5y2z_Pz_aa-2.0E0*1*I_KINETIC_H5y_Pz_a-2.0E0*4*I_KINETIC_H3y2z_Pz_a+4*1*I_KINETIC_F3y_Pz;
  abcd[311] = 4.0E0*I_KINETIC_K4y3z_Pz_aa-2.0E0*2*I_KINETIC_H4yz_Pz_a-2.0E0*3*I_KINETIC_H2y3z_Pz_a+3*2*I_KINETIC_F2yz_Pz;
  abcd[312] = 4.0E0*I_KINETIC_K3y4z_Pz_aa-2.0E0*3*I_KINETIC_H3y2z_Pz_a-2.0E0*2*I_KINETIC_Hy4z_Pz_a+2*3*I_KINETIC_Fy2z_Pz;
  abcd[313] = 4.0E0*I_KINETIC_K2y5z_Pz_aa-2.0E0*4*I_KINETIC_H2y3z_Pz_a-2.0E0*1*I_KINETIC_H5z_Pz_a+4*I_KINETIC_F3z_Pz;
  abcd[314] = 4.0E0*I_KINETIC_Ky6z_Pz_aa-2.0E0*5*I_KINETIC_Hy4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_P_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_P_aa
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_H_P_a
   * RHS shell quartet name: SQ_KINETIC_F_P
   ************************************************************/
  abcd[315] = 4.0E0*I_KINETIC_K5x2z_Px_aa-2.0E0*1*I_KINETIC_H5x_Px_a;
  abcd[316] = 4.0E0*I_KINETIC_K4xy2z_Px_aa-2.0E0*1*I_KINETIC_H4xy_Px_a;
  abcd[317] = 4.0E0*I_KINETIC_K4x3z_Px_aa-2.0E0*1*I_KINETIC_H4xz_Px_a-2.0E0*2*I_KINETIC_H4xz_Px_a;
  abcd[318] = 4.0E0*I_KINETIC_K3x2y2z_Px_aa-2.0E0*1*I_KINETIC_H3x2y_Px_a;
  abcd[319] = 4.0E0*I_KINETIC_K3xy3z_Px_aa-2.0E0*1*I_KINETIC_H3xyz_Px_a-2.0E0*2*I_KINETIC_H3xyz_Px_a;
  abcd[320] = 4.0E0*I_KINETIC_K3x4z_Px_aa-2.0E0*2*I_KINETIC_H3x2z_Px_a-2.0E0*3*I_KINETIC_H3x2z_Px_a+2*1*I_KINETIC_F3x_Px;
  abcd[321] = 4.0E0*I_KINETIC_K2x3y2z_Px_aa-2.0E0*1*I_KINETIC_H2x3y_Px_a;
  abcd[322] = 4.0E0*I_KINETIC_K2x2y3z_Px_aa-2.0E0*1*I_KINETIC_H2x2yz_Px_a-2.0E0*2*I_KINETIC_H2x2yz_Px_a;
  abcd[323] = 4.0E0*I_KINETIC_K2xy4z_Px_aa-2.0E0*2*I_KINETIC_H2xy2z_Px_a-2.0E0*3*I_KINETIC_H2xy2z_Px_a+2*1*I_KINETIC_F2xy_Px;
  abcd[324] = 4.0E0*I_KINETIC_K2x5z_Px_aa-2.0E0*3*I_KINETIC_H2x3z_Px_a-2.0E0*4*I_KINETIC_H2x3z_Px_a+3*2*I_KINETIC_F2xz_Px;
  abcd[325] = 4.0E0*I_KINETIC_Kx4y2z_Px_aa-2.0E0*1*I_KINETIC_Hx4y_Px_a;
  abcd[326] = 4.0E0*I_KINETIC_Kx3y3z_Px_aa-2.0E0*1*I_KINETIC_Hx3yz_Px_a-2.0E0*2*I_KINETIC_Hx3yz_Px_a;
  abcd[327] = 4.0E0*I_KINETIC_Kx2y4z_Px_aa-2.0E0*2*I_KINETIC_Hx2y2z_Px_a-2.0E0*3*I_KINETIC_Hx2y2z_Px_a+2*1*I_KINETIC_Fx2y_Px;
  abcd[328] = 4.0E0*I_KINETIC_Kxy5z_Px_aa-2.0E0*3*I_KINETIC_Hxy3z_Px_a-2.0E0*4*I_KINETIC_Hxy3z_Px_a+3*2*I_KINETIC_Fxyz_Px;
  abcd[329] = 4.0E0*I_KINETIC_Kx6z_Px_aa-2.0E0*4*I_KINETIC_Hx4z_Px_a-2.0E0*5*I_KINETIC_Hx4z_Px_a+4*3*I_KINETIC_Fx2z_Px;
  abcd[330] = 4.0E0*I_KINETIC_K5y2z_Px_aa-2.0E0*1*I_KINETIC_H5y_Px_a;
  abcd[331] = 4.0E0*I_KINETIC_K4y3z_Px_aa-2.0E0*1*I_KINETIC_H4yz_Px_a-2.0E0*2*I_KINETIC_H4yz_Px_a;
  abcd[332] = 4.0E0*I_KINETIC_K3y4z_Px_aa-2.0E0*2*I_KINETIC_H3y2z_Px_a-2.0E0*3*I_KINETIC_H3y2z_Px_a+2*1*I_KINETIC_F3y_Px;
  abcd[333] = 4.0E0*I_KINETIC_K2y5z_Px_aa-2.0E0*3*I_KINETIC_H2y3z_Px_a-2.0E0*4*I_KINETIC_H2y3z_Px_a+3*2*I_KINETIC_F2yz_Px;
  abcd[334] = 4.0E0*I_KINETIC_Ky6z_Px_aa-2.0E0*4*I_KINETIC_Hy4z_Px_a-2.0E0*5*I_KINETIC_Hy4z_Px_a+4*3*I_KINETIC_Fy2z_Px;
  abcd[335] = 4.0E0*I_KINETIC_K7z_Px_aa-2.0E0*5*I_KINETIC_H5z_Px_a-2.0E0*6*I_KINETIC_H5z_Px_a+5*4*I_KINETIC_F3z_Px;
  abcd[336] = 4.0E0*I_KINETIC_K5x2z_Py_aa-2.0E0*1*I_KINETIC_H5x_Py_a;
  abcd[337] = 4.0E0*I_KINETIC_K4xy2z_Py_aa-2.0E0*1*I_KINETIC_H4xy_Py_a;
  abcd[338] = 4.0E0*I_KINETIC_K4x3z_Py_aa-2.0E0*1*I_KINETIC_H4xz_Py_a-2.0E0*2*I_KINETIC_H4xz_Py_a;
  abcd[339] = 4.0E0*I_KINETIC_K3x2y2z_Py_aa-2.0E0*1*I_KINETIC_H3x2y_Py_a;
  abcd[340] = 4.0E0*I_KINETIC_K3xy3z_Py_aa-2.0E0*1*I_KINETIC_H3xyz_Py_a-2.0E0*2*I_KINETIC_H3xyz_Py_a;
  abcd[341] = 4.0E0*I_KINETIC_K3x4z_Py_aa-2.0E0*2*I_KINETIC_H3x2z_Py_a-2.0E0*3*I_KINETIC_H3x2z_Py_a+2*1*I_KINETIC_F3x_Py;
  abcd[342] = 4.0E0*I_KINETIC_K2x3y2z_Py_aa-2.0E0*1*I_KINETIC_H2x3y_Py_a;
  abcd[343] = 4.0E0*I_KINETIC_K2x2y3z_Py_aa-2.0E0*1*I_KINETIC_H2x2yz_Py_a-2.0E0*2*I_KINETIC_H2x2yz_Py_a;
  abcd[344] = 4.0E0*I_KINETIC_K2xy4z_Py_aa-2.0E0*2*I_KINETIC_H2xy2z_Py_a-2.0E0*3*I_KINETIC_H2xy2z_Py_a+2*1*I_KINETIC_F2xy_Py;
  abcd[345] = 4.0E0*I_KINETIC_K2x5z_Py_aa-2.0E0*3*I_KINETIC_H2x3z_Py_a-2.0E0*4*I_KINETIC_H2x3z_Py_a+3*2*I_KINETIC_F2xz_Py;
  abcd[346] = 4.0E0*I_KINETIC_Kx4y2z_Py_aa-2.0E0*1*I_KINETIC_Hx4y_Py_a;
  abcd[347] = 4.0E0*I_KINETIC_Kx3y3z_Py_aa-2.0E0*1*I_KINETIC_Hx3yz_Py_a-2.0E0*2*I_KINETIC_Hx3yz_Py_a;
  abcd[348] = 4.0E0*I_KINETIC_Kx2y4z_Py_aa-2.0E0*2*I_KINETIC_Hx2y2z_Py_a-2.0E0*3*I_KINETIC_Hx2y2z_Py_a+2*1*I_KINETIC_Fx2y_Py;
  abcd[349] = 4.0E0*I_KINETIC_Kxy5z_Py_aa-2.0E0*3*I_KINETIC_Hxy3z_Py_a-2.0E0*4*I_KINETIC_Hxy3z_Py_a+3*2*I_KINETIC_Fxyz_Py;
  abcd[350] = 4.0E0*I_KINETIC_Kx6z_Py_aa-2.0E0*4*I_KINETIC_Hx4z_Py_a-2.0E0*5*I_KINETIC_Hx4z_Py_a+4*3*I_KINETIC_Fx2z_Py;
  abcd[351] = 4.0E0*I_KINETIC_K5y2z_Py_aa-2.0E0*1*I_KINETIC_H5y_Py_a;
  abcd[352] = 4.0E0*I_KINETIC_K4y3z_Py_aa-2.0E0*1*I_KINETIC_H4yz_Py_a-2.0E0*2*I_KINETIC_H4yz_Py_a;
  abcd[353] = 4.0E0*I_KINETIC_K3y4z_Py_aa-2.0E0*2*I_KINETIC_H3y2z_Py_a-2.0E0*3*I_KINETIC_H3y2z_Py_a+2*1*I_KINETIC_F3y_Py;
  abcd[354] = 4.0E0*I_KINETIC_K2y5z_Py_aa-2.0E0*3*I_KINETIC_H2y3z_Py_a-2.0E0*4*I_KINETIC_H2y3z_Py_a+3*2*I_KINETIC_F2yz_Py;
  abcd[355] = 4.0E0*I_KINETIC_Ky6z_Py_aa-2.0E0*4*I_KINETIC_Hy4z_Py_a-2.0E0*5*I_KINETIC_Hy4z_Py_a+4*3*I_KINETIC_Fy2z_Py;
  abcd[356] = 4.0E0*I_KINETIC_K7z_Py_aa-2.0E0*5*I_KINETIC_H5z_Py_a-2.0E0*6*I_KINETIC_H5z_Py_a+5*4*I_KINETIC_F3z_Py;
  abcd[357] = 4.0E0*I_KINETIC_K5x2z_Pz_aa-2.0E0*1*I_KINETIC_H5x_Pz_a;
  abcd[358] = 4.0E0*I_KINETIC_K4xy2z_Pz_aa-2.0E0*1*I_KINETIC_H4xy_Pz_a;
  abcd[359] = 4.0E0*I_KINETIC_K4x3z_Pz_aa-2.0E0*1*I_KINETIC_H4xz_Pz_a-2.0E0*2*I_KINETIC_H4xz_Pz_a;
  abcd[360] = 4.0E0*I_KINETIC_K3x2y2z_Pz_aa-2.0E0*1*I_KINETIC_H3x2y_Pz_a;
  abcd[361] = 4.0E0*I_KINETIC_K3xy3z_Pz_aa-2.0E0*1*I_KINETIC_H3xyz_Pz_a-2.0E0*2*I_KINETIC_H3xyz_Pz_a;
  abcd[362] = 4.0E0*I_KINETIC_K3x4z_Pz_aa-2.0E0*2*I_KINETIC_H3x2z_Pz_a-2.0E0*3*I_KINETIC_H3x2z_Pz_a+2*1*I_KINETIC_F3x_Pz;
  abcd[363] = 4.0E0*I_KINETIC_K2x3y2z_Pz_aa-2.0E0*1*I_KINETIC_H2x3y_Pz_a;
  abcd[364] = 4.0E0*I_KINETIC_K2x2y3z_Pz_aa-2.0E0*1*I_KINETIC_H2x2yz_Pz_a-2.0E0*2*I_KINETIC_H2x2yz_Pz_a;
  abcd[365] = 4.0E0*I_KINETIC_K2xy4z_Pz_aa-2.0E0*2*I_KINETIC_H2xy2z_Pz_a-2.0E0*3*I_KINETIC_H2xy2z_Pz_a+2*1*I_KINETIC_F2xy_Pz;
  abcd[366] = 4.0E0*I_KINETIC_K2x5z_Pz_aa-2.0E0*3*I_KINETIC_H2x3z_Pz_a-2.0E0*4*I_KINETIC_H2x3z_Pz_a+3*2*I_KINETIC_F2xz_Pz;
  abcd[367] = 4.0E0*I_KINETIC_Kx4y2z_Pz_aa-2.0E0*1*I_KINETIC_Hx4y_Pz_a;
  abcd[368] = 4.0E0*I_KINETIC_Kx3y3z_Pz_aa-2.0E0*1*I_KINETIC_Hx3yz_Pz_a-2.0E0*2*I_KINETIC_Hx3yz_Pz_a;
  abcd[369] = 4.0E0*I_KINETIC_Kx2y4z_Pz_aa-2.0E0*2*I_KINETIC_Hx2y2z_Pz_a-2.0E0*3*I_KINETIC_Hx2y2z_Pz_a+2*1*I_KINETIC_Fx2y_Pz;
  abcd[370] = 4.0E0*I_KINETIC_Kxy5z_Pz_aa-2.0E0*3*I_KINETIC_Hxy3z_Pz_a-2.0E0*4*I_KINETIC_Hxy3z_Pz_a+3*2*I_KINETIC_Fxyz_Pz;
  abcd[371] = 4.0E0*I_KINETIC_Kx6z_Pz_aa-2.0E0*4*I_KINETIC_Hx4z_Pz_a-2.0E0*5*I_KINETIC_Hx4z_Pz_a+4*3*I_KINETIC_Fx2z_Pz;
  abcd[372] = 4.0E0*I_KINETIC_K5y2z_Pz_aa-2.0E0*1*I_KINETIC_H5y_Pz_a;
  abcd[373] = 4.0E0*I_KINETIC_K4y3z_Pz_aa-2.0E0*1*I_KINETIC_H4yz_Pz_a-2.0E0*2*I_KINETIC_H4yz_Pz_a;
  abcd[374] = 4.0E0*I_KINETIC_K3y4z_Pz_aa-2.0E0*2*I_KINETIC_H3y2z_Pz_a-2.0E0*3*I_KINETIC_H3y2z_Pz_a+2*1*I_KINETIC_F3y_Pz;
  abcd[375] = 4.0E0*I_KINETIC_K2y5z_Pz_aa-2.0E0*3*I_KINETIC_H2y3z_Pz_a-2.0E0*4*I_KINETIC_H2y3z_Pz_a+3*2*I_KINETIC_F2yz_Pz;
  abcd[376] = 4.0E0*I_KINETIC_Ky6z_Pz_aa-2.0E0*4*I_KINETIC_Hy4z_Pz_a-2.0E0*5*I_KINETIC_Hy4z_Pz_a+4*3*I_KINETIC_Fy2z_Pz;
  abcd[377] = 4.0E0*I_KINETIC_K7z_Pz_aa-2.0E0*5*I_KINETIC_H5z_Pz_a-2.0E0*6*I_KINETIC_H5z_Pz_a+5*4*I_KINETIC_F3z_Pz;
}
