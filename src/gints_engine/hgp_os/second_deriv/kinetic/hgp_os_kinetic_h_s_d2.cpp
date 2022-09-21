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
// BRA1 as redundant position, total RHS integrals evaluated as: 4518
// BRA2 as redundant position, total RHS integrals evaluated as: 1806
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

void hgp_os_kinetic_h_s_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_K7x_S_aa = 0.0E0;
  Double I_KINETIC_K6xy_S_aa = 0.0E0;
  Double I_KINETIC_K6xz_S_aa = 0.0E0;
  Double I_KINETIC_K5x2y_S_aa = 0.0E0;
  Double I_KINETIC_K5xyz_S_aa = 0.0E0;
  Double I_KINETIC_K5x2z_S_aa = 0.0E0;
  Double I_KINETIC_K4x3y_S_aa = 0.0E0;
  Double I_KINETIC_K4x2yz_S_aa = 0.0E0;
  Double I_KINETIC_K4xy2z_S_aa = 0.0E0;
  Double I_KINETIC_K4x3z_S_aa = 0.0E0;
  Double I_KINETIC_K3x4y_S_aa = 0.0E0;
  Double I_KINETIC_K3x3yz_S_aa = 0.0E0;
  Double I_KINETIC_K3x2y2z_S_aa = 0.0E0;
  Double I_KINETIC_K3xy3z_S_aa = 0.0E0;
  Double I_KINETIC_K3x4z_S_aa = 0.0E0;
  Double I_KINETIC_K2x5y_S_aa = 0.0E0;
  Double I_KINETIC_K2x4yz_S_aa = 0.0E0;
  Double I_KINETIC_K2x3y2z_S_aa = 0.0E0;
  Double I_KINETIC_K2x2y3z_S_aa = 0.0E0;
  Double I_KINETIC_K2xy4z_S_aa = 0.0E0;
  Double I_KINETIC_K2x5z_S_aa = 0.0E0;
  Double I_KINETIC_Kx6y_S_aa = 0.0E0;
  Double I_KINETIC_Kx5yz_S_aa = 0.0E0;
  Double I_KINETIC_Kx4y2z_S_aa = 0.0E0;
  Double I_KINETIC_Kx3y3z_S_aa = 0.0E0;
  Double I_KINETIC_Kx2y4z_S_aa = 0.0E0;
  Double I_KINETIC_Kxy5z_S_aa = 0.0E0;
  Double I_KINETIC_Kx6z_S_aa = 0.0E0;
  Double I_KINETIC_K7y_S_aa = 0.0E0;
  Double I_KINETIC_K6yz_S_aa = 0.0E0;
  Double I_KINETIC_K5y2z_S_aa = 0.0E0;
  Double I_KINETIC_K4y3z_S_aa = 0.0E0;
  Double I_KINETIC_K3y4z_S_aa = 0.0E0;
  Double I_KINETIC_K2y5z_S_aa = 0.0E0;
  Double I_KINETIC_Ky6z_S_aa = 0.0E0;
  Double I_KINETIC_K7z_S_aa = 0.0E0;
  Double I_KINETIC_H5x_S_a = 0.0E0;
  Double I_KINETIC_H4xy_S_a = 0.0E0;
  Double I_KINETIC_H4xz_S_a = 0.0E0;
  Double I_KINETIC_H3x2y_S_a = 0.0E0;
  Double I_KINETIC_H3xyz_S_a = 0.0E0;
  Double I_KINETIC_H3x2z_S_a = 0.0E0;
  Double I_KINETIC_H2x3y_S_a = 0.0E0;
  Double I_KINETIC_H2x2yz_S_a = 0.0E0;
  Double I_KINETIC_H2xy2z_S_a = 0.0E0;
  Double I_KINETIC_H2x3z_S_a = 0.0E0;
  Double I_KINETIC_Hx4y_S_a = 0.0E0;
  Double I_KINETIC_Hx3yz_S_a = 0.0E0;
  Double I_KINETIC_Hx2y2z_S_a = 0.0E0;
  Double I_KINETIC_Hxy3z_S_a = 0.0E0;
  Double I_KINETIC_Hx4z_S_a = 0.0E0;
  Double I_KINETIC_H5y_S_a = 0.0E0;
  Double I_KINETIC_H4yz_S_a = 0.0E0;
  Double I_KINETIC_H3y2z_S_a = 0.0E0;
  Double I_KINETIC_H2y3z_S_a = 0.0E0;
  Double I_KINETIC_Hy4z_S_a = 0.0E0;
  Double I_KINETIC_H5z_S_a = 0.0E0;
  Double I_KINETIC_F3x_S = 0.0E0;
  Double I_KINETIC_F2xy_S = 0.0E0;
  Double I_KINETIC_F2xz_S = 0.0E0;
  Double I_KINETIC_Fx2y_S = 0.0E0;
  Double I_KINETIC_Fxyz_S = 0.0E0;
  Double I_KINETIC_Fx2z_S = 0.0E0;
  Double I_KINETIC_F3y_S = 0.0E0;
  Double I_KINETIC_F2yz_S = 0.0E0;
  Double I_KINETIC_Fy2z_S = 0.0E0;
  Double I_KINETIC_F3z_S = 0.0E0;

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
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_S_vrr = PAY*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
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
     * totally 7 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_S_vrr = PAX*I_TWOBODYOVERLAP_H5x_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_S_vrr = PAY*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_S_vrr = PAZ*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_S_vrr = PAY*I_TWOBODYOVERLAP_H4xy_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_S_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_S_vrr = PAX*I_TWOBODYOVERLAP_H5y_S_vrr;
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
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_KINETIC_D2x_S_vrr = PAX*I_KINETIC_Px_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_Dxy_S_vrr = PAY*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_KINETIC_D2y_S_vrr = PAY*I_KINETIC_Py_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
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
     * shell quartet name: SQ_KINETIC_I_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 7 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_KINETIC_I6x_S_vrr = PAX*I_KINETIC_H5x_S_vrr+5*oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I5xy_S_vrr = PAY*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_KINETIC_I5xz_S_vrr = PAZ*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_KINETIC_I4x2y_S_vrr = PAY*I_KINETIC_H4xy_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_S_vrr-bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I4x2z_S_vrr = PAZ*I_KINETIC_H4xz_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_S_vrr-bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I3x3y_S_vrr = PAY*I_KINETIC_H3x2y_S_vrr+2*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_I3x2yz_S_vrr = PAZ*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_KINETIC_I3x3z_S_vrr = PAZ*I_KINETIC_H3x2z_S_vrr+2*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_I2x4y_S_vrr = PAX*I_KINETIC_Hx4y_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_S_vrr-bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I2x3yz_S_vrr = PAZ*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_KINETIC_I2xy3z_S_vrr = PAY*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_KINETIC_I2x4z_S_vrr = PAX*I_KINETIC_Hx4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_S_vrr-bdz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_Ix5y_S_vrr = PAX*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_S_vrr;
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
     * shell quartet name: SQ_KINETIC_K_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_K_S_aa_coefs = alpha*alpha;
    I_KINETIC_K7x_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K7x_S_vrr;
    I_KINETIC_K6xy_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K6xy_S_vrr;
    I_KINETIC_K6xz_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K6xz_S_vrr;
    I_KINETIC_K5x2y_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K5x2y_S_vrr;
    I_KINETIC_K5xyz_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K5xyz_S_vrr;
    I_KINETIC_K5x2z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K5x2z_S_vrr;
    I_KINETIC_K4x3y_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K4x3y_S_vrr;
    I_KINETIC_K4x2yz_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K4x2yz_S_vrr;
    I_KINETIC_K4xy2z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K4xy2z_S_vrr;
    I_KINETIC_K4x3z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K4x3z_S_vrr;
    I_KINETIC_K3x4y_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K3x4y_S_vrr;
    I_KINETIC_K3x3yz_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K3x3yz_S_vrr;
    I_KINETIC_K3x2y2z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K3x2y2z_S_vrr;
    I_KINETIC_K3xy3z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K3xy3z_S_vrr;
    I_KINETIC_K3x4z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K3x4z_S_vrr;
    I_KINETIC_K2x5y_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K2x5y_S_vrr;
    I_KINETIC_K2x4yz_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K2x4yz_S_vrr;
    I_KINETIC_K2x3y2z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K2x3y2z_S_vrr;
    I_KINETIC_K2x2y3z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K2x2y3z_S_vrr;
    I_KINETIC_K2xy4z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K2xy4z_S_vrr;
    I_KINETIC_K2x5z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K2x5z_S_vrr;
    I_KINETIC_Kx6y_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_Kx6y_S_vrr;
    I_KINETIC_Kx5yz_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_Kx5yz_S_vrr;
    I_KINETIC_Kx4y2z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_Kx4y2z_S_vrr;
    I_KINETIC_Kx3y3z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_Kx3y3z_S_vrr;
    I_KINETIC_Kx2y4z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_Kx2y4z_S_vrr;
    I_KINETIC_Kxy5z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_Kxy5z_S_vrr;
    I_KINETIC_Kx6z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_Kx6z_S_vrr;
    I_KINETIC_K7y_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K7y_S_vrr;
    I_KINETIC_K6yz_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K6yz_S_vrr;
    I_KINETIC_K5y2z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K5y2z_S_vrr;
    I_KINETIC_K4y3z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K4y3z_S_vrr;
    I_KINETIC_K3y4z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K3y4z_S_vrr;
    I_KINETIC_K2y5z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K2y5z_S_vrr;
    I_KINETIC_Ky6z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_Ky6z_S_vrr;
    I_KINETIC_K7z_S_aa += SQ_KINETIC_K_S_aa_coefs*I_KINETIC_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_H_S_a_coefs = alpha;
    I_KINETIC_H5x_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H5x_S_vrr;
    I_KINETIC_H4xy_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H4xy_S_vrr;
    I_KINETIC_H4xz_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H4xz_S_vrr;
    I_KINETIC_H3x2y_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H3x2y_S_vrr;
    I_KINETIC_H3xyz_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H3xyz_S_vrr;
    I_KINETIC_H3x2z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H3x2z_S_vrr;
    I_KINETIC_H2x3y_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H2x3y_S_vrr;
    I_KINETIC_H2x2yz_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H2x2yz_S_vrr;
    I_KINETIC_H2xy2z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H2xy2z_S_vrr;
    I_KINETIC_H2x3z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H2x3z_S_vrr;
    I_KINETIC_Hx4y_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_Hx4y_S_vrr;
    I_KINETIC_Hx3yz_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_Hx3yz_S_vrr;
    I_KINETIC_Hx2y2z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_Hx2y2z_S_vrr;
    I_KINETIC_Hxy3z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_Hxy3z_S_vrr;
    I_KINETIC_Hx4z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_Hx4z_S_vrr;
    I_KINETIC_H5y_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H5y_S_vrr;
    I_KINETIC_H4yz_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H4yz_S_vrr;
    I_KINETIC_H3y2z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H3y2z_S_vrr;
    I_KINETIC_H2y3z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H2y3z_S_vrr;
    I_KINETIC_Hy4z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_Hy4z_S_vrr;
    I_KINETIC_H5z_S_a += SQ_KINETIC_H_S_a_coefs*I_KINETIC_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_F3x_S += I_KINETIC_F3x_S_vrr;
    I_KINETIC_F2xy_S += I_KINETIC_F2xy_S_vrr;
    I_KINETIC_F2xz_S += I_KINETIC_F2xz_S_vrr;
    I_KINETIC_Fx2y_S += I_KINETIC_Fx2y_S_vrr;
    I_KINETIC_Fxyz_S += I_KINETIC_Fxyz_S_vrr;
    I_KINETIC_Fx2z_S += I_KINETIC_Fx2z_S_vrr;
    I_KINETIC_F3y_S += I_KINETIC_F3y_S_vrr;
    I_KINETIC_F2yz_S += I_KINETIC_F2yz_S_vrr;
    I_KINETIC_Fy2z_S += I_KINETIC_Fy2z_S_vrr;
    I_KINETIC_F3z_S += I_KINETIC_F3z_S_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_S_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_S_aa
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_F_S
   ************************************************************/
  abcd[0] = 4.0E0*I_KINETIC_K7x_S_aa-2.0E0*5*I_KINETIC_H5x_S_a-2.0E0*6*I_KINETIC_H5x_S_a+5*4*I_KINETIC_F3x_S;
  abcd[1] = 4.0E0*I_KINETIC_K6xy_S_aa-2.0E0*4*I_KINETIC_H4xy_S_a-2.0E0*5*I_KINETIC_H4xy_S_a+4*3*I_KINETIC_F2xy_S;
  abcd[2] = 4.0E0*I_KINETIC_K6xz_S_aa-2.0E0*4*I_KINETIC_H4xz_S_a-2.0E0*5*I_KINETIC_H4xz_S_a+4*3*I_KINETIC_F2xz_S;
  abcd[3] = 4.0E0*I_KINETIC_K5x2y_S_aa-2.0E0*3*I_KINETIC_H3x2y_S_a-2.0E0*4*I_KINETIC_H3x2y_S_a+3*2*I_KINETIC_Fx2y_S;
  abcd[4] = 4.0E0*I_KINETIC_K5xyz_S_aa-2.0E0*3*I_KINETIC_H3xyz_S_a-2.0E0*4*I_KINETIC_H3xyz_S_a+3*2*I_KINETIC_Fxyz_S;
  abcd[5] = 4.0E0*I_KINETIC_K5x2z_S_aa-2.0E0*3*I_KINETIC_H3x2z_S_a-2.0E0*4*I_KINETIC_H3x2z_S_a+3*2*I_KINETIC_Fx2z_S;
  abcd[6] = 4.0E0*I_KINETIC_K4x3y_S_aa-2.0E0*2*I_KINETIC_H2x3y_S_a-2.0E0*3*I_KINETIC_H2x3y_S_a+2*1*I_KINETIC_F3y_S;
  abcd[7] = 4.0E0*I_KINETIC_K4x2yz_S_aa-2.0E0*2*I_KINETIC_H2x2yz_S_a-2.0E0*3*I_KINETIC_H2x2yz_S_a+2*1*I_KINETIC_F2yz_S;
  abcd[8] = 4.0E0*I_KINETIC_K4xy2z_S_aa-2.0E0*2*I_KINETIC_H2xy2z_S_a-2.0E0*3*I_KINETIC_H2xy2z_S_a+2*1*I_KINETIC_Fy2z_S;
  abcd[9] = 4.0E0*I_KINETIC_K4x3z_S_aa-2.0E0*2*I_KINETIC_H2x3z_S_a-2.0E0*3*I_KINETIC_H2x3z_S_a+2*1*I_KINETIC_F3z_S;
  abcd[10] = 4.0E0*I_KINETIC_K3x4y_S_aa-2.0E0*1*I_KINETIC_Hx4y_S_a-2.0E0*2*I_KINETIC_Hx4y_S_a;
  abcd[11] = 4.0E0*I_KINETIC_K3x3yz_S_aa-2.0E0*1*I_KINETIC_Hx3yz_S_a-2.0E0*2*I_KINETIC_Hx3yz_S_a;
  abcd[12] = 4.0E0*I_KINETIC_K3x2y2z_S_aa-2.0E0*1*I_KINETIC_Hx2y2z_S_a-2.0E0*2*I_KINETIC_Hx2y2z_S_a;
  abcd[13] = 4.0E0*I_KINETIC_K3xy3z_S_aa-2.0E0*1*I_KINETIC_Hxy3z_S_a-2.0E0*2*I_KINETIC_Hxy3z_S_a;
  abcd[14] = 4.0E0*I_KINETIC_K3x4z_S_aa-2.0E0*1*I_KINETIC_Hx4z_S_a-2.0E0*2*I_KINETIC_Hx4z_S_a;
  abcd[15] = 4.0E0*I_KINETIC_K2x5y_S_aa-2.0E0*1*I_KINETIC_H5y_S_a;
  abcd[16] = 4.0E0*I_KINETIC_K2x4yz_S_aa-2.0E0*1*I_KINETIC_H4yz_S_a;
  abcd[17] = 4.0E0*I_KINETIC_K2x3y2z_S_aa-2.0E0*1*I_KINETIC_H3y2z_S_a;
  abcd[18] = 4.0E0*I_KINETIC_K2x2y3z_S_aa-2.0E0*1*I_KINETIC_H2y3z_S_a;
  abcd[19] = 4.0E0*I_KINETIC_K2xy4z_S_aa-2.0E0*1*I_KINETIC_Hy4z_S_a;
  abcd[20] = 4.0E0*I_KINETIC_K2x5z_S_aa-2.0E0*1*I_KINETIC_H5z_S_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_S_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_S_aa
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_F_S
   ************************************************************/
  abcd[21] = 4.0E0*I_KINETIC_K6xy_S_aa-2.0E0*5*I_KINETIC_H4xy_S_a;
  abcd[22] = 4.0E0*I_KINETIC_K5x2y_S_aa-2.0E0*1*I_KINETIC_H5x_S_a-2.0E0*4*I_KINETIC_H3x2y_S_a+4*1*I_KINETIC_F3x_S;
  abcd[23] = 4.0E0*I_KINETIC_K5xyz_S_aa-2.0E0*4*I_KINETIC_H3xyz_S_a;
  abcd[24] = 4.0E0*I_KINETIC_K4x3y_S_aa-2.0E0*2*I_KINETIC_H4xy_S_a-2.0E0*3*I_KINETIC_H2x3y_S_a+3*2*I_KINETIC_F2xy_S;
  abcd[25] = 4.0E0*I_KINETIC_K4x2yz_S_aa-2.0E0*1*I_KINETIC_H4xz_S_a-2.0E0*3*I_KINETIC_H2x2yz_S_a+3*1*I_KINETIC_F2xz_S;
  abcd[26] = 4.0E0*I_KINETIC_K4xy2z_S_aa-2.0E0*3*I_KINETIC_H2xy2z_S_a;
  abcd[27] = 4.0E0*I_KINETIC_K3x4y_S_aa-2.0E0*3*I_KINETIC_H3x2y_S_a-2.0E0*2*I_KINETIC_Hx4y_S_a+2*3*I_KINETIC_Fx2y_S;
  abcd[28] = 4.0E0*I_KINETIC_K3x3yz_S_aa-2.0E0*2*I_KINETIC_H3xyz_S_a-2.0E0*2*I_KINETIC_Hx3yz_S_a+2*2*I_KINETIC_Fxyz_S;
  abcd[29] = 4.0E0*I_KINETIC_K3x2y2z_S_aa-2.0E0*1*I_KINETIC_H3x2z_S_a-2.0E0*2*I_KINETIC_Hx2y2z_S_a+2*1*I_KINETIC_Fx2z_S;
  abcd[30] = 4.0E0*I_KINETIC_K3xy3z_S_aa-2.0E0*2*I_KINETIC_Hxy3z_S_a;
  abcd[31] = 4.0E0*I_KINETIC_K2x5y_S_aa-2.0E0*4*I_KINETIC_H2x3y_S_a-2.0E0*1*I_KINETIC_H5y_S_a+4*I_KINETIC_F3y_S;
  abcd[32] = 4.0E0*I_KINETIC_K2x4yz_S_aa-2.0E0*3*I_KINETIC_H2x2yz_S_a-2.0E0*1*I_KINETIC_H4yz_S_a+3*I_KINETIC_F2yz_S;
  abcd[33] = 4.0E0*I_KINETIC_K2x3y2z_S_aa-2.0E0*2*I_KINETIC_H2xy2z_S_a-2.0E0*1*I_KINETIC_H3y2z_S_a+2*I_KINETIC_Fy2z_S;
  abcd[34] = 4.0E0*I_KINETIC_K2x2y3z_S_aa-2.0E0*1*I_KINETIC_H2x3z_S_a-2.0E0*1*I_KINETIC_H2y3z_S_a+1*I_KINETIC_F3z_S;
  abcd[35] = 4.0E0*I_KINETIC_K2xy4z_S_aa-2.0E0*1*I_KINETIC_Hy4z_S_a;
  abcd[36] = 4.0E0*I_KINETIC_Kx6y_S_aa-2.0E0*5*I_KINETIC_Hx4y_S_a;
  abcd[37] = 4.0E0*I_KINETIC_Kx5yz_S_aa-2.0E0*4*I_KINETIC_Hx3yz_S_a;
  abcd[38] = 4.0E0*I_KINETIC_Kx4y2z_S_aa-2.0E0*3*I_KINETIC_Hx2y2z_S_a;
  abcd[39] = 4.0E0*I_KINETIC_Kx3y3z_S_aa-2.0E0*2*I_KINETIC_Hxy3z_S_a;
  abcd[40] = 4.0E0*I_KINETIC_Kx2y4z_S_aa-2.0E0*1*I_KINETIC_Hx4z_S_a;
  abcd[41] = 4.0E0*I_KINETIC_Kxy5z_S_aa;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_S_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_S_aa
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_F_S
   ************************************************************/
  abcd[42] = 4.0E0*I_KINETIC_K6xz_S_aa-2.0E0*5*I_KINETIC_H4xz_S_a;
  abcd[43] = 4.0E0*I_KINETIC_K5xyz_S_aa-2.0E0*4*I_KINETIC_H3xyz_S_a;
  abcd[44] = 4.0E0*I_KINETIC_K5x2z_S_aa-2.0E0*1*I_KINETIC_H5x_S_a-2.0E0*4*I_KINETIC_H3x2z_S_a+4*1*I_KINETIC_F3x_S;
  abcd[45] = 4.0E0*I_KINETIC_K4x2yz_S_aa-2.0E0*3*I_KINETIC_H2x2yz_S_a;
  abcd[46] = 4.0E0*I_KINETIC_K4xy2z_S_aa-2.0E0*1*I_KINETIC_H4xy_S_a-2.0E0*3*I_KINETIC_H2xy2z_S_a+3*1*I_KINETIC_F2xy_S;
  abcd[47] = 4.0E0*I_KINETIC_K4x3z_S_aa-2.0E0*2*I_KINETIC_H4xz_S_a-2.0E0*3*I_KINETIC_H2x3z_S_a+3*2*I_KINETIC_F2xz_S;
  abcd[48] = 4.0E0*I_KINETIC_K3x3yz_S_aa-2.0E0*2*I_KINETIC_Hx3yz_S_a;
  abcd[49] = 4.0E0*I_KINETIC_K3x2y2z_S_aa-2.0E0*1*I_KINETIC_H3x2y_S_a-2.0E0*2*I_KINETIC_Hx2y2z_S_a+2*1*I_KINETIC_Fx2y_S;
  abcd[50] = 4.0E0*I_KINETIC_K3xy3z_S_aa-2.0E0*2*I_KINETIC_H3xyz_S_a-2.0E0*2*I_KINETIC_Hxy3z_S_a+2*2*I_KINETIC_Fxyz_S;
  abcd[51] = 4.0E0*I_KINETIC_K3x4z_S_aa-2.0E0*3*I_KINETIC_H3x2z_S_a-2.0E0*2*I_KINETIC_Hx4z_S_a+2*3*I_KINETIC_Fx2z_S;
  abcd[52] = 4.0E0*I_KINETIC_K2x4yz_S_aa-2.0E0*1*I_KINETIC_H4yz_S_a;
  abcd[53] = 4.0E0*I_KINETIC_K2x3y2z_S_aa-2.0E0*1*I_KINETIC_H2x3y_S_a-2.0E0*1*I_KINETIC_H3y2z_S_a+1*I_KINETIC_F3y_S;
  abcd[54] = 4.0E0*I_KINETIC_K2x2y3z_S_aa-2.0E0*2*I_KINETIC_H2x2yz_S_a-2.0E0*1*I_KINETIC_H2y3z_S_a+2*I_KINETIC_F2yz_S;
  abcd[55] = 4.0E0*I_KINETIC_K2xy4z_S_aa-2.0E0*3*I_KINETIC_H2xy2z_S_a-2.0E0*1*I_KINETIC_Hy4z_S_a+3*I_KINETIC_Fy2z_S;
  abcd[56] = 4.0E0*I_KINETIC_K2x5z_S_aa-2.0E0*4*I_KINETIC_H2x3z_S_a-2.0E0*1*I_KINETIC_H5z_S_a+4*I_KINETIC_F3z_S;
  abcd[57] = 4.0E0*I_KINETIC_Kx5yz_S_aa;
  abcd[58] = 4.0E0*I_KINETIC_Kx4y2z_S_aa-2.0E0*1*I_KINETIC_Hx4y_S_a;
  abcd[59] = 4.0E0*I_KINETIC_Kx3y3z_S_aa-2.0E0*2*I_KINETIC_Hx3yz_S_a;
  abcd[60] = 4.0E0*I_KINETIC_Kx2y4z_S_aa-2.0E0*3*I_KINETIC_Hx2y2z_S_a;
  abcd[61] = 4.0E0*I_KINETIC_Kxy5z_S_aa-2.0E0*4*I_KINETIC_Hxy3z_S_a;
  abcd[62] = 4.0E0*I_KINETIC_Kx6z_S_aa-2.0E0*5*I_KINETIC_Hx4z_S_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_S_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_S_aa
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_F_S
   ************************************************************/
  abcd[63] = 4.0E0*I_KINETIC_K5x2y_S_aa-2.0E0*1*I_KINETIC_H5x_S_a;
  abcd[64] = 4.0E0*I_KINETIC_K4x3y_S_aa-2.0E0*1*I_KINETIC_H4xy_S_a-2.0E0*2*I_KINETIC_H4xy_S_a;
  abcd[65] = 4.0E0*I_KINETIC_K4x2yz_S_aa-2.0E0*1*I_KINETIC_H4xz_S_a;
  abcd[66] = 4.0E0*I_KINETIC_K3x4y_S_aa-2.0E0*2*I_KINETIC_H3x2y_S_a-2.0E0*3*I_KINETIC_H3x2y_S_a+2*1*I_KINETIC_F3x_S;
  abcd[67] = 4.0E0*I_KINETIC_K3x3yz_S_aa-2.0E0*1*I_KINETIC_H3xyz_S_a-2.0E0*2*I_KINETIC_H3xyz_S_a;
  abcd[68] = 4.0E0*I_KINETIC_K3x2y2z_S_aa-2.0E0*1*I_KINETIC_H3x2z_S_a;
  abcd[69] = 4.0E0*I_KINETIC_K2x5y_S_aa-2.0E0*3*I_KINETIC_H2x3y_S_a-2.0E0*4*I_KINETIC_H2x3y_S_a+3*2*I_KINETIC_F2xy_S;
  abcd[70] = 4.0E0*I_KINETIC_K2x4yz_S_aa-2.0E0*2*I_KINETIC_H2x2yz_S_a-2.0E0*3*I_KINETIC_H2x2yz_S_a+2*1*I_KINETIC_F2xz_S;
  abcd[71] = 4.0E0*I_KINETIC_K2x3y2z_S_aa-2.0E0*1*I_KINETIC_H2xy2z_S_a-2.0E0*2*I_KINETIC_H2xy2z_S_a;
  abcd[72] = 4.0E0*I_KINETIC_K2x2y3z_S_aa-2.0E0*1*I_KINETIC_H2x3z_S_a;
  abcd[73] = 4.0E0*I_KINETIC_Kx6y_S_aa-2.0E0*4*I_KINETIC_Hx4y_S_a-2.0E0*5*I_KINETIC_Hx4y_S_a+4*3*I_KINETIC_Fx2y_S;
  abcd[74] = 4.0E0*I_KINETIC_Kx5yz_S_aa-2.0E0*3*I_KINETIC_Hx3yz_S_a-2.0E0*4*I_KINETIC_Hx3yz_S_a+3*2*I_KINETIC_Fxyz_S;
  abcd[75] = 4.0E0*I_KINETIC_Kx4y2z_S_aa-2.0E0*2*I_KINETIC_Hx2y2z_S_a-2.0E0*3*I_KINETIC_Hx2y2z_S_a+2*1*I_KINETIC_Fx2z_S;
  abcd[76] = 4.0E0*I_KINETIC_Kx3y3z_S_aa-2.0E0*1*I_KINETIC_Hxy3z_S_a-2.0E0*2*I_KINETIC_Hxy3z_S_a;
  abcd[77] = 4.0E0*I_KINETIC_Kx2y4z_S_aa-2.0E0*1*I_KINETIC_Hx4z_S_a;
  abcd[78] = 4.0E0*I_KINETIC_K7y_S_aa-2.0E0*5*I_KINETIC_H5y_S_a-2.0E0*6*I_KINETIC_H5y_S_a+5*4*I_KINETIC_F3y_S;
  abcd[79] = 4.0E0*I_KINETIC_K6yz_S_aa-2.0E0*4*I_KINETIC_H4yz_S_a-2.0E0*5*I_KINETIC_H4yz_S_a+4*3*I_KINETIC_F2yz_S;
  abcd[80] = 4.0E0*I_KINETIC_K5y2z_S_aa-2.0E0*3*I_KINETIC_H3y2z_S_a-2.0E0*4*I_KINETIC_H3y2z_S_a+3*2*I_KINETIC_Fy2z_S;
  abcd[81] = 4.0E0*I_KINETIC_K4y3z_S_aa-2.0E0*2*I_KINETIC_H2y3z_S_a-2.0E0*3*I_KINETIC_H2y3z_S_a+2*1*I_KINETIC_F3z_S;
  abcd[82] = 4.0E0*I_KINETIC_K3y4z_S_aa-2.0E0*1*I_KINETIC_Hy4z_S_a-2.0E0*2*I_KINETIC_Hy4z_S_a;
  abcd[83] = 4.0E0*I_KINETIC_K2y5z_S_aa-2.0E0*1*I_KINETIC_H5z_S_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_S_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_S_aa
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_F_S
   ************************************************************/
  abcd[84] = 4.0E0*I_KINETIC_K5xyz_S_aa;
  abcd[85] = 4.0E0*I_KINETIC_K4x2yz_S_aa-2.0E0*1*I_KINETIC_H4xz_S_a;
  abcd[86] = 4.0E0*I_KINETIC_K4xy2z_S_aa-2.0E0*1*I_KINETIC_H4xy_S_a;
  abcd[87] = 4.0E0*I_KINETIC_K3x3yz_S_aa-2.0E0*2*I_KINETIC_H3xyz_S_a;
  abcd[88] = 4.0E0*I_KINETIC_K3x2y2z_S_aa-2.0E0*1*I_KINETIC_H3x2y_S_a-2.0E0*1*I_KINETIC_H3x2z_S_a+1*I_KINETIC_F3x_S;
  abcd[89] = 4.0E0*I_KINETIC_K3xy3z_S_aa-2.0E0*2*I_KINETIC_H3xyz_S_a;
  abcd[90] = 4.0E0*I_KINETIC_K2x4yz_S_aa-2.0E0*3*I_KINETIC_H2x2yz_S_a;
  abcd[91] = 4.0E0*I_KINETIC_K2x3y2z_S_aa-2.0E0*1*I_KINETIC_H2x3y_S_a-2.0E0*2*I_KINETIC_H2xy2z_S_a+2*1*I_KINETIC_F2xy_S;
  abcd[92] = 4.0E0*I_KINETIC_K2x2y3z_S_aa-2.0E0*2*I_KINETIC_H2x2yz_S_a-2.0E0*1*I_KINETIC_H2x3z_S_a+2*I_KINETIC_F2xz_S;
  abcd[93] = 4.0E0*I_KINETIC_K2xy4z_S_aa-2.0E0*3*I_KINETIC_H2xy2z_S_a;
  abcd[94] = 4.0E0*I_KINETIC_Kx5yz_S_aa-2.0E0*4*I_KINETIC_Hx3yz_S_a;
  abcd[95] = 4.0E0*I_KINETIC_Kx4y2z_S_aa-2.0E0*1*I_KINETIC_Hx4y_S_a-2.0E0*3*I_KINETIC_Hx2y2z_S_a+3*1*I_KINETIC_Fx2y_S;
  abcd[96] = 4.0E0*I_KINETIC_Kx3y3z_S_aa-2.0E0*2*I_KINETIC_Hx3yz_S_a-2.0E0*2*I_KINETIC_Hxy3z_S_a+2*2*I_KINETIC_Fxyz_S;
  abcd[97] = 4.0E0*I_KINETIC_Kx2y4z_S_aa-2.0E0*3*I_KINETIC_Hx2y2z_S_a-2.0E0*1*I_KINETIC_Hx4z_S_a+3*I_KINETIC_Fx2z_S;
  abcd[98] = 4.0E0*I_KINETIC_Kxy5z_S_aa-2.0E0*4*I_KINETIC_Hxy3z_S_a;
  abcd[99] = 4.0E0*I_KINETIC_K6yz_S_aa-2.0E0*5*I_KINETIC_H4yz_S_a;
  abcd[100] = 4.0E0*I_KINETIC_K5y2z_S_aa-2.0E0*1*I_KINETIC_H5y_S_a-2.0E0*4*I_KINETIC_H3y2z_S_a+4*1*I_KINETIC_F3y_S;
  abcd[101] = 4.0E0*I_KINETIC_K4y3z_S_aa-2.0E0*2*I_KINETIC_H4yz_S_a-2.0E0*3*I_KINETIC_H2y3z_S_a+3*2*I_KINETIC_F2yz_S;
  abcd[102] = 4.0E0*I_KINETIC_K3y4z_S_aa-2.0E0*3*I_KINETIC_H3y2z_S_a-2.0E0*2*I_KINETIC_Hy4z_S_a+2*3*I_KINETIC_Fy2z_S;
  abcd[103] = 4.0E0*I_KINETIC_K2y5z_S_aa-2.0E0*4*I_KINETIC_H2y3z_S_a-2.0E0*1*I_KINETIC_H5z_S_a+4*I_KINETIC_F3z_S;
  abcd[104] = 4.0E0*I_KINETIC_Ky6z_S_aa-2.0E0*5*I_KINETIC_Hy4z_S_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_S_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_K_S_aa
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_H_S_a
   * RHS shell quartet name: SQ_KINETIC_F_S
   ************************************************************/
  abcd[105] = 4.0E0*I_KINETIC_K5x2z_S_aa-2.0E0*1*I_KINETIC_H5x_S_a;
  abcd[106] = 4.0E0*I_KINETIC_K4xy2z_S_aa-2.0E0*1*I_KINETIC_H4xy_S_a;
  abcd[107] = 4.0E0*I_KINETIC_K4x3z_S_aa-2.0E0*1*I_KINETIC_H4xz_S_a-2.0E0*2*I_KINETIC_H4xz_S_a;
  abcd[108] = 4.0E0*I_KINETIC_K3x2y2z_S_aa-2.0E0*1*I_KINETIC_H3x2y_S_a;
  abcd[109] = 4.0E0*I_KINETIC_K3xy3z_S_aa-2.0E0*1*I_KINETIC_H3xyz_S_a-2.0E0*2*I_KINETIC_H3xyz_S_a;
  abcd[110] = 4.0E0*I_KINETIC_K3x4z_S_aa-2.0E0*2*I_KINETIC_H3x2z_S_a-2.0E0*3*I_KINETIC_H3x2z_S_a+2*1*I_KINETIC_F3x_S;
  abcd[111] = 4.0E0*I_KINETIC_K2x3y2z_S_aa-2.0E0*1*I_KINETIC_H2x3y_S_a;
  abcd[112] = 4.0E0*I_KINETIC_K2x2y3z_S_aa-2.0E0*1*I_KINETIC_H2x2yz_S_a-2.0E0*2*I_KINETIC_H2x2yz_S_a;
  abcd[113] = 4.0E0*I_KINETIC_K2xy4z_S_aa-2.0E0*2*I_KINETIC_H2xy2z_S_a-2.0E0*3*I_KINETIC_H2xy2z_S_a+2*1*I_KINETIC_F2xy_S;
  abcd[114] = 4.0E0*I_KINETIC_K2x5z_S_aa-2.0E0*3*I_KINETIC_H2x3z_S_a-2.0E0*4*I_KINETIC_H2x3z_S_a+3*2*I_KINETIC_F2xz_S;
  abcd[115] = 4.0E0*I_KINETIC_Kx4y2z_S_aa-2.0E0*1*I_KINETIC_Hx4y_S_a;
  abcd[116] = 4.0E0*I_KINETIC_Kx3y3z_S_aa-2.0E0*1*I_KINETIC_Hx3yz_S_a-2.0E0*2*I_KINETIC_Hx3yz_S_a;
  abcd[117] = 4.0E0*I_KINETIC_Kx2y4z_S_aa-2.0E0*2*I_KINETIC_Hx2y2z_S_a-2.0E0*3*I_KINETIC_Hx2y2z_S_a+2*1*I_KINETIC_Fx2y_S;
  abcd[118] = 4.0E0*I_KINETIC_Kxy5z_S_aa-2.0E0*3*I_KINETIC_Hxy3z_S_a-2.0E0*4*I_KINETIC_Hxy3z_S_a+3*2*I_KINETIC_Fxyz_S;
  abcd[119] = 4.0E0*I_KINETIC_Kx6z_S_aa-2.0E0*4*I_KINETIC_Hx4z_S_a-2.0E0*5*I_KINETIC_Hx4z_S_a+4*3*I_KINETIC_Fx2z_S;
  abcd[120] = 4.0E0*I_KINETIC_K5y2z_S_aa-2.0E0*1*I_KINETIC_H5y_S_a;
  abcd[121] = 4.0E0*I_KINETIC_K4y3z_S_aa-2.0E0*1*I_KINETIC_H4yz_S_a-2.0E0*2*I_KINETIC_H4yz_S_a;
  abcd[122] = 4.0E0*I_KINETIC_K3y4z_S_aa-2.0E0*2*I_KINETIC_H3y2z_S_a-2.0E0*3*I_KINETIC_H3y2z_S_a+2*1*I_KINETIC_F3y_S;
  abcd[123] = 4.0E0*I_KINETIC_K2y5z_S_aa-2.0E0*3*I_KINETIC_H2y3z_S_a-2.0E0*4*I_KINETIC_H2y3z_S_a+3*2*I_KINETIC_F2yz_S;
  abcd[124] = 4.0E0*I_KINETIC_Ky6z_S_aa-2.0E0*4*I_KINETIC_Hy4z_S_a-2.0E0*5*I_KINETIC_Hy4z_S_a+4*3*I_KINETIC_Fy2z_S;
  abcd[125] = 4.0E0*I_KINETIC_K7z_S_aa-2.0E0*5*I_KINETIC_H5z_S_a-2.0E0*6*I_KINETIC_H5z_S_a+5*4*I_KINETIC_F3z_S;
}
