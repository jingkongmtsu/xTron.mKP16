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
// BRA1 as redundant position, total RHS integrals evaluated as: 2377
// BRA2 as redundant position, total RHS integrals evaluated as: 2379
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: BRA1
//

//
// @@@@ derivative position-direction information
// BRA2
// X
// Y
// Z
// ####

void hgp_os_kinetic_g_d_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_G4x_F3x_b = 0.0E0;
  Double I_KINETIC_G3xy_F3x_b = 0.0E0;
  Double I_KINETIC_G3xz_F3x_b = 0.0E0;
  Double I_KINETIC_G2x2y_F3x_b = 0.0E0;
  Double I_KINETIC_G2xyz_F3x_b = 0.0E0;
  Double I_KINETIC_G2x2z_F3x_b = 0.0E0;
  Double I_KINETIC_Gx3y_F3x_b = 0.0E0;
  Double I_KINETIC_Gx2yz_F3x_b = 0.0E0;
  Double I_KINETIC_Gxy2z_F3x_b = 0.0E0;
  Double I_KINETIC_Gx3z_F3x_b = 0.0E0;
  Double I_KINETIC_G4y_F3x_b = 0.0E0;
  Double I_KINETIC_G3yz_F3x_b = 0.0E0;
  Double I_KINETIC_G2y2z_F3x_b = 0.0E0;
  Double I_KINETIC_Gy3z_F3x_b = 0.0E0;
  Double I_KINETIC_G4z_F3x_b = 0.0E0;
  Double I_KINETIC_G4x_F2xy_b = 0.0E0;
  Double I_KINETIC_G3xy_F2xy_b = 0.0E0;
  Double I_KINETIC_G3xz_F2xy_b = 0.0E0;
  Double I_KINETIC_G2x2y_F2xy_b = 0.0E0;
  Double I_KINETIC_G2xyz_F2xy_b = 0.0E0;
  Double I_KINETIC_G2x2z_F2xy_b = 0.0E0;
  Double I_KINETIC_Gx3y_F2xy_b = 0.0E0;
  Double I_KINETIC_Gx2yz_F2xy_b = 0.0E0;
  Double I_KINETIC_Gxy2z_F2xy_b = 0.0E0;
  Double I_KINETIC_Gx3z_F2xy_b = 0.0E0;
  Double I_KINETIC_G4y_F2xy_b = 0.0E0;
  Double I_KINETIC_G3yz_F2xy_b = 0.0E0;
  Double I_KINETIC_G2y2z_F2xy_b = 0.0E0;
  Double I_KINETIC_Gy3z_F2xy_b = 0.0E0;
  Double I_KINETIC_G4z_F2xy_b = 0.0E0;
  Double I_KINETIC_G4x_F2xz_b = 0.0E0;
  Double I_KINETIC_G3xy_F2xz_b = 0.0E0;
  Double I_KINETIC_G3xz_F2xz_b = 0.0E0;
  Double I_KINETIC_G2x2y_F2xz_b = 0.0E0;
  Double I_KINETIC_G2xyz_F2xz_b = 0.0E0;
  Double I_KINETIC_G2x2z_F2xz_b = 0.0E0;
  Double I_KINETIC_Gx3y_F2xz_b = 0.0E0;
  Double I_KINETIC_Gx2yz_F2xz_b = 0.0E0;
  Double I_KINETIC_Gxy2z_F2xz_b = 0.0E0;
  Double I_KINETIC_Gx3z_F2xz_b = 0.0E0;
  Double I_KINETIC_G4y_F2xz_b = 0.0E0;
  Double I_KINETIC_G3yz_F2xz_b = 0.0E0;
  Double I_KINETIC_G2y2z_F2xz_b = 0.0E0;
  Double I_KINETIC_Gy3z_F2xz_b = 0.0E0;
  Double I_KINETIC_G4z_F2xz_b = 0.0E0;
  Double I_KINETIC_G4x_Fx2y_b = 0.0E0;
  Double I_KINETIC_G3xy_Fx2y_b = 0.0E0;
  Double I_KINETIC_G3xz_Fx2y_b = 0.0E0;
  Double I_KINETIC_G2x2y_Fx2y_b = 0.0E0;
  Double I_KINETIC_G2xyz_Fx2y_b = 0.0E0;
  Double I_KINETIC_G2x2z_Fx2y_b = 0.0E0;
  Double I_KINETIC_Gx3y_Fx2y_b = 0.0E0;
  Double I_KINETIC_Gx2yz_Fx2y_b = 0.0E0;
  Double I_KINETIC_Gxy2z_Fx2y_b = 0.0E0;
  Double I_KINETIC_Gx3z_Fx2y_b = 0.0E0;
  Double I_KINETIC_G4y_Fx2y_b = 0.0E0;
  Double I_KINETIC_G3yz_Fx2y_b = 0.0E0;
  Double I_KINETIC_G2y2z_Fx2y_b = 0.0E0;
  Double I_KINETIC_Gy3z_Fx2y_b = 0.0E0;
  Double I_KINETIC_G4z_Fx2y_b = 0.0E0;
  Double I_KINETIC_G4x_Fxyz_b = 0.0E0;
  Double I_KINETIC_G3xy_Fxyz_b = 0.0E0;
  Double I_KINETIC_G3xz_Fxyz_b = 0.0E0;
  Double I_KINETIC_G2x2y_Fxyz_b = 0.0E0;
  Double I_KINETIC_G2xyz_Fxyz_b = 0.0E0;
  Double I_KINETIC_G2x2z_Fxyz_b = 0.0E0;
  Double I_KINETIC_Gx3y_Fxyz_b = 0.0E0;
  Double I_KINETIC_Gx2yz_Fxyz_b = 0.0E0;
  Double I_KINETIC_Gxy2z_Fxyz_b = 0.0E0;
  Double I_KINETIC_Gx3z_Fxyz_b = 0.0E0;
  Double I_KINETIC_G4y_Fxyz_b = 0.0E0;
  Double I_KINETIC_G3yz_Fxyz_b = 0.0E0;
  Double I_KINETIC_G2y2z_Fxyz_b = 0.0E0;
  Double I_KINETIC_Gy3z_Fxyz_b = 0.0E0;
  Double I_KINETIC_G4z_Fxyz_b = 0.0E0;
  Double I_KINETIC_G4x_Fx2z_b = 0.0E0;
  Double I_KINETIC_G3xy_Fx2z_b = 0.0E0;
  Double I_KINETIC_G3xz_Fx2z_b = 0.0E0;
  Double I_KINETIC_G2x2y_Fx2z_b = 0.0E0;
  Double I_KINETIC_G2xyz_Fx2z_b = 0.0E0;
  Double I_KINETIC_G2x2z_Fx2z_b = 0.0E0;
  Double I_KINETIC_Gx3y_Fx2z_b = 0.0E0;
  Double I_KINETIC_Gx2yz_Fx2z_b = 0.0E0;
  Double I_KINETIC_Gxy2z_Fx2z_b = 0.0E0;
  Double I_KINETIC_Gx3z_Fx2z_b = 0.0E0;
  Double I_KINETIC_G4y_Fx2z_b = 0.0E0;
  Double I_KINETIC_G3yz_Fx2z_b = 0.0E0;
  Double I_KINETIC_G2y2z_Fx2z_b = 0.0E0;
  Double I_KINETIC_Gy3z_Fx2z_b = 0.0E0;
  Double I_KINETIC_G4z_Fx2z_b = 0.0E0;
  Double I_KINETIC_G4x_F3y_b = 0.0E0;
  Double I_KINETIC_G3xy_F3y_b = 0.0E0;
  Double I_KINETIC_G3xz_F3y_b = 0.0E0;
  Double I_KINETIC_G2x2y_F3y_b = 0.0E0;
  Double I_KINETIC_G2xyz_F3y_b = 0.0E0;
  Double I_KINETIC_G2x2z_F3y_b = 0.0E0;
  Double I_KINETIC_Gx3y_F3y_b = 0.0E0;
  Double I_KINETIC_Gx2yz_F3y_b = 0.0E0;
  Double I_KINETIC_Gxy2z_F3y_b = 0.0E0;
  Double I_KINETIC_Gx3z_F3y_b = 0.0E0;
  Double I_KINETIC_G4y_F3y_b = 0.0E0;
  Double I_KINETIC_G3yz_F3y_b = 0.0E0;
  Double I_KINETIC_G2y2z_F3y_b = 0.0E0;
  Double I_KINETIC_Gy3z_F3y_b = 0.0E0;
  Double I_KINETIC_G4z_F3y_b = 0.0E0;
  Double I_KINETIC_G4x_F2yz_b = 0.0E0;
  Double I_KINETIC_G3xy_F2yz_b = 0.0E0;
  Double I_KINETIC_G3xz_F2yz_b = 0.0E0;
  Double I_KINETIC_G2x2y_F2yz_b = 0.0E0;
  Double I_KINETIC_G2xyz_F2yz_b = 0.0E0;
  Double I_KINETIC_G2x2z_F2yz_b = 0.0E0;
  Double I_KINETIC_Gx3y_F2yz_b = 0.0E0;
  Double I_KINETIC_Gx2yz_F2yz_b = 0.0E0;
  Double I_KINETIC_Gxy2z_F2yz_b = 0.0E0;
  Double I_KINETIC_Gx3z_F2yz_b = 0.0E0;
  Double I_KINETIC_G4y_F2yz_b = 0.0E0;
  Double I_KINETIC_G3yz_F2yz_b = 0.0E0;
  Double I_KINETIC_G2y2z_F2yz_b = 0.0E0;
  Double I_KINETIC_Gy3z_F2yz_b = 0.0E0;
  Double I_KINETIC_G4z_F2yz_b = 0.0E0;
  Double I_KINETIC_G4x_Fy2z_b = 0.0E0;
  Double I_KINETIC_G3xy_Fy2z_b = 0.0E0;
  Double I_KINETIC_G3xz_Fy2z_b = 0.0E0;
  Double I_KINETIC_G2x2y_Fy2z_b = 0.0E0;
  Double I_KINETIC_G2xyz_Fy2z_b = 0.0E0;
  Double I_KINETIC_G2x2z_Fy2z_b = 0.0E0;
  Double I_KINETIC_Gx3y_Fy2z_b = 0.0E0;
  Double I_KINETIC_Gx2yz_Fy2z_b = 0.0E0;
  Double I_KINETIC_Gxy2z_Fy2z_b = 0.0E0;
  Double I_KINETIC_Gx3z_Fy2z_b = 0.0E0;
  Double I_KINETIC_G4y_Fy2z_b = 0.0E0;
  Double I_KINETIC_G3yz_Fy2z_b = 0.0E0;
  Double I_KINETIC_G2y2z_Fy2z_b = 0.0E0;
  Double I_KINETIC_Gy3z_Fy2z_b = 0.0E0;
  Double I_KINETIC_G4z_Fy2z_b = 0.0E0;
  Double I_KINETIC_G4x_F3z_b = 0.0E0;
  Double I_KINETIC_G3xy_F3z_b = 0.0E0;
  Double I_KINETIC_G3xz_F3z_b = 0.0E0;
  Double I_KINETIC_G2x2y_F3z_b = 0.0E0;
  Double I_KINETIC_G2xyz_F3z_b = 0.0E0;
  Double I_KINETIC_G2x2z_F3z_b = 0.0E0;
  Double I_KINETIC_Gx3y_F3z_b = 0.0E0;
  Double I_KINETIC_Gx2yz_F3z_b = 0.0E0;
  Double I_KINETIC_Gxy2z_F3z_b = 0.0E0;
  Double I_KINETIC_Gx3z_F3z_b = 0.0E0;
  Double I_KINETIC_G4y_F3z_b = 0.0E0;
  Double I_KINETIC_G3yz_F3z_b = 0.0E0;
  Double I_KINETIC_G2y2z_F3z_b = 0.0E0;
  Double I_KINETIC_Gy3z_F3z_b = 0.0E0;
  Double I_KINETIC_G4z_F3z_b = 0.0E0;
  Double I_KINETIC_G4x_Px = 0.0E0;
  Double I_KINETIC_G3xy_Px = 0.0E0;
  Double I_KINETIC_G3xz_Px = 0.0E0;
  Double I_KINETIC_G2x2y_Px = 0.0E0;
  Double I_KINETIC_G2xyz_Px = 0.0E0;
  Double I_KINETIC_G2x2z_Px = 0.0E0;
  Double I_KINETIC_Gx3y_Px = 0.0E0;
  Double I_KINETIC_Gx2yz_Px = 0.0E0;
  Double I_KINETIC_Gxy2z_Px = 0.0E0;
  Double I_KINETIC_Gx3z_Px = 0.0E0;
  Double I_KINETIC_G4y_Px = 0.0E0;
  Double I_KINETIC_G3yz_Px = 0.0E0;
  Double I_KINETIC_G2y2z_Px = 0.0E0;
  Double I_KINETIC_Gy3z_Px = 0.0E0;
  Double I_KINETIC_G4z_Px = 0.0E0;
  Double I_KINETIC_G4x_Py = 0.0E0;
  Double I_KINETIC_G3xy_Py = 0.0E0;
  Double I_KINETIC_G3xz_Py = 0.0E0;
  Double I_KINETIC_G2x2y_Py = 0.0E0;
  Double I_KINETIC_G2xyz_Py = 0.0E0;
  Double I_KINETIC_G2x2z_Py = 0.0E0;
  Double I_KINETIC_Gx3y_Py = 0.0E0;
  Double I_KINETIC_Gx2yz_Py = 0.0E0;
  Double I_KINETIC_Gxy2z_Py = 0.0E0;
  Double I_KINETIC_Gx3z_Py = 0.0E0;
  Double I_KINETIC_G4y_Py = 0.0E0;
  Double I_KINETIC_G3yz_Py = 0.0E0;
  Double I_KINETIC_G2y2z_Py = 0.0E0;
  Double I_KINETIC_Gy3z_Py = 0.0E0;
  Double I_KINETIC_G4z_Py = 0.0E0;
  Double I_KINETIC_G4x_Pz = 0.0E0;
  Double I_KINETIC_G3xy_Pz = 0.0E0;
  Double I_KINETIC_G3xz_Pz = 0.0E0;
  Double I_KINETIC_G2x2y_Pz = 0.0E0;
  Double I_KINETIC_G2xyz_Pz = 0.0E0;
  Double I_KINETIC_G2x2z_Pz = 0.0E0;
  Double I_KINETIC_Gx3y_Pz = 0.0E0;
  Double I_KINETIC_Gx2yz_Pz = 0.0E0;
  Double I_KINETIC_Gxy2z_Pz = 0.0E0;
  Double I_KINETIC_Gx3z_Pz = 0.0E0;
  Double I_KINETIC_G4y_Pz = 0.0E0;
  Double I_KINETIC_G3yz_Pz = 0.0E0;
  Double I_KINETIC_G2y2z_Pz = 0.0E0;
  Double I_KINETIC_Gy3z_Pz = 0.0E0;
  Double I_KINETIC_G4z_Pz = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 20 integrals are omitted 
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
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_Px_vrr = PAX*I_TWOBODYOVERLAP_F3x_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Px_vrr = PAY*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Px_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Px_vrr = PAX*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Px_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
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
    Double I_TWOBODYOVERLAP_G2xyz_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Py_vrr = PAX*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Py_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
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
    Double I_TWOBODYOVERLAP_G2xyz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Pz_vrr = PAX*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Pz_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Pz_vrr;
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
     * totally 30 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_D2x_vrr = PBX*I_TWOBODYOVERLAP_G4x_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2x_vrr = PBX*I_TWOBODYOVERLAP_G3xy_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2x_vrr = PBX*I_TWOBODYOVERLAP_G3xz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_D2x_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
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
    Double I_TWOBODYOVERLAP_G2xyz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4y_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G4y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G3yz_Px_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_Px_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_Dxy_vrr = PBY*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_D2y_vrr = PBY*I_TWOBODYOVERLAP_G4x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2y_vrr = PBY*I_TWOBODYOVERLAP_G3xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2y_vrr = PBY*I_TWOBODYOVERLAP_G3xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2y_vrr = PBY*I_TWOBODYOVERLAP_G4y_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2y_vrr = PBY*I_TWOBODYOVERLAP_G3yz_Py_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2y_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_Py_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2y_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2y_vrr = PBY*I_TWOBODYOVERLAP_G4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G4x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_Pz_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_Pz_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_D2z_vrr = PBZ*I_TWOBODYOVERLAP_G4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_F3x_vrr = PBX*I_TWOBODYOVERLAP_G4x_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3x_vrr = PBX*I_TWOBODYOVERLAP_G3xy_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3x_vrr = PBX*I_TWOBODYOVERLAP_G3xz_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3x_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F3x_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3x_vrr = PBX*I_TWOBODYOVERLAP_G4y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_G3yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3x_vrr = PBX*I_TWOBODYOVERLAP_G4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G3xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G3xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G4y_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G3yz_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G4x_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G3xy_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G3xz_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G3yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_Dxy_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G4x_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G3xy_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G3xz_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G3yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4x_F3y_vrr = PBY*I_TWOBODYOVERLAP_G4x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3y_vrr = PBY*I_TWOBODYOVERLAP_G3xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3y_vrr = PBY*I_TWOBODYOVERLAP_G3xz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3y_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F3y_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3y_vrr = PBY*I_TWOBODYOVERLAP_G4y_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_G3yz_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3y_vrr = PBY*I_TWOBODYOVERLAP_G4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G3xy_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G3xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G4y_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G3yz_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4x_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G4x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_G4z_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G4z_Pz_vrr;

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
     * shell quartet name: SQ_KINETIC_F_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 20 integrals are omitted 
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
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_KINETIC_G4x_Px_vrr = PBX*I_KINETIC_G4x_S_vrr+4*oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_KINETIC_G3xy_Px_vrr = PBX*I_KINETIC_G3xy_S_vrr+3*oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_KINETIC_G3xz_Px_vrr = PBX*I_KINETIC_G3xz_S_vrr+3*oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_KINETIC_G2x2y_Px_vrr = PBX*I_KINETIC_G2x2y_S_vrr+2*oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_KINETIC_G2xyz_Px_vrr = PBX*I_KINETIC_G2xyz_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_KINETIC_G2x2z_Px_vrr = PBX*I_KINETIC_G2x2z_S_vrr+2*oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_KINETIC_Gx3y_Px_vrr = PBX*I_KINETIC_Gx3y_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_KINETIC_Gx2yz_Px_vrr = PBX*I_KINETIC_Gx2yz_S_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_KINETIC_Gxy2z_Px_vrr = PBX*I_KINETIC_Gxy2z_S_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_KINETIC_Gx3z_Px_vrr = PBX*I_KINETIC_Gx3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_KINETIC_G4y_Px_vrr = PBX*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_KINETIC_G3yz_Px_vrr = PBX*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_KINETIC_G2y2z_Px_vrr = PBX*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_KINETIC_Gy3z_Px_vrr = PBX*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_KINETIC_G4z_Px_vrr = PBX*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_KINETIC_G4x_Py_vrr = PBY*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_KINETIC_G3xy_Py_vrr = PBY*I_KINETIC_G3xy_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_KINETIC_G3xz_Py_vrr = PBY*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_KINETIC_G2x2y_Py_vrr = PBY*I_KINETIC_G2x2y_S_vrr+2*oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_KINETIC_G2xyz_Py_vrr = PBY*I_KINETIC_G2xyz_S_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_KINETIC_G2x2z_Py_vrr = PBY*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_KINETIC_Gx3y_Py_vrr = PBY*I_KINETIC_Gx3y_S_vrr+3*oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_KINETIC_Gx2yz_Py_vrr = PBY*I_KINETIC_Gx2yz_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_KINETIC_Gxy2z_Py_vrr = PBY*I_KINETIC_Gxy2z_S_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_KINETIC_Gx3z_Py_vrr = PBY*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_KINETIC_G4y_Py_vrr = PBY*I_KINETIC_G4y_S_vrr+4*oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_KINETIC_G3yz_Py_vrr = PBY*I_KINETIC_G3yz_S_vrr+3*oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_KINETIC_G2y2z_Py_vrr = PBY*I_KINETIC_G2y2z_S_vrr+2*oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_KINETIC_Gy3z_Py_vrr = PBY*I_KINETIC_Gy3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_KINETIC_G4z_Py_vrr = PBY*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_KINETIC_G4x_Pz_vrr = PBZ*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_KINETIC_G3xy_Pz_vrr = PBZ*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_KINETIC_G3xz_Pz_vrr = PBZ*I_KINETIC_G3xz_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_KINETIC_G2x2y_Pz_vrr = PBZ*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_KINETIC_G2xyz_Pz_vrr = PBZ*I_KINETIC_G2xyz_S_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Pz_vrr;
    Double I_KINETIC_G2x2z_Pz_vrr = PBZ*I_KINETIC_G2x2z_S_vrr+2*oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_KINETIC_Gx3y_Pz_vrr = PBZ*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_KINETIC_Gx2yz_Pz_vrr = PBZ*I_KINETIC_Gx2yz_S_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr;
    Double I_KINETIC_Gxy2z_Pz_vrr = PBZ*I_KINETIC_Gxy2z_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr;
    Double I_KINETIC_Gx3z_Pz_vrr = PBZ*I_KINETIC_Gx3z_S_vrr+3*oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_KINETIC_G4y_Pz_vrr = PBZ*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_KINETIC_G3yz_Pz_vrr = PBZ*I_KINETIC_G3yz_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_KINETIC_G2y2z_Pz_vrr = PBZ*I_KINETIC_G2y2z_S_vrr+2*oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_KINETIC_Gy3z_Pz_vrr = PBZ*I_KINETIC_Gy3z_S_vrr+3*oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_KINETIC_G4z_Pz_vrr = PBZ*I_KINETIC_G4z_S_vrr+4*oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 30 integrals are omitted 
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
    Double I_KINETIC_G2xyz_D2x_vrr = PBX*I_KINETIC_G2xyz_Px_vrr+2*oned2z*I_KINETIC_Fxyz_Px_vrr+oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_D2x_vrr-adz*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_D2x_vrr = PBX*I_KINETIC_G2x2z_Px_vrr+2*oned2z*I_KINETIC_Fx2z_Px_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2x_vrr-adz*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_KINETIC_Gx3y_D2x_vrr = PBX*I_KINETIC_Gx3y_Px_vrr+oned2z*I_KINETIC_F3y_Px_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2x_vrr-adz*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_D2x_vrr = PBX*I_KINETIC_Gx2yz_Px_vrr+oned2z*I_KINETIC_F2yz_Px_vrr+oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr-adz*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_D2x_vrr = PBX*I_KINETIC_Gxy2z_Px_vrr+oned2z*I_KINETIC_Fy2z_Px_vrr+oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
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
    Double I_KINETIC_G2xyz_Dxy_vrr = PBY*I_KINETIC_G2xyz_Px_vrr+oned2z*I_KINETIC_F2xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Dxy_vrr;
    Double I_KINETIC_G2x2z_Dxy_vrr = PBY*I_KINETIC_G2x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Dxy_vrr;
    Double I_KINETIC_Gx3y_Dxy_vrr = PBY*I_KINETIC_Gx3y_Px_vrr+3*oned2z*I_KINETIC_Fx2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Dxy_vrr;
    Double I_KINETIC_Gx2yz_Dxy_vrr = PBY*I_KINETIC_Gx2yz_Px_vrr+2*oned2z*I_KINETIC_Fxyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Dxy_vrr;
    Double I_KINETIC_Gxy2z_Dxy_vrr = PBY*I_KINETIC_Gxy2z_Px_vrr+oned2z*I_KINETIC_Fx2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Dxy_vrr;
    Double I_KINETIC_Gx3z_Dxy_vrr = PBY*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Dxy_vrr;
    Double I_KINETIC_G4y_Dxy_vrr = PBY*I_KINETIC_G4y_Px_vrr+4*oned2z*I_KINETIC_F3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Dxy_vrr;
    Double I_KINETIC_G3yz_Dxy_vrr = PBY*I_KINETIC_G3yz_Px_vrr+3*oned2z*I_KINETIC_F2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Dxy_vrr;
    Double I_KINETIC_G2y2z_Dxy_vrr = PBY*I_KINETIC_G2y2z_Px_vrr+2*oned2z*I_KINETIC_Fy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Dxy_vrr;
    Double I_KINETIC_Gy3z_Dxy_vrr = PBY*I_KINETIC_Gy3z_Px_vrr+oned2z*I_KINETIC_F3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Dxy_vrr;
    Double I_KINETIC_G4z_Dxy_vrr = PBY*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Dxy_vrr;
    Double I_KINETIC_G4x_D2y_vrr = PBY*I_KINETIC_G4x_Py_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2y_vrr-adz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_G3xy_D2y_vrr = PBY*I_KINETIC_G3xy_Py_vrr+oned2z*I_KINETIC_F3x_Py_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2y_vrr-adz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_D2y_vrr = PBY*I_KINETIC_G3xz_Py_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2y_vrr-adz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_D2y_vrr = PBY*I_KINETIC_G2x2y_Py_vrr+2*oned2z*I_KINETIC_F2xy_Py_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2y_vrr-adz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_G2xyz_D2y_vrr = PBY*I_KINETIC_G2xyz_Py_vrr+oned2z*I_KINETIC_F2xz_Py_vrr+oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_D2y_vrr-adz*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_D2y_vrr = PBY*I_KINETIC_G2x2z_Py_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2y_vrr-adz*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_KINETIC_Gx3y_D2y_vrr = PBY*I_KINETIC_Gx3y_Py_vrr+3*oned2z*I_KINETIC_Fx2y_Py_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2y_vrr-adz*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_D2y_vrr = PBY*I_KINETIC_Gx2yz_Py_vrr+2*oned2z*I_KINETIC_Fxyz_Py_vrr+oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr-adz*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_D2y_vrr = PBY*I_KINETIC_Gxy2z_Py_vrr+oned2z*I_KINETIC_Fx2z_Py_vrr+oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_KINETIC_Gx3z_D2y_vrr = PBY*I_KINETIC_Gx3z_Py_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2y_vrr-adz*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_D2y_vrr = PBY*I_KINETIC_G4y_Py_vrr+4*oned2z*I_KINETIC_F3y_Py_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2y_vrr-adz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_G3yz_D2y_vrr = PBY*I_KINETIC_G3yz_Py_vrr+3*oned2z*I_KINETIC_F2yz_Py_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2y_vrr-adz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_D2y_vrr = PBY*I_KINETIC_G2y2z_Py_vrr+2*oned2z*I_KINETIC_Fy2z_Py_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2y_vrr-adz*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_KINETIC_Gy3z_D2y_vrr = PBY*I_KINETIC_Gy3z_Py_vrr+oned2z*I_KINETIC_F3z_Py_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2y_vrr-adz*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_D2y_vrr = PBY*I_KINETIC_G4z_Py_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2y_vrr-adz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_G4x_D2z_vrr = PBZ*I_KINETIC_G4x_Pz_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_D2z_vrr-adz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_G3xy_D2z_vrr = PBZ*I_KINETIC_G3xy_Pz_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_D2z_vrr-adz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_D2z_vrr = PBZ*I_KINETIC_G3xz_Pz_vrr+oned2z*I_KINETIC_F3x_Pz_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_D2z_vrr-adz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_D2z_vrr = PBZ*I_KINETIC_G2x2y_Pz_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_D2z_vrr-adz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_G2xyz_D2z_vrr = PBZ*I_KINETIC_G2xyz_Pz_vrr+oned2z*I_KINETIC_F2xy_Pz_vrr+oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_D2z_vrr-adz*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_D2z_vrr = PBZ*I_KINETIC_G2x2z_Pz_vrr+2*oned2z*I_KINETIC_F2xz_Pz_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_D2z_vrr-adz*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_KINETIC_Gx3y_D2z_vrr = PBZ*I_KINETIC_Gx3y_Pz_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_D2z_vrr-adz*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_D2z_vrr = PBZ*I_KINETIC_Gx2yz_Pz_vrr+oned2z*I_KINETIC_Fx2y_Pz_vrr+oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr-adz*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_D2z_vrr = PBZ*I_KINETIC_Gxy2z_Pz_vrr+2*oned2z*I_KINETIC_Fxyz_Pz_vrr+oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_KINETIC_Gx3z_D2z_vrr = PBZ*I_KINETIC_Gx3z_Pz_vrr+3*oned2z*I_KINETIC_Fx2z_Pz_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_D2z_vrr-adz*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_D2z_vrr = PBZ*I_KINETIC_G4y_Pz_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_D2z_vrr-adz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_G3yz_D2z_vrr = PBZ*I_KINETIC_G3yz_Pz_vrr+oned2z*I_KINETIC_F3y_Pz_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_D2z_vrr-adz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_D2z_vrr = PBZ*I_KINETIC_G2y2z_Pz_vrr+2*oned2z*I_KINETIC_F2yz_Pz_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_D2z_vrr-adz*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_KINETIC_Gy3z_D2z_vrr = PBZ*I_KINETIC_Gy3z_Pz_vrr+3*oned2z*I_KINETIC_Fy2z_Pz_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_D2z_vrr-adz*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_D2z_vrr = PBZ*I_KINETIC_G4z_Pz_vrr+4*oned2z*I_KINETIC_F3z_Pz_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_D2z_vrr-adz*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_D
     * RHS shell quartet name: SQ_KINETIC_F_D
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_KINETIC_G4x_F3x_vrr = PBX*I_KINETIC_G4x_D2x_vrr+4*oned2z*I_KINETIC_F3x_D2x_vrr+2*oned2z*I_KINETIC_G4x_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_KINETIC_G3xy_F3x_vrr = PBX*I_KINETIC_G3xy_D2x_vrr+3*oned2z*I_KINETIC_F2xy_D2x_vrr+2*oned2z*I_KINETIC_G3xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_KINETIC_G3xz_F3x_vrr = PBX*I_KINETIC_G3xz_D2x_vrr+3*oned2z*I_KINETIC_F2xz_D2x_vrr+2*oned2z*I_KINETIC_G3xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_KINETIC_G2x2y_F3x_vrr = PBX*I_KINETIC_G2x2y_D2x_vrr+2*oned2z*I_KINETIC_Fx2y_D2x_vrr+2*oned2z*I_KINETIC_G2x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_KINETIC_G2xyz_F3x_vrr = PBX*I_KINETIC_G2xyz_D2x_vrr+2*oned2z*I_KINETIC_Fxyz_D2x_vrr+2*oned2z*I_KINETIC_G2xyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_KINETIC_G2x2z_F3x_vrr = PBX*I_KINETIC_G2x2z_D2x_vrr+2*oned2z*I_KINETIC_Fx2z_D2x_vrr+2*oned2z*I_KINETIC_G2x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_KINETIC_Gx3y_F3x_vrr = PBX*I_KINETIC_Gx3y_D2x_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+2*oned2z*I_KINETIC_Gx3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_KINETIC_Gx2yz_F3x_vrr = PBX*I_KINETIC_Gx2yz_D2x_vrr+oned2z*I_KINETIC_F2yz_D2x_vrr+2*oned2z*I_KINETIC_Gx2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_KINETIC_Gxy2z_F3x_vrr = PBX*I_KINETIC_Gxy2z_D2x_vrr+oned2z*I_KINETIC_Fy2z_D2x_vrr+2*oned2z*I_KINETIC_Gxy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_KINETIC_Gx3z_F3x_vrr = PBX*I_KINETIC_Gx3z_D2x_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+2*oned2z*I_KINETIC_Gx3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_KINETIC_G4y_F3x_vrr = PBX*I_KINETIC_G4y_D2x_vrr+2*oned2z*I_KINETIC_G4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_KINETIC_G3yz_F3x_vrr = PBX*I_KINETIC_G3yz_D2x_vrr+2*oned2z*I_KINETIC_G3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_KINETIC_G2y2z_F3x_vrr = PBX*I_KINETIC_G2y2z_D2x_vrr+2*oned2z*I_KINETIC_G2y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_KINETIC_Gy3z_F3x_vrr = PBX*I_KINETIC_Gy3z_D2x_vrr+2*oned2z*I_KINETIC_Gy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_KINETIC_G4z_F3x_vrr = PBX*I_KINETIC_G4z_D2x_vrr+2*oned2z*I_KINETIC_G4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_KINETIC_G4x_F2xy_vrr = PBY*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_KINETIC_G3xy_F2xy_vrr = PBY*I_KINETIC_G3xy_D2x_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2xy_vrr;
    Double I_KINETIC_G3xz_F2xy_vrr = PBY*I_KINETIC_G3xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2xy_vrr;
    Double I_KINETIC_G2x2y_F2xy_vrr = PBY*I_KINETIC_G2x2y_D2x_vrr+2*oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr;
    Double I_KINETIC_G2xyz_F2xy_vrr = PBY*I_KINETIC_G2xyz_D2x_vrr+oned2z*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F2xy_vrr;
    Double I_KINETIC_G2x2z_F2xy_vrr = PBY*I_KINETIC_G2x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2xy_vrr;
    Double I_KINETIC_Gx3y_F2xy_vrr = PBY*I_KINETIC_Gx3y_D2x_vrr+3*oned2z*I_KINETIC_Fx2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr;
    Double I_KINETIC_Gx2yz_F2xy_vrr = PBY*I_KINETIC_Gx2yz_D2x_vrr+2*oned2z*I_KINETIC_Fxyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F2xy_vrr;
    Double I_KINETIC_Gxy2z_F2xy_vrr = PBY*I_KINETIC_Gxy2z_D2x_vrr+oned2z*I_KINETIC_Fx2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F2xy_vrr;
    Double I_KINETIC_Gx3z_F2xy_vrr = PBY*I_KINETIC_Gx3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr;
    Double I_KINETIC_G4y_F2xy_vrr = PBY*I_KINETIC_G4y_D2x_vrr+4*oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_KINETIC_G3yz_F2xy_vrr = PBY*I_KINETIC_G3yz_D2x_vrr+3*oned2z*I_KINETIC_F2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2xy_vrr;
    Double I_KINETIC_G2y2z_F2xy_vrr = PBY*I_KINETIC_G2y2z_D2x_vrr+2*oned2z*I_KINETIC_Fy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2xy_vrr;
    Double I_KINETIC_Gy3z_F2xy_vrr = PBY*I_KINETIC_Gy3z_D2x_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2xy_vrr;
    Double I_KINETIC_G4z_F2xy_vrr = PBY*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_KINETIC_G4x_F2xz_vrr = PBZ*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2xz_vrr;
    Double I_KINETIC_G3xy_F2xz_vrr = PBZ*I_KINETIC_G3xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2xz_vrr;
    Double I_KINETIC_G3xz_F2xz_vrr = PBZ*I_KINETIC_G3xz_D2x_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2xz_vrr;
    Double I_KINETIC_G2x2y_F2xz_vrr = PBZ*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2xz_vrr;
    Double I_KINETIC_G2xyz_F2xz_vrr = PBZ*I_KINETIC_G2xyz_D2x_vrr+oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F2xz_vrr;
    Double I_KINETIC_G2x2z_F2xz_vrr = PBZ*I_KINETIC_G2x2z_D2x_vrr+2*oned2z*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2xz_vrr;
    Double I_KINETIC_Gx3y_F2xz_vrr = PBZ*I_KINETIC_Gx3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2xz_vrr;
    Double I_KINETIC_Gx2yz_F2xz_vrr = PBZ*I_KINETIC_Gx2yz_D2x_vrr+oned2z*I_KINETIC_Fx2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F2xz_vrr;
    Double I_KINETIC_Gxy2z_F2xz_vrr = PBZ*I_KINETIC_Gxy2z_D2x_vrr+2*oned2z*I_KINETIC_Fxyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F2xz_vrr;
    Double I_KINETIC_Gx3z_F2xz_vrr = PBZ*I_KINETIC_Gx3z_D2x_vrr+3*oned2z*I_KINETIC_Fx2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2xz_vrr;
    Double I_KINETIC_G4y_F2xz_vrr = PBZ*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2xz_vrr;
    Double I_KINETIC_G3yz_F2xz_vrr = PBZ*I_KINETIC_G3yz_D2x_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2xz_vrr;
    Double I_KINETIC_G2y2z_F2xz_vrr = PBZ*I_KINETIC_G2y2z_D2x_vrr+2*oned2z*I_KINETIC_F2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2xz_vrr;
    Double I_KINETIC_Gy3z_F2xz_vrr = PBZ*I_KINETIC_Gy3z_D2x_vrr+3*oned2z*I_KINETIC_Fy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2xz_vrr;
    Double I_KINETIC_G4z_F2xz_vrr = PBZ*I_KINETIC_G4z_D2x_vrr+4*oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2xz_vrr;
    Double I_KINETIC_G4x_Fx2y_vrr = PBX*I_KINETIC_G4x_D2y_vrr+4*oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_KINETIC_G3xy_Fx2y_vrr = PBX*I_KINETIC_G3xy_D2y_vrr+3*oned2z*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr;
    Double I_KINETIC_G3xz_Fx2y_vrr = PBX*I_KINETIC_G3xz_D2y_vrr+3*oned2z*I_KINETIC_F2xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fx2y_vrr;
    Double I_KINETIC_G2x2y_Fx2y_vrr = PBX*I_KINETIC_G2x2y_D2y_vrr+2*oned2z*I_KINETIC_Fx2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr;
    Double I_KINETIC_G2xyz_Fx2y_vrr = PBX*I_KINETIC_G2xyz_D2y_vrr+2*oned2z*I_KINETIC_Fxyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fx2y_vrr;
    Double I_KINETIC_G2x2z_Fx2y_vrr = PBX*I_KINETIC_G2x2z_D2y_vrr+2*oned2z*I_KINETIC_Fx2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr;
    Double I_KINETIC_Gx3y_Fx2y_vrr = PBX*I_KINETIC_Gx3y_D2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr;
    Double I_KINETIC_Gx2yz_Fx2y_vrr = PBX*I_KINETIC_Gx2yz_D2y_vrr+oned2z*I_KINETIC_F2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fx2y_vrr;
    Double I_KINETIC_Gxy2z_Fx2y_vrr = PBX*I_KINETIC_Gxy2z_D2y_vrr+oned2z*I_KINETIC_Fy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fx2y_vrr;
    Double I_KINETIC_Gx3z_Fx2y_vrr = PBX*I_KINETIC_Gx3z_D2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr;
    Double I_KINETIC_G4y_Fx2y_vrr = PBX*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_KINETIC_G3yz_Fx2y_vrr = PBX*I_KINETIC_G3yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fx2y_vrr;
    Double I_KINETIC_G2y2z_Fx2y_vrr = PBX*I_KINETIC_G2y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr;
    Double I_KINETIC_Gy3z_Fx2y_vrr = PBX*I_KINETIC_Gy3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr;
    Double I_KINETIC_G4z_Fx2y_vrr = PBX*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_KINETIC_G4x_Fxyz_vrr = PBZ*I_KINETIC_G4x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fxyz_vrr;
    Double I_KINETIC_G3xy_Fxyz_vrr = PBZ*I_KINETIC_G3xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr;
    Double I_KINETIC_G3xz_Fxyz_vrr = PBZ*I_KINETIC_G3xz_Dxy_vrr+oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fxyz_vrr;
    Double I_KINETIC_G2x2y_Fxyz_vrr = PBZ*I_KINETIC_G2x2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr;
    Double I_KINETIC_G2xyz_Fxyz_vrr = PBZ*I_KINETIC_G2xyz_Dxy_vrr+oned2z*I_KINETIC_F2xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fxyz_vrr;
    Double I_KINETIC_G2x2z_Fxyz_vrr = PBZ*I_KINETIC_G2x2z_Dxy_vrr+2*oned2z*I_KINETIC_F2xz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr;
    Double I_KINETIC_Gx3y_Fxyz_vrr = PBZ*I_KINETIC_Gx3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr;
    Double I_KINETIC_Gx2yz_Fxyz_vrr = PBZ*I_KINETIC_Gx2yz_Dxy_vrr+oned2z*I_KINETIC_Fx2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fxyz_vrr;
    Double I_KINETIC_Gxy2z_Fxyz_vrr = PBZ*I_KINETIC_Gxy2z_Dxy_vrr+2*oned2z*I_KINETIC_Fxyz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fxyz_vrr;
    Double I_KINETIC_Gx3z_Fxyz_vrr = PBZ*I_KINETIC_Gx3z_Dxy_vrr+3*oned2z*I_KINETIC_Fx2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr;
    Double I_KINETIC_G4y_Fxyz_vrr = PBZ*I_KINETIC_G4y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fxyz_vrr;
    Double I_KINETIC_G3yz_Fxyz_vrr = PBZ*I_KINETIC_G3yz_Dxy_vrr+oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fxyz_vrr;
    Double I_KINETIC_G2y2z_Fxyz_vrr = PBZ*I_KINETIC_G2y2z_Dxy_vrr+2*oned2z*I_KINETIC_F2yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr;
    Double I_KINETIC_Gy3z_Fxyz_vrr = PBZ*I_KINETIC_Gy3z_Dxy_vrr+3*oned2z*I_KINETIC_Fy2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr;
    Double I_KINETIC_G4z_Fxyz_vrr = PBZ*I_KINETIC_G4z_Dxy_vrr+4*oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fxyz_vrr;
    Double I_KINETIC_G4x_Fx2z_vrr = PBX*I_KINETIC_G4x_D2z_vrr+4*oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_KINETIC_G3xy_Fx2z_vrr = PBX*I_KINETIC_G3xy_D2z_vrr+3*oned2z*I_KINETIC_F2xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr;
    Double I_KINETIC_G3xz_Fx2z_vrr = PBX*I_KINETIC_G3xz_D2z_vrr+3*oned2z*I_KINETIC_F2xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fx2z_vrr;
    Double I_KINETIC_G2x2y_Fx2z_vrr = PBX*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_Fx2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr;
    Double I_KINETIC_G2xyz_Fx2z_vrr = PBX*I_KINETIC_G2xyz_D2z_vrr+2*oned2z*I_KINETIC_Fxyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fx2z_vrr;
    Double I_KINETIC_G2x2z_Fx2z_vrr = PBX*I_KINETIC_G2x2z_D2z_vrr+2*oned2z*I_KINETIC_Fx2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr;
    Double I_KINETIC_Gx3y_Fx2z_vrr = PBX*I_KINETIC_Gx3y_D2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr;
    Double I_KINETIC_Gx2yz_Fx2z_vrr = PBX*I_KINETIC_Gx2yz_D2z_vrr+oned2z*I_KINETIC_F2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fx2z_vrr;
    Double I_KINETIC_Gxy2z_Fx2z_vrr = PBX*I_KINETIC_Gxy2z_D2z_vrr+oned2z*I_KINETIC_Fy2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fx2z_vrr;
    Double I_KINETIC_Gx3z_Fx2z_vrr = PBX*I_KINETIC_Gx3z_D2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr;
    Double I_KINETIC_G4y_Fx2z_vrr = PBX*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_KINETIC_G3yz_Fx2z_vrr = PBX*I_KINETIC_G3yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fx2z_vrr;
    Double I_KINETIC_G2y2z_Fx2z_vrr = PBX*I_KINETIC_G2y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr;
    Double I_KINETIC_Gy3z_Fx2z_vrr = PBX*I_KINETIC_Gy3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr;
    Double I_KINETIC_G4z_Fx2z_vrr = PBX*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_KINETIC_G4x_F3y_vrr = PBY*I_KINETIC_G4x_D2y_vrr+2*oned2z*I_KINETIC_G4x_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_KINETIC_G3xy_F3y_vrr = PBY*I_KINETIC_G3xy_D2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+2*oned2z*I_KINETIC_G3xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_KINETIC_G3xz_F3y_vrr = PBY*I_KINETIC_G3xz_D2y_vrr+2*oned2z*I_KINETIC_G3xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_KINETIC_G2x2y_F3y_vrr = PBY*I_KINETIC_G2x2y_D2y_vrr+2*oned2z*I_KINETIC_F2xy_D2y_vrr+2*oned2z*I_KINETIC_G2x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_KINETIC_G2xyz_F3y_vrr = PBY*I_KINETIC_G2xyz_D2y_vrr+oned2z*I_KINETIC_F2xz_D2y_vrr+2*oned2z*I_KINETIC_G2xyz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_KINETIC_G2x2z_F3y_vrr = PBY*I_KINETIC_G2x2z_D2y_vrr+2*oned2z*I_KINETIC_G2x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_KINETIC_Gx3y_F3y_vrr = PBY*I_KINETIC_Gx3y_D2y_vrr+3*oned2z*I_KINETIC_Fx2y_D2y_vrr+2*oned2z*I_KINETIC_Gx3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_KINETIC_Gx2yz_F3y_vrr = PBY*I_KINETIC_Gx2yz_D2y_vrr+2*oned2z*I_KINETIC_Fxyz_D2y_vrr+2*oned2z*I_KINETIC_Gx2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_KINETIC_Gxy2z_F3y_vrr = PBY*I_KINETIC_Gxy2z_D2y_vrr+oned2z*I_KINETIC_Fx2z_D2y_vrr+2*oned2z*I_KINETIC_Gxy2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_KINETIC_Gx3z_F3y_vrr = PBY*I_KINETIC_Gx3z_D2y_vrr+2*oned2z*I_KINETIC_Gx3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_KINETIC_G4y_F3y_vrr = PBY*I_KINETIC_G4y_D2y_vrr+4*oned2z*I_KINETIC_F3y_D2y_vrr+2*oned2z*I_KINETIC_G4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_KINETIC_G3yz_F3y_vrr = PBY*I_KINETIC_G3yz_D2y_vrr+3*oned2z*I_KINETIC_F2yz_D2y_vrr+2*oned2z*I_KINETIC_G3yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_KINETIC_G2y2z_F3y_vrr = PBY*I_KINETIC_G2y2z_D2y_vrr+2*oned2z*I_KINETIC_Fy2z_D2y_vrr+2*oned2z*I_KINETIC_G2y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_KINETIC_Gy3z_F3y_vrr = PBY*I_KINETIC_Gy3z_D2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+2*oned2z*I_KINETIC_Gy3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_KINETIC_G4z_F3y_vrr = PBY*I_KINETIC_G4z_D2y_vrr+2*oned2z*I_KINETIC_G4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_KINETIC_G4x_F2yz_vrr = PBZ*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2yz_vrr;
    Double I_KINETIC_G3xy_F2yz_vrr = PBZ*I_KINETIC_G3xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2yz_vrr;
    Double I_KINETIC_G3xz_F2yz_vrr = PBZ*I_KINETIC_G3xz_D2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2yz_vrr;
    Double I_KINETIC_G2x2y_F2yz_vrr = PBZ*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2yz_vrr;
    Double I_KINETIC_G2xyz_F2yz_vrr = PBZ*I_KINETIC_G2xyz_D2y_vrr+oned2z*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F2yz_vrr;
    Double I_KINETIC_G2x2z_F2yz_vrr = PBZ*I_KINETIC_G2x2z_D2y_vrr+2*oned2z*I_KINETIC_F2xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2yz_vrr;
    Double I_KINETIC_Gx3y_F2yz_vrr = PBZ*I_KINETIC_Gx3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2yz_vrr;
    Double I_KINETIC_Gx2yz_F2yz_vrr = PBZ*I_KINETIC_Gx2yz_D2y_vrr+oned2z*I_KINETIC_Fx2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F2yz_vrr;
    Double I_KINETIC_Gxy2z_F2yz_vrr = PBZ*I_KINETIC_Gxy2z_D2y_vrr+2*oned2z*I_KINETIC_Fxyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F2yz_vrr;
    Double I_KINETIC_Gx3z_F2yz_vrr = PBZ*I_KINETIC_Gx3z_D2y_vrr+3*oned2z*I_KINETIC_Fx2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2yz_vrr;
    Double I_KINETIC_G4y_F2yz_vrr = PBZ*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2yz_vrr;
    Double I_KINETIC_G3yz_F2yz_vrr = PBZ*I_KINETIC_G3yz_D2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2yz_vrr;
    Double I_KINETIC_G2y2z_F2yz_vrr = PBZ*I_KINETIC_G2y2z_D2y_vrr+2*oned2z*I_KINETIC_F2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2yz_vrr;
    Double I_KINETIC_Gy3z_F2yz_vrr = PBZ*I_KINETIC_Gy3z_D2y_vrr+3*oned2z*I_KINETIC_Fy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2yz_vrr;
    Double I_KINETIC_G4z_F2yz_vrr = PBZ*I_KINETIC_G4z_D2y_vrr+4*oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2yz_vrr;
    Double I_KINETIC_G4x_Fy2z_vrr = PBY*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fy2z_vrr;
    Double I_KINETIC_G3xy_Fy2z_vrr = PBY*I_KINETIC_G3xy_D2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr;
    Double I_KINETIC_G3xz_Fy2z_vrr = PBY*I_KINETIC_G3xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fy2z_vrr;
    Double I_KINETIC_G2x2y_Fy2z_vrr = PBY*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_F2xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr;
    Double I_KINETIC_G2xyz_Fy2z_vrr = PBY*I_KINETIC_G2xyz_D2z_vrr+oned2z*I_KINETIC_F2xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fy2z_vrr;
    Double I_KINETIC_G2x2z_Fy2z_vrr = PBY*I_KINETIC_G2x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr;
    Double I_KINETIC_Gx3y_Fy2z_vrr = PBY*I_KINETIC_Gx3y_D2z_vrr+3*oned2z*I_KINETIC_Fx2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr;
    Double I_KINETIC_Gx2yz_Fy2z_vrr = PBY*I_KINETIC_Gx2yz_D2z_vrr+2*oned2z*I_KINETIC_Fxyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fy2z_vrr;
    Double I_KINETIC_Gxy2z_Fy2z_vrr = PBY*I_KINETIC_Gxy2z_D2z_vrr+oned2z*I_KINETIC_Fx2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fy2z_vrr;
    Double I_KINETIC_Gx3z_Fy2z_vrr = PBY*I_KINETIC_Gx3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr;
    Double I_KINETIC_G4y_Fy2z_vrr = PBY*I_KINETIC_G4y_D2z_vrr+4*oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fy2z_vrr;
    Double I_KINETIC_G3yz_Fy2z_vrr = PBY*I_KINETIC_G3yz_D2z_vrr+3*oned2z*I_KINETIC_F2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fy2z_vrr;
    Double I_KINETIC_G2y2z_Fy2z_vrr = PBY*I_KINETIC_G2y2z_D2z_vrr+2*oned2z*I_KINETIC_Fy2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr;
    Double I_KINETIC_Gy3z_Fy2z_vrr = PBY*I_KINETIC_Gy3z_D2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr;
    Double I_KINETIC_G4z_Fy2z_vrr = PBY*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fy2z_vrr;
    Double I_KINETIC_G4x_F3z_vrr = PBZ*I_KINETIC_G4x_D2z_vrr+2*oned2z*I_KINETIC_G4x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_KINETIC_G3xy_F3z_vrr = PBZ*I_KINETIC_G3xy_D2z_vrr+2*oned2z*I_KINETIC_G3xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_KINETIC_G3xz_F3z_vrr = PBZ*I_KINETIC_G3xz_D2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+2*oned2z*I_KINETIC_G3xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_KINETIC_G2x2y_F3z_vrr = PBZ*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_G2x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_KINETIC_G2xyz_F3z_vrr = PBZ*I_KINETIC_G2xyz_D2z_vrr+oned2z*I_KINETIC_F2xy_D2z_vrr+2*oned2z*I_KINETIC_G2xyz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G2xyz_Pz_vrr;
    Double I_KINETIC_G2x2z_F3z_vrr = PBZ*I_KINETIC_G2x2z_D2z_vrr+2*oned2z*I_KINETIC_F2xz_D2z_vrr+2*oned2z*I_KINETIC_G2x2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_KINETIC_Gx3y_F3z_vrr = PBZ*I_KINETIC_Gx3y_D2z_vrr+2*oned2z*I_KINETIC_Gx3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_KINETIC_Gx2yz_F3z_vrr = PBZ*I_KINETIC_Gx2yz_D2z_vrr+oned2z*I_KINETIC_Fx2y_D2z_vrr+2*oned2z*I_KINETIC_Gx2yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr;
    Double I_KINETIC_Gxy2z_F3z_vrr = PBZ*I_KINETIC_Gxy2z_D2z_vrr+2*oned2z*I_KINETIC_Fxyz_D2z_vrr+2*oned2z*I_KINETIC_Gxy2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr;
    Double I_KINETIC_Gx3z_F3z_vrr = PBZ*I_KINETIC_Gx3z_D2z_vrr+3*oned2z*I_KINETIC_Fx2z_D2z_vrr+2*oned2z*I_KINETIC_Gx3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_KINETIC_G4y_F3z_vrr = PBZ*I_KINETIC_G4y_D2z_vrr+2*oned2z*I_KINETIC_G4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_KINETIC_G3yz_F3z_vrr = PBZ*I_KINETIC_G3yz_D2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+2*oned2z*I_KINETIC_G3yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_KINETIC_G2y2z_F3z_vrr = PBZ*I_KINETIC_G2y2z_D2z_vrr+2*oned2z*I_KINETIC_F2yz_D2z_vrr+2*oned2z*I_KINETIC_G2y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_KINETIC_Gy3z_F3z_vrr = PBZ*I_KINETIC_Gy3z_D2z_vrr+3*oned2z*I_KINETIC_Fy2z_D2z_vrr+2*oned2z*I_KINETIC_Gy3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_KINETIC_G4z_F3z_vrr = PBZ*I_KINETIC_G4z_D2z_vrr+4*oned2z*I_KINETIC_F3z_D2z_vrr+2*oned2z*I_KINETIC_G4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_F_b
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_G_F_b_coefs = beta;
    I_KINETIC_G4x_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_F3x_vrr;
    I_KINETIC_G3xy_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_F3x_vrr;
    I_KINETIC_G3xz_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_F3x_vrr;
    I_KINETIC_G2x2y_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_F3x_vrr;
    I_KINETIC_G2xyz_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_F3x_vrr;
    I_KINETIC_G2x2z_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_F3x_vrr;
    I_KINETIC_Gx3y_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_F3x_vrr;
    I_KINETIC_Gx2yz_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_F3x_vrr;
    I_KINETIC_Gxy2z_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_F3x_vrr;
    I_KINETIC_Gx3z_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_F3x_vrr;
    I_KINETIC_G4y_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_F3x_vrr;
    I_KINETIC_G3yz_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_F3x_vrr;
    I_KINETIC_G2y2z_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_F3x_vrr;
    I_KINETIC_Gy3z_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_F3x_vrr;
    I_KINETIC_G4z_F3x_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_F3x_vrr;
    I_KINETIC_G4x_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_F2xy_vrr;
    I_KINETIC_G3xy_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_F2xy_vrr;
    I_KINETIC_G3xz_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_F2xy_vrr;
    I_KINETIC_G2x2y_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_F2xy_vrr;
    I_KINETIC_G2xyz_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_F2xy_vrr;
    I_KINETIC_G2x2z_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_F2xy_vrr;
    I_KINETIC_Gx3y_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_F2xy_vrr;
    I_KINETIC_Gx2yz_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_F2xy_vrr;
    I_KINETIC_Gxy2z_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_F2xy_vrr;
    I_KINETIC_Gx3z_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_F2xy_vrr;
    I_KINETIC_G4y_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_F2xy_vrr;
    I_KINETIC_G3yz_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_F2xy_vrr;
    I_KINETIC_G2y2z_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_F2xy_vrr;
    I_KINETIC_Gy3z_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_F2xy_vrr;
    I_KINETIC_G4z_F2xy_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_F2xy_vrr;
    I_KINETIC_G4x_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_F2xz_vrr;
    I_KINETIC_G3xy_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_F2xz_vrr;
    I_KINETIC_G3xz_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_F2xz_vrr;
    I_KINETIC_G2x2y_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_F2xz_vrr;
    I_KINETIC_G2xyz_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_F2xz_vrr;
    I_KINETIC_G2x2z_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_F2xz_vrr;
    I_KINETIC_Gx3y_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_F2xz_vrr;
    I_KINETIC_Gx2yz_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_F2xz_vrr;
    I_KINETIC_Gxy2z_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_F2xz_vrr;
    I_KINETIC_Gx3z_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_F2xz_vrr;
    I_KINETIC_G4y_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_F2xz_vrr;
    I_KINETIC_G3yz_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_F2xz_vrr;
    I_KINETIC_G2y2z_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_F2xz_vrr;
    I_KINETIC_Gy3z_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_F2xz_vrr;
    I_KINETIC_G4z_F2xz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_F2xz_vrr;
    I_KINETIC_G4x_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_Fx2y_vrr;
    I_KINETIC_G3xy_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_Fx2y_vrr;
    I_KINETIC_G3xz_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_Fx2y_vrr;
    I_KINETIC_G2x2y_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_Fx2y_vrr;
    I_KINETIC_G2xyz_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_Fx2y_vrr;
    I_KINETIC_G2x2z_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_Fx2y_vrr;
    I_KINETIC_Gx3y_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_Fx2y_vrr;
    I_KINETIC_Gx2yz_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_Fx2y_vrr;
    I_KINETIC_Gxy2z_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_Fx2y_vrr;
    I_KINETIC_Gx3z_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_Fx2y_vrr;
    I_KINETIC_G4y_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_Fx2y_vrr;
    I_KINETIC_G3yz_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_Fx2y_vrr;
    I_KINETIC_G2y2z_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_Fx2y_vrr;
    I_KINETIC_Gy3z_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_Fx2y_vrr;
    I_KINETIC_G4z_Fx2y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_Fx2y_vrr;
    I_KINETIC_G4x_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_Fxyz_vrr;
    I_KINETIC_G3xy_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_Fxyz_vrr;
    I_KINETIC_G3xz_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_Fxyz_vrr;
    I_KINETIC_G2x2y_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_Fxyz_vrr;
    I_KINETIC_G2xyz_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_Fxyz_vrr;
    I_KINETIC_G2x2z_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_Fxyz_vrr;
    I_KINETIC_Gx3y_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_Fxyz_vrr;
    I_KINETIC_Gx2yz_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_Fxyz_vrr;
    I_KINETIC_Gxy2z_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_Fxyz_vrr;
    I_KINETIC_Gx3z_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_Fxyz_vrr;
    I_KINETIC_G4y_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_Fxyz_vrr;
    I_KINETIC_G3yz_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_Fxyz_vrr;
    I_KINETIC_G2y2z_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_Fxyz_vrr;
    I_KINETIC_Gy3z_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_Fxyz_vrr;
    I_KINETIC_G4z_Fxyz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_Fxyz_vrr;
    I_KINETIC_G4x_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_Fx2z_vrr;
    I_KINETIC_G3xy_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_Fx2z_vrr;
    I_KINETIC_G3xz_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_Fx2z_vrr;
    I_KINETIC_G2x2y_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_Fx2z_vrr;
    I_KINETIC_G2xyz_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_Fx2z_vrr;
    I_KINETIC_G2x2z_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_Fx2z_vrr;
    I_KINETIC_Gx3y_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_Fx2z_vrr;
    I_KINETIC_Gx2yz_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_Fx2z_vrr;
    I_KINETIC_Gxy2z_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_Fx2z_vrr;
    I_KINETIC_Gx3z_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_Fx2z_vrr;
    I_KINETIC_G4y_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_Fx2z_vrr;
    I_KINETIC_G3yz_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_Fx2z_vrr;
    I_KINETIC_G2y2z_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_Fx2z_vrr;
    I_KINETIC_Gy3z_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_Fx2z_vrr;
    I_KINETIC_G4z_Fx2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_Fx2z_vrr;
    I_KINETIC_G4x_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_F3y_vrr;
    I_KINETIC_G3xy_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_F3y_vrr;
    I_KINETIC_G3xz_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_F3y_vrr;
    I_KINETIC_G2x2y_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_F3y_vrr;
    I_KINETIC_G2xyz_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_F3y_vrr;
    I_KINETIC_G2x2z_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_F3y_vrr;
    I_KINETIC_Gx3y_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_F3y_vrr;
    I_KINETIC_Gx2yz_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_F3y_vrr;
    I_KINETIC_Gxy2z_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_F3y_vrr;
    I_KINETIC_Gx3z_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_F3y_vrr;
    I_KINETIC_G4y_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_F3y_vrr;
    I_KINETIC_G3yz_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_F3y_vrr;
    I_KINETIC_G2y2z_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_F3y_vrr;
    I_KINETIC_Gy3z_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_F3y_vrr;
    I_KINETIC_G4z_F3y_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_F3y_vrr;
    I_KINETIC_G4x_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_F2yz_vrr;
    I_KINETIC_G3xy_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_F2yz_vrr;
    I_KINETIC_G3xz_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_F2yz_vrr;
    I_KINETIC_G2x2y_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_F2yz_vrr;
    I_KINETIC_G2xyz_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_F2yz_vrr;
    I_KINETIC_G2x2z_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_F2yz_vrr;
    I_KINETIC_Gx3y_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_F2yz_vrr;
    I_KINETIC_Gx2yz_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_F2yz_vrr;
    I_KINETIC_Gxy2z_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_F2yz_vrr;
    I_KINETIC_Gx3z_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_F2yz_vrr;
    I_KINETIC_G4y_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_F2yz_vrr;
    I_KINETIC_G3yz_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_F2yz_vrr;
    I_KINETIC_G2y2z_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_F2yz_vrr;
    I_KINETIC_Gy3z_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_F2yz_vrr;
    I_KINETIC_G4z_F2yz_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_F2yz_vrr;
    I_KINETIC_G4x_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_Fy2z_vrr;
    I_KINETIC_G3xy_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_Fy2z_vrr;
    I_KINETIC_G3xz_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_Fy2z_vrr;
    I_KINETIC_G2x2y_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_Fy2z_vrr;
    I_KINETIC_G2xyz_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_Fy2z_vrr;
    I_KINETIC_G2x2z_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_Fy2z_vrr;
    I_KINETIC_Gx3y_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_Fy2z_vrr;
    I_KINETIC_Gx2yz_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_Fy2z_vrr;
    I_KINETIC_Gxy2z_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_Fy2z_vrr;
    I_KINETIC_Gx3z_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_Fy2z_vrr;
    I_KINETIC_G4y_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_Fy2z_vrr;
    I_KINETIC_G3yz_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_Fy2z_vrr;
    I_KINETIC_G2y2z_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_Fy2z_vrr;
    I_KINETIC_Gy3z_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_Fy2z_vrr;
    I_KINETIC_G4z_Fy2z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_Fy2z_vrr;
    I_KINETIC_G4x_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4x_F3z_vrr;
    I_KINETIC_G3xy_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xy_F3z_vrr;
    I_KINETIC_G3xz_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3xz_F3z_vrr;
    I_KINETIC_G2x2y_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2y_F3z_vrr;
    I_KINETIC_G2xyz_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2xyz_F3z_vrr;
    I_KINETIC_G2x2z_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2x2z_F3z_vrr;
    I_KINETIC_Gx3y_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3y_F3z_vrr;
    I_KINETIC_Gx2yz_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx2yz_F3z_vrr;
    I_KINETIC_Gxy2z_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gxy2z_F3z_vrr;
    I_KINETIC_Gx3z_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gx3z_F3z_vrr;
    I_KINETIC_G4y_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4y_F3z_vrr;
    I_KINETIC_G3yz_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G3yz_F3z_vrr;
    I_KINETIC_G2y2z_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G2y2z_F3z_vrr;
    I_KINETIC_Gy3z_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_Gy3z_F3z_vrr;
    I_KINETIC_G4z_F3z_b += SQ_KINETIC_G_F_b_coefs*I_KINETIC_G4z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_P
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_G4x_Px += I_KINETIC_G4x_Px_vrr;
    I_KINETIC_G3xy_Px += I_KINETIC_G3xy_Px_vrr;
    I_KINETIC_G3xz_Px += I_KINETIC_G3xz_Px_vrr;
    I_KINETIC_G2x2y_Px += I_KINETIC_G2x2y_Px_vrr;
    I_KINETIC_G2xyz_Px += I_KINETIC_G2xyz_Px_vrr;
    I_KINETIC_G2x2z_Px += I_KINETIC_G2x2z_Px_vrr;
    I_KINETIC_Gx3y_Px += I_KINETIC_Gx3y_Px_vrr;
    I_KINETIC_Gx2yz_Px += I_KINETIC_Gx2yz_Px_vrr;
    I_KINETIC_Gxy2z_Px += I_KINETIC_Gxy2z_Px_vrr;
    I_KINETIC_Gx3z_Px += I_KINETIC_Gx3z_Px_vrr;
    I_KINETIC_G4y_Px += I_KINETIC_G4y_Px_vrr;
    I_KINETIC_G3yz_Px += I_KINETIC_G3yz_Px_vrr;
    I_KINETIC_G2y2z_Px += I_KINETIC_G2y2z_Px_vrr;
    I_KINETIC_Gy3z_Px += I_KINETIC_Gy3z_Px_vrr;
    I_KINETIC_G4z_Px += I_KINETIC_G4z_Px_vrr;
    I_KINETIC_G4x_Py += I_KINETIC_G4x_Py_vrr;
    I_KINETIC_G3xy_Py += I_KINETIC_G3xy_Py_vrr;
    I_KINETIC_G3xz_Py += I_KINETIC_G3xz_Py_vrr;
    I_KINETIC_G2x2y_Py += I_KINETIC_G2x2y_Py_vrr;
    I_KINETIC_G2xyz_Py += I_KINETIC_G2xyz_Py_vrr;
    I_KINETIC_G2x2z_Py += I_KINETIC_G2x2z_Py_vrr;
    I_KINETIC_Gx3y_Py += I_KINETIC_Gx3y_Py_vrr;
    I_KINETIC_Gx2yz_Py += I_KINETIC_Gx2yz_Py_vrr;
    I_KINETIC_Gxy2z_Py += I_KINETIC_Gxy2z_Py_vrr;
    I_KINETIC_Gx3z_Py += I_KINETIC_Gx3z_Py_vrr;
    I_KINETIC_G4y_Py += I_KINETIC_G4y_Py_vrr;
    I_KINETIC_G3yz_Py += I_KINETIC_G3yz_Py_vrr;
    I_KINETIC_G2y2z_Py += I_KINETIC_G2y2z_Py_vrr;
    I_KINETIC_Gy3z_Py += I_KINETIC_Gy3z_Py_vrr;
    I_KINETIC_G4z_Py += I_KINETIC_G4z_Py_vrr;
    I_KINETIC_G4x_Pz += I_KINETIC_G4x_Pz_vrr;
    I_KINETIC_G3xy_Pz += I_KINETIC_G3xy_Pz_vrr;
    I_KINETIC_G3xz_Pz += I_KINETIC_G3xz_Pz_vrr;
    I_KINETIC_G2x2y_Pz += I_KINETIC_G2x2y_Pz_vrr;
    I_KINETIC_G2xyz_Pz += I_KINETIC_G2xyz_Pz_vrr;
    I_KINETIC_G2x2z_Pz += I_KINETIC_G2x2z_Pz_vrr;
    I_KINETIC_Gx3y_Pz += I_KINETIC_Gx3y_Pz_vrr;
    I_KINETIC_Gx2yz_Pz += I_KINETIC_Gx2yz_Pz_vrr;
    I_KINETIC_Gxy2z_Pz += I_KINETIC_Gxy2z_Pz_vrr;
    I_KINETIC_Gx3z_Pz += I_KINETIC_Gx3z_Pz_vrr;
    I_KINETIC_G4y_Pz += I_KINETIC_G4y_Pz_vrr;
    I_KINETIC_G3yz_Pz += I_KINETIC_G3yz_Pz_vrr;
    I_KINETIC_G2y2z_Pz += I_KINETIC_G2y2z_Pz_vrr;
    I_KINETIC_Gy3z_Pz += I_KINETIC_Gy3z_Pz_vrr;
    I_KINETIC_G4z_Pz += I_KINETIC_G4z_Pz_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_D_dbx
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_G_F_b
   * RHS shell quartet name: SQ_KINETIC_G_P
   ************************************************************/
  abcd[0] = 2.0E0*I_KINETIC_G4x_F3x_b-2*I_KINETIC_G4x_Px;
  abcd[1] = 2.0E0*I_KINETIC_G3xy_F3x_b-2*I_KINETIC_G3xy_Px;
  abcd[2] = 2.0E0*I_KINETIC_G3xz_F3x_b-2*I_KINETIC_G3xz_Px;
  abcd[3] = 2.0E0*I_KINETIC_G2x2y_F3x_b-2*I_KINETIC_G2x2y_Px;
  abcd[4] = 2.0E0*I_KINETIC_G2xyz_F3x_b-2*I_KINETIC_G2xyz_Px;
  abcd[5] = 2.0E0*I_KINETIC_G2x2z_F3x_b-2*I_KINETIC_G2x2z_Px;
  abcd[6] = 2.0E0*I_KINETIC_Gx3y_F3x_b-2*I_KINETIC_Gx3y_Px;
  abcd[7] = 2.0E0*I_KINETIC_Gx2yz_F3x_b-2*I_KINETIC_Gx2yz_Px;
  abcd[8] = 2.0E0*I_KINETIC_Gxy2z_F3x_b-2*I_KINETIC_Gxy2z_Px;
  abcd[9] = 2.0E0*I_KINETIC_Gx3z_F3x_b-2*I_KINETIC_Gx3z_Px;
  abcd[10] = 2.0E0*I_KINETIC_G4y_F3x_b-2*I_KINETIC_G4y_Px;
  abcd[11] = 2.0E0*I_KINETIC_G3yz_F3x_b-2*I_KINETIC_G3yz_Px;
  abcd[12] = 2.0E0*I_KINETIC_G2y2z_F3x_b-2*I_KINETIC_G2y2z_Px;
  abcd[13] = 2.0E0*I_KINETIC_Gy3z_F3x_b-2*I_KINETIC_Gy3z_Px;
  abcd[14] = 2.0E0*I_KINETIC_G4z_F3x_b-2*I_KINETIC_G4z_Px;
  abcd[15] = 2.0E0*I_KINETIC_G4x_F2xy_b-1*I_KINETIC_G4x_Py;
  abcd[16] = 2.0E0*I_KINETIC_G3xy_F2xy_b-1*I_KINETIC_G3xy_Py;
  abcd[17] = 2.0E0*I_KINETIC_G3xz_F2xy_b-1*I_KINETIC_G3xz_Py;
  abcd[18] = 2.0E0*I_KINETIC_G2x2y_F2xy_b-1*I_KINETIC_G2x2y_Py;
  abcd[19] = 2.0E0*I_KINETIC_G2xyz_F2xy_b-1*I_KINETIC_G2xyz_Py;
  abcd[20] = 2.0E0*I_KINETIC_G2x2z_F2xy_b-1*I_KINETIC_G2x2z_Py;
  abcd[21] = 2.0E0*I_KINETIC_Gx3y_F2xy_b-1*I_KINETIC_Gx3y_Py;
  abcd[22] = 2.0E0*I_KINETIC_Gx2yz_F2xy_b-1*I_KINETIC_Gx2yz_Py;
  abcd[23] = 2.0E0*I_KINETIC_Gxy2z_F2xy_b-1*I_KINETIC_Gxy2z_Py;
  abcd[24] = 2.0E0*I_KINETIC_Gx3z_F2xy_b-1*I_KINETIC_Gx3z_Py;
  abcd[25] = 2.0E0*I_KINETIC_G4y_F2xy_b-1*I_KINETIC_G4y_Py;
  abcd[26] = 2.0E0*I_KINETIC_G3yz_F2xy_b-1*I_KINETIC_G3yz_Py;
  abcd[27] = 2.0E0*I_KINETIC_G2y2z_F2xy_b-1*I_KINETIC_G2y2z_Py;
  abcd[28] = 2.0E0*I_KINETIC_Gy3z_F2xy_b-1*I_KINETIC_Gy3z_Py;
  abcd[29] = 2.0E0*I_KINETIC_G4z_F2xy_b-1*I_KINETIC_G4z_Py;
  abcd[30] = 2.0E0*I_KINETIC_G4x_F2xz_b-1*I_KINETIC_G4x_Pz;
  abcd[31] = 2.0E0*I_KINETIC_G3xy_F2xz_b-1*I_KINETIC_G3xy_Pz;
  abcd[32] = 2.0E0*I_KINETIC_G3xz_F2xz_b-1*I_KINETIC_G3xz_Pz;
  abcd[33] = 2.0E0*I_KINETIC_G2x2y_F2xz_b-1*I_KINETIC_G2x2y_Pz;
  abcd[34] = 2.0E0*I_KINETIC_G2xyz_F2xz_b-1*I_KINETIC_G2xyz_Pz;
  abcd[35] = 2.0E0*I_KINETIC_G2x2z_F2xz_b-1*I_KINETIC_G2x2z_Pz;
  abcd[36] = 2.0E0*I_KINETIC_Gx3y_F2xz_b-1*I_KINETIC_Gx3y_Pz;
  abcd[37] = 2.0E0*I_KINETIC_Gx2yz_F2xz_b-1*I_KINETIC_Gx2yz_Pz;
  abcd[38] = 2.0E0*I_KINETIC_Gxy2z_F2xz_b-1*I_KINETIC_Gxy2z_Pz;
  abcd[39] = 2.0E0*I_KINETIC_Gx3z_F2xz_b-1*I_KINETIC_Gx3z_Pz;
  abcd[40] = 2.0E0*I_KINETIC_G4y_F2xz_b-1*I_KINETIC_G4y_Pz;
  abcd[41] = 2.0E0*I_KINETIC_G3yz_F2xz_b-1*I_KINETIC_G3yz_Pz;
  abcd[42] = 2.0E0*I_KINETIC_G2y2z_F2xz_b-1*I_KINETIC_G2y2z_Pz;
  abcd[43] = 2.0E0*I_KINETIC_Gy3z_F2xz_b-1*I_KINETIC_Gy3z_Pz;
  abcd[44] = 2.0E0*I_KINETIC_G4z_F2xz_b-1*I_KINETIC_G4z_Pz;
  abcd[45] = 2.0E0*I_KINETIC_G4x_Fx2y_b;
  abcd[46] = 2.0E0*I_KINETIC_G3xy_Fx2y_b;
  abcd[47] = 2.0E0*I_KINETIC_G3xz_Fx2y_b;
  abcd[48] = 2.0E0*I_KINETIC_G2x2y_Fx2y_b;
  abcd[49] = 2.0E0*I_KINETIC_G2xyz_Fx2y_b;
  abcd[50] = 2.0E0*I_KINETIC_G2x2z_Fx2y_b;
  abcd[51] = 2.0E0*I_KINETIC_Gx3y_Fx2y_b;
  abcd[52] = 2.0E0*I_KINETIC_Gx2yz_Fx2y_b;
  abcd[53] = 2.0E0*I_KINETIC_Gxy2z_Fx2y_b;
  abcd[54] = 2.0E0*I_KINETIC_Gx3z_Fx2y_b;
  abcd[55] = 2.0E0*I_KINETIC_G4y_Fx2y_b;
  abcd[56] = 2.0E0*I_KINETIC_G3yz_Fx2y_b;
  abcd[57] = 2.0E0*I_KINETIC_G2y2z_Fx2y_b;
  abcd[58] = 2.0E0*I_KINETIC_Gy3z_Fx2y_b;
  abcd[59] = 2.0E0*I_KINETIC_G4z_Fx2y_b;
  abcd[60] = 2.0E0*I_KINETIC_G4x_Fxyz_b;
  abcd[61] = 2.0E0*I_KINETIC_G3xy_Fxyz_b;
  abcd[62] = 2.0E0*I_KINETIC_G3xz_Fxyz_b;
  abcd[63] = 2.0E0*I_KINETIC_G2x2y_Fxyz_b;
  abcd[64] = 2.0E0*I_KINETIC_G2xyz_Fxyz_b;
  abcd[65] = 2.0E0*I_KINETIC_G2x2z_Fxyz_b;
  abcd[66] = 2.0E0*I_KINETIC_Gx3y_Fxyz_b;
  abcd[67] = 2.0E0*I_KINETIC_Gx2yz_Fxyz_b;
  abcd[68] = 2.0E0*I_KINETIC_Gxy2z_Fxyz_b;
  abcd[69] = 2.0E0*I_KINETIC_Gx3z_Fxyz_b;
  abcd[70] = 2.0E0*I_KINETIC_G4y_Fxyz_b;
  abcd[71] = 2.0E0*I_KINETIC_G3yz_Fxyz_b;
  abcd[72] = 2.0E0*I_KINETIC_G2y2z_Fxyz_b;
  abcd[73] = 2.0E0*I_KINETIC_Gy3z_Fxyz_b;
  abcd[74] = 2.0E0*I_KINETIC_G4z_Fxyz_b;
  abcd[75] = 2.0E0*I_KINETIC_G4x_Fx2z_b;
  abcd[76] = 2.0E0*I_KINETIC_G3xy_Fx2z_b;
  abcd[77] = 2.0E0*I_KINETIC_G3xz_Fx2z_b;
  abcd[78] = 2.0E0*I_KINETIC_G2x2y_Fx2z_b;
  abcd[79] = 2.0E0*I_KINETIC_G2xyz_Fx2z_b;
  abcd[80] = 2.0E0*I_KINETIC_G2x2z_Fx2z_b;
  abcd[81] = 2.0E0*I_KINETIC_Gx3y_Fx2z_b;
  abcd[82] = 2.0E0*I_KINETIC_Gx2yz_Fx2z_b;
  abcd[83] = 2.0E0*I_KINETIC_Gxy2z_Fx2z_b;
  abcd[84] = 2.0E0*I_KINETIC_Gx3z_Fx2z_b;
  abcd[85] = 2.0E0*I_KINETIC_G4y_Fx2z_b;
  abcd[86] = 2.0E0*I_KINETIC_G3yz_Fx2z_b;
  abcd[87] = 2.0E0*I_KINETIC_G2y2z_Fx2z_b;
  abcd[88] = 2.0E0*I_KINETIC_Gy3z_Fx2z_b;
  abcd[89] = 2.0E0*I_KINETIC_G4z_Fx2z_b;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_D_dby
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_G_F_b
   * RHS shell quartet name: SQ_KINETIC_G_P
   ************************************************************/
  abcd[90] = 2.0E0*I_KINETIC_G4x_F2xy_b;
  abcd[91] = 2.0E0*I_KINETIC_G3xy_F2xy_b;
  abcd[92] = 2.0E0*I_KINETIC_G3xz_F2xy_b;
  abcd[93] = 2.0E0*I_KINETIC_G2x2y_F2xy_b;
  abcd[94] = 2.0E0*I_KINETIC_G2xyz_F2xy_b;
  abcd[95] = 2.0E0*I_KINETIC_G2x2z_F2xy_b;
  abcd[96] = 2.0E0*I_KINETIC_Gx3y_F2xy_b;
  abcd[97] = 2.0E0*I_KINETIC_Gx2yz_F2xy_b;
  abcd[98] = 2.0E0*I_KINETIC_Gxy2z_F2xy_b;
  abcd[99] = 2.0E0*I_KINETIC_Gx3z_F2xy_b;
  abcd[100] = 2.0E0*I_KINETIC_G4y_F2xy_b;
  abcd[101] = 2.0E0*I_KINETIC_G3yz_F2xy_b;
  abcd[102] = 2.0E0*I_KINETIC_G2y2z_F2xy_b;
  abcd[103] = 2.0E0*I_KINETIC_Gy3z_F2xy_b;
  abcd[104] = 2.0E0*I_KINETIC_G4z_F2xy_b;
  abcd[105] = 2.0E0*I_KINETIC_G4x_Fx2y_b-1*I_KINETIC_G4x_Px;
  abcd[106] = 2.0E0*I_KINETIC_G3xy_Fx2y_b-1*I_KINETIC_G3xy_Px;
  abcd[107] = 2.0E0*I_KINETIC_G3xz_Fx2y_b-1*I_KINETIC_G3xz_Px;
  abcd[108] = 2.0E0*I_KINETIC_G2x2y_Fx2y_b-1*I_KINETIC_G2x2y_Px;
  abcd[109] = 2.0E0*I_KINETIC_G2xyz_Fx2y_b-1*I_KINETIC_G2xyz_Px;
  abcd[110] = 2.0E0*I_KINETIC_G2x2z_Fx2y_b-1*I_KINETIC_G2x2z_Px;
  abcd[111] = 2.0E0*I_KINETIC_Gx3y_Fx2y_b-1*I_KINETIC_Gx3y_Px;
  abcd[112] = 2.0E0*I_KINETIC_Gx2yz_Fx2y_b-1*I_KINETIC_Gx2yz_Px;
  abcd[113] = 2.0E0*I_KINETIC_Gxy2z_Fx2y_b-1*I_KINETIC_Gxy2z_Px;
  abcd[114] = 2.0E0*I_KINETIC_Gx3z_Fx2y_b-1*I_KINETIC_Gx3z_Px;
  abcd[115] = 2.0E0*I_KINETIC_G4y_Fx2y_b-1*I_KINETIC_G4y_Px;
  abcd[116] = 2.0E0*I_KINETIC_G3yz_Fx2y_b-1*I_KINETIC_G3yz_Px;
  abcd[117] = 2.0E0*I_KINETIC_G2y2z_Fx2y_b-1*I_KINETIC_G2y2z_Px;
  abcd[118] = 2.0E0*I_KINETIC_Gy3z_Fx2y_b-1*I_KINETIC_Gy3z_Px;
  abcd[119] = 2.0E0*I_KINETIC_G4z_Fx2y_b-1*I_KINETIC_G4z_Px;
  abcd[120] = 2.0E0*I_KINETIC_G4x_Fxyz_b;
  abcd[121] = 2.0E0*I_KINETIC_G3xy_Fxyz_b;
  abcd[122] = 2.0E0*I_KINETIC_G3xz_Fxyz_b;
  abcd[123] = 2.0E0*I_KINETIC_G2x2y_Fxyz_b;
  abcd[124] = 2.0E0*I_KINETIC_G2xyz_Fxyz_b;
  abcd[125] = 2.0E0*I_KINETIC_G2x2z_Fxyz_b;
  abcd[126] = 2.0E0*I_KINETIC_Gx3y_Fxyz_b;
  abcd[127] = 2.0E0*I_KINETIC_Gx2yz_Fxyz_b;
  abcd[128] = 2.0E0*I_KINETIC_Gxy2z_Fxyz_b;
  abcd[129] = 2.0E0*I_KINETIC_Gx3z_Fxyz_b;
  abcd[130] = 2.0E0*I_KINETIC_G4y_Fxyz_b;
  abcd[131] = 2.0E0*I_KINETIC_G3yz_Fxyz_b;
  abcd[132] = 2.0E0*I_KINETIC_G2y2z_Fxyz_b;
  abcd[133] = 2.0E0*I_KINETIC_Gy3z_Fxyz_b;
  abcd[134] = 2.0E0*I_KINETIC_G4z_Fxyz_b;
  abcd[135] = 2.0E0*I_KINETIC_G4x_F3y_b-2*I_KINETIC_G4x_Py;
  abcd[136] = 2.0E0*I_KINETIC_G3xy_F3y_b-2*I_KINETIC_G3xy_Py;
  abcd[137] = 2.0E0*I_KINETIC_G3xz_F3y_b-2*I_KINETIC_G3xz_Py;
  abcd[138] = 2.0E0*I_KINETIC_G2x2y_F3y_b-2*I_KINETIC_G2x2y_Py;
  abcd[139] = 2.0E0*I_KINETIC_G2xyz_F3y_b-2*I_KINETIC_G2xyz_Py;
  abcd[140] = 2.0E0*I_KINETIC_G2x2z_F3y_b-2*I_KINETIC_G2x2z_Py;
  abcd[141] = 2.0E0*I_KINETIC_Gx3y_F3y_b-2*I_KINETIC_Gx3y_Py;
  abcd[142] = 2.0E0*I_KINETIC_Gx2yz_F3y_b-2*I_KINETIC_Gx2yz_Py;
  abcd[143] = 2.0E0*I_KINETIC_Gxy2z_F3y_b-2*I_KINETIC_Gxy2z_Py;
  abcd[144] = 2.0E0*I_KINETIC_Gx3z_F3y_b-2*I_KINETIC_Gx3z_Py;
  abcd[145] = 2.0E0*I_KINETIC_G4y_F3y_b-2*I_KINETIC_G4y_Py;
  abcd[146] = 2.0E0*I_KINETIC_G3yz_F3y_b-2*I_KINETIC_G3yz_Py;
  abcd[147] = 2.0E0*I_KINETIC_G2y2z_F3y_b-2*I_KINETIC_G2y2z_Py;
  abcd[148] = 2.0E0*I_KINETIC_Gy3z_F3y_b-2*I_KINETIC_Gy3z_Py;
  abcd[149] = 2.0E0*I_KINETIC_G4z_F3y_b-2*I_KINETIC_G4z_Py;
  abcd[150] = 2.0E0*I_KINETIC_G4x_F2yz_b-1*I_KINETIC_G4x_Pz;
  abcd[151] = 2.0E0*I_KINETIC_G3xy_F2yz_b-1*I_KINETIC_G3xy_Pz;
  abcd[152] = 2.0E0*I_KINETIC_G3xz_F2yz_b-1*I_KINETIC_G3xz_Pz;
  abcd[153] = 2.0E0*I_KINETIC_G2x2y_F2yz_b-1*I_KINETIC_G2x2y_Pz;
  abcd[154] = 2.0E0*I_KINETIC_G2xyz_F2yz_b-1*I_KINETIC_G2xyz_Pz;
  abcd[155] = 2.0E0*I_KINETIC_G2x2z_F2yz_b-1*I_KINETIC_G2x2z_Pz;
  abcd[156] = 2.0E0*I_KINETIC_Gx3y_F2yz_b-1*I_KINETIC_Gx3y_Pz;
  abcd[157] = 2.0E0*I_KINETIC_Gx2yz_F2yz_b-1*I_KINETIC_Gx2yz_Pz;
  abcd[158] = 2.0E0*I_KINETIC_Gxy2z_F2yz_b-1*I_KINETIC_Gxy2z_Pz;
  abcd[159] = 2.0E0*I_KINETIC_Gx3z_F2yz_b-1*I_KINETIC_Gx3z_Pz;
  abcd[160] = 2.0E0*I_KINETIC_G4y_F2yz_b-1*I_KINETIC_G4y_Pz;
  abcd[161] = 2.0E0*I_KINETIC_G3yz_F2yz_b-1*I_KINETIC_G3yz_Pz;
  abcd[162] = 2.0E0*I_KINETIC_G2y2z_F2yz_b-1*I_KINETIC_G2y2z_Pz;
  abcd[163] = 2.0E0*I_KINETIC_Gy3z_F2yz_b-1*I_KINETIC_Gy3z_Pz;
  abcd[164] = 2.0E0*I_KINETIC_G4z_F2yz_b-1*I_KINETIC_G4z_Pz;
  abcd[165] = 2.0E0*I_KINETIC_G4x_Fy2z_b;
  abcd[166] = 2.0E0*I_KINETIC_G3xy_Fy2z_b;
  abcd[167] = 2.0E0*I_KINETIC_G3xz_Fy2z_b;
  abcd[168] = 2.0E0*I_KINETIC_G2x2y_Fy2z_b;
  abcd[169] = 2.0E0*I_KINETIC_G2xyz_Fy2z_b;
  abcd[170] = 2.0E0*I_KINETIC_G2x2z_Fy2z_b;
  abcd[171] = 2.0E0*I_KINETIC_Gx3y_Fy2z_b;
  abcd[172] = 2.0E0*I_KINETIC_Gx2yz_Fy2z_b;
  abcd[173] = 2.0E0*I_KINETIC_Gxy2z_Fy2z_b;
  abcd[174] = 2.0E0*I_KINETIC_Gx3z_Fy2z_b;
  abcd[175] = 2.0E0*I_KINETIC_G4y_Fy2z_b;
  abcd[176] = 2.0E0*I_KINETIC_G3yz_Fy2z_b;
  abcd[177] = 2.0E0*I_KINETIC_G2y2z_Fy2z_b;
  abcd[178] = 2.0E0*I_KINETIC_Gy3z_Fy2z_b;
  abcd[179] = 2.0E0*I_KINETIC_G4z_Fy2z_b;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_D_dbz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_G_F_b
   * RHS shell quartet name: SQ_KINETIC_G_P
   ************************************************************/
  abcd[180] = 2.0E0*I_KINETIC_G4x_F2xz_b;
  abcd[181] = 2.0E0*I_KINETIC_G3xy_F2xz_b;
  abcd[182] = 2.0E0*I_KINETIC_G3xz_F2xz_b;
  abcd[183] = 2.0E0*I_KINETIC_G2x2y_F2xz_b;
  abcd[184] = 2.0E0*I_KINETIC_G2xyz_F2xz_b;
  abcd[185] = 2.0E0*I_KINETIC_G2x2z_F2xz_b;
  abcd[186] = 2.0E0*I_KINETIC_Gx3y_F2xz_b;
  abcd[187] = 2.0E0*I_KINETIC_Gx2yz_F2xz_b;
  abcd[188] = 2.0E0*I_KINETIC_Gxy2z_F2xz_b;
  abcd[189] = 2.0E0*I_KINETIC_Gx3z_F2xz_b;
  abcd[190] = 2.0E0*I_KINETIC_G4y_F2xz_b;
  abcd[191] = 2.0E0*I_KINETIC_G3yz_F2xz_b;
  abcd[192] = 2.0E0*I_KINETIC_G2y2z_F2xz_b;
  abcd[193] = 2.0E0*I_KINETIC_Gy3z_F2xz_b;
  abcd[194] = 2.0E0*I_KINETIC_G4z_F2xz_b;
  abcd[195] = 2.0E0*I_KINETIC_G4x_Fxyz_b;
  abcd[196] = 2.0E0*I_KINETIC_G3xy_Fxyz_b;
  abcd[197] = 2.0E0*I_KINETIC_G3xz_Fxyz_b;
  abcd[198] = 2.0E0*I_KINETIC_G2x2y_Fxyz_b;
  abcd[199] = 2.0E0*I_KINETIC_G2xyz_Fxyz_b;
  abcd[200] = 2.0E0*I_KINETIC_G2x2z_Fxyz_b;
  abcd[201] = 2.0E0*I_KINETIC_Gx3y_Fxyz_b;
  abcd[202] = 2.0E0*I_KINETIC_Gx2yz_Fxyz_b;
  abcd[203] = 2.0E0*I_KINETIC_Gxy2z_Fxyz_b;
  abcd[204] = 2.0E0*I_KINETIC_Gx3z_Fxyz_b;
  abcd[205] = 2.0E0*I_KINETIC_G4y_Fxyz_b;
  abcd[206] = 2.0E0*I_KINETIC_G3yz_Fxyz_b;
  abcd[207] = 2.0E0*I_KINETIC_G2y2z_Fxyz_b;
  abcd[208] = 2.0E0*I_KINETIC_Gy3z_Fxyz_b;
  abcd[209] = 2.0E0*I_KINETIC_G4z_Fxyz_b;
  abcd[210] = 2.0E0*I_KINETIC_G4x_Fx2z_b-1*I_KINETIC_G4x_Px;
  abcd[211] = 2.0E0*I_KINETIC_G3xy_Fx2z_b-1*I_KINETIC_G3xy_Px;
  abcd[212] = 2.0E0*I_KINETIC_G3xz_Fx2z_b-1*I_KINETIC_G3xz_Px;
  abcd[213] = 2.0E0*I_KINETIC_G2x2y_Fx2z_b-1*I_KINETIC_G2x2y_Px;
  abcd[214] = 2.0E0*I_KINETIC_G2xyz_Fx2z_b-1*I_KINETIC_G2xyz_Px;
  abcd[215] = 2.0E0*I_KINETIC_G2x2z_Fx2z_b-1*I_KINETIC_G2x2z_Px;
  abcd[216] = 2.0E0*I_KINETIC_Gx3y_Fx2z_b-1*I_KINETIC_Gx3y_Px;
  abcd[217] = 2.0E0*I_KINETIC_Gx2yz_Fx2z_b-1*I_KINETIC_Gx2yz_Px;
  abcd[218] = 2.0E0*I_KINETIC_Gxy2z_Fx2z_b-1*I_KINETIC_Gxy2z_Px;
  abcd[219] = 2.0E0*I_KINETIC_Gx3z_Fx2z_b-1*I_KINETIC_Gx3z_Px;
  abcd[220] = 2.0E0*I_KINETIC_G4y_Fx2z_b-1*I_KINETIC_G4y_Px;
  abcd[221] = 2.0E0*I_KINETIC_G3yz_Fx2z_b-1*I_KINETIC_G3yz_Px;
  abcd[222] = 2.0E0*I_KINETIC_G2y2z_Fx2z_b-1*I_KINETIC_G2y2z_Px;
  abcd[223] = 2.0E0*I_KINETIC_Gy3z_Fx2z_b-1*I_KINETIC_Gy3z_Px;
  abcd[224] = 2.0E0*I_KINETIC_G4z_Fx2z_b-1*I_KINETIC_G4z_Px;
  abcd[225] = 2.0E0*I_KINETIC_G4x_F2yz_b;
  abcd[226] = 2.0E0*I_KINETIC_G3xy_F2yz_b;
  abcd[227] = 2.0E0*I_KINETIC_G3xz_F2yz_b;
  abcd[228] = 2.0E0*I_KINETIC_G2x2y_F2yz_b;
  abcd[229] = 2.0E0*I_KINETIC_G2xyz_F2yz_b;
  abcd[230] = 2.0E0*I_KINETIC_G2x2z_F2yz_b;
  abcd[231] = 2.0E0*I_KINETIC_Gx3y_F2yz_b;
  abcd[232] = 2.0E0*I_KINETIC_Gx2yz_F2yz_b;
  abcd[233] = 2.0E0*I_KINETIC_Gxy2z_F2yz_b;
  abcd[234] = 2.0E0*I_KINETIC_Gx3z_F2yz_b;
  abcd[235] = 2.0E0*I_KINETIC_G4y_F2yz_b;
  abcd[236] = 2.0E0*I_KINETIC_G3yz_F2yz_b;
  abcd[237] = 2.0E0*I_KINETIC_G2y2z_F2yz_b;
  abcd[238] = 2.0E0*I_KINETIC_Gy3z_F2yz_b;
  abcd[239] = 2.0E0*I_KINETIC_G4z_F2yz_b;
  abcd[240] = 2.0E0*I_KINETIC_G4x_Fy2z_b-1*I_KINETIC_G4x_Py;
  abcd[241] = 2.0E0*I_KINETIC_G3xy_Fy2z_b-1*I_KINETIC_G3xy_Py;
  abcd[242] = 2.0E0*I_KINETIC_G3xz_Fy2z_b-1*I_KINETIC_G3xz_Py;
  abcd[243] = 2.0E0*I_KINETIC_G2x2y_Fy2z_b-1*I_KINETIC_G2x2y_Py;
  abcd[244] = 2.0E0*I_KINETIC_G2xyz_Fy2z_b-1*I_KINETIC_G2xyz_Py;
  abcd[245] = 2.0E0*I_KINETIC_G2x2z_Fy2z_b-1*I_KINETIC_G2x2z_Py;
  abcd[246] = 2.0E0*I_KINETIC_Gx3y_Fy2z_b-1*I_KINETIC_Gx3y_Py;
  abcd[247] = 2.0E0*I_KINETIC_Gx2yz_Fy2z_b-1*I_KINETIC_Gx2yz_Py;
  abcd[248] = 2.0E0*I_KINETIC_Gxy2z_Fy2z_b-1*I_KINETIC_Gxy2z_Py;
  abcd[249] = 2.0E0*I_KINETIC_Gx3z_Fy2z_b-1*I_KINETIC_Gx3z_Py;
  abcd[250] = 2.0E0*I_KINETIC_G4y_Fy2z_b-1*I_KINETIC_G4y_Py;
  abcd[251] = 2.0E0*I_KINETIC_G3yz_Fy2z_b-1*I_KINETIC_G3yz_Py;
  abcd[252] = 2.0E0*I_KINETIC_G2y2z_Fy2z_b-1*I_KINETIC_G2y2z_Py;
  abcd[253] = 2.0E0*I_KINETIC_Gy3z_Fy2z_b-1*I_KINETIC_Gy3z_Py;
  abcd[254] = 2.0E0*I_KINETIC_G4z_Fy2z_b-1*I_KINETIC_G4z_Py;
  abcd[255] = 2.0E0*I_KINETIC_G4x_F3z_b-2*I_KINETIC_G4x_Pz;
  abcd[256] = 2.0E0*I_KINETIC_G3xy_F3z_b-2*I_KINETIC_G3xy_Pz;
  abcd[257] = 2.0E0*I_KINETIC_G3xz_F3z_b-2*I_KINETIC_G3xz_Pz;
  abcd[258] = 2.0E0*I_KINETIC_G2x2y_F3z_b-2*I_KINETIC_G2x2y_Pz;
  abcd[259] = 2.0E0*I_KINETIC_G2xyz_F3z_b-2*I_KINETIC_G2xyz_Pz;
  abcd[260] = 2.0E0*I_KINETIC_G2x2z_F3z_b-2*I_KINETIC_G2x2z_Pz;
  abcd[261] = 2.0E0*I_KINETIC_Gx3y_F3z_b-2*I_KINETIC_Gx3y_Pz;
  abcd[262] = 2.0E0*I_KINETIC_Gx2yz_F3z_b-2*I_KINETIC_Gx2yz_Pz;
  abcd[263] = 2.0E0*I_KINETIC_Gxy2z_F3z_b-2*I_KINETIC_Gxy2z_Pz;
  abcd[264] = 2.0E0*I_KINETIC_Gx3z_F3z_b-2*I_KINETIC_Gx3z_Pz;
  abcd[265] = 2.0E0*I_KINETIC_G4y_F3z_b-2*I_KINETIC_G4y_Pz;
  abcd[266] = 2.0E0*I_KINETIC_G3yz_F3z_b-2*I_KINETIC_G3yz_Pz;
  abcd[267] = 2.0E0*I_KINETIC_G2y2z_F3z_b-2*I_KINETIC_G2y2z_Pz;
  abcd[268] = 2.0E0*I_KINETIC_Gy3z_F3z_b-2*I_KINETIC_Gy3z_Pz;
  abcd[269] = 2.0E0*I_KINETIC_G4z_F3z_b-2*I_KINETIC_G4z_Pz;
}
