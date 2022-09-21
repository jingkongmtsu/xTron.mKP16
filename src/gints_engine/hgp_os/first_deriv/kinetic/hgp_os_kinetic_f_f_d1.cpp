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

void hgp_os_kinetic_f_f_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_G4x_F3x_a = 0.0E0;
  Double I_KINETIC_G3xy_F3x_a = 0.0E0;
  Double I_KINETIC_G3xz_F3x_a = 0.0E0;
  Double I_KINETIC_G2x2y_F3x_a = 0.0E0;
  Double I_KINETIC_G2xyz_F3x_a = 0.0E0;
  Double I_KINETIC_G2x2z_F3x_a = 0.0E0;
  Double I_KINETIC_Gx3y_F3x_a = 0.0E0;
  Double I_KINETIC_Gx2yz_F3x_a = 0.0E0;
  Double I_KINETIC_Gxy2z_F3x_a = 0.0E0;
  Double I_KINETIC_Gx3z_F3x_a = 0.0E0;
  Double I_KINETIC_G4y_F3x_a = 0.0E0;
  Double I_KINETIC_G3yz_F3x_a = 0.0E0;
  Double I_KINETIC_G2y2z_F3x_a = 0.0E0;
  Double I_KINETIC_Gy3z_F3x_a = 0.0E0;
  Double I_KINETIC_G4z_F3x_a = 0.0E0;
  Double I_KINETIC_G4x_F2xy_a = 0.0E0;
  Double I_KINETIC_G3xy_F2xy_a = 0.0E0;
  Double I_KINETIC_G3xz_F2xy_a = 0.0E0;
  Double I_KINETIC_G2x2y_F2xy_a = 0.0E0;
  Double I_KINETIC_G2xyz_F2xy_a = 0.0E0;
  Double I_KINETIC_G2x2z_F2xy_a = 0.0E0;
  Double I_KINETIC_Gx3y_F2xy_a = 0.0E0;
  Double I_KINETIC_Gx2yz_F2xy_a = 0.0E0;
  Double I_KINETIC_Gxy2z_F2xy_a = 0.0E0;
  Double I_KINETIC_Gx3z_F2xy_a = 0.0E0;
  Double I_KINETIC_G4y_F2xy_a = 0.0E0;
  Double I_KINETIC_G3yz_F2xy_a = 0.0E0;
  Double I_KINETIC_G2y2z_F2xy_a = 0.0E0;
  Double I_KINETIC_Gy3z_F2xy_a = 0.0E0;
  Double I_KINETIC_G4z_F2xy_a = 0.0E0;
  Double I_KINETIC_G4x_F2xz_a = 0.0E0;
  Double I_KINETIC_G3xy_F2xz_a = 0.0E0;
  Double I_KINETIC_G3xz_F2xz_a = 0.0E0;
  Double I_KINETIC_G2x2y_F2xz_a = 0.0E0;
  Double I_KINETIC_G2xyz_F2xz_a = 0.0E0;
  Double I_KINETIC_G2x2z_F2xz_a = 0.0E0;
  Double I_KINETIC_Gx3y_F2xz_a = 0.0E0;
  Double I_KINETIC_Gx2yz_F2xz_a = 0.0E0;
  Double I_KINETIC_Gxy2z_F2xz_a = 0.0E0;
  Double I_KINETIC_Gx3z_F2xz_a = 0.0E0;
  Double I_KINETIC_G4y_F2xz_a = 0.0E0;
  Double I_KINETIC_G3yz_F2xz_a = 0.0E0;
  Double I_KINETIC_G2y2z_F2xz_a = 0.0E0;
  Double I_KINETIC_Gy3z_F2xz_a = 0.0E0;
  Double I_KINETIC_G4z_F2xz_a = 0.0E0;
  Double I_KINETIC_G4x_Fx2y_a = 0.0E0;
  Double I_KINETIC_G3xy_Fx2y_a = 0.0E0;
  Double I_KINETIC_G3xz_Fx2y_a = 0.0E0;
  Double I_KINETIC_G2x2y_Fx2y_a = 0.0E0;
  Double I_KINETIC_G2xyz_Fx2y_a = 0.0E0;
  Double I_KINETIC_G2x2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Gx3y_Fx2y_a = 0.0E0;
  Double I_KINETIC_Gx2yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_Gxy2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Gx3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_G4y_Fx2y_a = 0.0E0;
  Double I_KINETIC_G3yz_Fx2y_a = 0.0E0;
  Double I_KINETIC_G2y2z_Fx2y_a = 0.0E0;
  Double I_KINETIC_Gy3z_Fx2y_a = 0.0E0;
  Double I_KINETIC_G4z_Fx2y_a = 0.0E0;
  Double I_KINETIC_G4x_Fxyz_a = 0.0E0;
  Double I_KINETIC_G3xy_Fxyz_a = 0.0E0;
  Double I_KINETIC_G3xz_Fxyz_a = 0.0E0;
  Double I_KINETIC_G2x2y_Fxyz_a = 0.0E0;
  Double I_KINETIC_G2xyz_Fxyz_a = 0.0E0;
  Double I_KINETIC_G2x2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Gx3y_Fxyz_a = 0.0E0;
  Double I_KINETIC_Gx2yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_Gxy2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Gx3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_G4y_Fxyz_a = 0.0E0;
  Double I_KINETIC_G3yz_Fxyz_a = 0.0E0;
  Double I_KINETIC_G2y2z_Fxyz_a = 0.0E0;
  Double I_KINETIC_Gy3z_Fxyz_a = 0.0E0;
  Double I_KINETIC_G4z_Fxyz_a = 0.0E0;
  Double I_KINETIC_G4x_Fx2z_a = 0.0E0;
  Double I_KINETIC_G3xy_Fx2z_a = 0.0E0;
  Double I_KINETIC_G3xz_Fx2z_a = 0.0E0;
  Double I_KINETIC_G2x2y_Fx2z_a = 0.0E0;
  Double I_KINETIC_G2xyz_Fx2z_a = 0.0E0;
  Double I_KINETIC_G2x2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Gx3y_Fx2z_a = 0.0E0;
  Double I_KINETIC_Gx2yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_Gxy2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Gx3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_G4y_Fx2z_a = 0.0E0;
  Double I_KINETIC_G3yz_Fx2z_a = 0.0E0;
  Double I_KINETIC_G2y2z_Fx2z_a = 0.0E0;
  Double I_KINETIC_Gy3z_Fx2z_a = 0.0E0;
  Double I_KINETIC_G4z_Fx2z_a = 0.0E0;
  Double I_KINETIC_G4x_F3y_a = 0.0E0;
  Double I_KINETIC_G3xy_F3y_a = 0.0E0;
  Double I_KINETIC_G3xz_F3y_a = 0.0E0;
  Double I_KINETIC_G2x2y_F3y_a = 0.0E0;
  Double I_KINETIC_G2xyz_F3y_a = 0.0E0;
  Double I_KINETIC_G2x2z_F3y_a = 0.0E0;
  Double I_KINETIC_Gx3y_F3y_a = 0.0E0;
  Double I_KINETIC_Gx2yz_F3y_a = 0.0E0;
  Double I_KINETIC_Gxy2z_F3y_a = 0.0E0;
  Double I_KINETIC_Gx3z_F3y_a = 0.0E0;
  Double I_KINETIC_G4y_F3y_a = 0.0E0;
  Double I_KINETIC_G3yz_F3y_a = 0.0E0;
  Double I_KINETIC_G2y2z_F3y_a = 0.0E0;
  Double I_KINETIC_Gy3z_F3y_a = 0.0E0;
  Double I_KINETIC_G4z_F3y_a = 0.0E0;
  Double I_KINETIC_G4x_F2yz_a = 0.0E0;
  Double I_KINETIC_G3xy_F2yz_a = 0.0E0;
  Double I_KINETIC_G3xz_F2yz_a = 0.0E0;
  Double I_KINETIC_G2x2y_F2yz_a = 0.0E0;
  Double I_KINETIC_G2xyz_F2yz_a = 0.0E0;
  Double I_KINETIC_G2x2z_F2yz_a = 0.0E0;
  Double I_KINETIC_Gx3y_F2yz_a = 0.0E0;
  Double I_KINETIC_Gx2yz_F2yz_a = 0.0E0;
  Double I_KINETIC_Gxy2z_F2yz_a = 0.0E0;
  Double I_KINETIC_Gx3z_F2yz_a = 0.0E0;
  Double I_KINETIC_G4y_F2yz_a = 0.0E0;
  Double I_KINETIC_G3yz_F2yz_a = 0.0E0;
  Double I_KINETIC_G2y2z_F2yz_a = 0.0E0;
  Double I_KINETIC_Gy3z_F2yz_a = 0.0E0;
  Double I_KINETIC_G4z_F2yz_a = 0.0E0;
  Double I_KINETIC_G4x_Fy2z_a = 0.0E0;
  Double I_KINETIC_G3xy_Fy2z_a = 0.0E0;
  Double I_KINETIC_G3xz_Fy2z_a = 0.0E0;
  Double I_KINETIC_G2x2y_Fy2z_a = 0.0E0;
  Double I_KINETIC_G2xyz_Fy2z_a = 0.0E0;
  Double I_KINETIC_G2x2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Gx3y_Fy2z_a = 0.0E0;
  Double I_KINETIC_Gx2yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_Gxy2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Gx3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_G4y_Fy2z_a = 0.0E0;
  Double I_KINETIC_G3yz_Fy2z_a = 0.0E0;
  Double I_KINETIC_G2y2z_Fy2z_a = 0.0E0;
  Double I_KINETIC_Gy3z_Fy2z_a = 0.0E0;
  Double I_KINETIC_G4z_Fy2z_a = 0.0E0;
  Double I_KINETIC_G4x_F3z_a = 0.0E0;
  Double I_KINETIC_G3xy_F3z_a = 0.0E0;
  Double I_KINETIC_G3xz_F3z_a = 0.0E0;
  Double I_KINETIC_G2x2y_F3z_a = 0.0E0;
  Double I_KINETIC_G2xyz_F3z_a = 0.0E0;
  Double I_KINETIC_G2x2z_F3z_a = 0.0E0;
  Double I_KINETIC_Gx3y_F3z_a = 0.0E0;
  Double I_KINETIC_Gx2yz_F3z_a = 0.0E0;
  Double I_KINETIC_Gxy2z_F3z_a = 0.0E0;
  Double I_KINETIC_Gx3z_F3z_a = 0.0E0;
  Double I_KINETIC_G4y_F3z_a = 0.0E0;
  Double I_KINETIC_G3yz_F3z_a = 0.0E0;
  Double I_KINETIC_G2y2z_F3z_a = 0.0E0;
  Double I_KINETIC_Gy3z_F3z_a = 0.0E0;
  Double I_KINETIC_G4z_F3z_a = 0.0E0;
  Double I_KINETIC_D2x_F3x = 0.0E0;
  Double I_KINETIC_Dxy_F3x = 0.0E0;
  Double I_KINETIC_Dxz_F3x = 0.0E0;
  Double I_KINETIC_D2y_F3x = 0.0E0;
  Double I_KINETIC_Dyz_F3x = 0.0E0;
  Double I_KINETIC_D2z_F3x = 0.0E0;
  Double I_KINETIC_D2x_F2xy = 0.0E0;
  Double I_KINETIC_Dxy_F2xy = 0.0E0;
  Double I_KINETIC_Dxz_F2xy = 0.0E0;
  Double I_KINETIC_D2y_F2xy = 0.0E0;
  Double I_KINETIC_Dyz_F2xy = 0.0E0;
  Double I_KINETIC_D2z_F2xy = 0.0E0;
  Double I_KINETIC_D2x_F2xz = 0.0E0;
  Double I_KINETIC_Dxy_F2xz = 0.0E0;
  Double I_KINETIC_Dxz_F2xz = 0.0E0;
  Double I_KINETIC_D2y_F2xz = 0.0E0;
  Double I_KINETIC_Dyz_F2xz = 0.0E0;
  Double I_KINETIC_D2z_F2xz = 0.0E0;
  Double I_KINETIC_D2x_Fx2y = 0.0E0;
  Double I_KINETIC_Dxy_Fx2y = 0.0E0;
  Double I_KINETIC_Dxz_Fx2y = 0.0E0;
  Double I_KINETIC_D2y_Fx2y = 0.0E0;
  Double I_KINETIC_Dyz_Fx2y = 0.0E0;
  Double I_KINETIC_D2z_Fx2y = 0.0E0;
  Double I_KINETIC_D2x_Fxyz = 0.0E0;
  Double I_KINETIC_Dxy_Fxyz = 0.0E0;
  Double I_KINETIC_Dxz_Fxyz = 0.0E0;
  Double I_KINETIC_D2y_Fxyz = 0.0E0;
  Double I_KINETIC_Dyz_Fxyz = 0.0E0;
  Double I_KINETIC_D2z_Fxyz = 0.0E0;
  Double I_KINETIC_D2x_Fx2z = 0.0E0;
  Double I_KINETIC_Dxy_Fx2z = 0.0E0;
  Double I_KINETIC_Dxz_Fx2z = 0.0E0;
  Double I_KINETIC_D2y_Fx2z = 0.0E0;
  Double I_KINETIC_Dyz_Fx2z = 0.0E0;
  Double I_KINETIC_D2z_Fx2z = 0.0E0;
  Double I_KINETIC_D2x_F3y = 0.0E0;
  Double I_KINETIC_Dxy_F3y = 0.0E0;
  Double I_KINETIC_Dxz_F3y = 0.0E0;
  Double I_KINETIC_D2y_F3y = 0.0E0;
  Double I_KINETIC_Dyz_F3y = 0.0E0;
  Double I_KINETIC_D2z_F3y = 0.0E0;
  Double I_KINETIC_D2x_F2yz = 0.0E0;
  Double I_KINETIC_Dxy_F2yz = 0.0E0;
  Double I_KINETIC_Dxz_F2yz = 0.0E0;
  Double I_KINETIC_D2y_F2yz = 0.0E0;
  Double I_KINETIC_Dyz_F2yz = 0.0E0;
  Double I_KINETIC_D2z_F2yz = 0.0E0;
  Double I_KINETIC_D2x_Fy2z = 0.0E0;
  Double I_KINETIC_Dxy_Fy2z = 0.0E0;
  Double I_KINETIC_Dxz_Fy2z = 0.0E0;
  Double I_KINETIC_D2y_Fy2z = 0.0E0;
  Double I_KINETIC_Dyz_Fy2z = 0.0E0;
  Double I_KINETIC_D2z_Fy2z = 0.0E0;
  Double I_KINETIC_D2x_F3z = 0.0E0;
  Double I_KINETIC_Dxy_F3z = 0.0E0;
  Double I_KINETIC_Dxz_F3z = 0.0E0;
  Double I_KINETIC_D2y_F3z = 0.0E0;
  Double I_KINETIC_Dyz_F3z = 0.0E0;
  Double I_KINETIC_D2z_F3z = 0.0E0;

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
     * totally 9 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_Px_vrr = PAX*I_TWOBODYOVERLAP_Px_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr+oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_Px_vrr = PAY*I_TWOBODYOVERLAP_Py_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Px_vrr+oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_Py_vrr = PAX*I_TWOBODYOVERLAP_Px_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_Py_vrr = PAY*I_TWOBODYOVERLAP_Py_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr+oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Py_vrr+oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_Pz_vrr = PAX*I_TWOBODYOVERLAP_Px_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_Py_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_S_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_S_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_S_F3x_vrr = PBX*I_TWOBODYOVERLAP_S_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_TWOBODYOVERLAP_S_F2xy_vrr = PBY*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_S_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_S_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_S_Fxyz_vrr = PBZ*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_S_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_S_F3y_vrr = PBY*I_TWOBODYOVERLAP_S_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_TWOBODYOVERLAP_S_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_S_Fy2z_vrr = PBY*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_S_F3z_vrr = PBZ*I_TWOBODYOVERLAP_S_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 18 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_D2x_vrr = PAX*I_TWOBODYOVERLAP_Px_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Px_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2x_vrr = PAY*I_TWOBODYOVERLAP_Py_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Px_Py_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_Py_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_Py_Px_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_Px_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_Py_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Px_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2y_vrr = PAX*I_TWOBODYOVERLAP_Px_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2y_vrr = PAY*I_TWOBODYOVERLAP_Py_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Py_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_Px_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_Py_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Py_Pz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Py_vrr;
    Double I_TWOBODYOVERLAP_D2x_D2z_vrr = PAX*I_TWOBODYOVERLAP_Px_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_D2z_vrr = PAY*I_TWOBODYOVERLAP_Py_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_Pz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_F3x_vrr = PAX*I_TWOBODYOVERLAP_S_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Py_F3x_vrr = PAY*I_TWOBODYOVERLAP_S_F3x_vrr;
    Double I_TWOBODYOVERLAP_Pz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_S_F3x_vrr;
    Double I_TWOBODYOVERLAP_Px_F2xy_vrr = PAX*I_TWOBODYOVERLAP_S_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Py_F2xy_vrr = PAY*I_TWOBODYOVERLAP_S_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Pz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_S_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Px_F2xz_vrr = PAX*I_TWOBODYOVERLAP_S_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Py_F2xz_vrr = PAY*I_TWOBODYOVERLAP_S_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Pz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_S_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_TWOBODYOVERLAP_Px_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_S_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Py_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_S_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Pz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_S_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Px_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_S_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Py_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_S_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Pz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_S_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Px_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_S_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Py_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_S_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Pz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_S_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Px_F3y_vrr = PAX*I_TWOBODYOVERLAP_S_F3y_vrr;
    Double I_TWOBODYOVERLAP_Py_F3y_vrr = PAY*I_TWOBODYOVERLAP_S_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Pz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_S_F3y_vrr;
    Double I_TWOBODYOVERLAP_Px_F2yz_vrr = PAX*I_TWOBODYOVERLAP_S_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Py_F2yz_vrr = PAY*I_TWOBODYOVERLAP_S_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Pz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_S_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_TWOBODYOVERLAP_Px_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_S_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Py_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_S_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_TWOBODYOVERLAP_Pz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_S_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Px_F3z_vrr = PAX*I_TWOBODYOVERLAP_S_F3z_vrr;
    Double I_TWOBODYOVERLAP_Py_F3z_vrr = PAY*I_TWOBODYOVERLAP_S_F3z_vrr;
    Double I_TWOBODYOVERLAP_Pz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_S_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_S_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_F3x_vrr = PAX*I_TWOBODYOVERLAP_Px_F3x_vrr+oned2z*I_TWOBODYOVERLAP_S_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F3x_vrr = PAY*I_TWOBODYOVERLAP_Px_F3x_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_Px_F3x_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3x_vrr = PAY*I_TWOBODYOVERLAP_Py_F3x_vrr+oned2z*I_TWOBODYOVERLAP_S_F3x_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_Py_F3x_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_Pz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_S_F3x_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2xy_vrr = PAX*I_TWOBODYOVERLAP_Px_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_S_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F2xy_vrr = PAY*I_TWOBODYOVERLAP_Px_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_Px_F2xy_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_Py_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_S_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_Py_F2xy_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_Pz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_S_F2xy_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2xz_vrr = PAX*I_TWOBODYOVERLAP_Px_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_S_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F2xz_vrr = PAY*I_TWOBODYOVERLAP_Px_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_Px_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_Py_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_S_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_Py_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_S_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_Px_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_S_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_Px_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_Px_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_Py_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_S_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_Py_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_S_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_Px_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_Px_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_Px_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_Py_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_Py_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_S_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_Px_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_S_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_Px_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_Px_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_Py_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_S_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_Py_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_S_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_D2x_F3y_vrr = PAX*I_TWOBODYOVERLAP_Px_F3y_vrr+oned2z*I_TWOBODYOVERLAP_S_F3y_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F3y_vrr = PAY*I_TWOBODYOVERLAP_Px_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_Px_F3y_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3y_vrr = PAY*I_TWOBODYOVERLAP_Py_F3y_vrr+oned2z*I_TWOBODYOVERLAP_S_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_Py_F3y_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_Pz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_S_F3y_vrr;
    Double I_TWOBODYOVERLAP_D2x_F2yz_vrr = PAX*I_TWOBODYOVERLAP_Px_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_S_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F2yz_vrr = PAY*I_TWOBODYOVERLAP_Px_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_Px_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_Py_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_S_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_Py_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_Pz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_S_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_TWOBODYOVERLAP_D2x_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_Px_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_S_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Dxy_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_Px_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dxz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_Px_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_TWOBODYOVERLAP_D2y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_Py_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_S_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_Dyz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_Py_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dyz_vrr;
    Double I_TWOBODYOVERLAP_D2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_Pz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_S_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_D2x_F3z_vrr = PAX*I_TWOBODYOVERLAP_Px_F3z_vrr+oned2z*I_TWOBODYOVERLAP_S_F3z_vrr;
    Double I_TWOBODYOVERLAP_Dxy_F3z_vrr = PAY*I_TWOBODYOVERLAP_Px_F3z_vrr;
    Double I_TWOBODYOVERLAP_Dxz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_Px_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2y_F3z_vrr = PAY*I_TWOBODYOVERLAP_Py_F3z_vrr+oned2z*I_TWOBODYOVERLAP_S_F3z_vrr;
    Double I_TWOBODYOVERLAP_Dyz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_Py_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_D2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_Pz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_S_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Pz_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2x_vrr = PAX*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2x_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2x_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2y_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxy_vrr = PAX*I_TWOBODYOVERLAP_D2z_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxy_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dxy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxy_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Dxy_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dxz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dxz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Px_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dxz_vrr = PAX*I_TWOBODYOVERLAP_D2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dxz_vrr = PAX*I_TWOBODYOVERLAP_D2z_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dxz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Px_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dxz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Dxz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dxz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Px_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2y_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2y_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2y_vrr = PAX*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2y_vrr = PAX*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2y_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2y_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Dyz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Dyz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Py_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Dyz_vrr = PAX*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Dyz_vrr = PAX*I_TWOBODYOVERLAP_D2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F3y_Dyz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Py_vrr;
    Double I_TWOBODYOVERLAP_F3z_Dyz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Dyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Dyz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Py_vrr;
    Double I_TWOBODYOVERLAP_F3x_D2z_vrr = PAX*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_D2z_vrr = PAY*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Pz_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_D2z_vrr = PAX*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_D2z_vrr = PAX*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_D2z_vrr = PAY*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_F3z_D2z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 20 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_F3x_vrr = PAX*I_TWOBODYOVERLAP_D2x_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F3x_vrr = PAY*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F3x_vrr = PAX*I_TWOBODYOVERLAP_D2y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F3x_vrr = PAX*I_TWOBODYOVERLAP_D2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3y_F3x_vrr = PAY*I_TWOBODYOVERLAP_D2y_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F3x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F3x_vrr;
    Double I_TWOBODYOVERLAP_F3x_F2xy_vrr = PAX*I_TWOBODYOVERLAP_D2x_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F2xy_vrr = PAY*I_TWOBODYOVERLAP_D2x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F2xy_vrr = PAX*I_TWOBODYOVERLAP_D2y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_D2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_D2y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_F3z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_F3x_F2xz_vrr = PAX*I_TWOBODYOVERLAP_D2x_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F2xz_vrr = PAY*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2x_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F2xz_vrr = PAX*I_TWOBODYOVERLAP_D2y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_D2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F3y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_D2y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F2xz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_D2x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_D2x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_D2y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_D2z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_D2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Fxyz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_D2x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_D2y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_D2z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_F3x_F3y_vrr = PAX*I_TWOBODYOVERLAP_D2x_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F3y_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F3y_vrr = PAY*I_TWOBODYOVERLAP_D2x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F3y_vrr = PAX*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F3y_vrr = PAX*I_TWOBODYOVERLAP_D2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3y_F3y_vrr = PAY*I_TWOBODYOVERLAP_D2y_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F3y_vrr;
    Double I_TWOBODYOVERLAP_F3x_F2yz_vrr = PAX*I_TWOBODYOVERLAP_D2x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F2yz_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F2yz_vrr = PAY*I_TWOBODYOVERLAP_D2x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2y_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F2yz_vrr = PAX*I_TWOBODYOVERLAP_D2y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_D2z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_F3y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_D2y_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_F3x_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_D2y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_D2z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_F3y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F3z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_D2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_F3x_F3z_vrr = PAX*I_TWOBODYOVERLAP_D2x_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2xy_F3z_vrr = PAY*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2xz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_D2x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_D2z_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_F3z_vrr = PAX*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_F3z_vrr = PAX*I_TWOBODYOVERLAP_D2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_F3y_F3z_vrr = PAY*I_TWOBODYOVERLAP_D2y_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_F3z_vrr;
    Double I_TWOBODYOVERLAP_F2yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_D2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_F3z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_D2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_F3x_vrr = PAX*I_TWOBODYOVERLAP_F3x_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3x_vrr = PAY*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3x_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3x_vrr = PAX*I_TWOBODYOVERLAP_F3y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F3x_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3x_vrr = PAX*I_TWOBODYOVERLAP_F3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3x_vrr = PAY*I_TWOBODYOVERLAP_F3y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3x_vrr = PAY*I_TWOBODYOVERLAP_F3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3x_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2xy_vrr = PAX*I_TWOBODYOVERLAP_F3x_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2xy_vrr = PAY*I_TWOBODYOVERLAP_F3x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2xy_vrr = PAX*I_TWOBODYOVERLAP_F3y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2xy_vrr = PAX*I_TWOBODYOVERLAP_F3z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2xy_vrr = PAY*I_TWOBODYOVERLAP_F3y_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2xy_vrr = PAY*I_TWOBODYOVERLAP_F3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2xy_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2xz_vrr = PAX*I_TWOBODYOVERLAP_F3x_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2xz_vrr = PAY*I_TWOBODYOVERLAP_F3x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2xz_vrr = PAX*I_TWOBODYOVERLAP_F3y_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2xz_vrr = PAX*I_TWOBODYOVERLAP_F3z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2xz_vrr = PAY*I_TWOBODYOVERLAP_F3y_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2xz_vrr = PAY*I_TWOBODYOVERLAP_F3z_F2xz_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2xz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_F3x_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_F3x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_F3y_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr = PAX*I_TWOBODYOVERLAP_F3z_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_F3y_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr = PAY*I_TWOBODYOVERLAP_F3z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fx2y_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_F3x_Fxyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_F3x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_F3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr = PAX*I_TWOBODYOVERLAP_F3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_F3y_Fxyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr = PAY*I_TWOBODYOVERLAP_F3z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fxyz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Fxyz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Fxyz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Dxy_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_F3x_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_F3x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_F3y_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr = PAX*I_TWOBODYOVERLAP_F3z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_F3y_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Dxz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr = PAY*I_TWOBODYOVERLAP_F3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fx2z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dxz_vrr;
    Double I_TWOBODYOVERLAP_G4x_F3y_vrr = PAX*I_TWOBODYOVERLAP_F3x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3y_vrr = PAY*I_TWOBODYOVERLAP_F3x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F3y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3y_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3y_vrr = PAX*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F3y_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3y_vrr = PAX*I_TWOBODYOVERLAP_F3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3y_vrr = PAY*I_TWOBODYOVERLAP_F3y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3y_vrr = PAY*I_TWOBODYOVERLAP_F3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3y_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_G4x_F2yz_vrr = PAX*I_TWOBODYOVERLAP_F3x_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F2yz_vrr = PAY*I_TWOBODYOVERLAP_F3x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F2yz_vrr = PAX*I_TWOBODYOVERLAP_F3y_F2yz_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F2yz_vrr = PAX*I_TWOBODYOVERLAP_F3z_F2yz_vrr;
    Double I_TWOBODYOVERLAP_G4y_F2yz_vrr = PAY*I_TWOBODYOVERLAP_F3y_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F2yz_vrr = PAY*I_TWOBODYOVERLAP_F3z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G4z_F2yz_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_G4x_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_F3x_Fy2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_F3x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3x_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_F3y_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr = PAX*I_TWOBODYOVERLAP_F3z_Fy2z_vrr;
    Double I_TWOBODYOVERLAP_G4y_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_F3y_Fy2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3y_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_Dyz_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr = PAY*I_TWOBODYOVERLAP_F3z_Fy2z_vrr+oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_G4z_Fy2z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_Fy2z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_Fy2z_vrr+2*oned2z*I_TWOBODYOVERLAP_F3z_Dyz_vrr;
    Double I_TWOBODYOVERLAP_G4x_F3z_vrr = PAX*I_TWOBODYOVERLAP_F3x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_TWOBODYOVERLAP_G3xy_F3z_vrr = PAY*I_TWOBODYOVERLAP_F3x_F3z_vrr;
    Double I_TWOBODYOVERLAP_G3xz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F3x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3x_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_F3z_vrr = PAY*I_TWOBODYOVERLAP_F2xy_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_F3z_vrr = PAX*I_TWOBODYOVERLAP_F3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_F3z_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_F3z_vrr = PAX*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_G4y_F3z_vrr = PAY*I_TWOBODYOVERLAP_F3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_G3yz_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_D2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_F3z_vrr = PAY*I_TWOBODYOVERLAP_F3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_G4z_F3z_vrr = PAZ*I_TWOBODYOVERLAP_F3z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_F3z_D2z_vrr;

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
     * totally 9 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_KINETIC_D2x_Px_vrr = PAX*I_KINETIC_Px_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_D2y_Px_vrr = PAY*I_KINETIC_Py_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_D2z_Px_vrr = PAZ*I_KINETIC_Pz_Px_vrr+oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Px_vrr-bdz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_D2x_Py_vrr = PAX*I_KINETIC_Px_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_D2y_Py_vrr = PAY*I_KINETIC_Py_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_D2z_Py_vrr = PAZ*I_KINETIC_Pz_Py_vrr+oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Py_vrr-bdz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_D2x_Pz_vrr = PAX*I_KINETIC_Px_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_KINETIC_D2y_Pz_vrr = PAY*I_KINETIC_Py_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;
    Double I_KINETIC_D2z_Pz_vrr = PAZ*I_KINETIC_Pz_Pz_vrr+oned2z*I_KINETIC_S_Pz_vrr+oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Pz_vrr-bdz*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_S_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_D
     * RHS shell quartet name: SQ_KINETIC_S_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_P
     ************************************************************/
    Double I_KINETIC_S_F3x_vrr = PBX*I_KINETIC_S_D2x_vrr+2*oned2z*I_KINETIC_S_Px_vrr+twoxi*I_TWOBODYOVERLAP_S_F3x_vrr-2*adz*I_TWOBODYOVERLAP_S_Px_vrr;
    Double I_KINETIC_S_F2xy_vrr = PBY*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_S_F2xy_vrr;
    Double I_KINETIC_S_F2xz_vrr = PBZ*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_S_F2xz_vrr;
    Double I_KINETIC_S_Fx2y_vrr = PBX*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_S_Fx2y_vrr;
    Double I_KINETIC_S_Fxyz_vrr = PBZ*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_S_Fxyz_vrr;
    Double I_KINETIC_S_Fx2z_vrr = PBX*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_S_Fx2z_vrr;
    Double I_KINETIC_S_F3y_vrr = PBY*I_KINETIC_S_D2y_vrr+2*oned2z*I_KINETIC_S_Py_vrr+twoxi*I_TWOBODYOVERLAP_S_F3y_vrr-2*adz*I_TWOBODYOVERLAP_S_Py_vrr;
    Double I_KINETIC_S_F2yz_vrr = PBZ*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_S_F2yz_vrr;
    Double I_KINETIC_S_Fy2z_vrr = PBY*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_S_Fy2z_vrr;
    Double I_KINETIC_S_F3z_vrr = PBZ*I_KINETIC_S_D2z_vrr+2*oned2z*I_KINETIC_S_Pz_vrr+twoxi*I_TWOBODYOVERLAP_S_F3z_vrr-2*adz*I_TWOBODYOVERLAP_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 18 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_KINETIC_S_D
     * RHS shell quartet name: SQ_KINETIC_P_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_D
     ************************************************************/
    Double I_KINETIC_D2x_D2x_vrr = PAX*I_KINETIC_Px_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+2*oned2z*I_KINETIC_Px_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_D2y_D2x_vrr = PAY*I_KINETIC_Py_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_D2z_D2x_vrr = PAZ*I_KINETIC_Pz_D2x_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2x_vrr-bdz*I_TWOBODYOVERLAP_S_D2x_vrr;
    Double I_KINETIC_D2x_Dxy_vrr = PAX*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_Px_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2y_Dxy_vrr = PAY*I_KINETIC_Py_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+oned2z*I_KINETIC_Py_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2z_Dxy_vrr = PAZ*I_KINETIC_Pz_Dxy_vrr+oned2z*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxy_vrr-bdz*I_TWOBODYOVERLAP_S_Dxy_vrr;
    Double I_KINETIC_D2x_Dxz_vrr = PAX*I_KINETIC_Px_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_Px_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_D2y_Dxz_vrr = PAY*I_KINETIC_Py_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_D2z_Dxz_vrr = PAZ*I_KINETIC_Pz_Dxz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+oned2z*I_KINETIC_Pz_Px_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dxz_vrr-bdz*I_TWOBODYOVERLAP_S_Dxz_vrr;
    Double I_KINETIC_D2x_D2y_vrr = PAX*I_KINETIC_Px_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_D2y_D2y_vrr = PAY*I_KINETIC_Py_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+2*oned2z*I_KINETIC_Py_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_D2z_D2y_vrr = PAZ*I_KINETIC_Pz_D2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2y_vrr-bdz*I_TWOBODYOVERLAP_S_D2y_vrr;
    Double I_KINETIC_D2x_Dyz_vrr = PAX*I_KINETIC_Px_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2y_Dyz_vrr = PAY*I_KINETIC_Py_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_Py_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2z_Dyz_vrr = PAZ*I_KINETIC_Pz_Dyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+oned2z*I_KINETIC_Pz_Py_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Dyz_vrr-bdz*I_TWOBODYOVERLAP_S_Dyz_vrr;
    Double I_KINETIC_D2x_D2z_vrr = PAX*I_KINETIC_Px_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_KINETIC_D2y_D2z_vrr = PAY*I_KINETIC_Py_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;
    Double I_KINETIC_D2z_D2z_vrr = PAZ*I_KINETIC_Pz_D2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+2*oned2z*I_KINETIC_Pz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_D2z_D2z_vrr-bdz*I_TWOBODYOVERLAP_S_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_F
     * RHS shell quartet name: SQ_KINETIC_S_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_F
     ************************************************************/
    Double I_KINETIC_Px_F3x_vrr = PAX*I_KINETIC_S_F3x_vrr+3*oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Px_F3x_vrr;
    Double I_KINETIC_Py_F3x_vrr = PAY*I_KINETIC_S_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Py_F3x_vrr;
    Double I_KINETIC_Pz_F3x_vrr = PAZ*I_KINETIC_S_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F3x_vrr;
    Double I_KINETIC_Px_F2xy_vrr = PAX*I_KINETIC_S_F2xy_vrr+2*oned2z*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Px_F2xy_vrr;
    Double I_KINETIC_Py_F2xy_vrr = PAY*I_KINETIC_S_F2xy_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Py_F2xy_vrr;
    Double I_KINETIC_Pz_F2xy_vrr = PAZ*I_KINETIC_S_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F2xy_vrr;
    Double I_KINETIC_Px_F2xz_vrr = PAX*I_KINETIC_S_F2xz_vrr+2*oned2z*I_KINETIC_S_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Px_F2xz_vrr;
    Double I_KINETIC_Py_F2xz_vrr = PAY*I_KINETIC_S_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Py_F2xz_vrr;
    Double I_KINETIC_Pz_F2xz_vrr = PAZ*I_KINETIC_S_F2xz_vrr+oned2z*I_KINETIC_S_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F2xz_vrr;
    Double I_KINETIC_Px_Fx2y_vrr = PAX*I_KINETIC_S_Fx2y_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Px_Fx2y_vrr;
    Double I_KINETIC_Py_Fx2y_vrr = PAY*I_KINETIC_S_Fx2y_vrr+2*oned2z*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Py_Fx2y_vrr;
    Double I_KINETIC_Pz_Fx2y_vrr = PAZ*I_KINETIC_S_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Fx2y_vrr;
    Double I_KINETIC_Px_Fxyz_vrr = PAX*I_KINETIC_S_Fxyz_vrr+oned2z*I_KINETIC_S_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Px_Fxyz_vrr;
    Double I_KINETIC_Py_Fxyz_vrr = PAY*I_KINETIC_S_Fxyz_vrr+oned2z*I_KINETIC_S_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Py_Fxyz_vrr;
    Double I_KINETIC_Pz_Fxyz_vrr = PAZ*I_KINETIC_S_Fxyz_vrr+oned2z*I_KINETIC_S_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Fxyz_vrr;
    Double I_KINETIC_Px_Fx2z_vrr = PAX*I_KINETIC_S_Fx2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Px_Fx2z_vrr;
    Double I_KINETIC_Py_Fx2z_vrr = PAY*I_KINETIC_S_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Py_Fx2z_vrr;
    Double I_KINETIC_Pz_Fx2z_vrr = PAZ*I_KINETIC_S_Fx2z_vrr+2*oned2z*I_KINETIC_S_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Fx2z_vrr;
    Double I_KINETIC_Px_F3y_vrr = PAX*I_KINETIC_S_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Px_F3y_vrr;
    Double I_KINETIC_Py_F3y_vrr = PAY*I_KINETIC_S_F3y_vrr+3*oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Py_F3y_vrr;
    Double I_KINETIC_Pz_F3y_vrr = PAZ*I_KINETIC_S_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F3y_vrr;
    Double I_KINETIC_Px_F2yz_vrr = PAX*I_KINETIC_S_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Px_F2yz_vrr;
    Double I_KINETIC_Py_F2yz_vrr = PAY*I_KINETIC_S_F2yz_vrr+2*oned2z*I_KINETIC_S_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Py_F2yz_vrr;
    Double I_KINETIC_Pz_F2yz_vrr = PAZ*I_KINETIC_S_F2yz_vrr+oned2z*I_KINETIC_S_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F2yz_vrr;
    Double I_KINETIC_Px_Fy2z_vrr = PAX*I_KINETIC_S_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Px_Fy2z_vrr;
    Double I_KINETIC_Py_Fy2z_vrr = PAY*I_KINETIC_S_Fy2z_vrr+oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Py_Fy2z_vrr;
    Double I_KINETIC_Pz_Fy2z_vrr = PAZ*I_KINETIC_S_Fy2z_vrr+2*oned2z*I_KINETIC_S_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Pz_Fy2z_vrr;
    Double I_KINETIC_Px_F3z_vrr = PAX*I_KINETIC_S_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Px_F3z_vrr;
    Double I_KINETIC_Py_F3z_vrr = PAY*I_KINETIC_S_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Py_F3z_vrr;
    Double I_KINETIC_Pz_F3z_vrr = PAZ*I_KINETIC_S_F3z_vrr+3*oned2z*I_KINETIC_S_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Pz_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_F
     * RHS shell quartet name: SQ_KINETIC_S_F
     * RHS shell quartet name: SQ_KINETIC_P_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_F
     ************************************************************/
    Double I_KINETIC_D2x_F3x_vrr = PAX*I_KINETIC_Px_F3x_vrr+oned2z*I_KINETIC_S_F3x_vrr+3*oned2z*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3x_vrr-bdz*I_TWOBODYOVERLAP_S_F3x_vrr;
    Double I_KINETIC_Dxy_F3x_vrr = PAY*I_KINETIC_Px_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F3x_vrr;
    Double I_KINETIC_Dxz_F3x_vrr = PAZ*I_KINETIC_Px_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F3x_vrr;
    Double I_KINETIC_D2y_F3x_vrr = PAY*I_KINETIC_Py_F3x_vrr+oned2z*I_KINETIC_S_F3x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3x_vrr-bdz*I_TWOBODYOVERLAP_S_F3x_vrr;
    Double I_KINETIC_Dyz_F3x_vrr = PAZ*I_KINETIC_Py_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F3x_vrr;
    Double I_KINETIC_D2z_F3x_vrr = PAZ*I_KINETIC_Pz_F3x_vrr+oned2z*I_KINETIC_S_F3x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_S_F3x_vrr;
    Double I_KINETIC_D2x_F2xy_vrr = PAX*I_KINETIC_Px_F2xy_vrr+oned2z*I_KINETIC_S_F2xy_vrr+2*oned2z*I_KINETIC_Px_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2xy_vrr-bdz*I_TWOBODYOVERLAP_S_F2xy_vrr;
    Double I_KINETIC_Dxy_F2xy_vrr = PAY*I_KINETIC_Px_F2xy_vrr+oned2z*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F2xy_vrr;
    Double I_KINETIC_Dxz_F2xy_vrr = PAZ*I_KINETIC_Px_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F2xy_vrr;
    Double I_KINETIC_D2y_F2xy_vrr = PAY*I_KINETIC_Py_F2xy_vrr+oned2z*I_KINETIC_S_F2xy_vrr+oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2xy_vrr-bdz*I_TWOBODYOVERLAP_S_F2xy_vrr;
    Double I_KINETIC_Dyz_F2xy_vrr = PAZ*I_KINETIC_Py_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F2xy_vrr;
    Double I_KINETIC_D2z_F2xy_vrr = PAZ*I_KINETIC_Pz_F2xy_vrr+oned2z*I_KINETIC_S_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_S_F2xy_vrr;
    Double I_KINETIC_D2x_F2xz_vrr = PAX*I_KINETIC_Px_F2xz_vrr+oned2z*I_KINETIC_S_F2xz_vrr+2*oned2z*I_KINETIC_Px_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2xz_vrr-bdz*I_TWOBODYOVERLAP_S_F2xz_vrr;
    Double I_KINETIC_Dxy_F2xz_vrr = PAY*I_KINETIC_Px_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F2xz_vrr;
    Double I_KINETIC_Dxz_F2xz_vrr = PAZ*I_KINETIC_Px_F2xz_vrr+oned2z*I_KINETIC_Px_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F2xz_vrr;
    Double I_KINETIC_D2y_F2xz_vrr = PAY*I_KINETIC_Py_F2xz_vrr+oned2z*I_KINETIC_S_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2xz_vrr-bdz*I_TWOBODYOVERLAP_S_F2xz_vrr;
    Double I_KINETIC_Dyz_F2xz_vrr = PAZ*I_KINETIC_Py_F2xz_vrr+oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F2xz_vrr;
    Double I_KINETIC_D2z_F2xz_vrr = PAZ*I_KINETIC_Pz_F2xz_vrr+oned2z*I_KINETIC_S_F2xz_vrr+oned2z*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_S_F2xz_vrr;
    Double I_KINETIC_D2x_Fx2y_vrr = PAX*I_KINETIC_Px_Fx2y_vrr+oned2z*I_KINETIC_S_Fx2y_vrr+oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_S_Fx2y_vrr;
    Double I_KINETIC_Dxy_Fx2y_vrr = PAY*I_KINETIC_Px_Fx2y_vrr+2*oned2z*I_KINETIC_Px_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fx2y_vrr;
    Double I_KINETIC_Dxz_Fx2y_vrr = PAZ*I_KINETIC_Px_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Fx2y_vrr;
    Double I_KINETIC_D2y_Fx2y_vrr = PAY*I_KINETIC_Py_Fx2y_vrr+oned2z*I_KINETIC_S_Fx2y_vrr+2*oned2z*I_KINETIC_Py_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_S_Fx2y_vrr;
    Double I_KINETIC_Dyz_Fx2y_vrr = PAZ*I_KINETIC_Py_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Fx2y_vrr;
    Double I_KINETIC_D2z_Fx2y_vrr = PAZ*I_KINETIC_Pz_Fx2y_vrr+oned2z*I_KINETIC_S_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_S_Fx2y_vrr;
    Double I_KINETIC_D2x_Fxyz_vrr = PAX*I_KINETIC_Px_Fxyz_vrr+oned2z*I_KINETIC_S_Fxyz_vrr+oned2z*I_KINETIC_Px_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_S_Fxyz_vrr;
    Double I_KINETIC_Dxy_Fxyz_vrr = PAY*I_KINETIC_Px_Fxyz_vrr+oned2z*I_KINETIC_Px_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fxyz_vrr;
    Double I_KINETIC_Dxz_Fxyz_vrr = PAZ*I_KINETIC_Px_Fxyz_vrr+oned2z*I_KINETIC_Px_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Fxyz_vrr;
    Double I_KINETIC_D2y_Fxyz_vrr = PAY*I_KINETIC_Py_Fxyz_vrr+oned2z*I_KINETIC_S_Fxyz_vrr+oned2z*I_KINETIC_Py_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_S_Fxyz_vrr;
    Double I_KINETIC_Dyz_Fxyz_vrr = PAZ*I_KINETIC_Py_Fxyz_vrr+oned2z*I_KINETIC_Py_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Fxyz_vrr;
    Double I_KINETIC_D2z_Fxyz_vrr = PAZ*I_KINETIC_Pz_Fxyz_vrr+oned2z*I_KINETIC_S_Fxyz_vrr+oned2z*I_KINETIC_Pz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_S_Fxyz_vrr;
    Double I_KINETIC_D2x_Fx2z_vrr = PAX*I_KINETIC_Px_Fx2z_vrr+oned2z*I_KINETIC_S_Fx2z_vrr+oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_S_Fx2z_vrr;
    Double I_KINETIC_Dxy_Fx2z_vrr = PAY*I_KINETIC_Px_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fx2z_vrr;
    Double I_KINETIC_Dxz_Fx2z_vrr = PAZ*I_KINETIC_Px_Fx2z_vrr+2*oned2z*I_KINETIC_Px_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Fx2z_vrr;
    Double I_KINETIC_D2y_Fx2z_vrr = PAY*I_KINETIC_Py_Fx2z_vrr+oned2z*I_KINETIC_S_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_S_Fx2z_vrr;
    Double I_KINETIC_Dyz_Fx2z_vrr = PAZ*I_KINETIC_Py_Fx2z_vrr+2*oned2z*I_KINETIC_Py_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Fx2z_vrr;
    Double I_KINETIC_D2z_Fx2z_vrr = PAZ*I_KINETIC_Pz_Fx2z_vrr+oned2z*I_KINETIC_S_Fx2z_vrr+2*oned2z*I_KINETIC_Pz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_S_Fx2z_vrr;
    Double I_KINETIC_D2x_F3y_vrr = PAX*I_KINETIC_Px_F3y_vrr+oned2z*I_KINETIC_S_F3y_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3y_vrr-bdz*I_TWOBODYOVERLAP_S_F3y_vrr;
    Double I_KINETIC_Dxy_F3y_vrr = PAY*I_KINETIC_Px_F3y_vrr+3*oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F3y_vrr;
    Double I_KINETIC_Dxz_F3y_vrr = PAZ*I_KINETIC_Px_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F3y_vrr;
    Double I_KINETIC_D2y_F3y_vrr = PAY*I_KINETIC_Py_F3y_vrr+oned2z*I_KINETIC_S_F3y_vrr+3*oned2z*I_KINETIC_Py_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3y_vrr-bdz*I_TWOBODYOVERLAP_S_F3y_vrr;
    Double I_KINETIC_Dyz_F3y_vrr = PAZ*I_KINETIC_Py_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F3y_vrr;
    Double I_KINETIC_D2z_F3y_vrr = PAZ*I_KINETIC_Pz_F3y_vrr+oned2z*I_KINETIC_S_F3y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_S_F3y_vrr;
    Double I_KINETIC_D2x_F2yz_vrr = PAX*I_KINETIC_Px_F2yz_vrr+oned2z*I_KINETIC_S_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F2yz_vrr-bdz*I_TWOBODYOVERLAP_S_F2yz_vrr;
    Double I_KINETIC_Dxy_F2yz_vrr = PAY*I_KINETIC_Px_F2yz_vrr+2*oned2z*I_KINETIC_Px_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F2yz_vrr;
    Double I_KINETIC_Dxz_F2yz_vrr = PAZ*I_KINETIC_Px_F2yz_vrr+oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F2yz_vrr;
    Double I_KINETIC_D2y_F2yz_vrr = PAY*I_KINETIC_Py_F2yz_vrr+oned2z*I_KINETIC_S_F2yz_vrr+2*oned2z*I_KINETIC_Py_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F2yz_vrr-bdz*I_TWOBODYOVERLAP_S_F2yz_vrr;
    Double I_KINETIC_Dyz_F2yz_vrr = PAZ*I_KINETIC_Py_F2yz_vrr+oned2z*I_KINETIC_Py_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F2yz_vrr;
    Double I_KINETIC_D2z_F2yz_vrr = PAZ*I_KINETIC_Pz_F2yz_vrr+oned2z*I_KINETIC_S_F2yz_vrr+oned2z*I_KINETIC_Pz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_S_F2yz_vrr;
    Double I_KINETIC_D2x_Fy2z_vrr = PAX*I_KINETIC_Px_Fy2z_vrr+oned2z*I_KINETIC_S_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_S_Fy2z_vrr;
    Double I_KINETIC_Dxy_Fy2z_vrr = PAY*I_KINETIC_Px_Fy2z_vrr+oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_Fy2z_vrr;
    Double I_KINETIC_Dxz_Fy2z_vrr = PAZ*I_KINETIC_Px_Fy2z_vrr+2*oned2z*I_KINETIC_Px_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_Fy2z_vrr;
    Double I_KINETIC_D2y_Fy2z_vrr = PAY*I_KINETIC_Py_Fy2z_vrr+oned2z*I_KINETIC_S_Fy2z_vrr+oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_S_Fy2z_vrr;
    Double I_KINETIC_Dyz_Fy2z_vrr = PAZ*I_KINETIC_Py_Fy2z_vrr+2*oned2z*I_KINETIC_Py_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_Fy2z_vrr;
    Double I_KINETIC_D2z_Fy2z_vrr = PAZ*I_KINETIC_Pz_Fy2z_vrr+oned2z*I_KINETIC_S_Fy2z_vrr+2*oned2z*I_KINETIC_Pz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_D2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_S_Fy2z_vrr;
    Double I_KINETIC_D2x_F3z_vrr = PAX*I_KINETIC_Px_F3z_vrr+oned2z*I_KINETIC_S_F3z_vrr+twoxi*I_TWOBODYOVERLAP_D2x_F3z_vrr-bdz*I_TWOBODYOVERLAP_S_F3z_vrr;
    Double I_KINETIC_Dxy_F3z_vrr = PAY*I_KINETIC_Px_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_F3z_vrr;
    Double I_KINETIC_Dxz_F3z_vrr = PAZ*I_KINETIC_Px_F3z_vrr+3*oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dxz_F3z_vrr;
    Double I_KINETIC_D2y_F3z_vrr = PAY*I_KINETIC_Py_F3z_vrr+oned2z*I_KINETIC_S_F3z_vrr+twoxi*I_TWOBODYOVERLAP_D2y_F3z_vrr-bdz*I_TWOBODYOVERLAP_S_F3z_vrr;
    Double I_KINETIC_Dyz_F3z_vrr = PAZ*I_KINETIC_Py_F3z_vrr+3*oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Dyz_F3z_vrr;
    Double I_KINETIC_D2z_F3z_vrr = PAZ*I_KINETIC_Pz_F3z_vrr+oned2z*I_KINETIC_S_F3z_vrr+3*oned2z*I_KINETIC_Pz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_D2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_S_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 12 integrals are omitted 
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
    Double I_KINETIC_Fx2z_D2x_vrr = PAX*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2x_vrr;
    Double I_KINETIC_F3y_D2x_vrr = PAY*I_KINETIC_D2y_D2x_vrr+2*oned2z*I_KINETIC_Py_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2x_vrr;
    Double I_KINETIC_F2yz_D2x_vrr = PAZ*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2x_vrr;
    Double I_KINETIC_F3z_D2x_vrr = PAZ*I_KINETIC_D2z_D2x_vrr+2*oned2z*I_KINETIC_Pz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2x_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2x_vrr;
    Double I_KINETIC_F3x_Dxy_vrr = PAX*I_KINETIC_D2x_Dxy_vrr+2*oned2z*I_KINETIC_Px_Dxy_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dxy_vrr;
    Double I_KINETIC_F2xy_Dxy_vrr = PAY*I_KINETIC_D2x_Dxy_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxy_vrr;
    Double I_KINETIC_F2xz_Dxy_vrr = PAZ*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxy_vrr;
    Double I_KINETIC_Fx2y_Dxy_vrr = PAX*I_KINETIC_D2y_Dxy_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxy_vrr;
    Double I_KINETIC_Fx2z_Dxy_vrr = PAX*I_KINETIC_D2z_Dxy_vrr+oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxy_vrr;
    Double I_KINETIC_F3y_Dxy_vrr = PAY*I_KINETIC_D2y_Dxy_vrr+2*oned2z*I_KINETIC_Py_Dxy_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dxy_vrr;
    Double I_KINETIC_F2yz_Dxy_vrr = PAZ*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxy_vrr;
    Double I_KINETIC_F3z_Dxy_vrr = PAZ*I_KINETIC_D2z_Dxy_vrr+2*oned2z*I_KINETIC_Pz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxy_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Dxy_vrr;
    Double I_KINETIC_F3x_Dxz_vrr = PAX*I_KINETIC_D2x_Dxz_vrr+2*oned2z*I_KINETIC_Px_Dxz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dxz_vrr;
    Double I_KINETIC_F2xy_Dxz_vrr = PAY*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dxz_vrr;
    Double I_KINETIC_F2xz_Dxz_vrr = PAZ*I_KINETIC_D2x_Dxz_vrr+oned2z*I_KINETIC_D2x_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dxz_vrr;
    Double I_KINETIC_Fx2y_Dxz_vrr = PAX*I_KINETIC_D2y_Dxz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dxz_vrr;
    Double I_KINETIC_Fx2z_Dxz_vrr = PAX*I_KINETIC_D2z_Dxz_vrr+oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dxz_vrr;
    Double I_KINETIC_F3y_Dxz_vrr = PAY*I_KINETIC_D2y_Dxz_vrr+2*oned2z*I_KINETIC_Py_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dxz_vrr;
    Double I_KINETIC_F2yz_Dxz_vrr = PAZ*I_KINETIC_D2y_Dxz_vrr+oned2z*I_KINETIC_D2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dxz_vrr;
    Double I_KINETIC_F3z_Dxz_vrr = PAZ*I_KINETIC_D2z_Dxz_vrr+2*oned2z*I_KINETIC_Pz_Dxz_vrr+oned2z*I_KINETIC_D2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dxz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Dxz_vrr;
    Double I_KINETIC_F3x_D2y_vrr = PAX*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_Px_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Px_D2y_vrr;
    Double I_KINETIC_F2xy_D2y_vrr = PAY*I_KINETIC_D2x_D2y_vrr+2*oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2y_vrr;
    Double I_KINETIC_F2xz_D2y_vrr = PAZ*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2y_vrr;
    Double I_KINETIC_Fx2y_D2y_vrr = PAX*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2y_vrr;
    Double I_KINETIC_Fx2z_D2y_vrr = PAX*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2y_vrr;
    Double I_KINETIC_F3y_D2y_vrr = PAY*I_KINETIC_D2y_D2y_vrr+2*oned2z*I_KINETIC_Py_D2y_vrr+2*oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2y_vrr;
    Double I_KINETIC_F2yz_D2y_vrr = PAZ*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2y_vrr;
    Double I_KINETIC_F3z_D2y_vrr = PAZ*I_KINETIC_D2z_D2y_vrr+2*oned2z*I_KINETIC_Pz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2y_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2y_vrr;
    Double I_KINETIC_F3x_Dyz_vrr = PAX*I_KINETIC_D2x_Dyz_vrr+2*oned2z*I_KINETIC_Px_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Dyz_vrr;
    Double I_KINETIC_F2xy_Dyz_vrr = PAY*I_KINETIC_D2x_Dyz_vrr+oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Dyz_vrr;
    Double I_KINETIC_F2xz_Dyz_vrr = PAZ*I_KINETIC_D2x_Dyz_vrr+oned2z*I_KINETIC_D2x_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Dyz_vrr;
    Double I_KINETIC_Fx2y_Dyz_vrr = PAX*I_KINETIC_D2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Dyz_vrr;
    Double I_KINETIC_Fx2z_Dyz_vrr = PAX*I_KINETIC_D2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Dyz_vrr;
    Double I_KINETIC_F3y_Dyz_vrr = PAY*I_KINETIC_D2y_Dyz_vrr+2*oned2z*I_KINETIC_Py_Dyz_vrr+oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Dyz_vrr;
    Double I_KINETIC_F2yz_Dyz_vrr = PAZ*I_KINETIC_D2y_Dyz_vrr+oned2z*I_KINETIC_D2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Dyz_vrr;
    Double I_KINETIC_F3z_Dyz_vrr = PAZ*I_KINETIC_D2z_Dyz_vrr+2*oned2z*I_KINETIC_Pz_Dyz_vrr+oned2z*I_KINETIC_D2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Dyz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Dyz_vrr;
    Double I_KINETIC_F3x_D2z_vrr = PAX*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_Px_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Px_D2z_vrr;
    Double I_KINETIC_F2xy_D2z_vrr = PAY*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_D2z_vrr;
    Double I_KINETIC_F2xz_D2z_vrr = PAZ*I_KINETIC_D2x_D2z_vrr+2*oned2z*I_KINETIC_D2x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_D2z_vrr;
    Double I_KINETIC_Fx2y_D2z_vrr = PAX*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_D2z_vrr;
    Double I_KINETIC_Fx2z_D2z_vrr = PAX*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_D2z_vrr;
    Double I_KINETIC_F3y_D2z_vrr = PAY*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_Py_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_D2z_vrr;
    Double I_KINETIC_F2yz_D2z_vrr = PAZ*I_KINETIC_D2y_D2z_vrr+2*oned2z*I_KINETIC_D2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_D2z_vrr;
    Double I_KINETIC_F3z_D2z_vrr = PAZ*I_KINETIC_D2z_D2z_vrr+2*oned2z*I_KINETIC_Pz_D2z_vrr+2*oned2z*I_KINETIC_D2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_F3z_D2z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 20 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_F
     * RHS shell quartet name: SQ_KINETIC_P_F
     * RHS shell quartet name: SQ_KINETIC_D_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_F
     ************************************************************/
    Double I_KINETIC_F3x_F3x_vrr = PAX*I_KINETIC_D2x_F3x_vrr+2*oned2z*I_KINETIC_Px_F3x_vrr+3*oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F3x_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F3x_vrr;
    Double I_KINETIC_F2xy_F3x_vrr = PAY*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F3x_vrr;
    Double I_KINETIC_F2xz_F3x_vrr = PAZ*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F3x_vrr;
    Double I_KINETIC_Fx2y_F3x_vrr = PAX*I_KINETIC_D2y_F3x_vrr+3*oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F3x_vrr;
    Double I_KINETIC_Fx2z_F3x_vrr = PAX*I_KINETIC_D2z_F3x_vrr+3*oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F3x_vrr;
    Double I_KINETIC_F3y_F3x_vrr = PAY*I_KINETIC_D2y_F3x_vrr+2*oned2z*I_KINETIC_Py_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F3x_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F3x_vrr;
    Double I_KINETIC_F2yz_F3x_vrr = PAZ*I_KINETIC_D2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F3x_vrr;
    Double I_KINETIC_F3z_F3x_vrr = PAZ*I_KINETIC_D2z_F3x_vrr+2*oned2z*I_KINETIC_Pz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F3x_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F3x_vrr;
    Double I_KINETIC_F3x_F2xy_vrr = PAX*I_KINETIC_D2x_F2xy_vrr+2*oned2z*I_KINETIC_Px_F2xy_vrr+2*oned2z*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F2xy_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F2xy_vrr;
    Double I_KINETIC_F2xy_F2xy_vrr = PAY*I_KINETIC_D2x_F2xy_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F2xy_vrr;
    Double I_KINETIC_F2xz_F2xy_vrr = PAZ*I_KINETIC_D2x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F2xy_vrr;
    Double I_KINETIC_Fx2y_F2xy_vrr = PAX*I_KINETIC_D2y_F2xy_vrr+2*oned2z*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F2xy_vrr;
    Double I_KINETIC_Fx2z_F2xy_vrr = PAX*I_KINETIC_D2z_F2xy_vrr+2*oned2z*I_KINETIC_D2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F2xy_vrr;
    Double I_KINETIC_F3y_F2xy_vrr = PAY*I_KINETIC_D2y_F2xy_vrr+2*oned2z*I_KINETIC_Py_F2xy_vrr+oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F2xy_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F2xy_vrr;
    Double I_KINETIC_F2yz_F2xy_vrr = PAZ*I_KINETIC_D2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F2xy_vrr;
    Double I_KINETIC_F3z_F2xy_vrr = PAZ*I_KINETIC_D2z_F2xy_vrr+2*oned2z*I_KINETIC_Pz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F2xy_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F2xy_vrr;
    Double I_KINETIC_F3x_F2xz_vrr = PAX*I_KINETIC_D2x_F2xz_vrr+2*oned2z*I_KINETIC_Px_F2xz_vrr+2*oned2z*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F2xz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F2xz_vrr;
    Double I_KINETIC_F2xy_F2xz_vrr = PAY*I_KINETIC_D2x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F2xz_vrr;
    Double I_KINETIC_F2xz_F2xz_vrr = PAZ*I_KINETIC_D2x_F2xz_vrr+oned2z*I_KINETIC_D2x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F2xz_vrr;
    Double I_KINETIC_Fx2y_F2xz_vrr = PAX*I_KINETIC_D2y_F2xz_vrr+2*oned2z*I_KINETIC_D2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F2xz_vrr;
    Double I_KINETIC_Fx2z_F2xz_vrr = PAX*I_KINETIC_D2z_F2xz_vrr+2*oned2z*I_KINETIC_D2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F2xz_vrr;
    Double I_KINETIC_F3y_F2xz_vrr = PAY*I_KINETIC_D2y_F2xz_vrr+2*oned2z*I_KINETIC_Py_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F2xz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F2xz_vrr;
    Double I_KINETIC_F2yz_F2xz_vrr = PAZ*I_KINETIC_D2y_F2xz_vrr+oned2z*I_KINETIC_D2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F2xz_vrr;
    Double I_KINETIC_F3z_F2xz_vrr = PAZ*I_KINETIC_D2z_F2xz_vrr+2*oned2z*I_KINETIC_Pz_F2xz_vrr+oned2z*I_KINETIC_D2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F2xz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F2xz_vrr;
    Double I_KINETIC_F3x_Fx2y_vrr = PAX*I_KINETIC_D2x_Fx2y_vrr+2*oned2z*I_KINETIC_Px_Fx2y_vrr+oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fx2y_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Fx2y_vrr;
    Double I_KINETIC_F2xy_Fx2y_vrr = PAY*I_KINETIC_D2x_Fx2y_vrr+2*oned2z*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fx2y_vrr;
    Double I_KINETIC_F2xz_Fx2y_vrr = PAZ*I_KINETIC_D2x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fx2y_vrr;
    Double I_KINETIC_Fx2y_Fx2y_vrr = PAX*I_KINETIC_D2y_Fx2y_vrr+oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fx2y_vrr;
    Double I_KINETIC_Fx2z_Fx2y_vrr = PAX*I_KINETIC_D2z_Fx2y_vrr+oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fx2y_vrr;
    Double I_KINETIC_F3y_Fx2y_vrr = PAY*I_KINETIC_D2y_Fx2y_vrr+2*oned2z*I_KINETIC_Py_Fx2y_vrr+2*oned2z*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fx2y_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fx2y_vrr;
    Double I_KINETIC_F2yz_Fx2y_vrr = PAZ*I_KINETIC_D2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fx2y_vrr;
    Double I_KINETIC_F3z_Fx2y_vrr = PAZ*I_KINETIC_D2z_Fx2y_vrr+2*oned2z*I_KINETIC_Pz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fx2y_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Fx2y_vrr;
    Double I_KINETIC_F3x_Fxyz_vrr = PAX*I_KINETIC_D2x_Fxyz_vrr+2*oned2z*I_KINETIC_Px_Fxyz_vrr+oned2z*I_KINETIC_D2x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Fxyz_vrr;
    Double I_KINETIC_F2xy_Fxyz_vrr = PAY*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fxyz_vrr;
    Double I_KINETIC_F2xz_Fxyz_vrr = PAZ*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_D2x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fxyz_vrr;
    Double I_KINETIC_Fx2y_Fxyz_vrr = PAX*I_KINETIC_D2y_Fxyz_vrr+oned2z*I_KINETIC_D2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fxyz_vrr;
    Double I_KINETIC_Fx2z_Fxyz_vrr = PAX*I_KINETIC_D2z_Fxyz_vrr+oned2z*I_KINETIC_D2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fxyz_vrr;
    Double I_KINETIC_F3y_Fxyz_vrr = PAY*I_KINETIC_D2y_Fxyz_vrr+2*oned2z*I_KINETIC_Py_Fxyz_vrr+oned2z*I_KINETIC_D2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fxyz_vrr;
    Double I_KINETIC_F2yz_Fxyz_vrr = PAZ*I_KINETIC_D2y_Fxyz_vrr+oned2z*I_KINETIC_D2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fxyz_vrr;
    Double I_KINETIC_F3z_Fxyz_vrr = PAZ*I_KINETIC_D2z_Fxyz_vrr+2*oned2z*I_KINETIC_Pz_Fxyz_vrr+oned2z*I_KINETIC_D2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fxyz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Fxyz_vrr;
    Double I_KINETIC_F3x_Fx2z_vrr = PAX*I_KINETIC_D2x_Fx2z_vrr+2*oned2z*I_KINETIC_Px_Fx2z_vrr+oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fx2z_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Fx2z_vrr;
    Double I_KINETIC_F2xy_Fx2z_vrr = PAY*I_KINETIC_D2x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fx2z_vrr;
    Double I_KINETIC_F2xz_Fx2z_vrr = PAZ*I_KINETIC_D2x_Fx2z_vrr+2*oned2z*I_KINETIC_D2x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fx2z_vrr;
    Double I_KINETIC_Fx2y_Fx2z_vrr = PAX*I_KINETIC_D2y_Fx2z_vrr+oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fx2z_vrr;
    Double I_KINETIC_Fx2z_Fx2z_vrr = PAX*I_KINETIC_D2z_Fx2z_vrr+oned2z*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fx2z_vrr;
    Double I_KINETIC_F3y_Fx2z_vrr = PAY*I_KINETIC_D2y_Fx2z_vrr+2*oned2z*I_KINETIC_Py_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fx2z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fx2z_vrr;
    Double I_KINETIC_F2yz_Fx2z_vrr = PAZ*I_KINETIC_D2y_Fx2z_vrr+2*oned2z*I_KINETIC_D2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fx2z_vrr;
    Double I_KINETIC_F3z_Fx2z_vrr = PAZ*I_KINETIC_D2z_Fx2z_vrr+2*oned2z*I_KINETIC_Pz_Fx2z_vrr+2*oned2z*I_KINETIC_D2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fx2z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Fx2z_vrr;
    Double I_KINETIC_F3x_F3y_vrr = PAX*I_KINETIC_D2x_F3y_vrr+2*oned2z*I_KINETIC_Px_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F3y_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F3y_vrr;
    Double I_KINETIC_F2xy_F3y_vrr = PAY*I_KINETIC_D2x_F3y_vrr+3*oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F3y_vrr;
    Double I_KINETIC_F2xz_F3y_vrr = PAZ*I_KINETIC_D2x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F3y_vrr;
    Double I_KINETIC_Fx2y_F3y_vrr = PAX*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F3y_vrr;
    Double I_KINETIC_Fx2z_F3y_vrr = PAX*I_KINETIC_D2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F3y_vrr;
    Double I_KINETIC_F3y_F3y_vrr = PAY*I_KINETIC_D2y_F3y_vrr+2*oned2z*I_KINETIC_Py_F3y_vrr+3*oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F3y_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F3y_vrr;
    Double I_KINETIC_F2yz_F3y_vrr = PAZ*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F3y_vrr;
    Double I_KINETIC_F3z_F3y_vrr = PAZ*I_KINETIC_D2z_F3y_vrr+2*oned2z*I_KINETIC_Pz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F3y_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F3y_vrr;
    Double I_KINETIC_F3x_F2yz_vrr = PAX*I_KINETIC_D2x_F2yz_vrr+2*oned2z*I_KINETIC_Px_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F2yz_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F2yz_vrr;
    Double I_KINETIC_F2xy_F2yz_vrr = PAY*I_KINETIC_D2x_F2yz_vrr+2*oned2z*I_KINETIC_D2x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F2yz_vrr;
    Double I_KINETIC_F2xz_F2yz_vrr = PAZ*I_KINETIC_D2x_F2yz_vrr+oned2z*I_KINETIC_D2x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F2yz_vrr;
    Double I_KINETIC_Fx2y_F2yz_vrr = PAX*I_KINETIC_D2y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F2yz_vrr;
    Double I_KINETIC_Fx2z_F2yz_vrr = PAX*I_KINETIC_D2z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F2yz_vrr;
    Double I_KINETIC_F3y_F2yz_vrr = PAY*I_KINETIC_D2y_F2yz_vrr+2*oned2z*I_KINETIC_Py_F2yz_vrr+2*oned2z*I_KINETIC_D2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F2yz_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F2yz_vrr;
    Double I_KINETIC_F2yz_F2yz_vrr = PAZ*I_KINETIC_D2y_F2yz_vrr+oned2z*I_KINETIC_D2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F2yz_vrr;
    Double I_KINETIC_F3z_F2yz_vrr = PAZ*I_KINETIC_D2z_F2yz_vrr+2*oned2z*I_KINETIC_Pz_F2yz_vrr+oned2z*I_KINETIC_D2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F2yz_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F2yz_vrr;
    Double I_KINETIC_F3x_Fy2z_vrr = PAX*I_KINETIC_D2x_Fy2z_vrr+2*oned2z*I_KINETIC_Px_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_Px_Fy2z_vrr;
    Double I_KINETIC_F2xy_Fy2z_vrr = PAY*I_KINETIC_D2x_Fy2z_vrr+oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_Fy2z_vrr;
    Double I_KINETIC_F2xz_Fy2z_vrr = PAZ*I_KINETIC_D2x_Fy2z_vrr+2*oned2z*I_KINETIC_D2x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_Fy2z_vrr;
    Double I_KINETIC_Fx2y_Fy2z_vrr = PAX*I_KINETIC_D2y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_Fy2z_vrr;
    Double I_KINETIC_Fx2z_Fy2z_vrr = PAX*I_KINETIC_D2z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_Fy2z_vrr;
    Double I_KINETIC_F3y_Fy2z_vrr = PAY*I_KINETIC_D2y_Fy2z_vrr+2*oned2z*I_KINETIC_Py_Fy2z_vrr+oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_Fy2z_vrr;
    Double I_KINETIC_F2yz_Fy2z_vrr = PAZ*I_KINETIC_D2y_Fy2z_vrr+2*oned2z*I_KINETIC_D2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_Fy2z_vrr;
    Double I_KINETIC_F3z_Fy2z_vrr = PAZ*I_KINETIC_D2z_Fy2z_vrr+2*oned2z*I_KINETIC_Pz_Fy2z_vrr+2*oned2z*I_KINETIC_D2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_F3z_Fy2z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_Fy2z_vrr;
    Double I_KINETIC_F3x_F3z_vrr = PAX*I_KINETIC_D2x_F3z_vrr+2*oned2z*I_KINETIC_Px_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3x_F3z_vrr-2*bdz*I_TWOBODYOVERLAP_Px_F3z_vrr;
    Double I_KINETIC_F2xy_F3z_vrr = PAY*I_KINETIC_D2x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_F3z_vrr;
    Double I_KINETIC_F2xz_F3z_vrr = PAZ*I_KINETIC_D2x_F3z_vrr+3*oned2z*I_KINETIC_D2x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_F3z_vrr;
    Double I_KINETIC_Fx2y_F3z_vrr = PAX*I_KINETIC_D2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_F3z_vrr;
    Double I_KINETIC_Fx2z_F3z_vrr = PAX*I_KINETIC_D2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_F3z_vrr;
    Double I_KINETIC_F3y_F3z_vrr = PAY*I_KINETIC_D2y_F3z_vrr+2*oned2z*I_KINETIC_Py_F3z_vrr+twoxi*I_TWOBODYOVERLAP_F3y_F3z_vrr-2*bdz*I_TWOBODYOVERLAP_Py_F3z_vrr;
    Double I_KINETIC_F2yz_F3z_vrr = PAZ*I_KINETIC_D2y_F3z_vrr+3*oned2z*I_KINETIC_D2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_F3z_vrr;
    Double I_KINETIC_F3z_F3z_vrr = PAZ*I_KINETIC_D2z_F3z_vrr+2*oned2z*I_KINETIC_Pz_F3z_vrr+3*oned2z*I_KINETIC_D2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_F3z_F3z_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_F
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_F
     * RHS shell quartet name: SQ_KINETIC_D_F
     * RHS shell quartet name: SQ_KINETIC_F_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
     ************************************************************/
    Double I_KINETIC_G4x_F3x_vrr = PAX*I_KINETIC_F3x_F3x_vrr+3*oned2z*I_KINETIC_D2x_F3x_vrr+3*oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3x_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_KINETIC_G3xy_F3x_vrr = PAY*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3x_vrr;
    Double I_KINETIC_G3xz_F3x_vrr = PAZ*I_KINETIC_F3x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3x_vrr;
    Double I_KINETIC_G2x2y_F3x_vrr = PAY*I_KINETIC_F2xy_F3x_vrr+oned2z*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3x_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_KINETIC_G2xyz_F3x_vrr = PAZ*I_KINETIC_F2xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F3x_vrr;
    Double I_KINETIC_G2x2z_F3x_vrr = PAZ*I_KINETIC_F2xz_F3x_vrr+oned2z*I_KINETIC_D2x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3x_vrr;
    Double I_KINETIC_Gx3y_F3x_vrr = PAX*I_KINETIC_F3y_F3x_vrr+3*oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3x_vrr;
    Double I_KINETIC_Gx2yz_F3x_vrr = PAZ*I_KINETIC_Fx2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F3x_vrr;
    Double I_KINETIC_Gxy2z_F3x_vrr = PAY*I_KINETIC_Fx2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F3x_vrr;
    Double I_KINETIC_Gx3z_F3x_vrr = PAX*I_KINETIC_F3z_F3x_vrr+3*oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3x_vrr;
    Double I_KINETIC_G4y_F3x_vrr = PAY*I_KINETIC_F3y_F3x_vrr+3*oned2z*I_KINETIC_D2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3x_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_KINETIC_G3yz_F3x_vrr = PAZ*I_KINETIC_F3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3x_vrr;
    Double I_KINETIC_G2y2z_F3x_vrr = PAZ*I_KINETIC_F2yz_F3x_vrr+oned2z*I_KINETIC_D2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3x_vrr-bdz*I_TWOBODYOVERLAP_D2y_F3x_vrr;
    Double I_KINETIC_Gy3z_F3x_vrr = PAY*I_KINETIC_F3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3x_vrr;
    Double I_KINETIC_G4z_F3x_vrr = PAZ*I_KINETIC_F3z_F3x_vrr+3*oned2z*I_KINETIC_D2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3x_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F3x_vrr;
    Double I_KINETIC_G4x_F2xy_vrr = PAX*I_KINETIC_F3x_F2xy_vrr+3*oned2z*I_KINETIC_D2x_F2xy_vrr+2*oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2xy_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_KINETIC_G3xy_F2xy_vrr = PAY*I_KINETIC_F3x_F2xy_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2xy_vrr;
    Double I_KINETIC_G3xz_F2xy_vrr = PAZ*I_KINETIC_F3x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2xy_vrr;
    Double I_KINETIC_G2x2y_F2xy_vrr = PAY*I_KINETIC_F2xy_F2xy_vrr+oned2z*I_KINETIC_D2x_F2xy_vrr+oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_KINETIC_G2xyz_F2xy_vrr = PAZ*I_KINETIC_F2xy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F2xy_vrr;
    Double I_KINETIC_G2x2z_F2xy_vrr = PAZ*I_KINETIC_F2xz_F2xy_vrr+oned2z*I_KINETIC_D2x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2xy_vrr;
    Double I_KINETIC_Gx3y_F2xy_vrr = PAX*I_KINETIC_F3y_F2xy_vrr+2*oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr;
    Double I_KINETIC_Gx2yz_F2xy_vrr = PAZ*I_KINETIC_Fx2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F2xy_vrr;
    Double I_KINETIC_Gxy2z_F2xy_vrr = PAY*I_KINETIC_Fx2z_F2xy_vrr+oned2z*I_KINETIC_Fx2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F2xy_vrr;
    Double I_KINETIC_Gx3z_F2xy_vrr = PAX*I_KINETIC_F3z_F2xy_vrr+2*oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr;
    Double I_KINETIC_G4y_F2xy_vrr = PAY*I_KINETIC_F3y_F2xy_vrr+3*oned2z*I_KINETIC_D2y_F2xy_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2xy_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_KINETIC_G3yz_F2xy_vrr = PAZ*I_KINETIC_F3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2xy_vrr;
    Double I_KINETIC_G2y2z_F2xy_vrr = PAZ*I_KINETIC_F2yz_F2xy_vrr+oned2z*I_KINETIC_D2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2xy_vrr-bdz*I_TWOBODYOVERLAP_D2y_F2xy_vrr;
    Double I_KINETIC_Gy3z_F2xy_vrr = PAY*I_KINETIC_F3z_F2xy_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2xy_vrr;
    Double I_KINETIC_G4z_F2xy_vrr = PAZ*I_KINETIC_F3z_F2xy_vrr+3*oned2z*I_KINETIC_D2z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2xy_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F2xy_vrr;
    Double I_KINETIC_G4x_F2xz_vrr = PAX*I_KINETIC_F3x_F2xz_vrr+3*oned2z*I_KINETIC_D2x_F2xz_vrr+2*oned2z*I_KINETIC_F3x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2xz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_KINETIC_G3xy_F2xz_vrr = PAY*I_KINETIC_F3x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2xz_vrr;
    Double I_KINETIC_G3xz_F2xz_vrr = PAZ*I_KINETIC_F3x_F2xz_vrr+oned2z*I_KINETIC_F3x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2xz_vrr;
    Double I_KINETIC_G2x2y_F2xz_vrr = PAY*I_KINETIC_F2xy_F2xz_vrr+oned2z*I_KINETIC_D2x_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2xz_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_KINETIC_G2xyz_F2xz_vrr = PAZ*I_KINETIC_F2xy_F2xz_vrr+oned2z*I_KINETIC_F2xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F2xz_vrr;
    Double I_KINETIC_G2x2z_F2xz_vrr = PAZ*I_KINETIC_F2xz_F2xz_vrr+oned2z*I_KINETIC_D2x_F2xz_vrr+oned2z*I_KINETIC_F2xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2xz_vrr;
    Double I_KINETIC_Gx3y_F2xz_vrr = PAX*I_KINETIC_F3y_F2xz_vrr+2*oned2z*I_KINETIC_F3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2xz_vrr;
    Double I_KINETIC_Gx2yz_F2xz_vrr = PAZ*I_KINETIC_Fx2y_F2xz_vrr+oned2z*I_KINETIC_Fx2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F2xz_vrr;
    Double I_KINETIC_Gxy2z_F2xz_vrr = PAY*I_KINETIC_Fx2z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F2xz_vrr;
    Double I_KINETIC_Gx3z_F2xz_vrr = PAX*I_KINETIC_F3z_F2xz_vrr+2*oned2z*I_KINETIC_F3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2xz_vrr;
    Double I_KINETIC_G4y_F2xz_vrr = PAY*I_KINETIC_F3y_F2xz_vrr+3*oned2z*I_KINETIC_D2y_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2xz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F2xz_vrr;
    Double I_KINETIC_G3yz_F2xz_vrr = PAZ*I_KINETIC_F3y_F2xz_vrr+oned2z*I_KINETIC_F3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2xz_vrr;
    Double I_KINETIC_G2y2z_F2xz_vrr = PAZ*I_KINETIC_F2yz_F2xz_vrr+oned2z*I_KINETIC_D2y_F2xz_vrr+oned2z*I_KINETIC_F2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2xz_vrr-bdz*I_TWOBODYOVERLAP_D2y_F2xz_vrr;
    Double I_KINETIC_Gy3z_F2xz_vrr = PAY*I_KINETIC_F3z_F2xz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2xz_vrr;
    Double I_KINETIC_G4z_F2xz_vrr = PAZ*I_KINETIC_F3z_F2xz_vrr+3*oned2z*I_KINETIC_D2z_F2xz_vrr+oned2z*I_KINETIC_F3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2xz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F2xz_vrr;
    Double I_KINETIC_G4x_Fx2y_vrr = PAX*I_KINETIC_F3x_Fx2y_vrr+3*oned2z*I_KINETIC_D2x_Fx2y_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fx2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_KINETIC_G3xy_Fx2y_vrr = PAY*I_KINETIC_F3x_Fx2y_vrr+2*oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr;
    Double I_KINETIC_G3xz_Fx2y_vrr = PAZ*I_KINETIC_F3x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fx2y_vrr;
    Double I_KINETIC_G2x2y_Fx2y_vrr = PAY*I_KINETIC_F2xy_Fx2y_vrr+oned2z*I_KINETIC_D2x_Fx2y_vrr+2*oned2z*I_KINETIC_F2xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_KINETIC_G2xyz_Fx2y_vrr = PAZ*I_KINETIC_F2xy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fx2y_vrr;
    Double I_KINETIC_G2x2z_Fx2y_vrr = PAZ*I_KINETIC_F2xz_Fx2y_vrr+oned2z*I_KINETIC_D2x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fx2y_vrr;
    Double I_KINETIC_Gx3y_Fx2y_vrr = PAX*I_KINETIC_F3y_Fx2y_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr;
    Double I_KINETIC_Gx2yz_Fx2y_vrr = PAZ*I_KINETIC_Fx2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fx2y_vrr;
    Double I_KINETIC_Gxy2z_Fx2y_vrr = PAY*I_KINETIC_Fx2z_Fx2y_vrr+2*oned2z*I_KINETIC_Fx2z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fx2y_vrr;
    Double I_KINETIC_Gx3z_Fx2y_vrr = PAX*I_KINETIC_F3z_Fx2y_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr;
    Double I_KINETIC_G4y_Fx2y_vrr = PAY*I_KINETIC_F3y_Fx2y_vrr+3*oned2z*I_KINETIC_D2y_Fx2y_vrr+2*oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fx2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_KINETIC_G3yz_Fx2y_vrr = PAZ*I_KINETIC_F3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fx2y_vrr;
    Double I_KINETIC_G2y2z_Fx2y_vrr = PAZ*I_KINETIC_F2yz_Fx2y_vrr+oned2z*I_KINETIC_D2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr-bdz*I_TWOBODYOVERLAP_D2y_Fx2y_vrr;
    Double I_KINETIC_Gy3z_Fx2y_vrr = PAY*I_KINETIC_F3z_Fx2y_vrr+2*oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr;
    Double I_KINETIC_G4z_Fx2y_vrr = PAZ*I_KINETIC_F3z_Fx2y_vrr+3*oned2z*I_KINETIC_D2z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fx2y_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Fx2y_vrr;
    Double I_KINETIC_G4x_Fxyz_vrr = PAX*I_KINETIC_F3x_Fxyz_vrr+3*oned2z*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_F3x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fxyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Fxyz_vrr;
    Double I_KINETIC_G3xy_Fxyz_vrr = PAY*I_KINETIC_F3x_Fxyz_vrr+oned2z*I_KINETIC_F3x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fxyz_vrr;
    Double I_KINETIC_G3xz_Fxyz_vrr = PAZ*I_KINETIC_F3x_Fxyz_vrr+oned2z*I_KINETIC_F3x_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fxyz_vrr;
    Double I_KINETIC_G2x2y_Fxyz_vrr = PAY*I_KINETIC_F2xy_Fxyz_vrr+oned2z*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_F2xy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fxyz_vrr;
    Double I_KINETIC_G2xyz_Fxyz_vrr = PAZ*I_KINETIC_F2xy_Fxyz_vrr+oned2z*I_KINETIC_F2xy_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fxyz_vrr;
    Double I_KINETIC_G2x2z_Fxyz_vrr = PAZ*I_KINETIC_F2xz_Fxyz_vrr+oned2z*I_KINETIC_D2x_Fxyz_vrr+oned2z*I_KINETIC_F2xz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fxyz_vrr;
    Double I_KINETIC_Gx3y_Fxyz_vrr = PAX*I_KINETIC_F3y_Fxyz_vrr+oned2z*I_KINETIC_F3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fxyz_vrr;
    Double I_KINETIC_Gx2yz_Fxyz_vrr = PAZ*I_KINETIC_Fx2y_Fxyz_vrr+oned2z*I_KINETIC_Fx2y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fxyz_vrr;
    Double I_KINETIC_Gxy2z_Fxyz_vrr = PAY*I_KINETIC_Fx2z_Fxyz_vrr+oned2z*I_KINETIC_Fx2z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fxyz_vrr;
    Double I_KINETIC_Gx3z_Fxyz_vrr = PAX*I_KINETIC_F3z_Fxyz_vrr+oned2z*I_KINETIC_F3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fxyz_vrr;
    Double I_KINETIC_G4y_Fxyz_vrr = PAY*I_KINETIC_F3y_Fxyz_vrr+3*oned2z*I_KINETIC_D2y_Fxyz_vrr+oned2z*I_KINETIC_F3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fxyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Fxyz_vrr;
    Double I_KINETIC_G3yz_Fxyz_vrr = PAZ*I_KINETIC_F3y_Fxyz_vrr+oned2z*I_KINETIC_F3y_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fxyz_vrr;
    Double I_KINETIC_G2y2z_Fxyz_vrr = PAZ*I_KINETIC_F2yz_Fxyz_vrr+oned2z*I_KINETIC_D2y_Fxyz_vrr+oned2z*I_KINETIC_F2yz_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fxyz_vrr-bdz*I_TWOBODYOVERLAP_D2y_Fxyz_vrr;
    Double I_KINETIC_Gy3z_Fxyz_vrr = PAY*I_KINETIC_F3z_Fxyz_vrr+oned2z*I_KINETIC_F3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fxyz_vrr;
    Double I_KINETIC_G4z_Fxyz_vrr = PAZ*I_KINETIC_F3z_Fxyz_vrr+3*oned2z*I_KINETIC_D2z_Fxyz_vrr+oned2z*I_KINETIC_F3z_Dxy_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fxyz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Fxyz_vrr;
    Double I_KINETIC_G4x_Fx2z_vrr = PAX*I_KINETIC_F3x_Fx2z_vrr+3*oned2z*I_KINETIC_D2x_Fx2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fx2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_KINETIC_G3xy_Fx2z_vrr = PAY*I_KINETIC_F3x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr;
    Double I_KINETIC_G3xz_Fx2z_vrr = PAZ*I_KINETIC_F3x_Fx2z_vrr+2*oned2z*I_KINETIC_F3x_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fx2z_vrr;
    Double I_KINETIC_G2x2y_Fx2z_vrr = PAY*I_KINETIC_F2xy_Fx2z_vrr+oned2z*I_KINETIC_D2x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_KINETIC_G2xyz_Fx2z_vrr = PAZ*I_KINETIC_F2xy_Fx2z_vrr+2*oned2z*I_KINETIC_F2xy_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fx2z_vrr;
    Double I_KINETIC_G2x2z_Fx2z_vrr = PAZ*I_KINETIC_F2xz_Fx2z_vrr+oned2z*I_KINETIC_D2x_Fx2z_vrr+2*oned2z*I_KINETIC_F2xz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fx2z_vrr;
    Double I_KINETIC_Gx3y_Fx2z_vrr = PAX*I_KINETIC_F3y_Fx2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr;
    Double I_KINETIC_Gx2yz_Fx2z_vrr = PAZ*I_KINETIC_Fx2y_Fx2z_vrr+2*oned2z*I_KINETIC_Fx2y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fx2z_vrr;
    Double I_KINETIC_Gxy2z_Fx2z_vrr = PAY*I_KINETIC_Fx2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fx2z_vrr;
    Double I_KINETIC_Gx3z_Fx2z_vrr = PAX*I_KINETIC_F3z_Fx2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr;
    Double I_KINETIC_G4y_Fx2z_vrr = PAY*I_KINETIC_F3y_Fx2z_vrr+3*oned2z*I_KINETIC_D2y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fx2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_KINETIC_G3yz_Fx2z_vrr = PAZ*I_KINETIC_F3y_Fx2z_vrr+2*oned2z*I_KINETIC_F3y_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fx2z_vrr;
    Double I_KINETIC_G2y2z_Fx2z_vrr = PAZ*I_KINETIC_F2yz_Fx2z_vrr+oned2z*I_KINETIC_D2y_Fx2z_vrr+2*oned2z*I_KINETIC_F2yz_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr-bdz*I_TWOBODYOVERLAP_D2y_Fx2z_vrr;
    Double I_KINETIC_Gy3z_Fx2z_vrr = PAY*I_KINETIC_F3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr;
    Double I_KINETIC_G4z_Fx2z_vrr = PAZ*I_KINETIC_F3z_Fx2z_vrr+3*oned2z*I_KINETIC_D2z_Fx2z_vrr+2*oned2z*I_KINETIC_F3z_Dxz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fx2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Fx2z_vrr;
    Double I_KINETIC_G4x_F3y_vrr = PAX*I_KINETIC_F3x_F3y_vrr+3*oned2z*I_KINETIC_D2x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3y_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_KINETIC_G3xy_F3y_vrr = PAY*I_KINETIC_F3x_F3y_vrr+3*oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3y_vrr;
    Double I_KINETIC_G3xz_F3y_vrr = PAZ*I_KINETIC_F3x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3y_vrr;
    Double I_KINETIC_G2x2y_F3y_vrr = PAY*I_KINETIC_F2xy_F3y_vrr+oned2z*I_KINETIC_D2x_F3y_vrr+3*oned2z*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3y_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_KINETIC_G2xyz_F3y_vrr = PAZ*I_KINETIC_F2xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F3y_vrr;
    Double I_KINETIC_G2x2z_F3y_vrr = PAZ*I_KINETIC_F2xz_F3y_vrr+oned2z*I_KINETIC_D2x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3y_vrr;
    Double I_KINETIC_Gx3y_F3y_vrr = PAX*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3y_vrr;
    Double I_KINETIC_Gx2yz_F3y_vrr = PAZ*I_KINETIC_Fx2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F3y_vrr;
    Double I_KINETIC_Gxy2z_F3y_vrr = PAY*I_KINETIC_Fx2z_F3y_vrr+3*oned2z*I_KINETIC_Fx2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F3y_vrr;
    Double I_KINETIC_Gx3z_F3y_vrr = PAX*I_KINETIC_F3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3y_vrr;
    Double I_KINETIC_G4y_F3y_vrr = PAY*I_KINETIC_F3y_F3y_vrr+3*oned2z*I_KINETIC_D2y_F3y_vrr+3*oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3y_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_KINETIC_G3yz_F3y_vrr = PAZ*I_KINETIC_F3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3y_vrr;
    Double I_KINETIC_G2y2z_F3y_vrr = PAZ*I_KINETIC_F2yz_F3y_vrr+oned2z*I_KINETIC_D2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3y_vrr-bdz*I_TWOBODYOVERLAP_D2y_F3y_vrr;
    Double I_KINETIC_Gy3z_F3y_vrr = PAY*I_KINETIC_F3z_F3y_vrr+3*oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3y_vrr;
    Double I_KINETIC_G4z_F3y_vrr = PAZ*I_KINETIC_F3z_F3y_vrr+3*oned2z*I_KINETIC_D2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3y_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F3y_vrr;
    Double I_KINETIC_G4x_F2yz_vrr = PAX*I_KINETIC_F3x_F2yz_vrr+3*oned2z*I_KINETIC_D2x_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F2yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_KINETIC_G3xy_F2yz_vrr = PAY*I_KINETIC_F3x_F2yz_vrr+2*oned2z*I_KINETIC_F3x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F2yz_vrr;
    Double I_KINETIC_G3xz_F2yz_vrr = PAZ*I_KINETIC_F3x_F2yz_vrr+oned2z*I_KINETIC_F3x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F2yz_vrr;
    Double I_KINETIC_G2x2y_F2yz_vrr = PAY*I_KINETIC_F2xy_F2yz_vrr+oned2z*I_KINETIC_D2x_F2yz_vrr+2*oned2z*I_KINETIC_F2xy_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F2yz_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_KINETIC_G2xyz_F2yz_vrr = PAZ*I_KINETIC_F2xy_F2yz_vrr+oned2z*I_KINETIC_F2xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F2yz_vrr;
    Double I_KINETIC_G2x2z_F2yz_vrr = PAZ*I_KINETIC_F2xz_F2yz_vrr+oned2z*I_KINETIC_D2x_F2yz_vrr+oned2z*I_KINETIC_F2xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_D2x_F2yz_vrr;
    Double I_KINETIC_Gx3y_F2yz_vrr = PAX*I_KINETIC_F3y_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F2yz_vrr;
    Double I_KINETIC_Gx2yz_F2yz_vrr = PAZ*I_KINETIC_Fx2y_F2yz_vrr+oned2z*I_KINETIC_Fx2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F2yz_vrr;
    Double I_KINETIC_Gxy2z_F2yz_vrr = PAY*I_KINETIC_Fx2z_F2yz_vrr+2*oned2z*I_KINETIC_Fx2z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F2yz_vrr;
    Double I_KINETIC_Gx3z_F2yz_vrr = PAX*I_KINETIC_F3z_F2yz_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F2yz_vrr;
    Double I_KINETIC_G4y_F2yz_vrr = PAY*I_KINETIC_F3y_F2yz_vrr+3*oned2z*I_KINETIC_D2y_F2yz_vrr+2*oned2z*I_KINETIC_F3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F2yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F2yz_vrr;
    Double I_KINETIC_G3yz_F2yz_vrr = PAZ*I_KINETIC_F3y_F2yz_vrr+oned2z*I_KINETIC_F3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F2yz_vrr;
    Double I_KINETIC_G2y2z_F2yz_vrr = PAZ*I_KINETIC_F2yz_F2yz_vrr+oned2z*I_KINETIC_D2y_F2yz_vrr+oned2z*I_KINETIC_F2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F2yz_vrr-bdz*I_TWOBODYOVERLAP_D2y_F2yz_vrr;
    Double I_KINETIC_Gy3z_F2yz_vrr = PAY*I_KINETIC_F3z_F2yz_vrr+2*oned2z*I_KINETIC_F3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F2yz_vrr;
    Double I_KINETIC_G4z_F2yz_vrr = PAZ*I_KINETIC_F3z_F2yz_vrr+3*oned2z*I_KINETIC_D2z_F2yz_vrr+oned2z*I_KINETIC_F3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F2yz_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F2yz_vrr;
    Double I_KINETIC_G4x_Fy2z_vrr = PAX*I_KINETIC_F3x_Fy2z_vrr+3*oned2z*I_KINETIC_D2x_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Fy2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_KINETIC_G3xy_Fy2z_vrr = PAY*I_KINETIC_F3x_Fy2z_vrr+oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Fy2z_vrr;
    Double I_KINETIC_G3xz_Fy2z_vrr = PAZ*I_KINETIC_F3x_Fy2z_vrr+2*oned2z*I_KINETIC_F3x_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Fy2z_vrr;
    Double I_KINETIC_G2x2y_Fy2z_vrr = PAY*I_KINETIC_F2xy_Fy2z_vrr+oned2z*I_KINETIC_D2x_Fy2z_vrr+oned2z*I_KINETIC_F2xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_KINETIC_G2xyz_Fy2z_vrr = PAZ*I_KINETIC_F2xy_Fy2z_vrr+2*oned2z*I_KINETIC_F2xy_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Fy2z_vrr;
    Double I_KINETIC_G2x2z_Fy2z_vrr = PAZ*I_KINETIC_F2xz_Fy2z_vrr+oned2z*I_KINETIC_D2x_Fy2z_vrr+2*oned2z*I_KINETIC_F2xz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_D2x_Fy2z_vrr;
    Double I_KINETIC_Gx3y_Fy2z_vrr = PAX*I_KINETIC_F3y_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Fy2z_vrr;
    Double I_KINETIC_Gx2yz_Fy2z_vrr = PAZ*I_KINETIC_Fx2y_Fy2z_vrr+2*oned2z*I_KINETIC_Fx2y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Fy2z_vrr;
    Double I_KINETIC_Gxy2z_Fy2z_vrr = PAY*I_KINETIC_Fx2z_Fy2z_vrr+oned2z*I_KINETIC_Fx2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Fy2z_vrr;
    Double I_KINETIC_Gx3z_Fy2z_vrr = PAX*I_KINETIC_F3z_Fy2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Fy2z_vrr;
    Double I_KINETIC_G4y_Fy2z_vrr = PAY*I_KINETIC_F3y_Fy2z_vrr+3*oned2z*I_KINETIC_D2y_Fy2z_vrr+oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Fy2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_Fy2z_vrr;
    Double I_KINETIC_G3yz_Fy2z_vrr = PAZ*I_KINETIC_F3y_Fy2z_vrr+2*oned2z*I_KINETIC_F3y_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Fy2z_vrr;
    Double I_KINETIC_G2y2z_Fy2z_vrr = PAZ*I_KINETIC_F2yz_Fy2z_vrr+oned2z*I_KINETIC_D2y_Fy2z_vrr+2*oned2z*I_KINETIC_F2yz_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Fy2z_vrr-bdz*I_TWOBODYOVERLAP_D2y_Fy2z_vrr;
    Double I_KINETIC_Gy3z_Fy2z_vrr = PAY*I_KINETIC_F3z_Fy2z_vrr+oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Fy2z_vrr;
    Double I_KINETIC_G4z_Fy2z_vrr = PAZ*I_KINETIC_F3z_Fy2z_vrr+3*oned2z*I_KINETIC_D2z_Fy2z_vrr+2*oned2z*I_KINETIC_F3z_Dyz_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Fy2z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_Fy2z_vrr;
    Double I_KINETIC_G4x_F3z_vrr = PAX*I_KINETIC_F3x_F3z_vrr+3*oned2z*I_KINETIC_D2x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G4x_F3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_KINETIC_G3xy_F3z_vrr = PAY*I_KINETIC_F3x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_F3z_vrr;
    Double I_KINETIC_G3xz_F3z_vrr = PAZ*I_KINETIC_F3x_F3z_vrr+3*oned2z*I_KINETIC_F3x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_F3z_vrr;
    Double I_KINETIC_G2x2y_F3z_vrr = PAY*I_KINETIC_F2xy_F3z_vrr+oned2z*I_KINETIC_D2x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_F3z_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_KINETIC_G2xyz_F3z_vrr = PAZ*I_KINETIC_F2xy_F3z_vrr+3*oned2z*I_KINETIC_F2xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_F3z_vrr;
    Double I_KINETIC_G2x2z_F3z_vrr = PAZ*I_KINETIC_F2xz_F3z_vrr+oned2z*I_KINETIC_D2x_F3z_vrr+3*oned2z*I_KINETIC_F2xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_D2x_F3z_vrr;
    Double I_KINETIC_Gx3y_F3z_vrr = PAX*I_KINETIC_F3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_F3z_vrr;
    Double I_KINETIC_Gx2yz_F3z_vrr = PAZ*I_KINETIC_Fx2y_F3z_vrr+3*oned2z*I_KINETIC_Fx2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_F3z_vrr;
    Double I_KINETIC_Gxy2z_F3z_vrr = PAY*I_KINETIC_Fx2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_F3z_vrr;
    Double I_KINETIC_Gx3z_F3z_vrr = PAX*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_F3z_vrr;
    Double I_KINETIC_G4y_F3z_vrr = PAY*I_KINETIC_F3y_F3z_vrr+3*oned2z*I_KINETIC_D2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_G4y_F3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_KINETIC_G3yz_F3z_vrr = PAZ*I_KINETIC_F3y_F3z_vrr+3*oned2z*I_KINETIC_F3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_F3z_vrr;
    Double I_KINETIC_G2y2z_F3z_vrr = PAZ*I_KINETIC_F2yz_F3z_vrr+oned2z*I_KINETIC_D2y_F3z_vrr+3*oned2z*I_KINETIC_F2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_F3z_vrr-bdz*I_TWOBODYOVERLAP_D2y_F3z_vrr;
    Double I_KINETIC_Gy3z_F3z_vrr = PAY*I_KINETIC_F3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_F3z_vrr;
    Double I_KINETIC_G4z_F3z_vrr = PAZ*I_KINETIC_F3z_F3z_vrr+3*oned2z*I_KINETIC_D2z_F3z_vrr+3*oned2z*I_KINETIC_F3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_G4z_F3z_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_F_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_G_F_a_coefs = alpha;
    I_KINETIC_G4x_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_F3x_vrr;
    I_KINETIC_G3xy_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_F3x_vrr;
    I_KINETIC_G3xz_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_F3x_vrr;
    I_KINETIC_G2x2y_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_F3x_vrr;
    I_KINETIC_G2xyz_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_F3x_vrr;
    I_KINETIC_G2x2z_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_F3x_vrr;
    I_KINETIC_Gx3y_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_F3x_vrr;
    I_KINETIC_Gx2yz_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_F3x_vrr;
    I_KINETIC_Gxy2z_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_F3x_vrr;
    I_KINETIC_Gx3z_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_F3x_vrr;
    I_KINETIC_G4y_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_F3x_vrr;
    I_KINETIC_G3yz_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_F3x_vrr;
    I_KINETIC_G2y2z_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_F3x_vrr;
    I_KINETIC_Gy3z_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_F3x_vrr;
    I_KINETIC_G4z_F3x_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_F3x_vrr;
    I_KINETIC_G4x_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_F2xy_vrr;
    I_KINETIC_G3xy_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_F2xy_vrr;
    I_KINETIC_G3xz_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_F2xy_vrr;
    I_KINETIC_G2x2y_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_F2xy_vrr;
    I_KINETIC_G2xyz_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_F2xy_vrr;
    I_KINETIC_G2x2z_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_F2xy_vrr;
    I_KINETIC_Gx3y_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_F2xy_vrr;
    I_KINETIC_Gx2yz_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_F2xy_vrr;
    I_KINETIC_Gxy2z_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_F2xy_vrr;
    I_KINETIC_Gx3z_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_F2xy_vrr;
    I_KINETIC_G4y_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_F2xy_vrr;
    I_KINETIC_G3yz_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_F2xy_vrr;
    I_KINETIC_G2y2z_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_F2xy_vrr;
    I_KINETIC_Gy3z_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_F2xy_vrr;
    I_KINETIC_G4z_F2xy_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_F2xy_vrr;
    I_KINETIC_G4x_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_F2xz_vrr;
    I_KINETIC_G3xy_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_F2xz_vrr;
    I_KINETIC_G3xz_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_F2xz_vrr;
    I_KINETIC_G2x2y_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_F2xz_vrr;
    I_KINETIC_G2xyz_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_F2xz_vrr;
    I_KINETIC_G2x2z_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_F2xz_vrr;
    I_KINETIC_Gx3y_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_F2xz_vrr;
    I_KINETIC_Gx2yz_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_F2xz_vrr;
    I_KINETIC_Gxy2z_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_F2xz_vrr;
    I_KINETIC_Gx3z_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_F2xz_vrr;
    I_KINETIC_G4y_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_F2xz_vrr;
    I_KINETIC_G3yz_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_F2xz_vrr;
    I_KINETIC_G2y2z_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_F2xz_vrr;
    I_KINETIC_Gy3z_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_F2xz_vrr;
    I_KINETIC_G4z_F2xz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_F2xz_vrr;
    I_KINETIC_G4x_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_Fx2y_vrr;
    I_KINETIC_G3xy_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_Fx2y_vrr;
    I_KINETIC_G3xz_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_Fx2y_vrr;
    I_KINETIC_G2x2y_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_Fx2y_vrr;
    I_KINETIC_G2xyz_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_Fx2y_vrr;
    I_KINETIC_G2x2z_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_Fx2y_vrr;
    I_KINETIC_Gx3y_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_Fx2y_vrr;
    I_KINETIC_Gx2yz_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_Fx2y_vrr;
    I_KINETIC_Gxy2z_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_Fx2y_vrr;
    I_KINETIC_Gx3z_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_Fx2y_vrr;
    I_KINETIC_G4y_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_Fx2y_vrr;
    I_KINETIC_G3yz_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_Fx2y_vrr;
    I_KINETIC_G2y2z_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_Fx2y_vrr;
    I_KINETIC_Gy3z_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_Fx2y_vrr;
    I_KINETIC_G4z_Fx2y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_Fx2y_vrr;
    I_KINETIC_G4x_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_Fxyz_vrr;
    I_KINETIC_G3xy_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_Fxyz_vrr;
    I_KINETIC_G3xz_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_Fxyz_vrr;
    I_KINETIC_G2x2y_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_Fxyz_vrr;
    I_KINETIC_G2xyz_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_Fxyz_vrr;
    I_KINETIC_G2x2z_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_Fxyz_vrr;
    I_KINETIC_Gx3y_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_Fxyz_vrr;
    I_KINETIC_Gx2yz_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_Fxyz_vrr;
    I_KINETIC_Gxy2z_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_Fxyz_vrr;
    I_KINETIC_Gx3z_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_Fxyz_vrr;
    I_KINETIC_G4y_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_Fxyz_vrr;
    I_KINETIC_G3yz_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_Fxyz_vrr;
    I_KINETIC_G2y2z_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_Fxyz_vrr;
    I_KINETIC_Gy3z_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_Fxyz_vrr;
    I_KINETIC_G4z_Fxyz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_Fxyz_vrr;
    I_KINETIC_G4x_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_Fx2z_vrr;
    I_KINETIC_G3xy_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_Fx2z_vrr;
    I_KINETIC_G3xz_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_Fx2z_vrr;
    I_KINETIC_G2x2y_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_Fx2z_vrr;
    I_KINETIC_G2xyz_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_Fx2z_vrr;
    I_KINETIC_G2x2z_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_Fx2z_vrr;
    I_KINETIC_Gx3y_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_Fx2z_vrr;
    I_KINETIC_Gx2yz_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_Fx2z_vrr;
    I_KINETIC_Gxy2z_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_Fx2z_vrr;
    I_KINETIC_Gx3z_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_Fx2z_vrr;
    I_KINETIC_G4y_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_Fx2z_vrr;
    I_KINETIC_G3yz_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_Fx2z_vrr;
    I_KINETIC_G2y2z_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_Fx2z_vrr;
    I_KINETIC_Gy3z_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_Fx2z_vrr;
    I_KINETIC_G4z_Fx2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_Fx2z_vrr;
    I_KINETIC_G4x_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_F3y_vrr;
    I_KINETIC_G3xy_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_F3y_vrr;
    I_KINETIC_G3xz_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_F3y_vrr;
    I_KINETIC_G2x2y_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_F3y_vrr;
    I_KINETIC_G2xyz_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_F3y_vrr;
    I_KINETIC_G2x2z_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_F3y_vrr;
    I_KINETIC_Gx3y_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_F3y_vrr;
    I_KINETIC_Gx2yz_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_F3y_vrr;
    I_KINETIC_Gxy2z_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_F3y_vrr;
    I_KINETIC_Gx3z_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_F3y_vrr;
    I_KINETIC_G4y_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_F3y_vrr;
    I_KINETIC_G3yz_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_F3y_vrr;
    I_KINETIC_G2y2z_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_F3y_vrr;
    I_KINETIC_Gy3z_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_F3y_vrr;
    I_KINETIC_G4z_F3y_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_F3y_vrr;
    I_KINETIC_G4x_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_F2yz_vrr;
    I_KINETIC_G3xy_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_F2yz_vrr;
    I_KINETIC_G3xz_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_F2yz_vrr;
    I_KINETIC_G2x2y_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_F2yz_vrr;
    I_KINETIC_G2xyz_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_F2yz_vrr;
    I_KINETIC_G2x2z_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_F2yz_vrr;
    I_KINETIC_Gx3y_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_F2yz_vrr;
    I_KINETIC_Gx2yz_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_F2yz_vrr;
    I_KINETIC_Gxy2z_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_F2yz_vrr;
    I_KINETIC_Gx3z_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_F2yz_vrr;
    I_KINETIC_G4y_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_F2yz_vrr;
    I_KINETIC_G3yz_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_F2yz_vrr;
    I_KINETIC_G2y2z_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_F2yz_vrr;
    I_KINETIC_Gy3z_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_F2yz_vrr;
    I_KINETIC_G4z_F2yz_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_F2yz_vrr;
    I_KINETIC_G4x_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_Fy2z_vrr;
    I_KINETIC_G3xy_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_Fy2z_vrr;
    I_KINETIC_G3xz_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_Fy2z_vrr;
    I_KINETIC_G2x2y_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_Fy2z_vrr;
    I_KINETIC_G2xyz_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_Fy2z_vrr;
    I_KINETIC_G2x2z_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_Fy2z_vrr;
    I_KINETIC_Gx3y_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_Fy2z_vrr;
    I_KINETIC_Gx2yz_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_Fy2z_vrr;
    I_KINETIC_Gxy2z_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_Fy2z_vrr;
    I_KINETIC_Gx3z_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_Fy2z_vrr;
    I_KINETIC_G4y_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_Fy2z_vrr;
    I_KINETIC_G3yz_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_Fy2z_vrr;
    I_KINETIC_G2y2z_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_Fy2z_vrr;
    I_KINETIC_Gy3z_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_Fy2z_vrr;
    I_KINETIC_G4z_Fy2z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_Fy2z_vrr;
    I_KINETIC_G4x_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4x_F3z_vrr;
    I_KINETIC_G3xy_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xy_F3z_vrr;
    I_KINETIC_G3xz_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3xz_F3z_vrr;
    I_KINETIC_G2x2y_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2y_F3z_vrr;
    I_KINETIC_G2xyz_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2xyz_F3z_vrr;
    I_KINETIC_G2x2z_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2x2z_F3z_vrr;
    I_KINETIC_Gx3y_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3y_F3z_vrr;
    I_KINETIC_Gx2yz_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx2yz_F3z_vrr;
    I_KINETIC_Gxy2z_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gxy2z_F3z_vrr;
    I_KINETIC_Gx3z_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gx3z_F3z_vrr;
    I_KINETIC_G4y_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4y_F3z_vrr;
    I_KINETIC_G3yz_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G3yz_F3z_vrr;
    I_KINETIC_G2y2z_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G2y2z_F3z_vrr;
    I_KINETIC_Gy3z_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_Gy3z_F3z_vrr;
    I_KINETIC_G4z_F3z_a += SQ_KINETIC_G_F_a_coefs*I_KINETIC_G4z_F3z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_F
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_D2x_F3x += I_KINETIC_D2x_F3x_vrr;
    I_KINETIC_Dxy_F3x += I_KINETIC_Dxy_F3x_vrr;
    I_KINETIC_Dxz_F3x += I_KINETIC_Dxz_F3x_vrr;
    I_KINETIC_D2y_F3x += I_KINETIC_D2y_F3x_vrr;
    I_KINETIC_Dyz_F3x += I_KINETIC_Dyz_F3x_vrr;
    I_KINETIC_D2z_F3x += I_KINETIC_D2z_F3x_vrr;
    I_KINETIC_D2x_F2xy += I_KINETIC_D2x_F2xy_vrr;
    I_KINETIC_Dxy_F2xy += I_KINETIC_Dxy_F2xy_vrr;
    I_KINETIC_Dxz_F2xy += I_KINETIC_Dxz_F2xy_vrr;
    I_KINETIC_D2y_F2xy += I_KINETIC_D2y_F2xy_vrr;
    I_KINETIC_Dyz_F2xy += I_KINETIC_Dyz_F2xy_vrr;
    I_KINETIC_D2z_F2xy += I_KINETIC_D2z_F2xy_vrr;
    I_KINETIC_D2x_F2xz += I_KINETIC_D2x_F2xz_vrr;
    I_KINETIC_Dxy_F2xz += I_KINETIC_Dxy_F2xz_vrr;
    I_KINETIC_Dxz_F2xz += I_KINETIC_Dxz_F2xz_vrr;
    I_KINETIC_D2y_F2xz += I_KINETIC_D2y_F2xz_vrr;
    I_KINETIC_Dyz_F2xz += I_KINETIC_Dyz_F2xz_vrr;
    I_KINETIC_D2z_F2xz += I_KINETIC_D2z_F2xz_vrr;
    I_KINETIC_D2x_Fx2y += I_KINETIC_D2x_Fx2y_vrr;
    I_KINETIC_Dxy_Fx2y += I_KINETIC_Dxy_Fx2y_vrr;
    I_KINETIC_Dxz_Fx2y += I_KINETIC_Dxz_Fx2y_vrr;
    I_KINETIC_D2y_Fx2y += I_KINETIC_D2y_Fx2y_vrr;
    I_KINETIC_Dyz_Fx2y += I_KINETIC_Dyz_Fx2y_vrr;
    I_KINETIC_D2z_Fx2y += I_KINETIC_D2z_Fx2y_vrr;
    I_KINETIC_D2x_Fxyz += I_KINETIC_D2x_Fxyz_vrr;
    I_KINETIC_Dxy_Fxyz += I_KINETIC_Dxy_Fxyz_vrr;
    I_KINETIC_Dxz_Fxyz += I_KINETIC_Dxz_Fxyz_vrr;
    I_KINETIC_D2y_Fxyz += I_KINETIC_D2y_Fxyz_vrr;
    I_KINETIC_Dyz_Fxyz += I_KINETIC_Dyz_Fxyz_vrr;
    I_KINETIC_D2z_Fxyz += I_KINETIC_D2z_Fxyz_vrr;
    I_KINETIC_D2x_Fx2z += I_KINETIC_D2x_Fx2z_vrr;
    I_KINETIC_Dxy_Fx2z += I_KINETIC_Dxy_Fx2z_vrr;
    I_KINETIC_Dxz_Fx2z += I_KINETIC_Dxz_Fx2z_vrr;
    I_KINETIC_D2y_Fx2z += I_KINETIC_D2y_Fx2z_vrr;
    I_KINETIC_Dyz_Fx2z += I_KINETIC_Dyz_Fx2z_vrr;
    I_KINETIC_D2z_Fx2z += I_KINETIC_D2z_Fx2z_vrr;
    I_KINETIC_D2x_F3y += I_KINETIC_D2x_F3y_vrr;
    I_KINETIC_Dxy_F3y += I_KINETIC_Dxy_F3y_vrr;
    I_KINETIC_Dxz_F3y += I_KINETIC_Dxz_F3y_vrr;
    I_KINETIC_D2y_F3y += I_KINETIC_D2y_F3y_vrr;
    I_KINETIC_Dyz_F3y += I_KINETIC_Dyz_F3y_vrr;
    I_KINETIC_D2z_F3y += I_KINETIC_D2z_F3y_vrr;
    I_KINETIC_D2x_F2yz += I_KINETIC_D2x_F2yz_vrr;
    I_KINETIC_Dxy_F2yz += I_KINETIC_Dxy_F2yz_vrr;
    I_KINETIC_Dxz_F2yz += I_KINETIC_Dxz_F2yz_vrr;
    I_KINETIC_D2y_F2yz += I_KINETIC_D2y_F2yz_vrr;
    I_KINETIC_Dyz_F2yz += I_KINETIC_Dyz_F2yz_vrr;
    I_KINETIC_D2z_F2yz += I_KINETIC_D2z_F2yz_vrr;
    I_KINETIC_D2x_Fy2z += I_KINETIC_D2x_Fy2z_vrr;
    I_KINETIC_Dxy_Fy2z += I_KINETIC_Dxy_Fy2z_vrr;
    I_KINETIC_Dxz_Fy2z += I_KINETIC_Dxz_Fy2z_vrr;
    I_KINETIC_D2y_Fy2z += I_KINETIC_D2y_Fy2z_vrr;
    I_KINETIC_Dyz_Fy2z += I_KINETIC_Dyz_Fy2z_vrr;
    I_KINETIC_D2z_Fy2z += I_KINETIC_D2z_Fy2z_vrr;
    I_KINETIC_D2x_F3z += I_KINETIC_D2x_F3z_vrr;
    I_KINETIC_Dxy_F3z += I_KINETIC_Dxy_F3z_vrr;
    I_KINETIC_Dxz_F3z += I_KINETIC_Dxz_F3z_vrr;
    I_KINETIC_D2y_F3z += I_KINETIC_D2y_F3z_vrr;
    I_KINETIC_Dyz_F3z += I_KINETIC_Dyz_F3z_vrr;
    I_KINETIC_D2z_F3z += I_KINETIC_D2z_F3z_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_F_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_G_F_a
   * RHS shell quartet name: SQ_KINETIC_D_F
   ************************************************************/
  abcd[0] = 2.0E0*I_KINETIC_G4x_F3x_a-3*I_KINETIC_D2x_F3x;
  abcd[1] = 2.0E0*I_KINETIC_G3xy_F3x_a-2*I_KINETIC_Dxy_F3x;
  abcd[2] = 2.0E0*I_KINETIC_G3xz_F3x_a-2*I_KINETIC_Dxz_F3x;
  abcd[3] = 2.0E0*I_KINETIC_G2x2y_F3x_a-1*I_KINETIC_D2y_F3x;
  abcd[4] = 2.0E0*I_KINETIC_G2xyz_F3x_a-1*I_KINETIC_Dyz_F3x;
  abcd[5] = 2.0E0*I_KINETIC_G2x2z_F3x_a-1*I_KINETIC_D2z_F3x;
  abcd[6] = 2.0E0*I_KINETIC_Gx3y_F3x_a;
  abcd[7] = 2.0E0*I_KINETIC_Gx2yz_F3x_a;
  abcd[8] = 2.0E0*I_KINETIC_Gxy2z_F3x_a;
  abcd[9] = 2.0E0*I_KINETIC_Gx3z_F3x_a;
  abcd[10] = 2.0E0*I_KINETIC_G4x_F2xy_a-3*I_KINETIC_D2x_F2xy;
  abcd[11] = 2.0E0*I_KINETIC_G3xy_F2xy_a-2*I_KINETIC_Dxy_F2xy;
  abcd[12] = 2.0E0*I_KINETIC_G3xz_F2xy_a-2*I_KINETIC_Dxz_F2xy;
  abcd[13] = 2.0E0*I_KINETIC_G2x2y_F2xy_a-1*I_KINETIC_D2y_F2xy;
  abcd[14] = 2.0E0*I_KINETIC_G2xyz_F2xy_a-1*I_KINETIC_Dyz_F2xy;
  abcd[15] = 2.0E0*I_KINETIC_G2x2z_F2xy_a-1*I_KINETIC_D2z_F2xy;
  abcd[16] = 2.0E0*I_KINETIC_Gx3y_F2xy_a;
  abcd[17] = 2.0E0*I_KINETIC_Gx2yz_F2xy_a;
  abcd[18] = 2.0E0*I_KINETIC_Gxy2z_F2xy_a;
  abcd[19] = 2.0E0*I_KINETIC_Gx3z_F2xy_a;
  abcd[20] = 2.0E0*I_KINETIC_G4x_F2xz_a-3*I_KINETIC_D2x_F2xz;
  abcd[21] = 2.0E0*I_KINETIC_G3xy_F2xz_a-2*I_KINETIC_Dxy_F2xz;
  abcd[22] = 2.0E0*I_KINETIC_G3xz_F2xz_a-2*I_KINETIC_Dxz_F2xz;
  abcd[23] = 2.0E0*I_KINETIC_G2x2y_F2xz_a-1*I_KINETIC_D2y_F2xz;
  abcd[24] = 2.0E0*I_KINETIC_G2xyz_F2xz_a-1*I_KINETIC_Dyz_F2xz;
  abcd[25] = 2.0E0*I_KINETIC_G2x2z_F2xz_a-1*I_KINETIC_D2z_F2xz;
  abcd[26] = 2.0E0*I_KINETIC_Gx3y_F2xz_a;
  abcd[27] = 2.0E0*I_KINETIC_Gx2yz_F2xz_a;
  abcd[28] = 2.0E0*I_KINETIC_Gxy2z_F2xz_a;
  abcd[29] = 2.0E0*I_KINETIC_Gx3z_F2xz_a;
  abcd[30] = 2.0E0*I_KINETIC_G4x_Fx2y_a-3*I_KINETIC_D2x_Fx2y;
  abcd[31] = 2.0E0*I_KINETIC_G3xy_Fx2y_a-2*I_KINETIC_Dxy_Fx2y;
  abcd[32] = 2.0E0*I_KINETIC_G3xz_Fx2y_a-2*I_KINETIC_Dxz_Fx2y;
  abcd[33] = 2.0E0*I_KINETIC_G2x2y_Fx2y_a-1*I_KINETIC_D2y_Fx2y;
  abcd[34] = 2.0E0*I_KINETIC_G2xyz_Fx2y_a-1*I_KINETIC_Dyz_Fx2y;
  abcd[35] = 2.0E0*I_KINETIC_G2x2z_Fx2y_a-1*I_KINETIC_D2z_Fx2y;
  abcd[36] = 2.0E0*I_KINETIC_Gx3y_Fx2y_a;
  abcd[37] = 2.0E0*I_KINETIC_Gx2yz_Fx2y_a;
  abcd[38] = 2.0E0*I_KINETIC_Gxy2z_Fx2y_a;
  abcd[39] = 2.0E0*I_KINETIC_Gx3z_Fx2y_a;
  abcd[40] = 2.0E0*I_KINETIC_G4x_Fxyz_a-3*I_KINETIC_D2x_Fxyz;
  abcd[41] = 2.0E0*I_KINETIC_G3xy_Fxyz_a-2*I_KINETIC_Dxy_Fxyz;
  abcd[42] = 2.0E0*I_KINETIC_G3xz_Fxyz_a-2*I_KINETIC_Dxz_Fxyz;
  abcd[43] = 2.0E0*I_KINETIC_G2x2y_Fxyz_a-1*I_KINETIC_D2y_Fxyz;
  abcd[44] = 2.0E0*I_KINETIC_G2xyz_Fxyz_a-1*I_KINETIC_Dyz_Fxyz;
  abcd[45] = 2.0E0*I_KINETIC_G2x2z_Fxyz_a-1*I_KINETIC_D2z_Fxyz;
  abcd[46] = 2.0E0*I_KINETIC_Gx3y_Fxyz_a;
  abcd[47] = 2.0E0*I_KINETIC_Gx2yz_Fxyz_a;
  abcd[48] = 2.0E0*I_KINETIC_Gxy2z_Fxyz_a;
  abcd[49] = 2.0E0*I_KINETIC_Gx3z_Fxyz_a;
  abcd[50] = 2.0E0*I_KINETIC_G4x_Fx2z_a-3*I_KINETIC_D2x_Fx2z;
  abcd[51] = 2.0E0*I_KINETIC_G3xy_Fx2z_a-2*I_KINETIC_Dxy_Fx2z;
  abcd[52] = 2.0E0*I_KINETIC_G3xz_Fx2z_a-2*I_KINETIC_Dxz_Fx2z;
  abcd[53] = 2.0E0*I_KINETIC_G2x2y_Fx2z_a-1*I_KINETIC_D2y_Fx2z;
  abcd[54] = 2.0E0*I_KINETIC_G2xyz_Fx2z_a-1*I_KINETIC_Dyz_Fx2z;
  abcd[55] = 2.0E0*I_KINETIC_G2x2z_Fx2z_a-1*I_KINETIC_D2z_Fx2z;
  abcd[56] = 2.0E0*I_KINETIC_Gx3y_Fx2z_a;
  abcd[57] = 2.0E0*I_KINETIC_Gx2yz_Fx2z_a;
  abcd[58] = 2.0E0*I_KINETIC_Gxy2z_Fx2z_a;
  abcd[59] = 2.0E0*I_KINETIC_Gx3z_Fx2z_a;
  abcd[60] = 2.0E0*I_KINETIC_G4x_F3y_a-3*I_KINETIC_D2x_F3y;
  abcd[61] = 2.0E0*I_KINETIC_G3xy_F3y_a-2*I_KINETIC_Dxy_F3y;
  abcd[62] = 2.0E0*I_KINETIC_G3xz_F3y_a-2*I_KINETIC_Dxz_F3y;
  abcd[63] = 2.0E0*I_KINETIC_G2x2y_F3y_a-1*I_KINETIC_D2y_F3y;
  abcd[64] = 2.0E0*I_KINETIC_G2xyz_F3y_a-1*I_KINETIC_Dyz_F3y;
  abcd[65] = 2.0E0*I_KINETIC_G2x2z_F3y_a-1*I_KINETIC_D2z_F3y;
  abcd[66] = 2.0E0*I_KINETIC_Gx3y_F3y_a;
  abcd[67] = 2.0E0*I_KINETIC_Gx2yz_F3y_a;
  abcd[68] = 2.0E0*I_KINETIC_Gxy2z_F3y_a;
  abcd[69] = 2.0E0*I_KINETIC_Gx3z_F3y_a;
  abcd[70] = 2.0E0*I_KINETIC_G4x_F2yz_a-3*I_KINETIC_D2x_F2yz;
  abcd[71] = 2.0E0*I_KINETIC_G3xy_F2yz_a-2*I_KINETIC_Dxy_F2yz;
  abcd[72] = 2.0E0*I_KINETIC_G3xz_F2yz_a-2*I_KINETIC_Dxz_F2yz;
  abcd[73] = 2.0E0*I_KINETIC_G2x2y_F2yz_a-1*I_KINETIC_D2y_F2yz;
  abcd[74] = 2.0E0*I_KINETIC_G2xyz_F2yz_a-1*I_KINETIC_Dyz_F2yz;
  abcd[75] = 2.0E0*I_KINETIC_G2x2z_F2yz_a-1*I_KINETIC_D2z_F2yz;
  abcd[76] = 2.0E0*I_KINETIC_Gx3y_F2yz_a;
  abcd[77] = 2.0E0*I_KINETIC_Gx2yz_F2yz_a;
  abcd[78] = 2.0E0*I_KINETIC_Gxy2z_F2yz_a;
  abcd[79] = 2.0E0*I_KINETIC_Gx3z_F2yz_a;
  abcd[80] = 2.0E0*I_KINETIC_G4x_Fy2z_a-3*I_KINETIC_D2x_Fy2z;
  abcd[81] = 2.0E0*I_KINETIC_G3xy_Fy2z_a-2*I_KINETIC_Dxy_Fy2z;
  abcd[82] = 2.0E0*I_KINETIC_G3xz_Fy2z_a-2*I_KINETIC_Dxz_Fy2z;
  abcd[83] = 2.0E0*I_KINETIC_G2x2y_Fy2z_a-1*I_KINETIC_D2y_Fy2z;
  abcd[84] = 2.0E0*I_KINETIC_G2xyz_Fy2z_a-1*I_KINETIC_Dyz_Fy2z;
  abcd[85] = 2.0E0*I_KINETIC_G2x2z_Fy2z_a-1*I_KINETIC_D2z_Fy2z;
  abcd[86] = 2.0E0*I_KINETIC_Gx3y_Fy2z_a;
  abcd[87] = 2.0E0*I_KINETIC_Gx2yz_Fy2z_a;
  abcd[88] = 2.0E0*I_KINETIC_Gxy2z_Fy2z_a;
  abcd[89] = 2.0E0*I_KINETIC_Gx3z_Fy2z_a;
  abcd[90] = 2.0E0*I_KINETIC_G4x_F3z_a-3*I_KINETIC_D2x_F3z;
  abcd[91] = 2.0E0*I_KINETIC_G3xy_F3z_a-2*I_KINETIC_Dxy_F3z;
  abcd[92] = 2.0E0*I_KINETIC_G3xz_F3z_a-2*I_KINETIC_Dxz_F3z;
  abcd[93] = 2.0E0*I_KINETIC_G2x2y_F3z_a-1*I_KINETIC_D2y_F3z;
  abcd[94] = 2.0E0*I_KINETIC_G2xyz_F3z_a-1*I_KINETIC_Dyz_F3z;
  abcd[95] = 2.0E0*I_KINETIC_G2x2z_F3z_a-1*I_KINETIC_D2z_F3z;
  abcd[96] = 2.0E0*I_KINETIC_Gx3y_F3z_a;
  abcd[97] = 2.0E0*I_KINETIC_Gx2yz_F3z_a;
  abcd[98] = 2.0E0*I_KINETIC_Gxy2z_F3z_a;
  abcd[99] = 2.0E0*I_KINETIC_Gx3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_F_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_G_F_a
   * RHS shell quartet name: SQ_KINETIC_D_F
   ************************************************************/
  abcd[100] = 2.0E0*I_KINETIC_G3xy_F3x_a;
  abcd[101] = 2.0E0*I_KINETIC_G2x2y_F3x_a-1*I_KINETIC_D2x_F3x;
  abcd[102] = 2.0E0*I_KINETIC_G2xyz_F3x_a;
  abcd[103] = 2.0E0*I_KINETIC_Gx3y_F3x_a-2*I_KINETIC_Dxy_F3x;
  abcd[104] = 2.0E0*I_KINETIC_Gx2yz_F3x_a-1*I_KINETIC_Dxz_F3x;
  abcd[105] = 2.0E0*I_KINETIC_Gxy2z_F3x_a;
  abcd[106] = 2.0E0*I_KINETIC_G4y_F3x_a-3*I_KINETIC_D2y_F3x;
  abcd[107] = 2.0E0*I_KINETIC_G3yz_F3x_a-2*I_KINETIC_Dyz_F3x;
  abcd[108] = 2.0E0*I_KINETIC_G2y2z_F3x_a-1*I_KINETIC_D2z_F3x;
  abcd[109] = 2.0E0*I_KINETIC_Gy3z_F3x_a;
  abcd[110] = 2.0E0*I_KINETIC_G3xy_F2xy_a;
  abcd[111] = 2.0E0*I_KINETIC_G2x2y_F2xy_a-1*I_KINETIC_D2x_F2xy;
  abcd[112] = 2.0E0*I_KINETIC_G2xyz_F2xy_a;
  abcd[113] = 2.0E0*I_KINETIC_Gx3y_F2xy_a-2*I_KINETIC_Dxy_F2xy;
  abcd[114] = 2.0E0*I_KINETIC_Gx2yz_F2xy_a-1*I_KINETIC_Dxz_F2xy;
  abcd[115] = 2.0E0*I_KINETIC_Gxy2z_F2xy_a;
  abcd[116] = 2.0E0*I_KINETIC_G4y_F2xy_a-3*I_KINETIC_D2y_F2xy;
  abcd[117] = 2.0E0*I_KINETIC_G3yz_F2xy_a-2*I_KINETIC_Dyz_F2xy;
  abcd[118] = 2.0E0*I_KINETIC_G2y2z_F2xy_a-1*I_KINETIC_D2z_F2xy;
  abcd[119] = 2.0E0*I_KINETIC_Gy3z_F2xy_a;
  abcd[120] = 2.0E0*I_KINETIC_G3xy_F2xz_a;
  abcd[121] = 2.0E0*I_KINETIC_G2x2y_F2xz_a-1*I_KINETIC_D2x_F2xz;
  abcd[122] = 2.0E0*I_KINETIC_G2xyz_F2xz_a;
  abcd[123] = 2.0E0*I_KINETIC_Gx3y_F2xz_a-2*I_KINETIC_Dxy_F2xz;
  abcd[124] = 2.0E0*I_KINETIC_Gx2yz_F2xz_a-1*I_KINETIC_Dxz_F2xz;
  abcd[125] = 2.0E0*I_KINETIC_Gxy2z_F2xz_a;
  abcd[126] = 2.0E0*I_KINETIC_G4y_F2xz_a-3*I_KINETIC_D2y_F2xz;
  abcd[127] = 2.0E0*I_KINETIC_G3yz_F2xz_a-2*I_KINETIC_Dyz_F2xz;
  abcd[128] = 2.0E0*I_KINETIC_G2y2z_F2xz_a-1*I_KINETIC_D2z_F2xz;
  abcd[129] = 2.0E0*I_KINETIC_Gy3z_F2xz_a;
  abcd[130] = 2.0E0*I_KINETIC_G3xy_Fx2y_a;
  abcd[131] = 2.0E0*I_KINETIC_G2x2y_Fx2y_a-1*I_KINETIC_D2x_Fx2y;
  abcd[132] = 2.0E0*I_KINETIC_G2xyz_Fx2y_a;
  abcd[133] = 2.0E0*I_KINETIC_Gx3y_Fx2y_a-2*I_KINETIC_Dxy_Fx2y;
  abcd[134] = 2.0E0*I_KINETIC_Gx2yz_Fx2y_a-1*I_KINETIC_Dxz_Fx2y;
  abcd[135] = 2.0E0*I_KINETIC_Gxy2z_Fx2y_a;
  abcd[136] = 2.0E0*I_KINETIC_G4y_Fx2y_a-3*I_KINETIC_D2y_Fx2y;
  abcd[137] = 2.0E0*I_KINETIC_G3yz_Fx2y_a-2*I_KINETIC_Dyz_Fx2y;
  abcd[138] = 2.0E0*I_KINETIC_G2y2z_Fx2y_a-1*I_KINETIC_D2z_Fx2y;
  abcd[139] = 2.0E0*I_KINETIC_Gy3z_Fx2y_a;
  abcd[140] = 2.0E0*I_KINETIC_G3xy_Fxyz_a;
  abcd[141] = 2.0E0*I_KINETIC_G2x2y_Fxyz_a-1*I_KINETIC_D2x_Fxyz;
  abcd[142] = 2.0E0*I_KINETIC_G2xyz_Fxyz_a;
  abcd[143] = 2.0E0*I_KINETIC_Gx3y_Fxyz_a-2*I_KINETIC_Dxy_Fxyz;
  abcd[144] = 2.0E0*I_KINETIC_Gx2yz_Fxyz_a-1*I_KINETIC_Dxz_Fxyz;
  abcd[145] = 2.0E0*I_KINETIC_Gxy2z_Fxyz_a;
  abcd[146] = 2.0E0*I_KINETIC_G4y_Fxyz_a-3*I_KINETIC_D2y_Fxyz;
  abcd[147] = 2.0E0*I_KINETIC_G3yz_Fxyz_a-2*I_KINETIC_Dyz_Fxyz;
  abcd[148] = 2.0E0*I_KINETIC_G2y2z_Fxyz_a-1*I_KINETIC_D2z_Fxyz;
  abcd[149] = 2.0E0*I_KINETIC_Gy3z_Fxyz_a;
  abcd[150] = 2.0E0*I_KINETIC_G3xy_Fx2z_a;
  abcd[151] = 2.0E0*I_KINETIC_G2x2y_Fx2z_a-1*I_KINETIC_D2x_Fx2z;
  abcd[152] = 2.0E0*I_KINETIC_G2xyz_Fx2z_a;
  abcd[153] = 2.0E0*I_KINETIC_Gx3y_Fx2z_a-2*I_KINETIC_Dxy_Fx2z;
  abcd[154] = 2.0E0*I_KINETIC_Gx2yz_Fx2z_a-1*I_KINETIC_Dxz_Fx2z;
  abcd[155] = 2.0E0*I_KINETIC_Gxy2z_Fx2z_a;
  abcd[156] = 2.0E0*I_KINETIC_G4y_Fx2z_a-3*I_KINETIC_D2y_Fx2z;
  abcd[157] = 2.0E0*I_KINETIC_G3yz_Fx2z_a-2*I_KINETIC_Dyz_Fx2z;
  abcd[158] = 2.0E0*I_KINETIC_G2y2z_Fx2z_a-1*I_KINETIC_D2z_Fx2z;
  abcd[159] = 2.0E0*I_KINETIC_Gy3z_Fx2z_a;
  abcd[160] = 2.0E0*I_KINETIC_G3xy_F3y_a;
  abcd[161] = 2.0E0*I_KINETIC_G2x2y_F3y_a-1*I_KINETIC_D2x_F3y;
  abcd[162] = 2.0E0*I_KINETIC_G2xyz_F3y_a;
  abcd[163] = 2.0E0*I_KINETIC_Gx3y_F3y_a-2*I_KINETIC_Dxy_F3y;
  abcd[164] = 2.0E0*I_KINETIC_Gx2yz_F3y_a-1*I_KINETIC_Dxz_F3y;
  abcd[165] = 2.0E0*I_KINETIC_Gxy2z_F3y_a;
  abcd[166] = 2.0E0*I_KINETIC_G4y_F3y_a-3*I_KINETIC_D2y_F3y;
  abcd[167] = 2.0E0*I_KINETIC_G3yz_F3y_a-2*I_KINETIC_Dyz_F3y;
  abcd[168] = 2.0E0*I_KINETIC_G2y2z_F3y_a-1*I_KINETIC_D2z_F3y;
  abcd[169] = 2.0E0*I_KINETIC_Gy3z_F3y_a;
  abcd[170] = 2.0E0*I_KINETIC_G3xy_F2yz_a;
  abcd[171] = 2.0E0*I_KINETIC_G2x2y_F2yz_a-1*I_KINETIC_D2x_F2yz;
  abcd[172] = 2.0E0*I_KINETIC_G2xyz_F2yz_a;
  abcd[173] = 2.0E0*I_KINETIC_Gx3y_F2yz_a-2*I_KINETIC_Dxy_F2yz;
  abcd[174] = 2.0E0*I_KINETIC_Gx2yz_F2yz_a-1*I_KINETIC_Dxz_F2yz;
  abcd[175] = 2.0E0*I_KINETIC_Gxy2z_F2yz_a;
  abcd[176] = 2.0E0*I_KINETIC_G4y_F2yz_a-3*I_KINETIC_D2y_F2yz;
  abcd[177] = 2.0E0*I_KINETIC_G3yz_F2yz_a-2*I_KINETIC_Dyz_F2yz;
  abcd[178] = 2.0E0*I_KINETIC_G2y2z_F2yz_a-1*I_KINETIC_D2z_F2yz;
  abcd[179] = 2.0E0*I_KINETIC_Gy3z_F2yz_a;
  abcd[180] = 2.0E0*I_KINETIC_G3xy_Fy2z_a;
  abcd[181] = 2.0E0*I_KINETIC_G2x2y_Fy2z_a-1*I_KINETIC_D2x_Fy2z;
  abcd[182] = 2.0E0*I_KINETIC_G2xyz_Fy2z_a;
  abcd[183] = 2.0E0*I_KINETIC_Gx3y_Fy2z_a-2*I_KINETIC_Dxy_Fy2z;
  abcd[184] = 2.0E0*I_KINETIC_Gx2yz_Fy2z_a-1*I_KINETIC_Dxz_Fy2z;
  abcd[185] = 2.0E0*I_KINETIC_Gxy2z_Fy2z_a;
  abcd[186] = 2.0E0*I_KINETIC_G4y_Fy2z_a-3*I_KINETIC_D2y_Fy2z;
  abcd[187] = 2.0E0*I_KINETIC_G3yz_Fy2z_a-2*I_KINETIC_Dyz_Fy2z;
  abcd[188] = 2.0E0*I_KINETIC_G2y2z_Fy2z_a-1*I_KINETIC_D2z_Fy2z;
  abcd[189] = 2.0E0*I_KINETIC_Gy3z_Fy2z_a;
  abcd[190] = 2.0E0*I_KINETIC_G3xy_F3z_a;
  abcd[191] = 2.0E0*I_KINETIC_G2x2y_F3z_a-1*I_KINETIC_D2x_F3z;
  abcd[192] = 2.0E0*I_KINETIC_G2xyz_F3z_a;
  abcd[193] = 2.0E0*I_KINETIC_Gx3y_F3z_a-2*I_KINETIC_Dxy_F3z;
  abcd[194] = 2.0E0*I_KINETIC_Gx2yz_F3z_a-1*I_KINETIC_Dxz_F3z;
  abcd[195] = 2.0E0*I_KINETIC_Gxy2z_F3z_a;
  abcd[196] = 2.0E0*I_KINETIC_G4y_F3z_a-3*I_KINETIC_D2y_F3z;
  abcd[197] = 2.0E0*I_KINETIC_G3yz_F3z_a-2*I_KINETIC_Dyz_F3z;
  abcd[198] = 2.0E0*I_KINETIC_G2y2z_F3z_a-1*I_KINETIC_D2z_F3z;
  abcd[199] = 2.0E0*I_KINETIC_Gy3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_F_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_G_F_a
   * RHS shell quartet name: SQ_KINETIC_D_F
   ************************************************************/
  abcd[200] = 2.0E0*I_KINETIC_G3xz_F3x_a;
  abcd[201] = 2.0E0*I_KINETIC_G2xyz_F3x_a;
  abcd[202] = 2.0E0*I_KINETIC_G2x2z_F3x_a-1*I_KINETIC_D2x_F3x;
  abcd[203] = 2.0E0*I_KINETIC_Gx2yz_F3x_a;
  abcd[204] = 2.0E0*I_KINETIC_Gxy2z_F3x_a-1*I_KINETIC_Dxy_F3x;
  abcd[205] = 2.0E0*I_KINETIC_Gx3z_F3x_a-2*I_KINETIC_Dxz_F3x;
  abcd[206] = 2.0E0*I_KINETIC_G3yz_F3x_a;
  abcd[207] = 2.0E0*I_KINETIC_G2y2z_F3x_a-1*I_KINETIC_D2y_F3x;
  abcd[208] = 2.0E0*I_KINETIC_Gy3z_F3x_a-2*I_KINETIC_Dyz_F3x;
  abcd[209] = 2.0E0*I_KINETIC_G4z_F3x_a-3*I_KINETIC_D2z_F3x;
  abcd[210] = 2.0E0*I_KINETIC_G3xz_F2xy_a;
  abcd[211] = 2.0E0*I_KINETIC_G2xyz_F2xy_a;
  abcd[212] = 2.0E0*I_KINETIC_G2x2z_F2xy_a-1*I_KINETIC_D2x_F2xy;
  abcd[213] = 2.0E0*I_KINETIC_Gx2yz_F2xy_a;
  abcd[214] = 2.0E0*I_KINETIC_Gxy2z_F2xy_a-1*I_KINETIC_Dxy_F2xy;
  abcd[215] = 2.0E0*I_KINETIC_Gx3z_F2xy_a-2*I_KINETIC_Dxz_F2xy;
  abcd[216] = 2.0E0*I_KINETIC_G3yz_F2xy_a;
  abcd[217] = 2.0E0*I_KINETIC_G2y2z_F2xy_a-1*I_KINETIC_D2y_F2xy;
  abcd[218] = 2.0E0*I_KINETIC_Gy3z_F2xy_a-2*I_KINETIC_Dyz_F2xy;
  abcd[219] = 2.0E0*I_KINETIC_G4z_F2xy_a-3*I_KINETIC_D2z_F2xy;
  abcd[220] = 2.0E0*I_KINETIC_G3xz_F2xz_a;
  abcd[221] = 2.0E0*I_KINETIC_G2xyz_F2xz_a;
  abcd[222] = 2.0E0*I_KINETIC_G2x2z_F2xz_a-1*I_KINETIC_D2x_F2xz;
  abcd[223] = 2.0E0*I_KINETIC_Gx2yz_F2xz_a;
  abcd[224] = 2.0E0*I_KINETIC_Gxy2z_F2xz_a-1*I_KINETIC_Dxy_F2xz;
  abcd[225] = 2.0E0*I_KINETIC_Gx3z_F2xz_a-2*I_KINETIC_Dxz_F2xz;
  abcd[226] = 2.0E0*I_KINETIC_G3yz_F2xz_a;
  abcd[227] = 2.0E0*I_KINETIC_G2y2z_F2xz_a-1*I_KINETIC_D2y_F2xz;
  abcd[228] = 2.0E0*I_KINETIC_Gy3z_F2xz_a-2*I_KINETIC_Dyz_F2xz;
  abcd[229] = 2.0E0*I_KINETIC_G4z_F2xz_a-3*I_KINETIC_D2z_F2xz;
  abcd[230] = 2.0E0*I_KINETIC_G3xz_Fx2y_a;
  abcd[231] = 2.0E0*I_KINETIC_G2xyz_Fx2y_a;
  abcd[232] = 2.0E0*I_KINETIC_G2x2z_Fx2y_a-1*I_KINETIC_D2x_Fx2y;
  abcd[233] = 2.0E0*I_KINETIC_Gx2yz_Fx2y_a;
  abcd[234] = 2.0E0*I_KINETIC_Gxy2z_Fx2y_a-1*I_KINETIC_Dxy_Fx2y;
  abcd[235] = 2.0E0*I_KINETIC_Gx3z_Fx2y_a-2*I_KINETIC_Dxz_Fx2y;
  abcd[236] = 2.0E0*I_KINETIC_G3yz_Fx2y_a;
  abcd[237] = 2.0E0*I_KINETIC_G2y2z_Fx2y_a-1*I_KINETIC_D2y_Fx2y;
  abcd[238] = 2.0E0*I_KINETIC_Gy3z_Fx2y_a-2*I_KINETIC_Dyz_Fx2y;
  abcd[239] = 2.0E0*I_KINETIC_G4z_Fx2y_a-3*I_KINETIC_D2z_Fx2y;
  abcd[240] = 2.0E0*I_KINETIC_G3xz_Fxyz_a;
  abcd[241] = 2.0E0*I_KINETIC_G2xyz_Fxyz_a;
  abcd[242] = 2.0E0*I_KINETIC_G2x2z_Fxyz_a-1*I_KINETIC_D2x_Fxyz;
  abcd[243] = 2.0E0*I_KINETIC_Gx2yz_Fxyz_a;
  abcd[244] = 2.0E0*I_KINETIC_Gxy2z_Fxyz_a-1*I_KINETIC_Dxy_Fxyz;
  abcd[245] = 2.0E0*I_KINETIC_Gx3z_Fxyz_a-2*I_KINETIC_Dxz_Fxyz;
  abcd[246] = 2.0E0*I_KINETIC_G3yz_Fxyz_a;
  abcd[247] = 2.0E0*I_KINETIC_G2y2z_Fxyz_a-1*I_KINETIC_D2y_Fxyz;
  abcd[248] = 2.0E0*I_KINETIC_Gy3z_Fxyz_a-2*I_KINETIC_Dyz_Fxyz;
  abcd[249] = 2.0E0*I_KINETIC_G4z_Fxyz_a-3*I_KINETIC_D2z_Fxyz;
  abcd[250] = 2.0E0*I_KINETIC_G3xz_Fx2z_a;
  abcd[251] = 2.0E0*I_KINETIC_G2xyz_Fx2z_a;
  abcd[252] = 2.0E0*I_KINETIC_G2x2z_Fx2z_a-1*I_KINETIC_D2x_Fx2z;
  abcd[253] = 2.0E0*I_KINETIC_Gx2yz_Fx2z_a;
  abcd[254] = 2.0E0*I_KINETIC_Gxy2z_Fx2z_a-1*I_KINETIC_Dxy_Fx2z;
  abcd[255] = 2.0E0*I_KINETIC_Gx3z_Fx2z_a-2*I_KINETIC_Dxz_Fx2z;
  abcd[256] = 2.0E0*I_KINETIC_G3yz_Fx2z_a;
  abcd[257] = 2.0E0*I_KINETIC_G2y2z_Fx2z_a-1*I_KINETIC_D2y_Fx2z;
  abcd[258] = 2.0E0*I_KINETIC_Gy3z_Fx2z_a-2*I_KINETIC_Dyz_Fx2z;
  abcd[259] = 2.0E0*I_KINETIC_G4z_Fx2z_a-3*I_KINETIC_D2z_Fx2z;
  abcd[260] = 2.0E0*I_KINETIC_G3xz_F3y_a;
  abcd[261] = 2.0E0*I_KINETIC_G2xyz_F3y_a;
  abcd[262] = 2.0E0*I_KINETIC_G2x2z_F3y_a-1*I_KINETIC_D2x_F3y;
  abcd[263] = 2.0E0*I_KINETIC_Gx2yz_F3y_a;
  abcd[264] = 2.0E0*I_KINETIC_Gxy2z_F3y_a-1*I_KINETIC_Dxy_F3y;
  abcd[265] = 2.0E0*I_KINETIC_Gx3z_F3y_a-2*I_KINETIC_Dxz_F3y;
  abcd[266] = 2.0E0*I_KINETIC_G3yz_F3y_a;
  abcd[267] = 2.0E0*I_KINETIC_G2y2z_F3y_a-1*I_KINETIC_D2y_F3y;
  abcd[268] = 2.0E0*I_KINETIC_Gy3z_F3y_a-2*I_KINETIC_Dyz_F3y;
  abcd[269] = 2.0E0*I_KINETIC_G4z_F3y_a-3*I_KINETIC_D2z_F3y;
  abcd[270] = 2.0E0*I_KINETIC_G3xz_F2yz_a;
  abcd[271] = 2.0E0*I_KINETIC_G2xyz_F2yz_a;
  abcd[272] = 2.0E0*I_KINETIC_G2x2z_F2yz_a-1*I_KINETIC_D2x_F2yz;
  abcd[273] = 2.0E0*I_KINETIC_Gx2yz_F2yz_a;
  abcd[274] = 2.0E0*I_KINETIC_Gxy2z_F2yz_a-1*I_KINETIC_Dxy_F2yz;
  abcd[275] = 2.0E0*I_KINETIC_Gx3z_F2yz_a-2*I_KINETIC_Dxz_F2yz;
  abcd[276] = 2.0E0*I_KINETIC_G3yz_F2yz_a;
  abcd[277] = 2.0E0*I_KINETIC_G2y2z_F2yz_a-1*I_KINETIC_D2y_F2yz;
  abcd[278] = 2.0E0*I_KINETIC_Gy3z_F2yz_a-2*I_KINETIC_Dyz_F2yz;
  abcd[279] = 2.0E0*I_KINETIC_G4z_F2yz_a-3*I_KINETIC_D2z_F2yz;
  abcd[280] = 2.0E0*I_KINETIC_G3xz_Fy2z_a;
  abcd[281] = 2.0E0*I_KINETIC_G2xyz_Fy2z_a;
  abcd[282] = 2.0E0*I_KINETIC_G2x2z_Fy2z_a-1*I_KINETIC_D2x_Fy2z;
  abcd[283] = 2.0E0*I_KINETIC_Gx2yz_Fy2z_a;
  abcd[284] = 2.0E0*I_KINETIC_Gxy2z_Fy2z_a-1*I_KINETIC_Dxy_Fy2z;
  abcd[285] = 2.0E0*I_KINETIC_Gx3z_Fy2z_a-2*I_KINETIC_Dxz_Fy2z;
  abcd[286] = 2.0E0*I_KINETIC_G3yz_Fy2z_a;
  abcd[287] = 2.0E0*I_KINETIC_G2y2z_Fy2z_a-1*I_KINETIC_D2y_Fy2z;
  abcd[288] = 2.0E0*I_KINETIC_Gy3z_Fy2z_a-2*I_KINETIC_Dyz_Fy2z;
  abcd[289] = 2.0E0*I_KINETIC_G4z_Fy2z_a-3*I_KINETIC_D2z_Fy2z;
  abcd[290] = 2.0E0*I_KINETIC_G3xz_F3z_a;
  abcd[291] = 2.0E0*I_KINETIC_G2xyz_F3z_a;
  abcd[292] = 2.0E0*I_KINETIC_G2x2z_F3z_a-1*I_KINETIC_D2x_F3z;
  abcd[293] = 2.0E0*I_KINETIC_Gx2yz_F3z_a;
  abcd[294] = 2.0E0*I_KINETIC_Gxy2z_F3z_a-1*I_KINETIC_Dxy_F3z;
  abcd[295] = 2.0E0*I_KINETIC_Gx3z_F3z_a-2*I_KINETIC_Dxz_F3z;
  abcd[296] = 2.0E0*I_KINETIC_G3yz_F3z_a;
  abcd[297] = 2.0E0*I_KINETIC_G2y2z_F3z_a-1*I_KINETIC_D2y_F3z;
  abcd[298] = 2.0E0*I_KINETIC_Gy3z_F3z_a-2*I_KINETIC_Dyz_F3z;
  abcd[299] = 2.0E0*I_KINETIC_G4z_F3z_a-3*I_KINETIC_D2z_F3z;
}
