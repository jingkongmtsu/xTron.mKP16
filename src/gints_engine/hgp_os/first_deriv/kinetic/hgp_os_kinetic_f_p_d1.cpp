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
// BRA1 as redundant position, total RHS integrals evaluated as: 1965
// BRA2 as redundant position, total RHS integrals evaluated as: 1359
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

void hgp_os_kinetic_f_p_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_G4x_Px_a = 0.0E0;
  Double I_KINETIC_G3xy_Px_a = 0.0E0;
  Double I_KINETIC_G3xz_Px_a = 0.0E0;
  Double I_KINETIC_G2x2y_Px_a = 0.0E0;
  Double I_KINETIC_G2xyz_Px_a = 0.0E0;
  Double I_KINETIC_G2x2z_Px_a = 0.0E0;
  Double I_KINETIC_Gx3y_Px_a = 0.0E0;
  Double I_KINETIC_Gx2yz_Px_a = 0.0E0;
  Double I_KINETIC_Gxy2z_Px_a = 0.0E0;
  Double I_KINETIC_Gx3z_Px_a = 0.0E0;
  Double I_KINETIC_G4y_Px_a = 0.0E0;
  Double I_KINETIC_G3yz_Px_a = 0.0E0;
  Double I_KINETIC_G2y2z_Px_a = 0.0E0;
  Double I_KINETIC_Gy3z_Px_a = 0.0E0;
  Double I_KINETIC_G4z_Px_a = 0.0E0;
  Double I_KINETIC_G4x_Py_a = 0.0E0;
  Double I_KINETIC_G3xy_Py_a = 0.0E0;
  Double I_KINETIC_G3xz_Py_a = 0.0E0;
  Double I_KINETIC_G2x2y_Py_a = 0.0E0;
  Double I_KINETIC_G2xyz_Py_a = 0.0E0;
  Double I_KINETIC_G2x2z_Py_a = 0.0E0;
  Double I_KINETIC_Gx3y_Py_a = 0.0E0;
  Double I_KINETIC_Gx2yz_Py_a = 0.0E0;
  Double I_KINETIC_Gxy2z_Py_a = 0.0E0;
  Double I_KINETIC_Gx3z_Py_a = 0.0E0;
  Double I_KINETIC_G4y_Py_a = 0.0E0;
  Double I_KINETIC_G3yz_Py_a = 0.0E0;
  Double I_KINETIC_G2y2z_Py_a = 0.0E0;
  Double I_KINETIC_Gy3z_Py_a = 0.0E0;
  Double I_KINETIC_G4z_Py_a = 0.0E0;
  Double I_KINETIC_G4x_Pz_a = 0.0E0;
  Double I_KINETIC_G3xy_Pz_a = 0.0E0;
  Double I_KINETIC_G3xz_Pz_a = 0.0E0;
  Double I_KINETIC_G2x2y_Pz_a = 0.0E0;
  Double I_KINETIC_G2xyz_Pz_a = 0.0E0;
  Double I_KINETIC_G2x2z_Pz_a = 0.0E0;
  Double I_KINETIC_Gx3y_Pz_a = 0.0E0;
  Double I_KINETIC_Gx2yz_Pz_a = 0.0E0;
  Double I_KINETIC_Gxy2z_Pz_a = 0.0E0;
  Double I_KINETIC_Gx3z_Pz_a = 0.0E0;
  Double I_KINETIC_G4y_Pz_a = 0.0E0;
  Double I_KINETIC_G3yz_Pz_a = 0.0E0;
  Double I_KINETIC_G2y2z_Pz_a = 0.0E0;
  Double I_KINETIC_Gy3z_Pz_a = 0.0E0;
  Double I_KINETIC_G4z_Pz_a = 0.0E0;
  Double I_KINETIC_D2x_Px = 0.0E0;
  Double I_KINETIC_Dxy_Px = 0.0E0;
  Double I_KINETIC_Dxz_Px = 0.0E0;
  Double I_KINETIC_D2y_Px = 0.0E0;
  Double I_KINETIC_Dyz_Px = 0.0E0;
  Double I_KINETIC_D2z_Px = 0.0E0;
  Double I_KINETIC_D2x_Py = 0.0E0;
  Double I_KINETIC_Dxy_Py = 0.0E0;
  Double I_KINETIC_Dxz_Py = 0.0E0;
  Double I_KINETIC_D2y_Py = 0.0E0;
  Double I_KINETIC_Dyz_Py = 0.0E0;
  Double I_KINETIC_D2z_Py = 0.0E0;
  Double I_KINETIC_D2x_Pz = 0.0E0;
  Double I_KINETIC_Dxy_Pz = 0.0E0;
  Double I_KINETIC_Dxz_Pz = 0.0E0;
  Double I_KINETIC_D2y_Pz = 0.0E0;
  Double I_KINETIC_Dyz_Pz = 0.0E0;
  Double I_KINETIC_D2z_Pz = 0.0E0;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_Px_vrr = PBX*I_TWOBODYOVERLAP_G4x_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Px_vrr = PBX*I_TWOBODYOVERLAP_G3xy_S_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Px_vrr = PBX*I_TWOBODYOVERLAP_G3xz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Px_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Px_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Px_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Px_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Px_vrr = PBX*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Px_vrr = PBX*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Px_vrr = PBX*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Py_vrr = PBY*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Py_vrr = PBY*I_TWOBODYOVERLAP_G3xy_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Py_vrr = PBY*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Py_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Py_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Py_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Py_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Py_vrr = PBY*I_TWOBODYOVERLAP_G4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Py_vrr = PBY*I_TWOBODYOVERLAP_G3yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Py_vrr = PBY*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

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
     * shell quartet name: SQ_KINETIC_G_P_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_G_P_a_coefs = alpha;
    I_KINETIC_G4x_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G4x_Px_vrr;
    I_KINETIC_G3xy_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G3xy_Px_vrr;
    I_KINETIC_G3xz_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G3xz_Px_vrr;
    I_KINETIC_G2x2y_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2x2y_Px_vrr;
    I_KINETIC_G2xyz_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2xyz_Px_vrr;
    I_KINETIC_G2x2z_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2x2z_Px_vrr;
    I_KINETIC_Gx3y_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gx3y_Px_vrr;
    I_KINETIC_Gx2yz_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gx2yz_Px_vrr;
    I_KINETIC_Gxy2z_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gxy2z_Px_vrr;
    I_KINETIC_Gx3z_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gx3z_Px_vrr;
    I_KINETIC_G4y_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G4y_Px_vrr;
    I_KINETIC_G3yz_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G3yz_Px_vrr;
    I_KINETIC_G2y2z_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2y2z_Px_vrr;
    I_KINETIC_Gy3z_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gy3z_Px_vrr;
    I_KINETIC_G4z_Px_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G4z_Px_vrr;
    I_KINETIC_G4x_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G4x_Py_vrr;
    I_KINETIC_G3xy_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G3xy_Py_vrr;
    I_KINETIC_G3xz_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G3xz_Py_vrr;
    I_KINETIC_G2x2y_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2x2y_Py_vrr;
    I_KINETIC_G2xyz_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2xyz_Py_vrr;
    I_KINETIC_G2x2z_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2x2z_Py_vrr;
    I_KINETIC_Gx3y_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gx3y_Py_vrr;
    I_KINETIC_Gx2yz_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gx2yz_Py_vrr;
    I_KINETIC_Gxy2z_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gxy2z_Py_vrr;
    I_KINETIC_Gx3z_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gx3z_Py_vrr;
    I_KINETIC_G4y_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G4y_Py_vrr;
    I_KINETIC_G3yz_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G3yz_Py_vrr;
    I_KINETIC_G2y2z_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2y2z_Py_vrr;
    I_KINETIC_Gy3z_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gy3z_Py_vrr;
    I_KINETIC_G4z_Py_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G4z_Py_vrr;
    I_KINETIC_G4x_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G4x_Pz_vrr;
    I_KINETIC_G3xy_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G3xy_Pz_vrr;
    I_KINETIC_G3xz_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G3xz_Pz_vrr;
    I_KINETIC_G2x2y_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2x2y_Pz_vrr;
    I_KINETIC_G2xyz_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2xyz_Pz_vrr;
    I_KINETIC_G2x2z_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2x2z_Pz_vrr;
    I_KINETIC_Gx3y_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gx3y_Pz_vrr;
    I_KINETIC_Gx2yz_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gx2yz_Pz_vrr;
    I_KINETIC_Gxy2z_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gxy2z_Pz_vrr;
    I_KINETIC_Gx3z_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gx3z_Pz_vrr;
    I_KINETIC_G4y_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G4y_Pz_vrr;
    I_KINETIC_G3yz_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G3yz_Pz_vrr;
    I_KINETIC_G2y2z_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G2y2z_Pz_vrr;
    I_KINETIC_Gy3z_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_Gy3z_Pz_vrr;
    I_KINETIC_G4z_Pz_a += SQ_KINETIC_G_P_a_coefs*I_KINETIC_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_P
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_D2x_Px += I_KINETIC_D2x_Px_vrr;
    I_KINETIC_Dxy_Px += I_KINETIC_Dxy_Px_vrr;
    I_KINETIC_Dxz_Px += I_KINETIC_Dxz_Px_vrr;
    I_KINETIC_D2y_Px += I_KINETIC_D2y_Px_vrr;
    I_KINETIC_Dyz_Px += I_KINETIC_Dyz_Px_vrr;
    I_KINETIC_D2z_Px += I_KINETIC_D2z_Px_vrr;
    I_KINETIC_D2x_Py += I_KINETIC_D2x_Py_vrr;
    I_KINETIC_Dxy_Py += I_KINETIC_Dxy_Py_vrr;
    I_KINETIC_Dxz_Py += I_KINETIC_Dxz_Py_vrr;
    I_KINETIC_D2y_Py += I_KINETIC_D2y_Py_vrr;
    I_KINETIC_Dyz_Py += I_KINETIC_Dyz_Py_vrr;
    I_KINETIC_D2z_Py += I_KINETIC_D2z_Py_vrr;
    I_KINETIC_D2x_Pz += I_KINETIC_D2x_Pz_vrr;
    I_KINETIC_Dxy_Pz += I_KINETIC_Dxy_Pz_vrr;
    I_KINETIC_Dxz_Pz += I_KINETIC_Dxz_Pz_vrr;
    I_KINETIC_D2y_Pz += I_KINETIC_D2y_Pz_vrr;
    I_KINETIC_Dyz_Pz += I_KINETIC_Dyz_Pz_vrr;
    I_KINETIC_D2z_Pz += I_KINETIC_D2z_Pz_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_P_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_D_P
   ************************************************************/
  abcd[0] = 2.0E0*I_KINETIC_G4x_Px_a-3*I_KINETIC_D2x_Px;
  abcd[1] = 2.0E0*I_KINETIC_G3xy_Px_a-2*I_KINETIC_Dxy_Px;
  abcd[2] = 2.0E0*I_KINETIC_G3xz_Px_a-2*I_KINETIC_Dxz_Px;
  abcd[3] = 2.0E0*I_KINETIC_G2x2y_Px_a-1*I_KINETIC_D2y_Px;
  abcd[4] = 2.0E0*I_KINETIC_G2xyz_Px_a-1*I_KINETIC_Dyz_Px;
  abcd[5] = 2.0E0*I_KINETIC_G2x2z_Px_a-1*I_KINETIC_D2z_Px;
  abcd[6] = 2.0E0*I_KINETIC_Gx3y_Px_a;
  abcd[7] = 2.0E0*I_KINETIC_Gx2yz_Px_a;
  abcd[8] = 2.0E0*I_KINETIC_Gxy2z_Px_a;
  abcd[9] = 2.0E0*I_KINETIC_Gx3z_Px_a;
  abcd[10] = 2.0E0*I_KINETIC_G4x_Py_a-3*I_KINETIC_D2x_Py;
  abcd[11] = 2.0E0*I_KINETIC_G3xy_Py_a-2*I_KINETIC_Dxy_Py;
  abcd[12] = 2.0E0*I_KINETIC_G3xz_Py_a-2*I_KINETIC_Dxz_Py;
  abcd[13] = 2.0E0*I_KINETIC_G2x2y_Py_a-1*I_KINETIC_D2y_Py;
  abcd[14] = 2.0E0*I_KINETIC_G2xyz_Py_a-1*I_KINETIC_Dyz_Py;
  abcd[15] = 2.0E0*I_KINETIC_G2x2z_Py_a-1*I_KINETIC_D2z_Py;
  abcd[16] = 2.0E0*I_KINETIC_Gx3y_Py_a;
  abcd[17] = 2.0E0*I_KINETIC_Gx2yz_Py_a;
  abcd[18] = 2.0E0*I_KINETIC_Gxy2z_Py_a;
  abcd[19] = 2.0E0*I_KINETIC_Gx3z_Py_a;
  abcd[20] = 2.0E0*I_KINETIC_G4x_Pz_a-3*I_KINETIC_D2x_Pz;
  abcd[21] = 2.0E0*I_KINETIC_G3xy_Pz_a-2*I_KINETIC_Dxy_Pz;
  abcd[22] = 2.0E0*I_KINETIC_G3xz_Pz_a-2*I_KINETIC_Dxz_Pz;
  abcd[23] = 2.0E0*I_KINETIC_G2x2y_Pz_a-1*I_KINETIC_D2y_Pz;
  abcd[24] = 2.0E0*I_KINETIC_G2xyz_Pz_a-1*I_KINETIC_Dyz_Pz;
  abcd[25] = 2.0E0*I_KINETIC_G2x2z_Pz_a-1*I_KINETIC_D2z_Pz;
  abcd[26] = 2.0E0*I_KINETIC_Gx3y_Pz_a;
  abcd[27] = 2.0E0*I_KINETIC_Gx2yz_Pz_a;
  abcd[28] = 2.0E0*I_KINETIC_Gxy2z_Pz_a;
  abcd[29] = 2.0E0*I_KINETIC_Gx3z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_P_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_D_P
   ************************************************************/
  abcd[30] = 2.0E0*I_KINETIC_G3xy_Px_a;
  abcd[31] = 2.0E0*I_KINETIC_G2x2y_Px_a-1*I_KINETIC_D2x_Px;
  abcd[32] = 2.0E0*I_KINETIC_G2xyz_Px_a;
  abcd[33] = 2.0E0*I_KINETIC_Gx3y_Px_a-2*I_KINETIC_Dxy_Px;
  abcd[34] = 2.0E0*I_KINETIC_Gx2yz_Px_a-1*I_KINETIC_Dxz_Px;
  abcd[35] = 2.0E0*I_KINETIC_Gxy2z_Px_a;
  abcd[36] = 2.0E0*I_KINETIC_G4y_Px_a-3*I_KINETIC_D2y_Px;
  abcd[37] = 2.0E0*I_KINETIC_G3yz_Px_a-2*I_KINETIC_Dyz_Px;
  abcd[38] = 2.0E0*I_KINETIC_G2y2z_Px_a-1*I_KINETIC_D2z_Px;
  abcd[39] = 2.0E0*I_KINETIC_Gy3z_Px_a;
  abcd[40] = 2.0E0*I_KINETIC_G3xy_Py_a;
  abcd[41] = 2.0E0*I_KINETIC_G2x2y_Py_a-1*I_KINETIC_D2x_Py;
  abcd[42] = 2.0E0*I_KINETIC_G2xyz_Py_a;
  abcd[43] = 2.0E0*I_KINETIC_Gx3y_Py_a-2*I_KINETIC_Dxy_Py;
  abcd[44] = 2.0E0*I_KINETIC_Gx2yz_Py_a-1*I_KINETIC_Dxz_Py;
  abcd[45] = 2.0E0*I_KINETIC_Gxy2z_Py_a;
  abcd[46] = 2.0E0*I_KINETIC_G4y_Py_a-3*I_KINETIC_D2y_Py;
  abcd[47] = 2.0E0*I_KINETIC_G3yz_Py_a-2*I_KINETIC_Dyz_Py;
  abcd[48] = 2.0E0*I_KINETIC_G2y2z_Py_a-1*I_KINETIC_D2z_Py;
  abcd[49] = 2.0E0*I_KINETIC_Gy3z_Py_a;
  abcd[50] = 2.0E0*I_KINETIC_G3xy_Pz_a;
  abcd[51] = 2.0E0*I_KINETIC_G2x2y_Pz_a-1*I_KINETIC_D2x_Pz;
  abcd[52] = 2.0E0*I_KINETIC_G2xyz_Pz_a;
  abcd[53] = 2.0E0*I_KINETIC_Gx3y_Pz_a-2*I_KINETIC_Dxy_Pz;
  abcd[54] = 2.0E0*I_KINETIC_Gx2yz_Pz_a-1*I_KINETIC_Dxz_Pz;
  abcd[55] = 2.0E0*I_KINETIC_Gxy2z_Pz_a;
  abcd[56] = 2.0E0*I_KINETIC_G4y_Pz_a-3*I_KINETIC_D2y_Pz;
  abcd[57] = 2.0E0*I_KINETIC_G3yz_Pz_a-2*I_KINETIC_Dyz_Pz;
  abcd[58] = 2.0E0*I_KINETIC_G2y2z_Pz_a-1*I_KINETIC_D2z_Pz;
  abcd[59] = 2.0E0*I_KINETIC_Gy3z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_F_P_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_D_P
   ************************************************************/
  abcd[60] = 2.0E0*I_KINETIC_G3xz_Px_a;
  abcd[61] = 2.0E0*I_KINETIC_G2xyz_Px_a;
  abcd[62] = 2.0E0*I_KINETIC_G2x2z_Px_a-1*I_KINETIC_D2x_Px;
  abcd[63] = 2.0E0*I_KINETIC_Gx2yz_Px_a;
  abcd[64] = 2.0E0*I_KINETIC_Gxy2z_Px_a-1*I_KINETIC_Dxy_Px;
  abcd[65] = 2.0E0*I_KINETIC_Gx3z_Px_a-2*I_KINETIC_Dxz_Px;
  abcd[66] = 2.0E0*I_KINETIC_G3yz_Px_a;
  abcd[67] = 2.0E0*I_KINETIC_G2y2z_Px_a-1*I_KINETIC_D2y_Px;
  abcd[68] = 2.0E0*I_KINETIC_Gy3z_Px_a-2*I_KINETIC_Dyz_Px;
  abcd[69] = 2.0E0*I_KINETIC_G4z_Px_a-3*I_KINETIC_D2z_Px;
  abcd[70] = 2.0E0*I_KINETIC_G3xz_Py_a;
  abcd[71] = 2.0E0*I_KINETIC_G2xyz_Py_a;
  abcd[72] = 2.0E0*I_KINETIC_G2x2z_Py_a-1*I_KINETIC_D2x_Py;
  abcd[73] = 2.0E0*I_KINETIC_Gx2yz_Py_a;
  abcd[74] = 2.0E0*I_KINETIC_Gxy2z_Py_a-1*I_KINETIC_Dxy_Py;
  abcd[75] = 2.0E0*I_KINETIC_Gx3z_Py_a-2*I_KINETIC_Dxz_Py;
  abcd[76] = 2.0E0*I_KINETIC_G3yz_Py_a;
  abcd[77] = 2.0E0*I_KINETIC_G2y2z_Py_a-1*I_KINETIC_D2y_Py;
  abcd[78] = 2.0E0*I_KINETIC_Gy3z_Py_a-2*I_KINETIC_Dyz_Py;
  abcd[79] = 2.0E0*I_KINETIC_G4z_Py_a-3*I_KINETIC_D2z_Py;
  abcd[80] = 2.0E0*I_KINETIC_G3xz_Pz_a;
  abcd[81] = 2.0E0*I_KINETIC_G2xyz_Pz_a;
  abcd[82] = 2.0E0*I_KINETIC_G2x2z_Pz_a-1*I_KINETIC_D2x_Pz;
  abcd[83] = 2.0E0*I_KINETIC_Gx2yz_Pz_a;
  abcd[84] = 2.0E0*I_KINETIC_Gxy2z_Pz_a-1*I_KINETIC_Dxy_Pz;
  abcd[85] = 2.0E0*I_KINETIC_Gx3z_Pz_a-2*I_KINETIC_Dxz_Pz;
  abcd[86] = 2.0E0*I_KINETIC_G3yz_Pz_a;
  abcd[87] = 2.0E0*I_KINETIC_G2y2z_Pz_a-1*I_KINETIC_D2y_Pz;
  abcd[88] = 2.0E0*I_KINETIC_Gy3z_Pz_a-2*I_KINETIC_Dyz_Pz;
  abcd[89] = 2.0E0*I_KINETIC_G4z_Pz_a-3*I_KINETIC_D2z_Pz;
}
