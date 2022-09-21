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
// BRA1 as redundant position, total RHS integrals evaluated as: 6456
// BRA2 as redundant position, total RHS integrals evaluated as: 3771
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

void hgp_os_kinetic_g_p_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_I6x_Px_aa = 0.0E0;
  Double I_KINETIC_I5xy_Px_aa = 0.0E0;
  Double I_KINETIC_I5xz_Px_aa = 0.0E0;
  Double I_KINETIC_I4x2y_Px_aa = 0.0E0;
  Double I_KINETIC_I4xyz_Px_aa = 0.0E0;
  Double I_KINETIC_I4x2z_Px_aa = 0.0E0;
  Double I_KINETIC_I3x3y_Px_aa = 0.0E0;
  Double I_KINETIC_I3x2yz_Px_aa = 0.0E0;
  Double I_KINETIC_I3xy2z_Px_aa = 0.0E0;
  Double I_KINETIC_I3x3z_Px_aa = 0.0E0;
  Double I_KINETIC_I2x4y_Px_aa = 0.0E0;
  Double I_KINETIC_I2x3yz_Px_aa = 0.0E0;
  Double I_KINETIC_I2x2y2z_Px_aa = 0.0E0;
  Double I_KINETIC_I2xy3z_Px_aa = 0.0E0;
  Double I_KINETIC_I2x4z_Px_aa = 0.0E0;
  Double I_KINETIC_Ix5y_Px_aa = 0.0E0;
  Double I_KINETIC_Ix4yz_Px_aa = 0.0E0;
  Double I_KINETIC_Ix3y2z_Px_aa = 0.0E0;
  Double I_KINETIC_Ix2y3z_Px_aa = 0.0E0;
  Double I_KINETIC_Ixy4z_Px_aa = 0.0E0;
  Double I_KINETIC_Ix5z_Px_aa = 0.0E0;
  Double I_KINETIC_I6y_Px_aa = 0.0E0;
  Double I_KINETIC_I5yz_Px_aa = 0.0E0;
  Double I_KINETIC_I4y2z_Px_aa = 0.0E0;
  Double I_KINETIC_I3y3z_Px_aa = 0.0E0;
  Double I_KINETIC_I2y4z_Px_aa = 0.0E0;
  Double I_KINETIC_Iy5z_Px_aa = 0.0E0;
  Double I_KINETIC_I6z_Px_aa = 0.0E0;
  Double I_KINETIC_I6x_Py_aa = 0.0E0;
  Double I_KINETIC_I5xy_Py_aa = 0.0E0;
  Double I_KINETIC_I5xz_Py_aa = 0.0E0;
  Double I_KINETIC_I4x2y_Py_aa = 0.0E0;
  Double I_KINETIC_I4xyz_Py_aa = 0.0E0;
  Double I_KINETIC_I4x2z_Py_aa = 0.0E0;
  Double I_KINETIC_I3x3y_Py_aa = 0.0E0;
  Double I_KINETIC_I3x2yz_Py_aa = 0.0E0;
  Double I_KINETIC_I3xy2z_Py_aa = 0.0E0;
  Double I_KINETIC_I3x3z_Py_aa = 0.0E0;
  Double I_KINETIC_I2x4y_Py_aa = 0.0E0;
  Double I_KINETIC_I2x3yz_Py_aa = 0.0E0;
  Double I_KINETIC_I2x2y2z_Py_aa = 0.0E0;
  Double I_KINETIC_I2xy3z_Py_aa = 0.0E0;
  Double I_KINETIC_I2x4z_Py_aa = 0.0E0;
  Double I_KINETIC_Ix5y_Py_aa = 0.0E0;
  Double I_KINETIC_Ix4yz_Py_aa = 0.0E0;
  Double I_KINETIC_Ix3y2z_Py_aa = 0.0E0;
  Double I_KINETIC_Ix2y3z_Py_aa = 0.0E0;
  Double I_KINETIC_Ixy4z_Py_aa = 0.0E0;
  Double I_KINETIC_Ix5z_Py_aa = 0.0E0;
  Double I_KINETIC_I6y_Py_aa = 0.0E0;
  Double I_KINETIC_I5yz_Py_aa = 0.0E0;
  Double I_KINETIC_I4y2z_Py_aa = 0.0E0;
  Double I_KINETIC_I3y3z_Py_aa = 0.0E0;
  Double I_KINETIC_I2y4z_Py_aa = 0.0E0;
  Double I_KINETIC_Iy5z_Py_aa = 0.0E0;
  Double I_KINETIC_I6z_Py_aa = 0.0E0;
  Double I_KINETIC_I6x_Pz_aa = 0.0E0;
  Double I_KINETIC_I5xy_Pz_aa = 0.0E0;
  Double I_KINETIC_I5xz_Pz_aa = 0.0E0;
  Double I_KINETIC_I4x2y_Pz_aa = 0.0E0;
  Double I_KINETIC_I4xyz_Pz_aa = 0.0E0;
  Double I_KINETIC_I4x2z_Pz_aa = 0.0E0;
  Double I_KINETIC_I3x3y_Pz_aa = 0.0E0;
  Double I_KINETIC_I3x2yz_Pz_aa = 0.0E0;
  Double I_KINETIC_I3xy2z_Pz_aa = 0.0E0;
  Double I_KINETIC_I3x3z_Pz_aa = 0.0E0;
  Double I_KINETIC_I2x4y_Pz_aa = 0.0E0;
  Double I_KINETIC_I2x3yz_Pz_aa = 0.0E0;
  Double I_KINETIC_I2x2y2z_Pz_aa = 0.0E0;
  Double I_KINETIC_I2xy3z_Pz_aa = 0.0E0;
  Double I_KINETIC_I2x4z_Pz_aa = 0.0E0;
  Double I_KINETIC_Ix5y_Pz_aa = 0.0E0;
  Double I_KINETIC_Ix4yz_Pz_aa = 0.0E0;
  Double I_KINETIC_Ix3y2z_Pz_aa = 0.0E0;
  Double I_KINETIC_Ix2y3z_Pz_aa = 0.0E0;
  Double I_KINETIC_Ixy4z_Pz_aa = 0.0E0;
  Double I_KINETIC_Ix5z_Pz_aa = 0.0E0;
  Double I_KINETIC_I6y_Pz_aa = 0.0E0;
  Double I_KINETIC_I5yz_Pz_aa = 0.0E0;
  Double I_KINETIC_I4y2z_Pz_aa = 0.0E0;
  Double I_KINETIC_I3y3z_Pz_aa = 0.0E0;
  Double I_KINETIC_I2y4z_Pz_aa = 0.0E0;
  Double I_KINETIC_Iy5z_Pz_aa = 0.0E0;
  Double I_KINETIC_I6z_Pz_aa = 0.0E0;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_I_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_Px_vrr = PBX*I_TWOBODYOVERLAP_I6x_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Px_vrr = PBX*I_TWOBODYOVERLAP_I5xy_S_vrr+5*oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Px_vrr = PBX*I_TWOBODYOVERLAP_I5xz_S_vrr+5*oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Px_vrr = PBX*I_TWOBODYOVERLAP_I4x2y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Px_vrr = PBX*I_TWOBODYOVERLAP_I4xyz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Px_vrr = PBX*I_TWOBODYOVERLAP_I4x2z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Px_vrr = PBX*I_TWOBODYOVERLAP_I3x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_I3x2yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_I3xy2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Px_vrr = PBX*I_TWOBODYOVERLAP_I3x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Px_vrr = PBX*I_TWOBODYOVERLAP_I2x4y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Px_vrr = PBX*I_TWOBODYOVERLAP_I2x3yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_I2x2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Px_vrr = PBX*I_TWOBODYOVERLAP_I2xy3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Px_vrr = PBX*I_TWOBODYOVERLAP_I2x4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix5y_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix3y2z_S_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix2y3z_S_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Px_vrr = PBX*I_TWOBODYOVERLAP_Ixy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_Px_vrr = PBX*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Px_vrr = PBX*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Px_vrr = PBX*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Px_vrr = PBX*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Px_vrr = PBX*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_Px_vrr = PBX*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_I6x_Py_vrr = PBY*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Py_vrr = PBY*I_TWOBODYOVERLAP_I5xy_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Py_vrr = PBY*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Py_vrr = PBY*I_TWOBODYOVERLAP_I4x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Py_vrr = PBY*I_TWOBODYOVERLAP_I4xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Py_vrr = PBY*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Py_vrr = PBY*I_TWOBODYOVERLAP_I3x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_I3x2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_I3xy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Py_vrr = PBY*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Py_vrr = PBY*I_TWOBODYOVERLAP_I2x4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Py_vrr = PBY*I_TWOBODYOVERLAP_I2x3yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_I2x2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Py_vrr = PBY*I_TWOBODYOVERLAP_I2xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Py_vrr = PBY*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix5y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix4yz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix3y2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix2y3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Py_vrr = PBY*I_TWOBODYOVERLAP_Ixy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_Py_vrr = PBY*I_TWOBODYOVERLAP_I6y_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Py_vrr = PBY*I_TWOBODYOVERLAP_I5yz_S_vrr+5*oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_I4y2z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Py_vrr = PBY*I_TWOBODYOVERLAP_I3y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Py_vrr = PBY*I_TWOBODYOVERLAP_I2y4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Py_vrr = PBY*I_TWOBODYOVERLAP_Iy5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_Py_vrr = PBY*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_I6x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I5xz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I4xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I4x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3xy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2x2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2xy3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2x4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix2y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ixy4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I5yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I4y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2y4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Iy5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I6z_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;

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
     * shell quartet name: SQ_KINETIC_I_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_I_S
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P
     ************************************************************/
    Double I_KINETIC_I6x_Px_vrr = PBX*I_KINETIC_I6x_S_vrr+6*oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Px_vrr;
    Double I_KINETIC_I5xy_Px_vrr = PBX*I_KINETIC_I5xy_S_vrr+5*oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Px_vrr;
    Double I_KINETIC_I5xz_Px_vrr = PBX*I_KINETIC_I5xz_S_vrr+5*oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Px_vrr;
    Double I_KINETIC_I4x2y_Px_vrr = PBX*I_KINETIC_I4x2y_S_vrr+4*oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Px_vrr;
    Double I_KINETIC_I4xyz_Px_vrr = PBX*I_KINETIC_I4xyz_S_vrr+4*oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Px_vrr;
    Double I_KINETIC_I4x2z_Px_vrr = PBX*I_KINETIC_I4x2z_S_vrr+4*oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Px_vrr;
    Double I_KINETIC_I3x3y_Px_vrr = PBX*I_KINETIC_I3x3y_S_vrr+3*oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Px_vrr;
    Double I_KINETIC_I3x2yz_Px_vrr = PBX*I_KINETIC_I3x2yz_S_vrr+3*oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Px_vrr;
    Double I_KINETIC_I3xy2z_Px_vrr = PBX*I_KINETIC_I3xy2z_S_vrr+3*oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Px_vrr;
    Double I_KINETIC_I3x3z_Px_vrr = PBX*I_KINETIC_I3x3z_S_vrr+3*oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Px_vrr;
    Double I_KINETIC_I2x4y_Px_vrr = PBX*I_KINETIC_I2x4y_S_vrr+2*oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Px_vrr;
    Double I_KINETIC_I2x3yz_Px_vrr = PBX*I_KINETIC_I2x3yz_S_vrr+2*oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Px_vrr;
    Double I_KINETIC_I2x2y2z_Px_vrr = PBX*I_KINETIC_I2x2y2z_S_vrr+2*oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Px_vrr;
    Double I_KINETIC_I2xy3z_Px_vrr = PBX*I_KINETIC_I2xy3z_S_vrr+2*oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Px_vrr;
    Double I_KINETIC_I2x4z_Px_vrr = PBX*I_KINETIC_I2x4z_S_vrr+2*oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Px_vrr;
    Double I_KINETIC_Ix5y_Px_vrr = PBX*I_KINETIC_Ix5y_S_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Px_vrr;
    Double I_KINETIC_Ix4yz_Px_vrr = PBX*I_KINETIC_Ix4yz_S_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Px_vrr;
    Double I_KINETIC_Ix3y2z_Px_vrr = PBX*I_KINETIC_Ix3y2z_S_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Px_vrr;
    Double I_KINETIC_Ix2y3z_Px_vrr = PBX*I_KINETIC_Ix2y3z_S_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Px_vrr;
    Double I_KINETIC_Ixy4z_Px_vrr = PBX*I_KINETIC_Ixy4z_S_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Px_vrr;
    Double I_KINETIC_Ix5z_Px_vrr = PBX*I_KINETIC_Ix5z_S_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Px_vrr;
    Double I_KINETIC_I6y_Px_vrr = PBX*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Px_vrr;
    Double I_KINETIC_I5yz_Px_vrr = PBX*I_KINETIC_I5yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Px_vrr;
    Double I_KINETIC_I4y2z_Px_vrr = PBX*I_KINETIC_I4y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Px_vrr;
    Double I_KINETIC_I3y3z_Px_vrr = PBX*I_KINETIC_I3y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Px_vrr;
    Double I_KINETIC_I2y4z_Px_vrr = PBX*I_KINETIC_I2y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Px_vrr;
    Double I_KINETIC_Iy5z_Px_vrr = PBX*I_KINETIC_Iy5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Px_vrr;
    Double I_KINETIC_I6z_Px_vrr = PBX*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Px_vrr;
    Double I_KINETIC_I6x_Py_vrr = PBY*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Py_vrr;
    Double I_KINETIC_I5xy_Py_vrr = PBY*I_KINETIC_I5xy_S_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Py_vrr;
    Double I_KINETIC_I5xz_Py_vrr = PBY*I_KINETIC_I5xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Py_vrr;
    Double I_KINETIC_I4x2y_Py_vrr = PBY*I_KINETIC_I4x2y_S_vrr+2*oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Py_vrr;
    Double I_KINETIC_I4xyz_Py_vrr = PBY*I_KINETIC_I4xyz_S_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Py_vrr;
    Double I_KINETIC_I4x2z_Py_vrr = PBY*I_KINETIC_I4x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Py_vrr;
    Double I_KINETIC_I3x3y_Py_vrr = PBY*I_KINETIC_I3x3y_S_vrr+3*oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Py_vrr;
    Double I_KINETIC_I3x2yz_Py_vrr = PBY*I_KINETIC_I3x2yz_S_vrr+2*oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Py_vrr;
    Double I_KINETIC_I3xy2z_Py_vrr = PBY*I_KINETIC_I3xy2z_S_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Py_vrr;
    Double I_KINETIC_I3x3z_Py_vrr = PBY*I_KINETIC_I3x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Py_vrr;
    Double I_KINETIC_I2x4y_Py_vrr = PBY*I_KINETIC_I2x4y_S_vrr+4*oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Py_vrr;
    Double I_KINETIC_I2x3yz_Py_vrr = PBY*I_KINETIC_I2x3yz_S_vrr+3*oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Py_vrr;
    Double I_KINETIC_I2x2y2z_Py_vrr = PBY*I_KINETIC_I2x2y2z_S_vrr+2*oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Py_vrr;
    Double I_KINETIC_I2xy3z_Py_vrr = PBY*I_KINETIC_I2xy3z_S_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Py_vrr;
    Double I_KINETIC_I2x4z_Py_vrr = PBY*I_KINETIC_I2x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Py_vrr;
    Double I_KINETIC_Ix5y_Py_vrr = PBY*I_KINETIC_Ix5y_S_vrr+5*oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Py_vrr;
    Double I_KINETIC_Ix4yz_Py_vrr = PBY*I_KINETIC_Ix4yz_S_vrr+4*oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Py_vrr;
    Double I_KINETIC_Ix3y2z_Py_vrr = PBY*I_KINETIC_Ix3y2z_S_vrr+3*oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Py_vrr;
    Double I_KINETIC_Ix2y3z_Py_vrr = PBY*I_KINETIC_Ix2y3z_S_vrr+2*oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Py_vrr;
    Double I_KINETIC_Ixy4z_Py_vrr = PBY*I_KINETIC_Ixy4z_S_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Py_vrr;
    Double I_KINETIC_Ix5z_Py_vrr = PBY*I_KINETIC_Ix5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Py_vrr;
    Double I_KINETIC_I6y_Py_vrr = PBY*I_KINETIC_I6y_S_vrr+6*oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Py_vrr;
    Double I_KINETIC_I5yz_Py_vrr = PBY*I_KINETIC_I5yz_S_vrr+5*oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Py_vrr;
    Double I_KINETIC_I4y2z_Py_vrr = PBY*I_KINETIC_I4y2z_S_vrr+4*oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Py_vrr;
    Double I_KINETIC_I3y3z_Py_vrr = PBY*I_KINETIC_I3y3z_S_vrr+3*oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Py_vrr;
    Double I_KINETIC_I2y4z_Py_vrr = PBY*I_KINETIC_I2y4z_S_vrr+2*oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Py_vrr;
    Double I_KINETIC_Iy5z_Py_vrr = PBY*I_KINETIC_Iy5z_S_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Py_vrr;
    Double I_KINETIC_I6z_Py_vrr = PBY*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Py_vrr;
    Double I_KINETIC_I6x_Pz_vrr = PBZ*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Pz_vrr;
    Double I_KINETIC_I5xy_Pz_vrr = PBZ*I_KINETIC_I5xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Pz_vrr;
    Double I_KINETIC_I5xz_Pz_vrr = PBZ*I_KINETIC_I5xz_S_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Pz_vrr;
    Double I_KINETIC_I4x2y_Pz_vrr = PBZ*I_KINETIC_I4x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Pz_vrr;
    Double I_KINETIC_I4xyz_Pz_vrr = PBZ*I_KINETIC_I4xyz_S_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Pz_vrr;
    Double I_KINETIC_I4x2z_Pz_vrr = PBZ*I_KINETIC_I4x2z_S_vrr+2*oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Pz_vrr;
    Double I_KINETIC_I3x3y_Pz_vrr = PBZ*I_KINETIC_I3x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Pz_vrr;
    Double I_KINETIC_I3x2yz_Pz_vrr = PBZ*I_KINETIC_I3x2yz_S_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Pz_vrr;
    Double I_KINETIC_I3xy2z_Pz_vrr = PBZ*I_KINETIC_I3xy2z_S_vrr+2*oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Pz_vrr;
    Double I_KINETIC_I3x3z_Pz_vrr = PBZ*I_KINETIC_I3x3z_S_vrr+3*oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Pz_vrr;
    Double I_KINETIC_I2x4y_Pz_vrr = PBZ*I_KINETIC_I2x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Pz_vrr;
    Double I_KINETIC_I2x3yz_Pz_vrr = PBZ*I_KINETIC_I2x3yz_S_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Pz_vrr;
    Double I_KINETIC_I2x2y2z_Pz_vrr = PBZ*I_KINETIC_I2x2y2z_S_vrr+2*oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Pz_vrr;
    Double I_KINETIC_I2xy3z_Pz_vrr = PBZ*I_KINETIC_I2xy3z_S_vrr+3*oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Pz_vrr;
    Double I_KINETIC_I2x4z_Pz_vrr = PBZ*I_KINETIC_I2x4z_S_vrr+4*oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Pz_vrr;
    Double I_KINETIC_Ix5y_Pz_vrr = PBZ*I_KINETIC_Ix5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Pz_vrr;
    Double I_KINETIC_Ix4yz_Pz_vrr = PBZ*I_KINETIC_Ix4yz_S_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Pz_vrr;
    Double I_KINETIC_Ix3y2z_Pz_vrr = PBZ*I_KINETIC_Ix3y2z_S_vrr+2*oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Pz_vrr;
    Double I_KINETIC_Ix2y3z_Pz_vrr = PBZ*I_KINETIC_Ix2y3z_S_vrr+3*oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Pz_vrr;
    Double I_KINETIC_Ixy4z_Pz_vrr = PBZ*I_KINETIC_Ixy4z_S_vrr+4*oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Pz_vrr;
    Double I_KINETIC_Ix5z_Pz_vrr = PBZ*I_KINETIC_Ix5z_S_vrr+5*oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Pz_vrr;
    Double I_KINETIC_I6y_Pz_vrr = PBZ*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Pz_vrr;
    Double I_KINETIC_I5yz_Pz_vrr = PBZ*I_KINETIC_I5yz_S_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Pz_vrr;
    Double I_KINETIC_I4y2z_Pz_vrr = PBZ*I_KINETIC_I4y2z_S_vrr+2*oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Pz_vrr;
    Double I_KINETIC_I3y3z_Pz_vrr = PBZ*I_KINETIC_I3y3z_S_vrr+3*oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Pz_vrr;
    Double I_KINETIC_I2y4z_Pz_vrr = PBZ*I_KINETIC_I2y4z_S_vrr+4*oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Pz_vrr;
    Double I_KINETIC_Iy5z_Pz_vrr = PBZ*I_KINETIC_Iy5z_S_vrr+5*oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Pz_vrr;
    Double I_KINETIC_I6z_Pz_vrr = PBZ*I_KINETIC_I6z_S_vrr+6*oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_P_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_I_P_aa_coefs = alpha*alpha;
    I_KINETIC_I6x_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I6x_Px_vrr;
    I_KINETIC_I5xy_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I5xy_Px_vrr;
    I_KINETIC_I5xz_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I5xz_Px_vrr;
    I_KINETIC_I4x2y_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4x2y_Px_vrr;
    I_KINETIC_I4xyz_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4xyz_Px_vrr;
    I_KINETIC_I4x2z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4x2z_Px_vrr;
    I_KINETIC_I3x3y_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3x3y_Px_vrr;
    I_KINETIC_I3x2yz_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3x2yz_Px_vrr;
    I_KINETIC_I3xy2z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3xy2z_Px_vrr;
    I_KINETIC_I3x3z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3x3z_Px_vrr;
    I_KINETIC_I2x4y_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x4y_Px_vrr;
    I_KINETIC_I2x3yz_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x3yz_Px_vrr;
    I_KINETIC_I2x2y2z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x2y2z_Px_vrr;
    I_KINETIC_I2xy3z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2xy3z_Px_vrr;
    I_KINETIC_I2x4z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x4z_Px_vrr;
    I_KINETIC_Ix5y_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix5y_Px_vrr;
    I_KINETIC_Ix4yz_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix4yz_Px_vrr;
    I_KINETIC_Ix3y2z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix3y2z_Px_vrr;
    I_KINETIC_Ix2y3z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix2y3z_Px_vrr;
    I_KINETIC_Ixy4z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ixy4z_Px_vrr;
    I_KINETIC_Ix5z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix5z_Px_vrr;
    I_KINETIC_I6y_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I6y_Px_vrr;
    I_KINETIC_I5yz_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I5yz_Px_vrr;
    I_KINETIC_I4y2z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4y2z_Px_vrr;
    I_KINETIC_I3y3z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3y3z_Px_vrr;
    I_KINETIC_I2y4z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2y4z_Px_vrr;
    I_KINETIC_Iy5z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Iy5z_Px_vrr;
    I_KINETIC_I6z_Px_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I6z_Px_vrr;
    I_KINETIC_I6x_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I6x_Py_vrr;
    I_KINETIC_I5xy_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I5xy_Py_vrr;
    I_KINETIC_I5xz_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I5xz_Py_vrr;
    I_KINETIC_I4x2y_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4x2y_Py_vrr;
    I_KINETIC_I4xyz_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4xyz_Py_vrr;
    I_KINETIC_I4x2z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4x2z_Py_vrr;
    I_KINETIC_I3x3y_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3x3y_Py_vrr;
    I_KINETIC_I3x2yz_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3x2yz_Py_vrr;
    I_KINETIC_I3xy2z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3xy2z_Py_vrr;
    I_KINETIC_I3x3z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3x3z_Py_vrr;
    I_KINETIC_I2x4y_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x4y_Py_vrr;
    I_KINETIC_I2x3yz_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x3yz_Py_vrr;
    I_KINETIC_I2x2y2z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x2y2z_Py_vrr;
    I_KINETIC_I2xy3z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2xy3z_Py_vrr;
    I_KINETIC_I2x4z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x4z_Py_vrr;
    I_KINETIC_Ix5y_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix5y_Py_vrr;
    I_KINETIC_Ix4yz_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix4yz_Py_vrr;
    I_KINETIC_Ix3y2z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix3y2z_Py_vrr;
    I_KINETIC_Ix2y3z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix2y3z_Py_vrr;
    I_KINETIC_Ixy4z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ixy4z_Py_vrr;
    I_KINETIC_Ix5z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix5z_Py_vrr;
    I_KINETIC_I6y_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I6y_Py_vrr;
    I_KINETIC_I5yz_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I5yz_Py_vrr;
    I_KINETIC_I4y2z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4y2z_Py_vrr;
    I_KINETIC_I3y3z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3y3z_Py_vrr;
    I_KINETIC_I2y4z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2y4z_Py_vrr;
    I_KINETIC_Iy5z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Iy5z_Py_vrr;
    I_KINETIC_I6z_Py_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I6z_Py_vrr;
    I_KINETIC_I6x_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I6x_Pz_vrr;
    I_KINETIC_I5xy_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I5xy_Pz_vrr;
    I_KINETIC_I5xz_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I5xz_Pz_vrr;
    I_KINETIC_I4x2y_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4x2y_Pz_vrr;
    I_KINETIC_I4xyz_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4xyz_Pz_vrr;
    I_KINETIC_I4x2z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4x2z_Pz_vrr;
    I_KINETIC_I3x3y_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3x3y_Pz_vrr;
    I_KINETIC_I3x2yz_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3x2yz_Pz_vrr;
    I_KINETIC_I3xy2z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3xy2z_Pz_vrr;
    I_KINETIC_I3x3z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3x3z_Pz_vrr;
    I_KINETIC_I2x4y_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x4y_Pz_vrr;
    I_KINETIC_I2x3yz_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x3yz_Pz_vrr;
    I_KINETIC_I2x2y2z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x2y2z_Pz_vrr;
    I_KINETIC_I2xy3z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2xy3z_Pz_vrr;
    I_KINETIC_I2x4z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2x4z_Pz_vrr;
    I_KINETIC_Ix5y_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix5y_Pz_vrr;
    I_KINETIC_Ix4yz_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix4yz_Pz_vrr;
    I_KINETIC_Ix3y2z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix3y2z_Pz_vrr;
    I_KINETIC_Ix2y3z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix2y3z_Pz_vrr;
    I_KINETIC_Ixy4z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ixy4z_Pz_vrr;
    I_KINETIC_Ix5z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Ix5z_Pz_vrr;
    I_KINETIC_I6y_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I6y_Pz_vrr;
    I_KINETIC_I5yz_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I5yz_Pz_vrr;
    I_KINETIC_I4y2z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I4y2z_Pz_vrr;
    I_KINETIC_I3y3z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I3y3z_Pz_vrr;
    I_KINETIC_I2y4z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I2y4z_Pz_vrr;
    I_KINETIC_Iy5z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_Iy5z_Pz_vrr;
    I_KINETIC_I6z_Pz_aa += SQ_KINETIC_I_P_aa_coefs*I_KINETIC_I6z_Pz_vrr;

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
   * shell quartet name: SQ_KINETIC_G_P_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_P_aa
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_D_P
   ************************************************************/
  abcd[0] = 4.0E0*I_KINETIC_I6x_Px_aa-2.0E0*4*I_KINETIC_G4x_Px_a-2.0E0*5*I_KINETIC_G4x_Px_a+4*3*I_KINETIC_D2x_Px;
  abcd[1] = 4.0E0*I_KINETIC_I5xy_Px_aa-2.0E0*3*I_KINETIC_G3xy_Px_a-2.0E0*4*I_KINETIC_G3xy_Px_a+3*2*I_KINETIC_Dxy_Px;
  abcd[2] = 4.0E0*I_KINETIC_I5xz_Px_aa-2.0E0*3*I_KINETIC_G3xz_Px_a-2.0E0*4*I_KINETIC_G3xz_Px_a+3*2*I_KINETIC_Dxz_Px;
  abcd[3] = 4.0E0*I_KINETIC_I4x2y_Px_aa-2.0E0*2*I_KINETIC_G2x2y_Px_a-2.0E0*3*I_KINETIC_G2x2y_Px_a+2*1*I_KINETIC_D2y_Px;
  abcd[4] = 4.0E0*I_KINETIC_I4xyz_Px_aa-2.0E0*2*I_KINETIC_G2xyz_Px_a-2.0E0*3*I_KINETIC_G2xyz_Px_a+2*1*I_KINETIC_Dyz_Px;
  abcd[5] = 4.0E0*I_KINETIC_I4x2z_Px_aa-2.0E0*2*I_KINETIC_G2x2z_Px_a-2.0E0*3*I_KINETIC_G2x2z_Px_a+2*1*I_KINETIC_D2z_Px;
  abcd[6] = 4.0E0*I_KINETIC_I3x3y_Px_aa-2.0E0*1*I_KINETIC_Gx3y_Px_a-2.0E0*2*I_KINETIC_Gx3y_Px_a;
  abcd[7] = 4.0E0*I_KINETIC_I3x2yz_Px_aa-2.0E0*1*I_KINETIC_Gx2yz_Px_a-2.0E0*2*I_KINETIC_Gx2yz_Px_a;
  abcd[8] = 4.0E0*I_KINETIC_I3xy2z_Px_aa-2.0E0*1*I_KINETIC_Gxy2z_Px_a-2.0E0*2*I_KINETIC_Gxy2z_Px_a;
  abcd[9] = 4.0E0*I_KINETIC_I3x3z_Px_aa-2.0E0*1*I_KINETIC_Gx3z_Px_a-2.0E0*2*I_KINETIC_Gx3z_Px_a;
  abcd[10] = 4.0E0*I_KINETIC_I2x4y_Px_aa-2.0E0*1*I_KINETIC_G4y_Px_a;
  abcd[11] = 4.0E0*I_KINETIC_I2x3yz_Px_aa-2.0E0*1*I_KINETIC_G3yz_Px_a;
  abcd[12] = 4.0E0*I_KINETIC_I2x2y2z_Px_aa-2.0E0*1*I_KINETIC_G2y2z_Px_a;
  abcd[13] = 4.0E0*I_KINETIC_I2xy3z_Px_aa-2.0E0*1*I_KINETIC_Gy3z_Px_a;
  abcd[14] = 4.0E0*I_KINETIC_I2x4z_Px_aa-2.0E0*1*I_KINETIC_G4z_Px_a;
  abcd[15] = 4.0E0*I_KINETIC_I6x_Py_aa-2.0E0*4*I_KINETIC_G4x_Py_a-2.0E0*5*I_KINETIC_G4x_Py_a+4*3*I_KINETIC_D2x_Py;
  abcd[16] = 4.0E0*I_KINETIC_I5xy_Py_aa-2.0E0*3*I_KINETIC_G3xy_Py_a-2.0E0*4*I_KINETIC_G3xy_Py_a+3*2*I_KINETIC_Dxy_Py;
  abcd[17] = 4.0E0*I_KINETIC_I5xz_Py_aa-2.0E0*3*I_KINETIC_G3xz_Py_a-2.0E0*4*I_KINETIC_G3xz_Py_a+3*2*I_KINETIC_Dxz_Py;
  abcd[18] = 4.0E0*I_KINETIC_I4x2y_Py_aa-2.0E0*2*I_KINETIC_G2x2y_Py_a-2.0E0*3*I_KINETIC_G2x2y_Py_a+2*1*I_KINETIC_D2y_Py;
  abcd[19] = 4.0E0*I_KINETIC_I4xyz_Py_aa-2.0E0*2*I_KINETIC_G2xyz_Py_a-2.0E0*3*I_KINETIC_G2xyz_Py_a+2*1*I_KINETIC_Dyz_Py;
  abcd[20] = 4.0E0*I_KINETIC_I4x2z_Py_aa-2.0E0*2*I_KINETIC_G2x2z_Py_a-2.0E0*3*I_KINETIC_G2x2z_Py_a+2*1*I_KINETIC_D2z_Py;
  abcd[21] = 4.0E0*I_KINETIC_I3x3y_Py_aa-2.0E0*1*I_KINETIC_Gx3y_Py_a-2.0E0*2*I_KINETIC_Gx3y_Py_a;
  abcd[22] = 4.0E0*I_KINETIC_I3x2yz_Py_aa-2.0E0*1*I_KINETIC_Gx2yz_Py_a-2.0E0*2*I_KINETIC_Gx2yz_Py_a;
  abcd[23] = 4.0E0*I_KINETIC_I3xy2z_Py_aa-2.0E0*1*I_KINETIC_Gxy2z_Py_a-2.0E0*2*I_KINETIC_Gxy2z_Py_a;
  abcd[24] = 4.0E0*I_KINETIC_I3x3z_Py_aa-2.0E0*1*I_KINETIC_Gx3z_Py_a-2.0E0*2*I_KINETIC_Gx3z_Py_a;
  abcd[25] = 4.0E0*I_KINETIC_I2x4y_Py_aa-2.0E0*1*I_KINETIC_G4y_Py_a;
  abcd[26] = 4.0E0*I_KINETIC_I2x3yz_Py_aa-2.0E0*1*I_KINETIC_G3yz_Py_a;
  abcd[27] = 4.0E0*I_KINETIC_I2x2y2z_Py_aa-2.0E0*1*I_KINETIC_G2y2z_Py_a;
  abcd[28] = 4.0E0*I_KINETIC_I2xy3z_Py_aa-2.0E0*1*I_KINETIC_Gy3z_Py_a;
  abcd[29] = 4.0E0*I_KINETIC_I2x4z_Py_aa-2.0E0*1*I_KINETIC_G4z_Py_a;
  abcd[30] = 4.0E0*I_KINETIC_I6x_Pz_aa-2.0E0*4*I_KINETIC_G4x_Pz_a-2.0E0*5*I_KINETIC_G4x_Pz_a+4*3*I_KINETIC_D2x_Pz;
  abcd[31] = 4.0E0*I_KINETIC_I5xy_Pz_aa-2.0E0*3*I_KINETIC_G3xy_Pz_a-2.0E0*4*I_KINETIC_G3xy_Pz_a+3*2*I_KINETIC_Dxy_Pz;
  abcd[32] = 4.0E0*I_KINETIC_I5xz_Pz_aa-2.0E0*3*I_KINETIC_G3xz_Pz_a-2.0E0*4*I_KINETIC_G3xz_Pz_a+3*2*I_KINETIC_Dxz_Pz;
  abcd[33] = 4.0E0*I_KINETIC_I4x2y_Pz_aa-2.0E0*2*I_KINETIC_G2x2y_Pz_a-2.0E0*3*I_KINETIC_G2x2y_Pz_a+2*1*I_KINETIC_D2y_Pz;
  abcd[34] = 4.0E0*I_KINETIC_I4xyz_Pz_aa-2.0E0*2*I_KINETIC_G2xyz_Pz_a-2.0E0*3*I_KINETIC_G2xyz_Pz_a+2*1*I_KINETIC_Dyz_Pz;
  abcd[35] = 4.0E0*I_KINETIC_I4x2z_Pz_aa-2.0E0*2*I_KINETIC_G2x2z_Pz_a-2.0E0*3*I_KINETIC_G2x2z_Pz_a+2*1*I_KINETIC_D2z_Pz;
  abcd[36] = 4.0E0*I_KINETIC_I3x3y_Pz_aa-2.0E0*1*I_KINETIC_Gx3y_Pz_a-2.0E0*2*I_KINETIC_Gx3y_Pz_a;
  abcd[37] = 4.0E0*I_KINETIC_I3x2yz_Pz_aa-2.0E0*1*I_KINETIC_Gx2yz_Pz_a-2.0E0*2*I_KINETIC_Gx2yz_Pz_a;
  abcd[38] = 4.0E0*I_KINETIC_I3xy2z_Pz_aa-2.0E0*1*I_KINETIC_Gxy2z_Pz_a-2.0E0*2*I_KINETIC_Gxy2z_Pz_a;
  abcd[39] = 4.0E0*I_KINETIC_I3x3z_Pz_aa-2.0E0*1*I_KINETIC_Gx3z_Pz_a-2.0E0*2*I_KINETIC_Gx3z_Pz_a;
  abcd[40] = 4.0E0*I_KINETIC_I2x4y_Pz_aa-2.0E0*1*I_KINETIC_G4y_Pz_a;
  abcd[41] = 4.0E0*I_KINETIC_I2x3yz_Pz_aa-2.0E0*1*I_KINETIC_G3yz_Pz_a;
  abcd[42] = 4.0E0*I_KINETIC_I2x2y2z_Pz_aa-2.0E0*1*I_KINETIC_G2y2z_Pz_a;
  abcd[43] = 4.0E0*I_KINETIC_I2xy3z_Pz_aa-2.0E0*1*I_KINETIC_Gy3z_Pz_a;
  abcd[44] = 4.0E0*I_KINETIC_I2x4z_Pz_aa-2.0E0*1*I_KINETIC_G4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_P_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_P_aa
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_D_P
   ************************************************************/
  abcd[45] = 4.0E0*I_KINETIC_I5xy_Px_aa-2.0E0*4*I_KINETIC_G3xy_Px_a;
  abcd[46] = 4.0E0*I_KINETIC_I4x2y_Px_aa-2.0E0*1*I_KINETIC_G4x_Px_a-2.0E0*3*I_KINETIC_G2x2y_Px_a+3*1*I_KINETIC_D2x_Px;
  abcd[47] = 4.0E0*I_KINETIC_I4xyz_Px_aa-2.0E0*3*I_KINETIC_G2xyz_Px_a;
  abcd[48] = 4.0E0*I_KINETIC_I3x3y_Px_aa-2.0E0*2*I_KINETIC_G3xy_Px_a-2.0E0*2*I_KINETIC_Gx3y_Px_a+2*2*I_KINETIC_Dxy_Px;
  abcd[49] = 4.0E0*I_KINETIC_I3x2yz_Px_aa-2.0E0*1*I_KINETIC_G3xz_Px_a-2.0E0*2*I_KINETIC_Gx2yz_Px_a+2*1*I_KINETIC_Dxz_Px;
  abcd[50] = 4.0E0*I_KINETIC_I3xy2z_Px_aa-2.0E0*2*I_KINETIC_Gxy2z_Px_a;
  abcd[51] = 4.0E0*I_KINETIC_I2x4y_Px_aa-2.0E0*3*I_KINETIC_G2x2y_Px_a-2.0E0*1*I_KINETIC_G4y_Px_a+3*I_KINETIC_D2y_Px;
  abcd[52] = 4.0E0*I_KINETIC_I2x3yz_Px_aa-2.0E0*2*I_KINETIC_G2xyz_Px_a-2.0E0*1*I_KINETIC_G3yz_Px_a+2*I_KINETIC_Dyz_Px;
  abcd[53] = 4.0E0*I_KINETIC_I2x2y2z_Px_aa-2.0E0*1*I_KINETIC_G2x2z_Px_a-2.0E0*1*I_KINETIC_G2y2z_Px_a+1*I_KINETIC_D2z_Px;
  abcd[54] = 4.0E0*I_KINETIC_I2xy3z_Px_aa-2.0E0*1*I_KINETIC_Gy3z_Px_a;
  abcd[55] = 4.0E0*I_KINETIC_Ix5y_Px_aa-2.0E0*4*I_KINETIC_Gx3y_Px_a;
  abcd[56] = 4.0E0*I_KINETIC_Ix4yz_Px_aa-2.0E0*3*I_KINETIC_Gx2yz_Px_a;
  abcd[57] = 4.0E0*I_KINETIC_Ix3y2z_Px_aa-2.0E0*2*I_KINETIC_Gxy2z_Px_a;
  abcd[58] = 4.0E0*I_KINETIC_Ix2y3z_Px_aa-2.0E0*1*I_KINETIC_Gx3z_Px_a;
  abcd[59] = 4.0E0*I_KINETIC_Ixy4z_Px_aa;
  abcd[60] = 4.0E0*I_KINETIC_I5xy_Py_aa-2.0E0*4*I_KINETIC_G3xy_Py_a;
  abcd[61] = 4.0E0*I_KINETIC_I4x2y_Py_aa-2.0E0*1*I_KINETIC_G4x_Py_a-2.0E0*3*I_KINETIC_G2x2y_Py_a+3*1*I_KINETIC_D2x_Py;
  abcd[62] = 4.0E0*I_KINETIC_I4xyz_Py_aa-2.0E0*3*I_KINETIC_G2xyz_Py_a;
  abcd[63] = 4.0E0*I_KINETIC_I3x3y_Py_aa-2.0E0*2*I_KINETIC_G3xy_Py_a-2.0E0*2*I_KINETIC_Gx3y_Py_a+2*2*I_KINETIC_Dxy_Py;
  abcd[64] = 4.0E0*I_KINETIC_I3x2yz_Py_aa-2.0E0*1*I_KINETIC_G3xz_Py_a-2.0E0*2*I_KINETIC_Gx2yz_Py_a+2*1*I_KINETIC_Dxz_Py;
  abcd[65] = 4.0E0*I_KINETIC_I3xy2z_Py_aa-2.0E0*2*I_KINETIC_Gxy2z_Py_a;
  abcd[66] = 4.0E0*I_KINETIC_I2x4y_Py_aa-2.0E0*3*I_KINETIC_G2x2y_Py_a-2.0E0*1*I_KINETIC_G4y_Py_a+3*I_KINETIC_D2y_Py;
  abcd[67] = 4.0E0*I_KINETIC_I2x3yz_Py_aa-2.0E0*2*I_KINETIC_G2xyz_Py_a-2.0E0*1*I_KINETIC_G3yz_Py_a+2*I_KINETIC_Dyz_Py;
  abcd[68] = 4.0E0*I_KINETIC_I2x2y2z_Py_aa-2.0E0*1*I_KINETIC_G2x2z_Py_a-2.0E0*1*I_KINETIC_G2y2z_Py_a+1*I_KINETIC_D2z_Py;
  abcd[69] = 4.0E0*I_KINETIC_I2xy3z_Py_aa-2.0E0*1*I_KINETIC_Gy3z_Py_a;
  abcd[70] = 4.0E0*I_KINETIC_Ix5y_Py_aa-2.0E0*4*I_KINETIC_Gx3y_Py_a;
  abcd[71] = 4.0E0*I_KINETIC_Ix4yz_Py_aa-2.0E0*3*I_KINETIC_Gx2yz_Py_a;
  abcd[72] = 4.0E0*I_KINETIC_Ix3y2z_Py_aa-2.0E0*2*I_KINETIC_Gxy2z_Py_a;
  abcd[73] = 4.0E0*I_KINETIC_Ix2y3z_Py_aa-2.0E0*1*I_KINETIC_Gx3z_Py_a;
  abcd[74] = 4.0E0*I_KINETIC_Ixy4z_Py_aa;
  abcd[75] = 4.0E0*I_KINETIC_I5xy_Pz_aa-2.0E0*4*I_KINETIC_G3xy_Pz_a;
  abcd[76] = 4.0E0*I_KINETIC_I4x2y_Pz_aa-2.0E0*1*I_KINETIC_G4x_Pz_a-2.0E0*3*I_KINETIC_G2x2y_Pz_a+3*1*I_KINETIC_D2x_Pz;
  abcd[77] = 4.0E0*I_KINETIC_I4xyz_Pz_aa-2.0E0*3*I_KINETIC_G2xyz_Pz_a;
  abcd[78] = 4.0E0*I_KINETIC_I3x3y_Pz_aa-2.0E0*2*I_KINETIC_G3xy_Pz_a-2.0E0*2*I_KINETIC_Gx3y_Pz_a+2*2*I_KINETIC_Dxy_Pz;
  abcd[79] = 4.0E0*I_KINETIC_I3x2yz_Pz_aa-2.0E0*1*I_KINETIC_G3xz_Pz_a-2.0E0*2*I_KINETIC_Gx2yz_Pz_a+2*1*I_KINETIC_Dxz_Pz;
  abcd[80] = 4.0E0*I_KINETIC_I3xy2z_Pz_aa-2.0E0*2*I_KINETIC_Gxy2z_Pz_a;
  abcd[81] = 4.0E0*I_KINETIC_I2x4y_Pz_aa-2.0E0*3*I_KINETIC_G2x2y_Pz_a-2.0E0*1*I_KINETIC_G4y_Pz_a+3*I_KINETIC_D2y_Pz;
  abcd[82] = 4.0E0*I_KINETIC_I2x3yz_Pz_aa-2.0E0*2*I_KINETIC_G2xyz_Pz_a-2.0E0*1*I_KINETIC_G3yz_Pz_a+2*I_KINETIC_Dyz_Pz;
  abcd[83] = 4.0E0*I_KINETIC_I2x2y2z_Pz_aa-2.0E0*1*I_KINETIC_G2x2z_Pz_a-2.0E0*1*I_KINETIC_G2y2z_Pz_a+1*I_KINETIC_D2z_Pz;
  abcd[84] = 4.0E0*I_KINETIC_I2xy3z_Pz_aa-2.0E0*1*I_KINETIC_Gy3z_Pz_a;
  abcd[85] = 4.0E0*I_KINETIC_Ix5y_Pz_aa-2.0E0*4*I_KINETIC_Gx3y_Pz_a;
  abcd[86] = 4.0E0*I_KINETIC_Ix4yz_Pz_aa-2.0E0*3*I_KINETIC_Gx2yz_Pz_a;
  abcd[87] = 4.0E0*I_KINETIC_Ix3y2z_Pz_aa-2.0E0*2*I_KINETIC_Gxy2z_Pz_a;
  abcd[88] = 4.0E0*I_KINETIC_Ix2y3z_Pz_aa-2.0E0*1*I_KINETIC_Gx3z_Pz_a;
  abcd[89] = 4.0E0*I_KINETIC_Ixy4z_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_P_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_P_aa
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_D_P
   ************************************************************/
  abcd[90] = 4.0E0*I_KINETIC_I5xz_Px_aa-2.0E0*4*I_KINETIC_G3xz_Px_a;
  abcd[91] = 4.0E0*I_KINETIC_I4xyz_Px_aa-2.0E0*3*I_KINETIC_G2xyz_Px_a;
  abcd[92] = 4.0E0*I_KINETIC_I4x2z_Px_aa-2.0E0*1*I_KINETIC_G4x_Px_a-2.0E0*3*I_KINETIC_G2x2z_Px_a+3*1*I_KINETIC_D2x_Px;
  abcd[93] = 4.0E0*I_KINETIC_I3x2yz_Px_aa-2.0E0*2*I_KINETIC_Gx2yz_Px_a;
  abcd[94] = 4.0E0*I_KINETIC_I3xy2z_Px_aa-2.0E0*1*I_KINETIC_G3xy_Px_a-2.0E0*2*I_KINETIC_Gxy2z_Px_a+2*1*I_KINETIC_Dxy_Px;
  abcd[95] = 4.0E0*I_KINETIC_I3x3z_Px_aa-2.0E0*2*I_KINETIC_G3xz_Px_a-2.0E0*2*I_KINETIC_Gx3z_Px_a+2*2*I_KINETIC_Dxz_Px;
  abcd[96] = 4.0E0*I_KINETIC_I2x3yz_Px_aa-2.0E0*1*I_KINETIC_G3yz_Px_a;
  abcd[97] = 4.0E0*I_KINETIC_I2x2y2z_Px_aa-2.0E0*1*I_KINETIC_G2x2y_Px_a-2.0E0*1*I_KINETIC_G2y2z_Px_a+1*I_KINETIC_D2y_Px;
  abcd[98] = 4.0E0*I_KINETIC_I2xy3z_Px_aa-2.0E0*2*I_KINETIC_G2xyz_Px_a-2.0E0*1*I_KINETIC_Gy3z_Px_a+2*I_KINETIC_Dyz_Px;
  abcd[99] = 4.0E0*I_KINETIC_I2x4z_Px_aa-2.0E0*3*I_KINETIC_G2x2z_Px_a-2.0E0*1*I_KINETIC_G4z_Px_a+3*I_KINETIC_D2z_Px;
  abcd[100] = 4.0E0*I_KINETIC_Ix4yz_Px_aa;
  abcd[101] = 4.0E0*I_KINETIC_Ix3y2z_Px_aa-2.0E0*1*I_KINETIC_Gx3y_Px_a;
  abcd[102] = 4.0E0*I_KINETIC_Ix2y3z_Px_aa-2.0E0*2*I_KINETIC_Gx2yz_Px_a;
  abcd[103] = 4.0E0*I_KINETIC_Ixy4z_Px_aa-2.0E0*3*I_KINETIC_Gxy2z_Px_a;
  abcd[104] = 4.0E0*I_KINETIC_Ix5z_Px_aa-2.0E0*4*I_KINETIC_Gx3z_Px_a;
  abcd[105] = 4.0E0*I_KINETIC_I5xz_Py_aa-2.0E0*4*I_KINETIC_G3xz_Py_a;
  abcd[106] = 4.0E0*I_KINETIC_I4xyz_Py_aa-2.0E0*3*I_KINETIC_G2xyz_Py_a;
  abcd[107] = 4.0E0*I_KINETIC_I4x2z_Py_aa-2.0E0*1*I_KINETIC_G4x_Py_a-2.0E0*3*I_KINETIC_G2x2z_Py_a+3*1*I_KINETIC_D2x_Py;
  abcd[108] = 4.0E0*I_KINETIC_I3x2yz_Py_aa-2.0E0*2*I_KINETIC_Gx2yz_Py_a;
  abcd[109] = 4.0E0*I_KINETIC_I3xy2z_Py_aa-2.0E0*1*I_KINETIC_G3xy_Py_a-2.0E0*2*I_KINETIC_Gxy2z_Py_a+2*1*I_KINETIC_Dxy_Py;
  abcd[110] = 4.0E0*I_KINETIC_I3x3z_Py_aa-2.0E0*2*I_KINETIC_G3xz_Py_a-2.0E0*2*I_KINETIC_Gx3z_Py_a+2*2*I_KINETIC_Dxz_Py;
  abcd[111] = 4.0E0*I_KINETIC_I2x3yz_Py_aa-2.0E0*1*I_KINETIC_G3yz_Py_a;
  abcd[112] = 4.0E0*I_KINETIC_I2x2y2z_Py_aa-2.0E0*1*I_KINETIC_G2x2y_Py_a-2.0E0*1*I_KINETIC_G2y2z_Py_a+1*I_KINETIC_D2y_Py;
  abcd[113] = 4.0E0*I_KINETIC_I2xy3z_Py_aa-2.0E0*2*I_KINETIC_G2xyz_Py_a-2.0E0*1*I_KINETIC_Gy3z_Py_a+2*I_KINETIC_Dyz_Py;
  abcd[114] = 4.0E0*I_KINETIC_I2x4z_Py_aa-2.0E0*3*I_KINETIC_G2x2z_Py_a-2.0E0*1*I_KINETIC_G4z_Py_a+3*I_KINETIC_D2z_Py;
  abcd[115] = 4.0E0*I_KINETIC_Ix4yz_Py_aa;
  abcd[116] = 4.0E0*I_KINETIC_Ix3y2z_Py_aa-2.0E0*1*I_KINETIC_Gx3y_Py_a;
  abcd[117] = 4.0E0*I_KINETIC_Ix2y3z_Py_aa-2.0E0*2*I_KINETIC_Gx2yz_Py_a;
  abcd[118] = 4.0E0*I_KINETIC_Ixy4z_Py_aa-2.0E0*3*I_KINETIC_Gxy2z_Py_a;
  abcd[119] = 4.0E0*I_KINETIC_Ix5z_Py_aa-2.0E0*4*I_KINETIC_Gx3z_Py_a;
  abcd[120] = 4.0E0*I_KINETIC_I5xz_Pz_aa-2.0E0*4*I_KINETIC_G3xz_Pz_a;
  abcd[121] = 4.0E0*I_KINETIC_I4xyz_Pz_aa-2.0E0*3*I_KINETIC_G2xyz_Pz_a;
  abcd[122] = 4.0E0*I_KINETIC_I4x2z_Pz_aa-2.0E0*1*I_KINETIC_G4x_Pz_a-2.0E0*3*I_KINETIC_G2x2z_Pz_a+3*1*I_KINETIC_D2x_Pz;
  abcd[123] = 4.0E0*I_KINETIC_I3x2yz_Pz_aa-2.0E0*2*I_KINETIC_Gx2yz_Pz_a;
  abcd[124] = 4.0E0*I_KINETIC_I3xy2z_Pz_aa-2.0E0*1*I_KINETIC_G3xy_Pz_a-2.0E0*2*I_KINETIC_Gxy2z_Pz_a+2*1*I_KINETIC_Dxy_Pz;
  abcd[125] = 4.0E0*I_KINETIC_I3x3z_Pz_aa-2.0E0*2*I_KINETIC_G3xz_Pz_a-2.0E0*2*I_KINETIC_Gx3z_Pz_a+2*2*I_KINETIC_Dxz_Pz;
  abcd[126] = 4.0E0*I_KINETIC_I2x3yz_Pz_aa-2.0E0*1*I_KINETIC_G3yz_Pz_a;
  abcd[127] = 4.0E0*I_KINETIC_I2x2y2z_Pz_aa-2.0E0*1*I_KINETIC_G2x2y_Pz_a-2.0E0*1*I_KINETIC_G2y2z_Pz_a+1*I_KINETIC_D2y_Pz;
  abcd[128] = 4.0E0*I_KINETIC_I2xy3z_Pz_aa-2.0E0*2*I_KINETIC_G2xyz_Pz_a-2.0E0*1*I_KINETIC_Gy3z_Pz_a+2*I_KINETIC_Dyz_Pz;
  abcd[129] = 4.0E0*I_KINETIC_I2x4z_Pz_aa-2.0E0*3*I_KINETIC_G2x2z_Pz_a-2.0E0*1*I_KINETIC_G4z_Pz_a+3*I_KINETIC_D2z_Pz;
  abcd[130] = 4.0E0*I_KINETIC_Ix4yz_Pz_aa;
  abcd[131] = 4.0E0*I_KINETIC_Ix3y2z_Pz_aa-2.0E0*1*I_KINETIC_Gx3y_Pz_a;
  abcd[132] = 4.0E0*I_KINETIC_Ix2y3z_Pz_aa-2.0E0*2*I_KINETIC_Gx2yz_Pz_a;
  abcd[133] = 4.0E0*I_KINETIC_Ixy4z_Pz_aa-2.0E0*3*I_KINETIC_Gxy2z_Pz_a;
  abcd[134] = 4.0E0*I_KINETIC_Ix5z_Pz_aa-2.0E0*4*I_KINETIC_Gx3z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_P_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_P_aa
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_D_P
   ************************************************************/
  abcd[135] = 4.0E0*I_KINETIC_I4x2y_Px_aa-2.0E0*1*I_KINETIC_G4x_Px_a;
  abcd[136] = 4.0E0*I_KINETIC_I3x3y_Px_aa-2.0E0*1*I_KINETIC_G3xy_Px_a-2.0E0*2*I_KINETIC_G3xy_Px_a;
  abcd[137] = 4.0E0*I_KINETIC_I3x2yz_Px_aa-2.0E0*1*I_KINETIC_G3xz_Px_a;
  abcd[138] = 4.0E0*I_KINETIC_I2x4y_Px_aa-2.0E0*2*I_KINETIC_G2x2y_Px_a-2.0E0*3*I_KINETIC_G2x2y_Px_a+2*1*I_KINETIC_D2x_Px;
  abcd[139] = 4.0E0*I_KINETIC_I2x3yz_Px_aa-2.0E0*1*I_KINETIC_G2xyz_Px_a-2.0E0*2*I_KINETIC_G2xyz_Px_a;
  abcd[140] = 4.0E0*I_KINETIC_I2x2y2z_Px_aa-2.0E0*1*I_KINETIC_G2x2z_Px_a;
  abcd[141] = 4.0E0*I_KINETIC_Ix5y_Px_aa-2.0E0*3*I_KINETIC_Gx3y_Px_a-2.0E0*4*I_KINETIC_Gx3y_Px_a+3*2*I_KINETIC_Dxy_Px;
  abcd[142] = 4.0E0*I_KINETIC_Ix4yz_Px_aa-2.0E0*2*I_KINETIC_Gx2yz_Px_a-2.0E0*3*I_KINETIC_Gx2yz_Px_a+2*1*I_KINETIC_Dxz_Px;
  abcd[143] = 4.0E0*I_KINETIC_Ix3y2z_Px_aa-2.0E0*1*I_KINETIC_Gxy2z_Px_a-2.0E0*2*I_KINETIC_Gxy2z_Px_a;
  abcd[144] = 4.0E0*I_KINETIC_Ix2y3z_Px_aa-2.0E0*1*I_KINETIC_Gx3z_Px_a;
  abcd[145] = 4.0E0*I_KINETIC_I6y_Px_aa-2.0E0*4*I_KINETIC_G4y_Px_a-2.0E0*5*I_KINETIC_G4y_Px_a+4*3*I_KINETIC_D2y_Px;
  abcd[146] = 4.0E0*I_KINETIC_I5yz_Px_aa-2.0E0*3*I_KINETIC_G3yz_Px_a-2.0E0*4*I_KINETIC_G3yz_Px_a+3*2*I_KINETIC_Dyz_Px;
  abcd[147] = 4.0E0*I_KINETIC_I4y2z_Px_aa-2.0E0*2*I_KINETIC_G2y2z_Px_a-2.0E0*3*I_KINETIC_G2y2z_Px_a+2*1*I_KINETIC_D2z_Px;
  abcd[148] = 4.0E0*I_KINETIC_I3y3z_Px_aa-2.0E0*1*I_KINETIC_Gy3z_Px_a-2.0E0*2*I_KINETIC_Gy3z_Px_a;
  abcd[149] = 4.0E0*I_KINETIC_I2y4z_Px_aa-2.0E0*1*I_KINETIC_G4z_Px_a;
  abcd[150] = 4.0E0*I_KINETIC_I4x2y_Py_aa-2.0E0*1*I_KINETIC_G4x_Py_a;
  abcd[151] = 4.0E0*I_KINETIC_I3x3y_Py_aa-2.0E0*1*I_KINETIC_G3xy_Py_a-2.0E0*2*I_KINETIC_G3xy_Py_a;
  abcd[152] = 4.0E0*I_KINETIC_I3x2yz_Py_aa-2.0E0*1*I_KINETIC_G3xz_Py_a;
  abcd[153] = 4.0E0*I_KINETIC_I2x4y_Py_aa-2.0E0*2*I_KINETIC_G2x2y_Py_a-2.0E0*3*I_KINETIC_G2x2y_Py_a+2*1*I_KINETIC_D2x_Py;
  abcd[154] = 4.0E0*I_KINETIC_I2x3yz_Py_aa-2.0E0*1*I_KINETIC_G2xyz_Py_a-2.0E0*2*I_KINETIC_G2xyz_Py_a;
  abcd[155] = 4.0E0*I_KINETIC_I2x2y2z_Py_aa-2.0E0*1*I_KINETIC_G2x2z_Py_a;
  abcd[156] = 4.0E0*I_KINETIC_Ix5y_Py_aa-2.0E0*3*I_KINETIC_Gx3y_Py_a-2.0E0*4*I_KINETIC_Gx3y_Py_a+3*2*I_KINETIC_Dxy_Py;
  abcd[157] = 4.0E0*I_KINETIC_Ix4yz_Py_aa-2.0E0*2*I_KINETIC_Gx2yz_Py_a-2.0E0*3*I_KINETIC_Gx2yz_Py_a+2*1*I_KINETIC_Dxz_Py;
  abcd[158] = 4.0E0*I_KINETIC_Ix3y2z_Py_aa-2.0E0*1*I_KINETIC_Gxy2z_Py_a-2.0E0*2*I_KINETIC_Gxy2z_Py_a;
  abcd[159] = 4.0E0*I_KINETIC_Ix2y3z_Py_aa-2.0E0*1*I_KINETIC_Gx3z_Py_a;
  abcd[160] = 4.0E0*I_KINETIC_I6y_Py_aa-2.0E0*4*I_KINETIC_G4y_Py_a-2.0E0*5*I_KINETIC_G4y_Py_a+4*3*I_KINETIC_D2y_Py;
  abcd[161] = 4.0E0*I_KINETIC_I5yz_Py_aa-2.0E0*3*I_KINETIC_G3yz_Py_a-2.0E0*4*I_KINETIC_G3yz_Py_a+3*2*I_KINETIC_Dyz_Py;
  abcd[162] = 4.0E0*I_KINETIC_I4y2z_Py_aa-2.0E0*2*I_KINETIC_G2y2z_Py_a-2.0E0*3*I_KINETIC_G2y2z_Py_a+2*1*I_KINETIC_D2z_Py;
  abcd[163] = 4.0E0*I_KINETIC_I3y3z_Py_aa-2.0E0*1*I_KINETIC_Gy3z_Py_a-2.0E0*2*I_KINETIC_Gy3z_Py_a;
  abcd[164] = 4.0E0*I_KINETIC_I2y4z_Py_aa-2.0E0*1*I_KINETIC_G4z_Py_a;
  abcd[165] = 4.0E0*I_KINETIC_I4x2y_Pz_aa-2.0E0*1*I_KINETIC_G4x_Pz_a;
  abcd[166] = 4.0E0*I_KINETIC_I3x3y_Pz_aa-2.0E0*1*I_KINETIC_G3xy_Pz_a-2.0E0*2*I_KINETIC_G3xy_Pz_a;
  abcd[167] = 4.0E0*I_KINETIC_I3x2yz_Pz_aa-2.0E0*1*I_KINETIC_G3xz_Pz_a;
  abcd[168] = 4.0E0*I_KINETIC_I2x4y_Pz_aa-2.0E0*2*I_KINETIC_G2x2y_Pz_a-2.0E0*3*I_KINETIC_G2x2y_Pz_a+2*1*I_KINETIC_D2x_Pz;
  abcd[169] = 4.0E0*I_KINETIC_I2x3yz_Pz_aa-2.0E0*1*I_KINETIC_G2xyz_Pz_a-2.0E0*2*I_KINETIC_G2xyz_Pz_a;
  abcd[170] = 4.0E0*I_KINETIC_I2x2y2z_Pz_aa-2.0E0*1*I_KINETIC_G2x2z_Pz_a;
  abcd[171] = 4.0E0*I_KINETIC_Ix5y_Pz_aa-2.0E0*3*I_KINETIC_Gx3y_Pz_a-2.0E0*4*I_KINETIC_Gx3y_Pz_a+3*2*I_KINETIC_Dxy_Pz;
  abcd[172] = 4.0E0*I_KINETIC_Ix4yz_Pz_aa-2.0E0*2*I_KINETIC_Gx2yz_Pz_a-2.0E0*3*I_KINETIC_Gx2yz_Pz_a+2*1*I_KINETIC_Dxz_Pz;
  abcd[173] = 4.0E0*I_KINETIC_Ix3y2z_Pz_aa-2.0E0*1*I_KINETIC_Gxy2z_Pz_a-2.0E0*2*I_KINETIC_Gxy2z_Pz_a;
  abcd[174] = 4.0E0*I_KINETIC_Ix2y3z_Pz_aa-2.0E0*1*I_KINETIC_Gx3z_Pz_a;
  abcd[175] = 4.0E0*I_KINETIC_I6y_Pz_aa-2.0E0*4*I_KINETIC_G4y_Pz_a-2.0E0*5*I_KINETIC_G4y_Pz_a+4*3*I_KINETIC_D2y_Pz;
  abcd[176] = 4.0E0*I_KINETIC_I5yz_Pz_aa-2.0E0*3*I_KINETIC_G3yz_Pz_a-2.0E0*4*I_KINETIC_G3yz_Pz_a+3*2*I_KINETIC_Dyz_Pz;
  abcd[177] = 4.0E0*I_KINETIC_I4y2z_Pz_aa-2.0E0*2*I_KINETIC_G2y2z_Pz_a-2.0E0*3*I_KINETIC_G2y2z_Pz_a+2*1*I_KINETIC_D2z_Pz;
  abcd[178] = 4.0E0*I_KINETIC_I3y3z_Pz_aa-2.0E0*1*I_KINETIC_Gy3z_Pz_a-2.0E0*2*I_KINETIC_Gy3z_Pz_a;
  abcd[179] = 4.0E0*I_KINETIC_I2y4z_Pz_aa-2.0E0*1*I_KINETIC_G4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_P_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_P_aa
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_D_P
   ************************************************************/
  abcd[180] = 4.0E0*I_KINETIC_I4xyz_Px_aa;
  abcd[181] = 4.0E0*I_KINETIC_I3x2yz_Px_aa-2.0E0*1*I_KINETIC_G3xz_Px_a;
  abcd[182] = 4.0E0*I_KINETIC_I3xy2z_Px_aa-2.0E0*1*I_KINETIC_G3xy_Px_a;
  abcd[183] = 4.0E0*I_KINETIC_I2x3yz_Px_aa-2.0E0*2*I_KINETIC_G2xyz_Px_a;
  abcd[184] = 4.0E0*I_KINETIC_I2x2y2z_Px_aa-2.0E0*1*I_KINETIC_G2x2y_Px_a-2.0E0*1*I_KINETIC_G2x2z_Px_a+1*I_KINETIC_D2x_Px;
  abcd[185] = 4.0E0*I_KINETIC_I2xy3z_Px_aa-2.0E0*2*I_KINETIC_G2xyz_Px_a;
  abcd[186] = 4.0E0*I_KINETIC_Ix4yz_Px_aa-2.0E0*3*I_KINETIC_Gx2yz_Px_a;
  abcd[187] = 4.0E0*I_KINETIC_Ix3y2z_Px_aa-2.0E0*1*I_KINETIC_Gx3y_Px_a-2.0E0*2*I_KINETIC_Gxy2z_Px_a+2*1*I_KINETIC_Dxy_Px;
  abcd[188] = 4.0E0*I_KINETIC_Ix2y3z_Px_aa-2.0E0*2*I_KINETIC_Gx2yz_Px_a-2.0E0*1*I_KINETIC_Gx3z_Px_a+2*I_KINETIC_Dxz_Px;
  abcd[189] = 4.0E0*I_KINETIC_Ixy4z_Px_aa-2.0E0*3*I_KINETIC_Gxy2z_Px_a;
  abcd[190] = 4.0E0*I_KINETIC_I5yz_Px_aa-2.0E0*4*I_KINETIC_G3yz_Px_a;
  abcd[191] = 4.0E0*I_KINETIC_I4y2z_Px_aa-2.0E0*1*I_KINETIC_G4y_Px_a-2.0E0*3*I_KINETIC_G2y2z_Px_a+3*1*I_KINETIC_D2y_Px;
  abcd[192] = 4.0E0*I_KINETIC_I3y3z_Px_aa-2.0E0*2*I_KINETIC_G3yz_Px_a-2.0E0*2*I_KINETIC_Gy3z_Px_a+2*2*I_KINETIC_Dyz_Px;
  abcd[193] = 4.0E0*I_KINETIC_I2y4z_Px_aa-2.0E0*3*I_KINETIC_G2y2z_Px_a-2.0E0*1*I_KINETIC_G4z_Px_a+3*I_KINETIC_D2z_Px;
  abcd[194] = 4.0E0*I_KINETIC_Iy5z_Px_aa-2.0E0*4*I_KINETIC_Gy3z_Px_a;
  abcd[195] = 4.0E0*I_KINETIC_I4xyz_Py_aa;
  abcd[196] = 4.0E0*I_KINETIC_I3x2yz_Py_aa-2.0E0*1*I_KINETIC_G3xz_Py_a;
  abcd[197] = 4.0E0*I_KINETIC_I3xy2z_Py_aa-2.0E0*1*I_KINETIC_G3xy_Py_a;
  abcd[198] = 4.0E0*I_KINETIC_I2x3yz_Py_aa-2.0E0*2*I_KINETIC_G2xyz_Py_a;
  abcd[199] = 4.0E0*I_KINETIC_I2x2y2z_Py_aa-2.0E0*1*I_KINETIC_G2x2y_Py_a-2.0E0*1*I_KINETIC_G2x2z_Py_a+1*I_KINETIC_D2x_Py;
  abcd[200] = 4.0E0*I_KINETIC_I2xy3z_Py_aa-2.0E0*2*I_KINETIC_G2xyz_Py_a;
  abcd[201] = 4.0E0*I_KINETIC_Ix4yz_Py_aa-2.0E0*3*I_KINETIC_Gx2yz_Py_a;
  abcd[202] = 4.0E0*I_KINETIC_Ix3y2z_Py_aa-2.0E0*1*I_KINETIC_Gx3y_Py_a-2.0E0*2*I_KINETIC_Gxy2z_Py_a+2*1*I_KINETIC_Dxy_Py;
  abcd[203] = 4.0E0*I_KINETIC_Ix2y3z_Py_aa-2.0E0*2*I_KINETIC_Gx2yz_Py_a-2.0E0*1*I_KINETIC_Gx3z_Py_a+2*I_KINETIC_Dxz_Py;
  abcd[204] = 4.0E0*I_KINETIC_Ixy4z_Py_aa-2.0E0*3*I_KINETIC_Gxy2z_Py_a;
  abcd[205] = 4.0E0*I_KINETIC_I5yz_Py_aa-2.0E0*4*I_KINETIC_G3yz_Py_a;
  abcd[206] = 4.0E0*I_KINETIC_I4y2z_Py_aa-2.0E0*1*I_KINETIC_G4y_Py_a-2.0E0*3*I_KINETIC_G2y2z_Py_a+3*1*I_KINETIC_D2y_Py;
  abcd[207] = 4.0E0*I_KINETIC_I3y3z_Py_aa-2.0E0*2*I_KINETIC_G3yz_Py_a-2.0E0*2*I_KINETIC_Gy3z_Py_a+2*2*I_KINETIC_Dyz_Py;
  abcd[208] = 4.0E0*I_KINETIC_I2y4z_Py_aa-2.0E0*3*I_KINETIC_G2y2z_Py_a-2.0E0*1*I_KINETIC_G4z_Py_a+3*I_KINETIC_D2z_Py;
  abcd[209] = 4.0E0*I_KINETIC_Iy5z_Py_aa-2.0E0*4*I_KINETIC_Gy3z_Py_a;
  abcd[210] = 4.0E0*I_KINETIC_I4xyz_Pz_aa;
  abcd[211] = 4.0E0*I_KINETIC_I3x2yz_Pz_aa-2.0E0*1*I_KINETIC_G3xz_Pz_a;
  abcd[212] = 4.0E0*I_KINETIC_I3xy2z_Pz_aa-2.0E0*1*I_KINETIC_G3xy_Pz_a;
  abcd[213] = 4.0E0*I_KINETIC_I2x3yz_Pz_aa-2.0E0*2*I_KINETIC_G2xyz_Pz_a;
  abcd[214] = 4.0E0*I_KINETIC_I2x2y2z_Pz_aa-2.0E0*1*I_KINETIC_G2x2y_Pz_a-2.0E0*1*I_KINETIC_G2x2z_Pz_a+1*I_KINETIC_D2x_Pz;
  abcd[215] = 4.0E0*I_KINETIC_I2xy3z_Pz_aa-2.0E0*2*I_KINETIC_G2xyz_Pz_a;
  abcd[216] = 4.0E0*I_KINETIC_Ix4yz_Pz_aa-2.0E0*3*I_KINETIC_Gx2yz_Pz_a;
  abcd[217] = 4.0E0*I_KINETIC_Ix3y2z_Pz_aa-2.0E0*1*I_KINETIC_Gx3y_Pz_a-2.0E0*2*I_KINETIC_Gxy2z_Pz_a+2*1*I_KINETIC_Dxy_Pz;
  abcd[218] = 4.0E0*I_KINETIC_Ix2y3z_Pz_aa-2.0E0*2*I_KINETIC_Gx2yz_Pz_a-2.0E0*1*I_KINETIC_Gx3z_Pz_a+2*I_KINETIC_Dxz_Pz;
  abcd[219] = 4.0E0*I_KINETIC_Ixy4z_Pz_aa-2.0E0*3*I_KINETIC_Gxy2z_Pz_a;
  abcd[220] = 4.0E0*I_KINETIC_I5yz_Pz_aa-2.0E0*4*I_KINETIC_G3yz_Pz_a;
  abcd[221] = 4.0E0*I_KINETIC_I4y2z_Pz_aa-2.0E0*1*I_KINETIC_G4y_Pz_a-2.0E0*3*I_KINETIC_G2y2z_Pz_a+3*1*I_KINETIC_D2y_Pz;
  abcd[222] = 4.0E0*I_KINETIC_I3y3z_Pz_aa-2.0E0*2*I_KINETIC_G3yz_Pz_a-2.0E0*2*I_KINETIC_Gy3z_Pz_a+2*2*I_KINETIC_Dyz_Pz;
  abcd[223] = 4.0E0*I_KINETIC_I2y4z_Pz_aa-2.0E0*3*I_KINETIC_G2y2z_Pz_a-2.0E0*1*I_KINETIC_G4z_Pz_a+3*I_KINETIC_D2z_Pz;
  abcd[224] = 4.0E0*I_KINETIC_Iy5z_Pz_aa-2.0E0*4*I_KINETIC_Gy3z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_G_P_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_P_aa
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_G_P_a
   * RHS shell quartet name: SQ_KINETIC_D_P
   ************************************************************/
  abcd[225] = 4.0E0*I_KINETIC_I4x2z_Px_aa-2.0E0*1*I_KINETIC_G4x_Px_a;
  abcd[226] = 4.0E0*I_KINETIC_I3xy2z_Px_aa-2.0E0*1*I_KINETIC_G3xy_Px_a;
  abcd[227] = 4.0E0*I_KINETIC_I3x3z_Px_aa-2.0E0*1*I_KINETIC_G3xz_Px_a-2.0E0*2*I_KINETIC_G3xz_Px_a;
  abcd[228] = 4.0E0*I_KINETIC_I2x2y2z_Px_aa-2.0E0*1*I_KINETIC_G2x2y_Px_a;
  abcd[229] = 4.0E0*I_KINETIC_I2xy3z_Px_aa-2.0E0*1*I_KINETIC_G2xyz_Px_a-2.0E0*2*I_KINETIC_G2xyz_Px_a;
  abcd[230] = 4.0E0*I_KINETIC_I2x4z_Px_aa-2.0E0*2*I_KINETIC_G2x2z_Px_a-2.0E0*3*I_KINETIC_G2x2z_Px_a+2*1*I_KINETIC_D2x_Px;
  abcd[231] = 4.0E0*I_KINETIC_Ix3y2z_Px_aa-2.0E0*1*I_KINETIC_Gx3y_Px_a;
  abcd[232] = 4.0E0*I_KINETIC_Ix2y3z_Px_aa-2.0E0*1*I_KINETIC_Gx2yz_Px_a-2.0E0*2*I_KINETIC_Gx2yz_Px_a;
  abcd[233] = 4.0E0*I_KINETIC_Ixy4z_Px_aa-2.0E0*2*I_KINETIC_Gxy2z_Px_a-2.0E0*3*I_KINETIC_Gxy2z_Px_a+2*1*I_KINETIC_Dxy_Px;
  abcd[234] = 4.0E0*I_KINETIC_Ix5z_Px_aa-2.0E0*3*I_KINETIC_Gx3z_Px_a-2.0E0*4*I_KINETIC_Gx3z_Px_a+3*2*I_KINETIC_Dxz_Px;
  abcd[235] = 4.0E0*I_KINETIC_I4y2z_Px_aa-2.0E0*1*I_KINETIC_G4y_Px_a;
  abcd[236] = 4.0E0*I_KINETIC_I3y3z_Px_aa-2.0E0*1*I_KINETIC_G3yz_Px_a-2.0E0*2*I_KINETIC_G3yz_Px_a;
  abcd[237] = 4.0E0*I_KINETIC_I2y4z_Px_aa-2.0E0*2*I_KINETIC_G2y2z_Px_a-2.0E0*3*I_KINETIC_G2y2z_Px_a+2*1*I_KINETIC_D2y_Px;
  abcd[238] = 4.0E0*I_KINETIC_Iy5z_Px_aa-2.0E0*3*I_KINETIC_Gy3z_Px_a-2.0E0*4*I_KINETIC_Gy3z_Px_a+3*2*I_KINETIC_Dyz_Px;
  abcd[239] = 4.0E0*I_KINETIC_I6z_Px_aa-2.0E0*4*I_KINETIC_G4z_Px_a-2.0E0*5*I_KINETIC_G4z_Px_a+4*3*I_KINETIC_D2z_Px;
  abcd[240] = 4.0E0*I_KINETIC_I4x2z_Py_aa-2.0E0*1*I_KINETIC_G4x_Py_a;
  abcd[241] = 4.0E0*I_KINETIC_I3xy2z_Py_aa-2.0E0*1*I_KINETIC_G3xy_Py_a;
  abcd[242] = 4.0E0*I_KINETIC_I3x3z_Py_aa-2.0E0*1*I_KINETIC_G3xz_Py_a-2.0E0*2*I_KINETIC_G3xz_Py_a;
  abcd[243] = 4.0E0*I_KINETIC_I2x2y2z_Py_aa-2.0E0*1*I_KINETIC_G2x2y_Py_a;
  abcd[244] = 4.0E0*I_KINETIC_I2xy3z_Py_aa-2.0E0*1*I_KINETIC_G2xyz_Py_a-2.0E0*2*I_KINETIC_G2xyz_Py_a;
  abcd[245] = 4.0E0*I_KINETIC_I2x4z_Py_aa-2.0E0*2*I_KINETIC_G2x2z_Py_a-2.0E0*3*I_KINETIC_G2x2z_Py_a+2*1*I_KINETIC_D2x_Py;
  abcd[246] = 4.0E0*I_KINETIC_Ix3y2z_Py_aa-2.0E0*1*I_KINETIC_Gx3y_Py_a;
  abcd[247] = 4.0E0*I_KINETIC_Ix2y3z_Py_aa-2.0E0*1*I_KINETIC_Gx2yz_Py_a-2.0E0*2*I_KINETIC_Gx2yz_Py_a;
  abcd[248] = 4.0E0*I_KINETIC_Ixy4z_Py_aa-2.0E0*2*I_KINETIC_Gxy2z_Py_a-2.0E0*3*I_KINETIC_Gxy2z_Py_a+2*1*I_KINETIC_Dxy_Py;
  abcd[249] = 4.0E0*I_KINETIC_Ix5z_Py_aa-2.0E0*3*I_KINETIC_Gx3z_Py_a-2.0E0*4*I_KINETIC_Gx3z_Py_a+3*2*I_KINETIC_Dxz_Py;
  abcd[250] = 4.0E0*I_KINETIC_I4y2z_Py_aa-2.0E0*1*I_KINETIC_G4y_Py_a;
  abcd[251] = 4.0E0*I_KINETIC_I3y3z_Py_aa-2.0E0*1*I_KINETIC_G3yz_Py_a-2.0E0*2*I_KINETIC_G3yz_Py_a;
  abcd[252] = 4.0E0*I_KINETIC_I2y4z_Py_aa-2.0E0*2*I_KINETIC_G2y2z_Py_a-2.0E0*3*I_KINETIC_G2y2z_Py_a+2*1*I_KINETIC_D2y_Py;
  abcd[253] = 4.0E0*I_KINETIC_Iy5z_Py_aa-2.0E0*3*I_KINETIC_Gy3z_Py_a-2.0E0*4*I_KINETIC_Gy3z_Py_a+3*2*I_KINETIC_Dyz_Py;
  abcd[254] = 4.0E0*I_KINETIC_I6z_Py_aa-2.0E0*4*I_KINETIC_G4z_Py_a-2.0E0*5*I_KINETIC_G4z_Py_a+4*3*I_KINETIC_D2z_Py;
  abcd[255] = 4.0E0*I_KINETIC_I4x2z_Pz_aa-2.0E0*1*I_KINETIC_G4x_Pz_a;
  abcd[256] = 4.0E0*I_KINETIC_I3xy2z_Pz_aa-2.0E0*1*I_KINETIC_G3xy_Pz_a;
  abcd[257] = 4.0E0*I_KINETIC_I3x3z_Pz_aa-2.0E0*1*I_KINETIC_G3xz_Pz_a-2.0E0*2*I_KINETIC_G3xz_Pz_a;
  abcd[258] = 4.0E0*I_KINETIC_I2x2y2z_Pz_aa-2.0E0*1*I_KINETIC_G2x2y_Pz_a;
  abcd[259] = 4.0E0*I_KINETIC_I2xy3z_Pz_aa-2.0E0*1*I_KINETIC_G2xyz_Pz_a-2.0E0*2*I_KINETIC_G2xyz_Pz_a;
  abcd[260] = 4.0E0*I_KINETIC_I2x4z_Pz_aa-2.0E0*2*I_KINETIC_G2x2z_Pz_a-2.0E0*3*I_KINETIC_G2x2z_Pz_a+2*1*I_KINETIC_D2x_Pz;
  abcd[261] = 4.0E0*I_KINETIC_Ix3y2z_Pz_aa-2.0E0*1*I_KINETIC_Gx3y_Pz_a;
  abcd[262] = 4.0E0*I_KINETIC_Ix2y3z_Pz_aa-2.0E0*1*I_KINETIC_Gx2yz_Pz_a-2.0E0*2*I_KINETIC_Gx2yz_Pz_a;
  abcd[263] = 4.0E0*I_KINETIC_Ixy4z_Pz_aa-2.0E0*2*I_KINETIC_Gxy2z_Pz_a-2.0E0*3*I_KINETIC_Gxy2z_Pz_a+2*1*I_KINETIC_Dxy_Pz;
  abcd[264] = 4.0E0*I_KINETIC_Ix5z_Pz_aa-2.0E0*3*I_KINETIC_Gx3z_Pz_a-2.0E0*4*I_KINETIC_Gx3z_Pz_a+3*2*I_KINETIC_Dxz_Pz;
  abcd[265] = 4.0E0*I_KINETIC_I4y2z_Pz_aa-2.0E0*1*I_KINETIC_G4y_Pz_a;
  abcd[266] = 4.0E0*I_KINETIC_I3y3z_Pz_aa-2.0E0*1*I_KINETIC_G3yz_Pz_a-2.0E0*2*I_KINETIC_G3yz_Pz_a;
  abcd[267] = 4.0E0*I_KINETIC_I2y4z_Pz_aa-2.0E0*2*I_KINETIC_G2y2z_Pz_a-2.0E0*3*I_KINETIC_G2y2z_Pz_a+2*1*I_KINETIC_D2y_Pz;
  abcd[268] = 4.0E0*I_KINETIC_Iy5z_Pz_aa-2.0E0*3*I_KINETIC_Gy3z_Pz_a-2.0E0*4*I_KINETIC_Gy3z_Pz_a+3*2*I_KINETIC_Dyz_Pz;
  abcd[269] = 4.0E0*I_KINETIC_I6z_Pz_aa-2.0E0*4*I_KINETIC_G4z_Pz_a-2.0E0*5*I_KINETIC_G4z_Pz_a+4*3*I_KINETIC_D2z_Pz;
}
