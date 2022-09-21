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

void hgp_os_kinetic_h_g(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
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
     * totally 30 integrals are omitted 
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
     * shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 45 integrals are omitted 
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
     * shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_Px_vrr = PAX*I_TWOBODYOVERLAP_G4x_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Px_vrr = PAY*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Px_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Px_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Px_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Px_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Px_vrr = PAX*I_TWOBODYOVERLAP_G4y_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Px_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Px_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Px_vrr = PAX*I_TWOBODYOVERLAP_G4z_Px_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_Px_vrr = PAY*I_TWOBODYOVERLAP_G4y_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Px_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Px_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Px_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Px_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Px_vrr = PAY*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_Px_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Px_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_Py_vrr = PAX*I_TWOBODYOVERLAP_G4x_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Py_vrr = PAY*I_TWOBODYOVERLAP_G4x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Py_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Py_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Py_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Py_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Py_vrr = PAX*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Py_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Py_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Py_vrr = PAX*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5y_Py_vrr = PAY*I_TWOBODYOVERLAP_G4y_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Py_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Py_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Py_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Py_vrr = PAY*I_TWOBODYOVERLAP_G4z_Py_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_Py_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Py_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5x_Pz_vrr = PAX*I_TWOBODYOVERLAP_G4x_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Pz_vrr = PAY*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G4x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Pz_vrr = PAY*I_TWOBODYOVERLAP_G3xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3x_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Pz_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Pz_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Pz_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Pz_vrr = PAX*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Pz_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Pz_vrr = PAX*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5y_Pz_vrr = PAY*I_TWOBODYOVERLAP_G4y_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G4y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3y_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Pz_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Pz_vrr = PAY*I_TWOBODYOVERLAP_G4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5z_Pz_vrr = PAZ*I_TWOBODYOVERLAP_G4z_Pz_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_Pz_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 30 integrals are omitted 
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
     * shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * expanding position: BRA1
     * code section is: VRR
     * totally 63 integrals are omitted 
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
     * shell quartet name: SQ_TWOBODYOVERLAP_H_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 42 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_F3x_vrr = PBX*I_TWOBODYOVERLAP_H5x_D2x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F3x_vrr = PBX*I_TWOBODYOVERLAP_H4xy_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F3x_vrr = PBX*I_TWOBODYOVERLAP_H4xz_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F3x_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F3x_vrr = PBX*I_TWOBODYOVERLAP_H3xyz_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2xyz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H3xyz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F3x_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2xy2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F3x_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx3yz_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx2y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hxy3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5y_F3x_vrr = PBX*I_TWOBODYOVERLAP_H5y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F3x_vrr = PBX*I_TWOBODYOVERLAP_H4yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F3x_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5z_F3x_vrr = PBX*I_TWOBODYOVERLAP_H5z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_TWOBODYOVERLAP_H5x_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H4xy_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5y_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H5y_D2x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H4yz_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_F2xy_vrr = PBY*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5y_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_D2x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_D2x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_D2x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_D2x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_F2xz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_D2x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H5x_D2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H4xy_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H4xz_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H3xyz_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2xyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fx2y_vrr = PBX*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H5x_D2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H4xy_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H4xz_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H3xyz_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2xyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H2xy2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Hx3yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Hx2y2z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Hxy3z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H4yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Fx2z_vrr = PBX*I_TWOBODYOVERLAP_H5z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5x_F3y_vrr = PBY*I_TWOBODYOVERLAP_H5x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F3y_vrr = PBY*I_TWOBODYOVERLAP_H4xy_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F3y_vrr = PBY*I_TWOBODYOVERLAP_H4xz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F3y_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F3y_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3xyz_Py_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F3y_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_Py_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2xy2z_Py_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F3y_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx3yz_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx2y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hxy3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5y_F3y_vrr = PBY*I_TWOBODYOVERLAP_H5y_D2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F3y_vrr = PBY*I_TWOBODYOVERLAP_H4yz_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H2y3z_Py_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F3y_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5z_F3y_vrr = PBY*I_TWOBODYOVERLAP_H5z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_TWOBODYOVERLAP_H5x_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_D2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_D2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_D2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_D2y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5z_F2yz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_D2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H5x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xy_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4xz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3xyz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3x2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2xy2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx3yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hxy3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5y_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H5y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_TWOBODYOVERLAP_H4yz_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_D2z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_Pz_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H3y2z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_D2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H2y3z_Pz_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_D2z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_Pz_vrr;
    Double I_TWOBODYOVERLAP_H5z_F3z_vrr = PBZ*I_TWOBODYOVERLAP_H5z_D2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_D2z_vrr+2*oned2z*I_TWOBODYOVERLAP_H5z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_G
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_G4x_vrr = PBX*I_TWOBODYOVERLAP_H5x_F3x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G4x_vrr = PBX*I_TWOBODYOVERLAP_H4xy_F3x_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G4x_vrr = PBX*I_TWOBODYOVERLAP_H4xz_F3x_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G4x_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G4x_vrr = PBX*I_TWOBODYOVERLAP_H3xyz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2xyz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H3xyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G4x_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G4x_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G4x_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G4x_vrr = PBX*I_TWOBODYOVERLAP_H2xy2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G4x_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G4x_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G4x_vrr = PBX*I_TWOBODYOVERLAP_Hx3yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G4x_vrr = PBX*I_TWOBODYOVERLAP_Hx2y2z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G4x_vrr = PBX*I_TWOBODYOVERLAP_Hxy3z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G4x_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5y_G4x_vrr = PBX*I_TWOBODYOVERLAP_H5y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G4x_vrr = PBX*I_TWOBODYOVERLAP_H4yz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G4x_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G4x_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H2y3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G4x_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_G4x_vrr = PBX*I_TWOBODYOVERLAP_H5z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H5x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H4xy_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H4xz_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_F3x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5y_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H5y_F3x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H4yz_F3x_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5z_G3xy_vrr = PBY*I_TWOBODYOVERLAP_H5z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5x_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_F3x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5y_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_F3x_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_F3x_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F3x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_F3x_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F3x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_F3x_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5z_G3xz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_F3x_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F3x_vrr;
    Double I_TWOBODYOVERLAP_H5x_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H5x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H4xy_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H4xz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H3xyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_F2xy_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5y_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H5y_F2xy_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H4yz_F2xy_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_G2x2y_vrr = PBY*I_TWOBODYOVERLAP_H5z_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_F2xy_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H5y_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_F2xy_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_F2xy_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_F2xy_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_F2xy_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H5z_G2xyz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_F2xy_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F2xy_vrr;
    Double I_TWOBODYOVERLAP_H5x_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H5x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H3xyz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_F2xz_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5y_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H5y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_F2xz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_F2xz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_D2x_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_F2xz_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5z_G2x2z_vrr = PBZ*I_TWOBODYOVERLAP_H5z_F2xz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F2xz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H5x_F3y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H4xy_F3y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H4xz_F3y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H3xyz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2xyz_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H2xy2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Hx3yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Hx2y2z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Hxy3z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H5y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H4yz_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gx3y_vrr = PBX*I_TWOBODYOVERLAP_H5z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_Fx2y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_Fx2y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_Fx2y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_Fx2y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_Fx2y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gx2yz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_Fx2y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_Fx2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H5x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H4xy_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H4xz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_Fx2z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H5y_Fx2z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H4yz_Fx2z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_Fx2z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_Fx2z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_Fx2z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gxy2z_vrr = PBY*I_TWOBODYOVERLAP_H5z_Fx2z_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H5x_F3z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H4xy_F3z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H4xz_F3z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H3xyz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2xyz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H2xy2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Hx3yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Hx2y2z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Hxy3z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H5y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H4yz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gx3z_vrr = PBX*I_TWOBODYOVERLAP_H5z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5x_G4y_vrr = PBY*I_TWOBODYOVERLAP_H5x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G4y_vrr = PBY*I_TWOBODYOVERLAP_H4xy_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G4y_vrr = PBY*I_TWOBODYOVERLAP_H4xz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H4xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G4y_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G4y_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H3xyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G4y_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G4y_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G4y_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G4y_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G4y_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G4y_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_F3y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G4y_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G4y_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G4y_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G4y_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_G4y_vrr = PBY*I_TWOBODYOVERLAP_H5y_F3y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G4y_vrr = PBY*I_TWOBODYOVERLAP_H4yz_F3y_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G4y_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G4y_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H2y3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G4y_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5z_G4y_vrr = PBY*I_TWOBODYOVERLAP_H5z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_F3y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5y_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_F3y_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_F3y_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F3y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_F3y_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F3y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_F3y_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5z_G3yz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_F3y_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F3y_vrr;
    Double I_TWOBODYOVERLAP_H5x_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H5x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H3xyz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_F2yz_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5y_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H5y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_F2yz_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_F2yz_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_D2y_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_F2yz_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5z_G2y2z_vrr = PBZ*I_TWOBODYOVERLAP_H5z_F2yz_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F2yz_vrr+oned2z*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_TWOBODYOVERLAP_H5x_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H5x_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H4xy_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H4xz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_F3z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5y_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H5y_F3z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_F3z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H4yz_F3z_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_F3z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_F3z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5z_Gy3z_vrr = PBY*I_TWOBODYOVERLAP_H5z_F3z_vrr;
    Double I_TWOBODYOVERLAP_H5x_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H5x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xy_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H4xy_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4xz_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4x_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H4xz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H3xyz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4y_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx3yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx2y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Hxy3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_F3z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5y_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H5y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_TWOBODYOVERLAP_H4yz_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_F3z_vrr+oned2z*I_TWOBODYOVERLAP_G4y_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H4yz_D2z_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_F3z_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H2y3z_D2z_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_F3z_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_Hy4z_D2z_vrr;
    Double I_TWOBODYOVERLAP_H5z_G4z_vrr = PBZ*I_TWOBODYOVERLAP_H5z_F3z_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_F3z_vrr+3*oned2z*I_TWOBODYOVERLAP_H5z_D2z_vrr;

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
     * totally 30 integrals are omitted 
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
     * shell quartet name: SQ_KINETIC_G_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 45 integrals are omitted 
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
     * shell quartet name: SQ_KINETIC_G_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 30 integrals are omitted 
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
     * shell quartet name: SQ_KINETIC_H_D
     * expanding position: BRA2
     * code section is: VRR
     * totally 63 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_P
     * RHS shell quartet name: SQ_KINETIC_G_P
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_KINETIC_H5x_D2x_vrr = PBX*I_KINETIC_H5x_Px_vrr+5*oned2z*I_KINETIC_G4x_Px_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2x_vrr-adz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_H4xy_D2x_vrr = PBX*I_KINETIC_H4xy_Px_vrr+4*oned2z*I_KINETIC_G3xy_Px_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2x_vrr-adz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_D2x_vrr = PBX*I_KINETIC_H4xz_Px_vrr+4*oned2z*I_KINETIC_G3xz_Px_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2x_vrr-adz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_D2x_vrr = PBX*I_KINETIC_H3x2y_Px_vrr+3*oned2z*I_KINETIC_G2x2y_Px_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2x_vrr-adz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_H3xyz_D2x_vrr = PBX*I_KINETIC_H3xyz_Px_vrr+3*oned2z*I_KINETIC_G2xyz_Px_vrr+oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_D2x_vrr-adz*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_KINETIC_H3x2z_D2x_vrr = PBX*I_KINETIC_H3x2z_Px_vrr+3*oned2z*I_KINETIC_G2x2z_Px_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2x_vrr-adz*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_KINETIC_H2x3y_D2x_vrr = PBX*I_KINETIC_H2x3y_Px_vrr+2*oned2z*I_KINETIC_Gx3y_Px_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2x_vrr-adz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_H2x2yz_D2x_vrr = PBX*I_KINETIC_H2x2yz_Px_vrr+2*oned2z*I_KINETIC_Gx2yz_Px_vrr+oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2xy2z_D2x_vrr = PBX*I_KINETIC_H2xy2z_Px_vrr+2*oned2z*I_KINETIC_Gxy2z_Px_vrr+oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr-adz*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_KINETIC_H2x3z_D2x_vrr = PBX*I_KINETIC_H2x3z_Px_vrr+2*oned2z*I_KINETIC_Gx3z_Px_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2x_vrr-adz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_Hx4y_D2x_vrr = PBX*I_KINETIC_Hx4y_Px_vrr+oned2z*I_KINETIC_G4y_Px_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx3yz_D2x_vrr = PBX*I_KINETIC_Hx3yz_Px_vrr+oned2z*I_KINETIC_G3yz_Px_vrr+oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_KINETIC_Hx2y2z_D2x_vrr = PBX*I_KINETIC_Hx2y2z_Px_vrr+oned2z*I_KINETIC_G2y2z_Px_vrr+oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_KINETIC_Hxy3z_D2x_vrr = PBX*I_KINETIC_Hxy3z_Px_vrr+oned2z*I_KINETIC_Gy3z_Px_vrr+oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_KINETIC_Hx4z_D2x_vrr = PBX*I_KINETIC_Hx4z_Px_vrr+oned2z*I_KINETIC_G4z_Px_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_D2x_vrr = PBX*I_KINETIC_H5y_Px_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2x_vrr-adz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_H4yz_D2x_vrr = PBX*I_KINETIC_H4yz_Px_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2x_vrr-adz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_D2x_vrr = PBX*I_KINETIC_H3y2z_Px_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2x_vrr-adz*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_KINETIC_H2y3z_D2x_vrr = PBX*I_KINETIC_H2y3z_Px_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2x_vrr-adz*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_KINETIC_Hy4z_D2x_vrr = PBX*I_KINETIC_Hy4z_Px_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2x_vrr-adz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_D2x_vrr = PBX*I_KINETIC_H5z_Px_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2x_vrr-adz*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_KINETIC_H5x_D2y_vrr = PBY*I_KINETIC_H5x_Py_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2y_vrr-adz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_H4xy_D2y_vrr = PBY*I_KINETIC_H4xy_Py_vrr+oned2z*I_KINETIC_G4x_Py_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2y_vrr-adz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_D2y_vrr = PBY*I_KINETIC_H4xz_Py_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2y_vrr-adz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_D2y_vrr = PBY*I_KINETIC_H3x2y_Py_vrr+2*oned2z*I_KINETIC_G3xy_Py_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2y_vrr-adz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_H3xyz_D2y_vrr = PBY*I_KINETIC_H3xyz_Py_vrr+oned2z*I_KINETIC_G3xz_Py_vrr+oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_D2y_vrr-adz*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_KINETIC_H3x2z_D2y_vrr = PBY*I_KINETIC_H3x2z_Py_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2y_vrr-adz*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_KINETIC_H2x3y_D2y_vrr = PBY*I_KINETIC_H2x3y_Py_vrr+3*oned2z*I_KINETIC_G2x2y_Py_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2y_vrr-adz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_H2x2yz_D2y_vrr = PBY*I_KINETIC_H2x2yz_Py_vrr+2*oned2z*I_KINETIC_G2xyz_Py_vrr+oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2xy2z_D2y_vrr = PBY*I_KINETIC_H2xy2z_Py_vrr+oned2z*I_KINETIC_G2x2z_Py_vrr+oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr-adz*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_KINETIC_H2x3z_D2y_vrr = PBY*I_KINETIC_H2x3z_Py_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2y_vrr-adz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_Hx4y_D2y_vrr = PBY*I_KINETIC_Hx4y_Py_vrr+4*oned2z*I_KINETIC_Gx3y_Py_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx3yz_D2y_vrr = PBY*I_KINETIC_Hx3yz_Py_vrr+3*oned2z*I_KINETIC_Gx2yz_Py_vrr+oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_KINETIC_Hx2y2z_D2y_vrr = PBY*I_KINETIC_Hx2y2z_Py_vrr+2*oned2z*I_KINETIC_Gxy2z_Py_vrr+oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_KINETIC_Hxy3z_D2y_vrr = PBY*I_KINETIC_Hxy3z_Py_vrr+oned2z*I_KINETIC_Gx3z_Py_vrr+oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_KINETIC_Hx4z_D2y_vrr = PBY*I_KINETIC_Hx4z_Py_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_D2y_vrr = PBY*I_KINETIC_H5y_Py_vrr+5*oned2z*I_KINETIC_G4y_Py_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2y_vrr-adz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_H4yz_D2y_vrr = PBY*I_KINETIC_H4yz_Py_vrr+4*oned2z*I_KINETIC_G3yz_Py_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2y_vrr-adz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_D2y_vrr = PBY*I_KINETIC_H3y2z_Py_vrr+3*oned2z*I_KINETIC_G2y2z_Py_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2y_vrr-adz*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_KINETIC_H2y3z_D2y_vrr = PBY*I_KINETIC_H2y3z_Py_vrr+2*oned2z*I_KINETIC_Gy3z_Py_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2y_vrr-adz*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_KINETIC_Hy4z_D2y_vrr = PBY*I_KINETIC_Hy4z_Py_vrr+oned2z*I_KINETIC_G4z_Py_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2y_vrr-adz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_D2y_vrr = PBY*I_KINETIC_H5z_Py_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2y_vrr-adz*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_KINETIC_H5x_D2z_vrr = PBZ*I_KINETIC_H5x_Pz_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_D2z_vrr-adz*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_KINETIC_H4xy_D2z_vrr = PBZ*I_KINETIC_H4xy_Pz_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_D2z_vrr-adz*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_D2z_vrr = PBZ*I_KINETIC_H4xz_Pz_vrr+oned2z*I_KINETIC_G4x_Pz_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_D2z_vrr-adz*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_D2z_vrr = PBZ*I_KINETIC_H3x2y_Pz_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_D2z_vrr-adz*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_KINETIC_H3xyz_D2z_vrr = PBZ*I_KINETIC_H3xyz_Pz_vrr+oned2z*I_KINETIC_G3xy_Pz_vrr+oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_D2z_vrr-adz*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_KINETIC_H3x2z_D2z_vrr = PBZ*I_KINETIC_H3x2z_Pz_vrr+2*oned2z*I_KINETIC_G3xz_Pz_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_D2z_vrr-adz*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_KINETIC_H2x3y_D2z_vrr = PBZ*I_KINETIC_H2x3y_Pz_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_D2z_vrr-adz*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_KINETIC_H2x2yz_D2z_vrr = PBZ*I_KINETIC_H2x2yz_Pz_vrr+oned2z*I_KINETIC_G2x2y_Pz_vrr+oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2xy2z_D2z_vrr = PBZ*I_KINETIC_H2xy2z_Pz_vrr+2*oned2z*I_KINETIC_G2xyz_Pz_vrr+oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_D2z_vrr-adz*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_KINETIC_H2x3z_D2z_vrr = PBZ*I_KINETIC_H2x3z_Pz_vrr+3*oned2z*I_KINETIC_G2x2z_Pz_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_D2z_vrr-adz*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_KINETIC_Hx4y_D2z_vrr = PBZ*I_KINETIC_Hx4y_Pz_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx3yz_D2z_vrr = PBZ*I_KINETIC_Hx3yz_Pz_vrr+oned2z*I_KINETIC_Gx3y_Pz_vrr+oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_KINETIC_Hx2y2z_D2z_vrr = PBZ*I_KINETIC_Hx2y2z_Pz_vrr+2*oned2z*I_KINETIC_Gx2yz_Pz_vrr+oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_KINETIC_Hxy3z_D2z_vrr = PBZ*I_KINETIC_Hxy3z_Pz_vrr+3*oned2z*I_KINETIC_Gxy2z_Pz_vrr+oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_KINETIC_Hx4z_D2z_vrr = PBZ*I_KINETIC_Hx4z_Pz_vrr+4*oned2z*I_KINETIC_Gx3z_Pz_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_D2z_vrr = PBZ*I_KINETIC_H5y_Pz_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_D2z_vrr-adz*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_KINETIC_H4yz_D2z_vrr = PBZ*I_KINETIC_H4yz_Pz_vrr+oned2z*I_KINETIC_G4y_Pz_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_D2z_vrr-adz*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_D2z_vrr = PBZ*I_KINETIC_H3y2z_Pz_vrr+2*oned2z*I_KINETIC_G3yz_Pz_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_D2z_vrr-adz*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_KINETIC_H2y3z_D2z_vrr = PBZ*I_KINETIC_H2y3z_Pz_vrr+3*oned2z*I_KINETIC_G2y2z_Pz_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_D2z_vrr-adz*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_KINETIC_Hy4z_D2z_vrr = PBZ*I_KINETIC_Hy4z_Pz_vrr+4*oned2z*I_KINETIC_Gy3z_Pz_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_D2z_vrr-adz*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_D2z_vrr = PBZ*I_KINETIC_H5z_Pz_vrr+5*oned2z*I_KINETIC_G4z_Pz_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_D2z_vrr-adz*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_F
     * expanding position: BRA2
     * code section is: VRR
     * totally 42 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_D
     * RHS shell quartet name: SQ_KINETIC_G_D
     * RHS shell quartet name: SQ_KINETIC_H_P
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_F
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     ************************************************************/
    Double I_KINETIC_H5x_F3x_vrr = PBX*I_KINETIC_H5x_D2x_vrr+5*oned2z*I_KINETIC_G4x_D2x_vrr+2*oned2z*I_KINETIC_H5x_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_KINETIC_H4xy_F3x_vrr = PBX*I_KINETIC_H4xy_D2x_vrr+4*oned2z*I_KINETIC_G3xy_D2x_vrr+2*oned2z*I_KINETIC_H4xy_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_KINETIC_H4xz_F3x_vrr = PBX*I_KINETIC_H4xz_D2x_vrr+4*oned2z*I_KINETIC_G3xz_D2x_vrr+2*oned2z*I_KINETIC_H4xz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_KINETIC_H3x2y_F3x_vrr = PBX*I_KINETIC_H3x2y_D2x_vrr+3*oned2z*I_KINETIC_G2x2y_D2x_vrr+2*oned2z*I_KINETIC_H3x2y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_KINETIC_H3xyz_F3x_vrr = PBX*I_KINETIC_H3xyz_D2x_vrr+3*oned2z*I_KINETIC_G2xyz_D2x_vrr+2*oned2z*I_KINETIC_H3xyz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H3xyz_Px_vrr;
    Double I_KINETIC_H3x2z_F3x_vrr = PBX*I_KINETIC_H3x2z_D2x_vrr+3*oned2z*I_KINETIC_G2x2z_D2x_vrr+2*oned2z*I_KINETIC_H3x2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_KINETIC_H2x3y_F3x_vrr = PBX*I_KINETIC_H2x3y_D2x_vrr+2*oned2z*I_KINETIC_Gx3y_D2x_vrr+2*oned2z*I_KINETIC_H2x3y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_KINETIC_H2x2yz_F3x_vrr = PBX*I_KINETIC_H2x2yz_D2x_vrr+2*oned2z*I_KINETIC_Gx2yz_D2x_vrr+2*oned2z*I_KINETIC_H2x2yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_KINETIC_H2xy2z_F3x_vrr = PBX*I_KINETIC_H2xy2z_D2x_vrr+2*oned2z*I_KINETIC_Gxy2z_D2x_vrr+2*oned2z*I_KINETIC_H2xy2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H2xy2z_Px_vrr;
    Double I_KINETIC_H2x3z_F3x_vrr = PBX*I_KINETIC_H2x3z_D2x_vrr+2*oned2z*I_KINETIC_Gx3z_D2x_vrr+2*oned2z*I_KINETIC_H2x3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_KINETIC_Hx4y_F3x_vrr = PBX*I_KINETIC_Hx4y_D2x_vrr+oned2z*I_KINETIC_G4y_D2x_vrr+2*oned2z*I_KINETIC_Hx4y_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_KINETIC_Hx3yz_F3x_vrr = PBX*I_KINETIC_Hx3yz_D2x_vrr+oned2z*I_KINETIC_G3yz_D2x_vrr+2*oned2z*I_KINETIC_Hx3yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Hx3yz_Px_vrr;
    Double I_KINETIC_Hx2y2z_F3x_vrr = PBX*I_KINETIC_Hx2y2z_D2x_vrr+oned2z*I_KINETIC_G2y2z_D2x_vrr+2*oned2z*I_KINETIC_Hx2y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Hx2y2z_Px_vrr;
    Double I_KINETIC_Hxy3z_F3x_vrr = PBX*I_KINETIC_Hxy3z_D2x_vrr+oned2z*I_KINETIC_Gy3z_D2x_vrr+2*oned2z*I_KINETIC_Hxy3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Hxy3z_Px_vrr;
    Double I_KINETIC_Hx4z_F3x_vrr = PBX*I_KINETIC_Hx4z_D2x_vrr+oned2z*I_KINETIC_G4z_D2x_vrr+2*oned2z*I_KINETIC_Hx4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_KINETIC_H5y_F3x_vrr = PBX*I_KINETIC_H5y_D2x_vrr+2*oned2z*I_KINETIC_H5y_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_KINETIC_H4yz_F3x_vrr = PBX*I_KINETIC_H4yz_D2x_vrr+2*oned2z*I_KINETIC_H4yz_Px_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_KINETIC_H3y2z_F3x_vrr = PBX*I_KINETIC_H3y2z_D2x_vrr+2*oned2z*I_KINETIC_H3y2z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_KINETIC_H2y3z_F3x_vrr = PBX*I_KINETIC_H2y3z_D2x_vrr+2*oned2z*I_KINETIC_H2y3z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_KINETIC_Hy4z_F3x_vrr = PBX*I_KINETIC_Hy4z_D2x_vrr+2*oned2z*I_KINETIC_Hy4z_Px_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_KINETIC_H5z_F3x_vrr = PBX*I_KINETIC_H5z_D2x_vrr+2*oned2z*I_KINETIC_H5z_Px_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F3x_vrr-2*adz*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_KINETIC_H5x_F2xy_vrr = PBY*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F2xy_vrr;
    Double I_KINETIC_H4xy_F2xy_vrr = PBY*I_KINETIC_H4xy_D2x_vrr+oned2z*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F2xy_vrr;
    Double I_KINETIC_H4xz_F2xy_vrr = PBY*I_KINETIC_H4xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F2xy_vrr;
    Double I_KINETIC_H3x2y_F2xy_vrr = PBY*I_KINETIC_H3x2y_D2x_vrr+2*oned2z*I_KINETIC_G3xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F2xy_vrr;
    Double I_KINETIC_H3xyz_F2xy_vrr = PBY*I_KINETIC_H3xyz_D2x_vrr+oned2z*I_KINETIC_G3xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F2xy_vrr;
    Double I_KINETIC_H3x2z_F2xy_vrr = PBY*I_KINETIC_H3x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F2xy_vrr;
    Double I_KINETIC_H2x3y_F2xy_vrr = PBY*I_KINETIC_H2x3y_D2x_vrr+3*oned2z*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F2xy_vrr;
    Double I_KINETIC_H2x2yz_F2xy_vrr = PBY*I_KINETIC_H2x2yz_D2x_vrr+2*oned2z*I_KINETIC_G2xyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F2xy_vrr;
    Double I_KINETIC_H2xy2z_F2xy_vrr = PBY*I_KINETIC_H2xy2z_D2x_vrr+oned2z*I_KINETIC_G2x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F2xy_vrr;
    Double I_KINETIC_H2x3z_F2xy_vrr = PBY*I_KINETIC_H2x3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F2xy_vrr;
    Double I_KINETIC_Hx4y_F2xy_vrr = PBY*I_KINETIC_Hx4y_D2x_vrr+4*oned2z*I_KINETIC_Gx3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F2xy_vrr;
    Double I_KINETIC_Hx3yz_F2xy_vrr = PBY*I_KINETIC_Hx3yz_D2x_vrr+3*oned2z*I_KINETIC_Gx2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F2xy_vrr;
    Double I_KINETIC_Hx2y2z_F2xy_vrr = PBY*I_KINETIC_Hx2y2z_D2x_vrr+2*oned2z*I_KINETIC_Gxy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F2xy_vrr;
    Double I_KINETIC_Hxy3z_F2xy_vrr = PBY*I_KINETIC_Hxy3z_D2x_vrr+oned2z*I_KINETIC_Gx3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F2xy_vrr;
    Double I_KINETIC_Hx4z_F2xy_vrr = PBY*I_KINETIC_Hx4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F2xy_vrr;
    Double I_KINETIC_H5y_F2xy_vrr = PBY*I_KINETIC_H5y_D2x_vrr+5*oned2z*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F2xy_vrr;
    Double I_KINETIC_H4yz_F2xy_vrr = PBY*I_KINETIC_H4yz_D2x_vrr+4*oned2z*I_KINETIC_G3yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F2xy_vrr;
    Double I_KINETIC_H3y2z_F2xy_vrr = PBY*I_KINETIC_H3y2z_D2x_vrr+3*oned2z*I_KINETIC_G2y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F2xy_vrr;
    Double I_KINETIC_H2y3z_F2xy_vrr = PBY*I_KINETIC_H2y3z_D2x_vrr+2*oned2z*I_KINETIC_Gy3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F2xy_vrr;
    Double I_KINETIC_Hy4z_F2xy_vrr = PBY*I_KINETIC_Hy4z_D2x_vrr+oned2z*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F2xy_vrr;
    Double I_KINETIC_H5z_F2xy_vrr = PBY*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F2xy_vrr;
    Double I_KINETIC_H5x_F2xz_vrr = PBZ*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F2xz_vrr;
    Double I_KINETIC_H4xy_F2xz_vrr = PBZ*I_KINETIC_H4xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F2xz_vrr;
    Double I_KINETIC_H4xz_F2xz_vrr = PBZ*I_KINETIC_H4xz_D2x_vrr+oned2z*I_KINETIC_G4x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F2xz_vrr;
    Double I_KINETIC_H3x2y_F2xz_vrr = PBZ*I_KINETIC_H3x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F2xz_vrr;
    Double I_KINETIC_H3xyz_F2xz_vrr = PBZ*I_KINETIC_H3xyz_D2x_vrr+oned2z*I_KINETIC_G3xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F2xz_vrr;
    Double I_KINETIC_H3x2z_F2xz_vrr = PBZ*I_KINETIC_H3x2z_D2x_vrr+2*oned2z*I_KINETIC_G3xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F2xz_vrr;
    Double I_KINETIC_H2x3y_F2xz_vrr = PBZ*I_KINETIC_H2x3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F2xz_vrr;
    Double I_KINETIC_H2x2yz_F2xz_vrr = PBZ*I_KINETIC_H2x2yz_D2x_vrr+oned2z*I_KINETIC_G2x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F2xz_vrr;
    Double I_KINETIC_H2xy2z_F2xz_vrr = PBZ*I_KINETIC_H2xy2z_D2x_vrr+2*oned2z*I_KINETIC_G2xyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F2xz_vrr;
    Double I_KINETIC_H2x3z_F2xz_vrr = PBZ*I_KINETIC_H2x3z_D2x_vrr+3*oned2z*I_KINETIC_G2x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F2xz_vrr;
    Double I_KINETIC_Hx4y_F2xz_vrr = PBZ*I_KINETIC_Hx4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F2xz_vrr;
    Double I_KINETIC_Hx3yz_F2xz_vrr = PBZ*I_KINETIC_Hx3yz_D2x_vrr+oned2z*I_KINETIC_Gx3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F2xz_vrr;
    Double I_KINETIC_Hx2y2z_F2xz_vrr = PBZ*I_KINETIC_Hx2y2z_D2x_vrr+2*oned2z*I_KINETIC_Gx2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F2xz_vrr;
    Double I_KINETIC_Hxy3z_F2xz_vrr = PBZ*I_KINETIC_Hxy3z_D2x_vrr+3*oned2z*I_KINETIC_Gxy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F2xz_vrr;
    Double I_KINETIC_Hx4z_F2xz_vrr = PBZ*I_KINETIC_Hx4z_D2x_vrr+4*oned2z*I_KINETIC_Gx3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F2xz_vrr;
    Double I_KINETIC_H5y_F2xz_vrr = PBZ*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F2xz_vrr;
    Double I_KINETIC_H4yz_F2xz_vrr = PBZ*I_KINETIC_H4yz_D2x_vrr+oned2z*I_KINETIC_G4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F2xz_vrr;
    Double I_KINETIC_H3y2z_F2xz_vrr = PBZ*I_KINETIC_H3y2z_D2x_vrr+2*oned2z*I_KINETIC_G3yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F2xz_vrr;
    Double I_KINETIC_H2y3z_F2xz_vrr = PBZ*I_KINETIC_H2y3z_D2x_vrr+3*oned2z*I_KINETIC_G2y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F2xz_vrr;
    Double I_KINETIC_Hy4z_F2xz_vrr = PBZ*I_KINETIC_Hy4z_D2x_vrr+4*oned2z*I_KINETIC_Gy3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F2xz_vrr;
    Double I_KINETIC_H5z_F2xz_vrr = PBZ*I_KINETIC_H5z_D2x_vrr+5*oned2z*I_KINETIC_G4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F2xz_vrr;
    Double I_KINETIC_H5x_Fx2y_vrr = PBX*I_KINETIC_H5x_D2y_vrr+5*oned2z*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fx2y_vrr;
    Double I_KINETIC_H4xy_Fx2y_vrr = PBX*I_KINETIC_H4xy_D2y_vrr+4*oned2z*I_KINETIC_G3xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fx2y_vrr;
    Double I_KINETIC_H4xz_Fx2y_vrr = PBX*I_KINETIC_H4xz_D2y_vrr+4*oned2z*I_KINETIC_G3xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fx2y_vrr;
    Double I_KINETIC_H3x2y_Fx2y_vrr = PBX*I_KINETIC_H3x2y_D2y_vrr+3*oned2z*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fx2y_vrr;
    Double I_KINETIC_H3xyz_Fx2y_vrr = PBX*I_KINETIC_H3xyz_D2y_vrr+3*oned2z*I_KINETIC_G2xyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Fx2y_vrr;
    Double I_KINETIC_H3x2z_Fx2y_vrr = PBX*I_KINETIC_H3x2z_D2y_vrr+3*oned2z*I_KINETIC_G2x2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fx2y_vrr;
    Double I_KINETIC_H2x3y_Fx2y_vrr = PBX*I_KINETIC_H2x3y_D2y_vrr+2*oned2z*I_KINETIC_Gx3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fx2y_vrr;
    Double I_KINETIC_H2x2yz_Fx2y_vrr = PBX*I_KINETIC_H2x2yz_D2y_vrr+2*oned2z*I_KINETIC_Gx2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fx2y_vrr;
    Double I_KINETIC_H2xy2z_Fx2y_vrr = PBX*I_KINETIC_H2xy2z_D2y_vrr+2*oned2z*I_KINETIC_Gxy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Fx2y_vrr;
    Double I_KINETIC_H2x3z_Fx2y_vrr = PBX*I_KINETIC_H2x3z_D2y_vrr+2*oned2z*I_KINETIC_Gx3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fx2y_vrr;
    Double I_KINETIC_Hx4y_Fx2y_vrr = PBX*I_KINETIC_Hx4y_D2y_vrr+oned2z*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fx2y_vrr;
    Double I_KINETIC_Hx3yz_Fx2y_vrr = PBX*I_KINETIC_Hx3yz_D2y_vrr+oned2z*I_KINETIC_G3yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Fx2y_vrr;
    Double I_KINETIC_Hx2y2z_Fx2y_vrr = PBX*I_KINETIC_Hx2y2z_D2y_vrr+oned2z*I_KINETIC_G2y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Fx2y_vrr;
    Double I_KINETIC_Hxy3z_Fx2y_vrr = PBX*I_KINETIC_Hxy3z_D2y_vrr+oned2z*I_KINETIC_Gy3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Fx2y_vrr;
    Double I_KINETIC_Hx4z_Fx2y_vrr = PBX*I_KINETIC_Hx4z_D2y_vrr+oned2z*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fx2y_vrr;
    Double I_KINETIC_H5y_Fx2y_vrr = PBX*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fx2y_vrr;
    Double I_KINETIC_H4yz_Fx2y_vrr = PBX*I_KINETIC_H4yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fx2y_vrr;
    Double I_KINETIC_H3y2z_Fx2y_vrr = PBX*I_KINETIC_H3y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fx2y_vrr;
    Double I_KINETIC_H2y3z_Fx2y_vrr = PBX*I_KINETIC_H2y3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fx2y_vrr;
    Double I_KINETIC_Hy4z_Fx2y_vrr = PBX*I_KINETIC_Hy4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fx2y_vrr;
    Double I_KINETIC_H5z_Fx2y_vrr = PBX*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fx2y_vrr;
    Double I_KINETIC_H5x_Fx2z_vrr = PBX*I_KINETIC_H5x_D2z_vrr+5*oned2z*I_KINETIC_G4x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Fx2z_vrr;
    Double I_KINETIC_H4xy_Fx2z_vrr = PBX*I_KINETIC_H4xy_D2z_vrr+4*oned2z*I_KINETIC_G3xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Fx2z_vrr;
    Double I_KINETIC_H4xz_Fx2z_vrr = PBX*I_KINETIC_H4xz_D2z_vrr+4*oned2z*I_KINETIC_G3xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Fx2z_vrr;
    Double I_KINETIC_H3x2y_Fx2z_vrr = PBX*I_KINETIC_H3x2y_D2z_vrr+3*oned2z*I_KINETIC_G2x2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Fx2z_vrr;
    Double I_KINETIC_H3xyz_Fx2z_vrr = PBX*I_KINETIC_H3xyz_D2z_vrr+3*oned2z*I_KINETIC_G2xyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Fx2z_vrr;
    Double I_KINETIC_H3x2z_Fx2z_vrr = PBX*I_KINETIC_H3x2z_D2z_vrr+3*oned2z*I_KINETIC_G2x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Fx2z_vrr;
    Double I_KINETIC_H2x3y_Fx2z_vrr = PBX*I_KINETIC_H2x3y_D2z_vrr+2*oned2z*I_KINETIC_Gx3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Fx2z_vrr;
    Double I_KINETIC_H2x2yz_Fx2z_vrr = PBX*I_KINETIC_H2x2yz_D2z_vrr+2*oned2z*I_KINETIC_Gx2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Fx2z_vrr;
    Double I_KINETIC_H2xy2z_Fx2z_vrr = PBX*I_KINETIC_H2xy2z_D2z_vrr+2*oned2z*I_KINETIC_Gxy2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Fx2z_vrr;
    Double I_KINETIC_H2x3z_Fx2z_vrr = PBX*I_KINETIC_H2x3z_D2z_vrr+2*oned2z*I_KINETIC_Gx3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Fx2z_vrr;
    Double I_KINETIC_Hx4y_Fx2z_vrr = PBX*I_KINETIC_Hx4y_D2z_vrr+oned2z*I_KINETIC_G4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Fx2z_vrr;
    Double I_KINETIC_Hx3yz_Fx2z_vrr = PBX*I_KINETIC_Hx3yz_D2z_vrr+oned2z*I_KINETIC_G3yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Fx2z_vrr;
    Double I_KINETIC_Hx2y2z_Fx2z_vrr = PBX*I_KINETIC_Hx2y2z_D2z_vrr+oned2z*I_KINETIC_G2y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Fx2z_vrr;
    Double I_KINETIC_Hxy3z_Fx2z_vrr = PBX*I_KINETIC_Hxy3z_D2z_vrr+oned2z*I_KINETIC_Gy3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Fx2z_vrr;
    Double I_KINETIC_Hx4z_Fx2z_vrr = PBX*I_KINETIC_Hx4z_D2z_vrr+oned2z*I_KINETIC_G4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Fx2z_vrr;
    Double I_KINETIC_H5y_Fx2z_vrr = PBX*I_KINETIC_H5y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Fx2z_vrr;
    Double I_KINETIC_H4yz_Fx2z_vrr = PBX*I_KINETIC_H4yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Fx2z_vrr;
    Double I_KINETIC_H3y2z_Fx2z_vrr = PBX*I_KINETIC_H3y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Fx2z_vrr;
    Double I_KINETIC_H2y3z_Fx2z_vrr = PBX*I_KINETIC_H2y3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Fx2z_vrr;
    Double I_KINETIC_Hy4z_Fx2z_vrr = PBX*I_KINETIC_Hy4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Fx2z_vrr;
    Double I_KINETIC_H5z_Fx2z_vrr = PBX*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Fx2z_vrr;
    Double I_KINETIC_H5x_F3y_vrr = PBY*I_KINETIC_H5x_D2y_vrr+2*oned2z*I_KINETIC_H5x_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_KINETIC_H4xy_F3y_vrr = PBY*I_KINETIC_H4xy_D2y_vrr+oned2z*I_KINETIC_G4x_D2y_vrr+2*oned2z*I_KINETIC_H4xy_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_KINETIC_H4xz_F3y_vrr = PBY*I_KINETIC_H4xz_D2y_vrr+2*oned2z*I_KINETIC_H4xz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_KINETIC_H3x2y_F3y_vrr = PBY*I_KINETIC_H3x2y_D2y_vrr+2*oned2z*I_KINETIC_G3xy_D2y_vrr+2*oned2z*I_KINETIC_H3x2y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_KINETIC_H3xyz_F3y_vrr = PBY*I_KINETIC_H3xyz_D2y_vrr+oned2z*I_KINETIC_G3xz_D2y_vrr+2*oned2z*I_KINETIC_H3xyz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H3xyz_Py_vrr;
    Double I_KINETIC_H3x2z_F3y_vrr = PBY*I_KINETIC_H3x2z_D2y_vrr+2*oned2z*I_KINETIC_H3x2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_KINETIC_H2x3y_F3y_vrr = PBY*I_KINETIC_H2x3y_D2y_vrr+3*oned2z*I_KINETIC_G2x2y_D2y_vrr+2*oned2z*I_KINETIC_H2x3y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_KINETIC_H2x2yz_F3y_vrr = PBY*I_KINETIC_H2x2yz_D2y_vrr+2*oned2z*I_KINETIC_G2xyz_D2y_vrr+2*oned2z*I_KINETIC_H2x2yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H2x2yz_Py_vrr;
    Double I_KINETIC_H2xy2z_F3y_vrr = PBY*I_KINETIC_H2xy2z_D2y_vrr+oned2z*I_KINETIC_G2x2z_D2y_vrr+2*oned2z*I_KINETIC_H2xy2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H2xy2z_Py_vrr;
    Double I_KINETIC_H2x3z_F3y_vrr = PBY*I_KINETIC_H2x3z_D2y_vrr+2*oned2z*I_KINETIC_H2x3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H2x3z_Py_vrr;
    Double I_KINETIC_Hx4y_F3y_vrr = PBY*I_KINETIC_Hx4y_D2y_vrr+4*oned2z*I_KINETIC_Gx3y_D2y_vrr+2*oned2z*I_KINETIC_Hx4y_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_KINETIC_Hx3yz_F3y_vrr = PBY*I_KINETIC_Hx3yz_D2y_vrr+3*oned2z*I_KINETIC_Gx2yz_D2y_vrr+2*oned2z*I_KINETIC_Hx3yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Hx3yz_Py_vrr;
    Double I_KINETIC_Hx2y2z_F3y_vrr = PBY*I_KINETIC_Hx2y2z_D2y_vrr+2*oned2z*I_KINETIC_Gxy2z_D2y_vrr+2*oned2z*I_KINETIC_Hx2y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Hx2y2z_Py_vrr;
    Double I_KINETIC_Hxy3z_F3y_vrr = PBY*I_KINETIC_Hxy3z_D2y_vrr+oned2z*I_KINETIC_Gx3z_D2y_vrr+2*oned2z*I_KINETIC_Hxy3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Hxy3z_Py_vrr;
    Double I_KINETIC_Hx4z_F3y_vrr = PBY*I_KINETIC_Hx4z_D2y_vrr+2*oned2z*I_KINETIC_Hx4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_KINETIC_H5y_F3y_vrr = PBY*I_KINETIC_H5y_D2y_vrr+5*oned2z*I_KINETIC_G4y_D2y_vrr+2*oned2z*I_KINETIC_H5y_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_KINETIC_H4yz_F3y_vrr = PBY*I_KINETIC_H4yz_D2y_vrr+4*oned2z*I_KINETIC_G3yz_D2y_vrr+2*oned2z*I_KINETIC_H4yz_Py_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_KINETIC_H3y2z_F3y_vrr = PBY*I_KINETIC_H3y2z_D2y_vrr+3*oned2z*I_KINETIC_G2y2z_D2y_vrr+2*oned2z*I_KINETIC_H3y2z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_KINETIC_H2y3z_F3y_vrr = PBY*I_KINETIC_H2y3z_D2y_vrr+2*oned2z*I_KINETIC_Gy3z_D2y_vrr+2*oned2z*I_KINETIC_H2y3z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H2y3z_Py_vrr;
    Double I_KINETIC_Hy4z_F3y_vrr = PBY*I_KINETIC_Hy4z_D2y_vrr+oned2z*I_KINETIC_G4z_D2y_vrr+2*oned2z*I_KINETIC_Hy4z_Py_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_KINETIC_H5z_F3y_vrr = PBY*I_KINETIC_H5z_D2y_vrr+2*oned2z*I_KINETIC_H5z_Py_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F3y_vrr-2*adz*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_KINETIC_H5x_F2yz_vrr = PBZ*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F2yz_vrr;
    Double I_KINETIC_H4xy_F2yz_vrr = PBZ*I_KINETIC_H4xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F2yz_vrr;
    Double I_KINETIC_H4xz_F2yz_vrr = PBZ*I_KINETIC_H4xz_D2y_vrr+oned2z*I_KINETIC_G4x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F2yz_vrr;
    Double I_KINETIC_H3x2y_F2yz_vrr = PBZ*I_KINETIC_H3x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F2yz_vrr;
    Double I_KINETIC_H3xyz_F2yz_vrr = PBZ*I_KINETIC_H3xyz_D2y_vrr+oned2z*I_KINETIC_G3xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F2yz_vrr;
    Double I_KINETIC_H3x2z_F2yz_vrr = PBZ*I_KINETIC_H3x2z_D2y_vrr+2*oned2z*I_KINETIC_G3xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F2yz_vrr;
    Double I_KINETIC_H2x3y_F2yz_vrr = PBZ*I_KINETIC_H2x3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F2yz_vrr;
    Double I_KINETIC_H2x2yz_F2yz_vrr = PBZ*I_KINETIC_H2x2yz_D2y_vrr+oned2z*I_KINETIC_G2x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F2yz_vrr;
    Double I_KINETIC_H2xy2z_F2yz_vrr = PBZ*I_KINETIC_H2xy2z_D2y_vrr+2*oned2z*I_KINETIC_G2xyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F2yz_vrr;
    Double I_KINETIC_H2x3z_F2yz_vrr = PBZ*I_KINETIC_H2x3z_D2y_vrr+3*oned2z*I_KINETIC_G2x2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F2yz_vrr;
    Double I_KINETIC_Hx4y_F2yz_vrr = PBZ*I_KINETIC_Hx4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F2yz_vrr;
    Double I_KINETIC_Hx3yz_F2yz_vrr = PBZ*I_KINETIC_Hx3yz_D2y_vrr+oned2z*I_KINETIC_Gx3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F2yz_vrr;
    Double I_KINETIC_Hx2y2z_F2yz_vrr = PBZ*I_KINETIC_Hx2y2z_D2y_vrr+2*oned2z*I_KINETIC_Gx2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F2yz_vrr;
    Double I_KINETIC_Hxy3z_F2yz_vrr = PBZ*I_KINETIC_Hxy3z_D2y_vrr+3*oned2z*I_KINETIC_Gxy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F2yz_vrr;
    Double I_KINETIC_Hx4z_F2yz_vrr = PBZ*I_KINETIC_Hx4z_D2y_vrr+4*oned2z*I_KINETIC_Gx3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F2yz_vrr;
    Double I_KINETIC_H5y_F2yz_vrr = PBZ*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F2yz_vrr;
    Double I_KINETIC_H4yz_F2yz_vrr = PBZ*I_KINETIC_H4yz_D2y_vrr+oned2z*I_KINETIC_G4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F2yz_vrr;
    Double I_KINETIC_H3y2z_F2yz_vrr = PBZ*I_KINETIC_H3y2z_D2y_vrr+2*oned2z*I_KINETIC_G3yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F2yz_vrr;
    Double I_KINETIC_H2y3z_F2yz_vrr = PBZ*I_KINETIC_H2y3z_D2y_vrr+3*oned2z*I_KINETIC_G2y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F2yz_vrr;
    Double I_KINETIC_Hy4z_F2yz_vrr = PBZ*I_KINETIC_Hy4z_D2y_vrr+4*oned2z*I_KINETIC_Gy3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F2yz_vrr;
    Double I_KINETIC_H5z_F2yz_vrr = PBZ*I_KINETIC_H5z_D2y_vrr+5*oned2z*I_KINETIC_G4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F2yz_vrr;
    Double I_KINETIC_H5x_F3z_vrr = PBZ*I_KINETIC_H5x_D2z_vrr+2*oned2z*I_KINETIC_H5x_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5x_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_KINETIC_H4xy_F3z_vrr = PBZ*I_KINETIC_H4xy_D2z_vrr+2*oned2z*I_KINETIC_H4xy_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_KINETIC_H4xz_F3z_vrr = PBZ*I_KINETIC_H4xz_D2z_vrr+oned2z*I_KINETIC_G4x_D2z_vrr+2*oned2z*I_KINETIC_H4xz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H4xz_Pz_vrr;
    Double I_KINETIC_H3x2y_F3z_vrr = PBZ*I_KINETIC_H3x2y_D2z_vrr+2*oned2z*I_KINETIC_H3x2y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H3x2y_Pz_vrr;
    Double I_KINETIC_H3xyz_F3z_vrr = PBZ*I_KINETIC_H3xyz_D2z_vrr+oned2z*I_KINETIC_G3xy_D2z_vrr+2*oned2z*I_KINETIC_H3xyz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H3xyz_Pz_vrr;
    Double I_KINETIC_H3x2z_F3z_vrr = PBZ*I_KINETIC_H3x2z_D2z_vrr+2*oned2z*I_KINETIC_G3xz_D2z_vrr+2*oned2z*I_KINETIC_H3x2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H3x2z_Pz_vrr;
    Double I_KINETIC_H2x3y_F3z_vrr = PBZ*I_KINETIC_H2x3y_D2z_vrr+2*oned2z*I_KINETIC_H2x3y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H2x3y_Pz_vrr;
    Double I_KINETIC_H2x2yz_F3z_vrr = PBZ*I_KINETIC_H2x2yz_D2z_vrr+oned2z*I_KINETIC_G2x2y_D2z_vrr+2*oned2z*I_KINETIC_H2x2yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr;
    Double I_KINETIC_H2xy2z_F3z_vrr = PBZ*I_KINETIC_H2xy2z_D2z_vrr+2*oned2z*I_KINETIC_G2xyz_D2z_vrr+2*oned2z*I_KINETIC_H2xy2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H2xy2z_Pz_vrr;
    Double I_KINETIC_H2x3z_F3z_vrr = PBZ*I_KINETIC_H2x3z_D2z_vrr+3*oned2z*I_KINETIC_G2x2z_D2z_vrr+2*oned2z*I_KINETIC_H2x3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H2x3z_Pz_vrr;
    Double I_KINETIC_Hx4y_F3z_vrr = PBZ*I_KINETIC_Hx4y_D2z_vrr+2*oned2z*I_KINETIC_Hx4y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_KINETIC_Hx3yz_F3z_vrr = PBZ*I_KINETIC_Hx3yz_D2z_vrr+oned2z*I_KINETIC_Gx3y_D2z_vrr+2*oned2z*I_KINETIC_Hx3yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Hx3yz_Pz_vrr;
    Double I_KINETIC_Hx2y2z_F3z_vrr = PBZ*I_KINETIC_Hx2y2z_D2z_vrr+2*oned2z*I_KINETIC_Gx2yz_D2z_vrr+2*oned2z*I_KINETIC_Hx2y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr;
    Double I_KINETIC_Hxy3z_F3z_vrr = PBZ*I_KINETIC_Hxy3z_D2z_vrr+3*oned2z*I_KINETIC_Gxy2z_D2z_vrr+2*oned2z*I_KINETIC_Hxy3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Hxy3z_Pz_vrr;
    Double I_KINETIC_Hx4z_F3z_vrr = PBZ*I_KINETIC_Hx4z_D2z_vrr+4*oned2z*I_KINETIC_Gx3z_D2z_vrr+2*oned2z*I_KINETIC_Hx4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_KINETIC_H5y_F3z_vrr = PBZ*I_KINETIC_H5y_D2z_vrr+2*oned2z*I_KINETIC_H5y_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5y_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_KINETIC_H4yz_F3z_vrr = PBZ*I_KINETIC_H4yz_D2z_vrr+oned2z*I_KINETIC_G4y_D2z_vrr+2*oned2z*I_KINETIC_H4yz_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H4yz_Pz_vrr;
    Double I_KINETIC_H3y2z_F3z_vrr = PBZ*I_KINETIC_H3y2z_D2z_vrr+2*oned2z*I_KINETIC_G3yz_D2z_vrr+2*oned2z*I_KINETIC_H3y2z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H3y2z_Pz_vrr;
    Double I_KINETIC_H2y3z_F3z_vrr = PBZ*I_KINETIC_H2y3z_D2z_vrr+3*oned2z*I_KINETIC_G2y2z_D2z_vrr+2*oned2z*I_KINETIC_H2y3z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H2y3z_Pz_vrr;
    Double I_KINETIC_Hy4z_F3z_vrr = PBZ*I_KINETIC_Hy4z_D2z_vrr+4*oned2z*I_KINETIC_Gy3z_D2z_vrr+2*oned2z*I_KINETIC_Hy4z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_Hy4z_Pz_vrr;
    Double I_KINETIC_H5z_F3z_vrr = PBZ*I_KINETIC_H5z_D2z_vrr+5*oned2z*I_KINETIC_G4z_D2z_vrr+2*oned2z*I_KINETIC_H5z_Pz_vrr+twoxi*I_TWOBODYOVERLAP_H5z_F3z_vrr-2*adz*I_TWOBODYOVERLAP_H5z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_G
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_F
     * RHS shell quartet name: SQ_KINETIC_G_F
     * RHS shell quartet name: SQ_KINETIC_H_D
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_G
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
     ************************************************************/
    Double I_KINETIC_H5x_G4x_vrr = PBX*I_KINETIC_H5x_F3x_vrr+5*oned2z*I_KINETIC_G4x_F3x_vrr+3*oned2z*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_KINETIC_H4xy_G4x_vrr = PBX*I_KINETIC_H4xy_F3x_vrr+4*oned2z*I_KINETIC_G3xy_F3x_vrr+3*oned2z*I_KINETIC_H4xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_KINETIC_H4xz_G4x_vrr = PBX*I_KINETIC_H4xz_F3x_vrr+4*oned2z*I_KINETIC_G3xz_F3x_vrr+3*oned2z*I_KINETIC_H4xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_KINETIC_H3x2y_G4x_vrr = PBX*I_KINETIC_H3x2y_F3x_vrr+3*oned2z*I_KINETIC_G2x2y_F3x_vrr+3*oned2z*I_KINETIC_H3x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_KINETIC_H3xyz_G4x_vrr = PBX*I_KINETIC_H3xyz_F3x_vrr+3*oned2z*I_KINETIC_G2xyz_F3x_vrr+3*oned2z*I_KINETIC_H3xyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H3xyz_D2x_vrr;
    Double I_KINETIC_H3x2z_G4x_vrr = PBX*I_KINETIC_H3x2z_F3x_vrr+3*oned2z*I_KINETIC_G2x2z_F3x_vrr+3*oned2z*I_KINETIC_H3x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_KINETIC_H2x3y_G4x_vrr = PBX*I_KINETIC_H2x3y_F3x_vrr+2*oned2z*I_KINETIC_Gx3y_F3x_vrr+3*oned2z*I_KINETIC_H2x3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_KINETIC_H2x2yz_G4x_vrr = PBX*I_KINETIC_H2x2yz_F3x_vrr+2*oned2z*I_KINETIC_Gx2yz_F3x_vrr+3*oned2z*I_KINETIC_H2x2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr;
    Double I_KINETIC_H2xy2z_G4x_vrr = PBX*I_KINETIC_H2xy2z_F3x_vrr+2*oned2z*I_KINETIC_Gxy2z_F3x_vrr+3*oned2z*I_KINETIC_H2xy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr;
    Double I_KINETIC_H2x3z_G4x_vrr = PBX*I_KINETIC_H2x3z_F3x_vrr+2*oned2z*I_KINETIC_Gx3z_F3x_vrr+3*oned2z*I_KINETIC_H2x3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_KINETIC_Hx4y_G4x_vrr = PBX*I_KINETIC_Hx4y_F3x_vrr+oned2z*I_KINETIC_G4y_F3x_vrr+3*oned2z*I_KINETIC_Hx4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_KINETIC_Hx3yz_G4x_vrr = PBX*I_KINETIC_Hx3yz_F3x_vrr+oned2z*I_KINETIC_G3yz_F3x_vrr+3*oned2z*I_KINETIC_Hx3yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr;
    Double I_KINETIC_Hx2y2z_G4x_vrr = PBX*I_KINETIC_Hx2y2z_F3x_vrr+oned2z*I_KINETIC_G2y2z_F3x_vrr+3*oned2z*I_KINETIC_Hx2y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr;
    Double I_KINETIC_Hxy3z_G4x_vrr = PBX*I_KINETIC_Hxy3z_F3x_vrr+oned2z*I_KINETIC_Gy3z_F3x_vrr+3*oned2z*I_KINETIC_Hxy3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr;
    Double I_KINETIC_Hx4z_G4x_vrr = PBX*I_KINETIC_Hx4z_F3x_vrr+oned2z*I_KINETIC_G4z_F3x_vrr+3*oned2z*I_KINETIC_Hx4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_KINETIC_H5y_G4x_vrr = PBX*I_KINETIC_H5y_F3x_vrr+3*oned2z*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_KINETIC_H4yz_G4x_vrr = PBX*I_KINETIC_H4yz_F3x_vrr+3*oned2z*I_KINETIC_H4yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_KINETIC_H3y2z_G4x_vrr = PBX*I_KINETIC_H3y2z_F3x_vrr+3*oned2z*I_KINETIC_H3y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H3y2z_D2x_vrr;
    Double I_KINETIC_H2y3z_G4x_vrr = PBX*I_KINETIC_H2y3z_F3x_vrr+3*oned2z*I_KINETIC_H2y3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H2y3z_D2x_vrr;
    Double I_KINETIC_Hy4z_G4x_vrr = PBX*I_KINETIC_Hy4z_F3x_vrr+3*oned2z*I_KINETIC_Hy4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_KINETIC_H5z_G4x_vrr = PBX*I_KINETIC_H5z_F3x_vrr+3*oned2z*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G4x_vrr-3*adz*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_KINETIC_H5x_G3xy_vrr = PBY*I_KINETIC_H5x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G3xy_vrr;
    Double I_KINETIC_H4xy_G3xy_vrr = PBY*I_KINETIC_H4xy_F3x_vrr+oned2z*I_KINETIC_G4x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G3xy_vrr;
    Double I_KINETIC_H4xz_G3xy_vrr = PBY*I_KINETIC_H4xz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G3xy_vrr;
    Double I_KINETIC_H3x2y_G3xy_vrr = PBY*I_KINETIC_H3x2y_F3x_vrr+2*oned2z*I_KINETIC_G3xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G3xy_vrr;
    Double I_KINETIC_H3xyz_G3xy_vrr = PBY*I_KINETIC_H3xyz_F3x_vrr+oned2z*I_KINETIC_G3xz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G3xy_vrr;
    Double I_KINETIC_H3x2z_G3xy_vrr = PBY*I_KINETIC_H3x2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G3xy_vrr;
    Double I_KINETIC_H2x3y_G3xy_vrr = PBY*I_KINETIC_H2x3y_F3x_vrr+3*oned2z*I_KINETIC_G2x2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G3xy_vrr;
    Double I_KINETIC_H2x2yz_G3xy_vrr = PBY*I_KINETIC_H2x2yz_F3x_vrr+2*oned2z*I_KINETIC_G2xyz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G3xy_vrr;
    Double I_KINETIC_H2xy2z_G3xy_vrr = PBY*I_KINETIC_H2xy2z_F3x_vrr+oned2z*I_KINETIC_G2x2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G3xy_vrr;
    Double I_KINETIC_H2x3z_G3xy_vrr = PBY*I_KINETIC_H2x3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G3xy_vrr;
    Double I_KINETIC_Hx4y_G3xy_vrr = PBY*I_KINETIC_Hx4y_F3x_vrr+4*oned2z*I_KINETIC_Gx3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G3xy_vrr;
    Double I_KINETIC_Hx3yz_G3xy_vrr = PBY*I_KINETIC_Hx3yz_F3x_vrr+3*oned2z*I_KINETIC_Gx2yz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G3xy_vrr;
    Double I_KINETIC_Hx2y2z_G3xy_vrr = PBY*I_KINETIC_Hx2y2z_F3x_vrr+2*oned2z*I_KINETIC_Gxy2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G3xy_vrr;
    Double I_KINETIC_Hxy3z_G3xy_vrr = PBY*I_KINETIC_Hxy3z_F3x_vrr+oned2z*I_KINETIC_Gx3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G3xy_vrr;
    Double I_KINETIC_Hx4z_G3xy_vrr = PBY*I_KINETIC_Hx4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G3xy_vrr;
    Double I_KINETIC_H5y_G3xy_vrr = PBY*I_KINETIC_H5y_F3x_vrr+5*oned2z*I_KINETIC_G4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G3xy_vrr;
    Double I_KINETIC_H4yz_G3xy_vrr = PBY*I_KINETIC_H4yz_F3x_vrr+4*oned2z*I_KINETIC_G3yz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G3xy_vrr;
    Double I_KINETIC_H3y2z_G3xy_vrr = PBY*I_KINETIC_H3y2z_F3x_vrr+3*oned2z*I_KINETIC_G2y2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G3xy_vrr;
    Double I_KINETIC_H2y3z_G3xy_vrr = PBY*I_KINETIC_H2y3z_F3x_vrr+2*oned2z*I_KINETIC_Gy3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G3xy_vrr;
    Double I_KINETIC_Hy4z_G3xy_vrr = PBY*I_KINETIC_Hy4z_F3x_vrr+oned2z*I_KINETIC_G4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G3xy_vrr;
    Double I_KINETIC_H5z_G3xy_vrr = PBY*I_KINETIC_H5z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G3xy_vrr;
    Double I_KINETIC_H5x_G3xz_vrr = PBZ*I_KINETIC_H5x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G3xz_vrr;
    Double I_KINETIC_H4xy_G3xz_vrr = PBZ*I_KINETIC_H4xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G3xz_vrr;
    Double I_KINETIC_H4xz_G3xz_vrr = PBZ*I_KINETIC_H4xz_F3x_vrr+oned2z*I_KINETIC_G4x_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G3xz_vrr;
    Double I_KINETIC_H3x2y_G3xz_vrr = PBZ*I_KINETIC_H3x2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G3xz_vrr;
    Double I_KINETIC_H3xyz_G3xz_vrr = PBZ*I_KINETIC_H3xyz_F3x_vrr+oned2z*I_KINETIC_G3xy_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G3xz_vrr;
    Double I_KINETIC_H3x2z_G3xz_vrr = PBZ*I_KINETIC_H3x2z_F3x_vrr+2*oned2z*I_KINETIC_G3xz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G3xz_vrr;
    Double I_KINETIC_H2x3y_G3xz_vrr = PBZ*I_KINETIC_H2x3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G3xz_vrr;
    Double I_KINETIC_H2x2yz_G3xz_vrr = PBZ*I_KINETIC_H2x2yz_F3x_vrr+oned2z*I_KINETIC_G2x2y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G3xz_vrr;
    Double I_KINETIC_H2xy2z_G3xz_vrr = PBZ*I_KINETIC_H2xy2z_F3x_vrr+2*oned2z*I_KINETIC_G2xyz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G3xz_vrr;
    Double I_KINETIC_H2x3z_G3xz_vrr = PBZ*I_KINETIC_H2x3z_F3x_vrr+3*oned2z*I_KINETIC_G2x2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G3xz_vrr;
    Double I_KINETIC_Hx4y_G3xz_vrr = PBZ*I_KINETIC_Hx4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G3xz_vrr;
    Double I_KINETIC_Hx3yz_G3xz_vrr = PBZ*I_KINETIC_Hx3yz_F3x_vrr+oned2z*I_KINETIC_Gx3y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G3xz_vrr;
    Double I_KINETIC_Hx2y2z_G3xz_vrr = PBZ*I_KINETIC_Hx2y2z_F3x_vrr+2*oned2z*I_KINETIC_Gx2yz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G3xz_vrr;
    Double I_KINETIC_Hxy3z_G3xz_vrr = PBZ*I_KINETIC_Hxy3z_F3x_vrr+3*oned2z*I_KINETIC_Gxy2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G3xz_vrr;
    Double I_KINETIC_Hx4z_G3xz_vrr = PBZ*I_KINETIC_Hx4z_F3x_vrr+4*oned2z*I_KINETIC_Gx3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G3xz_vrr;
    Double I_KINETIC_H5y_G3xz_vrr = PBZ*I_KINETIC_H5y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G3xz_vrr;
    Double I_KINETIC_H4yz_G3xz_vrr = PBZ*I_KINETIC_H4yz_F3x_vrr+oned2z*I_KINETIC_G4y_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G3xz_vrr;
    Double I_KINETIC_H3y2z_G3xz_vrr = PBZ*I_KINETIC_H3y2z_F3x_vrr+2*oned2z*I_KINETIC_G3yz_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G3xz_vrr;
    Double I_KINETIC_H2y3z_G3xz_vrr = PBZ*I_KINETIC_H2y3z_F3x_vrr+3*oned2z*I_KINETIC_G2y2z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G3xz_vrr;
    Double I_KINETIC_Hy4z_G3xz_vrr = PBZ*I_KINETIC_Hy4z_F3x_vrr+4*oned2z*I_KINETIC_Gy3z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G3xz_vrr;
    Double I_KINETIC_H5z_G3xz_vrr = PBZ*I_KINETIC_H5z_F3x_vrr+5*oned2z*I_KINETIC_G4z_F3x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G3xz_vrr;
    Double I_KINETIC_H5x_G2x2y_vrr = PBY*I_KINETIC_H5x_F2xy_vrr+oned2z*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_KINETIC_H4xy_G2x2y_vrr = PBY*I_KINETIC_H4xy_F2xy_vrr+oned2z*I_KINETIC_G4x_F2xy_vrr+oned2z*I_KINETIC_H4xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_KINETIC_H4xz_G2x2y_vrr = PBY*I_KINETIC_H4xz_F2xy_vrr+oned2z*I_KINETIC_H4xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_KINETIC_H3x2y_G2x2y_vrr = PBY*I_KINETIC_H3x2y_F2xy_vrr+2*oned2z*I_KINETIC_G3xy_F2xy_vrr+oned2z*I_KINETIC_H3x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_KINETIC_H3xyz_G2x2y_vrr = PBY*I_KINETIC_H3xyz_F2xy_vrr+oned2z*I_KINETIC_G3xz_F2xy_vrr+oned2z*I_KINETIC_H3xyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H3xyz_D2x_vrr;
    Double I_KINETIC_H3x2z_G2x2y_vrr = PBY*I_KINETIC_H3x2z_F2xy_vrr+oned2z*I_KINETIC_H3x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_KINETIC_H2x3y_G2x2y_vrr = PBY*I_KINETIC_H2x3y_F2xy_vrr+3*oned2z*I_KINETIC_G2x2y_F2xy_vrr+oned2z*I_KINETIC_H2x3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_KINETIC_H2x2yz_G2x2y_vrr = PBY*I_KINETIC_H2x2yz_F2xy_vrr+2*oned2z*I_KINETIC_G2xyz_F2xy_vrr+oned2z*I_KINETIC_H2x2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr;
    Double I_KINETIC_H2xy2z_G2x2y_vrr = PBY*I_KINETIC_H2xy2z_F2xy_vrr+oned2z*I_KINETIC_G2x2z_F2xy_vrr+oned2z*I_KINETIC_H2xy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr;
    Double I_KINETIC_H2x3z_G2x2y_vrr = PBY*I_KINETIC_H2x3z_F2xy_vrr+oned2z*I_KINETIC_H2x3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_KINETIC_Hx4y_G2x2y_vrr = PBY*I_KINETIC_Hx4y_F2xy_vrr+4*oned2z*I_KINETIC_Gx3y_F2xy_vrr+oned2z*I_KINETIC_Hx4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_KINETIC_Hx3yz_G2x2y_vrr = PBY*I_KINETIC_Hx3yz_F2xy_vrr+3*oned2z*I_KINETIC_Gx2yz_F2xy_vrr+oned2z*I_KINETIC_Hx3yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr;
    Double I_KINETIC_Hx2y2z_G2x2y_vrr = PBY*I_KINETIC_Hx2y2z_F2xy_vrr+2*oned2z*I_KINETIC_Gxy2z_F2xy_vrr+oned2z*I_KINETIC_Hx2y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr;
    Double I_KINETIC_Hxy3z_G2x2y_vrr = PBY*I_KINETIC_Hxy3z_F2xy_vrr+oned2z*I_KINETIC_Gx3z_F2xy_vrr+oned2z*I_KINETIC_Hxy3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr;
    Double I_KINETIC_Hx4z_G2x2y_vrr = PBY*I_KINETIC_Hx4z_F2xy_vrr+oned2z*I_KINETIC_Hx4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_KINETIC_H5y_G2x2y_vrr = PBY*I_KINETIC_H5y_F2xy_vrr+5*oned2z*I_KINETIC_G4y_F2xy_vrr+oned2z*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_KINETIC_H4yz_G2x2y_vrr = PBY*I_KINETIC_H4yz_F2xy_vrr+4*oned2z*I_KINETIC_G3yz_F2xy_vrr+oned2z*I_KINETIC_H4yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_KINETIC_H3y2z_G2x2y_vrr = PBY*I_KINETIC_H3y2z_F2xy_vrr+3*oned2z*I_KINETIC_G2y2z_F2xy_vrr+oned2z*I_KINETIC_H3y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H3y2z_D2x_vrr;
    Double I_KINETIC_H2y3z_G2x2y_vrr = PBY*I_KINETIC_H2y3z_F2xy_vrr+2*oned2z*I_KINETIC_Gy3z_F2xy_vrr+oned2z*I_KINETIC_H2y3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H2y3z_D2x_vrr;
    Double I_KINETIC_Hy4z_G2x2y_vrr = PBY*I_KINETIC_Hy4z_F2xy_vrr+oned2z*I_KINETIC_G4z_F2xy_vrr+oned2z*I_KINETIC_Hy4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_KINETIC_H5z_G2x2y_vrr = PBY*I_KINETIC_H5z_F2xy_vrr+oned2z*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G2x2y_vrr-adz*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_KINETIC_H5x_G2xyz_vrr = PBZ*I_KINETIC_H5x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G2xyz_vrr;
    Double I_KINETIC_H4xy_G2xyz_vrr = PBZ*I_KINETIC_H4xy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G2xyz_vrr;
    Double I_KINETIC_H4xz_G2xyz_vrr = PBZ*I_KINETIC_H4xz_F2xy_vrr+oned2z*I_KINETIC_G4x_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G2xyz_vrr;
    Double I_KINETIC_H3x2y_G2xyz_vrr = PBZ*I_KINETIC_H3x2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G2xyz_vrr;
    Double I_KINETIC_H3xyz_G2xyz_vrr = PBZ*I_KINETIC_H3xyz_F2xy_vrr+oned2z*I_KINETIC_G3xy_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G2xyz_vrr;
    Double I_KINETIC_H3x2z_G2xyz_vrr = PBZ*I_KINETIC_H3x2z_F2xy_vrr+2*oned2z*I_KINETIC_G3xz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G2xyz_vrr;
    Double I_KINETIC_H2x3y_G2xyz_vrr = PBZ*I_KINETIC_H2x3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G2xyz_vrr;
    Double I_KINETIC_H2x2yz_G2xyz_vrr = PBZ*I_KINETIC_H2x2yz_F2xy_vrr+oned2z*I_KINETIC_G2x2y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G2xyz_vrr;
    Double I_KINETIC_H2xy2z_G2xyz_vrr = PBZ*I_KINETIC_H2xy2z_F2xy_vrr+2*oned2z*I_KINETIC_G2xyz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G2xyz_vrr;
    Double I_KINETIC_H2x3z_G2xyz_vrr = PBZ*I_KINETIC_H2x3z_F2xy_vrr+3*oned2z*I_KINETIC_G2x2z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G2xyz_vrr;
    Double I_KINETIC_Hx4y_G2xyz_vrr = PBZ*I_KINETIC_Hx4y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G2xyz_vrr;
    Double I_KINETIC_Hx3yz_G2xyz_vrr = PBZ*I_KINETIC_Hx3yz_F2xy_vrr+oned2z*I_KINETIC_Gx3y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G2xyz_vrr;
    Double I_KINETIC_Hx2y2z_G2xyz_vrr = PBZ*I_KINETIC_Hx2y2z_F2xy_vrr+2*oned2z*I_KINETIC_Gx2yz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G2xyz_vrr;
    Double I_KINETIC_Hxy3z_G2xyz_vrr = PBZ*I_KINETIC_Hxy3z_F2xy_vrr+3*oned2z*I_KINETIC_Gxy2z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G2xyz_vrr;
    Double I_KINETIC_Hx4z_G2xyz_vrr = PBZ*I_KINETIC_Hx4z_F2xy_vrr+4*oned2z*I_KINETIC_Gx3z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G2xyz_vrr;
    Double I_KINETIC_H5y_G2xyz_vrr = PBZ*I_KINETIC_H5y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G2xyz_vrr;
    Double I_KINETIC_H4yz_G2xyz_vrr = PBZ*I_KINETIC_H4yz_F2xy_vrr+oned2z*I_KINETIC_G4y_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G2xyz_vrr;
    Double I_KINETIC_H3y2z_G2xyz_vrr = PBZ*I_KINETIC_H3y2z_F2xy_vrr+2*oned2z*I_KINETIC_G3yz_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G2xyz_vrr;
    Double I_KINETIC_H2y3z_G2xyz_vrr = PBZ*I_KINETIC_H2y3z_F2xy_vrr+3*oned2z*I_KINETIC_G2y2z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G2xyz_vrr;
    Double I_KINETIC_Hy4z_G2xyz_vrr = PBZ*I_KINETIC_Hy4z_F2xy_vrr+4*oned2z*I_KINETIC_Gy3z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G2xyz_vrr;
    Double I_KINETIC_H5z_G2xyz_vrr = PBZ*I_KINETIC_H5z_F2xy_vrr+5*oned2z*I_KINETIC_G4z_F2xy_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G2xyz_vrr;
    Double I_KINETIC_H5x_G2x2z_vrr = PBZ*I_KINETIC_H5x_F2xz_vrr+oned2z*I_KINETIC_H5x_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H5x_D2x_vrr;
    Double I_KINETIC_H4xy_G2x2z_vrr = PBZ*I_KINETIC_H4xy_F2xz_vrr+oned2z*I_KINETIC_H4xy_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H4xy_D2x_vrr;
    Double I_KINETIC_H4xz_G2x2z_vrr = PBZ*I_KINETIC_H4xz_F2xz_vrr+oned2z*I_KINETIC_G4x_F2xz_vrr+oned2z*I_KINETIC_H4xz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H4xz_D2x_vrr;
    Double I_KINETIC_H3x2y_G2x2z_vrr = PBZ*I_KINETIC_H3x2y_F2xz_vrr+oned2z*I_KINETIC_H3x2y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H3x2y_D2x_vrr;
    Double I_KINETIC_H3xyz_G2x2z_vrr = PBZ*I_KINETIC_H3xyz_F2xz_vrr+oned2z*I_KINETIC_G3xy_F2xz_vrr+oned2z*I_KINETIC_H3xyz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H3xyz_D2x_vrr;
    Double I_KINETIC_H3x2z_G2x2z_vrr = PBZ*I_KINETIC_H3x2z_F2xz_vrr+2*oned2z*I_KINETIC_G3xz_F2xz_vrr+oned2z*I_KINETIC_H3x2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H3x2z_D2x_vrr;
    Double I_KINETIC_H2x3y_G2x2z_vrr = PBZ*I_KINETIC_H2x3y_F2xz_vrr+oned2z*I_KINETIC_H2x3y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H2x3y_D2x_vrr;
    Double I_KINETIC_H2x2yz_G2x2z_vrr = PBZ*I_KINETIC_H2x2yz_F2xz_vrr+oned2z*I_KINETIC_G2x2y_F2xz_vrr+oned2z*I_KINETIC_H2x2yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_D2x_vrr;
    Double I_KINETIC_H2xy2z_G2x2z_vrr = PBZ*I_KINETIC_H2xy2z_F2xz_vrr+2*oned2z*I_KINETIC_G2xyz_F2xz_vrr+oned2z*I_KINETIC_H2xy2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H2xy2z_D2x_vrr;
    Double I_KINETIC_H2x3z_G2x2z_vrr = PBZ*I_KINETIC_H2x3z_F2xz_vrr+3*oned2z*I_KINETIC_G2x2z_F2xz_vrr+oned2z*I_KINETIC_H2x3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H2x3z_D2x_vrr;
    Double I_KINETIC_Hx4y_G2x2z_vrr = PBZ*I_KINETIC_Hx4y_F2xz_vrr+oned2z*I_KINETIC_Hx4y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Hx4y_D2x_vrr;
    Double I_KINETIC_Hx3yz_G2x2z_vrr = PBZ*I_KINETIC_Hx3yz_F2xz_vrr+oned2z*I_KINETIC_Gx3y_F2xz_vrr+oned2z*I_KINETIC_Hx3yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Hx3yz_D2x_vrr;
    Double I_KINETIC_Hx2y2z_G2x2z_vrr = PBZ*I_KINETIC_Hx2y2z_F2xz_vrr+2*oned2z*I_KINETIC_Gx2yz_F2xz_vrr+oned2z*I_KINETIC_Hx2y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Hx2y2z_D2x_vrr;
    Double I_KINETIC_Hxy3z_G2x2z_vrr = PBZ*I_KINETIC_Hxy3z_F2xz_vrr+3*oned2z*I_KINETIC_Gxy2z_F2xz_vrr+oned2z*I_KINETIC_Hxy3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Hxy3z_D2x_vrr;
    Double I_KINETIC_Hx4z_G2x2z_vrr = PBZ*I_KINETIC_Hx4z_F2xz_vrr+4*oned2z*I_KINETIC_Gx3z_F2xz_vrr+oned2z*I_KINETIC_Hx4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Hx4z_D2x_vrr;
    Double I_KINETIC_H5y_G2x2z_vrr = PBZ*I_KINETIC_H5y_F2xz_vrr+oned2z*I_KINETIC_H5y_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H5y_D2x_vrr;
    Double I_KINETIC_H4yz_G2x2z_vrr = PBZ*I_KINETIC_H4yz_F2xz_vrr+oned2z*I_KINETIC_G4y_F2xz_vrr+oned2z*I_KINETIC_H4yz_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H4yz_D2x_vrr;
    Double I_KINETIC_H3y2z_G2x2z_vrr = PBZ*I_KINETIC_H3y2z_F2xz_vrr+2*oned2z*I_KINETIC_G3yz_F2xz_vrr+oned2z*I_KINETIC_H3y2z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H3y2z_D2x_vrr;
    Double I_KINETIC_H2y3z_G2x2z_vrr = PBZ*I_KINETIC_H2y3z_F2xz_vrr+3*oned2z*I_KINETIC_G2y2z_F2xz_vrr+oned2z*I_KINETIC_H2y3z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H2y3z_D2x_vrr;
    Double I_KINETIC_Hy4z_G2x2z_vrr = PBZ*I_KINETIC_Hy4z_F2xz_vrr+4*oned2z*I_KINETIC_Gy3z_F2xz_vrr+oned2z*I_KINETIC_Hy4z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_Hy4z_D2x_vrr;
    Double I_KINETIC_H5z_G2x2z_vrr = PBZ*I_KINETIC_H5z_F2xz_vrr+5*oned2z*I_KINETIC_G4z_F2xz_vrr+oned2z*I_KINETIC_H5z_D2x_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G2x2z_vrr-adz*I_TWOBODYOVERLAP_H5z_D2x_vrr;
    Double I_KINETIC_H5x_Gx3y_vrr = PBX*I_KINETIC_H5x_F3y_vrr+5*oned2z*I_KINETIC_G4x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gx3y_vrr;
    Double I_KINETIC_H4xy_Gx3y_vrr = PBX*I_KINETIC_H4xy_F3y_vrr+4*oned2z*I_KINETIC_G3xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gx3y_vrr;
    Double I_KINETIC_H4xz_Gx3y_vrr = PBX*I_KINETIC_H4xz_F3y_vrr+4*oned2z*I_KINETIC_G3xz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gx3y_vrr;
    Double I_KINETIC_H3x2y_Gx3y_vrr = PBX*I_KINETIC_H3x2y_F3y_vrr+3*oned2z*I_KINETIC_G2x2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gx3y_vrr;
    Double I_KINETIC_H3xyz_Gx3y_vrr = PBX*I_KINETIC_H3xyz_F3y_vrr+3*oned2z*I_KINETIC_G2xyz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gx3y_vrr;
    Double I_KINETIC_H3x2z_Gx3y_vrr = PBX*I_KINETIC_H3x2z_F3y_vrr+3*oned2z*I_KINETIC_G2x2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gx3y_vrr;
    Double I_KINETIC_H2x3y_Gx3y_vrr = PBX*I_KINETIC_H2x3y_F3y_vrr+2*oned2z*I_KINETIC_Gx3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gx3y_vrr;
    Double I_KINETIC_H2x2yz_Gx3y_vrr = PBX*I_KINETIC_H2x2yz_F3y_vrr+2*oned2z*I_KINETIC_Gx2yz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gx3y_vrr;
    Double I_KINETIC_H2xy2z_Gx3y_vrr = PBX*I_KINETIC_H2xy2z_F3y_vrr+2*oned2z*I_KINETIC_Gxy2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gx3y_vrr;
    Double I_KINETIC_H2x3z_Gx3y_vrr = PBX*I_KINETIC_H2x3z_F3y_vrr+2*oned2z*I_KINETIC_Gx3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gx3y_vrr;
    Double I_KINETIC_Hx4y_Gx3y_vrr = PBX*I_KINETIC_Hx4y_F3y_vrr+oned2z*I_KINETIC_G4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gx3y_vrr;
    Double I_KINETIC_Hx3yz_Gx3y_vrr = PBX*I_KINETIC_Hx3yz_F3y_vrr+oned2z*I_KINETIC_G3yz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gx3y_vrr;
    Double I_KINETIC_Hx2y2z_Gx3y_vrr = PBX*I_KINETIC_Hx2y2z_F3y_vrr+oned2z*I_KINETIC_G2y2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gx3y_vrr;
    Double I_KINETIC_Hxy3z_Gx3y_vrr = PBX*I_KINETIC_Hxy3z_F3y_vrr+oned2z*I_KINETIC_Gy3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gx3y_vrr;
    Double I_KINETIC_Hx4z_Gx3y_vrr = PBX*I_KINETIC_Hx4z_F3y_vrr+oned2z*I_KINETIC_G4z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gx3y_vrr;
    Double I_KINETIC_H5y_Gx3y_vrr = PBX*I_KINETIC_H5y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gx3y_vrr;
    Double I_KINETIC_H4yz_Gx3y_vrr = PBX*I_KINETIC_H4yz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gx3y_vrr;
    Double I_KINETIC_H3y2z_Gx3y_vrr = PBX*I_KINETIC_H3y2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gx3y_vrr;
    Double I_KINETIC_H2y3z_Gx3y_vrr = PBX*I_KINETIC_H2y3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gx3y_vrr;
    Double I_KINETIC_Hy4z_Gx3y_vrr = PBX*I_KINETIC_Hy4z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gx3y_vrr;
    Double I_KINETIC_H5z_Gx3y_vrr = PBX*I_KINETIC_H5z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gx3y_vrr;
    Double I_KINETIC_H5x_Gx2yz_vrr = PBZ*I_KINETIC_H5x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gx2yz_vrr;
    Double I_KINETIC_H4xy_Gx2yz_vrr = PBZ*I_KINETIC_H4xy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gx2yz_vrr;
    Double I_KINETIC_H4xz_Gx2yz_vrr = PBZ*I_KINETIC_H4xz_Fx2y_vrr+oned2z*I_KINETIC_G4x_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gx2yz_vrr;
    Double I_KINETIC_H3x2y_Gx2yz_vrr = PBZ*I_KINETIC_H3x2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gx2yz_vrr;
    Double I_KINETIC_H3xyz_Gx2yz_vrr = PBZ*I_KINETIC_H3xyz_Fx2y_vrr+oned2z*I_KINETIC_G3xy_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gx2yz_vrr;
    Double I_KINETIC_H3x2z_Gx2yz_vrr = PBZ*I_KINETIC_H3x2z_Fx2y_vrr+2*oned2z*I_KINETIC_G3xz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gx2yz_vrr;
    Double I_KINETIC_H2x3y_Gx2yz_vrr = PBZ*I_KINETIC_H2x3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gx2yz_vrr;
    Double I_KINETIC_H2x2yz_Gx2yz_vrr = PBZ*I_KINETIC_H2x2yz_Fx2y_vrr+oned2z*I_KINETIC_G2x2y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gx2yz_vrr;
    Double I_KINETIC_H2xy2z_Gx2yz_vrr = PBZ*I_KINETIC_H2xy2z_Fx2y_vrr+2*oned2z*I_KINETIC_G2xyz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gx2yz_vrr;
    Double I_KINETIC_H2x3z_Gx2yz_vrr = PBZ*I_KINETIC_H2x3z_Fx2y_vrr+3*oned2z*I_KINETIC_G2x2z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gx2yz_vrr;
    Double I_KINETIC_Hx4y_Gx2yz_vrr = PBZ*I_KINETIC_Hx4y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gx2yz_vrr;
    Double I_KINETIC_Hx3yz_Gx2yz_vrr = PBZ*I_KINETIC_Hx3yz_Fx2y_vrr+oned2z*I_KINETIC_Gx3y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gx2yz_vrr;
    Double I_KINETIC_Hx2y2z_Gx2yz_vrr = PBZ*I_KINETIC_Hx2y2z_Fx2y_vrr+2*oned2z*I_KINETIC_Gx2yz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gx2yz_vrr;
    Double I_KINETIC_Hxy3z_Gx2yz_vrr = PBZ*I_KINETIC_Hxy3z_Fx2y_vrr+3*oned2z*I_KINETIC_Gxy2z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gx2yz_vrr;
    Double I_KINETIC_Hx4z_Gx2yz_vrr = PBZ*I_KINETIC_Hx4z_Fx2y_vrr+4*oned2z*I_KINETIC_Gx3z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gx2yz_vrr;
    Double I_KINETIC_H5y_Gx2yz_vrr = PBZ*I_KINETIC_H5y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gx2yz_vrr;
    Double I_KINETIC_H4yz_Gx2yz_vrr = PBZ*I_KINETIC_H4yz_Fx2y_vrr+oned2z*I_KINETIC_G4y_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gx2yz_vrr;
    Double I_KINETIC_H3y2z_Gx2yz_vrr = PBZ*I_KINETIC_H3y2z_Fx2y_vrr+2*oned2z*I_KINETIC_G3yz_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gx2yz_vrr;
    Double I_KINETIC_H2y3z_Gx2yz_vrr = PBZ*I_KINETIC_H2y3z_Fx2y_vrr+3*oned2z*I_KINETIC_G2y2z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gx2yz_vrr;
    Double I_KINETIC_Hy4z_Gx2yz_vrr = PBZ*I_KINETIC_Hy4z_Fx2y_vrr+4*oned2z*I_KINETIC_Gy3z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gx2yz_vrr;
    Double I_KINETIC_H5z_Gx2yz_vrr = PBZ*I_KINETIC_H5z_Fx2y_vrr+5*oned2z*I_KINETIC_G4z_Fx2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gx2yz_vrr;
    Double I_KINETIC_H5x_Gxy2z_vrr = PBY*I_KINETIC_H5x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gxy2z_vrr;
    Double I_KINETIC_H4xy_Gxy2z_vrr = PBY*I_KINETIC_H4xy_Fx2z_vrr+oned2z*I_KINETIC_G4x_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gxy2z_vrr;
    Double I_KINETIC_H4xz_Gxy2z_vrr = PBY*I_KINETIC_H4xz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gxy2z_vrr;
    Double I_KINETIC_H3x2y_Gxy2z_vrr = PBY*I_KINETIC_H3x2y_Fx2z_vrr+2*oned2z*I_KINETIC_G3xy_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gxy2z_vrr;
    Double I_KINETIC_H3xyz_Gxy2z_vrr = PBY*I_KINETIC_H3xyz_Fx2z_vrr+oned2z*I_KINETIC_G3xz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gxy2z_vrr;
    Double I_KINETIC_H3x2z_Gxy2z_vrr = PBY*I_KINETIC_H3x2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gxy2z_vrr;
    Double I_KINETIC_H2x3y_Gxy2z_vrr = PBY*I_KINETIC_H2x3y_Fx2z_vrr+3*oned2z*I_KINETIC_G2x2y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gxy2z_vrr;
    Double I_KINETIC_H2x2yz_Gxy2z_vrr = PBY*I_KINETIC_H2x2yz_Fx2z_vrr+2*oned2z*I_KINETIC_G2xyz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gxy2z_vrr;
    Double I_KINETIC_H2xy2z_Gxy2z_vrr = PBY*I_KINETIC_H2xy2z_Fx2z_vrr+oned2z*I_KINETIC_G2x2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gxy2z_vrr;
    Double I_KINETIC_H2x3z_Gxy2z_vrr = PBY*I_KINETIC_H2x3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gxy2z_vrr;
    Double I_KINETIC_Hx4y_Gxy2z_vrr = PBY*I_KINETIC_Hx4y_Fx2z_vrr+4*oned2z*I_KINETIC_Gx3y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gxy2z_vrr;
    Double I_KINETIC_Hx3yz_Gxy2z_vrr = PBY*I_KINETIC_Hx3yz_Fx2z_vrr+3*oned2z*I_KINETIC_Gx2yz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gxy2z_vrr;
    Double I_KINETIC_Hx2y2z_Gxy2z_vrr = PBY*I_KINETIC_Hx2y2z_Fx2z_vrr+2*oned2z*I_KINETIC_Gxy2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gxy2z_vrr;
    Double I_KINETIC_Hxy3z_Gxy2z_vrr = PBY*I_KINETIC_Hxy3z_Fx2z_vrr+oned2z*I_KINETIC_Gx3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gxy2z_vrr;
    Double I_KINETIC_Hx4z_Gxy2z_vrr = PBY*I_KINETIC_Hx4z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gxy2z_vrr;
    Double I_KINETIC_H5y_Gxy2z_vrr = PBY*I_KINETIC_H5y_Fx2z_vrr+5*oned2z*I_KINETIC_G4y_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gxy2z_vrr;
    Double I_KINETIC_H4yz_Gxy2z_vrr = PBY*I_KINETIC_H4yz_Fx2z_vrr+4*oned2z*I_KINETIC_G3yz_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gxy2z_vrr;
    Double I_KINETIC_H3y2z_Gxy2z_vrr = PBY*I_KINETIC_H3y2z_Fx2z_vrr+3*oned2z*I_KINETIC_G2y2z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gxy2z_vrr;
    Double I_KINETIC_H2y3z_Gxy2z_vrr = PBY*I_KINETIC_H2y3z_Fx2z_vrr+2*oned2z*I_KINETIC_Gy3z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gxy2z_vrr;
    Double I_KINETIC_Hy4z_Gxy2z_vrr = PBY*I_KINETIC_Hy4z_Fx2z_vrr+oned2z*I_KINETIC_G4z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gxy2z_vrr;
    Double I_KINETIC_H5z_Gxy2z_vrr = PBY*I_KINETIC_H5z_Fx2z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gxy2z_vrr;
    Double I_KINETIC_H5x_Gx3z_vrr = PBX*I_KINETIC_H5x_F3z_vrr+5*oned2z*I_KINETIC_G4x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gx3z_vrr;
    Double I_KINETIC_H4xy_Gx3z_vrr = PBX*I_KINETIC_H4xy_F3z_vrr+4*oned2z*I_KINETIC_G3xy_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gx3z_vrr;
    Double I_KINETIC_H4xz_Gx3z_vrr = PBX*I_KINETIC_H4xz_F3z_vrr+4*oned2z*I_KINETIC_G3xz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gx3z_vrr;
    Double I_KINETIC_H3x2y_Gx3z_vrr = PBX*I_KINETIC_H3x2y_F3z_vrr+3*oned2z*I_KINETIC_G2x2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gx3z_vrr;
    Double I_KINETIC_H3xyz_Gx3z_vrr = PBX*I_KINETIC_H3xyz_F3z_vrr+3*oned2z*I_KINETIC_G2xyz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gx3z_vrr;
    Double I_KINETIC_H3x2z_Gx3z_vrr = PBX*I_KINETIC_H3x2z_F3z_vrr+3*oned2z*I_KINETIC_G2x2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gx3z_vrr;
    Double I_KINETIC_H2x3y_Gx3z_vrr = PBX*I_KINETIC_H2x3y_F3z_vrr+2*oned2z*I_KINETIC_Gx3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gx3z_vrr;
    Double I_KINETIC_H2x2yz_Gx3z_vrr = PBX*I_KINETIC_H2x2yz_F3z_vrr+2*oned2z*I_KINETIC_Gx2yz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gx3z_vrr;
    Double I_KINETIC_H2xy2z_Gx3z_vrr = PBX*I_KINETIC_H2xy2z_F3z_vrr+2*oned2z*I_KINETIC_Gxy2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gx3z_vrr;
    Double I_KINETIC_H2x3z_Gx3z_vrr = PBX*I_KINETIC_H2x3z_F3z_vrr+2*oned2z*I_KINETIC_Gx3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gx3z_vrr;
    Double I_KINETIC_Hx4y_Gx3z_vrr = PBX*I_KINETIC_Hx4y_F3z_vrr+oned2z*I_KINETIC_G4y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gx3z_vrr;
    Double I_KINETIC_Hx3yz_Gx3z_vrr = PBX*I_KINETIC_Hx3yz_F3z_vrr+oned2z*I_KINETIC_G3yz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gx3z_vrr;
    Double I_KINETIC_Hx2y2z_Gx3z_vrr = PBX*I_KINETIC_Hx2y2z_F3z_vrr+oned2z*I_KINETIC_G2y2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gx3z_vrr;
    Double I_KINETIC_Hxy3z_Gx3z_vrr = PBX*I_KINETIC_Hxy3z_F3z_vrr+oned2z*I_KINETIC_Gy3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gx3z_vrr;
    Double I_KINETIC_Hx4z_Gx3z_vrr = PBX*I_KINETIC_Hx4z_F3z_vrr+oned2z*I_KINETIC_G4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gx3z_vrr;
    Double I_KINETIC_H5y_Gx3z_vrr = PBX*I_KINETIC_H5y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gx3z_vrr;
    Double I_KINETIC_H4yz_Gx3z_vrr = PBX*I_KINETIC_H4yz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gx3z_vrr;
    Double I_KINETIC_H3y2z_Gx3z_vrr = PBX*I_KINETIC_H3y2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gx3z_vrr;
    Double I_KINETIC_H2y3z_Gx3z_vrr = PBX*I_KINETIC_H2y3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gx3z_vrr;
    Double I_KINETIC_Hy4z_Gx3z_vrr = PBX*I_KINETIC_Hy4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gx3z_vrr;
    Double I_KINETIC_H5z_Gx3z_vrr = PBX*I_KINETIC_H5z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gx3z_vrr;
    Double I_KINETIC_H5x_G4y_vrr = PBY*I_KINETIC_H5x_F3y_vrr+3*oned2z*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_KINETIC_H4xy_G4y_vrr = PBY*I_KINETIC_H4xy_F3y_vrr+oned2z*I_KINETIC_G4x_F3y_vrr+3*oned2z*I_KINETIC_H4xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_KINETIC_H4xz_G4y_vrr = PBY*I_KINETIC_H4xz_F3y_vrr+3*oned2z*I_KINETIC_H4xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H4xz_D2y_vrr;
    Double I_KINETIC_H3x2y_G4y_vrr = PBY*I_KINETIC_H3x2y_F3y_vrr+2*oned2z*I_KINETIC_G3xy_F3y_vrr+3*oned2z*I_KINETIC_H3x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_KINETIC_H3xyz_G4y_vrr = PBY*I_KINETIC_H3xyz_F3y_vrr+oned2z*I_KINETIC_G3xz_F3y_vrr+3*oned2z*I_KINETIC_H3xyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H3xyz_D2y_vrr;
    Double I_KINETIC_H3x2z_G4y_vrr = PBY*I_KINETIC_H3x2z_F3y_vrr+3*oned2z*I_KINETIC_H3x2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H3x2z_D2y_vrr;
    Double I_KINETIC_H2x3y_G4y_vrr = PBY*I_KINETIC_H2x3y_F3y_vrr+3*oned2z*I_KINETIC_G2x2y_F3y_vrr+3*oned2z*I_KINETIC_H2x3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_KINETIC_H2x2yz_G4y_vrr = PBY*I_KINETIC_H2x2yz_F3y_vrr+2*oned2z*I_KINETIC_G2xyz_F3y_vrr+3*oned2z*I_KINETIC_H2x2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr;
    Double I_KINETIC_H2xy2z_G4y_vrr = PBY*I_KINETIC_H2xy2z_F3y_vrr+oned2z*I_KINETIC_G2x2z_F3y_vrr+3*oned2z*I_KINETIC_H2xy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr;
    Double I_KINETIC_H2x3z_G4y_vrr = PBY*I_KINETIC_H2x3z_F3y_vrr+3*oned2z*I_KINETIC_H2x3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H2x3z_D2y_vrr;
    Double I_KINETIC_Hx4y_G4y_vrr = PBY*I_KINETIC_Hx4y_F3y_vrr+4*oned2z*I_KINETIC_Gx3y_F3y_vrr+3*oned2z*I_KINETIC_Hx4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_KINETIC_Hx3yz_G4y_vrr = PBY*I_KINETIC_Hx3yz_F3y_vrr+3*oned2z*I_KINETIC_Gx2yz_F3y_vrr+3*oned2z*I_KINETIC_Hx3yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr;
    Double I_KINETIC_Hx2y2z_G4y_vrr = PBY*I_KINETIC_Hx2y2z_F3y_vrr+2*oned2z*I_KINETIC_Gxy2z_F3y_vrr+3*oned2z*I_KINETIC_Hx2y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr;
    Double I_KINETIC_Hxy3z_G4y_vrr = PBY*I_KINETIC_Hxy3z_F3y_vrr+oned2z*I_KINETIC_Gx3z_F3y_vrr+3*oned2z*I_KINETIC_Hxy3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr;
    Double I_KINETIC_Hx4z_G4y_vrr = PBY*I_KINETIC_Hx4z_F3y_vrr+3*oned2z*I_KINETIC_Hx4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_KINETIC_H5y_G4y_vrr = PBY*I_KINETIC_H5y_F3y_vrr+5*oned2z*I_KINETIC_G4y_F3y_vrr+3*oned2z*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_KINETIC_H4yz_G4y_vrr = PBY*I_KINETIC_H4yz_F3y_vrr+4*oned2z*I_KINETIC_G3yz_F3y_vrr+3*oned2z*I_KINETIC_H4yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_KINETIC_H3y2z_G4y_vrr = PBY*I_KINETIC_H3y2z_F3y_vrr+3*oned2z*I_KINETIC_G2y2z_F3y_vrr+3*oned2z*I_KINETIC_H3y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H3y2z_D2y_vrr;
    Double I_KINETIC_H2y3z_G4y_vrr = PBY*I_KINETIC_H2y3z_F3y_vrr+2*oned2z*I_KINETIC_Gy3z_F3y_vrr+3*oned2z*I_KINETIC_H2y3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H2y3z_D2y_vrr;
    Double I_KINETIC_Hy4z_G4y_vrr = PBY*I_KINETIC_Hy4z_F3y_vrr+oned2z*I_KINETIC_G4z_F3y_vrr+3*oned2z*I_KINETIC_Hy4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_KINETIC_H5z_G4y_vrr = PBY*I_KINETIC_H5z_F3y_vrr+3*oned2z*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G4y_vrr-3*adz*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_KINETIC_H5x_G3yz_vrr = PBZ*I_KINETIC_H5x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G3yz_vrr;
    Double I_KINETIC_H4xy_G3yz_vrr = PBZ*I_KINETIC_H4xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G3yz_vrr;
    Double I_KINETIC_H4xz_G3yz_vrr = PBZ*I_KINETIC_H4xz_F3y_vrr+oned2z*I_KINETIC_G4x_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G3yz_vrr;
    Double I_KINETIC_H3x2y_G3yz_vrr = PBZ*I_KINETIC_H3x2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G3yz_vrr;
    Double I_KINETIC_H3xyz_G3yz_vrr = PBZ*I_KINETIC_H3xyz_F3y_vrr+oned2z*I_KINETIC_G3xy_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G3yz_vrr;
    Double I_KINETIC_H3x2z_G3yz_vrr = PBZ*I_KINETIC_H3x2z_F3y_vrr+2*oned2z*I_KINETIC_G3xz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G3yz_vrr;
    Double I_KINETIC_H2x3y_G3yz_vrr = PBZ*I_KINETIC_H2x3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G3yz_vrr;
    Double I_KINETIC_H2x2yz_G3yz_vrr = PBZ*I_KINETIC_H2x2yz_F3y_vrr+oned2z*I_KINETIC_G2x2y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G3yz_vrr;
    Double I_KINETIC_H2xy2z_G3yz_vrr = PBZ*I_KINETIC_H2xy2z_F3y_vrr+2*oned2z*I_KINETIC_G2xyz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G3yz_vrr;
    Double I_KINETIC_H2x3z_G3yz_vrr = PBZ*I_KINETIC_H2x3z_F3y_vrr+3*oned2z*I_KINETIC_G2x2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G3yz_vrr;
    Double I_KINETIC_Hx4y_G3yz_vrr = PBZ*I_KINETIC_Hx4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G3yz_vrr;
    Double I_KINETIC_Hx3yz_G3yz_vrr = PBZ*I_KINETIC_Hx3yz_F3y_vrr+oned2z*I_KINETIC_Gx3y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G3yz_vrr;
    Double I_KINETIC_Hx2y2z_G3yz_vrr = PBZ*I_KINETIC_Hx2y2z_F3y_vrr+2*oned2z*I_KINETIC_Gx2yz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G3yz_vrr;
    Double I_KINETIC_Hxy3z_G3yz_vrr = PBZ*I_KINETIC_Hxy3z_F3y_vrr+3*oned2z*I_KINETIC_Gxy2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G3yz_vrr;
    Double I_KINETIC_Hx4z_G3yz_vrr = PBZ*I_KINETIC_Hx4z_F3y_vrr+4*oned2z*I_KINETIC_Gx3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G3yz_vrr;
    Double I_KINETIC_H5y_G3yz_vrr = PBZ*I_KINETIC_H5y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G3yz_vrr;
    Double I_KINETIC_H4yz_G3yz_vrr = PBZ*I_KINETIC_H4yz_F3y_vrr+oned2z*I_KINETIC_G4y_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G3yz_vrr;
    Double I_KINETIC_H3y2z_G3yz_vrr = PBZ*I_KINETIC_H3y2z_F3y_vrr+2*oned2z*I_KINETIC_G3yz_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G3yz_vrr;
    Double I_KINETIC_H2y3z_G3yz_vrr = PBZ*I_KINETIC_H2y3z_F3y_vrr+3*oned2z*I_KINETIC_G2y2z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G3yz_vrr;
    Double I_KINETIC_Hy4z_G3yz_vrr = PBZ*I_KINETIC_Hy4z_F3y_vrr+4*oned2z*I_KINETIC_Gy3z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G3yz_vrr;
    Double I_KINETIC_H5z_G3yz_vrr = PBZ*I_KINETIC_H5z_F3y_vrr+5*oned2z*I_KINETIC_G4z_F3y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G3yz_vrr;
    Double I_KINETIC_H5x_G2y2z_vrr = PBZ*I_KINETIC_H5x_F2yz_vrr+oned2z*I_KINETIC_H5x_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H5x_D2y_vrr;
    Double I_KINETIC_H4xy_G2y2z_vrr = PBZ*I_KINETIC_H4xy_F2yz_vrr+oned2z*I_KINETIC_H4xy_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H4xy_D2y_vrr;
    Double I_KINETIC_H4xz_G2y2z_vrr = PBZ*I_KINETIC_H4xz_F2yz_vrr+oned2z*I_KINETIC_G4x_F2yz_vrr+oned2z*I_KINETIC_H4xz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H4xz_D2y_vrr;
    Double I_KINETIC_H3x2y_G2y2z_vrr = PBZ*I_KINETIC_H3x2y_F2yz_vrr+oned2z*I_KINETIC_H3x2y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H3x2y_D2y_vrr;
    Double I_KINETIC_H3xyz_G2y2z_vrr = PBZ*I_KINETIC_H3xyz_F2yz_vrr+oned2z*I_KINETIC_G3xy_F2yz_vrr+oned2z*I_KINETIC_H3xyz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H3xyz_D2y_vrr;
    Double I_KINETIC_H3x2z_G2y2z_vrr = PBZ*I_KINETIC_H3x2z_F2yz_vrr+2*oned2z*I_KINETIC_G3xz_F2yz_vrr+oned2z*I_KINETIC_H3x2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H3x2z_D2y_vrr;
    Double I_KINETIC_H2x3y_G2y2z_vrr = PBZ*I_KINETIC_H2x3y_F2yz_vrr+oned2z*I_KINETIC_H2x3y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H2x3y_D2y_vrr;
    Double I_KINETIC_H2x2yz_G2y2z_vrr = PBZ*I_KINETIC_H2x2yz_F2yz_vrr+oned2z*I_KINETIC_G2x2y_F2yz_vrr+oned2z*I_KINETIC_H2x2yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H2x2yz_D2y_vrr;
    Double I_KINETIC_H2xy2z_G2y2z_vrr = PBZ*I_KINETIC_H2xy2z_F2yz_vrr+2*oned2z*I_KINETIC_G2xyz_F2yz_vrr+oned2z*I_KINETIC_H2xy2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H2xy2z_D2y_vrr;
    Double I_KINETIC_H2x3z_G2y2z_vrr = PBZ*I_KINETIC_H2x3z_F2yz_vrr+3*oned2z*I_KINETIC_G2x2z_F2yz_vrr+oned2z*I_KINETIC_H2x3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H2x3z_D2y_vrr;
    Double I_KINETIC_Hx4y_G2y2z_vrr = PBZ*I_KINETIC_Hx4y_F2yz_vrr+oned2z*I_KINETIC_Hx4y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Hx4y_D2y_vrr;
    Double I_KINETIC_Hx3yz_G2y2z_vrr = PBZ*I_KINETIC_Hx3yz_F2yz_vrr+oned2z*I_KINETIC_Gx3y_F2yz_vrr+oned2z*I_KINETIC_Hx3yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Hx3yz_D2y_vrr;
    Double I_KINETIC_Hx2y2z_G2y2z_vrr = PBZ*I_KINETIC_Hx2y2z_F2yz_vrr+2*oned2z*I_KINETIC_Gx2yz_F2yz_vrr+oned2z*I_KINETIC_Hx2y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Hx2y2z_D2y_vrr;
    Double I_KINETIC_Hxy3z_G2y2z_vrr = PBZ*I_KINETIC_Hxy3z_F2yz_vrr+3*oned2z*I_KINETIC_Gxy2z_F2yz_vrr+oned2z*I_KINETIC_Hxy3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Hxy3z_D2y_vrr;
    Double I_KINETIC_Hx4z_G2y2z_vrr = PBZ*I_KINETIC_Hx4z_F2yz_vrr+4*oned2z*I_KINETIC_Gx3z_F2yz_vrr+oned2z*I_KINETIC_Hx4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Hx4z_D2y_vrr;
    Double I_KINETIC_H5y_G2y2z_vrr = PBZ*I_KINETIC_H5y_F2yz_vrr+oned2z*I_KINETIC_H5y_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H5y_D2y_vrr;
    Double I_KINETIC_H4yz_G2y2z_vrr = PBZ*I_KINETIC_H4yz_F2yz_vrr+oned2z*I_KINETIC_G4y_F2yz_vrr+oned2z*I_KINETIC_H4yz_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H4yz_D2y_vrr;
    Double I_KINETIC_H3y2z_G2y2z_vrr = PBZ*I_KINETIC_H3y2z_F2yz_vrr+2*oned2z*I_KINETIC_G3yz_F2yz_vrr+oned2z*I_KINETIC_H3y2z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H3y2z_D2y_vrr;
    Double I_KINETIC_H2y3z_G2y2z_vrr = PBZ*I_KINETIC_H2y3z_F2yz_vrr+3*oned2z*I_KINETIC_G2y2z_F2yz_vrr+oned2z*I_KINETIC_H2y3z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H2y3z_D2y_vrr;
    Double I_KINETIC_Hy4z_G2y2z_vrr = PBZ*I_KINETIC_Hy4z_F2yz_vrr+4*oned2z*I_KINETIC_Gy3z_F2yz_vrr+oned2z*I_KINETIC_Hy4z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_Hy4z_D2y_vrr;
    Double I_KINETIC_H5z_G2y2z_vrr = PBZ*I_KINETIC_H5z_F2yz_vrr+5*oned2z*I_KINETIC_G4z_F2yz_vrr+oned2z*I_KINETIC_H5z_D2y_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G2y2z_vrr-adz*I_TWOBODYOVERLAP_H5z_D2y_vrr;
    Double I_KINETIC_H5x_Gy3z_vrr = PBY*I_KINETIC_H5x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Gy3z_vrr;
    Double I_KINETIC_H4xy_Gy3z_vrr = PBY*I_KINETIC_H4xy_F3z_vrr+oned2z*I_KINETIC_G4x_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Gy3z_vrr;
    Double I_KINETIC_H4xz_Gy3z_vrr = PBY*I_KINETIC_H4xz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Gy3z_vrr;
    Double I_KINETIC_H3x2y_Gy3z_vrr = PBY*I_KINETIC_H3x2y_F3z_vrr+2*oned2z*I_KINETIC_G3xy_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Gy3z_vrr;
    Double I_KINETIC_H3xyz_Gy3z_vrr = PBY*I_KINETIC_H3xyz_F3z_vrr+oned2z*I_KINETIC_G3xz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Gy3z_vrr;
    Double I_KINETIC_H3x2z_Gy3z_vrr = PBY*I_KINETIC_H3x2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Gy3z_vrr;
    Double I_KINETIC_H2x3y_Gy3z_vrr = PBY*I_KINETIC_H2x3y_F3z_vrr+3*oned2z*I_KINETIC_G2x2y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Gy3z_vrr;
    Double I_KINETIC_H2x2yz_Gy3z_vrr = PBY*I_KINETIC_H2x2yz_F3z_vrr+2*oned2z*I_KINETIC_G2xyz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Gy3z_vrr;
    Double I_KINETIC_H2xy2z_Gy3z_vrr = PBY*I_KINETIC_H2xy2z_F3z_vrr+oned2z*I_KINETIC_G2x2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Gy3z_vrr;
    Double I_KINETIC_H2x3z_Gy3z_vrr = PBY*I_KINETIC_H2x3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Gy3z_vrr;
    Double I_KINETIC_Hx4y_Gy3z_vrr = PBY*I_KINETIC_Hx4y_F3z_vrr+4*oned2z*I_KINETIC_Gx3y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Gy3z_vrr;
    Double I_KINETIC_Hx3yz_Gy3z_vrr = PBY*I_KINETIC_Hx3yz_F3z_vrr+3*oned2z*I_KINETIC_Gx2yz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Gy3z_vrr;
    Double I_KINETIC_Hx2y2z_Gy3z_vrr = PBY*I_KINETIC_Hx2y2z_F3z_vrr+2*oned2z*I_KINETIC_Gxy2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Gy3z_vrr;
    Double I_KINETIC_Hxy3z_Gy3z_vrr = PBY*I_KINETIC_Hxy3z_F3z_vrr+oned2z*I_KINETIC_Gx3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Gy3z_vrr;
    Double I_KINETIC_Hx4z_Gy3z_vrr = PBY*I_KINETIC_Hx4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Gy3z_vrr;
    Double I_KINETIC_H5y_Gy3z_vrr = PBY*I_KINETIC_H5y_F3z_vrr+5*oned2z*I_KINETIC_G4y_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Gy3z_vrr;
    Double I_KINETIC_H4yz_Gy3z_vrr = PBY*I_KINETIC_H4yz_F3z_vrr+4*oned2z*I_KINETIC_G3yz_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Gy3z_vrr;
    Double I_KINETIC_H3y2z_Gy3z_vrr = PBY*I_KINETIC_H3y2z_F3z_vrr+3*oned2z*I_KINETIC_G2y2z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Gy3z_vrr;
    Double I_KINETIC_H2y3z_Gy3z_vrr = PBY*I_KINETIC_H2y3z_F3z_vrr+2*oned2z*I_KINETIC_Gy3z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Gy3z_vrr;
    Double I_KINETIC_Hy4z_Gy3z_vrr = PBY*I_KINETIC_Hy4z_F3z_vrr+oned2z*I_KINETIC_G4z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Gy3z_vrr;
    Double I_KINETIC_H5z_Gy3z_vrr = PBY*I_KINETIC_H5z_F3z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Gy3z_vrr;
    Double I_KINETIC_H5x_G4z_vrr = PBZ*I_KINETIC_H5x_F3z_vrr+3*oned2z*I_KINETIC_H5x_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5x_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H5x_D2z_vrr;
    Double I_KINETIC_H4xy_G4z_vrr = PBZ*I_KINETIC_H4xy_F3z_vrr+3*oned2z*I_KINETIC_H4xy_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H4xy_D2z_vrr;
    Double I_KINETIC_H4xz_G4z_vrr = PBZ*I_KINETIC_H4xz_F3z_vrr+oned2z*I_KINETIC_G4x_F3z_vrr+3*oned2z*I_KINETIC_H4xz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H4xz_D2z_vrr;
    Double I_KINETIC_H3x2y_G4z_vrr = PBZ*I_KINETIC_H3x2y_F3z_vrr+3*oned2z*I_KINETIC_H3x2y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H3x2y_D2z_vrr;
    Double I_KINETIC_H3xyz_G4z_vrr = PBZ*I_KINETIC_H3xyz_F3z_vrr+oned2z*I_KINETIC_G3xy_F3z_vrr+3*oned2z*I_KINETIC_H3xyz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H3xyz_D2z_vrr;
    Double I_KINETIC_H3x2z_G4z_vrr = PBZ*I_KINETIC_H3x2z_F3z_vrr+2*oned2z*I_KINETIC_G3xz_F3z_vrr+3*oned2z*I_KINETIC_H3x2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H3x2z_D2z_vrr;
    Double I_KINETIC_H2x3y_G4z_vrr = PBZ*I_KINETIC_H2x3y_F3z_vrr+3*oned2z*I_KINETIC_H2x3y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H2x3y_D2z_vrr;
    Double I_KINETIC_H2x2yz_G4z_vrr = PBZ*I_KINETIC_H2x2yz_F3z_vrr+oned2z*I_KINETIC_G2x2y_F3z_vrr+3*oned2z*I_KINETIC_H2x2yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H2x2yz_D2z_vrr;
    Double I_KINETIC_H2xy2z_G4z_vrr = PBZ*I_KINETIC_H2xy2z_F3z_vrr+2*oned2z*I_KINETIC_G2xyz_F3z_vrr+3*oned2z*I_KINETIC_H2xy2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H2xy2z_D2z_vrr;
    Double I_KINETIC_H2x3z_G4z_vrr = PBZ*I_KINETIC_H2x3z_F3z_vrr+3*oned2z*I_KINETIC_G2x2z_F3z_vrr+3*oned2z*I_KINETIC_H2x3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H2x3z_D2z_vrr;
    Double I_KINETIC_Hx4y_G4z_vrr = PBZ*I_KINETIC_Hx4y_F3z_vrr+3*oned2z*I_KINETIC_Hx4y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Hx4y_D2z_vrr;
    Double I_KINETIC_Hx3yz_G4z_vrr = PBZ*I_KINETIC_Hx3yz_F3z_vrr+oned2z*I_KINETIC_Gx3y_F3z_vrr+3*oned2z*I_KINETIC_Hx3yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Hx3yz_D2z_vrr;
    Double I_KINETIC_Hx2y2z_G4z_vrr = PBZ*I_KINETIC_Hx2y2z_F3z_vrr+2*oned2z*I_KINETIC_Gx2yz_F3z_vrr+3*oned2z*I_KINETIC_Hx2y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Hx2y2z_D2z_vrr;
    Double I_KINETIC_Hxy3z_G4z_vrr = PBZ*I_KINETIC_Hxy3z_F3z_vrr+3*oned2z*I_KINETIC_Gxy2z_F3z_vrr+3*oned2z*I_KINETIC_Hxy3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Hxy3z_D2z_vrr;
    Double I_KINETIC_Hx4z_G4z_vrr = PBZ*I_KINETIC_Hx4z_F3z_vrr+4*oned2z*I_KINETIC_Gx3z_F3z_vrr+3*oned2z*I_KINETIC_Hx4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Hx4z_D2z_vrr;
    Double I_KINETIC_H5y_G4z_vrr = PBZ*I_KINETIC_H5y_F3z_vrr+3*oned2z*I_KINETIC_H5y_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5y_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H5y_D2z_vrr;
    Double I_KINETIC_H4yz_G4z_vrr = PBZ*I_KINETIC_H4yz_F3z_vrr+oned2z*I_KINETIC_G4y_F3z_vrr+3*oned2z*I_KINETIC_H4yz_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H4yz_D2z_vrr;
    Double I_KINETIC_H3y2z_G4z_vrr = PBZ*I_KINETIC_H3y2z_F3z_vrr+2*oned2z*I_KINETIC_G3yz_F3z_vrr+3*oned2z*I_KINETIC_H3y2z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H3y2z_D2z_vrr;
    Double I_KINETIC_H2y3z_G4z_vrr = PBZ*I_KINETIC_H2y3z_F3z_vrr+3*oned2z*I_KINETIC_G2y2z_F3z_vrr+3*oned2z*I_KINETIC_H2y3z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H2y3z_D2z_vrr;
    Double I_KINETIC_Hy4z_G4z_vrr = PBZ*I_KINETIC_Hy4z_F3z_vrr+4*oned2z*I_KINETIC_Gy3z_F3z_vrr+3*oned2z*I_KINETIC_Hy4z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_Hy4z_D2z_vrr;
    Double I_KINETIC_H5z_G4z_vrr = PBZ*I_KINETIC_H5z_F3z_vrr+5*oned2z*I_KINETIC_G4z_F3z_vrr+3*oned2z*I_KINETIC_H5z_D2z_vrr+twoxi*I_TWOBODYOVERLAP_H5z_G4z_vrr-3*adz*I_TWOBODYOVERLAP_H5z_D2z_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_G
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    abcd[0] += I_KINETIC_H5x_G4x_vrr;
    abcd[1] += I_KINETIC_H4xy_G4x_vrr;
    abcd[2] += I_KINETIC_H4xz_G4x_vrr;
    abcd[3] += I_KINETIC_H3x2y_G4x_vrr;
    abcd[4] += I_KINETIC_H3xyz_G4x_vrr;
    abcd[5] += I_KINETIC_H3x2z_G4x_vrr;
    abcd[6] += I_KINETIC_H2x3y_G4x_vrr;
    abcd[7] += I_KINETIC_H2x2yz_G4x_vrr;
    abcd[8] += I_KINETIC_H2xy2z_G4x_vrr;
    abcd[9] += I_KINETIC_H2x3z_G4x_vrr;
    abcd[10] += I_KINETIC_Hx4y_G4x_vrr;
    abcd[11] += I_KINETIC_Hx3yz_G4x_vrr;
    abcd[12] += I_KINETIC_Hx2y2z_G4x_vrr;
    abcd[13] += I_KINETIC_Hxy3z_G4x_vrr;
    abcd[14] += I_KINETIC_Hx4z_G4x_vrr;
    abcd[15] += I_KINETIC_H5y_G4x_vrr;
    abcd[16] += I_KINETIC_H4yz_G4x_vrr;
    abcd[17] += I_KINETIC_H3y2z_G4x_vrr;
    abcd[18] += I_KINETIC_H2y3z_G4x_vrr;
    abcd[19] += I_KINETIC_Hy4z_G4x_vrr;
    abcd[20] += I_KINETIC_H5z_G4x_vrr;
    abcd[21] += I_KINETIC_H5x_G3xy_vrr;
    abcd[22] += I_KINETIC_H4xy_G3xy_vrr;
    abcd[23] += I_KINETIC_H4xz_G3xy_vrr;
    abcd[24] += I_KINETIC_H3x2y_G3xy_vrr;
    abcd[25] += I_KINETIC_H3xyz_G3xy_vrr;
    abcd[26] += I_KINETIC_H3x2z_G3xy_vrr;
    abcd[27] += I_KINETIC_H2x3y_G3xy_vrr;
    abcd[28] += I_KINETIC_H2x2yz_G3xy_vrr;
    abcd[29] += I_KINETIC_H2xy2z_G3xy_vrr;
    abcd[30] += I_KINETIC_H2x3z_G3xy_vrr;
    abcd[31] += I_KINETIC_Hx4y_G3xy_vrr;
    abcd[32] += I_KINETIC_Hx3yz_G3xy_vrr;
    abcd[33] += I_KINETIC_Hx2y2z_G3xy_vrr;
    abcd[34] += I_KINETIC_Hxy3z_G3xy_vrr;
    abcd[35] += I_KINETIC_Hx4z_G3xy_vrr;
    abcd[36] += I_KINETIC_H5y_G3xy_vrr;
    abcd[37] += I_KINETIC_H4yz_G3xy_vrr;
    abcd[38] += I_KINETIC_H3y2z_G3xy_vrr;
    abcd[39] += I_KINETIC_H2y3z_G3xy_vrr;
    abcd[40] += I_KINETIC_Hy4z_G3xy_vrr;
    abcd[41] += I_KINETIC_H5z_G3xy_vrr;
    abcd[42] += I_KINETIC_H5x_G3xz_vrr;
    abcd[43] += I_KINETIC_H4xy_G3xz_vrr;
    abcd[44] += I_KINETIC_H4xz_G3xz_vrr;
    abcd[45] += I_KINETIC_H3x2y_G3xz_vrr;
    abcd[46] += I_KINETIC_H3xyz_G3xz_vrr;
    abcd[47] += I_KINETIC_H3x2z_G3xz_vrr;
    abcd[48] += I_KINETIC_H2x3y_G3xz_vrr;
    abcd[49] += I_KINETIC_H2x2yz_G3xz_vrr;
    abcd[50] += I_KINETIC_H2xy2z_G3xz_vrr;
    abcd[51] += I_KINETIC_H2x3z_G3xz_vrr;
    abcd[52] += I_KINETIC_Hx4y_G3xz_vrr;
    abcd[53] += I_KINETIC_Hx3yz_G3xz_vrr;
    abcd[54] += I_KINETIC_Hx2y2z_G3xz_vrr;
    abcd[55] += I_KINETIC_Hxy3z_G3xz_vrr;
    abcd[56] += I_KINETIC_Hx4z_G3xz_vrr;
    abcd[57] += I_KINETIC_H5y_G3xz_vrr;
    abcd[58] += I_KINETIC_H4yz_G3xz_vrr;
    abcd[59] += I_KINETIC_H3y2z_G3xz_vrr;
    abcd[60] += I_KINETIC_H2y3z_G3xz_vrr;
    abcd[61] += I_KINETIC_Hy4z_G3xz_vrr;
    abcd[62] += I_KINETIC_H5z_G3xz_vrr;
    abcd[63] += I_KINETIC_H5x_G2x2y_vrr;
    abcd[64] += I_KINETIC_H4xy_G2x2y_vrr;
    abcd[65] += I_KINETIC_H4xz_G2x2y_vrr;
    abcd[66] += I_KINETIC_H3x2y_G2x2y_vrr;
    abcd[67] += I_KINETIC_H3xyz_G2x2y_vrr;
    abcd[68] += I_KINETIC_H3x2z_G2x2y_vrr;
    abcd[69] += I_KINETIC_H2x3y_G2x2y_vrr;
    abcd[70] += I_KINETIC_H2x2yz_G2x2y_vrr;
    abcd[71] += I_KINETIC_H2xy2z_G2x2y_vrr;
    abcd[72] += I_KINETIC_H2x3z_G2x2y_vrr;
    abcd[73] += I_KINETIC_Hx4y_G2x2y_vrr;
    abcd[74] += I_KINETIC_Hx3yz_G2x2y_vrr;
    abcd[75] += I_KINETIC_Hx2y2z_G2x2y_vrr;
    abcd[76] += I_KINETIC_Hxy3z_G2x2y_vrr;
    abcd[77] += I_KINETIC_Hx4z_G2x2y_vrr;
    abcd[78] += I_KINETIC_H5y_G2x2y_vrr;
    abcd[79] += I_KINETIC_H4yz_G2x2y_vrr;
    abcd[80] += I_KINETIC_H3y2z_G2x2y_vrr;
    abcd[81] += I_KINETIC_H2y3z_G2x2y_vrr;
    abcd[82] += I_KINETIC_Hy4z_G2x2y_vrr;
    abcd[83] += I_KINETIC_H5z_G2x2y_vrr;
    abcd[84] += I_KINETIC_H5x_G2xyz_vrr;
    abcd[85] += I_KINETIC_H4xy_G2xyz_vrr;
    abcd[86] += I_KINETIC_H4xz_G2xyz_vrr;
    abcd[87] += I_KINETIC_H3x2y_G2xyz_vrr;
    abcd[88] += I_KINETIC_H3xyz_G2xyz_vrr;
    abcd[89] += I_KINETIC_H3x2z_G2xyz_vrr;
    abcd[90] += I_KINETIC_H2x3y_G2xyz_vrr;
    abcd[91] += I_KINETIC_H2x2yz_G2xyz_vrr;
    abcd[92] += I_KINETIC_H2xy2z_G2xyz_vrr;
    abcd[93] += I_KINETIC_H2x3z_G2xyz_vrr;
    abcd[94] += I_KINETIC_Hx4y_G2xyz_vrr;
    abcd[95] += I_KINETIC_Hx3yz_G2xyz_vrr;
    abcd[96] += I_KINETIC_Hx2y2z_G2xyz_vrr;
    abcd[97] += I_KINETIC_Hxy3z_G2xyz_vrr;
    abcd[98] += I_KINETIC_Hx4z_G2xyz_vrr;
    abcd[99] += I_KINETIC_H5y_G2xyz_vrr;
    abcd[100] += I_KINETIC_H4yz_G2xyz_vrr;
    abcd[101] += I_KINETIC_H3y2z_G2xyz_vrr;
    abcd[102] += I_KINETIC_H2y3z_G2xyz_vrr;
    abcd[103] += I_KINETIC_Hy4z_G2xyz_vrr;
    abcd[104] += I_KINETIC_H5z_G2xyz_vrr;
    abcd[105] += I_KINETIC_H5x_G2x2z_vrr;
    abcd[106] += I_KINETIC_H4xy_G2x2z_vrr;
    abcd[107] += I_KINETIC_H4xz_G2x2z_vrr;
    abcd[108] += I_KINETIC_H3x2y_G2x2z_vrr;
    abcd[109] += I_KINETIC_H3xyz_G2x2z_vrr;
    abcd[110] += I_KINETIC_H3x2z_G2x2z_vrr;
    abcd[111] += I_KINETIC_H2x3y_G2x2z_vrr;
    abcd[112] += I_KINETIC_H2x2yz_G2x2z_vrr;
    abcd[113] += I_KINETIC_H2xy2z_G2x2z_vrr;
    abcd[114] += I_KINETIC_H2x3z_G2x2z_vrr;
    abcd[115] += I_KINETIC_Hx4y_G2x2z_vrr;
    abcd[116] += I_KINETIC_Hx3yz_G2x2z_vrr;
    abcd[117] += I_KINETIC_Hx2y2z_G2x2z_vrr;
    abcd[118] += I_KINETIC_Hxy3z_G2x2z_vrr;
    abcd[119] += I_KINETIC_Hx4z_G2x2z_vrr;
    abcd[120] += I_KINETIC_H5y_G2x2z_vrr;
    abcd[121] += I_KINETIC_H4yz_G2x2z_vrr;
    abcd[122] += I_KINETIC_H3y2z_G2x2z_vrr;
    abcd[123] += I_KINETIC_H2y3z_G2x2z_vrr;
    abcd[124] += I_KINETIC_Hy4z_G2x2z_vrr;
    abcd[125] += I_KINETIC_H5z_G2x2z_vrr;
    abcd[126] += I_KINETIC_H5x_Gx3y_vrr;
    abcd[127] += I_KINETIC_H4xy_Gx3y_vrr;
    abcd[128] += I_KINETIC_H4xz_Gx3y_vrr;
    abcd[129] += I_KINETIC_H3x2y_Gx3y_vrr;
    abcd[130] += I_KINETIC_H3xyz_Gx3y_vrr;
    abcd[131] += I_KINETIC_H3x2z_Gx3y_vrr;
    abcd[132] += I_KINETIC_H2x3y_Gx3y_vrr;
    abcd[133] += I_KINETIC_H2x2yz_Gx3y_vrr;
    abcd[134] += I_KINETIC_H2xy2z_Gx3y_vrr;
    abcd[135] += I_KINETIC_H2x3z_Gx3y_vrr;
    abcd[136] += I_KINETIC_Hx4y_Gx3y_vrr;
    abcd[137] += I_KINETIC_Hx3yz_Gx3y_vrr;
    abcd[138] += I_KINETIC_Hx2y2z_Gx3y_vrr;
    abcd[139] += I_KINETIC_Hxy3z_Gx3y_vrr;
    abcd[140] += I_KINETIC_Hx4z_Gx3y_vrr;
    abcd[141] += I_KINETIC_H5y_Gx3y_vrr;
    abcd[142] += I_KINETIC_H4yz_Gx3y_vrr;
    abcd[143] += I_KINETIC_H3y2z_Gx3y_vrr;
    abcd[144] += I_KINETIC_H2y3z_Gx3y_vrr;
    abcd[145] += I_KINETIC_Hy4z_Gx3y_vrr;
    abcd[146] += I_KINETIC_H5z_Gx3y_vrr;
    abcd[147] += I_KINETIC_H5x_Gx2yz_vrr;
    abcd[148] += I_KINETIC_H4xy_Gx2yz_vrr;
    abcd[149] += I_KINETIC_H4xz_Gx2yz_vrr;
    abcd[150] += I_KINETIC_H3x2y_Gx2yz_vrr;
    abcd[151] += I_KINETIC_H3xyz_Gx2yz_vrr;
    abcd[152] += I_KINETIC_H3x2z_Gx2yz_vrr;
    abcd[153] += I_KINETIC_H2x3y_Gx2yz_vrr;
    abcd[154] += I_KINETIC_H2x2yz_Gx2yz_vrr;
    abcd[155] += I_KINETIC_H2xy2z_Gx2yz_vrr;
    abcd[156] += I_KINETIC_H2x3z_Gx2yz_vrr;
    abcd[157] += I_KINETIC_Hx4y_Gx2yz_vrr;
    abcd[158] += I_KINETIC_Hx3yz_Gx2yz_vrr;
    abcd[159] += I_KINETIC_Hx2y2z_Gx2yz_vrr;
    abcd[160] += I_KINETIC_Hxy3z_Gx2yz_vrr;
    abcd[161] += I_KINETIC_Hx4z_Gx2yz_vrr;
    abcd[162] += I_KINETIC_H5y_Gx2yz_vrr;
    abcd[163] += I_KINETIC_H4yz_Gx2yz_vrr;
    abcd[164] += I_KINETIC_H3y2z_Gx2yz_vrr;
    abcd[165] += I_KINETIC_H2y3z_Gx2yz_vrr;
    abcd[166] += I_KINETIC_Hy4z_Gx2yz_vrr;
    abcd[167] += I_KINETIC_H5z_Gx2yz_vrr;
    abcd[168] += I_KINETIC_H5x_Gxy2z_vrr;
    abcd[169] += I_KINETIC_H4xy_Gxy2z_vrr;
    abcd[170] += I_KINETIC_H4xz_Gxy2z_vrr;
    abcd[171] += I_KINETIC_H3x2y_Gxy2z_vrr;
    abcd[172] += I_KINETIC_H3xyz_Gxy2z_vrr;
    abcd[173] += I_KINETIC_H3x2z_Gxy2z_vrr;
    abcd[174] += I_KINETIC_H2x3y_Gxy2z_vrr;
    abcd[175] += I_KINETIC_H2x2yz_Gxy2z_vrr;
    abcd[176] += I_KINETIC_H2xy2z_Gxy2z_vrr;
    abcd[177] += I_KINETIC_H2x3z_Gxy2z_vrr;
    abcd[178] += I_KINETIC_Hx4y_Gxy2z_vrr;
    abcd[179] += I_KINETIC_Hx3yz_Gxy2z_vrr;
    abcd[180] += I_KINETIC_Hx2y2z_Gxy2z_vrr;
    abcd[181] += I_KINETIC_Hxy3z_Gxy2z_vrr;
    abcd[182] += I_KINETIC_Hx4z_Gxy2z_vrr;
    abcd[183] += I_KINETIC_H5y_Gxy2z_vrr;
    abcd[184] += I_KINETIC_H4yz_Gxy2z_vrr;
    abcd[185] += I_KINETIC_H3y2z_Gxy2z_vrr;
    abcd[186] += I_KINETIC_H2y3z_Gxy2z_vrr;
    abcd[187] += I_KINETIC_Hy4z_Gxy2z_vrr;
    abcd[188] += I_KINETIC_H5z_Gxy2z_vrr;
    abcd[189] += I_KINETIC_H5x_Gx3z_vrr;
    abcd[190] += I_KINETIC_H4xy_Gx3z_vrr;
    abcd[191] += I_KINETIC_H4xz_Gx3z_vrr;
    abcd[192] += I_KINETIC_H3x2y_Gx3z_vrr;
    abcd[193] += I_KINETIC_H3xyz_Gx3z_vrr;
    abcd[194] += I_KINETIC_H3x2z_Gx3z_vrr;
    abcd[195] += I_KINETIC_H2x3y_Gx3z_vrr;
    abcd[196] += I_KINETIC_H2x2yz_Gx3z_vrr;
    abcd[197] += I_KINETIC_H2xy2z_Gx3z_vrr;
    abcd[198] += I_KINETIC_H2x3z_Gx3z_vrr;
    abcd[199] += I_KINETIC_Hx4y_Gx3z_vrr;
    abcd[200] += I_KINETIC_Hx3yz_Gx3z_vrr;
    abcd[201] += I_KINETIC_Hx2y2z_Gx3z_vrr;
    abcd[202] += I_KINETIC_Hxy3z_Gx3z_vrr;
    abcd[203] += I_KINETIC_Hx4z_Gx3z_vrr;
    abcd[204] += I_KINETIC_H5y_Gx3z_vrr;
    abcd[205] += I_KINETIC_H4yz_Gx3z_vrr;
    abcd[206] += I_KINETIC_H3y2z_Gx3z_vrr;
    abcd[207] += I_KINETIC_H2y3z_Gx3z_vrr;
    abcd[208] += I_KINETIC_Hy4z_Gx3z_vrr;
    abcd[209] += I_KINETIC_H5z_Gx3z_vrr;
    abcd[210] += I_KINETIC_H5x_G4y_vrr;
    abcd[211] += I_KINETIC_H4xy_G4y_vrr;
    abcd[212] += I_KINETIC_H4xz_G4y_vrr;
    abcd[213] += I_KINETIC_H3x2y_G4y_vrr;
    abcd[214] += I_KINETIC_H3xyz_G4y_vrr;
    abcd[215] += I_KINETIC_H3x2z_G4y_vrr;
    abcd[216] += I_KINETIC_H2x3y_G4y_vrr;
    abcd[217] += I_KINETIC_H2x2yz_G4y_vrr;
    abcd[218] += I_KINETIC_H2xy2z_G4y_vrr;
    abcd[219] += I_KINETIC_H2x3z_G4y_vrr;
    abcd[220] += I_KINETIC_Hx4y_G4y_vrr;
    abcd[221] += I_KINETIC_Hx3yz_G4y_vrr;
    abcd[222] += I_KINETIC_Hx2y2z_G4y_vrr;
    abcd[223] += I_KINETIC_Hxy3z_G4y_vrr;
    abcd[224] += I_KINETIC_Hx4z_G4y_vrr;
    abcd[225] += I_KINETIC_H5y_G4y_vrr;
    abcd[226] += I_KINETIC_H4yz_G4y_vrr;
    abcd[227] += I_KINETIC_H3y2z_G4y_vrr;
    abcd[228] += I_KINETIC_H2y3z_G4y_vrr;
    abcd[229] += I_KINETIC_Hy4z_G4y_vrr;
    abcd[230] += I_KINETIC_H5z_G4y_vrr;
    abcd[231] += I_KINETIC_H5x_G3yz_vrr;
    abcd[232] += I_KINETIC_H4xy_G3yz_vrr;
    abcd[233] += I_KINETIC_H4xz_G3yz_vrr;
    abcd[234] += I_KINETIC_H3x2y_G3yz_vrr;
    abcd[235] += I_KINETIC_H3xyz_G3yz_vrr;
    abcd[236] += I_KINETIC_H3x2z_G3yz_vrr;
    abcd[237] += I_KINETIC_H2x3y_G3yz_vrr;
    abcd[238] += I_KINETIC_H2x2yz_G3yz_vrr;
    abcd[239] += I_KINETIC_H2xy2z_G3yz_vrr;
    abcd[240] += I_KINETIC_H2x3z_G3yz_vrr;
    abcd[241] += I_KINETIC_Hx4y_G3yz_vrr;
    abcd[242] += I_KINETIC_Hx3yz_G3yz_vrr;
    abcd[243] += I_KINETIC_Hx2y2z_G3yz_vrr;
    abcd[244] += I_KINETIC_Hxy3z_G3yz_vrr;
    abcd[245] += I_KINETIC_Hx4z_G3yz_vrr;
    abcd[246] += I_KINETIC_H5y_G3yz_vrr;
    abcd[247] += I_KINETIC_H4yz_G3yz_vrr;
    abcd[248] += I_KINETIC_H3y2z_G3yz_vrr;
    abcd[249] += I_KINETIC_H2y3z_G3yz_vrr;
    abcd[250] += I_KINETIC_Hy4z_G3yz_vrr;
    abcd[251] += I_KINETIC_H5z_G3yz_vrr;
    abcd[252] += I_KINETIC_H5x_G2y2z_vrr;
    abcd[253] += I_KINETIC_H4xy_G2y2z_vrr;
    abcd[254] += I_KINETIC_H4xz_G2y2z_vrr;
    abcd[255] += I_KINETIC_H3x2y_G2y2z_vrr;
    abcd[256] += I_KINETIC_H3xyz_G2y2z_vrr;
    abcd[257] += I_KINETIC_H3x2z_G2y2z_vrr;
    abcd[258] += I_KINETIC_H2x3y_G2y2z_vrr;
    abcd[259] += I_KINETIC_H2x2yz_G2y2z_vrr;
    abcd[260] += I_KINETIC_H2xy2z_G2y2z_vrr;
    abcd[261] += I_KINETIC_H2x3z_G2y2z_vrr;
    abcd[262] += I_KINETIC_Hx4y_G2y2z_vrr;
    abcd[263] += I_KINETIC_Hx3yz_G2y2z_vrr;
    abcd[264] += I_KINETIC_Hx2y2z_G2y2z_vrr;
    abcd[265] += I_KINETIC_Hxy3z_G2y2z_vrr;
    abcd[266] += I_KINETIC_Hx4z_G2y2z_vrr;
    abcd[267] += I_KINETIC_H5y_G2y2z_vrr;
    abcd[268] += I_KINETIC_H4yz_G2y2z_vrr;
    abcd[269] += I_KINETIC_H3y2z_G2y2z_vrr;
    abcd[270] += I_KINETIC_H2y3z_G2y2z_vrr;
    abcd[271] += I_KINETIC_Hy4z_G2y2z_vrr;
    abcd[272] += I_KINETIC_H5z_G2y2z_vrr;
    abcd[273] += I_KINETIC_H5x_Gy3z_vrr;
    abcd[274] += I_KINETIC_H4xy_Gy3z_vrr;
    abcd[275] += I_KINETIC_H4xz_Gy3z_vrr;
    abcd[276] += I_KINETIC_H3x2y_Gy3z_vrr;
    abcd[277] += I_KINETIC_H3xyz_Gy3z_vrr;
    abcd[278] += I_KINETIC_H3x2z_Gy3z_vrr;
    abcd[279] += I_KINETIC_H2x3y_Gy3z_vrr;
    abcd[280] += I_KINETIC_H2x2yz_Gy3z_vrr;
    abcd[281] += I_KINETIC_H2xy2z_Gy3z_vrr;
    abcd[282] += I_KINETIC_H2x3z_Gy3z_vrr;
    abcd[283] += I_KINETIC_Hx4y_Gy3z_vrr;
    abcd[284] += I_KINETIC_Hx3yz_Gy3z_vrr;
    abcd[285] += I_KINETIC_Hx2y2z_Gy3z_vrr;
    abcd[286] += I_KINETIC_Hxy3z_Gy3z_vrr;
    abcd[287] += I_KINETIC_Hx4z_Gy3z_vrr;
    abcd[288] += I_KINETIC_H5y_Gy3z_vrr;
    abcd[289] += I_KINETIC_H4yz_Gy3z_vrr;
    abcd[290] += I_KINETIC_H3y2z_Gy3z_vrr;
    abcd[291] += I_KINETIC_H2y3z_Gy3z_vrr;
    abcd[292] += I_KINETIC_Hy4z_Gy3z_vrr;
    abcd[293] += I_KINETIC_H5z_Gy3z_vrr;
    abcd[294] += I_KINETIC_H5x_G4z_vrr;
    abcd[295] += I_KINETIC_H4xy_G4z_vrr;
    abcd[296] += I_KINETIC_H4xz_G4z_vrr;
    abcd[297] += I_KINETIC_H3x2y_G4z_vrr;
    abcd[298] += I_KINETIC_H3xyz_G4z_vrr;
    abcd[299] += I_KINETIC_H3x2z_G4z_vrr;
    abcd[300] += I_KINETIC_H2x3y_G4z_vrr;
    abcd[301] += I_KINETIC_H2x2yz_G4z_vrr;
    abcd[302] += I_KINETIC_H2xy2z_G4z_vrr;
    abcd[303] += I_KINETIC_H2x3z_G4z_vrr;
    abcd[304] += I_KINETIC_Hx4y_G4z_vrr;
    abcd[305] += I_KINETIC_Hx3yz_G4z_vrr;
    abcd[306] += I_KINETIC_Hx2y2z_G4z_vrr;
    abcd[307] += I_KINETIC_Hxy3z_G4z_vrr;
    abcd[308] += I_KINETIC_Hx4z_G4z_vrr;
    abcd[309] += I_KINETIC_H5y_G4z_vrr;
    abcd[310] += I_KINETIC_H4yz_G4z_vrr;
    abcd[311] += I_KINETIC_H3y2z_G4z_vrr;
    abcd[312] += I_KINETIC_H2y3z_G4z_vrr;
    abcd[313] += I_KINETIC_Hy4z_G4z_vrr;
    abcd[314] += I_KINETIC_H5z_G4z_vrr;
  }
}
