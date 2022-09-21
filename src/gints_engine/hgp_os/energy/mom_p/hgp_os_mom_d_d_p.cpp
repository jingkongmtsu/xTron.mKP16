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

void hgp_os_mom_d_d_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_MOM_G4x_S_Px = 0.0E0;
  Double I_MOM_G3xy_S_Px = 0.0E0;
  Double I_MOM_G3xz_S_Px = 0.0E0;
  Double I_MOM_G2x2y_S_Px = 0.0E0;
  Double I_MOM_G2xyz_S_Px = 0.0E0;
  Double I_MOM_G2x2z_S_Px = 0.0E0;
  Double I_MOM_Gx3y_S_Px = 0.0E0;
  Double I_MOM_Gx2yz_S_Px = 0.0E0;
  Double I_MOM_Gxy2z_S_Px = 0.0E0;
  Double I_MOM_Gx3z_S_Px = 0.0E0;
  Double I_MOM_G4y_S_Px = 0.0E0;
  Double I_MOM_G3yz_S_Px = 0.0E0;
  Double I_MOM_G2y2z_S_Px = 0.0E0;
  Double I_MOM_Gy3z_S_Px = 0.0E0;
  Double I_MOM_G4z_S_Px = 0.0E0;
  Double I_MOM_G4x_S_Py = 0.0E0;
  Double I_MOM_G3xy_S_Py = 0.0E0;
  Double I_MOM_G3xz_S_Py = 0.0E0;
  Double I_MOM_G2x2y_S_Py = 0.0E0;
  Double I_MOM_G2xyz_S_Py = 0.0E0;
  Double I_MOM_G2x2z_S_Py = 0.0E0;
  Double I_MOM_Gx3y_S_Py = 0.0E0;
  Double I_MOM_Gx2yz_S_Py = 0.0E0;
  Double I_MOM_Gxy2z_S_Py = 0.0E0;
  Double I_MOM_Gx3z_S_Py = 0.0E0;
  Double I_MOM_G4y_S_Py = 0.0E0;
  Double I_MOM_G3yz_S_Py = 0.0E0;
  Double I_MOM_G2y2z_S_Py = 0.0E0;
  Double I_MOM_Gy3z_S_Py = 0.0E0;
  Double I_MOM_G4z_S_Py = 0.0E0;
  Double I_MOM_G4x_S_Pz = 0.0E0;
  Double I_MOM_G3xy_S_Pz = 0.0E0;
  Double I_MOM_G3xz_S_Pz = 0.0E0;
  Double I_MOM_G2x2y_S_Pz = 0.0E0;
  Double I_MOM_G2xyz_S_Pz = 0.0E0;
  Double I_MOM_G2x2z_S_Pz = 0.0E0;
  Double I_MOM_Gx3y_S_Pz = 0.0E0;
  Double I_MOM_Gx2yz_S_Pz = 0.0E0;
  Double I_MOM_Gxy2z_S_Pz = 0.0E0;
  Double I_MOM_Gx3z_S_Pz = 0.0E0;
  Double I_MOM_G4y_S_Pz = 0.0E0;
  Double I_MOM_G3yz_S_Pz = 0.0E0;
  Double I_MOM_G2y2z_S_Pz = 0.0E0;
  Double I_MOM_Gy3z_S_Pz = 0.0E0;
  Double I_MOM_G4z_S_Pz = 0.0E0;
  Double I_MOM_F3x_S_Px = 0.0E0;
  Double I_MOM_F2xy_S_Px = 0.0E0;
  Double I_MOM_F2xz_S_Px = 0.0E0;
  Double I_MOM_Fx2y_S_Px = 0.0E0;
  Double I_MOM_Fxyz_S_Px = 0.0E0;
  Double I_MOM_Fx2z_S_Px = 0.0E0;
  Double I_MOM_F3y_S_Px = 0.0E0;
  Double I_MOM_F2yz_S_Px = 0.0E0;
  Double I_MOM_Fy2z_S_Px = 0.0E0;
  Double I_MOM_F3z_S_Px = 0.0E0;
  Double I_MOM_F3x_S_Py = 0.0E0;
  Double I_MOM_F2xy_S_Py = 0.0E0;
  Double I_MOM_F2xz_S_Py = 0.0E0;
  Double I_MOM_Fx2y_S_Py = 0.0E0;
  Double I_MOM_Fxyz_S_Py = 0.0E0;
  Double I_MOM_Fx2z_S_Py = 0.0E0;
  Double I_MOM_F3y_S_Py = 0.0E0;
  Double I_MOM_F2yz_S_Py = 0.0E0;
  Double I_MOM_Fy2z_S_Py = 0.0E0;
  Double I_MOM_F3z_S_Py = 0.0E0;
  Double I_MOM_F3x_S_Pz = 0.0E0;
  Double I_MOM_F2xy_S_Pz = 0.0E0;
  Double I_MOM_F2xz_S_Pz = 0.0E0;
  Double I_MOM_Fx2y_S_Pz = 0.0E0;
  Double I_MOM_Fxyz_S_Pz = 0.0E0;
  Double I_MOM_Fx2z_S_Pz = 0.0E0;
  Double I_MOM_F3y_S_Pz = 0.0E0;
  Double I_MOM_F2yz_S_Pz = 0.0E0;
  Double I_MOM_Fy2z_S_Pz = 0.0E0;
  Double I_MOM_F3z_S_Pz = 0.0E0;
  Double I_MOM_D2x_S_Px = 0.0E0;
  Double I_MOM_Dxy_S_Px = 0.0E0;
  Double I_MOM_Dxz_S_Px = 0.0E0;
  Double I_MOM_D2y_S_Px = 0.0E0;
  Double I_MOM_Dyz_S_Px = 0.0E0;
  Double I_MOM_D2z_S_Px = 0.0E0;
  Double I_MOM_D2x_S_Py = 0.0E0;
  Double I_MOM_Dxy_S_Py = 0.0E0;
  Double I_MOM_Dxz_S_Py = 0.0E0;
  Double I_MOM_D2y_S_Py = 0.0E0;
  Double I_MOM_Dyz_S_Py = 0.0E0;
  Double I_MOM_D2z_S_Py = 0.0E0;
  Double I_MOM_D2x_S_Pz = 0.0E0;
  Double I_MOM_Dxy_S_Pz = 0.0E0;
  Double I_MOM_Dxz_S_Pz = 0.0E0;
  Double I_MOM_D2y_S_Pz = 0.0E0;
  Double I_MOM_Dyz_S_Pz = 0.0E0;
  Double I_MOM_D2z_S_Pz = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
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
    Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;
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
     * shell quartet name: SQ_MOM_D_S_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_MOM_P_S_S
     * RHS shell quartet name: SQ_MOM_S_S_S
     ************************************************************/
    Double I_MOM_D2x_S_S_vrr = PAX*I_MOM_Px_S_S_vrr+oned2z*I_MOM_S_S_S_vrr;
    Double I_MOM_Dxy_S_S_vrr = PAY*I_MOM_Px_S_S_vrr;
    Double I_MOM_D2y_S_S_vrr = PAY*I_MOM_Py_S_S_vrr+oned2z*I_MOM_S_S_S_vrr;
    Double I_MOM_D2z_S_S_vrr = PAZ*I_MOM_Pz_S_S_vrr+oned2z*I_MOM_S_S_S_vrr;

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
     * shell quartet name: SQ_MOM_F_S_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_MOM_D_S_S
     * RHS shell quartet name: SQ_MOM_P_S_S
     ************************************************************/
    Double I_MOM_F3x_S_S_vrr = PAX*I_MOM_D2x_S_S_vrr+2*oned2z*I_MOM_Px_S_S_vrr;
    Double I_MOM_F2xy_S_S_vrr = PAY*I_MOM_D2x_S_S_vrr;
    Double I_MOM_F2xz_S_S_vrr = PAZ*I_MOM_D2x_S_S_vrr;
    Double I_MOM_Fx2y_S_S_vrr = PAX*I_MOM_D2y_S_S_vrr;
    Double I_MOM_Fx2z_S_S_vrr = PAX*I_MOM_D2z_S_S_vrr;
    Double I_MOM_F3y_S_S_vrr = PAY*I_MOM_D2y_S_S_vrr+2*oned2z*I_MOM_Py_S_S_vrr;
    Double I_MOM_F2yz_S_S_vrr = PAZ*I_MOM_D2y_S_S_vrr;
    Double I_MOM_F3z_S_S_vrr = PAZ*I_MOM_D2z_S_S_vrr+2*oned2z*I_MOM_Pz_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_F_S_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_MOM_D_S_P
     * RHS shell quartet name: SQ_MOM_P_S_P
     * RHS shell quartet name: SQ_MOM_D_S_S
     ************************************************************/
    Double I_MOM_F3x_S_Px_vrr = PAX*I_MOM_D2x_S_Px_vrr+2*oned2z*I_MOM_Px_S_Px_vrr+oned2z*I_MOM_D2x_S_S_vrr;
    Double I_MOM_F2xy_S_Px_vrr = PAY*I_MOM_D2x_S_Px_vrr;
    Double I_MOM_F2xz_S_Px_vrr = PAZ*I_MOM_D2x_S_Px_vrr;
    Double I_MOM_Fx2y_S_Px_vrr = PAX*I_MOM_D2y_S_Px_vrr+oned2z*I_MOM_D2y_S_S_vrr;
    Double I_MOM_Fxyz_S_Px_vrr = PAZ*I_MOM_Dxy_S_Px_vrr;
    Double I_MOM_Fx2z_S_Px_vrr = PAX*I_MOM_D2z_S_Px_vrr+oned2z*I_MOM_D2z_S_S_vrr;
    Double I_MOM_F3y_S_Px_vrr = PAY*I_MOM_D2y_S_Px_vrr+2*oned2z*I_MOM_Py_S_Px_vrr;
    Double I_MOM_F2yz_S_Px_vrr = PAZ*I_MOM_D2y_S_Px_vrr;
    Double I_MOM_Fy2z_S_Px_vrr = PAY*I_MOM_D2z_S_Px_vrr;
    Double I_MOM_F3z_S_Px_vrr = PAZ*I_MOM_D2z_S_Px_vrr+2*oned2z*I_MOM_Pz_S_Px_vrr;
    Double I_MOM_F3x_S_Py_vrr = PAX*I_MOM_D2x_S_Py_vrr+2*oned2z*I_MOM_Px_S_Py_vrr;
    Double I_MOM_F2xy_S_Py_vrr = PAY*I_MOM_D2x_S_Py_vrr+oned2z*I_MOM_D2x_S_S_vrr;
    Double I_MOM_F2xz_S_Py_vrr = PAZ*I_MOM_D2x_S_Py_vrr;
    Double I_MOM_Fx2y_S_Py_vrr = PAX*I_MOM_D2y_S_Py_vrr;
    Double I_MOM_Fxyz_S_Py_vrr = PAZ*I_MOM_Dxy_S_Py_vrr;
    Double I_MOM_Fx2z_S_Py_vrr = PAX*I_MOM_D2z_S_Py_vrr;
    Double I_MOM_F3y_S_Py_vrr = PAY*I_MOM_D2y_S_Py_vrr+2*oned2z*I_MOM_Py_S_Py_vrr+oned2z*I_MOM_D2y_S_S_vrr;
    Double I_MOM_F2yz_S_Py_vrr = PAZ*I_MOM_D2y_S_Py_vrr;
    Double I_MOM_Fy2z_S_Py_vrr = PAY*I_MOM_D2z_S_Py_vrr+oned2z*I_MOM_D2z_S_S_vrr;
    Double I_MOM_F3z_S_Py_vrr = PAZ*I_MOM_D2z_S_Py_vrr+2*oned2z*I_MOM_Pz_S_Py_vrr;
    Double I_MOM_F3x_S_Pz_vrr = PAX*I_MOM_D2x_S_Pz_vrr+2*oned2z*I_MOM_Px_S_Pz_vrr;
    Double I_MOM_F2xy_S_Pz_vrr = PAY*I_MOM_D2x_S_Pz_vrr;
    Double I_MOM_F2xz_S_Pz_vrr = PAZ*I_MOM_D2x_S_Pz_vrr+oned2z*I_MOM_D2x_S_S_vrr;
    Double I_MOM_Fx2y_S_Pz_vrr = PAX*I_MOM_D2y_S_Pz_vrr;
    Double I_MOM_Fxyz_S_Pz_vrr = PAZ*I_MOM_Dxy_S_Pz_vrr+oned2z*I_MOM_Dxy_S_S_vrr;
    Double I_MOM_Fx2z_S_Pz_vrr = PAX*I_MOM_D2z_S_Pz_vrr;
    Double I_MOM_F3y_S_Pz_vrr = PAY*I_MOM_D2y_S_Pz_vrr+2*oned2z*I_MOM_Py_S_Pz_vrr;
    Double I_MOM_F2yz_S_Pz_vrr = PAZ*I_MOM_D2y_S_Pz_vrr+oned2z*I_MOM_D2y_S_S_vrr;
    Double I_MOM_Fy2z_S_Pz_vrr = PAY*I_MOM_D2z_S_Pz_vrr;
    Double I_MOM_F3z_S_Pz_vrr = PAZ*I_MOM_D2z_S_Pz_vrr+2*oned2z*I_MOM_Pz_S_Pz_vrr+oned2z*I_MOM_D2z_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_G_S_P
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_MOM_F_S_P
     * RHS shell quartet name: SQ_MOM_D_S_P
     * RHS shell quartet name: SQ_MOM_F_S_S
     ************************************************************/
    Double I_MOM_G4x_S_Px_vrr = PAX*I_MOM_F3x_S_Px_vrr+3*oned2z*I_MOM_D2x_S_Px_vrr+oned2z*I_MOM_F3x_S_S_vrr;
    Double I_MOM_G3xy_S_Px_vrr = PAY*I_MOM_F3x_S_Px_vrr;
    Double I_MOM_G3xz_S_Px_vrr = PAZ*I_MOM_F3x_S_Px_vrr;
    Double I_MOM_G2x2y_S_Px_vrr = PAY*I_MOM_F2xy_S_Px_vrr+oned2z*I_MOM_D2x_S_Px_vrr;
    Double I_MOM_G2xyz_S_Px_vrr = PAZ*I_MOM_F2xy_S_Px_vrr;
    Double I_MOM_G2x2z_S_Px_vrr = PAZ*I_MOM_F2xz_S_Px_vrr+oned2z*I_MOM_D2x_S_Px_vrr;
    Double I_MOM_Gx3y_S_Px_vrr = PAX*I_MOM_F3y_S_Px_vrr+oned2z*I_MOM_F3y_S_S_vrr;
    Double I_MOM_Gx2yz_S_Px_vrr = PAZ*I_MOM_Fx2y_S_Px_vrr;
    Double I_MOM_Gxy2z_S_Px_vrr = PAY*I_MOM_Fx2z_S_Px_vrr;
    Double I_MOM_Gx3z_S_Px_vrr = PAX*I_MOM_F3z_S_Px_vrr+oned2z*I_MOM_F3z_S_S_vrr;
    Double I_MOM_G4y_S_Px_vrr = PAY*I_MOM_F3y_S_Px_vrr+3*oned2z*I_MOM_D2y_S_Px_vrr;
    Double I_MOM_G3yz_S_Px_vrr = PAZ*I_MOM_F3y_S_Px_vrr;
    Double I_MOM_G2y2z_S_Px_vrr = PAZ*I_MOM_F2yz_S_Px_vrr+oned2z*I_MOM_D2y_S_Px_vrr;
    Double I_MOM_Gy3z_S_Px_vrr = PAY*I_MOM_F3z_S_Px_vrr;
    Double I_MOM_G4z_S_Px_vrr = PAZ*I_MOM_F3z_S_Px_vrr+3*oned2z*I_MOM_D2z_S_Px_vrr;
    Double I_MOM_G4x_S_Py_vrr = PAX*I_MOM_F3x_S_Py_vrr+3*oned2z*I_MOM_D2x_S_Py_vrr;
    Double I_MOM_G3xy_S_Py_vrr = PAY*I_MOM_F3x_S_Py_vrr+oned2z*I_MOM_F3x_S_S_vrr;
    Double I_MOM_G3xz_S_Py_vrr = PAZ*I_MOM_F3x_S_Py_vrr;
    Double I_MOM_G2x2y_S_Py_vrr = PAY*I_MOM_F2xy_S_Py_vrr+oned2z*I_MOM_D2x_S_Py_vrr+oned2z*I_MOM_F2xy_S_S_vrr;
    Double I_MOM_G2xyz_S_Py_vrr = PAZ*I_MOM_F2xy_S_Py_vrr;
    Double I_MOM_G2x2z_S_Py_vrr = PAZ*I_MOM_F2xz_S_Py_vrr+oned2z*I_MOM_D2x_S_Py_vrr;
    Double I_MOM_Gx3y_S_Py_vrr = PAX*I_MOM_F3y_S_Py_vrr;
    Double I_MOM_Gx2yz_S_Py_vrr = PAZ*I_MOM_Fx2y_S_Py_vrr;
    Double I_MOM_Gxy2z_S_Py_vrr = PAY*I_MOM_Fx2z_S_Py_vrr+oned2z*I_MOM_Fx2z_S_S_vrr;
    Double I_MOM_Gx3z_S_Py_vrr = PAX*I_MOM_F3z_S_Py_vrr;
    Double I_MOM_G4y_S_Py_vrr = PAY*I_MOM_F3y_S_Py_vrr+3*oned2z*I_MOM_D2y_S_Py_vrr+oned2z*I_MOM_F3y_S_S_vrr;
    Double I_MOM_G3yz_S_Py_vrr = PAZ*I_MOM_F3y_S_Py_vrr;
    Double I_MOM_G2y2z_S_Py_vrr = PAZ*I_MOM_F2yz_S_Py_vrr+oned2z*I_MOM_D2y_S_Py_vrr;
    Double I_MOM_Gy3z_S_Py_vrr = PAY*I_MOM_F3z_S_Py_vrr+oned2z*I_MOM_F3z_S_S_vrr;
    Double I_MOM_G4z_S_Py_vrr = PAZ*I_MOM_F3z_S_Py_vrr+3*oned2z*I_MOM_D2z_S_Py_vrr;
    Double I_MOM_G4x_S_Pz_vrr = PAX*I_MOM_F3x_S_Pz_vrr+3*oned2z*I_MOM_D2x_S_Pz_vrr;
    Double I_MOM_G3xy_S_Pz_vrr = PAY*I_MOM_F3x_S_Pz_vrr;
    Double I_MOM_G3xz_S_Pz_vrr = PAZ*I_MOM_F3x_S_Pz_vrr+oned2z*I_MOM_F3x_S_S_vrr;
    Double I_MOM_G2x2y_S_Pz_vrr = PAY*I_MOM_F2xy_S_Pz_vrr+oned2z*I_MOM_D2x_S_Pz_vrr;
    Double I_MOM_G2xyz_S_Pz_vrr = PAZ*I_MOM_F2xy_S_Pz_vrr+oned2z*I_MOM_F2xy_S_S_vrr;
    Double I_MOM_G2x2z_S_Pz_vrr = PAZ*I_MOM_F2xz_S_Pz_vrr+oned2z*I_MOM_D2x_S_Pz_vrr+oned2z*I_MOM_F2xz_S_S_vrr;
    Double I_MOM_Gx3y_S_Pz_vrr = PAX*I_MOM_F3y_S_Pz_vrr;
    Double I_MOM_Gx2yz_S_Pz_vrr = PAZ*I_MOM_Fx2y_S_Pz_vrr+oned2z*I_MOM_Fx2y_S_S_vrr;
    Double I_MOM_Gxy2z_S_Pz_vrr = PAY*I_MOM_Fx2z_S_Pz_vrr;
    Double I_MOM_Gx3z_S_Pz_vrr = PAX*I_MOM_F3z_S_Pz_vrr;
    Double I_MOM_G4y_S_Pz_vrr = PAY*I_MOM_F3y_S_Pz_vrr+3*oned2z*I_MOM_D2y_S_Pz_vrr;
    Double I_MOM_G3yz_S_Pz_vrr = PAZ*I_MOM_F3y_S_Pz_vrr+oned2z*I_MOM_F3y_S_S_vrr;
    Double I_MOM_G2y2z_S_Pz_vrr = PAZ*I_MOM_F2yz_S_Pz_vrr+oned2z*I_MOM_D2y_S_Pz_vrr+oned2z*I_MOM_F2yz_S_S_vrr;
    Double I_MOM_Gy3z_S_Pz_vrr = PAY*I_MOM_F3z_S_Pz_vrr;
    Double I_MOM_G4z_S_Pz_vrr = PAZ*I_MOM_F3z_S_Pz_vrr+3*oned2z*I_MOM_D2z_S_Pz_vrr+oned2z*I_MOM_F3z_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_G_S_P
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_MOM_G4x_S_Px += I_MOM_G4x_S_Px_vrr;
    I_MOM_G3xy_S_Px += I_MOM_G3xy_S_Px_vrr;
    I_MOM_G3xz_S_Px += I_MOM_G3xz_S_Px_vrr;
    I_MOM_G2x2y_S_Px += I_MOM_G2x2y_S_Px_vrr;
    I_MOM_G2xyz_S_Px += I_MOM_G2xyz_S_Px_vrr;
    I_MOM_G2x2z_S_Px += I_MOM_G2x2z_S_Px_vrr;
    I_MOM_Gx3y_S_Px += I_MOM_Gx3y_S_Px_vrr;
    I_MOM_Gx2yz_S_Px += I_MOM_Gx2yz_S_Px_vrr;
    I_MOM_Gxy2z_S_Px += I_MOM_Gxy2z_S_Px_vrr;
    I_MOM_Gx3z_S_Px += I_MOM_Gx3z_S_Px_vrr;
    I_MOM_G4y_S_Px += I_MOM_G4y_S_Px_vrr;
    I_MOM_G3yz_S_Px += I_MOM_G3yz_S_Px_vrr;
    I_MOM_G2y2z_S_Px += I_MOM_G2y2z_S_Px_vrr;
    I_MOM_Gy3z_S_Px += I_MOM_Gy3z_S_Px_vrr;
    I_MOM_G4z_S_Px += I_MOM_G4z_S_Px_vrr;
    I_MOM_G4x_S_Py += I_MOM_G4x_S_Py_vrr;
    I_MOM_G3xy_S_Py += I_MOM_G3xy_S_Py_vrr;
    I_MOM_G3xz_S_Py += I_MOM_G3xz_S_Py_vrr;
    I_MOM_G2x2y_S_Py += I_MOM_G2x2y_S_Py_vrr;
    I_MOM_G2xyz_S_Py += I_MOM_G2xyz_S_Py_vrr;
    I_MOM_G2x2z_S_Py += I_MOM_G2x2z_S_Py_vrr;
    I_MOM_Gx3y_S_Py += I_MOM_Gx3y_S_Py_vrr;
    I_MOM_Gx2yz_S_Py += I_MOM_Gx2yz_S_Py_vrr;
    I_MOM_Gxy2z_S_Py += I_MOM_Gxy2z_S_Py_vrr;
    I_MOM_Gx3z_S_Py += I_MOM_Gx3z_S_Py_vrr;
    I_MOM_G4y_S_Py += I_MOM_G4y_S_Py_vrr;
    I_MOM_G3yz_S_Py += I_MOM_G3yz_S_Py_vrr;
    I_MOM_G2y2z_S_Py += I_MOM_G2y2z_S_Py_vrr;
    I_MOM_Gy3z_S_Py += I_MOM_Gy3z_S_Py_vrr;
    I_MOM_G4z_S_Py += I_MOM_G4z_S_Py_vrr;
    I_MOM_G4x_S_Pz += I_MOM_G4x_S_Pz_vrr;
    I_MOM_G3xy_S_Pz += I_MOM_G3xy_S_Pz_vrr;
    I_MOM_G3xz_S_Pz += I_MOM_G3xz_S_Pz_vrr;
    I_MOM_G2x2y_S_Pz += I_MOM_G2x2y_S_Pz_vrr;
    I_MOM_G2xyz_S_Pz += I_MOM_G2xyz_S_Pz_vrr;
    I_MOM_G2x2z_S_Pz += I_MOM_G2x2z_S_Pz_vrr;
    I_MOM_Gx3y_S_Pz += I_MOM_Gx3y_S_Pz_vrr;
    I_MOM_Gx2yz_S_Pz += I_MOM_Gx2yz_S_Pz_vrr;
    I_MOM_Gxy2z_S_Pz += I_MOM_Gxy2z_S_Pz_vrr;
    I_MOM_Gx3z_S_Pz += I_MOM_Gx3z_S_Pz_vrr;
    I_MOM_G4y_S_Pz += I_MOM_G4y_S_Pz_vrr;
    I_MOM_G3yz_S_Pz += I_MOM_G3yz_S_Pz_vrr;
    I_MOM_G2y2z_S_Pz += I_MOM_G2y2z_S_Pz_vrr;
    I_MOM_Gy3z_S_Pz += I_MOM_Gy3z_S_Pz_vrr;
    I_MOM_G4z_S_Pz += I_MOM_G4z_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_F_S_P
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_MOM_F3x_S_Px += I_MOM_F3x_S_Px_vrr;
    I_MOM_F2xy_S_Px += I_MOM_F2xy_S_Px_vrr;
    I_MOM_F2xz_S_Px += I_MOM_F2xz_S_Px_vrr;
    I_MOM_Fx2y_S_Px += I_MOM_Fx2y_S_Px_vrr;
    I_MOM_Fxyz_S_Px += I_MOM_Fxyz_S_Px_vrr;
    I_MOM_Fx2z_S_Px += I_MOM_Fx2z_S_Px_vrr;
    I_MOM_F3y_S_Px += I_MOM_F3y_S_Px_vrr;
    I_MOM_F2yz_S_Px += I_MOM_F2yz_S_Px_vrr;
    I_MOM_Fy2z_S_Px += I_MOM_Fy2z_S_Px_vrr;
    I_MOM_F3z_S_Px += I_MOM_F3z_S_Px_vrr;
    I_MOM_F3x_S_Py += I_MOM_F3x_S_Py_vrr;
    I_MOM_F2xy_S_Py += I_MOM_F2xy_S_Py_vrr;
    I_MOM_F2xz_S_Py += I_MOM_F2xz_S_Py_vrr;
    I_MOM_Fx2y_S_Py += I_MOM_Fx2y_S_Py_vrr;
    I_MOM_Fxyz_S_Py += I_MOM_Fxyz_S_Py_vrr;
    I_MOM_Fx2z_S_Py += I_MOM_Fx2z_S_Py_vrr;
    I_MOM_F3y_S_Py += I_MOM_F3y_S_Py_vrr;
    I_MOM_F2yz_S_Py += I_MOM_F2yz_S_Py_vrr;
    I_MOM_Fy2z_S_Py += I_MOM_Fy2z_S_Py_vrr;
    I_MOM_F3z_S_Py += I_MOM_F3z_S_Py_vrr;
    I_MOM_F3x_S_Pz += I_MOM_F3x_S_Pz_vrr;
    I_MOM_F2xy_S_Pz += I_MOM_F2xy_S_Pz_vrr;
    I_MOM_F2xz_S_Pz += I_MOM_F2xz_S_Pz_vrr;
    I_MOM_Fx2y_S_Pz += I_MOM_Fx2y_S_Pz_vrr;
    I_MOM_Fxyz_S_Pz += I_MOM_Fxyz_S_Pz_vrr;
    I_MOM_Fx2z_S_Pz += I_MOM_Fx2z_S_Pz_vrr;
    I_MOM_F3y_S_Pz += I_MOM_F3y_S_Pz_vrr;
    I_MOM_F2yz_S_Pz += I_MOM_F2yz_S_Pz_vrr;
    I_MOM_Fy2z_S_Pz += I_MOM_Fy2z_S_Pz_vrr;
    I_MOM_F3z_S_Pz += I_MOM_F3z_S_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_MOM_D_S_P
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_MOM_D2x_S_Px += I_MOM_D2x_S_Px_vrr;
    I_MOM_Dxy_S_Px += I_MOM_Dxy_S_Px_vrr;
    I_MOM_Dxz_S_Px += I_MOM_Dxz_S_Px_vrr;
    I_MOM_D2y_S_Px += I_MOM_D2y_S_Px_vrr;
    I_MOM_Dyz_S_Px += I_MOM_Dyz_S_Px_vrr;
    I_MOM_D2z_S_Px += I_MOM_D2z_S_Px_vrr;
    I_MOM_D2x_S_Py += I_MOM_D2x_S_Py_vrr;
    I_MOM_Dxy_S_Py += I_MOM_Dxy_S_Py_vrr;
    I_MOM_Dxz_S_Py += I_MOM_Dxz_S_Py_vrr;
    I_MOM_D2y_S_Py += I_MOM_D2y_S_Py_vrr;
    I_MOM_Dyz_S_Py += I_MOM_Dyz_S_Py_vrr;
    I_MOM_D2z_S_Py += I_MOM_D2z_S_Py_vrr;
    I_MOM_D2x_S_Pz += I_MOM_D2x_S_Pz_vrr;
    I_MOM_Dxy_S_Pz += I_MOM_Dxy_S_Pz_vrr;
    I_MOM_Dxz_S_Pz += I_MOM_Dxz_S_Pz_vrr;
    I_MOM_D2y_S_Pz += I_MOM_D2y_S_Pz_vrr;
    I_MOM_Dyz_S_Pz += I_MOM_Dyz_S_Pz_vrr;
    I_MOM_D2z_S_Pz += I_MOM_D2z_S_Pz_vrr;
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
   * shell quartet name: SQ_MOM_D_P_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_MOM_F_S_P
   * RHS shell quartet name: SQ_MOM_D_S_P
   ************************************************************/
  Double I_MOM_D2x_Px_Px = I_MOM_F3x_S_Px+ABX*I_MOM_D2x_S_Px;
  Double I_MOM_Dxy_Px_Px = I_MOM_F2xy_S_Px+ABX*I_MOM_Dxy_S_Px;
  Double I_MOM_Dxz_Px_Px = I_MOM_F2xz_S_Px+ABX*I_MOM_Dxz_S_Px;
  Double I_MOM_D2y_Px_Px = I_MOM_Fx2y_S_Px+ABX*I_MOM_D2y_S_Px;
  Double I_MOM_Dyz_Px_Px = I_MOM_Fxyz_S_Px+ABX*I_MOM_Dyz_S_Px;
  Double I_MOM_D2z_Px_Px = I_MOM_Fx2z_S_Px+ABX*I_MOM_D2z_S_Px;
  Double I_MOM_D2x_Py_Px = I_MOM_F2xy_S_Px+ABY*I_MOM_D2x_S_Px;
  Double I_MOM_Dxy_Py_Px = I_MOM_Fx2y_S_Px+ABY*I_MOM_Dxy_S_Px;
  Double I_MOM_Dxz_Py_Px = I_MOM_Fxyz_S_Px+ABY*I_MOM_Dxz_S_Px;
  Double I_MOM_D2y_Py_Px = I_MOM_F3y_S_Px+ABY*I_MOM_D2y_S_Px;
  Double I_MOM_Dyz_Py_Px = I_MOM_F2yz_S_Px+ABY*I_MOM_Dyz_S_Px;
  Double I_MOM_D2z_Py_Px = I_MOM_Fy2z_S_Px+ABY*I_MOM_D2z_S_Px;
  Double I_MOM_D2x_Pz_Px = I_MOM_F2xz_S_Px+ABZ*I_MOM_D2x_S_Px;
  Double I_MOM_Dxy_Pz_Px = I_MOM_Fxyz_S_Px+ABZ*I_MOM_Dxy_S_Px;
  Double I_MOM_Dxz_Pz_Px = I_MOM_Fx2z_S_Px+ABZ*I_MOM_Dxz_S_Px;
  Double I_MOM_D2y_Pz_Px = I_MOM_F2yz_S_Px+ABZ*I_MOM_D2y_S_Px;
  Double I_MOM_Dyz_Pz_Px = I_MOM_Fy2z_S_Px+ABZ*I_MOM_Dyz_S_Px;
  Double I_MOM_D2z_Pz_Px = I_MOM_F3z_S_Px+ABZ*I_MOM_D2z_S_Px;
  Double I_MOM_D2x_Px_Py = I_MOM_F3x_S_Py+ABX*I_MOM_D2x_S_Py;
  Double I_MOM_Dxy_Px_Py = I_MOM_F2xy_S_Py+ABX*I_MOM_Dxy_S_Py;
  Double I_MOM_Dxz_Px_Py = I_MOM_F2xz_S_Py+ABX*I_MOM_Dxz_S_Py;
  Double I_MOM_D2y_Px_Py = I_MOM_Fx2y_S_Py+ABX*I_MOM_D2y_S_Py;
  Double I_MOM_Dyz_Px_Py = I_MOM_Fxyz_S_Py+ABX*I_MOM_Dyz_S_Py;
  Double I_MOM_D2z_Px_Py = I_MOM_Fx2z_S_Py+ABX*I_MOM_D2z_S_Py;
  Double I_MOM_D2x_Py_Py = I_MOM_F2xy_S_Py+ABY*I_MOM_D2x_S_Py;
  Double I_MOM_Dxy_Py_Py = I_MOM_Fx2y_S_Py+ABY*I_MOM_Dxy_S_Py;
  Double I_MOM_Dxz_Py_Py = I_MOM_Fxyz_S_Py+ABY*I_MOM_Dxz_S_Py;
  Double I_MOM_D2y_Py_Py = I_MOM_F3y_S_Py+ABY*I_MOM_D2y_S_Py;
  Double I_MOM_Dyz_Py_Py = I_MOM_F2yz_S_Py+ABY*I_MOM_Dyz_S_Py;
  Double I_MOM_D2z_Py_Py = I_MOM_Fy2z_S_Py+ABY*I_MOM_D2z_S_Py;
  Double I_MOM_D2x_Pz_Py = I_MOM_F2xz_S_Py+ABZ*I_MOM_D2x_S_Py;
  Double I_MOM_Dxy_Pz_Py = I_MOM_Fxyz_S_Py+ABZ*I_MOM_Dxy_S_Py;
  Double I_MOM_Dxz_Pz_Py = I_MOM_Fx2z_S_Py+ABZ*I_MOM_Dxz_S_Py;
  Double I_MOM_D2y_Pz_Py = I_MOM_F2yz_S_Py+ABZ*I_MOM_D2y_S_Py;
  Double I_MOM_Dyz_Pz_Py = I_MOM_Fy2z_S_Py+ABZ*I_MOM_Dyz_S_Py;
  Double I_MOM_D2z_Pz_Py = I_MOM_F3z_S_Py+ABZ*I_MOM_D2z_S_Py;
  Double I_MOM_D2x_Px_Pz = I_MOM_F3x_S_Pz+ABX*I_MOM_D2x_S_Pz;
  Double I_MOM_Dxy_Px_Pz = I_MOM_F2xy_S_Pz+ABX*I_MOM_Dxy_S_Pz;
  Double I_MOM_Dxz_Px_Pz = I_MOM_F2xz_S_Pz+ABX*I_MOM_Dxz_S_Pz;
  Double I_MOM_D2y_Px_Pz = I_MOM_Fx2y_S_Pz+ABX*I_MOM_D2y_S_Pz;
  Double I_MOM_Dyz_Px_Pz = I_MOM_Fxyz_S_Pz+ABX*I_MOM_Dyz_S_Pz;
  Double I_MOM_D2z_Px_Pz = I_MOM_Fx2z_S_Pz+ABX*I_MOM_D2z_S_Pz;
  Double I_MOM_D2x_Py_Pz = I_MOM_F2xy_S_Pz+ABY*I_MOM_D2x_S_Pz;
  Double I_MOM_Dxy_Py_Pz = I_MOM_Fx2y_S_Pz+ABY*I_MOM_Dxy_S_Pz;
  Double I_MOM_Dxz_Py_Pz = I_MOM_Fxyz_S_Pz+ABY*I_MOM_Dxz_S_Pz;
  Double I_MOM_D2y_Py_Pz = I_MOM_F3y_S_Pz+ABY*I_MOM_D2y_S_Pz;
  Double I_MOM_Dyz_Py_Pz = I_MOM_F2yz_S_Pz+ABY*I_MOM_Dyz_S_Pz;
  Double I_MOM_D2z_Py_Pz = I_MOM_Fy2z_S_Pz+ABY*I_MOM_D2z_S_Pz;
  Double I_MOM_D2x_Pz_Pz = I_MOM_F2xz_S_Pz+ABZ*I_MOM_D2x_S_Pz;
  Double I_MOM_Dxy_Pz_Pz = I_MOM_Fxyz_S_Pz+ABZ*I_MOM_Dxy_S_Pz;
  Double I_MOM_Dxz_Pz_Pz = I_MOM_Fx2z_S_Pz+ABZ*I_MOM_Dxz_S_Pz;
  Double I_MOM_D2y_Pz_Pz = I_MOM_F2yz_S_Pz+ABZ*I_MOM_D2y_S_Pz;
  Double I_MOM_Dyz_Pz_Pz = I_MOM_Fy2z_S_Pz+ABZ*I_MOM_Dyz_S_Pz;
  Double I_MOM_D2z_Pz_Pz = I_MOM_F3z_S_Pz+ABZ*I_MOM_D2z_S_Pz;

  /************************************************************
   * shell quartet name: SQ_MOM_F_P_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 15 integrals are omitted 
   * RHS shell quartet name: SQ_MOM_G_S_P
   * RHS shell quartet name: SQ_MOM_F_S_P
   ************************************************************/
  Double I_MOM_F3x_Px_Px = I_MOM_G4x_S_Px+ABX*I_MOM_F3x_S_Px;
  Double I_MOM_F2xy_Px_Px = I_MOM_G3xy_S_Px+ABX*I_MOM_F2xy_S_Px;
  Double I_MOM_F2xz_Px_Px = I_MOM_G3xz_S_Px+ABX*I_MOM_F2xz_S_Px;
  Double I_MOM_Fx2y_Px_Px = I_MOM_G2x2y_S_Px+ABX*I_MOM_Fx2y_S_Px;
  Double I_MOM_Fxyz_Px_Px = I_MOM_G2xyz_S_Px+ABX*I_MOM_Fxyz_S_Px;
  Double I_MOM_Fx2z_Px_Px = I_MOM_G2x2z_S_Px+ABX*I_MOM_Fx2z_S_Px;
  Double I_MOM_F3y_Px_Px = I_MOM_Gx3y_S_Px+ABX*I_MOM_F3y_S_Px;
  Double I_MOM_F2yz_Px_Px = I_MOM_Gx2yz_S_Px+ABX*I_MOM_F2yz_S_Px;
  Double I_MOM_Fy2z_Px_Px = I_MOM_Gxy2z_S_Px+ABX*I_MOM_Fy2z_S_Px;
  Double I_MOM_F3z_Px_Px = I_MOM_Gx3z_S_Px+ABX*I_MOM_F3z_S_Px;
  Double I_MOM_F2xy_Py_Px = I_MOM_G2x2y_S_Px+ABY*I_MOM_F2xy_S_Px;
  Double I_MOM_F2xz_Py_Px = I_MOM_G2xyz_S_Px+ABY*I_MOM_F2xz_S_Px;
  Double I_MOM_Fx2y_Py_Px = I_MOM_Gx3y_S_Px+ABY*I_MOM_Fx2y_S_Px;
  Double I_MOM_Fxyz_Py_Px = I_MOM_Gx2yz_S_Px+ABY*I_MOM_Fxyz_S_Px;
  Double I_MOM_Fx2z_Py_Px = I_MOM_Gxy2z_S_Px+ABY*I_MOM_Fx2z_S_Px;
  Double I_MOM_F3y_Py_Px = I_MOM_G4y_S_Px+ABY*I_MOM_F3y_S_Px;
  Double I_MOM_F2yz_Py_Px = I_MOM_G3yz_S_Px+ABY*I_MOM_F2yz_S_Px;
  Double I_MOM_Fy2z_Py_Px = I_MOM_G2y2z_S_Px+ABY*I_MOM_Fy2z_S_Px;
  Double I_MOM_F3z_Py_Px = I_MOM_Gy3z_S_Px+ABY*I_MOM_F3z_S_Px;
  Double I_MOM_F2xz_Pz_Px = I_MOM_G2x2z_S_Px+ABZ*I_MOM_F2xz_S_Px;
  Double I_MOM_Fxyz_Pz_Px = I_MOM_Gxy2z_S_Px+ABZ*I_MOM_Fxyz_S_Px;
  Double I_MOM_Fx2z_Pz_Px = I_MOM_Gx3z_S_Px+ABZ*I_MOM_Fx2z_S_Px;
  Double I_MOM_F2yz_Pz_Px = I_MOM_G2y2z_S_Px+ABZ*I_MOM_F2yz_S_Px;
  Double I_MOM_Fy2z_Pz_Px = I_MOM_Gy3z_S_Px+ABZ*I_MOM_Fy2z_S_Px;
  Double I_MOM_F3z_Pz_Px = I_MOM_G4z_S_Px+ABZ*I_MOM_F3z_S_Px;
  Double I_MOM_F3x_Px_Py = I_MOM_G4x_S_Py+ABX*I_MOM_F3x_S_Py;
  Double I_MOM_F2xy_Px_Py = I_MOM_G3xy_S_Py+ABX*I_MOM_F2xy_S_Py;
  Double I_MOM_F2xz_Px_Py = I_MOM_G3xz_S_Py+ABX*I_MOM_F2xz_S_Py;
  Double I_MOM_Fx2y_Px_Py = I_MOM_G2x2y_S_Py+ABX*I_MOM_Fx2y_S_Py;
  Double I_MOM_Fxyz_Px_Py = I_MOM_G2xyz_S_Py+ABX*I_MOM_Fxyz_S_Py;
  Double I_MOM_Fx2z_Px_Py = I_MOM_G2x2z_S_Py+ABX*I_MOM_Fx2z_S_Py;
  Double I_MOM_F3y_Px_Py = I_MOM_Gx3y_S_Py+ABX*I_MOM_F3y_S_Py;
  Double I_MOM_F2yz_Px_Py = I_MOM_Gx2yz_S_Py+ABX*I_MOM_F2yz_S_Py;
  Double I_MOM_Fy2z_Px_Py = I_MOM_Gxy2z_S_Py+ABX*I_MOM_Fy2z_S_Py;
  Double I_MOM_F3z_Px_Py = I_MOM_Gx3z_S_Py+ABX*I_MOM_F3z_S_Py;
  Double I_MOM_F2xy_Py_Py = I_MOM_G2x2y_S_Py+ABY*I_MOM_F2xy_S_Py;
  Double I_MOM_F2xz_Py_Py = I_MOM_G2xyz_S_Py+ABY*I_MOM_F2xz_S_Py;
  Double I_MOM_Fx2y_Py_Py = I_MOM_Gx3y_S_Py+ABY*I_MOM_Fx2y_S_Py;
  Double I_MOM_Fxyz_Py_Py = I_MOM_Gx2yz_S_Py+ABY*I_MOM_Fxyz_S_Py;
  Double I_MOM_Fx2z_Py_Py = I_MOM_Gxy2z_S_Py+ABY*I_MOM_Fx2z_S_Py;
  Double I_MOM_F3y_Py_Py = I_MOM_G4y_S_Py+ABY*I_MOM_F3y_S_Py;
  Double I_MOM_F2yz_Py_Py = I_MOM_G3yz_S_Py+ABY*I_MOM_F2yz_S_Py;
  Double I_MOM_Fy2z_Py_Py = I_MOM_G2y2z_S_Py+ABY*I_MOM_Fy2z_S_Py;
  Double I_MOM_F3z_Py_Py = I_MOM_Gy3z_S_Py+ABY*I_MOM_F3z_S_Py;
  Double I_MOM_F2xz_Pz_Py = I_MOM_G2x2z_S_Py+ABZ*I_MOM_F2xz_S_Py;
  Double I_MOM_Fxyz_Pz_Py = I_MOM_Gxy2z_S_Py+ABZ*I_MOM_Fxyz_S_Py;
  Double I_MOM_Fx2z_Pz_Py = I_MOM_Gx3z_S_Py+ABZ*I_MOM_Fx2z_S_Py;
  Double I_MOM_F2yz_Pz_Py = I_MOM_G2y2z_S_Py+ABZ*I_MOM_F2yz_S_Py;
  Double I_MOM_Fy2z_Pz_Py = I_MOM_Gy3z_S_Py+ABZ*I_MOM_Fy2z_S_Py;
  Double I_MOM_F3z_Pz_Py = I_MOM_G4z_S_Py+ABZ*I_MOM_F3z_S_Py;
  Double I_MOM_F3x_Px_Pz = I_MOM_G4x_S_Pz+ABX*I_MOM_F3x_S_Pz;
  Double I_MOM_F2xy_Px_Pz = I_MOM_G3xy_S_Pz+ABX*I_MOM_F2xy_S_Pz;
  Double I_MOM_F2xz_Px_Pz = I_MOM_G3xz_S_Pz+ABX*I_MOM_F2xz_S_Pz;
  Double I_MOM_Fx2y_Px_Pz = I_MOM_G2x2y_S_Pz+ABX*I_MOM_Fx2y_S_Pz;
  Double I_MOM_Fxyz_Px_Pz = I_MOM_G2xyz_S_Pz+ABX*I_MOM_Fxyz_S_Pz;
  Double I_MOM_Fx2z_Px_Pz = I_MOM_G2x2z_S_Pz+ABX*I_MOM_Fx2z_S_Pz;
  Double I_MOM_F3y_Px_Pz = I_MOM_Gx3y_S_Pz+ABX*I_MOM_F3y_S_Pz;
  Double I_MOM_F2yz_Px_Pz = I_MOM_Gx2yz_S_Pz+ABX*I_MOM_F2yz_S_Pz;
  Double I_MOM_Fy2z_Px_Pz = I_MOM_Gxy2z_S_Pz+ABX*I_MOM_Fy2z_S_Pz;
  Double I_MOM_F3z_Px_Pz = I_MOM_Gx3z_S_Pz+ABX*I_MOM_F3z_S_Pz;
  Double I_MOM_F2xy_Py_Pz = I_MOM_G2x2y_S_Pz+ABY*I_MOM_F2xy_S_Pz;
  Double I_MOM_F2xz_Py_Pz = I_MOM_G2xyz_S_Pz+ABY*I_MOM_F2xz_S_Pz;
  Double I_MOM_Fx2y_Py_Pz = I_MOM_Gx3y_S_Pz+ABY*I_MOM_Fx2y_S_Pz;
  Double I_MOM_Fxyz_Py_Pz = I_MOM_Gx2yz_S_Pz+ABY*I_MOM_Fxyz_S_Pz;
  Double I_MOM_Fx2z_Py_Pz = I_MOM_Gxy2z_S_Pz+ABY*I_MOM_Fx2z_S_Pz;
  Double I_MOM_F3y_Py_Pz = I_MOM_G4y_S_Pz+ABY*I_MOM_F3y_S_Pz;
  Double I_MOM_F2yz_Py_Pz = I_MOM_G3yz_S_Pz+ABY*I_MOM_F2yz_S_Pz;
  Double I_MOM_Fy2z_Py_Pz = I_MOM_G2y2z_S_Pz+ABY*I_MOM_Fy2z_S_Pz;
  Double I_MOM_F3z_Py_Pz = I_MOM_Gy3z_S_Pz+ABY*I_MOM_F3z_S_Pz;
  Double I_MOM_F2xz_Pz_Pz = I_MOM_G2x2z_S_Pz+ABZ*I_MOM_F2xz_S_Pz;
  Double I_MOM_Fxyz_Pz_Pz = I_MOM_Gxy2z_S_Pz+ABZ*I_MOM_Fxyz_S_Pz;
  Double I_MOM_Fx2z_Pz_Pz = I_MOM_Gx3z_S_Pz+ABZ*I_MOM_Fx2z_S_Pz;
  Double I_MOM_F2yz_Pz_Pz = I_MOM_G2y2z_S_Pz+ABZ*I_MOM_F2yz_S_Pz;
  Double I_MOM_Fy2z_Pz_Pz = I_MOM_Gy3z_S_Pz+ABZ*I_MOM_Fy2z_S_Pz;
  Double I_MOM_F3z_Pz_Pz = I_MOM_G4z_S_Pz+ABZ*I_MOM_F3z_S_Pz;

  /************************************************************
   * shell quartet name: SQ_MOM_D_D_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_MOM_F_P_P
   * RHS shell quartet name: SQ_MOM_D_P_P
   ************************************************************/
  abcd[0] = I_MOM_F3x_Px_Px+ABX*I_MOM_D2x_Px_Px;
  abcd[1] = I_MOM_F2xy_Px_Px+ABX*I_MOM_Dxy_Px_Px;
  abcd[2] = I_MOM_F2xz_Px_Px+ABX*I_MOM_Dxz_Px_Px;
  abcd[3] = I_MOM_Fx2y_Px_Px+ABX*I_MOM_D2y_Px_Px;
  abcd[4] = I_MOM_Fxyz_Px_Px+ABX*I_MOM_Dyz_Px_Px;
  abcd[5] = I_MOM_Fx2z_Px_Px+ABX*I_MOM_D2z_Px_Px;
  abcd[6] = I_MOM_F2xy_Px_Px+ABY*I_MOM_D2x_Px_Px;
  abcd[7] = I_MOM_Fx2y_Px_Px+ABY*I_MOM_Dxy_Px_Px;
  abcd[8] = I_MOM_Fxyz_Px_Px+ABY*I_MOM_Dxz_Px_Px;
  abcd[9] = I_MOM_F3y_Px_Px+ABY*I_MOM_D2y_Px_Px;
  abcd[10] = I_MOM_F2yz_Px_Px+ABY*I_MOM_Dyz_Px_Px;
  abcd[11] = I_MOM_Fy2z_Px_Px+ABY*I_MOM_D2z_Px_Px;
  abcd[12] = I_MOM_F2xz_Px_Px+ABZ*I_MOM_D2x_Px_Px;
  abcd[13] = I_MOM_Fxyz_Px_Px+ABZ*I_MOM_Dxy_Px_Px;
  abcd[14] = I_MOM_Fx2z_Px_Px+ABZ*I_MOM_Dxz_Px_Px;
  abcd[15] = I_MOM_F2yz_Px_Px+ABZ*I_MOM_D2y_Px_Px;
  abcd[16] = I_MOM_Fy2z_Px_Px+ABZ*I_MOM_Dyz_Px_Px;
  abcd[17] = I_MOM_F3z_Px_Px+ABZ*I_MOM_D2z_Px_Px;
  abcd[18] = I_MOM_F2xy_Py_Px+ABY*I_MOM_D2x_Py_Px;
  abcd[19] = I_MOM_Fx2y_Py_Px+ABY*I_MOM_Dxy_Py_Px;
  abcd[20] = I_MOM_Fxyz_Py_Px+ABY*I_MOM_Dxz_Py_Px;
  abcd[21] = I_MOM_F3y_Py_Px+ABY*I_MOM_D2y_Py_Px;
  abcd[22] = I_MOM_F2yz_Py_Px+ABY*I_MOM_Dyz_Py_Px;
  abcd[23] = I_MOM_Fy2z_Py_Px+ABY*I_MOM_D2z_Py_Px;
  abcd[24] = I_MOM_F2xz_Py_Px+ABZ*I_MOM_D2x_Py_Px;
  abcd[25] = I_MOM_Fxyz_Py_Px+ABZ*I_MOM_Dxy_Py_Px;
  abcd[26] = I_MOM_Fx2z_Py_Px+ABZ*I_MOM_Dxz_Py_Px;
  abcd[27] = I_MOM_F2yz_Py_Px+ABZ*I_MOM_D2y_Py_Px;
  abcd[28] = I_MOM_Fy2z_Py_Px+ABZ*I_MOM_Dyz_Py_Px;
  abcd[29] = I_MOM_F3z_Py_Px+ABZ*I_MOM_D2z_Py_Px;
  abcd[30] = I_MOM_F2xz_Pz_Px+ABZ*I_MOM_D2x_Pz_Px;
  abcd[31] = I_MOM_Fxyz_Pz_Px+ABZ*I_MOM_Dxy_Pz_Px;
  abcd[32] = I_MOM_Fx2z_Pz_Px+ABZ*I_MOM_Dxz_Pz_Px;
  abcd[33] = I_MOM_F2yz_Pz_Px+ABZ*I_MOM_D2y_Pz_Px;
  abcd[34] = I_MOM_Fy2z_Pz_Px+ABZ*I_MOM_Dyz_Pz_Px;
  abcd[35] = I_MOM_F3z_Pz_Px+ABZ*I_MOM_D2z_Pz_Px;
  abcd[36] = I_MOM_F3x_Px_Py+ABX*I_MOM_D2x_Px_Py;
  abcd[37] = I_MOM_F2xy_Px_Py+ABX*I_MOM_Dxy_Px_Py;
  abcd[38] = I_MOM_F2xz_Px_Py+ABX*I_MOM_Dxz_Px_Py;
  abcd[39] = I_MOM_Fx2y_Px_Py+ABX*I_MOM_D2y_Px_Py;
  abcd[40] = I_MOM_Fxyz_Px_Py+ABX*I_MOM_Dyz_Px_Py;
  abcd[41] = I_MOM_Fx2z_Px_Py+ABX*I_MOM_D2z_Px_Py;
  abcd[42] = I_MOM_F2xy_Px_Py+ABY*I_MOM_D2x_Px_Py;
  abcd[43] = I_MOM_Fx2y_Px_Py+ABY*I_MOM_Dxy_Px_Py;
  abcd[44] = I_MOM_Fxyz_Px_Py+ABY*I_MOM_Dxz_Px_Py;
  abcd[45] = I_MOM_F3y_Px_Py+ABY*I_MOM_D2y_Px_Py;
  abcd[46] = I_MOM_F2yz_Px_Py+ABY*I_MOM_Dyz_Px_Py;
  abcd[47] = I_MOM_Fy2z_Px_Py+ABY*I_MOM_D2z_Px_Py;
  abcd[48] = I_MOM_F2xz_Px_Py+ABZ*I_MOM_D2x_Px_Py;
  abcd[49] = I_MOM_Fxyz_Px_Py+ABZ*I_MOM_Dxy_Px_Py;
  abcd[50] = I_MOM_Fx2z_Px_Py+ABZ*I_MOM_Dxz_Px_Py;
  abcd[51] = I_MOM_F2yz_Px_Py+ABZ*I_MOM_D2y_Px_Py;
  abcd[52] = I_MOM_Fy2z_Px_Py+ABZ*I_MOM_Dyz_Px_Py;
  abcd[53] = I_MOM_F3z_Px_Py+ABZ*I_MOM_D2z_Px_Py;
  abcd[54] = I_MOM_F2xy_Py_Py+ABY*I_MOM_D2x_Py_Py;
  abcd[55] = I_MOM_Fx2y_Py_Py+ABY*I_MOM_Dxy_Py_Py;
  abcd[56] = I_MOM_Fxyz_Py_Py+ABY*I_MOM_Dxz_Py_Py;
  abcd[57] = I_MOM_F3y_Py_Py+ABY*I_MOM_D2y_Py_Py;
  abcd[58] = I_MOM_F2yz_Py_Py+ABY*I_MOM_Dyz_Py_Py;
  abcd[59] = I_MOM_Fy2z_Py_Py+ABY*I_MOM_D2z_Py_Py;
  abcd[60] = I_MOM_F2xz_Py_Py+ABZ*I_MOM_D2x_Py_Py;
  abcd[61] = I_MOM_Fxyz_Py_Py+ABZ*I_MOM_Dxy_Py_Py;
  abcd[62] = I_MOM_Fx2z_Py_Py+ABZ*I_MOM_Dxz_Py_Py;
  abcd[63] = I_MOM_F2yz_Py_Py+ABZ*I_MOM_D2y_Py_Py;
  abcd[64] = I_MOM_Fy2z_Py_Py+ABZ*I_MOM_Dyz_Py_Py;
  abcd[65] = I_MOM_F3z_Py_Py+ABZ*I_MOM_D2z_Py_Py;
  abcd[66] = I_MOM_F2xz_Pz_Py+ABZ*I_MOM_D2x_Pz_Py;
  abcd[67] = I_MOM_Fxyz_Pz_Py+ABZ*I_MOM_Dxy_Pz_Py;
  abcd[68] = I_MOM_Fx2z_Pz_Py+ABZ*I_MOM_Dxz_Pz_Py;
  abcd[69] = I_MOM_F2yz_Pz_Py+ABZ*I_MOM_D2y_Pz_Py;
  abcd[70] = I_MOM_Fy2z_Pz_Py+ABZ*I_MOM_Dyz_Pz_Py;
  abcd[71] = I_MOM_F3z_Pz_Py+ABZ*I_MOM_D2z_Pz_Py;
  abcd[72] = I_MOM_F3x_Px_Pz+ABX*I_MOM_D2x_Px_Pz;
  abcd[73] = I_MOM_F2xy_Px_Pz+ABX*I_MOM_Dxy_Px_Pz;
  abcd[74] = I_MOM_F2xz_Px_Pz+ABX*I_MOM_Dxz_Px_Pz;
  abcd[75] = I_MOM_Fx2y_Px_Pz+ABX*I_MOM_D2y_Px_Pz;
  abcd[76] = I_MOM_Fxyz_Px_Pz+ABX*I_MOM_Dyz_Px_Pz;
  abcd[77] = I_MOM_Fx2z_Px_Pz+ABX*I_MOM_D2z_Px_Pz;
  abcd[78] = I_MOM_F2xy_Px_Pz+ABY*I_MOM_D2x_Px_Pz;
  abcd[79] = I_MOM_Fx2y_Px_Pz+ABY*I_MOM_Dxy_Px_Pz;
  abcd[80] = I_MOM_Fxyz_Px_Pz+ABY*I_MOM_Dxz_Px_Pz;
  abcd[81] = I_MOM_F3y_Px_Pz+ABY*I_MOM_D2y_Px_Pz;
  abcd[82] = I_MOM_F2yz_Px_Pz+ABY*I_MOM_Dyz_Px_Pz;
  abcd[83] = I_MOM_Fy2z_Px_Pz+ABY*I_MOM_D2z_Px_Pz;
  abcd[84] = I_MOM_F2xz_Px_Pz+ABZ*I_MOM_D2x_Px_Pz;
  abcd[85] = I_MOM_Fxyz_Px_Pz+ABZ*I_MOM_Dxy_Px_Pz;
  abcd[86] = I_MOM_Fx2z_Px_Pz+ABZ*I_MOM_Dxz_Px_Pz;
  abcd[87] = I_MOM_F2yz_Px_Pz+ABZ*I_MOM_D2y_Px_Pz;
  abcd[88] = I_MOM_Fy2z_Px_Pz+ABZ*I_MOM_Dyz_Px_Pz;
  abcd[89] = I_MOM_F3z_Px_Pz+ABZ*I_MOM_D2z_Px_Pz;
  abcd[90] = I_MOM_F2xy_Py_Pz+ABY*I_MOM_D2x_Py_Pz;
  abcd[91] = I_MOM_Fx2y_Py_Pz+ABY*I_MOM_Dxy_Py_Pz;
  abcd[92] = I_MOM_Fxyz_Py_Pz+ABY*I_MOM_Dxz_Py_Pz;
  abcd[93] = I_MOM_F3y_Py_Pz+ABY*I_MOM_D2y_Py_Pz;
  abcd[94] = I_MOM_F2yz_Py_Pz+ABY*I_MOM_Dyz_Py_Pz;
  abcd[95] = I_MOM_Fy2z_Py_Pz+ABY*I_MOM_D2z_Py_Pz;
  abcd[96] = I_MOM_F2xz_Py_Pz+ABZ*I_MOM_D2x_Py_Pz;
  abcd[97] = I_MOM_Fxyz_Py_Pz+ABZ*I_MOM_Dxy_Py_Pz;
  abcd[98] = I_MOM_Fx2z_Py_Pz+ABZ*I_MOM_Dxz_Py_Pz;
  abcd[99] = I_MOM_F2yz_Py_Pz+ABZ*I_MOM_D2y_Py_Pz;
  abcd[100] = I_MOM_Fy2z_Py_Pz+ABZ*I_MOM_Dyz_Py_Pz;
  abcd[101] = I_MOM_F3z_Py_Pz+ABZ*I_MOM_D2z_Py_Pz;
  abcd[102] = I_MOM_F2xz_Pz_Pz+ABZ*I_MOM_D2x_Pz_Pz;
  abcd[103] = I_MOM_Fxyz_Pz_Pz+ABZ*I_MOM_Dxy_Pz_Pz;
  abcd[104] = I_MOM_Fx2z_Pz_Pz+ABZ*I_MOM_Dxz_Pz_Pz;
  abcd[105] = I_MOM_F2yz_Pz_Pz+ABZ*I_MOM_D2y_Pz_Pz;
  abcd[106] = I_MOM_Fy2z_Pz_Pz+ABZ*I_MOM_Dyz_Pz_Pz;
  abcd[107] = I_MOM_F3z_Pz_Pz+ABZ*I_MOM_D2z_Pz_Pz;
}
