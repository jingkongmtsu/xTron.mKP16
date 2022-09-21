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

void hgp_os_twobodyoverlap_g_f(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_K7x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K5xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4xy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3xy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x4yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x3y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x2y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2xy4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx5yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx4y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx3y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx2y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kxy5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K7y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K6yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K5y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2y5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ky6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K7z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I6x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H5x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G4x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double fbra  = ifac[ip2];
    Double onedz = iexp[ip2];
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;
    if(fabs(I_TWOBODYOVERLAP_S_S_vrr)<THRESHOLD_MATH) continue;


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
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_S_vrr = PAZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_S_vrr = PAX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_S_vrr = PAY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_S_vrr = PAX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_S_vrr = PAX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_S_vrr = PAY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2y_S_vrr;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_K7x_S += I_TWOBODYOVERLAP_K7x_S_vrr;
    I_TWOBODYOVERLAP_K6xy_S += I_TWOBODYOVERLAP_K6xy_S_vrr;
    I_TWOBODYOVERLAP_K6xz_S += I_TWOBODYOVERLAP_K6xz_S_vrr;
    I_TWOBODYOVERLAP_K5x2y_S += I_TWOBODYOVERLAP_K5x2y_S_vrr;
    I_TWOBODYOVERLAP_K5xyz_S += I_TWOBODYOVERLAP_K5xyz_S_vrr;
    I_TWOBODYOVERLAP_K5x2z_S += I_TWOBODYOVERLAP_K5x2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3y_S += I_TWOBODYOVERLAP_K4x3y_S_vrr;
    I_TWOBODYOVERLAP_K4x2yz_S += I_TWOBODYOVERLAP_K4x2yz_S_vrr;
    I_TWOBODYOVERLAP_K4xy2z_S += I_TWOBODYOVERLAP_K4xy2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3z_S += I_TWOBODYOVERLAP_K4x3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4y_S += I_TWOBODYOVERLAP_K3x4y_S_vrr;
    I_TWOBODYOVERLAP_K3x3yz_S += I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    I_TWOBODYOVERLAP_K3x2y2z_S += I_TWOBODYOVERLAP_K3x2y2z_S_vrr;
    I_TWOBODYOVERLAP_K3xy3z_S += I_TWOBODYOVERLAP_K3xy3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4z_S += I_TWOBODYOVERLAP_K3x4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5y_S += I_TWOBODYOVERLAP_K2x5y_S_vrr;
    I_TWOBODYOVERLAP_K2x4yz_S += I_TWOBODYOVERLAP_K2x4yz_S_vrr;
    I_TWOBODYOVERLAP_K2x3y2z_S += I_TWOBODYOVERLAP_K2x3y2z_S_vrr;
    I_TWOBODYOVERLAP_K2x2y3z_S += I_TWOBODYOVERLAP_K2x2y3z_S_vrr;
    I_TWOBODYOVERLAP_K2xy4z_S += I_TWOBODYOVERLAP_K2xy4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5z_S += I_TWOBODYOVERLAP_K2x5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6y_S += I_TWOBODYOVERLAP_Kx6y_S_vrr;
    I_TWOBODYOVERLAP_Kx5yz_S += I_TWOBODYOVERLAP_Kx5yz_S_vrr;
    I_TWOBODYOVERLAP_Kx4y2z_S += I_TWOBODYOVERLAP_Kx4y2z_S_vrr;
    I_TWOBODYOVERLAP_Kx3y3z_S += I_TWOBODYOVERLAP_Kx3y3z_S_vrr;
    I_TWOBODYOVERLAP_Kx2y4z_S += I_TWOBODYOVERLAP_Kx2y4z_S_vrr;
    I_TWOBODYOVERLAP_Kxy5z_S += I_TWOBODYOVERLAP_Kxy5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6z_S += I_TWOBODYOVERLAP_Kx6z_S_vrr;
    I_TWOBODYOVERLAP_K7y_S += I_TWOBODYOVERLAP_K7y_S_vrr;
    I_TWOBODYOVERLAP_K6yz_S += I_TWOBODYOVERLAP_K6yz_S_vrr;
    I_TWOBODYOVERLAP_K5y2z_S += I_TWOBODYOVERLAP_K5y2z_S_vrr;
    I_TWOBODYOVERLAP_K4y3z_S += I_TWOBODYOVERLAP_K4y3z_S_vrr;
    I_TWOBODYOVERLAP_K3y4z_S += I_TWOBODYOVERLAP_K3y4z_S_vrr;
    I_TWOBODYOVERLAP_K2y5z_S += I_TWOBODYOVERLAP_K2y5z_S_vrr;
    I_TWOBODYOVERLAP_Ky6z_S += I_TWOBODYOVERLAP_Ky6z_S_vrr;
    I_TWOBODYOVERLAP_K7z_S += I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_I6x_S += I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S += I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S += I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S += I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S += I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S += I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S += I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S += I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S += I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S += I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S += I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S += I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S += I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S += I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S += I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S += I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S += I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S += I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S += I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S += I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S += I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S += I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S += I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S += I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S += I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S += I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S += I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S += I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_H5x_S += I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S += I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S += I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S += I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S += I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S += I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S += I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S += I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S += I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S += I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S += I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S += I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S += I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S += I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S += I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S += I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S += I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S += I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S += I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S += I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S += I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_G4x_S += I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S += I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S += I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S += I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S += I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S += I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S += I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S += I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S += I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S += I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S += I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S += I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S += I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S += I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S += I_TWOBODYOVERLAP_G4z_S_vrr;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_G_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_Px = I_TWOBODYOVERLAP_H5x_S+ABX*I_TWOBODYOVERLAP_G4x_S;
  Double I_TWOBODYOVERLAP_G3xy_Px = I_TWOBODYOVERLAP_H4xy_S+ABX*I_TWOBODYOVERLAP_G3xy_S;
  Double I_TWOBODYOVERLAP_G3xz_Px = I_TWOBODYOVERLAP_H4xz_S+ABX*I_TWOBODYOVERLAP_G3xz_S;
  Double I_TWOBODYOVERLAP_G2x2y_Px = I_TWOBODYOVERLAP_H3x2y_S+ABX*I_TWOBODYOVERLAP_G2x2y_S;
  Double I_TWOBODYOVERLAP_G2xyz_Px = I_TWOBODYOVERLAP_H3xyz_S+ABX*I_TWOBODYOVERLAP_G2xyz_S;
  Double I_TWOBODYOVERLAP_G2x2z_Px = I_TWOBODYOVERLAP_H3x2z_S+ABX*I_TWOBODYOVERLAP_G2x2z_S;
  Double I_TWOBODYOVERLAP_Gx3y_Px = I_TWOBODYOVERLAP_H2x3y_S+ABX*I_TWOBODYOVERLAP_Gx3y_S;
  Double I_TWOBODYOVERLAP_Gx2yz_Px = I_TWOBODYOVERLAP_H2x2yz_S+ABX*I_TWOBODYOVERLAP_Gx2yz_S;
  Double I_TWOBODYOVERLAP_Gxy2z_Px = I_TWOBODYOVERLAP_H2xy2z_S+ABX*I_TWOBODYOVERLAP_Gxy2z_S;
  Double I_TWOBODYOVERLAP_Gx3z_Px = I_TWOBODYOVERLAP_H2x3z_S+ABX*I_TWOBODYOVERLAP_Gx3z_S;
  Double I_TWOBODYOVERLAP_G4y_Px = I_TWOBODYOVERLAP_Hx4y_S+ABX*I_TWOBODYOVERLAP_G4y_S;
  Double I_TWOBODYOVERLAP_G3yz_Px = I_TWOBODYOVERLAP_Hx3yz_S+ABX*I_TWOBODYOVERLAP_G3yz_S;
  Double I_TWOBODYOVERLAP_G2y2z_Px = I_TWOBODYOVERLAP_Hx2y2z_S+ABX*I_TWOBODYOVERLAP_G2y2z_S;
  Double I_TWOBODYOVERLAP_Gy3z_Px = I_TWOBODYOVERLAP_Hxy3z_S+ABX*I_TWOBODYOVERLAP_Gy3z_S;
  Double I_TWOBODYOVERLAP_G4z_Px = I_TWOBODYOVERLAP_Hx4z_S+ABX*I_TWOBODYOVERLAP_G4z_S;
  Double I_TWOBODYOVERLAP_G4x_Py = I_TWOBODYOVERLAP_H4xy_S+ABY*I_TWOBODYOVERLAP_G4x_S;
  Double I_TWOBODYOVERLAP_G3xy_Py = I_TWOBODYOVERLAP_H3x2y_S+ABY*I_TWOBODYOVERLAP_G3xy_S;
  Double I_TWOBODYOVERLAP_G3xz_Py = I_TWOBODYOVERLAP_H3xyz_S+ABY*I_TWOBODYOVERLAP_G3xz_S;
  Double I_TWOBODYOVERLAP_G2x2y_Py = I_TWOBODYOVERLAP_H2x3y_S+ABY*I_TWOBODYOVERLAP_G2x2y_S;
  Double I_TWOBODYOVERLAP_G2xyz_Py = I_TWOBODYOVERLAP_H2x2yz_S+ABY*I_TWOBODYOVERLAP_G2xyz_S;
  Double I_TWOBODYOVERLAP_G2x2z_Py = I_TWOBODYOVERLAP_H2xy2z_S+ABY*I_TWOBODYOVERLAP_G2x2z_S;
  Double I_TWOBODYOVERLAP_Gx3y_Py = I_TWOBODYOVERLAP_Hx4y_S+ABY*I_TWOBODYOVERLAP_Gx3y_S;
  Double I_TWOBODYOVERLAP_Gx2yz_Py = I_TWOBODYOVERLAP_Hx3yz_S+ABY*I_TWOBODYOVERLAP_Gx2yz_S;
  Double I_TWOBODYOVERLAP_Gxy2z_Py = I_TWOBODYOVERLAP_Hx2y2z_S+ABY*I_TWOBODYOVERLAP_Gxy2z_S;
  Double I_TWOBODYOVERLAP_Gx3z_Py = I_TWOBODYOVERLAP_Hxy3z_S+ABY*I_TWOBODYOVERLAP_Gx3z_S;
  Double I_TWOBODYOVERLAP_G4y_Py = I_TWOBODYOVERLAP_H5y_S+ABY*I_TWOBODYOVERLAP_G4y_S;
  Double I_TWOBODYOVERLAP_G3yz_Py = I_TWOBODYOVERLAP_H4yz_S+ABY*I_TWOBODYOVERLAP_G3yz_S;
  Double I_TWOBODYOVERLAP_G2y2z_Py = I_TWOBODYOVERLAP_H3y2z_S+ABY*I_TWOBODYOVERLAP_G2y2z_S;
  Double I_TWOBODYOVERLAP_Gy3z_Py = I_TWOBODYOVERLAP_H2y3z_S+ABY*I_TWOBODYOVERLAP_Gy3z_S;
  Double I_TWOBODYOVERLAP_G4z_Py = I_TWOBODYOVERLAP_Hy4z_S+ABY*I_TWOBODYOVERLAP_G4z_S;
  Double I_TWOBODYOVERLAP_G4x_Pz = I_TWOBODYOVERLAP_H4xz_S+ABZ*I_TWOBODYOVERLAP_G4x_S;
  Double I_TWOBODYOVERLAP_G3xy_Pz = I_TWOBODYOVERLAP_H3xyz_S+ABZ*I_TWOBODYOVERLAP_G3xy_S;
  Double I_TWOBODYOVERLAP_G3xz_Pz = I_TWOBODYOVERLAP_H3x2z_S+ABZ*I_TWOBODYOVERLAP_G3xz_S;
  Double I_TWOBODYOVERLAP_G2x2y_Pz = I_TWOBODYOVERLAP_H2x2yz_S+ABZ*I_TWOBODYOVERLAP_G2x2y_S;
  Double I_TWOBODYOVERLAP_G2xyz_Pz = I_TWOBODYOVERLAP_H2xy2z_S+ABZ*I_TWOBODYOVERLAP_G2xyz_S;
  Double I_TWOBODYOVERLAP_G2x2z_Pz = I_TWOBODYOVERLAP_H2x3z_S+ABZ*I_TWOBODYOVERLAP_G2x2z_S;
  Double I_TWOBODYOVERLAP_Gx3y_Pz = I_TWOBODYOVERLAP_Hx3yz_S+ABZ*I_TWOBODYOVERLAP_Gx3y_S;
  Double I_TWOBODYOVERLAP_Gx2yz_Pz = I_TWOBODYOVERLAP_Hx2y2z_S+ABZ*I_TWOBODYOVERLAP_Gx2yz_S;
  Double I_TWOBODYOVERLAP_Gxy2z_Pz = I_TWOBODYOVERLAP_Hxy3z_S+ABZ*I_TWOBODYOVERLAP_Gxy2z_S;
  Double I_TWOBODYOVERLAP_Gx3z_Pz = I_TWOBODYOVERLAP_Hx4z_S+ABZ*I_TWOBODYOVERLAP_Gx3z_S;
  Double I_TWOBODYOVERLAP_G4y_Pz = I_TWOBODYOVERLAP_H4yz_S+ABZ*I_TWOBODYOVERLAP_G4y_S;
  Double I_TWOBODYOVERLAP_G3yz_Pz = I_TWOBODYOVERLAP_H3y2z_S+ABZ*I_TWOBODYOVERLAP_G3yz_S;
  Double I_TWOBODYOVERLAP_G2y2z_Pz = I_TWOBODYOVERLAP_H2y3z_S+ABZ*I_TWOBODYOVERLAP_G2y2z_S;
  Double I_TWOBODYOVERLAP_Gy3z_Pz = I_TWOBODYOVERLAP_Hy4z_S+ABZ*I_TWOBODYOVERLAP_Gy3z_S;
  Double I_TWOBODYOVERLAP_G4z_Pz = I_TWOBODYOVERLAP_H5z_S+ABZ*I_TWOBODYOVERLAP_G4z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_Px = I_TWOBODYOVERLAP_I6x_S+ABX*I_TWOBODYOVERLAP_H5x_S;
  Double I_TWOBODYOVERLAP_H4xy_Px = I_TWOBODYOVERLAP_I5xy_S+ABX*I_TWOBODYOVERLAP_H4xy_S;
  Double I_TWOBODYOVERLAP_H4xz_Px = I_TWOBODYOVERLAP_I5xz_S+ABX*I_TWOBODYOVERLAP_H4xz_S;
  Double I_TWOBODYOVERLAP_H3x2y_Px = I_TWOBODYOVERLAP_I4x2y_S+ABX*I_TWOBODYOVERLAP_H3x2y_S;
  Double I_TWOBODYOVERLAP_H3xyz_Px = I_TWOBODYOVERLAP_I4xyz_S+ABX*I_TWOBODYOVERLAP_H3xyz_S;
  Double I_TWOBODYOVERLAP_H3x2z_Px = I_TWOBODYOVERLAP_I4x2z_S+ABX*I_TWOBODYOVERLAP_H3x2z_S;
  Double I_TWOBODYOVERLAP_H2x3y_Px = I_TWOBODYOVERLAP_I3x3y_S+ABX*I_TWOBODYOVERLAP_H2x3y_S;
  Double I_TWOBODYOVERLAP_H2x2yz_Px = I_TWOBODYOVERLAP_I3x2yz_S+ABX*I_TWOBODYOVERLAP_H2x2yz_S;
  Double I_TWOBODYOVERLAP_H2xy2z_Px = I_TWOBODYOVERLAP_I3xy2z_S+ABX*I_TWOBODYOVERLAP_H2xy2z_S;
  Double I_TWOBODYOVERLAP_H2x3z_Px = I_TWOBODYOVERLAP_I3x3z_S+ABX*I_TWOBODYOVERLAP_H2x3z_S;
  Double I_TWOBODYOVERLAP_Hx4y_Px = I_TWOBODYOVERLAP_I2x4y_S+ABX*I_TWOBODYOVERLAP_Hx4y_S;
  Double I_TWOBODYOVERLAP_Hx3yz_Px = I_TWOBODYOVERLAP_I2x3yz_S+ABX*I_TWOBODYOVERLAP_Hx3yz_S;
  Double I_TWOBODYOVERLAP_Hx2y2z_Px = I_TWOBODYOVERLAP_I2x2y2z_S+ABX*I_TWOBODYOVERLAP_Hx2y2z_S;
  Double I_TWOBODYOVERLAP_Hxy3z_Px = I_TWOBODYOVERLAP_I2xy3z_S+ABX*I_TWOBODYOVERLAP_Hxy3z_S;
  Double I_TWOBODYOVERLAP_Hx4z_Px = I_TWOBODYOVERLAP_I2x4z_S+ABX*I_TWOBODYOVERLAP_Hx4z_S;
  Double I_TWOBODYOVERLAP_H5y_Px = I_TWOBODYOVERLAP_Ix5y_S+ABX*I_TWOBODYOVERLAP_H5y_S;
  Double I_TWOBODYOVERLAP_H4yz_Px = I_TWOBODYOVERLAP_Ix4yz_S+ABX*I_TWOBODYOVERLAP_H4yz_S;
  Double I_TWOBODYOVERLAP_H3y2z_Px = I_TWOBODYOVERLAP_Ix3y2z_S+ABX*I_TWOBODYOVERLAP_H3y2z_S;
  Double I_TWOBODYOVERLAP_H2y3z_Px = I_TWOBODYOVERLAP_Ix2y3z_S+ABX*I_TWOBODYOVERLAP_H2y3z_S;
  Double I_TWOBODYOVERLAP_Hy4z_Px = I_TWOBODYOVERLAP_Ixy4z_S+ABX*I_TWOBODYOVERLAP_Hy4z_S;
  Double I_TWOBODYOVERLAP_H5z_Px = I_TWOBODYOVERLAP_Ix5z_S+ABX*I_TWOBODYOVERLAP_H5z_S;
  Double I_TWOBODYOVERLAP_H5x_Py = I_TWOBODYOVERLAP_I5xy_S+ABY*I_TWOBODYOVERLAP_H5x_S;
  Double I_TWOBODYOVERLAP_H4xy_Py = I_TWOBODYOVERLAP_I4x2y_S+ABY*I_TWOBODYOVERLAP_H4xy_S;
  Double I_TWOBODYOVERLAP_H4xz_Py = I_TWOBODYOVERLAP_I4xyz_S+ABY*I_TWOBODYOVERLAP_H4xz_S;
  Double I_TWOBODYOVERLAP_H3x2y_Py = I_TWOBODYOVERLAP_I3x3y_S+ABY*I_TWOBODYOVERLAP_H3x2y_S;
  Double I_TWOBODYOVERLAP_H3xyz_Py = I_TWOBODYOVERLAP_I3x2yz_S+ABY*I_TWOBODYOVERLAP_H3xyz_S;
  Double I_TWOBODYOVERLAP_H3x2z_Py = I_TWOBODYOVERLAP_I3xy2z_S+ABY*I_TWOBODYOVERLAP_H3x2z_S;
  Double I_TWOBODYOVERLAP_H2x3y_Py = I_TWOBODYOVERLAP_I2x4y_S+ABY*I_TWOBODYOVERLAP_H2x3y_S;
  Double I_TWOBODYOVERLAP_H2x2yz_Py = I_TWOBODYOVERLAP_I2x3yz_S+ABY*I_TWOBODYOVERLAP_H2x2yz_S;
  Double I_TWOBODYOVERLAP_H2xy2z_Py = I_TWOBODYOVERLAP_I2x2y2z_S+ABY*I_TWOBODYOVERLAP_H2xy2z_S;
  Double I_TWOBODYOVERLAP_H2x3z_Py = I_TWOBODYOVERLAP_I2xy3z_S+ABY*I_TWOBODYOVERLAP_H2x3z_S;
  Double I_TWOBODYOVERLAP_Hx4y_Py = I_TWOBODYOVERLAP_Ix5y_S+ABY*I_TWOBODYOVERLAP_Hx4y_S;
  Double I_TWOBODYOVERLAP_Hx3yz_Py = I_TWOBODYOVERLAP_Ix4yz_S+ABY*I_TWOBODYOVERLAP_Hx3yz_S;
  Double I_TWOBODYOVERLAP_Hx2y2z_Py = I_TWOBODYOVERLAP_Ix3y2z_S+ABY*I_TWOBODYOVERLAP_Hx2y2z_S;
  Double I_TWOBODYOVERLAP_Hxy3z_Py = I_TWOBODYOVERLAP_Ix2y3z_S+ABY*I_TWOBODYOVERLAP_Hxy3z_S;
  Double I_TWOBODYOVERLAP_Hx4z_Py = I_TWOBODYOVERLAP_Ixy4z_S+ABY*I_TWOBODYOVERLAP_Hx4z_S;
  Double I_TWOBODYOVERLAP_H5y_Py = I_TWOBODYOVERLAP_I6y_S+ABY*I_TWOBODYOVERLAP_H5y_S;
  Double I_TWOBODYOVERLAP_H4yz_Py = I_TWOBODYOVERLAP_I5yz_S+ABY*I_TWOBODYOVERLAP_H4yz_S;
  Double I_TWOBODYOVERLAP_H3y2z_Py = I_TWOBODYOVERLAP_I4y2z_S+ABY*I_TWOBODYOVERLAP_H3y2z_S;
  Double I_TWOBODYOVERLAP_H2y3z_Py = I_TWOBODYOVERLAP_I3y3z_S+ABY*I_TWOBODYOVERLAP_H2y3z_S;
  Double I_TWOBODYOVERLAP_Hy4z_Py = I_TWOBODYOVERLAP_I2y4z_S+ABY*I_TWOBODYOVERLAP_Hy4z_S;
  Double I_TWOBODYOVERLAP_H5z_Py = I_TWOBODYOVERLAP_Iy5z_S+ABY*I_TWOBODYOVERLAP_H5z_S;
  Double I_TWOBODYOVERLAP_H5x_Pz = I_TWOBODYOVERLAP_I5xz_S+ABZ*I_TWOBODYOVERLAP_H5x_S;
  Double I_TWOBODYOVERLAP_H4xy_Pz = I_TWOBODYOVERLAP_I4xyz_S+ABZ*I_TWOBODYOVERLAP_H4xy_S;
  Double I_TWOBODYOVERLAP_H4xz_Pz = I_TWOBODYOVERLAP_I4x2z_S+ABZ*I_TWOBODYOVERLAP_H4xz_S;
  Double I_TWOBODYOVERLAP_H3x2y_Pz = I_TWOBODYOVERLAP_I3x2yz_S+ABZ*I_TWOBODYOVERLAP_H3x2y_S;
  Double I_TWOBODYOVERLAP_H3xyz_Pz = I_TWOBODYOVERLAP_I3xy2z_S+ABZ*I_TWOBODYOVERLAP_H3xyz_S;
  Double I_TWOBODYOVERLAP_H3x2z_Pz = I_TWOBODYOVERLAP_I3x3z_S+ABZ*I_TWOBODYOVERLAP_H3x2z_S;
  Double I_TWOBODYOVERLAP_H2x3y_Pz = I_TWOBODYOVERLAP_I2x3yz_S+ABZ*I_TWOBODYOVERLAP_H2x3y_S;
  Double I_TWOBODYOVERLAP_H2x2yz_Pz = I_TWOBODYOVERLAP_I2x2y2z_S+ABZ*I_TWOBODYOVERLAP_H2x2yz_S;
  Double I_TWOBODYOVERLAP_H2xy2z_Pz = I_TWOBODYOVERLAP_I2xy3z_S+ABZ*I_TWOBODYOVERLAP_H2xy2z_S;
  Double I_TWOBODYOVERLAP_H2x3z_Pz = I_TWOBODYOVERLAP_I2x4z_S+ABZ*I_TWOBODYOVERLAP_H2x3z_S;
  Double I_TWOBODYOVERLAP_Hx4y_Pz = I_TWOBODYOVERLAP_Ix4yz_S+ABZ*I_TWOBODYOVERLAP_Hx4y_S;
  Double I_TWOBODYOVERLAP_Hx3yz_Pz = I_TWOBODYOVERLAP_Ix3y2z_S+ABZ*I_TWOBODYOVERLAP_Hx3yz_S;
  Double I_TWOBODYOVERLAP_Hx2y2z_Pz = I_TWOBODYOVERLAP_Ix2y3z_S+ABZ*I_TWOBODYOVERLAP_Hx2y2z_S;
  Double I_TWOBODYOVERLAP_Hxy3z_Pz = I_TWOBODYOVERLAP_Ixy4z_S+ABZ*I_TWOBODYOVERLAP_Hxy3z_S;
  Double I_TWOBODYOVERLAP_Hx4z_Pz = I_TWOBODYOVERLAP_Ix5z_S+ABZ*I_TWOBODYOVERLAP_Hx4z_S;
  Double I_TWOBODYOVERLAP_H5y_Pz = I_TWOBODYOVERLAP_I5yz_S+ABZ*I_TWOBODYOVERLAP_H5y_S;
  Double I_TWOBODYOVERLAP_H4yz_Pz = I_TWOBODYOVERLAP_I4y2z_S+ABZ*I_TWOBODYOVERLAP_H4yz_S;
  Double I_TWOBODYOVERLAP_H3y2z_Pz = I_TWOBODYOVERLAP_I3y3z_S+ABZ*I_TWOBODYOVERLAP_H3y2z_S;
  Double I_TWOBODYOVERLAP_H2y3z_Pz = I_TWOBODYOVERLAP_I2y4z_S+ABZ*I_TWOBODYOVERLAP_H2y3z_S;
  Double I_TWOBODYOVERLAP_Hy4z_Pz = I_TWOBODYOVERLAP_Iy5z_S+ABZ*I_TWOBODYOVERLAP_Hy4z_S;
  Double I_TWOBODYOVERLAP_H5z_Pz = I_TWOBODYOVERLAP_I6z_S+ABZ*I_TWOBODYOVERLAP_H5z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_D2x = I_TWOBODYOVERLAP_H5x_Px+ABX*I_TWOBODYOVERLAP_G4x_Px;
  Double I_TWOBODYOVERLAP_G3xy_D2x = I_TWOBODYOVERLAP_H4xy_Px+ABX*I_TWOBODYOVERLAP_G3xy_Px;
  Double I_TWOBODYOVERLAP_G3xz_D2x = I_TWOBODYOVERLAP_H4xz_Px+ABX*I_TWOBODYOVERLAP_G3xz_Px;
  Double I_TWOBODYOVERLAP_G2x2y_D2x = I_TWOBODYOVERLAP_H3x2y_Px+ABX*I_TWOBODYOVERLAP_G2x2y_Px;
  Double I_TWOBODYOVERLAP_G2xyz_D2x = I_TWOBODYOVERLAP_H3xyz_Px+ABX*I_TWOBODYOVERLAP_G2xyz_Px;
  Double I_TWOBODYOVERLAP_G2x2z_D2x = I_TWOBODYOVERLAP_H3x2z_Px+ABX*I_TWOBODYOVERLAP_G2x2z_Px;
  Double I_TWOBODYOVERLAP_Gx3y_D2x = I_TWOBODYOVERLAP_H2x3y_Px+ABX*I_TWOBODYOVERLAP_Gx3y_Px;
  Double I_TWOBODYOVERLAP_Gx2yz_D2x = I_TWOBODYOVERLAP_H2x2yz_Px+ABX*I_TWOBODYOVERLAP_Gx2yz_Px;
  Double I_TWOBODYOVERLAP_Gxy2z_D2x = I_TWOBODYOVERLAP_H2xy2z_Px+ABX*I_TWOBODYOVERLAP_Gxy2z_Px;
  Double I_TWOBODYOVERLAP_Gx3z_D2x = I_TWOBODYOVERLAP_H2x3z_Px+ABX*I_TWOBODYOVERLAP_Gx3z_Px;
  Double I_TWOBODYOVERLAP_G4y_D2x = I_TWOBODYOVERLAP_Hx4y_Px+ABX*I_TWOBODYOVERLAP_G4y_Px;
  Double I_TWOBODYOVERLAP_G3yz_D2x = I_TWOBODYOVERLAP_Hx3yz_Px+ABX*I_TWOBODYOVERLAP_G3yz_Px;
  Double I_TWOBODYOVERLAP_G2y2z_D2x = I_TWOBODYOVERLAP_Hx2y2z_Px+ABX*I_TWOBODYOVERLAP_G2y2z_Px;
  Double I_TWOBODYOVERLAP_Gy3z_D2x = I_TWOBODYOVERLAP_Hxy3z_Px+ABX*I_TWOBODYOVERLAP_Gy3z_Px;
  Double I_TWOBODYOVERLAP_G4z_D2x = I_TWOBODYOVERLAP_Hx4z_Px+ABX*I_TWOBODYOVERLAP_G4z_Px;
  Double I_TWOBODYOVERLAP_G4x_Dxy = I_TWOBODYOVERLAP_H4xy_Px+ABY*I_TWOBODYOVERLAP_G4x_Px;
  Double I_TWOBODYOVERLAP_G3xy_Dxy = I_TWOBODYOVERLAP_H3x2y_Px+ABY*I_TWOBODYOVERLAP_G3xy_Px;
  Double I_TWOBODYOVERLAP_G3xz_Dxy = I_TWOBODYOVERLAP_H3xyz_Px+ABY*I_TWOBODYOVERLAP_G3xz_Px;
  Double I_TWOBODYOVERLAP_G2x2y_Dxy = I_TWOBODYOVERLAP_H2x3y_Px+ABY*I_TWOBODYOVERLAP_G2x2y_Px;
  Double I_TWOBODYOVERLAP_G2xyz_Dxy = I_TWOBODYOVERLAP_H2x2yz_Px+ABY*I_TWOBODYOVERLAP_G2xyz_Px;
  Double I_TWOBODYOVERLAP_G2x2z_Dxy = I_TWOBODYOVERLAP_H2xy2z_Px+ABY*I_TWOBODYOVERLAP_G2x2z_Px;
  Double I_TWOBODYOVERLAP_Gx3y_Dxy = I_TWOBODYOVERLAP_Hx4y_Px+ABY*I_TWOBODYOVERLAP_Gx3y_Px;
  Double I_TWOBODYOVERLAP_Gx2yz_Dxy = I_TWOBODYOVERLAP_Hx3yz_Px+ABY*I_TWOBODYOVERLAP_Gx2yz_Px;
  Double I_TWOBODYOVERLAP_Gxy2z_Dxy = I_TWOBODYOVERLAP_Hx2y2z_Px+ABY*I_TWOBODYOVERLAP_Gxy2z_Px;
  Double I_TWOBODYOVERLAP_Gx3z_Dxy = I_TWOBODYOVERLAP_Hxy3z_Px+ABY*I_TWOBODYOVERLAP_Gx3z_Px;
  Double I_TWOBODYOVERLAP_G4y_Dxy = I_TWOBODYOVERLAP_H5y_Px+ABY*I_TWOBODYOVERLAP_G4y_Px;
  Double I_TWOBODYOVERLAP_G3yz_Dxy = I_TWOBODYOVERLAP_H4yz_Px+ABY*I_TWOBODYOVERLAP_G3yz_Px;
  Double I_TWOBODYOVERLAP_G2y2z_Dxy = I_TWOBODYOVERLAP_H3y2z_Px+ABY*I_TWOBODYOVERLAP_G2y2z_Px;
  Double I_TWOBODYOVERLAP_Gy3z_Dxy = I_TWOBODYOVERLAP_H2y3z_Px+ABY*I_TWOBODYOVERLAP_Gy3z_Px;
  Double I_TWOBODYOVERLAP_G4z_Dxy = I_TWOBODYOVERLAP_Hy4z_Px+ABY*I_TWOBODYOVERLAP_G4z_Px;
  Double I_TWOBODYOVERLAP_G4x_D2y = I_TWOBODYOVERLAP_H4xy_Py+ABY*I_TWOBODYOVERLAP_G4x_Py;
  Double I_TWOBODYOVERLAP_G3xy_D2y = I_TWOBODYOVERLAP_H3x2y_Py+ABY*I_TWOBODYOVERLAP_G3xy_Py;
  Double I_TWOBODYOVERLAP_G3xz_D2y = I_TWOBODYOVERLAP_H3xyz_Py+ABY*I_TWOBODYOVERLAP_G3xz_Py;
  Double I_TWOBODYOVERLAP_G2x2y_D2y = I_TWOBODYOVERLAP_H2x3y_Py+ABY*I_TWOBODYOVERLAP_G2x2y_Py;
  Double I_TWOBODYOVERLAP_G2xyz_D2y = I_TWOBODYOVERLAP_H2x2yz_Py+ABY*I_TWOBODYOVERLAP_G2xyz_Py;
  Double I_TWOBODYOVERLAP_G2x2z_D2y = I_TWOBODYOVERLAP_H2xy2z_Py+ABY*I_TWOBODYOVERLAP_G2x2z_Py;
  Double I_TWOBODYOVERLAP_Gx3y_D2y = I_TWOBODYOVERLAP_Hx4y_Py+ABY*I_TWOBODYOVERLAP_Gx3y_Py;
  Double I_TWOBODYOVERLAP_Gx2yz_D2y = I_TWOBODYOVERLAP_Hx3yz_Py+ABY*I_TWOBODYOVERLAP_Gx2yz_Py;
  Double I_TWOBODYOVERLAP_Gxy2z_D2y = I_TWOBODYOVERLAP_Hx2y2z_Py+ABY*I_TWOBODYOVERLAP_Gxy2z_Py;
  Double I_TWOBODYOVERLAP_Gx3z_D2y = I_TWOBODYOVERLAP_Hxy3z_Py+ABY*I_TWOBODYOVERLAP_Gx3z_Py;
  Double I_TWOBODYOVERLAP_G4y_D2y = I_TWOBODYOVERLAP_H5y_Py+ABY*I_TWOBODYOVERLAP_G4y_Py;
  Double I_TWOBODYOVERLAP_G3yz_D2y = I_TWOBODYOVERLAP_H4yz_Py+ABY*I_TWOBODYOVERLAP_G3yz_Py;
  Double I_TWOBODYOVERLAP_G2y2z_D2y = I_TWOBODYOVERLAP_H3y2z_Py+ABY*I_TWOBODYOVERLAP_G2y2z_Py;
  Double I_TWOBODYOVERLAP_Gy3z_D2y = I_TWOBODYOVERLAP_H2y3z_Py+ABY*I_TWOBODYOVERLAP_Gy3z_Py;
  Double I_TWOBODYOVERLAP_G4z_D2y = I_TWOBODYOVERLAP_Hy4z_Py+ABY*I_TWOBODYOVERLAP_G4z_Py;
  Double I_TWOBODYOVERLAP_G4x_D2z = I_TWOBODYOVERLAP_H4xz_Pz+ABZ*I_TWOBODYOVERLAP_G4x_Pz;
  Double I_TWOBODYOVERLAP_G3xy_D2z = I_TWOBODYOVERLAP_H3xyz_Pz+ABZ*I_TWOBODYOVERLAP_G3xy_Pz;
  Double I_TWOBODYOVERLAP_G3xz_D2z = I_TWOBODYOVERLAP_H3x2z_Pz+ABZ*I_TWOBODYOVERLAP_G3xz_Pz;
  Double I_TWOBODYOVERLAP_G2x2y_D2z = I_TWOBODYOVERLAP_H2x2yz_Pz+ABZ*I_TWOBODYOVERLAP_G2x2y_Pz;
  Double I_TWOBODYOVERLAP_G2xyz_D2z = I_TWOBODYOVERLAP_H2xy2z_Pz+ABZ*I_TWOBODYOVERLAP_G2xyz_Pz;
  Double I_TWOBODYOVERLAP_G2x2z_D2z = I_TWOBODYOVERLAP_H2x3z_Pz+ABZ*I_TWOBODYOVERLAP_G2x2z_Pz;
  Double I_TWOBODYOVERLAP_Gx3y_D2z = I_TWOBODYOVERLAP_Hx3yz_Pz+ABZ*I_TWOBODYOVERLAP_Gx3y_Pz;
  Double I_TWOBODYOVERLAP_Gx2yz_D2z = I_TWOBODYOVERLAP_Hx2y2z_Pz+ABZ*I_TWOBODYOVERLAP_Gx2yz_Pz;
  Double I_TWOBODYOVERLAP_Gxy2z_D2z = I_TWOBODYOVERLAP_Hxy3z_Pz+ABZ*I_TWOBODYOVERLAP_Gxy2z_Pz;
  Double I_TWOBODYOVERLAP_Gx3z_D2z = I_TWOBODYOVERLAP_Hx4z_Pz+ABZ*I_TWOBODYOVERLAP_Gx3z_Pz;
  Double I_TWOBODYOVERLAP_G4y_D2z = I_TWOBODYOVERLAP_H4yz_Pz+ABZ*I_TWOBODYOVERLAP_G4y_Pz;
  Double I_TWOBODYOVERLAP_G3yz_D2z = I_TWOBODYOVERLAP_H3y2z_Pz+ABZ*I_TWOBODYOVERLAP_G3yz_Pz;
  Double I_TWOBODYOVERLAP_G2y2z_D2z = I_TWOBODYOVERLAP_H2y3z_Pz+ABZ*I_TWOBODYOVERLAP_G2y2z_Pz;
  Double I_TWOBODYOVERLAP_Gy3z_D2z = I_TWOBODYOVERLAP_Hy4z_Pz+ABZ*I_TWOBODYOVERLAP_Gy3z_Pz;
  Double I_TWOBODYOVERLAP_G4z_D2z = I_TWOBODYOVERLAP_H5z_Pz+ABZ*I_TWOBODYOVERLAP_G4z_Pz;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 16 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_Px = I_TWOBODYOVERLAP_K7x_S+ABX*I_TWOBODYOVERLAP_I6x_S;
  Double I_TWOBODYOVERLAP_I5xy_Px = I_TWOBODYOVERLAP_K6xy_S+ABX*I_TWOBODYOVERLAP_I5xy_S;
  Double I_TWOBODYOVERLAP_I5xz_Px = I_TWOBODYOVERLAP_K6xz_S+ABX*I_TWOBODYOVERLAP_I5xz_S;
  Double I_TWOBODYOVERLAP_I4x2y_Px = I_TWOBODYOVERLAP_K5x2y_S+ABX*I_TWOBODYOVERLAP_I4x2y_S;
  Double I_TWOBODYOVERLAP_I4xyz_Px = I_TWOBODYOVERLAP_K5xyz_S+ABX*I_TWOBODYOVERLAP_I4xyz_S;
  Double I_TWOBODYOVERLAP_I4x2z_Px = I_TWOBODYOVERLAP_K5x2z_S+ABX*I_TWOBODYOVERLAP_I4x2z_S;
  Double I_TWOBODYOVERLAP_I3x3y_Px = I_TWOBODYOVERLAP_K4x3y_S+ABX*I_TWOBODYOVERLAP_I3x3y_S;
  Double I_TWOBODYOVERLAP_I3x2yz_Px = I_TWOBODYOVERLAP_K4x2yz_S+ABX*I_TWOBODYOVERLAP_I3x2yz_S;
  Double I_TWOBODYOVERLAP_I3xy2z_Px = I_TWOBODYOVERLAP_K4xy2z_S+ABX*I_TWOBODYOVERLAP_I3xy2z_S;
  Double I_TWOBODYOVERLAP_I3x3z_Px = I_TWOBODYOVERLAP_K4x3z_S+ABX*I_TWOBODYOVERLAP_I3x3z_S;
  Double I_TWOBODYOVERLAP_I2x4y_Px = I_TWOBODYOVERLAP_K3x4y_S+ABX*I_TWOBODYOVERLAP_I2x4y_S;
  Double I_TWOBODYOVERLAP_I2x3yz_Px = I_TWOBODYOVERLAP_K3x3yz_S+ABX*I_TWOBODYOVERLAP_I2x3yz_S;
  Double I_TWOBODYOVERLAP_I2x2y2z_Px = I_TWOBODYOVERLAP_K3x2y2z_S+ABX*I_TWOBODYOVERLAP_I2x2y2z_S;
  Double I_TWOBODYOVERLAP_I2xy3z_Px = I_TWOBODYOVERLAP_K3xy3z_S+ABX*I_TWOBODYOVERLAP_I2xy3z_S;
  Double I_TWOBODYOVERLAP_I2x4z_Px = I_TWOBODYOVERLAP_K3x4z_S+ABX*I_TWOBODYOVERLAP_I2x4z_S;
  Double I_TWOBODYOVERLAP_Ix5y_Px = I_TWOBODYOVERLAP_K2x5y_S+ABX*I_TWOBODYOVERLAP_Ix5y_S;
  Double I_TWOBODYOVERLAP_Ix4yz_Px = I_TWOBODYOVERLAP_K2x4yz_S+ABX*I_TWOBODYOVERLAP_Ix4yz_S;
  Double I_TWOBODYOVERLAP_Ix3y2z_Px = I_TWOBODYOVERLAP_K2x3y2z_S+ABX*I_TWOBODYOVERLAP_Ix3y2z_S;
  Double I_TWOBODYOVERLAP_Ix2y3z_Px = I_TWOBODYOVERLAP_K2x2y3z_S+ABX*I_TWOBODYOVERLAP_Ix2y3z_S;
  Double I_TWOBODYOVERLAP_Ixy4z_Px = I_TWOBODYOVERLAP_K2xy4z_S+ABX*I_TWOBODYOVERLAP_Ixy4z_S;
  Double I_TWOBODYOVERLAP_Ix5z_Px = I_TWOBODYOVERLAP_K2x5z_S+ABX*I_TWOBODYOVERLAP_Ix5z_S;
  Double I_TWOBODYOVERLAP_I5yz_Px = I_TWOBODYOVERLAP_Kx5yz_S+ABX*I_TWOBODYOVERLAP_I5yz_S;
  Double I_TWOBODYOVERLAP_I4y2z_Px = I_TWOBODYOVERLAP_Kx4y2z_S+ABX*I_TWOBODYOVERLAP_I4y2z_S;
  Double I_TWOBODYOVERLAP_I3y3z_Px = I_TWOBODYOVERLAP_Kx3y3z_S+ABX*I_TWOBODYOVERLAP_I3y3z_S;
  Double I_TWOBODYOVERLAP_I2y4z_Px = I_TWOBODYOVERLAP_Kx2y4z_S+ABX*I_TWOBODYOVERLAP_I2y4z_S;
  Double I_TWOBODYOVERLAP_Iy5z_Px = I_TWOBODYOVERLAP_Kxy5z_S+ABX*I_TWOBODYOVERLAP_Iy5z_S;
  Double I_TWOBODYOVERLAP_I5xy_Py = I_TWOBODYOVERLAP_K5x2y_S+ABY*I_TWOBODYOVERLAP_I5xy_S;
  Double I_TWOBODYOVERLAP_I4x2y_Py = I_TWOBODYOVERLAP_K4x3y_S+ABY*I_TWOBODYOVERLAP_I4x2y_S;
  Double I_TWOBODYOVERLAP_I4xyz_Py = I_TWOBODYOVERLAP_K4x2yz_S+ABY*I_TWOBODYOVERLAP_I4xyz_S;
  Double I_TWOBODYOVERLAP_I3x3y_Py = I_TWOBODYOVERLAP_K3x4y_S+ABY*I_TWOBODYOVERLAP_I3x3y_S;
  Double I_TWOBODYOVERLAP_I3x2yz_Py = I_TWOBODYOVERLAP_K3x3yz_S+ABY*I_TWOBODYOVERLAP_I3x2yz_S;
  Double I_TWOBODYOVERLAP_I3xy2z_Py = I_TWOBODYOVERLAP_K3x2y2z_S+ABY*I_TWOBODYOVERLAP_I3xy2z_S;
  Double I_TWOBODYOVERLAP_I2x4y_Py = I_TWOBODYOVERLAP_K2x5y_S+ABY*I_TWOBODYOVERLAP_I2x4y_S;
  Double I_TWOBODYOVERLAP_I2x3yz_Py = I_TWOBODYOVERLAP_K2x4yz_S+ABY*I_TWOBODYOVERLAP_I2x3yz_S;
  Double I_TWOBODYOVERLAP_I2x2y2z_Py = I_TWOBODYOVERLAP_K2x3y2z_S+ABY*I_TWOBODYOVERLAP_I2x2y2z_S;
  Double I_TWOBODYOVERLAP_I2xy3z_Py = I_TWOBODYOVERLAP_K2x2y3z_S+ABY*I_TWOBODYOVERLAP_I2xy3z_S;
  Double I_TWOBODYOVERLAP_Ix5y_Py = I_TWOBODYOVERLAP_Kx6y_S+ABY*I_TWOBODYOVERLAP_Ix5y_S;
  Double I_TWOBODYOVERLAP_Ix4yz_Py = I_TWOBODYOVERLAP_Kx5yz_S+ABY*I_TWOBODYOVERLAP_Ix4yz_S;
  Double I_TWOBODYOVERLAP_Ix3y2z_Py = I_TWOBODYOVERLAP_Kx4y2z_S+ABY*I_TWOBODYOVERLAP_Ix3y2z_S;
  Double I_TWOBODYOVERLAP_Ix2y3z_Py = I_TWOBODYOVERLAP_Kx3y3z_S+ABY*I_TWOBODYOVERLAP_Ix2y3z_S;
  Double I_TWOBODYOVERLAP_Ixy4z_Py = I_TWOBODYOVERLAP_Kx2y4z_S+ABY*I_TWOBODYOVERLAP_Ixy4z_S;
  Double I_TWOBODYOVERLAP_I6y_Py = I_TWOBODYOVERLAP_K7y_S+ABY*I_TWOBODYOVERLAP_I6y_S;
  Double I_TWOBODYOVERLAP_I5yz_Py = I_TWOBODYOVERLAP_K6yz_S+ABY*I_TWOBODYOVERLAP_I5yz_S;
  Double I_TWOBODYOVERLAP_I4y2z_Py = I_TWOBODYOVERLAP_K5y2z_S+ABY*I_TWOBODYOVERLAP_I4y2z_S;
  Double I_TWOBODYOVERLAP_I3y3z_Py = I_TWOBODYOVERLAP_K4y3z_S+ABY*I_TWOBODYOVERLAP_I3y3z_S;
  Double I_TWOBODYOVERLAP_I2y4z_Py = I_TWOBODYOVERLAP_K3y4z_S+ABY*I_TWOBODYOVERLAP_I2y4z_S;
  Double I_TWOBODYOVERLAP_Iy5z_Py = I_TWOBODYOVERLAP_K2y5z_S+ABY*I_TWOBODYOVERLAP_Iy5z_S;
  Double I_TWOBODYOVERLAP_I5xz_Pz = I_TWOBODYOVERLAP_K5x2z_S+ABZ*I_TWOBODYOVERLAP_I5xz_S;
  Double I_TWOBODYOVERLAP_I4xyz_Pz = I_TWOBODYOVERLAP_K4xy2z_S+ABZ*I_TWOBODYOVERLAP_I4xyz_S;
  Double I_TWOBODYOVERLAP_I4x2z_Pz = I_TWOBODYOVERLAP_K4x3z_S+ABZ*I_TWOBODYOVERLAP_I4x2z_S;
  Double I_TWOBODYOVERLAP_I3x2yz_Pz = I_TWOBODYOVERLAP_K3x2y2z_S+ABZ*I_TWOBODYOVERLAP_I3x2yz_S;
  Double I_TWOBODYOVERLAP_I3xy2z_Pz = I_TWOBODYOVERLAP_K3xy3z_S+ABZ*I_TWOBODYOVERLAP_I3xy2z_S;
  Double I_TWOBODYOVERLAP_I3x3z_Pz = I_TWOBODYOVERLAP_K3x4z_S+ABZ*I_TWOBODYOVERLAP_I3x3z_S;
  Double I_TWOBODYOVERLAP_I2x3yz_Pz = I_TWOBODYOVERLAP_K2x3y2z_S+ABZ*I_TWOBODYOVERLAP_I2x3yz_S;
  Double I_TWOBODYOVERLAP_I2x2y2z_Pz = I_TWOBODYOVERLAP_K2x2y3z_S+ABZ*I_TWOBODYOVERLAP_I2x2y2z_S;
  Double I_TWOBODYOVERLAP_I2xy3z_Pz = I_TWOBODYOVERLAP_K2xy4z_S+ABZ*I_TWOBODYOVERLAP_I2xy3z_S;
  Double I_TWOBODYOVERLAP_I2x4z_Pz = I_TWOBODYOVERLAP_K2x5z_S+ABZ*I_TWOBODYOVERLAP_I2x4z_S;
  Double I_TWOBODYOVERLAP_Ix4yz_Pz = I_TWOBODYOVERLAP_Kx4y2z_S+ABZ*I_TWOBODYOVERLAP_Ix4yz_S;
  Double I_TWOBODYOVERLAP_Ix3y2z_Pz = I_TWOBODYOVERLAP_Kx3y3z_S+ABZ*I_TWOBODYOVERLAP_Ix3y2z_S;
  Double I_TWOBODYOVERLAP_Ix2y3z_Pz = I_TWOBODYOVERLAP_Kx2y4z_S+ABZ*I_TWOBODYOVERLAP_Ix2y3z_S;
  Double I_TWOBODYOVERLAP_Ixy4z_Pz = I_TWOBODYOVERLAP_Kxy5z_S+ABZ*I_TWOBODYOVERLAP_Ixy4z_S;
  Double I_TWOBODYOVERLAP_Ix5z_Pz = I_TWOBODYOVERLAP_Kx6z_S+ABZ*I_TWOBODYOVERLAP_Ix5z_S;
  Double I_TWOBODYOVERLAP_I5yz_Pz = I_TWOBODYOVERLAP_K5y2z_S+ABZ*I_TWOBODYOVERLAP_I5yz_S;
  Double I_TWOBODYOVERLAP_I4y2z_Pz = I_TWOBODYOVERLAP_K4y3z_S+ABZ*I_TWOBODYOVERLAP_I4y2z_S;
  Double I_TWOBODYOVERLAP_I3y3z_Pz = I_TWOBODYOVERLAP_K3y4z_S+ABZ*I_TWOBODYOVERLAP_I3y3z_S;
  Double I_TWOBODYOVERLAP_I2y4z_Pz = I_TWOBODYOVERLAP_K2y5z_S+ABZ*I_TWOBODYOVERLAP_I2y4z_S;
  Double I_TWOBODYOVERLAP_Iy5z_Pz = I_TWOBODYOVERLAP_Ky6z_S+ABZ*I_TWOBODYOVERLAP_Iy5z_S;
  Double I_TWOBODYOVERLAP_I6z_Pz = I_TWOBODYOVERLAP_K7z_S+ABZ*I_TWOBODYOVERLAP_I6z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 48 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_D2x = I_TWOBODYOVERLAP_I6x_Px+ABX*I_TWOBODYOVERLAP_H5x_Px;
  Double I_TWOBODYOVERLAP_H4xy_D2x = I_TWOBODYOVERLAP_I5xy_Px+ABX*I_TWOBODYOVERLAP_H4xy_Px;
  Double I_TWOBODYOVERLAP_H4xz_D2x = I_TWOBODYOVERLAP_I5xz_Px+ABX*I_TWOBODYOVERLAP_H4xz_Px;
  Double I_TWOBODYOVERLAP_H3x2y_D2x = I_TWOBODYOVERLAP_I4x2y_Px+ABX*I_TWOBODYOVERLAP_H3x2y_Px;
  Double I_TWOBODYOVERLAP_H3xyz_D2x = I_TWOBODYOVERLAP_I4xyz_Px+ABX*I_TWOBODYOVERLAP_H3xyz_Px;
  Double I_TWOBODYOVERLAP_H3x2z_D2x = I_TWOBODYOVERLAP_I4x2z_Px+ABX*I_TWOBODYOVERLAP_H3x2z_Px;
  Double I_TWOBODYOVERLAP_H2x3y_D2x = I_TWOBODYOVERLAP_I3x3y_Px+ABX*I_TWOBODYOVERLAP_H2x3y_Px;
  Double I_TWOBODYOVERLAP_H2x2yz_D2x = I_TWOBODYOVERLAP_I3x2yz_Px+ABX*I_TWOBODYOVERLAP_H2x2yz_Px;
  Double I_TWOBODYOVERLAP_H2xy2z_D2x = I_TWOBODYOVERLAP_I3xy2z_Px+ABX*I_TWOBODYOVERLAP_H2xy2z_Px;
  Double I_TWOBODYOVERLAP_H2x3z_D2x = I_TWOBODYOVERLAP_I3x3z_Px+ABX*I_TWOBODYOVERLAP_H2x3z_Px;
  Double I_TWOBODYOVERLAP_Hx4y_D2x = I_TWOBODYOVERLAP_I2x4y_Px+ABX*I_TWOBODYOVERLAP_Hx4y_Px;
  Double I_TWOBODYOVERLAP_Hx3yz_D2x = I_TWOBODYOVERLAP_I2x3yz_Px+ABX*I_TWOBODYOVERLAP_Hx3yz_Px;
  Double I_TWOBODYOVERLAP_Hx2y2z_D2x = I_TWOBODYOVERLAP_I2x2y2z_Px+ABX*I_TWOBODYOVERLAP_Hx2y2z_Px;
  Double I_TWOBODYOVERLAP_Hxy3z_D2x = I_TWOBODYOVERLAP_I2xy3z_Px+ABX*I_TWOBODYOVERLAP_Hxy3z_Px;
  Double I_TWOBODYOVERLAP_Hx4z_D2x = I_TWOBODYOVERLAP_I2x4z_Px+ABX*I_TWOBODYOVERLAP_Hx4z_Px;
  Double I_TWOBODYOVERLAP_H5y_D2x = I_TWOBODYOVERLAP_Ix5y_Px+ABX*I_TWOBODYOVERLAP_H5y_Px;
  Double I_TWOBODYOVERLAP_H4yz_D2x = I_TWOBODYOVERLAP_Ix4yz_Px+ABX*I_TWOBODYOVERLAP_H4yz_Px;
  Double I_TWOBODYOVERLAP_H3y2z_D2x = I_TWOBODYOVERLAP_Ix3y2z_Px+ABX*I_TWOBODYOVERLAP_H3y2z_Px;
  Double I_TWOBODYOVERLAP_H2y3z_D2x = I_TWOBODYOVERLAP_Ix2y3z_Px+ABX*I_TWOBODYOVERLAP_H2y3z_Px;
  Double I_TWOBODYOVERLAP_Hy4z_D2x = I_TWOBODYOVERLAP_Ixy4z_Px+ABX*I_TWOBODYOVERLAP_Hy4z_Px;
  Double I_TWOBODYOVERLAP_H5z_D2x = I_TWOBODYOVERLAP_Ix5z_Px+ABX*I_TWOBODYOVERLAP_H5z_Px;
  Double I_TWOBODYOVERLAP_H4xz_Dxy = I_TWOBODYOVERLAP_I4xyz_Px+ABY*I_TWOBODYOVERLAP_H4xz_Px;
  Double I_TWOBODYOVERLAP_H3xyz_Dxy = I_TWOBODYOVERLAP_I3x2yz_Px+ABY*I_TWOBODYOVERLAP_H3xyz_Px;
  Double I_TWOBODYOVERLAP_H3x2z_Dxy = I_TWOBODYOVERLAP_I3xy2z_Px+ABY*I_TWOBODYOVERLAP_H3x2z_Px;
  Double I_TWOBODYOVERLAP_H2x2yz_Dxy = I_TWOBODYOVERLAP_I2x3yz_Px+ABY*I_TWOBODYOVERLAP_H2x2yz_Px;
  Double I_TWOBODYOVERLAP_H2xy2z_Dxy = I_TWOBODYOVERLAP_I2x2y2z_Px+ABY*I_TWOBODYOVERLAP_H2xy2z_Px;
  Double I_TWOBODYOVERLAP_H2x3z_Dxy = I_TWOBODYOVERLAP_I2xy3z_Px+ABY*I_TWOBODYOVERLAP_H2x3z_Px;
  Double I_TWOBODYOVERLAP_Hx3yz_Dxy = I_TWOBODYOVERLAP_Ix4yz_Px+ABY*I_TWOBODYOVERLAP_Hx3yz_Px;
  Double I_TWOBODYOVERLAP_Hx2y2z_Dxy = I_TWOBODYOVERLAP_Ix3y2z_Px+ABY*I_TWOBODYOVERLAP_Hx2y2z_Px;
  Double I_TWOBODYOVERLAP_Hxy3z_Dxy = I_TWOBODYOVERLAP_Ix2y3z_Px+ABY*I_TWOBODYOVERLAP_Hxy3z_Px;
  Double I_TWOBODYOVERLAP_Hx4z_Dxy = I_TWOBODYOVERLAP_Ixy4z_Px+ABY*I_TWOBODYOVERLAP_Hx4z_Px;
  Double I_TWOBODYOVERLAP_H4yz_Dxy = I_TWOBODYOVERLAP_I5yz_Px+ABY*I_TWOBODYOVERLAP_H4yz_Px;
  Double I_TWOBODYOVERLAP_H3y2z_Dxy = I_TWOBODYOVERLAP_I4y2z_Px+ABY*I_TWOBODYOVERLAP_H3y2z_Px;
  Double I_TWOBODYOVERLAP_H2y3z_Dxy = I_TWOBODYOVERLAP_I3y3z_Px+ABY*I_TWOBODYOVERLAP_H2y3z_Px;
  Double I_TWOBODYOVERLAP_Hy4z_Dxy = I_TWOBODYOVERLAP_I2y4z_Px+ABY*I_TWOBODYOVERLAP_Hy4z_Px;
  Double I_TWOBODYOVERLAP_H5z_Dxy = I_TWOBODYOVERLAP_Iy5z_Px+ABY*I_TWOBODYOVERLAP_H5z_Px;
  Double I_TWOBODYOVERLAP_H5x_D2y = I_TWOBODYOVERLAP_I5xy_Py+ABY*I_TWOBODYOVERLAP_H5x_Py;
  Double I_TWOBODYOVERLAP_H4xy_D2y = I_TWOBODYOVERLAP_I4x2y_Py+ABY*I_TWOBODYOVERLAP_H4xy_Py;
  Double I_TWOBODYOVERLAP_H4xz_D2y = I_TWOBODYOVERLAP_I4xyz_Py+ABY*I_TWOBODYOVERLAP_H4xz_Py;
  Double I_TWOBODYOVERLAP_H3x2y_D2y = I_TWOBODYOVERLAP_I3x3y_Py+ABY*I_TWOBODYOVERLAP_H3x2y_Py;
  Double I_TWOBODYOVERLAP_H3xyz_D2y = I_TWOBODYOVERLAP_I3x2yz_Py+ABY*I_TWOBODYOVERLAP_H3xyz_Py;
  Double I_TWOBODYOVERLAP_H3x2z_D2y = I_TWOBODYOVERLAP_I3xy2z_Py+ABY*I_TWOBODYOVERLAP_H3x2z_Py;
  Double I_TWOBODYOVERLAP_H2x3y_D2y = I_TWOBODYOVERLAP_I2x4y_Py+ABY*I_TWOBODYOVERLAP_H2x3y_Py;
  Double I_TWOBODYOVERLAP_H2x2yz_D2y = I_TWOBODYOVERLAP_I2x3yz_Py+ABY*I_TWOBODYOVERLAP_H2x2yz_Py;
  Double I_TWOBODYOVERLAP_H2xy2z_D2y = I_TWOBODYOVERLAP_I2x2y2z_Py+ABY*I_TWOBODYOVERLAP_H2xy2z_Py;
  Double I_TWOBODYOVERLAP_H2x3z_D2y = I_TWOBODYOVERLAP_I2xy3z_Py+ABY*I_TWOBODYOVERLAP_H2x3z_Py;
  Double I_TWOBODYOVERLAP_Hx4y_D2y = I_TWOBODYOVERLAP_Ix5y_Py+ABY*I_TWOBODYOVERLAP_Hx4y_Py;
  Double I_TWOBODYOVERLAP_Hx3yz_D2y = I_TWOBODYOVERLAP_Ix4yz_Py+ABY*I_TWOBODYOVERLAP_Hx3yz_Py;
  Double I_TWOBODYOVERLAP_Hx2y2z_D2y = I_TWOBODYOVERLAP_Ix3y2z_Py+ABY*I_TWOBODYOVERLAP_Hx2y2z_Py;
  Double I_TWOBODYOVERLAP_Hxy3z_D2y = I_TWOBODYOVERLAP_Ix2y3z_Py+ABY*I_TWOBODYOVERLAP_Hxy3z_Py;
  Double I_TWOBODYOVERLAP_Hx4z_D2y = I_TWOBODYOVERLAP_Ixy4z_Py+ABY*I_TWOBODYOVERLAP_Hx4z_Py;
  Double I_TWOBODYOVERLAP_H5y_D2y = I_TWOBODYOVERLAP_I6y_Py+ABY*I_TWOBODYOVERLAP_H5y_Py;
  Double I_TWOBODYOVERLAP_H4yz_D2y = I_TWOBODYOVERLAP_I5yz_Py+ABY*I_TWOBODYOVERLAP_H4yz_Py;
  Double I_TWOBODYOVERLAP_H3y2z_D2y = I_TWOBODYOVERLAP_I4y2z_Py+ABY*I_TWOBODYOVERLAP_H3y2z_Py;
  Double I_TWOBODYOVERLAP_H2y3z_D2y = I_TWOBODYOVERLAP_I3y3z_Py+ABY*I_TWOBODYOVERLAP_H2y3z_Py;
  Double I_TWOBODYOVERLAP_Hy4z_D2y = I_TWOBODYOVERLAP_I2y4z_Py+ABY*I_TWOBODYOVERLAP_Hy4z_Py;
  Double I_TWOBODYOVERLAP_H5z_D2y = I_TWOBODYOVERLAP_Iy5z_Py+ABY*I_TWOBODYOVERLAP_H5z_Py;
  Double I_TWOBODYOVERLAP_H5x_D2z = I_TWOBODYOVERLAP_I5xz_Pz+ABZ*I_TWOBODYOVERLAP_H5x_Pz;
  Double I_TWOBODYOVERLAP_H4xy_D2z = I_TWOBODYOVERLAP_I4xyz_Pz+ABZ*I_TWOBODYOVERLAP_H4xy_Pz;
  Double I_TWOBODYOVERLAP_H4xz_D2z = I_TWOBODYOVERLAP_I4x2z_Pz+ABZ*I_TWOBODYOVERLAP_H4xz_Pz;
  Double I_TWOBODYOVERLAP_H3x2y_D2z = I_TWOBODYOVERLAP_I3x2yz_Pz+ABZ*I_TWOBODYOVERLAP_H3x2y_Pz;
  Double I_TWOBODYOVERLAP_H3xyz_D2z = I_TWOBODYOVERLAP_I3xy2z_Pz+ABZ*I_TWOBODYOVERLAP_H3xyz_Pz;
  Double I_TWOBODYOVERLAP_H3x2z_D2z = I_TWOBODYOVERLAP_I3x3z_Pz+ABZ*I_TWOBODYOVERLAP_H3x2z_Pz;
  Double I_TWOBODYOVERLAP_H2x3y_D2z = I_TWOBODYOVERLAP_I2x3yz_Pz+ABZ*I_TWOBODYOVERLAP_H2x3y_Pz;
  Double I_TWOBODYOVERLAP_H2x2yz_D2z = I_TWOBODYOVERLAP_I2x2y2z_Pz+ABZ*I_TWOBODYOVERLAP_H2x2yz_Pz;
  Double I_TWOBODYOVERLAP_H2xy2z_D2z = I_TWOBODYOVERLAP_I2xy3z_Pz+ABZ*I_TWOBODYOVERLAP_H2xy2z_Pz;
  Double I_TWOBODYOVERLAP_H2x3z_D2z = I_TWOBODYOVERLAP_I2x4z_Pz+ABZ*I_TWOBODYOVERLAP_H2x3z_Pz;
  Double I_TWOBODYOVERLAP_Hx4y_D2z = I_TWOBODYOVERLAP_Ix4yz_Pz+ABZ*I_TWOBODYOVERLAP_Hx4y_Pz;
  Double I_TWOBODYOVERLAP_Hx3yz_D2z = I_TWOBODYOVERLAP_Ix3y2z_Pz+ABZ*I_TWOBODYOVERLAP_Hx3yz_Pz;
  Double I_TWOBODYOVERLAP_Hx2y2z_D2z = I_TWOBODYOVERLAP_Ix2y3z_Pz+ABZ*I_TWOBODYOVERLAP_Hx2y2z_Pz;
  Double I_TWOBODYOVERLAP_Hxy3z_D2z = I_TWOBODYOVERLAP_Ixy4z_Pz+ABZ*I_TWOBODYOVERLAP_Hxy3z_Pz;
  Double I_TWOBODYOVERLAP_Hx4z_D2z = I_TWOBODYOVERLAP_Ix5z_Pz+ABZ*I_TWOBODYOVERLAP_Hx4z_Pz;
  Double I_TWOBODYOVERLAP_H5y_D2z = I_TWOBODYOVERLAP_I5yz_Pz+ABZ*I_TWOBODYOVERLAP_H5y_Pz;
  Double I_TWOBODYOVERLAP_H4yz_D2z = I_TWOBODYOVERLAP_I4y2z_Pz+ABZ*I_TWOBODYOVERLAP_H4yz_Pz;
  Double I_TWOBODYOVERLAP_H3y2z_D2z = I_TWOBODYOVERLAP_I3y3z_Pz+ABZ*I_TWOBODYOVERLAP_H3y2z_Pz;
  Double I_TWOBODYOVERLAP_H2y3z_D2z = I_TWOBODYOVERLAP_I2y4z_Pz+ABZ*I_TWOBODYOVERLAP_H2y3z_Pz;
  Double I_TWOBODYOVERLAP_Hy4z_D2z = I_TWOBODYOVERLAP_Iy5z_Pz+ABZ*I_TWOBODYOVERLAP_Hy4z_Pz;
  Double I_TWOBODYOVERLAP_H5z_D2z = I_TWOBODYOVERLAP_I6z_Pz+ABZ*I_TWOBODYOVERLAP_H5z_Pz;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
   ************************************************************/
  abcd[0] = I_TWOBODYOVERLAP_H5x_D2x+ABX*I_TWOBODYOVERLAP_G4x_D2x;
  abcd[1] = I_TWOBODYOVERLAP_H4xy_D2x+ABX*I_TWOBODYOVERLAP_G3xy_D2x;
  abcd[2] = I_TWOBODYOVERLAP_H4xz_D2x+ABX*I_TWOBODYOVERLAP_G3xz_D2x;
  abcd[3] = I_TWOBODYOVERLAP_H3x2y_D2x+ABX*I_TWOBODYOVERLAP_G2x2y_D2x;
  abcd[4] = I_TWOBODYOVERLAP_H3xyz_D2x+ABX*I_TWOBODYOVERLAP_G2xyz_D2x;
  abcd[5] = I_TWOBODYOVERLAP_H3x2z_D2x+ABX*I_TWOBODYOVERLAP_G2x2z_D2x;
  abcd[6] = I_TWOBODYOVERLAP_H2x3y_D2x+ABX*I_TWOBODYOVERLAP_Gx3y_D2x;
  abcd[7] = I_TWOBODYOVERLAP_H2x2yz_D2x+ABX*I_TWOBODYOVERLAP_Gx2yz_D2x;
  abcd[8] = I_TWOBODYOVERLAP_H2xy2z_D2x+ABX*I_TWOBODYOVERLAP_Gxy2z_D2x;
  abcd[9] = I_TWOBODYOVERLAP_H2x3z_D2x+ABX*I_TWOBODYOVERLAP_Gx3z_D2x;
  abcd[10] = I_TWOBODYOVERLAP_Hx4y_D2x+ABX*I_TWOBODYOVERLAP_G4y_D2x;
  abcd[11] = I_TWOBODYOVERLAP_Hx3yz_D2x+ABX*I_TWOBODYOVERLAP_G3yz_D2x;
  abcd[12] = I_TWOBODYOVERLAP_Hx2y2z_D2x+ABX*I_TWOBODYOVERLAP_G2y2z_D2x;
  abcd[13] = I_TWOBODYOVERLAP_Hxy3z_D2x+ABX*I_TWOBODYOVERLAP_Gy3z_D2x;
  abcd[14] = I_TWOBODYOVERLAP_Hx4z_D2x+ABX*I_TWOBODYOVERLAP_G4z_D2x;
  abcd[15] = I_TWOBODYOVERLAP_H4xy_D2x+ABY*I_TWOBODYOVERLAP_G4x_D2x;
  abcd[16] = I_TWOBODYOVERLAP_H3x2y_D2x+ABY*I_TWOBODYOVERLAP_G3xy_D2x;
  abcd[17] = I_TWOBODYOVERLAP_H3xyz_D2x+ABY*I_TWOBODYOVERLAP_G3xz_D2x;
  abcd[18] = I_TWOBODYOVERLAP_H2x3y_D2x+ABY*I_TWOBODYOVERLAP_G2x2y_D2x;
  abcd[19] = I_TWOBODYOVERLAP_H2x2yz_D2x+ABY*I_TWOBODYOVERLAP_G2xyz_D2x;
  abcd[20] = I_TWOBODYOVERLAP_H2xy2z_D2x+ABY*I_TWOBODYOVERLAP_G2x2z_D2x;
  abcd[21] = I_TWOBODYOVERLAP_Hx4y_D2x+ABY*I_TWOBODYOVERLAP_Gx3y_D2x;
  abcd[22] = I_TWOBODYOVERLAP_Hx3yz_D2x+ABY*I_TWOBODYOVERLAP_Gx2yz_D2x;
  abcd[23] = I_TWOBODYOVERLAP_Hx2y2z_D2x+ABY*I_TWOBODYOVERLAP_Gxy2z_D2x;
  abcd[24] = I_TWOBODYOVERLAP_Hxy3z_D2x+ABY*I_TWOBODYOVERLAP_Gx3z_D2x;
  abcd[25] = I_TWOBODYOVERLAP_H5y_D2x+ABY*I_TWOBODYOVERLAP_G4y_D2x;
  abcd[26] = I_TWOBODYOVERLAP_H4yz_D2x+ABY*I_TWOBODYOVERLAP_G3yz_D2x;
  abcd[27] = I_TWOBODYOVERLAP_H3y2z_D2x+ABY*I_TWOBODYOVERLAP_G2y2z_D2x;
  abcd[28] = I_TWOBODYOVERLAP_H2y3z_D2x+ABY*I_TWOBODYOVERLAP_Gy3z_D2x;
  abcd[29] = I_TWOBODYOVERLAP_Hy4z_D2x+ABY*I_TWOBODYOVERLAP_G4z_D2x;
  abcd[30] = I_TWOBODYOVERLAP_H4xz_D2x+ABZ*I_TWOBODYOVERLAP_G4x_D2x;
  abcd[31] = I_TWOBODYOVERLAP_H3xyz_D2x+ABZ*I_TWOBODYOVERLAP_G3xy_D2x;
  abcd[32] = I_TWOBODYOVERLAP_H3x2z_D2x+ABZ*I_TWOBODYOVERLAP_G3xz_D2x;
  abcd[33] = I_TWOBODYOVERLAP_H2x2yz_D2x+ABZ*I_TWOBODYOVERLAP_G2x2y_D2x;
  abcd[34] = I_TWOBODYOVERLAP_H2xy2z_D2x+ABZ*I_TWOBODYOVERLAP_G2xyz_D2x;
  abcd[35] = I_TWOBODYOVERLAP_H2x3z_D2x+ABZ*I_TWOBODYOVERLAP_G2x2z_D2x;
  abcd[36] = I_TWOBODYOVERLAP_Hx3yz_D2x+ABZ*I_TWOBODYOVERLAP_Gx3y_D2x;
  abcd[37] = I_TWOBODYOVERLAP_Hx2y2z_D2x+ABZ*I_TWOBODYOVERLAP_Gx2yz_D2x;
  abcd[38] = I_TWOBODYOVERLAP_Hxy3z_D2x+ABZ*I_TWOBODYOVERLAP_Gxy2z_D2x;
  abcd[39] = I_TWOBODYOVERLAP_Hx4z_D2x+ABZ*I_TWOBODYOVERLAP_Gx3z_D2x;
  abcd[40] = I_TWOBODYOVERLAP_H4yz_D2x+ABZ*I_TWOBODYOVERLAP_G4y_D2x;
  abcd[41] = I_TWOBODYOVERLAP_H3y2z_D2x+ABZ*I_TWOBODYOVERLAP_G3yz_D2x;
  abcd[42] = I_TWOBODYOVERLAP_H2y3z_D2x+ABZ*I_TWOBODYOVERLAP_G2y2z_D2x;
  abcd[43] = I_TWOBODYOVERLAP_Hy4z_D2x+ABZ*I_TWOBODYOVERLAP_Gy3z_D2x;
  abcd[44] = I_TWOBODYOVERLAP_H5z_D2x+ABZ*I_TWOBODYOVERLAP_G4z_D2x;
  abcd[45] = I_TWOBODYOVERLAP_H5x_D2y+ABX*I_TWOBODYOVERLAP_G4x_D2y;
  abcd[46] = I_TWOBODYOVERLAP_H4xy_D2y+ABX*I_TWOBODYOVERLAP_G3xy_D2y;
  abcd[47] = I_TWOBODYOVERLAP_H4xz_D2y+ABX*I_TWOBODYOVERLAP_G3xz_D2y;
  abcd[48] = I_TWOBODYOVERLAP_H3x2y_D2y+ABX*I_TWOBODYOVERLAP_G2x2y_D2y;
  abcd[49] = I_TWOBODYOVERLAP_H3xyz_D2y+ABX*I_TWOBODYOVERLAP_G2xyz_D2y;
  abcd[50] = I_TWOBODYOVERLAP_H3x2z_D2y+ABX*I_TWOBODYOVERLAP_G2x2z_D2y;
  abcd[51] = I_TWOBODYOVERLAP_H2x3y_D2y+ABX*I_TWOBODYOVERLAP_Gx3y_D2y;
  abcd[52] = I_TWOBODYOVERLAP_H2x2yz_D2y+ABX*I_TWOBODYOVERLAP_Gx2yz_D2y;
  abcd[53] = I_TWOBODYOVERLAP_H2xy2z_D2y+ABX*I_TWOBODYOVERLAP_Gxy2z_D2y;
  abcd[54] = I_TWOBODYOVERLAP_H2x3z_D2y+ABX*I_TWOBODYOVERLAP_Gx3z_D2y;
  abcd[55] = I_TWOBODYOVERLAP_Hx4y_D2y+ABX*I_TWOBODYOVERLAP_G4y_D2y;
  abcd[56] = I_TWOBODYOVERLAP_Hx3yz_D2y+ABX*I_TWOBODYOVERLAP_G3yz_D2y;
  abcd[57] = I_TWOBODYOVERLAP_Hx2y2z_D2y+ABX*I_TWOBODYOVERLAP_G2y2z_D2y;
  abcd[58] = I_TWOBODYOVERLAP_Hxy3z_D2y+ABX*I_TWOBODYOVERLAP_Gy3z_D2y;
  abcd[59] = I_TWOBODYOVERLAP_Hx4z_D2y+ABX*I_TWOBODYOVERLAP_G4z_D2y;
  abcd[60] = I_TWOBODYOVERLAP_H4xz_Dxy+ABZ*I_TWOBODYOVERLAP_G4x_Dxy;
  abcd[61] = I_TWOBODYOVERLAP_H3xyz_Dxy+ABZ*I_TWOBODYOVERLAP_G3xy_Dxy;
  abcd[62] = I_TWOBODYOVERLAP_H3x2z_Dxy+ABZ*I_TWOBODYOVERLAP_G3xz_Dxy;
  abcd[63] = I_TWOBODYOVERLAP_H2x2yz_Dxy+ABZ*I_TWOBODYOVERLAP_G2x2y_Dxy;
  abcd[64] = I_TWOBODYOVERLAP_H2xy2z_Dxy+ABZ*I_TWOBODYOVERLAP_G2xyz_Dxy;
  abcd[65] = I_TWOBODYOVERLAP_H2x3z_Dxy+ABZ*I_TWOBODYOVERLAP_G2x2z_Dxy;
  abcd[66] = I_TWOBODYOVERLAP_Hx3yz_Dxy+ABZ*I_TWOBODYOVERLAP_Gx3y_Dxy;
  abcd[67] = I_TWOBODYOVERLAP_Hx2y2z_Dxy+ABZ*I_TWOBODYOVERLAP_Gx2yz_Dxy;
  abcd[68] = I_TWOBODYOVERLAP_Hxy3z_Dxy+ABZ*I_TWOBODYOVERLAP_Gxy2z_Dxy;
  abcd[69] = I_TWOBODYOVERLAP_Hx4z_Dxy+ABZ*I_TWOBODYOVERLAP_Gx3z_Dxy;
  abcd[70] = I_TWOBODYOVERLAP_H4yz_Dxy+ABZ*I_TWOBODYOVERLAP_G4y_Dxy;
  abcd[71] = I_TWOBODYOVERLAP_H3y2z_Dxy+ABZ*I_TWOBODYOVERLAP_G3yz_Dxy;
  abcd[72] = I_TWOBODYOVERLAP_H2y3z_Dxy+ABZ*I_TWOBODYOVERLAP_G2y2z_Dxy;
  abcd[73] = I_TWOBODYOVERLAP_Hy4z_Dxy+ABZ*I_TWOBODYOVERLAP_Gy3z_Dxy;
  abcd[74] = I_TWOBODYOVERLAP_H5z_Dxy+ABZ*I_TWOBODYOVERLAP_G4z_Dxy;
  abcd[75] = I_TWOBODYOVERLAP_H5x_D2z+ABX*I_TWOBODYOVERLAP_G4x_D2z;
  abcd[76] = I_TWOBODYOVERLAP_H4xy_D2z+ABX*I_TWOBODYOVERLAP_G3xy_D2z;
  abcd[77] = I_TWOBODYOVERLAP_H4xz_D2z+ABX*I_TWOBODYOVERLAP_G3xz_D2z;
  abcd[78] = I_TWOBODYOVERLAP_H3x2y_D2z+ABX*I_TWOBODYOVERLAP_G2x2y_D2z;
  abcd[79] = I_TWOBODYOVERLAP_H3xyz_D2z+ABX*I_TWOBODYOVERLAP_G2xyz_D2z;
  abcd[80] = I_TWOBODYOVERLAP_H3x2z_D2z+ABX*I_TWOBODYOVERLAP_G2x2z_D2z;
  abcd[81] = I_TWOBODYOVERLAP_H2x3y_D2z+ABX*I_TWOBODYOVERLAP_Gx3y_D2z;
  abcd[82] = I_TWOBODYOVERLAP_H2x2yz_D2z+ABX*I_TWOBODYOVERLAP_Gx2yz_D2z;
  abcd[83] = I_TWOBODYOVERLAP_H2xy2z_D2z+ABX*I_TWOBODYOVERLAP_Gxy2z_D2z;
  abcd[84] = I_TWOBODYOVERLAP_H2x3z_D2z+ABX*I_TWOBODYOVERLAP_Gx3z_D2z;
  abcd[85] = I_TWOBODYOVERLAP_Hx4y_D2z+ABX*I_TWOBODYOVERLAP_G4y_D2z;
  abcd[86] = I_TWOBODYOVERLAP_Hx3yz_D2z+ABX*I_TWOBODYOVERLAP_G3yz_D2z;
  abcd[87] = I_TWOBODYOVERLAP_Hx2y2z_D2z+ABX*I_TWOBODYOVERLAP_G2y2z_D2z;
  abcd[88] = I_TWOBODYOVERLAP_Hxy3z_D2z+ABX*I_TWOBODYOVERLAP_Gy3z_D2z;
  abcd[89] = I_TWOBODYOVERLAP_Hx4z_D2z+ABX*I_TWOBODYOVERLAP_G4z_D2z;
  abcd[90] = I_TWOBODYOVERLAP_H4xy_D2y+ABY*I_TWOBODYOVERLAP_G4x_D2y;
  abcd[91] = I_TWOBODYOVERLAP_H3x2y_D2y+ABY*I_TWOBODYOVERLAP_G3xy_D2y;
  abcd[92] = I_TWOBODYOVERLAP_H3xyz_D2y+ABY*I_TWOBODYOVERLAP_G3xz_D2y;
  abcd[93] = I_TWOBODYOVERLAP_H2x3y_D2y+ABY*I_TWOBODYOVERLAP_G2x2y_D2y;
  abcd[94] = I_TWOBODYOVERLAP_H2x2yz_D2y+ABY*I_TWOBODYOVERLAP_G2xyz_D2y;
  abcd[95] = I_TWOBODYOVERLAP_H2xy2z_D2y+ABY*I_TWOBODYOVERLAP_G2x2z_D2y;
  abcd[96] = I_TWOBODYOVERLAP_Hx4y_D2y+ABY*I_TWOBODYOVERLAP_Gx3y_D2y;
  abcd[97] = I_TWOBODYOVERLAP_Hx3yz_D2y+ABY*I_TWOBODYOVERLAP_Gx2yz_D2y;
  abcd[98] = I_TWOBODYOVERLAP_Hx2y2z_D2y+ABY*I_TWOBODYOVERLAP_Gxy2z_D2y;
  abcd[99] = I_TWOBODYOVERLAP_Hxy3z_D2y+ABY*I_TWOBODYOVERLAP_Gx3z_D2y;
  abcd[100] = I_TWOBODYOVERLAP_H5y_D2y+ABY*I_TWOBODYOVERLAP_G4y_D2y;
  abcd[101] = I_TWOBODYOVERLAP_H4yz_D2y+ABY*I_TWOBODYOVERLAP_G3yz_D2y;
  abcd[102] = I_TWOBODYOVERLAP_H3y2z_D2y+ABY*I_TWOBODYOVERLAP_G2y2z_D2y;
  abcd[103] = I_TWOBODYOVERLAP_H2y3z_D2y+ABY*I_TWOBODYOVERLAP_Gy3z_D2y;
  abcd[104] = I_TWOBODYOVERLAP_Hy4z_D2y+ABY*I_TWOBODYOVERLAP_G4z_D2y;
  abcd[105] = I_TWOBODYOVERLAP_H4xz_D2y+ABZ*I_TWOBODYOVERLAP_G4x_D2y;
  abcd[106] = I_TWOBODYOVERLAP_H3xyz_D2y+ABZ*I_TWOBODYOVERLAP_G3xy_D2y;
  abcd[107] = I_TWOBODYOVERLAP_H3x2z_D2y+ABZ*I_TWOBODYOVERLAP_G3xz_D2y;
  abcd[108] = I_TWOBODYOVERLAP_H2x2yz_D2y+ABZ*I_TWOBODYOVERLAP_G2x2y_D2y;
  abcd[109] = I_TWOBODYOVERLAP_H2xy2z_D2y+ABZ*I_TWOBODYOVERLAP_G2xyz_D2y;
  abcd[110] = I_TWOBODYOVERLAP_H2x3z_D2y+ABZ*I_TWOBODYOVERLAP_G2x2z_D2y;
  abcd[111] = I_TWOBODYOVERLAP_Hx3yz_D2y+ABZ*I_TWOBODYOVERLAP_Gx3y_D2y;
  abcd[112] = I_TWOBODYOVERLAP_Hx2y2z_D2y+ABZ*I_TWOBODYOVERLAP_Gx2yz_D2y;
  abcd[113] = I_TWOBODYOVERLAP_Hxy3z_D2y+ABZ*I_TWOBODYOVERLAP_Gxy2z_D2y;
  abcd[114] = I_TWOBODYOVERLAP_Hx4z_D2y+ABZ*I_TWOBODYOVERLAP_Gx3z_D2y;
  abcd[115] = I_TWOBODYOVERLAP_H4yz_D2y+ABZ*I_TWOBODYOVERLAP_G4y_D2y;
  abcd[116] = I_TWOBODYOVERLAP_H3y2z_D2y+ABZ*I_TWOBODYOVERLAP_G3yz_D2y;
  abcd[117] = I_TWOBODYOVERLAP_H2y3z_D2y+ABZ*I_TWOBODYOVERLAP_G2y2z_D2y;
  abcd[118] = I_TWOBODYOVERLAP_Hy4z_D2y+ABZ*I_TWOBODYOVERLAP_Gy3z_D2y;
  abcd[119] = I_TWOBODYOVERLAP_H5z_D2y+ABZ*I_TWOBODYOVERLAP_G4z_D2y;
  abcd[120] = I_TWOBODYOVERLAP_H4xy_D2z+ABY*I_TWOBODYOVERLAP_G4x_D2z;
  abcd[121] = I_TWOBODYOVERLAP_H3x2y_D2z+ABY*I_TWOBODYOVERLAP_G3xy_D2z;
  abcd[122] = I_TWOBODYOVERLAP_H3xyz_D2z+ABY*I_TWOBODYOVERLAP_G3xz_D2z;
  abcd[123] = I_TWOBODYOVERLAP_H2x3y_D2z+ABY*I_TWOBODYOVERLAP_G2x2y_D2z;
  abcd[124] = I_TWOBODYOVERLAP_H2x2yz_D2z+ABY*I_TWOBODYOVERLAP_G2xyz_D2z;
  abcd[125] = I_TWOBODYOVERLAP_H2xy2z_D2z+ABY*I_TWOBODYOVERLAP_G2x2z_D2z;
  abcd[126] = I_TWOBODYOVERLAP_Hx4y_D2z+ABY*I_TWOBODYOVERLAP_Gx3y_D2z;
  abcd[127] = I_TWOBODYOVERLAP_Hx3yz_D2z+ABY*I_TWOBODYOVERLAP_Gx2yz_D2z;
  abcd[128] = I_TWOBODYOVERLAP_Hx2y2z_D2z+ABY*I_TWOBODYOVERLAP_Gxy2z_D2z;
  abcd[129] = I_TWOBODYOVERLAP_Hxy3z_D2z+ABY*I_TWOBODYOVERLAP_Gx3z_D2z;
  abcd[130] = I_TWOBODYOVERLAP_H5y_D2z+ABY*I_TWOBODYOVERLAP_G4y_D2z;
  abcd[131] = I_TWOBODYOVERLAP_H4yz_D2z+ABY*I_TWOBODYOVERLAP_G3yz_D2z;
  abcd[132] = I_TWOBODYOVERLAP_H3y2z_D2z+ABY*I_TWOBODYOVERLAP_G2y2z_D2z;
  abcd[133] = I_TWOBODYOVERLAP_H2y3z_D2z+ABY*I_TWOBODYOVERLAP_Gy3z_D2z;
  abcd[134] = I_TWOBODYOVERLAP_Hy4z_D2z+ABY*I_TWOBODYOVERLAP_G4z_D2z;
  abcd[135] = I_TWOBODYOVERLAP_H4xz_D2z+ABZ*I_TWOBODYOVERLAP_G4x_D2z;
  abcd[136] = I_TWOBODYOVERLAP_H3xyz_D2z+ABZ*I_TWOBODYOVERLAP_G3xy_D2z;
  abcd[137] = I_TWOBODYOVERLAP_H3x2z_D2z+ABZ*I_TWOBODYOVERLAP_G3xz_D2z;
  abcd[138] = I_TWOBODYOVERLAP_H2x2yz_D2z+ABZ*I_TWOBODYOVERLAP_G2x2y_D2z;
  abcd[139] = I_TWOBODYOVERLAP_H2xy2z_D2z+ABZ*I_TWOBODYOVERLAP_G2xyz_D2z;
  abcd[140] = I_TWOBODYOVERLAP_H2x3z_D2z+ABZ*I_TWOBODYOVERLAP_G2x2z_D2z;
  abcd[141] = I_TWOBODYOVERLAP_Hx3yz_D2z+ABZ*I_TWOBODYOVERLAP_Gx3y_D2z;
  abcd[142] = I_TWOBODYOVERLAP_Hx2y2z_D2z+ABZ*I_TWOBODYOVERLAP_Gx2yz_D2z;
  abcd[143] = I_TWOBODYOVERLAP_Hxy3z_D2z+ABZ*I_TWOBODYOVERLAP_Gxy2z_D2z;
  abcd[144] = I_TWOBODYOVERLAP_Hx4z_D2z+ABZ*I_TWOBODYOVERLAP_Gx3z_D2z;
  abcd[145] = I_TWOBODYOVERLAP_H4yz_D2z+ABZ*I_TWOBODYOVERLAP_G4y_D2z;
  abcd[146] = I_TWOBODYOVERLAP_H3y2z_D2z+ABZ*I_TWOBODYOVERLAP_G3yz_D2z;
  abcd[147] = I_TWOBODYOVERLAP_H2y3z_D2z+ABZ*I_TWOBODYOVERLAP_G2y2z_D2z;
  abcd[148] = I_TWOBODYOVERLAP_Hy4z_D2z+ABZ*I_TWOBODYOVERLAP_Gy3z_D2z;
  abcd[149] = I_TWOBODYOVERLAP_H5z_D2z+ABZ*I_TWOBODYOVERLAP_G4z_D2z;
}
