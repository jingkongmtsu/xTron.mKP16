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

void hgp_os_twobodyoverlap_h_g(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_M9x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M8xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M8xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M7x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M7xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M7x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M6x3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M6x2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M6xy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M6x3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M5xy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x5y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x4yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x3y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x2y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M4xy4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x6y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x5yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x4y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x3y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x2y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M3xy5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x7y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x6yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x5y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x4y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x3y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x2y5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M2xy6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x7z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx8y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx7yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx6y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx5y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx4y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx3y5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx2y6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Mxy7z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx8z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M9y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M8yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M7y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M6y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M5y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M4y5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M3y6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M2y7z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_My8z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_M9z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L8x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L6xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L5xy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L4xy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x4yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x3y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x2y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L3xy4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x5yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x4y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x3y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x2y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L2xy5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx6yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx5y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx4y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx3y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx2y5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Lxy6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L8y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L7yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L6y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L5y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L4y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L3y5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L2y6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ly7z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_L8z_S = 0.0E0;
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
     * totally 4 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_S_vrr = PAX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_S_vrr = PAY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_S_vrr = PAY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2y_S_vrr;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_L_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_L8x_S_vrr = PAX*I_TWOBODYOVERLAP_K7x_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L7xy_S_vrr = PAY*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_L7xz_S_vrr = PAZ*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_L6x2y_S_vrr = PAY*I_TWOBODYOVERLAP_K6xy_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L6xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_K6xy_S_vrr;
    Double I_TWOBODYOVERLAP_L6x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K6xz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L5x3y_S_vrr = PAY*I_TWOBODYOVERLAP_K5x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_L5x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L5xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    Double I_TWOBODYOVERLAP_L5x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_K5x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_L4x4y_S_vrr = PAY*I_TWOBODYOVERLAP_K4x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L4x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    Double I_TWOBODYOVERLAP_L4x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L4xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    Double I_TWOBODYOVERLAP_L4x4z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_L3x5y_S_vrr = PAX*I_TWOBODYOVERLAP_K2x5y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K3x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_K3xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_L3xy4z_S_vrr = PAY*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    Double I_TWOBODYOVERLAP_L3x5z_S_vrr = PAX*I_TWOBODYOVERLAP_K2x5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x6y_S_vrr = PAX*I_TWOBODYOVERLAP_Kx6y_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K2x4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x3y3z_S_vrr = PAX*I_TWOBODYOVERLAP_Kx3y3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_K2xy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_L2xy5z_S_vrr = PAY*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x6z_S_vrr = PAX*I_TWOBODYOVERLAP_Kx6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx7y_S_vrr = PAX*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_Lx6yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    Double I_TWOBODYOVERLAP_Lx5y2z_S_vrr = PAX*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx4y3z_S_vrr = PAX*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx3y4z_S_vrr = PAX*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx2y5z_S_vrr = PAX*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    Double I_TWOBODYOVERLAP_Lxy6z_S_vrr = PAY*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx7z_S_vrr = PAX*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_L8y_S_vrr = PAY*I_TWOBODYOVERLAP_K7y_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L7yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_L6y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K6yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L5y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_K5y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_L4y4z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_L3y5z_S_vrr = PAY*I_TWOBODYOVERLAP_K2y5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2y6z_S_vrr = PAY*I_TWOBODYOVERLAP_Ky6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_Ly7z_S_vrr = PAY*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_L8z_S_vrr = PAZ*I_TWOBODYOVERLAP_K7z_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_M_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_M9x_S_vrr = PAX*I_TWOBODYOVERLAP_L8x_S_vrr+8*oned2z*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_M8xy_S_vrr = PAY*I_TWOBODYOVERLAP_L8x_S_vrr;
    Double I_TWOBODYOVERLAP_M8xz_S_vrr = PAZ*I_TWOBODYOVERLAP_L8x_S_vrr;
    Double I_TWOBODYOVERLAP_M7x2y_S_vrr = PAY*I_TWOBODYOVERLAP_L7xy_S_vrr+oned2z*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_M7xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_L7xy_S_vrr;
    Double I_TWOBODYOVERLAP_M7x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L7xz_S_vrr+oned2z*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_M6x3y_S_vrr = PAY*I_TWOBODYOVERLAP_L6x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_K6xy_S_vrr;
    Double I_TWOBODYOVERLAP_M6x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L6x2y_S_vrr;
    Double I_TWOBODYOVERLAP_M6xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_L6x2z_S_vrr;
    Double I_TWOBODYOVERLAP_M6x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_L6x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_K6xz_S_vrr;
    Double I_TWOBODYOVERLAP_M5x4y_S_vrr = PAY*I_TWOBODYOVERLAP_L5x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    Double I_TWOBODYOVERLAP_M5x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L5x3y_S_vrr;
    Double I_TWOBODYOVERLAP_M5x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L5x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    Double I_TWOBODYOVERLAP_M5xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_L5x3z_S_vrr;
    Double I_TWOBODYOVERLAP_M5x4z_S_vrr = PAZ*I_TWOBODYOVERLAP_L5x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    Double I_TWOBODYOVERLAP_M4x5y_S_vrr = PAX*I_TWOBODYOVERLAP_L3x5y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    Double I_TWOBODYOVERLAP_M4x4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L4x4y_S_vrr;
    Double I_TWOBODYOVERLAP_M4x3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L4x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    Double I_TWOBODYOVERLAP_M4x2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_L4xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    Double I_TWOBODYOVERLAP_M4xy4z_S_vrr = PAY*I_TWOBODYOVERLAP_L4x4z_S_vrr;
    Double I_TWOBODYOVERLAP_M4x5z_S_vrr = PAX*I_TWOBODYOVERLAP_L3x5z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    Double I_TWOBODYOVERLAP_M3x6y_S_vrr = PAX*I_TWOBODYOVERLAP_L2x6y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    Double I_TWOBODYOVERLAP_M3x5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L3x5y_S_vrr;
    Double I_TWOBODYOVERLAP_M3x4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L3x4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    Double I_TWOBODYOVERLAP_M3x3y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_L3x3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    Double I_TWOBODYOVERLAP_M3x2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_L3xy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    Double I_TWOBODYOVERLAP_M3xy5z_S_vrr = PAY*I_TWOBODYOVERLAP_L3x5z_S_vrr;
    Double I_TWOBODYOVERLAP_M3x6z_S_vrr = PAX*I_TWOBODYOVERLAP_L2x6z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    Double I_TWOBODYOVERLAP_M2x7y_S_vrr = PAX*I_TWOBODYOVERLAP_Lx7y_S_vrr+oned2z*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_M2x6yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L2x6y_S_vrr;
    Double I_TWOBODYOVERLAP_M2x5y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L2x5yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    Double I_TWOBODYOVERLAP_M2x4y3z_S_vrr = PAX*I_TWOBODYOVERLAP_Lx4y3z_S_vrr+oned2z*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    Double I_TWOBODYOVERLAP_M2x3y4z_S_vrr = PAX*I_TWOBODYOVERLAP_Lx3y4z_S_vrr+oned2z*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    Double I_TWOBODYOVERLAP_M2x2y5z_S_vrr = PAY*I_TWOBODYOVERLAP_L2xy5z_S_vrr+oned2z*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    Double I_TWOBODYOVERLAP_M2xy6z_S_vrr = PAY*I_TWOBODYOVERLAP_L2x6z_S_vrr;
    Double I_TWOBODYOVERLAP_M2x7z_S_vrr = PAX*I_TWOBODYOVERLAP_Lx7z_S_vrr+oned2z*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx8y_S_vrr = PAX*I_TWOBODYOVERLAP_L8y_S_vrr;
    Double I_TWOBODYOVERLAP_Mx7yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Lx7y_S_vrr;
    Double I_TWOBODYOVERLAP_Mx6y2z_S_vrr = PAX*I_TWOBODYOVERLAP_L6y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx5y3z_S_vrr = PAX*I_TWOBODYOVERLAP_L5y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx4y4z_S_vrr = PAX*I_TWOBODYOVERLAP_L4y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx3y5z_S_vrr = PAX*I_TWOBODYOVERLAP_L3y5z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx2y6z_S_vrr = PAX*I_TWOBODYOVERLAP_L2y6z_S_vrr;
    Double I_TWOBODYOVERLAP_Mxy7z_S_vrr = PAY*I_TWOBODYOVERLAP_Lx7z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx8z_S_vrr = PAX*I_TWOBODYOVERLAP_L8z_S_vrr;
    Double I_TWOBODYOVERLAP_M9y_S_vrr = PAY*I_TWOBODYOVERLAP_L8y_S_vrr+8*oned2z*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_M8yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L8y_S_vrr;
    Double I_TWOBODYOVERLAP_M7y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L7yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_M6y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_L6y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_K6yz_S_vrr;
    Double I_TWOBODYOVERLAP_M5y4z_S_vrr = PAZ*I_TWOBODYOVERLAP_L5y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    Double I_TWOBODYOVERLAP_M4y5z_S_vrr = PAY*I_TWOBODYOVERLAP_L3y5z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    Double I_TWOBODYOVERLAP_M3y6z_S_vrr = PAY*I_TWOBODYOVERLAP_L2y6z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    Double I_TWOBODYOVERLAP_M2y7z_S_vrr = PAY*I_TWOBODYOVERLAP_Ly7z_S_vrr+oned2z*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_My8z_S_vrr = PAY*I_TWOBODYOVERLAP_L8z_S_vrr;
    Double I_TWOBODYOVERLAP_M9z_S_vrr = PAZ*I_TWOBODYOVERLAP_L8z_S_vrr+8*oned2z*I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_M_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_M9x_S += I_TWOBODYOVERLAP_M9x_S_vrr;
    I_TWOBODYOVERLAP_M8xy_S += I_TWOBODYOVERLAP_M8xy_S_vrr;
    I_TWOBODYOVERLAP_M8xz_S += I_TWOBODYOVERLAP_M8xz_S_vrr;
    I_TWOBODYOVERLAP_M7x2y_S += I_TWOBODYOVERLAP_M7x2y_S_vrr;
    I_TWOBODYOVERLAP_M7xyz_S += I_TWOBODYOVERLAP_M7xyz_S_vrr;
    I_TWOBODYOVERLAP_M7x2z_S += I_TWOBODYOVERLAP_M7x2z_S_vrr;
    I_TWOBODYOVERLAP_M6x3y_S += I_TWOBODYOVERLAP_M6x3y_S_vrr;
    I_TWOBODYOVERLAP_M6x2yz_S += I_TWOBODYOVERLAP_M6x2yz_S_vrr;
    I_TWOBODYOVERLAP_M6xy2z_S += I_TWOBODYOVERLAP_M6xy2z_S_vrr;
    I_TWOBODYOVERLAP_M6x3z_S += I_TWOBODYOVERLAP_M6x3z_S_vrr;
    I_TWOBODYOVERLAP_M5x4y_S += I_TWOBODYOVERLAP_M5x4y_S_vrr;
    I_TWOBODYOVERLAP_M5x3yz_S += I_TWOBODYOVERLAP_M5x3yz_S_vrr;
    I_TWOBODYOVERLAP_M5x2y2z_S += I_TWOBODYOVERLAP_M5x2y2z_S_vrr;
    I_TWOBODYOVERLAP_M5xy3z_S += I_TWOBODYOVERLAP_M5xy3z_S_vrr;
    I_TWOBODYOVERLAP_M5x4z_S += I_TWOBODYOVERLAP_M5x4z_S_vrr;
    I_TWOBODYOVERLAP_M4x5y_S += I_TWOBODYOVERLAP_M4x5y_S_vrr;
    I_TWOBODYOVERLAP_M4x4yz_S += I_TWOBODYOVERLAP_M4x4yz_S_vrr;
    I_TWOBODYOVERLAP_M4x3y2z_S += I_TWOBODYOVERLAP_M4x3y2z_S_vrr;
    I_TWOBODYOVERLAP_M4x2y3z_S += I_TWOBODYOVERLAP_M4x2y3z_S_vrr;
    I_TWOBODYOVERLAP_M4xy4z_S += I_TWOBODYOVERLAP_M4xy4z_S_vrr;
    I_TWOBODYOVERLAP_M4x5z_S += I_TWOBODYOVERLAP_M4x5z_S_vrr;
    I_TWOBODYOVERLAP_M3x6y_S += I_TWOBODYOVERLAP_M3x6y_S_vrr;
    I_TWOBODYOVERLAP_M3x5yz_S += I_TWOBODYOVERLAP_M3x5yz_S_vrr;
    I_TWOBODYOVERLAP_M3x4y2z_S += I_TWOBODYOVERLAP_M3x4y2z_S_vrr;
    I_TWOBODYOVERLAP_M3x3y3z_S += I_TWOBODYOVERLAP_M3x3y3z_S_vrr;
    I_TWOBODYOVERLAP_M3x2y4z_S += I_TWOBODYOVERLAP_M3x2y4z_S_vrr;
    I_TWOBODYOVERLAP_M3xy5z_S += I_TWOBODYOVERLAP_M3xy5z_S_vrr;
    I_TWOBODYOVERLAP_M3x6z_S += I_TWOBODYOVERLAP_M3x6z_S_vrr;
    I_TWOBODYOVERLAP_M2x7y_S += I_TWOBODYOVERLAP_M2x7y_S_vrr;
    I_TWOBODYOVERLAP_M2x6yz_S += I_TWOBODYOVERLAP_M2x6yz_S_vrr;
    I_TWOBODYOVERLAP_M2x5y2z_S += I_TWOBODYOVERLAP_M2x5y2z_S_vrr;
    I_TWOBODYOVERLAP_M2x4y3z_S += I_TWOBODYOVERLAP_M2x4y3z_S_vrr;
    I_TWOBODYOVERLAP_M2x3y4z_S += I_TWOBODYOVERLAP_M2x3y4z_S_vrr;
    I_TWOBODYOVERLAP_M2x2y5z_S += I_TWOBODYOVERLAP_M2x2y5z_S_vrr;
    I_TWOBODYOVERLAP_M2xy6z_S += I_TWOBODYOVERLAP_M2xy6z_S_vrr;
    I_TWOBODYOVERLAP_M2x7z_S += I_TWOBODYOVERLAP_M2x7z_S_vrr;
    I_TWOBODYOVERLAP_Mx8y_S += I_TWOBODYOVERLAP_Mx8y_S_vrr;
    I_TWOBODYOVERLAP_Mx7yz_S += I_TWOBODYOVERLAP_Mx7yz_S_vrr;
    I_TWOBODYOVERLAP_Mx6y2z_S += I_TWOBODYOVERLAP_Mx6y2z_S_vrr;
    I_TWOBODYOVERLAP_Mx5y3z_S += I_TWOBODYOVERLAP_Mx5y3z_S_vrr;
    I_TWOBODYOVERLAP_Mx4y4z_S += I_TWOBODYOVERLAP_Mx4y4z_S_vrr;
    I_TWOBODYOVERLAP_Mx3y5z_S += I_TWOBODYOVERLAP_Mx3y5z_S_vrr;
    I_TWOBODYOVERLAP_Mx2y6z_S += I_TWOBODYOVERLAP_Mx2y6z_S_vrr;
    I_TWOBODYOVERLAP_Mxy7z_S += I_TWOBODYOVERLAP_Mxy7z_S_vrr;
    I_TWOBODYOVERLAP_Mx8z_S += I_TWOBODYOVERLAP_Mx8z_S_vrr;
    I_TWOBODYOVERLAP_M9y_S += I_TWOBODYOVERLAP_M9y_S_vrr;
    I_TWOBODYOVERLAP_M8yz_S += I_TWOBODYOVERLAP_M8yz_S_vrr;
    I_TWOBODYOVERLAP_M7y2z_S += I_TWOBODYOVERLAP_M7y2z_S_vrr;
    I_TWOBODYOVERLAP_M6y3z_S += I_TWOBODYOVERLAP_M6y3z_S_vrr;
    I_TWOBODYOVERLAP_M5y4z_S += I_TWOBODYOVERLAP_M5y4z_S_vrr;
    I_TWOBODYOVERLAP_M4y5z_S += I_TWOBODYOVERLAP_M4y5z_S_vrr;
    I_TWOBODYOVERLAP_M3y6z_S += I_TWOBODYOVERLAP_M3y6z_S_vrr;
    I_TWOBODYOVERLAP_M2y7z_S += I_TWOBODYOVERLAP_M2y7z_S_vrr;
    I_TWOBODYOVERLAP_My8z_S += I_TWOBODYOVERLAP_My8z_S_vrr;
    I_TWOBODYOVERLAP_M9z_S += I_TWOBODYOVERLAP_M9z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_L_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_L8x_S += I_TWOBODYOVERLAP_L8x_S_vrr;
    I_TWOBODYOVERLAP_L7xy_S += I_TWOBODYOVERLAP_L7xy_S_vrr;
    I_TWOBODYOVERLAP_L7xz_S += I_TWOBODYOVERLAP_L7xz_S_vrr;
    I_TWOBODYOVERLAP_L6x2y_S += I_TWOBODYOVERLAP_L6x2y_S_vrr;
    I_TWOBODYOVERLAP_L6xyz_S += I_TWOBODYOVERLAP_L6xyz_S_vrr;
    I_TWOBODYOVERLAP_L6x2z_S += I_TWOBODYOVERLAP_L6x2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3y_S += I_TWOBODYOVERLAP_L5x3y_S_vrr;
    I_TWOBODYOVERLAP_L5x2yz_S += I_TWOBODYOVERLAP_L5x2yz_S_vrr;
    I_TWOBODYOVERLAP_L5xy2z_S += I_TWOBODYOVERLAP_L5xy2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3z_S += I_TWOBODYOVERLAP_L5x3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4y_S += I_TWOBODYOVERLAP_L4x4y_S_vrr;
    I_TWOBODYOVERLAP_L4x3yz_S += I_TWOBODYOVERLAP_L4x3yz_S_vrr;
    I_TWOBODYOVERLAP_L4x2y2z_S += I_TWOBODYOVERLAP_L4x2y2z_S_vrr;
    I_TWOBODYOVERLAP_L4xy3z_S += I_TWOBODYOVERLAP_L4xy3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4z_S += I_TWOBODYOVERLAP_L4x4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5y_S += I_TWOBODYOVERLAP_L3x5y_S_vrr;
    I_TWOBODYOVERLAP_L3x4yz_S += I_TWOBODYOVERLAP_L3x4yz_S_vrr;
    I_TWOBODYOVERLAP_L3x3y2z_S += I_TWOBODYOVERLAP_L3x3y2z_S_vrr;
    I_TWOBODYOVERLAP_L3x2y3z_S += I_TWOBODYOVERLAP_L3x2y3z_S_vrr;
    I_TWOBODYOVERLAP_L3xy4z_S += I_TWOBODYOVERLAP_L3xy4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5z_S += I_TWOBODYOVERLAP_L3x5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6y_S += I_TWOBODYOVERLAP_L2x6y_S_vrr;
    I_TWOBODYOVERLAP_L2x5yz_S += I_TWOBODYOVERLAP_L2x5yz_S_vrr;
    I_TWOBODYOVERLAP_L2x4y2z_S += I_TWOBODYOVERLAP_L2x4y2z_S_vrr;
    I_TWOBODYOVERLAP_L2x3y3z_S += I_TWOBODYOVERLAP_L2x3y3z_S_vrr;
    I_TWOBODYOVERLAP_L2x2y4z_S += I_TWOBODYOVERLAP_L2x2y4z_S_vrr;
    I_TWOBODYOVERLAP_L2xy5z_S += I_TWOBODYOVERLAP_L2xy5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6z_S += I_TWOBODYOVERLAP_L2x6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7y_S += I_TWOBODYOVERLAP_Lx7y_S_vrr;
    I_TWOBODYOVERLAP_Lx6yz_S += I_TWOBODYOVERLAP_Lx6yz_S_vrr;
    I_TWOBODYOVERLAP_Lx5y2z_S += I_TWOBODYOVERLAP_Lx5y2z_S_vrr;
    I_TWOBODYOVERLAP_Lx4y3z_S += I_TWOBODYOVERLAP_Lx4y3z_S_vrr;
    I_TWOBODYOVERLAP_Lx3y4z_S += I_TWOBODYOVERLAP_Lx3y4z_S_vrr;
    I_TWOBODYOVERLAP_Lx2y5z_S += I_TWOBODYOVERLAP_Lx2y5z_S_vrr;
    I_TWOBODYOVERLAP_Lxy6z_S += I_TWOBODYOVERLAP_Lxy6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7z_S += I_TWOBODYOVERLAP_Lx7z_S_vrr;
    I_TWOBODYOVERLAP_L8y_S += I_TWOBODYOVERLAP_L8y_S_vrr;
    I_TWOBODYOVERLAP_L7yz_S += I_TWOBODYOVERLAP_L7yz_S_vrr;
    I_TWOBODYOVERLAP_L6y2z_S += I_TWOBODYOVERLAP_L6y2z_S_vrr;
    I_TWOBODYOVERLAP_L5y3z_S += I_TWOBODYOVERLAP_L5y3z_S_vrr;
    I_TWOBODYOVERLAP_L4y4z_S += I_TWOBODYOVERLAP_L4y4z_S_vrr;
    I_TWOBODYOVERLAP_L3y5z_S += I_TWOBODYOVERLAP_L3y5z_S_vrr;
    I_TWOBODYOVERLAP_L2y6z_S += I_TWOBODYOVERLAP_L2y6z_S_vrr;
    I_TWOBODYOVERLAP_Ly7z_S += I_TWOBODYOVERLAP_Ly7z_S_vrr;
    I_TWOBODYOVERLAP_L8z_S += I_TWOBODYOVERLAP_L8z_S_vrr;

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
   * shell quartet name: SQ_TWOBODYOVERLAP_I_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
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
  Double I_TWOBODYOVERLAP_I6y_Px = I_TWOBODYOVERLAP_Kx6y_S+ABX*I_TWOBODYOVERLAP_I6y_S;
  Double I_TWOBODYOVERLAP_I5yz_Px = I_TWOBODYOVERLAP_Kx5yz_S+ABX*I_TWOBODYOVERLAP_I5yz_S;
  Double I_TWOBODYOVERLAP_I4y2z_Px = I_TWOBODYOVERLAP_Kx4y2z_S+ABX*I_TWOBODYOVERLAP_I4y2z_S;
  Double I_TWOBODYOVERLAP_I3y3z_Px = I_TWOBODYOVERLAP_Kx3y3z_S+ABX*I_TWOBODYOVERLAP_I3y3z_S;
  Double I_TWOBODYOVERLAP_I2y4z_Px = I_TWOBODYOVERLAP_Kx2y4z_S+ABX*I_TWOBODYOVERLAP_I2y4z_S;
  Double I_TWOBODYOVERLAP_Iy5z_Px = I_TWOBODYOVERLAP_Kxy5z_S+ABX*I_TWOBODYOVERLAP_Iy5z_S;
  Double I_TWOBODYOVERLAP_I6z_Px = I_TWOBODYOVERLAP_Kx6z_S+ABX*I_TWOBODYOVERLAP_I6z_S;
  Double I_TWOBODYOVERLAP_I6x_Py = I_TWOBODYOVERLAP_K6xy_S+ABY*I_TWOBODYOVERLAP_I6x_S;
  Double I_TWOBODYOVERLAP_I5xy_Py = I_TWOBODYOVERLAP_K5x2y_S+ABY*I_TWOBODYOVERLAP_I5xy_S;
  Double I_TWOBODYOVERLAP_I5xz_Py = I_TWOBODYOVERLAP_K5xyz_S+ABY*I_TWOBODYOVERLAP_I5xz_S;
  Double I_TWOBODYOVERLAP_I4x2y_Py = I_TWOBODYOVERLAP_K4x3y_S+ABY*I_TWOBODYOVERLAP_I4x2y_S;
  Double I_TWOBODYOVERLAP_I4xyz_Py = I_TWOBODYOVERLAP_K4x2yz_S+ABY*I_TWOBODYOVERLAP_I4xyz_S;
  Double I_TWOBODYOVERLAP_I4x2z_Py = I_TWOBODYOVERLAP_K4xy2z_S+ABY*I_TWOBODYOVERLAP_I4x2z_S;
  Double I_TWOBODYOVERLAP_I3x3y_Py = I_TWOBODYOVERLAP_K3x4y_S+ABY*I_TWOBODYOVERLAP_I3x3y_S;
  Double I_TWOBODYOVERLAP_I3x2yz_Py = I_TWOBODYOVERLAP_K3x3yz_S+ABY*I_TWOBODYOVERLAP_I3x2yz_S;
  Double I_TWOBODYOVERLAP_I3xy2z_Py = I_TWOBODYOVERLAP_K3x2y2z_S+ABY*I_TWOBODYOVERLAP_I3xy2z_S;
  Double I_TWOBODYOVERLAP_I3x3z_Py = I_TWOBODYOVERLAP_K3xy3z_S+ABY*I_TWOBODYOVERLAP_I3x3z_S;
  Double I_TWOBODYOVERLAP_I2x4y_Py = I_TWOBODYOVERLAP_K2x5y_S+ABY*I_TWOBODYOVERLAP_I2x4y_S;
  Double I_TWOBODYOVERLAP_I2x3yz_Py = I_TWOBODYOVERLAP_K2x4yz_S+ABY*I_TWOBODYOVERLAP_I2x3yz_S;
  Double I_TWOBODYOVERLAP_I2x2y2z_Py = I_TWOBODYOVERLAP_K2x3y2z_S+ABY*I_TWOBODYOVERLAP_I2x2y2z_S;
  Double I_TWOBODYOVERLAP_I2xy3z_Py = I_TWOBODYOVERLAP_K2x2y3z_S+ABY*I_TWOBODYOVERLAP_I2xy3z_S;
  Double I_TWOBODYOVERLAP_I2x4z_Py = I_TWOBODYOVERLAP_K2xy4z_S+ABY*I_TWOBODYOVERLAP_I2x4z_S;
  Double I_TWOBODYOVERLAP_Ix5y_Py = I_TWOBODYOVERLAP_Kx6y_S+ABY*I_TWOBODYOVERLAP_Ix5y_S;
  Double I_TWOBODYOVERLAP_Ix4yz_Py = I_TWOBODYOVERLAP_Kx5yz_S+ABY*I_TWOBODYOVERLAP_Ix4yz_S;
  Double I_TWOBODYOVERLAP_Ix3y2z_Py = I_TWOBODYOVERLAP_Kx4y2z_S+ABY*I_TWOBODYOVERLAP_Ix3y2z_S;
  Double I_TWOBODYOVERLAP_Ix2y3z_Py = I_TWOBODYOVERLAP_Kx3y3z_S+ABY*I_TWOBODYOVERLAP_Ix2y3z_S;
  Double I_TWOBODYOVERLAP_Ixy4z_Py = I_TWOBODYOVERLAP_Kx2y4z_S+ABY*I_TWOBODYOVERLAP_Ixy4z_S;
  Double I_TWOBODYOVERLAP_Ix5z_Py = I_TWOBODYOVERLAP_Kxy5z_S+ABY*I_TWOBODYOVERLAP_Ix5z_S;
  Double I_TWOBODYOVERLAP_I6y_Py = I_TWOBODYOVERLAP_K7y_S+ABY*I_TWOBODYOVERLAP_I6y_S;
  Double I_TWOBODYOVERLAP_I5yz_Py = I_TWOBODYOVERLAP_K6yz_S+ABY*I_TWOBODYOVERLAP_I5yz_S;
  Double I_TWOBODYOVERLAP_I4y2z_Py = I_TWOBODYOVERLAP_K5y2z_S+ABY*I_TWOBODYOVERLAP_I4y2z_S;
  Double I_TWOBODYOVERLAP_I3y3z_Py = I_TWOBODYOVERLAP_K4y3z_S+ABY*I_TWOBODYOVERLAP_I3y3z_S;
  Double I_TWOBODYOVERLAP_I2y4z_Py = I_TWOBODYOVERLAP_K3y4z_S+ABY*I_TWOBODYOVERLAP_I2y4z_S;
  Double I_TWOBODYOVERLAP_Iy5z_Py = I_TWOBODYOVERLAP_K2y5z_S+ABY*I_TWOBODYOVERLAP_Iy5z_S;
  Double I_TWOBODYOVERLAP_I6z_Py = I_TWOBODYOVERLAP_Ky6z_S+ABY*I_TWOBODYOVERLAP_I6z_S;
  Double I_TWOBODYOVERLAP_I6x_Pz = I_TWOBODYOVERLAP_K6xz_S+ABZ*I_TWOBODYOVERLAP_I6x_S;
  Double I_TWOBODYOVERLAP_I5xy_Pz = I_TWOBODYOVERLAP_K5xyz_S+ABZ*I_TWOBODYOVERLAP_I5xy_S;
  Double I_TWOBODYOVERLAP_I5xz_Pz = I_TWOBODYOVERLAP_K5x2z_S+ABZ*I_TWOBODYOVERLAP_I5xz_S;
  Double I_TWOBODYOVERLAP_I4x2y_Pz = I_TWOBODYOVERLAP_K4x2yz_S+ABZ*I_TWOBODYOVERLAP_I4x2y_S;
  Double I_TWOBODYOVERLAP_I4xyz_Pz = I_TWOBODYOVERLAP_K4xy2z_S+ABZ*I_TWOBODYOVERLAP_I4xyz_S;
  Double I_TWOBODYOVERLAP_I4x2z_Pz = I_TWOBODYOVERLAP_K4x3z_S+ABZ*I_TWOBODYOVERLAP_I4x2z_S;
  Double I_TWOBODYOVERLAP_I3x3y_Pz = I_TWOBODYOVERLAP_K3x3yz_S+ABZ*I_TWOBODYOVERLAP_I3x3y_S;
  Double I_TWOBODYOVERLAP_I3x2yz_Pz = I_TWOBODYOVERLAP_K3x2y2z_S+ABZ*I_TWOBODYOVERLAP_I3x2yz_S;
  Double I_TWOBODYOVERLAP_I3xy2z_Pz = I_TWOBODYOVERLAP_K3xy3z_S+ABZ*I_TWOBODYOVERLAP_I3xy2z_S;
  Double I_TWOBODYOVERLAP_I3x3z_Pz = I_TWOBODYOVERLAP_K3x4z_S+ABZ*I_TWOBODYOVERLAP_I3x3z_S;
  Double I_TWOBODYOVERLAP_I2x4y_Pz = I_TWOBODYOVERLAP_K2x4yz_S+ABZ*I_TWOBODYOVERLAP_I2x4y_S;
  Double I_TWOBODYOVERLAP_I2x3yz_Pz = I_TWOBODYOVERLAP_K2x3y2z_S+ABZ*I_TWOBODYOVERLAP_I2x3yz_S;
  Double I_TWOBODYOVERLAP_I2x2y2z_Pz = I_TWOBODYOVERLAP_K2x2y3z_S+ABZ*I_TWOBODYOVERLAP_I2x2y2z_S;
  Double I_TWOBODYOVERLAP_I2xy3z_Pz = I_TWOBODYOVERLAP_K2xy4z_S+ABZ*I_TWOBODYOVERLAP_I2xy3z_S;
  Double I_TWOBODYOVERLAP_I2x4z_Pz = I_TWOBODYOVERLAP_K2x5z_S+ABZ*I_TWOBODYOVERLAP_I2x4z_S;
  Double I_TWOBODYOVERLAP_Ix5y_Pz = I_TWOBODYOVERLAP_Kx5yz_S+ABZ*I_TWOBODYOVERLAP_Ix5y_S;
  Double I_TWOBODYOVERLAP_Ix4yz_Pz = I_TWOBODYOVERLAP_Kx4y2z_S+ABZ*I_TWOBODYOVERLAP_Ix4yz_S;
  Double I_TWOBODYOVERLAP_Ix3y2z_Pz = I_TWOBODYOVERLAP_Kx3y3z_S+ABZ*I_TWOBODYOVERLAP_Ix3y2z_S;
  Double I_TWOBODYOVERLAP_Ix2y3z_Pz = I_TWOBODYOVERLAP_Kx2y4z_S+ABZ*I_TWOBODYOVERLAP_Ix2y3z_S;
  Double I_TWOBODYOVERLAP_Ixy4z_Pz = I_TWOBODYOVERLAP_Kxy5z_S+ABZ*I_TWOBODYOVERLAP_Ixy4z_S;
  Double I_TWOBODYOVERLAP_Ix5z_Pz = I_TWOBODYOVERLAP_Kx6z_S+ABZ*I_TWOBODYOVERLAP_Ix5z_S;
  Double I_TWOBODYOVERLAP_I6y_Pz = I_TWOBODYOVERLAP_K6yz_S+ABZ*I_TWOBODYOVERLAP_I6y_S;
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
   * totally 63 integrals are omitted 
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
   * shell quartet name: SQ_TWOBODYOVERLAP_K_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_K7x_Px = I_TWOBODYOVERLAP_L8x_S+ABX*I_TWOBODYOVERLAP_K7x_S;
  Double I_TWOBODYOVERLAP_K6xy_Px = I_TWOBODYOVERLAP_L7xy_S+ABX*I_TWOBODYOVERLAP_K6xy_S;
  Double I_TWOBODYOVERLAP_K6xz_Px = I_TWOBODYOVERLAP_L7xz_S+ABX*I_TWOBODYOVERLAP_K6xz_S;
  Double I_TWOBODYOVERLAP_K5x2y_Px = I_TWOBODYOVERLAP_L6x2y_S+ABX*I_TWOBODYOVERLAP_K5x2y_S;
  Double I_TWOBODYOVERLAP_K5xyz_Px = I_TWOBODYOVERLAP_L6xyz_S+ABX*I_TWOBODYOVERLAP_K5xyz_S;
  Double I_TWOBODYOVERLAP_K5x2z_Px = I_TWOBODYOVERLAP_L6x2z_S+ABX*I_TWOBODYOVERLAP_K5x2z_S;
  Double I_TWOBODYOVERLAP_K4x3y_Px = I_TWOBODYOVERLAP_L5x3y_S+ABX*I_TWOBODYOVERLAP_K4x3y_S;
  Double I_TWOBODYOVERLAP_K4x2yz_Px = I_TWOBODYOVERLAP_L5x2yz_S+ABX*I_TWOBODYOVERLAP_K4x2yz_S;
  Double I_TWOBODYOVERLAP_K4xy2z_Px = I_TWOBODYOVERLAP_L5xy2z_S+ABX*I_TWOBODYOVERLAP_K4xy2z_S;
  Double I_TWOBODYOVERLAP_K4x3z_Px = I_TWOBODYOVERLAP_L5x3z_S+ABX*I_TWOBODYOVERLAP_K4x3z_S;
  Double I_TWOBODYOVERLAP_K3x4y_Px = I_TWOBODYOVERLAP_L4x4y_S+ABX*I_TWOBODYOVERLAP_K3x4y_S;
  Double I_TWOBODYOVERLAP_K3x3yz_Px = I_TWOBODYOVERLAP_L4x3yz_S+ABX*I_TWOBODYOVERLAP_K3x3yz_S;
  Double I_TWOBODYOVERLAP_K3x2y2z_Px = I_TWOBODYOVERLAP_L4x2y2z_S+ABX*I_TWOBODYOVERLAP_K3x2y2z_S;
  Double I_TWOBODYOVERLAP_K3xy3z_Px = I_TWOBODYOVERLAP_L4xy3z_S+ABX*I_TWOBODYOVERLAP_K3xy3z_S;
  Double I_TWOBODYOVERLAP_K3x4z_Px = I_TWOBODYOVERLAP_L4x4z_S+ABX*I_TWOBODYOVERLAP_K3x4z_S;
  Double I_TWOBODYOVERLAP_K2x5y_Px = I_TWOBODYOVERLAP_L3x5y_S+ABX*I_TWOBODYOVERLAP_K2x5y_S;
  Double I_TWOBODYOVERLAP_K2x4yz_Px = I_TWOBODYOVERLAP_L3x4yz_S+ABX*I_TWOBODYOVERLAP_K2x4yz_S;
  Double I_TWOBODYOVERLAP_K2x3y2z_Px = I_TWOBODYOVERLAP_L3x3y2z_S+ABX*I_TWOBODYOVERLAP_K2x3y2z_S;
  Double I_TWOBODYOVERLAP_K2x2y3z_Px = I_TWOBODYOVERLAP_L3x2y3z_S+ABX*I_TWOBODYOVERLAP_K2x2y3z_S;
  Double I_TWOBODYOVERLAP_K2xy4z_Px = I_TWOBODYOVERLAP_L3xy4z_S+ABX*I_TWOBODYOVERLAP_K2xy4z_S;
  Double I_TWOBODYOVERLAP_K2x5z_Px = I_TWOBODYOVERLAP_L3x5z_S+ABX*I_TWOBODYOVERLAP_K2x5z_S;
  Double I_TWOBODYOVERLAP_Kx6y_Px = I_TWOBODYOVERLAP_L2x6y_S+ABX*I_TWOBODYOVERLAP_Kx6y_S;
  Double I_TWOBODYOVERLAP_Kx5yz_Px = I_TWOBODYOVERLAP_L2x5yz_S+ABX*I_TWOBODYOVERLAP_Kx5yz_S;
  Double I_TWOBODYOVERLAP_Kx4y2z_Px = I_TWOBODYOVERLAP_L2x4y2z_S+ABX*I_TWOBODYOVERLAP_Kx4y2z_S;
  Double I_TWOBODYOVERLAP_Kx3y3z_Px = I_TWOBODYOVERLAP_L2x3y3z_S+ABX*I_TWOBODYOVERLAP_Kx3y3z_S;
  Double I_TWOBODYOVERLAP_Kx2y4z_Px = I_TWOBODYOVERLAP_L2x2y4z_S+ABX*I_TWOBODYOVERLAP_Kx2y4z_S;
  Double I_TWOBODYOVERLAP_Kxy5z_Px = I_TWOBODYOVERLAP_L2xy5z_S+ABX*I_TWOBODYOVERLAP_Kxy5z_S;
  Double I_TWOBODYOVERLAP_Kx6z_Px = I_TWOBODYOVERLAP_L2x6z_S+ABX*I_TWOBODYOVERLAP_Kx6z_S;
  Double I_TWOBODYOVERLAP_K7y_Px = I_TWOBODYOVERLAP_Lx7y_S+ABX*I_TWOBODYOVERLAP_K7y_S;
  Double I_TWOBODYOVERLAP_K6yz_Px = I_TWOBODYOVERLAP_Lx6yz_S+ABX*I_TWOBODYOVERLAP_K6yz_S;
  Double I_TWOBODYOVERLAP_K5y2z_Px = I_TWOBODYOVERLAP_Lx5y2z_S+ABX*I_TWOBODYOVERLAP_K5y2z_S;
  Double I_TWOBODYOVERLAP_K4y3z_Px = I_TWOBODYOVERLAP_Lx4y3z_S+ABX*I_TWOBODYOVERLAP_K4y3z_S;
  Double I_TWOBODYOVERLAP_K3y4z_Px = I_TWOBODYOVERLAP_Lx3y4z_S+ABX*I_TWOBODYOVERLAP_K3y4z_S;
  Double I_TWOBODYOVERLAP_K2y5z_Px = I_TWOBODYOVERLAP_Lx2y5z_S+ABX*I_TWOBODYOVERLAP_K2y5z_S;
  Double I_TWOBODYOVERLAP_Ky6z_Px = I_TWOBODYOVERLAP_Lxy6z_S+ABX*I_TWOBODYOVERLAP_Ky6z_S;
  Double I_TWOBODYOVERLAP_K7z_Px = I_TWOBODYOVERLAP_Lx7z_S+ABX*I_TWOBODYOVERLAP_K7z_S;
  Double I_TWOBODYOVERLAP_K6xy_Py = I_TWOBODYOVERLAP_L6x2y_S+ABY*I_TWOBODYOVERLAP_K6xy_S;
  Double I_TWOBODYOVERLAP_K6xz_Py = I_TWOBODYOVERLAP_L6xyz_S+ABY*I_TWOBODYOVERLAP_K6xz_S;
  Double I_TWOBODYOVERLAP_K5x2y_Py = I_TWOBODYOVERLAP_L5x3y_S+ABY*I_TWOBODYOVERLAP_K5x2y_S;
  Double I_TWOBODYOVERLAP_K5xyz_Py = I_TWOBODYOVERLAP_L5x2yz_S+ABY*I_TWOBODYOVERLAP_K5xyz_S;
  Double I_TWOBODYOVERLAP_K5x2z_Py = I_TWOBODYOVERLAP_L5xy2z_S+ABY*I_TWOBODYOVERLAP_K5x2z_S;
  Double I_TWOBODYOVERLAP_K4x3y_Py = I_TWOBODYOVERLAP_L4x4y_S+ABY*I_TWOBODYOVERLAP_K4x3y_S;
  Double I_TWOBODYOVERLAP_K4x2yz_Py = I_TWOBODYOVERLAP_L4x3yz_S+ABY*I_TWOBODYOVERLAP_K4x2yz_S;
  Double I_TWOBODYOVERLAP_K4xy2z_Py = I_TWOBODYOVERLAP_L4x2y2z_S+ABY*I_TWOBODYOVERLAP_K4xy2z_S;
  Double I_TWOBODYOVERLAP_K4x3z_Py = I_TWOBODYOVERLAP_L4xy3z_S+ABY*I_TWOBODYOVERLAP_K4x3z_S;
  Double I_TWOBODYOVERLAP_K3x4y_Py = I_TWOBODYOVERLAP_L3x5y_S+ABY*I_TWOBODYOVERLAP_K3x4y_S;
  Double I_TWOBODYOVERLAP_K3x3yz_Py = I_TWOBODYOVERLAP_L3x4yz_S+ABY*I_TWOBODYOVERLAP_K3x3yz_S;
  Double I_TWOBODYOVERLAP_K3x2y2z_Py = I_TWOBODYOVERLAP_L3x3y2z_S+ABY*I_TWOBODYOVERLAP_K3x2y2z_S;
  Double I_TWOBODYOVERLAP_K3xy3z_Py = I_TWOBODYOVERLAP_L3x2y3z_S+ABY*I_TWOBODYOVERLAP_K3xy3z_S;
  Double I_TWOBODYOVERLAP_K3x4z_Py = I_TWOBODYOVERLAP_L3xy4z_S+ABY*I_TWOBODYOVERLAP_K3x4z_S;
  Double I_TWOBODYOVERLAP_K2x5y_Py = I_TWOBODYOVERLAP_L2x6y_S+ABY*I_TWOBODYOVERLAP_K2x5y_S;
  Double I_TWOBODYOVERLAP_K2x4yz_Py = I_TWOBODYOVERLAP_L2x5yz_S+ABY*I_TWOBODYOVERLAP_K2x4yz_S;
  Double I_TWOBODYOVERLAP_K2x3y2z_Py = I_TWOBODYOVERLAP_L2x4y2z_S+ABY*I_TWOBODYOVERLAP_K2x3y2z_S;
  Double I_TWOBODYOVERLAP_K2x2y3z_Py = I_TWOBODYOVERLAP_L2x3y3z_S+ABY*I_TWOBODYOVERLAP_K2x2y3z_S;
  Double I_TWOBODYOVERLAP_K2xy4z_Py = I_TWOBODYOVERLAP_L2x2y4z_S+ABY*I_TWOBODYOVERLAP_K2xy4z_S;
  Double I_TWOBODYOVERLAP_K2x5z_Py = I_TWOBODYOVERLAP_L2xy5z_S+ABY*I_TWOBODYOVERLAP_K2x5z_S;
  Double I_TWOBODYOVERLAP_Kx6y_Py = I_TWOBODYOVERLAP_Lx7y_S+ABY*I_TWOBODYOVERLAP_Kx6y_S;
  Double I_TWOBODYOVERLAP_Kx5yz_Py = I_TWOBODYOVERLAP_Lx6yz_S+ABY*I_TWOBODYOVERLAP_Kx5yz_S;
  Double I_TWOBODYOVERLAP_Kx4y2z_Py = I_TWOBODYOVERLAP_Lx5y2z_S+ABY*I_TWOBODYOVERLAP_Kx4y2z_S;
  Double I_TWOBODYOVERLAP_Kx3y3z_Py = I_TWOBODYOVERLAP_Lx4y3z_S+ABY*I_TWOBODYOVERLAP_Kx3y3z_S;
  Double I_TWOBODYOVERLAP_Kx2y4z_Py = I_TWOBODYOVERLAP_Lx3y4z_S+ABY*I_TWOBODYOVERLAP_Kx2y4z_S;
  Double I_TWOBODYOVERLAP_Kxy5z_Py = I_TWOBODYOVERLAP_Lx2y5z_S+ABY*I_TWOBODYOVERLAP_Kxy5z_S;
  Double I_TWOBODYOVERLAP_Kx6z_Py = I_TWOBODYOVERLAP_Lxy6z_S+ABY*I_TWOBODYOVERLAP_Kx6z_S;
  Double I_TWOBODYOVERLAP_K7y_Py = I_TWOBODYOVERLAP_L8y_S+ABY*I_TWOBODYOVERLAP_K7y_S;
  Double I_TWOBODYOVERLAP_K6yz_Py = I_TWOBODYOVERLAP_L7yz_S+ABY*I_TWOBODYOVERLAP_K6yz_S;
  Double I_TWOBODYOVERLAP_K5y2z_Py = I_TWOBODYOVERLAP_L6y2z_S+ABY*I_TWOBODYOVERLAP_K5y2z_S;
  Double I_TWOBODYOVERLAP_K4y3z_Py = I_TWOBODYOVERLAP_L5y3z_S+ABY*I_TWOBODYOVERLAP_K4y3z_S;
  Double I_TWOBODYOVERLAP_K3y4z_Py = I_TWOBODYOVERLAP_L4y4z_S+ABY*I_TWOBODYOVERLAP_K3y4z_S;
  Double I_TWOBODYOVERLAP_K2y5z_Py = I_TWOBODYOVERLAP_L3y5z_S+ABY*I_TWOBODYOVERLAP_K2y5z_S;
  Double I_TWOBODYOVERLAP_Ky6z_Py = I_TWOBODYOVERLAP_L2y6z_S+ABY*I_TWOBODYOVERLAP_Ky6z_S;
  Double I_TWOBODYOVERLAP_K7z_Py = I_TWOBODYOVERLAP_Ly7z_S+ABY*I_TWOBODYOVERLAP_K7z_S;
  Double I_TWOBODYOVERLAP_K6xy_Pz = I_TWOBODYOVERLAP_L6xyz_S+ABZ*I_TWOBODYOVERLAP_K6xy_S;
  Double I_TWOBODYOVERLAP_K6xz_Pz = I_TWOBODYOVERLAP_L6x2z_S+ABZ*I_TWOBODYOVERLAP_K6xz_S;
  Double I_TWOBODYOVERLAP_K5x2y_Pz = I_TWOBODYOVERLAP_L5x2yz_S+ABZ*I_TWOBODYOVERLAP_K5x2y_S;
  Double I_TWOBODYOVERLAP_K5xyz_Pz = I_TWOBODYOVERLAP_L5xy2z_S+ABZ*I_TWOBODYOVERLAP_K5xyz_S;
  Double I_TWOBODYOVERLAP_K5x2z_Pz = I_TWOBODYOVERLAP_L5x3z_S+ABZ*I_TWOBODYOVERLAP_K5x2z_S;
  Double I_TWOBODYOVERLAP_K4x3y_Pz = I_TWOBODYOVERLAP_L4x3yz_S+ABZ*I_TWOBODYOVERLAP_K4x3y_S;
  Double I_TWOBODYOVERLAP_K4x2yz_Pz = I_TWOBODYOVERLAP_L4x2y2z_S+ABZ*I_TWOBODYOVERLAP_K4x2yz_S;
  Double I_TWOBODYOVERLAP_K4xy2z_Pz = I_TWOBODYOVERLAP_L4xy3z_S+ABZ*I_TWOBODYOVERLAP_K4xy2z_S;
  Double I_TWOBODYOVERLAP_K4x3z_Pz = I_TWOBODYOVERLAP_L4x4z_S+ABZ*I_TWOBODYOVERLAP_K4x3z_S;
  Double I_TWOBODYOVERLAP_K3x4y_Pz = I_TWOBODYOVERLAP_L3x4yz_S+ABZ*I_TWOBODYOVERLAP_K3x4y_S;
  Double I_TWOBODYOVERLAP_K3x3yz_Pz = I_TWOBODYOVERLAP_L3x3y2z_S+ABZ*I_TWOBODYOVERLAP_K3x3yz_S;
  Double I_TWOBODYOVERLAP_K3x2y2z_Pz = I_TWOBODYOVERLAP_L3x2y3z_S+ABZ*I_TWOBODYOVERLAP_K3x2y2z_S;
  Double I_TWOBODYOVERLAP_K3xy3z_Pz = I_TWOBODYOVERLAP_L3xy4z_S+ABZ*I_TWOBODYOVERLAP_K3xy3z_S;
  Double I_TWOBODYOVERLAP_K3x4z_Pz = I_TWOBODYOVERLAP_L3x5z_S+ABZ*I_TWOBODYOVERLAP_K3x4z_S;
  Double I_TWOBODYOVERLAP_K2x5y_Pz = I_TWOBODYOVERLAP_L2x5yz_S+ABZ*I_TWOBODYOVERLAP_K2x5y_S;
  Double I_TWOBODYOVERLAP_K2x4yz_Pz = I_TWOBODYOVERLAP_L2x4y2z_S+ABZ*I_TWOBODYOVERLAP_K2x4yz_S;
  Double I_TWOBODYOVERLAP_K2x3y2z_Pz = I_TWOBODYOVERLAP_L2x3y3z_S+ABZ*I_TWOBODYOVERLAP_K2x3y2z_S;
  Double I_TWOBODYOVERLAP_K2x2y3z_Pz = I_TWOBODYOVERLAP_L2x2y4z_S+ABZ*I_TWOBODYOVERLAP_K2x2y3z_S;
  Double I_TWOBODYOVERLAP_K2xy4z_Pz = I_TWOBODYOVERLAP_L2xy5z_S+ABZ*I_TWOBODYOVERLAP_K2xy4z_S;
  Double I_TWOBODYOVERLAP_K2x5z_Pz = I_TWOBODYOVERLAP_L2x6z_S+ABZ*I_TWOBODYOVERLAP_K2x5z_S;
  Double I_TWOBODYOVERLAP_Kx6y_Pz = I_TWOBODYOVERLAP_Lx6yz_S+ABZ*I_TWOBODYOVERLAP_Kx6y_S;
  Double I_TWOBODYOVERLAP_Kx5yz_Pz = I_TWOBODYOVERLAP_Lx5y2z_S+ABZ*I_TWOBODYOVERLAP_Kx5yz_S;
  Double I_TWOBODYOVERLAP_Kx4y2z_Pz = I_TWOBODYOVERLAP_Lx4y3z_S+ABZ*I_TWOBODYOVERLAP_Kx4y2z_S;
  Double I_TWOBODYOVERLAP_Kx3y3z_Pz = I_TWOBODYOVERLAP_Lx3y4z_S+ABZ*I_TWOBODYOVERLAP_Kx3y3z_S;
  Double I_TWOBODYOVERLAP_Kx2y4z_Pz = I_TWOBODYOVERLAP_Lx2y5z_S+ABZ*I_TWOBODYOVERLAP_Kx2y4z_S;
  Double I_TWOBODYOVERLAP_Kxy5z_Pz = I_TWOBODYOVERLAP_Lxy6z_S+ABZ*I_TWOBODYOVERLAP_Kxy5z_S;
  Double I_TWOBODYOVERLAP_Kx6z_Pz = I_TWOBODYOVERLAP_Lx7z_S+ABZ*I_TWOBODYOVERLAP_Kx6z_S;
  Double I_TWOBODYOVERLAP_K6yz_Pz = I_TWOBODYOVERLAP_L6y2z_S+ABZ*I_TWOBODYOVERLAP_K6yz_S;
  Double I_TWOBODYOVERLAP_K5y2z_Pz = I_TWOBODYOVERLAP_L5y3z_S+ABZ*I_TWOBODYOVERLAP_K5y2z_S;
  Double I_TWOBODYOVERLAP_K4y3z_Pz = I_TWOBODYOVERLAP_L4y4z_S+ABZ*I_TWOBODYOVERLAP_K4y3z_S;
  Double I_TWOBODYOVERLAP_K3y4z_Pz = I_TWOBODYOVERLAP_L3y5z_S+ABZ*I_TWOBODYOVERLAP_K3y4z_S;
  Double I_TWOBODYOVERLAP_K2y5z_Pz = I_TWOBODYOVERLAP_L2y6z_S+ABZ*I_TWOBODYOVERLAP_K2y5z_S;
  Double I_TWOBODYOVERLAP_Ky6z_Pz = I_TWOBODYOVERLAP_Ly7z_S+ABZ*I_TWOBODYOVERLAP_Ky6z_S;
  Double I_TWOBODYOVERLAP_K7z_Pz = I_TWOBODYOVERLAP_L8z_S+ABZ*I_TWOBODYOVERLAP_K7z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 84 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_D2x = I_TWOBODYOVERLAP_K7x_Px+ABX*I_TWOBODYOVERLAP_I6x_Px;
  Double I_TWOBODYOVERLAP_I5xy_D2x = I_TWOBODYOVERLAP_K6xy_Px+ABX*I_TWOBODYOVERLAP_I5xy_Px;
  Double I_TWOBODYOVERLAP_I5xz_D2x = I_TWOBODYOVERLAP_K6xz_Px+ABX*I_TWOBODYOVERLAP_I5xz_Px;
  Double I_TWOBODYOVERLAP_I4x2y_D2x = I_TWOBODYOVERLAP_K5x2y_Px+ABX*I_TWOBODYOVERLAP_I4x2y_Px;
  Double I_TWOBODYOVERLAP_I4xyz_D2x = I_TWOBODYOVERLAP_K5xyz_Px+ABX*I_TWOBODYOVERLAP_I4xyz_Px;
  Double I_TWOBODYOVERLAP_I4x2z_D2x = I_TWOBODYOVERLAP_K5x2z_Px+ABX*I_TWOBODYOVERLAP_I4x2z_Px;
  Double I_TWOBODYOVERLAP_I3x3y_D2x = I_TWOBODYOVERLAP_K4x3y_Px+ABX*I_TWOBODYOVERLAP_I3x3y_Px;
  Double I_TWOBODYOVERLAP_I3x2yz_D2x = I_TWOBODYOVERLAP_K4x2yz_Px+ABX*I_TWOBODYOVERLAP_I3x2yz_Px;
  Double I_TWOBODYOVERLAP_I3xy2z_D2x = I_TWOBODYOVERLAP_K4xy2z_Px+ABX*I_TWOBODYOVERLAP_I3xy2z_Px;
  Double I_TWOBODYOVERLAP_I3x3z_D2x = I_TWOBODYOVERLAP_K4x3z_Px+ABX*I_TWOBODYOVERLAP_I3x3z_Px;
  Double I_TWOBODYOVERLAP_I2x4y_D2x = I_TWOBODYOVERLAP_K3x4y_Px+ABX*I_TWOBODYOVERLAP_I2x4y_Px;
  Double I_TWOBODYOVERLAP_I2x3yz_D2x = I_TWOBODYOVERLAP_K3x3yz_Px+ABX*I_TWOBODYOVERLAP_I2x3yz_Px;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2x = I_TWOBODYOVERLAP_K3x2y2z_Px+ABX*I_TWOBODYOVERLAP_I2x2y2z_Px;
  Double I_TWOBODYOVERLAP_I2xy3z_D2x = I_TWOBODYOVERLAP_K3xy3z_Px+ABX*I_TWOBODYOVERLAP_I2xy3z_Px;
  Double I_TWOBODYOVERLAP_I2x4z_D2x = I_TWOBODYOVERLAP_K3x4z_Px+ABX*I_TWOBODYOVERLAP_I2x4z_Px;
  Double I_TWOBODYOVERLAP_Ix5y_D2x = I_TWOBODYOVERLAP_K2x5y_Px+ABX*I_TWOBODYOVERLAP_Ix5y_Px;
  Double I_TWOBODYOVERLAP_Ix4yz_D2x = I_TWOBODYOVERLAP_K2x4yz_Px+ABX*I_TWOBODYOVERLAP_Ix4yz_Px;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2x = I_TWOBODYOVERLAP_K2x3y2z_Px+ABX*I_TWOBODYOVERLAP_Ix3y2z_Px;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2x = I_TWOBODYOVERLAP_K2x2y3z_Px+ABX*I_TWOBODYOVERLAP_Ix2y3z_Px;
  Double I_TWOBODYOVERLAP_Ixy4z_D2x = I_TWOBODYOVERLAP_K2xy4z_Px+ABX*I_TWOBODYOVERLAP_Ixy4z_Px;
  Double I_TWOBODYOVERLAP_Ix5z_D2x = I_TWOBODYOVERLAP_K2x5z_Px+ABX*I_TWOBODYOVERLAP_Ix5z_Px;
  Double I_TWOBODYOVERLAP_I6y_D2x = I_TWOBODYOVERLAP_Kx6y_Px+ABX*I_TWOBODYOVERLAP_I6y_Px;
  Double I_TWOBODYOVERLAP_I5yz_D2x = I_TWOBODYOVERLAP_Kx5yz_Px+ABX*I_TWOBODYOVERLAP_I5yz_Px;
  Double I_TWOBODYOVERLAP_I4y2z_D2x = I_TWOBODYOVERLAP_Kx4y2z_Px+ABX*I_TWOBODYOVERLAP_I4y2z_Px;
  Double I_TWOBODYOVERLAP_I3y3z_D2x = I_TWOBODYOVERLAP_Kx3y3z_Px+ABX*I_TWOBODYOVERLAP_I3y3z_Px;
  Double I_TWOBODYOVERLAP_I2y4z_D2x = I_TWOBODYOVERLAP_Kx2y4z_Px+ABX*I_TWOBODYOVERLAP_I2y4z_Px;
  Double I_TWOBODYOVERLAP_Iy5z_D2x = I_TWOBODYOVERLAP_Kxy5z_Px+ABX*I_TWOBODYOVERLAP_Iy5z_Px;
  Double I_TWOBODYOVERLAP_I6z_D2x = I_TWOBODYOVERLAP_Kx6z_Px+ABX*I_TWOBODYOVERLAP_I6z_Px;
  Double I_TWOBODYOVERLAP_I6x_D2y = I_TWOBODYOVERLAP_K6xy_Py+ABY*I_TWOBODYOVERLAP_I6x_Py;
  Double I_TWOBODYOVERLAP_I5xy_D2y = I_TWOBODYOVERLAP_K5x2y_Py+ABY*I_TWOBODYOVERLAP_I5xy_Py;
  Double I_TWOBODYOVERLAP_I5xz_D2y = I_TWOBODYOVERLAP_K5xyz_Py+ABY*I_TWOBODYOVERLAP_I5xz_Py;
  Double I_TWOBODYOVERLAP_I4x2y_D2y = I_TWOBODYOVERLAP_K4x3y_Py+ABY*I_TWOBODYOVERLAP_I4x2y_Py;
  Double I_TWOBODYOVERLAP_I4xyz_D2y = I_TWOBODYOVERLAP_K4x2yz_Py+ABY*I_TWOBODYOVERLAP_I4xyz_Py;
  Double I_TWOBODYOVERLAP_I4x2z_D2y = I_TWOBODYOVERLAP_K4xy2z_Py+ABY*I_TWOBODYOVERLAP_I4x2z_Py;
  Double I_TWOBODYOVERLAP_I3x3y_D2y = I_TWOBODYOVERLAP_K3x4y_Py+ABY*I_TWOBODYOVERLAP_I3x3y_Py;
  Double I_TWOBODYOVERLAP_I3x2yz_D2y = I_TWOBODYOVERLAP_K3x3yz_Py+ABY*I_TWOBODYOVERLAP_I3x2yz_Py;
  Double I_TWOBODYOVERLAP_I3xy2z_D2y = I_TWOBODYOVERLAP_K3x2y2z_Py+ABY*I_TWOBODYOVERLAP_I3xy2z_Py;
  Double I_TWOBODYOVERLAP_I3x3z_D2y = I_TWOBODYOVERLAP_K3xy3z_Py+ABY*I_TWOBODYOVERLAP_I3x3z_Py;
  Double I_TWOBODYOVERLAP_I2x4y_D2y = I_TWOBODYOVERLAP_K2x5y_Py+ABY*I_TWOBODYOVERLAP_I2x4y_Py;
  Double I_TWOBODYOVERLAP_I2x3yz_D2y = I_TWOBODYOVERLAP_K2x4yz_Py+ABY*I_TWOBODYOVERLAP_I2x3yz_Py;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2y = I_TWOBODYOVERLAP_K2x3y2z_Py+ABY*I_TWOBODYOVERLAP_I2x2y2z_Py;
  Double I_TWOBODYOVERLAP_I2xy3z_D2y = I_TWOBODYOVERLAP_K2x2y3z_Py+ABY*I_TWOBODYOVERLAP_I2xy3z_Py;
  Double I_TWOBODYOVERLAP_I2x4z_D2y = I_TWOBODYOVERLAP_K2xy4z_Py+ABY*I_TWOBODYOVERLAP_I2x4z_Py;
  Double I_TWOBODYOVERLAP_Ix5y_D2y = I_TWOBODYOVERLAP_Kx6y_Py+ABY*I_TWOBODYOVERLAP_Ix5y_Py;
  Double I_TWOBODYOVERLAP_Ix4yz_D2y = I_TWOBODYOVERLAP_Kx5yz_Py+ABY*I_TWOBODYOVERLAP_Ix4yz_Py;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2y = I_TWOBODYOVERLAP_Kx4y2z_Py+ABY*I_TWOBODYOVERLAP_Ix3y2z_Py;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2y = I_TWOBODYOVERLAP_Kx3y3z_Py+ABY*I_TWOBODYOVERLAP_Ix2y3z_Py;
  Double I_TWOBODYOVERLAP_Ixy4z_D2y = I_TWOBODYOVERLAP_Kx2y4z_Py+ABY*I_TWOBODYOVERLAP_Ixy4z_Py;
  Double I_TWOBODYOVERLAP_Ix5z_D2y = I_TWOBODYOVERLAP_Kxy5z_Py+ABY*I_TWOBODYOVERLAP_Ix5z_Py;
  Double I_TWOBODYOVERLAP_I6y_D2y = I_TWOBODYOVERLAP_K7y_Py+ABY*I_TWOBODYOVERLAP_I6y_Py;
  Double I_TWOBODYOVERLAP_I5yz_D2y = I_TWOBODYOVERLAP_K6yz_Py+ABY*I_TWOBODYOVERLAP_I5yz_Py;
  Double I_TWOBODYOVERLAP_I4y2z_D2y = I_TWOBODYOVERLAP_K5y2z_Py+ABY*I_TWOBODYOVERLAP_I4y2z_Py;
  Double I_TWOBODYOVERLAP_I3y3z_D2y = I_TWOBODYOVERLAP_K4y3z_Py+ABY*I_TWOBODYOVERLAP_I3y3z_Py;
  Double I_TWOBODYOVERLAP_I2y4z_D2y = I_TWOBODYOVERLAP_K3y4z_Py+ABY*I_TWOBODYOVERLAP_I2y4z_Py;
  Double I_TWOBODYOVERLAP_Iy5z_D2y = I_TWOBODYOVERLAP_K2y5z_Py+ABY*I_TWOBODYOVERLAP_Iy5z_Py;
  Double I_TWOBODYOVERLAP_I6z_D2y = I_TWOBODYOVERLAP_Ky6z_Py+ABY*I_TWOBODYOVERLAP_I6z_Py;
  Double I_TWOBODYOVERLAP_I6x_D2z = I_TWOBODYOVERLAP_K6xz_Pz+ABZ*I_TWOBODYOVERLAP_I6x_Pz;
  Double I_TWOBODYOVERLAP_I5xy_D2z = I_TWOBODYOVERLAP_K5xyz_Pz+ABZ*I_TWOBODYOVERLAP_I5xy_Pz;
  Double I_TWOBODYOVERLAP_I5xz_D2z = I_TWOBODYOVERLAP_K5x2z_Pz+ABZ*I_TWOBODYOVERLAP_I5xz_Pz;
  Double I_TWOBODYOVERLAP_I4x2y_D2z = I_TWOBODYOVERLAP_K4x2yz_Pz+ABZ*I_TWOBODYOVERLAP_I4x2y_Pz;
  Double I_TWOBODYOVERLAP_I4xyz_D2z = I_TWOBODYOVERLAP_K4xy2z_Pz+ABZ*I_TWOBODYOVERLAP_I4xyz_Pz;
  Double I_TWOBODYOVERLAP_I4x2z_D2z = I_TWOBODYOVERLAP_K4x3z_Pz+ABZ*I_TWOBODYOVERLAP_I4x2z_Pz;
  Double I_TWOBODYOVERLAP_I3x3y_D2z = I_TWOBODYOVERLAP_K3x3yz_Pz+ABZ*I_TWOBODYOVERLAP_I3x3y_Pz;
  Double I_TWOBODYOVERLAP_I3x2yz_D2z = I_TWOBODYOVERLAP_K3x2y2z_Pz+ABZ*I_TWOBODYOVERLAP_I3x2yz_Pz;
  Double I_TWOBODYOVERLAP_I3xy2z_D2z = I_TWOBODYOVERLAP_K3xy3z_Pz+ABZ*I_TWOBODYOVERLAP_I3xy2z_Pz;
  Double I_TWOBODYOVERLAP_I3x3z_D2z = I_TWOBODYOVERLAP_K3x4z_Pz+ABZ*I_TWOBODYOVERLAP_I3x3z_Pz;
  Double I_TWOBODYOVERLAP_I2x4y_D2z = I_TWOBODYOVERLAP_K2x4yz_Pz+ABZ*I_TWOBODYOVERLAP_I2x4y_Pz;
  Double I_TWOBODYOVERLAP_I2x3yz_D2z = I_TWOBODYOVERLAP_K2x3y2z_Pz+ABZ*I_TWOBODYOVERLAP_I2x3yz_Pz;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2z = I_TWOBODYOVERLAP_K2x2y3z_Pz+ABZ*I_TWOBODYOVERLAP_I2x2y2z_Pz;
  Double I_TWOBODYOVERLAP_I2xy3z_D2z = I_TWOBODYOVERLAP_K2xy4z_Pz+ABZ*I_TWOBODYOVERLAP_I2xy3z_Pz;
  Double I_TWOBODYOVERLAP_I2x4z_D2z = I_TWOBODYOVERLAP_K2x5z_Pz+ABZ*I_TWOBODYOVERLAP_I2x4z_Pz;
  Double I_TWOBODYOVERLAP_Ix5y_D2z = I_TWOBODYOVERLAP_Kx5yz_Pz+ABZ*I_TWOBODYOVERLAP_Ix5y_Pz;
  Double I_TWOBODYOVERLAP_Ix4yz_D2z = I_TWOBODYOVERLAP_Kx4y2z_Pz+ABZ*I_TWOBODYOVERLAP_Ix4yz_Pz;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2z = I_TWOBODYOVERLAP_Kx3y3z_Pz+ABZ*I_TWOBODYOVERLAP_Ix3y2z_Pz;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2z = I_TWOBODYOVERLAP_Kx2y4z_Pz+ABZ*I_TWOBODYOVERLAP_Ix2y3z_Pz;
  Double I_TWOBODYOVERLAP_Ixy4z_D2z = I_TWOBODYOVERLAP_Kxy5z_Pz+ABZ*I_TWOBODYOVERLAP_Ixy4z_Pz;
  Double I_TWOBODYOVERLAP_Ix5z_D2z = I_TWOBODYOVERLAP_Kx6z_Pz+ABZ*I_TWOBODYOVERLAP_Ix5z_Pz;
  Double I_TWOBODYOVERLAP_I6y_D2z = I_TWOBODYOVERLAP_K6yz_Pz+ABZ*I_TWOBODYOVERLAP_I6y_Pz;
  Double I_TWOBODYOVERLAP_I5yz_D2z = I_TWOBODYOVERLAP_K5y2z_Pz+ABZ*I_TWOBODYOVERLAP_I5yz_Pz;
  Double I_TWOBODYOVERLAP_I4y2z_D2z = I_TWOBODYOVERLAP_K4y3z_Pz+ABZ*I_TWOBODYOVERLAP_I4y2z_Pz;
  Double I_TWOBODYOVERLAP_I3y3z_D2z = I_TWOBODYOVERLAP_K3y4z_Pz+ABZ*I_TWOBODYOVERLAP_I3y3z_Pz;
  Double I_TWOBODYOVERLAP_I2y4z_D2z = I_TWOBODYOVERLAP_K2y5z_Pz+ABZ*I_TWOBODYOVERLAP_I2y4z_Pz;
  Double I_TWOBODYOVERLAP_Iy5z_D2z = I_TWOBODYOVERLAP_Ky6z_Pz+ABZ*I_TWOBODYOVERLAP_Iy5z_Pz;
  Double I_TWOBODYOVERLAP_I6z_D2z = I_TWOBODYOVERLAP_K7z_Pz+ABZ*I_TWOBODYOVERLAP_I6z_Pz;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 42 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_F3x = I_TWOBODYOVERLAP_I6x_D2x+ABX*I_TWOBODYOVERLAP_H5x_D2x;
  Double I_TWOBODYOVERLAP_H4xy_F3x = I_TWOBODYOVERLAP_I5xy_D2x+ABX*I_TWOBODYOVERLAP_H4xy_D2x;
  Double I_TWOBODYOVERLAP_H4xz_F3x = I_TWOBODYOVERLAP_I5xz_D2x+ABX*I_TWOBODYOVERLAP_H4xz_D2x;
  Double I_TWOBODYOVERLAP_H3x2y_F3x = I_TWOBODYOVERLAP_I4x2y_D2x+ABX*I_TWOBODYOVERLAP_H3x2y_D2x;
  Double I_TWOBODYOVERLAP_H3xyz_F3x = I_TWOBODYOVERLAP_I4xyz_D2x+ABX*I_TWOBODYOVERLAP_H3xyz_D2x;
  Double I_TWOBODYOVERLAP_H3x2z_F3x = I_TWOBODYOVERLAP_I4x2z_D2x+ABX*I_TWOBODYOVERLAP_H3x2z_D2x;
  Double I_TWOBODYOVERLAP_H2x3y_F3x = I_TWOBODYOVERLAP_I3x3y_D2x+ABX*I_TWOBODYOVERLAP_H2x3y_D2x;
  Double I_TWOBODYOVERLAP_H2x2yz_F3x = I_TWOBODYOVERLAP_I3x2yz_D2x+ABX*I_TWOBODYOVERLAP_H2x2yz_D2x;
  Double I_TWOBODYOVERLAP_H2xy2z_F3x = I_TWOBODYOVERLAP_I3xy2z_D2x+ABX*I_TWOBODYOVERLAP_H2xy2z_D2x;
  Double I_TWOBODYOVERLAP_H2x3z_F3x = I_TWOBODYOVERLAP_I3x3z_D2x+ABX*I_TWOBODYOVERLAP_H2x3z_D2x;
  Double I_TWOBODYOVERLAP_Hx4y_F3x = I_TWOBODYOVERLAP_I2x4y_D2x+ABX*I_TWOBODYOVERLAP_Hx4y_D2x;
  Double I_TWOBODYOVERLAP_Hx3yz_F3x = I_TWOBODYOVERLAP_I2x3yz_D2x+ABX*I_TWOBODYOVERLAP_Hx3yz_D2x;
  Double I_TWOBODYOVERLAP_Hx2y2z_F3x = I_TWOBODYOVERLAP_I2x2y2z_D2x+ABX*I_TWOBODYOVERLAP_Hx2y2z_D2x;
  Double I_TWOBODYOVERLAP_Hxy3z_F3x = I_TWOBODYOVERLAP_I2xy3z_D2x+ABX*I_TWOBODYOVERLAP_Hxy3z_D2x;
  Double I_TWOBODYOVERLAP_Hx4z_F3x = I_TWOBODYOVERLAP_I2x4z_D2x+ABX*I_TWOBODYOVERLAP_Hx4z_D2x;
  Double I_TWOBODYOVERLAP_H5y_F3x = I_TWOBODYOVERLAP_Ix5y_D2x+ABX*I_TWOBODYOVERLAP_H5y_D2x;
  Double I_TWOBODYOVERLAP_H4yz_F3x = I_TWOBODYOVERLAP_Ix4yz_D2x+ABX*I_TWOBODYOVERLAP_H4yz_D2x;
  Double I_TWOBODYOVERLAP_H3y2z_F3x = I_TWOBODYOVERLAP_Ix3y2z_D2x+ABX*I_TWOBODYOVERLAP_H3y2z_D2x;
  Double I_TWOBODYOVERLAP_H2y3z_F3x = I_TWOBODYOVERLAP_Ix2y3z_D2x+ABX*I_TWOBODYOVERLAP_H2y3z_D2x;
  Double I_TWOBODYOVERLAP_Hy4z_F3x = I_TWOBODYOVERLAP_Ixy4z_D2x+ABX*I_TWOBODYOVERLAP_Hy4z_D2x;
  Double I_TWOBODYOVERLAP_H5z_F3x = I_TWOBODYOVERLAP_Ix5z_D2x+ABX*I_TWOBODYOVERLAP_H5z_D2x;
  Double I_TWOBODYOVERLAP_H5x_F2xy = I_TWOBODYOVERLAP_I5xy_D2x+ABY*I_TWOBODYOVERLAP_H5x_D2x;
  Double I_TWOBODYOVERLAP_H4xy_F2xy = I_TWOBODYOVERLAP_I4x2y_D2x+ABY*I_TWOBODYOVERLAP_H4xy_D2x;
  Double I_TWOBODYOVERLAP_H4xz_F2xy = I_TWOBODYOVERLAP_I4xyz_D2x+ABY*I_TWOBODYOVERLAP_H4xz_D2x;
  Double I_TWOBODYOVERLAP_H3x2y_F2xy = I_TWOBODYOVERLAP_I3x3y_D2x+ABY*I_TWOBODYOVERLAP_H3x2y_D2x;
  Double I_TWOBODYOVERLAP_H3xyz_F2xy = I_TWOBODYOVERLAP_I3x2yz_D2x+ABY*I_TWOBODYOVERLAP_H3xyz_D2x;
  Double I_TWOBODYOVERLAP_H3x2z_F2xy = I_TWOBODYOVERLAP_I3xy2z_D2x+ABY*I_TWOBODYOVERLAP_H3x2z_D2x;
  Double I_TWOBODYOVERLAP_H2x3y_F2xy = I_TWOBODYOVERLAP_I2x4y_D2x+ABY*I_TWOBODYOVERLAP_H2x3y_D2x;
  Double I_TWOBODYOVERLAP_H2x2yz_F2xy = I_TWOBODYOVERLAP_I2x3yz_D2x+ABY*I_TWOBODYOVERLAP_H2x2yz_D2x;
  Double I_TWOBODYOVERLAP_H2xy2z_F2xy = I_TWOBODYOVERLAP_I2x2y2z_D2x+ABY*I_TWOBODYOVERLAP_H2xy2z_D2x;
  Double I_TWOBODYOVERLAP_H2x3z_F2xy = I_TWOBODYOVERLAP_I2xy3z_D2x+ABY*I_TWOBODYOVERLAP_H2x3z_D2x;
  Double I_TWOBODYOVERLAP_Hx4y_F2xy = I_TWOBODYOVERLAP_Ix5y_D2x+ABY*I_TWOBODYOVERLAP_Hx4y_D2x;
  Double I_TWOBODYOVERLAP_Hx3yz_F2xy = I_TWOBODYOVERLAP_Ix4yz_D2x+ABY*I_TWOBODYOVERLAP_Hx3yz_D2x;
  Double I_TWOBODYOVERLAP_Hx2y2z_F2xy = I_TWOBODYOVERLAP_Ix3y2z_D2x+ABY*I_TWOBODYOVERLAP_Hx2y2z_D2x;
  Double I_TWOBODYOVERLAP_Hxy3z_F2xy = I_TWOBODYOVERLAP_Ix2y3z_D2x+ABY*I_TWOBODYOVERLAP_Hxy3z_D2x;
  Double I_TWOBODYOVERLAP_Hx4z_F2xy = I_TWOBODYOVERLAP_Ixy4z_D2x+ABY*I_TWOBODYOVERLAP_Hx4z_D2x;
  Double I_TWOBODYOVERLAP_H5y_F2xy = I_TWOBODYOVERLAP_I6y_D2x+ABY*I_TWOBODYOVERLAP_H5y_D2x;
  Double I_TWOBODYOVERLAP_H4yz_F2xy = I_TWOBODYOVERLAP_I5yz_D2x+ABY*I_TWOBODYOVERLAP_H4yz_D2x;
  Double I_TWOBODYOVERLAP_H3y2z_F2xy = I_TWOBODYOVERLAP_I4y2z_D2x+ABY*I_TWOBODYOVERLAP_H3y2z_D2x;
  Double I_TWOBODYOVERLAP_H2y3z_F2xy = I_TWOBODYOVERLAP_I3y3z_D2x+ABY*I_TWOBODYOVERLAP_H2y3z_D2x;
  Double I_TWOBODYOVERLAP_Hy4z_F2xy = I_TWOBODYOVERLAP_I2y4z_D2x+ABY*I_TWOBODYOVERLAP_Hy4z_D2x;
  Double I_TWOBODYOVERLAP_H5z_F2xy = I_TWOBODYOVERLAP_Iy5z_D2x+ABY*I_TWOBODYOVERLAP_H5z_D2x;
  Double I_TWOBODYOVERLAP_H5x_F2xz = I_TWOBODYOVERLAP_I5xz_D2x+ABZ*I_TWOBODYOVERLAP_H5x_D2x;
  Double I_TWOBODYOVERLAP_H4xy_F2xz = I_TWOBODYOVERLAP_I4xyz_D2x+ABZ*I_TWOBODYOVERLAP_H4xy_D2x;
  Double I_TWOBODYOVERLAP_H4xz_F2xz = I_TWOBODYOVERLAP_I4x2z_D2x+ABZ*I_TWOBODYOVERLAP_H4xz_D2x;
  Double I_TWOBODYOVERLAP_H3x2y_F2xz = I_TWOBODYOVERLAP_I3x2yz_D2x+ABZ*I_TWOBODYOVERLAP_H3x2y_D2x;
  Double I_TWOBODYOVERLAP_H3xyz_F2xz = I_TWOBODYOVERLAP_I3xy2z_D2x+ABZ*I_TWOBODYOVERLAP_H3xyz_D2x;
  Double I_TWOBODYOVERLAP_H3x2z_F2xz = I_TWOBODYOVERLAP_I3x3z_D2x+ABZ*I_TWOBODYOVERLAP_H3x2z_D2x;
  Double I_TWOBODYOVERLAP_H2x3y_F2xz = I_TWOBODYOVERLAP_I2x3yz_D2x+ABZ*I_TWOBODYOVERLAP_H2x3y_D2x;
  Double I_TWOBODYOVERLAP_H2x2yz_F2xz = I_TWOBODYOVERLAP_I2x2y2z_D2x+ABZ*I_TWOBODYOVERLAP_H2x2yz_D2x;
  Double I_TWOBODYOVERLAP_H2xy2z_F2xz = I_TWOBODYOVERLAP_I2xy3z_D2x+ABZ*I_TWOBODYOVERLAP_H2xy2z_D2x;
  Double I_TWOBODYOVERLAP_H2x3z_F2xz = I_TWOBODYOVERLAP_I2x4z_D2x+ABZ*I_TWOBODYOVERLAP_H2x3z_D2x;
  Double I_TWOBODYOVERLAP_Hx4y_F2xz = I_TWOBODYOVERLAP_Ix4yz_D2x+ABZ*I_TWOBODYOVERLAP_Hx4y_D2x;
  Double I_TWOBODYOVERLAP_Hx3yz_F2xz = I_TWOBODYOVERLAP_Ix3y2z_D2x+ABZ*I_TWOBODYOVERLAP_Hx3yz_D2x;
  Double I_TWOBODYOVERLAP_Hx2y2z_F2xz = I_TWOBODYOVERLAP_Ix2y3z_D2x+ABZ*I_TWOBODYOVERLAP_Hx2y2z_D2x;
  Double I_TWOBODYOVERLAP_Hxy3z_F2xz = I_TWOBODYOVERLAP_Ixy4z_D2x+ABZ*I_TWOBODYOVERLAP_Hxy3z_D2x;
  Double I_TWOBODYOVERLAP_Hx4z_F2xz = I_TWOBODYOVERLAP_Ix5z_D2x+ABZ*I_TWOBODYOVERLAP_Hx4z_D2x;
  Double I_TWOBODYOVERLAP_H5y_F2xz = I_TWOBODYOVERLAP_I5yz_D2x+ABZ*I_TWOBODYOVERLAP_H5y_D2x;
  Double I_TWOBODYOVERLAP_H4yz_F2xz = I_TWOBODYOVERLAP_I4y2z_D2x+ABZ*I_TWOBODYOVERLAP_H4yz_D2x;
  Double I_TWOBODYOVERLAP_H3y2z_F2xz = I_TWOBODYOVERLAP_I3y3z_D2x+ABZ*I_TWOBODYOVERLAP_H3y2z_D2x;
  Double I_TWOBODYOVERLAP_H2y3z_F2xz = I_TWOBODYOVERLAP_I2y4z_D2x+ABZ*I_TWOBODYOVERLAP_H2y3z_D2x;
  Double I_TWOBODYOVERLAP_Hy4z_F2xz = I_TWOBODYOVERLAP_Iy5z_D2x+ABZ*I_TWOBODYOVERLAP_Hy4z_D2x;
  Double I_TWOBODYOVERLAP_H5z_F2xz = I_TWOBODYOVERLAP_I6z_D2x+ABZ*I_TWOBODYOVERLAP_H5z_D2x;
  Double I_TWOBODYOVERLAP_H5x_Fx2y = I_TWOBODYOVERLAP_I6x_D2y+ABX*I_TWOBODYOVERLAP_H5x_D2y;
  Double I_TWOBODYOVERLAP_H4xy_Fx2y = I_TWOBODYOVERLAP_I5xy_D2y+ABX*I_TWOBODYOVERLAP_H4xy_D2y;
  Double I_TWOBODYOVERLAP_H4xz_Fx2y = I_TWOBODYOVERLAP_I5xz_D2y+ABX*I_TWOBODYOVERLAP_H4xz_D2y;
  Double I_TWOBODYOVERLAP_H3x2y_Fx2y = I_TWOBODYOVERLAP_I4x2y_D2y+ABX*I_TWOBODYOVERLAP_H3x2y_D2y;
  Double I_TWOBODYOVERLAP_H3xyz_Fx2y = I_TWOBODYOVERLAP_I4xyz_D2y+ABX*I_TWOBODYOVERLAP_H3xyz_D2y;
  Double I_TWOBODYOVERLAP_H3x2z_Fx2y = I_TWOBODYOVERLAP_I4x2z_D2y+ABX*I_TWOBODYOVERLAP_H3x2z_D2y;
  Double I_TWOBODYOVERLAP_H2x3y_Fx2y = I_TWOBODYOVERLAP_I3x3y_D2y+ABX*I_TWOBODYOVERLAP_H2x3y_D2y;
  Double I_TWOBODYOVERLAP_H2x2yz_Fx2y = I_TWOBODYOVERLAP_I3x2yz_D2y+ABX*I_TWOBODYOVERLAP_H2x2yz_D2y;
  Double I_TWOBODYOVERLAP_H2xy2z_Fx2y = I_TWOBODYOVERLAP_I3xy2z_D2y+ABX*I_TWOBODYOVERLAP_H2xy2z_D2y;
  Double I_TWOBODYOVERLAP_H2x3z_Fx2y = I_TWOBODYOVERLAP_I3x3z_D2y+ABX*I_TWOBODYOVERLAP_H2x3z_D2y;
  Double I_TWOBODYOVERLAP_Hx4y_Fx2y = I_TWOBODYOVERLAP_I2x4y_D2y+ABX*I_TWOBODYOVERLAP_Hx4y_D2y;
  Double I_TWOBODYOVERLAP_Hx3yz_Fx2y = I_TWOBODYOVERLAP_I2x3yz_D2y+ABX*I_TWOBODYOVERLAP_Hx3yz_D2y;
  Double I_TWOBODYOVERLAP_Hx2y2z_Fx2y = I_TWOBODYOVERLAP_I2x2y2z_D2y+ABX*I_TWOBODYOVERLAP_Hx2y2z_D2y;
  Double I_TWOBODYOVERLAP_Hxy3z_Fx2y = I_TWOBODYOVERLAP_I2xy3z_D2y+ABX*I_TWOBODYOVERLAP_Hxy3z_D2y;
  Double I_TWOBODYOVERLAP_Hx4z_Fx2y = I_TWOBODYOVERLAP_I2x4z_D2y+ABX*I_TWOBODYOVERLAP_Hx4z_D2y;
  Double I_TWOBODYOVERLAP_H5y_Fx2y = I_TWOBODYOVERLAP_Ix5y_D2y+ABX*I_TWOBODYOVERLAP_H5y_D2y;
  Double I_TWOBODYOVERLAP_H4yz_Fx2y = I_TWOBODYOVERLAP_Ix4yz_D2y+ABX*I_TWOBODYOVERLAP_H4yz_D2y;
  Double I_TWOBODYOVERLAP_H3y2z_Fx2y = I_TWOBODYOVERLAP_Ix3y2z_D2y+ABX*I_TWOBODYOVERLAP_H3y2z_D2y;
  Double I_TWOBODYOVERLAP_H2y3z_Fx2y = I_TWOBODYOVERLAP_Ix2y3z_D2y+ABX*I_TWOBODYOVERLAP_H2y3z_D2y;
  Double I_TWOBODYOVERLAP_Hy4z_Fx2y = I_TWOBODYOVERLAP_Ixy4z_D2y+ABX*I_TWOBODYOVERLAP_Hy4z_D2y;
  Double I_TWOBODYOVERLAP_H5z_Fx2y = I_TWOBODYOVERLAP_Ix5z_D2y+ABX*I_TWOBODYOVERLAP_H5z_D2y;
  Double I_TWOBODYOVERLAP_H5x_Fx2z = I_TWOBODYOVERLAP_I6x_D2z+ABX*I_TWOBODYOVERLAP_H5x_D2z;
  Double I_TWOBODYOVERLAP_H4xy_Fx2z = I_TWOBODYOVERLAP_I5xy_D2z+ABX*I_TWOBODYOVERLAP_H4xy_D2z;
  Double I_TWOBODYOVERLAP_H4xz_Fx2z = I_TWOBODYOVERLAP_I5xz_D2z+ABX*I_TWOBODYOVERLAP_H4xz_D2z;
  Double I_TWOBODYOVERLAP_H3x2y_Fx2z = I_TWOBODYOVERLAP_I4x2y_D2z+ABX*I_TWOBODYOVERLAP_H3x2y_D2z;
  Double I_TWOBODYOVERLAP_H3xyz_Fx2z = I_TWOBODYOVERLAP_I4xyz_D2z+ABX*I_TWOBODYOVERLAP_H3xyz_D2z;
  Double I_TWOBODYOVERLAP_H3x2z_Fx2z = I_TWOBODYOVERLAP_I4x2z_D2z+ABX*I_TWOBODYOVERLAP_H3x2z_D2z;
  Double I_TWOBODYOVERLAP_H2x3y_Fx2z = I_TWOBODYOVERLAP_I3x3y_D2z+ABX*I_TWOBODYOVERLAP_H2x3y_D2z;
  Double I_TWOBODYOVERLAP_H2x2yz_Fx2z = I_TWOBODYOVERLAP_I3x2yz_D2z+ABX*I_TWOBODYOVERLAP_H2x2yz_D2z;
  Double I_TWOBODYOVERLAP_H2xy2z_Fx2z = I_TWOBODYOVERLAP_I3xy2z_D2z+ABX*I_TWOBODYOVERLAP_H2xy2z_D2z;
  Double I_TWOBODYOVERLAP_H2x3z_Fx2z = I_TWOBODYOVERLAP_I3x3z_D2z+ABX*I_TWOBODYOVERLAP_H2x3z_D2z;
  Double I_TWOBODYOVERLAP_Hx4y_Fx2z = I_TWOBODYOVERLAP_I2x4y_D2z+ABX*I_TWOBODYOVERLAP_Hx4y_D2z;
  Double I_TWOBODYOVERLAP_Hx3yz_Fx2z = I_TWOBODYOVERLAP_I2x3yz_D2z+ABX*I_TWOBODYOVERLAP_Hx3yz_D2z;
  Double I_TWOBODYOVERLAP_Hx2y2z_Fx2z = I_TWOBODYOVERLAP_I2x2y2z_D2z+ABX*I_TWOBODYOVERLAP_Hx2y2z_D2z;
  Double I_TWOBODYOVERLAP_Hxy3z_Fx2z = I_TWOBODYOVERLAP_I2xy3z_D2z+ABX*I_TWOBODYOVERLAP_Hxy3z_D2z;
  Double I_TWOBODYOVERLAP_Hx4z_Fx2z = I_TWOBODYOVERLAP_I2x4z_D2z+ABX*I_TWOBODYOVERLAP_Hx4z_D2z;
  Double I_TWOBODYOVERLAP_H5y_Fx2z = I_TWOBODYOVERLAP_Ix5y_D2z+ABX*I_TWOBODYOVERLAP_H5y_D2z;
  Double I_TWOBODYOVERLAP_H4yz_Fx2z = I_TWOBODYOVERLAP_Ix4yz_D2z+ABX*I_TWOBODYOVERLAP_H4yz_D2z;
  Double I_TWOBODYOVERLAP_H3y2z_Fx2z = I_TWOBODYOVERLAP_Ix3y2z_D2z+ABX*I_TWOBODYOVERLAP_H3y2z_D2z;
  Double I_TWOBODYOVERLAP_H2y3z_Fx2z = I_TWOBODYOVERLAP_Ix2y3z_D2z+ABX*I_TWOBODYOVERLAP_H2y3z_D2z;
  Double I_TWOBODYOVERLAP_Hy4z_Fx2z = I_TWOBODYOVERLAP_Ixy4z_D2z+ABX*I_TWOBODYOVERLAP_Hy4z_D2z;
  Double I_TWOBODYOVERLAP_H5z_Fx2z = I_TWOBODYOVERLAP_Ix5z_D2z+ABX*I_TWOBODYOVERLAP_H5z_D2z;
  Double I_TWOBODYOVERLAP_H5x_F3y = I_TWOBODYOVERLAP_I5xy_D2y+ABY*I_TWOBODYOVERLAP_H5x_D2y;
  Double I_TWOBODYOVERLAP_H4xy_F3y = I_TWOBODYOVERLAP_I4x2y_D2y+ABY*I_TWOBODYOVERLAP_H4xy_D2y;
  Double I_TWOBODYOVERLAP_H4xz_F3y = I_TWOBODYOVERLAP_I4xyz_D2y+ABY*I_TWOBODYOVERLAP_H4xz_D2y;
  Double I_TWOBODYOVERLAP_H3x2y_F3y = I_TWOBODYOVERLAP_I3x3y_D2y+ABY*I_TWOBODYOVERLAP_H3x2y_D2y;
  Double I_TWOBODYOVERLAP_H3xyz_F3y = I_TWOBODYOVERLAP_I3x2yz_D2y+ABY*I_TWOBODYOVERLAP_H3xyz_D2y;
  Double I_TWOBODYOVERLAP_H3x2z_F3y = I_TWOBODYOVERLAP_I3xy2z_D2y+ABY*I_TWOBODYOVERLAP_H3x2z_D2y;
  Double I_TWOBODYOVERLAP_H2x3y_F3y = I_TWOBODYOVERLAP_I2x4y_D2y+ABY*I_TWOBODYOVERLAP_H2x3y_D2y;
  Double I_TWOBODYOVERLAP_H2x2yz_F3y = I_TWOBODYOVERLAP_I2x3yz_D2y+ABY*I_TWOBODYOVERLAP_H2x2yz_D2y;
  Double I_TWOBODYOVERLAP_H2xy2z_F3y = I_TWOBODYOVERLAP_I2x2y2z_D2y+ABY*I_TWOBODYOVERLAP_H2xy2z_D2y;
  Double I_TWOBODYOVERLAP_H2x3z_F3y = I_TWOBODYOVERLAP_I2xy3z_D2y+ABY*I_TWOBODYOVERLAP_H2x3z_D2y;
  Double I_TWOBODYOVERLAP_Hx4y_F3y = I_TWOBODYOVERLAP_Ix5y_D2y+ABY*I_TWOBODYOVERLAP_Hx4y_D2y;
  Double I_TWOBODYOVERLAP_Hx3yz_F3y = I_TWOBODYOVERLAP_Ix4yz_D2y+ABY*I_TWOBODYOVERLAP_Hx3yz_D2y;
  Double I_TWOBODYOVERLAP_Hx2y2z_F3y = I_TWOBODYOVERLAP_Ix3y2z_D2y+ABY*I_TWOBODYOVERLAP_Hx2y2z_D2y;
  Double I_TWOBODYOVERLAP_Hxy3z_F3y = I_TWOBODYOVERLAP_Ix2y3z_D2y+ABY*I_TWOBODYOVERLAP_Hxy3z_D2y;
  Double I_TWOBODYOVERLAP_Hx4z_F3y = I_TWOBODYOVERLAP_Ixy4z_D2y+ABY*I_TWOBODYOVERLAP_Hx4z_D2y;
  Double I_TWOBODYOVERLAP_H5y_F3y = I_TWOBODYOVERLAP_I6y_D2y+ABY*I_TWOBODYOVERLAP_H5y_D2y;
  Double I_TWOBODYOVERLAP_H4yz_F3y = I_TWOBODYOVERLAP_I5yz_D2y+ABY*I_TWOBODYOVERLAP_H4yz_D2y;
  Double I_TWOBODYOVERLAP_H3y2z_F3y = I_TWOBODYOVERLAP_I4y2z_D2y+ABY*I_TWOBODYOVERLAP_H3y2z_D2y;
  Double I_TWOBODYOVERLAP_H2y3z_F3y = I_TWOBODYOVERLAP_I3y3z_D2y+ABY*I_TWOBODYOVERLAP_H2y3z_D2y;
  Double I_TWOBODYOVERLAP_Hy4z_F3y = I_TWOBODYOVERLAP_I2y4z_D2y+ABY*I_TWOBODYOVERLAP_Hy4z_D2y;
  Double I_TWOBODYOVERLAP_H5z_F3y = I_TWOBODYOVERLAP_Iy5z_D2y+ABY*I_TWOBODYOVERLAP_H5z_D2y;
  Double I_TWOBODYOVERLAP_H5x_F2yz = I_TWOBODYOVERLAP_I5xz_D2y+ABZ*I_TWOBODYOVERLAP_H5x_D2y;
  Double I_TWOBODYOVERLAP_H4xy_F2yz = I_TWOBODYOVERLAP_I4xyz_D2y+ABZ*I_TWOBODYOVERLAP_H4xy_D2y;
  Double I_TWOBODYOVERLAP_H4xz_F2yz = I_TWOBODYOVERLAP_I4x2z_D2y+ABZ*I_TWOBODYOVERLAP_H4xz_D2y;
  Double I_TWOBODYOVERLAP_H3x2y_F2yz = I_TWOBODYOVERLAP_I3x2yz_D2y+ABZ*I_TWOBODYOVERLAP_H3x2y_D2y;
  Double I_TWOBODYOVERLAP_H3xyz_F2yz = I_TWOBODYOVERLAP_I3xy2z_D2y+ABZ*I_TWOBODYOVERLAP_H3xyz_D2y;
  Double I_TWOBODYOVERLAP_H3x2z_F2yz = I_TWOBODYOVERLAP_I3x3z_D2y+ABZ*I_TWOBODYOVERLAP_H3x2z_D2y;
  Double I_TWOBODYOVERLAP_H2x3y_F2yz = I_TWOBODYOVERLAP_I2x3yz_D2y+ABZ*I_TWOBODYOVERLAP_H2x3y_D2y;
  Double I_TWOBODYOVERLAP_H2x2yz_F2yz = I_TWOBODYOVERLAP_I2x2y2z_D2y+ABZ*I_TWOBODYOVERLAP_H2x2yz_D2y;
  Double I_TWOBODYOVERLAP_H2xy2z_F2yz = I_TWOBODYOVERLAP_I2xy3z_D2y+ABZ*I_TWOBODYOVERLAP_H2xy2z_D2y;
  Double I_TWOBODYOVERLAP_H2x3z_F2yz = I_TWOBODYOVERLAP_I2x4z_D2y+ABZ*I_TWOBODYOVERLAP_H2x3z_D2y;
  Double I_TWOBODYOVERLAP_Hx4y_F2yz = I_TWOBODYOVERLAP_Ix4yz_D2y+ABZ*I_TWOBODYOVERLAP_Hx4y_D2y;
  Double I_TWOBODYOVERLAP_Hx3yz_F2yz = I_TWOBODYOVERLAP_Ix3y2z_D2y+ABZ*I_TWOBODYOVERLAP_Hx3yz_D2y;
  Double I_TWOBODYOVERLAP_Hx2y2z_F2yz = I_TWOBODYOVERLAP_Ix2y3z_D2y+ABZ*I_TWOBODYOVERLAP_Hx2y2z_D2y;
  Double I_TWOBODYOVERLAP_Hxy3z_F2yz = I_TWOBODYOVERLAP_Ixy4z_D2y+ABZ*I_TWOBODYOVERLAP_Hxy3z_D2y;
  Double I_TWOBODYOVERLAP_Hx4z_F2yz = I_TWOBODYOVERLAP_Ix5z_D2y+ABZ*I_TWOBODYOVERLAP_Hx4z_D2y;
  Double I_TWOBODYOVERLAP_H5y_F2yz = I_TWOBODYOVERLAP_I5yz_D2y+ABZ*I_TWOBODYOVERLAP_H5y_D2y;
  Double I_TWOBODYOVERLAP_H4yz_F2yz = I_TWOBODYOVERLAP_I4y2z_D2y+ABZ*I_TWOBODYOVERLAP_H4yz_D2y;
  Double I_TWOBODYOVERLAP_H3y2z_F2yz = I_TWOBODYOVERLAP_I3y3z_D2y+ABZ*I_TWOBODYOVERLAP_H3y2z_D2y;
  Double I_TWOBODYOVERLAP_H2y3z_F2yz = I_TWOBODYOVERLAP_I2y4z_D2y+ABZ*I_TWOBODYOVERLAP_H2y3z_D2y;
  Double I_TWOBODYOVERLAP_Hy4z_F2yz = I_TWOBODYOVERLAP_Iy5z_D2y+ABZ*I_TWOBODYOVERLAP_Hy4z_D2y;
  Double I_TWOBODYOVERLAP_H5z_F2yz = I_TWOBODYOVERLAP_I6z_D2y+ABZ*I_TWOBODYOVERLAP_H5z_D2y;
  Double I_TWOBODYOVERLAP_H5x_F3z = I_TWOBODYOVERLAP_I5xz_D2z+ABZ*I_TWOBODYOVERLAP_H5x_D2z;
  Double I_TWOBODYOVERLAP_H4xy_F3z = I_TWOBODYOVERLAP_I4xyz_D2z+ABZ*I_TWOBODYOVERLAP_H4xy_D2z;
  Double I_TWOBODYOVERLAP_H4xz_F3z = I_TWOBODYOVERLAP_I4x2z_D2z+ABZ*I_TWOBODYOVERLAP_H4xz_D2z;
  Double I_TWOBODYOVERLAP_H3x2y_F3z = I_TWOBODYOVERLAP_I3x2yz_D2z+ABZ*I_TWOBODYOVERLAP_H3x2y_D2z;
  Double I_TWOBODYOVERLAP_H3xyz_F3z = I_TWOBODYOVERLAP_I3xy2z_D2z+ABZ*I_TWOBODYOVERLAP_H3xyz_D2z;
  Double I_TWOBODYOVERLAP_H3x2z_F3z = I_TWOBODYOVERLAP_I3x3z_D2z+ABZ*I_TWOBODYOVERLAP_H3x2z_D2z;
  Double I_TWOBODYOVERLAP_H2x3y_F3z = I_TWOBODYOVERLAP_I2x3yz_D2z+ABZ*I_TWOBODYOVERLAP_H2x3y_D2z;
  Double I_TWOBODYOVERLAP_H2x2yz_F3z = I_TWOBODYOVERLAP_I2x2y2z_D2z+ABZ*I_TWOBODYOVERLAP_H2x2yz_D2z;
  Double I_TWOBODYOVERLAP_H2xy2z_F3z = I_TWOBODYOVERLAP_I2xy3z_D2z+ABZ*I_TWOBODYOVERLAP_H2xy2z_D2z;
  Double I_TWOBODYOVERLAP_H2x3z_F3z = I_TWOBODYOVERLAP_I2x4z_D2z+ABZ*I_TWOBODYOVERLAP_H2x3z_D2z;
  Double I_TWOBODYOVERLAP_Hx4y_F3z = I_TWOBODYOVERLAP_Ix4yz_D2z+ABZ*I_TWOBODYOVERLAP_Hx4y_D2z;
  Double I_TWOBODYOVERLAP_Hx3yz_F3z = I_TWOBODYOVERLAP_Ix3y2z_D2z+ABZ*I_TWOBODYOVERLAP_Hx3yz_D2z;
  Double I_TWOBODYOVERLAP_Hx2y2z_F3z = I_TWOBODYOVERLAP_Ix2y3z_D2z+ABZ*I_TWOBODYOVERLAP_Hx2y2z_D2z;
  Double I_TWOBODYOVERLAP_Hxy3z_F3z = I_TWOBODYOVERLAP_Ixy4z_D2z+ABZ*I_TWOBODYOVERLAP_Hxy3z_D2z;
  Double I_TWOBODYOVERLAP_Hx4z_F3z = I_TWOBODYOVERLAP_Ix5z_D2z+ABZ*I_TWOBODYOVERLAP_Hx4z_D2z;
  Double I_TWOBODYOVERLAP_H5y_F3z = I_TWOBODYOVERLAP_I5yz_D2z+ABZ*I_TWOBODYOVERLAP_H5y_D2z;
  Double I_TWOBODYOVERLAP_H4yz_F3z = I_TWOBODYOVERLAP_I4y2z_D2z+ABZ*I_TWOBODYOVERLAP_H4yz_D2z;
  Double I_TWOBODYOVERLAP_H3y2z_F3z = I_TWOBODYOVERLAP_I3y3z_D2z+ABZ*I_TWOBODYOVERLAP_H3y2z_D2z;
  Double I_TWOBODYOVERLAP_H2y3z_F3z = I_TWOBODYOVERLAP_I2y4z_D2z+ABZ*I_TWOBODYOVERLAP_H2y3z_D2z;
  Double I_TWOBODYOVERLAP_Hy4z_F3z = I_TWOBODYOVERLAP_Iy5z_D2z+ABZ*I_TWOBODYOVERLAP_Hy4z_D2z;
  Double I_TWOBODYOVERLAP_H5z_F3z = I_TWOBODYOVERLAP_I6z_D2z+ABZ*I_TWOBODYOVERLAP_H5z_D2z;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_L_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_M_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_L8x_Px = I_TWOBODYOVERLAP_M9x_S+ABX*I_TWOBODYOVERLAP_L8x_S;
  Double I_TWOBODYOVERLAP_L7xy_Px = I_TWOBODYOVERLAP_M8xy_S+ABX*I_TWOBODYOVERLAP_L7xy_S;
  Double I_TWOBODYOVERLAP_L7xz_Px = I_TWOBODYOVERLAP_M8xz_S+ABX*I_TWOBODYOVERLAP_L7xz_S;
  Double I_TWOBODYOVERLAP_L6x2y_Px = I_TWOBODYOVERLAP_M7x2y_S+ABX*I_TWOBODYOVERLAP_L6x2y_S;
  Double I_TWOBODYOVERLAP_L6xyz_Px = I_TWOBODYOVERLAP_M7xyz_S+ABX*I_TWOBODYOVERLAP_L6xyz_S;
  Double I_TWOBODYOVERLAP_L6x2z_Px = I_TWOBODYOVERLAP_M7x2z_S+ABX*I_TWOBODYOVERLAP_L6x2z_S;
  Double I_TWOBODYOVERLAP_L5x3y_Px = I_TWOBODYOVERLAP_M6x3y_S+ABX*I_TWOBODYOVERLAP_L5x3y_S;
  Double I_TWOBODYOVERLAP_L5x2yz_Px = I_TWOBODYOVERLAP_M6x2yz_S+ABX*I_TWOBODYOVERLAP_L5x2yz_S;
  Double I_TWOBODYOVERLAP_L5xy2z_Px = I_TWOBODYOVERLAP_M6xy2z_S+ABX*I_TWOBODYOVERLAP_L5xy2z_S;
  Double I_TWOBODYOVERLAP_L5x3z_Px = I_TWOBODYOVERLAP_M6x3z_S+ABX*I_TWOBODYOVERLAP_L5x3z_S;
  Double I_TWOBODYOVERLAP_L4x4y_Px = I_TWOBODYOVERLAP_M5x4y_S+ABX*I_TWOBODYOVERLAP_L4x4y_S;
  Double I_TWOBODYOVERLAP_L4x3yz_Px = I_TWOBODYOVERLAP_M5x3yz_S+ABX*I_TWOBODYOVERLAP_L4x3yz_S;
  Double I_TWOBODYOVERLAP_L4x2y2z_Px = I_TWOBODYOVERLAP_M5x2y2z_S+ABX*I_TWOBODYOVERLAP_L4x2y2z_S;
  Double I_TWOBODYOVERLAP_L4xy3z_Px = I_TWOBODYOVERLAP_M5xy3z_S+ABX*I_TWOBODYOVERLAP_L4xy3z_S;
  Double I_TWOBODYOVERLAP_L4x4z_Px = I_TWOBODYOVERLAP_M5x4z_S+ABX*I_TWOBODYOVERLAP_L4x4z_S;
  Double I_TWOBODYOVERLAP_L3x5y_Px = I_TWOBODYOVERLAP_M4x5y_S+ABX*I_TWOBODYOVERLAP_L3x5y_S;
  Double I_TWOBODYOVERLAP_L3x4yz_Px = I_TWOBODYOVERLAP_M4x4yz_S+ABX*I_TWOBODYOVERLAP_L3x4yz_S;
  Double I_TWOBODYOVERLAP_L3x3y2z_Px = I_TWOBODYOVERLAP_M4x3y2z_S+ABX*I_TWOBODYOVERLAP_L3x3y2z_S;
  Double I_TWOBODYOVERLAP_L3x2y3z_Px = I_TWOBODYOVERLAP_M4x2y3z_S+ABX*I_TWOBODYOVERLAP_L3x2y3z_S;
  Double I_TWOBODYOVERLAP_L3xy4z_Px = I_TWOBODYOVERLAP_M4xy4z_S+ABX*I_TWOBODYOVERLAP_L3xy4z_S;
  Double I_TWOBODYOVERLAP_L3x5z_Px = I_TWOBODYOVERLAP_M4x5z_S+ABX*I_TWOBODYOVERLAP_L3x5z_S;
  Double I_TWOBODYOVERLAP_L2x6y_Px = I_TWOBODYOVERLAP_M3x6y_S+ABX*I_TWOBODYOVERLAP_L2x6y_S;
  Double I_TWOBODYOVERLAP_L2x5yz_Px = I_TWOBODYOVERLAP_M3x5yz_S+ABX*I_TWOBODYOVERLAP_L2x5yz_S;
  Double I_TWOBODYOVERLAP_L2x4y2z_Px = I_TWOBODYOVERLAP_M3x4y2z_S+ABX*I_TWOBODYOVERLAP_L2x4y2z_S;
  Double I_TWOBODYOVERLAP_L2x3y3z_Px = I_TWOBODYOVERLAP_M3x3y3z_S+ABX*I_TWOBODYOVERLAP_L2x3y3z_S;
  Double I_TWOBODYOVERLAP_L2x2y4z_Px = I_TWOBODYOVERLAP_M3x2y4z_S+ABX*I_TWOBODYOVERLAP_L2x2y4z_S;
  Double I_TWOBODYOVERLAP_L2xy5z_Px = I_TWOBODYOVERLAP_M3xy5z_S+ABX*I_TWOBODYOVERLAP_L2xy5z_S;
  Double I_TWOBODYOVERLAP_L2x6z_Px = I_TWOBODYOVERLAP_M3x6z_S+ABX*I_TWOBODYOVERLAP_L2x6z_S;
  Double I_TWOBODYOVERLAP_Lx7y_Px = I_TWOBODYOVERLAP_M2x7y_S+ABX*I_TWOBODYOVERLAP_Lx7y_S;
  Double I_TWOBODYOVERLAP_Lx6yz_Px = I_TWOBODYOVERLAP_M2x6yz_S+ABX*I_TWOBODYOVERLAP_Lx6yz_S;
  Double I_TWOBODYOVERLAP_Lx5y2z_Px = I_TWOBODYOVERLAP_M2x5y2z_S+ABX*I_TWOBODYOVERLAP_Lx5y2z_S;
  Double I_TWOBODYOVERLAP_Lx4y3z_Px = I_TWOBODYOVERLAP_M2x4y3z_S+ABX*I_TWOBODYOVERLAP_Lx4y3z_S;
  Double I_TWOBODYOVERLAP_Lx3y4z_Px = I_TWOBODYOVERLAP_M2x3y4z_S+ABX*I_TWOBODYOVERLAP_Lx3y4z_S;
  Double I_TWOBODYOVERLAP_Lx2y5z_Px = I_TWOBODYOVERLAP_M2x2y5z_S+ABX*I_TWOBODYOVERLAP_Lx2y5z_S;
  Double I_TWOBODYOVERLAP_Lxy6z_Px = I_TWOBODYOVERLAP_M2xy6z_S+ABX*I_TWOBODYOVERLAP_Lxy6z_S;
  Double I_TWOBODYOVERLAP_Lx7z_Px = I_TWOBODYOVERLAP_M2x7z_S+ABX*I_TWOBODYOVERLAP_Lx7z_S;
  Double I_TWOBODYOVERLAP_L6x2y_Py = I_TWOBODYOVERLAP_M6x3y_S+ABY*I_TWOBODYOVERLAP_L6x2y_S;
  Double I_TWOBODYOVERLAP_L6xyz_Py = I_TWOBODYOVERLAP_M6x2yz_S+ABY*I_TWOBODYOVERLAP_L6xyz_S;
  Double I_TWOBODYOVERLAP_L5x3y_Py = I_TWOBODYOVERLAP_M5x4y_S+ABY*I_TWOBODYOVERLAP_L5x3y_S;
  Double I_TWOBODYOVERLAP_L5x2yz_Py = I_TWOBODYOVERLAP_M5x3yz_S+ABY*I_TWOBODYOVERLAP_L5x2yz_S;
  Double I_TWOBODYOVERLAP_L5xy2z_Py = I_TWOBODYOVERLAP_M5x2y2z_S+ABY*I_TWOBODYOVERLAP_L5xy2z_S;
  Double I_TWOBODYOVERLAP_L4x4y_Py = I_TWOBODYOVERLAP_M4x5y_S+ABY*I_TWOBODYOVERLAP_L4x4y_S;
  Double I_TWOBODYOVERLAP_L4x3yz_Py = I_TWOBODYOVERLAP_M4x4yz_S+ABY*I_TWOBODYOVERLAP_L4x3yz_S;
  Double I_TWOBODYOVERLAP_L4x2y2z_Py = I_TWOBODYOVERLAP_M4x3y2z_S+ABY*I_TWOBODYOVERLAP_L4x2y2z_S;
  Double I_TWOBODYOVERLAP_L4xy3z_Py = I_TWOBODYOVERLAP_M4x2y3z_S+ABY*I_TWOBODYOVERLAP_L4xy3z_S;
  Double I_TWOBODYOVERLAP_L3x5y_Py = I_TWOBODYOVERLAP_M3x6y_S+ABY*I_TWOBODYOVERLAP_L3x5y_S;
  Double I_TWOBODYOVERLAP_L3x4yz_Py = I_TWOBODYOVERLAP_M3x5yz_S+ABY*I_TWOBODYOVERLAP_L3x4yz_S;
  Double I_TWOBODYOVERLAP_L3x3y2z_Py = I_TWOBODYOVERLAP_M3x4y2z_S+ABY*I_TWOBODYOVERLAP_L3x3y2z_S;
  Double I_TWOBODYOVERLAP_L3x2y3z_Py = I_TWOBODYOVERLAP_M3x3y3z_S+ABY*I_TWOBODYOVERLAP_L3x2y3z_S;
  Double I_TWOBODYOVERLAP_L3xy4z_Py = I_TWOBODYOVERLAP_M3x2y4z_S+ABY*I_TWOBODYOVERLAP_L3xy4z_S;
  Double I_TWOBODYOVERLAP_L2x6y_Py = I_TWOBODYOVERLAP_M2x7y_S+ABY*I_TWOBODYOVERLAP_L2x6y_S;
  Double I_TWOBODYOVERLAP_L2x5yz_Py = I_TWOBODYOVERLAP_M2x6yz_S+ABY*I_TWOBODYOVERLAP_L2x5yz_S;
  Double I_TWOBODYOVERLAP_L2x4y2z_Py = I_TWOBODYOVERLAP_M2x5y2z_S+ABY*I_TWOBODYOVERLAP_L2x4y2z_S;
  Double I_TWOBODYOVERLAP_L2x3y3z_Py = I_TWOBODYOVERLAP_M2x4y3z_S+ABY*I_TWOBODYOVERLAP_L2x3y3z_S;
  Double I_TWOBODYOVERLAP_L2x2y4z_Py = I_TWOBODYOVERLAP_M2x3y4z_S+ABY*I_TWOBODYOVERLAP_L2x2y4z_S;
  Double I_TWOBODYOVERLAP_L2xy5z_Py = I_TWOBODYOVERLAP_M2x2y5z_S+ABY*I_TWOBODYOVERLAP_L2xy5z_S;
  Double I_TWOBODYOVERLAP_Lx7y_Py = I_TWOBODYOVERLAP_Mx8y_S+ABY*I_TWOBODYOVERLAP_Lx7y_S;
  Double I_TWOBODYOVERLAP_Lx6yz_Py = I_TWOBODYOVERLAP_Mx7yz_S+ABY*I_TWOBODYOVERLAP_Lx6yz_S;
  Double I_TWOBODYOVERLAP_Lx5y2z_Py = I_TWOBODYOVERLAP_Mx6y2z_S+ABY*I_TWOBODYOVERLAP_Lx5y2z_S;
  Double I_TWOBODYOVERLAP_Lx4y3z_Py = I_TWOBODYOVERLAP_Mx5y3z_S+ABY*I_TWOBODYOVERLAP_Lx4y3z_S;
  Double I_TWOBODYOVERLAP_Lx3y4z_Py = I_TWOBODYOVERLAP_Mx4y4z_S+ABY*I_TWOBODYOVERLAP_Lx3y4z_S;
  Double I_TWOBODYOVERLAP_Lx2y5z_Py = I_TWOBODYOVERLAP_Mx3y5z_S+ABY*I_TWOBODYOVERLAP_Lx2y5z_S;
  Double I_TWOBODYOVERLAP_Lxy6z_Py = I_TWOBODYOVERLAP_Mx2y6z_S+ABY*I_TWOBODYOVERLAP_Lxy6z_S;
  Double I_TWOBODYOVERLAP_L8y_Py = I_TWOBODYOVERLAP_M9y_S+ABY*I_TWOBODYOVERLAP_L8y_S;
  Double I_TWOBODYOVERLAP_L7yz_Py = I_TWOBODYOVERLAP_M8yz_S+ABY*I_TWOBODYOVERLAP_L7yz_S;
  Double I_TWOBODYOVERLAP_L6y2z_Py = I_TWOBODYOVERLAP_M7y2z_S+ABY*I_TWOBODYOVERLAP_L6y2z_S;
  Double I_TWOBODYOVERLAP_L5y3z_Py = I_TWOBODYOVERLAP_M6y3z_S+ABY*I_TWOBODYOVERLAP_L5y3z_S;
  Double I_TWOBODYOVERLAP_L4y4z_Py = I_TWOBODYOVERLAP_M5y4z_S+ABY*I_TWOBODYOVERLAP_L4y4z_S;
  Double I_TWOBODYOVERLAP_L3y5z_Py = I_TWOBODYOVERLAP_M4y5z_S+ABY*I_TWOBODYOVERLAP_L3y5z_S;
  Double I_TWOBODYOVERLAP_L2y6z_Py = I_TWOBODYOVERLAP_M3y6z_S+ABY*I_TWOBODYOVERLAP_L2y6z_S;
  Double I_TWOBODYOVERLAP_Ly7z_Py = I_TWOBODYOVERLAP_M2y7z_S+ABY*I_TWOBODYOVERLAP_Ly7z_S;
  Double I_TWOBODYOVERLAP_L6xyz_Pz = I_TWOBODYOVERLAP_M6xy2z_S+ABZ*I_TWOBODYOVERLAP_L6xyz_S;
  Double I_TWOBODYOVERLAP_L6x2z_Pz = I_TWOBODYOVERLAP_M6x3z_S+ABZ*I_TWOBODYOVERLAP_L6x2z_S;
  Double I_TWOBODYOVERLAP_L5x2yz_Pz = I_TWOBODYOVERLAP_M5x2y2z_S+ABZ*I_TWOBODYOVERLAP_L5x2yz_S;
  Double I_TWOBODYOVERLAP_L5xy2z_Pz = I_TWOBODYOVERLAP_M5xy3z_S+ABZ*I_TWOBODYOVERLAP_L5xy2z_S;
  Double I_TWOBODYOVERLAP_L5x3z_Pz = I_TWOBODYOVERLAP_M5x4z_S+ABZ*I_TWOBODYOVERLAP_L5x3z_S;
  Double I_TWOBODYOVERLAP_L4x3yz_Pz = I_TWOBODYOVERLAP_M4x3y2z_S+ABZ*I_TWOBODYOVERLAP_L4x3yz_S;
  Double I_TWOBODYOVERLAP_L4x2y2z_Pz = I_TWOBODYOVERLAP_M4x2y3z_S+ABZ*I_TWOBODYOVERLAP_L4x2y2z_S;
  Double I_TWOBODYOVERLAP_L4xy3z_Pz = I_TWOBODYOVERLAP_M4xy4z_S+ABZ*I_TWOBODYOVERLAP_L4xy3z_S;
  Double I_TWOBODYOVERLAP_L4x4z_Pz = I_TWOBODYOVERLAP_M4x5z_S+ABZ*I_TWOBODYOVERLAP_L4x4z_S;
  Double I_TWOBODYOVERLAP_L3x4yz_Pz = I_TWOBODYOVERLAP_M3x4y2z_S+ABZ*I_TWOBODYOVERLAP_L3x4yz_S;
  Double I_TWOBODYOVERLAP_L3x3y2z_Pz = I_TWOBODYOVERLAP_M3x3y3z_S+ABZ*I_TWOBODYOVERLAP_L3x3y2z_S;
  Double I_TWOBODYOVERLAP_L3x2y3z_Pz = I_TWOBODYOVERLAP_M3x2y4z_S+ABZ*I_TWOBODYOVERLAP_L3x2y3z_S;
  Double I_TWOBODYOVERLAP_L3xy4z_Pz = I_TWOBODYOVERLAP_M3xy5z_S+ABZ*I_TWOBODYOVERLAP_L3xy4z_S;
  Double I_TWOBODYOVERLAP_L3x5z_Pz = I_TWOBODYOVERLAP_M3x6z_S+ABZ*I_TWOBODYOVERLAP_L3x5z_S;
  Double I_TWOBODYOVERLAP_L2x5yz_Pz = I_TWOBODYOVERLAP_M2x5y2z_S+ABZ*I_TWOBODYOVERLAP_L2x5yz_S;
  Double I_TWOBODYOVERLAP_L2x4y2z_Pz = I_TWOBODYOVERLAP_M2x4y3z_S+ABZ*I_TWOBODYOVERLAP_L2x4y2z_S;
  Double I_TWOBODYOVERLAP_L2x3y3z_Pz = I_TWOBODYOVERLAP_M2x3y4z_S+ABZ*I_TWOBODYOVERLAP_L2x3y3z_S;
  Double I_TWOBODYOVERLAP_L2x2y4z_Pz = I_TWOBODYOVERLAP_M2x2y5z_S+ABZ*I_TWOBODYOVERLAP_L2x2y4z_S;
  Double I_TWOBODYOVERLAP_L2xy5z_Pz = I_TWOBODYOVERLAP_M2xy6z_S+ABZ*I_TWOBODYOVERLAP_L2xy5z_S;
  Double I_TWOBODYOVERLAP_L2x6z_Pz = I_TWOBODYOVERLAP_M2x7z_S+ABZ*I_TWOBODYOVERLAP_L2x6z_S;
  Double I_TWOBODYOVERLAP_Lx6yz_Pz = I_TWOBODYOVERLAP_Mx6y2z_S+ABZ*I_TWOBODYOVERLAP_Lx6yz_S;
  Double I_TWOBODYOVERLAP_Lx5y2z_Pz = I_TWOBODYOVERLAP_Mx5y3z_S+ABZ*I_TWOBODYOVERLAP_Lx5y2z_S;
  Double I_TWOBODYOVERLAP_Lx4y3z_Pz = I_TWOBODYOVERLAP_Mx4y4z_S+ABZ*I_TWOBODYOVERLAP_Lx4y3z_S;
  Double I_TWOBODYOVERLAP_Lx3y4z_Pz = I_TWOBODYOVERLAP_Mx3y5z_S+ABZ*I_TWOBODYOVERLAP_Lx3y4z_S;
  Double I_TWOBODYOVERLAP_Lx2y5z_Pz = I_TWOBODYOVERLAP_Mx2y6z_S+ABZ*I_TWOBODYOVERLAP_Lx2y5z_S;
  Double I_TWOBODYOVERLAP_Lxy6z_Pz = I_TWOBODYOVERLAP_Mxy7z_S+ABZ*I_TWOBODYOVERLAP_Lxy6z_S;
  Double I_TWOBODYOVERLAP_Lx7z_Pz = I_TWOBODYOVERLAP_Mx8z_S+ABZ*I_TWOBODYOVERLAP_Lx7z_S;
  Double I_TWOBODYOVERLAP_L6y2z_Pz = I_TWOBODYOVERLAP_M6y3z_S+ABZ*I_TWOBODYOVERLAP_L6y2z_S;
  Double I_TWOBODYOVERLAP_L5y3z_Pz = I_TWOBODYOVERLAP_M5y4z_S+ABZ*I_TWOBODYOVERLAP_L5y3z_S;
  Double I_TWOBODYOVERLAP_L4y4z_Pz = I_TWOBODYOVERLAP_M4y5z_S+ABZ*I_TWOBODYOVERLAP_L4y4z_S;
  Double I_TWOBODYOVERLAP_L3y5z_Pz = I_TWOBODYOVERLAP_M3y6z_S+ABZ*I_TWOBODYOVERLAP_L3y5z_S;
  Double I_TWOBODYOVERLAP_L2y6z_Pz = I_TWOBODYOVERLAP_M2y7z_S+ABZ*I_TWOBODYOVERLAP_L2y6z_S;
  Double I_TWOBODYOVERLAP_Ly7z_Pz = I_TWOBODYOVERLAP_My8z_S+ABZ*I_TWOBODYOVERLAP_Ly7z_S;
  Double I_TWOBODYOVERLAP_L8z_Pz = I_TWOBODYOVERLAP_M9z_S+ABZ*I_TWOBODYOVERLAP_L8z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_K_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 111 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P
   ************************************************************/
  Double I_TWOBODYOVERLAP_K7x_D2x = I_TWOBODYOVERLAP_L8x_Px+ABX*I_TWOBODYOVERLAP_K7x_Px;
  Double I_TWOBODYOVERLAP_K6xy_D2x = I_TWOBODYOVERLAP_L7xy_Px+ABX*I_TWOBODYOVERLAP_K6xy_Px;
  Double I_TWOBODYOVERLAP_K6xz_D2x = I_TWOBODYOVERLAP_L7xz_Px+ABX*I_TWOBODYOVERLAP_K6xz_Px;
  Double I_TWOBODYOVERLAP_K5x2y_D2x = I_TWOBODYOVERLAP_L6x2y_Px+ABX*I_TWOBODYOVERLAP_K5x2y_Px;
  Double I_TWOBODYOVERLAP_K5xyz_D2x = I_TWOBODYOVERLAP_L6xyz_Px+ABX*I_TWOBODYOVERLAP_K5xyz_Px;
  Double I_TWOBODYOVERLAP_K5x2z_D2x = I_TWOBODYOVERLAP_L6x2z_Px+ABX*I_TWOBODYOVERLAP_K5x2z_Px;
  Double I_TWOBODYOVERLAP_K4x3y_D2x = I_TWOBODYOVERLAP_L5x3y_Px+ABX*I_TWOBODYOVERLAP_K4x3y_Px;
  Double I_TWOBODYOVERLAP_K4x2yz_D2x = I_TWOBODYOVERLAP_L5x2yz_Px+ABX*I_TWOBODYOVERLAP_K4x2yz_Px;
  Double I_TWOBODYOVERLAP_K4xy2z_D2x = I_TWOBODYOVERLAP_L5xy2z_Px+ABX*I_TWOBODYOVERLAP_K4xy2z_Px;
  Double I_TWOBODYOVERLAP_K4x3z_D2x = I_TWOBODYOVERLAP_L5x3z_Px+ABX*I_TWOBODYOVERLAP_K4x3z_Px;
  Double I_TWOBODYOVERLAP_K3x4y_D2x = I_TWOBODYOVERLAP_L4x4y_Px+ABX*I_TWOBODYOVERLAP_K3x4y_Px;
  Double I_TWOBODYOVERLAP_K3x3yz_D2x = I_TWOBODYOVERLAP_L4x3yz_Px+ABX*I_TWOBODYOVERLAP_K3x3yz_Px;
  Double I_TWOBODYOVERLAP_K3x2y2z_D2x = I_TWOBODYOVERLAP_L4x2y2z_Px+ABX*I_TWOBODYOVERLAP_K3x2y2z_Px;
  Double I_TWOBODYOVERLAP_K3xy3z_D2x = I_TWOBODYOVERLAP_L4xy3z_Px+ABX*I_TWOBODYOVERLAP_K3xy3z_Px;
  Double I_TWOBODYOVERLAP_K3x4z_D2x = I_TWOBODYOVERLAP_L4x4z_Px+ABX*I_TWOBODYOVERLAP_K3x4z_Px;
  Double I_TWOBODYOVERLAP_K2x5y_D2x = I_TWOBODYOVERLAP_L3x5y_Px+ABX*I_TWOBODYOVERLAP_K2x5y_Px;
  Double I_TWOBODYOVERLAP_K2x4yz_D2x = I_TWOBODYOVERLAP_L3x4yz_Px+ABX*I_TWOBODYOVERLAP_K2x4yz_Px;
  Double I_TWOBODYOVERLAP_K2x3y2z_D2x = I_TWOBODYOVERLAP_L3x3y2z_Px+ABX*I_TWOBODYOVERLAP_K2x3y2z_Px;
  Double I_TWOBODYOVERLAP_K2x2y3z_D2x = I_TWOBODYOVERLAP_L3x2y3z_Px+ABX*I_TWOBODYOVERLAP_K2x2y3z_Px;
  Double I_TWOBODYOVERLAP_K2xy4z_D2x = I_TWOBODYOVERLAP_L3xy4z_Px+ABX*I_TWOBODYOVERLAP_K2xy4z_Px;
  Double I_TWOBODYOVERLAP_K2x5z_D2x = I_TWOBODYOVERLAP_L3x5z_Px+ABX*I_TWOBODYOVERLAP_K2x5z_Px;
  Double I_TWOBODYOVERLAP_Kx6y_D2x = I_TWOBODYOVERLAP_L2x6y_Px+ABX*I_TWOBODYOVERLAP_Kx6y_Px;
  Double I_TWOBODYOVERLAP_Kx5yz_D2x = I_TWOBODYOVERLAP_L2x5yz_Px+ABX*I_TWOBODYOVERLAP_Kx5yz_Px;
  Double I_TWOBODYOVERLAP_Kx4y2z_D2x = I_TWOBODYOVERLAP_L2x4y2z_Px+ABX*I_TWOBODYOVERLAP_Kx4y2z_Px;
  Double I_TWOBODYOVERLAP_Kx3y3z_D2x = I_TWOBODYOVERLAP_L2x3y3z_Px+ABX*I_TWOBODYOVERLAP_Kx3y3z_Px;
  Double I_TWOBODYOVERLAP_Kx2y4z_D2x = I_TWOBODYOVERLAP_L2x2y4z_Px+ABX*I_TWOBODYOVERLAP_Kx2y4z_Px;
  Double I_TWOBODYOVERLAP_Kxy5z_D2x = I_TWOBODYOVERLAP_L2xy5z_Px+ABX*I_TWOBODYOVERLAP_Kxy5z_Px;
  Double I_TWOBODYOVERLAP_Kx6z_D2x = I_TWOBODYOVERLAP_L2x6z_Px+ABX*I_TWOBODYOVERLAP_Kx6z_Px;
  Double I_TWOBODYOVERLAP_K7y_D2x = I_TWOBODYOVERLAP_Lx7y_Px+ABX*I_TWOBODYOVERLAP_K7y_Px;
  Double I_TWOBODYOVERLAP_K6yz_D2x = I_TWOBODYOVERLAP_Lx6yz_Px+ABX*I_TWOBODYOVERLAP_K6yz_Px;
  Double I_TWOBODYOVERLAP_K5y2z_D2x = I_TWOBODYOVERLAP_Lx5y2z_Px+ABX*I_TWOBODYOVERLAP_K5y2z_Px;
  Double I_TWOBODYOVERLAP_K4y3z_D2x = I_TWOBODYOVERLAP_Lx4y3z_Px+ABX*I_TWOBODYOVERLAP_K4y3z_Px;
  Double I_TWOBODYOVERLAP_K3y4z_D2x = I_TWOBODYOVERLAP_Lx3y4z_Px+ABX*I_TWOBODYOVERLAP_K3y4z_Px;
  Double I_TWOBODYOVERLAP_K2y5z_D2x = I_TWOBODYOVERLAP_Lx2y5z_Px+ABX*I_TWOBODYOVERLAP_K2y5z_Px;
  Double I_TWOBODYOVERLAP_Ky6z_D2x = I_TWOBODYOVERLAP_Lxy6z_Px+ABX*I_TWOBODYOVERLAP_Ky6z_Px;
  Double I_TWOBODYOVERLAP_K7z_D2x = I_TWOBODYOVERLAP_Lx7z_Px+ABX*I_TWOBODYOVERLAP_K7z_Px;
  Double I_TWOBODYOVERLAP_K6xy_D2y = I_TWOBODYOVERLAP_L6x2y_Py+ABY*I_TWOBODYOVERLAP_K6xy_Py;
  Double I_TWOBODYOVERLAP_K6xz_D2y = I_TWOBODYOVERLAP_L6xyz_Py+ABY*I_TWOBODYOVERLAP_K6xz_Py;
  Double I_TWOBODYOVERLAP_K5x2y_D2y = I_TWOBODYOVERLAP_L5x3y_Py+ABY*I_TWOBODYOVERLAP_K5x2y_Py;
  Double I_TWOBODYOVERLAP_K5xyz_D2y = I_TWOBODYOVERLAP_L5x2yz_Py+ABY*I_TWOBODYOVERLAP_K5xyz_Py;
  Double I_TWOBODYOVERLAP_K5x2z_D2y = I_TWOBODYOVERLAP_L5xy2z_Py+ABY*I_TWOBODYOVERLAP_K5x2z_Py;
  Double I_TWOBODYOVERLAP_K4x3y_D2y = I_TWOBODYOVERLAP_L4x4y_Py+ABY*I_TWOBODYOVERLAP_K4x3y_Py;
  Double I_TWOBODYOVERLAP_K4x2yz_D2y = I_TWOBODYOVERLAP_L4x3yz_Py+ABY*I_TWOBODYOVERLAP_K4x2yz_Py;
  Double I_TWOBODYOVERLAP_K4xy2z_D2y = I_TWOBODYOVERLAP_L4x2y2z_Py+ABY*I_TWOBODYOVERLAP_K4xy2z_Py;
  Double I_TWOBODYOVERLAP_K4x3z_D2y = I_TWOBODYOVERLAP_L4xy3z_Py+ABY*I_TWOBODYOVERLAP_K4x3z_Py;
  Double I_TWOBODYOVERLAP_K3x4y_D2y = I_TWOBODYOVERLAP_L3x5y_Py+ABY*I_TWOBODYOVERLAP_K3x4y_Py;
  Double I_TWOBODYOVERLAP_K3x3yz_D2y = I_TWOBODYOVERLAP_L3x4yz_Py+ABY*I_TWOBODYOVERLAP_K3x3yz_Py;
  Double I_TWOBODYOVERLAP_K3x2y2z_D2y = I_TWOBODYOVERLAP_L3x3y2z_Py+ABY*I_TWOBODYOVERLAP_K3x2y2z_Py;
  Double I_TWOBODYOVERLAP_K3xy3z_D2y = I_TWOBODYOVERLAP_L3x2y3z_Py+ABY*I_TWOBODYOVERLAP_K3xy3z_Py;
  Double I_TWOBODYOVERLAP_K3x4z_D2y = I_TWOBODYOVERLAP_L3xy4z_Py+ABY*I_TWOBODYOVERLAP_K3x4z_Py;
  Double I_TWOBODYOVERLAP_K2x5y_D2y = I_TWOBODYOVERLAP_L2x6y_Py+ABY*I_TWOBODYOVERLAP_K2x5y_Py;
  Double I_TWOBODYOVERLAP_K2x4yz_D2y = I_TWOBODYOVERLAP_L2x5yz_Py+ABY*I_TWOBODYOVERLAP_K2x4yz_Py;
  Double I_TWOBODYOVERLAP_K2x3y2z_D2y = I_TWOBODYOVERLAP_L2x4y2z_Py+ABY*I_TWOBODYOVERLAP_K2x3y2z_Py;
  Double I_TWOBODYOVERLAP_K2x2y3z_D2y = I_TWOBODYOVERLAP_L2x3y3z_Py+ABY*I_TWOBODYOVERLAP_K2x2y3z_Py;
  Double I_TWOBODYOVERLAP_K2xy4z_D2y = I_TWOBODYOVERLAP_L2x2y4z_Py+ABY*I_TWOBODYOVERLAP_K2xy4z_Py;
  Double I_TWOBODYOVERLAP_K2x5z_D2y = I_TWOBODYOVERLAP_L2xy5z_Py+ABY*I_TWOBODYOVERLAP_K2x5z_Py;
  Double I_TWOBODYOVERLAP_Kx6y_D2y = I_TWOBODYOVERLAP_Lx7y_Py+ABY*I_TWOBODYOVERLAP_Kx6y_Py;
  Double I_TWOBODYOVERLAP_Kx5yz_D2y = I_TWOBODYOVERLAP_Lx6yz_Py+ABY*I_TWOBODYOVERLAP_Kx5yz_Py;
  Double I_TWOBODYOVERLAP_Kx4y2z_D2y = I_TWOBODYOVERLAP_Lx5y2z_Py+ABY*I_TWOBODYOVERLAP_Kx4y2z_Py;
  Double I_TWOBODYOVERLAP_Kx3y3z_D2y = I_TWOBODYOVERLAP_Lx4y3z_Py+ABY*I_TWOBODYOVERLAP_Kx3y3z_Py;
  Double I_TWOBODYOVERLAP_Kx2y4z_D2y = I_TWOBODYOVERLAP_Lx3y4z_Py+ABY*I_TWOBODYOVERLAP_Kx2y4z_Py;
  Double I_TWOBODYOVERLAP_Kxy5z_D2y = I_TWOBODYOVERLAP_Lx2y5z_Py+ABY*I_TWOBODYOVERLAP_Kxy5z_Py;
  Double I_TWOBODYOVERLAP_Kx6z_D2y = I_TWOBODYOVERLAP_Lxy6z_Py+ABY*I_TWOBODYOVERLAP_Kx6z_Py;
  Double I_TWOBODYOVERLAP_K7y_D2y = I_TWOBODYOVERLAP_L8y_Py+ABY*I_TWOBODYOVERLAP_K7y_Py;
  Double I_TWOBODYOVERLAP_K6yz_D2y = I_TWOBODYOVERLAP_L7yz_Py+ABY*I_TWOBODYOVERLAP_K6yz_Py;
  Double I_TWOBODYOVERLAP_K5y2z_D2y = I_TWOBODYOVERLAP_L6y2z_Py+ABY*I_TWOBODYOVERLAP_K5y2z_Py;
  Double I_TWOBODYOVERLAP_K4y3z_D2y = I_TWOBODYOVERLAP_L5y3z_Py+ABY*I_TWOBODYOVERLAP_K4y3z_Py;
  Double I_TWOBODYOVERLAP_K3y4z_D2y = I_TWOBODYOVERLAP_L4y4z_Py+ABY*I_TWOBODYOVERLAP_K3y4z_Py;
  Double I_TWOBODYOVERLAP_K2y5z_D2y = I_TWOBODYOVERLAP_L3y5z_Py+ABY*I_TWOBODYOVERLAP_K2y5z_Py;
  Double I_TWOBODYOVERLAP_Ky6z_D2y = I_TWOBODYOVERLAP_L2y6z_Py+ABY*I_TWOBODYOVERLAP_Ky6z_Py;
  Double I_TWOBODYOVERLAP_K7z_D2y = I_TWOBODYOVERLAP_Ly7z_Py+ABY*I_TWOBODYOVERLAP_K7z_Py;
  Double I_TWOBODYOVERLAP_K6xy_D2z = I_TWOBODYOVERLAP_L6xyz_Pz+ABZ*I_TWOBODYOVERLAP_K6xy_Pz;
  Double I_TWOBODYOVERLAP_K6xz_D2z = I_TWOBODYOVERLAP_L6x2z_Pz+ABZ*I_TWOBODYOVERLAP_K6xz_Pz;
  Double I_TWOBODYOVERLAP_K5x2y_D2z = I_TWOBODYOVERLAP_L5x2yz_Pz+ABZ*I_TWOBODYOVERLAP_K5x2y_Pz;
  Double I_TWOBODYOVERLAP_K5xyz_D2z = I_TWOBODYOVERLAP_L5xy2z_Pz+ABZ*I_TWOBODYOVERLAP_K5xyz_Pz;
  Double I_TWOBODYOVERLAP_K5x2z_D2z = I_TWOBODYOVERLAP_L5x3z_Pz+ABZ*I_TWOBODYOVERLAP_K5x2z_Pz;
  Double I_TWOBODYOVERLAP_K4x3y_D2z = I_TWOBODYOVERLAP_L4x3yz_Pz+ABZ*I_TWOBODYOVERLAP_K4x3y_Pz;
  Double I_TWOBODYOVERLAP_K4x2yz_D2z = I_TWOBODYOVERLAP_L4x2y2z_Pz+ABZ*I_TWOBODYOVERLAP_K4x2yz_Pz;
  Double I_TWOBODYOVERLAP_K4xy2z_D2z = I_TWOBODYOVERLAP_L4xy3z_Pz+ABZ*I_TWOBODYOVERLAP_K4xy2z_Pz;
  Double I_TWOBODYOVERLAP_K4x3z_D2z = I_TWOBODYOVERLAP_L4x4z_Pz+ABZ*I_TWOBODYOVERLAP_K4x3z_Pz;
  Double I_TWOBODYOVERLAP_K3x4y_D2z = I_TWOBODYOVERLAP_L3x4yz_Pz+ABZ*I_TWOBODYOVERLAP_K3x4y_Pz;
  Double I_TWOBODYOVERLAP_K3x3yz_D2z = I_TWOBODYOVERLAP_L3x3y2z_Pz+ABZ*I_TWOBODYOVERLAP_K3x3yz_Pz;
  Double I_TWOBODYOVERLAP_K3x2y2z_D2z = I_TWOBODYOVERLAP_L3x2y3z_Pz+ABZ*I_TWOBODYOVERLAP_K3x2y2z_Pz;
  Double I_TWOBODYOVERLAP_K3xy3z_D2z = I_TWOBODYOVERLAP_L3xy4z_Pz+ABZ*I_TWOBODYOVERLAP_K3xy3z_Pz;
  Double I_TWOBODYOVERLAP_K3x4z_D2z = I_TWOBODYOVERLAP_L3x5z_Pz+ABZ*I_TWOBODYOVERLAP_K3x4z_Pz;
  Double I_TWOBODYOVERLAP_K2x5y_D2z = I_TWOBODYOVERLAP_L2x5yz_Pz+ABZ*I_TWOBODYOVERLAP_K2x5y_Pz;
  Double I_TWOBODYOVERLAP_K2x4yz_D2z = I_TWOBODYOVERLAP_L2x4y2z_Pz+ABZ*I_TWOBODYOVERLAP_K2x4yz_Pz;
  Double I_TWOBODYOVERLAP_K2x3y2z_D2z = I_TWOBODYOVERLAP_L2x3y3z_Pz+ABZ*I_TWOBODYOVERLAP_K2x3y2z_Pz;
  Double I_TWOBODYOVERLAP_K2x2y3z_D2z = I_TWOBODYOVERLAP_L2x2y4z_Pz+ABZ*I_TWOBODYOVERLAP_K2x2y3z_Pz;
  Double I_TWOBODYOVERLAP_K2xy4z_D2z = I_TWOBODYOVERLAP_L2xy5z_Pz+ABZ*I_TWOBODYOVERLAP_K2xy4z_Pz;
  Double I_TWOBODYOVERLAP_K2x5z_D2z = I_TWOBODYOVERLAP_L2x6z_Pz+ABZ*I_TWOBODYOVERLAP_K2x5z_Pz;
  Double I_TWOBODYOVERLAP_Kx6y_D2z = I_TWOBODYOVERLAP_Lx6yz_Pz+ABZ*I_TWOBODYOVERLAP_Kx6y_Pz;
  Double I_TWOBODYOVERLAP_Kx5yz_D2z = I_TWOBODYOVERLAP_Lx5y2z_Pz+ABZ*I_TWOBODYOVERLAP_Kx5yz_Pz;
  Double I_TWOBODYOVERLAP_Kx4y2z_D2z = I_TWOBODYOVERLAP_Lx4y3z_Pz+ABZ*I_TWOBODYOVERLAP_Kx4y2z_Pz;
  Double I_TWOBODYOVERLAP_Kx3y3z_D2z = I_TWOBODYOVERLAP_Lx3y4z_Pz+ABZ*I_TWOBODYOVERLAP_Kx3y3z_Pz;
  Double I_TWOBODYOVERLAP_Kx2y4z_D2z = I_TWOBODYOVERLAP_Lx2y5z_Pz+ABZ*I_TWOBODYOVERLAP_Kx2y4z_Pz;
  Double I_TWOBODYOVERLAP_Kxy5z_D2z = I_TWOBODYOVERLAP_Lxy6z_Pz+ABZ*I_TWOBODYOVERLAP_Kxy5z_Pz;
  Double I_TWOBODYOVERLAP_Kx6z_D2z = I_TWOBODYOVERLAP_Lx7z_Pz+ABZ*I_TWOBODYOVERLAP_Kx6z_Pz;
  Double I_TWOBODYOVERLAP_K6yz_D2z = I_TWOBODYOVERLAP_L6y2z_Pz+ABZ*I_TWOBODYOVERLAP_K6yz_Pz;
  Double I_TWOBODYOVERLAP_K5y2z_D2z = I_TWOBODYOVERLAP_L5y3z_Pz+ABZ*I_TWOBODYOVERLAP_K5y2z_Pz;
  Double I_TWOBODYOVERLAP_K4y3z_D2z = I_TWOBODYOVERLAP_L4y4z_Pz+ABZ*I_TWOBODYOVERLAP_K4y3z_Pz;
  Double I_TWOBODYOVERLAP_K3y4z_D2z = I_TWOBODYOVERLAP_L3y5z_Pz+ABZ*I_TWOBODYOVERLAP_K3y4z_Pz;
  Double I_TWOBODYOVERLAP_K2y5z_D2z = I_TWOBODYOVERLAP_L2y6z_Pz+ABZ*I_TWOBODYOVERLAP_K2y5z_Pz;
  Double I_TWOBODYOVERLAP_Ky6z_D2z = I_TWOBODYOVERLAP_Ly7z_Pz+ABZ*I_TWOBODYOVERLAP_Ky6z_Pz;
  Double I_TWOBODYOVERLAP_K7z_D2z = I_TWOBODYOVERLAP_L8z_Pz+ABZ*I_TWOBODYOVERLAP_K7z_Pz;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 85 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_D
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_F3x = I_TWOBODYOVERLAP_K7x_D2x+ABX*I_TWOBODYOVERLAP_I6x_D2x;
  Double I_TWOBODYOVERLAP_I5xy_F3x = I_TWOBODYOVERLAP_K6xy_D2x+ABX*I_TWOBODYOVERLAP_I5xy_D2x;
  Double I_TWOBODYOVERLAP_I5xz_F3x = I_TWOBODYOVERLAP_K6xz_D2x+ABX*I_TWOBODYOVERLAP_I5xz_D2x;
  Double I_TWOBODYOVERLAP_I4x2y_F3x = I_TWOBODYOVERLAP_K5x2y_D2x+ABX*I_TWOBODYOVERLAP_I4x2y_D2x;
  Double I_TWOBODYOVERLAP_I4xyz_F3x = I_TWOBODYOVERLAP_K5xyz_D2x+ABX*I_TWOBODYOVERLAP_I4xyz_D2x;
  Double I_TWOBODYOVERLAP_I4x2z_F3x = I_TWOBODYOVERLAP_K5x2z_D2x+ABX*I_TWOBODYOVERLAP_I4x2z_D2x;
  Double I_TWOBODYOVERLAP_I3x3y_F3x = I_TWOBODYOVERLAP_K4x3y_D2x+ABX*I_TWOBODYOVERLAP_I3x3y_D2x;
  Double I_TWOBODYOVERLAP_I3x2yz_F3x = I_TWOBODYOVERLAP_K4x2yz_D2x+ABX*I_TWOBODYOVERLAP_I3x2yz_D2x;
  Double I_TWOBODYOVERLAP_I3xy2z_F3x = I_TWOBODYOVERLAP_K4xy2z_D2x+ABX*I_TWOBODYOVERLAP_I3xy2z_D2x;
  Double I_TWOBODYOVERLAP_I3x3z_F3x = I_TWOBODYOVERLAP_K4x3z_D2x+ABX*I_TWOBODYOVERLAP_I3x3z_D2x;
  Double I_TWOBODYOVERLAP_I2x4y_F3x = I_TWOBODYOVERLAP_K3x4y_D2x+ABX*I_TWOBODYOVERLAP_I2x4y_D2x;
  Double I_TWOBODYOVERLAP_I2x3yz_F3x = I_TWOBODYOVERLAP_K3x3yz_D2x+ABX*I_TWOBODYOVERLAP_I2x3yz_D2x;
  Double I_TWOBODYOVERLAP_I2x2y2z_F3x = I_TWOBODYOVERLAP_K3x2y2z_D2x+ABX*I_TWOBODYOVERLAP_I2x2y2z_D2x;
  Double I_TWOBODYOVERLAP_I2xy3z_F3x = I_TWOBODYOVERLAP_K3xy3z_D2x+ABX*I_TWOBODYOVERLAP_I2xy3z_D2x;
  Double I_TWOBODYOVERLAP_I2x4z_F3x = I_TWOBODYOVERLAP_K3x4z_D2x+ABX*I_TWOBODYOVERLAP_I2x4z_D2x;
  Double I_TWOBODYOVERLAP_Ix5y_F3x = I_TWOBODYOVERLAP_K2x5y_D2x+ABX*I_TWOBODYOVERLAP_Ix5y_D2x;
  Double I_TWOBODYOVERLAP_Ix4yz_F3x = I_TWOBODYOVERLAP_K2x4yz_D2x+ABX*I_TWOBODYOVERLAP_Ix4yz_D2x;
  Double I_TWOBODYOVERLAP_Ix3y2z_F3x = I_TWOBODYOVERLAP_K2x3y2z_D2x+ABX*I_TWOBODYOVERLAP_Ix3y2z_D2x;
  Double I_TWOBODYOVERLAP_Ix2y3z_F3x = I_TWOBODYOVERLAP_K2x2y3z_D2x+ABX*I_TWOBODYOVERLAP_Ix2y3z_D2x;
  Double I_TWOBODYOVERLAP_Ixy4z_F3x = I_TWOBODYOVERLAP_K2xy4z_D2x+ABX*I_TWOBODYOVERLAP_Ixy4z_D2x;
  Double I_TWOBODYOVERLAP_Ix5z_F3x = I_TWOBODYOVERLAP_K2x5z_D2x+ABX*I_TWOBODYOVERLAP_Ix5z_D2x;
  Double I_TWOBODYOVERLAP_I6y_F3x = I_TWOBODYOVERLAP_Kx6y_D2x+ABX*I_TWOBODYOVERLAP_I6y_D2x;
  Double I_TWOBODYOVERLAP_I5yz_F3x = I_TWOBODYOVERLAP_Kx5yz_D2x+ABX*I_TWOBODYOVERLAP_I5yz_D2x;
  Double I_TWOBODYOVERLAP_I4y2z_F3x = I_TWOBODYOVERLAP_Kx4y2z_D2x+ABX*I_TWOBODYOVERLAP_I4y2z_D2x;
  Double I_TWOBODYOVERLAP_I3y3z_F3x = I_TWOBODYOVERLAP_Kx3y3z_D2x+ABX*I_TWOBODYOVERLAP_I3y3z_D2x;
  Double I_TWOBODYOVERLAP_I2y4z_F3x = I_TWOBODYOVERLAP_Kx2y4z_D2x+ABX*I_TWOBODYOVERLAP_I2y4z_D2x;
  Double I_TWOBODYOVERLAP_Iy5z_F3x = I_TWOBODYOVERLAP_Kxy5z_D2x+ABX*I_TWOBODYOVERLAP_Iy5z_D2x;
  Double I_TWOBODYOVERLAP_I6z_F3x = I_TWOBODYOVERLAP_Kx6z_D2x+ABX*I_TWOBODYOVERLAP_I6z_D2x;
  Double I_TWOBODYOVERLAP_I5xy_F2xy = I_TWOBODYOVERLAP_K5x2y_D2x+ABY*I_TWOBODYOVERLAP_I5xy_D2x;
  Double I_TWOBODYOVERLAP_I5xz_F2xy = I_TWOBODYOVERLAP_K5xyz_D2x+ABY*I_TWOBODYOVERLAP_I5xz_D2x;
  Double I_TWOBODYOVERLAP_I4x2y_F2xy = I_TWOBODYOVERLAP_K4x3y_D2x+ABY*I_TWOBODYOVERLAP_I4x2y_D2x;
  Double I_TWOBODYOVERLAP_I4xyz_F2xy = I_TWOBODYOVERLAP_K4x2yz_D2x+ABY*I_TWOBODYOVERLAP_I4xyz_D2x;
  Double I_TWOBODYOVERLAP_I4x2z_F2xy = I_TWOBODYOVERLAP_K4xy2z_D2x+ABY*I_TWOBODYOVERLAP_I4x2z_D2x;
  Double I_TWOBODYOVERLAP_I3x3y_F2xy = I_TWOBODYOVERLAP_K3x4y_D2x+ABY*I_TWOBODYOVERLAP_I3x3y_D2x;
  Double I_TWOBODYOVERLAP_I3x2yz_F2xy = I_TWOBODYOVERLAP_K3x3yz_D2x+ABY*I_TWOBODYOVERLAP_I3x2yz_D2x;
  Double I_TWOBODYOVERLAP_I3xy2z_F2xy = I_TWOBODYOVERLAP_K3x2y2z_D2x+ABY*I_TWOBODYOVERLAP_I3xy2z_D2x;
  Double I_TWOBODYOVERLAP_I3x3z_F2xy = I_TWOBODYOVERLAP_K3xy3z_D2x+ABY*I_TWOBODYOVERLAP_I3x3z_D2x;
  Double I_TWOBODYOVERLAP_I2x4y_F2xy = I_TWOBODYOVERLAP_K2x5y_D2x+ABY*I_TWOBODYOVERLAP_I2x4y_D2x;
  Double I_TWOBODYOVERLAP_I2x3yz_F2xy = I_TWOBODYOVERLAP_K2x4yz_D2x+ABY*I_TWOBODYOVERLAP_I2x3yz_D2x;
  Double I_TWOBODYOVERLAP_I2x2y2z_F2xy = I_TWOBODYOVERLAP_K2x3y2z_D2x+ABY*I_TWOBODYOVERLAP_I2x2y2z_D2x;
  Double I_TWOBODYOVERLAP_I2xy3z_F2xy = I_TWOBODYOVERLAP_K2x2y3z_D2x+ABY*I_TWOBODYOVERLAP_I2xy3z_D2x;
  Double I_TWOBODYOVERLAP_I2x4z_F2xy = I_TWOBODYOVERLAP_K2xy4z_D2x+ABY*I_TWOBODYOVERLAP_I2x4z_D2x;
  Double I_TWOBODYOVERLAP_Ix5y_F2xy = I_TWOBODYOVERLAP_Kx6y_D2x+ABY*I_TWOBODYOVERLAP_Ix5y_D2x;
  Double I_TWOBODYOVERLAP_Ix4yz_F2xy = I_TWOBODYOVERLAP_Kx5yz_D2x+ABY*I_TWOBODYOVERLAP_Ix4yz_D2x;
  Double I_TWOBODYOVERLAP_Ix3y2z_F2xy = I_TWOBODYOVERLAP_Kx4y2z_D2x+ABY*I_TWOBODYOVERLAP_Ix3y2z_D2x;
  Double I_TWOBODYOVERLAP_Ix2y3z_F2xy = I_TWOBODYOVERLAP_Kx3y3z_D2x+ABY*I_TWOBODYOVERLAP_Ix2y3z_D2x;
  Double I_TWOBODYOVERLAP_Ixy4z_F2xy = I_TWOBODYOVERLAP_Kx2y4z_D2x+ABY*I_TWOBODYOVERLAP_Ixy4z_D2x;
  Double I_TWOBODYOVERLAP_Ix5z_F2xy = I_TWOBODYOVERLAP_Kxy5z_D2x+ABY*I_TWOBODYOVERLAP_Ix5z_D2x;
  Double I_TWOBODYOVERLAP_I6y_F2xy = I_TWOBODYOVERLAP_K7y_D2x+ABY*I_TWOBODYOVERLAP_I6y_D2x;
  Double I_TWOBODYOVERLAP_I5yz_F2xy = I_TWOBODYOVERLAP_K6yz_D2x+ABY*I_TWOBODYOVERLAP_I5yz_D2x;
  Double I_TWOBODYOVERLAP_I4y2z_F2xy = I_TWOBODYOVERLAP_K5y2z_D2x+ABY*I_TWOBODYOVERLAP_I4y2z_D2x;
  Double I_TWOBODYOVERLAP_I3y3z_F2xy = I_TWOBODYOVERLAP_K4y3z_D2x+ABY*I_TWOBODYOVERLAP_I3y3z_D2x;
  Double I_TWOBODYOVERLAP_I2y4z_F2xy = I_TWOBODYOVERLAP_K3y4z_D2x+ABY*I_TWOBODYOVERLAP_I2y4z_D2x;
  Double I_TWOBODYOVERLAP_Iy5z_F2xy = I_TWOBODYOVERLAP_K2y5z_D2x+ABY*I_TWOBODYOVERLAP_Iy5z_D2x;
  Double I_TWOBODYOVERLAP_I6z_F2xy = I_TWOBODYOVERLAP_Ky6z_D2x+ABY*I_TWOBODYOVERLAP_I6z_D2x;
  Double I_TWOBODYOVERLAP_I5xz_F2xz = I_TWOBODYOVERLAP_K5x2z_D2x+ABZ*I_TWOBODYOVERLAP_I5xz_D2x;
  Double I_TWOBODYOVERLAP_I4xyz_F2xz = I_TWOBODYOVERLAP_K4xy2z_D2x+ABZ*I_TWOBODYOVERLAP_I4xyz_D2x;
  Double I_TWOBODYOVERLAP_I4x2z_F2xz = I_TWOBODYOVERLAP_K4x3z_D2x+ABZ*I_TWOBODYOVERLAP_I4x2z_D2x;
  Double I_TWOBODYOVERLAP_I3x2yz_F2xz = I_TWOBODYOVERLAP_K3x2y2z_D2x+ABZ*I_TWOBODYOVERLAP_I3x2yz_D2x;
  Double I_TWOBODYOVERLAP_I3xy2z_F2xz = I_TWOBODYOVERLAP_K3xy3z_D2x+ABZ*I_TWOBODYOVERLAP_I3xy2z_D2x;
  Double I_TWOBODYOVERLAP_I3x3z_F2xz = I_TWOBODYOVERLAP_K3x4z_D2x+ABZ*I_TWOBODYOVERLAP_I3x3z_D2x;
  Double I_TWOBODYOVERLAP_I2x3yz_F2xz = I_TWOBODYOVERLAP_K2x3y2z_D2x+ABZ*I_TWOBODYOVERLAP_I2x3yz_D2x;
  Double I_TWOBODYOVERLAP_I2x2y2z_F2xz = I_TWOBODYOVERLAP_K2x2y3z_D2x+ABZ*I_TWOBODYOVERLAP_I2x2y2z_D2x;
  Double I_TWOBODYOVERLAP_I2xy3z_F2xz = I_TWOBODYOVERLAP_K2xy4z_D2x+ABZ*I_TWOBODYOVERLAP_I2xy3z_D2x;
  Double I_TWOBODYOVERLAP_I2x4z_F2xz = I_TWOBODYOVERLAP_K2x5z_D2x+ABZ*I_TWOBODYOVERLAP_I2x4z_D2x;
  Double I_TWOBODYOVERLAP_Ix4yz_F2xz = I_TWOBODYOVERLAP_Kx4y2z_D2x+ABZ*I_TWOBODYOVERLAP_Ix4yz_D2x;
  Double I_TWOBODYOVERLAP_Ix3y2z_F2xz = I_TWOBODYOVERLAP_Kx3y3z_D2x+ABZ*I_TWOBODYOVERLAP_Ix3y2z_D2x;
  Double I_TWOBODYOVERLAP_Ix2y3z_F2xz = I_TWOBODYOVERLAP_Kx2y4z_D2x+ABZ*I_TWOBODYOVERLAP_Ix2y3z_D2x;
  Double I_TWOBODYOVERLAP_Ixy4z_F2xz = I_TWOBODYOVERLAP_Kxy5z_D2x+ABZ*I_TWOBODYOVERLAP_Ixy4z_D2x;
  Double I_TWOBODYOVERLAP_Ix5z_F2xz = I_TWOBODYOVERLAP_Kx6z_D2x+ABZ*I_TWOBODYOVERLAP_Ix5z_D2x;
  Double I_TWOBODYOVERLAP_I5yz_F2xz = I_TWOBODYOVERLAP_K5y2z_D2x+ABZ*I_TWOBODYOVERLAP_I5yz_D2x;
  Double I_TWOBODYOVERLAP_I4y2z_F2xz = I_TWOBODYOVERLAP_K4y3z_D2x+ABZ*I_TWOBODYOVERLAP_I4y2z_D2x;
  Double I_TWOBODYOVERLAP_I3y3z_F2xz = I_TWOBODYOVERLAP_K3y4z_D2x+ABZ*I_TWOBODYOVERLAP_I3y3z_D2x;
  Double I_TWOBODYOVERLAP_I2y4z_F2xz = I_TWOBODYOVERLAP_K2y5z_D2x+ABZ*I_TWOBODYOVERLAP_I2y4z_D2x;
  Double I_TWOBODYOVERLAP_Iy5z_F2xz = I_TWOBODYOVERLAP_Ky6z_D2x+ABZ*I_TWOBODYOVERLAP_Iy5z_D2x;
  Double I_TWOBODYOVERLAP_I6z_F2xz = I_TWOBODYOVERLAP_K7z_D2x+ABZ*I_TWOBODYOVERLAP_I6z_D2x;
  Double I_TWOBODYOVERLAP_I5xz_Fx2y = I_TWOBODYOVERLAP_K6xz_D2y+ABX*I_TWOBODYOVERLAP_I5xz_D2y;
  Double I_TWOBODYOVERLAP_I4xyz_Fx2y = I_TWOBODYOVERLAP_K5xyz_D2y+ABX*I_TWOBODYOVERLAP_I4xyz_D2y;
  Double I_TWOBODYOVERLAP_I4x2z_Fx2y = I_TWOBODYOVERLAP_K5x2z_D2y+ABX*I_TWOBODYOVERLAP_I4x2z_D2y;
  Double I_TWOBODYOVERLAP_I3x2yz_Fx2y = I_TWOBODYOVERLAP_K4x2yz_D2y+ABX*I_TWOBODYOVERLAP_I3x2yz_D2y;
  Double I_TWOBODYOVERLAP_I3xy2z_Fx2y = I_TWOBODYOVERLAP_K4xy2z_D2y+ABX*I_TWOBODYOVERLAP_I3xy2z_D2y;
  Double I_TWOBODYOVERLAP_I3x3z_Fx2y = I_TWOBODYOVERLAP_K4x3z_D2y+ABX*I_TWOBODYOVERLAP_I3x3z_D2y;
  Double I_TWOBODYOVERLAP_I2x3yz_Fx2y = I_TWOBODYOVERLAP_K3x3yz_D2y+ABX*I_TWOBODYOVERLAP_I2x3yz_D2y;
  Double I_TWOBODYOVERLAP_I2x2y2z_Fx2y = I_TWOBODYOVERLAP_K3x2y2z_D2y+ABX*I_TWOBODYOVERLAP_I2x2y2z_D2y;
  Double I_TWOBODYOVERLAP_I2xy3z_Fx2y = I_TWOBODYOVERLAP_K3xy3z_D2y+ABX*I_TWOBODYOVERLAP_I2xy3z_D2y;
  Double I_TWOBODYOVERLAP_I2x4z_Fx2y = I_TWOBODYOVERLAP_K3x4z_D2y+ABX*I_TWOBODYOVERLAP_I2x4z_D2y;
  Double I_TWOBODYOVERLAP_Ix4yz_Fx2y = I_TWOBODYOVERLAP_K2x4yz_D2y+ABX*I_TWOBODYOVERLAP_Ix4yz_D2y;
  Double I_TWOBODYOVERLAP_Ix3y2z_Fx2y = I_TWOBODYOVERLAP_K2x3y2z_D2y+ABX*I_TWOBODYOVERLAP_Ix3y2z_D2y;
  Double I_TWOBODYOVERLAP_Ix2y3z_Fx2y = I_TWOBODYOVERLAP_K2x2y3z_D2y+ABX*I_TWOBODYOVERLAP_Ix2y3z_D2y;
  Double I_TWOBODYOVERLAP_Ixy4z_Fx2y = I_TWOBODYOVERLAP_K2xy4z_D2y+ABX*I_TWOBODYOVERLAP_Ixy4z_D2y;
  Double I_TWOBODYOVERLAP_Ix5z_Fx2y = I_TWOBODYOVERLAP_K2x5z_D2y+ABX*I_TWOBODYOVERLAP_Ix5z_D2y;
  Double I_TWOBODYOVERLAP_I5yz_Fx2y = I_TWOBODYOVERLAP_Kx5yz_D2y+ABX*I_TWOBODYOVERLAP_I5yz_D2y;
  Double I_TWOBODYOVERLAP_I4y2z_Fx2y = I_TWOBODYOVERLAP_Kx4y2z_D2y+ABX*I_TWOBODYOVERLAP_I4y2z_D2y;
  Double I_TWOBODYOVERLAP_I3y3z_Fx2y = I_TWOBODYOVERLAP_Kx3y3z_D2y+ABX*I_TWOBODYOVERLAP_I3y3z_D2y;
  Double I_TWOBODYOVERLAP_I2y4z_Fx2y = I_TWOBODYOVERLAP_Kx2y4z_D2y+ABX*I_TWOBODYOVERLAP_I2y4z_D2y;
  Double I_TWOBODYOVERLAP_Iy5z_Fx2y = I_TWOBODYOVERLAP_Kxy5z_D2y+ABX*I_TWOBODYOVERLAP_Iy5z_D2y;
  Double I_TWOBODYOVERLAP_I6z_Fx2y = I_TWOBODYOVERLAP_Kx6z_D2y+ABX*I_TWOBODYOVERLAP_I6z_D2y;
  Double I_TWOBODYOVERLAP_I5xy_Fx2z = I_TWOBODYOVERLAP_K6xy_D2z+ABX*I_TWOBODYOVERLAP_I5xy_D2z;
  Double I_TWOBODYOVERLAP_I4x2y_Fx2z = I_TWOBODYOVERLAP_K5x2y_D2z+ABX*I_TWOBODYOVERLAP_I4x2y_D2z;
  Double I_TWOBODYOVERLAP_I4xyz_Fx2z = I_TWOBODYOVERLAP_K5xyz_D2z+ABX*I_TWOBODYOVERLAP_I4xyz_D2z;
  Double I_TWOBODYOVERLAP_I3x3y_Fx2z = I_TWOBODYOVERLAP_K4x3y_D2z+ABX*I_TWOBODYOVERLAP_I3x3y_D2z;
  Double I_TWOBODYOVERLAP_I3x2yz_Fx2z = I_TWOBODYOVERLAP_K4x2yz_D2z+ABX*I_TWOBODYOVERLAP_I3x2yz_D2z;
  Double I_TWOBODYOVERLAP_I3xy2z_Fx2z = I_TWOBODYOVERLAP_K4xy2z_D2z+ABX*I_TWOBODYOVERLAP_I3xy2z_D2z;
  Double I_TWOBODYOVERLAP_I2x4y_Fx2z = I_TWOBODYOVERLAP_K3x4y_D2z+ABX*I_TWOBODYOVERLAP_I2x4y_D2z;
  Double I_TWOBODYOVERLAP_I2x3yz_Fx2z = I_TWOBODYOVERLAP_K3x3yz_D2z+ABX*I_TWOBODYOVERLAP_I2x3yz_D2z;
  Double I_TWOBODYOVERLAP_I2x2y2z_Fx2z = I_TWOBODYOVERLAP_K3x2y2z_D2z+ABX*I_TWOBODYOVERLAP_I2x2y2z_D2z;
  Double I_TWOBODYOVERLAP_I2xy3z_Fx2z = I_TWOBODYOVERLAP_K3xy3z_D2z+ABX*I_TWOBODYOVERLAP_I2xy3z_D2z;
  Double I_TWOBODYOVERLAP_Ix5y_Fx2z = I_TWOBODYOVERLAP_K2x5y_D2z+ABX*I_TWOBODYOVERLAP_Ix5y_D2z;
  Double I_TWOBODYOVERLAP_Ix4yz_Fx2z = I_TWOBODYOVERLAP_K2x4yz_D2z+ABX*I_TWOBODYOVERLAP_Ix4yz_D2z;
  Double I_TWOBODYOVERLAP_Ix3y2z_Fx2z = I_TWOBODYOVERLAP_K2x3y2z_D2z+ABX*I_TWOBODYOVERLAP_Ix3y2z_D2z;
  Double I_TWOBODYOVERLAP_Ix2y3z_Fx2z = I_TWOBODYOVERLAP_K2x2y3z_D2z+ABX*I_TWOBODYOVERLAP_Ix2y3z_D2z;
  Double I_TWOBODYOVERLAP_Ixy4z_Fx2z = I_TWOBODYOVERLAP_K2xy4z_D2z+ABX*I_TWOBODYOVERLAP_Ixy4z_D2z;
  Double I_TWOBODYOVERLAP_I6y_Fx2z = I_TWOBODYOVERLAP_Kx6y_D2z+ABX*I_TWOBODYOVERLAP_I6y_D2z;
  Double I_TWOBODYOVERLAP_I5yz_Fx2z = I_TWOBODYOVERLAP_Kx5yz_D2z+ABX*I_TWOBODYOVERLAP_I5yz_D2z;
  Double I_TWOBODYOVERLAP_I4y2z_Fx2z = I_TWOBODYOVERLAP_Kx4y2z_D2z+ABX*I_TWOBODYOVERLAP_I4y2z_D2z;
  Double I_TWOBODYOVERLAP_I3y3z_Fx2z = I_TWOBODYOVERLAP_Kx3y3z_D2z+ABX*I_TWOBODYOVERLAP_I3y3z_D2z;
  Double I_TWOBODYOVERLAP_I2y4z_Fx2z = I_TWOBODYOVERLAP_Kx2y4z_D2z+ABX*I_TWOBODYOVERLAP_I2y4z_D2z;
  Double I_TWOBODYOVERLAP_Iy5z_Fx2z = I_TWOBODYOVERLAP_Kxy5z_D2z+ABX*I_TWOBODYOVERLAP_Iy5z_D2z;
  Double I_TWOBODYOVERLAP_I6x_F3y = I_TWOBODYOVERLAP_K6xy_D2y+ABY*I_TWOBODYOVERLAP_I6x_D2y;
  Double I_TWOBODYOVERLAP_I5xy_F3y = I_TWOBODYOVERLAP_K5x2y_D2y+ABY*I_TWOBODYOVERLAP_I5xy_D2y;
  Double I_TWOBODYOVERLAP_I5xz_F3y = I_TWOBODYOVERLAP_K5xyz_D2y+ABY*I_TWOBODYOVERLAP_I5xz_D2y;
  Double I_TWOBODYOVERLAP_I4x2y_F3y = I_TWOBODYOVERLAP_K4x3y_D2y+ABY*I_TWOBODYOVERLAP_I4x2y_D2y;
  Double I_TWOBODYOVERLAP_I4xyz_F3y = I_TWOBODYOVERLAP_K4x2yz_D2y+ABY*I_TWOBODYOVERLAP_I4xyz_D2y;
  Double I_TWOBODYOVERLAP_I4x2z_F3y = I_TWOBODYOVERLAP_K4xy2z_D2y+ABY*I_TWOBODYOVERLAP_I4x2z_D2y;
  Double I_TWOBODYOVERLAP_I3x3y_F3y = I_TWOBODYOVERLAP_K3x4y_D2y+ABY*I_TWOBODYOVERLAP_I3x3y_D2y;
  Double I_TWOBODYOVERLAP_I3x2yz_F3y = I_TWOBODYOVERLAP_K3x3yz_D2y+ABY*I_TWOBODYOVERLAP_I3x2yz_D2y;
  Double I_TWOBODYOVERLAP_I3xy2z_F3y = I_TWOBODYOVERLAP_K3x2y2z_D2y+ABY*I_TWOBODYOVERLAP_I3xy2z_D2y;
  Double I_TWOBODYOVERLAP_I3x3z_F3y = I_TWOBODYOVERLAP_K3xy3z_D2y+ABY*I_TWOBODYOVERLAP_I3x3z_D2y;
  Double I_TWOBODYOVERLAP_I2x4y_F3y = I_TWOBODYOVERLAP_K2x5y_D2y+ABY*I_TWOBODYOVERLAP_I2x4y_D2y;
  Double I_TWOBODYOVERLAP_I2x3yz_F3y = I_TWOBODYOVERLAP_K2x4yz_D2y+ABY*I_TWOBODYOVERLAP_I2x3yz_D2y;
  Double I_TWOBODYOVERLAP_I2x2y2z_F3y = I_TWOBODYOVERLAP_K2x3y2z_D2y+ABY*I_TWOBODYOVERLAP_I2x2y2z_D2y;
  Double I_TWOBODYOVERLAP_I2xy3z_F3y = I_TWOBODYOVERLAP_K2x2y3z_D2y+ABY*I_TWOBODYOVERLAP_I2xy3z_D2y;
  Double I_TWOBODYOVERLAP_I2x4z_F3y = I_TWOBODYOVERLAP_K2xy4z_D2y+ABY*I_TWOBODYOVERLAP_I2x4z_D2y;
  Double I_TWOBODYOVERLAP_Ix5y_F3y = I_TWOBODYOVERLAP_Kx6y_D2y+ABY*I_TWOBODYOVERLAP_Ix5y_D2y;
  Double I_TWOBODYOVERLAP_Ix4yz_F3y = I_TWOBODYOVERLAP_Kx5yz_D2y+ABY*I_TWOBODYOVERLAP_Ix4yz_D2y;
  Double I_TWOBODYOVERLAP_Ix3y2z_F3y = I_TWOBODYOVERLAP_Kx4y2z_D2y+ABY*I_TWOBODYOVERLAP_Ix3y2z_D2y;
  Double I_TWOBODYOVERLAP_Ix2y3z_F3y = I_TWOBODYOVERLAP_Kx3y3z_D2y+ABY*I_TWOBODYOVERLAP_Ix2y3z_D2y;
  Double I_TWOBODYOVERLAP_Ixy4z_F3y = I_TWOBODYOVERLAP_Kx2y4z_D2y+ABY*I_TWOBODYOVERLAP_Ixy4z_D2y;
  Double I_TWOBODYOVERLAP_Ix5z_F3y = I_TWOBODYOVERLAP_Kxy5z_D2y+ABY*I_TWOBODYOVERLAP_Ix5z_D2y;
  Double I_TWOBODYOVERLAP_I6y_F3y = I_TWOBODYOVERLAP_K7y_D2y+ABY*I_TWOBODYOVERLAP_I6y_D2y;
  Double I_TWOBODYOVERLAP_I5yz_F3y = I_TWOBODYOVERLAP_K6yz_D2y+ABY*I_TWOBODYOVERLAP_I5yz_D2y;
  Double I_TWOBODYOVERLAP_I4y2z_F3y = I_TWOBODYOVERLAP_K5y2z_D2y+ABY*I_TWOBODYOVERLAP_I4y2z_D2y;
  Double I_TWOBODYOVERLAP_I3y3z_F3y = I_TWOBODYOVERLAP_K4y3z_D2y+ABY*I_TWOBODYOVERLAP_I3y3z_D2y;
  Double I_TWOBODYOVERLAP_I2y4z_F3y = I_TWOBODYOVERLAP_K3y4z_D2y+ABY*I_TWOBODYOVERLAP_I2y4z_D2y;
  Double I_TWOBODYOVERLAP_Iy5z_F3y = I_TWOBODYOVERLAP_K2y5z_D2y+ABY*I_TWOBODYOVERLAP_Iy5z_D2y;
  Double I_TWOBODYOVERLAP_I6z_F3y = I_TWOBODYOVERLAP_Ky6z_D2y+ABY*I_TWOBODYOVERLAP_I6z_D2y;
  Double I_TWOBODYOVERLAP_I5xz_F2yz = I_TWOBODYOVERLAP_K5x2z_D2y+ABZ*I_TWOBODYOVERLAP_I5xz_D2y;
  Double I_TWOBODYOVERLAP_I4xyz_F2yz = I_TWOBODYOVERLAP_K4xy2z_D2y+ABZ*I_TWOBODYOVERLAP_I4xyz_D2y;
  Double I_TWOBODYOVERLAP_I4x2z_F2yz = I_TWOBODYOVERLAP_K4x3z_D2y+ABZ*I_TWOBODYOVERLAP_I4x2z_D2y;
  Double I_TWOBODYOVERLAP_I3x2yz_F2yz = I_TWOBODYOVERLAP_K3x2y2z_D2y+ABZ*I_TWOBODYOVERLAP_I3x2yz_D2y;
  Double I_TWOBODYOVERLAP_I3xy2z_F2yz = I_TWOBODYOVERLAP_K3xy3z_D2y+ABZ*I_TWOBODYOVERLAP_I3xy2z_D2y;
  Double I_TWOBODYOVERLAP_I3x3z_F2yz = I_TWOBODYOVERLAP_K3x4z_D2y+ABZ*I_TWOBODYOVERLAP_I3x3z_D2y;
  Double I_TWOBODYOVERLAP_I2x3yz_F2yz = I_TWOBODYOVERLAP_K2x3y2z_D2y+ABZ*I_TWOBODYOVERLAP_I2x3yz_D2y;
  Double I_TWOBODYOVERLAP_I2x2y2z_F2yz = I_TWOBODYOVERLAP_K2x2y3z_D2y+ABZ*I_TWOBODYOVERLAP_I2x2y2z_D2y;
  Double I_TWOBODYOVERLAP_I2xy3z_F2yz = I_TWOBODYOVERLAP_K2xy4z_D2y+ABZ*I_TWOBODYOVERLAP_I2xy3z_D2y;
  Double I_TWOBODYOVERLAP_I2x4z_F2yz = I_TWOBODYOVERLAP_K2x5z_D2y+ABZ*I_TWOBODYOVERLAP_I2x4z_D2y;
  Double I_TWOBODYOVERLAP_Ix4yz_F2yz = I_TWOBODYOVERLAP_Kx4y2z_D2y+ABZ*I_TWOBODYOVERLAP_Ix4yz_D2y;
  Double I_TWOBODYOVERLAP_Ix3y2z_F2yz = I_TWOBODYOVERLAP_Kx3y3z_D2y+ABZ*I_TWOBODYOVERLAP_Ix3y2z_D2y;
  Double I_TWOBODYOVERLAP_Ix2y3z_F2yz = I_TWOBODYOVERLAP_Kx2y4z_D2y+ABZ*I_TWOBODYOVERLAP_Ix2y3z_D2y;
  Double I_TWOBODYOVERLAP_Ixy4z_F2yz = I_TWOBODYOVERLAP_Kxy5z_D2y+ABZ*I_TWOBODYOVERLAP_Ixy4z_D2y;
  Double I_TWOBODYOVERLAP_Ix5z_F2yz = I_TWOBODYOVERLAP_Kx6z_D2y+ABZ*I_TWOBODYOVERLAP_Ix5z_D2y;
  Double I_TWOBODYOVERLAP_I5yz_F2yz = I_TWOBODYOVERLAP_K5y2z_D2y+ABZ*I_TWOBODYOVERLAP_I5yz_D2y;
  Double I_TWOBODYOVERLAP_I4y2z_F2yz = I_TWOBODYOVERLAP_K4y3z_D2y+ABZ*I_TWOBODYOVERLAP_I4y2z_D2y;
  Double I_TWOBODYOVERLAP_I3y3z_F2yz = I_TWOBODYOVERLAP_K3y4z_D2y+ABZ*I_TWOBODYOVERLAP_I3y3z_D2y;
  Double I_TWOBODYOVERLAP_I2y4z_F2yz = I_TWOBODYOVERLAP_K2y5z_D2y+ABZ*I_TWOBODYOVERLAP_I2y4z_D2y;
  Double I_TWOBODYOVERLAP_Iy5z_F2yz = I_TWOBODYOVERLAP_Ky6z_D2y+ABZ*I_TWOBODYOVERLAP_Iy5z_D2y;
  Double I_TWOBODYOVERLAP_I6z_F2yz = I_TWOBODYOVERLAP_K7z_D2y+ABZ*I_TWOBODYOVERLAP_I6z_D2y;
  Double I_TWOBODYOVERLAP_I6x_F3z = I_TWOBODYOVERLAP_K6xz_D2z+ABZ*I_TWOBODYOVERLAP_I6x_D2z;
  Double I_TWOBODYOVERLAP_I5xy_F3z = I_TWOBODYOVERLAP_K5xyz_D2z+ABZ*I_TWOBODYOVERLAP_I5xy_D2z;
  Double I_TWOBODYOVERLAP_I5xz_F3z = I_TWOBODYOVERLAP_K5x2z_D2z+ABZ*I_TWOBODYOVERLAP_I5xz_D2z;
  Double I_TWOBODYOVERLAP_I4x2y_F3z = I_TWOBODYOVERLAP_K4x2yz_D2z+ABZ*I_TWOBODYOVERLAP_I4x2y_D2z;
  Double I_TWOBODYOVERLAP_I4xyz_F3z = I_TWOBODYOVERLAP_K4xy2z_D2z+ABZ*I_TWOBODYOVERLAP_I4xyz_D2z;
  Double I_TWOBODYOVERLAP_I4x2z_F3z = I_TWOBODYOVERLAP_K4x3z_D2z+ABZ*I_TWOBODYOVERLAP_I4x2z_D2z;
  Double I_TWOBODYOVERLAP_I3x3y_F3z = I_TWOBODYOVERLAP_K3x3yz_D2z+ABZ*I_TWOBODYOVERLAP_I3x3y_D2z;
  Double I_TWOBODYOVERLAP_I3x2yz_F3z = I_TWOBODYOVERLAP_K3x2y2z_D2z+ABZ*I_TWOBODYOVERLAP_I3x2yz_D2z;
  Double I_TWOBODYOVERLAP_I3xy2z_F3z = I_TWOBODYOVERLAP_K3xy3z_D2z+ABZ*I_TWOBODYOVERLAP_I3xy2z_D2z;
  Double I_TWOBODYOVERLAP_I3x3z_F3z = I_TWOBODYOVERLAP_K3x4z_D2z+ABZ*I_TWOBODYOVERLAP_I3x3z_D2z;
  Double I_TWOBODYOVERLAP_I2x4y_F3z = I_TWOBODYOVERLAP_K2x4yz_D2z+ABZ*I_TWOBODYOVERLAP_I2x4y_D2z;
  Double I_TWOBODYOVERLAP_I2x3yz_F3z = I_TWOBODYOVERLAP_K2x3y2z_D2z+ABZ*I_TWOBODYOVERLAP_I2x3yz_D2z;
  Double I_TWOBODYOVERLAP_I2x2y2z_F3z = I_TWOBODYOVERLAP_K2x2y3z_D2z+ABZ*I_TWOBODYOVERLAP_I2x2y2z_D2z;
  Double I_TWOBODYOVERLAP_I2xy3z_F3z = I_TWOBODYOVERLAP_K2xy4z_D2z+ABZ*I_TWOBODYOVERLAP_I2xy3z_D2z;
  Double I_TWOBODYOVERLAP_I2x4z_F3z = I_TWOBODYOVERLAP_K2x5z_D2z+ABZ*I_TWOBODYOVERLAP_I2x4z_D2z;
  Double I_TWOBODYOVERLAP_Ix5y_F3z = I_TWOBODYOVERLAP_Kx5yz_D2z+ABZ*I_TWOBODYOVERLAP_Ix5y_D2z;
  Double I_TWOBODYOVERLAP_Ix4yz_F3z = I_TWOBODYOVERLAP_Kx4y2z_D2z+ABZ*I_TWOBODYOVERLAP_Ix4yz_D2z;
  Double I_TWOBODYOVERLAP_Ix3y2z_F3z = I_TWOBODYOVERLAP_Kx3y3z_D2z+ABZ*I_TWOBODYOVERLAP_Ix3y2z_D2z;
  Double I_TWOBODYOVERLAP_Ix2y3z_F3z = I_TWOBODYOVERLAP_Kx2y4z_D2z+ABZ*I_TWOBODYOVERLAP_Ix2y3z_D2z;
  Double I_TWOBODYOVERLAP_Ixy4z_F3z = I_TWOBODYOVERLAP_Kxy5z_D2z+ABZ*I_TWOBODYOVERLAP_Ixy4z_D2z;
  Double I_TWOBODYOVERLAP_Ix5z_F3z = I_TWOBODYOVERLAP_Kx6z_D2z+ABZ*I_TWOBODYOVERLAP_Ix5z_D2z;
  Double I_TWOBODYOVERLAP_I6y_F3z = I_TWOBODYOVERLAP_K6yz_D2z+ABZ*I_TWOBODYOVERLAP_I6y_D2z;
  Double I_TWOBODYOVERLAP_I5yz_F3z = I_TWOBODYOVERLAP_K5y2z_D2z+ABZ*I_TWOBODYOVERLAP_I5yz_D2z;
  Double I_TWOBODYOVERLAP_I4y2z_F3z = I_TWOBODYOVERLAP_K4y3z_D2z+ABZ*I_TWOBODYOVERLAP_I4y2z_D2z;
  Double I_TWOBODYOVERLAP_I3y3z_F3z = I_TWOBODYOVERLAP_K3y4z_D2z+ABZ*I_TWOBODYOVERLAP_I3y3z_D2z;
  Double I_TWOBODYOVERLAP_I2y4z_F3z = I_TWOBODYOVERLAP_K2y5z_D2z+ABZ*I_TWOBODYOVERLAP_I2y4z_D2z;
  Double I_TWOBODYOVERLAP_Iy5z_F3z = I_TWOBODYOVERLAP_Ky6z_D2z+ABZ*I_TWOBODYOVERLAP_Iy5z_D2z;
  Double I_TWOBODYOVERLAP_I6z_F3z = I_TWOBODYOVERLAP_K7z_D2z+ABZ*I_TWOBODYOVERLAP_I6z_D2z;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_G
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_F
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_F
   ************************************************************/
  abcd[0] = I_TWOBODYOVERLAP_I6x_F3x+ABX*I_TWOBODYOVERLAP_H5x_F3x;
  abcd[1] = I_TWOBODYOVERLAP_I5xy_F3x+ABX*I_TWOBODYOVERLAP_H4xy_F3x;
  abcd[2] = I_TWOBODYOVERLAP_I5xz_F3x+ABX*I_TWOBODYOVERLAP_H4xz_F3x;
  abcd[3] = I_TWOBODYOVERLAP_I4x2y_F3x+ABX*I_TWOBODYOVERLAP_H3x2y_F3x;
  abcd[4] = I_TWOBODYOVERLAP_I4xyz_F3x+ABX*I_TWOBODYOVERLAP_H3xyz_F3x;
  abcd[5] = I_TWOBODYOVERLAP_I4x2z_F3x+ABX*I_TWOBODYOVERLAP_H3x2z_F3x;
  abcd[6] = I_TWOBODYOVERLAP_I3x3y_F3x+ABX*I_TWOBODYOVERLAP_H2x3y_F3x;
  abcd[7] = I_TWOBODYOVERLAP_I3x2yz_F3x+ABX*I_TWOBODYOVERLAP_H2x2yz_F3x;
  abcd[8] = I_TWOBODYOVERLAP_I3xy2z_F3x+ABX*I_TWOBODYOVERLAP_H2xy2z_F3x;
  abcd[9] = I_TWOBODYOVERLAP_I3x3z_F3x+ABX*I_TWOBODYOVERLAP_H2x3z_F3x;
  abcd[10] = I_TWOBODYOVERLAP_I2x4y_F3x+ABX*I_TWOBODYOVERLAP_Hx4y_F3x;
  abcd[11] = I_TWOBODYOVERLAP_I2x3yz_F3x+ABX*I_TWOBODYOVERLAP_Hx3yz_F3x;
  abcd[12] = I_TWOBODYOVERLAP_I2x2y2z_F3x+ABX*I_TWOBODYOVERLAP_Hx2y2z_F3x;
  abcd[13] = I_TWOBODYOVERLAP_I2xy3z_F3x+ABX*I_TWOBODYOVERLAP_Hxy3z_F3x;
  abcd[14] = I_TWOBODYOVERLAP_I2x4z_F3x+ABX*I_TWOBODYOVERLAP_Hx4z_F3x;
  abcd[15] = I_TWOBODYOVERLAP_Ix5y_F3x+ABX*I_TWOBODYOVERLAP_H5y_F3x;
  abcd[16] = I_TWOBODYOVERLAP_Ix4yz_F3x+ABX*I_TWOBODYOVERLAP_H4yz_F3x;
  abcd[17] = I_TWOBODYOVERLAP_Ix3y2z_F3x+ABX*I_TWOBODYOVERLAP_H3y2z_F3x;
  abcd[18] = I_TWOBODYOVERLAP_Ix2y3z_F3x+ABX*I_TWOBODYOVERLAP_H2y3z_F3x;
  abcd[19] = I_TWOBODYOVERLAP_Ixy4z_F3x+ABX*I_TWOBODYOVERLAP_Hy4z_F3x;
  abcd[20] = I_TWOBODYOVERLAP_Ix5z_F3x+ABX*I_TWOBODYOVERLAP_H5z_F3x;
  abcd[21] = I_TWOBODYOVERLAP_I5xy_F3x+ABY*I_TWOBODYOVERLAP_H5x_F3x;
  abcd[22] = I_TWOBODYOVERLAP_I4x2y_F3x+ABY*I_TWOBODYOVERLAP_H4xy_F3x;
  abcd[23] = I_TWOBODYOVERLAP_I4xyz_F3x+ABY*I_TWOBODYOVERLAP_H4xz_F3x;
  abcd[24] = I_TWOBODYOVERLAP_I3x3y_F3x+ABY*I_TWOBODYOVERLAP_H3x2y_F3x;
  abcd[25] = I_TWOBODYOVERLAP_I3x2yz_F3x+ABY*I_TWOBODYOVERLAP_H3xyz_F3x;
  abcd[26] = I_TWOBODYOVERLAP_I3xy2z_F3x+ABY*I_TWOBODYOVERLAP_H3x2z_F3x;
  abcd[27] = I_TWOBODYOVERLAP_I2x4y_F3x+ABY*I_TWOBODYOVERLAP_H2x3y_F3x;
  abcd[28] = I_TWOBODYOVERLAP_I2x3yz_F3x+ABY*I_TWOBODYOVERLAP_H2x2yz_F3x;
  abcd[29] = I_TWOBODYOVERLAP_I2x2y2z_F3x+ABY*I_TWOBODYOVERLAP_H2xy2z_F3x;
  abcd[30] = I_TWOBODYOVERLAP_I2xy3z_F3x+ABY*I_TWOBODYOVERLAP_H2x3z_F3x;
  abcd[31] = I_TWOBODYOVERLAP_Ix5y_F3x+ABY*I_TWOBODYOVERLAP_Hx4y_F3x;
  abcd[32] = I_TWOBODYOVERLAP_Ix4yz_F3x+ABY*I_TWOBODYOVERLAP_Hx3yz_F3x;
  abcd[33] = I_TWOBODYOVERLAP_Ix3y2z_F3x+ABY*I_TWOBODYOVERLAP_Hx2y2z_F3x;
  abcd[34] = I_TWOBODYOVERLAP_Ix2y3z_F3x+ABY*I_TWOBODYOVERLAP_Hxy3z_F3x;
  abcd[35] = I_TWOBODYOVERLAP_Ixy4z_F3x+ABY*I_TWOBODYOVERLAP_Hx4z_F3x;
  abcd[36] = I_TWOBODYOVERLAP_I6y_F3x+ABY*I_TWOBODYOVERLAP_H5y_F3x;
  abcd[37] = I_TWOBODYOVERLAP_I5yz_F3x+ABY*I_TWOBODYOVERLAP_H4yz_F3x;
  abcd[38] = I_TWOBODYOVERLAP_I4y2z_F3x+ABY*I_TWOBODYOVERLAP_H3y2z_F3x;
  abcd[39] = I_TWOBODYOVERLAP_I3y3z_F3x+ABY*I_TWOBODYOVERLAP_H2y3z_F3x;
  abcd[40] = I_TWOBODYOVERLAP_I2y4z_F3x+ABY*I_TWOBODYOVERLAP_Hy4z_F3x;
  abcd[41] = I_TWOBODYOVERLAP_Iy5z_F3x+ABY*I_TWOBODYOVERLAP_H5z_F3x;
  abcd[42] = I_TWOBODYOVERLAP_I5xz_F3x+ABZ*I_TWOBODYOVERLAP_H5x_F3x;
  abcd[43] = I_TWOBODYOVERLAP_I4xyz_F3x+ABZ*I_TWOBODYOVERLAP_H4xy_F3x;
  abcd[44] = I_TWOBODYOVERLAP_I4x2z_F3x+ABZ*I_TWOBODYOVERLAP_H4xz_F3x;
  abcd[45] = I_TWOBODYOVERLAP_I3x2yz_F3x+ABZ*I_TWOBODYOVERLAP_H3x2y_F3x;
  abcd[46] = I_TWOBODYOVERLAP_I3xy2z_F3x+ABZ*I_TWOBODYOVERLAP_H3xyz_F3x;
  abcd[47] = I_TWOBODYOVERLAP_I3x3z_F3x+ABZ*I_TWOBODYOVERLAP_H3x2z_F3x;
  abcd[48] = I_TWOBODYOVERLAP_I2x3yz_F3x+ABZ*I_TWOBODYOVERLAP_H2x3y_F3x;
  abcd[49] = I_TWOBODYOVERLAP_I2x2y2z_F3x+ABZ*I_TWOBODYOVERLAP_H2x2yz_F3x;
  abcd[50] = I_TWOBODYOVERLAP_I2xy3z_F3x+ABZ*I_TWOBODYOVERLAP_H2xy2z_F3x;
  abcd[51] = I_TWOBODYOVERLAP_I2x4z_F3x+ABZ*I_TWOBODYOVERLAP_H2x3z_F3x;
  abcd[52] = I_TWOBODYOVERLAP_Ix4yz_F3x+ABZ*I_TWOBODYOVERLAP_Hx4y_F3x;
  abcd[53] = I_TWOBODYOVERLAP_Ix3y2z_F3x+ABZ*I_TWOBODYOVERLAP_Hx3yz_F3x;
  abcd[54] = I_TWOBODYOVERLAP_Ix2y3z_F3x+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F3x;
  abcd[55] = I_TWOBODYOVERLAP_Ixy4z_F3x+ABZ*I_TWOBODYOVERLAP_Hxy3z_F3x;
  abcd[56] = I_TWOBODYOVERLAP_Ix5z_F3x+ABZ*I_TWOBODYOVERLAP_Hx4z_F3x;
  abcd[57] = I_TWOBODYOVERLAP_I5yz_F3x+ABZ*I_TWOBODYOVERLAP_H5y_F3x;
  abcd[58] = I_TWOBODYOVERLAP_I4y2z_F3x+ABZ*I_TWOBODYOVERLAP_H4yz_F3x;
  abcd[59] = I_TWOBODYOVERLAP_I3y3z_F3x+ABZ*I_TWOBODYOVERLAP_H3y2z_F3x;
  abcd[60] = I_TWOBODYOVERLAP_I2y4z_F3x+ABZ*I_TWOBODYOVERLAP_H2y3z_F3x;
  abcd[61] = I_TWOBODYOVERLAP_Iy5z_F3x+ABZ*I_TWOBODYOVERLAP_Hy4z_F3x;
  abcd[62] = I_TWOBODYOVERLAP_I6z_F3x+ABZ*I_TWOBODYOVERLAP_H5z_F3x;
  abcd[63] = I_TWOBODYOVERLAP_I5xy_F2xy+ABY*I_TWOBODYOVERLAP_H5x_F2xy;
  abcd[64] = I_TWOBODYOVERLAP_I4x2y_F2xy+ABY*I_TWOBODYOVERLAP_H4xy_F2xy;
  abcd[65] = I_TWOBODYOVERLAP_I4xyz_F2xy+ABY*I_TWOBODYOVERLAP_H4xz_F2xy;
  abcd[66] = I_TWOBODYOVERLAP_I3x3y_F2xy+ABY*I_TWOBODYOVERLAP_H3x2y_F2xy;
  abcd[67] = I_TWOBODYOVERLAP_I3x2yz_F2xy+ABY*I_TWOBODYOVERLAP_H3xyz_F2xy;
  abcd[68] = I_TWOBODYOVERLAP_I3xy2z_F2xy+ABY*I_TWOBODYOVERLAP_H3x2z_F2xy;
  abcd[69] = I_TWOBODYOVERLAP_I2x4y_F2xy+ABY*I_TWOBODYOVERLAP_H2x3y_F2xy;
  abcd[70] = I_TWOBODYOVERLAP_I2x3yz_F2xy+ABY*I_TWOBODYOVERLAP_H2x2yz_F2xy;
  abcd[71] = I_TWOBODYOVERLAP_I2x2y2z_F2xy+ABY*I_TWOBODYOVERLAP_H2xy2z_F2xy;
  abcd[72] = I_TWOBODYOVERLAP_I2xy3z_F2xy+ABY*I_TWOBODYOVERLAP_H2x3z_F2xy;
  abcd[73] = I_TWOBODYOVERLAP_Ix5y_F2xy+ABY*I_TWOBODYOVERLAP_Hx4y_F2xy;
  abcd[74] = I_TWOBODYOVERLAP_Ix4yz_F2xy+ABY*I_TWOBODYOVERLAP_Hx3yz_F2xy;
  abcd[75] = I_TWOBODYOVERLAP_Ix3y2z_F2xy+ABY*I_TWOBODYOVERLAP_Hx2y2z_F2xy;
  abcd[76] = I_TWOBODYOVERLAP_Ix2y3z_F2xy+ABY*I_TWOBODYOVERLAP_Hxy3z_F2xy;
  abcd[77] = I_TWOBODYOVERLAP_Ixy4z_F2xy+ABY*I_TWOBODYOVERLAP_Hx4z_F2xy;
  abcd[78] = I_TWOBODYOVERLAP_I6y_F2xy+ABY*I_TWOBODYOVERLAP_H5y_F2xy;
  abcd[79] = I_TWOBODYOVERLAP_I5yz_F2xy+ABY*I_TWOBODYOVERLAP_H4yz_F2xy;
  abcd[80] = I_TWOBODYOVERLAP_I4y2z_F2xy+ABY*I_TWOBODYOVERLAP_H3y2z_F2xy;
  abcd[81] = I_TWOBODYOVERLAP_I3y3z_F2xy+ABY*I_TWOBODYOVERLAP_H2y3z_F2xy;
  abcd[82] = I_TWOBODYOVERLAP_I2y4z_F2xy+ABY*I_TWOBODYOVERLAP_Hy4z_F2xy;
  abcd[83] = I_TWOBODYOVERLAP_Iy5z_F2xy+ABY*I_TWOBODYOVERLAP_H5z_F2xy;
  abcd[84] = I_TWOBODYOVERLAP_I5xz_F2xy+ABZ*I_TWOBODYOVERLAP_H5x_F2xy;
  abcd[85] = I_TWOBODYOVERLAP_I4xyz_F2xy+ABZ*I_TWOBODYOVERLAP_H4xy_F2xy;
  abcd[86] = I_TWOBODYOVERLAP_I4x2z_F2xy+ABZ*I_TWOBODYOVERLAP_H4xz_F2xy;
  abcd[87] = I_TWOBODYOVERLAP_I3x2yz_F2xy+ABZ*I_TWOBODYOVERLAP_H3x2y_F2xy;
  abcd[88] = I_TWOBODYOVERLAP_I3xy2z_F2xy+ABZ*I_TWOBODYOVERLAP_H3xyz_F2xy;
  abcd[89] = I_TWOBODYOVERLAP_I3x3z_F2xy+ABZ*I_TWOBODYOVERLAP_H3x2z_F2xy;
  abcd[90] = I_TWOBODYOVERLAP_I2x3yz_F2xy+ABZ*I_TWOBODYOVERLAP_H2x3y_F2xy;
  abcd[91] = I_TWOBODYOVERLAP_I2x2y2z_F2xy+ABZ*I_TWOBODYOVERLAP_H2x2yz_F2xy;
  abcd[92] = I_TWOBODYOVERLAP_I2xy3z_F2xy+ABZ*I_TWOBODYOVERLAP_H2xy2z_F2xy;
  abcd[93] = I_TWOBODYOVERLAP_I2x4z_F2xy+ABZ*I_TWOBODYOVERLAP_H2x3z_F2xy;
  abcd[94] = I_TWOBODYOVERLAP_Ix4yz_F2xy+ABZ*I_TWOBODYOVERLAP_Hx4y_F2xy;
  abcd[95] = I_TWOBODYOVERLAP_Ix3y2z_F2xy+ABZ*I_TWOBODYOVERLAP_Hx3yz_F2xy;
  abcd[96] = I_TWOBODYOVERLAP_Ix2y3z_F2xy+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F2xy;
  abcd[97] = I_TWOBODYOVERLAP_Ixy4z_F2xy+ABZ*I_TWOBODYOVERLAP_Hxy3z_F2xy;
  abcd[98] = I_TWOBODYOVERLAP_Ix5z_F2xy+ABZ*I_TWOBODYOVERLAP_Hx4z_F2xy;
  abcd[99] = I_TWOBODYOVERLAP_I5yz_F2xy+ABZ*I_TWOBODYOVERLAP_H5y_F2xy;
  abcd[100] = I_TWOBODYOVERLAP_I4y2z_F2xy+ABZ*I_TWOBODYOVERLAP_H4yz_F2xy;
  abcd[101] = I_TWOBODYOVERLAP_I3y3z_F2xy+ABZ*I_TWOBODYOVERLAP_H3y2z_F2xy;
  abcd[102] = I_TWOBODYOVERLAP_I2y4z_F2xy+ABZ*I_TWOBODYOVERLAP_H2y3z_F2xy;
  abcd[103] = I_TWOBODYOVERLAP_Iy5z_F2xy+ABZ*I_TWOBODYOVERLAP_Hy4z_F2xy;
  abcd[104] = I_TWOBODYOVERLAP_I6z_F2xy+ABZ*I_TWOBODYOVERLAP_H5z_F2xy;
  abcd[105] = I_TWOBODYOVERLAP_I5xz_F2xz+ABZ*I_TWOBODYOVERLAP_H5x_F2xz;
  abcd[106] = I_TWOBODYOVERLAP_I4xyz_F2xz+ABZ*I_TWOBODYOVERLAP_H4xy_F2xz;
  abcd[107] = I_TWOBODYOVERLAP_I4x2z_F2xz+ABZ*I_TWOBODYOVERLAP_H4xz_F2xz;
  abcd[108] = I_TWOBODYOVERLAP_I3x2yz_F2xz+ABZ*I_TWOBODYOVERLAP_H3x2y_F2xz;
  abcd[109] = I_TWOBODYOVERLAP_I3xy2z_F2xz+ABZ*I_TWOBODYOVERLAP_H3xyz_F2xz;
  abcd[110] = I_TWOBODYOVERLAP_I3x3z_F2xz+ABZ*I_TWOBODYOVERLAP_H3x2z_F2xz;
  abcd[111] = I_TWOBODYOVERLAP_I2x3yz_F2xz+ABZ*I_TWOBODYOVERLAP_H2x3y_F2xz;
  abcd[112] = I_TWOBODYOVERLAP_I2x2y2z_F2xz+ABZ*I_TWOBODYOVERLAP_H2x2yz_F2xz;
  abcd[113] = I_TWOBODYOVERLAP_I2xy3z_F2xz+ABZ*I_TWOBODYOVERLAP_H2xy2z_F2xz;
  abcd[114] = I_TWOBODYOVERLAP_I2x4z_F2xz+ABZ*I_TWOBODYOVERLAP_H2x3z_F2xz;
  abcd[115] = I_TWOBODYOVERLAP_Ix4yz_F2xz+ABZ*I_TWOBODYOVERLAP_Hx4y_F2xz;
  abcd[116] = I_TWOBODYOVERLAP_Ix3y2z_F2xz+ABZ*I_TWOBODYOVERLAP_Hx3yz_F2xz;
  abcd[117] = I_TWOBODYOVERLAP_Ix2y3z_F2xz+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F2xz;
  abcd[118] = I_TWOBODYOVERLAP_Ixy4z_F2xz+ABZ*I_TWOBODYOVERLAP_Hxy3z_F2xz;
  abcd[119] = I_TWOBODYOVERLAP_Ix5z_F2xz+ABZ*I_TWOBODYOVERLAP_Hx4z_F2xz;
  abcd[120] = I_TWOBODYOVERLAP_I5yz_F2xz+ABZ*I_TWOBODYOVERLAP_H5y_F2xz;
  abcd[121] = I_TWOBODYOVERLAP_I4y2z_F2xz+ABZ*I_TWOBODYOVERLAP_H4yz_F2xz;
  abcd[122] = I_TWOBODYOVERLAP_I3y3z_F2xz+ABZ*I_TWOBODYOVERLAP_H3y2z_F2xz;
  abcd[123] = I_TWOBODYOVERLAP_I2y4z_F2xz+ABZ*I_TWOBODYOVERLAP_H2y3z_F2xz;
  abcd[124] = I_TWOBODYOVERLAP_Iy5z_F2xz+ABZ*I_TWOBODYOVERLAP_Hy4z_F2xz;
  abcd[125] = I_TWOBODYOVERLAP_I6z_F2xz+ABZ*I_TWOBODYOVERLAP_H5z_F2xz;
  abcd[126] = I_TWOBODYOVERLAP_I6x_F3y+ABX*I_TWOBODYOVERLAP_H5x_F3y;
  abcd[127] = I_TWOBODYOVERLAP_I5xy_F3y+ABX*I_TWOBODYOVERLAP_H4xy_F3y;
  abcd[128] = I_TWOBODYOVERLAP_I5xz_F3y+ABX*I_TWOBODYOVERLAP_H4xz_F3y;
  abcd[129] = I_TWOBODYOVERLAP_I4x2y_F3y+ABX*I_TWOBODYOVERLAP_H3x2y_F3y;
  abcd[130] = I_TWOBODYOVERLAP_I4xyz_F3y+ABX*I_TWOBODYOVERLAP_H3xyz_F3y;
  abcd[131] = I_TWOBODYOVERLAP_I4x2z_F3y+ABX*I_TWOBODYOVERLAP_H3x2z_F3y;
  abcd[132] = I_TWOBODYOVERLAP_I3x3y_F3y+ABX*I_TWOBODYOVERLAP_H2x3y_F3y;
  abcd[133] = I_TWOBODYOVERLAP_I3x2yz_F3y+ABX*I_TWOBODYOVERLAP_H2x2yz_F3y;
  abcd[134] = I_TWOBODYOVERLAP_I3xy2z_F3y+ABX*I_TWOBODYOVERLAP_H2xy2z_F3y;
  abcd[135] = I_TWOBODYOVERLAP_I3x3z_F3y+ABX*I_TWOBODYOVERLAP_H2x3z_F3y;
  abcd[136] = I_TWOBODYOVERLAP_I2x4y_F3y+ABX*I_TWOBODYOVERLAP_Hx4y_F3y;
  abcd[137] = I_TWOBODYOVERLAP_I2x3yz_F3y+ABX*I_TWOBODYOVERLAP_Hx3yz_F3y;
  abcd[138] = I_TWOBODYOVERLAP_I2x2y2z_F3y+ABX*I_TWOBODYOVERLAP_Hx2y2z_F3y;
  abcd[139] = I_TWOBODYOVERLAP_I2xy3z_F3y+ABX*I_TWOBODYOVERLAP_Hxy3z_F3y;
  abcd[140] = I_TWOBODYOVERLAP_I2x4z_F3y+ABX*I_TWOBODYOVERLAP_Hx4z_F3y;
  abcd[141] = I_TWOBODYOVERLAP_Ix5y_F3y+ABX*I_TWOBODYOVERLAP_H5y_F3y;
  abcd[142] = I_TWOBODYOVERLAP_Ix4yz_F3y+ABX*I_TWOBODYOVERLAP_H4yz_F3y;
  abcd[143] = I_TWOBODYOVERLAP_Ix3y2z_F3y+ABX*I_TWOBODYOVERLAP_H3y2z_F3y;
  abcd[144] = I_TWOBODYOVERLAP_Ix2y3z_F3y+ABX*I_TWOBODYOVERLAP_H2y3z_F3y;
  abcd[145] = I_TWOBODYOVERLAP_Ixy4z_F3y+ABX*I_TWOBODYOVERLAP_Hy4z_F3y;
  abcd[146] = I_TWOBODYOVERLAP_Ix5z_F3y+ABX*I_TWOBODYOVERLAP_H5z_F3y;
  abcd[147] = I_TWOBODYOVERLAP_I5xz_Fx2y+ABZ*I_TWOBODYOVERLAP_H5x_Fx2y;
  abcd[148] = I_TWOBODYOVERLAP_I4xyz_Fx2y+ABZ*I_TWOBODYOVERLAP_H4xy_Fx2y;
  abcd[149] = I_TWOBODYOVERLAP_I4x2z_Fx2y+ABZ*I_TWOBODYOVERLAP_H4xz_Fx2y;
  abcd[150] = I_TWOBODYOVERLAP_I3x2yz_Fx2y+ABZ*I_TWOBODYOVERLAP_H3x2y_Fx2y;
  abcd[151] = I_TWOBODYOVERLAP_I3xy2z_Fx2y+ABZ*I_TWOBODYOVERLAP_H3xyz_Fx2y;
  abcd[152] = I_TWOBODYOVERLAP_I3x3z_Fx2y+ABZ*I_TWOBODYOVERLAP_H3x2z_Fx2y;
  abcd[153] = I_TWOBODYOVERLAP_I2x3yz_Fx2y+ABZ*I_TWOBODYOVERLAP_H2x3y_Fx2y;
  abcd[154] = I_TWOBODYOVERLAP_I2x2y2z_Fx2y+ABZ*I_TWOBODYOVERLAP_H2x2yz_Fx2y;
  abcd[155] = I_TWOBODYOVERLAP_I2xy3z_Fx2y+ABZ*I_TWOBODYOVERLAP_H2xy2z_Fx2y;
  abcd[156] = I_TWOBODYOVERLAP_I2x4z_Fx2y+ABZ*I_TWOBODYOVERLAP_H2x3z_Fx2y;
  abcd[157] = I_TWOBODYOVERLAP_Ix4yz_Fx2y+ABZ*I_TWOBODYOVERLAP_Hx4y_Fx2y;
  abcd[158] = I_TWOBODYOVERLAP_Ix3y2z_Fx2y+ABZ*I_TWOBODYOVERLAP_Hx3yz_Fx2y;
  abcd[159] = I_TWOBODYOVERLAP_Ix2y3z_Fx2y+ABZ*I_TWOBODYOVERLAP_Hx2y2z_Fx2y;
  abcd[160] = I_TWOBODYOVERLAP_Ixy4z_Fx2y+ABZ*I_TWOBODYOVERLAP_Hxy3z_Fx2y;
  abcd[161] = I_TWOBODYOVERLAP_Ix5z_Fx2y+ABZ*I_TWOBODYOVERLAP_Hx4z_Fx2y;
  abcd[162] = I_TWOBODYOVERLAP_I5yz_Fx2y+ABZ*I_TWOBODYOVERLAP_H5y_Fx2y;
  abcd[163] = I_TWOBODYOVERLAP_I4y2z_Fx2y+ABZ*I_TWOBODYOVERLAP_H4yz_Fx2y;
  abcd[164] = I_TWOBODYOVERLAP_I3y3z_Fx2y+ABZ*I_TWOBODYOVERLAP_H3y2z_Fx2y;
  abcd[165] = I_TWOBODYOVERLAP_I2y4z_Fx2y+ABZ*I_TWOBODYOVERLAP_H2y3z_Fx2y;
  abcd[166] = I_TWOBODYOVERLAP_Iy5z_Fx2y+ABZ*I_TWOBODYOVERLAP_Hy4z_Fx2y;
  abcd[167] = I_TWOBODYOVERLAP_I6z_Fx2y+ABZ*I_TWOBODYOVERLAP_H5z_Fx2y;
  abcd[168] = I_TWOBODYOVERLAP_I5xy_Fx2z+ABY*I_TWOBODYOVERLAP_H5x_Fx2z;
  abcd[169] = I_TWOBODYOVERLAP_I4x2y_Fx2z+ABY*I_TWOBODYOVERLAP_H4xy_Fx2z;
  abcd[170] = I_TWOBODYOVERLAP_I4xyz_Fx2z+ABY*I_TWOBODYOVERLAP_H4xz_Fx2z;
  abcd[171] = I_TWOBODYOVERLAP_I3x3y_Fx2z+ABY*I_TWOBODYOVERLAP_H3x2y_Fx2z;
  abcd[172] = I_TWOBODYOVERLAP_I3x2yz_Fx2z+ABY*I_TWOBODYOVERLAP_H3xyz_Fx2z;
  abcd[173] = I_TWOBODYOVERLAP_I3xy2z_Fx2z+ABY*I_TWOBODYOVERLAP_H3x2z_Fx2z;
  abcd[174] = I_TWOBODYOVERLAP_I2x4y_Fx2z+ABY*I_TWOBODYOVERLAP_H2x3y_Fx2z;
  abcd[175] = I_TWOBODYOVERLAP_I2x3yz_Fx2z+ABY*I_TWOBODYOVERLAP_H2x2yz_Fx2z;
  abcd[176] = I_TWOBODYOVERLAP_I2x2y2z_Fx2z+ABY*I_TWOBODYOVERLAP_H2xy2z_Fx2z;
  abcd[177] = I_TWOBODYOVERLAP_I2xy3z_Fx2z+ABY*I_TWOBODYOVERLAP_H2x3z_Fx2z;
  abcd[178] = I_TWOBODYOVERLAP_Ix5y_Fx2z+ABY*I_TWOBODYOVERLAP_Hx4y_Fx2z;
  abcd[179] = I_TWOBODYOVERLAP_Ix4yz_Fx2z+ABY*I_TWOBODYOVERLAP_Hx3yz_Fx2z;
  abcd[180] = I_TWOBODYOVERLAP_Ix3y2z_Fx2z+ABY*I_TWOBODYOVERLAP_Hx2y2z_Fx2z;
  abcd[181] = I_TWOBODYOVERLAP_Ix2y3z_Fx2z+ABY*I_TWOBODYOVERLAP_Hxy3z_Fx2z;
  abcd[182] = I_TWOBODYOVERLAP_Ixy4z_Fx2z+ABY*I_TWOBODYOVERLAP_Hx4z_Fx2z;
  abcd[183] = I_TWOBODYOVERLAP_I6y_Fx2z+ABY*I_TWOBODYOVERLAP_H5y_Fx2z;
  abcd[184] = I_TWOBODYOVERLAP_I5yz_Fx2z+ABY*I_TWOBODYOVERLAP_H4yz_Fx2z;
  abcd[185] = I_TWOBODYOVERLAP_I4y2z_Fx2z+ABY*I_TWOBODYOVERLAP_H3y2z_Fx2z;
  abcd[186] = I_TWOBODYOVERLAP_I3y3z_Fx2z+ABY*I_TWOBODYOVERLAP_H2y3z_Fx2z;
  abcd[187] = I_TWOBODYOVERLAP_I2y4z_Fx2z+ABY*I_TWOBODYOVERLAP_Hy4z_Fx2z;
  abcd[188] = I_TWOBODYOVERLAP_Iy5z_Fx2z+ABY*I_TWOBODYOVERLAP_H5z_Fx2z;
  abcd[189] = I_TWOBODYOVERLAP_I6x_F3z+ABX*I_TWOBODYOVERLAP_H5x_F3z;
  abcd[190] = I_TWOBODYOVERLAP_I5xy_F3z+ABX*I_TWOBODYOVERLAP_H4xy_F3z;
  abcd[191] = I_TWOBODYOVERLAP_I5xz_F3z+ABX*I_TWOBODYOVERLAP_H4xz_F3z;
  abcd[192] = I_TWOBODYOVERLAP_I4x2y_F3z+ABX*I_TWOBODYOVERLAP_H3x2y_F3z;
  abcd[193] = I_TWOBODYOVERLAP_I4xyz_F3z+ABX*I_TWOBODYOVERLAP_H3xyz_F3z;
  abcd[194] = I_TWOBODYOVERLAP_I4x2z_F3z+ABX*I_TWOBODYOVERLAP_H3x2z_F3z;
  abcd[195] = I_TWOBODYOVERLAP_I3x3y_F3z+ABX*I_TWOBODYOVERLAP_H2x3y_F3z;
  abcd[196] = I_TWOBODYOVERLAP_I3x2yz_F3z+ABX*I_TWOBODYOVERLAP_H2x2yz_F3z;
  abcd[197] = I_TWOBODYOVERLAP_I3xy2z_F3z+ABX*I_TWOBODYOVERLAP_H2xy2z_F3z;
  abcd[198] = I_TWOBODYOVERLAP_I3x3z_F3z+ABX*I_TWOBODYOVERLAP_H2x3z_F3z;
  abcd[199] = I_TWOBODYOVERLAP_I2x4y_F3z+ABX*I_TWOBODYOVERLAP_Hx4y_F3z;
  abcd[200] = I_TWOBODYOVERLAP_I2x3yz_F3z+ABX*I_TWOBODYOVERLAP_Hx3yz_F3z;
  abcd[201] = I_TWOBODYOVERLAP_I2x2y2z_F3z+ABX*I_TWOBODYOVERLAP_Hx2y2z_F3z;
  abcd[202] = I_TWOBODYOVERLAP_I2xy3z_F3z+ABX*I_TWOBODYOVERLAP_Hxy3z_F3z;
  abcd[203] = I_TWOBODYOVERLAP_I2x4z_F3z+ABX*I_TWOBODYOVERLAP_Hx4z_F3z;
  abcd[204] = I_TWOBODYOVERLAP_Ix5y_F3z+ABX*I_TWOBODYOVERLAP_H5y_F3z;
  abcd[205] = I_TWOBODYOVERLAP_Ix4yz_F3z+ABX*I_TWOBODYOVERLAP_H4yz_F3z;
  abcd[206] = I_TWOBODYOVERLAP_Ix3y2z_F3z+ABX*I_TWOBODYOVERLAP_H3y2z_F3z;
  abcd[207] = I_TWOBODYOVERLAP_Ix2y3z_F3z+ABX*I_TWOBODYOVERLAP_H2y3z_F3z;
  abcd[208] = I_TWOBODYOVERLAP_Ixy4z_F3z+ABX*I_TWOBODYOVERLAP_Hy4z_F3z;
  abcd[209] = I_TWOBODYOVERLAP_Ix5z_F3z+ABX*I_TWOBODYOVERLAP_H5z_F3z;
  abcd[210] = I_TWOBODYOVERLAP_I5xy_F3y+ABY*I_TWOBODYOVERLAP_H5x_F3y;
  abcd[211] = I_TWOBODYOVERLAP_I4x2y_F3y+ABY*I_TWOBODYOVERLAP_H4xy_F3y;
  abcd[212] = I_TWOBODYOVERLAP_I4xyz_F3y+ABY*I_TWOBODYOVERLAP_H4xz_F3y;
  abcd[213] = I_TWOBODYOVERLAP_I3x3y_F3y+ABY*I_TWOBODYOVERLAP_H3x2y_F3y;
  abcd[214] = I_TWOBODYOVERLAP_I3x2yz_F3y+ABY*I_TWOBODYOVERLAP_H3xyz_F3y;
  abcd[215] = I_TWOBODYOVERLAP_I3xy2z_F3y+ABY*I_TWOBODYOVERLAP_H3x2z_F3y;
  abcd[216] = I_TWOBODYOVERLAP_I2x4y_F3y+ABY*I_TWOBODYOVERLAP_H2x3y_F3y;
  abcd[217] = I_TWOBODYOVERLAP_I2x3yz_F3y+ABY*I_TWOBODYOVERLAP_H2x2yz_F3y;
  abcd[218] = I_TWOBODYOVERLAP_I2x2y2z_F3y+ABY*I_TWOBODYOVERLAP_H2xy2z_F3y;
  abcd[219] = I_TWOBODYOVERLAP_I2xy3z_F3y+ABY*I_TWOBODYOVERLAP_H2x3z_F3y;
  abcd[220] = I_TWOBODYOVERLAP_Ix5y_F3y+ABY*I_TWOBODYOVERLAP_Hx4y_F3y;
  abcd[221] = I_TWOBODYOVERLAP_Ix4yz_F3y+ABY*I_TWOBODYOVERLAP_Hx3yz_F3y;
  abcd[222] = I_TWOBODYOVERLAP_Ix3y2z_F3y+ABY*I_TWOBODYOVERLAP_Hx2y2z_F3y;
  abcd[223] = I_TWOBODYOVERLAP_Ix2y3z_F3y+ABY*I_TWOBODYOVERLAP_Hxy3z_F3y;
  abcd[224] = I_TWOBODYOVERLAP_Ixy4z_F3y+ABY*I_TWOBODYOVERLAP_Hx4z_F3y;
  abcd[225] = I_TWOBODYOVERLAP_I6y_F3y+ABY*I_TWOBODYOVERLAP_H5y_F3y;
  abcd[226] = I_TWOBODYOVERLAP_I5yz_F3y+ABY*I_TWOBODYOVERLAP_H4yz_F3y;
  abcd[227] = I_TWOBODYOVERLAP_I4y2z_F3y+ABY*I_TWOBODYOVERLAP_H3y2z_F3y;
  abcd[228] = I_TWOBODYOVERLAP_I3y3z_F3y+ABY*I_TWOBODYOVERLAP_H2y3z_F3y;
  abcd[229] = I_TWOBODYOVERLAP_I2y4z_F3y+ABY*I_TWOBODYOVERLAP_Hy4z_F3y;
  abcd[230] = I_TWOBODYOVERLAP_Iy5z_F3y+ABY*I_TWOBODYOVERLAP_H5z_F3y;
  abcd[231] = I_TWOBODYOVERLAP_I5xz_F3y+ABZ*I_TWOBODYOVERLAP_H5x_F3y;
  abcd[232] = I_TWOBODYOVERLAP_I4xyz_F3y+ABZ*I_TWOBODYOVERLAP_H4xy_F3y;
  abcd[233] = I_TWOBODYOVERLAP_I4x2z_F3y+ABZ*I_TWOBODYOVERLAP_H4xz_F3y;
  abcd[234] = I_TWOBODYOVERLAP_I3x2yz_F3y+ABZ*I_TWOBODYOVERLAP_H3x2y_F3y;
  abcd[235] = I_TWOBODYOVERLAP_I3xy2z_F3y+ABZ*I_TWOBODYOVERLAP_H3xyz_F3y;
  abcd[236] = I_TWOBODYOVERLAP_I3x3z_F3y+ABZ*I_TWOBODYOVERLAP_H3x2z_F3y;
  abcd[237] = I_TWOBODYOVERLAP_I2x3yz_F3y+ABZ*I_TWOBODYOVERLAP_H2x3y_F3y;
  abcd[238] = I_TWOBODYOVERLAP_I2x2y2z_F3y+ABZ*I_TWOBODYOVERLAP_H2x2yz_F3y;
  abcd[239] = I_TWOBODYOVERLAP_I2xy3z_F3y+ABZ*I_TWOBODYOVERLAP_H2xy2z_F3y;
  abcd[240] = I_TWOBODYOVERLAP_I2x4z_F3y+ABZ*I_TWOBODYOVERLAP_H2x3z_F3y;
  abcd[241] = I_TWOBODYOVERLAP_Ix4yz_F3y+ABZ*I_TWOBODYOVERLAP_Hx4y_F3y;
  abcd[242] = I_TWOBODYOVERLAP_Ix3y2z_F3y+ABZ*I_TWOBODYOVERLAP_Hx3yz_F3y;
  abcd[243] = I_TWOBODYOVERLAP_Ix2y3z_F3y+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F3y;
  abcd[244] = I_TWOBODYOVERLAP_Ixy4z_F3y+ABZ*I_TWOBODYOVERLAP_Hxy3z_F3y;
  abcd[245] = I_TWOBODYOVERLAP_Ix5z_F3y+ABZ*I_TWOBODYOVERLAP_Hx4z_F3y;
  abcd[246] = I_TWOBODYOVERLAP_I5yz_F3y+ABZ*I_TWOBODYOVERLAP_H5y_F3y;
  abcd[247] = I_TWOBODYOVERLAP_I4y2z_F3y+ABZ*I_TWOBODYOVERLAP_H4yz_F3y;
  abcd[248] = I_TWOBODYOVERLAP_I3y3z_F3y+ABZ*I_TWOBODYOVERLAP_H3y2z_F3y;
  abcd[249] = I_TWOBODYOVERLAP_I2y4z_F3y+ABZ*I_TWOBODYOVERLAP_H2y3z_F3y;
  abcd[250] = I_TWOBODYOVERLAP_Iy5z_F3y+ABZ*I_TWOBODYOVERLAP_Hy4z_F3y;
  abcd[251] = I_TWOBODYOVERLAP_I6z_F3y+ABZ*I_TWOBODYOVERLAP_H5z_F3y;
  abcd[252] = I_TWOBODYOVERLAP_I5xz_F2yz+ABZ*I_TWOBODYOVERLAP_H5x_F2yz;
  abcd[253] = I_TWOBODYOVERLAP_I4xyz_F2yz+ABZ*I_TWOBODYOVERLAP_H4xy_F2yz;
  abcd[254] = I_TWOBODYOVERLAP_I4x2z_F2yz+ABZ*I_TWOBODYOVERLAP_H4xz_F2yz;
  abcd[255] = I_TWOBODYOVERLAP_I3x2yz_F2yz+ABZ*I_TWOBODYOVERLAP_H3x2y_F2yz;
  abcd[256] = I_TWOBODYOVERLAP_I3xy2z_F2yz+ABZ*I_TWOBODYOVERLAP_H3xyz_F2yz;
  abcd[257] = I_TWOBODYOVERLAP_I3x3z_F2yz+ABZ*I_TWOBODYOVERLAP_H3x2z_F2yz;
  abcd[258] = I_TWOBODYOVERLAP_I2x3yz_F2yz+ABZ*I_TWOBODYOVERLAP_H2x3y_F2yz;
  abcd[259] = I_TWOBODYOVERLAP_I2x2y2z_F2yz+ABZ*I_TWOBODYOVERLAP_H2x2yz_F2yz;
  abcd[260] = I_TWOBODYOVERLAP_I2xy3z_F2yz+ABZ*I_TWOBODYOVERLAP_H2xy2z_F2yz;
  abcd[261] = I_TWOBODYOVERLAP_I2x4z_F2yz+ABZ*I_TWOBODYOVERLAP_H2x3z_F2yz;
  abcd[262] = I_TWOBODYOVERLAP_Ix4yz_F2yz+ABZ*I_TWOBODYOVERLAP_Hx4y_F2yz;
  abcd[263] = I_TWOBODYOVERLAP_Ix3y2z_F2yz+ABZ*I_TWOBODYOVERLAP_Hx3yz_F2yz;
  abcd[264] = I_TWOBODYOVERLAP_Ix2y3z_F2yz+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F2yz;
  abcd[265] = I_TWOBODYOVERLAP_Ixy4z_F2yz+ABZ*I_TWOBODYOVERLAP_Hxy3z_F2yz;
  abcd[266] = I_TWOBODYOVERLAP_Ix5z_F2yz+ABZ*I_TWOBODYOVERLAP_Hx4z_F2yz;
  abcd[267] = I_TWOBODYOVERLAP_I5yz_F2yz+ABZ*I_TWOBODYOVERLAP_H5y_F2yz;
  abcd[268] = I_TWOBODYOVERLAP_I4y2z_F2yz+ABZ*I_TWOBODYOVERLAP_H4yz_F2yz;
  abcd[269] = I_TWOBODYOVERLAP_I3y3z_F2yz+ABZ*I_TWOBODYOVERLAP_H3y2z_F2yz;
  abcd[270] = I_TWOBODYOVERLAP_I2y4z_F2yz+ABZ*I_TWOBODYOVERLAP_H2y3z_F2yz;
  abcd[271] = I_TWOBODYOVERLAP_Iy5z_F2yz+ABZ*I_TWOBODYOVERLAP_Hy4z_F2yz;
  abcd[272] = I_TWOBODYOVERLAP_I6z_F2yz+ABZ*I_TWOBODYOVERLAP_H5z_F2yz;
  abcd[273] = I_TWOBODYOVERLAP_I5xy_F3z+ABY*I_TWOBODYOVERLAP_H5x_F3z;
  abcd[274] = I_TWOBODYOVERLAP_I4x2y_F3z+ABY*I_TWOBODYOVERLAP_H4xy_F3z;
  abcd[275] = I_TWOBODYOVERLAP_I4xyz_F3z+ABY*I_TWOBODYOVERLAP_H4xz_F3z;
  abcd[276] = I_TWOBODYOVERLAP_I3x3y_F3z+ABY*I_TWOBODYOVERLAP_H3x2y_F3z;
  abcd[277] = I_TWOBODYOVERLAP_I3x2yz_F3z+ABY*I_TWOBODYOVERLAP_H3xyz_F3z;
  abcd[278] = I_TWOBODYOVERLAP_I3xy2z_F3z+ABY*I_TWOBODYOVERLAP_H3x2z_F3z;
  abcd[279] = I_TWOBODYOVERLAP_I2x4y_F3z+ABY*I_TWOBODYOVERLAP_H2x3y_F3z;
  abcd[280] = I_TWOBODYOVERLAP_I2x3yz_F3z+ABY*I_TWOBODYOVERLAP_H2x2yz_F3z;
  abcd[281] = I_TWOBODYOVERLAP_I2x2y2z_F3z+ABY*I_TWOBODYOVERLAP_H2xy2z_F3z;
  abcd[282] = I_TWOBODYOVERLAP_I2xy3z_F3z+ABY*I_TWOBODYOVERLAP_H2x3z_F3z;
  abcd[283] = I_TWOBODYOVERLAP_Ix5y_F3z+ABY*I_TWOBODYOVERLAP_Hx4y_F3z;
  abcd[284] = I_TWOBODYOVERLAP_Ix4yz_F3z+ABY*I_TWOBODYOVERLAP_Hx3yz_F3z;
  abcd[285] = I_TWOBODYOVERLAP_Ix3y2z_F3z+ABY*I_TWOBODYOVERLAP_Hx2y2z_F3z;
  abcd[286] = I_TWOBODYOVERLAP_Ix2y3z_F3z+ABY*I_TWOBODYOVERLAP_Hxy3z_F3z;
  abcd[287] = I_TWOBODYOVERLAP_Ixy4z_F3z+ABY*I_TWOBODYOVERLAP_Hx4z_F3z;
  abcd[288] = I_TWOBODYOVERLAP_I6y_F3z+ABY*I_TWOBODYOVERLAP_H5y_F3z;
  abcd[289] = I_TWOBODYOVERLAP_I5yz_F3z+ABY*I_TWOBODYOVERLAP_H4yz_F3z;
  abcd[290] = I_TWOBODYOVERLAP_I4y2z_F3z+ABY*I_TWOBODYOVERLAP_H3y2z_F3z;
  abcd[291] = I_TWOBODYOVERLAP_I3y3z_F3z+ABY*I_TWOBODYOVERLAP_H2y3z_F3z;
  abcd[292] = I_TWOBODYOVERLAP_I2y4z_F3z+ABY*I_TWOBODYOVERLAP_Hy4z_F3z;
  abcd[293] = I_TWOBODYOVERLAP_Iy5z_F3z+ABY*I_TWOBODYOVERLAP_H5z_F3z;
  abcd[294] = I_TWOBODYOVERLAP_I5xz_F3z+ABZ*I_TWOBODYOVERLAP_H5x_F3z;
  abcd[295] = I_TWOBODYOVERLAP_I4xyz_F3z+ABZ*I_TWOBODYOVERLAP_H4xy_F3z;
  abcd[296] = I_TWOBODYOVERLAP_I4x2z_F3z+ABZ*I_TWOBODYOVERLAP_H4xz_F3z;
  abcd[297] = I_TWOBODYOVERLAP_I3x2yz_F3z+ABZ*I_TWOBODYOVERLAP_H3x2y_F3z;
  abcd[298] = I_TWOBODYOVERLAP_I3xy2z_F3z+ABZ*I_TWOBODYOVERLAP_H3xyz_F3z;
  abcd[299] = I_TWOBODYOVERLAP_I3x3z_F3z+ABZ*I_TWOBODYOVERLAP_H3x2z_F3z;
  abcd[300] = I_TWOBODYOVERLAP_I2x3yz_F3z+ABZ*I_TWOBODYOVERLAP_H2x3y_F3z;
  abcd[301] = I_TWOBODYOVERLAP_I2x2y2z_F3z+ABZ*I_TWOBODYOVERLAP_H2x2yz_F3z;
  abcd[302] = I_TWOBODYOVERLAP_I2xy3z_F3z+ABZ*I_TWOBODYOVERLAP_H2xy2z_F3z;
  abcd[303] = I_TWOBODYOVERLAP_I2x4z_F3z+ABZ*I_TWOBODYOVERLAP_H2x3z_F3z;
  abcd[304] = I_TWOBODYOVERLAP_Ix4yz_F3z+ABZ*I_TWOBODYOVERLAP_Hx4y_F3z;
  abcd[305] = I_TWOBODYOVERLAP_Ix3y2z_F3z+ABZ*I_TWOBODYOVERLAP_Hx3yz_F3z;
  abcd[306] = I_TWOBODYOVERLAP_Ix2y3z_F3z+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F3z;
  abcd[307] = I_TWOBODYOVERLAP_Ixy4z_F3z+ABZ*I_TWOBODYOVERLAP_Hxy3z_F3z;
  abcd[308] = I_TWOBODYOVERLAP_Ix5z_F3z+ABZ*I_TWOBODYOVERLAP_Hx4z_F3z;
  abcd[309] = I_TWOBODYOVERLAP_I5yz_F3z+ABZ*I_TWOBODYOVERLAP_H5y_F3z;
  abcd[310] = I_TWOBODYOVERLAP_I4y2z_F3z+ABZ*I_TWOBODYOVERLAP_H4yz_F3z;
  abcd[311] = I_TWOBODYOVERLAP_I3y3z_F3z+ABZ*I_TWOBODYOVERLAP_H3y2z_F3z;
  abcd[312] = I_TWOBODYOVERLAP_I2y4z_F3z+ABZ*I_TWOBODYOVERLAP_H2y3z_F3z;
  abcd[313] = I_TWOBODYOVERLAP_Iy5z_F3z+ABZ*I_TWOBODYOVERLAP_Hy4z_F3z;
  abcd[314] = I_TWOBODYOVERLAP_I6z_F3z+ABZ*I_TWOBODYOVERLAP_H5z_F3z;
}
